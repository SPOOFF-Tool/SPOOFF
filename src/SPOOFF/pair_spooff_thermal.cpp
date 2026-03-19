// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
 *
 * *** Smoothed Particle Open-source in-Operation Fuel Fragmentation ***
 *
 * The SPOOFF package is a modification of the MACHDYN package specifically 
 * designed for modelling fragmentation within nuclear fuels. This was 
 * created in collaboration betweein UKNNL and University of Bangor Nuclear 
 * Futures Centre under the OperaHPC project.
 * (2026) matthew.horton@uknnl.com
 * 
 * This file is modified from files part of the MACHDYN package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 https://www.lammps.org/, Sandia National Laboratories
 LAMMPS development team: developers@lammps.org

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mike Parks (SNL)
------------------------------------------------------------------------- */

#include "pair_spooff_thermal.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "spooff_kernels.h"
#include "comm.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace SPOOFF_Kernels;

#define SQRT2 1.414213562e0

/* ---------------------------------------------------------------------- */

PairThermal::PairThermal(LAMMPS *lmp) :
                Pair(lmp) {

        onerad_dynamic = onerad_frozen = maxrad_dynamic = maxrad_frozen = nullptr;
        comm_forward = 4; // this pair style communicates 20 doubles to ghost atoms : PK1 tensor + F tensor + shepardWeight   //MH added cpsph,ksph,hgsph,fracsph to forward comm
        cut_comm = MAX(neighbor->cutneighmax, comm->cutghostuser); // cutoff radius within which ghost atoms are communicated.
}

/* ---------------------------------------------------------------------- */

PairThermal::~PairThermal() {

        if (allocated) {
                memory->destroy(setflag);
                memory->destroy(cutsq);
                //memory->destroy(impact_velocity);

                delete[] onerad_dynamic;
                delete[] onerad_frozen;
                delete[] maxrad_dynamic;
                delete[] maxrad_frozen;
        }
}

/* ---------------------------------------------------------------------- */

void PairThermal::compute(int eflag, int vflag) {
        int i, j, ii, jj, inum, jnum, itype, jtype;
        double xtmp, ytmp, ztmp, delx, dely, delz;
        double vxtmp, vytmp, vztmp, delvx, delvy, delvz;
        double rsq, r, evdwl, fpair;
        int *ilist, *jlist, *numneigh, **firstneigh;
        double rcut, r_geom, delta, dt_crit;
        double *rmass = atom->rmass;

        evdwl = 0.0;
        ev_init(eflag, vflag);

        tagint *mol = atom->molecule;
        double **f = atom->f;
        double **x = atom->x;
        double **v = atom->vest;
        double **x0 = atom->x0;
        double *esph = atom->esph;
        double *hgsph = atom->hgsph;
        double *tsph = atom->tsph;
        double *rhosph = atom->rhosph;
        double *ksph = atom->ksph;
        double *cpsph = atom->cpsph;
        double *desph = atom->desph;
        double **qsph = atom->qsph;
        double *vfrac = atom->vfrac;
        int *type = atom->type;
        int nlocal = atom->nlocal;
        double *radius = atom->contact_radius;
        double *sph_radius = atom->radius;
        double rcutSq;
        double delx0, dely0, delz0, rSq0, sphCut, h, wf, wfd, deltaHE, deltaQ, E_reduced;

        int newton_pair = force->newton_pair;
        int periodic = (domain->xperiodic || domain->yperiodic || domain->zperiodic);

        inum = list->inum;
        ilist = list->ilist;
        numneigh = list->numneigh;
        firstneigh = list->firstneigh;

        stable_time_increment = 1.0e22;

        comm->forward_comm(this);

        // loop over neighbors of my atoms
        for (ii = 0; ii < inum; ii++) {
                i = ilist[ii];
                xtmp = x[i][0];
                ytmp = x[i][1];
                ztmp = x[i][2];
                vxtmp = v[i][0];
                vytmp = v[i][1];
                vztmp = v[i][2];
                itype = type[i];
                jlist = firstneigh[i];
                jnum = numneigh[i];

                for (jj = 0; jj < jnum; jj++) {
                        j = jlist[jj];
                        j &= NEIGHMASK;

                        jtype = type[j];

                        delx = xtmp - x[j][0];
                        dely = ytmp - x[j][1];
                        delz = ztmp - x[j][2];
                        delvx = vxtmp - v[j][0];
                        delvy = vytmp - v[j][1];
                        delvz = vztmp - v[j][2];

                        rsq = delx * delx + dely * dely + delz * delz;
                        
                        r = sqrt(rsq);
                        h = (sph_radius[i] + sph_radius[j]);
                        // heat flow
                        if (r < h) {
                                spiky_kernel_and_derivative(h, r, domain->dimension, wf, wfd);

                                if(!first){ // crude way to check the system has initialized
                                    
                                deltaHE =  2.0*rmass[i] * rmass[j] / (rmass[i]+rmass[j]);
                                deltaHE *=  (rhosph[i] + rhosph[j]) / (rhosph[i] * rhosph[j]);
                                deltaHE *=  2/((1/ksph[i])+(1/ksph[j]));  //harmonic mean
                                deltaHE *=  (cpsph[i] + cpsph[j])/(2*cpsph[i] * cpsph[j]);  //harmonic mean
                                deltaHE *= (1/(rhosph[i])) * (esph[i] - esph[j]) * wfd/r; 

                                if(false){
                                printf("i,j: %d %d \n",i,j);  //added for testing
                                printf("r: %f \n",r);  //added for testing
                                printf("cpsph: %f %f \n",cpsph[i],cpsph[j]);  //added for testing
                                printf("ksph: %f %f \n",ksph[i],ksph[j]);  //added for testing
                                printf("rhosph: %f %f \n",rhosph[i],rhosph[j]);  //added for testing
                                printf("wfd: %f \n",wfd);  //added for testing
                                printf("qsph: %f %f \n",qsph[i][0],qsph[j][0]);  //added for testing
                                printf("tsph: %f %f \n",tsph[i],tsph[j]);  //added for testing
                                printf("deltaQ: %f \n",deltaQ);  //added for testing
                                printf("deltaHE: %f \n",deltaHE);  //added for testing
                                printf("------------------------------------");  //added for testing
                                }
                                // calculate heatflux
                                deltaQ = rmass[i] * rmass[j] / (rmass[i]+rmass[j]);
                                deltaQ *=  (rhosph[i] + rhosph[j]) / (rhosph[i] * rhosph[j]);
                                deltaQ *= (-2/((1/ksph[i])+(1/ksph[j]))) * (tsph[i] - tsph[j]) * wfd/r; 


                                //update rate of change
                                desph[i] += deltaHE;
                                qsph[i][0] += deltaQ*delx; //missing correction as not calculated for thermal pair type
                                qsph[i][1] += deltaQ*dely;
                                qsph[i][2] += deltaQ*delz;

                                if (newton_pair || j < nlocal) {
                                        desph[j] -= deltaHE;
                                        qsph[j][0] -= deltaQ*delx; //missing correction as not calculated for thermal pair type
                                        qsph[j][1] -= deltaQ*dely;
                                        qsph[j][2] -= deltaQ*delz;
                                }
                        }


                        }
                                
                }
                desph[i] += (hgsph[i]*rmass[i]);
                
        }
        first = false;

        //for (i = 0; i < nlocal; i++)
                //desph[i] += (hgsph[i]*rmass[i]);

//      double stable_time_increment_all = 0.0;
//      MPI_Allreduce(&stable_time_increment, &stable_time_increment_all, 1, MPI_DOUBLE, MPI_MIN, world);
//      if (comm->me == 0) {
//              printf("stable time step for pair spooff/thermal is %f\n", stable_time_increment_all);
//      }
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairThermal::allocate() {
        allocated = 1;
        int n = atom->ntypes;

        memory->create(setflag, n + 1, n + 1, "pair:setflag");
        for (int i = 1; i <= n; i++)
                for (int j = i; j <= n; j++)
                        setflag[i][j] = 0;

        memory->create(cutsq, n + 1, n + 1, "pair:cutsq"); // always needs to be allocated, even with granular neighborlist


        onerad_dynamic = new double[n + 1];
        onerad_frozen = new double[n + 1];
        maxrad_dynamic = new double[n + 1];
        maxrad_frozen = new double[n + 1];

        //int nmax = atom->nmax;
        //int *numneigh, ii, i, inum, *ilist;
        //numneigh = list->numneigh;
        //inum = list->inum;
        //ilist = list->ilist;
        //int maxpartner = 1;
        //for (ii = 0; ii < inum; ii++) {
        //        i = ilist[ii];
        //        maxpartner = MAX(numneigh[i],maxpartner);
        //}
        //int maxall;
        //MPI_Allreduce(&maxpartner, &maxall, 1, MPI_INT, MPI_MAX, world);
        //maxpartner = maxall;
        //memory->create(impact_velocity, nmax, maxpartner, "pair:impact_velocity");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairThermal::settings(int narg, char **arg) {
        if (narg != 0)
                error->all(FLERR, "Illegal number of args for pair_style thermal");

}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairThermal::coeff(int narg, char **arg) {
        if (narg != 2)
                error->all(FLERR, "Incorrect args for pair coefficients for Thermal style, Need only pair coeficients");
        if (!allocated)
                allocate();

        int ilo, ihi, jlo, jhi;
        utils::bounds(FLERR,arg[0], 1, atom->ntypes, ilo, ihi, error);
        utils::bounds(FLERR,arg[1], 1, atom->ntypes, jlo, jhi, error);


        int count = 0;
        for (int i = ilo; i <= ihi; i++) {
                for (int j = MAX(jlo, i); j <= jhi; j++) {
                        setflag[i][j] = 1;
                        count++;
                }
        }

        if (count == 0)
                error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairThermal::init_one(int i, int j) {

        if (!allocated)
                allocate();

        if (setflag[i][j] == 0)
                error->all(FLERR, "All pair coeffs are not set");


        // cutoff = sum of max I,J radii for
        // dynamic/dynamic & dynamic/frozen interactions, but not frozen/frozen

        double cutoff = maxrad_dynamic[i] + maxrad_dynamic[j];
        cutoff = MAX(cutoff, maxrad_frozen[i] + maxrad_dynamic[j]);
        cutoff = MAX(cutoff, maxrad_dynamic[i] + maxrad_frozen[j]);

        if (comm->me == 0) {
                printf("cutoff for pair spooff/thermal = %f\n", cutoff);
        }
        return cutoff;
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairThermal::init_style() {
        int i;

        // error checks

        if (!atom->contact_radius_flag)
                error->all(FLERR, "Pair style spooff/thermal requires atom style with contact_radius");

        neighbor->add_request(this, NeighConst::REQ_SIZE);

        // set maxrad_dynamic and maxrad_frozen for each type
        // include future Fix pour particles as dynamic

        for (i = 1; i <= atom->ntypes; i++)
                onerad_dynamic[i] = onerad_frozen[i] = 0.0;

        double *radius = atom->radius;
        int *type = atom->type;
        int nlocal = atom->nlocal;

        for (i = 0; i < nlocal; i++) {
                onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]], radius[i]);
        }

        MPI_Allreduce(&onerad_dynamic[1], &maxrad_dynamic[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);
        MPI_Allreduce(&onerad_frozen[1], &maxrad_frozen[1], atom->ntypes, MPI_DOUBLE, MPI_MAX, world);
}

/* ----------------------------------------------------------------------
 neighbor callback to inform pair style of neighbor list to use
 optional granular history list
 ------------------------------------------------------------------------- */

void PairThermal::init_list(int id, NeighList *ptr) {
        if (id == 0)
                list = ptr;
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based arrays
 ------------------------------------------------------------------------- */

double PairThermal::memory_usage() {
        return 1.0;
}

void *PairThermal::extract(const char *str, int &/*i*/) {
        //printf("in PairTriSurf::extract\n");
        if (strcmp(str, "spooff/thermal/stable_time_increment_ptr") == 0) {
                return (void *) &stable_time_increment;
        }

        return nullptr;

}

/* ----------------------------------------------------------------------
 Communication of ghost atoms for mpi
 ------------------------------------------------------------------------- */

int PairThermal::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/) {
  int i, j, m;
  double *cpsph = atom->cpsph;
  double *ksph = atom->ksph;
  double *rhosph = atom->rhosph;
  tagint *mol = atom->molecule;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = cpsph[j]; //1
    buf[m++] = ksph[j]; //2
    buf[m++] = rhosph[j]; //3
    buf[m++] = mol[j]; //4
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairThermal::unpack_forward_comm(int n, int first, double *buf) {
  int i, m, last;
  double *cpsph = atom->cpsph;
  double *ksph = atom->ksph;
  double *rhosph = atom->rhosph;
  tagint *mol = atom->molecule;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    cpsph[i] = buf[m++]; //1
    ksph[i] = buf[m++];  //2
    rhosph[i] = buf[m++];  //3
    mol[i] = static_cast<int>(buf[m++]); //4
  }
}