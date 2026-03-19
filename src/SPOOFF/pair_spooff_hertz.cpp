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

#include "pair_spooff_hertz.h"

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

PairHertz::PairHertz(LAMMPS *lmp) :
                Pair(lmp) {

        onerad_dynamic = onerad_frozen = maxrad_dynamic = maxrad_frozen = nullptr;
        bulkmodulus = nullptr;
        damping_scaling = nullptr;
        restitution = nullptr;
        kn = nullptr;
        //impact_velocity = nullptr;
        scale = 1.0;
        comm_forward = 3; // this pair style communicates 20 doubles to ghost atoms : PK1 tensor + F tensor + shepardWeight   //MH added cpsph,ksph,hgsph,fracsph to forward comm
        cut_comm = MAX(neighbor->cutneighmax, comm->cutghostuser); // cutoff radius within which ghost atoms are communicated.
}

/* ---------------------------------------------------------------------- */

PairHertz::~PairHertz() {

        if (allocated) {
                memory->destroy(setflag);
                memory->destroy(cutsq);
                memory->destroy(bulkmodulus);
                memory->destroy(restitution);
                memory->destroy(damping_scaling);
                memory->destroy(kn);
                //memory->destroy(impact_velocity);

                delete[] onerad_dynamic;
                delete[] onerad_frozen;
                delete[] maxrad_dynamic;
                delete[] maxrad_frozen;
        }
}

/* ---------------------------------------------------------------------- */

void PairHertz::compute(int eflag, int vflag) {
        int i, j, ii, jj, inum, jnum, itype, jtype;
        double xtmp, ytmp, ztmp, delx, dely, delz;
        double vxtmp, vytmp, vztmp, delvx, delvy, delvz;
        double rsq, r, evdwl, fpair;
        int *ilist, *jlist, *numneigh, **firstneigh;
        double rcut, r_geom, delta, ri, rj, dt_crit;
        double *rmass = atom->rmass;

        evdwl = 0.0;
        ev_init(eflag, vflag);

        tagint *mol = atom->molecule;
        double **f = atom->f;
        double **x = atom->x;
        double **v = atom->vest;
        double **x0 = atom->x0;
        double *youngs = atom->youngs;
        double *poissons = atom->poissons;
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
                ri = scale * radius[i];
                jlist = firstneigh[i];
                jnum = numneigh[i];

                for (jj = 0; jj < jnum; jj++) {
                        if(first){
                                //impact_velocity[i][jj]=-1.0;
                        }
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

                        rj = scale * radius[j];
                        rcut = ri + rj;
                        rcutSq = rcut * rcut;
                        if (rsq < rcutSq) {



                                /*
                                 * self contact option:
                                 * if pair of particles was initially close enough to interact via a bulk continuum mechanism (e.g. SPH), exclude pair from contact forces.
                                 * this approach should work well if no updates of the reference configuration are performed.
                                 */

                                if (itype == jtype) {
                                        if (mol[i] > 0 & mol[j] > 0) { //only considers non deleted particles in reference
                                        delx0 = x0[j][0] - x0[i][0];
                                        dely0 = x0[j][1] - x0[i][1];
                                        delz0 = x0[j][2] - x0[i][2];
                                        if (periodic) {
                                                domain->minimum_image(delx0, dely0, delz0);
                                        }
                                        rSq0 = delx0 * delx0 + dely0 * dely0 + delz0 * delz0; // initial distance
                                        sphCut = sph_radius[i] + sph_radius[j];
                                        if (rSq0 < sphCut * sphCut) {
                                                rcut = 0.5 * rcut;
                                                rcutSq = rcut * rcut;
                                                if (rsq > rcutSq) {
                                                        continue;
                                                }
                                        }
                                        }
                                }

                                r = sqrt(rsq);
                                //printf("hertz interaction, r=%f, cut=%f, h=%f\n", r, rcut, sqrt(rSq0));

                                // Hertzian short-range forces
                                delta = rcut - r; // overlap distance
                                r_geom = ri * rj / rcut;
                                //assuming poisson ratio = 1/4 for 3d
                                E_reduced = bulkmodulus[itype][jtype] * 1.0/(((1-pow(poissons[i],2))/youngs[i])+((1-pow(poissons[j],2))/youngs[j]));
                                fpair = 1.066666667e0 * E_reduced * delta * sqrt(delta * r_geom); //  units: N
                                //damping
                                double S, damping_coeff, velocity;
                                velocity = sqrt(delx * delx + dely * dely + delz * delz);
                                S = delx * delvx + dely * delvy + delz * delvz;
                                //if(impact_velocity[i][jj]<0.0){
                                //  impact_velocity[i][j] = velocity;
                                //}
                                if(velocity==0.0){
                                  damping_coeff = 0.0;
                                }else{
                                  damping_coeff = damping_scaling[itype][jtype]*0.587*(1-restitution[itype][jtype]);
                                  if(S<0){
                                          fpair += damping_coeff * 1.066666667e0 *  sqrt(delta * r_geom) * velocity; //moving towards each other
                                  }else{
                                          fpair -= damping_coeff * 1.066666667e0 *  sqrt(delta * r_geom) * velocity; //moving away from each other
                                  }
                                }
                                if(fpair!=fpair){
                                  printf("damping_coeff: %f, delta: %f, r_geom: %f , velocity: %f \n",damping_coeff,delta,r_geom,velocity);
                                  printf("DAMPING_SCALING: %f , RESTITUTION: %f , E_reduced: %f  \n",damping_scaling[itype][jtype],restitution[itype][jtype],E_reduced);
                                  printf("------------------------------------------- \n");
                                }
                                fpair = MAX(fpair,0.0);
                                evdwl = fpair * 0.4e0 * delta; // GCG 25 April: this expression conserves total energy
                                if(fpair==fpair){ //protects agains nan
                                dt_crit = 3.14 * sqrt(0.5 * (rmass[i] + rmass[j]) / (fpair / delta));
                                }

                                stable_time_increment = MIN(stable_time_increment, dt_crit);
                                if (r > 2.0e-16) {
                                        fpair /= r; // divide by r and multiply with non-normalized distance vector
                                } else {
                                        fpair = 0.0;
                                }

                                /*
                                 * contact viscosity -- needs to be done, see GRANULAR package for normal & shear damping
                                 * for now: no damping and thus no viscous energy deltaE
                                 */

                                if (evflag) {
                                        ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
                                }

                                f[i][0] += delx * fpair;
                                f[i][1] += dely * fpair;
                                f[i][2] += delz * fpair;

                                if (newton_pair || j < nlocal) {
                                        f[j][0] -= delx * fpair;
                                        f[j][1] -= dely * fpair;
                                        f[j][2] -= delz * fpair;
                                }


                        }
                                
                }
        }
        first = false;

//      double stable_time_increment_all = 0.0;
//      MPI_Allreduce(&stable_time_increment, &stable_time_increment_all, 1, MPI_DOUBLE, MPI_MIN, world);
//      if (comm->me == 0) {
//              printf("stable time step for pair spooff/hertz is %f\n", stable_time_increment_all);
//      }
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairHertz::allocate() {
        allocated = 1;
        int n = atom->ntypes;

        memory->create(setflag, n + 1, n + 1, "pair:setflag");
        for (int i = 1; i <= n; i++)
                for (int j = i; j <= n; j++)
                        setflag[i][j] = 0;

        memory->create(bulkmodulus, n + 1, n + 1, "pair:kspring");
        memory->create(restitution, n + 1, n + 1, "pair:retitution");
        memory->create(damping_scaling, n + 1, n + 1, "pair:damping_scaling");
        memory->create(kn, n + 1, n + 1, "pair:kn");

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

void PairHertz::settings(int narg, char **arg) {
        if (narg != 1)
                error->all(FLERR, "Illegal number of args for pair_style hertz");

        scale = utils::numeric(FLERR, arg[0],false,lmp);
        if (comm->me == 0) {
                printf("\n>>========>>========>>========>>========>>========>>========>>========>>========\n");
                printf("SPOOFF/HERTZ CONTACT SETTINGS:\n");
                printf("... effective contact radius is scaled by %f\n", scale);
                printf(">>========>>========>>========>>========>>========>>========>>========>>========\n");
        }

}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairHertz::coeff(int narg, char **arg) {
        if (narg != 5)
                error->all(FLERR, "Incorrect args for pair coefficients for Herts style, Need: Force_Scaling, Damping_Scaling, Coeficient_Restitution");
        if (!allocated)
                allocate();

        int ilo, ihi, jlo, jhi;
        utils::bounds(FLERR,arg[0], 1, atom->ntypes, ilo, ihi, error);
        utils::bounds(FLERR,arg[1], 1, atom->ntypes, jlo, jhi, error);

        double bulkmodulus_one = utils::numeric(FLERR,arg[2],false,lmp);
        double damping_scaling_one = utils::numeric(FLERR,arg[3],false,lmp);
        double restitution_one = utils::numeric(FLERR,arg[4],false,lmp);

        // set short-range force constant
        double kn_one = 0.0;
        if (domain->dimension == 3) {
                kn_one = (16. / 15.) * bulkmodulus_one; //assuming poisson ratio = 1/4 for 3d
        } else {
                kn_one = 0.251856195 * (2. / 3.) * bulkmodulus_one; //assuming poisson ratio = 1/3 for 2d
        }

        int count = 0;
        for (int i = ilo; i <= ihi; i++) {
                for (int j = MAX(jlo, i); j <= jhi; j++) {
                        bulkmodulus[i][j] = bulkmodulus_one;
                        damping_scaling[i][j] = damping_scaling_one;
                        restitution[i][j] = restitution_one;
                        kn[i][j] = kn_one;
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

double PairHertz::init_one(int i, int j) {

        if (!allocated)
                allocate();

        if (setflag[i][j] == 0)
                error->all(FLERR, "All pair coeffs are not set");

        bulkmodulus[j][i] = bulkmodulus[i][j];
        restitution[j][i] = restitution[i][j];
        damping_scaling[j][i] = damping_scaling[i][j];
        kn[j][i] = kn[i][j];

        // cutoff = sum of max I,J radii for
        // dynamic/dynamic & dynamic/frozen interactions, but not frozen/frozen

        double cutoff = maxrad_dynamic[i] + maxrad_dynamic[j];
        cutoff = MAX(cutoff, maxrad_frozen[i] + maxrad_dynamic[j]);
        cutoff = MAX(cutoff, maxrad_dynamic[i] + maxrad_frozen[j]);

        if (comm->me == 0) {
                printf("cutoff for pair spooff/hertz = %f\n", cutoff);
        }
        return cutoff;
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairHertz::init_style() {
        int i;

        // error checks

        if (!atom->contact_radius_flag)
                error->all(FLERR, "Pair style spooff/hertz requires atom style with contact_radius");

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

void PairHertz::init_list(int id, NeighList *ptr) {
        if (id == 0)
                list = ptr;
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based arrays
 ------------------------------------------------------------------------- */

double PairHertz::memory_usage() {
        return 1.0;
}

void *PairHertz::extract(const char *str, int &/*i*/) {
        //printf("in PairTriSurf::extract\n");
        if (strcmp(str, "spooff/hertz/stable_time_increment_ptr") == 0) {
                return (void *) &stable_time_increment;
        }

        return nullptr;

}

/* ----------------------------------------------------------------------
 Communication of ghost atoms for mpi
 ------------------------------------------------------------------------- */

int PairHertz::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/) {
  int i, j, m;
  double *youngs = atom->youngs;
  double *poissons = atom->poissons;
  tagint *mol = atom->molecule;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = mol[j]; //1
    buf[m++] = youngs[j]; //2
    buf[m++] = poissons[j]; //3
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairHertz::unpack_forward_comm(int n, int first, double *buf) {
  int i, m, last;
  double *youngs = atom->youngs;
  double *poissons = atom->poissons;
  tagint *mol = atom->molecule;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    mol[i] = static_cast<int>(buf[m++]); //1
    youngs[i] = buf[m++];  //2
    poissons[i] = buf[m++];  //3
  }
}