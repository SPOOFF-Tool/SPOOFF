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

#include "fix_spooff_setconvection.h"

#include "fix_spooff_tlsph_reference_configuration.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "update.h"
#include "variable.h"

#include <cstring>
#include <Eigen/Eigen>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NONE, CONSTANT, EQUAL, ATOM };

/* ---------------------------------------------------------------------- */

FixSPOOFFSetConvection::FixSPOOFFSetConvection(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), estr(nullptr), idregion(nullptr),
    region(nullptr), sforce(nullptr)
{
  if (narg < 5) error->all(FLERR, "Illegal fix setconvection command 1");

  dynamic_group_allow = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;

  if (strstr(arg[3], "v_") == arg[3]) {
    estr = utils::strdup(&arg[3][2]);
  } else if (strcmp(arg[3], "NULL") == 0) {
    estyle = NONE;
  } else {
    hvalue = utils::numeric(FLERR, arg[3], false, lmp);
    T_inf = utils::numeric(FLERR, arg[4], false, lmp);
    estyle = CONSTANT;
  }

  fix_tlsph_reference_configuration = nullptr;
 


  // optional args

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "region") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix setconvection command 2");
      region = domain->get_region_by_id(arg[iarg + 1]);
      if (!region) error->all(FLERR, "Region {} for fix setconvection does not exist", arg[iarg + 1]);
      idregion = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix setconvection command 3");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;

  maxatom = atom->nmax;
  memory->create(sforce, maxatom, 3, "setconvection:sforce");
}

/* ---------------------------------------------------------------------- */

FixSPOOFFSetConvection::~FixSPOOFFSetConvection()
{
  delete[] estr;
  delete[] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixSPOOFFSetConvection::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPOOFFSetConvection::init()
{
  // check variables

  if (estr) {
    evar = input->variable->find(estr);
    if (evar < 0) error->all(FLERR, "Variable name for fix setconvection does not exist");
    if (input->variable->equalstyle(evar))
      estyle = EQUAL;
    else if (input->variable->atomstyle(evar))
      estyle = ATOM;
    else
      error->all(FLERR, "Variable for fix setconvection is invalid style");
  }
  

  // set index and check validity of region

  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix setconvection does not exist", idregion);
  }

  if (estyle == ATOM)
    varflag = ATOM;
  else if (estyle == EQUAL)
    varflag = EQUAL;
  else
    varflag = CONSTANT;

  // cannot use non-zero forces for a minimization since no energy is integrated
  // use fix addforce instead

  int flag = 0;
  if (update->whichflag == 2) {
    if (estyle == EQUAL || estyle == ATOM) flag = 1;
    if (estyle == CONSTANT && hvalue != 0.0) flag = 1;
  }
  if (flag) error->all(FLERR, "Cannot use non-zero forces in an energy minimization");


  // find associated SPOOFF_TLSPH_NEIGHBORS fix that must exist
// could have changed locations in fix list since created

  int i;
  ifix_tlsph = -1;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style, "SPOOFF_TLSPH_NEIGHBORS") == 0)
      ifix_tlsph = i;
  if (ifix_tlsph == -1)
    error->all(FLERR, "Fix SPOOFF_TLSPH_NEIGHBORS does not exist");

}

/* ---------------------------------------------------------------------- */

void FixSPOOFFSetConvection::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else
    error->all(FLERR, "Fix spooff/setconvection does not support RESPA");
}

/* ---------------------------------------------------------------------- */

void FixSPOOFFSetConvection::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSPOOFFSetConvection::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double *tsph = atom->tsph;
  double *ksph = atom->ksph;
  double **qsph = atom->qsph;
  double *desph = atom->desph;
  double *esph = atom->esph;
  double *cpsph = atom->cpsph;
  double *rmass = atom->rmass;
  double **vest = atom->vest;
  double *vfrac = atom->vfrac;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  tagint *mol = atom->molecule;
  int jnum, jj, i, j, jnearest, idim, first;
  double nearestT, nearestr, r, fh, length;
  //double ClosestR[nlocal] ={0};
  //double ClosestT[nlocal] ={0};
  //int First[nlocal] ={true};


  tagint **partner = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->partner;
  int *npartner = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->npartner;
  Eigen::Vector3d xi, xj, dx;

  // update region if necessary

  if (region) region->prematch();


  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
        if (estyle) {
          //need to first collect nearest system neighbour information (assumes convective boundary is at edge of simulation)
          jnum = npartner[i];
          nearestT = 0;
          nearestr= 0;
          jnearest = 0;
          first = true;
          for (jj = 0; jj < jnum; jj++) {

        if (partner[i][jj] == 0)
          continue;
        j = atom->map(partner[i][jj]);
        if (j < 0) { //                 // check if lost a partner without first breaking bond
          partner[i][jj] = 0;
          continue;
        }

        if (mol[j] < 0) { // particle has failed. do not include it for computing any property
          continue;
        }

        if (mol[i] != mol[j]) {
          continue;
        }

        if (region && region->match(x[j][0], x[j][1], x[j][2])){ 
        //ignore partners in convection region  (should update this to be just particles in given boundary region)
        }else{
        for (idim = 0; idim < 3; idim++) {
          xi(idim) = x[i][idim];
          xj(idim) = x[j][idim];
        }

        dx = xj - xi;
        r = dx.norm(); // current distance

        if(first==true){
          nearestT = tsph[j];
          nearestr = r;
          jnearest = j;
          first=false;
        }else{
          if(r<nearestr){
            nearestT = tsph[j];
            nearestr = r;
            jnearest = j;
          }
        }


      }}   
      if(first==false){
        if(hvalue==0){
          idim = domain->dimension;
          length = pow(vfrac[i],(1)/idim);
          tsph[i] = nearestT; 
          esph[i] = nearestT*cpsph[i]*rmass[i];
        }else{
          fh = (-hvalue*nearestr/ksph[i]);
          //tsph[i] = (nearestT*(1-(hvalue/ksph[i])*(nearestr-(length/2))) + T_inf*(nearestr*hvalue/ksph[i]))/(1+(length*hvalue/ksph[i])*0.5);
          //esph[i] = tsph[i]*cpsph[i]*rmass[i];
          tsph[i] = ((nearestT - T_inf*fh)/(1-fh));
          esph[i] = ((nearestT - T_inf*fh)/(1-fh))*cpsph[i]*rmass[i];
          //printf("------------------------------------------------------ \n");
          //printf("Boudary %i, x %f, y %f \n", i, x[i][0], x[i][1]);
          //printf("Nearest %i, x %f, y %f \n", jnearest, x[jnearest][0], x[jnearest][1]);
          //printf("Nearest Temperature %f, ", nearestT);
          //printf("Nearest Distance %f, ", nearestr);
          //printf("fh %f, ", fh);
          //printf("New Temperature %f \n ", tsph[i]);
        }
      }else{
          tsph[i] = T_inf;
          esph[i] = T_inf*cpsph[i]*rmass[i];
          //printf("------------------------------------------------------ \n");
          //printf("ID %i, ", i);
          //printf("No pair found %f, ", nearestT);
          //printf("New Temperature %f \n ", tsph[i]);
          //desph[i] = -1*(hvalue*(tsph[i]-T_inf) -qsph[i][0])*cpsph[i]*rmass[i];  // old style just uses heatflux which isn't good enough and makes boundary geometry hard (needs updateing to handel heatflux direction)
      }
        }
      }

  } 
}


/* ----------------------------------------------------------------------
 return components of total force on fix group before force was changed
 ------------------------------------------------------------------------- */

double FixSPOOFFSetConvection::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal, foriginal_all, 3, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return foriginal_all[n];
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double FixSPOOFFSetConvection::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax * 3 * sizeof(double);
  return bytes;
}
