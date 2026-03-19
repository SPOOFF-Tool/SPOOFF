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

#include "fix_spooff_setconductivity.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NONE, CONSTANT, EQUAL, ATOM };

/* ---------------------------------------------------------------------- */

FixSPOOFFSetConductivity::FixSPOOFFSetConductivity(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), estr(nullptr), idregion(nullptr),
    region(nullptr), sforce(nullptr)
{
  if (narg < 5) error->all(FLERR, "Illegal fix setconductivity command");

  dynamic_group_allow = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;

  if (strstr(arg[3], "v_") == arg[3]) {
    estr = utils::strdup(&arg[3][2]);
    ktype = utils::numeric(FLERR, arg[4], false, lmp);  //different conductivity types are currently: 0 constant, 1 is UO2 correlation, 2 zirrconium correlation UNITS ARE W/mK - use the scaling factor to adjust the units
  } else if (strcmp(arg[3], "NULL") == 0) {
    estyle = NONE;
  } else {
    scale = utils::numeric(FLERR, arg[3], false, lmp);
    ktype = utils::numeric(FLERR, arg[4], false, lmp);  //different conductivity types are currently: 0 constant, 1 is UO2 correlation, 2 zirrconium correlation
    estyle = CONSTANT;
  }



  // optional args

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "region") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix setconductivity command");
      region = domain->get_region_by_id(arg[iarg + 1]);
      if (!region) error->all(FLERR, "Region {} for fix setconductivity does not exist", arg[iarg + 1]);
      idregion = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix setconductivity command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;

  maxatom = atom->nmax;
  memory->create(sforce, maxatom, 3, "setconductivity:sforce");
}

/* ---------------------------------------------------------------------- */

FixSPOOFFSetConductivity::~FixSPOOFFSetConductivity()
{
  delete[] estr;
  delete[] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixSPOOFFSetConductivity::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPOOFFSetConductivity::init()
{
  // check variables

  if (estr) {
    evar = input->variable->find(estr);
    if (evar < 0) error->all(FLERR, "Variable name for fix setconductivity does not exist");
    if (input->variable->equalstyle(evar))
      estyle = EQUAL;
    else if (input->variable->atomstyle(evar))
      estyle = ATOM;
    else
      error->all(FLERR, "Variable for fix setconductivity is invalid style");
  }
  

  // set index and check validity of region

  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix setconductivity does not exist", idregion);
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
    if (estyle == CONSTANT && scale != 0.0) flag = 1;
  }
  if (flag) error->all(FLERR, "Cannot use non-zero forces in an energy minimization");
}

/* ---------------------------------------------------------------------- */

void FixSPOOFFSetConductivity::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else
    error->all(FLERR, "Fix spooff/setconductivity does not support RESPA");
}

/* ---------------------------------------------------------------------- */

void FixSPOOFFSetConductivity::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSPOOFFSetConductivity::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  double *ksph = atom->ksph;
  double *tsph = atom->tsph;
  double **vest = atom->vest;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double porosity = 1; //hardcoded for now

  // update region if necessary

  if (region) region->prematch();

  // reallocate sforce array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce, maxatom, 3, "setconductivity:sforce");
  }


  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
        if (estyle) {
          if(ktype==0){
            ksph[i] = scale;
          }else if(ktype==1){
            if(tsph[i]<(1650+273)){
              ksph[i] = scale*porosity*100*((40.4/(464+(tsph[i]-273))) + 1.216e-4*exp(1.867e-3*(tsph[i]-273)));  //100 factor to change to W/mK standard
            }else if(tsph[i]<(2840+273)){
              ksph[i] = scale*porosity*100*(0.0191 + 1.216e-4*exp(1.867e-3*(tsph[i]-273)));;  //100 factor to change to W/mK standard
            }else{
              ksph[i] = scale*porosity*0.04351656;
            }
          }else if(ktype==2){
            ksph[i] = scale*(7.51 + 2.09e-2*(tsph[i]-273) -1.45e-5*pow((tsph[i]-273),2) + 7.67e-9*pow((tsph[i]-273),3));
          }          
        }
      }

    // variable force, wrap with clear/add

  } else {

    modify->clearstep_compute();

    if (estyle == EQUAL)
      scale = input->variable->compute_equal(evar);
    else if (estyle == ATOM)
      input->variable->compute_atom(evar, igroup, &sforce[0][0], 3, 0);
    

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
        if (estyle == ATOM) {
          if(ktype==0){
            ksph[i] = scale;
          }else if(ktype==1){
            if(tsph[i]<(1650+273)){
              ksph[i] = scale*porosity*100*((40.4/(464+(tsph[i]-273))) + 1.216e-4*exp(1.867e-3*(tsph[i]-273)));  //100 factor to change to W/mK standard
            }else if(tsph[i]<(2840+273)){
              ksph[i] = scale*porosity*100*(0.0191 + 1.216e-4*exp(1.867e-3*(tsph[i]-273)));;  //100 factor to change to W/mK standard
            }else{
              ksph[i] = scale*porosity*0.04351656;
            }
          }else if(ktype==2){
            ksph[i] = scale*(7.51 + 2.09e-2*(tsph[i]-273) -1.45e-5*pow((tsph[i]-273),2) + 7.67e-9*pow((tsph[i]-273),3));
          }   
        } else if (estyle) {
          if(ktype==0){
            ksph[i] = scale;
          }else if(ktype==1){
            if(tsph[i]<(1650+273)){
              ksph[i] = scale*porosity*100*((40.4/(464+(tsph[i]-273))) + 1.216e-4*exp(1.867e-3*(tsph[i]-273)));  //100 factor to change to W/mK standard
            }else if(tsph[i]<(2840+273)){
              ksph[i] = scale*porosity*100*(0.0191 + 1.216e-4*exp(1.867e-3*(tsph[i]-273)));;  //100 factor to change to W/mK standard
            }else{
              ksph[i] = scale*porosity*0.04351656;
            }
          }else if(ktype==2){
            ksph[i] = scale*(7.51 + 2.09e-2*(tsph[i]-273) -1.45e-5*pow((tsph[i]-273),2) + 7.67e-9*pow((tsph[i]-273),3));
          }     
        }
      }
  }
}

/* ----------------------------------------------------------------------
 return components of total force on fix group before force was changed
 ------------------------------------------------------------------------- */

double FixSPOOFFSetConductivity::compute_vector(int n)
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

double FixSPOOFFSetConductivity::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax * 3 * sizeof(double);
  return bytes;
}
