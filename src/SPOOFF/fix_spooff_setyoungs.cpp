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

#include "fix_spooff_setyoungs.h"

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

FixSPOOFFSetYoungs::FixSPOOFFSetYoungs(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), estr(nullptr), idregion(nullptr),
    region(nullptr), sforce(nullptr)
{
  if (narg < 5) error->all(FLERR, "Illegal fix setyoungs command");

  dynamic_group_allow = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;

  if (strstr(arg[3], "v_") == arg[3]) {
    estr = utils::strdup(&arg[3][2]);
    ktype = utils::numeric(FLERR, arg[4], false, lmp);  //different youngs types are currently: 0 constant, 1 is UO2 correlation, 2 zirrconium beta phase correlation UNITS ARE N/m2 - use the scaling factor to adjust the units
  } else if (strcmp(arg[3], "NULL") == 0) {
    estyle = NONE;
  } else {
    scale = utils::numeric(FLERR, arg[3], false, lmp);
    ktype = utils::numeric(FLERR, arg[4], false, lmp);  //different youngs types are currently: 0 constant, 1 is UO2 correlation, 2 zirrconium beta phase correlation
    estyle = CONSTANT;
  }



  // optional args

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "region") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix setyoungs command");
      region = domain->get_region_by_id(arg[iarg + 1]);
      if (!region) error->all(FLERR, "Region {} for fix setyoungs does not exist", arg[iarg + 1]);
      idregion = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix setyoungs command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;

  maxatom = atom->nmax;
  memory->create(sforce, maxatom, 3, "setyoungs:sforce");
}

/* ---------------------------------------------------------------------- */

FixSPOOFFSetYoungs::~FixSPOOFFSetYoungs()
{
  delete[] estr;
  delete[] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixSPOOFFSetYoungs::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSPOOFFSetYoungs::init()
{
  // check variables

  if (estr) {
    evar = input->variable->find(estr);
    if (evar < 0) error->all(FLERR, "Variable name for fix setyoungs does not exist");
    if (input->variable->equalstyle(evar))
      estyle = EQUAL;
    else if (input->variable->atomstyle(evar))
      estyle = ATOM;
    else
      error->all(FLERR, "Variable for fix setyoungs is invalid style");
  }
  

  // set index and check validity of region

  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix setyoungs does not exist", idregion);
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

void FixSPOOFFSetYoungs::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else
    error->all(FLERR, "Fix spooff/setyoungs does not support RESPA");
}

/* ---------------------------------------------------------------------- */

void FixSPOOFFSetYoungs::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSPOOFFSetYoungs::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  double *youngs = atom->youngs;
  double *tsph = atom->tsph;
  double **vest = atom->vest;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *porosity = atom->porosity;

  // update region if necessary

  if (region) region->prematch();

  // reallocate sforce array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce, maxatom, 3, "setyoungs:sforce");
  }


  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
        if (estyle) {
          if(ktype==0){
            youngs[i] = scale;
          }else if(ktype==1){
              youngs[i] = scale*2.334E11*(1-(2.752*(1-porosity[i])))*(1-(1.0915E-4 * tsph[i]));  //UO2
          }else if(ktype==2){
            youngs[i] = scale*(9.21E11 - (4.05E7*tsph[i])); //Zircalloy Beta Phase
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
            youngs[i] = scale;
          }else if(ktype==1){
              youngs[i] = scale*2.334E11*(1-(2.752*(1-porosity[i])))*(1-(1.0915E-4 * tsph[i]));  //UO2
          }else if(ktype==2){
            youngs[i] = scale*(9.21E11 - (4.05E7*tsph[i])); //Zircalloy Beta Phase
          }   
        } else if (estyle) {
          if(ktype==0){
            youngs[i] = scale;
          }else if(ktype==1){
              youngs[i] = scale*2.334E11*(1-(2.752*(1-porosity[i])))*(1-(1.0915E-4 * tsph[i]));  //UO2
          }else if(ktype==2){
            youngs[i] = scale*(9.21E11 - (4.05E7*tsph[i])); //Zircalloy Beta Phase
          }    
        }
      }
  }
}

/* ----------------------------------------------------------------------
 return components of total force on fix group before force was changed
 ------------------------------------------------------------------------- */

double FixSPOOFFSetYoungs::compute_vector(int n)
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

double FixSPOOFFSetYoungs::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax * 3 * sizeof(double);
  return bytes;
}
