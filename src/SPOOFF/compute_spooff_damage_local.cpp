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

#include "fix_spooff_tlsph_reference_configuration.h"
#include "compute_spooff_damage_local.h"

#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"
#include "modify.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NONE, NEIGH, PAIR, BOND, ANGLE, DIHEDRAL, IMPROPER };
enum { TYPE, RADIUS };

#define DELTA 10000

/* ---------------------------------------------------------------------- */

ComputeDamageLocal::ComputeDamageLocal(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), vlocal(nullptr), alocal(nullptr), indices(nullptr),
    pack_choice(nullptr)
{
  if (narg < 4) error->all(FLERR, "Illegal compute damage/local command");

  local_flag = 1;
  nvalues = narg - 3;
  pack_choice = new FnPtrPack[nvalues];

  kindflag = NONE;

  int i;
  nvalues = 0;
  int iarg = 3;
  while (iarg < narg) {
    i = iarg - 3;


    if (strcmp(arg[iarg], "patom1") == 0) {
      pack_choice[i] = &ComputeDamageLocal::pack_patom1;
      if (kindflag != NONE && kindflag != PAIR)
        error->all(FLERR, "Compute damage/local cannot use these inputs together");
      kindflag = PAIR;
    } else if (strcmp(arg[iarg], "patom2") == 0) {
      pack_choice[i] = &ComputeDamageLocal::pack_patom2;
      if (kindflag != NONE && kindflag != PAIR)
        error->all(FLERR, "Compute damage/local cannot use these inputs together");
      kindflag = PAIR;
    } else if (strcmp(arg[iarg], "damage") == 0) {
      pack_choice[i] = &ComputeDamageLocal::pack_damage;
      if (kindflag != NONE && kindflag != PAIR)
        error->all(FLERR, "Compute damage/local cannot use these inputs together");
      kindflag = PAIR;
      } else if (strcmp(arg[iarg], "degredation") == 0) {
      pack_choice[i] = &ComputeDamageLocal::pack_degredation;
      if (kindflag != NONE && kindflag != PAIR)
        error->all(FLERR, "Compute damage/local cannot use these inputs together");
      kindflag = PAIR;
      } else if (strcmp(arg[iarg], "energy") == 0) {
      pack_choice[i] = &ComputeDamageLocal::pack_energy;
      if (kindflag != NONE && kindflag != PAIR)
        error->all(FLERR, "Compute damage/local cannot use these inputs together");
      kindflag = PAIR;

    } else
      break;

    nvalues++;
    iarg++;
  }

  if (nvalues == 1)
    size_local_cols = 0;
  else
    size_local_cols = nvalues;

  // optional args

  cutstyle = TYPE;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "cutoff") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal compute damage/local command");
      if (strcmp(arg[iarg + 1], "type") == 0)
        cutstyle = TYPE;
      else if (strcmp(arg[iarg + 1], "radius") == 0)
        cutstyle = RADIUS;
      else
        error->all(FLERR, "Illegal compute damage/local command");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal compute damage/local command");
  }

  // error check

  if (atom->molecular == 2 &&
      (kindflag == BOND || kindflag == ANGLE || kindflag == DIHEDRAL || kindflag == IMPROPER))
    error->all(FLERR,
               "Compute damage/local does not (yet) work "
               "with atom_style template");

  if (kindflag == BOND && atom->avec->bonds_allow == 0)
    error->all(FLERR, "Compute damage/local for damage that isn't allocated");
  if (kindflag == ANGLE && atom->avec->angles_allow == 0)
    error->all(FLERR, "Compute damage/local for damage that isn't allocated");
  if (kindflag == DIHEDRAL && atom->avec->dihedrals_allow == 0)
    error->all(FLERR, "Compute damage/local for damage that isn't allocated");
  if (kindflag == IMPROPER && atom->avec->impropers_allow == 0)
    error->all(FLERR, "Compute damage/local for damage that isn't allocated");
  if (cutstyle == RADIUS && !atom->radius_flag)
    error->all(FLERR, "Compute damage/local requires atom attribute radius");

  nmax = 0;
  vlocal = nullptr;
  alocal = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeDamageLocal::~ComputeDamageLocal()
{
  delete[] pack_choice;
  memory->destroy(vlocal);
  memory->destroy(alocal);
  memory->destroy(indices);
}

/* ---------------------------------------------------------------------- */

void ComputeDamageLocal::init()
{
  if (kindflag == NEIGH || kindflag == PAIR) {
    if (force->pair == nullptr)
      error->all(FLERR, "No pair style is defined for compute damage/local");
    if (force->pair->single_enable == 0)
      error->all(FLERR, "Pair style does not support compute damage/local");
  }

  // for NEIGH/PAIR need an occasional half neighbor list
  // set size to same value as request made by force->pair
  // this should enable it to always be a copy list  (e.g. for granular pstyle)

  if (kindflag == NEIGH || kindflag == PAIR) {
    int neighflags = NeighConst::REQ_OCCASIONAL;
    auto pairrequest = neighbor->find_request(force->pair);
    if (pairrequest && pairrequest->get_size()) neighflags |= NeighConst::REQ_SIZE;
    neighbor->add_request(this, neighflags);
  }

  // do initial memory allocation so that memory_usage() is correct
  // cannot be done yet for NEIGH/PAIR, since neigh list does not exist

  if (kindflag == NEIGH)
    ncount = 0;
  else if (kindflag == PAIR)
    ncount = 0;


  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;

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

void ComputeDamageLocal::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeDamageLocal::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and generate list of indices

  if (kindflag == NEIGH)
    ncount = count_pairs(0, 0);
  else if (kindflag == PAIR)
    ncount = count_pairs(0, 1);


  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;

  if (kindflag == NEIGH)
    ncount = count_pairs(1, 0);
  else if (kindflag == PAIR)
    ncount = count_pairs(1, 1);


  // fill vector or array with local values

  if (nvalues == 1) {
    buf = vlocal;
    (this->*pack_choice[0])(0);
  } else {
    if (alocal) buf = &alocal[0][0];
    for (int n = 0; n < nvalues; n++) (this->*pack_choice[n])(n);
  }
}

/* ----------------------------------------------------------------------
   count pairs and compute pair info on this proc
   only count pair once if newton_pair is off
   both atom I,J must be in group
   if allflag is set, compute requested info about pair
   if forceflag = 1, pair must be within force cutoff, else neighbor cutoff
------------------------------------------------------------------------- */

int ComputeDamageLocal::count_pairs(int allflag, int forceflag)
{
  int i, j, m, ii, jj, inum, jnum, itype, jtype;
  tagint itag, jtag;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, radsum;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double **x = atom->x;
  tagint *mol = atom->molecule;
  double *radius = atom->radius;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  tagint **partnerfirst = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->partnerfirst;
  int *npartner = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->npartner;

  // invoke half neighbor list (will copy or build if necessary)  

  inum = list->inum;

  // loop over neighbors of my atoms
  // skip if I or J are not in group
  // for newton = 0 and J = ghost atom,
  //   need to ensure I,J pair is only output by one proc
  //   use same itag,jtag logic as in Neighbor::neigh_half_nsq()

  double **cutsq = force->pair->cutsq;

  m = 0;
  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itag = tag[i];
    itype = type[i];
    jnum = npartner[i];

    for (jj = 0; jj < jnum; jj++) {
      if (partnerfirst[i][jj] == 0)
        continue;
      j = atom->map(partnerfirst[i][jj]);
      if (mol[i] < 0 | mol[j] < 0) 
        continue;
      //j &= NEIGHMASK;

      if (!(mask[j] & groupbit)) continue;

      // itag = jtag is possible for long cutoffs that include images of self

      if (j >= nlocal) {
        jtag = partnerfirst[i][jj];
        if (itag > jtag) {
          if ((itag + jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag + jtag) % 2 == 1) continue;
        } else {
          if (x[j][2] < ztmp) continue;
          if (x[j][2] == ztmp) {
            if (x[j][1] < ytmp) continue;
            if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
          }
        }
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      if (forceflag) {
        if (cutstyle == TYPE) {
          if (rsq >= cutsq[itype][jtype]) continue;
        } else {
          radsum = radius[i] + radius[j];
          if (rsq >= radsum * radsum) continue;
        }
      }

      if (allflag) {
        indices[m][0] = i;
        indices[m][1] = j;  // this is the actual indicies j
        indices[m][2] = jj;   // this is the array indicies for the reference configuration
      }
      m++;
    }
  }

  return m;
}








/* ---------------------------------------------------------------------- */

void ComputeDamageLocal::reallocate(int n)
{
  // grow vector_local or array_local, also indices

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vlocal);
    memory->create(vlocal, nmax, "damage/local:vector_local");
    vector_local = vlocal;
  } else {
    memory->destroy(alocal);
    memory->create(alocal, nmax, nvalues, "damage/local:array_local");
    array_local = alocal;
  }

  memory->destroy(indices);
  memory->create(indices, nmax, 3, "damage/local:indices");
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeDamageLocal::memory_usage()
{
  double bytes = (double) nmax * nvalues * sizeof(double);
  bytes += (double) nmax * 2 * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword compute damage/local can output
   the atom damage is packed into buf starting at n with stride nvalues
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void ComputeDamageLocal::pack_patom1(int n)
{
  int i;
  tagint *tag = atom->tag;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    buf[n] = tag[i];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeDamageLocal::pack_patom2(int n)
{
  int i, j, jj;
  tagint *tag = atom->tag;
  tagint **partnerfirst = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->partnerfirst;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    jj = indices[m][2];
    buf[n] = partnerfirst[i][jj];
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeDamageLocal::pack_damage(int n)
{

  int i, j, jj;
  float **degradation_ij = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->degradation_ij;
  

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    jj = indices[m][2];
    buf[n] = MIN(degradation_ij[i][jj],1.0);
    n += nvalues;
  }
}

void ComputeDamageLocal::pack_degredation(int n)
{

  int i, j, jj;
  float **degradation_ij = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->degradation_ij;
  

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    jj = indices[m][2];
    buf[n] = degradation_ij[i][jj];
    n += nvalues;
  }
}

void ComputeDamageLocal::pack_energy(int n)
{

  int i, j, jj;
  float **energy_per_bond = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->energy_per_bond;

  for (int m = 0; m < ncount; m++) {
    i = indices[m][0];
    j = indices[m][1];
    jj = indices[m][2];
    buf[n] = energy_per_bond[i][jj];
    n += nvalues;
  }
}

