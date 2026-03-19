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

#include "fix_spooff_tlsph_reference_configuration.h"

#include "compute_spooff_damage.h"
#include <cstring>
#include "atom.h"
#include "fix.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "pair.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

ComputeSPOOFFDamage::ComputeSPOOFFDamage(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute spooff/damage command");
  if (atom->damage_flag != 1) error->all(FLERR,"compute spooff/damage command requires atom_style with damage (e.g. spooff)");

  peratom_flag = 1;
  size_peratom_cols = 4;

  nmax = 0;
  damage_vector = nullptr;
  fix_tlsph_reference_configuration = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSPOOFFDamage::~ComputeSPOOFFDamage()
{
  memory->sfree(damage_vector);
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFDamage::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"spooff/damage") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute spooff/damage");

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

void ComputeSPOOFFDamage::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow Vector array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(damage_vector);
                nmax = atom->nmax;
                memory->create(damage_vector, nmax, size_peratom_cols, "atom:damage_vector");
                array_atom = damage_vector;
  }
  


  double *damage = atom->damage;  
  double *starting_neighs = atom->starting_neighs; 
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  tagint *mol = atom->molecule;
  int jnum, jj, i, j, idim;
  int itmp = 0;
  int *numNeighsRefConfig = (int *) force->pair->extract("spooff/tlsph/numNeighsRefConfig_ptr", itmp);
    if (numNeighsRefConfig == nullptr) {
        error->all(FLERR, "compute spooff/tlsph_num_neighs failed to access numNeighsRefConfig array");
    }

  tagint **partner = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->partner;
  int *npartner = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->npartner;
  float **degradation_ij = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->degradation_ij;
  float **energy_per_bond = (dynamic_cast<FixSPOOFF_TLSPH_ReferenceConfiguration *>(modify->fix[ifix_tlsph]))->energy_per_bond;
  int lostNeighs;


    for (int i = 0; i < nlocal; i++) {
      jnum = npartner[i];
      //floating point error protection for first time step
      if(starting_neighs[i] == 0){
        starting_neighs[i]=1.0;
      }
      if(numNeighsRefConfig[i] == 0){
        numNeighsRefConfig[i]=1.0;
      }
      lostNeighs = starting_neighs[i]-numNeighsRefConfig[i];
      if (mask[i] & groupbit) {
              damage_vector[i][0] = 1.0 - (numNeighsRefConfig[i]/starting_neighs[i]);  //most USEFUL measure (just percentage of neighbours broken)              
              damage_vector[i][1] = 0.0;              //used to track average degredation of all bonds - this is not that useful to asses material damage but very useful for model development as are the next 2
              damage_vector[i][2] = 0.0;              //used to tracke average energy per bond REMAINING unbroken
              damage_vector[i][3] = damage[i];        //used to track associated dame value for each atom 
              
              //need to first collect system neighbour information to average over it
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
        damage_vector[i][1] += MIN(degradation_ij[i][jj],1.0); //degradation over 1 casues failure anyway and gives a better understanding of 'average' if we limit it to 1
        damage_vector[i][2] += energy_per_bond[i][jj];
      
    }
    damage_vector[i][1] = damage_vector[i][1] + lostNeighs*1.0;
    damage_vector[i][1] = damage_vector[i][1]/starting_neighs[i];

    damage_vector[i][2] = damage_vector[i][2]/numNeighsRefConfig[i];
      }
      else {
              damage_vector[i][0] = 0.0;
              damage_vector[i][1] = 0.0;
              damage_vector[i][2] = 0.0;
              damage_vector[i][3] = 0.0;
      }
      
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeSPOOFFDamage::memory_usage()
{
  double bytes = (double)size_peratom_cols * nmax * sizeof(double);
  return bytes;
}
