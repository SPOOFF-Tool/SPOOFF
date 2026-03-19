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

#include <cstring>
#include "compute_spooff_internal_energy.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSPOOFFInternalEnergy::ComputeSPOOFFInternalEnergy(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute spooff/internal_energy command");
  if (atom->esph_flag != 1) error->all(FLERR,"compute spooff/internal_energy command requires atom_style with internal_energy (e.g. spooff)");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  internal_energy_vector = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSPOOFFInternalEnergy::~ComputeSPOOFFInternalEnergy()
{
  memory->sfree(internal_energy_vector);
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFInternalEnergy::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"spooff/internal_energy") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute spooff/internal_energy");
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFInternalEnergy::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow rhoVector array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(internal_energy_vector);
    nmax = atom->nmax;
    internal_energy_vector = (double *) memory->smalloc(nmax*sizeof(double),"atom:internal_energy_vector");
    vector_atom = internal_energy_vector;
  }

  double *esph = atom->esph;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
              internal_energy_vector[i] = esph[i];
      }
      else {
              internal_energy_vector[i] = 0.0;
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeSPOOFFInternalEnergy::memory_usage()
{
  double bytes = (double)nmax * sizeof(double);
  return bytes;
}
