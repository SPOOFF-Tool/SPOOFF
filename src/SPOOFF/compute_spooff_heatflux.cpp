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
#include "compute_spooff_heatflux.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSPOOFFHeatflux::ComputeSPOOFFHeatflux(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute spooff/heatflux command");
  if (atom->esph_flag != 1) error->all(FLERR,"compute spooff/heatflux command requires atom_style with internal_energy (e.g. spooff)");

  peratom_flag = 1;
  size_peratom_cols = 3;


  nmax = 0;
  heatflux_array = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSPOOFFHeatflux::~ComputeSPOOFFHeatflux()
{
  memory->sfree(heatflux_array);
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFHeatflux::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"spooff/heatflux") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute spooff/heatflux");
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFHeatflux::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow rhoVector array if necessary

  if (atom->nmax > nmax) {
                memory->destroy(heatflux_array);
                nmax = atom->nmax;
                memory->create(heatflux_array, nmax, size_peratom_cols, "atom:heatflux_array");
                array_atom = heatflux_array;
  }


  double **qsph = atom->qsph;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
              heatflux_array[i][0] = qsph[i][0];
              heatflux_array[i][1] = qsph[i][1];
              heatflux_array[i][2] = qsph[i][2];
      }
      else {
              heatflux_array[i][0] = 0.0;
              heatflux_array[i][1] = 0.0;
              heatflux_array[i][2] = 0.0;
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeSPOOFFHeatflux::memory_usage()
{
  double bytes = (double)size_peratom_cols * nmax * sizeof(double);
  return bytes;

}
