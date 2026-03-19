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

#include "compute_spooff_rho.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSPOOFFRho::ComputeSPOOFFRho(LAMMPS *lmp, int narg, char **arg) :
                Compute(lmp, narg, arg) {
        if (narg != 3)
                error->all(FLERR, "Illegal compute spooff/rho command");
        if (atom->vfrac_flag != 1)
                error->all(FLERR, "compute spooff/rho command requires atom_style with volume (e.g. spooff)");

        peratom_flag = 1;
        size_peratom_cols = 0;

        nmax = 0;
        rhoVector = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSPOOFFRho::~ComputeSPOOFFRho() {
        memory->sfree(rhoVector);
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFRho::init() {

        int count = 0;
        for (int i = 0; i < modify->ncompute; i++)
                if (strcmp(modify->compute[i]->style, "spooff/rho") == 0)
                        count++;
        if (count > 1 && comm->me == 0)
                error->warning(FLERR, "More than one compute spooff/rho");
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFRho::compute_peratom() {
        invoked_peratom = update->ntimestep;

        // grow rhoVector array if necessary

        if (atom->nmax > nmax) {
                memory->sfree(rhoVector);
                nmax = atom->nmax;
                rhoVector = (double *) memory->smalloc(nmax * sizeof(double), "atom:rhoVector");
                vector_atom = rhoVector;
        }

        double *rhosph = atom->rhosph;
        int *mask = atom->mask;
        int nlocal = atom->nlocal;

        for (int i = 0; i < nlocal; i++) {
                if (mask[i] & groupbit) {
                        rhoVector[i] = rhosph[i];
                } else {
                        rhoVector[i] = 0.0;
                }
        }

}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSPOOFFRho::memory_usage() {
        double bytes = (double)nmax * sizeof(double);
        return bytes;
}
