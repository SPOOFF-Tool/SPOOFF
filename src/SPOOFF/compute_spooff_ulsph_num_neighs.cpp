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
#include "compute_spooff_ulsph_num_neighs.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "pair.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSPOOFFULSPHNumNeighs::ComputeSPOOFFULSPHNumNeighs(LAMMPS *lmp, int narg, char **arg) :
        Compute(lmp, narg, arg) {
    if (narg != 3)
        error->all(FLERR, "Illegal compute spooff/ulsph_num_neighs command");

    peratom_flag = 1;
    size_peratom_cols = 0;

    nmax = 0;
    numNeighsOutput = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSPOOFFULSPHNumNeighs::~ComputeSPOOFFULSPHNumNeighs() {
    memory->destroy(numNeighsOutput);
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFULSPHNumNeighs::init() {
    int count = 0;
    for (int i = 0; i < modify->ncompute; i++)
        if (strcmp(modify->compute[i]->style, "spooff/ulsph_num_neighs") == 0)
            count++;
    if (count > 1 && comm->me == 0)
        error->warning(FLERR, "More than one compute spooff/ulsph_num_neighs");
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFULSPHNumNeighs::compute_peratom() {
    invoked_peratom = update->ntimestep;

    if (atom->nmax > nmax) {
        memory->destroy(numNeighsOutput);
        nmax = atom->nmax;
        memory->create(numNeighsOutput, nmax, "ulsph/num_neighs:numNeighsRefConfigOutput");
        vector_atom = numNeighsOutput;
    }

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    int itmp = 0;
    int *numNeighs = (int *) force->pair->extract("spooff/ulsph/numNeighs_ptr", itmp);
    if (numNeighs == nullptr) {
        error->all(FLERR, "compute spooff/ulsph_num_neighs failed to access numNeighs array");
    }

    for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
            numNeighsOutput[i] = numNeighs[i];
        } else {
            numNeighsOutput[i] = 0.0;
        }
    }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSPOOFFULSPHNumNeighs::memory_usage() {
    double bytes = (double)nmax * sizeof(double);
    return bytes;
}
