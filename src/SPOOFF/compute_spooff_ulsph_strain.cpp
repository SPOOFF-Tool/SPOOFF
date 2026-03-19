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

#include "compute_spooff_ulsph_strain.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace std;
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSPOOFFULSPHstrain::ComputeSPOOFFULSPHstrain(LAMMPS *lmp, int narg, char **arg) :
                Compute(lmp, narg, arg) {
        if (narg != 3)
                error->all(FLERR, "Illegal compute spooff/tlsph_strain command");

        peratom_flag = 1;
        size_peratom_cols = 6;

        nmax = 0;
        strainVector = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSPOOFFULSPHstrain::~ComputeSPOOFFULSPHstrain() {
        memory->sfree(strainVector);
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFULSPHstrain::init() {

        int count = 0;
        for (int i = 0; i < modify->ncompute; i++)
                if (strcmp(modify->compute[i]->style, "spooff/tlsph_strain") == 0)
                        count++;
        if (count > 1 && comm->me == 0)
                error->warning(FLERR, "More than one compute spooff/tlsph_strain");
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFULSPHstrain::compute_peratom() {
        double **atom_data9 = atom->smd_data_9; // ULSPH strain is stored in the first 6 entries of this data field

        invoked_peratom = update->ntimestep;

        // grow vector array if necessary

        if (atom->nmax > nmax) {
                memory->destroy(strainVector);
                nmax = atom->nmax;
                memory->create(strainVector, nmax, size_peratom_cols, "strainVector");
                array_atom = strainVector;
        }

        int *mask = atom->mask;
        int nlocal = atom->nlocal;

        for (int i = 0; i < nlocal; i++) {
                if (mask[i] & groupbit) {

                        strainVector[i][0] = atom_data9[i][0];
                        strainVector[i][1] = atom_data9[i][1];
                        strainVector[i][2] = atom_data9[i][2];
                        strainVector[i][3] = atom_data9[i][3];
                        strainVector[i][4] = atom_data9[i][4];
                        strainVector[i][5] = atom_data9[i][5];
                } else {
                        for (int j = 0; j < size_peratom_cols; j++) {
                                strainVector[i][j] = 0.0;
                        }
                }
        }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSPOOFFULSPHstrain::memory_usage() {
        double bytes = (double)size_peratom_cols * nmax * sizeof(double);
        return bytes;
}
