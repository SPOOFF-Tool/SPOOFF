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

#include "compute_spooff_tlsph_strain.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "update.h"

#include <cstring>
#include <Eigen/Eigen>          // IWYU pragma: export

using namespace Eigen;
using namespace std;
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSPOOFFTLSPHstrain::ComputeSPOOFFTLSPHstrain(LAMMPS *lmp, int narg, char **arg) :
                Compute(lmp, narg, arg) {
        if (narg != 3)
                error->all(FLERR, "Illegal compute spooff/tlsph_strain command");

        peratom_flag = 1;
        size_peratom_cols = 6;

        nmax = 0;
        strainVector = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSPOOFFTLSPHstrain::~ComputeSPOOFFTLSPHstrain() {
        memory->sfree(strainVector);
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFTLSPHstrain::init() {

        int count = 0;
        for (int i = 0; i < modify->ncompute; i++)
                if (strcmp(modify->compute[i]->style, "spooff/tlsph_strain") == 0)
                        count++;
        if (count > 1 && comm->me == 0)
                error->warning(FLERR, "More than one compute spooff/tlsph_strain");
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFTLSPHstrain::compute_peratom() {
        double **defgrad0 = atom->smd_data_9;

        invoked_peratom = update->ntimestep;

        // grow vector array if necessary

        if (atom->nmax > nmax) {
                memory->destroy(strainVector);
                nmax = atom->nmax;
                memory->create(strainVector, nmax, size_peratom_cols, "strainVector");
                array_atom = strainVector;
        }

        // copy data to output array
        int itmp = 0;
        auto Fincr = (Matrix3d *) force->pair->extract("spooff/tlsph/Fincr_ptr", itmp);
        if (Fincr == nullptr) {
                error->all(FLERR, "compute spooff/tlsph_strain failed to access Fincr array");
        }

        int *mask = atom->mask;
        int nlocal = atom->nlocal;
        Matrix3d E, eye, Ftotal, F0;
        eye.setIdentity();

        for (int i = 0; i < nlocal; i++) {
                if (mask[i] & groupbit) {

                        // old deformation gradient
                        F0(0, 0) = defgrad0[i][0];
                        F0(0, 1) = defgrad0[i][1];
                        F0(0, 2) = defgrad0[i][2];
                        F0(1, 0) = defgrad0[i][3];
                        F0(1, 1) = defgrad0[i][4];
                        F0(1, 2) = defgrad0[i][5];
                        F0(2, 0) = defgrad0[i][6];
                        F0(2, 1) = defgrad0[i][7];
                        F0(2, 2) = defgrad0[i][8];

                        // compute current total deformation gradient
                        Ftotal = F0 * Fincr[i]; // this is the total deformation gradient: reference deformation times incremental deformation


                        E = 0.5 * (Ftotal.transpose() * Ftotal - eye); // Green-Lagrange strain
                        strainVector[i][0] = E(0, 0);
                        strainVector[i][1] = E(1, 1);
                        strainVector[i][2] = E(2, 2);
                        strainVector[i][3] = E(0, 1);
                        strainVector[i][4] = E(0, 2);
                        strainVector[i][5] = E(1, 2);
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

double ComputeSPOOFFTLSPHstrain::memory_usage() {
        double bytes = (double)size_peratom_cols * nmax * sizeof(double);
        return bytes;
}
