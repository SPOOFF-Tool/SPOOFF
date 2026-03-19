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

#include "compute_spooff_tlsph_stress_cylindrical.h"
#include <cmath>
#include <cstring>
#include <Eigen/Eigen>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "pair.h"

using namespace Eigen;
using namespace LAMMPS_NS;


/*
 * deviator of a tensor
 */
static Matrix3d Deviator(const Matrix3d& M) {
        Matrix3d eye;
        eye.setIdentity();
        eye *= M.trace() / 3.0;
        return M - eye;
}

/* ---------------------------------------------------------------------- */

ComputeSPOOFFTLSPHStressCylindrical::ComputeSPOOFFTLSPHStressCylindrical(LAMMPS *lmp, int narg, char **arg) :
                Compute(lmp, narg, arg) {
        if (narg != 3)
                error->all(FLERR, "Illegal compute spooff/tlsph_stress_cylindrical command");

        peratom_flag = 1;
        size_peratom_cols = 7;

        nmax = 0;
        stress_cylindrical_array = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSPOOFFTLSPHStressCylindrical::~ComputeSPOOFFTLSPHStressCylindrical() {
        memory->sfree(stress_cylindrical_array);
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFTLSPHStressCylindrical::init() {

        int count = 0;
        for (int i = 0; i < modify->ncompute; i++)
                if (strcmp(modify->compute[i]->style, "spooff/tlsph_stress_cylindrical") == 0)
                        count++;
        if (count > 1 && comm->me == 0)
                error->warning(FLERR, "More than one compute spooff/tlsph_stress_cylindrical");
}

/* ---------------------------------------------------------------------- */

void ComputeSPOOFFTLSPHStressCylindrical::compute_peratom() {
        invoked_peratom = update->ntimestep;
        Matrix3d stress_deviator;
        double von_mises_stress;
        double **x = atom->x;
        double theta;

        // grow vector array if necessary

        if (atom->nmax > nmax) {
                memory->destroy(stress_cylindrical_array);
                nmax = atom->nmax;
                memory->create(stress_cylindrical_array, nmax, size_peratom_cols, "stress_cylindricaltensorVector");
                array_atom = stress_cylindrical_array;
        }

        int itmp = 0;
        auto T = (Matrix3d *) force->pair->extract("spooff/tlsph/stressTensor_ptr", itmp);
        if (T == nullptr) {
                error->all(FLERR, "compute spooff/tlsph_stress_cylindrical could not access stress tensors. Are the matching pair styles present?");
        }
        int nlocal = atom->nlocal;
        int *mask = atom->mask;

        for (int i = 0; i < nlocal; i++) {
                if (mask[i] & groupbit) {
                        theta = atan2(x[i][1],x[i][0]);
                        stress_deviator = Deviator(T[i]);
                        von_mises_stress = sqrt(3. / 2.) * stress_deviator.norm();
                        stress_cylindrical_array[i][0] = T[i](0, 0)*pow(cos(theta),2) + T[i](1, 1)*pow(sin(theta),2) + 2*T[i](0, 1)*sin(theta)*cos(theta); // rr
                        stress_cylindrical_array[i][1] = T[i](0, 0)*pow(sin(theta),2) + T[i](1, 1)*pow(cos(theta),2) - 2*T[i](0, 1)*sin(theta)*cos(theta); // tt
                        stress_cylindrical_array[i][2] = T[i](2, 2); // zz
                        stress_cylindrical_array[i][3] = (T[i](1, 1)-T[i](0, 0))*sin(theta)*cos(theta) + T[i](0, 1)*(pow(cos(theta),2)-pow(sin(theta),2)); // rt
                        stress_cylindrical_array[i][4] = T[i](0, 2)*cos(theta) + T[i](1, 2)*sin(theta); // rz
                        stress_cylindrical_array[i][5] = -T[i](0, 2)*sin(theta) + T[i](1, 2)*cos(theta); // tz
                        stress_cylindrical_array[i][6] = von_mises_stress;
                } else {
                        for (int j = 0; j < size_peratom_cols; j++) {
                                stress_cylindrical_array[i][j] = 0.0;
                        }
                }
        }
}

/* ----------------------------------------------------------------------
 memory usage of local atom-based array
 ------------------------------------------------------------------------- */

double ComputeSPOOFFTLSPHStressCylindrical::memory_usage() {
        double bytes = (double)size_peratom_cols * nmax * sizeof(double);
        return bytes;
}
