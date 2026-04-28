#!/bin/bash

cd Heat\=40Wperg/
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd Heat\=50Wperg/
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd Heat\=60Wperg/
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd Heat\=70Wperg/
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd Heat\=80Wperg/
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd Heat\=90Wperg/
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd Heat\=100Wperg/
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd Heat\=110Wperg/
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd Heat\=120Wperg/
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'

