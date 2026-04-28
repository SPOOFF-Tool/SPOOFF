#!/bin/bash
echo "Hello, World"
cd dp\=2
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd dp\=1.5
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd dp\=1
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd dp\=0.8
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd dp\=0.5
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd dp\=0.4
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd dp\=0.3
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
cd dp\=0.2
mpirun -np 7 ../../../src/lmp_mpi -in 'run.lmp'
cd ..
echo "Done"