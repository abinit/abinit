#!/bin/bash -x
(cd ../../../../tmp/; make mj4)
#mpirun -n 4 time abinit < t03_4_x.files  > log 2> err
#set OMP_NUM_THREADS=4
#export OMP_NUM_THREADS=4
mpirun -n 4 -genv OMP_NUM_THREADS=4 -genv I_MPI_PIN_DOMAIN=omp time abinit < t03_4_x.files > log 2> err  | tail -f log  | grep -e ETOT 



