#!/bin/bash -x
(cd ../../../../tmp/; make mj4)
#mpirun -n 2 time abinit < t03_2_x.files > log 2> err
export OMP_NUM_THREADS=6
mpirun -n 2 time abinit < t03_2_x.files > log 2> err  | tail -f log  | grep -e ETOT 

