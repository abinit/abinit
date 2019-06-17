#!/bin/bash -x
(cd ../../../../tmp/; make mj4)
#mpirun -n 1 time abinit < t03_x.files > log 2> err | tail -f log
#mpirun -n 1 time abinit < t03_x.files > log 2> err
export OMP_NUM_THREADS=12
mpirun -n 1 time abinit < t03_x.files > log 2> err | tail -f log  | grep -e ETOT 
