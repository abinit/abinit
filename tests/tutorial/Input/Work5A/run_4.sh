#!/bin/bash -x
(cd ../../../../tmp/; make mj4)
export OMP_NUM_THREADS=3
mpirun -n 4 time abinit < t03_4_x.files  > log 2> err
#mpirun -n 4 time abinit < t03_4_x.files > log 2> err  | tail -f log  | grep -e ETOT 



