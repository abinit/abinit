#!/bin/bash -x
export OMP_NUM_THREADS=1
mpirun -n 1 time abinit < t03_x.files > log_CB1_1node_1120bands 2> err_CB1_1node_1120bands
