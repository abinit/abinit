#!/bin/bash -x
export OMP_NUM_THREADS=1
mpirun -n 2 time abinit < t03_2_x.files > log_CB1_2node_1120bands 2> err_CB1_2node_1120bands


