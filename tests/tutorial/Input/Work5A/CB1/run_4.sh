#!/bin/bash -x
export OMP_NUM_THREADS=1
mpirun -n 4 time abinit < t03_4_x.files > log_CB1_4node_1120bands 2> err_CB1_4node_1120bands




