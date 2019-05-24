#!/bin/bash -x
mpirun -n 1 time abinit < t03_x.files > log 2> err
mpirun -n 2 time abinit < t03_x.files > log 2> err

