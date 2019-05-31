#!/bin/bash -x
(cd ../../../../tmp/; make mj4)
mpirun -n 1 time abinit < t03_x.files > log 2> err

