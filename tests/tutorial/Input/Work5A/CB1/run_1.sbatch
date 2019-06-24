#!/bin/bash -x
#MSUB -r test_abinit                # Request name
#MSUB -n 1                         # Total number of tasks to use
#MSUB -c 1                        # Number of threads per task to use
#MSUB -N 1                         # Number of nodes
#MSUB -T 1800                      # Elapsed time limit in seconds
#MSUB -o abinit_%I.o              # Standard output. %I is the job id
#MSUB -e abinit_%I.e              # Error output. %I is the job id
##MSUB -q haswell                   # Choosing nehalem nodes
#MSUB -q sklb                   # Choosing nehalem nodes

export OMP_NUM_THREADS=1
export ABI_HOME=$WORKDIR/_apps/abinit
export PATH=$ABI_HOME/tmp/src/98_main/:$PATH
export KMP_AFFINITY=compact

ccc_mprun -n 1 time abinit < t03_x.files #> log_CB1_1node_1120bands 2> err_CB1_1node_1120bands
#mpirun -n 1 time abinit < t03_x.files > log_CB1_1node_1120bands 2> err_CB1_1node_1120bands
