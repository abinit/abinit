#!/bin/bash


MPIRUN='mpirun -n 4'
ABINIT='abinit'

$MPIRUN $ABINIT < calc.files &> calc.log 2> stderr

