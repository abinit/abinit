#!/bin/bash


MPIRUN='mpirun -n 4'
ABINIT='abinit'

ln -nfs ../../../../Den/out_data/odat_DEN input_data/idat_DEN

$MPIRUN $ABINIT < calc.files &> calc.log 2> stderr

