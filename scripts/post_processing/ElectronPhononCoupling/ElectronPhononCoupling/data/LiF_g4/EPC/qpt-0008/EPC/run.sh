#!/bin/bash


MPIRUN='mpirun -n 4'
ABINIT='abinit'

ln -nfs ../../../../Den/out_data/odat_DEN input_data/idat_DEN
ln -nfs ../../../WFK/out_data/odat_WFK input_data/idat_WFK
ln -nfs ../../WFQ/out_data/odat_WFQ input_data/idat_WFQ
ln -nfs ../../../../DVSCF/qpt-0008/DVSCF/out_data/odat_DEN1 input_data/idat_DEN1
ln -nfs ../../../../DVSCF/qpt-0008/DVSCF/out_data/odat_DEN2 input_data/idat_DEN2
ln -nfs ../../../../DVSCF/qpt-0008/DVSCF/out_data/odat_DEN3 input_data/idat_DEN3
ln -nfs ../../../../DVSCF/qpt-0008/DVSCF/out_data/odat_DEN4 input_data/idat_DEN4
ln -nfs ../../../../DVSCF/qpt-0008/DVSCF/out_data/odat_DEN5 input_data/idat_DEN5
ln -nfs ../../../../DVSCF/qpt-0008/DVSCF/out_data/odat_DEN6 input_data/idat_DEN6

$MPIRUN $ABINIT < calc.files &> calc.log 2> stderr

