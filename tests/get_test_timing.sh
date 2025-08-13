#! /bin/bash

for fil in */Refs/*
do
 timeline=`tail  $fil | grep 'Proc.   0 individual time'`
 optdrivline=`awk '/optdriver[0-9]*      [ ]*[1-9]/ {print $0, " "}' $fil`
 libline=`awk '/libxc|wannier90|psml|libxml|xmlf90/ {print $0, " "}' $fil`
 execline=`awk '/mrgddb|anaddb|cut3d|mrgscr|multibinit|optic/ {print $0, " "}' $fil`
 echo $fil $timeline $optdrivline $libline $execline
done
