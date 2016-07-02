#!/usr/bin/env bash
# This script looks for source files not present in object lists

[ -s configure.ac ] || exit 1

echo ""
echo "Looking for OpenMP directives"
echo "-----------------------------"
echo ""

for d in src/[0-9]*
do
  omp_files=""

  echo -n "Directory $d:"

  cd ${d} > /dev/null
  for f in `ls *.f *.F90 2> /dev/null`
  do
    grep '!$OMP' ${f} > /dev/null && omp_files="${omp_files} ${f}"
  done
  cd - > /dev/null 2>&1

  echo "${omp_files:- none}"
done
