#!/bin/sh
#
# Copyright (C) 2005-2020 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

set -e

# Init
my_name="run-basic-tests"
my_cnffile="tests.env"

# Check arguments
if test "${#}" -lt "2"; then
 echo "Usage: ${my_name} test_dir test_number"
 exit 0
fi

# Check config file
if test -s "${my_cnffile}"; then
 . "${my_cnffile}"
else
 echo "${my_name}: config file ${my_cnffile} not found - aborting now."
 exit 1
fi

# Finish init
my_string="tmp-`hostname`_`date '+%Y%m%d'`"
test_dir="${1}"
test_number="${2}"

mkdir -p "${abinit_outdir}/${test_dir}/${my_string}"
my_output="${abinit_outdir}/${test_dir}/${my_string}/test${test_number}"

# Clean-up
rm -f "${my_output}.in" "${my_output}.out" "${my_output}.log"
rm -f "${my_output}i_*" "${my_output}o_*" "${my_output}_*"
cp "${abinit_inpdir}/${test_dir}/Input/test${test_number}.in" "${my_output}.in"
cp "${abinit_inpdir}/${test_dir}/Input/testin_wannier90o_w90.win" "${my_dir}/testin_wannier90o_w90.win"

# Write abinit.files for test
cat > "${my_output}.files" <<EOF
test${test_number}.in
test${test_number}.out
test${test_number}i
test${test_number}o
test${test_number}
EOF

case "${test_number}" in

 in_fast)
  echo "${abinit_pspdir}/01h.pspgth" >> "${my_output}.files"
  ;;

 in_v1)
  echo "${abinit_pspdir}/70yb.pspnc" >> "${my_output}.files"
  ;;

 in_v5)
  echo "${abinit_pspdir}/01h.pspgth" >> "${my_output}.files"
  echo "${abinit_pspdir}/04be.pspgth" >> "${my_output}.files"
  ;;

 in_bigdft)
  echo "${abinit_pspdir}/01h.pspgth" >> "${my_output}.files"
  ;;

 in_etsf_io)
  echo "${abinit_pspdir}/20ca.paw" >> "${my_output}.files"
  ;;

 in_libxc)
  echo "${abinit_pspdir}/83bi.psphgh" >> "${my_output}.files"
  ;;

 in_wannier90)
  echo "${abinit_pspdir}/31ga.pspnc" >> "${my_output}.files"
  echo "${abinit_pspdir}/33as.pspnc" >> "${my_output}.files"
  ;;

 *)
  echo "${my_name}: unknown test number ${test_number} - aborting now."
  rm -f "${my_output}.files"
  exit 2
  ;;

esac

cd "${abinit_outdir}/${test_dir}/${my_string}"
"${abinit_bindir}/abinit" < "${my_output}.files" > "${my_output}.log"

if test -s "test${test_number}_STATUS"; then
 cat "test${test_number}_STATUS"
fi
