# Script for testing abinit in charge (TC = Test Charge)
# Before running this script, read README in this directory.

# Copyright (C) 2000-2016 ABINIT group (XG,LSi)
# This file is distributed under the terms of the
# GNU General Public License, see ~abinit/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
#
# For the same input file, test the execution of 1 to 4
# concurrent executions, then analyze the results.

# Usage under unix:
# C-shell: ( RunTC machine_name [First_test] >& log_file )
# bash: RunTC machine_name [First_test] > log_file 2>&1 
#

#*****************************************************************************

set -e

# Init
my_name="run-charge-tests"
my_cnffile="../tests.env"
my_outdir=`hostname`

# Set-up environment
if test -s "${my_cnffile}"; then
 . ${my_cnffile}
else
 echo "${my_name}: ${my_cnffile} not found - aborting now"
 exit 1
fi

# Set codename for sequential version
CODE_SEQ="${abinit_bindir}/abinit"
# Set pseudopotential directory
PSPS="${abinit_pspdir}"
# Define the pseudopotential
PSP="${PSPS}/14si.psp"
# Set utilities directory
UTIL="${abinit_rundir}"
# Set machine name
MACH=$1

echo "testing the Abinit code on the $MACH platform"

#*******************************************************************************

#################################################
# set a date flag to make a new directory today #
#################################################

YYMMDD=`date "+DATE %Y%m%d" | cut -d" " -f2`

WORK_DIR="tmp-${MACH}_${YYMMDD}"
if test ! -d $WORK_DIR
then
  mkdir $WORK_DIR
fi

# **************************************** 
# Make preliminary operations in the working directory
cd $WORK_DIR

rm -f ab.files* 
tests="0 1 2 3 4 5 6 7 8 9"
for Tnum in $tests
do
cat >ab.files$Tnum <<End1
${abinit_inpdir}/paral/Input/t9.in
t9$Tnum.out
t9$Tnum.i
t9$Tnum.o
t9$Tnum
$PSP
End1
done
# ****************************************
# Jump to test to be run
#Case_1
#if test $2 -eq 1    This test does not work under C-shell ?!
#then
rm -f t90*
echo Case1
time $CODE_SEQ < ab.files0 2>&1 1> t90.log 
sleep 5
#exit
#fi

#Case_2
#if test $2 -le 2
#then
echo Case2
rm -f t91* t92*
time $CODE_SEQ < ab.files1 2>&1 1> t91.log &
time $CODE_SEQ < ab.files2 2>&1 1> t92.log 
sleep 5
#exit
#fi

#Case_3
#if test $2 -le 3
#then
echo Case3
rm -f t93* t94* t95*
time $CODE_SEQ < ab.files3 2>&1 1> t93.log &
time $CODE_SEQ < ab.files4 2>&1 1> t94.log &
time $CODE_SEQ < ab.files5 2>&1 1> t95.log
sleep 5
#exit
#fi

#Case_4
#if test $2 -le 4
#then
echo Case4
rm -f t96* t97* t98* t99*
time $CODE_SEQ < ab.files6 2>&1 1> t96.log &
time $CODE_SEQ < ab.files7 2>&1 1> t97.log &
time $CODE_SEQ < ab.files8 2>&1 1> t98.log &
time $CODE_SEQ < ab.files9 2>&1 1> t99.log 
sleep 5
#exit
#fi

#Case_5
#if test $2 -le 5
#then
#Analyze the output files
rm -f fldiff.set9.report
echo " " > fldiff.set9.report
for Tnum in $tests
do
diff -b t9$Tnum.out ${abinit_inpdir}/paral/Refs/t90.out > diff.t9$Tnum
echo " " >> fldiff.set9.report
echo "Case_9$Tnum" 
echo "Case_9$Tnum" >> fldiff.set9.report
$PERL $UTIL/fldiff.pl -ignore t9$Tnum.out ${abinit_inpdir}/paral/Refs/t90.out >> fldiff.set9.report
done
#fi
