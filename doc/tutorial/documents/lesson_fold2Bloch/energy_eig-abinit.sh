#!/bin/sh

# This program will process case_EIG file of ABINIT and generate a row
# of eigenvalies and k-points
#
# Usage:
#   $ ./energy_eig.sh  case_EIG [-ry2ev] > output.dat
#   read file and analyse 3 lines
#
# (c) Oleg Rubel	Jan. 20, 2010

# debuging option: debug=0 - no debug; debug>1 - print extended output
debug=1

# check if the file name is speciefied
if [ ! -n "$1" ]
then
  echo 'please specify the file name: "./energy_eig.sh case.output1 #lines"'
  echo 'EXIT!'
  exit
else
  filein=$1
  if [ $debug -ge 1 ]; then echo filein=$filein; fi
fi

# create output file
fileout=$filein.dat
echo fileout=$fileout
if [ ! -e $fileout ]
then
  if [ $debug -ge 1 ]; then echo file $fileout will be created; fi
  touch $fileout
else
  if [ $debug -ge 1 ]; then echo file $fileout exists and will be removed; fi
  rm $fileout
  touch $fileout
fi

# determine total number of k-points
kmax=`grep 'nkpt=' $filein | cut -d "=" -f 2 | cut -d "k" -f 1`
if [ $debug -ge 1 ]; then echo kmax=$kmax; fi

# determine line number for each k-point
linelist=`grep -n ', kpt=' $filein | cut -d ":" -f 1`
if [ $debug -ge 1 ]; then echo linelist=$linelist; fi

# check how many lines are allocated for eigenvalues
line1=`echo $linelist | cut -d " " -f 1`
line2=`echo $linelist | cut -d " " -f 2`
let nlines=line2-line1-1
if [ $debug -ge 1 ]; then echo line1=$line1, line2=$line2, nlines=$nlines; fi

# main loop while $lastline < $linemax
for line in $linelist; do

  if [ $debug -ge 1 ]; then echo line=$line; fi

  # determine the k-point number
  knum=`sed -n ${line}p $filein | awk -F "kpt#" '{ print $2 }' | cut -d "," -f 1`
  if [ $debug -ge 1 ]; then echo knum=$knum; fi

  # determine weight for the k-point
  weight=`sed -n ${line}p $filein | awk -F "wtk=" '{ print $2 }' | cut -d "," -f 1`
  if [ $debug -ge 1 ]; then echo weight=$weight; fi

  # determine coordinates for the k-point
  kptc=`sed -n ${line}p $filein | awk -F "kpt=" '{ print $2 }' | cut -d "(" -f 1`
  if [ $debug -ge 1 ]; then echo kptc=$kptc; fi

  # read eigenvalues
  let linestart=line+1
  let lineend=line+nlines
  if [ $debug -ge 1 ]; then echo linestart=$linestart lineend=$lineend; fi
  eiglist=`sed -n ${linestart},${lineend}p $filein`
  if [ $debug -ge 1 ]; then echo eiglist=$eiglist; fi

  # append each eigenvalues to the output file
  echo $knum $weight $kptc $eiglist >> $fileout
#  for eig in $eiglist; do
#    if [ $ry2ev -eq 1 ]; then # Ry -> eV conversion
#      eig=`echo $eig | awk '{ printf "%14.8f", $1 * 13.6056923}'`
#    fi
#    echo $knum $weight $eig >> $fileout
#  done

done

# END
echo 'DONE, exit'
