---
authors: RShaltaf
---

# the Mrgscr utility  

This file explains the i/o parameters needed for the screening files (SCR) merging code (MRGSCR).  

The user is advised to be familiar with the main [[help:abinit]] especially
the GW part before reading the present file.

## 1 Introduction
  
The mrgscr is a utility that comes with the ABINIT code. It is used to merge
partial screening files. These files contain the screening calculated on some
selected q points generated using the input variables 'nqptdm' and 'qptdm'.
The mrgscr utility merge these files into a single file that contain the
screening on the full list of q points.

Like other utilities within ABINIT (e.g. mrgddb, mrgkk), the input is very
simple, and could be given directly at the screen, or more conveniently, piped from a file.

The user should give the number of SCRs that will be merged in the first line,
then the name of the output file in the second line, after what he/she shall
write the whole set of filenames for the SCRs to be merged, one on each line.

## 2 Why such utility can be useful?
  
The GW part of ABINIT code in its present version (5.x) is not yet
parallelized. The input variables 'nqptdm' and 'qptdm' are meant to be used as
a kind of manual parallelization by splitting the whole screening calculation
into several smaller jobs which can be submitted into several machines
simultaneously. mrgscr later can be used to merge the several output SCRs into a single file.

## 3 How does the mrgscr code work ?
  
First the code reads the header of the first partial file, then it calculates
the full list of q points that should exist in the full screening file. Then
it checks the consistency of the other files and the existing q points in each
file. Last it merges the files into a single file, but it worth to say here
that merging WILL NOT be successful and no output file will be generated if:

1) the files are not consistent, i.e, different PPs, different DM size etc.;
usually the files MUST be generated under same conditions for every single
input variable, except of course for 'nqptdm' and =20 'qptdm'.

2) the various input partial screening files can not form a complete file via
merging. This often happen if the user still have one or =20 more q points on
which screening still need to be calculated.

## 4 What if there is a repetition of one or more q points through different input partial screening files?
  
No problem, the code will include only one q point from every two repeated q
points, and will report it in the log file.

## 5 How to check the status of the resulting screening file, or partial files, what is there, what is needed, etc.?
  
mrgscr can do it, the user only need to use a two-line input file. The first
line of the input file should be 0, then the second line should contain the
name of the file to be checked. The result of the checking is reported through the log file.

Examples
    
    merging case
    3                    <== the number of files to be merged
    out                  <== name of out put file
    1_SCR                <== start:name of the files, to be merged
    2_SCR
    3_SCR
    1                    <== to merge q-points (2 = to merge frequencies) 
    
    
    checking case
    
    0                    <== just write zero
    1_SCR                <== the name of the file to be checked
