---
authors: XG, DCA, JWZ
---

# The abitk utility  

The user is advised to be familiar with the main [[help:abinit]] help file
before reading the present file.

## 1 Introduction

This simple utility, which is built along with ABINIT during the compilation phase,
simply reads output files and displays selected information. It is designed primarily
to read output files saved in `netcdf` format. The options can be seen by running

    abitk -h
  
## 2 Typical usage

In a typical ABINIT run, an input file `run.abi` produces a number of output files,
including `run.abo` and `run.abo_GSR.nc`. The latter of these contains information
about the Ground State Run in `netcdf` format. One can examine this file with
  
    abitk crystal_print run.abo_GSR.nc

for example, which will extract and print information about the crystal structure used
in the run.
