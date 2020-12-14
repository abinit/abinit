---
authors: MG
---

# abitk

Abitk is a Fortran executable located *~abinit/src/89_main/* that provides
a command line interface to perform very basic post-processing of Abinit output files. 
The tool supports both Fortran binary files as well as netcdf files
The syntax is:

    abitk COMMAND FILE [--option1 foo --option2 bar ...]

where COMMAND defines to operation to to performed, FILE is the output file from which data is extracted followed by command line
options depending on COMMAND.
Use `abitk --help` to get list of possible commands.
It support (mainly netcdf files).
The idea is to expose a simplified interface to some of the Fortran routines thus replacing

## Why a new EPH code?

