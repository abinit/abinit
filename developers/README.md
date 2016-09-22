ABINIT developers and maintainers scripts
=========================================

Contains several scripts to ease development.


- bzr_helper
    Various script to facilitate the use of bzr.
    See the documentation within.

- mkroutine.sh
    Make a new F90 routine, with the correct robodoc header,
    and  structure that follows ABINIT rules.

- mkmodule.sh
    Make a new F90 module.


- FT77to90 and fixed_to_free
    `FT77to90` is a perl script that is able to translate a file written
    in Fortran77 fixed format to Fortran90 free format. What it does
    is relatively well explained. The csh script `fixed_to_free`
    is the driver of `FT77to90`, and slightly change its output.

- abirules.pl
    The source of the abirules script. The latter is able to enforce 
    automatically some of the ABINIT rules in F90 files.
 
- change
    A bash script, to change automatically some expression into another
    is a whole set of files, while making back-up copies of the old version.

- parents
    Locate all parents and children of a F90 routines, and write
    them in the routine.

- var-file-index.py
    Build the file `Infos/varfileindex.html` which refers all the
    input variables inside the input files of tests.