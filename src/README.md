Adding a new directory
======================

Instructions can be adapted to move a directory inside src (see point 4)
Do this only if you have the autotools.

1. create your directory, fill it, and add the files with git

2. in `config/specs/corelibs.conf`, declare the new directory,
   its attributes and its dependencies upon the external libraries

3. in `config/specs/binaries.conf`, specify the main codes that
   depends on the routines contained in this new directory

Then, in order for the make to work, issue:

    ~abinit/config/scripts/makemake

then `configure` and `make`

In case of a move, check that the doc (or other files) does not need any update.
Advice: find the occurrences thanks to:

    grep name_of_directory *
    grep name_of_directory */*
