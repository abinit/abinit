## Source tree : Directories and files
## Introduction 

THIS FILE IS OUTDATED !!

This page provides a description of the organisation of the abinit package, in terms of subdirectories and their content.

The main directory of the abinit package is referred to as ~abinit hereafter, independently of its absolute location. From the Web site, a file *abinit-x.y.z.tar.gz can be dumped and, following the installation notes, a ~abinit directory will be created for this version.

See the files ~abinit/doc/tutorial/welcome.html and ~abinit/doc/users/new_user_guide.html for an introduction to the abinit package. Instructions to install the code, make the executable and run some tests can be found in the ~abinit/doc/install_notes directory. These files are also available on the Web :
http://www.abinit.org/documentation/helpfiles .

## There are numerous subdirectories in ~abinit, grouped by purpose:

  *  [A. doc](#A.doc/): documentation.
  *  [B. src](#B.The src/directories): core source of 
     ABINIT.
  *  [C. tests](#developers:srctree_The_tests_Tutorial_and_Psps_for_tests_directories): tests cases, tutorial input/output, and pseudopotentials.
  *  D. config/: the build system of ABINIT. 

We will now describe the content and use of the doc/, src/, and tests/ directories.

## A.doc/

This subdirectory contains all the files that describe information related to the ABINIT package. There are several <nowiki>HTML</nowiki> files, accessible directly from the Web, the other ones being plain-text (or in the markdown format). Many of these files are also ported to the Web frequently.

**A.1 doc/users**

Here is the place any user of ABINIT should look for documentation in. In particular, we would like to attract your attention to the following help files:

  *  new_user_guide.html: a guide for the new user.
  *  abinis_help.html: main code (sequential use).
  *  aim_help.html: AIM utility of the abinit package (Bader Atom-In-Molecule charge density analysis).
  *  anaddb_help.html: for the "Analysis of Derivative DataBases" utility.
  *  mrgddb_help.html: for the "Merge of Derivative DataBases" utility.
  *  respfn_help.html: a complementary help file for the response features of the core part of ABINIT.
  *  cut3d_help.html: for the "Cut3D" (cut in three dimensions) utility.

Other <nowiki>HTML</nowiki> files:

  *  acknowledgments.html: suggested acknowledgments and references to be inserted in scientific papers the results have been obtained thanks to ABINIT.
  *  bibliography.html: a list of papers that provide the theoretical framework the code is built upon.

Some of the other <nowiki>(non-HTML)</nowiki> files are worth to mention (in alphabetical order, and forgetting about possible version numbers or dates):

   * band2eps_help: the help file for band2eps (utility to produce encapsulated postscript graph of phonon band structures).
   * cut3d_help: the help file for cut3d (the analyser of 3D-files, like _DEN or _POT files)
   * ddbs_upgrade.txt: how to upgrade DDBs produced prior to ABINIT.
   * gwmanual.txt: a rough guide to the GW part of ABINIT.
   * newsp_help.txt: a brief help file for newsp, the wavefunction translator.
   * paral_use.txt: describes how to use the parallel version of abinit.
   * tuning.txt: describes how to reduce the memory needs of the code, and describes briefly the most time-consuming routines of the code.

# A.2 doc/tutorial

    * welcome.html: welcome to the new user.

# A.3 doc/input_variables

  *  many other <nowiki>HTML</nowiki> files describe the input variables (varbas.html, vardev.html, ..., varrlx.html), or constitute an index to these input variables (keyhr.html).

# A.4 doc/developers

Some of the other <nowiki>(non-HTML)</nowiki> files are worth to mention (in alphabetical order, and forgetting about possible version numbers or dates):

  *  context: the context of development of the ABINIT project.
  *  contributors.txt: the list of contributors to the ABINIT project.
  *  contributing.html: a description of the procedures followed for the development of the ABINIT package through collaboration of different groups of persons, based at different places around the world.
  *  FFT_in_parallel.txt: a work document for the FFT parallellisation.
  *  rules_coding.txt: the rules for Fortran90 coding, adopted by the ABINIT group.
  *  programmer_guide.txt: a guide for the programmer.
  *  rules_paral.txt: the rules followed for the parallelization within ABINIT.
  *  use_cpp.txt: how to use cpp in the ABINIT context.
  *  checklist.txt: a few things to remember, when you plan to give a contribution.

# A.5 doc/help_make

  *  make_help: a brief help file for the make in the ~abinit directory, also called automatically when make is issued without option.

# A.6 Other subdirectories

In addition to these files, several subdirectories are based in the doc/ directory:

  *  features: a directory that contains the cumulative lists of features, for different version of ABINIT.
  *  install_notes: for different versions of the code, the procedure to be followed to install them (compilation and execution of tests). These files are written in HTML, and available on the web site at ABINIT.
  *  macroave: documentation about the macroave utility of the ABINIT package.

  *  psp_infos: information on different pseudopotentials that can be used with the ABINIT package.\\ Presently, it contains:
     *   psp1.info: describes the format 1 of pseudopotentials.
     *   psp1.data: data about the Troullier-Martins pseudopotentials available on the site.
     *   psp3.info: describes the way to use a HGH pseudopotential file.
     *   psp45.info: describes the way to modify a psp file generated by the '92 code of M. Teter, or other even older files.
     *   psp5spinorbitinfo: the use of spin-orbit for format 5 psps.
     *   psp6.info: describes the way to use a fhi pseudopotential file.

  *  release_notes: a directory that contains the release notes for different versions of ABINIT.
  *  tutorial: the ABINIT tutorial, to be coupled with the content of the *~abinit/tests/tutorial directory.
  *  theory: information about the implementations inside ABINIT, see the README file for details.
  *  presentation: the source of the presentation of ABINIT, available on the web site, or as a poster.
  *  misc: a directory that contains miscellaneous information, not directly linked to ABINIT, but which might interest ABINIT users.


# B.The src/directories 

These directories contain the source of the routines in the ABINIT package, as well as their lists (abinit_src files). 

Subdirectories are labelled with a long name "xy_shortname", where xy is a two-digit number, e.g. 32_util .

The digit in the long name has the meaning of a level (in the link procedure). Routines in the level 01 (in the directory src/01_gsl_ext) do not use any other-level routine. The routines of level 68 use only routines of directories xy_name where xy is lower than 68. The level ordering continues until level 98, with *src/98_main*, which contains main programs.

There are many 30 src/* directories, so it is rather difficult to memorize in which directory one can find one particular file. Note that the name of all subroutines are different from each other in these src/* directories. This means that, in order to edit a particular routine, one ought not know in which directory it is located! Suppose that you want to modify the routine kpgsph.F90 using the vi editor. The easiest way is to simply issue vi */kpgsph.F90 and modify the file. Just type ls */kpgsph.F90 if you really need to know in which directory the file is located.

# B.1. src/10_defs

In particular, this directory contains:

  *  the defs_basis module, with the definitions of global parameters of the ABINIT package;

  *  constants that define data types (e.g. dp , is used in real(dp));
  *  numerical constants (e.g. zero , one, half, pi, ...);
  *  physical constants (e.g. Ha_eV=27.2113961_dp! 1 Hartree, in eV).

The other defs_* modules presently contain the definitions of derived datatypes.


# B.2. src/32_util

This directory contains small utility routines that are used in different places in the ABINIT package. These utility can be related to the treatment of characters, to basic mathematical operations, or to basic chemistry and physics data.

# B.3. src/98_main

This directory contains the Fortran programs top routines, e.g. abinit.F90, newsp.F90, mrgddb.F90, anaddb.F90, lwf.F90, conducti.F90, aim.F90, cut3d.F90, mrggkk.F90, macroave.F90, optic.F90 ...


# C. The tests/*, Tutorial, and Psps_for_tests directories 

This directory contains a set of tests of the code, aimed at testing whether the code is coherent in time (successive versions), and exercising many parts of the code. However, they do not examine its accuracy on physical problems, mainly because the number of plane waves used is too small, and some tests are not run to
convergence.

In order to understand the structure of this test directory, one might start by looking at
http://www.abinit.org/developers/dev-environment/buildbot/howto-add-a-new-test-in-the-test-suite 