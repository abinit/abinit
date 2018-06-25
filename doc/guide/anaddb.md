---
authors: XG, DCA
---

# The anaddb utility  

This file explains the use and i/o parameters needed for the "Analysis of Derivative DataBase" code.

This code is able to compute interatomic force constants, but also, more
generally, many different physical properties from databases containing
derivatives of the total energy (Derivative DataBases - DDB).  
The user is not supposed to know how the Derivative DataBase (DBB) has been
generated. He/she should simply know what material is described by the DDB he/she wants to use.  

If he/she is interested in the generation of DDB, and wants to know more about
this topic, he/she will read different help files of the ABINIT package,
related to the [[help:abinit|main ABINIT executable]], to the
[[help:respfn|DFPT features of ABINIT]], and to the [[help:mrgddb|DDB merge tool]].

It will be easier to discover the present file with the help of the [tutorials](../tutorial),
especially the tutorials on [DFPT1](../tutorial/rf1) and [DFPT2](../tutorial/rf2).  

## 1 Introduction
  
In short, a Derivative DataBase contains a list of derivatives of the total
energy with respect to three kind of perturbations: phonons, electric field
and stresses. The present code analyses the DDB, and directly gives properties
of the material under investigation, like phonon spectrum, frequency-dependent
dielectric tensor, thermal properties.

Given an input file (parameters described below), the user must create a
"files" file which lists names for the files the job will require, including
the main input file, the main output file, the name of the DDB, and some other
file names optionally used for selected capabilities of the code.

The files file (called for example ab.files) could look like:
    
      anaddb.in  
      anaddb.out  
      ddb  
      band_eps  
      gkk  
      anaddb.ep  
      ddk  
     
In this example:  

  * the main input file is called "anaddb.in",   
  * the main output will be put into the file called "anaddb.out",   
  * the input DDB file is called "ddb",   
  * information to draw phonon band structures will go to band_eps  
  * the input GKK file is called "gkk" (used only for electron-phonon interactions)  
  * the base filename for electron-phonon output "anaddb.ep" (used only for electron-phonon interactions)  
  * the file name for ddk reference files: these are the GKK files generated in k-point derivative runs, 
    using the [[prtgkk]] abinit input variable (used only for electron-phonon transport calculations)

Other examples are given in the ~abinit/test/v2 directory. The latter three
filename information is often not used by anaddb. The maximal length of names
for the main input or output files is presently 264 characters.

The main executable file is called anaddb. Supposing that the "files" file is
called anaddb.files, and that the executable is placed in your working
directory, anaddb is run interactively (in Unix) with the command:

    anaddb < anaddb.files >& log
  
or, in the background, with the command

    anaddb < anaddb.files >& log &

where standard out and standard error are piped to the log file called "log"

!!! tip
    Piping the standard error, thanks to the '&' sign placed after '>' is
    **really important** for the analysis of eventual failures, when not due to
    ABINIT, but to other sources, like disk full problem.

The user can specify any names he/she wishes for any of these files. Variations of the
above commands could be needed, depending on the flavor of UNIX that is used
on the platform that is considered for running the code.

The syntax of the input file is strictly similar to the syntax of the main
abinit input files: the file is parsed, keywords are identified, comments are
also identified. However, the multidataset mode is not available.

## 2 Input variables
  
This ANADDB utility is able to perform many different tasks, each governed by
a selected set of input variables, with also some input variables common to
many of the different tasks. The 'flag' variables activates the different tasks 
e.g. [[dieflag@anaddb]], [[thmflag@anaddb]], [[elphflag@anaddb]]

The list of input variables for the anaddb input file are presented in the
[[varset:anaddb]] variable set. In order to discover them, it is easier to use
the different tutorials: start with the [second DFPT tutorial](../tutorial/rf2), then follow 
the [tutorial on elasticity and
piezoelectricity](../tutorial/elastic), the [tutorial on electron-phonon
interaction](../tutorial/eph), and the [tutorial on non-linear properties](../tutorial/nlo). 
