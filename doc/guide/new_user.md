---
authors: DCA, XG, RC
---

# New user help file  

This page gives a beginner's introduction to the ABINIT resources, 
the package, and the main ABINIT applications.  

## Foreword
  
The ABINIT project is a group effort of dozens of people worldwide, whose
central outcome is the main ABINIT application, delivered with many other
files in the ABINIT package. The ABINIT project includes also resources
provided on the [ABINIT Web site](https://www.abinit.org) and 
the [github organization](https://github.com/abinit).

Before reading the present page, and get some grasp about the main ABINIT
application, you should get some theoretical background. In case you have
already used another electronic structure code, or a quantum chemistry code,
it might be sufficient to read the introduction of [[cite:Payne1992]].
If you have never used another electronic structure code or a Quantum
Chemistry package, you should complete such reading by going (at your own
pace) through the Chaps. 1 to 13 , and appendices L and M of R.M. Martin's book [[cite:Martin2004]].

After having gone through the present New User's Guide, you should follow the
[[tutorial:index|ABINIT tutorial]].

## Introduction
  
ABINIT is a package whose main program allows to find the total energy, charge
density and electronic structure of systems made of electrons and nuclei
(molecules and periodic solids) within Density Functional Theory, using
pseudopotentials and a planewave basis, or augmented plane waves, or even wavelets. 

Some possibilities of ABINIT go beyond Density Functional Theory,
i.e. the many-body perturbation theory (GW approximation the Bethe-Salpether
equation), Time-Dependent Density Functional Theory, Dynamical Mean-Field
Theory, the Allen-Heine-Cardona theory to find temperature-dependent electronic structure. 

ABINIT also includes options to optimize the geometry
according to the DFT forces and stresses, or to perform molecular dynamics
simulation using these forces, or to generate dynamical (vibrations - phonons)
properties, dielectric properties, mechanical properties, thermodynamical
properties, etc. In addition to the main ABINIT code, different utility
programs are provided.

We suppose that you have downloaded the ABINIT package from the Web site,
unpacked it and installed it. If not, you might nevertheless continue reading
the present Web page, just to get an overview, but it might prove more
fruitful to have first downloaded the ABINIT package and at least unpacked it,
see the [installation notes](../installation).

!!! note

    We will use the name "~abinit" to refer to the directory that contains the
    ABINIT package after download. In practice, a version number is appended to
    this name, to give for example: abinit-8.8.0. The ABINIT package versioning
    scheme is explained later in this file.

~abinit contains different subdirectories. For example, the present file, as
well as other descriptive files, should be found in ~abinit/doc/.
Other subdirectories will be described later.

## The main executable: abinit
  
After compilation, the main code will be present in the package as
~abinit/src/98_main/abinit (or perhaps at another place, depending on your installation).

To run abinit you need four things:

1. Access to the executable, abinit. 
2. An input file. 
3. A files file (list of file names in a file). 
4. A pseudopotential input file for each kind of element in the unit cell. 

With these items a job can be run.

The full list of input variables, all of which are provided in the single
input file, is given in the ABINIT [[varset:allvars|list of all variable]].
The detailed description of input variables is given in many "Variable Set" files, including:

  * Basic variables, [[varset:basic]]
  * Ground-state calculation variables, [[varset:gstate]]
  * GW variables, [[varset:gw]]
  * Files handling variables, [[varset:files]]
  * Parallelisation variables, [[varset:paral]]
  * Density Functional Perturbation Theory variables, [[varset:dfpt]]

A set of examples aimed at guiding the beginner is available in the [[tutorial:index|tutorials]].

Other test cases (more than 800 input files) can be found in the ~abinit/test
subdirectories, e.g. "fast", the "vX" series (v1, v2, ... v67mbpt, v7, v8),
"libxc", "paral", the "tutoX" series ...

Many different sorts of pseudopotentials can be used with ABINIT. 
Most of them can be found on the [atomic data files](https://www.abinit.org/downloads/atomic-data-files) 
page of the ABINIT web site. 
There are official recommended pseudopotentials tables 
(the PAW JTH table, and the norm-conserving table from ONCVPSP), and also some older sets of pseudopotentials. 
Information on pseudopotential files can be found in the [[help:abinit#5|ABINIT help file]],
the [[theory:pseudopotentials|Pseudopotential theory document]], on the [ABINIT wiki](https://wiki.abinit.org/doku.php?id=developers:pseudos),
and in the [[topic:PseudosPAW|PseudosPAW]] topics.

!!! warning

    A subset of existing pseudopotentials are used for test cases, and are located in the 
    ~abinit/tests/Psps_for_tests directory but they **are not recommended** for production. 

## Other programs in the package
  
In addition to abinit, there are utility programs provided in the package.
Some utility programs are written in Fortran (like the main abinit program), and
their sources is also in ~abinit/src/98_main. 
These include:

mrgddb and anaddb 
:   They allow one to post-process responses to atomic
    displacements and/or to homogeneous electric field, and/or to strain
    perturbation, as generated by abinit, to produce full phonon band structures,
    thermodynamical functions, piezoelectric properties, superconducting
    properties, to name a few. `mrgddb` is for "Merge of Derivative DataBases",
    while `anaddb` is for "Analysis of Derivative DataBases".

cut3d 
:   It can be used to post-process the three-dimensional density (or
    potential) files generated by abinit. It allows one to deduce charge density
    in selected planes (for isodensity plots), along selected lines, or at
    selected points. It allows one also to make the Hirshfeld decomposition of the
    charge density in "atomic" contributions.

fold2Bloch
:   It is used for unfolding first-principle electronic band structures

aim
:   It is also a post-processor of the three-dimensional density files
    generated by abinit. It performs the Bader Atom-In-Molecule decomposition of
    the charge density in "atomic" contributions.

conducti
:   It allows one to compute the frequency-dependent optical conductivity.

Some utility programs are not written in Fortran, but in Python. They are
contained in ~abinit/scripts, where post-processing (numerous tools) and pre-
processing scripts are distinguished. Some allows one to visualize ABINIT
outputs, like abinit_eignc_to_bandstructure.py.

## Other resources outside the ABINIT package
  
In addition to the ABINIT package, other resources can be obtained from the
[ABINIT github site](https://github.com/abinit). The sources of the latest
version of the ABINIT package are actually mirrored on this site, but for
other resources (not in the package) this is the only download point.

[AbiPy](https://github.com/abinit/abipy)
:   is an open-source library for analyzing the results produced by
    ABINIT (including visualisation), and for preparing input files and workflows
    to automate calculations (so-called high-throughput calculations).
    It provides interface with [pymatgen](http://pymatgen.org/), 
    developed by the [Materials Project](http://materialsproject.org/).
    Abinit tutorials based on AbiPy are available in the [abitutorials repository](https://github.com/abinit/abitutorials).

[PseudoDojo](http://www.pseudo-dojo.org/)
:   is a Python framework for generating and validating
    pseudopotentials (or PAW atomic data files). Normal ABINIT users benefit a lot
    from this project, since the ABINIT recommended table of norm-conserving
    pseudopotentials has been generated thanks to it. 
    The recommended PAW table is also provided via the pseudo-dojo interface.

[abiconfig](https://github.com/abinit/abiconfig)
:   is a holding area for configuration files used to
    configure/compile Abinit on clusters. You might benefit from it if you are 
    installing Abinit on a cluster.

[abiflows](https://github.com/abinit/abiflows)
:   provides flows for high-throughput calculations with ABINIT.


[abiconda](https://github.com/abinit/abiconda)
:   contains conda recipes to build Abinit-related packages (like AbiPy). 
    You might benefit from it if you install Abipy on your machine.


In addition to the resources that the ABINIT developer provide to the
community through the ABINIT packages, portal and Github, many ABINIT-independent 
commercial or free applications can be used to visualize ABINIT
outputs or interact with ABINIT. We provide a (not very well maintained) list
of links in [the last section of http://www.abinit.org/sponsors](http://www.abinit.org/sponsors).
Of course, you might get more by browsing the Web.

## Input variables to abinit
  
As an overview, the most important input variables, to be provided in the
input file, are listed below:

**Specification of the geometry of the problem, and types of atoms:**

[[natom]]
:       total number of atoms in unit cell

[[ntypat]]
:   number of types of atoms

[[typat]]([[natom]]):
:   sequence of integers, specifying the type of each atom.
    NOTE: the atomic coordinates ([[xangst]], [[xcart]] or [[xred]])
    must be specified in the same order

[[rprim]](3,3)
:   unscaled primitive translations of periodic cell;
    each COLUMN of this array is one primitive translation

[[xangst]](3,[[natom]])
:   cartesian coordinates (Angstrom) of atoms in unit cell
    NOTE: only used when [[xred]] and [[xcart]] are absent

[[xcart]](3,[[natom]])
:   cartesian coordinates (Bohr) of atoms in unit cell
    NOTE: only used when [[xred]] and [[xangst]] are absent

[[xred]](3,[[natom]])
:   fractional coordinates for atomic locations;
    NOTE: leave out if [[xangst]] or [[xcart]] is used

[[znucl]]([[ntypat]])
:   Nuclear charge of each type of element; must agree with
    nuclear charge found in psp file.

**Specification of the planewave basis set, Brillouin zone wavevector sampling, and occupation of the bands:**

[[ecut]]            
:       planewave kinetic energy cutoff in Hartree

[[kptopt]]
:       option for specifying the k-point grid
        if [[kptopt]]=1, automatic generation, using ngkpt and shiftk.

[[ngkpt]](3)
:       dimensions of the three-dimensional grid of k-points

[[occopt]]
:       set the occupation of electronic levels:
        =1 for semiconductors
        =3 ... 7  for metals

**Specification of the type of calculation to be done:**

[[ionmov]]
:       when [[ionmov]] = 0: the ions and cell shape are fixed
                        = 2: search for the equilibrium geometry
                        = 6: molecular dynamics

[[iscf]]
:       either a positive number for defining self-consistent
        algorithm (usual), or -2 for band structure in fixed potential

[[optdriver]]
:       when == 3 and 4: will do GW calculations (many-body perturbation theory)

[[rfelfd]]
:       when /= 0: will do response calculation to electric field

[[rfphon]]
:       when = 1: will do response calculation to atomic displacements

**Specification of the numerical convergency of the calculation:**

[[nstep]]
:    maximal number of self-consistent cycles (on the order of 20)

[[tolvrs]]
:    tolerance on self-consistent convergence

[[ntime]]
:       number of molecular dynamics or relaxation steps

[[tolmxf]]
:       force tolerance for structural relaxation in Hartree/Bohr

## Output files
  
Output from an abinit run shows up in several files and in the standard output. 
Usually one runs the command with a pipe of standard output to a log
file, which can be inspected for warnings or error messages if anything goes
wrong or otherwise can be discarded at the end of a run. The more easily
readable formatted output goes to the output file whose name is given in the
"files" file, i.e. you provide the name of the formatted output file. No error
message is reported in the latter file. On the other hand, this is the file
that is usually kept for archival purposes.

In addition, wavefunctions can be input (starting point) or output (result of
the calculation), and possibly, charge density and/or electrostatic potential,
if they have been asked for. These three sets of data are stored in unformatted binary files (native Fortran),
or in NetCDF format.  

The Density Of States (DOS) can also be an output as a formatted (readable) file.
An analysis of geometry can also be provided (GEO file).
The name of these files is constructed from a "root" name, that must be
different for input files and output files, and that is provided by the user,
to which the code will append a descriptor, like WFK for wavefunctions, DEN
for the density, POT for the potential, DOS for the density of states...

There are also different temporary files. A "root" name should be provided by
the user, from which the code generates a full name. Amongst these files,
there is a "status" file, summarizing the current status of advancement of the
code, in long jobs. The [[help:abinit|ABINIT help file]] contains more details.

## What does the code do?
  
The simplest sort of job computes an electronic structure for a fixed set of
atomic positions within a periodic unit cell. By electronic structure, we mean
a set of eigenvalues and wavefunctions which achieve the lowest DFT energy
possible for that basis set (that number of planewaves). 

The code takes the description of the unit cell and atomic positions and assembles a crystal
potential from the input atomic pseudopotentials, then uses either an input
wavefunction or simple gaussians to generate the initial charge density and
screening potential, then uses a self-consistent algorithm to iteratively
adjust the planewave coefficients until a sufficient convergence is reached in the energy.

Analytic derivatives of the energy with respect to atomic positions and unit
cell primitive translations yield atomic forces and the stress tensor. The
code can optionally adjust atomic positions to move the forces toward zero and
adjust unit cell parameters to move toward zero stress. It can performs
molecular dynamics. It can also be used to find responses to atomic
displacements and homogeneous electric field, so that the full phonon band
structure can be constructed.

## Versioning logic
  
We finish this "help for new user" with a brief explanation of the logic of ABINIT version releases.

The full name of a version has three digits (for example, 8.8.3). The first
digit is the slowly varying one (in average, it is changed after two or three
years). It indicates the major efforts and trends in that version. At the
level of 1.x.y ABINIT (before 2000 !), the major effort was placed on the
"ground-state" properties (total energy, forces, geometry optimisation,
molecular dynamics ...). With version 2.x.y , response-function features
(phonons, dielectric response, effective charges, interatomic force constants
...) were included. The main additional characteristics of version 3.x.y were
the distribution under the GNU General Public Licence, the set-up of the
documentation and help to the user through the Web site in html format, and
the availability of GW capabilities. The version 4.x.y put a lot of effort in
the speed of ABINIT (e.g. PAW), and its parallelisation. These historical
developments explain why the tests are gathered in directories "v1", "v2",
"v3", etc. Every 4 to 8 months, we release a "production version" of ABINIT in
which the second digit, an even number, is incremented, which usually goes
with additional features. A [release notes document](about/release-notes) is issued, with the list of
additional capabilities, and other information with respect to modifications
with the previous release. The odd second digits are used for internal
management only, so-called "development versions" of ABINIT (for example
8.9.0). Two versions differing by the last (third) digit have the same
capabilities, but the one with the largest last digit is more debugged than
the other: version 8.8.3 is more debugged than 8.8.2, but no new features has
been added (so likely, no additional bug!).

In order to start using ABINIT, please follow [[tutorial:index|the tutorial]]
