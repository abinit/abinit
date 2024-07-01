---
authors: DRH
---

# Tutorial on elastic properties

## Elastic and piezoelectric properties.

This tutorial shows how to calculate physical properties related to strain, for
an insulator and a metal:

  * the rigid-atom elastic tensor
  * the rigid-atom piezoelectric tensor (insulators only)
  * the internal strain tensor
  * the atomic relaxation corrections to the elastic and piezoelectric tensor

You should complete tutorials [RF1](/tutorial/rf1) and [RF2](/tutorial/rf2)
to introduce the density-functional perturbation theory (DFPT) features of
ABINIT before starting this tutorial. You will learn to use additional DFPT
features of ABINIT, and to use relevant parts of the associated codes Mrgddb and Anaddb.

Visualisation tools are NOT covered in this tutorial.
Powerful visualisation procedures have been developed in the Abipy context,
relying on matplotlib. See the README of [Abipy](https://github.com/abinit/abipy)
and the [Abipy tutorials](https://github.com/abinit/abitutorials).

This tutorial should take about two hours.

[TUTORIAL_README]

## 1 The ground-state geometry of (hypothetical) wurtzite AlP

*Before beginning, you might consider working in a different subdirectory as for the other tutorials.
Why not create Work_elast in \$ABI_TESTS/tutorespfn/Input?*
You should also copy the file *telast_1.abi* from *\$ABI_TESTS/tutorespfn/Input* to *Work_elast*.

```sh
cd $ABI_TESTS/tutorespfn/Input
mkdir Work_elast
cd Work_elast
cp ../telast_1.abi .
```

You may wish to start the calculation (less than
one minute on a standard 3GHz machine) before you read the following.

    abinit telast_1.abi > telast_1.log &

Then, you should open your input file *telast_1.abi* with an editor and examine it as you read this discussion.

{% dialog tests/tutorespfn/Input/telast_1.abi %}

The hypothetical wurtzite structure for AlP retains the tetrahedral
coordination of the atoms of the actual zincblende structure of AlP, but has a hexagonal lattice.
It was chosen for this tutorial because the atomic
positions are not completely determined by symmetry.
Both the atomic positions and the lattice constants should be optimized before beginning DFPT
calculations, especially those related to strain properties.
While GS structural optimization was treated in tutorials 1-3, we are introducing a few
new features here, and you should look at the following new input variables which will be discussed below:

  * [[getxred]]
  * [[iatfix]]
  * [[natfix]]
  * [[strfact]]

There are two datasets specified in *telast_1.abi*. First, let us examine the
common input data. We specify a starting guess for [[acell]], and give an
accurate decimal specification for [[rprim]]. The definition of the atom
types and atoms follows [tutorial DFPT1](/tutorial/rf1). The reduced atomic
positions [[xred]] are a starting approximation, and will be replaced by our
converged results in the remaining input files, as will [[acell]].

We will work with a fixed plane wave cutoff [[ecut]] (=6 Ha), but introduce
[[ecutsm]] (0.5 Ha) as in [tutorial 3](/tutorial/base3) to smear the cutoff,
which produces smoothly varying stresses as the lattice parameters are
optimized. We will keep the same value of [[ecutsm]] for the DFPT calculations
as well, since changing it from the optimization run value could reintroduce
non-zero forces and stresses. For the k-point grid, we must explicitly specify
[[shiftk]] since the default value results in a grid shifted so as to break
hexagonal symmetry. The RF strain calculations check this, and will exit with
an error message if the grid does not have the proper symmetry.
The self-consistency procedures follow [tutorial RF1](/tutorial/rf1).

Dataset 1 optimizes the atomic positions keeping the lattice parameters fixed,
setting [[ionmov]]=2 as in [tutorial 1](/tutorial/base1). The optimization
steps proceed until the maximum force component on any atom is less than
[[tolmxf]]. It is always advised to relax the forces before beginning the
lattice parameter optimization. Dataset 2 optimizes the lattice parameters
with [[optcell]]=2 as in [tutorial 3](/tutorial/base3). However, tutorial 3
treats cubic Si, and the atom positions in reduced coordinates remained
fixed. In the present, more general case, the reduced atomic coordinates must
be reoptimized as the lattice parameters are optimized. Note that it is
necessary to include [[getxred]] = -1 so that the second dataset is
initialized with the relaxed coordinates. Coordinate and lattice parameter
optimizations actually take place simultaneously, with the computed stresses
at each step acting as forces on the lattice parameters. We have introduced
[[strfact]] which scales the stresses so that they may be compared with the
same [[tolmxf]] convergence test that is applied to the forces. The default
value of 100 is probably a good choice for many systems, but you should be
aware of what is happening.

From the hexagonal symmetry, we know that the positions of the atoms in the
a-b basal plane are fixed. However, a uniform translation along the c axis of
all the atoms leaves the structure invariant. Only the relative displacement
of the Al and As planes along the c axis is physically relevant. We will fix
the Al positions to be at reduced c-axis coordinates 0 and 1/2 (these are
related by symmetry) by introducing [[natfix]] and [[iatfix]] to constrain the
structural optimization. This is really just for cosmetic purposes, since
letting them all slide an arbitrary amount (as they otherwise would) won't
change any results. However, you probably wouldn't want to publish the results
that way, so we may as well develop good habits.

Now we shall examine the results of the structural optimization run. As
always, we should first examine the log file to make sure the run has
terminated cleanly. There are a number of warnings, but none of them are
apparently serious. Next, let us edit the output file, *telast_1.abo*.

{% dialog tests/tutorespfn/Refs/telast_1.abo %}

The first thing to look for is to see whether Abinit recognized the symmetry of the
system. In setting up a new data file, it's easy to make mistakes, so this is a valuable check.
We see

    DATASET    1 : space group P6_3 m c (#186); Bravais hP (primitive hexag.)

which is correct. Next, we confirm that the structural optimization converged.
The following lines tell us that things are OK.
From dataset 1 :

    At Broyd/MD step   4, gradients are converged :
    max grad (force/stress) = 2.7264E-08 < tolmxf= 1.0000E-06 ha/bohr (free atoms)

and from dataset 2 :

    At Broyd/MD step  14, gradients are converged :
    max grad (force/stress) = 5.3969E-07 < tolmxf= 1.0000E-06 ha/bohr (free atoms)

We can also confirm that the stresses are relaxed:

    Cartesian components of stress tensor (hartree/bohr^3)
     sigma(1 1)=  1.12504637E-09  sigma(3 2)=  0.00000000E+00
     sigma(2 2)=  1.12504637E-09  sigma(3 1)=  0.00000000E+00
     sigma(3 3)= -2.49308610E-09  sigma(2 1)=  0.00000000E+00

Now would be a good time to copy *telast_2.abi* into your
working directory, since we will use the present output to start the next run.
Locate the optimized lattice parameters and reduced atomic coordinates near
the end of *telast_1.abo*:

     acell2    7.2488954246E+00  7.2488954246E+00  1.1879499870E+01 Bohr

     xred2     3.3333333333E-01  6.6666666667E-01  0.0000000000E+00
               6.6666666667E-01  3.3333333333E-01  5.0000000000E-01
               3.3333333333E-01  6.6666666667E-01  3.7517446813E-01
               6.6666666667E-01  3.3333333333E-01  8.7517446813E-01

With your editor, copy and paste these into *telast_2.abi* at the indicated
places in the "Common input data" area. Be sure to change acell2 and xred2 to
acell and xred since these common values will apply to all datasets in the next set of calculations.

## 2 Response-function calculations of several second derivatives of the total energy

We will now compute second derivatives of the total energy (2DTE's) with
respect to all the perturbations we need to compute elastic and piezoelectric
properties. You may want to review the first paragraphs of the [[help:respfn]]
which you studied in [tutorial RF1](/tutorial/rf1).
We will introduce only one new input variable for the strain perturbation,

  * [[rfstrs]]

The treatment of strain as a perturbation has some subtle aspects. It would be
a good idea to read  Metric tensor formulation of strain in density-functional
perturbation theory, by D. R. Hamann, Xifan Wu, Karin M. Rabe, and David Vanderbilt [[cite:Hamann2005]]
especially Sec. II and Sec. IV. We will do all the RF calculations you learned in [tutorial RF1](/tutorial/rf1) together
with strain, so you should review the variables

  * [[rfphon]]
  * [[rfatpol]]
  * [[rfdir]]
  * [[rfelfd]]

If not yet done (see end of previous section), copy *telast_2.abi* into *Work_elast* and start the
calculation while you read (less than 2 minutes on a standard 3GHz machine). We now
assume that you know which command to use to launch ABINIT.
Look at *telast_2.abi* in your editor to follow the discussion.

{% dialog tests/tutorespfn/Input/telast_2.abi %}

This has been set up as a self-contained calculation with three datasets. The
first is simply a GS run to obtain the GS wave functions we will need for the
DFPT calculations. We have removed the convergence test from the common input
data to remind ourselves that different tests are needed for different
datasets. We set a tight limit on the convergence of the self-consistent
potential with [[tolvrs]]. Since we have specified [[nband]]=8, all the bands
are occupied and the potential test also assures us that all the wave
functions are well converged. This issue will come up again in the section on
metals. We could have used the output wave functions *telast_1o_DS2_WFK* as
input for our RF calculations and skipped dataset 1, but redoing the GS
calculation takes relatively little time for this simple system.

Dataset 2 involves the calculation of the derivatives of the wave functions
with respect to the Brillouin-zone wave vector, the so-called ddk wave
functions. Recall that these are auxiliary quantities needed to compute the
response to the [[lesson:rf1#5|electric field perturbation]] and
introduced in [tutorial RF1](/tutorial/rf1). It would be a good idea to review the relevant
parts of [[help:respfn#1|section 1]] of the respfn_help file.

Examining this section of *telast_2.abi*, note that electric
field as well as strain are uniform perturbations, so that they are defined only for
[[qpt]] = 0 0 0. [[rfelfd]] = 2 specifies that we want the ddk calculation to be
performed, which requires [[iscf]] = -3. The ddk wave functions will be used to
calculate both the piezoelectric tensor and the Born effective charges, and in
general we need them for **k** derivatives in all three (reduced) directions,
[[rfdir]] = 1 1 1 (that is the default). Since there is no potential self-consistency in the ddk
calculations, we must specify convergence in terms of the wave function
residuals using [[tolwfr]].

Finally, dataset 3 performs the actual calculations of the needed 2DTE's for
the elastic and piezoelectric tensors. Setting [[rfphon]] = 1 turns on the
atomic displacement perturbation, which we need for all atoms (see [[rfatpol]])
and all directions (see [[rfdir]]). Abinit will calculate first-order
wave functions for each atom and direction in turn, and use those to calculate
2DTE's with respect to all pairs of atomic displacements and with respect to
one atomic displacement and one component of electric field. These quantities,
the interatomic force constants (at $\Gamma$) and the Born effective charges will
be used later to compute the atomic relaxation contribution to the elastic and
piezoelectric tensor.

First-order wave functions for the strain perturbation are computed next.
Setting [[rfstrs]] = 3 specifies that we want both uniaxial and shear strains
to be treated, and [[rfdir]] = 1 1 1 cycles through strains xx, yy, and zz for
uniaxial and yz, xz, and xy for shear. We note that while other perturbations
in Abinit are treated in reduced coordinates, strain is better dealt with in
Cartesian coordinates for reasons discussed in the reference cited above.
These wave functions are used to compute three types of 2DTE's. Derivatives
with respect to two strain components give us the so-called rigid-ion elastic
tensor. Derivatives with respect to one strain and one electric field
component give us the rigid-ion piezoelectric tensor. Finally, derivatives
with respect to one strain and one atomic displacement yield the internal-strain 
force-response tensor, an intermediate quantity that will be necessary
to compute the atomic relaxation corrections to the rigid-ion quantities. As
in [tutorial DFPT1](/tutorial/rf1), we specify convergence in terms of the residual of the
potential (here the first-order potential) using [[tolvrs]].

Your run should have completed by now. Abinit should have created quite a few files,
among which:

  * telast_2.abo (main output file)
  * telast_2o_DS1_DDB (first derivatives of the energy from GS calculation)
  * telast_2o_DS3_DDB (second derivatives from the RF calculation)
  * telast_2o_DS1_WFK (GS wave functions)
  * telast_2o_DS2_1WF* (ddk wave functions)
  * telast_2o_DS3_1WF* (RF first-order wave functions from various perturbations)

The derivative database DDB files are ascii and
readable, but primarily for subsequent analysis by anaddb which we will
undertake in the next section. Finally, the various wave function binary files
are primarily of use for subsequent calculations, where they could cut the
number of needed iterations in, for example, convergence testing. 
File names have been generated according to the following convention.
After the "root" name (which is by default taken from the name of the .abi file),
follows the number of the dataset producing the file.
Finally, the first-order wave function 1WF files have a final "pertcase"
number described in [[help:respfn#1|section 1]] of the [[help:respfn|respfn_help file]].
While *telast_2.abi* specifies all atomic displacements, only the 
symmetry-inequivalent perturbations are treated, so the "pertcase" list is incomplete.
All cases specified in the input data are treated for the strain perturbation.

First, take a look at the end of the \*log file to make sure the run
has completed without error. You might wish to take a look at the WARNING's,
but they all appear to be harmless. Next, edit your *telast_2.abo* file.
Searching backwards for ETOT you will find

         iter   2DEtotal(Ha)        deltaE(Ha) residm    vres2
    -ETOT  1   2.6864034257676     -1.394E+01 7.628E+01 3.415E+02
     ETOT  2   1.5949143903393     -1.091E+00 3.669E-02 4.465E+00
     ETOT  3   1.5767914745018     -1.812E-02 5.026E-05 3.449E-01
     ETOT  4   1.5758067303228     -9.847E-04 9.600E-07 6.237E-03
     ETOT  5   1.5757929703648     -1.376E-05 6.251E-08 3.303E-05
     ETOT  6   1.5757928924185     -7.795E-08 1.128E-09 3.010E-07
     ETOT  7   1.5757928914135     -1.005E-09 9.367E-11 6.841E-09
     ETOT  8   1.5757928913731     -4.046E-11 1.281E-12 1.312E-10
     ETOT  9   1.5757928913723     -7.105E-13 1.065E-13 8.070E-12

     At SCF step    9       vres2   =  8.07E-12 < tolvrs=  1.00E-10 =>converged.

Abinit is solving a set of Schroedinger-like equations for the first-order
wave functions, and these functions minimize a variational expression for the
2DTE (technically, they are called self-consistent Sternheimer equations).
The  energy  convergence looks similar to that of GS calculations.  The fact
that *vres2*, the residual of the self-consistent first-order potential, has
reached [[tolvrs]] well within [[nstep]] (40) iterations indicates that the
2DTE calculation for this perturbation (xy strain) has converged. It would
pay to examine a few more cases for different perturbations (unless you have
looked through all the warnings in the log).

Another convergence item to examine in your .abo file is

     Seventeen components of 2nd-order total energy (hartree) are
     1,2,3: 0th-order hamiltonian combined with 1st-order wavefunctions
         kin0=   1.67341861E+01 eigvalue=  -4.11015159E-01  local=  -3.22857213E+00
     4,5,6,7: 1st-order hamiltonian combined with 1st and 0th-order wfs
     loc psp =  -9.55499839E+00  Hartree=   4.67777679E+00     xc=  -7.47227780E-01
         kin1=  -2.31525775E+01
     8,9,10: eventually, occupation + non-local contributions
        edocc=   0.00000000E+00     enl0=   6.39144174E-01   enl1=  -6.46481544E-03
     1-10 gives the relaxation energy (to be shifted if some occ is /=2.0)
       erelax=  -1.50497487E+01
     11,12,13 Non-relaxation  contributions : frozen-wavefunctions and Ewald
      fr.hart=  -1.37404595E-01   fr.kin=   1.27503317E+01 fr.loc=   4.84873306E-01
     14,15,16 Non-relaxation  contributions : frozen-wavefunctions and Ewald
      fr.nonl=   4.27603780E-01    fr.xc=  -3.16560587E-02  Ewald=   3.13179353E+00
     17 Non-relaxation  contributions : pseudopotential core energy
      pspcore=   0.00000000E+00
     Resulting in :
     2DEtotal=    0.1575792891E+01 Ha. Also 2DEtotal=    0.428795052510E+02 eV
        (2DErelax=   -1.5049748718E+01 Ha. 2DEnonrelax=    1.6625541610E+01 Ha)
        (  non-var. 2DEtotal :    1.5757929293E+00 Ha)

This detailed breakdown of the contributions to 2DTE is probably of limited
interest, but you should compare "2DEtotal" and "non-var. 2DEtotal" from the
last three lines. While the first-order wave function for the present
perturbation minimizes a variational expression for the second derivative
with respect to this perturbation as we just saw, the various 2DTE given as
elastic tensors, etc. in the output and in the DDB file are all computed using
non-variational expressions. Using the non-variational expressions, mixed
second derivatives with respect to the present perturbation and all other
perturbations of interest can be computed directly from the present 
first-order wave functions. The disadvantage is that the non-variational result
has errors which are linearly proportional to convergence errors in the GS and
first-order wave functions. Since errors in the variational 2DEtotal are
second-order in wave-function convergence errors, comparing this to the non-variational
result for the diagonal second derivative will give an idea of the
accuracy of the latter and perhaps indicate the need for tighter convergence
tolerances for both the GS and RF wave functions.
This is discussed in X. Gonze and C. Lee, Phys. Rev. B 55, 10355 (1997) [[cite:Gonze1997a]], Sec. II.
For an atomic-displacement perturbation, the corresponding breakdown of the 2DTE is headed
"Thirteen components."

Now let us take a look at the results we want, the various 2DTE's. They begin by

     ==> Compute Derivative Database <==

      2nd-order matrix (non-cartesian coordinates, masses not included,
       asr not included )
      cartesian coordinates for strain terms (1/ucvol factor
       for elastic tensor components not included)
         j1       j2             matrix element
      dir pert dir pert     real part     imaginary part
 
       1    1   1    1         5.7740048299         0.0000000000
       1    1   2    1        -2.8870024150         0.0000000000
       1    1   3    1         0.0000000000         0.0000000000
                .....

These are the "raw" 2DTE's, in reduced coordinates for atom-displacement and
electric-field perturbations, but cartesian coordinates for strain
perturbations. The same results with the same organization appear in the file
*telast_2_DS3_DDB* which will be used later as input for automated analysis and
converted to more useful notation and units by anaddb. A breakout of various
types of 2DTE's follows (all converted to Cartesian coordinates and in atomic units):

      Dynamical matrix, in cartesian coordinates,
       if specified in the inputs, asr has been imposed
         j1       j2             matrix element
      dir pert dir pert     real part    imaginary part

       1    1   1    1         0.1098826675         0.0000000000
       1    1   2    1         0.0000000000         0.0000000000
       1    1   3    1         0.0000000000         0.0000000000
                .....

This contains the interatomic force constant data that will be used later to
include atomic relaxation effects. "asr" refers to the acoustic sum rule,
which basically is a way of making sure that forces sum to zero when an atom is displaced.

      Effective charges, in cartesian coordinates,
      (from phonon response)
       if specified in the inputs, charge neutrality has been imposed
         j1       j2             matrix element
      dir pert dir pert     real part    imaginary part

       1    6   1    1         2.0670263917         0.0000000000
       2    6   1    1         0.0000000000         0.0000000000
       3    6   1    1         0.0000000000         0.0000000000
                .....

The Born effective charges will be used to find the atomic relaxation
contributions of the piezoelectric tensor.

      Rigid-atom elastic tensor , in cartesian coordinates,
         j1       j2             matrix element
      dir pert dir pert     real part    imaginary part

       1    7   1    7         0.0078004875         0.0000000000
       1    7   2    7         0.0019706468         0.0000000000
       1    7   3    7         0.0011850679         0.0000000000
                .....

The rigid-atom elastic tensor is the 2DTE with respect to a pair of strains.
We recall that "pert" = natom+3 and natom+4 for unaxial and shear strains, respectively.

      Internal strain coupling parameters, in cartesian coordinates,
       zero average net force deriv. has been imposed
         j1       j2             matrix element
      dir pert dir pert     real part    imaginary part

       1    1   1    7         0.1464077350         0.0000000000
       1    1   2    7        -0.1464077350         0.0000000000
       1    1   3    7         0.0000000000         0.0000000000
                .....

These 2DTE's with respect to one strain and one atomic displacement are needed
for atomic relaxation corrections to both the elastic tensor and piezoelectric
tensor. While this set of parameters is of limited direct interest, it should
be examined in cases when you think that high symmetry may eliminate the need
for these corrections. You are probably wrong, and any non-zero term indicates a correction.

      Rigid-atom proper piezoelectric tensor, in cartesian coordinates,
      (from strain response)
         j1       j2             matrix element
      dir pert dir pert     real part    imaginary part

       1    6   1    7         0.0000000000         0.0000000000
       1    6   2    7         0.0000000000         0.0000000000
       1    6   3    7        -0.0000000000         0.0000000000
       1    6   1    8         0.0000000000         0.0000000000
       1    6   2    8         0.0070195611         0.0000000000
       1    6   3    8         0.0000000000         0.0000000000

Finally, we have the piezoelectric tensor, the 2DTE with respect to one strain
and one uniform electric field component. (Yes, there are non-zero elements.)

## 3 ANADDB calculation of atom-relaxation effects

In this section, we will run the program anaddb, which analyzes DDB files
generated in prior RF calculations. You should copy *telast_3.abi* 
in your *Work_elast* directory. You should now go to the [[help:anaddb|anaddb help file]]
introduction. The bulk of the material in this help file is contained in the
description of the variables. You should read the descriptions of

  * [[anaddb:elaflag]],
  * [[anaddb:piezoflag]],
  * [[anaddb:instrflag]],
  * [[anaddb:chneut]].

For the theory underlying the incorporation of atom-relaxation corrections, it
is recommended you see  X. Wu, D. Vanderbilt, and D. R. Hamann [[cite:Wu2005]].

Anaddb can do lots of other things, such as calculate the frequency-dependent
dielectric tensor, interpolate the phonon spectrum to make nice phonon
dispersion plots, calculate Raman spectra, etc., but we are focusing on the
minimum needed for the elastic and piezoelectric constants at zero electric field.

We also mention that [[help:mrgddb|mrgddb]]
is another utility program that can be used to combine DDB files generated in
several different datasets or in different runs into a single DDB file that
can be analyzed by anaddb. One particular usage would be to combine the DDB
file produced by the GS run, which contains first-derivative information such
as stresses and forces with the RF DDB. It is anticipated that anaddb in a
future release will implement the finite-stress corrections to the elastic
tensor discussed in [notes by A. R. Oganov](/theory/documents/elasticity-oganov.pdf) .

Now would be a good time to edit *telast_3.abi* and observe that it is very
simple, consisting of nothing more than the four variables listed above set to
appropriate values.

{% dialog tests/tutorespfn/Input/telast_3.abi %}

A *telast_3.files* file is also needed for anaddb. Copy it from the tests/tutorespfn/Input directory.

{% dialog tests/tutorespfn/Input/telast_3.files %}

The first two lines specify the .abi and .abo files, the third line specifies the DDB file, and the
last lines are dummy names which would be used in connection with other
capabilities of anaddb. Now you should run the calculation:

    anaddb <telast_3.files >& telast_3.log

This calculation should only take a few seconds. You should edit the log file,
go to the end, and make sure the calculation terminated without error. Next,
examine *telast_3.abo*. After some header information, we come to tables giving
the "force-response" and "displacement-response" internal strain tensors.
These represent, respectively, the force on each atom and the displacement of
each atom in response to a unit strain of the specified type. These numbers
are of limited interest to us, but represent important intermediate quantities
in the treatment of atomic relaxation (see the X. Wu [[cite:Wu2005]] paper cited above).

Next, we come to the elastic tensor output:

     Elastic Tensor (clamped ion) (unit:10^2GP):

       2.2949823   0.5797842   0.3486590  -0.0000000  -0.0000000   0.0000001
       0.5797842   2.2949822   0.3486590  -0.0000000   0.0000000   0.0000001
       0.3486589   0.3486589   2.4696020  -0.0000000   0.0000000   0.0000001
       0.0000000   0.0000000   0.0000000   0.5821756   0.0000000   0.0000000
      -0.0000000   0.0000000   0.0000000   0.0000000   0.5821756   0.0000000
       0.0000000   0.0000000   0.0000000   0.0000000   0.0000000   0.8575990

     Elastic Tensor (relaxed ion) (unit:10^2GP):
      (at fixed electric field boundary condition)

       1.8986339   0.7520966   0.5652395  -0.0000000  -0.0000000   0.0000001
       0.7520966   1.8986339   0.5652395  -0.0000000   0.0000000   0.0000001
       0.5652395   0.5652395   2.0508555  -0.0000000   0.0000000   0.0000001
       0.0000000   0.0000000   0.0000000   0.4487908   0.0000000   0.0000000
       0.0000000  -0.0000000  -0.0000000   0.0000000   0.4487909  -0.0000000
       0.0000000  -0.0000000   0.0000000   0.0000000  -0.0000000   0.5732684


While not labeled, the rows and columns 1-6 here represent xx, yy, zz, yz, xz,
xy strains and stresses in the conventional Voigt notation.
The clamped-ion results were calculated in the telast_2 RF run, and are simply
converted to standard GPa units by anaddb (the terms "clamped ion," "clamped
atom," and "rigid atom" used in various places are interchangeable, similarly for "relaxed.")
The relaxed-ion result was calculated by anaddb by combining 2DTE's for
internal strain and interatomic force constants which are stored in the input
DDB file. Comparing the clamped and relaxed results, we see that all the
diagonal elastic constants have decreased in value.
This is plausible, since allowing the internal degrees of freedom to relax
should make a material less stiff. These tensors should be symmetric, and
certain tensor elements should be zero or identical by symmetry.
It's a good idea to check these properties against a standard text such as  J.
F. Nye, Physical Properties of Crystals (Oxford U. P., Oxford 1985) [[cite:Nye1985]].
Departures from expected symmetries (there are a few in the last decimal place
here) are due to either convergence errors or, if large, incorrectly specified
geometry (however, see the final comments on symmetry  below).

Later, we will obtain the (3,3) component of the Elastic Tensor (clamped ion)
from finite differences. Its value is 246.96020 GPa (note the unit:10^2GP indication above).

Next in *telast_3.abo we find the piezoelectric tensor results:

     Proper piezoelectric constants (clamped ion) (unit:c/m^2)

          0.00000000      0.00000000      0.36031593
          0.00000000      0.00000000      0.36031593
         -0.00000000     -0.00000000     -0.69614902
          0.00000000      0.40162256      0.00000000
          0.40162251      0.00000000      0.00000000
          0.00000000      0.00000000     -0.00000001

     ddb_piezo : WARNING -
      Acoustic sum rule violation met : the eigenvalues of accoustic mode
      are too large at Gamma point
      Increase cutoff energy or k-points sampling.
      The three eigenvalues are:   -1.016026E-04   -1.389885E-05   -2.054240E-05

     Proper piezoelectric constants (relaxed ion) (unit:c/m^2)

          0.00000000     -0.00000000     -0.12745531
         -0.00000000      0.00000000     -0.12745533
          0.00000000      0.00000000      0.24692958
          0.00000000     -0.15011587      0.00000000
         -0.15011568      0.00000000      0.00000000
         -0.00000000      0.00000000     -0.00000002

The 3 columns here represent x, y, and z electric polarization, and the 6 rows
the Voigt strains. The clamped-ion result was calculated in the telast_2 RF
run, and is simply scaled to conventional units by anaddb. The ion relaxation
contributions are based on 2DTE's for internal strain, interatomic force
constants, and Born effective charges, and typically constitute much larger
corrections to the piezoelectric tensor than to the elastic tensor. Once
again, symmetries should be checked (The slight discrepancies seen here can
be removed by setting tolvrs3 = 1.0d-18 in *telast_2.abi*). One should be aware
that the piezoelectric tensor is identically zero in any material which has a center of symmetry.

There is also a WARNING message. For the purpose of a tutorial, the [[ecut]]
has not been selected large enough. Changing [[ecut]] from 6 Ha to 15 Ha will eliminate this problem.

Since we are dealing with a hypothetical material, there is no experimental
data with which to compare our results. In the next section, we will calculate
a few of these numbers by a finite-difference method to gain confidence in the RF approach.

## 4 Finite-difference calculation of elastic and piezoelectric constants

You should copy *telast_4.abi* into your *Work_elast* directory.

{% dialog tests/tutorespfn/Input/telast_4.abi %}

Editing *telast_4.abi*, you will see that it has four datasets, the first two
with the c-axis contracted 0.01% and the second two with it expanded 0.01%,
which we specified by changing the third row of [[rprim]]. The common data is
essentially the same as telast_2.abi, and the relaxed [[acell]] values and
[[xred]] from telast_1.abo have already been included. Datasets 1 and 3 do the
self-consistent convergence of the GS wave functions for the strained lattices
and compute the stress. Datasets 2 and 4 introduce a new variable.

  * [[berryopt]]

Electric polarization in solids is a subtle topic which has been
rigorously resolved thirty years ago. It is now understood to be a bulk property, and to be
quantitatively described by a Berry phase formulation introduced by R. D.
King-Smith and D. Vanderbilt, Phys. Ref. B 47, 1651(1993) [[cite:Kingsmith1993]]. It can be
calculated in a GS calculation by integrating the gradient with respect to
**k** of the GS wave functions over the Brillouin zone. In GS calculations,
the gradients are approximated by finite-difference expressions constructed
from neighboring points in the **k** mesh. These are closely related to the
ddk wave functions used in RF calculations in section 2 and introduced in
[[lesson:rf1#5|tutorial DFPT1, section 5]]. We will use [[berryopt]] = -1,
which utilizes an improved coding of the calculation.
With the default [[rfdir]] = 1 1 1 all the cartesian components of the polarization are computed.

Now, run the *telast_4.abi* calculation, which should only take a minute or two, and
edit *telast_4.abo*. To calculate the elastic constants, we need to find the
stresses  sigma(1 1) and sigma(3 3). We see that each of the four datasets
have stress results, but that there are slight differences between those from,
for example dataset 1 and dataset 2, which should be identical. Despite our
tight limit, this is still a convergence issue. Look at the following convergence results,

    Dataset 1:
     At SCF step   13       vres2   =  6.68E-21 < tolvrs=  1.00E-18 =>converged.

    Dataset 2:
     At SCF step    1       vres2   =  5.31E-22 < tolvrs=  1.00E-18 =>converged.

Since dataset 2 has better convergence, we will use it.  Coherently, 
we will use also the dataset 4 results, choosing those in GPa units,

    -Cartesian components of stress tensor (GPa)         [Pressure=  9.1716E-03 GPa]
    - sigma(1 1)= -1.39848817E-03  sigma(3 2)=  0.00000000E+00
    - sigma(2 2)= -1.39848817E-03  sigma(3 1)=  0.00000000E+00
    - sigma(3 3)= -2.47179326E-02  sigma(2 1)=  0.00000000E+00

    -Cartesian components of stress tensor (GPa)         [Pressure= -1.1941E-02 GPa]
    - sigma(1 1)=  5.57418718E-03  sigma(3 2)=  0.00000000E+00
    - sigma(2 2)=  5.57418718E-03  sigma(3 1)=  0.00000000E+00
    - sigma(3 3)=  2.46737283E-02  sigma(2 1)=  0.00000000E+00

Let us now compute the numerical derivative of sigma(3 3) and compare to our
RF result. Recalling that our dimensionless strains were ±0.0001, we find
246.9583 GPa. This compares very well with the value 246.9602 GPa, the 3,3
element of the Rigid-ion elastic tensor we found from our anaddb calculation
in section 3. (Recall that our strains and stresses were both 3,3 or z,z or Voigt
3.) Similarly, the numerical derivative of  sigma(3 1) is 34.8634 GPa
compared to 34.8658 GPa, the 3,1 elastic-tensor element computed above.

The good agreement we found from this simple numerical differentiation
required that we had accurately relaxed the lattice so that the stress of the
unstrained structure was very small. Similar numerical-derivative comparisons
for systems with finite stress are more complicated, as discussed in
[notes by A. R. Oganov](/theory/documents/elasticity-oganov.pdf).
Numerical-derivative comparisons for the relaxed-ion results are extremely challenging
since they require relaxing atomic forces to exceedingly small limits.

Now let us examine the electric polarizations found in datasets 2 and 4,
focusing on the C/m^2 results,

           Polarization    -3.695942387E-15 C/m^2
           Polarization    -8.071946282E-16 C/m^2
           Polarization    -3.244653656E-01 C/m^2

           Polarization     5.396151684E-16 C/m^2
           Polarization     9.218443808E-18 C/m^2
           Polarization    -3.246052331E-01 C/m^2

While not labeled as such, these are the Cartesian x, y, and z components,
respectively, and the x and y components are zero within numerical accuracy as
they must be from symmetry. Numerical differentiation of the z component
yields -0.699337 C/m$^2$. This is to be compared with the z,3 element of our
rigid-ion piezoelectric tensor from section 3, -0.696149 C/m$^2$, and the two
results do not compare as well as we might wish.

What is wrong? There are two possibilities. The first is that the RF
calculation produces the proper piezoelectric tensor, while numerical
differentiation of the polarization produces the improper piezoelectric
tensor. This is a subtle point, for which you are referred to  D. Vanderbilt,
J. Phys. Chem. Solids 61, 147 (2000) [[cite:Vanderbilt2000]]. The improper-to-proper transformation
only effects certain tensor elements, however, and for our particular
combination of crystal symmetry and choice of strain there is no correction.
The second possibility is the subject of the next section.

## 5 Alternative DFPT calculation of some piezoelectric constants

Our GS calculation of the polarization in section 4 used, in effect, a finite-
difference approximation to ddk wave functions, while our RF calculations in
section 2 used analytic results based on the RF approach. Since the **k** grid
determined by [[ngkpt]] = 4 4 4 and [[nshiftk]] = 1 is rather coarse, this is
a probable source of discrepancy. Since this issue was noted previously in
connection with the calculation of Born effective charges by  Na Sai, K. M.
Rabe, and D. Vanderbilt, Phys. Rev. B 66, 104108 (2002) [[cite:Sai2002]], Abinit has
incorporated the ability to use finite-difference ddk wave functions from GS
calculations in RF calculations of electric-field-related 2DTE's. 

Copy
*telast_5.abi* into *Work_elast*, and edit *telast_5.abi*.

{% dialog tests/tutorespfn/Input/telast_5.abi %}

You should compare this with our previous RF data, *telast_2.abi*, and note that
dataset1 and the Common data (after entering relaxed structural results) are
essentially identical. Dataset 2 has been replaced by a non-self-consistent GS
calculation with [[berryopt]] = -2 specified to perform the finite-difference
ddk wave function calculation. (The finite-difference first-order wave
functions are implicit but not actually calculated in the GS polarization
calculation.) We have restricted [[rfdir]] to 0 0 1 since we are only
interested in the 3,3 piezoelectric constant. Now compare dataset 3 with that
in *telast_2.abi*. Can you figure out what we have dropped and why? Run the
*telast_5.abi* calculation, which will only take about a minute with our simplifications.

Now edit *telast_5.abo*, looking for the piezoelectric tensor,

      Rigid-atom proper piezoelectric tensor, in cartesian coordinates,
      (from strain response)
         j1       j2             matrix element
      dir pert dir pert     real part    imaginary part

       3    6   3    7        -0.0122230317         0.0000000000

{% dialog tests/tutorespfn/Refs/telast_5.abo %}

We immediately see a problem -- this output, like most of the .out file, is in
atomic units, while we computed our numerical derivative in conventional C/m$^2$
units. While you might think to simply run anaddb to do the conversion as
before, its present version is not happy with such an incomplete DDB file as
telast_5 has generated and will not produce the desired result. While it
should be left as an exercise to the student to dig the conversion factor out
of the literature, or better yet out of the source code, we will cheat and
tell you that 1 a.u.=57.2147606 C/m$^2$. Thus the new RF result for the 3,3 rigid-
ion piezoelectric constant is -0.699338  C/m$^2$ compared to the result found in
section 4 by a completely-GS finite difference calculation, -0.699337 C/m$^2$.
The agreement is now excellent!

The fully RF calculation in section 2 in fact will converge much more rapidly
with **k** sample than the partial-finite-difference method introduced here.
Is it worthwhile to have learned how to do this? We believe that is always
pays to have alternative ways to test results, and besides, this didn't take
much time. (Have you found the conversion factor on your own yet?)

A final word about AlP calculations: in this study, we have not used converged parameters for the sake of speed.
It is always the duty of the user to check the precision of his numerical calculations with respect to the
parameters of the calculation.


## 6 Response-function calculation of the elastic constants of Al metal

For metals, the existence of partially occupied bands is a complicating
feature for RF as well as GS calculations.
Now would be a good time to review [tutorial 4](/tutorial/base4) which dealt in detail with the interplay between
**k**-sample convergence and Fermi-surface broadening, especially [[lesson:base4#3|section 3 of tutorial 4]].
You should copy *telast_6.abi* into *Work_elast*, and begin your run
while you read on, since it involves a convergence study with multiple datasets and may take about two minutes.

{% dialog tests/tutorespfn/Input/telast_6.abi %}

While the run is in progress, edit *telast_6.abi*. As in *tbase4_3.abi*, we will set
[[udtset]] to specify a double loop.  In the present case, however, the outer
loop will be over 3 successively larger meshes of **k** points, while the
inner loop will be successively

  1. GS self-consistent runs with optimization of acell.
  2. GS density-generating run for the next step.
  3. Non-self-consistent GS run to converge unoccupied or slightly-occupied bands.
  4. RF run for symmetry-inequivalent elastic constants.

In Section 1, we did a separate GS structural optimization run and
transferred the results by hand to RF run section 2.  Because we are doing a
convergence test here, we have combined these steps, and use [[getcell]] to
transfer the optimized coordinates from the first dataset of the inner loop
forward to the rest.  If we were doing a more complicated structure with
internal coordinates that were also optimized, we would need to use both this
and [[getxred]] to transfer these, as in telast_1.abi.

The specific data for inner-loop dataset 1 is very similar to that for *telast_1.abi*.
Inner-loop dataset 2 is a bit of a hack.  We need the density
for inner-loop dataset 3, and while we could set [[prtden]] = 1 in dataset 1,
this would produce a separate density file for every step in the structural
optimization, and it isn't clear how to automatically pick out the last one.
So, dataset 2 picks up the wave functions from dataset 1 (only one file of
these is produced, for the optimized structure), does one more iteration with
fixed geometry, and writes a density file.

Inner-loop dataset 3 is a non-self-consistent run whose purpose is to ensure
that all the wave functions specified by [[nband]] are well converged. For
metals, we have to specify enough bands to make sure that the Fermi surface is
properly calculated.  Bands above the  Fermi level which have small occupancy
or near-zero occupancy if their energies exceed the Fermi energy by more than
a few times [[tsmear]], will have very little effect on the self-consistent
potential, so the [[tolvrs]] test in dataset 1 doesn't ensure their
convergence.  Using [[tolwfr]] in inner-loop dataset 3 does.  Partially-
occupied or unoccupied bands up to [[nband]]   play a different role in
constructing the first-order wave functions than do the many unoccupied bands
beyond [[nband]] which aren't explicitly treated in Abinit, as discussed in
S. de Gironcoli, Phys. Rev. B 51, 6773 (1995) [[cite:DeGironcoli1995]].  By setting [[nband]] exactly
equal to the number of occupied bands for RF calculations for semiconductors
and insulators, we avoid having to deal with the issue of converging
unoccupied bands.  Could we avoid the extra steps by simply using [[tolwfr]]
instead of [[tolvrs]] in dataset 1?  Perhaps, but experience has shown that
this does not necessarily lead to as well-converged a potential, and it is not
recommended.  These same considerations apply to phonon calculations for
metals, or in particular to [[qpt]]= 0 0 0 phonon calculations for the
interatomic force constants needed to find atom-relaxation contributions to
the elastic constants for non-trivial structures as in section 2 and section 3.

The data specific to the elastic-tensor RF calculation in inner-loop dataset 4
should by now be familiar.  We take advantage of the fact that for cubic
symmetry the only symmetry-inequivalent elastic constants are C$_{11}$, C$_{12}$, and
C$_{44}$.  Abinit, unfortunately, does not do this analysis automatically, so we
specify [[rfdir]] = 1 0 0 to avoid duplicate calculations.  (Note that if atom
relaxation is to be taken into account  for a more complex structure, the full
set of directions must be used.)

When the telast_6 calculations finish, first look at *telast_6.log* as usual to
make sure they have run to completion without error. Next, it would be a good
idea to look at the band occupancies occ?? (where ?? is a dual-loop dataset
index) reported at the end (following  `==END DATASET(S)==`).  The highest band,
the fourth in this case, should have zero or very small occupation, or you
need to increase [[nband]] or decrease [[tsmear]] .  Now, use your newly
perfected knowledge of the Abinit perturbation indexing conventions to scan
through *telast_6.abo* and find C$_{11}$ , C$_{12}$ , and C$_{44}$ for each of the three
**k** -sample choices, which will be  under the " Rigid-atom elastic tensor"
heading.  Also find the lattice constants for each case, whose convergence you
studied in [tutorial 4](/tutorial/base4).
You should be able to cut-and-paste these into a table like the following,

                C_11        C_12        C_44        acell
    ngkpt=3*6   0.003844    0.002294    0.001377    7.585323
    ngkpt=3*8   0.004409    0.002088    0.001355    7.583261 
    ngkpt=3*10  0.004392    0.002092    0.001354    7.583710

We can immediately see that the lattice constant converges considerably more
rapidly with **k** sample than the elastic constants.   For [[ngkpt]] =3*6,
acell is converged to 0.02%, while the C's have up to 15% error.  For [[ngkpt]]
=3*8, the C's are converged to better than 0.5%, even for the largest,
C$_{11}$, which should be acceptable.

As in [tutorial 4](/tutorial/base4), the [[ngkpt]] convergence is controlled by [[tsmear]].  The
smaller the broadening, the denser the **k** sample that is needed to get a
smooth variation of occupancy, and presumably stress, with strain.   While we
will not explore [[tsmear]] convergence in this tutorial, you may wish to do so
on your own.  

Also, even more than for the lattice parameter, the type of smearing
function plays an important role. The preferred smearing, [[occopt]]=7,
apparently performs worse than the the standard Fermi-Dirac broadening
[[occopt]]=3 that we have used above. Indeed, with [[occopt]]=7
and the same [[tsmear]]=0.02, one obtains:


                C_11        C_12        C_44        acell
    ngkpt=3*6   0.003113    0.002598    0.001575    7.543949
    ngkpt=3*8   0.004447    0.001948    0.001441    7.542962
    ngkpt=3*10  0.004781    0.001816    0.001233    7.544365
    ngkpt=3*12  0.004257    0.002108    0.001244    7.545547
    ngkpt=3*14  0.004056    0.002210    0.001300    7.545783
    ngkpt=3*16  0.004230    0.002120    0.001316    7.545453
    ngkpt=3*20  0.004271    0.002100    0.001311    7.545363


The reasons that this supposedly superior smoothing function performs poorly in this context
has not been investigated. Of course, [[tsmear]] has a different meaning
in both cases, and perhaps the value of 0.02Ha (=6315 Kelvin) for [[occopt]]=3 might yield a large error
anyhow (see the predicted [[acell]]).

The main thing to be learned is that checking
convergence with respect to all relevant parameters is **always** the user's
responsibility.   Simple systems that include the main physical features of a
complex system of interest will usually suffice for this testing.  Don't get
caught publishing a result that another researcher refutes on convergence
grounds, and do not blame such a mistake on Abinit!

Now we make a comparison with experiment.  Using the above computed values, converting the C's to standard
units (Ha/Bohr$^3$ = 2.94210119E+04 GPa) and using zero-temperature extrapolated
experimental results from P. M. Sutton, Phys. Rev. 91, 816 (1953) [[cite:Sutton1953]], we find

                                      C_11(GPa)  C_12(GPa)  C_44(GPa)
    Calculated (T=6315K)                129.2      61.5       39.8 
    Calculated (T=0K, Gaussian 0.02Ha)  125.6      61.8       38.6
    Experiment (T=0K)                   123.0      70.8       30.9


Is this good agreement?  The numerical values are not yet definitive,
as one should do more convergence studies anyhow. This being said, there isn't much literature on DFT calculations of
full sets of elastic constants.  Many calculations of the bulk modulus
(K=(C$_{11}$+2C$_{12}$ )/3 in the cubic case) typically are within 10% of experiment
for the LDA.  Running telast_6 with ixc=11, the Perdew-Burke-Enzerhof GGA,
increases the calculated C's by 1-2%, and wouldn't be expected to make a large
difference for a nearly-free-electron metal.

### Comment on symmetry

It is important to bear in mind that the way a tensor like the elastic tensor
appears is a function of the frame used. Thus for the aluminum fcc case
considered above, the nonzero elements are _C$_{11}$_, _C$_{12}$_, and _C$_{44}$_,
_provided that the crystal axes are aligned with the laboratory frame._
For an arbitrary alignment of the crystal axes, many more _C$_{ij}$_ elements will be
non-zero, and this can be confusing.

It is easy to see why this happens if you
imagine actually measuring the elastic tensor elements. If you start with the
conventional cubic cell, and apply pressure to one face, you can measure _C$_{11}$_.
 But if you turn the cell to some random angle, you'll measure a response
that is a mixture of _C$_{11}$_ and _C$_{12}$_.

Within ABINIT, if the aluminum fcc
cell is described using [[angdeg]] and [[acell]], then an axis of the
primitive cell will be aligned along the laboratory _z_ axis but this will not
lead to a (conventional) cell alignment with the laboratory frame. The
resulting elastic tensor will be correct but will appear to be more
complicated than in the illustration above. It can be rotated back to a simple
frame by hand (bearing in mind that all four indices of the fourth-rank
elastic tensor have to be rotated!) but it's easier to start with a more
conventional alignment of the unit cell.

If you use a standard text like
Bradley and Cracknell, The Mathematical Theory of Symmetry in Solids, Oxford [[cite:Bradley1972]]
you can find the standard primitive cell descriptions for the Bravais lattice
types and these are aligned as much as possible with a standard laboratory frame.
