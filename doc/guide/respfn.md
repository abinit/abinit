---
authors: XG, DCA
---

# The DFPT (respfn) code

This page complements the main [[help:abinit]], for matters related
to responses to perturbations computed with DFPT.
It will be easier to discover the present file with the help of the [[tutorial:rf1|DFPT1 tutorial]].  

<a id="intro"></a> 
## 0 Introducing the computation of responses
  
ABINIT can compute the response to different perturbations, and provide access
to quantities that are second derivatives of total energy (2DTE) with respect
to these perturbations. 
Presently, they can be of four types: 

1. phonons 
2. static homogeneous electric field
3. strain  
4. magnetic field (coupling to the spin, not the orbital motion)

The physical properties connected to 2DTE with respect to perturbations (1) and (2) are the phonon
dynamical matrices, the dielectric tensor, and the Born effective charges,
while the additional strain perturbation (3), mixed with phonon and electric
field leads to elastic constant, internal strain, and piezoelectricity.
The magnetic field perturbation is a recent addition to ABINIT, and will not be detailed at present.


More functionalities of the computation of responses should be implemented
sooner or later. Some third derivatives of the
total energy (3DTE) are also implemented. The 3DTE might give phonon-phonon
coupling, non-linear electric response, anharmonic elastic constants, Gruneisen parameters,...

The basic quantities that ABINIT will compute are the **first-order** derivatives
of the wavefunctions (1WF) with respect to these perturbations. The later
calculation of the 2DTE and 3DTE from these 1WF is an easy computational task:
the construction of the 2DTE with respect to perturbations j1 and j2
involves mainly evaluating matrix elements between the 1WF of j1 and/or the 1WF of j2. 
More details on this technique can be found in [[cite:Gonze1997]] and [[cite:Gonze1997a]].

The calculation of the 1WF for a particular perturbation is done using a
variational principle and an algorithm rather similar to the one used to find
the unperturbed ground-state wavefunctions. Thus, a lot of technical details
and parameters are the same for both ground-state and response-function
calculations. This justifies the development of one unique code for these two
classes of properties: many of the routines of abinit are common in these
calculations, or had to be duplicated, but with relatively small modifications.

The ABINIT code performs a rather primitive analysis of the calculated 2DTEs.
For example, it gives the phonon frequencies, electronic dielectric tensor and
effective charges. But the main output of the code is the Derivative DataBase
(DDB): a file that contains the set of all 2DTEs and 3DTEs calculated by the
code. This DDB can be manipulated by the MRGDDB code, and fully analyzed by
the Anaddb code. See the corresponding [[help:mrgddb]] and [[help:anaddb]].

<a id="1"></a>
## 1 Description of perturbations
  
The perturbation of the **phonon** type is the displacement of one atom along
one of the axis of the unit cell, by a unit of length (in reduced coordinates).
It is characterized by two integer numbers and one wavevector.
The two integer numbers are the number of the moved atom, which will be noted
**ipert**, and the number of the axis of the unit cell, which will be noted **idir**. 

!!! important

    *ipert* for phonon perturbation can have values between 1 and [[natom]],
    *idir* can have values between 1 and 3.

The set of all possible phonon perturbations for one wavevector has thus 3 * [[natom]] elements. 
From this basis set, all phonons can be constructed by linear combination. 
The set of atoms to be moved in one dataset of ABINIT is determined by the input
variable [[rfatpol]]. The set of directions to be considered in one dataset of
ABINIT is determined by the input variable [[rfdir]]. The wavevector to be
considered in one dataset of ABINIT is determined by the input variables
[[nqpt]], [[qpt]], and [[qptnrm]].

The perturbations of the **electric field** type are

  * the application of the homogeneous electric field along the axes of the reciprocal lattice 
  * the derivative of the Hamiltonian with respect to the electronic wavevector along 
    the axes of the reciprocal lattice (which allows to compute derivatives of wavefunctions with respect to their wavevector), 
    an **auxiliary** quantity needed **before** the application of the homogeneous electric field. 
    The perturbation is the change of wavevector by dk in the Hamiltonian, hence the perturbation 
    is referred to as the derivative dk perturbation ("ddk" perturbation). 

Note 1
:   the ddk perturbation is defined as the derivative with respect to k
    in reduced coordinates; this is equivalent to applying a linear perturbation
    of strength $2\pi$ along the conjugate direction in real space. This statement
    comes from the derivation of the phase factor $\exp(i2\pi kr)$ with respect to
    k in reduced coordinates.

Note 2
:   in case the electric field type perturbation is computed inside a
    finite electric field, the derivative of the Hamiltonian with respect to the
    electronic wavevector is not computed: everything is done by a finite-
    difference technique. Also, the definition of the homogeneous electric field
    perturbation is not along the axes of the reciprocal lattice, but in cartesian
    coordinates. Sorry for the possible confusion ...

These electric-type perturbations are also characterized by two numbers:
*ipert* being natom + 1 for the ddk perturbation and natom + 2 for the electric
field, and *idir* being 1, 2 or 3, as for phonon perturbations. Although the
possibility of electric field characterized by a non-zero wavevector is
envisioned for a future version of the code, at present only homogeneous
fields are considered. So the wavevector of the electric field type
perturbations is $\Gamma$ (q=0).

The perturbations of the **strain** type are either an uniaxial strain or a
shear strain. The strain perturbations are considered in cartesian coordinates
(x,y,z). They are characterized by two numbers, with *ipert* being natom + 3 for
the uniaxial strains, and natom + 4 for the shear strains, and *idir* describes
the particular component. 

Explicitly, for uniaxial strains:

* idir = 1 gives the xx strain perturbation, 
* idir = 2 gives the yy strain perturbation, 
* idir = 3 gives the zz strain perturbation, 

while for shear strains:

* idir=1 gives the yz strain perturbation, 
* idir=2 gives the xz perturbation, 
* idir=3 gives the xy perturbation.  

Note that the "rigid-atom" elastic constants, as output of ABINIT, are those
obtained with **frozen** internal coordinates. The internal coordinate
relaxation, needed to give "physical" elastic constants can be handled through
the knowledge of the internal strain and dynamical matrix at $\Gamma$, in ANADDB.
Of course, if all the internal coordinate are fixed by symmetry, all the
internal strains vanish, and the "rigid-atom" and "physical" elastic constants are equal.  
Limitations of the present implementation (as of v5.7):

  * Symmetry is presently used to skip redundant k points in the BZ sum, 
    but not to skip redundant strain perturbations.

We also define the index of the perturbation, called *pertcase*, equal to idir + 3*ipert. 
Accordingly, pertcase runs from 1 to 3 * (natom + 4), and will be
needed to identify output and input files, see section 6.

!!! summary

    To summarize, the perturbations are characterized by two numbers, **ipert** from
    1 to natom + 4, and **idir**, from 1 to 3, as well as one wavevector (that is
    gamma when a non-phonon perturbation is considered). A number called
    **pertcase** combines *ipert* and *idir*, and runs from 1 to 3 * (natom + 4).

The 2DTE, being derivative of the total energy with respect to two
perturbations, will be characterized by two sets of (idir,ipert), or by two
pertcase numbers, while 3DTE will need three such sets or pertcase numbers.
In addition they will depend on one wavevector (for 2DTE) or two wavevectors (for 3DTE).

In the non-stationary implementation of the 2DTE, used for off-diagonal elements in ABINIT, the first pertcase
corresponds to the perturbation that gives the derivative of the potential,
and the second pertcase corresponds to the perturbation that gives the
derivative of the wavefunctions.

<a id="2"></a>
## 2 Filenames and input of ground-state wavefunctions
  
The **same** 'files' file as for GS calculations is used for RF calculations.
Actually, in the multi-dataset mode, one will be able to make in one ABINIT
run, ground-state computations as well as response-function computations, so
that the 'files' file must be the same.... The 'input' file will have many
common input variables for these different cases, but also some separate ones.

Two ground-state wavefunction files might be needed:

  * the file of ground-state wavefunctions at a set of wavevectors, named k-points
  * the file of ground-state wavefunctions at the corresponding k+q, where q is the wavevector of the perturbation

These files also contain the corresponding eigenvalues.

Note that the second file (k+q) will be identical to the first one (k), in the
case of zero-wavevector perturbation. Also, if the q-wavevector maps the
k-point grid onto itself (q being the difference between points that belong to
the grid), all the information on the k+q grid is contained in the k grid, a
fact that ABINIT is able to exploit.

One should have a look at the input variables [[irdwfk]], and [[irdwfq]], for
a description of the ground-state wavefunction file names generated from the
root names provided in the 'files' file. In the multi-dataset mode, the
following input variables will be relevant: [[getwfk]], and [[getwfq]]. The
file names of the ground-state wavefunction file follow the same convention as
for the ground-state case. Thus, read the 
[[help:abinit#files-file|corresponding section]] of the abinit help file, if needed.

In the case of an electric field perturbation, the output 1WF of the
corresponding ddk perturbation is needed as input. If the option [[rfelfd]]=1
is used, then the code will take care of doing first the derivative dk
perturbation calculation, then write the 1WF at the correct place, as an
output file, then begin the homogeneous field perturbation calculation.
Usually, the use of [[rfelfd]]=1 is not recommended, as the ddk computation is
the most often done with different parameters as the electric field perturbation,
a dataset with [[rfelfd]]=2 being followed with a dataset with [[rfelfd]]=3.

The nomenclature for first-order wavefunction files is also given in the
[[help:abinit#files-file|abinit help]] file, but it is worth to specify it in more
detail here. The root name is formed from the string of character in the third
line of the 'files' file (for an input file) or the fourth line of the 'files'
file (for an output file), that is complemented, in the multi-dataset mode, by
' **_DS** ' and the dataset number, and then,
the string ' **_1WF** ' is added, followed by pertcase (the index of the perturbation). This gives, e.g., for
a 'files' file whose fourth line is '/tmp/o', for the dataset number 3, and a
perturbation corresponding to the displacement of the second atom in the x
direction (pertcase=4), the following name of the corresponding 1st-order
wavefunction output file:
    
     /tmp/o_DS3_1WF4

Such a file might be used as input of another computation, or of another dataset.

The relevant input variables are: [[ird1wf]], and [[irdddk]], as well as
[[get1wf]], and [[getddk]], in the multi-dataset mode.

<a id="symmetries"></a>
## 3 The use of symmetries
  
In order to understand correctly the behaviour of response-function runs, 
some information on the use of symmetries must be given.

Some perturbations (including their wavevector) may be invariant for some
symmetries. The code is able to use symmetries to skip perturbations of which a
symmetric has already been calculated (except in the case of strain
perturbations). ABINIT is also able to use the symmetries that keeps the
perturbations invariant, to reduce the number of k points needed for the
sampling of electronic wavefunctions in the Brillouin zone (although this
feature is not optimal yet). There is one exception to this, the ddk
perturbation, for which the spatial symmetries cannot be used yet.

In any case, unlike for the ground-state, the input k-point set for response
function should NOT have been decreased by using spatial symmetries, prior to
the loop over perturbations (see section 4). Only the time-reversal symmetry,
retained by calculations at q=0, ought to be used to decrease this input
k-point set. Considering each perturbation in turn, ABINIT will be able to
select from this non-reduced set of k-points, the proper k-point set,
automatically, by using the symmetries that leave each perturbation invariant.

Accordingly, the preferred way to generate the k-point grid is of course to
use the [[ngkpt]] or [[kptrlatt]] input variables, with different values of [[kptopt]]:

  * kptopt = 1 for the ground state
  * kptopt = 2 for response functions at q=0 
  * kptopt = 3 for response functions at non-zero q 

<a id="4"></a>
## 4 Organisation of response-function computations
  
In agreement with the information provided in the previous sections, different
cases can be distinguished.

When one considers the response to an atomic displacement with q=0, the
following procedure is suggested:

  * first, a self-consistent ground-state computation with the restricted set of k-points 
    in the Irreducible Brillouin Zone (with [[kptopt]]=1)

  * second, a self-consistent response-function computation with the atomic displacement perturbation, 
    with the half set of k-points (with [[kptopt]]=2)

When one considers the response to an electric field (with q=0), the following
procedure is suggested:

  * first, a self-consistent ground-state computation with the restricted set of k-points 
    in the Irreducible Brillouin Zone (with [[kptopt]]=1)

  * second, a non-self-consistent response-function computation of the d/dk perturbation, 
    with the half set of k-points (with [[kptopt]]=2, and [[iscf]]=-3)

  * third, a self-consistent response-function computation of the electric field perturbation, 
    with the half set of k-points (with [[kptopt]]=2)

When one considers the response to an atomic displacement in the special case
where q connects k-points that both belong to the special k-point grid, the
following procedure is suggested:

  * first, a self-consistent ground-state computation with the restricted set of k-points 
    in the Irreducible Brillouin Zone (with [[kptopt]]=1)

  * second, a self-consistent response-function computation of the atomic displacement perturbation, 
    with the full set of k-points (with [[kptopt]]=3)

When one considers the response to an atomic displacement for a general q
point, the following procedure is suggested:

  * first, a self-consistent ground-state computation with the restricted set of k-points 
    in the Irreducible Brillouin Zone (with [[kptopt]]=1)

  * second, a non-self-consistent ground-state run with the set of k+q points, that might be 
    reduced thanks to symmetries (with [[kptopt]]=1)

  * third, a self-consistent response-function computation of the atomic displacement perturbation, 
    with the full set of k-points (with [[kptopt]]=3)

Of course, these different steps can be combined when a set of responses is
looked for. In particular, the computations of responses at gamma, in the case
where the full dynamical matrix as well as the dielectric tensor and the Born
effective charges are needed, can be combined as follows:

  * first, a self-consistent ground-state computation with the restricted set of k-points 
    in the Irreducible Brillouin Zone (with [[kptopt]]=1)

  * second, the three non-self-consistent response-function computations (one for each direction) 
    of the d/dk perturbation, with the half set of k-points (with [[kptopt]]=2, and [[iscf]]=-3)

  * third, all the self-consistent response-function computations of the electric field perturbations 
    and of the atomic displacements, with the half set of k-points (with [[kptopt]]=2)

Still, computations of perturbations at different q wavevectors cannot be
mixed. But they can follow the other computations. Supposing that
perturbations at q=0 and a general q point are to be performed, they will be combined as follows:

  * first, a self-consistent ground-state computation with the restricted set of k-points 
    in the Irreducible Brillouin Zone (with [[kptopt]]=1)

  * second, the three non-self-consistent response-function computations (one for each direction) 
    of the d/dk perturbation, with the half set of k-points (with [[kptopt]]=2, and [[iscf]]=-3)

  * third, all q=0 self-consistent response-function computations of the electric field perturbations 
    and of the atomic displacements, with the half set of k-points (with [[kptopt]]=2)

  * fourth, a non-self-consistent ground-state computation with the set of k+q points, 
    that might be reduced thanks to symmetries (with [[kptopt]]=1)

  * fifth, the self-consistent response-function computations of the atomic displacement perturbations 
    with a q wavevector, with the full set of k-points (with [[kptopt]]=3)

Note that the error in the 2DTE is **linear** in the **ground-state**
wavefunction error (unlike the error due to the 1WFs). Moreover, a large
prefactor is associated with this source of error (it can even cause cause the
instability of the SCF procedure). As a consequence, the convergence of the
ground-state wavefunction should be very good. The same is true at the level
of the ddk wavefunctions.

As mentioned in the introduction, inside the response-function part of the
code, the calculation proceeds in two steps: first the calculation of the
first-order derivative of the wavefunctions (1WF), then the combinations of
these 1WF to build the 2DTE and 3DTE.

In an initialisation part, the input file and all the ground-state files are
read, and a few basic quantities are constructed.

In the first part, every perturbation is examined, one at a time, separately:

  * A file containing previous RF wavefunctions is eventually read.

  * The minimisation of the variational expression is performed, and this procedure generates 
    the 1WF as well as the first-order density and potential.

  * It is possible, knowing these quantities for the perturbation ipert1, to construct all the 2DTE 
    with respect to this perturbation and any ipert2, except for ipert2 being an homogeneous electric field, 
    in which case the derivative of the ground-state wavefunctions with respect to their wavevector (ddk) is also needed. 
    This feature has been implemented for ipert2 being any phonon (of the same wavevector than ipert1), 
    or an homogeneous electric field.

  * The first-order wavefunctions (1WF) are written as well as the first-order density or potential (if requested).

<a id="5"></a>
## 5 List of relevant input variables
  
A subset of the ABINIT input variables have a modified meaning or a modified
behaviour in case of RF calculations. Here is the list of these input
variables, together with the variables that applies only to RF computations.
Note that the code will do a RF calculation ([[optdriver]]=1) when one of
[[rfphon]] or [[rfelfd]] is non-zero.

  * [[amu]]
  * [[getwfk]], [[getwfq]], [[get1wf]], [[getddk]] 
  * [[irdwfk]], [[irdwfq]], [[ird1wf]], [[irdddk]] 
  * [[iscf]]
  * [[nkpt]]
  * [[nqpt]], [[qpt]], [[qptnrm]] 
  * [[nsym]]
  * [[rfasr]]
  * [[rfatpol]]
  * [[rfdir]]
  * [[rfelfd]]
  * [[rfphon]]
  * [[dfpt_sciss]]
  * [[tolwfr]], [[toldfe]], [[tolvrs]]

<a id="output"></a>
## 6 The different output files
  
Output from the code goes to several places listed below.

**6.1. The log file**

This file is the same as the log file of the abinit code when computing ground
state (GS) results. Possibly, the output of datasets related to response
functions will be intertwined with those concerned with ground-state case. The
purpose of this file is the same as in the GS case, and the use of error messages is unchanged.

<a id="6.2"></a>
**6.2. The main output file**

This file is the same as the main output file of the abinit code when
computing ground state (GS) results. Possibly, the output of datasets
related to response functions will be intertwined with those concerned with
ground-state case. We explain here the parts related to the RF computation.

The initialisation part is the same as for the GS. So, the reader is advised
to read the [[help:abinit#outputfile|section 6.2]] of the abinit help file,
as well as the first paragraph of the section [[help:abinit#6.3|6.3]] of
this file. Afterwards, the content of the main output file differs a bit...

The main output file reports on the initialisation of the ground-state
wavefunctions at k+q, then the loop on the perturbations begins. For each
perturbation, there is:

  * a short description of the perturbation
  * the timing information
  * the report on the initialisation of the 1WF
  * then the iterations for the minimisation begin, and the output file describes 
    the number of the iteration, the second derivative of the total energy obtained (2DTEnergy in Ha), 
    the change in 2DTEnergy since last iteration, the maximum residual over all bands and k points, 
    and the square of the potential residual.
  * after the iterations are completed, the residuals are reported
  * in case of the derivative dk perturbation, ek2 (to be explained) and the f-sum rule ratio value are printed. 
    The f-sum rule ratio should be close to 1 (not when ecutsm/=0, however).
  * then the components of the 2DTEnergy, broken in at most 14 pieces, depending on the perturbation
  * then the 2DTEnergy in Hartree and in eV
  * then the relaxation contribution (caused by changes in wavefunctions), and the non-relaxation contribution 
    (Ewald and frozen-wavefunction contribution)
  * then the 2DTEnergy evaluated using a non-variational expression.

After the computation of each perturbation, the output file reports on the
2DTE matrix elements. This part is not executed if the only perturbation is
the derivative dk perturbation. It will give the following information:

  * if [[prtvol]]=1 or bigger, the full detail of every contribution to the 2DTE in reduced coordinates.
  * the 2DTE in reduced coordinates.
  * then it will use the 2DTE to perform already some analysis of the data without use the Mrgddb and Anaddb codes, namely: 
    the full dynamical matrix (cartesian coordinates, masses included) the effective charges, and the dielectric tensor, 
    the phonon frequencies (including the analysis of the non- analyticity if we are at $\Gamma$). 
    Note that phonon frequencies lower than 1.0d-8Ha (absolute value) are automatically set to zero, 
    while imaginary phonon frequencies (square roots of negative eigenvalues - indicating an instability) 
    are printed as negative (this facilitates the post-processing).

For this last analysis, the code assumes that the whole set of perturbations
in a class has been calculated, either all the phonon perturbations or the
homogeneous electric field perturbation, or both. If this is not true, the
code will give results that may be wrong in the case that the reduced system
of coordinates is not cartesian (for the dynamical matrix, the effective
charge tensor of the dielectric matrix), or in all case wrong (the phonon
frequencies); also the code will put zero in the matrix elements that have not been calculated. 
A Warning message is issued if the above information cannot be trusted.

Finally, the code provides the timing information.

<a id="6.3"></a>
**6.3. The first-order wavefunction (1WF) files**

These are unformatted data files containing the planewaves coefficients of all
the wavefunctions, written in the following format. First, the header (see
[[help:abinit#header|section 6.4]] of the abinit help file), followed by
    
```fortran
     bantot=0                                    <-- counts over all bands
     index=0                                     <-- index for the wavefunctions
     do isppol=1,nsppol
      do ikpt=1,nkpt
       write(unit) npw1,nspinor,nband                    <-- for each k point
       write(unit) kg(1:3,1:npw1)                        <-- plane wave reduced coordinates
       do iband=1,nband1
        write(unit) (eigen1(jband+(iband-1)*nband+bantot),jband=1,2*nband)  <-- column of eigenvalue matrix
        write(unit) (cg(ii+index),ii=1,2*npw1*nspinor)   <-- wavefunction coefficients
       enddo                                            for a single band and k point
       bantot=bantot+nband
       index=index+2*npw1*nspinor*nband1
      enddo
     enddo
```
  
In this code section, npw1(ikpt) is the number of planewaves in the basis at
the k+q point, nband1(ikpt) is likewise the number of bands at the k point
(which may vary among k points depending on occopt), and the factor of 2 in
writing the wavefunction results from the fact that it is complex.  
eigen1 is the array that contains the matrix element of the first-order
Hamiltonian between the different ground-state wavefunctions. It could be used
to build the electron-phonon coupling and deformation potentials.  
Note that the format for first-order WF file differs from the format used for
the ground-state WF file by the fact that eigen1 is now an array, and no more
a vector, and is written with the corresponding wf, and no more before the
writing of all wf for one k point.

**6.4. The first-order density files**

They consist of the header lines, followed by
    
```fortran
    write(unit) (rhor1(ir),ir=1,cplex*ngfft(1)*ngfft(2)*ngfft(3))
```

Here, rhor1 is the change of electron density in electrons/Bohr^3. The
parameter cplex is 1 when q=0 and 2 when q/=0 . Indeed, for q=0, the density
change is a real quantity, while it is complex in general when q/=0 .

<a id="ddb"></a>
**6.5. The derivative database (DDB)**

It is made of two parts. The first should allow one to unambiguously identify
the run that has generated the DDB, while the second part contains the 2DTE,
grouped by blocks of data.

Note that the DDB output of the ABINIT code can be merged with other DDBs as
described in the [[help:mrgddb|Mrgddb help file]].

The first part contains:

  * the DDB version number (that defines the structure of the DDB)
  * seven parameters needed for the dimensionning of the DDB file 
    ([[natom]], [[nkpt]], [[nsppol]], [[nsym]], [[ntypat]], [[occopt]], and [[nband]] - 
    or the array [[nband]] ([[nkpt]]* [[nsppol]]) if [[occopt]]=2)
  * different information on the run that generated the 2DTE 
    ([[acell]],[[amu]],[[ecut]],[[iscf]],[[ixc]],[[kpt]],[[kptnrm]], 
     [[ngfft]],[[occ]],[[rprim]],[[dfpt_sciss]],[[symrel]],[[xred]],[[tnons]],[[typat]],[[tolwfr]],[[wtk]],[[ziontypat]], 
    as well as information on the pseudopotentials by means of their Kleinman-Bylander energies). 
    These values are simply a transcription of the input data, or other simple internal parameters.

Note: the format and content of this first part of the DDBs have to be updated in the future ...

The second part contains:

  * the number of data blocks
  * for each data block, the type of the block, its number of elements, and the list of elements.

The elements of a 2DTE block are described by 4 integers and 2 real numbers.
The 2 first integers define a first perturbation in the form (idir1,ipert1),
the two next define a second perturbation in the form (idir2,ipert2). The
matrix element corresponds to the derivative of the total energy with respect
to the parameters corresponding to these perturbations. The real numbers are
the real and imaginary parts of the 2DTE. Sometimes, the code uses spatial
symmetries, the time-reversal symmetry, or even the permutation of first and
second perturbations to deduce the value of non-computed matrix elements. This
behaviour might be improved, as it is sometimes confusing ...

<a id="numerical-quality"></a>
## 7 Numerical quality of the calculations
  
It is possible to get from the RF calculations essentially EXACT derivatives
of the total energy with respect to perturbations. There is a published
account of this fact in [[cite:Gonze1995]]. An agreement of 8 digits
or more was obtained in the comparison with finite-difference derivatives of GS data.

The accuracy of these calculations are thus entirely determined by the input
parameters the user choose for the RF run, and the preparatory GS runs.

We will now review the convergence parameters, usually the same as for the GS
calculations, indicating only the specific features related to RF calculations.

Input parameters that could influence the accuracy of the calculation are:

  * [[ecut]] (the energy cut-off, that depends strongly on the pseudopotential)
  * [[ixc]] (describing the exchange-correlation functional)
  * [[nkpt]](or, more accurately, the Brillouin zone sampling, that can be determined 
     alternatively by the inputs variables [[ngkpt]] or [[kptrlatt]])
  * one of the self-consistent convergence tolerance parameters, [[toldfe]], [[tolvrs]], or [[tolwfr]].

The input parameters [[ecut]], [[ecutsm]], [[ixc]], [[intxc]], the set of
k-wavevectors, as well as the related variables have to be the SAME in the
ground-state calculations that go before a RF run and this RF run.

Namely: do not try to use ground-state wavefunction files produced with
[[ecut]]=10Ha for a RF run with [[ecut]]=20Ha. In some cases, the code will complain
and stop, but in other cases, it might simply produce garbage !

If the value of [[ngfft]] is input by hand, its value must also be equal in
the GS and RF cases. ALWAYS use [[ngfft]] large enough to have boxcut=2 or
larger, in order to avoid any FFT filter error. In the GS case, boxcut as
small as 1.5 could be allowed in some cases. It is not allowed with RF
calculations, because they are more sensitive to that error.

The convergence tests with respect to [[ecut]], and the k-point grid should be
done carefully. This was already emphasized in the [[help:abinit]] and is re-
emphasized here. The user should test the convergence DIRECTLY on the property
he or she is interested in. For example, if the user wants a phonon frequency
accurate to 10 cm^-1, he/she could be lead to do a full calculation (GS+RF) of
phonons at 30Ha, then another full calculation at 35Ha, then another at
40Ha... It is an error to rely on tolerance on the total energy (for example
1mHa/atom) or geometry (accuracy of one part per thousand on the bond lengths)
to draw 'a priori' conclusions on the convergence of other quantities, and not
monitor the convergence of these directly. To be clear: if phonon frequencies
are needed, check the convergence of phonon frequencies !

The user should note that for bands with very small occupancy in the metallic
case as well as unoccupied bands for insulators, the ground state run
preceeding response function runs will not necessarily converge these
wavefunctions using usual ground-state tests such as [[toldfe]] or (better)
[[tolvrs]]. To be sure that inaccuracies are not introduced into the response
function calculation by poorly converged unoccupied bands, a separate run
starting from a saved charge density ([[prtden]]=1 in the self-consistent run)
and using [[iscf]]=-2 and [[tolwfr]] may be needed.
