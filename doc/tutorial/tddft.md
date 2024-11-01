---
authors: XG
---

# TDDFT (Casida)

## Time-Dependent Density Functional Theory, Casida's approach.

This tutorial aims at showing how to get the following physical properties, for finite systems:

  * Excitation energies
  * Associated oscillator strengths
  * Frequency-dependent polarizability and optical spectra

in the Casida approach, within Time-Dependent Density Functional Theory.

This tutorial should take about 30 minutes.

[TUTORIAL_README]

## Brief theoretical introduction

In order to perform Time-Dependent Density Functional Theory calculations (TDDFT)
of electronic excitations and oscillator strengths, in the Casida's approach,
you should first have some theoretical background.
TDDFT was first developed in the eighties, but the direct calculation of
electronic excitations was introduced much later, by Casida and co-workers.
A comprehensive description of the underlying formalism is given in

    "Time-Dependent Density Functional Response Theory of Molecular Systems:
     Theory, Computational Methods, and Functionals"
     M. E. Casida
     in Recent Development and Applications of Modern Density Functional Theory,
     edited by J.M. Seminario (Elsevier, Amsterdam, 1996), p. 391.
     <http://dx.doi.org/10.1016/S1380-7323(96)80093-8>

However this reference might be hard to get, that is why we have based the
tutorial instead on the following (also early) papers:
[[cite:Casida1998]], [[cite:Casida1998a]], and [[cite:Vasiliev1998]].

The first of these papers, [[cite:Casida1998]], will be used as main reference for our tutorial.

From these papers, you will learn that a TDDFT calculation of electronic
excitation energies start first from a usual ground-state calculation, with a
chosen XC functional. Such a calculation produces a spectrum of Kohn-Sham
electronic energies. It is widely known that differences between occupied and
unoccupied Kohn-Sham electronic energies resemble excitation energies (the
difference in energy between an excited state and the ground state), although
there is no real theoretical justification for this similarity.

These differences between Kohn-Sham electronic energies are the starting point
of Casida's approach: in the framework of TDDFT, their square give the main
contribution to the diagonal part of a matrix, whose eigenvalues will be the
square of the excitation energies. One has to add to the diagonal matrix made
from the squares of Kohn-Sham energy differences, a coupling matrix, whose
elements are four-wavefunction integrals of the Coulomb and exchange-correlation kernel.
The exchange-correlation kernel contribution will differ
in the spin-singlet and in the spin-triplet states, this being the only
difference between spin-singlet and spin-triplet states. See Eqs.(1.3) and
(1.4) of [[cite:Casida1998a]], and Eqs.(1-2) of [[cite:Vasiliev1998]].

The construction of the coupling matrix can be done on the basis of an
exchange-correlation kernel that is derived from the exchange-correlation
functional used for the ground-state, but this is not a requirement of the
theory, since such a correspondence only holds for the exact functional. In
practice, the approximation to the XC potential and the one to the XC kernel
are often different. See section III of [[cite:Casida1998]].

A big drawback of the currently known XC potentials and XC kernels is observed
when the system is infinite in at least one direction (e.g. polymers, slabs,
or solids). In this case, the addition of the coupling matrix is unable to
shift the edges of the Kohn-Sham band structure (each four-wavefunction
integral becomes too small). There is only a redistribution of the oscillator
strengths. In particular, the DFT band gap problem is NOT solved by TDDFT.
Also, the Casida's approach relies on the discreteness of the Kohn-Sham spectrum.

Thus, the TDDFT approach to electronic excitation energies in ABINIT is ONLY
valid for finite systems (atoms, molecules, clusters). Actually, only one
k-point can be used, and a "box center" must be defined, close to the center
of gravity of the system.

The Casida formalism also gives access to the oscillator strengths, needed to
obtain the frequency-dependent polarizability, and corresponding optical
spectrum. In the ABINIT implementation, the oscillators strengths are given as
a second-rank tensor, in cartesian coordinates, as well as the average over
all directions usually used for molecules and clusters. It is left to the user
to generate the polarisability spectrum, according to e.g. Eq.(1.2) of [[cite:Casida1998a]].

One can also combine the ground state total energy with the electronic
excitation energies to obtain Born-Oppenheimer potential energy curves for
excited states. This is illustrated for formaldehyde in [[cite:Casida1998a]].

Given its simplicity, and the relatively modest CPU cost of this type of
calculation, Casida's approach enjoys a wide popularity. There has been
hundreds of papers published on the basis of methodology. Still, its accuracy
might be below the expectations, as you will see. As often, test this method
to see if it suits your needs, and read the recent literature ...

## A first computation of electronic excitation energies and oscillator strengths, for N$_2$

We will now compute and analyse the excitation energies of the diatomic molecule N$_2$.
This is a rather simple system, with cylindrical symmetry,
allowing interesting understanding. Although we will suppose that you are
familiarized with quantum numbers for diatomic molecules, this should not play
an important role in the understanding of the way to use Abinit
implementation of Casida's formalism.

*Before beginning, you might consider to work in a different subdirectory as
for the other tutorials. Why not Work_tddft?*
Copy the file *ttddft_1.abi* in *Work_tddft*:

```sh
cd $ABI_TESTS/tutorial/Input
mkdir Work_tddft
cd Work_tddft
cp ../ttddft_1.abi .
```

So, issue now:

    abinit ttddft_1.abi > log &

The computation is quite fast: about 3 secs on a 2.8 GHz PC.
Let's examine the input file *ttddft_1.abi*.

{% dialog tests/tutorial/Input/ttddft_1.abi %}

There are two datasets: the first one corresponds to a typical ground-state
calculation, with only occupied bands. The density and wavefunctions are
written, for use in the second data set. The second dataset is the one where
the TDDFT calculation is done. Moreover, the non-self-consistent calculation
of the occupied eigenfunctions and corresponding eigenenergies is also
accomplished. This is obtained by setting [[iscf]] to -1.

Please, take now some time to read the information about this value of [[iscf]], and the few
input variables that acquire some meaning in this context (namely,
[[boxcenter]], [[td_mexcit]], and [[td_maxene]]). Actually, this is most of
the information that should be known to use the TDDFT in ABINIT!

You will note that we have 5 occupied bands (defined for dataset 1), and that
we add 7 unoccupied bands in the dataset 2, to obtain a total of 12 bands. The
box is not very large (6x5x5 Angstrom), the cutoff is quite reasonable, 25
Hartree), and as requested for the Casida's formalism, only one k point is
used. We are using the Perdew-Wang 92 LDA functional for both the self-
consistent and non-self-consistent calculations ([[ixc]] = -1012), as deduced
by abinit by looking at the pseudopotential file..

We can now examine the output file *ttddft_1.abo.*

{% dialog tests/tutorial/Refs/ttddft_1.abo %}

One can jump to the second dataset section, and skip a few non-interesting
information, in order to reach the following information:

     *** TDDFT : computation of excited states ***
     Splitting of  12 bands in   5 occupied bands, and   7 unoccupied bands,
     giving    35 excitations.

The matrix that is diagonalized, in the Casida's formalism, is thus a 35x35 matrix.
It will give 35 excitation energies.
Then, follows the list of excitation energies, obtained from the difference of
Kohn-Sham eigenvalues (occupied and unoccupied), for further reference. They
are ordered by increasing energy. In order to analyze the TDDFT as well as
experimental data in the next section, let us mention that the Kohn-Sham
eigenfunctions in this simulation have the following characteristics:

  * The first and fifth states are (non-degenerate) occupied $\sigma$ states (m=0), with even parity
  * The second state is a (non-degenerate) occupied $\sigma$ state (m=0), with odd parity
  * The third and fourth states are degenerate occupied $\pi$ states (m=+1,-1), with odd parity
  * The sixth and seventh states are degenerate unoccupied $\pi$ states (m=+1,-1), with even parity
  * The state 8 is a (non-degenerate) unoccupied $\sigma$ state (m=0), with even parity

Combining states 3,4 and 5 with 6, 7 and 8, give the first nine Kohn-Sham energy differences:

      Transition  (Ha)  and   (eV)   Tot. Ene. (Ha)  Aver     XX       YY       ZZ
       5->  6 3.06494E-01 8.34013E+00 -2.03684E+01 0.0000E+00 0.00E+00 0.00E+00 0.00E+00
       5->  7 3.06494E-01 8.34013E+00 -2.03684E+01 0.0000E+00 0.00E+00 0.00E+00 0.00E+00
       5->  8 3.50400E-01 9.53488E+00 -2.03245E+01 0.0000E+00 0.00E+00 0.00E+00 0.00E+00
       4->  6 3.60026E-01 9.79682E+00 -2.03149E+01 5.5534E-01 1.67E+00 0.00E+00 0.00E+00
       3->  6 3.60026E-01 9.79682E+00 -2.03149E+01 3.7821E-03 1.13E-02 0.00E+00 0.00E+00
       4->  7 3.60026E-01 9.79682E+00 -2.03149E+01 3.7822E-03 1.13E-02 0.00E+00 0.00E+00
       3->  7 3.60026E-01 9.79682E+00 -2.03149E+01 5.5534E-01 1.67E+00 0.00E+00 0.00E+00
       4->  8 4.03933E-01 1.09916E+01 -2.02710E+01 3.7976E-02 0.00E+00 1.14E-01 0.00E+00
       3->  8 4.03933E-01 1.09916E+01 -2.02710E+01 3.7976E-02 0.00E+00 0.00E+00 1.14E-01

Without the coupling matrix, these would be the excitation energies, for both
the spin-singlet and spin-triplet states. The coupling matrix modifies the
eigenenergies, by mixing different electronic excitations, and also lift some
degeneracies, e.g. the quadruplet formed by the combination of the degenerate
states 3-4 and 6-7 that gives the (four-fold degenerate) excitation energies 
with 3.60026E-01 Ha in the above table.

Indeed, concerning the spin-singlet, the following excitation energies are
obtained (see the next section of the output file):

      TDDFT singlet excitation energies (at most 20 of them are printed),
      and corresponding total energies.
      Excit#   (Ha)    and    (eV)    total energy (Ha)    major contributions
       1    3.45362E-01   9.39779E+00   -2.032952E+01    1.00(  5->  7)  0.00(  1->  7)
       2    3.45434E-01   9.39975E+00   -2.032945E+01    1.00(  5->  6)  0.00(  1->  6)
       3    3.60026E-01   9.79681E+00   -2.031486E+01    0.50(  3->  6)  0.50(  4->  7)
       4    3.68693E-01   1.00326E+01   -2.030619E+01    0.99(  5->  8)  0.00(  2-> 10)
       5    3.83765E-01   1.04428E+01   -2.029112E+01    0.50(  4->  7)  0.50(  3->  6)
       6    3.83798E-01   1.04437E+01   -2.029108E+01    0.50(  4->  6)  0.50(  3->  7)
       7    4.03285E-01   1.09740E+01   -2.027160E+01    0.99(  3->  8)  0.01(  4->  8)
       8    4.03304E-01   1.09745E+01   -2.027158E+01    0.99(  4->  8)  0.01(  3->  8)
       9    4.59051E-01   1.24914E+01   -2.021583E+01    0.91(  2->  8)  0.04(  3->  7)
       ...

The excitation energies are numbered according to increasing energies, in Ha
as well as in eV. The total energy is also given (adding excitation energy to
the ground-state energy), and finally, the two major contributions to each
of these excitations are mentioned (size of the contribution then identification).

It is seen that the first and second excitations are degenerate (numerical
inaccuracies accounts for the meV difference), and mainly comes from the first
and second Kohn-Sham energy differences (between occupied state 5 and
unoccupied states 6 and 7). This is also true for the fourth excitation, that
comes from the third Kohn-Sham energy difference (between occupied state 5 and
unoccupied state 8). The quadruplet of Kohn-Sham energy differences, that was
observed at 3.60026E-01 Ha, has been split into one doublet and two singlets,
with numbers 3 (the lowest singlet), 5-6 (the doublet) while the last singlet
is not present in the 20 lowest excitations.

The list of oscillator strength is then provided.

      Oscillator strengths :  (elements smaller than 1.e-6 are set to zero)
      Excit#   (Ha)   Average    XX        YY        ZZ         XY        XZ        YZ
       1 3.45362E-01 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
       2 3.45434E-01 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
       3 3.60026E-01 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
       4 3.68693E-01 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
       5 3.83765E-01 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
       6 3.83798E-01 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
       7 4.03285E-01 5.881E-02 0.000E+00 1.409E-03 1.750E-01  0.00E+00  0.00E+00 -1.57E-02
       8 4.03304E-01 5.685E-02 0.000E+00 1.692E-01 1.361E-03  0.00E+00  0.00E+00  1.52E-02
       9 4.59051E-01 8.613E-02 2.584E-01 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
      10 4.61447E-01 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
       ...

The first six transitions are forbidden, with zero oscillator strength. The
seventh and eighth transitions are allowed, with sizeable YY, YZ and ZZ components.

Next, one finds the excitation energies for the spin-triplet states:

      TDDFT triplet excitation energies (at most 20 of them are printed),
      and corresponding total energies.
      Excit#   (Ha)    and    (eV)    total energy (Ha)    major contributions
       1    2.84779E-01   7.74923E+00   -2.039010E+01    1.00(  5->  7)  0.00(  5->  6)
       2    2.84781E-01   7.74928E+00   -2.039010E+01    1.00(  5->  6)  0.00(  5->  7)
       3    2.97188E-01   8.08689E+00   -2.037770E+01    0.50(  3->  7)  0.50(  4->  6)
       4    3.30296E-01   8.98780E+00   -2.034459E+01    0.50(  4->  7)  0.50(  3->  6)
       5    3.30344E-01   8.98913E+00   -2.034454E+01    0.50(  4->  6)  0.50(  3->  7)
       6    3.43692E-01   9.35234E+00   -2.033119E+01    1.00(  5->  8)  0.00(  2-> 10)
       7    3.60026E-01   9.79681E+00   -2.031486E+01    0.50(  3->  6)  0.50(  4->  7)
       8    3.85278E-01   1.04839E+01   -2.028961E+01    0.91(  2->  7)  0.09(  3->  8)
       9    3.85310E-01   1.04848E+01   -2.028957E+01    0.91(  2->  6)  0.08(  4->  8)
       ...

Spin-triplet energies are markedly lower than the corresponding spin-singlet
energies. Also, the highest singlet derived from the Kohn-Sham quadruplet is
now the excitation number 7. The oscillator strengths also follow. At this
stage, we are in position to compare with experimental data, and try to
improve the quality of our calculation.

To summarize our results, we obtain the following five lowest-lying spin-
singlet excitation energies, with corresponding quantum numbers (that we
derive from the knowledge of the Kohn-Sham states quantum numbers):

    9.40 eV   m=+1,-1  even parity (Pi_g state)
    9.80 eV   m=0      odd parity  (Sigma_u state)
    10.03 eV  m=0      even parity (Sigma_g state)
    10.44 eV  m=+2,-2  odd parity  (Delta_u state)
    10.97 eV  m=+1,-1  odd parity  (Pi_u state)

and the following five lowest-lying spin-triplet excitations energies, with
corresponding quantum numbers:

    7.75 eV   m=+1,-1  even parity (Pi_g state)
    8.09 eV   m=0      odd parity  (Sigma_u state)
    8.99 eV   m=+2,-2  odd parity  (Delta_u state)
    9.35 eV   m=0      even parity (Sigma_g state)
    9.80 eV   m=0      odd parity  (Sigma_u state)

The quantum number related to the effect of a mirror plane, needed for $\Sigma$
states, could not be attributed on the sole basis of the knowledge of Kohn-
Sham orbitals quantum numbers.

The lowest-lying experimental spin-singlet excitation energies, see table III
of [[cite:Casida1998]], are as follows:

    9.31 eV   m=+1,-1  even parity (Pi_g state)
    9.92 eV   m=0      odd parity  (Sigma_u- state)
    10.27 eV  m=+2,-2  odd parity  (Delta_u state)

and the lowest-lying experimental spin-triplet excitations energies are:

    7.75 eV   m=0      odd parity  (Sigma_u+ state)
    8.04 eV   m=+1,-1  even parity (Pi_g state)
    8.88 eV   m=+2,-2  odd parity  (Delta_u state)
    9.67 eV   m=0      odd parity  (Sigma_u- state)

In several cases, the agreement is quite satisfactory, on the order of 0.1-0.2
eV. However, there are also noticeable discrepancies. Indeed, we have to understand, in our simulation:

  * The appearance of the spin-singlet $^1\Sigma_g$ state at 10.03 eV (Spin-singlet state 3)
  * The inversion between the spin-triplet $^3\Pi_g$ and $^3\Sigma_u$ states (Spin-triplet states 1 and 2)
  * The appearance of the spin-triplet $^3\Sigma_g$ state at 9.35 eV (Spin-triplet state 4)

Still, the agreement between these TDDFT values and the experimental values is
much better than anything that can be done on the sole basis of Kohn-Sham
energy differences, that are (for spin-singlet and -triplet):

    8.34 eV   m=+1,-1   even parity (Pi_g state)
    9.53 eV   m=0       odd parity  (Sigma_u state)
    9.80 eV   m=0(twice),+2,-2 odd parity  (Sigma_u and Delta_u states)
    10.99 eV  m=+1,-1   odd parity  (Pi_u state)

## Convergence studies

There are several parameters subject to convergence studies in this context:
the energy cut-off, the box size, and the number of unoccupied bands.

We will start with the number of unoccupied states. The only input parameter
to be changed in the input file is the value of nband2. The following results
are obtained, for nband2 = 12, 30, 60, 100 and 150 (Energies given in eV):

    Singlet 1 :  9.40   9.37   9.33   9.30   9.28
    Singlet 2 :  9.80   9.80   9.80   9.79   9.80
    Singlet 3 : 10.03   9.91   9.85   9.83   9.83
    Singlet 4 : 10.44  10.43  10.43  10.42  10.42
    Singlet 5 : 10.97  10.97  10.97  10.97  10.97
    Triplet 1 :  7.75   7.75   7.73   7.73   7.72
    Triplet 2 :  8.06   8.02   7.98   7.96   7.95
    Triplet 3 :  8.99   8.97   8.96   8.96   8.95
    Triplet 4 :  9.35   9.34   9.34   9.34   9.34
    Triplet 5 :  9.80   9.80   9.80   9.80   9.80

You might try to obtain one of these...

The computation with nband2 = 100 takes a bit more than 1 minute,
and gives a result likely converged within 0.01 eV.
Let us have a look at these data. Unfortunately, none of the above-
mentioned discrepancies with experimental data is resolved, although the
difference between the first and second spin-triplet states decreases
significantly. Although we see that at least 60 bands are needed to obtain
results converged within 0.05 eV, we will continue to rely on 12 bands to try
to understand the most important discrepancies, while keeping the CPU time to a low level.

We next try to increase the cut-off energy. Again, this is fairly easy. One
can e.g. set up a double dataset loop. The following results are obtained, for
ecut = 25, 30, 35 and 45 Ha:

    Singlet 1 :  9.40   9.37   9.36   9.36 
    Singlet 2 :  9.80   9.78   9.77   9.77 
    Singlet 3 : 10.03  10.03  10.04  10.04 
    Singlet 4 : 10.46  10.37  10.32  10.30
    Singlet 5 : 10.97  10.98  10.98  10.98
    Triplet 1 :  7.75   7.72   7.72   7.72 
    Triplet 2 :  8.09   8.06   8.06   8.06
    Triplet 3 :  8.99   8.97   8.96   8.96
    Triplet 4 :  9.35   9.35   9.35   9.35
    Triplet 5 :  9.80   9.78   9.77   9.77  

You might try to obtain one of these...

The computation with ecut=30 takes a couple of seconds
and gives a result likely converged within 0.01 eV. 
The modifications with respect to the results with ecut=25 Ha are quite small.

We finally examine the effect of the cell size. Again, this is fairly easy.
e keep 12 bands, but stick to ecut=30 Ha.
One can e.g. set up a double dataset loop. The following results are obtained,
for acell = (6 5 5), (8 7 7), (10 9 9), (12 11 11), (14 13 13), (16 15 15) and (20 19 19):

    Singlet 1 :  9.37   9.25   9.24   9.24   9.25   9.25   9.25
    Singlet 2 :  9.78   9.73   9.72   9.72   9.72   9.72   9.72
    Singlet 3 : 10.03  10.04  10.18  10.26  10.29  10.32  10.34
    Singlet 4 : 10.42  10.35  10.34  10.34  10.34  10.34  10.35
    Singlet 5 : 10.98  11.34  11.40  11.12  10.95  10.84  10.70
    Triplet 1 :  7.72   7.60   7.60   7.60   7.60   7.60   7.60
    Triplet 2 :  8.06   8.07   8.07   8.08   8.08   8.08   8.08
    Triplet 3 :  8.97   8.94   8.94   8.94   8.94   8.94   8.94
    Triplet 4 :  9.35   9.73   9.72   9.72   9.72   9.72   9.72
    Triplet 5 :  9.78   9.75  10.00  10.12  10.18  10.22  10.27

Obviously, the cell size plays an important role in the spurious appearance of
the states, that was remarked when comparing against experimental data.
Indeed, the singlet 3 and triplet 4 states energy increases strongly with the
cell size, while all other states quickly stabilize (except the still higher singlet 5 state).
Actually, the singlet 3 and 4 switch between (16 15 15) and (20 19 19),
and the triplet 4 and 5 switch between (6 5 5) and (8 7 7)  ...
The presence of the singlet 3 state (Sigma_g) in the three lowest ones,
an the triplet 4 state (Sigma_g) were two of our discrepancies.

There is one lesson to be learned from that convergence study: the
convergence of different states can be quite different. Usually, converging
the lower excited states do not require too much effort, while it is quite
difficult, especially concerning the supercell size, to converge higher states.

At this stage, we will simply stop this convergence study, and give the
results of an ABINIT calculation using ecut 30 Hartree, acell 10 9 9, and 100
bands, focusing only on those states already converged, then compare the results with other
LDA/TDLDA results (from [[cite:Casida1998]]) and experimental results:

                       present  Casida experimental
    Singlet Pi_g      :  9.20    9.05     9.31
    Singlet Sigma_u-  :  9.72    9.63     9.92
    Singlet Delta_u   : 10.33   10.22    10.27
    Triplet Sigma_u+  :  7.99    7.85     7.75
    Triplet Pi_g      :  7.59    7.54     8.04
    Triplet Delta_u   :  8.92    8.82     8.88
    Triplet Sigma_u-  :  9.72    9.63     9.67

Our calculation is based on pseudopotentials, while Casida's calculation is an
all-electron one. This fact might account for the discrepancy
between both calculations (maximal 0.15 eV - it is of course the user's responsibility to test
the influence of different pseudopotentials on his/her calculations). The
agreement with experimental data is on the order of 0.2 eV, with the exception
of the $^3\Pi_g$ state (0.45 eV). In particular, we note that LDA/TDLDA is
not able to get the correct ordering of the lowest two triplet states ... One
of our problems was intrinsic to the LDA/TDLDA approximation ...

## The choice of the exchange-correlation potential and kernel

As emphasized in [[cite:Casida1998]], choosing a different functional for the self-consistent part
(XC potential) and the generation of the coupling matrix (XC kernel) 
can give a better description of the higher-lying states. 
Indeed, a potential with a -1/r tail (unlike the LDA or GGA) like the van Leeuwen-Baerends one [[cite:VanLeeuwen1994]],
can reproduce fairly well the ionisation energy, giving a much
better description of the Rydberg states. Still, the LDA kernel works pretty well.

In order to activate this procedure, set the value of [[ixc]] in dataset 1 to the
SCF functional, and the value of ixc in dataset 2 to the XC functional to be
used for the kernel. Use pseudopotentials that agree with the SCF functional in dataset 1.
As of writing (end of 2020), the ABINIT implementation has a strong restriction on the XC kernels that can be used. They
must be of LDA-type. See the list of allowed [[ixc]] values in the description of [[iscf]]=-1.
