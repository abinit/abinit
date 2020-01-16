---
authors: XG
---

# Tutorial on TDDFT

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
This is a rather simple s.ystem, with cylindrical symmetry,
allowing interesting understanding. Although we will suppose that you are
familiarized with quantum numbers for diatomic molecules, this should not play
an important role in the understanding of the way to use Abinit
implementation of Casida's formalism.

*Before beginning, you might consider to work in a different subdirectory as
for the other tutorials. Why not Work_tddft?*
Copy the files *ttddft_x.files* and *ttddft_1.in* in *Work_tddft*:

```sh
cd $ABI_TESTS/tutorial/Input
mkdir Work_tddft
cd Work_tddft
cp ../ttddft_x.files .
cp ../ttddft_1.in .
```

So, issue now:

    abinit < ttddft_x.files > log 2> err &

The computation is quite fast: about 15 secs on a 2.8 GHz PC.
Let's examine the input file *ttddft_1.in*.

{% dialog tests/tutorial/Input/ttddft_1.in %}

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
used. We have chosen the Perdew-Wang 92 LDA functional for both the self-
consistent and non-self-consistent calculations ([[ixc]] = 7).

We can now examine the output file *ttddft_1.out.*

{% dialog tests/tutorial/Refs/ttddft_1.out %}

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
       5->  6 3.10888E-01 8.45969E+00 -1.92741E+01 0.0000E+00 0.00E+00 0.00E+00 0.00E+00
       5->  7 3.10888E-01 8.45969E+00 -1.92741E+01 0.0000E+00 0.00E+00 0.00E+00 0.00E+00
       5->  8 3.44036E-01 9.36171E+00 -1.92409E+01 0.0000E+00 0.00E+00 0.00E+00 0.00E+00
       4->  6 3.64203E-01 9.91046E+00 -1.92207E+01 1.4463E-01 4.34E-01 0.00E+00 0.00E+00
       3->  6 3.64203E-01 9.91046E+00 -1.92207E+01 4.2299E-01 1.27E+00 0.00E+00 0.00E+00
       4->  7 3.64203E-01 9.91046E+00 -1.92207E+01 4.2299E-01 1.27E+00 0.00E+00 0.00E+00
       3->  7 3.64203E-01 9.91046E+00 -1.92207E+01 1.4463E-01 4.34E-01 0.00E+00 0.00E+00
       4->  8 3.97351E-01 1.08125E+01 -1.91876E+01 4.0028E-02 0.00E+00 1.20E-01 0.00E+00
       3->  8 3.97351E-01 1.08125E+01 -1.91876E+01 4.0028E-02 0.00E+00 0.00E+00 1.20E-01

Without the coupling matrix, these would be the excitation energies, for both
the spin-singlet and spin-triplet states. The coupling matrix modifies the
eigenenergies, by mixing different electronic excitations, and also lift some
degeneracies, e.g. the quadruplet formed by the combination of the degenerate
states 3-4 and 6-7 that gives the excitation energies with 3.64203E-01 Ha in the above table.

Indeed, concerning the spin-singlet, the following excitation energies are
obtained (see the next section of the output file):

      TDDFT singlet excitation energies (at most 20 of them are printed),
      and corresponding total energies.
      Excit#   (Ha)    and    (eV)    total energy (Ha)    major contributions
       1    3.47952E-01   9.46826E+00   -1.923699E+01    0.83(  5->  6)  0.17(  5->  7)
       2    3.48006E-01   9.46971E+00   -1.923693E+01    0.83(  5->  7)  0.17(  5->  6)
       3    3.62425E-01   9.86208E+00   -1.922251E+01    0.99(  5->  8)  0.00(  2-> 10)
       4    3.64202E-01   9.91043E+00   -1.922074E+01    0.37(  3->  7)  0.37(  4->  6)
       5    3.84223E-01   1.04553E+01   -1.920072E+01    0.37(  4->  6)  0.37(  3->  7)
       6    3.84236E-01   1.04556E+01   -1.920070E+01    0.37(  4->  7)  0.37(  3->  6)
       7    3.96699E-01   1.07947E+01   -1.918824E+01    0.99(  3->  8)  0.01(  4->  8)
       8    3.96723E-01   1.07954E+01   -1.918822E+01    0.99(  4->  8)  0.01(  3->  8)
       9    4.54145E-01   1.23579E+01   -1.913079E+01    1.00(  5->  9)  0.00(  3-> 12)
       ...

The excitation energies are numbered according to increasing energies, in Ha
as well as in eV. The total energy is also given (adding excitation energy to
the ground-state energy), and finally, the two major contributions to each
of these excitations are mentioned (size of the contribution then identification).

It is seen that the first and second excitations are degenerate (numerical
inaccuracies accounts for the meV difference), and mainly comes from the first
and second Kohn-Sham energy differences (between occupied state 5 and
unoccupied states 6 and 7). This is also true for the third excitation, that
comes from the third Kohn-Sham energy difference (between occupied state 5 and
unoccupied state 8). The quadruplet of Kohn-Sham energy differences, that was
observed at 3.64203E-01 Ha, has been split into one doublet and two singlets,
with numbers 4 (the lowest singlet), 5-6 (the doublet) while the last singlet
is not present in the 20 lowest excitations.

The list of oscillator strength is then provided.

      Oscillator strengths :  (elements smaller than 1.e-6 are set to zero)
      Excit#   (Ha)   Average    XX        YY        ZZ         XY        XZ        YZ
       1 3.47952E-01 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
       2 3.48006E-01 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
       3 3.62425E-01 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
       4 3.64202E-01 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
       5 3.84223E-01 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
       6 3.84236E-01 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
       7 3.96699E-01 5.759E-02 0.000E+00 1.928E-03 1.709E-01  0.00E+00  0.00E+00 -1.82E-02
       8 3.96723E-01 5.544E-02 0.000E+00 1.645E-01 1.855E-03  0.00E+00  0.00E+00  1.75E-02
       9 4.54145E-01 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
      10 4.60223E-01 9.496E-02 2.849E-01 0.000E+00 0.000E+00  0.00E+00  0.00E+00  0.00E+00
       ...

The first six transitions are forbidden, with zero oscillator strength. The
seventh and eighth transitions are allowed, with sizeable YY, YZ and ZZ components.

Next, one finds the excitation energies for the spin-triplet states:

      TDDFT triplet excitation energies (at most 20 of them are printed),
      and corresponding total energies.
      Excit#   (Ha)    and    (eV)    total energy (Ha)    major contributions
       1    2.88423E-01   7.84838E+00   -1.929652E+01    0.82(  5->  6)  0.18(  5->  7)
       2    2.88424E-01   7.84842E+00   -1.929652E+01    0.82(  5->  7)  0.18(  5->  6)
       3    2.99762E-01   8.15693E+00   -1.928518E+01    0.37(  3->  6)  0.37(  4->  7)
       4    3.33749E-01   9.08177E+00   -1.925119E+01    0.37(  4->  6)  0.37(  3->  7)
       5    3.33809E-01   9.08339E+00   -1.925113E+01    0.37(  4->  7)  0.37(  3->  6)
       6    3.36922E-01   9.16812E+00   -1.924802E+01    1.00(  5->  8)  0.00(  2-> 10)
       7    3.64202E-01   9.91045E+00   -1.922074E+01    0.37(  3->  7)  0.37(  4->  6)
       8    3.90779E-01   1.06336E+01   -1.919416E+01    0.67(  3->  8)  0.27(  2->  6)
       9    3.90834E-01   1.06351E+01   -1.919411E+01    0.67(  4->  8)  0.27(  2->  7)
       ...

Spin-triplet energies are markedly lower than the corresponding spin-singlet
energies. Also, the highest singlet derived from the Kohn-Sham quadruplet is
now the excitation number 7. The oscillator strengths also follow. At this
stage, we are in position to compare with experimental data, and try to
improve the quality of our calculation.

To summarize our results, we obtain the following five lowest-lying spin-
singlet excitation energies, with corresponding quantum numbers (that we
derive from the knowledge of the Kohn-Sham states quantum numbers):

    9.47 eV   m=+1,-1  even parity (Pi_g state)
    9.86 eV   m=0      even parity (Sigma_g state)
    9.91 eV   m=0      odd parity  (Sigma_u state)
    10.46 eV  m=+2,-2  odd parity  (Delta_u state)
    10.79 eV  m=+1,-1  odd parity  (Pi_u state)

and the following five lowest-lying spin-triplet excitations energies, with
corresponding quantum numbers:

    7.85 eV   m=+1,-1  even parity (Pi_g state)
    8.16 eV   m=0      odd parity  (Sigma_u state)
    9.08 eV   m=+2,-2  odd parity  (Delta_u state)
    9.16 eV   m=0      even parity (Sigma_g state)
    9.91 eV   m=0      odd parity  (Sigma_u state)

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

  * The appearance of the spin-singlet $^1\Sigma_g$ state at 9.86 eV (Spin-singlet state 2)
  * The inversion between the spin-triplet $^3\Pi_g$ and $^3\Sigma_u$ states (Spin-triplet states 1 and 2)
  * The appearance of the spin-triplet $^3\Sigma_g$ state at 9.16 eV (Spin-triplet state 4)

Still, the agreement between these TDDFT values and the experimental values is
much better than anything that can be done on the sole basis of Kohn-Sham
energy differences, that are (for spin-singlet and -triplet):

    8.46 eV   m=+1,-1   even parity (Pi_g state)
    9.36 eV   m=0       odd parity  (Sigma_u state)
    9.91 eV   m=0(twice),+2,-2 odd parity  (Sigma_u and Delta_u states)
    10.81 eV  m=+1,-1   odd parity  (Pi_u state)

## Convergence studies

There are several parameters subject to convergence studies in this context:
the energy cut-off, the box size, and the number of unoccupied bands.

We will start with the number of unoccupied states. The only input parameter
to be changed in the input file is the value of nband2. The following results
are obtained, for nband2 = 12, 30, 60, 100 and 150 (Energies given in eV):

    Singlet 1 :  9.47   9.44   9.39   9.36   9.35
    Singlet 2 :  9.86   9.74   9.68   9.66   9.66
    Singlet 3 :  9.91   9.91   9.91   9.91   9.91
    Singlet 4 : 10.46  10.45  10.44  10.44  10.43
    Singlet 5 : 10.79  10.79  10.79  10.79  10.79
    Triplet 1 :  7.85   7.84   7.83   7.82   7.82
    Triplet 2 :  8.16   8.08   8.03   8.00   8.00
    Triplet 3 :  9.08   9.07   9.05   9.05   9.04
    Triplet 4 :  9.16   9.16   9.15   9.15   9.15
    Triplet 5 :  9.91   9.91   9.91   9.91   9.91

You might try to obtain one of these...

The computation with nband2 = 100 takes
about 7 minutes on a 2.8 GHz PC, and gives a result likely converged within
0.01 eV. Let's have a look at these data. Unfortunately, none of the above-
mentioned discrepancies with experimental data is resolved, although the
difference between the first and second spin-triplet states decreases
significantly. Although we see that at least 60 bands are needed to obtain
results converged within 0.05 eV, we will continue to rely on 12 bands to try
to understand the most important discrepancies, while keeping the CPU time to a low level.

We next try to increase the cut-off energy. Again, this is fairly easy. One
can e.g. set up a double dataset loop. The following results are obtained, for
ecut = 25, 35, 45, 55, 65, and 75 Ha:

    Singlet 1 :  9.47   9.41   9.39   9.36   9.36
    Singlet 2 :  9.86   9.83   9.78   9.76   9.76
    Singlet 3 :  9.91   9.97  10.01  10.02  10.03
    Singlet 4 : 10.46  10.37  10.32  10.30  10.29
    Singlet 5 : 10.79  10.87  10.90  10.91  10.92
    Triplet 1 :  7.85   7.79   7.76   7.74   7.73
    Triplet 2 :  8.16   8.02   7.94   7.92   7.91
    Triplet 3 :  9.08   8.98   8.92   8.90   8.89
    Triplet 4 :  9.16   9.28   9.33   9.34   9.34
    Triplet 5 :  9.91   9.83   9.78   9.77   9.76

You might try to obtain one of these...

The computation with ecut=75 takes
about 90 secs on a 2.8 GHz PC, and gives a result likely converged within 0.01
eV. Let us have a look at these data. Concerning the discrepancies with the
experimental results, we see that the position of the second spin-singlet
state has even worsened, the difference between the first and second spin-
triplet states decreases, so that, together with an increase of nband, their
order might become the correct one, and the fourth spin-triplet state energy
has increased, but not enough.

We finally examine the effect of the cell size. Again, this is fairly easy.
One can e.g. set up a double dataset loop. The following results are obtained,
for acell = (6 5 5), (7 6 6), (8 7 7), (9 8 8), (10 9 9) and (12 11 11):

    Singlet 1 :  9.47   9.37   9.33   9.33   9.33   9.33
    Singlet 2 :  9.86   9.78   9.84   9.91   9.96  10.03
    Singlet 3 :  9.91   9.88   9.85   9.85   9.85   9.85
    Singlet 4 : 10.46  10.41  10.38  10.37  10.37  10.37
    Singlet 5 : 10.79  10.98  11.14  11.27  11.19  11.04
    Triplet 1 :  7.85   7.75   7.72   7.71   7.72   7.72
    Triplet 2 :  8.16   8.18   8.18   8.18   8.19   8.20
    Triplet 3 :  9.08   9.07   9.06   9.06   9.06   9.06
    Triplet 4 :  9.16   9.36   9.55   9.68   9.78   9.90
    Triplet 5 :  9.91   9.88   9.85   9.85   9.85   9.85

Obviously, the cell size plays an important role in the spurious appearance of
the states, that was remarked when comparing against experimental data.
Indeed, the singlet 2 and triplet 4 states energy increases strongly with the
cell size, while all other states quickly stabilize (except the still higher singlet 5 state).

There is one lesson to be learned from that convergence study: the
convergence of different states can be quite different. Usually, converging
the lower excited states do not require too much effort, while it is quite
difficult, especially concerning the supercell size, to converge higher states.

At this stage, we will simply stop this convergence study, and give the
results of an ABINIT calculation using ecut 45 Hartree, acell 12 11 11, and 30
bands (not fully converged, though!), then compare the results with other
LDA/TDLDA results (from [[cite:Casida1998]]) and experimental results:

                          present Casida experimental
    Singlet Pi_g      :  9.25    9.05     9.31
    Singlet Sigma_u-  :  9.72    9.63     9.92
    Singlet Delta_u   : 10.22   10.22    10.27
    Triplet Sigma_u+  :  7.95    7.85     7.75
    Triplet Pi_g      :  7.64    7.54     8.04
    Triplet Delta_u   :  8.89    8.82     8.88
    Triplet Sigma_u-  :  9.72    9.63     9.67

Our calculation is based on pseudopotentials, while Casida's calculation is an
all-electron one. This fact might account for the 0.1-0.2 eV discrepancy
between both calculations (it is of course the user's responsibility to test
the influence of different pseudopotentials on his/her calculations). The
agreement with experimental data is on the order of 0.2 eV, with the exception
of the $^3\Pi_g$ state (0.4 eV). In particular, we note that LDA/TDLDA is
not able to get the correct ordering of the lower two triplet states ... One
of our problems was intrinsic to the LDA/TDLDA approximation ...

## The choice of the exchange-correlation potential

As emphasized in [[cite:Casida1998]], choosing a different functional for the self-consistent part
(XC potential) and the generation of the coupling matrix (XC
kernel) can give a better description of the higher-lying states. Indeed, a
potential with a -1/r tail (unlike the LDA or GGA) like the van Leeuwen-Baerends one [[cite:VanLeeuwen1994]],
can reproduce fairly well the ionisation energy, giving a much
better description of the Rydberg states. Still, the LDA kernel works pretty well.

In order to activate this procedure, set the value of [[ixc]] in dataset 1 to the
SCF functional, and the value of ixc in dataset 2 to the XC functional to be
used for the kernel. Use pseudopotentials that agree with the SCF functional.
