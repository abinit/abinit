---
authors: SPesant, MCote, XG, BAmadon
---

# DFT+U

## The projected density of states of NiO.

This tutorial aims at showing how to perform a DFT+U calculation using Abinit (see also [[cite:Amadon2008a]])

You will learn what is a DFT+U calculation and what are the main input
variables controlling this type of calculation.

It is supposed that you already know how to do PAW calculations using ABINIT.
Please follow the two tutorials on PAW in ABINIT ([PAW1](/tutorial/paw1), [PAW2](/tutorial/paw2)), if this is not the case.

This tutorial should take about 1 hour to complete.

[TUTORIAL_README]

## 0 Short summary of the DFT+U method

The standard Local Density Approximation (LDA), where the exchange and
correlation energy is fit to homogeneous electron gas results, is a functional
that works well for a vast number of compounds. But, for some crystals, the
interactions between electrons are so important that they cannot be
represented by the LDA alone. Generally, these highly correlated materials
contain rare-earth metals or transition metals, which have partially filled *d*
or *f* bands and thus localized electrons.

The LDA tends to delocalize electrons over the crystal, and each electron
feels an average of the Coulombic potential. For highly correlated materials,
the large Coulombic repulsion between localized electrons might not be well
represented by a functional such as the LDA. A way to avoid this problem is to
add a Hubbard-like, localised term, to the LDA density functional. This
approach is known as DFT+U (formerly referred to as LDA+U). In the actual implementation, we
separate localized d or f electrons, on which the Hubbard term will act, from
the delocalized ones (*s* and *p* electrons). The latter are correctly described
by the usual LDA calculation. In order to avoid the double counting of the
correlation part for localized electrons (already included in the LDA,
although in an average manner), another term - called the double-counting
correction - is subtracted from the Hamiltonian.

In Abinit, two double-counting corrections are currently implemented:

-The Full localized limit (FLL) [[cite:Liechtenstein1995]] ([[usepawu]]=1)

-The Around Mean Field (AMF) [[cite:Czyzyk1994]]  ([[usepawu]]=2)

For some systems, the result might depend on the choice of the double-counting method.
However, the two methods generally give similar results.

## 1 Ground state calculation of NiO using LDA

*Before continuing, you might consider to work in a different subdirectory as
for the other tutorials. Why not Work_dftu?
In what follows, the names of files will be mentioned as if you were in this subdirectory.*

Copy the file *tdftu_1.abi* from *\$ABI_TESTS/tutorial/Input* to your *Work_dftu* directory with:

```sh
cd $ABI_TESTS/tutorial/Input
mkdir Work_dftu
cd Work_dftu
cp ../tdftu_1.abi .
```

{% dialog tests/tutorial/Input/tdftu_1.abi %}

Now run the code as usual.
The job should take less than 20 seconds on a laptop. It calculates the LDA
ground state of the NiO crystal. A low cutoff and a small number of k-points
are used in order to speed up the calculation. During this time you can take a
look at the input file.

The NiO crystallizes in the rocksalt structure, with one Ni and one O atom in
the primitive cell (the crystallographic primitive cell). However, NiO is
known to exhibit an antiferromagnetic ordering at low temperature (along the
<111> direction). From the electronic point of view, the true unit cell has
two Ni and two O atoms: the local magnetic moment around the first Ni atom
will have a sign opposite to the one of the other Ni atom.

You should take some time to examine the values used for the input variables
[[xred]], [[rprim]] (note the last line!), [[typat]], [[spinat]], [[nsppol]],
and [[nspden]], that define this antiferromagnetic ordering along the <111>
direction (of a conventional cubic cell).

If you take a look at the output file (tdftu_1.out), you can see the
integrated total density in the PAW spheres (see the [PAW1](/tutorial/paw1)
and [PAW2](/tutorial/paw2) tutorials on PAW formalism). This value roughly
estimates the magnetic moment of NiO:

     Integrated electronic and magnetization densities in atomic spheres:
     ---------------------------------------------------------------------
     Radius=ratsph(iatom), smearing ratsm=  0.0000. Diff(up-dn)=approximate z local magnetic moment.
     Atom    Radius    up_density   dn_density  Total(up+dn)  Diff(up-dn)
        1   1.81432     8.564383     7.188016     15.752398     1.376367
        2   1.81432     7.188016     8.564383     15.752398    -1.376367
        3   1.41465     2.260902     2.260902      4.521804    -0.000000
        4   1.41465     2.260902     2.260902      4.521804     0.000000
     

The atoms in the output file, are listed as in the [[typat]] variable (the
first two are nickel atoms and the last two are oxygen atoms). The results
indicate that spins are located in each nickel atom of the doubled primitive
cell. Fortunately, the LDA succeeds to give an antiferromagnetic ground state
for the NiO. But the result does not agree with the experimental data.

The magnetic moment (the difference between up and down spin on the nickel atom)
range around 1.6-1.9 according to experiments  ([[cite:Cheetham1983]],[[cite:Neubeck1999]],[[cite:Sawatzky1984]],
[[cite:Hufner1984]])
Also, as the Fermi level is at 0.33748 Ha (see the *tdftu_1.abo* file), one
can see (on the *tdftu_1.o_EIG* file that contains eigenvalues for the three k-point of this calculation) that the band gap obtained between the last (24th) occupied band (0.31537 Ha, at k
point 3) and the first (25th) unoccupied band (0.35671 Ha, at kpoint 3) is
approximately 1.1 eV which is lower than the measured value of 4.0-4.3 eV
(This value could be modified using well-converged parameters but would still
be much lower than what is expected). A easier and graphical way to evaluate the gap would be to plot the density
of states (see last section of this tutorial).

Making abstraction of the effect of insufficiently convergence parameters, the
reason for the discrepancy between the DFT-LDA data and the experiments is
first the fact the DFT is a theory for the ground state and second, the lack
of correlation of the LDA. Alone, the homogeneous electron gas cannot
correctly represent the interactions among $d$ electrons of the Ni atom. That is
why we want to improve our functional, and be able to manage the strong correlation in NiO.

## 2 DFT+U with the FLL double-counting

As seen previously, the LDA does not gives good results for the magnetization
and band gap compared to experiments.
At this stage, we will try to improve the correspondence between calculation
and experimental data. First, we will use the DFT(LDA)+U with the Full
localized limit (FLL) double-counting method.

FLL and AMF double-counting expressions are given in the papers listed above,
and use the adequate number of electrons for each spin. For the Hubbard term,
the rotationally invariant interaction is used.

!!! note

    It is important to notice that in order to use DFT+U in Abinit, you must
    employ PAW pseudopotentials.

You should run abinit with the *tdftu_2.abi* input file. This calculation takes
also less than 20 seconds on a laptop.
During the calculation, you can take a look at the input file.

{% dialog tests/tutorial/Input/tdftu_2.abi %}

Some variable describing the DFT+U parameters have been added to the previous file. All
other parameters were kept constant from the preceding calculation. First, you
must set the variable [[usepawu]] to one (for the FLL method) and two (for the
AMF method) in order to enable the DFT+U calculation. Then, with [[lpawu]] you
give for each atomic species ([[znucl]]) the values of angular momentum (l) for
which the DFT+U correction will be applied. The choices are 1 for *p*-orbitals, 2 for *d*-orbitals
and 3 for *f*-orbitals. You cannot treat s orbitals with DFT+U in the
present version of ABINIT. Also, if you do not want to apply DFT+U correction
on a species, you can set the variable to -1. For the case of NiO, we put
[[lpawu]] to 2 for Ni and -1 for O.


!!! note

    The current implementation applies DFT+U correction only inside atomic sphere. To check if this
    approximation is realistic, relaunch the calculation with [[pawprtvol]] equal to three.
    Then search for ph0phiint in the log file:

        pawpuxinit: icount, ph0phiint(icount)= 1  0.90467

    This line indicates that the norm of atomic wavefunctions inside atomic sphere is 0.90, rather close
    to one. In the case of nickel, the approximation is thus realistic. The case where the norm is too small
    (close to 0.5) is discussed in [[cite:Geneste2017]].

Finally, as described in the article cited above for FLL and AMF, we must
define the screened Coulomb interaction between electrons that are treated in
DFT+U, with the help of the variable [[upawu]] and the screened exchange
interaction, with [[jpawu]]. Note that you can choose the energy unit by
indicating at the end of the line the unit abbreviation (e.g. eV or Ha). For
NiO, we will use variables that are generally accepted for this type of compound:

    upawu  8.0 0.0 eV
    jpawu  0.8 0.0 eV


You can take a look at the result of the calculation. The magnetic moment is now:


     Integrated electronic and magnetization densities in atomic spheres:
     ---------------------------------------------------------------------
     Radius=ratsph(iatom), smearing ratsm=  0.0000. Diff(up-dn)=approximate z local magnetic moment.
     Atom    Radius    up_density   dn_density  Total(up+dn)  Diff(up-dn)
        1   1.81432     8.749919     6.987384     15.737302     1.762535
        2   1.81432     6.987384     8.749919     15.737302    -1.762535
        3   1.41465     2.290397     2.290397      4.580793    -0.000000
        4   1.41465     2.290397     2.290397      4.580793    -0.000000


NiO is found antiferromagnetic, with a moment that is in reasonable agreement
with experimental results. Moreover, the system is a large gap insulator with
about 5.3 eV band gap (the 24th band at k point 3 has an eigenenergy of
0.26699 Ha, much lower than the eigenenergy of the 25th band at k point 1,
namely 0.46243 Ha, see the *tdftu_2.o_EIG* file). This number is very approximative, since the very rough
sampling of k points is not really appropriate to evaluate a band gap, still
one obtains the right physics.

A word of caution is in order here. It is NOT the case that one obtain
systematically a good result with the DFT+U method at the first trial. Indeed,
due to the nature of the modification of the energy functional, the landscape
of this energy functional might present numerous local minima (see for examples
[[cite:Jomard2008]] or [[cite:Dorado2009]]).

Unlike DFT+U, for the simple LDA (without U), in the non-spin-polarized case,
there is usually only one minimum, that is the global minimum. So, if it
converges, the self-consistency algorithm always find the same solution,
namely, the global minimum. This is already not true in the case of spin-
polarized calculations (where there might be several stable solutions of the
SCF cycles, like ferromagnetic and ferromagnetic), but usually, there are not
many local minima, and the use of the [[spinat]] input variables allows one to
adequately select the global physical characteristics of the sought solution.

By contrast, with the U, the [[spinat]] input variable is too primitive, and
one needs to be able to initialize a spin-density matrix on each atomic site
where a U is present, in order to guide the SCF algorithm.

The fact that [[spinat]] works for NiO comes from the relative simplicity of this system.

## 3 Initialization of the density matrix

*You should begin by running the tdftu_3.abi file before continuing.*

In order to help the DFT+U find the ground state, you can define the initial
density matrix for correlated orbitals with [[dmatpawu]]. For $d$ orbitals, this variable
must contains $5\times5$ square matrices. There should be one square matrix per nsppol and atom.
So in our case, there are 2 square matrices.
Also, to enable this
feature, [[usedmatpu]] must be set to a non-zero value (default is 0). When
positive, the density matrix is kept to the [[dmatpawu]] value for the
[[usedmatpu]] value steps. For our calculation(tdftu_3.abi) , [[usedmatpu]] is 5,
thus the spin-density matrix is kept constant for 5 SCF steps.
Let's examinates the input dmatpawu

{% dialog tests/tutorial/Input/tdftu_3.abi %}

To understand the density matrix used in the variable [[dmatpawu]] in this input file, have a look
to the section on this variable [[dmatpawu]]. This section show the order to orbitals in the density matrix. With
the help of this section, one can understand that the density matrix corresponds to all orbitals filled except
$e_g$ orbitals for one spin.

In the log file (not the usual output file), you will find for each step, the
calculated density matrix, followed by the imposed density matrix. After the
first 5 SCF steps, the initial density matrix is no longer imposed. Here is a
section of the log file, in which the imposed occupation matrices are echoed:

    -------------------------------------------------------------------------

    Occupation matrix for correlated orbitals is kept constant
    and equal to dmatpawu from input file !
    ----------------------------------------------------------

    == Atom   1 == Imposed occupation matrix for spin 1 ==
         0.90036    0.00000   -0.00003    0.00000    0.00000
         0.00000    0.90036   -0.00001    0.00000    0.00002
        -0.00003   -0.00001    0.91309   -0.00001    0.00000
         0.00000    0.00000   -0.00001    0.90036   -0.00002
         0.00000    0.00002    0.00000   -0.00002    0.91309

    == Atom   1 == Imposed occupation matrix for spin 2 ==
         0.89677   -0.00001    0.00011   -0.00001    0.00000
        -0.00001    0.89677    0.00006    0.00001   -0.00010
         0.00011    0.00006    0.11580    0.00006    0.00000
        -0.00001    0.00001    0.00006    0.89677    0.00010
         0.00000   -0.00010    0.00000    0.00010    0.11580

    == Atom   2 == Imposed occupation matrix for spin 1 ==
         0.89677   -0.00001    0.00011   -0.00001    0.00000
        -0.00001    0.89677    0.00006    0.00001   -0.00010
         0.00011    0.00006    0.11580    0.00006    0.00000
        -0.00001    0.00001    0.00006    0.89677    0.00010
         0.00000   -0.00010    0.00000    0.00010    0.11580

    == Atom   2 == Imposed occupation matrix for spin 2 ==
         0.90036    0.00000   -0.00003    0.00000    0.00000
         0.00000    0.90036   -0.00001    0.00000    0.00002
        -0.00003   -0.00001    0.91309   -0.00001    0.00000
         0.00000    0.00000   -0.00001    0.90036   -0.00002
         0.00000    0.00002    0.00000   -0.00002    0.91309

Generally, the DFT+U functional meets the problem of multiple local minima,
much more than the usual LDA or GGA functionals. One often gets trapped in a
local minimum. Trying different starting points might be important...

## 4 AMF double-counting method

Now we will use the other implementation for the double-counting term in DFT+U
(in Abinit), known as AMF. As the FLL method, this method uses the number of
electrons for each spin independently and the complete interactions $U(m_1,m_2,m_3,m_4)$ and $J(m_1,m_2,m_3,m_4)$.

As in the preceding run, we will start with a fixed density matrix for d
orbitals. You might now start your calculation, with the *tdftu_4.abi*, or skip the calculation, and rely on the reference file
provided in the *\$ABI_TESTS/tutorial/Refs* directory. Examine the *tdftu_4.abi* file.

{% dialog tests/tutorial/Input/tdftu_4.abi %}

The only difference in the input file compared to *tdftu_3.abi* is the
value of [[usepawu]] = 2. We obtain a band gap of 4.75 eV. The value of the
band gap with AMF and FLL is different. However, we have to remember that
these results are not well converged. By contrast, the magnetization,

      Integrated electronic and magnetization densities in atomic spheres:
      ---------------------------------------------------------------------
      Radius=ratsph(iatom), smearing ratsm=  0.0000. Diff(up-dn)=approximate z local magnetic moment.
      Atom    Radius    up_density   dn_density  Total(up+dn)  Diff(up-dn)
         1   1.81432     8.675718     6.993823     15.669541     1.681895
         2   1.81432     6.993823     8.675718     15.669541    -1.681895
         3   1.41465     2.288681     2.288681      4.577361    -0.000000
         4   1.41465     2.288681     2.288681      4.577361     0.000000


is very similar to the DFT+U FLL. 
For other systems, the difference can be more important. FLL is designed
to work well for systems in which occupations of orbitals are 0 or 1 for each
spin. The AMF should be used when orbital occupations are near the average occupancies.

## 5 Projected density of states in DFT+U

Using [[prtdos]] 3, you can now compute the projected d and f density of states.
For more information about projected density of states, for more details see the [PAW1](/tutorial/paw1) tutorial.
