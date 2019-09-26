---
authors: SPesant, MCote, XG, BAmadon
---

# Tutorial on DFT+U  

## The projected density of states of NiO.  

This tutorial aims at showing how to perform a DFT+U calculation using Abinit (see also [[cite:Amadon2008a]])

You will learn what is a DFT+U calculation and what are the main input
variables controlling this type of calculation.

It is supposed that you already know how to do PAW calculations using ABINIT.
Please follow the two tutorials on PAW in ABINIT ([PAW1](paw1), [PAW2](paw2)), if this is not the case.

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
approach is known as LDA+U (actually DFT+U). In the actual implementation, we
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

Copy the files *tdftu_1.in* and *tdftu_x.files* from *\$ABI_TESTS/tutorial/Input* to your *Work_dftu* directory with:

```sh
cd $ABI_TESTS/tutorial/Input
mkdir Work_dftu
cd Work_dftu
cp ../tdftu_x.files .  # You will need to edit this file.
cp ../tdftu_1.in .
```

{% dialog tests/tutorial/Input/tdftu_x.files tests/tutorial/Input/tdftu_1.in %}

Now run the code as usual.
The job should take less than 30 seconds on a PC 3 GHz. It calculates the LDA
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
integrated total density in the PAW spheres (see the [PAW1](paw1)
and [PAW2](paw2) tutorials on PAW formalism). This value roughly
estimate the magnetic moment of NiO:
    
     Integrated total density in atomic spheres:
     -------------------------------------------
     Atom  Sphere radius  Integrated_up_density  Integrated_dn_density  Total(up+dn)   Diff(up-dn)
        1        2.30000             9.05536980             7.85243738   16.90780718    1.20293241
        2        2.30000             7.85243738             9.05536980   16.90780718   -1.20293241
        3        1.21105             1.82716080             1.82716080    3.65432159   -0.00000000
        4        1.21105             1.82716080             1.82716080    3.65432159    0.00000000
     Note: Diff(up-dn) can be considered as a rough approximation of a local magnetic moment.

The atoms in the output file, are listed as in the [[typat]] variable (the
first two are nickel atoms and the last two are oxygen atoms). The results
indicate that spins are located in each nickel atom of the doubled primitive
cell. Fortunately, the LDA succeeds to give an antiferromagnetic ground state
for the NiO. But the result does not agree with the experimental data. 

The magnetic moment (the difference between up and down spin on the nickel atom)
range around 1.6-1.9 according to experiments  ([[cite:Cheetham1983]],[[cite:Neubeck1999]],[[cite:Sawatzky1984]],
[[cite:Hufner1984]])
Also, as the Fermi level is at 0.22347 Ha, one
can see that the band gap obtained between the last occupied (0.20672 Ha, at k
point 2) and the first unoccupied band (0.23642 Ha, at kpoint 3) is
approximately 0.8 eV which is lower than the measured value of 4.0-4.3 eV
(This value could be modified using well-converged parameters but would still
be much lower than what is expected).

Making abstraction of the effect of insufficiently convergence parameters, the
reason for the discrepancy between the DFT-LDA data and the experiments is
first the fact the DFT is a theory for the ground state and second, the lack
of correlation of the LDA. Alone, the homogeneous electron gas cannot
correctly represent the interactions among d electrons of the Ni atom. That is
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

    It is important to notice that in order to use LDA+U in Abinit, you must
    employ PAW pseudopotentials.

You should run abinit with the *tdftu_2.in* input file. This calculation takes
less than 30 seconds on a PC 3.0 GHz
During the calculation, you can take a look at the input file. 

{% dialog tests/tutorial/Input/tdftu_2.in %}

Some variable describing the LDA+U parameters have been added to the previous file. All
other parameters were kept constant from the preceding calculation. First, you
must set the variable [[usepawu]] to one (for the FLL method) and two (for the
AMT method) in order to enable the LDA+U calculation. Then, with [[lpawu]] you
give for each atomic species ([[znucl]]) the values of angular momentum (l) for
which the LDA+U correction will be applied. The choices are 2 for d-orbitals
and 3 for *f*-orbitals. You cannot treat s and p orbitals with LDA+U in the
present version of ABINIT. Also, if you do not want to apply LDA+U correction
on a species, you can set the variable to -1. For the case of NiO, we put
[[lpawu]] to 2 for Ni and -1 for O.

Finally, as described in the article cited above for FLL and AMF, we must
define the screened Coulomb interaction between electrons that are treated in
LDA+U, with the help of the variable [[upawu]] and the screened exchange
interaction, with [[jpawu]]. Note that you can choose the energy unit by
indicating at the end of the line the unit abbreviation (e.g. eV or Ha). For
NiO, we will use variables that are generally accepted for this type of compound:

    upawu 8.0 eV
    jpawu 0.8 eV 

You can take a look at the result of the calculation. The magnetic moment is now:
    
     Integrated total density in atomic spheres:
     -------------------------------------------
     Atom  Sphere radius  Integrated_up_density  Integrated_dn_density  Total(up+dn)   Diff(up-dn)
        1        2.30000             9.28514439             7.53721910   16.82236349    1.74792528
        2        2.30000             7.53721910             9.28514439   16.82236349   -1.74792528
        3        1.21105             1.84896670             1.84896670    3.69793339    0.00000000
        4        1.21105             1.84896670             1.84896670    3.69793339    0.00000000
     Note: Diff(up-dn) can be considered as a rough approximation of a local magnetic moment.

NiO is found antiferromagnetic, with a moment that is in reasonable agreement
with experimental results. Moreover, the system is a large gap insulator with
about 5.0 eV band gap (the 24th band at k point 3 has an eigenenergy of
0.15896 Ha, much lower than the eigenenergy of the 25th band at k point 1,
namely 0.24296 Ha). This number is very approximative, since the very rough
sampling of k points is not really appropriate to evaluate a band gap, still
one obtains the right physics.

A word of caution is in order here. It is NOT the case that one obtain
systematically a good result with the LDA+U method at the first trial. Indeed,
due to the nature of the modification of the energy functional, the landscape
of this energy functional might present numerous local minima.

Unlike LDA+U, for the simple LDA (without U), in the non-spin-polarized case,
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
  
*You should begin by running the tdftu_3.in file before continuing.*

In order to help the LDA+U find the ground state, you can define the initial
density matrix for correlated orbitals with [[dmatpawu]] To enable this
feature, [[usedmatpu]] must be set to a non-zero value (default is 0). When
positive, the density matrix is kept to the [[dmatpawu]] value for the
[[usedmatpu]] value steps. For our calculation(tdftu_3.in) , [[usedmatpu]] is 5, 
thus the spin-density matrix is kept constant for 5 SCF steps.

{% dialog tests/tutorial/Input/tdftu_3.in %}

In the log file (not the usual output file), you will find for each step, the
calculated density matrix, followed by the imposed density matrix. After the
first 5 SCF steps, the initial density matrix is no longer imposed. Here is a
section of the log file, in which the imposed occupation matrices are echoed:
    
    -------------------------------------------------------------------------
    
    Occupation matrix for correlated orbitals is kept constant
    and equal to initial one !
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
    
Generally, the LDA+U functional meets the problem of multiple local minima,
much more than the usual LDA or GGA functionals. One often gets trapped in a
local minimum. Trying different starting points might be important...

## 4 AMF double-counting method
  
Now we will use the other implementation for the double-counting term in LDA+U
(in Abinit), known as AMF. As the FLL method, this method uses the number of
electrons for each spin independently and the complete interactions $U(m_1,m_2,m_3,m_4)$ and $J(m_1,m_2,m_3,m_4)$.

As in the preceding run, we will start with a fixed density matrix for d
orbitals. You might now start your calculation, with the *tdftu_4.in* and
*tdftu_4.files*, or skip the calculation, and rely on the reference file
provided in the *\$ABI_TESTS/tutorial/Refs* directory. Examine the *tdftu_4.in* file. 

{% dialog tests/tutorial/Input/tdftu_4.in %}

The only difference in the input file compared to *tdftu_3.in* is the
value of [[usepawu]] = 2. We obtain a band gap of 4.3 eV. The value of the
band gap with AMF and FLL is different. However, we have to remember that
these results are not well converged. By contrast, the magnetization,
    
     Atom  Sphere radius  Integrated_up_density  Integrated_dn_density  Total(up+dn)   Diff(up-dn)
        1        2.30000             9.24026835             7.56013140   16.80039975    1.68013694
        2        2.30000             7.56013140             9.24026835   16.80039975   -1.68013694
        3        1.21105             1.84683993             1.84683993    3.69367986   -0.00000000
        4        1.21105             1.84683993             1.84683993    3.69367986    0.00000000
     Note: Diff(up-dn) can be considered as a rough approximation of a local magnetic moment.

is very similar to the LDA+U FLL. In fact, this system is not very complex.
But for other systems, the difference can be more important. FLL is designed
to work well for crystal with diagonal occupation matrix with 0 or 1 for each
spin. The AMF should be used when orbital occupations are near the average occupancies.

## 5 Projected density of states in LDA+U
  
Using [[prtdos]] 3, you can now compute the projected d and f density of states.
For more information about projected density of states, for more details see the [PAW1](paw1) tutorial.
