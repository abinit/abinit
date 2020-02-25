---
authors: DJA
---

# Calculation of *U* and J using Cococcioni's approach

## How to determine *U* for DFT+*U* in ABINIT? Cococcioni's approach.

This tutorial aims to show how you can determine *U* for further DFT+*U*
calculations consistently and in a fast an easy way. You will learn how to prepare
the input files for the determination and how to use the main parameters implemented for this aim.
It is supposed that you already know how to run ABINIT with PAW  (tutorial [PAW1](paw1)).
Obviously, you should also read the tutorial [DFT+U](dftu), and likely the tutorial [PAW2](paw2),
to generate PAW atomic data.

This tutorial should take about 1/2 hour.

[TUTORIAL_README]

## Summary of linear response method to determine *U*

The linear response method has been introduced by several authors
[[cite:Cococcioni2002]], [[cite:Cococcioni2005]],[[cite:Dederichs1984]],[[cite:Hybertsen1989]],
[[cite:Anisimov1991]],[[cite:Pickett1998]].
It is based on the fact that *U* corresponds to the energy to localize an additional
electron on the same site: *U=E[n+1]+E[n-1]-2E[n]* [[cite:Hybertsen1989]]. This can be reformulated
as the response to an infinitesimal change of of occupation of the orbital by
the electrons dn. Then *U* is the second derivative of the energy with respect
to the occupation $U=\frac{\delta^2 E}{\delta^2 n}$. The first method fixed the occupation by
cutting the hopping terms of localized orbitals. Later propositions
constrained the occupation through Lagrange multipliers [[cite:Dederichs1984]],[[cite:Anisimov1991]]. The Lagrange
multiplier $\alpha$ corresponds to a local potential that has to be applied to
increase or decrease the occupation by +/-1 electron. Note that the occupation
need not to vary by 1 electron, but the occupation shift can be infinitesimal.

It is recommended to read the following papers to understand the basic
concepts of the linear response calculations to calculate *U*:

[1] "A LDA+*U* study of selected iron compounds ", M. Cococcioni, Ph.D. thesis,
International School for Advanced Studies (SISSA), Trieste (2002)  [[cite:Cococcioni2002]]

[2] "Linear response approach to the calculation of the effective interaction
parameters in the LDA + *U* method", M. Cococcioni and S. de Gironcoli, Physical
Review B 71, 035105 (2005)  [[cite:Cococcioni2005]]

Some further reading:

[3] "Ground States of Constrained Systems: Application to Cerium Impurities",
P. H. Dederichs, S. Blugel, R. Zeller, and H. Akai, Phys. Rev. Lett. 53, 2512 (1984)   [[cite:Dederichs1984]]

[4] "Calculation of Coulomb-interaction parameters for La2CuO4 using a
constrained-density-functional approach", M. S. Hybertsen, M. Schluter, and N.
E. Christensen, Phys. Rev. B 39, 9028 (1989)   [[cite:Hybertsen1989]]

[5] "Density-functional calculation of effective Coulomb interactions in
metals", V. I. Anisimov and O. Gunnarsson, Phys. Rev. B42, 7570 (1991)  [[cite:Anisimov1991]]

[6] "Reformulation of the LDA+*U* method for a local-orbital basis", W. E.
Pickett, S. C. Erwin, and E. C. Ethridge, Phys. Rev. B58, 1201 (1998)  [[cite:Pickett1998]]

The implementation of the determination of *U* in ABINIT is briefly discussed in  [[cite:Gonze2016]].

## How to determine *U* in ABINIT

*Before continuing, you may consider to work in a different subdirectory as
for the other tutorials. Why not Work_udet?*

!!! important

    In what follows, the name of files are mentioned as if you were in this subdirectory.
    All the input files can be found in the $ABI_TESTS/tutorial/Input directory
     You can compare your results with reference output files located in
     $ABI_TESTS/tutorial/Refs directory (for the present tutorial they are named tudet*.out).

The input file *tudet_1.in* is an example of a file to prepare a wave function
for further processing. You might use the file *tudet_1.files* as a "files"
file, and get the corresponding output file ../Refs/tudet_1.out).

Copy the files *tudet_1.in* and *tudet_1.files* in your work directory, and run ABINIT:

```sh
cd $ABI_TESTS/tutorial/Input
mkdir Work_udet
cd Work_udet
cp ../tudet_1.files .
cp ../tudet_1.in .

abinit < tudet_1.files > log 2> err &
```

In the meantime, you can read the input file and see that this is a usual
DFT+U calculation, with *U*=0.

{% dialog tests/tutorial/Input/tudet_1.files tests/tutorial/Input/tudet_1.in %}

This setting allows us to read the occupations of
the Fe 3d orbitals ([[lpawu]] 2). The cell contains 2 atoms. This is the
minimum to get reasonable response matrices. We converge the electronic
structure to a high accuracy ([[tolvrs]] 10d-12), which usually allows one to
determine occupations with a precision of 10d-10. The [[ecut]] is chosen very low,
in order to accelerate calculations.
We do not suppress the writing of the *WFK* file, because this is the input for
the calculations of U.

Once this calculation has finished, run the second one:
Copy the files *tudet_2.in and tudet_2.files* in your work directory, and run ABINIT:

    abinit < tudet_2.files > tudet_1.log

{% dialog tests/tutorial/Input/tudet_2.files tests/tutorial/Input/tudet_2.in %}

As you can see from the *tudet_2.files* file, this run uses the *tudet_1o_WFK* as
an input. In the *tudet_2.in* all the symmetry relations are specified
explicitly. In the *tudet_2.log* you can verify that none of the symmetries
connects atoms 1 with atom 2:

    symatm: atom number    1 is reached starting at atom

       1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1

       1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1

     symatm: atom number    2 is reached starting at atom

       2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2

       2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2

This is important. Otherwise the occupation numbers have no freedom to evolve
separately on the atoms surrounding the atom on which you apply the perturbation.

You can generate these symmetries, in a separate run, where you specify the
atom where the perturbation is done as a different species. From the output
you read the number of symmetries ([[nsym]]), the symmetry operations
([[symrel]]) and the non-symmorphic vectors ([[tnons]]). This is already
done here and inserted in the *tudet_2.in* file. Note that you can alternatively
disable all symmetries with [[nsym]] = 1, or break specific symmetries by
displacing the impurity atom in the preliminary run. However, for the
determination of *U*, the positions should be the ideal positions and only the
symmetry should be reduced.

For the rest, usually it is enough to set [[macro_uj]] 1 to run the
calculation of *U*. Note also, that the [[irdwfk]] 1 and the [[tolvrs]] 1d-8
need not be set explicitly, because they are the defaults with [[macro_uj]] 1.

Once the calculation tudet_2 is converged, you can have look at the output.
You can see, that the atomic shift ([[atvshift]]) is automatically set:

             atvshift      0.00367    0.00367    0.00367    0.00367    0.00367
                           0.00367    0.00367    0.00367    0.00367    0.00367
                           0.00000    0.00000    0.00000    0.00000    0.00000
                           0.00000    0.00000    0.00000    0.00000    0.00000

This means, that all the 10 3d spin-spin orbitals on the first Fe atom where
shifted by 0.1 eV (=0.00367 Ha). On the second atom no shift was applied.
Self-consistency was reached twice: Once for a positive shift, once for the negative shift:

    grep SCF  tudet_2.out

The lines starting with URES

     URES      ii    nat       r_max    U(J)[eV]   U_ASA[eV]   U_inf[eV]
     URES       1      2     4.69390     4.74555     3.67983     3.20150
     URES       2     16     9.38770     8.77694     6.80588     5.92122
     URES       3     54    14.08160     9.17082     7.11130     6.18694
     URES       4    128    18.77540     9.25647     7.17772     6.24472
     URES       5    250    23.46930     9.28509     7.19991     6.26403

contain *U* for different supercells. The column "nat" indicates how many atoms
were involved in the supercell, r_max indicates the maximal distance of the
impurity atoms in that supercell. The column *U* indicates the actual *U* you
calculated and should use in your further calculations. U_ASA is an estimate
of *U* for more extended projectors and U_\inf is the estimate for a projector
extended even further.

Although it is enough to set [[macro_uj]] 1, you can further tune your runs.
As a standard, the potential shift to the 1st atom treated in DFT+*U*, with a
potential shift of 0.1 eV. If you wish to determine *U* on the second atom you
put [[pawujat]] 2. To change the size of the potential shift use e.g.
[[pawujv]] 0.05 eV. Our tests show that 0.1 eV is the optimal value, but the
linear response is linear in a wide range (1-0.001 eV).

## The ujdet utility

In general the calculation of *U* with abinit as described above is sufficient.
For some post-treatment that goes beyond the standard applications, a separate
executable *ujdet* was created. The output of abinit is formatted so that you
can easily "cut" the part with the ujdet input variables: you can generate
the standard input file for the ujdet utility by typing:

    sed -n "/MARK/,/MARK/p" tudet_2.out  > ujdet.in

Note that the input for the ujdet utility is always called *ujdet.in*

It contains the potential shifts applied vsh (there are 4 shifts: vsh1, vsh3
for non-selfconsistent calculations that allows to extract the contribution to
*U* originating from a non-interacting electron gas, and vsh2, vsh4 for positive
and negative potential shift). The same applies for the occupations occ[1-4].

We now calculate *U* for an even larger supercell: Uncomment the line scdim in *ujdet.in* and add

     scdim 6 6 6

to specify a 6 6 6 supercell or

     scdim 700 0 0

to specify the maximum total number of atoms in the supercell. Then, run *ujdet*:

    rm ujdet.[ol]* ; ujdet > ujdet.log

    grep URES ujdet.out

     URES      ii    nat       r_max    U(J)[eV]   U_ASA[eV]   U_inf[eV]
     URES       1      2     4.69390     4.74555     3.67983     3.20150
     URES       2     16     9.38770     8.77694     6.80588     5.92122
     URES       3     54    14.08160     9.17082     7.11130     6.18694
     URES       4    128    18.77540     9.25647     7.17772     6.24472
     URES       5    250    23.46930     9.28509     7.19991     6.26403
     URES       6    432    28.16310     9.29738     7.20944     6.27232


As you can see, *U* has now been extrapolated to a supercell containing 432 atoms.

The value of *U* depends strongly on the extension of the projectors used in the
calculation. If you want to use *U* in LMTO-ASA calculations you can use the
keyword [[pawujrad]] in the *ujdet.in* file to get grips of the *U* you want to
use there. Just uncomment the line and add the ASA-radius of the specific atom e.g.

    pawujrad 2.5

Running

    rm ujdet.[ol]* ; ujdet > ujdet.log

gives now higher values in the column U_ASA than in the runs before (8.07 eV
compared to 7.21 eV): For more localized projectors the *U* value has to be bigger.
