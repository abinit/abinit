---
authors: XG, RC
---

# Second (basic) tutorial

## The H<sub>2</sub> molecule, with convergence studies.

This tutorial aims at showing how to get converged values for the following physical properties:

* the bond length
* the atomisation energy

You will learn about the numerical quality of the calculations, then make
convergence studies with respect to the number of planewaves and the size of
the supercell, and finally consider the effect of the XC functional. The
problems related to the use of different pseudopotential are not examined.
You will also finish to read the [[help:abinit]].

This tutorial should take about 1 hour.

[TUTORIAL_READMEV9]

## Summary of the previous tutorial

We studied the H$_2$ molecule in a big box.
We used 10 Ha as cut-off energy, a 10x10x10 Bohr$^3$ supercell, the local-density approximation
(as well as the local-spin-density approximation) in the Teter parametrization ([[ixc]] = 1, the
default), and a pseudopotential from the Goedecker-Hutter-Teter table.

At this stage, we compared our results:

* bond length: 1.522 Bohr
* atomisation energy at that bond length: 0.1656 Ha = 4.506 eV

with the experimental data (as well as theoretical data using a much more accurate technique than DFT)

  * bond length: 1.401 Bohr
  * atomisation energy: 4.747 eV

The bond length is awful (nearly 10% off), and the atomisation energy is a bit too low, 5% off.

## 2 The convergence in ecut (I)

**2.1.a** **Computing the bond length and corresponding atomisation energy in one run.**

*Before beginning, you might consider to work in a different subdirectory as for tutorial 1.
Why not Work2?*

Because we will compute many times the bond length and atomisation energy, it
is worth to make a single input file that will do all the associated operations.
You should try to use 2 datasets (try to combine *\$ABI_TESTS/tutorial/Input/tbase1_3.in* with *tbase1_5.in*).
Do not try to have the same position of the H atom as one of the H$_2$ atoms in the optimized geometry.

```sh
cd $ABI_TESTS/tutorial/Input
mkdir Work2
cd Work2
cp ../tbase2_x.files .   # You will need to edit this file.
cp ../tbase2_1.in .
```

The input file *tbase2_1.in* is an example of file that will do the job,

{% dialog tests/tutorial/Input/tbase2_1.in %}

while *tbase2_1.out* is an example of output file:

{% dialog tests/tutorial/Refs/tbase2_1.out %}

You might use *$ABI_TESTS/tutorial/Input/tbase2_x.files* as *files* file
(do not forget to modify it, like in [[lesson:base1|tutorial 1]],
although it does not differ from *tbase1_x.files*.

{% dialog tests/tutorial/Input/tbase2_x.files %}

Execute the code with:

    abinit < tbase2_x.files > log 2> err &

The run should take less than one minute.

You should obtain the values:

        etotal1  -1.1058360644E+00
        etotal2  -4.7010531489E-01

and

        xcart1  -7.6091015760E-01  0.0000000000E+00  0.0000000000E+00
                 7.6091015760E-01  0.0000000000E+00  0.0000000000E+00

These are similar to those determined in [tutorial 1](base1),
although they have been obtained in one run.
You can also check that the residual forces are lower than `5.0d-4`.
Convergence issues are discussed in [[help:abinit#numerical-quality|section 6]] of the abinit help file, on numerical quality.
You should read it.
By the way, you have read many parts of the abinit help file!
You are missing the sections [[help:abinit#2|2]], [[help:abinit#pseudopotential-files|5]],  [[help:abinit#numerical-quality|6]].

You are also missing the description of many input variables.
We suggest that you finish reading entirely the abinit help file now, while
the knowledge of the input variables will come in the long run.

**2.1.b** Many convergence parameters have already been identified. We will
focus only on [[ecut]] and [[acell]]. This is because

* the convergence of the SCF cycle and geometry determination are well
   under control thanks to [[toldfe]], [[toldff]] and [[tolmxf]]
   (this might not be the case for other physical properties)

* there is no k point convergence study to be done for an isolated system in a big box:
  no additional information is gained by adding a k-point beyond one

* the boxcut value is automatically chosen larger than 2 by ABINIT, see the determination of the
  input variable [[ngfft]] by preprocessing

* we are using [[ionmov]] = 2 for the determination of the geometry.

## 3 The convergence in ecut (II)

For the check of convergence with respect to [[ecut]], you have the choice
between doing different runs of the *tbase2_1.in* file with different values of
[[ecut]], or doing a double loop of datasets, as proposed in *$ABI_TESTS/tutorial/Input/tbase2_2.in*.
The values of [[ecut]] have been chosen between 10 Ha and 35 Ha, by step of 5 Ha.
If you want to make a double loop, you might benefit of reading again the
[[help:abinit#loop|double-loop section]] of the abinit_help file.

**2.2.a** You have likely seen a big increase of the CPU time needed to do the
calculation. You should also look at the increase of the memory needed to do
the calculation (go back to the beginning of the output file).
The output data are as follows:

        etotal11 -1.1058360644E+00
        etotal12 -4.7010531489E-01
        etotal21 -1.1218716100E+00
        etotal22 -4.7529731401E-01
        etotal31 -1.1291943792E+00
        etotal32 -4.7773586424E-01
        etotal41 -1.1326879404E+00
        etotal42 -4.7899908214E-01
        etotal51 -1.1346739190E+00
        etotal52 -4.7972721653E-01
        etotal61 -1.1359660026E+00
        etotal62 -4.8022016459E-01

         xcart11 -7.6091015760E-01  0.0000000000E+00  0.0000000000E+00
                  7.6091015760E-01  0.0000000000E+00  0.0000000000E+00
         xcart12  0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
         xcart21 -7.5104912643E-01  0.0000000000E+00  0.0000000000E+00
                  7.5104912643E-01  0.0000000000E+00  0.0000000000E+00
         xcart22  0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
         xcart31 -7.3977108831E-01  0.0000000000E+00  0.0000000000E+00
                  7.3977108831E-01  0.0000000000E+00  0.0000000000E+00
         xcart32  0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
         xcart41 -7.3304273322E-01  0.0000000000E+00  0.0000000000E+00
                  7.3304273322E-01  0.0000000000E+00  0.0000000000E+00
         xcart42  0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
         xcart51 -7.3001570260E-01  0.0000000000E+00  0.0000000000E+00
                  7.3001570260E-01  0.0000000000E+00  0.0000000000E+00
         xcart52  0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
         xcart61 -7.2955902118E-01  0.0000000000E+00  0.0000000000E+00
                  7.2955902118E-01  0.0000000000E+00  0.0000000000E+00
         xcart62  0.0000000000E+00  0.0000000000E+00  0.0000000000E+00

The corresponding atomisation energies and interatomic distances are:

| ecut (Ha)  |  atomisation energy (Ha)  | interatomic distance (Bohr)
| :--        | :--                       | :--
  10         | .1656                     |  1.522
  15         | .1713                     |  1.502
  20         | .1737                     |  1.480
  25         | .1747                     |  1.466
  30         | .1753                     |  1.460
  35         | .1756                     |  1.459

In order to obtain 0.2% relative accuracy on the bond length or atomisation
energy, one should use a kinetic cut-off energy of 30 Ha.
We will keep in mind this value for the final run.

## 4 The convergence in acell

The same technique as for [[ecut]] should be now used for the convergence in [[acell]].
We will explore [[acell]] starting from `8 8 8` to `18 18 18`, by step of `2 2 2`.
We keep [[ecut]] 10 for this study. Indeed, it is a rather general rule that there is
little cross-influence between the convergence of [[ecut]] and the convergence of [[acell]].

The file *$ABI_TESTS/tutorial/Input/tbase2_3.in* can be used as an example.

{% dialog tests/tutorial/Input/tbase2_3.in %}

The output results in *$ABI_TESTS/tutorial/Refs/tbase2_3.out* are as follows:

        etotal11   -1.1188124709E+00
        etotal12   -4.8074164402E-01
        etotal21   -1.1058360838E+00
        etotal22   -4.7010531489E-01
        etotal31   -1.1039109527E+00
        etotal32   -4.6767804802E-01
        etotal41   -1.1039012868E+00
        etotal42   -4.6743724199E-01
        etotal51   -1.1041439411E+00
        etotal52   -4.6735895176E-01
        etotal61   -1.1042058281E+00
        etotal62   -4.6736729718E-01

         xcart11   -7.8330751426E-01  0.0000000000E+00  0.0000000000E+00
                    7.8330751426E-01  0.0000000000E+00  0.0000000000E+00
         xcart12    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
         xcart21   -7.6024281092E-01  0.0000000000E+00  0.0000000000E+00
                    7.6024281092E-01  0.0000000000E+00  0.0000000000E+00
         xcart22    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
         xcart31   -7.5428234893E-01  0.0000000000E+00  0.0000000000E+00
                    7.5428234893E-01  0.0000000000E+00  0.0000000000E+00
         xcart32    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
         xcart41   -7.5446921004E-01  0.0000000000E+00  0.0000000000E+00
                    7.5446921004E-01  0.0000000000E+00  0.0000000000E+00
         xcart42    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
         xcart51   -7.5384974520E-01  0.0000000000E+00  0.0000000000E+00
                    7.5384974520E-01  0.0000000000E+00  0.0000000000E+00
         xcart52    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
         xcart61   -7.5373336127E-01  0.0000000000E+00  0.0000000000E+00
                    7.5373336127E-01  0.0000000000E+00  0.0000000000E+00
         xcart62    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00

The corresponding atomisation energies and interatomic distances are:

| acell (Bohr) | atomisation energy (Ha) | interatomic distance (Bohr)
| :--  | :-- | :--
8  | .1574 | 1.568
10 | .1656 | 1.522
12 | .1686 | 1.509
14 | .1691 | 1.510
16 | .1694 | 1.508
18 | .1695 | 1.508

In order to reach 0.2% convergence on the interatomic distance, one needs `acell 12 12 12`.
The atomisation energy needs `acell 14 14 14` to be converged at that level.
At `12 12 12`, the difference is .0009 Ha = 0.024 eV, which is sufficiently small for practical purposes.
We will use `acell 12 12 12` for the final run.

For most solids the size of the unit cell will be smaller than that.
We are treating a lot of vacuum in this supercell!
So, the H$_2$ study, with this pseudopotential, turns out to be not really easy.
Of course, the number of states to be treated is minimal!
This allows to have reasonable CPU time still.

## 5 The final calculation in Local (Spin) Density Approximation

We now use the correct values of both [[ecut]] and [[acell]].
Well, you should modify the *tbase2_3.in* file to make a calculation with `acell 12 12 12` and `ecut 30`.
You can still use the double loop feature with `udtset 1 2`
(which reduces to a single loop), to minimize the modifications to the file.

The file *$ABI_TESTS/tutorial/Input/tbase2_4.in* can be taken as an example of input file:

{% dialog tests/tutorial/Input/tbase2_4.in %}

while *$ABI_TESTS/tutorial/Refs/tbase2_4.out* is as an example of output file:

{% dialog tests/tutorial/Refs/tbase2_4.out %}

Since we are doing the calculation at a single ([[ecut]], [[acell]]) pair, the
total CPU time is not as much as for the previous determinations of optimal
values through series calculations.
However, the memory needs have still increased a bit.

The output data are:

        etotal11 -1.1329369190E+00
        etotal12 -4.7765320721E-01

         xcart11 -7.2594741339E-01  0.0000000000E+00  0.0000000000E+00
                  7.2594741339E-01  0.0000000000E+00  0.0000000000E+00
         xcart12  0.0000000000E+00  0.0000000000E+00  0.0000000000E+00

* The corresponding atomisation energy is `0.1776 Ha = 4.833 eV`
* The interatomic distance is 1.452 Bohr.
* These are our final data for the local (spin) density approximation.

We have used [[ixc]] = 1.
Other expressions for the local (spin) density approximation [2, 3 ... 7] are possible.
The values 1, 2, 3 and 7 should give about the same results, since they all start
from the XC energy of the homogeneous electron gas, as determined by Quantum Monte Carlo calculations.
Other possibilities (*ixc* = 4, 5, 6) are older local density functionals, that could not rely on these data.

## 6 The use of the Generalized Gradient Approximation

We will use the Perdew-Burke-Ernzerhof functional proposed in [[cite:Perdew1996]]

In principle, for GGA, one should use another pseudopotential than for LDA.
However, for the special case of Hydrogen, and in general pseudopotentials
with a very small core (including only the 1s orbital), pseudopotentials
issued from the LDA and from the GGA are very similar.
So, we will not change our pseudopotential.
This will save us lot of time, as we should not redo an [[ecut]] convergence test

!!! important

    *ecut* is often characteristic of the pseudopotentials that are used in a calculation.

Independently of the pseudopotential, an [[acell]] convergence test should not
be done again, since the vacuum is treated similarly in LDA or GGA.

So, our final values within GGA will be easily obtained by setting [[ixc]] to 11 in *tbase2_4.in*.

{% dialog tests/tutorial/Input/tbase2_5.in %}

        etotal11 -1.1621428376E+00
        etotal12 -4.9869631917E-01

         xcart11 -7.1190611804E-01  0.0000000000E+00  0.0000000000E+00
                  7.1190611804E-01  0.0000000000E+00  0.0000000000E+00
         xcart12  0.0000000000E+00  0.0000000000E+00  0.0000000000E+00

* The corresponding atomisation energy is 0.1648 Ha = 4.483 eV
* The interatomic distance is 1.424 Bohr.
* These are our final data for the generalized gradient approximation.

Once more, here are the experimental data:

* bond length: 1.401 Bohr
* atomisation energy: 4.747 eV

In GGA, we are within 2% of the experimental bond length, but 5% of the experimental atomisation energy.
In LDA, we were within 4% of the experimental bond length, and within 2% of the experimental atomisation energy.

!!! important

    Do not forget that the typical accuracy of LDA and GGA varies with the class of materials studied.
