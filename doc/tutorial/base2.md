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

[TUTORIAL_README]

## Summary of the previous tutorial

We studied the H$_2$ molecule in a big box.
We used 10 Ha as cut-off energy, a 10x10x10 Bohr$^3$ supercell, the local-density approximation
(as well as the local-spin-density approximation) in the Perdew-Wang parametrization ([[ixc]] = -1012)
and a pseudopotential from the pseudodojo <http://www.pseudo-dojo.org/>.

At this stage, we compared our results:

* bond length: 1.486 Bohr
* atomisation energy at that bond length: 0.1704 Ha = 4.635 eV

with the experimental data (as well as theoretical data using a much more accurate technique than DFT)

  * bond length: 1.401 Bohr
  * atomisation energy: 4.747 eV

The bond length is rather bad (about 6% off), and the atomisation energy is a bit too low, 2.5% off.

## 2 The convergence in ecut (I)

**2.1.a** **Computing the bond length and corresponding atomisation energy in one run.**

*Before beginning, you might consider to work in a different subdirectory as for tutorial 1.
Why not Work2?*

Because we will compute many times the bond length and atomisation energy, it
is worth to make a single input file that will do all the associated operations.
You should try to use 2 datasets (try to combine *\$ABI_TESTS/tutorial/Input/tbase1_3.abi* with *tbase1_5.abi*).
Do not try to have the same position of the H atom as one of the H$_2$ atoms in the optimized geometry.

```sh
cd $ABI_TESTS/tutorial/Input
mkdir Work2
cd Work2
cp ../tbase2_1.abi .
```

The input file *tbase2_1.abi* is an example of file that will do the job,

{% dialog tests/tutorial/Input/tbase2_1.abi %}

while *tbase2_1.abo* is an example of output file:

{% dialog tests/tutorial/Refs/tbase2_1.abo %}

Execute the code with:

    abinit tbase2_1.abi > log &

The run should take less than one minute.

You should obtain the values:

           etotal1    -1.1182883137E+00
           etotal2    -4.7393103688E-01

and

            xcart1    -7.4307169181E-01  0.0000000000E+00  0.0000000000E+00
                       7.4307169181E-01  0.0000000000E+00  0.0000000000E+00

These are similar to those determined in [tutorial 1](/tutorial/base1),
although they have been obtained in one run.
You can also check that the residual forces are lower than `5.0d-4`.
Convergence issues are discussed in [[help:abinit#numerical-quality|section 6]] of the abinit help file, on numerical quality.
You should read it.
By the way, you have read many parts of the abinit help file!
You are missing the sections (or part of) [[help:abinit#2|2]], [[help:abinit#pseudopotential-files|5]],  [[help:abinit#numerical-quality|6]].

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

* the boxcut value (see [[boxcutmin]]) is automatically chosen larger than 2 by ABINIT, see the determination of the
  input variable [[ngfft]] by preprocessing

* we are using [[ionmov]] = 2 for the determination of the geometry.

## 3 The convergence in ecut (II)

For the check of convergence with respect to [[ecut]], you have the choice
between doing different runs of the *tbase2_1.abi* file with different values of
[[ecut]], or doing a double loop of datasets, as proposed in *$ABI_TESTS/tutorial/Input/tbase2_2.abi*.
The values of [[ecut]] have been chosen between 10 Ha and 35 Ha, by step of 5 Ha.
If you want to make a double loop, you might benefit of reading again the
[[help:abinit#loop|double-loop section]] of the abinit_help file.

**2.2.a** You have likely seen a big increase of the CPU time needed to do the
calculation. You should also look at the increase of the memory needed to do
the calculation (go back to the beginning of the output file).
The output data are as follows:

           etotal11   -1.1182883137E+00
           etotal12   -4.7393103688E-01
           etotal21   -1.1325211211E+00
           etotal22   -4.7860857539E-01
           etotal31   -1.1374371581E+00
           etotal32   -4.8027186429E-01
           etotal41   -1.1389569555E+00
           etotal42   -4.8083335144E-01
           etotal51   -1.1394234882E+00
           etotal52   -4.8101478048E-01
           etotal61   -1.1395511942E+00
           etotal62   -4.8107063412E-01

            xcart11   -7.4307169181E-01  0.0000000000E+00  0.0000000000E+00
                       7.4307169181E-01  0.0000000000E+00  0.0000000000E+00
            xcart12    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
            xcart21   -7.3687974546E-01  0.0000000000E+00  0.0000000000E+00
                       7.3687974546E-01  0.0000000000E+00  0.0000000000E+00
            xcart22    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
            xcart31   -7.3014665027E-01  0.0000000000E+00  0.0000000000E+00
                       7.3014665027E-01  0.0000000000E+00  0.0000000000E+00
            xcart32    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
            xcart41   -7.2642579309E-01  0.0000000000E+00  0.0000000000E+00
                       7.2642579309E-01  0.0000000000E+00  0.0000000000E+00
            xcart42    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
            xcart51   -7.2563260546E-01  0.0000000000E+00  0.0000000000E+00
                       7.2563260546E-01  0.0000000000E+00  0.0000000000E+00
            xcart52    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
            xcart61   -7.2554339763E-01  0.0000000000E+00  0.0000000000E+00
                       7.2554339763E-01  0.0000000000E+00  0.0000000000E+00
            xcart62    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00

The corresponding atomisation energies and interatomic distances are:

| ecut (Ha)  |  atomisation energy (Ha)  | interatomic distance (Bohr)
| :--        | :--                       | :--
  10         | .1704                     |  1.486
  15         | .1753                     |  1.474
  20         | .1769                     |  1.460
  25         | .1773                     |  1.453
  30         | .1774                     |  1.451
  35         | .1774                     |  1.451

In order to obtain 0.2% relative accuracy on the bond length or atomisation
energy, one should use a kinetic cut-off energy of 25 Ha.
We will keep in mind this value for the final run.

## 4 The convergence in acell

The same technique as for [[ecut]] should be now used for the convergence in [[acell]].
We will explore [[acell]] starting from `8 8 8` to `18 18 18`, by step of `2 2 2`.
We keep [[ecut]] 10 for this study. Indeed, it is a rather general rule that there is
little cross-influence between the convergence of [[ecut]] and the convergence of [[acell]].

The file *$ABI_TESTS/tutorial/Input/tbase2_3.abi* can be used as an example.

{% dialog tests/tutorial/Input/tbase2_3.abi %}

The output results in *$ABI_TESTS/tutorial/Refs/tbase2_3.abo* are as follows:

           etotal11   -1.1305202335E+00
           etotal12   -4.8429570903E-01
           etotal21   -1.1182883137E+00
           etotal22   -4.7393103688E-01
           etotal31   -1.1165450484E+00
           etotal32   -4.7158917506E-01
           etotal41   -1.1165327748E+00
           etotal42   -4.7136118536E-01
           etotal51   -1.1167740301E+00
           etotal52   -4.7128698890E-01
           etotal61   -1.1168374331E+00
           etotal62   -4.7129589330E-01

            xcart11   -7.6471149217E-01  0.0000000000E+00  0.0000000000E+00
                       7.6471149217E-01  0.0000000000E+00  0.0000000000E+00
            xcart12    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
            xcart21   -7.4307169181E-01  0.0000000000E+00  0.0000000000E+00
                       7.4307169181E-01  0.0000000000E+00  0.0000000000E+00
            xcart22    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
            xcart31   -7.3778405090E-01  0.0000000000E+00  0.0000000000E+00
                       7.3778405090E-01  0.0000000000E+00  0.0000000000E+00
            xcart32    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
            xcart41   -7.3794243127E-01  0.0000000000E+00  0.0000000000E+00
                       7.3794243127E-01  0.0000000000E+00  0.0000000000E+00
            xcart42    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
            xcart51   -7.3742475720E-01  0.0000000000E+00  0.0000000000E+00
                       7.3742475720E-01  0.0000000000E+00  0.0000000000E+00
            xcart52    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
            xcart61   -7.3733248368E-01  0.0000000000E+00  0.0000000000E+00
                       7.3733248368E-01  0.0000000000E+00  0.0000000000E+00
            xcart62    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00

The corresponding atomisation energies and interatomic distances are:

| acell (Bohr) | atomisation energy (Ha) | interatomic distance (Bohr)
| :--  | :-- | :--
8  | .1619 | 1.529
10 | .1704 | 1.486
12 | .1734 | 1.476
14 | .1738 | 1.478
16 | .1742 | 1.475
18 | .1742 | 1.475

In order to reach 0.2% convergence on the atomisation energy and interatomic distance one needs `acell 16 16 16`.
We will use `acell 16 16 16` for the final run.

For most solids the size of the unit cell will be smaller than that.
We are treating a lot of vacuum in this supercell !
So, the H$_2$ study, with this pseudopotential, turns out to be not really easy.
Of course, the number of states to be treated is minimal!
This allows to have reasonable CPU time still.

## 5 The final calculation in Local (Spin) Density Approximation

We now use the correct values of both [[ecut]] and [[acell]].
Well, you should modify the *tbase2_3.abi* file to make a calculation with `acell 16 16 16` and `ecut 25`.
You can still use the double loop feature with `udtset 1 2`
(which reduces to a single loop), to minimize the modifications to the file.

The file *$ABI_TESTS/tutorial/Input/tbase2_4.abi* can be taken as an example of input file:

{% dialog tests/tutorial/Input/tbase2_4.abi %}

while *$ABI_TESTS/tutorial/Refs/tbase2_4.abo* is as an example of output file:

{% dialog tests/tutorial/Refs/tbase2_4.abo %}

Since we are doing the calculation at a single ([[ecut]], [[acell]]) pair, the
total CPU time is not as much as for the previous determinations of optimal
values through series calculations.
However, the memory needs have still increased a bit.

The output data are:

           etotal11   -1.1369766875E+00
           etotal12   -4.7827555035E-01

            xcart11   -7.2259811794E-01  0.0000000000E+00  0.0000000000E+00
                       7.2259811794E-01  0.0000000000E+00  0.0000000000E+00
            xcart12    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00

* The corresponding atomisation energy is `0.1804 Ha = 4.910 eV`
* The interatomic distance is 1.445 Bohr.
* These are our final data for the local (spin) density approximation.

Without our choice of pseudopotential, the value of [[ixc]] was -1012, corresponding
to the Perdew-Wang [[cite:Perdew1992a]] parameterization of the LDA XC functional.
It is in principle the same as using [[ixc]] = 1.
Other expressions for the local (spin) density approximation [[ixc]]=[2, 3 ... 7] are possible.
The values 1, 2, 3 and 7 should give about the same results, since they all start
from the XC energy of the homogeneous electron gas, as determined by Quantum Monte Carlo calculations.
Other possibilities ([[ixc]] = 4, 5, 6) are older local density functionals, that could not rely on these data.

## 6 The use of the Generalized Gradient Approximation

We will use the Perdew-Burke-Ernzerhof functional proposed in [[cite:Perdew1996]]

For GGA, we use another pseudopotential than for LDA.
In principle we should redo the [[ecut]] convergence test, possibly coming to the conclusion
that another value of [[ecut]] should be used.
However, for the special case of Hydrogen, and in general pseudopotentials
with a very small core (including only the 1s orbital), pseudopotentials
issued from the LDA and from the GGA are very similar.
So, we will not redo an [[ecut]] convergence test.

!!! important

    *ecut* is often characteristic of the pseudopotentials that are used in a calculation.

Independently of the pseudopotential, an [[acell]] convergence test should not
be done again, since the vacuum is treated similarly in LDA or GGA.

So, our final values within GGA will be easily obtained by changing the pseudopotential with respect to the one used in *tbase2_4.abi*.

{% dialog tests/tutorial/Input/tbase2_5.abi %}

           etotal11   -1.1658082573E+00
           etotal12   -4.9940910146E-01

            xcart11   -7.0870975055E-01 -2.3009273766E-32 -6.6702161522E-32
                       7.0870975055E-01  2.3009273766E-32  6.6702161522E-32
            xcart12    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00

* The corresponding atomisation energy is 0.1670 Ha = 4.544 eV
* The interatomic distance is 1.417 Bohr.
* These are our final data for the generalized gradient approximation.

Once more, here are the experimental data:

* bond length: 1.401 Bohr
* atomisation energy: 4.747 eV

In GGA, we are within 2% of the experimental bond length, but 5% of the experimental atomisation energy.
In LDA, we were within 3% of the experimental bond length, and about 3.5% of the experimental atomisation energy.

!!! important

    Do not forget that the typical accuracy of LDA and GGA varies with the class of materials studied.
    Usually, LDA gives too small lattice parameters, by 1...3%, while GGA gives too large lattice parameters, by 1...3% as well,
    but there might be classes of materials for which the deviation is larger.
    See e.g. [[cite:Lejaeghere2014]].
