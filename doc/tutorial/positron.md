---
authors: JW, MT
---

# Electron-positron annihilation

## Calculation of positron lifetime and momentum distribution in silicon

This tutorial aims at showing how to perform **Two-Component Density-Functional
Theory (TCDFT)** calculations in the PAW framework to obtain the following
physical properties:

  * the positron lifetime in the perfect material,
  * the lifetime of a positron localized in a vacancy,
  * the electron-positron momentum distribution.

For the description of the implementation of TCDFT in ABINIT see [[cite:Wiktor2015]].

The user should be familiar with the four basic tutorials of ABINIT and the [first PAW tutorial](paw1).

This tutorial should take about 2 hours.

[TUTORIAL_READMEV9]

## Computing the positron lifetime in Si lattice

*Before beginning, you might consider to work in a different subdirectory as
for the other tutorials. Why not Work_positron?*

```sh
cd $ABI_TESTS/tutorial/Input
mkdir Work_positron
cd Work_positron
cp ../tpositron_x.files .   # You will need to edit this file.
cp ../tpositron_1.in .
```

The tutorial begins with a calculation of the positron lifetime in a silicon lattice.
In a perfect material the positron is delocalized. We can assume that its
density approaches zero and that it cannot affect the electron density. We will
perform a calculation in only two steps:

  1. Calculation of the ground-state electron density without the positron.
  2. Calculation of the ground-state positron density in the presence of the electron density from step 1.

The two densities are used to calculate the positron lifetime, which is
proportional to the inverse of the overlap of the electron and positron
densities. This 2-step calculation, considering the _zero-positron density
limit_, corresponds to the conventional scheme (CONV).

In the `tpositron_1.in` file, you will find two datasets.

{% dialog tests/tutorial/Input/tpositron_1.in %}

The first dataset is a standard ground-state calculation. The second one introduces a positron into
the system. You can see that in this case we set:

      positron2  1      ! Dataset 2 is a positronic GS calculation
      getden2    1      !   in presence of the previous electronic density
      kptopt2    0      !   Use only k=gamma point

      ixcpositron2  1   ! We are using the Boronski and Nieminen parametrization

Here we set [[positron]]=1, which corresponds to a positronic ground-state
calculation, considering that the electrons are not perturbed by the presence
of the positron (_zero-positron density limit_). The electron density is read from the file resulting from
dataset 1. As we consider the positron to be completely delocalized, we only
consider the &Gamma; point in the _Brillouin_ zone. The keyword [[ixcpositron]] selects the electron-
positron correlation functional and enhancement factor. In this calculation we
use the functional parametrized by Boronski and Nieminen [[cite:Boronski1986]], using the data provided by Arponen and Pajanne [[cite:Arponen1979]].

We can now run the calculation. In the directory
`~abinit/tests/tutorial/Input/Work_positron`, copy the files
`~abinit/tests/tutorial/Input/tpositron_x.files` and `tpositron_1.in`.

Then, issue:

    abinit < tpositron_x.files > log 2> err &

This calculation should only take a few seconds.

You can look at the *tpositron_1.out* file.
We find the positron lifetime calculated in the RPA limit:

    ########## Lifetime computation 2

     # Zero-positron density limit of Arponen and Pajanne provided by Boronski & Nieminen
       Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)
     # Enhancement factor of Boronski & Nieminen IN THE RPA LIMIT
       Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)

     Positron lifetime                         (ps)   =  2.22891945E+02

The lifetime of 223 ps agrees well with the value of 225 ps calculated with
the same number of valence electrons in [[cite:Wiktor2015]] and
with the experimental value of about 219 ps [[cite:Panda1997]].

!!! important

    If we had not used the "zero positron density limit" approximation
    (using, for example, another value of [[ixcpositron]]), we would theoretically
    have needed a box of infinite size for the positron to completely delocalise
    itself in the crystal. This can be avoided by the use of the [[posocc]] input
    parameter. Numerically, it is equivalent to calculate the density of 1
    positron in a box of size N and that of x positron in a box of size N * x.
    Thus, we can calculate the lifetime of a positron in a primitive cell by
    setting [[posocc]] to a small value (0.0001 ...). This value must obviously be tested...


## Computing the positron lifetime in a Si monovacancy

We will now perform a positron lifetime calculation for a monovacancy in
silicon in the conventional scheme (which we applied to the perfect lattice
previously). Note that when the positron localizes inside a vacancy, the _zero-
positron density limit_ does not apply anymore. However, in some cases, the
conventional scheme proved to yield results in agreement with experiments.

For the purpose of this tutorial, we generate a defect in a cell containing
only 16 atoms. This supercell is too small to get converged results, but the
calculation is relatively fast.

{% dialog tests/tutorial/Input/tpositron_2.in %}

You can now, issue (after having replaced *tpositron_1* by *tpositron_2* in the
*tpositron_x.files* file):

    abinit < tpositron_x.files > log 2> err &

Once the calculation is finished, look at the *tpositron_2.out* file.
Again, we look at the reported lifetime:

    ########## Lifetime computation 2

     # Zero-positron density limit of Arponen and Pajanne provided by Boronski & Nieminen
       Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)
     # Enhancement factor of Boronski & Nieminen IN THE RPA LIMIT
       Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)

     Positron lifetime                         (ps)   =  2.46936401E+02

We observe that when the positron localizes inside the vacancy, its lifetime
increases from 223 to 247 ps. This is because now the majority of the positron
density is localized in the vacancy region, where the electron density is
small. The overlap of the electron and positron densities is reduced, and the lifetime increased.

In the *Work_positron* directory, you will also find a *tpositron_2o_DS2_DEN_POSITRON*
file, containing the positron density. Visualizing this file (using e.g.
_cut3d_ and _XcrysDen_ or _VMD_ ) you can see that the positron is localized
inside the vacancy. You can see below how the positron (in red, isodensity at
30% of the maximum density) localized the silicon monovacancy looks like:

![](positron_assets/posdensity.png)


## Performing a self-consistent electron-positron calculation for a Si vacancy

We will now perform a self-consistent calculation of the positron and electron
densities. As this calculation will take a few minutes, you can already issue
(putting *tpositron_3.in* in *tpositron_x.files*):

    abinit < tpositron_x.files > log 2> err &

{% dialog tests/tutorial/Input/tpositron_3.in %}

This calculation is significantly longer than the previous one, because the
electron and positron steps will be repeated until the convergence criterion is reached.

In *tpositron_3.in* we only have one dataset and we set
[[positron]] = -10 to perform an automatic calculation of electrons and positron
densities. The convergence is controlled by [[postoldfe]] = 1d-5. This means
that we will repeat the electron and positron steps until the energy
difference between them is lower than 1d-5 Ha. This value should always be
larger than [[toldfe]]. In this calculation we still use [[ixcpositron]] = 1,
which means that we are using the GGGC scheme (see [[cite:Gilgien1994]] and [[cite:Wiktor2015]]

Once the calculation is finished, look at the positron lifetime in *tpositron_3.out*.

    ########## Lifetime computation 2

     # Zero-positron density limit of Arponen and Pajanne provided by Boronski & Nieminen
       Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)
     # Enhancement factor of Boronski & Nieminen IN THE RPA LIMIT
       Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)

     Positron lifetime                         (ps)   =  2.55617112E+02

Including the self-consistency increases the positron lifetime, because its
localization inside the vacancy becomes stronger when the positron and the
electron densities are allowed to relax.

## Relaxing the vacancy according to forces due to electrons and the positron

In addition to the self-consistency, the lifetime of a positron inside a
vacancy can be strongly affected by the relaxation of the atoms due to the
forces coming from both the electrons and the positron. You can already start
the relaxation of the vacancy by issuing:

    abinit < tpositron_4.files > log 2> err &

!!! important

    Don't forget to put *tpositron_4.in* in *tpositron_x.files*.

In this calculation we switched on the atomic relaxation by setting
[[ionmov]] = 2. We need to calculate forces to be able to move the atoms, so we
set [[optforces]] = 1. In the provided *tpositron_4.in* file, we only perform 4
relaxation steps ([[ntime]] = 4) to save time, but more steps would be needed to
converge the positron lifetime.

{% dialog tests/tutorial/Input/tpositron_3.in %}

Look at the positron lifetime in the RPA limit after each ionic step:

     Positron lifetime                         (ps)   =  2.55617112E+02
     Positron lifetime                         (ps)   =  2.56981105E+02
     Positron lifetime                         (ps)   =  2.81986785E+02
     Positron lifetime                         (ps)   =  2.82826327E+02

As the vacancy relaxes outwards, the positron lifetime increases. 4 steps were
not enough to relax the defect completely, as the lifetime still changes.
Indeed, setting [[ntime]] to 10 delivers:

    Positron lifetime                         (ps)   =  2.55617112E+02
    Positron lifetime                         (ps)   =  2.56981106E+02
    Positron lifetime                         (ps)   =  2.81986782E+02
    Positron lifetime                         (ps)   =  2.82826326E+02
    Positron lifetime                         (ps)   =  2.86660064E+02
    Positron lifetime                         (ps)   =  2.87040831E+02
    Positron lifetime                         (ps)   =  2.87284438E+02
    Positron lifetime                         (ps)   =  2.87360829E+02
    Positron lifetime                         (ps)   =  2.87302206E+02

Although the results at ionic steps 3 and 4 differ from each other by less
than one percent, they differ by more from the final result. The one percent
convergence is only reached at ionic step 5 and after.

Also, remember that the 16-atom supercell is not large enough to get converged
results. In Table IV of [[cite:Wiktor2015]] you can see converged
results of the positron lifetime of Si monovacancy within various methods.

## Computing the electron-positron momentum distribution (Doppler spectrum) of a Si lattice

In the last part of the tutorial we will calculate the electron-positron
momentum distribution (_Doppler spectrum_) of a silicon lattice in the conventional
scheme. This type of calculation is much more time and memory consuming than
the _lifetime_ calculation, as it is using the electron and positron
_wavefunctions_ (not only _densities_).

You can already issue (putting *tpositron_5.in* in *tpositron_x.files*):

    abinit < tpositron_5.files > log 2> err &

Now take a look at the input file *tpositron_5.in*.

{% dialog tests/tutorial/Input/tpositron_5.in %}

The momentum distribution calculation is activated by [[posdoppler]] = 1. You can also notice that instead
of having two datasets as in the first part of this tutorial, we now use the
automatic electron-positron loop and set [[posnstep]] = 2. This is done because
we need to have the full electron and positron wavefunctions in memory, which
is only the case when [[positron]] <= -10. Additionally, the momentum
distribution calculations require using a full k-point grid. In the input file we set:

    kptopt 0
    istwfk *1
    nkpt 8   # This corresponds to a 2x2x2 grid, denser grids may be needed to get converged spectra
    kpt
    0 0 0
    0 0 0.5
    0 0.5 0
    0.5 0 0
    0 0.5 0.5
    0.5 0 0.5
    0.5 0.5 0
    0.5 0.5 0.5

This grid is used in both electron and positron calculations, but only the
positron _wavefunction_ at the first point is taken in the momentum distribution
calculation, so the $\Gamma$ point should always be given first.

In the calculation of the momentum distribution, we need to include both _core_
and _valence_ electrons. The _wavefunctions_ of the core electrons are read from a
file (one per atom type), which needs to be provided. This _core WF file_ should
be named `<psp_file_name>.corewf` (where `<psp_file_name>` is the name of the
pseudo-potential (or PAW) file) or `corewf.abinit<ityp>` (where `<ityp>` is the
index of the atom type). _Core WF files_ can be obtained with the `atompaw` tool
(see [the tutorial on generating PAW datasets (PAW2)](paw2) ) by the use of the
`prtcorewf` keyword. You will find the core wavefunction file used in this calculation in
*$ABI_PSPDIR/Si.LDA-PW-paw.abinit.corewf*.

Once the calculation is complete, you can find a *tpositron_5o_DOPPLER* file
containing the _momentum distribution_ on the FFT grid. You can use the
*~abinit/scripts/post_processing/posdopspectra.F90* tool to generate 1D
projections (_Doppler spectra_) in (001), (011) and (111) directions and to
calculate the low- and high-momentum contributions to the
momentum distribution (so called `S` and `W` parameters, see [[cite:Wiktor2015]]).

## Studying the effect of the PAW dataset completeness

The positron lifetime and momentum distribution calculations within the PAW
method are very sensitive to the number of valence electrons in the **PAW
dataset**. It is due to the fact that it is not easy to describe the positron
_wavefunction_, tending to zero at the nucleus, using the electron atomic
orbitals. The **PAW basis set** in this case needs to be more complete than only
for describing the electron _wavefunctions_.

The simplest way to make the **PAW
dataset** more complete is to include `semicore electrons`. It is also possible to
add the `partial waves` corresponding to the `semicore electrons` in the basis
used only for the positron wave function description, while keeping the
initial number of valence electrons (as done in [[cite:Wiktor2015]]). However, this second method is less straightforward.

The previous calculations were done with only **4 valence electrons** (`3s` and `3p`).
We will now see what happens if we include the `2s` and `2p` states in the **PAW dataset**.
In *tpositron_12el_x.files* we have replaced the *Si.LDA-PW-paw.abinit*
dataset with *Si.12el.LDA-PW-paw.abinit*. We can now rerun the lifetime calculation:

    abinit < tpositron_12el_x.files > log 2> err

We now find the positron lifetime calculated in the RPA limit:

    ########## Lifetime computation 2

     # Zero-positron density limit of Arponen and Pajanne provided by Boronski & Nieminen
       Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)
     # Enhancement factor of Boronski & Nieminen IN THE RPA LIMIT
       Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)

     Positron lifetime                         (ps)   =  2.11481560E+02

This value is significantly lower than 223 ps achieved with 4 valence
electrons in the first step. **It is, therefore, very important to always test
the PAW dataset completeness for positron calculations**.

The PAW dataset completeness is even more important in the _Doppler spectra_
calculations. We will now recalculate the momentum distribution including 12
_valence electrons_ (using *tpositron_7.in* in *tpositron_12el_x.files*):

    abinit < tpositron_12el_x.files > log 2> err

Before processing the new *tpositron_7o_DOPPLER file*, you should copy files
`rho_001`, `rho_011`, `rho_111` from the fifth step to for instance `si4el_001`, `si4el_011` and `si4el_111`.
By plotting the _Doppler spectra_ in the (001) direction calculated with 4 and
12 valence electrons, you should obtain a figure like this:

![](positron_assets/doppler.png)

The dataset with 4 valence electrons is **not complete enough** to describe the
positron `wavefunction` around the nucleus. This is reflected in the
unphysically high probability at high momenta in the spectrum.

Further explanation of the influence of the PAW dataset on the _Doppler spectra_
can be found in [[cite:Wiktor2015]]. In case you need to generate
your own dataset for momentum distribution calculations, you can follow the [tutorial on generating PAW datasets (PAW2)](paw2).
