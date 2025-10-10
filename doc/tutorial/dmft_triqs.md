---
author: ECastiel
---

# Tutorial on DFT+DMFT with TRIQS/CT-HYB

## Introduction

In this tutorial, you will learn how to run a DFT+DMFT calculation using the interface
between Abinit and TRIQS/CT-HYB.
TRIQS/CT-HYB is used as the impurity solver,
while Abinit drives the calculation and connects to TRIQS/CT-HYB
as an external library.

Before starting, you should already know how to do a basic DFT calculation with Abinit
(see [basic1](/tutorial/base1), [basic2](/tutorial/base2), [basic3](/tutorial/base3),
and [basic4](/tutorial/base4)). It also helps if you're familiar with PAW
(see [PAW1](/tutorial/paw1) and [PAW2](/tutorial/paw2)), [DFT+U](/tutorial/dftu),
and Abinit's internal DMFT solvers (see [DMFT](/tutorial/dmft)).

We also recommend that you're comfortable with TRIQS/CT-HYB, and we'll assume you already
have at least a basic idea of how the DFT+DMFT method works.

The whole tutorial should take a few hours to complete.

[TUTORIAL_README]

## Presentation

Abinit comes with its own internal DMFT code and solvers, but they aren't numerically exact
and rely on the density-density approximation. TRIQS/CT-HYB, on the other hand, can handle the
full interaction Hamiltonian.

Our interface makes it easy to run charge self-consistent DFT+DMFT calculations on real materials
with just a single input file. It is designed to give the accuracy and rigor
needed for computing structural properties, including a stationary DFT+DMFT implementation,
with the computation of the Baym-Kadanoff functional and the use of the exact double counting formula.
That said, it is not meant for toy models - TRIQS's Python API and built-in tools
already handle those efficiently.

Some input variables are shared with Abinit's internal DMFT, while others are specific
to either the TRIQS interface or the internal implementation. You can find the full list in the
[[varset:dmft|DMFT]] and
[[varset:dmft_triqs|DMFT_TRIQS]] sections of the variable glossary.
Most TRIQS/CT-HYB specific variables directly correspond to TRIQS/CT-HYB's own inputs, so if you're
already familiar with it, you shouldn't get lost. Finally, we encourage you to read the
[[topic:DMFT|topic on DMFT]].

## Installation

First, you need to install [[https://triqs.github.io/triqs/latest/|TRIQS]] and
[[https://triqs.github.io/cthyb/latest/|TRIQS/CT-HYB]]. If you run into problems during
this step, please contact them, not us.

CT-HYB calculations can be very computationally demanding, so take your time when
compiling and pick an efficient BLAS implementation, as it can make a big difference
in performance.

After that, you need to reconfigure Abinit to link it with TRIQS. Add the following line to your
configuration file:

    with_triqs="yes"

The library should be found automatically if the paths to the TRIQS and
TRIQS∕CT-HYB libraries are included in your environment variables. If you are using
the complex version of TRIQS/CT-HYB, also add:

    enable_triqs_complex="yes"

Next, make sure your options are recognized correctly and that Abinit can
successfully link to TRIQS and compile a test code. When configuring, pay attention
to these lines:

    checking whether to enable the complex version of TRIQS... [yes|no]

and:

    checking whether the TRIQS library works... [yes|no]

If something goes wrong, carefully check the `config.log` file. Ask yourself: Why did the
test code fail to compile ? Was the library found and linked properly ? Is the test code
using syntax or functions from a more recent TRIQS/CT-HYB version ?

Remember, you can override any automatically added compilation or linker flags using the
configuration files variables: `TRIQS_CPPFLAGS`, `TRIQS_CFLAGS`, `TRIQS_CXXFLAGS`,
`TRIQS_LDFLAGS`, and `TRIQS_LIBS`.

As a last resort, if you still can't solve the issue, post a question on the
[[https://discourse.abinit.org/|Abinit user forum]].

Once everything is set, simply recompile Abinit.

## Basic parameters for DFT+DMFT

Several key parameters are needed for a DFT+DMFT calculation. You can check out the
[[topic:DMFT|DMFT topic]] later, especially the "compulsory" and "basic" sections, which
cover the main input variables you need to understand.

A DFT+DMFT charge self-consistent calculation usually involves two datasets: the first
is a well-converged DFT calculation, and the second a DFT+DMFT calculation that uses the
DFT output density from the first dataset as its starting point.

We won't go over DFT parameters or how to converge them here since that's already covered in the
DFT tutorials. Still, keep in mind that you may need to adjust your values and redo convergence
studies once DMFT is added. For example, you may need to converge your Kohn-Sham wavefunctions
better by increasing [[nline]] and [[nnsclo]]. Also, the electronic temperature (set via [[tsmear]])
becomes more important when DMFT is included, as temperature effects are explicitly treated
in the impurity solver, while they are often neglected in the DFT exchange-correlation functional.

Now let's go through the main DMFT variables ; you'll find more details in the variable glossary.

### Activate DMFT

To turn on DMFT, just set [[usedmft]]=1 in your input.

### Correlated electrons

Next, choose which atom types will be treated with DMFT and which angular momentum channel is correlated.
Currently, only one correlated angular momentum channel per atom type is allowed. This works like in
DFT+U using the variable [[lpawu]].

For each atom type, set [[lpawu]] to the correlated channel (e.g., 2 for $d$ orbitals, 3 for $f$ orbitals). Set [[lpawu]] = -1 if you don't want DMFT on that atom type.

Our interface lets you apply DMFT on multiple atoms, but each atom is treated separately, and interatomic
interactions are not handled dynamically. Besides, only [[lpawu]] $\le$ 3 is supported currently.

If you want DMFT on a subset of orbitals, use [[dmft_t2g]] (for $t_{2g}$ orbitals only) or [[dmft_x2my2d]] (for the $d_{x^2-y^2}$ orbital only).

### Correlated orbital

After choosing the correlated channel (angular part), you need to define the radial part (multiplied
by $r$), noted $u_l(r)$, of the local orbitals $\frac{u_l(r)}{r} \, Y_{lm}(\hat{r})$. This is controlled
by the variable [[dmft_orbital]].

Pick an orbital that captures as many electrons as possible for DMFT, focusing on those near half-filling
around the Fermi level where correlations matter most.

Since DMFT usually applies to localized electrons (like $3d$ or $4f$), an atomic orbital is often best.
By default, [[dmft_orbital]] corresponds to the most bound atomic orbital in your PAW dataset,
but you can also select a different one or provide a custom orbital from a file.

For technical reasons, this orbital must then be projected onto a finite set of Kohn-Sham bands,
defined by [[dmftbandi]] and [[dmftbandf]]. A wider energy window gives more localized orbitals,
while an infinite window brings you back to the original [[dmft_orbital]]. Make sure the window is
large enough to avoid overlaps with orbitals of neighboring atoms.

Because the orbitals are projected onto a finite energy window, they're no longer orthonormal.
To fix this, Wannier functions are built by orthonormalizing the projected orbitals,
following the schemes specified by [[dmft_wanorthnorm]].

### Interaction tensor

In DMFT, electrons don't interact through the bare Coulomb potential - instead, they feel a screened
interaction, which you'll need to specify. We assume spherical symmetry (the Slater parametrization),
so the full interaction tensor can be defined using just the Slater integrals $F^k$.

In practice, it is more common to use the average interaction $U$ and Hund's exchange $J$, set with
[[upawu]] and [[jpawu]] for each atom type. These are directly related to the Slater integrals.

By default, the higher-order Slater integrals are determined using the atomic ratios $F^k / F^2$,
but you can manually adjust them with [[f4of2_sla]] and [[f6of2_sla]] if needed.

### Self-consistent cycle

A charge self-consistent DFT+DMFT calculation involves two nested loops:

  * The outer DFT+DMFT loop, controlled by [[nstep]]
  * The inner DMFT loop, controlled by [[dmft_iter]]

This means the impurity solver is called a total of [[nstep]] x [[dmft_iter]] times. Each DMFT loop
runs at fixed electronic density, and that density gets updated [[nstep]] times in total - once every
[[dmft_iter]] DMFT iterations.

You can adjust how the self-energy is mixed using [[dmft_mxsf]].

Here's a schematic showing how the DFT+DMFT self-consistent cycle works:

![dmft_cycle_tuto](dmft_triqs_assets/dmft_cycle_tuto.png)

By default, the initial DMFT self-energy is set to the DFT double counting value.

If you are working with a magnetic system, you can speed up convergence by applying an initial
static shift between the two spin channels using [[dmft_shiftself]].

You can also restart from a previous DMFT run by loading an existing self-energy file with
[[getself]].

### Magnetism

Next, you will need to decide what kind of magnetism you want to include. If you set [[nsppol]]=1
and [[nspden]]=1, we enforce a paramagnetic solution by symmetrizing the two spin channels of the
Green's function. Even so, the impurity solver is still solved with all spin channels, so local
magnetic moments can still form.

If you set [[nsppol]]=2 and [[nspden]]=2, you enable collinear magnetism.

In some cases, it is better to keep the DFT exchange-correlation potential non-magnetic and let
magnetism emerge purely from DMFT. You can do that by setting [[usepawu]].

Finally, non-collinear calculations are also supported by our interface ([[nspinor]]=2 and [[nspden]]=4).
Here, we only treat spin-orbit coupling (SOC) as a perturbation at the DFT level - it is not included
in the DMFT model itself. SOC introduces spin off-diagonal terms, sometimes with an imaginary part,
and our interface can handle those.
However, a few features are missing in this case - for example, measuring the density matrix in the
solver ([[dmft_triqs_measure_density_matrix]]) isn't supported, which makes structural property
calculations impractical unless you neglect the off-diagonal terms. We emphasize that this limitation
comes from TRIQS∕CT-HYB itself, not from our implementation.

### Double counting

When combining DFT and DMFT, you run into the double counting problem - some of the local interactions
are already partly included at the DFT level, so they need to be subtracted and replaced by their exact
DMFT values.

Abinit provides the exact double counting formula, which removes these contributions exactly. It is a
bit more involved to use, so we will cover it later in its own section.

There are also several simpler, commonly used approximate formulas, which you can select using the
variable [[dmft_dc]].

### Impurity Solver

We haven't yet told Abinit to use the TRIQS/CT-HYB interface instead of its internal DMFT solvers.
This is done with [[dmft_solv]], which sets the impurity solver to use.

For TRIQS/CT-HYB, there are two relevant options:

  * [[dmft_solv]]=6: Uses TRIQS/CT-HYB but keeps only the density-density terms in the interaction
    tensor. This is much faster, though if that's what you need, a segment solver would be more
    efficient.
  * [[dmft_solv]]=7: Uses the full rotationally invariant Slater Hamiltonian, giving the most accurate
    treatment of interactions.

If instead you want to use Abinit's internal DMFT implementation with its built-in solvers, check out the
dedicated [DMFT tutorial](/tutorial/dmft).

## Electronic structure of Fe within DFT

Before jumping into the DFT+DMFT calculation for Fe, it is useful to first get some physical insight into
its electronic structure at the DFT level. This will help later when choosing some of the key DMFT
parameters.

We will start by computing the density of states (DOS), the band structure, and the fatbands. A fatband
plot shows how much of a given angular character (like $s$, $p$, or $d$) each Kohn-Sham wavefunction
has for every band and $k$-point - the thicker the line, the stronger that orbital's contribution.

Go to the tutorial input directory and make a new working folder:

```sh
cd $ABI_TESTS/tutoparal/Input
mkdir Work_dmft_triqs
cd Work_dmft_triqs
cp ../tdmft_triqs_1.abi .
```

{% dialog tests/tutoparal/Input/tdmft_triqs_1.abi %}

Launch the input file with Abinit. For the rest of this tutorial, we will run everything on 4 CPUs,
but feel free to use more if you can:

    mpirun -n 4 abinit tdmft_triqs_1.abi > log 2>&1

While it runs, take a look at the input file. It has two datasets:

  1. A DFT dataset to compute the DOS and the density.
  2. A non self-consistent dataset to compute the band structure along a specific $k$-path in the
     Brillouin zone (defined by [[kptbounds]]).

The fatbands are computed using the variable [[pawfatbnd]]. Here, it is set to 1, meaning fatbands are
computed for each angular momentum $l$. If you want orbitally resolved fatbands, set it to 2 instead.

Once the calculation finishes, you will find:

  * The DOS in the file: `tdmft_triqs_1o_DS1_DOS_AT0001`
  * The fatbands in files named like: `tdmft_triqs_1o_DS2_FATBANDS_at0001_Fe_is1_l000x` - one for
    each angular momentum.

For Fe, we are mainly interested in the $3d$ electrons, since they are the most strongly correlated.

You can plot the results using any plotting software you like - or simply use the provided Python
scripts:

```sh
python3 ../dmft_triqs_scripts/plot_dos.py
python3 ../dmft_triqs_scripts/plot_fatbands.py
```

You should get something like this for the DOS:

![dmft_cycle_tuto](dmft_triqs_assets/dos.png)

and this for the $d$ fatbands:

![dmft_cycle_tuto](dmft_triqs_assets/fatbands.png)

From the DOS, you can see that most of the $d$-electron spectral weight lies close to the Fermi level,
within an energy window of about 5 eV.

From the fatbands, you can tell that while the strongest $d$ character is indeed around the Fermi level,
it is not confined to the $3d$ bands (bands 6-10). There is significant hybridization with the $4s$
(band 5) and $4p$ (bands 11-13) states, and even some weaker but noticeable hybridization with higher
bands.

So when you move to the DMFT calculation, you will need to be careful in selecting which bands to
include, to make sure you capture as much of the $d$ spectral weight as possible.

** How to perform a charge self-consistent DFT+DMFT calculation: an example on Fe **

We will now show how to perform a basic charge self-consistent DFT+DMFT calculation on iron, applying DMFT to the 3d orbitals (lpawu=2).

Copy the second input, and launch it with Abinit. This might take quite some time (a few dozens minutes).

cp ../tdmft_triqs_2.abi .
mpirun -n 4 abinit tdmft_triqs_2.abi > log 2>&1

In the meantime, have a look at the input and get familiar with the input variables. You should already know some of them.
We will explain the ones that haven't been explained yet shortly, and show you how to choose their values. You can have more details by looking at the variable glossary.

We perform a high temperature calculation (tsmear = 3000 K), and choose the density-density Hamiltonian (dmft_solv=6) to make the calculation easier.
When performing a charge self-consistent calculation, we advise to set dmft_iter to 1 and only tweak nstep in order to speed up the convergence, as
achieving self-consistency on the DMFT self-energy when the density is not converged takes time for nothing.

Once the calculation is complete, you will see that a lot of output files specific to DMFT have been created.

 * tdmft_triqs_2o_DS2_CTQMC_DATA_iatom0001.h5 -> This file contains contains the last visited configuration of the last CT-QMC run. This file is automatically read
   at the beginning of each new call to the solver in order to restart the simulation from the previous configuration rather and speed up warmup. This behavior can be disabled
   by dmft_triqs_read_ctqmcdata if this causes some issues (e.g. if you have a violent change of behavior from one iteration to the other and the previous configuration is no
   longer suited). You can also use this file to initialize a new calculation with getctqmcdata.
 * tdmft_triqs_2o_DS2_CTQMC_HISTOGRAM_iatom0001.dat -> This file is printed whenever dmft_triqs_measure_density_matrix is activated. This prints the weight of each Fock state
   in the many-body state, directly sampled from the Monte-Carlo.
 * tdmft_triqs_2o_DS2_DMFTOCCND -> This file prints the DFT+DMFT occupation factors in the Kohn-Sham basis. Note that off-diagonal and imaginary components might appear.
 * tdmft_triqs_2o_DS2_DMFTORBITAL_itypat0001.dat -> This file prints the correlated orbital u_l(r) specified by dmft_orbital. Note that this is before projection, so this
   is not yet the Wannier function that is used in practice.
 * tdmft_triqs_2o_DS2_Gtau_diag_DLR_iatom0001.dat -> This file corresponds to the diagonal elements of the imaginary-time DLR representation of the Green's function
 * tdmft_triqs_2o_DS2_Gtau_diag_iatom0001.dat -> This file corresponds to the diagonal elements of the binned imaginary-time Green's function directly sampled from the Monte-Carlo.
 * tdmft_triqs_2o_DS2_Hybridization_diag_iatom0001.dat -> This file corresponds to the diagonal elements of the imaginary-time hybridization function.
 * tdmft_triqs_2o_DS2_Self-omega_iatom0001_isppol1 -> This file corresponds to the imaginary-frequency self-energy, and is read at the beginning of each new DMFT loop in order
   to perform a charge self-consistent calculation and start from the previous self-energy. This behavior can be disabled via dmft_rslf, though we don't see any point in doing that.
   This file can be used for restart and be read to initialize a brand new calculation via getself.

Let's now browse across the log file, and briefly explain how the code works and how to choose the relevant computational parameters. Note that we are voluntarily very agressive
towards the user, and we choose not to set a default value for the parameters we deem important. Indeed, we feel that most of the discrepancies that can be found in the various DFT+DMFT
results in the literature are due to poorly converged results, thus we want to force you to think and choose the values of your parameters wisely, according to your system.















