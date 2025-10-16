---
author: ECastiel
---

# Tutorial on DFT+DMFT with TRIQS/CT-HYB

## Introduction

In this tutorial, you will learn how to run a DFT+DMFT calculation using the interface
between ABINIT and TRIQS/CT-HYB. We show an example on iron.
TRIQS/CT-HYB is used as the impurity solver,
while ABINIT drives the calculation and connects to TRIQS/CT-HYB
as an external library.

Before starting, you should already know how to do a basic DFT calculation with ABINIT
(see [basic1](/tutorial/base1), [basic2](/tutorial/base2), [basic3](/tutorial/base3),
and [basic4](/tutorial/base4)). It also helps if you're familiar with PAW
(see [PAW1](/tutorial/paw1) and [PAW2](/tutorial/paw2)), [DFT+U](/tutorial/dftu),
and ABINIT's internal DMFT solvers (see [DMFT](/tutorial/dmft)).

We also recommend that you're comfortable with TRIQS/CT-HYB, and we'll assume you already
have at least a basic idea of how the DFT+DMFT method works.

The whole tutorial should take a few hours to complete.

[TUTORIAL_README]

## Presentation

ABINIT comes with its own internal DMFT code and solvers, but they aren't numerically exact
and rely on the density-density approximation. TRIQS/CT-HYB, on the other hand, can handle the
full interaction Hamiltonian.

Our interface makes it easy to run charge self-consistent DFT+DMFT calculations on real materials
with just a single input file. It is designed to give the accuracy and rigor
needed for computing structural properties, including a stationary DFT+DMFT implementation,
with the computation of the Baym-Kadanoff functional, the use of the exact double counting formula,
and an analytical evaluation of the high-frequency moments of the Green's function in a
self-consistent fashion.
That said, our interface does not handle toy models - TRIQS's Python API and built-in tools
already handle those efficiently.

Even though our interface with TRIQS/CT-HYB roughly follows the internal DMFT implementation
of ABINIT, there are a few important differences in input variables or available features, and
it is more than just swapping out the impurity solver.

You can find all the input variables related to DMFT listed [[varset:dmft|here]].

The input variables that are specific to the TRIQS/CT-HYB solver start with `dmft_triqs_`, and
many of them are direct copies of TRIQS/CT-HYB's own API inputs - so if you have used TRIQS
before, you should not get lost.

On the other hand, variables related to ABINIT's internal segment solver start with
`dmftctqmc_` or `dmftqmc_`. These don't apply here, so you can safely ignore them for this
tutorial.

Finally, variables that control the self-consistent cycle itself start with `dmft_`. Most of
these are shared between the internal DMFT and the TRIQS/CT-HYB interface, but some features
are available only in one or the other. The glossary clearly points out when that is the case.

If you want a good overview of the input variables relevant for the TRIQS/CT-HYB interface, we
recommend checking out the dedicated [[topic:DmftTriqsCthyb|topic]].

## Installation

First, you need to install [[https://triqs.github.io/triqs/latest/|TRIQS]] and
[[https://triqs.github.io/cthyb/latest/|TRIQS/CT-HYB]]. If you run into problems during
this step, please contact them, not us.

CT-HYB calculations can be very computationally demanding, so take your time when
compiling and pick an efficient BLAS implementation, as it can make a big difference
in performance.

After that, you need to reconfigure ABINIT to link it with TRIQS. Add the following line to your
configuration file:

    with_triqs="yes"

The library should be found automatically if the paths to the TRIQS and
TRIQS∕CT-HYB libraries are included in your environment variables. If you are using
the complex version of TRIQS/CT-HYB, also add:

    enable_triqs_complex="yes"

Next, make sure your options are recognized correctly and that ABINIT can
successfully link to TRIQS and compile a test code. When configuring, pay attention
to these lines:

    checking whether to enable the complex version of TRIQS... [yes|no]

and:

    checking whether the TRIQS library works... [yes|no]

If something goes wrong, carefully check the `config.log` file. Ask yourself: Why did the
test code fail to compile ? Was the library found and linked properly ? Is the test code
using syntax or functions from a more recent TRIQS/CT-HYB version ? This tutorial has
been written with version 4.0.x of TRIQS/CT-HYB.

Remember, you can override any automatically added compilation or linker flags using the
configuration files variables: `TRIQS_CPPFLAGS`, `TRIQS_CFLAGS`, `TRIQS_CXXFLAGS`,
`TRIQS_LDFLAGS`, and `TRIQS_LIBS`.

As a last resort, if you still can't solve the issue, post a question on the
[[https://discourse.abinit.org/|ABINIT user forum]].

Once everything is set, simply recompile ABINIT.

## Fatbands of Fe within DFT

Before jumping into the DFT+DMFT calculation for Fe, it is useful to first get some physical insight into
its electronic structure at the DFT level. This will help later when choosing some of the key DMFT
parameters.

We will start by computing the density of states (DOS), the band structure, and the fatbands. A fatband
plot shows how much of a given angular character (like $s$, $p$, or $d$) each Kohn-Sham wavefunction
has for every band and $k$-point - the thicker the line, the stronger that orbital's contribution.

Go to the tutorial input directory, make a new working folder, and copy the first input file:

```sh
cd $ABI_TESTS/tutoparal/Input
mkdir Work_dmft_triqs
cd Work_dmft_triqs
cp ../tdmft_triqs_1.abi .
```

{% dialog tests/tutoparal/Input/tdmft_triqs_1.abi %}

Launch the input file with ABINIT. For the rest of this tutorial, we will run everything on 4 CPUs,
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

Next, choose the correlated angular momentum $l$=0 ($s$), 1 ($p$), 2 ($d$)... for each atom type that will
be treated within DMFT. This works like in DFT+U using the variable [[lpawu]].
Currently, you can only choose one value per type of atom, and we only support [[lpawu]] $\le$ 3.

Set [[lpawu]]=-1 if you don't want DMFT on that atom type.

If you want DMFT on a subset of $d$ orbitals, use [[dmft_t2g]] (for $t_{2g}$ orbitals only) or
[[dmft_x2my2d]] (for the $d_{x^2-y^2}$ orbital only).

Our interface will then solve the impurity problem for each atom of this type separately, and interatomic
interactions are therefore not handled dynamically.

### Correlated orbital

After choosing the correlated channel (angular part), you need to define the radial part $u_l(r)$ of the
correlated orbitals $\frac{u_l(r)}{r} \, Y_{lm}(\hat{r})$. This is controlled
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
following the scheme specified by [[dmft_wanorthnorm]].

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

### Magnetism

Next, you will need to decide what kind of magnetism you want to include. If you set [[nsppol]]=1
and [[nspden]]=1, a paramagnetic solution is enforced by symmetrizing the two spin channels of the
Green's function. Even so, the impurity solver is solved with all spin channels, so local
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

ABINIT provides the exact double counting formula, which removes these contributions exactly. It is a
bit more involved to use, so we will cover it later in its own section.

There are also several simpler, commonly used approximate formulas, which you can select using the
variable [[dmft_dc]].

### Impurity Solver

We haven't yet told ABINIT to use the TRIQS/CT-HYB interface instead of its internal DMFT solvers.
This is done with [[dmft_solv]], which sets the impurity solver to use.

For TRIQS/CT-HYB, there are two relevant options:

  * [[dmft_solv]]=6: Uses TRIQS/CT-HYB but keeps only the density-density terms in the interaction
    tensor. This is much faster, though if that's what you need, a segment solver would be more
    efficient.
  * [[dmft_solv]]=7: Uses the full rotationally invariant Slater Hamiltonian, giving the most accurate
    treatment of interactions.

If instead you want to use ABINIT's internal DMFT implementation with its built-in solvers, check out the
dedicated [DMFT tutorial](/tutorial/dmft).

** How to run a charge self-consistent DFT+DMFT calculation: an example on Fe **

Let's now go through a simple charge self-consistent DFT+DMFT calculation on iron, applying DMFT to the
3d orbitals ([[lpawu]]=2).

Start by copying the second input file and running it with ABINIT. Note that this calculation can take
a while - expect at least a few minutes.

```sh
cp ../tdmft_triqs_2.abi .
mpirun -n 4 abinit tdmft_triqs_2.abi > log 2>&1
```

While it runs, open the input file and get familiar with the variables - many of them should already
look familiar.
We will go over the new ones shortly and explain how to choose their values. For even more details,
check out the variable glossary.

In this example, we run at high temperature ([[tsmear]]=3000 K) and use the density-density Hamiltonian
([[dmft_solv]]=6) to make things a bit lighter computationally.

When doing a charge self-consistent calculation, it is usually best to set [[dmft_iter]]=1 and only
adjust [[nstep]] to control convergence speed. That is because trying to reach self-consistency in the
DMFT loop before the charge density converges just wastes time.

We also perform a magnetic calculation ([[nsppol]]=2), but with a non-magnetic XC potential
([[usepawu]]=14), to avoid artificial magnetization from double counting. To help convergence, we
apply an initial shift between the two spin channels of the DMFT self-energy using [[dmft_shiftself]].

For the double-counting correction, we use the AMF scheme ([[dmft_dc]]=6), which works better for metals.
The correlated orbital is an atomic 3d orbital ([[dmft_orbital]]=1), projected onto Kohn-Sham bands 5 to
13 ([[dmftbandi]] and [[dmftbandf]]), which - as we saw earlier - carry most of the 3d spectral weight.

Once the calculation is finished, let's go through the log file and take a quick look at how the code
actually works. This will help you understand how to set the convergence parameters and make sense of
all the output files.
















