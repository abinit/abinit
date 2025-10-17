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
within an energy window of about $5$ eV.

From the fatbands, you can tell that while the strongest $d$ character is indeed around the Fermi level,
it is not confined to the $3d$ bands (bands $6-10$). There is significant hybridization with the $4s$
(band $5$) and $4p$ (bands $11-13$) states, and even some weaker but noticeable hybridization with higher
bands.

So when you move to the DMFT calculation, you will need to be careful in selecting which bands to
include, to make sure you capture as much of the $d$ spectral weight as possible.

## Charge self-consistent DFT+DMFT calculation on Fe

We will now walk you through the basics of running a charge self-consistent DFT+DMFT calculation,
using iron as an example. Along the way, we will introduce the main input variables you need to know,
show you how our implementation works, and go over the output files. We will also look at the main
convergence parameters and share some tips on how to choose their values wisely.

Before we start, an important note: do not use the values from our example input file as reference for
your own convergence parameters. We intentionally picked small values here to keep the runtime
reasonable. Our implementation is purposely strict about this - we do not provide default values for
many parameters, because we want users to choose them consciously.

Indeed, we feel that most of the inconsistencies found in DFT+DMFT results across the literature come
from poorly chosen convergence parameters. So, always perform proper convergence tests, and redo them
for each system you study - there is no one-size-fits-all setup.

Alright, let's dive in. Start by copying the second input file and running it with ABINIT. Feel free
to use more CPUs if you have them available:

```sh
cp ../tdmft_triqs_2.abi .
mpirun -n 4 abinit tdmft_triqs_2.abi > log 2>&1
```

{% dialog tests/tutoparal/Input/tdmft_triqs_2.abi %}

This calculation might take a few minutes, so while it runs, let's take a look at the input file and
see what's going on.

A charge self-consistent DFT+DMFT calculation usually consists of two datasets:

  1. The first one is a well-converged DFT calculation.
  2. The second one is a DFT+DMFT calculation, which starts from the DFT output density of
     the first dataset.

We won't revisit the DFT parameters or how to converge them here - that is already covered in the
DFT tutorials. However, remember that once you include DMFT, you'll likely need to readjust your
parameters and redo your convergence studies. For instance, you will often need to converge your
Kohn-Sham wavefunctions more tightly by increasing [[nline]] and [[nnsclo]].

Also, the electronic temperature (set via [[tsmear]]$=1/\beta$) becomes much more important when
DMFT is active, since temperature effects are explicitly included in the impurity solver -
unlike in the DFT exchange-correlation functional, where they are usually neglected.

Here, we use a CT-HYB solver, which represents the Green's function $G(\tau)$ on the imaginary
time axis, with $\tau$ ranging from $0$ to $\beta$. Because the impurity problem becomes more
complex as $\beta$ increases, we are choosing a high temperature ([[tsmear]]$=3000K$) to keep
things simpler and faster for this example.

For the second dataset, we start by activating DMFT by setting [[usedmft]]$=1$.

### Correlated Electrons

Next, we need to choose which angular momentum channel - that is, which type of orbital
($s$,$p$,$d$,...) - we want to treat as correlated within DMFT. Each channel corresponds
to a value of $l=0,1,2$..., and we select it using the variable [[lpawu]], just like in
DFT+U.

For each atom type, you specify the value of $l$ you want to apply DMFT to. If you do not
want to apply DMFT on a given atom type, simply set [[lpawu]]$=-1$.

Each correlated channel will then correspond to $2 \times (2l+1)$ orbitals per correlated
atom - the factor of $2$ accounts for spin degeneracy, and the $2l+1$ comes from the different
magnetic quantum numbers $m = -l,...l$.
Keep in mind that the complexity of the impurity problem increases exponentially with this
number of orbitals, so choosing your correlated channel wisely is important.

At the moment, only values of [[lpawu]] $\le 3$ are supported - anything higher would be
computationally impractical anyway.

In our iron example, the localized orbitals that show strong correlation effects are the $3d$
orbitals. Since there is only one atom type (Fe), we simply set [[lpawu]]$=2$.

### Correlated Orbital

Now that we have defined which angular momentum channel we are treating as correlated, we need
to specify the radial part of the correlated orbital - the reduced radial wavefunction $u_l(r)$.

Our convention for the local orbitals follows ABINIT's internal PAW convention:
$\frac{u_l(r)}{r} Y_{lm}(\hat{r})$ where $Y_{lm}(\hat{r})$ are real spherical harmonics.

The choice of the radial wavefunction is made through the variable [[dmft_orbital]]. For iron, the
$3d$ orbitals remain quite localized and keep their atomic character, so we usually just take the
lowest-energy atomic orbital from the PAW dataset (truncated at the PAW radius). That is the default
choice, [[dmft_orbital]]=1, but you can pick any orbital from your PAW dataset or even provide a
custom one from a file if needed.

When the calculation starts, you can see how the radial orbital is defined in the `log` file, right
at the beginning of the self-consistent cycle:

```sh
   =====  Build DMFT radial orbital for atom type 1 ========

   Using atomic orbital number 1 from PAW dataset
   Squared norm of the DMFT orbital: 0.8514
```

This lets you check the normalization of the orbital. The radial function $u_l(r)$ is also written
to the file `tdmft_triqs_2o_DS2_DMFTORBITAL_itypat0001.dat`. You can plot it with your favorite tool
 - for instance, all DMFT-related outputs are compatible with xmgrace. If you plot this file, you
should clearly see the $3d$ atomic orbital from your PAW dataset:

![dmft_cycle_tuto](dmft_triqs_assets/3d_orbital.png)
















