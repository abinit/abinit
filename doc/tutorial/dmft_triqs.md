---
author: ECastiel
---

# Tutorial on DFT+DMFT with TRIQS/CT-HYB

## Introduction

In this tutorial, you will learn how to run a DFT+DMFT calculation on iron using the interface
between ABINIT and TRIQS/CT-HYB.
TRIQS/CT-HYB is used as the impurity solver,
while ABINIT drives the calculation and connects to TRIQS/CT-HYB
as an external library.

Before starting, you should already know how to do a basic DFT calculation with ABINIT
(see [basic1](/tutorial/base1), [basic2](/tutorial/base2), [basic3](/tutorial/base3),
and [basic4](/tutorial/base4)). It also helps if you are familiar with PAW
(see [PAW1](/tutorial/paw1) and [PAW2](/tutorial/paw2)), [DFT+U](/tutorial/dftu),
and ABINIT's internal DMFT solvers (see [DMFT](/tutorial/dmft)).

We also recommend that you are comfortable with TRIQS/CT-HYB, and we will assume you already
have at least a basic idea of how the DFT+DMFT method works.

The whole tutorial should take a few hours to complete.

[TUTORIAL_README]

## Presentation

ABINIT comes with its own internal DMFT code and solvers, but they are not numerically exact
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

The fatbands are computed using the variable [[pawfatbnd]]. Here, it is set to $1$, meaning fatbands are
computed for each angular momentum $l$. If you want orbitally resolved fatbands, set it to $2$ instead.

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
for each system you study.

Alright, let's dive in. Start by copying the second input file and running it with ABINIT. Feel free
to use more CPUs if you have them available:

```sh
cp ../tdmft_triqs_2.abi .
mpirun -n 4 abinit tdmft_triqs_2.abi > log2 2>&1
```

{% dialog tests/tutoparal/Input/tdmft_triqs_2.abi %}

This calculation might take a few minutes, so while it runs, let's take a look at the input file and
see what's going on.

A charge self-consistent DFT+DMFT calculation usually consists of two datasets:

  1. The first one is a DFT calculation (that should be in principle well converged).
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
complex as $\beta$ increases, we are choosing a high temperature ([[tsmear]]$=3000 \, \mathrm{K}$) to keep
things simpler and faster for this example.

For the second dataset, we start by activating DMFT, setting [[usedmft]]$=1$.

### Correlated Electrons

Next, we need to choose which angular momentum channel ($s$,$p$,$d$,...) we want to treat as
correlated within DMFT. Each channel corresponds
to a value of $l=0,1,2$..., and we select it using the variable [[lpawu]], just like in
DFT+U.

For each atom type, you specify the value of $l$ you want to apply DMFT to. If you do not
want to apply DMFT on a given atom type, simply set [[lpawu]]$=-1$.

Each correlated channel will then correspond to $2 \times (2l+1)$ orbitals per correlated
atom - the factor of $2$ accounts for spin degeneracy, and the $2l+1$ comes from the different
magnetic quantum numbers $m = -l,...l$.
Keep in mind that the complexity of the impurity problem increases exponentially with this
number of orbitals, so choosing your correlated channel wisely is important.

!!! note

    At the moment, only values of [[lpawu]] $\le 3$ are supported - anything higher would be computationally impractical anyway.

In our iron example, the localized orbitals that show strong correlation effects are the $3d$
orbitals. Since there is only one atom type (Fe), we simply set [[lpawu]]$=2$.

### Correlated Orbital

Now that we have defined which angular momentum channel we are treating as correlated, we need
to specify the radial part of the correlated orbital - the reduced radial wavefunction $u_l(r)$.

Our convention for the local orbitals follows ABINIT's internal PAW convention:
$\frac{u_l(r)}{r} Y_{lm}(\hat{r})$ where $Y_{lm}(\hat{r})$ are real spherical harmonics.

The choice of the radial wavefunction is made through the variable [[dmft_orbital]]. For iron, the
$3d$ orbitals remain quite localized and keep their atomic character, so we just take the
lowest-energy atomic orbital from the PAW dataset (truncated at the PAW radius). That is the default
choice, [[dmft_orbital]]$=1$, but you can pick any orbital from your PAW dataset or even provide a
custom one from a file if needed.

!!! tip

    More generally, you should pick a radial orbital that captures as much spectral weight as possible in DMFT, especially near the Fermi level, since that's where electron correlations have the biggest impact.

When the calculation starts, you can see how the radial orbital is defined in the `log2` file, right
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

### Energy window

Unfortunately, we cannot directly compute the lattice Green's function over the full Hilbert space - it
would be far too computationally demanding. Instead, we restrict ourselves to a smaller subspace
$\mathcal{H}$, made up of low-energy basis elements. We then downfold this Green's function, which is
projected onto $\mathcal{H}$, onto the local orbitals defined by [[dmft_orbital]] in order to retrieve
the local Green's function.

But notice: this is actually the same as downfolding the full lattice Green's function onto the
projection of [[dmft_orbital]] onto $\mathcal{H}$.

Therefore, the DMFT local orbitals that are effectively handled in the code are not exactly the ones
defined by [[dmft_orbital]], but rather their projection on $\mathcal{H}$.

In our implementation, the lattice Green's function is represented in the Kohn-Sham basis, and
our low-energy Hilbert space is built from all Kohn-Sham wavefunctions whose band indices fall between
[[dmftbandi]] and [[dmftbandf]] - the range that defines the correlated energy window.

!!! note

    We would like to emphasize, however, that when we compute the trace of the Green's function - for instance, to locate the Fermi level - we include all bands (from $1$ to [[nband]]). This differs from many implementations that limit themselves strictly to the correlated energy window.

Since Kohn-Sham wavefunctions form a complete basis, you could in theory recover the original orbital
defined by [[dmft_orbital]] by setting [[dmftbandi]] to $1$ and converging [[dmftbandf]] to infinity.
In practice, however, this is not feasible - especially when your atomic orbital is truncated -
because it would involve an unmanageable number of bands.

The choice of your energy window is therefore critical. It directly governs both the spatial spread of
your orbital and your computation time.

  * A wider energy window yields a more localized orbital, but also increases hybridization with the
    bath. Since we use a hybridization-expansion impurity solver, a larger hybridization significantly
    slows down the solver and the evaluation of the Green's function.
  * A narrower energy window may leave out part of the spectral weight of your target orbital
    [[dmft_orbital]]. In this case, the projected orbital may become too extended, potentially
    overlapping with neighboring atoms. That would violate a key approximation of the method - that
    interatomic elements of the Green's function can be neglected.

This is extremely important since if this approximation breaks, the very important identity
Downfold(Upfold) = Identity might break down, which can lead to inaccuracies and unstability.
You can check the deviation in the `log2` file:

```sh
  == Check downfold(upfold)=identity ==

--- !WARNING
src_file: m_matlu.F90
src_line: 1373
message: |

    Differences between Downfold(Upfold) and
       Identity is too large:
         0.1284E-03 is larger than  0.1000E-03
...
```

As an exercise, you can check that this deviation becomes lower whenever you increase the energy window.

!!! tip

    Because Kohn-Sham wavefunctions vary between systems, the projection of [[dmft_orbital]] onto your subspace $\mathcal{H}$ will also vary, especially if the energy window is too narrow. To reliably compare energies between systems, you have to use the same orbital for all calculations, and converge your results with respect to [[dmftbandi]] and [[dmftbandf]]. See the section on the exact double counting formula for more guidance.

As shown earlier from the fatband analysis, the main spectral weight for the $d$ orbitals lie between
bands $5$ and $13$, so we use this as our energy window.

Once the orbitals have been projected onto this finite Hilbert space, they are no longer orthonormal.
You can inspect this in the `log2` file:

```sh
 == The DMFT orbitals are now projected on the correlated bands

 == Check: Downfolded Occupations and Norm of unnormalized projected orbitals

  ------ Symmetrized Occupations

   -------> For Correlated Atom 1

          -- polarization spin component 1
        0.50385   0.00000   0.00000   0.00000  -0.00000
        0.00000   0.50385  -0.00000   0.00000   0.00000
        0.00000   0.00000   0.45679   0.00000   0.00000
        0.00000   0.00000   0.00000   0.50385   0.00000
       -0.00000   0.00000   0.00000   0.00000   0.45679

  ------ Symmetrized Norm

   -------> For Correlated Atom 1

          -- polarization spin component 1
        0.80454   0.00000   0.00000   0.00000   0.00000
        0.00000   0.80454   0.00000   0.00000  -0.00000
        0.00000   0.00000   0.80154   0.00000   0.00000
       -0.00000   0.00000   0.00000   0.80454  -0.00000
        0.00000  -0.00000   0.00000   0.00000   0.80154
```

Here, the occupation and norm matrices (prior to orthonormalization) are
shown. All local operators are printed in matrix form in the `log` file,
sorted by spin and increasing magnetic quantum number $m$.

By default, all quantities are output in the real spherical harmonics (cubic) basis. The only exception
occurs when data is passed to the CT-HYB solver, which requires a rotation to the CT-QMC basis (see the
corresponding subsection for details). Wherever a basis change occurs, it is clearly indicated in the
`log`.

As shown above, the norms of the projected orbitals are slightly below the original atomic value
(which wasn't normalized to begin with). Try increasing [[dmftbandf]], and you will see the norms
approach the true atomic value of $0.8514$.

Finally, we promote the projected orbitals to proper Wannier functions by orthonormalizing them via the scheme specified by [[dmft_wanorthnorm]].

### Interaction tensor

In a solid, the correlated electrons do not interact via the bare Coulomb potential - their interaction is screened by all the other electrons in the system. Because of this,
you need to provide the screened interaction tensor $U_{ijkl}$.

In our implementation, we use the Slater parametrization, which expresses the full interaction tensor in terms of just a few Slater integrals $F^k$ ($k = 2i$, with $i = 0, l$).

By default, the ratios $F^k / F^2$ (for $k \ge 4$) are kept at their atomic values, which is usually a good approximation for $3d$ transition metals and rare-earth elements.
If you want to adjust them manually, you can do so via the variables [[f4of2_sla]] and [[f6of2_sla]].

Once these ratios are fixed, the two remaining Slater integrals can be determined uniquely from the average screened interaction $U$ and the screened Hund’s exchange $J$.
These parameters have clear physical meanings — $U$ sets the overall interaction strength, while $J$ controls the tendency of electrons to align their spins.
You can set them for each atom type using [[upawu]] and [[jpawu]].

There are several ways to compute these parameters directly within ABINIT — for example, using [cRPA](/tutorial/ucalc_crpa) or [linear response](/tutorial/lruj) - but those
are beyond the scope of this tutorial. Here, we will simply treat $U$ and $J$ as input parameters.

Keep in mind that these are matrix elements of the screened potential, meaning their values depend on the shape of the local orbitals.
So, if you change your energy window, you may need to readjust $U$ and $J$.

For this example, we will use the same parameters as in [[cite:Han2018]] —
$U = 5.5$ eV and $J = 0.84$ eV — since their chosen energy window is similar to ours.

### Self-consistent cycle

A charge self-consistent DFT+DMFT calculation involves two nested loops:

  * The outer DFT+DMFT loop, whose number of iterations is controlled by [[nstep]]
  * The inner DMFT loop, whose number of iterations is controlled by [[dmft_iter]]

This means the impurity solver is called a total of [[nstep]] $\times$ [[dmft_iter]] times. Each DMFT loop is performed at fixed electronic density, and that density gets updated
[[nstep]] times in total - once every [[dmft_iter]] DMFT iterations.

Here’s a schematic showing how the DFT+DMFT self-consistent cycle works:

![dmft_cycle_tuto](dmft_triqs_assets/dmft_cycle_tuto.png)

When doing a charge self-consistent calculation, it is usually best to set [[dmft_iter]]$=1$ and only adjust [[nstep]] to control convergence speed.
That is because trying to reach self-consistency in the DMFT loop before the charge density converges just wastes time.

### Magnetism

The next step is to decide what type of magnetism you want to include. Here, we choose [[nsppol]]$=$[[nspden]]$=1$, which enforces a collinear paramagnetic solution by
symmetrizing the two spin channels of the Green’s function. Even so, the impurity solver is solved with all spin channels, meaning local magnetic moments can still form.

Collinear magnetism can be enabled with [[nsppol]]$=$[[nspden]]$=2$. In this case, it can sometimes be better to keep the DFT exchange-correlation potential non-magnetic and let magnetism emerge purely from DMFT. You can do that by setting [[usepawu]]$=10$. Otherwise, set [[usepawu]]$=14$ for a magnetic XC potential.

Finally, non-collinear calculations are also supported by our interface ([[nspinor]]$=2$ and [[nspden]]$=4$).

### Double counting

When combining DFT and DMFT, you run into the double counting problem - some of the local interactions are already partly included at the DFT level, so they need to be subtracted and
replaced by their exact DMFT values.

ABINIT provides the exact double counting formula, which removes these contributions exactly. It is a bit more involved to use, so we will cover it later in its own section.

There are also several simpler, commonly used approximate formulas, which you can select using the variable [[dmft_dc]]. Here, we choose [[dmft_dc]]$=6$, corresponding to the AMF correction, better suited for metals.

### Impurity solver

We haven’t yet told ABINIT to use the TRIQS/CT-HYB interface instead of its internal DMFT solvers. This is done with [[dmft_solv]], which sets the impurity solver to be used.

For TRIQS/CT-HYB, there are two relevant options:

  * [[dmft_solv]]$=6$: Uses TRIQS/CT-HYB but keeps only the density-density terms in the interaction tensor. This is much faster, though if that’s what you need, a segment solver would be
      more efficient.
  * [[dmft_solv]]$=7$: Uses the full rotationally invariant Slater Hamiltonian, giving the most accurate treatment of interactions.

If instead you want to use ABINIT’s internal DMFT implementation with its built-in solvers, check out the dedicated [DMFT tutorial](/tutorial/dmft).

### Matsubara Mesh

In our implementation, the Green’s function is represented in imaginary frequencies as $G(i\omega_n)$, and is computed exactly on the first [[dmft_triqs_n_iw]] = $n_{i\omega}$ Matsubara frequencies, defined as
\begin{equation*}
    \omega_n = (2n + 1)\frac{\pi}{\beta}
\end{equation*}

During the calculation, several steps require integrating over all Matsubara frequencies.
However, directly summing:
\begin{equation*}
    \sum_{n=-n_{i \omega}}^{n_{i \omega}-1} G(i\omega_n) e^{i 0^{+}}
\end{equation*}
(with $e^{i0^+}$ a small damping factor) would be impractical — it would take far too long to converge.

To make this efficient, we use a moment expansion of the Green’s function up to fifth order:

\begin{equation*}
    \sum_{n=-\infty}^{\infty} G(i\omega_n) e^{i 0^{+}} \approx \sum_{n=-n_{i\omega}}^{n_{i\omega}-1} \left[ G(i\omega_n) - \sum_{j=1}^{5} \frac{G_j}{(i\omega_n)^j} \right] + \sum_{n=-\infty}^{\infty} \frac{G_j}{(i\omega_n)^j} e^{i 0^{+}}
\end{equation*}
where $G_j$ are the moments of the Green’s function.

The first term is computed numerically, while the second term is handled analytically.
Since the moments $G_j$ are obtained analytically in our implementation,
this approach is both efficient and guaranteed to converge to the exact result.

To ensure consistency, the code automatically compares this summation to the corresponding value for the DFT Green’s function, where the sum reduces to a Fermi–Dirac occupation factor.
If the difference is too large, the code stops and raises an error — this is a helpful way to check whether your [[dmft_triqs_n_iw]] value is properly chosen.

You can find this check in your `log2` file. If everything is working correctly, you’ll see lines similar to:

```sh
** Differences between occupations from DFT Green's function and
   Fermi-Dirac occupations are small enough:
   0.1460E-11 is lower than  0.1000E-03
```

and:

```sh
== Compute DFT Band Energy terms
     Differences between band energy with Fermi-Dirac occupations
     and occupations from DFT Green's function is:   -0.000000
     which is smaller than          0.00001
```

As you can see, this integration method is highly efficient — in this example, only about $500$ frequencies were needed to reach excellent accuracy.

!!! tip

    When you change the temperature, remember to increase the number of Matsubara frequencies linearly with $\beta$ to maintain the same effective energy range.

### CT-QMC basis

In the SCF cycle, we work in the cubic basis.
However, even though [[dmft_solv]]$=7$ is fully rotationally invariant and
basis independent, the choice of solver basis can in practice have a
huge impact on computation time.

Strong off-diagonal components in the hybridization function can make the sign problem worse. On the other hand, TRIQS/CT-HYB relies on conserved quantities (like angular momentum
or spin) to partition the local Hilbert space of size $2^{2 \times (2l+1)}$, which reduces the size of the matrices that need to be handled during the simulation. These conserved
quantities are detected automatically, and their number depends strongly on the choice of basis.

TRIQS/CT-HYB prints the number of subspaces it finds:

```sh
  Using autopartition algorithm to partition the local Hilbert space
  Found 1024 subspaces.
```

A higher number of subspaces (thus with smaller dimensions) means smaller matrices, which reduces computation time. In this example, the local Hilbert space size is $2^{10}=1024$
($l=2$), and we see $1024$ subspaces, each of dimension $1$. Thus, we only handle scalars, and do not even need a matrix solver, as we make the density-density approximation in this
example. You can check as an exercise that the matrix sizes will increase by setting [[dmft_solv]]$=7$, depending again on your basis choice.

So, there is a trade-off: you want to reduce off-diagonal components and maximize the number of subspaces, which isn't always easy. Fortunately, our interface provides many basis
options that can be set via [[dmft_triqs_basis]] - check the glossary for all available choices.

!!! note

    Even though our interface handles the full off-diagonal hybridization and Green's function, some features are not available in this case. For instance, density matrix sampling ([[dmft_triqs_measure_density_matrix]]) cannot be performed, making the sampling of static observables (like energy or electron number) much noisier. We wish to emphasize that this limitation comes from TRIQS/CT-HYB itself, not our interface.

To simplify things, you can set the off-diagonal components to zero in the CT-QMC basis via [[dmft_triqs_off_diag]]$=0$. Be careful, as the calculation is no longer numerically exact in this case.

In this example, we choose to stay in the cubic basis ([[dmft_triqs_basis]]$=0$) and neglect the off-diagonal components, since they vanish by symmetry anyway. The code then
automatically detects and prints the most compact block structure of the hybridization function and electronic levels:

```sh
   == Searching for the most optimal block structure of the electronic levels and the hybridization

   == Solving impurity model for atom 1, where there are 10 blocks

  --> Block 0 contains flavor 0
  --> Block 1 contains flavor 1
  --> Block 2 contains flavor 2
  --> Block 3 contains flavor 3
  --> Block 4 contains flavor 4
  --> Block 5 contains flavor 5
  --> Block 6 contains flavor 6
  --> Block 7 contains flavor 7
  --> Block 8 contains flavor 8
  --> Block 9 contains flavor 9

   == Schematic of the block structure

      0  .  .  .  .  .  .  .  .  .
      .  1  .  .  .  .  .  .  .  .
      .  .  2  .  .  .  .  .  .  .
      .  .  .  3  .  .  .  .  .  .
      .  .  .  .  4  .  .  .  .  .
      .  .  .  .  .  5  .  .  .  .
      .  .  .  .  .  .  6  .  .  .
      .  .  .  .  .  .  .  7  .  .
      .  .  .  .  .  .  .  .  8  .
      .  .  .  .  .  .  .  .  .  9
```

This allows TRIQS/CT-HYB to discard some irrelevant moves, and to optimize the computation of the determinant of the hybridization matrix. In our case, there are
$10$ blocks of size $1$, as there are no off-diagonal components.

### CT-HYB parameters

The key controls for the TRIQS/CT-HYB solver are [[dmft_triqs_n_cycles]], [[dmft_triqs_length_cycle]],
[[dmft_triqs_n_warmup_cycles_init]] and [[dmft_triqs_n_warmup_cycles_restart]]. These parameters are
crucial for both performance and accuracy and should be carefully converged.

Every CT-QMC simulation begins with a warmup phase, where the Markov chain hasn't reached equilibrium
yet, so no measurements are taken. Each warmup phase consists of a number of cycles, and each cycle
contains [[dmft_triqs_length_cycle]] sweeps.

There are two kinds of warmup:

  1. The initial warmup, triggered at the first DMFT iteration, when we start from scratch (i.e. an empty
     configuration) and no initial configuration is provided by the user via [[getctqmcdata]]. It is
     controlled by [[dmft_triqs_n_warmup_cycles_init]], can be quite long, and must be carefully
     converged.
  2. The restart warmup, used from the second DMFT iteration onward, and at the first iteration if
     a configuration file is provided via [[getctqmcdata]]. We automatically restart from the last saved
     CT-QMC configuration (`tdmft_triqs_2o_DS2_CTQMC_DATA_iatom0001.h5`), and thus the warmup can be
     much shorter, driven by [[dmft_triqs_n_warmup_cycles_restart]] - which can often be set to $0$.

The restart behavior can be disabled by [[dmft_triqs_read_ctqmcdata]] as it can cause some issues in
some rare occasions - at very low temperature for instance, when the configuration weights become
extremely low.

After warmup, the solver moves into the measurement phase, which lasts [[dmft_triqs_n_cycles]] cycles.
During this phase, observables are measured at the end of each cycle.

One key aspect to keep in mind is that only decorrelated/independent measurements matter for statistics.
That is why [[dmft_triqs_length_cycle]] should in theory be set to match the autocorrelation time.

If it is too small, you are measuring correlated samples - wasting time on measurement - and if it is
too big, you are discarding valid statistics - getting fewer meaningful samples.

In practice, to set this value, we advise you to keep increasing it until the time spent in measurement printed here:

```sh
  Total measure time                         | 0.0910871
```

is negligible compared to the simulation time. Then, check that you haven't lost statistics on your quantity of interest.

!!! tip

    We do not advise to check the autocorrelation time printed by TRIQS/CT-HYB (unlike what is said in this tutorial), as an automatic measurement of this quantity is highly unreliable.

At the end of a CT-HYB simulation, you will see a summary of the different Monte-Carlo moves and their
acceptance rates:

```sh
  Move set Insert two operators: 0.401734
  Move set Remove two operators: 0.404321
  Move  Shift one operator: 0.518682
```

These numbers give important clues about the sampling efficiency. If a move has a very low acceptance
rate and is not essential for ergodicity, you can disable it, but be sure of what you are doing.

You will find detailed explanations of the moves on the TRIQS/CT-HYB
[[https://triqs.github.io/cthyb/latest/basicnotions/moves.html|website]].

The key parameters to control the proposal distribution are [[dmft_triqs_move_double]], [[dmft_triqs_move_shift]] and [[dmft_triqs_pauli_prob]].

















