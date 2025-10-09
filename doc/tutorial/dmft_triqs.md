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

The whole tutorial should take about one hour to complete
(faster if you can use multiple processors).

[TUTORIAL_README]

## Presentation

Abinit comes with its own internal DMFT code and solvers, but they aren't numerically exact
and use the density-density approximation. TRIQS/CT-HYB, on the other hand, can handle the
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

A DFT+DMFT charge self-consistent calculation usually consists of two datasets: the first one
is a well-converged DFT calculation, and the second one a DFT+DMFT calculation, using the DFT
output density of the first dataset as a starting point.

We won't go over DFT parameters or how to converge them since that's already covered in the
DFT tutorials. Still, keep in mind that you may need to adjust your values and redo convergence
studies once DMFT is added. For example, you may need to converge your Kohn-Sham wavefunctions
better by increasing [[nline]] and [[nnsclo]]. Also, the electronic temperature (set via [[tsmear]])
becomes more important when DMFT is included, as temperature effects are directly taken into account
in the impurity solver, while they are often neglected in the DFT exchange-correlation functional.

Now let's go through the main DMFT variables for the second dataset. More details are in the variable
glossary.

### Activate DMFT

To turn on DMFT, just set:

    usedmft 1

in your input.

### Correlated electrons

Next, choose which atom types will be treated with DMFT and which angular momentum channel is correlated.
Currently, only one correlated angular momentum channel per atom type is allowed. This works like in
DFT+U using the variable [[lpawu]].

For each atom type, set [[lpawu]] to the correlated channel (e.g., 2 for $d$ orbitals, 3 for $f$ orbitals). Set [[lpawu]] = -1 if you don't want DMFT on that atom type.

Our interface lets you apply DMFT on multiple atoms, but each atom is treated separately, and interatomic
interactions are not handled dynamically. Besides, only [[lpawu]] $\le$ 3 is supported.

If you want DMFT on a subset of orbitals, use [[dmft_t2g]] (for $t_{2g}$ orbitals only) or [[dmft_x2my2d]] (for the $x^2-y^2$ orbital only).

### Correlated orbital

After choosing the correlated channel (angular part), you need to define the radial part (multiplied
by $r$) $u_l(r)$ of the local orbitals $\frac{u_l(r)}{r} \, Y_{lm}(\hat{r})$. This is set with
[[dmft_orbital]].

Pick an orbital that captures as many electrons as possible for DMFT, focusing on those near half-filling
around the Fermi level where correlations matter most.

Since DMFT usually applies to localized electrons, an atomic orbital is often best. By default,
[[dmft_orbital]] is the most bound atomic orbital in your PAW dataset, but you can choose any orbital
by providing one from a file.

For technical reasons, the orbital must then be projected onto a Hilbert space of finite dimension.
In our case, this is chosen as a finite set of Kohn-Sham bands (comprised between
[[dmftbandi]] and [[dmftbandf]]). A larger energy window gives more localized orbitals. An infinite
window recovers the original [[dmft_orbital]]. Make sure the window is large enough for DMFT to be
valid and to avoid overlap with neighboring atoms.

Projected orbitals are no longer orthonormal, so Wannier functions are built by orthonormalizing them
according to the scheme set by [[dmft_wanorthnorm]].

### Interaction tensor

DMFT electrons don't feel the bare Coulomb interaction - they interact via a screened potential, which
you must specify. We assume spherical symmetry (Slater parametrization) and thus only need the Slater
integrals $F_k$ to define the full interaction tensor.

Typically, we use instead the average charge interaction $U$ and Hund's exchange $J$ ([[upawu]] and
[[jpawu]] for each atom type), which are related to the Slater integrals.

The higher-order Slater integrals are set by default using the atomic values for $F_k / F_2$, but you
can override them with [[f4of2_sla]] and [[f6of2_sla]].

### Self-consistent cycle

A charge self-consistent calculation is comprised of two loops: the DFT+DMFT loop, whose number of
iterations is controlled by [[nstep]], and the DMFT inner loop, whose number of steps is controlled
by [[dmft_iter]]. Thus, the total number of calls to the impurity solver will be
[[nstep]] * [[dmft_iter]]. Each DMFT loop is performed at fixed density, and the density is updated
[[nstep]] times in total, every [[dmft_iter]] iterations. The mixing of the self-energy is controlled
by [[dmft_mxsf]].

### Magnetism

Next, you need to choose the type of magnetism. With [[nsppol]]=1 and [[nspden]]=1,
a paramagnetic solution is enforced by symmetrizing the two spin channels of the Green's function,
though the impurity problem is still solved with all spin channels, allowing the description of
local magnetic moments. With [[nsppol]]=2 and [[nspden]]=2, collinear magnetism is possible.

In some cases, it is preferential to have a non-magnetic exchange-correlation potential in DFT,
and let magnetism arise from DMFT only. You can do that with the variable [[usepawu]].

Finally, non-collinear calculations are handled by our interface ([[nspinor]]=2 and [[nspden]]=4),
but we only consider spin-orbit coupling (SOC) as a perturbation at the DFT level, and we do not
include it in the model. SOC leads to the appearance of spin off-diagonal components,
with possibly an imaginary part, which are fully handled by our interface, although some features
are lacking in this case, such as the measurement of the density matrix in the solver
([[dmft_triqs_measure_density_matrix]]), which is crucial to compute structural properties.
We emphasize this is not a limitation of our implementation, but rather of TRIQS/CT-HYB itself.

### Double counting

When combining DFT and DMFT arise the double counting issue. The local interactions already treated
to some extent at the DFT level need to be subtracted in order to be
replaced by their exact values given by DMFT. We provide the exact formula, which subtract them exactly,
though it is quite harder to use and will thus presented later in its
own section. We allow several commonly used and simpler approximate formulas, controlled by the variable
[[dmft_dc]].

### Impurity Solver

We still haven't told Abinit to use the interface with TRIQS/CT-HYB rather than the internal
implementation. This is done via [[dmft_solv]], which controls the impurity solver.
The two relevant values here for TRIQS/CT-HYB are [[dmft_solv]]=6 and [[dmft_solv]]=7.
The first one uses TRIQS/CT-HYB and keep only the density-density components in the interaction
tensor, drastically reducing the computation time (though in this case, you would be better off
using a segment solver), while [[dmft_solv]]=7 uses the full rotationally invariant Slater
Hamiltonian. If you want to use the internal Abinit implementation
with the internal solvers, please look the relevant tutorial.









