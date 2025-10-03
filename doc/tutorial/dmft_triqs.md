---
author: ECastiel
---

# Tutorial on DFT+DMFT with TRIQS/CT-HYB

## Introduction

In this tutorial, you'll learn how to run a DFT+DMFT calculation using the interface
between Abinit and TRIQS/CT-HYB.
TRIQS/CT-HYB is used as the impurity solver,
while Abinit handles the rest of the calculation and connects to TRIQS/CT-HYB
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
[[https://docs.abinit.org/variables/dmft/|DMFT]] and
[[https://docs.abinit.org/variables/dmft_triqs/|DMFT_TRIQS]] sections of the variable glossary.
Most TRIQS/CT-HYB specific variables directly correspond to TRIQS/CT-HYB's own inputs, so if you're
already familiar with it, you shouldn't get lost. Finally, we encourage to read the
[[topic:dmft|topic on DMFT]].

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











