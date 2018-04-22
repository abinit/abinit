---
description: How to select the SCF algorithm
authors: XG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to select the SCF algorithm with the ABINIT package.

## Introduction

Self-Consistent Field calculations allow to determine the solution of the
Kohn-Sham equations, ending with converged "self-consistent" wavefunctions,
density, and Kohn-Sham potentials. Different algorithms can be chosen to
converge to the solution of this set of equations, governed by the input
variable [[iscf]] and [[wfoptalg]]. [[iscf]] focuses on the density/potential
self-consistency algorithms, while [[wfoptalg]] focuses on the determination
of the wavefunction through the solution of the Shrodinger equation with fixed
Kohn-Sham potential.

The algorithm selected by [[iscf]] are the iterative kind, among which Pulay
mixing is one of the most efficient. Also, an efficient preconditioner will
speed up the convergence. Among different choices, a generalized Kerker
preconditioner is implemented, see [[diemac]], [[diemix]] and [[dielng]].  
In order to perform a non-self-consistent calculations of wavefunctions and
corresponding eigenvalues in a fixed potential, as for representing a full
band structure, the loop over density/potentials self-consistency must be
disabled, for which [[iscf]]=-2 must be chosen.

Among the algorithms to find the wavefunctions, selected by [[wfoptalg]], the
conjugate-gradient and the LOBPCG ones are the favourite. Use the Chebyshev
filtering for massive parallel runs.

Inner electronic eigenvalues can be computed thanks to the minimisation of the
residual with respect to a target energy value, see [[eshift]].


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

