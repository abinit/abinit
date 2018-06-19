---
description: How to perform calculation with constrained atomic magnetic moments
authors: EB
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform calculation with constrained atomic magnetic moments with the ABINIT package.

## Introduction

A complementary magnetic constraint method has been implemented in the ABINIT
code, wherein the magnetization around each atom is pushed to a desired
(vectorial) value. The constraint can either be on the full vector quantity,
$\vec{m}$, or only on the direction **m**. This is mainly useful for non
collinear systems, where the direction and amplitude of the magnetic moment
can change. The method follows that used in the Quantum Espresso
[[cite:Moscaconte2007]] and VASP [[cite:Ma2015]] codes: a Lagrangian
constraint is applied to the energy, and works through a resulting term in the
potential, which acts on the different spin components. The magnetization in a
sphere  Ωi around atom i at position **R** i is calculated as:

**m** = ∫Ωi **m** ( **r** ) d **r**

and the corresponding potential for spin component α is written as:

Vα = 2 λ f(| **r** - **R** i| / rs) **c** α

The function f(x) = x2(3+x(1+x(-6+3x))), is applied to smooth the transition
near the edge of the sphere around **R** i, over a thickness rs (by default
0.05 bohr, and f is set to 0 for | **r** - **R** i|> rs). This minimizes
discontinuous variations of the potential from iteration to iteration.

The constraint is managed by the keyword [[magconon]]. Value 1 gives a
constraint on the direction ( **c** = **m** \- **s** i ( **s** i. **m** ),
value 2 gives a full constraint on the vector ( **c** = **m** \- **s** i),
with respect to the keyword [[spinat]] ( **s** i above), giving a 3-vector for
each atom. The latter is quite a stringent constraint, and often may not
converge. The former value usually works, provided sufficient precision is
given for the calculation of the magnetic moment (kinetic energy cutoff in
particular).

The strength of the constraint is given by the keyword [[magcon_lambda]] (λ
above - real valued). Typical values are 10-2 but vary strongly with system
type: this value should be started small (here the constraint may not be
enforced fully) and increased. A too large value leads to oscillations of the
magnetization (the equivalent of charge sloshing) which do not converge. A
corresponding Lagrange penalty term is added to the total energy, and is
printed to the log file, along with the effective magnetic field being
applied. In an ideal case the energy penalty term should go to 0 (the
constraint is fully satisfied).



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

