---
description: How to perform calculation with constrained atomic magnetic moments
authors: EB
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform calculation with constrained atomic magnetic moments with the ABINIT package.

## Introduction

A complementary magnetic constraint method has been implemented in the ABINIT
code, wherein the magnetization around each atom $I$ is pushed to a desired
(vectorial) value. The constraint can either be on the full vector quantity,
$\vec{M}_I$, or only on the direction $\vec{e}_I$. This is mainly useful for non
collinear systems, where the direction and amplitude of the magnetic moment
can change. The method follows that used in
<!--- Quantum Espresso [[cite:Moscaconte2007]] and -->
VASP [[cite:Ma2015]]: a Lagrangian
constraint is applied to the energy, and works through a resulting term in the
potential, which acts on the different spin components. The magnetization in a
sphere  $\Omega_I$ around atom $I$ at position $\vec{R}_I$ is calculated as:

$$ \vec{M}_I = \int_{\Omega_I} \vec{m}(\vec{r}) F_I(|\vec{r}-\vec{R_I}|) d\vec{r} $$

and the corresponding potential for spin component $\alpha$  is written as:

$$V_{\alpha}(\vec{r}) = -2 \lambda F_I(|\vec{r}-\vec{R_I}|) \vec{c}  \vec{\sigma}_{\alpha}. $$

The function $F_I$ is zero outside of a sphere of radius [[ratsph]] for atom $I$,
and inside the sphere,  the function $f(x) = x^2(3+x(1+x(-6+3x)))$ (see Eq. (B4) in [[cite:Laflamme2016]]),
is applied to smooth the transition in a region of thickness $r_s$ (fixed to 0.05 bohr), otherwise it is 1.
This minimizes discontinuous variations of the potential from iteration to iteration.

The constraint is managed by the keyword [[magconon]]. Value 1 gives a
constraint on the direction:

$$\vec{c}  = (\vec{M}_I/|\vec{M}_I|)-(\vec{M}^{spinat}_I / |\vec{M}^{spinat}_I|),$$

while value 2 gives a full constraint on the vector 

$$\vec{c}  = \vec{M}_I - \vec{M}^{spinat}_I,$$

where in both cases $\vec{M}^{spinat}_I$ defined by [[spinat]], giving a 3-vector magnetic potential for
each atom. The latter is quite a stringent constraint, and often may not
converge. The former value usually works, provided sufficient precision is
given for the calculation of the magnetic moment (kinetic energy cutoff in
particular).

The strength of the constraint is given by the keyword [[magcon_lambda]] ($\lambda$
above - real valued). Typical values are $10^{-2}$ but vary strongly with system
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

