---
description: Linear-response calculations with constrained magnetic moments (Constrained DFPT).
authors: MR and MS
---
<!--- This is the source file for this topic. Can be edited. -->

This page describes how to perform noncollinear DFPT calculations in magnetic insulators 
using the constrained magnetic moments formalism implemented in ABINIT. This approach, referred to 
as **Constrained DFPT**, allows one to compute static and dynamic (at finite frequency) response functions 
--spanning interatomic force constatns, dielectric tensor, Born charges, etc.. but also magnetic susceptibilities, 
magnetic Born charges and magnetoelectric tensors-- while enforcing parametric control over the local magnetic moments.

The implementation follows the theoretical framework introduced in [[cite:Royo2025]].

It is intended for medium to advanced users familiar with noncollinear magnetism and DFPT.

---

## Introduction

In systems where time-reversal (TR) symmetry is broken due to the presence of internal ordering of spins,
a standar, static linear-response calculation yields unphysical response functions that are invariant under TR symmetry. 
In the context of phonons, this was first noticed by Mead and Truhlar [[cite:Mead1979]] who, by carefully 
considering the phases of the nuclear and electron wave functions when imposing the Born-Oppenheimer 
approximation, introduced a vector-potential term in the effective Hamiltonian of the ions. This new 
term enters the Born-Oppenheimer equations of motion as a Berry curvature and restores the expected 
magnetic symmetries of  in the parameter space of the ionic displacements. 




Noncollinear magnetic systems pose severe convergence difficulties within standard density-functional 
perturbation theory (DFPT). The presence of low-energy magnon excitations leads to a poorly conditioned 
linear-response problem, especially:

- At the Brillouin zone center (acoustic magnon as Goldstone mode),
- In frequency-dependent calculations near magnon resonances,
- In coupled spin–phonon problems.

Constrained DFPT resolves these issues by introducing a penalty functional that stiffens the magnetic 
degrees of freedom during the linear-response calculation. The magnetic moments are constrained to remain 
close to their ground-state configuration, thereby eliminating problematic low-energy resonances from the 
self-consistent loop.

The physically meaningful (relaxed-spin) response functions are subsequently reconstructed via exact 
linear-algebra relations derived from Legendre transformations.

This strategy:

- Dramatically improves convergence,
- Preserves full formal equivalence with the unconstrained theory,
- Enables robust access to dynamical spin–phonon response functions.

---

## Theoretical background

### Magnetic functionals and Legendre transforms

The formalism is based on introducing a penalty-modified Kohn–Sham functional,

\[
\tilde U = E_{\mathrm{KS}} + \frac{\alpha}{2} \sum_l (B_l - m_l)^2,
\]

where:

- \( m_l \) are local magnetic moments,
- \( B_l \) are target magnetic moments,
- \( \alpha \) is a penalty parameter.

This defines a **constrained-B internal energy functional**.  

From it, one can construct a family of equivalent magnetic thermodynamic functionals:

- Internal energy \( U(M_l) \) (constrained magnetic moments),
- Magnetic enthalpy \( F(H_l) \) (fixed local Zeeman fields),
- Penalty-based internal energy \( \tilde U(B_l) \),
- Modified enthalpy \( \tilde F(H_l) \).

All response functions can be expressed as second derivatives of one of these functionals. 
The various representations are related via exact Legendre transformations.

In practice, ABINIT performs the linear-response calculation using the penalty-based functional 
\( \tilde U \), while the physically meaningful response tensors are reconstructed at post-processing 
level.

---

### Spin susceptibilities

The physical spin susceptibility matrix is defined as

\[
\chi_{jl} = \frac{\partial m_j}{\partial H_l}.
\]

Within the constrained formalism, it is reconstructed from the second derivatives of the penalty functional:

\[
\chi = \tilde U^{-1} - \frac{1}{\alpha} I.
\]

Thus, the unconstrained susceptibility is obtained by simple matrix inversion and subtraction of the 
penalty contribution.

---

### Spin–phonon coupling

The method generalizes naturally to the coupled spin–phonon problem.

The Hessian matrix in the extended parameter space takes block form:

\[
U =
\begin{pmatrix}
U^{(ss)} & U^{(sp)} \\
U^{(ps)} & U^{(pp)}
\end{pmatrix}
\]

where:

- \( ss \): spin–spin sector,
- \( pp \): phonon–phonon sector,
- \( sp \), \( ps \): mixed spin–phonon couplings.

The relaxed-spin interatomic force constants are obtained as

\[
\Phi^{RS} = \Phi^{FS} - U^{(sp)} \chi U^{(ps)}.
\]

This partition clearly separates:

- Frozen-spin (non-resonant) contributions,
- Resonant corrections mediated by spin canting.

---

### Dynamical regime

At finite frequency \( \omega \), the formalism is formulated within time-dependent DFPT 
(using the Kohn–Sham action functional).

A key advantage of the constrained approach is that the frozen-spin internal energy matrix 
\( U(\omega) \) is smooth and weakly frequency-dependent. This enables a controlled 
adiabatic expansion:

\[
U(\omega) = K + i\omega G - \omega^2 M + \dots
\]

Two levels of approximation are available:

- **First-order adiabatic approximation (FOA)**: neglect electronic mass corrections.
- **Second-order adiabatic approximation (SOA)**: includes renormalization of magnon and phonon masses.

The SOA yields modified equations of motion including finite magnon inertia.

---

## Practical implementation in ABINIT

The constrained DFPT implementation builds upon:

- Noncollinear magnetism,
- Local Zeeman perturbations,
- Dynamical DFPT (frequency-dependent Sternheimer equation).

The main features are:

1. Addition of a penalty term to the magnetic kernel,
2. Modified interaction kernel in the SCF cycle,
3. Post-processing reconstruction of relaxed-spin response functions.

No modification of the potential-mixing scheme is required.

---

## General workflow

A typical constrained DFPT calculation proceeds as follows:

### 1. Ground-state calculation

- Perform a noncollinear magnetic ground-state calculation.
- Ensure convergence of local magnetic moments.
- Define Wigner–Seitz spheres for magnetic moment integration.

---

### 2. Constrained linear-response calculation

- Activate constrained magnetic moment mode.
- Choose an appropriate penalty parameter `alpha`.
- Perform DFPT calculations for desired perturbations:
  - Local Zeeman fields,
  - Atomic displacements,
  - Electric fields (if required),
  - Finite frequency (if dynamical response is desired).

The response is computed in the frozen-spin representation.

---

### 3. Reconstruction of physical response functions

After obtaining second derivatives of the constrained functional:

- Invert spin–spin block,
- Apply Legendre relations,
- Reconstruct:
  - Relaxed-spin force constants,
  - Spin susceptibilities,
  - Dielectric tensor,
  - Magnetoelectric response,
  - Spin–phonon Green's functions.

---

## Choice of penalty parameter

The parameter `alpha` must be chosen carefully:

- Too small: magnon resonances remain in the low-energy window.
- Too large: numerical instabilities may appear.

Guidelines:

- Increase `alpha` until spin–spin correlation poles move above the optical gap.
- Verify stability of reconstructed relaxed-spin quantities with respect to `alpha`.

Final reported physical results should be independent of `alpha`.

---

## Limitations

Current implementation is restricted to:

- Insulating systems,
- Noncollinear magnetism,
- Adiabatic local/semi-local exchange-correlation functionals,
- Time-reversal broken magnetic states.

Metals are not supported.

---

## Related Input Variables

{{ related_variables }}

---

## Selected Input Files

{{ selected_input_files }}

---

## Tutorials

A dedicated tutorial on coupled spin–phonon response and electromagnons 
is in preparation.
