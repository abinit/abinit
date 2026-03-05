---
description: Linear-response calculations with constrained magnetic moments (Constrained DFPT).
authors: MR and MS
---
<!--- This is the source file for this topic. Can be edited. -->

This page describes how to perform noncollinear DFPT calculations in magnetic insulators 
using the constrained magnetic moments formalism implemented in ABINIT. This approach, referred to 
as **Constrained DFPT**, allows one to compute static and dynamic (at finite frequency) response functions 
--spanning interatomic force constants, dielectric tensor, Born charges, etc.. but also magnetic susceptibilities, 
magnetic Born charges and magnetoelectric tensors-- while enforcing parametric control over the local magnetic moments.

The implementation follows the theoretical framework introduced in [[cite:Royo2026]].

It is intended for medium to advanced users familiar with noncollinear magnetism and DFPT.

---

## Introduction

In systems where time-reversal (TR) symmetry is broken due to the presence of internal ordering of spins,
a standard, static linear-response calculation yields unphysical response functions that are invariant under TR symmetry. 
In the context of phonons, this was first noticed by Mead and Truhlar [[cite:Mead1979]] who, by carefully 
considering the phases of the nuclei and electrons wave functions when imposing the Born-Oppenheimer 
approximation, introduced a vector-potential term in the effective Hamiltonian of the nuclei. This new 
contribution enters the phonon equations of motion as a Berry curvature in the parameter space of the nuclei displacements, 
it has the physical meaning of a force induced by a velocity, and restores the expected magnetic symmetries of the crystal [[cite:Bonini2023]].

Magnetic materials can also host spin-wave excitations (magnons) which typically overlap in energy with phonons
and introduce additional complications in the linear-response regime. On the one hand, magnons and phonons can interact
mutually influencing each other's spectra and, therefore, need to be simultaneously treated. This was solved in 
[[cite:Ren2025]] by working with a set of Hessians and Berry curvatures defined in an extended parameter space of 
atomic displacements, local spin cantings and interactions thereof. The resulting generalized equation of motion
provides the eigenfrequencies and eigenvectors of the coupled magnon-phonon system. On the other hand, the so-called
acoustic magnons typically have very low frequencies at the center of the Brillouin zone, a fact that causes severe numerical 
instabilities in the self-consistent linear-response calculation whenever a given perturbation couples with these magnon excitations.

The constrained DFPT method implemented in ABINIT resolves these convergence issues by introducing a penalty functional that stiffens 
the magnetic degrees of freedom during the linear-response calculation. The magnetic moments are constrained to remain 
close to their ground-state configuration, thereby eliminating problematic low-energy resonances from the self-consistent loop.
The physically meaningful response functions, i.e., those without the constraints, are subsequently reconstructed in ANADDB via exact 
linear-algebra relations derived from Legendre transformations.

---
## Theoretical background

The theoretical formalism is described in detail in [[cite:Royo2026]]. Below, we summarize its key aspects. 

### Magnetic functionals and Legendre transforms

The formalism is based on four different magnetic energy functionals:

- Penalty-based internal energy \( \tilde U(B_l) \) , with \( B_l \) as target magnetic moments ( /(l=\kappa,\alpha /) is a composite index that runs over the magnetic sites and Cartesian directions),
- Modified enthalpy \( \tilde F(H_l) \), with \( H_l \) as local Zeeman fields.
- Internal energy based on holonomic constraints \( U(M_l) \), with \( M_l \) as target magnetic moments,
- Magnetic enthalpy  .

The two internal-energy functionals correspond to those introduced by Ma and Dudarev [[cite:Ma2015]] and by Gonze \emph{et al.} [[cite:Gonze2022]] and mentioned in [[topic:ConstrainedDFT]]. The corresponding enthalpies are in turn obtained via Legendre transforms. In [[cite:Royo2026]] it was demonstrated that the general response functions (second-order derivatives of the total energy with respect to two arbitrary perturbations) calculated using these four functionals are related via trivial linear-algebra equations. This means that all four functionals provide the same information, however, in practice it is more convenient to perform the linear-response calculation on the
internal energies, as they are free from the problematic low-energy magnons. The physical, spin-relaxed response functions on the magnetic enthalpy \( F(H_l) \) are subsequently 

In the current ABINIT implementation, the linear-response calculation is done with the penalty-based functional  \( \tilde U(B_l) \). The output of such a calculation is 
written in DDB files that are subsequently used by ANADDB to obtain the physical, spin-relaxed response functions defined on the magnetic enthalpy \( F(H_l) \).

---

### Dynamical regime

Regarding the dynamical linear-response regime required to study magnets, the implementation offers two routes.

---

#### First-order adiabatic approximation

The first one is the so-called first-order adiabatic approximation (FOA) in Ref. [[cite:Royo2026]]. It emerges from adopting an
adiabatic expansion on the frequency dependence of the second-order internal energies: 

\[
U_{\lambda_1,\lambda_2}(\omega) = K_{\lambda_1,\lambda_2} + i\omega G_{\lambda_1,\lambda_2} - \omega^2 M_{\lambda_1,\lambda_2} + \dots
\]

where, \( \lambda_1, \lambda_2 \) indicate two possible perturbations (atomic displacements, electric or Zeeman fields).
The FOA stops the above expansion at first order in the frequency and, therefore, requires the calculation 
of static Hessians \( {\bf U} \) plus Berry curvatures \( {\bf G} \). The Hessians are obtained via static 
constrained DFPT calculations. The Berry curvatures, in turn, correspond to frequency derivatives of second-order 
energies in the static \( \omega \rightarrow 0 \) limit whose calculation, in the most general case, boils down to 
obtain:

\[
\frac{d E_{ab}({\bf q},\omega)}{d\omega}=
\int [d^3 k] \sum_m f_{n\bf k} \left( \langle u_{m\bf k,-q}^{\lambda_2} | {u}_{m\bf k,-q}^{\lambda_1} \rangle - 
\langle {u}_{m\bf k,q}^{\lambda_1}|u_{m\bf k,q}^{\lambda_2} \rangle \right).
\]

The above equation --related with the calculation of a time-dispersion property-- has been implemented in the longwave
driver of ABINIT --originally devoted to calculate spatial-dispersion properties and now generalized to the case of time-- 
for any arbitrary pair of perturbations. The longwave driver reads the first-order wave functions pre-calculated in a 
constrained DFPT run and computes a set of constrained Berry curvatures. 

Both Hessians and Berry curvatures are calculated with the penalty-based functional \( \tilde U(B_l) \) and written in 
DDB files as second- and third-order total-energy derivatives, 
respectively. These DDB files are then used by ANADDB to switch between the different Legendre related magnetic functionals
[[cite:Royo2026]]. For instance, this allows one to obtain \( U_{\lambda_1,\lambda_2}(\omega) \) via a lineal (in this case) interpolation 
in frequency to subsequently convert it to the physical spin- and ion-relaxed enthalpies: frequency dependent susceptibilites
(dielectric, magnetic, magnetoelelectric, etc...) including coupled phonon and magnon resonances. 

---

#### Frequency-dependent DFPT





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
