---
description: Linear-response calculations with constrained magnetic moments (Constrained DFPT).
authors: MR and MS
---

<!--- This is the source file for this topic. Can be edited. -->

This page describes how to perform noncollinear DFPT calculations in magnetic insulators 
using the constrained magnetic moments formalism implemented in ABINIT. This approach, referred to 
as **Constrained DFPT**, allows one to compute static and dynamic (finite-frequency) response functions 
--including interatomic force constants, dielectric tensors, Born effective charges, magnetic susceptibilities, 
magnetic Born charges, and magnetoelectric tensors-- while enforcing parametric control over the local magnetic moments.

The implementation follows the theoretical framework introduced in [[cite:Royo2026]].

It is intended for medium to advanced users familiar with noncollinear magnetism and DFPT.

---

# Introduction

In systems where time-reversal (TR) symmetry is broken due to internal spin ordering,
a standard static linear-response calculation yields unphysical response functions that are invariant under TR symmetry.

In the context of phonons, this issue was first noticed by Mead and Truhlar [[cite:Mead1979]], who, by carefully 
considering the phases of the nuclear and electronic wavefunctions within the Born–Oppenheimer 
approximation, introduced a vector-potential term in the effective Hamiltonian of the nuclei. 
This additional contribution enters the phonon equations of motion as a Berry curvature in the 
parameter space of nuclear displacements. Physically, it represents a velocity-induced force and 
restores the expected magnetic symmetries of the crystal [[cite:Bonini2023]].

Magnetic materials can also host spin-wave excitations (magnons), which typically overlap in energy with phonons 
and introduce additional complications in the linear-response regime. On the one hand, magnons and phonons can interact,
mutually influencing each other's spectra, and therefore must be treated simultaneously. This problem was addressed in 
[[cite:Ren2024]] by working with a set of Hessians and Berry curvatures defined in an extended parameter space of 
atomic displacements, local spin cantings, and their mutual interactions. The resulting generalized equations of motion 
provide the eigenfrequencies and eigenvectors of the coupled magnon–phonon system.

On the other hand, acoustic magnons typically have very low frequencies at the center of the Brillouin zone. 
This causes severe numerical instabilities in the self-consistent linear-response calculation whenever a perturbation 
couples to these magnon excitations.

The constrained DFPT method implemented in ABINIT resolves these convergence issues by introducing a penalty functional 
that stiffens the magnetic degrees of freedom during the linear-response calculation. The magnetic moments are constrained 
to remain close to their ground-state configuration, thereby eliminating problematic low-energy resonances from the 
self-consistent loop.

The physically meaningful response functions, i.e. those without the constraints, are subsequently reconstructed in ANADDB 
via exact linear-algebra relations derived from Legendre transformations.

The Constrained DFPT implementation is still under active development. Users are strongly encouraged to contact the developers (Miquel Royo and Massimiliano Stengel) before using it for production calculations.

---

# Theoretical background

The theoretical formalism is described in detail in [[cite:Royo2026]]. Below we summarize its key aspects.

## Magnetic functionals and Legendre transforms

The formalism is based on four different magnetic energy functionals:

- Penalty-based internal energy \( \tilde{U}(B_l) \), with \( B_l \) the target magnetic moments  
  (where \( l = (\kappa,\alpha) \) is a composite index running over magnetic sites and Cartesian directions).

- Modified enthalpy \( \tilde{F}(H_l) \), with \( H_l \) the local Zeeman fields.

- Internal energy based on holonomic constraints \( U(M_l) \), with \( M_l \) the target magnetic moments.

- Magnetic enthalpy \( F(H_l) \).

The two internal-energy functionals correspond to those introduced by Ma and Dudarev [[cite:Ma2015]] and by Gonze *et al.* [[cite:Gonze2022]], 
and are discussed in [[topic:ConstrainedDFT]]. The corresponding enthalpies are obtained via Legendre transformations.

In [[cite:Royo2026]] it was demonstrated that the general response functions (second derivatives of the total energy with respect to two arbitrary perturbations) calculated using these four functionals are related by simple linear-algebra relations. Therefore, all four functionals contain the same physical information. In practice, however, it is more convenient to perform the linear-response calculation using the internal energies, as they are free from the problematic low-energy magnon resonances.

In the current ABINIT implementation of constrained DFPT, the linear-response calculation is performed with the penalty-based functional \( \tilde{U}(B_l) \). The output of this calculation is written in DDB files that are subsequently used by ANADDB to extract the physical, spin-relaxed response functions defined in terms of the magnetic enthalpy \( F(H_l) \).

---

## Dynamical regime

Two different approaches are available to treat dynamical linear-response effects in magnetic systems.

---

### First-order adiabatic approximation

The first approach is the **first-order adiabatic approximation (FOA)** described in [[cite:Royo2026]]. 
It arises from an adiabatic expansion of the frequency dependence of the second-order internal energies:

\[
U_{\lambda_1,\lambda_2}(\omega)
=
K_{\lambda_1,\lambda_2}
+ i\omega G_{\lambda_1,\lambda_2}
- \omega^2 M_{\lambda_1,\lambda_2}
+ \dots
\label{eq:adiab}
\]

where \( \lambda_1 \) and \( \lambda_2 \) denote two possible perturbations (atomic displacements, electric or Zeeman fields, etc...).

The FOA truncates the expansion at first order in frequency. Consequently, it requires the calculation of:

- static Hessians \( \mathbf{U} \)
- Berry curvatures \( \mathbf{G} \)

The Hessians are obtained via static constrained DFPT calculations.  
The Berry curvatures correspond to frequency derivatives of the second-order energies in the static limit \( \omega \to 0 \). In the most general case, their calculation reduces to

\[
\frac{dE_{ab}(\mathbf{q},\omega)}{d\omega}
=
\int [d^3 k] \sum_n f_{n\mathbf{k}}
\left(
\langle u_{n\mathbf{k},-q}^{\lambda_2} | u_{n\mathbf{k},-q}^{\lambda_1} \rangle
-
\langle u_{n\mathbf{k},q}^{\lambda_1} | u_{n\mathbf{k},q}^{\lambda_2} \rangle
\right).
\]

This expression, related to the calculation of a time-dispersion property, has been implemented in the **longwave driver** of ABINIT. This driver, originally designed to compute spatial-dispersion properties, has been generalized to treat time dispersion as well. It reads the first-order wavefunctions obtained in a constrained DFPT calculation and computes a set of constrained Berry curvatures.

Both Hessians and Berry curvatures are calculated using the penalty-based functional \( \tilde{U}(B_l) \) and written in DDB files as second- and third-order total-energy derivatives, respectively. These files are then used by ANADDB to convert between the different Legendre-related magnetic functionals [[cite:Royo2026]].

For example, ANADDB can reconstruct \( U_{\lambda_1,\lambda_2}(\omega) \) through a linear interpolation in frequency using Eq. \ref{eq:adiab} truncated at first order. One can also enforce causality by replacing \( \omega \rightarrow \omega + i\eta \), where \( \eta \) is a positive infinitesimal.

The interpolated internal energies are then converted into the physical spin- and ion-relaxed enthalpies, yielding frequency-dependent susceptibilities (dielectric, magnetic, magnetoelectric, etc...), including coupled phonon and magnon resonances.

---

### Frequency-dependent DFPT

The DFPT linear-response driver in ABINIT was originally developed to treat static perturbations [[cite:Gonze1997]] [[cite:Gonze1997a]] within the adiabatic approximation. Recently, it has been extended to handle time-dependent perturbations oscillating at a given frequency \( \omega \).

Such perturbations intrinsically break time-reversal symmetry --even in nonmagnetic materials-- and therefore lift Kramers degeneracy. In practice, this means that the first-order wavefunctions at \( (\omega,\mathbf{q}) \) and \( (-\omega,-\mathbf{q}) \) are no longer related by simple complex conjugation. Instead, they must be computed independently by solving two coupled Sternheimer equations (or equivalently minimizing two second-order energy functionals, as implemented in ABINIT). These equations are coupled through the self-consistent fields (see Section III.A of [[cite:Royo2026]] for details).

The current conjugate-gradient minimization algorithm implemented in ABINIT enables frequency-dependent DFPT calculations in frequency ranges where the perturbations do not couple to resonances of the system. This essentially corresponds to transparent regimes in insulators. For magnetic insulators, however, frequency-dependent DFPT calculations must be combined with the constrained functional \( \tilde{U}(B_l) \) to avoid numerical instabilities caused by nearby magnon resonances.

Because the frequency dependence of the second-order internal energies \( U_{\lambda_1,\lambda_2}(\omega) \) is smooth and free of resonances in the sub-gap regime, it is not necessary to perform DFPT calculations on a dense frequency grid. Instead, a small number of calculations --typically fewer than five frequency points-- across the frequency range of interest is sufficient. The resulting DDB files are then merged with **MRGDDB** into a single DDB file containing different blocks corresponding to the computed \( (\omega,\mathbf{q}) \) cases. ANADDB then performs a fast polynomial interpolation in frequency of \( U_{\lambda_1,\lambda_2}(\omega) \), optionally including causality, and converts the results into spin- and ion-relaxed enthalpies.

---

# Limitations

The current implementation is restricted to:

- Magnetic insulating systems (metals are not supported)
- Local exchange-correlation functionals
- Norm-conserving pseudopotentials
- Atomic displacements, electric fields, Zeeman fields (uniform and local), and scalar potential perturbations

Some contributions are not yet included in the implementation (corresponding to the second term of Eq. (13) and Eq. (69) of [[cite:Royo2026]]). These terms may become important when the atomic spheres depend on the perturbation --for example when an atomic displacement moves the integration sphere along the atomic trajectory. These contributions vanish by symmetry in the systems studied in [[cite:Royo2026]], but they may become finite in other magnetic systems, such as those with noncollinear magnetic ground states.

The current constrained DFPT implementation cannot yet be used on top of a constrained ground-state calculation corresponding to any of the flavors described in [[topic:ConstrainedDFT]]. In practice, the linear-response calculation must start from a ground state with a stable magnetic configuration.

Finally, as mentioned above, ANADDB allows interpolation of the calculated response functions in frequency. However, combined interpolation in both frequency and momentum --needed, for instance, to obtain the Brillouin-zone dispersion of nonadiabatic coupled spin-phonon susceptibilities-- has not yet been implemented. Development in this direction is currently ongoing. 


---

# Related Input Variables

{{ related_variables }}

---

# Selected Input Files

{{ selected_input_files }}

---

# Tutorials

A dedicated tutorial on constrained DFPT calculations is currently in preparation.
