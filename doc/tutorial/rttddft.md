---
authors: FBrieuc
---

# Tutorial on real-time time-dependent DFT (RT-TDDFT)

This tutorial aims at showing how to perform a RT-TDDFT calculation using Abinit.

!!! warning

    RT-TDDFT is under active development and still experimental. It should thus be used with caution!

If you have not yet done any of the basic tutorials it might be useful to do 
at least the first three basic ones: [basic1](/tutorial/base1), [basic2](/tutorial/base2) and
[basic3](/tutorial/base3).
It could also be useful that you know how to do PAW calculations using Abinit. 
If you are interested you can follow the first two tutorials on PAW: 
[PAW1](/tutorial/paw1) and [PAW2](/tutorial/paw2).

This tutorial should take about XXX hour to complete, depending on the amount of 
computing power you have access to and how much of the convergence study you decide to do.

[TUTORIAL_README]

## 1. Basics of RT-TDDFT

Time-dependent DFT (TDDFT) has been developed in the 1980s as a generalization of 
DFT to the time-dependent case and is based on similar theorems than the standard 
Hohenberg and Kohn theorems. They were first developed, in the general case, by 
Runge and Gross in 1984 [[cite:Runge1984]] establishing a one-to-one correspondence
between the time-dependent potential and the time-dependent density $n(t)$.
In the spirit of Kohn and Sham, Runge and Gross also introduced a non-interacting
system evolving in an effective time-dependent potential such that its density 
is the same as the interacting system of interest.
Further work by Van Leeuwen [[cite:VanLeeuwen1999]] proved the existence of this 
effective potential under some rather mild approximations. 

Within this framework, one ends up with the so-called time-dependent Kohn-Sham (TDKS)
equations that have a similar form as the time-dependent Shrödinger equation
\begin{equation}
    i\hat{S}\partial_t \psi_{nk}(\vec{r},t) = \hat{H} \psi_{nk}(\vec{r},t),
\end{equation}
but with the Hamiltonian being now a functional of the density and the initial 
state $\psi^0$, $\hat{H} \equiv \hat{H}[n,\psi^0]$.
Note that formally the Hamiltonian is a functional of the density at all previous times 
ie. it depends on the whole history. This is quite complicated to account for as
it requires a functional that is non-local both in time and space. 
In the following, we will thus always work in the so-called _adiabatic approximation_
which discards this dependence on the past history and consider that the Hamiltonian 
is a functional of the density at time $t$ only.
The operator $\hat{S}$ is the overlap operator that appears in PAW. 
For norm-conserving pseudopotentials this operator is just the unity.

In the real-time formulation of TDDFT (RT-TDDFT), the TDKS equations are numerically 
integrated in _real-time_ by discretizing the time with a time step $\Delta t$. 
The propagator to go from time $t$ to $t+\Delta t$ can be approximated as
\begin{equation}
    \hat{U}(t,t+\Delta t) = \exp\left(-i\hat{S}^{-1}\hat{H}\Delta t\right).
\end{equation}
In the present implementation only two simple propagators are available: 
the so-called exponential rule (ER) and exponential midpoint rule (EMR)
propagators. In both of them, the exponential in the propagator is approximated 
by a Taylor expansion (usually up to fourth order as it was found to be the best
choice to ensure good stability while minimizing computation time [[cite:Castro2004]]).
These two propagators differ only by the choice of the Hamiltonian used in the exponential.
The simplest ER propagator uses the Hamiltonian at time $t$, while the EMR propagator
uses the Hamiltonian at time $t+\Delta t/2$. This last choice ensure the correct 
time reversal symmetry and is thus supposed to be better than the simple ER propagator. 
However, this choice introduces a self-consistency since one needs the density at 
time $t+\Delta t/2$ to evolve the orbitals at time $t$. This is resolved by a
predictor-corrector scheme. See [[cite:Castro2004]] for more information 
on propagators for the TDKS equations including the ER and EMR propagators.

It is possible to apply an impulse external electric field in order to 
compute the associated response functions: the electronic conductivity 
and the dielectric function. This tutorial concerns the computation of 
the dielectric function of diamond from the response to an impulse electric
field computed using RT-TDDFT.

!!! note

    The application of an external electric field in RT-TDDFT is only available in PAW.
    We will thus work in this formalism throughout this tutorial.

## 2. Ground state of diamond
_Before beginning, you might consider working in a different subdirectory, as for the other tutorials.
Why not Work_rttddft? In what follows, the name of files are mentioned as if you were in such a
 subdirectory in $ABI_TESTS/tutorial.
All the input files can be found in the `$ABI_TESTS/tutorial/Input` directory._

!!! important

    You can compare your results with reference output files located in
    `$ABI_TESTS/tutorial/Refs` directories (for the present tutorial they 
     are named `trttddft_*.abo`).

The first thing that a real-time TDDFT (RT-TDDFT) calculation needs is a set of initial orbitals.
In the following we will start our calculations from the ground state KS orbitals. 
Thus let us start by computing the ground state orbitals of diamond. 
This system was actually treated in details in the first PAW tutorial and we will use the same description here. 
If you have not done this first tutorial on PAW [PAW1](/tutorial/paw1) it might be a good idea to do it first.
