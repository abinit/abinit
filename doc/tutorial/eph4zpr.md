---
authors: MG
---

# Zero-point renormalization of the band gaps and temperature-dependent band structures

This tutorial explains how to compute the electron self-energy due to phonons, evaluate the zero-point
renormalization (ZPR) of the band gap, obtain spectral functions and temperature-dependent band structures.
We start with a very brief overview of the many-body formalism in the context of the e-ph interaction, 
then we discuss how to compute the ZPR, temperature-dependent band gaps and spectral functions.

## Formalism

The electron-phonon (e-ph) self-energy, $\Sigma^{\text{e-ph}}$, describes the renormalization of the single-particle
energy as well as the lifetime of the charged excitation due to the interaction with phonons.
This term should be added to the electron-electron (e-e) self-energy, $\Sigma^{\text{e-e}}$, that describes
many-body effects induced by the Coulomb interaction.
The e-e contribution can be estimated using, for instance, the $GW$ approximation but
in this tutorial we mainly focus on the e-ph self-energy ($\Sigma$ for brevity) and its temperature dependence.
Indeed, in semiconductors and insulators most of temperature dependence at low T originates from the e-ph interaction 
(and the thermal expansion).
Corrections due to $\Sigma^{\text{e-e}}$ are obviously important but their temperature dependence is rather small 
provided the system has an energy gap much larger than $kT$.

In state-of-the-art \abinitio methods, the e-ph coupling is described
within DFT by expanding the KS effective potential in the nuclear displacements,
and the vibrational properties are obtained with  (DFPT) [[cite:Gonze1997]], [[cite:Baroni2001]].
The e-ph self-energy consists of two terms: the Fan-Migdal (FM) self-energy 
and the Debye-Waller (DW) part [[cite:Giustino2017]]:

\begin{equation}
\Sigma^\eph = \Sigma^\FM(\ww) + \Sigma_^{\DW}
\end{equation}

The diagonal matrix elements of the (retarded) FM part in the KS basis set are given by

\begin{equation}
\begin{split}
    \Sigma^\FM_{n\kk}(\omega,\ef,T) =
                & \sum_{m,\nu} \int_\BZ \frac{d\qq}{\Omega_\BZ} |\gkq|^2 \\
                & \times \left[
                    \frac{n_\qnu(T) + f_{m\kk+\qq}(\ef,T)}
                         {\omega - \emkq  + \wqnu + i \eta} \right.\\
                & \left. +
                    \frac{n_\qnu(T) + 1 - f_{m\kk+\qq}(\ef,T)}
                         {\omega - \emkq  - \wqnu + i \eta} \right] ,
\end{split}
\label{eq:fan_selfen}
\end{equation}

where $f_{m\kk+\qq}(\ef,T)$ and $n_\qnu(T)$ correspond to the Fermi-Dirac and Bose-Einstein occupation functions
with $T$ the temperature and $\ef$ the Fermi energy.
<!--
For the sake of simplicity, the temperature and Fermi level are considered as parameters, and the dependence on $T$ and $\ef$ will be omitted in the following.
-->
The integral is performed over the $\qq$-points in the BZ of volume $\Omega_\BZ$ and $\eta$
is a positive real infinitesimal.
From a mathematical point of view, one should take the limit $\eta \rightarrow 0^+$.
At the level of the implementation, $\eta$ is set to a (small) finite value given by the [[zcut]] variable.

!!! tips

    In practice, one can use symmetries to reduce the integration to an appropriate irreducible zone defined by the
    operations of the little group of the $\kk$-point.
    This symmetrization trick is activated by default and can be deactivated by setting [[symsigma]] to zero.

<!--
First-principles calculations of the EPH self-energy are therefore crucial to understand the temperature-dependence
of band gaps, including the correction due to zero-point motion, as well as for computing phonon-limited mobilities
within the Boltzmann transport equation.
In ABINIT v9, it will be possible to compute the EPH self-energy in the Kohn-Sham representation using the EPH matrix elements.
The code employs optimized algorithms to compute either the full self-energy
(needed for QP corrections and spectral functions)
When computing the full self-energy, it is possible to reduce the number of empty states required for convergence
by using the first-order wavefunctions obtained by solving the relevant Sternheimer equation.
-->

The **static** DW term involves the second order derivative of
the KS potential with respect to the nuclear displacements.
State-of-the-art implementations approximate the DW contribution with

\begin{equation}
\label{eq:dw_selfen}
\Sigma_{n\kk}^{\DW} = \sum_{\qq\nu m} (2 n_{\qq\nu} + 1) \dfrac{g_{mn\nu}^{2,DW}(\kk, \qq)}{\ee_{n\kk} - \ee_{m\kk}},
\end{equation}

where $g_{mn\nu}^{2,\DW}(\kk,\qq)$ is an effective matrix element that, within the rigid-ion approximation,
can be expressed in terms of the $\gkq$ matrix elements [[cite:Giustino2017]].

Both the FM and the DW term converge slowly with respect to the number of bands [[nband]] and the $\qq$-sampling.
At the level of the implementation, the number of bands in the FM/DW expression is governed by [[nband]].
Note that, in principle, one should sum an infinite number of states (read from the external WFK file) 
and that the convergence is rather slow.
QP energy differences usually converge faster than absolute QP energies (the same behaviour is observed in $GW$).
The EPH code implements numerical tricks to accelerate the convergence with respect to [[nband]] that are discussed 
in more detail in the next sections.
Also the convergence with the $\qq$-sampling is rather slow and delicate.
This is especially true in polar materials due the divergence of the e-ph matrix elements 
in the $\qq \rightarrow 0$ limit [[cite:Vogl1976]].

### QP equation and spectral function

In principle, the quasi-particle (QP) excitations are defined by the solution(s)
in the complex plane of the equation $$ z = \ee_\nk + \Sigma_\nk^{\text{e-ph}}(z) $$.
provided that the non-diagonal components of the self-energy can be neglected
In practice, the problem is usually simplified by seeking approximated solutions along the real axis 
following two different approaches.

In the on-the-mass-shell approximation, the QP energy is given by the real part
of the self-energy evaluated at the bare KS eigenvalue

\begin{equation}
\ee^\QP_\nk = \ee_\nk + \Re\, \Sigma_\nk^{\text{e-ph}}(\ee_\nk).
\end{equation}

The second approach linearizes the self-energy near the KS eigenvalue and evaluates the QP correction using

\begin{equation}
  \ee^\QP_\nk = \ee_\nk
      + Z_\nk\,\Re  \Sigma^\text{e-ph}_\nk(\ee_\nk)
\end{equation}

with the renormalization factor $Z$ given by

\begin{equation}
  Z_\nk= \left(1 - \Re\left[ \frac{\partial\Sigma^\text{e-ph}_{\nk}}{\partial\ee}\right]\Bigg|_{\ee=\ee_\nk} \right)^{-1}.
\end{equation}

Both approaches are implemented in ABINIT although it should be noted that, according to recent works, 
the on-the-mass-shell approach provides results that are closer to those obtained with more advanced techniques
based on the cumulant expansion [[cite:Nery2018]]

<!--
Further details concerning the implementation are given in Ref.[[cite:Gonze2019]].
In order to accelerate the convergence with the number of empty states, ABINIT replaces the contributions
given by the high-energy states above a certain band index $M$ with the solution
of a non-self-consistent Sternheimer equation in which only the first $M$ states are required.
The methodology, proposed in [[cite:Gonze2011]], is based on a quasi-static approximation
in which the phonon frequencies in the denominator of Eq.~(\ref{eq:fan_selfen}) are neglected and
the frequency dependence of $\Sigma$ is approximated with the value computed at $\omega = \enk$.
This approximation is justified when the bands above $M$ are sufficiently high in energy with respect
to the $n\kk$ states that must be corrected.
The code can compute the spectral functions
-->

## Typical workflow for ZPR and T-dep calculations 

A typical workflow requires the following steps:

1. GS calculation to obtain the WFK and the DEN file.
   Remember to set [[prtpot]] to 1 to write the file with the KS potential required for the Sternheimer trick.

2. DFPT calculations for all the IBZ $\qq$-points corresponding to the *ab-initio* [[ddb_ngqpt]]
   used to perform the Fourier interpolation of the dynamical matrix and of the DFPT potentials.
   Remember to compute BECS and Q*

3. NSCF computation of a WFK file on a dense $\kk$-mesh containing the wavevectors 
   for which phonon-induced QP corrections are wanted. Use the DEN file from step #1. 
   Remember to include enough empty states so that it is possible to perform convergence studies wrt [[nband]]

4. Merge the partial DDB and POT files with *mrgddb* and *mrgdvdb*, respectively

5. Use the full DDB file, the DVDB file and the WFK file obtained in step 3 to start a ZPR calculation.

## Calculation of the ZPR for XXX

In this tutorial, we prefer to focus on the EPH part so we skip completely the DFPT computation 
and provide pre-computed DDB and POT1 files obtained with a xxx $\qq$-mesh.
We also provide a DEN.nc file that can be used to perform a NSCF calculation including empty states.

First of all, let's merge the DDB files using the following input file:

Now we can merge the DFPT potential with *mrgdv* and the input file:

For our first NSCF calculation, we use a xxx $\kk$-mesh and XXX bands.
The number of bands is sufficiently large so that we can perform initial convergence studies.
We also perform a NSCF calculation on a high-symmetry $\kk$-path to locate the position of the KS band edges
as these are the states we want to correct.
The input file is ...
We use [[getden_path]] to read the DEN.nc file instead of [[getden]] or [[irdden]] variables.

TODO: Mention ncdump and the input_string

In this particular case the CBM and the VMB ...

!!! important

  Note that the EPH code can compute QP corrections only for $\kk$-points belonging to the $\kk$-mesh
  associated to WFK file.
  In some cases, it is not possible to generate a $\kk$-mesh containing the band edges.
  In this case, you will have to select the closest wavevectors in the grid.

### Our first ZPR calculation

TODO: use and explaine [[structure]]

For our first ZPR we decided to use a rather simple input file that allows us to 
discuss the basic input variables.

[[optdriver]] [[eph_task]] is set to 4
[[dvdb_ngqpt]] is set to XXX
We use [[getddb_filepath]], [[getwfk_filepath]], [[getdvdb_filepath]] to specify the paths to the external files.
[[ngkpt]] 

[[nkptgw]]
[[kptgw]]

### Convergence wrt nband

At this point it should be not so difficult to write an input file to perform ZPR calculations for different
values of [[nband]].
An example of input file is given in 

grep something gives

So XXX bands are needed to convergece the ZPR

### How to reduce the number of bands with the Sternheimer approach

In this section, we discuss how to take advantange of the Sternheimer equation 
to accelerate the convergence with [[nband]].
To activate the Sternheimer method, use [[eph_stern]] = 1 
and pass the file to the GS KS potential file via [[getpot_filepath]].
The non-self-consistent Sternheimer equation is solved on-the-fly
using max [[nline]] iterations. 
The NSCF iteration stops when the solution is convergenced within [[tolwfr]].
In this example, we use the default values. You might need to adjust the values for more complicaed systems.

Note that the [[nband]]

### How to compute the spectral function

### Convergence of the ZPR wrt the q-mesh

### Estimating the ZPR with a generalized Fr\"ohlich model 

Last but not least, one can estimate the correction to the zero-point renormalization
of the band gap in polar materials using a generalized Fr\"ohlich model based on
*ab initio** effective masses computed with DFPT [[cite:Laflamme2016]]
For formalism is detailed in XXX.
An example input file is available in [[test:v7_88]].
