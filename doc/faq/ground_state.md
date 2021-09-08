---
authors: MG
---

# Ground-state FAQs

This page collects FAQs related to the GS part of Abinit.

## What are the dimensions defining the workload of a standard GS run?

In a GS run, we have to find the set of KS orbitals $\psi_\nks$ that solve the self-consistent problem:

\begin{align*}
H_\KS[n] \psi_\nks &= \epsilon_\nks \psi_\nks \\
n(\rr) &= \sum_\nks f_\nks \psi_\nks(\rr)
\end{align*}

where the number of bands $n$ is given by [[nband]], the number of **k**-points (in the IBZ) is [[nkpt]]
and the number of spin components $\sigma$ is given by [[nsppol]]  (note that we are assuming scalar wavefunctions,
that are either spin up or spin down. 
Two-component spinors require [[nspinor]] == 2).
The periodic part of the Bloch state is expanded in a planewave basis set according to:

$$
u_\nks(\rr) = \sum_{\GG}^{\text{npw}_\kk} u_\nks(\GG) e^{i\GG.\rr}
$$

For a given $\kk$ point only the $\text{npw}_\kk$ $\GG$-vectors with kinetic energy
$|\kG|^2/2 <= \text{ecut}_\text{eff}$ are included in the Fourier series
where $\text{ecut}_\text{eff} =$ [[ecut]] * [[dilatmx]] ** 2.
The maximum number of $\GG$-vectors over all the **k**-points is called [[mpw]].
Time-reversal symmetry is automatically exploited to halve the number of planewaves for
particular **k**-points, see the [[istwfk]] input variable.

The GS code stores in memory the Fourier components $u_\nks(\GG)$ and Fast Fourier Transforms
are used to go from $u_\nks(\GG)$ to $u_\nks(\rr)$.
The FFT mesh is defined by the first three entries of the [[ngfft]] array that are automatically
computed from [[ecut]] and [[boxcutmin]].
There are also other additional arrays that depend on [[natom]], especially if PAW is used.
All these variables are reported at the beginning of the main output file.

Keep in mind that there is a strong interplay between these dimensions and the MPI parallelization.
The parallelism over **k**-points and (collinear) spins is rather efficient
whereas the parallelism over bands and $\GG$-vectors is more network intensive.

!!! warning

    The conjugate gradient solver (default algorithm for GS calculations)
    **cannot use** more than [[nkpt]] * [[nsppol]] MPI processes.
    To take advantage of the band + $\GG$ parallelim use [[paral_kgb]] == 1 in the input file,
    possibly with [[autoparal]] = 1

## Is it a good idea to set the number of bands to the highest occupied state to accelerate a GS calculation?

Usually, the answer is **NO** and the explanation requires some technical discussion about the KS solver.
In principle, the workload as well as the memory footprint are proportional to [[nband]] so, at first sight,
reducing the number of KS states seems to be a smart choice.
On the other hand, a larger number of bands allows the eigenvalue solver to operate on a larger Hilbert subspace
and this can make a huge difference at the level of the number of iterations required to convergence.
This is especially true for algorithms such as the RMM-DIIS method ([[rmm_diis]] = 1).

In other words, more bands means more wall-time per SCF iteration but, on the other hand,
the SCF cycle is expected to reach converge in less iterations.

!!! important

    In the case of metals, one needs partially occupied states to account for the tail of the occupation
    function [[occopt]] with the broadening parameter [[tsmear]].
    Decreasing the number of bands is therefore a **very bad idea** and ABINIT will stop if there are not enough
    partially occupied states.

## Which convergence criterion do you recommend for the SCF cycle?

Well, the answer depends on the physical properties you are interested in.

Using [[toldfe]] to stop the SCF cycle is an excellent choice if you are interested in total energies,
but the other convergence criteria are much better if you are interested in other physical properties.

[[tolvrs]] is an excellent choice when computing the DEN file needed for a NSCF band structure calculation
since we want to make sure that our KS potential is well converged before using it in the NSCF run.

For structural relaxations or MD runs, one can use [[toldff]] or [[tolrff]] to obtain accurate forces
so that atoms are moved in the right direction at each relaxation step.
Note, however, that in high-symmetry systems the forces may be zero by symmetry.
To avoid this problem, use [[tolvrs]]

[[tolwfr]] is the most stringent convergence criterion since it takes into account the convergence of
all the [[nband]] - [[nbdbuf]] states (including **empty** states, if any).
This is the recommended approach when computing GS wavefunctions to be used in the DFPT code
or in the many-body part.
Keep in mind, however, that the last states usually require many more iterations to converge than the
low-energy states so we strongly recommend using [[tolwfr]] in conjunction with a [[nbdbuf]].

In the case of NSCF run, [[tolwfr]] is the only available convergence criterion

## The SCF cycle does not converge

The default SCF algorithm for SCF (see [[iscf]]) is an excellent compromise between speed and robustness.
It is referred to as Pulay II in Fig. 10 of
N. D. Woods and M. C. Payne and P. J. Hasnip, J. Phys.: Condens. Matter 31, 453001 (2019),
and appears on the Pareto frontier.
Some parameters in the Kerker preconditioning have also an effect on the balance between speed and robustness,
and the default values are good for metallic systems.
As thoroughly discussed in this reference, there is no fool-proof fast algorithm.

If your SCF cycle do not converge (see nres2 or vres2 column in the output file), the reasons can be:
(1) Insufficient underlying accuracy in the solution of the Schrödinger at fixed potential
(2) Transient non-linear behaviour of the SCF, due either to
    (2a) sudden change of occupation numbers (usually only for metals), or
    (2b) long-wavelength fluctuations of potential, bigger than the gap.

Start first to address (2), by some tuning which can come without significantly slowing ABINIT.
Test different values of [[diemac]].
The default is very large, and suitable fo metals.
Try running with a small value (perhaps 5) if your system is a bulk semiconductor or a molecule.
If it is a doped solid, some intermediate value (perhaps 50) might be appropriate.

If this does not work, try to address (1), set [[tolrde]] to 0.001, increase the value of [[nline]]
(e.g. 6 or 8 instead of the default 4), as well as of [[nnsclo]]  (e.g. 2 instead of the default).
Your residm should be significantly lower than without such modifications, and perhaps your problem will be solved.

If this still does not work, but your residm did not look bad after all before addressing (1),
then revert back to the default values of [[tolrde]], [[nline]] and [[nnsclo]], as this indeed slows down ABINIT.

Then, try to use [[iscf]] 2.
This is a very slow but unconditionally convergent algorithm provided [[diemix]] is small enough.
At some stage of convergence, you might restart with the obtained DEN a better algorithm,
as non-linear effects should have been eliminated at some stage.
