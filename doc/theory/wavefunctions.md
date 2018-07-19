---
authors: X. Gonze, Y. Suzukawa, M. Mikami, MG
---

$$
\newcommand{\mcR}{{\mathcal{R}}}
\newcommand{\omcR}{{\hat\mcR}}
\newcommand{\Atm}{{\bf A}}
\newcommand{\AFMomcR}{{\underline{\omcR}}}
\newcommand{\AFMmcR}{{\underline{\mcR}}}
\newcommand{\HH}{{\hat H}}
\newcommand{\nk}{{n\kk}}
\newcommand{\RRm}{\mcR^{-1}} 
\newcommand{\Rit}{\mcR^{-1T}}
\newcommand{\rmt}{{\rr-\tt}}
\newcommand{\aa}{\alpha}
\newcommand{\bb}{\beta}
\newcommand{\ss}{\sigma}
\newcommand{\mcJ}{{\mathcal{J}}}
\newcommand{\omcA}{{\hat\mcA}}
\newcommand{\thalf}{{\tfrac{1}{2}}}
\newcommand{\mcE}{{\mathcal{E}}}
\newcommand{\mcA}{{\mathcal{A}}}
$$

Representation and convension of one wavefunction
ABINIT data structures and their theoretical justifications.

## Notations and theoretical considerations

A Bloch wavefunction characterized by a wavevector $\kk$ is such that

$$ \psi_{\bf k}({\bf r}) = e^{i{\bf k}\cdot{\bf r}} u_{\bf k}({\bf r}) \,, $$ 

where $u_{\bf k}({\bf r})$ is periodic, that is

$$ u_{\bf k}({\bf r}+{\bf R})=u_{\bf k}({\bf r}) \,, $$ 

where ${\bf R}$ is a vector of the real space lattice.

Its representation by plane waves reads:

\begin{eqnarray*}
u_{\bf k}({\bf r})&=&\sum_{\bf G}u_{\bf k}({\bf G})e^{i{\bf G}\cdot{\bf r}}  \,, \\
\psi_{\bf k}({\bf r})&=&\sum_{\bf G}u_{\bf k}({\bf G})
e^{i ({\bf k}+{\bf G})\cdot{\bf r}} \,.
\end{eqnarray*} 

Normalization of ${u_k}$ reads:

$$ \sum_{\bf G}|u_{\bf k}({\bf G})|^2 = 1 \,.$$

For a **spinor** wavefunction, there is an additional variable, 
the spin $\sigma$ that can take two values, that is $\sigma=\uparrow$ (spin up) 
or $\sigma=\downarrow$ (spin down).
The following relations hold:

$$
u_{\bf k}({\bf r},\sigma) = \sum_{\bf G}u_{\bf k}({\bf G},\sigma) e^{i{\bf G} \cdot {\bf r}} \,,
$$

$$
\psi_{\bf k}({\bf r},\sigma) = \sum_{\bf G}u_{\bf k}({\bf G},\sigma) e^{i({\bf k}+{\bf G})\cdot{\bf r}} \,,
$$

$$
\sum_{\sigma}\sum_{\bf G}|u_{\bf k}({\bf G},\sigma)|^2 = 1 \,.
$$

## Properties of the wavefunctions (scalar case)

For ground-state wavefunctions, we have to solve the Schrödinger equation

$$ H |\psi_{n{\bf k}}\ra = \varepsilon_{n{\bf k}}|\psi_{n{\bf k}}\ra \,, $$

where $H$ is the Hamiltonian operator, $n$ labels the state (or the band), $\varepsilon_{n{\bf k}}$ is the eigenvalue.

As the wavevector labelling of an eigenstate comes from the property

$$ \psi_{\bf k}({\bf r}+{\bf R}) = e^{i{\bf k}{\bf R}} \psi_{\bf k}({\bf r}) \,, $$

in which $\kk$ can be replaced by $\kk + \GG$ where $\GG$ is any reciprocal space lattice vector, we can
*choose* the wavefunctions at $\kk$ and $\kk + \GG$
to be equal, or to make a linear combination of wavefunctions with the same energy. 
This is a choice of **gauge** that does not affect the value of physical observables.
In what follows we prefer to work with the gauge 

$$ \psi_{\kk + \GG}(\rr) = \psi_\kk(\rr) $$

to keep the notation as simple as possible,

The time-reversal symmetry (non-magnetic case) of the Hamiltonian gives the following relation:

$$
\psi_{n{\bf k}}({\bf r}) = \psi^{*}_{n-\bf k}({\bf r}) \,.
$$

For the ${\bf k}$ wavevectors that are half a reciprocal lattice vector
$(2{\bf k}={\bf G}_{0})$, there is a special relationship between 
the Fourier coefficients of the wavefunction:

$$
u_{n{\bf k}}({\bf G}) =
u_{n{\bf k}-{\bf G}_{0}}({\bf G}+{\bf G}_{0}) = 
u_{n-{\bf k}}({\bf G}+{\bf G}_{0}) =
u^{*}_{n{\bf k}}(-{\bf G}-{\bf G}_{0}) \,.
$$

That is, coefficients at $\GG$ and $-{\bf G}-{\bf G}_{0}$ are related.
This will allow to decrease by a factor of 2 the storage space for these
specific ${\bf k}$ points and accelerate CPU-intensive parts such as the 
application of the non-local part and the Fourier transform.

## Properties of the wavefunctions (spinor case)

One must distinguish two classes of Hamiltonians:

1. the Hamiltonian is spin-diagonal,
2. the Hamiltonian mixes the spin components.

In the first class, one finds usual non-spin-polarized, non-spin-orbit
Hamiltonians, in which case the spin up-spin up and spin down-spin down parts of the Hamiltonian
are equal, as well as spin-polarized
Hamiltonian when the magnetic field varies in strength but *not* in direction.
In the second class, one finds Hamiltonians that include the
spin-orbit splitting as well as non-collinear spin systems.

In the first class, the wavefunctions can be made entirely of *either* spin-up components *or* spin-down
components, and treated independently of those made of opposite spin.
This corresponds to [[nsppol]] = 2.
In the second class, one must stay with spinor wavefunctions.
This corresponds to [[nspinor ]] = 2.

These two classes are mutually exclusive. The possibilities are thus:

 nsppol      | nspinor     | wavefunction type
------------ | ----------- | -----------------
   1         |     1       |   scalar wavefunctions
   2         |     1       |   spin-polarized wavefunctions
   1         |     2       |   spinor wavefunctions


The inclusion of spin-orbit coupling in the Hamiltonian requires [[nspinor]] = 2.

## Plane wave basis set sphere

In order to avoid dealing with an infinite number of plane waves to represent Bloch wavefunctions,
one selects those with a kinetic energy lower than some cutoff energy $E_{\rm cut}$.
The set of allowed ${\bf G}$ vectors will be denoted by $\mcS_{\kk}(E_{\rm cut})$ such that

$$
\GG\,\in \mcS_{\kk}(E_{\rm cut}) \;\mbox{if}\; \dfrac{\lvert{\bf k + G}\rvert^{2}}{2} \leq E_{\rm cut} \,.
$$

The kinetic energy cutoff is computed from the input variables [[ecut]] and [[dilatmx]]
to give the *effective* value

$$ {\text{ecut_eff}} = {\text ecut} * ({\text dilatmx})^2 $$

For special $\kk$-points satisfying the condition $2 \kk = \GG_0$, not all coefficients must be stored. 
A specific storage mode, governed by the input variable [[istwfk]] has been
introduced for the following $\kk$ points:

$$
\Bigl(0, 0, 0\Bigr), 
\left(0, 0, \frac{1}{2}\right),
\left(0, \frac{1}{2}, 0\right),
\left(0, \frac{1}{2}, \frac{1}{2}\right),
\left(\frac{1}{2}, 0, 0,\right),
\left(\frac{1}{2}, 0, \frac{1}{2}\right),
\left(\frac{1}{2}, \frac{1}{2}, 0\right),
\left(\frac{1}{2}, \frac{1}{2}, \frac{1}{2}\right)
$$

For these points, the number of $\GG$ vectors to be taken into account, is decreased by about a factor of 2.
For the $\GG$'s that are not treated, the coefficients
$u_{n{\bf k}}({\bf G})$ can be recovered from those that are treated, thanks to

$$ u_{n{\bf k}}({\bf G}) = u^{*}_{n{\bf k}}(-{\bf G}-{\bf G}_0) $$

The value of [[istwfk]] is automatically computed by the code
on the basis of the k-point coordinates and the treatment of time-reversal symmetry as specified by [[kptopt]].
One can disable the time-reversal trick in the input file by setting explicitly the value of [[istwfk]]
with the syntax:

```
istwfk *1
```

<!--
The number of plane waves is {\text npw}
For ${\text ipw}=1\cdots {\text npw}$, the reduced coordinates of ${\bf G}$
are contained in the array {\text kg}:

$$
\mbox{these are integer numbers}
\cases{
  {\bf G}^{red}_{1}=& {\text kg}(1,{\text ipw}) \cr
  {\bf G}^{red}_{2}=& {\text kg}(2,{\text ipw}) \cr
  {\bf G}^{red}_{3}=& {\text kg}(3,{\text ipw}) \cr
}
$$
This list of $\GG$ vectors is computed in the routine {\text kpgsph.f}.
-->

The maximum number of $\GG$-vectors over all k-points is called [[mpw]] inside the code. 

## FFT grid and FFT box

For the generation of the density from wavefunctions, as well as for
the application of the local part of the potential, one needs to be
able to compute $u_{n\kk}(\rr)$
for a 3D-mesh of $\rr$-points, extremely fast, from the values $u_{n\kk}(\GG)$.

!!! note

    spin up and spin down parts can be treated separately in this
    operation, so they do not need to be specified otherwise in this section.

The FFT algorithm starts from values of a function

$$
\begin{aligned}
  z (j_{1},j_{2},j_{3})  \, \mbox{for} \, &j_{1}=0\cdots(N_{1}-1) \,, \\
               &j_{2}=0\cdots(N_{2}-1) \,,
               &j_{3}=0\cdots(N_{3}-1)
\end{aligned}
$$

and compute fast the transformed

$$
\begin{aligned}
\tilde{z}(l_{1},l_{2},l_{3}) \, \mbox{for} \, 
                                &l_{1}=0\cdots(N_{1}-1)\,, \\
				&l_{2}=0\cdots(N_{2}-1)\,, \\
				&l_{3}=0\cdots(N_{3}-1)  \\
\end{aligned}
$$

with

$$
\tilde{z}(l_{1},l_{2},l_{3})=\sum_{j_{1},j_{2},j_{3}} z(j_{1},j_{2},j_{3})
e^{i2\pi\left(\frac{j_{1}l_{1}}{N_{1}}+\frac{j_{2}l_{2}}{N_{2}}+\frac{j_{3}l_{3}}{N_{3}}\right)} \,.
$$

We want the values of $u_{\bf k}({\bf r})$ on a FFT grid with: 

\begin{eqnarray*}
r^{red}_{1}&=&\frac{0}{N_{1}},\frac{1}{N_{1}},\cdots
\frac{N_{1}-1}{N_{1}}\left(=\frac{l_{1}}{N_{1}}\right) \,, \\
r^{red}_{2}&=&\frac{0}{N_{2}},\frac{1}{N_{1}},\cdots
\frac{N_{2}-1}{N_{2}}\left(=\frac{l_{2}}{N_{2}}\right) \,, \\
r^{red}_{3}&=&\frac{0}{N_{3}},\frac{1}{N_{3}},\cdots
\frac{N_{3}-1}{N_{3}}\left(=\frac{l_{3}}{N_{3}}\right) \,. 
\end{eqnarray*}

The FFT algorithm has a computational cost that scales as $N\log(N)$ where $N$ is the 
total number of points in the mesh.
Note, however, that we cannot obtain $u_\kk(\rr)$ *everywhere* but only
on the discrete set of FFT points given in the above equations.
The choice of $N_{1},N_{2},N_{3}$ is not discussed here.

<!--
The effect of $G^{red}_{1}$ or $G^{red}_{1}+N_{1}$ (or any value of $G^{red}_{1}$ modulo $N$) will be similar.

\begin{eqnarray*}
u_{{\bf k}}({\bf r})&=&\sum_{\bf G} u_{\bf k}({\bf G})
e^{i2\pi {\bf G} \cdot {\bf r}} \\
 &=&\sum_{\bf G} u_{\bf k}({\bf G}) e^{i2\pi(G^{red}_{1}r^{red}_{1} +
G^{red}_{2}r^{red}_{2} + G^{red}_{3}r^{red}_{3})}
\end{eqnarray*}

Let us represent $u_{\kk}(\rr)$ by the segment

```fortran
wf_real (1:2,1:N_{1},1:N_{2},1:N_{3})
```

where the first index refer to the real or imaginary part and the
three others to the integer values $l_{1}+1,l_{2}+1,l_{3}+1$

Let us map the $c_{\bf k}({\bf G})$ coefficients on a similar segment

```fortran
wf_reciprocal (1:2,1:N_{1},1:N_{2},1:N_{3})
```

with a similar meaning of {\tt wf\_reciprocal}$(1:2,j_{1}+1,j_{2}+1,j_{3}+1)$:

\begin{eqnarray*}
j_{1}&=&{\tt mod}({\bf G}^{red}_{1},N_{1}) [\Rightarrow j_{1}\in[0,N_{1}-1]]\\
j_{2}&=&{\tt mod}({\bf G}^{red}_{2},N_{2}) \\
j_{3}&=&{\tt mod}({\bf G}^{red}_{3},N_{3})
\end{eqnarray*}

Then:

\begin{eqnarray*}
\lefteqn{
{\tt wf\_real}(\cdot ,l_{1}+1,l_{2}+1,l_{3}+1)}  \\
&=& \sum^{N_{1}-1}_{j_{1}=0}
\sum^{N_{2}-1}_{j_{2}=0} \sum^{N_{3}-1}_{j_{3}=0} {\tt wf\_reciprocal}
(\cdot ,j_{1}+1,j_{2}+1,j_{3}+1) \times
e^{i2\pi(\frac{j_{1}l_{1}}{N_{1}}+\frac{j_{2}l_{2}}{N_{2}}+\frac{j_{3}l_{3}}{N_{3}})}
\end{eqnarray*}

This is, up to the array indexing convention, precisely the operation done by the FFT algorithm.
-->

For FFT efficiency (minimisation of cache conflicts), the arrays
are not dimensioned with $(N_1, N_2, N_3)$ but with $(N_4, N_5, N_6)$ where:

\begin{eqnarray*}
N_4 = N_1 + 1\;{\text{if}}\; N_1\; {\text{is even else}}\; N_1, \\
N_5 = N_2 + 1\;{\text{if}}\; N_2\; {\text{is even else}}\; N_2, \\
N_6 = N_3 + 1\;{\text{if}}\; N_3\; {\text{is even else}}\; N_3.
\end{eqnarray*}

The FFT mesh is given by [[ngfft]]. 
PAW requires an additional **dense** FFT mesh for densities and potentials called [[ngfftdg]].

<!--
## Wavefunctions and spatial symmetries

If some spatial symmetry operation commutes with the Hamiltonian:

$$ [H,S_{\bf t}]=0 $$

then

\begin{eqnarray*}
H|\psi\ra = \varepsilon|\psi>&\Rightarrow&
S_{\bf t}H|\psi\ra = \varepsilon S_{t}|\psi\ra \\
 &\Rightarrow& H[S_{\bf t}|\psi\ra] = \varepsilon[S_{\bf t}|\psi]
\end{eqnarray*}

$S_{\bf t}|\psi\ra$ is also an eigenvector, with the same eigenvalue as $|\psi\ra$.

However its wavevector is different:

\begin{eqnarray*}
\psi_{n{\bf k}}({\bf r}+{\bf R}) &=&
 e^{i2\pi {\bf k} {\bf R}} \psi_{n{\bf k}}({\bf r}) \\
\Rightarrow (S_{\bf t} \psi_{n{\bf k}})({\bf r}+{\bf R})&=&
\psi_{n{\bf k}}((S_{\bf t})^{-1}({\bf r}+{\bf R})) \\
 &=&\psi_{n{\bf k}}(\sum_{\beta}S^{-1}_{\alpha\beta}(r_{\beta}+R_{\beta}-t_{\beta})) \\
 &=&\psi_{n{\bf k}}(\sum_{\beta}S^{-1}_{\alpha\beta}(r_{\beta}-t_{\beta})+\sum_{\beta}S^{-1}_{\alpha\beta}R_{\beta}) \\
 &=&\psi_{n{\bf k}}((S_{t})^{-1}({\bf r})+\sum_{\beta}S^{-1}_{\alpha\beta}R_{\beta}) \\
\noalign{\hbox{($S^{-1}_{\alpha\beta} R_{\beta}$ must be a vector of the real space lattice if $S_{t}$ leaves the lattice invariant)}}
 &=&e^{i2\pi \sum_{\alpha\beta} k_{\alpha}
S^{-1}_{\alpha\beta} R_{\beta}} \psi_{n{\bf k}}((S_{t})^{-1}({\bf r})) \\
 &=&e^{i2\pi {\bf k}'\cdot{\bf R}}(S_{\bf t}\psi_{n{\bf k}})({\bf r})
\end{eqnarray*}

where $({\bf k}')_{\alpha} = \sum_{\beta} S^{-1}_{\beta\alpha} k_{\beta}$

For a vector in the reciprocal space

$$ ({\bf k}')_{\beta} = (S_{\bf t}({\bf k}))_{\beta} = \sum_{\beta} S^{-1}_{\beta\alpha} k_{\beta} $$

i.e. the inverse transpose of $S_{\alpha\beta}$ is used.

The preceeding result means

\begin{eqnarray*}
\psi_{n(S^{-1,t}{\bf k})}
 &\stackrel{\rm L.C.}{=}& (S_{t}\psi_{n{\bf k}})({\bf r}) \\
 &\stackrel{\rm L.C.}{=}& \psi_{n{\bf k}} (\sum_{\beta}
 S^{-1}_{\alpha\beta}(r_{\beta}-t_{\beta}))
\end{eqnarray*}

\begin{eqnarray*}
&\Longrightarrow& u_{n(S^{-1,t} k)}({\bf r}) e^{i2\pi \sum_{\alpha\beta}
S^{-1,t}_{\alpha\beta} k_{\beta} r_{\alpha}} \stackrel{\rm L.C.}{=}
e^{i2\pi \sum_{\alpha\beta} k_{\alpha}
S^{-1}_{\alpha\beta}(r_{\beta}-t_{\beta})} \times u_{n{\bf k}}(\sum_{\beta}
S^{-1}_{\alpha\beta}(r_{\beta}-t_{\beta})) \\
&\Longrightarrow& u_{n(S^{-1,t} k)}({\bf r}) \stackrel{\rm L.C.}{=}
e^{-i2\pi \sum_{\alpha\beta} k_{\alpha} S^{-1}_{\alpha\beta}
t_{\beta}} u_{nk}(\sum_{\beta}
S^{-1}_{\alpha\beta}(r_{\beta}-t_{\beta})) \\
&\Longrightarrow& \sum_{{\bf G}} c_{n(S^{-1,t}k )}({\bf G})
 e^{i2\pi{\bf G}\cdot{\bf r}}
\stackrel{\rm L.C.}{=} e^{-i2\pi \sum_{\alpha\beta} k_{\alpha}
S^{-1}_{\alpha\beta} t_{\beta}} \sum_{{\bf G}'} c_{n{\bf k}}({\bf G}')
 e^{i2\pi \sum_{\alpha\beta} G'_{\alpha}S^{-1}_{\alpha\beta}
 (r_{\beta}-t_{\beta})} \\
&\Longrightarrow& c_{n(S^{-1,t} k)}(\sum_{\alpha} G'_{\alpha}
S^{-1}_{\alpha\beta}) \stackrel{\rm L.C.}{=} e^{-i2\pi
\sum_{\alpha\beta}(k_{\alpha}+G'_{\alpha}) S^{-1}_{\alpha\beta}
t_{\beta}} c_{n{\bf k}}({\bf G}') \\
\end{eqnarray*}

This formula allows to derive coefficients $c_{n}$ at one ${\bf k}$ point
from these at a symmetric ${\bf k}$ point.
-->

## Symmetry Properties

### Effect of space group symmetries on electron energies and wavefunctions

All the symmetry operations which leave the crystal unchanged constitute a space group.
Besides the translation symmetry operations, a space group contains proper or 
improper rotations followed by an appropriate fractional displacement.
There are 230 in total [[cite:Bassani1975]]. 

A generic element of the space group will be denoted in the following with $\hat\mcR_\tt$
where $\mcR$ is a real orthonormal matrix associated to
a proper or improper rotation while $\tt$ is the corresponding fractional translation.
If all the fractional translations are zero the space group is said to be *symmorphic*.
In Abinit, the rotations in real space are called [[symrel]], while the corresponding 
fractional translation are stored in the [[tnons]] array.
Note that both *symrel* and *tnons* are given in reduced coordinates. 

The application of a symmetry operation $\hat\mcR_\tt$ to a vector $\Atm_i$ 
defining the position of an atom in the unit cell gives:

\begin{equation}
\omcR_\tt\, \Atm_i \equiv \mcR^{-1} (\Atm_i -\tt) = \Atm_j + \RR,
\end{equation}

where $\Atm_j$ indicates an atom in the same unit cell of the same type as 
$\Atm_i$ which may be coincident with $\Atm_i$, and $\RR$ is a suitable lattice translation (possibly zero).
The application of the symmetry operation $\hat \mcR_\tt$ on a generic function $F(\rr)$ 
of the three spatial coordinates is conventionally defined by:

\begin{equation}
\omcR_\tt\, F (\rr) \equiv F(\Ri (\rr-\tt)).  
\end{equation}

Since $\omcR_\tt$ commutes with the Hamiltonian $\HH$ of the crystal, it readily
follows that, given $\Psi_{n\kk}(\rr)$ being eigenstate of $\HH$,
 $\omcR_\tt\, \Psi_{n\kk}(\rr)$ is also eigenstate of the Schrödinger problem with the same eigenvalue:

\begin{equation}
 \begin{cases}
  \bigl[
   \hat\mcR_\tt, \HH \bigr]  =0 &  \\
   \\
   \HH \Psi_{n\kk}(\rr) = \ee_{n\kk}\Psi_{n\kk}(\rr) &
 \end{cases}
 \Longrightarrow \HH\hat\mcR_\tt\, \Psi_{n\kk}= 
 \omcR_\tt\HH\,\Psi_{n\kk} = \ee_{n\kk}\,\omcR_\tt\Psi_{n\kk}.
\end{equation}

Although $\omcR_\tt\,\Psi_{n\kk}(\rr)$ has the same eigenenergy as $\Psi_{n\kk}(\rr)$, its 
crystalline momentum is different.
The operation $\hat\mcR_\tt$ of the space group transforms a Bloch function with vector $\kk$
into a new Bloch wave function of crystalline momentum $\mcR \kk$. 
This important property can be seen as follows:

\begin{equation}
\label{eq:Rotation_of_psi}
 \begin{aligned}
 \Bigl[ \omcR_\tt\Psi_\nk \Bigr] (\rr+\RR) =  \quad 
   & \Psi_\nk \bigl( \Ri (\rr+\RR-\tt) \bigr) \\ 
 = & e^{i \kk\cdot\Ri (\rr+\RR-\tt)}\, u_\nk \bigl( \Ri(\rr-\tt) \bigr) \\
 = & \quad e^{i\Rit \kk\cdot\rr}\, \Psi_\nk \bigl( \Ri (\rr-\tt)\bigr ) \\
 = & e^{i\,\mcR \kk \cdot \rr}\, \omcR_\tt \Psi_\nk(\rr),
 \end{aligned}
\end{equation}

where $\RR$ is an arbitrary vector of the direct Bravais lattice and the invariance of the periodic part of 
the Bloch wave function has been exploited.

!!! note 

    The last equality in follows from the orthogonality
    of the $\mcR$ matrix when referred to in a Cartesian frame of reference. 
    The orthogonality of the rotation matrix does not hold anymore if, as usually done in the ABINIT code,
    the symmetry operations are expressed in reduced coordinates.
    In this case, the correct matrix to use for operations in reciprocal space is given by the transpose of $\Ri$.
    (called *symrec* in the Fortran code)

For a nondegenerate state, one obtains:

\begin{equation}
\label{eq:rotation_wfn}
\Psi_{\mcR\kk}(\rr) = \Psi_{\kk} (\Ri(\rr-\tt)).
%u_{\mcR^{-1}\kk}(\GG)    & = u_\kk(\mcR^{-1}\GG)
\end{equation}

The above equation can be used to reduce the number of $\kk$-points and 
matrix elements that have to be evaluated and stored since the information in
the full Brillouin zone can be reconstructed by symmetry from an appropriate irreducible wedge.
The set of equations below summarizes the most useful relationships commonly employed:

\begin{eqnarray}
\label{eq:space_group_symmetry}
\begin{cases}
\ee_{\mcR\kk}    & =  \quad \ee_{\kk} 
\\
u_{\mcR\kk}(\rr) & =  \quad e^{-i \mcR\kk     \cdot\tt}\, u_{\kk}\bigl(\RRm(\rmt)\bigr) 
\\
u_{\mcR\kk}(\GG) & =  \quad e^{-i(\mcR\kk+\GG)\cdot\tt}\, u_{\kk}(\RRm\GG).
\end{cases}
\end{eqnarray}

The time invariance of the Hamiltonian, $\HH^\* = \HH$, might introduce additional constraints:

\begin{eqnarray}
\label{eq:time_reversal_symmetry}
\begin{cases}
\ee_\nk   & =  \quad \ee_{-\kk}, 
\\
u_\nk(\rr)& =  \quad u^\*_{n-\kk}(\rr), 
\\
u_\nk(\GG)& =  \quad u_{n-\kk}^\*(-\GG).
\end{cases}
\end{eqnarray}

It is important to stress that the set of equations in \ref{eq:space_group_symmetry}
hold only in case of nondegenerate states.
In the presence of degeneracy, the application of a symmetry operation on a Bloch 
wave function with momentum $\kk$ belonging to the set of degenerate states $\mcC_{\nk}$, 
produces a new eigenstate with same energy and crystalline momentum $\mcR \kk$. 
The new eigenstate is given by an appropriate linear combination of the degenerate states with wave vector $\mcR\kk$. 
More explicitly:

\begin{equation}
\label{eq:symmetry_for_degenerate_states}
\omcR_\tt\Psi_\nk (\rr)  = \Psi_\nk(\mcR^{-1}(\rr-\tt))
= \sum_{\aa \in\,\mcC_\nk} D_{\aa n}(\mcR)\,\Psi_{\aa\mcR\kk}(\rr),
%\psi_{n\mcR\kk}(\rr)% u_{\mcR^{-1}\kk}(\GG)    & = u_\kk(\mcR^{-1}\GG)
\end{equation}

where $D_{\aa\bb}(\mcR)$ is the unitary transformation associated with $\mcR$ [[cite:Bassani1975]]. 
Equation \ref{eq:symmetry_for_degenerate_states} reduces to \ref{eq:rotation_wfn}
if the state $\Psi_\nk(\rr)$ is nondegenerate since, in this particular case, $D_{\aa\bb}(\mcR) = \delta_{\aa\bb}$.

### Magnetic Space Groups

An extension of the concept of symmetry is needed in order to explain the magnetic properties 
of crystals [[cite:Landau1984]].
In addition to the spatial arrangement of atoms, the orientation of the atomic magnetic moment 
also becomes important in such systems.
It may turn out that the usual spatial operation, while leaving the crystal unchanged in 
regard to its geometrical structure, will reverse the orientation of spins.

Let the symbol $\mcJ$  indicate the operation of reversing all spins, and 
let $\mcE$ denote the identity operator. 
A combined operation consisting of an ordinary symmetry operation followed by $\mcJ$
is a new type of symmetry operation called a complementary operation. 
The rules for operator composition in the standard nonmagnetic group can be easily extended to include $\mcJ$:

\begin{equation}
\text{If}\;\; \omcA_1 \cdot \omcA_2 = \omcA_3 \;\;\text{then}
\begin{cases}
(\mcJ \omcA_1) \cdot (\mcJ\omcA_2) = \omcA_3
\\\\
(\mcJ \omcA_1) \cdot \omcA_2 = \mcJ\omcA_3
\\\\
\omcA_1 \cdot (\mcJ\omcA_2) = \mcJ\omcA_3.
\end{cases}
\end{equation}

Magnetic groups are obtained by replacing some of the symmetry elements of the initial non-magnetic 
space group by their complementary operations so that the resulting ensemble forms a group 
with respect to the new algebra.
Magnetic space groups are sometimes referred to as Shubnikov groups, and can be 
classified according to three different categories [[cite:Bhagavantam1964]].

Shubnikov type IV 
:   Groups which include $\mcJ$ explicitly, so-called grey groups or Shubnikov type IV (230 in number).
    Each group can be obtained by taking the direct product of each of the conventional space 
    groups with the group $(\mcE,\mcJ)$

Shubnikov type III
:   Groups which do not include $\mcJ$ explicitly, but which contain complementary symmetry 
    operations. Also called mixed groups or Shubnikov type III. There are 1191 in total.

Colorless 
:   Groups which do not include $\mcJ$ either explicitly or in conjunction with a conventional 
    symmetry operation, also named colorless groups. There are 230 and they are indistinguishable from 
    the conventional space groups

Henceforth, in order to keep the notation as simple as possible, a symmetry operation containing $\mcJ$ 
will be denoted with an underlined symbol $\AFMomcR_\tt$, and will be said to have anti-ferromagnetic character.
On the contrary, symmetry operations which preserve the orientation of the spin projection will be denoted using the 
standard notation $\omcR_\tt$, and will be said to have ferromagnetic character.
The action of a magnetic symmetry on a nondegenerate Bloch state has to be generalized according to:

\begin{eqnarray}
\label{eq:AFM_symmetries_wfs}
\omcR_\tt\,    \Psi^\ss_\kk(\rr)  \equiv & \Psi^\ss_\kk(\Ri(\rr-\tt))    & = \Psi_{\mcR\kk}^\ss(\rr),
\\
\AFMomcR_\tt\, \Psi^\ss_\kk(\rr)  \equiv & \Psi^{-\ss}_\kk(\Ri(\rr-\tt)) & = \Psi_{\mcR\kk}^{-\ss}(\rr),
\end{eqnarray}

where $\ss = \pm \thalf$. 
Using the above symmetry relations, 
one can verify that the two spin components of the electron density transform according to:

\begin{eqnarray}
\label{eq:rho_vxc_AFM_symmetry}
\omcR_\tt\,    n^\ss(\rr) =  & n^\ss \bigl(\mcR^{-1}(\rr-\tt)\bigr) & =  n^\ss   (\rr),
\\
\AFMomcR_\tt\, n^\ss(\rr) =  & n^\ss \bigl(\mcR^{-1}(\rr-\tt)\bigr) & =  n^{-\ss}(\rr),
\end{eqnarray}

from which it follows that the total charge and, consequently, the Hartree potential are invariant
under the application of any operation of the magnetic group.
Similar symmetry relationships hold for the two components of the exchange-correlation potential.

The magnetic character of the symmetry operation is stored in the [[symafm]] array.
