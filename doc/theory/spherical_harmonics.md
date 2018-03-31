---
authors: MG
---

$$
\newcommand{\lm}{{lm}}
\newcommand{\Ylm}{{Y_m^l}}
\newcommand{\rY}{{\mathcal Y}}
\newcommand{\rYlm}{{\mathcal Y}_m^l}
\newcommand{\rYLM}{{\mathcal Y}_M^L}
\newcommand{\lp}{{l'}}
\newcommand{\rrhat}{{\widehat\rr}}
\newcommand{\rrphat}{{{\widehat\rr}'}}
\newcommand{\mcP}{{\mathcal{P}}}
\newcommand{\df}{\equiv} %this required mathtools :=  
\newcommand{\li}{{l_i}}
\newcommand{\mi}{{m_i}}
\newcommand{\nj}{{n_j}}
\newcommand{\lj}{{l_j}}
\newcommand{\mj}{{m_j}}
\renewcommand{\mp}{{m'}}
\newcommand{\mcG}{{\mathcal G}}
\newcommand{\omcR}{{\hat\mcR}}
\newcommand{\mcR}{{\mathcal{R}}}
\newcommand{\ddO}{{\dd\Omega}}
\newcommand{\LM}{{LM}}
\newcommand{\mcRm}{\mcR^{-1}} 
\newcommand{\kpGhat}{\widehat{\kpG}}
\newcommand{\rrp}{{\bf r'}}
\newcommand{\Gaunt}{{\mathcal G}}
$$ 

# Complex and Real Spherical Harmonics

## Complex Spherical Harmonics

Complex spherical harmonics, $\Ylm$, are defined as the eigenfunctions of the orbital 
angular momentum operators, $\hat L^2$ and $\hat L_z$. 
They can be labeled by two quantum numbers, $l$ and $m$, related to the eigenvalues of 
$\hat L^2$ and $\hat L_z$: 

\begin{equation}
\hat L^2 \Ylm(\theta,\phi) = l(l+1) \Ylm(\theta,\phi),
\end{equation}

and 

\begin{equation}
\hat L_z \Ylm(\theta,\phi) = m \Ylm(\theta,\phi).
\end{equation}

Their explicit expression is:

\begin{equation}
\Ylm(\theta,\phi) = (-1)^m N_\lm P_m^l(\cos\theta) e^{im\phi},
\end{equation}

where $l$ takes non-negative integer values and the possible values for $m$ are 
integers from $-l$ to $l$. 
The two symbols $\theta$ and $\phi$ are used to denote angular spherical coordinates.
$N_\lm$ is the normalization factor:

\begin{equation}
N_\lm = \sqrt{\frac{2l+1}{4\pi} \dfrac{(l-m)!}{(l+m)!}}
\end{equation}

with $P_m^l$ being the associated Legendre polynomials.
Complex spherical harmonics are orthogonal functions on the sphere, i.e.:

\begin{equation}
\la\Ylm|Y_\mp^\lp\ra_\Omega = \delta_{l\lp}\delta_{m\mp},
\end{equation}

where $\la.|.\ra_\Omega$ denotes the integration over the unit sphere.
Furthermore, they form a complete basis set since any square-integrable 
function $f(\theta,\phi)$ can be expanded in series according to:

\begin{equation}
f(\theta,\phi) = \sum_{l=0}^\infty \sum_{m=-l}^l a_{\lm}\Ylm(\theta,\phi).
\end{equation}

Complex spherical harmonics satisfy two useful properties commonly employed in practical applications.
The parity properties:

\begin{equation}\label{eq:parity_property_Ylm}
{Y_m^l}^\dagger(\rrhat) = (-1)^m Y_{-m}^l (\rrhat),
\end{equation}

\begin{equation}
\Ylm(-\rrhat) = (-1)^l \Ylm(\rrhat),
\end{equation}

and the addition theorem:

\begin{equation}\label{eq:addition_theorem_CSH}
\sum_m \Ylm(\rrhat)\, \Ylm^\*(\rrphat) = 
\dfrac{2l+1}{4\pi}\,\mcP_l (\rrhat\cdot\rrphat),
\end{equation}

where $\mcP_l(q)$ are Legendre polynomials.
Equation \ref{eq:addition_theorem_CSH} is generally employed in the 
application of the nonlocal part of the Hamiltonian in the case of 
norm-conserving pseudopotentials when [[useylm]] = 0. 

## Real spherical harmonics <a id="RSH"><a>

Real spherical harmonics (RSH) are obtained by combining complex conjugate 
functions associated to opposite values of $m$. 
RSH are the most adequate basis functions for calculations in which atomic 
symmetry is important since they can be directly related to the irreducible
representations of the subgroups of $D_3$ [[cite:Blanco1997]].
Moreover, being real, they have half the memory requirement of complex spherical harmonics. 
This is clearly an advantage if high angular momenta are needed or several RHS values have to be stored in memory.
A possible definition for RHS is [[cite:Blanco1997]]:

\begin{equation}\label{eq:Definition_real_harmonics}
 \rYlm \df 
 \begin{cases}
  \quad \dfrac{(-1)^m}{\sqrt{2}}\,(\Ylm + \Ylm^\*) \quad  m > 0
  \\\\ 
  \quad Y_0^l  \quad m = 0 
  \\\\
  \quad \dfrac{(-1)^m}{i\sqrt{2}}\,(Y_{\lvert{m}\rvert}^l - {Y_{\lvert{m}\rvert}^l}^\*)  \quad  m < 0.
 \end{cases}
\end{equation}

Equation \ref{eq:Definition_real_harmonics} can be rewritten in matrix notation as

\begin{equation}\label{eq:Unitaty_transformation_Y_Re2Cmplx}
\underline\rY^l = U^l \underline{Y}^l,
\end{equation}

where $\underline\rY^l$ and $\underline{Y}^l$ are column vectors containing real and 
spherical harmonics ordered by increasing $m$ values.
$U^l$ is a unitary matrix whose explicit expression can be found in [[cite:Blanco1997]].

The basic properties of RSH can be easily derived from 
the properties of complex spherical harmonics by means of \ref{eq:Definition_real_harmonics}.

Reality
:
    \begin{equation}
    {\rYlm}^\*(\rrhat) = \rYlm(\rrhat),
    \end{equation}

Parity property
:
    \begin{equation}
    \rYlm(-\rrhat) = (-1)^l \rYlm (\rrhat),
    \end{equation}

Orthonormality
:
    \begin{equation}
    \la\rYlm|\rY_{\mp}^\lp\ra_\Omega = \delta_{l\lp} \delta_{m\mp}.
    \end{equation}

Certain on-site matrix elements include the integration of three real spherical harmonics, which 
can be evaluated exactly:

\begin{equation}\label{eq:Gaunt_def}
\Gaunt_{\li\mi\lj\mj}^{\LM} \df \int \rY_\mi^\li\,\rY_M^L\,\rY_\mj^\lj\ddO.
\end{equation}

The above integral is known in the literature under the name of Gaunt coefficient.
Due to selection rules, the integral in equation \ref{eq:Gaunt_def} is nonzero only if $\mi+\mj = M$
and $\lvert{\li+\lj}\rvert \ge (L,L+2,L+4, \dots) \ge \lvert{\li-\lj}\rvert$.

!!! note

    The expression for $\mcG$ follows the internal conventions 
    used inside the ABINIT code and slightly differs from the one reported in standard text-books.

## Symmetry transformation of RSH

The use of symmetry properties is of fundamental importance to simplify the
determination of electronic structure as it greatly 
reduces the number of non-equivalent integrals that have to be  computed,
as well as the size of any tables stored in memory.

Frequently one has to evaluate the value of a RSH at a rotated point after 
the application of a symmetry operation belonging to the space group of the crystal.
This usually occurs when $\kk$-dependent matrix elements have to be symmetrized in the full 
Brillouin zone starting from the knowledge of their symmetrical image in the irreducible wedge. 
The effect of a proper or improper rotation on a RSH can be deduced from:

\begin{equation}
\label{eq:RSH_rotation}
\omcR\,\rYlm(\rrhat) = \rYlm(\mcRm \rrhat) = 
\sum_\alpha D^l_{\alpha m}(\mcR)\,\rYlm(\rrhat).
\end{equation}

That is, spherical harmonics of given $l$ are transformed into a linear combination
of RHS of same $l$ where the coefficients are given by:

\begin{equation}
D^l_{\alpha m}(\mcR) \df \la\rY^l_\alpha|\hat\mcR|\rYlm\ra.
\end{equation}

## Useful expansions in terms of RSH

By means of the unitary transformation given in \ref{eq:Definition_real_harmonics}, 
it is possible to rewrite the Rayleigh expansion of a plane wave in terms of RSH as

\begin{equation}
\label{eq:Rayleigh_expansion_real_Ylm}
e^{i(\kpG)\cdot\rr} = 
4\pi \sum_\lm i^l j_l(\lvert{\kpG}\rvert r)\,\rYlm(\kpGhat)\,\rYlm(\rrhat)
\end{equation}

with $j_l(x)$ the spherical Bessel function of order $l$. 
In a similar way, the electrostatic potentials can be expanded in a real spherical harmonics basis set.
The final expression is very similar to the one obtained in the case of complex spherical harmonics:

\begin{equation}\label{eq:Expansion_Coulombian_RSH}
\dfrac{1}{\lvert{\rr-\rrp}\rvert} = \sum_{l=0}^\infty 
\dfrac{4\pi}{2l+1} 
\sum_{m=-l}^{m=l}
\dfrac{r_<^l}{r_>^{l+1}} \rYlm(\rrhat)\rYlm(\rrphat),
\end{equation}

with $r_< \df \min(r,r')$ and $r_> \df \max(r,r')$.
