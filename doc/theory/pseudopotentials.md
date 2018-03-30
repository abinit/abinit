---
authors: MG
---

$
\newcommand{\aa}{\alpha}
\newcommand{\PS}{{\text{PS}}}
\newcommand{\Ylm}{{Y_m^l}}
\newcommand{\Vloc}{V_{\text{loc}}}
\newcommand{\vv}{\hat v}
\newcommand{\vloc}{\vv_{\text{loc}}}
\newcommand{\Vnl}{V_{\text{nl}}}
\newcommand{\lm}{{lm}}
\newcommand{\KK}{{\bf K}}
\newcommand{\KKp}{{\bf{K'}}}
\newcommand{\KKhat}{\widehat \KK}
\newcommand{\KKphat}{\widehat{\KK}'}
\newcommand{\jl}{j_l}
\newcommand{\rrhat}{{\widehat\rr}}
$

## Norm-conserving Pseudopotentials in the Kleynman Bylander form

This page reports basic definitions and equations needed 
to evalute the matrix elements of the commutator between the position operator, $\rr$,
and the nonlocal part of the Hamiltonian in the particular case in which norm-conserving pseudopotentials are used.
Providing a consistent introduction to the theory of pseudopotentials is beyond the purposes of this section,
a more complete discussion of the topic can be found in specialized 
articles available in the literature ~\cite{Hamann1979,Bachelet1982,Kleinman1982,Hamann1989,Troullier1991,Gonze1991,Fuchs1999}.

For our purposes, we can summarize by saying that a pseudopotential in constructed in order to 
replace the atomic all-electron potential such that core states are eliminated 
and valence electrons are described by nodeless pseudo wavefunctions whose 
representation in the Fourier space decays rapidly.

Modern norm-conserving pseudopotentials are constructed following the procedure
described in~\cite{Hamann1979,Hamann1989,Troullier1991} such that the scattering properties of the 
all-electron atom are reproduced around the reference energy configuration, up to first order in energy.
In modern ab initio codes, the fully separable form proposed by Kleynman and Bylander 
is usually employed as it drastically reduces the number of operations 
required for the application of the nonlocal part of the Hamiltonian 
as well as the memory required to store the operator~\cite{ Kleinman1982}.
In the KB form, the interaction between a valence electron and an ion
of type $\aa$ is described by means of the operator:

\begin{equation}\label{eq:KBsingleatom}
\vv_\PS^\aa = 
\vloc^\aa(r) + \sum_\lm |\chi_\lm^\aa\ra\,E_l^\aa\,\la\chi_\lm^\aa|,
% \vloc_\aa(\rr) + \sum_\lm |\Ylm\,u_{l\aa}^\KB\ra E_{l\aa}^\KB \lau_{l\aa}^\KB\,\Ylm|
\end{equation}

where $\vloc^\aa(r)$ is a purely local potential with a Coulomb tail $\gamma/r$.

!!! note

    For simplicity, the discussion is limited to pseudopotentials with a single projector per angular channel.
    The generalization to multi-projector pseudopotentials requirens introducing an additional $n$ index
    in the KB energies i.e. $E_{nl}^\aa$ and in the projectors $\chi_{n\lm}^\aa(\rr)$.

The so-called Kleinman-Bylander energies, $E_l^\aa$,
measure the strength of the nonlocal component with respect to the local part.
The projectors $\chi_\lm^\aa(\rr)$ are 
short-ranged functions, expressed in terms of a complex spherical harmonic $\Ylm(\theta,\phi)$
multiplied by an $l$- and atom-dependent radial function 

\begin{equation}
\label{eq:KB_projectors}
\chi_\lm^\aa(\rr) = \dfrac{1}{r}\,\chi_l^\aa(r)\,\Ylm(\theta,\phi).
\end{equation}

<!--
where 

\begin{equation}\label{eq:KBfunctionU}
u_{l\aa}^\KB (\rr) = 
\frac{\Delta v_{l\aa}(\rr)\,u_{l\aa}^\PS (\rr)}
     {\norm{u_{l\aa}(\rr)\,\Delta_{l\aa}}^{1/2}}
\end{equation}
are localized functions defined in terms of the short-ranged ...
and the pseudo eigenfuncions of the reference atom.
\begin{equation}\label{eq:KBenergy}
E_{l\aa}^\KB = 
\frac{\la u_{l\aa}^\PS \Delta v_{l\aa} | \Delta v_{l\aa} u_{l\aa}^\PS \ra }
     {\la u_{l\aa}^\PS | \Delta v_{l\aa} | u_{l\aa}^\PS \ra}
\end{equation}
-->

The total nonlocal part of the Hamiltonian is obtained 
by summing the different atom-centered contributions.
The final expression in real space reads:

\begin{equation}
\Vnl(\rr_1,\rr_2)= 
\sum_{\substack{\RR \\ \aa\tt_\aa}}
\sum_\lm \la\rr_1-\RR-\tt_\aa|\chi_\lm^\aa\ra E_l^\aa\la\chi_\lm^\aa|\rr_2-\RR-\tt_\aa\ra,
\end{equation}

where $\RR$ runs over the sites of the Bravais lattice, and $\tt_\aa$ over the 
positions of the atoms with the same type $\aa$ located inside the unit cell.
Due to the nonlocality of the operator, its Fourier transform depends
on two separate indices, $\kk+\GG_1$ and $\kk+\GG_2$, instead of the simple difference $\GG_1-\GG_2$.
The shorthand notation $\KK = \kk + \GG$ is used in the following to simplify the derivation.

Using the Rayleigh expansion of planewaves in terms of spherical 
harmonics $\Ylm(\theta,\phi)$ and spherical Bessel functions $\jl(Kr)$:

\begin{equation}\label{eq:PWinYlmPl}
e^{i\KK\cdot\rr} = 
 4\pi\,\sum_\lm i^l \jl(Kr) \Ylm^\*(\KKhat)\,\Ylm(\rrhat),
\end{equation}

the Fourier representation of the projector defined by \ref{eq:KB_projectors} can be expressed as:

\begin{equation}
\la\KK|\chi_\lm^\aa\ra = 
 \frac{4\pi}{\Omega^{1/2}}\, (-1)^l\,\Ylm(\KKhat) 
 \int_0^\infty\, r\jl(Kr)\,\chi_l^\aa(r)r\,\dd  =
 \frac{4\pi}{\Omega^{1/2}}\, (-1)^l\,\Ylm(\KKhat) F_l^\aa(K).
\end{equation}

<!--
%where the form factors $F_l^\aa(K)$ related to the atom of type $\aa$ is defined by
%\begin{equation}\label{wq:defformfactors}
% F_l^\aa(K) \df 
% \frac{\int_0^\infty r\,j_l (Kr)\,u_{l\aa} \Delta v_{l\aa}\,\dd r}
%      {\norm{u_{l\aa} \Delta v_{l\aa}}^{1/2} }
%\end{equation}
-->

Finally, the expression for the total nonlocal operator in reciprocal space is:

\begin{equation}
\label{eq:VnlKBmatrixelements}
\Vnl(\KK,\KKp)= \frac{(4\pi)^2}{\Omega} \sum_{\aa \tt_\aa} \sum_\lm
 \,e^{-i(\KK-\KKp)\cdot\tt_\aa}\,
 \Ylm(\KKhat)\Ylm^\*(\KKphat)\, E_l^\aa F_l^\aa(K) F_l^\aa(K').
\end{equation}
