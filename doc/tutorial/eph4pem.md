---
authors: BG
---

# Generalized Frohlich model: polaron formation energy and effective mass

This tutorial explains how to compute the polaron formation energy based on the generalized Frohlich formalism.
Additionaly one will learn how to compute polaron effective masses in the case of triply degenerate electronic bands in cubic materials.
For the overview of generalized Frohlich model please refer to [[cite:Miglio2020]].
We start by describing the Frohlich formalism upon which we base our calculations for the polaron formation energy, 
also refered to as Zero Point Renormalization (ZPR).
We continue by analysing the generalized Frohlich model applied in particular to cubic materials, 
for which we are able to provide the polaron effective mass. 
Further details concerning the formalism and implementation related to cubic materials are given in **cite_paper**.

It is assumed the user has already completed the two tutorials [RF1](/tutorial/rf1) and [RF2](/tutorial/rf2),
and that he/she is familiar with the calculation of ground state (GS) and response properties
in particular phonons, Born effective charges and the high-frequency dielectric tensor.
The user should have read the [introduction tutorial for the EPH code](/tutorial/eph_intro)
before running these examples.

This lesson should take about xx hour.

## Formalism 

The generalized Frohlich Hamiltonian is expressed as follows:

\begin{equation}
\begin{split}
    \hat{H}^{gFr} = \sum_{\kk,n} \frac{\sigma \kk^2}{2m^{*}_{n}(\hat{k})} \hat{c}^{+}_{kk,n} \hat{c}^{\phantom{+}}_{\kk,n} + \sum_{\qq,j} \omega_{jLO}(\hat{\qq}) \hat{a}^{+}_{\qq,j} \hat{a}^{\phantom{+}}_{\qq,j} \\
    + \sum_{\kk nn',\qq j} g^{gFr}(\kk nn',\qq j) \hat{c}^{+}_{\kk+\qq,n'} \hat{c}^{\phantom{+}}_{\kk,n} (\hat{a}^{\phantom{+}}_{\qq,j} + \hat{a}^{+}_{-\qq,j} ),
\end{split}
\label{eq:gfr}
\end{equation}

where 

The generalized Frohlich Hamiltonian in the cubic case is expressed as follows:

\begin{equation}
\begin{split}
    \hat{H}^{cFr} = \sum_{\kk,n} \frac{\sigma \kk^2}{2m^{*}_{n}(\hat{k})} \hat{c}^{+}_{kk,n} \hat{c}^{\phantom{+}}_{\kk,n} + \sum_{\qq,j} \omega_{jLO} \hat{a}^{+}_{\qq,j} \hat{a}^{\phantom{+}}_{\qq,j} \\
    + \sum_{\kk nn',\qq j} g^{gFr}(\kk nn',\qq j) \hat{c}^{+}_{\kk+\qq,n'} \hat{c}^{\phantom{+}}_{\kk,n} (\hat{a}^{\phantom{+}}_{\qq,j} + \hat{a}^{+}_{-\qq,j} ),
\end{split}
\label{eq:cfr}
\end{equation}

Taking the extremum eigenvalue as reference,
the bare electronic energy dispersion is described by

\begin{equation}
    E_{n}(\kk)= \sigma\frac{k^{2}}{2m_{n}^*(\kk)},    
\end{equation}

up to quadratic order, with a direction- and band-dependent effective mass $m_{n}(\hat{k})$. Such effective mass fulfills the following eigenvalue equation,

\begin{equation}
    H_{LK}(\kk) \ee_{n}(\hat{k}) = \frac{k^{2}}{2m_{n}^*(\hat{k})}\ee_{n}(\hat{k}),
\label{eqn:eigen_eq_HLK}
\end{equation}


Specifically, in the cubic three-band degenerate case, the Luttinger-Kohn Hamiltonian writes

\begin{equation}
       H_{LK} (\kk ) = 
       \begin{pmatrix}
            Ak_{x}^2 + B(k_{y}^2 + k_{z}^2 )    &                       Ck_{x}k_{y} &                       Ck_{x}k_{z} \\
            Ck_{x}k_{y}                         & Ak_{y}^2 + B(k_{z}^2 + k_{x}^2 )  &                       Ck_{y}k_{z} \\
            Ck_{x}k_{z}                         & Ck_{y}k_{z}                       & Ak_{z}^2 + B(k_{x}^2 + k_{y}^2 )  \\
       \end{pmatrix}            
\label{eq:lhh}
\end{equation}

with three parameters $A$, $B$, and $C$.

Turning on the electron-phonon coupling, the polaron dispersion relation in the degenerate bands case will have the same behaviour, including the same symmetry characteristics, which gives:

\begin{equation}
\begin{split}
    E_{P,n}(\kk) &  = \sigma\frac{k^{2}}{2m_{n}^*(\kk)} + \Sigma_{n}(\kk) \\
                     & \approx \Sigma_{n}(\kk = 0) + \sigma\frac{k^2}{2m_{P,n}^*(\hat{k})} 
\end{split}     
\end{equation}

with

\begin{equation}
    \frac{1}{m^{*}_{P,n}(\hat{k})} = \frac{1}{m^{*}_n(\hat{k})} + \sigma\frac{ d^{2} \Sigma_n(\kk,E_n(\kk) )}{ dk^{2}} \Big|_{k=0,\kk=k.\hat{k}}.
    \label{eqn:m_P_deg}
\end{equation}


\begin{equation}
\Sigma_n(\kk,E_n(\kk)) = 
\frac{1}{\pi} \sum_{j}
\int d^{3}q \frac{\omega_{jLO}}{4\pi q^{2}\epsilon^{*}_{j}} 
\ee_{n}^{\phantom{.}T}(\hat{k}) 
\Bigg( 
\frac{k^{2}}{2m_{n}(\kk)}\mathbb{I}
- H(\kk+\qq)
- \sigma\omega_{jLO} \mathbb{I}
\Bigg)^{-1}
\ee_{n}(\hat{k}).
%\end{split} 
    \label{eqn:Sigma_k_Ek}
\end{equation}


In order to calculate the self-energy in the degenerate bands case, we perform a numerical integration over $\qq$ expressed in spherical coordinates:

\begin{equation}
\Sigma_n(\kk,E_n(\kk)) = 
\frac{1}{\pi} \sum_{j}
\frac{\omega_{jLO}}{\epsilon^{*}_{j}}
\int_0^{\infty}dq
\frac{1}{4\pi q^{2}}
\int_{4\pi}d\hat{q}
\,
\ee_{n}^{\phantom{.}T}(\hat{k}) 
\Bigg( 
\frac{k^{2}}{2m_{n}(\kk)}\mathbb{I}
- H(\kk+\qq)
- \sigma\omega_{jLO} \mathbb{I}
\Bigg)^{-1}
\ee_{n}(\hat{k}).
%\end{split} 
\label{eqn:Sigma_z_k_expanded}
\end{equation}

with 

\begin{equation}
    \qq = q \hat{q}= q 
    \begin{pmatrix}
     \sin{\theta}\cos{\phi} \\
     \sin{\theta}\sin{\phi} \\
     \cos{\theta}
    \end{pmatrix}.
\end{equation}

Regarding the integral over $q \in [0,\infty)$, one could use a homogeneous grid integration method, with the maximal value 
$q_{max}$ tending to infinity. 

However, such homogeneous grid integration approach converges slowly. 
Instead, the semi-infinite domain $q \in [0,\infty)$ can be mapped onto a finite one given that the behaviour of the integrand follows f(q) --> $\frac{1}{\alpha^{2} + \beta^{2}q^{2} }$, for large $q$. 
One performs the following change of variable

\begin{equation}
    q = \Big( \frac{\omega_{LO}}{\gamma}\Big)^{1/2} \tan{\xi} 
\end{equation}

giving

\begin{equation}
\frac{dq}{d\xi} = \Big( \frac{\omega_{LO}}{\gamma} \Big)^{1/2} \frac{1}{\cos^{2}{\xi}}
=\Big( \frac{\omega_{LO}}{\gamma} \Big)^{1/2}
\Big(\frac{\gamma q^2}{\omega_{LO}}+1\Big)
\end{equation}

and

\begin{equation}
    \int_{0}^{\infty} dq f(q) = \int_{0}^{\pi/2} d\xi \frac{dq}{d\xi} f(q(\xi)),
\end{equation}


## Typical workflow for Polaron Effective Mass

A typical workflow for ZPR calculations involves the following steps 
(see the [introductory e-ph tutorial](/tutorial/eph_intro)):

1. **GS calculation** to obtain the WFK and the DEN file.
   The $\kk$-mesh should be dense enough to converge both electronic and vibrational properties.

2. **DFPT calculations** for all the IBZ $\qq$-points corresponding to the *ab-initio* [[ddb_ngqpt]] mesh
   that will be used to perform the Fourier interpolation of the dynamical matrix and of the DFPT potentials.
   In the simplest case, one uses a $\qq$-mesh that is equal to the GS $\kk$-mesh (sub-meshes are also fine)
   and the DFPT calculations can directly start from the WFK produced in step #1.
   Remember to compute $\bm{\epsilon}^{\infty}$, $\bm{Z}^*$ (polar materials) and the dynamical quadrupoles
   $\bm{Q}^*$ as these quantities are needed for an accurate interpolation of phonon frequencies and DFPT potentials.

## Getting started

[TUTORIAL_README]
