---
authors: M. Giantomassi, X. Gonze, Y. Suzukawa, M. Mikami
---

# Geometric considerations

## Real space

The three primitive translation vectors are $\RR_1$, $\RR_2$, $\RR_3$.
Their representation in Cartesian coordinates (atomic units) is:

$$ \RR_1 \rightarrow {rprimd(:, 1)} $$
$$ \RR_2 \rightarrow {rprimd(:, 2)} $$
$$ \RR_3 \rightarrow {rprimd(:, 3)} $$

Related input variables: [[acell]], [[rprim]], [[angdeg]].

The atomic positions are specified by the coordinates ${\bf x}_{\tau}$
for $\tau=1 \dots N_{atom}$ where $N_{atom}$ is the number of atoms [[natom]].

Representation in reduced coordinates:

\begin{eqnarray*}
{\bf x}_{\tau} &=& x^{red}_{1\tau} \cdot {\bf R}_{1}
 + x^{red}_{2\tau} \cdot {\bf R}_{2}
 + x^{red}_{3\tau} \cdot {\bf R}_{3} \\
 \tau &\rightarrow& {iatom} \\
 N_{atom} &\rightarrow& {natom} \\
 x^{red}_{1\tau} &\rightarrow& {xred(1,iatom)} \\
 x^{red}_{2\tau} &\rightarrow& {xred(2,iatom)} \\
 x^{red}_{3\tau} &\rightarrow& {xred(3,iatom)}
\end{eqnarray*}

Related input variables: [[xangst]], [[xcart]], [[xred]]

The volume of the primitive unit cell (called *ucvol* in the code) is

\begin{eqnarray*}
\Omega &=& {\bf R}_1 \cdot ({\bf R}_2 \times {\bf R}_3)
\end{eqnarray*}

The scalar products in the reduced representation are valuated thanks to

$$ 
{\bf r} \cdot {\bf r'} =\left(
\begin{array}{ccc}
r^{red}_{1} & r^{red}_{2} & r^{red}_{1}
\end{array}
\right)
\left(
\begin{array}{ccc}
{\bf R}_{1} \cdot {\bf R}_{1} & {\bf R}_{1} \cdot {\bf R}_{2} &
{\bf R}_{1} \cdot {\bf R}_{3} \\
{\bf R}_{2} \cdot {\bf R}_{1} & {\bf R}_{2} \cdot {\bf R}_{2} &
{\bf R}_{2} \cdot {\bf R}_{3} \\
{\bf R}_{3} \cdot {\bf R}_{1} & {\bf R}_{3} \cdot {\bf R}_{2} &
{\bf R}_{3} \cdot {\bf R}_{3}
\end{array}
\right)
\left(
\begin{array}{c}
r^{red \prime}_{1} \\
r^{red \prime}_{2} \\
r^{red \prime}_{3}
\end{array}
\right) 
$$

that is

$$ {\bf r} \cdot {\bf r'} = \sum_{ij} r^{red}_{i} {\bf R}^{met}_{ij} r^{red \prime}_{j} $$

where ${\bf R}^{met}_{ij}$ is the metric tensor in real space:

$$ {\bf R}^{met}_{ij} \rightarrow {rmet(i,j)} $$

## Reciprocal space

The three primitive translation vectors in reciprocal space are
$\GG_1$, $\GG_2$,$\GG_3$

\begin{eqnarray*}
{\bf G}_{1}&=&\frac{1}{\Omega}({\bf R}_{2}\times{\bf R}_{3}) \rightarrow {gprimd(:,1)} \\
{\bf G}_{2}&=&\frac{1}{\Omega}({\bf R}_{3}\times{\bf R}_{1}) \rightarrow {gprimd(:,2)} \\
{\bf G}_{3}&=&\frac{1}{\Omega}({\bf R}_{1}\times{\bf R}_{2}) \rightarrow {gprimd(:,3)}
\end{eqnarray*}

This definition is such that $\GG_i \cdot \RR_j = \delta_{ij}$

!!! warning
    Often, a factor of $2\pi$ is present in definition of
    ${\bf G}_{i}$, but not here, for historical reasons.

Reduced representation of vectors (K) in reciprocal space

$$
{\bf K}=K^{red}_{1}{\bf G}_{1}+K^{red}_{2}{\bf G}_{2}
+K^{red}_{3}{\bf G}^{red}_{3} \rightarrow
(K^{red}_{1},K^{red}_{2},K^{red}_{3}) 
$$

e.g. the reduced representation of ${\bf G}_{1}$ is (1, 0, 0).

!!! important

    The reduced representation of the vectors of the reciprocal space
    lattice is made of triplets of integers.

The scalar products in the reduced representation are evaluated thanks to

$$ 
{\bf K} \cdot {\bf K'}=\left(
\begin{array}{ccc}
K^{red}_{1} & K^{red}_{2} & K^{red}_{1}
\end{array}
\right)
\left(
\begin{array}{ccc}
{\bf G}_{1} \cdot {\bf G}_{1} & {\bf G}_{1} \cdot {\bf G}_{2}
& {\bf G}_{1} \cdot {\bf G}_{3} \\
{\bf G}_{2} \cdot {\bf G}_{1} & {\bf G}_{2} \cdot {\bf G}_{2}
& {\bf G}_{2} \cdot {\bf G}_{3} \\
{\bf G}_{3} \cdot {\bf G}_{1} & {\bf G}_{3} \cdot {\bf G}_{2}
& {\bf G}_{3} \cdot {\bf G}_{3}
\end{array}
\right)
\left(
\begin{array}{c}
K^{red \prime}_{1} \\
K^{red \prime}_{2} \\
K^{red \prime}_{3}
\end{array}
\right) 
$$

that is 

$$ {\bf K} \cdot {\bf K'} = \sum_{ij} K^{red}_{i}{\bf G}^{met}_{ij}K^{red \prime}_{j} $$

where ${\bf G}^{met}_{ij}$ is the metric tensor in reciprocal space:

$$ {\bf G}^{met}_{ij} \rightarrow {gmet(i,j)} $$

## Fourier series for periodic lattice quantities

Any function with the periodicity of the lattice i.e any function fullfilling the property:

$$ u(\rr + \RR) = u(\rr) $$

can be represented with the discrete Fourier series:

$$  u(\rr)= \sum_\GG u(\GG)e^{i\GG\cdot\rr} $$ 

where the Fourier coefficient, $u(\GG)$, is given by:

$$ u(\GG) = \frac{1}{\Omega} \int_\Omega u(\rr)e^{-i\GG\cdot\rr}\dd\rr $$

<!--
This appendix reports the conventions used in this work for the Fourier 
representation in frequency- and momentum-space. 
The volume of the unit cell is denoted with $\Omega$, while $V$ is the
total volume of the crystal simulated employing Born-von K\'arman periodic boundary condition~\cite{Ashcroft1976}.

\begin{equation}\label{eq:IFT_2points_convention}
 f(\rr_1,\rr_2)= \frac{1}{V} \sum_{\substack{\qq \\ \GG_1 \GG_2}}  
 e^{i (\qq +\GG_1) \cdot \rr_1}\,f_{\GG_1 \GG_2}(\qq)\,e^{-i (\qq+\GG_2) \cdot \rr_2} 
\end{equation}

\begin{equation}\label{eq:FT_2points_convention}
 f_{\GG_1\GG_2}(\qq) = \frac{1}{V} \iint_V 
 e^{-i(\qq+\GG_1) \cdot \rr_1}\,f(\rr_1, \rr_2)\,e^{i (\qq+\GG_2) \cdot \rr_2}\dd\rr_1\dd\rr_2 
\end{equation}
-->
