---
authors: X. Gonze, Y. Suzukawa, M. Mikami
---

# Geometric considerations

## Real space

The three primitive translation vectors are ${\bf R}_{1p}$, ${\bf R}_{2p}$, ${\bf R}_{3p}$.
Their representation in Cartesian coordinates (atomic units) is:

$$ {\bf R}_{1p} \rightarrow {rprimd(1:3, 1)} $$
$$ {\bf R}_{2p} \rightarrow {rprimd(1:3, 2)} $$
$$ {\bf R}_{3p} \rightarrow {rprimd(1:3, 3)} $$

Related input variables: [[acell]], [[rprim]], [[angdeg]].

The atomic positions are specified by the coordinates ${\bf x}_{\tau}$
for $\tau=1 \dots N_{atom}$ where $N_{atom}$ is the number of atoms [[natom]].

Representation in reduced coordinates:

\begin{eqnarray*}
{\bf x}_{\tau} &=& x^{red}_{1\tau} \cdot {\bf R}_{1p}
 + x^{red}_{2\tau} \cdot {\bf R}_{2p}
 + x^{red}_{3\tau} \cdot {\bf R}_{3p} \\
 \tau &\rightarrow& {iatom} \\
 N_{atom} &\rightarrow& {natom} \\
 x^{red}_{1\tau} &\rightarrow& {xred(1,iatom)} \\
 x^{red}_{2\tau} &\rightarrow& {xred(2,iatom)} \\
 x^{red}_{3\tau} &\rightarrow& {xred(3,iatom)}
\end{eqnarray*}

Related input variables: [[xangst]], [[xcart]], [[xred]]

The volume of the primitive unit cell (ucvol) is

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
{\bf R}_{1p} \cdot {\bf R}_{1p} & {\bf R}_{1p} \cdot {\bf R}_{2p} &
{\bf R}_{1p} \cdot {\bf R}_{3p} \\
{\bf R}_{2p} \cdot {\bf R}_{1p} & {\bf R}_{2p} \cdot {\bf R}_{2p} &
{\bf R}_{2p} \cdot {\bf R}_{3p} \\
{\bf R}_{3p} \cdot {\bf R}_{1p} & {\bf R}_{3p} \cdot {\bf R}_{2p} &
{\bf R}_{3p} \cdot {\bf R}_{3p}
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
${\bf G}_{1p}$,${\bf G}_{2p}$,${\bf G}_{3p}$

\begin{eqnarray*}
{\bf G}_{1p}&=&\frac{1}{\Omega}({\bf R}_{2p}\times{\bf R}_{3p}) \rightarrow {gprimd(1:3,1)} \\
{\bf G}_{2p}&=&\frac{1}{\Omega}({\bf R}_{3p}\times{\bf R}_{1p}) \rightarrow {gprimd(1:3,2)} \\
{\bf G}_{3p}&=&\frac{1}{\Omega}({\bf R}_{1p}\times{\bf R}_{2p}) \rightarrow {gprimd(1:3,3)}
\end{eqnarray*}

This definition is such that ${\bf G}_{ip}\cdot{\bf R}_{jp}=\delta_{ij}$

!!! warning
    Often, a factor of $2\pi$ is present in definition of
    ${\bf G}_{ip}$, but not here, for historical reasons.

Reduced representation of vectors (K) in reciprocal space

$$
{\bf K}=K^{red}_{1}{\bf G}_{1p}+K^{red}_{2}{\bf G}_{2p}
+K^{red}_{3}{\bf G}^{red}_{3p} \rightarrow
(K^{red}_{1},K^{red}_{2},K^{red}_{3}) 
$$

e.g. the reduced representation of ${\bf G}_{1p}$ is (1,0,0).

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
{\bf G}_{1p} \cdot {\bf G}_{1p} & {\bf G}_{1p} \cdot {\bf G}_{2p}
& {\bf G}_{1p} \cdot {\bf G}_{3p} \\
{\bf G}_{2p} \cdot {\bf G}_{1p} & {\bf G}_{2p} \cdot {\bf G}_{2p}
& {\bf G}_{2p} \cdot {\bf G}_{3p} \\
{\bf G}_{3p} \cdot {\bf G}_{1p} & {\bf G}_{3p} \cdot {\bf G}_{2p}
& {\bf G}_{3p} \cdot {\bf G}_{3p}
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

## Symmetries

A symmetry operation in real space sends the point ${\bf r}$ to the point
${\bf r'}={\bf S_t}\{{\bf r}\}$ whose coordinates are
$({\bf r'})_{\alpha}=\sum_{\beta}S_{\alpha \beta}r_{\beta} + t_{\alpha}$ (Cartesian coordinates).

The symmetry operations that preserves the crystalline structure are
those that send every atom location on an atom location with the same atomic type.

The application of a symmetry operation to a function of spatial coordinates $f$ is such that:

$$
({\bf S_t}f)({\bf r}) = f(({\bf S_t})^{-1}{\bf r})
$$

$$
({\bf S_t})^{-1}({\bf r}) =
\sum_{\beta} S^{-1}_{\alpha \beta}(r_{\beta}-t_{\beta})
$$

For each symmetry operation,$isym=1 \dots nsym$, the $3\times3$
${\bf S}^{red}$ matrix is stored in {symrel(:,:,isym)}.

[in reduced coordinates:
 $r'^{red}_{\alpha}=\sum_{\beta}S^{red}_{\alpha \beta}
 r^{red}_{\beta}+t^{red}_{\beta}$ ]

and the vector ${\bf t}^{red}$ is stored in {tnons(:,isym)}.

The conversion between reduced coordinates and Cartesian coordinates is

$$
r'_{\gamma}=\sum_{\alpha \beta}(R_{\alpha p})_{\gamma}
[S^{red}_{\alpha \beta} r^{red}_{\beta}+t^{red}_{\alpha}]
$$

with $[{\rm as} \ G_{ip} \cdot R_{jp}=\delta_{ij}]$

$$
r_{\delta}=\sum_{\alpha}(R_{\alpha p})_\delta
r^{red}_{\alpha} \rightarrow
\sum_{\beta}(G_{\beta p})_{\delta} r_{\delta}=r^{red}_{\beta}
$$.

So 

$$
S_{\gamma \delta}=\sum_{\alpha \beta}(R_{\alpha p})_{\gamma}
S^{red}_{\alpha \beta}(G_{\beta p})_{\gamma}
$$
