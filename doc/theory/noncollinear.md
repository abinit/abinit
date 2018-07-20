---
authors: G. Zerah
---

## Non colinear magnetism

## Notations and theoretical considerations

We will denote the spinor by $\Psi^{\alpha\beta}$, ${\alpha, \beta}$ being the two spin indexes.
The magnetic properties are well represented by introducing the spin density matrix:

$$ \rho^{\alpha\beta}(\rr) = \sum_n f_n \la \rr|\Psi_n^\alpha\ra \la\Psi_n^\beta|\rr\ra $$

where the sum runs over all states and $f_n$ is the occupation of state $n$.

With $\rho^{\alpha\beta}(\rr)$, we can express the scalar density by

$$ \rho(\rr)=\sum_{\alpha} \rho^{\alpha\alpha}(\rr) $$

and the magnetization density $\vec m(\rr)$ (in units of $\hbar /2$) whose components are:

$$ m_i(\rr) = \sum_{\alpha\beta} \rho^{\alpha\beta}(\rr) \sigma_i^{\alpha\beta}, $$

where the $\sigma_i$ are the Pauli matrices.

In general, $E_{xc}$ is a functional of $\rho^{\alpha\beta}(\rr)$, or equivalently of $\vec m(\rr)$ and $\rho(\rr)$. 
It is therefore denoted as $E_{xc}[n(\rr), \vec m(\rr)]$.

The expression of $V_{xc}$ taking into account the above expression of $E_{xc}$ is:

$$
V_{xc}^{\alpha\beta}(\rr)={\delta E_{xc} \over \delta \rho (\rr)} \delta_{\alpha\beta} +
\sum_{i=1}^3 {\delta E_{xc} \over \delta m_i (\rr) }\sigma_i^{\alpha\beta}
$$

In the LDA approximation, due to its rotational invariance, $E_{xc}$ is indeed a functional of $n(\rr)$ and $|m(\rr)|$ only.
In the GGA approximation, on the contrary, we **assume** that it is a functional of $n(\rr)$ and $|m(\rr)|$ and their gradients.
(This is not the most general functional of $\vec m(\rr)$ dependent upon first order derivatives, and rotationally invariant.)
We therefore use exactly the same functional as in the spin polarized situation, using the local direction
of $\vec m(\rr)$ as polarization direction.

We then have

$$ 
{\delta E_{xc} \over \delta m_i (\rr) }={\delta E_{xc} \over \delta |m_i (\rr)| } \widehat {m(\rr)},
$$

where $\widehat {m(\rr)} = {m(\rr) \over |m(\rr)|}$.
Now, in the LDA-GGA formulations, $n_\uparrow + n_\downarrow =n$ and $|n_\uparrow-n_\downarrow|=|m|$
and therefore, if we set $n_\uparrow = (n+m)/2$ and $n_\downarrow=(n-n_\uparrow)$, we have:

$$
{\delta E_{xc} \over \delta \rho (\rr)} = {1 \over 2} \Bigl(
{\delta E_{xc} \over \delta n_\uparrow(\rr)}+
{\delta E_{xc} \over \delta n_\downarrow(\rr)}
\Bigr )
$$

and

$$
{\delta E_{xc} \over \delta |m (\rr)| }={1 \over 2} \Bigl ( 
{\delta E_{xc} \over \delta n_\uparrow(\rr)} -
{\delta E_{xc} \over \delta n_\downarrow(\rr)}
\Bigr )
$$

This makes the connection with the more usual spin polarized case.

Expression of $V_{xc}$ in LDA-GGA

$$
V_{xc}(\rr) = {\delta E_{xc} \over \delta \rho (\rr)} \delta_{\alpha\beta}+ {\delta E_{xc} \over \delta |m (\rr)| }
 {\widehat m(\rr)}.\sigma
$$

## Implementation

Computation of $\rho^{\alpha\beta}(\rr) = \sum_n f_n \la \rr|\Psi^\alpha\ra \la\Psi^\beta|\rr\ra$

One would like to use the routine *mkrho* which does precisely this
but this routine transforms only real quantities, whereas
$\rho^{\alpha\beta}(\rr)$ is hermitian and can have complex elements.
The *trick* is to use only the real quantities:

\begin{eqnarray*}
\rho^{11}(\rr)& = &\sum_n f_n \la \rr|\Psi^1\ra \la\Psi^1\ra \\
\rho^{22}(\rr)&=&\sum_n f_n \la \rr|\Psi^2\ra \la\Psi^2\ra \\
\rho(\rr)+m_x(\rr)&=&\sum_{n} f_n (\Psi^{1}+\Psi^{2})^*_n (\Psi^{1}+\Psi^{2})_n \\
\rho(\rr)+m_y(\rr)&=&\sum_{n} f_n (\Psi^{1}-i \Psi^{2})^*_n (\Psi^{1}-i \Psi^{2})_n
\end{eqnarray*}

and compute $\rho(\rr)$ and $\vec m(\rr)$ with the help of:

\begin{eqnarray*}
\rho(\rr)&=&\rho^{11}(\rr)+\rho^{22}(\rr) \\
m_z(\rr)&=&\rho^{11}(\rr) - \rho^{22}(\rr)
\end{eqnarray*}


For more information about noncollinear magnetism see [[cite:Hobbs2000]] 
and [[cite:Perdew1992]] for the xc functional.
