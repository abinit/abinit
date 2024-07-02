---
authors: SP, MR, MS
---

# Dynamical quadrupoles with DFPT

This tutorial discusses how to compute long-wavelength dynamical quadrupoles using
density functional perturbation theory, taking the specific case of AlAs as an example.  

It is assumed the user has already completed the two tutorials [RF1](/tutorial/rf1) and [RF2](/tutorial/rf2),
and that they are familiar with the calculation of ground state and response properties,
in particular phonons, Born effective charges and dielectric tensor.

The first-order long-range multipolar expansion are dynamical Born effective charge while 
the second-order are called dynamical quadrupoles. 

The main use of computing dynamical quadrupoles is to subtract the long-range electrostatic when doing 
Fourier interpolation of the dynamical matrix [[cite:Royo2020]], perturbed potential [[cite:Brunin2020]], 
and electron-phonon matrix elements [[cite:Ponce2021]].

Extention to the case of 2D materials exists for dynamical matrices [[cite:Royo2021]] and electron-phonon matrices [[cite:Ponce2023]].

This lesson should take about 1.0 hour.

## Formalism

The formalism has been derived in [[cite:Royo2019]] which allows to obtain the quadrupole
tensor as the second q-gradients in the direction $\gamma$ of the long-wavelength multipolar expansion
for an applied electric field in the direction $\delta$ and an atomic displacement of atom $\kappa$ in the 
direction $\beta$: 

\begin{equation}
Q_{\kappa\beta}^{(2,\gamma\delta)}= 2 i ( E_{\gamma}^{\mathcal{E}_\delta^* \tau_{\kappa\beta}}
 + E_{\delta}^{\mathcal{E}_\gamma^* \tau_{\kappa\beta}}), 
\end{equation}

where the first $\bf q$-gradient of the mixed response to an electric field and an atomic displacement is obtained for spin $s$ as:

\begin{equation}
 E_{\gamma}^{\mathcal{E}_\delta^* \tau_{\kappa\beta}} = s \int_{\rm BZ} [d^3k] \sum_{m} 
 E_{m{\bf k},\gamma}^{\mathcal{E}_\delta^* \tau_{\kappa\beta}}  +  \frac{1}{2} \int_{\Omega} \int K_{\gamma}({\bf r},{\bf r}') n^{\mathcal{E}_\delta} ({\bf r})
 n^{\tau_{\kappa\beta}}({\bf r}') d^3r d^3r',
\end{equation}

where $K_{\gamma}({\bf r},{\bf r}')$ is the $\bf q$-derivative of the Hartree and exchange-correlation kernel in the direction $\gamma$, $n^\lambda$ is the first-order perturbed density due to a generic pertubation $\lambda$:

\begin{equation}
n^{\lambda}({\bf r}) = 2s \int_{\rm BZ} [d^3 k] \, \sum_m
       \langle  u^{(0)}_{m {\bf k}} | {\bf r} \rangle \langle {\bf r} | u^{\lambda}_{m {\bf k}} \rangle.
\end{equation}

 and where 

\begin{equation}
  E_{m{\bf k},\gamma}^{\mathcal{E}_\delta^* \tau_{\kappa\beta}} = 
 \langle u_{m{\bf k}}^{\mathcal{E}_{\delta}} | \partial_{\gamma} \hat{H}^{(0)}_{{\bf k}} | u_{m{\bf k}}^{\tau_{\kappa\beta}} \rangle + \langle u_{m{\bf k}}^{\mathcal{E}_{\delta}}| \partial_{\gamma} \hat{Q}_{{\bf k}}  \hat{\mathcal{H}}_{{\bf k}}^{\tau_{\kappa\beta}}  | 
  u_{m{\bf k}}^{(0)} \rangle +
  \langle u_{m{\bf k}}^{(0)} |  V^{\mathcal{E}_{\delta}} \partial_{\gamma} \hat{Q}_{{\bf k}}| u_{m{\bf k}}^{\tau_{\kappa\beta}} \rangle + \langle  u_{m{\bf k}}^{\mathcal{E}_{\delta}}|  \hat{H}_{{\bf k},\gamma}^{\tau_{\kappa\beta}}  | u_{m{\bf k}}^{(0)} \rangle  -\frac{i}{2}  \langle u_{m{\bf k}}^{(0)}| \partial_{\gamma\delta} \hat{P}_{{\bf k}} |  u_{m{\bf k}}^{\tau_{\kappa\beta}} \rangle .
\end{equation}

where $\partial_{\gamma\delta} \hat{P}_{{\bf k}}$ is the second ${\bf k}$-gradient of the valence-band projector. 


[TUTORIAL_README]

## Quadrupole calculation in fcc AlAs

*Before beginning, you might consider creating a different subdirectory to work in.
Why not create Work_quad ?*

The file *tquad_1.abi* is the input file for the first step
(GS + DFPT perturbations for all the $\qq$-points in the IBZ).
Copy it to the working directory with:

```sh
cd $ABI_TESTS/tutorespfn/Input
mkdir Work_quad
cd Work_quad
cp ../tquad_1.abi .
```

{% dialog tests/tutorespfn/Input/tquad_1.abi %}

This step might be quite time-consuming so you may want to immediately start the job in background with:

```sh
abinit tquad_1.abi > tquad_1.log 2> err &
```

Open the output and look for the Quadrupole data block:

```sh
 Quadrupole tensor, in cartesian coordinates,
 efidir   atom   atdir    qgrdir          real part        imaginary part
    1       1       1       1           -0.0000000566        -0.0000000000
    1       1       2       1           -0.0000000157        -0.0000000000
    1       1       3       1            0.0000000157        -0.0000000000
    1       2       1       1           -0.0000001675        -0.0000000000
    1       2       2       1           -0.0000000042        -0.0000000000
    1       2       3       1            0.0000000042        -0.0000000000
    2       1       1       1           -0.0000000362        -0.0000000000
    2       1       2       1           -0.0000000205        -0.0000000000
    2       1       3       1           13.4866068621        -0.0000000000
    2       2       1       1           -0.0000000859        -0.0000000000
    2       2       2       1           -0.0000000816        -0.0000000000
    2       2       3       1           -6.0008872938        -0.0000000000
    3       1       1       1           -0.0000000362        -0.0000000000
    3       1       2       1           13.4866068621        -0.0000000000
    3       1       3       1           -0.0000000205        -0.0000000000
    3       2       1       1           -0.0000000859        -0.0000000000
    3       2       2       1           -6.0008872938        -0.0000000000
    3       2       3       1           -0.0000000816        -0.0000000000
```

Since we are in a binary FCC solid, there are only two independent quadrupoles values given by 
$Q_{\kappa\beta}^{\gamma\delta} = Q_\kappa |\varepsilon_{\beta\gamma\delta}|$ where $\varepsilon_{\beta\gamma\delta}$ is the Levi-Cevita symbol. 


