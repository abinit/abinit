---
authors: FB & JB
---

# First tutorial on aTDEP

![ATDEP1TOP](atdep1_assets/tatdep1top.gif)

## The 2$^{nd}$ order **effective** Interatomic Force Constants (IFC)

This tutorial shows how to capture anharmonicities by means of an harmonic **Temperature Dependent Effective Potential** (TDEP) by using the ABINIT package. In practice, this requires to obtain the $2^{nd}$ order **effective** IFC. Once obtained, almost all the dynamic (phonons...), elastic (constants, moduli...) and thermodynamic (entropy, free energy...) desired quantities can be derived therefrom.

You will learn:

1. how to launch aTDEP just after an ABINIT simulation, 
2. the meaning and effects of the main input variables, and 
3. how to exploit the data provided in the output files.

You are not supposed to know how to use ABINIT, but you are strongly encouraged to read the following documents:

* User guide: [[pdf:aTDEP_Guide| aTDEP guide]]  
* Theory: [[pdf:aTDEP_Paper|aTDEP paper]] corresponding to the article [[cite:Bottin2020]]

This tutorial should take about 1.5 hour.

[TUTORIAL_READMEV9]

## 1. Summary of the aTDEP method

The **Temperature Dependent Effective Potential** approach has been introduced by O. Hellman *et al.* [[cite:Hellman2011]] in 2011. The purpose of this method is to capture the anharmonic effects in an **effective** way. 

Let us consider that the potential energy $U$ of a crystal can be rewritten as a Taylor expansion around the equilibrium:

$$
U=   U_0 + \sum_{p\ge 1} \frac{1}{p\ !} \sum_{\substack{\alpha_1...\alpha_p \\ i_1...i_p}}\overset{(p)}{\Phi}\vphantom{\Phi}_{i_1...i_p}^{\alpha_1...\alpha_p}\prod_{k=1}^p u_{i_k}^{\alpha_k}
$$

> In this equation, and in the following, the Latin letters in subscripts $i, j, k...$ and the Greek letters in superscripts $\alpha, \beta, \gamma$... will define the atoms and the cartesian directions, respectively.

with U$_0$ the minimum of the potential energy, $u$ the displacement with respect to equilibrium, and $\overset{(p)}{\Phi}$ the $p^{th}$ order IFC, respectively. The aim of the TDEP method is to extract the IFC in an effective way. 

Now, let us assume that :

1. the previous equation is truncated at the $2^{nd}$ order such as : 
$$
U=   U_0 + \frac{1}{2!}\sum_{ij,\alpha\beta} \overset{(2)}{\Theta}\vphantom{\Theta}_{ij}^{\alpha\beta} u_i^\alpha u_j^\beta + 0(u^3)
$$

2. a set of $N_t$ forces $\mathbf{F}_{AIMD}(t)$ and displacements $\mathbf{u}_{AIMD}(t)$ is obtained using *ab initio* molecular dynamic (AIMD) simulations, leading to the following system of equations ($F=-\nabla U$): 

$$
F_{i,AIMD}^\alpha(t)=  - \sum_{j,\beta} \overset{(2)}{\Theta}\vphantom{\Theta}_{ij}^{\alpha\beta} u_{j,AIMD}^\beta(t)
$$

It is then possible to obtain the 2nd order **effective** IFC $\overset{(2)}{\Theta}$ by using a least squares method : $\mathbf{\overset{(2)}{\Theta}} = \mathbf{F} . \mathbf{u}^{-1}$. This fitting procedure modifies the 2nd order IFC by including (in an effective way) the anharmonic contributions coming from the terms above the truncation. Therefore, the IFC are no longer constant and become temperature dependent. That is the reason why we change the notation: in the following, the $\Phi$ will be referred to as the ‘‘true IFC’’ and the $\Theta$ as the ‘‘effective IFC’’.


## 2. A simple case : Al-fcc

Let us begin with the face centered cubic (fcc) phase of aluminum. This one is very simple for many reasons : 

- There is only one atom in the unitcell so we will have only three phonon branches in the spectrum.
- At 0 GPa, the fcc phase is stable from 0 K up to the melting so we do not expect any trouble coming from phonon instabilities.

*Before beginning, you might consider to work in a different subdirectory as for the other tutorials. Why not create Work_atdep1 in \$ABI_TESTS/tutoatdep/Input? You can copy all the input files within.*

```sh
cd $ABI_TESTS/tutoatdep/Input
mkdir Work_atdep1
cd Work_atdep1
cp ../tatdep1_1.* . 
```

### 	2.1 The input files

Let us discuss the meaning of these five files :

#### 		^^2.1.1 The data files *tatdep1_1xred.dat*, *tatdep1_1fcart.dat* and *tatdep1_1etotal.dat*^^

These ones store some data coming from the AIMD simulations : the reduced coordinates, the cartesian forces and the total energy of all the atoms in the supercell, respectively. In the present example, only 20 snapshots are extracted from a very long trajectory with thousands molecular dynamic time steps.

#### 		^^2.1.2 The input file *tatdep1_1.abi*^^ 

{% dialog tests/tutoatdep/Input/tatdep1_1.abi %}

This one lists (in a formated way) all the input variables needed. In the following we will comment each variable and the value used for the present calculation. Many of them have the same meaning as in the ABINIT main code.

- For the unitcell :

Input variable | Meaning
---------------|--------
[[atdep:brav]] | Defines the BRAVais lattice (as defined in the ABINIT code). For the present calculation (fcc) : 7 (cubic) and -3 (face centered).
[[atdep:natom_unitcell]] | Defines the Number of ATOMs in the UNITCELL. For the present calculation : 1 
[[atdep:xred_unitcell]] | Defines the Xyz REDuced coordinates in the UNITCELL. For the present calculation : 0.0 0.0 0.0
[[atdep:typat_unitcell]] | Defines the TYPes of AToms in the UNITCELL. For the present calculation : 1
[[atdep:ntypat]] | Defines the Number of TYPes of AToms. For the present calculation : 1
[[atdep:amu]] | Defines the Atomic masses in Mass Units. For the present calculation (Al) : 26.981539

- For the supercell :

Input variable | Meaning
---------------|--------
[[atdep:rprimd]] | Defines the Dimensional Real space PRMitive vectors of the SUPERCELL. For the present calculation : $\begin{pmatrix} 22.9089998 & 0.0 & 0.0 \\ 0.0 & 22.9089998 & 0.0 \\ 0.0 & 0.0 & 22.9089998 \end{pmatrix}$ 
[[atdep:multiplicity]] | Defines the MULTIPLICITY of the SUPERCELL with respect to the primitive UNICELL. For the present calculation : $\begin{pmatrix} -3 & 3 & 3 \\ 3 & -3 & 3 \\ 3 & 3 & -3 \end{pmatrix}$
[[atdep:natom]] | Defines the Number of ATOMs in the SUPERCELL. For the present calculation : 108
[[atdep:typat]] | Defines the TYPe of AToms in the SUPERCELL. For the present calculation : 108 * 1
[[atdep:temperature]] | Defines the TEMPERATURE of the system. For the present calculation : 900 K

- For the calculation :

Input variable | Meaning
---------------|--------
[[atdep:nstep_max]] | Defines the upper limit in the range of configurations that one wants to use. For the present calculation : 20
[[atdep:nstep_min]] | Defines the lower limit in the range of configurations that one wants to use. For the present calculation : 1
[[atdep:rcut]] | Defines the CUToff Radius used to compute the second order IFCs. For the present calculation : 11.45 ($\approx \frac{22.9089998}{2}$) 

- Optional :

Input variable | Meaning
---------------|--------
[[atdep:enunit]] | Defines the ENergy UNIT used for the phonon spectrum. For the present calculation : 3 (in THz)

#### 		^^2.1.3 The files file *tatdep1_1.files*^^

This one lists the input file name and the root of input and output files :

{% dialog tests/tutoatdep/Input/tatdep1_1.files %}

You can now execute `atdep`:

```sh
atdep < tatdep1_1.files > log 2> err &
```

The code should run very quickly.

### 	2.2 The output files

The ATDEP code writes **many** output files (some of them are available in *$ABI_TESTS/tutoatdep/Refs/). The reason is twofold : to remove all the "details" of the calculations from the main output file and to give all the thermodynamic data in an handable format. Let us detail these output files in the following :

#### 		^^2.2.1 The main output file *tatdep1_1.abo*^^

{% dialog tests/tutoatdep/Refs/tatdep1_1.abo %}

This file reproduces all the steps encountered during the execution of ATDEP. You are strongly adviced to detect all the following sequences in the main output file :

1. begins with the common header of the ABINIT output files
2. echoes all the input variables included in the input file
3. computes some quantities using the available data (the acell of the unitcell,...)
4. generates the **q**-point mesh
5. establishes a correspondence between the atoms in the unitcell, the multiplicity, the symmetries and the atoms in the supercell
6. finds the symmetry operations between atoms and pairs of atoms.
7. computes the number of non-zero independent IFC coefficients at the 1st and 2nd order, for each shell of coordination.
8. computes the constraints (for the IFC), builds the pseudo-inverse and solves the problem
9. lists all the IFC coefficients for each shell, at the 1st and 2nd order
10. writes in some files the dynamical matrix, the phonon spectrum and the vibrational density of states (vDOS)
11. echoes the elastic constants and some elastic moduli
12. computes the energy of the model (the temperature effective potential) and some convergence parameters
13. writes in a file the thermodynamics of the system (see below)
14. finishes with the standard aknowlegment section of ABINIT output files 

#### 		^^2.2.2 The phonon frequencies file *tatdep1_1omega.dat*^^

{% dialog tests/tutoatdep/Refs/tatdep1_1omega.dat %}

You can plot the phonon spectrum. If you use the |xmgrace| tool, launch:

      xmgrace -nxy tatdep1_1omega.dat

You should get this picture :

![Al phonon spectrum](atdep1_assets/tatdep1_1omega.png)

On the Y-axis, you have the frequencies (in THz, see the input file). On the X-axis, you have the **q**-points along a path in the Brillouin Zone (BZ). This one is defined by default and depends on the Bravais lattice.

!!! note
    The path along the BZ can be changed using the [[atdep:bzpath]] input variable. 

The BZ boundaries and all the **q**-points included in the path are available in the *tatdep1_1qpt.dat* file :

{% dialog tests/tutoatdep/Refs/tatdep1_1qpt.dat %}

As you can see, in the present calculation the path is as follows : $\Gamma$ - X - W - X' - K - $\Gamma$ - L. Concerning the results, you can compare the phonon spectrum obtained in this tutorial with the first figure of [this paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.91.135501). As you can see, the overall agreement is very good but not perfect due the too small number of atomic configurations (20 time steps) and the difference between the experimental (80 K) and theoretical (900 K) temperatures. If you perform AIMD simulations at various temperatures and store more than 20 atomic configurations, you will obtain the following picture : 

![Al experimental spectrum](atdep1_assets/tatdep1_1exp.jpg)

!!! warning 
    The *tatdep1_1qpt.dat* file do not be confused with the *tatdep1_1qbz.dat* file which defines the Monkhorst-Pack **q**-point mesh used to compute the vDOS : $g(\omega)=\frac{1}{3N_a}\sum_{s=1}^{3N_a}\sum_{\mathbf{q}\in BZ} \delta(\omega-\omega_s(\mathbf{q}))$ such as $\int_0^{\omega_{max}} g(\omega)d\omega =1$, with $\omega_{max}$ the highest phonon frequency of the system. The vDOS is written in the *tatdep1_1vdos.dat* file. You may plot it to verify that the vDOS is consistent with the phonon spectrum.


#### 		^^2.2.3 The thermodynamic file *tatdep1_1thermo.dat*^^

{% dialog tests/tutoatdep/Refs/tatdep1_1thermo.dat %}

In this file, we print all the thermodynamic data that we can compute by using the phonon spectrum and/or the vDOS. The main quantity is the free energy $\mathcal{F}$. This one can be splitted in two parts:

$$
\mathcal{F}(V,T)=U_0(V)+\mathcal{F}_{\rm vib}(V,T) 
$$

The first part is the cold contribution (at T = 0 K) whereas the second one is the vibrational contribution (with T $\neq$ 0). The cold contribution can be computed using a ground state specific calculation or using the following formulation :

$$
U_0(V)=< U_{\rm MD}(t)-\frac{1}{2!}\sum_{ij,\alpha\beta} \overset{(2)}{\Theta}\vphantom{\Theta}_{ij}^{\alpha\beta} u_i^\alpha u_j^\beta >
$$

The vibrational contributions (free energy $F_{\rm vib}$, internal energy $U_{\rm vib}$, entropy $S_{\rm vib}$ and heat capacity $C_{\rm vib,V}$) can be computed using the vDOS $g(\omega)$ in the harmonic approximation (see the paper of Lee & Gonze [[cite:Lee1995]]) :

$$
  F_{\rm vib} = \frac{3N_a}{\beta} \int_0^{\omega_{max}} \ln \left(2\sinh(\frac{\beta\hbar\omega}{2})\right) g(\omega)d\omega \label{eq:Fvib}
$$

$$
  U_{\rm vib} = 3N_a \int_0^{\omega_{max}} \frac{\hbar\omega}{2}\coth(\frac{\beta\hbar\omega}{2}) g(\omega)d\omega 
$$

$$
C_{\rm vib,V} = 3N_a k_B \int_0^{\omega_{max}} \left(\frac{\beta\hbar\omega}{2\sinh(\frac{\beta\hbar\omega}{2})}\right)^2 g(\omega)d\omega \label{eq:Cv}
$$

$$
  S_{\rm vib} = 3N_a k_B \int_0^{\omega_{max}} \left[\frac{\beta\hbar\omega}{2}\coth(\frac{\beta\hbar\omega}{2}) - \ln \left(2\sinh(\frac{\beta\hbar\omega}{2})\right)\right] g(\omega)d\omega
$$

All these thermodynamic data are computed and written in the *tatdep1_1thermo.dat*. Note that this file is divided in two parts :
 
- the first one is dedicated to the thermodynamic data obtained at the temperature defined by input variable [[atdep:temperature]].
- whereas in the second one, the thermodynamic data are **extrapolated** at all the temperatures using a fixed vDOS.

!!! note
    This last point corresponds to the harmonic approximation (HA). In this case, the phonon frequencies do not depend on the temperature but only on the volume $V$, so we have $\omega_{\rm HA} = \omega(V)$. Using a constant vDOS, it's then possible to compute all the thermodynamic data, whathever the temperature. In this case, the temperature variation of the thermodynamic quantities only comes from the filling of phononic states using the Bose-Einstein statistics. To go beyond, and capture the thermal expansion for example, we can assume that the temperature effects are **implicit** through the variation of the volume. This is the quasi-harmonic approximations (QHA) : $\omega_{\rm QHA}=\omega(V(T))$. If in many cases the QHA gives excellent results, it fails to reproduce an **explicit** variation of the thermodynamic data with respect to the temperature (by definition, the phonon frequencies cannot vary at constant volume, along an isochore). This **explicit** variation comes from anharmonic effects and only be captured by going beyond the second order in the energy expansion. That is the work done by aTDEP, by recasting all the 3rd, 4th... terms of the energy expansion within the 2nd order, in an **effective** way. Since the 2nd order effective IFC now takes into account all these terms, it capture the temperature effects and we have $\omega_{\rm Anh}=\omega(T,V(T)$.

In the *tatdep1_1thermo.dat* file corresponding to the present calculation, several remarks can be done. You can see that the specific heat $C_{\rm vib,V}$ is equal to 2.972 (in $k_B$ units) at $T$ = 900 K. In the second part of this file, you see that this quantity converges towards 3 at high temperature, as expected by the Dulong-Petit law (in this part we are in HA, so this law is fulfilled). This result is consistent with the experimental Debye temperature $\Theta_D \approx$ 400 K ; at $T$ = 900 K the behaviour of aluminum is classic and no longer quantic, since all the phononic states are filled. This can be seen also for another quantity. Plot the vibrational internal energy $U_{\rm vib}^{ \rm HA}$ as the function of temperature (see the second part of the file). And plot also $U_{\rm vib}^{\rm Classic}=3k_B T$ corresponding to the classic formulation (in eV, so use the conversion factor 1 eV = 11 604 K). You will see that the classic limit is achieved between 400 and 600 K, as expected.

### 	2.3 Numerical convergence (accuracy and precision)

Several input variables have strong impact on the convergence of the effective IFC, phonon frequencies and thermodynamic data. Two of them are in the *tatdep1_1.abi* input file (in the "DEFINE_COMPUTATIONAL_DETAILS" section) and others comes from the AIMD simulations. 

####		^^2.3.1 The cutoff radius [[atdep:rcut]]^^

The first one is the number of coordination shells included in the calculation of the 2nd order effective IFC. 

![Shells at 2nd order](atdep1_assets/tatdep1_1shell.jpg)

Let us see again the *tatdep1_1.abo* output file and go to the "SECOND ORDER" section. You will see the list of the five shells included in the present calculation and sorted as a function of the shell radius : 0.0000000000 (the onsite interaction), 5.3997030363 (the 2nd shell), 7.6363332667 (the 3rd), 9.3525600046 (the 4th) and 10.7994060725 a.u. (the 5th).  

In the *tatdep1_1.abi* input file the cutoff radius [[atdep:rcut]] equals to 11.45 (a.u.). Now, we will change this value to 6.0, 8.0 and 10.0 in order to have 2, 3 and 4 shells in the calculation, respectively. To do that, you can replace the third line of the *tatdep1_1.files* file by "Rcut6", then set "rcut 6.0" in the input file *tatdep1_1.abi" and finally launch `atdep`. Repeat this process for "Rcut8" and "Rcut10" and plot all the phonon spectra :  

      xmgrace -nxy Rcut6omega.dat -nxy Rcut8omega.dat -nxy Rcut10omega.dat -nxy tatdep1_1omega.dat

You should get the following picture :

![Al phonon spectrum wrt rcut](atdep1_assets/tatdep1_1rcut.png)

Concerning this very simple case, the frequencies are almost converged with only two shells (the onsite interaction and the 1st shell of coordination). In most situations, this is not the case. Here, we can see that some differences remain for rcut = 6.0 and 8.0 a.u. with respect to higher shell radii. With 4 shells and rcut = 10.0, the phonon spectrum seems to be converged and almost equal to 5 shells and rcut = 11.45. This is confirmed by AIMD simulations with 216 atoms in the supercell and a higher shell radius (not shown here).

!!! warning
    The cutoff radius [[atdep:rcut]] cannot be greater than the shortest dimension of the supercell. Otherwise, the shell will include spurious atomic vibrations. The only way to have a larger cutoff radius is to perform AIMD simulations with a larger supercell/number of atoms.    

####		^^2.3.2 The number of atomic configurations^^

Another key quantity is the number of atomic configurations used in the calculation. This one is defined by the difference between two input variables : [[atdep:nstep_max]] - [[atdep:nstep_min]]. For simplicity, we generally use as input data files (*etot.dat*, *xred.dat* and *fcart.dat*) the whole trajectory coming from the AIMD simulations, with thousands of atomic configurations. So, for an AIMD trajectory with 10 000 time steps including a thermalization over 5 000 time steps, we can set [[atdep:nstep_max]] to 10 000 and [[atdep:nstep_min]] to 5 000. However, all the AIMD time steps are not uncorrelated and 99% of the information coming from them is in general useless.

The number of uncorrelated configurations needed for the calculation is direcly related to the number of non-zero and independent IFC coefficients which has to be computed. At the 2nd order, the whole effective IFC $\mathbf{\Theta}$ is a $3N_a\times 3N_a$ matrix. For instance, in the present calculation with $N_a$ = 108 atoms, the whole IFC has 104 976 coefficients. So, if one wants to obtain them (using a least square method $\mathbf{\Theta} = \mathbf{F} . \mathbf{u}^{-1}$), it would require tens of thousands time steps, which is out of reach (see the seminal article of Hellman and coworkers [[cite:Hellman2011]]), 

Thanks to crystal symmetries, tensor symmetries (of the IFC, of the dynamical matrix, of the elastic tensor...) and invariances (translational and rotational) this huge number can be drastically reduced. For example, in the present calculation, we only need to compute 12 IFC coefficients (see "Total number of coefficients at the second order" in the *tatdep1_1.abo* output file) : 0 for the 1st shell then 3, 2, 4 and 3 coefficients for the higher shells. You can see their value in the output file (have a look at "List of (second order) IFC"). In fact, many of them are zero, symmetric or anti-symmetric, which gives the following picture of the whole IFC :

$$ 
\begin{pmatrix}
\begin{pmatrix} \theta_{21} & 0      & 0       \\  0      & \theta_{21} &  0     \\  0      & 0      & \theta_{21} \end{pmatrix} \quad
\begin{pmatrix} \theta_{22} & \theta_{23} & \theta_{24}  \\ -\theta_{23} & \theta_{22} & \theta_{26} \\  \theta_{25} & \theta_{26} & \theta_{22} \end{pmatrix} \quad  
\begin{matrix}  ... \\ ... \\ ... \end{matrix} \\
\begin{pmatrix} \theta_{22} &-\theta_{23} & \theta_{25}  \\ \theta_{23}  & \theta_{22} & \theta_{26} \\  \theta_{24} & \theta_{26} & \theta_{22} \end{pmatrix} \quad
\begin{pmatrix} \theta_{27} & 0      & 0       \\  0      & \theta_{27} & 0      \\  0      & 0      & \theta_{27} \end{pmatrix} \quad  
\begin{matrix}  ... \\ ... \\ ... \end{matrix} \\
\quad
\begin{matrix}  ... & ... & ... \end{matrix} \quad \quad \quad  
\begin{matrix}  ... & ... & ... \end{matrix} \quad \quad
\begin{matrix}  \ddots \end{matrix}
\end{pmatrix}
$$     

Thanks to this drastic reduction of the IFC coefficients, only 50 to 100 atomic uncorrelated configurations are generally needed to obtain converged properties at the 2nd order (in this example, and in the whole ABINIT package, we only propose examples with a maximum of 20 uncorrelated configurations in order to avoid a too huge amount of data). As previously for the cutoff radius, we can study the convergence of the calculation with respect to the number of uncorrelated atomic configurations. Set [[atdep:nstep_max]] equal to 5 in *tatdep1_1.abi*, replace the root of the ouput file name by "Ntep5" in the *tatdep1_1.files* and launch `atdep`. Do it again for 8 time steps then plot :

      xmgrace -nxy Nstep5omega.dat -nxy Nstep8omega.dat -nxy tatdep1_1omega.dat

You should get the following picture :

![Al phonon spectrum wrt nstep](atdep1_assets/tatdep1_1nstep.png)

In conclusion, a too small number of uncorrelated atomic configurations lead to a large error in the phonon spectrum. Therefore, do not hesitate to pursue the AIMD simulation (in order to accumulate a larger number of configurations) until achieving the convergence of the phonon spectrum.   

!!! note
    Another input variable impact the number of atomic configurations : [[atdep:slice]]. This optional variable select one configuration over [[atdep:slice]], so the calculation will have ([[atdep:nstep_max]]-[[atdep:nstep_min]])/[[atdep:slice]] configurations at all. To test its utility, you can add a line "slice 4" in the section "optional input variable", change the root of the output file name by "slice" and launch `atdep`. The value of this optional variable is now echoed at the begining of the output file and you can find that the "real number of time steps" is now 5 (and no longer 20). Finally, you can plot the phonon spectrum and see the differences with respect to have the 20 configurations (*tatdep1_1omega.dat*) or only the 5 first (*Nstep5omega.dat*). 

####		^^2.3.3 Other important parameters^^

- **k**-points [[cite:Anzellini2020]]
- size of the supercell

## 3. Temperature dependency of a soft mode : U-$\alpha$

This calculation is similar to the one performed in the following article [[cite:Bouchet2015]].

###	3.1 Failure of the QHA

###	3.2 Elastic moduli

###	3.3 Another question?

## 4. Dynamic stabilization due to anharmonic effects : U-$\gamma$

This calculation is similar to the ones performed in the following articles [[cite:Anzellini2020]] and [[cite:Bottin2020]].

###	4.1 Strong anharmonicity

###	4.2 Thermodynamics

###	4.3 Another question?

