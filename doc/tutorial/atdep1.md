---
authors: FB & JB
---

# First tutorial on aTDEP

![INP SIMaP](https://simap.grenoble-inp.fr/medias/photo/thumbnail4-b_1521620102788-jpg)

## The 2$^{nd}$ order **effective** Interatomic Force Constants (IFC)

This tutorial shows how to capture anharmonicities by means of an harmonic **Temperature Dependent Effective Potential** (TDEP) using the ABINIT package.

In practice, this means to obtain the $2^{nd}$ order **effective** IFC. Almost all the dynamic (phonons...), elastic (constants, moduli...) and thermodynamic (entropy, free energy...) desired quantities can be derived therefrom.

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

Let us consider that the potential energy of a crystal can be rewritten using a Taylor expansion around the equilibrium:

$$
U=   U_0 + \sum_{p\ge 1} \frac{1}{p\ !} \sum_{\substack{\alpha_1...\alpha_p \\ i_1...i_p}}\overset{(p)}{\Phi}\vphantom{\Phi}_{i_1...i_p}^{\alpha_1...\alpha_p}\prod_{k=1}^p u_{i_k}^{\alpha_k}
$$

> In this equation, and in the following, the Latin letters in subscripts $i, j, k...$ and the Greek letters in superscripts $\alpha, \beta, \gamma$... will define the atoms and the cartesian directions, respectively.

with U$_0$ the minimum of the potential energy, $u$ the displacement with respect to equilibrium, and $\overset{(p)}{\Phi}\vphantom{\Phi}$ the $p^{th}$ order IFC, respectively. The aim of the TDEP method is to extract the IFC in an effective way. 

Let us assume that :

1. the previous equation is truncated at the $2^{nd}$ order such as : 
$$
U=   U_0 + \frac{1}{2!}\sum_{ij,\alpha\beta} \overset{(2)}{\Theta}\vphantom{\Theta}_{ij}^{\alpha\beta} u_i^\alpha u_j^\beta + 0(u^3)
$$

2. a set of $N_t$ forces $\mathcal{F}_{AIMD}(t)$ and displacements $u_{AIMD}(t)$ is obtained using *ab initio* molecular dynamic (AIMD) simulations, leading to the following system of equations ($\mathcal{F}=-\nabla U$): 

$$
\mathcal{F}_{i,AIMD}^\alpha(t)=  - \sum_{j,\beta} \overset{(2)}{\Theta}\vphantom{\Theta}_{ij}^{\alpha\beta} u_{j,AIMD}^\beta(t)
$$

It is then possible to obtain the second order **effective** IFC $\overset{(2)}{\Theta}\vphantom{\Theta}$ by using a least squares method. This fitting procedure modifies the terms below the truncation by including (in an effective way) the anharmonic contributions coming from the terms above the truncation. Therefore, the IFC are no longer constant and become temperature dependent. That is the reason why we change the notation: in the following, the $\Phi$ will be referred to as the ‘‘true IFC’’ and the $\Theta$ as the ‘‘effective IFC’’.

![Harmonic Potential](https://upload.wikimedia.org/wikipedia/commons/6/60/Potential_approximation.png)

## 2. A simple case : Al-fcc

*Before beginning, you might consider to work in a different subdirectory as for the other tutorials. Why not create Work_atdep1 in \$ABI_TESTS/tutoatdep/Input? You can copy all the input files within.*

```sh
cd $ABI_TESTS/tutoatdep/Input
mkdir Work_atdep1
cd Work_atdep1
cp ../tatdep1_1.files .   # You will need to edit this file.
cp ../tatdep1_1.in .
cp ../HIST1_1.nc .
```

### 2.1 The input files

Let us discuss the meaning of these three files :

* The input file *tatdep1_1.in* :

{% dialog tests/tutoatdep/Input/tatdep1_1.in %}

* The NetCDF file *HIST1_1.nc* ...
* The files file *tatdep1_1.files* lists the file names and root names :

{% dialog tests/tutoatdep/Input/tatdep1_1.files %}

You can now execute `atdep`:

```sh
atdep < tatdep1_1.files > log 2> err &
```

The code should run very quickly.

### 2.2 The output files

The corresponding output files (available in *$ABI_TESTS/tutorial/Refs/) :

* The main output file *tatdep1_1.out* :

{% dialog tests/tutoatdep/Refs/tatdep1_1.out %}

* The phonon frequencies file *tatdep1_1omega.dat* ... You can plot the phonon spectrum. If you use the |xmgrace| tool, launch:

      xmgrace -nxy tatdep1_1omega.dat

  You should get this:

  ![Al phonon spectrum](atdep1_assets/xxx.png)

  At this stage, you have

* The thermodynamic file *tatdep1_1thermo.dat* :

{% dialog tests/tutoatdep/Refs/tatdep1_1thermo.dat %}

### 2.3 Another question?

## 3. Two atoms in the unitcell : Si-d

This calculation is similar to the one performed in the following article [[cite:Bottin2020]].

###	3.1 Convergence with respect to Rcut

###	3.2 Effect of the number of configurations

###	3.3 Another question?



## 4. Temperature dependency of a soft mode : U-$\alpha$

This calculation is similar to the one performed in the following article [[cite:Bouchet2015]].

###	4.1 Convergence with respect to Rcut

###	4.2 Another question?

## 5. Dynamic stabilization due to anharmonic effects : Zr-$\beta$

This calculation is similar to the ones performed in the following articles [[cite:Anzellini2020]] and [[cite:Bottin2020]].

###	4.1 Convergence with respect to Rcut

###	4.2 Another question?

## 6. Description of an alloy using VCA : UMo-$\gamma$

This calculation is similar to the one performed in the following article [[cite:Castellano2020]].

###	4.1 Convergence with respect to Rcut

###	4.2 Another question?
