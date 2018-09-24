---
authors: xuhe
---

# Spin model in Multibinit

## Build a spin model and run spin dynamics in Multibinit

This lesson aims at showing how to build a spin model and run spin dynamics.

**Before beginning, it is important to know the background of spin dynamics, which can be found in literatures (e.g.[[cite: Evans2014]] )**

With this lesson, you will learn to:

  * Generate the XML for the model 
  * Run a spin dynamics within MULTIBINIT

The TB2J python package, which can be used to generate spin model, can be found on the abinit gitlab website at https://gitlab.abinit.org/xuhe/TB2J. This package will be included in the Abinit package in the future. 



[TUTORIAL_README]

## 1 Heisenberg Model formalism

The spin model is defined as a classical Heisenberg model. In the current version of Multibinit, the interactions can be considered are: the exchange interaction, single ion anisotropy (SIA),  Dzyaloshinski-Moriya (DM) interaction, external magnetic field. 

$$E = E^{exc}+E^{SIA} + E^{DM}+E^{ext}$$

Here we outline the Heisenberg model formalism as implemented in Multibinit.

 The energy of exchange interaction can be written as

$$E^{exc} =- \sum_{i\neq j} J_{ij} \vec{S_i}\cdot{\vec{S_j}}$$ 

A few conventions used in Multibinit should be noted:

- all the $\vec{S}$ are normalized to 1. 
- Both $J_{ij}$ and $J_{ji}$ are in the Hamiltonian.
- There is a minus sign, which means that the interaction is ferromagnetic for $J_{ij} >0$.

As the sites $i$ and $j$ defined in the model are in a finite cell (often the primitive cell), there are interactions between sites inside (eg. site $i$) and outside the cell (site $j'$). For a site at the position $\vec{r_k}=\vec{r_{j'}}+\vec{R}$, here $r_j$ is the position of $j$ inside the cell, and $\vec{R}$ is an integer multiples of the lattice vectors. For site inside the cell,its $\vec{R}$ vector is (0,0,0).  Therefore, $J$ between site $i$ and $k$ are labeled with $i$, $j$ and $R$, denoted as $J_{ij}(\vec{0}, \vec{R})$. Note that due to the translation symmetry, $J_{ij}(\vec{R}, 0)= J_{ij}(\vec{0}, -\vec{R})$, thus we can alway use $\vec{0}$ for the first index, and $J_{ij}(\vec{0}, \vec{R})$ can be simplified as $J_{ij}(\vec{R})$. 

The SIA term can be written as 

$$E^{SIA}=-k_u \sum_i (\vec{S_i}\cdot \vec{e})^2$$

, where $k_u$ and $\vec{e}$ are the amplitude and direction of the single ion anisotropy.

The DM term can be written as

$$ E^{DM} = \sum_{i\neq j} \vec{D}_{ij}\cdot \vec{S_i}\times{\vec{S_j}}$$

, where $\vec{D_{ij}}$ is the amplitude of DM interaction.

The external magnetic field term can be written as

$$E^{ext}=- \sum_i   m_{i}   \vec{S_i}\cdot \vec{H}$$

, where $m_i$ is the magnetic moment of site $i$, and $\vec{H}$ is the magnetic field.

## 2. Format of xml file for spin mode (TODO: remove?)

## 3. Build spin model file

One way to calculate the Heisenberg model parameters is to use the spin force theorem (see [[cite:Liechtenstein1983]], [[cite:Katsnelson2000]]), which take rotations of localized spin as perturbations. In Abinit, the hamiltonian uses plane-waves as basis set, thus the localized spin is not directly accessible. We can construct localized Wannier functions as basis set and rewrite the hamiltonian. Then the exchange parameters can be calculated from this hamiltonian ( [[cite:Korotin2015]] ). 

For building the Wannier function Hamiltonian from Abinit, see tutorial [wannier90](wannier90). Other DFT codes interfaced with [Wannier90](http://www.wannier.org) can also be used. Then [TB2J](https://gitlab.abinit.org/xuhe/TB2J) package can be used to get the Heisenberg model parameters and generate the input model for Multibinit. 



## 4. Run spin dynamics

Now that we have the spin model xml file, we can run spin dynamics from multibinit.  An example input files can be found at ~abinit/tests/tutomultibinit/Input/tmulti_spin1.* .  Three files are there: 

* "tmulti_spin1.files" is the files file, which gives the names of the input and output files for  Multibinit.

* "tmulti_spin1.xml" is the file containing the Heisenberg model parameters.

* "tmulti_spin1.in" is the main input file, where you can put the parameters for the spin dynamics simulation.

You can copy these three files into a directory. 

## 5. Postprocessing










