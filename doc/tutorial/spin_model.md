---
authors: xuhe
---

# Spin model in Multibinit

## Build a spin model and run spin dynamics in Multibinit

This lesson aims at showing how to build a spin model and run spin dynamics.

**Before beginning, it is very important to read the reference [[cite:Wojdel2013]].**

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

$$E^{exc} =- \sum_{i\neq j} J_{ij} \vec{S_i}\cdot{\vec{S_j}}$$ .

A few conventions used in Multibinit should be noted:

- all the $\vec{S}$ are normalized to 1. 
- Both $J_{ij}$ and $J_{ji}$ are in the Hamiltonian.
- There is a minus sign, which means that the interaction is ferromagnetic for $J_{ij} >0$.

As the sites $i$ and $j$ defined in the model are in a finite cell (often the primitive cell), there are interactions between sites inside (eg. site $i$) and outside the cell (site $j'$). For a site at the position $\vec{r_k}=\vec{r_{j'}}+\vec{R}$, here $r_j$ is the position of $j$ inside the cell, and $\vec{R}$ is an integer multiples of the lattice vectors. For site inside the cell,its $\vec{R}$ vector is (0,0,0).  Therefore, $J$ between site $i$ and $k$ are labeled with $i$, $j$ and $R$, denoted as $J_{ij}(\vec{0}, \vec{R})$. Note that due to the translation symmetry, $J_{ij}(\vec{R}, 0)= J_{ij}(\vec{0}, -\vec{R})$, thus we can alway use $\vec{0}$ for the first index, and $J_{ij}(\vec{0}, \vec{R})$ can be simplified as $J_{ij}(\vec{R})$. 

## 2. Format of xml file for spin mode

## 3. Build spin model files

## 4. Run spin dynamics

Now that we have the spin model file, we can run spin dynamics

## 5. Postprocessing









