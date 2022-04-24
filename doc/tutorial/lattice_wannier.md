---
authors: HeXu
---

# Tutorials on constructing Lattice Wannier functions (LWF's). 

This tutorial introduces the methods for constructing lattice Wannier funtions (LWF's). You'll learn how to use the SCDM-k (selected columns of the density matrix in k-space) method for constructing the LWF's. And how to plot the band structure of the LWF's and compare it with the full phonon band structures. 

The LWF is a close analogy to the electron Wannier functions.  Each LWF is formed by a group of atomic displacement within a certain range. Together they form a localized basis set for the atomic distortions, which could span the (sub)space for the atomic vibrations (phonons). One typical use case is to build effective Hamiltonian of atomic displacement. In many phenomenons, only a few distortion modes are important. For example,  in a structural phase transition the soft modes are often related, whereas the phonons with much higher frequency are probably less relevant. Thus if the LWF's can well represent the relevant phonon modes, it could reduce the number of degrees of freedom in the study of such phenomenons. 

In anaddb, LWF's could be constructed with the selected columns of the density matrix (SCDM-k) algorithm or the projected Wannier function (PWF) algorithm. Both methods were initially proposed for electron Wannier functions and were adapted for the LWF's. In this tutoral, we'll first learn about the algorithms and the practical usage of the methods. 


## SCDM-k method

The SCDM-k method is based on the observation that the columns of the density matrices (DM) are often localized, given that the DM is expressed in local basis set. In the case of a periodic crystal, the columns density matrices for all wave vector is localized and  are continuous in reciprocal space. Therefore by transform these columns to the real space, one could get localized functions. Then if the columns of the DM can approximately span the space of the lattice distortions of interest, they can effectively form the basis for the lattice distortions. 




## Projected-Wannier function method


## Practical guide on to build LWF's


### Building LWF's with the SCDM-k method


### Building LWF's with the PWF method


