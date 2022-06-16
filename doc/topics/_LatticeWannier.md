---
description: How to build and run Lattice Wannier function Models 
authors: HeXu
---


This page gives hints on how to build and run lattice Wannier function (LWF's) models, including:
* How to build LWF's from phonon calculations.
<!-- * How to run dynamics of LWF's with MULTIBINIT. -->

## Introduction
The lattice Wannier functions (LWF's) are localized functions of atomic displacements in crystal structures which are the analogue of electronic Wannier functions. They could be used as basis set of effective Hamiltonians for the description of lattice distortions. 

In anaddb, we can use the "selected cloumns of the density matrix" (SCDM-k) algorithm for the construction of the LWF's, and the Hamiltonian of them within the harmonic approximation. Starting from that, we can build structures with various amplitudes of LWF's and compute the energies of them. Then the anharmonic interaction between the LWF's can be parameterized and added to the Hamiltonian. We can then run dynamics for the LWF's with the Hamiltonian. Since the LWF's can be constructed to give good description for a few phonon branches relevant to the phenomenons under study with the number of degrees of freedom largely reduced, it allows for simplification of the model without dropping the essential physics. The difficulty to fit of the anharmonic interaction, which grows with the number of degrees of freedoms, is thus much more easier. It also reduces the requirement for computational power when doing the dynamics. 


## Tutorials
[[tutorial:lattice_wannier][How to build Lattice Wannier functions]] : The tutorial on how to build LWF's with anaddb from phonon calculations.

  
## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

