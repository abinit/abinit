---
authors: FJ
description: GSintroduction Abinit topic
---

This page gives hints on how to build an input file for a ground state calculation with the ABINIT package.

## Introduction

The computation of the ground state energy of an assembly of nuclei and
electrons placed in a repeated cell can be done using (1) plane waves and
norm-conserving pseudopotentials, or, (2) so-called "Projector-Augmented
Waves" (PAW method), with appropriate pseudoatomic data, or (3) wavelets. The
wavelet framework is described [[topic:Wavelets|here]].  
In the plane wave framework, the program admits many different types of
pseudopotentials. There are several complete sets of norm-conserving
pseudopotentials available for most elements of the periodic table. The
recommended one (GGA) comes from the ONCVPSP generator (with spin-orbit
coupling). For PAW calculation,the recommended one (GGA and LDA) is the JTH
table in the PAW XML format. The choice between norm-conserving
pseudopotentials or PAW is deduced automatically by the choice of the
pseudopotential in the "files" file. An input file must specify the following
items: the [[topic:crystal|crystalline structure and symmetries]]. the set of
[[topic:k-points|k-points]] used. the [[topic:xc|exchange and correlation
functional]]. the convergency settings. possibly [[topic:PAW|PAW]] special
settings. possibly, input variables for [[topic:spinpolarisation|spin-
polarized systems and spin orbit coupling]] calculations.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[lesson:base1|The lesson 1]] deals with the H2 molecule: get the total energy, the electronic energies, the charge density, the bond length, the atomisation energy 
* [[lesson:base2|The lesson 2]] deals again with the H2 molecule: convergence studies, LDA versus GGA 
* [[lesson:base3|The lesson 3]] deals with crystalline silicon (an insulator): the definition of a k-point grid, the smearing of the cut-off energy, the computation of a band structure, and again, convergence studies ...
* [[lesson:base4|The lesson 4]]] deals with crystalline aluminum (a metal), and its surface: occupation numbers, smearing the Fermi-Dirac distribution, the surface energy, and again, convergence studies ...

