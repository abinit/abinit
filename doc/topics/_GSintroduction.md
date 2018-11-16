---
description: How to build an input file for a ground state calculation
authors: FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to build an input file for a ground state calculation with the ABINIT package.

## Introduction

The computation of the ground state energy of an assembly of nuclei and
electrons placed in a repeated cell can be done using (1) plane waves and
norm-conserving pseudopotentials, or, (2) so-called "Projector-Augmented
Waves" (PAW method), with appropriate pseudoatomic data, or (3) wavelets. The
wavelet framework is described [[topic:Wavelets|here]].  

In the plane wave framework, the program admits many different types of
pseudopotentials. There are several complete sets of norm-conserving
pseudopotentials available for most elements of the periodic table. 

The recommended tables (GGA-PBE, GGA-PBEsol and LDA) come from the |pseudodojo| project 
with ONCVPSP pseudopotentials ([[cite:Hamann2013]]) 
both in scalar-relativistic format and fully-relativistic version with spin-orbit coupling. 
For PAW calculation,the recommended one (GGA-PBE and LDA) is the JTH
table in the PAW XML format (([cite:Jollet2014]]). 

The choice between norm-conserving
pseudopotentials or PAW is deduced automatically by the choice of the
pseudopotential in the "files" file. An input file must specify the following
items: 

* the [[topic:crystal|crystalline structure and symmetries]]
* the set of [[topic:k-points|k-points]] used
* the [[topic:xc|exchange and correlation functional]]
* the convergence settings
* possibly [[topic:PAW|PAW]] special settings
* possibly, input variables for [[topic:spinpolarisation|spin-polarized systems and spin orbit coupling]] calculations.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:base1|The tutorial 1]] deals with the H2 molecule: get the total energy, the electronic energies, the charge density, the bond length, the atomisation energy 
* [[tutorial:base2|The tutorial 2]] deals again with the H2 molecule: convergence studies, LDA versus GGA 
* [[tutorial:base3|The tutorial 3]] deals with crystalline silicon (an insulator): the definition of a k-point grid, the smearing of the cut-off energy, the computation of a band structure, and again, convergence studies ...
* [[tutorial:base4|The tutorial 4]]] deals with crystalline aluminum (a metal), and its surface: occupation numbers, smearing the Fermi-Dirac distribution, the surface energy, and again, convergence studies ...

