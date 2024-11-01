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
table in the PAW XML format ([[cite:Jollet2014]]). 

The choice between norm-conserving
pseudopotentials or PAW is deduced automatically by the choice of the
pseudopotential in the "files" file. An input file must specify the following items: 

* the [[topic:crystal|crystalline structure and symmetries]]
* the set of [[topic:k-points|k-points]] used
* the [[topic:xc|exchange and correlation functional]]
* the convergence settings
* possibly [[topic:PAW|PAW]] special settings
* possibly, input variables for [[topic:spinpolarisation|spin-polarized systems and spin orbit coupling]] calculations.


An example of a minimal input file to calculate the ground state of crystalline aluminium is given here:


```
# Crystalline aluminum. Calculation of the total energy
# at fixed number of k points and broadening.

#Definition of occupation numbers
occopt 4
tsmear 0.05

#Definition of the unit cell
acell 3*7.60           # This is equivalent to   7.60 7.60 7.60
rprim  0.0  0.5  0.5   # FCC primitive vectors (to be scaled by acell)
       0.5  0.0  0.5
       0.5  0.5  0.0

#Definition of the atom types
ntypat 1          # There is only one type of atom
znucl 13          # The keyword "znucl" refers to the atomic number of the
                  # possible type(s) of atom. The pseudopotential(s)
                  # mentioned in the "files" file must correspond
                  # to the type(s) of atom. Here, the only type is Aluminum

#Definition of the atoms
natom 1           # There is only one atom per cell
typat 1           # This atom is of type 1, that is, Aluminum
xred  0.0  0.0  0.0 # This keyword indicate that the location of the atoms
                    # will follow, one triplet of number for each atom
                    # Triplet giving the REDUCED coordinate of atom 1.

#Definition of the planewave basis set
ecut  6.0         # Maximal kinetic energy cut-off, in Hartree
pawecutdg  10.0   #Maximal kinetic energy cut-off, in Hartree for the fine grid in case of PAW calculation

#Definition of the k-point grid
ngkpt 2 2 2       # This is a 2x2x2 FCC grid, based on the primitive vectors
chksymbreak 0
#Definition of the SCF procedure
nstep 10          # Maximal number of SCF cycles
toldfe 1.0d-6     # Will stop when, twice in a row, the difference
                  # between two consecutive evaluations of total energy
                  # differ by less than toldfe (in Hartree)
                  # This value is way too large for most realistic studies of materials       
```

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:base1|The tutorial 1]] deals with the H2 molecule: get the total energy, the electronic energies, the charge density, the bond length, the atomisation energy 
* [[tutorial:base2|The tutorial 2]] deals again with the H2 molecule: convergence studies, LDA versus GGA 
* [[tutorial:base3|The tutorial 3]] deals with crystalline silicon (an insulator): the definition of a k-point grid, the smearing of the cut-off energy, the computation of a band structure, and again, convergence studies ...
* [[tutorial:base4|The tutorial 4]]] deals with crystalline aluminum (a metal), and its surface: occupation numbers, smearing the Fermi-Dirac distribution, the surface energy, and again, convergence studies ...

