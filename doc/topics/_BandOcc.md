---
description: How to to specify bands and occupation numbers, for metals or insulators
authors: FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to to specify bands and occupation numbers, 
for metals or insulators with the ABINIT package.

## Introduction

  
Metallic as well as insulating systems can be treated, depending on the value
of [[occopt]]. The default value of [[occopt]] corresponds to an insulator (or
finite molecule): the number of bands (or states for a molecule) is deduced
from the number of electrons brought by each pseudopotential ion, and then all
the bands are occupied (by two electrons in case of a non-spin-polarized
system, or by 1 electron in the case of a spin-polarized system), and a small
number of empty bands are added, e.g. to obtain the band gap.

For a metallic system, use a value of [[occopt]] between 3 and 7. ABINIT will
compute a default number of bands, including some nearly unoccupied ones, and
find the occupation numbers. The different values of [[occopt]] correspond to
different smearing schemes (smearning defined by [[tsmear]] for defining the
occupation numbers, e.g. Fermi broadening, the Gaussian broadening, the
Gaussian-Hermite broadening, as well as the modifications proposed by Marzari.
Finite temperatures can also be treated thanks to a smearing scheme
(Verstraete scheme) using [[tphysel]].

It is possible to define manually the number of bands (input variable
[[nband]]) as well as the occupation numbers (input variable [[occ]]). This
might be useful to perform a Δ-SCF calculation for excited states.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:base1|The tutorial 1]] deals with the H2 molecule: get the total energy, the electronic energies, the charge density, the bond length, the atomisation energy 
* [[tutorial:base2|The tutorial 2]] deals again with the H2 molecule: convergence studies, LDA versus GGA 
* [[tutorial:base3|The tutorial 3]] deals with crystalline silicon (an insulator): the definition of a k-point grid, the smearing of the cut-off energy, the computation of a band structure, and again, convergence studies ...
* [[tutorial:base4|The tutorial 4]] deals with crystalline aluminum (a metal), and its surface: occupation numbers, smearing the Fermi-Dirac distribution, the surface energy, and again, convergence studies ...

