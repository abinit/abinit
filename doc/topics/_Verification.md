---
description: How to become convinced that results are numerically correct
authors: XG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to become convinced that results are numerically correct with the ABINIT package.

## Introduction

One can never be sure that an implementation of a complicated formalism is
numerically correct. However, several studies have helped to gain confidence
in selected properties computed within ABINIT, by making detailed comparisons
with independently developed software applications.

The ground state total energy, geometry relaxation, phonon frequencies,
electron-phonon matrix elements, temperature dependence of the electronic gap
has been cross checked with Quantum Espresso and YAMBO. See [[cite:Ponce2014]].

It is possible to choose the same convention for the definition of the average
eletrostatic potential, than Quantum Espresso, allowing verification of
results for charged systems. See [[cite:Bruneval2014]].

See the "Validation" section in [[topic:PseudosPAW]].



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

