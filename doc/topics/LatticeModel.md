---
description: How to fit build a lattice model in Multibinit
authors: AM
---

This page gives hints on how to a lattice model in multibinit.

## Introduction

The MULTIBINIT software is using a second-principles approach for lattice dynamics simulations based on atomic potentials fitted on first-principles calculations [[cite:Wojdel2013]].


[[topic:FitProcess | Topic for the fit process]]

[[topic:BoundProcess | Topic for the bound process]]

[[topic:DynamicsMultibinit | Topic to run a dynamics]]  
  
## Related Input Variables

*basic:*

- [[multibinit:dipdip]]  DIPole-DIPole interaction
- [[multibinit:prt_model]]  Effective potential XML output
 
*useful:*

- [[multibinit:coefficients]]  values of the COEFFICIENTS
- [[multibinit:energy_reference]]  Energy of the refences structure
- [[multibinit:ncoeff]]  Number of anharmonic COEFFicients
- [[multibinit:ngqpt]]  Number of Grids points for Q PoinTs
- [[multibinit:nqshft]]  Number of Q SHiFTs
 
*expert:*

- [[multibinit:dipdip_prt]]  DIPole-DIPole PRinT
- [[multibinit:dipdip_range]]  Dipole-Dipole interaction
 

## Selected Input Files

*paral:*

- [[tests/paral/Input/t100.in]]
- [[tests/paral/Input/t101.in]]
- [[tests/paral/Input/t102.in]]
- [[tests/paral/Input/t95.in]]
 
*v8:*

- [[tests/v8/Input/t06.in]]
- [[tests/v8/Input/t07.in]]
- [[tests/v8/Input/t09.in]]
- [[tests/v8/Input/t10.in]]
- [[tests/v8/Input/t11.in]]
- [[tests/v8/Input/t13.in]]
- [[tests/v8/Input/t14.in]]
- [[tests/v8/Input/t15.in]]
 

