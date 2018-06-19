---
description: How to perform a LDA-1/2 calculation
authors: F. Jollet, G. Zerah
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform a LDA-1/2 calculation with the ABINIT package.

## Introduction

This feature is available only in PAW. The LDA-1/2 framework is described in
[[cite:Ferreira2008]]. In ABINIT, the LDA-1/2 approximation needs to use PAW
data files that contain a small potential vminushalf that is added to the
local potential inside the code. The LDA-1/2 approach is based on Slater's
half occupation technique to restore the band gap of semi-conductors and
insulators.


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

