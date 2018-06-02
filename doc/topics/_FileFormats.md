---
description: How to manage file formats, and the interfacing with other applications outside of the ABINIT organisation
authors: XG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to manage file formats, and the interfacing with other applications outside of
the ABINIT organisation with the ABINIT package.

## Introduction

For a long time, ABINIT has produced files in either text format or in Fortran
binary format. This has the drawback of being difficult to read by other
codes, and also being difficult to maintain in the long run, as any change
must be propagated to the reading routines.

Following the advent of the [Nanoquanta/ETSF file format](http://www.etsf.eu/fileformats), 
ABINIT has gradually shifted toward
the use of file formats that are addressable by content, especially NetCDF and XML.

NetCDF is used to file that typically store a large amount of data, like
wavefunctions, density, potentials, etc . One should even shift to HDF5 in due
time for these files.

XML is used for smaller files. E.g. PSML pseudopotentials can be used, these
being common with the SIESTA code. YAML files have allso recently appeared,
for the ABINIT documentation.


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

