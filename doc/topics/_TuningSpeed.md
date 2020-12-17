---
description: How to tune the speed and memory usage
authors: XG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to tune the speed and memory usage with the ABINIT package.

## Introduction

The major factors governing the speed of an ABINIT run, for a given physical system, are:

  * the size of the planewave basis set, see [[topic:Planewaves]];
  * the size of the wavevector grid sampling, see [[topic:k-points]];
  * and the parallelism, see [[topic:parallelism]].

For the two first factors, there is a trade-off between CPU time and precision
of the computation, while for the third factor, there is some limit on the
maximal speed-up that can be achieved (and also, the resources must be available.

Beyond these major factors, there is still room for some adjustment. The
needed planewave basis set will depend on the pseudopotential (or PAW atomic
dataset) that is used. Some might be softer than others and need a smaller
planewave basis set. They might possibly be less accurate as well ...

If one is only interested in ground-state properties and forces, one might
also get some speed up by using a real-space representation of density and
potential on a real space FFT grid that does not allow their fine details to
be taken into account (actually, filtering such quantities). This is achieved
by lowering [[boxcutmin]] below its theoretically needed value of 2.0.

The choice of the FFT algorithm implementation, see [[fftalg]] might also lead
to significant speed up, on specific machines.

For specific k-points, time-reversal symmetry can be used to represent the
wavefunctions with their real part, instead of both their real and complex parts. 
This allows halving the memory needs, as well as the CPU time. 
See [[istwfk]].

Other input variables related to tuning the speed or the memory usage are for expert users only.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

