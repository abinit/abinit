# GPU support in Abinit

  IMPORTANT: GPU support is currently highly *EXPERIMENTAL* and should
             be used by experienced developers only. If you encounter
             any problem, please report it to Yann Pouillon before doing
             anything else.

## GPU-related parameters

GPU support is activated by the --enable-gpu option of configure.
Another option of importance is the --with-gpu-flavor one, which selects
the kind of GPU support that will be activated. A convenience option,
codename --with-gpu-prefix, is also provided, in order to set
automatically all relevant parameters whenever possible. A few other
options are available as well, mainly for fine-tuning of the build
parameters and testing purposes.

Full descriptions of all these options can be found in the
~abinit/doc/build/config-template.ac9 file. Do not hesitate to ask
questions on https://forum.abinit.org/.

In addition, the permitted GPU-related preprocessiong options are:

  * HAVE_GPU        : generic use;
  * HAVE_GPU_SERIAL : serial GPU support;
  * HAVE_GPU_MPI    : MPI-aware GPU support.


## Cuda support

At present it is possible to ask for single- or double-precision Cuda
support. The configure script will check that the Cuda libraries are
properly working, but however not whether double-precision is actually
supported by your version of Cuda (this might be added in the future).

All calls to Cuda routines should be carefully embedded within
'#if defined HAVE_GPU_CUDA ... #else ... #endif' preprocessing blocks.
When a feature does require Cuda and will not work without it, the
corresponding '#else' part should display an error and cause Abinit to
abort.

The permitted Cuda-related preprocessing options are :

  * HAVE_GPU_CUDA    : generic use;
  * HAVE_GPU_CUDA_SP : single-precision calculations;
  * HAVE_GPU_CUDA_DP : double-precision calculations.

All high-level routines directly accessing Cuda features have to be put
in ~abinit/src/52_manage_cuda/, and low-level ones in
~abinit/shared/common/src/17_gpu_toolbox/. 

All files belonging to nVidia must *not* be distributed with Abinit.
Please discuss with Yann Pouillon if you need them inside the Abinit
source tree during the build.

In any case, all Cuda-related developments should be done in good
coordination with:

  * Marc Torrent
  * Yann Pouillon

## Cuda version

To take advantage of the multiple FFT in cuda (FFT in batch), ABINIT
have to be compiled with a Cuda version>=3.0.
It is possible to build with previous versions (>2.1 tested) but you
make some changes. 
cuda implementation support devices with capabilty (revision)>1.0

## Magma support

The MAGMA project aims to develop a dense linear algebra library similar
to LAPACK but for heterogeneous/hybrid architectures, starting
with current "Multicore+GPU" systems.
It is recommended to take advantage of MAGMA when using ABINIT with
Cuda.
Magma is not distributed within ABINIT package; it has to be preliminary
installed. To activate MAGMA support during building process, use
--wih-linalg-flavor="...+magma" at configure level.

## OpenCL support

OpenCL support is currently under discussion. More info will come once
decisions have been taken.


## S_GPU support

The S_GPU library provides higher performance and better load balancing
when each GPU of a hybrid computer is shared by several processes, e.g.
MPI tasks.

It will be supported in Abinit in the future, from its version 2.

See http://ligforge.imag.fr/projects/sgpu/ for details.
