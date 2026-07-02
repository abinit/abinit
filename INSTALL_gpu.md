# GPU support in ABINIT

## Introduction

This page explains how to compile ABINIT with GPU support enabled and provides a brief overview of the current state of GPU support in ABINIT.

Current GPU support primarily relies on OpenMP GPU offloading and is actively developed to support more features of ABINIT.
Legacy GPU implementations also remain in the code and can still be enabled, but their scope is limited to ground-state calculations on NVIDIA GPUs.

## State of GPU support

Since ABINIT v10.0 (April 2024), a new GPU implementation of ABINIT has been made available. It is based on the OpenMP GPU offload paradigm and supports both NVIDIA and AMD GPUs.

GPU vendor libraries are also used:

 - cuBLAS/cuFFT/cuSOLVER from CUDA SDK for NVIDIA GPUs,
 - hipBLAS/hipFFT/hipSOLVER from ROCm SDK for AMD GPUs,

As of version 10.8, the only FULLY validated compilers for enabling GPU support in ABINIT through OpenMP GPU offload are:

- [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk), or simply NVHPC, exclusively for NVIDIA GPUs
- [Cray Compiler Environment](https://cpe.ext.hpe.com/docs/latest/cce/index.html), or simply Cray or CCE, validated for AMD GPUs

Support for other compilers such as GCC or LLVM is currently being explored but isn't production-ready yet.


## How to enable GPU during compilation

ABINIT can be built using two build systems, Autotools and CMake, the latter being under developement.

The three following settings are required to enable GPU support:

- GPU architecture identifier (example : `80` for NVIDIA A100, `gfx90a` for AMD MI250)
- CUDA or ROCm usage and installation root path (path auto-detected in CMake)
- OpenMP offload activation (auto-detected in CMake)


### Using Autotools

The following flags are REQUIRED for enabling GPU support:

 - `--enable-openmp-offload`
 - `--with-cuda=<path>` or `--with-rocm=<path>` for CUDA/ROCm root installation path
 - `GPU_ARCH` variable must be set to the architecture identifier of the targeted GPU

The following flags may be used for additional tuning:

 - `--enable-mpi-gpu-aware` enables passing GPU buffers directly to selected MPI calls, to leverage GPU-Direct capabilities
 - `--enable-gpu-nvidia-unified-memory` is specific to NVHPC compiler and enables unified-memory optimizations

Finally, the flag `--with-gpu-markers=yes` may be passed to enable NVTX/rocTX markers calls in ABINIT code for profiling purposes (using tools such as [NVIDIA NSight Systems](https://developer.nvidia.com/nsight-systems) or [AMD ROCprofiler](https://rocm.docs.amd.com/projects/rocprofiler-sdk/en/latest/index.html)).


The traditional way to use Autotools with ABINIT is usually by storing flags inside a `.ac9` configuration file.

Three example configuration files are provided below, covering both NVIDIA and AMD GPUs.
Adapt them as needed for your hardware and software environment.

- [**Architecture 1**](build/GPU_Nvidia_A100+EPYC7763.ac9.md)
   CPU: AMD Milan EPYC 7763
   GPU: Nvidia A100 SXM4 80 GB
   Compilers: Nvidia HPC compilers
   Libraries: MKL, CUDA

- [**Architecture 2**](build/GPU_Nvidia_GH200.ac9.md)
   CPU: NVIDIA Grace ARM Neoverse-V2
   GPU: Nvidia H200 90 GB
   Compilers: Nvidia HPC compilers
   Libraries: NVPL, CUDA

- [**Architecture 3**](build/GPU_AMD_InstinctMI250X+EPYC7453.ac9.md)
   CPU: AMD Trento EPYC 7A53
   GPU: AMD Instinct MI250X 2*64 GB
   Compilers: Cray Compiler Environment (CCE)
   Libraries: FFTW3, libSCI, HIP+ROCm


### Using CMake

The following flags are REQUIRED for enabling GPU support:

| GPU vendor | NVIDIA | AMD |
|------------|--------|-----|
| Enable GPU |`-DABINIT_ENABLE_GPU_CUDA=ON` | `-DABINIT_ENABLE_GPU_HIP=ON` |
| GPU architecture code | `-DCMAKE_CUDA_ARCHITECTURES=XYZ` | `-DCMAKE_HIP_ARCHITECTURES=gfxXYZ` |


OpenMP offload support should be detected and enabled automatically.

The following flags may be used for additional tuning:

 - `-DABINIT_ENABLE_GPU_AWARE_MPI=ON` enables passing GPU buffers directly to selected MPI calls, to leverage GPU-Direct capabilities
 - `-DABINIT_ENABLE_NVIDIA_UNIFIED_MEM=ON` is specific to NVHPC compiler and tells to enable unified memory tuning

Finally, the flag `-DABINIT_ENABLE_GPU_MARKERS=ON` may be passed for enabling NVTX/rocTX markers calls in ABINIT code for profiling purposes (using tools such as [NVIDIA NSight Systems](https://developer.nvidia.com/nsight-systems) or [AMD ROCprofiler](https://rocm.docs.amd.com/projects/rocprofiler-sdk/en/latest/index.html)).


## Running on GPUs

In order to effectively use GPU at runtime, a specific parameter must be set in the ABINIT input file, [[gpu_option]].

Simply add the following line in your file :
```
 gpu_option "GPU_OPENMP"
```
or the less verbose:
```
 gpu_option 2
```

For Ground-State calculations, it is best to use GPU with Chebyshev Filtering as the wavefunction optimisation algorithm :
```
wfoptalg 111

# second choice, less efficient, and requires tuning with nblock_lobpcg
# wfoptalg 114
```

In addition to ground-state calculations, GPU acceleration is also available for:

- DFPT (response function), with [[rfphon]], [[rfelfd]], [[rfstrs]] or [[rfddk]] perturbations
- DMFT (many-body theory), for accelerating Green's function computations, and parts of Hubbard-One and CT-QMC solvers



Other values for [[gpu_option]] serve to access previous GPU implementations and are restricted to GS calculation.


## Legacy GPU implementations

Previous GPU implementations have been deprecated but kept in the code.
They can still be used for testing purposes but are no longer supported.

The first implementation from 2013 was based on CUDA and MAGMA and is now fully obsolete.
It is restricted to NVIDIA GPUs and only requires flags related to CUDA support.
However, it remains available for testing purposes only, provided that [[wfoptalg]] is set to 14.

A subsequent GPU implementation using [Kokkos](https://github.com/kokkos/kokkos) and [YAKL](https://github.com/mrnorman/YAKL) and partially leveraging previous CUDA implementation was also developed.

This Kokkos implementation supports the same `Chebyshev filtering` algorithm ([[wfoptalg]] set to `111`) as current OpenMP implementation and should offer similar performance on ground-state calculations.
However, its development was discontinued, so it does not support additional use cases.

To enable Kokkos support in the Autotools build system, one must provide the following flags:

- `--with-kokkos=<path>` for Kokkos installation directory
- `--with-yakl=<path>` for YAKL installation directory
- `--with-cuda=<path>` for CUDA installation directory

CMake requires similar flags but should auto-detect installation directories:

- `-DABINIT_KOKKOS_WANTED=ON` to enable Kokkos
- `-DABINIT_YAKL_WANTED=ON` to enable YAKL
- `-DABINIT_ENABLE_GPU_CUDA` to enable CUDA

Both Kokkos and YAKL can be built on the fly if not available:

- `-DABINIT_KOKKOS_BUILD=ON` to build Kokkos (must be installed otherwise)
- `-DABINIT_YAKL_BUILD=ON` to build YAKL (must be installed otherwise)
