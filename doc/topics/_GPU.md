---
description: How to set parameters for a parallel calculation with GPU
authors: MT,FJ,MS
---
<!--- This is the source file for this topics. Can be edited. -->

This page provides guidance on how to configure parallel calculations using GPU acceleration in ABINIT.

It assumes that you already know the basics of [parallelism](parallelism.md) in ABINIT.

GPU support must have been enabled when compiling ABINIT (details [here](../INSTALL_gpu.md)).

> IMPORTANT:
Unless stated otherwise, this page only covers the OpenMP GPU Offload implementation ([[gpu_option]]=GPU_OPENMP).


## Supported use cases

The following use cases are supported when using [[gpu_option]]=`GPU_OPENMP`:

- Basic Ground-State calculations, with ChebFI or LOBPCG as SCF solvers ([[wfoptalg]]={111,114}, described in [[cite:Lygatsika2026]]), including LDA, GGA, meta-GGA and hybrid functionals

- Response Function calculations (DFPT), with respect to atom displacement, electric fields and strain perturbation ([[rfphon]], [[rfelfd]] and [[rfstrs]], along [[rfddk]])

- DMFT calculations, with acceleration of Green's function calculation, and partial optimisation of Hubbard-One and CT-QMC solvers ([[dmft_solv]]={2,5})

For other values of [[gpu_option]], GPU acceleration is limited to basic ground-state calculations, using ChebFI ([[gpu_option]]=`GPU_KOKKOS`) and the legacy LOBPCG implementation ([[gpu_option]]=`GPU_LEGACY`), respectively.


## Basic usage

Running ABINIT in parallel using GPU requires setting a specific input parameter : [[gpu_option]].
As stated earlier, preferred value is `GPU_OPENMP` (other choices aren't actively maintained).

You can then run ABINIT with MPI as usual using `mpirun`. In most cases, using as many MPI tasks as GPUs generally provides good performance, **each MPI task being bound to a single GPU**.

For example, if you run on a single node having 4 GPUs, you would use 4 MPI tasks:

```
mpirun -n 4 abinit run.abi > log 2> err
```

You can also combine GPU acceleration with OpenMP thread parallelism:

```
export OMP_NUM_THREADS=8
mpirun -n 4 abinit run.abi > log 2> err
```

If [[paral_kgb]] is enabled and you have one K-point, [[npband]] should be only used to distribute MPI tasks.

As for SCF solver, Chebyshev Filtering algorithm ([[wfoptalg]]=111) should be preferred in most cases, while LOBPCG ([[wfoptalg]]=114) serves as a fallback.

In some cases, such as ground-state calculations with many k-points, it may be beneficial to run multiple MPI tasks per GPU, and distributing using [[np_spkpt]].

## Tuning for GPU

Some use case or computing environments may require additional tuning to achieve optimal performance.

### Tuning options at compilation

Some compilation flags should be enabled to take advantage of some hardware features of your GPU nodes:

- GPU-aware MPI: allows the MPI library to access GPU buffers directly, avoiding additional data transfers between CPU and GPU
- Unified memory (NVHPC): this flag is specific to NVHPC and is only useful on NVIDIA so-called Superchips such as GH200 or GB200 that feature unified memory

Usually, it is always a good idea to enable GPU-aware MPI directives.
As for Unified memory, it may improve performance in some cases but is currently known to be penalizing on huge systems.

See the installation procedure for [enabling GPU in ABINIT](../INSTALL_gpu.md#how-to-enable-gpu-during-compilation) where those flags are documented.


### MPI tasks placement/binding on GPUs

Like CPUs, GPUs are assigned numerical identifiers and are addressed using those identifiers.. When using ABINIT, each MPI task will select a GPU to work with for the entirety of the computation. In a simple configuration, each MPI task selects the GPU whose ID satisfies `modulo(mpi_rank/gpu_count)`.

WWhile explicit MPI task binding is often beneficial on CPU-only workloads, it becomes non-trivial on some GPU configurations, such as AMD Instinct MI250X nodes where GPU numbering does not match CPU topology.

Running multiple MPI tasks per GPU may also require a dedicated binding strategy to achieve good performance.

Please refer to your cluster documentation for advice on binding on GPU nodes.


### Running multiple MPI tasks per GPU on NVIDIA GPUs

Some use cases may benefit from running using many MPI tasks per GPU.
This approach is most useful for relatively small calculations where a single MPI task does not generate enough work to fully occupy the GPU.

On NVIDIA GPUs, a feature called [MPS (Multi-Process Service)](https://docs.nvidia.com/deploy/mps/latest/index.html) exists to accelerate this setting. MPS allows a GPU to be shared efficiently among multiple MPI tasks. A common approach is to allocate a fraction of the GPU resources to each task by specifying a percentage of the available GPU capacity.

For example, if you run with 4 MPI per GPU, you would specify a GPU share of 25 % to NVIDIA MPS, so that each MPI task would use its dedicated share of resources on the GPU without competing with others.

Because MPS deployment is system-dependent, you may need to refer to your cluster documentation for enabling this tuning feature.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:basepar|An introduction on ABINIT in Parallel]] should be read before going to the next tutorials about parallelism. One simple example of parallelism in ABINIT will be shown.
* [[tutorial:paral_bandpw|Parallelism over bands and plane waves]] presents the combined k-point (K), plane-wave (G), band (B), spin/spinor parallelism of ABINIT (so, the "KGB" parallelism), for the computation of total energy, density, and ground state properties
* [[tutorial:paral_dfpt|Parallelism of response-function calculations]] - you need to be familiarized with the calculation of linear-response properties within ABINIT, see the tutorial [[tutorial:rf1|Response-Function 1 (RF1)]]

