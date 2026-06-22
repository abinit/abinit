---
title: Parallelization of ground-state using GPU
authors: MS
---

# Parallelization of ground-state calculations using GPU

## Explore the *k-points/plane waves/bands* parallelization wih GPU enabled

This tutorial discusses how to accelerate ground-state calculations with GPUs using ABINIT.

It is assumed that you followed previous [[tutorial:paral_bandpw|tutorial exploring k-points/plane waves/bands parallelism]], as
the same concept is leveraged here.

This tutorial should take less than an hour and assumes that you are running on a system with 128 CPU cores and 4 NVIDIA GPUs.

You may run on fewer resources, but be aware that consumer-grade GPUs (such as NVIDIA GeForce RTX and AMD Radeon models) typically provide limited double-precision performance.
You need to run on a server-grade GPU such as those usually found in computing clusters otherwise GPU execution may be significantly slower.

You are supposed to know already some basics of ABINIT.
Some useful references: [[cite:Levitt2015]], [[cite:Bottin2008]], [[cite:Knyazev2001]]

If you are interested in details of the GPU implementation and its performance, please refer to [[cite:Lygatsika2026]].

[TUTORIAL_README]

## Introduction

*Before continuing you might work in a different subdirectory as for the other
tutorials. Why not work_paral_gpu?*

!!! important 

    In what follows, the names of files are mentioned as if you were in this subdirectory.
    All the input files can be found in the `\$ABI_TESTS/tutoparal/Input` directory.
    You can compare your results with reference output files located in `\$ABI_TESTS/tutoparal/Refs`.

    In the following, when "run ABINIT over _nn_ CPU cores" appears, you have to use
    a specific command line according to the operating system and architecture of
    the computer you are using. This can be for instance: `mpirun -n nn abinit input.abi`
    or the use of a specific submission file.


Conceptually, a GPU (Graphics Processing Unit) can be thought of as a chip containing thousands of relatively simple cores, while a CPU contains a smaller number of more sophisticated and powerful cores.

In other words, GPUs provide massive data parallelism through thousands of lightweight execution units, while CPUs are optimized for low-latency execution and complex control flow.

GPUs act as accelerators attached to a CPU node, allowing a running application to offload
execution of selected parts of the calculation that can take advantage of the massive parallelism capabilities of GPUs.

A GPU typically has its own dedicated memory, meaning data transfers happen between CPU and GPU for computation on GPU to occur.
However, this is not true on some recent systems such as NVIDIA Grace Superchips (ex: GH200, GB200) that provide a unified memory paradigm.
In such systems, CPU and GPU can access a shared memory space, reducing or eliminating some explicit data transfers.

In ABINIT, GPU porting efforts have been made to accelerate performance-critical parts of ground-state calculations.
Enabling GPU support is straightforward, but achieving good performance also requires choosing appropriate parallelization input parameters.

In this tutorial, we will discuss:

* How to enable GPU for a ground-state calculation
* Which parallelism settings to enable for better performance with GPU

## GPU acceleration on a huge case

The first system we're going to accelerate is a moderately big system with thousands of electronic bands.

In the GPU implmentation of ABINIT, bands are typically one of the dimensions along which GPU parallelism is exploited.

The following parallelism variables are still relevant in GPU-enabled calculations:

 - [[np_spkpt]] (number of processes for spin and k points),
 - [[npband]] (number of processes for bands),

And of course, [[paral_kgb]], to enable those variables.
Compared to a CPU execution, only [[npfft]] is unused for GPU calculation.

Another important parameter is the wavefield optimisation algorithm:

 - [[wfoptalg]]=111 : Chebyshev Filtering, aka **CHEBFI**, *favoured on big systems*
 - [[wfoptalg]]=114 : Locally Optimised Block Processing Conjugate Gradient, aka **LOBPCG**, *historic*
 - [[rmm_diis]]=1 : Residual Minimization Method, Direct Inversion in the Iterative Subspace, aka **RMM-DIIS**, *works with GPU but not optimised*

!!! important

    When [[wfoptalg]] isn't set, default algorithm Conjugate Gradient ([[wfoptalg]]=0) is used.
    **This algorithm only runs on CPU**, so be sure to set [[wfoptalg]] to use either CHEBFI or LOBPCG

Out of the 3, CHEBFI is the more recent and was designed to run efficiently on larger systems by reducing calls to costly Rayleigh-Ritz method, compared to LOBPCG.
Therefore, CHEBFI is generally the preferred algorithm on GPUs, while LOBPCG is more of a fallback.

Only ChebFi will be considered in this tutorial.

### First run using CPU

First, we will not distribute the workload over k-points and only rely on band parallelism.

The first case here is a 128-atom titanium system, moderately big, with few SCF steps.

As a first step, try running this case on CPU, by using [[npband]] to distribute MPI tasks on bands and [[wfoptalg]] to use ChebFi algorithm.
Let's use 32 MPI tasks and 4 OpenMP threads.

It should take about a minute to complete.

###################################################################################################################
###################################################################################################################
Copy the `tparal_gpu_01.abi` file the tutorial directory into your working directory.

{% dialog tests/tutoparal/Input/tparal_gpu_01.abi %}
###################################################################################################################
###################################################################################################################

It is rather slow, so let's see how to accelerate execution using GPU.


### Accelerate using GPU

The first crucial step is to add the following input parameter :

```
 gpu_option "GPU_OPENMP"
```
or (less verbose) :
```
 gpu_option 2
```

As for MPI parallelism, we're going to use as many tasks as we have GPUs.
With 4 GPUs, we use 4 MPI tasks... so set [[npband]] accordingly.

###################################################################################################################
###################################################################################################################
Copy the `tparal_gpu_02.abi` file the tutorial directory into your working directory.

{% dialog tests/tutoparal/Input/tparal_gpu_02.abi %}
###################################################################################################################
###################################################################################################################

It should run significantly faster.


### Information about GPU in ABINIT output

If you inspect ABINIT output, you will see 2 informative sections related to GPU.

The first one is located after the input parameter summary and displays a summary of the GPU selected by MPI task 0:
```
 ________________________ Graphic Card Properties _______________________________
~
    Device                0 : NVIDIA A100-SXM4-40GB
 Revision number:                   8.0
 Total amount of global memory:     40440 Mbytes
 Clock rate:                        1.4 GHz
 Number of processors/cores:          108/       32
 Max FP64 GFLOPS:                  4872 GFP
 Total  constant memory:            65536 bytes
 Shared memory per block:           49152 bytes
 Number of registers per block:     65536
 UUID:                              000000CB-00000000-00000000-00000000
~
 ________________________________________________________________________________
 Using CUDA version: 12.4
```

The second one is located at the first SCF step and is a summary of GPU memory consumption estimate:
```
GPU memory consumption estimate per MPI task for K-point   1:
 Considered available memory             :  36396.336 MiB

|                 Buffers governed by blocking/slicing                 |
|:--------------------------|---------:|-------------:|---------------:|
|  gemm_nonlop_projectors   |    1 blk |  npw,*natom* |     299.32 MiB |
|  fourwf (fofr work array) |    1 blk | npw,*bandpp* |     615.09 MiB |
|  xFFT~internal buffers    |       NA |           NA |       0.00 MiB |

|   Static buffers, computed once and permanently on card    |
|:-------------------------|--------------:|----------------:|
|  invovl (mkinvovl)       |        natom  |       81.01 MiB |
|  chebfi2                 |          npw  |      299.39 MiB |
|  hamiltonian arrays      |          npw  |       19.16 MiB |

|  Work buffers (mostly sized after bandpp or nblock_lobpcg) |
|:-------------------------|--------------:|----------------:|
|  mkrho~vtowfk_extra      |   npw,bandpp  |       33.26 MiB |
|  hegvd                   |       bandpp  |       48.00 MiB |
|  prep_nonlop             |   npw,bandpp  |       99.77 MiB |
|  invovl                  | natom,bandpp  |       45.00 MiB |
|  chebfi2 (RR buffers)    |        nband  |       32.01 MiB |
|  chebfiwf (cg,resid,eig) |    npw,nband  |       33.29 MiB |

Sum                                      :    1605.31 MiB
```

This summary is in Markdown format and can be display as tables on some platforms.

Its purpose is to outline which arrays are consuming the most GPU memory.

Normally, ABINIT tries to split some part of its computation in blocks to fit in GPU memory if necessary but may fail sometimes.

In case of failure:

- if **buffers governed by blocking/slicing** are too big, try setting [[gpu_nl_splitsize]] or [[gpu_nfft_blocks]] yourself
- if **static buffers** are too big, your system is too big (too many atoms and/or planewaves)
- if **work buffers** are too big, try running on more nodes



## GPU acceleration on a case with many k-points

The second case is a 31-atoms of gold system, and has a small number of bands but many k-points

As explained earlier, GPU massive parallelism capabilities are mostly exploited on the number of bands.
This means that with a small number of bands, the GPU may not be fully utilized.

In that case, it may be advantageous to run using more than one MPI task per GPU, especially if you have many k-points.

### First step with CPU

As a first step, run this case on CPU, still using 32 MPI tasks and 4 OpenMP threads.

Instead of using [[npband]], use [[np_spkpt]] to distribute MPI tasks among k-points, as it is more efficient than [[npband]].

###################################################################################################################
###################################################################################################################
Copy the `tparal_gpu_03.abi` file the tutorial directory into your working directory.

{% dialog tests/tutoparal/Input/tparal_gpu_03.abi %}
###################################################################################################################
###################################################################################################################

### Accelerate using GPU

Then, for enabling GPU, proceed as usual : add the parameter [[gpu_option]] and set it to `GPU_OPENMP`.

Considering 4 GPUs on your node, you will also use as many MPI tasks, and consequently set [[np_spkpt]] accordingly.

###################################################################################################################
###################################################################################################################
Copy the `tparal_gpu_04.abi` file the tutorial directory into your working directory.

{% dialog tests/tutoparal/Input/tparal_gpu_04.abi %}
###################################################################################################################
###################################################################################################################

### Using many MPI tasks per GPU

While previous run was accelerated by enabling GPU, we can do better by using many MPI tasks per GPU.

The reason is that this system is relatively small, with only a limited number of bands and plane waves. As a result, a single MPI rank may not generate enough work to fully utilize the GPU.

In the specific case of NVIDIA GPUs, you'll need to enable a background service named Multiple-Process Service, a.k.a **MPS**.

It consists of initializing a single CUDA context and splitting GPU resources so that each MPI task using the GPU won't compete for the same hardware resources.

!!! important

    Enablement of NVIDIA MPS may depend on your cluster rules or managed through a job scheduler like SLURM, so check your local documentation.


Considering a local execution, you can enable MPS and split GPU resources using :
```
# Check all GPUs compute mode are set to DEFAULT (EXCLUSIVE means it won't work)
nvidia-smi -q -d COMPUTE

# Launch MPS daemon
nvidia-cuda-mps-control -d

# Split GPUs in 4 shares (1/4 = 25 %)
echo "set_default_active_thread_percentage 25.0" | nvidia-cuda-mps-control
```

By splitting GPUs in 4 shares, we can now allocate 4 MPI tasks per GPU, or 16 in total, so we need to change [[np_spkpt]] accordingly.

###################################################################################################################
###################################################################################################################
Copy the `tparal_gpu_05.abi` file the tutorial directory into your working directory.

{% dialog tests/tutoparal/Input/tparal_gpu_05.abi %}
###################################################################################################################
###################################################################################################################

Performance should have been improved again.
Note that 4 MPI ranks per GPU is only an example; the optimal value depends on the system size, GPU model, and workload characteristics.


## 

