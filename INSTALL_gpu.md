# GPU support in Abinit

> IMPORTANT:
GPU support is currently highly *EXPERIMENTAL* and should be used by experienced developers only.

Since April 2024, a new GPU implementation of ABINIT has become available.
A new user manual for ABINIT on GPU is currently being written and will be made available soon. In the meantime, you may refer to the description of the [[gpu_option]] variable for guidance, as entry point.

*The previous implementation from 2013 (based on `cuda`/`magma`) is now obsolete. However, it remains available for testing purposes only.*

To use ABINIT with GPUs, two distinct programming models are available:
- [OpenMP offload](https://www.openmp.org/specifications/)
  Available for mostly all ground-state properties (DFPT is ongoing). `LOBPCG`and `Chebyshev filtering` iterative eigensolvers available (see [[wfoptalg]]).
-  [cuda](https://docs.nvidia.com/cuda)+[Kokkos](https://github.com/kokkos/kokkos)+[YAKL](https://github.com/mrnorman/YAKL)
  Available only for the calculation of the ground-state total energy, using the `Chebyshev Filtering` (ChebFi) iterative eigensolver.

To utilize either of these programming models, it is essential to compile ABINIT with specific configurations by enabling the relevant options during the setup stage and selecting the correct value for the input parameter [[gpu_option]].

If you wish to conduct tests on a CPU+GPU computing architecture, you can look to the three (at present) example configuration files for inspiration,
corresponding to three different tested installations. 
They are available in the directory doc/build of the package, and also available from
<https://github.com/abinit/abinit/tree/master/doc/build> .
Simply adapt them to suit your specific CPU+GPU hardware/software:

- **Architecture 1**
   CPU: AMD Milan EPYC 7763
   GPU: Nvidia A100 SXM4 80 Go
   Programming model: `openMP offload`
   Compilers: Nvidia HPC compilers
   Libraries: MKL, Cuda
    
- **Architecture 2**
   CPU: AMD Milan EPYC 7543
   GPU: Nvidia A100 SXM4 80 Go
   Programming model: `openMP offload`
   Compilers: Nvidia HPC compilers
   Libraries: MKL, Cuda
      
- **Architecture 3**
   CPU: AMD Trendo EPYC 7A53
   GPU: AMD Instinct MI250X 64 GB
   Programming model: `openMP offload`
   Compilers: CRAY
   Libraries: FTTW3, libSCI, HIP+ROCm
