## How to profile code

Profiling is useful for identifying potential code optimizations in general and becomes an important diagnostic tool 
when porting ABINIT to GPUs in particular. An internal timer is implemented in ABINIT and its report is enabled using 
the keyword [*timopt*](https://docs.abinit.org/variables/gstate/#timopt). Moreover, memory profiling on CPU is 
activated internally using the *enable_memory_profiling* compilation flag. In order to profile memory on GPU or to 
perform a more advanced analysis, we need to use a dedicated tool. Here we focus on the ones provided by GPU manifacturers NVIDIA and AMD.

### How to profile MPI jobs using NVTX/Nsight

The first tool is not limited to GPU and presents an interest for timeline traces of parallel processes on CPU as well.
We make a general presentation for CPU while an analogous procedure applies to GPU.

Among tools that analyze MPI usage and performance, it is possible to annotate source code using 
[NVTX](https://nvidia.github.io/NVTX/) (NVIDIA Tools Extension Library) and trace an MPI job execution as a timeline of 
API events per process using the profiler [NVIDIA Nsight Systems](https://developer.nvidia.com/nsight-systems). Code 
annotation libraries are activated by default in ABINIT when compiled on GPU, since NVIDIA CUDA Toolkit includes NVTX 
(and AMD includes ROCTX). As of ABINIT version 10.3.6, NVTX annotation is also supported on CPU (see *with_gpu_markers* 
input variable). A guide for minimal installation of selected NVIDIA developer tools required to use NVTX on CPUs is 
provided here. Note that we don't need to install the entire NVIDIA CUDA Toolkit to profile code on CPU.

The NVIDIA Tools Extensions (NVTX) API can be installed on Linux with:

    sudo dnf install cuda-nvtx-12

You can check the location of the installed library using the command `locate nvToolsExt`.
Environment variables must be set by adding the following lines to your .bashrc (assuming here version 12.8 has been 
installed): 

    export PATH=/usr/local/cuda-12.8/bin:$PATH
    export LD_LIBRARY_PATH=/usr/local/cuda-12.8/lib64:$LD_LIBRARY_PATH

Nsight Systems installer can be directly downloaded from NVIDIA
[website](https://developer.nvidia.com/nsight-systems/get-started). We recommend to download the Full Version for 
profiling from the GUI. A minimal installation for profiling from the CLI only is also available but not used here. 
To use NVTX API you need to link ABINIT with `nvToolsExt` .so library and activate the NVTX markers macros. This can be
achieved by adding the following lines to your `autoconf` configuration file when compiling ABINIT on CPU with 
autotools:

    with_gpu_markers="yes"
    abi_gpu_nvtx_v3="yes"
    GPU_LIBS="-L/usr/local/cuda-12.8/lib64 -lnvToolsExt"

Note that these lines are not needed when GPU is enabled because marker libraries are configured automatically 
from the CUDA library root. The serial profiler can be attached to individual MPI processes for 
generating a timeline of selected API events using the following command, for example here tracing MPI and NVTX events 
of a parallel ABINIT calculation: 

    mpirun -n 4 nsys profile --trace=mpi,nvtx path_to_abinit_exe path_to_input_abi

Once a report has been generated on your machine for each process, it can be opened (in single or multi-report view) 
in Nsight GUI with: 

    nsys-ui report1.nsys-rep

An output example of tracing NVTX/MPI API events on a single MPI process using NVIDIA Nsight Systems:

![nsight_screenshot](nsight.png)

### How to profile GPU kernels using the roofline model

The roofline model can be used to visualize the achieved performance and arithmetic intensity of a GPU kernel. This 
information allows to assess whether a GPU kernel execution makes effective use of the available compute capabilities 
of a GPU architecture, or whether it underutilizes the resources and may benefit from optimization. For instance, for
CUDA kernels, the roofline can reveal if Tensor cores are used or not. 

The roofline is generated in two steps. First, a set of required metrics is produced using 
[Nsight Compute](https://docs.nvidia.com/nsight-compute/) for CUDA kernels and 
[rocPROFv3](https://rocm.docs.amd.com/projects/rocprofiler-sdk/en/latest/how-to/using-rocprofv3.html) for 
ROCm kernels. Then, the metrics are post-processed using user-defined scripts in order to produce the roofline chart 
depending on the GPU specifications namely the peak performance and memory bandwidth. We use the Python scripts 
provided by the Adastra supercomputing center of [GENCI](https://www.genci.fr/en), equipped with AMD GPUs, to 
demonstrate how to produce rooflines for GPU kernels in ABINIT. For more information, we refer to the extensive Adastra 
user documentation on [GPU roofline](https://dci.dci-gitlab.cines.fr/webextranet/software_stack/tools/index.html#id1).

A production build can be used to generate roofline metrics, there is no need for debug build. 
An example of configuration file can be found in *abinit/doc/build
/GPU_InstinctMI250X+EPYC7453.ac9*. ABINIT implements ROCm annotations and markers that are included in the profiler reports for the ease of analysis per
separated code parts. It also implements useful macros for profiler start/stop calls, for both NVTX and ROCm. These allow 
to limit the profiler execution to specific parts of the code and to reduce overhead. We use them for readability of the
roofline in order to avoid overlapping points. Do not forget the compilation option *with_gpu_markers="yes"* in order to be able to use 
annotations and profiler switches for ROCm.

In this example, we sandwich a critical part of the code that is computationally demanding, namely the Rayleigh-Ritz 
procedure. This is possible by using the 

first adding a line at the beginning of *src/98_main/abinit.F90* to switch off the 
profiler as early as possible 

    #if defined(HAVE_GPU_MARKERS)
        NVTX_INIT()
        NVTX_PROFILER_STOP()
    #endif
    ! .. rest of the code ..

Note that NVTX init simply fills the internal nvtx_names and nvtx_ids arrays used to label regions so putting the
stopped before or after essentially makes no difference.

then switch on the profiler in the interesting part in *src/48_diago/m_chebfi2.F90*: 

    NVTX_PROFILER_START()                     ! start profiler

    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_RR)     ! start annotation
    call xg_RayleighRitz(chebfi%X,chebfi%AX%self,chebfi%BX%self,eigen,ierr,0,tim_RR,& 
        chebfi%gpu_option,solve_ax_bx=.true.)
    ABI_NVTX_END_RANGE()                      ! end annotation

    NVTX_PROFILER_STOP()                      ! end profiler

For this test we run a ground-state calculation for 320 atoms of Ga2O3 with 1536 bands, a cutoff of 18 Hartree and 
1 k-point. We allocate an entire compute node (8 GPUs for Adastra MI25OX), deploy 8 parallel tasks with 8 logical cores
each and make sure to export multi-threaded *OMP_NUM_THREADS=8* so that the used solver is rocSOLVER and not ScaLAPACK. 
We use a wrapper to call rocPROFv3 for the rank 0 process only and avoid overhead: 

    $ cat rocprofv3_wrapper.sh
    #!/bin/bash

    if [ "${SLURM_PROCID}" == "0" ]; then
	    exec -- rocprofv3 --stats --kernel-trace --marker-trace \
		    --input=counters.txt --output-format=csv -o "${OUTPUT_FILE}" \
		    -- "${@}"
    else
	    exec -- "${@}"
    fi

where the *counters.txt* file contains the counters for binary64 profiling recommended in the Adastra
documentation (see [Performance
counters](https://dci.dci-gitlab.cines.fr/webextranet/software_stack/tools/index.html#performance-counters) section),
that we copy here for self-consistency:

    $ cat counters.txt
    pmc: TCC_EA_RDREQ_32B_sum TCC_EA_RDREQ_sum TCC_EA_WRREQ_sum TCC_EA_WRREQ_64B_sum SQ_INSTS_VALU_ADD_F64 SQ_INSTS_VALU_MUL_F64 SQ_INSTS_VALU_FMA_F64 SQ_INSTS_VALU_TRANS_F64 SQ_INSTS_VALU_MFMA_MOPS_F64

Note that without *--selected-regions*, rocprofv3 starts profiling immediately.

We then make the profiler call in our SLURM job as

    srun rocprofv3_wrapper.sh ${ABINIT} ga2o3.abi

The output is a counter collection stored in CSV. We post-process it in order to plot the roofline using the 
Python scripts provided in the Adastra documentation 
(see [Roofline](https://dci.dci-gitlab.cines.fr/webextranet/software_stack/tools/index.html#roofline) section).

The kernels of interest are then diagnosed using the roofline: GEMM operation by rocBLAS is compute-bound, while the 
HEGV eigensolver of rocSOLVER is memory-bound.

![roofline_screenshot](roofline.png)

Note that for the moment we cannot use *-selected-regions* option as ROCm in the Adastra toolchain is version 6.4.3 while
the option is available > 7.2.0. For now we used the fallback described in the official ROCm documentation
[here](https://rocm.docs.amd.com/projects/rocprofiler-sdk/en/develop/how-to/using-rocprofiler-sdk-roctx.html#profiler-control-with-selected-regions). AMD explicitly states that counter collection for selected ROCTx regions was added in ROCm 7.14. ROCm 7.14 release notes. Citing the documentation, "Counter collection for selected regions is available in ROCm 7.14.0." see more [here](https://rocm.docs.amd.com/en/docs-7.14.0/about/release-notes.html#selective-roctx-region-profiling-with-counter-collection).
