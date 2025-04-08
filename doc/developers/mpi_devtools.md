## How to debug MPI jobs using gdb

Serial debuggers attached to individual processes in an MPI job are useful.
Run an MPI job (4 processes in this example) that opens one `xterm` terminal window for each process and then 
loads the ABINIT executable in the GNU debugger for that process using the command:

    mpirun -n 4 xterm -e gdb path_to_abinit_exe

Run ABINIT on each process separately by typing the following command on every terminal window opened during the
previous step:

    (gdb) run path_to_input_abi

Next, gdb can be used as in the serial case
[here](https://docs.abinit.org/developers/developers_howto/#how-to-debug-with-gdb). As an example, we wait for the
error in one of the processes, then print the *backtrace* on that process, select a stack frame (the first one here) 
and display additional information, such as values of function arguments before the error, with the commands:

    (gdb) bt
    (gdb) select-frame 1
    (gdb) info args

Remember to compile the code with the -g option in order to produce debugging information. 

## How to profile MPI jobs using NVTX/Nsight

Among tools that analyze MPI usage and performance, it is possible to annotate source code using [NVTX](https://nvidia.github.io/NVTX/) (NVIDIA Tools Extension Library) and trace the execution of an MPI job as a timeline of events per process using the profiler [NVIDIA Nsight Systems](https://developer.nvidia.com/nsight-systems). Code annotation is activated by default in ABINIT when compiled on GPUs, using NVIDIA (includes NVTX) or AMD (ROCTX) annotation libraries. As of ABINIT version 10.3.5, NVTX annotation is also supported for CPUs (see *with_gpu_markers* input variable). A guide for minimal installation of selected NVIDIA developer tools for CPUs is provided here.

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
achieved by adding the following lines to your ac9 configuration file when compiling ABINIT on CPU with autotools:

    with_gpu_markers="yes"
    abi_gpu_nvtx_v3="yes"
    GPU_LIBS="-L/usr/local/cuda-12.8/lib64 -lnvToolsExt"

Note that these lines are not needed when GPU is enabled because markers are configured automatically 
from the CUDA library root. The serial profiler can be attached to individual MPI processes for 
generating a timeline of selected API events using the following command, for example here tracing MPI and NVTX events of a parallel ABINIT calculation: 
 
    mpirun -n 4 nsys profile --trace=mpi,nvtx path_to_abinit_exe path_to_input_abi

Once a report has been generated on your machine for each process, it can be opened (in single or multi-report view) 
in Nsight GUI with: 

    nsys-ui report1.nsys-rep

An output example of tracing NVTX/MPI API events on a single MPI process using NVIDIA Nsight Systems:

![nsight_screenshot](nsight.png)
