---
authors: J-M Beuken
---

# How to install ABINIT v9 on CentOS

This step-by-step recipe describes how to build ABINIT on Fedora/RHEL/CentOS *nix distributions.
 
Explicit testing was done with __CentOS 8.2__

[Quick Guide for the Impatient](#quick-guide-for-the-impatient)<br />
[Quick Guide for the Impatient (MKL version)](#quick-guide-for-the-impatient-mkl-version)

## Prerequisites

1.  __A fortran compiler__<br />
   Possible options:       
     - gfortran, the GNU compiler. ([https://gcc.gnu.org/](https://gcc.gnu.org/))
     - ifort, the intel compiler. This is a commercial compiler, slightly more complicated to use but more optimized for intel architecture.
2.  __A Python interpreter__ (**v3.8 recommended**)
3.  __A MPI library__ (If you want to benefit from parallelism; **recommended**).<br />
   Possible options:
     - Open MPI  ([http://www.open-mpi.org](http://www.open-mpi.org))
     - MPICH  ([http://www.mpich.org](http://www.mpich.org))
4.  __A Linear Algebra library__<br />
   Possible options:
     - MKL (Intel® Math Kernel Library) ([Free Download](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library/choose-download/linux.html)) (**recommended**)
     - OpenBLAS (An optimized BLAS library) ([https://www.openblas.net/](https://www.openblas.net/)) (**recommended with GNU**)
     - Netlib (blas,lapack,scalapack) ([https://www.netlib.org/](https://www.netlib.org/))
     - ATLAS (Automatically Tuned Linear Algebra Software) ([http://math-atlas.sourceforge.net/](http://math-atlas.sourceforge.net/))
5.  __Mandatory libraries__<br />
     - HDF5 (High-performance data management and storage suite) ([https://www.hdfgroup.org/solutions/hdf5/](https://www.hdfgroup.org/solutions/hdf5/))
     - NetCDF (Network Common Data Form) ([https://www.unidata.ucar.edu/software/netcdf/](https://www.unidata.ucar.edu/software/netcdf/))
     - libXC (library of exchange-correlation functionals) ([https://tddft.org/programs/libxc/download/](https://tddft.org/programs/libxc/download/))
6.  __Optional supported librairies__<br />
     - FFTW (library for computing the discrete Fourier transform) ([http://www.fftw.org/](http://www.fftw.org/)) (**recommended with GNU**)
     - libxml2 (XML C parser) ([http://xmlsoft.org/](http://xmlsoft.org/downloads.html)) 
     - Wannier90  ([http://www.wannier.org/](http://www.wannier.org/)) (**recommended**)
     - LibPSML ([https://esl.cecam.org/PSML](https://esl.cecam.org/PSML)) + xmlf90 ([https://gitlab.com/siesta-project/libraries/xmlf90](https://gitlab.com/siesta-project/libraries/xmlf90))

## Installation of tools and libraries

All mandatory libraries are installed through the DNF package manager.<br />
Concerning some optional libraries, it will be necessary to compile them from the sources.

Following step installs a possible selection of required packages from dnf <br />
for __a simple parallel ABINIT__ compilation ( mandatory libraries, MPICH, fftw3 and OpenBLAS)

1. __compiler__<br />
    `sudo dnf install gcc-gfortran`

2. __MPI libraries__ - choice for MPICH<br />
    `sudo dnf install mpich mpich-devel`

3. __linear algebra library__ - choice for OpenBLAS<br />
    `sudo dnf install openblas`

4. __mandatory libraries__ - choice for hdf5 parallel version<br />
    `sudo dnf install hdf5-mpich hdf5-mpich-devel`<br />
    `sudo dnf install netcdf-mpich-devel netcdf-fortran-mpich-devel`<br />
    `sudo dnf install libxc libxc-devel`

5. __fftw3__<br />
    `sudo dnf install fftw fftw-devel`

6. __python interpreter__<br />
    `sudo dnf install python3`

__!!! Important:__

Before continuing, it is important to test if your development environment is properly configured.<br />
For example, to check if the MPICH package is functional, you can execute the following command:

    mpif90 --version

If the output is:

    bash: mpif90: command not found...

then, you need to find out where the MPI wrappers are installed...<br />
When you installed the MPICH package via dnf, all directories can be displayed by using this command:

    rpm -ql mpich-devel | grep mpif90

the output is

    /usr/lib64/mpich/bin/mpif90

 __The $PATH variable need to be updated__:

     export PATH=/usr/lib64/mpich/bin:$PATH

     mpif90 --version

     GNU Fortran (GCC) 8.3.1 20191121 (Red Hat 8.3.1-5)
     /usr/lib64/mpich/bin/mpif90


## Compiling, testing and installing ABINIT

__Download ABINIT__.<br />
For normal users, it is advised to get the newest version from our [website](https://www.abinit.org/packages) (replace 9.0.4 by the newest version available).

    wget https://www.abinit.org/sites/default/files/packages/abinit-9.0.4.tar.gz
    tar xzf abinit-9.0.4.tar.gz
    cd abinit-9.0.4

__Create a working directory__:

    mkdir build && cd build

__To configure__:

    ../configure -q --with-config-file='my_config_file.ac'

where `my_config_file.ac` is either a self made config file.<br />
More on the configure options is presented in [next section](#the-config-file).

__Compile with__:

    make -j 4

__!!! Tip:__

    '-j 4'  means that 4 cores are used to compile but you can use more

__Testing with__:

    cd tests
    ./runtests.py fast -j 4

__!!! Important:__

The end of the test output __must__ be something like:

    Suite   failed  passed  succeeded  skipped  disabled  run_etime  tot_etime
    fast         0       0         11        0         0      27.72      27.98

    Completed in 9.95 [s]. Average time for test=2.52 [s], stdev=2.94 [s]
    Summary: failed=0, succeeded=11, passed=0, skipped=0, disabled=0

otherwise there is a __problem__ with the compilation: see [Troubleshooting](#troubleshooting)

__Install__ (optional):

    make install

## The config file

The configure command takes as input some variables and flags.<br />
For example:

    ../configure --with-mpi="yes"

tells ABINIT to enable MPI support.<br />
Any variable or flag can be found by typing:

    ../configure --help

Most options are detected automatically by ABINIT.<br />
For example, with the option `--with-mpi="yes"`, ABINIT will try to use the parallel fortran compiler (mpifort)
and detect directories with useful library and header files for MPI support.<br />

When a lot of options are used, it is advised to use a config file.<br />

__!!! Important !!!__

    The name of the options in the `.ac` files is in normalized form that is
    the initial `--` is removed from the option name and all the other `-` characters
    in the string are replaced by an underscore `_`.
    Following these simple rules, the configure option `--with-mpi` becomes `with_mpi`
    in the ac file.

 The `.ac` file for __our simple parallel ABINIT__ using OpenBLAS can be written as:

    # installation location
    prefix=$HOME/local

    # MPI settings
    with_mpi="yes"
    enable_mpi_io="yes"

    # linear algebra settings
    with_linalg_flavor="openblas"
    LINALG_LIBS="-L/usr/lib64 -lopenblas"

    # mandatory libraries
    with_hdf5="yes"
    with_netcdf="yes"
    with_netcdf_fortran="yes"
    with_libxc="yes"

    # FFT flavor
    with_fft_flavor="fftw3"
    FFTW3_LIBS="-L/usr/lib64 -lfftw3 -lfftw3f"

    # Enable Netcdf mode in Abinit (use netcdf as default I/O library)
    enable_netcdf_default="yes"


## To go further

- compiling optional libraries with fallback project: Wannier90, libPSML/XMLF90,...
- enabling OpenMP
- using libxml2


__-------------------------------------------------------------------------------------------------__


## Quick Guide for the Impatient

We will build ABINIT with the following components:

- GNU compilers
- MPICH
- OpenBLAS
- FFTW3

### Installing needed packages

    sudo dnf install gcc-gfortran
    sudo dnf install mpich mpich-devel
    sudo dnf install openblas
    sudo dnf install hdf5-mpich hdf5-mpich-devel
    sudo dnf install netcdf-mpich-devel netcdf-fortran-mpich-devel
    sudo dnf install libxc libxc-devel
    sudo dnf install fftw fftw-devel
    sudo dnf install python3

### Getting sources of ABINIT

    wget https://www.abinit.org/sites/default/files/packages/abinit-9.0.4.tar.gz
    tar xzf abinit-9.0.4.tar.gz
    cd abinit-9.0.4
    mkdir build && cd build
    export PATH=/usr/lib64/mpich/bin:$PATH

### Creating a config file

Edit a `config.ac` file:

    # installation location
    prefix=$HOME/local

    # MPI settings
    with_mpi="yes"
    enable_mpi_io="yes"

    # linear algebra settings
    with_linalg_flavor="openblas"
    LINALG_LIBS="-L/usr/lib64 -lopenblas"

    # mandatory libraries
    with_hdf5="yes"
    with_netcdf="yes"
    with_netcdf_fortran="yes"
    with_libxc="yes"

    # FFT flavor
    with_fft_flavor="fftw3"
    FFTW3_LIBS="-L/usr/lib64 -lfftw3 -lfftw3f"

    # Enable Netcdf mode in Abinit (use netcdf as default I/O library)
    enable_netcdf_default="yes"

### Compiling ABINIT

     ../configure -q --with-config-file='config.ac'
     make -j 8

### Testing ABINIT

    cd tests
    export OPENBLAS_NUM_THREADS=1
    ./runtest.py fast -j 8 --no-logo

### Installing ABINIT

    make install


__-------------------------------------------------------------------------------------------------__


## Quick Guide for the Impatient (MKL version)

We will build ABINIT with the following components:

- GNU compilers
- MPICH
- MKL

### Installing needed packages

    sudo dnf install gcc-gfortran
    sudo dnf install mpich mpich-devel
    sudo dnf install hdf5-mpich hdf5-mpich-devel
    sudo dnf install netcdf-mpich-devel netcdf-fortran-mpich-devel
    sudo dnf install libxc libxc-devel
    sudo dnf install python3

### Getting sources of ABINIT

    wget https://www.abinit.org/sites/default/files/packages/abinit-9.0.4.tar.gz
    tar xzf abinit-9.0.4.tar.gz
    cd abinit-9.0.4
    mkdir build && cd build
    export PATH=/usr/lib64/mpich/bin:$PATH

### Creating a config file

 Edit a `config.ac` file:

    # installation location
    prefix=$HOME/local

    # MPI settings
    with_mpi="yes"
    enable_mpi_io="yes"

    # linear algebra settings
    with_linalg_flavor="mkl"
    LINALG_CPPFLAGS="-I${MKLROOT}/include"
    LINALG_FCFLAGS="-I${MKLROOT}/include"
    LINALG_LIBS="-L${MKLROOT}/lib/intel64 -Wl,--start-group  -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group"

    # mandatory libraries
    with_hdf5="yes"
    with_netcdf="yes"
    with_netcdf_fortran="yes"
    with_libxc="yes"

    # FFT flavor
    with_fft_flavor="dfti"
    FFT_FCFLAGS="-I${MKLROOT}/include"

    # Enable Netcdf mode in Abinit (use netcdf as default I/O library)
    enable_netcdf_default="yes"

### Compiling ABINIT

     ../configure -q --with-config-file='config.ac'
     make -j 8

### Testing ABINIT

    cd tests
    export MKL_NUM_THREADS=1
    ./runtest.py fast -j 8 --no-logo

### Installing ABINIT

    make install


__-------------------------------------------------------------------------------------------------__


## Troubleshooting

-   

