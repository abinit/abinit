---
authors: J-M Beuken
---

# How to install ABINIT v9 on CentOS

This step-by-step recipe describes how to build ABINIT on Fedora/RHEL/CentOS *nix distributions.
Tested with __CentOS 8.2__

[Quick Guide for the Impatient](#quick-guide-for-the-impatient)

[Quick Guide for the Impatient (MKL version)](#quick-guide-for-the-impatient-mkl-version)

## Prerequisites

1.  __Fortran compiler__

    Possible options:

    - gfortran, the GNU compiler. ([https://gcc.gnu.org/](https://gcc.gnu.org))
    - ifort, the intel compiler. This is a commercial compiler, slightly more complicated
      to use but more optimized for intel architecture.

2.  __Python interpreter__ (**v3.7+ recommended**)

3.  __MPI library__ (if you want to benefit from parallelism; **recommended**).

     Possible options:

     - [Open MPI](http://www.open-mpi.org)
     - [MPICH](http://www.mpich.org)

4.  __Linear Algebra library__

    Possible options:

    - MKL (Intel® Math Kernel Library): [Free Download](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library/choose-download/linux.html), **recommended** for performance
    - [OpenBLAS](https://www.openblas.net): An optimized BLAS library, **recommended with GNU**.
    - [Netlib](https://www.netlib.org): blas, lapack, scalapack
    - [ATLAS](http://math-atlas.sourceforge.net/): Automatically Tuned Linear Algebra Software

5.  __Mandatory libraries:__

    - [HDF5](https://www.hdfgroup.org/solutions/hdf5/): High-performance data management and storage suite
    - [NetCDF](https://www.unidata.ucar.edu/software/netcdf): Network Common Data Form
    - [libXC](https://tddft.org/programs/libxc/download/): Library of exchange-correlation functionals

6.  __Optional supported librairies:__

    - [FFTW3](http://www.fftw.org/): Library for computing the discrete Fourier transform, **recommended with GNU**
    - [libxml2](http://xmlsoft.org/downloads.html): XML C parser, recommended for multibinit
    - [Wannier90](http://www.wannier.org)
    - [LibPSML](https://esl.cecam.org/PSML) + [xmlf90](https://gitlab.com/siesta-project/libraries/xmlf90)
      to read pseudopotentials in psml format

## Installation of tools and libraries

All mandatory libraries are installed through the DNF package manager.
For other optional libraries, compilation from source is needed.

The steps required to install MPICH, fftw3 and OpenBLAS with dnf and compile
__a relatively simple parallel version of ABINIT__ are summarized below:

1. __Install the compiler__

    `sudo dnf install gcc-gfortran`

2. __Install MPI library__ (MPICH)

    `sudo dnf install mpich mpich-devel`

3. __Install linear algebra library__  (OpenBLAS)

    `sudo dnf install openblas`

4. __Install other mandatory libraries__  (use hdf5 with support for parallel MPI-IO)

    `sudo dnf install hdf5-mpich hdf5-mpich-devel`

    `sudo dnf install netcdf-mpich-devel netcdf-fortran-mpich-devel`

    `sudo dnf install libxc libxc-devel`

5. __Install fftw3__

    `sudo dnf install fftw fftw-devel`

6. __Install python interpreter__

    `sudo dnf install python3`

!!! important

    Before continuing, it is important to test whether your development environment is properly configured.
    To check whether the MPICH package is installed, execute the following command:

    ```sh
    mpif90 --version
    ```

    If the output is:

    ```sh
    bash: mpif90: command not found...
    ```

    then, you need to find out where the MPI wrappers are installed.

    If you installed the MPICH package via dnf, the installation directories can be obtained by using e.g.

    ```sh
    rpm -ql mpich-devel | grep mpif90
    ```

    that should print

    ```sh
    /usr/lib64/mpich/bin/mpif90
    ```

    The $PATH variable needs to be updated:

    ```sh
    export PATH=/usr/lib64/mpich/bin:$PATH

    mpif90 --version

    GNU Fortran (GCC) 8.3.1 20191121 (Red Hat 8.3.1-5)
    /usr/lib64/mpich/bin/mpif90
    ```

## Compiling, testing and installing ABINIT

__Download ABINIT__.

For normal users, it is advised to get the latest stable version
from our [website](https://www.abinit.org/packages) (replace 9.0.4 by the newest version available).

```sh
wget https://www.abinit.org/sites/default/files/packages/abinit-9.0.4.tar.gz
tar xzf abinit-9.0.4.tar.gz
cd abinit-9.0.4
```

__Create a working directory__:

```sh
mkdir build && cd build
```

__Configure with__:

```sh
../configure --with-config-file='my_config_file.ac'
```

where `my_config_file.ac` is an external file providing all the configuration flags and options.
More on the configure options is presented in [next section](#the-config-file).

__Compile with__:

```
make -j 4
```

where `-j 4` means that 4 cores are used to compile. Adjust this value according to number of
physical cores available on your machine.

To run the test suite, issue:

```sh
cd tests
./runtests.py fast -j 4
```

!!! Important

    At the end of the test, one should get something like:

    ```
    Suite   failed  passed  succeeded  skipped  disabled  run_etime  tot_etime
    fast         0       0         11        0         0      27.72      27.98

    Completed in 9.95 [s]. Average time for test=2.52 [s], stdev=2.94 [s]
    Summary: failed=0, succeeded=11, passed=0, skipped=0, disabled=0
    ```

otherwise there is a __problem__ with the compilation: see [Troubleshooting](#troubleshooting)

__Install__ (optional):

    make install

## The configuration file

The configure command takes in input variables and flags.
For example:

```sh
../configure --with-mpi="yes"
```

tells ABINIT to enable MPI support.
All the variables and flags supported by the script can be found by typing:

```sh
../configure --help
```

Some options are detected automatically by the script.
For example, with the option `--with-mpi="yes"`, ABINIT will try to use the parallel fortran compiler
found in $PATH (e.g. mpifort) and will try to detect the directories containing the associated libraries
and the header files required by MPI.

When a lot of options are needed, it is advised to use a config file.

The `.ac` file for __our simple parallel ABINIT__ build based on OpenBLAS is:

```sh
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
```

!!! Important

    The name of the options in the `.ac` files is in normalized form that is
    the initial `--` is removed from the option name and all the other `-` characters
    in the string are replaced by an underscore `_`.
    Following these simple rules, the configure option `--with-mpi` becomes `with_mpi`
    in the ac file.

<!--
## To go further

- compiling optional libraries with the fallback project: Wannier90, libPSML/XMLF90.
- enabling OpenMP
- using libxml2
-->

## Quick Guide for the impatient

We will build ABINIT with the following components:

- GNU compilers
- MPICH
- OpenBLAS
- FFTW3

### Installing required packages

```sh
sudo dnf install gcc-gfortran
sudo dnf install mpich mpich-devel
sudo dnf install openblas
sudo dnf install hdf5-mpich hdf5-mpich-devel
sudo dnf install netcdf-mpich-devel netcdf-fortran-mpich-devel
sudo dnf install libxc libxc-devel
sudo dnf install fftw fftw-devel
sudo dnf install python3
```

### Getting the ABINIT tarball

```sh
wget https://www.abinit.org/sites/default/files/packages/abinit-9.0.4.tar.gz
tar xzf abinit-9.0.4.tar.gz
cd abinit-9.0.4
mkdir build && cd build
export PATH=/usr/lib64/mpich/bin:$PATH
```

### Creating a config file

Edit a `config.ac` file:

```sh
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
```

### Compiling ABINIT

```sh
../configure -q --with-config-file='config.ac'
make -j 8
```

### Testing ABINIT

```sh
cd tests
export OPENBLAS_NUM_THREADS=1
./runtest.py fast -j 8 --no-logo
```

### Installing ABINIT

```sh
make install
```

## Quick Guide for the impatient (MKL version)

We will build ABINIT with the following components:

- GNU compilers
- MPICH
- MKL

### Installing needed packages

```sh
sudo dnf install gcc-gfortran
sudo dnf install mpich mpich-devel
sudo dnf install hdf5-mpich hdf5-mpich-devel
sudo dnf install netcdf-mpich-devel netcdf-fortran-mpich-devel
sudo dnf install libxc libxc-devel
sudo dnf install python3
```

### Getting sources of ABINIT

```sh
wget https://www.abinit.org/sites/default/files/packages/abinit-9.0.4.tar.gz
tar xzf abinit-9.0.4.tar.gz
cd abinit-9.0.4
mkdir build && cd build
export PATH=/usr/lib64/mpich/bin:$PATH
```

### Creating a configuration file

Edit a `config.ac` file:

```sh
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
```

### Compiling ABINIT

```sh
../configure --with-config-file='config.ac'
make -j 8
```

### Testing ABINIT

```sh
cd tests
export MKL_NUM_THREADS=1
./runtest.py fast -j 8 --no-logo
```

### Installing ABINIT

```sh
make install
```
