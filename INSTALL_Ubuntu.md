---
authors: J. Van Bever
---

# How to install ABINIT on Ubuntu

This file describes how to install ABINIT on Ubuntu.
The procedure for other Linux-based systems is similar, but the package manager and
the location of the libraries may differ.
Explicit testing was done with Ubuntu 19.10.

This page discusses how to install ABINIT with two different approaches:

- Using the apt package manager to install a precompiled version of abinit
- Compiling abinit from source using external libraries installed via apt.

Other possibilities (such as installation using the binary files provided by *conda*) are not discussed.

## Using apt

To install a precompiled version of ABINIT, just type:

    sudo apt install abinit

then enter your user password and ABINIT should install smoothly with all its dependencies.

!!! warning

    This version of ABINIT is most likely an old one.
    Use `abinit -v` to print the version number.

## Compiling from source under Ubuntu

### Introduction

The source code of the latest stable release is available at [https://www.abinit.org/packages](https://www.abinit.org/packages).
Compiling from source enables one to optimize the build, activate MPI parallelism, customize linear algebra libraries etc.
All the required packages are available via *apt* and can be installed with the syntax:

    sudo apt install [package]

Use the commands `apt search` and `apt show` to search for a package or to get info about a particular package.

Ubuntu packages install their executables in the bin folder (e.g. /usr/bin).
The location of the executable can be found via the unix `which` command.
For example

    which gfortran

will show the location of the gfortran executable provided the binary is in $PATH.

Other packages of interest when compiling code from source
are the so-called develop packages (denoted by the `-dev` suffix) that contain among others:

  - header files (extension `.h`) providing the declaration of prototypes and named constants.
    These files are usually installed in an include folder (e.g. /usr/include).
  - shared library files (extension `.so`).
    that are usually installed in a lib folder (e.g. /usr/lib).

To obtain the location of such files, use the command

    dpkg -L [package]

from the *dpkg* package

### Prerequisites

The prerequisites are first discussed qualitatively, because the installation may
depend on the linux distribution. Then we discuss how to compile the source code.

A possible list of prerequisites (tested for Ubuntu 19.10) is:

1. Ubuntu or similar Ubuntu-based distributions

2. A Fortran compiler. Possible options:

      - gfortran, the GNU compiler. It is open source and available via e.g. apt (gfortran package).
      - ifort, the intel compiler. This is a commercial compiler, slightly more complicated
        to use but more optimized for intel architecture.

3. A MPI library installed (If you want to benefit from parallelism; recommended).
   Possible options:

      - Open MPI from apt (`libopenmpi-dev` package) or [http://www.open-mpi.org](http://www.open-mpi.org)
      - MPICH from apt (`libmpich-dev` package) or [http://www.mpich.org](http://www.mpich.org)

    Depending on your distribution, you might need to manually add the `mpi-default-dev` package,
    a metapackage for both MPI libraries.

4. A Linear Algebra library installed.
   A fallback (see next point) is available inside ABINIT (basic version of lapack),
   but you might want to install a math library yourself, especially for parallel computations:
   You can choose among:

     - `blas` (libblas-dev)
     - `lapack` (liblapack-dev)
     - `scalapack` (libscalapack-...-dev)
     - `atlas` (libatlas-base-dev),
     - `mkl` from Intel (or you might try the libmkl-full-dev package).

5. Some mandatory libraries:

      - HDF5, NetCDF and NetCDF-Fortran, libraries to write/read binary files in netcdf4 format.
       These libraries are available via the `libhdf5-dev`, `libnetcdf-dev` and `libnetcdff-dev` packages from apt.
       For parallel IO, the `libpnetcdf-dev` is required.

      - LIBXC, a library containing exchange-correlation potentials, from the `libxc-dev` package.

   Note that it is also possible to generate these libraries via the ABINIT fallbacks:

   ```sh
   cd fallbacks
   ./build-abinit-fallbacks.sh
   ```

   In the latter case, the ABINIT configuration file (see later) should contain `with_fallbacks="yes"`.

These are the commands required to install the required packages from apt 
assuming a relatively simple ABINIT build with MPI support.
The list of commands may change depending on your linux distribution,
the exact ABINIT version you want to compile and the libraries you want to use.

```sh
# 1 # compiler
sudo apt install gfortran

# 2 # MPI libraries - choice for Open MPI
sudo apt install mpi-default-dev libopenmpi-dev

# 3 # math libraries - choice for lapack and blas
sudo apt install liblapack-dev libblas-dev

# 4 # mandatory libraries
sudo apt install libhdf5-dev libnetcdf-dev libnetcdff-dev libpnetcdf-dev libxc-dev
```

### Compiling ABINIT

For normal users it is advised to download the latest stable version from our website
(replace 9.0.4 by the newest version available).

```sh
wget https://www.abinit.org/sites/default/files/packages/abinit-9.0.4.tar.gz
tar xzf abinit-9.0.4.tar.gz
cd abinit-9.0.4
```

Create a working directory:

```sh
mkdir build && cd build
```

To configure, use:

```sh
../configure --with-config-file='my_config_file.ac9'
```

where 'my_config_file.ac9' is a configuration file that is discussed in more details in the next section.

Compile with:

```sh
make -j4
```

Install (optional):

```sh
make install
```

### The configuration file

The configure script accepts variables and flags to customize the configuration.
For example

```sh
../configure FC=mpifort --with-mpi="yes"
```

tells ABINIT to configure for a MPI build with the *mpifort* MPI wrapper.
To obtain the documentation for the different variables and flags, use:

```sh
../configure --help
```

Most configuration options are detected automatically by configure.
For example, if `with_mpi` is set to 'yes', configure will try to use the parallel fortran compiler (mpifort)
and automatically detect the MPI installation with libraries (.so) and header (.h) files.
When you install the Open MPI package via apt, these directories can be printed to terminal 
by using `dpkg -L 'libopenmpi-dev'`.

When a lot of options must be passed to configure, it is advised to use an external configuration file with
the syntax:

```sh
../configure --with-config-file='my_config_file.ac9'
```

An example of configuration file to build abinit with MPI, lapack and blas and automatic 
detection for libxc, hdf5, and netcdf:

```sh
# MPI settings
with_mpi="yes"
enable_mpi_io="yes"

# linear algebra settings
with_linalg_flavor="netlib"
LINALG_LIBS="-L/usr/lib/x86_64-linux-gnu -llapack -lblas"

# mandatory libraries
with_libxc="yes"
with_hdf5="yes"
with_netcdf="yes"
with_netcdf_fortran="yes"
```

Note that:

  - one uses '-' when typing a flag but '_' inside the config file, e.g. `--with-mpi="yes"` becomes `with_mpi="yes"`.

  - the LINALG_LIBS variable was explicitly set for this linux distrubution.
    The directory was extracted via `dpkg -L liblapack-dev` and `dpkg -L libblas-dev`.

  - when fine tuning variables and flags for a particular linux distribution, it is advised to
    take a look at the template file `~abinit/doc/build/config-template.ac9`.
    For example, the setting of `LINALG_LIBS` in this template file is given
    by the line `#LINALG_LIBS="-L/usr/local/lib -llapack -lblas"`.

More specialized libraries might be harder to detect.
For example, the following section was added to the config file to detect customized FFT and XML libraries.
These libraries are available via apt (`libfftw3-dev `and `libxml2-dev`).
The directories for the corresponding library and header files can be found by using `dpkg -L [package]`
and other flags can be extracted from the `~abinit/doc/build/config-template.ac9` template

```sh
# fast fourier settings
with_fft_flavor="fftw3"
FFTW3_CPPFLAGS="-I/usr/include"
FFTW3_FCFLAGS="-I/usr/include"
FFTW3_LIBS="-L/usr/lib/x86_64-linux-gnu -lfftw3 -lfftw3f"

# XML2 library (used by multibinit)
with_libxml2="yes"
LIBXML2_FCFLAGS="-I/usr/lib/x86_64-linux-gnu/"
LIBXML2_LIBS="-L/usr/include/libxml2/libxml/ -lxmlf90"
```
