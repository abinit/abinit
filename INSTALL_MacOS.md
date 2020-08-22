---
authors: M. Torrent
---

# How to install ABINIT on Mac OS

This file describes how to install ABINIT on MacOS:

 - Using the Homebrew package manager
 - Using the MacPorts package manager
 - Compiling from source 

## Using [homebrew](http://brew.sh)

Tested with:

    1. mac OS X v10.12 (Sierra)
    2. 10.13 (High Sierra)
    3. 10.14 (Mojave)

A Homebrew official formula for ABINIT is available in this [github repository](https://github.com/abinit/homebrew-tap).

### Prerequisites

1. Homebrew installed (see: <http://brew.sh/#install>)

### Installing ABINIT

Before the first installation, type:

    brew tap abinit/tap

To install ABINIT just type:

    brew install abinit

and ABINIT should install smoothly with its dependencies.

Note:

* the LibXC and netCDF fallbacks (plugins) are used by default.
  Wannier90 and BigDFT are not available in Homebrew.
  AtomPAW can be installed as a separate package (formula).

* The following options are available for the ABINIT formula:

      - --with-testsuite    --> Run full test suite (time consuming)
      - --without-check     --> Skip build-time tests (not recommended)
      - --without-openmp    --> Disable openMP multithreading
      - --without-netcdf    --> Build without netcdf support
      - --without-libxc     --> Build without libXC support
      - --without-fftw      --> Build without fftw support
      - --without-scalapack --> Build without scalapack support

## Using [macports](http://www.macports.org)

ABINIT is available on the MacPorts project, not necessarily in its latest version.
Tested with mac OS X v10.8 (Mountain Lion) --> v10.14(Mojave)

### Prerequisites:

1. MacPorts installed (see <https://www.macports.org/install.php>)

2. Some basic ports already installed:

    1. gcc (last version) with Fortran variant (Fortran compiler),
    2. mpich or openmpi (MPI)

3. Before starting, it is preferable to update MacPorts system:

        sudo port selfupdate
        sudo port upgrade outdated

### Installing ABINIT

To install ABINIT just type:

    sudo port install abinit

### ABINIT port variants

By default, ABINIT is installed with the following plugins/fallbacks:

    libXC, Wannier90

To activate support for the FFTW3 library:

     sudo port install abinit @X.Y.Z +fftw3

To link ABINIT with the parallel Linear Algebra ScaLapack:

     sudo port install abinit @X.Y.Z +scalapack

To install a multi-threaded (openMP) version of ABINIT:

     sudo port install abinit @X.Y.Z +threads

It is possible to mix all previous variants: 

    sudo port install abinit @X.Y.Z +fftw3+threads+scalapack

Other options available, see:

    port info abinit

## Compiling from source under MacOS

### Prerequisites

1. Mac OS

2. Xcode installed with "Xcode command line tools"; just type:

        xcode-select --install

3. A Fortran compiler installed. Possible options:

      - [gfortran-for-macOS project](https://github.com/fxcoudert/gfortran-for-macOS/releases).
      - gfortran binary from: [http://hpc.sourceforge.net](http://hpc.sourceforge.net)
      - gfortran binary from: [https://gcc.gnu.org/wiki/GFortranBinaries#MacOS](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS)
      - gfortran installed via a package manager (MacPorts, Homebrew, Fink)
      - Intel Fortran compiler

4. A MPI library installed (if you want to benefit from parallelism; recommended).
   Possible options:

      - mpich from [http://www.mpich.org](http://www.mpich.org), or via a package manager
      - openmpi from [http://www.open-mpi.org](http://www.open-mpi.org), or via a package manager

5. A Linear Algebra library installed.
  By default the `accelerate` Framework is installed on MacOS
  and ABINIT build system should find it.
  But you might want to install a parallel library: `scalapack`, `atlas`, `mkl`, ...

### Installing ABINIT

Download ABINIT. 
For normal users it is advised to get the newest version from our website (replace 9.0.4 by the newest version available).

    wget https://www.abinit.org/sites/default/files/packages/abinit-9.0.4.tar.gz
    tar xzf abinit-9.0.4.tar.gz
    cd abinit-9.0.4

Create a working directory:

    mkdir build && cd build

To configure the sequential version:

    ../configure FC=gfortran CC=clang FCFLAGS_EXTRA="-ffree-line-length-none"

For the parallel version (only with MPI installed):

    ../configure FC=mpif90 CC=mpicc FCFLAGS_EXTRA="-ffree-line-length-none" \
        --enable-mpi  --enable-mpi-io

Compile with:

    make -j4

Install (optional):

    make install

## Comments

To benefit from the "fallbacks" (`libXC`, `netCDF`, ...), consult the configure script help: 

    ../configure --help

The `--with-dft-flavor` and `--with-trio-flavor` have to be adjusted.
