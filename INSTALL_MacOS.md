---
authors: M. Torrent
---

# How to install ABINIT on macOS

This file describes how to install ABINIT on macOS using one of the following approaches:

 - Homebrew package manager
 - MacPorts package manager
 - Compilation from source
 - Example : MacPorts + ABINIT source

## Using [homebrew](http://brew.sh)

Tested with:

1. macOS 10.13 (High Sierra)
2. macOS 10.14 (Mojave)
3. macOS 10.15 (Catalina)

A Homebrew official formula for ABINIT is available in this [github repository](https://github.com/abinit/homebrew-tap).

### Prerequisites

- Homebrew installed (see: <http://brew.sh/#install>)

### Installing ABINIT

Before the first installation, type:

    brew tap abinit/tap

To install ABINIT just type

    brew install abinit

and ABINIT should install smoothly with all its dependencies.

Note that:

* the LibXC and netCDF fallbacks are enabled by default.
  Wannier90 and BigDFT are not available in Homebrew.
  AtomPAW can be installed as a separate formula.

* the following extra options are available for the ABINIT formula:

      * --with-testsuite    --> Run full test suite (time consuming)
      * --without-check     --> Skip build-time tests (not recommended)
      * --without-openmp    --> Disable openMP multithreading
      * --without-netcdf    --> Build without netcdf support
      * --without-libxc     --> Build without libXC support
      * --without-fftw      --> Build without fftw support
      * --without-scalapack --> Build without scalapack support

## Using [macports](http://www.macports.org)

ABINIT is available on the MacPorts project, but it is not necessarily the latest version.
The procedure has been tested with Mac OS X v10.8 (Mountain Lion) and v10.15 (Catalina)

### Prerequisites:

1. MacPorts installed (see <https://www.macports.org/install.php>)

2. Some basic ports already installed:

    1. gcc (last version) with Fortran variant (Fortran compiler),
    2. mpich or openmpi (MPI)

3. Before starting, it is preferable to update the MacPorts system:

        sudo port selfupdate
        sudo port upgrade outdated

### Installing ABINIT

To install ABINIT, just type:

    sudo port install abinit

### ABINIT port variants

By default, ABINIT is installed with libXC and Wannier90

To activate support for the FFTW3 library:

     sudo port install abinit @X.Y.Z +fftw3

where X.Y.X is to be replaced by the version of ABINIT that has been installed (you might skip X.Y.Z if there is only one version installed).

To link ABINIT with the parallel Linear Algebra ScaLapack:

     sudo port install abinit @X.Y.Z +scalapack

To install a multi-threaded (openMP) version of ABINIT:

     sudo port install abinit @X.Y.Z +threads

It is possible to mix all previous variants:

    sudo port install abinit @X.Y.Z +fftw3+threads+scalapack

Other options available, see:

    port info abinit

## Compiling from source under macOS

### Prerequisites

1. macOS

2. Xcode installed with "Xcode command line tools"; just type:

        xcode-select --install

3. A Fortran compiler installed. Possible options:

      - [gfortran-for-macOS project](https://github.com/fxcoudert/gfortran-for-macOS/releases).
      - gfortran binary from [http://hpc.sourceforge.net](http://hpc.sourceforge.net)
      - gfortran binary from [https://gcc.gnu.org/wiki/GFortranBinaries#MacOS](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS)
      - gfortran installed via a package manager (MacPorts, Homebrew, Fink)
      - intel Fortran compiler

4. Mandatory libraries.

      - [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
        (High-performance data management and storage suite)

      - [NetCDF](https://www.unidata.ucar.edu/software/netcdf/)
        (Network Common Data Form) 
      - [libXC](https://tddft.org/programs/libxc/download/) (library of exchange-correlation functionals)

5. A MPI library installed (if you want to benefit from parallelism; recommended).
   Possible options:

      - mpich from [http://www.mpich.org](http://www.mpich.org), or via package manager
      - open-mpi from [http://www.open-mpi.org](http://www.open-mpi.org), or via package manager

6. A Linear Algebra library installed.
   By default the `accelerate` Framework is installed on macOS and the ABINIT build system should find it.
   But you might want to install a parallel library: `scalapack`, `atlas`, `mkl`, etc.
   If ABINIT is linked with the accelerate library, make sure the code is configured with
   `--enable-zdot-bugfix="yes"`. 
   Use `otool -L abinit` to print the shared libraries required by the program.

### Installing ABINIT from source

For normal users it is advised to get the newest version from our website (replace 9.4.2 by the newest version available).

    wget https://www.abinit.org/sites/default/files/packages/abinit-9.4.2.tar.gz
    tar xzf abinit-9.4.2.tar.gz
    cd abinit-9.4.2

Create a working directory:

    mkdir build && cd build

To configure the sequential version, use:

    ../configure FC=gfortran CC=clang FCFLAGS_EXTRA="-ffree-line-length-none"

For the parallel version (only if MPI installed):

    ../configure FC=mpif90 CC=mpicc FCFLAGS_EXTRA="-ffree-line-length-none" \
                 --enable-mpi --enable-mpi-io

For the GNU compiler version 10 or beyond, "-ffree-line-length-none" needs to be replaced
by "-ffree-line-length-none -fallow-argument-mismatch".

Compile with:

    make -j4

Install (optional):

    make install

Remember that on MacOs, the environment variant LD_LIBRARY_PATH should be replaced by DYLD_LIBRARY_PATH.

## Example : [macports](http://www.macports.org) + ABINIT source

We give an example of using ABINIT sources (to get the latest ABINIT version) combined with the ease of installation of Macports.
Indeed, the libraries (mandatory or optional) used by ABINIT are available from the MacPorts project.
The procedure has been tested with Mac OS X v11.5 (Big Sur) 

Start with the prerequisite of the macports approach, as explained in the section [using macports](#using-macports). Then,
install the libraries needed by ABINIT. There are different choices for the linear algebra and mpi,
as mentioned in the section [Compiling from source](#compiling-from-source-under-macos). Focusing on the choice
of OpenBLAS and openmpi, and supposing that GNU compiler version 11 is to be used, one might issue 

     sudo port install gcc11
     sudo port install OpenBLAS +gcc11+fortran
     sudo port install openmpi-gcc11 +gfortran
     sudo port install fftw-3 +gfortran
     sudo port install fftw-3-single +gfortran
     sudo port install fftw-3-long +gfortran
     sudo port install hdf5 +cxx+gcc11+hl+openmpi
     sudo port install netcdf +cdf5+dap+x+gcc11+netcdf4+openmpi 
     sudo port install netcdf-fortran +gcc11+openmpi 
     sudo port install libxc4 +gcc11

and optionally

     sudo port install wannier90 +accelerate+gcc11
     sudo port install atompaw +accelerate+gcc11+libxc

After this step, one can follow the steps described in the subsection "Installing ABINIT from source" of
[Compiling from source](#compiling-from-source-under-macos).

## Troubleshooting

When switching from GNU compiler version 10 to GNU compiler version 11 on the same machine, with MacOs Big Sur v11.5, it has been seen that Xcode needs to be reinstalled.
Indeed, at the level of the abinit compilation, the system library (-lSystem) was not found, which was quite hard to debug.
Actually, the reinstallation of Xcode might be good practice to solve other problems.
To do this, use

    sudo rm -rf /Library/Developer/CommandLineTools
    sudo xcode-select --install

You might have to accept the licence on the screen. Also, possibly first using
 
    sudo xcodebuild -license


## Comments

To benefit from the optional "fallbacks" (`Wannier90`, `libPSML`, ...),
consult the [abinit-fallbacks Project](https://gitlab.abinit.org/buildbot/abinit-fallbacks)
