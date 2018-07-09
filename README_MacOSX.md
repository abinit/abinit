---
authors: M. Torrent
---

# How to install ABINIT on Mac OSX

This file describes how to install ABINIT on Mac OS X:

 - Using the Homebrew package manager
 - Using the MacPorts package manager
 - Compiling from source 

## Using [homebrew](http://brew.sh)

Tested with mac OS X v10.9 (Mavericks), v10.10 (Yosemite), 10.11 (El Capitan).
A Homebrew official formula for ABINIT is available (authors: D. Davydov + M. Torrent).

### Prerequesites

1. Homebrew installed (see: <http://brew.sh/#install>)
2. Homebrew needs Xcode and "Xcode command line tools" to be installed. 
   Just type:

        xcode-select --install

### Installing ABINIT

Just type:

    brew install homebrew/science/abinit

ABINIT should install smoothly with its dependencies.

Note:

* At the first installation, this can take time, because of the
  required compilation of some dependencies (netcdf, fftw ...).

* LibXC, ETSF_IO and netCDF fallbacks (plugins) are used by default.
  FoX, Wannier90 and BigDFT are not available in Homebrew.
  AtomPAW can be installed as a separate package (formula).

* The following options are available for ABINIT formula:

      - --with-testsuite --> Run full test suite (time consuming)
      - --without-check --> Skip build-time tests (not recommended)
      - --without-etsf_io --> Build without etsf_io support
      - --without-netcdf --> Build without netcdf support
      - --without-fftw --> Build without fftw support
      - --without-gsl --> Build without gsl support
      - --without-scalapack --> Build without scalapack support

## Using [macports](http://www.macports.org)

Available from ABINIT v7.4.3
Tested with mac OS X v10.8 (Mountain Lion), v10.9 (Mavericks), 
v10.10 (Yosemite) v10.11 (El Capitan)

Prerequesites:

* MacPorts installed (see <https://www.macports.org/install.php>)

* Some basic ports already installed:

    1. gcc (last version) with Fortran variant (Fortran compiler),
    2. mpich or openmpi (MPI)

* Before starting, it is preferable to update MacPorts system:

        sudo port selfupdate
        sudo port upgrade outdated

Notes:

  - It is recommended to completely reinstall MacPorts after a MacOS upgrade.
  - MacPorts needs Xcode and "Xcode command line tools" to be installed; just type:

        xcode-select --install

##  Using the official ABINIT port

There are available ports in the MacPorts system for some ABINIT versions.
They have been originally created by ABINIT developers.
But they are not necessarily maintained for all ABINIT versions.

Try first to install official ABINIT port:

    sudo port install abinit

If the lastest ABINIT is successfully installed, it's OK for you.
If not, follow the procedure below.


## Using ABINIT port from ftp.abinit.org

### Preparing your MacPorts installation

You need to create a local repository designed to receive local ports.
Let's name it LOCAL_REPOSITORY in the following. This might be:

    LOCAL_REPOSITORY=/Users/my_login/ports   
    
or:

    /opt/local/localports

Replace LOCAL_REPOSITORY by the correct string in the following...

1. Edit the file /opt/local/etc/macports/sources.conf (with administrator privileges)
   and add these two lines at the end:

        file://LOCAL_REPOSITORY
        rsync://rsync.macports.org/release/ports [default]

2. Create the directory: LOCAL_REPOSITORY/science/abinit
  
### Installing ABINIT

Download the "portfile" corresponding to your ABINIT version
and copy it in the LOCAL_REPOSITORY/science/abinit directory.
The portfile for ABINIT vX.Y.Z is available at: http://ftp.abinit.org/MacPorts/X.Y.Z/Portfile.tar.gz

Untar the file in the local repository:

    cd LOCAL_REPOSITORY/science/abinit
    tar -xvzf Portfile.tar.gz

Update your local repository index:

    cd LOCAL_REPOSITORY
    sudo portindex

Install ABINIT port:

    sudo port install abinit @X.Y.Z

### ABINIT port variants

By default, ABINIT is installed with the following plugins/fallbacks:

    libXC, ETSF_IO, Wannier90

Linking ABINIT to FFTW3 library:

     sudo port install abinit @X.Y.Z +fftw3

Linking ABINIT to parallel Linear Algebra ScaLapack:

     sudo port install abinit @X.Y.Z +scalapack

Installing a multi-threaded (openMP) version of ABINIT:

     sudo port install abinit @X.Y.Z +threads

Installing atompaw PAW atomic dataset generator:

     sudo port install abinit @X.Y.Z +atompaw
or 

    sudo port install atompaw

It is possible to mix all previous variants: 

    sudo port install abinit @X.Y.Z +fftw3+threads+scalapack+atompaw

Other options available by typing:

    port info abinit

## Compiling from source under MacOSX

Prerequesites

1. Mac OSX

2. Xcode installed with "Xcode command line tools"; just type:

        xcode-select --install

3. A Fortran compiler installed. Possible options:

      - gfortran binary from: http://hpc.sourceforge.net
      - gfortran binary from: https://gcc.gnu.org/wiki/GFortranBinaries#MacOS
      - gfortran installed via a package manager (MacPorts, Homebrew, Fink)
      - Intel Fortran compiler

4. A MPI library installed  (If you want to benefit from parallelism; recommended).
   Possible options:

      - mpich from http://www.mpich.org, or via a package manager
      - openmpi from http://www.open-mpi.org, or via a package manager

5. A Linear Algebra library installed.
  By default the 'accelerate' Framework is installed on MacOSX
  and ABINIT build system should find it.
  But you might want to install a parallel library: scalapack, atlas, mkl, ...

### Installing ABINIT

Create a working directory:

    cd abinit_src_dir
    mkdir build && cd build

To configure the sequential version:

    ../configure FC=gfortran CC=clang FCFLAGS_EXTRA="-ffree-line-length-none"

for the parallel version (only with MPI installed):

    ../configure FC=mpif90 CC=mpicc FCFLAGS_EXTRA="-ffree-line-length-none" \
        --enable-mpi  --enable-mpi-io

Compile with:
    
    make mj4

Install (optional):

    make install

## Comments

To benefit from the "fallbacks" (libXC, netCDF, ...), consult the configure
script help: 

    ../configure --help

The --with-dft-flavor and --with-trio-flavor have to be adjusted.
