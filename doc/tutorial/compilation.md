---
authors: MG
---

# How to compile ABINIT

This tutorial explains how to compile ABINIT and the required dependencies
without relying on pre-compiled libraries, package managers and root privileges.
You will learn how to use the standard approach based on **configure** and **GNU make** 
to build and install your own software stack.
It is assumed that you already have a standard Linux installation 
that already provides some basic tools to build software from source.
The changes required for MacOsX are briefly mentioned when needed.
In the last part of the tutorial, we briefly explain how to use modules 
to load the external dependencies, how to link against the intel MKL library and, optionally, 
how to activate support for  OpenMP threads. 

<!--
We will make extensive use of the shell hence familiarity with the Unix terminal is assumed.
Actually, you are supposed to perform all the steps of this tutorial without any interaction with the mouse or
any fancy drag-and-drop interface because the environment typically found on super-computing centers is not so
user-friendly and it's better that you get used to it.
Last but not least, forget about `sudo` and installation in system-wide directiories such as */usr/local*.
Our goal is to teach you how to create an environment to compile Abinit with minimal privileges.
A list of useful commands and options is also available 
[here](https://www3.ntu.edu.sg/home/ehchua/programming/cpp/gcc_make.html)
-->

If this is the first time you use "configure && make" to build software, 
we strongly recommended to read this 
[quick introduction](https://www.codecoffee.com/software-installation-configure-make-install/)
before proceeding with the next steps.

<!--
## Useful concepts
### Semantic version
### Dynamic linkage
Before starting, it is worth checking whether your Linux installation provides some fundamental tools 
needed to build application from source.
-->

## Getting started

Since ABINIT is written in Fortran, we need a modern Fortran compiler that supports the **F2003 specifications**
as well as a C compiler (the C++ compiler is optional and required for advanced features such as DMFT).
In what follows, we will be focusing on the GNU toolchain: *gcc* for C and *gfortran* for Fortran.

In the terminal, issue:

```sh
gfortran --version
```

to get the version of the compiler (must be >= ??).

!!! tip

    If gfortran is not installed, use the package manager provided by your Linux distribution to install it.
    On Ubuntu, use:

    ```sh
    sudo apt-get install gfortran
    ```

    Life is much more complicated if you are a Mac-OsX users since Mac does not officialy 
    support Fortran and you will need to install a Fortran compiler.


The *which* command, reports the absolute path of the executable.

```sh
which gfortran
```

This point is very important as we will use the *which* command as lot in the rest of this tutorial 
to find the location of the different dependencies.

Now check whether GNU **make** is installed with

```sh
which make
```

MaxOsX users will need to install **make** via Xcode.

## How to compile BLAS and LAPACK

BLAS/LAPACK represents the workhorse of any scientific code.
<!--
In principle this step can be skipped since any decent Linux distribution already provides
pre-compiled versions but, as already mentioned in the introduction, we prefer 
to compile everything from source and BLAS/LAPACK represents an excellent exercise.
-->
The compilation of these two libraries is indeed straightforward and this allows us 
to introduce some basic concepts that will reveal useful in the other parts of the tutorial.

First of all, let's create a directory in your `$HOME` directory (let's call it *local*).
The *src* directory will be used to store the different source files while artifacts 
will be installed in `$HOME/local/bin`, `$HOME/local/lib` etc.

Let's create the *local* directory with:

```sh
cd $HOME && mkdir local
```

!!! tip

    "&&" is used to chain commands together, such that the next command is executed if and only 
    if the preceding command exited without errors (or, more accurately, exits with a return code of 0).

Now create the `src` subirectory with:

```sh
cd local && mkdir src && cd src
```

Download the openblas tarball from the [openblas website](https://www.openblas.net/) with:

```sh
wget https://github.com/xianyi/OpenBLAS/archive/v0.3.7.tar.gz
```

If `wget` is not available, use

```
curl
```

Now unpack the tarball with:

```sh
tar -xvf v0.3.7.tar.gz
```

then cd to the directory with

```sh
cd OpenBLAS-0.3.7/
```

and execute the *configure* script with:

```sh
./configure --prefix=$HOME/local/
make -j2
```

There are two things worth noticing here:

* The `--prefix` option specifies the location where all the libraries, executables, include files, modules, 
  documentation pages, etc will be installed when `make install` is issued.
* The `-j2` option tells make to use 2 processes to build the package. 
  Adjust this value according to the number of CPUs available on your machine.

At the end of make, you should get the following output:

```md
 OpenBLAS build complete. (BLAS CBLAS LAPACK LAPACKE)

  OS               ... Linux
  Architecture     ... x86_64
  BINARY           ... 64bit
  C compiler       ... GCC  (command line : cc)
  Fortran compiler ... GFORTRAN  (command line : gfortran)
  Library Name     ... libopenblas_haswellp-r0.3.7.a (Multi threaded; Max num-threads is 12)

To install the library, you can run "make PREFIX=/path/to/your/installation install".
```

Before installing the library, it is good common practice to run the test suite to validate the build.
Many packages allows users to run the test suite by issuing:

```sh
make tests
```

other packages provide `make check`. 
As usual, use `make --help` to list the available options.

```sh
make PREFIX=$HOME/local/ install
```

So far so good, we managed to compile and install our own version of BLAS/LAPACK.
At this point, you should have the following libraries in $HOME/local/lib:

```sh
ls $HOME/local/lib/

cmake          libopenblas.so    libopenblas_haswellp-r0.3.7.a   pkgconfig
libopenblas.a  libopenblas.so.0  libopenblas_haswellp-r0.3.7.so
```

<!--
Use

```sh
nm $HOME/local/lib/libopenblas.so
```

to list the symbols presented in the library.
-->

Since we have installed the package inside a non-standard location ($HOME/local), 
we need to update two important shell variables so that bash and 
Add these two lines at the end of your $HOME/.bash_profile file 

```sh
export PATH=$HOME/local/bin:$PATH
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
```

and then execute:

```sh
source $HOME/.bash_profile file
```

to activate these changes.
Now use

```sh
echo $PATH
echo $LD_LIBRARY_PATH
```

to print the value of the variables.

If this is the first time you hear about $PATH and $LD_LIBRARY_PATH, please take some time to read
about the meaning of this environment variables.
Also remember that one can use `env` to print all the environment variables.

!!! tip

    MaxOsx users should replace `LD_LIBRARY_PATH` with `DYLD_LIBRARY_PATH`


## How to compile Libxc

At this point, it should not be so difficult to compile and install the libxc library for the XC functional.
Libxc is written in C, it supports the standard `configure && make` approach and, most importantly, 
the library does not have any external dependecies other than basic libraries that are available on every decent
Unix-like installations.
The only thing worth noticing is that you will need to activate the Fortran interface

Also in this case, you are supposed to configure the package with the *--prefix* option, run the tests
and finally issue `make install`.

```sh
# Get the tarball
cd $HOME/local/src
wget http://www.tddft.org/programs/libxc/down.php?file=4.3.4/libxc-4.3.4.tar.gz
tar -zxvf libxc-4.3.4.tar.gz

# Configure
cd libxc-4.3.4
./configure --prefix=$HOME/local

# Build + Tests + Installation
make -j2
make check
make install
```

Let's have a look at the libraires we have just installed:

```sh
[gmatteo@bob libxc-4.3.4]$ ls ~/local/lib/libxc*
/home/gmatteo/local/lib/libxc.a   /home/gmatteo/local/lib/libxcf03.a   /home/gmatteo/local/lib/libxcf90.a
/home/gmatteo/local/lib/libxc.la  /home/gmatteo/local/lib/libxcf03.la  /home/gmatteo/local/lib/libxcf90.la
```

## Installing FFTW

FFTW is a C library for computing the Fast Fourier transform in one or more dimensions.
Abinit alread provides an internal implementation of the FFT algorithm in Fortran 
so FFTW is an optional dependency.
If you care about performance, FFTW (or MKL) is a nice thing to have because it is usually faster 
than the internal version.

The FFTW source code can be dowloaded from [fftw.org](http://www.fftw.org/), 
the tarball of the latest version is available at this URL: <http://www.fftw.org/fftw-3.3.8.tar.gz>.

The compilation procedure is very similar to the one previously used for libxc. 
Note, however, that Abinit needs both the **single-precision** and the **double-precision** version.
This means that you need to configure, build and install the package twice:

At the end, you should have the following libraries installed 

```sh
[gmatteo@bob libxc-4.3.4]$ ls ~/local/lib/libxc*
/home/gmatteo/local/lib/libxc.a   /home/gmatteo/local/lib/libxcf03.a   /home/gmatteo/local/lib/libxcf90.a
/home/gmatteo/local/lib/libxc.la  /home/gmatteo/local/lib/libxcf03.la  /home/gmatteo/local/lib/libxcf90.la
```

## Installing MPI

<!--
This sections is optional although highly recommended as in 99.9% of the cases, 
you will be using a parallel version of ABINIT to perform some serious calculation.
-->

First of all, keep in mind that there are several MPI implementations available 
(openmpi, mpich, intel mpi, etc) so you must choose one and stick to it.
Each MPI library provides a `mpif90` script wrapping the Fortran compiler 
and a `mpirun` (`mpiexec`) launcher to run applications with multiple processes.

All the libraries (executables) requiring MPI should be compiled, linked and executed 
with the same MPI library.
In other words, don't try to link a library compiled with mpich if you are building the source code with openmpi 
and don't try to execute an executable compiled with intel mpi with the mpirun provided by the openmpi implementation!
Also in this case, the `which` command allows you

In this tutorial, we use the mpich implementation 
that can be downloaded from [this page](https://www.mpich.org/downloads/).

```sh
wget http://www.mpich.org/static/downloads/3.3.2/mpich-3.3.2.tar.gz
tar -zxvf mpich-3.3.2.tar.gz
cd mpich-3.3.2/
./configure --prefix=$HOME/local

make -j4
make check
make install
```

```sh
[gmatteo@bob mpich-3.3.2]$ ls $HOME/local/bin/mpi*
/home/gmatteo/local/bin/mpic++        /home/gmatteo/local/bin/mpiexec        /home/gmatteo/local/bin/mpifort
/home/gmatteo/local/bin/mpicc         /home/gmatteo/local/bin/mpiexec.hydra  /home/gmatteo/local/bin/mpirun
/home/gmatteo/local/bin/mpichversion  /home/gmatteo/local/bin/mpif77         /home/gmatteo/local/bin/mpivars
/home/gmatteo/local/bin/mpicxx        /home/gmatteo/local/bin/mpif90


[gmatteo@bob mpich-3.3.2]$ ls $HOME/local/include/mpi*
/home/gmatteo/local/include/mpi.h              /home/gmatteo/local/include/mpicxx.h
/home/gmatteo/local/include/mpi.mod            /home/gmatteo/local/include/mpif.h
/home/gmatteo/local/include/mpi_base.mod       /home/gmatteo/local/include/mpio.h
/home/gmatteo/local/include/mpi_constants.mod  /home/gmatteo/local/include/mpiof.h
/home/gmatteo/local/include/mpi_sizeofs.mod
```

!!! important

    The `.mod` files are Fortran modules generated by the Fortran compiler
    These modules will be "used" by other Fortran codes such as Abinit during the compilation process.
    Unfortunately, `.mod` files are compiler-depedent and the format may also depend on the compiler version.
    In a nutshell, one cannot use this MPI library to compile the code with another Fortran compiler.


## Installing HDF5 and netcdf4 (C and Fortran bindings)

Netcdf4 is built on top of HDF5 and consists of two different parts: the low-level layer implemented in C
and the Fortran bindings.
To create the libraries required by ABINIT, we will build these three different pieces 
in bottom-up fashion starting from HDF5.

Download the HDF5 tarball from https://www.hdfgroup.org/downloads/hdf5

Configure the package with:

```sh
    ./configure --prefix=$HOME/local/ CC=$HOME/local/bin/mpicc --enable-parallel
```


```sh
                     AM C Flags:
               Shared C Library: yes
               Static C Library: yes


                        Fortran: no

                            C++: no

                           Java: no


Features:
---------
                   Parallel HDF5: yes
Parallel Filtered Dataset Writes: yes
              Large Parallel I/O: yes
              High-level library: yes
                    Threadsafety: no
             Default API mapping: v110
  With deprecated public symbols: yes
          I/O filters (external): deflate(zlib)
                             MPE:
                      Direct VFD: no
                         dmalloc: no
  Packages w/ extra debug output: none
                     API tracing: no
            Using memory checker: no
 Memory allocation sanity checks: no
          Function stack tracing: no
       Strict file format checks: no
    Optimization instrumentation: no
```


Again make check (Warning: it may take some time!) and make install

Now let's move to netcdf.
Download the C version and the Fortran bindings from the 
[netcdf website](https://www.unidata.ucar.edu/downloads/netcdf/) with:


```sh
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-c-4.7.3.tar.gz
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.5.2.tar.gz
```

and unpack the tarball files as usual.
To compile the C library, use:

```sh
cd netcdf-c-4.7.3/
./configure --prefix=$HOME/local/ CC=$HOME/local/bin/mpicc
```

where we have used the `CC` variable to specify that we want to use the `mpicc` wrapper installed previously.
This step is important if we want to activate support for parallel IO.

```sh
# NetCDF C Configuration Summary
==============================

# General
-------
NetCDF Version:		4.7.3
Dispatch Version:       1
Configured On:		Sun Feb  9 02:43:05 CET 2020
Host System:		x86_64-pc-linux-gnu
Build Directory: 	/home/gmatteo/local/src/netcdf-c-4.7.3
Install Prefix:         /home/gmatteo/local

# Compiling Options
-----------------
C Compiler:		/home/gmatteo/local/bin/mpicc
CFLAGS:
CPPFLAGS:
LDFLAGS:
AM_CFLAGS:
AM_CPPFLAGS:
AM_LDFLAGS:
Shared Library:		yes
Static Library:		yes
Extra libraries:	-lhdf5_hl -lhdf5 -lm -ldl -lz -lcurl

# Features
--------
NetCDF-2 API:		yes
HDF4 Support:		no
HDF5 Support:		yes
NetCDF-4 API:		yes
NC-4 Parallel Support:	yes
PnetCDF Support:	no
DAP2 Support:		yes
DAP4 Support:		yes
Byte-Range Support:	no
Diskless Support:	yes
MMap Support:		no
JNA Support:		no
CDF5 Support:		yes
ERANGE Fill Support:	no
Relaxed Boundary Check:	yes
```

Now use the standard sequence of commands to compile and install HDF5:

    make -j4
    make check
    make install

Once your are done with the installation

```sh
[gmatteo@bob netcdf-c-4.7.3]$ which nc-config
~/local/bin/nc-config
```

nc-config

```sh
[gmatteo@bob netcdf-c-4.7.3]$ nc-config --all

This netCDF 4.7.3 has been built with the following features:

  --cc            -> /home/gmatteo/local/bin/mpicc
  --cflags        -> -I/home/gmatteo/local/include
  --libs          -> -L/home/gmatteo/local/lib -lnetcdf
  --static        -> -lhdf5_hl -lhdf5 -lm -ldl -lz -lcurl
```

To compile the Fortran bindings, use:

```sh
    ./configure --prefix=$HOME/local/ FC=$HOME/local/bin/mpif90
```

where `FC` points to out mpif90 wrapper.
Again:

    make -j4
    make check
    make install


nf-config


## Get the Abinit source and compile the code

From this point on, we assume that gcc, gfortran and make are installed 
and we discuss how to compile and install the hard-dependecies

There are two ways of getting the source code of ABINIT:

  * directly from the ABINIT web site ([abinit.org/](https://www.abinit.org/)) by downloading the latest production tarball;
  * from the ABINIT gitlab git repository. This is favored, as it allows easier integration and merging, testing, etc...

Once you have got the tarball, uncompress it by typing:

    tar xvzf abinit-<version> .tar.gz

where _<version>_ is the version number you downloaded, e.g. "8.6.3". Then go
into the newly-created _abinit- <version>_ directory and have a look at it.
Then answer the following questions:

Now you can try to build ABINIT. Information on how to do it is stored inside
the INSTALL file. Please read it now.

Before actually starting the compilation, type:

    ./configure --help

and read carefully the output. You might then find useful to have a look at
the template for config files stored in _~abinit/doc/build/config-template.ac9_
which will provide you with more details on the configuration. Other example
config files in that subdirectory can be used to set up your build more
quickly.  If you have a configuration file called _~/.abinit/build/hostname.ac_
(with hostname equal to the $HOSTNAME shell variable for your machine) ABINIT's
configure will load it at runtime.

The compilation will likely take more than 10 minutes. 


## Compiling Abinit on a cluster with the Intel toolchain

On intel-based machines, we suggest to compile ABINIT with the intel compilers and the intel mkl library 
The mkl package, in particular, implements highly-optimized versions of BLAS/LAPACK, FFT routines
as well as SCALAPACK.
that are usually more performant that 
usually 

To list the available modules, use:

    module avail

At this point, you should decide the toolchain to be used to compile Abinit:
the compiler e.g. intel ifort, GNU gfortran
This step is strongly cluster-dependent


### How to compile ABINIT with OpenMP support (threads)

Compiling ABINIT with OpenMP is not that difficult as everything boils down to:

* Using a threaded version for BLAS/LAPACK/FFT
* Passing the correct options and libraries to the ABINIT configure script

On the contrary, answering the questions:

* When and why should I use threads?
* How many threads should I use and what is the parallel speedup I should expect?

is more difficult because there are several factors that should be considered.

To keep a long story short, one should start to look into threads 
when we start to trigger limitations and bottlenecks in the MPI implementation, 
especially at the level of the memory requirements.
These problems are usually observed in large calculations

For a gentle introduction to threads, please consult this tutorial
