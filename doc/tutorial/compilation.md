---
authors: MG
---

# How to compile ABINIT

This tutorial explains how to compile ABINIT including the required external dependencies
without relying on pre-compiled libraries, package managers and root privileges.
You will learn how to use the standard **configure** and **make** Unix tools
to build and install your own software stack including the MPI library and the associated 
*mpif90* and *mpicc* wrappers required to compile parallel MPI applications.

It is assumed that you already have a standard Unix-like installation 
that provides the basic tools needed to build software from source (Fortran/C compilers and *make*).
The changes required for MacOsX are briefly mentioned when needed.
Windows users should install [cygwin](https://cygwin.com/index.html) that 
provides a POSIX-compatible environment 
or, alternatively, use a [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about).
Note that the procedure descrived in this tutorial has been tested using Linux/MacOsX installations so
we welcome feedback from Windows users.

!!! important

    In the last part of the tutorial, we discuss more advanced topics such as using modules in supercomputing centers, 
    compiling and linking the with the intel compiler and mkl, activating support for OpenMP threads.
    You may want to jump to this section if you are already familiar with software compilation.
<!--
    In the last part of the tutorial, we treat more advanced topics related to the usage of modules in supercomputing centers.
    We also explain how to link ABINIT with the intel MKL library and how to activate support for OpenMP threads. 
-->

Note that we will make extensive use of the bash shell hence familiarity with the terminal is assumed.
For a quick introduction to the command line, please consult 
this [Ubuntu tutorial](https://ubuntu.com/tutorials/command-line-for-beginners#1-overview).
If this is the first time you use the **configure && make** approach to build software, 
we **strongly** recommended to read this 
[guide](https://www.codecoffee.com/software-installation-configure-make-install/)
before proceeding with the next steps.

If you are not interested in compiling everything from source, you may want to consider the following alternatives:

* Internal fallbacks

* Compilation with external libraries provided by apt-get (Linux users)

* Precompiled versions provided by conda-forge (for Linux and MacOsX)

* Homebrew bottles (MacOsX)

!!! important 

    The aim of this tutorial is to teach you how to compile code from source but we cannot guarantee
    that these recipes will work out of the box on every possible architecture.
    We will do our best to explain how to setup your environment and how to avoid the typical pitfalls but 
    we cannot cover all the possible cases.
    Fortunately, the internet provides lots of resources.
    Search engines are your best friends and in some cases one can find the solution by just copying 
    the error message in the search bar.
    For more complicated situations, you can ask for help on the Abinit forum but keep in mind but 
    one should always provide enough information about the problem.


## Getting started

Since ABINIT is written in Fortran, we need a **recent** Fortran compiler 
that supports the **F2003 specifications** as well as a C compiler.
At the time of writing ( |today| ), the C++ compiler is optional and required only for advanced features 
such as the interface with the TRIQS library that won't be treated in this lesson.

In what follows, we will be focusing on the GNU toolchain i.e. *gcc* for C and *gfortran* for Fortran.
These "sequential" compilers are adequate if you don't need to compile parallel MPI applications.
The compilation of MPI code, indeed, requires the installation of **additional libraries**
and **specialized wrappers** (*mpif90*, *mpicc*) that "replace" the "sequential" compilers.
This very important case scenario is covered in more details in the next sections.
For the time being, we mainly focus on the compilation of sequential applications/libraries.

First of all, let's make sure the **gfortran** compiler is installed on your machine
by issuing in the terminal the following command:

```sh
which gfortran
/usr/bin/gfortran
```

The **which** command, returns the absolute path of the executable.
This Unix tool is extremely useful to pinpoint possible problems and we will use it 
a lot in the rest of this tutorial.

To get the version of the compiler, use the `--version` option:

```sh
gfortran --version
GNU Fortran (GCC) 5.3.1 20160406 (Red Hat 5.3.1-6)
Copyright (C) 2015 Free Software Foundation, Inc.
```

At present, ABINIT requires version >= ??.

In our case, we are lucky that the system administrator installed *gfortran* in */usr/bin* and we can start to use
the compiler to build our software stack.
If gfortran is not installed, you may want to use the package manager provided by your Linux distribution to install it.
On Ubuntu, for instance, one can use:

```sh
sudo apt-get install gfortran
```

Now let's check whether **make** is already installed with:

```sh
which make
/usr/bin/make
```

!!! tip

    Things get more complicated if you are a Mac-OsX user as Apple does not officially 
    support Fortran so you will need to install gfortran either via 
    [homebrew](https://brew.sh/) or [macport](https://www.macports.org/).
    Alternatively, one can install gfortran using one of the standalone DMG installers
    provided by the [gfortran-for-macOS project](https://github.com/fxcoudert/gfortran-for-macOS/releases).
    Note also that MaxOsX users will need to install **make** via [Xcode](https://developer.apple.com/xcode/).

## How to compile BLAS and LAPACK 

BLAS/LAPACK represent the workhorse of many scientific codes and an optimized implementation 
is crucial for achieving **good performance**.
In principle this step can be skipped as any decent Linux distribution already provides
pre-compiled versions but, as already mentioned in the introduction, we are geeks and we 
prefer to compile everything from source.
Moreover the compilation of BLAS/LAPACK represents an excellent exercise 
that gives us the opportunity to discuss some basic concepts that 
will reveal useful in the other parts of this tutorial.

First of all, let's create a new directory inside your `$HOME` (let's call it **local**) using the command:

```sh    
cd $HOME && mkdir local
```

!!! tip

    $HOME is a standard shell variable that stores the absolute path to your home directory.
    Use:

    ```sh
    echo My home directory is $HOME
    ```

    to print the value of the variable.

    The **&&** syntax is used to chain commands together, such that the next command is executed if and only 
    if the preceding command exited without errors (or, more accurately, exits with a return code of 0).
    We wil use this trick a lot in the other examples to reduce the number of lines we have to type 
    in the terminal so that one can easily cut and paste the examples.


Now create the `src` subdirectory inside $HOME/local with:

```sh
cd $HOME/local && mkdir src && cd src
```

Note that the *src* directory will be used to store the packages with the source files and compile the code, 
whereas executables and libraries will be installed in `$HOME/local/bin` and`$HOME/local/lib`, respectively.

Download the tarball from the [openblas website](https://www.openblas.net/) with:

```sh
wget https://github.com/xianyi/OpenBLAS/archive/v0.3.7.tar.gz
```

If *wget* is not available, use *curl* with the `-o` option to specify the name of the output file as in:

```sh
curl https://github.com/xianyi/OpenBLAS/archive/v0.3.7.tar.gz -o v0.3.7.tar.gz 
```

!!! tip

    To get the URL associated to a HTML link inside the browser, hover the mouse pointer over the link, 
    press the right mouse button and then select `Copy Link Address` to copy the link to the system clipboard. 
    Then paste the text in the terminal by selecting the `Copy` action in the menu
    activated by clicking on the right button.
    Alternatively, one can press the central button (mouse wheel) or use CMD + V on MacOsX.
    This trick is quite handy to fetch tarballs from the web inside the terminal.


Now uncompress the tarball with:

```sh
tar -xvf v0.3.7.tar.gz
```

then `cd` to the directory with:

```sh
cd OpenBLAS-0.3.7
```

and execute 

```sh
make -j2 USE_THREAD=0 USE_LOCKING=1 
```

to build the single thread version.
By default, openblas activates threads (see [FAQ page](https://github.com/xianyi/OpenBLAS/wiki/Faq#multi-threaded))
but in our case we prefer to use the sequential version as Abinit is mainly optimized for MPI.

<!--
OpenMP threads are supported by Abinit but, for the time being, we prefer to focus 
and execute the *configure* script with:

```sh
./configure --prefix=$HOME/local/
```

where we used the `--prefix` option to specify the location where all the libraries, executables, include files, 
Fortran modules, man pages, etc will be installed when we execute `make install`
This point is very important because the default is `/usr/local`, a directory
in which normal users don't have write privileges. 
This means that that you won't be able to install software in this directory
unless you perform this step as root with *sudo*. 
The reason why we are installing in $HOME/local/ is that we want to keep our software stack well separated 
from to the libraries installed by your Linux distribution so that we can easily test new libraries and/or 
different versions without affecting the OS installation.

Now issue:
-->
The `-j2` option tells *make* to use 2 processes to build the code in order to speed up the compilation. 
Adjust this value according to the number of **physical cores** available on your machine.

At the end, you should get the following output:

```md
 OpenBLAS build complete. (BLAS CBLAS LAPACK LAPACKE)

  OS               ... Linux
  Architecture     ... x86_64
  BINARY           ... 64bit
  C compiler       ... GCC  (command line : cc)
  Fortran compiler ... GFORTRAN  (command line : gfortran)
  Library Name     ... libopenblas_haswell-r0.3.7.a (Single threaded)

To install the library, you can run "make PREFIX=/path/to/your/installation install".
```

A compilation with plain *make* would give:

```md
  Library Name     ... libopenblas_haswellp-r0.3.7.a (Multi threaded; Max num-threads is 12)
```

that indicates that our library supports threads.

You may have noticed that, in this particular case, *make* is not just building the library but is also 
running unit tests to validate the build.
This means that if *make* completes successfully, we can be confident that the build is OK 
and we can proceed  with the installation.
Other packages use a different philosophy and provide a `make check` option that should be executed after *make* 
in order to run the test suite before instaling the package.

To install openblas in $HOME/local, issue:

```sh
make PREFIX=$HOME/local/ install
```

At this point, one should have the following libraries installed in $HOME/local/lib:

```sh
ls $HOME/local/lib/libopenblas*

/home/gmatteo/local/lib/libopenblas.a     /home/gmatteo/local/lib/libopenblas_haswell-r0.3.7.a
/home/gmatteo/local/lib/libopenblas.so    /home/gmatteo/local/lib/libopenblas_haswell-r0.3.7.so
/home/gmatteo/local/lib/libopenblas.so.0
```

Files ending with `.so` are shared libraries (`.so` stands for **shared object**) whereas 
`.a` files are static libraries.
In our examples, we will mainly use dynamic linking as this is the most common scenario.

<!--
but this also means that we need to set the value of the LD_LIBRARY_PATH environment variable.
So far so good, we managed to compile and install our version of BLAS/LAPACK.
Now use

```sh
nm $HOME/local/lib/libopenblas.so
```

to list the symbols presented in the library.
-->

Since we have installed the package inside a non-standard directory ($HOME/local), 
we need to update two important shell variables: `$PATH` and `$LD_LIBRARY_PATH`.

Add these two lines at the end of your `$HOME/.bash_profile` file 

```sh
export PATH=$HOME/local/bin:$PATH
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH

# This for OpenMP
export OMP_NUM_THREADS=1
```

then execute:

```sh
source $HOME/.bash_profile
```

to activate these changes without having to start a new terminal session.

!!! important

    Remember the 0-th rule in software developement: if something does not work as expected 
    and you have followed all the instructions step by step, press restart!

Now use:

```sh
echo $PATH
echo $LD_LIBRARY_PATH
```

to print the value of these variables to the terminal.
In my case, I get:

```sh
echo $PATH
/home/gmatteo/local/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin

echo $LD_LIBRARY_PATH
/home/gmatteo/local/lib:
```

If this is the first time you hear about $PATH and $LD_LIBRARY_PATH, please take some time to learn
about the meaning of these environment variables.
More information about `$PATH` is available [here](http://www.linfo.org/path_env_var.html)

!!! tip

    MaxOsx users should replace `LD_LIBRARY_PATH` with `DYLD_LIBRARY_PATH`

    Remember that one can use `env` to print all the environment variables defined 
    in your session and pipe the results to other Unix tools. 
    Try e.g.:

    ```sh
    env | grep LD_
    ```

    to print only the variables starting with **_LD**


## How to compile libxc

At this point, it should not be so difficult to compile and install the libxc library that provides 
many useful XC functionals such as MGGA and hybrid functionals.
Libxc is written in C and can be built using the standard `configure && make` approach.
No external dependency is needed, except for basic C libraries that are available
on any decent Linux distribution.

Also in this case, we are supposed to configure the package using the *--prefix* option, 
run the tests to validate the build and finally execute `make install`.
The required commands are reported below:

```sh
# Get the tarball. 
# Mind the -O option used in wget to specify the name of the output file

cd $HOME/local/src
wget http://www.tddft.org/programs/libxc/down.php?file=4.3.4/libxc-4.3.4.tar.gz -O libxc.tar.gz
tar -zxvf libxc.tar.gz
```

Configure the package with:

```sh
cd libxc-4.3.4 && ./configure --prefix=$HOME/local
```

Finally, build, run the tests and install libxc with:

```sh
make -j2
make check && make install
```

Now let's have a look at the libraries we have just installed:

```sh
ls ~/local/lib/libxc*
/home/gmatteo/local/lib/libxc.a   /home/gmatteo/local/lib/libxcf03.a   /home/gmatteo/local/lib/libxcf90.a
/home/gmatteo/local/lib/libxc.la  /home/gmatteo/local/lib/libxcf03.la  /home/gmatteo/local/lib/libxcf90.la
```

Note the following:

  * libxc is the C library.
  * libxcf90 is the F90 library 
  * libxcf03 is the F2003 library

At present, the Fortran interface is not required by Abinit, only the C library.
We made this choice, because one can easily use the C library with Abinit compiled with different Fortran compilers/version
without having to rebuild libxc with a different Fortran compiler.
Note, however, that for other libraries (fftw3, netcdf) we will use Fortran bindings, in this case it is important that
the bindings are built with the **same compiler as the one used to compile Abinit**. 
In a nutshell, Fortran/C++ libraries are compiler-dependent because these two languages are relatively high-level so they need 
to rely on implementation details.
C libraries, on the other hand, are usually more portable as the C language has less abstractions 
and, last but not least, the Linux-OS is written in C, the Linux kernel in 99.9% of the cases 
is compiled with GNU-GCC and all the other vendors need to maintain compatibility with the GCC ABI.

In principle, one can reuse libraries as long as the major version of the Fortran compiler is the same 
but experience has  shown that it's always a good idea to require strict matching.

TODO: Discuss FCFLAGS and FCDFLAGS

## Compiling and installing FFTW

FFTW is a C library for computing the Fast Fourier transform in one or more dimensions.
ABINIT already provides an internal implementation of the FFT algorithm implemented in Fortran
hence FFTW is considered an optional dependency. 
Nevertheless, **we do not recommend the internal implementation if you really care about performance**.
The reason is that FFTW (or, even better, the DFTI library provided by intel MKL) 
is usually much faster than the internal version.

!!! important

    FFTW is very easy to install on Linux machines once you have *gcc* and *gfortran*.
    The [[fftalg]] variable defines the implementation to be used and 312 corresponds to the FFTW implementation.
    The default value of [[fftalg]] is automatically set by the *configure* script via pre-preprocessing options.
    In other words, if you activate support for FFTW (DFTI) at configure time, 
    ABINIT will use [[fftalg]] 312 (512) as default.

The FFTW source code can be downloaded from [fftw.org](http://www.fftw.org/), 
and the tarball of the latest version is available at <http://www.fftw.org/fftw-3.3.8.tar.gz>.

```
cd $HOME/local/src
wget http://www.fftw.org/fftw-3.3.8.tar.gz
tar -zxvf fftw-3.3.8.tar.gz && cd fftw-3.3.8
```

The compilation procedure is very similar to the one already used for the libxc package. 
Note, however, that ABINIT needs both the **single-precision** and the **double-precision** version.
This means that we need to configure, build and install the fftw package twice.

To build the single precision version, use:

```sh
./configure --prefix=$HOME/local --enable-single
make -j2 
make check && make install
```

Let's have a look at the libraries we've just installed with:

```sh
ls $HOME/local/lib/libfftw3*
/home/gmatteo/local/lib/libfftw3f.a  /home/gmatteo/local/lib/libfftw3f.la
```

the `f` at the end stands for `float` (C jargon for single precision)
Now we configure for the double precision version (default behaviour)

```sh
./configure --prefix=$HOME/local
make -j2 
make check && make install
```

After this step, you should have both the single (`f`) and the double (`l`) version of the library:

```sh
ls $HOME/local/lib/libfftw3*
/home/gmatteo/local/lib/libfftw3.a   /home/gmatteo/local/lib/libfftw3f.a 
/home/gmatteo/local/lib/libfftw3.la  /home/gmatteo/local/lib/libfftw3f.la
```

To make sure we have the Fortran API bundled in the library, one can use the `nm` tool
to get the list of symbols provided by the library and then use *grep* to search for the relevant API:


```sh
[gmatteo@bob fftw-3.3.8]$ nm $HOME/local/lib/libfftw3f.a | grep sfftw_plan_many_dft
0000000000000400 T sfftw_plan_many_dft_
0000000000003570 T sfftw_plan_many_dft__
0000000000001a90 T sfftw_plan_many_dft_c2r_
0000000000004c00 T sfftw_plan_many_dft_c2r__
0000000000000f60 T sfftw_plan_many_dft_r2c_
00000000000040d0 T sfftw_plan_many_dft_r2c__

[gmatteo@bob fftw-3.3.8]$ nm $HOME/local/lib/libfftw3.a | grep dfftw_plan_many_dft
0000000000000400 T dfftw_plan_many_dft_
0000000000003570 T dfftw_plan_many_dft__
0000000000001a90 T dfftw_plan_many_dft_c2r_
0000000000004c00 T dfftw_plan_many_dft_c2r__
0000000000000f60 T dfftw_plan_many_dft_r2c_
00000000000040d0 T dfftw_plan_many_dft_r2c__
```

!!! note 

    At present, there's no need to compile FFTW with MPI support because ABINIT implements its own
    version of the MPI-FFT algorithm based on the sequential FFTW version.
    The MPI-algorithm implemented in ABINIT is optimized for plane-waves codes 
    as it supports zero-padding and composite transforms for the applications of the local part of the KS potential. 


## Installing MPI

In this section, we discuss how to compile and install the MPI library.
This step is required if you want to run ABINIT with multiple processes and/or you
need to compile MPI-based libraries such as PBLAS/Scalapack or the HDF5 library with support for parallel IO.

It is worth stressing that the MPI installation provides two scripts (**mpif90** and **mpicc**)
that act as a sort of wrapper around the sequential Fortran and the C compilers, respectively.
These scripts **must be used** to compile parallel software using MPI instead 
of the "sequential" compilers e.g. `gfortran` and `gcc`. 
The MPI library also provides launcher scripts installed in the *bin* directory (*mpirun* or *mpiexec*)
that must be used to execute MPI applications with NUM_PROCS MPI processes with the syntax:

```sh
mpirun -n NUM_PROCS EXECUTABLE [ARGS]
```

Keep in mind that there are several MPI implementations available around
(e.g. *openmpi*, *mpich*, *intel mpi*, etc) so you must **choose one implementation and stick to it** 
when building your software stack.
In other words, all the libraries and executables requiring MPI must be compiled, linked and executed 
**with the same MPI library**.
Don't try to link a library compiled with e.g. *mpich* if you are building the code with 
the *mpif90* wrapper provided by e.g. *openmpi*.
By the same token, don't try to run executables compiled with e.g. *intel mpi* with the 
*mpirun* launcher provided by *openmpi*!

In this tutorial, we employ the *mpich* implementation whose source code can be downloaded 
from this [webpage](https://www.mpich.org/downloads/).
In the terminal, issue:

```sh
cd $HOME/local/src
wget http://www.mpich.org/static/downloads/3.3.2/mpich-3.3.2.tar.gz
tar -zxvf mpich-3.3.2.tar.gz
cd mpich-3.3.2/
```

to download and uncompress the tarball.
Then configure/compile/test/install the library with:

```sh
./configure --prefix=$HOME/local
make -j2
make check && make install
```

```sh
----------------------------------------------------------------------
Libraries have been installed in:
   /home/gmatteo/local/lib

If you ever happen to want to link against installed libraries
in a given directory, LIBDIR, you must either use libtool, and
specify the full pathname of the library, or use the '-LLIBDIR'
flag during linking and do at least one of the following:
   - add LIBDIR to the 'LD_LIBRARY_PATH' environment variable
     during execution
   - add LIBDIR to the 'LD_RUN_PATH' environment variable
     during linking
   - use the '-Wl,-rpath -Wl,LIBDIR' linker flag
   - have your system administrator add LIBDIR to '/etc/ld.so.conf'

See any operating system documentation about shared libraries for
more information, such as the ld(1) and ld.so(8) manual pages.
----------------------------------------------------------------------
```

Let's have a look at the MPI executables we've just installed in $HOME/local/bin:

```sh
ls $HOME/local/bin/mpi*
/home/gmatteo/local/bin/mpic++        /home/gmatteo/local/bin/mpiexec        /home/gmatteo/local/bin/mpifort
/home/gmatteo/local/bin/mpicc         /home/gmatteo/local/bin/mpiexec.hydra  /home/gmatteo/local/bin/mpirun
/home/gmatteo/local/bin/mpichversion  /home/gmatteo/local/bin/mpif77         /home/gmatteo/local/bin/mpivars
/home/gmatteo/local/bin/mpicxx        /home/gmatteo/local/bin/mpif90
```

Since we added $HOME/local/bin to our $PATH, we should see that *mpi90* is actually 
pointing to the version we have just installed:

```sh
which mpif90
~/local/bin/mpif90
```

As already mentioned, *mpif90* is a wrapper around the sequential Fortran compiler. 
To show the Fortran compiler invoked by *mpif90*, one can use:

```sh
mpif90 -v
mpifort for MPICH version 3.3.2
Using built-in specs.
COLLECT_GCC=gfortran
COLLECT_LTO_WRAPPER=/usr/libexec/gcc/x86_64-redhat-linux/5.3.1/lto-wrapper
Target: x86_64-redhat-linux
Configured with: ../configure --enable-bootstrap --enable-languages=c,c++,objc,obj-c++,fortran,ada,go,lto --prefix=/usr --mandir=/usr/share/man --infodir=/usr/share/info --with-bugurl=http://bugzilla.redhat.com/bugzilla --enable-shared --enable-threads=posix --enable-checking=release --enable-multilib --with-system-zlib --enable-__cxa_atexit --disable-libunwind-exceptions --enable-gnu-unique-object --enable-linker-build-id --with-linker-hash-style=gnu --enable-plugin --enable-initfini-array --disable-libgcj --with-isl --enable-libmpx --enable-gnu-indirect-function --with-tune=generic --with-arch_32=i686 --build=x86_64-redhat-linux
Thread model: posix
gcc version 5.3.1 20160406 (Red Hat 5.3.1-6) (GCC)
```

The C include files (*.h*) and the Fortran modules (*.mod*) have been installed in `$HOME/local/include`:

```sh
ls $HOME/local/include/mpi*
/home/gmatteo/local/include/mpi.h              /home/gmatteo/local/include/mpicxx.h
/home/gmatteo/local/include/mpi.mod            /home/gmatteo/local/include/mpif.h
/home/gmatteo/local/include/mpi_base.mod       /home/gmatteo/local/include/mpio.h
/home/gmatteo/local/include/mpi_constants.mod  /home/gmatteo/local/include/mpiof.h
/home/gmatteo/local/include/mpi_sizeofs.mod
```

!!! important

    The `.mod` files are Fortran modules produced by the Fortran compiler.
    The location of these modules must be *passed* to the Fortran compiler when compiling ABINIT.
    Note also that these `.mod` files are **compiler- and version-dependent**.
    In other words, one cannot use these `.mod` files to compile code with a different Fortran compiler.
    Moreover, you should not expect to be able to use modules compiled with a different version of the different compiler.
    This is one of the reasons why the version of the Fortran compiler used to build the software stack is very important.

## Installing HDF5 and netcdf4

In Abinit v9, HDF5 and netcdf5 have become hard-requirements.
This means that these libraries are now mandatory to compile Abinit. 
The reason is that the Abinit developers are moving away from Fortran binary files as this 
format is not portable and it is difficult to read when using high-level languages such as python. 

Netcdf4 is built on top of HDF5 and consists of two different layers: 

* The low-level C library

* The Fortran bindings i.e. Fortran routines calling the C implementation.
  This is the high-level API/ABI used by ABINIT to perform all the IO operations on netcdf files.

To build the libraries required by ABINIT, we will compile the different layers
in a bottom-up fashion starting from the HDF5 package.
Since we want to activate support for parallel IO, we need to compile the library using the wrappers
provided by our MPI installation instead of using *gcc* or *gfortran* directly.

Let's start by downloading the HDF5 tarball from this [download page](https://www.hdfgroup.org/downloads/hdf5/source-code/)
Uncompress the archive with *tar* as usual, then configure the package with:

```sh
./configure --prefix=$HOME/local/ CC=$HOME/local/bin/mpicc --enable-parallel --enable-shared
```

where we've used the *CC* variable to specify the C compiler.
This step is important in order to enable support for parallel IO.

At the end of the configuration step, you should get the following output:

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

The line with:

```sh
Parallel HDF5: yes
```

tells us that our HDF5 build will support parallel IO (because we used CC=mpicc during the configuration step).
Also note that, as far as ABINIT is concerned, Fortran support is **optional**
as ABINIT will be interfaced with HDF5 through the Fortran bindings provided by netcdf-fortran.
In other words, ABINIT requires netcdf-fortran and not the Fortran bindings for HDF5.

Again, issue `make -j NUM` followed by 
`make check` (**Warning: it may take some time!**) and finally `make install`.

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
cd netcdf-c-4.7.3
./configure --prefix=$HOME/local/ CC=$HOME/local/bin/mpicc LDFLAGS=-L$HOME/local/lib CPPFLAGS=-I$HOME/local/include
```

where we have used the `CC` variable to specify that we want to use the `mpicc` wrapper installed previously.
At the end of the configuration step, one should have

```sh
# NetCDF C Configuration Summary
==============================

# General
-------
NetCDF Version:		4.7.3
Dispatch Version:       1
Configured On:		Wed Apr  8 00:53:19 CEST 2020
Host System:		x86_64-pc-linux-gnu
Build Directory: 	/home/gmatteo/local/src/netcdf-c-4.7.3
Install Prefix:         /home/gmatteo/local

# Compiling Options
-----------------
C Compiler:		/home/gmatteo/local/bin/mpicc
CFLAGS:
CPPFLAGS:		-I/home/gmatteo/local/include
LDFLAGS:		-L/home/gmatteo/local/lib
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

Now use the standard sequence of commands to compile and install the package:

```sh
make -j2
make check && make install
```

Once your are done with the installation, use the `nc-config` executable to 
inspect the features provided by the library we've just installed.

```sh
[gmatteo@bob netcdf-c-4.7.3]$ which nc-config
~/local/bin/nc-config
```

To get a summary of the options used to build the C layer and the available features, use

```sh
[gmatteo@bob netcdf-c-4.7.3]$ nc-config --all

This netCDF 4.7.3 has been built with the following features:

  --cc            -> /home/gmatteo/local/bin/mpicc
  --cflags        -> -I/home/gmatteo/local/include
  --libs          -> -L/home/gmatteo/local/lib -lnetcdf
  --static        -> -lhdf5_hl -lhdf5 -lm -ldl -lz -lcurl
```

To compile the Fortran bindings, execute the configure script with the options:

```sh
./configure --prefix=$HOME/local/ FC=$HOME/local/bin/mpif90 LDFLAGS=-L$HOME/local/lib CPPFLAGS=-I$HOME/local/include
```

where **FC** points to our *mpif90* wrapper.
Then issue:

```sh
make -j2
make check && make install
```

To inspect the features activated in our Fortran library, use `nf-config` instead of **nc-config**:

```sh
which nf-config
```

To get a summary of the options used to build the Fortran bindings and the list of available features, use

```sh
[gmatteo@bob netcdf-c-4.7.3]$ nf-config --all
```

See also <https://www.unidata.ucar.edu/software/netcdf/docs/building_netcdf_fortran.html>

## How to compile ABINIT

In this section, we discuss how to compile and install ABINIT 
using the MPI compilers and the libraries installed previously.
First of all, download the ABINIT tarball from [this page](https://www.abinit.org/packages) using

```sh
wget https://www.abinit.org/sites/default/files/packages/abinit-9.0.2.tar.gz
```

In this case, we are using version "9.0.2" but you may want to download the
latest production version to take advantage of new features and benefit from bug fixes.
Once you got the tarball, uncompress it by typing:

```sh
tar -xvzf abinit-9.0.2.tar.gz
```

Then cd into the newly created *abinit-9.0.2* directory and take some time to read the `INSTALL` file.

Now let's try to build ABINIT from source.
Before actually starting the compilation, type:

```sh
./configure --help
```

and read carefully the documentaton of the different options. 
You might then find useful to have a look at the template for the configuration files 
stored in _~abinit/doc/build/config-template.ac9_
which will provide you with more details on the configuration. 
Other example of config files in that subdirectory can be used to set up your build more quickly.  
If you have a configuration file called _~/.abinit/build/hostname.ac_
(with hostname equal to the $HOSTNAME shell variable for your machine) ABINIT's
configure will load it at runtime.

Instead of passing options to configure from the command line, we'll be using an external file 
to gather all our options.
Note that one can use shell variables and reuse the output of external tools using
[backtick syntax](https://unix.stackexchange.com/questions/48392/understanding-backtick/48393) 
as is *\`nf-config --flibs`* to reduce the amount of typing 
and have a configuration file that can be reused in other contexts.

https://www.cprogramming.com/tutorial/shared-libraries-linux-gcc.html

In this example, we will be taking advantage of the high-level interface provided by the *with_XXX* options
to tell the build system where external dependencies are located instead of passing options explicitly.
This is the easiest approach but if configure cannot detect the dependency properly, you may need to set 
some options manually.

```sh
# -------------------------------------------------------------------------- #
# MPI support                                                                #
# -------------------------------------------------------------------------- #

# Determine whether to build parallel code (default is auto)

# Note:
#
#   * the build system expects to find subdirectories named bin/, lib/,
#     include/ under the prefix.
#
with_mpi=$HOME/local/

# Flavor of linear algebra libraries to use (default is netlib)
#
with_linalg_flavor="openblas"

# Library flags for linear algebra (default is unset)
#
LINALG_LIBS="-L$HOME/local/lib -lopenblas"

# -------------------------------------------------------------------------- #
# Optimized FFT support                                                      #
# -------------------------------------------------------------------------- #

# Flavor of FFT framework to support (default is auto)
#
with_fft_flavor="fftw3"

FFTW3_LIBS="-L$HOME/local/lib -lfftw3f -lfftw3"

# LibXC
# Website: http://www.tddft.org/programs/octopus/wiki/index.php/Libxc

# Trigger and install prefix for LibXC (default is unset)
#
#with_libxc="$HOME/local/lib"
LIBXC_LIBS="-L$HOME/local/lib -lxcf90 -lxc"
LIBXC_FCFLAGS="-L$HOME/local/include"

# NetCDF
# Website: http://www.unidata.ucar.edu/software/netcdf/

# Trigger and install prefix for NetCDF (default is unset)
#
#with_netcdf=$(nc-config --prefix)
with_netcdf="yes"

#with_netcdf_fortran=$(nf-config --prefix)
NETCDF_FORTRAN_LIBS=`nf-config --flibs`
NETCDF_FORTRAN_FCFLAGS=`nf-config --fflags`

HDF5_LIBS=`nf-config --flibs`
HDF5_FCFLAGS=`nf-config --fflags`

# Enable OpenMP (gmatteo, torrent)
enable_openmp="no"
```

!!! important

    The name of the options in the `.ac9` files is in **normalized form** that is
    the initial `--` is removed from the option name and all the other `-` characters in the string 
    are replaced by an underscore `_`.
    Following these simple rules, the  *configure* option `--with-mpi` becomes `with_mpi` in the ac file.

Save all these options in the *myconf.ac9* file and pass it to the *configure* script using the syntax:

```sh
./configure --with-config-file="myconf.ac9"
```

Remember to crosscheck the summary produced by the configure script to make sure you are 
getting what you expect.

The configure script has generated several files required by make and the *config.h* include file 
containing all the pre-processing options used to build ABINIT.
Let's have a look at this file. 

```cpp
/* Define to 1 if you have a working MPI installation. */
#define HAVE_MPI 1

/* Define to 1 if you have a MPI-1 implementation (obsolete, broken). */
/* #undef HAVE_MPI1 */

/* Define to 1 if you have a MPI-2 implementation. */
#define HAVE_MPI2 1

/* Define to 1 if you want MPI I/O support. */
#define HAVE_MPI_IO 1
```

To get the summary of options activated during the build, run *abinit* with the `-b` option 
(or `--build` if you prefer the verbose version)

```sh
./src/98_main/abinit -b
```

Use `abinit --help` for the full list of options available.

So far so good, we managed

#### Dynamic libraries 

```sh
ldd abinit
```

On MacOsX, replace *ldd* with *otool* and the syntax:

```sh
otool -L abinit
```

To understand why LD_LIBRARY_PATH (or DYLD_LIBRARY_PATH), let's try to reset the value 
of this variable with

```sh
unset LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH
```

now rerun *ldd* (or *otool*) again. 
Do you understand what's happening here? 
Can you execute the abinit executable with an empty $LD_LIBRARY_PATH?
How would you fix the problem?

## How to compile ABINIT on a cluster with the intel toolchain

On intel-based clusters, we suggest to compile ABINIT with the intel compilers and the MKL library 
in order to achieve better performance.
The MKL library, indeed, provides highly-optimized implementations for BLAS, LAPACK, FFT, and SCALAPACK
that can lead to a significant speedup while simplifying considerably the compilation process.

In what follows, we assume a cluster in which the sysadmin has already installed all the modules 
(compilers, MPI and libs) required to compile ABINIT.
If some of the required libraries are lacking, it should not be that difficult to reuse the expertise acquired 
in this tutorial to build and install your own libraries inside $HOME/local although the best solution would be
to ask the syadmin to provide modules with the dependencies required by Abinit or, even better, 
a module that loads the ABINIT executable and all it dependencies (again, feel free to ask your sysadmin to
provide such a module).

For a quick introduction to the environment modules, please consult
[this documentation](https://support.ceci-hpc.be/doc/_contents/UsingSoftwareAndLibraries/UsingPreInstalledSoftware/index.html).

To list the available modules, use:

```sh
module avail
```

Hopefully, the sysadmin has already installed for you all the dependencies required by ABINIT using a consistent tool chain 
i.e. same compiler version, same MPI version. 
If this is not the case, ask the sysadmin to add the missing modules.
If, on the other hand, the sysadmin does not respond to your mails because is busy watching kittens videos on youtube, you can 
always compile the library from source and install it using the `make && make install` procedure discussed in this tutorial.

```sh
module load
```

To list the modules currently loaded use:

```sh
module list
```

At this point, you should decide the toolchain to be used to compile ABINIT:
the compiler e.g. intel ifort, GNU gfortran
This step is strongly cluster-dependent

To see what libraries are recommended for a particular use case, specify the parameters in the drop down lists below.
[mkl-link-line-advisor](https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor)


### How to compile ABINIT with support for OpenMP threads

For a gentle introduction to OpenMP threads, please consult this tutorial.

Compiling ABINIT with OpenMP is not that difficult as everything boils down to:

* Using a **threaded version** for the BLAS/LAPACK/FFT libraries
* Passing the correct options and libraries to the ABINIT configure script

On the contrary, answering the questions:

* When and why should I use OpenMP threads for my calculations?
* How many threads should I use and what is the parallel speedup I should expect?

is much more difficult to answer as there are several factors that should be taken into account.

To keep a long story short, one should use OpenMP threads 
when one starts to trigger limitations or bottlenecks in the MPI implementation, 
especially at the level of the memory requirements or of the parallel scalability.
These problems are usually observed in large calculations that is large [[natom]], [[mpw]], [[nband]])  
Note that it does not make any sense to compile ABINIT with OpenMP if your calculations are relatively small.
ABINIT, indeed, is mainly designed with MPI-parallelism in mind.
Calculations with a relatively large number of $\kk$-points will benefit more of MPI than OpenMP, 
especially if the number of MPI processes divides exactly the number of $\kk$-points.
Even worse, do not compile the code with OpenMP support if you do not plan to use threads because the OpenMP
version will have an additional overhead due to the creation of threaded sections.


!!! Important

    When using threaded libraries remember to set explicitly the number of threads with e.g.

    ```sh
    export OMP_NUM_THREADS=1
    ```

    either in your bash_profile or in the submission script.
    By default, indeed, OpenMP uses all the cores available on the system so it's very easy to overload
    the system especially when one starts to use threads in conjunction with MPI processes.

    Remember also that increasing the number of threads does not necessarily leads to faster calculations.
    There's always an optimal value for the number of threads beyond which the parallel efficiency start to decrease.
    Unfortunately, this value is strongly hardware and software dependent so you will need to benchmark the code
    before running production calculations.

### Trouble shooting 

Remember to use *which*

Inspect *config.log* for possible errors messages, especially the error messages 
associated to the feature you are trying to activate.
Note, indeed, there configure is trying to detect your architecture at runtime.
Tests performed at the begining may fail so not all the
error messages reported in *config.log* are necessarily relevant.

Last but not least, remember to attach the log file when asking for help on the ABINIT forum

How to use gdb


How to use LLDB

$ lldb ../../../src/98_main/abinit
(lldb) target create "../../../src/98_main/abinit"
Current executable set to '../../../src/98_main/abinit' (x86_64).
(lldb) settings set target.input-path t85.in
(lldb) run

Avoid VECLIB on MacOsx

	/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libLAPACK.dylib (compatibility version 1.0.0, current version 1.0.0)
	/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib (compatibility version 1.0.0, current version 1.0.0)

```sh
# Enable ZDOTC and ZDOTU bugfix (gmatteo)
#
enable_zdot_bugfix="yes"
```

### Additional resources


