---
authors: MG
---

# Configuration and compilation

This page collects FAQs related to

## I'm new to configure. Can you explain how to use it?

Please consult <https://wiki.abinit.org/doku.php?id=build:configure>

To obtain the full list of available options and their documentation, use

```
./configure --help
```

As an example, the option to specify the Fortran compiler vendor appears as with_fc_vendor in the config file,
while it is typed --with-fc-vendor on the command line.

You will find detailed instructions on how to set the various options of the configure script
in the template file ~abinit/doc/build/config-template.ac, including how to have your config file
included automatically each time you build Abinit.

There are only 2 undocumented options in this template:

the --with-config-file option, which lets you use arbitrary names for your config files,
as long as you specify it each time on the command line;
the --disable-config-file option, to force Abinit not to use any config file.
For obvious reasons, they are the only options that cannot be stored in a config file.

## Configure complains that Abinit cannot be built!

Well, there are many reasons why configure can fail.

First of all, make sure that your LD_LIBRARY_PATH (DYLD_LIBRARY_PATH on MacOs), PATH environment variables
are properly set **before** running configure.
By the same token, one should load all the modules before configuring.
If your environment is propertly set and configure keeps on failing, you will need to search 
for inside *config.log*.

## Compilation stops with "internal compiler error"

This usually indicates a bug in the compiler that should be reported to the compiler vendor.
To bypass the problem, reduce the optimization level from e.g. -O2 (default) to -O1 or even -O0.
If you don't like the idea of compiling all the source file with suboptimal options, you may try
to disable optimation only for the problematic source files using this quick and dirty recipe:

- Run *make* till you trigger the internal compiler error
- Copy the compilation command that triggers the internal error
- Add `-O0` to the compiler options and execute this command line to compile the module.

## Cannot run Abinit due to "undefined reference" error

This means that your Abinit executable is **dynamically linked** with external libraries
but the OS cannot find these libraries at runtime.
Please make sure that LD_LIBRARY_PATH (DYLD_LIBRARY_PATH on MacOs) is properly set that is you have
the same value as the one used when configuring/compiling ABINIT.
The same error message may show up if you forgot to load modules before trying to run the executable.
Again, you are supposed to **load the same modules** that were used for configuring/compiling Abinit.

## Is there any precompiled package for Abinit?

Yes, please consult the different subsections available in the [installation](/installation) page.

## Where can I find examples of configuration files for HPC clusters?

The |abiconfig| package provides configuration files used to configure and compile Abinit on HPC clusters.
Each configuration file contains a header with metadata in json format followed by the configure options
supported  by the Abinit build system.
Examples for Abinit v9 are available [here](https://github.com/abinit/abiconfig/tree/master/abiconfig/clusters).

## Is there any EasyBuild recipe for Abinit?

Yes. An HowTo tutorial is also available [here](/INSTALL_EasyBuild)

## Is there any Spack recipe for Abinit?

Yes but, at present, only Abinit8 is supported.

## Is there any conda package for Abinit?

Yes. conda install abinit --channel conda-forge

Note, however, we do not reccomend using precompiled version on HPC clusters.
The version provided by Spack/EasyBuild pa

## Why do we provide fallbacks?

Abinit Fallbacks is a package builder for the external dependencies of ABINIT in environments lacking these components.
Note, indeed, that **all** the external libraries should be compiled with the same MPI wrapper/compiler
used to build ABINIT.
Unfortunately, not all HPC centers provide a complete set modules covering the external dependencies required by ABINIT.
<!--
In the worst case scenario, the modules are available but they do not work properly
They do not provide full support for the advanced features of Abinit nor HPC-grade calculation capabilities.
They are designed to let developers quickly test new versions of these external dependencies in various situations
before proposing upgrades, as well as to allow end-users to run calculations on their favorite PCs.
In case some dependencies are missing on your computers,
Abinit provides fallback libraries that you can build and install from their sources before compiling Abinit itself.
-->

!!! important

    Feel free to contact the sysadmin of your cluster to ask him/her to install these libraries.

## How can I compile Abinit with OpenMP support?

Compiling ABINIT with OpenMP support is not that difficult once we understand that OpenMP can be activated
at two different levels.

1. The Abinit Fortran routines
2. External libraries for BLAS/Lapack and FFTs.

To activate OpenMP in 1), add to the .ac file

```
enable_openmp="yes"
FCFLAGS_OPENMP="-fopenmp"    # for the GFORTRAN compiler
#FCFLAGS_OPENMP="-qopenmp"   # for the INTEL FORTRAN compiler
```

Finally, remember to set the environment variable

```
export OMP_NUM_THREADS=1
```

!!! important

    Do not mix MKL with FFTW3, especially when activating OpenMP.
