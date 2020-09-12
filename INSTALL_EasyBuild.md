---
authors: M. Giantomassi
---

# How to install ABINIT with EasyBuild

This page describes how to compile and install ABINIT with [EasyBuild](https://github.com/easybuilders/easybuild),
a python framework for managing scientific software on clusters
that allows one to build an entire software stack and the associated modules with a single command.

For further information about EasyBuild, please consult the [official documentation](https://easybuild.readthedocs.io/en/latest/).
We also recommended to read this [tutorial](https://ulhpc-tutorials.readthedocs.io/en/latest/tools/easybuild/) in which
the different steps required to build and install software are discussed in more detail.
For an introduction to Lmod and Environment Modules,
see the [Lmod user guide](https://lmod.readthedocs.io/en/latest/index.html).

Note also that EasyBuild configuration files (*easyconfigs*) for different ABINIT versions and different toolchains
are available on [github](https://github.com/easybuilders/easybuild-easyconfigs/tree/develop/easybuild/easyconfigs/a/ABINIT).
If you are not interested in learning how to use EasyBuild,
feel free to ask your sysadmin to build an ABINIT module using one of the *easyconfigs* already available.
If, on the other hand, you love learning how to use new tools to facilitate your work
and you like to organize your own software stack in a clean way with automatically-generated modules, continue reading.

## Getting started with EasyBuild

First of all, let's use the `module spider` command to check whether an ABINIT module is already installed:

```sh
$ module spider abinit

------------------------------------------------------------------------------------------------------------------------
  ABINIT:
------------------------------------------------------------------------------------------------------------------------
    Description:
      ABINIT is a package whose main program allows one to find the total energy, charge density and electronic
      structure of systems made of electrons and nuclei (molecules and periodic solids) within Density Functional
      Theory (DFT), using pseudopotentials and a planewave or wavelet basis.

     Versions:
        ABINIT/8.2.2-foss-2016b
        ABINIT/8.4.4-intel-2017b
        ABINIT/8.6.3-intel-2018a
        ABINIT/8.10.2-intel-2017b
        ABINIT/8.10.3-intel-2018b

------------------------------------------------------------------------------------------------------------------------
  For detailed information about a specific "ABINIT" package (including how to load the modules) use the module's full name.
  Note that names that have a trailing (E) are extensions provided by other modules.
  For example:

     $ module spider ABINIT/8.10.3-intel-2018b
------------------------------------------------------------------------------------------------------------------------
```

The output shows that there are five ABINIT versions already installed.
The good news is that we can directly load one of these modules with

```sh
$ module load ABINIT/8.10.3-intel-2018b
```

to activate the 8.10.3 executable compiled with the intel toolchain 2018b.
The bad news is that all the available versions are rather old so
we have to use the EasyBuild command line interface to build a more recent version.

To build software with EasyBuild, we need the `eb` python script.
If `eb` is not already in $PATH, you will get the following error message:

```sh
$ which eb
bash: eb: command not found...
```

In this case, issue `module spider EasyBuild` and follow the instructions printed to the terminal
to load the EasyBuild module, e.g.:

```sh
$ module load EasyBuild/4.2.2

$ which eb
/opt/sw/arch/easybuild/software/EasyBuild/4.2.2/bin/eb
```

Now execute `eb` with the `-S` option to search for ABINIT easyconfigs:

```sh
$ eb -S abinit

== found valid index for /usr/easybuild/easyconfigs, so using it...
== found valid index for /usr/easybuild/easyconfigs, so using it...
CFGS1=/usr/easybuild/easyconfigs
CFGS2=/usr/easybuild/easyconfigs
 * $CFGS1/a/ABINIT/ABINIT-8.0.8-intel-2016a.eb
 * $CFGS1/a/ABINIT/ABINIT-8.0.8b-foss-2016b.eb
 * $CFGS1/a/ABINIT/ABINIT-8.0.8b-intel-2016b.eb
 * $CFGS1/a/ABINIT/ABINIT-8.10.2-intel-2018b.eb
 * $CFGS1/a/ABINIT/ABINIT-8.10.3-intel-2018b.eb
 * $CFGS1/a/ABINIT/ABINIT-8.2.2-foss-2016b.eb
 * $CFGS1/a/ABINIT/ABINIT-8.2.2-intel-2016b.eb
 * $CFGS1/a/ABINIT/ABINIT-8.4.4-intel-2017b.eb
 * $CFGS1/a/ABINIT/ABINIT-8.6.3-intel-2018a.eb
 * $CFGS2/a/ABINIT/ABINIT-8.0.8-intel-2016a.eb
 * $CFGS2/a/ABINIT/ABINIT-8.0.8b-foss-2016b.eb
 * $CFGS2/a/ABINIT/ABINIT-8.0.8b-intel-2016b.eb
 * $CFGS2/a/ABINIT/ABINIT-8.10.2-intel-2018b.eb
 * $CFGS2/a/ABINIT/ABINIT-8.10.3-intel-2018b.eb
 * $CFGS2/a/ABINIT/ABINIT-8.2.2-foss-2016b.eb
 * $CFGS2/a/ABINIT/ABINIT-8.2.2-intel-2016b.eb
 * $CFGS2/a/ABINIT/ABINIT-8.4.4-intel-2017b.eb
 * $CFGS2/a/ABINIT/ABINIT-8.6.3-intel-2018a.eb
 * $CFGS2/a/ABINIT/ABINIT-9.0.4-foss-2019b.eb
 * $CFGS2/a/ABINIT/ABINIT-9.0.4-intel-2019b.eb
```

Various *easyconfig* files are found for different ABINIT versions
and different toolchains: `foss 2016b`, `intel 2018a`, etc.
Note that the output may change depending on the version of EasyBuild installed on your machine as
these `eb` files are shipped with the EasyBuild installation.
More recent *easyconfigs* files can be found in the
[github repository](https://github.com/easybuilders/easybuild-easyconfigs/tree/develop/easybuild/easyconfigs/a/ABINIT).

!!! important

    `eb` files typically follow the naming scheme `<name>-<version>[-<toolchain>][<versionsuffix>].eb`.
    The name reflects the fact each easyconfig file is associated to a particular version of the library/application
    and a particular [toolchain](https://easybuild.readthedocs.io/en/latest/Common-toolchains.html).
    For example, `foss` stands for Free and Open Source Software toolchain based on 
    GCC, OpenMPI, OpenBLAS/LAPACK, ScaLAPACK, and FFTW3
    whereas the `intel` toolchain is based on the intel compilers, MKL and intel MPI.

    Applications based on `foss` are relatively easy to build since all the basic building blocks
    (including the GCC compiler suite) can be compiled from source if needed.
    On the contrary, the `intel` toolchain requires licensed software and valid licenses must be available to install and use it.

    It is worth to mention that applications compiled with the `intel` toolchain are usually more performant
    than the `foss` version especially on intel hardware.
    Unfortunately, compiling and linking with the intel toolchain requires a pre-existent licensed installation
    that is usually done by the sysadmin.
    As a consequence, building software with an intel easyconfig is not that easy if the sysadmin hasn't
    already installed the intel toolchain.

We are lucky as our installation kindly provides two configuration files for ABINIT 9.0.4:

```sh
* $CFGS2/a/ABINIT/ABINIT/9.0.4-foss-2019b.eb
* $CFGS2/a/ABINIT/ABINIT/9.0.4-intel-2019b.eb
```

but before discussing the build process, let's have a look at the `Meta Modules` section that lists the
releases available on our machine:

```sh
$ module avail

------------------ Meta Modules ------------------
   releases/elic-2017b        releases/2016b (S,L)
   releases/2016a      (S)    releases/2017b (S)
   releases/2018a (S)         releases/2019b (S,D)    use.own
   releases/2018b (S)         tis/2018.01    (S,L)

<snip>
```

The `(D)` after the name of the module stands for Default.
It means that EasyBuild uses this release by default.

Since we want to build `9.0.4-foss-2019b.eb` and `2019b` is the default release,
we can proceed directly with the next steps.
If you want to build a recipe associated to a different release e.g. `2016a`,
you need to load the associated module before building the code using e.g.:

```sh
module load releases/2016a`
```

Now invoke `eb` with the name of the easyconfig file and the two options `--dry-run` and `--robot`
(or just `-Dr` if you prefer the short version):

```sh
eb ABINIT-9.1.0-foss-2019b.eb --robot --dry-run
```

`--dry-run` tells `eb` to print all the required dependencies **without actually building software**.
We **strongly** suggest to use this option to understand what's going on before performing a full installation.
The `--robot` option tells `eb` to automatically build and install all dependencies
while searching for easyconfigs files in a set of pre-defined directories.
To prepend additional directories to search for eb files (like the current directory $PWD),
use the syntax `--robot-paths=$PWD`.
Multiple directories can be specified using `--robot-paths=PATH1:PATH2`

On my machine, I get the following results:

```sh
eb ABINIT-9.1.0-foss-2019b.eb --robot --dry-run

== temporary log file in case of crash /tmp/eb-ttdRNl/easybuild-hVtLoS.log
== found valid index for /usr/easybuild/easyconfigs, so using it...
Dry run: printing build status of easyconfigs and dependencies
 * [x] /usr/easybuild/easyconfigs/m/M4/M4-1.4.18.eb (module: M4/1.4.18)
 * [x] /usr/easybuild/easyconfigs/b/Bison/Bison-3.3.2.eb (module: Bison/3.3.2)
 * [x] /usr/easybuild/easyconfigs/z/zlib/zlib-1.2.11.eb (module: zlib/1.2.11)
 * [x] /usr/easybuild/easyconfigs/h/help2man/help2man-1.47.4.eb (module: help2man/1.47.4)
 * [x] /usr/easybuild/easyconfigs/f/flex/flex-2.6.4.eb (module: flex/2.6.4)
 * [x] /usr/easybuild/easyconfigs/b/binutils/binutils-2.32.eb (module: binutils/2.32)
 * [x] /usr/easybuild/easyconfigs/g/GCCcore/GCCcore-8.3.0.eb (module: GCCcore/8.3.0)
 * [x] /usr/easybuild/easyconfigs/h/help2man/help2man-1.47.8-GCCcore-8.3.0.eb (module: help2man/1.47.8-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/m/M4/M4-1.4.18-GCCcore-8.3.0.eb (module: M4/1.4.18-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/z/zlib/zlib-1.2.11-GCCcore-8.3.0.eb (module: zlib/1.2.11-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/b/Bison/Bison-3.3.2-GCCcore-8.3.0.eb (module: Bison/3.3.2-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/f/flex/flex-2.6.4-GCCcore-8.3.0.eb (module: flex/2.6.4-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/b/binutils/binutils-2.32-GCCcore-8.3.0.eb (module: binutils/2.32-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/g/GCC/GCC-8.3.0.eb (module: GCC/8.3.0)
 * [x] /usr/easybuild/easyconfigs/s/Szip/Szip-2.1.1-GCCcore-8.3.0.eb (module: Szip/2.1.1-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/n/ncurses/ncurses-6.1-GCCcore-8.3.0.eb (module: ncurses/6.1-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/e/expat/expat-2.2.7-GCCcore-8.3.0.eb (module: expat/2.2.7-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/b/bzip2/bzip2-1.0.8-GCCcore-8.3.0.eb (module: bzip2/1.0.8-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/l/libreadline/libreadline-8.0-GCCcore-8.3.0.eb (module: libreadline/8.0-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/c/cURL/cURL-7.66.0-GCCcore-8.3.0.eb (module: cURL/7.66.0-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/p/pkg-config/pkg-config-0.29.2-GCCcore-8.3.0.eb (module: pkg-config/0.29.2-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/l/libtool/libtool-2.4.6-GCCcore-8.3.0.eb (module: libtool/2.4.6-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/o/OpenBLAS/OpenBLAS-0.3.7-GCC-8.3.0.eb (module: OpenBLAS/0.3.7-GCC-8.3.0)
 * [x] /usr/easybuild/easyconfigs/c/CMake/CMake-3.15.3-GCCcore-8.3.0.eb (module: CMake/3.15.3-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/p/Perl/Perl-5.30.0-GCCcore-8.3.0.eb (module: Perl/5.30.0-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/d/Doxygen/Doxygen-1.8.16-GCCcore-8.3.0.eb (module: Doxygen/1.8.16-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/a/Autoconf/Autoconf-2.69-GCCcore-8.3.0.eb (module: Autoconf/2.69-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/a/Automake/Automake-1.16.1-GCCcore-8.3.0.eb (module: Automake/1.16.1-GCCcore-8.3.0)
 * [ ] /usr/easybuild/easyconfigs/l/libxc/libxc-4.3.4-GCC-8.3.0.eb (module: libxc/4.3.4-GCC-8.3.0)
 * [x] /usr/easybuild/easyconfigs/a/Autotools/Autotools-20180311-GCCcore-8.3.0.eb (module: Autotools/20180311-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/n/numactl/numactl-2.0.12-GCCcore-8.3.0.eb (module: numactl/2.0.12-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/x/xorg-macros/xorg-macros-1.19.2-GCCcore-8.3.0.eb (module: xorg-macros/1.19.2-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/l/libpciaccess/libpciaccess-0.14-GCCcore-8.3.0.eb (module: libpciaccess/0.14-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/n/ncurses/ncurses-6.0.eb (module: ncurses/6.0)
 * [x] /usr/easybuild/easyconfigs/g/gettext/gettext-0.19.8.1.eb (module: gettext/0.19.8.1)
 * [x] /usr/easybuild/easyconfigs/x/XZ/XZ-5.2.4-GCCcore-8.3.0.eb (module: XZ/5.2.4-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/l/libxml2/libxml2-2.9.9-GCCcore-8.3.0.eb (module: libxml2/2.9.9-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/h/hwloc/hwloc-1.11.12-GCCcore-8.3.0.eb (module: hwloc/1.11.12-GCCcore-8.3.0)
 * [x] /usr/easybuild/easyconfigs/o/OpenMPI/OpenMPI-3.1.4-GCC-8.3.0.eb (module: OpenMPI/3.1.4-GCC-8.3.0)
 * [x] /usr/easybuild/easyconfigs/g/gompi/gompi-2019b.eb (module: gompi/2019b)
 * [x] /usr/easybuild/easyconfigs/f/FFTW/FFTW-3.3.8-gompi-2019b.eb (module: FFTW/3.3.8-gompi-2019b)
 * [x] /usr/easybuild/easyconfigs/s/ScaLAPACK/ScaLAPACK-2.0.2-gompi-2019b.eb (module: ScaLAPACK/2.0.2-gompi-2019b)
 * [x] /usr/easybuild/easyconfigs/h/HDF5/HDF5-1.10.5-gompi-2019b.eb (module: HDF5/1.10.5-gompi-2019b)
 * [x] /usr/easybuild/easyconfigs/n/netCDF/netCDF-4.7.1-gompi-2019b.eb (module: netCDF/4.7.1-gompi-2019b)
 * [x] /usr/easybuild/easyconfigs/n/netCDF-Fortran/netCDF-Fortran-4.5.2-gompi-2019b.eb (module: netCDF-Fortran/4.5.2-gompi-2019b)
 * [x] /usr/easybuild/easyconfigs/f/foss/foss-2019b.eb (module: foss/2019b)
 * [ ] /home/users/g/m/gmatteo/try_eb/ABINIT-9.1.0-foss-2019b.eb (module: ABINIT/9.0.4-foss-2019b)
== Temporary log file(s) /tmp/eb-ttdRNl/easybuild-hVtLoS.log* have been removed.
== Temporary directory /tmp/eb-ttdRNl has been removed.
```

The output indicates that most of the dependencies are already installed (**checked boxes**).
Only *libxc/4.3.4-GCC-8.3.0* and *ABINIT-9.1.0-foss-2019b.eb* must be built from source (**unchecked boxes**).

Now we can finally build and install ABINIT v9.1.0 with *foss-2019b*
by executing the same command without `--dry-run`:

```sh
eb ABINIT-9.1.0-foss-2019b.eb --robot

== temporary log file in case of crash /tmp/eb-4IOAjW/easybuild-vBwXik.log
== found valid index for /usr/easybuild/easyconfigs, so using it...
== resolving dependencies ...
== processing EasyBuild easyconfig /usr/easybuild/easyconfigs/l/libxc/libxc-4.3.4-GCC-8.3.0.eb
== building and installing libxc/4.3.4-GCC-8.3.0...
== fetching files...
== creating build dir, resetting environment...
== starting iteration #0 ...
== unpacking...
== patching...
== preparing...
== configuring...
== building...
== testing...
== installing...
== taking care of extensions...
== creating build dir, resetting environment...
== starting iteration #1 ...
== unpacking...
== patching...
== preparing...
== configuring...
== building...
== testing...
== installing...
== taking care of extensions...
== restore after iterating...
== postprocessing...
== sanity checking...
== cleaning up...
== creating module...
== permissions...
== packaging...
== COMPLETED: Installation ended successfully (took 18 min 29 sec)
== Results of the build can be found in the log file(s) /home/ucl/modl/gmatteo/.local/easybuild/software/libxc/4.3.4-GCC-8.3.0/easybuild/easybuild-libxc-4.3.4-20200904.010629.log
== processing EasyBuild easyconfig /home/users/g/m/gmatteo/try_eb/ABINIT-9.1.0-foss-2019b.eb
== building and installing ABINIT/9.0.4-foss-2019b...
== fetching files...
== creating build dir, resetting environment...
== unpacking...
== patching...
== preparing...
== configuring...
== building...
== testing...
== installing...
== taking care of extensions...
== restore after iterating...
== postprocessing...
== sanity checking...
== cleaning up...
== creating module...
== permissions...
== packaging...
== COMPLETED: Installation ended successfully (took 7 min 4 sec)
== Results of the build can be found in the log file(s) /home/ucl/modl/gmatteo/.local/easybuild/software/ABINIT/9.0.4-foss-2019b/easybuild/easybuild-ABINIT-9.0.4-20200904.011333.log
== Build succeeded for 2 out of 2
== Temporary log file(s) /tmp/eb-4IOAjW/easybuild-vBwXik.log* have been removed.
== Temporary directory /tmp/eb-4IOAjW has been removed.
```

As you can see from the log messages, *eb* has compiled *libxc* and finally *abinit*.
Executables, libraries and modules are now installed in our `~/.local/easybuild/` directory.

```sh
$ ls ~/.local/easybuild/
build  ebfiles_repo  modules  software  sources

$ tree  ~/.local/easybuild/modules/
/home/ucl/modl/gmatteo/.local/easybuild/modules/
|-- all
|   |-- ABINIT
|   |   `-- 9.0.4-foss-2019b.lua
|   `-- libxc
|       `-- 4.3.4-GCC-8.3.0.lua
`-- chem
    |-- ABINIT
    |   `-- 9.0.4-foss-2019b.lua -> /home/ucl/modl/gmatteo/.local/easybuild/modules/all/ABINIT/9.0.4-foss-2019b.lua
    `-- libxc
        `-- 4.3.4-GCC-8.3.0.lua -> /home/ucl/modl/gmatteo/.local/easybuild/modules/all/libxc/4.3.4-GCC-8.3.0.lua
```

Since we need to load modules installed in `~/.local/easybuild/modules/all`,
we must add this directory to $MODULEPATH.
With a recent version of EasyBuild, issue:

```sh
module use use.own
```

If the `use.own` module is not available, execute the command:

```sh
module use ~/.local/easybuild/modules/all
```

At this point, we should have our version of ABINIT and the dependencies available as `User Modules`

```sh
$ module avail

------------- User Modules -------------
   ABINIT/9.0.4-foss-2019b    libxc/4.3.4-GCC-8.3.0
```

We can finally load our brand-new ABINIT module with:

```sh
module load ABINIT/9.0.4-foss-2019b
```

and check that the abinit executable is in $PATH

```sh
$ which abinit
~/.local/easybuild/software/ABINIT/9.0.4-foss-2019b/bin/abinit
```

and that we have the correct version:

```sh
$ abinit -v
9.0.4
```

The output of `ldd` shows that the executable is dynamically linked to our version of libxc
whereas *openblas*, *fftw3*, *netcdf*, *hdf5* and MPI are provided by the `2019b` release of `foss` as expected:

```sh
$ ldd `which abinit`

libxc.so.5 => /home/ucl/modl/gmatteo/.local/easybuild/software/libxc/4.3.4-GCC-8.3.0/lib64/libxc.so.5 (0x00002aaaaaccf000)
libfftw3f.so.3 => /opt/sw/arch/easybuild/2019b/software/FFTW/3.3.8-gompi-2019b/lib/libfftw3f.so.3 (0x00002aaaab2d9000)
libfftw3.so.3 => /opt/sw/arch/easybuild/2019b/software/FFTW/3.3.8-gompi-2019b/lib/libfftw3.so.3 (0x00002aaaab5c0000)
libopenblas.so.0 => /opt/sw/arch/easybuild/2019b/software/OpenBLAS/0.3.7-GCC-8.3.0/lib/libopenblas.so.0 (0x00002aaaab856000)
libnetcdff.so.7 => /opt/sw/arch/easybuild/2019b/software/netCDF-Fortran/4.5.2-gompi-2019b/lib/libnetcdff.so.7 (0x00002aaaaaad3000)
libnetcdf.so.15 => /opt/sw/arch/easybuild/2019b/software/netCDF/4.7.1-gompi-2019b/lib64/libnetcdf.so.15 (0x00002aaaaab59000)
libhdf5_hl.so.100 => /opt/sw/arch/easybuild/2019b/software/HDF5/1.10.5-gompi-2019b/lib/libhdf5_hl.so.100 (0x00002aaaaac95000)
libhdf5.so.103 => /opt/sw/arch/easybuild/2019b/software/HDF5/1.10.5-gompi-2019b/lib/libhdf5.so.103 (0x00002aaaac607000)
libmpi.so.40 => /opt/sw/arch/easybuild/2019b/software/OpenMPI/3.1.4-GCC-8.3.0/lib/libmpi.so.40 (0x00002aaaad432000)
```

The binary works out of the box because the ABINIT module is automatically loading all the dependencies we need
and our local installation is prepended to `LD_LIBRARY_PATH` and `PATH`.
To see this, use `module show`:

```sh
module show ABINIT/9.0.4-foss-2019b

whatis("Description: ABINIT is a package whose main program allows one to find the total energy,
 charge density and  electronic structure of systems made of electrons and nuclei (molecules
 and periodic solids) within Density Functional  Theory (DFT), using pseudopotentials and a
 planewave or wavelet basis.")
whatis("Homepage: https://www.abinit.org/")
whatis("URL: https://www.abinit.org/")
conflict("ABINIT")
load("foss/2019b")
load("libxc/4.3.4-GCC-8.3.0")
load("HDF5/1.10.5-gompi-2019b")
load("netCDF/4.7.1-gompi-2019b")
load("netCDF-Fortran/4.5.2-gompi-2019b")
prepend_path("CMAKE_PREFIX_PATH","/home/ucl/modl/gmatteo/.local/easybuild/software/ABINIT/9.0.4-foss-2019b")
prepend_path("LD_LIBRARY_PATH","/home/ucl/modl/gmatteo/.local/easybuild/software/ABINIT/9.0.4-foss-2019b/lib")
prepend_path("LIBRARY_PATH","/home/ucl/modl/gmatteo/.local/easybuild/software/ABINIT/9.0.4-foss-2019b/lib")
prepend_path("PATH","/home/ucl/modl/gmatteo/.local/easybuild/software/ABINIT/9.0.4-foss-2019b/bin")
prepend_path("PKG_CONFIG_PATH","/home/ucl/modl/gmatteo/.local/easybuild/software/ABINIT/9.0.4-foss-2019b/lib/pkgconfig")
setenv("EBROOTABINIT","/home/ucl/modl/gmatteo/.local/easybuild/software/ABINIT/9.0.4-foss-2019b")
setenv("EBVERSIONABINIT","9.0.4")
setenv("EBDEVELABINIT","/home/ucl/modl/gmatteo/.local/easybuild/software/ABINIT/9.0.4-foss-2019b/easybuild/ABINIT-9.0.4-foss
-2019b-easybuild-devel")
```

Q: What happens if you run

```sh
module load releases/2016a
eb ABINIT-9.1.0-foss-2019b.eb -Dr
```

## Using a customized easyconfig

In the previous example we used one of the official easyconfig files shipped with the EasyBuild package.
There are cases, however, in which we need to perform a customized build.
For example, we may want to compile the latest stable version of ABINIT recently released,
or use a different release, or maybe we just need to activate configuration options that are not present in the official file.
In this case, we can simply copy the original easyconfig file,
change it according to our needs and finally pass this customized file to the `eb` script.

First of all, let's have a look at `ABINIT-9.0.4-foss-2019b.eb`:

```python
easyblock = 'ConfigureMake'

name = 'ABINIT'
version = '9.0.4'

homepage = 'https://www.abinit.org/'
description = """ABINIT is a package whose main program allows one to find the total energy,
 charge density and  electronic structure of systems made of electrons and nuclei (molecules
 and periodic solids) within Density Functional  Theory (DFT), using pseudopotentials and a
 planewave or wavelet basis."""

toolchain = {'name': 'foss', 'version': '2019b'}
toolchainopts = {'usempi': True, 'pic': True}

source_urls = ['https://www.abinit.org/sites/default/files/packages/']
sources = [SOURCELOWER_TAR_GZ]
#checksums = ['82e8d071088ab8dc1b3a24380e30b68c544685678314df1213180b449c84ca65']

dependencies = [
    ('libxc', '4.3.4'),
    ('netCDF-Fortran', '4.5.2'),
]

# Ensure MPI.
configopts = '--with-mpi="yes" --enable-openmp="no" '

# BLAS/Lapack
configopts += '--with-linalg-flavor="openblas"  LINALG_LIBS="-L${EBROOTOPENBLAS}/lib -lopenblas" '

# FFTW3 support
configopts += '--with-fft-flavor=fftw3 FFTW3_LIBS="-L${EBROOTFFTW} -lfftw3f -lfftw3" '

# libxc support
configopts += '--with-libxc=${EBROOTLIBXC} '

# hdf5/netcdf4.
configopts += 'with_netcdf="${EBROOTNETCDF}" '
configopts += 'with_netcdf_fortran="${EBROOTNETCDFMINFORTRAN}" '
configopts += 'with_hdf5="${EBROOTHDF5}" '

# make sure --free-line-length-none is added to FCFLAGS
configopts += 'FCFLAGS="${FCFLAGS} --free-line-length-none" '

runtest = 'check'

sanity_check_paths = {
    'files': ['bin/%s' % x for x in ['abinit', 'aim', 'cut3d', 'conducti', 'mrgddb', 'mrgscr', 'optic']],
    'dirs': ['lib/pkgconfig'],
}

moduleclass = 'chem'
```

The meaning of the different easyconfig parameters
is explained [here](https://easybuild.readthedocs.io/en/latest/Writing_easyconfig_files.html).

In most of the cases, a customized build requires changing one of the following entries:

- `version`: string with the version number. Used to download the tarball from `source_urls`
- toolchain version
- dependencies
- `configopts`: string with the options passed to the ABINIT configure script

Changing `version` and/or toolchain version is straightforward.
To use a different version of ABINIT (e.g. 9.2.0), simply change the value of `version`.

```diff
- version = '9.0.4'
+ version = '9.2.0'
```

save the new file as `ABINIT-9.2.0-foss-2019b.eb` and execute:

```
eb ./ABINIT-9.2.0-foss-2019b.eb --robot
```

to install a new module for version 9.2.0.
Everything should work out of the box provided version 9.2.0 is still supporting
the version of the libraries listed in `dependencies`.
By the same token, one can compile with version `2020a` of the `foss` toolchain by just changing

```diff
- toolchain = {'name': 'foss', 'version': '2019b'}
+ toolchain = {'name': 'foss', 'version': '2020a'}
```

provided `releases/2020a` is available on our cluster.

Changing `dependencies` and `configopts` requires some basic understanding of the ABINIT build system
and of the logic used by EasyBuild to automate the compilation.
If you need a specialized `eb` file, feel free to contact the ABINIT developers on the forum or the EasyBuild 
developers to ask for support.

<!--
Unfortunately, we cannot cover all the possible scenarios so we focus on the simplest case in which
we want to **extend** the recipe that is we just want to activate support for an optional library
without changing the other easyconfigs parameters.
Let's assume, for instance, that we want to build an Abinit version with Wannier90 support.

Expert users may want to activate support for optional features that require additional external libraries
e.g. wannier90, atompaw, etc.
In this case, one has to list the external library in the `dependencies` section, add the correct configuration
options to `configopts` and make sure that `eb` can find an easyconfig file that can be used to build the library.
-->
