# Installation notes for ABINIT

This file provides a description of the operations needed to install the
ABINIT package, to generate the executable and to make the tests. 
It provides also the description of simple modifications of the package, for developers.

Examples of configuration files to configure and compile Abinit on clusters are available
in the |abiconfig| package on github while the configuration files 
used for our buildbot testfarm are available in the [autoconf_examples section](developers/autoconf_examples/).

<!--
See a recent version of the [new user's guide](..),
for an introduction to the abinit package. 
See a recent version of the [[help:abinit]] file for learning how to use the code. 
Both of them can be found either on the Web, or in the doc subdirectory of the package.
-->

Any comment or suggestion to improve the procedure will be welcome! 
Simply contact the ABINIT group <https://forum.abinit.org/>

## Short Howto

For the vast majority of people willing to use ABINIT (Unix/Linux, not developers, but users), 
here follows a short list of instructions needed to install it:

  1. Download, gunzip and untar the 
    [latest version of the ABINIT package](https://www.abinit.org/packages) (abinit-x.y.z.tar.gz) 
    from the abinit Web site, then change the current directory to the top of the directory that was created. 
    If you do not know what all this means, go to [How to get a version of ABINIT? (Section 1B Normal users)](#1b-normal-user)
  2. Issue configure or ./configure (or first create a tmp directory, then cd tmp, then ../configure)
  3. Issue make (or make mj4 , or make multi multi_nprocs=n for using "n" processors on a SMP machine 
    where you have to replace "n" by its value) 
  4. Issue (optionally) "make install"

Well, it might also be that information on the Fortran compiler is needed, in
which case something like:

    ./configure FC=your_F90_compiler 

instead of the bare "configure", might work, where "your_F90_compiler" has to
be replaced by the location of your F90 compiler, such as /usr/local/gcc472/bin/gfortran.

If you succeeded to download, gunzip and untar the ABINIT package, but failed
with the next steps, please go to **2\. How to make the executables?**

If you succeeded to make the executables, but would like to check whether
ABINIT has been installed correctly, please go to **3\. How to make the
internal tests?** and the following sections.

If you want to have a much better handling on ABINIT than normal users, or if
you download ABINIT from Gitlab or GitHub, then go to **1\. How to get a
version of ABINIT? (Section 1A Expert users of developer)**.

To build the executables you will need at least 400 MBytes of free disk space.
Running all tests will require more than 2.5 GBytes.

## How to get a version of ABINIT?

We will distinguish two cases:

  * 1.A. Expert user or developer. 
    You have a F90 compiler under UNIX/Linux or MacOS X, as well as (free) **software applications 
    like git, automake, autoconf, libtool, perl, python, and you want to have a full handle** on the package 
    (compilation, modification of files, writing of scripts ... this is the preferred case for developers). 
    This is also needed if you download ABINIT from the primary ABINIT GitLab repository or its GitHub mirror. 
    We will sometimes refer to this case as the "autotools working mode".

  * 1.B. Normal user. You have a F90 compiler under UNIX/Linux or MacOS X and you want 
    **simply to compile the source code**, and, from time to time, **modify and/or add a new file** 
    (this is the case of most users, system managers, and also occasional developers). 
    We will sometimes refer to this case as being "without autotools".

You should read only the appropriate section (you can safely ignore the other one ...).

### 1.A. Expert user or developer

You have a F90 compiler under UNIX/Linux or MacOS X, as well as (free)
software like git, automake, autoconf, libtool, perl, python, and you want to
have a full handle on the package (compilation, modification of files, writing
of scripts). This is the preferred case for developers... We will sometimes
refer to this case as the "autotools development mode".

If you do not have these tools, and would like to have them, please consult
your local computer guru, and the following pages:

  * [An overview of ABINIT development](https://wiki.abinit.org/doku.php?id=developers:overview)
  * [10 steps to hike ABINIT](https://wiki.abinit.org/doku.php?id=developers:hike)
  * [Buildbot and the test farm](https://wiki.abinit.org/doku.php?id=bb:overview)

If you want to develop on a regular basis, please have a Git(lab) access
created for you by contacting Jean-Michel Beuken, as described in these pages.
If it is only very occasional, you might as well rely on the [ABINIT Github Web site](https://github.com/abinit).

It is strongly advised to subscribe to the [ABINIT forum](https://forum.abinit.org/) 
if you want to be able to get the latest information concerning the autotools development mode.

After having installed git, and obtained a branch on the ABINIT worldwide
repository, create an automomous copy of the source code, on top of which you
have to make your development. 
This is explained in the ABINIT wiki 
[git(lab): ABINIT specificities](https://wiki.abinit.org/doku.php?id=developers:specificities_git_abinit)

For your branches on the ABINIT worldwide repository, you will have the
permission not only to clone/fetch/pull, but also to commit/push your
modifications. You might alternatively download other branches of the
archives, but you will not be able to commit in these branches. So, do not
start to modify these, you will not be able to include them afterwards in the archive.

Working with `git clone` creates a local archive for your daily work, this
archive being linked to the main ABINIT archive. This very efficient technique
is recommended, as it makes you more independent for the management of your
work (you will be able to create new branches). One big advantage of this
technique is that people working with a laptop can develop and commit safely
without a network connection.

Now, cd to the newly created abinit directory, and issue:

    ./config/scripts/makemake

This command initializes a whole set of files and scripts, needed for the
autotools, as well as for the global work on ABINIT sources. This
initialization might take up to two minutes.
After this initialisation, you can proceed to the generation of the
executables, as described in section 2.

### 1.B. Normal user

You have a F90 compiler under UNIX/Linux or MacOS X and you want simply to
**compile the source files**, and, from time to time, **modify and/or add a
new file**. This is the case of most users, system managers, and also many
developers. If you want to modify and/or add a new file, please consult the
[section "For developers: How to modify the code?"](#for-developers-how-to-modify-the-code).
after reading the present
section. In what follows, _x.y.z_ represents the ABINIT version.

In order to get the ABINIT package, you have first to download the file
**abinit-_x.y.z_.tar.gz** from the ABINIT Web site, see 
[the following page](http://www.abinit.org/packages) then issue:

    gunzip abinit-_x.y.z_.tar.gz | tar -xvf -

That's it.

The abinit-_x.y.z_.tar.gz gzipped tar file contains all the needed files, including:

  * the sources of the abinit code (also, the files needed for generating the blas and lapack libraries), 
    in the directories "src" and "lib";
  * the documentation, in the directory "doc";
  * the complete set of tests, and the pseudopotentials needed for the tests, in the directory "tests";
  * all the scripts and information needed to produce makefiles, in other directories, especially "config".

The package does not contain the object files and the binary executable files.

You can go to the next section, to generate the executable files, if this worked.
If this did not work, here are more detailed explanations ...

So, execute the following actions:

  * under mac os x, open a terminal session, so you can work as if it were a unix platform.
  * transfer the above-mentioned file(s) to your machine, in a directory suitable for 
    the installation of the present version of abinit, and subsequent ones.  
    you should have about 250 mb of disk space to install the code, maybe more,
    depending on the version, and the number of tests that you will do.

  * gunzip (on some machine you need gzip -d) and untar the file abinit-_x.y.z_.tar.gz:  
    gunzip abinit-_x.y.z_.tar.gz | tar -xvf -

If correctly done, a main directory, denoted ~abinit in the present document
(usually, its real name will be **abinit-_x.y.z_**) and a whole set of
subdirectories should have been created. One of them is called 'doc/'. It
contains many important informations. many of the files it contains are .html
files, that are placed on the web site. However, many other files are not
available in .html format, and are not found on the web site. In the future,
keep in mind that the information that you are looking for (but that you
cannot find on the web site) might be in the subdirectories of doc/, esp.
doc/theory/, doc/users/, doc/psp_infos/.

You might now go to the section 2.

Do not forget: if you want to modify and/or add a new file, please consult the
section 7\. for developers: how to modify the code? after reading the present section.

## 2\. How to make the executables? (=How to compile the executables?)

We now suppose that you have a F90 compiler and you want to compile the source files.

In most cases, you will have to provide to the 'make' utility some
information: the location of the F90 compiler (and sometimes even the C
compiler) on your machine, the adequate compiler options, and, if you want to
produce the parallel binaries, the location of the MPI library on your machine.

Although the presently implemented building tools should be powerful enough to
succeed to make the binaries without you giving such information, it has been
seen that on a significant number of platforms, it is still better to give
them. Indeed, you might generate a clearly suboptimal executable, executing
slowly, or with downgraded capabilities.

Supposing that you are in the lucky case where the build system is able to
find all the information, then the build of ABINIT is very simple. Issue:

  * configure or ./configure (or first create a tmp directory, then cd tmp, then ../configure)
  * make (or make mj4 , or make multi multi_nprocs=n for using "n" processors on a SMP machine 
    where you have to replace "n" by its value) 
  * (optionally) make install

Well, it might also be that only one additional information is needed, in which case something like:

    configure FC=gcc
    make 

might work. In both cases, let's explain a bit what is done, and the further possibilities.

The 'configure' step produces the set of Makefile files (among other things),
taking into account information about your machine and the hostname.ac file.
It takes three minute long, or less. The 'make' step compiles everything,
according to the Makefile files produced in the prior step. The time to make
everything is highly dependent on the compiler and platform. On a 2.8 GHz
quad-proc machine (using make mj4), the whole compilation is about 5 minutes.
On some other platforms, with only one processor, it might be more than one hour.

The executables will be located in the subdirectory ~abinit/src/98_main, if
you have chosen to issue ./configure in the ~abinit directory. If you have
issued ./configure in another directory, it will be placed accordingly.

The 'make' command can also be used in many different ways, by mentioning one
or more targets. A (partial) list of targets for users can be obtained by typing:
   
    make help

Additional targets, for developers, can be obtained by typing:
    
    make help_dev

It is possible to compile only one of the executable, just after the configure
step by typing:
  
    make name_of_the_binary

(where name_of_the_binary can be abinit, cut3d, anaddb, ...).

These are only a tiny fraction of the things that you can realize thanks to
'make'. Moreover, there are also 'Makefile' files in most of the
subdirectories of ABINIT, with often their own (partial) list of targets for
users (and also sometimes for developers). To obtain these lists, go to the
directory, and type:
  
    make help

or
    
    make help_dev

Finally:
    
    make install

will install abinit in the /usr/local directory.

In case you want to go further, please consult files in ~abinit/doc/build .

Let's come back to the case where the build system needs some more
information. This information should be stored in a file named hostname.ac,
where "hostname" is the result of executing the command `hostname` on your
machine, e.g. abiref.pcpm.ucl.ac.be or my_machine ... , and taking the first
word of the returned chain of character, e.g. abiref or my_machine ...

There is a template for such a file, located in ~abinit/doc/config/. Its name
is config-template.ac8. Examples of such files, that have been used for testing
the package, can be found in ~abinit/doc/build/config-examples/. By the way,
the known problems observed for these different tests are mentioned in the
~abinit/KNOWN_PROBLEMS file, and the hostname.ac files are correspondingly
indicated at the beginning of this file.

Most of the examples provided in the ~abinit/doc/build/config-examples/
directory contain about five definitions: F90 and C locations, F90 and C
options, MPI library location (or the indication that MPI is not enabled). On
the other hand, there are many other possible control flags ("with_XYZ"),
needed for advanced use. In case you have trouble with some library (LibXC,
WANNIER90, ETSF_IO ...), you may disable its build.

Your hostname.ac file might be placed in your home directory in a new
directory that you will name ~/.abinit/build/. At that location, everytime you
install a new version of ABINIT, the needed information will be found by
ABINIT, so you do not have to care anymore about this file after the first installation.

On the other hand, if you need to play with several computers, you can place
the hostname.ac file directory in the ~abinit directory, where such a
hostname.ac file will be also seen by the build system (and preferred over the
one located in ~/.abinit/build/) or in your build directory (like
~abinit/tmp/). As mentioned above, you might even type at the terminal the
definitions contained in the hostname.ac file.

Note the order of precedence for the location of the hostname.ac file (or
command-line information), in case more than one possibility is used,
(decreasing order of precedence):

  * Command line (overcome all other information)
  * Your build directory (~abinit/tmp/)
  * The ABINIT top source directory (~abinit/)
  * ~/.abinit/build/
  * /etc/abinit/build/

When the hostname.ac file is ready, you have to issue, in the ~abinit directory:

  * configure or ./configure (or first create a tmp directory, then cd tmp, create a hostname.ac file, then ../configure)
  * make or make mj4 (or make multi for using several processors on a SMP machine)
  * (optionally) make install

## How to make the internal tests?

In case you are running under Unix (Linux or another flavour), the abinit code
has several small internal tests (three basic ones, called "fast", "v1" and
"v5", and then one for each of the libraries "bigdft", "etsf_io", "libxc",
"wannier90"), that can be issued automatically, and that check themselves
whether the results that have been obtained are right or wrong. These tests
are available whether you have got the package from the Web or from the ABINIT
archive. Of course, you need to have compiled abinit in order to run the
internal tests. Moreover, the simple implementation procedure assumes that the
executable is located in ~abinit/src/98_main (the standard location after
issuing "make").

You can begin with the test "fast". Simply issue the command:
    
    make test_fast

It will run during a few seconds. It should print:

```text
Status file, reporting on built-in test fast

==> The run finished cleanly.
    Moreover, comparison of the total energy, and other (few) relevant quantities with reference values has been successful.
    This does not mean that no problem is present, however.
    Please run the complete set of ABINIT tests to gain a better confidence in your installation.
```

This means that the internal test "fast" ran successfully. If you do not get
this message, then the executables were not properly generated, or there is a
problem with the makefile that drives the internal test. In this case, after
having tried to solve the problem by yourself, you should contact somebody in
the ABINIT group.

In addition to this small message, you can have access to all generated files,
that are located inside the tests/built-in/Input subdirectory.

Supposing test "fast" was OK, then you might issue the command:
    
    make tests_in

The test "fast" will be done once more, followed by the other internal tests.
Again, we hope that you will get the positive diagnostics for the other tests.
Of course, the "bigdft", "etsf_io", "libxc", and "wannier90" needs the
appropriate library to be installed in order to work properly.

For further information on these internal tests, see the ~abinit/tests/built-in/README file.

You might now read the [new user's guide](..), in
order to learn how to use the code, and then to follow the four basic
tutorials, see the [entry page for the tutorials](tutorial/index.md).
This is useful if you consider that the installation has been successful. Or
you might continue to read the present Web page, and try to perform the speed
tests, as well as the other tests.

## How to make the other tests?

Let us pursue with the testing procedure. Go to the ~abinit/tests directory, and issue:
    
    make help

The workhorse script to run the tests is called runtests.py . It is very
flexible. A reasonable set of tests (those contained in the fast and v"x"
directories), can be run automatically by:
    
    ./runtests.py

or e.g.
    
    ./runtests.py -j4

(if you have 4 cores on your computer)

This is the recommended procedure for developers. In order to execute these
tests, you will need a larger disk space than for the simple installation of
the code (the total additional disk space required is on the order of 1GB).

The video below gives an overwiew of the command line options of `runtests.py`

[![asciicast](https://asciinema.org/a/40324.png)](https://asciinema.org/a/40324)

Let us now examine the different subdirectories.

~ABINIT/tests/fast/ is the simplest, and its content will be described in some
detail below. For tests of the parallel version see the directory
tests/paral/, as well as the ~abinit/doc/users/paral_use text file. For tests
of the response function features of abinit, and for tests of mrgddb and
anaddb, see the subdirectories tests/v2. The other directories tests/v3,
tests/v4, ... presents further tests of recently implemented features of ABINIT.

**1) tests/fast/** (for the sequential version only)

This subdirectory contains a basic set of tests of the code, aimed at testing
whether the code is coherent in time (successive versions), and exercising
several parts of the code. However, they do not examine its accuracy on
physical problems, mainly because the number of plane waves used is too small,
and some tests are not run to self-consistent convergence. 32 MB of memory
should be enough for these tests.

The input files for each of the tests can be found in the
~abinit/tests/fast/Input directory. At the bottom of the each of the input
file, some metadata is present. Such metadata mentions the executable to be
used (automatically), possibly the other input files (like pseudopotentials),
the output files to be analyzed, the admitted tolerances with respect to
reference output files, the author of the test, and a brief description of the test.

To run only the tests in this directory, simply issue, in the ~abinit/tests/ directory:
    
    ./runtests.py fast

It will create a directory named Test_suite. All the results will be in that
directory. The output files will be automatically compared, thanks to a 'diff'
command, to a set of reference files (in ~abinit/tests/fast/Refs/). The
corresponding difference files are prefixed by 'diff.'.

In addition to 'diff', there are two other levels of automatic analysis: one
based on a comparing tool called 'fldiff', producing 'fldiff.report' files,
and another where the output of 'fldiff' is further analyzed, and produce a
brief report called 'report'. The latter step is only performed in case all
the tests cases of one directory are performed (including the case where tests
of different directories are performed)

The one-line summaries produced by fldiff (see later) are compared with the
tolerances indicated in the input file (metadata added at the end of the input
file). This procedure produces a file called "report", in which there is a one
line assessment of the behaviour of each test: succeeded (everything is OK),
passed (the test is OK for users in production), passed marginally (the test
is within 1.5 of the usually accepted deviation, which is likely OK for most
applications - still to be improved by the development team, though), failed
(there is a problem, the deviation is usually not accepted). This is by far
the most convenient tool to analyze the automatic tests of abinit.

The vast majority of tests cases succeed or pass on all platforms that are
used by the developer team in Louvain-la-neuve. Some problems are mentioned in
the file ~abinit/KNOWN_PROBLEMS , Additionally, there might be specific
problems for some test case for some platforms, also mentioned in
~abinit/KNOWN_PROBLEMS. So, apart of the known problems, mentioned in this
file, the "report" file should mention, for each test case, only "succeeded" or "passed".

The comparing tool 'fldiff' -for 'floating diff'- performs in a more detailed
way the comparison of floating numbers between the output files and the
reference files than in the case of a 'diff' command. As used presently, for
each run inside one directory, one single file, called 'fldiff.report', will
be produced, and gather the analysis for all tests in that directory.

If for one test case, the two files differ by the number of lines, the
'fldiff.report' file will report that it cannot compare the two files. Usually
this problem will be seen at the level of 'command signs' appearing sometimes
in the first column of the output files, so a typical error message
(announcing something went wrong) will be:
    
    Case_1
    22
    The diff analysis cannot be pursued: the command sign differ.

By contrast, it will identify the floating numbers and ignore their
differences if they are within some prescribed tolerance, or if the difference
is not relevant. For example, it is able to ignore the differences in timings.
If everything goes fine for a test, fldiff should identify only the
differences in:

  * the dates of execution (possibly);
  * the version numbers (possibly);
  * the platform description (possibly);
  * the overall execution time (this is ALWAYS printed, even without differences).

So, a successful execution of one test case may be announced as follows in the fldiff.report file:
    
    Case_1
    2
    <  Version 8.0.6  of ABINIT
    >  Version 8.0.3  of ABINIT
    5
    <  Starting date: Mon 23 May 2016.
    >  Starting date: Mon  4 Apr 2016
    202
    < +overall time at end (sec): cpu=          7.1  wall=          8.0
    > +Overall time at end (sec): cpu=          7.3  wall=          8.0
    Summary Case_1: no significant difference has been found
    

The fldiff.report file will have one such section for each test_case that was
run. It begins with the number of the test case, then includes a few blocks of
three lines: the number of the line where something is happening, followed by
the content of the two lines. Finally, there is a one-line summary for each test case.

More information on the fldiff script can be found in the ~abinit/tests/Scripts/fldiff.pl file.

**2) tests/v1**

This directory contains tests built in the same spirit as those in the
test/fast directory, but that exercise other basic features, like the
treatment of metals, the GGA, the new pseudopotentials, the multi-dataset
mode, the cell parameters optimization, and the spatial symmetry groups
database. These were developed during the development time of the version 1 of
ABINIT. Of course, the automatic difference procedure only compares to recent
runs of the ABINIT code. As for the 'fast' test cases, the fldiff.report and
report files are also available. 64 MB of memory should be enough for these tests.

**3) tests/v2**

This directory contains tests built in the same spirit as those in the
tests/fast/ or tests/v1 directory, but that exercise features not present in
version 1 of the ABINIT package, mainly the response function features, and
the use of the mrgddb and anaddb codes. Again, the automatic difference
procedure only compares to recent runs of the ABINIT code. As for the 'fast'
test cases, the fldiff.report and report files are also available. 64 MB of
memory should be enough for these tests.

**4) tests/v3, tests/v4, tests/v5, tests/v6, tests/v67mbpt, tests/v7, tests/v8, tests/bigdft, 
    tests/etsf_io, tests/libxc, tests/wannier90**

These directories contain tests built in the same spirit as those in the
tests/fast/ directory, but that exercise features not present in the
Plane_Wave code, nor in version 1 or 2 of the ABINIT package, noticeably the
use of the GW code, the utilities Cut3d, AIM, .., the PAW ... . Or the
interfacing with fallbacks. Again, the automatic difference procedure only
compares to recent runs of the ABINIT code. As for the 'fast' test cases, the
fldiff.report and report files are also available. 64 MB of memory should be
enough for these tests.

**5) tests/paral/ and tests/mpiio/** (need MPI support)

This directory contains tests built in the same spirit as those in the
test/fast/ directory, but that exercise the parallel version of the ABINIT code.

The script runtests.py considers one of the different input files, and for
this file, it will use the parallel code with one processing node, then
perform different parallel runs with an increasing number of processing nodes,
as specified in the metadata contained in the input file. As for the other
series of test, the diff and the fldiff utilities are used automatically, and
fldiff.report and report files are produced automatically.

**6) tests/cpu ** (for the sequential version only)

This subdirectory contains the scripts, and input files needed for testing the
cpu time, either on progressively finer real space grids, or on progressively
bigger unit cells. Please read the README file of this directory. Also for
this suite of tests is activated with:
    
    make tests_cpu

Unlike in the previous case, many directories will be created (more than 10 in
the present version). Their name begins with the test name (A1, A2, A3, A4,
B1, B2, B3, B4, C3, D3), and is followed by the machine name and the date.
Inside these directories, many runs are done. There is a 'report' file that
summarizes the timing of the different runs, and there is a 'diff' file, that
compares these timings with the reference (output files from a PIV at 2.8 MHz, usually).

The structure of these tests is more complex than that of the test/fast/ and
test/v1/ directories. The tools used are the 'serie' scripts (serieA,serieB,
serieC and serieD) as well as the 'getrep' script. For an explanation, contact
the ABINIT group. For the largest tests (B and D series), up to 200MB of
central memory are required.

## Things that are NOT in the installation packages.

  * Many pseudopotentials are not in the installation package:
    The package contain several dozen pseudopotentials, for testing purposes, see
    ~abinit/tests/Psps_for_tests/. However, the largest set of pseudopotentials
    and PAW atomic datafiles can be found on the [ABINIT Atomic data pseudopoentials and PAW datasets Web
    page](http://www.abinit.org/psp-tables). Other
    pseudopotentials have been generated by many different users, and might be
    shared, but you might have to contact them.

  * The "fallbacks" are **not** in the installation package. They are downloaded from 
    the Web automatically at configure time. If you do not have the command "wget" 
    install on your machine, or if you do not have an access to the internet, you should disable all the fallbacks. 

  * The Web site <https://www.abinit.org> contains many other things, including links to the forum, 
    the mailing list, the ABINIT events ...

## For developers: how to modify the code?

### To modify a file (arguments unchanged).

If you want simply to modify the content of an existing file (e.g. one of the
Fortran files, located in one of the src/* directories), without changing its
list of arguments, and recompile the code, the procedure is quite simple.
Modify that file, and reissue:
    
    make

in ~abinit.

### To modify a file (arguments changed).

If you want to modify the content of an existing Fortran file as well as the
list of arguments, and recompile the code, the procedure is to modify that
file, and then issue, in the ~abinit directory:
    
    ./config/scripts/abilint . .
    make

!!! important

    Do not forget the two dots in the abilint command.

### To add a Fortran file.

If you want to add a new Fortran file, there is a difference whether you work
with or without the autotools.

In both cases, you must **add the file name in the abinit.src file of the
corresponding directory**. As an example, suppose that you want to add a new
routine named **blork.F90** in the directory ~abinit/src/65_nonlocal/. Then,
you must also declare it in the ~abinit/src/65_nonlocal/abinit.src file.

Note that the choice of directory is important. You should choose to put your
new routine in a directory that has a prefix number (e.g. 42 for 42_nlstrain)
higher than all the directories that contain routines you will use (except if
the called routines are in the same directory), and smaller than all the
directories that call your routine (except for routines that are in the same
directory as yours).

If you work with the autotools, you have now to reissue:
    
    ./config/scripts/makemake
    ./configure
    make

(see the section 2).

If you work without the autotools, **you must also modify by hand** (i.e. find
the places where the Fortran files are listed, and complete the list) the
Makefile.in of the directory where the new file has been added (e.g.
~abinit/src/65_nonlocal/Makefile.in).

Then issue:
    
    ./config/scripts/abilint . .
    make

in ~abinit/.

### To generate the source package.

If you want to produce the source package abinit-_x.y.z_.tar.gz, type:
    
    make dist

in ~abinit/.

Do not forget to change its name (e.g. add your name after _**x.y.z**_, to
identify that this is a modified version of ABINIT).
