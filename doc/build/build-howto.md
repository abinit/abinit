How to build ABINIT
===================

Introduction
------------

This document describes how to build ABINIT using the autotools-based
build framework.

If you find something out-of-date or incorrect, please send a message
to Yann Pouillon.



Build requirements
------------------

Before you try to build ABINIT, please make sure that you have installed
the correct tools, and they are correctly setup too.

### Unix/Linux ###

#### Hardware ####

 * 256 MB RAM (128 may be sufficient if you're very patient). For static
   builds, 512MB RAM minimum.
 * For debug builds: at least 2.5 GB free disk space (4 GB recommended).
 * For optimized builds: at least 400 MB free disk space (600 MB recommended).

#### Build Tools ####

 * A recent POSIX Shell
   The default /bin/sh that comes with some older unices (notably OSF/1)
   is known to be deficient. In that case, you should use 'make SHELL=ksh'
   for instance.
 * A C compiler
   GCC 3.2 or higher is recommended; egcs 1.0.3 (or higher), gcc 2.95.2
   (or higher) will also work; or your platform's native C compiler.
   Redhat 7.0 users, the compiler distributed with RH 7.0 is buggy, and it
   is recommended that you upgrade to the latest gcc 2.9x compiler
   (2.96-77 or later). You need the packages named *gcc* and *cpp*.
 * Perl 5.6 or higher
   Older Perl versions may work, though they have not been tested.
 * Python 2.2 or higher
   Older Python versions may work, though they have not been tested. Please
   note that this will become a true requirement soon.
 * GNU make 3.79.1 or higher
   Other flavours of make may work, yet we do not support them.

#### Optional Software ####

The following software are used only by developers. If you have no idea
what their purposes are, then don't worry. You don't need them to
build ABINIT.

 * Autoconf 2.59 or higher
   Autoconf is necessary if you want to hack on configure.in. Earlier versions
   might work, though they have not been tested.
 * Automake 1.9 or higher
   Automake is necessary to generate a *Makefile.in* from a *Makefile.am*.
   Autmoake 1.6 and earlier will likely not work.
 * Bazaar 1.4 or TLA 1.3.2 (or higher)
   TLA will allow you to access the latest version of the source.



### Windows ###

### Mac OS X ###

ABINIT builds for Mac OS X are Mach-O builds, built with gcc and Makefiles.
Doing builds on Mac OS X is therefore very similar to doing Unix builds.
You should be comfortable with the Terminal app (in /Applications/Utilities).
To install the software requirements, you will need to have admin privileges
on the machine.



#### Hardware ####

 * Mac OS X 10.1 or later
 * 2 GB or more free disk space



#### Build tools ####

 * The latest Developer Tools from Apple, if you don't have them already.
   You'll need a (free) ADC membership to download the tools. The December
   2002 tools or later are required to build ABINIT.

 * The current version of Fink (0.7.0 at time of writing). Note that there
   is conflict between Fink and Virex (both use /sw), so if you have Virex
   installed, you need to update to the latest version of Virex (7.2.1) and
   resolve the conflict as described in the Fink news (2003-04-16).

Earlier versions of Mac OS are not supported.



Get the source
--------------

The source package of ABINIT is usually available through the ABINIT website,
https://www.abinit.org/. Developers are strongly advised to use Bazaar
instead, as it will accelerate a lot the integration of their contributions.
For more details on Bazaar, please consult the ABINIT website
(https://www.abinit.org/bzr/) and Bazaar's website (http://bazaar-vcs.org/).



Configure build options
-----------------------

Running configure and make with the default options will not give you
a fully optimized build. In order to take full benefit from your machine,
you should use a config file. Please read these directions carefully before
building ABINIT.

### Using a per-user config file ###

Though it is possible to manually call "configure" with command-line
options, there is a better and more permanent way to tune the build
parameters. This is typically done by placing a config file in your
home directory: <tt>~/.abinit/build/*<hostname>*.ac</tt> .

This file is actually a Bourne shell-script which will be executed by
*configure* to set-up or override its default values. Since it is named
after the host it is stored on, you will be able to build ABINIT on several
machines sharing the same home directory.

You will find a self-documented tunable configuration for ABINIT in
~abinit/doc/build/config-template.ac8.



### Building with an OBJDIR ###

It is **highly recommended** that you use a separate object directory (OBJDIR)
when building ABINIT. This means that the source code and object files are not
intermingled, and that you can build ABINIT for multiple architectures from
the same source tree.

To build ABINIT using an OBJDIR, just create it, go inside and run *configure*
from there, e.g.:

<pre>
mkdir tmp-builddir ; cd tmp-builddir ; ../configure ; make
</pre>



### Other options ###

There are many other options recognized by the configure script. Some of
them are special-purpose options intended for embedders or other users,
and should not be used to build the full suite. The full list
can be obtained by running *./configure --help*.

Advice: if you can't figure out what a configure option does, don't use it!



Build & install
---------------

Building ABINIT is done by following the now traditional
*configure - make - make install* trilogy. To be more precise:

 * preliminary step: <tt>mkdir tmp-builddir ; cd tmp-builddir</tt>
 * 1. <tt>../configure</tt>
 * 2. <tt>make</tt>
 * 3. <tt>make install</tt>

A <tt>make check</tt> may be performed after <tt>make</tt>, at your option.



