The current configuration of Abinit requires the presence of
fallbacks, i.e. external packages that can be built and
bundled together for development and/or testing purposes when
full versions of these external packages are not available
on the system.

If you intend to perform production calculations with
optional features enabled, we highly recommend you to use
the full version of these external packages instead. For
more information on how to have them installed for Abinit,
please consult the Abinit Wiki at https://wiki.abinit.org/.

If you are currently working on a development/unstable
version of Abinit or making test calculations, you may use
the Abinit Fallbacks.

The Abinit Fallbacks are developed as an independent
and slowly-evolving project. You can download them from
https://launchpad.net/abinit-fallbacks, which is the
recommended way to access and build them, since it lets you
fully tune their build. However, for your convenience, we
provide them within the Abinit source tree, under the
fallbacks/ subdirectory.

In order to help you, we have created a script named
build-abinit-fallbacks.sh in the fallbacks/ subdirectory of
the Abinit build directory. To build the fallbacks required
for the current Abinit configuration, you can run this
script through the following commands:

cd fallbacks
/bin/sh ./build-abinit-fallbacks.sh

The configure script will now fail. You will need to run it
again with the --with-fallbacks option after having
properly installed the required fallbacks. You will then be
able to reuse them for many subsequent builds of Abinit.
