


This tutorial describes how to produce test-coverage metrics for C/C++ projects.  We start by instrumenting an autotools build to produce `gcov` output with `gcc`; we then generate test coverage artifacts by running our test suite, and finally we explore our results using `lcov`.

### Introduction

We'd like to assess the quality of the existing test suite for each Product Strategy projects.  A measurement of *test coverage* will tell us what part of the project's code is "covered" or exercised by its tests.  50% is a start; 80% would be a good goal for a neglected project; one rarely encounters 100% test coverage, but we'd like to get as close as we can.  Initially we'll use these findings to gain an overview of the test-quality of each project; ultimately these metrics will guide the improvement of our codebase, and enable us to monitor our progress using Jenkins and associated open-source tools.

### A Three-Part Process in Several Steps

We'll enable test-coverage of a particular C/C++ project in a three-part process:

* enabling a special build
* running the tests
* studying the output

Our first step will be to enable a special build.  Ultimately this just means adding a few flags to our `gcc` invocation, but it never seems so straightforward with autotools projects ;) .  The Product Strategy Quality team has made available a set of files to facilitate this build--please check these out now if you'd like to follow along.

    bzr co lp:~allanlesage/coverage-tutorial

Inspecting this archive you'll find 

* __`gcov.m4`__: an `autoconf` macro which will check for relevant tools
* __`Makefile.am.coverage`__, which includes our coverage-enabled `automake` targets
* an old revision of __`dbus-test-runner`__
* a copy of this tutorial

Before we start let's make sure we have `lcov` installed.  `lcov` incoporates the GNU tool `gcov` and produces pretty HTML reports of our coverage results.

    sudo apt-get lcov

(Note that the Debian `lcov` package includes `genhtml`.)  We're ready to begin.


##### 1.  Branch a project of interest

I heard you like testing, so for this example we'll test __dbus-test-runner__ which is a runner which runs tests against the dbus messaging system.  An old revision is included in our coverage-tutorial archive; if you want to start fresh (or on a different project), you would

    bzr branch lp:dbus-test-runner


##### 2.  Install the `gcov.m4` `autoconf` macro and invoke it

Open `gcov.m4`: here's where we're defining the compiler flags and checking for necessary tools.  Let's put this file where `autoconf` can find it.

Your project may already have a directory for `m4` macros:

    $ grep AC_CONFIG_MACRO_DIR *
    configure.ac:AC_CONFIG_MACRO_DIR([m4])

If you find such a directive, simply copy the `gcov.m4` file into the named directory.  If not, create an `m4` directory and copy the `gcov.m4` file there--and don't forget to include a reference to the directory in the body of `configure.ac` (it'll look like the `grep` result above).


##### 3.  Massage `configure.ac` to include our necessary coverage flags.

This is the essential move of our special build.  `gcc` supports coverage reporting with the addition of a few compiler flags.  Here they are in the `gcov.m4` file:

    # Remove all optimization flags from CFLAGS
    changequote({,})
    CFLAGS=`echo "$CFLAGS" | $SED -e 's/-O[0-9]*//g'`
    changequote([,])
    # Add the special gcc flags
    COVERAGE_CFLAGS="-O0 -fprofile-arcs -ftest-coverage"
    COVERAGE_CXXFLAGS="-O0 -fprofile-arcs -ftest-coverage"
    COVERAGE_LDFLAGS="-lgcov"

When the tests are run, `--fprofile-arcs` produces a tally of the execution of each *arc* of the code--one `.gcda` file for each source file.  The `--ftest-coverage` flag produces `.gcno` files, which link an *arc* to a source line so that we can see which lines were touched.  We'll watch these files being produced in a little bit; you can read about the flags at length in the [GNU documentation](http://gcc.gnu.org/onlinedocs/gcc/Debugging-Options.html).  Note that we also remove optimization with the `-O0` flag to get more precise results.

Now here's the step which requires knowledge both of your code and little bit of `autoconf`.  We've included the flags above (in step 3), but we'll now need to edit our `configure.ac` to make the flags available to our build.  For dbus-test-runner the diff looks like this:

    === modified file configure.ac
    --- configure.ac  2010-12-08 02:35:12 +0000
    +++ configure.ac  2011-12-06 21:42:04 +0000
    @@ -45,6 +45,16 @@
     AM_GLIB_GNU_GETTEXT
     
     ###########################
    +# gcov coverage reporting
    +###########################
    +
    +m4_include([m4/gcov.m4])
    +AC_TDD_GCOV
    +AC_SUBST(COVERAGE_CFLAGS)
    +AC_SUBST(COVERAGE_CXXFLAGS)
    +AC_SUBST(COVERAGE_LDFLAGS)
    +
    +###########################
     # Files
     ###########################

And then having added these flags to the build process, we need to actually *actually* add them to the build proper:

    === modified file src/Makefile.am
    --- src/Makefile.am     2009-12-07 21:00:43 +0000
    +++ src/Makefile.am     2011-12-06 21:42:04 +0000
    @@ -3,6 +3,8 @@
     
     dbus_test_runner_SOURCES = dbus-test-runner.c
     dbus_test_runner_CFLAGS  = $(DBUS_TEST_RUNNER_CFLAGS) \
    +                         $(COVERAGE_CFLAGS) \
                              -DDEFAULT_SESSION_CONF="\"$(datadir)/dbus-test-runner/session.conf\"" \
                              -Wall -Werror
     dbus_test_runner_LDADD   = $(DBUS_TEST_RUNNER_LIBS)
    +dbus_test_runner_LDFLAGS = $(COVERAGE_LDFLAGS)

Here we risk running afoul of the autotools "ancient ones".  This is an uncomplicated project--if you're not getting the results you want in the next step, I can recommend [*The Goat Book*](http://sourceware.org/autobook/autobook/autobook_40.html#SEC40) as a decent tutorial on how these macros work.

##### 4.  Verify that `autoconf.sh` `--enable-gcov` generates `.gcno` files

Having enabled our special build with the `gcc` compiler flags, let's verify that `autoconf` is generating the `.gcno` files we've asked for:

    $ ./autogen.sh --enable-gcov
    ...
    $ find -name *.gcno
    ./src/dbus_test_runner-dbus-test-runner.gcno

If all goes well we'll find a `.gcno` for each `.o` file compiled.  If not, then our flags aren't being respected--maybe something's wrong with our `configure.ac` addition in step 3.  It's important to have these checkpoints along the way to divide the autotools process into testable pieces.

##### 5.  Install `Makefile.am.coverage`

So now that we have the instrumentation we need, we'll use a special build to run the tests.  `Copy Makefile.am.coverage` into the top-level directory.

    cp ../Makefile.am.coverage .

Inspection shows that this file defines some extra make targets to generate coverage results using `lcov`.  We'll alert automake to this by simply including the file--here's the diff for dbus-test-runner:

    === modified file Makefile.am
    --- Makefile.am   2009-04-23 21:19:56 +0000
    +++ Makefile.am   2011-12-19 18:00:38 +0000
    @@ -1,1 +1,3 @@
     SUBDIRS = data src tests po
    +
    +include $(top_srcdir)/Makefile.am.coverage


##### 6.  Verify that make check generates `.gcda` files

Now we'll actually run the tests.

    $ make coverage-html
    ...
    $ find -name *.gcda
    ./src/dbus_test_runner-dbus-test-runner.gcda

Again these `.gcda` files represent a tally of 'arcs' through the code.  Note that we can generate these files not just while running tests, but also during normal execution--you can see the kinship with profiling.  `lcov` (via `gcov`) has used these tallies with our `.gcno` files to produce line-by-line results which we'll study in a moment.  Here's what the test coverage looks like for dbus-test-runner:

    Overall coverage rate:
      lines......: 78.8% (215 of 273 lines)
      functions..: 86.4% (19 of 22 functions)
      branches...: 61.2% (60 of 98 branches)

For a small project these figures show that we have a good test-suite to build on.  We're curious about what we're missing, though, so we'll investigate below.  Before you register a comment about "no coverage results", recognize that you may actually have zero code coverage.  But at least you're able to measure it, riiight?

##### 7.  Have a look at the `lcov` pages

`lcov` has produced a coverage-results directory at the top-level of our project; open index.html with a web browser and explore.

![lcov index](./lcov_index.jpg "lcov index")

At a glance our line coverage is 78.8%, and our "test-coverage progress bar" is yellow.  Be aware of the difference between the offered metrics:

* __line coverage__: how many lines have our tests touched
* __function coverage__: how many functions have our tests touched 
* __branch coverage__: for the graph which describes all possible paths of control through this file--especially through conditionals, e.g.--what percentage has our tests touched

In our opinion the killer feature of the `lcov` output is the source-file display, which shows which lines weren't touched by tests.  Drill into the source directory to see the results for a particular file.

![lcov detail](./lcov_detail.jpg "lcov detail")

Here in `dbus-test-runner.c` line 93 we've hit the line that tests the status of `G_IO_STATUS_NORMAL` 61 times, but never taken the path in the line below, for which `G_IO_STATUS_NORMAL` is false.


### Conclusion

We've presented an introduction to test coverage, featuring the `gcc` toolchain and autotools--as demonstrated with a set of utilities we've found useful here in the Quality group.

I hope it's obvious that this tutorial is just the beginning--not just of applying these methods to our code (we and you), but also of understanding what test coverage means to the improvement of Quality.  Our team has witnessed that having an observable measure like the `lcov` green bars encourages developers to increase test coverage.  However while chasing into a particularly narrow corner of their code, our developers have produced some confounding results which have caused us to doubt the accuracy of the available tools--I hope that Thomas Voss will follow with a post on some of his findings in the future.

Regardless, test coverage will be an important metric for our group for the coming cycles, and the bars will continue to grow greener as Spring approaches. . . .  Meanwhile I'll be interested to learn about your autotools/`gcc`/coverage experiences:

* What `lcov` alternatives have you tried?
* Which hallowed autotools rituals have I profaned in preparing this tutorial?
* What's the highest coverage percentage you've witnessed in a GPL project?
