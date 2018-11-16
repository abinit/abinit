---
authors: XG
---


# How to contribute to ABINIT

!!! Warning
    This file is obsolete, and should be reexamined/rewritten.
    By contrast, the documentation for developers on the [ABINIT Wiki](https://wiki.abinit.org/doku.php?id=developers:overview) 
    is up-to-date.

This page provides a description of the procedures followed for development of
the ABINIT package through collaboration of different groups of persons, based
in different places in the world. 

Any comment or suggestion to improve these procedures will be welcome!

[TOC]

* * *

## Introduction

The ABINIT package is aimed at being used by different groups of people,
without mandatory control by the main contributors of the ABINIT group. In the
same way, the ABINIT development project is fundamentally open to the
contributions of different persons, not located in Louvain-la-neuve or
Corning. These contributing persons are members _de facto_ of the ABINIT
group.

People using the code might consider adding their personal subroutines,
without trying to make them part of the official ABINIT package. However, this
has two drawbacks for them: in subsequent versions, their modifications will
not be incorporated, so that they might have to check and modify the interface
for each new version ; moreover, their contribution is not tested by other
users of the code, so that bugs might remain unnoticed. It is also nicer to
share the result of their coding efforts with other users of the code.

Of course, a collaborative effort has also some drawbacks. In particular, the
collaboration between distant developers should be carefully planned, since
orthogonal modifications of the same piece of code by two different people at
the same time is very likely to happen, generating "negative progress", i.e. a
large waste of time when synchronization is to be done. Also, it is required
to use a well-defined coding style, to provide test case files, and to comment
the modifications and additions as much as possible, in order to facilitate
the maintenance and the future modifications.

This document aims at defining the protocol to be followed to avoid "negative
progress" due to lack of synchronization. The analogy with the procedures to
be used for the parallelization of a code is obvious. The aim is that each
external 'node' does not waste its time, that communications are kept at the
lowest level possible, and that the final result is correct ! We will need
barriers for synchronization, and so on ...

The ability to incorporate the contributions of different groups in a
harmonious way might become a noticeable strength of the ABINIT project.

## Code repositories

Managing different versions of ABINIT is done thanks to a utility called
[**Bazaar**](http://www.bazaar-vcs.org/). For an introduction to this powerful
and versatile tool, you can have a look at our dedicated pages in the
Developers' corner of the ABINIT web site.

Thanks to Bazaar, the development of a project becomes completely transparent,
since all the changes in the files are registered: the latest version is of
course available, but one can come backwards in time to track bugs easily or
to know what anybody did. The development of different branches is also
managed, as well as the subsequent merging procedure. This feature is
important to allow development by many different people. The place where all
the files are stored, including their history, is called a **repository**.

For developers, using the repository is a privileged way of managing their
contributions, since the coordinator is automatically and instantly informed
of their progresses. They also benefit from regular backups of all the
contributions.

The ABINIT repository is divided into several **main categories**:

  * **abinit-release**, where go all contributions for the official release packages, these contributions consisting mainly in fixes and documentation;
  * **abinit-devel**, where all new and ongoing developments are managed;
  * **abinit-doc**, containing new and incomplete documentation, until it is consistent enough to be moved to the source code.

Each category is divided into **versions**, 3 digits for _abinit-release_ and
2 digits for _abinit-devel_. Then, in each category and for each version,
there is one **reference branch**, codename _merge_, which is the backbone of
the development effort. All concerned developers are supposed to use this
branch as a starting point for their tasks and to keep permanently in sync
with it. They may have at least one branch of their own.

## Basic philosophy

In the following, we will distinguish between **"debugging" contributions**
and **"development" contributions**.

_Debugging contributions_ are typically modifications of a few lines in one or
relatively few routines, needed for the code to work properly or to be
properly documented. Sometimes they are related to comments within routines or
corrections to documents. Usually, such modifications do not need any
synchronization and can be sent directly to the coordinator via email, in the
form of a patch. It is however possible to have them synchronized by recording
them inside the **abinit-release** category of the repository.

For the time being, Xavier
[&lt;xavier.gonze@uclouvain.be&gt;](mailto:xavier.gonze@uclouvain.be) should be
contacted.

_Development contributions_ usually involve the addition of new capabilities
to the code. Despite the use of Bazaar, synchronization with the coordinator
is **ALWAYS** needed: one has to make sure that nobody else is already in
charge of a similar project! The development contributions might be quite
local (basically adding one routine, called by a few lines from an existing
routine), or, on the contrary, involve modifications of many existing
routines. Even for the local type of modifications, discussion with the
coordinator is mandatory. Though Bazaar now deals nicely with conflicts at the
file level, it is of great importance to avoid semantic conflicts from the
very beginning.

The developer will be allocated a _development task_. Related to this task,
they will be free to code, experiment, debug, and check the result of their
work without any communication with the rest of the ABINIT group. The
allocation of the task has obviously to be done in coordination with the rest
of the group prior to the work, while the result of the development has to be
incorporated into a new official version. The task is thus also limited in
time.

The prior allocation and subsequent incorporation, taking into account the
possibility that many different developers work independently, must be done
centrally by one or more coordinators (at present Xavier but this might change
in the future), in order to guarantee the harmony, relevancy and consistency
of all contributions.

It will be the responsibility of the developer to make enough checks of the
correctness of their modifications or additions. The developer should provide
adequate documentation: basically the description of the input variables and
output data in the _abinit_help_ file, as well as a possible update of the
bibliography. They should also provide one or more tests to be added to the
standard suite of tests. This is needed to ensure that the transfer to the
official version of the code has been done properly, and also that the new
capabilities will be preserved by the subsequent modifications of the code.
Finally, they have to show that their modifications have not suppressed
existing capabilities of the codes, by running the set of already existing
automatic tests.

It will be the responsibility of the coordinator to transfer the result of the
development effort of each developer to an official version, in such a way
that the test is reproduced in a satisfactory way. The subsequent maintenance
will be automatically done by checking that the corresponding test is still
working despite modifications.

## Detailed protocol for the developer

**a.** The developer proposes a modification or addition to the code to the coordinator. In most cases, such a proposal will be warmly welcome! There might be some discussions possibly involving other scientists to improve the original proposal. The developers mailing list has been set-up for that purpose (see the "Community section" on the ABINIT web site). 

**b.** With respect to the latest official version, the developers define their task and make a list of routines they wish to heavily modify (and those they would like to add, but this is not really needed at that time). Following the current coding practice, adding new input variables would need rather small modifications of a very large number of routines (invars, abinit, gstate, brdmin, mover, scfcv_core, vtorho ...). In order to avoid the allocation of too large an area of development within the source code, which would make other development efforts impossible, it is assumed that the developer will make appropriate use of the different 'user' input variables. This appropriate use means that none of the aforementioned routines will be modified for that purpose. The easiest case to handle is when a routine or a set of routines are added to the package, and the interface to the rest of the code involves only a set of called subroutines the interface (list of arguments) should not change, and just one (modified) calling routine. Many development projects of this type can coexist. By contrast, if one of the development efforts involves major changes to major routines, it will prohibit the execution of another development project of the same type, taking into account the current structure of the code. 

**c.** Having gathered the different suggestions made on the basis of the current official version, the coordinator releases the next version, adequately prepared to allow the proposed developments by different developer groups, with a list of routines allocated to each group, the number of the test cases allocated, and, for information, the tasks that will be done. 

**d.** After installation on the development machine, and prior to any modification, the developer runs the internal tests, as well as the tests series _v1_ to _v5_ and the _paral_ (sequential cases) suite of tests. The "`make tests_dev`" command is perfectly suited for that purpose. Checking against previous reference files that the results are OK is always nice enough. The developer can replace the reference files provided with the official version by those produced on their machines, to ease further checking of the results. Then the development effort can begin. It is expected _a priori_ that the following files will be modified in any case: 

  * _abinit.src_ files in some _src/*_ directories, because the new routines must be listed there;
  * the _allvariables.html_ and some of the _var*.html_ files, because the new input variables must be listed and documented there (the 'user' variables are now used, but a final, aesthetic, name should be proposed);
  * the _README_ and _tests.cnf_ files in the subdirectory _v5_, as well as an additional input file (at least).

**IMPORTANT**  
  
During the development, only the allocated routines should be modified. This
is very important. Many others can be created.

In both cases, the developers **MUST** follow the current ABINIT coding style,
presented in the latest version of the document _rules_coding_, in the
_~abinit/doc/developers_ subdirectory. In particular, they should mention
their initials in the header of the new or modified routines. This will be
useful if somebody needs information about this routine some time later.

At the end of the development effort, it is mandatory that the developer runs
again the _tests_dev_ series, to be sure that the developments have not
spoiled some other feature of the code. This is easily done by issuing **`make
tests_dev`**. This command will run automatically the required tests, and
produce a file _summary_tests.tar.gz_ that should be send with the updated
routines and separate new tests. This is important, since the coordinator will
have to run the suite of tests as well, but at that time, having trouble with
a modification done by some developer would mean an important delay in the
delivery of the new official version, and thus a large waste of time. It is
desirable that the tests of the new feature do not last more than 1 minute on
Pentium III 1GHz, or about 20 seconds on other (faster) processors.

**e.** After the development has occurred, the developer prepares a gzipped (compressed) tar file with all the needed files (additional routines, modified routines, makemake, abinit_help, README, tests.cnf, ...), the _summary_tests.tar.gz_ file, and makes it available to the coordinator (by email, FTP or SSH). 

## Example

a. David proposes a nice addition to the ABINIT package. Doug and Xavier are
enthusiastic about it.

b. On the basis of version 5.3.4, David proposes a list of routines that
should be allocated to him.

c. Xavier delivers a 5.4.x version that has been adequately prepared for the
independent development that David wants to do.

d. David implements his addition and modifications on the basis of version
5.4.x, checks whether the suite of tests is still OK, and makes the files
available to Xavier.

e. Xavier thanks David very much, transfers the work to the official version
5.5.y and performs the final modification of the names of the (new) input
variables, including the propagation through the routines that were not
allocated to the developper.

