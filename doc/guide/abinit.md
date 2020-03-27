---
authors: DCA,  XG
---

# ABINIT, the main code  

This document explains the i/o parameters and format needed for the main code (abinit) in the ABINIT package.  

The new user is advised to read first the [new user's guide](..),
before reading the present file. 
It will be easier to discover the present file with the help of the [[tutorial:index|tutorial]].

When the user will be sufficiently familiarized with ABINIT, reading the
~abinit/doc/users/tuning.txt file might be useful (this file, as many
additional documentation files, is not available on the Web, but is available in the package). 
For calculating response properties using abinit, the complementary [[help:respfn]] is needed.

<a id="intro"></a>
## 1 How to run the code?

### 1.1 Introducing the 'files' file
  
Given an input file (parameters described below) and the required
pseudopotential files, the user must create a "files" file which lists names
for the files the job will require, including the main input file, the main
output file, root names for other input, output, or temporary files, and the
names of different pseudopotential or PAW atomic data files (one per line).
The files file (called for example ab.files) could look like:
    
       ab_in
       ab_out
       abi
       abo
       tmp
       Si-GGA.psp8
       O-GGA.psp8

In this example:  

* The main input file is called "ab_in".  
* The main output will be put into the file called "ab_out".  
* The name of input wavefunctions (if any) will be built from the root "abi"
    (namely abi_WFK, see later).  
* The output wavefunctions will be written to abo_WFK. Other output files
  might be build from this root.  
* The temporary files will have a name that use the root "tmp". (for example tmp_STATUS).  
* The pseudopotentials needed for this job are "Si-GGA.psp8" and "O-GGA.psp8".  
  
Other examples are given in the subdirectories of the ~abinit/tests directory.

!!! important 

    The maximal length of names for the main input or output files is presently
    132 characters. It is 112 characters for the root strings, since they will be
    supplemented by different character strings.

<a id="exec"></a>
### 1.2 Running the code
  
The main executable file is called abinit. Supposing that the "files" file is
called ab.files, and that the executable is placed in your working directory,
abinit is run interactively (in Unix) with the command

    abinit < ab.files >& log
  
or, in the background, with the command

    abinit < ab.files >& log &

where standard out and standard error are piped to the log file called "log"
(piping the standard error, thanks to the '&' sign placed after '>' is
**really important** for the analysis of eventual failures, when not due to
ABINIT, but to other sources, like disk full problem...). The user can
specify any names he/she wishes for any of these files. Variations of the
above commands could be needed, depending on the flavor of UNIX that is used
on the platform that is considered for running the code.  

<a id="2"></a>
## 2 The underlying theoretical framework and algorithms

The methods employed in this computer code to solve the electronic structure
problem are described in part in different review papers as well as research
papers. The code is an implementation of the Local Density Approximation to
the Density Functional Theory, based upon a plane wave basis set and separable
pseudopotentials. The iterative minimization algorithm is a combination of
fixed potential preconditioned conjugate gradient optimization of wavefunction
and a choice of different algorithms for the update of the potential, one of
which is a potential-based conjugate gradient algorithm.

The representation of potential, density and wavefunctions in real space will
be done on a regular 3D grid of points. Its spacing will be determined by the
cut-off energy (see the input variable [[ecut]]) of the planewave basis in
reciprocal space. This grid of points will also be the starting point of Fast
Fourier Transforms between real and reciprocal space. The number of such
points, called [[ngfft]], should be sufficiently large for adequate
representation of the functions, but not too large, for reasons of
computational efficiency. The trade-off between accuracy and computational
efficiency is present in many places of the code, and addressed briefly at the
end of the present help file.

We recommend a good introduction to many different concepts valid for this
code, available in this [[cite:Payne1992|Reviews of Modern Physics article]].
Note that this paper does NOT reflect the present status of the code. 
ABINIT is closer in spirit to the [[cite:Kresse1996|paper]] of Kresse and Furthmuller.
If you have never used another electronic structure code or a Quantum
Chemistry package, you should browse through the Chaps. 1 to 13, and
appendices L and M of [[cite:Martin2004|this book]] by R. M. Martin

<a id="input"></a>
## 3 The input file

### 3.1 Format of the input file
  
Note that this input file was called ab_in in the example of section 1.1.  
We first explain the content of the input file without use of the 
"multi-dataset" possibility (that will be explained in section 3.3).

The parameters are input to the code from a single input file. Each parameter
value is provided by giving the name of the input variable and then placing
the numerical value(s) beside the name, separated by one or more spaces, or
even an equal sign (equal signs are replaced by blanks by the parser).
Depending on the input variable, the numerical value may be an integer or a
real number (internal representation as double precision number), and may
actually represent an array of values. If it represents an array, the next set
of numbers separated by spaces are taken as the values for the array.

  * Do NOT separate a minus sign from the number to which it applies. 
  * Do NOT use tabs. 
  * NOTE THAT NO LINE OF THE INPUT FILE MAY EXCEED 132 CHARACTERS. 
    That is, only the first 132 characters of each line of the input file 
    will be read and parsed for input variables and their values. 

The names of all the parameters can be found in the [[varset:allvars|input variable database]]. 
The list of input variables present in the latter file links them to their definitions, contained 
in different "variable set" files, some of which are listed here:

  * Basic variables, [[varset:basic]]
  * Files handling variables, [[varset:files]]
  * Ground-state calculation variables, [[varset:gstate]]
  * GW variables, [[varset:gw]]
  * Parallelisation variables, [[varset:paral]]
  * Density Functional Perturbation Theory variables, [[varset:dfpt]]

In the actual input file, these parameters may be given in any desired order,
and more than one may be given per line. Spaces are used to separate values
and additional spaces are ignored.  
An as example of input, the parameter for length scales is called [[acell]]
and is an array [[acell]](3) for the lengths of the primitive translations in
Bohr atomic units. To input a typical Si diamond lattice one would have the line

    acell 10.25311 10.25311 10.25311

in the input file. This may equivalently be written

    acell 3*10.25311

and will still be parsed correctly: it is equivalent to the above line. Even

    acell *10.25311

will work. In the latter case the '*' sign means that the parser should use
the given value to fill the array, by repeating it as many time as needed.  
Multiple spaces are ignored, as is any text which does not contain the
character strings which correspond to some input parameters. In case of
arrays, only the needed numbers will be considered, and the eventual numbers
after those needed will also be ignored. For example,

    natom 3           # This gives the number of atoms  
    typat 1 1 2 2 3   # typat(1:natom) gives the type of each atom: only  
                      # the first three data are read, since [[natom]]=3 `  

A given variable is identified by the parser by having at least one blank
before it and after it (again, multiple blanks are irrelevant).  
ABINIT has also some (very limited) interpretor capabilities:

  * It can identify one slash sign (/) being placed between two numbers 
    (without a separating blank) as being the definition of a fraction 
    (e.g. 1/3 will be interpreted as 0.33333333333333d0); 

  * It can identify sqrt(...) or -sqrt(...) as being the definition of a square root, 
    when applied to one valid number - also without a separating blank - 
    (e.g. -sqrt(0.75) will be interpreted as -0.8660254038d0); 

  * Note, however, that these capabilities are NOT recursive. 
    At most, a sqrt identifier can contain an expression that uses a fraction (e.g. sqrt(3/4) is OK), 
    but two fractions (or two sqrt) cannot be used in one expression, and a sqrt cannot be present 
    in the numerator or denominator of a fraction. 

Comments should be placed to the right of the comment characters # or !;
anything to the right of a "#" or a "!" on any line is simply ignored by the
parser. Additional text, not preceded by a "#" or a "!" would not otherwise
cause trouble unless the text inadvertently contained character strings which
were the same as variable names (e.g. [[acell]]). The characters "#" or "!"
can also be used to "store" old values of variables or place anything else of
convenience into the file in such a way as to be ignored by the parser when
the data is read.  

Case is irrelevant as the entire input string is mapped to upper case before
parsing, to remove case sensitivity.  
More than one parameter per line may be given. If a given parameter name is
given more than once in the input file, an error message is printed, and the code stops.

External input files can be included with the syntax:

    include "geometry.inc"

where geometry.inc gives the crystalline structure in the Abinit format:

    cat geometry.in

    # Si in diamond structure
    acell 3*10.25
    rprim   
      0.0 0.5 0.5  
      0.5 0.0 0.5  
      0.5 0.5 0.0
    natom  2
    ntypat 1
    typat  2*1
    xred   0.00  0.00  0.00
           0.25  0.25  0.25
    znucl 14.0


<a id="parameters"></a>
### 3.2 More about ABINIT input variables
  
In each section of the ABINIT input variables files, a generic information on
the input variable is given: a **mnemonics**, possibly some
**characteristics**, the **variable type** (integer, real, string), and the
**default value**. Then, follows the description of the variable.

The **characteristics** can be one of the following: 

-   **DEVELOP**,
-   **NO_MULTI**
-   **INTERNAL_ONLY**
-   **INPUT_ONLY**
-   **EVOLVING**
-   **ENERGY**
-   **LENGTH**
-   **MAGNETIC FIELD**

#### Physical information

The **ENERGY**, **LENGTH** and **MAGNETIC FIELD** characteristics indicate
that the physical meaning of the variable is known by ABINIT, so that ABINIT
can treat its physical dimensions with some units.

The use of the atomic unit system (e.g. the Hartree for energy, about 27.211
eV, and the Bohr for lengths about 0.529 Angstroms) is strictly enforced
within the code. However, the dimension of some input variables can be
specified and read correctly. 

At present, this applies to three types of
variables: those that have the dimension of an energy, those that have a
dimension of length, and those that have a dimension of magnetic field. The
first class of variables have the characteristics **ENERGY**, and can be
specified in atomic units (Hartree), or electron-volts, or Rydbergs, or even Kelvin. 

The second class of variables have the characteristics **LENGTH**,
and can be specified in atomic units (Bohr), nm (nanometer) and angstrom. 
The third class of
variables have the characteristics **MAGNETIC FIELD**, and can be
specified in atomic units and Tesla. The abinit parser recognize a dimension
if it is specified after the list of numbers following the input variable
keyword, in the input file. The specification can be upper or lower case, or a
mix thereof. Here is the list of recognized chains of characters:

  * Ry or Rydberg or Rydbergs --> Rydberg (for energies) 
  * eV --> electron-volts (for energies) 
  * K  --> Kelvin (for energies) 
  * Angstr --> Angstrom (for lengths) 
  * nm --> nanometer (for lengths) 
  * T or Tesla --> Tesla (for magnetic fields) 
  * S or Sec or Second --> second (for time) 

Other character chains, like "au" (for
atomic units) or "Hartree", or "Bohr" are not recognized, but make the parser
choose (by default) atomic units, which is the correct behaviour. Example:
    
        acell 8 8 8 angstrom
        ecut 8 Ry
        tsmear 1000 K

or
    
         acell 3*10 Bohr  ecut 270 eV  tsmear 0.01

The use of the atomic units is mandatory for other dimensioned input
variables, like the tolerance on forces ([[toldff]]), parameters that define
an 'object' ([[objaax]], [[objaax]], [[objbax]], [[objatr]], [[objbtr]]), and
the initial velocity of atoms ([[vel]] if needed).

The initial atomic positions can be input in Bohr or Angstrom through
[[xcart]], but also, independently, in Angstrom through [[xangst]], or
even in reduced coordinates, through [[xred]].

#### Flow information

Most of the variables can be used in the multi-dataset mode (see section 3.3),
but those that must have a unique value throughout all the datasets are
signaled with the indication **NO_MULTI**.

Some of the input variables, with characteristics **INPUT_ONLY** are only
used by the parser, to initialize other input variables, but are not
transmitted inside the code, beyond the parser. In particular, they are not
echoed in the output file.

At variance, some internal variables, with characteristics **INTERNAL_ONLY** 
are documented in the help files, but are not accessible as input variables.
The documentation is provided because such variables are sometimes mentioned
in the output file.

Most of the input variables do not change while a run is performed. Some of
them, by contrast, may evolve, like the atomic positions, the atomic
velocities, the cell shape, and the occupation numbers. Their echo, after the
run has proceeded, will of course differ from their input value. They are
signaled by the indication **EVOLVING**.

#### Other information

**DEVELOP** refers to input variables that are not used in production
runs, but have been introduced during development time, of a feature that is
likely not finalized. For non ABINIT developers, it is strongly advised to skip them.

In addition to giving the input variables, the input file can be useful for
another purpose: placing the word " **exit** " on the top line will cause the
job to end smoothly on the very next iteration, if the [[chkexit]] input
variable is non-zero. This functions because the program closes and reopens
the input file on every iteration and checks the top line for the keyword
"exit". THE WORD MUST BE PLACED WITH SPACES (BLANKS) ON BOTH SIDES. Thus
placing exit on the top line of the input file WHILE THE JOB IS ALREADY
RUNNING will force the job to end smoothly on the very next iteration. On some
machines, this does not work always (we do not know why...). Another
possibility is offered: one can create a file named "abinit.exit" in the
directory where the job was started. The code should also smoothly end. In
both cases, the stop is not immediate. It can take a significant fraction
(about 20% at most) of one SCF step to execute properly the instruction still needed.

<a id="multidatasets"></a>
### 3.3 The multi-dataset mode
  
Until now, we have assumed that the user wants to make computations
corresponding to one set of data: for example, determination of the total
energy for some geometry, with some set of plane waves and some set of k-points.

It is often needed to redo the calculations for different values of some
parameter, letting all the other things equal. As typical examples, we have
convergence studies needed to determine which cut-off energy gives the needed
accuracy. In other cases, one makes chains of calculations in order to compute
the band structure: first a self-consistent calculation of the density and
potential, then the eigenenergy computation along different lines.

For that purpose, the **multi-dataset mode** has been implemented.

It allows the code to treat, in one run, different sets of data, and to chain
them. The number of datasets to be treated is specified by the variable
[[ndtset]], while the indices of the datasets (by default 1, 2, 3, and so on)
can be eventually provided by the array [[jdtset]].

For each dataset to be treated, characterized by some index, each input
variable will determined by the following **rules** (actually, it is easier to
understand when one looks at examples, see below):

  * (1) ABINIT looks whether the variable name (e.g. [[ecut]] ), appended with the index 
    of the dataset (e.g. [[jdtset]]=2), exists (e.g. "ecut2" ). It will take the data that follows this keyword, if it exists.

  * (2) If this modified variable name does not exist, it will look whether a metacharacter, 
    a series or a double-loop data set has been defined, see sections 3.4 or 3.5.

  * (3) If the variable name appended with the index of the dataset does not exist, and 
    if there is no series nor double-loop dataset for this keyword, it looks for an occurrence 
    of the variable name without any index appended, and take the corresponding data. (This corresponds to the single dataset mode)

  * (4) If such occurrences do not exist, it takes the default value. (Also, similar to the single dataset mode)
    
        ---------------
    
        1st example.
    
        ndtset   2
         acell   8 8 8
          ecut1  10
          ecut2  15

means that there are 2 datasets: a first in which
    
         acell 8 8 8  ecut 10 

has to be used, and a second in which
    
         acell 8 8 8  ecut 15

has to be used.
    
        ------------------
    
        2nd example
    
        ndtset 2     jdtset 4 5
    
        acell   8 8 8
        acell5 10 10 10
        ecut1  10
        ecut2  15
        ecut3  20
        ecut4  25
        ecut5  30
    
this means that there are still two datasets, but now characterized by the
indices 4 and 5, so that the first run will use the generic "acell", and "ecut4":
    
         acell 8 8 8 ecut 25

and the second run will use "acell5" and "ecut5":
    
         acell 10 10 10 ecut 30 

Note that ecut1, ecut2 and ecut3 are not used.

<a id="series"></a>
### 3.4 Defining a series
  
Rule (2) is split in three parts: (2a), (2b) and (2c). Series relate with (2b):

(2b) If the variable name appended with the index of the dataset does not
exist, the code looks whether a series has been defined for this keyword.

There are two kinds of series:

  * arithmetic series (constant _increment_ between terms of the series) 

  * geometric series (constant _ratio_ between terms of the series) 

The first term of the series is defined by the keyword appended with a colon
(e.g. **ecut:** ), while the increment of an arithmetic series is defined by
the keyword appended with a plus (e.g. **ecut+** ), and the factor of a
geometric series is defined by the keyword appended with a times (e.g. **ecut*** ).

If the index of the dataset is 1, the first term of the series is used, while
for index N, the appropriate input data is obtained by considering the Nth
term of the series.
    
     ------------------
    
     3rd example
    
       ndtset 6
       ecut1 10
       ecut2 15
       ecut3 20
       ecut4 25
       ecut5 30
       ecut6 35

is equivalent to
    
        ndtset 6 ecut: 10 ecut+ 5

In both cases, there are six datasets, with increasing values of [[ecut]].

<a id="loop"></a>
### 3.5 Defining a double loop dataset
  
To define a double loop dataset, one has first to define the upper limit of
two loop counters, thanks to the variable [[udtset]]. The inner loop will
execute from 1 to [[udtset]](2), and the outer loop will execute from 1 to
[[udtset]](1). Note that the largest value for [[udtset]](1) is presently 999,
while it is 9 for [[udtset]](2) (so, only the last digit for the inner loop).

The value of [[ndtset]] must be coherent with [[udtset]] (it must equal the
product `udtset(1) * udtset(2)`).

A dataset index is created by the concatenation of the outer loop index and
the inner loop index.  
For example, if [[udtset]](1) is 2 and [[udtset]](2) is 4, the index will
assume the following values: `11, 12, 13, 14, 21, 22, 23, and 24`.

Independently of the use of [[udtset]], rules (2a) and (2c) will be used to
define the value of an input variable:

(2a) The question mark " **?** " can be used as a metacharacter, replacing any
digit from 1 to 9, to define an index of a dataset.  
For example, **ecut?** 1 means that the input value that follows it can be
used for [[ecut]] for the datasets `01, 11, 21, 31, 41, 51, 61, 71, 81, and 91`.

(2c) If the variable name appended with the index of the dataset does not
exist, the code looks whether a double-loop series has been defined for this
keyword. Series can be defined for the inner loop index or the outer loop
index. Two signs will be appended to the variable name (instead of one in the
simple series case). One of these signs must be a question mark " **?** ",
again used as a metacharacter able to assume the values 1 to 9.  
If it is found in the first of the two positions, it means that the series
does not care about the outer loop index (so the values generated are equal
for all outer loop index values). If it is found in the second of the two
positions, the series does not care about the inner loop index. The other sign
can be a colon, a plus or a times, as in the case of the series defined in
(2a), with the same meaning.

Rule (1) has precedence over them, they have precedence over rules (3) or (4),
rule (2a) has precedence over rules (2b) or (2c) and the two latter cannot be
used simultaneously for the same variable.
    
        ------------------
    
        4th example
        ndtset 6    udtset 2 3
        acell1?  10 10 10
        acell2?  15 15 15
        ecut?: 5    ecut?+ 1
     
is equivalent to
    
        ndtset 6     jdtset 11 12 13  21 22 23
        acell11  10 10 10     ecut11 5
        acell12  10 10 10     ecut12 6
        acell13  10 10 10     ecut13 7
        acell21  15 15 15     ecut21 5
        acell22  15 15 15     ecut22 6
        acell23  15 15 15     ecut23 7
     

!!! tip

    More examples can be found in the directory ~abinit/tests/v1, cases 59 and later.

<a id="filenames-multidataset"></a>
### 3.6 File names in the multi-dataset mode
  
The root names for input and output files (potential, density, wavefunctions
and so on) will receive an appendix: **_DS** followed by the index of the
dataset. See section 4.

The  **get** variables can be used to chain the calculations.

Let us mention a few of them: [[getwfk]], [[getwfq]], [[getddk]], [[get1wf]],
[[getden]], [[getcell]], [[getxred]] and [[getxcart]].

  * [[getwfk]] allows to take the output wavefunctions of a previous dataset and use them as input wavefunctions 
  * [[getwfq]], [[getddk]] and [[get1wf]] do similar things for response function calculations 
  * [[getden]] does the same for the density; [[getcell]] does the same for [[acell]] and [[rprim]] 
  * [[getxred]] and [[getxcart]] do the same for the atomic positions, either in reduced coordinates, or in cartesian coordinates. 

The different variables corresponding to each dataset are echoed using the
same indexing convention as for the input step. For the last echo of the code
variables, some output variables are also summarized, using the same conventions:

  * **etotal** (total energy) 
  * **fcart** (cartesian forces) 
  * **strten** (the stress tensor). 

<a id="files-file"></a>
## 4 More detailed presentation of the files file
  
Note: _This "files" file is called _ab.files_ in section 1._

As mentioned in section 1 (you might read it again if needed), the "files"
file contains the file names or root names needed to build file names. These
are listed below: there are 5 names or root names for input, output and
temporaries, and then a list of pseudopotentials (one per line). These names
may be provided from unit 05 interactively during the run but are more
typically provided by piping from a file in Unix (the "files" file).

**ab_in**  
Filename of file containing the input data, described in the preceding sections.

**ab_out**  
Filename of the main file in which formatted output will be placed (the main
output file). Error messages and other diagnostics will NOT be placed in this
file, but sent to unit 06 (terminal or log file); the unit 06 output can be
ignored unless something goes wrong. The code repeats a lot of information to
both unit 06 and to the main output file. The unit 06 output is intended to be
discarded if the run completes successfully, with the main output file keeping
the record of the run in a nicer looking format.

**abi**  
The other files READ by the code will have a name that is constructed from the
root "abi". This apply to optionally read wavefunction, density or potential
files. In the multi-dataset mode, this root will be complemented by **_DS**
and the dataset index. The list of possible input files, with their name
created from the root 'abi', is the following (a similar list exist when 
**_DS** and the dataset index are appended to 'abi'):

  * **abi_WFK**   
filename of file containing input wavefunction coefficients created from an
earlier run (with [[nqpt]]=0). Will be opened and read if [[irdwfk]] is 1.
The wavefunction file is unformatted and can be very large. **Warning**: in
the multi dataset mode, if getwfk is non-zero, a wavefunction file build from
**abo** will be read.

  * **abi_WFQ**   
filename of file containing input wavefunction coefficients created from an
earlier run (with [[nqpt]]=1), as needed for response function calculations.
The wavefunction file is unformatted and can be very large. **Warning**: in
the multi dataset mode, if getwfk is non-zero, a wavefunction file build from
**abo** will be read.

  * **abi_1WFxx**   
filename of file containing input first-order wavefunctions created from an
earlier RF run. xx is the index of the perturbation

  * **abi_DEN**   
filename of file containing density created from an earlier run. See
explanations related to negative values of [[iscf]]. This file is also
unformatted. **Warning**: in the multi dataset mode, if getwfk is non-zero, a
density file build from **abo** will be read.

  * **abi_HES**   
filename of file containing an approximate hessian, for eventual
(re)initialisation of Broyden minimisation. See brdmin.F90 routine. The use of
[[restartxf]] is preferred.

**abo**  
Except "ab_out" and "log", the other files WRITTEN by the code will have a
name that is constructed from the root "abo". This apply to optionally written
wavefunction, density, potential, or density of states files. In the multi-
dataset mode, this root will be complemented by **_DS** and the dataset
index. Also in the multi-dataset mode, the root "abo" can be used to build the
name of **input** files, thanks to the 'get' variables. The list of possible
output files, with their name created from the root 'abo' is the following (a
similar list exists when **_DS** and the dataset index are appended to 'abo'):

  * **abo_WFK**   
Filename of file containing output wavefunction coefficients, if [[nqpt]]=0.
The wavefunction file is unformatted and can be very large.

  * **abo_WFQ**   
Same as **abo_WFK**, but for the case [[nqpt]]=1. The wavefunctions are
always output, either with the name **abo_WFK**, or with the name
**abo_WFQ**.

  * **abo_1WFxx**   
Same as **abo_WFK**, but for first-order wavefunctions, xx is the index of
the perturbation, see the section [[help:respfn#6.3|section 6.3]] of the [[help:respfn]].

  * **abo_DDB**   
The derivative database, produced by a response-function dataset, see [[help:respfn#ddb|this section]]
of the respfn help file.

  * **abo_DEN**   
filename of file containing density, in the case [[ionmov]]=0. See the keyword
[[prtden]]. This file is unformatted, but can be read by cut3d.

  * **abo_TIMx_DEN**   
filenames of files containing density, in the case [[ionmov]]/=0. The value of
"x" after " **TIM** " is described hereafter. See the keyword [[prtden]]. This
file is unformatted, but can be read by cut3d.

  * **abo_POT**   
filename of file containing Kohn-Sham potential See the keyword [[prtpot]].
This file is unformatted, but can be read by cut3d.

  * **abo_TIMx_POT**   
filenames of files containing Kohn-Sham potential in the case [[ionmov]]/=0.
The value of "x" after "TIM" is described hereafter. See the keyword
[[prtpot]]. This file is unformatted, but can be read by cut3d.

  * **abo_DOS**   
filename of file containing density of states. See the keyword [[prtdos]].
This file is formatted.

  * **abo_TIMx_DOS**   
filenames of files containing the density of states in the case [[prtdos]]=2
and [[ionmov]]=1 or 2. The value of "x" after "TIM" is described hereafter.
See also the keyword [[prtdos]]. This file is formatted.

  * **abo_GEO**   
filename of file containing the geometrical analysis (bond lengths and bond
angles) in the case [[ionmov]]=0. See the keyword [[prtgeo]]. This file is
formatted.

  * **abo_TIMx_GEO**   
filenames of files containing the geometrical analysis (bond lengths and bond
angles) in the case [[ionmov]]=1 or 2. The value of "x" after "TIM" is
described hereafter. See also the keyword [[prtgeo]]. This file is formatted.

  * **abo_KSS**   
filename of file containing output wavefunction coefficients, if
[[nbandkss]]/=0. This wavefunction file is unformatted and can be very large.
Its purpose is to start a GW calculation using M.Torrent's code. A different
format than for **abo_WFK** is used, see the file
~abinit/doc/developers/format_KSS.txt.

  * **abo_EIG**   
A file containing the electronic eigenvalues, for subsequent plotting of band
structure.

When [[ionmov]]/=0, the **POT**, **DEN**, or **GEO** files are output each
time that a SCF cycle is finished. The " **x** " of **TIMx** aims at giving
each of these files a different name. It is attributed as follows:  
\- case ionmov==1: there is an initialization phase, that takes 4 calls to
the SCF calculation. The value of x will be A, B, C, and D. Then, x will be 1,
2, 3 ..., actually in agreement with the value of itime (see the keyword
[[ntime]])  
\- other ionmov cases: the initialisation phase take only one SCF call. The
value of x will be 0 for that call. Then, the value of x is 1, 2, 3... in
agreement with the value of itime (see the keyword [[ntime]])

**tmp**  
The temporary files created by the codes will have a name that is constructed
from the root " **tmp** ". tmp should usually be chosen such as to give access
to a disk of the machine that is running the job, not a remote (NFS) disk.
Under Unix, the name might be something like `/tmp/user_name/temp`. As an
example, **tmp_STATUS**  
gives the status of advancement of the calculation, and is updated very
frequently

**psp1**  
filename of first pseudopotential input file. The pseudopotential data files
are formatted. There must be as many filenames provided sequentially here as
there are types of atoms in the system, and the order in which the names are
given establishes the identity of the atoms in the unit cell. (psp2, psp3, ...)

<a id="5"></a>
## 5 The pseudopotential files and PAW atomic data files
  
Actually, no real understanding of these files is needed to run the code. The
recommended pseudopotentials can be downloaded from the ABINIT Web site at
[[https://www.abinit.org/psp-tables]]. Documentation is
provided there as well as in the dedicated [[topic:PseudosPAW]]. Note that it
is not possible to mix norm-conserving pseudopotentials and PAW atomic data
sets in the same run. Also, every such file has been generated for a
particular choice of the exchange-correlation functional [[ixc]]. It is in
principle incorrect to use a pseudopotential (or PAW data) with another
exchange-correlation functional than the one it has been generated for, but
ABINIT will only send a warning.

For different other reasons, it might nevertheless be useful to be able to
grasp some information from the file. For norm-conserving pseudopotentials
different format are possible (labelled 1 to 8 presently). The associated
internal variable is called pspcod. Information on the header of these
pseudopotential files can be found in the abinit wiki at
[[https://wiki.abinit.org/doku.php?id=developers:pseudos]], that you should
read now (do not pursue with the description of each format, though).

## 6 The different output files
  
Explanation of the output from the code

Output from the code goes to several places listed below.

<a id="logfile"></a>
### 6.1 The log file
  
The "log" file (this is the standard UNIX output file, and corresponds to
Fortran unit number 06): a file which echoes the values of the input
parameters and describes various steps of the calculation, typically in much
more detail than is desired as a permanent record of the run. This log file is
intended to be informative in case of an error or for a fuller description of
the run. For a successful run the user will generally delete the log file
afterwards. There are four types of exception messages: **ERROR**, **BUG**,
**WARNING** and **COMMENT** messages.

**ERROR** and **BUG** messages cause the code to stop, immediately, or after a
very small delay. An **ERROR** is attributed to the user, while a **BUG** is
attributed to the developer.

A **WARNING** message indicates that something happened that is not as
expected, but this something is not so important as to make the code stop. A
**COMMENT** message gives some information to the user, concerning something
unusual. None of them should appear when the run is completely normal.

After a run is completed, always have a look at the end of the log file, to
see whether an **ERROR** or a **BUG** occurred.  

Also, the code gives the number of **WARNING** or **COMMENT** it issued. It is
advised to read at least the **WARNING** messages, during the first month of
ABINIT use.

<a id="outputfile"></a>
### 6.2 The main output file
  
The **main output file** is a formatted output file to be kept as the permanent record of the run.

Note that it is expected **not** to exist at the beginning of the run:  
If a file with the name specified in the "files" file already exists, the code
will generate, from the given one, another name, appended with **.A**. If
this new name already exists, it will try to append **.B**, and so on, until **.Z**.  
Then, the code stops, and asks you to clean the directory.

The **main output file** starts with a heading:

  * version number and specified platform 
  * copyright notice and distribution licence 
  * date 
  * echo of "files" file (except pseudopotential name) 

Then, for each dataset, it reports the point symmetry group and Bravais
lattice, and the expected memory needs. It echoes the input data, and report
on checks of data consistency for each dataset.

<a id="6.3"></a>
### 6.3 More on the main output file
  
Then, for each dataset, the real computation is done, and the code will report
on some initialisations, the SCF convergence, and the final analysis of
results for this dataset. Each of these phases is now described in more
details.

The code reports:

  * the real and reciprocal space translation vectors ( _Note_: the definition of the reciprocal vector is such that Ri.Gj= deltaij)
  * the volume of the unit cell
  * the ratio between linear dimension of the FFT box and the sphere of plane waves, called boxcut

It must be above 2 for exact treatment of convolutions by FFT. 
[[ngfft]] has been automatically chosen to give a boxcut value larger than 2, but not much
larger, since more CPU time is needed for larger FFT grids;

  * the code also mention that for the same FFT grid you might treat (slightly) larger [[ecut]] 
    (so, with a rather small increase of CPU time) 

  * the heading for each pseudopotential which has been input 

  * from the inwffil subroutine, a description of the wavefunction initialization 
  (random number initialization or input from a disk file), that is, a report 
  of the number of planewaves (npw) in the basis at each k point

  * from the setup2 subroutine, the average number of planewaves over all k points is reported in two forms, 
  arithmetic average and geometric average. 

Until here, the output of a ground-state computation is identical to the one
of a response-function calculation. See the [[help:respfn]] for the latter,
especially [[help:respfn#6.2|section 6.2]].

Next the code reports information for each SCF iteration:

  * the iteration number 
  * the (pseudo) total energy (Etot) in Hartree [This is not the total energy of the system, 
    since the pseudopotential approximation has been made: a constant energy (in the frozen-core approximation) 
    should be added to the present pseudo total energy in order to obtain a total energy, 
    that includes the contributions from the core electrons. Since only differences of total energy matter 
    (except is extremely rare cases), one can work with this pseudo energy like if it were the true total energy, 
    except that the missing constant depends on the pseudopotential that has been used. 
    Thus one has to perform differences of pseudo energies between simulations that use the same pseudopotentials]. 
  * the change in Etot since last iteration (deltaE) 
  * the maximum squared residual residm over all bands and k points 
    (residm - the residual measures the quality of the wavefunction convergence) 
  * the squared residual of the potential in the SCF procedure (vres2) 
  * the maximum change in the gradients of Etot with respect to fractional coordinates (diffor, in Hartree) 
  * the rms value of the gradients of Etot with respect to fractional coordinates (maxfor, in Hartree).  
    The latter two are directly related to forces on each atom.

  * Then comes an assessment of the SCF convergence: the criterion for fulfillment of the SCF criterion 
    (defined by [[toldfe]], [[toldff]], [[tolwfr]] or [[tolvrs]]) might be satisfied or not... 
  * Then the stresses are reported. 

This ends the content of a fixed atomic position calculation.

Many such blocks can follow.

When the atomic positions have been eventually relaxed, according to the value
of [[ntime]], the code output more information:

  * The squared residuals for each band are reported, k point by k point. 
  * Then the fractional or reduced coordinates are given, 
  * followed by the energy gradients, 
  * followed by the cartesian coordinates in Angstroms, 
  * followed by the cartesian forces in Hartree/Bohr and eV/Angstrom. 
  * Also are given the rms force ( **frms** ) and the maximum absolute value of any force component ( **max** ). 
  * Next are the length scales of the unit cell in Bohr and in Angstroms. 
  * Next are the eigenvalues of each band for each k point, in eV or Hartree or both depending on the choice of [[enunit]].   

<a id="averagepot"></a>
NOTE that the average electrostatic potential of a periodically repeated cell is UNDEFINED.  
In the present implementation, the average Hartree potential and local
potential are imposed to be zero, but not the average exchange-correlation
potential. This definition gives some meaning to the absolute values of
eigenenergies, thanks to Janak's theorem: they are derivatives of the total
energy with respect to occupation number. Indeed, the G=0 contributions of the
Hartree, local potential and ion-ion to the total energy is independent of the
occupation number in the present implementation. With this noticeable
exception, one should always work with **differences** in eigenenergies, as
well as **differences** between eigenenergies and the potential. For example,
the absolute eigenenergies of a bulk cell should not be used to try to predict
a work function. The latter quantity should be obtained in a supercell
geometry, by comparing the Fermi energy in a slab and the potential in the
vacuum in the same supercell.

  * Next are the minimum and maximum values for charge density, and next smaller or larger values (in order to see degeneracies). 

  * Next are the total energy (Ha and eV) and its components: 
    * kinetic, 
    * Hartree, 
    * exchange and correlation (xc), 
    * Ewald (ion-ion energy), 
    * *core correction* to the local pseudopotential, 
    * local pseudopotential, and 
    * nonlocal pseudopotential.  The sum of the Kohn-Sham energies (termed "band energy") is also given. 

  * Next is the stress tensor, (1/ucvol) d(Etot)/d(strain(a,b))
    for Etot=total energy per unit cell and **a**, **b** are **x**, **y**, or **z** components of strain.
    The stress tensor is given in cartesian coordinates in Hartree/Bohr 3 and GPa.
    The basics of the stress tensor are described in [[cite:Nielsen1985]] and [[cite:Nielsen1985a]]. 

Having finished all the calculations for the different datasets, the code
echoes the parameters listed in the input file, using the latest values e.g.
for [[xred]], [[vel]], and [[xcart]], and supplement them with the values
obtained for the total energy, the forces and stresses, as well as occupation numbers.
The latter echoes are very convenient for a quick look at the result of calculation!

This is followed finally by the timing output: both "cpu" time and "wall
clock" time as provided by calls within the code.
The total cpu and wall clock times are reported first, in seconds, minutes,
and hours for convenient checking at a glance.
Next are the cpu and wall times for the principal time-consuming subroutine
calls, each of which is independent of the others. The sum of these times
usually accounts for about 90% of the run time.

The main subroutines, for BIG jobs, are

1. fourwf: the subroutine which performs the fast Fourier transform for the wavefunctions: 
2. fourdp: the subroutine which performs the fast Fourier transform related to density and potential 
3. rhohxc: computes the Hartree and exchange-correlation energy and potential and sometimes 
   derivative of potential; only the XC timing is reported, excluding time connected to the FFTs: `xc:pot/=fourdp.`
4. nonlop: computes the matrix elements of the nonlocal pseudopotential: $\langle G|V_{non-local}|C \rangle$
5. projbd: Gram-Schmidt orthogonalisation 

In case of small jobs, other (initialisation) routines may take a larger
share, and the sum of the times for the principal time-consuming subroutine
calls will not make 90% of the run time.

If the long printing option has been selected ([[prtvol]]=1), the code gives
much more information in the whole output file. These should be rather 
self-explanatory, usually. Some need more explanation.  
In particular the cpu and wall times for major subroutines which are NOT
independent of each other; for example vtorho conducts the loop over k points
and calls practically everything else. In case of a ground state calculation,
at fixed atomic positions, these subroutines are:

1. **abinit**: the main routine 
2. **driver**: select ground state or response calculations 
3. **gstate**: the driver of the ground state calculations 
4. **scfcv_core**: the SCF cycle driver 
5. **vtorho**: compute the density from the potential (it includes a loop over spins and k-points) 
6. **vtowfk**: compute the wavefunctions at one particular k-point (includes a non self consistent loop, and a loop over bands) 
7. **cgwf**: optimize one wavefunction in a fixed potential 
8. **getghc**: computes $\langle G|H|C\rangle$, that is, applies the Hamiltonian operator to an input vector. 

<a id="header"></a>
### 6.4 The header
  
The **wavefunction files**, **density files**, and **potential files** all
begin with the same records, called the "header".
This header is treated using the hdr_type Fortran data structure inside ABINIT. 
There are dedicated routines inside ABINIT for initializing a header, updating it,
reading the header of an unformatted disk file, writing a header to an unformatted disk file,
echoing a header to a formatted disk file, cleaning a header data structure.

The header is made of 4 + [[ntypat]] unformatted records, obtained by the
following Fortran90 instructions (format 5.7):
    
```fortran
     write(unit=header) codvsn,headform,fform
     write(unit=header) bantot,date,intxc,ixc,natom,ngfft(1:3),&
    & nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase,usepaw,&
    & ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3),rprimd(1:3,1:3),stmbias,tphysel,tsmear,usewvl

     write(unit=header) istwfk(1:nkpt),nband(1:nkpt*nsppol),&
    & npwarr(1:nkpt),so_psp(1:npsp),symafm(1:nsym),symrel(1:3,1:3,1:nsym),typat(1:natom),&
    & kpt(1:3,1:nkpt),occ(1:bantot),tnons(1:3,1:nsym),znucltypat(1:ntypat),wtk(1:nkpt)
     do ipsp=1,npsp
    ! (npsp lines, 1 for each pseudopotential; npsp=ntypat, except if alchemical pseudo-atoms)
      write(unit=unit) title,znuclpsp,zionpsp,pspso,pspdat,pspcod,pspxc,lmn_size
     enddo
    !(in case of usepaw==0, final record: residm, coordinates, total energy, Fermi energy)
     write(unit=unit) residm,xred(1:3,1:natom),etotal,fermie
    !(in case of usepaw==1, there are some additional records)
     if (usepaw==1)then
      write(unit=unit)( pawrhoij(iatom)%nrhoijsel(1:nspden),iatom=1,natom), cplex, nspden
      write(unit=unit)((pawrhoij(iatom)%rhoijselect(1:      nrhoijsel(ispden),ispden),ispden=1,nspden),iatom=1,natom),&
    &                 ((pawrhoij(iatom)%rhoijp     (1:cplex*nrhoijsel(ispden),ispden),ispden=1,nspden),iatom=1,natom)
     endif
```

where the type of the different variables is:
    
```fortran
    character*6 :: codvsn
    integer :: headform,fform
    integer :: bantot,date,intxc,ixc,natom,ngfft(3),nkpt,npsp,
     nspden,nspinor,nsppol,nsym,ntypat,occopt,pertcase,usepaw
     integer :: usewvl, cplex, nspden
    double precision :: acell(3),ecut,ecutdg,ecutsm,ecut_eff,qptn(3),rprimd(3,3),stmbias,tphysel,tsmear
    integer :: istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt),so_psp(npsp),&
    & symafm(nsym),symrel(3,3,nsym),typat(natom),nrhoijsel(nspden),rhoijselect(*,nspden)
    double precision :: kpt(3,nkpt),occ(bantot),tnons(3,nsym),znucltypat(ntypat),wtk(nkpt)
    character*132 :: title
    double precision :: znuclpsp,zionpsp
    integer :: pspso,pspdat,pspcod,pspxc,lmax,lloc,mmax=integers
    double precision :: residm,xred(3,natom),etotal,fermie,rhoij(*,nspden)
```

NOTE: _etotal is set to its true value only for density and potential files.
For other files, it is set to 1.0d20_  
NOTE: _ecut_eff= [[ecut]]* [[dilatmx]] 2_  
NOTE: _For all cases where occupation numbers are defined (that is, positive
iscf, and iscf=-3), and for non-metallic occupation numbers, the Fermi energy
is set to the highest occupied eigenenergy. This might not correspond to the
expected Fermi energy for a later non-self-consistent calculation (e.g. the band structure)_

The header might differ for different versions of ABINIT. One pre-v5.3 format
is described below. Note however, that the current version of ABINIT should be
able to read all the previous formats (not to write them), with the exception
of wavefunction files for which the [[ecutsm]] value was non-zero (there has
been a change of definition of the smearing function in v4.4).

The format for version 4.4, 4.5, 4.6, 5.0, 5.1 and 5.2 was:
   
```fortran
     write(unit=header) codvsn,headform,fform
     write(unit=header) bantot,date,intxc,ixc,natom,ngfft(1:3),&
    & nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase,usepaw,&
    & ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3),rprimd(1:3,1:3),stmbias,tphysel,tsmear
     write(unit=header) istwfk(1:nkpt),nband(1:nkpt*nsppol),&
    & npwarr(1:nkpt),so_typat(1:ntypat),symafm(1:nsym),symrel(1:3,1:3,1:nsym),typat(1:natom),&
    & kpt(1:3,1:nkpt),occ(1:bantot),tnons(1:3,1:nsym),znucltypat(1:ntypat)
     do ipsp=1,npsp
    ! (npsp lines, 1 for each pseudopotential; npsp=ntypat, except if alchemical pseudo-atoms)
      write(unit=unit) title,znuclpsp,zionpsp,pspso,pspdat,pspcod,pspxc,lmn_size
     enddo
    !(in case of usepaw==0, final record: residm, coordinates, total energy, Fermi energy)
     write(unit=unit) residm,xred(1:3,1:natom),etotal,fermie
    !(in case of usepaw==1, there are some additional records)
     if (usepaw==1)then
      write(unit=unit)(pawrhoij(iatom)%nrhoijsel(1:nspden),iatom=1,natom)
      write(unit=unit)((pawrhoij(iatom)%rhoijselect(1:nrhoijsel(ispden),ispden),ispden=1,nspden),iatom=1,natom),&
    &                 ((pawrhoij(iatom)%rhoijp     (1:nrhoijsel(ispden),ispden),ispden=1,nspden),iatom=1,natom)
     endif
```

<a id="denfile"></a>
### 6.5 The density output file
  
This is an unformatted data file containing the electron density on the real
space FFT grid. It consists of the header records followed by

```fortran
    do ispden=1,nspden
     write(unit) (rhor(ir),ir=1,cplex*ngfft(1)*ngfft(2)*ngfft(3))
    enddo
```

where **rhor** is the electron density in electrons/Bohr^3, and cplex is the
number of complex components of the density (cplex=1 for GS calculations -the
density is real-, and cplex=1 or 2 for RF). The input variable [[nspden]]
describes the number of components of the density. The first component (the
only one present when [[nspden]]=1) is always the total charge density. When
[[nspden]]=2, the second component is the density associated with spin-up
electrons. When [[nspden]]=4, the second, third and fourth components
correspond to the x, y and z projections of the local magnetization, in units
of hbar/2. Note that the meaning of the different components of the density
differs for the density array (rhor) and for the different potential arrays
(vxc...), see section  6.6.

To identify the points in real space which correspond with the index "ir"
above, consider the following.  
The first array value (ir=1) corresponds with the first grid point which is at
the origin of the unit cell, (x=0, y=0, z=0).  
The next grid point (ir=2) lies along the first primitive translation at the
next fft grid point, which is (1/[[ngfft]](1))*[[acell]](1)*[[rprim]](mu,1).
This is 1/[[ngfft]](1) of the way along the first primitive translation.  
The rest of the values up to ir=[[ngfft]](1) lie along this vector, at
(ir-1)/[[ngfft]](1) of the way along the first primitive translation. The
point at ir=[[ngfft]](1)+1 lies at 1/[[ngfft]](2) along the second primitive translation.  
The next points up to ir=[[ngfft]](1)+[[ngfft]](1) are displaced in the
direction of the second primitive translation by 1/[[ngfft]](2) and in the
first translation by (ir-[[ngfft]](1)-1)/[[ngfft]](1).  
This pattern continues until ir=[[ngfft]](1)*[[ngfft]](2).  
The next point after that is displaced along the third primitive translation
by 1/ngfft(3), and so forth until ir varies all the way from 1 to
[[ngfft]](1)*[[ngfft]](2)*[[ngfft]](3). This last point is in the corner
diagonally opposite from the origin, or right alongside the origin if the
whole grid is viewed as being periodically repeated.

<a id="localpotfile"></a>
### 6.6 The potential files
  
Also unformatted files consisting of the header records and
    
```fortran
    do ispden=1,nspden
     write(unit) (potential(ir),ir=1,cplex*ngfft(1)*ngfft(2)*ngfft(3))
    enddo
```

where **potential** can be either the sum of the Hartree potential, exchange-
correlation and local pseudopotential (see [[prtpot]]), the Hartree potential
(see [[prtvha]]), the Hartree+XC potential (see [[prtvhxc]]), the local
pseudopotential (see [[prtvpsp]]) or the XC potential (see [[prtvxc]]), These
are defined on the real space grid in Hartree energy units. The underlying
grid is as described above. If [[nspden]]=2, the different components are the
spin-up potential and the spin-down potential. In the case [[nspden]]=4, the
components correspond to the up-up potential, the down-down potential, the
real part of the up-down potential, and the imaginary part of the up-down
potential. Note that the Hartree potential is NOT spin-dependent, but in order
to use the same format as for the other potential files, the spin-independent
array is written twice, once for spin-up and one for spin-down.

<a id="wfkfile"></a>
**6.7. The wavefunction output file**

This is an unformatted data file containing the planewaves coefficients of all
the wavefunctions, and different supplementary data.

The **ground-state** wf file consists of the header records, and data written
with the following lines of FORTRAN (version 4.0 and more recent versions):

```fortran    
          bantot=0                                    <-- counts over all bands
          index=0                                     <-- index for the wavefunction location
          do isppol=1,nsppol
           do ikpt=1,nkpt
            write(unit) npw,nspinor,nband                    <-- for each k point
            write(unit) kg(1:3,1:npw)                        <-- plane wave reduced coordinates
            write(unit) eigen(1+bantot:nband+bantot),        <-- eigenvalues for this k point
                        occ(1+bantot:nband+bantot)           <-- occupation numbers for this k point
            do iband=1,nband
             write(unit) (cg(ii+index),ii=1,2*npw*nspinor)   <-- wavefunction coefficients
            enddo                                            for a single band and k point
            bantot=bantot+nband
            index=index+2*npw*nspinor*nband
           enddo
          enddo
```

If the job ended without problem, a few supplementary lines are added, in
order to give the history of atomic positions and corresponding forces. The
integer nxfh gives the number of pairs (x,f) of positions and forces in reduced coordinates:
    
```fortran
     write(unit)nxfh
     do ixfh=1,nxfh
      write(unit) xred(1:3,1:natom,ixfh),dummy(1:3,1:4),&
    &             fred(1:3,1:natom,ixfh),dummy(1:3,1:4)
     enddo
```

The dummy variables might contain, in the future, the description of the unit
cell, and the stresses. The type of the different variables is:
    
```fortran
    integer :: kg,nband,npw,nspinor,nxfh
    double precision :: cg,dummy,eigen,fred,occ,xred
```

The **response-function** wf file consists of the header records, and data
written with the following lines of FORTRAN (version 4.0 and more recent versions):
    
```fortran
    bantot=0                                    <-- counts over all bands
    do isppol=1,nsppol
     do ikpt=1,nkpt
      write(unit) npw,nspinor,nband                    <-- for each k point
      write(unit) kg(1:3,1:npw)                        <-- plane wave reduced coordinates
      do iband=1,nband
       write(unit) (eigen(jband+(iband-1)*nband+bantot),jband=1,2*nband)  <-- column of eigenvalue matrix
       write(unit) (cg(ii+index),ii=1,2*npw*nspinor)     <-- wavefunction coefficients
      enddo                                            for a single band and k point
      bantot=bantot+nband
     enddo
    enddo
```

In version previous to 4.0, npw and nspinor were combined:
    
```fortran
    write(unit) npw*nspinor,nband
```

while the planewave coordinate record was not present (in both GS and RF cases).

Note that there is an alternative format (_KSS) for the output of the
wavefunction coefficients, activated by a non-zero value of [[nbandkss]].

### 6.8 Other output files

  
There are many other output files, optionally written, all formatted files at
present. Their use is usually governed by a specific input variable. Please
consult the description of this input variable, in order to have more
information on such files:

  * [[prtdos]] to print a file with the electronic Density-Of-States
  * [[prteig]] to print a file with the list of k points and eigenenergies
  * [[prtgeo]] to print a file with a geometrical analysis (bond lengths and bond angles), that also contains an XMOL section
  * [[prt1dm]] to print a one-dimensional projection of potential and density, for the three axes.

### 6.9 Control of output in the parallel case
  
For massively parallel runs, one cannot afford to have some of the output
files that are usually created. Explicitly, the log file and also the status
file become problematic. By default, with less than N processors, they are
created, but beyond N processors, they are deactivated except for the main log
file (master processor).

This default behaviour can be changed as follows. If a file named "_NOLOG"
exists in the current directory, then no log file and no status file will be
created, even with less than N processors. By contrast, if a file "_LOG"
exists in the current directory, then a log file and the status files will be
created, even with more than N processors. Alternatively, if a file named
"_MAINLOG" exists and there are less than N processors, only the master
processor writes the log and status files (this mimic the default behavior
when using more than N processors but with less than N processors)

In ABINITv7, N was set at N=100. However, with ABINITv8, N has been switched
to 2. It can be changed "by hand", though: modify NPROC_NO_EXTRA_LOG in
src/10_defs/defs_basis.F90 and recompile. See src/95_drive/iofn1.F90 for more explanation.

<a id="7"></a>
## 7 Numerical quality of the calculations
  
The following section describes various parameters which affect convergence
and the numerical quality of calculations.

The list of these input parameters is

  * (1) [[ecut]] 
  * (2) [[toldfe]], [[toldff]], [[tolwfr]], and [[tolvrs]], as well as [[nstep]] 
  * (3) [[nkpt]] 
  * (4) [[ngfft]] 
  * (5) [[tolmxf]], as well as [[amu]], [[dtion]], [[vis]], [[ntime]] 
  * (6) [[acell]] and [[rprim]] 
  
The technical design of the pseudopotential also affects the quality of the results.

(1) The first issue regarding convergence is the number of planewaves in the
basis for a given set of atoms. Some atoms (notably those in the first row or
first transition series row) have relatively deep pseudopotentials which
require many planewaves for convergence. In contrast are atoms like Si for
which fewer planewaves are needed. A typical value of [[ecut]] for silicon
might be 5-10 Hartree for quite good convergence, while the value for oxygen
might be 25-35 Hartree or more depending on the convergence desired and the
design of the pseudo- potential.

NOTE: It is necessary in every new problem to **TEST** the convergence by
**RAISING** [[ecut]] for a given calculation until the results being computed
are constant to within some tolerance. This is up to the user and is very
important. For a given [[acell]] and [[rprim]], [[ecut]] is the parameter
which controls the number of planewaves. Of course if [[rprim]] or [[acell]]
is varied then the number of planewaves will also change.

Let us reiterate that extremely careful pseudopotential design can optimize
the convergence of _e.g._ the total energy within some range of planewave
number or [[ecut]]. It is appropriate to attempt to optimize this convergence,
especially for difficult atoms like oxygen or copper, as long as one does not
significantly compromise the quality or transferability of the
pseudopotential. There are many people working on new techniques for
optimizing convergence.

For information on extended norm conservation, see E. L. Shirley, D. C. Allan,
R. M. Martin, and J. D. Joannopoulos, Phys. Rev. B 40, 3652 (1989).

For information on optimizing the convergence of pseudopotentials, see A. M.
Rappe, K. M. Rabe, E. Kaxiras, and J. D. Joannopoulos, Phys. Rev. B 41, 1227 (1990).

(2) In addition to achieving convergence in the number of planewaves in the
basis, one must ensure that the SCF iterations which solve the electronic
structure for a given set of atomic coordinates are also converged. This
convergence is controlled by the parameters [[toldfe]], [[toldff]],
[[tolwfr]], and [[tolvrs]], as well as the parameter [[nstep]]. One of the
"tolerance" parameters must be chosen, and, when the required level of
tolerance is fulfilled, the SCF cycles will stop. The [[nstep]] variable also
controls convergence in preconditioned conjugate gradient iterations by
forcing the calculation to stop whenever the number of such iterations exceeds
nstep. Usually one wants nstep to be set larger than needed to reach a given
tolerance, or else one wants to restart insufficiently converged calculations
until the required tolerance is reached.

Note that, if the gap in the system closes (e.g. due to defect formation or if
the system is metallic in the first place), the presently coded algorithm will
be slower to converge than for insulating materials. Convergence trouble
during iterations usually signals closure of the gap. The code will suggest to
treat at least one unoccupied state (or band) in order to be able to monitor such a closure.

(3) For self consistent calculations ([[iscf]] positive) it is important to
test the adequacy of the k point integration. If symmetry is used then one
usually tests a set of "special point" grids. Otherwise one tests the addition
of more and more k points, presumably on uniform grids, to ensure that a
sufficient number has been included for good k point integration. The
parameter nkpt indicates how many k points are being used, and their
coordinates are given by kpt and kptnrm, described above. The weight given to
each k point is provided by input variable [[wtk]]. Systematic tests of k
point integration are much more difficult than tests of the adequacy of the
number of planewaves. The difficulty I refer to is simply the lack of a very
systematic method for generating k point grids for tests.

(4) It is possible to run calculations for which the fft box is not quite
large enough to avoid aliasing error in fft convolutions. An aliasing error,
or a Fourier filter approximation, is occurring when the output variable "
**boxcut** " is less than 2. boxcut is the smallest ratio of the fft box side
to the planewave basis sphere diameter. If this ratio is 2 or larger then e.g.
the calculation of the Hartree potential from the charge density is done
without approximation.  

NOTE: the values of [[ngfft]](1:3) are chosen automatically by the code to
give boxcut > 2, if [[ngfft]] has not been set by hand. At ratios smaller than
2, certain of the highest Fourier components are corrupted in the convolution.
If the basis is nearly complete, this Fourier filter can be an excellent
approximation. In this case values of boxcut can be as small as about 1.5
without incurring significant error. For a given [[ecut]], [[acell]], and
[[rprim]], one should run tests for which [[ngfft]] is large enough to give
boxcut >= 2, and then one may try smaller values of [[ngfft]] if the results
are not significantly altered. See the descriptions of these variables above.

(5) If you are running calculations to relax or equilibrate structures, i.e.
with [[ionmov]]=1 and possibly [[vis]]>0, then the quality of your molecular
dynamics or relaxation will be affected by the parameters [[amu]], [[dtion]],
[[vis]], [[ntime]], [[tolmxf]]. Clearly if you want a relaxed structure you
must either run long enough or make repeated runs until the largest force in
the problem (output as fmax) is smaller than what you will tolerate (see
[[tolmxf]]).  
If [[dtion]] is too large for the given values of masses ([[amu]]) and
viscosity ([[vis]]) then the molecular dynamics will be unstable. If [[dtion]]
is too small, then the molecular dynamics will move inefficiently slowly. A
consensus exists in the community that forces larger than about 0.1
eV/Angstrom are really too large to consider the relaxation to be converged.
It is best for the user to get experience with this in his/her own
application.  
The option [[ionmov]]=2, 3 or 7 are also available This uses the Broyden
(BFGS) scheme for structural optimization and is much more efficient than
viscous damping for structural relaxation.

(6) If you are running supercell calculations (i.e. an isolated atom or
molecule in a big box, or a defect in a solid, or a slab calculation) you must
check the convergence of your calculation with respect to the supercell and
system size.

  * For an isolated molecule in a big box: increase concurrently the three dimensions 
    of your supercell ([[acell]]), and check the convergence of your physical property. 
  * For a defect in a solid: your supercell must be a multiple of the primitive cell of the bulk solid, 
    so you have less freedom. Still, be sure that your supercell is large enough for your properties 
    of interest to be accurate at the level you want it to be. 
  * For a slab calculation: you must increase the vacuum in the cell, but also the thickness of your slab systematically... 
