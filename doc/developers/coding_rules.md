## ABINIT style for Fortran programming
(revised many times from the original draft in 1991, BUT NOT REALLY UP_TO_DATE ... SHOULD BE REVISED ONCE MORE)

### 1. Foreword 

The ABINIT code should conform to the ANSI Fortran norm (Fortran 2003).
Still, most of the following rules can already used with Fortran95. This norm is abbreviated F95, while the older norm, Fortran 77, is abbreviated F77.
The following set of rules complements the F95 (or Fortran2003) standard.
The ABINIT code conforms to most of these additional rules already.
Each modification, or new routine is expected to adopt all
the standards, that will be also enforced gradually in existing routines.
The developer is supposed to have a good knowledge of the Fortran constructs. 
For an introduction to them, consult the following URL: (to be filled)

The rules described afterwards do NOT apply to the routines that are in the *_noabirules directories.
Their source is voluntarily not modified,
to keep coherence with the original versions, provided by somebody who followed other coding rules.
Names having a special meaning in F95 are written here in uppercase.
They are usually NOT in uppercase in the code, however.

  * Could be nice to enforce that ?

There must be basic rules that are enforced: lowercase for routine names, variable names in English, comments in English, explicative names! Avoid confusing names.

The basic utility of this document is to provide a common style to all the routines. 
This is important,in order for ABINIT developers to understand quickly what a routine does: 
they expect some information to be always present in some part of each routine, or, in the
opposite sense, they should know where to place each information (and not spend time to evaluate
where each information should be placed).
This document is also a way to enforce a kind of "locality of information" principle: all the information needed to understand what a routine does should be contained in the same routine, and if possible, the information about a small part of the routine should be contained within a screen page of this part of the routine.

In general, the developer will conform to the style by taking an existing routine as an example, and
use it to generate a new routine. 
However, it is worth reading the present rules from time to time ...
The new rules will be enforced gradually.
This means that not only the format, but also the content of present document is IN A TRANSITIONAL state.

### 2. Declarations

  * Use ''IMPLICIT NONE'', forcing explicit declaration of all variables, in all subroutines or modules

  * ''COMPLEX*16'' or ''REAL*8'' are not ANSI standards. In the ''defs_basis'' module, one defines the integer variable dp <code>INTEGER, PARAMETER :: dp = kind(1.0d0)</code> Accordingly, define floating point numbers as<code>

          REAL(dp) :: real_variable 
          !(the double precision is standard, avoid single precision). 
          For the declaration of complex numbers, use two real numbers:
      
          REAL(dp) :: complex_variable(2) 
          ! That is, place complex numbers in a single array with 
          ! real parts in odd array components 
          ! and imaginary parts in even array components.</code>
       
  * When setting values of constants, include ''_dp'' on the end to indicate double precision constants, e.g. 0.0_dp. For integers, one has the choice between:<code>
   <code>     
      INTEGER :: ii
      
      ! and          
      
      INTEGER(i4b) :: ii
      
      ! where i4b has been defined in the 'defs_basis' module</code>
   </code>
  * Mention all the attributes of a variable in a single-line declaration. This means that ''PARAMETER'', ''SAVE'', ''ALLOCATABLE'', ''INTENT(IN)'', etc. should not be declared later. However, avoid the ''DIMENSION'' attribute, and declare the dimension after each array. This is to save numerous line in the declaration section. The use of the double colon is the standard. 
  (JB:For arrays, I feel like one per line is enough since it becomes harder to have a glance at the dimensions afterwards) For example<code>

        REAL(dp), PARAMETER :: aa=1.0d0, bb=1.0d0/3.0d0
        
        INTEGER, ALLOCATABLE :: cc(:)</code>

  * For character strings, use:

        CHARACTER(len=10) :: name(5)

  * For a string of characters whose length is known only in the calling subroutine, use<code>
  
        CHARACTER(len=*) :: names(8)</code>

  * Logical variables are allowed, especially for very big arrays, where the memory gain is important, or for temporary logical arrays (work arrays). However, when a variable has to describe a choice, use an integer instead of a logical: this paves the way to the introduction of further choices. (JB: discussion about magic numbers ?)

  * This should not be enforced but we should provide guidelines: if we suspect that the size of the array may be very large, use allocatable arrays.... Avoid automatic arrays (except for small arrays such as ''local_arr(small_num)''), use allocatable arrays.
 
     - The handling of automatic and allocatable variables in memory is different: automatic variables are placed in the stack, while allocatable arrays are placed in the heap. The stack size is usually limited.
     - This also allows the checking of the parts of the program where memory is allocated, through the use of a grep allocatable * or  grep allocate * command. In order for such a search to be efficient, never continue an allocate or allocatable command on the next line. (see also rule 5.f)

  * Should be used only if the "pointer" power is needed to avoid waste of memory, never use pointer for dynamic (temporary) arrays)). Use pointer to point something. Coding rule: nullify pointer when you are done with init, set pointers to null() in the datypes. Use allocatable Example:<code>

        REAL(dp),DIMENSION(:),POINTER :: array
        nullify(array)
        CALL init(array)
        (using array)
        DEALLOCATE(array)

        ! where the subroutine init contains:

        INTEGER :: size_array
        REAL(dp),DIMENSION(:),POINTER :: array
        size_array=...(some expression)
        ALLOCATE(array(size_array))  
        (initialize array) 
        END SUBROUTINE</code>

  * Assumed shape arrays are encouraged when the shape of the array is not known and/or efficiency is not important. Avoid passing the arrays descriptot when you are calling F77 routines e.g. BLAS/LAPACK)
assumed shape arrays Continue to use explicit shape  the use of assumed shape array is allowed Always declare explicitly in the subroutine the dimension of a array. This is to allow to understand locally exactly what array operations are doing. Try to  High-level routines with many arguments should use datatypes to gather the information so that we can simplify the API. We cannot continue to use optional arguments and/or zillions of options and arrays to handle all the possible cases. Look for example at the XC routines and the implementation of MGGA!JB: YES ! Actually the number of arguments should be limited to a defined value like let say 10 and be enforced !)  So, use:<code>

       INTEGER :: zat(natom)
       
       ! and avoid
       INTEGER :: zat(*) ! to be avoided</code>

  * JB : Sometimes it is usefull to not declare explicitely the size of an array but just it shape. For instance: Remember that an explicit interface is required (declare your procedure in a module if possible)<code>
      REAL :: gw(:,:) ! freq, flavors
      INTEGER :: nfreq, nflavors

      nfreq=SIZE(gw,1)
      nflavor=SIZE(gw,2)</code>
      This allows to work on a given subset of an array and the routine is then adaptative
    

  
### 3. Variables

  * Never use a one-letter name for a variable, because they are too difficult to find systematically. Even two-letter variables (except duplicated letters, easy to find, or names including j, k, q, x, or z, or some special variables like "dp") should be avoided. So, use ii,jj,kk instead of i,j,k. Note that F90 allows for names longer than 6 characters, unlike F77. (enforced)

  * Naming of looping variables (integers):

          Use i... for the current (dummy) value of a given index
          Use n... for the number of a given quantity, in general
                   also equal to its maximum (unlike in the old standard)
          Use m... for the maximum number of a given quantity,
                   when it is not equal to the actual number of this given quantity
          For example: iatom, natom, matom, respectively (partially enforced).
          
  * Try to name your variables following the standard names used in Abinit (e.g. cg) 

  * For global variables, use explicative names, possibly with a prefix that is constructed from the name of the module (e.g. XMPI_WORLD). Avoid global variables that can be changed at run-time or provide a setter/getter method 

### 4. File format

  * Use the Fortran90 free format. (132 column maximum, continuation sign & ... &, comment sign !, no  restriction on the use of the first columns). Note that doubling the comment sign (!!) will be used to allow processing by Robodoc.

  * Always use the continuation sign on both the continued line (at its end), and on the continuation line (at its beginning) (enforced)

  * Use indentation of two columns for each DO-loop or each IF-block, or any other of structure (usually enforced)

  *  Don't break the indentation even when you insert CPP sections 


The header of each routine must use a special format, so that it can be processed by the 
Robodoc documentation program (see http://www.xs4all.nl/~rfsber/Robo/robodoc.html and http://www.xs4all.nl/~rfsber/Robo/manual.html).

This header of the file will have to mention:

    - the name of the subroutine
    - the function of the subroutine
    - the copyright
    - the complete set of input arguments with explanation
    - the complete set of output arguments with explanation
    - the complete set of arguments that are both input and output, with explanation
    - the list of calling routines (produced automatically: do not worry about it).
    - eventually some notes.

Then, the body of the routine will follow (with CPP F90 instructions, see rule 3.f).
For Robodoc treating Fortran 90, there is a special meaning attached
to the lines of the source files that begin with a double bang !!
There are also keywords (''NAME'', ''FUNCTION'', ''COPYRIGHT'', ''INPUTS'', ''OUTPUT'', ''SIDE EFFECTS'', ''NOTES'', ''SOURCE'' ...) to be used.

(MG: We should standardize the format used to write the documentation inside the fields. RST is an interesting possibility...)

The header format is the most easily understood by a commented example:

<code>
  !!****f* ABINIT/name_of_subroutine
  !! NAME
  !! name_of_subroutine
  !!
  !! FUNCTION
  !! explanation of "name_of_subroutine"
  !! will be described here ...
  !!
  !! COPYRIGHT
  !! Copyright (C) yyy1[-yyy2] ABINIT group (initials)
  !! This file is distributed under the terms of the
  !! GNU General Public License, see ~abinit/COPYING
  !! or http://www.gnu.org/copyleft/gpl.txt .
  !! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
  !!
  !! INPUTS
  !!  arg1 = ...
  !!  arg2 = ...
  !!
  !! OUTPUT
  !!  arg3 = ...
  !!  arg4 = ...
  !!
  !! SIDE EFFECTS
  !! Input/Output   ! if this line is found in a routine, use "SIDE EFFECT".
  !!  arg5 = ...
  !!
  !! PARENTS        ! This paragraph will be produced automatically
  !!                ! It will mention the routines that call the present one.
  !!
  !! NOTES
  !!  (miscellaneous note, warning, etc., if any)
  !!
  !! SOURCE
</code>

The string "name_of_subroutine" must be replaced by the name of the subroutine.
In the copyright notice, "yyy1" and "yyy2" must be replaced by years of modifications,
"initials" must be replaced by the initials of the contributors, for example

<code>
   !! Copyright (C) 2007-2013 ABINIT group (DCA, XG).
   "arg1=..." 
</code>

and similar strings must be replaced by the actual description
of input or output variables. They must appear in the same order in which
they appear as arguments to the subroutine. Note that the Robodoc keywords
are "INPUTS" (with a "S") and "OUTPUT" (without a "S").

If physical quantities are being passed as arguments, give the appropriate dimension for the quantity
(e.g. Hartrees or Bohr). No work space private to each routine should be an argument
of the routine, use ALLOCATABLE, ALLOCATE, and DEALLOCATE inside the subroutine.

After the line

<code>
  !! SOURCE
</code>

follows the CPP and F90 lines (enforced).

3.f. After the header, and before the F90 declaration of variables, one must mention, in order:

  * eventually, global CPP definitions
  * the subroutine name and arguments
  * eventually, "include" files

As an example:<code>


  #define FNLEN 132    /* defines the length of file name variables */

  subroutine name_of_subroutine (arg1,arg2,arg3,arg4,arg5)

  include 'file'
</code>

3.g. 
Every routine should then have the declaration IMPLICIT NONE.
This should then be followed by other declarations as follows.
*Very important note*: the order of declarations is automatically enforced by the ''abirules'' script. 
The developer can perfectly leave to the script the duty of ordering his new declarations,
at the time of the merge with the rest of the code
Note that within each list, the variable names should be in alphabetic order.
(JB: INTENT should be enforced for all arguments, shouldn't it ?)

A1. Arguments: scalars


  INTENT(IN)
  (1) integers
  (2) double precision
  (3) other variables in argument list
  INTENT(OUT)
  (4) integers
  (5) double precision
  (6) other variables in argument list
  INTENT(INOUT)
  (7) integers
  (8) double precision
  (9) other variables in argument list


A2. Arguments: arrays same as (1)-(9à) above

B1. Local variables (sequential) : scalars

  INTEGER
  (1) integer parameters
  (2) integer with save attribute
  (3) other integers
  DOUBLE PRECISION
  (4) double precision parameters
  (5) double precision with save attribute
  (6) other double precision
  (7) other scalar variables

B2. Local variables (sequential) : arrays

  INTEGER
  (1) integer arrays
  (2) allocatable integer arrays
  (3) pointers to integer arrays
  DOUBLE PRECISION
  (4) real arrays
  (5) allocatable real arrays
  (6) pointers to real arrays
  (7) other array variables

C1. Local variables (parallel case: MPI code section) : scalars

C2. Local variables (parallel case: MPI code section) : arrays

See the description of the sequential variables

3.h. Do coding in lower case. Lower case is suggested because it looks prettier. 
Use standard English grammar for comments (start sentences with uppercase). 
The keywords DEBUG and ENDDEBUG, in uppercase, should emphasize each part of the code, commented 
or uncommented, that is written for the purpose of debugging. (enforced)

3.i. Use comments frequently to explain manipulations and make them in English. (enforced)

3.l. Equations mentioned in comment lines should be written in Latex, for processing using Src2tex.
There are several ways to describe latex equations in the F90 comments
(these might be introduced by one "!" or two) :MOVE this section to a separate section


  1)  !! $ equation $
  1') !  $ equation $
  2)  !! $\displaystyle equation$
  2') !  $\displaystyle equation$
  3)  !!{{\ \begin{equation}
      !!  ... equation ...
      !!\end{equation} }}
  3') !{{\ \begin{equation}
      !  ... equation ...
      !\end{equation} }}
  4)  !!{{\ \begin{eqnarray}
      !!  ... equation ...
      !!\end{eqnarray} }}
  4') !{{\ \begin{eqnarray}
      !  ... equation ...
      !\end{eqnarray} }}

For example:

  !! epsatm=$(4 \pi) \int_0^\infty [r^2 (V(r)+Zv/r) dr]$ (hartree)

If a long fraction sentence is found, "!! $\displaystyle equation $" may be better for visualization.

In some cases, 

  !!\begin{equation} ...!!\end{equation}

or

  !!\begin{eqnarray} ...!!\end{eqnarray}

format might be better, especially for very long equation lines.

3.m. Additional information about Src2tex and Robodoc.

As mentioned in rule 3.e , the content of the PARENT paragraph will be created or updated automatically.
Also, the source files will be preprocessed automatically before applying src2tex:

   * the very first line will be added: !{\src2tex{texfont=tt}}
   * the very last line will also be added: %%!!***%%

(usually enforced)

3.o. ***MG_NEIN. Let's fix the scripts so that we can write readable code! I have a fortran parser written in python taken from f2py. Detecting functions is not a piece of cake but it's feasible!*** In order to ease the writing of the script that generates "PARENT" routines before 
the self-documentation call to Robodoc, it was needed to impose two additional rules:

  * subroutines should be used instead of functions
  * the instruction CALL should always be the first of an instruction line. 

In particular, the construct:

      IF (condition) CALL subroutine  !to be avoided

should be avoided, and replaced by:

      IF(condition) THEN
         CALL subroutine
      ENDIF

(enforced)


3.q. 
The INTENT must always be specified, except for pointers (which have cannot have intent in F90 (MG: F2003 allows it and the intent refers to the pointer association status)) and subprogram names given as formal arguments (rare).
(JB:This rule should be upper with the argument declaration of a routine ... so my previous remark should be read as "this rule 3.q come to late")
Nevertheless pointers should be placed in list of formal arguments according their 'quasi-intent'
in source code (read/modified/affected data).
On the other hand, a subprogram given as formal argument behaves always like 'intent(in)'. 
In that case, the formal argument consists of the explicit interface of the subprogram.
A formal argument belonging to a structures type must have an 'intent(in)' (respectively 'intent(out)')
if all fields of the structure are read-only (resp. write-only). 
In all other cases, specify an 'intent(inout)' and detail (with comments) the status for each field.

### 5. Constructs for flow control

5.a. Do not use numeral labels. These can be avoided through the use of the construct

   DO ...
      ...
   ENDDO

as well as the instruction (new to F90) CYCLE or EXIT, or the constructs based on IF, CASE, WHERE or FORALL.
Note: Numeral labels can also be avoided for FORMAT statements by placing the format descriptor in a character string, like in:

   WRITE (unit, '(i3)' ) nn

or:

   READ (text, '(i6, f6.2)' ) nn, xx

A character string can be explicitly declared, then initialized, if it is to be used more than once. (enforced)

5.b. Literal labels or comment lines must be used for each construct that is longer than one page of a screen (about 20 lines), both at the beginning and at the end of the construct. (enforced)

### 6. Use of arrays 

6.a. Distinguish explicitly a vector operation from its scalar counterpart, by specifying colons:

suppose <code>
  REAL(dp) :: aa(10), bb(10), cc(10)</code>
then use
<code>
  aa(:)=bb(:)+cc(:)</code>
and not
<code>
  aa=bb+cc ! to be avoided
</code>  

However, this rule can lead to problems with the stack size, see rules 5c-e. (enforced) (JB:I think the stack size has really increased and for a small routine that need temporary array of small size this is better that spending time for allocation/deallocation. )

6.b. In order to improve the readability of the code, use vector operations instead of explicit DO-loop constructs. Use also the intrinsic subroutines ''DOT_PRODUCT'', or ''MATMUL'', as well as
functions like ''SUM'' or ''MAXVAL''. However, avoid using segments as arguments. 
Use directly the array name: you do not need to distinguish this operation from an hypothetical scalar counterpart.
For reason of efficiency, very heavily used parts of codes could be coded differently.
Use also BLAS and LAPACK routines for higher-level constructs, although they
do not lead to as good a readability (see later). (partially enforced)

6.c DO NOT pass as an argument of a subroutine, a segment of an array. For example :

<code>   call routine(aa(1:naa))         ! to be avoided</code>

Indeed, some compilers copy all the values of the segment into
a temporary variable. (**MG: And they are wrong because the compile shall pass the array descriptor since the routine must have a F90 interface. Here I'm strongly in favor of a F90 approach, especially if the routine if not very CPU-critical**)
(JB:This is not a problem if 1) the so called routine use the size function to know the size of the array instead of explicit declaration of the size, 2) the segment is aligned in memory : call routine(aa(1:naa,5) is fine but call routine (aa(2,1:nbb)) will creat a temporary array that slows down the code. This implies to correctly allocate the arrays...)
Thus, this is memory and CPU time consuming.
Moreover, the zone of memory that is used for that purpose can be limited
in size, and thus the code is constrained not by the full physical
memory of the machine, but by the size of that zone (stack size on the DECs).
For character strings, this restriction is not too important,
since character strings are rather short in any case, so the rule
is not to be followed for character strings. The same argument
applies when the size of the segment is know to be
limited in all cases (like in nonlop.f, setsym.f, xctotl.f or pseudopotential routines for example).
(enforced, with the restrictions mentioned)

6.d Avoid the following construct, except when the segment size is sufficiently small:

<code>  aa(:,ii)=aa(:,jj)</code>

  (JB: Is there really a problem there ? I think I do that and I don't have trouble)
  (JB:NB: all the implicit loop are good for new cpus since then can vectorize the instruction. So loops should really be though carefully)
On the DEC machines, the compiler thinks it must first store the segment aa(:,jj) in the stack, then copy the stack in aa(:,ii).
This not only doubles the time for copying, but can lead to stack size problems. 
This likely happens each time the same array appears both at the left-hand side and right-hand side of an assignment expression.
However, there is no problem with <code> aa(:,ii)=bb(:,jj) (enforced)</code>

6.e Avoid the use of large segments in read or write operations (**MG: not easy to accomplish, especially for binary IO**)

<code>  read(6,*)aa(:,ii)</code>

Again, on the DEC machines, the stack is used, and there are potential problems with the declared stack size.
Use instead :

<code>  read(6,*)(aa(jj,ii),jj=jmin,jmax)</code>

A full array, not presented as a segment, also has no problem :

<code>  read(6,*)bigarray</code>

But avoid :

<code>  read(6,*)bigarray(:)</code>

(enforced)

6.f For allocation and deallocation of arrays use the macros ABI_ALLOCATE(array, (dim1, dim2))) and ABI_DEALLOCATE(array) defined in ~abinit/src/incs/abi_common.h. Use ABI_DATATYPE_ALLOCATE(foo, shape) and ABI_DATATYPE_DEALLOCATE(foo) if foo is a user-defined datatype

### 7. Coding practice

  * Run innermost loop fastest in looping. This is simply less CPU time consuming in Fortran.(JB:I'm sad it cannot be enforced :o()

  * Use generic function names instead of specializing to double precision function calls for Fortran intrinsic functions. E.g. use ''cos'' or ''sin'', not ''dcos'' or ''dsin'' (enforced).

  * ***MG_TOO_COMPLICATED_USEFASTLIBS_INSTEAD*** On many machines cache memory conflicts, which slow down performance, are avoided by not using dimension of arrays in multiples of 256. Thus it is recommended to avoid array lengths which are multiples of 256. This becomes especially important in the use of Fast Fourier Transforms. (usually enforced)

### 8. Exception handling and I/Os

MG: This part must be rewritten. We should stop the code via macros such as MSG_ERROR
but, at the same time, we should define some kind of protocol for critical Warnings
and unrecoverable failures... Low level routines should never stop but return an exit status
to be checked by the caller!
JB: Yes, we could print a stack of the error but this seem complicated with the actuel high level routine that alread have billions of arg and do zillions of thing

The code issues four types of special messages:

  * ERROR
  * BUG
  * WARNING
  * COMMENT

ERROR and BUG messages cause the code to stop, immediately, or after a very small delay. 
An ERROR is attributed to the user, while a BUG is attributed to the developer.
A WARNING message indicates that something happened that is not as expected,
but this something is not so important as to make the code stop.
A COMMENT message gives some information to the user, about something unusual. 
None of them should appear when the run is completely normal, and runs as expected (enforced).

Some checking should be provided for input variables for physical reasonableness or internal consistency.
This is intended to maintain a user-friendliness already present in the code, and will greatly improve the possibility of widespread use by new users. (usually enforced)

Exception diagnostics from routines should be directed to the standard output (std_out), and report:

  - the name of the routine which is reporting the error
  - the sort of exception that is occurring (ERROR, BUG, WARNING or COMMENT)
  - the present status
  - the expected status
  - in the case of an ERROR message or a BUG message, what is wrong, like "The volume of the unit cell should not be zero"
  - in the case of ERROR message, what action the user should take

(note that BUG messages are automatically appended with the directive "Action : contact ABINIT group").

Experience shows that the trouble spent putting good descriptive error diagnostics in place is 
well repaid by minimizing the time spent fixing trivial bugs. (usually enforced)

  -  (MG: Old way. Now we should always use the macros MSG_ERROR, MSG_BUG) Exiting should always be done through the LEAVE_NEW subroutine, except for temporary debugging purposes. The corresponding error message must be printed through the wrtout routine, since this routine can be used to modify the error message in a standard way when needed. See also the rules_parallel file, since exiting is much more complicated in the parallel case. (exit through leave_new is enforced, systematic use of wrtout not yet)

  - Never write a character in the first column of the main output file (ab_out). This is because the fldiff tool (see programmer_guide) uses the
signs in this column to perform comparison of test case output files with reference ones. (enforced)

  - When you output a message using the wrtout routine, do not use the Fortran90 / (Slash) editing. It is not treated correctly. Use instead  char(10), that will be handled in a special way
(note that the abbreviation 'ch10' for 'char(10)' is defined in basis_defs.f )
By contrast, do not use  char(10) when you output something using the Fortran90 write intrinsic :
char(10) is the end-of-line sign ONLY under UNIX.
For example MacOS needs char(13) ... So, use the Fortran90 / (Slash) editing in that case ... F or example:

<code>  write(std_out,'(a,/,a,/)') ' hello',' world !'</code>

  - When you write to the main output file, using the 'f' format of Fortran (like in f12.4), in case there
is a possibility of different rounding on different computers, add a small number (on the order of 100000 times smaller than the precision of the f format) to the number you want to write.
For example, be the number  numb=0.56250000000_dp , to be written using the format f11.3 :

<code>  write(ab_out,'(f11.3)')numb      ! will give portability problems</code>

On some computers, it will be written 0.562 , and on some others, it will be written 0.563 . Instead, use:

<code>  write(ab_out,'(f11.3)')numb+tol10</code>

as tol10 is defined 0.0000000001_dp in defs_basis.f. In this way, it will print as 0.563 always.

  - Use always the string "std_out" to write to the standard output file  (the "log" file, with unit number 6). 
Use always the string "ab_out" to write to the main output file (with unit number 7). 

    - Use always the string "std_err" to write to the standard error file (with unit number 0). 
      This is to allow correct redirection in case of massively  parallel run (beyond 100 processors): the writing can be globally suppressed.  Also, writing to the standard error file should only be done when the code is about to stop  gain to avoid problems with numerous write with massively parallel machines.
(enforced)

### 9. To be avoided 

  * Use

     <  instead of the old-fashioned .LT.
     <= instead of .LE.
     >  instead of .GT.
     >= instead of .GE.
     == instead of .EQ.
     /= instead of .NE.

  * Avoid RETURN at the end of each subroutine, use only END SUBROUTINE
. 
  * Avoid RETURN at the end of each function, use only END FUNCTION

  * Avoid DATA statements, initialize variables in their declaration line :

<code>       INTEGER :: nline = 6</code>
       (note: initialized variables get automatically the SAVE status)
       (note: there are arrays constructors to initialize arrays in their declaration line)

  * Avoid the use of COMMON and EQUIVALENCE

  * Do not use an isolated quote (') in a comment, like in the following :

 <code>      ! This isn't a good thing
       ! Indeed, the CPP (preprocessor) does not treat properly 
       ! such occurrences on some machines. Instead replace it simply with
       !
       ! This is not a good thing
       ! or
       ! This isn t a good thing
       ! Everybody will understand the latter, 
       ! although grammatically incorrect. (usually enforced)</code>
       
### 10. Use of BLAS and LAPACK subroutines 

BLAS and LAPACK subroutines are public domain subroutines, gathered
in libraries, that can be optimized by each hardware vendor on their platform. 
BLAS means Basic Linear Algebra Subroutines. There are different levels of BLAS:

  * BLAS 1: basic scalar and vector subroutines
  * BLAS 2: matrix-vector subroutines
  * BLAS 3: matrix-matrix subroutines
  * Sparse BLAS: sparse vector operations

LAPACK is a set of subroutine that supersedes the old LINPACK and EISPACK
libraries, and addresses Linear equation solution and Eigensystem solution
(matrix inversion, linear least squares problems, orthogonal factorization,
ordinary eigenvalue problems, generalized eigenvalue problems, singular value decomposition).

The Web addresses:

  * http://www.netlib.org/lapack/index.html

  * http://www.netlib.org/blas/index.html

For inversion or diagonalisation of small matrices (however larger than 3*3), 
LAPACK should be used preferentially to any other subroutines.

In general, the capabilities of BLAS subroutines can be obtained by the array capabilities of F90. 
The corresponding code is usually easier to read. For example,

  REAL(dp) :: aa(npw),bb(npw)
  CALL DCOPY(npw,aa,1,bb,1)

can be replaced by

  REAL(dp) :: aa(npw),bb(npw)
  bb(:)=aa(:)

So, the rule will be to use F90 instead of BLAS subroutines, with some exceptions (see below). 
In general, the following can be avoided : DSCAL, DCOPY, DAXPY, DZERO, DSUM, DMAX ..., as well as
their corresponding COMPLEX variant. (partially enforced)
  (JB: I don't know if it is somewhere so I put it here: Avoid kind specific routines/functions and use generic ones :<code>
    a=EXP(b)</code>
   instead of <code>
   a=DEXP(b)</code>
  This is more convenient
  )

There are exceptions :

  - constructs that would still be easier to read with BLAS calls,
  - constructs that are critical for speed, and faster with BLAS calls.

In the first class of routines, one find the scalar product between complex numbers, in our particular representation by

  REAL(dp) :: aa(2,size),bb(2,size)
  
Indeed, the F90 function DOT_PRODUCT(aa,bb) cannot be used, since it would not lead to the evaluation of a complex dot product, but a real dot product. So, use call DDOT(2*size,aa,1,bb,1)
(at present, no other BLAS routines than those called by LAPACK are used in ABINIT, so these rules are indications for the future)

### 11. Modules 

MG: This section should be rewritten from scratch!
JB: This section should be highlited somewhere so people use module instead of unpacked routines. And that would help every developers for debugging !

The architecture of ABINIT is modular, albeit not principally relying on F90 modules in its present status.
Routines are partitioned in different subdirectories, with a related hierarchy of usage. 
A more extensive use of modules at that level is the topics of current reflexion. 
See HM5_document, section 5 of the proposal by H. Matthis, in the present directory, ~abinit/doc/developers .

The module named 'defs_basis' contains global variables of ABINIT, and is used in every subroutine.
It contains the definitions of:

  * subtypes (e.g. kind of floating numbers)
  * I/O units
  * selected real constants or trigonometric constants
  * physical constants (conversion factors)

Modules are also useful for the definition of F90 interfaces (for optional arguments), or F90 derived datatypes (structured types).
Although they are allowed, no rules has been set up for the use of optional arguments, and calls based about keywords. 
This is a topics of current reflexion. 
Be cautious when you introduce such a feature, and mention it to the ABINIT coordinator.

### 12. Derived datatypes

12.a. Derived datatypes should be declared in the adequate module (MG: And this rule is not followed in many places, e.g dataset_type, header_type ...)
These are powerful F90 constructs, but the information about them is not local to the subroutine 
where they are used, so they should be introduced in a controlled way, in order for the programmers to become sufficiently easily familiarized with them: the introduction of a new derived datatype must be made in agreement with the coordinator of ABINIT.
The introduction of appropriate derived datatypes in ABINIT is one of the central issues of v3.2 and later versions.

12.b. Suffix a type identifier with '_type'. (MG: What about ''obj_t'' At present we use both!) For example:

    TYPE ipe_type              ! information on program execution
     INTEGER :: nb_proc
     CHARACTER*63 :: fname
    END TYPE ipe_type

which allows to define in subroutines::

   TYPE(ipe_type) :: ipe

12.c. Pros and cons.

Grouping connected variables into structured types is interesting for readability (it avoids too long
lists of formal arguments for instance) and may facilitate code changes inside subprograms.
Imagine you want to add a new general information and you want an access to this data in many subprograms :
if you insert a new field in 'ipe_type' structure, it will not be necessary to change each subprogram
interface.

However, source code itself may become less readable. Also, remember that the use of structured types is never more efficient for CPU: complex declarations should be avoided.

====== 13. Other topics ======

For the time being, pointers are only allowed when an array has to be allocated in a subprogram, and
then used and deallocated outside it. See rule 1.f. They should be nullified as soon as possible.

F95 constructs are not yet allowed 
It has been checked for ABINITv3.1 that not all compilers are able to deal with them.
MG_NEWRULE: The reference standard is F95, F2003 extension are allowed provided that:

   * all the slaves of the testfarm compile and pass the tests
   * It's possible to use these extensions either with wrapper functions or CPP macros so that the semantics is not changed when we have a compiler that do not support that particular F2003 feature. Example:<code>
 
#ifdef HAVE_FC_IOMSG
   open(file=trim(file),unit=unit, iostat=iostat, iomsg=iomsg)
#else
   open(file=trim(file),unit=unit, iostat=iostat)
   iomsg = "IOMSG not supported by the compiler"
#endif</code>

Object-oriented programming in F90 (definition of operators acting on objects...) is the subject 
of current reflexion. See for example the URL: http://www.cs.rpi.edu/~szymansk/oof90.html

### 13. Useful links 

On the WEB page 

http://www.cs.rpi.edu/~szymansk/oof90.html): "Introduction to Object-Oriented Concepts using Fortran90":

http://www.cs.rpi.edu/~szymansk/OOF90/Forum.html (shorter version)

http://www.Subsubsection title..
^^^^^^^^^^^^^^^^^^^^^cs.rpi.edu/~szymansk/OOF90/F90_Objects.html (PS file available)

(e.g. the data encapsulation is explained, taking a FFT routine as such an example. This might be related to the port using fftw...)
The other PS file ("How to Express C++ Concepts in Fortran 90") seems similar to the above notes.
Also, to the link "Gocha ! 
Click here for Fortran90 bug bites" on this page. This page warns the usage of the new functions in Fortran90.

"Compiler and tools tricks"
http://www.fortran-2000.com/ArnaudRecipes/CompilerTricks.html

"Arnaud's advice to develop and debug numerical applications"
http://www.fortran-2000.com/ArnaudRecipes/ArnaudAdvice.html


IDA: An interprocedural analyser for Fortran 90 and High Performance Fortran programs
http://www.vcpc.univie.ac.at/information/software/ida

Diagramf: a simple diagrammer for Fortran Language Programs

These were found from: http://www.fortranlib.com or http://www.uni-comp.com/fortran

Benchmarks of F90 compilers for Linux/Windows-boxes. See http://www.polyhedron.com

Summary tables for F90/95, and other Fortran information
http://www.owlnet.rice.edu/~mech517/F90_docs/tables.pdf
http://www.owlnet.rice.edu/~mech517/F90_OOP.html
http://www.owlnet.rice.edu/~mech517/books_new.html

### MG Additional comments

  * Use assumed length strings in low-level routines whenever possible. One can always add run-time checks to avoid string truncation. Example 
 
  * Prefer allocatable arrays over pointers, especially if you need a dynamic entity in a datatype. Allocatable arrays, unlike pointer, are contiguous and more efficient than pointers

  * Use => null() to nullify pointers in the declaration of the datatype (the initial status of a pointer is *undefined*

  * objects and methods should be defined in the same module

  * use capital case for global variables (what about global variables that are exported, e.g dp?)
    
  * Avoid magic numbers. Use if (optdriver==RUNL_GSTATE) instead of if (optdriver==1). Magic numbers render the code less readable and more difficult to modify/maintain
 
  * module names should start with m_ (what about the defs_ modules?)

  * Import a module with "use module only: ..." whenever possible, unless you are designing a low-level modules that export many constants and variables (e.g. m_xmpi). In the latter case, all the public entities of the module should start with the same prefix e.g xmpi_

  * Proposed conventions for standard methods: obj_init, obj_free, obj_nullify, obj_malloc, obj_print, obj_bcast, obj_ncwrite

  * new convention for types: foo_type -> foo_t

  * Rename the macros used for allocation and deallocation: ABI_ALLOCATE --> ABI_MALLOC, ABI_DEALLOCATE --> ABI_FREE, ABI_DATATYPE_ALLOCATE --> ABI_DT_MALLOC (the names presently used are TOOOOOOOOOOOOOOOOOOO long and readability matters!)
 
  * One COPYRIGHT section per file (we all know that ABINIT is a GPL code). If a file contains more that one procedure, the name of the author of the procedure can be reported in the ROBODOC section of the procedure or in the COPYRIGHT section reported at the beginning of the file.
 
  * Fortran include files should have the ".finc" file extension so that the abirules python scripts can easily locate them. Fortran90 files should use ".F90" as file extension. Other extensions are not accepted!
 
  * Avoid cluttered code:

       * Interlay code with comments and empty lines. Avoid huge blocks of statements without empty lines. 
       * Use helper functions to make the code more readable   

  * Use inlined comments wisely. A typical example where inlined comments are better than a separated comment line: <code>
         do spin=1,sppol
            ... a lot of code
         ! end spin loop <<<< This is bad!
         end do
           
         do spin=1,nsppol
            .... a lot of code
         end do ! end spin loop    <<< This is OK</code>

  * Avoid writing routines in which the *real* intent of the dummy arguments depend on some options. Long time ago we had to call xredxcart(red,cart,option) to convert from reduced to cartesian or vice versa depending on option. As a consequence both red and cart had to be defined with intent(inout). Now we have two routines for the different cases: xred2cart and xcart2red. The new approach is much more readable and we don't have to use intent(inout) for both arguments anymore. Similarly, avoid writing procedures that perform different kind of IO depending on the value of some input option Example:<code> 
  
    subroutine read_write(data, rdwr)
      real(dp),intent(inout) :: data(:)
      ! Very bad as any caller of read_write that 
      ! wants to write data to unit  must declare data 
      ! with intent(inout) and this renders 
      ! the code difficult to read and understand
     
      if (rdwr==1) then
         read(unit) data
       else if (rdwr==2) then
         write(unit) data
       end if  
       
     ! use
     subroutine read_data(data)
     subroutine write_data(data) 
     </code>
     
  * When writing modules, try to import all the required entities at the module level instead of using different use statements for the different module procedures. This approach reduces the amount of boilerplate code and clarifies the dependencies of the module. Use private as default attribute and export only those procedures or those entities that form the public API of the module. Example<code>

     module foo 
        use defs_basis
        use bar, only: bar_func
         
        implicit none ! Procedures *contained* in this module will inherit implicit none
        
        private  ! Everything is private unless explicitly stated
        public :: list_of_entities_exported_by_this_module
        
     contains
        ! All the procedures in this module will have access to defs_basis, bar_func1
        ! implicit none is automatically enforced   
     end module foo</code>
     
   * Avoid passing dataset and mpi_enreg to low-level routines in which you only need few parameters to control the algorithm. Passing dataset and mpi_enreg makes the code less reusable since you have to create a new big object from scratch if you want to reuse the code in another context. 

   * Write reusable code, split the algorithm in independent steps especially in the algorithm consists of a computational part and the IO section. Example: suppose you want to write a routine to compute the DOS and you want to save the results in formatted file. It doesn't make sense to perform these two tasks in the same procedure: write a routines that computes and returns the DOS and then wrap this procedure in another subroutine that calls the former one to get the DOS and then writes the results on file. The routine for the DOS can be reused in different contexts unlike getnel.F90

Work in progress:

Logging and Event handling: we distinguish between two classes of events: Critical Events and Normal Events:

Normal event: Something unexpected happened but the software is still working as expected (e.g. COMMENTs and WARNINGs)

Critical event: A serious error or condition, indicating that the program itself may be unable to continue running or that the results may be inaccurate or wrong (e.g. ERROR or an SCF cycle that is not converged)

Normal events are logged to the log file (std_out) while critical events require a more elaborate protocol so that users and/or an external software driving the calculation are **always** informed about the occurrence of such events. Note that this requirement is not easy to fulfill, especially when we have a MPI calculation on many processors or a run that uses MPI pools.

Normal Warnings, Comments are redirected to the std_out of the processor via the macros MSG_WARNING(msg), MSG_COMMENT(warning), respectively.

If a critical warning occurs on a MPI node, the process shall call the new macro ABI_CRITICAL(message, event_id) where message is a string and event_id is a unique identifier (integer > 0, the list of ids is defined in an include file and made available via the m_errors.F90 module) 

ABI_CRITICAL will try to report the event to ab_out. If the MPI process is not connected to the main output file, the logging will be redirected to the file __warnings__P[MPI_rank]__.txt and a file __abinit_warning__, containing the same message,is created (if the file already exists, no message is written, see below)
A similar approach is used in the case of Errors, the only difference being that the message is redirected to __errors__[MPI_RANK].txt and MPI_ABORT is called after the logging.

The file __abinit_warning__.txt is used in order to avoid problems with calculations done on many MPI nodes. In this case, indeed, the logging is disabled and no processors (except for master) can perform formatted IO. __abinit_warnings__ is mainly used for handling this particular case. Note that abinit should remove all the files named __abinit_warnings__ and __abinit_errors__ located in the working directory before starting the run (this approach won't work if the user is running several instances of abinit in the same directory and I won't do any attempt to address this particular case) 
### ABINIT style for MPI programming 

(Was valid before the band-by-band parallelism. Should be heavily updated!! In v4.2, see FFT_in_parallel )

This file is still under development. It should describe the rules
followed in the abinit code specifically related with the features that make the MPI parallel execution possible.

In the following, we will distinguish two different modes for the processes:

   * they are doing the same things: COLLECTIVE mode (i.e. same instructions and same data)

   * they are doing different things: PERSONAL mode (i.e. different instructions or different data)

### Basic rules for MPI parallelization 

  * Avoid defining variable whose name starts with MPI_ or mpi_ (Fortran is case insensitive). This prefix has been reserved by the MPI specifications for future additions.

  * Avoid explicit calls to the MPI routines as much as possible. Use the wrappers and the constants provided my m_xmpi.F90

  * Try to make the MPI algorithm consistent with the sequential version as much as possible. Think parallel since the very start (MPI is the common case, sequential is the exception). Avoid using HAVE_MPI or xmpi_paral==1 in high-level procedures as much as possible.

  * Simple MPI algorithms can be easily implemented once the MPI communicator is known.  Avoid passing mpi_enreg to your routine especially if you only need MPI_COMM, the number of processors and the CPU rank (use xcom_size, xcom_rank) NOTE: MPI_enreg clearly breaks the first rule in several procedures!

  * Use MPI barriers only when are really needed. Many of the MPI calls used in ABINIT are collective, thus a collective synchronization is implicitly performed. Remember: MPI barriers are expensive, especially when many CPUs are involved.

### Dealing with IO 

MG: The distinction between COLL and PERS is too complicated. Now each MPI node has its own log file
(logging can be easily disabled if nprocs > THRESHOLD), Hence the only subtle point is when we call wrtout(ab_out,"this goes to the main output file"). Well written code should insert these calls inside a if (rank==master) condition so that only the master node (i.e. the node who has access to ab_out) writes the important results to ab_out. I propose to remove "COLL and  "PERS" from wrtout) 

* Use subroutine wrtout(unit,message,mode) when writing to std_out or ab_out

* In the case of a COLLECTIVE mode, the 0-process will write the message in the
unit file, while the other processes will do nothing.

In the case of PERSONAL mode, we distinguish writing to the log file (standard output) from writing to other files.

  * In the case of a writing to the log file, the process will simply write the message preceeded by is proc-number.

  * The case of writing to other files should be avoided if possible. If not, the problem is still to be addressed.

### Dealing with stop 

  * MG: There's only one way to stop a MPI code when something unexpected happens i.e. MPI_ABORT. MPI_FINALIZE should be called only before leaving the main program.

  * Never use the Fortran stop explicitly (again MPI is the common case)

subroutine leave_new(unit,message,mode) + subroutine test_leave
In the case of a COLLECTIVE mode, the 0-process will write the message in the
unit file, while the other processes will do nothing. Then all processes will stop.

### ABINIT style for OpenMP programming 

The following sections are covered:

  * Basics
  * Portability
  * Efficiency

### Basics 

  * Information about OpenMP can be found at the following URL: (to be completed)
 
  * OMP is case insensitive, however the use of upper case keywords is strongly advised
 
  * !$OMP END PARALLEL DO is facultative in many cases. Avoid it whenever possible.
 
  * The default OMP name scope for variables is SHARED. In many cases one only needs to specify which variables are PRIVATE 

  * Do not use OMP to parallelize sections that may call MPI routines. OMP should be used for fine grained parallelization
      
  * Do not use OMP to parallelize loops that can be easily rewritten in terms of BLAS calls. Let the BLAS library parallelize the calculation for you, as a side effect the code will run faster even when threads are not used

### Portability 

  * OpenMp version 3.0 is assumed.

  * (IBM) named constants are not permitted in SHARED clause<code>

       integer, parameter :: mfac=11
       ! not accepted : mfac in a parameter
       !$OMP PARALLEL DO SHARED(mfac)      

       ! not accepted : two_pi is a parameter defined in defs_basis
       !$OMP PARALLEL DO SHARED(two_pi) </code>
                                            
  * Do not forget the continuation sign &<code>

      ! incorrect version: the continuation sign & is missing
      !$OMP PARALLEL DO DEFAULT(PRIVATE)   
      !$OMP&SHARED(aft,bef,gbound,g3max,ic,ind,len3,lot)

      ! correct version: the continuation sign & is present
      !$OMP PARALLEL DO DEFAULT(PRIVATE)&   
      !$OMP&SHARED(aft,bef,gbound,g3max,ic,ind,len3,lot) </code>
      
  * The intel compiler (version??) miscompiles COLLAPSE. Use the configure option --enable-omp-collapes to disable it. 
### Efficiency

(to be completed)