---
authors: JB
---

Documentation of `Xg_t` and `XgBlock_t` for abinit
==================================================

Written on: 2018/09/07  
Updated on:

Introduction and motivation
---------------------------

For the last decade, computers and more precisely processors have evolved a lot. 
To answer the issue of renewal, a new set of functions has been written to ease and adapt abinit to new architectures. 
This *abstract layer* is still a prototype and proof of concept. 
The true library should be written in an Oriented Object way to allow to switch during runtime between different libraries, or accelerators.

This current state of the code is stable. 
It is designed for **2D arrays** of **complex** or **real**
It is `Fortran2003` compliant and without any OOP.
Here is a list of what you can do

* Memory management

    - Allocate memory
    - Free memory
    - Block and sub-blocks
    - Map an object to an already allocated memory space
    - Get a `Fortran90` array from an object
    - Get/Set values from/to `Fortran90` array
    - Copy an object to an other one.
    - Pack a matrix to Upper/Lower Triangular Matrix
    - Reshape

* BLAS operations

    - POTRF
    - TRSM
    - GEMM
    - ADD
    - Scale

* LAPACK operations
    - HEEV
    - HEEVD
    - HPEV
    - HPEVD
    - HEGV
    - HEGVX
    - HPGVX
    - HPGVD

* ScaLAPACK
    - PHEEV
    - PHEGV

* Custom extensions
    - Array Shift
    - Colwise Norm 2
    - Colwise $Y-aX$
    - Set to 0
    - Set to Identity
    - Set the diagonal
    - Set the diagonal and zero off-diagonal terms
    - Get average
    - Get deviation
    - Print

* Tools
    - Get size
    - Get space

Each function is explained in detail below.

Memory management
-----------------
In this section, two types are explained. The first one, `xg_t`, is for allocating memory whereas the second one, `xgBlock_t` works on an already allocated object or pointer. I did not find an other way of doing thing in Fortran. Of course, in `C`, both types would be the same.

## Allocate memory
To create a new `xg_t`, one has to declare a `xg_t` and allocate its memory through the `xg_init` function:

```fortran
type(xg_t) :: xg1
call xg_init(xg1,space,rows,cols,comm)
```

`space` can either be `SPACE_C`, `SPACE_R` or `SPACE_CR`.
`SPACE_C` is for true complex numbers, `SPACE_R` is for real numbers taken from only the real part of complex numbers whereas `SPACE_CR` is consider both the real and imaginary parts of complex numbers as independent real numbers.
The difference between `SPACE_R` and `SPACE_CR` is only visible for getting and setting the values.
See the *Get/Set values* section for more details.  

`rows` would be the number of rows to consider and `cols` the number of columns.
`comm` is optional and is use if reductions are needed for algebra operations, and even ScaLAPACK if available.

## Free memory
Once you do not need an array anymore, or you do not need the memory anymore, you can free it with `xg_free`:

```fortran
call xg_free(xg1)
```

## Block and sub block

It may happen that you do not want to work on the full array but only a part of the array.
This is the main motivation to introduce the *block* notion.
A block is a part of an array on which you can work.
The main thing to remember is that a block work on an already allocated memory.
Therefore there is no check allocation nor freeing of block.
The `Fortran90` type for this object is `xgBlock_t` and you can create it with

```fortran
type(xgBlock_t) :: xgblockA
```

Basically, this is the main object you will use.
Even if you created a `xg_t` object to allocate memory, you cannot work on `xg_t` but rather on its own block
that you can access with `xg%self`.
`xg%self` is of type `xgBlock_t` and represents the full array you allocated.

To work on a subarray, you will need a `xgBlock_t` object, and build it with the `xg_setBlock` function

```fortran
type(xg_t)::xg1
type(xgBlock_t)::xgblockA
call xg_init(xg1,SPACE_C,maxRows,maxCols)
call xg_setBlock(xg1,xgblockA,firstColumn,rows,cols)
```

`xgblockA` represents a subarray of `xg1` starting from the first column `firstColumn`
(the very first one is 1) of `xg1` and has `rows` rows and `cols` columns.

## Map an object to an already allocated memory space

Usually, you will have a `Fortran90` array that you may want to plug in an `xgBlock_t` or *vice versa* you have a `xgBlock_t` object and want to get a `Fortran90` array variable.
You can achieve this operation with what is called *mapping* and *reverse mapping*.
From `Fortran90` array to `xgBlock_t` it is  *mapping*
On a regular cpu, this should be free, no allocations, no memory transfer, nothing.
But keep in mind that if you `xgBlock_t` represent an array on a GPU or other coprocessor unit, 
this will allocate memory and transfer it!

```fortran
double precision :: array(2,10)
type(xgBlock_t) :: xgblockA
! if array represents complex numbers, it has 10 complex numbers
call xgBlock_map(xgblockA,array,SPACE_C,rows,cols,comm)
! if array represents real numbers, it has 20 real numbers (SPACE_CR)
call xgBlock_map(xgblockA,array,SPACE_CR,rows,cols,comm)
```
`rows` and `cols` are the number of rows and columns for the representation inside the `xgBlock_t`.
It takes into account the space: regardless the space you chose, row*cols is the size of elements.
Example:
`array` if `Fortran90`

| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10|
|---|---|---|---|---|---|---|---|---|---|
| a | c | e | g | i | k | m | o | q | s |
| b | d | f | h | j | l | n | p | r | t |

```fortran
call xgBlock_map(xgblockA,array,SPACE_C,5,2,comm)
```

gives

|   1   |   2   |
|-------|-------|
| (a,b) | (k,l) | 
| (c,d) | (m,n) | 
| (e,f) | (o,p) | 
| (g,h) | (q,r) | 
| (i,j) | (s,t) | 

`(a,b)` is the complex number at position (1,1)

```fortran
call xgBlock_map(xgblockA,array,SPACE_CR,5,2,comm)
```

gives

|   1   |   2   |
|-------|-------|
|   a   |   f   | 
|   b   |   g   | 
|   c   |   h   | 
|   d   |   i   | 
|   e   |   j   | 

`a` is the real number at position (1,1)

## Get a `Fortran90` array from an object

From `xgBlock_t` to `Fortran90` array is a reverse mapping.
This is exactly the opposite operation as the mapping.

```fortran
double precision :: array(2,10)
type(xgBlock_t) :: xgblockA
! if array represents complex numbers, we take only the 10 first complex numbers
call xgBlock_reverseMap(xgblockA,array,rows,cols)
! if array represents real numbers, it has 20 real numbers (SPACE_CR) to put in a (2,10) array.
call xgBlock_reverseMap(xgblockA,array,rows,cols)
```

Be careful! If you work with complex, the array(2,10) has 10 complex numbers, only **one** row of complex numbers, not 2!!!  
Example:
`xgBlock_t` is (complex)

|   1   |   2   |
|-------|-------|
| (a,b) | (k,l) | 
| (c,d) | (m,n) | 
| (e,f) | (o,p) | 
| (g,h) | (q,r) | 
| (i,j) | (s,t) | 

```fortran
call xgBlock_reverseMap(xgblockA,array,1,10)
```

gives

| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10|
|---|---|---|---|---|---|---|---|---|---|
| a | c | e | g | i | k | m | o | q | s |
| b | d | f | h | j | l | n | p | r | t |

and
`xgBlock_t` is (real)

|   1   |   2   |
|-------|-------|
|   a   |   f   | 
|   b   |   g   | 
|   c   |   h   | 
|   d   |   i   | 
|   e   |   j   | 

```fortran
call xgBlock_reverseMap(xgblockA,array,2,5)
```

gives

| 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|
| a | c | e | g | i |
| b | d | f | h | j |

```fortran
call xgBlock_reverseMap(xgblockA,array,1,10)
```

gives

| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10|
|---|---|---|---|---|---|---|---|---|---|
| a | b | c | d | e | f | g | h | i | j |

You can access the space of an `xgBlock_t` and the dimensions with 

```fortran
! depending on the version it can be space(xgblockA)
space = getSpace(xgblockA)
! depending on the version it can be cols(xgblockA)
cols = getCols(xgBlock1)
! In a next version there will be a getRows(xgblockA) function too
call xgBlock_getSize(xgBlock1,rows,cols)
```

## Get/Set values from/to `Fortran90` array

Those are old functions that should not be used that much.
Before the mapping and reverse mapping processes, on had to copy values from `Fortran90` array to `xgBlock_t` and *vice versa*
Here it is very import to take care of the declared space for the `xgBlock_t`.
Furthermore, for abinit compliance, a `Fortran90` array is here always a 2D array with the first dimension equals to 2 `double precision :: array(2,:)`.
Each setter/getter exist for `xg_t` and `xgBlock_t`. 
Only the `xgBlock_t` version follows.

```fortran
double precision :: array(2,10)
type(xgBlock_t) :: xgblockA
... Manage memory
call xgBlock_set(xgblockA,array,shiftCol,rows)
```

Let's say array is

| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10|
|---|---|---|---|---|---|---|---|---|---|
| a | c | e | g | i | k | m | o | q | s |
| b | d | f | h | j | l | n | p | r | t |

`shiftCol = 1` and `rows = 5`.  
The result array in `xgblockA` will be
* if `SPACE_C`

|   1   |   2   |   3   |
|-------|-------|-------|
| xxxxx | (a,b) | (k,l) |
| xxxxx | (c,d) | (m,n) |
| xxxxx | (e,f) | (o,p) |
| xxxxx | (g,h) | (q,r) |
| xxxxx | (i,j) | (s,t) |

* if `SPACE_CR`
Here we set 5 rows of complex numbers so there are 5*2=10 rows of reals

|   1   |   2   |   3   |
|-------|-------|-------|
| xxxxx |   a   |   b   |
| xxxxx |   c   |   d   |
| xxxxx |   e   |   f   |
| xxxxx |   g   |   h   |
| xxxxx |   i   |   j   |
| xxxxx |   k   |   l   |
| xxxxx |   m   |   n   |
| xxxxx |   o   |   p   |
| xxxxx |   q   |   r   |
| xxxxx |   s   |   t   |

Just to be clear, the first half columns are the real parts and the second half parts are the imaginary parts.
This is absolutely not important since both are decoupled.

* if `SPACE_R`

Here we only take the reals parts of the number so there is actually 5 rows.

|   1   |   2   |   3   |
|-------|-------|-------|
| xxxxx |   a   |   k   |
| xxxxx |   c   |   m   |
| xxxxx |   e   |   o   |
| xxxxx |   g   |   q   |
| xxxxx |   i   |   s   |

The getter works to reverse the `set` operation.

## Copy an object to an other one.

If you want to copy a `xgBlock_t` into an other one, then you can use this functions.
As in the BLAS `copy` function, you can specify the incrementation for each object

```fortran
call xgBlock_copy(xgblockA,xgblockB,inc1,inc2)
```

This will copy `xgblockA` by incrementation of `inc1` into `xgblockB` by incrementation of `inc2`

## Pack a matrix to Upper/Lower Triangular Matrix
Packing is only usefull if you intent to use LAPACK functions that need packing.
Both upper and lower triangular packing are implemented.

```fortran
call xgBlock_pack(xgblockA,xgblockB,uplo)
```

`xgblockA` will be packed into `xgblockB` according to the value of `uplo`: *l* or *L* for lower packing and *u* or *U* for upper packing (See BLAS documentation for more explanations).

## Reshape

This tool is only here if the shape of the block should be changed.

```fortran
call xgBlock_reshape(xgblockA,newShape)
```

`newShape` is a 2 values array, the first one for the number of rows and the second one for the number of columns.

Here were presented the tools to allocate, free and manage memory blocks.
The next part will focus on BLAS functions and ease of use compare to bare BLAS functions.

BLAS operations
---------------

All Blas functions have input variable names identical to real BLAS calls. 
Missing arguments are not needed with this library.
Please refer to blas documentation at
[Netlib documentation](http://www.netlib.org/lapack/explore-html/dir_8ff1c41005a6034efef54f08b5a85d05.html)

## POTRF
## TRSM
## GEMM
## ADD
## Scale

LAPACK operations
-----------------
All Lapack functions have input variable names identical to real Lapack calls. 
Missing arguments are not needed with this library.
Please refer to blas documentation at
[Netlib documentation](http://www.netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen.html)

## HEEV
## HEEVD
## HPEV
## HPEVD
## HEGV
## HEGVX
## HPGVX
## HPGVD

ScaLAPACK
---------

!!! warning

    At the present time (2018/10), none of scalapack implementation is thread safe.
    ELPA should work but requires `MPI_THREAD_MULTIPLE` to work which is not yet available on our architectures! 
    So it is not tested and for security, scalapack/elpa is disabled at runtime with threads.

_______________________________

Before using a ScaLAPACK function, one has to check if it is possible or not (or event efficient).
The `xgScalapack_init` function is design for this purpose and is very similar to what was already done in the old implementation but in a more flexible way.

```fortran
use m_xgScalapack
type(xgScalapack_t) :: xgScalapack
integer :: comm
integer :: maxDim
integer :: verbosity
logical :: usable

call xgScalapack_init(xgScalapack,comm,maxDim,verbosity,usable)
```

`comm` is the MPI communicator to use. 
Note that not all of the processes available will be used.
`maxDim` is the maximum dimension of a matrix that will be used with this `xgScalapack` object.
`verbosity` is 0 or positive value.
If 0 then nothing is printed. 
If positive, then the number of processes inside the communicator `comm` is displayed.

To configure the number of processus that will be selected, you must invoke the function `xgScalapack_config` before calling `xgScalapack_init`.

```fortran
call xgScalapack_config(config,rankpp)
```

where `config` is one of `{SLK_AUTO,SLK_DISABLED}` or any positive value. 
A positive value will force scalapack to use exactly `config` processes.
In `SLK_AUTO` mode, `maxDim` will drive the number of processes to use on the fly.
Basically it is the maximum size of a matrix that a processes can hold.
For instance, if the full matrix is `2000x2000` and `maxDim=500`, then scalapack will use at most 4 processes.
When `config=SLK_DISABLED`, scalapack will not be used and `usable` will be set to `false`.
In addition, if the matrix is too small with respect to the configuration of the moddule, `usable` is set to `false`.
`usable` is set to true only if scalapack can be used.
In this case, subcommincators are created in the `xgScalapack` object.
Do not forget to call `xgScalapack_free(xgScalapack)` to free memory and destroy communicators.

With threads, `usable` is always set to `false`

Calling the following functions when `usable=false` is not tested and will crash (I think) or give trash results !

## PHEEV

This is equivalent to `HEEV` to solve `matrixA*X=eigenvalues*X`

```fortran
call xgScalapack_heev(xgScalapack,matrixA,eigenvalues)
```

`xgScalapack` is a previously initialized `xgScalapac_t` object, `matrixA` is the matrix to diagonalize and `eigenvalues` the eigenvalues.
Eigenvectors are placed in `matrixA` at the end.

## PHEGV

This is equivalent to `HEGV` to solve `matrixA*X=eigenvalues*matrixB*X`

```fortran
call xgScalapack_hegv(xgScalapack,matrixA,matrixB,eigenvalues)
```

Eigenvectors are placed in `matrixA` at the end.

Custom extensions
-----------------

## Array Shift

This function is the same as the standard fortran `cshift` function.
It is just an interface.

```fortran
  integer :: nshift
  integer :: shiftdim
  call xgBlock_cshift(xgBlock,nshift,shiftdim)
```

`nshift` is the shift to perform along the `shiftdim` dimension.
See fortran documentation for more information.

## Colwise Norm 2

This function is usefull to compute the square L2 norm -- meaning dot product -- of each column.
The function also provide information on min and max if there are present.

```fortran
  call xgBlock_colwiseNorm2(xgBlock,dot,max_val,max_elt,min_val,min_elt)
```

`dot` is the resulting dot product (BLAS call) for each column of `xgBlock`.
Therefore `dot` is a `SPACE_R` column vector with `xgBlock%cols` rows.
`max_val`, `min_val` are the maximal and minimal values of the `dot` vector and `max_elt`, `min_elt` are the indices of maximal and minimal values.
All those 4 variables are optional.

## Colwise $$Y-aX$$

The colwise operation done here is done column by columns.
It is a useful function to compute residuals for eigen problems.
Let suppose you solve HX=eX or HX=eBX, then calling this function will compute HX-eX or HX-eBX for each eigenvector of the matrix X.

```fortran
  call xgBlock_colwiseCymax(xgBlockA, da, xgBlockB,xgBlockW)
```

Here the operation is `xgBlockA(:,i) = xgBlockW(:,i)-da(i)*xgBlockB(:,i)`.
`da` is a `xgBlock_t` colum vector of `xgBlockA%cols` rows. 
In our example, `xgBlockA` is HX and `xgBlockW` is X or BX and `da` is the eigenvalues vector.
*Note:* this function is very similar to the BLAS `saxpy` function by the last one is not used
because the minus operation would require an array copy that I want to avoid.

## Set to 0

Nullify --set to 0-- all coefficient of the block.

```fortran
  call xgBlock_zero(xgBlock)
```

## Set to Identity

This function just set the diagonal to ones. 
The off-diagonal coefficients are not modified.
If you want to build the identity matrix, first call the `xgBlock_zero` function.

```fortran
  call xgBlock_one(xgBlock)
```

## Set the diagonal

Impose the diagonal of a matrix.

```fortran
  call xgBlock_diagonal(xgBlock,diag)
```

`diag` is of type `xgBlock_t` is is a column vector.
It is the diagonal to impose to `xgBlock` object.

## Set the diagonal and zero off-diagonal terms

Only keep the diagonal termes and set the off-diagonal terms to 0.

```fortran
  call xgBlock_diagonalOnly(xgBlock)
```

## Get average

Compute the average of the full block

```fortran
  call xgBlock_average(xgBlock,average)
```

`average` is the resulting real number.

## Get deviation

Compute the deviation of the full block.

```fortran
  call xgBlock_deviation(xgBlock,deviation)
```

`deviation` is the resulting real number.

## Print

As obvious as it seems, this function just print in a convenient way the data inside the block.

```fortran
  call xgBlock_print(xgBlock,outunit)
```

`xgBlock` is the block to print and `outunit` is the unit number in which the data will be written.
