---
authors: AB,  OR
---

# the Fold2Bloch utility  

This file explains the use and i/o parameters needed for the "fold2Bloch" post-processor of the ABINIT package.  

This program generates an unfolded spectrum in a small cell, from a folded
spectrum of eigenvalues in a supercell. The new k-wavevectors and the weights
for the range of energies (standard postscript file **_fold2Bloch.out**_ ) are
based on the user's input of number and direction of folds as one of the
execution arguments (x:y:z). For the calculations program uses coefficients
and their vectors, K values, and eigenvalues found in the **__WFK_** file.
Fold2Bloch uses **_wfk_open_read, wfk_read_band_block and hdr_io_** routines
from the ABINIT package to read the _WFK file header and information about the
element structure which is used for our calculations.

It will be easier to discover the present file with the help of the tutorial ([[tutorial:fold2bloch]]).  

## 1 Compiling instructions
  
Compiled as part of abinit package. Fold2Bloch does not differentiate between
parallel and serial configuration so it always compiles for serial operation.

## 2 Execution
  
The program fold2Bloch is executed by invoking the following command _from the directory of the case file_ :

**fold2Bloch case_WFK x:y:z**

Where:

**case** is the name of the file that comes before _WFK identification  
 **x:y:z** are integer values( >0) that represent a multiplicity in the
corresponding direction used when constructing the super-cell; separated by a
colon (:).

Ex.

**fold2Bloch Ga8As7Bi1_WFK 2:2:2**

NOTE:  
If cell structure is spin polarized then there are two output files:  
**fold2Block_up.out** and **fold2Bloch_down.out**

## 3 Output sample
  
_**Output Sample**_

Below are three randomly selected parts from an output file. 
The data is presented as follows (from left to right columns):

New K Values (x, y, z) Eigenvalue(Ht) Weight

0.000000 0.000000 0.000000 -0.884864 1.000000

0.000000 0.000000 1.000000 -0.884864 0.000000

0.000000 0.000000 -1.000000 -0.884864 0.000000

0.000000 1.000000 0.000000 -0.884864 0.000000

0.000000 1.000000 1.000000 -0.884864 0.000000

0.000000 1.000000 -1.000000 -0.884864 0.000000

0.000000 0.000000 0.000000 -0.432062 0.000000

0.000000 0.000000 1.000000 -0.432062 0.500000

0.000000 0.000000 -1.000000 -0.432062 0.500000

0.000000 1.000000 0.000000 -0.432062 0.000000

0.000000 1.000000 1.000000 -0.432062 0.000000

0.000000 1.000000 -1.000000 -0.432062 0.000000

0.000000 0.000000 0.000000 -0.432062 0.000000

0.000000 0.000000 1.000000 -0.432062 0.500000

0.000000 0.000000 -1.000000 -0.432062 0.500000

0.000000 1.000000 0.000000 -0.432062 0.000000

0.000000 1.000000 1.000000 -0.432062 0.000000

0.000000 1.000000 -1.000000 -0.432062 0.000000

  

-0.333333 0.666667 0.500000 1.574788 0.499999

-0.333333 0.666667 1.500000 1.574788 0.000000

-0.333333 0.666667 -0.500000 1.574788 0.500001

-0.333333 -0.333333 0.500000 1.574788 0.000000

-0.333333 -0.333333 1.500000 1.574788 0.000000

-0.333333 -0.333333 -0.500000 1.574788 0.000000

-0.333333 0.666667 0.500000 1.574788 0.500001

-0.333333 0.666667 1.500000 1.574788 0.000000

-0.333333 0.666667 -0.500000 1.574788 0.499999

-0.333333 -0.333333 0.500000 1.574788 0.000000

-0.333333 -0.333333 1.500000 1.574788 0.000000

-0.333333 -0.333333 -0.500000 1.574788 0.000000

  

-0.333333 -0.166667 0.000000 -0.403296 1.000000

-0.333333 -0.166667 1.000000 -0.403296 0.000000

-0.333333 -0.166667 -1.000000 -0.403296 0.000000

-0.333333 0.833333 0.000000 0.040617 0.000000

-0.333333 0.833333 1.000000 0.040617 0.000000

-0.333333 0.833333 -1.000000 0.040617 0.000000

-0.333333 -0.166667 0.000000 0.040617 0.000000

-0.333333 -0.166667 1.000000 0.040617 0.500003

-0.333333 -0.166667 -1.000000 0.040617 0.499997

-0.333333 0.833333 0.000000 0.040617 0.000000

-0.333333 0.833333 1.000000 0.040617 0.000000

-0.333333 0.833333 -1.000000 0.040617 0.000000

-0.333333 -0.166667 0.000000 0.040617 0.000000

-0.333333 -0.166667 1.000000 0.040617 0.499997

-0.333333 -0.166667 -1.000000 0.040617 0.500003

## 4 Subroutines and functions

  
_**Subroutines and Functions**_

**Subroutine Getargs (getargs.F90)**

This routine read in the arguments from the command line, checks them for
correct format, and checks whether the input file exists. All the values are
assigned to their appropriate variables and are passed back to **fold2Bloch**.

**Parameters**

_Input:_

folds, fname

folds : empty array(1,3)  
fname : empty string

_Output:_

folds, fname

folds : array containing number of folds in x, y, and z directions. from
command argument  
fname : string containing the name of the input file. from the command
arguments

_Local:_

num_args, argcount, ii, ios, args, argfolds, dir

num_args : integer containing the number of command arguments sent  
argcount : integer counter used to iterate through all arguments  
ii : index used to indicate ':' position in the folds argument  
ios : integer used with iostat to check if folds are digits  
args : array containing all the arguments  
argfolds : temp string containing folds argument  
dir : boolean use to check if input file exists

**Subroutine Progress (Progress.F90)**

Progress compares K point currently being process with the total number of K
points and calculates percent complete. It then outputs the results and the
current K point that is being processed to the screen for the user to know.

**Parameters**

_Input:_

ikpt, nkpt, kpt

ikpt : current K point number  
nkpt : total number of K points  
kpt : current K point being processed

_Output:_

NONE

_Local:_

NONE

**Subroutine NewK (NewK.F90)**

This subroutine finds new K values for the unfolded spectrum. The K values are
determined based on how many times and in what direction the function needs to
be unfolded, which gives us the shape and size of the unfolded Brillouin zone.

**Parameters**

_Input:_

XX, YY, ZZ, FX, FY, FZ

XX : Original K point coordinate in the x-plane  
YY : Original K point coordinate in the y-plane  
ZZ : Original K point coordinate in the z-plane  
FX : Number of folds in the X direction  
FY : Number of folds in the Y direction  
FZ : Number of folds in the Z direction

_Output:_

NKVal

NKVal : 2D array containing x-y-z coordinates of new K points.

_Local:_

field_x, field_y, field_z, TX, TY, TZ, loop, ii, jj, kk, size

field_x : size of unfolded Brillouin zone in X direction  
field_y : size of unfolded Brillouin zone in Y direction  
field_z : size of unfolded Brillouin zone in Z direction  
TX : temporary holds the X coordinate of a new K point  
TY : temporary holds the Y coordinate of a new K point  
TZ : temporary holds the Z coordinate of a new K point  
loop : points to a place in NKVal to store a new K point  
I, j, k : counters for x, y, and z directions respectively  
size : number of new K points, used for NKVal allocation

**Subroutine SortC (SortC.F90)**

This subroutine sorts energy coefficients into appropriate groups and finds
each group's weight on the unfolded spectrum. Number of groups corresponds to
the number of new K points found in **NewK**. The coefficients are sorted
based on their relative position to the Gamma point.

**Parameters**

_Input:_

FX, FY, FZ, Vector, CoefC, NV

FX, FY, FZ : number of folds in each directions (see **NewK** )  
Vector : array containing x-y-z coordinates of each energy coefficient  
Coef : array containing energy coefficients for a specific Eigenvalue  
NV : number of vectors, and coefficients (one vector per coefficient)

_Output_

Weights

Weights : an array containing weights of each group

_Local_

TGroupc, Sums, counter, sumtot, remainder_x, remainder_y, remainder_z, jj, kk, ll, pp, el

TGroupC : an array containing groups of sorted coefficients  
Sums : an array containing sums of squared coefficients in each group  
counter : an array that keeps track of number of coefficients in each group  
sumtot : total sum of sums of each group  
reminader_x : relative position of each coefficient in x-axis  
reminader_y : relative position of each coefficient in y-axis  
reminader_z : relative position of each coefficient in z-axis  
jj, kk, ll : counters cycling through each group  
p : counter for elements in each group  
el : points to a location in Sums to store the next value

**MAIN (fol2Bloch.F90)**

**Parameters**

_Input_

fname, folds, nkval, nkpt, kpts, cg, kg, eig, weight, nsppol, nspinor, npwarr,
nband.

fname : case name inputted as part of the initial arguments by user, from
getargs()  
folds : initial argument (x:y:z), indicating number and direction of folds,
from getargs()  
nkval : array of new k points after the unfolding, from newk()  
nkpt : total number of K points, from hdr%nkpt  
kpts : K points coordinates, from hdr%kptns  
cg :  
kg : array of vectors per K point, from wfk_read_band_bloch() **  
**eig : array of eigenvalues per K point, from wfk_read_band_bloch() **  
**weight : weights calculate by sortc()  
nsppol : integer indicating up/down polarization (1 if none, 2 if polarized)  
nspinor : number of local orbitals, from hdr%nspinor  
npwarr : number of vectors in each kpoint, from hdr%npwarr  
nband : number of bands in K point, from hdr%nband  

_Output_

nkval, eig, weights  
(written to fold2Bloch.out or fold2Bloch_up.out and fold2Bloch_down.out if spin polarized)

nkval : new K points after unfolding  
eigval : original eigenvalues  
weights : weight of each new K point for a certain eigenvalue

_Local_

ikpt, iband, csppol, cg_b, count, infile, outfile, coefc, mband, mcg, outname

ikpt : integer number of current K point being processed  
iband : integer current band being processed  
csppol : integer current spin polarization being processed  
infile : integer input file number  
outfile : integer output file number  
coefc : array containing coefficient for the band being processed  
mband : integer maximum number of bands  
mcg : integer maximum possible number of coefficients to allocate them all  
 **  
**  
