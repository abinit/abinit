Tutorial on Many-Body calculations in parallel

alpha-quartz SiO2 
-----------------------------------------

First test: Generation of the WFK file 
 
Standard GS + NSCF calculation to produce the WFK file
with a Gamma-centered 4x4x3 mesh (converged calculations
would require a much denser sampling e.g. 6x6x5 should be ok).

This test does not scale very well with the number of processors since we have
only 9 k-points in the IBZ and we don't use the paralkgb option. 

Also the number of bands in the WFK file is not enough to converge the GW 
calculations (1200 bands for the screening, 800 for sigma would be needed).

How to run the test:

(mpirun ...)  abinit < tmbt_1.files > log_1

then rename the WFK file

mv tmbt_1o_DS2_WFK 443_WFK

In the next steps, we will be using this file to start the GW calculations.

-------------------------------------------

Second step: Screening calculation with the Adler-Wiser formulation and gwpara=2.

First create a symbolic link for the WFK file

ln -s 443_WFK tmbt_2i_WFK

change the files file, then run the code as usual 

(mpirun ...)  abinit < tmbt_2.files > log_2

The RPA polarizabilty includes 26 empty bands therefore
this test should scale very well up to 26 processors.
Mind that, for optimal speedup, the number of processors should divide 26.

At the end of the calculation we have to rename the output SCR file 

mv tmbt_2o_SCR  443_ppm_SCR

The prefix _ppm means that this SCR file can be used for GW calculations
within the plasmon-pole approximation (only two frequencies have been computed).

-------------------------------------------


Third step: Screening calculation with the Hilbert transform method and gwpara=2

First create a symbolic link for the WFK file

ln -s 443_WFK tmbt_3i_WFK

change the files file specifying the correct index associated to the test,
then run the code as usual 

(mpirun ...)  abinit < tmbt_3.files > log_3

This test uses the same number of empty states as the previous step 
hence it will scale up to 26 processors.

At the end of the calculation we have to rename the output SCR file 

mv tmbt_3o_SCR  443_cd_SCR

The prefix _ppm means that this SCR file can be used for GW calculations
withing the contour-defomation technique.

-------------------------------------------


Fourth step: GW calculation with the plasmon-pole approximation and gwpara=2.

First create two symbolic links for the WFK and the SCR file

ln -s 443_WFK     tmbt_4i_WFK
ln -s 443_ppm_SCR tmbt_4i_SCR

then change the files as usual and issue

(mpirun ...)  abinit < tmbt_4.files > log_4

The number of bands is 50, hence the correlation part will scale up to 50 nodes.

The same input file can be used for performing CD calculations provided that we link the correct SCR file using

ln -s 443_cd_SCR tmbt_4i_SCR

and we replace the line

gwcalctyp 0  ppmodel

with

gwcalctyp 2  

in the input file

-------------------------------------------

Fifth step: GW calculation with the plasmon-pole approximation and gwpara=2.

First create a symbolic link for the WFK and the SCR file

ln -s 443_WFK     tmbt_5i_WFK
ln -s 443_ppm_SCR tmbt_5i_SCR

change the files file than run the code as usual 

(mpirun ...)  abinit < tmbt_5.files > log_5

This input file should scale up ... processors.
