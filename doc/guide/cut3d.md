---
authors: XG, RC, GMR, JFB, MCote
---

# the Cut3D utility  

This file explains the use and i/o parameters needed for the
"Cut 3-Dimensional files" post-processor of the ABINIT package.

This code is able to analyse the files produced by ABINIT, that contain
3-Dimensional real space data, like all types of potential files, density
files. Wavefunction files data can also be analysed: first, a k-point number,
and the band number must be given, then, the corresponding wavefunction is transformed to real space.

In all these cases, thanks to Cut3D, one can obtain 2-Dimensional data
corresponding to a cut by a plane, or 1-Dimensional data along a line. One can
also translate the original formatting into many different ones.

Finally, one can also perform angular momenta analysis of wavefunctions with
respect to any given atom, computation of the Hirshfeld atomic charge (starting from a density file).

## 1 How to run cut3d
  
To run cut3d, simply type:
    
    cut3d
  
then, provide answer to the questions. You will have to give first the name of the unformatted file.
For example, t1xo_DEN.

## 2 Analyze most types of files, excluding wavefunction files
  
Supposing that you are not treating a wavefunction file,
you will have to choose between different possibilities:

1. computation of data for a point, to be specified
2. computation of data along a line, to be specified
3. computation of data on a 2D grid, to be specified
4. computation of data on a 3D grid, to be specified
5. conversion to formatted file
6. conversion to indexed formatted
7. conversion to Molekel format
8. conversion to 3D data with coordinates (tecplot ASCII format)
9. output .xsf file for XCrysDen
11. compute atomic charge using the Hirshfeld method

For option 1) you will have the possibility to specify a point in reduced or cartesian coordinates.
For option 2) you will have the possibility to specify a line by its two 
end-points in reduced or cartesian coordinates, or by it being perpendicular to some plane.
For option 3) and 4) many possibilities are offered, including specifications
thanks to points defined in reduced coordinates, cartesian coordinates, or atomic positions.  
To continue the analysis, simply answer the questions of the code, that should
be sufficiently self-explanatory.

## 3 Analyze wavefunction files
  
Instead, supposing that you are treating a wavefunction file, you will be able
to perform the analysis of one wave function. You will have to define the k
point, the number of the band, and possibly the spin-polarization or the spinor component.
Then, you will be asked whether you want to perform the angular component
analysis. You will have to provide the radius of the sphere(s) around each
atom, for which the angular analysis will be performed.
Finally, you will be given the choice between different formatting of the
wavefunction real-space data, including bare files, or XCrysDen formatted file.
