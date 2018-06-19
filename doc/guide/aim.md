---
authors: PCasek, FF, XG
---

# the AIM utility  

This file explains the use and i/o parameters needed for the Atom-In-Molecule 
(AIM - Bader analysis) utility of the ABINIT package.  
The AIM utility allows to analyse charge densities produced by the ABINIT
code. The AIM analysis (Atom-In-Molecule) has been proposed by Bader. Thanks
to topological properties of the charge density, the space is partitioned in
non-overlapping regions, each containing a nucleus. The charge density of each
region is attributed to the corresponding nucleus, hence the concept of Atom- In-Molecule.

## 1 Introduction
  
The Bader technique allows to partition the space in attraction regions. 
Each of these regions is centered on one atom. The only input for this technique is
the total charge density: the density gradient line starting from one point
in space leads to one unique attracting atom. 
(References to the relevant literature are to be provided).

Around each atom, the basin of attraction forms a irregular, curved
polyhedron. Different polyhedra might have faces, vertices of apices in
common. Altogether, these polyhedra span the whole space.

The points where the density gradient vanishes are called "critical points" (CP). 
They all belong to the surface of some Bader polyhedra. According to the
number of positive eigenvalues of the Hessian of the density, one will distinguish:

  * "Bonding CPs" (a minimum along one eigenvector of the Hessian, but a maximum along the two other eigenvectors of the Hessian), 
     that belong to a face of the Bader volume (one BCP per face)

  * "Ring CPs" (a minimum along two eigenvectors of the Hessian, but a maximum along the remaining eigenvector of the Hessian), 
     that belong to a vertex of the Bader volume (one RCP per vertex)

  * "Cage CPs" (a global minimum, defining an apex of the Bader volume).

The Euler relation must be fulfilled: the number of BCPs minus the number of
RCPs plus the number of CCPs must equal 2. 
The Bader polyhedra might be convex (this is the usual case), but might as well not be convex.

In the present implementation, the user is required to specify one atom for
which he/she wants to compute the Bader volume, surfaces or critical points.
Different runs are needed for different atoms.

In case of the search for critical points, one start from the middle of the
segment with a neighbouring atom (all neighbouring atoms are examined in
turn), and evolves towards a nearby point with zero gradient. 

Then, in case [[aim:crit]] equals 2, one checks that the CP that has been found belongs to
the attraction region of the desired atom. This last step is by no means
trivial. In the present implementation, the check is done by means of the
straight line (radius) connecting the point with the atom. In case the Bader
volume is not convex, it might be that a correctly identified CP of the Bader
volume defines a radius that goes through a region that does not belong to the
Bader volume: the CP is "hidden" from the atom defining the attraction
region. In this case, the CP is considered as NOT being part of the Bader
volume, unfortunately. The reader is advised to look at the automatic test of
the MgO molecule to see such a pathology : no cage critical point is found for
the Mg atom. By chance, this problem is not a severe one, when the user is
interested by other aspects of the Bader analysis, as described below.

In case of the search for the Bader surface, or the integral of the charge
within the Bader surface, the user should define a series of radii of
integration, equally spread over the theta and phi angular variables.
Symmetries can be used to decrease the angular region to be sampled. Along
each of these radii, the algorithm will determine at which distance the radius
crosses the Bader surface. The integral of the charge will be done up to this
distance. For this search of the Bader surface, the information needed from
the critical points analysis is rather restricted : only an estimation of the
minimal and maximal radii of the Bader surface. This is why the fact that not
all CP have been determined is rather unimportant. On the other hand, the fact
that some part of the Bader surface might not be "seen" by the defining atom
must be known by the user. There might be a small amount of "hidden" charge as
well. Numerical tests show that the amount of hidden charge is quite small,
likely less than 0.01 electron.

The determination of density, gradient of density and hessian of the density
is made thanks to an interpolation scheme, within each (small) parallelepiped
of the FFT grid (there are n1*n2*n3 such parallelepiped). Numerical subtleties
associated with such a finite element scheme are more delicate than for the
usual treatment of the density within the main ABINIT code ! There are many
more parameters to be adjusted, defining the criteria for stopping the search
for a CP, or the distance beyond which CPs are considered different. The user
is strongly advised to experiment with the different parameters, in order to
have some feeling about the robustness of his/her calculations against
variation of these. Examples from the automatic tests should help him/her, as
well as the associated comments in the corresponding README files.

Note that the AIM program can also determine the Bader distance for one given
angular direction, or determine the density and laplacian at several, given,
points in space, according to the user will.

## 2 Input and output files
  
To run the program one needs to prepare two files:

  * [files-file] file which contains the name of the input file, the root for the names 
    of the different output files, and the name of other data files. 

  * [input-file] file which gives the values of the input variables 

Except these files you need the valence density in real space (*_DEN file,
output of ABINIT) and the core density files (*.fc file, output of the FHI
pseudopotential generator code, actually available from the ABINIT Web page)

The files file (called for example aim.files) could look like:
    
      aim.in    # input-file
      abo_DEN   # valence density (output of ABINIT)
      aim       # the root of the different output files
      at1.fc    # core density files (in the same order as
      at2.fc    # in the ABINIT files-file )
      ...
    
About the _DEN file:  
Usually, the grid in the real space for the valence density should be finer
than the one proposed by ABINIT. (For example for the lattice parameter
7-8~a.u. , ngfft at least 64 gives the precision of the Bader charge estimated
to be better than 0.003~electrons).

About the core density files:  
LDA core density files have been generated for the whole periodic table, and
are available on the ABINIT web site. Since the densities are weakly dependent
on the choice of the XC functional, and moreover, the charge density analysis
is mostly a qualitative tool, these files can be used for other functionals.
Still, if really accurate values for the Bader charge analysis are needed, one
should generate core density files with the same XC functional as for the valence density.

The main executable file is called aim. Supposing that the "files" file is
called aim.files, and that the executable is placed in your working directory,
aim is run interactively (in Unix) with the command

    aim < aim.files >& log

or, in the background, with the command

    aim < aim.files >& log &

where standard out and standard error are piped to the log file called "log"
(piping the standard error, thanks to the '&' sign placed after '>' is
**really important** for the analysis of eventual failures, when not due to
AIM, but to other sources, like disk full problem ...). 

The user can specify any names he/she wishes for any of these files. 
Variations of the above commands could be needed, depending on the flavor of UNIX that is used on the
platform that is considered for running the code.

The syntax of the input file is quite similar to the syntax of the main abinit
input files: the file is parsed, keywords are identified, comments are also
identified. However, the multidataset mode is not available.

Note that you have to specify what you want to calculate (default = nothing).
An example of the simple input file for Oxygen in bulk MgO is given in
~abinit/test/v3/Input/t57.in. 
There are also corresponding output files in this directory.

Before giving the description of the input variables for the aim input file,
we give some explanation concerning the output file.

The atomic units and cartesian coordinates are used for all output parameters.
The names of the output files are of the form root+suffix. There are different output files:

* [*.out] - the central output file - there are many informations which are clearly described. 
  Here is some additional information on the integration - there is separated integrations of the core and the valence density. 
  In both cases the radial integration is performed using cubic splines and the angular ones by Gauss quadrature. 
  However the principal part of the core density **of the considered atom** is integrated in the sphere 
  of minimal Bader radius using spherical symmetry. 
  The rest of the core density (of this atom, out of this sphere) together with all the core contributions 
  of the neighbors are added to the valence density integration. In the output file, there is the result 
  of the **complete** integration of the core density of the atom, then the two contributions 
  (spherical integration vs. the others) and then the total Bader charge. 

* [*.crit] - the file with the critical points (CPs) - they are listed in the order 
  Bond CPs (BCPs), Ring CPs (RCPs) and Cage CPs (CCPs). 
  The line of the output contains these informations: 
  $$ position \quad \frac{eigenvalues}{of\ Hessian} \quad \frac{index\ of}{bonded\ atom} \quad \Delta \rho_c \quad \rho_c $$ 
  (note: position is evaluated with respect to the considered atom, while the index of bonded atom is listed only for BCPs). 

* [*.surf] - the file with the Bader surface - there is a head of the form (in latex): 
  
```latex
      \begin{tabular}{ccc}
          index of the atom & \multicolumn{2}{c}{position} \\
          ntheta & thetamin & thetamax \\
          nphi & phimin & phimax \\
      \end{tabular}
```

and the list of the Bader surface radii: $$ \theta \quad \phi \quad
r(\theta,\phi) \quad W(\theta,\phi)$$ (note : $W(\theta,\phi)$ is the weight
for the Gauss quadrature). The minimal and maximal radii are given at the last line.

  * [*.log] - the log file - a lot of informations but it is not very nice actually. 
  * [*.gp] - gnuplot script showing the calculated part of the Bader surface with lines. 

The gnuplot scripts are made in the manner that one needs type only
    
    load  'file'

(quotes are necessary). Note, that this isn't considered as the visualization (it is only for working purpose)!

## 3 List of input variables
  
The list of input variables for the aim input file is presented in the [[varset:aim|aim set of input variables]].
