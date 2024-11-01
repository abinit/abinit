# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

executable = "aim"

from abimkdocs.variables import ValueWithUnit, MultipleValue, Range
#from abipy.abio.abivar_database.variables import ValueWithUnit, MultipleValue, Range, ValueWithConditions
ValueWithConditions = dict
Variable=dict

variables = [
Variable(
    abivarname="atom@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_basic'],
    dimensions="scalar",
    defaultval=1,
    mnemonics="index of ATOM",
    added_in_version="before_v9",
    text=r"""
Index of the investigated atom.
""",
),

Variable(
    abivarname="atrad@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_expert'],
    dimensions="scalar",
    defaultval=1.0,
    mnemonics="bader ATomic RADius",
    added_in_version="before_v9",
    text=r"""
A first estimation of the Bader radius (not too important - it is used only
two times)
""",
),

Variable(
    abivarname="coff1@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_expert'],
    dimensions="scalar",
    defaultval=0.98,
    mnemonics="COeFFicient 1",
    added_in_version="before_v9",
    text=r"""
See the input variable [[ratmin@aim]].
""",
),

Variable(
    abivarname="coff2@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_expert'],
    dimensions="scalar",
    defaultval=0.95,
    mnemonics="COeFFicient 2",
    added_in_version="before_v9",
    text=r"""
See the input variable [[ratmin@aim]].
""",
),

Variable(
    abivarname="crit@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="computation of CRITical points",
    added_in_version="before_v9",
    text=r"""
Drives the computation of critical points.

  * [0] not
  * [-1] reading from the file ``root''.crit
  * [1] calculated (simplified version)
  * [2] calculated (standard version - recommended)
  * [3] calculated (the original version)

The original version searches all critical points (CPs) starting from the
center between two and three atoms (atom - neighbor(s)) by Newton-Raphson
algorithm - without tests (not recommended) - don't use together with surface
analysis !

The simplified and standard versions search CP(3,-1) starting from the center
of the pairs~atom-neighbor; then CP(3,1) from the center between two CP(3,-1)
and finally CP(3,3) from the center between two CP(3,1). The robust
Popeliers's algorithm is used. The difference between the two is based in the
fact that the standard version makes the test if the CP is really on the Bader
surface of the calculated atom for each CP, while the simplified version does
this only for CP(3,-1). When CP analysis is rather fast (with respect to
surface determination), 2 is recommended. In all cases the number of neighbors
considered is limited by distance cutoff (variable [[maxatd@aim]])
""",
),

Variable(
    abivarname="denout@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="electronic DENsity OUTput",
    added_in_version="before_v9",
    text=r"""
Output of the electronic density. The specification of the line (plane) in the
real space must be given in the input variable [[vpts@aim]] and grid in
[[ngrid@aim]]. It is also possible to get only the valence density or the core
density (see [[dltyp@aim]]).

  * 0, no output
  * 1, 1D distribution
  * 2, 2D distribution
""",
),

Variable(
    abivarname="dltyp@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Density or Laplacian TYP output",
    added_in_version="before_v9",
    text=r"""
Specification of the contribution of the electronic density corresponding to
the density and/or laplacian output (see [[denout@aim]] and [[lapout@aim]])

  * 0, total electronic density
  * 1, only the valence density
  * 2, only the core density
""",
),

Variable(
    abivarname="dpclim@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_useful'],
    dimensions="scalar",
    defaultval="1.d-2",
    mnemonics="DPCLIM",
    added_in_version="before_v9",
    text=r"""
If two "numerically different" critical points are separated by less than
**dpclim**, they are considered to be the same critical point. This often
happens because of numerical inaccuracies: one CP might be "seen" by two
different finite elements. The default should be OK when the ecut is quite
large, on the order of 60 Hartree. For less accurate calculations of the
density, increase the default value to 5.d-2, let's say.
""",
),

Variable(
    abivarname="foldep@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_expert'],
    dimensions=[3],
    defaultval="3*0.0",
    mnemonics="FOLlow DEParture",
    added_in_version="before_v9",
    text=r"""
Needed in the case [[aim:follow]]=1 only. Defines the starting point.
""",
),

Variable(
    abivarname="follow@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="FOLLOW the gradient path",
    added_in_version="before_v9",
    text=r"""
Follow the gradient path to the corresponding atom starting from the position
specified in the input variable [[aim:foldep]].
""",
),

Variable(
    abivarname="folstp@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_expert'],
    dimensions="scalar",
    defaultval=0.5,
    mnemonics="FOLlow STeP",
    added_in_version="before_v9",
    text=r"""
The first step for following the gradient path.
""",
),

Variable(
    abivarname="gpsurf@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="GraPhic output for the bader SURFace",
    added_in_version="before_v9",
    text=r"""
Drives the graphic output (gnuplot script) of the irreducible part of the
calculated Bader surface.

  * 0, not output
  * 1, output
""",
),

Variable(
    abivarname="inpt@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_expert'],
    dimensions="scalar",
    defaultval=100,
    mnemonics="numer of INtegration PoinTs",
    added_in_version="before_v9",
    text=r"""
Number of radial points used for integration of the Bader charge (not too
sensitive).
""",
),

Variable(
    abivarname="irho@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integration of the charge density RHO",
    added_in_version="before_v9",
    text=r"""
Drives the integration of the charge of the Bader atom.

  * 0, not calculated
  * 1, calculated (usual mode)
""",
),

Variable(
    abivarname="ivol@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="Integration of the VOLume",
    added_in_version="before_v9",
    text=r"""
Drives the integration of the volume of the Bader atom.

  * 0, not calculated
  * 1, calculated
""",
),

Variable(
    abivarname="lapout@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="electronic density LAPlacian OUTput",
    added_in_version="before_v9",
    text=r"""
Output of the laplacian of electronic density. The specification of the line
(plane) in the real space must be given in the input variable [[aim:vpts]] and
grid in [[aim:ngrid]]. It is also possible to get only the valence density or
the core density (see [[aim:dltyp]]).

  * 0, no output
  * 1, 1D distribution
  * 2, 2D distribution
""",
),

Variable(
    abivarname="lgrad@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_useful'],
    dimensions="scalar",
    defaultval="1.d-12",
    mnemonics="Low GRADient criterion",
    added_in_version="before_v9",
    text=r"""
The search for one particular CP is decided to be successful when either the
norm of the gradient of the electron density is smaller than **lgrad** or when
the length of the planned search step is smaller than [[aim:lstep]]. If the
number of search step becomes larger than an internal limit (presently set to
100), one will allow a weaker criteria for satisfaction, based on
[[aim:lgrad2]] and [[aim:lstep2]]. If the internal limit is reached, and the
criteria on [[aim:lgrad2]] and [[aim:lstep2]] are not satisfied, then the
searching procedure continues with the next seed.
""",
),

Variable(
    abivarname="lgrad2@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_useful'],
    dimensions="scalar",
    defaultval="1.d-5",
    mnemonics="Low GRADient criterion 2",
    added_in_version="before_v9",
    text=r"""
Determines the criterion for deciding that a CP has been found. See
[[aim:lgrad]] for more details.
""",
),

Variable(
    abivarname="lstep@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_useful'],
    dimensions="scalar",
    defaultval="1.d-10",
    mnemonics="Length of the planned search STEP",
    added_in_version="before_v9",
    text=r"""
Determines the criterion for deciding a CP has been found. See [[aim:lgrad]]
for more details.
""",
),

Variable(
    abivarname="lstep2@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_useful'],
    dimensions="scalar",
    defaultval="1.d-5",
    mnemonics="Length of the planned search STEP 2",
    added_in_version="before_v9",
    text=r"""
Determines the criterion for deciding that a CP has been found. See
[[aim:lgrad]] for more details.
""",
),

Variable(
    abivarname="maxatd@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_useful'],
    dimensions="scalar",
    defaultval=10.0,
    mnemonics="MAXimal ATomic Distance",
    added_in_version="before_v9",
    text=r"""
Atoms within this maximal distance are considered in order to start the search
of a CP.

Note that the supercell, determined by [[aim:nsa]], [[aim:nsb]], and
[[aim:nsc]] might be too small to actually lead to the consideration of all
the desired atoms.
""",
),

Variable(
    abivarname="maxcpd@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_useful'],
    dimensions="scalar",
    defaultval=5.0,
    mnemonics="MAXimal CP Distance",
    added_in_version="before_v9",
    text=r"""
The CPs are searched for within this maximal distance.

Note that the supercell, determined by [[aim:nsa]], [[aim:nsb]], and
[[aim:nsc]] might be too small to actually lead to the consideration of all
the critical points.
""",
),

Variable(
    abivarname="ngrid@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_expert'],
    dimensions=[2],
    defaultval="2*30",
    mnemonics="Number of GRID points",
    added_in_version="before_v9",
    text=r"""
Defines the grid in real space, for the density and laplacian outputs,
governed by [[aim:denout]] and [[aim:lapout]].
""",
),

Variable(
    abivarname="nphi@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_basic'],
    dimensions="scalar",
    defaultval=48,
    mnemonics="Number of PHI angle",
    added_in_version="before_v9",
    text=r"""
With [[aim:ntheta]], this variable defines the angular grid for the
integration within the Bader volume, in particular, the number of phi angles,
to be used between [[aim:phimin]] and [[aim:phimax]]. When the difference
between these two variables is 2 $\pi$, the recommended value of **nphi** is 48.
When it is $\pi$ (for symmetry reasons), the recommended value is 32. When it is
$\pi$/2 (for symmetry reasons), the recommended value is 20.
""",
),

Variable(
    abivarname="nsa@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_expert'],
    dimensions="scalar",
    defaultval=3,
    mnemonics="Number of Supercell points in direction A",
    added_in_version="before_v9",
    text=r"""
These variables define a "supercell", from the primitive cell repeated along
each primitive direction. This supercell is build as follows:



      do isa=-nsa,nsa
       do isb=-nsb,nsb
        do isc=-nsc,nsc
          -> here, the cell is translated by the vector
          -> (isa,isb,isc) in crystallographic coordinates
          -> and accumulated, to give the supercell
        enddo
       enddo
      enddo
""",
),

Variable(
    abivarname="nsb@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_expert'],
    dimensions="scalar",
    defaultval=3,
    mnemonics="Number of Supercell points in direction B",
    added_in_version="before_v9",
    text=r"""
These variables define a "supercell", from the primitive cell repeated along
each primitive direction. This supercell is build as follows:



      do isa=-nsa,nsa
       do isb=-nsb,nsb
        do isc=-nsc,nsc
          -> here, the cell is translated by the vector
          -> (isa,isb,isc) in crystallographic coordinates
          -> and accumulated, to give the supercell
        enddo
       enddo
      enddo
""",
),

Variable(
    abivarname="nsc@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_expert'],
    dimensions="scalar",
    defaultval=3,
    mnemonics="Number of Supercell points in direction C",
    added_in_version="before_v9",
    text=r"""
These variables define a "supercell", from the primitive cell repeated along
each primitive direction. This supercell is build as follows:



      do isa=-nsa,nsa
       do isb=-nsb,nsb
        do isc=-nsc,nsc
          -> here, the cell is translated by the vector
          -> (isa,isb,isc) in crystallographic coordinates
          -> and accumulated, to give the supercell
        enddo
       enddo
      enddo
""",
),

Variable(
    abivarname="ntheta@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_basic'],
    dimensions="scalar",
    defaultval=32,
    mnemonics="Number of THETA angles",
    added_in_version="before_v9",
    text=r"""
With [[aim:nphi]], this variable defines the angular grid for the integration
within the Bader volume, in particular, the number of theta angles, to be used
between [[aim:thetamin]] and [[aim:thetamax]]. When the difference between
these two variables is $\pi$, the recommended value of **ntheta** is 32. When it
is $\pi$/2 (for symmetry reasons), the recommended value is 20.
""",
),

Variable(
    abivarname="phimax@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_basic'],
    dimensions="scalar",
    defaultval=2.0,
    mnemonics="PHI MAXimal angle",
    added_in_version="before_v9",
    text=r"""
Angular limits of integration of the Bader volume for the phi variables. The
number of integration points is given by [[aim:nphi]]. The range of
integration can be decreased if there are symmetry reasons for doing this.
""",
),

Variable(
    abivarname="phimin@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_expert'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="PHI MINimal angle",
    added_in_version="before_v9",
    text=r"""
Angular limits of integration of the Bader volume for the phi variables. The
number of integration points is given by [[aim:nphi]]. The range of
integration can be decreased if there are symmetry reasons for doing this.
""",
),

Variable(
    abivarname="radstp@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_expert'],
    dimensions="scalar",
    defaultval=0.05,
    mnemonics="RADial STeP",
    added_in_version="before_v9",
    text=r"""
The length of the first step in the search of the exact Bader radius.
""",
),

Variable(
    abivarname="ratmin@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_expert'],
    dimensions="scalar",
    defaultval=1.0,
    mnemonics="Radius Atomic MINimal",
    added_in_version="before_v9",
    text=r"""
The first estimation of the smallest radius of the basin of the atom (the
distance at which the procedure that follows the gradient path announces that
the gradient path finishes in the corresponding atom) This parameter is very
important for the speed of the calculation, but this first estimation is not
usually used because the program makes a new one based on the knowledge of
CPs. In fact after the CP analysis, the new estimation is done by the product
of the ad hoc parameter [[aim:coff1]] (default 0.98) by the distance of the
nearest bonding CP. If there is a problem later, [[aim:coff2]] (default 0.95)
is used instead.
""",
),

Variable(
    abivarname="rsurdir@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_expert'],
    dimensions=[2],
    defaultval="2*0.0",
    mnemonics="Radius SURface DIRection",
    added_in_version="before_v9",
    text=r"""
In the case [[aim:rsurf]]=1, gives the direction (angular coordinates
theta,phi) along which the radius of the Bader surface is to be determined.
""",
),

Variable(
    abivarname="rsurf@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="computation of the Radius bader SURFace",
    added_in_version="before_v9",
    text=r"""
Drive the computation of the radius of the Bader surface for the angles
specified in the input variable [[aim:rsurdir]]

  * 0, not calculated
  * 1, calculated
""",
),

Variable(
    abivarname="scal@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_expert'],
    dimensions=[3],
    defaultval="1.0 1.0 1.0",
    mnemonics="SCALing of the cartesian coordinates",
    added_in_version="before_v9",
    text=r"""
SCALing of the cartesian coordinates.
""",
),

Variable(
    abivarname="surf@aim",
    varset="aim",
    vartype="integer",
    topics=['Bader_compulsory'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="computation of the bader SURFace",
    added_in_version="before_v9",
    text=r"""
Drive the computation of the full Bader surface.

  * 0, not calculated
  * 1, calculated
""",
),

Variable(
    abivarname="thetamax@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_basic'],
    dimensions="scalar",
    defaultval=r"$\pi$",
    mnemonics="THETA MAXimal angle",
    added_in_version="before_v9",
    text=r"""
Angular limits of integration of the Bader volume for the theta variables. The
number of integration points is given by [[aim:ntheta]]. The range of
integration can be decreased if there are symmetry reasons for doing this.
""",
),

Variable(
    abivarname="thetamin@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_expert'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="THETA MINimal angle",
    added_in_version="before_v9",
    text=r"""
Angular limits of integration of the Bader volume for the theta variables. The
number of integration points is given by [[aim:ntheta]]. The range of
integration can be decreased if there are symmetry reasons for doing this.
""",
),

Variable(
    abivarname="vpts@aim",
    varset="aim",
    vartype="real",
    topics=['Bader_useful'],
    dimensions=[6],
    defaultval="6*0.0",
    mnemonics="Vectors defining the PoinTS of the surface",
    commentdims="6 for 1D, 9 for 2D",
    added_in_version="before_v9",
    text=r"""
Basic vectors of the line or rectangle in real space, defining the points for
which the density or laplacian will be computed, thanks to [[aim:denout]] or
[[aim:lapout]]
""",
),

]
