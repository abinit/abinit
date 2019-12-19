# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

executable = "optic"

from abimkdocs.variables import ValueWithUnit, MultipleValue, Range
#from abipy.abio.abivar_database.variables import ValueWithUnit, MultipleValue, Range, ValueWithConditions
ValueWithConditions = dict
Variable=dict

variables = [
Variable(
    abivarname="broadening@optic",
    varset="optic",
    vartype="real",
    topics=['Optic_basic'],
    dimensions="scalar",
    defaultval="1.d-3 Ha",
    mnemonics="BROADENING",
    text=r"""
This parameter applies a broadening to the spectrum and is used to avoid
divergences in the sum-over-states approach.
The sum-over-states approach to the linear and nonlinear susceptibilities
inevitably include resonant denominators of the form
$(\omega - \omega_{nm})^{-1}$, see for example Eq. 46 of [[cite:Sharma2004]].
Numerically these denominators lead to infinities. In order to avoid them
one could do one of two things.
One could change the sum over k-points to integration, and then use the
linear tetrahedron method (see [[cite:Hughes1996]] for details).
Another way to get around the problem, as we will do in the present case,
is to avoid the singularities by adding a small imaginary contribution to the
denominator. This addition prevents the denominator from ever going to 0,
and acts as a broadening to the spectrum. The broadening should not be too
large, as this would wash out the features in the spectrum.
""",
),

Variable(
    abivarname="ddkfile@optic",
    varset="optic",
    vartype="string",
    topics=['Optic_basic'],
    dimensions="scalar",
    mnemonics="DDK FILE",
    commentdefault="no default",
    text=r"""
This parameter specifies the name of the file containing the matrix elements of the
$d/dk$ operator in direction X, as the string ddkfile_X. This file should have been
produced by a preparatory Abinit run.
This file must not contain the first-order wavefunctions, and may be generated
using [[prtwf]] 3. Make sure that the number of bands, spin channels, and
k-points are the same in all the files.
""",
),

Variable(
    abivarname="domega@optic",
    varset="optic",
    vartype="real",
    topics=['Optic_basic'],
    dimensions="scalar",
    defaultval="1.d-3 Ha",
    mnemonics="Delta OMEGA",
    text=r"""
This parameter specifies the step size $\Delta\omega$ for the grid over which the
optic utility computes the susceptibilities. The maximum energy is set by the
companion parameter [[optic:maxomega]]. The susceptibilities are thus computed at
[[optic:maxomega]]/[[optic:domega]] energy points (zero excluded). In order
to capture more features, decrease the step size to get a finer energy
grid. In order to go to higher frequency, increase the maximum.
""",
),

Variable(
    abivarname="lin_comp@optic",
    varset="optic",
    vartype="integer",
    topics=['Optic_basic'],
    dimensions=[['num_lin_comp']],
    defaultval=0,
    mnemonics="LINear COMPonents",
    text=r"""
This parameter specifies the directions of the [[optic:num_lin_comp]] requested components
of the dielectric tensor. The components are specified in
cartesian coordinates, where 1, 2, and 3 represent x, y, and z respectively. For
example, 11 represents the xx component, and 32 represents zy. There should be
[[optic:num_lin_comp]] entries. Note that these directions are denoted by
$a$ and $b$ in [[cite:Sharma2004]].
""",
),

Variable(
    abivarname="maxomega@optic",
    varset="optic",
    vartype="real",
    topics=['Optic_basic'],
    dimensions="scalar",
    defaultval="1 Ha",
    mnemonics="MAXimum value of OMEGA",
    text=r"""
This parameter specifies the maximum energy for the grid over which the
optic utility computes the susceptibilities. The grid step size is set by the
companion parameter [[optic:domega]]. The susceptibilities are thus computed at
[[optic:maxomega]]/[[optic:domega]] energy points (zero excluded). In order
to capture more features, decrease the step size to get a finer energy
grid. In order to go to higher frequency, increase the maximum.
""",
),

Variable(
    abivarname="nonlin_comp@optic",
    varset="optic",
    vartype="integer",
    topics=['Optic_basic'],
    dimensions=[['num_nonlin_comp']],
    defaultval=0,
    mnemonics="NON-LINear COMPonents",
    text=r"""
This parameter specifies the directions of the [[optic:num_nonlin_comp]] requested components
of the second-order nonlinear dielectric tensor. The components are specified in
cartesian coordinates, where 1, 2, and 3 represent x, y, and z respectively. For
example, 111 represents the xxx component, and 321 represents zyx. There should be
[[optic:num_nonlin_comp]] entries. Note that these directions are denoted by
$a$ and $b$ in [[cite:Sharma2004]].
""",
),

Variable(
    abivarname="num_lin_comp@optic",
    varset="optic",
    vartype="integer",
    topics=['Optic_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="NUMber of LINear COMPonents",
    text=r"""
This parameter specifies how many components (out of 9 possible)
of the linear optical dielectric tensor to calculate.
Some of these may be either equal to each other, or zero, depending upon the
symmetry of the material (for detail see [[cite:Draxl2006]]).
The directions of the requested components are specified by the parameter
[[optic:lin_comp]].
""",
),

Variable(
    abivarname="num_nonlin_comp@optic",
    varset="optic",
    vartype="integer",
    topics=['Optic_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="NUMber of NON-LINear COMPonents",
    text=r"""
This parameter specifies how many components (out of 27 possible)
of the second-order nonlinear optical dielectric tensor to calculate.
Some of these may be either equal to each other, or zero, depending upon the
symmetry of the material (for detail see [[cite:Draxl2006]]).
The directions of the requested components are specified by the parameter
[[optic:nonlin_comp]].
""",
),

Variable(
    abivarname="scissor@optic",
    varset="optic",
    vartype="real",
    topics=['Optic_basic'],
    dimensions="scalar",
    defaultval=0.0,
    mnemonics="SCISSOR operator",
    commentdefault="in Ha",
    text=r"""
This parameter provides a fixed shift to all the conduction bands. As
LDA/GGA are known to underestimate the band-gap by a significant amount in
some cases, in order to obtain a reasonable optical spectrum and make a realistic
comparison with experiments one needs to correct for this.
The scissors shift is normally chosen to be the difference between the experimental and
theoretical band-gap, and simply shifts the conduction bands. Alternatively, one may
determine the self energy using the [[tutorial:gw1|GW approach]], in which case
the opening of the gap due to the GW correction can be used as the scissor shift.
""",
),

Variable(
    abivarname="tolerance@optic",
    varset="optic",
    vartype="real",
    topics=['Optic_basic'],
    dimensions="scalar",
    defaultval="1.d-3 Ha",
    mnemonics="TOLERANCE",
    text=r"""
This parameter sets a scale for discarding small energy denominators.
When energy denominators are smaller than **tolerance** , the term is discarded from the sum.
See also [[optic:broadening]].
""",
),

Variable(
    abivarname="wfkfile@optic",
    varset="optic",
    vartype="string",
    topics=['Optic_basic'],
    dimensions="scalar",
    mnemonics="WaveFunction K FILE",
    commentdefault="no default",
    text=r"""
This parameter specifies the name of the ground state wavefunction file, which
should have been produced in a preparatory Abinit run. It should include both
the valence and conduction states to be used in the optic calculation
(see [[help:optic]]).
""",
),

]
