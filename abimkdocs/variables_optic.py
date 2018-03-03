from __future__ import print_function, division, unicode_literals, absolute_import

executable = "optic"

from abimkdocs.variables import ValueWithUnit, MultipleValue, Range
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
    text="""
In Eq. 46 of [[cite:Sharma2004], it is clear that when ever wnm(k) is equal to w, there is a resonance.
Numerically this would lead to an infinity. In order to avoid this one could do two things.

You could change the sum over k-points to integration and then use
linear tetrahedron method (see [[cite:Hughes1996]] for details).
Another way to get around the problem is, like we do in the present case, avoid this
singularity by adding a small complex number to the denominator. This prevents
the denominator from ever going to 0 and acts as a broadening to the spectrum.
The broadening should not be too large as this would wash out the features in the spectrum.
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
    text="""
Specify the filename that has been produced by the preparatory Abinit run.
This file must contain the matrix elements of the d/dk operator along direction X.
It must not contain the first-order wavefunctions and may be generated using [[prtwf]] 3.
You should make sure that the number of bands, of spin channels and of
k-points are the same in all the files.

use as string with the filename: ddkfile_X, where X is the file number.
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
    text="""
The step and maximum sets your energy grid for the calculation using the
formula number of energy mesh points=maximum/step (zero excluded). So in order
to capture more features you can decrease the step size to get a finer energy
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
    text="""
This tells which component of the dielectric tensor you want to calculate.
These numbers are called a and b Eqs. 46 in [[cite:Sharma2004]].
1 2 3 represent x y and z respectively. For example 11 would be xx and 32 would mean zy.
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
    text="""
The step and maximum sets your energy grid for the calculation using the
formula number of energy mesh points=maximum/step (zero excluded). So in order
to capture more features you can decrease the step size to get a finer energy grid.
In order to go to higher frequency, increase the maximum.
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
    text="""
This tells which component of the dielectric tensor you want to calculate.
These numbers are called a, b and c in [[cite:Sharma2004]].
1 2 3 represent x y and z respectively. For example 111 would be xxx and 321 would mean zyx.
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
    text="""
How many components out of 9 of the linear optical dielectric tensor do you
want to calculate. Most of these are either equal or zero depending upon the
symmetry of the material (for detail see [[cite:Draxl2006]]).
Note that the directions are along the Cartesian axis.
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
    text="""
How many components out of 27 of the non-linear optical dielectric tensor do
you want to calculate. Most of these are either equal or zero depending upon
the symmetry of the material (for detail see [[cite:Draxl2006]]).
Note that the directions are along the Cartesian axis.
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
    text="""
LDA/GGA are known to underestimate the band-gap by up to 100%. In order
to get the optical spectrum and make a realistic comparison with experiments
one needs to correct for this. This can be achieved in two ways.

The scissors shift is normally chosen to be the difference between the experimental and
theoretical band-gap and is used to shift the conduction bands only. Another
way in which you do not have to rely on experimental data is to determine the
self energy using the [[lesson:gw1|GW approach]].
In this case the opening of the gap due to the GW correction can be used as scissor shift.
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
    text="""
When energy denominators are smaller than **tolerance** , the term is discarded from the sum.
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
    text="""
Specify the filename that has been produced by the preparatory Abinit run.
This file must contain the matrix elements of the d/dk operator along
direction X. It must not contain the first-order wavefunctions and may be
generated using [[prtwf]] 3.
You should make sure that the number of bands, of spin channels and of
k-points are the same in all the files.
""",
),

]
