"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
]

#: List of keywords that are automatically added to all the tests of this suite.
keywords = [
]

#: List of input files
inp_files = [
# Constrained DFT
"t01.in",
"t02.in",
"t03.in",
# Structure variable
"t04.in",
# Optic
"t05.in", # GGA-PBE nsppol=1
"t06.in", # optic
"t07.in", # GGA-PBE nsppol=2
"t08.in", # optic
"t09.in", # GGA-PBE PAW nsppol=1
"t10.in", # optic
"t11.in", # GGA-PBE PAW nsppol=2
"t12.in", # optic
"t13.in", # metallic iron GGA-PBE PAW nsppol=2
"t14.in", # optic
"t15.in", # check slight misalignment of rprim, with large tolsym
"t16.in", # check slightly incorrect xred, with large tolsym
"t17.in", # check slightly incorrect rprim and xred, yielding correction to tnons, although stil not tolerated.
# GW/BSE
"t21.in", # HF exchange checking q->0 terms
"t22.in", # AC GW
# DFPT
"t41.in",
"t42.in",
# Electron-phonon
"t50.in",
"t51.in",
"t52.in",
"t53.in",
"t54.in",
"t55.in",
"t56.in",
"t57.in",
"t58.in",
"t59.in",
"t60.in",
"t61.in",
# Spin dynamics in multibinit
"t81.in", # set initial spin using rotation q-modulation
"t82.in", # damping
# New lattice mover in multibinit
"t83.in", # Langevin NVT 
"t84.in", # Berendsen NVT
"t85.in", # Velocity Verlet NVE
"t86.in", # Spin lattice coupling, Only Oiju term activated.
"t87.in", # Spin lattice coupling, Oiju and Tijuv. 
]
