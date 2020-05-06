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
"t05.in",
"t06.in",
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
