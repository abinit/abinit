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
"t01.abi",
"t02.abi",
"t03.abi",
# Structure variable
"t04.abi",
# Optic
"t05.abi", # GGA-PBE nsppol=1
"t06.abi", # optic
"t07.abi", # GGA-PBE nsppol=2
"t08.abi", # optic
"t09.abi", # GGA-PBE PAW nsppol=1
"t10.abi", # optic
"t11.abi", # GGA-PBE PAW nsppol=2
"t12.abi", # optic
"t13.abi", # metallic iron GGA-PBE PAW nsppol=2
"t14.abi", # optic
"t15.abi", # check slight misalignment of rprim, with large tolsym
"t16.abi", # check slightly incorrect xred, with large tolsym
"t17.abi", # check slightly incorrect rprim and xred, yielding correction to tnons, although stil not tolerated.
"t18.abi", # check slightly incorrect rprim and xred, yielding correction to tnons, tolerated.
"t19.abi", # disable all checks thanks to expert_user
"t20.abi", # test treatment of inaccurate POSCAR file
"t21.abi", # test treatment of inaccurate POSCAR file
"t29.abi", # RMM-DIIS eigsolver for NC.
"t30.abi", # RMM-DIIS eigsolver for PAW.
# GW/BSE
"t31.abi", # HF exchange checking q->0 terms
"t32.abi", # AC GW
# DFPT
"t41.abi",
"t42.abi",
"t43.abi",
# Electron-phonon
"t50.abi",
"t51.abi",
"t52.abi",
"t53.abi",
"t54.abi",
"t55.abi",
"t56.abi",
"t57.abi",
"t58.abi",
"t59.abi",
"t60.abi",
"t61.abi",
"t62.abi",
"t63.abi",
"t64.abi",
"t65.abi",
# Spin dynamics in multibinit
"t81.abi", # set initial spin using rotation q-modulation
"t82.abi", # damping
# New lattice mover in multibinit
"t83.abi", # Langevin NVT
"t84.abi", # Berendsen NVT
"t85.abi", # Velocity Verlet NVE
"t86.abi", # Spin lattice coupling, Only Oiju term activated.
"t87.abi", # Spin lattice coupling, Oiju and Tijuv.
# GS Coulomb cut-off
"t90.abi", # checkG Coulomb cut-off, large tolerance a.t.m.
]
