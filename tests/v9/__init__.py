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
"t18.in", # check slightly incorrect rprim and xred, yielding correction to tnons, tolerated.
"t19.in", # disable all checks thanks to expert_user
"t20.in", # test treatment of inaccurate POSCAR file
"t21.in", # test treatment of inaccurate POSCAR file
"t29.in", # RMM-DIIS eigsolver for NC.
"t30.in", # RMM-DIIS eigsolver for PAW.
# GW/BSE
"t31.in", # HF exchange checking q->0 terms
"t32.in", # AC GW
"t33.in", # GW 1RDM and related quantities 
"t34.in", # Same as t33.in but reading checkpoints
"t35.in", # GW 1RDM and related quantities (using only Sigma_x)
"t36.in", # GW 1RDM and related quantities but using Silicon
# DFPT
"t41.in",
"t42.in",
"t43.in",
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
"t62.in",
"t63.in",
"t64.in",
"t65.in",
# Spin dynamics in multibinit
"t81.in", # set initial spin using rotation q-modulation
"t82.in", # damping
# New lattice mover in multibinit
"t83.in", # Langevin NVT
"t84.in", # Berendsen NVT
"t85.in", # Velocity Verlet NVE
"t86.in", # Spin lattice coupling, Only Oiju term activated.
"t87.in", # Spin lattice coupling, Oiju and Tijuv.
# GS Coulomb cut-off
"t90.in", # checkG Coulomb cut-off, large tolerance a.t.m.
]
