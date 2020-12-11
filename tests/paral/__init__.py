"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
]

#: List of keywords that are automatically added to all the tests of this suite.
keywords = [
]

#: This suite contains tests executed with different numbers of MPI processes.
is_multi_parallel = True


#: List of input files
inp_files = [
"t01.abi",
"t02.abi",
"t03.abi",
"t05.abi",
"t06.abi",
"t07.abi",
"t08.abi",
"t09.abi",
"t21.abi",
"t22.abi",
"-t23.abi", # disabled
"t24.abi",
"t25.abi",
"t26.abi",
"t27.abi",
"t28.abi",
"t29.abi",
"t30.abi",
"t31.abi",
"t32.abi",
"t41.abi",
"-t48.abi", # MPI-IO, not tested here!
"-t49.abi", # MPI-IO, not tested here!
"t51.abi",
"t52.abi",
"t53.abi",
"t54.abi",
"t55.abi",
"t56.abi",
"t57.abi",
"-t58.abi", # disabled
"t59.abi",
"t60.abi",
"-t61.abi", # disabled for now
"t62.abi",
"t63.abi",
"t64.abi",
"t71.abi",
"t72.abi",
"t73.abi",
"t74.abi",
"t75.abi",
"t76.abi",
"t77.abi",
"t80.abi",
"t81.abi",
"t82.abi",
"t83.abi",
"t84.abi",
"-t85.abi", # disabled for now
"t86.abi",
"t91.abi",
"t92.abi",
"t93.abi",
"t94.abi",
"t95.abi",
"t96.abi",
"t97.abi",
"t98.abi",
"t99.abi"
]
