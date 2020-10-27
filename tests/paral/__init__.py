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
"t01.in",
"t02.in",
"t03.in",
"t05.in",
"t06.in",
"t07.in",
"t08.in",
"t09.in",
"t21.in",
"t22.in",
"-t23.in", # disabled
"t24.in",
"t25.in",
"t26.in",
"t27.in",
"t28.in",
"t29.in",
"t30.in",
"t31.in",
"t41.in",
"-t48.in", # MPI-IO, not tested here!
"-t49.in", # MPI-IO, not tested here!
"t51.in",
"t52.in",
"t53.in",
"t54.in",
"t55.in",
"t56.in",
"t57.in",
"-t58.in", # disabled
"t59.in",
"t60.in",
"-t61.in", # disabled for now
"t62.in",
"t63.in",
"t64.in",
"t71.in",
"t72.in",
"t73.in",
"t74.in",
"t75.in",
"t76.in",
"t77.in",
"t80.in",
"t81.in",
"t82.in",
"t83.in",
"t84.in",
"-t85.in", # disabled for now
"t86.in",
"t91.in",
"t92.in",
"t93.in",
"t94.in",
"t95.in",
"t96.in",
"t97.in",
"t98.in",
"t99.in"
]
