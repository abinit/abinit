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
"t33.abi", # test npband>1 for Norm-Conserving
"t35.abi", # test xg_nonlop_option=1 (PAW,istwfk=1)
"t36.abi", # test cprj_in_memory=1 (PAW,istwfk>1,no DFTI)
"t37.abi", # test cprj_in_memory=1 (NC,istwfk>1,no DFTI)
"t38.abi", # test cprj_in_memory=1 (PAW,nspinor=2,compare xg_nonlop_option=0 and 1)
"t39.abi", # test cprj_in_memory=1 (PAW,istwfk=1,npkpt>1)
"t40.abi", # test cprj_in_memory=1 (PAW,nspinor=2)
"t41.abi",
"t43.abi", # test chebfi2 (paw+soc)
"t44.abi", # test chebfi2 (istwfk>1,paw,DFTI)
"t45.abi", # test chebfi2 (istwfk>1,norm conserving,DFTI)
"t46.abi", # test chebfi2 (istwfk>1,paw,not DFTI)
"t47.abi", # test chebfi2 (istwfk>1,norm conserving,not DFTI)
"t48.abi", # test cprj_in_memory=1 (PAW,istwfk>1,DFTI)
"t49.abi", # test cprj_in_memory=1 (NC,istwfk>1,DFTI)
#"t50.abi", # test cprj_in_memory=1 (NC,nspinor=2) ! not working yet
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
"t65.abi",
"t66.abi",
"t71.abi",
"t72.abi",
"t73.abi",
"t74.abi",
"t75.abi",
"t76.abi",
"t77.abi",
"t78.abi", # GWR
"t79.abi", # GWR
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
"t99.abi",
"t100.abi",
"t101.abi",
]
