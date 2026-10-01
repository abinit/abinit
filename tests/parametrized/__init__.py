"""
Global variables associated to the parametrized test suite.

This suite is optional: it is registered only with `runtests.py --with-parametrized` or `with_parametrized`
in testbot.json, and requires the pyyaml and pydantic packages.
Input/ contains templates (*.tmpl.abi) that are expanded on the fly into tests:
see tests/plans/parametrized_tests_plan.md.
"""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
]

#: List of keywords that are automatically added to all the tests of this suite.
keywords = [
]

#: True as soon as a template sweeps nproc, so that testbot runs it in the multi-parallel loop (plan 9.1).
is_multi_parallel = True

#: List of templates. The fixtures (*.fixture.abi) used by t07, and the prerequisites
#: of a chain (setup.test_chain), are not listed.
inp_files = [
"t01.tmpl.abi",
"t02.tmpl.abi",
"t03.tmpl.abi",
"t04.tmpl.abi",
"t05.tmpl.abi",
"t06.tmpl.abi",
"t07.tmpl.abi",
]
