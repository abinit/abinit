"""
Run a script as if pydantic were not installed (non-regression checks, plan section 11.3):

    python tests/pymods/parametrized/without_pydantic.py tests/runtests.py v1 -j6
"""
import runpy
import sys

sys.modules["pydantic"] = None  # every `import pydantic` now raises ImportError
script, sys.argv = sys.argv[1], sys.argv[1:]
runpy.run_path(script, run_name="__main__")
