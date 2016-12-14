# coding: utf-8

import os
from glob import glob

from setuptools import setup, find_packages

#---------------------------------------------------------------------------
# Basic project information
#---------------------------------------------------------------------------

# release.py contains version, authors, license, url, keywords, etc.
release_file = os.path.join('ElectronPhononCoupling', 'config', 'release.py')

with open(release_file) as f:
    code = compile(f.read(), release_file, 'exec')
    exec(code)


#---------------------------------------------------------------------------
# Find scripts
#---------------------------------------------------------------------------

def find_scripts():
    """Find the scripts."""
    scripts = []
    pyfiles = glob(os.path.join('ElectronPhononCoupling', 'scripts', "*"))
    scripts.extend(pyfiles)
    return scripts


def get_long_desc():
    with open("README.rst") as fh:
        return fh.read()

def write_manifest():
    content = """\
include *.rst
recursive-include ElectronPhononCoupling *.py
recursive-include ElectronPhononCoupling/scripts *
recursive-include ElectronPhononCoupling/data *
prune ElectronPhononCoupling/data/inputs-for-tests/output
graft Examples
exclude Examples/Calculations/*/odat*
exclude Examples/Calculations/*/*.out*
exclude Examples/Calculations/*/*.log*
exclude Examples/Calculations/*/*fort*
exclude Examples/Out/*
"""
    with open('MANIFEST.in', 'write') as f:
        f.write(content)


#---------------------------------------------------------------------------
# Setup
#---------------------------------------------------------------------------

write_manifest()

install_requires = [
    'numpy >=1.8',
    'mpi4py >=2.0',
    'netCDF4 >=1.2',
    ]

my_package_data = {
        'ElectronPhononCoupling.data' : ['data_*/*', 'inputs_for_tests/*.*', 'outputs_of_tests/*'],
    }

my_exclude_package_data = {
        'ElectronPhononCoupling.data' : ['inputs_for_tests/ouput/*'],
    }


setup_args = dict(
      name             = name,
      version          = __version__,
      description      = description,
      long_description = long_description,
      author           = author,
      author_email     = author_email,
      url              = url,
      license          = license,
      packages         = find_packages(),
      scripts          = find_scripts(),
      package_data     = my_package_data,
      exclude_package_data = my_exclude_package_data,
      install_requires = install_requires,
      )

if __name__ == "__main__":
    setup(**setup_args)

