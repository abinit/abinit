
ElectronPhononCoupling
======================

Python module to analyze electron-phonon related quantities from ABINIT.


Istallation
-----------

Simply issue

    >$ python setup.py install

This should install the module somewhere in your $PYTHONPATH
and the script "ElectronPhononCoupling/scripts/pp-temperature" in your $PATH

requires

    numpy >= 1.8.1
    mpi4py >= 2.0.0

Usage
-----

Interactive usage:

    >$ electron-phonon-coupling

or in parallel, e.g.:

    >$ mpirun -n 4 electron-phonon-coupling

As a python module:

    from ElectronPhononCoupling import compute_epc

    ...

You can run a python script that calls the function 'compute_epc' 
in serial or in parallel with e.g.:

    mpirun -n 4 python my_script.py

See the examples in ElectronPhononCoupling/data/inputs_for_tests/
for how to use this module.
The generation of the data with Abinit is explained in the
input file located in ElectronPhononCoupling/data/data_LiF/.

