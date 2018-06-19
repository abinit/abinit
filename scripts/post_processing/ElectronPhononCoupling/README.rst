
ElectronPhononCoupling
======================

ElectronPhononCoupling (EPC) is a python module
to analyze electron-phonon related quantities computed with Abinit.


Istallation
-----------

Issue

    >$ python setup.py install

Requires

    * numpy >= 1.8.1
    * mpi4py >= 2.0.0
    * netCDF4 >= 1.2.1

Building the netCDF4 dependency is sometimes delicate. On certain systems,
it might be necessary to set the CC environment variable to the same compiler
that was used to build python. E.g. CC=gcc

Usage
-----

Example:

    import ElectronPhononCoupling as epc

    epc.compute(
        renormalization=True,
        broadening=True,
        self_energy=True,
        spectral_function=True,
        temperature=True,
        ...


You can run such python script in parallel with, e.g.:

    mpirun -n 4 python myscript.py

Documentation
-------------
 
* For how to use this module, see the Examples directory.

* For the theory pertaining the electronic self-energy
    due to electron-phonon coupling, and temperature dependence
    of electronic structure, see [PRB 92, 085137 (2015)].

* For the advanced user and developer, see the Doc directory.


