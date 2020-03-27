#! /usr/bin/python

#
#    Copyright (C) 2003-2020 ABINIT group
#
#    Written by Gabriel Antonius in python (compatible v2.7).
#    This is free software, and you are welcome to redistribute it
#    under certain conditions (GNU General Public License,
#    see ~abinit/COPYING or http://www.gnu.org/copyleft/gpl.txt).
#
#    ABINIT is a project of the Universite Catholique de Louvain,
#    Corning Inc. and other collaborators, see ~abinit/doc/developers/contributors.txt.
#    Please read ~abinit/doc/biblio/generated_files/bib_acknow.html for suggested
#    acknowledgments of the ABINIT effort.
#
#    For more information, see https://www.abinit.org .

"""
This script can be run interactively,
but it is recommended to import it as a module:

    >>> from merge_ddb_nc import merge_ddb_nc
    >>> merge_ddb_nc(out_fname, fnames)

"""

from __future__ import print_function
import numpy as np
import netCDF4 as nc

__version__ = '1.0.0'

def merge_ddb_nc(out_fname, fnames):
    """
    Merge a list of DDB.nc files containing different elements of the same qpoint.

    Arguments
    ---------

    out_fname: Name for the merged file (will overwrite any existing file).
    fnames: List of DDB.nc files.

    """
    if not fnames:
        raise Exception('Empty list of files given for merge')

    fname0 = fnames.pop(0)

    with nc.Dataset(out_fname, 'w') as dsout:
        with nc.Dataset(fname0, 'r') as dsin:
            nc_copy(dsin, dsout)
            q0 = dsin.variables[u'q_point_reduced_coord'][...]

        for fname in fnames:
            with nc.Dataset(fname, 'r') as dsin:

                # Check that the qpoints are the same
                q = dsin.variables[u'q_point_reduced_coord'][...]
                if not all(np.isclose(q0, q)):
                    raise Exception('Cannot merge DDB.nc at different q-points.')

                # Merge dynamical matrix
                dynmat = dsin.variables[u'second_derivative_of_energy'][...]
                dynmat_mask = dsin.variables[u'second_derivative_of_energy_mask'][...]

                out_dynmat = dsin.variables[u'second_derivative_of_energy']
                out_dynmat_mask = dsin.variables[u'second_derivative_of_energy_mask']

                ni,nj,nk,nl = dynmat_mask.shape

                for i in range(ni):
                 for j in range(nj):
                  for k in range(nk):
                   for l in range(nl):
                    if dynmat_mask[i,j,k,l]:
                        dsout.variables[u'second_derivative_of_energy'][i,j,k,l,:] = (
                            dynmat[i,j,k,l,:])

                        dsout.variables[u'second_derivative_of_energy_mask'][i,j,k,l] = (
                            dynmat_mask[i,j,k,l])

                # Born effective charge tensor
                BECT = dsin.variables[u'born_effective_charge_tensor'][...]
                BECT_mask = dsin.variables[u'born_effective_charge_tensor_mask'][...]

                ni,nj,nk = BECT_mask.shape

                for i in range(ni):
                 for j in range(nj):
                  for k in range(nk):
                    if BECT_mask[i,j,k]:
                        dsout.variables[u'born_effective_charge_tensor'][i,j,k] = (
                            BECT[i,j,k])

                        dsout.variables[u'born_effective_charge_tensor_mask'][i,j,k] = (
                            BECT_mask[i,j,k])


def nc_copy(dsin, dsout):
    """
    Copy all dimensions and variable of one nc.Dataset instance into another.
    """

    #Copy dimensions
    for dname, dim in dsin.dimensions.iteritems():
        dsout.createDimension(dname, len(dim))

    #Copy variables
    for vname, varin in dsin.variables.iteritems():
        outVar = dsout.createVariable(vname, varin.datatype, varin.dimensions)
        outVar[...] = varin[...]


def interactive_merge_ddb_nc():
    """Get inputs from the user and run merge_ddb_nc."""

    program_name = 'merge_ddb_nc'
    description = """Merge several DDB.nc files, belonging to the same q-point."""

    def get_user(s):
        return raw_input(s.rstrip() + '\n').split('#')[0]


    print(program_name)
    print(len(program_name) * '-')
    print(description + '\n')

    ui = get_user('Enter a name for the output file in which to merge (will overwrite any existing file):')
    out_fname = str(ui)

    ui = get_user('Enter the number of files to merge:')
    nfiles = int(ui)

    fnames = list()
    for i in range(nfiles):
        ui = get_user('Enter the name of file {}:'.format(i+1))
        fname = str(ui)
        fnames.append(fname)

    # Main execution
    print('Executing...')
    merge_ddb_nc(out_fname, fnames)

    print('All done.')


# =========================================================================== #
# Run interactive program
# =========================================================================== #

if __name__ == '__main__':
    interactive_merge_ddb_nc()

