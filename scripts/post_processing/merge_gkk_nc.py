#! /usr/bin/env python

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

    >>> from merge_gkk_nc import merge_gkk_nc
    >>> merge_gkk_nc(out_fname, fnames)

"""

from __future__ import print_function
import numpy as np
import netCDF4 as nc

__version__ = '1.0.0'

def merge_gkk_nc(out_fname, fnames):
    """
    Merge a list of GKK<i>.nc files containing different elements of the same qpoint.

    Arguments
    ---------

    out_fname: Name for the merged file (will overwrite any existing file).
    fnames: List of GKK<i>.nc files.

    """
    if not fnames:
        raise Exception('Empty list of files given for merge')

    fname0 = fnames[0]

    with nc.Dataset(out_fname, 'w') as dsout:
        with nc.Dataset(fname0, 'r') as dsin:

            nc_copy(dsin, dsout,
                except_dimensions=['number_of_atoms_for_gkk', 'number_of_cartesian_directions_for_gkk'],
                except_variables=['second_derivative_eigenenergies_actif'],
                )

            q0 = dsin.variables['current_q_point'][...]
            natom = len(dsin.dimensions[u'number_of_atoms'])
            ncart = len(dsin.dimensions[u'number_of_cartesian_directions'])

        dsout.createDimension('number_of_atoms_for_gkk', natom)
        dsout.createDimension('number_of_cartesian_directions_for_gkk', ncart)

        gkk = dsout.createVariable('second_derivative_eigenenergies_actif', np.dtype('float64'),
                                    ('max_number_of_states', 'number_of_atoms_for_gkk',
                                     'number_of_cartesian_directions_for_gkk', 'number_of_kpoints',
                                     'product_mband_nsppol2'))

        for i, fname in enumerate(fnames):

            iatom = i // ncart
            icart = i % ncart

            with nc.Dataset(fname, 'r') as dsin:

                # Check that the qpoints are the same
                q = dsin.variables['current_q_point'][...]
                if not all(np.isclose(q0, q)):
                    raise Exception('Cannot merge GKK.nc at different q-points.')

                gkki = dsin.variables[u'second_derivative_eigenenergies_actif'][:,0,0,...]
                gkk[:,iatom,icart,...] = gkki


def nc_copy(dsin, dsout, except_dimensions=None, except_variables=None):
    """
    Copy all dimensions and variable of one nc.Dataset instance into another.
    """

    #Copy dimensions
    for dname, dim in dsin.dimensions.iteritems():
        if except_dimensions and dname in except_dimensions:
            continue
        dsout.createDimension(dname, len(dim))

    #Copy variables
    for vname, varin in dsin.variables.iteritems():
        if except_variables and vname in except_variables:
            continue
        outVar = dsout.createVariable(vname, varin.datatype, varin.dimensions)
        outVar[...] = varin[...]


def interactive_merge_gkk_nc():
    """Get inputs from the user and run merge_gkk_nc."""

    program_name = 'merge_gkk_nc'
    description = """Merge several GKK<i>.nc files, belonging to the same q-point."""

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
    merge_gkk_nc(out_fname, fnames)

    print('All done.')


# =========================================================================== #
# Run interactive program
# =========================================================================== #

if __name__ == '__main__':
    interactive_merge_gkk_nc()

