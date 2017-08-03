import numpy as np
import netCDF4 as nc

from .ncutil import nc_copy

def reduce_kpoints(subset, fname_in, fname_out):
    """
    Reduce the number of kpoints by selecting a subset of indices.
    Read the file fname_in and write the reduced data in fname_out.
    Use the name of fname_in to determine the type of file, with possible
    extentions: EIG.nc, EIGR2D.nc, EIGI2D.nc, GKK.nc, FAN.nc
    """
    if fname_in.endswith('EIG.nc'):
        return reduce_kpoints_eig(subset, fname_in, fname_out)
    elif fname_in.endswith('EIGR2D.nc') or fname_in.endswith('EIGI2D.nc'):
        return reduce_kpoints_eigr2d(subset, fname_in, fname_out)
    elif fname_in.endswith('GKK.nc'):
        return reduce_kpoints_gkk(subset, fname_in, fname_out)
    elif fname_in.endswith('FAN.nc'):
        return reduce_kpoints_fan(subset, fname_in, fname_out)
    else:
        raise Exception('Unrecognized file type for file {}'.format(fname_in))


def reduce_kpoints_eig(subset, fname_in, fname_out):
    """
    Operates on a EIG.nc file.
    Reduce the number of kpoints by selecting a subset of indices.
    Read the file fname_in and write the reduced data in fname_out.
    """

    # Name of the dimensions for nkpt in the nc file
    dname = 'nkpt'

    # Name of the variables with a nkpt dimension
    varnames = (
        'Eigenvalues',
        'Kptns',
        'NBandK',
        )

    reduce_dim(subset, dname, varnames, fname_in, fname_out)


def reduce_kpoints_eigr2d(subset, fname_in, fname_out):
    """
    Operates on a EIGR2D.nc or EIGI2D file.
    Reduce the number of kpoints by selecting a subset of indices.
    Read the file fname_in and write the reduced data in fname_out.
    """

    # Name of the dimensions for nkpt in the nc file
    dname = 'number_of_kpoints'

    # Name of the variables with a nkpt dimension
    varnames = (
        'reduced_coordinates_of_kpoints',
        'kpoint_weights',
        'number_of_states',
        'eigenvalues',
        'occupations',
        'istwfk',
        'second_derivative_eigenenergies',
        )

    reduce_dim(subset, dname, varnames, fname_in, fname_out)


def reduce_kpoints_gkk(subset, fname_in, fname_out):
    """
    Operates on a GKK.nc file.
    Reduce the number of kpoints by selecting a subset of indices.
    Read the file fname_in and write the reduced data in fname_out.
    """

    # Name of the dimensions for nkpt in the nc file
    dname = 'number_of_kpoints'

    # Name of the variables with a nkpt dimension
    varnames = (
        'reduced_coordinates_of_kpoints',
        'kpoint_weights',
        'number_of_states',
        'eigenvalues',
        'occupations',
        'istwfk',
        'second_derivative_eigenenergies_actif',
        )

    reduce_dim(subset, dname, varnames, fname_in, fname_out)


def reduce_kpoints_fan(subset, fname_in, fname_out):
    """
    Operates on a FAN.nc file.
    Reduce the number of kpoints by selecting a subset of indices.
    Read the file fname_in and write the reduced data in fname_out.
    """

    # Name of the dimensions for nkpt in the nc file
    dname = 'number_of_kpoints'

    # Name of the variables with a nkpt dimension
    varnames = (
        'reduced_coordinates_of_kpoints',
        'kpoint_weights',
        'number_of_states',
        'eigenvalues',
        'occupations',
        'istwfk',
        'second_derivative_eigenenergies_actif',
        )

    reduce_dim(subset, dname, varnames, fname_in, fname_out)


def reduce_dim(subset, dname, varnames, fname_in, fname_out):
    """
    Copy a netcdf file from fname_in into fname_out,
    but reduce the dimension 'dname' 
    """

    newdim = len(subset)

    with nc.Dataset(fname_in, 'r') as dsin:
        with nc.Dataset(fname_out, 'w') as dsout:
            nc_copy(dsin, dsout,
                    except_dimensions=[dname],
                    except_variables=varnames)

            dsout.createDimension(dname, newdim)

            for varname in varnames:

                datatype = dsin.variables[varname].datatype
                dimensions = dsin.variables[varname].dimensions
                axis = dimensions.index(dname)

                data = dsout.createVariable(varname, datatype, dimensions)
                data[...] = np.take(dsin.variables[varname][...], subset, axis)

