
import netCDF4 as nc

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

