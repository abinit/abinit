#!/usr/bin/env python
"""
This module provides functions and objects to validate netcdf files written in the ETSF-IO file format.
For a quick reference to the etsf specs see: http://esl.cecam.org/mediawiki/index.php/ETSF_File_Format_Specifications
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import re
import logging
logger = logging.getLogger(__name__)

from .termcolor import cprint

try:
    import numpy as np
    import netCDF4
except ImportError as exc:
    errmsg = str(exc) + "\nCannot import numpy or netCDF4. Use `anaconda or pip install netcdf`\n"
    logger.warning(errmsg)
    raise ImportError(errmsg)


def all_subclasses(cls):
    """
    Given a class `cls`, this recursive function returns a list with
    all subclasses, subclasses of subclasses, and so on.
    """
    subclasses = cls.__subclasses__() 
    return subclasses + [g for s in subclasses for g in all_subclasses(s)]


class EtsfObject(object):
    """
    Base class for netcdf Dimensions, Variables, Attributes.
    Subclasses implement a `validate` method that receives a nc dataset 
    and return a list of errors (strings).
    """
    def __str__(self):
        return "<%s: %s>" % (self.__class__.__name__, self.name)


class EtsfDimension(EtsfObject):
    """A dimension has a name, a type and, optionally, a list of allowed values."""
    def __init__(self, name, xtype, allowed=None):
        self.name = name
        self.xtype = xtype
        self.allowed = allowed

    def validate(self, ncdata):
        """Validate the content in the nc dataset, return list of errors."""
        errors = []
        eapp = errors.append

        # Dimension should be present.
        if self.name not in ncdata.dimensions:
            eapp("<%s> Not present in ncdata" % self.name)
            return errors

        # Check if value is allowed.
        if self.allowed is not None:
            n = len(ncdata.dimensions[self.name])
            if n not in self.allowed:
                eapp.append("Value %s not in allowed %s" % (n, self.allowed))

        return errors


class EtsfAttribute(EtsfObject):
    """A dimension has a name, a type and, optionally, a shape and list of allowed values."""
    def __init__(self, name, xtype, shape=None, allowed=None):
        self.name = name 
        self.allowed = allowed

    def validate(self, ncdata):
        """Validate the content in the nc dataset, return list of errors."""
        errors = []
        eapp = errors.append
        nc_attrs = ncdata.ncattrs()

        # Dimension should be present.
        if self.name not in nc_attrs:
            eapp("<%s> not present in ncdata" % self.name)

        # Check if value is allowed.
        if self.allowed is not None:
            value = nc_attrs[self.name] 
            if value not in self.allowed:
                eapp("Value %s is not allowed" % str(value))

        return errors


class EtsfVariable(EtsfObject):
    """
    A variable has a name, a type and a list of dimensions.
    The list of allowed values and the attributes that must be specified are optional.
    """
    # Stores all the instances we are gonna create.
    all_variables = []

    def __init__(self, name, xtype, dimensions, allowed=None, reqattrs=None):
        self.name = name
        self.xtype = xtype
        #if xtype == "char": assert shape
        self.dimensions = dimensions
        #self.shape = tuple(shape) if shape else ()
        self.allowed = allowed
        self.reqattrs = reqattrs
        # TODO dtype

        self.ndim = len(self.dimensions)

        # Store the instance in cls.all_variables
        self.__class__.all_variables.append(self)

    def validate(self, ncdata):
        """Validate the content in the nc dataset, return list of errors."""
        errors = []
        eapp = errors.append

        # Variable should be present.
        if self.name not in ncdata.variables:
            eapp("<%s> not present in ncdata" % self.name)
            return errors

        # If variable is present, check that dimension names agree with the specifications.
        ncvar = ncdata.variables[self.name]

        # Number of dimensions
        if self.ndim != ncvar.ndim:
            eapp("Has different ndim:\n\tFile: %s\n\tSpecs: %s" % (self.ndim, ncvar.ndim))
            return errors

        # Test shape
        shape = np.array([len(ncdata.dimensions[dim.name]) for dim in self.dimensions])
        if np.any(shape != ncvar.shape):
            eapp("Has different shape" % (shape, ncvar.shape))
            return errors

        # Test dimension names.
        dim_names = np.array([dim.name for dim in self.dimensions])
        ncdim_names = np.array(ncvar.dimensions)
        #print("dim_names: ", dim_names); print("ncdim_names: ", ncdim_names)
        if any(dim_names != ncdim_names):
            eapp("Wrong dimension names.\n\tFile: %s\n\tSpecs: %s" % (dim_names, ncdim_names))

        # Test allowed values.
        if self.allowed is not None:
            assert not ncvar.shape 
            try:
                value = ncvar.getValue()[0] if not ncvar.shape else ncvar[:]
            except IndexError:
                value = ncvar.getValue() if not ncvar.shape else ncvar[:]

            if value not in self.allowed:
                eapp("Value %s not in allowed %s" % (value, self.allowed))

        # Test attributes.
        if self.reqattrs is not None:
            nc_attrs = ncvar.ncattrs()
            for attr in self.reqattrs:
                aname = attr.name
                if aname not in nc_attrs:
                    eapp("Attribute %s is missing for variable: %s" % (aname, self.name))

        return errors


class VariableWithUnits(EtsfVariable):
    """A Variable that requires the specification of units."""
    def __init__(self, name, xtype, dimensions, allowed=None, reqattrs=None):
        if reqattrs is None: reqattrs = []
        reqattrs = reqattrs[:]
        reqattrs.append(units)

        super(VariableWithUnits, self).__init__(name, xtype, dimensions, allowed=allowed, reqattrs=reqattrs)

    def validate(self, ncdata):
        """Validate the content in the nc dataset, return list of errors."""
        errors = super(VariableWithUnits, self).validate(ncdata)

        if self.name not in ncdata.variables:
            assert errors
            return errors

        ncvar = ncdata.variables[self.name]
        nc_attrs = ncvar.ncattrs()
        if "scale_to_atomic_units" not in nc_attrs:
            errors.append("scale_to_atomic_units attribute must be specified") 

        return errors

# AND NOW, LADIES & GENTLEMEN, THE ETSF-IO SPECIFICATIONS

######################################################################################
# Mandatory attributes
# This table gathers specifications for required attributes in any ETSF NetCDF files.
######################################################################################
file_format = EtsfAttribute("file_format", "char", [80])
file_format_version = EtsfAttribute("file_format_version", "real")
Conventions = EtsfAttribute("Conventions", "char", [80])

mandatory_attributes = [file_format, file_format_version, Conventions]

#################################################################
# Optional attributes
# This table presents optional attributes for ETSF NetCDF files.
#################################################################
history	= EtsfAttribute("history", "char", [1024])	
title = EtsfAttribute("title", "char", [80])	

k_dependent = EtsfAttribute("k_dependent", "char", [80], allowed=["yes", "no"])	
# Attribute of number_of_states, flag-type, see Flag-like attributes.

symmorphic = EtsfAttribute("symmorphic", "char", [80], allowed=["yes", "no"])
# flag-type attribute, see Flag-like attributes.

##################################################################
# Generic attributes of variables
# A few attributes might apply to a large number of variables. 
# The following table presents the generic attributes that might 
# be mandatory for selected variables in ETSF NetCDF files.
##################################################################
units = EtsfAttribute("units", "char", [80])
scale_to_atomic_units = EtsfAttribute("scale_to_atomic_units", "double")

##############################################################################
# Dimensions that cannot be split
# This table list the dimensions that are not supposed to lead to a splitting.
##############################################################################
character_string_length = EtsfDimension("character_string_length", "integer", allowed=[80])
real_or_complex_coefficients = EtsfDimension("real_or_complex_coefficients", "integer", allowed=[1, 2])
real_or_complex_density	= EtsfDimension("real_or_complex_density", "integer", allowed=[1, 2])
real_or_complex_gw_corrections = EtsfDimension("real_or_complex_gw_corrections", "integer", allowed=[1, 2])
real_or_complex_potential = EtsfDimension("real_or_complex_potential", "integer", allowed=[1, 2])
real_or_complex_wavefunctions = EtsfDimension("real_or_complex_wavefunctions", "integer", allowed=[1, 2])
number_of_cartesian_directions = EtsfDimension("number_of_cartesian_directions", "integer", allowed=[3])
number_of_reduced_dimensions = EtsfDimension("number_of_reduced_dimensions", "integer", [3])
number_of_vectors = EtsfDimension("number_of_vectors", "integer", allowed=[3])
number_of_symmetry_operations = EtsfDimension("number_of_symmetry_operations", "integer")
number_of_atoms	= EtsfDimension("number_of_atoms", "integer")
number_of_atom_species = EtsfDimension("number_of_atom_species", "integer")
symbol_length = EtsfDimension("symbol_length", "integer", allowed=[2])

##############################################################################
# Dimensions that can be split
# This table list the dimensions that might be used to define a splitting 
# (e.g. in case of parallelism). 
##############################################################################
max_number_of_states = EtsfDimension("max_number_of_states", "integer")
number_of_kpoints = EtsfDimension("number_of_kpoints", "integer")
number_of_spins = EtsfDimension("number_of_spins", "integer", allowed=[1, 2])
number_of_spinor_components = EtsfDimension("number_of_spinor_components", "integer", allowed=[1, 2])
number_of_components = EtsfDimension("number_of_components", "integer", allowed=[1, 2, 4])
max_number_of_coefficients = EtsfDimension("max_number_of_coefficients", "integer")
number_of_grid_points_vector1 = EtsfDimension("number_of_grid_points_vector1", "integer")
number_of_grid_points_vector2 = EtsfDimension("number_of_grid_points_vector2", "integer")
number_of_grid_points_vector3 = EtsfDimension("number_of_grid_points_vector3", "integer")
max_number_of_basis_grid_points	= EtsfDimension("max_number_of_basis_grid_points", "integer")
number_of_localisation_regions = EtsfDimension("number_of_localisation_regions", "integer", allowed=[1])

####################
# Atomic information
####################
valence_charges = EtsfVariable("valence_charges", "double",  [number_of_atom_species])
pseudopotential_types = EtsfVariable("pseudopotential_types", "char", [number_of_atom_species, character_string_length])

######################
# Electronic structure
######################
number_of_electrons = EtsfVariable("number_of_electrons", "integer", [])	
exchange_functional = EtsfVariable("exchange_functional", "char", [character_string_length])
correlation_functional = EtsfVariable("correlation_functional", "char", [character_string_length])
fermi_energy = VariableWithUnits("fermi_energy", "double", [], reqattrs=[units])
# Units attribute required. The attribute scale to atomic units might also be mandatory
smearing_scheme	= EtsfVariable("smearing_scheme", "char", [character_string_length])
smearing_width = VariableWithUnits("smearing_width", "double", [], reqattrs=[units])
# Units attribute required. The attribute scale to atomic units might also be mandatory

##################
# Reciprocal space
##################
kinetic_energy_cutoff = VariableWithUnits("kinetic_energy_cutoff", "double", [], reqattrs=[units])
# Units attribute required. The attribute scale to atomic units might also be mandatory
kpoint_grid_shift = EtsfVariable("kpoint_grid_shift", "double", [number_of_reduced_dimensions])
kpoint_grid_vectors = EtsfVariable("kpoint_grid_vectors", "double", [number_of_vectors, number_of_reduced_dimensions])
monkhorst_pack_folding = EtsfVariable("monkhorst_pack_folding", "integer", [number_of_vectors])

###################################################################################
# Atomic structure and symmetry operations
# Variables and attributes to specify the atomic structure and symmetry operations.
###################################################################################
primitive_vectors = EtsfVariable("primitive_vectors", "double", [number_of_vectors, number_of_cartesian_directions])

reduced_symmetry_matrices = EtsfVariable("reduced_symmetry_matrices", "integer", 
  [number_of_symmetry_operations, number_of_reduced_dimensions, number_of_reduced_dimensions],
  reqattrs=[symmorphic]) # The "symmorphic" attribute is needed.

reduced_symmetry_translations = EtsfVariable("reduced_symmetry_translations", "double", 
    [number_of_symmetry_operations, number_of_reduced_dimensions])
# The "symmorphic" attribute is needed.

# In principle: allowed=range(1, 233)) but I usually use 0 when the space_group is not avaiable
space_group = EtsfVariable("space_group", "integer", [], allowed=range(0, 233)) 
atom_species = EtsfVariable("atom_species", "integer", [number_of_atoms]) # Between 1 and number_of_atom_species.

reduced_atom_positions = EtsfVariable("reduced_atom_positions", "double", [number_of_atoms, number_of_reduced_dimensions])
atomic_numbers = EtsfVariable("atomic_numbers", "double", [number_of_atom_species])
atom_species_names = EtsfVariable("atom_species_names", "char", [number_of_atom_species, character_string_length])
chemical_symbols = EtsfVariable("chemical_symbol", "char", [number_of_atom_species, symbol_length])

##########
# K-points
##########
reduced_coordinates_of_kpoints = EtsfVariable("reduced_coordinates_of_kpoints", "double", 
    [number_of_kpoints, number_of_reduced_dimensions])
kpoint_weights = EtsfVariable("kpoint_weights", "double", [number_of_kpoints])

########
# States
########
number_of_states = EtsfVariable("number_of_states", "integer", [number_of_spins, number_of_kpoints], 
    reqattrs=[k_dependent])  # The attribute "k_dependent" must be defined.

eigenvalues = VariableWithUnits("eigenvalues", "double", [number_of_spins, number_of_kpoints, max_number_of_states])
occupations = EtsfVariable("occupations", "double", [number_of_spins, number_of_kpoints, max_number_of_states])

# Density
# A density in such a format (represented on a 3D homogeneous grid) is suited for the representation 
# of smooth densities, as obtained naturally from pseudopotential calculations using plane waves.
# This specification for a density can also accommodate the response densities of Density-Functional Perturbation Theory.

density = VariableWithUnits("density", "double", [number_of_components, number_of_grid_points_vector3, 
            number_of_grid_points_vector2, number_of_grid_points_vector1, real_or_complex_density],
            reqattrs=[units])
# By default, the density is given in atomic units, that is, number of electrons per Bohr3. 
# The "units" attribute is required. The attribute "scale_to_atomic_units" might also be mandatory

# Exchange and correlation
correlation_potential = VariableWithUnits("correlation_potential", "double", [number_of_components, 
  number_of_grid_points_vector3, number_of_grid_points_vector2, number_of_grid_points_vector1, real_or_complex_potential],
  reqattrs=[units])
#Units attribute required. The attribute scale to atomic units might also be mandatory

exchange_potential = VariableWithUnits("exchange_potential", "double", [number_of_components, 
  number_of_grid_points_vector3, number_of_grid_points_vector2, number_of_grid_points_vector1, real_or_complex_potential],
  reqattrs=[units])
# Units attribute required. The attribute scale to atomic units might also be mandatory

exchange_correlation_potential = VariableWithUnits("exchange_correlation_potential", "double", [number_of_components, 
  number_of_grid_points_vector3, number_of_grid_points_vector2, number_of_grid_points_vector1, real_or_complex_potential],
  reqattrs=[units])
# Units attribute required. The attribute "scale to atomic units" might also be mandatory


class EtsfGroup(object):
    """"
    This object is essentially a container of variables
    A Group can contain other subgroups.
    """
    attributes = []
    dimensions = []
    variables = []
    at_least_one_variable_in = []
    subgroups = []

    @classmethod
    def validate_file(cls, path):
        ncdata = netCDF4.Dataset(path, mode="r")
        errors = cls.validate(ncdata)
        ncdata.close()
        return errors

    @classmethod
    def validate(cls, ncdata):
        errors = []

        # Test attributes.
        for attr in cls.attributes:
            errors += attr.validate(ncdata)

        # Test dimensions
        for dim in cls.dimensions:
            errors += dim.validate(ncdata)

        # Test variables
        for var in cls.variables:
            errors += var.validate(ncdata)

        # Test variables.
        ok, elist = 0, []
        for var in cls.at_least_one_variable_in:
            if var.name in ncdata.variables:
                ok += 1
            else:
                elist.append("<%s> not present in ncdata" % var.name)

        if not ok:
            errors.extend(elist)
                                                                                                
        # Test subgroups.
        for group in cls.subgroups:
            errors += group.validate(ncdata)

        return errors


class CrystalGroup(EtsfGroup):
    """
    A ETSF NetCDF file for crystallographic data should contain the following set of mandatory information:

    The three attributes defined in Mandatory attributes

    The following dimensions from Dimensions that cannot be split:

        number_of_cartesian_directions
        number_of_vectors
        number_of_atoms
        number_of_atom_species
        number_of_symmetry_operations

    The following variables defined in Atomic structure and symmetry operations:

        primitive_vectors
        reduced_symmetry_matrices
        reduced_symmetry_translations
        space_group
        atom_species
        reduced_atom_positions

    At least one of the following variables defined in Atomic structure and symmetry operations, 
    to specify the kind of atoms:

        atomic_numbers
        atom_species_names
        chemical_symbols

    The use of atomic_numbers is preferred. If atomic_numbers is not available, atom_species_names 
    will be preferred over chemical_symbols. In case more than one such variables are present in a file, 
    the same order of preference should be followed by the reading program.
    """
    attributes = mandatory_attributes

    dimensions = [
        number_of_cartesian_directions,
        number_of_vectors,
        number_of_atoms,
        number_of_atom_species,
        number_of_symmetry_operations,
    ]

    variables = [
        primitive_vectors,
        reduced_symmetry_matrices,
        reduced_symmetry_translations,
        space_group,
        atom_species,
        reduced_atom_positions,
    ]

    at_least_one_variable_in = [
        atomic_numbers,
        atom_species_names,
        chemical_symbols,
    ]


class KpointsGroup(EtsfGroup):
    variables = [
        reduced_coordinates_of_kpoints,
        kpoint_weights,
    ]

class StatesGroup(EtsfGroup):
    variables = [
        number_of_states,
        eigenvalues,
        occupations,
    ]


class DenPotGroup(EtsfGroup):
    """
    A ETSF NetCDF file for a density should contain the following set of mandatory information:

    The three attributes defined in Mandatory attributes

    The following dimensions from Dimensions that cannot be split:

        number_of_cartesian_directions
        number_of_vectors
        real_or_complex_density and/or real_or_complex_potential

    The following dimensions from Dimensions that can be split:

        number_of_components
        number_of_grid_points_vector1
        number_of_grid_points_vector2
        number_of_grid_points_vector3

    The primitive vectors of the cell, as defined in Atomic structure and symmetry operations.

    The density or potential, as defined in Density or Exchange and correlation. 
    This variable must be the last, in order not to be limited to 4 GB.

    As mentioned in General considerations and General specifications, such file might contain 
    additional information agreed within ETSF, such as any of the variables specified in General specifications. 
    It might even contain enough information to be declared a ETSF NetCDF file "containing crystallographic data or 
    "containing the wavefunctions", or both. Such file might also contain additional information specific 
    to the software that generated the file. It is not expected that this other software-specific information 
    be used by another software.

    A ETSF NetCDF exchange, correlation, or exchange-correlation potential file should contain at least 
    one variable among the three presented in Exchange and correlation in replacement of the specification 
    of the density. The type and size of such variables are similar to the one of the density. 
    The other variables required for a density are also required for a potential file. 
    Additional ETSF or software-specific information might be added, as described previously.
    The information might distributed among different files, thanks to the use of splitting of data for variables:

        number_of_components
        number_of_grid_points_vector1
        number_of_grid_points_vector2
        number_of_grid_points_vector3

    In case the splitting related to one of these variables is activated, then the corresponding variables 
    in Auxiliary variables for splitting must be defined. Accordingly, the dimensions of the variables 
    in Density and/or Exchange and correlation will be changed, to accommodate only the segment of data 
    effectively contained in the file.
    """
    attributes = mandatory_attributes

    dimensions = [
        number_of_cartesian_directions,
        number_of_vectors,
        #real_or_complex_density and/or real_or_complex_potential,
        number_of_components,
        number_of_grid_points_vector1,
        number_of_grid_points_vector2,
        number_of_grid_points_vector3,
        # Atomic structure and symmetry operations
    ]

    variables = [
        density,
    ]

    subgroups = [
        CrystalGroup,
    ]


class WavefunctionGroup(EtsfGroup):
    """
    Specification for files containing the wavefunctions

    Specification

    A ETSF NetCDF file "containing the wavefunctions" should contain at least the information needed to build 
    the density from this file. Also, since the eigenvalues are intimately linked to eigenfunctions, 
    it is expected that such a file contain eigenvalues. Of course, files might contain less information 
    than the one required, but still follow the naming convention of ETSF. It might also contain more information, 
    of the kind specified in other tables of the present document.

    A ETSF NetCDF file "containing the wavefunctions" should contain the following set of mandatory information:

    The three attributes defined in Mandatory attributes

    The following dimensions from Dimensions that cannot be split:

        character_string_length
        number_of_cartesian_directions
        number_of_vectors
        real_or_complex_coefficients and/or real_or_complex_wavefunctions
        number_of_symmetry_operations
        number_of_reduced_dimensions

    The following dimensions from Dimensions that can be split:

        max_number_of_states
        number_of_kpoints
        number_of_spins
        number_of_spinor_components

    In case of a real-space wavefunctions, the following dimensions from Dimensions that can be split:

        number_of_grid_points_vector1
        number_of_grid_points_vector2
        number_of_grid_points_vector3

    In case of a wavefunction given in terms of a basis set, the following dimensions from Dimensions that can be split:

        max_number_of_coefficients

    In case of a wavefunction given in terms of a Daubechies wavelet basis set, the following dimensions 
    from Dimensions that can be split:

        max_number_of_basis_grid_points
        number_of_localization_regions

    The primitive vectors of the cell, as defined in Atomic structure and 
    symmetry operations (variable primitive_vectors)

    The symmetry operations, as defined in Atomic structure and symmetry operations 
    (given by the variables reduced_symmetry_translations and reduced_symmetry_matrices)

    The information related to each kpoint, as defined in K-points.

    The information related to each state (including eigenenergies and occupation numbers), as defined in States.

    In case of basis set representation, the information related to the basis set, and the variable 
    coefficients_of_wavefunctions, as defined in Wavefunctions.

    For basis_set equal to "plane_waves", the following variable is required from Wavefunctions:

        reduced_coordinates_of_plane_waves

    For basis_set equal "daubechies_wavelets", the following variables are required from Wavefunctions:

        coordinates_of_basis_grid_points
        number_of_coefficients_per_grid_point

    In case of real-space representation, the variable following variable from Wavefunctions:

        real_space_wavefunctions

    As mentioned in General considerations and General specifications, such a file might contain 
    additional information agreed on within ETSF, such as any of the variables specified in General specifications. 
    It might even contain enough information to be declared a ETSF NetCDF file "containing crystallographic data" 
    or "containing the density", or both. Such a file might also contain additional information specific 
    to the software that generated the file. It is not expected that this other software-specific information 
    be used by another software.

    The information might be distributed among different files, thanks to the use of splitting of data for variables:

        max_number_of_states
        number_of_kpoints
        number_of_spins

    And, either

        number_of_grid_points_vector1
        number_of_grid_points_vector2
        number_of_grid_points_vector3

    or
        max_number_of_coefficients

    In case the splitting related to one of these variables is activated, then the corresponding variables 
    in Split wavefunctions must be defined. Accordingly, the dimensions of the variables in K-points, States, 
    Wavefunctions, and BSE/GW might have to be changed, to accommodate only the segment of data effectively 
    contained in the file.
    """
    dimensions = [
        character_string_length,
        number_of_cartesian_directions,
        number_of_vectors,
        #real_or_complex_coefficients and/or real_or_complex_wavefunctions
        number_of_symmetry_operations,
        number_of_reduced_dimensions,
        #
        max_number_of_states,
        number_of_kpoints,
        number_of_spins,
        number_of_spinor_components,
        #
        number_of_grid_points_vector1,
        number_of_grid_points_vector2,
        number_of_grid_points_vector3,

        #coefficients_of_wavefunctions,
        #reduced_coordinates_of_plane_waves
    ]

    subgroups = [
        CrystalGroup,
        KpointsGroup,
        StatesGroup,
    ]

class DensityGroup(DenPotGroup):
    pass

class PotentialGroup(DenPotGroup):
    pass


def validate_vars(path):
    """
    Validate the etsf variables declared in file `path`. 
    Return list of errors.
    """
    ncdata = netCDF4.Dataset(path, mode="r")
    evars, all_errors = [] , []
    d = {}
    for var in EtsfVariable.all_variables:
       if var.name in ncdata.variables:
            elist = var.validate(ncdata)
            if elist:
                all_errors.append(elist)
                evars.append(var)
    ncdata.close()

    if not all_errors: return all_errors

    # Print errors
    try:
        title = "=== ERRORS IN FILE %s ===" % os.path.relpath(path)
    except Exception as exc:
        # I've seen this kind of error.
        logger.warning("relpath returned %s" % str(exc))
        title = "=== ERRORS IN FILE %s ===" % path

    bar = "=" * len(title)
    cprint(bar, "magenta")
    cprint(title, "magenta")
    cprint(bar, "magenta")

    for var, elist in zip(evars, all_errors):
        cprint("%s: Found [%d] error(s):" % (var, len(elist)), "red")
        for i, e in enumerate(elist): 
            print("[%d] " % i, e)

    return all_errors


def find_groups(path):
    """Find the etsf groups present in path and validate them."""
    ncdata = netCDF4.Dataset(path, mode="r")
    groups = []
    for g in all_subclasses(EtsfGroup):
        if g.validate(ncdata):
            groups.append(g)
            print("group:", g, "present and validated")
        else:
            print("group:", g, "present but not ok")
    ncdata.close()
    return groups


def validate_groups(path, groups):
    """
    Validate the presence and the consistency of a list of groups..
    Return list of errors
    """
    ncdata = netCDF4.Dataset(path, mode="r")

    wrong_datamodel = ""
    # TODO: Activate this test.
    #if ncdata.data_model != "NETCDF4":
    #    wrong_datamodel = "Found data_model %s while it should be NETCDF4" % ncdata.data_model

    egroups, errors = [], []
    for group in groups:
        elist = group.validate(ncdata)
        if elist: 
            egroups.append(group)
            errors.append(elist)
    ncdata.close()

    if wrong_datamodel or errors:
        # Print errors
        try:
            title = "=== ERRORS IN FILE %s ===" % os.path.relpath(path)
        except Exception as exc:
            # I've seen this kind of error.
            logger.warning("relpath returned %s" % str(exc))
            title = "=== ERRORS IN FILE %s ===" % path
                                                                        
        bar = "=" * len(title)
        cprint(bar, "magenta")
        cprint(title, "magenta")
        cprint(bar, "magenta")

        for group, elist in zip(egroups, errors):
            print(group)
            for i, e in enumerate(elist): 
                print("[%d]" % i, e)

        if wrong_datamodel: 
            errors.append(wrong_datamodel)
            print(wrong_datamodel)

    return errors


def validate_ncfile(path):
    """
    This function validates the netcdf files produced by Abinit.

    Return:
        List of strings with error messages.
    """
    # Every netcdf file produced by Abinit must be listed in this dictionary.
    # that maps the file extension to the list of Groups contained in the file.
    # Legacy files that do not have etsfio groups are tagged with None and ignored.
    # New netcdf files introduced in Abinit should try to use the etsfio specs as much 
    # as possible and extend the format to treat new cases.
    # Remember that gmatteo is watching you so don't try to hack this code.
    ext2groups = {
        # GS files.
        "GSR.nc": [CrystalGroup, KpointsGroup, StatesGroup],
        "DEN.nc": [DensityGroup],
        "EIG.nc": None,
        "OUT.nc": None,
        "WFK.nc": [WavefunctionGroup],
        #"POT.nc": [PotentialGroup],
        #"VH.nc": [PotentialGroup],
        #"VXC.nc": [PotentialGroup],
        "VHXC.nc": [PotentialGroup],
        #"KDEN.nc": [DensityGroup],
        "HIST.nc": None,
        #"PAWDEN.nc": [DensityGroup],
        #"ELF.nc": [DensityGroup],
        #"ELF_DOWN.nc": [DensityGroup],
        #"ELF_UP.nc": [DensityGroup],
        "PSPS.nc": None,
        #"STM.nc": [],
        #"GDEN1.nc":
        #"GDEN2.nc":
        #"GDEN3.nc":
        #"KDEN.nc":
        #"LDEN.nc":

        # DFPT files.
        "DDB.nc": [CrystalGroup],
        "WFQ.nc": [WavefunctionGroup],
        #"DEN1.nc": [DensityGroup],
        #"POT1.nc": [PotentialGroup],
        #"EIGR2D.nc": 
        #"EIGI2D.nc":

        # Anaddb files
        "anaddb.nc": [CrystalGroup],
        "EC.nc": [CrystalGroup],
        "PHBST.nc": [CrystalGroup],

        # GW files.
        #"KSS.nc": [WavefunctionGroup],
        "SCR.nc": [CrystalGroup, KpointsGroup, StatesGroup],
        "SUS.nc": [CrystalGroup, KpointsGroup, StatesGroup],
        "SIGRES.nc": [CrystalGroup, KpointsGroup, StatesGroup],
        #"QP_DEN.nc": [DensityGroup],

        # BSE files.
        "MDF.nc": [CrystalGroup, KpointsGroup, StatesGroup],
    }

    groups = None
    fname = os.path.basename(path)

    # DFPT files require regular expressions e.g. 
    # _DS2_1WF6.nc, _DS2_DEN6, _DEN_POSITRON
    re_1wf = re.compile("(\w+_)1WF(\d+)(\.nc)$")
    re_1den = re.compile("(\w+_)DEN(\d+)(\.nc)$")
    #re_1den = re.compile("(\w+_)POT(\d+)(\.nc)$")
    if re_1wf.match(fname): groups = [WavefunctionGroup]
    if re_1den.match(fname): groups = [DensityGroup]
    if groups is not None: return validate_groups(path, groups)

    ext = fname.split("_")[-1]
    try:
        groups = ext2groups[ext]
    except KeyError as exc:
        errors = ["Unknown file extension in file %s" % fname]
        print(errors)
        return errors

    if groups is None: return []
    return validate_groups(path, groups)
