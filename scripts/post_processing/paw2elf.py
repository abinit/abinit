# SPDX-License-Identifier: GPL-3.0-only
#
# Copyright (C) 2026 Josef Zwanziger
# jzwanzig@dal.ca
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# Portions of this program were initially drafted with assistance from
# ChatGPT (OpenAI). The scientific design, integration, review,
# modification, validation, and responsibility for the final code are
# those of the human author.

# TL;DR
# paw2elf.py generates an ELF netcdf file from PAW density and kinetic energy density
# netcdf files.

#  ===================================================================================
#  OVERVIEW of paw2elf.py script
#  ===================================================================================
#  Currently Abinit does not generate the electron localization function, ELF,
#  during a PAW run, as it does for norm-conserving pseudopotentials. This program
#  constructs ELF as a post-processing routine, using the computed density and
#  kinetic energy density files.
#
#  The preparatory Abinit calculation must include output of the DEN and KDEN files,
#  that is, density and kinetic energy density. Therefore you must include the 
#  keywords "prtden 1" or "prtden 2", "prtkden 1", and "usekden 1". Also, the output
#  must be in netcdf form ("iomode 3"). Note that "prtden 1" generates the DEN.nc file,
#  while "prtden 2" generates both DEN.nc and PAWDEN.nc. The program can use either 
#  one, and this is selected by input option -d. 
#
#  The program will also read the PAW data files used, in order to extract the core
#  density and kinetic energy density. It will determine where to find them based on
#  the content of the DEN.nc or PAWDEN.nc files. Therefore it is probably safest to
#  run this program in the same directory where you generated the DEN files; if you
#  move these files to some other computer the program probably won't be able to 
#  find the PAW files.
#
#  The output will be a new netcdf file, named DEN_ELF.nc or PAWDEN_ELF.nc, which will
#  include both the original DEN information as well as the new ELF information.
 
import xml.etree.ElementTree as ET
import argparse
import textwrap
import numpy as np
import math
from scipy.integrate import simpson
from scipy.interpolate import PchipInterpolator
from itertools import product
import netCDF4
import warnings
# Silence the specific NumPy 2.5 shape deprecation warning
warnings.filterwarnings("ignore", category=DeprecationWarning, message="Setting the shape on a NumPy array")
import os
import sys
import hashlib
import shutil

def write_elf_like_pawrhor(
        input_filename,
        output_filename,
        elf,
        paw_density_variable="density"):
    """
    Copy a PAWDEN NetCDF file and add ELF with the same storage ordering
    as pawrhor.

    elf is supplied in Python order (n1,n2,n3).
    """
    elf = np.asarray(elf, dtype=np.float64)

    shutil.copyfile(input_filename, output_filename)

    with netCDF4.Dataset(output_filename, mode="a") as nc:
        pawrhor = nc.variables[paw_density_variable]

        #print("pawrhor dimensions:", pawrhor.dimensions)
        #print("pawrhor shape:", pawrhor.shape)
        #print("input ELF shape:", elf.shape)

        if pawrhor.ndim != 5:
            raise ValueError(
                "Expected pawrhor to have five dimensions."
            )

        expected_spatial_shape = (
            pawrhor.shape[3],
            pawrhor.shape[2],
            pawrhor.shape[1]
        )

        if elf.shape != expected_spatial_shape:
            raise ValueError(
                f"ELF has shape {elf.shape}, but the expected Python "
                f"spatial shape is {expected_spatial_shape}."
            )

        # Convert:
        #
        # (n1,n2,n3) -> (component,n3,n2,n1,real_or_complex)
        elf_stored = elf.transpose(2, 1, 0)[None, :, :, :, None]

        if elf_stored.shape != pawrhor.shape:
            raise ValueError(
                f"Converted ELF shape {elf_stored.shape} does not "
                f"match pawrhor shape {pawrhor.shape}."
            )

        if "pawelf" in nc.variables:
            elf_var = nc.variables["pawelf"]

            if elf_var.shape != elf_stored.shape:
                raise ValueError(
                    "An ELF variable exists but has an incompatible shape."
                )
        else:
            elf_var = nc.createVariable(
                "pawelf",
                "f8",
                pawrhor.dimensions,
                zlib=True,
                complevel=4,
                shuffle=True
            )

        elf_var[:] = elf_stored

        elf_var.long_name = "Electron localization function"
        elf_var.units = "dimensionless"
        elf_var.valid_min = 0.0
        elf_var.valid_max = 1.0
        elf_var.spatial_axis_order = (
            "Stored as component,n3,n2,n1,real_or_complex"
        )

    #print(f"Wrote {output_filename}")

def calculate_md5(file_path):
    md5_hash = hashlib.md5()

    # Open the file in read-binary mode
    with open(file_path, "rb") as f:
        # Read the file in 4KB chunks
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)

    return md5_hash.hexdigest()

def minimum_image_vector(delta_frac, lattice, image_range=1):
    """
    Find the shortest periodic displacement vector in Cartesian coordinates.

    Parameters
    ----------
    delta_frac : ndarray, shape (..., 3)
        Fractional displacement grid_position - atom_position.
    lattice : ndarray, shape (3, 3)
        Direct-lattice vectors stored as rows.
    image_range : int
        Number of periodic images searched in each direction.

    Returns
    -------
    displacement_cart : ndarray, shape (..., 3)
        Shortest Cartesian displacement vector.
    distance : ndarray, shape (...)
        Length of the displacement vector.
    """
    delta0 = delta_frac - np.rint(delta_frac)

    translations = np.array(
        list(product(
            range(-image_range, image_range + 1),
            repeat=3
        )),
        dtype=float
    )

    best_r2 = np.full(delta0.shape[:-1], np.inf)
    best_vector = np.zeros_like(delta0, dtype=float)

    for translation in translations:
        candidate_frac = delta0 + translation
        candidate_cart = candidate_frac @ lattice

        candidate_r2 = np.einsum(
            "...i,...i->...",
            candidate_cart,
            candidate_cart
        )

        replace = candidate_r2 < best_r2

        best_r2[replace] = candidate_r2[replace]
        best_vector[replace] = candidate_cart[replace]

    return best_vector, np.sqrt(best_r2)

def core_density_gradient(
        grid_shape,
        lattice,
        atoms,
        chunk_size=200_000,
        image_range=1):
    """
    Construct the Cartesian gradient of atom-centered spherical core densities.

    Parameters
    ----------
    grid_shape : tuple
        Shape (n1,n2,n3) of the real-space grid.

    lattice : ndarray, shape (3,3)
        Direct-lattice vectors stored as rows.

    atoms : sequence of dictionaries
        Each dictionary contains:

            "frac_position" : atomic fractional coordinates
            "r"             : radial grid
            "core_density"  : spherical 3D core density n_c(r)

    Returns
    -------
    grad_core : ndarray, shape (3,n1,n2,n3)
        Cartesian components of the core-density gradient.
    """
    lattice = np.asarray(lattice, dtype=float)

    n1, n2, n3 = grid_shape
    npoints = n1 * n2 * n3

    prepared_atoms = []

    for atom in atoms:
        r = np.asarray(atom["r"], dtype=float)
        rho = np.asarray(atom["core_density"], dtype=float)

        order = np.argsort(r)
        r = r[order]
        rho = rho[order]

        r, unique_indices = np.unique(r, return_index=True)
        rho = rho[unique_indices]

        rho_interpolator = PchipInterpolator(
            r,
            rho,
            extrapolate=False
        )

        prepared_atoms.append({
            "frac_position": np.asarray(
                atom["frac_position"],
                dtype=float
            ),
            "drho_dr": rho_interpolator.derivative(),
            "rmin": float(r[0]),
            "rmax": float(r[-1]),
        })

    # Stored as (Cartesian component, flattened grid index).
    grad_core_flat = np.zeros((3, npoints), dtype=float)

    for start in range(0, npoints, chunk_size):
        stop = min(start + chunk_size, npoints)
        flat_indices = np.arange(start, stop)

        i, j, k = np.unravel_index(
            flat_indices,
            grid_shape,
            order="C"
        )

        grid_frac = np.column_stack((
            i / n1,
            j / n2,
            k / n3
        ))

        grad_block = np.zeros((stop - start, 3), dtype=float)

        for atom in prepared_atoms:
            delta_frac = grid_frac - atom["frac_position"]

            displacement, radius = minimum_image_vector(
                delta_frac,
                lattice,
                image_range=image_range
            )

            # Avoid division by zero at the nucleus.
            inside = (
                (radius >= atom["rmin"])
                & (radius <= atom["rmax"])
                & (radius > 1.0e-14)
            )

            radial_derivative = atom["drho_dr"](
                radius[inside]
            )

            grad_block[inside] += (
                radial_derivative[:, None]
                * displacement[inside]
                / radius[inside, None]
            )

            # At r=0, the gradient of a smooth spherical density is zero.

        grad_core_flat[:, start:stop] = grad_block.T

    return grad_core_flat.reshape(
        (3, n1, n2, n3),
        order="C"
    )

def density_gradient_magnitude(rho_ae, lattice, return_components=False):
    """
    Calculate |grad rho| on a periodic triclinic real-space grid using
    spectral FFT differentiation.

    Parameters
    ----------
    rho_ae : ndarray, shape (n1, n2, n3)
        Real-space electron density on a uniform fractional-coordinate grid.

    lattice : ndarray, shape (3, 3)
        Direct-lattice vectors stored as rows:

            lattice[0, :] = a1
            lattice[1, :] = a2
            lattice[2, :] = a3

        If the lattice vectors are in bohr and rho is in electrons/bohr^3,
        the returned gradient is in electrons/bohr^4.

    return_components : bool, optional
        If True, also return the Cartesian gradient components.

    Returns
    -------
    grad_magnitude : ndarray, shape (n1, n2, n3)
        Magnitude |grad rho|.

    grad_cart : ndarray, shape (3, n1, n2, n3), optional
        Cartesian components of the gradient. Returned only when
        return_components=True. The first index corresponds to x, y, z.
    """
    rho_ae = np.asarray(rho_ae, dtype=float)
    lattice = np.asarray(lattice, dtype=float)

    if rho_ae.ndim != 3:
        raise ValueError("rho_ae must be a three-dimensional array.")

    if lattice.shape != (3, 3):
        raise ValueError("lattice must have shape (3, 3).")

    n1, n2, n3 = rho_ae.shape

    # Integer Fourier indices associated with the fractional coordinates.
    #
    # np.fft.fftfreq(N, d=1/N) gives:
    #   0, 1, ..., N/2-1, -N/2, ..., -1
    m1 = np.fft.fftfreq(n1, d=1.0 / n1)
    m2 = np.fft.fftfreq(n2, d=1.0 / n2)
    m3 = np.fft.fftfreq(n3, d=1.0 / n3)

    M1, M2, M3 = np.meshgrid(
        m1, m2, m3,
        indexing="ij",
        sparse=True
    )

    # If r = s @ lattice, a Fourier mode is
    #
    #   exp(i G.r) = exp(i 2*pi*m.s),
    #
    # and therefore
    #
    #   G = 2*pi*m @ lattice^{-T}.
    #
    inverse_lattice_transpose = np.linalg.inv(lattice).T

    gx = 2.0 * np.pi * (
        M1 * inverse_lattice_transpose[0, 0]
        + M2 * inverse_lattice_transpose[1, 0]
        + M3 * inverse_lattice_transpose[2, 0]
    )

    gy = 2.0 * np.pi * (
        M1 * inverse_lattice_transpose[0, 1]
        + M2 * inverse_lattice_transpose[1, 1]
        + M3 * inverse_lattice_transpose[2, 1]
    )

    gz = 2.0 * np.pi * (
        M1 * inverse_lattice_transpose[0, 2]
        + M2 * inverse_lattice_transpose[1, 2]
        + M3 * inverse_lattice_transpose[2, 2]
    )

    rho_fourier = np.fft.fftn(rho_ae)

    grad_x = np.fft.ifftn(1j * gx * rho_fourier).real
    grad_y = np.fft.ifftn(1j * gy * rho_fourier).real
    grad_z = np.fft.ifftn(1j * gz * rho_fourier).real

    grad_magnitude = np.sqrt(
        grad_x**2 + grad_y**2 + grad_z**2
    )

    if return_components:
        grad_cart = np.stack((grad_x, grad_y, grad_z), axis=0)
        return grad_magnitude, grad_cart

    return grad_magnitude

def radial_interpolator(r, f):
    """
    Construct a shape-preserving interpolation of a radial function.

    Values outside the tabulated radial range are set to zero.
    """
    r = np.asarray(r, dtype=float)
    f = np.asarray(f, dtype=float)

    order = np.argsort(r)
    r = r[order]
    f = f[order]

    # Remove duplicate radial points, should they occur.
    r, unique_indices = np.unique(r, return_index=True)
    f = f[unique_indices]

    return PchipInterpolator(
        r,
        f,
        extrapolate=False
    )

def triclinic_minimum_image_distance(delta_frac, lattice,
                                      image_range=1):
    """
    Return the shortest Cartesian distance between periodically
    equivalent points in a triclinic cell.

    Parameters
    ----------
    delta_frac : ndarray, shape (..., 3)
        Fractional displacement grid_position - atom_position.
    lattice : ndarray, shape (3, 3)
        Lattice vectors stored as rows.
    image_range : int
        Number of neighboring images tested around the fractionally
        wrapped image. image_range=1 tests 27 images.

    Returns
    -------
    rmin : ndarray, shape delta_frac.shape[:-1]
        Minimum Cartesian distance.
    """
    lattice = np.asarray(lattice, dtype=float)

    # First bring the fractional displacement near the central cell.
    delta0 = delta_frac - np.rint(delta_frac)

    translations = np.array(
        list(product(
            range(-image_range, image_range + 1),
            repeat=3
        )),
        dtype=float
    )

    r2_min = np.full(delta0.shape[:-1], np.inf)

    for translation in translations:
        candidate_frac = delta0 + translation

        # Row-vector convention:
        # r_cart = s_frac @ lattice
        candidate_cart = candidate_frac @ lattice

        r2 = np.einsum(
            "...i,...i->...",
            candidate_cart,
            candidate_cart
        )

        r2_min = np.minimum(r2_min, r2)

    return np.sqrt(r2_min)

def add_spherical_core_fields(
        valence_density,
        valence_tau,
        lattice,
        atoms,
        chunk_size=200_000,
        image_range=1):
    """
    Add atom-centered spherical core density and kinetic-energy density
    to fields defined on a periodic triclinic FFT grid.

    Parameters
    ----------
    valence_density : ndarray, shape (n1, n2, n3)
        PAW valence electron density.
    valence_tau : ndarray, shape (n1, n2, n3)
        PAW valence kinetic-energy density.
    lattice : ndarray, shape (3, 3)
        Lattice vectors in Cartesian coordinates, stored as rows.
    atoms : sequence of dictionaries
        Each dictionary must contain

            "frac_position" : array-like, length 3
            "r"             : radial grid
            "core_density"  : core density on radial grid
            "core_tau"      : core kinetic-energy density on radial grid

    chunk_size : int
        Number of real-space grid points processed at once.
    image_range : int
        Periodic images tested in each lattice direction.

    Returns
    -------
    total_density, total_tau : ndarrays
        Valence plus core fields.
    """
    valence_density = np.asarray(valence_density, dtype=float)
    valence_tau = np.asarray(valence_tau, dtype=float)

    if valence_density.shape != valence_tau.shape:
        raise ValueError("Density and tau grids must have the same shape.")

    shape = valence_density.shape
    n1, n2, n3 = shape
    npoints = n1 * n2 * n3

    core_density = np.zeros(npoints, dtype=float)
    core_tau = np.zeros(npoints, dtype=float)

    # Prepare the radial interpolators once per atom.
    prepared_atoms = []

    for atom in atoms:
        prepared_atoms.append({
            "frac_position":
                np.asarray(atom["frac_position"], dtype=float),

            "density_interp":
                radial_interpolator(
                    atom["r"],
                    atom["core_density"]
                ),

            "tau_interp":
                radial_interpolator(
                    atom["r"],
                    atom["core_tau"]
                ),

            "rmax":
                float(np.max(atom["r"]))
        })

    # Process flattened FFT-grid indices in blocks to limit memory use.
    for start in range(0, npoints, chunk_size):
        stop = min(start + chunk_size, npoints)
        flat_index = np.arange(start, stop)

        # C-order flattening: k is the fastest-running index.
        i, j, k = np.unravel_index(flat_index, shape, order="C")

        grid_frac = np.column_stack((
            i / n1,
            j / n2,
            k / n3
        ))

        rho_block = np.zeros(stop - start)
        tau_block = np.zeros(stop - start)

        for atom in prepared_atoms:
            delta_frac = grid_frac - atom["frac_position"]

            r = triclinic_minimum_image_distance(
                delta_frac,
                lattice,
                image_range=image_range
            )

            inside = r <= atom["rmax"]

            # PchipInterpolator returns NaN outside its domain because
            # extrapolate=False, so evaluate only inside the radial grid.
            rho_block[inside] += atom["density_interp"](r[inside])
            tau_block[inside] += atom["tau_interp"](r[inside])

        core_density[start:stop] = rho_block
        core_tau[start:stop] = tau_block

    core_density = core_density.reshape(shape, order="C")
    core_tau = core_tau.reshape(shape, order="C")

    total_density = valence_density + core_density
    total_tau = valence_tau + core_tau

    return total_density, total_tau

def get_core_densities_from_atompaw_file(pawfile):

    # evidently atompaw works with densities as n_rad = 4*Pi*r^2*n rather than n directly
    # sqrt(4*Pi) probably comes from implicit use of real spherical harmonics which would
    # occur in all the wavefunction routines
    # ae_core from paw xml file contains sqrt(4*Pi)*n_rad(r)/(4*Pi*r^2)
    # tau from paw xml file contains 1/2 * sqrt(4*Pi)*coretau(r)/(4*Pi*r^2) = sqrt(4*pi)*tau
    sphharm=math.sqrt(4.0*math.pi)
    
    try:
        tree = ET.parse(pawfile)
    except ET.ParseError as e:
        print(f"reading {pawfile}, XML parsing error: {e}")

    root = tree.getroot()

    grid = root.find('radial_grid').find('values').text
    grid_data = np.fromstring(grid,sep=' ')
    
    ae_core = root.find('ae_core_density').text
    ae_core_data = np.fromstring(ae_core,sep=' ')/sphharm

    kden_core = root.find('ae_core_kinetic_energy_density').text
    kden_core_data = np.fromstring(kden_core,sep=' ')/sphharm

    return grid_data,ae_core_data,kden_core_data

def make_atom_list(xred,typat,pseudofiles):
    atom_list=[]
    for iatom,atom in enumerate(xred):
        pseudofile=pseudofiles[typat[iatom]-1]
        grid,aecore,kdencore=get_core_densities_from_atompaw_file(pseudofile)

        atom_list.append({ "frac_position":xred[iatom],"r":grid, "core_density":aecore, "core_tau":kdencore })
    return atom_list

def red2cart(red,rprimd):
    return red[0]*rprimd[0]+red[1]*rprimd[1]+red[2]*rprimd[2]

def get_atom_positions(xred,rprimd):
    xcart=np.zeros(xred.shape)
    for iatom,atom in enumerate(xred):
        xcart[iatom]=atom[0]*rprimd[0]+atom[1]*rprimd[1]+atom[2]*rprimd[2]
       
    return xcart

def get_pseudos_from_nc(input_string,pawdenfile,md5_pseudos):
    for iword,word in enumerate(input_string):
        if 'pp_dirpath' in word:
            pp_dirpath=input_string[iword+1]
        if 'pseudos' in word:
            pseudos=input_string[iword+1]
    pp_dirpath = os.path.expandvars(pp_dirpath.strip('"'))
    pseudos = pseudos.strip('"').split(',')
    pseudo_list=[os.path.join(pp_dirpath, f) for f in pseudos]
    for ipseudo,pseudo_file in enumerate(pseudo_list):
        pseudo_md5=calculate_md5(pseudo_file)
        if not pseudo_md5 in md5_pseudos[ipseudo]:
            print("\n  PAW data file "+pseudo_file+" does not match md5 referenced in "+pawdenfile+", exiting \n")
            sys.exit()

    return pseudo_list

def get_elf(rho_ae,grad_rho_ae,tau_ae):
    
    # Becke and Edgecombe JCP 92, 5397 (1990)
    # Uniform electron gas at n(r), Eq. 13
    #UniformEG=(3.0/5.0)*(6.0*math.pi**2)**(2.0/3.0)*rho_ae**(5.0/3.0)

    # ELOC, Eq. 10
    #ELOC=tau_ae - (1.0/4.0)*grad_rho_ae**2/rho_ae

    # Savin et al. using total densities
    UniformEG=(3.0/10.0)*(3.0*math.pi**2)**(2.0/3.0)*rho_ae**(5.0/3.0)
    ELOC=tau_ae - (1.0/8.0)*grad_rho_ae**2/rho_ae
    
    return 1.0/(1.0+(ELOC/UniformEG)**2)

def get_ncfiles(filepath,dentype):
    # return a dictionary of field:{filename,fieldname}
    # for input

    dictout={}
    if dentype == "DEN":
        denfile=filepath+'_DEN.nc'
        if not os.path.exists(denfile):
            print("Desired DEN file "+denfile+" not found, exiting \n")
            sys.exit()
        dictout["density"]={"denfile":denfile,"denfield":"density"}
        elffile=filepath+'_DEN_ELF.nc'
        dictout["elf"]={"elffile":elffile,"elffield":"pawelf"}
    if dentype == "PAWDEN":
        pawdenfile=filepath+'_PAWDEN.nc'
        if not os.path.exists(pawdenfile):
            print("Desired PAWDEN file "+pawdenfile+" not found, exiting \n")
            sys.exit()
        dictout["density"]={"denfile":pawdenfile,"denfield":"pawrhor"}
        elffile=filepath+'_PAWDEN_ELF.nc'
        dictout["elf"]={"elffile":elffile,"elffield":"pawelf"}
    kdenfile=filepath+'_KDEN.nc'
    if not os.path.exists(kdenfile):
        print("Required KDEN file "+kdenfile+" not found, exiting \n")
        sys.exit()
    dictout["tau"]={"taufile":kdenfile,"taufield":"kinedr"}

    return dictout


long_instructions=textwrap.dedent('''
  ===================================================================================
  OVERVIEW of paw2elf.py script
  ===================================================================================
  Currently Abinit does not generate the electron localization function, ELF,
  during a PAW run, as it does for norm-conserving pseudopotentials. This program
  constructs ELF as a post-processing routine, using the computed density and
  kinetic energy density files.

  The preparatory Abinit calculation must include output of the DEN and KDEN files,
  that is, density and kinetic energy density. Therefore you must include the 
  keywords "prtden 1" or "prtden 2", "prtkden 1", and "usekden 1". Also, the output
  must be in netcdf form ("iomode 3"). Note that "prtden 1" generates the DEN.nc file,
  while "prtden 2" generates both DEN.nc and PAWDEN.nc. The program can use either 
  one, and this is selected by input option -d. 

  The program will also read the PAW data files used, in order to extract the core
  density and kinetic energy density. It will determine where to find them based on
  the content of the DEN.nc or PAWDEN.nc files. Therefore it is probably safest to
  run this program in the same directory where you generated the DEN files; if you
  move these files to some other computer the program probably won't be able to 
  find the PAW files.

  The output will be a new netcdf file, named DEN_ELF.nc or PAWDEN_ELF.nc, which will
  include both the original DEN information as well as the new ELF information.
        ''')

parser=argparse.ArgumentParser(
        description=long_instructions,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
parser.add_argument("-f","--filepath",help="file path, default gso, _DEN.nc etc will be added",default="gso")
parser.add_argument("-d","--dentype",choices=["DEN","PAWDEN"],default="DEN",help="density file DEN (default) or PAWDEN")
args=parser.parse_args()

file_dict=get_ncfiles(args.filepath,args.dentype)

with netCDF4.Dataset(file_dict['density']['denfile'],'r') as ncroot:
    n4=ncroot.dimensions['number_of_grid_points_vector1'].size
    n5=ncroot.dimensions['number_of_grid_points_vector2'].size
    n6=ncroot.dimensions['number_of_grid_points_vector3'].size
    rprimd=ncroot.variables['primitive_vectors'][:]
    xred=ncroot.variables['reduced_atom_positions'][:]
    atom_species=ncroot.variables['atom_species'][:]
    atom_species_names=netCDF4.chartostring(ncroot.variables['atom_species_names'][:])
    atomic_numbers=ncroot.variables['atomic_numbers'][:]
    byte_string=ncroot.variables['md5_pseudos'][:]
    md5_pseudos=[b"".join(x).decode("utf-8") for x in byte_string]
    byte_string=ncroot.variables['input_string'][:]
    s=b"".join(byte_string).decode("utf-8")
    input_string=s.split()
    pseudofiles=get_pseudos_from_nc(input_string,file_dict['density']['denfile'],md5_pseudos)

dims = (n4,n5,n6)
scalar_field=np.zeros((n6,n5,n4))
voxel_scale=np.linalg.det(rprimd)/(n4*n5*n6)

with netCDF4.Dataset(file_dict['density']['denfile'],'r') as ncroot:
    evalence=(ncroot.variables[file_dict['density']['denfield']][0,:,:,:,0]).T
with netCDF4.Dataset(file_dict['tau']['taufile'],'r') as ncroot:
    tauvalence=(ncroot.variables[file_dict['tau']['taufield']][0,:,:,:,0]).T

atom_list=make_atom_list(xred,atom_species,pseudofiles)

rho_ae, tau_ae = add_spherical_core_fields(
    evalence,
    tauvalence,
    rprimd,
    atom_list
)

grad_rho_ae_fft,grad_components_fft=density_gradient_magnitude(evalence,rprimd,return_components=True)
grad_core=core_density_gradient(evalence.shape,rprimd,atom_list)
grad_ae=grad_components_fft+grad_core
grad_rho_ae = np.sqrt(np.sum(grad_ae**2,axis=0))

# leave out core
#grad_pawvalence = np.sqrt(np.sum(grad_components_fft**2,axis=0))

elf=get_elf(rho_ae,grad_rho_ae,tau_ae)

write_elf_like_pawrhor(
    file_dict['density']['denfile'],
    file_dict['elf']['elffile'],
    elf,
    paw_density_variable=file_dict['density']['denfield']
)

