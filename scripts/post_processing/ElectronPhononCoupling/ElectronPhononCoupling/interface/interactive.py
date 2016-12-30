from __future__ import print_function

import numpy as N

from ..config import __version__
from ..core.mpi import comm, mpi_abort_if_exception, i_am_master

from .compute_epc import compute_epc

__all__ = ['run_interactive', 'get_user_input']


def run_interactive():
    """Run the calculation after getting inputs interactively from the user."""
    with mpi_abort_if_exception():
        if i_am_master:
            arguments = get_user_input()
        else:
            arguments = None
        arguments = comm.bcast(arguments, root=0)

    epc = compute_epc(**arguments)



def get_user_input():
    """Get all inputs for the calculation interactively."""

    arguments = {
        'calc_type' : 1,
        'nqpt' : 1,
        'wtq' : [1.0],
        'output' : '',
        'smearing_eV' : 3.6749e-4,
        'temperature' : False,
        'temp_range' : [0, 0, 1],
        'omega_range' : [-0.1, 0.1, 0.001],
        'lifetime' : False,
        'DDB_fnames' : list(),
        'eigq_fnames' : list(),
        'EIGR2D_fnames' : list(),
        'EIGI2D_fnames' : list(),
        'FAN_fnames' : list(),
        'GKK_fnames' : list(),
        'eig0_fname' : '',
        'verbose' : True,
        }

    # Interaction with the user
    logo = """

  _____ _           _                         ____  _                                ____                  _ _             
 | ____| | ___  ___| |_ _ __ ___  _ __       |  _ \| |__   ___  _ __   ___  _ __    / ___|___  _   _ _ __ | (_)_ __   __ _ 
 |  _| | |/ _ \/ __| __| '__/ _ \| '_ \ _____| |_) | '_ \ / _ \| '_ \ / _ \| '_ \  | |   / _ \| | | | '_ \| | | '_ \ / _` |
 | |___| |  __/ (__| |_| | | (_) | | | |_____|  __/| | | | (_) | | | | (_) | | | | | |__| (_) | |_| | |_) | | | | | | (_| |
 |_____|_|\___|\___|\__|_|  \___/|_| |_|     |_|   |_| |_|\___/|_| |_|\___/|_| |_|  \____\___/ \__,_| .__/|_|_|_| |_|\__, |
                                                                                                    |_|              |___/ 

                                                                                                            Version {} 

    """.format(__version__)

    description = """
    This module compute the renormalization and broadening (lifetime)
    of electronic energies due to electron-phonon interaction,
    using either the static or dynamical AHC theory at zero and finite temperatures.
    Also computes the self-energy and spectral function.
    """
    print(logo)
    print(description)

    def get_user(s):
        return raw_input(s.rstrip('\n') + '\n').split('#')[0]

    # Type of calculation the user want to perform
    ui = get_user("""
Please select the calculation type:
    1  Static AHC calculation.
    2  Dynamic AHC calculation.
    3  Static AHC calculation with control over active space.
    4  Frequency-dependent self-energy and spectral function.

Note that option 2, 3 and 4 requires _FAN.nc files obtained
through ABINIT option 'ieig2rf 4'
""")
    calc_type = N.int(ui)
    arguments.update(calc_type=calc_type)
    
    # Define the output file name
    ui = get_user('Enter name of the output file')
    output = ui.strip()
    arguments.update(output=output)
    
    # Enter the value of the smearing parameter for dynamic AHC
    if (calc_type > 1):
      ui = get_user('Enter value of the smearing parameter (in eV)')
      smearing_eV = N.float(ui)
    else:
      smearing_eV = None
    arguments.update(smearing_eV=smearing_eV)
    
    # Temperature dependence analysis?
    ui = get_user('Do you want to compute the change of eigenergies with temperature? [y/n]')
    temperature = ui.split()[0]
    if temperature.lower() == 'y':
      temperature = True
      arguments.update(temperature=temperature)
    else:
      temperature = False

    if temperature:
      ui = get_user('Introduce the starting temperature, max temperature and steps. e.g. 0 2000 100')
      temp_range = map(float, ui.split())
      arguments.update(temp_range=temp_range)


    # frequency range
    if calc_type == 4:
      ui = get_user('Introduce the starting frequency, max frequency and steps. e.g. -1 1 0.05')
      omega_range = map(float, ui.split())
      arguments.update(omega_range=omega_range)
    
    # Broadening lifetime of the electron
    ui = get_user('Do you want to compute the lifetime of the electrons? [y/n]')
    tmp =ui.split()[0]
    if tmp == 'y':
      lifetime = True
    else:
      lifetime = False
    arguments.update(lifetime=lifetime)

    # Get the nb of random Q-points from user 
    ui = get_user('Enter the number of Q-points')
    try:
      nqpt = int(ui)
    except ValueError:
      raise Exception('The value you enter is not an integer!')
    arguments.update(nqpt=nqpt)

    # Get the q-points weights
    wtq = list()
    for ii in N.arange(nqpt):
      ui = get_user('Enter the weight of the %s q-point' %ii)
      wtq.append(float(ui.split()[0]))
    arguments.update(wtq=wtq)
    
    # Get the path of the DDB files from user
    DDB_fnames = []
    for ii in N.arange(nqpt):
      ui = get_user('Enter the name of the %s DDB file' %ii)
      if len(ui.split()) != 1:
        raise Exception("You should provide only 1 file")
      else: # Append and TRIM the input string with STRIP
        DDB_fnames.append(ui.strip(' \t\n\r'))
    arguments.update(DDB_fnames=DDB_fnames)

    # Get the path of the eigq files from user
    eigq_fnames = []
    for ii in N.arange(nqpt):
      ui = get_user('Enter the name of the %s eigq file' %ii)
      if len(ui.split()) != 1:
        raise Exception("You should provide only 1 file")
      else:
        eigq_fnames.append(ui.strip(' \t\n\r'))
    arguments.update(eigq_fnames=eigq_fnames)
    
    # Get the path of the EIGR2D files from user
    EIGR2D_fnames = []
    for ii in N.arange(nqpt):
      ui = get_user('Enter the name of the %s EIGR2D file' %ii)
      if len(ui.split()) != 1:
        raise Exception("You should provide only 1 file")
      else:
        EIGR2D_fnames.append(ui.strip(' \t\n\r'))
    arguments.update(EIGR2D_fnames=EIGR2D_fnames)
    
    # Get the path of the EIGI2D files from user
    if lifetime:
      EIGI2D_fnames = []
      for ii in N.arange(nqpt):
        ui = get_user('Enter the name of the %s EIGI2D file' %ii)
        if len(ui.split()) != 1:
          raise Exception("You should provide only 1 file")
        else:
          EIGI2D_fnames.append(ui.strip(' \t\n\r'))
        arguments.update(EIGI2D_fnames=EIGI2D_fnames)
    
    # Get the path of the FAN files from user if dynamical calculation
    if (calc_type > 1):
      fnames = []
      for ii in N.arange(nqpt):
        ui = get_user('Enter the name of the %s GKK or FAN file' %ii)
        if len(ui.split()) != 1:
          raise Exception("You should provide only 1 file")
        else:
          fnames.append(ui.strip(' \t\n\r'))
        if fnames:
            fname0 = fnames[0]
            if fname0.endswith('FAN.nc'):
                arguments.update(FAN_fnames=fnames)
            elif fname0.endswith('GKK.nc'):
                arguments.update(GKK_fnames=fnames)
            else:
                raise Exception('Extension not recognized. '
                                'Expected GKK.nc or FAN.nc file.')
    
    # Take the EIG at Gamma
    ui = get_user('Enter the name of the unperturbed EIG.nc file at Gamma')
    if len(ui.split()) != 1:
      raise Exception("You sould only provide 1 file")
    else:
      eig0_fname=ui.strip(' \t\n\r')
    arguments.update(eig0_fname=eig0_fname)

    return arguments


