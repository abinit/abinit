from __future__ import print_function

from ..core.constants import Ha2eV
from ..core.util import report_runtime
from ..core import EpcAnalyzer

__all__ = ['compute_epc']


@report_runtime
def compute_epc(
        calc_type = 1,
        nqpt = 1,
        wtq = [1.0],
        write = True,
        output='epc.out',
        temperature = False,
        lifetime = False,
        smearing_eV = 0.01,
        temp_range = [0, 0, 1],
        omega_range=[-0.1,0.1,0.001],
        eig0_fname = '',
        eigq_fnames = list(),
        DDB_fnames = list(),
        EIGR2D_fnames = list(),
        EIGI2D_fnames = list(),
        FAN_fnames = list(),
        GKK_fnames = list(),
        verbose=False,
        **kwargs):
    """
    Compute electron-phonon coupling related quantities, such as:
        - the zero-point renormalization
        - the temperature dependance of eigenvalues
        - the quasiparticle lifetime from the el-ph self-energy


    Keyword arguments
    -----------------

    calc_type: (1)
        Governs the type of calculation performed
            1  -  static AHC calculation
            2  -  dynamic AHC calculation
            3  -  static AHC calculation with control over active space
            4  -  frequency-dependent self-energy and spectral function

    write: (True)
        Write the results on the disk.

    output: ('epc.out')
        Rootname for the output files.

    smearing_eV: (0.01)
        Imaginary parameter for ZPR contributions and broadening.

    temperature: (False)
        Compute temperature dependance of eigenvalues corrections / broadening

    temp_range: [0,0,1]
        Minimum, maximum and step temperature for eigenvalues dependance.

    omega_range: [0,0,1]
        Minimum, maximum and step frequency for the self-energy.

    lifetime: (False)
        Compute broadening

    nqpt: (1)
        Number of q-points

    DDB_fnames: ([])
        Names of _DDB files

    eigq_fnames: ([])
        Names of _EIG.nc files at k+q

    EIGR2D_fnames: ([])
        Names of _EIGR2D.nc files

    EIGI2D_fnames: ([])
        Names of _EIGI2D.nc files

    FAN_fnames: ([])
        Names of _FAN.nc files

    GKK_fnames: ([])
        Names of _GKK.nc files

    eig0_fname: ('')
        Name of the _EIG.nc file for the eigenvalues being corrected


    Returns
    -------

    epc: EpcAnalyzer
        Object containing the response function data
    """

    if smearing_eV is None:
        smearing_Ha = None
    else:
        smearing_Ha = smearing_eV / Ha2eV

    # Initialize epc
    epc = EpcAnalyzer(nqpt=nqpt, 
                       wtq=wtq,
                       eig0_fname=eig0_fname,
                       eigq_fnames=eigq_fnames,
                       DDB_fnames=DDB_fnames,
                       EIGR2D_fnames=EIGR2D_fnames,
                       EIGI2D_fnames=EIGI2D_fnames,
                       FAN_fnames=FAN_fnames,
                       GKK_fnames=GKK_fnames,
                       temp_range=temp_range,
                       omega_range=omega_range,
                       smearing=smearing_Ha,
                       output=output,
                       verbose=verbose,
                       **kwargs)

    # Compute renormalization
    if calc_type == 1:

        if temperature:
            epc.compute_static_td_renormalization()
        else:
            epc.compute_static_zp_renormalization()

    elif calc_type == 2:

        if temperature:
            epc.compute_dynamical_td_renormalization()
        else:
            epc.compute_dynamical_zp_renormalization()

    elif calc_type == 3:

        if temperature:
            epc.compute_static_control_td_renormalization()
        else:
            epc.compute_static_control_zp_renormalization()

    elif calc_type == 4:
        if temperature:
            epc.compute_td_self_energy()
            epc.compute_td_spectral_function()
        else:
            epc.compute_zp_self_energy()
            epc.compute_zp_spectral_function()

    else:
        raise Exception('Calculation type must be 1, 2, 3 or 4')

    # Compute lifetime
    if lifetime:

        if calc_type == 1:

            if temperature:
                epc.compute_static_td_broadening()
            else:
                epc.compute_static_zp_broadening()

        elif calc_type == 2:

            if temperature:
                epc.compute_dynamical_td_broadening()
            else:
                epc.compute_dynamical_zp_broadening()

        elif calc_type == 3:

            if temperature:
                epc.compute_static_control_td_broadening()
            else:
                epc.compute_static_control_zp_broadening()


    # Write the files
    if write:
        epc.write_netcdf()
        epc.write_renormalization()
        if lifetime:
            epc.write_broadening()

    return epc

