import numpy as np
from numpy import zeros
from .constants import tol6, kb_HaK

def get_bose(omega, temperatures):
    """
    Get the Bose-Einstein occupations on a range of temperatures
    for a given boson frequency omega (in Hartree).

    Arguments:
        omega: float
        temperatures: [ntemp]

    Returns: bose [ntemp]
    """

    #if not temperatures:
    #    ntemp = 1
    #    temperatures = [0.]
    #else:
    ntemp = len(temperatures)
    bose = zeros(ntemp)

    omega = np.float(omega)
    if omega < tol6:
        return bose

    upper_limit = 50.
    
    for tt, T in enumerate(temperatures):
    
        if T < tol6:
            continue
    
        x = omega / (kb_HaK * T)
    
        if x > upper_limit:
            continue
    
        bose[tt] = 1. / (np.exp(x) - 1)

    return bose

