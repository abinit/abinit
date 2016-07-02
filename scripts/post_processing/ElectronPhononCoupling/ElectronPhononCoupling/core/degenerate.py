
from copy import copy
import numpy as N
from numpy import zeros

def iter_spin_band_eig(eigens, ikpt):
    """
    Iterator over spin index and eigenvalues.
    Yields tuples (ispin, iband, eig) in order of increasing eig at a given k-point.
    """
    nspin, nkpt, nband = eigens.shape
    sbe = [(ispin, 0, eigens[ispin,ikpt,0]) for ispin in range(nspin)]
    cmp_sbe = lambda sbe1, sbe2: cmp(sbe1[2], sbe2[2])
    while sbe:
        min_sbe = sorted(sbe, cmp=cmp_sbe)[0]
        yield min_sbe

        i = sbe.index(min_sbe)
        ispin, iband, eig = min_sbe

        if iband == nband - 1:
            del sbe[i]
        else:
            sbe[i] = (ispin, iband+1, eigens[ispin, ikpt, iband+1])
    

def get_degen(eigens):
    """
    Compute the degeneracy of the bands.

    Arguments
    ---------

    eigens: numpy.ndarray of shape (nspin, nkpt, nband)
        Eigenvalues.

    Returns
    -------
    degen: 2D list (nkpt, )
        For each k-point, contains a list of groups of (s, n) tuples
        which are degenerated.
        For example, if there is only a single spin (spin unpolarized case)
        and two k-points, and at the first k-point there are two triplets,
        then the output is
            [[[(0,1), (0,2), (0,3)],  [(0,4), (0,5), (0,6)]], []]

    """
    nspin, nkpt, nband = eigens.shape

    degen = list()
    for ikpt in range(nkpt):

        kpt_degen = list()
        group = list()
        last_ispin, last_iband, last_eig = 0, 0, -float('inf')

        for sbe in iter_spin_band_eig(eigens, ikpt):
            ispin, iband, eig = sbe

            if N.isclose(last_eig, eig, rtol=1e-12, atol=1e-5):
                if not group:
                    group.append((last_ispin, last_iband))
                group.append((ispin, iband))

            else:
                if group:
                    kpt_degen.append(group)
                    group = list()

            last_ispin, last_iband, last_eig = ispin, iband, eig

        degen.append(kpt_degen)

    return degen


def make_average(arr, degen):
    """ 
    Average a quantity over degenerated states.
    Does not work with spin yet.

    Arguments
    ---------

    arr: numpy.ndarray(..., nkpt, nband)
        An array of any dimension, of which the two last indicies are
        the kpoint and the band.

    degen: list(nkpt, )
        Degeneracy of the states at each k-point, obtained from get_degen.

    Returns
    -------

    arr: numpy.ndarray(..., nkpt, nband)
        The array with the values of degenerated bands averaged.

    """
    nkpt, nband = arr.shape[-2:]

    for ikpt in range(nkpt):
        for group in degen[ikpt]:
            average = copy(arr[...,ikpt,group[0][1]])
            for ispin, iband in group[1:]:
                average += arr[...,ikpt,iband]

            average /= len(group)
            for ispin, iband in group:
                arr[...,ikpt,iband] = average

    return arr


def symmetrize_fan_degen(fan_epc, degen):
    """
    Enforce coupling terms to be zero on the diagonal
    and in degenerate states subset.

    Arguments
    ---------
    fan_epc: np.ndarray, shape=(nkpt,nband,nband,nmode)
        The coupling matrix V_ij V_ji

    degen: list(nkpt, )
        Degeneracy of the states at each k-point, obtained from get_degen.

    Returns
    -------

    fan_epc_symmetrized: np.ndarray, shape=(nkpt,nband,nband,nmode)
    """
    nkpt, nband, mband, nmode = fan_epc.shape
 
    offdiag = N.zeros((nkpt, nband, nband))
    offdiag[:] =  N.ones((nband, nband)) - N.identity(nband)

    for ikpt in range(nkpt):
        for group in degen[ikpt]:
            for degi in group:
                for degj in group:
                    ieig, jeig = degi[1], degj[1]
                    offdiag[ikpt][ieig][jeig] = 0

    fan_epc_sym = N.einsum('ijkl,ijk->ijkl', fan_epc, offdiag)

    return fan_epc_sym
  



