from __future__ import print_function

import numpy as np
from numpy import zeros, ones, einsum

from .constants import tol6, tol8, tol12, Ha2eV, kb_HaK

from .mathutil import delta_lorentzian
from . import EigFile, Eigr2dFile, FanFile, DdbFile, GkkFile

__author__ = "Gabriel Antonius"

__all__ = ['QptAnalyzer']


class QptAnalyzer(object):

    def __init__(self,
                 ddb_fname=None,
                 eigq_fname=None,
                 eigk_fname=None,
                 eigr2d_fname=None,
                 eigr2d0_fname=None,
                 eigi2d_fname=None,
                 fan_fname=None,
                 fan0_fname=None,
                 gkk_fname=None,
                 gkk0_fname=None,
                 wtq=1.0,
                 smearing=0.00367,
                 temperatures=None,
                 omegase=None,
                 amu=None,
                 asr=True,
                 mu=None,
                 double_smearing = False,
                 smearing_width = 0.0367,
                 smearing_above = 0.00367,
                 smearing_below = 0.00367,
                 ):

        # Files
        self.ddb = DdbFile(ddb_fname, read=False, asr=asr)
        self.eigq = EigFile(eigq_fname, read=False)
        self.eigr2d = Eigr2dFile(eigr2d_fname, read=False)
        self.eigi2d = Eigr2dFile(eigi2d_fname, read=False)
        self.fan = FanFile(fan_fname, read=False)
        self.eig0 = EigFile(eigk_fname, read=False)
        self.eigr2d0 = Eigr2dFile(eigr2d0_fname, read=False)
        self.fan0 = FanFile(fan0_fname, read=False)
        self.gkk = GkkFile(gkk_fname, read=False)
        self.gkk0 = GkkFile(gkk0_fname, read=False)

        self.wtq = wtq
        self.smearing = smearing
        self.omegase = omegase if omegase else list()
        self.temperatures = temperatures if temperatures else list()
        self.mu = mu
        self.amu = amu

        self.double_smearing = double_smearing
        self.smearing_width = smearing_width
        self.smearing_above = smearing_above
        self.smearing_below = smearing_below

    @property
    def nkpt(self):
        if self.eigr2d.fname:
            return self.eigr2d.nkpt
        elif self.fan.fname:
            return self.fan.nkpt
        elif self.gkk.fname:
            return self.gkk.nkpt
        else:
            raise Exception("Don't know nkpt. No files to read.")

    @property
    def kpts(self):
        if self.eigr2d.fname:
            return self.eigr2d.kpt
        elif self.fan.fname:
            return self.fan.kpt
        elif self.gkk.fname:
            return self.gkk.kpt
        else:
            raise Exception("Don't know kpt. No files to read.")

    @property
    def nband(self):
        if self.eigr2d.fname:
            return self.eigr2d.nband
        elif self.fan.fname:
            return self.fan.nband
        elif self.gkk.fname:
            return self.gkk.nband
        else:
            raise Exception("Don't know nband. No files to read.")

    @property
    def natom(self):
        return self.ddb.natom

    @property
    def nmode(self):
        return self.ddb.nmode

    @property
    def is_gamma(self):
        return self.ddb.is_gamma

    @property
    def qred(self):
        return self.ddb.qred

    @property
    def omega(self):
        return self.ddb.omega

    @property
    def nomegase(self):
        return len(self.omegase)

    @property
    def ntemp(self):
        return len(self.temperatures)

    @property
    def use_gkk(self):
        return (bool(self.gkk.fname) and bool(self.gkk0.fname))

    @property
    def has_active(self):
        return (bool(self.fan.fname) and bool(self.fan0.fname)) or self.use_gkk

    def read_nonzero_files(self):
        """Read all nc files that are not specifically related to q=0."""
        for f in (self.ddb, self.eigq, self.eigr2d, self.eigi2d,
                  self.fan, self.gkk):
            if f.fname:
                f.read_nc()

        if self.amu is not None:
            self.ddb.set_amu(self.amu)

        self.ddb.compute_dynmat()

    def read_ddb(self):
        """Read the ddb and diagonalize the matrix, setting omega."""
        self.ddb.read_nc()
        if self.amu is not None:
            self.ddb.set_amu(self.amu)
        self.ddb.compute_dynmat()

    def read_zero_files(self):
        """Read all nc files that are related to q=0."""
        for f in (self.eig0, self.eigr2d0, self.fan0, self.gkk0):
            if f.fname:
                f.read_nc()

    def broadcast_zero_files(self):
        """Broadcast the data related to q=0 from master to all workers."""

        if self.eig0.fname:
            self.eig0.broadcast()
            self.eig0.get_degen()

        if self.eigr2d0.fname:
            self.eigr2d0.broadcast()

        if self.fan0.fname:
            self.fan0.broadcast()

        if self.gkk0.fname:
            self.gkk0.broadcast()

    def get_occ_kq_nospin(self):
        """
        Get the occupations, being either 0 or 1, regardless of spinor.
        Assumes a gapped system, where occupations are the same at all kpts. 
        Returns: occ[nband]
        """
        if self.eigr2d.fname:
            occ = self.eigr2d.occ[0,0,:]
        elif self.fan.fname:
            occ = self.fan.occ[0,0,:]
        elif self.gkk.fname:
            occ = self.gkk.occ[0,0,:]
        else:
            raise Exception("Don't know nband. No files to read.")

        if any(occ == 2.0):
            occ = occ / 2.0

        return occ

    def get_max_val(self):
        """Get the maximum valence band energy."""
        occ0 = self.get_occ_kq_nospin()
        eig = self.eigq.EIG[0,0,:]

        E_last = eig[0]
        for f, E in zip(occ0, eig):
            if f < 0.5:
                break
            E_last = E

        return E_last

    def get_min_cond(self):
        """Get the minimum conduction band energy."""
        occ0 = self.get_occ_kq_nospin()
        eig = self.eigq.EIG[0,0,:]

        for f, E in zip(occ0, eig):
            if f <= 0.5:
                break

        return E

    def find_fermi_level(self):
        """
        Find the Fermi level locally, using the eigenvalues
        at all k+q points available. Assuming a gapped system.
        """
        return (self.get_max_val() + self.get_min_cond()) / 2.0

    def get_eta(self, omega_se):
        """
        Get the imaginary broadening parameter,
        which may depend on frequency. 
        
        Returns: eta[nkpt, nband, nomegase]
        """

        # nkpt, nband
        occ0 = self.eig0.get_fermi_function_T0(self.mu)[0,:,:]

        # nkpt, nband
        signs = (2 * occ0 - 1)

        nomegase = len(omega_se)
        values = np.zeros(nomegase, dtype=float)
        #abso = np.abs(omega_se)

        if self.double_smearing:

            eta1 = self.smearing_below
            eta2 = self.smearing
            eta3 = self.smearing_above
            w = self.smearing_width

            for i, o in enumerate(omega_se):

                if o < -w:
                    values[i] = eta1
                elif o <= 0.:
                    values[i] = eta2 + (eta1-eta2) * abs(o) / w
                elif o < w:
                    values[i] = eta2 + (eta3-eta2) * abs(o) / w
                else:
                    values[i] = eta3

        else:
            values[:] = self.smearing

        eta = einsum('kn,l->knl', signs, values) * 1j 

        return eta

    @staticmethod
    def get_se_indices(mode=False, temperature=False, omega=False):
        """Get the indices string corresponding to the self-energy dimensions"""
        indices = ''
        if mode:
            indices += 'o'
        if temperature:
            indices += 't'
        if omega:
            indices += 'l'
        indices += 'kn'
        return indices

    @staticmethod
    def reduce_array(arr, mode=False, temperature=False, omega=False,
                     shape=None):
        """
        Eliminate dimensions from an array of shape
        (nmode, ntemp, nomegase, nkpt, nband)
        by summing over any or all of the first three dimension.

        mode:
            Keep the first dimension
        temperature:
            Keep the second dimension
        omega:
            Keep the third dimension
        shape: 'otlkn'
            Any subset of the initial indices, in any order.
            Reshaping the array should be avoided if unnecessary.
        """

        initial_indices = 'otlkn'

        final_indices = QptAnalyzer.get_se_indices(
            mode=mode, temperature=temperature, omega=omega)

        if shape is not None:

            if bool(shape.strip(initial_indices)):
                raise Exception('shape ({}) must be a new combination of ' +
                                'the intitial indices ({}).'.format(
                                 shape, initial_indices))

            final_indices = shape

        summation = initial_indices + '->' + final_indices

        return einsum(summation, arr)

    def get_fan_ddw_sternheimer(self,
        mode=False, omega=False, temperature=False, shape=None):
        """
        Compute the fan and ddw contribution to the self-energy
        obtained from the Sternheimer equation,
        that is, the contribution of the upper bands.

        Do not include the q-point weight.

        Returns: fan, ddw

        The return arrays vary in dimensions, depending on the input arguments.
        These arrays are at most of dimension 5, as

            fan[nmode, ntemp, nomegase, nkpt, nband]
            ddw[nmode, ntemp, nomegase, nkpt, nband]

        Depending on the truth value of the input arguments,
        the dimension nomegase (omega) and ntemp (temperature)
        will be eliminated. 
        The dimension nmode will be summed over in case mode=False.

        In the semi-static approximation, these quantities do not actually
        depend on omega, so the arrays are simply repated along the omega axis.
        """

        nkpt = self.nkpt
        nband = self.nband
        natom = self.natom
        nmode = self.nmode
        nomegase = self.nomegase
        ntemp = self.ntemp

        # Get reduced displacement (scaled with frequency)
        displ_red_FAN2, displ_red_DDW2 = self.ddb.get_reduced_displ_squared()
    
        # FIXME this will not work for nsppol=2
        # nmode, nkpt, nband
        fan = einsum('knabij,objai->okn', self.eigr2d.EIG2D, displ_red_FAN2)
        ddw = einsum('knabij,objai->okn', self.eigr2d0.EIG2D, displ_red_DDW2)

        # Temperature dependence factor
        n_B = self.ddb.get_bose(self.temperatures)
        tdep = 2 * n_B + 1 if temperature else ones((nmode,1))

        # Omega dependence factor
        odep = ones(nomegase) if omega else ones(1)

        # nmode, ntemp, nkpt, nband
        fan = einsum('okn,ot->otkn', fan, tdep)
        ddw = einsum('okn,ot->otkn', ddw, tdep)

        # nmode, ntemp, nomega, nkpt, nband
        fan = einsum('otkn,l->otlkn', fan, odep)
        ddw = einsum('otkn,l->otlkn', ddw, odep)

        # Reduce the arrays
        fan = self.reduce_array(fan, mode=mode, temperature=temperature,
                                omega=omega, shape=shape)
        ddw = self.reduce_array(ddw, mode=mode, temperature=temperature,
                                omega=omega, shape=shape)
        return fan, ddw
    
    def get_fan_ddw_gkk2_active(self):
        """
        Compute the squared gkk elements for the fan ddw terms.

        Returns:
            fan[nkpt, nband, nband, nmode]
            ddw[nkpt, nband, nband, nmode]
        """

        if not self.has_active:
            raise Exception('You should provide GKK files or FAN files '
                            'to compute active space contribution.')

        if self.use_gkk:

            #gkk2 = self.gkk.get_gkk_squared()
            #gkk02 = self.gkk0.get_gkk_squared()

            fan = np.abs(self.gkk.get_gkk_mode(self.ddb)) ** 2
            ddw = self.gkk0.get_gkk2_DW_mode(self.ddb)

        else:
            # Get reduced displacement (scaled with frequency)
            displ_red_FAN2, displ_red_DDW2 = self.ddb.get_reduced_displ_squared()

            gkk2 = self.fan.FAN
            gkk02 = self.fan0.FAN

            # nkpt, nband, nband, nmode
            fan = einsum('kniajbm,oabij->knmo', gkk2, displ_red_FAN2)
            ddw = einsum('kniajbm,oabij->knmo', gkk02, displ_red_DDW2)

        # Enforce the diagonal coupling terms to be zero at Gamma
        ddw = self.eig0.symmetrize_fan_degen(ddw)
        if self.is_gamma:
            fan = self.eig0.symmetrize_fan_degen(fan)
      
        return fan, ddw
    
    def get_fan_ddw_active(self, mode=False, omega=False, temperature=False,
                           dynamical=True, shape=None):
        """
        Compute the fan and ddw contributions to the self-energy
        from the active space, that is, the the lower bands.

        Do not include the q-point weight.

        Returns: fan, ddw

        The return arrays vary in dimensions, depending on the input arguments.
        These arrays are at most of dimension 5, as

            fan[nmode, ntemp, nomegase, nkpt, nband]
            ddw[nmode, ntemp, nomegase, nkpt, nband]

        Depending on the truth value of the input arguments,
        the dimension nomegase (omega) and ntemp (temperature)
        will be eliminated. 
        The dimension nmode will be summed over in case mode=False.

        The Debye-Waller term does not actually depends on omega,
        but this dimension is kept anyway.
        """

        nkpt = self.nkpt
        nband = self.nband
        nmode = self.nmode

        if temperature:
            ntemp = self.ntemp
            temperatures = self.temperatures

            # Bose-Enstein occupation number
            # nmode, ntemp
            n_B = self.ddb.get_bose(temperatures)

        else:
            ntemp = 1
            temperatures = zeros(1)
            n_B = zeros((nmode,1))

        if omega:
            nomegase = self.nomegase
            omega_se = self.omegase
        else:
            # omega_se is measured from the bare eigenvalues
            nomegase = 1
            omega_se = zeros(1)

        if dynamical:
            omega_q = self.ddb.omega[:].real
        else:
            omega_q = zeros(nmode)


        # Fermi-Dirac occupation number
        # nspin, nkpt, nband, ntemp
        occ = self.eigq.get_fermi_function(self.mu, temperatures)

        # G^2
        # nkpt, nband, nband, nmode
        fan_g2, ddw_g2 = self.get_fan_ddw_gkk2_active()


        # DDW term
        # --------

        # nkpt, nband
        occ0 = self.eig0.get_fermi_function_T0(self.mu)[0,:,:]
    
        # nkpt, nband, nband
        delta_E_ddw = (einsum('kn,m->knm', self.eig0.EIG[0,:,:].real, ones(nband))
                     - einsum('kn,m->kmn', self.eig0.EIG[0,:,:].real, ones(nband))
                     - einsum('m,kn->knm', ones(nband), (2*occ0-1)) * self.smearing * 1j)

        # nmode, nkpt, nband
        ddw = einsum('knmo,knm->okn', ddw_g2, 1.0 / delta_E_ddw)

        # nmode, ntemp
        tdep = 2 * n_B + 1

        # FIXME This is not optimal: The mode indices will be summed
        #       so there is no need to create an array this big.
        #       in case omega=True and mode=False
        # nmode, ntemp, nkpt, nband
        ddw = einsum('okn,ot->otkn', ddw, tdep)

        odep = ones(nomegase) if omega else ones(0)

        # ntemp, nomega, nkpt, nband
        ddw = einsum('otkn,l->otlkn', ddw, ones(nomegase))

        # Reduce the arrays
        ddw = self.reduce_array(ddw, mode=mode,
                                temperature=temperature, omega=omega,
                                shape=shape)


        # Fan term
        # --------

        # nmode, ntemp, nomegase, nkpt, nband
        fan = zeros((nmode, ntemp, nomegase, nkpt, nband), dtype=complex)

        # n + 1 - f
        # nkpt, nband, nmode, ntemp
        num1 = (einsum('ot,kn->knot', n_B, ones((nkpt,nband)))
              + 1. - einsum('knt,o->knot', occ[0,:,:,:], ones(nmode)))

        # n + f
        # nkpt, nband, nmode, ntemp
        num2 = (einsum('ot,kn->knot', n_B, ones((nkpt,nband)))
              + einsum('knt,o->knot', occ[0,:,:,:], ones(nmode)))

        # nkpt, nband
        #eta = (2 * occ0 - 1) * self.smearing * 1j

        # nkpt, nband, nomegase
        eta = self.get_eta(omega_se)

        for jband in range(nband):

            # nkpt, nband
            delta_E = (
                self.eig0.EIG[0,:,:].real
              - einsum('k,n->kn', self.eigq.EIG[0,:,jband].real, ones(nband))
              )
    
            # nkpt, nband, nomegase
            delta_E_omega = (
                  einsum('kn,l->knl', delta_E, ones(nomegase))
                + einsum('kn,l->knl', ones((nkpt,nband)), omega_se)
                - eta
                )
    
            # Emission term
            # nkpt, nband, nomegase, nmode
            deno1 = (einsum('knl,o->knlo', delta_E_omega, ones(nmode))
                   - einsum('knl,o->knlo', ones((nkpt,nband,nomegase)), omega_q))

            # nmode, nkpt, nband, nomegase, ntemp
            div1 = einsum('kot,knlo->oknlt', num1[:,jband,:,:], 1.0 / deno1)
    
            del deno1
    
            # Absorption term
            # nkpt, nband, nomegase, nmode
            deno2 = (einsum('knl,o->knlo', delta_E_omega, ones(nmode))
                   + einsum('knl,o->knlo', ones((nkpt,nband,nomegase)), omega_q))
    
            # nmode, nkpt, nband, nomegase, ntemp
            div2 = einsum('kot,knlo->oknlt', num2[:,jband,:,:], 1.0 / deno2)

            del deno2
    
            # FIXME This is not optimal: The mode indices will be summed
            #       so there is no need to create an array this big.
            # in case omega=True and mode=False

            # nmode, ntemp, nomegase, nkpt, nband
            fan += einsum('kno,oknlt->otlkn', fan_g2[:,:,jband,:], div1 + div2)
    
            del div1, div2
      
        # Reduce the arrays
        fan = self.reduce_array(fan, mode=mode, temperature=temperature,
                                omega=omega, shape=shape)

        return fan, ddw

    def get_fan_ddw(self, mode=False, temperature=False,
                    omega=False, dynamical=False, shape=None):
        """
        Compute the sum of the Fan and the Diagonal Debye-Waller term.

        Keyword arguments
        -----------------

        mode:
            Keep mode decomposition (additional dimension)
        temperature:
            Inlcude temperature depdencente (additional dimension).
        omega:
            Inlcude frequency depdencente (additional dimension).
        dynamical:
            Use the full dynamical theory by including the phonon frequencies
            in the location of the poles of the self-energy.

        """

        kwargs = dict(
            mode=mode,
            temperature=temperature,
            omega=omega,
            shape=shape,
            )

        fan_stern, ddw_stern = self.get_fan_ddw_sternheimer(**kwargs)
        fan_active, ddw_active = self.get_fan_ddw_active(dynamical=dynamical,
                                                         **kwargs)

        fan = fan_active + fan_stern
        ddw = ddw_active + ddw_stern

        return fan, ddw


    def get_self_energy(self,
                        mode=False,
                        temperature=False,
                        omega=False,
                        dynamical=True,
                        only_sternheimer=False,
                        only_active=False,
                        only_fan=False,
                        only_ddw=False,
                        real=False,
                        imag=False,
                        shape=None,
                        ):
        """
        Compute the self energy with various options to separate the terms
        into separate fan / ddw, or active / sternheimer.

        Keyword arguments
        -----------------

        mode:
            Keep mode decomposition (additional dimension)
        temperature:
            Inlcude temperature depdencente (additional dimension).
        omega:
            Inlcude frequency depdencente (additional dimension).
        dynamical:
            Use the full dynamical theory by including the phonon frequencies
            in the location of the poles of the self-energy.

        only_sternheimer:
            Only include the (frequency independent) Sternheimer part.
        only_active:
            Only include the "active" part contribution.
        only_fan:
            Only include the (frequency dependent) Fan term.
        only_ddw:
            Only include the (frequency independent) diagonal Debye-Waller term.

        """

        kwargs = dict(
            mode=mode,
            temperature=temperature,
            omega=omega,
            #shape=shape,
            )

        if only_sternheimer and only_active:
            raise Exception(
            'only_sternheimer and only_active cannot be True at the same time')

        elif only_sternheimer:
            fan, ddw = self.get_fan_ddw_sternheimer(**kwargs)

        elif only_active:
            fan, ddw = self.get_fan_ddw_active(dynamical=dynamical, **kwargs)

        else:
            fan, ddw = self.get_fan_ddw(dynamical=dynamical, **kwargs)

        if only_fan:
            se_q = fan
        elif only_ddw:
            se_q = - ddw
        else:
            se_q = fan - ddw

        se = self.wtq * se_q
        se = self.eig0.make_average(se)

        if real:
            se = se.real
        elif imag:
            se = se.imag

        # Note that we must not reshape before calling eig0.make_average
        if shape is not None:
            initial_indices = self.get_se_indices(
                mode=mode, temperature=temperature, omega=omega)
            summation = initial_indices + '->' + shape
            se = einsum(summation, se)

        return se


    def get_broadening(self, mode=False, temperature=False,
                       omega=False, dynamical=True, shape=None):
        """
        Compute the zp broadening contribution from one q-point in a dynamical scheme.
        Only take the active space contribution.
        """

        nkpt = self.nkpt
        nband = self.nband
        nmode = self.nmode

        if temperature:
            ntemp = self.ntemp
            temperatures = self.temperatures

            # Bose-Enstein occupation number
            # nmode, ntemp
            n_B = self.ddb.get_bose(temperatures)

        else:
            ntemp = 1
            temperatures = zeros(1)
            n_B = zeros((nmode,1))

        if omega:
            nomegase = self.nomegase
            omega_se = self.omegase
        else:
            # omega_se is measured from the bare eigenvalues
            nomegase = 1
            omega_se = zeros(1)

        if dynamical:
            omega_q = self.ddb.omega[:].real
        else:
            omega_q = zeros(nmode)

        # Fermi-Dirac occupation number
        # nspin, nkpt, nband, ntemp
        occ = self.eigq.get_fermi_function(self.mu, temperatures)
    
        # nkpt, nband, ntemp
        f = occ[0]

        # nkpt, nband
        occ0 = self.eig0.get_fermi_function_T0(self.mu)[0,:,:]
    
        # nkpt, nband, nband, nmode
        fan_g2, ddw_g2 = self.get_fan_ddw_gkk2_active()
      
        # nmode, ntemp, nkpt
        n_B = einsum('ot,q->otq', n_B, ones(nkpt))

        # nmode, ntemp, nkpt,nband
        f = einsum('qmt,o->otqm', f, ones(nmode))

        # nkpt, nband
        sign = np.sign(- (2 * occ0 - 1.))

        broadening = zeros((nmode,ntemp,nomegase,nkpt,nband))

        for jband in range(nband):

            # nmode, ntemp, nkpt
            num1 = (n_B + f[...,jband])
            num2 = (n_B + 1 - f[...,jband])

            # nkpt, nband
            delta_E = (
                self.eig0.EIG[0,:,:].real
                - einsum('q,n->qn', self.eigq.EIG[0,:,jband].real, ones(nband))
                )

            # nkpt, nband, nomegase
            delta_E_omega = (einsum('kn,l->knl', delta_E, ones(nomegase))
                           + einsum('kn,l->knl', ones((nkpt,nband)), omega_se))

            # nmode, nkpt, nband, nomegase
            deno1 = (
                einsum('knl,o->oknl', delta_E_omega, ones(nmode))
              + einsum('o,knl->oknl', omega_q, ones((nkpt,nband,nomegase)))
                )

            # nmode, nkpt, nband, nomegase
            deno2 = (
                einsum('knl,o->oknl', delta_E_omega, ones(nmode))
              - einsum('o,knl->oknl', omega_q, ones((nkpt,nband,nomegase)))
                )

            # nmode, nkpt, nband, nomegase
            delta1 = np.pi * delta_lorentzian(deno1, self.smearing)
            delta2 = np.pi * delta_lorentzian(deno2, self.smearing)
    
            # nmode, ntemp, nomegase, nkpt, nband
            term1 = einsum('otk,oknl->otlkn', num1, delta1)
            term2 = einsum('otk,oknl->otlkn', num2, delta2)

            deltas = einsum('kn,otlkn->otlkn', sign, term1 + term2)
            broadening_j = einsum('kno,otlkn->otlkn', fan_g2[:,:,jband,:], deltas)

            broadening += broadening_j.real

        # Reduce the arrays
        broadening = self.reduce_array(broadening, mode=mode,
                                       temperature=temperature,
                                       omega=omega)

        broadening *= self.wtq
        broadening = self.eig0.make_average(broadening)

        # Note that we must not reshape before calling eig0.make_average
        if shape is not None:
            initial_indices = self.get_se_indices(
                mode=mode, temperature=temperature, omega=omega)
            summation = initial_indices + '->' + shape
            broadening = einsum(summation, broadening)

        return broadening


    def get_zp_self_energy(self):
        """
        Compute the zp frequency-dependent dynamical self-energy
        from one q-point.
    
        The self-energy is evaluated on a frequency mesh 'omegase'
        that is shifted by the bare energies, such that, what is retured is
    
            Simga'_kn(omega) = Sigma_kn(omega + E^0_kn)

        Returns: sigma[nkpt,nband,nomegase]
        """
        self.sigma = self.get_self_energy(
            mode=False,
            temperature=False,
            omega=True,
            dynamical=True,
            only_sternheimer=False,
            only_active=False,
            shape='knl'
            )
        return self.sigma

    def get_td_self_energy(self):
        """
        Compute the temperature depended and frequency-dependent
        dynamical self-energy from one q-point.
    
        The self-energy is evaluated on a frequency mesh 'omegase'
        that is shifted by the bare energies, such that, what is retured is
    
            Simga'_kn(omega,T) = Sigma_kn(omega + E^0_kn, T)
    
        Returns: sigma[nkpt,nband,nomegase,ntemp]
        """
        self.sigma = self.get_self_energy(
            mode=False,
            temperature=True,
            omega=True,
            dynamical=True,
            only_sternheimer=False,
            only_active=False,
            shape='knlt',
            )
        # nkpt, nband, nomegase, nband
        return self.sigma

    def get_zp_self_energy_active(self):
        """
        Compute the zp frequency-dependent dynamical self-energy
        from one q-point.
        Only include the active space contribution.
    
        The self-energy is evaluated on a frequency mesh 'omegase'
        that is shifted by the bare energies, such that, what is retured is
    
            Simga'_kn(omega) = Sigma_kn(omega + E^0_kn)

        Returns: sigma[nkpt,nband,nomegase]
        """
        self.sigma = self.get_self_energy(
            mode=False,
            temperature=False,
            omega=True,
            dynamical=True,
            only_sternheimer=False,
            only_active=True,
            shape='knl',
            )
        return self.sigma

    def get_zp_self_energy_sternheimer(self):
        """
        Compute the zp frequency-dependent dynamical self-energy
        from one q-point.
        Only include the Sternheimer contribution.
    
        The self-energy is evaluated on a frequency mesh 'omegase'
        that is shifted by the bare energies, such that, what is retured is
    
            Simga'_kn(omega) = Sigma_kn(omega + E^0_kn)

        Returns: sigma[nkpt,nband,nomegase]
        """
        self.sigma = self.get_self_energy(
            mode=False,
            temperature=False,
            omega=True,
            dynamical=True,
            only_sternheimer=True,
            only_active=False,
            shape='knl',
            )
        return self.sigma

    def get_td_self_energy_active(self):
        """
        Compute the temperature depended and frequency-dependent
        dynamical self-energy from one q-point.
        Only include the active space contribution.
    
        The self-energy is evaluated on a frequency mesh 'omegase'
        that is shifted by the bare energies, such that, what is retured is
    
            Simga'_kn(omega,T) = Sigma_kn(omega + E^0_kn, T)
    
        Returns: sigma[nkpt,nband,nomegase,ntemp]
        """
        self.sigma = self.get_self_energy(
            mode=False,
            temperature=True,
            omega=True,
            dynamical=True,
            only_sternheimer=False,
            only_active=True,
            shape='knlt',
            )
        return self.sigma

    def get_td_self_energy_sternheimer(self):
        """
        Compute the temperature depended and frequency-dependent
        dynamical self-energy from one q-point.
        Only include the Sternheimer contribution.
    
        The self-energy is evaluated on a frequency mesh 'omegase'
        that is shifted by the bare energies, such that, what is retured is
    
            Simga'_kn(omega,T) = Sigma_kn(omega + E^0_kn, T)
    
        Returns: sigma[nkpt,nband,nomegase,ntemp]
        """
        self.sigma = self.get_self_energy(
            mode=False,
            temperature=True,
            omega=True,
            dynamical=True,
            only_sternheimer=True,
            only_active=False,
            shape='knlt',
            )
        return self.sigma

    def get_zpr_static_sternheimer(self):
        """Compute the q-point zpr contribution in a static scheme."""

        self.zpr = self.get_self_energy(
            mode=False,
            temperature=False,
            omega=False,
            dynamical=False,
            only_sternheimer=True,
            only_active=False,
            real=True,
            )
        return self.zpr

    def get_zpr_static_sternheimer_modes(self):
        """Compute the q-point zpr contribution in a static scheme."""

        self.zpr = self.get_self_energy(
            mode=True,
            temperature=False,
            omega=False,
            dynamical=False,
            only_sternheimer=True,
            only_active=False,
            ).real
        return self.zpr

    def get_zpr_static(self):
        """
        Compute the q-point zpr contribution in a static scheme,
        with the transitions split between active and sternheimer.
        """
        self.zpr = self.get_self_energy(
            mode=False,
            temperature=False,
            omega=False,
            dynamical=False,
            only_sternheimer=False,
            only_active=False,
            real=True,
            )
        return self.zpr

    def get_zpr_dynamical(self):
        """
        Compute the q-point zpr contribution in a static scheme
        with the transitions split between active and sternheimer.
        """
        self.zpr = self.get_self_energy(
            mode=False,
            temperature=False,
            omega=False,
            dynamical=True,
            only_sternheimer=False,
            only_active=False,
            real=True,
            )
        return self.zpr

    def get_tdr_static(self):
        """
        Compute the q-point contribution to the temperature-dependent
        renormalization in a static scheme,
        with the transitions split between active and sternheimer.
        """
        self.tdr = self.get_self_energy(
            mode=False,
            temperature=True,
            omega=False,
            dynamical=False,
            only_sternheimer=False,
            only_active=False,
            real=True,
            shape='knt',
            )
        return self.tdr
    
    def get_tdr_dynamical(self):
        """
        Compute the q-point contribution to the temperature-dependent
        renormalization in a dynamical scheme.
        """
        self.tdr = self.get_self_energy(
            mode=False,
            temperature=True,
            omega=False,
            dynamical=True,
            only_sternheimer=False,
            only_active=False,
            real=True,
            shape='knt',
            )
        return self.tdr

    def get_tdr_static_nosplit(self):
        """
        Compute the q-point contribution to the temperature-dependent
        renormalization in a static scheme.
        """
        self.tdr = self.get_self_energy(
            mode=False,
            temperature=True,
            omega=False,
            dynamical=False,
            only_sternheimer=True,
            only_active=False,
            real=True,
            shape='knt',
            )
        return self.tdr

    def get_tdr_dynamical_active(self):
        """
        Compute the q-point contribution to the temperature-dependent
        renormalization in a dynamical scheme,
        taking only the active space contribution.
        """
        self.tdr = self.get_self_energy(
            mode=False,
            temperature=True,
            omega=False,
            dynamical=True,
            only_sternheimer=False,
            only_active=True,
            real=True,
            shape='knt',
            )
        return self.tdr

    def get_tdr_static_sternheimer(self):
        """Compute the q-point zpr contribution in a static scheme."""

        self.tdr = self.get_self_energy(
            mode=False,
            temperature=True,
            omega=False,
            dynamical=False,
            only_sternheimer=True,
            only_active=False,
            real=True,
            shape='knt',
            )
        return self.tdr

    def get_zpr_dynamical_active(self):
        """
        Compute the q-point contribution to the zero point
        renormalization in a dynamical scheme,
        taking only the active space contribution.
        """
        self.zpr = self.get_self_energy(
            mode=False,
            temperature=False,
            omega=False,
            dynamical=True,
            only_sternheimer=False,
            only_active=True,
            real=True,
            )
        # nkpt, nband, ntemp
        return self.zpr

    def get_zpr_dynamical_active_modes(self):
        """
        Compute the q-point contribution to the zero point
        renormalization in a dynamical scheme,
        taking only the active space contribution.
        """
        self.zpr = self.get_self_energy(
            mode=True,
            temperature=False,
            omega=False,
            dynamical=True,
            only_sternheimer=False,
            only_active=True,
            ).real
        # nkpt, nband, ntemp
        return self.zpr

    def get_zpr_static_modes(self):
        """
        Compute the q-point zpr contribution in a static scheme,
        with the transitions split between active and sternheimer.
        Retain the mode decomposition of the zpr.
        """
        self.zpr = self.get_self_energy(
            mode=True,
            temperature=False,
            omega=False,
            dynamical=False,
            only_sternheimer=False,
            only_active=False,
            real=True,
            )
        # nmode, nkpt, nband
        return self.zpr

    def get_zpr_dynamical_modes(self):
        """
        Compute the q-point zpr contribution in a dynamical scheme,
        with the transitions split between active and sternheimer.
        Retain the mode decomposition of the zpr.
        """
        self.zpr = self.get_self_energy(
            mode=True,
            temperature=False,
            omega=False,
            dynamical=True,
            only_sternheimer=False,
            only_active=False,
            real=True,
            )
        # nmode, nkpt, nband
        return self.zpr

    def get_zpb_dynamical(self):
        """
        Compute the zp broadening contribution from one q-point in a dynamical scheme.
        Only take the active space contribution.
        Returns: zpb[nkpt,nband]
        """
        self.zpb = self.get_broadening(mode=False, temperature=False,
                                       omega=False, dynamical=True)
        return self.zpb
    
    def get_tdb_dynamical(self):
        """
        Compute the td broadening contribution from one q-point in a dynamical scheme.
        Only take the active space contribution.
        Returns: zpb[nkpt,nband,ntemp]
        """
        self.tdb = self.get_broadening(
            mode=False,
            temperature=True,
            omega=False,
            dynamical=True,
            shape='knt',
            )
        return self.tdb
    
    def get_zpb_static(self):
        """
        Compute the zp broadening contribution from one q-point in a static scheme.
        Only take the active space contribution.
        Returns: zpb[nkpt,nband]
        """
        self.zpb = self.get_broadening(mode=False, temperature=False,
                                       omega=False, dynamical=False)
        return self.zpb

    def get_tdb_static(self):
        """
        Compute the td broadening contribution from one q-point in a static scheme.
        Only take the active space contribution.
        Returns: zpb[nkpt,nband,ntemp]
        """
        self.tdb = self.get_broadening(
            mode=False,
            temperature=True,
            omega=False,
            dynamical=False,
            shape='knt',
            )
        return self.tdb

    def get_tdb_static_nosplit(self):
        """
        Compute the q-point contribution to the temperature-dependent broadening
        in a static scheme from the EIGI2D files.
        """
    
        nkpt = self.nkpt
        nband = self.nband
        natom = self.natom
        ntemp = self.ntemp
          
        # These indicies be swapped at the end
        self.tdb = zeros((ntemp, nkpt, nband), dtype=complex)
    
        # Get reduced displacement (scaled with frequency)
        displ_red_FAN2, displ_red_DDW2 = self.ddb.get_reduced_displ_squared()
    
        bose = self.ddb.get_bose(self.temperatures)
    
        fan_corrQ = einsum('ijklmn,olnkm->oij', self.eigi2d.EIG2D, displ_red_FAN2)
    
        for imode in np.arange(3*natom):
          for tt, T in enumerate(self.temperatures):
            self.tdb[tt,:,:] += np.pi * fan_corrQ[imode,:,:] * (2*bose[imode,tt] + 1.)
    
        self.tdb = self.tdb * self.wtq
    
        self.tdb = self.eig0.make_average(self.tdb)
    
        # nkpt, nband, ntemp
        self.tdb = np.einsum('tkn->knt', self.tdb)

        return self.tdb

    def get_zpb_static_nosplit(self):
        """
        Compute the zp broadening contribution from one q-point in a static scheme
        from the EIGI2D files.
        """
    
        nkpt = self.nkpt
        nband = self.nband
        natom = self.natom
    
        self.zpb = zeros((nkpt, nband), dtype=complex)
    
        # Get reduced displacement (scaled with frequency)
        displ_red_FAN2, displ_red_DDW2 = self.ddb.get_reduced_displ_squared()
        
        fan_corrQ = einsum('ijklmn,olnkm->oij', self.eigi2d.EIG2D, displ_red_FAN2)
    
        self.zpb += np.pi * np.sum(fan_corrQ, axis=0)
        self.zpb = self.zpb * self.wtq
    
        if np.any(self.zpb[:,:].imag > tol12):
          warnings.warn("The real part of the broadening is non zero: {}".format(broadening))
    
        self.zpb = self.eig0.make_average(self.zpb)
    
        return self.zpb

    def get_zpr_ddw_active(self):
        """
        Compute the q-point zpr contribution in a static scheme
        with the transitions split between active and sternheimer.
        """
        self.zpr = self.get_self_energy(
            mode=False,
            temperature=False,
            omega=False,
            dynamical=False,
            only_sternheimer=False,
            only_active=True,
            only_ddw=True,
            real=True,
            )
        return self.zpr

    def get_tdr_ddw_active(self):
        """
        Compute the q-point zpr contribution in a static scheme
        with the transitions split between active and sternheimer.
        """
        self.tdr = self.get_self_energy(
            mode=False,
            temperature=True,
            omega=False,
            dynamical=False,
            only_sternheimer=False,
            only_active=True,
            only_ddw=True,
            real=True,
            shape='knt',
            )
        return self.tdr

    def get_zp_fan_active(self):
        """
        Compute only the active part of fan term
        """
        self.zpaf = self.get_self_energy(
            mode=False,
            temperature=False,
            omega=True,
            dynamical=True,
            only_sternheimer=False,
            only_active=True,
            only_fan=True,
            )
        # nkpt, nband, nomegase, nband
        self.zpaf = einsum('lkn->knl', self.zpaf) # FIXME why??

        return self.zpaf

