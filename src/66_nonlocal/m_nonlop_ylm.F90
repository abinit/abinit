!!****m* ABINIT/m_nonlop_ylm
!! NAME
!!  m_nonlop_ylm
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_nonlop_ylm

 use defs_basis
 use m_xmpi
 use m_abicore
 use m_errors

 use defs_abitypes, only : MPI_type
 use m_geometry,    only : strconv
 use m_kg,          only : ph1d3d, mkkpg
 use m_pawcprj,     only : pawcprj_type
 use m_opernla_ylm, only : opernla_ylm
 use m_opernlb_ylm, only : opernlb_ylm
 use m_opernlc_ylm, only : opernlc_ylm
 use m_opernld_ylm, only : opernld_ylm
 use m_kg,          only : mkkpgcart

 implicit none

 private
!!***

 public :: nonlop_ylm
!!***

contains
!!***

!!****f* ABINIT/nonlop_ylm
!! NAME
!! nonlop_ylm
!!
!! FUNCTION
!! * Compute application of a nonlocal operator Vnl in order to get:
!!    - contracted elements (energy, forces, stresses, ...), if signs=1
!!    - a function in reciprocal space (|out> = Vnl|in>), if signs=2
!!   Operator Vnl, as the following general form:
!!    $Vnl=sum_{R,lmn,l''m''n''} {|P_{Rlmn}> Enl^{R}_{lmn,l''m''n''} <P_{Rl''m''n''}|}$
!!   Operator Vnl is -- in the typical case -- the nonlocal potential.
!!   - With norm-conserving pseudopots, $Enl^{R}_{lmn,l''m''n''}$ is the
!!     Kleinmann-Bylander energy $Ekb^{R}_{ln}$.
!!   - In a PAW calculation, $Enl^{R}_{lmn,l''m''n''}$ are the nonlocal
!!     coefficients to connect projectors $D_{ij}$.
!!   - The |P_{Rlmn}> are the projector functions.
!! * Optionnaly, in case of PAW calculation, compute:
!!   - Application of the overlap matrix in reciprocal space
!!     (<in|S|in> or (I+S)|in>).
!!   - Application of (Vnl-lambda.S) in reciprocal space
!!     (<in|Vnl-lambda.S|in> and derivatives or (Vnl-lambda.S)|in>).
!! * This routine uses spherical harmonics Ylm to express Vnl.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  choice: chooses possible output:
!!    choice=0 => do nothing (only compute WF projected with NL projectors)
!!          =1 => non-local energy contribution
!!          =2 => 1st derivative(s) with respect to atomic position(s)
!!          =3 => 1st derivative(s) with respect to strain(s)
!!          =22=> mixed 2nd derivative(s) with respect to atomic pos. and q vector (at q=0)
!!          =25=> mixed 3rd derivative(s) with respect to atomic pos. and two q vectors (at q=0)
!!          =23=> 1st derivative(s) with respect to atomic pos. and
!!                1st derivative(s) with respect to atomic pos. and strains
!!          =4 => 2nd derivative(s) with respect to 2 atomic pos.
!!          =24=> 1st derivative(s) with respect to atm. pos. and
!!          =33=> mixed 2nd derivative(s) with respect to strain and q vector (at q=0)
!!                2nd derivative(s) with respect to 2 atomic pos.
!!          =5 => 1st derivative(s) with respect to k wavevector, typically
!!                sum_ij [ |p_i> D_ij <dp_j/dk| + |dp_i/dk> D_ij < p_j| ]
!!          =6 => 2nd derivative(s) with respect to 2 strains and
!!                mixed 2nd derivative(s) with respect to strains & atomic pos.
!!          =51 =>right 1st derivative(s) with respect to k wavevector, typically
!!                sum_ij [ |p_i> D_ij <dp_j/dk| ]
!!          =52 =>left 1st derivative(s) with respect to k wavevector, typically
!!                sum_ij [ |dp_i/dk> D_ij < p_j| ]
!!          =53 =>twist 1st derivative(s) with respect to k, typically
!!                sum_ij [ |dp_i/dk_(idir+1)> D_ij <dp_j//dk_(idir-1)|
!!                        -|dp_i/dk_(idir-1)> D_ij <dp_j//dk_(idir+1)|]
!!          =54=> mixed 2nd derivative(s) with respect to atomic pos. and left k wavevector
!!          =55=> mixed 2nd derivative(s) with respect to strain and right k wavevector
!!          =7 => apply operator $\sum_i [ |p_i> <p_i| ],
!!                same as overlap operator with s_ij=identity (paw_opt==3 only)
!!          =8 => 2nd derivatives with respect to 2 k wavevectors
!!          =81=> partial 2nd derivatives with respect to 2 k wavevectors,
!!                full derivative with respect to k1, right derivative with respect to k2,
!!                (derivative with respect to k of choice 51), typically
!!                sum_ij [ |dp_i/dk1> D_ij <dp_j/dk2| + |p_i> D_ij < d2p_j/dk1dk2| ]
!!    Only choices 1,2,3,23,4,5,6 are compatible with useylm=0.
!!    Only choices 1,2,22,25,3,5,33,51,52,53,7,8,81 are compatible with signs=2
!!  cpopt=flag defining the status of cprjin%cp(:)=<Proj_i|Cnk> scalars (see below, side effects)
!!  dimenl1,dimenl2=dimensions of enl (see enl)
!!  dimekbq=1 if enl factors do not contain a exp(-iqR) phase, 2 is they do
!!  dimffnlin=second dimension of ffnlin (1+number of derivatives)
!!  dimffnlout=second dimension of ffnlout (1+number of derivatives)
!!  enl(cplex_enl*dimenl1,dimenl2,nspinortot**2,dimekbq)=
!!  ->Norm conserving : ==== when paw_opt=0 ====
!!                      (Real) Kleinman-Bylander energies (hartree)
!!                      dimenl1=lmnmax  -  dimenl2=ntypat
!!                      dimekbq is 2 if Enl contains a exp(-iqR) phase, 1 otherwise
!!  ->PAW :             ==== when paw_opt=1, 2 or 4 ====
!!                      (Real or complex, hermitian) Dij coefs to connect projectors
!!                      dimenl1=cplex_enl*lmnmax*(lmnmax+1)/2  -  dimenl2=natom
!!                      These are complex numbers if cplex_enl=2
!!                        enl(:,:,1) contains Dij^up-up
!!                        enl(:,:,2) contains Dij^dn-dn
!!                        enl(:,:,3) contains Dij^up-dn (only if nspinor=2)
!!                        enl(:,:,4) contains Dij^dn-up (only if nspinor=2)
!!                      dimekbq is 2 if Dij contains a exp(-iqR) phase, 1 otherwise
!!  ffnlin(npwin,dimffnlin,lmnmax,ntypat)=nonlocal form factors to be used
!!          for the application of the nonlocal operator to the |in> vector
!!  ffnlout(npwout,dimffnlout,lmnmax,ntypat)=nonlocal form factors to be used
!!          for the application of the nonlocal operator to the |out> vector
!!  -----------------------------------------------------------
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2) or (choice=22,signs=2)
!!                        - k point direction in the case (choice=5,signs=2S)
!!                          for choice 53, twisted derivative involves idir+1 and idir-1
!!                        - strain component (1:6) in the case (choice=3,signs=2) or (choice=6,signs=1)
!!                        - strain component (1:9) in the case (choice=33,signs=2)
!!                        - (1:9) components to specify the atom to be moved and the second q-gradient
!!                          direction in the case (choice=25,signs=2)
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=lmn
!!  istwf_k=option parameter that describes the storage of wfs
!!  kgin(3,npwin)=integer coords of planewaves in basis sphere, for the |in> vector
!!  kgout(3,npwout)=integer coords of planewaves in basis sphere, for the |out> vector
!!  kpgin(npw,npkgin)= (k+G) components and related data, for the |in> vector
!!  kpgout(npw,nkpgout)=(k+G) components and related data, for the |out> vector
!!  kptin(3)=k point in terms of recip. translations, for the |in> vector
!!  kptout(3)=k point in terms of recip. translations, for the |out> vector
!!  lambda=factor to be used when computing (Vln-lambda.S) - only for paw_opt=2
!!         Typically lambda is the eigenvalue (or its guess)
!!  lmnmax=max. number of (l,m,n) components over all types of atoms
!!  matblk=dimension of the arrays ph3din and ph3dout
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell
!!  nattyp(ntypat)=number of atoms of each type
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpgin,nkpgout=second sizes of arrays kpgin/kpgout
!!  nloalg(3)=governs the choice of the algorithm for nonlocal operator
!!  nnlout=dimension of enlout (when signs=1 and choice>0):
!!         ==== if paw_opt=0, 1 or 2 ====
!!         choice   nnlout     |  choice   nnlout
!!              1   1          |      51   6 (complex)
!!              2   3*natom    |      52   6 (complex)
!!              3   6          |      53   6
!!              4   6*natom    |      54   9*natom
!!             23   6+3*natom  |      55   36 (complex)
!!             24   9*natom    |       6   36+18*natom
!!              5   3          |       8   6
!!                             |      81   18 (complex)
!!         ==== if paw_opt=3 ====
!!         choice   nnlout
!!              1   1
!!              2   3*natom
!!              5   3
!!             51   3
!!             52   3
!!             54   9*natom
!!             55   36
!!              7   1
!!              8   6
!!             81   9
!!         ==== if paw_opt=4 ====
!!         not available
!!  npwin=number of planewaves for given k point, for the |in> vector
!!  npwout=number of planewaves for given k point, for the |out> vector
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nspinortot=total number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in cell
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!  phkxredin(2,natom)=phase factors exp(2 pi kptin.xred)
!!  phkxredout(2,natom)=phase factors exp(2 pi kptout.xred)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1D structure factors phase information
!!  ph3din(2,npwin,matblk)=3D structure factors, for each atom and plane wave (in)
!!  ph3dout(2,npwout,matblk)=3-dim structure factors, for each atom and plane wave (out)
!!  [qdir]= optional,direction of the q-gradient (only for choice=22 choice=25 and choice=33)
!!  signs= if 1, get contracted elements (energy, forces, stress, ...)
!!         if 2, applies the non-local operator to a function in reciprocal space
!!  sij(dimenl1,ntypat*(paw_opt/3))=overlap matrix components (only if paw_opt=2, 3 or 4)
!!  ucvol=unit cell volume (bohr^3)
!!  vectin(2,npwin*nspinor)=input cmplx wavefunction coefficients <G|in>
!!  [cprjin_left(natom,nspinor)]=The projected input wave function <p_nlm|in_left>
!!    for the left wavefunction. Data are assumed to be in memory, they are NOT recalculated here.
!!
!! OUTPUT
!! ==== if (signs==1) ====
!! --If (paw_opt==0, 1 or 2)
!!    enlout(nnlout)= contribution to the non-local part of the following properties:
!!      if choice=1 : enlout(1)             -> the energy
!!      if choice=2 : enlout(3*natom)       -> 1st deriv. of energy wrt atm. pos (forces)
!!      if choice=3 : enlout(6)             -> 1st deriv. of energy wrt strain (stresses)
!!      if choice=4 : enlout(6*natom)       -> 2nd deriv. of energy wrt 2 atm. pos (dyn. mat.)
!!      if choice=23: enlout(6+3*natom)     -> 1st deriv. of energy wrt atm. pos (forces) and
!!                                             1st deriv. of energy wrt strain (stresses)
!!      if choice=24: enlout(9*natom)       -> 1st deriv. of energy wrt atm. pos (forces) and
!!                                             2nd deriv. of energy wrt 2 atm. pos (dyn. mat.)
!!      if choice=5 : enlout(3)             -> 1st deriv. of energy wrt k
!!      if choice=51: enlout(3)             -> 1st deriv. (right) of energy wrt k
!!      if choice=52: enlout(3)             -> 1st deriv. (left) of energy wrt k
!!      if choice=53: enlout(3)             -> 1st deriv. (twist) of energy wrt k
!!      if choice=54: enlout(18*natom)      -> 2nd deriv. of energy wrt atm. pos and right k (Born eff. charge)
!!      if choice=55: enlout(36)            -> 2nd deriv. of energy wrt strain and right k (piezoelastic tensor)
!!      if choice=6 : enlout(36+18*natom)   -> 2nd deriv. of energy wrt 2 strains (elast. tensor) and
!!                                             2nd deriv. of energy wrt to atm. pos and strain (internal strain)
!!      if choice=8 : enlout(6)             -> 2nd deriv. of energy wrt 2 k
!!      if choice=81: enlout(9)             -> 2nd deriv.of E: full derivative w.r.t. k1, right derivative w.r.t k2
!! --If (paw_opt==3)
!!      if choice=1 : enlout(1)             -> contribution to <c|S|c> (note: not including <c|c>)
!!      if choice=2 : enlout(3*natom)       -> contribution to <c|dS/d_atm.pos|c>
!!      if choice=51: enlout(3)             -> contribution to <c|d(right)S/d_k|c>
!!      if choice=52: enlout(3)             -> contribution to <c|d(left)S/d_k|c>
!!      if choice=54: enlout(18*natom)      -> 2nd deriv. of energy wrt atm. pos and right k (Born eff. charge)
!!      if choice=55: enlout(36)            -> 2nd deriv. of energy wrt strain and right k (piezoelastic tensor)
!!      if choice=7 : enlout(1)             -> contribution to <c|sum_i[p_i><p_i]|c>
!!      if choice=8 : enlout(6)             -> contribution to <c|d2S/d_k1d_k2|c>
!!      if choice=81: enlout(9)             -> contribution to <c|dS/d_k1[d(right)d_k2]|c>
!! --If (paw_opt==4)
!!      not available
!! ==== if (signs==2) ====
!! --if (paw_opt=0)
!!    vectout(2,npwout*my_nspinor*ndat)=result of the aplication of the concerned operator
!!                or one of its derivatives to the input vect.
!!      if (choice=22) <G|d2V_nonlocal/d(atm. pos)dq|vect_in> (at q=0)
!!      if (choice=25) <G|d3V_nonlocal/d(atm. pos)dqdq|vect_in> (at q=0)
!!      if (choice=33) <G|d2V_nonlocal/d(strain)dq|vect_in> (at q=0)
!! --if (paw_opt=0, 1 or 4)
!!    vectout(2,npwout*my_nspinor*ndat)=result of the aplication of the concerned operator
!!                or one of its derivatives to the input vect.:
!!      if (choice=1)  <G|V_nonlocal|vect_in>
!!      if (choice=2)  <G|dV_nonlocal/d(atm. pos)|vect_in>
!!      if (choice=3)  <G|dV_nonlocal/d(strain)|vect_in>
!!      if (choice=5)  <G|dV_nonlocal/d(k)|vect_in>
!!      if (choice=51) <G|d(right)V_nonlocal/d(k)|vect_in>
!!      if (choice=52) <G|d(left)V_nonlocal/d(k)|vect_in>
!!      if (choice=53) <G|d(twist)V_nonlocal/d(k)|vect_in>
!!      if (choice=8)  <G|d2V_nonlocal/d(k)d(k)|vect_in>
!!      if (choice=81) <G|d[d(right)V_nonlocal/d(k)]/d(k)|vect_in>
!! --if (paw_opt=2)
!!    vectout(2,npwout*my_nspinor*ndat)=final vector in reciprocal space:
!!      if (choice=1)  <G|V_nonlocal-lamdba.(I+S)|vect_in>
!!      if (choice=2)  <G|d[V_nonlocal-lamdba.(I+S)]/d(atm. pos)|vect_in>
!!      if (choice=3)  <G|d[V_nonlocal-lamdba.(I+S)]/d(strain)|vect_in>
!!      if (choice=5)  <G|d[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!!      if (choice=51) <G|d(right)[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!!      if (choice=52) <G|d(left)[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!!      if (choice=53) <G|d(twist)[V_nonlocal-lamdba.(I+S)]/d(k)|vect_in>
!!      if (choice=8)  <G|d2[V_nonlocal-lamdba.(I+S)]/d(k)d(k)|vect_in>
!!      if (choice=81) <G|d[d(right[V_nonlocal-lamdba.(I+S)]/d(k)]/d(k)|vect_in>
!! --if (paw_opt=3 or 4)
!!    svectout(2,npwout*my_nspinor*ndat)=result of the aplication of Sij (overlap matrix)
!!                  or one of its derivatives to the input vect.:
!!      if (choice=1)  <G|I+S|vect_in>
!!      if (choice=2)  <G|dS/d(atm. pos)|vect_in>
!!      if (choice=3)  <G|dS/d(strain)|vect_in>
!!      if (choice=5)  <G|dS/d(k)|vect_in>
!!      if (choice=51) <G|d(right)S/d(k)|vect_in>
!!      if (choice=52) <G|d(left)S/d(k)|vect_in>
!!      if (choice=53) <G|d(twist)S/d(k)|vect_in>
!!      if (choice=3)  <G|d[V_nonlocal-lamdba.(I+S)]/d(strain)|vect_in>
!!      if (choice=7)  <G|sum_i[p_i><p_i]|vect_in>
!!      if (choice=8)  <G|d2S/d(k)d(k)|vect_in>
!!      if (choice=81) <G|d[d(right)S/d(k)]/d(k)|vect_in>
!!
!! SIDE EFFECTS
!!  cprjin(natom,nspinor) <type(pawcprj_type)>=projected input wave function |in> on non-local projectors
!!                                            =<p_lmn|in> and derivatives
!!                    Treatment depends on cpopt parameter:
!!                     if cpopt=-1, <p_lmn|in> (and derivatives)
!!                                  are computed here (and not saved)
!!                     if cpopt= 0, <p_lmn|in> are computed here and saved
!!                                  derivatives are eventually computed but not saved
!!                     if cpopt= 1, <p_lmn|in> and first derivatives are computed here and saved
!!                                  other derivatives are eventually computed but not saved
!!                     if cpopt= 2  <p_lmn|in> are already in memory;
!!                                  first (and 2nd) derivatives are computed here and not saved
!!                     if cpopt= 3  <p_lmn|in> are already in memory;
!!                                  first derivatives are computed here and saved
!!                                  other derivatives are eventually computed but not saved
!!                     if cpopt= 4  <p_lmn|in> and first derivatives are already in memory;
!!                                  other derivatives are not computed (except when choice=8 or 81)
!!                                  This option is not compatible with choice=4,24 or 6
!!                     Warning: for cpopt= 1 or 3, derivatives wrt strains do not contain
!!                              the contribution due to the volume change;
!!                              i.e. <dp_lmn/dEps|in> are incomplete.
!!
!! NOTES
!! This application of the nonlocal operator is programmed using a direct
!! implementation of spherical harmonics (Ylm). Abinit used historically
!! Legendre polynomials for the application of nonlocal operator; but the
!! implementation of PAW algorithm enforced the use of Ylm.
!!
!! In the case signs=1, the array vectout is not used, nor modified
!! so that the same array as vectin can be used as a dummy argument;
!! the same is true for the pairs npwin-npwout, ffnlin-ffnlout,
!! kgin-kgout, ph3din-ph3dout, phkredin-phkxredout).
!!
!! Notes about choice==33:
!!  **Since the 2nd derivative w.r.t q-vector is calculated along cartesian
!!    directions, the 1/twopi**2 factor (that in the rest of the code is applied
!!    in the reduced to cartesian derivative conversion process) is here
!!    explicictly included in the formulas.
!!
!!  **Notice that idir=1-9, in contrast to the strain perturbation (idir=1-6),
!!    because this term is not symmetric w.r.t permutations of the two strain
!!    indices.(Also applies for choice=25)
!!
!!  **A -i factor has been factorized out in all the contributions of the second
!!    q-gradient of the metric Hamiltonian and in the first and second q-gradients
!!    of the atomic displacement Hamiltonian. This is lately included in the
!!    matrix element calculation.
!!
!! TODO
!! * Complete implementation of spin-orbit
!!
!! PARENTS
!!      m_nonlop
!!
!! CHILDREN
!!      mkkpg,mkkpgcart,opernla_ylm,opernlb_ylm,opernlc_ylm,opernld_ylm,ph1d3d
!!      strconv,xmpi_sum
!!
!! SOURCE

 subroutine nonlop_ylm(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,&
&                      enl,enlout,ffnlin,ffnlout,gprimd,idir,indlmn,istwf_k,&
&                      kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                      mpi_enreg,natom,nattyp,ngfft,nkpgin,nkpgout,nloalg,nnlout,&
&                      npwin,npwout,nspinor,nspinortot,ntypat,paw_opt,phkxredin,phkxredout,ph1d,&
&                      ph3din,ph3dout,signs,sij,svectout,ucvol,vectin,vectout,cprjin_left,qdir)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,idir
 integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,natom,nkpgin,nkpgout,nnlout
 integer,intent(in) :: npwin,npwout,nspinor,nspinortot,ntypat,paw_opt,signs
 integer,intent(in),optional :: qdir
 real(dp),intent(in) :: lambda,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),kgin(3,npwin)
 integer,intent(in),target :: indlmn(6,lmnmax,ntypat)
 integer,intent(in) :: kgout(3,npwout),nattyp(ntypat),ngfft(18),nloalg(3)
 real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2,dimekbq)
 real(dp),intent(in),target :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
 real(dp),intent(in),target :: ffnlout(npwout,dimffnlout,lmnmax,ntypat)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in),target :: kpgin(npwin,nkpgin),kpgout(npwout,nkpgout)
 real(dp),intent(in) :: kptin(3),kptout(3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: phkxredin(2,natom),phkxredout(2,natom)
 real(dp),intent(in) :: sij(dimenl1,ntypat*((paw_opt+1)/3))
 real(dp),intent(inout) :: ph3din(2,npwin,matblk),ph3dout(2,npwout,matblk)
 real(dp),intent(inout) :: vectin(:,:)
 real(dp),intent(out) :: enlout(:)
 real(dp),intent(out) :: svectout(:,:)
 real(dp),intent(inout) :: vectout (:,:)
 type(pawcprj_type),intent(inout) :: cprjin(:,:)
 type(pawcprj_type),optional,intent(in) :: cprjin_left(:,:)

!Local variables-------------------------------
!scalars
 integer :: choice_a,choice_b,cplex,cplex_enl,cplex_fac,ia,ia1,ia2,ia3,ia4,ia5
 integer :: iatm,ic,idir1,idir2,ii,ierr,ilmn,ishift,ispinor,itypat,jc,mincat,mu,mua,mub,mu0
 integer :: n1,n2,n3,nd2gxdt,ndgxdt,ndgxdt_stored,nd2gxdtfac,ndgxdtfac
 integer :: nincat,nkpgin_,nkpgout_,nlmn,nu,nua1,nua2,nub1,nub2,optder
 real(dp) :: enlk
 logical :: check,testnl
 character(len=500) :: message
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer,parameter :: gamma(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))
 integer,allocatable :: cplex_dgxdt(:),cplex_d2gxdt(:)
 integer,ABI_CONTIGUOUS pointer :: indlmn_typ(:,:)
 real(dp),allocatable :: d2gxdt(:,:,:,:,:),d2gxdtfac(:,:,:,:,:),d2gxdtfac_sij(:,:,:,:,:)
 real(dp),allocatable :: dgxdt(:,:,:,:,:),dgxdtfac(:,:,:,:,:),dgxdtfac_sij(:,:,:,:,:)
 real(dp),allocatable :: ddkk(:),fnlk(:),gmet(:,:)
 real(dp),allocatable :: gx(:,:,:,:),gxfac(:,:,:,:),gxfac_sij(:,:,:,:),gx_left(:,:,:,:)
 real(dp),allocatable :: sij_typ(:),strnlk(:)
 real(dp),allocatable :: work1(:),work2(:),work3(:,:),work4(:,:),work5(:,:,:),work6(:,:,:),work7(:,:,:)
 real(dp),ABI_CONTIGUOUS pointer :: ffnlin_typ(:,:,:),ffnlout_typ(:,:,:),kpgin_(:,:),kpgout_(:,:)

! **********************************************************************

 DBG_ENTER("COLL")

!Check consistency of arguments
!==============================================================

!signs=1, almost all choices
 if (signs==1) then
   if(paw_opt<3) then
     check=(choice==0 .or.choice==1 .or.choice==2 .or.choice==3 .or.choice==4 .or.&
&     choice==23.or.choice==24.or.choice==5 .or.choice==51.or.choice==52.or.&
&     choice==53.or.choice==54.or.choice==55.or.&
&     choice==6 .or.choice==8 .or.choice==81)
   else if (paw_opt==3) then
     check=(choice== 0.or.choice== 1.or.choice== 2.or.choice==3.or.choice==5.or.&
&     choice==23.or.choice==51.or.choice==52.or.choice==54.or.choice==55.or.&
&     choice== 8.or.choice==81)
   else
     check = .false.
   end if
   ABI_CHECK(check,'BUG: choice not compatible (for signs=1)')
 end if

!signs=2, less choices
 if (signs==2) then
   check=(choice==0.or.choice==1.or.choice==2.or.choice==22.or.choice==25.or.choice==3 .or.&
&   choice==5.or.choice==33.or.choice==51.or.choice==52.or.choice==53.or.choice==54.or.&
&   choice==7.or.choice==8.or.choice==81)
   ABI_CHECK(check,'BUG: choice not compatible (for signs=2)')
 end if
!1<=idir<=6 is required when choice=3 and signs=2
 if (choice==3.and.signs==2) then
   check=(idir>=1.and.idir<=6)
   ABI_CHECK(check,'BUG: choice=3 and signs=2 requires 1<=idir<=6')
!1<=idir<=9 is required when choice= 25 or 33 and signs=2
 else if ((choice==25.or.choice==33).and.signs==2) then
   check=(idir>=1.and.idir<=9)
   ABI_CHECK(check,'BUG: choice= 25 or 33 and signs=2 requires 1<=idir<=9')
!1<=idir<=9 is required when choice==8/81 and signs=2
 else if ((choice==8.or.choice==81.or.choice==54).and.signs==2) then
   check=(idir>=1.and.idir<=9)
   ABI_CHECK(check,'BUG: choice=8/81 and signs=2 requires 1<=idir<=9')
 else
!  signs=2 requires 1<=idir<=3 when choice>1
   check=(signs/=2.or.choice<=1.or.choice==7.or.(idir>=1.and.idir<=3))
   ABI_CHECK(check,'BUG: signs=2 requires 1<=idir<=3')
 end if
!1<=qdir<=3 is required when choice==22 or choice==25 or choice==33 and signs=2
 if ((choice==22.or.choice==25.or.choice==33).and.signs==2) then
   check=(qdir>=1.and.qdir<=3)
   ABI_CHECK(check,'BUG: choice=22,25 or 33 and signs=2 requires 1<=qdir<=3')
 end if

!check allowed values for cpopt
 check=(cpopt>=-1.and.cpopt<=4)
 ABI_CHECK(check,'bad value for cpopt')
 check=(cpopt/=4.or.(choice/=4.and.choice/=24.and.choice/=6))
 ABI_CHECK(check,'BUG: cpopt=4 not allowed for 2nd derivatives')
 check=(cpopt/=2.or.(choice/=8.and.choice/=81))
 ABI_CHECK(check,'BUG: cpopt=2 not allowed for choice=8,81, use cpopt=4 instead')
!check conditions for optional arguments
 check=((.not.present(cprjin_left)).or.(signs==1.and.choice==1))
 ABI_CHECK(check,'BUG: when cprjin_left is present, must have choice=1,signs=1')
!protect special case choice==7
 check=(choice/=7.or.paw_opt==3)
 ABI_CHECK(check,'BUG: when choice=7, paw_opt must be 3')
!spin-orbit not yet allowed
 check=(maxval(indlmn(6,:,:))<=1)
 ABI_CHECK(check,'BUG: spin-orbit with Yml for nonlop not yet allowed')

!Test: size of blocks of atoms
 mincat=min(NLO_MINCAT,maxval(nattyp))
 if (nloalg(2)<=0.and.mincat>matblk) then
   write(message, '(a,a,a,i4,a,i4,a)' ) &
&   'With nloc_mem<=0, mincat must be less than matblk.',ch10,&
&   'Their value is ',mincat,' and ',matblk,'.'
   MSG_BUG(message)
 end if

!Define dimensions of projected scalars
!==============================================================

!Define some useful variables
 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3)
 choice_a=merge(choice,1,choice/=7);choice_b=choice_a
 if (cpopt>=2) choice_a=-choice_a
 cplex=2;if (istwf_k>1) cplex=1 !Take into account TR-symmetry
 cplex_enl=1;if (paw_opt>0) cplex_enl=2*dimenl1/(lmnmax*(lmnmax+1))
 cplex_fac=max(cplex,dimekbq)
 if ((nspinortot==2.or.cplex_enl==2).and.paw_opt>0.and.choice/=7) cplex_fac=2

!Define dimensions of projected scalars
 ndgxdt=0;ndgxdtfac=0;nd2gxdt=0;nd2gxdtfac=0
 if (choice==2) then
   if (signs==1) ndgxdt=3
   if (signs==2) ndgxdt=1
   if (signs==2) ndgxdtfac=1
 end if
 if (choice==22) then
   if (signs==2) ndgxdt=1
   if (signs==2) ndgxdtfac=1
 end if
 if (choice==25) then
   if (signs==2) ndgxdt=1
   if (signs==2) ndgxdtfac=1
 end if
 if (choice==3) then
   if (signs==1) ndgxdt=6
   if (signs==2) ndgxdt=1
   if (signs==2) ndgxdtfac=1
 end if
 if (choice==23) then
   if (signs==1) ndgxdt=9
 end if
 if (choice==4) then
   if(signs==1) ndgxdt=3
   if(signs==1) ndgxdtfac=3
   if(signs==1) nd2gxdt=6
 end if
 if (choice==24) then
   if(signs==1) ndgxdt=3
   if(signs==1) ndgxdtfac=3
   if(signs==1) nd2gxdt=6
 end if
 if (choice==5) then
   if(signs==1) ndgxdt=3
   if(signs==2) ndgxdt=1
   if(signs==2) ndgxdtfac=1
 end if
 if (choice==33) then
   if(signs==2) ndgxdt=2
   if(signs==2) ndgxdtfac=2
   if(signs==2) nd2gxdt=3
   if(signs==2) nd2gxdtfac=3
 end if
 if (choice==51) then
   if(signs==1) ndgxdt=3
   if(signs==2) ndgxdt=1
   if(signs==2) ndgxdtfac=1
 end if
 if (choice==52) then
   if(signs==1) ndgxdt=3
   if(signs==2) ndgxdt=1
   if(signs==2) ndgxdtfac=1
 end if
 if (choice==53) then
   if(signs==1) ndgxdt=3
   if(signs==1) ndgxdtfac=3
   if(signs==2) ndgxdt=2
   if(signs==2) ndgxdtfac=2
 end if
 if (choice==54) then
   if(signs==1) ndgxdt=6
   if(signs==1) ndgxdtfac=6
   if(signs==1) nd2gxdt=9
   if(signs==2) ndgxdt=1
   if(signs==2) nd2gxdt=1
   if(signs==2) ndgxdtfac=1
   if(signs==2) nd2gxdtfac=1
 end if
 if (choice==55) then
   if(signs==1) ndgxdt=9
   if(signs==1) ndgxdtfac=9
   if(signs==1) nd2gxdt=18
 end if
 if (choice==6) then
   if(signs==1) ndgxdt=9
   if(signs==1) ndgxdtfac=9
   if(signs==1) nd2gxdt=54
 end if
 if (choice==8) then
   if(signs==1) ndgxdt=3
   if(signs==1) ndgxdtfac=3
   if(signs==1) nd2gxdt=6
   if(signs==2) ndgxdt=2
   if(signs==2) ndgxdtfac=2
   if(signs==2) nd2gxdt=1
   if(signs==2) nd2gxdtfac=1
 end if
 if (choice==81) then
   if(signs==1) ndgxdt=3
   if(signs==1) ndgxdtfac=3
   if(signs==1) nd2gxdt=6
   if(signs==2) ndgxdt=1
   if(signs==2) ndgxdtfac=1
   if(signs==2) nd2gxdt=1
   if(signs==2) nd2gxdtfac=1
 end if
 ABI_CHECK(ndgxdtfac<=ndgxdt,"BUG: ndgxdtfac>ndgxdt!")
 optder=0;if (ndgxdtfac>0) optder=1
 if (nd2gxdtfac>0) optder=2

!Consistency tests
 if (cpopt==4) then
   if (ndgxdt>0.and.cprjin(1,1)%ncpgr<=0) then
     message='cprjin%ncpgr=0 not allowed with cpopt=4 and these (choice,signs) !'
     MSG_BUG(message)
   end if
 end if
 if (cpopt==1.or.cpopt==3) then
   if (cprjin(1,1)%ncpgr<ndgxdt) then
     message='should have cprjin%ncpgr>=ndgxdt with cpopt=1 or 3 !'
     MSG_BUG(message)
   end if
 end if

!Additional steps before calculation
!==============================================================

!Initialize output arrays
 if (signs==1) then
   ABI_ALLOCATE(fnlk,(3*natom))
   ABI_ALLOCATE(ddkk,(6))
   ABI_ALLOCATE(strnlk,(6))
   enlk=zero;fnlk=zero;ddkk=zero;strnlk=zero
   enlout(:)=zero
 end if
 if (signs==2) then
   if (paw_opt==0.or.paw_opt==1.or.paw_opt==4) vectout(:,:)=zero
   if (paw_opt==2.and.choice==1) vectout(:,:)=-lambda*vectin(:,:)
   if (paw_opt==2.and.choice> 1) vectout(:,:)=zero
   if (paw_opt==3.or.paw_opt==4) then
     if (choice==1) svectout(:,:)=vectin(:,:)
     if (choice> 1) svectout(:,:)=zero
   end if
 end if

!Eventually re-compute (k+G) vectors (and related data)
 nkpgin_=0
 if (choice==2.or.choice==22.or.choice==25.or.choice==33.or.choice==54) nkpgin_=3
 if (signs==1) then
   if (choice==4.or.choice==24) nkpgin_=9
   if (choice==3.or.choice==23.or.choice==6) nkpgin_=3
   if (choice==55) nkpgin_=3
 end if
 if (nkpgin<nkpgin_) then
   ABI_ALLOCATE(kpgin_,(npwin,nkpgin_))

   !For the metric derivatives we need kpg in Cartesian coordinates
   if (choice==33) then
     call mkkpgcart(gprimd,kgin,kpgin_,kptin,nkpgin_,npwin)
   else
     call mkkpg(kgin,kpgin_,kptin,nkpgin_,npwin)
   end if

 else
   nkpgin_ = nkpgin
   kpgin_  => kpgin
 end if

 nkpgout_=0
 if ((choice==2.or.choice==22.or.choice==25.or.choice==3.or.choice==33.or.choice==54).and.signs==2) nkpgout_=3
 if (nkpgout<nkpgout_) then
   ABI_ALLOCATE(kpgout_,(npwout,nkpgout_))

   !For the metric derivatives we need kpg in Cartesian coordinates
   if (choice==33) then
     call mkkpgcart(gprimd,kgout,kpgout_,kptout,nkpgout_,npwout)
   else
     call mkkpg(kgout,kpgout_,kptout,nkpgout_,npwout)
   end if

 else
   nkpgout_ = nkpgout
   kpgout_ => kpgout
 end if

!Big loop on atom types
!==============================================================

 ia1=1;iatm=0
 do itypat=1,ntypat

!  Get atom loop indices for different types:
   ia2=ia1+nattyp(itypat)-1;ia5=1

!  Select quantities specific to the current type of atom
   nlmn=count(indlmn(3,:,itypat)>0)

!  Test on local part
   testnl=(paw_opt/=0)
   if (paw_opt==0) testnl=any(enl(:,:,:,:)>tol10)

!  Some non-local part is to be applied for that type of atom
   if (testnl) then

!    Store some quantities depending only of the atom type
     ffnlin_typ => ffnlin(:,:,:,itypat)
     indlmn_typ => indlmn(:,:,itypat)
     if (signs==2) then
       ffnlout_typ => ffnlout(:,:,:,itypat)
     end if
     if (paw_opt>=2) then
       ABI_ALLOCATE(sij_typ,(nlmn*(nlmn+1)/2))
       if (cplex_enl==1) then
         do ilmn=1,nlmn*(nlmn+1)/2
           sij_typ(ilmn)=sij(ilmn,itypat)
         end do
       else
         do ilmn=1,nlmn*(nlmn+1)/2
           sij_typ(ilmn)=sij(2*ilmn-1,itypat)
         end do
       end if
     else
       ABI_ALLOCATE(sij_typ,(0))
     end if

!    Loop over atoms of the same type
!    ==============================================================

!    Cut the sum on different atoms in blocks, to allow memory saving.
!    Inner summations on atoms will be done from ia3 to ia4.
!    Note: the maximum range from ia3 to ia4 is mincat (max. increment of atoms).

     do ia3=ia1,ia2,mincat
       ia4=min(ia2,ia3+mincat-1)
!      Give the increment of number of atoms in this subset.
       nincat=ia4-ia3+1

!      Prepare the phase factors if they were not already computed
       if (nloalg(2)<=0) then
         call ph1d3d(ia3,ia4,kgin,matblk,natom,npwin,n1,n2,n3,phkxredin,ph1d,ph3din)
       end if

!      Allocate memory for projected scalars
       ABI_ALLOCATE(gx,(cplex,nlmn,nincat,nspinor))
       ABI_ALLOCATE(dgxdt,(cplex,ndgxdt,nlmn,nincat,nspinor))
       ABI_ALLOCATE(d2gxdt,(cplex,nd2gxdt,nlmn,nincat,nspinor))
       ABI_ALLOCATE(d2gxdtfac,(cplex_fac,nd2gxdtfac,nlmn,nincat,nspinor))
       ABI_ALLOCATE(dgxdtfac,(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor))
       ABI_ALLOCATE(gxfac,(cplex_fac,nlmn,nincat,nspinor))
       gx(:,:,:,:)=zero;gxfac(:,:,:,:)=zero
       if (ndgxdt>0) dgxdt(:,:,:,:,:)=zero
       if (ndgxdtfac>0) dgxdtfac(:,:,:,:,:)=zero
       if (nd2gxdt>0) d2gxdt(:,:,:,:,:)=zero
       if (nd2gxdtfac>0) d2gxdtfac(:,:,:,:,:)=zero
       if (paw_opt>=3) then
         ABI_ALLOCATE(gxfac_sij,(cplex,nlmn,nincat,nspinor))
         ABI_ALLOCATE(dgxdtfac_sij,(cplex,ndgxdtfac,nlmn,nincat,nspinor))
         ABI_ALLOCATE(d2gxdtfac_sij,(cplex,nd2gxdtfac,nlmn,nincat,nspinor))
         gxfac_sij(:,:,:,:)=zero
         if (ndgxdtfac>0) dgxdtfac_sij(:,:,:,:,:)=zero
         if (nd2gxdtfac>0) d2gxdtfac_sij(:,:,:,:,:) = zero
       else
         ABI_ALLOCATE(gxfac_sij,(0,0,0,0))
         ABI_ALLOCATE(dgxdtfac_sij,(0,0,0,0,0))
         ABI_ALLOCATE(d2gxdtfac_sij,(0,0,0,0,0))
       end if

!      When istwf_k > 1, gx derivatives can be real or pure imaginary
!      cplex_dgxdt(i)  = 1 if dgxdt(1,i,:,:)  is real, 2 if it is pure imaginary
!      cplex_d2gxdt(i) = 1 if d2gxdt(1,i,:,:) is real, 2 if it is pure imaginary
       ABI_ALLOCATE(cplex_dgxdt,(ndgxdt))
       ABI_ALLOCATE(cplex_d2gxdt,(nd2gxdt))
       cplex_dgxdt(:) = 1 ; cplex_d2gxdt(:) = 1
       if(ndgxdt > 0) then
         if (choice==5.or.choice==51.or.choice==52.or.choice==53.or. &
&         choice==8.or.choice==81) cplex_dgxdt(:) = 2
         if (choice==54.and.signs==1) cplex_dgxdt(4:6) = 2
         if (choice==54.and.signs==2) cplex_dgxdt(:)   = 2
         if (choice==55.and.signs==1) cplex_dgxdt(7:9) = 2
       end if
       if(nd2gxdt > 0) then
         if (choice==54) cplex_d2gxdt(:) = 2
         if (choice==55.and.signs==1) cplex_d2gxdt(1:18)= 2
       end if

!      Compute projection of current wave function |c> on each
!      non-local projector: <p_lmn|c>
!      ==============================================================

!      Retrieve eventually <p_lmn|c> coeffs (and derivatives)
       if (cpopt>=2) then
         do ispinor=1,nspinor
           do ia=1,nincat
             gx(1:cplex,1:nlmn,ia,ispinor)=cprjin(iatm+ia,ispinor)%cp(1:cplex,1:nlmn)
           end do
         end do
       end if
       if (cpopt==4.and.ndgxdt>0) then
         ndgxdt_stored = cprjin(1,1)%ncpgr
         ishift=0
         if (((choice==2).or.(choice==3)).and.(ndgxdt_stored>ndgxdt).and.(signs==2)) ishift=idir-ndgxdt
         if ((choice==2).and.(ndgxdt_stored==9).and.(signs==2)) ishift=ishift+6
         if (choice==2.and.(ndgxdt_stored>ndgxdt).and.(signs==1)) ishift=ndgxdt_stored-ndgxdt
         if(cplex == 2) then
           do ispinor=1,nspinor
             do ia=1,nincat
               if (ndgxdt_stored==ndgxdt.or.(ndgxdt_stored>ndgxdt.and.((choice==2).or.(choice==3)))) then
                 dgxdt(1:2,1:ndgxdt,1:nlmn,ia,ispinor)=cprjin(iatm+ia,ispinor)%dcp(1:2,1+ishift:ndgxdt+ishift,1:nlmn)
               else if (signs==2.and.ndgxdt_stored==3) then
                 if (choice==5.or.choice==51.or.choice==52) then ! ndgxdt=1
                   dgxdt(1:2,1,1:nlmn,ia,ispinor)=cprjin(iatm+ia,ispinor)%dcp(1:2,idir,1:nlmn)
                 else if (choice==8) then ! ndgxdt=2
                   idir1=(idir-1)/3+1; idir2=mod((idir-1),3)+1
                   dgxdt(1:2,1,1:nlmn,ia,ispinor)=cprjin(iatm+ia,ispinor)%dcp(1:2,idir1,1:nlmn)
                   dgxdt(1:2,2,1:nlmn,ia,ispinor)=cprjin(iatm+ia,ispinor)%dcp(1:2,idir2,1:nlmn)
                 else if (choice==81) then ! ndgxdt=1
                   idir1=(idir-1)/3+1; idir2=mod((idir-1),3)+1
                   dgxdt(1:2,1,1:nlmn,ia,ispinor)=cprjin(iatm+ia,ispinor)%dcp(1:2,idir2,1:nlmn)
                 end if
               end if
             end do
           end do
         else
           do ispinor=1,nspinor
             do ia=1,nincat
               do ilmn=1,nlmn
                 if (ndgxdt_stored==ndgxdt.or.(ndgxdt_stored>ndgxdt.and.((choice==2).or.(choice==3)))) then
                   do ii=1,ndgxdt
                     ic = cplex_dgxdt(ii)
                     dgxdt(1,ii,ilmn,ia,ispinor)=cprjin(iatm+ia,ispinor)%dcp(ic,ii+ishift,ilmn)
                   end do
                 else if (signs==2.and.ndgxdt_stored==3) then
                   if (choice==5.or.choice==51.or.choice==52) then ! ndgxdt=1
                     dgxdt(1,1,ilmn,ia,ispinor)=cprjin(iatm+ia,ispinor)%dcp(cplex_dgxdt(1),idir,ilmn)
                   else if (choice==8) then ! ndgxdt=2
                     idir1=(idir-1)/3+1; idir2=mod((idir-1),3)+1
                     dgxdt(1,1,ilmn,ia,ispinor)=cprjin(iatm+ia,ispinor)%dcp(cplex_dgxdt(1),idir1,ilmn)
                     dgxdt(1,2,ilmn,ia,ispinor)=cprjin(iatm+ia,ispinor)%dcp(cplex_dgxdt(2),idir2,ilmn)
                   else if (choice==81) then ! ndgxdt=1
                     idir1=(idir-1)/3+1; idir2=mod((idir-1),3)+1
                     dgxdt(1,1,ilmn,ia,ispinor)=cprjin(iatm+ia,ispinor)%dcp(cplex_dgxdt(1),idir2,ilmn)
                   end if
                 end if
               end do
             end do
           end do
         end if
       end if

!      Computation or <p_lmn|c> (and derivatives) for this block of atoms
       if ((cpopt<4.and.choice_a/=-1).or.choice==8.or.choice==81) then
         call opernla_ylm(choice_a,cplex,cplex_dgxdt,cplex_d2gxdt,dimffnlin,d2gxdt,dgxdt,ffnlin_typ,gx,&
&         ia3,idir,indlmn_typ,istwf_k,kpgin_,matblk,mpi_enreg,nd2gxdt,ndgxdt,nincat,nkpgin_,nlmn,&
&         nloalg,npwin,nspinor,ph3din,signs,ucvol,vectin,qdir=qdir)
       end if

!      Transfer result to output variable cprj (if requested)
!      cprj(:)%cp receive the <p_i|Psi> factors (p_i: non-local projector)
!      Be careful: cprj(:)%dcp does not exactly contain the derivative of cprj(:)%cp.
!                  - Volume contributions (in case of strain derivative) are not included,
!                  - Global coordinate transformation are not applied,
!                  cprj(:)%dcp is meant to be a restart argument of the present nonlop routine.
       if (cpopt==0.or.cpopt==1) then
         do ispinor=1,nspinor
           do ia=1,nincat
             cprjin(iatm+ia,ispinor)%nlmn=nlmn
             cprjin(iatm+ia,ispinor)%cp(1:cplex,1:nlmn)=gx(1:cplex,1:nlmn,ia,ispinor)
             if (cplex==1) cprjin(iatm+ia,ispinor)%cp(2,1:nlmn)=zero
           end do
         end do
       end if
       if ((cpopt==1.or.cpopt==3).and.ndgxdt>0) then
         ishift=0
         if ((choice==2).and.(cprjin(1,1)%ncpgr>ndgxdt)) ishift=cprjin(1,1)%ncpgr-ndgxdt
         if(cplex==2)then
           do ispinor=1,nspinor
             do ia=1,nincat
               cprjin(iatm+ia,ispinor)%dcp(1:2,1+ishift:ndgxdt+ishift,1:nlmn)=dgxdt(1:2,1:ndgxdt,1:nlmn,ia,ispinor)
!               cprjin(iatm+ia,ispinor)%dcp(1:2,1:ndgxdt,1:nlmn)=dgxdt(1:2,1:ndgxdt,1:nlmn,ia,ispinor)
             end do
           end do
         else
           do ispinor=1,nspinor
             do ia=1,nincat
               do ilmn=1,nlmn
                 do ii=1,ndgxdt
                   ic = cplex_dgxdt(ii) ; jc = 3 - ic
                   cprjin(iatm+ia,ispinor)%dcp(ic,ii+ishift,ilmn)=dgxdt(1,ii,ilmn,ia,ispinor)
                   cprjin(iatm+ia,ispinor)%dcp(jc,ii+ishift,ilmn)=zero
!                   cprjin(iatm+ia,ispinor)%dcp(ic,ii,ilmn)=dgxdt(1,ii,ilmn,ia,ispinor)
!                   cprjin(iatm+ia,ispinor)%dcp(jc,ii,ilmn)=zero
                 end do
               end do
             end do
           end do
         end if
       end if

!      If choice==0, that's all for these atoms !
       if (choice>0) then
         if(choice/=7) then
!          Contraction from <p_i|c> to Sum_j[Dij.<p_j|c>] (and derivatives)
           call opernlc_ylm(atindx1,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_enl,cplex_fac,dgxdt,dgxdtfac,dgxdtfac_sij,&
&           d2gxdt,d2gxdtfac,d2gxdtfac_sij,dimenl1,dimenl2,dimekbq,enl,gx,gxfac,gxfac_sij,&
&           iatm,indlmn_typ,itypat,lambda,mpi_enreg,natom,ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,&
&           nincat,nlmn,nspinor,nspinortot,optder,paw_opt,sij_typ)
         else
           gxfac_sij=gx
         end if

!        Operate with the non-local potential on the projected scalars,
!        in order to get contributions to energy/forces/stress/dyn.mat
!        ==============================================================
         if (signs==1) then
           if (.not.present(cprjin_left)) then
             call opernld_ylm(choice_b,cplex,cplex_fac,ddkk,dgxdt,dgxdtfac,dgxdtfac_sij,d2gxdt,&
&             enlk,enlout,fnlk,gx,gxfac,gxfac_sij,ia3,natom,nd2gxdt,ndgxdt,ndgxdtfac,&
&             nincat,nlmn,nnlout,nspinor,paw_opt,strnlk)
           else
             ABI_ALLOCATE(gx_left,(cplex,nlmn,nincat,nspinor))
!            Retrieve <p_lmn|c> coeffs
             do ispinor=1,nspinor
               do ia=1,nincat
                 gx_left(1:cplex,1:nlmn,ia,ispinor)=cprjin_left(iatm+ia,ispinor)%cp(1:cplex,1:nlmn)
               end do
             end do
!            TODO
!            if (cpopt==4.and.ndgxdt>0) then
!            do ispinor=1,nspinor
!            do ia=1,nincat
!            dgxdt_left(1:cplex,1:ndgxdt,1:nlmn,ia,ispinor)=cprjin_left(iatm+ia,ispinor)%dcp(1:cplex,1:ndgxdt,1:nlmn)
!            end do
!            end do
!            end if
             call opernld_ylm(choice_b,cplex,cplex_fac,ddkk,dgxdt,dgxdtfac,dgxdtfac_sij,d2gxdt,&
&             enlk,enlout,fnlk,gx_left,gxfac,gxfac_sij,ia3,natom,nd2gxdt,ndgxdt,ndgxdtfac,&
&             nincat,nlmn,nnlout,nspinor,paw_opt,strnlk)
             ABI_DEALLOCATE(gx_left)
           end if
         end if

!        Operate with the non-local potential on the projected scalars,
!        in order to get matrix element
!        ==============================================================
         if (signs==2) then
!          Prepare the phase factors if they were not already computed
           if(nloalg(2)<=0) then
             call ph1d3d(ia3,ia4,kgout,matblk,natom,npwout,n1,n2,n3,phkxredout,ph1d,ph3dout)
           end if
           call opernlb_ylm(choice_b,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_fac,&
&           d2gxdtfac,d2gxdtfac_sij,dgxdtfac,dgxdtfac_sij,dimffnlout,ffnlout_typ,gxfac,gxfac_sij,ia3,&
&           idir,indlmn_typ,kpgout_,matblk,ndgxdtfac,nd2gxdtfac,nincat,nkpgout_,nlmn,&
&           nloalg,npwout,nspinor,paw_opt,ph3dout,svectout,ucvol,vectout,qdir=qdir)
         end if

       end if ! choice==0

!      Deallocate temporary projected scalars
       ABI_DEALLOCATE(gx)
       ABI_DEALLOCATE(gxfac)
       ABI_DEALLOCATE(dgxdt)
       ABI_DEALLOCATE(dgxdtfac)
       ABI_DEALLOCATE(d2gxdt)
       ABI_DEALLOCATE(d2gxdtfac)
       ABI_DEALLOCATE(dgxdtfac_sij)
       ABI_DEALLOCATE(d2gxdtfac_sij)
       ABI_DEALLOCATE(gxfac_sij)
       ABI_DEALLOCATE(cplex_dgxdt)
       ABI_DEALLOCATE(cplex_d2gxdt)

!      End sum on atom subset loop
       iatm=iatm+nincat;ia5=ia5+nincat
     end do
     !if (paw_opt>=2)  then
     ABI_DEALLOCATE(sij_typ)
     !end if

!    End condition of existence of a non-local part
   else
     if (cpopt==0.or.cpopt==1) then
       do ispinor=1,nspinor
         do ia=1,nattyp(itypat)
           cprjin(iatm+ia,ispinor)%cp(:,1:nlmn)=zero
         end do
       end do
     end if
     if ((cpopt==1.or.cpopt==3).and.ndgxdt>0) then
       ishift=0
       if ((choice==2).and.(cprjin(1,1)%ncpgr>ndgxdt)) ishift=cprjin(1,1)%ncpgr-ndgxdt
       do ispinor=1,nspinor
         do ia=1,nattyp(itypat)
           cprjin(iatm+ia,ispinor)%dcp(:,1+ishift:ndgxdt+ishift,1:nlmn)=zero
!           cprjin(iatm+ia,ispinor)%dcp(:,1:ndgxdt,1:nlmn)=zero
         end do
       end do
     end if
     iatm=iatm+nattyp(itypat)
   end if

!  End atom type loop
   ia1=ia2+1
 end do


!Reduction in case of parallelism
!==============================================================

 if (signs==1.and.mpi_enreg%paral_spinor==1) then
   if (nnlout/=0) then
     call xmpi_sum(enlout,mpi_enreg%comm_spinor,ierr)
   end if
   if (choice==3.or.choice==6.or.choice==23) then
     call xmpi_sum(enlk,mpi_enreg%comm_spinor,ierr)
   end if
   if (choice==6) then
     call xmpi_sum(fnlk,mpi_enreg%comm_spinor,ierr)
     call xmpi_sum(strnlk,mpi_enreg%comm_spinor,ierr)
   end if
   if (choice==55) then
     call xmpi_sum(ddkk,mpi_enreg%comm_spinor,ierr)
   end if
 end if

!Coordinate transformations
!==============================================================

!Need sometimes gmet
 if ((signs==1.and.paw_opt<=3).and. &
& (choice==5 .or.choice==51.or.choice==52.or.choice==53.or.&
& choice==54.or.choice==55)) then
   ABI_ALLOCATE(gmet,(3,3))
   gmet = MATMUL(TRANSPOSE(gprimd),gprimd)
 end if

!1st derivative wrt to strain (stress tensor):
! - convert from reduced to cartesian coordinates
! - substract volume contribution
 if ((choice==3.or.choice==23).and.signs==1.and.paw_opt<=3) then
   mu0=0 ! Shift to be applied in enlout array
   ABI_ALLOCATE(work1,(6))
   work1(1:6)=enlout(mu0+1:mu0+6)
   call strconv(work1,gprimd,work1)
   enlout(mu0+1:mu0+3)=(work1(1:3)-enlk)
   enlout(mu0+4:mu0+6)= work1(4:6)
   ABI_DEALLOCATE(work1)
 end if

!1st derivative wrt to k wave vector (ddk):
! - convert from cartesian to reduced coordinates
 if ((choice==5.or.choice==53).and.signs==1.and.paw_opt<=3) then
   mu0=0 ! Shift to be applied in enlout array
   ABI_ALLOCATE(work1,(3))
   work1(:)=enlout(mu0+1:mu0+3)
   enlout(mu0+1:mu0+3)=gmet(:,1)*work1(1)+gmet(:,2)*work1(2)+gmet(:,3)*work1(3)
   ABI_DEALLOCATE(work1)
 end if
 if ((choice==51.or.choice==52).and.signs==1.and.paw_opt<=3) then
   mu0=0 ! Shift to be applied in enlout array
   ABI_ALLOCATE(work1,(3))
   do mu=1,2 ! Loop for Re,Im
     work1(1:3)=(/enlout(mu0+1),enlout(mu0+3),enlout(mu0+5)/)
     enlout(mu0+1)=gmet(1,1)*work1(1)+gmet(1,2)*work1(2)+gmet(1,3)*work1(3)
     enlout(mu0+3)=gmet(2,1)*work1(1)+gmet(2,2)*work1(2)+gmet(2,3)*work1(3)
     enlout(mu0+5)=gmet(3,1)*work1(1)+gmet(3,2)*work1(2)+gmet(3,3)*work1(3)
     mu0=mu0+1
   end do
   ABI_DEALLOCATE(work1)
 end if

!2nd derivative wrt to k wave vector and atomic position (effective charges):
! - convert from cartesian to reduced coordinates
 if (choice==54.and.signs==1.and.paw_opt<=3) then
   mu0=0 ! Shift to be applied in enlout array
   ABI_ALLOCATE(work1,(3))
   ABI_ALLOCATE(work2,(3))
   do mu=1,3*natom
!    First, real part
     work1(1)=enlout(mu0+1);work1(2)=enlout(mu0+3);work1(3)=enlout(mu0+5)
     work2(:)=gmet(:,1)*work1(1)+gmet(:,2)*work1(2)+gmet(:,3)*work1(3)
     enlout(mu0+1)=work2(1);enlout(mu0+3)=work2(2);enlout(mu0+5)=work2(3)
!    Then imaginary part
     work1(1)=enlout(mu0+2);work1(2)=enlout(mu0+4);work1(3)=enlout(mu0+6)
     work2(:)=gmet(:,1)*work1(1)+gmet(:,2)*work1(2)+gmet(:,3)*work1(3)
     enlout(mu0+2)=work2(1);enlout(mu0+4)=work2(2);enlout(mu0+6)=work2(3)
     mu0=mu0+6
   end do
   ABI_DEALLOCATE(work1)
   ABI_DEALLOCATE(work2)
 end if

!2nd derivative wrt to k wave vector and strain (piezoelectric tensor):
! - convert from cartesian to reduced coordinates (k point)
! - convert from reduced to cartesian coordinates (strain)
! - substract volume contribution
! - symetrize strain components
 if (choice==55.and.signs==1.and.paw_opt<=3) then
   ABI_ALLOCATE(work3,(2,3))
   ABI_ALLOCATE(work4,(2,3))
   ABI_ALLOCATE(work5,(2,3,6))
   ABI_ALLOCATE(work7,(2,3,6))
   ABI_ALLOCATE(work6,(2,3,3))
   do ic=1,3 ! gamma
     work5=zero
     do jc=1,3 ! nu
       do ii=1,3 ! lambda
         mu=(gamma(jc,ii)-1)*3+1
         work5(1,jc,ii)=gmet(ic,1)*enlout(2*mu-1)+gmet(ic,2)*enlout(2*mu+1) &
&         +gmet(ic,3)*enlout(2*mu+3)
         work5(2,jc,ii)=gmet(ic,1)*enlout(2*mu  )+gmet(ic,2)*enlout(2*mu+2) &
&         +gmet(ic,3)*enlout(2*mu+4)
       end do
     end do
     work6=zero
     do jc=1,3 ! nu
       do ii=1,3 ! beta
         work6(1:cplex,ii,jc)=gprimd(ii,1)*work5(1:cplex,jc,1)+gprimd(ii,2)*work5(1:cplex,jc,2) &
&         +gprimd(ii,3)*work5(1:cplex,jc,3)
       end do
     end do
     do jc=1,3 ! alpha
       do ii=1,3 ! beta
         mu=gamma(jc,ii)
         work7(1:cplex,ic,mu)=gprimd(jc,1)*work6(1:cplex,ii,1)+gprimd(jc,2)*work6(1:cplex,ii,2) &
&         +gprimd(jc,3)*work6(1:cplex,ii,3)
       end do
     end do
   end do ! gamma

   do ii=1,3 ! alpha
     work3(1,ii)=gprimd(ii,1)*ddkk(2*1-1)+gprimd(ii,2)*ddkk(2*2-1) &
&     +gprimd(ii,3)*ddkk(2*3-1)
     work3(2,ii)=gprimd(ii,1)*ddkk(2*1  )+gprimd(ii,2)*ddkk(2*2  ) &
&     +gprimd(ii,3)*ddkk(2*3  )
   end do
   do ii=1,3 ! gamma
     work4(1,ii)=gmet(ii,1)*ddkk(2*1-1)+gmet(ii,2)*ddkk(2*2-1) &
&     +gmet(ii,3)*ddkk(2*3-1)
     work4(2,ii)=gmet(ii,1)*ddkk(2*1  )+gmet(ii,2)*ddkk(2*2  ) &
&     +gmet(ii,3)*ddkk(2*3  )
   end do

   do mu=1,6
     ii=alpha(mu) ! alpha
     ic=beta(mu) ! beta
     do jc=1,3 ! gamma
       work7(1:cplex,jc,mu)=work7(1:cplex,jc,mu)-half &
&       *(gprimd(ic,jc)*work3(1:cplex,ii)+gprimd(ii,jc)*work3(1:cplex,ic))
       if (ii==ic) work7(1:cplex,jc,mu)=work7(1:cplex,jc,mu)-work4(1:cplex,jc)
     end do
   end do
   do mu=1,6 ! alpha,beta
     do nu=1,3 ! gamma
       mu0=3*(mu-1)+nu
       enlout(2*mu0-1)=work7(1,nu,mu)
       enlout(2*mu0  )=work7(2,nu,mu)
     end do
   end do
   ABI_DEALLOCATE(gmet)
   ABI_DEALLOCATE(work3)
   ABI_DEALLOCATE(work4)
   ABI_DEALLOCATE(work5)
   ABI_DEALLOCATE(work6)
   ABI_DEALLOCATE(work7)
 end if

!2nd derivative wrt to 2 k wave vectors (effective mass):
! - convert from cartesian to reduced coordinates
 if ((choice==8.or.choice==81).and.signs==1.and.paw_opt<=3) then
   mu0=0 ! Shift to be applied in enlout array
   ABI_ALLOCATE(work3,(3,3))
   ABI_ALLOCATE(work4,(3,3))
   mua=1;if (choice==81) mua=2
   do ii=1,mua ! Loop Re,Im
     if (choice==8) then ! enlout is real in Voigt notation
       work3(1,1)=enlout(mu0+1) ; work3(1,2)=enlout(mu0+6) ; work3(1,3)=enlout(mu0+5)
       work3(2,1)=enlout(mu0+6) ; work3(2,2)=enlout(mu0+2) ; work3(2,3)=enlout(mu0+4)
       work3(3,1)=enlout(mu0+5) ; work3(3,2)=enlout(mu0+4) ; work3(3,3)=enlout(mu0+3)
     else                ! enlout is complex in matrix notation
       work3(1,1)=enlout(mu0+1 ) ; work3(1,2)=enlout(mu0+3 ) ; work3(1,3)=enlout(mu0+5 )
       work3(2,1)=enlout(mu0+7 ) ; work3(2,2)=enlout(mu0+9 ) ; work3(2,3)=enlout(mu0+11)
       work3(3,1)=enlout(mu0+13) ; work3(3,2)=enlout(mu0+15) ; work3(3,3)=enlout(mu0+17)
     end if
     do mu=1,3
       work4(:,mu)=gprimd(:,1)*work3(mu,1)+gprimd(:,2)*work3(mu,2)+gprimd(:,3)*work3(mu,3)
     end do
     do mu=1,3
       work3(:,mu)=gprimd(:,1)*work4(mu,1)+gprimd(:,2)*work4(mu,2)+gprimd(:,3)*work4(mu,3)
     end do
     do mu=1,3
       work4(:,mu)=gprimd(1,:)*work3(mu,1)+gprimd(2,:)*work3(mu,2)+gprimd(3,:)*work3(mu,3)
     end do
     do mu=1,3
       work3(:,mu)=gprimd(1,:)*work4(mu,1)+gprimd(2,:)*work4(mu,2)+gprimd(3,:)*work4(mu,3)
     end do
     if (choice==8) then ! enlout is real in Voigt notation
       enlout(mu0+1) = work3(1,1) ; enlout(mu0+2) = work3(2,2) ; enlout(mu0+3) = work3(3,3)
       enlout(mu0+4) = work3(3,2) ; enlout(mu0+5) = work3(1,3) ; enlout(mu0+6) = work3(2,1)
     else                ! enlout is complex in matrix notation
       enlout(mu0+1 )=work3(1,1) ; enlout(mu0+3 )=work3(1,2) ; enlout(mu0+5 )=work3(1,3)
       enlout(mu0+7 )=work3(2,1) ; enlout(mu0+9 )=work3(2,2) ; enlout(mu0+11)=work3(2,3)
       enlout(mu0+13)=work3(3,1) ; enlout(mu0+15)=work3(3,2) ; enlout(mu0+17)=work3(3,3)
     end if
     mu0=mu0+1
   end do
   ABI_DEALLOCATE(work3)
   ABI_DEALLOCATE(work4)
 end if

!2nd derivative wrt to 2 strains (elastic tensor):
! - convert from reduced to cartesian coordinates
! - substract volume contribution
 if (choice==6.and.signs==1.and.paw_opt<=3) then
   mu0=0 ! Shift to be applied in enlout array
   ABI_ALLOCATE(work1,(6))
   ABI_ALLOCATE(work2,(6))
   ABI_ALLOCATE(work3,(6+3*natom,6))
   work3(:,:)=reshape(enlout(mu0+1:mu0+6*(6+3*natom)),(/6+3*natom,6/))
   do mu=1,6
     call strconv(work3(1:6,mu),gprimd,work3(1:6,mu))
   end do
   do mu=1,6+3*natom
     work1(1:6)=work3(mu,1:6)
     call strconv(work1,gprimd,work2)
     work3(mu,1:6)=work2(1:6)
   end do
   enlout(mu0+1:mu0+6*(6+3*natom))=reshape(work3(:,:),(/6*(6+3*natom)/))
   ABI_DEALLOCATE(work1)
   ABI_DEALLOCATE(work2)
   ABI_DEALLOCATE(work3)
   call strconv(strnlk,gprimd,strnlk)
   do mub=1,6
     nub1=alpha(mub);nub2=beta(mub)
     do mua=1,6
       mu=mu0+mua+(3*natom+6)*(mub-1)
       nua1=alpha(mua);nua2=beta(mua)
       if (mua<=3.and.mub<=3) enlout(mu)=enlout(mu)+enlk
       if (mua<=3) enlout(mu)=enlout(mu)-strnlk(mub)
       if (mub<=3) enlout(mu)=enlout(mu)-strnlk(mua)
       if (nub1==nua2) enlout(mu)=enlout(mu)-0.25d0*strnlk(gamma(nua1,nub2))
       if (nub2==nua2) enlout(mu)=enlout(mu)-0.25d0*strnlk(gamma(nua1,nub1))
       if (nub1==nua1) enlout(mu)=enlout(mu)-0.25d0*strnlk(gamma(nua2,nub2))
       if (nub2==nua1) enlout(mu)=enlout(mu)-0.25d0*strnlk(gamma(nua2,nub1))
     end do
     if (mub<=3) then
       do nua1=1,natom
         nua2=3*(nua1-1);mu=mu0+nua2+6+(3*natom+6)*(mub-1)
         enlout(mu+1:mu+3)=enlout(mu+1:mu+3)-fnlk(nua2+1:nua2+3)
       end do
     end if
   end do
 end if

 if (allocated(gmet)) then
   ABI_DEALLOCATE(gmet)
 end if

!Final deallocations
!==============================================================

 if (signs==1)  then
   ABI_DEALLOCATE(fnlk)
   ABI_DEALLOCATE(ddkk)
   ABI_DEALLOCATE(strnlk)
 end if

 if (nkpgin<nkpgin_) then
   ABI_DEALLOCATE(kpgin_)
 end if
 if (nkpgout<nkpgout_) then
   ABI_DEALLOCATE(kpgout_)
 end if

 DBG_EXIT("COLL")

end subroutine nonlop_ylm
!!***

end module m_nonlop_ylm
!!***
