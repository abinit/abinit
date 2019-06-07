!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_nonlop
!! NAME
!!  m_nonlop
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (MT)
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

module m_nonlop

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_abicore
 use m_xmpi
 use m_cgtools
 use m_gemm_nonlop

 use m_time,        only : timab
 use m_hamiltonian, only : gs_hamiltonian_type, KPRIME_H_K, K_H_KPRIME, K_H_K, KPRIME_H_KPRIME
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_copy
 use m_nonlop_pl,   only : nonlop_pl
 use m_nonlop_ylm,  only : nonlop_ylm

#if defined HAVE_GPU_CUDA
 use m_manage_cuda
#endif

 implicit none

 private
!!***

 public :: nonlop
!!***

contains
!!***

!!****f* ABINIT/nonlop
!! NAME
!! nonlop
!!
!! FUNCTION
!! This routine is a driver to compute:
!! * Application of a nonlocal operator Vnl_k_k^prime in order to get:
!!    - contracted elements (energy, forces, stresses, ...), if signs=1
!!    - a function in reciprocal space (|out> = Vnl|in>),    if signs=2
!! * Optionally, in case of PAW calculation:
!!   - Application of the overlap matrix in reciprocal space
!!     (<in|S|in> or (I+S)|in>).
!!   - Application of (Vnl-lambda.S) in reciprocal space
!! According to user s choice, the routine calls a subroutine, computing all quantities:
!!   - using Legendre Polynomials Pl (Norm-conserving psps only)
!!   - using Spherical Harmonics Ylm (N-conserving or PAW ; compulsory for PAW)
!!   - using GPUs (N-conserving or PAW)
!!
!! INPUTS
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
!!                2nd derivative(s) with respect to 2 atomic pos.
!!          =33=> mixed 2nd derivative(s) with respect to strain and q vector (at q=0)
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
!!  [enl]=optional (if not present, use hamk%ekb); non-local coeffs connecting projectors
!!        see hamk%ekb description
!!  hamk <type(gs_hamiltonian_type)>=data defining the Hamiltonian at a given k (NL part involved here)
!!     | atindx1(natom)=index table for atoms, inverse of atindx
!!     | dimekb1,dimekb2=dimensions of ekb (see ham%ekb)
!!     | dimekbq=1 if enl factors do not contain a exp(-iqR) phase, 2 is they do
!!     | ekb(dimekb1,dimekb2,nspinor**2,dimekbq)=
!!     |   ->NC psps (paw_opt=0) : Kleinman-Bylander energies (hartree)
!!     |                           dimekb1=lmnmax, dimekb2=ntypat
!!     |   ->PAW (paw_opt=1 or 4): Dij coefs connecting projectors (ij symmetric)
!!     |                           dimekb1=cplex_ekb*lmnmax*(lmnmax+1)/2, dimekb2=natom
!!     |                           Complex numbers if cplex_ekb=2
!!     |                           ekb(:,:,1)= Dij^up-up, ekb(:,:,2)= Dij^dn-dn
!!     |                           ekb(:,:,3)= Dij^up-dn, ekb(:,:,4)= Dij^dn-up (only if nspinor=2)
!!     | ffnl_k(npw_k,dimffnl_k,lmnmax,ntypat)=nonlocal form factors at k
!!     | ffnl_kp(npw_kp,dimffnl_kp,lmnmax,ntypat)=nonlocal form factors at k^prime
!!     | gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!     | gprimd(3,3)=dimensional reciprocal space primitive translations
!!     | indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=ln (useylm=0) or i=lmn (useylm=1)
!!     | istwf_k=option parameter that describes the storage of wfs at k
!!     | istwf_kp=option parameter that describes the storage of wfs at k^prime
!!     | lmnmax=max. number of (l,m,n) components over all types of atoms
!!     | matblk=dimension of the arrays ph3d_k and ph3d_kp
!!     | mgfft=maximum size of 1D FFTs
!!     | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!     | mpssoang= 1+max(spin*angular momentum) for nonlocal pseudopotentials
!!     | natom=number of atoms in cell
!!     | nattyp(ntypat)=number of atoms of each type
!!     | ngfft(18)=contain all needed information about 3D FFT
!!     ! kg_k(3,npw_k)=integer coords of planewaves in basis sphere, for k
!!     ! kg_kp(3,npw_kp)=integer coords of planewaves in basis sphere, for k^prime
!!     ! kpg_k(npw_k,:)= (k+G) components and related data
!!     ! kpg_kp(npw_kp,:)=(k^prime+G) components and related data,
!!     ! kpt(3)=k point in terms of recip. translations
!!     ! kptp(3)=k^prime point in terms of recip. translations
!!     | nloalg(3)=governs the choice of the algorithm for nonlocal operator
!!     ! npw_k=number of (k+G) planewaves
!!     ! npw_kp=number of (k^prime+G) planewaves
!!     | ntypat=number of types of atoms in cell
!!     | nspinor=total number of spinorial components of the wavefunctions
!!     | ph1d(2,3*(2*mgfft+1)*natom)=1D structure factors phase information
!!     | ph3d_k(2,npw_k,matblk)=3D structure factors, for each atom and (k+g) plane wave
!!     | ph3d_kp(2,npw_kp,matblk)=3-dim structure factors, for each atom and (k^prime+g) plane wave
!!     | phkxred(2,natom)=phase factors exp(2 pi k.xred)
!!     | phkpxred(2,natom)=phase factors exp(2 pi k^prime.xred)
!!     | sij(dimekb1,ntypat)=overlap matrix components (only if paw_opt=2, 3 or 4)
!!     | ucvol=unit cell volume (bohr^3)
!!     | use_gpu_cuda=governs wheter we do the hamiltonian calculation on gpu or not
!!     | useylm=how the NL operator is to be applied: 1=using Ylm, 0=using Legendre polynomials
!!  [iatom_only]=optional. If present (and >0), only projectors related to atom of index iatom_only
!!          will be applied. (used fi to apply derivative of NL operator wrt an atomic displacement)
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2) or (choice=22,signs=2) 
!!                        - k point direction in the case (choice=5,51,or 52)
!!                          for choice 53 signs=2, cross derivatives are in idir-1 and idir+1 directions
!!                        - strain component (1:6) in the case (choice=3,signs=2) or (choice=6,signs=1)
!!                        - strain component (1:9) in the case (choice=33,signs=2)
!!                        - (1:9) components to specify the atom to be moved and the second q-gradient 
!!                          direction in the case (choice=25,signs=2)
!!  lambda=factor to be used when computing (Vln-lambda.S) - only for paw_opt=2
!!         Typically lambda is the eigenvalue (or its guess)
!!  mpi_enreg=informations about MPI parallelization
!!  ndat=number of wavefunctions on which to apply nonlop
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
!!  [only_SO]=optional, flag to calculate only the SO part in nonlop
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!  [qdir]= optional, direction of the q-gradient (only for choice=22, choice=25 and choice=33)
!!  [select_k]=optional, option governing the choice of k points to be used.
!!             hamk datastructure contains quantities needed to apply NL operator
!!             in reciprocal space between 2 kpoints, k and k^prime (equal in most cases);
!!             if select_k=1, <k^prime|Vnl|k>       is applied [default]
!!             if select_k=2, <k|Vnl|k^prime>       is applied
!!             if select_k=3, <k|Vnl|k>             is applied
!!             if select_k=4, <k^prime|Vnl|k^prime> is applied
!!  signs= if 1, get contracted elements (energy, forces, stress, ...)
!!         if 2, applies the non-local operator to a function in reciprocal space
!!  tim_nonlop=timing code of the calling routine (can be set to 0 if not attributed)
!!  vectin(2,npwin*my_nspinor*ndat)=input cmplx wavefunction coefficients <G|Cnk>
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
!!  ==== ONLY IF useylm=1
!!  cprjin(natom,my_nspinor*ndat) <type(pawcprj_type)>=projected input wave function |in> on non-local projectors
!!                                  =<p_lmn|in> and derivatives
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
!!                                  other derivatives are not computed
!!                                  This option is not compatible with choice=4,24 or 6
!!                     If useylm=0, must have cpopt=-1!
!!                     Warning: for cpopt= 1 or 3, derivatives wrt strains do not contain
!!                              the contribution due to the volume change;
!!                              i.e. <dp_lmn/dEps|in> are incomplete.
!!
!! NOTES
!! * See nonlop_pl and nonlop_ylm to have more comments...
!! * In the case signs=1, the array vectout is not used.
!!
!! PARENTS
!!      d2frnl,dfpt_nsteltwf,dfptnl_resp,energy,fock_getghc,forstrnps,getgh1c
!!      getgh1dq,getgh2c,getghc,getgsc,m_invovl,m_lobpcgwf,make_grad_berry,
!!      nonlop_test,prep_nonlop,vtowfk,wfd_vnlpsi
!!
!! CHILDREN
!!      gemm_nonlop,nonlop_gpu,nonlop_pl,nonlop_ylm,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,timab
!!
!! SOURCE

subroutine nonlop(choice,cpopt,cprjin,enlout,hamk,idir,lambda,mpi_enreg,ndat,nnlout,&
&                 paw_opt,signs,svectout,tim_nonlop,vectin,vectout,&
&                 enl,iatom_only,only_SO,qdir,select_k) !optional arguments

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cpopt,idir,ndat,nnlout,paw_opt,signs,tim_nonlop
 integer,intent(in),optional :: iatom_only,only_SO,qdir,select_k
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(in),target :: hamk
!arrays
 real(dp),intent(in) :: lambda(ndat)
 real(dp),intent(in),target,optional :: enl(:,:,:,:)
 real(dp),intent(inout),target :: vectin(:,:)
 real(dp),intent(out),target :: enlout(:),svectout(:,:)
 real(dp),intent(inout),target :: vectout(:,:)
 type(pawcprj_type),intent(inout),target :: cprjin(:,:)

!Local variables-------------------------------
!scalars
 integer :: dimenl1,dimenl2,dimenl2_,dimekbq,dimffnlin,dimffnlout,dimsij,iatm,iatom_only_,idat
 integer :: ii,ispden,ispinor,istwf_k,itypat,jspinor,matblk_,my_nspinor,n1,n2,n3,natom_,ncpgr_atm
 integer :: nkpgin,nkpgout,npwin,npwout,ntypat_,only_SO_,select_k_,shift1,shift2,shift3
 logical :: atom_pert,force_recompute_ph3d,kpgin_allocated,kpgout_allocated,use_gemm_nonlop
 character(len=500) :: msg
!arrays
 integer :: nlmn_atm(1),nloalg_(3)
 integer,pointer :: kgin(:,:),kgout(:,:)
 integer, ABI_CONTIGUOUS pointer :: atindx1_(:),indlmn_(:,:,:),nattyp_(:)
 real(dp) :: tsec(2)
 real(dp),pointer :: enl_ptr(:,:,:,:)
 real(dp),pointer :: ffnlin(:,:,:,:),ffnlin_(:,:,:,:),ffnlout(:,:,:,:),ffnlout_(:,:,:,:)
 real(dp),pointer :: kpgin(:,:),kpgout(:,:)
 real(dp) :: kptin(3),kptout(3)
 real(dp),pointer :: ph3din(:,:,:),ph3din_(:,:,:),ph3dout(:,:,:),ph3dout_(:,:,:)
 real(dp),pointer :: phkxredin(:,:),phkxredin_(:,:),phkxredout(:,:),phkxredout_(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: ph1d_(:,:),sij_(:,:)
 real(dp), pointer :: enl_(:,:,:,:)
 type(pawcprj_type),pointer :: cprjin_(:,:)
  integer :: b0,b1,b2,b3,b4,e0,e1,e2,e3,e4

! **********************************************************************

 DBG_ENTER("COLL")

!Keep track of time spent in this routine (selection of different slots for different choices)
 call timab(220+tim_nonlop,1,tsec)

 only_SO_=0; if (present(only_SO)) only_SO_=only_SO
 my_nspinor=max(1,hamk%nspinor/mpi_enreg%nproc_spinor)

 force_recompute_ph3d=.false.

!Error(s) on incorrect input
 if (hamk%useylm==0) then
   if (paw_opt>0) then
     MSG_BUG('When paw_opt>0 you must use ylm version of nonlop! Set useylm 1.')
   end if
   if (cpopt/=-1) then
     MSG_BUG('If useylm=0, ie no PAW, then cpopt/=-1 is not allowed !')
   end if
   if (hamk%dimekbq/=1) then
     MSG_BUG('If useylm=0, ie no PAW, then dimekbq/=-1 is not allowed !')
   end if
   if (hamk%use_gpu_cuda/=0) then
     msg = 'When use_gpu_cuda/=0 you must use ylm version of nonlop! Set useylm 1.'
     MSG_BUG(msg)
   end if
 end if
 if (hamk%use_gpu_cuda/=0.and.hamk%dimekbq/=1) then
   msg = 'GPU version of nonlop not compatible with a exp(-iqR) phase!'
   MSG_BUG(msg)
 end if
 if ((.not.associated(hamk%kg_k)).or.(.not.associated(hamk%kg_kp))) then
   MSG_BUG('kg_k/kg_kp should be associated!')
 end if
 if ((.not.associated(hamk%ffnl_k)).or.(.not.associated(hamk%ffnl_kp))) then
   MSG_BUG('ffnl_k/ffnl_kp should be associated!')
 end if
!if (hamk%istwf_k/=hamk%istwf_kp) then
!  msg = 'istwf has to be the same for both k-points.'
!  MSG_BUG(msg)
!end if

!Select k-dependent objects according to select_k input parameter
 select_k_=1;if (present(select_k)) select_k_=select_k
 nkpgin=0;nkpgout=0;nullify(kpgin);nullify(kpgout)
 nullify(ph3din);nullify(ph3dout)
 if (select_k_==KPRIME_H_K) then
!  ===== <k^prime|Vnl|k> =====
   kptin = hamk%kpt_k ; kptout = hamk%kpt_kp
   npwin=hamk%npw_fft_k ; npwout=hamk%npw_fft_kp
   kgin => hamk%kg_k ; kgout => hamk%kg_kp
   if (associated(hamk%kpg_k)) then
     kpgin => hamk%kpg_k ; nkpgin=size(kpgin,2)
   end if
   if (associated(hamk%kpg_kp)) then
     kpgout => hamk%kpg_kp ; nkpgout=size(kpgout,2)
   end if
   phkxredin => hamk%phkxred ; phkxredout => hamk%phkpxred
   ffnlin => hamk%ffnl_k ; ffnlout => hamk%ffnl_kp
   if (associated(hamk%ph3d_k )) ph3din  => hamk%ph3d_k
   if (associated(hamk%ph3d_kp)) ph3dout => hamk%ph3d_kp
   force_recompute_ph3d=(.not.(associated(hamk%ph3d_k).and.associated(hamk%ph3d_kp)))
   istwf_k=hamk%istwf_k
 else if (select_k_==K_H_KPRIME) then
!  ===== <k|Vnl|k^prime> =====
   kptin = hamk%kpt_kp ; kptout = hamk%kpt_k
   npwin=hamk%npw_fft_kp ; npwout=hamk%npw_fft_k
   kgin => hamk%kg_kp ; kgout => hamk%kg_k
   if (associated(hamk%kpg_kp)) then
     kpgin => hamk%kpg_kp ; nkpgin=size(kpgin,2)
   end if
   if (associated(hamk%kpg_k)) then
     kpgout => hamk%kpg_k ; nkpgout=size(kpgout,2)
   end if
   phkxredin => hamk%phkpxred ; phkxredout => hamk%phkxred
   ffnlin => hamk%ffnl_kp ; ffnlout => hamk%ffnl_k
   if (associated(hamk%ph3d_kp)) ph3din  => hamk%ph3d_kp
   if (associated(hamk%ph3d_k )) ph3dout => hamk%ph3d_k
   force_recompute_ph3d=(.not.(associated(hamk%ph3d_kp).and.associated(hamk%ph3d_k)))
   istwf_k=hamk%istwf_kp
 else if (select_k_==K_H_K) then
!  ===== <k|Vnl|k> =====
   kptin = hamk%kpt_k ; kptout = hamk%kpt_k
   npwin=hamk%npw_fft_k ; npwout=hamk%npw_fft_k
   kgin => hamk%kg_k ; kgout => hamk%kg_k
   if (associated(hamk%kpg_k)) then
     kpgin => hamk%kpg_k ; nkpgin=size(kpgin,2)
   end if
   if (associated(hamk%kpg_k)) then
     kpgout => hamk%kpg_k ; nkpgout=size(kpgout,2)
   end if
   phkxredin => hamk%phkxred ; phkxredout => hamk%phkxred
   ffnlin => hamk%ffnl_k ; ffnlout => hamk%ffnl_k
   if (associated(hamk%ph3d_k)) ph3din  => hamk%ph3d_k
   if (associated(hamk%ph3d_k)) ph3dout => hamk%ph3d_k
   force_recompute_ph3d=(.not.(associated(hamk%ph3d_k)))
   istwf_k=hamk%istwf_k
 else if (select_k_==KPRIME_H_KPRIME) then
!  ===== <k^prime|Vnl|k^prime> =====
   kptin = hamk%kpt_kp ; kptout = hamk%kpt_kp
   npwin=hamk%npw_fft_kp ; npwout=hamk%npw_fft_kp
   kgin => hamk%kg_kp ; kgout => hamk%kg_kp
   if (associated(hamk%kpg_kp)) then
     kpgin => hamk%kpg_kp ; nkpgin=size(kpgin,2)
   end if
   if (associated(hamk%kpg_kp)) then
     kpgout => hamk%kpg_kp ; nkpgout=size(kpgout,2)
   end if
   phkxredin => hamk%phkpxred ; phkxredout => hamk%phkpxred
   ffnlin => hamk%ffnl_kp ; ffnlout => hamk%ffnl_kp
   if (associated(hamk%ph3d_kp)) ph3din  => hamk%ph3d_kp
   if (associated(hamk%ph3d_kp)) ph3dout => hamk%ph3d_kp
   force_recompute_ph3d=(.not.(associated(hamk%ph3d_kp)))
   istwf_k=hamk%istwf_kp
 end if

 if (npwin==0.or.npwout==0) return
 dimffnlin=size(ffnlin,2);dimffnlout=size(ffnlout,2)
 kpgin_allocated=(.not.associated(kpgin))
 if (kpgin_allocated) then
   ABI_ALLOCATE(kpgin,(npwin,0))
 end if
 kpgout_allocated=(.not.associated(kpgout))
 if (kpgout_allocated) then
   ABI_ALLOCATE(kpgout,(npwout,0))
 end if

!Check some sizes for safety
!if (paw_opt==0.or.cpopt<2.or.((cpopt==2.or.cpopt==3).and.choice>1)) then
 if (size(ffnlin,1)/=npwin.or.size(ffnlin,3)/=hamk%lmnmax) then
   msg = 'Incorrect size for ffnlin!'
!   MSG_BUG(msg)
 end if
 if(signs==2) then
   if (size(ffnlout,1)/=npwout.or.size(ffnlout,3)/=hamk%lmnmax) then
     msg = 'Incorrect size for ffnlout!'
     MSG_BUG(msg)
   end if
 end if
!This test is OK only because explicit sizes are passed to nonlop_* routines
 if (size(vectin)<2*npwin*my_nspinor*ndat) then
   msg = 'Incorrect size for vectin!'
   MSG_BUG(msg)
 end if
 if(choice/=0.and.signs==2) then
   if(paw_opt/=3) then
!    This test is OK only because explicit sizes are passed to nonlop_* routines
     if (size(vectout)<2*npwout*my_nspinor*ndat) then
       msg = 'Incorrect size for vectout!'
       MSG_BUG(msg)
     end if
   end if
   if(paw_opt>=3) then
     if (size(svectout)<2*npwout*my_nspinor*ndat) then
       msg = 'Incorrect size for svectout!'
       MSG_BUG(msg)
     end if
   end if
 end if
 if(cpopt>=0) then
   if (size(cprjin)/=hamk%natom*my_nspinor*ndat) then
     msg = 'Incorrect size for cprjin!'
     MSG_BUG(msg)
   end if
 end if

!Non-local coefficients connecting projectors:
!If enl is present in the arg list, use it; instead use hamk%ebk
 if (present(enl)) then
   enl_ptr => enl
   dimenl1=size(enl,1);dimenl2=size(enl,2);dimekbq=size(enl,4)
 else
   enl_ptr => hamk%ekb
   dimenl1=hamk%dimekb1;dimenl2=hamk%dimekb2;dimekbq=1
 end if

!In the case of a derivative with respect to an atomic displacement,
!and if <g|dVnl/dR|c> is required (signs=2), we only need to compute the
!derivatives of the projectors associated with the displaced atom.
 iatom_only_=-1;if (present(iatom_only)) iatom_only_=iatom_only
 atom_pert=((signs==2).and.(choice==2.or.choice==4.or.choice==22.or.choice==24.or.choice==25.or.choice==54))

 if (iatom_only_>0.and.atom_pert) then
!   We consider only atom with index iatom_only
   iatm=hamk%atindx(iatom_only_);itypat=hamk%typat(iatom_only_)
   natom_=1 ; ntypat_=1 ; dimenl2_=1 ; matblk_=1
   nloalg_(:)=hamk%nloalg(:)
   ABI_ALLOCATE(atindx1_,(1))
   ABI_ALLOCATE(nattyp_,(1))
   atindx1_(1)=1 ; nattyp_(1)=1
!  Store at the right place the 1d phases
   n1=hamk%ngfft(1);n2=hamk%ngfft(2);n3=hamk%ngfft(3)
   ABI_ALLOCATE(ph1d_,(2,(2*n1+1)+(2*n2+1)+(2*n3+1)))
   shift1=(iatm-1)*(2*n1+1)
   ph1d_(:,1:2*n1+1)=hamk%ph1d(:,1+shift1:2*n1+1+shift1)
   shift2=(iatm-1)*(2*n2+1)+hamk%natom*(2*n1+1)
   ph1d_(:,1+2*n1+1:2*n2+1+2*n1+1)=hamk%ph1d(:,1+shift2:2*n2+1+shift2)
   shift3=(iatm-1)*(2*n3+1)+hamk%natom*(2*n1+1+2*n2+1)
   ph1d_(:,1+2*n1+1+2*n2+1:2*n3+1+2*n2+1+2*n1+1)=hamk%ph1d(:,1+shift3:2*n3+1+shift3)
   ABI_ALLOCATE(phkxredin_,(2,1))
   ABI_ALLOCATE(phkxredout_,(2,1))
   phkxredin_(:,1)=phkxredin(:,iatm)
   phkxredout_(:,1)=phkxredout(:,iatm)
   ABI_ALLOCATE(ph3din_,(2,npwin,1))
   ABI_ALLOCATE(ph3dout_,(2,npwout,1))
   if (force_recompute_ph3d.or.hamk%matblk<hamk%natom) then
     nloalg_(2)=-abs(nloalg_(2)) !Will compute the 3D phase factors inside nonlop
   else
     ph3din_(:,1:npwin,1)=ph3din(:,1:npwin,iatm)
     ph3dout_(:,1:npwout,1)=ph3dout(:,1:npwout,iatm)
   end if
   ABI_ALLOCATE(ffnlin_,(npwin,dimffnlin,hamk%lmnmax,1))
   ABI_ALLOCATE(ffnlout_,(npwout,dimffnlout,hamk%lmnmax,1))
   ffnlin_(:,:,:,1)=ffnlin(:,:,:,itypat)
   ffnlout_(:,:,:,1)=ffnlout(:,:,:,itypat)
   ABI_DATATYPE_ALLOCATE(cprjin_,(1,my_nspinor*((cpopt+5)/5)))
   if (cpopt>=0) then
     nlmn_atm(1)=cprjin(iatm,1)%nlmn
     ncpgr_atm=cprjin(iatm,1)%ncpgr
     call pawcprj_alloc(cprjin_,ncpgr_atm,nlmn_atm)
     do idat=1,ndat
       do ispinor=1,my_nspinor
         jspinor=ispinor+(idat-1)*my_nspinor
         call pawcprj_copy(cprjin(iatm:iatm,jspinor:jspinor),cprjin_(1:1,jspinor:jspinor))
       end do
     end do
   end if
   if (size(enl_ptr)>0) then
     ABI_ALLOCATE(enl_,(size(enl_ptr,1),1,hamk%nspinor**2,size(enl_ptr,4)))
     do ii=1,size(enl_ptr,4)
       do ispden=1,hamk%nspinor**2
         if (dimenl2==hamk%natom .and. hamk%usepaw==1) then
           enl_(:,1,ispden,ii)=enl_ptr(:,iatom_only_,ispden,ii)
         else if (dimenl2==hamk%ntypat) then
           enl_(:,1,ispden,ii)=enl_ptr(:,itypat,ispden,ii)
         else
           enl_(:,1,ispden,ii)=enl_ptr(:,1,ispden,ii)
         end if
       end do
     end do
   else
     ABI_ALLOCATE(enl_,(0,0,0,0))
   end if
   if (allocated(hamk%sij)) then
     dimsij=size(hamk%sij,1)
     ABI_ALLOCATE(sij_,(dimsij,1))
     if (size(hamk%sij,2)==hamk%ntypat) then
       sij_(:,1)=hamk%sij(:,itypat)
     else if (size(hamk%sij)>0) then
       sij_(:,1)=hamk%sij(:,1)
     end if
   end if
   ABI_ALLOCATE(indlmn_,(6,hamk%lmnmax,1))
   indlmn_(:,:,1)=hamk%indlmn(:,:,itypat)

 else
!  Usual case: all atoms are processed
   natom_  =hamk%natom; ntypat_=hamk%ntypat
   dimenl2_=dimenl2   ; matblk_=hamk%matblk
   nloalg_(:)  = hamk%nloalg(:)
   atindx1_    => hamk%atindx1
   nattyp_     => hamk%nattyp
   ph1d_       => hamk%ph1d
   phkxredin_  => phkxredin
   phkxredout_ => phkxredout
   ffnlin_     => ffnlin
   ffnlout_    => ffnlout
   cprjin_     => cprjin
   enl_        => enl_ptr
   sij_        => hamk%sij
   indlmn_     => hamk%indlmn
   if (force_recompute_ph3d) then
     nloalg_(2)=-abs(nloalg_(2)) !Will compute the 3D phase factors inside nonlop
     ABI_ALLOCATE(ph3din_,(2,npwin,hamk%matblk))
     ABI_ALLOCATE(ph3dout_,(2,npwout,hamk%matblk))
   else
     ph3din_     => ph3din
     ph3dout_    => ph3dout
   end if
 end if

!A specific version of nonlop based on BLAS3 can be used
!But there are several restrictions

! use_gemm_nonlop= ( gemm_nonlop_use_gemm .and. &
!& signs == 2 .and. paw_opt /= 2 .and. hamk%nspinor == 1 .and. &
!& cpopt < 3 .and. hamk%useylm /= 0 .and. &
!& (choice < 2 .or. choice == 7) )

 use_gemm_nonlop= ( gemm_nonlop_use_gemm .and. &
& signs == 2 .and. paw_opt /= 2 .and. &
& cpopt < 3 .and. hamk%useylm /= 0 .and. &
& (choice < 2 .or. choice == 7) )

 if(use_gemm_nonlop) then
   call gemm_nonlop(atindx1_,choice,cpopt,cprjin_,dimenl1,dimenl2_,dimekbq,&
&   dimffnlin,dimffnlout,enl_,enlout,ffnlin_,ffnlout_,hamk%gmet,hamk%gprimd,&
&   idir,indlmn_,istwf_k,kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,&
&   hamk%lmnmax,matblk_,hamk%mgfft,mpi_enreg,hamk%mpsang,hamk%mpssoang,&
&   natom_,nattyp_,ndat,hamk%ngfft,nkpgin,nkpgout,nloalg_,&
&   nnlout,npwin,npwout,my_nspinor,hamk%nspinor,ntypat_,only_SO_,paw_opt,&
&   phkxredin_,phkxredout_,ph1d_,ph3din_,ph3dout_,signs,sij_,svectout,&
&   tim_nonlop,hamk%ucvol,hamk%useylm,vectin,vectout,hamk%use_gpu_cuda)

 else
   !$omp parallel do default(shared), &
   !$omp& firstprivate(ndat,npwin,my_nspinor,choice,signs,paw_opt,npwout,cpopt,nnlout), &
   !$omp& private(b0,b1,b2,b3,b4,e0,e1,e2,e3,e4)
   !!$omp& schedule(static), if(hamk%use_gpu_cuda==0)
   do idat=1, ndat
     !vectin_idat => vectin(:,1+npwin*my_nspinor*(idat-1):npwin*my_nspinor*idat)
     b0 = 1+npwin*my_nspinor*(idat-1)
     e0 = npwin*my_nspinor*idat
     if (choice/=0.and.signs==2.and.paw_opt/=3) then
       !vectout_idat => vectout(:,1+npwout*my_nspinor*(idat-1):npwout*my_nspinor*idat)
       b1 = 1+npwout*my_nspinor*(idat-1)
       e1 = npwout*my_nspinor*idat
     else
       !vectout_idat => vectout
       b1 = lbound(vectout,dim=2)
       e1 = ubound(vectout,dim=2)
     end if
     if (choice/=0.and.signs==2.and.paw_opt>=3) then
       !svectout_idat => svectout(:,1+npwout*my_nspinor*(idat-1):npwout*my_nspinor*idat)
       b2 = 1+npwout*my_nspinor*(idat-1)
       e2 = npwout*my_nspinor*idat
     else
       !svectout_idat => svectout
       b2 = lbound(svectout,dim=2)
       e2 = ubound(svectout,dim=2)
     end if

     if (cpopt>=0) then
       !cprjin_idat => cprjin_(:,my_nspinor*(idat-1)+1:my_nspinor*(idat))
       b3 = my_nspinor*(idat-1)+1
       e3 = my_nspinor*(idat)
     else
       !cprjin_idat => cprjin_
       b3 = lbound(cprjin_,dim=2)
       e3 = ubound(cprjin_,dim=2)
     end if
     if (nnlout>0) then
      !enlout_idat => enlout((idat-1)*nnlout+1:(idat*nnlout))
       b4 = (idat-1)*nnlout+1
       e4 = (idat*nnlout)
     else
      !enlout_idat => enlout
       b4 = lbound(enlout,dim=1)
       e4 = ubound(enlout,dim=1)
     end if

!    Legendre Polynomials version
     if (hamk%useylm==0) then
       call nonlop_pl(choice,dimenl1,dimenl2_,dimffnlin,dimffnlout,enl_,&
&       enlout(b4:e4),ffnlin_,ffnlout_,hamk%gmet,hamk%gprimd,idir,indlmn_,istwf_k,&
&       kgin,kgout,kpgin,kpgout,kptin,kptout,hamk%lmnmax,matblk_,hamk%mgfft,&
&       mpi_enreg,hamk%mpsang,hamk%mpssoang,natom_,nattyp_,hamk%ngfft,&
&       nkpgin,nkpgout,nloalg_,npwin,npwout,my_nspinor,hamk%nspinor,&
&       ntypat_,only_SO_,phkxredin_,phkxredout_,ph1d_,ph3din_,ph3dout_,signs,hamk%ucvol,&
&       vectin(:,b0:e0),vectout(:,b1:e1))
!    Spherical Harmonics version
     else if (hamk%use_gpu_cuda==0) then
       call nonlop_ylm(atindx1_,choice,cpopt,cprjin_(:,b3:e3),dimenl1,dimenl2_,dimekbq,&
&       dimffnlin,dimffnlout,enl_,enlout(b4:e4),ffnlin_,ffnlout_,hamk%gprimd,idir,&
&       indlmn_,istwf_k,kgin,kgout,kpgin,kpgout,kptin,kptout,lambda(idat),&
&       hamk%lmnmax,matblk_,hamk%mgfft,mpi_enreg,natom_,nattyp_,hamk%ngfft,&
&       nkpgin,nkpgout,nloalg_,nnlout,npwin,npwout,my_nspinor,hamk%nspinor,&
&       ntypat_,paw_opt,phkxredin_,phkxredout_,ph1d_,ph3din_,ph3dout_,signs,sij_,&
&       svectout(:,b2:e2),hamk%ucvol,vectin(:,b0:e0),vectout(:,b1:e1),qdir=qdir)
!    GPU version
     else
       call nonlop_gpu(atindx1_,choice,cpopt,cprjin(:,b3:e3),dimenl1,dimenl2_,&
&       dimffnlin,dimffnlout,enl_,enlout(b4:e4),ffnlin_,ffnlout_,hamk%gprimd,idir,&
&       indlmn_,istwf_k,kgin,kgout,kpgin,kpgout,kptin,kptout,lambda(idat),&
&       hamk%lmnmax,matblk_,hamk%mgfft,mpi_enreg,natom_,nattyp_,hamk%ngfft,&
&       nkpgin,nkpgout,nloalg_,nnlout,npwin,npwout,my_nspinor,hamk%nspinor,&
&       ntypat_,paw_opt,phkxredin_,phkxredout_,ph1d_,ph3din_,ph3dout_,signs,sij_,&
&       svectout(:,b2:e2),hamk%ucvol,vectin(:,b0:e0),vectout(:,b1:e1))
     end if
   end do
   !$omp end parallel do
 end if

!Release temporary storage
 if (iatom_only_>0.and.atom_pert) then
   if (cpopt>=0) then
     call pawcprj_free(cprjin_)
   end if
   ABI_DEALLOCATE(atindx1_)
   ABI_DEALLOCATE(nattyp_)
   ABI_DEALLOCATE(ph1d_)
   ABI_DEALLOCATE(ph3din_)
   ABI_DEALLOCATE(ph3dout_)
   ABI_DEALLOCATE(phkxredin_)
   ABI_DEALLOCATE(phkxredout_)
   ABI_DEALLOCATE(ffnlin_)
   ABI_DEALLOCATE(ffnlout_)
   ABI_DEALLOCATE(enl_)
   ABI_DEALLOCATE(indlmn_)
   ABI_DATATYPE_DEALLOCATE(cprjin_)
   if (allocated(hamk%sij)) then
     ABI_DEALLOCATE(sij_)
   end if
 else if (force_recompute_ph3d) then
   ABI_DEALLOCATE(ph3din_)
   ABI_DEALLOCATE(ph3dout_)
 end if
 if (kpgin_allocated) then
   ABI_DEALLOCATE(kpgin)
 end if
 if (kpgout_allocated) then
   ABI_DEALLOCATE(kpgout)
 end if

 call timab(220+tim_nonlop,2,tsec)

 DBG_EXIT("COLL")

end subroutine nonlop
!!***

!!****f* ABINIT/nonlop_gpu
!! NAME
!! nonlop_gpu
!!
!! FUNCTION
!!  Compute application of a nonlocal operator, using GPU (NVidia Cuda)
!!  This routine is an interface to Cuda Kernel gpu_nonlop.cu
!!
!! COPYRIGHT
!! Copyright (C) 2011-2019 ABINIT group (FDahm, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  choice: chooses possible output:
!!    choice=0 => do nothing (only compute WF projected with NL projectors)
!!          =1 => a non-local energy contribution
!!          =2 => a gradient with respect to atomic position(s)
!!          =3 => a gradient with respect to strain(s)
!!          =23=> a gradient with respect to atm. pos. and strain(s)
!!  cpopt=flag defining the status of cprjin%cp(:)=<Proj_i|Cnk> scalars (see below, side effects)
!!  dimenl1,dimenl2=dimensions of enl (see enl)
!!  dimffnlin=second dimension of ffnlin (1+number of derivatives)
!!  dimffnlout=second dimension of ffnlout (1+number of derivatives)
!!  enl(dimenl1,dimenl2,nspinortot**2)=
!!  ->Norm conserving : ==== when paw_opt=0 ====
!!                      (Real) Kleinman-Bylander energies (hartree)
!!                      dimenl1=lmnmax  -  dimenl2=ntypat
!!  ->PAW :             ==== when paw_opt=1, 2 or 4 ====
!!                      (Real or complex, hermitian) Dij coefs to connect projectors
!!                      dimenl1=cplex_enl*lmnmax*(lmnmax+1)/2  -  dimenl2=natom
!!  ffnlin(npwin,dimffnlin,lmnmax,ntypat)=nonlocal form factors to be used
!!          for the application of the nonlocal operator to the |in> vector
!!  ffnlout(npwout,dimffnlout,lmnmax,ntypat)=nonlocal form factors to be used
!!          for the application of the nonlocal operator to the |out> vector
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2),
!!                        - k point direction in the case (choice=5,signs=2)
!!                          for choice 53, twisted derivative involves idir+1 and idir-1
!!                        - strain component (1:6) in the case (choice=3,signs=2) or (choice=6,signs=1)
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
!!         choice=1=>nnlout=1   choice=2=>nnlout=3*natom    choice=3=>nnlout=6
!!         ==== if paw_opt=3 ====
!!         choice=1 =>nnlout=1
!!         ==== if paw_opt=4 ====
!!         not available
!!  npwin=number of planewaves for given k point, for the |in> vector
!!  npwout=number of planewaves for given k point, for the |out> vector
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nspinortot=number of spinorial components of the wavefunctions on current proc
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
!!  signs= if 1, get contracted elements (energy, forces, stress, ...)
!!         if 2, applies the non-local operator to a function in reciprocal space
!!  sij(dimenl1,ntypat*(paw_opt/3))=overlap matrix components (only if paw_opt=2, 3 or 4)
!!  ucvol=unit cell volume (bohr^3)
!!  vectin(2,npwin*nspinor)=input cmplx wavefunction coefficients <G|Cnk>
!!  [cprjin_left(natom,nspinor)]=The projected input wave function <p_nlm|in_left>
!!    for the left wavefunction. Data are assumed to be in memory, they are NOT recalculated here.
!!    Only signs==1 and choice==1 are supported.
!!
!! OUTPUT
!! ==== if (signs==1) ====
!! --If (paw_opt==0, 1 or 2)
!!    enlout(nnlout)= contribution to the non-local part of the following properties:
!!      if choice=1 : enlout(1)               -> the energy
!!      if choice=2 : enlout(1:3*natom)       -> the forces
!!      if choice=3 : enlout(1:6)             -> the stresses
!!      if choice=23: enlout(1:6+3*natom)     -> the forces and the stresses
!! --If (paw_opt==3)
!!    if choice=1 : enlout(nnlout)= contribution to <c|S|c>  (nnlout=1)
!! --If (paw_opt==4)
!!    not available
!! ==== if (signs==2) ====
!! --if (paw_opt=0, 1 or 4)
!!    vectout(2,npwout*nspinor)=result of the aplication of the concerned operator
!!                or one of its derivatives to the input vect.:
!!      if (choice=1) <G|V_nonlocal|vect_start>
!!      if (choice=2) <G|dV_nonlocal/d(atm coord)|vect_start>
!!      if (choice=3) <G|dV_nonlocal/d(strain)|vect_start>
!!  if (paw_opt=2)
!!    vectout(2,npwout*nspinor)=final vector in reciprocal space:
!!      if (choice=1) <G|V_nonlocal-lamdba.(I+S)|vect_start>
!!      if (choice=2) <G|d[V_nonlocal-lamdba.(I+S)]/d(atm coord)|vect_start>
!!      if (choice=3) <G|d[V_nonlocal-lamdba.(I+S)]/d(strain)|vect_start>
!! --if (paw_opt=3 or 4)
!!    svectout(2,npwout*nspinor)=result of the aplication of Sij (overlap matrix)
!!                  or one of its derivatives to the input vect.:
!!      if (choice=1) <G|I+S|vect_start>
!!      if (choice=2) <G|dS/d(atm coord)|vect_start>
!!      if (choice=3) <G|dS/d(strain)|vect_start>
!!
!! SIDE EFFECTS
!!  cprjin(natom,nspinor) <type(pawcprj_type)>=projected input wave function |in> on non-local projectors
!!                                  =<p_lmn|in> and derivatives
!!                    Treatment depends on cpopt parameter:
!!                     if cpopt=-1, <p_lmn|in> (and derivatives)
!!                                  are computed here (and not saved)
!!                     if cpopt= 0, <p_lmn|in> are computed here and saved
!!                                  derivatives are eventually computed but not saved
!!                     if cpopt= 1, <p_lmn|in> and first derivatives are computed here and saved
!!                                  other derivatives are eventually computed but not saved
!!
!! TODO
!! * Implementation for spinorial wave functions (nspinor=2)
!! * Implementation for response function (phonons, ddk, elastic tensor, ...)
!!
!! PARENTS
!!      nonlop
!!
!! CHILDREN
!!      dotprod_g,gpu_nonlop
!!
!! SOURCE

 subroutine nonlop_gpu(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimffnlin,dimffnlout,&
&                      enl,enlout,ffnlin,ffnlout,gprimd,idir,indlmn,istwf_k,&
&                      kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                      mpi_enreg,natom,nattyp,ngfft,nkpgin,nkpgout,nloalg,nnlout,&
&                      npwin,npwout,nspinor,nspinortot,ntypat,paw_opt,phkxredin,phkxredout,ph1d,&
&                      ph3din,ph3dout,signs,sij,svectout,ucvol,vectin,vectout)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimffnlin,dimffnlout,idir
 integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,natom,nkpgin,nkpgout,nnlout
 integer,intent(in) :: npwin,npwout,nspinor,nspinortot,ntypat,paw_opt,signs
 real(dp),intent(in) :: lambda,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),indlmn(6,lmnmax,ntypat),kgin(3,npwin)
 integer,intent(in) :: kgout(3,npwout),nattyp(ntypat),ngfft(18),nloalg(3)
 real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2)
 real(dp),intent(in) :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
 real(dp),intent(in) :: ffnlout(npwout,dimffnlout,lmnmax,ntypat) !,gmet(3,3)
 real(dp),intent(in) :: gprimd(3,3),kpgin(npwin,nkpgin),kpgout(npwout,nkpgout)
 real(dp),intent(in) :: kptin(3),kptout(3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: phkxredin(2,natom),phkxredout(2,natom)
 real(dp),intent(in) :: sij(dimenl1,ntypat*((paw_opt+1)/3))
 real(dp),intent(inout) :: ph3din(2,npwin,matblk),ph3dout(2,npwout,matblk)
 real(dp),intent(inout) :: vectin(:,:)
 real(dp),intent(out) :: enlout(:)
 real(dp),intent(out),target :: svectout(:,:)
 real(dp),intent(out),target :: vectout (:,:)
 type(pawcprj_type),intent(inout) :: cprjin(:,:)

!Local variables-------------------------------
!scalars
 integer :: ia,iatom,ilmn,iproj,ispinor,itypat,signs_
 real(dp) :: doti
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: proj(:,:)
 real(dp),pointer :: svectout_(:,:),vectout_(:,:)

! **********************************************************************

 DBG_ENTER("COLL")

!Error on bad choice
 if ((choice<0 .or. choice>3).and. choice/=23 .and. choice/=24) then
   write(msg,'(a,i0,a)')'Does not presently support this choice=',choice,'.'
   MSG_BUG(msg)
 end if
 if (cpopt<-1.or.cpopt>1) then
   msg='  Bad value for cpopt !'
   MSG_BUG(msg)
 end if
 if (nspinor==2) then
   msg='  nspinor=2 (spinorial WF) not yet allowed !'
   MSG_ERROR(msg)
 end if

 if ((cpopt==0).or.(cpopt==1))  then
   ABI_ALLOCATE(proj,(2,lmnmax*natom))
   proj=zero;
 end if

!Workaround to get choice=1/signs=1 working
 if (choice==1.and.signs==1) then
   signs_=2
   ABI_ALLOCATE(vectout_,(2,npwin*nspinor))
   ABI_ALLOCATE(svectout_,(2,npwin*nspinor*(paw_opt/3)))
 else
   signs_=signs;vectout_=>vectout;svectout_=>svectout
 end if

#if defined HAVE_GPU_CUDA
 call gpu_nonlop(atindx1,choice,cpopt,proj,dimenl1,dimenl2,dimffnlin,dimffnlout,&
& enl,enlout,ffnlin,ffnlout,gprimd,idir,indlmn,istwf_k,&
& kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
& mpi_enreg%me_g0,natom,nattyp,ngfft,nkpgin,nkpgout,nloalg,nnlout,&
& npwin,npwout,nspinor,ntypat,paw_opt,phkxredin,phkxredout,ph1d,&
& ph3din,ph3dout,signs_,sij,svectout_,pi,ucvol,vectin,vectout_)
#else
 ABI_UNUSED(nnlout)
#endif

 if (choice==1.and.signs==1) then
   if (paw_opt/=3) then
     call dotprod_g(enlout(1),doti,istwf_k,npwin*nspinor,1,vectin,vectout_,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   else
     call dotprod_g(enlout(1),doti,istwf_k,npwin*nspinor,1,vectin,svectout_,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   end if
   ABI_DEALLOCATE(vectout_)
   ABI_DEALLOCATE(svectout_)
 else
   nullify(vectout_,svectout_)
 end if

 if ((cpopt==0).or.(cpopt==1)) then
   iproj=0
   do ispinor=1,nspinor
     iatom=0
     do itypat=1,ntypat
       do ia=1,nattyp(itypat)
         iatom=iatom+1;cprjin(iatom,1)%nlmn=count(indlmn(3,:,itypat)>0)
         do ilmn=1,cprjin(iatom,1)%nlmn
           iproj=iproj+1;cprjin(iatom,1)%cp(:,ilmn)=proj(:,iproj)
         end do
       end do
     end do
   end do
   ABI_DEALLOCATE(proj)
 end if

 DBG_EXIT("COLL")

!Fake statements to satisfy ABI rules
#if ! defined HAVE_GPU_CUDA
 if (.false.) then
   write(std_out,*) atindx1,enl,ffnlin,ffnlout,gprimd
   write(std_out,*) idir,istwf_k,kgin,kgout,kpgin,kpgout
   write(std_out,*) kptin,kptout,lambda,mpi_enreg%me
   write(std_out,*) ngfft,nloalg,ph1d,ph3din,ph3dout
   write(std_out,*) phkxredin,phkxredout,signs,sij
   write(std_out,*) ucvol,vectin
 end if
#endif

end subroutine nonlop_gpu
!!***

end module m_nonlop
!!***
