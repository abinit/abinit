!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ddb_omega_interpol
!! NAME
!!  m_ddb_omega_interpol
!!
!! FUNCTION
!! Interpolate the nonadiabatic second-order susceptibilities 
!! onto a fine frequency grid and incorporate the lattice-mediated contributions.
!!
!! COPYRIGHT
!!  Copyright (C) 2024 ABINIT group (MR and MS)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

module m_ddb_omega_interpol
    
 use defs_basis
 use m_abicore
 use m_profiling_abi
 use m_errors
 use m_ddb
 use m_ddb_magpen,      only : local_spinsus, magmom, mp_d2etot, asrw0
 use m_fstrings,        only : itoa, sjoin
 use m_macroave,        only : POLINT
 use m_io_tools,        only : open_file
 use m_cgtools,         only : fxphas_seq
 use m_dynmat,          only : pheigvec_normalize,phdispl_from_eigvec
 use m_numeric_tools,   only : polcoe

 implicit none

 public :: ddb_omega_interpol ! Perform an interpolation of the derivatives calculated at constrained magnetic moments
                              ! and later on convert them into the physically relevant (spin-relaxed) ones at each value of interpolated omega.

 private

! *************************************************************************

contains 
!!***

!!****f* m_ddb_omega_interpol/ddb_omega_interpol
!! NAME
!! ddb_omega_interpol
!!
!! FUNCTION
!! Interpolate over frequency the secon-order derivatives calculated at 
!! constrained magnetic moments, latter on convert them into susceptibilities 
!! calculated at fixed/relaxed ions/spins.
!!
!! INPUTS
!! ddb (INOUT) = ddb block datastructure
!! magpen = amplitude (in Ha) of the applied magnetic penalty 
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert = maximum number of ipert
!! mpopt = 1 calculate the frozen-magnetic second-order quantities
!!         2 calculate the spin-relaxed second-order quantities 
!! natom= number of atoms in unit cell
!! ntypat= number of atom types
!! ucvol= unit cell volume
!!
!! OUTPUT
!! ddb= ddb%val updated with the corrected second-order derivatives
!!
!! SOURCE

 subroutine ddb_omega_interpol(amu,ddb,ddb_lw,delta_asrw0,delta_asrw0_fm, & 
& dissip,eta,outfilename_radix,magpen,mpatpol,mpdir,mpert,mpopt,natom, &
& nomega,ntypat,omegaflag,omegamax,omegamin,prtvol,rftyp,typat,ucvol,xred)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: dissip,mpert,mpopt,natom,nomega,ntypat,omegaflag,prtvol,rftyp
 real(dp),intent(in) :: eta,magpen,omegamax,omegamin,ucvol
 character(len=*),intent(in) :: outfilename_radix
!arrays
 type(ddb_type),intent(inout) :: ddb,ddb_lw
 integer,intent(in) :: mpatpol(2),mpdir(3),typat(natom)
 real(dp),intent(in) :: amu(ntypat)
 real(dp), intent(inout) :: delta_asrw0(3*natom,3), delta_asrw0_fm(3*natom,3)
 real(dp),intent(in) :: xred(3,natom)

!Local variables -------------------------
!scalars
 integer :: alpha_unit,diel_unit,fs2rs,i,iblok,ifound
 integer :: ii,imode,ipert1,ipert2,iw,j,jblok,jw,kblok,lblok
 integer :: mmag_unit,mmom_unit,mmspec_unit,nblok,ndim 
 integer :: nmat,nmdir,nwcalc,optgb,phon_unit,prtopt
 integer :: locmagsus_unit,zeff_unit,zeffspec_unit,zfield_unit
 real(dp) :: omegastp
 character(len=5000) :: msg,pfmt
 character(len=fnlen) :: alpha_filename
 character(len=fnlen) :: diel_filename,mmag_filename,locmagsus_filename,mmom_filename,mmspec_filename
 character(len=fnlen) :: phon_filename,zeff_filename,zeffspec_filename,zfield_filename
 complex(dpc) :: cplxvar,cplx_weta
 logical :: qeq0
!arrays
 real(dp) :: qphnrm(3),qphon(3,3)
 real(dp), allocatable :: dint_fsddb(:,:),int_fsddb(:,:,:),int_rsddb(:,:,:)
 real(dp), allocatable :: omega(:),omegacalc(:)
 real(dp), allocatable :: w0hessian(:,:),w0berry(:,:)
 real(dp), allocatable :: displ(:),eigvec(:),eigvec_fm(:,:,:,:,:),phfrq(:,:)
 real(dp), allocatable :: mode_phonspec(:,:),phonspec(:)
 real(dp), allocatable :: coeffs(:,:,:)
 complex(dpc), allocatable :: dummysus(:,:)
 complex(dpc), allocatable :: invmagsus(:,:,:), lm_magsus(:,:,:), magsus(:,:,:), invhmat(:,:)
 complex(dpc), allocatable :: dummymom(:,:),dummymom_tr(:,:),mmom(:,:,:), mmom_tr(:,:,:)
 complex(dpc), allocatable :: ri_mmom(:,:,:)
 complex(dpc), allocatable :: lm_zfield(:,:,:),zfield(:,:,:),zfield_tr(:,:)
 complex(dpc), allocatable :: bc_barmagsus(:,:),bc_ss(:,:),bc_sp(:,:)
 complex(dpc), allocatable :: ci_alpha(:,:,:),lm_alpha(:,:,:),lm_alpha_nm(:,:,:,:)
 complex(dpc), allocatable :: ci_localpha(:,:,:),lm_localpha(:,:,:)
 complex(dpc), allocatable :: ci_epsilon(:,:,:),lm_epsilon(:,:,:),lm_epsilon_nm(:,:,:,:)
 complex(dpc), allocatable :: ci_mchi(:,:,:),lm_mchi(:,:,:),lm_mchi_nm(:,:,:,:)
 complex(dpc), allocatable :: modemm(:,:,:),modedisp(:,:,:),modezf(:,:,:)
 complex(dpc), allocatable :: zeff(:,:),zeff_tr(:,:),modemeff(:,:,:),modezeff(:,:,:)
 complex(dpc), allocatable :: fmzeff(:,:),fmzeff_tr(:,:)
 complex(dpc), allocatable :: phongreen(:,:),phongreen_fm(:,:)
 complex(dpc), allocatable :: macmagsus(:,:,:)
 complex(dpc), allocatable :: genzeff_tr(:,:), ri_genelsus(:,:,:)

!TMP: CrI3 varaibles:
 complex(dpc) :: work(4,4)
 complex(dpc),parameter :: ure=(1.d0,0.d0),uim=(0.d0,1.d0)
 real(dp) :: totnorm
 
! *********************************************************************

 write(msg, '(2a,(80a),4a)' ) ch10,('=',ii=1,80),ch10,ch10,&
 ' Omega interpolation of magnetic penalty quantities section ',ch10
 call wrtout([std_out, ab_out], msg)

!Identify the calculated omegas
 nwcalc=ddb%nblok
 ABI_MALLOC(omegacalc,(nwcalc))
 omegacalc(:)=ddb%omega(1,:)

!Define the omega discretization
 omegastp=(omegamax-omegamin)/(nomega-1) 

 prtopt=0
 if (magpen<zero) then
   nmat= 1
 else if (magpen>zero) then
   nmat= mpatpol(2) - mpatpol(1) + 1
 end if
 nmdir=sum(mpdir(:))
 ndim=nmat*nmdir
 fs2rs=1
 qphon=zero
 qphon(:,1)=ddb%qpt(1:3,1)
 qeq0=(sqrt(sum(qphon(:,1)**2))<tol8)
 qphnrm(:)=ddb%nrm(1,1)
 optgb=0
 ABI_MALLOC(omega,(nomega))
 ABI_MALLOC(phfrq,(3*natom,nomega))
 ABI_MALLOC(phonspec,(nomega))
 ABI_MALLOC(phongreen,(3*natom,3*natom))
 ABI_MALLOC(phongreen_fm,(3*natom,3*natom))
 ABI_MALLOC(mode_phonspec,(3*natom,nomega))
 ABI_MALLOC(displ,(2*3*natom*3*natom))
 ABI_MALLOC(eigvec,(2*3*natom*3*natom))
 ABI_MALLOC(eigvec_fm,(2,3,natom,3,natom))
 ABI_MALLOC(modemm,(ndim,3*natom,nomega))
 ABI_MALLOC(modedisp,(3*natom,3*natom,nomega))
 ABI_MALLOC(modezf,(ndim,3*natom,nomega))
 ABI_MALLOC(fmzeff_tr,(3*natom,3))
 ABI_MALLOC(zeff,(3,3*natom))
 ABI_MALLOC(zeff_tr,(3*natom,3))
 ABI_MALLOC(genzeff_tr,(3*natom+ndim,3))
 ABI_MALLOC(ri_genelsus,(3*natom+ndim,3,nomega))
 ABI_MALLOC(ci_alpha,(3,3,nomega))
 ABI_MALLOC(lm_alpha,(3,3,nomega))
 ABI_MALLOC(lm_alpha_nm,(3,3,3*natom,nomega))
 ABI_MALLOC(ci_localpha,(ndim,3,nomega))
 ABI_MALLOC(lm_localpha,(ndim,3,nomega))
 ABI_MALLOC(ci_epsilon,(3,3,nomega))
 ABI_MALLOC(lm_epsilon,(3,3,nomega))
 ABI_MALLOC(lm_epsilon_nm,(3,3,3*natom,nomega))
 ABI_MALLOC(ci_mchi,(3,3,nomega))
 ABI_MALLOC(lm_mchi,(3,3,nomega))
 ABI_MALLOC(lm_mchi_nm,(3,3,3*natom,nomega))
 ABI_MALLOC(modemeff,(3,3*natom,nomega))
 ABI_MALLOC(modezeff,(3,3*natom,nomega))
 ABI_MALLOC(dummysus,(ndim,ndim))
 ABI_MALLOC(magsus,(ndim,ndim,nomega))
 ABI_MALLOC(lm_magsus,(ndim,ndim,nomega))
 ABI_MALLOC(invmagsus,(ndim,ndim,nomega))
 ABI_MALLOC(dummymom,(ndim,(natom+5)*3))
 ABI_MALLOC(dummymom_tr,(ndim,(natom+5)*3))
 ABI_MALLOC(mmom,(ndim,(natom+5)*3,nomega))
 ABI_MALLOC(mmom_tr,((natom+5)*3,ndim,nomega))
 ABI_MALLOC(lm_zfield,(ndim,3,nomega))
 ABI_MALLOC(zfield,(ndim,(natom+5)*3,nomega))
 ABI_MALLOC(zfield_tr,((natom+5)*3,ndim))
 ABI_MALLOC(macmagsus,(3,3,nomega))
 ABI_MALLOC(dint_fsddb,(2,ddb%msize))
 ABI_MALLOC(int_fsddb,(2,ddb%msize,1))
 ABI_MALLOC(int_rsddb,(2,ddb%msize,1))
 ABI_MALLOC(ri_mmom,(ndim,3,nomega))

!For linear interpolation detect the w=0 Hessians and Berry curvatures
 if (omegaflag == 1) then

   ABI_MALLOC(w0hessian,(2,ddb%msize))
   nblok= ddb%nblok
   ifound= 0
   do iblok= 1, nblok
     if (abs(ddb%omega(1,iblok)) < tol12) then
       w0hessian(:,:)= ddb%val_fs(:,:,iblok)
       ifound= 1
     end if
   end do
   if (ifound==0) then
     write(msg, '(3a)' )' No omega=0 block with second-order derivatives', &
   & ' found in the DDB file. This is necessary if omegaflag=1 ',ch10
     ABI_ERROR(msg)
   end if

   ABI_MALLOC(w0berry,(2,ddb_lw%msize))
   nblok= ddb_lw%nblok
   ifound= 0
   do iblok= 1, nblok
     if (abs(ddb_lw%omega(1,iblok)) < tol12) then
       w0berry(:,:)= ddb_lw%val_fs(:,:,iblok)
       ifound= 1
     end if
   end do
   if (ifound==0) then
     write(msg, '(3a)' )' No omega=0 block with third-order derivatives', &
   & ' found in the DDB file. This is necessary if omegaflag=1 ',ch10
     ABI_ERROR(msg)
   end if
 end if

!For nwcalc-1 Taylor-expansion interpolation precalculate the coefficients
 if (omegaflag == 3) then
   ABI_MALLOC(coeffs,(2,nwcalc,ddb%msize))
   do ii=1,ddb%msize
     if (all(ddb%flg(ii,:)==1)) then
       call polcoe(omegacalc,ddb%val_fs(1,ii,:),nwcalc,coeffs(1,:,ii))
       call polcoe(omegacalc,ddb%val_fs(2,ii,:),nwcalc,coeffs(2,:,ii))
     end if
   end do
 end if

!Loop over the frequency
 do iw=1,nomega
   omega(iw)=omegamin+omegastp*(iw-1)
   if (nomega==1) omega(iw)=omegamin

   !Perform the different interpolations
   !Lineal (with analytic Berry curvature) with dissipation if eta/=0
   if (omegaflag==1) then
     do ii=1,ddb%msize
       if (all(ddb%flg(ii,:)==1)) then
         int_fsddb(1,ii,1)= w0hessian(1,ii) + omega(iw)*w0berry(1,ii) - eta*w0berry(2,ii)
         int_fsddb(2,ii,1)= w0hessian(2,ii) + omega(iw)*w0berry(2,ii) + eta*w0berry(1,ii)
       end if
     end do

   !Polynomial with no dissipation 
   else if (omegaflag==2) then
     do ii=1,ddb%msize
       if (all(ddb%flg(ii,:)==1)) then
         call POLINT(omegacalc,ddb%val_fs(1,ii,:),nwcalc,omega(iw),int_fsddb(1,ii,1),dint_fsddb(1,ii)) 
         call POLINT(omegacalc,ddb%val_fs(2,ii,:),nwcalc,omega(iw),int_fsddb(2,ii,1),dint_fsddb(2,ii)) 
       else if (count(ddb%flg(ii,:)==0)/=nwcalc) then
         write(msg,'(a,a,a)')&
         'ddb_omega_interpol detects differences between the DDB bloks for each frequency.',ch10,&
       & ' The interpolation has been stopped.' 
         ABI_ERROR(msg)
       end if
     end do 

   !Taylor-expansion around w=0 with dissipation if eta/=0
   else if (omegaflag==3) then
     cplx_weta=cmplx(omega(iw),eta,16)
     do ii=1,ddb%msize
       if (all(ddb%flg(ii,:)==1)) then
         cplxvar=cmplx(zero,zero,16)
         do jw= 1, nwcalc
           cplxvar= cplxvar + cplx_weta**(jw-1)*cmplx(coeffs(1,jw,ii),coeffs(2,jw,ii),16)
         end do
         int_fsddb(1,ii,1)=real(cplxvar)
         int_fsddb(2,ii,1)=aimag(cplxvar) 
       else if (count(ddb%flg(ii,:)==0)/=nwcalc) then
         write(msg,'(a,a,a)')&
         'ddb_omega_interpol detects differences between the DDB bloks for each frequency.',ch10,&
       & ' The interpolation has been stopped.' 
         ABI_ERROR(msg)
       end if
     end do 
   end if

   !Calculate the local spin susceptibilities
   call local_spinsus(dummysus,ddb,1,dummysus,invmagsus(:,:,iw),&
 & dummysus,magpen,magsus(:,:,iw),mpatpol,mpdir,mpert,natom,1,ndim,nmdir,prtopt,prtvol,&
 & fs2rs=fs2rs,blkval_fs=int_fsddb)

   !Calculate the 1st-order magnetic moments
   call magmom(dummymom,dummymom_tr,ddb,dummysus,dummysus,1,1,1,magpen,&
 & magsus(:,:,iw),mmom(:,:,iw),mmom_tr(:,:,iw),mpatpol,mpdir,mpert,natom,&
 & 1,ndim,nmdir,prtopt,prtvol,qphon,xred,zfield(:,:,iw),zfield_tr,&
 & fs2rs=fs2rs,blkval_fs=int_fsddb)
  
   !Now calculate the non-magnetic second-order quantities
   call ddb%to_d2etot(int_fsddb,1,0,qeq0,qphon,qphnrm,ucvol,optgb,omega=omega(iw))

   call mp_d2etot(dummysus,ddb,1,dummysus,magsus(:,:,iw),&
 & magpen,mpert,mpopt,natom,1,ndim,qphon,xred,zfield(:,:,iw),zfield_tr, &
 & fs2rs=fs2rs,blkval_fs=int_fsddb,blkval_rs=int_rsddb)

   call ddb%to_d2etot(int_fsddb,1,1,qeq0,qphon,qphnrm,ucvol,optgb,omega=omega(iw))
   if (mpopt==2) call ddb%to_d2etot(int_rsddb,1,1,qeq0,qphon,qphnrm,ucvol,optgb,omega=omega(iw))

   !Calculate the phonon propagator (Green's function) and spectral function
   !and the lattice-mediated contributions to the different susceptibilities.
   if (mpopt==1) then
     call phonon_green(amu,displ,eigvec,eta,int_fsddb,& 
   & mode_phonspec(:,iw),mpert,natom,ntypat,omega(iw),&
   & phfrq(:,iw),phongreen,phonspec(iw),typat)

     call ri_d2etot(int_fsddb,ci_alpha(:,:,iw),ci_epsilon(:,:,iw),ci_localpha(:,:,iw),ci_mchi(:,:,iw),&
   & lm_alpha(:,:,iw),lm_epsilon(:,:,iw),lm_localpha(:,:,iw),lm_magsus(:,:,iw),lm_mchi(:,:,iw),&
   & magsus(:,:,iw),mpert,mmom(:,:,iw),mmom_tr(:,:,iw),natom,ndim,phongreen,ucvol,zeff)
   else if (mpopt==2) then
     call phonon_green(amu,displ,eigvec,eta,int_rsddb,& 
   & mode_phonspec(:,iw),mpert,natom,ntypat,omega(iw),&
   & phfrq(:,iw),phongreen,phonspec(iw),typat)

     call ri_d2etot(int_rsddb,ci_alpha(:,:,iw),ci_epsilon(:,:,iw),ci_localpha(:,:,iw),ci_mchi(:,:,iw),&
   & lm_alpha(:,:,iw),lm_epsilon(:,:,iw),lm_localpha(:,:,iw),lm_magsus(:,:,iw),lm_mchi(:,:,iw),& 
   & magsus(:,:,iw),mpert,mmom(:,:,iw),mmom_tr(:,:,iw),natom,ndim,phongreen,ucvol,zeff)

     call lm_normal_modes(int_rsddb,displ,eta,lm_alpha_nm(:,:,:,iw),lm_epsilon_nm(:,:,:,iw),lm_mchi_nm(:,:,:,iw), &
   & mmom(:,:,iw),modemm(:,:,iw),modedisp(:,:,iw),modemeff(:,:,iw),modezeff(:,:,iw),modezf(:,:,iw),&
   & mpert,natom,ndim,omega(iw),phfrq(:,iw),ucvol,zfield(:,:,iw))

   end if

 end do

!!!  Print results of interpolation
!Local magnetic susceptibilities
 locmagsus_filename=trim(outfilename_radix)//"_LOCMAGSUS"
 if (open_file(locmagsus_filename, msg, newunit=locmagsus_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(locmagsus_unit,*) '#'
 write(locmagsus_unit,*) '#  Local magnetic susceptibilities calculated and interpolated by ANADDB'
 write(locmagsus_unit,*) '#'
 write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' )  ndim**2

 write(locmagsus_unit,*) '#  Real part of clamped-ion local magnetic susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(locmagsus_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(magsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(locmagsus_unit,msg,'COLL')
 end do

 write(locmagsus_unit,*) ' '
 write(locmagsus_unit,*) '#  Imaginary part of clamped-ion local magnetic susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(locmagsus_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(magsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(locmagsus_unit,msg,'COLL')
 end do

 write(locmagsus_unit,*) ' '
 write(locmagsus_unit,*) '#  Real part of relaxed-ion local magnetic susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(locmagsus_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(magsus(i,j,iw)+lm_magsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(locmagsus_unit,msg,'COLL')
 end do

 write(locmagsus_unit,*) ' '
 write(locmagsus_unit,*) '#  Imaginary part of relaxed-ion local magnetic susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(locmagsus_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(magsus(i,j,iw)+lm_magsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(locmagsus_unit,msg,'COLL')
 end do

 write(locmagsus_unit,*) ' '
 write(locmagsus_unit,*) '#  Real part of the inverse of the clamped-ion local magnetic susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X^{-1}_11     X^{-1}_12     ...     X^{-1}_21     X^{-1}_22     ...'
 call wrtout(locmagsus_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(invmagsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(locmagsus_unit,msg,'COLL')
 end do

 write(locmagsus_unit,*) ' '
 write(locmagsus_unit,*) '#  Imaginary part of the inverse of the clamped-ion local magnetic susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X^{-1}_11     X^{-1}_12     ...     X^{-1}_21     X^{-1}_22     ...'
 call wrtout(locmagsus_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(invmagsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(locmagsus_unit,msg,'COLL')
 end do

 close (locmagsus_unit)

!Zfields
 zfield_filename=trim(outfilename_radix)//"_LOCZFIELDS"
 if (open_file(zfield_filename, msg, newunit=zfield_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  ndim
 do imode= 1, 3*natom
   write(zfield_unit,*) ' '
   write(zfield_unit,'(a,i3)') '#  Real part of local Zeeman fields (at. units) induced by phonon mode:', imode
   write(msg,'(a,a,a)') ch10,&
 &           ' # At  hw     Z_{mat_1,1}     Z_{mat_1,2}     ...     Z_{mat_2,1}     Z_{mat_2,2}'
   call wrtout(zfield_unit,msg,'COLL')
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (real(modezf(i,imode,iw)),i=1,ndim)
     call wrtout(zfield_unit,msg,'COLL')
   end do
   write(zfield_unit,*) ' '
   write(zfield_unit,'(a,i3)') '#  Imaginary part of local Zeeman fields (at. units) induced by phonon mode:', imode
   write(msg,'(a,a,a)') ch10,&
 &           ' # At  hw     Z_{mat_1,1}     Z_{mat_1,2}     ...     Z_{mat_2,1}     Z_{mat_2,2}'
   call wrtout(zfield_unit,msg,'COLL')
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (aimag(modezf(i,imode,iw)),i=1,ndim)
     call wrtout(zfield_unit,msg,'COLL')
   end do
 end do

 write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' )  ndim*3
 write(zfield_unit,*) ' '
 write(zfield_unit,*) '#  Real part of clamped-ion local Zeeman fields induced by electric field (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     Z_11     Z_12      Z_13    ...     Z_21     Z_22     ...'
 call wrtout(zfield_unit,msg,'COLL')
 do iw=1,nomega
   write(msg,pfmt) omega(iw), ((real(zfield(i,(natom+1)*3+j,iw)),j=1,3),i=1,ndim)
   call wrtout(zfield_unit,msg,'COLL')
 end do

 write(zfield_unit,*) ' '
 write(zfield_unit,*) '#  Imag part of clamped-ion local Zeeman fields induced by electric field (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     Z_11     Z_12      Z_13    ...     Z_21     Z_22     ...'
 call wrtout(zfield_unit,msg,'COLL')
 do iw=1,nomega
   write(msg,pfmt) omega(iw), ((aimag(zfield(i,(natom+1)*3+j,iw)),j=1,3),i=1,ndim)
   call wrtout(zfield_unit,msg,'COLL')
 end do

 write(zfield_unit,*) ' '
 write(zfield_unit,*) '#  Real part of clamped-ion local Zeeman fields induced by macroscopic Zeeman field (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     Z_11     Z_12      Z_13    ...     Z_21     Z_22     ...'
 call wrtout(zfield_unit,msg,'COLL')
 do iw=1,nomega
   write(msg,pfmt) omega(iw), ((real(zfield(i,(natom+4)*3+j,iw)),j=1,3),i=1,ndim)
   call wrtout(zfield_unit,msg,'COLL')
 end do

 write(zfield_unit,*) ' '
 write(zfield_unit,*) '#  Imag part of clamped-ion local Zeeman fields induced by macroscopic Zeeman field (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     Z_11     Z_12      Z_13    ...     Z_21     Z_22     ...'
 call wrtout(zfield_unit,msg,'COLL')
 do iw=1,nomega
   write(msg,pfmt) omega(iw), ((aimag(zfield(i,(natom+4)*3+j,iw)),j=1,3),i=1,ndim)
   call wrtout(zfield_unit,msg,'COLL')
 end do
 close (zfield_unit)

!Magnetic moments
 mmom_filename=trim(outfilename_radix)//"_LOCMAGMOM"
 if (open_file(mmom_filename, msg, newunit=mmom_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(mmom_unit,*) '#'
 write(mmom_unit,*) '#  Local magnetic moments calculated and interpolated by ANADDB'
 write(mmom_unit,*) '#'

 write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  ndim
 do imode= 1, 3*natom
   write(mmom_unit,*) ' '
   write(mmom_unit,'(a,i3)') '#  Real part of local magnetic moments (at. units) induced by phonon mode:', imode
   write(msg,'(a,a,a)') ch10,&
 &           ' # At  hw     m_{mat_1,1}     m_{mat_1,2}     ...     m_{mat_2,1}     m_{mat_2,2}'
   call wrtout(mmom_unit,msg,'COLL')
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (real(modemm(i,imode,iw)),i=1,ndim)
     call wrtout(mmom_unit,msg,'COLL')
   end do
   write(mmom_unit,*) ' '
   write(mmom_unit,'(a,i3)') '#  Imaginary part of local magnetic moments (at. units) induced by phonon mode:', imode
   write(msg,'(a,a,a)') ch10,&
 &           ' # At  hw     m_{mat_1,1}     m_{mat_1,2}     ...     m_{mat_2,1}     m_{mat_2,2}'
   call wrtout(mmom_unit,msg,'COLL')
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (aimag(modemm(i,imode,iw)),i=1,ndim)
     call wrtout(mmom_unit,msg,'COLL')
   end do
 end do

 write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  ndim*3
 write(mmom_unit,*) ' '
 write(mmom_unit,*) '#  Real part of the clamped-ion local magnetoelectric tensor(at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(ci_localpha(i,j,iw)),j=1,3),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

 write(mmom_unit,*) ' '
 write(mmom_unit,*) '#  Imaginary part of the clamped-ion local magnetoelectric tensor(at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(ci_localpha(i,j,iw)),j=1,3),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

 write(mmom_unit,*) ' '
 write(mmom_unit,*) '#  Real part of relaxed-ion local magnetoelectric tensor (at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(ci_localpha(i,j,iw)+lm_localpha(i,j,iw)),j=1,3),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

 write(mmom_unit,*) ' '
 write(mmom_unit,*) '#  Imaginary part of relaxed-ion local magnetoelectric tensor (at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(ci_localpha(i,j,iw)+lm_localpha(i,j,iw)),j=1,3),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

 close(mmom_unit)

!Dielectric susceptibility
 diel_filename=trim(outfilename_radix)//"_DIELTENS"
 if (open_file(diel_filename, msg, newunit=diel_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 if (mpopt==1) then
   write(diel_unit,*) '#  Fixed-spin dielectric tensor calculated and interpolated by ANADDB'
 else if (mpopt==2) then
   write(diel_unit,*) '#  Relaxed-spin dielectric tensor calculated and interpolated by ANADDB'
 else
   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
   ABI_ERROR(msg)
 end if
 
 write(diel_unit,*) '#'
 write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' ) 9 

 write(diel_unit,*) '#  Real part of clamped-ion dielectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
 call wrtout(diel_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(ci_epsilon(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do

 write(diel_unit,*) ' '
 write(diel_unit,*) '#  Imaginary part of clamped-ion dielectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
 call wrtout(diel_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(ci_epsilon(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do

 write(diel_unit,*) ' '
 write(diel_unit,*) '#  Real part of relaxed-ion dielectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
 call wrtout(diel_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(lm_epsilon(i,j,iw)+ci_epsilon(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do

 write(diel_unit,*) ' '
 write(diel_unit,*) '#  Imaginary part of relaxed-ion dielectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
 call wrtout(diel_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(lm_epsilon(i,j,iw)+ci_epsilon(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do
 write(diel_unit,*) ' '

 write(diel_unit,*) '#'
 write(diel_unit,*) '# Phonon modes contribution to dielectric tensor' 
 write(diel_unit,*) '#'

 do imode= 1, 3*natom
   write(diel_unit,*) ' '
   write(diel_unit,*) '#  Real part of dielectric tensor due to phonon mode:', imode
   write(msg,'(a,a)') ch10,&
 &           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
   call wrtout(diel_unit,msg,'COLL')
   do iw=1,nomega
      write(msg,pfmt) &
   &  omega(iw), ((real(lm_epsilon_nm(i,j,imode,iw)),j=1,3),i=1,3)
      call wrtout(diel_unit,msg,'COLL')
   end do
   write(diel_unit,*) ' '
   write(diel_unit,*) '#  Imaginary part of dielectric tensor due to phonon mode:', imode
   write(msg,'(a,a)') ch10,&
 &           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
   call wrtout(diel_unit,msg,'COLL')
   do iw=1,nomega
      write(msg,pfmt) &
   &  omega(iw), ((aimag(lm_epsilon_nm(i,j,imode,iw)),j=1,3),i=1,3)
      call wrtout(diel_unit,msg,'COLL')
   end do
 end do 

 close(diel_unit)

!Magnetoelectric tensor
 alpha_filename=trim(outfilename_radix)//"_MAGNETOELTENS"
 if (open_file(alpha_filename, msg, newunit=alpha_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 if (mpopt==1) then
   write(alpha_unit,*) '#  Fixed-spin magnetoelectric tensor calculated and interpolated by ANADDB'
 else if (mpopt==2) then
   write(alpha_unit,*) '#  Relaxed-spin magnetoelectric tensor calculated and interpolated by ANADDB'
 else
   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
   ABI_ERROR(msg)
 end if
 
 write(alpha_unit,*) '#'
 write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' ) 9 

 write(alpha_unit,*) '#  Real part of clamped-ion magnetoelectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     alpha_11     alpha_12     ...     alpha_21     alpha_22     ...'
 call wrtout(alpha_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(ci_alpha(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do

 write(alpha_unit,*) ' '
 write(alpha_unit,*) '#  Imaginary part of clamped-ion magnetoelectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     alpha_11     alpha_12     ...     alpha_21     alpha_22     ...'
 call wrtout(alpha_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(ci_alpha(i,j,iw)),j=1,3),i=1,3)
    call wrtout(alpha_unit,msg,'COLL')
 end do

 write(alpha_unit,*) ' '
 write(alpha_unit,*) '#  Real part of relaxed-ion magnetoelectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     alpha_11     alpha_12     ...     alpha_21     alpha_22     ...'
 call wrtout(alpha_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(lm_alpha(i,j,iw)+ci_alpha(i,j,iw)),j=1,3),i=1,3)
    call wrtout(alpha_unit,msg,'COLL')
 end do

 write(alpha_unit,*) ' '
 write(alpha_unit,*) '#  Imaginary part of relaxed-ion magnetoelectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     alpha_11     alpha_12     ...     alpha_21     alpha_22     ...'
 call wrtout(alpha_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(lm_alpha(i,j,iw)+ci_alpha(i,j,iw)),j=1,3),i=1,3)
    call wrtout(alpha_unit,msg,'COLL')
 end do
 write(alpha_unit,*) ''

 write(alpha_unit,*) '#'
 write(alpha_unit,*) '# Phonon modes contribution to magnetoelectric tensor' 
 write(alpha_unit,*) '#'

! do iw=1,nomega
!   write(200,*) omega(iw),(real(lm_alpha_nm(1,1,i,iw)),i=1,3*natom)
!   write(201,*) omega(iw),(aimag(lm_alpha_nm(1,1,i,iw)),i=1,3*natom)
! end do

 do imode= 1, 3*natom
   write(alpha_unit,*) ' '
   write(alpha_unit,*) '#  Real part of magnetoelectric tensor due to phonon mode:', imode
   write(msg,'(a,a)') ch10,&
 &           ' # At  hw     alpha_11     alpha_12     ...     alpha_21     alpha_22     ...'
   call wrtout(alpha_unit,msg,'COLL')
   do iw=1,nomega
      write(msg,pfmt) &
   &  omega(iw), ((real(lm_alpha_nm(i,j,imode,iw)),j=1,3),i=1,3)
      call wrtout(alpha_unit,msg,'COLL')
   end do
   write(alpha_unit,*) ' '
   write(alpha_unit,*) '#  Imaginary part of magnetoelectric tensor due to phonon mode:', imode
   write(msg,'(a,a)') ch10,&
 &           ' # At  hw     alpha_11     alpha_12     ...     alpha_21     alpha_22     ...'
   call wrtout(alpha_unit,msg,'COLL')
   do iw=1,nomega
      write(msg,pfmt) &
   &  omega(iw), ((aimag(lm_alpha_nm(i,j,imode,iw)),j=1,3),i=1,3)
      call wrtout(alpha_unit,msg,'COLL')
   end do
 end do 

 close(alpha_unit)

!Magnetic susceptibility
 mmag_filename=trim(outfilename_radix)//"_MAGSUS"
 if (open_file(mmag_filename, msg, newunit=mmag_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 if (mpopt==1) then
   write(mmag_unit,*) '#  Fixed-spin magnetic susceptibility calculated and interpolated by ANADDB'
 else if (mpopt==2) then
   write(mmag_unit,*) '#  Relaxed-spin magnetic susceptibility calculated and interpolated by ANADDB'
 else
   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
   ABI_ERROR(msg)
 end if
 
 write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' ) 9 

 write(mmag_unit,*) ' '
 write(mmag_unit,*) '#  Real part of clamped-ion magnetic susceptibility (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(mmag_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(ci_mchi(i,j,iw)),j=1,3),i=1,3)
    call wrtout(mmag_unit,msg,'COLL')
 end do

 write(mmag_unit,*) ' '
 write(mmag_unit,*) '#  Imaginary part of clamped-ion magnetic susceptibility (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(mmag_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(ci_mchi(i,j,iw)),j=1,3),i=1,3)
    call wrtout(mmag_unit,msg,'COLL')
 end do

 write(mmag_unit,*) ' '
 write(mmag_unit,*) '#  Real part of relaxed-ion magnetic susceptibility (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(mmag_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(ci_mchi(i,j,iw)+lm_mchi(i,j,iw)),j=1,3),i=1,3)
    call wrtout(mmag_unit,msg,'COLL')
 end do

 write(mmag_unit,*) ' '
 write(mmag_unit,*) '#  Imaginary part of relaxed-ion magnetic susceptibility (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(mmag_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(ci_mchi(i,j,iw)+lm_mchi(i,j,iw)),j=1,3),i=1,3)
    call wrtout(mmag_unit,msg,'COLL')
 end do

 write(mmag_unit,*) ''

 write(mmag_unit,*) '#'
 write(mmag_unit,*) '# Phonon modes contribution to magnetoelectric tensor' 
 write(mmag_unit,*) '#'

 do imode= 1, 3*natom
   write(mmag_unit,*) ' '
   write(mmag_unit,*) '#  Real part of magnetic susceptibility due to phonon mode:', imode
   write(msg,'(a,a)') ch10,&
 &           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
   call wrtout(mmag_unit,msg,'COLL')
   do iw=1,nomega
      write(msg,pfmt) &
   &  omega(iw), ((real(lm_mchi_nm(i,j,imode,iw)),j=1,3),i=1,3)
      call wrtout(mmag_unit,msg,'COLL')
   end do
   write(mmag_unit,*) ' '
   write(mmag_unit,*) '#  Imaginary part of magnetic susceptibility due to phonon mode:', imode
   write(msg,'(a,a)') ch10,&
 &           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
   call wrtout(mmag_unit,msg,'COLL')
   do iw=1,nomega
      write(msg,pfmt) &
   &  omega(iw), ((aimag(lm_mchi_nm(i,j,imode,iw)),j=1,3),i=1,3)
      call wrtout(mmag_unit,msg,'COLL')
   end do
 end do 

 close(mmag_unit)

!Phonon spectral function
 phon_filename=trim(outfilename_radix)//"_SPECTRAL_PHONON"
 if (open_file(phon_filename, msg, newunit=phon_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(phon_unit,*) '#'
 if (mpopt==1) then
   write(phon_unit,*) '#  Fixed-spin phonon spectral function calculated and interpolated by ANADDB'
 else if (mpopt==2) then
   write(phon_unit,*) '#  Relaxed-spin phonon spectral function calculated and interpolated by ANADDB'
 else
   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
   ABI_ERROR(msg)
 end if
 
 write(msg,'(a,a)') ch10, ' # At  hw               Phonon SF      '
 call wrtout(phon_unit,msg,'COLL')
 do iw=1,nomega
   write(msg,*) omega(iw), phonspec(iw)
   call wrtout(phon_unit,msg,'COLL')
 end do

 close(phon_unit)

!Phonon frequencies
 phon_filename=trim(outfilename_radix)//"_PHFRQ"
 if (open_file(phon_filename, msg, newunit=phon_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(phon_unit,*) '#'
 if (mpopt==1) then
   write(phon_unit,*) '#  Frozen-spin phonon frequencies calculated and interpolated by ANADDB'
 else if (mpopt==2) then
   write(phon_unit,*) '#  Relaxed-Spin phonon frequencies calculated and interpolated by ANADDB'
 else
   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
   ABI_ERROR(msg)
 end if

 write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' ) natom*3
 write(msg,'(a,a)') ch10,&
&           ' # At  hw    eval(1)     eval(2) ...'
 call wrtout(phon_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) omega(iw), phfrq(:,iw)
    call wrtout(phon_unit,msg,'COLL')
 end do
 
 close(phon_unit)

!Phonon eigendisplacements
 phon_filename=trim(outfilename_radix)//"_PHDISP"
 if (open_file(phon_filename, msg, newunit=phon_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(phon_unit,*) '#'
 if (mpopt==1) then
   write(phon_unit,*) '#  Frozen-spin phonon eigendisplacements calculated and interpolated by ANADDB'
 else if (mpopt==2) then
   write(phon_unit,*) '#  Relaxed-spin phonon eigendisplacements calculated and interpolated by ANADDB'
 else
   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
   ABI_ERROR(msg)
 end if

 write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  3*natom
 do imode= 1, 3*natom
   write(phon_unit,*) ' '
   write(phon_unit,'(a,i3)') '#  Real part of phonon eigendisplacement:', imode
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (real(modedisp(i,imode,iw)),i=1,3*natom)
     call wrtout(phon_unit,msg,'COLL')
   end do

   write(phon_unit,*) ' '
   write(phon_unit,'(a,i3)') '#  Imaginary part of phonon eigendisplacement:', imode
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (aimag(modedisp(i,imode,iw)),i=1,3*natom)
     call wrtout(phon_unit,msg,'COLL')
   end do
 end do

!Born effective charges
 zeff_filename=trim(outfilename_radix)//"_ZEFF"
 if (open_file(zeff_filename, msg, newunit=zeff_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(zeff_unit,*) '#'
 write(zeff_unit,*) '#  Electric Born effective charges calculated and interpolated by ANADDB'
 write(zeff_unit,*) '#'

 write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' ) 3 
 do imode= 1, 3*natom
   write(zeff_unit,*) ' '
   write(zeff_unit,'(a,i3)') '#  Real part of electric Born charge (at. units) induced by phonon mode:', imode
   write(msg,'(a,a)') ch10,&
 &           ' # At  hw     Z^x_{n}     Z^y_{n}     Z^z_{n}'
   call wrtout(zeff_unit,msg,'COLL')
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (real(modezeff(i,imode,iw)),i=1,3)
     call wrtout(zeff_unit,msg,'COLL')
   end do
   write(zeff_unit,*) ' '
   write(zeff_unit,'(a,i3)') '#  Imaginary part of electric Born charge (at. units) induced by phonon mode:', imode
   write(msg,'(a,a)') ch10,&
 &           ' # At  hw     Z^x_{n}     Z^y_{n}     Z^z_{n}'
   call wrtout(zeff_unit,msg,'COLL')
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (aimag(modezeff(i,imode,iw)),i=1,3)
     call wrtout(zeff_unit,msg,'COLL')
   end do
 end do

 do imode= 1, 3*natom
   write(zeff_unit,*) ' '
   write(zeff_unit,'(a,i3)') '#  Real part of magnetic Born charge (at. units) induced by phonon mode:', imode
   write(msg,'(a,a)') ch10,&
 &           ' # At  hw     M^x_{n}     M^y_{n}     M^z_{n}'
   call wrtout(zeff_unit,msg,'COLL')
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (real(modemeff(i,imode,iw)),i=1,3)
     call wrtout(zeff_unit,msg,'COLL')
   end do
   write(zeff_unit,*) ' '
   write(zeff_unit,'(a,i3)') '#  Imaginary part of magnetic Born charge (at. units) induced by phonon mode:', imode
   write(msg,'(a,a)') ch10,&
 &           ' # At  hw     M^x_{n}     M^y_{n}     M^z_{n}'
   call wrtout(zeff_unit,msg,'COLL')
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (aimag(modemeff(i,imode,iw)),i=1,3)
     call wrtout(zeff_unit,msg,'COLL')
   end do
 end do

 close(zeff_unit)

 ABI_FREE(dint_fsddb)
 ABI_FREE(int_fsddb)
 ABI_FREE(int_rsddb)
 ABI_FREE(dummysus)
 ABI_FREE(magsus)
 ABI_FREE(lm_magsus)
 ABI_FREE(invmagsus)
 ABI_FREE(mmom)
 ABI_FREE(mmom_tr)
 ABI_FREE(zfield)
 ABI_FREE(zfield_tr)
 ABI_FREE(ci_epsilon)
 ABI_FREE(lm_epsilon)
 ABI_FREE(lm_epsilon_nm)
 ABI_FREE(ci_mchi)
 ABI_FREE(lm_mchi)
 ABI_FREE(ci_alpha)
 ABI_FREE(lm_alpha)
 ABI_FREE(lm_alpha_nm)
 ABI_FREE(ci_localpha)
 ABI_FREE(lm_localpha)
 ABI_FREE(macmagsus)
 ABI_FREE(omega)
 ABI_FREE(phfrq)
 ABI_FREE(phongreen)
 ABI_FREE(phongreen_fm)
 ABI_FREE(phonspec)
 ABI_FREE(mode_phonspec)
 ABI_FREE(displ)
 ABI_FREE(eigvec)
 ABI_FREE(eigvec_fm)
 ABI_FREE(modemm)
 ABI_FREE(modedisp)
 ABI_FREE(modezf)
 ABI_FREE(fmzeff_tr)
 ABI_FREE(zeff)
 ABI_FREE(zeff_tr)
 ABI_FREE(genzeff_tr)
 ABI_FREE(ri_genelsus)
 ABI_FREE(modezeff)
 ABI_FREE(modemeff)
 ABI_SFREE(w0hessian)
 ABI_SFREE(w0berry)
 ABI_SFREE(coeffs)

 end subroutine ddb_omega_interpol
!!***

!!****f* ABINIT/phonon_green
!! NAME
!!  phonon_green
!!
!! FUNCTION
!!  Computes the phonon Green's function and spectral function at 
!!  a given value of frequency and imaginary damping
!!
!! COPYRIGHT
!!  Copyright (C) 2024 ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amu(ntypat)= atomic masses
!!  ifc(3*natom,3*natom)= Interatomic-force constants calculated at a given omega
!!  natom= number of atoms in the cell
!!  ntypat= number of atom types in the cell
!!  omega= frequency at which the IFCs have been calculated
!!  typat(natom)= array with the type of atoms in the cell
!!
!! OUTPUT
!!  phonspec= phonon spectral function at the input omega
!!  phfrq= phonon frequencies calculated with the IFCs at the input omega
!!  phfrq= phonon eigenvectors calculated with the IFCs at the input omega
!!  
!!
!! SIDE EFFECTS
!!
!! NOTES
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


subroutine phonon_green(amu,displ,eigvec,eta,blkval,& 
& mode_phonspec,mpert,natom,ntypat,omega,&
& phfrq,phongreen,phonspec,typat)

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: mpert,natom,ntypat 
 real(dp), intent(in) :: eta,omega
 real(dp), intent(out) :: phonspec
!arrays
 integer, intent(in)  :: typat(natom)
 real(dp), intent(in) :: amu(ntypat)
 real(dp), intent(in)  :: blkval(2,3,mpert,3,mpert,1)
 real(dp), intent(out) :: displ(2*3*natom*3*natom)
 real(dp), intent(out) :: eigvec(2*3*natom*3*natom)
 real(dp), intent(out) :: phfrq(3*natom)
 real(dp), intent(out) :: mode_phonspec(3*natom)
 complex(dpc), intent(out) :: phongreen(3*natom,3*natom)

!Local variables-------------------------------
!scalars
 integer :: iat1,iat2,idir1,idir2,ipert1,ipert2
 integer :: icol,ier,imode,info,irow,lwork,mpdim,pdim
 real(dp) :: mfac1, mfac2
 complex(dpc) :: cplx_eta
!arrays
 integer, allocatable :: ipiv(:)
 real(dp) :: dum(2,0) 
 real(dp), allocatable :: invmassfac(:,:)
 real(dp), allocatable :: matrx(:,:),zhpev1(:,:),zhpev2(:)
 real(dp), allocatable :: eigval(:)
 real(dp), allocatable, save :: delta_asrw0(:,:)
 complex(dpc), allocatable :: ifc(:,:)
 complex(dpc), allocatable :: dynmat(:,:),ifc_w2mass(:,:)
 complex(dpc),allocatable :: work(:),work1(:,:)
 complex(dpc),allocatable :: mass_phongreen(:,:)
!character(len=500) :: msg                   

! *************************************************************************

 DBG_ENTER("COLL")

!Extract the IFCs.
 pdim=3*natom
 ABI_MALLOC(ifc,(pdim,pdim))
 do ipert2= 1, natom
   do idir2= 1, 3
     icol= (ipert2-1)*3 + idir2
     do ipert1= 1, natom
       do idir1= 1, 3
         irow= (ipert1-1)*3 + idir1
         ifc(irow,icol)= &
       & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,1), &
       & blkval(2,idir1,ipert1,idir2,ipert2,1),16)
       end do
     end do
   end do
 end do

!Apply ASR: it has weird consequences on the intensities of the spectral function
!better not applied.
! ABI_MALLOC_IFNOT(delta_asrw0,(3*natom,3))
! if (omega < tol14) then
!   call asrw0(delta_asrw0,ifc,natom,0) 
! else 
!   call asrw0(delta_asrw0,ifc,natom,1) 
! end if

!Build an array with the inverse mass factors
 ABI_MALLOC(invmassfac,(natom,natom))
 do iat2= 1, natom
   do iat1= 1, natom
     invmassfac(iat1,iat2)=one/sqrt(amu(typat(iat1))*amu(typat(iat2)))/amu_emass
   end do
 end do

!Build the dynamical and (Phi-M(w+eta)**2) matrices
 ABI_MALLOC(dynmat,(pdim,pdim))
 ABI_MALLOC(ifc_w2mass,(pdim,pdim))

 cplx_eta=cmplx(0.0_dp,eta)
 do iat2= 1, natom
   do idir2= 1, 3
     icol= (iat2-1)*3 + idir2
     do iat1= 1, natom
       do idir1= 1, 3
         irow= (iat1-1)*3 + idir1
         dynmat(irow,icol)= invmassfac(iat1,iat2)*ifc(irow,icol)
         ifc_w2mass(irow,icol)= ifc(irow,icol)
         if (irow==icol) then
           ifc_w2mass(irow,icol)= ifc_w2mass(irow,icol) - &
         & amu(typat(iat1))*amu_emass*(omega+cplx_eta)**2 
         end if
       end do
     end do
   end do
 end do

!Invert to obtain the phonon Green's function
 ABI_MALLOC(work1,(pdim,pdim))
 work1=ifc_w2mass

 ABI_MALLOC(ipiv,(pdim))
 call zgetrf( pdim, pdim, work1, pdim, ipiv, info )
 ABI_CHECK(info == 0, sjoin('zgetrf returned:', itoa(info)))

 ABI_MALLOC(work,(2))
 call zgetri( pdim, work1, pdim, ipiv, work, -1, info )
 ABI_CHECK(info == 0, sjoin('zgetri returned:', itoa(info)))
 lwork=int(work(1))

 ABI_REMALLOC(work,(lwork))
 call zgetri( pdim, work1, pdim, ipiv, work, lwork, info )
 ABI_CHECK(info == 0, sjoin('zgetri returned:', itoa(info)))

 phongreen=-work1

!Now apply the mass factors
 ABI_MALLOC(mass_phongreen,(pdim,pdim))
 do icol= 1, pdim
   iat2= ceiling(icol/three)
   mfac2= sqrt(amu(typat(iat2))*amu_emass)
   do irow= 1, pdim
     iat1= ceiling(irow/three)
     mfac1= sqrt(amu(typat(iat1))*amu_emass)
     mass_phongreen(irow,icol)= -mfac1*work1(irow,icol)*mfac2
   end do
 end do

!Finally extract the spectral function from the trace
 do irow= 1, pdim
   mode_phonspec(irow)= -one/pi * aimag(two*cmplx(omega,eta,16)*mass_phongreen(irow,irow))
 end do
 phonspec= sum(mode_phonspec(:))

!Diagonalize the Dynamical matrix
 ABI_MALLOC(matrx,(2,(3*natom*(3*natom+1))/2))
 ABI_MALLOC(eigval,(pdim))
 do icol= 1, pdim
   do irow= 1, icol
     matrx(1,irow + (icol-1)*icol/2)=real(dynmat(irow,icol))
     matrx(2,irow + (icol-1)*icol/2)=aimag(dynmat(irow,icol))
   end do
 end do 

 ABI_MALLOC(zhpev1,(2,2*3*natom-1))
 ABI_MALLOC(zhpev2,(3*3*natom-2))

 call ZHPEV ('V','U',3*natom,matrx,eigval,eigvec,3*natom,zhpev1,zhpev2,ier)
 ABI_CHECK(ier == 0, sjoin('zhpev returned:', itoa(ier)))

 ABI_FREE(matrx)
 ABI_FREE(zhpev1)
 ABI_FREE(zhpev2)

!Get the phonon frequencies (negative by convention, if the eigenvalue of the dynamical matrix is negative)
 do imode=1,3*natom
   if(eigval(imode)>=1.0d-16)then
     phfrq(imode)=sqrt(eigval(imode))
   else if(eigval(imode)>=-1.0d-16)then
     phfrq(imode)=zero
   else
     phfrq(imode)=-sqrt(-eigval(imode))
   end if
 end do

!Fix the phase of the eigenvectors
 call fxphas_seq(eigvec,dum, 0, 0, 1, 3*natom*3*natom, 0, 3*natom, 3*natom, 0)

!Normalise the eigenvectors
 call pheigvec_normalize(natom, eigvec)

 ! Get the phonon displacements
 call phdispl_from_eigvec(natom, ntypat, typat, amu, eigvec, displ)

 ABI_FREE(ifc)
 ABI_FREE(ifc_w2mass)
 ABI_FREE(ipiv)
 ABI_FREE(work1)
 ABI_FREE(mass_phongreen)
 ABI_FREE(invmassfac)
! ABI_FREE(delta_asrw0)

 DBG_EXIT("COLL")

end subroutine phonon_green
!!***

!!****f* ABINIT/lm_normal_modes
!! NAME
!!  lm_normal_modes
!!
!! FUNCTION
!!  Calculates the lattice-mediated contributions of the different type of
!!  susceptibilities by projecting the calculation on the contributions of 
!!  the phonon modes. 
!!
!! COPYRIGHT
!!  Copyright (C) 2024 ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amu(ntypat)= atomic masses
!!  blkval(2,3,mpert,3,mpert,1)= array with second-order derivatives
!!  displ(2,3,natom,3,natom)= phonon eigendisplacements
!!  mmom(ndim,(natom+2)*3)= first-order magnetic moments
!!  natom= number of atoms in the cell
!!  ndim= dimension of the penalized degrees of freedom
!!  ntypat= number of atom types in the cell
!!  typat(natom)= array with the type of atoms in the cell
!!
!! OUTPUT
!!  modemm(ndim,3*natom)= mode-resolved local magnetic moments
!!  modezf(ndim,3*natom)= mode-resolved local Zeeman fields
!!  modemeff(3,3*natom)= mode-resolved magnetic Born charges
!!  modezeff(3,3*natom)= mode-resolved electric Born charges
!!
!! SIDE EFFECTS
!!
!! NOTES
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

subroutine lm_normal_modes(blkval,displ,eta,lm_alpha_nm,lm_epsilon_nm,lm_mchi_nm, &
& mmom,modemm,modedisp,modemeff,modezeff,modezf,mpert,natom,ndim,omega,phfrq,ucvol,zfield)

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: mpert,natom,ndim 
 real(dp), intent(in) :: eta,omega,ucvol
!arrays
 real(dp), intent(in) :: blkval(2,3,mpert,3,mpert,1)
 real(dp), intent(in) :: displ(2*3*natom*3*natom)
 real(dp), intent(in) :: phfrq(3*natom)
 complex(dpc), intent(in) :: mmom(ndim,(natom+5)*3)
 complex(dpc), intent(in) :: zfield(ndim,(natom+5)*3)
 complex(dpc), intent(out) :: lm_alpha_nm(3,3,3*natom)
 complex(dpc), intent(out) :: lm_epsilon_nm(3,3,3*natom)
 complex(dpc), intent(out) :: lm_mchi_nm(3,3,3*natom)
 complex(dpc), intent(out) :: modemm(ndim,3*natom)
 complex(dpc), intent(out) :: modedisp(3*natom,3*natom)
 complex(dpc), intent(out) :: modemeff(3,3*natom)
 complex(dpc), intent(out) :: modezeff(3,3*natom)
 complex(dpc), intent(out) :: modezf(ndim,3*natom)

!Local variables-------------------------------
!scalars
 integer :: i1,iat1,iat2,idir1,idir2,icol,im,imode,index,ipert,ipert1,ipert2,irow,jmode,jpert
 integer :: pdim
 real(dp) :: fac
 complex(dpc) :: cplx_eta,cplx_w2
!arrays
 real(dp), allocatable :: mass(:)
 complex(dpc),allocatable :: c_blkval(:,:,:,:),norm(:)
 complex(dpc),allocatable :: zeff(:,:), zeff_tr(:,:), modezeff_tr(:,:)
 complex(dpc),allocatable :: meff(:,:), meff_tr(:,:), modemeff_tr(:,:)
!character(len=500) :: msg                   

! *************************************************************************

 DBG_ENTER("COLL")

!Define the complex eigendisplacementes array
do imode=1,3*natom
   do idir1=1,3
     do ipert1=1,natom
       i1=idir1+(ipert1-1)*3
       index=i1+3*natom*(imode-1)
       modedisp(i1,imode)= cmplx(displ(2*index-1),displ(2*index),16)
     end do
   end do
 end do

 ABI_MALLOC(norm,(3*natom))
 do imode= 1, 3*natom
   norm(imode)= sqrt(dot_product(modedisp(:,imode),modedisp(:,imode)))
 end do 

!Compute the mode-resolved macroscopic quantities
!(Born and magnetic charges)
 pdim= 3*natom
 ABI_MALLOC(c_blkval,(3,mpert,3,mpert))
 c_blkval= cmplx(blkval(1,:,:,:,:,1),blkval(2,:,:,:,:,1),16)

!Born charges
 ABI_MALLOC(zeff,(3,pdim))
 ABI_MALLOC(zeff_tr,(pdim,3))
 ABI_MALLOC(modezeff_tr,(pdim,3))
 ABI_MALLOC(meff,(3,pdim))
 ABI_MALLOC(meff_tr,(pdim,3))
 ABI_MALLOC(modemeff_tr,(pdim,3))
 modezeff(:,:)=(zero,zero)
 modemeff(:,:)=(zero,zero)
 modezeff_tr(:,:)=(zero,zero)
 modemeff_tr(:,:)=(zero,zero)
 ipert= natom + 2
 jpert= natom + 5
 do im= 1, 3
   do iat2= 1, natom
     do idir2= 1, 3
       imode= (iat2-1)*3 + idir2
       do iat1= 1, natom
         do idir1= 1, 3
           irow= (iat1-1)*3 + idir1

           !Electric Born charges
           zeff(im,irow)= c_blkval(im,ipert,idir1,iat1)
           zeff_tr(irow,im)= c_blkval(idir1,iat1,im,ipert)
           modezeff(im,imode)= modezeff(im,imode) + zeff(im,irow)* &
         & modedisp(irow,imode)
           modezeff_tr(imode,im)= modezeff_tr(imode,im) + zeff_tr(irow,im)* &
         & conjg(modedisp(irow,imode))

           !Magnetic Born charges
           meff(im,irow)= c_blkval(im,jpert,idir1,iat1)
           meff_tr(irow,im)= c_blkval(idir1,iat1,im,jpert)
           modemeff(im,imode)= modemeff(im,imode) + meff(im,irow)* &
         & modedisp(irow,imode)
           modemeff_tr(imode,im)= modemeff_tr(imode,im) + meff_tr(irow,im)* &
         & conjg(modedisp(irow,imode))

         end do
       end do
     end do
   end do
 end do
 ABI_FREE(zeff)
 ABI_FREE(zeff_tr)
 ABI_FREE(meff)
 ABI_FREE(meff_tr)

!Compute the mode-resolved local magnetic moments and fields
 modemm(:,:)=(zero,zero)
 modezf(:,:)=(zero,zero)
 do im= 1, ndim
   do iat2= 1, natom
     do idir2= 1, 3
       imode= (iat2-1)*3 + idir2
       do iat1= 1, natom
         do idir1= 1, 3
           irow= (iat1-1)*3 + idir1
           modemm(im,imode)= modemm(im,imode) +  mmom(im,irow)*modedisp(irow,imode)
           modezf(im,imode)= modezf(im,imode) +  zfield(im,irow)*modedisp(irow,imode)
         end do
       end do
     end do
   end do
 end do

!Compute the normal modes contribution to the susceptibilities
 cplx_eta= cmplx(0.0_dp,eta)
 cplx_w2= (omega+cplx_eta)**2

!Dielectric tensor
 fac= -four_pi/ucvol
 do idir1= 1, 3
   do idir2= 1, 3
     do imode= 1, pdim
       lm_epsilon_nm(idir1,idir2,imode)= fac*modezeff(idir1,imode)*modezeff_tr(imode,idir2)/ &
     & (cplx_w2 - phfrq(imode)**2)
     end do
   end do
 end do

!Magnetoelectric susceptibility
 fac= -one/ucvol
 do idir1= 1, 3
   do idir2= 1, 3
     do imode= 1, pdim
       lm_alpha_nm(idir1,idir2,imode)= fac*modemeff(idir1,imode)*modezeff_tr(imode,idir2)/ &
     & (cplx_w2 - phfrq(imode)**2)
     end do
   end do
 end do

!Magnetic susceptibility
 fac= -one/ucvol
 do idir1= 1, 3
   do idir2= 1, 3
     do imode= 1, pdim
       lm_mchi_nm(idir1,idir2,imode)= fac*modemeff(idir1,imode)*modemeff_tr(imode,idir2)/ &
     & (cplx_w2 - phfrq(imode)**2)
     end do
   end do
 end do

!Normalize mode-projected quantities
 do imode= 1, natom*3
   modezeff(:,imode)= modezeff(:,imode) / norm(imode)
   modemeff(:,imode)= modemeff(:,imode) / norm(imode)
   modemm(:,imode)= modemm(:,imode) / norm(imode)
   modezf(:,imode)= modezf(:,imode) / norm(imode)
 end do

 ABI_FREE(norm)
 ABI_FREE(c_blkval)
 ABI_FREE(modezeff_tr)
 ABI_FREE(modemeff_tr)

 DBG_EXIT("COLL")

end subroutine lm_normal_modes
!!***

!!****f* ABINIT/mode_mmom
!! NAME
!!  mode_mmom
!!
!! FUNCTION
!!  Projects the magnetic moments and Zeeman fields on the eigenmodes 
!!  of the dynamical matrix calculated at each value of omega
!!
!! COPYRIGHT
!!  Copyright (C) 2024 ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amu(ntypat)= atomic masses
!!  eigvec(2,3,natom,3,natom)= dynamical matrix eigenvectors
!!  mmom(ndim,(natom+2)*3)= first-order magnetic moments
!!  natom= number of atoms in the cell
!!  ndim= dimension of the penalized degrees of freedom
!!  ntypat= number of atom types in the cell
!!  typat(natom)= array with the type of atoms in the cell
!!
!! OUTPUT
!!  modemm(ndim,3*natom)= mode-resolved magnetic moments
!!  modezf(ndim,3*natom)= mode-resolved Zeeman fields
!!
!! SIDE EFFECTS
!!
!! NOTES
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


subroutine mode_mmom(amu,eigvec,mmom,modemm,modedisp,modezf,natom,ndim,ntypat,typat,zfield)

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: natom,ndim,ntypat 
!arrays
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: amu(ntypat)
 real(dp), intent(in) :: eigvec(2,3,natom,3,natom)
 complex(dpc), intent(in) :: mmom(ndim,(natom+5)*3)
 complex(dpc), intent(in) :: zfield(ndim,(natom+5)*3)
 complex(dpc), intent(out) :: modemm(ndim,3*natom)
 complex(dpc), intent(out) :: modedisp(3*natom,3*natom)
 complex(dpc), intent(out) :: modezf(ndim,3*natom)

!Local variables-------------------------------
!scalars
 integer :: iat1,iat2,idir1,idir2,icol,im,imode,irow
 real(dp) :: mcell
!arrays
 real(dp), allocatable :: mass(:)
!character(len=500) :: msg                   

! *************************************************************************

 DBG_ENTER("COLL")

!Define the mass factors
 ABI_MALLOC(mass,(natom))
 mcell=zero
 do iat1= 1, natom
   mass(iat1)= amu(typat(iat1))
   mcell= mcell + amu(typat(iat1))
 end do
 mass(:)=sqrt(mcell/mass(:))

!Reshape the eigenvector array
 do iat2= 1, natom
   do idir2= 1, 3
     imode= (iat2-1)*3 + idir2
     do iat1= 1, natom
       do idir1= 1, 3
         irow= (iat1-1)*3 + idir1
         modedisp(irow,imode)= cmplx(eigvec(1,idir1,iat1,idir2,iat2),eigvec(2,idir1,iat1,idir2,iat2),16)
       end do
     end do
   end do
 end do

!Compute the mode-resolved moments
 modemm(:,:)=(zero,zero)
 modezf(:,:)=(zero,zero)
 do im= 1, ndim
   do iat2= 1, natom
     do idir2= 1, 3
       imode= (iat2-1)*3 + idir2
       do iat1= 1, natom
         do idir1= 1, 3
           irow= (iat1-1)*3 + idir1
           modemm(im,imode)= modemm(im,imode) +  mass(iat1)*mmom(im,irow)*modedisp(irow,imode)
           modezf(im,imode)= modezf(im,imode) +  mass(iat1)*zfield(im,irow)*modedisp(irow,imode)
         end do
       end do
     end do
   end do
 end do

 ABI_FREE(mass)


 DBG_EXIT("COLL")

end subroutine mode_mmom
!!***

!!****f* ABINIT/mode_zeff
!! NAME
!!  mode_zeff
!!
!! FUNCTION
!!  Projects the Born effective charges on the eigenmodes of the dynamical 
!!  matrix calculated at each value of omega
!!
!! COPYRIGHT
!!  Copyright (C) 2024 ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amu(ntypat)= atomic masses
!!  eigvec(2,3,natom,3,natom)= dynamical matrix eigenvectors
!!  zeff(ndim,(natom+2)*3)= first-order magnetic moments
!!  natom= number of atoms in the cell
!!  ntypat= number of atom types in the cell
!!  typat(natom)= array with the type of atoms in the cell
!!  zeff(3,3*natom)= atomic Born effective charges at a given omega
!!
!! OUTPUT
!!  modezeff(3,3*natom)= mode-resolved magnetic moments
!!
!! SIDE EFFECTS
!!
!! NOTES
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


subroutine mode_zeff(amu,eigvec,modezeff,natom,ntypat,typat,zeff)

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: natom,ntypat 
!arrays
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: amu(ntypat)
 real(dp), intent(in) :: eigvec(2,3,natom,3,natom)
 complex(dpc), intent(in) :: zeff(3,natom*3)
 complex(dpc), intent(out) :: modezeff(3,3*natom)

!Local variables-------------------------------
!scalars
 integer :: iat1,iat2,idir1,idir2,icol,im,imode,irow
 real(dp) :: mcell
!arrays
 real(dp), allocatable :: mass(:)
 complex(dpc), allocatable :: mode_zeffspec(:,:)
!character(len=500) :: msg                   

! *************************************************************************

 DBG_ENTER("COLL")

!Define the mass factors
 ABI_MALLOC(mass,(natom))
 mcell=zero
 do iat1= 1, natom
   mass(iat1)= amu(typat(iat1))
   mcell= mcell + amu(typat(iat1))
 end do
 mass(:)=sqrt(mcell/mass(:))

!Compute the mode-resolved Born charges
 modezeff(:,:)=(zero,zero)
 do im= 1, 3
   do iat2= 1, natom
     do idir2= 1, 3
       imode= (iat2-1)*3 + idir2
       do iat1= 1, natom
         do idir1= 1, 3
           irow= (iat1-1)*3 + idir1
           modezeff(im,imode)= modezeff(im,imode) +  mass(iat1)*zeff(im,irow)* &
         & cmplx(eigvec(1,idir1,iat1,idir2,iat2),eigvec(2,idir1,iat1,idir2,iat2),16)
         end do
       end do
     end do
   end do
 end do

 ABI_FREE(mass)

 DBG_EXIT("COLL")

end subroutine mode_zeff
!!***

!!****f* ABINIT/me_altcalc
!! NAME
!!  me_altcalc
!!
!! FUNCTION
!!  Calculates the relaxed-ion magnetic moments induced 
!!  by an electric field
!!
!! COPYRIGHT
!!  Copyright (C) 2024 ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amu(ntypat)= atomic masses
!!  eigvec(2,3,natom,3,natom)= dynamical matrix eigenvectors
!!  natom= number of atoms in the cell
!!  ntypat= number of atom types in the cell
!!  typat(natom)= array with the type of atoms in the cell
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
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


subroutine me_altcalc(amu,eigvec,lm_magsus,lm_zfield,phongreen_fm,magsus,natom,ndim,ntypat,omega,ri_mmom,typat,&
& fmzeff_tr,zfield)

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: natom,ndim,ntypat 
 real(dp), intent(in) :: omega
!arrays
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: amu(ntypat)
 real(dp), intent(in) :: eigvec(2,3,natom,3,natom)
 complex(dpc), intent(in) :: lm_magsus(ndim,ndim)
 complex(dpc), intent(out) :: lm_zfield(ndim,3)
 complex(dpc), intent(in) :: magsus(ndim,ndim)
 complex(dpc), intent(in) :: phongreen_fm(3*natom,3*natom)
 complex(dpc), intent(out) :: ri_mmom(ndim,3)
 complex(dpc), intent(in) :: fmzeff_tr(3*natom,3)
 complex(dpc), intent(in) :: zfield(ndim,(natom+2)*3)

!Local variables-------------------------------
!scalars
 integer :: i,iat1,iat2,idir1,idir2,icol,im,imode1,imode2,irow,j,k,l
 integer :: zfield_unit
 real(dp) :: mcell,mfac1,mfac2,norm
 complex*16,parameter :: ure=(1.d0,0.d0),uim=(0.d0,1.d0)
!arrays
 real(dp), allocatable :: mass(:)
 complex(dpc), allocatable :: me_altcalcspec(:,:)
 complex(dpc),allocatable :: mass_phongreen(:,:)
 complex(dpc) :: ri_zfield(ndim,3),ri_magsus(ndim,ndim)
 complex(dpc) :: eigdisp(3*natom,3*natom)
 complex(dpc) :: nm_zfield(ndim,3*natom)
 complex(dpc) :: nm_fmzeff_tr(3*natom,3)
 complex(dpc) :: nm_phongreen(3*natom,3*natom)
 complex(dpc) :: basein(2,2),baseout(2,2),rmat(2,2),hmat(2,2)
 complex(dpc) :: vecin(2),vecout(2),totbasein(3*natom,2),totbaseout(3*natom,2)

! *************************************************************************

 DBG_ENTER("COLL")

!Define the mass factors
 ABI_MALLOC(mass,(natom))
 mcell=zero
 do iat1= 1, natom
   mass(iat1)= amu(typat(iat1))
   mcell= mcell + amu(typat(iat1))
 end do
 mass(:)=sqrt(mcell/mass(:))

 ABI_MALLOC(mass_phongreen,(3*natom,3*natom))
 do icol= 1, 3*natom
   iat2= ceiling(icol/three)
   mfac2= sqrt(amu(typat(iat2))*amu_emass)
   do irow= 1, 3*natom
     iat1= ceiling(irow/three)
     mfac1= sqrt(amu(typat(iat1))*amu_emass)
     mass_phongreen(irow,icol)= mfac1*phongreen_fm(irow,icol)*mfac2
   end do
 end do

!Different formula on the sublattice space
 !Spin part
 ri_magsus= magsus + lm_magsus
 ri_mmom=-matmul(ri_magsus,zfield(:,(natom+1)*3+1:(natom+2)*3))
 !Lattice part
 lm_zfield= matmul(zfield(:,1:3*natom),matmul(phongreen_fm,fmzeff_tr))
! ri_mmom= ri_mmom + matmul(ri_magsus,lm_zfield)

!Project on the normal-modes space the lattice part
 do iat2= 1, natom
   do idir2= 1, 3
     imode2= (iat2-1)*3 + idir2
     do iat1= 1, natom
       do idir1= 1, 3
         imode1= (iat1-1)*3 + idir1
         eigdisp(imode1,imode2)= &
       & cmplx(eigvec(1,idir1,iat1,idir2,iat2),eigvec(2,idir1,iat1,idir2,iat2),16)
       end do
     end do
   end do
 end do

!tmp change the phase of mode 21
! do i=1, natom*3
!   eigdisp(i,21)= (0.d0,1.d0)*eigdisp(i,21)
! end do

!enforce a circularly polarized pair of modes 20 and 21
!B
 baseout(:,1)=1.d0/sqrt(2.d0)*(/ure,uim/)
 baseout(:,2)=1.d0/sqrt(2.d0)*(/ure,-uim/)

!A
! norm=dot_product(eigdisp(1:2,20),eigdisp(1:2,20))
! basein(:,1)=1.d0/sqrt(norm)*eigdisp(1:2,20)    

! norm=dot_product(eigdisp(1:2,21),eigdisp(1:2,21))
! basein(:,2)=1.d0/sqrt(norm)*eigdisp(1:2,21) 

  basein(:,1)=eigdisp(1:2,20)
  basein(:,2)=eigdisp(1:2,21)
  hmat=matmul(transpose(conjg(basein)),basein)

  basein=basein/sqrt(hmat(1,1))

! rmat=matmul(baseout,transpose(conjg(basein)))
 rmat=matmul(transpose(conjg(baseout)),basein)
!  rmat=matmul(basein,transpose(conjg(baseout)))

 totbasein(:,1)=eigdisp(:,20)
 totbasein(:,2)=eigdisp(:,21)
 totbaseout=matmul(totbasein,transpose(conjg(rmat)))
 eigdisp(:,20)=totbaseout(:,1)
 eigdisp(:,21)=totbaseout(:,2)

! do i= 1, 3*natom, 3
!   vecin(1:2)= eigdisp(i:i+1,20)
!   vecout(:)=matmul(transpose(conjg(rmat)),vecin)
!   eigdisp(i:i+1,20)= vecout(1:2)
!
!   vecin(1:2)= eigdisp(i:i+1,21)
!   vecout(:)= matmul(vecin,transpose(conjg(rmat)))
!   vecout(:)=matmul(transpose(conjg(rmat)),vecin)
!   eigdisp(i:i+1,21)= vecout(1:2)
! end do

!now do the calculation
 do i=1, ndim
!   nm_zfield(i,1:3*natom)=matmul(conjg(zfield(i,1:3*natom)),eigdisp)
    nm_zfield(i,1:3*natom)=matmul(zfield(i,1:3*natom),eigdisp)
 end do
 nm_fmzeff_tr=matmul(transpose(conjg(eigdisp)),fmzeff_tr)
 nm_phongreen=matmul(transpose(conjg(eigdisp)),matmul(phongreen_fm,eigdisp))
 
 lm_zfield= matmul(nm_zfield(:,1:3*natom),matmul(nm_phongreen,nm_fmzeff_tr))

 !restrict only to a set of modes
! lm_zfield=cmplx(zero,zero,16)
! do i= 1, ndim
!   do j= 1, 3
!     do k= 21,21
!       do l= 21,21
!         lm_zfield(i,j)= lm_zfield(i,j) + nm_zfield(i,k)*nm_phongreen(k,l)*nm_fmzeff_tr(l,j)
!       end do 
!     end do 
!   end do
! end do

 write(120,*) omega,real(nm_zfield(1:4,20))
 write(220,*) omega,aimag(nm_zfield(1:4,20))
 write(121,*) omega,real(nm_zfield(1:4,21))
 write(221,*) omega,aimag(nm_zfield(1:4,21))

 write(420,*) omega,real(nm_fmzeff_tr(20,1:3))
 write(421,*) omega,real(nm_fmzeff_tr(21,1:3))
 write(520,*) omega,aimag(nm_fmzeff_tr(20,1:3))
 write(521,*) omega,aimag(nm_fmzeff_tr(21,1:3))


 do i=1, 3*natom
   write(320,*) real(eigdisp(i,20)),aimag(eigdisp(i,20))
   write(321,*) real(eigdisp(i,21)),aimag(eigdisp(i,21))
 end do 
 write(320,*) 
 write(321,*) 

 ri_mmom= ri_mmom + matmul(ri_magsus,lm_zfield)

 ABI_FREE(mass)
 ABI_FREE(mass_phongreen)

 DBG_EXIT("COLL")

end subroutine me_altcalc
!!***

!!****f* ABINIT/ri_d2etot
!! NAME
!!  ri_d2etot
!!
!! FUNCTION
!!  Extracts the CI 2nd-order susceptibilities and calculates the corresponding 
!!  lattice-mediated contributions
!!
!! COPYRIGHT
!!  Copyright (C) 2024 ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  blkval(2,3,mpert,3,mpert,1)= array with second-order derivatives
!!  magsus(ndim,ndim)= local magnetic susceptibility (RS) or its inverse (FS)
!!  mpert= maximum number of perturbations
!!  mcoup(ndim,natom+5)= magnetic Zeman fields (FS) or moments (RS)
!!  mcoup_tr(natom+5,ndim)= hermitian conjugate of mmom 
!!  natom= number of atoms in the cell
!!  ndim= dimension of the local magnetic degrees of freedom
!!  phongreen(3*natom,3*natom)= Phonon Green's function
!!
!! OUTPUT
!!  blkval_lm(2,3,mpert,3,mpert,1)= array with the lattice-mediated 
!!    second-order derivatives
!!
!! SIDE EFFECTS
!!
!! NOTES
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


 subroutine ri_d2etot(blkval,ci_alpha,ci_epsilon,ci_localpha,ci_mchi,&
& lm_alpha,lm_epsilon,lm_localpha,lm_magsus,lm_mchi,magsus,mpert,mcoup, &
& mcoup_tr,natom,ndim,phongreen,ucvol,zeff)

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: mpert,natom,ndim
 real(dp), intent(in) :: ucvol
!arrays
 real(dp), intent(in) :: blkval(2,3,mpert,3,mpert,1)
 complex(dpc), intent(out) :: ci_alpha(3,3)
 complex(dpc), intent(out) :: lm_alpha(3,3)
 complex(dpc), intent(out) :: ci_localpha(ndim,3)
 complex(dpc), intent(out) :: lm_localpha(ndim,3)
 complex(dpc), intent(out) :: ci_epsilon(3,3)
 complex(dpc), intent(out) :: lm_epsilon(3,3)
 complex(dpc), intent(out) :: ci_mchi(3,3)
 complex(dpc), intent(out) :: lm_mchi(3,3)
 complex(dpc), intent(inout) :: magsus(ndim,ndim)
 complex(dpc), intent(out) :: lm_magsus(ndim,ndim)
 complex(dpc), intent(in) :: mcoup(ndim,(natom+5)*3)
 complex(dpc), intent(in) :: mcoup_tr((natom+5)*3,ndim)
 complex(dpc), intent(in) :: phongreen(3*natom,3*natom)
 complex(dpc), intent(out) :: zeff(3,3*natom)

!Local variables-------------------------------
!scalars
 integer :: i,iat1,iat2,icol,idir1,idir2,ipert1,ipert2,irow,j,k,l
 integer :: pdim
 real(dp) :: fac
!arrays
 complex(dpc),allocatable :: c_blkval(:,:,:,:)
 complex(dpc),allocatable :: coup(:,:),coup_tr(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 pdim= 3*natom
 ABI_MALLOC(c_blkval,(3,mpert,3,mpert))
 c_blkval= cmplx(blkval(1,:,:,:,:,1),blkval(2,:,:,:,:,1),16)

 !Dielectric tensor
 fac= -four_pi/ucvol
 ABI_MALLOC(coup,(3,pdim))
 ABI_MALLOC(coup_tr,(pdim,3))
 ipert1= natom + 2
 do idir1= 1, 3
   irow= idir1
   do ipert2= 1, natom
     do idir2= 1, 3
       icol= (ipert2-1)*3 + idir2
       coup(irow,icol)= c_blkval(idir1,ipert1,idir2,ipert2)
       coup_tr(icol,irow)= c_blkval(idir2,ipert2,idir1,ipert1)
       zeff(irow,icol)= coup(irow,icol)
     end do
   end do
 end do
 ipert2= natom + 2
 lm_epsilon= fac*matmul(coup(:,:),matmul(phongreen,coup_tr(:,:)))
 
 do idir1= 1, 3
   do idir2= 1, 3
     ci_epsilon(idir1,idir2)= c_blkval(idir1,ipert1,idir2,ipert2)
   end do
 end do 

 !Magnetoelectric susceptibility
 fac= -one/ucvol
 ipert1= natom + 5
 do idir1= 1, 3
   irow= idir1
   do ipert2= 1, natom
     do idir2= 1, 3
       icol= (ipert2-1)*3 + idir2
       coup(irow,icol)= c_blkval(idir1,ipert1,idir2,ipert2)
     end do
   end do
 end do

 ipert1= natom + 2
 do idir1= 1, 3
   irow= idir1
   do ipert2= 1, natom
     do idir2= 1, 3
       icol= (ipert2-1)*3 + idir2
       coup_tr(icol,irow)= c_blkval(idir2,ipert2,idir1,ipert1)
     end do
   end do
 end do

 lm_alpha= fac*matmul(coup(:,:),matmul(phongreen,coup_tr(:,:)))
 
 ipert1= natom + 5
 ipert2= natom + 2
 do idir1= 1, 3
   do idir2= 1, 3
     ci_alpha(idir1,idir2)= c_blkval(idir1,ipert1,idir2,ipert2)/ucvol
   end do
 end do 

 !Local magnetoelectric susceptibilty
  ci_localpha(:,:)= mcoup(:,(natom+1)*3+1:(natom+2)*3)
  lm_localpha(:,:)= -matmul(mcoup(:,1:natom*3),matmul(phongreen,coup_tr(:,:))) 

 !Magnetic susceptibility 
 fac= -one/ucvol
 ipert1= natom + 5
 do idir1= 1, 3
   irow= idir1
   do ipert2= 1, natom
     do idir2= 1, 3
       icol= (ipert2-1)*3 + idir2
       coup(irow,icol)= c_blkval(idir1,ipert1,idir2,ipert2)
       coup_tr(icol,irow)= c_blkval(idir2,ipert2,idir1,ipert1)
     end do
   end do
 end do
 ipert2= natom + 5
 lm_mchi= fac*matmul(coup(:,:),matmul(phongreen,coup_tr(:,:)))
 
 do idir1= 1, 3
   do idir2= 1, 3
     ci_mchi(idir1,idir2)= c_blkval(idir1,ipert1,idir2,ipert2)/ucvol
   end do
 end do 
 
 !Local magnetic susceptibility 
 lm_magsus(:,:)= -matmul(mcoup(:,1:natom*3),matmul(phongreen,mcoup_tr(1:natom*3,:)))

 ABI_FREE(c_blkval)
 ABI_FREE(coup)
 ABI_FREE(coup_tr)

 DBG_EXIT("COLL")

end subroutine ri_d2etot
!!***

end module m_ddb_omega_interpol
!!***
