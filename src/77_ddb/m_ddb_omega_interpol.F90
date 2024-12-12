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
 use m_ddb_magpen,      only : local_spinsus, magmom
 use m_fstrings,        only : itoa, sjoin
 use m_macroave,        only : POLINT
 use m_io_tools,        only : open_file
 use m_cgtools,         only : fxphas_seq
 use m_dynmat,          only : pheigvec_normalize
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
 integer :: diel_unit,fs2rs,i,iblok,ifound,ii,imode,ipert1,ipert2,iw,j,jblok,jw,kblok,lblok,mmom_unit,mmspec_unit,nblok,ndim 
 integer :: nmat,nmdir,nwcalc,phon_unit,prtopt,spin_unit,zeff_unit,zeffspec_unit,zfield_unit
 real(dp) :: omegastp
 character(len=5000) :: msg,pfmt
 character(len=fnlen) :: diel_filename,spin_filename,mmom_filename,mmspec_filename
 character(len=fnlen) :: phon_filename,zeff_filename,zeffspec_filename,zfield_filename
 complex(dpc) :: cplxvar,cplx_weta
!arrays
 real(dp) :: qphnrm(3),qphon(3,3)
 real(dp), allocatable :: dint_fsddb(:,:),int_fsddb(:,:,:),omega(:),omegacalc(:)
 real(dp), allocatable :: w0hessian(:,:),w0berry(:,:)
 real(dp), allocatable :: eigvec(:,:,:,:,:),eigvec_fm(:,:,:,:,:),phfrq(:,:)
 real(dp), allocatable :: magphonspec(:),mode_magphonspec(:,:),mode_phonspec(:,:),phonspec(:)
 real(dp), allocatable :: coeffs(:,:,:)
 complex(dpc), allocatable :: dummysus(:,:)
 complex(dpc), allocatable :: invmagsus(:,:,:), lm_magsus(:,:,:), magsus(:,:,:), invhmat(:,:)
 complex(dpc), allocatable :: dummymom(:,:),dummymom_tr(:,:),mmom(:,:,:), mmom_tr(:,:,:)
 complex(dpc), allocatable :: ri_mmom(:,:,:)
 complex(dpc), allocatable :: lm_zfield(:,:,:),zfield(:,:,:),zfield_tr(:,:)
 complex(dpc), allocatable :: bc_barmagsus(:,:),bc_ss(:,:),bc_sp(:,:)
 complex(dpc), allocatable :: barepsilon(:,:,:),epsilon(:,:,:),ifcmat(:,:),ifcmat_fm(:,:)
 complex(dpc), allocatable :: modemm(:,:,:),zeff(:,:),zeff_tr(:,:),modezeff(:,:,:)
 complex(dpc), allocatable :: fmzeff(:,:),fmzeff_tr(:,:)
 complex(dpc), allocatable :: zeffspec(:,:),mmomspec(:,:),magphongreen(:,:),phongreen(:,:),phongreen_fm(:,:)
 complex(dpc), allocatable :: lm_epsilon(:,:,:),ri_magelsus(:,:),lm_magelsus(:,:,:)
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
 ABI_MALLOC(omega,(nomega))
 ABI_MALLOC(phfrq,(3*natom,nomega))
 ABI_MALLOC(phonspec,(nomega))
 ABI_MALLOC(magphonspec,(nomega))
 ABI_MALLOC(magphongreen,(3*natom+ndim,3*natom+ndim))
 ABI_MALLOC(phongreen,(3*natom,3*natom))
 ABI_MALLOC(phongreen_fm,(3*natom,3*natom))
 ABI_MALLOC(mode_phonspec,(3*natom,nomega))
 ABI_MALLOC(mode_magphonspec,(3*natom+ndim,nomega))
 ABI_MALLOC(zeffspec,(3,nomega))
 ABI_MALLOC(mmomspec,(ndim,nomega))
 ABI_MALLOC(eigvec,(2,3,natom,3,natom))
 ABI_MALLOC(eigvec_fm,(2,3,natom,3,natom))
 ABI_MALLOC(modemm,(ndim,3*natom,nomega))
 ABI_MALLOC(fmzeff,(3,3*natom))
 ABI_MALLOC(fmzeff_tr,(3*natom,3))
 ABI_MALLOC(zeff,(3,3*natom))
 ABI_MALLOC(zeff_tr,(3*natom,3))
 ABI_MALLOC(genzeff_tr,(3*natom+ndim,3))
 ABI_MALLOC(ri_genelsus,(3*natom+ndim,3,nomega))
 ABI_MALLOC(lm_epsilon,(3,3,nomega))
 ABI_MALLOC(modezeff,(3,3*natom,nomega))
 ABI_MALLOC(dummysus,(ndim,ndim))
 ABI_MALLOC(magsus,(ndim,ndim,nomega))
 ABI_MALLOC(lm_magsus,(ndim,ndim,nomega))
 ABI_MALLOC(ri_magelsus,(ndim,3))
 ABI_MALLOC(lm_magelsus,(ndim,3,nomega))
 ABI_MALLOC(invmagsus,(ndim,ndim,nomega))
 ABI_MALLOC(dummymom,(ndim,(natom+2)*3))
 ABI_MALLOC(dummymom_tr,(ndim,(natom+2)*3))
 ABI_MALLOC(mmom,(ndim,(natom+2)*3,nomega))
 ABI_MALLOC(mmom_tr,((natom+2)*3,ndim,nomega))
 ABI_MALLOC(lm_zfield,(ndim,3,nomega))
 ABI_MALLOC(zfield,(ndim,(natom+5)*3,nomega))
 ABI_MALLOC(zfield_tr,((natom+5)*3,ndim))
 ABI_MALLOC(barepsilon,(3,3,nomega))
 ABI_MALLOC(epsilon,(3,3,nomega))
 ABI_MALLOC(macmagsus,(3,3,nomega))
 ABI_MALLOC(ifcmat,(3*natom,3*natom))
 ABI_MALLOC(ifcmat_fm,(3*natom,3*natom))
 ABI_MALLOC(dint_fsddb,(2,ddb%msize))
 ABI_MALLOC(int_fsddb,(2,ddb%msize,1))
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

   !Perform the different interpolations
   !Lineal (with analytic Berry curvature) with dissipation if eta/=0
   if (omegaflag==1) then
!     call lineal_omega_interp(w0hessian,w0berry,eta,ifcmat_fm, &
!   & invmagsus(:,:,iw),magsus(:,:,iw),mpatpol,mpdir,mpert,ddb%msize, &
!   & natom,ndim,nmdir,int_fsddb,omega(iw),zfield(:,:,iw),zfield_tr)
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

   !Calculate the magnetic moments
   call magmom(dummymom,dummymom_tr,ddb,dummysus,dummysus,1,1,1,magpen,&
   & magsus(:,:,iw),mmom(:,:,iw),mmom_tr(:,:,iw),mpatpol,mpdir,mpert,natom,&
   &1,ndim,nmdir,prtopt,prtvol,qphon,xred,zfield(:,:,iw),zfield_tr,&
   &fs2rs=fs2rs,blkval_fs=int_fsddb)
  
!     !Calculate the dielectric susceptibility
!     call mp_diel(barepsilon(:,:,iw),barmagsus(:,:,iw),barmmom,barmmom_tr,&
!   & int_fsddb,dissip,epsilon(:,:,iw),1,invhmat,magpen,magsus(:,:,iw),mpert,mpopt,&
!   & natom,1,ndim,prtopt,prtvol,ucvol,zfield(:,:,iw),zfield_tr)
!  
!     !Calculate the interatomic force constants
!     call mp_ifc(barmagsus(:,:,iw),barmmom,barmmom_tr,int_fsddb,dissip,1,ifcmat,&
!   & ifcmat_fm,invhmat,magsus(:,:,iw),magpen,mpert,mpopt,&
!   & natom,1,ndim,omega(iw),omegaflag,prtopt,prtvol,qphon,xred,zfield(:,:,iw),zfield_tr)
!  
!     !Calculate the macroscopic magnetic susceptibility
!     call mp_macmagsus(barmagsus(:,:,iw),int_fsddb,&
!   & dissip,1,invbarmagsus(:,:,iw),invhmat, macmagsus(:,:,iw),magpen,&
!   & magsus(:,:,iw),mpatpol,mpdir,mpert,mpopt,&
!   & natom,nblok,ndim,nmdir,prtopt,prtvol,ucvol)

   !Apply ASR
!   if (omega(1) < tol12) then
!     call asrw0(delta_asrw0,ifcmat,natom,0) 
!     call asrw0(delta_asrw0_fm,ifcmat_fm,natom,0) 
!   else 
!     call asrw0(delta_asrw0,ifcmat,natom,1) 
!     call asrw0(delta_asrw0_fm,ifcmat_fm,natom,1) 
!   end if

   !Calculate the phonon and magnon-phonon Green's functions and spectral functions
!   call phonon_green(amu,eigvec,eigvec_fm,eta,ifcmat,ifcmat_fm,invmagsus(:,:,iw),& 
! & magphongreen,magphonspec(iw),mode_magphonspec(:,iw),mode_phonspec(:,iw),natom,ndim,ntypat,omega(iw),&
! & phfrq(:,iw),phongreen,phongreen_fm,phonspec(iw),typat,zfield(:,:,iw),zfield_tr)

!   if (omegaflag==2.or.omegaflag==3) then
     !Calculate the mode-resolved magnetic moments
!     call mode_mmom(amu,eigvec,mmom(:,:,iw),mmomspec(:,iw),modemm(:,:,iw),mode_phonspec(:,iw),natom,ndim,ntypat,typat)

!     !Calculate the Born effective charges
!     call mp_zeff(barmagsus(:,:,iw),barmmom,barmmom_tr,int_fsddb,&
!   & dissip,fmzeff,fmzeff_tr,1,invhmat,lm_epsilon(:,:,iw),magpen,magsus(:,:,iw),mpert,mpopt,&
!   & natom,1,ndim,phongreen,prtopt,prtvol,ucvol,zeff,zeff_tr,zfield(:,:,iw),zfield_tr)

     !Calculate the mode-resolved Born effective charges
!     call mode_zeff(amu,eigvec,mode_phonspec(:,iw),modezeff(:,:,iw),natom,ntypat,&
!   & typat,zeff(:,:,iw),zeffspec(:,iw))

!     !Calculate here the phonons contribution to the spin susceptibility
!     if (dissip==0) then
!       lm_magsus(:,:,iw)=-matmul(mmom(:,1:natom*3,iw),matmul(phongreen,transpose(conjg(mmom(:,1:natom*3,iw)))))
!     else if (dissip==1) then
!       lm_magsus(:,:,iw)=-matmul(mmom(:,1:natom*3,iw),matmul(phongreen,mmom_tr(1:natom*3,:,iw)))
!     end if
!
!     !Calclate here the lattice-mediated magnetic moments induced by an electric field
!     if (dissip==0) then
!       lm_magelsus(:,:,iw)=-matmul(mmom(:,1:natom*3,iw),matmul(phongreen,transpose(conjg(zeff(:,:)))))
!     else if (dissip==1) then
!       lm_magelsus(:,:,iw)=-matmul(mmom(:,1:natom*3,iw),matmul(phongreen,zeff_tr(:,:)))
!     end if
!
!     !Alternative calculation
!     genzeff_tr(1:natom*3,:)= fmzeff_tr(:,:)
!     genzeff_tr(natom*3:natom*3+ndim,:)= zfield(:,(natom+2)*3-2:(natom+2)*3,iw)
!     ri_genelsus(:,:,iw)=-matmul(magphongreen,genzeff_tr)
!
!     !Another alternative
!
!!TMP: FM eigvecs
!     call me_altcalc(amu,eigvec_fm,lm_magsus(:,:,iw),lm_zfield(:,:,iw),phongreen_fm,magsus(:,:,iw),natom,ndim,ntypat,&
!   & omega(iw),ri_mmom(:,:,iw),typat,fmzeff_tr,zfield(:,:,iw))
!
!!TMP: SR eigvecs
!!     call me_altcalc(amu,eigvec,lm_magsus(:,:,iw),lm_zfield(:,:,iw),phongreen_fm,magsus(:,:,iw),natom,ndim,ntypat,&
!!   & omega(iw),ri_mmom(:,:,iw),typat,fmzeff_tr,zfield(:,:,iw))
!
!     !Convert magnetic susceptibilities to the magnon basis
!!     work(:,:)=magsus(:,:,iw)
!!     magsus(:,:,iw)=matmul(transpose(conjg(magbasis)),matmul(work,magbasis))
!!     work(:,:)=lm_magsus(:,:,iw)
!!     lm_magsus(:,:,iw)=matmul(transpose(conjg(magbasis)),matmul(work,magbasis))
!   end if
 end do

!!Calculate the norm of the phonon spectral function
! totnorm=zero
! do imode=1, 3*natom
!   write(msg,'(a,i3,a,es15.7)') 'Phonon mode: ', imode, & 
! & '. Norm of the spectral function: ', sum(mode_phonspec(imode,:))*omegastp
!   call wrtout([ab_out,std_out],msg,'COLL')
!   totnorm= totnorm + sum(mode_phonspec(imode,:))*omegastp
! end do
! write(msg,'(a,es15.7)') & 
! & ' Total norm of the spectral function: ', totnorm
!   call wrtout([ab_out,std_out],msg,'COLL')
!
!!Calculate the norm of the magnon phonon spectral function
! totnorm=zero
! do imode=1, 3*natom+ndim
!   write(msg,'(a,i3,a,es15.7)') 'Magnon-Phonon mode: ', imode, & 
! & '. Norm of the spectral function: ', sum(mode_magphonspec(imode,:))*omegastp
!   call wrtout([ab_out,std_out],msg,'COLL')
!   totnorm= totnorm + sum(mode_magphonspec(imode,:))*omegastp
! end do
! write(msg,'(a,es15.7)') & 
! & ' Total norm of the spectral function: ', totnorm
!   call wrtout([ab_out,std_out],msg,'COLL')

!!!  Print results of interpolation
!Spin susceptibilities
 spin_filename=trim(outfilename_radix)//"_SPINSUS"
 if (open_file(spin_filename, msg, newunit=spin_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(spin_unit,*) '#'
 write(spin_unit,*) '#  Spin susceptibilities calculated and interpolated by ANADDB'
 write(spin_unit,*) '#'
 write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' )  ndim**2

 write(spin_unit,*) '#  Real part of clamped-ion local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(magsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

 write(spin_unit,*) ' '
 write(spin_unit,*) '#  Imaginary part of clamped-ion local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(magsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

! write(spin_unit,*) ' '
! write(spin_unit,*) '#  Real part of lattice-mediated local spin susceptibility tensor (at. units)'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
! call wrtout(spin_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((real(lm_magsus(i,j,iw)),j=1,ndim),i=1,ndim)
!    call wrtout(spin_unit,msg,'COLL')
! end do
!
! write(spin_unit,*) ' '
! write(spin_unit,*) '#  Imaginary part of lattice-mediated local spin susceptibility tensor (at. units)'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
! call wrtout(spin_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((aimag(lm_magsus(i,j,iw)),j=1,ndim),i=1,ndim)
!    call wrtout(spin_unit,msg,'COLL')
! end do
!
! write(spin_unit,*) ' '
! write(spin_unit,*) '#  Real part of relaxed-ion local spin susceptibility tensor (at. units)'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
! call wrtout(spin_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((real(magsus(i,j,iw)+lm_magsus(i,j,iw)),j=1,ndim),i=1,ndim)
!    call wrtout(spin_unit,msg,'COLL')
! end do
!
! write(spin_unit,*) ' '
! write(spin_unit,*) '#  Imaginary part of relaxed-ion local spin susceptibility tensor (at. units)'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
! call wrtout(spin_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((aimag(magsus(i,j,iw)+lm_magsus(i,j,iw)),j=1,ndim),i=1,ndim)
!    call wrtout(spin_unit,msg,'COLL')
! end do

 write(spin_unit,*) ' '
 write(spin_unit,*) '#  Real part of the inverse of the clamped-ion local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X^{-1}_11     X^{-1}_12     ...     X^{-1}_21     X^{-1}_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(invmagsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

 write(spin_unit,*) ' '
 write(spin_unit,*) '#  Imaginary part of the inverse of the clamped-ion local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X^{-1}_11     X^{-1}_12     ...     X^{-1}_21     X^{-1}_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(invmagsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

! write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' ) 9 
!
! write(spin_unit,*) ' '
! write(spin_unit,*) '#  Real part of clamped-ion macroscopic spin susceptibility tensor (at. units)'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
! call wrtout(spin_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((real(macmagsus(i,j,iw)),j=1,3),i=1,3)
!    call wrtout(spin_unit,msg,'COLL')
! end do
!
! write(spin_unit,*) ' '
! write(spin_unit,*) '#  Imaginary part of clamped-ion macroscopic spin susceptibility tensor (at. units)'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
! call wrtout(spin_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((aimag(macmagsus(i,j,iw)),j=1,3),i=1,3)
!    call wrtout(spin_unit,msg,'COLL')
! end do

 close (spin_unit)

!Zfields
 write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' )  ndim*3
 zfield_filename=trim(outfilename_radix)//"_ZFIELDS"
 if (open_file(zfield_filename, msg, newunit=zfield_unit) /= 0) then
   ABI_ERROR(msg)
 end if

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
! write(zfield_unit,*) '#  Real part of lattice-mediated local Zeeman fields induced by electric field (at. units)'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     Z_11     Z_12      Z_13    ...     Z_21     Z_22     ...'
! call wrtout(zfield_unit,msg,'COLL')
! do iw=1,nomega
!   write(msg,pfmt) omega(iw), ((real(lm_zfield(i,j,iw)),j=1,3),i=1,ndim)
!   call wrtout(zfield_unit,msg,'COLL')
! end do
!
! write(zfield_unit,*) ' '
! write(zfield_unit,*) '#  Imag part of lattice-mediated local Zeeman fields induced by electric field (at. units)'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     Z_11     Z_12      Z_13    ...     Z_21     Z_22     ...'
! call wrtout(zfield_unit,msg,'COLL')
! do iw=1,nomega
!   write(msg,pfmt) omega(iw), ((aimag(lm_zfield(i,j,iw)),j=1,3),i=1,ndim)
!   call wrtout(zfield_unit,msg,'COLL')
! end do
! 
 close (zfield_unit)

!Magnetic moments
 mmom_filename=trim(outfilename_radix)//"_MAGMOM"
 if (open_file(mmom_filename, msg, newunit=mmom_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(mmom_unit,*) '#'
 write(mmom_unit,*) '#  Magnetic moments calculated and interpolated by ANADDB'
 write(mmom_unit,*) '#'

! write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  ndim
! do imode= 1, 3*natom
!   write(mmom_unit,*) ' '
!   write(mmom_unit,'(a,i3)') '#  Real part of magnetic moments (at. units) induced by phonon mode:', imode
!   write(msg,'(a,a,a)') ch10,&
! &           ' # At  hw     m_{mat_1,1}     m_{mat_1,2}     ...     m_{mat_2,1}     m_{mat_2,2}'
!   call wrtout(mmom_unit,msg,'COLL')
!   do iw=1,nomega
!     write(msg,pfmt) &
!   & omega(iw), (real(modemm(i,imode,iw)),i=1,ndim)
!     call wrtout(mmom_unit,msg,'COLL')
!   end do
!   write(mmom_unit,*) ' '
!   write(mmom_unit,'(a,i3)') '#  Imaginary part of magnetic moments (at. units) induced by phonon mode:', imode
!   write(msg,'(a,a,a)') ch10,&
! &           ' # At  hw     m_{mat_1,1}     m_{mat_1,2}     ...     m_{mat_2,1}     m_{mat_2,2}'
!   call wrtout(mmom_unit,msg,'COLL')
!   do iw=1,nomega
!     write(msg,pfmt) &
!   & omega(iw), (aimag(modemm(i,imode,iw)),i=1,ndim)
!     call wrtout(mmom_unit,msg,'COLL')
!   end do
! end do

 write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  ndim*3
 write(mmom_unit,*) ' '
 write(mmom_unit,*) '#  Real part of the clamped-ion local magnetoelectric tensor(at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(mmom(i,j,iw)),j=3*(natom+1)+1,3*(natom+2)),i=1,ndim)
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
 &  omega(iw), ((aimag(mmom(i,j,iw)),j=3*(natom+1)+1,3*(natom+2)),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

! write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  ndim*3
! write(mmom_unit,*) ' '
! write(mmom_unit,*) '#  Real part of the phonon-modes contribution to the local magnetoelectric tensor(at. units)'
! write(msg,'(a,a,a)') ch10,&
!&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
!&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
! call wrtout(mmom_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((real(lm_magelsus(i,j,iw)),j=1,3),i=1,ndim)
!    call wrtout(mmom_unit,msg,'COLL')
! end do
!
! write(mmom_unit,*) ' '
! write(mmom_unit,*) '#  Imaginary part of the phonon-modes contribution to the local magnetoelectric tensor(at. units)'
! write(msg,'(a,a,a)') ch10,&
!&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
!&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
! call wrtout(mmom_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((aimag(lm_magelsus(i,j,iw)),j=1,3),i=1,ndim)
!    call wrtout(mmom_unit,msg,'COLL')
! end do
!
! write(mmom_unit,*) ' '
! write(mmom_unit,*) '#  Real part of relaxed-ion local magnetoelectric tensor (at. units)'
! write(msg,'(a,a,a)') ch10,&
!&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
!&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
! call wrtout(mmom_unit,msg,'COLL')
! do iw=1,nomega
!    ri_magelsus(:,:)=mmom(:,3*(natom+1)+1:3*(natom+2),iw)+lm_magelsus(:,:,iw)
!    write(msg,pfmt) &
! &  omega(iw), ((real(ri_magelsus(i,j)),j=1,3),i=1,ndim)
!    call wrtout(mmom_unit,msg,'COLL')
! end do
!
! write(mmom_unit,*) ' '
! write(mmom_unit,*) '#  Imaginary part of relaxed-ion local magnetoelectric tensor (at. units)'
! write(msg,'(a,a,a)') ch10,&
!&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
!&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
! call wrtout(mmom_unit,msg,'COLL')
! do iw=1,nomega
!    ri_magelsus(:,:)=mmom(:,3*(natom+1)+1:3*(natom+2),iw)+lm_magelsus(:,:,iw)
!    write(msg,pfmt) &
! &  omega(iw), ((aimag(ri_magelsus(i,j)),j=1,3),i=1,ndim)
!    call wrtout(mmom_unit,msg,'COLL')
! end do
!
! write(mmom_unit,*) ' '
! write(mmom_unit,*) '# [Alternative calc] Real part of relaxed-ion local magnetoelectric tensor (at. units)'
! write(msg,'(a,a,a)') ch10,&
!&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
!&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
! call wrtout(mmom_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((real(ri_genelsus(3*natom+i,j,iw)),j=1,3),i=1,ndim)
!    call wrtout(mmom_unit,msg,'COLL')
! end do
!
! write(mmom_unit,*) ' '
! write(mmom_unit,*) '#  [Alternative calc] Imaginary part of relaxed-ion local magnetoelectric tensor (at. units)'
! write(msg,'(a,a,a)') ch10,&
!&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
!&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
! call wrtout(mmom_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((aimag(ri_genelsus(3*natom+i,j,iw)),j=1,3),i=1,ndim)
!    call wrtout(mmom_unit,msg,'COLL')
! end do
!
! write(mmom_unit,*) ' '
! write(mmom_unit,*) '# [Another alternative calc] Real part of relaxed-ion local magnetoelectric tensor (at. units)'
! write(msg,'(a,a,a)') ch10,&
!&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
!&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
! call wrtout(mmom_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((real(ri_mmom(i,j,iw)),j=1,3),i=1,ndim)
!    call wrtout(mmom_unit,msg,'COLL')
! end do
!
! write(mmom_unit,*) ' '
! write(mmom_unit,*) '#  [Another alternative calc] Imaginary part of relaxed-ion local magnetoelectric tensor (at. units)'
! write(msg,'(a,a,a)') ch10,&
!&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
!&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
! call wrtout(mmom_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((aimag(ri_mmom(i,j,iw)),j=1,3),i=1,ndim)
!    call wrtout(mmom_unit,msg,'COLL')
! end do
!
! close(mmom_unit)
!
!! mmspec_filename=trim(outfilename_radix)//"_SPECTRAL_MAGMOM"
!! if (open_file(mmspec_filename, msg, newunit=mmspec_unit) /= 0) then
!!   ABI_ERROR(msg)
!! end if
!!
!! write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  ndim
!! write(mmspec_unit,*) ' '
!! write(mmspec_unit,'(a)') '#  Real part of magnetic moments spectral function:'
!! write(mmspec_unit,*) ' '
!! write(msg,'(a,a)') ch10,&
!! &           ' # At  hw     m_{mat_1,1}     m_{mat_1,2}     ...     m_{mat_2,1}     m_{mat_2,2}'
!! call wrtout(mmspec_unit,msg,'COLL')
!! do iw=1,nomega
!!   write(msg,pfmt) omega(iw), real(mmomspec(:,iw))
!!   call wrtout(mmspec_unit,msg,'COLL')
!! end do
!!
!! write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  ndim
!! write(mmspec_unit,*) ' '
!! write(mmspec_unit,'(a)') '#  Imaginary part of magnetic moments spectral function:'
!! write(mmspec_unit,*) ' '
!! write(msg,'(a,a)') ch10,&
!! &           ' # At  hw     m_{mat_1,1}     m_{mat_1,2}     ...     m_{mat_2,1}     m_{mat_2,2}'
!! call wrtout(mmspec_unit,msg,'COLL')
!! do iw=1,nomega
!!   write(msg,pfmt) omega(iw), aimag(mmomspec(:,iw))
!!   call wrtout(mmspec_unit,msg,'COLL')
!! end do
!!
!! close(mmspec_unit)
!
!!Dielectric susceptibility
! diel_filename=trim(outfilename_radix)//"_DIELSUS"
! if (open_file(diel_filename, msg, newunit=diel_unit) /= 0) then
!   ABI_ERROR(msg)
! end if
!
! write(diel_unit,*) '#'
! write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' ) 9
!
! write(diel_unit,*) '#  Real part of clamped-ion penalized dielectric tensor'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
! call wrtout(diel_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((real(barepsilon(i,j,iw)),j=1,3),i=1,3)
!    call wrtout(diel_unit,msg,'COLL')
! end do
!
! write(diel_unit,*) ' '
! write(diel_unit,*) '#  Imaginary part of clamped-ion penalized dielectric tensor'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
! call wrtout(diel_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((aimag(barepsilon(i,j,iw)),j=1,3),i=1,3)
!    call wrtout(diel_unit,msg,'COLL')
! end do
!
! write(diel_unit,*) '#'
! if (mpopt==1) then
!   write(diel_unit,*) '#  Frozen-magnetic clamped-ion dielectric tensor calculated and interpolated by ANADDB'
! else if (mpopt==2) then
!   write(diel_unit,*) '#  Spin-relaxed clamped-ion dielectric tensor calculated and interpolated by ANADDB'
! else
!   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
!   ABI_ERROR(msg)
! end if
! 
! write(diel_unit,*) '#'
! write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' ) 9 
!
! write(diel_unit,*) '#  Real part of clamped-ion dielectric tensor'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
! call wrtout(diel_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((real(epsilon(i,j,iw)),j=1,3),i=1,3)
!    call wrtout(diel_unit,msg,'COLL')
! end do
!
! write(diel_unit,*) ' '
! write(diel_unit,*) '#  Imaginary part of clamped-ion dielectric tensor'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
! call wrtout(diel_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((aimag(epsilon(i,j,iw)),j=1,3),i=1,3)
!    call wrtout(diel_unit,msg,'COLL')
! end do
!
! write(diel_unit,*) '#'
! if (mpopt==1) then
!   write(diel_unit,*) '#  Frozen-magnetic lattice-mediated dielectric tensor calculated and interpolated by ANADDB'
! else if (mpopt==2) then
!   write(diel_unit,*) '#  Spin-relaxed lattice-mediated dielectric tensor calculated and interpolated by ANADDB'
! else
!   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
!   ABI_ERROR(msg)
! end if
! 
! write(diel_unit,*) '#'
! write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' ) 9 
!
! write(diel_unit,*) '#  Real part of lattice-mediated dielectric tensor'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
! call wrtout(diel_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((real(lm_epsilon(i,j,iw)),j=1,3),i=1,3)
!    call wrtout(diel_unit,msg,'COLL')
! end do
!
! write(diel_unit,*) ' '
! write(diel_unit,*) '#  Imaginary part of lattice-mediated dielectric tensor'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
! call wrtout(diel_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((aimag(lm_epsilon(i,j,iw)),j=1,3),i=1,3)
!    call wrtout(diel_unit,msg,'COLL')
! end do
!
! write(diel_unit,*) ' '
! write(diel_unit,*) '#  Real part of relaxed-ion dielectric tensor'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
! call wrtout(diel_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((real(lm_epsilon(i,j,iw)+epsilon(i,j,iw)),j=1,3),i=1,3)
!    call wrtout(diel_unit,msg,'COLL')
! end do
!
! write(diel_unit,*) ' '
! write(diel_unit,*) '#  Imaginary part of relaxed-ion dielectric tensor'
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
! call wrtout(diel_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) &
! &  omega(iw), ((aimag(lm_epsilon(i,j,iw)+epsilon(i,j,iw)),j=1,3),i=1,3)
!    call wrtout(diel_unit,msg,'COLL')
! end do
!
! close(diel_unit)
!
!!Phonon spectral function
! phon_filename=trim(outfilename_radix)//"_SPECTRAL_PHONON"
! if (open_file(phon_filename, msg, newunit=phon_unit) /= 0) then
!   ABI_ERROR(msg)
! end if
!
! write(phon_unit,*) '#'
! if (mpopt==1) then
!   write(phon_unit,*) '#  Frozen-magnetic phonon spectral function calculated and interpolated by ANADDB'
! else if (mpopt==2) then
!   write(phon_unit,*) '#  Spin-relaxed phonon spectral function calculated and interpolated by ANADDB'
! else
!   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
!   ABI_ERROR(msg)
! end if
! 
! write(msg,'(a,a)') ch10, ' # At  hw               Phonon SF           Magnon-phonon SF'
! call wrtout(phon_unit,msg,'COLL')
! do iw=1,nomega
!   write(msg,*) omega(iw), phonspec(iw), magphonspec(iw)
!   call wrtout(phon_unit,msg,'COLL')
! end do
!
! close(phon_unit)
!
!!Phonon frequencies
! phon_filename=trim(outfilename_radix)//"_PHFRQ"
! if (open_file(phon_filename, msg, newunit=phon_unit) /= 0) then
!   ABI_ERROR(msg)
! end if
!
! write(phon_unit,*) '#'
! if (mpopt==1) then
!   write(phon_unit,*) '#  Frozen-magnetic phonon frequencies calculated and interpolated by ANADDB'
! else if (mpopt==2) then
!   write(phon_unit,*) '#  Spin-relaxed phonon frequencies calculated and interpolated by ANADDB'
! else
!   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
!   ABI_ERROR(msg)
! end if
!
! write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' ) natom*3
! write(msg,'(a,a)') ch10,&
!&           ' # At  hw    eval(1)     eval(2) ...'
! call wrtout(phon_unit,msg,'COLL')
! do iw=1,nomega
!    write(msg,pfmt) omega(iw), phfrq(:,iw)
!    call wrtout(phon_unit,msg,'COLL')
! end do
! 
! close(phon_unit)
!
!!!Born effective charges
!! zeff_filename=trim(outfilename_radix)//"_ZEFF"
!! if (open_file(zeff_filename, msg, newunit=zeff_unit) /= 0) then
!!   ABI_ERROR(msg)
!! end if
!!
!! write(zeff_unit,*) '#'
!! write(zeff_unit,*) '#  Born effective charges calculated and interpolated by ANADDB'
!! write(zeff_unit,*) '#'
!!
!! write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' ) 3 
!! do imode= 1, 3*natom
!!   write(zeff_unit,*) ' '
!!   write(zeff_unit,'(a,i3)') '#  Real part of Born charge (at. units) induced by phonon mode:', imode
!!   write(msg,'(a,a)') ch10,&
!! &           ' # At  hw     Z^x_{n}     Z^y_{n}     Z^z_{n}'
!!   call wrtout(zeff_unit,msg,'COLL')
!!   do iw=1,nomega
!!     write(msg,pfmt) &
!!   & omega(iw), (real(modezeff(i,imode,iw)),i=1,3)
!!     call wrtout(zeff_unit,msg,'COLL')
!!   end do
!!   write(zeff_unit,*) ' '
!!   write(zeff_unit,'(a,i3)') '#  Imaginary part of Born charge (at. units) induced by phonon mode:', imode
!!   write(msg,'(a,a)') ch10,&
!! &           ' # At  hw     Z^x_{n}     Z^y_{n}     Z^z_{n}'
!!   call wrtout(zeff_unit,msg,'COLL')
!!   do iw=1,nomega
!!     write(msg,pfmt) &
!!   & omega(iw), (aimag(modezeff(i,imode,iw)),i=1,3)
!!     call wrtout(zeff_unit,msg,'COLL')
!!   end do
!! end do
!!
!! close(zeff_unit)
!
!! zeffspec_filename=trim(outfilename_radix)//"_SPECTRAL_ZEFF"
!! if (open_file(zeffspec_filename, msg, newunit=zeffspec_unit) /= 0) then
!!   ABI_ERROR(msg)
!! end if
!! write(zeffspec_unit,*) ' '
!! write(zeffspec_unit,'(a)') '#  Real part of Born charges spectral function:'
!! write(zeffspec_unit,*) ' '
!! write(msg,'(a,a)') ch10,&
!! &           ' # At  hw     Z^x     Z^y     Z^z'
!! call wrtout(zeffspec_unit,msg,'COLL')
!! do iw=1,nomega
!!   write(msg,'(4es15.7)') omega(iw), real(zeffspec(:,iw))
!!   call wrtout(zeffspec_unit,msg,'COLL')
!! end do
!!
!! write(zeffspec_unit,*) ' '
!! write(zeffspec_unit,'(a)') '#  Imaginary part of Born charges spectral function:'
!! write(zeffspec_unit,*) ' '
!! write(msg,'(a,a)') ch10,&
!! &           ' # At  hw     Z^x     Z^y     Z^z'
!! call wrtout(zeffspec_unit,msg,'COLL')
!! do iw=1,nomega
!!   write(msg,'(4es15.7)') omega(iw), aimag(zeffspec(:,iw))
!!   call wrtout(zeffspec_unit,msg,'COLL')
!! end do
!!
!! close(zeffspec_unit)

 ABI_FREE(dint_fsddb)
 ABI_FREE(int_fsddb)
 ABI_FREE(dummysus)
 ABI_FREE(magsus)
 ABI_FREE(lm_magsus)
 ABI_FREE(invmagsus)
 ABI_FREE(mmom)
 ABI_FREE(mmom_tr)
 ABI_FREE(zfield)
 ABI_FREE(zfield_tr)
 ABI_FREE(barepsilon)
 ABI_FREE(epsilon)
 ABI_FREE(macmagsus)
 ABI_FREE(ifcmat)
 ABI_FREE(ifcmat_fm)
 ABI_FREE(omega)
 ABI_FREE(phfrq)
 ABI_FREE(magphongreen)
 ABI_FREE(phongreen)
 ABI_FREE(phongreen_fm)
 ABI_FREE(magphonspec)
 ABI_FREE(phonspec)
 ABI_FREE(zeffspec)
 ABI_FREE(mmomspec)
 ABI_FREE(mode_magphonspec)
 ABI_FREE(mode_phonspec)
 ABI_FREE(eigvec)
 ABI_FREE(eigvec_fm)
 ABI_FREE(modemm)
 ABI_FREE(fmzeff)
 ABI_FREE(fmzeff_tr)
 ABI_FREE(zeff)
 ABI_FREE(zeff_tr)
 ABI_FREE(genzeff_tr)
 ABI_FREE(lm_epsilon)
 ABI_FREE(ri_magelsus)
 ABI_FREE(lm_magelsus)
 ABI_FREE(ri_genelsus)
 ABI_FREE(modezeff)
 ABI_SFREE(w0hessian)
 ABI_SFREE(w0berry)
 ABI_SFREE(coeffs)

 end subroutine ddb_omega_interpol
!!***

!!****f* ABINIT/lineal_omega_interp
!! NAME
!!  lineal_omega_interp
!!
!! FUNCTION
!!  Calgulates a omega dependent matrix of penalized second-order energies via
!!  \Phi(w)=K+iwG, where K and G are the w=0 penalized Hessians and Berry curvatures.    
!!
!! COPYRIGHT
!!  Copyright (C) 2024 ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert)=  Second-order w=0 derivative matrices
!! blkval_lw(2,3*mpert*3*mpert)=  Third-order w=0 derivative matrices
!!  mpert= maximum number of perturbations
!!  natom= number of atoms in the cell
!!  ndim= dimension of the local spin susceptibilities 
!!  omega= frequency at which the IFCs have been calculated
!!
!! OUTPUT
!!  int_fsddb(2,ddb%msize,1)= interpolated hessian at omega
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

subroutine lineal_omega_interp(blkval,blkval_lw,eta,ifcmat_fm, &
& invmagsus,magsus,mpatpol,mpdir,mpert,msize,natom,ndim,nmdir,int_fsddb, &
& omega,zfield,zfield_tr)

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: mpert,msize,natom,ndim,nmdir 
 real(dp), intent(in) :: eta,omega
!arrays
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp), intent(in) :: blkval(2,3,mpert,3,mpert)
 real(dp), intent(in) :: blkval_lw(2,3,mpert,3,mpert,3,mpert)
 real(dp), intent(out) :: int_fsddb(2,msize,1)
 complex(dpc), intent(out) :: ifcmat_fm(3*natom,3*natom)
 complex(dpc), intent(out) :: invmagsus(ndim,ndim)
 complex(dpc), intent(out) :: magsus(ndim,ndim)
 complex(dpc), intent(out) :: zfield(ndim,(natom+2)*3)
 complex(dpc), intent(out) :: zfield_tr((natom+2)*3,ndim)
 
!Local variables -------------------------
!scalars
 integer :: idir1,idir2,idir3,info,ipert1,ipert2,ipert3
 integer :: iat1,iat2,icol,irow,idir1_red,idir2_red,ipert1_red,ipert2_red
 integer :: lwork
!arrays
 real(dp), allocatable :: lhess(:,:,:,:,:)
 integer, allocatable :: ipiv(:)
 complex(dpc),allocatable :: work(:),work1(:,:)
 
! *********************************************************************

 ABI_MALLOC(lhess,(2,3,mpert,3,mpert))
 ipert3= natom + 9
 idir3= 1
 do ipert2= 1, mpert
   do idir2= 1, 3
     do ipert1= 1, mpert
       do idir1= 1, 3
         lhess(1,idir1,ipert1,idir2,ipert2)= blkval(1,idir1,ipert1,idir2,ipert2) + &
       & omega*blkval_lw(1,idir1,ipert1,idir2,ipert2,idir3,ipert3) - &
       & eta*blkval_lw(2,idir1,ipert1,idir2,ipert2,idir3,ipert3)
         lhess(2,idir1,ipert1,idir2,ipert2)= blkval(2,idir1,ipert1,idir2,ipert2) + &
       & omega*blkval_lw(2,idir1,ipert1,idir2,ipert2,idir3,ipert3) + &
       & eta*blkval_lw(1,idir1,ipert1,idir2,ipert2,idir3,ipert3)
!         lhess(:,idir1,ipert1,idir2,ipert2)= blkval(:,idir1,ipert1,idir2,ipert2) + &
!       & omega*blkval_lw(:,idir1,ipert1,idir2,ipert2,idir3,ipert3) 
       end do
     end do
   end do
 end do
 int_fsddb(1,:,1)= reshape( lhess(1,:,:,:,:), shape = (/msize/) )
 int_fsddb(2,:,1)= reshape( lhess(2,:,:,:,:), shape = (/msize/) )

!Store the relevant matrices for the magnon-phonon Green's function
!First phonon-phonon
 do ipert2= 1, natom
   do idir2= 1, 3
     icol=( ipert2-1)*3 + idir2
     do ipert1= 1, natom
       do idir1= 1, 3
         irow=( ipert1-1)*3 + idir1
         ifcmat_fm(irow,icol)= cmplx(lhess(1,idir1,ipert1,idir2,ipert2), &
       & lhess(2,idir1,ipert1,idir2,ipert2),16)
       end do
     end do
   end do
 end do 

!Then, spin-spin
 ipert2_red= 0
 invmagsus=cmplx(zero,zero,16)
 do iat2= mpatpol(1), mpatpol(2)
   ipert2= natom + 11 + iat2
   ipert2_red= ipert2_red + 1
   idir2_red= 0
   do idir2= 1, 3
     if (mpdir(idir2)==0) cycle
     idir2_red= idir2_red + 1
     icol=idir2_red+(ipert2_red-1)*nmdir
     ipert1_red=0
     do iat1= mpatpol(1), mpatpol(2)
       ipert1= natom + 11 + iat1
       ipert1_red= ipert1_red + 1
       idir1_red= 0
       do idir1= 1, 3
         if (mpdir(idir1)==0) cycle
         idir1_red=idir1_red+1
         irow=idir1_red+(ipert1_red-1)*nmdir
         invmagsus(irow,icol)= cmplx(lhess(1,idir1,ipert1,idir2,ipert2), &
       & lhess(2,idir1,ipert1,idir2,ipert2),16)
       end do
     end do
   end do
 end do

!Invert to obtain magsus
 ABI_MALLOC(work1,(ndim,ndim))
 work1=invmagsus
 
 ABI_MALLOC(ipiv,(ndim))
 call zgetrf( ndim, ndim, work1, ndim, ipiv, info )
 ABI_CHECK(info == 0, sjoin('zgetrf returned:', itoa(info)))

 ABI_MALLOC(work,(2))
 call zgetri( ndim, work1, ndim, ipiv, work, -1, info )
 ABI_CHECK(info == 0, sjoin('zgetri returned:', itoa(info)))
 lwork=int(work(1))

 ABI_REMALLOC(work,(lwork))
 call zgetri( ndim, work1, ndim, ipiv, work, lwork, info )
 ABI_CHECK(info == 0, sjoin('zgetri returned:', itoa(info)))

 magsus=work1
 ABI_FREE(ipiv)
 ABI_FREE(work)
 ABI_FREE(work1)

!Finally, the spin-phonon
 do ipert2=1,natom
   do idir2=1,3
     icol=idir2+(ipert2-1)*3
     ipert1_red= 0
     do iat1= mpatpol(1), mpatpol(2)
       ipert1= natom + 11 + iat1
       ipert1_red= ipert1_red + 1
       idir1_red= 0
       do idir1= 1, 3
         if (mpdir(idir1)==0) cycle
         idir1_red= idir1_red + 1
         irow=idir1_red+(ipert1_red-1)*nmdir
         zfield(irow,icol)= cmplx(lhess(1,idir1,ipert1,idir2,ipert2), &
       & lhess(2,idir1,ipert1,idir2,ipert2),16)
         zfield_tr(icol,irow)= cmplx(lhess(1,idir2,ipert2,idir1,ipert1), &
       & lhess(2,idir2,ipert2,idir1,ipert1),16)
       end do
     end do
   end do
 end do

 ABI_FREE(lhess)

end subroutine lineal_omega_interp
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


subroutine phonon_green(amu,eigvec,eigvec_fm,eta,ifc,ifc_fm,invmagsus,& 
& magphongreen,magphonspec,mode_magphonspec,mode_phonspec,natom,ndim,ntypat,omega,&
& phfrq,phongreen,phongreen_fm,phonspec,typat,zfield,zfield_tr)

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: natom,ndim,ntypat 
 real(dp), intent(in) :: eta,omega
 real(dp), intent(out) :: magphonspec,phonspec
!arrays
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: amu(ntypat)
 real(dp),intent(out) :: eigvec(2*3*natom*3*natom)
 real(dp),intent(out) :: eigvec_fm(2*3*natom*3*natom)
 real(dp),intent(out) :: phfrq(3*natom)
 real(dp), intent(out) :: mode_phonspec(3*natom)
 real(dp), intent(out) :: mode_magphonspec(3*natom+ndim)
 complex(dpc), intent(in) :: ifc(3*natom,3*natom)
 complex(dpc), intent(in) :: ifc_fm(3*natom,3*natom)
 complex(dpc), intent(in) :: invmagsus(ndim,ndim)
 complex(dpc), intent(out) :: magphongreen(3*natom+ndim,3*natom+ndim)
 complex(dpc), intent(out) :: phongreen(3*natom,3*natom)
 complex(dpc), intent(out) :: phongreen_fm(3*natom,3*natom)
 complex(dpc), intent(in) :: zfield(ndim,(natom+2)*3)
 complex(dpc), intent(in) :: zfield_tr((natom+2)*3,ndim)

!Local variables-------------------------------
!scalars
 integer :: iat1,iat2,idir1,idir2,icol,ier,imode,info,irow,lwork,mpdim,pdim
 real(dp) :: mfac1, mfac2
 complex(dpc) :: cplx_eta
!arrays
 integer, allocatable :: ipiv(:)
 real(dp) :: dum(2,0) 
 real(dp), allocatable :: invmassfac(:,:)
 real(dp), allocatable :: matrx(:,:),zhpev1(:,:),zhpev2(:)
 real(dp), allocatable :: eigval(:)
 complex(dpc), allocatable :: dynmat(:,:),w2dynmat(:,:)
 complex(dpc),allocatable :: work(:),work1(:,:)
 complex(dpc),allocatable :: mass_magphongreen(:,:)
!character(len=500) :: msg                   

! *************************************************************************

 DBG_ENTER("COLL")

!Build an array with the inverse mass factors
 ABI_MALLOC(invmassfac,(natom,natom))
 do iat2= 1, natom
   do iat1= 1, natom
     invmassfac(iat1,iat2)=one/sqrt(amu(typat(iat1))*amu(typat(iat2)))/amu_emass
   end do
 end do

!Build the ((w+eta)**2 - D(w)) matrix 
 pdim=3*natom
 ABI_MALLOC(dynmat,(pdim,pdim))
 ABI_MALLOC(w2dynmat,(pdim,pdim))

 cplx_eta=cmplx(0.0_dp,eta)
 do iat2= 1, natom
   do idir2= 1, 3
     icol= (iat2-1)*3 + idir2
     do iat1= 1, natom
       do idir1= 1, 3
         irow= (iat1-1)*3 + idir1
         dynmat(irow,icol)= invmassfac(iat1,iat2)*ifc(irow,icol)
         w2dynmat(irow,icol)= -one*dynmat(irow,icol)
         if (irow==icol) then
           w2dynmat(irow,icol)= (omega+cplx_eta)**2 + w2dynmat(irow,icol)
         end if
       end do
     end do
   end do
 end do

!Invert to obtain the phonon Green's function
 ABI_MALLOC(work1,(pdim,pdim))
 work1=w2dynmat

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

 phongreen=work1

!Finally extract the spectral function from the trace
 do irow= 1, pdim
   mode_phonspec(irow)= -two*omega/pi * aimag(work1(irow,irow))
!   mode_phonspec(irow)= -one/pi * aimag(two*cmplx(omega,eta)*work1(irow,irow))
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

 ! Get the phonon frequencies (negative by convention, if the eigenvalue of the dynamical matrix is negative)
 do imode=1,3*natom
   if(eigval(imode)>=1.0d-16)then
     phfrq(imode)=sqrt(eigval(imode))
   else if(eigval(imode)>=-1.0d-16)then
     phfrq(imode)=zero
   else
     phfrq(imode)=-sqrt(-eigval(imode))
   end if
 end do

 ! Fix the phase of the eigenvectors
 call fxphas_seq(eigvec,dum, 0, 0, 1, 3*natom*3*natom, 0, 3*natom, 3*natom, 0)

 ! Normalise the eigenvectors
 call pheigvec_normalize(natom, eigvec)

 ! Apply mass factos to Green's function to use it later in the calculation of the 
 ! phonon contributions to the susceptibilities.
 do iat2= 1, natom
   do idir2= 1, 3
     icol= (iat2-1)*3 + idir2
     do iat1= 1, natom
       do idir1= 1, 3
         irow= (iat1-1)*3 + idir1
         phongreen(irow,icol)= invmassfac(iat1,iat2)*phongreen(irow,icol)
       end do
     end do
   end do
 end do

!Calcuate also the frozen-spin phonon Greens function
 do iat2= 1, natom
   do idir2= 1, 3
     icol= (iat2-1)*3 + idir2
     do iat1= 1, natom
       do idir1= 1, 3
         irow= (iat1-1)*3 + idir1
         dynmat(irow,icol)= invmassfac(iat1,iat2)*ifc_fm(irow,icol)
         w2dynmat(irow,icol)= -one*dynmat(irow,icol)
         if (irow==icol) then
           w2dynmat(irow,icol)= (omega+cplx_eta)**2 + w2dynmat(irow,icol)
         end if
       end do
     end do
   end do
 end do

!Invert to obtain the phonon Green's function
 work1=w2dynmat

 call zgetrf( pdim, pdim, work1, pdim, ipiv, info )
 ABI_CHECK(info == 0, sjoin('zgetrf returned:', itoa(info)))

 ABI_REMALLOC(work,(2))
 call zgetri( pdim, work1, pdim, ipiv, work, -1, info )
 ABI_CHECK(info == 0, sjoin('zgetri returned:', itoa(info)))
 lwork=int(work(1))

 ABI_REMALLOC(work,(lwork))
 call zgetri( pdim, work1, pdim, ipiv, work, lwork, info )
 ABI_CHECK(info == 0, sjoin('zgetri returned:', itoa(info)))

 phongreen_fm=work1

!Diagonalize the Dynamical matrix
 ABI_MALLOC(matrx,(2,(3*natom*(3*natom+1))/2))
 do icol= 1, pdim
   do irow= 1, icol
     matrx(1,irow + (icol-1)*icol/2)=real(dynmat(irow,icol))
     matrx(2,irow + (icol-1)*icol/2)=aimag(dynmat(irow,icol))
   end do
 end do 
 ABI_FREE(dynmat)

 ABI_MALLOC(zhpev1,(2,2*3*natom-1))
 ABI_MALLOC(zhpev2,(3*3*natom-2))

 call ZHPEV ('V','U',3*natom,matrx,eigval,eigvec_fm,3*natom,zhpev1,zhpev2,ier)
 ABI_CHECK(ier == 0, sjoin('zhpev returned:', itoa(ier)))

 ABI_FREE(matrx)
 ABI_FREE(zhpev1)
 ABI_FREE(zhpev2)

 ! Fix the phase of the eigenvectors
 call fxphas_seq(eigvec_fm,dum, 0, 0, 1, 3*natom*3*natom, 0, 3*natom, 3*natom, 0)

 ! Normalise the eigenvectors
 call pheigvec_normalize(natom, eigvec_fm)


 ! Apply mass factos to Green's function to use it later in the calculation of the 
 ! phonon contributions to the susceptibilities.
 do iat2= 1, natom
   do idir2= 1, 3
     icol= (iat2-1)*3 + idir2
     do iat1= 1, natom
       do idir1= 1, 3
         irow= (iat1-1)*3 + idir1
         phongreen_fm(irow,icol)= invmassfac(iat1,iat2)*phongreen_fm(irow,icol)
       end do
     end do
   end do
 end do

!Now calculate the generalized magnon-phonon Green's function
!Build the (M(w+eta)**2 - C(w)) matrix 
 mpdim= pdim + ndim
 ABI_REMALLOC(w2dynmat,(mpdim,mpdim))
 w2dynmat= (zero,zero)

!First the phonon-phonon sector
 do iat2= 1, natom
   do idir2= 1, 3
     icol= (iat2-1)*3 + idir2
     do iat1= 1, natom
       do idir1= 1, 3
         irow= (iat1-1)*3 + idir1
         w2dynmat(irow,icol)= -one*ifc_fm(irow,icol)
         if (irow==icol) then
           w2dynmat(irow,icol)= amu(typat(iat1))*amu_emass* &
         & (omega+cplx_eta)**2 + w2dynmat(irow,icol)
         end if
       end do
     end do
   end do
 end do

!Next the magnon-magnon sector
 do icol= 1, ndim
   do irow= 1, ndim
     w2dynmat(pdim+irow,pdim+icol)= -invmagsus(irow,icol)
   end do
 end do

!Finally the phonon-magnon and magnon-phonon sector
 do icol= 1, ndim
   do iat1= 1, natom
     do idir1= 1, 3
       irow= (iat1-1)*3 + idir1
       w2dynmat(irow,pdim+icol)= -zfield_tr(irow,icol)
       w2dynmat(pdim+icol,irow)= -zfield(icol,irow)
     end do
   end do
 end do
 
!Invert to obtain the generalized Green's function
 ABI_REMALLOC(work1,(mpdim,mpdim))
 work1=w2dynmat

 ABI_REMALLOC(ipiv,(mpdim))
 call zgetrf( mpdim, mpdim, work1, mpdim, ipiv, info )
 ABI_CHECK(info == 0, sjoin('zgetrf returned:', itoa(info)))

 ABI_REMALLOC(work,(2))
 call zgetri( mpdim, work1, mpdim, ipiv, work, -1, info )
 ABI_CHECK(info == 0, sjoin('zgetri returned:', itoa(info)))
 lwork=int(work(1))

 ABI_REMALLOC(work,(lwork))
 call zgetri( mpdim, work1, mpdim, ipiv, work, lwork, info )
 ABI_CHECK(info == 0, sjoin('zgetri returned:', itoa(info)))

 magphongreen=work1

!Now apply the mass factors
 ABI_MALLOC(mass_magphongreen,(mpdim,mpdim))
 do icol= 1, mpdim
   if (icol <= pdim) then
     iat2= ceiling(icol/three)
     mfac2= sqrt(amu(typat(iat2))*amu_emass)
   else
      mfac2=zero
   end if
   do irow= 1, mpdim
     if (irow <= pdim) then
       iat1= ceiling(irow/three)
       mfac1= sqrt(amu(typat(iat1))*amu_emass)
     else
       mfac1=zero
     end if
     mass_magphongreen(irow,icol)= mfac1*work1(irow,icol)*mfac2
   end do
 end do

!Finally extract the generalized spectral function
 do irow= 1, mpdim
   mode_magphonspec(irow)= -two*omega/pi * aimag(mass_magphongreen(irow,irow))
 end do
 magphonspec= sum(mode_magphonspec(:))

 ABI_FREE(w2dynmat)
 ABI_FREE(ipiv)
 ABI_FREE(work1)
 ABI_FREE(mass_magphongreen)
 ABI_FREE(invmassfac)

 DBG_EXIT("COLL")

end subroutine phonon_green
!!***

!!****f* ABINIT/mode_mmom
!! NAME
!!  mode_mmom
!!
!! FUNCTION
!!  Projects the magnetic moments on the eigenmodes of the dynamical 
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
!!  mmom(ndim,(natom+2)*3)= first-order magnetic moments
!!  natom= number of atoms in the cell
!!  ndim= dimension of the penalized degrees of freedom
!!  ntypat= number of atom types in the cell
!!  typat(natom)= array with the type of atoms in the cell
!!
!! OUTPUT
!!  modemm(ndim,3*natom)= mode-resolved magnetic moments
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


subroutine mode_mmom(amu,eigvec,mmom,mmomspec,modemm,mode_phonspec,natom,ndim,ntypat,typat)

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
 real(dp), intent(in) :: mode_phonspec(3*natom)
 complex(dpc), intent(in) :: mmom(ndim,(natom+2)*3)
 complex(dpc), intent(out) :: modemm(ndim,3*natom)
 complex(dpc), intent(out) :: mmomspec(ndim)

!Local variables-------------------------------
!scalars
 integer :: iat1,iat2,idir1,idir2,icol,im,imode,irow
 real(dp) :: mcell
!arrays
 real(dp), allocatable :: mass(:)
 complex(dpc), allocatable :: mode_mmomspec(:,:)
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

!Compute the mode-resolved moments
 modemm(:,:)=(zero,zero)
 do im= 1, ndim
   do iat2= 1, natom
     do idir2= 1, 3
       imode= (iat2-1)*3 + idir2
       do iat1= 1, natom
         do idir1= 1, 3
           irow= (iat1-1)*3 + idir1
           modemm(im,imode)= modemm(im,imode) +  mass(iat1)*mmom(im,irow)* &
         & cmplx(eigvec(1,idir1,iat1,idir2,iat2),eigvec(2,idir1,iat1,idir2,iat2),16)
         end do
       end do
     end do
   end do
 end do

 ABI_FREE(mass)

!Compute the magnetic moments weighted by the phonon spectral function
 ABI_MALLOC(mode_mmomspec,(ndim,3*natom))
 do irow= 1, 3*natom
   mode_mmomspec(:,irow)= modemm(:,irow)*mode_phonspec(irow)
 end do
 do im= 1, ndim
   mmomspec(im)=sum(mode_mmomspec(im,:))
 end do 
 ABI_FREE(mode_mmomspec)
 

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


subroutine mode_zeff(amu,eigvec,mode_phonspec,modezeff,natom,ntypat,typat,zeff,zeffspec)

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
 real(dp), intent(in) :: mode_phonspec(3*natom)
 complex(dpc), intent(in) :: zeff(3,natom*3)
 complex(dpc), intent(out) :: modezeff(3,3*natom)
 complex(dpc), intent(out) :: zeffspec(3)

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

!Compute the Born charges weighted by the phonon spectral function
 ABI_MALLOC(mode_zeffspec,(3,3*natom))
 do irow= 1, 3*natom
   mode_zeffspec(:,irow)= modezeff(:,irow)*mode_phonspec(irow)
 end do
 zeffspec(1)=sum(mode_zeffspec(1,:))
 zeffspec(2)=sum(mode_zeffspec(2,:))
 zeffspec(3)=sum(mode_zeffspec(3,:))

 ABI_FREE(mode_zeffspec)

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


end module m_ddb_omega_interpol
!!***
