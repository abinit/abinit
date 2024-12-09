!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ddb_magpen
!! NAME
!!  m_ddb_magpen
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2023 ABINIT group (FIXME: add author)
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

module m_ddb_magpen
    
 use defs_basis
 use m_abicore
 use m_profiling_abi
 use m_errors
 use m_ddb
 use m_fstrings,        only : itoa, sjoin
 use m_macroave,        only : POLINT
 use m_io_tools,        only : open_file
 use m_cgtools,         only : fxphas_seq
 use m_dynmat,          only : pheigvec_normalize
 use m_numeric_tools,   only : polcoe

 implicit none

 public :: ddb_magpen       ! Convert the derivatives calculated with the magnetic penalty into the physically relevant ones.
 public :: ddb_omega_interpol ! Perform an interpolation of the he derivatives calculated with the magnetic
                              ! penalty and later on convert them into the physically relevant ones at each value of interpolated omega.


 private

! *************************************************************************

contains 
!!***

!!****f* m_ddb_magpen/ddb_magpen
!! NAME
!! ddb_magpen
!!
!! FUNCTION
!! Convert the second- and possibly third-order derivatives calculated with the magnetic penalty
!! into physically relevant quantities.
!!
!! INPUTS
!! ddb (INOUT) = ddb block datastructure
!! ddb_lw (INOUT) = ddb_lw block datastructure
!! magpen = amplitude (in Ha) of the applied magnetic penalty 
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert = maximum number of ipert
!! mpopt = 1 calculate the frozen-magnetic second-order quantities
!!         2 calculate the spin-relaxed second-order quantities 
!! natom= number of atoms in unit cell
!! ntypat= number of atom types
!! timdisp= 1 calculate the third-order Berry curvatures
!! ucvol= unit cell volume
!!
!! OUTPUT
!! ddb= ddb%val updated with the corrected second-order derivatives
!!
!! SOURCE

 subroutine ddb_magpen(ddb,ddb_lw,delta_asrw0,delta_asrw0_fm,dissip,& 
& magpen,mpatpol,mpdir,mpert,mpopt,natom, &
& ntypat,omegaflag,prtvol,rftyp,ucvol,timdisp,xred)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: dissip,mpert,mpopt,natom,ntypat,omegaflag,prtvol,rftyp,timdisp
 real(dp),intent(in) :: magpen,ucvol
!arrays
 type(ddb_type),intent(inout) :: ddb,ddb_lw
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp), intent(inout) :: delta_asrw0(3*natom,3), delta_asrw0_fm(3*natom,3)
 real(dp),intent(in) :: xred(3,natom)

!Local variables -------------------------
!scalars
 integer :: iblok,ii,ipert1,ipert2,jblok,kblok,lblok,nblok,ndim 
 integer :: nmat,nmdir,prtopt
 character(len=500) :: msg
 logical :: qeq0
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4),rfmagn(4),rffreq(4)
 real(dp) :: omega(3),qphnrm(3),qphon(3,3)
 complex(dpc) :: barepsilon(3,3),epsilon(3,3), macmagsus(3,3)
 complex(dpc), allocatable :: barmagsus(:,:),invbarmagsus(:,:)
 complex(dpc), allocatable :: invmagsus(:,:), magsus(:,:), invhmat(:,:)
 complex(dpc), allocatable :: barmmom(:,:),barmmom_tr(:,:),mmom(:,:),mmom_tr(:,:)
 complex(dpc), allocatable :: zfield(:,:), zfield_tr(:,:)
 complex(dpc), allocatable :: bc_barmagsus(:,:),bc_ss(:,:),bc_sp(:,:)
 complex(dpc), allocatable :: ifcmat(:,:),ifcmat_fm(:,:),zeff(:,:),zeff_tr(:,:)
 complex(dpc), allocatable :: fmzeff(:,:),fmzeff_tr(:,:)
 complex(dpc), allocatable :: lm_epsilon(:,:),dum_phongreen(:,:)

! *********************************************************************
 write(msg, '(2a,(80a),4a)' ) ch10,('=',ii=1,80),ch10,ch10,&
 ' Linear-response magnetic penalty section ',ch10
 call wrtout([std_out, ab_out], msg)

 prtopt=1
 if (magpen<zero) then
   nmat= 1
 else if (magpen>zero) then
   nmat= mpatpol(2) - mpatpol(1) + 1
 end if
 nmdir=sum(mpdir(:))
 ndim=nmat*nmdir
 ABI_MALLOC(barmagsus,(ndim,ndim))
 ABI_MALLOC(magsus,(ndim,ndim))
 ABI_MALLOC(invbarmagsus,(ndim,ndim))
 ABI_MALLOC(invmagsus,(ndim,ndim))
 ABI_MALLOC(invhmat,(ndim,ndim))
 ABI_MALLOC(barmmom,(ndim,(natom+5)*3))
 ABI_MALLOC(barmmom_tr,((natom+5)*3,ndim))
 ABI_MALLOC(mmom,(ndim,(natom+5)*3))
 ABI_MALLOC(mmom_tr,((natom+5)*3,ndim))
 ABI_MALLOC(zfield,(ndim,(natom+5)*3))
 ABI_MALLOC(zfield_tr,((natom+5)*3,ndim))

 ABI_MALLOC(ddb%val_fs,(2,ddb%msize,ddb%nblok))
 if (mpopt==2) ABI_MALLOC(ddb%val_rs,(2,ddb%msize,ddb%nblok))

 ABI_MALLOC(ifcmat,(3*natom,3*natom))
 ABI_MALLOC(ifcmat_fm,(3*natom,3*natom))
 ABI_MALLOC(fmzeff,(3,3*natom))
 ABI_MALLOC(fmzeff_tr,(3*natom,3))
 ABI_MALLOC(zeff,(3,3*natom))
 ABI_MALLOC(zeff_tr,(3*natom,3))
 ABI_MALLOC(lm_epsilon,(3,3))
 ABI_MALLOC(dum_phongreen,(3*natom,3*natom))

 nblok=ddb%nblok
 do kblok=1,nblok

   ! Look for the local spin-susceptibility block in the DDB
   omega=zero
   qphon=zero
   qphon(:,1)=ddb%qpt(1:3,kblok)
   qeq0=(sqrt(sum(qphon(:,1)**2))<tol8)
   qphnrm(:)=ddb%nrm(1,kblok)
   omega(1)=ddb%omega(1,kblok)
   rfphon(1:2)=0
   rfelfd(1:2)=0
   rfstrs(1:2)=0
   if (magpen<zero) then
     rfmagn(1:2)= 1
   else if (magpen>zero) then
     rfmagn(1:2)= 2
   end if

   write(msg, '(1a,(80a),2a,3f16.8,2a,f16.8,a)' ) ch10,('-',ii=1,80),ch10, &
   ' q point  ', qphon(:,1),ch10,&
   ' frequency', ddb%omega(1,kblok), ch10
   call wrtout([std_out, ab_out], msg)

   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, &
  & mpatpol=mpatpol,mpdir=mpdir,omega=omega,rfmagn=rfmagn)

   ! Calculate and write the spin-susceptibility matrices
   if (iblok /= 0) then
     if (prtvol>1) then
       if (magpen<zero) then
         write(msg, '(2a)' ) ' Spin susceptibility (Uniform Zeeman) ',ch10
         call wrtout([std_out, ab_out], msg)
       else if (magpen>zero) then
         write(msg, '(2a)' ) ' Spin susceptibility (Local Zeeman) ',ch10
         call wrtout([std_out, ab_out], msg)
       end if
     end if

     call local_spinsus(barmagsus,ddb,iblok,invbarmagsus,invmagsus,invhmat,magpen,magsus,&
   & mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,omega(1),omegaflag,prtopt,prtvol)

   end if

   ! Calculate and write the induced magnetic moments 
   if (prtvol>1) then
     write(msg, '(2a,(80a),4a)' ) ch10,('-',ii=1,80),ch10,ch10,&
     ' First-order magnetic moments ',ch10
     call wrtout([std_out, ab_out], msg)
   end if

   ! First atomic-displacement
   rfphon(2)=1
   rfelfd(1:2)=0
   rfstrs(1:2)=0
   rfmagn(:)=0
   if (magpen<zero) then
     rfmagn(1)= 1
   else if (magpen>zero) then
     rfmagn(1)= 2
   end if

   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, &
  & mpatpol=mpatpol,mpdir=mpdir,omega=omega,rfmagn=rfmagn)

   ! Then electric field
   ! Look for the induced magnetic moments block in the DDB
   jblok=0
   if (qeq0) then
     rfphon(2)=0
     rfelfd(2)=2
     rfstrs(1:2)=0
     rfmagn(:)=0
     if (magpen<zero) then
       rfmagn(1)= 1
     else if (magpen>zero) then
       rfmagn(1)= 2
     end if

     call ddb%get_block(jblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, &
   & mpatpol=mpatpol,mpdir=mpdir,omega=omega,rfmagn=rfmagn)
   end if

   ! Then macroscopic Zeeman field
   ! Look for the induced magnetic moments block in the DDB
   lblok=0
   if (qeq0) then
     rfphon(:)=0
     rfelfd(:)=0
     rfstrs(1:2)=0
     rfmagn(2)=1
     if (magpen<zero) then
       rfmagn(1)= 1
     else if (magpen>zero) then
       rfmagn(1)= 2
     end if

     call ddb%get_block(lblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, &
   & mpatpol=mpatpol,mpdir=mpdir,omega=omega,rfmagn=rfmagn)
   end if

   if (iblok /= 0 .or. jblok /=0 .or. lblok/=0) then
     call magmom(barmmom,barmmom_tr,ddb,invbarmagsus,invhmat,iblok,jblok,lblok,magpen,magsus,mmom,mmom_tr,&
   & mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,omega(1),omegaflag,prtopt,prtvol,qphon,xred,zfield,zfield_tr)
   end if

   !Now calculate the non-magnetic second-order quantities
   write(msg, '(2a,(80a),4a)' ) ch10,('-',ii=1,80),ch10,ch10,&
   ' Second-order linear-response tensors ',ch10
   call wrtout([std_out, ab_out], msg)

   !Convert ddb%val to second-order energies
   call ddb%to_d2etot(ddb%val,kblok,0,qeq0,qphon,qphnrm,ucvol,omega=omega)

   !Convert second-order derivatives to diferent magnetic boundary conditions
   call mp_d2etot(barmagsus,ddb,0,kblok,&
 & invhmat,magsus,magpen,mpert,mpopt,natom,nblok,ndim,qphon,xred,zfield,zfield_tr)

   !Convert second-order energies to the physical quantities of ddb%val
   call ddb%to_d2etot(ddb%val,kblok,1,qeq0,qphon,qphnrm,ucvol,omega=omega)
   call ddb%to_d2etot(ddb%val_fs,kblok,1,qeq0,qphon,qphnrm,ucvol,omega=omega)
   if (mpopt==2) call ddb%to_d2etot(ddb%val_rs,kblok,1,qeq0,qphon,qphnrm,ucvol,omega=omega)

   !Print the physical quantities in the new magnetic boundary conditions
   if (prtopt==1) then
     if (mpopt==1) then
       call mp_d2etot_print(ddb,ddb%val_fs,kblok,mpert,natom,nblok,1,omega,prtvol,qeq0,qphnrm,qphon,ucvol)
     else if (mpopt==2) then
       call mp_d2etot_print(ddb,ddb%val_rs,kblok,mpert,natom,nblok,2,omega,prtvol,qeq0,qphnrm,qphon,ucvol)
     end if
   end if

   rfmagn(:)=0
   rfelfd(:)=0
   rfphon(:)=0
   
!   !IFCs block
!   rfphon(1:2)=1
!   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega)
!   if (iblok /= 0) then
!     call mp_ifc(barmagsus,barmmom,barmmom_tr,ddb%val,0,iblok,ifcmat,&
!   & ifcmat_fm,invhmat,magsus,magpen,mpert,mpopt,&
!   & natom,nblok,ndim,omega(1),omegaflag,prtopt,prtvol,qphon,xred,zfield,zfield_tr)
!     
!     !Apply ASR
!     if (omega(1) < tol12) then
!       call asrw0(delta_asrw0,ifcmat,natom,0) 
!       call asrw0(delta_asrw0_fm,ifcmat_fm,natom,0) 
!     else 
!       call asrw0(delta_asrw0,ifcmat,natom,1) 
!       call asrw0(delta_asrw0_fm,ifcmat_fm,natom,1) 
!     end if
!   end if
!
!   !Born effective charges block
!   jblok=0
!   if (qeq0) then
!     rfphon(1:2)=1
!     rfelfd(1:2)=2
!     call ddb%get_block(jblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega)
!   end if
!   if (jblok /= 0 ) then
!     call mp_zeff(barmagsus,barmmom,barmmom_tr,ddb%val,&
!   & 0,fmzeff,fmzeff_tr,jblok,invhmat,lm_epsilon,magpen,magsus,mpert,mpopt,&
!   & natom,nblok,ndim,dum_phongreen,prtopt,prtvol,ucvol,zeff,zeff_tr,zfield,zfield_tr)
!   end if
!
!   !Dielectric susceptibility block
!   lblok=0
!   if (qeq0) then
!     rfphon(:)=0
!     rfelfd(1:2)=2
!     call ddb%get_block(lblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega)
!   end if
!   if (lblok /= 0 ) then
!     call mp_diel(barepsilon,barmagsus,barmmom,barmmom_tr,ddb%val,&
!   & 0,epsilon,lblok,invhmat,magpen,magsus,mpert,mpopt,&
!   & natom,nblok,ndim,prtopt,prtvol,ucvol,zfield,zfield_tr)
!   end if
!
!   !Magnetic susceptibility block
!   iblok=0
!   rfphon(:)=0
!   rfelfd(1:2)=0
!   rfstrs(1:2)=0
!   rfmagn(1:2)=1
!   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega)
!   if (lblok /= 0 ) then
!     call mp_macmagsus(barmagsus,ddb%val,&
!   & 0,iblok,invbarmagsus,invhmat,macmagsus,magpen,magsus,mpatpol,mpdir,mpert,mpopt,&
!   & natom,nblok,ndim,nmdir,prtopt,prtvol,ucvol)
!   end if

 end do

 ! BERRY CURVATURES
 if (timdisp==1) then

!   if (dissip==1) then
!     ABI_BUG("Berry curvatures calculation is not implemented with dissipation, set dissip=0")
!   end if

   ABI_MALLOC(bc_barmagsus,(ndim,ndim))
   ABI_MALLOC(bc_ss,(ndim,ndim))
   ABI_MALLOC(bc_sp,(ndim,(natom+2)*3))

   write(msg, '(2a,(80a),4a)' ) ch10,('=',ii=1,80),ch10,ch10,&
   ' Frequency-derivatives magnetic penalty section ',ch10
   call wrtout([std_out, ab_out], msg)

   rffreq(:)=0
   nblok=ddb_lw%nblok
   do kblok=1,nblok

     if (ddb_lw%typ(kblok)/=33) cycle

     !Berry curvature of the penalized spin-susceptibility
     qphon=zero
     qphon(:,1)=ddb_lw%qpt(1:3,kblok)
     qphnrm(:)=ddb_lw%nrm(1,kblok)
     omega(:)=ddb_lw%omega(:,kblok)
     rfphon(1:3)=0
     rfelfd(1:3)=0
     rfstrs(1:3)=0
     rffreq(3)=1
     iblok=0
     if (magpen<zero) then
       rfmagn(1:2)= 1
     else if (magpen>zero) then
       rfmagn(1:2)= 2
     end if

     write(msg, '(a,(80a),a,3(a,3f16.8,a),3(a,f16.8,a))' ) ch10,('-',ii=1,80),ch10, &
     ' q point 1  ', qphon(:,1),ch10,&
     ' q point 2  ', qphon(:,2),ch10,&  
     ' q point 3  ', qphon(:,3),ch10,&  
     ' frequency 1', ddb_lw%omega(1,kblok), ch10,&
     ' frequency 2', ddb_lw%omega(2,kblok), ch10,&
     ' frequency 3', ddb_lw%omega(3,kblok), ch10
     call wrtout([std_out, ab_out], msg)

     call ddb_lw%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, 33, &
   & mpatpol=mpatpol,mpdir=mpdir,omega=omega,rfmagn=rfmagn,rffreq=rffreq)

     if (iblok /= 0) then
       call berrycurv_ss(bc_barmagsus,bc_ss,ddb_lw%val,iblok,invbarmagsus,mpatpol,mpdir,mpert,&
     & natom,nblok,ndim,nmdir,omega(1),omegaflag,prtvol)
     end if

     !Berry curvature of the induced Zeeman fields

     !First atomic-displacement
     iblok=0
     rfphon(2)=1
     rfmagn(:)=0
     if (magpen<zero) then
       rfmagn(1)= 1
     else if (magpen>zero) then
       rfmagn(1)= 2
     end if

     call ddb_lw%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, 33, &
   & mpatpol=mpatpol,mpdir=mpdir,omega=omega,rfmagn=rfmagn,rffreq=rffreq)

     ! Then electric field
     jblok=0
     if (qeq0) then
       rfphon(2)=0
       rfelfd(2)=2

       call ddb_lw%get_block(jblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, 33, &
     & mpatpol=mpatpol,mpdir=mpdir,omega=omega,rfmagn=rfmagn,rffreq=rffreq)
     end if

     if (iblok /= 0 .or. jblok /=0) then
       call berrycurv_sp(barmmom,bc_sp,bc_ss,ddb_lw%val,iblok,invbarmagsus,jblok, &
     & mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,omega(1),omegaflag,prtvol,qphon,xred)
     end if

     !Berry curvature of other second-order quantites

     !IFCs
     rfphon(1:2)=1
     rfelfd(:)=0
     rfmagn(:)=0
     call ddb_lw%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, 33, &
   & omega=omega,rffreq=rffreq)
     if (iblok /= 0 ) then
       call berrycurv_pp(barmagsus,bc_barmagsus,bc_sp,ddb_lw%val,iblok, &
     & mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,omega(1),omegaflag,prtvol,qphon,xred,zfield)
     end if
   end do 

   ABI_FREE(bc_ss)
 end if

!Deallocations
 ABI_FREE(barmagsus)
 ABI_FREE(barmmom)
 ABI_FREE(barmmom_tr)
 ABI_FREE(magsus)
 ABI_FREE(invbarmagsus)
 ABI_FREE(invmagsus)
 ABI_FREE(mmom)
 ABI_FREE(mmom_tr)
 ABI_FREE(zfield)
 ABI_FREE(zfield_tr)
 ABI_FREE(ifcmat)
 ABI_FREE(ifcmat_fm)
 ABI_FREE(zeff)
 ABI_FREE(zeff_tr)
 ABI_FREE(fmzeff)
 ABI_FREE(fmzeff_tr)
 ABI_FREE(lm_epsilon)
 ABI_FREE(dum_phongreen)

 end subroutine ddb_magpen
!!***

!!****f* m_ddb_magpen/local_spinsus
!! NAME
!! local_spinsus
!!
!! FUNCTION
!! Calculate the spin-susceptibility matrix and its inverse
!!
!! INPUTS
!! ddb=  Second-order derivative arrais
!! iblok= index of the current block
!! magpen = amplitude (in Ha) of the applied magnetic penalty 
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! nmdir= number of directions along which the magnetic penalty was applied
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!!
!! OUTPUT
!! barmagsus(ndim,ndim)= Penalized spin-sussceptibility tensor
!! invbarmagsus(ndim,ndim)= Inverse of the penalized spin-sussceptibility tensor
!! magsus(ndim,ndim)= Spin-sussceptibility tensor
!! invmagsus(ndim,ndim)= Inverse of spin-sussceptibility tensor
!!
!! SOURCE

 subroutine local_spinsus(barmagsus,ddb,iblok,invbarmagsus,invmagsus,invhmat,magpen,magsus,&
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,omega,omegaflag,prtopt,prtvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok,ndim,nmdir
 integer,intent(in) :: omegaflag,prtopt,prtvol
 real(dp),intent(in) :: magpen,omega
!arrays
 type(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: mpatpol(2),mpdir(3)
 complex(dpc),intent(out) :: barmagsus(ndim,ndim)
 complex(dpc),intent(out) :: invbarmagsus(ndim,ndim)
 complex(dpc),intent(out) :: invhmat(ndim,ndim)
 complex(dpc),intent(out) :: magsus(ndim,ndim)
 complex(dpc),intent(out) :: invmagsus(ndim,ndim)

!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,index,info,ipert1,ipert2,irow,lwork
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
 real(dp) :: fac
 character(len=1000) :: msg
!arrays
 complex(dpc) :: idty(ndim,ndim)
 integer(dp) :: indexat(ndim),indexdir(ndim)
 integer, allocatable :: ipiv(:)
 complex(dpc),allocatable :: work(:),work1(:,:),work2(:,:)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

 !Extract the penalized susceptibility
 idty=(zero,zero)
 ipert2_red= 0
 do iat2= mpatpol(1), mpatpol(2)
   ipert2= natom + 11 + iat2
   ipert2_red= ipert2_red + 1
   idir2_red= 0
   do idir2= 1, 3
     if (mpdir(idir2)==0) cycle
     idir2_red= idir2_red + 1
     icol=idir2_red+(ipert2_red-1)*nmdir
     indexat(icol)=iat2
     indexdir(icol)=idir2
     idty(icol,icol)=(one,zero)
     ipert1_red=0
     do iat1= mpatpol(1), mpatpol(2)
       ipert1= natom + 11 + iat1
       ipert1_red= ipert1_red + 1
       idir1_red= 0
       do idir1= 1, 3
         if (mpdir(idir1)==0) cycle
         idir1_red=idir1_red+1
         irow=idir1_red+(ipert1_red-1)*nmdir
         index= idir1 + 3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))

         barmagsus(irow,icol)= &
       & cmplx(ddb%val(1,index,iblok),ddb%val(2,index,iblok),16)

       end do
     end do
   end do
 end do

!Use magsus to store the intermediate array
 magsus=idty-magpen*barmagsus
     
!Invert the arrays
 ABI_MALLOC(work1,(ndim,ndim))
 ABI_MALLOC(work2,(ndim,ndim))
 work1=barmagsus
 work2=magsus

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


 call zgetrf( ndim, ndim, work2, ndim, ipiv, info )
 ABI_CHECK(info == 0, sjoin('zgetrf returned:', itoa(info)))

 call zgetri( ndim, work2, ndim, ipiv, work, -1, info )
 ABI_CHECK(info == 0, sjoin('zgetri returned:', itoa(info)))
 lwork=int(work(1))

 ABI_REMALLOC(work,(lwork))
 call zgetri( ndim, work2, ndim, ipiv, work, lwork, info )
 ABI_CHECK(info == 0, sjoin('zgetri returned:', itoa(info)))
 ABI_FREE(work)

 !Write the results in meaningfull arrays
 !use idty here to store an intermediate array
 invbarmagsus=work1
 
 !At last, calculate the susceptibility and its inverse
 invhmat=work2
 magsus=matmul(work2,barmagsus)
 invmagsus=invbarmagsus-magpen*idty

 ABI_FREE(ipiv)
 ABI_FREE(work1)
 ABI_FREE(work2)

 if (prtopt>1.and.prtvol>1) then
!TODO: remove
!  ! fac=2.714943600699**2*27.2114/four
!   fac=27.2114/four
!   open(10,file='k_ss.txt')
!     do irow=1, ndim
!       write(10,*) invmagsus(irow,1:ndim)*fac
!     end do 
!   close(10)
  
   !Write results in output
   if (magpen > zero) then
     call wrtout([ab_out,std_out], ' Local spin susceptibility ')
     call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
     do irow=1, ndim
       do icol=1, ndim
         write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
       & indexat(irow), cart(indexdir(irow)), indexat(icol), cart(indexdir(icol)), &
       & real(magsus(irow,icol)), aimag(magsus(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
  
  
     call wrtout([ab_out,std_out], ' Inverse of local spin susceptibility ')
     call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
     do irow=1, ndim
       do icol=1, ndim
         write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
       & indexat(irow), cart(indexdir(irow)), indexat(icol), cart(indexdir(icol)), &
       & real(invmagsus(irow,icol)), aimag(invmagsus(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
  
     if (prtvol > 2) then
       call wrtout([ab_out,std_out], ' Penalized local spin susceptibility ')
       call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
       do irow=1, ndim
         do icol=1, ndim
           write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
         & indexat(irow), cart(indexdir(irow)), indexat(icol), cart(indexdir(icol)), &
         & real(barmagsus(irow,icol)), aimag(barmagsus(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
    
       call wrtout([ab_out,std_out], ' Inverse of penalized local spin susceptibility ')
       call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
       do irow=1, ndim
         do icol=1, ndim
           write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
         & indexat(irow), cart(indexdir(irow)), indexat(icol), cart(indexdir(icol)), &
         & real(invbarmagsus(irow,icol)), aimag(invbarmagsus(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
  
       call wrtout([ab_out,std_out], ' Inverse of H matrix (I-\alpha \barchi)^{-1} ')
       call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
       do irow=1, ndim
         do icol=1, ndim
           write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
         & indexat(irow), cart(indexdir(irow)), indexat(icol), cart(indexdir(icol)), &
         & real(invhmat(irow,icol)), aimag(invhmat(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
     end if
  
   end if !magpen>zero

 end if !prtopt

!Store the FS and RS flavors on the DDB array
 ipert2_red= 0
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
         index= idir1 + 3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))
 
         ddb%val_fs(1,index,iblok)=real(invmagsus(irow,icol))
         ddb%val_fs(2,index,iblok)=aimag(invmagsus(irow,icol))

         ddb%val_rs(1,index,iblok)=real(magsus(irow,icol))
         ddb%val_rs(2,index,iblok)=aimag(magsus(irow,icol))
 
       end do
     end do
   end do
 end do
 
 end subroutine local_spinsus
!!***

!!****f* m_ddb_magpen/magmom
!! NAME
!! magmom
!!
!! FUNCTION
!! Calculate the first-order magnetic moments and the constrained
!! Zeeman fields
!!
!! INPUTS
!! ddb=  Second-order derivative arrais
!! invbarmagsus(ndim,ndim)= Inverse of the penalized spin-sussceptibility tensor
!! iblok= index of the atomic displacement block
!! jblok= index of the electric field block
!! magsus(ndim,ndim)= Spin-sussceptibility tensor
!! magpen = amplitude (in Ha) of the applied magnetic penalty 
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! nmdir= number of directions along which the magnetic penalty was applied
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!!
!! OUTPUT
!! barmmom(ndim,(natom+2)*3)= penalized first order magnetic moments on the atoms and 
!!  directions of the penalty induced by atomic displacements and/or electric fields.
!! mmom(ndim,(natom+2)*3)= first order magnetic moments on the atoms and 
!!  directions of the penalty induced by atomic displacements and/or electric fields.
!! zfield(ndim,(natom+2)*3)= Zeeman fields at constrained magnetic moments.
!!
!! SOURCE

 subroutine magmom(barmmom,barmmom_tr,ddb,invbarmagsus,invhmat,iblok,jblok,lblok,magpen,magsus,mmom,mmom_tr,&
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,omega,omegaflag,prtopt,prtvol,qphon,xred,zfield,zfield_tr)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,jblok,lblok,mpert,natom,nblok,ndim,nmdir,omegaflag,prtopt,prtvol
 real(dp),intent(in) :: omega,magpen
!arrays
 type(ddb_type),intent(inout) :: ddb
 real(dp),intent(in) :: qphon(3),xred(3,natom)
 integer,intent(in) :: mpatpol(2),mpdir(3)
 complex(dpc),intent(out) :: barmmom(ndim,(natom+5)*3)
 complex(dpc),intent(out) :: barmmom_tr((natom+5)*3,ndim)
 complex(dpc),intent(in) :: invbarmagsus(ndim,ndim)
 complex(dpc),intent(in) :: invhmat(ndim,ndim)
 complex(dpc),intent(in) :: magsus(ndim,ndim)
 complex(dpc),intent(out) :: mmom(ndim,(natom+5)*3)
 complex(dpc),intent(out) :: mmom_tr((natom+5)*3,ndim)
 complex(dpc),intent(out) :: zfield(ndim,(natom+5)*3)
 complex(dpc),intent(out) :: zfield_tr((natom+5)*3,ndim)
!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,index,ipert1,ipert2,irow
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red,jndex,zblok
 real(dp) :: fac
 character(len=1000) :: msg
!arrays
 integer(dp) :: indexat1(ndim),indexdir1(ndim)
 integer(dp) :: indexat2((natom+5)*3),indexdir2((natom+5)*3)
 complex(dpc) :: mmom_alt(ndim,(natom+5)*3)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

!Extract the penalized moments
 do ipert2=1,natom+5
   !exclude strain perturbation
   if (ipert2==natom+3.or.ipert2==natom+4) cycle
   do idir2=1,3
     icol=idir2+(ipert2-1)*3
     indexat2(icol)=ipert2
     indexdir2(icol)=idir2

     ipert1_red= 0
     do iat1= mpatpol(1), mpatpol(2)
       ipert1= natom + 11 + iat1
       ipert1_red= ipert1_red + 1
       idir1_red= 0
       do idir1= 1, 3
         if (mpdir(idir1)==0) cycle
         idir1_red= idir1_red + 1
         irow=idir1_red+(ipert1_red-1)*nmdir
         indexat1(irow)=iat1
         indexdir1(irow)=idir1
         index= idir1 + 3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))
         jndex= idir2 + 3*((ipert2-1)+mpert*((idir1-1)+3*(ipert1-1)))

         if (iblok /=0 .and. ipert2 <= natom) then
           barmmom(irow,icol)= cmplx(ddb%val(1,index,iblok),ddb%val(2,index,iblok),16)
           barmmom_tr(icol,irow)= cmplx(ddb%val(1,jndex,iblok),ddb%val(2,jndex,iblok),16)
         else if (jblok /=0 .and. ipert2 == natom+2) then
           barmmom(irow,icol)= cmplx(ddb%val(1,index,jblok),ddb%val(2,index,jblok),16)
           barmmom_tr(icol,irow)= cmplx(ddb%val(1,jndex,jblok),ddb%val(2,jndex,jblok),16)
         else if (lblok /=0 .and. ipert2 == natom+5) then
           barmmom(irow,icol)= cmplx(ddb%val(1,index,lblok),ddb%val(2,index,lblok),16)
           barmmom_tr(icol,irow)= cmplx(ddb%val(1,jndex,lblok),ddb%val(2,jndex,lblok),16)
         end if

       end do
     end do
   end do
 end do

!Compute the Zeeman fields 
 zfield=-matmul(invbarmagsus,barmmom)
 zfield_tr=-matmul(barmmom_tr,invbarmagsus)

!Compute the moments
 mmom=-matmul(magsus,zfield)
 mmom_alt=matmul(invhmat,barmmom)
 mmom_tr=matmul(barmmom_tr,invhmat)

 if (prtopt==1.and.prtvol>1) then
    !TODO: the change of phase should rather be done on the magnetic variables, for them
    !to follow the same criterion as the atomic displacement ones. 
!  ! fac=2.714943600699/two*27.2114/0.529177
!   fac=27.2114/0.529177/two
!  
!   open(10,file='k_ps.txt')
!   do iat1= 1, natom
!     do idir1= 1, 3
!       icol= (iat1-1)*3 + idir1
!       !MR: caution, this conjg might be incorrect in presence of dissipation
!       write(10,*) conjg(zfield(1:ndim,icol)*fac* &
!     & exp(two_pi*(0.d0,1.d0)* dot_product(qphon,xred(:,iat1))))
!     end do
!   end do 
!   close(10)
  
  !Write the results
   if (magpen > zero) then

     !Atomic displacements
     if (iblok /= 0) then
       call wrtout([ab_out,std_out], ' Local Zeeman fields induced by atomic displacements (at constrained magnetic moments)')
       call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
       do irow=1, ndim
         do icol=1, natom*3
           write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
         & indexat1(irow), cart(indexdir1(irow)), indexat2(icol), cart(indexdir2(icol)), &
         & real(zfield(irow,icol)), aimag(zfield(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
  
       call wrtout([ab_out,std_out], ' Local magnetic moments induced by atomic displacements (from induced Zeeman fields)')
       call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
       do irow=1, ndim
         do icol=1, natom*3
           write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
         & indexat1(irow), cart(indexdir1(irow)), indexat2(icol), cart(indexdir2(icol)), &
         & real(mmom(irow,icol)), aimag(mmom(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
  
       call wrtout([ab_out,std_out], ' Local magnetic moments induced by atomic displacements (from induced penalized moments)')
       call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
       do irow=1, ndim
         do icol=1, natom*3
           write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
         & indexat1(irow), cart(indexdir1(irow)), indexat2(icol), cart(indexdir2(icol)), &
         & real(mmom_alt(irow,icol)), aimag(mmom_alt(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
     end if
     
     !Electric field
     if (jblok /= 0) then
       call wrtout([ab_out,std_out], ' Local Zeeman fields induced by electric field (at constrained magnetic moments)')
       call wrtout([ab_out,std_out], '  atom1  dir  efld.dir         Real              Imag')
       do irow=1, ndim
         do icol=(natom+1)*3+1, (natom+2)*3
           write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
         & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
         & real(zfield(irow,icol)), aimag(zfield(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
       call wrtout([ab_out,std_out], ' Local magnetic moments induced by electric field (from induced Zeeman fields)')
       call wrtout([ab_out,std_out], '  atom1  dir  efld.dir         Real              Imag')
       do irow=1, ndim
         do icol=(natom+1)*3+1, (natom+2)*3
           write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
         & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
         & real(mmom(irow,icol)), aimag(mmom(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
       call wrtout([ab_out,std_out], ' Local magnetic moments induced by electric field (from induced penalized moments)')
       call wrtout([ab_out,std_out], '  atom1  dir  efld.dir         Real              Imag')
       do irow=1, ndim
         do icol=(natom+1)*3+1, (natom+2)*3
           write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
         & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
         & real(mmom_alt(irow,icol)), aimag(mmom_alt(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
     end if
     
     !Macroscopic Zeeman
     if (lblok /= 0) then
       call wrtout([ab_out,std_out], ' Local Zeeman fields induced by macroscopic Zeeman field (at constrained magnetic moments)')
       call wrtout([ab_out,std_out], '  atom1  dir  Bfld.dir         Real              Imag')
       do irow=1, ndim
         do icol=(natom+4)*3+1, (natom+5)*3
           write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
         & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
         & real(zfield(irow,icol)), aimag(zfield(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
       call wrtout([ab_out,std_out], ' Local magnetic moments induced by macroscopic Zeeman field (from induced Zeeman fields)')
       call wrtout([ab_out,std_out], '  atom1  dir  Bfld.dir         Real              Imag')
       do irow=1, ndim
         do icol=(natom+4)*3+1, (natom+5)*3
           write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
         & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
         & real(mmom(irow,icol)), aimag(mmom(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
       call wrtout([ab_out,std_out], ' Local magnetic moments induced by macroscopic Zeeman field (from induced penalized moments)')
       call wrtout([ab_out,std_out], '  atom1  dir  Bfld.dir         Real              Imag')
       do irow=1, ndim
         do icol=(natom+4)*3+1, (natom+5)*3
           write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
         & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
         & real(mmom_alt(irow,icol)), aimag(mmom_alt(irow,icol))
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], '   ')
     end if

     if (prtvol > 2) then
       if (iblok /= 0) then
         call wrtout([ab_out,std_out], ' Penalized local magnetic moments induced by atomic displacements ')
         call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
         do irow=1, ndim
           do icol=1, natom*3
             write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
           & indexat1(irow), cart(indexdir1(irow)), indexat2(icol), cart(indexdir2(icol)), &
           & real(barmmom(irow,icol)), aimag(barmmom(irow,icol))
             call wrtout([ab_out,std_out], msg)
           end do
         end do
         call wrtout([ab_out,std_out], '   ')
       end if
       if (jblok /= 0) then
         call wrtout([ab_out,std_out], ' Penalized local magnetic moments induced by electric field ')
         call wrtout([ab_out,std_out], '  atom1  dir  efld.dir         Real              Imag')
         do irow=1, ndim
           do icol=(natom+1)*3+1, (natom+2)*3
             write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
           & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
           & real(barmmom(irow,icol)), aimag(barmmom(irow,icol))
             call wrtout([ab_out,std_out], msg)
           end do
         end do
         call wrtout([ab_out,std_out], '   ')
       end if
       if (lblok /= 0) then
         call wrtout([ab_out,std_out], ' Penalized local magnetic moments induced by macroscopic Zeeman field ')
         call wrtout([ab_out,std_out], '  atom1  dir  Bfld.dir         Real              Imag')
         do irow=1, ndim
           do icol=(natom+4)*3+1, (natom+5)*3
             write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
           & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
           & real(barmmom(irow,icol)), aimag(barmmom(irow,icol))
             call wrtout([ab_out,std_out], msg)
           end do
         end do
         call wrtout([ab_out,std_out], '   ')
       end if
     end if
   end if
 end if 
  
!Store the FS and RS flavors on the DDB array
 do ipert2=1,natom+2
   do idir2=1,3
     icol=idir2+(ipert2-1)*3
     indexat2(icol)=ipert2
     indexdir2(icol)=idir2

     ipert1_red= 0
     do iat1= mpatpol(1), mpatpol(2)
       ipert1= natom + 11 + iat1
       ipert1_red= ipert1_red + 1
       idir1_red= 0
       do idir1= 1, 3
         if (mpdir(idir1)==0) cycle
         idir1_red= idir1_red + 1
         irow=idir1_red+(ipert1_red-1)*nmdir
         indexat1(irow)=iat1
         indexdir1(irow)=idir1
         index= idir1 + 3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))
         jndex= idir2 + 3*((ipert2-1)+mpert*((idir1-1)+3*(ipert1-1)))
         
         zblok= 0
         if (iblok /=0 .and. ipert2 <= natom) zblok= iblok
         if (jblok /=0 .and. ipert2 == natom+2) zblok= jblok
         if (lblok /=0 .and. ipert2 == natom+5) zblok= lblok
         if (zblok /= 0) then
           ddb%val_fs(1,index,zblok)=real(zfield(irow,icol))
           ddb%val_fs(2,index,zblok)=aimag(zfield(irow,icol))
           ddb%val_fs(1,jndex,zblok)=real(zfield_tr(icol,irow))
           ddb%val_fs(2,jndex,zblok)=aimag(zfield_tr(icol,irow))

           ddb%val_rs(1,index,zblok)=real(mmom(irow,icol))
           ddb%val_rs(2,index,zblok)=aimag(mmom(irow,icol))
           ddb%val_rs(1,jndex,zblok)=real(mmom_tr(icol,irow))
           ddb%val_rs(2,jndex,zblok)=aimag(mmom_tr(icol,irow))
         end if

       end do
     end do
   end do
 end do

 end subroutine magmom
!!***

!!****f* m_ddb_magpen/mp_d2etot
!! NAME
!! mp_d2etot
!!
!! FUNCTION
!! Calculate the different magnetic flavors (see mpopt below) of the second-
!! order derivatives of total energy
!!
!! INPUTS
!! barmagsus(ndim,ndim)= Penalized spin-sussceptibility tensor (\bar{\chi})
!! (equal to barmom^{\dagger} in the nondissipative regime)
!! ddb= the ddb object
!! dissip= if 0 a nondissipative regime is assumed
!!         if 1 a dissipative regime is assumed with a finite \eta introduced at the interpolation in omega regime
!! iblok= index of the IFCs block
!! invhmat(ndim,ndim)= (I-\alpha\bar{\chi})^-1 matrix
!! magsus(ndim,ndim)= Spin-sussceptibility tensor
!! magpen= magnetic penalty amplitude 
!! mpert =maximum number of ipert
!! mpopt = 1 calculate the frozen-spin second-order quantities
!!         2 calculate the relaxed-spin second-order quantities 
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! ndim= number of local magnetic degres of freedom
!! qphon= momentum wave-vector
!! xred(3,natom)= reduced atomic coordinates
!! zfield(ndim,(natom+2)*3)= First-order induced Zeeman fields
!! zfield_tr(natom+2)*3,ndim)= Linear-responses to external Zeeman fields
!! (equal to zfield^{\dagger} in the nondissipative regime)
!!
!! OUTPUT
!! ddb%val_fs(2,msize,nblok)= second-order derivatives at fixed spin.
!! ddb%val_rs(2,msize,nblok)= second-order derivatives at relaxed spin.
!!
!! SOURCE

 subroutine mp_d2etot(barmagsus,ddb,dissip,&
& iblok,invhmat,magsus,magpen,mpert,mpopt,&
& natom,nblok,ndim,qphon,xred,zfield,zfield_tr)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,dissip,mpert,mpopt,natom,nblok,ndim
 real(dp),intent(in) :: magpen
!arrays
 type(ddb_type),intent(inout) :: ddb
 real(dp),intent(in) :: qphon(3),xred(3,natom)
 complex(dpc),intent(in) :: barmagsus(ndim,ndim)
 complex(dpc),intent(in) :: invhmat(ndim,ndim)
 complex(dpc),intent(in) :: magsus(ndim,ndim)
 complex(dpc),intent(in) :: zfield(ndim,(natom+5)*3)
 complex(dpc),intent(in) :: zfield_tr((natom+5)*3,ndim)
!Local variables -------------------------
!scalars
 integer :: idir1,idir2,ipert1,ipert2,index,irow,icol
 complex(dpc) :: val_ps,val_fs, val_rs
!arrays
 
! *********************************************************************

 do ipert2= 1, natom+5
   do idir2= 1, 3
     icol= (ipert2-1)*3 + idir2
     do ipert1= 1, natom+5
       do idir1= 1, 3
         irow= (ipert1-1)*3 + idir1
         index= idir1 + 3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))

         !Extract the penalized second-order derivatives 
         val_ps= cmplx(ddb%val(1,index,iblok),ddb%val(2,index,iblok),16)
       
         !Calculate the fixed-spin flavor
         val_fs= val_ps + &
       & sum( zfield_tr(irow,:) * matmul( barmagsus,zfield(:,icol) ) ) 
         ddb%val_fs(1,index,iblok)= real(val_fs)
         ddb%val_fs(2,index,iblok)= aimag(val_fs)

         if (mpopt==2) then
           !Calculate the relaxed-spin flavor
           val_rs= val_fs - &
         & sum( zfield_tr(irow,:) * matmul( magsus,zfield(:,icol) ) ) 
           ddb%val_rs(1,index,iblok)= real(val_rs)
           ddb%val_rs(2,index,iblok)= aimag(val_rs)
         end if

       end do
     end do
   end do
 end do 

!TODO:This should be applied to the magnetic variables instead
! !Adopt the same phase convention as for the local Zeeman perturbation
! do ipert2= 1, natom
!   do idir2= 1, 3
!     icol=( ipert2-1)*3 + idir2
!     do ipert1= 1, natom
!       do idir1= 1, 3
!         irow=( ipert1-1)*3 + idir1
!         fmifc_sf(irow,icol)=fmifc(irow,icol)* &
!       & exp(two_pi*(0.d0,1.d0)* dot_product(qphon,xred(:,ipert2)-xred(:,ipert1)))
!       end do
!     end do
!   end do
! end do 

 end subroutine mp_d2etot
!!***

!!****f* m_ddb_magpen/mp_d2etot_print
!! NAME
!! mp_d2etot_print
!!
!! FUNCTION
!! Write on output file the fixed- and relaxed-spin susceptibilities
!! Only macroscopic quantities are printed if prtvol=1
!!
!! INPUTS
!! blkval= 2nd-order susceptibilities matrix 
!! kblok= index of the current block
!! opt= 1 write the frozen-spin second-order quantities
!!      2 write the relaxed-spin second-order quantities 
!! omega= frequency of the perturbation
!! qphon= momentum wave-vector
!! ucvol= unit-cell volume
!!
!! OUTPUT
!!
!! SOURCE

 subroutine mp_d2etot_print(ddb,blkval,kblok,mpert,natom,nblok,opt,omega,prtvol,qeq0,qphnrm,qphon,ucvol)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(in) :: ddb
 integer,intent(in) :: kblok,mpert,natom,nblok,opt,prtvol
 logical,intent(in) :: qeq0
 real(dp),intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: omega(3)
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 real(dp),intent(inout) :: qphnrm(3),qphon(3,3)

!Local variables -------------------------
!scalars
 integer :: iblok,idir1,idir2,ipert1,ipert2,index,irow,icol
 integer :: rftyp
 character(len=1000) :: msg
!arrays
 integer :: rfelfd(4),rfmagn(4),rfphon(4),rfstrs(4)
 real(dp) :: val(2)
 character(len=1) :: cart(3)=(/'x','y','z'/)
 
! *********************************************************************

 rfelfd(:)=0
 rfphon(:)=0
 rfstrs(:)=0
 rfmagn(:)=0
 rftyp = 1

 !IFCs
 if (prtvol >1) then
   rfphon(1:2)=1
   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega)
   if (iblok/=0.and.iblok==kblok) then
     if (opt==1) then
       call wrtout([ab_out,std_out], ' Frozen-spin interatomic force constants')
     else if (opt==2) then
       call wrtout([ab_out,std_out], ' Relaxed-spin interatomic force constants')
     end if
     call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
     do ipert1= 1, natom
       do idir1= 1, 3
         irow=( ipert1-1)*3 + idir1
         do ipert2= 1, natom
           do idir2= 1, 3
             icol=( ipert2-1)*3 + idir2
             val(:)=blkval(:,idir1,ipert1,idir2,ipert2,kblok)
             write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)') &
           & ipert1, cart(idir1), ipert2, cart(idir2), val(1), val(2)
             call wrtout([ab_out,std_out], msg)
           end do
         end do
         call wrtout([ab_out,std_out], ' ')
       end do
     end do
   end if
 end if

 if (qeq0) then

   !Born charges 
   rfphon(1:2)=1
   rfelfd(1:2)=2
   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega)
   if (iblok/=0.and.iblok==kblok) then
     if (opt==1) then
       call wrtout([ab_out,std_out], ' Frozen-spin Born effective charges')
     else if (opt==2) then
       call wrtout([ab_out,std_out], ' Relaxed-spin Born effective charges')
     end if
     call wrtout([ab_out,std_out], ' efld.dir   atom   dir        Real              Imag')
     ipert1= natom + 2
     do idir1= 1, 3
       do ipert2= 1, natom
         do idir2= 1, 3
           val(:)=blkval(:,idir1,ipert1,idir2,ipert2,kblok)
           write(msg,'(3x,a2,7x,i3,4x,a2,2x,2es18.9)') &
         & cart(idir1), ipert2, cart(idir2), val(1), val(2)
           call wrtout([ab_out,std_out], msg)
         end do
       end do
       call wrtout([ab_out,std_out], ' ')
     end do
   end if

   !Dielectric tensor
   iblok=0
   rfphon(:)=0
   rfelfd(1:2)=2
   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega)
   if (iblok/=0.and.iblok==kblok) then
     if (opt==1) then
       call wrtout([ab_out,std_out], ' Frozen-spin clamped-ion dielectric tensor')
     else if (opt==2) then
       call wrtout([ab_out,std_out], ' Relaxed-spin clamped-ion dielectric tensor')
     end if
     call wrtout([ab_out,std_out], '  dir  dir        Real              Imag')
     ipert1= ddb%natom + 2
     ipert2= ddb%natom + 2
     do idir2= 1, 3
       do idir1= 1, 3
         val(:)=blkval(:,idir1,ipert1,idir2,ipert2,kblok)
         write(msg,'(2x,a2,3x,a2,2x,2es18.9)' ) cart(idir1), cart(idir2), &
       & val(1), val(2)
         call wrtout([ab_out,std_out], msg)
       end do
       call wrtout([ab_out,std_out], ' ')
     end do
   end if

   !Magnetoelectric susceptibility
   iblok=0
   rfphon(:)=0
   rfelfd(1)=0
   rfelfd(2)=2
   rfmagn(1)=1
   rfmagn(2)=0
   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega)
   if (iblok/=0.and.iblok==kblok) then
     if (opt==1) then
       call wrtout([ab_out,std_out], ' Frozen-spin clamped-ion magnetoelectric susceptibility')
     else if (opt==2) then
       call wrtout([ab_out,std_out], ' Relaxed-spin clamped-ion magnetoelectric susceptibility')
     end if
     call wrtout([ab_out,std_out], ' M-dir E-dir        Real              Imag')
     ipert1= ddb%natom + 5
     ipert2= ddb%natom + 2
     do idir2= 1, 3
       do idir1= 1, 3
         val(:)=blkval(:,idir1,ipert1,idir2,ipert2,kblok)/ucvol
         write(msg,'(2x,a2,3x,a2,2x,2es18.9)' ) cart(idir1), cart(idir2), &
       & val(1), val(2)
         call wrtout([ab_out,std_out], msg)
       end do
       call wrtout([ab_out,std_out], ' ')
     end do
   end if

 end if

 !Magnetic susceptibility
 iblok=0
 rfphon(:)=0
 rfelfd(:)=0
 rfmagn(1)=1
 rfmagn(2)=1
 call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega)
 if (iblok/=0.and.iblok==kblok) then
   if (opt==1) then
     call wrtout([ab_out,std_out], ' Frozen-spin clamped-ion magnetic susceptibility')
   else if (opt==2) then
     call wrtout([ab_out,std_out], ' Relaxed-spin clamped-ion magnetic susceptibility')
   end if
   call wrtout([ab_out,std_out], '  dir  dir        Real              Imag')
   ipert1= ddb%natom + 5
   ipert2= ddb%natom + 5
   do idir2= 1, 3
     do idir1= 1, 3
       val(:)=blkval(:,idir1,ipert1,idir2,ipert2,kblok)/ucvol
       write(msg,'(2x,a2,3x,a2,2x,2es18.9)' ) cart(idir1), cart(idir2), &
     & val(1), val(2)
       call wrtout([ab_out,std_out], msg)
     end do
     call wrtout([ab_out,std_out], ' ')
   end do
 end if

 !Magnetic Born effective charges
 iblok=0
 rfphon(1)=1
 rfelfd(:)=0
 rfmagn(1)=0
 rfmagn(2)=1
 call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega)
 if (iblok/=0.and.iblok==kblok) then
   if (opt==1) then
     call wrtout([ab_out,std_out], ' Frozen-spin magnetic Born effective charges')
   else if (opt==2) then
     call wrtout([ab_out,std_out], ' Relaxed-spin magnetic Born effective charges')
   end if
   call wrtout([ab_out,std_out], ' atom   dir     B-dir        Real              Imag')

   ipert2= ddb%natom + 5
   do ipert1= 1, natom
     do idir1= 1, 3
       do idir2= 1, 3
         val(:)=blkval(:,idir1,ipert1,idir2,ipert2,kblok)
         write(msg,'(i3,4x,a2,7x,a2,2x,2es18.9)') &
       & ipert1, cart(idir1), cart(idir2), val(1), val(2)
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], ' ')
   end do
 end if

 end subroutine mp_d2etot_print
!!***

!!****f* m_ddb_magpen/asrw0
!! NAME
!! asrw0
!!
!! FUNCTION
!! Impose the Acoustic Sum Rule from the w=0 IFCs
!!
!! INPUTS
!! natom= namber of atoms
!! option= if 0, this is the w=0 case, calculate delta_asrw0
!          if 1, use the previously calculated delta_asrw0
!!
!! OUTPUT
!! ifcmat(3*natom,3*natom)= IFC matrix after ASR has been applied.
!! delta_asrw0(3*natom,3)= Amount to remove in order to enforce ASR.
!!
!! SOURCE

 subroutine asrw0(delta_asrw0,ifcmat,natom,option)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,option
!arrays
 real(dp),intent(inout) :: delta_asrw0(3*natom,3)
 complex(dpc),intent(inout) :: ifcmat(3*natom,3*natom)
!Local variables -------------------------
!scalars
 integer :: icol,idir1,idir2,ipert1,ipert2,irow
 real(dp) :: fac
 character(len=1000) :: msg
!arrays

! *********************************************************************

!Calculate the ASR correction 
 if (option == 0) then
   delta_asrw0= zero
   do idir1= 1, 3
     do ipert1= 1, natom
       irow= (ipert1-1)*3 + idir1
       do idir2= 1, 3
         do ipert2= 1, natom
           icol= (ipert2-1)*3 + idir2
           delta_asrw0(irow,idir2)=delta_asrw0(irow,idir2) + &
         & real(ifcmat(irow,icol))
         end do
       end do
     end do
   end do
 end if
           
!Apply the ASR
 do idir1= 1, 3
   do ipert1= 1, natom
     irow= (ipert1-1)*3 + idir1
     do idir2= 1, 3
       icol= (ipert1-1)*3 + idir2
       ifcmat(irow,icol)= ifcmat(irow,icol) - delta_asrw0(irow,idir2)
     end do
   end do
 end do
  
 end subroutine asrw0
!!***

!!****f* m_ddb_magpen/berrycurv_ss
!! NAME
!! berrycurv_ss
!!
!! FUNCTION
!! Calculate the Berry curvature of the inverse magnetic susceptibility
!! (this is equivalent to the G^(ss) matrix of S.Ren et al.)
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert*3*mpert,nblok)=  Third-order derivative matrices
!!  In our case, the nblok is restricted to iblok
!! iblok= index of the current block
!! invbarmagsus(ndim,ndim)= Inverse of the penalized spin-sussceptibility tensor
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! nmdir= number of directions along which the magnetic penalty was applied
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!!
!! OUTPUT
!! bc_ss(ndim,ndim)= Berry-curvature of the inverse local-spin susceptibility
!!
!! SOURCE

 subroutine berrycurv_ss(bc_barmagsus,bc_ss,blkval,iblok,invbarmagsus,&
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,omega,omegaflag,prtvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok,ndim,nmdir,omegaflag,prtvol
 real(dp),intent(in) :: omega
!arrays
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(inout) :: blkval(2,3,mpert,3,mpert,3,mpert,nblok)
 complex(dpc),intent(in) :: invbarmagsus(ndim,ndim)
 complex(dpc),intent(out) :: bc_ss(ndim,ndim)
 complex(dpc),intent(out) :: bc_barmagsus(ndim,ndim)
!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,idir3,info,ipert1,ipert2,ipert3,irow,lwork
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
 real(dp) :: fac
 complex(dpc), parameter :: ione=(0.d0,1.d0)
 character(len=1000) :: msg
!arrays
 complex(dpc) :: idty(ndim,ndim)
 integer(dp) :: indexat(ndim),indexdir(ndim)
 integer, allocatable :: ipiv(:)
 complex(dpc),allocatable :: work(:),work1(:,:),work2(:,:)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

 !Extract the Berry-curvature of the penalized susceptibility
 idty=(zero,zero)
 ipert3= natom + 9
 idir3= 1
 ipert2_red= 0
 do iat2= mpatpol(1), mpatpol(2)
   ipert2= natom + 11 + iat2
   ipert2_red= ipert2_red + 1
   idir2_red= 0
   do idir2= 1, 3
     if (mpdir(idir2)==0) cycle
     idir2_red= idir2_red + 1
     icol=idir2_red+(ipert2_red-1)*nmdir
     indexat(icol)=iat2
     indexdir(icol)=idir2
     idty(icol,icol)=(one,zero)
     ipert1_red=0
     do iat1= mpatpol(1), mpatpol(2)
       ipert1= natom + 11 + iat1
       ipert1_red= ipert1_red + 1
       idir1_red= 0
       do idir1= 1, 3
         if (mpdir(idir1)==0) cycle
         idir1_red=idir1_red+1
         irow=idir1_red+(ipert1_red-1)*nmdir

         bc_barmagsus(irow,icol)= -one* &
       & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok), &
       & blkval(2,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok),16)

       end do
     end do
   end do
 end do

!Calculate the Berry-curvature of the inverse magnetic susceptibility
 bc_ss=-matmul(invbarmagsus,matmul(bc_barmagsus,invbarmagsus)) 

! fac=2.714943600699**2/four !TMP
 fac=one/four !TMP
 open(10,file='g_ss.txt')
   do irow=1, ndim
     write(10,*) -ione*bc_ss(irow,1:ndim)*fac
   end do 
 close(10)


 call wrtout([ab_out,std_out], ' Berry curvature of the inverse spin susceptibility ')
 call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
 do irow=1, ndim
   do icol=1, ndim
     write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
   & indexat(irow), cart(indexdir(irow)), indexat(icol), cart(indexdir(icol)), &
   & real(bc_ss(irow,icol)), aimag(bc_ss(irow,icol))
     call wrtout([ab_out,std_out], msg)
   end do
 end do
 call wrtout([ab_out,std_out], '   ')

!For linear interpolation of Hessians substitute FM Berry curvature 
!into the ddb_lw object.
 if (omegaflag==2.and.abs(omega)<tol12) then
   ipert2_red= 0
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
           blkval(1,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok)= &
         & real(bc_ss(irow,icol))
           blkval(2,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok)= &
         & aimag(bc_ss(irow,icol))
         end do
       end do
     end do
   end do
 end if

 end subroutine berrycurv_ss
!!***

!!****f* m_ddb_magpen/berrycurv_sp
!! NAME
!! berrycurv_sp
!!
!! FUNCTION
!! Calculate the Berry curvature of the spin-phonon Hessian
!! (equivalent to the Berry curvature of the induced Zeeman fields at constrained
!! magnetic moments)
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert*3*mpert,nblok)=  Third-order derivative matrices
!!  In our case, the nblok is restricted to iblok
!! iblok= index of the current block
!! invbarmagsus(ndim,ndim)= Inverse of the penalized spin-sussceptibility tensor
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! nmdir= number of directions along which the magnetic penalty was applied
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!!
!! OUTPUT
!!
!! SOURCE

 subroutine berrycurv_sp(barmmom,bc_sp,bc_ss,blkval,iblok,invbarmagsus,&
& jblok,mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,omega,omegaflag,prtvol,qphon,xred)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,jblok,mpert,natom,nblok,ndim,nmdir,omegaflag,prtvol
 real(dp),intent(in) :: omega
!arrays
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(inout) :: blkval(2,3,mpert,3,mpert,3,mpert,nblok)
 real(dp),intent(in) :: qphon(3),xred(3,natom)
 complex(dpc),intent(in) :: barmmom(ndim,(natom+2)*3)
 complex(dpc),intent(in) :: bc_ss(ndim,ndim)
 complex(dpc),intent(in) :: invbarmagsus(ndim,ndim)
 complex(dpc),intent(out) :: bc_sp(ndim,(natom+2)*3)
!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,idir3,ipert1,ipert2,ipert3,irow
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
 real(dp) :: fac,re,im
 complex(dpc), parameter :: ione=(0.d0,1.d0)
 character(len=1000) :: msg
!arrays
 integer(dp) :: indexat1(ndim),indexdir1(ndim)
 integer(dp) :: indexat2((natom+2)*3),indexdir2((natom+2)*3)
 complex(dpc) :: bc_barsp(ndim,(natom+2)*3)
 complex(dpc) :: bc_ps((natom+2)*3,ndim)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************


!Extract the berry curvature of the penalized moments
 bc_barsp=(zero,zero)
 ipert3= natom + 9
 idir3= 1
 do ipert2=1,natom+2
   do idir2=1,3
     icol=idir2+(ipert2-1)*3
     indexat2(icol)=ipert2
     indexdir2(icol)=idir2

     ipert1_red= 0
     do iat1= mpatpol(1), mpatpol(2)
       ipert1= natom + 11 + iat1
       ipert1_red= ipert1_red + 1
       idir1_red= 0
       do idir1= 1, 3
         if (mpdir(idir1)==0) cycle
         idir1_red= idir1_red + 1
         irow=idir1_red+(ipert1_red-1)*nmdir
         indexat1(irow)=iat1
         indexdir1(irow)=idir1
         
         if (iblok /=0 .and. ipert2 <= natom) then
           bc_barsp(irow,icol)= -one* &
         & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok), &
         & blkval(2,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok),16)
         else if (jblok /=0 .and. ipert2 == natom+2) then
           bc_barsp(irow,icol)= -one* &
         & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,idir3,ipert3,jblok), &
         & blkval(2,idir1,ipert1,idir2,ipert2,idir3,ipert3,jblok),16)
         end if

       end do
     end do
   end do
 end do

 !Calculate the Berry curvature of the induced Zeeman fields
 bc_sp= -matmul(bc_ss,barmmom) - matmul(invbarmagsus,bc_barsp)

 do irow=1,ndim
   do iat1= 1, natom
     do idir1= 1, 3
       icol= (iat1-1)*3 + idir1
       !MR: Caution, this conjg might be wrong in presence of dissipation
       bc_ps(icol,irow)=conjg(bc_sp(irow,icol)*exp(two_pi*(0.d0,1.d0)* dot_product(qphon,xred(:,iat1))))
     end do
   end do
 end do 

! fac=2.714943600699/two/0.52917 !TMP
 fac=one/two/0.52917 !TMP
 open(10,file='g_ps.txt')
 do irow=1,natom*3
   write(10,*) (0.d0,-1.d0)*bc_ps(irow,1:ndim)*fac
 end do 
 close(10)

 if (iblok /= 0) then
   call wrtout([ab_out,std_out], ' Berry curvature of the Zeeman fields induced by atomic displacements (at constrained magnetic moments)')
   call wrtout([ab_out,std_out], '  atom1  dir  efld.dir         Real              Imag')
   do irow=1, ndim
     do icol=1, natom*3
       write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)' ) &
     & indexat1(irow), cart(indexdir1(irow)), indexat2(icol), cart(indexdir2(icol)), &
     & real(bc_sp(irow,icol)), aimag(bc_sp(irow,icol))
       call wrtout([ab_out,std_out], msg)
     end do
   end do
   call wrtout([ab_out,std_out], '   ')
 end if
 if (jblok /= 0) then
   call wrtout([ab_out,std_out], ' Berry curvature of the Zeeman fields induced by electric field (at constrained magnetic moments)')
   call wrtout([ab_out,std_out], '  atom1  dir  efld.dir         Real              Imag')
   do irow=1, ndim
     do icol=(natom+1)*3+1, (natom+2)*3
       write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
     & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
     & real(bc_sp(irow,icol)), aimag(bc_sp(irow,icol))
       call wrtout([ab_out,std_out], msg)
     end do
   end do
   call wrtout([ab_out,std_out], '   ')
 end if

!For linear interpolation of Hessians substitute sp  Berry curvature 
!into the ddb_lw object.
 if (omegaflag==2.and.abs(omega)<tol12) then
   do ipert2=1,natom+2
     do idir2=1,3
       icol=idir2+(ipert2-1)*3
       indexat2(icol)=ipert2
       indexdir2(icol)=idir2
  
       ipert1_red= 0
       do iat1= mpatpol(1), mpatpol(2)
         ipert1= natom + 11 + iat1
         ipert1_red= ipert1_red + 1
         idir1_red= 0
         do idir1= 1, 3
           if (mpdir(idir1)==0) cycle
           idir1_red= idir1_red + 1
           irow=idir1_red+(ipert1_red-1)*nmdir
           
           if (iblok /=0 .and. ipert2 <= natom) then
           blkval(1,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok)= &
         & real(bc_sp(irow,icol))
           blkval(2,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok)= &
         & aimag(bc_sp(irow,icol))
           blkval(1,idir2,ipert2,idir1,ipert1,idir3,ipert3,iblok)= &
         & real(bc_sp(irow,icol))
           blkval(2,idir2,ipert2,idir1,ipert1,idir3,ipert3,iblok)= &
         & -aimag(bc_sp(irow,icol))
           end if
  
         end do
       end do
     end do
   end do
 end if

 end subroutine berrycurv_sp
!!***

!!****f* m_ddb_magpen/berrycurv_pp
!! NAME
!! berrycurv_pp
!!
!! FUNCTION
!! Calculate the Berry curvature of the phonon-phonon Hessian
!! (at constrained magnetic moments)
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert*3*mpert,nblok)=  Third-order derivative matrices
!!  In our case, the nblok is restricted to iblok
!! iblok= index of the current block
!! invbarmagsus(ndim,ndim)= Inverse of the penalized spin-sussceptibility tensor
!! mpatpol(2) = Atoms on which the magnetic penalty has been applied
!! mpdir(3) = Directions along which the spin-degrees of freedom have been stiffened 
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! nblok= number of blocks in the DDB
!! nmdir= number of directions along which the magnetic penalty was applied
!! ndim= dimension of the square susceptibilities 
!! prtvol= control the volume of information written on output
!!
!! OUTPUT
!!
!! SOURCE

 subroutine berrycurv_pp(barmagsus,bc_barmagsus,bc_sp,blkval,iblok,&
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,omega,omegaflag,prtvol,qphon,xred,zfield)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok,ndim,nmdir,omegaflag,prtvol
 real(dp),intent(in) :: omega
!arrays
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(inout) :: blkval(2,3,mpert,3,mpert,3,mpert,nblok)
 real(dp),intent(in) :: qphon(3),xred(3,natom)
 complex(dpc),intent(in) :: barmagsus(ndim,ndim)
 complex(dpc),intent(in) :: bc_barmagsus(ndim,ndim)
 complex(dpc),intent(in) :: bc_sp(ndim,(natom+2)*3)
 complex(dpc),intent(in) :: zfield(ndim,(natom+2)*3)
!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,idir3,ipert1,ipert2,ipert3,irow
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
 real(dp) :: fac
 complex(dpc), parameter :: ione=(0.d0,1.d0)
 character(len=1000) :: msg
!arrays
 integer(dp) :: indexat1(ndim),indexdir1(ndim)
 integer(dp) :: indexat2((natom+2)*3),indexdir2((natom+2)*3)
 complex(dpc) :: bc_barpp(natom*3,natom*3), bc_pp(natom*3,natom*3) 
 complex(dpc) :: bc_pp_sf(natom*3,natom*3), term(natom*3,natom*3,3)
 complex(dpc) :: ifc_bc_sp(ndim,natom*3), ifc_zfield(ndim,natom*3)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

!Extract the frequency derivative of the penalized IFCs 
 ipert3= natom + 9
 idir3= 1
 do ipert2= 1, natom
   do idir2= 1, 3
     icol=( ipert2-1)*3 + idir2
     do ipert1= 1, natom
       do idir1= 1, 3
         irow=( ipert1-1)*3 + idir1
         bc_barpp(irow,icol)= &
       & cmplx(blkval(1,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok), &
       &       blkval(2,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok),16)
       end do
     end do
   end do
 end do 

!Calculate the different terms entering the Berry curvature
 ifc_zfield(:,:)=zfield(:,1:natom*3)
 ifc_bc_sp(:,:)=bc_sp(:,1:natom*3)
 term(:,:,1)=matmul(transpose(conjg(ifc_bc_sp)),matmul(barmagsus,ifc_zfield))
 term(:,:,2)=matmul(transpose(conjg(ifc_zfield)),matmul(bc_barmagsus,ifc_zfield))
 term(:,:,3)=matmul(transpose(conjg(ifc_zfield)),matmul(barmagsus,ifc_bc_sp))

 bc_pp(:,:)= bc_barpp(:,:) + term(:,:,1) + term(:,:,2) + term(:,:,3)

 do ipert2= 1, natom
   do idir2= 1, 3
     icol=( ipert2-1)*3 + idir2
     do ipert1= 1, natom
       do idir1= 1, 3
         irow=( ipert1-1)*3 + idir1
         bc_pp_sf(irow,icol)=bc_pp(irow,icol)* &
       & exp(two_pi*(0.d0,1.d0)* dot_product(qphon,xred(:,ipert2)-xred(:,ipert1)))
       end do
     end do
   end do
 end do 

 fac=one/(0.52917)**2
 open(10,file='g_pp.txt')
 do irow=1,natom*3
   write(10,*) -ione*bc_pp_sf(irow,:)*fac
 end do 
 close(10)

 !Write the results
 call wrtout([ab_out,std_out], ' Berry curvature of interatomic force constants (at constrained magnetic moments)')
 call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
 do ipert1= 1, natom
   do idir1= 1, 3
     irow=( ipert1-1)*3 + idir1
     do ipert2= 1, natom
       do idir2= 1, 3
         icol=( ipert2-1)*3 + idir2
         write(msg,'(2(i4,4x,a2,2x),2x,2es18.9)') &
       & ipert1, cart(idir1), ipert2, cart(idir2), &
       & real(bc_pp(irow,icol)), aimag(bc_pp(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
   end do
 end do 

!For linear interpolation of Hessians substitute pp Berry curvature 
!into the ddb_lw object.
 if (omegaflag==2.and.abs(omega)<tol12) then
   do ipert2= 1, natom
     do idir2= 1, 3
       icol=( ipert2-1)*3 + idir2
       do ipert1= 1, natom
         do idir1= 1, 3
           irow=( ipert1-1)*3 + idir1
           blkval(1,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok)= &
         & real(bc_pp(irow,icol))
           blkval(2,idir1,ipert1,idir2,ipert2,idir3,ipert3,iblok)= &
         & aimag(bc_pp(irow,icol))
         end do
       end do
     end do
   end do 

 end if

 end subroutine berrycurv_pp

!!****f* m_ddb_omega_interpol/ddb_omega_interpol
!! NAME
!! ddb_omega_interpol
!!
!! FUNCTION
!! Interpolate over frequency the secon-order derivatives calculated with 
!! the magnetic penalty and latter on convert them into physically relevant 
!! quantities.
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
& dissip,eta,eta_phongreen,outfilename_radix,magpen,mpatpol,mpdir,mpert,mpopt,natom, &
& nomega,ntypat,omegaflag,omegamax,omegamin,prtvol,rftyp,typat,ucvol,xred)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: dissip,mpert,mpopt,natom,nomega,ntypat,omegaflag,prtvol,rftyp
 real(dp),intent(in) :: eta,eta_phongreen,magpen,omegamax,omegamin,ucvol
 character(len=*),intent(in) :: outfilename_radix
!arrays
 type(ddb_type),intent(inout) :: ddb,ddb_lw
 integer,intent(in) :: mpatpol(2),mpdir(3),typat(natom)
 real(dp),intent(in) :: amu(ntypat)
 real(dp), intent(inout) :: delta_asrw0(3*natom,3), delta_asrw0_fm(3*natom,3)
 real(dp),intent(in) :: xred(3,natom)

!Local variables -------------------------
!scalars
 integer :: diel_unit,i,iblok,ifound,ii,imode,ipert1,ipert2,iw,j,jblok,jw,kblok,lblok,mmom_unit,mmspec_unit,nblok,ndim 
 integer :: nmat,nmdir,nwcalc,phon_unit,prtopt,spin_unit,zeff_unit,zeffspec_unit,zfield_unit
 real(dp) :: omegastp
 character(len=5000) :: msg,pfmt
 character(len=fnlen) :: diel_filename,spin_filename,mmom_filename,mmspec_filename
 character(len=fnlen) :: phon_filename,zeff_filename,zeffspec_filename,zfield_filename
 complex(dpc) :: cplxvar,cplx_weta
!arrays
 real(dp) :: qphnrm(3),qphon(3,3)
 real(dp), allocatable :: dint_barddb(:,:),int_barddb(:,:,:),omega(:),omegacalc(:)
 real(dp), allocatable :: w0hessian(:,:),w0berry(:,:)
 real(dp), allocatable :: eigvec(:,:,:,:,:),eigvec_fm(:,:,:,:,:),phfrq(:,:)
 real(dp), allocatable :: magphonspec(:),mode_magphonspec(:,:),mode_phonspec(:,:),phonspec(:)
 real(dp), allocatable :: coeffs(:,:)
 complex(dpc), allocatable :: barmagsus(:,:,:),invbarmagsus(:,:,:)
 complex(dpc), allocatable :: invmagsus(:,:,:), lm_magsus(:,:,:), magsus(:,:,:), invhmat(:,:)
 complex(dpc), allocatable :: barmmom(:,:),barmmom_tr(:,:),mmom(:,:,:), mmom_tr(:,:,:)
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
 complex(dpc) :: magbasis(4,4),work(4,4)
 complex(dpc),parameter :: ure=(1.d0,0.d0),uim=(0.d0,1.d0)
 real(dp) :: totnorm
 
! *********************************************************************

!TMP: Complete magnon basis
 magbasis(1,:)=0.5d0*(/ure,uim,ure,uim/)
 magbasis(2,:)=0.5d0*(/ure,-uim,ure,-uim/)
 magbasis(3,:)=0.5d0*(/ure,uim,-ure,-uim/)
 magbasis(4,:)=0.5d0*(/ure,-uim,-ure,uim/)

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
 ABI_MALLOC(barmagsus,(ndim,ndim,nomega))
 ABI_MALLOC(magsus,(ndim,ndim,nomega))
 ABI_MALLOC(lm_magsus,(ndim,ndim,nomega))
 ABI_MALLOC(ri_magelsus,(ndim,3))
 ABI_MALLOC(lm_magelsus,(ndim,3,nomega))
 ABI_MALLOC(invbarmagsus,(ndim,ndim,nomega))
 ABI_MALLOC(invmagsus,(ndim,ndim,nomega))
 ABI_MALLOC(invhmat,(ndim,ndim))
 ABI_MALLOC(barmmom,(ndim,(natom+2)*3))
 ABI_MALLOC(barmmom_tr,(ndim,(natom+2)*3))
 ABI_MALLOC(mmom,(ndim,(natom+2)*3,nomega))
 ABI_MALLOC(mmom_tr,((natom+2)*3,ndim,nomega))
 ABI_MALLOC(lm_zfield,(ndim,3,nomega))
 ABI_MALLOC(zfield,(ndim,(natom+2)*3,nomega))
 ABI_MALLOC(zfield_tr,((natom+2)*3,ndim))
 ABI_MALLOC(barepsilon,(3,3,nomega))
 ABI_MALLOC(epsilon,(3,3,nomega))
 ABI_MALLOC(macmagsus,(3,3,nomega))
 ABI_MALLOC(ifcmat,(3*natom,3*natom))
 ABI_MALLOC(ifcmat_fm,(3*natom,3*natom))
 ABI_MALLOC(dint_barddb,(2,ddb%msize))
 ABI_MALLOC(int_barddb,(2,ddb%msize,1))
 ABI_MALLOC(ri_mmom,(ndim,3,nomega))
 if (omegaflag == 3) then
   ABI_MALLOC(coeffs,(2,nwcalc))
 end if

!For linear interpolation detect the w=0 Hessians and Berry curvatures
 if (omegaflag == 2) then

   ABI_MALLOC(w0hessian,(2,ddb%msize))
   nblok= ddb%nblok
   ifound= 0
   do iblok= 1, nblok
     if (abs(ddb%omega(1,iblok)) < tol12) then
       w0hessian(:,:)= ddb%val(:,:,iblok)
       ifound= 1
     end if
   end do
   if (ifound==0) then
     write(msg, '(3a)' )' No omega=0 block with second-order derivatives', &
   & ' found in the DDB file. This is necessary if omegaflag=2 ',ch10
     ABI_ERROR(msg)
   end if

   ABI_MALLOC(w0berry,(2,ddb_lw%msize))
   nblok= ddb_lw%nblok
   ifound= 0
   do iblok= 1, nblok
     if (abs(ddb_lw%omega(1,iblok)) < tol12) then
       w0berry(:,:)= ddb_lw%val(:,:,iblok)
       ifound= 1
     end if
   end do
   if (ifound==0) then
     write(msg, '(3a)' )' No omega=0 block with third-order derivatives', &
   & ' found in the DDB file. This is necessary if omegaflag=2 ',ch10
     ABI_ERROR(msg)
   end if
 end if

!Loop over the frequency
 do iw=1,nomega
   omega(iw)=omegamin+omegastp*(iw-1)

   if (omegaflag==1) then
     do ii=1,ddb%msize
       if (all(ddb%flg(ii,:)==1)) then
         call POLINT(omegacalc,ddb%val(1,ii,:),nwcalc,omega(iw),int_barddb(1,ii,1),dint_barddb(1,ii)) 
         call POLINT(omegacalc,ddb%val(2,ii,:),nwcalc,omega(iw),int_barddb(2,ii,1),dint_barddb(2,ii)) 
       else if (count(ddb%flg(ii,:)==0)/=nwcalc) then
         write(msg,'(a,a,a)')&
         'ddb_omega_interpol detects differences between the DDB bloks for each frequency.',ch10,&
       & ' The interpolation has been stopped.' 
         ABI_ERROR(msg)
       end if
     end do 
   else if (omegaflag==2) then
     call lineal_omega_interp(w0hessian,w0berry,eta,ifcmat_fm, &
   & invmagsus(:,:,iw),magsus(:,:,iw),mpatpol,mpdir,mpert,ddb%msize, &
   & natom,ndim,nmdir,int_barddb,omega(iw),zfield(:,:,iw),zfield_tr)
   else if (omegaflag==3) then
     cplx_weta=cmplx(omega(iw),eta,16)
     do ii=1,ddb%msize
       if (all(ddb%flg(ii,:)==1)) then
         call polcoe(omegacalc,ddb%val(1,ii,:),nwcalc,coeffs(1,:))
         call polcoe(omegacalc,ddb%val(2,ii,:),nwcalc,coeffs(2,:))
         cplxvar=cmplx(zero,zero,16)
         do jw= 1, nwcalc
           cplxvar= cplxvar + cplx_weta**(jw-1)*cmplx(coeffs(1,jw),coeffs(2,jw),16)
         end do
         int_barddb(1,ii,1)=real(cplxvar)
         int_barddb(2,ii,1)=aimag(cplxvar) 
       else if (count(ddb%flg(ii,:)==0)/=nwcalc) then
         write(msg,'(a,a,a)')&
         'ddb_omega_interpol detects differences between the DDB bloks for each frequency.',ch10,&
       & ' The interpolation has been stopped.' 
         ABI_ERROR(msg)
       end if
     end do 
   end if

!To deactivate interpolation
!   int_barddb(:,:,1)=ddb%val(:,:,1)

   if (omegaflag==1.or.omegaflag==3) then
     !Calculate the local spin susceptibilities
     call local_spinsus(barmagsus(:,:,iw),ddb,1,invbarmagsus(:,:,iw),invmagsus(:,:,iw),&
   & invhmat,magpen,magsus(:,:,iw),mpatpol,mpdir,mpert,natom,1,ndim,nmdir,omega(iw),omegaflag,prtopt,prtvol)

     !Calculate the magnetic moments
     call magmom(barmmom,barmmom_tr,ddb,invbarmagsus(:,:,iw),invhmat,1,1,1,magpen,magsus(:,:,iw),mmom(:,:,iw),mmom_tr(:,:,iw),&
   & mpatpol,mpdir,mpert,natom,1,ndim,nmdir,omega(iw),omegaflag,prtopt,prtvol,qphon,xred,zfield(:,:,iw),zfield_tr)
  
!     !Calculate the dielectric susceptibility
!     call mp_diel(barepsilon(:,:,iw),barmagsus(:,:,iw),barmmom,barmmom_tr,&
!   & int_barddb,dissip,epsilon(:,:,iw),1,invhmat,magpen,magsus(:,:,iw),mpert,mpopt,&
!   & natom,1,ndim,prtopt,prtvol,ucvol,zfield(:,:,iw),zfield_tr)
!  
!     !Calculate the interatomic force constants
!     call mp_ifc(barmagsus(:,:,iw),barmmom,barmmom_tr,int_barddb,dissip,1,ifcmat,&
!   & ifcmat_fm,invhmat,magsus(:,:,iw),magpen,mpert,mpopt,&
!   & natom,1,ndim,omega(iw),omegaflag,prtopt,prtvol,qphon,xred,zfield(:,:,iw),zfield_tr)
!  
!     !Calculate the macroscopic magnetic susceptibility
!     call mp_macmagsus(barmagsus(:,:,iw),int_barddb,&
!   & dissip,1,invbarmagsus(:,:,iw),invhmat, macmagsus(:,:,iw),magpen,&
!   & magsus(:,:,iw),mpatpol,mpdir,mpert,mpopt,&
!   & natom,nblok,ndim,nmdir,prtopt,prtvol,ucvol)
   end if

   !Apply ASR
!   if (omega(1) < tol12) then
!     call asrw0(delta_asrw0,ifcmat,natom,0) 
!     call asrw0(delta_asrw0_fm,ifcmat_fm,natom,0) 
!   else 
!     call asrw0(delta_asrw0,ifcmat,natom,1) 
!     call asrw0(delta_asrw0_fm,ifcmat_fm,natom,1) 
!   end if

   !Calculate the phonon and magnon-phonon Green's functions and spectral functions
   call phonon_green(amu,eigvec,eigvec_fm,eta_phongreen,ifcmat,ifcmat_fm,invmagsus(:,:,iw),& 
 & magphongreen,magphonspec(iw),mode_magphonspec(:,iw),mode_phonspec(:,iw),natom,ndim,ntypat,omega(iw),&
 & phfrq(:,iw),phongreen,phongreen_fm,phonspec(iw),typat,zfield(:,:,iw),zfield_tr)

   if (omegaflag==1.or.omegaflag==3) then
     !Calculate the mode-resolved magnetic moments
     call mode_mmom(amu,eigvec,mmom(:,:,iw),mmomspec(:,iw),modemm(:,:,iw),mode_phonspec(:,iw),natom,ndim,ntypat,typat)

!     !Calculate the Born effective charges
!     call mp_zeff(barmagsus(:,:,iw),barmmom,barmmom_tr,int_barddb,&
!   & dissip,fmzeff,fmzeff_tr,1,invhmat,lm_epsilon(:,:,iw),magpen,magsus(:,:,iw),mpert,mpopt,&
!   & natom,1,ndim,phongreen,prtopt,prtvol,ucvol,zeff,zeff_tr,zfield(:,:,iw),zfield_tr)

     !Calculate the mode-resolved Born effective charges
!     call mode_zeff(amu,eigvec,mode_phonspec(:,iw),modezeff(:,:,iw),natom,ntypat,&
!   & typat,zeff(:,:,iw),zeffspec(:,iw))

     !Calculate here the phonons contribution to the spin susceptibility
     if (dissip==0) then
       lm_magsus(:,:,iw)=-matmul(mmom(:,1:natom*3,iw),matmul(phongreen,transpose(conjg(mmom(:,1:natom*3,iw)))))
     else if (dissip==1) then
       lm_magsus(:,:,iw)=-matmul(mmom(:,1:natom*3,iw),matmul(phongreen,mmom_tr(1:natom*3,:,iw)))
     end if

     !Calclate here the lattice-mediated magnetic moments induced by an electric field
     if (dissip==0) then
       lm_magelsus(:,:,iw)=-matmul(mmom(:,1:natom*3,iw),matmul(phongreen,transpose(conjg(zeff(:,:)))))
     else if (dissip==1) then
       lm_magelsus(:,:,iw)=-matmul(mmom(:,1:natom*3,iw),matmul(phongreen,zeff_tr(:,:)))
     end if

     !Alternative calculation
     genzeff_tr(1:natom*3,:)= fmzeff_tr(:,:)
     genzeff_tr(natom*3:natom*3+ndim,:)= zfield(:,(natom+2)*3-2:(natom+2)*3,iw)
     ri_genelsus(:,:,iw)=-matmul(magphongreen,genzeff_tr)

     !Another alternative

!TMP: FM eigvecs
     call me_altcalc(amu,eigvec_fm,lm_magsus(:,:,iw),lm_zfield(:,:,iw),phongreen_fm,magsus(:,:,iw),natom,ndim,ntypat,&
   & omega(iw),ri_mmom(:,:,iw),typat,fmzeff_tr,zfield(:,:,iw))

!TMP: SR eigvecs
!     call me_altcalc(amu,eigvec,lm_magsus(:,:,iw),lm_zfield(:,:,iw),phongreen_fm,magsus(:,:,iw),natom,ndim,ntypat,&
!   & omega(iw),ri_mmom(:,:,iw),typat,fmzeff_tr,zfield(:,:,iw))

     !Convert magnetic susceptibilities to the magnon basis
!     work(:,:)=magsus(:,:,iw)
!     magsus(:,:,iw)=matmul(transpose(conjg(magbasis)),matmul(work,magbasis))
!     work(:,:)=lm_magsus(:,:,iw)
!     lm_magsus(:,:,iw)=matmul(transpose(conjg(magbasis)),matmul(work,magbasis))
   end if
 end do

!Calculate the norm of the phonon spectral function
 totnorm=zero
 do imode=1, 3*natom
   write(msg,'(a,i3,a,es15.7)') 'Phonon mode: ', imode, & 
 & '. Norm of the spectral function: ', sum(mode_phonspec(imode,:))*omegastp
   call wrtout([ab_out,std_out],msg,'COLL')
   totnorm= totnorm + sum(mode_phonspec(imode,:))*omegastp
 end do
 write(msg,'(a,es15.7)') & 
 & ' Total norm of the spectral function: ', totnorm
   call wrtout([ab_out,std_out],msg,'COLL')

!Calculate the norm of the magnon phonon spectral function
 totnorm=zero
 do imode=1, 3*natom+ndim
   write(msg,'(a,i3,a,es15.7)') 'Magnon-Phonon mode: ', imode, & 
 & '. Norm of the spectral function: ', sum(mode_magphonspec(imode,:))*omegastp
   call wrtout([ab_out,std_out],msg,'COLL')
   totnorm= totnorm + sum(mode_magphonspec(imode,:))*omegastp
 end do
 write(msg,'(a,es15.7)') & 
 & ' Total norm of the spectral function: ', totnorm
   call wrtout([ab_out,std_out],msg,'COLL')

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

 write(spin_unit,*) '#  Real part of local spin susceptibility tensor (at. units)'
! write(spin_unit,*) '#  Real part of penalized local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(magsus(i,j,iw)),j=1,ndim),i=1,ndim)
! &  omega(iw), ((real(barmagsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

 write(spin_unit,*) ' '
 write(spin_unit,*) '#  Imaginary part of local spin susceptibility tensor (at. units)'
! write(spin_unit,*) '#  Imaginary part of penalized local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(magsus(i,j,iw)),j=1,ndim),i=1,ndim)
! &  omega(iw), ((aimag(barmagsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

 write(spin_unit,*) ' '
 write(spin_unit,*) '#  Real part of lattice-mediated local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(lm_magsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

 write(spin_unit,*) ' '
 write(spin_unit,*) '#  Imaginary part of lattice-mediated local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(lm_magsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

 write(spin_unit,*) ' '
 write(spin_unit,*) '#  Real part of relaxed-ion local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(magsus(i,j,iw)+lm_magsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

 write(spin_unit,*) ' '
 write(spin_unit,*) '#  Imaginary part of relaxed-ion local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(magsus(i,j,iw)+lm_magsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

 write(spin_unit,*) ' '
! write(spin_unit,*) '#  Real part of the inverse of the penalized local spin susceptibility tensor (at. units)'
 write(spin_unit,*) '#  Real part of the inverse of the local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X^{-1}_11     X^{-1}_12     ...     X^{-1}_21     X^{-1}_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
! &  omega(iw), ((real(invbarmagsus(i,j,iw)),j=1,ndim),i=1,ndim)
 &  omega(iw), ((real(invmagsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

 write(spin_unit,*) ' '
! write(spin_unit,*) '#  Imaginary part of the inverse of the penalized local spin susceptibility tensor (at. units)'
 write(spin_unit,*) '#  Imaginary part of the inverse of the local spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X^{-1}_11     X^{-1}_12     ...     X^{-1}_21     X^{-1}_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
! &  omega(iw), ((aimag(invbarmagsus(i,j,iw)),j=1,ndim),i=1,ndim)
 &  omega(iw), ((aimag(invmagsus(i,j,iw)),j=1,ndim),i=1,ndim)
    call wrtout(spin_unit,msg,'COLL')
 end do

 write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' ) 9 

 write(spin_unit,*) ' '
 write(spin_unit,*) '#  Real part of clamped-ion macroscopic spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(macmagsus(i,j,iw)),j=1,3),i=1,3)
    call wrtout(spin_unit,msg,'COLL')
 end do

 write(spin_unit,*) ' '
 write(spin_unit,*) '#  Imaginary part of clamped-ion macroscopic spin susceptibility tensor (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     X_11     X_12     ...     X_21     X_22     ...'
 call wrtout(spin_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(macmagsus(i,j,iw)),j=1,3),i=1,3)
    call wrtout(spin_unit,msg,'COLL')
 end do

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
 write(zfield_unit,*) '#  Real part of lattice-mediated local Zeeman fields induced by electric field (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     Z_11     Z_12      Z_13    ...     Z_21     Z_22     ...'
 call wrtout(zfield_unit,msg,'COLL')
 do iw=1,nomega
   write(msg,pfmt) omega(iw), ((real(lm_zfield(i,j,iw)),j=1,3),i=1,ndim)
   call wrtout(zfield_unit,msg,'COLL')
 end do

 write(zfield_unit,*) ' '
 write(zfield_unit,*) '#  Imag part of lattice-mediated local Zeeman fields induced by electric field (at. units)'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     Z_11     Z_12      Z_13    ...     Z_21     Z_22     ...'
 call wrtout(zfield_unit,msg,'COLL')
 do iw=1,nomega
   write(msg,pfmt) omega(iw), ((aimag(lm_zfield(i,j,iw)),j=1,3),i=1,ndim)
   call wrtout(zfield_unit,msg,'COLL')
 end do
 
 close (zfield_unit)

!Magnetic moments
 mmom_filename=trim(outfilename_radix)//"_MAGMOM"
 if (open_file(mmom_filename, msg, newunit=mmom_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(mmom_unit,*) '#'
 write(mmom_unit,*) '#  Magnetic moments calculated and interpolated by ANADDB'
 write(mmom_unit,*) '#'

 write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  ndim
 do imode= 1, 3*natom
   write(mmom_unit,*) ' '
   write(mmom_unit,'(a,i3)') '#  Real part of magnetic moments (at. units) induced by phonon mode:', imode
   write(msg,'(a,a,a)') ch10,&
 &           ' # At  hw     m_{mat_1,1}     m_{mat_1,2}     ...     m_{mat_2,1}     m_{mat_2,2}'
   call wrtout(mmom_unit,msg,'COLL')
   do iw=1,nomega
     write(msg,pfmt) &
   & omega(iw), (real(modemm(i,imode,iw)),i=1,ndim)
     call wrtout(mmom_unit,msg,'COLL')
   end do
   write(mmom_unit,*) ' '
   write(mmom_unit,'(a,i3)') '#  Imaginary part of magnetic moments (at. units) induced by phonon mode:', imode
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

 write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  ndim*3
 write(mmom_unit,*) ' '
 write(mmom_unit,*) '#  Real part of the phonon-modes contribution to the local magnetoelectric tensor(at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(lm_magelsus(i,j,iw)),j=1,3),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

 write(mmom_unit,*) ' '
 write(mmom_unit,*) '#  Imaginary part of the phonon-modes contribution to the local magnetoelectric tensor(at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(lm_magelsus(i,j,iw)),j=1,3),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

 write(mmom_unit,*) ' '
 write(mmom_unit,*) '#  Real part of relaxed-ion local magnetoelectric tensor (at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    ri_magelsus(:,:)=mmom(:,3*(natom+1)+1:3*(natom+2),iw)+lm_magelsus(:,:,iw)
    write(msg,pfmt) &
 &  omega(iw), ((real(ri_magelsus(i,j)),j=1,3),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

 write(mmom_unit,*) ' '
 write(mmom_unit,*) '#  Imaginary part of relaxed-ion local magnetoelectric tensor (at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    ri_magelsus(:,:)=mmom(:,3*(natom+1)+1:3*(natom+2),iw)+lm_magelsus(:,:,iw)
    write(msg,pfmt) &
 &  omega(iw), ((aimag(ri_magelsus(i,j)),j=1,3),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

 write(mmom_unit,*) ' '
 write(mmom_unit,*) '# [Alternative calc] Real part of relaxed-ion local magnetoelectric tensor (at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(ri_genelsus(3*natom+i,j,iw)),j=1,3),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

 write(mmom_unit,*) ' '
 write(mmom_unit,*) '#  [Alternative calc] Imaginary part of relaxed-ion local magnetoelectric tensor (at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(ri_genelsus(3*natom+i,j,iw)),j=1,3),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

 write(mmom_unit,*) ' '
 write(mmom_unit,*) '# [Another alternative calc] Real part of relaxed-ion local magnetoelectric tensor (at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(ri_mmom(i,j,iw)),j=1,3),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

 write(mmom_unit,*) ' '
 write(mmom_unit,*) '#  [Another alternative calc] Imaginary part of relaxed-ion local magnetoelectric tensor (at. units)'
 write(msg,'(a,a,a)') ch10,&
&           ' # At  hw     m_{mat_1,1}^{Ex}     m_{mat_1,1}^{Ey}',&
&           '      ...     m_{mat_1,2}^{Ex}     ...     m_{mat_2,1}^{Ex}     ...'
 call wrtout(mmom_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(ri_mmom(i,j,iw)),j=1,3),i=1,ndim)
    call wrtout(mmom_unit,msg,'COLL')
 end do

 close(mmom_unit)

! mmspec_filename=trim(outfilename_radix)//"_SPECTRAL_MAGMOM"
! if (open_file(mmspec_filename, msg, newunit=mmspec_unit) /= 0) then
!   ABI_ERROR(msg)
! end if
!
! write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  ndim
! write(mmspec_unit,*) ' '
! write(mmspec_unit,'(a)') '#  Real part of magnetic moments spectral function:'
! write(mmspec_unit,*) ' '
! write(msg,'(a,a)') ch10,&
! &           ' # At  hw     m_{mat_1,1}     m_{mat_1,2}     ...     m_{mat_2,1}     m_{mat_2,2}'
! call wrtout(mmspec_unit,msg,'COLL')
! do iw=1,nomega
!   write(msg,pfmt) omega(iw), real(mmomspec(:,iw))
!   call wrtout(mmspec_unit,msg,'COLL')
! end do
!
! write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' )  ndim
! write(mmspec_unit,*) ' '
! write(mmspec_unit,'(a)') '#  Imaginary part of magnetic moments spectral function:'
! write(mmspec_unit,*) ' '
! write(msg,'(a,a)') ch10,&
! &           ' # At  hw     m_{mat_1,1}     m_{mat_1,2}     ...     m_{mat_2,1}     m_{mat_2,2}'
! call wrtout(mmspec_unit,msg,'COLL')
! do iw=1,nomega
!   write(msg,pfmt) omega(iw), aimag(mmomspec(:,iw))
!   call wrtout(mmspec_unit,msg,'COLL')
! end do
!
! close(mmspec_unit)

!Dielectric susceptibility
 diel_filename=trim(outfilename_radix)//"_DIELSUS"
 if (open_file(diel_filename, msg, newunit=diel_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(diel_unit,*) '#'
 write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' ) 9

 write(diel_unit,*) '#  Real part of clamped-ion penalized dielectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
 call wrtout(diel_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(barepsilon(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do

 write(diel_unit,*) ' '
 write(diel_unit,*) '#  Imaginary part of clamped-ion penalized dielectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
 call wrtout(diel_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(barepsilon(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do

 write(diel_unit,*) '#'
 if (mpopt==1) then
   write(diel_unit,*) '#  Frozen-magnetic clamped-ion dielectric tensor calculated and interpolated by ANADDB'
 else if (mpopt==2) then
   write(diel_unit,*) '#  Spin-relaxed clamped-ion dielectric tensor calculated and interpolated by ANADDB'
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
 &  omega(iw), ((real(epsilon(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do

 write(diel_unit,*) ' '
 write(diel_unit,*) '#  Imaginary part of clamped-ion dielectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
 call wrtout(diel_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(epsilon(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do

 write(diel_unit,*) '#'
 if (mpopt==1) then
   write(diel_unit,*) '#  Frozen-magnetic lattice-mediated dielectric tensor calculated and interpolated by ANADDB'
 else if (mpopt==2) then
   write(diel_unit,*) '#  Spin-relaxed lattice-mediated dielectric tensor calculated and interpolated by ANADDB'
 else
   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
   ABI_ERROR(msg)
 end if
 
 write(diel_unit,*) '#'
 write(pfmt, '( "(es15.7, ", I4, "(es15.7))" )' ) 9 

 write(diel_unit,*) '#  Real part of lattice-mediated dielectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
 call wrtout(diel_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(lm_epsilon(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do

 write(diel_unit,*) ' '
 write(diel_unit,*) '#  Imaginary part of lattice-mediated dielectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
 call wrtout(diel_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(lm_epsilon(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do

 write(diel_unit,*) ' '
 write(diel_unit,*) '#  Real part of relaxed-ion dielectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
 call wrtout(diel_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((real(lm_epsilon(i,j,iw)+epsilon(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do

 write(diel_unit,*) ' '
 write(diel_unit,*) '#  Imaginary part of relaxed-ion dielectric tensor'
 write(msg,'(a,a)') ch10,&
&           ' # At  hw     eps_11     eps_12     ...     eps_21     eps_22     ...'
 call wrtout(diel_unit,msg,'COLL')
 do iw=1,nomega
    write(msg,pfmt) &
 &  omega(iw), ((aimag(lm_epsilon(i,j,iw)+epsilon(i,j,iw)),j=1,3),i=1,3)
    call wrtout(diel_unit,msg,'COLL')
 end do

 close(diel_unit)

!Phonon spectral function
 phon_filename=trim(outfilename_radix)//"_SPECTRAL_PHONON"
 if (open_file(phon_filename, msg, newunit=phon_unit) /= 0) then
   ABI_ERROR(msg)
 end if

 write(phon_unit,*) '#'
 if (mpopt==1) then
   write(phon_unit,*) '#  Frozen-magnetic phonon spectral function calculated and interpolated by ANADDB'
 else if (mpopt==2) then
   write(phon_unit,*) '#  Spin-relaxed phonon spectral function calculated and interpolated by ANADDB'
 else
   write(msg,'(a)') 'ddb_omega_interpol: variable mpopt just can be 1 or 2'
   ABI_ERROR(msg)
 end if
 
 write(msg,'(a,a)') ch10, ' # At  hw               Phonon SF           Magnon-phonon SF'
 call wrtout(phon_unit,msg,'COLL')
 do iw=1,nomega
   write(msg,*) omega(iw), phonspec(iw), magphonspec(iw)
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
   write(phon_unit,*) '#  Frozen-magnetic phonon frequencies calculated and interpolated by ANADDB'
 else if (mpopt==2) then
   write(phon_unit,*) '#  Spin-relaxed phonon frequencies calculated and interpolated by ANADDB'
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

!!Born effective charges
! zeff_filename=trim(outfilename_radix)//"_ZEFF"
! if (open_file(zeff_filename, msg, newunit=zeff_unit) /= 0) then
!   ABI_ERROR(msg)
! end if
!
! write(zeff_unit,*) '#'
! write(zeff_unit,*) '#  Born effective charges calculated and interpolated by ANADDB'
! write(zeff_unit,*) '#'
!
! write(pfmt, '( "(es15.7, ", I2, "(es17.7))" )' ) 3 
! do imode= 1, 3*natom
!   write(zeff_unit,*) ' '
!   write(zeff_unit,'(a,i3)') '#  Real part of Born charge (at. units) induced by phonon mode:', imode
!   write(msg,'(a,a)') ch10,&
! &           ' # At  hw     Z^x_{n}     Z^y_{n}     Z^z_{n}'
!   call wrtout(zeff_unit,msg,'COLL')
!   do iw=1,nomega
!     write(msg,pfmt) &
!   & omega(iw), (real(modezeff(i,imode,iw)),i=1,3)
!     call wrtout(zeff_unit,msg,'COLL')
!   end do
!   write(zeff_unit,*) ' '
!   write(zeff_unit,'(a,i3)') '#  Imaginary part of Born charge (at. units) induced by phonon mode:', imode
!   write(msg,'(a,a)') ch10,&
! &           ' # At  hw     Z^x_{n}     Z^y_{n}     Z^z_{n}'
!   call wrtout(zeff_unit,msg,'COLL')
!   do iw=1,nomega
!     write(msg,pfmt) &
!   & omega(iw), (aimag(modezeff(i,imode,iw)),i=1,3)
!     call wrtout(zeff_unit,msg,'COLL')
!   end do
! end do
!
! close(zeff_unit)

! zeffspec_filename=trim(outfilename_radix)//"_SPECTRAL_ZEFF"
! if (open_file(zeffspec_filename, msg, newunit=zeffspec_unit) /= 0) then
!   ABI_ERROR(msg)
! end if
! write(zeffspec_unit,*) ' '
! write(zeffspec_unit,'(a)') '#  Real part of Born charges spectral function:'
! write(zeffspec_unit,*) ' '
! write(msg,'(a,a)') ch10,&
! &           ' # At  hw     Z^x     Z^y     Z^z'
! call wrtout(zeffspec_unit,msg,'COLL')
! do iw=1,nomega
!   write(msg,'(4es15.7)') omega(iw), real(zeffspec(:,iw))
!   call wrtout(zeffspec_unit,msg,'COLL')
! end do
!
! write(zeffspec_unit,*) ' '
! write(zeffspec_unit,'(a)') '#  Imaginary part of Born charges spectral function:'
! write(zeffspec_unit,*) ' '
! write(msg,'(a,a)') ch10,&
! &           ' # At  hw     Z^x     Z^y     Z^z'
! call wrtout(zeffspec_unit,msg,'COLL')
! do iw=1,nomega
!   write(msg,'(4es15.7)') omega(iw), aimag(zeffspec(:,iw))
!   call wrtout(zeffspec_unit,msg,'COLL')
! end do
!
! close(zeffspec_unit)

 ABI_FREE(dint_barddb)
 ABI_FREE(int_barddb)
 ABI_FREE(barmagsus)
 ABI_FREE(magsus)
 ABI_FREE(lm_magsus)
 ABI_FREE(invbarmagsus)
 ABI_FREE(invmagsus)
 ABI_FREE(invhmat)
 ABI_FREE(barmmom)
 ABI_FREE(barmmom_tr)
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
!!  int_barddb(2,ddb%msize,1)= interpolated hessian at omega
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
& invmagsus,magsus,mpatpol,mpdir,mpert,msize,natom,ndim,nmdir,int_barddb, &
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
 real(dp), intent(out) :: int_barddb(2,msize,1)
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
 int_barddb(1,:,1)= reshape( lhess(1,:,:,:,:), shape = (/msize/) )
 int_barddb(2,:,1)= reshape( lhess(2,:,:,:,:), shape = (/msize/) )

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


subroutine phonon_green(amu,eigvec,eigvec_fm,eta_phongreen,ifc,ifc_fm,invmagsus,& 
& magphongreen,magphonspec,mode_magphonspec,mode_phonspec,natom,ndim,ntypat,omega,&
& phfrq,phongreen,phongreen_fm,phonspec,typat,zfield,zfield_tr)

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: natom,ndim,ntypat 
 real(dp), intent(in) :: eta_phongreen,omega
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

 cplx_eta=cmplx(0.0_dp,eta_phongreen)
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
!   mode_phonspec(irow)= -one/pi * aimag(two*cmplx(omega,eta_phongreen)*work1(irow,irow))
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

end module m_ddb_magpen
!!***
