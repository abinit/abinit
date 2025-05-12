!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ddb_magpen
!! NAME
!!  m_ddb_magpen
!!
!! FUNCTION
!!  Convert second- and -third (Berry curvatures) order total energy derivatives
!!  calculated with the magnetic penalty (constrained-B functional) into the 
!!  corresponding quantities of different magnetic functionals: 
!!  --constrained-M (fixed-spin)
!!  --constrained-H (relaxed spin)
!!  Calculate and write the ensuing clamped-ion susceptibilities. Lattice-mediated
!!  contributions are incorporated in m_ddb_omega_interpol.
!!
!! COPYRIGHT
!!  Copyright (C) 2023 ABINIT group (MR and MS)
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
 public :: local_spinsus    ! Treat local spin susceptibility (2nd-order magnetic derivatives)
 public :: magmom           ! Treat first-order magnetic moments (2nd-order mixed derivatives)
 public :: mp_d2etot        ! Treat 2nd-order nonmagnetic derivatives
 public :: asrw0            ! Apply the ASR correction calculated at w=0 at any value of w

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
 integer :: nmat,nmdir,optgb,prtopt
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
 optgb=1
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
! if (mpopt==2) ABI_MALLOC(ddb%val_rs,(2,ddb%msize,ddb%nblok))
 ABI_MALLOC(ddb%val_rs,(2,ddb%msize,ddb%nblok))

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
   & mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtopt,prtvol)

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
   & mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtopt,prtvol,qphon,xred,zfield,zfield_tr)
   end if

   !Now calculate the non-magnetic second-order quantities
   write(msg, '(2a,(80a),4a)' ) ch10,('-',ii=1,80),ch10,ch10,&
   ' Second-order linear-response tensors ',ch10
   call wrtout([std_out, ab_out], msg)

   !Convert ddb%val to second-order energies
   call ddb%to_d2etot(ddb%val,kblok,0,qeq0,qphon,qphnrm,ucvol,optgb,omega=omega)

   !Convert second-order derivatives to diferent magnetic boundary conditions
   call mp_d2etot(barmagsus,ddb,kblok,invhmat,magsus,magpen,mpert,mpopt,natom, &
 & nblok,ndim,qphon,xred,zfield,zfield_tr)

   !Convert second-order energies to the physical quantities of ddb%val
   call ddb%to_d2etot(ddb%val,kblok,1,qeq0,qphon,qphnrm,ucvol,optgb,omega=omega)
   call ddb%to_d2etot(ddb%val_fs,kblok,1,qeq0,qphon,qphnrm,ucvol,optgb,omega=omega)
   if (mpopt==2) call ddb%to_d2etot(ddb%val_rs,kblok,1,qeq0,qphon,qphnrm,ucvol,optgb,omega=omega)

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
   
 end do

 ! BERRY CURVATURES
 if (timdisp==1) then

!   if (dissip==1) then
!     ABI_BUG("Berry curvatures calculation is not implemented with dissipation, set dissip=0")
!   end if
   ABI_MALLOC(ddb_lw%val_fs,(2,ddb_lw%msize,ddb_lw%nblok))

   ABI_MALLOC(bc_barmagsus,(ndim,ndim))
   ABI_MALLOC(bc_ss,(ndim,ndim))
   ABI_MALLOC(bc_sp,(ndim,(natom+5)*3))

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
       call berrycurv_ss(bc_barmagsus,bc_ss,ddb_lw,iblok,invbarmagsus,mpatpol,mpdir,mpert,&
     & natom,nblok,ndim,nmdir,prtvol)
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

     ! Then macroscopic Zeeman field
     lblok=0
     if (qeq0) then
       rfelfd(2)=0
       rfmagn(2)=1

       call ddb_lw%get_block(lblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, 33, &
     & mpatpol=mpatpol,mpdir=mpdir,omega=omega,rfmagn=rfmagn,rffreq=rffreq)
     end if

     if (iblok /= 0 .or. jblok /=0 .or. lblok /= 0) then
       call berrycurv_sp(barmmom,bc_sp,bc_ss,ddb_lw,iblok,invbarmagsus,jblok,lblok, &
     & mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol,qphon,xred)
     end if

     !Berry curvature of other second-order quantites
     call berrycurv_pp(barmagsus,bc_barmagsus,bc_sp,ddb_lw,kblok, &
   & mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol,qeq0,ucvol,xred,zfield)

     !Print them
     call mp_d3etot_print(ddb_lw,ddb_lw%val_fs,kblok,mpert,natom,nblok,1,omega,&
   & prtvol,qeq0,qphnrm,qphon,ucvol)

   end do
   ABI_FREE(bc_ss)
   ABI_FREE(bc_sp)
   ABI_FREE(bc_barmagsus)
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
!! fs2rs= (optional) if 1, the routine starts from a precalculated blkval_fs 
!! blkval_fs(2,3,mpert,3,mpert)= fixed-spin 2nd-order derivatives
!!
!! OUTPUT
!! barmagsus(ndim,ndim)= Penalized spin-sussceptibility tensor
!! invbarmagsus(ndim,ndim)= Inverse of the penalized spin-sussceptibility tensor
!! magsus(ndim,ndim)= Spin-sussceptibility tensor
!! invmagsus(ndim,ndim)= Inverse of spin-sussceptibility tensor
!!
!! SOURCE

 subroutine local_spinsus(barmagsus,ddb,iblok,invbarmagsus,invmagsus,invhmat,magpen,magsus, &
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtopt,prtvol, &
& fs2rs,blkval_fs) !optional

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok,ndim,nmdir
 integer,intent(in) :: prtopt,prtvol
 integer,intent(in),optional :: fs2rs
 real(dp),intent(in) :: magpen
!arrays
 type(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(in),optional :: blkval_fs(2,3,mpert,3,mpert,1)
 complex(dpc),intent(out) :: barmagsus(ndim,ndim)
 complex(dpc),intent(out) :: invbarmagsus(ndim,ndim)
 complex(dpc),intent(out) :: invhmat(ndim,ndim)
 complex(dpc),intent(out) :: magsus(ndim,ndim)
 complex(dpc),intent(out) :: invmagsus(ndim,ndim)

!Local variables -------------------------
!scalars
 integer :: fs2rs_
 integer :: iat1,iat2,icol,idir1,idir2,index,info,ipert1,ipert2,irow,lwork
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
 real(dp) :: fac
 character(len=1000) :: msg
!arrays
 complex(dpc) :: idty(ndim,ndim)
 integer :: indexat(ndim),indexdir(ndim)
 integer, allocatable :: ipiv(:)
 complex(dpc),allocatable :: work(:),work1(:,:),work2(:,:)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

!If fixed-spin case has been precalculated do less stuff
 fs2rs_=0; if (present(fs2rs)) fs2rs_=fs2rs
 if (fs2rs_==1) then
   if (.not.present(blkval_fs)) then
     write(msg, '(3a)' )' No fixed-spin array has been passed to local_spinsus', &
   & ' but fs2rs=1 ',ch10
     ABI_ERROR(msg)
   end if 
 end if

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

         if (fs2rs_==0) then
           barmagsus(irow,icol)= &
         & cmplx(ddb%val(1,index,iblok),ddb%val(2,index,iblok),16)
         else if (fs2rs_==1) then
           invmagsus(irow,icol)= &
         & cmplx(blkval_fs(1,idir1,ipert1,idir2,ipert2,iblok), &
         & blkval_fs(2,idir1,ipert1,idir2,ipert2,iblok),16)
         end if

       end do
     end do
   end do
 end do

!Use magsus to store the intermediate array
 magsus=idty-magpen*barmagsus
     
!Invert the arrays
 ABI_MALLOC(work1,(ndim,ndim))
 ABI_MALLOC(work2,(ndim,ndim))
 if (fs2rs_==0) then
   work1=barmagsus
 else if (fs2rs_==1) then
   work1=invmagsus
 end if
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

 if (fs2rs_==1) then
   magsus=work1
   return
 end if

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
 invbarmagsus=work1
 
 !At last, calculate the susceptibility and its inverse
 invhmat=work2
 magsus=matmul(work2,barmagsus)
 invmagsus=invbarmagsus-magpen*idty

 ABI_FREE(ipiv)
 ABI_FREE(work1)
 ABI_FREE(work2)

!!TMP shift of invmagsus
! do irow= 1, ndim
!   invmagsus(irow,irow)= invmagsus(irow,irow) - 2.072d-4 
! end do 
! invmagsus(1,3)= invmagsus(1,3) + 1.924d-4
! invmagsus(2,4)= invmagsus(2,4) + 1.924d-4
! invmagsus(3,1)= invmagsus(3,1) + 1.924d-4
! invmagsus(4,2)= invmagsus(4,2) + 1.924d-4

 if (prtopt==1.and.prtvol>1) then
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
!! fs2rs= (optional) if 1, the routine starts from a precalculated blkval_fs 
!! blkval_fs(2,3,mpert,3,mpert)= fixed-spin 2nd-order derivatives
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
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtopt,prtvol,qphon,xred,zfield,zfield_tr, &
& fs2rs,blkval_fs) !optional

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,jblok,lblok,mpert,natom,nblok,ndim,nmdir,prtopt,prtvol
 integer,intent(in),optional :: fs2rs
 real(dp),intent(in) :: magpen
!arrays
 type(ddb_type),intent(inout) :: ddb
 real(dp),intent(in) :: qphon(3),xred(3,natom)
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(in),optional :: blkval_fs(2,3,mpert,3,mpert,1)
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
 integer :: fs2rs_
 integer :: iat1,iat2,icol,idir1,idir2,index,ipert1,ipert2,irow
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red,jndex,zblok
 real(dp) :: fac
 character(len=1000) :: msg
!arrays
 integer :: indexat1(ndim),indexdir1(ndim)
 integer :: indexat2((natom+5)*3),indexdir2((natom+5)*3)
 complex(dpc) :: mmom_alt(ndim,(natom+5)*3)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

!If fixed-spin case has been precalculated do less stuff
 fs2rs_=0; if (present(fs2rs)) fs2rs_=fs2rs
 if (fs2rs_==1) then
   if (.not.present(blkval_fs)) then
     write(msg, '(3a)' )' No fixed-spin array has been passed to magmom', &
   & ' but fs2rs=1 ',ch10
     ABI_ERROR(msg)
   end if 
 end if

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

         if (fs2rs_==0) then
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
         else if (fs2rs_==1) then
           zfield(irow,icol)= &
         & cmplx(blkval_fs(1,idir1,ipert1,idir2,ipert2,iblok), &
         & blkval_fs(2,idir1,ipert1,idir2,ipert2,iblok),16)
           zfield_tr(icol,irow)= &
         & cmplx(blkval_fs(1,idir2,ipert2,idir1,ipert1,iblok), &
         & blkval_fs(2,idir2,ipert2,idir1,ipert1,iblok),16)
         end if

       end do
     end do
   end do
 end do

!Compute the Zeeman fields 
 if (fs2rs_==0) then
   zfield=-matmul(invbarmagsus,barmmom)
   zfield_tr=-matmul(barmmom_tr,invbarmagsus)
 end if

!Compute the moments
 if (fs2rs_==0) then
   mmom=-matmul(magsus,zfield)
   mmom_alt=matmul(invhmat,barmmom)
   mmom_tr=matmul(barmmom_tr,invhmat)
 else if (fs2rs_==1) then
   mmom=-matmul(magsus,zfield)
   mmom_tr=-matmul(zfield_tr,magsus)
   return
 end if

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
       call wrtout([ab_out,std_out], '  atom1  dir  E-dir            Real              Imag')
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
       call wrtout([ab_out,std_out], '  atom1  dir  E-dir            Real              Imag')
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
       call wrtout([ab_out,std_out], '  atom1  dir  E-dir            Real              Imag')
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
       call wrtout([ab_out,std_out], '  atom1  dir  B-dir            Real              Imag')
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
       call wrtout([ab_out,std_out], '  atom1  dir  B-dir            Real              Imag')
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
       call wrtout([ab_out,std_out], '  atom1  dir  B-dir            Real              Imag')
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
         call wrtout([ab_out,std_out], '  atom1  dir  E-dir            Real              Imag')
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
         call wrtout([ab_out,std_out], '  atom1  dir  B-dir            Real              Imag')
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
 do ipert2=1,natom+5
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
!! fs2rs= (optional) if 1, the routine starts from a precalculated blkval_fs 
!! blkval_fs(2,3,mpert,3,mpert,1)= fixed-spin 2nd-order derivatives
!!
!! OUTPUT
!! ddb%val_fs(2,msize,nblok)= second-order derivatives at fixed spin.
!! ddb%val_rs(2,msize,nblok)= second-order derivatives at relaxed spin.
!! blkval_rs(2,3,mpert,3,mpert,1)= (optional) relaxed-spin 2nd-order derivatives
!!
!! SOURCE

 subroutine mp_d2etot(barmagsus,ddb,&
& iblok,invhmat,magsus,magpen,mpert,mpopt,&
& natom,nblok,ndim,qphon,xred,zfield,zfield_tr,&
& fs2rs,blkval_fs,blkval_rs) !optional

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,mpopt,natom,nblok,ndim
 integer,intent(in),optional :: fs2rs
 real(dp),intent(in) :: magpen
!arrays
 type(ddb_type),intent(inout) :: ddb
 real(dp),intent(in) :: qphon(3),xred(3,natom)
 real(dp),intent(in),optional :: blkval_fs(2,3,mpert,3,mpert,1)
 real(dp),intent(out),optional :: blkval_rs(2,3,mpert,3,mpert,1)
 complex(dpc),intent(in) :: barmagsus(ndim,ndim)
 complex(dpc),intent(in) :: invhmat(ndim,ndim)
 complex(dpc),intent(in) :: magsus(ndim,ndim)
 complex(dpc),intent(in) :: zfield(ndim,(natom+5)*3)
 complex(dpc),intent(in) :: zfield_tr((natom+5)*3,ndim)
!Local variables -------------------------
!scalars
 integer :: fs2rs_
 integer :: idir1,idir2,ipert1,ipert2,index,irow,icol
 complex(dpc) :: val_ps,val_fs, val_rs
 character(len=1000) :: msg
!arrays
 
! *********************************************************************

!If fixed-spin case has been precalculated do less stuff
 fs2rs_=0; if (present(fs2rs)) fs2rs_=fs2rs
 if (fs2rs_==1) then
   if (.not.present(blkval_fs)) then
     write(msg, '(3a)' )' No fixed-spin array has been passed to mp_d2etot', &
   & ' but fs2rs=1 ',ch10
     ABI_ERROR(msg)
   end if 
   if (.not.present(blkval_rs)) then
     write(msg, '(3a)' )' No relaxed-spin array has been passed to mp_d2etot', &
   & ' but fs2rs=1 ',ch10
     ABI_ERROR(msg)
   end if 
 end if

!Extract the penalized/constrained quantities
 do ipert2= 1, natom+5
   do idir2= 1, 3
     icol= (ipert2-1)*3 + idir2
     do ipert1= 1, natom+5
       do idir1= 1, 3
         irow= (ipert1-1)*3 + idir1
         index= idir1 + 3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))

         if (fs2rs_==0) then
           !Extract the penalized second-order derivatives 
           val_ps= cmplx(ddb%val(1,index,iblok),ddb%val(2,index,iblok),16)
           !Calculate the fixed-spin flavor
           val_fs= val_ps + &
         & sum( zfield_tr(irow,:) * matmul( barmagsus,zfield(:,icol) ) ) 
           ddb%val_fs(1,index,iblok)= real(val_fs)
           ddb%val_fs(2,index,iblok)= aimag(val_fs)
         else if (fs2rs_==1) then
           val_fs= &
         & cmplx(blkval_fs(1,idir1,ipert1,idir2,ipert2,iblok), &
         & blkval_fs(2,idir1,ipert1,idir2,ipert2,iblok),16)
         end if

         if (mpopt==2) then
           !Calculate the relaxed-spin flavor
           val_rs= val_fs - &
         & sum( zfield_tr(irow,:) * matmul( magsus,zfield(:,icol) ) ) 
           if (fs2rs_==0) then
             ddb%val_rs(1,index,iblok)= real(val_rs)
             ddb%val_rs(2,index,iblok)= aimag(val_rs)
           else if (fs2rs_==1) then
             blkval_rs(1,idir1,ipert1,idir2,ipert2,iblok)=real(val_rs)
             blkval_rs(2,idir1,ipert1,idir2,ipert2,iblok)=aimag(val_rs)
           end if
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
     call wrtout([ab_out,std_out], ' E-dir      atom   dir        Real              Imag')
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
   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega, &
 & rfmagn=rfmagn)
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
 call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega, &
& rfmagn=rfmagn)
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
 call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega, &
& rfmagn=rfmagn)
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

 subroutine berrycurv_ss(bc_barmagsus,bc_ss,ddb_lw,iblok,invbarmagsus,&
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok,ndim,nmdir,prtvol
!arrays
 type(ddb_type),intent(inout) :: ddb_lw
 integer,intent(in) :: mpatpol(2),mpdir(3)
 complex(dpc),intent(in) :: invbarmagsus(ndim,ndim)
 complex(dpc),intent(out) :: bc_ss(ndim,ndim)
 complex(dpc),intent(out) :: bc_barmagsus(ndim,ndim)
!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,idir3,index
 integer :: info,ipert1,ipert2,ipert3,irow,lwork
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
 real(dp) :: fac
 complex(dpc), parameter :: ione=(0.d0,1.d0)
 character(len=1000) :: msg
!arrays
 complex(dpc) :: idty(ndim,ndim)
 integer :: indexat(ndim),indexdir(ndim)
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
         index = idir1 + &
       & 3*((ipert1 - 1) + mpert*((idir2 - 1) + &
       & 3*((ipert2 -1 ) + mpert*((idir3 - 1) + 3*(ipert3 - 1)))))

         bc_barmagsus(irow,icol)= -one* &
       & cmplx(ddb_lw%val(1,index,iblok),ddb_lw%val(2,index,iblok),16)

       end do
     end do
   end do
 end do

!Calculate the Berry-curvature of the inverse magnetic susceptibility
 bc_ss=-matmul(invbarmagsus,matmul(bc_barmagsus,invbarmagsus)) 

!TODO: Adapt the phase of this quantity
! fac=2.714943600699**2/four !TMP
! fac=one/four !TMP
! open(10,file='g_ss.txt')
!   do irow=1, ndim
!     write(10,*) -ione*bc_ss(irow,1:ndim)*fac
!   end do 
! close(10)


 if (prtvol>1) then
   call wrtout([ab_out,std_out], ' Fixed-spin Berry curvature of the inverse spin susceptibility ')
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
 end if

!Store the FM flavor in the DDB file
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
         index = idir1 + &
       & 3*((ipert1 - 1) + mpert*((idir2 - 1) + &
       & 3*((ipert2 -1 ) + mpert*((idir3 - 1) + 3*(ipert3 - 1)))))

         ddb_lw%val_fs(1,index,iblok)= real(bc_ss(irow,icol))
         ddb_lw%val_fs(2,index,iblok)= aimag(bc_ss(irow,icol))
       end do
     end do
   end do
 end do

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
!! ddb_lw=  Third-order derivative ddb
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

 subroutine berrycurv_sp(barmmom,bc_sp,bc_ss,ddb_lw,iblok,invbarmagsus,&
& jblok,lblok,mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol,qphon,xred)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,jblok,lblok,mpert,natom,nblok,ndim,nmdir,prtvol
!arrays
 type(ddb_type),intent(inout) :: ddb_lw
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(in) :: qphon(3),xred(3,natom)
 complex(dpc),intent(in) :: barmmom(ndim,(natom+5)*3)
 complex(dpc),intent(in) :: bc_ss(ndim,ndim)
 complex(dpc),intent(in) :: invbarmagsus(ndim,ndim)
 complex(dpc),intent(out) :: bc_sp(ndim,(natom+5)*3)
!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,idir3,index,ipert1,ipert2,ipert3,irow
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red,jndex
 real(dp) :: fac,re,im
 complex(dpc), parameter :: ione=(0.d0,1.d0)
 character(len=1000) :: msg
!arrays
 integer :: indexat1(ndim),indexdir1(ndim)
 integer :: indexat2((natom+5)*3),indexdir2((natom+5)*3)
 complex(dpc) :: bc_barsp(ndim,(natom+5)*3)
! complex(dpc) :: bc_ps((natom+2)*3,ndim)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

!Extract the berry curvature of the penalized moments
 bc_barsp=(zero,zero)
 ipert3= natom + 9
 idir3= 1
 do ipert2=1,natom+5
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
         index = idir1 + &
       & 3*((ipert1 - 1) + mpert*((idir2 - 1) + &
       & 3*((ipert2 -1 ) + mpert*((idir3 - 1) + 3*(ipert3 - 1)))))
         
         if (iblok /=0 .and. ipert2 <= natom) then
           bc_barsp(irow,icol)= -one* &
         & cmplx(ddb_lw%val(1,index,iblok),ddb_lw%val(2,index,iblok),16)
         else if (jblok /=0 .and. ipert2 == natom+2) then
           bc_barsp(irow,icol)= -one* &
         & cmplx(ddb_lw%val(1,index,jblok),ddb_lw%val(2,index,jblok),16)
         else if (lblok /=0 .and. ipert2 == natom+5) then
           bc_barsp(irow,icol)= -one* &
         & cmplx(ddb_lw%val(1,index,lblok),ddb_lw%val(2,index,lblok),16)
         end if

       end do
     end do
   end do
 end do

 !Calculate the Berry curvature of the induced Zeeman fields
 bc_sp= -matmul(bc_ss,barmmom) - matmul(invbarmagsus,bc_barsp)

!TO DO: apply the q-depenent phase on the magnetic variables
! do irow=1,ndim
!   do iat1= 1, natom
!     do idir1= 1, 3
!       icol= (iat1-1)*3 + idir1
!       !MR: Caution, this conjg might be wrong in presence of dissipation
!       bc_ps(icol,irow)=conjg(bc_sp(irow,icol)*exp(two_pi*(0.d0,1.d0)* dot_product(qphon,xred(:,iat1))))
!     end do
!   end do
! end do 
!
!! fac=2.714943600699/two/0.52917 !TMP
! fac=one/two/0.52917 !TMP
! open(10,file='g_ps.txt')
! do irow=1,natom*3
!   write(10,*) (0.d0,-1.d0)*bc_ps(irow,1:ndim)*fac
! end do 
! close(10)

 if (prtvol > 1) then
   if (iblok /= 0) then
     call wrtout([ab_out,std_out], ' Fixed-spin Berry curvature of the Zeeman fields induced by atomic displacements')
     call wrtout([ab_out,std_out], '  atom1  dir  E-dir            Real              Imag')
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
     call wrtout([ab_out,std_out], ' Fixed-spin Berry curvature of the Zeeman fields induced by electric field')
     call wrtout([ab_out,std_out], '  atom1  dir  E-dir            Real              Imag')
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
   if (lblok /= 0) then
     call wrtout([ab_out,std_out], ' Fixed-spin Berry curvature of the Zeeman fields induced by macroscopic Zeeman field')
     call wrtout([ab_out,std_out], '  atom1  dir  B-dir            Real              Imag')
     do irow=1, ndim
       do icol=(natom+4)*3+1, (natom+5)*3
         write(msg,'(i4,4x,a2,4x,a2,6x,2es18.9)' ) &
       & indexat1(irow), cart(indexdir1(irow)), cart(indexdir2(icol)), &
       & real(bc_sp(irow,icol)), aimag(bc_sp(irow,icol))
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], '   ')
   end if
 end if

!Store the FM flavor in the DDB file
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
         index = idir1 + &
       & 3*((ipert1 - 1) + mpert*((idir2 - 1) + &
       & 3*((ipert2 -1 ) + mpert*((idir3 - 1) + 3*(ipert3 - 1)))))
         jndex = idir2 + &
       & 3*((ipert2 - 1) + mpert*((idir1 - 1) + &
       & 3*((ipert1 -1 ) + mpert*((idir3 - 1) + 3*(ipert3 - 1)))))
         
         if (iblok /=0 .and. ipert2 <= natom) then
           ddb_lw%val_fs(1,index,iblok)=real(bc_sp(irow,icol))
           ddb_lw%val_fs(2,index,iblok)=aimag(bc_sp(irow,icol))
           ddb_lw%val_fs(1,jndex,iblok)=real(bc_sp(irow,icol))
           ddb_lw%val_fs(2,jndex,iblok)=-aimag(bc_sp(irow,icol))
         else if (jblok /=0 .and. ipert2 == natom+2) then
           ddb_lw%val_fs(1,index,jblok)=real(bc_sp(irow,icol))
           ddb_lw%val_fs(2,index,jblok)=aimag(bc_sp(irow,icol))
           ddb_lw%val_fs(1,jndex,jblok)=real(bc_sp(irow,icol))
           ddb_lw%val_fs(2,jndex,jblok)=-aimag(bc_sp(irow,icol))
         else if (lblok /=0 .and. ipert2 == natom+5) then
           ddb_lw%val_fs(1,index,lblok)=real(bc_sp(irow,icol))
           ddb_lw%val_fs(2,index,lblok)=aimag(bc_sp(irow,icol))
           ddb_lw%val_fs(1,jndex,lblok)=real(bc_sp(irow,icol))
           ddb_lw%val_fs(2,jndex,lblok)=-aimag(bc_sp(irow,icol))
         end if

       end do
     end do
   end do
 end do

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
!! kblok= index of the current block
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

 subroutine berrycurv_pp(barmagsus,bc_barmagsus,bc_sp,ddb_lw,kblok,&
& mpatpol,mpdir,mpert,natom,nblok,ndim,nmdir,prtvol,qeq0,ucvol,xred,zfield)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: kblok,mpert,natom,nblok,ndim,nmdir,prtvol
 real(dp),intent(in) :: ucvol
 logical,intent(in) :: qeq0
!arrays
 type(ddb_type),intent(inout) :: ddb_lw
 integer,intent(in) :: mpatpol(2),mpdir(3)
 real(dp),intent(in) :: xred(3,natom)
 complex(dpc),intent(in) :: barmagsus(ndim,ndim)
 complex(dpc),intent(in) :: bc_barmagsus(ndim,ndim)
 complex(dpc),intent(in) :: bc_sp(ndim,(natom+5)*3)
 complex(dpc),intent(in) :: zfield(ndim,(natom+5)*3)
!Local variables -------------------------
!scalars
 integer :: iat1,iat2,icol,idir1,idir2,idir3,index,ipert1,ipert2,ipert3,irow
 integer :: ipert1_red,ipert2_red,idir1_red,idir2_red
 real(dp) :: fac
 complex(dpc), parameter :: ione=(0.d0,1.d0)
 complex(dpc) :: cval
 character(len=1000) :: msg
!arrays
 complex(dpc) :: bc_barpp((natom+5)*3,(natom+5)*3), bc_pp((natom+5)*3,(natom+5)*3) 
 complex(dpc) :: term((natom+5)*3,(natom+5)*3,3)
 character(len=1) :: cart(3)=(/'x','y','z'/)

! *********************************************************************

!Extract the frequency derivative of the penalized 2nd order nonmagnetic quantities
 ipert3= natom + 9
 idir3= 1
 do ipert2= 1, natom+5
   do idir2= 1, 3
     icol=( ipert2-1)*3 + idir2
     do ipert1= 1, natom+5
       do idir1= 1, 3
         irow=( ipert1-1)*3 + idir1
         index = idir1 + &
       & 3*((ipert1 - 1) + mpert*((idir2 - 1) + &
       & 3*((ipert2 -1 ) + mpert*((idir3 - 1) + 3*(ipert3 - 1)))))

         bc_barpp(irow,icol)= cmplx(ddb_lw%val(1,index,kblok),ddb_lw%val(2,index,kblok),16)

       end do
     end do
   end do
 end do 

!Calculate the different terms entering the Berry curvature
 term(:,:,1)=matmul(transpose(conjg(bc_sp)),matmul(barmagsus,zfield))
 term(:,:,2)=matmul(transpose(conjg(zfield)),matmul(bc_barmagsus,zfield))
 term(:,:,3)=matmul(transpose(conjg(zfield)),matmul(barmagsus,bc_sp))

 bc_pp(:,:)= bc_barpp(:,:) + term(:,:,1) + term(:,:,2) + term(:,:,3)

!Apply factors to convert derivatives of energy into susceptibilities
!Born charges
 if (qeq0) then
   ipert1= natom + 2
   do ipert2= 1, natom
     do idir2= 1, 3
       icol=(ipert2-1)*3 + idir2
       do idir1= 1, 3
         irow=(ipert1-1)*3 + idir1
         cval=bc_pp(irow,icol)
         bc_pp(irow,icol)=-one*cval
         cval=bc_pp(icol,irow)
         bc_pp(icol,irow)=-one*cval
       end do
     end do
   end do
 end if

!Magnetic charges induced by atomic displacement
 ipert1= natom + 5
 do ipert2= 1, natom
   do idir2= 1, 3
     icol=( ipert2-1)*3 + idir2
     do idir1= 1, 3
       irow=( ipert1-1)*3 + idir1
       cval=bc_pp(irow,icol)
       bc_pp(irow,icol)=-cval
       cval=bc_pp(icol,irow)
       bc_pp(icol,irow)=-cval
     end do
   end do
 end do

!Dielectric tensor
 if (qeq0) then
   ipert1= natom + 2
   ipert2= natom + 2
   do idir2= 1, 3
     icol=( ipert2-1)*3 + idir2
     do idir1= 1, 3
       irow=( ipert1-1)*3 + idir1
       cval=bc_pp(irow,icol)
       bc_pp(irow,icol)=-four_pi/ucvol*cval
     end do
   end do
 end if

!TODO: the next two susceptibilities miss a 1/ucvol factor that needs to be first
!incorporated in the 2nd-order susceptibilities of ABINIT.

!Magnetoelectric susceptibility
 if (qeq0) then
   ipert1= natom + 5
   ipert2= natom + 2
   do idir2= 1, 3
     icol=( ipert2-1)*3 + idir2
     do idir1= 1, 3
       irow=( ipert1-1)*3 + idir1
       cval=bc_pp(irow,icol)
       bc_pp(irow,icol)=-cval
       cval=bc_pp(icol,irow)
       bc_pp(icol,irow)=-cval
     end do
   end do
 end if

!Magnetic susceptibility
  ipert1= natom + 5
  ipert2= natom + 5
  do idir2= 1, 3
    icol=( ipert2-1)*3 + idir2
    do idir1= 1, 3
      irow=( ipert1-1)*3 + idir1
      cval=bc_pp(irow,icol)
      bc_pp(irow,icol)=-cval
    end do
  end do

!Store the FM flavor in the DDB file
 do ipert2= 1, natom+5
   do idir2= 1, 3
     icol=( ipert2-1)*3 + idir2
     do ipert1= 1, natom+5
       do idir1= 1, 3
         irow=( ipert1-1)*3 + idir1
         index = idir1 + &
       & 3*((ipert1 - 1) + mpert*((idir2 - 1) + &
       & 3*((ipert2 -1 ) + mpert*((idir3 - 1) + 3*(ipert3 - 1)))))
         ddb_lw%val_fs(1,index,kblok)=real(bc_pp(irow,icol))
         ddb_lw%val_fs(2,index,kblok)=aimag(bc_pp(irow,icol))
       end do
     end do
   end do
 end do 

 end subroutine berrycurv_pp
!!***

!!****f* m_ddb_magpen/mp_d3etot_print
!! NAME
!! mp_d3etot_print
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

 subroutine mp_d3etot_print(ddb_lw,blkval,kblok,mpert,natom,nblok,opt,omega,prtvol,qeq0,qphnrm,qphon,ucvol)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(in) :: ddb_lw
 integer,intent(in) :: kblok,mpert,natom,nblok,opt,prtvol
 logical,intent(in) :: qeq0
 real(dp),intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: omega(3)
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,3,mpert,nblok)
 real(dp),intent(inout) :: qphnrm(3),qphon(3,3)

!Local variables -------------------------
!scalars
 integer :: iblok,idir1,idir2,idir3,ipert1,ipert2,ipert3,index,irow,icol
 integer :: rftyp
 character(len=1000) :: msg
!arrays
 integer :: rfelfd(4),rfmagn(4),rfphon(4),rfstrs(4),rffreq(4)
 real(dp) :: val(2)
 character(len=1) :: cart(3)=(/'x','y','z'/)
 
! *********************************************************************

 rfelfd(:)=0
 rfphon(:)=0
 rfstrs(:)=0
 rfmagn(:)=0
 rffreq(:)=0
 rffreq(3)=1
 rftyp = 33

 ipert3= natom + 9
 idir3= 1

 !IFCs
 if (prtvol>1) then
   rfphon(1:2)=1
   call ddb_lw%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega, &
  & rffreq=rffreq)
   if (iblok/=0.and.iblok==kblok) then
     if (opt==1) then
       call wrtout([ab_out,std_out], ' Frozen-spin Berry curvature of interatomic force constants')
     else if (opt==2) then
       call wrtout([ab_out,std_out], ' Relaxed-spin Berry curvature of interatomic force constants')
     end if
     call wrtout([ab_out,std_out], '  atom1  dir  atom2  dir        Real              Imag')
     do ipert1= 1, natom
       do idir1= 1, 3
         irow=( ipert1-1)*3 + idir1
         do ipert2= 1, natom
           do idir2= 1, 3
             icol=( ipert2-1)*3 + idir2
             val(:)=blkval(:,idir1,ipert1,idir2,ipert2,idir3,ipert3,kblok)
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
   call ddb_lw%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega, &
 & rffreq=rffreq)
   if (iblok/=0.and.iblok==kblok) then
     if (opt==1) then
       call wrtout([ab_out,std_out], ' Frozen-spin Berry curvature of Born effective charges')
     else if (opt==2) then
       call wrtout([ab_out,std_out], ' Relaxed-spin Berry curvature of Born effective charges')
     end if
     call wrtout([ab_out,std_out], ' E-dir      atom   dir        Real              Imag')
     ipert1= natom + 2
     do idir1= 1, 3
       do ipert2= 1, natom
         do idir2= 1, 3
           val(:)=blkval(:,idir1,ipert1,idir2,ipert2,idir3,ipert3,kblok)
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
   call ddb_lw%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega, &
 & rffreq=rffreq)
   if (iblok/=0.and.iblok==kblok) then
     if (opt==1) then
       call wrtout([ab_out,std_out], ' Frozen-spin Berry curvature of clamped-ion dielectric tensor')
     else if (opt==2) then
       call wrtout([ab_out,std_out], ' Relaxed-spin Berry curvature of clamped-ion dielectric tensor')
     end if
     call wrtout([ab_out,std_out], '  dir  dir        Real              Imag')
     ipert1= ddb_lw%natom + 2
     ipert2= ddb_lw%natom + 2
     do idir2= 1, 3
       do idir1= 1, 3
         val(:)=blkval(:,idir1,ipert1,idir2,ipert2,idir3,ipert3,kblok)
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
   call ddb_lw%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega, &
 & rfmagn=rfmagn,rffreq=rffreq)
   if (iblok/=0.and.iblok==kblok) then
     if (opt==1) then
       call wrtout([ab_out,std_out], ' Frozen-spin Berry curvature of clamped-ion magnetoelectric susceptibility')
     else if (opt==2) then
       call wrtout([ab_out,std_out], ' Relaxed-spin Berry curvature of clamped-ion magnetoelectric susceptibility')
     end if
     call wrtout([ab_out,std_out], ' M-dir E-dir        Real              Imag')
     ipert1= ddb_lw%natom + 5
     ipert2= ddb_lw%natom + 2
     do idir2= 1, 3
       do idir1= 1, 3
         val(:)=blkval(:,idir1,ipert1,idir2,ipert2,idir3,ipert3,kblok)/ucvol
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
 call ddb_lw%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega, &
& rfmagn=rfmagn,rffreq=rffreq)
 if (iblok/=0.and.iblok==kblok) then
   if (opt==1) then
     call wrtout([ab_out,std_out], ' Frozen-spin Berry curvature of clamped-ion magnetic susceptibility')
   else if (opt==2) then
     call wrtout([ab_out,std_out], ' Relaxed-spin Berry curvature of clamped-ion magnetic susceptibility')
   end if
   call wrtout([ab_out,std_out], '  dir  dir        Real              Imag')
   ipert1= ddb_lw%natom + 5
   ipert2= ddb_lw%natom + 5
   do idir2= 1, 3
     do idir1= 1, 3
       val(:)=blkval(:,idir1,ipert1,idir2,ipert2,idir3,ipert3,kblok)/ucvol
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
 call ddb_lw%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, omega=omega, &
& rfmagn=rfmagn,rffreq=rffreq)
 if (iblok/=0.and.iblok==kblok) then
   if (opt==1) then
     call wrtout([ab_out,std_out], ' Frozen-spin Berry curvature of magnetic Born effective charges')
   else if (opt==2) then
     call wrtout([ab_out,std_out], ' Relaxed-spin Berry curvature of magnetic Born effective charges')
   end if
   call wrtout([ab_out,std_out], ' atom   dir     B-dir        Real              Imag')

   ipert2= ddb_lw%natom + 5
   do ipert1= 1, natom
     do idir1= 1, 3
       do idir2= 1, 3
         val(:)=blkval(:,idir1,ipert1,idir2,ipert2,idir3,ipert3,kblok)
         write(msg,'(i3,4x,a2,7x,a2,2x,2es18.9)') &
       & ipert1, cart(idir1), cart(idir2), val(1), val(2)
         call wrtout([ab_out,std_out], msg)
       end do
     end do
     call wrtout([ab_out,std_out], ' ')
   end do
 end if

 end subroutine mp_d3etot_print
!!***

end module m_ddb_magpen
!!***
