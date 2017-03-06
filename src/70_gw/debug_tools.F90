!{\src2tex{textfont=tt}}
!!****f* ABINIT/check_zarot
!! NAME
!! check_zarot
!!
!! FUNCTION
!!  Debugging routine used to test zarot.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2017 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3,get_bz_item,getcprj,kdata_free,kdata_init
!!      mkkpg,paw_symcprj,pawcprj_alloc,pawcprj_free,wfd_get_cprj,wfd_sym_ur
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine check_zarot(npwvec,Cryst,ngfft,gvec,psps,pawang,grottb,grottbm1)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 
 use m_fft_mesh,     only : rotate_FFT_mesh
 use m_geometry,     only : normv
 use m_crystal,      only : crystal_t
 use m_paw_sphharm,      only : initylmr
 use m_mpinfo,       only : destroy_mpi_enreg
 use m_pawang,       only : pawang_type
 use m_pawtab,       only : pawtab_type
 use m_pawcprj,      only : pawcprj_type, pawcprj_alloc, pawcprj_free

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'check_zarot'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwvec
 type(crystal_t),intent(in) :: Cryst
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: grottb(npwvec,Cryst%timrev,Cryst%nsym),grottbm1(npwvec,Cryst%timrev,Cryst%nsym)
 integer,intent(in) :: gvec(3,npwvec)

!Local variables-------------------------------
!scalars
 integer :: aa,ig,ig_sym,iginv_sym,ii,ilpa,ilpm,isym,itim,jj,ll,lmax,mm,mqmem_
 integer :: nqpt_,optder,option,normchoice,npts,ix,iy,iz,ir_sym,ir !,nx,ny,nz
 real(dp) :: err,max_diff,test,tmp,ylm_sym,rx,ry,rz
 logical :: found !,iscompatibleFFT
 character(len=500) :: message
 type(MPI_type) :: Fake_MPI_enreg
!arrays
 integer :: toinv(Cryst%nsym),trial(3,3),rm1(3,3)
 integer,allocatable :: nband(:),npwarr(:),irottb(:,:)
 real(dp),allocatable :: DS_mmpl(:,:,:),DSinv_mmpl(:,:,:),qptns(:,:),ylm_q(:,:)
 real(dp),allocatable :: ylmgr_q(:,:,:)
 real(dp),allocatable :: ylmr(:,:),ylmr_gr(:,:,:),nrm(:),rr(:,:) !,dum_tnons(:,:)
 real(dp) :: search(3)

! *************************************************************************

 write(message,'(a)')' check_zarot  : enter '
 call wrtout(std_out,message,'COLL')

 do jj=1,Cryst%nsym
  found=.FALSE.
  do ii=1,Cryst%nsym
   call mati3inv(Cryst%symrec(:,:,ii),trial)
   trial=transpose(trial)
   if (ALL(trial==Cryst%symrec(:,:,jj))) then
    toinv(jj)=ii
    found=.TRUE.
    exit
   end if
  end do
  if (.not. found) then 
    MSG_ERROR("inverse not found!")
  end if
 end do

 mqmem_=1 ; nqpt_=1 ; optder=0
 ABI_MALLOC(npwarr,(mqmem_))
 ABI_MALLOC(qptns,(3,mqmem_))
 npwarr(:)=npwvec ; qptns(:,:)=zero

 lmax=psps%mpsang-1
 write(std_out,*)'lmax= ',lmax
 ABI_MALLOC(ylm_q,(npwvec*mqmem_,(lmax+1)**2))
 ABI_MALLOC(ylmgr_q,(npwvec*mqmem_,3+6*(optder/2),(lmax+1)**2))
 call initmpi_seq(Fake_MPI_enreg)
 ABI_MALLOC(nband,(1))
 nband=0

 ! Note: dtset%nband and dtset%nsppol are not used in sequential mode
 call initylmg(Cryst%gprimd,gvec,qptns,mqmem_,Fake_MPI_enreg,Psps%mpsang,npwvec,nband,nqpt_,npwarr,0,optder,&
& Cryst%rprimd,ylm_q,ylmgr_q)

 call destroy_mpi_enreg(Fake_MPI_enreg)

 ABI_MALLOC(DS_mmpl,(2*lmax+1,2*lmax+1,lmax+1))
 ABI_MALLOC(DSinv_mmpl,(2*lmax+1,2*lmax+1,lmax+1))
 max_diff=zero ; test=zero

 do ig=1,npwvec
  if (ig==1) cycle

  do isym=1,Cryst%nsym
   do itim=1,Cryst%timrev

    ig_sym=grottb(ig,itim,isym) !index of IS G
    DS_mmpl(:,:,:)=pawang%zarot(:,:,:,isym)

    iginv_sym=grottbm1(ig,itim,isym) !index of (IS)^-1 G
    DSinv_mmpl(:,:,:)=pawang%zarot(:,:,:,toinv(isym))

    do ll=0,lmax
     do mm=1,2*ll+1
      ilpm=1+ll**2+ll+(mm-1-ll)
      ylm_sym=ylm_q(ig_sym,ilpm)     !Ylm(IS   G)
      !ylm_sym=ylm_q(iginv_sym,ilpm) !Ylm(IS^-1G)
      !
      ! here we calculate the symmetric
      tmp=zero
      do aa=1,2*ll+1
       test=MAX(test,ABS(DS_mmpl(aa,mm,ll+1)-DSinv_mmpl(mm,aa,ll+1)))
       ilpa=1+ll**2+ll+(aa-1-ll)
       tmp= tmp+ ylm_q(ig,ilpa)*DS_mmpl(aa,mm,ll+1)
      end do
      if (itim==2) tmp=tmp*(-1)**ll
      err=ABS(tmp-ylm_sym) !Ylm(IS G) = D_am Yma(S) (-1)**l

      if (err > tol6) then
       write(std_out,*)'WARNING check fort 77'
       write(77,'(6(a,i3),a)')' -- ig: ',ig,' igsym: ',ig_sym,' isym ',isym,' itim:',itim,' ll: ',ll,' mm: ',(mm-1-ll)," --"
       write(77,*)tmp,ylm_sym,ABS(tmp-ylm_sym)
      end if
      max_diff=MAX(max_diff,err)

     end do
    end do !itim

   end do  !isym
  end do !sym
 end do !ig

 write(std_out,*)"MAX DIFF ",max_diff
 write(std_out,*)"MAX TEST ",test


 ABI_FREE(nband)
 ABI_FREE(npwarr)
 ABI_FREE(qptns)
 ABI_FREE(ylm_q)
 ABI_FREE(ylmgr_q)

 npts = PRODUCT(ngfft(1:3))

 ABI_MALLOC(irottb,(npts,Cryst%nsym))
 !allocate(dum_tnons(3,Cryst%nsym)); dum_tnons=zero
 !call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,dum_tnons,ngfft,irottb,iscompatibleFFT)
 !if (.not.iscompatibleFFT) then
 !  MSG_ERROR("Uncompatible FFT mesh")
 !end if
 !deallocate(dum_tnons)

 ABI_MALLOC(rr,(3,npts))
 ABI_MALLOC(nrm,(npts))
 ii = 0
 do iz=0,ngfft(3)-1
   do iy=0,ngfft(2)-1
     do ix=0,ngfft(1)-1
       ii = ii + 1
       if (ix <= ngfft(1)/2) then
         rx = DBLE(ix)/ngfft(1)
       else
         rx = DBLE(ix-ngfft(1))/ngfft(1)
       end if
       if (iy <= ngfft(2)/2) then
         ry = DBLE(iy)/ngfft(2)
       else
         ry = DBLE(iy-ngfft(2))/ngfft(2)
       end if
       if (iz <= ngfft(3)/2) then
         rz = DBLE(iz)/ngfft(3)
       else
         rz = DBLE(iz-ngfft(3))/ngfft(3)
       end if
       rr(:,ii) = (/rx,ry,rz/)
       nrm(ii) = normv(rr(:,ii),Cryst%rmet,"R")
     end do
   end do
 end do

 irottb = HUGE(0)
 do isym=1,Cryst%nsym
   call mati3inv(Cryst%symrel(:,:,isym),rm1)
   rm1 = transpose(rm1)
   do ii=1,npts
     search = MATMUL(rm1,rr(:,ii))
     do jj=1,npts
       if (ALL (ABS(search-rr(:,jj)) < tol6)) irottb(ii,isym) = jj
     end do
   end do
 end do

 option=1; normchoice=1
 ABI_MALLOC(ylmr,(Psps%mpsang**2,npts))
 ABI_MALLOC(ylmr_gr,(3*(option/2)+6*(option/3),Psps%mpsang**2,npts))

 call initylmr(Psps%mpsang,normchoice,npts,nrm,option,rr,ylmr,ylmr_gr)

 max_diff=zero ; test=zero

 do isym=1,Cryst%nsym
   do ir=1,npts
     ir_sym = irottb(ir,isym) ! idx of R^{-1} (r-\tau)
     if (ir_sym == HUGE(0)) then
       write(std_out,*)"Got HUGE"
       CYCLE
     end if

     do ll=0,lmax
       do mm=1,2*ll+1
         ilpm=1+ll**2+ll+(mm-1-ll)
         ylm_sym=ylmr(ir_sym,ilpm)      !Ylm(R^{-1}(r-t))
         !ylm_sym=ylm_q(iginv_sym,ilpm) !Ylm(IS^-1G)
         !
         ! here we calculate the symmetric
         tmp=zero
         do aa=1,2*ll+1
           test=MAX(test,ABS(DS_mmpl(aa,mm,ll+1)-DSinv_mmpl(mm,aa,ll+1)))
           ilpa=1+ll**2+ll+(aa-1-ll)
           tmp= tmp+ ylmr(ir,ilpa)*DS_mmpl(aa,mm,ll+1)
         end do
         !if (itim==2) tmp=tmp*(-1)**ll
         err=ABS(tmp-ylm_sym) ! Ylm(R^{1}(r-t)) = D_am Yma(r)

         if (err > tol6) then
           write(std_out,*)'WARNING check fort 78'
           write(77,'(5(a,i3),a)')' -- ir: ',ir,' ir_sym: ',ir_sym,' isym ',isym,' ll: ',ll,' mm: ',(mm-1-ll)," --"
           write(77,*)tmp,ylm_sym,ABS(tmp-ylm_sym)
         end if
         max_diff=MAX(max_diff,err)

       end do ! ll
     end do ! mm

   end do ! ir
 end do ! isym

 write(std_out,*)"MAX DIFF REAL SPACE ",max_diff
 write(std_out,*)"MAX TEST REAL SPACE ",test

 ABI_FREE(ylmr)
 ABI_FREE(ylmr_gr)
 ABI_FREE(irottb)
 ABI_FREE(rr)
 ABI_FREE(nrm)

 ABI_FREE(DS_mmpl)
 ABI_FREE(DSinv_mmpl)

end subroutine check_zarot
!!***

!!****f* ABINIT/paw_check_symcprj
!! NAME
!! paw_check_symcprj
!!
!! FUNCTION
!!   Test the routines used to symmetrize PAW CPRJ
!!
!! COPYRIGHT
!! Copyright (C) 2010-2017 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      fftbox_execute,fftbox_plan3,get_bz_item,getcprj,kdata_free,kdata_init
!!      mkkpg,paw_symcprj,pawcprj_alloc,pawcprj_free,wfd_get_cprj,wfd_sym_ur
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine paw_check_symcprj(Wfd,ik_bz,band,spin,sym_mode,Cryst,Kmesh,Psps,Pawtab,Pawang,Cprj_bz)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_fft            

 use m_pawang,         only : pawang_type
 use m_pawtab,         only : pawtab_type
 use m_pawcprj,        only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_crystal,        only : crystal_t
 use m_bz_mesh,        only : kmesh_t, get_BZ_item
 use m_wfd,            only : wfd_t, wfd_get_cprj, kdata_init, kdata_free, kdata_t, wfd_sym_ur

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_check_symcprj'
 use interfaces_65_paw
 use interfaces_66_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,band,spin,sym_mode
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfd_t),intent(inout) :: Wfd
!arrays
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 type(pawcprj_type),intent(out) :: Cprj_bz(Cryst%natom,Wfd%nspinor)

!Local variables ------------------------------
!scalars
 integer :: k_sym,k_tim,ik_ibz,ig,fft_idx
 integer :: cpopt,choice,npw_k,istwf_k,nkpg
 integer :: iatom,iatm,isp
 complex(dpc) :: k_eimkt
 logical :: k_isirred
 type(Kdata_t) :: Gdata
 type(fftbox_plan3_t) :: plan
!arrays
 integer :: k_umklp(3)
 real(dp) :: k_bz(3)
 real(dp),allocatable :: kpg_k(:,:),vectin(:,:)
 complex(dpc) :: ur1_dpc(Wfd%nfft*Wfd%nspinor)
 complex(gwpc) :: ur1(Wfd%nfft*Wfd%nspinor)
 type(pawcprj_type),allocatable :: Cprj_srt(:,:)

!************************************************************************

 call get_BZ_item(Kmesh,ik_bz,k_bz,ik_ibz,k_sym,k_tim,k_eimkt,k_umklp,k_isirred)

 if (k_isirred) then  ! Symmetrization is not needed. Retrieve Cprj_ibz from Wfd and return immediately.
   call wfd_get_cprj(Wfd,band,ik_ibz,spin,Cryst,Cprj_bz,sorted=.FALSE.)
   RETURN
 end if

 select case (sym_mode)

 case (1) ! Faster Symmetrization in reciprocal space.

   call wfd_get_cprj(Wfd,band,ik_ibz,spin,Cryst,Cprj_bz,sorted=.FALSE.)
   call paw_symcprj(ik_bz,Wfd%nspinor,1,Cryst,Kmesh,Pawtab,Pawang,Cprj_bz)

 case (2) ! Symmetrize u(r) in reciprocal space, FFT from r to G then call getcprj to obtain the symmetrized cprj.

   ! Symmetrization in real space on the FFT BOX.
   call wfd_sym_ur(Wfd,Cryst,Kmesh,band,ik_bz,spin,ur1)

   istwf_k = 1

   ! Init k_data associated to the G-sphere centered at k_bz.
   call kdata_init(Gdata,Cryst,Psps,k_bz,istwf_k,Wfd%ngfft,Wfd%MPI_enreg,ecut=Wfd%ecut)

   npw_k = Gdata%npw
   !
   ! Compute (k+G) vectors
   nkpg=0
   ABI_MALLOC(kpg_k,(npw_k,nkpg))
   if (nkpg>0) then
     call mkkpg(Gdata%kg_k,kpg_k,k_bz,nkpg,npw_k)
   end if

   ABI_MALLOC(vectin,(2,npw_k*Wfd%nspinor))
   !ABI_CHECK(npw_k==Wfd%npwwfn,"Wrong npw")
   !
   ! FFT R -> G TODO Fix issue with double precision complex.
   ur1_dpc = ur1

   call fftbox_plan3(plan,Wfd%ngfft(1:3),Wfd%ngfft(7),-1)
   call fftbox_execute(plan,ur1_dpc)

   do ig=1,npw_k ! FFT box to G-sphere.
     fft_idx = Gdata%igfft0(ig)
     if (fft_idx/=0) then ! G-G0 belong to the FFT mesh.
       vectin(1,ig) = DBLE (ur1_dpc(fft_idx))
       vectin(2,ig) = AIMAG(ur1_dpc(fft_idx))
     else
       vectin(:,ig) = zero ! Set this component to zero.
     end if
   end do
   !
   ! Calculate SORTED cprj.
   cpopt   = 0 ! Nothing is already calculated.
   choice  = 1

   ABI_DT_MALLOC(Cprj_srt,(Wfd%natom,Wfd%nspinor))
   call pawcprj_alloc(Cprj_srt,0,Wfd%nlmn_sort)

   call getcprj(choice,cpopt,vectin,Cprj_srt,Gdata%fnl_dir0der0,&
&    0,Wfd%indlmn,istwf_k,Gdata%kg_k,kpg_k,k_bz,Wfd%lmnmax,Wfd%mgfft,Wfd%MPI_enreg,&
&    Cryst%natom,Cryst%nattyp,Wfd%ngfft,Wfd%nloalg,npw_k,Wfd%nspinor,Cryst%ntypat,&
&    Gdata%phkxred,Wfd%ph1d,Gdata%ph3d,Cryst%ucvol,1)

   ABI_FREE(vectin)
   ABI_FREE(kpg_k)
   !
   ! Reorder cprj (sorted --> unsorted)
   do iatom=1,Cryst%natom
     iatm=Cryst%atindx(iatom)
     do isp=1,Wfd%nspinor
       Cprj_bz(iatom,isp)%cp=Cprj_srt(iatm,isp)%cp
     end do
   end do
   call pawcprj_free(Cprj_srt)
   ABI_DT_FREE(Cprj_srt)

   call kdata_free(Gdata)

 case default
   MSG_ERROR("Wrong sym_mode")
 end select

end subroutine paw_check_symcprj
!!***

