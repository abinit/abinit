!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_npert_rbz
!! NAME
!! get_npert_rbz
!!
!! FUNCTION
!! Get the number of effective pertubation done in looper3, nkpt_rbz, nband_rbz
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!
!! OUTPUT
!!  npert=number of effective pertubation done in looper3
!!  nkpt_rbz= nkpt in the reduced brillouin zone
!!  nband_rbz= nband in the reduced brillouin zone
!!
!! PARENTS
!!      finddistrproc,initmpi_pert,mpi_setup
!!
!! CHILDREN
!!      irreducible_set_pert,littlegroup_pert,littlegroup_q,mati3inv,metric
!!      mkrdim,symatm,symkpt,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine get_npert_rbz(dtset,nband_rbz,nkpt_rbz,npert)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_npert_rbz'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: npert
!arrays
 integer,pointer :: nkpt_rbz(:)
 real(dp),pointer :: nband_rbz(:,:)
 type(dataset_type),intent(in) :: dtset

!Local variables-------------------------------
!scalars
 integer :: icase,idir,ikpt,ikpt1,ipert,isppol,isym,maxidir,mpert,nband_k,nsym1,timrev,timrev_pert
 integer :: to_compute_this_pert
 real(dp) :: tolsym8,ucvol
 character(len=500) :: message
!arrays
 integer :: rfdir(9),rf2dir(9),rf2_dir1(3),rf2_dir2(3)
 integer,allocatable :: indkpt1(:,:),indsym(:,:,:),pertsy(:,:),rfpert(:),symq(:,:,:),symrec(:,:,:)
 integer, allocatable :: pert_tmp(:,:), pert_calc(:,:)
 integer,allocatable :: symaf1(:),symrc1(:,:,:),symrl1(:,:,:),symrl1_tmp(:,:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: tnons1_tmp(:,:),wtk_folded(:)

! *************************************************************************

!Define the set of admitted perturbations
 mpert=dtset%natom+6
 if(dtset%natom+10/=0.or.dtset%natom+11/=0) mpert=dtset%natom+11

 ABI_ALLOCATE(symrec,(3,3,dtset%nsym))
!Get the symmetry matrices in terms of reciprocal basis
 do isym=1,dtset%nsym
   call mati3inv(dtset%symrel(:,:,isym),symrec(:,:,isym))
 end do

 ABI_ALLOCATE(indsym,(4,dtset%nsym,dtset%natom))
!Obtain a list of rotated atom labels:
 tolsym8=tol8
 call symatm(indsym,dtset%natom,dtset%nsym,symrec,dtset%tnons,tolsym8,dtset%typat,dtset%xred_orig)

 ABI_ALLOCATE(symq,(4,2,dtset%nsym))
 timrev=1
 call littlegroup_q(dtset%nsym,dtset%qptn,symq,symrec,dtset%symafm,timrev)

!Initialize the list of perturbations rfpert
 ABI_ALLOCATE(rfpert,(mpert))
 rfpert(:)=0
 if(dtset%rfphon==1)rfpert(dtset%rfatpol(1):dtset%rfatpol(2))=1

 if(dtset%rfddk==1)rfpert(dtset%natom+1)=1
 if(dtset%rfddk==2)rfpert(dtset%natom+6)=1

 if(dtset%rf2_dkdk/=0) rfpert(dtset%natom+10)=1
 if(dtset%rf2_dkde/=0) rfpert(dtset%natom+11)=1

 if(dtset%rfelfd==1.or.dtset%rfelfd==2)rfpert(dtset%natom+1)=1
 if(dtset%rfelfd==1.or.dtset%rfelfd==3)rfpert(dtset%natom+2)=1

 if(dtset%rfstrs==1.or.dtset%rfstrs==3)rfpert(dtset%natom+3)=1
 if(dtset%rfstrs==2.or.dtset%rfstrs==3)rfpert(dtset%natom+4)=1

 if(dtset%rfuser==1.or.dtset%rfuser==3)rfpert(dtset%natom+6)=1
 if(dtset%rfuser==2.or.dtset%rfuser==3)rfpert(dtset%natom+7)=1

 if(dtset%rfmagn==1) rfpert(dtset%natom+5)=1 

 ABI_ALLOCATE(pertsy,(3,mpert))
 call irreducible_set_pert(indsym,mpert,dtset%natom,dtset%nsym,pertsy,dtset%rfdir,rfpert,symq,symrec,dtset%symrel)
 npert=0
! ABI_ALLOCATE(pert_tmp,(3*mpert))

! do ipert=1,mpert
!   do idir=1,3
!     if( rfpert(ipert)==1 .and. dtset%rfdir(idir) == 1 )then
!       if (pertsy(idir,ipert)==1.or.&
!&       (dtset%prepanl == 1.and.ipert == dtset%natom+2).or.&
!&       (dtset%prepgkk == 1.and.ipert <= dtset%natom)  ) then
!         npert = npert+1;
!         pert_tmp(npert) = idir+(ipert-1)*3;
!       else
!         write(message, '(a,a,i0,a,i0,a,a,a,a,a,a)' )ch10,&
!&         'The perturbation idir=',idir,'  ipert=',ipert,' is',ch10,&
!&         'symmetric of a previously calculated perturbation.',ch10,&
!&         'So, its SCF calculation is not needed.',ch10
!         call wrtout(std_out,message,'COLL')
!       end if ! Test of existence of symmetry of perturbation
!     end if ! Test of existence of perturbation
!   end do
! end do

!Initialize rf2dir :
 rf2_dir1(1:3)=dtset%rf2_pert1_dir(1:3)
 rf2_dir2(1:3)=dtset%rf2_pert2_dir(1:3)
!Diagonal terms :
 rf2dir(1) = rf2_dir1(1)*rf2_dir2(1)
 rf2dir(2) = rf2_dir1(2)*rf2_dir2(2)
 rf2dir(3) = rf2_dir1(3)*rf2_dir2(3)
!Upper triangular terms :
 rf2dir(4) = rf2_dir1(2)*rf2_dir2(3)
 rf2dir(5) = rf2_dir1(1)*rf2_dir2(3)
 rf2dir(6) = rf2_dir1(1)*rf2_dir2(2)
!Lower triangular terms :
 rf2dir(7) = rf2_dir1(3)*rf2_dir2(2)
 rf2dir(8) = rf2_dir1(3)*rf2_dir2(1)
 rf2dir(9) = rf2_dir1(2)*rf2_dir2(1)

!Determine existence of pertubations and of pertubation symmetries
!Create array with pertubations which have to be calculated
 ABI_ALLOCATE(pert_tmp,(2,3*(dtset%natom+6)+18))

 do ipert=1,mpert
   if (ipert<dtset%natom+10) then
     maxidir = 3
     rfdir(1:3) = dtset%rfdir(:)
     rfdir(4:9) = 0
   else
     maxidir = 9
     rfdir(1:9) = rf2dir(:)
   end if
   do idir=1,maxidir
     if( rfpert(ipert)==1 .and. rfdir(idir) == 1 )then
       to_compute_this_pert = 0
       if (ipert>=dtset%natom+10) then
         to_compute_this_pert = 1
       else if ((pertsy(idir,ipert)==1).or.&
&         ((dtset%prepanl == 1).and.(ipert == dtset%natom+2)).or.&
&         ((dtset%prepgkk == 1).and.(ipert <= dtset%natom))  ) then
         to_compute_this_pert = 1
       end if
       if (to_compute_this_pert /= 0) then
         npert = npert+1;
         pert_tmp(1,npert) = ipert
         pert_tmp(2,npert) = idir
       else
         write(message, '(a,a,i4,a,i4,a,a,a,a,a,a)' )ch10,&
&         ' The perturbation idir=',idir,'  ipert=',ipert,' is',ch10,&
&         ' symmetric of a previously calculated perturbation.',ch10,&
&         ' So, its SCF calculation is not needed.',ch10
         call wrtout(std_out,message,'COLL')
       end if ! Test of existence of symmetry of perturbation
     end if ! Test of existence of perturbation
   end do
 end do
 ABI_ALLOCATE(pert_calc,(2,npert))
 do icase=1,npert
   pert_calc(:,icase)=pert_tmp(:,icase)
 end do
 ABI_DEALLOCATE(pert_tmp)
 ABI_DEALLOCATE(pertsy)
 ABI_DEALLOCATE(rfpert)

! Write YAML doc with the list of irreducible perturbations. Example.
!
!--- !IrredPerts
!# List of irreducible perturbations
!irred_perts:
!  - qpt: [ 0.0000000000000000,  0.0000000000000000,  0.0000000000000000]
!    ipert : 1
!    idir  : 1
!  - qpt: [ 0.0000000000000000,  0.0000000000000000,  0.0000000000000000]
!    ipert : 2
!    idir  : 1
!..
 write(std_out,'(a)')"--- !IrredPerts"
 write(std_out,'(a)')'# List of irreducible perturbations'
 write(std_out,'(a)')'irred_perts:'

 do icase=1,npert
!   pert = pert_tmp(icase)

!   if (pert <= dtset%natom*3) then
!     idir = mod(pert, 3)
!     if (idir==0) idir=3
!     ipert=((pert-idir) / 3 + 1)
!   else
!     idir = mod(pert, 3)
!     if (idir==0) idir=3
!     ipert = dtset%natom + ((pert - 3*dtset%natom - 1) / 3) + 1
!   end if
   ipert = pert_calc(1,icase)
   idir = pert_calc(2,icase)

   write(std_out,'(a,3(f20.16,a))')"   - qpt: [ ",dtset%qptn(1),", ", dtset%qptn(2),", ", dtset%qptn(3),"]"
   write(std_out,'(a,i0)')"     ipert: ",ipert
   write(std_out,'(a,i0)')"     idir: ",idir
 end do

 write(std_out,'(a)')"..."

! ABI_ALLOCATE(pert_calc,(npert))
! do icase=1,npert
!   pert_calc(icase) = pert_tmp(icase)
! end do

 call mkrdim(dtset%acell_orig(1:3,1),dtset%rprim_orig(1:3,1:3,1),rprimd)
 call metric(gmet,gprimd,std_out,rmet,rprimd,ucvol)

 ABI_ALLOCATE(nkpt_rbz,(npert))
 ABI_ALLOCATE(indkpt1,(dtset%nkpt,npert))
 indkpt1=0

 do icase=1,npert
!   if (pert_calc(icase) <= dtset%natom*3) then
!     idir = mod(pert_calc(icase),3)
!     if (idir==0) idir=3
!     ipert=( (pert_calc(icase)-idir) / 3 + 1)
!   else
!     ipert = dtset%natom + ((pert_calc(icase) - 3*dtset%natom - 1) / 3) + 1
!     idir = mod(pert_calc(icase),3)
!     if (idir==0) idir=3
!   end if
   ipert = pert_calc(1,icase)
   idir = pert_calc(2,icase)

   ABI_ALLOCATE(symrl1_tmp,(3,3,dtset%nsym))
   ABI_ALLOCATE(symaf1,(dtset%nsym))
   ABI_ALLOCATE(tnons1_tmp,(3,dtset%nsym))
!  MJV TODO: check whether prepgkk should be used here
   if (dtset%prepanl /= 1 .and. dtset%berryopt /=4 .and. dtset%berryopt /=6 .and. dtset%berryopt /=7 .and. &
&   dtset%berryopt /=14 .and. dtset%berryopt /=16 .and. dtset%berryopt /=17) then   !!HONG
     call littlegroup_pert(gprimd,idir,indsym,std_out,ipert,dtset%natom,dtset%nsym,nsym1,2,&
&     dtset%symafm,symaf1,symq,symrec,&
&     dtset%symrel,symrl1_tmp,0,dtset%tnons,tnons1_tmp)
   else
     nsym1 = 1
   end if
   ABI_DEALLOCATE(tnons1_tmp)
   ABI_DEALLOCATE(symaf1)

   ABI_ALLOCATE(symrc1,(3,3,nsym1))
   ABI_ALLOCATE(symrl1,(3,3,nsym1))
   symrl1(:,:,1:nsym1)=symrl1_tmp(:,:,1:nsym1)
   ABI_DEALLOCATE(symrl1_tmp)
   do isym=1,nsym1
     call mati3inv(symrl1(:,:,isym),symrc1(:,:,isym))
   end do
   ABI_DEALLOCATE(symrl1)

   ABI_ALLOCATE(wtk_folded,(dtset%nkpt))
   timrev_pert=timrev
   if(dtset%ieig2rf>0) then
     call symkpt(0,gmet,indkpt1(:,icase),std_out,dtset%kptns,dtset%nkpt,nkpt_rbz(icase),&
&     1,symrc1,0,dtset%wtk,wtk_folded)
   else
!    For the time being, the time reversal symmetry is not used
!    for ddk, elfd, mgfd perturbations.
     if(ipert==dtset%natom+1 .or. ipert==dtset%natom+10 .or. ipert==dtset%natom+11 .or. &
&     ipert==dtset%natom+2 .or. dtset%berryopt==4 .or. dtset%berryopt==6 .or. dtset%berryopt==7  &
&     .or. dtset%berryopt==14 .or. dtset%berryopt==16 .or. dtset%berryopt==17 )timrev_pert=0  !!HONG
     call symkpt(0,gmet,indkpt1(:,icase),std_out,dtset%kptns,dtset%nkpt,nkpt_rbz(icase),&
     nsym1,symrc1,timrev_pert,dtset%wtk,wtk_folded)
   end if
   ABI_DEALLOCATE(wtk_folded)
   ABI_DEALLOCATE(symrc1)
 end do

 ABI_ALLOCATE(nband_rbz,(maxval(nkpt_rbz)*dtset%nsppol,npert))
 nband_rbz=zero
 do icase=1,npert
   do isppol=1,dtset%nsppol
     ikpt1=1
     do ikpt=1,dtset%nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
!      Must test against ikpt1/=nkpt_rbz+1, before evaluate indkpt1(ikpt1)
       if(ikpt1/=nkpt_rbz(icase)+1)then
         if(ikpt==indkpt1(ikpt1,icase))then
           nband_rbz(ikpt1+(isppol-1)*nkpt_rbz(icase),icase)=nband_k
         end if
       end if
     end do
   end do

 end do

 ABI_DEALLOCATE(indkpt1)
 ABI_DEALLOCATE(symq)
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(indsym)
 ABI_DEALLOCATE(pert_calc)

end subroutine get_npert_rbz
!!***
