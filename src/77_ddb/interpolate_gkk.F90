!{\src2tex{textfont=tt}}
!!****f* ABINIT/interpolate_gkk
!!
!! NAME
!! interpolate_gkk
!!
!! FUNCTION
!! This routine interpolates the gkk matrices for all q vectors
!!  between points on the full kpt_phon grid.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2016 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!   kpt_phon = coordinates of all kpoints close to the FS
!!
!! OUTPUT
!!   elph_ds = modified gkq
!!
!! NOTES
!!  inspired to some extent by epcouple.f from the DecAFT package by J. Kay Dewhurst
!!  most inputs taken from mkifc.f
!!  in anaddb set ifcflag 1 such that the IFC are calculated in atmfrc prior to calling elphon
!!
!! PARENTS
!!      get_all_gkk2
!!
!! CHILDREN
!!      ftgkk,ifc_fourq,wrap2_pmhalf,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine interpolate_gkk(crystal,ifc,elph_ds,kpt_phon)

 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_errors

 use m_numeric_tools,   only : wrap2_pmhalf
 use m_crystal,         only : crystal_t
 use m_ifc,             only : ifc_type, ifc_fourq

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'interpolate_gkk'
 use interfaces_77_ddb, except_this_one => interpolate_gkk
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: crystal
 type(ifc_type),intent(in) :: ifc
 type(elph_type),intent(inout) :: elph_ds
!arrays
 real(dp),intent(in) :: kpt_phon(3,elph_ds%k_phon%nkpt)

!Local variables-------------------------------
  ! output variables for dfpt_phfrq
! variables for zhpev
! variables for phonon interpolation
!scalars
 integer :: i1,i2,ikpt_phon2,iFSqpt,ib1,ib2,ier,ii
 integer :: iost,isppol,qtor,natom
 integer :: sz1,sz2,sz3,sz4,unit_gkkp
 real(dp) :: qphnrm,res
 !character(len=500) :: msg
!arrays
 real(dp) :: gprim(3,3)
 real(dp) :: displ(2,elph_ds%nbranch,elph_ds%nbranch),eigval(3*crystal%natom)
 real(dp) :: eigvec(3*3*crystal%natom*3*crystal%natom)
 real(dp) :: pheigvec(2*elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: phfrq_tmp(elph_ds%nbranch),qphon(3),redkpt(3)
 real(dp),allocatable :: gkk2_diag_tmp(:,:,:,:),gkk2_tmp(:,:,:,:,:,:,:)
 real(dp),allocatable :: matrx(:,:),zhpev1(:,:)
 real(dp),allocatable :: zhpev2(:)

! *************************************************************************

!
!NOTE: mjv 18/5/2008 reverted to old style of ftgkk with all kpt done together.
!may want to modify this later to use the new cleaner format with 1 FT at a
!time.
!
 write(std_out,*) 'interpolate_gkk : enter'

 natom = crystal%natom
 gprim = ifc%gprim

 if (elph_ds%nsppol /= 1) then
   MSG_ERROR("interpolate_gkk not coded with nsppol>1 yet")
 end if
 isppol = 1


!------------------------------------------------------
!complete dynamical matrices for all qpts between points
!on full kpt grid (interpolation from IFC)
!------------------------------------------------------

 sz1=elph_ds%ngkkband;sz2=elph_ds%nbranch
 sz3=elph_ds%k_phon%nkpt;sz4=elph_ds%nFSband
!allocate (gkk_tmp(2,sz1,sz1,sz2,sz2,1,1))
!DEBUG
!allocate (gkk_tmp_full(2,sz1,sz1,sz2,elph_ds%nFSband,sz3))
!allocate (gkk_tmp_full(2,s2,sz4,sz4,sz3))
!ENDDEBUG
 ABI_ALLOCATE(gkk2_tmp,(2,sz1,sz1,sz2,sz2,sz3,1))
 ABI_ALLOCATE(gkk2_diag_tmp,(sz1,sz1,sz2,sz3))
 ABI_ALLOCATE(zhpev1,(2,2*3*natom-1))
 ABI_ALLOCATE(zhpev2,(3*3*natom-2))
 ABI_ALLOCATE(matrx,(2,(3*natom*(3*natom+1))/2))

 qphnrm = one
!in this part use the inverse Fourier transform to get 1 (arbitrary) qpt at a
!time
 ii = 0
 qtor = 0
 unit_gkkp = 150
 open (unit=unit_gkkp,file='gkkp_file_ascii',form='formatted',status='unknown',iostat=iost)
 if (iost /= 0) then
   MSG_ERROR("error opening gkkpfile as new")
 end if

!loop over all FS pairs.
!do ikpt1=1,elph_ds%k_phon%nkptirr
!do iFSqpt=1,elph_ds%k_phon%nkpt

!
!this should run through the sparse mesh of 2x2x2 kpoints
!
 do iFSqpt=1,elph_ds%k_phon%nkpt
   res = 2.0_dp*(kpt_phon(1,iFSqpt)+one)
   if (abs(res-int(res)) > tol10) cycle
   res = 2.0_dp*(kpt_phon(2,iFSqpt)+one)
   if (abs(res-int(res)) > tol10) cycle
   res = 2.0_dp*(kpt_phon(3,iFSqpt)+one)
   if (abs(res-int(res)) > tol10) cycle

!  do ikpt1=1,1
!  
!  NOTE: should be very easy to parallelize!
!  
!  write(std_out,*) ' interpolate_gkk : ikpt1 = ',ikpt1, ' / ', elph_ds%k_phon%nkptirr
   write(std_out,*) ' interpolate_gkk : ikpt1 = ',iFSqpt, ' / ', elph_ds%k_phon%nkpt

!  DEBUG
!  write(std_out,*) ' interpolate_gkk : Warning debug version'
!  cycle
!  ENDDEBUG

   gkk2_tmp(:,:,:,:,:,:,:) = zero

!  qphon = 1 - 2    ie.  1 = 2+qphon
   qphon(:) = kpt_phon(:,iFSqpt)

!  shouldnt be necessary here, but oh well
   call wrap2_pmhalf(qphon(1),redkpt(1),res)
   call wrap2_pmhalf(qphon(2),redkpt(2),res)
   call wrap2_pmhalf(qphon(3),redkpt(3),res)

   qphon(:) = redkpt(:)
   redkpt(1) = qphon(1)*gprim(1,1)+qphon(2)*gprim(1,2)+qphon(3)*gprim(1,3)
   redkpt(2) = qphon(1)*gprim(2,1)+qphon(2)*gprim(2,2)+qphon(3)*gprim(2,3)
   redkpt(3) = qphon(1)*gprim(3,1)+qphon(2)*gprim(3,2)+qphon(3)*gprim(3,3)
   write (unit_gkkp,*) 'qp= ', redkpt

   call ifc_fourq(ifc,crystal,qphon,phfrq_tmp,displ,out_eigvec=pheigvec)
   write (unit_gkkp,*) phfrq_tmp(:)*Ha_cmm1

   ii = ii+1
!  if(ii > 0 .and. ii < 1000) write(std_out,'(a,i5,3E16.6,2x)') &
!  &   ' wrote phfrq_tmp for time ', ii, phfrq_tmp
!  end if

!  phonon eigenvectors are in eigvec
!  real and imaginary parts
!  phonon displacements = eigvec/sqrt(M_i) are in displ
!  real and imaginary parts

!  DEBUG
!  test: uniform phonon frequency
!  phfrq_tmp(:) = 0.0001_dp
!  ENDDEBUG

!  FT gamma matrices for all kpt_phon points, and
!  for qpoint = qphon(:) = kpt_phon(ikpt_phon)

   call ftgkk(ifc%wghatm,gkk2_tmp,elph_ds%gkk_rpt,elph_ds%gkqwrite,&
&   elph_ds%gkk_rptwrite,gprim,1,&
&   natom,elph_ds%k_phon%nkpt,elph_ds%ngkkband,elph_ds%k_phon%nkpt,1,ifc%nrpt,elph_ds%nsppol,&
&   qtor,ifc%rpt,qphon,elph_ds%unit_gkk_rpt,elph_ds%unitgkq)

!  NOTE: Normally the eigenvectors of the gkk2_tmp should be the same as eigvec

!  Diagonalize gamma matrices at qpoint (complex matrix) for all kpt_phon.
!  Copied from dfpt_phfrq
   do ikpt_phon2=1,elph_ds%k_phon%nkpt
     res = 8.0_dp*(kpt_phon(1,ikpt_phon2)+one)
     if (abs(res-int(res)) > tol10) cycle
     res = 8.0_dp*(kpt_phon(2,ikpt_phon2)+one)
     if (abs(res-int(res)) > tol10) cycle
     res = 8.0_dp*(kpt_phon(3,ikpt_phon2)+one)
     if (abs(res-int(res)) > tol10) cycle

     write (unit_gkkp,*) 'kp= ', kpt_phon(:,ikpt_phon2)

     do ib1=1,elph_ds%ngkkband
       do ib2=1,elph_ds%ngkkband
         ier=0
         ii=1
         do i2=1,3*natom
           do i1=1,i2
             matrx(1,ii)=gkk2_tmp(1,ib1,ib2,i1,i2,ikpt_phon2,1)
             matrx(2,ii)=gkk2_tmp(2,ib1,ib2,i1,i2,ikpt_phon2,1)
             ii=ii+1
           end do
         end do
         call ZHPEV ('N','U',3*natom,matrx,eigval,eigvec,3*natom,zhpev1,&
&         zhpev2,ier)

         gkk2_diag_tmp(ib2,ib1,:,ikpt_phon2) = eigval(:)
         do i1=1,3*natom
           write (unit_gkkp,*) elph_ds%minFSband-1+ib1,elph_ds%minFSband-1+ib2,i1,&
&           eigval(i1)
         end do
       end do
     end do
   end do

   if (elph_ds%gkk2write == 1) then
     write(std_out,*) 'WARNING COMMENTED WRITE TO BINARY FILE!!!'
!    write (elph_ds%unit_gkk2,REC=iFSqpt) gkk2_diag_tmp(:,:,:,:)
     write(std_out,'(a,i4,4(2E16.6,2x))') ' gkk2 loop ', &
&     iFSqpt,gkk2_diag_tmp(1,1,:,1:2),gkk2_diag_tmp(1,1,:,elph_ds%k_phon%nkpt-1:elph_ds%k_phon%nkpt)
!    &    ikpt1,gkk2_tmp(:,1,1,1,1,1:2),gkk2_tmp(:,1,1,elph_ds%k_phon%nkpt-1:elph_ds%k_phon%nkpt)
   else if (elph_ds%gkk2write == 0) then
     elph_ds%gkk2(:,:,:,:,iFSqpt,isppol) = gkk2_diag_tmp(:,:,:,:)
!    elph_ds%gkk2(:,:,:,:,ikpt1) = gkk2_tmp
     write(std_out,*) ' interpolate_gkk : gkk2(b=1,b=1,:,kpt=1,iFSqpt) = '
     write(std_out,*) gkk2_diag_tmp(1,1,:,1)
   end if

 end do
!end do on iFSqpt

 ABI_DEALLOCATE(matrx)
 ABI_DEALLOCATE(zhpev1)
 ABI_DEALLOCATE(zhpev2)

end subroutine interpolate_gkk
!!***
