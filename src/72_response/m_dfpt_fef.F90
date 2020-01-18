!!****m* ABINIT/m_dfpt_fef
!! NAME
!!  m_dfpt_fef
!!
!! FUNCTION
!!  Response calculations in finite electric field.
!!
!! COPYRIGHT
!!  Copyright (C) 2004-2020 ABINIT group (XW).
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

module m_dfpt_fef

 use defs_basis
 use m_abicore
 use m_errors
 use m_efield
 use m_dtset

 use defs_abitypes, only : MPI_type
 use m_kg,        only : kpgio
 use m_cgtools,   only : overlap_g
 use m_mpinfo,    only : proc_distrb_cycle

 implicit none

 private
!!***

 public :: dfptff_initberry
 public :: dfptff_gradberry
 public :: dfptff_gbefd
 public :: dfptff_edie
 public :: dfptff_ebp
 public :: dfptff_die
 public :: dfptff_bec
 public :: qmatrix
!!***

contains
!!***

!!****f* ABINIT/dfptff_initberry
!! NAME
!! dfptff_initberry
!!
!! FUNCTION
!! Initialization of response calculations in finite electric field.
!!
!! INPUTS
!! dtset <type(dataset_type)> = all input variables in this dataset
!! gmet(3,3) = reciprocal space metric tensor in bohr**-2
!! kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!! kg1(3,mpw1*mkmem) = reduced (integer) coordinates of G vecs for response wfs
!! mband =  maximum number of bands
!! mkmem = maximum number of k-points in core memory
!! mpw = maximum number of plane waves
!! mpw1 = maximum number of plane waves for response wavefunctions
!! nkpt = number of k points
!! npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!! npwar1(nkpt) = number of planewaves in basis and boundary for response wfs
!! nsppol = 1 for unpolarized, 2 for spin-polarized
!! occ(mband*nkpt*nsppol) = occup number for each band at each k point
!! rprimd(3,3) = dimensional primitive vectors
!!
!! OUTPUT
!! dtefield=variables related to response Berry-phase calculation
!! pwindall(max(mpw,mpw1)*mkmem,8,3) = array used to compute the overlap matrices
!! pwindall(:,1,:) <- <u^(0)_i|u^(0)_i+1>
!! pwindall(:,2,:) <- <u^(0)_i|u^(0)_i-1>
!! pwindall(:,3,:) <- <u^(1)_i|u^(1)_i+1>
!! pwindall(:,4,:) <- <u^(1)_i|u^(1)_i-1>
!! pwindall(:,5,:) <- <u^(1)_i|u^(0)_i+n+1>
!! pwindall(:,6,:) <- <u^(1)_i|u^(0)_i+n-1>
!! pwindall(:,7,:) <- <u^(0)_i|u^(1)_i-n+1>
!! pwindall(:,8,:) <- <u^(0)_i|u^(1)_i-n-1>
!!
!! NOTES
!!      this duplicates in part initberry for the init of the dtefield - should be made
!!      into a common constructor in m_dtefield or somethin
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      kpgio,wrtout
!!
!! SOURCE

subroutine dfptff_initberry(dtefield,dtset,gmet,kg,kg1,mband,mkmem,mpi_enreg,&
&                mpw,mpw1,nkpt,npwarr,npwar1,nsppol,occ,pwindall,rprimd)

!Arguments ----------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,mpw,mpw1,nkpt,nsppol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(efield_type),intent(inout) :: dtefield !vz_i needs efield2
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),kg1(3,mpw1*mkmem),npwar1(nkpt)
 integer,intent(in) :: npwarr(nkpt)
 integer,intent(out) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
 real(dp),intent(in) :: gmet(3,3),occ(mband*nkpt*nsppol),rprimd(3,3)

!Local variables ----------------------------------
!scalars
 integer :: flag,iband,icg,idir,ifor,ikg,ikg1,ikpt,ikpt1,ikpt2,ikstr
 integer :: index,ipw,isppol,istr,iunmark,jpw,nband_k,mband_occ_k,nkstr,npw_k
 integer :: npw_k1,orig
 real(dp) :: ecut_eff,occ_val
 character(len=500) :: message
!arrays
 integer :: dg(3)
 integer,allocatable :: kg_tmp(:,:),kpt_mark(:),npwarr_tmp(:),npwtot(:)
 real(dp) :: diffk(3),dk(3)
 real(dp),allocatable :: kpt1(:,:)

! *************************************************************************

!-----------------------------------------------------------
!---------------- initilize dtefield //---------------------
!-----------------------------------------------------------

!Initialization of efield_type variables
 dtefield%efield_dot(:) = zero
 dtefield%dkvecs(:,:) = zero
 dtefield%maxnstr = 0    ; dtefield%maxnkstr  = 0
 dtefield%nstr(:) = 0    ; dtefield%nkstr(:) = 0
 ABI_ALLOCATE(dtefield%ikpt_dk,(nkpt,9,3))
 ABI_ALLOCATE(dtefield%cgindex,(nkpt,nsppol*2))
 ABI_ALLOCATE(dtefield%kgindex,(nkpt))
 dtefield%ikpt_dk(:,:,:) = 0
 dtefield%cgindex(:,:) = 0
 dtefield%mband_occ = 0
 ABI_ALLOCATE(dtefield%nband_occ,(nsppol))
 dtefield%nband_occ = 0
 pwindall(:,:,:) = 0

!Compute the number of occupied bands and check that--------
!it is the same for each k-point----------------------------

 occ_val = two/(dtset%nsppol*one)

 index = 0
 do isppol = 1, nsppol
   do ikpt = 1, nkpt

     mband_occ_k = 0
     nband_k = dtset%nband(ikpt)

     do iband = 1, nband_k
       index = index + 1
       if (abs(occ(index) - occ_val) < tol8) mband_occ_k = mband_occ_k + 1
     end do

     if (ikpt > 1) then
       if (dtefield%nband_occ(isppol) /= mband_occ_k) then
         message = ' The number of valence bands is not the same for every k-point for present spin'
         MSG_ERROR(message)
       end if
     else
       dtefield%mband_occ = max(dtefield%mband_occ,mband_occ_k)
       dtefield%nband_occ(isppol) = mband_occ_k
     end if

   end do
 end do

!DEBUG
!write(std_out,*)'dfptff_initberry:nkpt',nkpt
!ENDDEBUG

!Compute the location of zero-order wavefunction --------------
 icg = 0
 do isppol = 1, nsppol
   do ikpt = 1, nkpt

     dtefield%cgindex(ikpt,isppol) = icg
     nband_k = dtset%nband(ikpt)
     npw_k = npwarr(ikpt)
     icg = icg + dtset%nspinor*npw_k*nband_k

   end do
 end do

!Compute the location of kg --------------
 ikg = 0
 do ikpt = 1, nkpt

   dtefield%kgindex(ikpt) = ikg
   npw_k = npwarr(ikpt)
   ikg = ikg + npw_k

 end do

!Compute the location of first-order wavefunction -------------------
 icg = 0
 do isppol = 1, nsppol
   do ikpt = 1, nkpt

     dtefield%cgindex(ikpt,isppol+nsppol) = icg
     nband_k = dtset%nband(ikpt)
     npw_k = npwar1(ikpt)
     icg = icg + dtset%nspinor*npw_k*nband_k

   end do
 end do

!Compute the reciprocal lattice coordinates of the electric field-----

 dtefield%efield_dot(1) = dot_product(dtset%efield(:),rprimd(:,1))
 dtefield%efield_dot(2) = dot_product(dtset%efield(:),rprimd(:,2))
 dtefield%efield_dot(3) = dot_product(dtset%efield(:),rprimd(:,3))

 write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
& ' initberry: Reciprocal lattice coordinates of the electric field',ch10,&
& '  efield_dot(1:3) = ',dtefield%efield_dot(1:3),ch10
 call wrtout(std_out,message,'COLL')

!find the related k points to every k point in full BZ-----------------

!loop over three reciprocal directions
 do idir = 1, 3

   if (dtset%rfdir(idir) == 1) then

!    Compute dk(:), the vector between a k-point and its nearest
!    neighbour along the direction idir

     dk(:) = zero
     dk(idir) = 1000_dp   ! initialize with a big number
     do ikpt = 2, nkpt
       diffk(:) = abs(dtset%kptns(:,ikpt) - dtset%kptns(:,1))
       if ((diffk(1) < dk(1)+tol8).and.(diffk(2) < dk(2)+tol8).and.&
&       (diffk(3) < dk(3)+tol8)) dk(:) = diffk(:)
     end do
     dtefield%dkvecs(:,idir) = dk(:)

!    For each k point, find k_prim such that k_prim= k + dk mod(G)
!    where G is a vector of the reciprocal lattice

     do ikpt = 1, nkpt

!      First: k + dk
       do ikpt1 = 1, nkpt
         diffk(:) = abs(dtset%kptns(:,ikpt1) - dtset%kptns(:,ikpt) - dk(:))
         if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
           dtefield%ikpt_dk(ikpt,1,idir) = ikpt1
           exit
         end if
       end do

!      Second: k - dk
       do ikpt1 = 1, nkpt
         diffk(:) = abs(dtset%kptns(:,ikpt1) - dtset%kptns(:,ikpt) + dk(:))
         if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
           dtefield%ikpt_dk(ikpt,2,idir) = ikpt1
           exit
         end if
       end do

!      new
!      3rd: k + (n+1)dk

       do ikpt1 = 1, nkpt
         diffk(:) = abs(dtset%kptns(:,ikpt1) - dtset%kptns(:,ikpt) - dk(:) - dtset%qptn(:))
         if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
           dtefield%ikpt_dk(ikpt,3,idir) = ikpt1
           exit
         end if
       end do

!      6th: k - (n+1)dk

       do ikpt1 = 1, nkpt
         diffk(:) = abs(dtset%kptns(:,ikpt1) - dtset%kptns(:,ikpt) + dk(:) + dtset%qptn(:))
         if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
           dtefield%ikpt_dk(ikpt,6,idir) = ikpt1
           exit
         end if
       end do

!      4th: k + (n-1)dk

       do ikpt1 = 1, nkpt
         diffk(:) = abs(dtset%kptns(:,ikpt1) - dtset%kptns(:,ikpt) + dk(:) - dtset%qptn(:))
         if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
           dtefield%ikpt_dk(ikpt,4,idir) = ikpt1
           exit
         end if
       end do

!      5th: k - (n-1)dk

       do ikpt1 = 1, nkpt
         diffk(:) = abs(dtset%kptns(:,ikpt1) - dtset%kptns(:,ikpt) - dk(:) + dtset%qptn(:))
         if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
           dtefield%ikpt_dk(ikpt,5,idir) = ikpt1
           exit
         end if
       end do

!      7th: k+n dk

       do ikpt1 = 1, nkpt
         diffk(:) = abs(dtset%kptns(:,ikpt1) - dtset%kptns(:,ikpt) - dtset%qptn(:))
         if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
           dtefield%ikpt_dk(ikpt,7,idir) = ikpt1
           exit
         end if
       end do

!      8th: k-n dk

       do ikpt1 = 1, nkpt
         diffk(:) = abs(dtset%kptns(:,ikpt1) - dtset%kptns(:,ikpt) + dtset%qptn(:))
         if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
           dtefield%ikpt_dk(ikpt,8,idir) = ikpt1
           exit
         end if
       end do

!      9th: -k

       do ikpt1 = 1, nkpt
         diffk(:) = abs(dtset%kptns(:,ikpt1) + dtset%kptns(:,ikpt))
         if(sum(diffk(:))  < 3*tol8) then
           dtefield%ikpt_dk(ikpt,9,idir) = ikpt1
           exit
         end if
       end do

     end do     ! ikpt



!    Find the string length, starting from k point 1
!    (all strings must have the same number of points)

     nkstr = 1
     ikpt1 = 1
     do ikpt = 1, nkpt
       ikpt1 = dtefield%ikpt_dk(ikpt1,1,idir)
       if (ikpt1 == 1) exit
       nkstr = nkstr + 1
     end do

!    Check that the string length is a divisor of nkpt
     if(mod(nkpt,nkstr) /= 0) then
       write(message,'(a,i0,a,i0)')' The string length = ',nkstr,', is not a divisor of nkpt =',nkpt
       MSG_BUG(message)
     end if
     dtefield%nkstr(idir) = nkstr
     dtefield%nstr(idir)  = nkpt/nkstr

   end if      ! dtset%rfdir(idir) == 1

   write(message,'(a,i1,a,i3,a,i3)')&
&   '  dfptff_initberry: for direction ',idir,', nkstr = ',dtefield%nkstr(idir),&
&   ', nstr = ',dtefield%nstr(idir)
   call wrtout(std_out,message,'COLL')

 end do     ! close loop over idir

 dtefield%maxnstr  = maxval(dtefield%nstr(:))
 dtefield%maxnkstr = maxval(dtefield%nkstr(:))
 ABI_ALLOCATE(dtefield%idxkstr,(dtefield%maxnkstr,dtefield%maxnstr,3))
 dtefield%idxkstr(:,:,:) = 0


!Build the different strings------------------------------------------

 ABI_ALLOCATE(kpt_mark,(nkpt))
 do idir = 1, 3

   if (dtset%rfdir(idir) == 1) then

     iunmark = 1
     kpt_mark(:) = 0
     do istr = 1, dtefield%nstr(idir)

       do while(kpt_mark(iunmark) /= 0)
         iunmark = iunmark + 1
       end do
       dtefield%idxkstr(1,istr,idir) = iunmark
       kpt_mark(iunmark) = 1
       do ikstr = 2, dtefield%nkstr(idir)
         ikpt1 = dtefield%idxkstr(ikstr-1,istr,idir)
         ikpt2 = dtefield%ikpt_dk(ikpt1,1,idir)
         dtefield%idxkstr(ikstr,istr,idir) = ikpt2
         kpt_mark(ikpt2) = 1
       end do

     end do    ! istr

   end if         ! rfdir(idir) == 1

 end do           ! close loop over idir

 ABI_DEALLOCATE(kpt_mark)


!Build the array pwindall that is needed to compute the different overlap matrices
!at k +- dk

 ABI_ALLOCATE(kg_tmp,(3,max(mpw,mpw1)*mkmem))
 ABI_ALLOCATE(kpt1,(3,nkpt))
 ABI_ALLOCATE(npwarr_tmp,(nkpt))
 ABI_ALLOCATE(npwtot,(nkpt))
 ecut_eff = dtset%ecut*(dtset%dilatmx)**2

 do idir = 1, 3

   if (dtset%rfdir(idir) == 1) then

     dk(:) = dtefield%dkvecs(:,idir)

     do ifor = 1, 2

       if (ifor == 2) dk(:) = -1_dp*dk(:)

!      Build the array kpt1 = kptns + dk
!      all k-poins of kpt1 must be in the same BZ as those of kptns
       kpt1(:,:) = zero
       do ikpt = 1, nkpt
         ikpt1 = dtefield%ikpt_dk(ikpt,ifor,idir)
         kpt1(:,ikpt) = dtset%kptns(:,ikpt1)
       end do  ! close loop over ikpt

!      Set up the basis sphere of plane waves at kpt1
       kg_tmp(:,:) = 0
       call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg_tmp,&
&       kpt1,dtset%mkmem,dtset%nband,nkpt,'PERS',mpi_enreg,mpw,&
&       npwarr_tmp,npwtot,dtset%nsppol)

       ikg = 0 ; ikg1 = 0
       do ikpt = 1, nkpt

         nband_k = dtset%nband(ikpt)
         ikpt1 = dtefield%ikpt_dk(ikpt,ifor,idir)


         dg(:) = nint(dtset%kptns(:,ikpt) + dk(:) + tol10 - &
&         dtset%kptns(:,ikpt1))

         flag = 0; orig = 1
         if (dg(1)*dg(1) + dg(2)*dg(2) + dg(3)*dg(3) > 0) flag = 1

         ikpt1 = dtefield%ikpt_dk(ikpt,ifor,idir)
         npw_k  = npwarr(ikpt)
         npw_k1 = npwarr_tmp(ikpt)

         do ipw = 1, npw_k
           do jpw = orig, npw_k1
             if ((kg(1,ikg + ipw) == kg_tmp(1,ikg1 + jpw) - dg(1)).and. &
             (kg(2,ikg + ipw) == kg_tmp(2,ikg1 + jpw) - dg(2)).and. &
             (kg(3,ikg + ipw) == kg_tmp(3,ikg1 + jpw) - dg(3))) then
               pwindall((ikpt-1)*max(mpw,mpw1) + ipw,ifor,idir) = jpw
               if (flag == 0) orig = jpw
               exit
             end if
           end do
         end do

         ikg  = ikg + npw_k
         ikg1 = ikg1 + npw_k1

       end do    ! close loop over ikpt

     end do    ! close loop over ifor

   end if      ! rfdir(idir) == 1

 end do        ! close loop over idir

!------------------------------------------------------------------------------
!<u_q|u_q>
!at k +- dk

 do idir = 1, 3

   if (dtset%rfdir(idir) == 1) then

     dk(:) = dtefield%dkvecs(:,idir)

     do ifor = 1, 2

       if (ifor == 2) dk(:) = -1_dp*dk(:)

!      Build the array kpt1 = kptns + qptn + dk
!      all k-poins of kpt1 must be in the same BZ as those of kptns
       kpt1(:,:) = zero
       do ikpt = 1, nkpt
         ikpt1 = dtefield%ikpt_dk(ikpt,ifor,idir)
         kpt1(:,ikpt) = dtset%kptns(:,ikpt1)+dtset%qptn(:)
       end do  ! close loop over ikpt

!      Set UP THE BASIS SPHERE OF PLANE waves at kpt1
       kg_tmp(:,:) = 0
       call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg_tmp,&
&       kpt1,dtset%mkmem,dtset%nband,nkpt,'PERS',mpi_enreg,mpw1,&
&       npwarr_tmp,npwtot,dtset%nsppol)


       ikg = 0 ; ikg1 = 0
       do ikpt = 1, nkpt

         nband_k = dtset%nband(ikpt)
         ikpt1 = dtefield%ikpt_dk(ikpt,ifor,idir)

         if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,1,mpi_enreg%me)) cycle

         dg(:) = nint(dtset%kptns(:,ikpt) + dk(:) + tol10 - &
&         dtset%kptns(:,ikpt1))

         flag = 0; orig = 1
         if (dg(1)*dg(1) + dg(2)*dg(2) + dg(3)*dg(3) > 0) flag = 1

         ikpt1 = dtefield%ikpt_dk(ikpt,ifor,idir)
         npw_k  = npwar1(ikpt)
         npw_k1 = npwarr_tmp(ikpt)

         do ipw = 1, npw_k
           do jpw = orig, npw_k1
             if ((kg1(1,ikg + ipw) == kg_tmp(1,ikg1 + jpw) - dg(1)).and. &
             (kg1(2,ikg + ipw) == kg_tmp(2,ikg1 + jpw) - dg(2)).and. &
             (kg1(3,ikg + ipw) == kg_tmp(3,ikg1 + jpw) - dg(3))) then
               pwindall((ikpt-1)*max(mpw,mpw1) +ipw,ifor+2,idir) = jpw
               if (flag == 0) orig = jpw
               exit
             end if
           end do
         end do

         ikg  = ikg + npw_k
         ikg1 = ikg1 + npw_k1

       end do    ! close loop over ikpt
     end do    ! close loop over ifor
   end if      ! rfdir(idir) == 1
 end do        ! close loop over idir


!---------------------------------------------------------------------------

 do idir = 1, 3

   if (dtset%rfdir(idir) == 1) then

     dk(:) = dtset%qptn(:) + dtefield%dkvecs(:,idir)

     do ifor = 1, 2

       if (ifor == 2) dk(:) = dtset%qptn(:) - dtefield%dkvecs(:,idir)

!      Build the array kpt1 = kptns + qptn + dk
!      all k-poins of kpt1 must be in the same BZ as those of kptns
       kpt1(:,:) = zero
       do ikpt = 1, nkpt
         ikpt1 = dtefield%ikpt_dk(ikpt,ifor+2,idir)
         kpt1(:,ikpt) = dtset%kptns(:,ikpt1)
       end do  ! close loop over ikpt

!      Set UP THE BASIS SPHERE OF PLANE waves at kpt1
       kg_tmp(:,:) = 0
       call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg_tmp,&
&       kpt1,dtset%mkmem,dtset%nband,nkpt,'PERS',mpi_enreg,mpw,&
&       npwarr_tmp,npwtot,dtset%nsppol)

       ikg = 0 ; ikg1 = 0
       do ikpt = 1, nkpt

         nband_k = dtset%nband(ikpt)
         ikpt1 = dtefield%ikpt_dk(ikpt,ifor+2,idir)

         if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,1,mpi_enreg%me)) cycle

         dg(:) = nint(dtset%kptns(:,ikpt) + dk(:) + tol10 - &
&         dtset%kptns(:,ikpt1))

         flag = 0; orig = 1
         if (dg(1)*dg(1) + dg(2)*dg(2) + dg(3)*dg(3) > 0) flag = 1

         ikpt1 = dtefield%ikpt_dk(ikpt,ifor+2,idir)
         npw_k  = npwar1(ikpt)
         npw_k1 = npwarr_tmp(ikpt)

         do ipw = 1, npw_k
           do jpw = orig, npw_k1
             if ((kg1(1,ikg + ipw) == kg_tmp(1,ikg1 + jpw) - dg(1)).and. &
             (kg1(2,ikg + ipw) == kg_tmp(2,ikg1 + jpw) - dg(2)).and. &
             (kg1(3,ikg + ipw) == kg_tmp(3,ikg1 + jpw) - dg(3))) then
               pwindall((ikpt-1)*max(mpw,mpw1)+ipw,ifor+4,idir) = jpw
               if (flag == 0) orig = jpw
               exit
             end if
           end do
         end do
         ikg  = ikg + npw_k
         ikg1 = ikg1 + npw_k1

       end do    ! close loop over ikpt
     end do    ! close loop over ifor
   end if      ! rfdir(idir) == 1
 end do        ! close loop over idir===============================================================


!Build the array pwind3 that is needed to compute the overlap matrices===============================

 do idir = 1, 3

   if (dtset%rfdir(idir) == 1) then

     dk(:) = - dtset%qptn(:) + dtefield%dkvecs(:,idir)

     do ifor = 1, 2

       if (ifor == 2) dk(:) = - dtset%qptn(:) -  dtefield%dkvecs(:,idir)

!      Build the array kpt1 = kptns + qptn + dk
!      all k-poins of kpt1 must be in the same BZ as those of kptns
       kpt1(:,:) = zero
       do ikpt = 1, nkpt
         ikpt1 = dtefield%ikpt_dk(ikpt,ifor+4,idir)
         kpt1(:,ikpt) = dtset%kptns(:,ikpt1)+ dtset%qptn(:)
       end do  ! close loop over ikpt

!      Set UP THE BASIS SPHERE OF PLANE waves at kpt1
       kg_tmp(:,:) = 0
       call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg_tmp,&
&       kpt1,dtset%mkmem,dtset%nband,nkpt,'PERS',mpi_enreg,mpw1,&
&       npwarr_tmp,npwtot,dtset%nsppol)

       ikg = 0 ; ikg1 = 0
       do ikpt = 1, nkpt

         nband_k = dtset%nband(ikpt)
         ikpt1 = dtefield%ikpt_dk(ikpt,ifor+4,idir)

         if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,1,mpi_enreg%me)) cycle

         dg(:) = nint(dtset%kptns(:,ikpt) + dk(:) + tol10 - &
&         dtset%kptns(:,ikpt1))

         flag = 0; orig = 1
         if (dg(1)*dg(1) + dg(2)*dg(2) + dg(3)*dg(3) > 0) flag = 1

         ikpt1 = dtefield%ikpt_dk(ikpt,ifor+4,idir)
         npw_k  = npwarr(ikpt)
         npw_k1 = npwarr_tmp(ikpt)

         do ipw = 1, npw_k
           do jpw = orig, npw_k1
             if ((kg(1,ikg + ipw) == kg_tmp(1,ikg1 + jpw) - dg(1)).and. &
             (kg(2,ikg + ipw) == kg_tmp(2,ikg1 + jpw) - dg(2)).and. &
             (kg(3,ikg + ipw) == kg_tmp(3,ikg1 + jpw) - dg(3))) then
               pwindall((ikpt-1)*max(mpw,mpw1) + ipw,ifor+6,idir) = jpw
               if (flag == 0) orig = jpw
               exit
             end if
           end do
         end do
         ikg  = ikg + npw_k
         ikg1 = ikg1 + npw_k1

       end do    ! close loop over ikpt
     end do    ! close loop over ifor
   end if      ! rfdir(idir) == 1
 end do        ! close loop over idir====================================================================

 ABI_DEALLOCATE(kg_tmp)
 ABI_DEALLOCATE(kpt1)
 ABI_DEALLOCATE(npwarr_tmp)
 ABI_DEALLOCATE(npwtot)

end subroutine dfptff_initberry
!!***

!!****f* ABINIT/dfptff_gradberry
!! NAME
!! dfptff_gradberry
!!
!! FUNCTION
!! Calculation of the gradient of Berry-phase term in finite electric field.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2020 ABINIT group (XW).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cg(2,mpw*nspinor*mband*mkmem*nsppol) = planewave coefficients of wavefunctions
!! cg1(2,mpw1*nspinor*mband*mk1mem*nsppol) = pw coefficients of
!! RF wavefunctions at k,q.
!! dtefield = variables related to finite electric field calculation
!! ikpt = the index of the current k point
!! isppol=1 for unpolarized, 2 for spin-polarized
!! mband =  maximum number of bands
!! mkmem = maximum number of k-points in core memory
!! mpw = maximum number of plane waves
!! mpw1 = maximum number of plane waves for response wavefunctions
!! nkpt = number of k points
!! npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!! npwar1(nkpt) = number of planewaves in basis and boundary for response wfs
!! nspinor = 1 for scalar wfs, 2 for spinor wfs
!! nsppol = 1 for unpolarized, 2 for spin-polarized
!! qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3) =
!! inverse of the overlap matrix
!! pwindall(max(mpw,mpw1)*mkmem,8,3) = array used to compute the overlap matrices
!! pwindall(:,1,:) <- <u^(0)_i|u^(0)_i+1>
!! pwindall(:,2,:) <- <u^(0)_i|u^(0)_i-1>
!! pwindall(:,3,:) <- <u^(1)_i|u^(1)_i+1>
!! pwindall(:,4,:) <- <u^(1)_i|u^(1)_i-1>
!! pwindall(:,5,:) <- <u^(1)_i|u^(0)_i+n+1>
!! pwindall(:,6,:) <- <u^(1)_i|u^(0)_i+n-1>
!! pwindall(:,7,:) <- <u^(0)_i|u^(1)_i-n+1>
!! pwindall(:,8,:) <- <u^(0)_i|u^(1)_i-n-1>
!!
!! OUTPUT
!! grad_berry(2,mpw1,dtefield%mband_occ) = the gradient of the Berry phase term
!!
!! PARENTS
!!      dfpt_vtorho
!!
!! CHILDREN
!!      overlap_g
!!
!! SOURCE

subroutine dfptff_gradberry(cg,cg1,dtefield,grad_berry,ikpt,isppol,mband,mpw,mpw1,mkmem,mk1mem,nkpt,&
&                     npwarr,npwar1,nspinor,nsppol,qmat,pwindall)

!Arguments ----------------------------------------
!scalars
 integer,intent(in) :: ikpt,isppol,mband,mk1mem,mkmem,mpw,mpw1,nkpt,nspinor
 integer,intent(in) :: nsppol
 type(efield_type),intent(in) :: dtefield
!arrays
 integer,intent(in) :: npwar1(nkpt),npwarr(nkpt)
 integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3)
 real(dp),intent(out) :: grad_berry(2,mpw1,dtefield%mband_occ)

!Local variables -------------------------
!scalars
 integer :: iband,icg,icg1,idir,ikpt1
 integer :: ikpt1m,ikptn,ikptnm,ikptnp1,ipw,jband,jpw,kband
 integer :: mpw_tmp,npw_k1,npw_k2,pwmax,pwmin
 real(dp) :: doti,dotr,fac,wfi,wfr
!arrays
 integer,allocatable :: pwind_tmp(:)
 real(dp) :: z1(2),z2(2)
 real(dp),allocatable :: Amat(:,:,:),Bmat(:,:,:),s1mat(:,:,:),vect1(:,:)
 real(dp),allocatable :: vect2(:,:)

! *************************************************************************

 mpw_tmp=max(mpw,mpw1)
 ABI_ALLOCATE(vect1,(2,0:mpw_tmp))
 ABI_ALLOCATE(vect2,(2,0:mpw_tmp))
 ABI_ALLOCATE(s1mat,(2,dtefield%mband_occ,dtefield%mband_occ))
 ABI_ALLOCATE(pwind_tmp,(mpw_tmp))
 ABI_ALLOCATE(Amat,(2,dtefield%mband_occ,dtefield%mband_occ))
 ABI_ALLOCATE(Bmat,(2,dtefield%mband_occ,dtefield%mband_occ))
 vect1(:,0) = zero ; vect2(:,0) = zero
 s1mat(:,:,:)=zero
 grad_berry(:,:,:) = zero

 do idir=1,3
   fac = dtefield%efield_dot(idir)*dble(nkpt)/&
&   (dble(dtefield%nstr(idir))*four_pi)

!  prepare
   ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
   icg1 = dtefield%cgindex(ikpt1,isppol+nsppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwar1(ikpt1)
   pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,3,idir)

   do ipw = 1, npw_k1
     jpw = pwind_tmp(ipw)

     if (jpw > 0) then
       do iband = 1, dtefield%mband_occ
         wfr = cg1(1,icg1 + (iband - 1)*npw_k2*nspinor + jpw)
         wfi = cg1(2,icg1 + (iband - 1)*npw_k2*nspinor + jpw)

         do  jband = 1, dtefield%mband_occ

           grad_berry(1,ipw,jband) = &
&           grad_berry(1,ipw,jband) + &
&           fac*qmat(1,iband,jband,ikpt,1,idir)*wfr - fac*qmat(2,iband,jband,ikpt,1,idir)*wfi

           grad_berry(2,ipw,jband) = &
&           grad_berry(2,ipw,jband) + &
&           fac*qmat(1,iband,jband,ikpt,1,idir)*wfi + fac*qmat(2,iband,jband,ikpt,1,idir)*wfr

         end do
       end do
     end if

   end do

!  compute <u^(0)_{k_j+n}|u^(1)_{k_j+1,q}> matrix----------------------------------------------------

!  prepare to calculate overlap matrix
   ikptn = dtefield%ikpt_dk(ikpt,7,idir)
   ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
   icg = dtefield%cgindex(ikptn,isppol)
   icg1 = dtefield%cgindex(ikpt1,isppol+nsppol)
   npw_k1 = npwarr(ikptn)
   npw_k2 = npwar1(ikpt1)
   pwind_tmp(1:npw_k1) = pwindall((ikptn-1)*mpw_tmp+1:(ikptn-1)*mpw_tmp+npw_k1,7,idir)


   vect1(:,0) = zero ; vect2(:,0) = zero
   do jband = 1, dtefield%mband_occ
     vect2(:,1:npw_k2) = &
&     cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
     if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero

     do iband = 1, dtefield%mband_occ

       pwmin = (iband-1)*npw_k1*nspinor
       pwmax = pwmin + npw_k1*nspinor
       vect1(:,1:npw_k1) = &
&       cg(:,icg + 1 + pwmin:icg + pwmax)
       if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
       call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&       vect1,vect2)
       s1mat(1,iband,jband) = dotr
       s1mat(2,iband,jband) = doti

     end do    ! iband
   end do    !jband

!  compute <u^(0)_{-k_j+1}|u^(1)_{-k_j+n},q> matrix--------------------

!  prepare to calculate overlap matrix
   ikptn = dtefield%ikpt_dk(ikpt,7,idir)
   ikptnm= dtefield%ikpt_dk(ikptn,9,idir)
   ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
   ikpt1m= dtefield%ikpt_dk(ikpt1,9,idir)

   icg = dtefield%cgindex(ikpt1m,isppol)
   icg1 = dtefield%cgindex(ikptnm,isppol+nsppol)
   npw_k1 = npwarr(ikpt1m)
   npw_k2 = npwar1(ikptnm)
   pwind_tmp(1:npw_k1) = pwindall((ikpt1m-1)*mpw_tmp+1:(ikpt1m-1)*mpw_tmp+npw_k1,7,idir)

   vect1(:,0) = zero ; vect2(:,0) = zero
   do jband = 1, dtefield%mband_occ
     vect2(:,1:npw_k2) = &
&     cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
     if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero

     do iband = 1, dtefield%mband_occ

       pwmin = (iband-1)*npw_k1*nspinor
       pwmax = pwmin + npw_k1*nspinor
       vect1(:,1:npw_k1) = &
&       cg(:,icg + 1 + pwmin:icg + pwmax)
       if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
       call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&       vect1,vect2)

       s1mat(1,jband,iband) = s1mat(1,jband,iband) + dotr
       s1mat(2,jband,iband) = s1mat(2,jband,iband) + doti

     end do    ! iband
   end do    !jband

   Amat(:,:,:)=zero

!  calculate Amat
   do iband=1, dtefield%mband_occ
     do jband=1, dtefield%mband_occ
       do kband=1, dtefield%mband_occ
         Amat(1,iband,jband) = Amat(1,iband,jband) + s1mat(1,iband,kband)*&
&         qmat(1,kband,jband,ikpt,1,idir)&
&         - s1mat(2,iband,kband)*qmat(2,kband,jband,ikpt,1,idir)
         Amat(2,iband,jband) = Amat(2,iband,jband) + s1mat(1,iband,kband)*&
&         qmat(2,kband,jband,ikpt,1,idir)&
&         + s1mat(2,iband,kband)*qmat(1,kband,jband,ikpt,1,idir)
       end do
     end do
   end do

   Bmat(:,:,:)=zero

!  calculate Bmat
   ikptn = dtefield%ikpt_dk(ikpt,7,idir)
   do iband=1, dtefield%mband_occ
     do jband=1, dtefield%mband_occ
       do kband=1, dtefield%mband_occ
         Bmat(1,jband,kband) = Bmat(1,jband,kband) + Amat(1,iband,kband)*&
&         qmat(1,jband,iband,ikptn,1,idir)&
&         - Amat(2,iband,kband)*qmat(2,jband,iband,ikptn,1,idir)
         Bmat(2,jband,kband) = Bmat(2,jband,kband) + Amat(1,iband,kband)*&
&         qmat(2,jband,iband,ikptn,1,idir)&
&         + Amat(2,iband,kband)*qmat(1,jband,iband,ikptn,1,idir)
       end do
     end do
   end do

!  calc. the second term of gradient------------------------------

!  preparation

   ikptnp1 = dtefield%ikpt_dk(ikpt,3,idir)
   icg = dtefield%cgindex(ikptnp1,isppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwarr(ikptnp1)
   pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,5,idir)

   z1(:) = zero
   z2(:) = zero

   do ipw = 1, npw_k1

     jpw = pwind_tmp(ipw)

     if (jpw > 0) then

       do iband = 1, dtefield%mband_occ
         wfr = cg(1,icg + (iband - 1)*npw_k2*nspinor + jpw)
         wfi = cg(2,icg + (iband - 1)*npw_k2*nspinor + jpw)

         do jband=1, dtefield%mband_occ

           grad_berry(1,ipw,jband) = grad_berry(1,ipw,jband) &
&           - fac*(Bmat(1,iband,jband)*wfr - Bmat(2,iband,jband)*wfi)
           grad_berry(2,ipw,jband) = grad_berry(2,ipw,jband) &
&           - fac*(Bmat(1,iband,jband)*wfi + Bmat(2,iband,jband)*wfr)

         end do
       end do
     end if
   end do

!  Second part of gradient of Berry phase++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   vect1(:,0) = zero ; vect2(:,0) = zero

!  prepare
   ikpt1 = dtefield%ikpt_dk(ikpt,2,idir)
   icg1 = dtefield%cgindex(ikpt1,isppol+nsppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwar1(ikpt1)
   pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,4,idir)

!  write(std_out,*)'dfpt_cgwf:pwind_tmp',pwind_tmp
!  stop

   do ipw = 1, npw_k1

     jpw = pwind_tmp(ipw)

     if (jpw > 0) then
       do iband = 1, dtefield%mband_occ
         wfr = cg1(1,icg1 + (iband - 1)*npw_k2*nspinor + jpw)
         wfi = cg1(2,icg1 + (iband - 1)*npw_k2*nspinor + jpw)

         do  jband = 1, dtefield%mband_occ

           grad_berry(1,ipw,jband) = &
&           grad_berry(1,ipw,jband) - &
&           fac*qmat(1,iband,jband,ikpt,2,idir)*wfr + fac*qmat(2,iband,jband,ikpt,2,idir)*wfi

           grad_berry(2,ipw,jband) = &
&           grad_berry(2,ipw,jband) - &
&           fac*qmat(1,iband,jband,ikpt,2,idir)*wfi - fac*qmat(2,iband,jband,ikpt,2,idir)*wfr

         end do
       end do
     end if
   end do


!  compute <u^(0)_{k_j+n}|u^(1)_{k_j-1,q}> matrix----------------------------------------------------

!  prepare to calculate overlap matrix
   ikptn = dtefield%ikpt_dk(ikpt,7,idir)
   ikpt1 = dtefield%ikpt_dk(ikpt,2,idir)
   icg = dtefield%cgindex(ikptn,isppol)
   icg1 = dtefield%cgindex(ikpt1,isppol+nsppol)
   npw_k1 = npwarr(ikptn)
   npw_k2 = npwar1(ikpt1)
   pwind_tmp(1:npw_k1) =pwindall((ikptn-1)*mpw_tmp+1:(ikptn-1)*mpw_tmp+npw_k1,8,idir)

   vect1(:,0) = zero ; vect2(:,0) = zero
   do jband = 1, dtefield%mband_occ
     vect2(:,1:npw_k2) = &
&     cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
     if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero

     do iband = 1, dtefield%mband_occ

       pwmin = (iband-1)*npw_k1*nspinor
       pwmax = pwmin + npw_k1*nspinor
       vect1(:,1:npw_k1) = &
&       cg(:,icg + 1 + pwmin:icg + pwmax)
       if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
       call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&       vect1,vect2)
       s1mat(1,iband,jband) = dotr
       s1mat(2,iband,jband) = doti

     end do    ! iband
   end do    !jband

!  compute <u^(0)_{-k_j-1}|u^(1)_{-k_j+n,q}> matrix-----------------------------------------------------

!  prepare to calculate overlap matrix
   ikptn = dtefield%ikpt_dk(ikpt,7,idir)
   ikptnm= dtefield%ikpt_dk(ikptn,9,idir)
   ikpt1 = dtefield%ikpt_dk(ikpt,2,idir)
   ikpt1m= dtefield%ikpt_dk(ikpt1,9,idir)
   icg = dtefield%cgindex(ikpt1m,isppol)
   icg1 = dtefield%cgindex(ikptnm,isppol+nsppol)
   npw_k1 = npwarr(ikpt1m)
   npw_k2 = npwar1(ikptnm)
   pwind_tmp(1:npw_k1) =pwindall((ikpt1m-1)*mpw_tmp+1:(ikpt1m-1)*mpw_tmp+npw_k1,8,idir)



   vect1(:,0) = zero ; vect2(:,0) = zero
   do jband = 1, dtefield%mband_occ
     vect2(:,1:npw_k2) = &
&     cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
     if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
     do iband = 1, dtefield%mband_occ
       pwmin = (iband-1)*npw_k1*nspinor
       pwmax = pwmin + npw_k1*nspinor
       vect1(:,1:npw_k1) = &
&       cg(:,icg + 1 + pwmin:icg + pwmax)
       if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
       call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&       vect1,vect2)
       s1mat(1,jband,iband) = s1mat(1,jband,iband) + dotr
       s1mat(2,jband,iband) = s1mat(2,jband,iband) + doti

     end do    ! iband
   end do    !jband

   Amat(:,:,:)=zero

!  calculate Amat
   do iband=1, dtefield%mband_occ
     do jband=1, dtefield%mband_occ
       do kband=1, dtefield%mband_occ
         Amat(1,iband,jband) = Amat(1,iband,jband) + s1mat(1,iband,kband)*&
&         qmat(1,kband,jband,ikpt,2,idir)&
&         - s1mat(2,iband,kband)*qmat(2,kband,jband,ikpt,2,idir)
         Amat(2,iband,jband) = Amat(2,iband,jband) + s1mat(1,iband,kband)*&
&         qmat(2,kband,jband,ikpt,2,idir)&
&         + s1mat(2,iband,kband)*qmat(1,kband,jband,ikpt,2,idir)
       end do
     end do
   end do

   Bmat(:,:,:)=zero

!  calculate Bmat
   ikptn = dtefield%ikpt_dk(ikpt,7,idir)
   do iband=1, dtefield%mband_occ
     do jband=1, dtefield%mband_occ
       do kband=1, dtefield%mband_occ
         Bmat(1,jband,kband) = Bmat(1,jband,kband) + Amat(1,iband,kband)*&
&         qmat(1,jband,iband,ikptn,2,idir)&
&         - Amat(2,iband,kband)*qmat(2,jband,iband,ikptn,2,idir)
         Bmat(2,jband,kband) = Bmat(2,jband,kband) + Amat(1,iband,kband)*&
&         qmat(2,jband,iband,ikptn,2,idir)&
         + Amat(2,iband,kband)*qmat(1,jband,iband,ikptn,2,idir)
       end do
     end do
   end do

!  calc. the second term of gradient------------------------------

!  preparation

   ikptnp1 = dtefield%ikpt_dk(ikpt,4,idir)
   icg = dtefield%cgindex(ikptnp1,isppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwarr(ikptnp1)
   pwind_tmp(1:npw_k1) =pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,6,idir)

   z1(:) = zero
   z2(:) = zero
   do ipw = 1, npw_k1

     jpw = pwind_tmp(ipw)

     if (jpw > 0) then

       do iband = 1, dtefield%mband_occ
         wfr = cg(1,icg + (iband - 1)*npw_k2*nspinor + jpw)
         wfi = cg(2,icg + (iband - 1)*npw_k2*nspinor + jpw)

         do jband=1, dtefield%mband_occ

           grad_berry(1,ipw,jband) = grad_berry(1,ipw,jband) + fac*(Bmat(1,iband,jband)*wfr - Bmat(2,iband,jband)*wfi)
           grad_berry(2,ipw,jband) = grad_berry(2,ipw,jband) + fac*(Bmat(1,iband,jband)*wfi + Bmat(2,iband,jband)*wfr)

         end do
       end do
     end if
   end do

 end do !idir

 ABI_DEALLOCATE(vect1)
 ABI_DEALLOCATE(vect2)
 ABI_DEALLOCATE(s1mat)
 ABI_DEALLOCATE(Amat)
 ABI_DEALLOCATE(Bmat)
 ABI_DEALLOCATE(pwind_tmp)

end subroutine dfptff_gradberry
!!***

!!****f* ABINIT/dfptff_gbefd
!! NAME
!! dfptff_gbefd
!!
!! FUNCTION
!! calculate the gradient of the second order \Omega E \cdot P
!! term, Eq.(23) in PRB 75, 115116(2007) [[cite:Wang2007]].
!!
!! INPUTS
!! cg(2,mpw*nspinor*mband*mkmem*nsppol) = planewave coefficients of wavefunctions
!! cg1(2,mpw1*nspinor*mband*mk1mem*nsppol) = pw coefficients of
!! RF wavefunctions at k,q.
!! dtefield = variables related to response Berry-phase calculation
!! ikpt = the index of the current k point
!! isppol = the index of the spin component
!! mband =  maximum number of bands
!! mkmem = maximum number of k-points in core memory
!! mpw = maximum number of plane waves
!! mpw1 = maximum number of plane waves for response wavefunctions
!! nkpt = number of k points
!! npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!! npwar1(nkpt) = number of planewaves in basis and boundary for response wfs
!! nspinor = 1 for scalar wfs, 2 for spinor wfs
!! nsppol = 1 for unpolarized, 2 for spin-polarized
!! qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3) =
!! inverse of the overlap matrix
!! pwindall(max(mpw,mpw1)*mkmem,8,3) = array used to compute the overlap matrices
!! pwindall(:,1,:) <- <u^(0)_i|u^(0)_i+1>
!! pwindall(:,2,:) <- <u^(0)_i|u^(0)_i-1>
!! pwindall(:,3,:) <- <u^(1)_i|u^(1)_i+1>
!! pwindall(:,4,:) <- <u^(1)_i|u^(1)_i-1>
!! pwindall(:,5,:) <- <u^(1)_i|u^(0)_i+n+1>
!! pwindall(:,6,:) <- <u^(1)_i|u^(0)_i+n-1>
!! pwindall(:,7,:) <- <u^(0)_i|u^(1)_i-n+1>
!! pwindall(:,8,:) <- <u^(0)_i|u^(1)_i-n-1>
!!
!! OUTPUT
!! grad_berry = the gradient of the Berry phase term
!!
!! PARENTS
!!      dfpt_vtorho
!!
!! CHILDREN
!!      overlap_g
!!
!! SOURCE

subroutine dfptff_gbefd(cg,cg1,dtefield,grad_berry,idir_efield,ikpt,isppol,mband,mpw,mpw1,mkmem,mk1mem,nkpt,&
&                 npwarr,npwar1,nspinor,nsppol,qmat,pwindall,rprimd)

!Arguments ----------------------------------------
!scalars
 integer,intent(in) :: idir_efield,ikpt,isppol,mband,mk1mem,mkmem,mpw,mpw1,nkpt
 integer,intent(in) :: nspinor,nsppol
 type(efield_type),intent(in) :: dtefield
!arrays
 integer,intent(in) :: npwar1(nkpt),npwarr(nkpt)
 integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(out) :: grad_berry(2,mpw1,dtefield%mband_occ)

!Local variables -------------------------
!scalars
 integer :: iband,icg,icg1,idir,ikpt1
 integer :: ikptn,ikptnp1,ipw,jband,jpw,kband
 integer :: mpw_tmp,npw_k1,npw_k2,pwmax,pwmin
 real(dp) :: doti,dotr,fac,wfi,wfr
!arrays
 integer,allocatable :: pwind_tmp(:)
 real(dp) :: z1(2),z2(2)
 real(dp),allocatable :: Amat(:,:,:),Bmat(:,:,:),s1mat(:,:,:),vect1(:,:)
 real(dp),allocatable :: vect2(:,:)

! *************************************************************************

 mpw_tmp=max(mpw,mpw1)
 ABI_ALLOCATE(vect1,(2,0:mpw_tmp))
 ABI_ALLOCATE(vect2,(2,0:mpw_tmp))
 ABI_ALLOCATE(s1mat,(2,dtefield%mband_occ,dtefield%mband_occ))
 ABI_ALLOCATE(pwind_tmp,(mpw_tmp))
 ABI_ALLOCATE(Amat,(2,dtefield%mband_occ,dtefield%mband_occ))
 ABI_ALLOCATE(Bmat,(2,dtefield%mband_occ,dtefield%mband_occ))
 vect1(:,0) = zero ; vect2(:,0) = zero
 s1mat(:,:,:)=zero
 grad_berry(:,:,:) = zero

 do idir=1,3
   fac = dtefield%efield_dot(idir)*dble(nkpt)/&
&   (dble(dtefield%nstr(idir))*four_pi)

!  prepare
   ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
   icg1 = dtefield%cgindex(ikpt1,isppol+nsppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwar1(ikpt1)
   pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,3,idir)

   do ipw = 1, npw_k1
     jpw = pwind_tmp(ipw)
     if (jpw > 0) then
       do iband = 1, dtefield%mband_occ
         wfr = cg1(1,icg1 + (iband - 1)*npw_k2*nspinor + jpw)
         wfi = cg1(2,icg1 + (iband - 1)*npw_k2*nspinor + jpw)
         do  jband = 1, dtefield%mband_occ
           grad_berry(1,ipw,jband) = &
&           grad_berry(1,ipw,jband) + &
&           fac*qmat(1,iband,jband,ikpt,1,idir)*wfr - fac*qmat(2,iband,jband,ikpt,1,idir)*wfi

           grad_berry(2,ipw,jband) = &
&           grad_berry(2,ipw,jband) + &
&           fac*qmat(1,iband,jband,ikpt,1,idir)*wfi + fac*qmat(2,iband,jband,ikpt,1,idir)*wfr
         end do
       end do
     end if
   end do

!  compute <u^(0)_{k_j}|u^(1)_{k_j+1}> matrix----------------------------------------------------

!  prepare to calculate overlap matrix
   ikptn = ikpt
   ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
   icg = dtefield%cgindex(ikptn,isppol)
   icg1 = dtefield%cgindex(ikpt1,isppol+nsppol)
   npw_k1 = npwarr(ikptn)
   npw_k2 = npwar1(ikpt1)
   pwind_tmp(1:npw_k1) = pwindall((ikptn-1)*mpw_tmp+1:(ikptn-1)*mpw_tmp+npw_k1,7,idir)

   vect1(:,0) = zero ; vect2(:,0) = zero
   do jband = 1, dtefield%mband_occ
     vect2(:,1:npw_k2) = &
&     cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
     if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
     do iband = 1, dtefield%mband_occ
       pwmin = (iband-1)*npw_k1*nspinor
       pwmax = pwmin + npw_k1*nspinor
       vect1(:,1:npw_k1) = &
&       cg(:,icg + 1 + pwmin:icg + pwmax)
       if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
       call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&       vect1,vect2)
       s1mat(1,iband,jband) = dotr
       s1mat(2,iband,jband) = doti
     end do    ! iband
   end do    !jband

!  compute <u^(1)_{k_j}|u^(0)_{k_j+1}> matrix-----------------------------------------------------

!  prepare to calculate overlap matrix
   ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
   icg = dtefield%cgindex(ikpt,isppol+nsppol)
   icg1 = dtefield%cgindex(ikpt1,isppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwarr(ikpt1)
   pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,5,idir)

   vect1(:,0) = zero ; vect2(:,0) = zero
   do jband = 1, dtefield%mband_occ
     vect2(:,1:npw_k2) = &
&     cg(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
     if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
     do iband = 1, dtefield%mband_occ
       pwmin = (iband-1)*npw_k1*nspinor
       pwmax = pwmin + npw_k1*nspinor
       vect1(:,1:npw_k1) = &
&       cg1(:,icg + 1 + pwmin:icg + pwmax)
       if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
       call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&       vect1,vect2)
       s1mat(1,jband,iband) = s1mat(1,jband,iband) + dotr
       s1mat(2,jband,iband) = s1mat(2,jband,iband) + doti
     end do    ! iband
   end do    !jband

   Amat(:,:,:)=zero

!  calculate Amat
   do iband=1, dtefield%mband_occ
     do jband=1, dtefield%mband_occ
       do kband=1, dtefield%mband_occ
         Amat(1,iband,jband) = Amat(1,iband,jband) + s1mat(1,iband,kband)*qmat(1,kband,jband,ikpt,1,idir)&
&         - s1mat(2,iband,kband)*qmat(2,kband,jband,ikpt,1,idir)
         Amat(2,iband,jband) = Amat(2,iband,jband) + s1mat(1,iband,kband)*qmat(2,kband,jband,ikpt,1,idir)&
&         + s1mat(2,iband,kband)*qmat(1,kband,jband,ikpt,1,idir)
       end do
     end do
   end do

   Bmat(:,:,:)=zero

!  calculate Bmat
   ikptn = dtefield%ikpt_dk(ikpt,7,idir)
   do iband=1, dtefield%mband_occ
     do jband=1, dtefield%mband_occ
       do kband=1, dtefield%mband_occ
         Bmat(1,jband,kband) = Bmat(1,jband,kband) + Amat(1,iband,kband)*qmat(1,jband,iband,ikptn,1,idir)&
&         - Amat(2,iband,kband)*qmat(2,jband,iband,ikptn,1,idir)
         Bmat(2,jband,kband) = Bmat(2,jband,kband) + Amat(1,iband,kband)*qmat(2,jband,iband,ikptn,1,idir)&
&         + Amat(2,iband,kband)*qmat(1,jband,iband,ikptn,1,idir)
       end do
     end do
   end do

!  calc. the second term of gradient------------------------------

!  preparation

   ikptnp1 = dtefield%ikpt_dk(ikpt,3,idir)
   icg = dtefield%cgindex(ikptnp1,isppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwarr(ikptnp1)
   pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,5,idir)

   z1(:) = zero
   z2(:) = zero

   do ipw = 1, npw_k1
     jpw = pwind_tmp(ipw)
     if (jpw > 0) then
       do iband = 1, dtefield%mband_occ
         wfr = cg(1,icg + (iband - 1)*npw_k2*nspinor + jpw)
         wfi = cg(2,icg + (iband - 1)*npw_k2*nspinor + jpw)
         do jband=1, dtefield%mband_occ
           grad_berry(1,ipw,jband) = grad_berry(1,ipw,jband) - fac*(Bmat(1,iband,jband)*wfr - Bmat(2,iband,jband)*wfi)
           grad_berry(2,ipw,jband) = grad_berry(2,ipw,jband) - fac*(Bmat(1,iband,jband)*wfi + Bmat(2,iband,jband)*wfr)
         end do
       end do
     end if
   end do

!  Second part of gradient of Berry phase++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   vect1(:,0) = zero ; vect2(:,0) = zero

!  prepare
   ikpt1 = dtefield%ikpt_dk(ikpt,2,idir)
   icg1 = dtefield%cgindex(ikpt1,isppol+nsppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwar1(ikpt1)
   pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,4,idir)

   do ipw = 1, npw_k1
     jpw = pwind_tmp(ipw)
     if (jpw > 0) then
       do iband = 1, dtefield%mband_occ
         wfr = cg1(1,icg1 + (iband - 1)*npw_k2*nspinor + jpw)
         wfi = cg1(2,icg1 + (iband - 1)*npw_k2*nspinor + jpw)
         do  jband = 1, dtefield%mband_occ
           grad_berry(1,ipw,jband) = &
&           grad_berry(1,ipw,jband) - &
&           fac*qmat(1,iband,jband,ikpt,2,idir)*wfr + fac*qmat(2,iband,jband,ikpt,2,idir)*wfi

           grad_berry(2,ipw,jband) = &
&           grad_berry(2,ipw,jband) - &
&           fac*qmat(1,iband,jband,ikpt,2,idir)*wfi - fac*qmat(2,iband,jband,ikpt,2,idir)*wfr
         end do
       end do
     end if
   end do

!  compute <u^(0)_{k_j}|u^(1)_{k_j-1}> matrix----------------------------------------------------

!  prepare to calculate overlap matrix
   ikptn = ikpt
   ikpt1 = dtefield%ikpt_dk(ikpt,2,idir)
   icg = dtefield%cgindex(ikptn,isppol)
   icg1 = dtefield%cgindex(ikpt1,isppol+nsppol)
   npw_k1 = npwarr(ikptn)
   npw_k2 = npwar1(ikpt1)
   pwind_tmp(1:npw_k1) =pwindall((ikptn-1)*mpw_tmp+1:(ikptn-1)*mpw_tmp+npw_k1,8,idir)

   vect1(:,0) = zero ; vect2(:,0) = zero
   do jband = 1, dtefield%mband_occ
     vect2(:,1:npw_k2) = &
&     cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
     if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
     do iband = 1, dtefield%mband_occ
       pwmin = (iband-1)*npw_k1*nspinor
       pwmax = pwmin + npw_k1*nspinor
       vect1(:,1:npw_k1) = &
&       cg(:,icg + 1 + pwmin:icg + pwmax)
       if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
       call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&       vect1,vect2)
       s1mat(1,iband,jband) = dotr
       s1mat(2,iband,jband) = doti
     end do    ! iband
   end do    !jband

!  compute <u^(1)_{k_j}|u^(0)_{k_j-1}> matrix-----------------------------------------------------

!  prepare to calculate overlap matrix
   ikpt1 = dtefield%ikpt_dk(ikpt,2,idir)
   icg = dtefield%cgindex(ikpt,isppol+nsppol)
   icg1 = dtefield%cgindex(ikpt1,isppol)
   npw_k1 = npwarr(ikpt)
   npw_k2 = npwar1(ikpt1)
   pwind_tmp(1:npw_k1) =pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,6,idir)

   vect1(:,0) = zero ; vect2(:,0) = zero
   do jband = 1, dtefield%mband_occ
     vect2(:,1:npw_k2) = &
&     cg(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
     if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
     do iband = 1, dtefield%mband_occ
       pwmin = (iband-1)*npw_k1*nspinor
       pwmax = pwmin + npw_k1*nspinor
       vect1(:,1:npw_k1) = &
&       cg1(:,icg + 1 + pwmin:icg + pwmax)
       if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
       call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&       vect1,vect2)
       s1mat(1,jband,iband) = s1mat(1,jband,iband) + dotr
       s1mat(2,jband,iband) = s1mat(2,jband,iband) + doti
     end do    ! iband
   end do    !jband

   Amat(:,:,:)=zero

!  calculate Amat
   do iband=1, dtefield%mband_occ
     do jband=1, dtefield%mband_occ
       do kband=1, dtefield%mband_occ
         Amat(1,iband,jband) = Amat(1,iband,jband) + s1mat(1,iband,kband)*qmat(1,kband,jband,ikpt,2,idir)&
&         - s1mat(2,iband,kband)*qmat(2,kband,jband,ikpt,2,idir)
         Amat(2,iband,jband) = Amat(2,iband,jband) + s1mat(1,iband,kband)*qmat(2,kband,jband,ikpt,2,idir)&
&         + s1mat(2,iband,kband)*qmat(1,kband,jband,ikpt,2,idir)
       end do
     end do
   end do

   Bmat(:,:,:)=zero

!  calculate Bmat
   ikptn = dtefield%ikpt_dk(ikpt,7,idir)
   do iband=1, dtefield%mband_occ
     do jband=1, dtefield%mband_occ
       do kband=1, dtefield%mband_occ
         Bmat(1,jband,kband) = Bmat(1,jband,kband) + Amat(1,iband,kband)*qmat(1,jband,iband,ikptn,2,idir)&
&         - Amat(2,iband,kband)*qmat(2,jband,iband,ikptn,2,idir)
         Bmat(2,jband,kband) = Bmat(2,jband,kband) + Amat(1,iband,kband)*qmat(2,jband,iband,ikptn,2,idir)&
         + Amat(2,iband,kband)*qmat(1,jband,iband,ikptn,2,idir)
       end do
     end do
   end do

!  calc. the second term of gradient------------------------------

!  preparation

   ikptnp1 = dtefield%ikpt_dk(ikpt,4,idir)
   icg = dtefield%cgindex(ikptnp1,isppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwarr(ikptnp1)
   pwind_tmp(1:npw_k1) =pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,6,idir)
   z1(:) = zero
   z2(:) = zero
   do ipw = 1, npw_k1
     jpw = pwind_tmp(ipw)
     if (jpw > 0) then
       do iband = 1, dtefield%mband_occ
         wfr = cg(1,icg + (iband - 1)*npw_k2*nspinor + jpw)
         wfi = cg(2,icg + (iband - 1)*npw_k2*nspinor + jpw)
         do jband=1, dtefield%mband_occ
           grad_berry(1,ipw,jband) = grad_berry(1,ipw,jband) + fac*(Bmat(1,iband,jband)*wfr - Bmat(2,iband,jband)*wfi)
           grad_berry(2,ipw,jband) = grad_berry(2,ipw,jband) + fac*(Bmat(1,iband,jband)*wfi + Bmat(2,iband,jband)*wfr)
         end do
       end do
     end if
   end do

 end do !idir

!!----------------------------------------third part of gradient------------------------------------------------------
 do idir=1,3
   fac = rprimd(idir_efield,idir)*dble(nkpt)/&
&   (dble(dtefield%nstr(idir))*four_pi)

!  prepare
   ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
   icg1 = dtefield%cgindex(ikpt1,isppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwarr(ikpt1)
   pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,5,idir)
   do ipw = 1, npw_k1
     jpw = pwind_tmp(ipw)
     if (jpw > 0) then
       do iband = 1, dtefield%mband_occ
         wfr = cg(1,icg1 + (iband - 1)*npw_k2*nspinor + jpw)
         wfi = cg(2,icg1 + (iband - 1)*npw_k2*nspinor + jpw)
         do  jband = 1, dtefield%mband_occ
           grad_berry(1,ipw,jband) = &
&           grad_berry(1,ipw,jband) + &
&           fac*qmat(1,iband,jband,ikpt,1,idir)*wfr - fac*qmat(2,iband,jband,ikpt,1,idir)*wfi
           grad_berry(2,ipw,jband) = &
&           grad_berry(2,ipw,jband) + &
&           fac*qmat(1,iband,jband,ikpt,1,idir)*wfi + fac*qmat(2,iband,jband,ikpt,1,idir)*wfr
         end do
       end do
     end if
   end do

!  prepare
   ikpt1 = dtefield%ikpt_dk(ikpt,2,idir)
   icg1 = dtefield%cgindex(ikpt1,isppol)
   npw_k1 = npwar1(ikpt)
   npw_k2 = npwarr(ikpt1)
   pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,6,idir)

   do ipw = 1, npw_k1
     jpw = pwind_tmp(ipw)
     if (jpw > 0) then
       do iband = 1, dtefield%mband_occ
         wfr = cg(1,icg1 + (iband - 1)*npw_k2*nspinor + jpw)
         wfi = cg(2,icg1 + (iband - 1)*npw_k2*nspinor + jpw)
         do  jband = 1, dtefield%mband_occ
           grad_berry(1,ipw,jband) = &
&           grad_berry(1,ipw,jband) - &
&           fac*qmat(1,iband,jband,ikpt,2,idir)*wfr + fac*qmat(2,iband,jband,ikpt,2,idir)*wfi

           grad_berry(2,ipw,jband) = &
&           grad_berry(2,ipw,jband) - &
&           fac*qmat(1,iband,jband,ikpt,2,idir)*wfi - fac*qmat(2,iband,jband,ikpt,2,idir)*wfr

         end do
       end do
     end if
   end do

 end do !idir

 ABI_DEALLOCATE(vect1)
 ABI_DEALLOCATE(vect2)
 ABI_DEALLOCATE(s1mat)
 ABI_DEALLOCATE(Amat)
 ABI_DEALLOCATE(Bmat)
 ABI_DEALLOCATE(pwind_tmp)

end subroutine dfptff_gbefd
!!***

!!****f* ABINIT/dfptff_edie
!! NAME
!! dfptff_edie
!!
!! FUNCTION
!! calculate the second order energy from the contribution of \Omega E \cdot P
!! term, Eq.(6) in PRB 75, 115116(2007) [[cite:Wang2007]].
!!
!! INPUTS
!! cg(2,mpw*nspinor*mband*mkmem*nsppol) = planewave coefficients of wavefunctions
!! cg1(2,mpw1*nspinor*mband*mk1mem*nsppol) = pw coefficients of
!! RF wavefunctions at k,q.
!! dtefield = variables related to response Berry-phase calculation
!! ikpt = the index of the current k point
!! isppol = the index of the spin component
!! mband =  maximum number of bands
!! mkmem = maximum number of k-points in core memory
!! mpw = maximum number of plane waves
!! mpw1 = maximum number of plane waves for response wavefunctions
!! nkpt = number of k points
!! npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!! npwar1(nkpt) = number of planewaves in basis and boundary for response wfs
!! nspinor = 1 for scalar wfs, 2 for spinor wfs
!! nsppol = 1 for unpolarized, 2 for spin-polarized
!! qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3) =
!! inverse of the overlap matrix
!! pwindall(max(mpw,mpw1)*mkmem,8,3) = array used to compute the overlap matrices
!! pwindall(:,1,:) <- <u^(0)_i|u^(0)_i+1>
!! pwindall(:,2,:) <- <u^(0)_i|u^(0)_i-1>
!! pwindall(:,3,:) <- <u^(1)_i|u^(1)_i+1>
!! pwindall(:,4,:) <- <u^(1)_i|u^(1)_i-1>
!! pwindall(:,5,:) <- <u^(1)_i|u^(0)_i+n+1>
!! pwindall(:,6,:) <- <u^(1)_i|u^(0)_i+n-1>
!! pwindall(:,7,:) <- <u^(0)_i|u^(1)_i-n+1>
!! pwindall(:,8,:) <- <u^(0)_i|u^(1)_i-n-1>
!!
!! OUTPUT
!! eberry = the energy of the Berry phase term
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      overlap_g
!!
!! SOURCE

subroutine dfptff_edie(cg,cg1,dtefield,eberry,idir_efield,mband,mkmem,&
&                mpw,mpw1,nkpt,npwarr,npwar1,nsppol,nspinor,pwindall,qmat,rprimd)

!Arguments ----------------------------------------
!scalars
 integer,intent(in) :: idir_efield,mband,mkmem,mpw,mpw1,nkpt,nspinor,nsppol
 real(dp),intent(out) :: eberry
 type(efield_type),intent(in) :: dtefield
!arrays
 integer,intent(in) :: npwar1(nkpt),npwarr(nkpt)
 integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
 real(dp),intent(in) :: cg(2,mpw*mband*mkmem*nspinor*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*mband*mkmem*nspinor*nsppol)
 real(dp),intent(in) :: qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables ----------------------------------
!scalars
 integer :: iband,icg,icg1,idir
 integer :: ikpt,ikpt1,ikptn,ikptnm
 integer :: jband,kband,mpw_tmp,npw_k1,npw_k2,pwmax,pwmin
 real(dp) :: doti,dotr,e0,fac
!arrays
 integer,allocatable :: pwind_tmp(:)
 real(dp) :: z1(2)
 real(dp),allocatable :: Amat(:,:,:),umat(:,:,:,:),vect1(:,:),vect2(:,:)

! *************************************************************************

!calculate 4 matrices -----------------------------
 mpw_tmp=max(mpw,mpw1)
 ABI_ALLOCATE(umat,(2,dtefield%mband_occ,dtefield%mband_occ,4))
 ABI_ALLOCATE(vect1,(2,0:mpw_tmp))
 ABI_ALLOCATE(vect2,(2,0:mpw_tmp))
 ABI_ALLOCATE(pwind_tmp,(mpw_tmp))
 ABI_ALLOCATE(Amat,(2,dtefield%mband_occ,dtefield%mband_occ))
 vect1(:,0) = zero ; vect2(:,0) = zero
 eberry=zero

 do ikpt=1,nkpt
   do idir=1,3
     fac = dtefield%efield_dot(idir)/&
&     (dble(dtefield%nstr(idir))*four_pi)

!    compute <u^(1)_{k_j,q}|u^(1)_{k_j+1,q}> matrix----------------------------------------------------

!    prepare to calculate overlap matrix
     ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
     icg = dtefield%cgindex(ikpt,1+nsppol)
     icg1 = dtefield%cgindex(ikpt1,1+nsppol)
     npw_k1 = npwar1(ikpt)
     npw_k2 = npwar1(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,3,idir)
     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
       do iband = 1, dtefield%mband_occ
         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg1(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         umat(1,iband,jband,1) = dotr
         umat(2,iband,jband,1) = doti
       end do    ! iband
     end do    !jband

!    compute <u^(0)_{k_j}|u^(1)_{k_j-n+1,q}> matrix----------------------------------------------------

!    prepare to calculate overlap matrix
     ikpt1 = dtefield%ikpt_dk(ikpt,5,idir)
     icg = dtefield%cgindex(ikpt,1)
     icg1 = dtefield%cgindex(ikpt1,1+nsppol)
     npw_k1 = npwarr(ikpt)
     npw_k2 = npwar1(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,7,idir)

     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
       do iband = 1, dtefield%mband_occ
         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         umat(1,iband,jband,2) = dotr
         umat(2,iband,jband,2) = doti
       end do    ! iband
     end do    !jband

!    compute <u^(1)_{k_j-n,q}|u^(0)_{k_j+1}> matrix----------------------------------------------------

!    prepare to calculate overlap matrix
     ikptn  = dtefield%ikpt_dk(ikpt,8,idir)
     ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
     icg = dtefield%cgindex(ikptn,1+nsppol)
     icg1 = dtefield%cgindex(ikpt1,1)
     npw_k1 = npwar1(ikptn)
     npw_k2 = npwarr(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikptn-1)*mpw_tmp+1:(ikptn-1)*mpw_tmp+npw_k1,5,idir)

     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
       do iband = 1, dtefield%mband_occ
         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg1(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         umat(1,iband,jband,3) = dotr
         umat(2,iband,jband,3) = doti
       end do    ! iband
     end do    !jband

!    compute <u^(0)_{-k_j-n+1}|u^(1)_{-k_j,q}> matrix----------------------------------------------------

!    prepare to calculate overlap matrix
     ikptn = dtefield%ikpt_dk(ikpt,5,idir)
     ikptnm = dtefield%ikpt_dk(ikptn,9,idir)
     ikpt1 = dtefield%ikpt_dk(ikpt,9,idir)
     icg = dtefield%cgindex(ikptnm,1)
     icg1 = dtefield%cgindex(ikpt1,1+nsppol)
     npw_k1 = npwarr(ikptnm)
     npw_k2 = npwar1(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikptnm-1)*mpw_tmp+1:(ikptnm-1)*mpw_tmp+npw_k1,7,idir)
     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
       do iband = 1, dtefield%mband_occ
         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         umat(1,iband,jband,4) = dotr
         umat(2,iband,jband,4) = doti
       end do    ! iband
     end do    !jband

!    sum over the whole------------------------------------------------------------

     e0=zero
     do iband=1,dtefield%mband_occ
       do jband=1,dtefield%mband_occ
         e0 = e0 + 2_dp*(umat(1,iband,jband,1)*qmat(2,jband,iband,ikpt,1,idir)&
&         +       umat(2,iband,jband,1)*qmat(1,jband,iband,ikpt,1,idir))
       end do
     end do
     eberry = eberry - e0*fac
     e0=zero
     ikptn=dtefield%ikpt_dk(ikpt,8,idir)
     Amat(:,:,:)=zero
!    calculate Amat
     do iband=1, dtefield%mband_occ
       do jband=1, dtefield%mband_occ
         do kband=1, dtefield%mband_occ
           Amat(1,iband,jband) = Amat(1,iband,jband) + (umat(1,iband,kband,3))*qmat(1,kband,jband,ikpt,1,idir)&
&           - (umat(2,iband,kband,3))*qmat(2,kband,jband,ikpt,1,idir)
           Amat(2,iband,jband) = Amat(2,iband,jband) + (umat(1,iband,kband,3))*qmat(2,kband,jband,ikpt,1,idir)&
&           + (umat(2,iband,kband,3))*qmat(1,kband,jband,ikpt,1,idir)
         end do
       end do
     end do

     do iband=1, dtefield%mband_occ
       do jband=1, dtefield%mband_occ
         do kband=1, dtefield%mband_occ
           z1(1) = (umat(1,jband,iband,4)+umat(1,iband,jband,2))*qmat(1,jband,kband,ikptn,1,idir)&
&           - (umat(2,jband,iband,4)+umat(2,iband,jband,2))*qmat(2,jband,kband,ikptn,1,idir)
           z1(2) = (umat(1,jband,iband,4)+umat(1,iband,jband,2))*qmat(2,jband,kband,ikptn,1,idir)&
&           + (umat(2,jband,iband,4)+umat(2,iband,jband,2))*qmat(1,jband,kband,ikptn,1,idir)

           e0 = e0 - 4_dp*(z1(1)*Amat(2,kband,iband)+z1(2)*Amat(1,kband,iband))
         end do
       end do
     end do

     eberry = eberry - e0*fac

!    !---------------------------------last part---------------------------------------------

     fac = rprimd(idir_efield,idir)/&
&     (dble(dtefield%nstr(idir))*two_pi)

!    compute <u^(1)_{k_j-n,q}|u^(0)_{k_j+1}> matrix----------------------------------------------------

!    prepare to calculate overlap matrix
     ikptn  = dtefield%ikpt_dk(ikpt,8,idir)
     ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
     icg = dtefield%cgindex(ikptn,1+nsppol)
     icg1 = dtefield%cgindex(ikpt1,1)
     npw_k1 = npwar1(ikptn)
     npw_k2 = npwarr(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikptn-1)*mpw_tmp+1:(ikptn-1)*mpw_tmp+npw_k1,5,idir)

     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
       do iband = 1, dtefield%mband_occ
         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg1(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         umat(1,iband,jband,1) = dotr
         umat(2,iband,jband,1) = doti
       end do    ! iband
     end do    !jband

!    compute <u^(0)_{k_j}|u^(1)_{k_j-n+1,q}> matrix----------------------------------------------------

!    prepare to calculate overlap matrix
     ikpt1 = dtefield%ikpt_dk(ikpt,5,idir)
     icg = dtefield%cgindex(ikpt,1)
     icg1 = dtefield%cgindex(ikpt1,1+nsppol)
     npw_k1 = npwarr(ikpt)
     npw_k2 = npwar1(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,7,idir)
     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
       do iband = 1, dtefield%mband_occ
         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         umat(1,iband,jband,1) = umat(1,iband,jband,1) + dotr
         umat(2,iband,jband,1) = umat(2,iband,jband,1) + doti
       end do    ! iband
     end do    !jband

     e0=zero

     do iband=1,dtefield%mband_occ
       do jband=1,dtefield%mband_occ
         e0 = e0 +    (umat(1,iband,jband,1)*qmat(2,jband,iband,ikpt,1,idir)&
&         +       umat(2,iband,jband,1)*qmat(1,jband,iband,ikpt,1,idir))
       end do
     end do

     eberry = eberry - e0*fac

   end do !end idir
 end do !end ikpt

 ABI_DEALLOCATE(umat)
 ABI_DEALLOCATE(vect1)
 ABI_DEALLOCATE(vect2)
 ABI_DEALLOCATE(pwind_tmp)
 ABI_DEALLOCATE(Amat)

end subroutine dfptff_edie
!!***

!!****f* ABINIT/dfptff_ebp
!! NAME
!! dfptff_ebp
!!
!! FUNCTION
!! calculation of the energy from the term \Omega E \cdot P
!!
!! INPUTS
!! cg(2,mpw*nspinor*mband*mkmem*nsppol) = planewave coefficients of wavefunctions
!! cg1(2,mpw1*nspinor*mband*mk1mem*nsppol) = pw coefficients of
!! RF wavefunctions at k,q.
!! dtefield = variables related to response Berry-phase calculation
!! ikpt = the index of the current k point
!! isppol = the index of the spin component
!! mband =  maximum number of bands
!! mkmem = maximum number of k-points in core memory
!! mpw = maximum number of plane waves
!! mpw1 = maximum number of plane waves for response wavefunctions
!! nkpt = number of k points
!! npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!! npwar1(nkpt) = number of planewaves in basis and boundary for response wfs
!! nspinor = 1 for scalar wfs, 2 for spinor wfs
!! nsppol = 1 for unpolarized, 2 for spin-polarized
!! qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3) =
!! inverse of the overlap matrix
!! pwindall(max(mpw,mpw1)*mkmem,8,3) = array used to compute the overlap matrices
!! pwindall(:,1,:) <- <u^(0)_i|u^(0)_i+1>
!! pwindall(:,2,:) <- <u^(0)_i|u^(0)_i-1>
!! pwindall(:,3,:) <- <u^(1)_i|u^(1)_i+1>
!! pwindall(:,4,:) <- <u^(1)_i|u^(1)_i-1>
!! pwindall(:,5,:) <- <u^(1)_i|u^(0)_i+n+1>
!! pwindall(:,6,:) <- <u^(1)_i|u^(0)_i+n-1>
!! pwindall(:,7,:) <- <u^(0)_i|u^(1)_i-n+1>
!! pwindall(:,8,:) <- <u^(0)_i|u^(1)_i-n-1>
!!
!! OUTPUT
!! grad_berry(2,mpw1,dtefield%mband_occ) = the gradient of the Berry phase term
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      overlap_g
!!
!! SOURCE

subroutine dfptff_ebp(cg,cg1,dtefield,eberry,mband,mkmem,&
&               mpw,mpw1,nkpt,npwarr,npwar1,nsppol,nspinor,pwindall,qmat)

!Arguments ----------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,mpw,mpw1,nkpt,nspinor,nsppol
 real(dp),intent(out) :: eberry
 type(efield_type),intent(in) :: dtefield
!arrays
 integer,intent(in) :: npwar1(nkpt),npwarr(nkpt)
 integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
 real(dp),intent(in) :: cg(2,mpw*mband*mkmem*nspinor*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*mband*mkmem*nspinor*nsppol)
 real(dp),intent(in) :: qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3)

!Local variables ----------------------------------
!scalars
 integer :: iband,icg,icg1,idir
 integer :: ikpt,ikpt1,ikptn,ikptnm
 integer :: jband,kband,mpw_tmp,npw_k1,npw_k2,pwmax,pwmin
 real(dp) :: doti,dotr,e0,fac
!arrays
 integer,allocatable :: pwind_tmp(:)
 real(dp) :: z1(2)
 real(dp),allocatable :: Amat(:,:,:),umat(:,:,:,:),vect1(:,:),vect2(:,:)

! *************************************************************************

!calculate 4 matrices -----------------------------
 mpw_tmp=max(mpw,mpw1)
 ABI_ALLOCATE(umat,(2,dtefield%mband_occ,dtefield%mband_occ,4))
 ABI_ALLOCATE(vect1,(2,0:mpw_tmp))
 ABI_ALLOCATE(vect2,(2,0:mpw_tmp))
 ABI_ALLOCATE(pwind_tmp,(mpw_tmp))
 ABI_ALLOCATE(Amat,(2,dtefield%mband_occ,dtefield%mband_occ))
 vect1(:,0) = zero ; vect2(:,0) = zero
 eberry=zero

 do ikpt=1,nkpt

   do idir=1,3

     fac = dtefield%efield_dot(idir)/&
&     (dble(dtefield%nstr(idir))*four_pi)

!    compute <u^(1)_{k_j,q}|u^(1)_{k_j+1,q}> matrix---------

!    prepare to calculate overlap matrix
     ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
     icg = dtefield%cgindex(ikpt,1+nsppol)
     icg1 = dtefield%cgindex(ikpt1,1+nsppol)
     npw_k1 = npwar1(ikpt)
     npw_k2 = npwar1(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,3,idir)
     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)

       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero

       do iband = 1, dtefield%mband_occ

         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg1(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         umat(1,iband,jband,1) = dotr
         umat(2,iband,jband,1) = doti

       end do    ! iband
     end do    !jband

!    compute <u^(0)_{k_j}|u^(1)_{k_j-n+1,q}> matrix----------------------------------------------------

!    prepare to calculate overlap matrix
     ikpt1 = dtefield%ikpt_dk(ikpt,5,idir)
     icg = dtefield%cgindex(ikpt,1)
     icg1 = dtefield%cgindex(ikpt1,1+nsppol)
     npw_k1 = npwarr(ikpt)
     npw_k2 = npwar1(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,7,idir)

     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero

       do iband = 1, dtefield%mband_occ

         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         umat(1,iband,jband,2) = dotr
         umat(2,iband,jband,2) = doti

       end do    ! iband
     end do    !jband

!    compute <u^(1)_{k_j-n,q}|u^(0)_{k_j+1}> matrix----------------------------------------------------

!    prepare to calculate overlap matrix
     ikptn  = dtefield%ikpt_dk(ikpt,8,idir)
     ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
     icg = dtefield%cgindex(ikptn,1+nsppol)
     icg1 = dtefield%cgindex(ikpt1,1)
     npw_k1 = npwar1(ikptn)
     npw_k2 = npwarr(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikptn-1)*mpw_tmp+1:(ikptn-1)*mpw_tmp+npw_k1,5,idir)

     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero

       do iband = 1, dtefield%mband_occ

         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg1(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         umat(1,iband,jband,3) = dotr
         umat(2,iband,jband,3) = doti

       end do    ! iband
     end do    !jband

!    compute <u^(0)_{-k_j-n+1}|u^(1)_{-k_j,q}> matrix----------------------------------------------------

!    prepare to calculate overlap matrix
     ikptn = dtefield%ikpt_dk(ikpt,5,idir)
     ikptnm = dtefield%ikpt_dk(ikptn,9,idir)
     ikpt1 = dtefield%ikpt_dk(ikpt,9,idir)
     icg = dtefield%cgindex(ikptnm,1)
     icg1 = dtefield%cgindex(ikpt1,1+nsppol)
     npw_k1 = npwarr(ikptnm)
     npw_k2 = npwar1(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikptnm-1)*mpw_tmp+1:(ikptnm-1)*mpw_tmp+npw_k1,7,idir)

     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero

       do iband = 1, dtefield%mband_occ

         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         umat(1,iband,jband,4) = dotr
         umat(2,iband,jband,4) = doti

       end do    ! iband
     end do    !jband

!    sum over the whole------------------------------------------------------------

     e0=zero
     do iband=1,dtefield%mband_occ
       do jband=1,dtefield%mband_occ
         e0 = e0 + 4_dp*(umat(1,iband,jband,1)*qmat(2,jband,iband,ikpt,1,idir)&
&         +       umat(2,iband,jband,1)*qmat(1,jband,iband,ikpt,1,idir))

       end do
     end do

     eberry = eberry - e0*fac

     e0=zero

     ikptn=dtefield%ikpt_dk(ikpt,8,idir)

     Amat(:,:,:)=zero

!    calculate Amat
     do iband=1, dtefield%mband_occ
       do jband=1, dtefield%mband_occ
         do kband=1, dtefield%mband_occ
           Amat(1,iband,jband) = Amat(1,iband,jband) + (umat(1,iband,kband,3))*&
&           qmat(1,kband,jband,ikpt,1,idir)&
&           - (umat(2,iband,kband,3))*qmat(2,kband,jband,ikpt,1,idir)
           Amat(2,iband,jband) = Amat(2,iband,jband) + (umat(1,iband,kband,3))*&
&           qmat(2,kband,jband,ikpt,1,idir)&
&           + (umat(2,iband,kband,3))*qmat(1,kband,jband,ikpt,1,idir)
         end do
       end do
     end do

     do iband=1, dtefield%mband_occ
       do jband=1, dtefield%mband_occ
         do kband=1, dtefield%mband_occ

           z1(1) = (umat(1,jband,iband,4)+umat(1,iband,jband,2))*&
&           qmat(1,jband,kband,ikptn,1,idir)&
&           -    (umat(2,jband,iband,4)+umat(2,iband,jband,2))*&
&           qmat(2,jband,kband,ikptn,1,idir)
           z1(2) = (umat(1,jband,iband,4)+umat(1,iband,jband,2))*&
&           qmat(2,jband,kband,ikptn,1,idir)&
&           +    (umat(2,jband,iband,4)+umat(2,iband,jband,2))*&
&           qmat(1,jband,kband,ikptn,1,idir)

           e0 = e0 - 4_dp*(z1(1)*Amat(2,kband,iband)+z1(2)*Amat(1,kband,iband))

         end do
       end do
     end do

     eberry = eberry - e0*fac

   end do !end idir
 end do !end ikpt

 ABI_DEALLOCATE(umat)
 ABI_DEALLOCATE(vect1)
 ABI_DEALLOCATE(vect2)
 ABI_DEALLOCATE(pwind_tmp)
 ABI_DEALLOCATE(Amat)

end subroutine dfptff_ebp
!!***

!!****f* ABINIT/dfptff_die
!! NAME
!! dfptff_die
!!
!! FUNCTION
!! calculate electric susceptibility tensor in Eq.(28) in PRB 75, 115116(2007) [[cite:Wang2007]].
!!
!! INPUTS
!! cg(2,mpw*nspinor*mband*mkmem*nsppol) = planewave coefficients of wavefunctions
!! cg1(2,mpw1*nspinor*mband*mk1mem*nsppol) = pw coefficients of
!! RF wavefunctions at k,q.
!! dtefield = variables related to response Berry-phase calculation
!! idirpert = the current coloumn of the dielectric permittivity tensor
!! mband =  maximum number of bands
!! mkmem = maximum number of k-points in core memory
!! mpw = maximum number of plane waves
!! mpw1 = maximum number of plane waves for response wavefunctions
!! nkpt = number of k points
!! npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!! npwar1(nkpt) = number of planewaves in basis and boundary for response wfs
!! nspinor = 1 for scalar wfs, 2 for spinor wfs
!! nsppol = 1 for unpolarized, 2 for spin-polarized
!! qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3) =
!! inverse of the overlap matrix
!! pwindall(max(mpw,mpw1)*mkmem,8,3) = array used to compute the overlap matrices
!! pwindall(:,1,:) <- <u^(0)_i|u^(0)_i+1>
!! pwindall(:,2,:) <- <u^(0)_i|u^(0)_i-1>
!! pwindall(:,3,:) <- <u^(1)_i|u^(1)_i+1>
!! pwindall(:,4,:) <- <u^(1)_i|u^(1)_i-1>
!! pwindall(:,5,:) <- <u^(1)_i|u^(0)_i+n+1>
!! pwindall(:,6,:) <- <u^(1)_i|u^(0)_i+n-1>
!! pwindall(:,7,:) <- <u^(0)_i|u^(1)_i-n+1>
!! pwindall(:,8,:) <- <u^(0)_i|u^(1)_i-n-1>
!! rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! OUTPUT
!! diet = electric susceptibility tensor
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      overlap_g
!!
!! SOURCE

subroutine dfptff_die(cg,cg1,dtefield,d2lo,idirpert,ipert,mband,mkmem,&
&               mpw,mpw1,mpert,nkpt,npwarr,npwar1,nsppol,nspinor,pwindall,qmat,rprimd)

!Arguments ----------------------------------------
!scalars
 integer,intent(in) :: idirpert,ipert,mband,mkmem,mpert,mpw,mpw1,nkpt,nspinor
 integer,intent(in) :: nsppol
 type(efield_type),intent(in) :: dtefield
!arrays
 integer,intent(in) :: npwar1(nkpt),npwarr(nkpt)
 integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
 real(dp),intent(in) :: cg(2,mpw*mband*mkmem*nspinor*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*mband*mkmem*nspinor*nsppol)
 real(dp),intent(in) :: qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: d2lo(2,3,mpert,3,mpert) !vz_i

!Local variables ----------------------------------
!scalars
 integer :: ialpha,iband,icg,icg1,idir,ikpt,ikpt1,jband,mpw_tmp,npw_k1
 integer :: npw_k2,pwmax,pwmin
 real(dp) :: doti,dotr,e0,fac
!arrays
 integer,allocatable :: pwind_tmp(:)
 real(dp) :: edir(3)
 real(dp),allocatable :: s1mat(:,:,:),vect1(:,:),vect2(:,:)

! *************************************************************************

!calculate s1 matrices -----------------------------
 mpw_tmp=max(mpw,mpw1)
 ABI_ALLOCATE(s1mat,(2,dtefield%mband_occ,dtefield%mband_occ))
 ABI_ALLOCATE(vect1,(2,0:mpw_tmp))
 ABI_ALLOCATE(vect2,(2,0:mpw_tmp))
 ABI_ALLOCATE(pwind_tmp,(mpw_tmp))
 vect1(:,0) = zero ; vect2(:,0) = zero

 edir(:)=zero

 do ikpt=1,nkpt
   do idir=1,3
!    compute <u^(0)_{k_j}|u^(1)_{k_j+1,q}> matrix--- q=0 ----------------------------------------

!    prepare to calculate overlap matrix
     ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
     icg = dtefield%cgindex(ikpt,1)
     icg1 = dtefield%cgindex(ikpt1,1+nsppol)
     npw_k1 = npwarr(ikpt)
     npw_k2 = npwar1(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,7,idir)

     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
       do iband = 1, dtefield%mband_occ
         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         s1mat(1,iband,jband) = dotr
         s1mat(2,iband,jband) = doti
       end do    ! iband
     end do    !jband

!    compute <u^(1)_{k_j,q}|u^(0)_{k_j+1}> matrix-- q=0 -------------------------------------

!    prepare to calculate overlap matrix
     ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
     icg = dtefield%cgindex(ikpt,1+nsppol)
     icg1 = dtefield%cgindex(ikpt1,1)
     npw_k1 = npwar1(ikpt)
     npw_k2 = npwarr(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,5,idir)
     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
       do iband = 1, dtefield%mband_occ
         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg1(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         s1mat(1,iband,jband) = s1mat(1,iband,jband) + dotr
         s1mat(2,iband,jband) = s1mat(2,iband,jband) + doti
       end do    ! iband
     end do    !jband

!    sum over the whole------------------------------------------------------------

     e0=zero

     do iband=1,dtefield%mband_occ
       do jband=1,dtefield%mband_occ
         e0 = e0 + (s1mat(1,iband,jband)*qmat(2,jband,iband,ikpt,1,idir)&
&         +    s1mat(2,iband,jband)*qmat(1,jband,iband,ikpt,1,idir))

       end do
     end do

     do ialpha=1,3
       fac = rprimd(ialpha,idir)/&
&       (dble(dtefield%nstr(idir))*pi)
       edir(ialpha)=edir(ialpha)+ e0*fac
     end do

   end do !idir
 end do !ikpt

 d2lo(1,1:3,ipert,idirpert,ipert)=edir(:)

 ABI_DEALLOCATE(s1mat)
 ABI_DEALLOCATE(vect1)
 ABI_DEALLOCATE(vect2)
 ABI_DEALLOCATE(pwind_tmp)

end subroutine dfptff_die
!!***

!!****f* ABINIT/dfptff_bec
!! NAME
!! dfptff_bec
!!
!! FUNCTION
!! calculate Born effective charge tensor in Eq.(33) in PRB 75, 115116(2007) [[cite:Wang2007]].
!!
!! INPUTS
!! cg(2,mpw*nspinor*mband*mkmem*nsppol) = planewave coefficients of wavefunctions
!! cg1(2,mpw1*nspinor*mband*mk1mem*nsppol) = pw coefficients of
!! RF wavefunctions at k,q.
!! dtefield = variables related to response Berry-phase calculation
!! idirpert = the current coloumn of the dielectric permittivity tensor
!! mband =  maximum number of bands
!! mkmem = maximum number of k-points in core memory
!! mpw = maximum number of plane waves
!! mpw1 = maximum number of plane waves for response wavefunctions
!! nkpt = number of k points
!! npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!! npwar1(nkpt) = number of planewaves in basis and boundary for response wfs
!! nspinor = 1 for scalar wfs, 2 for spinor wfs
!! nsppol = 1 for unpolarized, 2 for spin-polarized
!! qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3) =
!! inverse of the overlap matrix
!! pwindall(max(mpw,mpw1)*mkmem,8,3) = array used to compute the overlap matrices
!! pwindall(:,1,:) <- <u^(0)_i|u^(0)_i+1>
!! pwindall(:,2,:) <- <u^(0)_i|u^(0)_i-1>
!! pwindall(:,3,:) <- <u^(1)_i|u^(1)_i+1>
!! pwindall(:,4,:) <- <u^(1)_i|u^(1)_i-1>
!! pwindall(:,5,:) <- <u^(1)_i|u^(0)_i+n+1>
!! pwindall(:,6,:) <- <u^(1)_i|u^(0)_i+n-1>
!! pwindall(:,7,:) <- <u^(0)_i|u^(1)_i-n+1>
!! pwindall(:,8,:) <- <u^(0)_i|u^(1)_i-n-1>
!! rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! OUTPUT
!! d2lo(1,1:3,natom+5,1:3,1:natom) = Born effective charge tensor
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      overlap_g
!!
!! SOURCE

subroutine dfptff_bec(cg,cg1,dtefield,natom,d2lo,idirpert,ipert,mband,mkmem,&
&               mpw,mpw1,mpert,nkpt,npwarr,npwar1,nsppol,nspinor,pwindall,qmat,rprimd)

!Arguments ----------------------------------------
!scalars
 integer,intent(in) :: idirpert,ipert,mband,mkmem,mpert,mpw,mpw1,natom,nkpt
 integer,intent(in) :: nspinor,nsppol
 type(efield_type),intent(in) :: dtefield
!arrays
 integer,intent(in) :: npwar1(nkpt),npwarr(nkpt)
 integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
 real(dp),intent(in) :: cg(2,mpw*mband*mkmem*nspinor*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*mband*mkmem*nspinor*nsppol)
 real(dp),intent(in) :: qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: d2lo(2,3,mpert,3,mpert) !vz_i

!Local variables ----------------------------------
!scalars
 integer :: ialpha,iband,icg,icg1,idir,ikpt,ikpt1,jband,mpw_tmp,npw_k1
 integer :: npw_k2,pwmax,pwmin
 real(dp) :: doti,dotr,e0,fac
!arrays
 integer,allocatable :: pwind_tmp(:)
 real(dp) :: edir(3)
 real(dp),allocatable :: s1mat(:,:,:),vect1(:,:),vect2(:,:)

! *************************************************************************

!calculate s1 matrices -----------------------------
 mpw_tmp=max(mpw,mpw1)
 ABI_ALLOCATE(s1mat,(2,dtefield%mband_occ,dtefield%mband_occ))
 ABI_ALLOCATE(vect1,(2,0:mpw_tmp))
 ABI_ALLOCATE(vect2,(2,0:mpw_tmp))
 ABI_ALLOCATE(pwind_tmp,(mpw_tmp))
 vect1(:,0) = zero ; vect2(:,0) = zero

 edir(:)=zero

 do ikpt=1,nkpt
   do idir=1,3

!    compute <u^(0)_{k_j}|u^(1)_{k_j+1,q}> matrix--- q=0 ----------------------------------------

!    prepare to calculate overlap matrix
     ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
     icg = dtefield%cgindex(ikpt,1)
     icg1 = dtefield%cgindex(ikpt1,1+nsppol)
     npw_k1 = npwarr(ikpt)
     npw_k2 = npwar1(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,7,idir)
     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg1(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
       do iband = 1, dtefield%mband_occ
         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         s1mat(1,iband,jband) = dotr
         s1mat(2,iband,jband) = doti
       end do    ! iband
     end do    !jband

!    compute <u^(1)_{k_j,q}|u^(0)_{k_j+1}> matrix-- q=0 -------------------------------------

!    prepare to calculate overlap matrix
     ikpt1 = dtefield%ikpt_dk(ikpt,1,idir)
     icg = dtefield%cgindex(ikpt,1+nsppol)
     icg1 = dtefield%cgindex(ikpt1,1)
     npw_k1 = npwar1(ikpt)
     npw_k2 = npwarr(ikpt1)
     pwind_tmp(1:npw_k1) = pwindall((ikpt-1)*mpw_tmp+1:(ikpt-1)*mpw_tmp+npw_k1,5,idir)
     vect1(:,0) = zero ; vect2(:,0) = zero
     do jband = 1, dtefield%mband_occ
       vect2(:,1:npw_k2) = &
&       cg(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
       if (npw_k2 < mpw_tmp) vect2(:,npw_k2+1:mpw_tmp) = zero
       do iband = 1, dtefield%mband_occ
         pwmin = (iband-1)*npw_k1*nspinor
         pwmax = pwmin + npw_k1*nspinor
         vect1(:,1:npw_k1) = &
&         cg1(:,icg + 1 + pwmin:icg + pwmax)
         if (npw_k1 < mpw_tmp) vect1(:,npw_k1+1:mpw_tmp) = zero
         call overlap_g(doti,dotr,mpw_tmp,npw_k1,npw_k2,nspinor,pwind_tmp,&
&         vect1,vect2)
         s1mat(1,iband,jband) = s1mat(1,iband,jband) + dotr
         s1mat(2,iband,jband) = s1mat(2,iband,jband) + doti
       end do    ! iband
     end do    !jband

!    sum over the whole------------------------------------------------------------

     e0=zero

     do iband=1,dtefield%mband_occ
       do jband=1,dtefield%mband_occ
         e0 = e0 + (s1mat(1,iband,jband)*qmat(2,jband,iband,ikpt,1,idir)&
&         +    s1mat(2,iband,jband)*qmat(1,jband,iband,ikpt,1,idir))
       end do
     end do

     do ialpha=1,3
       fac = rprimd(ialpha,idir)/&
&       (dble(dtefield%nstr(idir))*pi)

       edir(ialpha)=edir(ialpha)+ e0*fac
     end do

   end do
 end do

 d2lo(1,1:3,natom+2,idirpert,ipert)=edir(:)

 ABI_DEALLOCATE(s1mat)
 ABI_DEALLOCATE(vect1)
 ABI_DEALLOCATE(vect2)
 ABI_DEALLOCATE(pwind_tmp)

end subroutine dfptff_bec
!!***

!!****f* ABINIT/qmatrix
!! NAME
!! qmatrix
!!
!! FUNCTION
!! calculation of the inverse of the overlap matrix
!!
!! INPUTS
!! cg(2,mpw*nspinor*mband*mkmem*nsppol) = planewave coefficients of wavefunctions
!! RF wavefunctions at k,q.
!! dtefield = variables related to response Berry-phase calculation
!! ikpt = the index of the current k point
!! mband =  maximum number of bands
!! mkmem = maximum number of k-points in core memory
!! mpw = maximum number of plane waves
!! mpw1 = maximum number of plane waves for response wavefunctions
!! nkpt = number of k points
!! npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!! npwar1(nkpt) = number of planewaves in basis and boundary for response wfs
!! nspinor = 1 for scalar wfs, 2 for spinor wfs
!! nsppol = 1 for unpolarized, 2 for spin-polarized
!! pwindall(max(mpw,mpw1)*mkmem,8,3) = array used to compute the overlap matrices
!!
!! OUTPUT
!! qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3) = inverse of the overlap matrix
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      dzgedi,dzgefa,overlap_g
!!
!! SOURCE

subroutine qmatrix(cg,dtefield,qmat,mpw,mpw1,mkmem,mband,npwarr,nkpt,nspinor,nsppol,pwindall)

 use m_hide_lapack, only : dzgedi, dzgefa

!Arguments ----------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,mpw,mpw1,nkpt,nspinor,nsppol
 type(efield_type),intent(in) :: dtefield
!arrays
 integer,intent(in) :: npwarr(nkpt),pwindall(max(mpw,mpw1)*mkmem,8,3)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(out) :: qmat(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3)

!Local variables -------------------------
!scalars
 integer :: iband,icg,icg1,idir,ifor,ikpt,ikpt2,info,jband,job
 integer :: npw_k1,npw_k2,pwmax,pwmin
 integer :: isppol
 real(dp) :: doti,dotr
!arrays
 integer,allocatable :: ipvt(:),pwind_k(:)
 real(dp) :: det(2,2)
 real(dp),allocatable :: sinv(:,:,:),smat_k(:,:,:),vect1(:,:),vect2(:,:)
 real(dp),allocatable :: zgwork(:,:)

! *************************************************************************

 ABI_ALLOCATE(ipvt,(dtefield%mband_occ))
 ABI_ALLOCATE(sinv,(2,dtefield%mband_occ,dtefield%mband_occ))
 ABI_ALLOCATE(zgwork,(2,dtefield%mband_occ))
 ABI_ALLOCATE(vect1,(2,0:mpw))
 ABI_ALLOCATE(vect2,(2,0:mpw))
 ABI_ALLOCATE(smat_k,(2,dtefield%mband_occ,dtefield%mband_occ))
 ABI_ALLOCATE(pwind_k,(max(mpw,mpw1)))
 vect1(:,0) = zero ; vect2(:,0) = zero

 job = 11

!loop over k points
 do isppol = 1, nsppol
   do ikpt = 1, nkpt
     npw_k1 = npwarr(ikpt)
     icg  = dtefield%cgindex(ikpt,1)
     do idir = 1, 3
       do ifor = 1, 2

         ikpt2 = dtefield%ikpt_dk(ikpt,ifor,idir)
         npw_k2 = npwarr(ikpt2)
         icg1 = dtefield%cgindex(ikpt2,1)
         pwind_k(1:npw_k1) = pwindall((ikpt-1)*max(mpw,mpw1)+1:(ikpt-1)*max(mpw,mpw1)+npw_k1,ifor,idir)

         do jband = 1, dtefield%nband_occ(isppol)
           vect2(:,1:npw_k2) = cg(:,icg1 + 1 + (jband-1)*npw_k2*nspinor:icg1 + jband*npw_k2*nspinor)
           if (npw_k2 < mpw) vect2(:,npw_k2+1:mpw) = zero

           do iband = 1, dtefield%nband_occ(isppol)

             pwmin = (iband-1)*npw_k1*nspinor
             pwmax = pwmin + npw_k1*nspinor
             vect1(:,1:npw_k1) = cg(:,icg + 1 + pwmin:icg + pwmax)
             if (npw_k1 < mpw) vect1(:,npw_k1+1:mpw) = zero
             call overlap_g(doti,dotr,mpw,npw_k1,npw_k2,nspinor,pwind_k,vect1,vect2)
             smat_k(1,iband,jband) = dotr
             smat_k(2,iband,jband) = doti
           end do    ! iband
         end do    !jband

         sinv(:,:,:) = smat_k(:,:,:)

         call dzgefa(sinv,dtefield%mband_occ,dtefield%nband_occ(isppol),ipvt,info)
         call dzgedi(sinv,dtefield%mband_occ,dtefield%nband_occ(isppol),ipvt,det,zgwork,job)

         qmat(:,:,:,ikpt,ifor,idir) = sinv(:,:,:)
       end do
     end do
   end do  !end loop over k
 end do

 ABI_DEALLOCATE(ipvt)
 ABI_DEALLOCATE(sinv)
 ABI_DEALLOCATE(zgwork)
 ABI_DEALLOCATE(vect1)
 ABI_DEALLOCATE(vect2)
 ABI_DEALLOCATE(smat_k)
 ABI_DEALLOCATE(pwind_k)

end subroutine qmatrix
!!***

end module m_dfpt_fef
!!***
