!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfptff_initberry
!! NAME
!! dfptff_initberry
!!
!! FUNCTION
!! Initialization of response calculations in finite electric field. 
!!
!! COPYRIGHT
!! Copyright (C) 2004-2016 ABINIT group (XW).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      kpgio,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine  dfptff_initberry(dtefield,dtset,gmet,kg,kg1,mband,mkmem,mpi_enreg,&
&                mpw,mpw1,nkpt,npwarr,npwar1,nsppol,occ,pwindall,rprimd)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_efield

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptff_initberry'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ----------------------------------------
! scalars
! arrays
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
! scalars
! arrays
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
 pwindall(:,:,:) = 0

!Compute the number of occupied bands and check that--------
!it is the same for each k-point----------------------------

 occ_val = two/(dtset%nsppol*one)

 index = 0
 do ikpt = 1, nkpt

   mband_occ_k = 0
   nband_k = dtset%nband(ikpt)

   do iband = 1, nband_k
     index = index + 1
     if (abs(occ(index) - occ_val) < tol8) mband_occ_k = mband_occ_k + 1
   end do

   if (ikpt > 1) then
     if (dtefield%mband_occ /= mband_occ_k) then
       message = ' The number of valence bands is not the same for every k-point'
       MSG_ERROR(message)
     end if
   else
     dtefield%mband_occ = mband_occ_k
   end if

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
