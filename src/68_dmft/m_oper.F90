!!****m* ABINIT/m_oper
!! NAME
!!  m_oper
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

! nvtx related macro definition
#include "nvtx_macros.h"

MODULE m_oper

 use, intrinsic :: iso_c_binding, only: c_size_t, c_loc
 use defs_basis
 use m_abicore
 use m_errors
 use m_xomp
 use m_abi_linalg

 !use m_abi_linalg, only : abi_xgemm
 use m_hide_lapack, only : xginv
 use m_matlu, only : copy_matlu,destroy_matlu,diff_matlu,identity_matlu,init_matlu, &
               & inverse_matlu,matlu_type,print_matlu,prod_matlu,trace_matlu,zero_matlu
 use m_paw_dmft, only : mpi_distrib_dmft_type,paw_dmft_type
 use m_xmpi, only : xmpi_allgatherv,xmpi_gatherv,xmpi_sum,xmpi_sum_master

#ifdef HAVE_GPU_MARKERS
 use m_nvtx_data
#endif

 implicit none

 private

 public :: init_oper
 public :: init_oper_ndat
 public :: diff_oper
 public :: destroy_oper
 public :: print_oper
 public :: inverse_oper
 public :: downfold_oper
 public :: identity_oper
 public :: copy_oper
 public :: copy_oper_from_ndat
 public :: copy_oper_to_ndat
 public :: trace_oper
 public :: upfold_oper
 public :: prod_oper
 public :: trace_prod_oper
 public :: gather_oper
 public :: gather_oper_ks
!!***

!!****t* m_oper/oper_type
!! NAME
!!  oper_type
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! SOURCE

 type, public :: oper_type ! for each atom

!  integer :: maxlpawu         ! Number of correlated atoms
!
!  integer :: mband
!  ! Number of bands

  ! Wether ks and matlu are stored on GPU
  integer :: gpu_option

  integer :: ndat

  integer :: has_operks
  ! Is the operator allocated in the KS basis ?

  integer :: has_opermatlu
  ! Is the operator allocated in the local basis ?
!
  integer :: mbandc
  ! Total number of correlated bands

  integer :: natom
  ! Number of atoms

  integer :: nkpt
  ! Number of k-point in the IBZ.
!
  integer :: nspinor
  ! Number of spinors
!
  integer :: nsppol
  ! Number of spin polarizations

  integer :: paral
  ! =1 if the operator has been memory-parallelized over kpt, 0 otherwise

  integer :: shiftk
  ! Shift to get the physical kpt index (when the operator is memory-parallelized over kpt)

  !character(len=12) :: whichoper
  ! describe the type of operator computed (DFT, DMFT, KS..)

!  ! Polarisation
  type(matlu_type), allocatable :: matlu(:)
  ! Local projection on correlated orbitals

  complex(dp), allocatable :: ks(:,:,:,:)
  ! In the KS basis  (mbandc,mbandc,nkpt,nsppol)

  real(dp), ABI_CONTIGUOUS pointer :: wtk(:) => null()
  ! Weights for each kpt

 end type oper_type
!!***

!----------------------------------------------------------------------


CONTAINS  !========================================================================================
!!***

!!****f* m_oper/init_oper
!! NAME
!! init_oper
!!
!! FUNCTION
!!  Allocate variables used in type oper_type.
!!
!! INPUTS
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  nkpt = number of k-pts
!!  wtk = weights for each k-pt
!!  shiftk = shift for the kpt index
!!  opt_ksloc = 1: initialize in KS space only
!!            = 2: initialize in local space only
!!            = 3: initialize in both KS and local space
!!
!! OUTPUTS
!!  oper <type(oper_type)>= operator
!!
!! SOURCE

subroutine init_oper(paw_dmft,oper,nkpt,wtk,shiftk,opt_ksloc)

!Arguments ------------------------------------
 integer, optional, intent(in) :: nkpt,opt_ksloc,shiftk
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(oper_type), intent(inout) :: oper
 real(dp), target, optional :: wtk(paw_dmft%nkpt)
!Local variables ------------------------------------
 integer :: optksloc
!************************************************************************

 DBG_ENTER("COLL")

 optksloc = 3
 if (present(opt_ksloc)) optksloc = opt_ksloc

 !if(optksloc/=3) then
    ! FIXME: empty line!
 !endif

 oper%has_operks    = 0
 oper%has_opermatlu = 0
 oper%gpu_option    = ABI_GPU_DISABLED

! ===================
!  Integers
! ===================
 oper%mbandc  = paw_dmft%mbandc
 oper%natom   = paw_dmft%natom
 oper%nspinor = paw_dmft%nspinor
 oper%nsppol  = paw_dmft%nsppol
 oper%paral   = 0
 oper%shiftk  = 0
 oper%ndat    = 1

 if (present(shiftk)) oper%shiftk = shiftk

 oper%nkpt = paw_dmft%nkpt
 if (present(nkpt)) oper%nkpt = nkpt

 if (present(shiftk) .or. oper%nkpt /= paw_dmft%nkpt) oper%paral = 1

! allocate(oper%wtk(oper%nkpt))
 if (present(wtk)) then
   oper%wtk => wtk(:)
 else
   oper%wtk => paw_dmft%wtk(:)
 end if ! present(wtk)

! ===================
!  KS variables
! ===================
 if (optksloc == 1 .or. optksloc == 3) then

   ABI_MALLOC(oper%ks,(oper%mbandc,oper%mbandc,oper%nkpt,oper%nsppol))
   oper%has_operks  = 1
   oper%ks(:,:,:,:) = czero

 end if ! optksloc=1 or optksloc=3

! ===================
!  matlu variables
! ===================
 if (optksloc == 2 .or. optksloc == 3) then
   ABI_MALLOC(oper%matlu,(oper%natom))
   oper%has_opermatlu = 1
   call init_matlu(oper%natom,oper%nspinor,oper%nsppol,paw_dmft%lpawu(:),oper%matlu(:))
 end if ! optksloc=2 or optksloc=3

 DBG_EXIT("COLL")

end subroutine init_oper
!!***

!!****f* m_oper/init_oper_ndat
!! NAME
!! init_oper_ndat
!!
!! FUNCTION
!!  Allocate variables used in type oper_type.
!!
!! INPUTS
!!
!! OUTPUTS
!! oper  = operator of type oper_type
!!
!! SOURCE

subroutine init_oper_ndat(paw_dmft,oper,ndat,nkpt,wtk,shiftk,opt_ksloc,gpu_option)

 use m_matlu, only : init_matlu
 use m_paw_dmft, only : paw_dmft_type

!Arguments ------------------------------------
 integer, optional, intent(in) :: nkpt,opt_ksloc,shiftk,gpu_option
 integer, intent(in) :: ndat
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(oper_type), intent(inout) :: oper
 real(dp), target, optional :: wtk(paw_dmft%nkpt)
!Local variables ------------------------------------
 integer :: optksloc,ndat_,l_gpu_option
!************************************************************************

 DBG_ENTER("COLL")

 optksloc = 3
 if (present(opt_ksloc)) optksloc = opt_ksloc
 l_gpu_option=ABI_GPU_DISABLED; if(present(gpu_option)) l_gpu_option=gpu_option

 !if(optksloc/=3) then
    ! FIXME: empty line!
 !endif
 ndat_=ndat;

 oper%gpu_option    = l_gpu_option
 oper%has_operks    = 0
 oper%has_opermatlu = 0

! ===================
!  Integers
! ===================
 oper%mbandc  = paw_dmft%mbandc
 oper%natom   = paw_dmft%natom
 oper%nspinor = paw_dmft%nspinor
 oper%nsppol  = paw_dmft%nsppol
 oper%paral   = 0
 oper%shiftk  = 0
 oper%ndat    = ndat_

 if (present(shiftk)) oper%shiftk = shiftk

 oper%nkpt = paw_dmft%nkpt
 if (present(nkpt)) oper%nkpt = nkpt

 if (present(shiftk) .or. oper%nkpt /= paw_dmft%nkpt) oper%paral = 1

! allocate(oper%wtk(oper%nkpt))
 if (present(wtk)) then
   oper%wtk => wtk(:)
 else
   oper%wtk => paw_dmft%wtk(:)
 end if ! present(wtk)

! ===================
!  KS variables
! ===================
 if (optksloc == 1 .or. optksloc == 3) then

   ABI_MALLOC(oper%ks,(oper%mbandc,oper%mbandc*ndat_,oper%nkpt,oper%nsppol))
   oper%has_operks  = 1
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(alloc:oper%ks) IF(l_gpu_option==ABI_GPU_OPENMP)
#endif
   if(gpu_option==ABI_GPU_OPENMP) then
     call gpu_set_to_zero_complex(oper%ks, int(oper%nsppol,c_size_t)*ndat_*oper%mbandc*oper%mbandc*oper%nkpt)
   else
     oper%ks(:,:,:,:) = czero
   end if

 end if ! optksloc=1 or optksloc=3

! ===================
!  matlu variables
! ===================
 if (optksloc == 2 .or. optksloc == 3) then
   ABI_MALLOC(oper%matlu,(oper%natom))
   oper%has_opermatlu = 1
   call init_matlu(oper%natom,oper%nspinor,oper%nsppol*ndat_,paw_dmft%lpawu(:),oper%matlu(:),gpu_option=l_gpu_option)
 end if ! optksloc=2 or optksloc=3

 DBG_EXIT("COLL")

end subroutine init_oper_ndat
!!***

!!****f* m_oper/destroy_oper
!! NAME
!! destroy_oper
!!
!! FUNCTION
!!  deallocate oper
!!
!! INPUTS
!!  oper <type(oper_type)>= operator
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_oper(oper)

!Arguments ------------------------------------
 type(oper_type), intent(inout) :: oper
!Local variables-------------------------------
!! *********************************************************************

 DBG_ENTER("COLL")

 if (oper%has_opermatlu == 1) then
   call destroy_matlu(oper%matlu(:),oper%natom)
 !else
 !  message = " Operator is not defined to be used in destroy_oper"
 !  ABI_ERROR(message)
 end if ! has_opermatlu=1

 if (allocated(oper%matlu)) then
   ABI_FREE(oper%matlu)
   oper%has_opermatlu = 0
 end if

 if (allocated(oper%ks)) then
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET EXIT DATA MAP(delete:oper%ks) IF(oper%gpu_option==ABI_GPU_OPENMP)
#endif
   ABI_FREE(oper%ks)
   oper%has_operks = 0
 end if

 oper%wtk => null()
!  no deallocation for wtk: wtk is an explicit pointer

 DBG_EXIT("COLL")

end subroutine destroy_oper
!!***

!!****f* m_oper/copy_oper
!! NAME
!! copy_oper
!!
!! FUNCTION
!!  Copy oper1 into oper2
!!
!! INPUTS
!!  oper1 <type(oper_type)>= operator
!!
!! OUTPUT
!!  oper2 <type(oper_type)>= operator
!!
!! SOURCE

subroutine copy_oper(oper1,oper2)

!Arguments ------------------------------------
 type(oper_type), intent(in) :: oper1
 type(oper_type), intent(inout) :: oper2 !vz_i
!Local variables-------------------------------
! *********************************************************************

 DBG_ENTER("COLL")

 if (oper1%has_opermatlu == 1 .and. oper2%has_opermatlu == 1) then
   call copy_matlu(oper1%matlu(:),oper2%matlu(:),oper1%natom)
 end if

 if (oper1%has_operks == 1 .and. oper2%has_operks == 1) &
    & oper2%ks(:,:,:,:) = oper1%ks(:,:,:,:)

 DBG_EXIT("COLL")

end subroutine copy_oper
!!***

!!****f* m_oper/copy_oper_from_ndat
!! NAME
!! copy_oper_from_ndat
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine copy_oper_from_ndat(oper1,oper2,ndat,nw,proct,me_freq,copy_ks)

 use defs_basis
 use m_matlu, only : copy_matlu_from_ndat
 use m_errors

!Arguments ------------------------------------
!type
 integer,intent(in) :: nw,ndat,me_freq
 logical,intent(in) :: copy_ks
 integer,intent(in) :: proct(nw)
 type(oper_type),target,intent(in) :: oper1
 type(oper_type),intent(inout) :: oper2(nw) !vz_i

!oper variables-------------------------------
 integer ::  ikpt, isppol, idat, iw, iatom, mbandc
 complex(dp), ABI_CONTIGUOUS pointer :: mat(:,:,:)
! *********************************************************************
 DBG_ENTER("COLL")
 ABI_CHECK(oper1%ndat==ndat, "Bad value for ndat!")
 mbandc=oper1%mbandc
 if(oper1%has_opermatlu==1 .and. oper1%gpu_option==ABI_GPU_OPENMP) then
   do iatom=1,oper1%natom
     mat => oper1%matlu(iatom)%mat ! array of structs in OpenMP loosely supported
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET UPDATE FROM(mat)
#endif
   end do
 end if
 idat=1
 do iw=1,nw
   if (proct(iw) /= me_freq) cycle
   if(oper1%has_opermatlu==1.and.oper2(iw)%has_opermatlu==1)  then
     call copy_matlu_from_ndat(oper1%matlu,oper2(iw)%matlu,oper1%natom,ndat,idat)
     idat=idat+1
   endif
 enddo

 if(allocated(oper1%ks) .and. copy_ks) then
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET UPDATE FROM(oper1%ks) IF(oper1%gpu_option==ABI_GPU_OPENMP)
#endif
   idat=1
   do iw=1,nw
     if (proct(iw) /= me_freq) cycle
     do isppol=1,oper1%nsppol
       do ikpt=1,oper1%nkpt
         oper2(iw)%ks(:,:,ikpt,isppol)=oper1%ks(:,1+(idat-1)*mbandc:idat*mbandc,ikpt,isppol)
       enddo
     enddo
     idat=idat+1
   enddo
 endif

 DBG_EXIT("COLL")
end subroutine copy_oper_from_ndat
!!***

!!****f* m_oper/copy_oper_to_ndat
!! NAME
!! copy_oper_to_ndat
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine copy_oper_to_ndat(oper1,oper2,ndat,nw,proct,me_freq,copy_ks)

 use defs_basis
 use m_matlu, only : copy_matlu_to_ndat
 use m_errors

!Arguments ------------------------------------
!type
 integer,intent(in) :: nw,ndat,me_freq
 logical,intent(in) :: copy_ks
 integer,intent(in) :: proct(nw)
 type(oper_type),intent(in) :: oper1(nw)
 type(oper_type),target,intent(inout) :: oper2 !vz_i

!oper variables-------------------------------
 integer :: ikpt, isppol, idat, iw, iatom, mbandc
 complex(dp), ABI_CONTIGUOUS pointer :: mat(:,:,:)
! *********************************************************************
 DBG_ENTER("COLL")
 ABI_CHECK(oper2%ndat==ndat, "Bad value for ndat!")
 ABI_CHECK(oper2%mbandc==oper1(1)%mbandc, "Bad value for mbandc!")
 mbandc=oper2%mbandc
 idat=1
 do iw=1,nw
   if (proct(iw) /= me_freq) cycle
   if(oper1(iw)%has_opermatlu==1.and.oper2%has_opermatlu==1)  then
     call copy_matlu_to_ndat(oper1(iw)%matlu,oper2%matlu,oper2%natom,ndat,idat)
     idat=idat+1
   endif
 enddo
 if(oper2%has_opermatlu==1 .and. oper2%gpu_option==ABI_GPU_OPENMP) then
   do iatom=1,oper2%natom
     mat => oper2%matlu(iatom)%mat ! array of structs in OpenMP loosely supported
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET UPDATE TO(mat)
#endif
   end do
 end if

 if(allocated(oper2%ks) .and. copy_ks) then
   ABI_CHECK(size(oper2%ks,dim=2) == mbandc*ndat, "well?")
   ABI_CHECK(size(oper2%ks,dim=1) == mbandc, "uh?")
   idat=1
   do iw=1,nw
     if (proct(iw) /= me_freq) cycle
     do isppol=1,oper2%nsppol
       do ikpt=1,oper2%nkpt
         oper2%ks(:,1+(idat-1)*mbandc:idat*mbandc,ikpt,isppol)=oper1(iw)%ks(:,:,ikpt,isppol)
       enddo
     enddo
     idat=idat+1
   enddo
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET UPDATE TO(oper2%ks) IF(oper2%gpu_option==ABI_GPU_OPENMP)
#endif
 endif

 DBG_EXIT("COLL")
end subroutine copy_oper_to_ndat
!!***

!!****f* m_oper/print_oper
!! NAME
!! print_oper
!!
!! FUNCTION
!!
!! INPUTS
!! oper <type(oper_type)>= operator
!! option= < 5: write diagonal part of KS occupation matrix
!!         > 8: write all elements of KS occup. matrix.
!! paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!! prtopt= in local space: print option for print_matlu
!!         in KS space: only prints if abs(prtopt)>=3
!!                      print off-diagonal elements if abs(prtopt)>=4
!!
!! OUTPUT
!!
!! SOURCE

subroutine print_oper(oper,option,paw_dmft,prtopt)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(oper_type), intent(in) :: oper
 integer, intent(in) :: option,prtopt
!Local variables-------------------------------
 integer :: ib,ib1,iband1,iband2,ikpt,isppol,mbandc,nkpt,nkptr
 character(len=50000) :: message
 logical  :: ximag
 real(dp) :: maximag
! *********************************************************************

 DBG_ENTER("COLL")

 if (oper%has_opermatlu == 1) then
   write(message,'(2a)') ch10,'   = In the atomic basis'
   call wrtout(std_out,message,'COLL')
   call print_matlu(oper%matlu(:),oper%natom,prtopt)
 end if ! has_opermatlu=1

 if (oper%has_operks == 1) then
   write(message,'(2a)') ch10,'   = In the Kohn-Sham basis'
   call wrtout(std_out,message,'COLL')

!todo_ba complete print_out
   mbandc  = oper%mbandc
   iband1  = 1
   iband2  = mbandc
   maximag = zero
   nkpt    = oper%nkpt
   ximag   = .false.
!   do ib=1,oper%mbandc
!     if(-(paw_dmft%eigen_dft(1,1,ib)+paw_dmft%fermie).ge.0.3) iband1=ib
!     if( (paw_dmft%eigen_dft(1,1,ib)-paw_dmft%fermie).le.0.3) iband2=ib
!   enddo

   if (abs(prtopt) >= 3 .and. ((option < 5) .or. (option > 8))) then
!     write(message,'(x,a,a,i4,2x,a)') ch10,'  -KS states'
!     call wrtout(std_out,message,'COLL')
     nkptr = min(nkpt,4)
     do isppol=1,paw_dmft%nsppol
       write(message,'(a,3x,a,1x,i1)') ch10,"--isppol--",isppol
       call wrtout(std_out,message,'COLL')
       write(message,'(2a)') ch10,&
         & "   - (in the following only the values for the correlated bands and the first k-points are printed)"
       call wrtout(std_out,message,'COLL')
       do ikpt=1,nkptr
         write(message,'(2a,i4,2x,f14.5,a)') ch10,&
             & "   -k-pt--",ikpt,oper%wtk(ikpt),"(<-weight(k-pt))"
         call wrtout(std_out,message,'COLL')
         if (option < 5) then
           write(message,'(19x,a,6x,a)') "Eigenvalues","Occupations"
           call wrtout(std_out,message,'COLL')
         else if (abs(prtopt) >= 4 .or. option > 8) then
           write(message,'(a,10x,2000(i5,12x))') ch10,(paw_dmft%include_bands(ib),ib=iband1,iband2)
           call wrtout(std_out,message,'COLL')
         end if ! option
         do ib=1,mbandc
           if (option < 5) then
             if (abs(aimag(oper%ks(ib,ib,ikpt,isppol))) >= tol10) then
               write(message,'(a,i5,e14.5,3x,e14.5,3x,e21.14)') "   -iband--",paw_dmft%include_bands(ib),&
                 & paw_dmft%eigen_dft(ib,ikpt,isppol),oper%ks(ib,ib,ikpt,isppol)
             else
               write(message,'(a,i5,e14.5,3x,e14.5)') "   -iband--",paw_dmft%include_bands(ib),&
                 & paw_dmft%eigen_dft(ib,ikpt,isppol),dble(oper%ks(ib,ib,ikpt,isppol))
             end if ! imaginary part
             call wrtout(std_out,message,'COLL')
           end if ! option<5
           if (abs(prtopt) >= 4 .or. option > 8 .and. ib >= iband1 .and. ib <= iband2) then

             write(message,'(i5,1x,2000(2f7.3,3x))') paw_dmft%include_bands(ib),(dble(oper%ks(ib,ib1,ikpt,isppol)), &
                 & aimag(oper%ks(ib,ib1,ikpt,isppol)),ib1=iband1,iband2)
             call wrtout(std_out,message,'COLL')

!   to write imaginary part
!             write(message, '(1000(2f9.3,2x))') &
!&               (real(oper%ks(isppol,ikpt,ib,ib1)),imag(oper%ks(isppol,ikpt,ib,ib1)),ib1=iband1,iband2)
!             call wrtout(std_out,message,'COLL')
           end if ! prtopt>=20
           if (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7) then ! no sense to perform this check on off-diagonal elements
             if (abs(aimag(oper%ks(ib,ib,ikpt,isppol))) > max(tol10,maximag)) then
               ximag   = .true.
               maximag = aimag(oper%ks(ib,ib,ikpt,isppol))
             end if
           else
             do ib1=1,mbandc
               if (abs(aimag(oper%ks(ib1,ib,ikpt,isppol))) > max(tol10,maximag)) then
                 ximag   = .true.
                 maximag = aimag(oper%ks(ib1,ib,ikpt,isppol))
               end if
             end do ! ib1
           end if
         end do ! ib
       end do ! ikpt
     end do ! isppol
   else
    write(message,'(5x,a,i10,a)') '(not written)'
    call wrtout(std_out,message,'COLL')
   end if ! abs(prtopt)>=3 and (option<5 or option>8)
   if (ximag) then
     write(message,'(3a,e12.4,a)') "Occupations are imaginary !",ch10,&
        & "  Maximal value is ",maximag,ch10
     ABI_WARNING(message)
   end if ! ximag
 else if (abs(prtopt) >= 3 .and. ((option < 5) .or. (option > 8))) then
   write(message, '(2a)') ch10," Prb with options and has_operks in print_oper"
   call wrtout(std_out,message,'COLL')
 end if ! if oper%has_operks
! write(message, '(2a)') ch10," end print_oper"
!     call wrtout(std_out,message,'COLL')

 DBG_EXIT("COLL")

end subroutine print_oper
!!***

!!****f* m_oper/inverse_oper
!! NAME
!! inverse_oper
!!
!! FUNCTION
!!  Compute the inverse of the operator either in the KS space or in the
!!  correlated subspace.
!!
!! INPUTS
!!  oper <type(oper_type)>= operator
!!  option=1 do inversion in KS band space
!!        =2 do inversion in local space
!!        =3 do both
!!  procb(ikpt)=for kpt parallelization; gives the rank (in the kpt communicator) of the CPU handling each ikpt
!!  iproc=rank of the current process in the kpt communicator
!!
!! OUTPUT
!!  oper <type(oper_type)>= operator inverted
!!
!! SOURCE

subroutine inverse_oper(oper,option,procb,iproc,gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: option
 type(oper_type), target, intent(inout) :: oper
 integer, optional, intent(in) :: iproc,gpu_option
 integer, optional, intent(in) :: procb(oper%nkpt)
!Local variables-------------------------------
 integer :: ikpt,isppol,idat,paral,mbandc
 integer :: l_gpu_option
 !integer :: blk,iatom,ib
 complex(dp), ABI_CONTIGUOUS pointer :: ks(:,:,:,:)
#ifdef HAVE_OPENMP_OFFLOAD
 complex(dp), allocatable :: work(:,:)
 complex(dp), ABI_CONTIGUOUS pointer :: mat(:,:,:)
#endif
!todo_ba: prb with gwpc here: necessary for matcginv but should be dp
! *********************************************************************

 DBG_ENTER("COLL")
 ABI_NVTX_START_RANGE(NVTX_DMFT_INVERSE_OPER)

 l_gpu_option=ABI_GPU_DISABLED; if(present(gpu_option)) l_gpu_option=gpu_option
 paral = 0
 ks => oper%ks
 mbandc = oper%mbandc
 if (present(procb) .and. present(iproc) .and. oper%paral == 0) paral = 1

 !if (((option == 1 .or. option == 3) .and. (oper%has_operks == 0)) .or. &
 !  & ((option == 2 .or. option == 3) .and. (oper%has_opermatlu == 0))) then
 !  message = " Options are not coherent with definitions of this operator"
 !  ABI_ERROR(message)
 !end if

 if (option == 2 .or. option == 3) then
   call inverse_matlu(oper%matlu(:),oper%natom)
 end if

#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET ENTER DATA MAP(alloc:ks) IF(l_gpu_option==ABI_GPU_OPENMP .and. oper%gpu_option/=ABI_GPU_OPENMP)
 !$OMP TARGET UPDATE TO(ks) IF(l_gpu_option==ABI_GPU_OPENMP .and. oper%gpu_option/=ABI_GPU_OPENMP)
#endif
 if (option == 1 .or. option == 3) then
   do isppol=1,oper%nsppol
     do ikpt=1,oper%nkpt
       if (paral == 1) then
         if (procb(ikpt) /= iproc) cycle
       end if
       if(l_gpu_option==ABI_GPU_DISABLED) then
         do idat=1,oper%ndat
  !          write(std_out,*) "isppol,ikpt",isppol,ikpt,m
  !          write(std_out,*) "isppol,ikpt",matrix
           !call matcginv_dpc(matrix,oper%mbandc,oper%mbandc)
           call xginv(oper%ks(:,1+(idat-1)*oper%mbandc:idat*oper%mbandc,ikpt,isppol),oper%mbandc)
         end do ! idat
       else if(l_gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD

         ABI_MALLOC(work, (mbandc,mbandc*oper%ndat))
         !$OMP TARGET ENTER DATA MAP(alloc:work)
         !do idat=1,oper%ndat,32
           !blk=min(32,oper%ndat-idat+1)
           !$OMP TARGET DATA USE_DEVICE_ADDR(ks,work)
           !call gpu_xginv_strided(2,mbandc,ks(:,1+(idat-1)*mbandc:idat*mbandc,ikpt,isppol),mbandc,mbandc*mbandc,1)
           !call gpu_xginv_strided(2,mbandc,ks(:,1+(idat-1)*mbandc:(idat+blk-1)*mbandc,ikpt,isppol),mbandc,mbandc*mbandc,blk,work)
           call gpu_xginv_strided(2,mbandc,ks(:,:,ikpt,isppol),mbandc,mbandc*mbandc,oper%ndat,work)
           !$OMP END TARGET DATA
         !end do ! idat
         !$OMP TARGET EXIT DATA MAP(delete:work)
         ABI_FREE(work)

#endif
       end if
     end do ! ikpt
   end do ! isppol
 end if ! option
#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET EXIT DATA MAP(from:ks) IF(l_gpu_option==ABI_GPU_OPENMP .and. oper%gpu_option/=ABI_GPU_OPENMP)
#endif

 ABI_NVTX_END_RANGE()
 DBG_EXIT("COLL")

end subroutine inverse_oper
!!***

!!****f* m_oper/downfold_oper
!! NAME
!! downfold_oper
!!
!! FUNCTION
!!  Downfold an operator from KS space to local space.
!!
!! INPUTS
!!  oper <type(oper_type)>= operator
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  procb(ikpt)=for kpt parallelization; gives the rank (in the kpt communicator) of the CPU handling each ikpt
!!  iproc=rank of the current process in the kpt communicator
!!  option = 1 (default) : downfold an operator represented by a matrix in KS space
!!         = 2 : downfold the identity
!!         = 3 : downfold a diagonal KS operator
!!         = 4 : computes downfold(upfold)
!!  op_ks_diag = when option=3, you can provide the diagonal KS operator in this variable
!!               with the format mbandc*nkpt*nsppol, instead of storing it in the mband*mband matrix
!!               of oper%ks
!!
!! OUTPUT
!!
!! SOURCE

subroutine downfold_oper(oper,paw_dmft,procb,iproc,option,op_ks_diag,gpu_option)

!Arguments ------------------------------------
 type(oper_type),target,intent(inout) :: oper
 type(paw_dmft_type),target,intent(in) :: paw_dmft
 integer, optional, intent(in) :: iproc,option,gpu_option
 integer, optional, intent(in) :: procb(oper%nkpt)
 real(dp), optional, intent(in) :: op_ks_diag(oper%mbandc,oper%nkpt,oper%nsppol)
!oper variables-------------------------------
 integer :: iatom,ib,ik,ikpt,isppol,im,idat,lpawu,mbandc,ndim
 integer :: ndim_max,nspinor,ndat,opt,paral,shift
 integer :: l_gpu_option
 complex(dp) :: alpha
 complex(dp), ABI_CONTIGUOUS pointer :: ks(:,:,:,:),mat(:,:,:),chipsi(:,:,:,:,:)
 real(dp), ABI_CONTIGUOUS pointer :: wtk(:)
 character(len=500) :: message
 complex(dp), allocatable :: mat_temp(:,:,:),mat_temp2(:,:,:),mat_temp3(:,:)
! *********************************************************************

 DBG_ENTER("COLL")
 ABI_NVTX_START_RANGE(NVTX_DMFT_DOWNFOLD_OPER)

#ifndef HAVE_OPENMP_OFFLOAD
 ABI_UNUSED(alpha); ABI_UNUSED(im)
#endif

 if (oper%has_opermatlu == 0) then
   message = " Operator is not defined to be used in downfold_oper"
   ABI_ERROR(message)
 end if

 l_gpu_option = ABI_GPU_DISABLED; if(present(gpu_option)) l_gpu_option = gpu_option
 paral        = 0; if (present(procb) .and. present(iproc) .and. oper%paral == 0) paral = 1
 opt          = 1; if (present(option)) opt = option

 if(l_gpu_option==ABI_GPU_OPENMP) then
   ABI_CHECK(opt==1 .or. opt==3, "Incompatible codepath with OpenMP GPU")
 end if

 mbandc   = oper%mbandc
 nspinor  = oper%nspinor
 ndim_max = nspinor * (2*paw_dmft%maxlpawu+1)
 shift    = oper%shiftk
 ndat     = oper%ndat
 if(l_gpu_option==ABI_GPU_OPENMP) then
   ks => oper%ks
   wtk => oper%wtk
   chipsi => paw_dmft%chipsi
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(to:chipsi,wtk)
   !$OMP TARGET ENTER DATA MAP(to:ks) IF(oper%gpu_option/=ABI_GPU_OPENMP)
   if (present(op_ks_diag)) then
     !$OMP TARGET ENTER DATA MAP(to:op_ks_diag)
   end if
#endif
 end if

 do iatom=1,oper%natom
   lpawu = oper%matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   mat => oper%matlu(iatom)%mat
   ndim = nspinor * (2*lpawu+1)
   if(oper%gpu_option==ABI_GPU_DISABLED) then
     mat(:,:,:) = czero
   else if(oper%gpu_option==ABI_GPU_OPENMP) then
     call gpu_set_to_zero_complex(mat, int(oper%nsppol,c_size_t)*ndat*ndim*ndim)
   end if
   ABI_MALLOC(mat_temp,(ndim,mbandc,ndat))
   ABI_MALLOC(mat_temp2,(ndim,ndim,ndat))
   ABI_MALLOC(mat_temp3,(ndim,ndim))
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(alloc:mat_temp,mat_temp2,mat_temp3) IF(l_gpu_option==ABI_GPU_OPENMP)
   !$OMP TARGET ENTER DATA MAP(to:mat) IF(l_gpu_option==ABI_GPU_OPENMP .and. oper%gpu_option/=ABI_GPU_OPENMP)
#endif
   do isppol=1,oper%nsppol
     do ikpt=1,oper%nkpt ! index of kpt on the current CPU

       if (paral == 1) then
         if (procb(ikpt) /= iproc) cycle
       end if

       ik = ikpt + shift ! true kpt index (needed for chipsi)

       if (opt == 1 .or. opt == 3) then


         if (opt == 1) then

           if(l_gpu_option == ABI_GPU_DISABLED) then
             call abi_zgemm_2dd("n","n",ndim,mbandc*ndat,mbandc,cone,paw_dmft%chipsi(:,:,ik,isppol,iatom),&
             &    ndim_max,oper%ks(:,:,ikpt,isppol),mbandc,czero,mat_temp(:,:,:),ndim)
           else if(l_gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
             !$OMP TARGET DATA USE_DEVICE_ADDR(mat_temp,chipsi,ks)
             call abi_gpu_xgemm(2,"n","n",ndim,mbandc*ndat,mbandc,cone,c_loc(chipsi(:,:,ik,isppol,iatom)),&
             &    ndim_max,c_loc(ks(:,:,ikpt,isppol)),mbandc,czero,c_loc(mat_temp(:,:,:)),ndim)
             !$OMP END TARGET DATA
#endif
           end if

         else if (opt == 3) then

           if(l_gpu_option == ABI_GPU_DISABLED) then
             do idat=1,ndat
               do ib=1,mbandc
                 if (present(op_ks_diag)) then
                   mat_temp(:,ib,idat) = paw_dmft%chipsi(1:ndim,ib,ik,isppol,iatom) * op_ks_diag(ib,ikpt,isppol)
                 else
                   mat_temp(:,ib,idat) = paw_dmft%chipsi(1:ndim,ib,ik,isppol,iatom) * oper%ks(ib,ib+(idat-1)*mbandc,ikpt,isppol)
                 end if ! present(op_ks_diag)
               end do ! ib
             end do ! ndat
           else if(l_gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
             if (present(op_ks_diag)) then
               !$OMP TARGET TEAMS DISTRIBUTE MAP(to:chipsi,op_ks_diag,mat_temp) PRIVATE(idat)
               do idat=1,ndat
                 !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(ib,im)
                 do ib=1,mbandc
                   do im=1,ndim
                     mat_temp(im,ib,idat) = chipsi(im,ib,ik,isppol,iatom) * op_ks_diag(ib,ikpt,isppol)
                   end do
                 end do
               end do
             else
               !$OMP TARGET TEAMS DISTRIBUTE MAP(to:chipsi,ks,mat_temp) PRIVATE(idat)
               do idat=1,ndat
                 !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(ib,im)
                 do ib=1,mbandc
                   do im=1,ndim
                     mat_temp(im,ib,idat) = chipsi(im,ib,ik,isppol,iatom) * ks(ib,ib+(idat-1)*mbandc,ikpt,isppol)
                   end do
                 end do
               end do
             end if ! present(op_ks_diag)
#endif
           end if

         end if ! opt=1 or 3

         if(l_gpu_option == ABI_GPU_DISABLED) then
           do idat=1,ndat
             call abi_xgemm("n","c",ndim,ndim,mbandc,cone,mat_temp(:,:,idat),ndim,&
             &    paw_dmft%chipsi(:,:,ik,isppol,iatom),ndim_max,czero,mat_temp2(:,:,idat),ndim)
           end do ! ndat
         else if(l_gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
           !$OMP TARGET DATA USE_DEVICE_ADDR(mat_temp,chipsi,mat_temp2)
           call abi_gpu_xgemm_strided(2,'n','c',ndim,ndim,mbandc,cone,c_loc(mat_temp(:,:,:)),ndim,ndim*mbandc,&
           &    c_loc(chipsi(:,:,ik,isppol,iatom)),ndim_max,0,czero,c_loc(mat_temp2(:,:,:)),ndim,ndim*ndim,ndat)
           !$OMP END TARGET DATA
#endif
         end if

       else if (opt == 2) then

         call abi_xgemm("n","c",ndim,ndim,mbandc,cone,paw_dmft%chipsi(:,:,ik,isppol,iatom),&
                      & ndim_max,paw_dmft%chipsi(:,:,ik,isppol,iatom),ndim_max,czero,mat_temp2(:,:,1),ndim)

       else if (opt == 4) then

         call abi_xgemm("n","c",ndim,ndim,mbandc,cone,paw_dmft%chipsi(:,:,ik,isppol,iatom),&
                      & ndim_max,paw_dmft%chipsi(:,:,ik,isppol,iatom),ndim_max,czero,mat_temp3(:,:),ndim)

         call abi_xgemm("n","n",ndim,ndim,ndim,cone,mat_temp3(:,:),ndim,&
                      & mat_temp3(:,:),ndim,czero,mat_temp2(:,:,1),ndim)

       end if ! opt

       if(l_gpu_option == ABI_GPU_DISABLED) then
         do idat=1,ndat
           oper%matlu(iatom)%mat(:,:,idat+(isppol-1)*ndat) = &
           &    oper%matlu(iatom)%mat(:,:,idat+(isppol-1)*ndat) + mat_temp2(:,:,idat)*oper%wtk(ik)
         end do ! ndat
       else if(l_gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
         alpha = dcmplx(wtk(ik), 0.0_dp)
         !$OMP TARGET DATA USE_DEVICE_ADDR(mat,mat_temp2)
         call abi_gpu_xaxpy(2, ndim*ndim*ndat, alpha, &
         &    c_loc(mat_temp2), 1, c_loc(mat(:,:,1+(isppol-1)*ndat:isppol*ndat)), 1)
         !$OMP END TARGET DATA
#endif
       end if

     end do ! ikpt
   end do ! isppol
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET EXIT DATA MAP(from:mat) IF(l_gpu_option==ABI_GPU_OPENMP .and. oper%gpu_option/=ABI_GPU_OPENMP)
   !$OMP TARGET EXIT DATA MAP(delete:mat_temp,mat_temp2,mat_temp3) IF(l_gpu_option==ABI_GPU_OPENMP)
#endif
   ABI_FREE(mat_temp)
   ABI_FREE(mat_temp2)
   ABI_FREE(mat_temp3)
 end do ! iatom

#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET EXIT DATA MAP(delete:ks) IF(l_gpu_option==ABI_GPU_OPENMP .and. oper%gpu_option/=ABI_GPU_OPENMP)
 !$OMP TARGET EXIT DATA MAP(delete:chipsi,wtk) IF(l_gpu_option==ABI_GPU_OPENMP)
 if (present(op_ks_diag)) then
   !$OMP TARGET EXIT DATA MAP(delete:op_ks_diag) IF(l_gpu_option==ABI_GPU_OPENMP)
 end if
#endif
!do isppol=1,nsppol
 ! do ikpt=1,nkpt
 !  ikpt1=ikpt
 !  if(present(jkpt)) ikpt1=jkpt
 !  lvz=paral==0  !vz_d
 !  if(present(iproc)) lvz=lvz.or.(paral==1.and.(procb2(ikpt1)==iproc))  !vz_d
!!  if ((paral==1.and.(procb2(ikpt1)==iproc)).or.(paral==0)) then    !vz_d
 !  if(lvz) then !vz_d
 !  do ib=1,mbandc
 !   do ib1=1,mbandc
 !    do iatom=1,natom
 !     if(oper%matlu(iatom)%lpawu.ne.-1) then
 !     ndim=2*oper%matlu(iatom)%lpawu+1
 !      do im=1,ndim
 !       do im1=1,ndim
 !        do ispinor=1,nspinor
 !         do ispinor1=1,nspinor
 !           if (im1 == im .and. im1 == 1 .and. ib1 == ib .and. iatom == 1) then

  !          end if
  !          oper%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=     &
!&            oper%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)+    &
!&            paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)*        &
!&            conjg(paw_dmft%psichi(isppol,ikpt1,ib1,ispinor1,iatom,im1))* &
!false&            paw_dmft%psichi(isppol,ikpt1,ib1,ispinor1,iatom,im1)*
!&
!false&            conjg(paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im))* &
!&            oper%ks(isppol,ikpt,ib,ib1)*oper%wtk(ikpt)
! one  could suppress wtk here if present(jkpt)
! ks(ib,ib1)=ks(ib1,ib) -> ib and ib1 can be underchanged !
!          enddo ! ispinor1
!         enddo ! ispinor
!        enddo ! im1
!       enddo ! im
!      endif
!     enddo ! iatom
 !   enddo ! ib
 !  enddo ! ib
 !  endif
 ! enddo ! ikpt
 !enddo ! isppol

 DBG_EXIT("COLL")

 ABI_NVTX_END_RANGE()
end subroutine downfold_oper
!!***


!!****f* m_oper/upfold_oper
!! NAME
!! upfold_oper
!!
!! FUNCTION
!!  Upfold an operator from local space to KS space
!!
!! INPUTS
!!  oper <type(oper_type)>= operator
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  procb(ikpt)=for kpt parallelization; gives the rank (in the kpt communicator) of the CPU handling each ikpt
!!  iproc=rank of the current process in the kpt communicator
!!
!! OUTPUT
!!
!! SOURCE

subroutine upfold_oper(oper,paw_dmft,procb,iproc,gpu_option)

!Arguments ------------------------------------
 type(oper_type),target, intent(inout)  :: oper
 type(paw_dmft_type),target, intent(in) :: paw_dmft
 integer, optional, intent(in)   :: iproc,gpu_option
 integer, optional, intent(in)   :: procb(oper%nkpt)
!Local variables-------------------------------
 integer :: iatom,ik,ikpt,isppol,idat,lpawu,mbandc,l_gpu_option
 integer :: ndim,ndim_max,ndat,nspinor,paral,shift
 complex(dp), ABI_CONTIGUOUS pointer :: ks(:,:,:,:),mat(:,:,:),chipsi(:,:,:,:,:)
 complex(dp), allocatable :: mat_temp(:,:),mat_temp2(:,:)
! *********************************************************************

 l_gpu_option=ABI_GPU_DISABLED; if(present(gpu_option)) l_gpu_option=gpu_option

 DBG_ENTER("COLL")
 ABI_NVTX_START_RANGE(NVTX_DMFT_UPFOLD_OPER)

 !if ((oper%has_opermatlu == 0) .or. (oper%has_operks == 0)) then
 !  message = " Operator is not defined to be used in upfold_oper"
 !  ABI_ERROR(message)
 !end if

 mbandc   = paw_dmft%mbandc
 nspinor  = paw_dmft%nspinor
 ndim_max = nspinor * (2*paw_dmft%maxlpawu+1)
 paral    = 0
 shift    = oper%shiftk
 ndat     = oper%ndat

 if (present(procb) .and. present(iproc) .and. oper%paral == 0) paral = 1

 if(l_gpu_option==ABI_GPU_DISABLED) oper%ks(:,:,:,:) = czero

 ABI_MALLOC(mat_temp,(mbandc,ndim_max*ndat))
 ABI_MALLOC(mat_temp2,(mbandc,mbandc*ndat))
 if(l_gpu_option == ABI_GPU_OPENMP) then
   ks => oper%ks
   chipsi => paw_dmft%chipsi

#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(alloc:ks) IF(l_gpu_option==ABI_GPU_OPENMP .and. oper%gpu_option/=ABI_GPU_OPENMP)
   !$OMP TARGET ENTER DATA MAP(alloc:chipsi,mat_temp,mat_temp2) IF(l_gpu_option==ABI_GPU_OPENMP)
   !$OMP TARGET UPDATE TO(chipsi) IF(l_gpu_option==ABI_GPU_OPENMP)
   call gpu_set_to_zero_complex(ks, int(oper%nsppol,c_size_t)*ndat*mbandc*mbandc*oper%nkpt)
#endif
 end if

 do iatom=1,oper%natom
   lpawu = oper%matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim = (2*lpawu+1) * nspinor
   mat => oper%matlu(iatom)%mat
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(to:mat) IF(l_gpu_option==ABI_GPU_OPENMP .and. oper%gpu_option/=ABI_GPU_OPENMP)
#endif
   do isppol=1,oper%nsppol
     do ikpt=1,oper%nkpt ! index of kpt on the current CPU

       if (paral == 1) then
         if (procb(ikpt) /= iproc) cycle
       end if

       ik = ikpt + shift ! true kpt index (needed for chipsi)

       if(l_gpu_option == ABI_GPU_DISABLED) then

         call abi_zgemm_2dd("c","n",mbandc,ndat*ndim,ndim,cone,paw_dmft%chipsi(:,:,ik,isppol,iatom),&
                      & ndim_max,oper%matlu(iatom)%mat(:,:,(isppol-1)*ndat+1:isppol*ndat),ndim,czero,mat_temp(:,:),mbandc)

         do idat=1,ndat

           call abi_xgemm("n","n",mbandc,mbandc,ndim,cone,mat_temp(:,1+(idat-1)*ndim:idat*ndim),mbandc,&
                        & paw_dmft%chipsi(:,:,ik,isppol,iatom),ndim_max,czero,mat_temp2(:,1+(idat-1)*mbandc:idat*mbandc),mbandc)

         end do ! idat

         !oper%ks(:,:,ikpt,isppol) = oper%ks(:,:,ikpt,isppol) + mat_temp2(:,:)
         call zaxpy(mbandc*mbandc*ndat, cone, mat_temp2, 1, oper%ks(:,:,ikpt,isppol), 1)

       else if(l_gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET DATA USE_DEVICE_ADDR(mat_temp,chipsi,mat)
         call abi_gpu_xgemm(2,"c","n",mbandc,ndat*ndim,ndim,cone,c_loc(chipsi(:,:,ik,isppol,iatom)),&
         &    ndim_max,c_loc(mat(:,:,(isppol-1)*ndat+1:isppol*ndat)),ndim,czero,c_loc(mat_temp(:,:)),mbandc)
         !$OMP END TARGET DATA

         !$OMP TARGET DATA USE_DEVICE_ADDR(mat_temp,chipsi,mat_temp2)
         call abi_gpu_xgemm_strided(2,'n','n',mbandc,mbandc,ndim,cone,c_loc(mat_temp(:,:)),mbandc,ndim*mbandc,&
         &    c_loc(chipsi(:,:,ik,isppol,iatom)),ndim_max,0,czero,c_loc(mat_temp2(:,:)),mbandc,mbandc*mbandc,ndat)
         !$OMP END TARGET DATA
         !$OMP TARGET DATA USE_DEVICE_ADDR(ks,mat_temp2)
         call abi_gpu_xaxpy(1, 2*mbandc*mbandc*ndat, cone, &
         &    c_loc(mat_temp2), 1, c_loc(ks(:,:,ikpt,isppol)), 1)
         !$OMP END TARGET DATA
#endif
       end if

     end do ! ikpt
   end do ! isppol
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET EXIT DATA MAP(delete:mat) IF(l_gpu_option==ABI_GPU_OPENMP .and. oper%gpu_option/=ABI_GPU_OPENMP)
#endif
 end do ! iatom

#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET UPDATE FROM(ks) IF(l_gpu_option==ABI_GPU_OPENMP .and. oper%gpu_option/=ABI_GPU_OPENMP)
 !$OMP TARGET EXIT DATA MAP(delete:ks) IF(l_gpu_option==ABI_GPU_OPENMP .and. oper%gpu_option/=ABI_GPU_OPENMP)
 !$OMP TARGET EXIT DATA MAP(delete:chipsi,mat_temp,mat_temp2) IF(l_gpu_option==ABI_GPU_OPENMP)
#endif
 ABI_FREE(mat_temp)
 ABI_FREE(mat_temp2)

 ABI_NVTX_END_RANGE()
 DBG_EXIT("COLL")

end subroutine upfold_oper
!!***

!!****f* m_oper/identity_oper
!! NAME
!! identity_oper
!!
!! FUNCTION
!!  Construct the identity operator
!!
!! INPUTS
!!  oper <type(oper_type)>= operator
!!  option = 1: in KS space
!!         = 2: in local space
!!         = 3: both
!!
!! OUTPUT
!!
!! SOURCE

subroutine identity_oper(oper,option)

!Arguments ------------------------------------
 integer, intent(in) :: option
 type(oper_type), intent(inout) :: oper
!Local variables-------------------------------
 integer :: ib,natom
 character(len=500) :: message
! *********************************************************************

 DBG_ENTER("COLL")

 if (((option == 1 .or. option == 3) .and. (oper%has_operks == 0)) .or. &
   & ((option == 2 .or. option == 3) .and. (oper%has_opermatlu == 0))) then
   message = " Options in identity_oper are not coherent with definitions of this operator"
   ABI_ERROR(message)
 end if

 if (option == 1 .or. option == 3) then

   oper%ks(:,:,:,:) = czero
   do ib=1,oper%mbandc
     oper%ks(ib,ib,:,:) = cone
   end do ! ib

 end if ! option=1 or 3

 if (option == 2 .or. option == 3) then
   natom = oper%natom
   call zero_matlu(oper%matlu(:),natom)
   call identity_matlu(oper%matlu(:),natom)
 end if ! option=2 or 3

 DBG_EXIT("COLL")

end subroutine identity_oper
!!***

!!****f* m_oper/diff_oper
!! NAME
!! diff_oper
!!
!! FUNCTION
!! Compute a norm of the differences between two occupations matrices.
!!
!! INPUTS
!!  char1 = character describing occup1
!!  char2 = character describing occup2
!!  occup1 <type(oper_type)>= occupations
!!  occup2 <type(oper_type)>= occupations
!!  option : option for printing (if 1 assume data are related to lda only)
!!  toldiff : tolerance for the difference
!!
!! OUTPUT
!!
!! SOURCE

subroutine diff_oper(char1,char2,occup1,occup2,option,toldiff)

!Arguments ------------------------------------
 type(oper_type), intent(in) :: occup1,occup2
 integer, intent(in) :: option
 real(dp), intent(in) :: toldiff
 character(len=*), intent(in) :: char1,char2
!Local variables-------------------------------
 character(len=500) :: message
! *********************************************************************

 DBG_ENTER("COLL")

 if (occup1%has_opermatlu == 0 .or. occup2%has_opermatlu == 0) then
   message = " Operators are not defined to be used in diff_oper"
   ABI_ERROR(message)
 end if

 if (occup1%nkpt /= occup2%nkpt) then
   write(message,'(a,2x,2i9)')' Operators are not equal',occup1%nkpt,occup2%nkpt
   ABI_ERROR(message)
 end if

 call diff_matlu(char1,char2,occup1%matlu(:),occup2%matlu(:),occup1%natom,option,toldiff)
! if(option==1) then
!  toldiff=tol4
!  if( matludiff < toldiff ) then
!   write(message,'(6a,e12.4,a,e12.4)') ch10,&
!&   '   Differences between ',trim(char1),' and ',trim(char2),' is small enough:',&
!&   matludiff,'is lower than',toldiff
!   call wrtout(std_out,message,'COLL')
!  else
!   write(message,'(6a,e12.4,a,e12.4)') ch10,&
!&   '   Error: Differences between ',trim(char1),' and ',trim(char2),' is too large:',&
!&   matludiff,'is large than',toldiff
!   call wrtout(std_out,message,'COLL')
!   call abi_abort('COLL')
!  endif
! endif
! call abi_abort('COLL')

 DBG_EXIT("COLL")

end subroutine diff_oper
!!***

!!****f* m_oper/trace_oper
!! NAME
!! trace_oper
!!
!! FUNCTION
!!  Computes the trace of an operator
!!
!! INPUTS
!!  oper <type(oper_type)>= operator
!!  opt_ksloc = 1: trace in KS space
!!            = 2: trace in local space
!!            = 3: both
!!
!! OUTPUT
!!  trace_ks  :: trace in KS space
!!  trace_loc :: trace in local space
!!  trace_ks_cmplx :: complex trace in KS space
!!
!! SOURCE

subroutine trace_oper(oper,trace_ks,trace_loc,opt_ksloc,trace_ks_cmplx)

!Arguments ------------------------------------
 type(oper_type), intent(in) :: oper
 real(dp), intent(out) :: trace_ks  !vz_i
 real(dp), intent(inout) :: trace_loc(oper%nsppol+1,oper%natom) !vz_i
 integer, intent(in) :: opt_ksloc
 complex(dp), optional, intent(out) :: trace_ks_cmplx
!Local variables-------------------------------
 integer :: ib,ikpt,isppol
 complex(dp) :: trace
 character(len=500) :: message
! *********************************************************************

 DBG_ENTER("COLL")

 if (((opt_ksloc == 1 .or. opt_ksloc == 3) .and. (oper%has_operks == 0)) .or. &
   & ((opt_ksloc == 2 .or. opt_ksloc == 3) .and. (oper%has_opermatlu == 0))) then
   message = " Options in trace_oper are not coherent with definitions of this operator"
   ABI_ERROR(message)
 end if

 if (opt_ksloc == 1 .or. opt_ksloc == 3) then
   trace = czero
   !temp1=zero
   do isppol=1,oper%nsppol
     do ikpt=1,oper%nkpt
       do ib=1,oper%mbandc
         trace = trace + oper%ks(ib,ib,ikpt,isppol)*oper%wtk(ikpt+oper%shiftk)
         !temp1=temp1+oper%wtk(ikpt)
       end do ! ib
     end do ! ikpt
   end do ! isppol
   if (oper%nsppol == 1 .and. oper%nspinor == 1) trace = two * trace
   if (present(trace_ks_cmplx)) trace_ks_cmplx = trace
   trace_ks = dble(trace)
!   write(std_out,*) "temp1",temp1
 end if ! opt_ksloc

 if (opt_ksloc == 2 .or. opt_ksloc == 3) then
   call trace_matlu(oper%matlu(:),oper%natom,trace_loc=trace_loc(:,:))
 end if

 DBG_EXIT("COLL")

end subroutine trace_oper
!!***

!!****f* m_oper/prod_oper
!! NAME
!! prod_oper
!!
!! FUNCTION
!!  Computes the matrix product of oper1 and oper2
!!
!! INPUTS
!!  oper1,oper2 <type(oper_type)>= operator
!!  opt_ksloc = 1 : in KS space
!!            = 2 : in local space
!!  opt_diag = 1 if oper1 and oper2 are diagonal in KS space, 0 otherwise (default)
!!
!! OUTPUT
!!  oper3 <type(oper_type)>= matrix product of oper1 and oper2
!!
!! SOURCE

subroutine prod_oper(oper1,oper2,oper3,opt_ksloc,opt_diag)

!Arguments ------------------------------------
 type(oper_type), intent(in) :: oper1,oper2
 type(oper_type), intent(inout) :: oper3
 integer, intent(in) :: opt_ksloc
 integer, optional, intent(in) :: opt_diag
!Local variables-------------------------------
 integer :: ib,ikpt,isppol,mbandc
 logical :: diag
! *********************************************************************

 DBG_ENTER("COLL")

 if (opt_ksloc == 2 .and. oper1%has_opermatlu == 1 .and. &
   & oper2%has_opermatlu == 1 .and. oper3%has_opermatlu == 1) then
   call prod_matlu(oper1%matlu(:),oper2%matlu(:),oper3%matlu(:),oper1%natom)
 end if

 if (opt_ksloc == 1 .and. oper1%has_operks == 1 .and. &
    & oper2%has_operks == 1 .and. oper3%has_operks == 1) then
   mbandc = oper1%mbandc
   diag   = .false.
   if (present(opt_diag)) then
     if (opt_diag == 1) diag = .true.
   end if
   if (diag) then
     do ib=1,mbandc
       oper3%ks(ib,ib,:,:) = oper1%ks(ib,ib,:,:) * oper2%ks(ib,ib,:,:)
     end do ! ib
   else
     do isppol=1,oper1%nsppol
       do ikpt=1,oper1%nkpt
         call abi_xgemm("n","n",mbandc,mbandc,mbandc,cone,oper1%ks(:,:,ikpt,isppol),mbandc,&
                      & oper2%ks(:,:,ikpt,isppol),mbandc,czero,oper3%ks(:,:,ikpt,isppol),mbandc)
       end do ! ikpt
     end do ! isppol
   end if ! diag
 end if ! opt_ksloc=1

 DBG_EXIT("COLL")

end subroutine prod_oper
!!***

!!****f* m_oper/trace_prod_oper
!! NAME
!! trace_prod_oper
!!
!! FUNCTION
!!  Computes Tr(oper1*oper2) in KS space
!!
!! INPUTS
!!  oper1,oper2 <type(oper_type)>= operator
!!
!! OUTPUT
!!  trace = Tr(op1*op2)
!!
!! SOURCE

subroutine trace_prod_oper(oper1,oper2,trace)

!Arguments ------------------------------------
 type(oper_type), intent(in) :: oper1,oper2
 complex(dp), intent(out) :: trace
!Local variables-------------------------------
 integer :: ikpt,isppol
 character(len=500) :: message
! *********************************************************************

 if (oper1%shiftk /= oper2%shiftk) then
   message = "Inconsistency in trace_prod_oper: oper1%shiftk should be equal to oper2%shiftk"
   ABI_ERROR(message)
 end if

 trace = czero

 do isppol=1,oper1%nsppol
   do ikpt=1,oper1%nkpt
     trace = trace + sum(oper1%ks(:,:,ikpt,isppol)*transpose(oper2%ks(:,:,ikpt,isppol)))*&
       & oper1%wtk(ikpt+oper1%shiftk)
   end do ! ikpt
 end do ! isppol

 if (oper1%nsppol == 1 .and. oper1%nspinor == 1) trace = trace * two

end subroutine trace_prod_oper
!!***

!!****f* m_oper/gather_oper
!! NAME
!! gather_oper
!!
!! FUNCTION
!!  Gather the contributions from all CPUs, for a frequency-dependent
!!  operator, and for both levels of parallelization (kpt and then frequency,
!!  and frequency only).
!!
!! INPUTS
!!  oper <type(oper_type)>= operator for each frequency
!!  distrib <type(mpi_distrib_dmft_type)> = mpi related data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  opt_ksloc = 1 : gather the KS operator on the kpt and frequency communicator
!!              2 : gather the local operator (the exact behavior can be defined via opt_commkpt)
!!  master = if present, only gather on the master node
!!  opt_diag = 1 if the operator is diagonal in KS space, 0 (default) otherwise
!!  opt_commkpt (only meaningful for the local quantity)
!!              = 0 (default) : frequency-only parallelization
!!                              -> xmpi_allgatherv on the whole communicator
!!              = 1 : kpt and then frequency parallelization (CAREFUL: here, frequencies
!!                     are not distributed in the same way as the frequency-only parallelization scheme)
!!                    -> xmpi_sum on the kpt-communicator, and then xmpi_allgatherv
!!                      on the frequency communicator (useful after a downfold for instance)
!!
!! OUTPUT
!!
!! SOURCE

subroutine gather_oper(oper,distrib,paw_dmft,opt_ksloc,master,opt_diag,opt_commkpt)

!Arguments ------------------------------------
 type(mpi_distrib_dmft_type), target, intent(in) :: distrib
 type(oper_type), intent(inout) :: oper(distrib%nw)
 type(paw_dmft_type) :: paw_dmft
 integer, intent(in) :: opt_ksloc
 integer, optional, intent(in) :: master,opt_commkpt,opt_diag
!Local variables-------------------------------
 integer :: comm,iatom,ib1,ibuf,ierr,ifreq,ikpt,im1
 integer :: irank,irank1,irank2,isppol,lpawu,mbandc,myproc
 integer :: myproc2,natom,ndim,nkpt,nproc,nproc_freq,nproc_kpt
 integer :: nproc2,nspinor,nsppol,nw,optcommkpt,siz_buf
 logical :: diag
 integer, allocatable :: displs(:),recvcounts(:)
 complex(dp), allocatable :: buffer(:),buffer_tot(:)
! *********************************************************************

 comm    = paw_dmft%spacecomm
 mbandc  = paw_dmft%mbandc
 myproc  = paw_dmft%myproc
 natom   = paw_dmft%natom
 nkpt    = paw_dmft%nkpt
 nproc   = paw_dmft%nproc
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 nw      = distrib%nw

 nproc_kpt  = min(nkpt,nproc)

 optcommkpt = 0
 if (present(opt_commkpt)) optcommkpt = opt_commkpt

 if (opt_ksloc == 1) then

   nproc_freq = max(1,nproc/nkpt)

   ABI_MALLOC(recvcounts,(nproc))
   ABI_MALLOC(displs,(nproc))

   diag = .false.
   if (present(opt_diag)) then
     if (opt_diag == 1) diag = .true.
   end if

   irank2 = 1
   do irank=0,nproc_kpt-1
     do irank1=0,nproc_freq-1
       recvcounts(irank2) = distrib%nkpt_mem(irank+1) * distrib%nw_mem_kptparal(irank1+1)
       irank2 = irank2 + 1
     end do ! irank1
   end do ! irank
   if (nproc > nproc_freq*nproc_kpt) recvcounts(nproc_freq*nproc_kpt+1:nproc) = 0

   recvcounts(:) = recvcounts(:) * merge(mbandc,mbandc**2,diag)
   displs(1) = 0
   do irank=2,nproc
     displs(irank) = displs(irank-1) + recvcounts(irank-1)
   end do ! irank

   ABI_MALLOC(buffer,(recvcounts(myproc+1)))
   ABI_MALLOC(buffer_tot,(displs(nproc)+recvcounts(nproc)))

   do isppol=1,nsppol

     ibuf = 0
     do ikpt=1,nkpt
       if (distrib%procb(ikpt) /= distrib%me_kpt) cycle
       do ifreq=1,nw
         if (distrib%proct(ifreq) /= distrib%me_freq) cycle
         do ib1=1,mbandc
           if (diag) then
             ibuf = ibuf + 1
             buffer(ibuf) = oper(ifreq)%ks(ib1,ib1,ikpt,isppol)
           else
             buffer(ibuf+1:ibuf+mbandc) = oper(ifreq)%ks(:,ib1,ikpt,isppol)
             ibuf = ibuf + mbandc
           end if ! diag
         end do ! ib1
       end do ! ifreq
     end do ! ikpt

     if (present(master)) then
       call xmpi_gatherv(buffer(:),recvcounts(myproc+1),buffer_tot(:),recvcounts(:),displs(:),master,comm,ierr)
     else
       call xmpi_allgatherv(buffer(:),recvcounts(myproc+1),buffer_tot(:),recvcounts(:),displs(:),comm,ierr)
     end if  ! present(master)

     ibuf = 0
     do ikpt=1,nkpt
       do ifreq=1,nw
         do ib1=1,mbandc
           if (diag) then
             ibuf = ibuf + 1
             oper(ifreq)%ks(ib1,ib1,ikpt,isppol) = buffer_tot(ibuf)
           else
             oper(ifreq)%ks(:,ib1,ikpt,isppol) = buffer_tot(ibuf+1:ibuf+mbandc)
             ibuf = ibuf + mbandc
           end if ! diag
         end do ! ib1
       end do ! ifreq
     end do ! ikpt

   end do ! isppol

 else if (opt_ksloc == 2) then

   nproc_freq = nproc / nkpt
   myproc2 = merge(distrib%me_freq,myproc,optcommkpt==1)
   nproc2  = merge(nproc_freq+1,nproc,optcommkpt==1)

   ABI_MALLOC(recvcounts,(nproc2))
   ABI_MALLOC(displs,(nproc2))

   siz_buf = 0

   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     siz_buf = siz_buf + (2*lpawu+1)**2
   end do ! iatom

   siz_buf = siz_buf * (nspinor**2) * nsppol
   if (optcommkpt == 1) then
     recvcounts(:) = siz_buf * distrib%nw_mem_kptparal(:)
   else
     recvcounts(:) = siz_buf * distrib%nw_mem(:)
   end if
   displs(1) = 0
   do irank=2,nproc2
     displs(irank) = displs(irank-1) + recvcounts(irank-1)
   end do ! irank

   if (optcommkpt == 1 .and. recvcounts(myproc2+1) == 0) then
     siz_buf = siz_buf * merge(nw,distrib%nw_mem_kptparal(mod(paw_dmft%myproc,nproc_freq)+1),nproc_freq<=1)
   else
     siz_buf = recvcounts(myproc2+1)
   end if

   ABI_MALLOC(buffer,(siz_buf))
   ABI_MALLOC(buffer_tot,(recvcounts(nproc2)+displs(nproc2)))

   buffer(:) = czero

   ibuf = 0
   do ifreq=1,nw
     if (optcommkpt == 1) then
       if (distrib%proct(ifreq) /= distrib%me_freq) cycle
     else if (optcommkpt == 0) then
       if (distrib%procf(ifreq) /= myproc) cycle
     end if ! optcommkpt
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       ndim = (2*lpawu+1) * nspinor
       do isppol=1,nsppol
         do im1=1,ndim
           buffer(ibuf+1:ibuf+ndim) = oper(ifreq)%matlu(iatom)%mat(:,im1,isppol)
           ibuf = ibuf + ndim
         end do ! im1
       end do ! isppol
     end do ! iatom
   end do ! ifreq

   if (optcommkpt == 1) comm = distrib%comm_freq
   if (present(master)) then
     if (optcommkpt == 1) then
       call xmpi_sum_master(buffer(:),master,distrib%comm_kpt,ierr)
     end if
     call xmpi_gatherv(buffer(:),recvcounts(myproc2+1),buffer_tot(:),recvcounts(:),displs(:),master,comm,ierr)
   else
     if (optcommkpt == 1) then
       call xmpi_sum(buffer(:),distrib%comm_kpt,ierr)
     end if
     call xmpi_allgatherv(buffer(:),recvcounts(myproc2+1),buffer_tot(:),recvcounts(:),displs(:),comm,ierr)
   end if  ! present(master)

   ibuf = 0
   do ifreq=1,nw
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       ndim = (2*lpawu+1) * nspinor
       do isppol=1,nsppol
         do im1=1,ndim
           oper(ifreq)%matlu(iatom)%mat(:,im1,isppol) = buffer_tot(ibuf+1:ibuf+ndim)
           ibuf = ibuf + ndim
         end do ! im1
       end do ! isppol
     end do ! iatom
   end do ! ifreq

 end if ! opt_ksloc

 ABI_FREE(recvcounts)
 ABI_FREE(displs)
 ABI_FREE(buffer)
 ABI_FREE(buffer_tot)

end subroutine gather_oper
!!***

!!****f* m_oper/gather_oper_ks
!! NAME
!! gather_oper_ks
!!
!! FUNCTION
!!  For a single KS operator, performs a xmpi_sum on the frequency communicator,
!!  and a xmpi_allgatherv on the kpt-communicator
!!
!! INPUTS
!!  oper <type(oper_type)>= operator
!!  distrib <type(mpi_distrib_dmft_type)> = mpi related data
!!  opt_diag = 1 if the operator is diagonal in KS space, 0 (default) otherwise
!!
!! OUTPUT
!!
!! SOURCE

subroutine gather_oper_ks(oper,distrib,paw_dmft,opt_diag)

!Arguments ------------------------------------
 type(oper_type), intent(inout) :: oper
 type(mpi_distrib_dmft_type), intent(in) :: distrib
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, optional, intent(in) :: opt_diag
!Local variables-------------------------------
 integer :: ib1,ibuf,ierr,ikpt,irank,isppol,mbandc
 integer :: me_kpt,nkpt,nproc,nproc_freq,nsppol,siz_buf
 logical :: diag
 integer, allocatable :: displs(:),recvcounts(:)
 complex(dp), allocatable :: buffer(:),buffer_tot(:)
! *********************************************************************

 mbandc = paw_dmft%mbandc
 me_kpt = distrib%me_kpt
 nkpt   = paw_dmft%nkpt
 nproc  = paw_dmft%nproc
 nsppol = paw_dmft%nsppol

 nproc_freq = nproc / nkpt

 diag = .false.
 if (present(opt_diag)) then
   if (opt_diag == 1) diag = .true.
 end if

 ABI_MALLOC(recvcounts,(nproc))
 ABI_MALLOC(displs,(nproc))

 recvcounts(:) = distrib%nkpt_mem(:) * merge(mbandc,mbandc**2,diag)

 displs(1) = 0
 do irank=2,nproc
   displs(irank) = displs(irank-1) + recvcounts(irank-1)
 end do ! irank

 siz_buf = recvcounts(me_kpt+1)
 if (siz_buf == 0) then
   siz_buf = mbandc
   if (nproc_freq <= 1) siz_buf = siz_buf * distrib%nkpt_mem(mod(paw_dmft%myproc,min(nkpt,nproc))+1)
   if (.not. diag) siz_buf = siz_buf * mbandc
 end if
 ABI_MALLOC(buffer,(siz_buf))
 ABI_MALLOC(buffer_tot,(recvcounts(nproc)+displs(nproc)))

 do isppol=1,nsppol

   buffer(:) = czero

   ibuf = 0
   do ikpt=1,nkpt
     if (distrib%procb(ikpt) /= me_kpt) cycle
     do ib1=1,mbandc
       if (diag) then
         ibuf = ibuf + 1
         buffer(ibuf) = oper%ks(ib1,ib1,ikpt,isppol)
       else
         buffer(ibuf+1:ibuf+mbandc) = oper%ks(:,ib1,ikpt,isppol)
         ibuf = ibuf + mbandc
       end if ! diag
     end do ! ib1
   end do ! ikpt

   call xmpi_sum(buffer(:),distrib%comm_freq,ierr)

   call xmpi_allgatherv(buffer(:),recvcounts(me_kpt+1),&
      & buffer_tot(:),recvcounts(:),displs(:),distrib%comm_kpt,ierr)

   ibuf = 0
   do ikpt=1,nkpt
     do ib1=1,mbandc
       if (diag) then
         ibuf = ibuf + 1
         oper%ks(ib1,ib1,ikpt,isppol) = buffer_tot(ibuf)
       else
         oper%ks(:,ib1,ikpt,isppol) = buffer_tot(ibuf+1:ibuf+mbandc)
         ibuf = ibuf + mbandc
       end if ! diag
     end do ! ib
   end do ! ikpt

 end do ! isppol

 ABI_FREE(buffer)
 ABI_FREE(buffer_tot)
 ABI_FREE(recvcounts)
 ABI_FREE(displs)

end subroutine gather_oper_ks
!!***

END MODULE m_oper
!!***
