!!****m* ABINIT/m_matlu
!! NAME
!!  m_matlu
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
!! NOTES
!!  subroutines in this module must never call
!!   a subroutine of m_oper, m_green, m_self
!!   in order to avoid circular dependancies
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

MODULE m_matlu

 use defs_basis
 use m_errors
 use m_abicore
 use, intrinsic :: iso_c_binding, only: c_size_t, c_loc
 use m_abi_linalg

#ifdef HAVE_GPU
 use m_gpu_toolbox
#endif

 use m_abi_linalg, only : abi_xgemm
 use m_fstrings, only : int2char4
 use m_hide_lapack, only : xginv
 use m_io_tools, only : flush_unit
 use m_matrix, only : blockdiago_fordsyev
 use m_paw_dmft, only : paw_dmft_type
 use m_xmpi, only : xmpi_bcast,xmpi_sum

 implicit none

 private

 public :: init_matlu
 public :: inverse_matlu
 public :: destroy_matlu
 public :: diff_matlu
 public :: add_matlu
 public :: print_matlu
 public :: sym_matlu
 public :: copy_matlu
 public :: copy_matlu_from_ndat
 public :: copy_matlu_to_ndat
 public :: gather_matlu
 public :: zero_matlu
 public :: trace_matlu
 public :: diag_matlu
 public :: rotate_matlu
 public :: shift_matlu
 public :: checkdiag_matlu
 public :: checkreal_matlu
 public :: prod_matlu
 public :: conjg_matlu
 public :: ln_matlu
 public :: slm2ylm_matlu
 public :: fac_matlu
 public :: printplot_matlu
 public :: identity_matlu
 public :: magmomforb_matlu
 public :: magmomfspin_matlu
 public :: magmomfzeeman_matlu
 public :: chi_matlu
 public :: trace_prod_matlu
 public :: xmpi_matlu
 public :: symmetrize_matlu
 public :: ylm2jmj_matlu
 public :: magnfield_matlu
 public :: magmomjmj_matlu
!!***

!!****t* m_matlu/matlu_type
!! NAME
!!  matlu_type
!!
!! FUNCTION
!!  This structured datatype contains a matrix for the correlated subspace
!!
!! SOURCE

 type, public :: matlu_type ! for each atom

  integer :: lpawu
  ! Value of the angular momentum for each correlated electrons

!  integer :: natom
   ! number of atoms (given for each atom, not useful..could be changed)
!
!  integer :: mband
!  ! Number of bands
!
!  integer :: mbandc
!  ! Total number of bands in the Kohn-Sham Basis for PAW+DMFT
!
!  integer :: natpawu         ! Number of correlated atoms
!
!  integer :: nkpt
!  ! Number of k-point in the IBZ.
  !character(len=12) :: whichmatlu
  ! describe the type of local matrix computed (greenDFT, etc..)
!
  integer :: gpu_option
  ! Wether ks and matlu are stored on GPU
!
  integer :: ndat
  ! Number of elements computed in batch
!
  integer :: nspinor
  ! Number of spinorial components
!
  integer :: nsppol
  ! Number of polarizations

  complex(dp), allocatable :: mat(:,:,:)
  ! Local quantity

 end type matlu_type

!----------------------------------------------------------------------


CONTAINS  !========================================================================================
!!***

!!****f* m_matlu/init_matlu
!! NAME
!! init_matlu
!!
!! FUNCTION
!!  Allocate variables used in type matlu_type.
!!
!! INPUTS
!!  natom   = number of atoms
!!  nspinor = number of spinorial components
!!  nsppol  = number of polarisation components
!!  lpawu_natom(natom) = value of lpawu for every atom
!!  matlu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!
!! OUTPUTS
!!  matlu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!
!! SOURCE

subroutine init_matlu(natom,nspinor,nsppol,lpawu_natom,matlu,gpu_option,ndat)

!Arguments ------------------------------------
 integer, intent(in) :: natom,nspinor,nsppol
 integer, intent(in) :: lpawu_natom(natom)
 integer, intent(in), optional :: gpu_option,ndat
 type(matlu_type), target, intent(inout) :: matlu(natom)
!Local variables ------------------------------------
 integer :: iatom,lpawu,ndim,l_gpu_option,l_ndat
 complex(dp), ABI_CONTIGUOUS pointer :: mat(:,:,:)
!************************************************************************

 l_gpu_option=ABI_GPU_DISABLED; if(present(gpu_option)) l_gpu_option=gpu_option
 l_ndat=1; if(present(ndat)) l_ndat=ndat
! matlu%mband       = mband
! matlu%dmftbandf   = dmftbandf
! matlu%dmftbandi   = dmftbandi
! matlu%nkpt        = nkpt
! matlu%mbandc  = 0
 do iatom=1,natom

   lpawu = lpawu_natom(iatom)
   matlu(iatom)%lpawu   = lpawu
   matlu(iatom)%nspinor = nspinor
   matlu(iatom)%nsppol  = nsppol
   matlu(iatom)%ndat  = l_ndat
   matlu(iatom)%gpu_option  = l_gpu_option
   if (lpawu == -1) cycle
   ndim = (2*lpawu+1) * nspinor
   ABI_MALLOC(matlu(iatom)%mat,(ndim,ndim,nsppol))
   if(l_gpu_option==ABI_GPU_DISABLED) then
     matlu(iatom)%mat(:,:,:) = czero
   else if(l_gpu_option==ABI_GPU_OPENMP) then
     matlu(iatom)%mat(:,:,:) = czero
     mat => matlu(iatom)%mat ! array of structs in OpenMP loosely supported
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET ENTER DATA MAP(alloc:mat)
#endif
     call gpu_set_to_zero_complex(matlu(iatom)%mat, int(nsppol,c_size_t)*ndim*ndim)
   end if

 end do ! iatom

end subroutine init_matlu
!!***

!!****f* m_matlu/zero_matlu
!! NAME
!! zero_matlu
!!
!! FUNCTION
!!  Set the elements of matlu to 0.
!!
!! INPUTS
!!  matlu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  natom = number of atoms
!!  onlynondiag = set all the off-diagonal elements to 0
!!  onlyimag = set the imaginary part to 0
!!
!! OUTPUT
!!  matlu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  err = maximal off-diagonal/imaginary element that is neglected
!!
!! SOURCE

subroutine zero_matlu(matlu,natom,onlynondiag,onlyimag,err)

!Arguments ------------------------------------
 integer, intent(in) :: natom
 type(matlu_type), intent(inout) :: matlu(natom)
 integer, optional, intent(in) :: onlyimag,onlynondiag
 real(dp), optional, intent(out) :: err
!Local variables-------------------------------
 integer :: iatom,im,im1,isppol
 integer :: lpawu,ndim,nspinor,nsppol,tndim
 real(dp) :: err_
!*********************************************************************

 nspinor = matlu(1)%nspinor
 nsppol  = matlu(1)%nsppol

 if (present(err)) err = zero

 do iatom=1,natom
   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   if (present(onlynondiag)) then
     ndim  = 2*lpawu + 1
     tndim = nspinor * ndim
     do isppol=1,nsppol
       do im1=1,tndim
         do im=1,tndim
           if (im /= im1) then
             if (present(err)) then
               err_ = abs(matlu(iatom)%mat(im,im1,isppol))
               if (err_ > err) err = err_
             end if
             matlu(iatom)%mat(im,im1,isppol) = czero
           end if ! im/=im1
         end do ! im
       end do ! im1
     end do ! isppol
   else if (present(onlyimag)) then
     if (present(err)) err = maxval(abs(aimag(matlu(iatom)%mat(:,:,:))))
     matlu(iatom)%mat(:,:,:) = cmplx(dble(matlu(iatom)%mat(:,:,:)),zero,kind=dp)
   else
     matlu(iatom)%mat(:,:,:) = czero
   end if ! onlynondiag
 end do ! iatom

end subroutine zero_matlu
!!***

!!****f* m_matlu/destroy_matlu
!! NAME
!! destroy_matlu
!!
!! FUNCTION
!!  Deallocate matlu
!!
!! INPUTS
!!  matlu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  natom = number of atoms
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_matlu(matlu,natom)

!Arguments ------------------------------------
 integer, intent(in) :: natom
 type(matlu_type),target, intent(inout) :: matlu(natom)
!Local variables-------------------------------
 integer :: iatom
 complex(dp), ABI_CONTIGUOUS pointer :: mat(:,:,:)
! *********************************************************************

 do iatom=1,natom
   mat => matlu(iatom)%mat ! array of structs in OpenMP loosely supported
   if(matlu(iatom)%gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET EXIT DATA MAP(delete:mat)
#endif
   end if
   ABI_SFREE(matlu(iatom)%mat)
 end do ! iatom

end subroutine destroy_matlu
!!***

!!****f* m_matlu/copy_matlu
!! NAME
!! copy_matlu
!!
!! FUNCTION
!!  Copy mat1 into mat2
!!
!! INPUTS
!!  mat1 <type(matlu_type)>= density matrix nmat1 in the local orbital basis and related variables
!!  natom = number of atoms
!!  opt_diag = if present, only copy the diagonal elements (the off-diagonal elements are not set to 0)
!!  opt_non_diag = if present, only copy the off-diagonal elements
!!  opt_re = if present, only copy the real part
!!
!! OUTPUT
!!  mat2 <type(matlu_type)>= density matrix nmat2 in the local orbital basis and related variables
!!
!! SOURCE

subroutine copy_matlu(mat1,mat2,natom,opt_diag,opt_non_diag,opt_re)

!Arguments ------------------------------------
 integer, intent(in) :: natom
 type(matlu_type), intent(in) :: mat1(natom)
 type(matlu_type), intent(inout) :: mat2(natom) !vz_i
 integer, optional, intent(in) :: opt_diag,opt_non_diag,opt_re
!Local variables-------------------------------
 integer :: iatom,im,im1,isppol,lpawu,ndim,nspinor,nsppol
! *********************************************************************

 nspinor = mat1(1)%nspinor
 nsppol  = mat1(1)%nsppol

 do iatom=1,natom

   lpawu = mat1(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim = (2*lpawu+1) * nspinor

   !if both matrix are on GPU and no opt is provided, perform copy on GPU
   !Other cases can be handled through OpenMP kernels on GPU but no use case exists.
   if(mat1(iatom)%gpu_option==ABI_GPU_OPENMP .and. mat2(iatom)%gpu_option==ABI_GPU_OPENMP &
   &    .and. .not. (present(opt_diag) .or. present(opt_non_diag) .or. present(opt_re))) then
#ifdef HAVE_OPENMP_OFFLOAD
     call gpu_copy_complex(mat2(iatom)%mat, mat1(iatom)%mat, int(nsppol, c_size_t)*ndim*ndim)
#endif
     cycle
   end if


#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET UPDATE FROM(mat1(iatom)%mat) IF(mat1(iatom)%gpu_option==ABI_GPU_OPENMP)
   !$OMP TARGET UPDATE FROM(mat2(iatom)%mat) IF(mat2(iatom)%gpu_option==ABI_GPU_OPENMP)
#endif
   if (present(opt_diag)) then
     do isppol=1,nsppol
       do im=1,ndim
         mat2(iatom)%mat(im,im,isppol) = mat1(iatom)%mat(im,im,isppol)
       end do ! im
     end do ! isppol
   else if (present(opt_non_diag)) then
     do isppol=1,nsppol
       do im1=1,ndim
         do im=1,ndim
           if (im /= im1) mat2(iatom)%mat(im,im1,isppol) = mat1(iatom)%mat(im,im1,isppol)
         end do ! im
       end do ! im1
     end do ! isppol
   else if (present(opt_re)) then
     mat2(iatom)%mat(:,:,:) = cmplx(dble(mat1(iatom)%mat(:,:,:)),zero,kind=dp)
   else
     mat2(iatom)%mat(:,:,:) = mat1(iatom)%mat(:,:,:)
   end if ! opt
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET UPDATE TO(mat2(iatom)%mat) IF(mat2(iatom)%gpu_option==ABI_GPU_OPENMP)
#endif

 end do ! iatom

!   do iatom=1,natom
!    lpawu=nmat1(iatom)%lpawu
!    if(lpawu.ne.-1) then
!     nmat2(iatom)%mat=nmat1(iatom)%mat
!    endif
!   enddo

end subroutine copy_matlu
!!***

!!****f* m_matlu/copy_matlu_from_ndat
!! NAME
!! copy_matlu_from_ndat
!!
!! FUNCTION
!!  Copy matlu1 into matlu2
!!
!! INPUTS
!!  maltu1 <type(matlu_type)>= density matrix matlu1 in the local orbital basis and related variables
!!  natom = number of atoms
!!
!! OUTPUT
!!  maltu2 <type(matlu_type)>= density matrix matlu2 in the local orbital basis and related variables
!!
!! SOURCE

subroutine copy_matlu_from_ndat(mat1,mat2,natom,ndat,idat,opt_diag,opt_non_diag,opt_re)

!Arguments ------------------------------------
!type
 integer, intent(in) :: natom,ndat,idat
 type(matlu_type),intent(in) :: mat1(natom)
 type(matlu_type),intent(inout) :: mat2(natom) !vz_i
 integer, optional, intent(in) :: opt_diag,opt_non_diag,opt_re

!Local variables-------------------------------
 integer :: iatom,isppol,im1,im,ndim,nsppol,nspinor,lpawu
! *********************************************************************

 ABI_CHECK(mat1(1)%nsppol==mat2(1)%nsppol*ndat, "bad ndat value")
 nspinor = mat1(1)%nspinor
 nsppol  = mat2(1)%nsppol
 do iatom=1,natom
   lpawu = mat1(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim = (2*lpawu+1) * nspinor

   if (present(opt_diag)) then
     do isppol=1,nsppol
       do im=1,ndim
         mat2(iatom)%mat(im,im,isppol) = mat1(iatom)%mat(im,im,(isppol-1)*ndat+idat)
       end do ! im
     end do ! isppol
   else if (present(opt_non_diag)) then
     do isppol=1,nsppol
       do im1=1,ndim
         do im=1,ndim
           if (im /= im1) mat2(iatom)%mat(im,im1,isppol) = mat1(iatom)%mat(im,im1,(isppol-1)*ndat+idat)
         end do ! im
       end do ! im1
     end do ! isppol

   else if (present(opt_re)) then
     do isppol=1,nsppol
       mat2(iatom)%mat(:,:,isppol) = cmplx(dble(mat1(iatom)%mat(:,:,(isppol-1)*ndat+idat)),zero,kind=dp)
     end do ! isppol
   else
     do isppol=1,nsppol
       mat2(iatom)%mat(:,:,isppol) = mat1(iatom)%mat(:,:,(isppol-1)*ndat+idat)
     end do ! isppol
   end if ! opt
 enddo ! iatom

end subroutine copy_matlu_from_ndat
!!***

!!****f* m_matlu/copy_matlu_to_ndat
!! NAME
!! copy_matlu_to_ndat
!!
!! FUNCTION
!!  Copy matlu1 into matlu2
!!
!! INPUTS
!!  maltu1 <type(matlu_type)>= density matrix matlu1 in the local orbital basis and related variables
!!  natom = number of atoms
!!
!! OUTPUT
!!  maltu2 <type(matlu_type)>= density matrix matlu2 in the local orbital basis and related variables
!!
!! SOURCE

subroutine copy_matlu_to_ndat(mat1,mat2,natom,ndat,idat,opt_diag,opt_non_diag,opt_re)

!Arguments ------------------------------------
!type
 integer, intent(in) :: natom,ndat,idat
 type(matlu_type),intent(in) :: mat1(natom)
 type(matlu_type),intent(inout) :: mat2(natom) !vz_i
 integer, optional, intent(in) :: opt_diag,opt_non_diag,opt_re

!Local variables-------------------------------
 integer :: iatom,isppol,im1,im,ndim,nsppol,nspinor,lpawu
! *********************************************************************


 ABI_CHECK(mat1(1)%nsppol*ndat==mat2(1)%nsppol, "bad ndat value")
 nspinor = mat1(1)%nspinor
 nsppol  = mat1(1)%nsppol
 do iatom=1,natom
   lpawu = mat1(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim = (2*lpawu+1) * nspinor

   if (present(opt_diag)) then
     do isppol=1,nsppol
       do im=1,ndim
         mat2(iatom)%mat(im,im,(isppol-1)*ndat+idat) = mat1(iatom)%mat(im,im,isppol)
       end do ! im
     end do ! isppol
   else if (present(opt_non_diag)) then
     do isppol=1,nsppol
       do im1=1,ndim
         do im=1,ndim
           if (im /= im1) mat2(iatom)%mat(im,im1,(isppol-1)*ndat+idat) = mat1(iatom)%mat(im,im1,isppol)
         end do ! im
       end do ! im1
     end do ! isppol

   else if (present(opt_re)) then
     do isppol=1,nsppol
       mat2(iatom)%mat(:,:,(isppol-1)*ndat+idat) = cmplx(dble(mat1(iatom)%mat(:,:,isppol)),zero,kind=dp)
     end do ! isppol
   else
     do isppol=1,nsppol
       mat2(iatom)%mat(:,:,(isppol-1)*ndat+idat) = mat1(iatom)%mat(:,:,isppol)
     end do ! isppol
   end if ! opt
 enddo ! iatom

end subroutine copy_matlu_to_ndat
!!***

!!****f* m_matlu/print_matlu
!! NAME
!! print_matlu
!!
!! FUNCTION
!!  Print matlu
!!
!! INPUTS
!!  matlu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  natom= number of atoms
!!  prtopt=   /=0   print matlu
!!            >=5   print matlu in n,mx,my,mz representation
!!  opt_diag=   0   print non diagonal matrix (real or complex according to nspinor)
!!             -1   print non diagonal complex matrix
!!            >=1   print diagonal matrix (real or complex according to nspinor)
!!  opt_ab_out=  0  print matrix on std_out
!!             /=0  print matrix on ab_out
!!  opt_exp=        write in exponent format if present
!!  argout=         output unit
!!  compl=      1   print complex matrix
!!            /=1   only print complex matrix if nspinor=2
!!
!! OUTPUT
!!
!! SOURCE

subroutine print_matlu(matlu,natom,prtopt,opt_diag,opt_ab_out,opt_exp,argout,compl)

!Arguments ------------------------------------
 integer, intent(in):: natom,prtopt
 type(matlu_type), intent(in) :: matlu(natom)
 integer, optional, intent(in) :: opt_diag,opt_ab_out,opt_exp,argout,compl
!Local variables-------------------------------
 integer :: arg_out,iatom,im,im1,ispinor,ispinor1,isppol,lpawu
 integer :: ndim,nspinor,nsppol,optab_out,optdiag
 logical :: testcmplx,testcmplx_
 complex(dp), allocatable :: mat_nmrep(:,:)
 character(len=500) :: message
 character(len=4) :: mode_paral,tag_at
 character(len=9), parameter :: dspinm(2,2) = RESHAPE((/"n        ","mx       ","my       ","mz       "/),(/2,2/))
! *********************************************************************

 arg_out    = ab_out
 mode_paral = 'COLL'
 optab_out  = 0
 optdiag    = 0

 if (present(opt_diag)) optdiag = opt_diag
 if (present(opt_ab_out)) optab_out = opt_ab_out
 if (optab_out == 0) arg_out = std_out

 if (present(argout)) then
  arg_out    = argout
  mode_paral = 'PERS'
 end if

 nspinor    = matlu(1)%nspinor
 nsppol     = matlu(1)%nsppol
 testcmplx_ = (nspinor == 2)
 if (present(compl)) testcmplx_ = (nspinor == 2) .or. (compl == 1)

 do iatom=1,natom

   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim = 2*lpawu + 1

   write(tag_at,'(i4)') iatom
   write(message,'(3a)') ch10,'   -------> For Correlated Atom ',adjustl(tag_at)
   call wrtout(arg_out,message,mode_paral)

   testcmplx = testcmplx_
   if (maxval(abs(aimag(matlu(iatom)%mat(:,:,:)))) > tol5) testcmplx = .true.

   !do isppol=1,nsppol
   !  if (present(opt_ab_out) .and. nsppol == 2) then
   !    noccspin = zero
   !    do im=1,ndim
   !      noccspin = noccspin + REAL(matlu(iatom)%mat(im,im,isppol))
   !    end do
       !write(message,fmt='(7x,a,i3,a,f10.5)') ". Occ. for lpawu and for spin",isppol," =",noccspin
       !call wrtout(arg_out, message,mode_paral)
   !  end if
   !end do ! isppol

   do isppol=1,nsppol
     if (nspinor == 1) then
       write(message,'(a,10x,a,1x,i1)') ch10,'-- polarization spin component',isppol
       call wrtout(arg_out,message,mode_paral)
     end if ! nspinor=1
     do ispinor=1,nspinor
       do ispinor1=1,nspinor
         if (nspinor == 2) then
           write(message,'(a,10x,a,i1,1x,i1)') ch10,'-- spin components ',ispinor,ispinor1
           call wrtout(arg_out,message,mode_paral)
         end if
         if (optdiag <= 0) then
           do im1=1,ndim
             if (optdiag == 0) then
               if ((.not. testcmplx) .and. (abs(prtopt) > 0)) then
                 if (present(opt_exp)) then
                   write(message,'(5x,20e24.14)') (dble(matlu(iatom)%mat(im1+(ispinor-1)*ndim,im+(ispinor1-1)*ndim,isppol)),im=1,ndim)
!                  call wrtout(arg_out,  message,mode_paral)
!                  write(message,'(5x,20e20.14)') (REAL(sqrt(matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1))),m=1,2*lpawu+1)
!                  call wrtout(arg_out,  message,mode_paral)
!                  write(message,'(5x,20e20.14)') (REAL(1.d0/sqrt(matlu(iatom)%mat(m,m,isppol,ispinor,ispinor1))),m=1,2*lpawu+1)
                 else
                   write(message,'(5x,20f10.5)') (dble(matlu(iatom)%mat(im1+(ispinor-1)*ndim,im+(ispinor1-1)*ndim,isppol)),im=1,ndim)
                 end if ! opt_exp
               else if (testcmplx .and. (abs(prtopt) > 0)) then
                 if (present(opt_exp)) then
                   if (opt_exp == 2) then
                     write(message,'(5x,14(2e18.10,1x))') ((matlu(iatom)%mat(im1+(ispinor-1)*ndim,im+(ispinor1-1)*ndim,isppol)),im=1,ndim)
                   else
                     write(message,'(5x,14(2e14.4,2x))') ((matlu(iatom)%mat(im1+(ispinor-1)*ndim,im+(ispinor1-1)*ndim,isppol)),im=1,ndim)
                   end if ! opt_exp=2
                 else
                   write(message,'(5x,14(2f9.5,2x))') ((matlu(iatom)%mat(im1+(ispinor-1)*ndim,im+(ispinor1-1)*ndim,isppol)),im=1,ndim)
                 end if ! opt_exp
!&               write(message,'(5x,14(2f15.11,2x))')((matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
               end if ! testcmplx
             else if (optdiag == -1) then
               write(message,'(5x,14(2f10.5,2x))') ((matlu(iatom)%mat(im1+(ispinor-1)*ndim,im+(ispinor1-1)*ndim,isppol)),im=1,ndim)
             end if ! optdiag
             call wrtout(arg_out,message,mode_paral)
           end do ! im1
         else if (optdiag >= 1) then
           if ((.not. testcmplx) .and. (abs(prtopt) > 0)) write(message,'(5x,20f10.5)') &
               & (dble(matlu(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor1-1)*ndim,isppol)),im=1,ndim)
           if (testcmplx .and. (abs(prtopt) > 0)) write(message,'(5x,14(2f9.5,2x))') &
               & ((matlu(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor1-1)*ndim,isppol)),im=1,ndim)
!            write(std_out,'(5x,14(2f9.5,2x))')((matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
           call wrtout(arg_out,message,mode_paral)
         end if ! optdiag
       end do ! ispinor1
     end do ! ispinor
     if (nspinor == 2 .and. prtopt >= 5) then
       ABI_MALLOC(mat_nmrep,(2*ndim,2*ndim)) ! Put matlu in n,mx,my,mz representation
       do im=1,ndim
         do im1=1,ndim
           mat_nmrep(im1,im) = matlu(iatom)%mat(im1,im,isppol) + matlu(iatom)%mat(im1+ndim,im+ndim,isppol)  ! n
           mat_nmrep(im1+ndim,im+ndim) = matlu(iatom)%mat(im1,im,isppol) - matlu(iatom)%mat(im1+ndim,im+ndim,isppol) ! mz
           mat_nmrep(im1+ndim,im) = matlu(iatom)%mat(im1,im+ndim,isppol) + matlu(iatom)%mat(im1+ndim,im,isppol)  ! mx
           mat_nmrep(im1,im+ndim) = (matlu(iatom)%mat(im1,im+ndim,isppol)-matlu(iatom)%mat(im+ndim,im,isppol)) * j_dpc ! my
         end do ! im1
       end do ! im
       do ispinor=1,nspinor
         do ispinor1=1,nspinor
           write(message,'(a,10x,2a)') ch10,'-- spin components',dspinm(ispinor1,ispinor)
           call wrtout(arg_out,message,mode_paral)
           do im1=1,ndim
             write(message,'(5x,14(2f9.5,2x))') ((mat_nmrep(im1+(ispinor1-1)*ndim,im+(ispinor-1)*ndim)),im=1,ndim)
             call wrtout(arg_out,message,mode_paral)
           end do ! im1
         end do ! ispinor1
       end do ! ispinor
       ABI_FREE(mat_nmrep)
     end if ! nspinor=2 and prtopt >=5
   end do ! isppol
!     if(nsppol==1.and.nspinor==1) then
!       write(message,'(a,10x,a,i3,a)')  ch10,'-- polarization spin component',isppol+1,' is identical'
!       call wrtout(arg_out,  message,mode_paral)
!     endif
 end do ! iatom

end subroutine print_matlu
!!***

!!****f* m_matlu/sym_matlu
!! NAME
!! sym_matlu
!!
!! FUNCTION
!! Symmetrize local quantity.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gloc(natom) <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!  gloc(natom) <type(matlu_type)>= density matrix symmetrized in the local orbital basis and related variables
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine sym_matlu(gloc,paw_dmft)

!Arguments ------------------------------------
 type(paw_dmft_type), target, intent(in) :: paw_dmft
 type(matlu_type), target, intent(inout) :: gloc(paw_dmft%natom)
!Local variables-------------------------------
 integer :: at_indx,iatom,irot,isppol,lpawu,m1,m2,mu,natom
 integer :: ndim,ndim_max,nspinor,nsppol,nsym,nu,gpu_option
 complex(dp), allocatable :: gloc_tmp(:,:,:),gloc_tmp2(:,:,:)
 complex(dp), allocatable :: gloc_tmp3(:,:,:,:),gloc_tmp4(:,:,:,:)
 type(matlu_type), allocatable, target :: gloc_nmrep(:),glocsym(:)
 complex(dp), ABI_CONTIGUOUS pointer :: zarot(:,:,:,:),gloc_mat(:,:,:),glocsym_mat(:,:,:)
 real(dp), ABI_CONTIGUOUS pointer :: symrec_cart(:,:,:)
 complex(dp) :: ratio

 natom    = paw_dmft%natom
 ndim_max = 2*paw_dmft%maxlpawu + 1
 nspinor  = paw_dmft%nspinor
 nsppol   = gloc(1)%nsppol
 nsym     = paw_dmft%nsym
 gpu_option = gloc(1)%gpu_option
 zarot    => paw_dmft%zarot
 ratio = dcmplx(1.0_dp/nsym,0.0_dp)

 !zarot       => paw_dmft%zarot(:,1:ndim,irot,lpawu+1)
 ABI_MALLOC(glocsym,(natom))
#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET ENTER DATA MAP(alloc:zarot) IF(gpu_option==ABI_GPU_OPENMP)
 !$OMP TARGET UPDATE TO(zarot) IF(gpu_option==ABI_GPU_OPENMP)
#else
  ABI_UNUSED((/m1,m2/))
#endif

!=========  Case nspinor ==1 ========================

 if (nspinor == 1) then

   ABI_MALLOC(gloc_tmp,(ndim_max,ndim_max*nsppol,nsym))
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(alloc:gloc_tmp) IF(gpu_option==ABI_GPU_OPENMP)
#endif

   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu(:),glocsym(:),gpu_option=gpu_option)

   do iatom=1,natom

     lpawu = gloc(iatom)%lpawu
     if (lpawu == -1) cycle
     ndim = 2*lpawu + 1
     glocsym_mat => glocsym(iatom)%mat

     ABI_MALLOC(gloc_tmp2,(ndim,ndim*nsppol,nsym))
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET ENTER DATA MAP(alloc:gloc_tmp2) IF(gpu_option==ABI_GPU_OPENMP)
#endif


     if(gpu_option==ABI_GPU_DISABLED) then

       do irot=1,nsym
         at_indx = paw_dmft%indsym(irot,iatom)
         gloc_mat    => gloc(at_indx)%mat

         do isppol=1,nsppol
           call abi_xgemm("n","n",ndim,ndim,ndim,cone,gloc_mat(:,:,isppol),ndim,&
                        & zarot(:,1:ndim,irot,lpawu+1),ndim_max,czero,gloc_tmp(:,1+ndim*(isppol-1):ndim*isppol,irot),ndim_max)
         end do ! isppol

         call abi_xgemm("t","n",ndim,ndim*nsppol,ndim,cone,zarot(:,1:ndim,irot,lpawu+1),ndim_max,&
                      & gloc_tmp(:,:,irot),ndim_max,czero,gloc_tmp2(:,:,irot),ndim)

         do isppol=1,nsppol
           glocsym_mat(:,:,isppol) = glocsym_mat(:,:,isppol) + gloc_tmp2(:,1+ndim*(isppol-1):ndim*isppol,irot)
         end do ! isppol

       end do ! irot

       glocsym_mat(:,:,:) = glocsym_mat(:,:,:) / dble(nsym)

     else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD

       do irot=1,nsym
         at_indx = paw_dmft%indsym(irot,iatom)
         gloc_mat    => gloc(at_indx)%mat

         !$OMP TARGET DATA USE_DEVICE_ADDR(gloc_tmp,zarot,gloc_mat)
         call abi_gpu_xgemm_strided(2,"n","n",ndim,ndim,ndim,cone,&
         &    c_loc(gloc_mat(:,:,:)),ndim,ndim*ndim,&
         &    c_loc(zarot(:,1:ndim,irot,lpawu+1)),ndim_max,0,czero,&
         &    c_loc(gloc_tmp(:,:,irot)),ndim_max,ndim*ndim,nsppol,async=.true.,stream_id=irot)
         !$OMP END TARGET DATA
       end do ! irot
       call gpu_device_synchronize()

       !$OMP TARGET DATA USE_DEVICE_ADDR(gloc_tmp,zarot,gloc_tmp2)
       call abi_gpu_xgemm_strided(2,"t","n",ndim,ndim*nsppol,ndim,cone,&
       &    c_loc(zarot(:,:,:,lpawu+1)),ndim_max,ndim_max*ndim_max,&
       &    c_loc(gloc_tmp(:,:,:)),ndim_max,ndim_max*ndim_max*nsppol,czero,&
       &    c_loc(gloc_tmp2(:,:,:)),ndim,ndim*ndim*nsppol,nsym)
       !$OMP END TARGET DATA

       !$OMP TARGET TEAMS DISTRIBUTE MAP(to:glocsym_mat,gloc_tmp2) PRIVATE(isppol)
       do isppol=1,nsppol
         !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(m1,m2,irot)
         do m2=1,ndim
           do m1=1,ndim
             do irot=1,nsym
               glocsym_mat(m1,m2,isppol) = glocsym_mat(m1,m2,isppol) + gloc_tmp2(m1,m2+ndim*(isppol-1),irot)
             end do ! irot
           end do ! m1
         end do  ! m2
       end do ! isppol

       !$OMP TARGET DATA USE_DEVICE_ADDR(glocsym_mat)
       call abi_gpu_xscal(2, ndim*ndim*nsppol, ratio, c_loc(glocsym_mat), 1)
       !$OMP END TARGET DATA

#endif
     end if

#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET EXIT DATA MAP(delete:gloc_tmp2) IF(gpu_option==ABI_GPU_OPENMP)
#endif
     ABI_FREE(gloc_tmp2)

   end do ! iatom

   !==  Put glocsym into gloc
   call copy_matlu(glocsym(:),gloc(:),natom)

#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET EXIT DATA MAP(delete:gloc_tmp) IF(gpu_option==ABI_GPU_OPENMP)
#endif
   ABI_FREE(gloc_tmp)
!=========  Case nspinor ==2 ========================

 else

   symrec_cart => paw_dmft%symrec_cart(:,:,:)
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(to:symrec_cart) IF(gpu_option==ABI_GPU_OPENMP)
#endif
   !== Allocate temporary arrays
   ABI_MALLOC(gloc_nmrep,(natom))
   call init_matlu(natom,1,4*nsppol,paw_dmft%lpawu(:),glocsym(:),gpu_option=gpu_option)
   call init_matlu(natom,1,4*nsppol,paw_dmft%lpawu(:),gloc_nmrep(:),gpu_option=gpu_option)

   ! Put gloc into gloc_nmrep (density and magnetization representation)
   ! gloc_nmrep(iatom)%mat(:,:,i) = n,mx,my,mz for i=1,2,3,4 respectively
   call chg_repr_matlu(gloc(:),gloc_nmrep(:),natom,1,1)

  !==  Do the sum over symmetrized density matrix (in n,m repr)
   do iatom=1,natom

     lpawu = gloc(iatom)%lpawu
     if (lpawu == -1) cycle
     ndim = 2*lpawu + 1
     glocsym_mat => glocsym(iatom)%mat


     ABI_MALLOC(gloc_tmp3,(ndim,ndim,4*nsppol,nsym))
     ABI_MALLOC(gloc_tmp4,(ndim,ndim,4*nsppol,nsym))

     if(gpu_option==ABI_GPU_DISABLED) then
       do irot=1,nsym

         at_indx = paw_dmft%indsym(irot,iatom)
         gloc_mat    => gloc_nmrep(at_indx)%mat

         do isppol=1,nsppol
           do mu=1,4 ! Symmetrize density and magnetization

             call abi_xgemm("n","n",ndim,ndim,ndim,cone,gloc_mat(:,:,mu+(isppol-1)*4),ndim, &
                          & zarot(:,1:ndim,irot,lpawu+1),ndim_max,czero,gloc_tmp3(:,:,mu+(isppol-1)*4,irot),ndim)

           end do ! mu
         end do ! isppol
       end do ! irot

       do irot=1,nsym
         call abi_zgemm_2dd("t","n",ndim,ndim*4*nsppol,ndim,cone,zarot(:,1:ndim,irot,lpawu+1),ndim_max,&
                      & gloc_tmp3(:,:,:,irot),ndim,czero,gloc_tmp4(:,:,:,irot),ndim)

       end do ! irot

       do irot=1,nsym
         do isppol=1,nsppol
           glocsym_mat(:,:,1+(isppol-1)*4) = glocsym_mat(:,:,1+(isppol-1)*4) + gloc_tmp4(:,:,1+(isppol-1)*4,irot)
         end do ! isppol
       end do ! irot

         ! Symmetrize magnetization

       do irot=1,nsym
         do isppol=1,nsppol
           do nu=2,4
             do mu=2,4
               glocsym_mat(:,:,mu+(isppol-1)*4) = glocsym_mat(:,:,mu+(isppol-1)*4) + &
                 &    symrec_cart(mu-1,nu-1,irot)*gloc_tmp4(:,:,nu+(isppol-1)*4,irot)
             end do ! mu
           end do ! nu
         end do ! isppol
       end do ! irot

    !  ==  Normalize sum
       glocsym_mat(:,:,:) = glocsym_mat(:,:,:) / dble(nsym)

     else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET ENTER DATA MAP(alloc:gloc_tmp3,gloc_tmp4)

       do irot=1,nsym

         at_indx = paw_dmft%indsym(irot,iatom)
         gloc_mat    => gloc_nmrep(at_indx)%mat

         !$OMP TARGET DATA USE_DEVICE_ADDR(gloc_tmp3,zarot,gloc_mat)
         call abi_gpu_xgemm_strided(2,"n","n",ndim,ndim,ndim,cone,&
         &    c_loc(gloc_mat(:,:,:)),ndim,ndim*ndim,&
         &    c_loc(zarot(:,1:ndim,irot,lpawu+1)),ndim_max,0,czero,&
         &    c_loc(gloc_tmp3(:,:,:,irot)),ndim,ndim*ndim,4*nsppol,async=.true.,stream_id=irot)
         !$OMP END TARGET DATA
       end do ! irot
       call gpu_device_synchronize()


       !$OMP TARGET DATA USE_DEVICE_ADDR(gloc_tmp3,zarot,gloc_tmp4)
       call abi_gpu_xgemm_strided(2,"t","n",ndim,ndim*4*nsppol,ndim,cone,&
       &    c_loc(zarot(:,:,:,lpawu+1)),ndim_max,ndim_max*ndim_max,&
       &    c_loc(gloc_tmp3(:,:,:,:)),ndim,ndim*ndim*4*nsppol,czero,&
       &    c_loc(gloc_tmp4(:,:,:,:)),ndim,ndim*ndim*4*nsppol,nsym)
       !$OMP END TARGET DATA


       !$OMP TARGET TEAMS DISTRIBUTE MAP(to:glocsym_mat,gloc_tmp4) PRIVATE(isppol)
       do isppol=1,nsppol
         !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(m1,m2,irot)
         do m2=1,ndim
           do m1=1,ndim
             do irot=1,nsym
               glocsym_mat(m1,m2,1+(isppol-1)*4) = glocsym_mat(m1,m2,1+(isppol-1)*4) + gloc_tmp4(m1,m2,1+(isppol-1)*4,irot)
             end do ! irot
           end do ! m1
         end do  ! m2
       end do ! isppol

         ! Symmetrize magnetization

       !$OMP TARGET TEAMS DISTRIBUTE MAP(to:glocsym_mat,gloc_tmp4,symrec_cart) PRIVATE(isppol,irot)
       do isppol=1,nsppol
         !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(nu,mu,m2,m1)
         do mu=2,4
           do m2=1,ndim
             do m1=1,ndim
               do irot=1,nsym
                 do nu=2,4
                   glocsym_mat(m1,m2,mu+(isppol-1)*4) = glocsym_mat(m1,m2,mu+(isppol-1)*4) + &
                   &    symrec_cart(mu-1,nu-1,irot)*gloc_tmp4(m1,m2,nu+(isppol-1)*4,irot)
                 end do ! m1
               end do  ! m2
             end do ! mu
           end do ! nu
         end do ! isppol
       end do ! irot

    !  ==  Normalize sum
       !$OMP TARGET DATA USE_DEVICE_ADDR(glocsym_mat)
       call abi_gpu_xscal(2, ndim*ndim*4*nsppol, ratio, c_loc(glocsym_mat), 1)
       !$OMP END TARGET DATA

       !$OMP TARGET EXIT DATA MAP(delete:gloc_tmp3,gloc_tmp4)
#endif
     end if ! gpu_option

     ABI_FREE(gloc_tmp3)
     ABI_FREE(gloc_tmp4)

   end do ! iatom

!==  Compute back density matrix in upup dndn updn dnup representation
   call chg_repr_matlu(gloc(:),glocsym(:),natom,-1,1)

   call destroy_matlu(gloc_nmrep(:),natom)
   ABI_FREE(gloc_nmrep)

#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET EXIT DATA MAP(delete:symrec_cart) IF(gpu_option==ABI_GPU_OPENMP)
#endif
  !==============end of nspinor=2 case ===========
 end if ! nspinor

 call destroy_matlu(glocsym(:),natom)
 ABI_FREE(glocsym)
#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET EXIT DATA MAP(delete:zarot) IF(gpu_option==ABI_GPU_OPENMP)
#endif

 !mt2g(1)=1
 !mt2g(2)=2
 !mt2g(3)=4
 !mx2my2d=5
 !t2g=paw_dmft%dmftqmc_t2g
 !x2my2d=paw_dmft%dmftqmc_x2my2d

! DBG_ENTER("COLL")

 !ci=cone
 !nspinor=gloc(1)%nspinor
 !nsppol=gloc(1)%nsppol
 !natom=cryst_struc%natom

 !ABI_MALLOC(glocnm,(natom))
 !ABI_MALLOC(glocnms,(natom))
 !ABI_MALLOC(glocsym,(natom))
 !ABI_MALLOC(lpawu_natom,(natom))
 !lpawu_natom(1:natom)=gloc(1:natom)%lpawu ! If gloc(1:natom)%lpawu is directly used in the next three lines, warnings are generated by some compilers.
 !call init_matlu(natom,nspinor,nsppol,lpawu_natom,glocnm)
 !call init_matlu(natom,nspinor,nsppol,lpawu_natom,glocnms)
 !call init_matlu(natom,nspinor,nsppol,lpawu_natom,glocsym)
 !ABI_FREE(lpawu_natom)

!=========  Case nspinor ==1 ========================

 !if (nspinor==1) then
 ! ispinor=1
 ! ispinor1=1
 ! do iatom=1,cryst_struc%natom
 !  do isppol=1,nsppol
 !   if(gloc(iatom)%lpawu/=-1) then
 !    lpawu=gloc(iatom)%lpawu
 !    do m1=1, 2*lpawu+1
 !     do m2=1, 2*lpawu+1
 !      do irot=1,cryst_struc%nsym
 !       at_indx=cryst_struc%indsym(4,irot,iatom)
 !       do m3=1, 2*lpawu+1
 !        do m4=1, 2*lpawu+1
 !         if(t2g==1) then
 !          m1s=mt2g(m1)
 !          m2s=mt2g(m2)
 !          m3s=mt2g(m3)
 !          m4s=mt2g(m4)
 !          lpawu_zarot=2
 !         else if (x2my2d==1) then
 !          m1s=mx2my2d
 !          m2s=mx2my2d
 !          m3s=mx2my2d
 !          m4s=mx2my2d
 !          lpawu_zarot=2
 !         else
 !          m1s=m1
 !          m2s=m2
 !          m3s=m3
 !          m4s=m4
 !          lpawu_zarot=lpawu
 !         endif
 !         zarot2=pawang%zarot(m3s,m1s,lpawu_zarot+1,irot)*pawang%zarot(m4s,m2s,lpawu_zarot+1,irot)
 !         glocsym(iatom)%mat(m1,m2,isppol,ispinor,ispinor1)=&
!&          glocsym(iatom)%mat(m1,m2,isppol,ispinor,ispinor1)&
!&          +gloc(at_indx)%mat(m3,m4,isppol,ispinor,ispinor1)*zarot2
 !        end do  ! m3
 !       end do  ! m4
 !      end do  ! irot
 !      glocsym(iatom)%mat(m1,m2,isppol,ispinor,ispinor1)=&
!&       glocsym(iatom)%mat(m1,m2,isppol,ispinor,ispinor1)/real(cryst_struc%nsym,kind=dp)
!      end do ! m2
!     end do ! m1
!    endif ! lpawu/=-1
!   end do ! isppol
!  end do ! iatom
!==  Put glocsym into gloc
!  do iatom=1,cryst_struc%natom
!    if(gloc(iatom)%lpawu/=-1) then
!      gloc(iatom)%mat=glocsym(iatom)%mat
!      gloc(iatom)%mat(:,:,1,:,:)=(glocsym(iatom)%mat(:,:,1,:,:) &
!&      + glocsym(iatom)%mat(:,:,2,:,:))/two
!      gloc(iatom)%mat(:,:,2,:,:)= gloc(iatom)%mat(:,:,1,:,:)
!      write(std_out,*) "WARNING: SYM non mag"
!      write(ab_out,*) "WARNING: SYM non mag"
!    endif
!  end do ! iatom

!=========  Case nspinor ==2 ========================

! else if (nspinor==2) then

!== Allocate temporary arrays
!  do iatom=1,cryst_struc%natom
!   if(gloc(iatom)%lpawu/=-1) then
!    ndim=2*gloc(iatom)%lpawu+1
!    ABI_FREE(glocnm(iatom)%mat)
!    ABI_FREE(glocnms(iatom)%mat)
!    ABI_FREE(glocsym(iatom)%mat)
!    ABI_MALLOC(glocnm(iatom)%mat,(ndim,ndim,nsppol,4,1))
!    ABI_MALLOC(glocnms(iatom)%mat,(ndim,ndim,nsppol,4,1))
!    ABI_MALLOC(glocsym(iatom)%mat,(ndim,ndim,nsppol,2,2))
!   endif
!  enddo
!  ABI_MALLOC(symrec_cart,(3,3,cryst_struc%nsym))

!==  Compute symrec_cart
!  do irot=1,cryst_struc%nsym
!   call symredcart(cryst_struc%gprimd,cryst_struc%rprimd,symrec_cart(:,:,irot),cryst_struc%symrec(:,:,irot))
!  end do

!==  Compute density matrix in density and magnetization representation
 ! call chg_repr_matlu(gloc,glocnm,cryst_struc%natom,option=1,prtopt=1)

!==  Do the sum over symetrized density matrix (in n,m repr)
!  isppol=1
!  do iatom=1,cryst_struc%natom
!   if(gloc(iatom)%lpawu/=-1) then
!    lpawu=gloc(iatom)%lpawu
!    ndim=2*gloc(iatom)%lpawu+1
!    do m1=1, 2*lpawu+1
!     do m2=1, 2*lpawu+1
!      sumrho=czero
!      rotmag=czero
!      do irot=1,cryst_struc%nsym
!       summag=czero
!       at_indx=cryst_struc%indsym(4,irot,iatom)
!       do m3=1, 2*lpawu+1
!        do m4=1, 2*lpawu+1
!          if(t2g==1) then
!           m1s=mt2g(m1)
!           m2s=mt2g(m2)
!           m3s=mt2g(m3)
!           m4s=mt2g(m4)
!           lpawu_zarot=2
!!          else if (x2my2d==1) then
!           m1s=mx2my2d
!           m2s=mx2my2d
!           m3s=mx2my2d
!           m4s=mx2my2d
!           lpawu_zarot=2
!          else
!           m1s=m1
!           m2s=m2
!           m3s=m3
!           m4s=m4
!           lpawu_zarot=lpawu
!          endif
!         zarot2=pawang%zarot(m3s,m2s,lpawu_zarot+1,irot)*pawang%zarot(m4s,m1s,lpawu_zarot+1,irot)
!         sumrho=sumrho +  glocnm(at_indx)%mat(m4,m3,isppol,1,1)  * zarot2
!         do mu=1,3
!          summag(mu)=summag(mu) + glocnm(at_indx)%mat(m4,m3,isppol,mu+1,1) * zarot2
!         enddo
!        end do ! m3
!       end do !m4

!       ==  special case of magnetization
 !      do nu=1,3
 !       do mu=1,3
 !        rotmag(mu)=rotmag(mu)+symrec_cart(mu,nu,irot)*summag(nu)
 !       end do
 !      end do
!      write(std_out,'(a,3i4,2x,3(2f10.5,2x))') "rotmag",irot,m1,m2,(rotmag(mu),mu=1,3)
 !     end do ! irot

!       ==  Normalizes sum
 !     sumrho=sumrho/real(cryst_struc%nsym,kind=dp)
!        sumrho=glocnm(isppol,1,iatom,m1,m2) ! test without sym
 !     glocnms(iatom)%mat(m1,m2,isppol,1,1)=sumrho
 !     do mu=1,3
 !      rotmag(mu)=rotmag(mu)/real(cryst_struc%nsym,kind=dp)
!          rotmag(mu)=glocnm(isppol,mu+1,iatom,m1,m2) ! test without sym
 !      glocnms(iatom)%mat(m1,m2,isppol,mu+1,1)=rotmag(mu)
 !     enddo
 !    end do  ! m2
 !   end do ! m1
 !  endif ! lpawu/=-1
 ! end do ! iatom

!==  Compute back density matrix in upup dndn updn dnup representation
 ! call chg_repr_matlu(glocsym,glocnms,cryst_struc%natom,option=-1,prtopt=1)

!==  Put glocsym into gloc
!  do iatom=1,cryst_struc%natom
!    if(gloc(iatom)%lpawu/=-1) then
!      gloc(iatom)%mat=glocsym(iatom)%mat
!      gloc(iatom)%mat(:,:,1,:,:)=(glocsym(iatom)%mat(:,:,1,:,:) &
!&      + glocsym(iatom)%mat(:,:,2,:,:))/two
!      gloc(iatom)%mat(:,:,2,:,:)= gloc(iatom)%mat(:,:,1,:,:)
!      write(std_out,*) "WARNING: SYM non mag"
!      write(ab_out,*) "WARNING: SYM non mag"
!    endif
!  end do ! iatom

!  ABI_FREE(symrec_cart)
! endif

! call destroy_matlu(glocnm,cryst_struc%natom)
! call destroy_matlu(glocnms,cryst_struc%natom)
! call destroy_matlu(glocsym,cryst_struc%natom)
! ABI_FREE(glocnm)
! ABI_FREE(glocnms)
! ABI_FREE(glocsym)

 end subroutine sym_matlu
!!***

!!****f* m_matlu/inverse_matlu
!! NAME
!! inverse_matlu
!!
!! FUNCTION
!! Inverse local quantity.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom) :: input quantity to inverse
!!  natom=number of atoms in cell.
!!
!! OUTPUT
!!  matlu(natom) :: inverse of input matrix
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine inverse_matlu(matlu,natom)

!Arguments ------------------------------------
 integer, intent(in) :: natom
 type(matlu_type), intent(inout) :: matlu(natom)
!Local variables-------------------------------
 integer :: iatom,isppol,lpawu,ndim,nspinor,nsppol
 !************************************************************************

 nspinor = matlu(1)%nspinor
 nsppol  = matlu(1)%nsppol

 do iatom=1,natom
   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim = nspinor * (2*lpawu+1)
   do isppol=1,nsppol
     call xginv(matlu(iatom)%mat(:,:,isppol),ndim)
   end do ! isppol
 end do ! iatom

 !if(prtopt>0) then
 !endif
 !ABI_MALLOC(gathermatlu,(natom))
 !do iatom=1,natom
 !  if(matlu(iatom)%lpawu.ne.-1) then
 !    tndim=nsppol*nspinor*(2*matlu(iatom)%lpawu+1)
 !    ABI_MALLOC(gathermatlu(iatom)%value,(tndim,tndim))
 !    gathermatlu(iatom)%value=czero
 !  endif
 !enddo

 !call gather_matlu(matlu,gathermatlu,natom,option=1,prtopt=1)
 !do iatom=1,natom
 !  if(matlu(iatom)%lpawu.ne.-1) then
 !    tndim=nsppol*nspinor*(2*matlu(iatom)%lpawu+1)
     !call matcginv_dpc(gathermatlu(iatom)%value,tndim,tndim)
 !    call xginv(gathermatlu(iatom)%value,tndim)
 !  endif
 !enddo
 !call gather_matlu(matlu,gathermatlu,natom,option=-1,prtopt=1)

 !do iatom=1,natom
 !  if(matlu(iatom)%lpawu.ne.-1) then
 !    ABI_FREE(gathermatlu(iatom)%value)
 !  endif
 !enddo
 !ABI_FREE(gathermatlu)

 end subroutine inverse_matlu
!!***

!!****f* m_matlu/diff_matlu
!! NAME
!! diff_matlu
!!
!! FUNCTION
!!
!! INPUTS
!!  char1 = character describing matlu1
!!  char2 = character describing matlu2
!!  matlu1(natom) <type(matlu_type)>= density matrix 1 in the local orbital basis and related variables
!!  matlu2(natom) <type(matlu_type)>= density matrix 2 in the local orbital basis and related variables
!!  natom = number of atoms
!!  option =1      if diff > toldiff , stop
!!          0      print diff and toldiff
!!          else   do not test and do not print
!!  toldiff = threshold for the difference between matlu1 and matlu2
!!  zero_or_one = useful when comparing the overlap of Wannier functions at one kpt with the identity
!!
!! OUTPUT
!!  ierr = 0 if diff < toldiff
!!       = -1 otherwise
!!
!! SOURCE

subroutine diff_matlu(char1,char2,matlu1,matlu2,natom,option,toldiff,ierr,zero_or_one)

!Arguments ------------------------------------
 integer, intent(in) :: natom,option
 type(matlu_type), intent(in) :: matlu1(natom),matlu2(natom)
 character(len=*), intent(in) :: char1,char2
 real(dp), intent(in) :: toldiff
 integer, optional, intent(out) :: ierr
 integer, optional, intent(in) :: zero_or_one
!Local variables-------------------------------
 integer  :: iatom,idiff,lpawu,nspinor,nsppol
 real(dp) :: matludiff
 character(len=500) :: message
! *********************************************************************

 if (option /= 1 .and. option /= 0) return

 idiff     = 0
 matludiff = zero
 nspinor   = matlu1(1)%nspinor
 nsppol    = matlu1(1)%nsppol

 do iatom=1,natom
   lpawu = matlu1(iatom)%lpawu
   if (lpawu == -1) cycle
   matludiff = matludiff + sum(abs(matlu1(iatom)%mat(:,:,:)-matlu2(iatom)%mat(:,:,:)))
   idiff = idiff + (2*lpawu+1)**2
 end do ! iatom
 idiff = idiff * (nspinor**2) * nsppol

 if (.not. present(zero_or_one)) matludiff = matludiff / dble(idiff)

 if (matludiff < toldiff) then
   write(message,'(5a,6x,3a,4x,e12.4,a,e12.4)') ch10,'   ** Differences between ',trim(char1),' and ',&
     & ch10,trim(char2),' are small enough:',ch10,matludiff,' is lower than',toldiff
   call wrtout(std_out,message,'COLL')
   if (present(ierr)) ierr = 0
 else
   write(message,'(5a,3x,3a,3x,e12.4,a,e12.4)') ch10,'Differences between ',trim(char1),' and ',&
     & ch10,trim(char2),' is too large:',ch10,matludiff,' is larger than',toldiff
   ABI_WARNING(message)
!    write(message,'(8a,4x,e12.4,a,e12.4)') ch10,"  Matrix for ",trim(char1)
   write(message,'(a,3x,a)') ch10,trim(char1)
   call wrtout(std_out,message,'COLL')
   call print_matlu(matlu1(:),natom,1,opt_diag=-1)
   write(message,'(a,3x,a)') ch10,trim(char2)
   call wrtout(std_out,message,'COLL')
   call print_matlu(matlu2(:),natom,1,opt_diag=-1)
   if (present(zero_or_one) .and. (mod(matludiff,one) < toldiff)) then
     write(message,'(a,3x,a)') ch10," The norm is not identity for this k-point &
        & but is compatible with a high symmetry point"
     call wrtout(std_out,message,'COLL')
   else if (present(zero_or_one)) then
     write(message,'(a,3x,a)') ch10," The norm is not identity for this k-point but might be compatible &
       & with a high symmetry point: it should be checked"
     call wrtout(std_out,message,'COLL')
   else if (option == 1) then
     call flush_unit(std_out)
     write(message,'(5a,6x,3a,4x,e12.4,a,e12.4)') ch10,'   ** Differences between ',trim(char1),' and ',&
       & ch10,trim(char2),' are too high:',ch10,matludiff,' is greater than',toldiff
     ABI_ERROR(message)
   end if ! zero_or_one
   if (present(ierr)) ierr = -1
 end if ! matludiff < toldiff

end subroutine diff_matlu
!!***

!!****f* m_matlu/add_matlu
!! NAME
!! add_matlu
!!
!! FUNCTION
!!
!! INPUTS
!!  matlu1 <type(matlu_type)>= density matrix matlu1 in the local orbital basis and related variables
!!  matlu2 <type(matlu_type)>= density matrix matlu2 in the local orbital basis and related variables
!!  natom = number of atoms
!!  sign_matlu2= 1 add matlu1 and matlu2
!!              -1 substract matlu2 to matlu1
!!
!! OUTPUT
!!  matlu3 <type(matlu_type)>= density matrix matlu3, sum/substract matlu1 and matlu2
!!
!! SOURCE

subroutine add_matlu(matlu1,matlu2,matlu3,natom,sign_matlu2,idat,ndat)

!Arguments ------------------------------------
 integer, intent(in) :: natom,sign_matlu2
 integer, optional,intent(in) :: idat,ndat
 type(matlu_type), intent(in) :: matlu1(natom),matlu2(natom)
 type(matlu_type), intent(inout) :: matlu3(natom) !vz_i
!Local variables-------------------------------
 integer :: iatom,lpawu,isppol
! *********************************************************************

 if(present(idat) .and. present(ndat)) then
   do iatom=1,natom
     lpawu = matlu1(iatom)%lpawu
     if (lpawu == -1) cycle
     do isppol=1,matlu1(iatom)%nsppol
       if (sign_matlu2 == 1) then
         matlu3(iatom)%mat(:,:,(isppol-1)*ndat+idat) = matlu1(iatom)%mat(:,:,isppol) + matlu2(iatom)%mat(:,:,isppol)
       else if (sign_matlu2 == -1) then
         matlu3(iatom)%mat(:,:,(isppol-1)*ndat+idat) = matlu1(iatom)%mat(:,:,isppol) - matlu2(iatom)%mat(:,:,isppol)
       end if
     end do ! isppol
   end do ! iatom
 else
   do iatom=1,natom
     lpawu = matlu1(iatom)%lpawu
     if (lpawu == -1) cycle
     if (sign_matlu2 == 1) then
       matlu3(iatom)%mat(:,:,:) = matlu1(iatom)%mat(:,:,:) + matlu2(iatom)%mat(:,:,:)
     else if (sign_matlu2 == -1) then
       matlu3(iatom)%mat(:,:,:) = matlu1(iatom)%mat(:,:,:) - matlu2(iatom)%mat(:,:,:)
     end if
   end do ! iatom
 end if

end subroutine add_matlu
!!***

!!****f* m_matlu/chg_repr_matlu
!! NAME
!! chg_repr_matlu
!!
!! FUNCTION
!! Change representation of density matrix (useful for nspinor=2)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  glocspsp(natom) :: density matrix in the spin spin representation
!!  glocnm(natom) :: density matrix in the magnetization representation
!!  natom=number of atoms in cell.
!!  option=  1 glocspsp is input, glocnm is computed
!!        = -1 glocspsp is computed, glocnm is input
!!  prtopt= abs(prtopt) >= 3 : print in magnetization representation
!!
!! OUTPUT
!!  glocspsp(natom) :: density matrix in the spin spin representation
!!  glocnm(natom) :: density matrix in the magnetization representation
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine chg_repr_matlu(glocspsp,glocnm,natom,option,prtopt)

!Arguments ------------------------------------
 integer, intent(in) :: natom,option,prtopt
 type(matlu_type), intent(inout), target :: glocnm(natom),glocspsp(natom)
!Local variables-------------------------------
 integer :: iatom,lpawu,m1,m2,mu,ndim,nsppol,isppol,gpu_option
 complex(dp), ABI_CONTIGUOUS pointer :: glocnm_mat(:,:,:),glocspsp_mat(:,:,:)
 character(len=500) :: message

! DBG_ENTER("COLL")
 ABI_CHECK(glocspsp(1)%nsppol*4==glocnm(1)%nsppol, "Mismatch in nsppol between glocspsp and glocnm")
 ABI_CHECK(glocspsp(1)%gpu_option==glocnm(1)%gpu_option, "Mismatch in gpu_option between glocspsp and glocnm")
 nsppol=glocspsp(1)%nsppol
 gpu_option=glocspsp(1)%gpu_option

!==  Compute density matrix in density magnetization representation
 if (option == 1) then
   do iatom=1,natom
     lpawu = glocspsp(iatom)%lpawu
     if (lpawu == -1) cycle
     ndim = 2*lpawu + 1
     glocnm_mat => glocnm(iatom)%mat
     glocspsp_mat => glocspsp(iatom)%mat

#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET TEAMS DISTRIBUTE PRIVATE(isppol) MAP(to:glocnm,glocspsp) &
     !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
     do isppol=1,nsppol
       !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(m1,m2)
       do m2=1,ndim
         do m1=1,ndim
           glocnm_mat(m1,m2,1+(isppol-1)*4) = glocspsp_mat(m1,m2,     isppol)  + glocspsp_mat(m1+ndim,m2+ndim,isppol)
           glocnm_mat(m1,m2,4+(isppol-1)*4) = glocspsp_mat(m1,m2,     isppol)  - glocspsp_mat(m1+ndim,m2+ndim,isppol)
           glocnm_mat(m1,m2,2+(isppol-1)*4) = glocspsp_mat(m1,m2+ndim,isppol)  + glocspsp_mat(m1+ndim,m2,isppol)
           glocnm_mat(m1,m2,3+(isppol-1)*4) = (glocspsp_mat(m1,m2+ndim,isppol) - glocspsp_mat(m1+ndim,m2,isppol)) * j_dpc
         end do  ! m1
       end do ! m2
     end do ! isppol

     if (abs(prtopt) >= 3) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET UPDATE FROM(glocnm) IF(gpu_option==ABI_GPU_OPENMP)
#endif
       write(message,'(a)') "        -- in n, m repr "
       call wrtout(std_out,message,'COLL')
       do mu=1,4
         do m1=1,ndim
           write(message,'(8x,(14(2f9.5,2x)))') (glocnm(iatom)%mat(m1,m2,isppol+mu-1),m2=1,ndim)
           call wrtout(std_out,message,'COLL')
         end do ! m1
         write(message,'(a)') ch10
         call wrtout(std_out,message,'COLL')
       end do ! mu
     end if ! prtopt >3
   end do ! iatom

!==  Compute back density matrix in upup dndn updn dnup representation
 else if (option == -1) then

   do iatom=1,natom
     lpawu = glocnm(iatom)%lpawu
     if (lpawu == -1) cycle
     ndim = 2*lpawu + 1
     glocnm_mat => glocnm(iatom)%mat
     glocspsp_mat => glocspsp(iatom)%mat

#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET TEAMS DISTRIBUTE PRIVATE(isppol) MAP(to:glocnm,glocspsp) &
     !$OMP& IF(gpu_option==ABI_GPU_OPENMP)
#endif
     do isppol=1,nsppol
       !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(m1,m2)
       do m2=1,ndim
         do m1=1,ndim
           glocspsp_mat(m1,m2,isppol) = &
           &    half * (glocnm_mat(m1,m2,1+(isppol-1)*4)+glocnm_mat(m1,m2,4+(isppol-1)*4))
           glocspsp_mat(m1+ndim,m2+ndim,isppol) = &
           &    half * (glocnm_mat(m1,m2,1+(isppol-1)*4)-glocnm_mat(m1,m2,4+(isppol-1)*4))
           glocspsp_mat(m1,m2+ndim,isppol) = &
           &    half * (glocnm_mat(m1,m2,2+(isppol-1)*4)-j_dpc*glocnm_mat(m1,m2,3+(isppol-1)*4))
           glocspsp_mat(m1+ndim,m2,isppol) = &
           &    half * (glocnm_mat(m1,m2,2+(isppol-1)*4)+j_dpc*glocnm_mat(m1,m2,3+(isppol-1)*4))
         end do ! m1
       end do ! m2
     end do ! isppol
     !if (abs(prtopt) > 6) then
     !  write(message,'(a)') "        -- in spin spin repr "
     !  call wrtout(std_out,message,'COLL')
     !  do mu=1,4
     !    do m1=1,ndim
     !      write(message,'(8x,14(2f9.5,2x))') (glocspsp(iatom)%mat(m1,m2,isppol,mu,1),m2=1,ndim)
     !      call wrtout(std_out,  message,'COLL')
     !    end do ! m1
     !    write(message,'(a)') ch10
     !    call wrtout(std_out,message,'COLL')
     !  end do
     !end if ! prtopt>6
   end do ! iatom
 else
   message = "stop in chg_repr_matlu"
   ABI_ERROR(message)
 end if ! option

! DBG_EXIT("COLL")

 end subroutine chg_repr_matlu
!!***

!!****f* m_matlu/trace_matlu
!! NAME
!! trace_matlu
!!
!! FUNCTION
!!  Compute the trace of the matlu matrix
!!
!! INPUTS
!!  matlu(natom) <type(matlu_type)>= density matrix in the
!!               local orbital basis and related variables
!!  natom = number of atoms
!!  itau = flag for print
!!      not present (default) : occupations from G(iw)
!!       = 1  : occupations from G(tau)
!!       = -1 : occupations from G0(tau)
!!       = 4  : trace of matlu
!!
!! OUTPUT
!!  trace_loc(nsppol+1,natom)= trace for each atom and each polarization,
!!                             trace_loc(iatom,nsppol+1) is
!!                             the full trace over all polarizations
!!  trace= trace over all correlated atoms
!!
!! SOURCE

 subroutine trace_matlu(matlu,natom,trace_loc,itau,trace)

!Arguments ------------------------------------
 integer, intent(in) :: natom
 type(matlu_type), intent(in) :: matlu(natom)
 real(dp), target, optional, intent(inout) :: trace_loc(matlu(1)%nsppol+1,natom)
 integer, optional, intent(in) :: itau
 complex(dp), optional, intent(out) :: trace
!Local variables-------------------------------
 integer :: iatom,im,isppol,lpawu,ndim,nspinor,nsppol
 complex(dp) :: trace_tmp,trace_tmp2
 real(dp), ABI_CONTIGUOUS pointer :: traceloc(:,:) => null()
 character(len=4) :: tag
 character(len=12) :: tag_nb_elec
 character(len=500) :: message
! *********************************************************************

 nspinor = matlu(1)%nspinor
 nsppol  = matlu(1)%nsppol

 if (present(trace_loc)) then
   traceloc => trace_loc(:,:)
 else
   ABI_MALLOC(traceloc,(nsppol+1,natom))
 end if

 trace_tmp = czero
 traceloc(:,:) = zero

 do iatom=1,natom

   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle

   ndim = nspinor * (2*lpawu+1)
   write(tag,'(i4)') iatom
   write(message,'(3a)') ch10,'   -------> For Correlated Atom ',adjustl(tag)
   if (.not. present(itau)) then
     call wrtout(std_out,message,'COLL')
   end if

   if (present(itau)) then
     if (itau > 0) then
       call wrtout(std_out,message,'COLL')
     end if
   end if ! present(itau)

   do isppol=1,nsppol
     trace_tmp2 = czero
     do im=1,ndim
       trace_tmp2 = trace_tmp2 + matlu(iatom)%mat(im,im,isppol)
     end do ! im
     trace_tmp = trace_tmp + trace_tmp2
     traceloc(isppol,iatom) = dble(trace_tmp2)
     traceloc(nsppol+1,iatom) = traceloc(nsppol+1,iatom) + traceloc(isppol,iatom)
   end do ! isppol
   if (nsppol == 1 .and. nspinor == 1) traceloc(nsppol+1,iatom) = traceloc(nsppol+1,iatom) * two
   write(tag_nb_elec,'(f12.6)') traceloc(nsppol+1,iatom)
   tag_nb_elec = adjustl(tag_nb_elec)
   if (.not. present(itau)) then
     write(message,'(8x,2a)') 'Nb of Corr. elec. from G(iw) is: ',tag_nb_elec
     call wrtout(std_out,message,'COLL')
   end if ! not present(itau)
   if (present(itau)) then
     if (itau == 1) then
       write(message,'(8x,2a)') 'Nb of Corr. elec. from G(tau=0-) is: ',tag_nb_elec
       call wrtout(std_out,message,'COLL')
     else if (itau == -1) then
       write(message,'(8x,2a)') 'Nb: Sum of the values of G0(tau=0-) is: ',tag_nb_elec
       call wrtout(std_out,message,'COLL')
     else if (itau == 4) then
       write(message,'(8x,2a)') 'Trace of matlu matrix is: ',tag_nb_elec
       call wrtout(std_out,message,'COLL')
     end if ! itau
   end if ! present(itau)
 end do ! iatom

 if (present(trace)) then
   if (nsppol == 1 .and. nspinor == 1) trace_tmp = trace_tmp * two
   trace = trace_tmp
 end if ! present(trace)

 if (nsppol > 1 .and. (.not. present(trace))) then
   do iatom=1,natom
     lpawu = matlu(iatom)%lpawu
     if (lpawu == -1) cycle
     !! MAG
    ! if(nsppol>1.and.present(itau)) then
    !   if(itau==1) then
     write(tag_nb_elec,'(f12.6)') traceloc(2,iatom) - traceloc(1,iatom)
     write(message,'(8x,2a)') 'DMFT Corr. Elec. Mag.: ',adjustl(tag_nb_elec)
     call wrtout(std_out,message,'COLL')
    !   endif
   end do ! iatom
 end if ! nsppol>1

 if (.not. present(trace_loc)) then
   ABI_FREE(traceloc)
 end if
 traceloc => null()

 end subroutine trace_matlu
!!***

!!****f* m_matlu/gather_matlu
!! NAME
!! gather_matlu
!!
!! FUNCTION
!! Create new array from matlu
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gloc(natom) <type(matlu_type)>        = density matrix in the spin spin representation
!!  gatherloc(natom) <type(coeff2c_type)> = density matrix where spin and angular momentum are gathered in the same index
!!  natom=number of atoms in cell.
!!  option= 1 go from gloc to gathergloc
!!  option= -1 go from gathergloc to gloc
!!  prtopt= option to define level of printing
!!
!! OUTPUT
!!  gloc(natom) <type(matlu_type)>        = density matrix in the spin spin representation
!!  gatherloc(natom) <type(coeff2c_type)> = density matrix where spin and angular momentum are gathered in the same index
!!
!! SIDE EFFECTS
!!
!! SOURCE
 subroutine gather_matlu(gloc,gathergloc,natom,option,prtopt)

 use defs_wvltypes
 use m_crystal, only : crystal_t

! type  matlus_type
!  SEQUENCE
!  complex(dp), pointer :: mat(:,:)
! end type matlus_type

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,option,prtopt
 type(coeff2c_type), intent(inout) :: gathergloc(natom)
 type(matlu_type),intent(inout) :: gloc(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im1,im2,ispinor,ispinor1,isppol,isppol1
 integer :: jc1,jc2,ml1,ml2,ndim,nspinor,nsppol,tndim
 character(len=500) :: message

! DBG_ENTER("COLL")
 nsppol=gloc(1)%nsppol
 nspinor=gloc(1)%nspinor

 do iatom=1,natom
   if(gloc(iatom)%lpawu.ne.-1) then
!==-------------------------------------

     ndim=2*gloc(iatom)%lpawu+1
     tndim=nsppol*nspinor*ndim

!== Put norm into array "gathergloc"
     jc1=0
     do isppol=1,nsppol
       do ispinor=1,nspinor
         do ml1=1,ndim
           jc1=jc1+1
           jc2=0
           do isppol1=1,nsppol
             do ispinor1=1,nspinor
               do ml2=1,ndim
                 jc2=jc2+1
                 if(option==1) then
                   if(isppol==isppol1) then
                     gathergloc(iatom)%value(jc1,jc2)=gloc(iatom)%mat(ml1+(ispinor-1)*ndim,ml2+(ispinor1-1)*ndim,isppol)
                   endif
                 else if(option==-1) then
                   if(isppol==isppol1) then
                     gloc(iatom)%mat(ml1+(ispinor-1)*ndim,ml2+(ispinor1-1)*ndim,isppol)=gathergloc(iatom)%value(jc1,jc2)
                   endif
                 endif
               enddo
             enddo ! ispinor1
           enddo ! isppol1
         enddo
       enddo !ispinor
     enddo ! isppol
   endif
 enddo ! iatom
 if(option==1.and.prtopt==3) then
   do iatom=1,natom
     if(gloc(iatom)%lpawu.ne.-1) then
       tndim=nsppol*nspinor*(2*gloc(iatom)%lpawu+1)
       write(message,'(2a,i5)') ch10,' (gathermatlu:) For atom', iatom
       call wrtout(std_out,message,'COLL')
       do im1=1,tndim
         write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
&         (gathergloc(iatom)%value(im1,im2),im2=1,tndim)
         call wrtout(std_out,message,'COLL')
       end do
     endif
   enddo ! iatom
 else if(option==-1.and.prtopt==3) then
   call print_matlu(gloc,natom,prtopt)
 endif



! DBG_EXIT("COLL")

 end subroutine gather_matlu
!!***

!!****f* m_matlu/diag_matlu
!! NAME
!! diag_matlu
!!
!! FUNCTION
!! Diagonalize hermitian matlu matrix
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom) :: input quantity to diagonalize (careful, matlu must be hermitian)
!!  natom=number of atoms
!!  prtopt: >=3 : print eigenvectors
!!          >=4 : print matlu before diagonalization
!!  nsppol_imp= if 1, one can diagonalize with the same matrix the Up
!!   and Dn matlu matrix. It is convenient because one can thus have the
!!   same interaction matrix for up and dn spins. Default is nsppol.
!!  checkstop= if true (default), print the matrix for spin down in the diagonalization basis of spin up
!!             (useful when nsppol=2 and nsppol_imp=1)
!!  optreal= diagonalize the real matrix if max(imag(matlu)) < 1e-6
!!  test= if 8, use the block diagonalization algorithm (only when the real matrix is diagonalized)
!!
!! OUTPUT
!!  matlu_diag(natom) :: diagonalized density matrix
!!  eigvectmatlu(natom) = Eigenvectors corresponding to the diagonalization
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine diag_matlu(matlu,matlu_diag,natom,prtopt,eigvectmatlu,nsppol_imp,checkstop,opt_real,test)

!Arguments ------------------------------------
 integer, intent(in) :: natom,prtopt
 type(matlu_type), intent(in) :: matlu(natom)
 type(matlu_type), intent(inout) :: eigvectmatlu(natom),matlu_diag(natom) !vz_i
 integer, optional, intent(in) :: nsppol_imp,opt_real,test
 logical, optional, intent(in) :: checkstop
!Local variables-------------------------------
 integer :: iatom,im1,im2,info,isppol,lpawu,lwork,lworkr
 integer :: nspinor,nsppol,nsppolimp,optreal,tndim
 logical :: blockdiag,checkstop_in,print_temp_mat2
 character(len=4) :: tag
 character(len=500) :: message
 real(dp), allocatable :: eig(:),rwork(:),valuer(:,:),work(:)!,valuer2(:,:)
 !real(dp),allocatable :: valuer3(:,:),valuer4(:,:)
! real(dp),allocatable :: eigvec(:,:)
 complex(dp), allocatable :: temp_mat(:,:),zwork(:)
!debug complex(dp),allocatable :: temp_mat3(:,:)
!************************************************************************

 blockdiag    = .false.
 checkstop_in = .true.
 nspinor      = matlu(1)%nspinor
 nsppol       = matlu(1)%nsppol
 nsppolimp    = nsppol
 optreal      = 0

 if (present(nsppol_imp)) nsppolimp = nsppol_imp
 if (present(checkstop)) checkstop_in = checkstop
 if (present(test)) blockdiag = (test == 8)
 if (present(opt_real)) optreal = opt_real

 call zero_matlu(matlu_diag(:),natom)
 call copy_matlu(matlu(:),eigvectmatlu(:),natom)

 !donotdiag=.true.
 !donotdiag=.false.
! ===========================
! Check is diagonalization is necessary and how
! ===========================
 !do isppol=1,matlu(1)%nsppol
 !  do iatom=1,natom
 !     if(matlu(iatom)%lpawu.ne.-1) then
 !      tndim=(2*matlu(iatom)%lpawu+1)
 !       do im1=1,tndim
 !         do im2=1,tndim
 !           do ispinor=1,nspinor
 !             do ispinor1=1,nspinor
 !               if(abs(matlu(iatom)%mat(im1,im2,isppol,ispinor,ispinor1))>tol8.and.&
!&                 (im1/=im2.or.ispinor/=ispinor1)) then
!                 if matrix is diagonal: do not diagonalize
 !                 donotdiag=.false.
 !                 exit
 !               endif
 !             enddo
 !           enddo
 !         enddo
 !       enddo
 !     endif
 !  enddo
 !enddo

 !if(donotdiag) then
 !  do isppol=1,matlu(1)%nsppol
 !    do iatom=1,natom
 !       if(matlu(iatom)%lpawu.ne.-1) then
 !         tndim=(2*matlu(iatom)%lpawu+1)
 !         eigvectmatlu(iatom,isppol)%value(:,:)=czero
 !         do im1=1,tndim
 !           eigvectmatlu(iatom,isppol)%value(im1,im1)=cone
 !         enddo
 !       endif
 !    enddo
 !  enddo
 !  call copy_matlu(matlu,matlu_diag,natom)
 !  write(message,'(a)')  "   Diagonalisation of matlu will not be performed"
 !  call wrtout(std_out,message,'COLL')
 !  return
 !endif

! For nsppol=2, and if nsppolimp=1, the eigenvectors are computed for isppol=1, and applied through
! rotate_matlu to isppol=2. It is the reason why the sum below is only from 1 to nsppolimp !

 do iatom=1,natom

   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   write(tag,'(i4)') iatom
   tndim = nspinor * (2*lpawu+1)
   lwork = 2*tndim - 1
   lworkr = tndim * (tndim+2) * 2
   ABI_MALLOC(eig,(tndim))

! ===========================
! Define gathermatlu
! ===========================
   !ABI_MALLOC(gathermatlu,(natom))
   !do iatom=1,natom
   !  if(matlu(iatom)%lpawu.ne.-1) then
   !    tndim=nspinor*(2*matlu(iatom)%lpawu+1)
   !    ABI_MALLOC(gathermatlu(iatom)%value,(tndim,tndim))
   !    gathermatlu(iatom)%value=czero
   !  endif
   !enddo
   !if(nsppol==1.and.nspinor==2) then
   !  call gather_matlu(matlu,gathermatlu,natom,option=1,prtopt=1)
   !else if((nsppol==2.or.nsppol==1).and.nspinor==1) then
   !  do iatom=1,natom
   !    if(matlu(iatom)%lpawu.ne.-1) then
   !      tndim=nspinor*(2*matlu(iatom)%lpawu+1)
   !      do im1=1,tndim
   !        do im2=1,tndim
   !          gathermatlu(iatom)%value(im1,im2)=matlu(iatom)%mat(im1,im2,isppol,1,1)
   !        enddo
   !      enddo
   !    endif
   !  enddo
   !endif

   ! ===========================
   ! Diagonalize
   ! ===========================
   do isppol=1,nsppolimp
!debug       allocate(temp_mat2(tndim,tndim))
!debug       temp_mat2=zero
!         ABI_MALLOC(valuer2,(tndim,tndim))
!         ABI_MALLOC(valuer3,(tndim,tndim))
!         ABI_MALLOC(valuer4,(tndim,tndim))
!         valuer2=zero
!         valuer3=zero
!         valuer4=zero
     if (prtopt >= 4) then
       write(message,'(a,i4,a,i4)') "       BEFORE DIAGONALIZATION for atom",iatom,"  and isppol",isppol
       call wrtout(std_out,message,'COLL')
       do im1=1,tndim
         write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))') (matlu(iatom)%mat(im1,im2,isppol),im2=1,tndim)
         call wrtout(std_out,message,'COLL')
       end do ! im1
     end if ! prtopt>=4
!debug       temp_mat2(:,:)=gathermatlu(iatom)%value(:,:)
!           write(std_out,*)"diag"

     if (optreal == 1 .and. maxval(abs(aimag(matlu(iatom)%mat(:,:,isppol)))) < tol6) then
       write(message,'(a,2x,a,e9.3,a)') ch10,"Imaginary part of Local Hamiltonian is lower than ",&
         & tol6,": the real matrix is used"
       call wrtout(std_out,message,'COLL')
       ABI_MALLOC(valuer,(tndim,tndim))
       valuer(:,:) = dble(matlu(iatom)%mat(:,:,isppol))
!           write(message,'(a)') ch10
!           call wrtout(std_out,message,'COLL')
!           write(message,'(a,i4,a,i4)')  "BEFORE valuer for atom",iatom,"  and isppol",isppol
!           call wrtout(std_out,message,'COLL')
!           do im1=1,tndim
!             write(message,'(2(1x,18(1x,"(",f20.15,",",f20.15,")")))')&
!&             (valuer(im1,im2),im2=1,tndim)
!             call wrtout(std_out,message,'COLL')
!           end do
!           do im1=1,tndim
!             valuer(im1,im1)=real(im1,kind=dp)*0.00000000001_dp+valuer(im1,im1)
!           enddo
!           write(message,'(a)') ch10
!           call wrtout(std_out,message,'COLL')
!           write(message,'(a,i4,a,i4)')  "BEFORE valuer for atom",iatom,"  and isppol",isppol
!           call wrtout(std_out,message,'COLL')
!           do im1=1,tndim
!             write(message,'(2(1x,18(1x,f20.15,f20.15)))')&
!&             (valuer(im1,im2),im2=1,tndim)
!             call wrtout(std_out,message,'COLL')
!           end do
           !call dsyev('v','u',tndim,valuer,tndim,eig,work,lworkr,info)
       if (blockdiag) then
         call blockdiago_fordsyev(valuer(:,:),tndim,eig(:))
       else
         ABI_MALLOC(work,(lworkr))
         call dsyev('v','u',tndim,valuer(:,:),tndim,eig(:),work(:),lworkr,info)
         ABI_FREE(work)
       end if ! blockdiag
!!       For reproductibility
!           ! valuer2: eigenvector for the perturb matrix
!           valuer2=real(gathermatlu(iatom)%value,kind=dp)
!           do im1=1,tndim
!             valuer2(im1,im1)=float(im1)*0.00000000001+valuer2(im1,im1)
!           enddo
!           call dsyev('v','u',tndim,valuer2,tndim,eig,work,lworkr,info)
!           write(message,'(a)') ch10
!           call wrtout(std_out,message,'COLL')
!           write(message,'(a,i4,a,i4)')  "       valuer2 for atom",iatom,"  and isppol",isppol
!           call wrtout(std_out,message,'COLL')
!           do im1=1,tndim
!             write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
!&             (valuer2(im1,im2),im2=1,tndim)
!             call wrtout(std_out,message,'COLL')
!           end do
!           call dgemm('n','n',tndim,tndim,tndim,cone,valuer,tndim,&
!&            valuer2,tndim,czero,valuer3                ,tndim)
!           call dgemm('c','n',tndim,tndim,tndim,cone,valuer2,tndim,&
!&            valuer3                   ,tndim,czero,valuer4,tndim)
!           ! valuer4: compute unpert matrix in the basis of the
!           ! perturb basis
!           write(message,'(a)') ch10
!           call wrtout(std_out,message,'COLL')
!           write(message,'(a,i4,a,i4)')  "BEFORE valuer4 for atom",iatom,"  and isppol",isppol
!           call wrtout(std_out,message,'COLL')
!           do im1=1,tndim
!             write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
!&             (valuer4(im1,im2),im2=1,tndim)
!             call wrtout(std_out,message,'COLL')
!           end do
!           call dsyev('v','u',tndim,valuer4,tndim,eig,work,lworkr,info)
!           ! valuer4: Diago valuer4 (nearly diag)
!           write(message,'(a)') ch10
!           call wrtout(std_out,message,'COLL')
!           write(message,'(a,i4,a,i4)')  "AFTER  valuer4 for atom",iatom,"  and isppol",isppol
!           call wrtout(std_out,message,'COLL')
!           do im1=1,tndim
!             write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
!&             (valuer4(im1,im2),im2=1,tndim)
!             call wrtout(std_out,message,'COLL')
!           end do
!           call dgemm('n','n',tndim,tndim,tndim,cone,valuer2,tndim,&
!&            valuer4,tndim,czero,valuer                ,tndim)
           !write(6,*) "INFO",info
       eigvectmatlu(iatom)%mat(:,:,isppol) = cmplx(valuer(:,:),zero,kind=dp)
       ABI_FREE(valuer)
!           write(message,'(a,i4,a,i4)')  "AFTER valuer for atom",iatom,"  and isppol",isppol
!           call wrtout(std_out,message,'COLL')
!           do im1=1,tndim
!             write(message,'(2(1x,18(1x,"(",f20.15,",",f20.15,")")))')&
!&             (valuer(im1,im2),im2=1,tndim)
!             call wrtout(std_out,message,'COLL')
!           end do
     else
       if (optreal == 1 .and. maxval(abs(aimag(matlu(iatom)%mat(:,:,isppol)))) > tol8) then
         write(message,'(a)') " Local hamiltonian in correlated basis is complex"
         ABI_COMMENT(message)
       end if
       ABI_MALLOC(zwork,(lwork))
       ABI_MALLOC(rwork,(3*tndim-2))
       call zheev('v','u',tndim,eigvectmatlu(iatom)%mat(:,:,isppol),tndim,eig(:),zwork(:),lwork,rwork(:),info)
       ABI_FREE(zwork)
       ABI_FREE(rwork)
           !call blockdiago_forzheev(gathermatlu(iatom)%value,tndim,eig)
     end if ! present(optreal)
     if (prtopt >= 3) then
       write(message,'(a)') ch10
       call wrtout(std_out,message,'COLL')
       write(message,'(3a,i1)') "       EIGENVECTORS for atom ",trim(adjustl(tag))," and isppol ",isppol
       call wrtout(std_out,message,'COLL')
       do im1=1,tndim
         write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))') (eigvectmatlu(iatom)%mat(im1,im2,isppol),im2=1,tndim)
         call wrtout(std_out,message,'COLL')
       end do ! im1
          ! do im1=1,tndim
          !   xcheck=czero
          !   do im3=1,tndim
          !     do im2=1,tndim
          !       xcheck=xcheck+gathermatlu(iatom)%value(im1,im2)*conjg(gathermatlu(iatom)%value(im2,im3))
          !     end do
          !   end do
          !   write(6,*) "check",im3,im1,xcheck
          ! end do
     end if ! prtopt>=3
!       write(std_out,*) "eig",eig
! ===========================
! Put eigenvalue in matlu_diag
! ===========================
     do im1=1,tndim
       matlu_diag(iatom)%mat(im1,im1,isppol) = cmplx(eig(im1),zero,kind=dp)
     end do ! im1
         !if(prtopt>=2) then
!            write(message,'(a,12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
!&            ch10,(eig(im1),im1=1,tndim)
!             call wrtout(std_out,message,'COLL')
           !call wrtout(std_out,message,'COLL')
            !write(std_out,*) "EIG", eig
         !endif
!         ABI_FREE(valuer2)
!         ABI_FREE(valuer3)
!         ABI_FREE(valuer4)
!     endif
!   enddo
! ===========================
! Keep eigenvectors gathermatlu
! ===========================
         !if (present(eigvectmatlu)) then
         !  tndim=nspinor*(2*matlu(iatom)%lpawu+1)
         !  eigvectmatlu(iatom,isppol)%value(:,:)=gathermatlu(iatom)%value(:,:)
!           write(std_out,*) "eigvect in diag_matlu"
!           do im1=1,tndim
!             write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
!&             (gathermatlu(iatom)%value(im1,im2),im2=1,tndim)
!             call wrtout(std_out,message,'COLL')
!           end do
         !endif

     if (nsppolimp == 1 .and. nsppol == 2) then
! ==================================================================
! If necessary rotate levels for this other spin, assuming the same
! rotation matrix: it has to be checked afterwards that the matrix is
! diagonal
! ===================================================================

!        input matrix: gathermatlu
!        rotation matrix: eigvectmatlu
!        intermediate matrix: temp_mat
!        result matrix: temp_mat2
         !do im1=1,tndim
         !  do im2=1,tndim
         !    gathermatlu(iatom)%value(im1,im2)=matlu(iatom)%mat(im1,im2,2,1,1)
         !  enddo
         !enddo

       if (prtopt >= 3) then
         write(message,'(a,i4,a,i4)') "       MATLU for atom",iatom," inside if nsppolimp==1"
         call wrtout(std_out,message,'COLL')
         do im1=1,tndim
           write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))') (matlu(iatom)%mat(im1,im2,2),im2=1,tndim)
           call wrtout(std_out,message,'COLL')
         end do ! im1
       end if ! prtopt>=3

       ABI_MALLOC(temp_mat,(tndim,tndim))

       call abi_xgemm('n','n',tndim,tndim,tndim,cone,matlu(iatom)%mat(:,:,2),tndim,&
                    & eigvectmatlu(iatom)%mat(:,:,1),tndim,czero,temp_mat(:,:),tndim)

       call abi_xgemm('c','n',tndim,tndim,tndim,cone,eigvectmatlu(iatom)%mat(:,:,1),tndim,&
                    & temp_mat(:,:),tndim,czero,matlu_diag(iatom)%mat(:,:,2),tndim)

       ABI_FREE(temp_mat)

       eigvectmatlu(iatom)%mat(:,:,2) = eigvectmatlu(iatom)%mat(:,:,1)
       print_temp_mat2 = .false.
       do im2=1,tndim
         do im1=1,tndim
           if (im1 /= im2 .and. abs(matlu_diag(iatom)%mat(im1,im2,2)) > tol5) then
             write(message,'(3a,i4,2f16.4)') ch10,'diag_matlu= Matrix for spin number 2 obtained with', &
              & ' eigenvectors from diagonalization for spin nb 1 is non diagonal for atom:',iatom,&
              & abs(matlu_diag(iatom)%mat(im1,im2,2)),tol5
             call wrtout(std_out,message,'COLL')
             if (abs(matlu_diag(iatom)%mat(im1,im2,2)) > tol1 .or. checkstop_in) print_temp_mat2 = .true.
           end if
         end do ! im1
       end do ! im2

       if (print_temp_mat2 .and. prtopt >= 3) then
         write(message,'(a)') "       temp_mat2"
         call wrtout(std_out,message,'COLL')
         do im1=1,tndim
           write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))') (matlu_diag(iatom)%mat(im1,im2,2),im2=1,tndim)
           call wrtout(std_out,message,'COLL')
         end do ! im1
         if (iatom == 2) ABI_ERROR("iatom==2")
       end if ! print_temp_mat2
     end if ! nsppol_imp=1 and nsppol=2
   end do ! isppol

   ABI_FREE(eig)

!!  for check only
!debug     if(matlu(iatom)%lpawu.ne.-1) then
!debug       allocate(temp_mat(tndim,tndim))
!debug       allocate(temp_mat3(tndim,tndim))
!debug           do im1=1,tndim
!debug             do im2=1,tndim
!debug!               rot_mat(iatom,isppol)%value(im1,im2)=rot_mat_orig(iatom,isppol)%value(im1,im2)
!debug               temp_mat3(im1,im2)=conjg(gathermatlu(iatom)%value(im2,im1))
!debug             enddo
!debug           enddo
!debug       temp_mat(:,:)=czero
!debug!      input matrix: temp_mat2
!debug!      rotation matrix: gathermatlu
!debug!      intermediate matrix: temp_mat
!debug!      result matrix: temp_mat2
!debug       call zgemm('n','c',tndim,tndim,tndim,cone,temp_mat2   ,tndim,&
!debug&        temp_mat3,tndim,czero,temp_mat                ,tndim)
!debug       call zgemm('n','n',tndim,tndim,tndim,cone,temp_mat3,tndim,&
!debug&        temp_mat                   ,tndim,czero,temp_mat2,tndim)
!debug!       call zgemm('n','c',tndim,tndim,tndim,cone,temp_mat2   ,tndim,&
!debug!&        gathermatlu(iatom)%value,tndim,czero,temp_mat                ,tndim)
!debug!       call zgemm('n','n',tndim,tndim,tndim,cone,gathermatlu(iatom)%value,tndim,&
!debug!&        temp_mat                   ,tndim,czero,temp_mat2,tndim)
!debug         write(std_out,*) "result"
!debug         do im1=1,tndim
!debug           write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
!debug&           (temp_mat2(im1,im2),im2=1,tndim)
!debug           call wrtout(std_out,message,'COLL')
!debug         end do
!debug       deallocate(temp_mat)
!debug       deallocate(temp_mat3)
!debug     endif ! lpawu

 end do  ! iatom
! End loop over atoms
! ===========================

 end subroutine diag_matlu
!!***

!!****f* m_matlu/rotate_matlu
!! NAME
!! rotate_matlu
!!
!! FUNCTION
!! Rotate matlu matrix
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu_inp(natom) :: input quantity to rotate
!!  rot_mat(natom) :: Rotation matrix (usually from diag_matlu)
!!  natom=number of atoms in cell.
!!  inverse=   1: rot_mat^H * matlu * rot_mat (from original basis to diagonal basis)
!!           /=1: rot_mat * matlu * rot_mat^H (from diagonal basis to original basis)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine rotate_matlu(matlu_inp,rot_mat,natom,inverse)

!Arguments ------------------------------------
 integer, intent(in) :: inverse,natom
 type(matlu_type), intent(inout) :: matlu_inp(natom)
 type(matlu_type), intent(in) :: rot_mat(natom)
!Local variables-------------------------------
 integer :: iatom,isppol,lpawu,nspinor,nsppol,tndim
 complex(dp), allocatable :: temp_mat(:,:)
 character(len=1) :: c1,c2
!************************************************************************

 nspinor = matlu_inp(1)%nspinor
 nsppol  = matlu_inp(1)%nsppol

 if (inverse == 1) then
   c1 = "n" ; c2 = "c"
 else
   c1 = "c" ; c2 = "n"
 end if ! inverse

 do iatom=1,natom
   lpawu = matlu_inp(iatom)%lpawu
   if (lpawu == -1) cycle
   tndim = nspinor * (2*lpawu+1)
   ABI_MALLOC(temp_mat,(tndim,tndim))
   do isppol=1,nsppol
     call abi_xgemm('n',c1,tndim,tndim,tndim,cone,matlu_inp(iatom)%mat(:,:,isppol),tndim,&
                  & rot_mat(iatom)%mat(:,:,isppol),tndim,czero,temp_mat(:,:),tndim)
     call abi_xgemm(c2,'n',tndim,tndim,tndim,cone,rot_mat(iatom)%mat(:,:,isppol),tndim,&
                  & temp_mat(:,:),tndim,czero,matlu_inp(iatom)%mat(:,:,isppol),tndim)
   end do ! isppol
   ABI_FREE(temp_mat)
 end do ! iatom

 !do isppol=1,nsppol

! ===========================
! Define gathermatlu and rot_mat_orig and allocate
! ===========================
   !ABI_MALLOC(rot_mat_orig,(natom))
   !ABI_MALLOC(gathermatlu,(natom))
   !do iatom=1,natom
   !  if(matlu(iatom)%lpawu.ne.-1) then
   !    tndim=nspinor*(2*matlu(iatom)%lpawu+1)
   !    ABI_MALLOC(gathermatlu(iatom)%value,(tndim,tndim))
   !    gathermatlu(iatom)%value=czero
!       ABI_MALLOC(rot_mat_orig(iatom,isppol)%value,(tndim,tndim))
!       rot_mat_orig(iatom,isppol)%value(:,:)=rot_mat(iatom,isppol)%value(:,:)
   !    ABI_MALLOC(rot_mat_orig(iatom)%value,(tndim,tndim))
   !    rot_mat_orig(iatom)%value(:,:)=rot_mat(iatom,isppol)%value(:,:)
   !  endif
   !enddo
   !if(nsppol==1.and.nspinor==2) then
   !  call gather_matlu(matlu,gathermatlu,natom,option=1,prtopt=1)
   !else if((nsppol==2.or.nsppol==1).and.nspinor==1) then
   !  do iatom=1,natom
   !    if(matlu(iatom)%lpawu.ne.-1) then
   !      tndim=nspinor*(2*matlu(iatom)%lpawu+1)
   !      do im1=1,tndim
   !        do im2=1,tndim
   !          gathermatlu(iatom)%value(im1,im2)=matlu(iatom)%mat(im1,im2,isppol,1,1)
   !        enddo
   !      enddo
   !    endif
   !  enddo
   !endif
        ! write(std_out,*) "gathermatlu in rotate matlu"
        ! do im1=1,tndim
        !   write(message,'(12(1x,18(1x,"(",e17.10,",",e17.10,")")))')&
        !    (gathermatlu(1)%value(im1,im2),im2=1,tndim)
        !   call wrtout(std_out,message,'COLL')
        ! end do

! ===========================
! If necessary, invert rot_mat
! ===========================
   !if(inverse==1) then
   !  do iatom=1,natom
   !    if(matlu(iatom)%lpawu.ne.-1) then
   !      tndim=nspinor*(2*matlu(iatom)%lpawu+1)
   !        do im1=1,tndim
   !          do im2=1,tndim
!               rot_mat(iatom,isppol)%value(im1,im2)=conjg(rot_mat_orig(iatom,isppol)%value(im2,im1))
   !            rot_mat(iatom,isppol)%value(im1,im2)=conjg(rot_mat_orig(iatom)%value(im2,im1))
   !          enddo
   !        enddo
   !    endif ! lpawu
   !  enddo ! iatom
   !endif
        ! write(std_out,*) "rot_mat_orig "
        ! do im1=1,tndim
        !   write(message,'(12(1x,18(1x,"(",e18.10,",",e18.10,")")))')&
        ! &   (rot_mat_orig(1)%value(im1,im2),im2=1,tndim)
        !   call wrtout(std_out,message,'COLL')
        ! end do
        ! write(std_out,*) "rot_mat "
        ! do im1=1,tndim
        !   write(message,'(12(1x,18(1x,"(",e18.10,",",e18.10,")")))')&
        ! &   (rot_mat(1,1)%value(im1,im2),im2=1,tndim)
        !   call wrtout(std_out,message,'COLL')
        ! end do

! ===========================
! Rotate
! ===========================
   !ABI_MALLOC(temp_mat,(tndim,tndim))
   !do iatom=1,natom
   !  if(matlu(iatom)%lpawu.ne.-1) then
   !    tndim=nspinor*(2*matlu(iatom)%lpawu+1)
   !    temp_mat(:,:)=czero
!      input matrix: gathermatlu
!      rotation matrix: rot_mat
!      intermediate matrix: temp_mat
!      result matrix: gathermatlu
       ! temp_mat = gathermatlu * conjg(rot_mat)
   !    call zgemm('n','c',tndim,tndim,tndim,cone,gathermatlu(iatom)%value
   !    ,tndim,&
!&        rot_mat(iatom,isppol)%value,tndim,czero,temp_mat
!,tndim)
       ! gathermatlu = rot_mat * temp_mat = rot_mat * gathermatlu *
       ! conjg(rot_mat)
   !    call
   !    zgemm('n','n',tndim,tndim,tndim,cone,rot_mat(iatom,isppol)%value,tndim,&
!&        temp_mat
!,tndim,czero,gathermatlu(iatom)%value,tndim)
  !   endif ! lpawu
  ! enddo ! iatom
   !do iatom=1,natom
   !  if(matlu(iatom)%lpawu.ne.-1) then
   !    write(std_out,*) "temp_mat in rotate_matlu 2"
   !    do im1=1,tndim
   !      write(message,'(12(1x,18(1x,"(",f17.10,",",f17.10,")")))')&
   !&       (temp_mat(im1,im2),im2=1,tndim)
   !      call wrtout(std_out,message,'COLL')
   !    end do
   !  endif
   !enddo
   !do iatom=1,natom
   !  if(matlu(iatom)%lpawu.ne.-1) then
   !    write(std_out,*) "gathermatlu in rotate_matlu 2"
   !    do im1=1,tndim
   !      write(message,'(12(1x,18(1x,"(",f17.10,",",f17.10,")")))')&
   !&       (gathermatlu(iatom)%value(im1,im2),im2=1,tndim)
   !      call wrtout(std_out,message,'COLL')
   !    end do
   !  endif
   !enddo
  ! ABI_FREE(temp_mat)
     !ABI_ERROR("Aborting now")

! Choose inverse rotation: reconstruct correct rot_mat from rot_mat_orig
! ========================================================================
   !if(inverse==1) then
   !  do iatom=1,natom
   !    if(matlu(iatom)%lpawu.ne.-1) then
   !      tndim=nspinor*(2*matlu(iatom)%lpawu+1)
   !        do im1=1,tndim
   !          do im2=1,tndim
!  !
!  rot_mat(iatom,isppol)%value(im1,im2)=rot_mat_orig(iatom,isppol)%value(im1,im2)
   !            rot_mat(iatom,isppol)%value(im1,im2)=rot_mat_orig(iatom)%value(im1,im2)
   !          enddo
   !        enddo
   !    endif ! lpawu
   !  enddo ! iatom
   !endif

! ===========================
! Put data into matlu(iatom)
! ===========================
   !if(nsppol==1.and.nspinor==2) then
   !  call gather_matlu(matlu,gathermatlu,natom,option=-1,prtopt=1)
   !else if((nsppol==2.or.nsppol==1).and.nspinor==1) then
   !  do iatom=1,natom
   !    if(matlu(iatom)%lpawu.ne.-1) then
   !      tndim=nspinor*(2*matlu(iatom)%lpawu+1)
   !      do im1=1,tndim
   !        do im2=1,tndim
   !          matlu(iatom)%mat(im1,im2,isppol,1,1)=
   !          gathermatlu(iatom)%value(im1,im2)
   !        enddo
   !      enddo
   !    endif
   !  enddo
   !endif ! test nsppol/nspinor
! ===========================
! Deallocations
! ===========================
   !do iatom=1,natom
   !  if(matlu(iatom)%lpawu.ne.-1) then
   !    ABI_FREE(gathermatlu(iatom)%value)
!       ABI_FREE(rot_mat_orig(iatom,isppol)%value)
   !    ABI_FREE(rot_mat_orig(iatom)%value)
   !  endif
   !enddo
   !ABI_FREE(gathermatlu)
   !ABI_FREE(rot_mat_orig)
 !enddo ! isppol

 end subroutine rotate_matlu
!!***

!!****f* m_matlu/shift_matlu
!! NAME
!! shift_matlu
!!
!! FUNCTION
!! Add/subtract a scalar to the diagonal part
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom) :: input quantity to shift
!!  natom=number of atoms in cell.
!!  shift= shift of the diagonal part.
!!  signe= 1 (default) : the shift is substracted
!!       = -1 : the shift is added
!!
!! OUTPUT
!!  matlu(natom) :: shifted matrix
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine shift_matlu(matlu,natom,shift,signe)

!Arguments ------------------------------------
 integer, intent(in) :: natom
 type(matlu_type), intent(inout) :: matlu(natom)
 complex(dp), intent(in) :: shift(natom)
 integer, optional, intent(in) :: signe
!Local variables-------------------------------
 integer :: iatom,im,lpawu,ndim,nspinor,nsppol,signe_used
! character(len=500) :: message
!************************************************************************

 nspinor    = matlu(1)%nspinor
 nsppol     = matlu(1)%nsppol
 signe_used = 1

 if (present(signe)) then
   if (signe == -1) signe_used = -1
 end if ! present(signe)

 do iatom=1,natom
   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim = nspinor * (2*lpawu+1)
   do im=1,ndim
     matlu(iatom)%mat(im,im,:) = matlu(iatom)%mat(im,im,:) + merge(-shift(iatom),shift(iatom),signe_used==1)
   end do ! im
 end do ! iatom

 end subroutine shift_matlu
!!***

!!****f* m_matlu/checkreal_matlu
!! NAME
!! checkreal_matlu
!!
!! FUNCTION
!! Check that matlu is real and diagonal with given precision
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dmft_solv :: impurity solver
!!  matlu(natom) :: input quantity to check
!!  natom=number of atoms in cell.
!!  tol : threshold. Print a warning if max(abs(imag(off diagonal elements))) > tol or
!!        max(abs(off diagonal elements)) > tol, and throws an error if max(abs(imag(diagonal elements))) > tol
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine checkreal_matlu(matlu,natom,tol)

!Arguments ------------------------------------
 real(dp), intent(in) :: tol
 integer, intent(in)  :: natom
 type(matlu_type), intent(in) :: matlu(natom)
!Local variables-------------------------------
 integer :: iatom,im,im1,isppol,lpawu,ndim,nspinor,nsppol
 character(len=500) :: message
 real(dp) :: elem,maximag,maximagdiag,maxoffdiag
!************************************************************************

 maximag     = zero
 maximagdiag = zero
 maxoffdiag  = zero
 nspinor     = matlu(1)%nspinor
 nsppol      = matlu(1)%nsppol

 do iatom=1,natom
   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim = nspinor * (2*lpawu+1)
   do isppol=1,nsppol
     do im1=1,ndim
       do im=1,ndim
         elem = abs(aimag(matlu(iatom)%mat(im,im1,isppol)))
         if (elem > maximag) maximag = elem
         if (im == im1) then
           if (elem > maximagdiag) maximagdiag = elem
         else
           elem = abs(matlu(iatom)%mat(im,im1,isppol))
           if (elem > maxoffdiag) maxoffdiag = elem
         end if ! im/=im1
       end do ! im
     end do ! im1
   end do ! isppol
 end do ! iatom

 if (maximagdiag > tol) then
   write(message,'(3x,2a,e12.4,a,e12.4,2a)') ch10,&
     & ' Diagonal part of the occupation matrix is complex: the imaginary part ',&
     & maximagdiag,' is larger than',tol,ch10,  &
     & "The calculation cannot handle it : check that your calculation is meaningful"
   ABI_ERROR(message)
 end if ! maximagdiag > tol
 if (maximag > tol) then
   write(message,'(3x,2a,e12.4,a,e12.4,2a)') ch10,&
     & ' The off diagonal occupation matrix is complex: the imaginary part ',maximag,' is larger than',tol,ch10,&
     & "Check that your calculation is meaningful"
   ABI_WARNING(message)
 end if ! maximag > tol
 if (maxoffdiag > tol) then
   write(message,'(3x,2a,e12.4,a,e12.4,6a)') ch10,&
        & ' Occupation matrix is non diagonal : the maximum off-diag part ',maxoffdiag,' is larger than',tol,ch10,&
        & "The corresponding non diagonal elements will be neglected in the Weiss/Hybridization functions",ch10,&
        & "(Except if dmft_solv=8,9 where these elements are taken into account)",ch10,"This is an approximation."
   ABI_WARNING(message)
 else
   write(message,'(3x,2a,e12.4,a,e12.4,2a)') ch10,' Occupation matrix is diagonal : the off-diag part ',&
     & maxoffdiag,' is lower than',tol
   ABI_COMMENT(message)
 end if ! maxoffdiag > tol

 end subroutine checkreal_matlu
!!***

!!****f* m_matlu/checkdiag_matlu
!! NAME
!! checkdiag_matlu
!!
!! FUNCTION
!! Check that matlu is diagonal in the orbital index with given precision
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom) :: input quantity to check
!!  natom=number of atoms in cell.
!!  tol : precision
!!
!! OUTPUT
!!  nondiag= true if max(abs(off diagonal elements)) > tol
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine checkdiag_matlu(matlu,natom,tol,nondiag)

!Arguments ------------------------------------
 real(dp), intent(in) :: tol
 integer, intent(in)  :: natom
 logical, intent(out) :: nondiag
 type(matlu_type), intent(in) :: matlu(natom)
!Local variables-------------------------------
 integer :: iatom,im,im1,isppol,lpawu,ndim,nsppol,nspinor
!************************************************************************

 nondiag = .false.
 nspinor = matlu(1)%nspinor
 nsppol  = matlu(1)%nsppol

 do iatom=1,natom
   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim = nspinor * (2*lpawu+1)
   do isppol=1,nsppol
     do im1=1,ndim
       do im=1,ndim
!              if(im/=im1) write(std_out,*) "im,im1",im,im1,matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor)
 !            if(present(nondiag).eqv..false.) then
 !              if(im/=im1.and.(abs(matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor))>tol))  then
 !                write(message,'(5i5)') im,im1,isppol,ispinor,ispinor
 !                call wrtout(std_out,message,'COLL')
 !                write(message,'(a,3e16.5)')" checkdiag_matlu: Warning ",matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor),tol
 !                call wrtout(std_out,message,'COLL')
 !                if(.not.present(opt)) ABI_ERROR("not present(opt)")
 !                if(matlu(1)%nspinor==1) ABI_ERROR("matlu%nspinor==1")
 !              endif
!             endif
 !              if(present(nondiag)) then
         if ((im /= im1) .and. abs(dble(matlu(iatom)%mat(im,im1,isppol))) > tol) nondiag = .true.
                  ! write(6,*) "NONDIAG", matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)
               !if(ispinor/=ispinor1.and.(abs(matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))>tol))  then
               !  write(message,'(a,3e16.5)')" checkdiag_matlu :i Warning ",matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1),tol
               !  call wrtout(std_out,message,'COLL')
               !  write(message,'(5i5)') im,im1,isppol,ispinor,ispinor
               !  call wrtout(std_out,message,'COLL')
               !  if(matlu(1)%nspinor==1) ABI_ERROR("matlu%nspinor==1")
               !endif
       end do ! im
     end do ! im1
   end do ! isppol
 end do ! iatom

 end subroutine checkdiag_matlu
!!***

!!****f* m_matlu/prod_matlu
!! NAME
!! prod_matlu
!!
!! FUNCTION
!! Matrix product of two matlus
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu1(natom) :: input quantity
!!  matlu2(natom) :: input quantity
!!
!! OUTPUT
!!  matlu3(natom) :: output quantity
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine prod_matlu(matlu1,matlu2,matlu3,natom)

!Arguments ------------------------------------
 integer, intent(in) :: natom
 type(matlu_type), intent(in) :: matlu1(natom),matlu2(natom)
 type(matlu_type), intent(inout) :: matlu3(natom)
!Local variables-------------------------------
 integer :: iatom,isppol,lpawu,ndim,nspinor,nsppol
!************************************************************************

 nspinor = matlu1(1)%nspinor
 nsppol  = matlu1(1)%nsppol

 do iatom=1,natom
   lpawu = matlu1(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim = nspinor * (2*lpawu+1)
   do isppol=1,nsppol
     call abi_xgemm("n","n",ndim,ndim,ndim,cone,matlu1(iatom)%mat(:,:,isppol),ndim,&
                  & matlu2(iatom)%mat(:,:,isppol),ndim,czero,matlu3(iatom)%mat(:,:,isppol),ndim)
   end do ! isppol
 end do ! iatom

  !call zero_matlu(matlu3,natom)
 !do iatom=1,natom
 !  lpawu=matlu1(iatom)%lpawu
 !  if(lpawu.ne.-1) then
 !    do isppol=1,matlu1(1)%nsppol
 !      do ispinor1=1,matlu1(1)%nspinor
 !        do ispinor2=1,matlu1(1)%nspinor
 !          do ispinor3=1,matlu1(1)%nspinor
 !            do im1=1,2*lpawu+1
 !              do im2=1,2*lpawu+1
 !                do im3=1,2*lpawu+1
 !                  matlu3(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)= &
!&                    matlu3(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)+ &
!&                    matlu1(iatom)%mat(im1,im3,isppol,ispinor1,ispinor3)*&
!&                    matlu2(iatom)%mat(im3,im2,isppol,ispinor3,ispinor2)
 !                enddo ! im3
 !              enddo ! im2
 !            enddo ! im1
 !          enddo ! ispinor3
 !        enddo ! ispinor2
 !      enddo ! ispinor1
 !    enddo ! isppol
 !  endif ! lpawu
 !enddo ! iatom

 end subroutine prod_matlu
!!***

!!****f* m_matlu/conjg_matlu
!! NAME
!! conjg_matlu
!!
!! FUNCTION
!! conjugate of input matlu
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu1(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE
 subroutine conjg_matlu(matlu1,natom)
 use defs_wvltypes

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
!arrays
 type(matlu_type), intent(inout) :: matlu1(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im1,im2,ispinor2,ispinor1,isppol
 integer :: lpawu,ndim
!arrays
!************************************************************************
 do iatom=1,natom
   lpawu=matlu1(iatom)%lpawu
   if(lpawu.ne.-1) then
     ndim=2*lpawu+1
     do isppol=1,matlu1(1)%nsppol
       do ispinor1=1,matlu1(1)%nspinor
         do ispinor2=1,matlu1(1)%nspinor
           do im1=1,2*lpawu+1
             do im2=1,2*lpawu+1
               matlu1(iatom)%mat(im1+(ispinor1-1)*ndim,im2+(ispinor2-1)*ndim,isppol)= &
&               conjg(matlu1(iatom)%mat(im1+(ispinor1-1)*ndim,im2+(ispinor2-1)*ndim,isppol))
             enddo ! im2
           enddo ! im1
         enddo ! ispinor2
       enddo ! ispinor1
     enddo ! isppol
   endif ! lpawu
 enddo ! iatom

 end subroutine conjg_matlu
!!***

!!****f* m_matlu/ln_matlu
!! NAME
!! ln_matlu
!!
!! FUNCTION
!! Compute the logarithm of matlu (only if diagonal for the moment)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu1(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity

!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE
 subroutine ln_matlu(matlu1,natom)
 use defs_wvltypes

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
!arrays
 type(matlu_type), intent(inout) :: matlu1(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im,ispinor,isppol
 integer :: lpawu,ndim
 character(len=500) :: message
!arrays
!************************************************************************
 !call checkdiag_matlu(matlu1,natom,tol8)
 do iatom=1,natom
   lpawu=matlu1(iatom)%lpawu
   if(lpawu.ne.-1) then
     ndim=2*lpawu+1
     do isppol=1,matlu1(1)%nsppol
       do ispinor=1,matlu1(1)%nspinor
         do im=1,2*lpawu+1
           if( real(matlu1(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol))<zero) then
             write(message,'(2a,2es13.5,a)') ch10," ln_matlu: PROBLEM " &
&             , matlu1(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol)
             ABI_ERROR(message)
           endif
           matlu1(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol)= &
&           log(matlu1(iatom)%mat(im+(ispinor-1)*ndim,im+(ispinor-1)*ndim,isppol))
         enddo ! im
       enddo ! ispinor
     enddo ! isppol
   endif ! lpawu
 enddo ! iatom

 end subroutine ln_matlu
!!***

!!****f* m_matlu/slm2ylm_matlu
!! NAME
!! slm2ylm_matlu
!!
!! FUNCTION
!! Transform mat from Slm to Ylm basis or vice versa
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu1(natom) :: input quantity
!!  natom :: number of atoms
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  option=1 go from Slm to Ylm basis
!!  option=2 go from Ylm to Slm basis
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine slm2ylm_matlu(matlu,natom,paw_dmft,option,optprt)

!Arguments ------------------------------------
 integer, intent(in) :: natom,option,optprt
 type(matlu_type), target, intent(inout) :: matlu(natom)
 type(paw_dmft_type), target, intent(in) :: paw_dmft
!Local variables-------------------------------
 integer :: iatom,im1,im2,ispin,ispinor1,ispinor2,isppol
 integer :: lpawu,ndim,ndim_max,nspin,nspinor,nsppol
 complex(dp), pointer :: mat_out(:,:) => null(), slm2ylm(:,:) => null()
 complex(dp), allocatable :: mat_inp(:,:),mat_tmp(:,:)
 complex(dp), target, allocatable :: mat_tmp2(:,:)
 character(len=1) :: c1,c2
 character(len=500) :: message
!************************************************************************

 ndim_max = 2*paw_dmft%maxlpawu + 1
 nspinor  = paw_dmft%nspinor
 nsppol   = paw_dmft%nsppol
 nspin    = nsppol * (nspinor**2)

 if (option == 1) then
   c1 = "n" ; c2 = "c"
 else if (option == 2) then
   c1 = "c" ; c2 = "n"
 end if

 do iatom=1,natom

   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   ndim = 2*lpawu + 1
   slm2ylm => paw_dmft%slm2ylm(:,:,lpawu+1)

   if (optprt > 2) then
     write(message,'(2a)') ch10,"SLM2YLM matrix"
     call wrtout(std_out,message,'COLL')
     do im1=1,ndim
       write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))') (slm2ylm(im1,im2),im2=1,ndim)
       call wrtout(std_out,message,'COLL')
     end do ! im1
   end if ! optprt>2

   ABI_MALLOC(mat_inp,(ndim,ndim))
   ABI_MALLOC(mat_tmp,(ndim,nspin*ndim))
   ABI_MALLOC(mat_tmp2,(ndim,nspin*ndim))

   ispin = 0

   do isppol=1,nsppol
     do ispinor2=1,nspinor
       do ispinor1=1,nspinor

         ispin = ispin + 1

         ! Make copy here instead of creating a temporary when calling zgemm in order to please -fcheck
         mat_inp(:,:) = matlu(iatom)%mat(1+(ispinor1-1)*ndim:ndim*ispinor1,1+(ispinor2-1)*ndim:ispinor2*ndim,isppol)

         if (optprt > 2) then
           write(message,'(2a,i2,a,i2,a,i2)') ch10,"SLM input matrix,&
             & isppol=",isppol,", ispinor1=",ispinor1,", ispinor2=",ispinor2
           call wrtout(std_out,message,'COLL')
           do im1=1,ndim
             write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))') &
               & (mat_inp(im1,im2),im2=1,ndim)
             call wrtout(std_out,message,'COLL')
           end do ! im1
         end if ! optprt>2

         call abi_xgemm("n",c2,ndim,ndim,ndim,cone,mat_inp(:,:),ndim,slm2ylm(:,1:ndim), &
                      & ndim_max,czero,mat_tmp(:,1+(ispin-1)*ndim:ndim*ispin),ndim)

       end do ! ispinor2
     end do ! ispinor1
   end do ! isppol

   call abi_xgemm(c1,"n",ndim,ndim*nspin,ndim,cone,slm2ylm(:,1:ndim),ndim_max,mat_tmp(:,:),ndim, &
                & czero,mat_tmp2(:,:),ndim)

   ispin = 0

   do isppol=1,nsppol
     do ispinor2=1,nspinor
       do ispinor1=1,nspinor

         ispin = ispin + 1

         mat_out => mat_tmp2(:,1+(ispin-1)*ndim:ndim*ispin)

         matlu(iatom)%mat(1+(ispinor1-1)*ndim:ndim*ispinor1,1+(ispinor2-1)*ndim:ispinor2*ndim,isppol) = mat_out(:,:)

         if (optprt > 2) then
           write(message,'(2a,i2,a,i2,a,i2)') ch10,"YLM output matrix, isppol=",isppol,", ispinor=",ispinor1,&
              & ", ispinor2=",ispinor2
           call wrtout(std_out,message,'COLL')
           do im1=1,ndim
             write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))') (mat_out(im1,im2),im2=1,ndim)
             call wrtout(std_out,message,'COLL')
           end do ! im1
         end if ! optprt>2

       end do ! ispinor1
     end do ! ispinor2
   end do ! isppol

   ABI_FREE(mat_inp)
   ABI_FREE(mat_tmp)
   ABI_FREE(mat_tmp2)

 end do ! iatom

 mat_out => null()
 slm2ylm => null()

 !do iatom=1,natom
 !  lpawu=matlu(iatom)%lpawu
 !  if(lpawu.ne.-1) then
 !    ndim=2*lpawu+1
 !    ll=lpawu
 !    ABI_MALLOC(slm2ylm,(2*ll+1,2*ll+1))
 !    slm2ylm=czero
 !    do im=1,2*ll+1
 !      mm=im-ll-1;jm=-mm+ll+1
 !      onem=dble((-1)**mm)
 !      if (mm> 0) then
 !        slm2ylm(im,im)= cmplx(onem*invsqrt2,zero,kind=dp)
 !        slm2ylm(jm,im)= cmplx(invsqrt2,     zero,kind=dp)
 !      end if
 !      if (mm==0) then
 !        slm2ylm(im,im)=cone
 !      end if
 !      if (mm< 0) then
 !        slm2ylm(im,im)= cmplx(zero,     invsqrt2,kind=dp)
 !        slm2ylm(jm,im)=-cmplx(zero,onem*invsqrt2,kind=dp)
 !      end if
 !    end do
 !    if(optprt>2) then
 !      write(message,'(2a)') ch10,"SLM2YLM matrix"
 !      call wrtout(std_out,message,'COLL')
 !      do im1=1,ll*2+1
 !        write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
!&         (slm2ylm(im1,im2),im2=1,ll*2+1)
 !        call wrtout(std_out,message,'COLL')
 !      end do
 !    endif
 !    do isppol=1,matlu(1)%nsppol
 !      do ispinor=1,matlu(1)%nspinor
 !        do ispinor2=1,matlu(1)%nspinor
 !          ABI_MALLOC(mat_out_c,(2*ll+1,2*ll+1))
 !          ABI_MALLOC(mat_inp_c,(2*ll+1,2*ll+1))
 !          mat_inp_c(:,:) = matlu(iatom)%mat(1+(ispinor-1)*ndim:ndim+(ispinor-1)*ndim,1+(ispinor2-1)*ndim:ndim+(ispinor2-1)*ndim,isppol)
 !          mat_out_c=czero

 !          if(optprt>2) then
 !            write(message,'(2a, i2, a, i2, a, i2)') ch10,"SLM input matrix, isppol=", isppol, ", ispinor=", ispinor,&
!&             ", ispinor2=", ispinor2
 !            call wrtout(std_out,message,'COLL')
 !            do im1=1,ll*2+1
 !              write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
!&               (mat_inp_c(im1,im2),im2=1,ll*2+1)
 !              call wrtout(std_out,message,'COLL')
 !            end do
 !          endif

  !         do jm=1,2*ll+1
  !           do im=1,2*ll+1
  !             tmp2=czero
  !             do ii=1,2*ll+1
  !               do jj=1,2*ll+1
  !                 if(option==1) then
  !                   tmp2=tmp2+mat_inp_c(ii,jj)*(slm2ylm(im,ii))*CONJG(slm2ylm(jm,jj))
  !                 else if(option==2) then
  !                   tmp2=tmp2+mat_inp_c(ii,jj)*CONJG(slm2ylm(ii,im))*(slm2ylm(jj,jm))
  !                 end if
  !               end do
  !             end do
  !             mat_out_c(im,jm)=tmp2
  !           end do
  !         end do

  !         if(optprt>2) then
  !           write(message,'(2a, i2, a, i2, a, i2)') ch10,"YLM output matrix, isppol=", isppol, ", ispinor=", ispinor,&
!&             ", ispinor2=", ispinor2
  !           call wrtout(std_out,message,'COLL')
  !           do im1=1,ll*2+1
  !             write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
  !    &         (mat_out_c(im1,im2),im2=1,ll*2+1)
  !             call wrtout(std_out,message,'COLL')
  !           end do
  !         endif

  !         matlu(iatom)%mat(1+(ispinor-1)*ndim:ndim+(ispinor-1)*ndim,1+(ispinor2-1)*ndim:ndim+(ispinor2-1)*ndim,isppol)=mat_out_c(:,:)
  !         ABI_FREE(mat_out_c)
  !         ABI_FREE(mat_inp_c)
  !      enddo ! im
 !      enddo ! ispinor
 !    enddo ! isppol
 !    ABI_FREE(slm2ylm)
 !  endif ! lpawu
 !enddo ! iatom

 end subroutine slm2ylm_matlu
!!***

!!****f* m_matlu/fac_matlu
!! NAME
!! fac_matlu
!!
!! FUNCTION
!! Multiply matlu by a scalar
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom) :: input quantity
!!  natom=number of atoms in cell.
!!  fac= factor
!!
!! OUTPUT
!!  matlu(natom) :: fac * matlu
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine fac_matlu(matlu,natom,fac)

!Arguments ------------------------------------
 integer, intent(in) :: natom
 type(matlu_type), intent(inout) :: matlu(natom)
 complex(dp), intent(in) :: fac
!Local variables-------------------------------
 integer :: iatom,lpawu
! character(len=500) :: message
!************************************************************************

 do iatom=1,natom
   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   matlu(iatom)%mat(:,:,:) = fac * matlu(iatom)%mat(:,:,:)
 end do ! iatom

 end subroutine fac_matlu
!!***

!!****f* m_matlu/printplot_matlu
!! NAME
!! printplot_matlu
!!
!! FUNCTION
!! Write matlu for a given frequency
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom) :: input quantity to write
!!  natom=number of atoms in cell.
!!  freq :: frequency
!!  char1 :: name of the file on which to write
!!  units :: unit of the file
!!  imre :: if present, write real and imaginary parts in two different files
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine printplot_matlu(matlu,natom,freq,char1,units,imre)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,units
 integer, optional, intent(in) :: imre
 real(dp), intent(in) :: freq
!arrays
 type(matlu_type), intent(in) :: matlu(natom)
 character(len=*), intent(in) :: char1
!Local variables-------------------------------
!scalars
 integer :: iatom,im,im1,ispinor,ispinor1,isppol
 integer :: lpawu,ndim,nspinor,nsppol,unitnb
 character(len=4) :: tag_at
 character(len=fnlen) :: tmpfil,tmpfilim,tmpfilre
! character(len=500) :: message
!arrays
!************************************************************************

 nspinor = matlu(1)%nspinor
 nsppol  = matlu(1)%nsppol

 ! not yet tested and used
 do iatom=1,natom
   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim   = 2*lpawu + 1
   unitnb = units + iatom
   call int2char4(iatom,tag_at)
   if (present(imre)) then
     tmpfilre = trim(char1)//tag_at//"re"
     tmpfilim = trim(char1)//tag_at//"im"
     open(unit=unitnb+10,file=trim(tmpfilre),status='unknown',form='formatted')
     open(unit=unitnb+20,file=trim(tmpfilim),status='unknown',form='formatted')
     write(unitnb+10,'(400e26.16)') freq,(((((dble(matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,isppol)),&
      & im=1,ndim),ispinor=1,nspinor),im1=1,ndim),ispinor1=1,nspinor),isppol=1,nsppol)
     write(unitnb+20,'(400e26.16)') freq,(((((aimag(matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,isppol)),&
      & im=1,ndim),ispinor=1,nspinor),im1=1,ndim),ispinor1=1,nspinor),isppol=1,nsppol)
   else
     tmpfil = trim(char1)//tag_at
     open(unit=unitnb,file=trim(tmpfil),status='unknown',form='formatted')
     write(unitnb,'(400e26.16)') freq,(((((matlu(iatom)%mat(im+(ispinor-1)*ndim,im1+(ispinor1-1)*ndim,isppol),&
      & im=1,ndim),ispinor=1,nspinor),im1=1,ndim),ispinor1=1,nspinor),isppol=1,nsppol)
   end if ! present(imre)
 end do ! iatom

 end subroutine printplot_matlu
!!***

!!****f* m_matlu/identity_matlu
!! NAME
!! identity_matlu
!!
!! FUNCTION
!!  Set the diagonal elements to 1 (the off-diagonal are not set to 0)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom) :: input quantity
!!  natom=number of atoms in cell.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine identity_matlu(matlu,natom)

!Arguments ------------------------------------
 integer, intent(in) :: natom
 type(matlu_type), intent(inout) :: matlu(natom)
!Local variables-------------------------------
 integer :: iatom,im,lpawu,ndim,nspinor
! character(len=500) :: message
!arrays
!************************************************************************

 nspinor = matlu(1)%nspinor

 do iatom=1,natom
   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim = nspinor * (2*lpawu+1)
   do im=1,ndim
     matlu(iatom)%mat(im,im,:) = cone
   end do ! im
 end do ! iatom

 end subroutine identity_matlu
!!***

!!***
!!****f* m_matlu/magmomforb_matlu
!! NAME
!! magmomforb_matlu
!!
!! FUNCTION
!! return the product of occupation matrix of dimension [(2*ll+1)]**4 in the Ylm basis
!! with the matrix of orbital angular momentum element for the x,y and z direction.
!! Option gives the direction of the magnetic moment x==1, y==2 and z==3.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (FGendron)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! matlu1(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity in Ylm basis
!! natom :: number of atoms
!! option = 1 :: x axis
!!	  = 2 :: y axis
!!        = 3 :: z axis
!! optptr > 2 :: print orbital angular matrix elements and resulting product
!!
!! OUTPUT
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: product
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE
 subroutine magmomforb_matlu(matlu,mu,natom,option,optprt)
 use defs_wvltypes

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,option,optprt
 complex(dp), allocatable, intent(inout) :: mu(:)
!arrays
 type(matlu_type), intent(inout) :: matlu(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im,ispinor,isppol,ispinor2
 integer :: lpawu,ll,jm,ml1,jc1,ms1,lcor,im1,im2,tndim,ndim
 character(len=500) :: message
 real(dp) :: xj
!arrays
 complex(dp),allocatable :: mat_out_c(:,:)
! integer, allocatable :: ind_msml(:,:)
 complex(dp), allocatable :: temp_mat(:,:)
 type(coeff2c_type), allocatable :: gathermatlu(:)
 type(coeff2c_type), allocatable :: muorb(:)
!************************************************************************

 !=====================================
 ! Allocate matrices
 !=====================================

 ABI_MALLOC(gathermatlu,(natom))
 ABI_MALLOC(muorb,(natom))
 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     tndim=2*(2*matlu(iatom)%lpawu+1)
     ABI_MALLOC(gathermatlu(iatom)%value,(tndim,tndim))
     gathermatlu(iatom)%value=czero
     ABI_MALLOC(muorb(iatom)%value,(tndim,tndim))
     muorb(iatom)%value=czero
   end if
 end do

 do iatom=1,natom
   lpawu=matlu(iatom)%lpawu
   if(lpawu.ne.-1) then

     ll=lpawu
     lcor=lpawu
     !=====================================
     !build orbital angular momentum matrix along x axis
     !=====================================
     if(option==1) then

       jc1=0
       do ms1 =-1,1
         xj=float(ms1)+half
         do ml1 = -ll,ll
            jc1=jc1+1
            if(jc1 == 1) then
               muorb(iatom)%value(jc1,jc1+1) = sqrt(float((lcor*(lcor+1)) - ml1*(ml1 + 1)))*0.5
            endif
            if(jc1 > 1 .and. jc1 < 2*(2*ll+1)) then
               muorb(iatom)%value(jc1,jc1+1) = sqrt(float((lcor*(lcor+1)) - ml1*(ml1 + 1)))*0.5
               muorb(iatom)%value(jc1,jc1-1) = sqrt(float((lcor*(lcor+1)) - ml1*(ml1 - 1)))*0.5
            else if(jc1== 2*(2*ll+1)) then
               muorb(iatom)%value(jc1,jc1-1) = sqrt(float((lcor*(lcor+1)) - ml1*(ml1 - 1)))*0.5
            end if
          end do
        end do

     !=====================================
     !build orbital angular momentum matrix along y axis
     !=====================================
     else if(option==2) then

       jc1=0
       do ms1 =-1,1
         xj=float(ms1)+half
         do ml1 = -ll,ll
            jc1=jc1+1
            if(jc1 == 1) then
               muorb(iatom)%value(jc1,jc1+1) = cmplx(zero,sqrt(float((lcor*(lcor+1)) - ml1*(ml1 + 1)))*0.5,kind=dp)
            endif
            if(jc1 > 1 .and. jc1 < 2*(2*ll+1)) then
               muorb(iatom)%value(jc1,jc1+1) = cmplx(zero,sqrt(float((lcor*(lcor+1)) - ml1*(ml1 + 1)))*0.5,kind=dp)
               muorb(iatom)%value(jc1,jc1-1) = cmplx(zero,-sqrt(float((lcor*(lcor+1)) - ml1*(ml1 - 1)))*0.5,kind=dp)
            else if(jc1 == 2*(2*ll+1)) then
               muorb(iatom)%value(jc1,jc1-1) = cmplx(zero,-sqrt(float((lcor*(lcor+1)) - ml1*(ml1 - 1)))*0.5,kind=dp)
            end if
          end do
        end do

     !=====================================
     !build orbital angular momentum matrix along z axis
     !=====================================
     else if(option==3) then
       jc1=0
       do ms1=-1,1
         do ml1=-ll,ll
            jc1=jc1+1
            if(jc1 < tndim+1) then
              muorb(iatom)%value(jc1,jc1) = ml1
            endif
          end do
        end do
     end if

     if(optprt>2) then
        write(message,'(a,i4)') "Orbital angular momentum matrix elements in |m_l,m_s> basis for axis=", option
        call wrtout(std_out,message,"COLL")
        do im=1,2*(ll*2+1)
          write(message,'(6(1x,9(1x,f4.1,",",f4.1)))') (muorb(iatom)%value(im,jm),jm=1,2*(ll*2+1))
          call wrtout(std_out,message,"COLL")
        end do
     end if

   end if !lpawu
 end do !atom

     !=====================================
     ! Reshape input Ylm matlu in one 14x14 matrix
     !=====================================

 call gather_matlu(matlu,gathermatlu,natom,option=1,prtopt=1)

!!printing for debug
! write(std_out,*) "gathermatlu in magmomforb"
! do im1=1,tndim
!   write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
!        (gathermatlu(1)%value(im1,im2),im2=1,tndim)
!   call wrtout(std_out,message,'COLL')
! end do

     !=====================================
     ! Matrix product of Occ and muorb
     !=====================================
 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     tndim=2*(2*matlu(iatom)%lpawu+1)
     ABI_MALLOC(temp_mat,(tndim,tndim))

     call zgemm('n','n',tndim,tndim,tndim,cone,gathermatlu(iatom)%value,tndim,muorb(iatom)%value,tndim,czero,temp_mat,tndim)

     gathermatlu(iatom)%value=temp_mat
     ABI_FREE(temp_mat)

     !=====================================
     ! Trace of matrix product
     !=====================================

   do im1=1,tndim
     do im2=1,tndim
       if(im1==im2) then
         mu(iatom) = mu(iatom) + gathermatlu(iatom)%value(im1,im2)
       end if
     end do
   end do


     !=====================================
     ! Reshape product matrix into matlu format
     !=====================================

 !call gather_matlu(matlu,gathermatlu,natom,option=-1,prtopt=1)


     if(optprt>2) then
       ABI_MALLOC(mat_out_c,(2*ll+1,2*ll+1))
       ndim = 2*ll+1
       do isppol=1,matlu(1)%nsppol
         do ispinor=1,matlu(1)%nspinor
           do ispinor2=1,matlu(1)%nspinor
             mat_out_c(:,:) = matlu(iatom)%mat(1+(ispinor-1)*ndim:ndim+(ispinor-1)*ndim,1+(ispinor2-1)*ndim:ndim+(ispinor2-1)*ndim,isppol)

             write(message,'(2a, i2, a, i2, a, i2)') ch10,"Orbital angular momentum matrix, isppol=", isppol, ", ispinor=",&
&            ispinor,", ispinor2=", ispinor2
             call wrtout(std_out,message,'COLL')
             do im1=1,ll*2+1
               write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
      &         (mat_out_c(im1,im2),im2=1,ll*2+1)
               call wrtout(std_out,message,'COLL')
             end do

           end do ! ispinor2
         end do ! ispinor
       end do ! isppol
       ABI_FREE(mat_out_c)
     endif

   end if !lpawu
 end do !atom

     !=====================================
     ! Deallocate gathermatlu
     !=====================================

 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     ABI_FREE(gathermatlu(iatom)%value)
     ABI_FREE(muorb(iatom)%value)
   end if
 end do
 ABI_FREE(gathermatlu)
 ABI_FREE(muorb)

 end subroutine magmomforb_matlu

!!***


!!***
!!****f* m_matlu/magmomfspin_matlu
!! NAME
!! magmomfspin_matlu
!!
!! FUNCTION
!! return the product of occupation matrix of dimension [(2*ll+1)]**4 in the Ylm basis
!! with the matrix of spin angular momentum element for the x,y and z direction.
!! Option gives the direction of the magnetic moment x==1, y==2 and z==3.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (FGendron)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! matlu1(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity in Ylm basis
!! natom :: number of atoms
!! option = 1 :: x axis
!!	  = 2 :: y axis
!!        = 3 :: z axis
!! optptr > 2 :: print spin angular matrix elements and resulting product
!!
!! OUTPUT
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: product
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE
 subroutine magmomfspin_matlu(matlu,mu,natom,option,optprt)
 use defs_wvltypes

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,option,optprt
 complex(dp), allocatable, intent(inout) :: mu(:)
!arrays
 type(matlu_type), intent(inout) :: matlu(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im,ispinor,isppol,ispinor2
 integer :: lpawu,ll,jm,ml1,jc1,ms1,lcor,im1,im2,tndim,ndim
 character(len=500) :: message
 real(dp) :: xj
!arrays
 complex(dp),allocatable :: mat_out_c(:,:)
 integer, allocatable :: ind_msml(:,:)
 complex(dp), allocatable :: temp_mat(:,:)
 type(coeff2c_type), allocatable :: gathermatlu(:)
 type(coeff2c_type), allocatable :: muspin(:)
!************************************************************************

 !=====================================
 ! Allocate matrices
 !=====================================

 ABI_MALLOC(gathermatlu,(natom))
 ABI_MALLOC(muspin,(natom))
 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     tndim=2*(2*matlu(iatom)%lpawu+1)
     ABI_MALLOC(gathermatlu(iatom)%value,(tndim,tndim))
     gathermatlu(iatom)%value=czero
     ABI_MALLOC(muspin(iatom)%value,(tndim,tndim))
     muspin(iatom)%value=czero
   end if
 end do

 do iatom=1,natom
   lpawu=matlu(iatom)%lpawu
   if(lpawu.ne.-1) then
     ll=lpawu
     lcor=lpawu
     ABI_MALLOC(ind_msml,(2,-ll:ll))
     ind_msml=czero
     !=====================================
     !build spin angular momentum matrix along x axis
     !=====================================
     if(option==1) then
       jc1=0
       do ms1=1,2
         do ml1=-ll,ll
          jc1=jc1+1
          ind_msml(ms1,ml1)=jc1
         end do
       end do

       jc1=0
       do ms1 =-1,1
         xj=float(ms1)+half
         do ml1 = -ll,ll
            jc1=jc1+1
            if(xj < 0.0 ) then
               muspin(iatom)%value(ind_msml(2,ml1),jc1) = 0.5
            else if(xj > 0.0) then
               muspin(iatom)%value(ind_msml(1,ml1),ind_msml(2,ml1)) = 0.5
            end if
          end do
        end do

     !=====================================
     !build spin angular momentum matrix along y axis
     !up spin is first
     !=====================================
     else if(option==2) then
        jc1=0
       do ms1=1,2
         do ml1=-ll,ll
          jc1=jc1+1
          ind_msml(ms1,ml1)=jc1
         end do
       end do

       jc1=0
       do ms1 =-1,1
         xj=float(ms1)+half
         do ml1 = -ll,ll
            jc1=jc1+1
            if(xj < 0.0 ) then
               muspin(iatom)%value(ind_msml(2,ml1),jc1) = cmplx(zero,0.5,kind=dp)
            else if(xj > 0.0) then
               muspin(iatom)%value(ind_msml(1,ml1),ind_msml(2,ml1)) = cmplx(zero,-0.5,kind=dp)
            end if
          end do
        end do

     !=====================================
     !build spin angular momentum matrix along z axis
     !up spin is first
     !=====================================
     else if(option==3) then
       jc1=0
       do ms1=-1,1
         xj=float(ms1)+half
         do ml1=-ll,ll
            jc1=jc1+1
            if(jc1 < tndim+1) then
              if(xj < 0.0 ) then
                 muspin(iatom)%value(jc1,jc1) = -xj
              else if(xj > 0.0) then
                 muspin(iatom)%value(jc1,jc1) = -xj
              end if
            endif
         end do
       end do
     end if

     if(optprt>2) then
        write(message,'(a,i4)') "Spin angular momentum matrix elements in |m_l,m_s> basis for axis", option
        call wrtout(std_out,message,"COLL")
        do im=1,2*(ll*2+1)
          write(message,'(6(1x,9(1x,f4.1,",",f4.1)))') (muspin(iatom)%value(im,jm),jm=1,2*(ll*2+1))
          call wrtout(std_out,message,"COLL")
        end do
     end if

   ABI_FREE(ind_msml)
   end if !lpawu
 end do !atom

     !=====================================
     ! Reshape input Ylm matlu in one 14x14 matrix
     !=====================================

 call gather_matlu(matlu,gathermatlu,natom,option=1,prtopt=1)

!!printing for debug
!! write(std_out,*) "gathermatlu in magmomfspin"
!! do im1=1,tndim
!!   write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
!!        (gathermatlu(1)%value(im1,im2),im2=1,tndim)
!!   call wrtout(std_out,message,'coll')
!! end do

     !=====================================
     ! Matrix product of Occ and muspin
     !=====================================
 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     tndim=2*(2*matlu(iatom)%lpawu+1)
     ABI_MALLOC(temp_mat,(tndim,tndim))

     call zgemm('n','n',tndim,tndim,tndim,cone,gathermatlu(iatom)%value,tndim,muspin(iatom)%value,tndim,czero,temp_mat,tndim)

     gathermatlu(iatom)%value=temp_mat
     ABI_FREE(temp_mat)

     !!printing for debug
     !!write(std_out,*) "gathermatlu in magmomfspin after product"
     !!do im1=1,tndim
     !!   write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
     !!        (gathermatlu(1)%value(im1,im2),im2=1,tndim)
     !!   call wrtout(std_out,message,'coll')
     !!end do

     !=====================================
     ! Trace of matrix product
     !=====================================

   do im1=1,tndim
     do im2=1,tndim
       if(im1==im2) then
         mu(iatom) = mu(iatom) + gathermatlu(iatom)%value(im1,im2)
       end if
     end do
   end do

     !=====================================
     ! Reshape product matrix into matlu format
     !=====================================

    !call gather_matlu(matlu,gathermatlu,natom,option=-1,prtopt=1)


     !=====================================
     ! Print matlu
     !=====================================
     if(optprt>2) then
       ABI_MALLOC(mat_out_c,(2*ll+1,2*ll+1))
       ndim = 2*ll+1
       do isppol=1,matlu(1)%nsppol
         do ispinor=1,matlu(1)%nspinor
           do ispinor2=1,matlu(1)%nspinor
             mat_out_c(:,:) = matlu(iatom)%mat(1+(ispinor-1)*ndim:ndim+(ispinor-1)*ndim,1+(ispinor2-1)*ndim:ndim+(ispinor2-1)*ndim,isppol)

             write(message,'(2a, i2, a, i2, a, i2)') ch10,"Spin angular momentum matrix, isppol=", isppol, ", ispinor=", ispinor,&
&             ", ispinor2=", ispinor2
             call wrtout(std_out,message,'COLL')
             do im1=1,ll*2+1
               write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
      &         (mat_out_c(im1,im2),im2=1,ll*2+1)
               call wrtout(std_out,message,'COLL')
             end do

           end do ! im
         end do ! ispinor
       end do ! isppol
       ABI_FREE(mat_out_c)
     endif

   end if !lpawu
 end do !atom

     !=====================================
     ! Deallocate gathermatlu
     !=====================================

 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     ABI_FREE(gathermatlu(iatom)%value)
     ABI_FREE(muspin(iatom)%value)
   end if
 end do
 ABI_FREE(gathermatlu)
 ABI_FREE(muspin)

 end subroutine magmomfspin_matlu

!!***


!!***
!!****f* m_matlu/magmomfzeeman_matlu
!! NAME
!! magmomfspin_matlu
!!
!! FUNCTION
!! return the product of occupation matrix of dimension [(2*ll+1)]**4 in the Ylm basis
!! with the matrix of Zeeman angular momentum element (L_u + 2*S_u) for the x,y and z direction.
!! Option gives the direction of the magnetic moment x==1, y==2 and z==3.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (FGendron)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! matlu1(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity in Ylm basis
!! natom :: number of atoms
!! option = 1 :: x axis
!!	  = 2 :: y axis
!!        = 3 :: z axis
!! optptr > 2 :: print Zeeman angular matrix elements and resulting product
!!
!! OUTPUT
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: product
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE
 subroutine magmomfzeeman_matlu(matlu,mu,natom,option,optprt)
 use defs_wvltypes

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,option,optprt
 complex(dp), allocatable, intent(inout) :: mu(:)
!arrays
 type(matlu_type), intent(inout) :: matlu(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im,ispinor,isppol,ispinor2
 integer :: lpawu,ll,jm,ml1,jc1,ms1,lcor,im1,im2,tndim,ndim
 character(len=500) :: message
 real(dp) :: xj
!arrays
 complex(dp),allocatable :: mat_out_c(:,:)
 integer, allocatable :: ind_msml(:,:)
 complex(dp), allocatable :: temp_mat(:,:)
 type(coeff2c_type), allocatable :: gathermatlu(:)
 type(coeff2c_type), allocatable :: muzeeman(:)
!************************************************************************

 !=====================================
 ! Allocate matrices
 !=====================================

 ABI_MALLOC(gathermatlu,(natom))
 ABI_MALLOC(muzeeman,(natom))
 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     tndim=2*(2*matlu(iatom)%lpawu+1)
     ABI_MALLOC(gathermatlu(iatom)%value,(tndim,tndim))
     gathermatlu(iatom)%value=czero
     ABI_MALLOC(muzeeman(iatom)%value,(tndim,tndim))
     muzeeman(iatom)%value=czero
   end if
 end do

 do iatom=1,natom
   lpawu=matlu(iatom)%lpawu
   if(lpawu.ne.-1) then
     ll=lpawu
     lcor=lpawu
     ABI_MALLOC(ind_msml,(2,-ll:ll))
     ind_msml=czero
     !=====================================
     !build Zeeman angular momentum matrix along x axis
     !=====================================
     if(option==1) then
       jc1=0
       do ms1=1,2
         do ml1=-ll,ll
          jc1=jc1+1
          ind_msml(ms1,ml1)=jc1
         end do
       end do

       jc1=0
       do ms1 =-1,1
         xj=float(ms1)+half
         do ml1 = -ll,ll
            jc1=jc1+1
            if(xj < 0.0 ) then
               muzeeman(iatom)%value(ind_msml(2,ml1),jc1) = 2*0.5
            else if(xj > 0.0) then
               muzeeman(iatom)%value(ind_msml(1,ml1),ind_msml(2,ml1)) = 2*0.5
            end if
            if(jc1 == 1) then
               muzeeman(iatom)%value(jc1,jc1+1) = sqrt(float((lcor*(lcor+1)) - ml1*(ml1 + 1)))*0.5
            endif
            if(jc1 > 1 .and.jc1 < 2*(2*ll+1)) then
               muzeeman(iatom)%value(jc1,jc1+1) = sqrt(float((lcor*(lcor+1)) - ml1*(ml1 + 1)))*0.5
               muzeeman(iatom)%value(jc1,jc1-1) = sqrt(float((lcor*(lcor+1)) - ml1*(ml1 - 1)))*0.5
            else if(jc1 == 2*(2*ll+1)) then
               muzeeman(iatom)%value(jc1,jc1-1) = sqrt(float((lcor*(lcor+1)) - ml1*(ml1 - 1)))*0.5
            end if
          end do
        end do

     !=====================================
     !build Zeeman angular momentum matrix along y axis
     !up spin is first
     !=====================================
     else if(option==2) then
        jc1=0
       do ms1=1,2
         do ml1=-ll,ll
          jc1=jc1+1
          ind_msml(ms1,ml1)=jc1
         end do
       end do

       jc1=0
       do ms1 =-1,1
         xj=float(ms1)+half
         do ml1 = -ll,ll
            jc1=jc1+1
            if(xj < 0.0 ) then
               muzeeman(iatom)%value(ind_msml(2,ml1),jc1) = 2*cmplx(zero,0.5,kind=dp)
            else if(xj > 0.0) then
               muzeeman(iatom)%value(ind_msml(1,ml1),ind_msml(2,ml1)) = 2*cmplx(zero,-0.5,kind=dp)
            end if
            if(jc1 == 1) then
               muzeeman(iatom)%value(jc1,jc1+1) = cmplx(zero,sqrt(float((lcor*(lcor+1)) - ml1*(ml1 + 1)))*0.5,kind=dp)
            endif
            if(jc1 > 1 .and. jc1 < 2*(2*ll+1)) then
               muzeeman(iatom)%value(jc1,jc1+1) = cmplx(zero,sqrt(float((lcor*(lcor+1)) - ml1*(ml1 + 1)))*0.5,kind=dp)
               muzeeman(iatom)%value(jc1,jc1-1) = cmplx(zero,-sqrt(float((lcor*(lcor+1)) - ml1*(ml1 - 1)))*0.5,kind=dp)
            else if(jc1== 2*(2*ll+1)) then
               muzeeman(iatom)%value(jc1,jc1-1) = cmplx(zero,-sqrt(float((lcor*(lcor+1)) - ml1*(ml1 - 1)))*0.5,kind=dp)
            end if
          end do
        end do


     !=====================================
     !build Zeeman angular momentum matrix along z axis
     !up spin is first
     !=====================================
     else if(option==3) then
       jc1=0
       do ms1=-1,1
         xj=float(ms1)+half
         do ml1=-ll,ll
            jc1=jc1+1
            if(jc1 < tndim+1) then
              if(xj < 0.0 ) then
                 muzeeman(iatom)%value(jc1,jc1) = ml1-2*xj
              else if(xj > 0.0) then
                 muzeeman(iatom)%value(jc1,jc1) = ml1-2*xj
              end if
            endif
         end do
       end do
     end if


     if(optprt>2) then
        write(message,'(a,i4)') "Zeeman angular momentum matrix elements in |m_l,m_s> basis for axis", option
        call wrtout(std_out,message,"COLL")
        do im=1,2*(ll*2+1)
          write(message,'(6(1x,9(1x,f4.1,",",f4.1)))') (muzeeman(iatom)%value(im,jm),jm=1,2*(ll*2+1))
          call wrtout(std_out,message,"COLL")
        end do
     end if

   ABI_FREE(ind_msml)
   end if !lpawu
 end do !atom

     !=====================================
     ! Reshape input Ylm matlu in one 14x14 matrix
     !=====================================

 !ABI_MALLOC(gathermatlu,(natom))
 !do iatom=1,natom
 !  if(matlu(iatom)%lpawu.ne.-1) then
 !    tndim=2*(2*matlu(iatom)%lpawu+1)
 !    ABI_MALLOC(gathermatlu(iatom)%value,(tndim,tndim))
 !    gathermatlu(iatom)%value=czero
 !  end if
 !end do

 call gather_matlu(matlu,gathermatlu,natom,option=1,prtopt=1)

!!printing for debug
!! write(std_out,*) "gathermatlu in magmomfspin"
!! do im1=1,tndim
!!   write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
!!        (gathermatlu(1)%value(im1,im2),im2=1,tndim)
!!   call wrtout(std_out,message,'coll')
!! end do

     !=====================================
     ! Matrix product of Occ and muzeeman
     !=====================================
 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     tndim=2*(2*matlu(iatom)%lpawu+1)
     ABI_MALLOC(temp_mat,(tndim,tndim))

     call zgemm('n','n',tndim,tndim,tndim,cone,gathermatlu(iatom)%value,tndim,muzeeman(iatom)%value,tndim,czero,temp_mat,tndim)

     gathermatlu(iatom)%value=temp_mat
     ABI_FREE(temp_mat)

     !!printing for debug
     !!write(std_out,*) "gathermatlu in magmomfspin after product"
     !!do im1=1,tndim
     !!   write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
     !!        (gathermatlu(1)%value(im1,im2),im2=1,tndim)
     !!   call wrtout(std_out,message,'coll')
     !!end do

     !=====================================
     ! Trace of matrix product
     !=====================================

   do im1=1,tndim
     do im2=1,tndim
       if(im1==im2) then
         mu(iatom) = mu(iatom) + gathermatlu(iatom)%value(im1,im2)
       end if
     end do
   end do

     !=====================================
     ! Reshape product matrix into matlu format
     !=====================================

    !call gather_matlu(matlu,gathermatlu,natom,option=-1,prtopt=1)


     !=====================================
     ! Print matlu
     !=====================================
     if(optprt>2) then
       ABI_MALLOC(mat_out_c,(2*ll+1,2*ll+1))
       ndim = 2*ll+1
       do isppol=1,matlu(1)%nsppol
         do ispinor=1,matlu(1)%nspinor
           do ispinor2=1,matlu(1)%nspinor
             mat_out_c(:,:) = matlu(iatom)%mat(1+(ispinor-1)*ndim:ndim+(ispinor-1)*ndim,1+(ispinor2-1)*ndim:ndim+(ispinor2-1)*ndim,isppol)

             write(message,'(2a, i2, a, i2, a, i2)') ch10,"Zeeman angular momentum matrix, isppol=", isppol, ", ispinor=",&
&            ispinor,", ispinor2=", ispinor2
             call wrtout(std_out,message,'COLL')
             do im1=1,ll*2+1
               write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
      &         (mat_out_c(im1,im2),im2=1,ll*2+1)
               call wrtout(std_out,message,'COLL')
             end do

           end do ! im
         end do ! ispinor
       end do ! isppol
       ABI_FREE(mat_out_c)
     endif ! optprt

   end if !lpawu
 end do !atom

     !=====================================
     ! Deallocate gathermatlu
     !=====================================

 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     ABI_FREE(gathermatlu(iatom)%value)
     ABI_FREE(muzeeman(iatom)%value)
   end if
 end do
 ABI_FREE(gathermatlu)
 ABI_FREE(muzeeman)
 end subroutine magmomfzeeman_matlu

!!***

!!***
!!****f* m_matlu/chi_matlu
!! NAME
!! chi_matlu
!!
!! FUNCTION
!! return the matrix of dimension [(2*ll+1)]**4 in the Ylm basis
!! with the matrix elements of the orbital (option=1), spin (option=2) and
!! total (option=3) angular momentum for z direction. It is used by the
!! QMC for the correlation function of the magnetic moment
!!
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (FGendron)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! matlu1(natom)%(nsppol,nspinor,nspinor,ndim,ndim) ::
!! natom :: number of atoms
!! option = 1 :: Orbital angular momentum along z axis
!!	  = 2 :: 2*Spin angluar momentum alonf z axis
!!        = 3 :: total angular momentum along z axis
!! optptr > 2 :: print angular matrix elements
!!
!! OUTPUT
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: quantity in Ylm basis
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE
 subroutine chi_matlu(matlu,natom,option,optprt)
 use defs_wvltypes

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,option,optprt
!arrays
 type(matlu_type), intent(inout) :: matlu(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im
 integer :: lpawu,ll,jm,ml1,jc1,ms1,lcor,tndim
 character(len=500) :: message
 real(dp) :: xj
!arrays
 integer, allocatable :: ind_msml(:,:)
! type(coeff2c_type), allocatable :: gathermatlu(:)
 type(coeff2c_type), allocatable :: muchi(:)
!************************************************************************

 !=====================================
 ! Allocate matrices
 !=====================================

 ABI_MALLOC(muchi,(natom))
 do iatom=1,natom
   !if(matlu(iatom)%lpawu.ne.-1) then
     tndim=2*(2*matlu(iatom)%lpawu+1)
     ABI_MALLOC(muchi(iatom)%value,(tndim,tndim))
     muchi(iatom)%value=czero
   !end if
 end do

 do iatom=1,natom
   lpawu=matlu(iatom)%lpawu
   if(lpawu.ne.-1) then
     ll=lpawu
     lcor=lpawu
     ABI_MALLOC(ind_msml,(2,-ll:ll))
     ind_msml=czero
     !=====================================
     !Orbital angular momentum matrix along z axis
     !=====================================
     if(option==1) then
       jc1=0
       do ms1=1,2
         do ml1=-ll,ll
          jc1=jc1+1
          ind_msml(ms1,ml1)=jc1
         end do
       end do

       jc1=0
       do ms1=-1,1
         do ml1=-ll,ll
            jc1=jc1+1
            if(jc1 <tndim+1) then
                muchi(iatom)%value(jc1,jc1) = ml1
            endif
          end do
        end do

     !=====================================
     !Spin angular momentum matrix along z axis
     !up spin is first
     !=====================================
     else if(option==2) then
        jc1=0
       do ms1=1,2
         do ml1=-ll,ll
          jc1=jc1+1
          ind_msml(ms1,ml1)=jc1
         end do
       end do

       jc1=0
       do ms1=-1,1
         xj=float(ms1)+half
         do ml1=-ll,ll
            jc1=jc1+1
            if(jc1<tndim+1) then
                if(xj < 0.0 ) then
                        muchi(iatom)%value(jc1,jc1) = -2*xj
                else if(xj > 0.0) then
                        muchi(iatom)%value(jc1,jc1) = -2*xj

                end if
             endif
         end do
       end do

     !=====================================
     !Total angular momentum matrix along z axis
     !up spin is first
     !=====================================
     else if(option==3) then
       jc1=0
       do ms1=-1,1
         xj=float(ms1)+half
         do ml1=-ll,ll
            jc1=jc1+1
            if(jc1<tndim+1) then
                if(xj < 0.0 ) then
                        muchi(iatom)%value(jc1,jc1) = ml1-2*xj
                else if(xj > 0.0) then
                        muchi(iatom)%value(jc1,jc1) = ml1-2*xj
                end if
            endif
         end do
       end do
     end if

     if(optprt>2) then
        if(option==1) then
          write(message,'(a)') "Orbital angular momentum matrix elements in |m_l,m_s> basis"
        else if(option==2) then
          write(message,'(a)') "Spin angular momentum matrix elements in |m_l,m_s> basis"
        else if(option==3) then
          write(message,'(a)') "Zeeman angular momentum matrix elements in |m_l,m_s> basis"
        end if
        call wrtout(std_out,message,"COLL")
        do im=1,2*(ll*2+1)
          write(message,'(6(1x,9(1x,f4.1,",",f4.1)))') (muchi(iatom)%value(im,jm),jm=1,2*(ll*2+1))
          call wrtout(std_out,message,"COLL")
        end do
     end if

   ABI_FREE(ind_msml)

     !=====================================
     ! Reshape matrix into matlu format
     !=====================================

   ! call gather_matlu(matlu,muchi(iatom),natom=1,option=-1,prtopt=1)


   end if !lpawu
 end do !atom

  !=====================================
  ! Reshape matrix into matlu format
  !=====================================

  call gather_matlu(matlu,muchi,natom,option=-1,prtopt=1)

     !=====================================
     ! Deallocate gathermatlu
     !=====================================

 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     ABI_FREE(muchi(iatom)%value)
   end if
 end do
 ABI_FREE(muchi)

 end subroutine chi_matlu

!!***

!!****f* m_matlu/trace_prod_matlu
!! NAME
!! trace_prod_matlu
!!
!! FUNCTION
!! Computes Tr(matlu1*matlu2) for each atom. It is NOT assumed
!! that either matlu1 or matlu2 is symmetric, so this routine is
!! suboptimal if this is the case.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu1(natom),matlu2(natom) :: input quantity
!!  natom=number of atoms in cell.
!!  opt_add = if present, add the result to trace
!!  trace_tot = if present, computes the trace over all atoms
!!  iatom = if present, only computes the contribution of iatom
!!          if <=0, add all contributions as usual
!!
!! OUTPUT
!!  trace(natom) :: Tr(matlu1*matlu2) for each atom
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine trace_prod_matlu(matlu1,matlu2,natom,trace,trace_tot,iatom)

!Arguments ------------------------------------
 integer, intent(in) :: natom
 type(matlu_type), intent(in) :: matlu1(natom),matlu2(natom)
 complex(dp), intent(inout) :: trace(natom)
 complex(dp), optional, intent(out) :: trace_tot
 integer, optional, intent(in) :: iatom
!Local variables-------------------------------
 integer :: ia1,ia2,iatom_,isppol,lpawu,nspinor,nsppol
!************************************************************************

 nspinor = matlu1(1)%nspinor
 nsppol  = matlu1(1)%nsppol

 trace(:)  = czero

 ia1 = 1 ; ia2 = natom
 if (present(iatom)) then
   if (iatom > 0) then
     ia1 = iatom ; ia2 = iatom
   end if
 end if ! present(iatom)

 do iatom_=ia1,ia2
   lpawu = matlu1(iatom_)%lpawu
   if (lpawu == -1) cycle
   do isppol=1,nsppol
     trace(iatom_) = trace(iatom_) + sum(matlu1(iatom_)%mat(:,:,isppol)*transpose(matlu2(iatom_)%mat(:,:,isppol)))
   end do ! isppol
   if (nsppol == 1 .and. nspinor == 1) trace(iatom_) = trace(iatom_) * two
 end do ! iatom

 if (present(trace_tot)) trace_tot = sum(trace(:))

 end subroutine trace_prod_matlu
!!***

!!****f* m_matlu/xmpi_matlu
!! NAME
!! xmpi_sum_matlu
!!
!! FUNCTION
!!  Put matlu into a buffer and perform the required MPI operation
!!  (xmpi_sum or xmpi_bcast) on the input communicator.
!!
!! INPUTS
!!  matlu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  natom = number of atoms
!!  comm = communicator on which to perform the required communication
!!  master = master node (only used in the case of xmpi_bcast, default is 0)
!!  option = 1 (default) : xmpi_sum
!!         = 2 : xmpi_bcast
!!
!! OUTPUT
!!
!! SOURCE

 subroutine xmpi_matlu(matlu,natom,comm,master,option)

!Arguments ------------------------------------
 integer, intent(in) :: comm,natom
 type(matlu_type), intent(inout) :: matlu(natom)
 integer, optional, intent(in) :: master,option
!Local variables-------------------------------
 integer :: iatom,ibuf,im1,ierr,isppol,lpawu
 integer :: master_node,ndim,nspinor,nsppol,opt,siz_buf
 complex(dp), allocatable :: buffer(:)
!************************************************************************

 nspinor = matlu(1)%nspinor
 nsppol  = matlu(1)%nsppol
 opt = 1
 if (present(option)) opt = option
 master_node = 0
 if (present(master)) master_node = master

 siz_buf = 0
 do iatom=1,natom
   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   siz_buf = siz_buf + (2*lpawu+1)**2
 end do ! iatom

 siz_buf = siz_buf * (nspinor**2) * nsppol

 ABI_MALLOC(buffer,(siz_buf))

 ibuf = 0
 do iatom=1,natom
   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim = nspinor * (2*lpawu+1)
   do isppol=1,nsppol
     do im1=1,ndim
       buffer(ibuf+1:ibuf+ndim) = matlu(iatom)%mat(:,im1,isppol)
       ibuf = ibuf + ndim
     end do ! im1
   end do ! isppol
 end do ! iatom

 if (opt == 1) then
   call xmpi_sum(buffer(:),comm,ierr)
 else if (opt == 2) then
   call xmpi_bcast(buffer(:),master_node,comm,ierr)
 end if ! opt

 ibuf = 0
 do iatom=1,natom
   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim = nspinor * (2*lpawu+1)
   do isppol=1,nsppol
     do im1=1,ndim
       matlu(iatom)%mat(:,im1,isppol) = buffer(ibuf+1:ibuf+ndim)
       ibuf = ibuf + ndim
     end do ! im1
   end do ! isppol
 end do ! iatom

 ABI_FREE(buffer)

 end subroutine xmpi_matlu
!!***

!!****f* m_matlu/symmetrize_matlu
!! NAME
!! symmetrize_matlu
!!
!! FUNCTION
!!  Symmetrizes matlu (A = (A+A^T)/2
!!
!! INPUTS
!!  matlu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  natom = number of atoms
!!
!! OUTPUT
!!  err = difference between A and symmetrized A
!!
!! SOURCE

 subroutine symmetrize_matlu(matlu,natom,err)

!Arguments ------------------------------------
 integer, intent(in) :: natom
 type(matlu_type), intent(inout) :: matlu(natom)
 real(dp), optional, intent(out) :: err
!Local variables-------------------------------
 integer :: iatom,isppol,lpawu,nspinor,nsppol,tndim
 real(dp) :: err_
 complex(dp), allocatable :: mat_tmp(:,:)
!************************************************************************

 nspinor = matlu(1)%nspinor
 nsppol  = matlu(1)%nsppol

 if (present(err)) err = zero

 do iatom=1,natom
   lpawu = matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   tndim = nspinor * (2*lpawu+1)
   ABI_MALLOC(mat_tmp,(tndim,tndim))
   do isppol=1,nsppol
     mat_tmp(:,:) = half * (matlu(iatom)%mat(:,:,isppol)+ &
                  & transpose(matlu(iatom)%mat(:,:,isppol)))
     if (present(err)) then
       err_ = sum(abs(mat_tmp(:,:)-matlu(iatom)%mat(:,:,isppol)))
       if (err_ > err) err = err_
     end if ! present(err)
     matlu(iatom)%mat(:,:,isppol) = mat_tmp(:,:)
   end do ! isppol
   ABI_FREE(mat_tmp)
 end do ! iatom

 end subroutine symmetrize_matlu
!!***

!!****f* m_matlu/ylm2jmj_matlu
!! NAME
!! ylm2jmj_matlu
!!
!! FUNCTION
!! Transform mat from Ylm to JmJ basis or vice versa
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom) :: input quantity
!!  natom :: number of atoms
!!  option=1 go from Ylm to JmJ basis
!!  option=2 go from JmJ to Ylm basis
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE

 subroutine ylm2jmj_matlu(matlu,natom,option,paw_dmft)

!Arguments ------------------------------------
 integer, intent(in) :: natom,option
 type(matlu_type), intent(inout) :: matlu(natom)
 type(paw_dmft_type), intent(in) :: paw_dmft
!Local variables-------------------------------
 integer :: iatom,lpawu,nspinor,tndim,tndim_max
 complex(dp), allocatable :: mat_tmp(:,:)
 character(len=1) :: c1,c2
!************************************************************************

 nspinor = paw_dmft%nspinor
 if (nspinor == 1) ABI_BUG("nspinor should be equal to 2")

 tndim_max = nspinor * (2*paw_dmft%maxlpawu+1)

 if (option == 1) then
   c1 = "c" ; c2 = "n"
 else
   c1 = "n" ; c2 = "c"
 end if

 ABI_MALLOC(mat_tmp,(tndim_max,tndim_max))

 do iatom=1,natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   if (lpawu == 0) ABI_BUG("l should not be equal to 0")
   tndim = nspinor * (2*lpawu+1)

   call abi_xgemm("n",c2,tndim,tndim,tndim,cone,matlu(iatom)%mat(:,:,1),tndim, &
                & paw_dmft%jmj2ylm(:,1:tndim,lpawu+1),tndim_max,czero,mat_tmp(:,1:tndim),tndim_max)

   call abi_xgemm(c1,"n",tndim,tndim,tndim,cone,paw_dmft%jmj2ylm(:,1:tndim,lpawu+1), &
                & tndim_max,mat_tmp(:,1:tndim),tndim_max,czero,matlu(iatom)%mat(:,:,1),tndim)

 end do ! iatom

 ABI_FREE(mat_tmp)

 end subroutine ylm2jmj_matlu
!!***

!!***
!!****f* m_matlu/magnfield_matlu
!! NAME
!! magnfield_matlu
!!
!! FUNCTION
!! return the matrix of magnetic moment mz times Bz
!!
!!
!! COPYRIGHT
!! Copyright (C) 2005-2026 ABINIT group (FGendron)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! matlu1(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity in Ylm basis
!! natom :: number of atoms
!! bfield :: value of magnetic field in Tesla
!! option = 1 :: scalar spin angular momentum along z axis
!! option = 2 :: SOC total angular momentum (L+2S) along z axis
!!
!! OUTPUT
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: product
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! SOURCE
 subroutine magnfield_matlu(matlu,natom,bfield,option)
 use defs_basis
 use defs_wvltypes
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,option
 real(dp) :: bfield
!arrays
 type(matlu_type), intent(inout) :: matlu(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im,ndim,isppol
 integer :: ll,ml1,jc1,ms1,tndim
 real(dp) :: xj
!arrays
 type(coeff2c_type), allocatable :: magnmatb(:)
!************************************************************************

 !================================
 ! Allocate matrices
 !================================

 ABI_MALLOC(magnmatb,(natom))
 do iatom=1,natom
   if(matlu(iatom)%lpawu .ne. -1) then
     tndim=2*(2*matlu(iatom)%lpawu+1)
     ABI_MALLOC(magnmatb(iatom)%value,(tndim,tndim))
     magnmatb(iatom)%value=czero
   endif
 enddo

 if(option .eq. 1) then

 !================================
 ! Scalar magnetism (Spin only case)
 ! H = mu_B*g_e*S_Z*B_z
 !================================

   do iatom=1,natom
     if(matlu(iatom)%lpawu .ne. -1) then
       ndim=2*matlu(iatom)%lpawu+1
       do isppol=1,matlu(iatom)%nsppol
         do im=1,ndim
           if (isppol .eq. 1) then
             matlu(iatom)%mat(im,im,isppol) = half*bfield
           else
             matlu(iatom)%mat(im,im,isppol) = -half*bfield
           endif
         enddo ! im
       enddo ! isppol
     endif ! lpawu
   enddo ! natom


 elseif(option .eq. 2) then

 !================================
 ! Spin-orbit magnetism
 ! H = mu_B*(L_z+g_e*S_Z)*B_z
 !================================

   do iatom=1,natom
     if(matlu(iatom)%lpawu .ne. -1) then
       tndim=2*(2*matlu(iatom)%lpawu+1)
       ll=matlu(iatom)%lpawu

       jc1=0
       do ms1=-1,1
         xj=float(ms1)+half
         do ml1=-ll,ll
           jc1=jc1+1
           if(jc1 < tndim+1) then
             if (xj < 0.0) then
               magnmatb(iatom)%value(jc1,jc1) = half*(ml1-2*xj)*bfield
             elseif(xj > 0.0) then
               magnmatb(iatom)%value(jc1,jc1) = half*(ml1-2*xj)*bfield
             endif
           endif
         enddo !ml1
       enddo ! ms1
     endif !lpawu
   enddo !natom
 endif !option

 !=======================
 ! reshape matrix
 !=======================

 if(option .eq. 2) then
   call gather_matlu(matlu,magnmatb(natom),natom,option=-1,prtopt=1)
 endif

 !================================
 ! Deallocate matrices
 !================================

 do iatom=1,natom
   if(matlu(iatom)%lpawu .ne. -1) then
     ABI_FREE(magnmatb(iatom)%value)
   endif
 enddo

 ABI_FREE(magnmatb)

 end subroutine magnfield_matlu
!!***

!!****f* m_matlu/magmomjmj_matlu                                                                         
!! NAME                                                                                                  
!! magmomjmj_matlu                                                                                       
!!                                                                                                       
!! FUNCTION                                                                                              
!! return the matrix of magnetic moments in the Jmj basis                  
!!                                                                                                       
!!                                                                                                       
!! COPYRIGHT                                                                                             
!! Copyright (C) 2005-2026 ABINIT group (FGendron)                                                       
!! This file is distributed under the terms of the                                                       
!! GNU General Public License, see ~abinit/COPYING                                                       
!! or http://www.gnu.org/copyleft/gpl.txt .                                                              
!!                                                                                                       
!! INPUTS                                                                                                
!!                                                                                                       
!! OUTPUT                                                                                                
!!                                                                                                       
!! SIDE EFFECTS                                                                                          
!!                                                                                                       
!! NOTES                                                                                                  
!!                                                                                                        
!! SOURCE                                                                                                 
 subroutine magmomjmj_matlu(matlu,natom)                                                                  
 use defs_basis                                                                                           
 use defs_wvltypes                                                                                        
 implicit none                                                                                            
                                                                                                          
!Arguments ------------------------------------                                                           
!scalars                                                                                                  
 integer, intent(in) :: natom                                                                             
!arrays                                                                                                   
 type(matlu_type), intent(inout) :: matlu(natom)                                                          
!Local variables-------------------------------                                                           
!scalars                                                                                                  
 integer :: iatom,lpawu,ll,ml1,ms1,jm,jc1,tndim,jj                                                        
 real(dp) :: xj,xmj                                                                                       
!arrays                                                                                                   
 integer, allocatable :: ind_msml(:,:)                                                                    
 type(coeff2c_type), allocatable :: gathermatlu(:)                                                        
 complex(dpc),allocatable :: mlms2jmj(:,:)                                                                
!************************************************************************                                 
                                                                                                          
 !=====================================                                                                   
 ! Allocate Matrices                                                                                      
 !=====================================                                                                   
                                                                                                          
 ABI_MALLOC(gathermatlu,(natom))                                                                          
                                                                                                          
 do iatom=1,natom                                                                                         
   lpawu=matlu(iatom)%lpawu                                                                               
   if(lpawu.ne.-1) then                                                                                   
     ll=lpawu                                                                                             
     tndim=2*(2*ll+1)                                                                                     
                                                                                                          
     ABI_MALLOC(gathermatlu(iatom)%value,(tndim,tndim))                                                   
     gathermatlu(iatom)%value=czero                                                                       
     ABI_MALLOC(mlms2jmj,(tndim,tndim))                                                                   
     mlms2jmj=czero                                                                                       
     ABI_MALLOC(ind_msml,(2,-ll:ll))                                                                      
     mlms2jmj=czero                                                                                       
                                                                                                          
 !=====================================                                                                   
 ! Build J,M_J matrix                                                                                     
 !=====================================                                                                   
                                                                                                          
    jc1=0                                                                                                 
    do ms1=1,2                                                                                            
      do ml1=-ll,ll                                                                                       
        jc1=jc1+1                                                                                         
        ind_msml(ms1,ml1)=jc1                                                                             
      end do                                                                                              
    end do                                                                                                
                                                                                                          
    jc1=0                                                                                                 
    do jj=ll,ll+1                                                                                         
      xj=float(jj)-half !  xj is in {ll-0.5, ll+0.5}                                                      
      do jm=-jj,jj-1                                                                                      
        xmj=float(jm)+half  ! xmj is in {-xj,xj}                                                          
        jc1=jc1+1           ! Global index for JMJ                                                        
        if(nint(xj+0.5)==ll+1) then  ! if xj=ll+0.5                                                       
          mlms2jmj(jc1,jc1)=xmj   !  J=L+0.5 and m_J=L+0.5                                                
        else if(nint(xj-0.5)==ll-1) then                                                                  
          mlms2jmj(jc1,jc1)=xmj   !  J=L+0.5 and m_J=-L-0.5                                               
        end if                                                                                            
      end do                                                                                              
    end do                                                                                                
                                                                                                          
    !print to debug                                                                    
    !write(message,'(3a)') ch10,"JMJ Matrix"                                           
    !call wrtout(std_out,message,"COLL")                                               
    !do im=1,2*(ll*2+1)                                                                
    !  write(message,'(12(1x,18(1x,f5.2,f5.2)))') (mlms2jmj(im,jm),jm=1,2*(ll*2+1))    
    !  call wrtout(std_out,message,"COLL")                                             
    !end do                                                                            
                                                                                        
  !=====================================                                                
  ! Put back into matlu format                                                          
  !=====================================                                                
                                                                                        
   gathermatlu(iatom)%value=mlms2jmj                                                  
                                                                                       
   call gather_matlu(matlu,gathermatlu(iatom),natom=1,option=-1,prtopt=0)             
                                                                                        
  !=====================================                                                
  ! Deallocate Matrices                                                                 
  !=====================================                                                
                                                                                        
   ABI_FREE(gathermatlu(iatom)%value)                                                
    end if !lpawu                                                                       
  end do !natom                                                                         

ABI_FREE(mlms2jmj) 
ABI_FREE(ind_msml)
ABI_FREE(gathermatlu)                                                                 
                                                                                        
end subroutine magmomjmj_matlu                                                        


END MODULE m_matlu
!!***
