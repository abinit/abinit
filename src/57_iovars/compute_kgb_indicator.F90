!{\src2tex{textfont=tt}}
!!****f* ABINIT/compute_kgb_indicator
!! NAME
!! compute_kgb_indicator
!!
!! FUNCTION
!! Only for "KGB" parallelism (LOBPCG algorithm for Ground-state):
!!  Give an indicator of performance for a given distribution of processors
!!  (npband, npfft and bandpp).
!!  Determine best choice of parameters for Scalapack and/or Magma Linear Algebra routines.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (FD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  bandpp=internal lobpcg optimization variable
!!  glb_comm=communicator for global MPI communications
!!  mband=maximum number of bands.
!!  mband=maximum number of plane waves
!!  npband=number of processor 'band'
!!  npfft = number of processor 'fft'
!!  uselinalggpu=indicate if we also test the gpu linear algebra
!!
!! OUTPUT
!!  acc_kgb = indicator of performance
!!  npslk = number of process to used in communicators
!!
!! SIDE EFFECTS
!! This routine can be used to find an indicator in order to refine automatic process distribution.
!!   This indicator is returned in acc_kgb
!! This routine can be used to find the optimal values of np_slk parameter (ScaLapack)
!!   and wheter or not we should use Magma for Linear Algebra in lobpcgwf
!!
!! PARENTS
!!      finddistrproc
!!
!! CHILDREN
!!      abi_linalg_finalize,abi_linalg_init,abi_xhegv,abi_xorthonormalize
!!      wrtout,xmpi_bcast,xmpi_comm_free
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine compute_kgb_indicator(acc_kgb,bandpp,glb_comm,mband,mpw,npband,npfft,npslk,&
&                                uselinalggpu)

 use defs_basis
 use m_abi_linalg
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_time,  only : abi_wtime

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'compute_kgb_indicator'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bandpp,glb_comm,mband,mpw,npband,npfft
 integer,intent(inout) :: npslk,uselinalggpu
 real(dp),intent(inout) :: acc_kgb
!arrays

!Local variables-------------------------------
!scalars
 integer,parameter :: max_number_of_npslk=10,max_number_of_iter=10
 integer :: blocksize,bigorder,ierr,ii,islk,islk1,iter,jj,keep_gpu
 integer :: kgb_comm,my_rank,np_slk,np_slk_max,np_slk_best,nranks
 integer :: use_lapack_gpu,use_slk,vectsize
 real(dp) :: min_eigen,min_ortho,time_xeigen,time_xortho
 character(len=500) :: message
!arrays
 integer,allocatable :: ranks(:),val_npslk(:)
 real(dp),allocatable :: eigen(:),grama(:,:),gramb(:,:)
 complex(dpc),allocatable :: blockvectorbx(:,:),blockvectorx(:,:),sqgram(:,:)

!******************************************************************

 DBG_ENTER("COLL")

#ifdef DEBUG_MODE
 write(message,'(a,3i3)') 'compute_kgb_indicator : (bpp,npb,npf) = ', bandpp, npband, npfft
 call wrtout(std_out,message,'PERS')
#endif

!Create local communicator for test
 if (xmpi_paral==1) then
   nranks=npfft*npband
   ABI_ALLOCATE(ranks,(nranks))
   ranks=(/((my_rank-1),my_rank=1,nranks)/)
   kgb_comm=xmpi_subcomm(glb_comm,nranks,ranks,my_rank_in_group=my_rank)
   ABI_DEALLOCATE(ranks)
 else
   kgb_comm=xmpi_comm_self
   my_rank=0
 end if

!Only for process in the new subgroup
 if (my_rank/=xmpi_undefined) then

!  We enforce vectsize >=blocksize  (This is not true in lobpcg but
!  these are rare cases and this simplify the matrix constructions below...)
   blocksize=npband*bandpp
   vectsize=max(1+mpw/(npband*npfft),blocksize)
   bigorder=3*blocksize

   ABI_ALLOCATE(blockvectorx,(vectsize,blocksize))
   ABI_ALLOCATE(blockvectorbx,(vectsize,blocksize))
   ABI_ALLOCATE(sqgram,(blocksize,blocksize))
   ABI_ALLOCATE(grama,(2*bigorder,bigorder))
   ABI_ALLOCATE(gramb,(2*bigorder,bigorder))
   ABI_ALLOCATE(eigen,(bigorder))
   ABI_ALLOCATE(val_npslk,(max_number_of_npslk)) ! not too much values tested

   min_eigen=greatest_real
   min_ortho=greatest_real
   np_slk_best=-1 ; np_slk_max=0
#ifdef HAVE_LINALG_SCALAPACK
   np_slk_max=min(max_number_of_npslk,npband*npfft)
#endif

!  Preselect a range of available np_slk values
   val_npslk(1)=0 ; val_npslk(2)=1
   do islk=3,np_slk_max
     np_slk=val_npslk(islk-1)*2
     do while ((modulo(npband*npfft,np_slk)>0).and.(np_slk<(npband*npfft)))
       np_slk=np_slk+1
     end do
     if(np_slk>(npband*npfft)) exit
     val_npslk(islk)=np_slk
   end do
   np_slk_max=islk-1

!  Loop over np_slk values
   islk1=1
#ifdef HAVE_LINALG_MAGMA
   islk1=1-uselinalggpu
#endif
   do islk=islk1,np_slk_max

     time_xortho=zero ; time_xeigen=zero

     use_slk=0
     if (islk==0) then
!      This is the test for the GPU
       use_lapack_gpu=1 ; np_slk=0
     else
       use_lapack_gpu=0 ; np_slk=val_npslk(islk)
       if (np_slk>0) use_slk=1
     end if

!    Initialize linalg parameters for this np_slk value
!    For the first np_slk value, everything is initialized
!    For the following np_slk values, only Scalapack parameters are updated
     call abi_linalg_init(kgb_comm,np_slk,bigorder,my_rank,only_scalapack=(islk>islk1))

!    We could do mband/blocksize iter as in lobpcg but it's too long
     do iter=1,max_number_of_iter

!      Build matrixes
       do ii=1,vectsize
         do jj=1,blocksize
           if (ii>jj) then
             blockvectorx(ii,jj) =czero
             blockvectorbx(ii,jj)=czero
           else
             blockvectorx(ii,jj) =cone
             blockvectorbx(ii,jj)=cone
           end if
         end do
       end do
       grama=zero;gramb=zero
       do jj=1,bigorder
         do ii=jj,bigorder
           if (ii==jj) then
             grama(2*ii-1,jj)=one
             gramb(2*ii-1,jj)=one
           else
             grama(2*ii-1:2*ii,jj)=one
             grama(2*jj-1,ii)= one
             grama(2*jj  ,ii)=-one
           end if
         end do
       end do

!      Call to abi_xorthonormalize
       time_xortho=time_xortho-abi_wtime()
       call abi_xorthonormalize(blockvectorx,blockvectorbx,blocksize,kgb_comm,sqgram,vectsize)
       time_xortho = time_xortho + abi_wtime()

!      Call to abi_xhegv
       time_xeigen=time_xeigen-abi_wtime()
       call abi_xhegv(1,'v','u',bigorder,grama,gramb,eigen,&
&       x_cplx=2,use_slk=use_slk,use_gpu=use_lapack_gpu)
       time_xeigen=time_xeigen+abi_wtime()

     end do ! iter

!    Finalize linalg parameters for this np_slk value
!    For the last np_slk value, everything is finalized
!    For the previous np_slk values, only Scalapack parameters are updated
     call abi_linalg_finalize(only_scalapack=(islk<np_slk_max))

     time_xortho= time_xortho*mband/blocksize
     time_xeigen= time_xeigen*mband/blocksize
     if (time_xortho<min_ortho) min_ortho=time_xortho
     if (time_xeigen<min_eigen) then
       min_eigen=time_xeigen
       np_slk_best=np_slk
       keep_gpu=use_lapack_gpu
     end if

   end do ! np_slk

#ifdef DEBUG_MODE
   write(message,'(2(a,es15.3),a,i3)') ' In the best case, xortho took ',min_ortho,&
&   ' and xeigen took ',min_eigen,' for np_slk=',np_slk_best
   call wrtout(std_out,message,'PERS')
#endif

!  Final values to be sent to others process
   acc_kgb=min_ortho+four*min_eigen
   npslk=np_slk_best
   uselinalggpu=keep_gpu

   ABI_DEALLOCATE(blockvectorx)
   ABI_DEALLOCATE(blockvectorbx)
   ABI_DEALLOCATE(sqgram)
   ABI_DEALLOCATE(grama)
   ABI_DEALLOCATE(gramb)
   ABI_DEALLOCATE(eigen)
   ABI_DEALLOCATE(val_npslk)

 end if ! my_rank in group

!Free local MPI communicator
 call xmpi_comm_free(kgb_comm)

!Broadcast of results to be sure every process has them
 call xmpi_bcast(acc_kgb,0,glb_comm,ierr)
 call xmpi_bcast(npslk,0,glb_comm,ierr)
 call xmpi_bcast(uselinalggpu,0,glb_comm,ierr)

#ifndef DEBUG_MODE
 ABI_UNUSED(message)
#endif

 DBG_EXIT("COLL")

end subroutine compute_kgb_indicator
!!***
