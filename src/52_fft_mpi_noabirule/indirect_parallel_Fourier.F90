!{\src2tex{textfont=tt}}
!!****f* ABINIT/indirect_parallel_Fourier
!! NAME
!! indirect_parallel_Fourier
!!
!! FUNCTION
!! The purpose of this routine is to transfer data from right to left
!! right(:,index(i))=left(:,i)
!! The difficulty is that right and left are distributed among processors
!! We will suppose that the distribution is done as a density in Fourier space
!! We first order the right hand side data according to the processor
!! in which they are going to be located in the left hand side.
!! This is done is a way such that  a mpi_alltoall put the data on the correct processor.
!! We also transfer their future adress.
!! A final ordering put everything in place
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (GZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  index(sizeindex)= global adress for the transfer from right to left
!!  left(2,nleft)=left hand side
!!  mpi_enreg=information about MPI parallelization
!!  ngleft(18)=contain all needed information about 3D FFT for the left hand side
!!  see ~abinit/doc/variables/vargs.htm#ngfft
!!  ngright(18)=contain all needed information about 3D FFT for the right hand side
!!  see ~abinit/doc/variables/vargs.htm#ngfft
!!  nleft=second dimension of left array (for this processor)
!!  nright=second dimension of right array (for this processor)
!!  sizeindex=size of the index array (different form nright, because it is global to all proccessors)
!!
!! OUTPUT
!!  left(2,nleft)=the elements of the right hand side, at the correct palce in the correct processor
!!
!! NOTES
!!  A lot of things to improve.
!!
!! PARENTS
!!      prcref,prcref_PMA,transgrid
!!
!! CHILDREN
!!      mpi_alltoall,ptabs_fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine indirect_parallel_Fourier(index,left,mpi_enreg,ngleft,ngright,nleft,nright,paral_kgb,right,sizeindex)

 use m_profiling_abi
 use m_errors
 use defs_basis
 use defs_abitypes

 use m_mpinfo,   only : ptabs_fourdp

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'indirect_parallel_Fourier'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ngleft(18),ngright(18),nleft,nright,paral_kgb,sizeindex
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: index(sizeindex)
 real(dp),intent(in) :: right(2,nright)
 real(dp),intent(inout) :: left(2,nleft)

!Local variables ---------------------------------------
!scalars
 integer :: ierr,i_global,ileft,iright,iright_global
 integer :: j,j1,j2,j3,j_global,jleft_global
 integer :: jleft_local,me_fft,n1l,n2l,n3l,n1r,n2r,n3r,nd2l,nd2r
 integer :: nproc_fft,proc_dest,r2,siz_slice_max
!arrays
 integer,allocatable :: index_recv(:),index_send(:),siz_slice(:), ffti2r_global(:)
 integer, ABI_CONTIGUOUS pointer :: fftn2l_distrib(:),ffti2l_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3l_distrib(:),ffti3l_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn2r_distrib(:),ffti2r_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3r_distrib(:),ffti3r_local(:)
 real(dp),allocatable :: right_send(:,:),right_recv(:,:)

! *************************************************************************
 n1r=ngright(1);n2r=ngright(2);n3r=ngright(3)
 n1l=ngleft(1) ;n2l=ngleft(2) ;n3l=ngleft(3)
 nproc_fft=mpi_enreg%nproc_fft; me_fft=mpi_enreg%me_fft
 nd2r=n2r/nproc_fft; nd2l=n2l/nproc_fft

 !Get the distrib associated with the left fft_grid 
 call ptabs_fourdp(mpi_enreg,n2l,n3l,fftn2l_distrib,ffti2l_local,fftn3l_distrib,ffti3l_local)

 !Get the distrib associated with the right fft_grid 
 call ptabs_fourdp(mpi_enreg,n2r,n3r,fftn2r_distrib,ffti2r_local,fftn3r_distrib,ffti3r_local)

 !Precompute local --> global corespondance
 ABI_ALLOCATE(ffti2r_global,(nd2r))
 ffti2r_global(:) = -1
 do j2=1,n2r
    if( fftn2r_distrib(j2) == me_fft ) then 
       ffti2r_global( ffti2r_local(j2) ) = j2 
    end if
 end do
 

 ABI_ALLOCATE(siz_slice,(nproc_fft))
 siz_slice(:)=0
 do i_global=1,sizeindex !look for the maximal size of slice of data
  j_global=index(i_global)!; write(std_out,*) j_global,i_global
  if(j_global /=0) then
    !use the fact that (j-1)=i1 + n1l*(j2l-1 + n2l*(j3l-1))
   proc_dest= fftn2l_distrib( modulo((j_global-1)/n1l,n2l) + 1) 
   siz_slice(proc_dest+1)=siz_slice(proc_dest+1)+1
!DEBUG
!write(std_out,*) 'in indirect proc',proc_dest,siz_slice(proc_dest+1)
!ENDDEBUG
  end if
 end do
 siz_slice_max=maxval(siz_slice) !This value could be made smaller by looking locally
!and performing a allgather with a max
!DEBUG
!write(std_out,*) 'siz_slice,sizeindex,siz_slice',siz_slice(:),sizeindex,siz_slice_max
!write(std_out,*) 'sizeindex,nright,nleft',sizeindex,nright,nleft
!ENDDEBUG
 ABI_ALLOCATE(right_send,(2,nproc_fft*siz_slice_max))
 ABI_ALLOCATE(index_send,(nproc_fft*siz_slice_max))
 siz_slice(:)=0; index_send(:)=0; right_send(:,:)=zero
 do iright=1,nright
  j=iright-1;j1=modulo(j,n1r);j2=modulo(j/n1r,nd2r);j3=j/(n1r*nd2r)
  j2 = ffti2r_global(j2+1) - 1 
  iright_global=n1r*(n2r*j3+j2)+j1+1
  jleft_global=index(iright_global)
  if(jleft_global/=0)then
     j=jleft_global-1;j1=modulo(j,n1l);j2=modulo(j/n1l,n2l);j3=j/(n1l*n2l); r2=ffti2l_local(j2+1)-1
   jleft_local=n1l*(nd2l*j3+r2)+j1+1
   proc_dest=fftn2l_distrib(j2+1) 
   siz_slice(proc_dest+1)=siz_slice(proc_dest+1)+1
   right_send(:,proc_dest*siz_slice_max+siz_slice(proc_dest+1))=right(:,iright)
   index_send(proc_dest*siz_slice_max+siz_slice(proc_dest+1))=jleft_local
!DEBUG
!   write(std_out,*) 'loop ir',jleft_local,jleft_global,iright_global,iright
!ENDDEBUG
  end if
 end do
 ABI_ALLOCATE(right_recv,(2,nproc_fft*siz_slice_max))
 ABI_ALLOCATE(index_recv,(nproc_fft*siz_slice_max))
#if defined HAVE_MPI
  if(paral_kgb == 1) then
    call mpi_alltoall (right_send,2*siz_slice_max, &
&                          MPI_double_precision, &
&                          right_recv,2*siz_slice_max, &
&                          MPI_double_precision,mpi_enreg%comm_fft,ierr)
    call mpi_alltoall (index_send,siz_slice_max, &
&                          MPI_integer, &
&                          index_recv,siz_slice_max, &
&                          MPI_integer,mpi_enreg%comm_fft,ierr)
  endif
#endif
 do ileft=1,siz_slice_max*nproc_fft
!DEBUG
!write(std_out,*)index_recv(ileft)
!ENDEBUG
 if(index_recv(ileft) /=0 ) left(:,index_recv(ileft))=right_recv(:,ileft)
 end do
 ABI_DEALLOCATE(right_recv)
 ABI_DEALLOCATE(index_recv)
 ABI_DEALLOCATE(right_send)
 ABI_DEALLOCATE(index_send)
 ABI_DEALLOCATE(siz_slice)
 ABI_DEALLOCATE(ffti2r_global)

end subroutine indirect_parallel_Fourier
!!***
