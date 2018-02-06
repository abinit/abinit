!{\src2tex{textfont=tt}}
!!****f* ABINIT/zerosym
!! NAME
!! zerosym
!!
!! FUNCTION
!! Symmetrize an array on the FFT grid by vanishing some term on the boundaries.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (GZ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cplex= if 1, input array is REAL, if 2, input array is COMPLEX
!!  mpi_enreg=informations about MPI parallelization
!!  n1,n2,n3=FFT dimensions nfft=n1*n2*n3
!!  ig1,ig2,ig3=optional arguments= indexes of unbalanced g-vectors to cancel
!!              if not present, ig1=1+n1/2, ig2=1+n2/2, ig3=1+n3/2 for even n1,n2,n3
!!              if igj=-1, nothing is done in direction j
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  array(cplex,n1*n2*n3)=complex array to be symetrized
!!
!! PARENTS
!!      atm2fft,dfpt_atm2fft,dfpt_dyfro,forces,hartre,initro,m_fock,pawmknhat
!!      pawmknhat_psipsi,prcref,prcref_PMA,stress,transgrid
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine zerosym(array,cplex,n1,n2,n3,&
&                  ig1,ig2,ig3,comm_fft,distribfft) ! Optional arguments

 use defs_basis
 use m_profiling_abi
 use m_distribfft
 use m_errors
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'zerosym'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,n1,n2,n3
 integer,optional,intent(in) :: ig1,ig2,ig3,comm_fft
 type(distribfft_type),intent(in),target,optional :: distribfft
!arrays
 real(dp),intent(inout) :: array(cplex,n1*n2*n3)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,ifft,ifft_proc,index,j,j1,j2,j3,me_fft,nd2
 integer :: nproc_fft,n1sel,nn12,n2sel,n3sel,r2
 !arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)

! **********************************************************************

 DBG_ENTER("COLL")

 me_fft=0;nproc_fft=1
 if (present(comm_fft)) then
   me_fft=xmpi_comm_rank(comm_fft)
   nproc_fft=xmpi_comm_size(comm_fft)
 end if
 nd2=(n2-1)/nproc_fft+1
 nn12=n1*n2

!Get the distrib associated with this fft_grid 
 if (present(distribfft)) then
   if (n2== distribfft%n2_coarse) then
     fftn2_distrib => distribfft%tab_fftdp2_distrib
     ffti2_local => distribfft%tab_fftdp2_local
   else if(n2 == distribfft%n2_fine) then
     fftn2_distrib => distribfft%tab_fftdp2dg_distrib
     ffti2_local => distribfft%tab_fftdp2dg_local
   else
     MSG_BUG("Unable to find an allocated distrib for this fft grid")
   end if
 else
   ABI_ALLOCATE(fftn2_distrib,(n2))
   ABI_ALLOCATE(ffti2_local,(n2))
   fftn2_distrib=0;ffti2_local=(/(i2,i2=1,n2)/)
 end if

 if (present(ig1)) then
   n1sel=ig1
 else if (mod(n1,2)==0) then
   n1sel=1+n1/2
 else
   n1sel=-1
 end if
 if (present(ig2)) then
   n2sel=ig2
 else if (mod(n2,2)==0) then
   n2sel=1+n2/2
 else
   n2sel=-1
 end if
 if (present(ig3)) then
   n3sel=ig3
 else if (mod(n3,2)==0) then
   n3sel=1+n3/2
 else
   n3sel=-1
 end if

 if (n1sel>0) then
   index=n1sel-nn12-n1
   do i3=1,n3
     index=index+nn12;ifft=index
     do i2=1,n2
       ifft=ifft+n1
       if (nproc_fft>1) then
!        MPIWF: consider ifft only if it is treated by the current proc and compute its adress
         j=ifft-1;j1=modulo(j,n1);j2=modulo(j/n1,n2);j3=j/(n1*n2) !;r2=modulo(j2,nd2)
         if(fftn2_distrib(j2+1)==me_fft) then ! MPIWF this ifft is to be treated by me_fft
           r2= ffti2_local(j2+1) - 1 
           ifft_proc=n1*(nd2*j3+r2)+j1+1 !this is ifft in the current proc
           array(:,ifft_proc)=zero
         end if
       else
         array(:,ifft)=zero
       end if
     end do
   end do
 end if

 if (n2sel>0) then
   index=n1*n2sel-nn12-n1
   do i3=1,n3
     index=index+nn12;ifft=index
     do i1=1,n1
       ifft=ifft+1
       if (nproc_fft>1) then
!        MPIWF: consider ifft only if it is treated by the current proc and compute its adress
         j=ifft-1;j1=modulo(j,n1);j2=modulo(j/n1,n2);j3=j/(n1*n2);
         if(fftn2_distrib(j2+1)==me_fft) then ! MPIWF this ifft is to be treated by me_fft
           r2= ffti2_local(j2+1) - 1 
           ifft_proc=n1*(nd2*j3+r2)+j1+1 !this is ifft in the current proc
           array(:,ifft_proc)=zero
         end if
       else
         array(:,ifft)=zero
       end if
     end do
   end do
 end if

 if (n3sel>0) then
   index=nn12*n3sel-nn12-n1
   do i2=1,n2
     index=index+n1;ifft=index
     do i1=1,n1
       ifft=ifft+1
       if (nproc_fft>1) then
!        MPIWF: consider ifft only if it is treated by the current proc and compute its adress
         j=ifft-1;j1=modulo(j,n1);j2=modulo(j/n1,n2);j3=j/(n1*n2)
         if(fftn2_distrib(j2+1)==me_fft) then ! MPIWF this ifft is to be treated by me_fft
           r2= ffti2_local(j2+1) - 1 
           ifft_proc=n1*(nd2*j3+r2)+j1+1 !this is ifft in the current proc
           array(:,ifft_proc)=zero
         end if
       else
         array(:,ifft)=zero
       end if
     end do
   end do
 end if

 if (.not.present(distribfft)) then
   ABI_DEALLOCATE(fftn2_distrib)
   ABI_DEALLOCATE(ffti2_local)
 end if

 DBG_EXIT("COLL")

end subroutine zerosym
!!***
