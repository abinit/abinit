!{\src2tex{textfont=tt}}
!!****f* ABINIT/kgindex
!! NAME
!! kgindex
!!
!! FUNCTION
!! Compute the index of each plane wave on a FFT grid.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  kg_k(3,npw_k)=dimensionless coords of G vecs (integer)
!!  mpi_enreg=informations about MPI parallelization
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  npw_k=number of planewaves
!!
!! OUTPUT
!!  indpw_k(npw_k)=linear list number (in fft box) of given G vector for the current processor (local adress)
!!                =0 if kg_k(ipw) is not treated by this procesor
!!  mask(npw_k)=True if  kg_k(ipw) belongs to this procesor, false otherwise.
!!
!! NOTES
!!   mpi_enreg is not necessary in this case (the info is also in ngfft), but much more easy to read...
!!
!! PARENTS
!!      m_fft_prof,m_gsphere,m_screening,m_shirley,m_wfd,prcref,prcref_PMA
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine kgindex(indpw_k,kg_k,mask,mpi_enreg,ngfft,npw_k)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kgindex'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: kg_k(3,npw_k),ngfft(18)
 integer,intent(out) :: indpw_k(npw_k)
 logical,intent(out) :: mask(npw_k)
!Local variables-------------------------------
!scalars
 integer :: ig,ig1,ig2,ig3,me_fft,n1,n2,n3,nd2
 character(len=500) :: msg
 !arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 !integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)

! *************************************************************************

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)

!Use the following indexing (N means ngfft of the adequate direction)
!0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= kg
!1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index

 me_fft=mpi_enreg%me_fft
 nd2=(n2-1)/mpi_enreg%nproc_fft+1

 if (n2== mpi_enreg%distribfft%n2_coarse) then
   fftn2_distrib => mpi_enreg%distribfft%tab_fftdp2_distrib
   ffti2_local => mpi_enreg%distribfft%tab_fftdp2_local
 else if (n2 == mpi_enreg%distribfft%n2_fine) then
   fftn2_distrib => mpi_enreg%distribfft%tab_fftdp2dg_distrib
   ffti2_local => mpi_enreg%distribfft%tab_fftdp2dg_local
 else
   MSG_BUG("Unable to find an allocated distrib for this fft grid")
 end if

 !call ptabs_fourwf(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 do ig=1,npw_k
   ig1=modulo(kg_k(1,ig),n1)
   ig2=modulo(kg_k(2,ig),n2)
   ig3=modulo(kg_k(3,ig),n3)
   if(me_fft==fftn2_distrib(ig2+1)) then
     ig2=ffti2_local(ig2+1) - 1
     indpw_k(ig)=ig1+1+n1*(ig2+nd2*ig3)
     mask(ig)=.true.
   else
     indpw_k(ig)=0
     mask(ig)=.false.
   end if
   if ( ANY(kg_k(:,ig)>ngfft(1:3)/2) .or. ANY(kg_k(:,ig)<-(ngfft(1:3)-1)/2) ) then
     write(msg,'(a,i0,a)')" The G-vector with ig: ",ig," falls outside the FFT box."
     MSG_ERROR(msg)
   end if
 end do

end subroutine kgindex
!!***
