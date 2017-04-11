!{\src2tex{textfont=tt}}
!!****f* ABINIT/multipoles_fftr
!! NAME
!! multipoles_fftr
!!
!! FUNCTION
!!  Compute spatial multipole moments of input array on FFT grid
!!  Namely, the electrical dipole, quadrupole, etc... of the electron density
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT group (MJV, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  arraysp(nfft,nspden)=the array whose average has to be computed
!!  [distribfft<type(distribfft_type)>]= -optional- infos related to FFT parallelism
!!  [mpi_comm_grid]= -optional- MPI communicator over the grid
!!  origin(3)=vector defining the origin of the dipole (point of observation, in reduced coordinates)
!!  nfft=number of FFT points stored by one proc
!!  nfftot=total number of FFT points
!!  ngfft =number of subdivisions along each lattice vector
!!  nspden=number of spin-density components
!!  rprimd = dimensionful lattice vectors
!!
!! OUTPUT
!!  dipole(nspden)=mean value of the dipole of input array, for each nspden component
!!
!! PARENTS
!!      multipoles_fftr
!!
!! CHILDREN
!!      multipoles_fftr,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine multipoles_fftr(arraysp,dipole,nfft,ngfft,nspden,rprimd,origin,&
&                          distribfft,mpi_comm_grid)

 use m_profiling_abi
 use defs_basis
 use m_errors
 use defs_abitypes
 use m_distribfft
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'multipoles_fftr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nspden
 integer,intent(in),optional :: mpi_comm_grid
 type(distribfft_type),intent(in),optional,target :: distribfft
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: arraysp(nfft,nspden),origin(3),rprimd(3,3)
 real(dp),intent(out) :: dipole(3,nspden)
!Local variables-------------------------------
!scalars
 integer,parameter :: ishift=5
 integer :: ierr,ifft1,ifft2,ifft3,ifft0,ifft,ispden,i1,i2,i3,n1,n2,n3
 integer :: me_fft,my_mpi_comm,nfftot
 logical :: fftgrid_found
 real(dp) :: invn1,invn2,invn3,invnfftot,wrapfft1,wrapfft2,wrapfft3
 character(len=500) :: msg
 type(distribfft_type),pointer :: my_distribfft
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: dipole_tmp(3,nspden)

! *************************************************************************


!Several initializations
 my_mpi_comm=xmpi_comm_self;if (present(mpi_comm_grid)) my_mpi_comm=mpi_comm_grid
 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3)
 nfftot=n1*n2*n3;invnfftot=one/(dble(nfftot))
 invn1=one/real(n1,kind=dp);invn2=one/real(n2,kind=dp);invn3=one/real(n3,kind=dp)
 me_fft=xmpi_comm_rank(my_mpi_comm)

 dipole(:,:)=zero

!Get the distrib associated with the FFT grid
 if (present(distribfft)) then
   my_distribfft => distribfft
 else
   ABI_DATATYPE_ALLOCATE(my_distribfft,)
   call init_distribfft_seq(my_distribfft,'f',n2,n3,'fourdp')
 end if
 fftgrid_found=.false.
 if (n2 == my_distribfft%n2_coarse ) then
   if (n3 == size(my_distribfft%tab_fftdp3_distrib)) then
     fftn3_distrib => my_distribfft%tab_fftdp3_distrib
     ffti3_local => my_distribfft%tab_fftdp3_local
     fftgrid_found=.true.
   end if
 end if
 if (n2 == my_distribfft%n2_fine ) then
   if (n3 == size(my_distribfft%tab_fftdp3dg_distrib)) then
     fftn3_distrib => my_distribfft%tab_fftdp3dg_distrib
     ffti3_local => my_distribfft%tab_fftdp3dg_local
     fftgrid_found = .true.
   end if
 end if
 if (.not.(fftgrid_found)) then
   msg='Unable to find an allocated distrib for the FFT grid!'
   MSG_BUG(msg)
 end if

!Loop over FFT grid points
!$OMP PARALLEL PRIVATE(ifft,i1,i2,ifft1,ifft2,ispden,wrapfft1,wrapfft2)
 do ifft3=1,n3
   i3=mod(ifft3-1+ishift*n3,n3)
   if(fftn3_distrib(1+i3)==me_fft) then
     wrapfft3=mod(real(ifft3-1,kind=dp)*invn3-origin(3)+1.5_dp,one)-half
     ifft0=1+n1*n2*(ffti3_local(1+i3)-1)
!$OMP SINGLE
     dipole_tmp=zero
!$OMP END SINGLE
!$OMP DO COLLAPSE(2), REDUCTION(+:dipole_tmp)
     do ifft2=1,n2
       do ifft1=1,n1
         i2=mod(ifft2-1+ishift*n2,n2)
         i1=mod(ifft1-1+ishift*n1,n1)
         wrapfft2=mod(real(ifft2-1,kind=dp)*invn2-origin(2)+1.5_dp,one)-half
         wrapfft1=mod(real(ifft1-1,kind=dp)*invn1-origin(1)+1.5_dp,one)-half
         ifft=ifft0+n1*i2+i1

!        Computation of integral(s)
         do ispden=1,nspden
           dipole_tmp(1,ispden)=dipole_tmp(1,ispden)+wrapfft1*arraysp(ifft,ispden)
           dipole_tmp(2,ispden)=dipole_tmp(2,ispden)+wrapfft2*arraysp(ifft,ispden)
           dipole_tmp(3,ispden)=dipole_tmp(3,ispden)+wrapfft3*arraysp(ifft,ispden)
         end do

       end do
     end do
!$OMP END DO
!$OMP SINGLE
     dipole=dipole+dipole_tmp
!$OMP END SINGLE
   end if
 end do
!$OMP END PARALLEL

!MPI parallelization
 if (xmpi_comm_size(my_mpi_comm)>1) then
   call xmpi_sum(dipole,my_mpi_comm,ierr)
 end if

!From reduced to cartesian coordinates
 do ispden=1,nspden
   dipole(:,ispden)=matmul(rprimd,dipole(:,ispden))*invnfftot
 end do

 if (.not.present(distribfft)) then
   call destroy_distribfft(my_distribfft)
   ABI_DATATYPE_DEALLOCATE(my_distribfft)
 end if

end subroutine multipoles_fftr
!!***
