!{\src2tex{textfont=tt}}
!!****f* ABINIT/fftpac
!! NAME
!! fftpac
!!
!! FUNCTION
!! Allow for data copying to modify the stride (dimensioning) of a three-
!! dimensional array, for more efficient three dimensional fft.
!! NOTE that the arrays are in REAL space.
!!
!! Note that arrays aa and bb may be the same array (start at the same address).
!! The array aa also incorporate a spin variable.
!! MG FIXME: THIS IS **VERY BAD** AS FORTRAN DOES NOT ALLOW FOR ALIASING
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, MF, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ispden=actual spin-density of interest
!!  nspden=number of spin-density components
!!  n1,n2,n3=actual data dimensions, dimensions of complex array a
!!  nd1,nd2,nd3=array dimensions of (larger) array b
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  option= see description of side effects
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  aa & bb arrays are treated as input or output depending on option: 
!!  option=1  aa(n1*n2*n3,ispden) <-- bb(nd1,nd2,nd3) real case
!!  option=2  aa(n1*n2*n3,ispden) --> bb(nd1,nd2,nd3) real case
!!  option=10 aa(n1*n2*n3,ispden) <-- bb(nd1,nd2,nd3) complex case like option 1 real part
!!  option=11 aa(n1*n2*n3,ispden) <-- bb(nd1,nd2,nd3) complex case like option 1 imag part
!!
!! PARENTS
!!      dfpt_mkrho,dfpt_nstpaw,dfpt_rhofermi,dfpt_vtorho,dfptnl_resp,energy
!!      fock_getghc,getgh1c,gwls_hamiltonian,ks_ddiago,m_epjdos,m_io_kss,mkrho
!!      suscep_stat,vtorho
!!
!! CHILDREN
!!      ptabs_fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fftpac(ispden,mpi_enreg,nspden,n1,n2,n3,nd1,nd2,nd3,ngfft,aa,bb,option)

 use defs_basis
 use m_errors

 use defs_abitypes,  only : MPI_type
 use m_mpinfo,       only : ptabs_fourdp

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftpac'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ispden,n1,n2,n3,nd1,nd2,nd3,nspden,option
 type(mpi_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(inout) :: aa(n1*n2*n3/ngfft(10),nspden),bb(nd1,nd2,nd3)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,index,me_fft,nproc_fft
 character(len=500) :: message
 !arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)

! *************************************************************************

 me_fft=ngfft(11); nproc_fft=ngfft(10)

 if (option==1.or.option==2) then
   if (nd1<n1.or.nd2<n2.or.nd3<n3) then
     write(message,'(a,3i0,2a,3i0,a)')&
&     'Each of nd1,nd2,nd3=',nd1,nd2,nd3,ch10,&
&     'must be >=      n1, n2, n3 =',n1,n2,n3,'.'
     MSG_BUG(message)
   end if
 else
   if (2*nd1<n1.or.nd2<n2.or.nd3<n3) then
     write(message,'(a,3i0,2a,3i0,a)')&
&     'Each of 2*nd1,nd2,nd3=',2*nd1,nd2,nd3,ch10,&
&     'must be >= (n1, n2, n3) =',n1,n2,n3,'.'
     MSG_BUG(message)
   end if
 end if

 ! Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 if (option==1) then
   do i3=1,n3
     if (me_fft==fftn3_distrib(i3)) then
       do i2=1,n2
         do i1=1,n1
           aa(i1+n1*(i2-1+n2*(ffti3_local(i3)-1)),ispden)=bb(i1,i2,i3)
         end do
       end do
     end if
   end do

 else if (option==2) then
   !  Here we avoid corrupting the data in a while writing to b in the
   !  case in which a and b are same array.
   !  Also: replace "trash" data with 0 s to avoid floating point
   !  exceptions when this data is actually manipulated in fft.
   do i3=nd3,n3+1,-1
     do i2=nd2,1,-1
       do i1=nd1,1,-1
         bb(i1,i2,i3)=0.d0
       end do
     end do
   end do
   do i3=n3,1,-1
     if (me_fft==fftn3_distrib(i3)) then
       do i2=nd2,n2+1,-1
         do i1=nd1,1,-1
           bb(i1,i2,i3)=0.d0
         end do
       end do
       do i2=n2,1,-1
         do i1=nd1,n1+1,-1
           bb(i1,i2,i3)=0.d0
         end do
         do i1=n1,1,-1
           bb(i1,i2,i3)=aa(i1+n1*(i2-1+n2*(ffti3_local(i3) - 1)),ispden)
         end do
       end do
     end if
   end do
!  MF
 else if (option==10 .or. option==11) then
   index=1
   if(option==11) index=2
   do i3=1,n3
     do i2=1,n2
       do i1=1,n1/2
         aa(index,ispden)=bb(i1,i2,i3)
         index=index+2
       end do
     end do
   end do
!  MF
 else
   write(message,'(a,i0,a)')' Bad option =',option,'.'
   MSG_BUG(message)
 end if

end subroutine fftpac
!!***
