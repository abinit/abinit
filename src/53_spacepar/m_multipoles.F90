!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_multipoles
!! NAME
!!  m_multipoles
!!
!! FUNCTION
!!  Compute spatial multipole moments of input array on FFT grid
!!
!! COPYRIGHT
!!  Copyright (C) 2003-2018 ABINIT group (MJV, MT)
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

module m_multipoles

 use defs_basis
 use m_errors
 use m_profiling_abi
 use defs_abitypes
 use m_distribfft
 use m_xmpi


 implicit none

 private
!!***

 public :: multipoles_out    ! Output multipole moments of input array on FFT grid, calculated with multipoles_fftr
!!***

contains
!!***

!!****f* m_multipoles/multipoles_fftr
!! NAME
!! multipoles_fftr
!!
!! FUNCTION
!!  Compute spatial multipole moments of input array on FFT grid
!!  Namely, the electrical dipole, quadrupole, etc... of the electron density
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
!!      multipoles_out
!!
!! CHILDREN
!!      destroy_distribfft,init_distribfft_seq,xmpi_sum
!!
!! SOURCE

subroutine multipoles_fftr(arraysp,dipole,nfft,ngfft,nspden,rprimd,origin,&
&                          distribfft,mpi_comm_grid)


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
!$OMP DO COLLAPSE(2) REDUCTION(+:dipole_tmp)
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

!!****f* m_multipoles/multipoles_out
!! NAME
!! multipoles_out
!!
!! FUNCTION
!!  Output multipole moments of input array on FFT grid, calculated with multipoles_fftr
!!  Namely, the electrical dipole, quadrupole, etc... of the electron density
!!
!! INPUTS
!!  rhor(nfft,nspden)=electronic density
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms
!!  nfft=number of FFT points stored by one proc
!!  ngfft =number of subdivisions along each lattice vector
!!  nspden=number of spin-density components
!!  ntypat=number of atom types
!!  rprimd(3,3)=dimensionful lattice vectors
!!  typat(ntypat)=types of atoms
!!  ucvol=unit cell volume
!!  unit_out=file unit to print out
!!  ziontypat(ntypat)=ionic charge of each type of atom
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      multipoles_fftr,wrtout
!!
!! SOURCE

subroutine multipoles_out(rhor,mpi_enreg,natom,nfft,ngfft,nspden,&
&                         ntypat,rprimd,typat,ucvol,unit_out,xred,ziontypat)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'multipoles_out'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,ntypat,unit_out
 real(dp), intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer, intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3),xred(3,natom),ziontypat(ntypat)
!Local variables ------------------------------
!scalars
 integer :: iatom,nspden_updn
 real(dp) :: ziontotal
 character(len=500) :: message
!arrays
 real(dp) :: center_of_charge(3),dipole_el(3,2),dipole_ions_cart(3)
 real(dp) :: dipole_ions_red(3),dipole_tot(3),tmp(3)

! *************************************************************************

!Separate spins only for nspden=2
 nspden_updn=merge(1,2,nspden/=2)

!Title

!Get nuclear part of dipole
 dipole_ions_red(:) = zero ; ziontotal = zero
 do iatom = 1, natom
   dipole_ions_red(:) = dipole_ions_red(:) + xred(:,iatom)*ziontypat(typat(iatom))
   ziontotal = ziontotal + ziontypat(typat(iatom))
 end do
 dipole_ions_cart(:) = matmul(rprimd,dipole_ions_red(:))

!Find coordinates of center of charge on FFT grid
 center_of_charge(1:3) = dipole_ions_red(1:3)/ziontotal

!Get electronic part of dipole with respect to center of charge of ions (in cart. coord.)
 dipole_el = zero
 call multipoles_fftr(rhor(:,1:nspden_updn),dipole_el(:,1:nspden_updn),nfft,ngfft,nspden_updn,&
& rprimd,center_of_charge,distribfft=mpi_enreg%distribfft,mpi_comm_grid=mpi_enreg%comm_fft)
 dipole_el(1:3,1:nspden_updn)=dipole_el(1:3,1:nspden_updn)*ucvol

!Take into account storage of rhor (up+dn,up)
 if (nspden==2) then
   tmp(1:3)=dipole_el(1:3,1)
   dipole_el(1:3,1)=dipole_el(1:3,2)
   dipole_el(1:3,2)=tmp(1:3)-dipole_el(1:3,2)
 end if

!Compute total dipole
!NOTE: wrt center of charge, dipole_ions is 0
 dipole_tot(1) = - sum(dipole_el(1,1:nspden_updn))
 dipole_tot(2) = - sum(dipole_el(2,1:nspden_updn))
 dipole_tot(3) = - sum(dipole_el(3,1:nspden_updn))

!Output
 write (message, '(2a)') ch10,' ----- Electric nuclear dipole wrt the center of ionic charge ----- '
 call wrtout(unit_out, message, 'COLL')
 write (message, '(a,3(1x,ES12.5))') &
& ' Center of charge for ionic distribution (red. coord.): ',center_of_charge(1:3)
 call wrtout(unit_out, message, 'COLL')
 write (message, '(3a,3(1x,E16.6),3a,3(1x,E16.6),a)') ' -----',ch10,&
& ' Ionic dipole (cart. coord.)     = ',dipole_ions_cart, ' (a.u.)',ch10, &
& '                                 = ',dipole_ions_cart/dipole_moment_debye,' (D)'
 call wrtout(unit_out, message, 'COLL')
 if (nspden/=2) then
   !This is compatible with nspden=4
   write (message, '(3a,3(1x,E16.6),3a,3(1x,E16.6),a)') ' -----',ch10,&
&   ' Electronic dipole (cart. coord.)= ',dipole_el(:,1),' (a.u.)',ch10,&
&   '                                 = ',dipole_el(:,1)/dipole_moment_debye,' (D)'
 else
   write (message, '(3a,3(1x,E16.6),a,3(2a,3(1x,E16.6),a))') ' -----',ch10,&
&   ' Electronic dipole (cart. coord.)= ',dipole_el(:,1),' up (a.u.)',ch10,&
&   '                                 = ',dipole_el(:,2),' dn (a.u.)',ch10,&
&   '                                 = ',dipole_el(:,1)/dipole_moment_debye,' up (D)',ch10,&
&   '                                 = ',dipole_el(:,2)/dipole_moment_debye,' dn (D)'
 end if
 call wrtout(unit_out, message, 'COLL')
 write (message, '(3a,3(1x,E16.6),3a,3(1x,E16.6),a)') ' -----',ch10,&
& ' Total dipole (cart. coord.)     = ',dipole_tot,' (a.u.)',ch10,&
& '                                 = ',dipole_tot/dipole_moment_debye,' (D)'
 call wrtout(unit_out, message, 'COLL')
 call wrtout(unit_out, ' ', 'COLL')

end subroutine multipoles_out
!!***

end module m_multipoles
!!***
