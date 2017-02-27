!{\src2tex{textfont=tt}}
!!****f* ABINIT/multipoles_fftr
!! NAME
!! multipoles_fftr
!!
!! FUNCTION
!!  Compute spatial multipole moments of input array on FFT grid
!!  Namely, the electrical dipole, quadrupole, etc... of the electron density
!!  call mean_fftr to deal with the averaging over several MPI OMP processors
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT group (MJV, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  arraysp(nfft,nspden)=the array whose average has to be computed
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

subroutine multipoles_fftr(arraysp,dipole,mpi_enreg,nfft,ngfft,nspden,rprimd,neworigin)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_cgtools

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'multipoles_fftr'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,ngfft(3),nspden
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(in) :: neworigin(3)
 real(dp),intent(in) :: arraysp(nfft,nspden)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(out) :: dipole(3,nspden)

!Local variables-------------------------------
!scalars
 integer :: ifft1, ifft2, ifft3, ifft, ispden, nfftot
 character(len=500) :: message
 real(dp) :: invn1, invn2, invn3
!arrays
 real(dp) :: meansp(nspden)
 real(dp),allocatable :: wrapfft(:)
 real(dp), allocatable :: tmpsp(:,:)

! *************************************************************************

 ABI_ALLOCATE(tmpsp,(nfft,nspden))

 nfftot = ngfft(1)*ngfft(2)*ngfft(3)
 invn1 = one / ngfft(1)
 invn2 = one / ngfft(2)
 invn3 = one / ngfft(3)


!for the moment impose no fft parallelization
!FIXME: needs to know somehow which fft points (i1i2i3 triplets) are on this processor...
 if (nfft /= nfftot) then 
   write (message,'(3a)') 'Error: fft parallelization of multipoles_fftr is not coded yet.',ch10,&
&   ' return from routine and continue as if nothing happened'
   call wrtout(std_out, message, 'COLL')
   return
 end if

!multiply value at each point by xred1
!FIXME: make this more efficient than looping explicitly over n2 n3
!also check if ordering is correct wrt fft grid: z cycles fastest?
 ABI_ALLOCATE(wrapfft,(ngfft(1)))
 do ifft1 = 1, ngfft(1)
   wrapfft(ifft1) = mod((ifft1-1-neworigin(1)+half*ngfft(1)), dble(ngfft(1))) - half*ngfft(1)
 end do
 ifft = 1
 do ifft3 = 1, ngfft(3)
   do ifft2 = 1, ngfft(2)
     do ifft1 = 1, ngfft(1)
       tmpsp(ifft,:) = arraysp(ifft,:) * wrapfft(ifft1) * invn1
       ifft = ifft + 1
     end do
   end do
 end do
 ABI_DEALLOCATE(wrapfft)

!average over cell
 call mean_fftr(tmpsp,meansp,nfft,nfftot,nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
 dipole(1,:) = meansp

!multiply value at each point by xred2
 ABI_ALLOCATE(wrapfft,(ngfft(2)))
 do ifft2 = 1, ngfft(2)
   wrapfft(ifft2) = mod((ifft2-1-neworigin(2)+half*ngfft(2)), dble(ngfft(2))) - half*ngfft(2)
 end do
 ifft = 1
 do ifft3 = 1, ngfft(3)
   do ifft2 = 1, ngfft(2)
     do ifft1 = 1, ngfft(1)
       tmpsp(ifft,:) = arraysp(ifft,:) * wrapfft(ifft2) * invn2
       ifft = ifft + 1
     end do
   end do
 end do
 ABI_DEALLOCATE(wrapfft)

!average over cell
 call mean_fftr(tmpsp,meansp,nfft,nfftot,nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
 dipole(2,:) = meansp

!multiply value at each point by xred3
 ABI_ALLOCATE(wrapfft,(ngfft(3)))
 do ifft3 = 1, ngfft(3)
   wrapfft(ifft3) = mod((ifft3-1-neworigin(3)+half*ngfft(3)), dble(ngfft(3))) - half*ngfft(3)
 end do
 ifft = 1
 do ifft3 = 1, ngfft(3)
   do ifft2 = 1, ngfft(2)
     do ifft1 = 1, ngfft(1)
       tmpsp(ifft,:) = arraysp(ifft,:) * wrapfft(ifft3) * invn3
       ifft = ifft + 1
     end do
   end do
 end do
 ABI_DEALLOCATE(wrapfft)

!average over cell
 call mean_fftr(tmpsp,meansp,nfft,nfftot,nspden)
 dipole(3,:) = meansp

!Turn xred to cartesian space for output
 do ispden = 1, nspden
   dipole(:,ispden) = matmul(rprimd,dipole(:,ispden))
 end do


 ABI_DEALLOCATE(tmpsp)

end subroutine multipoles_fftr
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/multipoles_out
!! NAME
!! multipoles_out
!!
!! FUNCTION
!!  Output multipole moments of input array on FFT grid, calculated with multipoles_fftr
!!  Namely, the electrical dipole, quadrupole, etc... of the electron density
!!
!! COPYRIGHT
!! Copyright (C) 2010-2017 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  arraysp(nfft,nspden)=the array whose average has to be computed
!!  mpi_enreg=information about MPI parallelization
!!  nfft=number of FFT points stored by one proc
!!  nfftot=total number of FFT points
!!  ngfft =number of subdivisions along each lattice vector
!!  nspden=number of spin-density components
!!  rprimd = dimensionful lattice vectors
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  print to str_out
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine multipoles_out(arraysp,mpi_enreg,natom,nfft,ngfft,nspden,&
&  ntypat,rprimd,typat,ucvol,xred,ziontypat)

 use m_profiling_abi
 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'multipoles_out'
 use interfaces_14_hidewrite
 use interfaces_53_spacepar, except_this_one => multipoles_out
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom, ntypat
 integer,intent(in) :: nfft,ngfft(3),nspden
 type(MPI_type),intent(in) :: mpi_enreg
 real(dp), intent(in) :: ucvol
!arrays
 integer, intent(in) :: typat(natom)
 real(dp),intent(in) :: arraysp(nfft,nspden)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(in) :: ziontypat(ntypat)

! local vars
 integer :: iatom, ispden
 character(len=500) :: message

 real(dp) :: center_of_charge(3)
 real(dp), allocatable :: dipole_el(:,:)
 real(dp) :: dipole_ions(3), ziontotal, dipole_tot(3)

! *************************************************************************

!nuclear part of dipole
 dipole_ions = zero
 ziontotal = zero
 do iatom = 1, natom
   dipole_ions = dipole_ions + xred(:,iatom)*ziontypat(typat(iatom))
   ziontotal = ziontotal + ziontypat(typat(iatom))
 end do
 write (message, '(a,3E18.6,3a,3E18.6,2a)') ' Ionic dipole =        ',  dipole_ions, ' (a.u.) or ', ch10, &
& '              =        ', dipole_ions/dipole_moment_debye, ' (D)',ch10
 call wrtout(std_out, message, 'COLL')

!find coordinates of center of charge on fft grid
!NOTE: wrt center of charge, dipole_ions is 0
 center_of_charge(:) = dipole_ions(:)/ziontotal * ngfft(:)
 write (message, '(a,3E20.10)') ' Center of charge for ionic distribution (red coordinates): ', dipole_ions(:)/ziontotal
 call wrtout(std_out, message, 'COLL')

!get electronic part of dipole with respect to center of charge of ions
 ABI_ALLOCATE(dipole_el,(3,nspden))

 call multipoles_fftr(arraysp,dipole_el,mpi_enreg,nfft,ngfft,nspden,rprimd,center_of_charge)
 dipole_el = dipole_el * ucvol

 dipole_tot(1) = -sum(dipole_el(1,:))
 dipole_tot(2) = -sum(dipole_el(2,:))
 dipole_tot(3) = -sum(dipole_el(3,:))

!output
 write (message, '(2a)') ch10,' Dipole in cartesian coord  with respect to the center of ionic charge'
 call wrtout(std_out, message, 'COLL')
 write (message, '(a)') ' Electronic part : '
 call wrtout(std_out, message, 'COLL')
 do ispden = 1, nspden
   write (message, '(a,I5,a,3E18.6,3a,3E18.6,a)') '   density component ', ispden, ' dipole = ',  dipole_el(:,ispden),&
&   ' (a.u.) ', ch10, '                                  = ', dipole_el(:,ispden)/dipole_moment_debye, ' (D)' 
   call wrtout(std_out, message, 'COLL')
 end do

 write (message, '(a,3E18.6,3a,3E18.6,2a)') ' Total dipole =        ',  dipole_tot, ' (a.u.) or ', ch10, &
& '              =        ', dipole_tot/dipole_moment_debye, ' (D)',ch10
 call wrtout(std_out, message, 'COLL')

 ABI_DEALLOCATE(dipole_el)
 
end subroutine multipoles_out
!!***
