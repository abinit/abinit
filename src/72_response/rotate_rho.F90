!{\src2tex{textfont=tt}}
!!****f* ABINIT/rotate_rho
!! NAME
!! rotate_rho
!!
!! FUNCTION
!! rotate density in real and reciprocal space 
!!
!! COPYRIGHT
!! Copyright (C) 2013-2018 ABINIT group (MVerst)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!!  ngfft=array of dimensions for different FFT grids
!!  nspden=number of spin-density components
!!  rhor1(cplex*nfft,nspden)=array for Fourier transform of RF electron density
!!  === if psps%usepaw==1 TODO: extend to PAW
!!    pawrhoij1(natom) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!  symrel1=single symmetry operation in real space to apply to rho
!!  tnon = eventual translation associated to symrel1
!!
!! OUTPUT
!!  rhog1_eq(2,nfft)= symmetric density in reciprocal space for equivalent perturbation
!!  rhor1_eq(cplex*nfft,nspden) = symmetric density in real space for equivalent perturbation
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!      fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine rotate_rho(cplex, itirev, mpi_enreg, nfft, ngfft, nspden, &
&   rhor1, rhog1_eq, rhor1_eq, symrel1, tnon)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rotate_rho'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!args
 integer,intent(in) :: cplex, nfft, nspden, itirev
 integer,intent(in) :: ngfft(18)

 integer, intent(in) :: symrel1(3,3)
 real(dp),intent(in) :: tnon(3)
 real(dp),intent(inout) :: rhor1(cplex*nfft,nspden)

 real(dp),intent(out) :: rhog1_eq(2,nfft)
 real(dp),intent(out) :: rhor1_eq(cplex*nfft,nspden)

 type(MPI_type),intent(in) :: mpi_enreg

! local vars
 integer :: id1, id2, id3
 integer :: n1, n2, n3, nd2
 integer :: l1, l2, l3
 integer :: i1, i2, i3, ind1, ind2
 integer :: j1, j2, j3
 integer :: k1, k2, k3
 integer :: nproc_fft, ispden, me_fft
 real(dp) :: arg
 logical :: t_tnon_nonzero

 real(dp) :: phnon1(2)
 real(dp), allocatable :: workg(:,:), workg_eq(:,:)
 character(len=500) :: message

! *************************************************************************

 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3);nproc_fft=ngfft(10);me_fft=ngfft(11);nd2=n2/nproc_fft

 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2

 rhog1_eq = zero
 rhor1_eq = zero

 if (itirev == 2) then
   write (message,'(3a,9I4,1a)') 'using time reversal. ',ch10,'Symrel1 = ', symrel1, ch10
 else
   write (message,'(3a,9I4,1a)') 'no time reversal. ',ch10,'Symrel1 = ', symrel1, ch10
 end if
 !call wrtout(std_out, message, 'COLL')
 
 t_tnon_nonzero = (any(abs(tnon) > tol12))

! eventually, for FFT parallelization
! call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 ABI_ALLOCATE(workg,(2,nfft))
 ABI_ALLOCATE(workg_eq,(2,nfft))
 do ispden = 1, nspden

! fft input rhor1 to reciprocal space: uses work* as a buffer
   call fourdp(cplex,workg,rhor1(:,ispden),-1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,0)

! below taken from irrzg and setsym
!  Loop over reciprocal space grid points:
!  loop over local points in workg, and get back transform from rhog1,
!  which is presumed complete on each proc
   ind1=0
   do i3=1,n3
     do i2=1,n2
!       if(fftn2_distrib(i2)/=me_fft)  cycle ! this ind is not to be treated by me_fft
       do i1=1,n1

         ind1=ind1+1
!       r2=ffti2_local(i2+1) - 1
!       ind=n1*(nd2*i3+r2)+i1+1 !this is ind in the current proc

!      Get location of G vector (grid point) centered at 0 0 0
         l1=i1-(i1/id1)*n1-1
         l2=i2-(i2/id2)*n2-1
         l3=i3-(i3/id3)*n3-1

!      Get rotated G vector Gj for each symmetry element
!      -- here we use the TRANSPOSE of symrel1; assuming symrel1 expresses
!      the rotation in real space, the transpose is then appropriate
!      for G space symmetrization (p. 1172d,e of notes, 2 June 1995).
         j1=symrel1(1,1)*l1+symrel1(2,1)*l2+symrel1(3,1)*l3
         j2=symrel1(1,2)*l1+symrel1(2,2)*l2+symrel1(3,2)*l3
         j3=symrel1(1,3)*l1+symrel1(2,3)*l2+symrel1(3,3)*l3

!      Map into [0,n-1] and then add 1 for array index in [1,n]
         k1=1+mod(n1+mod(j1,n1),n1)
         k2=1+mod(n2+mod(j2,n2),n2)
         k3=1+mod(n3+mod(j3,n3),n3)

!      Get linear index of rotated point Gj
         ind2=k1+n1*((k2-1)+n2*(k3-1))
!       r2=ffti2_local(j2+1) - 1
!       ind=n1*(nd2*j3+r2)+j1+1 !this is ind may be in another proc!!

         phnon1(1) = one
         phnon1(2) = zero
         if (t_tnon_nonzero) then
!        compute exp(-2*Pi*I*G dot tau) using original G
! NB: this phase is same as that in irrzg and phnons1, and corresponds to complex conjugate of phase from G to Gj; 
! we use it immediately below, to go _to_ workg(ind1)
! TODO : replace this with complex powers of exp(2pi tnon(1)) etc...
           arg=two_pi*(dble(l1)*tnon(1)+dble(l2)*tnon(2)+dble(l3)*tnon(3))
           phnon1(1) = cos(arg)
           phnon1(2) =-sin(arg)
         end if

!      rho(Strans*G)=exp(2*Pi*I*(G) dot tau_S) rho(G)
         workg_eq (1, ind1) = phnon1(1) * workg(1, ind2) &
&         - phnon1(2) * workg(2, ind2)
         workg_eq (2, ind1) = phnon1(1) * workg(2, ind2) &
&         + phnon1(2) * workg(1, ind2)

       end do
     end do
   end do

! accumulate rhog1_eq
   if (ispden == 1) rhog1_eq = workg_eq

! FFT back to real space to get rhor1_eq
!    Pull out full or spin up density, now symmetrized
   call fourdp(cplex,workg_eq,rhor1_eq(:,ispden),1,mpi_enreg,nfft,ngfft,mpi_enreg%paral_kgb,0)

 end do !nspden

 ABI_DEALLOCATE(workg)
 ABI_DEALLOCATE(workg_eq)

end subroutine rotate_rho
!!***
