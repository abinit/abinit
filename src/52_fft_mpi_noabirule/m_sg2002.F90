!!****m* ABINIT/m_sg2002
!! NAME
!!  m_sg2002
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2002-2007 Stefan Goedecker, CEA Grenoble
!!  Copyright (C) 2014-2020 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

module m_sg2002

 use defs_basis
 use defs_fftdata
 use m_abicore
 use m_errors
 use m_xmpi

 use m_time,         only : timab
 use m_fstrings,     only : itoa
 use m_fftcore,      only : sphere_fft1, fill, scramble, switchreal, switch, mpiswitch,&
&                           unfill, unscramble, unswitchreal, unswitch, unmpiswitch,&
&                           fill_cent, switch_cent, switchreal_cent, mpiswitch_cent, multpot, addrho,&
&                           unfill_cent, unswitchreal_cent, unswitch_cent, unmpiswitch_cent, unscramble,&
&                           mpifft_fg2dbox, mpifft_dbox2fr, mpifft_fr2dbox, mpifft_dbox2fg

 implicit none

 private

 ! Public API:
 !public :: sg2002_seqfourdp   ! seq-FFT of densities and potentials.
 public :: sg2002_mpifourdp    ! MPI-FFT of densities and potentials.
 !public :: mpi_fourwf
 !public :: sg2002_seqfourwf   ! seq-FFT of wavefunctions.
 !public :: sg2002_mpifourwf   ! MPI-FFT of wavefunctions.

! Low-level tools.
! These procedure shouls be accessed via a wrapper that selected the library via fftalg
 public :: sg2002_back           ! G --> R for densities and potentials
 public :: sg2002_forw           ! R --> G for densities and potentials
 public :: sg2002_mpiback_wf     ! G --> R for wavefunctions
 public :: sg2002_mpiforw_wf     ! R --> G for wavefunctions
 public :: sg2002_applypot       ! Compute <G|vloc|u> where u is given in reciprocal space.
 public :: sg2002_applypot_many  ! Compute <G|vloc|u> where u is given in reciprocal space.
 public :: sg2002_accrho         ! Compute rho = weigth_r*Re(u(r))**2 + weigth_i*Im(u(r))**2

contains
!!***

!!****f* m_sg2002/sg2002_back
!! NAME
!!  sg2002_back
!!
!! FUNCTION
!!   CALCULATES THE DISCRETE FOURIER TRANSFORM  in parallel using MPI/OpenMP
!!
!!   ZR(I1,I2,I3)= \sum_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) ZF(j1,j3,j2)
!!
!! Adopt standard convention that isign=1 for backward transform
!!
!! INPUTS:
!!    cplex=1 for real --> complex, 2 for complex --> complex
!!    ZF: input array in G-space (note the switch of i2 and i3)
!!
!!         real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!         imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!
!!         i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!! OUTPUTS:
!!    ZR: output array in R space.
!!
!!         ZR(1,i1,i2,i3,idat)=real(R(i1,i2,i3,idat))
!!         ZR(2,i1,i2,i3,idat)=imag(R(i1,i2,i3,idat))
!!
!!         i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!!
!!    nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!    me_fft: [0:nproc_fft-1] number of processor as returned by MPI_COMM_RANK
!!    n1,n2,n3: logical dimension of the transform. As transform lengths
!!              most products of the prime factors 2,3,5 are allowed.
!!              The detailed table with allowed transform lengths can
!!              be found in subroutine CTRIG
!!    nd1,nd2,nd3: Dimension of ZF and ZR
!!    nd2proc=((nd2-1)/nproc_fft)+1 maximal number of 2nd dim slices
!!    nd3proc=((nd3-1)/nproc_fft)+1 maximal number of 3rd dim slices
!!
!! NOTES:
!!   The maximum number of processors that can reasonably be used is max(n2,n3)
!!   It is very important to find the optimal
!!   value of NCACHE. NCACHE determines the size of the work array ZW, that
!!   has to fit into cache. It has therefore to be chosen to equal roughly
!!    half the size of the physical cache in units of real*8 numbers.
!!   The optimal value of ncache can easily be determined by numerical
!!   experimentation. A too large value of ncache leads to a dramatic
!!   and sudden decrease of performance, a too small value to a to a
!!   slow and less dramatic decrease of performance. If NCACHE is set
!!   to a value so small, that not even a single one dimensional transform
!!   can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fft,m_sg2002
!!
!! CHILDREN
!!
!! SOURCE

subroutine sg2002_back(cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd1eff,nd2proc,nd3proc,option,zf,zr,comm_fft)

 implicit none

!Arguments ------------------------------------
! real space input
 integer,intent(in) :: cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd1eff,nd2proc,nd3proc,option,comm_fft
 real(dp),intent(in) :: zf(2,nd1,nd3,nd2proc,ndat)
 real(dp),intent(out) :: zr(2,nd1eff,nd2,nd3proc,ndat)

!Local variables-------------------------------
!scalars
 integer :: i,j,i1,ic1,ic2,ic3,idat,ierr,includelast,inzee,j2,j2st,j3,jeff,jp2st,lot,lzt
 integer :: ma,mb,n1dfft,n1eff,n2eff,n1zt,ncache,nnd3,nproc_fft,me_fft
 character(len=500) :: msg
!arrays
 real(dp), allocatable :: zt(:,:,:)  ! work arrays for transpositions
 real(dp), allocatable :: zmpi1(:,:,:,:),zmpi2(:,:,:,:) ! work arrays for MPI
 real(dp), allocatable :: zw(:,:,:) ! cache work array
 real(dp) :: tsec(2)
! FFT work arrays
 real(dp), allocatable, dimension(:,:) :: trig1,trig2,trig3
 integer, allocatable, dimension(:) :: after1,now1,before1,after2,now2,before2,after3,now3,before3

! *************************************************************************

 nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)

 ! find cache size that gives optimal performance on machine
 ncache=4*max(n1,n2,n3,1024)

 if (ncache/(4*max(n1,n2,n3))<1) then
   write(msg,'(5a)') &
&    'ncache has to be enlarged to be able to hold at',ch10, &
&    'least one 1-d FFT of each size even though this will',ch10,&
&    'reduce the performance for shorter transform lengths'
   ABI_ERROR(msg)
 end if

! check input
 if (nd1<n1 .or. nd2<n2 .or. nd3<n3) then
   ABI_ERROR("nd1<n1 .or. nd2<n2 .or. nd3<n3")
 end if

 ! Effective n1 and n2 (complex-to-complex or real-to-complex)
 n1eff=n1; n2eff=n2; n1zt=n1
 if (cplex==1) then
   n1eff=(n1+1)/2 ; n2eff=n2/2+1 ; n1zt=2*(n1/2+1)
 end if

 lzt=n2eff
 if (mod(n2eff,2) == 0) lzt=lzt+1
 if (mod(n2eff,4) == 0) lzt=lzt+1

! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ABI_ALLOCATE(trig1,(2,n1))
 ABI_ALLOCATE(after1,(mdata))
 ABI_ALLOCATE(now1,(mdata))
 ABI_ALLOCATE(before1,(mdata))
 ABI_ALLOCATE(trig2,(2,n2))
 ABI_ALLOCATE(after2,(mdata))
 ABI_ALLOCATE(now2,(mdata))
 ABI_ALLOCATE(before2,(mdata))
 ABI_ALLOCATE(trig3,(2,n3))
 ABI_ALLOCATE(after3,(mdata))
 ABI_ALLOCATE(now3,(mdata))
 ABI_ALLOCATE(before3,(mdata))
 ABI_ALLOCATE(zw,(2,ncache/4,2))
 ABI_ALLOCATE(zt,(2,lzt,n1zt))
 ABI_ALLOCATE(zmpi2,(2,n1,nd2proc,nnd3))
 if (nproc_fft>1)  then
   ABI_ALLOCATE(zmpi1,(2,n1,nd2proc,nnd3))
 end if

 call ctrig(n3,trig3,after3,before3,now3,1,ic3)
 call ctrig(n1,trig1,after1,before1,now1,1,ic1)
 call ctrig(n2,trig2,after2,before2,now2,1,ic2)

!DEBUG
! write(std_out,'(a,3i4)' )'sg2002_back,zf n1,n2,n3',n1,n2,n3
! write(std_out,'(a,3i4)' )'nd1,nd2,nd3proc',nd1,nd2,nd3proc
! write(std_out,'(a,3i4)' )'m1,m2,m3',m1,m2,m3
! write(std_out,'(a,3i4)' )'max1,max2,max3',max1,max2,max3
! write(std_out,'(a,3i4)' )'md1,md2proc,md3',md1,md2proc,md3
! write(std_out,'(a,3i4)' )'n1eff,m2eff,m1zt',n1eff,m2eff,m1zt
!ENDDEBUG

 do idat=1,ndat
   ! transform along z axis
   ! input: I1,I3,J2,(Jp2)
   lot=ncache/(4*n3)

   do j2=1,nd2proc
     if (me_fft*nd2proc+j2 <= n2eff) then

       do i1=1,n1,lot
         ma=i1
         mb=min(i1+(lot-1),n1)
         n1dfft=mb-ma+1

         ! input: G1,G3,G2,(Gp2)
         call fill(nd1,nd3,lot,n1dfft,n3,zf(1,i1,1,j2,idat),zw(1,1,1))

         inzee=1
         do i=1,ic3
           call fftstp(lot,n1dfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
&                      trig3,after3(i),now3(i),before3(i),1)
           inzee=3-inzee
         end do

         ! input:  G1,R3,G2,(Gp2)
         ! output: G1,G2,R3,(Gp2)
         call scramble(i1,j2,lot,n1dfft,n1,n3,nd2proc,nd3,zw(1,1,inzee),zmpi2)
       end do
     end if
   end do

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Rp3,(Gp2)
   ! output: G1,G2,G3,Gp2,(Rp3)
   if (nproc_fft>1) then
     call timab(543,1,tsec)
     call xmpi_alltoall(zmpi2,2*n1*nd2proc*nd3proc, &
&                       zmpi1,2*n1*nd2proc*nd3proc,comm_fft,ierr)
     call timab(543,2,tsec)
   end if

   do j3=1,nd3proc
     if (me_fft*nd3proc+j3 <= n3) then
       Jp2st=1
       J2st=1

       ! transform along x axis
       lot=ncache/(4*n1)

       do j=1,n2eff,lot
         ma=j
         mb=min(j+(lot-1),n2eff)
         n1dfft=mb-ma+1

         ! input:  G1,G2,R3,Gp2,(Rp3)
         ! output: G2,G1,R3,Jp2,(Rp3)
         if (nproc_fft == 1) then
           call mpiswitch(j3,n1dfft,Jp2st,J2st,lot,n1,nd2proc,nd3proc,nproc_fft,option,zmpi2,zw(1,1,1))
         else
           call mpiswitch(j3,n1dfft,Jp2st,J2st,lot,n1,nd2proc,nd3proc,nproc_fft,option,zmpi1,zw(1,1,1))
         end if

         ! input:  G2,G1,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         inzee=1
         do i=1,ic1-1
           call fftstp(lot,n1dfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
&                       trig1,after1(i),now1(i),before1(i),1)
           inzee=3-inzee
         end do

         i=ic1
         call fftstp(lot,n1dfft,n1,lzt,n1zt,zw(1,1,inzee),zt(1,j,1), &
&                    trig1,after1(i),now1(i),before1(i),1)
       end do

       ! transform along y axis
       lot=ncache/(4*n2)

       do j=1,n1eff,lot
         ma=j
         mb=min(j+(lot-1),n1eff)
         n1dfft=mb-ma+1
         includelast=1

         if (cplex==1) then
          jeff=2*j-1
          includelast=1
          if (mb==n1eff .and. n1eff*2/=n1) includelast=0
         end if

         ! input:  G2,R1,R3,(Rp3)
         ! output: R1,G2,R3,(Rp3)
         if (cplex==2) then
           call switch(n1dfft,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
         else
           call switchreal(includelast,n1dfft,n2,n2eff,lot,n1zt,lzt,zt(1,1,jeff),zw(1,1,1))
         end if

         inzee=1
         do i=1,ic2-1
           call fftstp(lot,n1dfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
&                       trig2,after2(i),now2(i),before2(i),1)
           inzee=3-inzee
         end do

         i=ic2
         call fftstp(lot,n1dfft,n2,nd1eff,nd2,zw(1,1,inzee),zr(1,j,1,j3,idat), &
&                    trig2,after2(i),now2(i),before2(i),1)

       end do
       ! output: R1,R2,R3,(Rp3)

     end if
   end do
 end do ! idat

 ABI_DEALLOCATE(trig1)
 ABI_DEALLOCATE(after1)
 ABI_DEALLOCATE(now1)
 ABI_DEALLOCATE(before1)
 ABI_DEALLOCATE(trig2)
 ABI_DEALLOCATE(after2)
 ABI_DEALLOCATE(now2)
 ABI_DEALLOCATE(before2)
 ABI_DEALLOCATE(trig3)
 ABI_DEALLOCATE(after3)
 ABI_DEALLOCATE(now3)
 ABI_DEALLOCATE(before3)
 ABI_DEALLOCATE(zmpi2)
 ABI_DEALLOCATE(zw)
 ABI_DEALLOCATE(zt)
 if (nproc_fft>1)  then
   ABI_DEALLOCATE(zmpi1)
 end if

end subroutine sg2002_back
!!***

!----------------------------------------------------------------------

!!****f* m_sg2002/sg2002_forw
!! NAME
!!  sg2002_forw
!!
!! FUNCTION
!!   Adopt standard convention that isign=-1 for forward transform
!!   CALCULATES THE DISCRETE FOURIERTRANSFORM ZF(I1,I3,I2)=
!!   S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) ZR(j1,j2,j3)
!!   in parallel using MPI/OpenMP and BLAS library calls.
!!
!! INPUTS
!!    ZR: input array
!!         ZR(1,i1,i2,i3,idat)=real(R(i1,i2,i3,idat))
!!         ZR(2,i1,i2,i3,idat)=imag(R(i1,i2,i3,idat))
!!         i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!! OUTPUTS
!!    ZF: output array (note the switch of i2 and i3)
!!         real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!         imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!         i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!!    nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!    me_fft: [0:nproc_fft-1] number of processor as returned by MPI_COMM_RANK
!!     n1,n2,n3: logical dimension of the transform. As transform lengths
!!               most products of the prime factors 2,3,5 are allowed.
!!              The detailed table with allowed transform lengths can
!!              be found in subroutine CTRIG
!!     nd1,nd2,nd3: Dimension of ZR and ZF
!!    nd2proc=((nd2-1)/nproc_fft)+1 maximal number of 2nd dim slices
!!    nd3proc=((nd3-1)/nproc_fft)+1 maximal number of 3rd dim slices
!!
!! NOTES
!!  SHOULD describe nd1eff
!!  SHOULD put cplex and nd1eff in OMP declarations
!!  SHOULD describe the change of value of nd2prod
!!
!!  The maximum number of processors that can reasonably be used is max(n2,n3)
!!
!!  It is very important to find the optimal
!!  value of NCACHE. NCACHE determines the size of the work array ZW, that
!!  has to fit into cache. It has therefore to be chosen to equal roughly
!!   half the size of the physical cache in units of real*8 numbers.
!!  The optimal value of ncache can easily be determined by numerical
!!  experimentation. A too large value of ncache leads to a dramatic
!!  and sudden decrease of performance, a too small value to a to a
!!  slow and less dramatic decrease of performance. If NCACHE is set
!!  to a value so small, that not even a single one dimensional transform
!!  can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fft,m_sg2002
!!
!! CHILDREN
!!
!! SOURCE

subroutine sg2002_forw(cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd1eff,nd2proc,nd3proc,option,zr,zf,comm_fft)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,comm_fft
 integer,intent(in) :: ndat,n1,n2,n3,nd1,nd2,nd3,nd1eff,nd2proc,nd3proc,option
!arrays
 real(dp),intent(in) :: zr(2,nd1eff,nd2,nd3proc,ndat)
 real(dp),intent(out) :: zf(2,nd1,nd3,nd2proc,ndat)

!Local variables-------------------------------
!scalars
 integer :: i,j,i1,ic1,ic2,ic3,idat,ierr,inzee,j2,j2st,j3,jp2st,lot,lzt
 integer :: ma,mb,n1dfft,n1eff,n2eff,n1zt,ncache,nnd3,nproc_fft,me_fft
 character(len=500) :: msg
!arrays
 real(dp), allocatable :: zt(:,:,:) ! work arrays for transpositions
 real(dp), allocatable :: zmpi1(:,:,:,:),zmpi2(:,:,:,:) ! work arrays for MPI
 real(dp), allocatable :: zw(:,:,:) ! cache work array
 real(dp) :: tsec(2)
! FFT work arrays
 real(dp), allocatable, dimension(:,:) :: trig1,trig2,trig3
 integer, allocatable, dimension(:) :: after1,now1,before1,after2,now2,before2,after3,now3,before3

! *************************************************************************

 nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)

 ! find cache size that gives optimal performance on machine
 ncache=4*max(n1,n2,n3,1024)
 if (ncache/(4*max(n1,n2,n3))<1) then
   write(msg,'(5a)')&
&     'ncache has to be enlarged to be able to hold at',ch10, &
&     'least one 1-d FFT of each size even though this will',ch10,&
&     'reduce the performance for shorter transform lengths'
   ABI_ERROR(msg)
 end if

 ! check input
 if (nd1<n1 .or. nd2<n2 .or. nd3<n3) then
   ABI_ERROR("nd1<n1 .or. nd2<n2 .or. nd3<n3")
 end if

!Effective n1 and n2 (complex-to-complex or real-to-complex)
 n1eff=n1; n2eff=n2; n1zt=n1
 if (cplex==1) then
   n1eff=(n1+1)/2; n2eff=n2/2+1; n1zt=2*(n1/2+1)
 end if

 lzt=n2eff
 if (mod(n2eff,2) == 0) lzt=lzt+1
 if (mod(n2eff,4) == 0) lzt=lzt+1

 ! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ABI_ALLOCATE(trig1,(2,n1))
 ABI_ALLOCATE(after1,(mdata))
 ABI_ALLOCATE(now1,(mdata))
 ABI_ALLOCATE(before1,(mdata))
 ABI_ALLOCATE(trig2,(2,n2))
 ABI_ALLOCATE(after2,(mdata))
 ABI_ALLOCATE(now2,(mdata))
 ABI_ALLOCATE(before2,(mdata))
 ABI_ALLOCATE(trig3,(2,n3))
 ABI_ALLOCATE(after3,(mdata))
 ABI_ALLOCATE(now3,(mdata))
 ABI_ALLOCATE(before3,(mdata))
 ABI_ALLOCATE(zw,(2,ncache/4,2))
 ABI_ALLOCATE(zt,(2,lzt,n1zt))
 ABI_ALLOCATE(zmpi2,(2,n1,nd2proc,nnd3))
 if (nproc_fft>1)  then
   ABI_ALLOCATE(zmpi1,(2,n1,nd2proc,nnd3))
 end if

 call ctrig(n2,trig2,after2,before2,now2,-1,ic2)
 call ctrig(n1,trig1,after1,before1,now1,-1,ic1)
 call ctrig(n3,trig3,after3,before3,now3,-1,ic3)

 do idat=1,ndat
   do j3=1,nd3proc
     if (me_fft*(nd3proc)+j3 <= n3) then
       Jp2st=1; J2st=1

       ! transform along y axis
       ! input: R1,R2,R3,(Rp3)
       lot=ncache/(4*n2)

       do j=1,n1eff,lot
         ma=j
         mb=min(j+(lot-1),n1eff)
         n1dfft=mb-ma+1
         i=1
         call fftstp(nd1eff,n1dfft,nd2,lot,n2,zr(1,j,1,j3,idat),zw(1,1,1), &
&                    trig2,after2(i),now2(i),before2(i),-1)

         inzee=1
         do i=2,ic2
           call fftstp(lot,n1dfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
&                       trig2,after2(i),now2(i),before2(i),-1)
            inzee=3-inzee
         end do

         !  input: R1,G2,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         if(cplex==2)then
           call unswitch(n1dfft,n2,lot,n1zt,lzt,zw(1,1,inzee),zt(1,1,j))
         else
           call unswitchreal(n1dfft,n2,n2eff,lot,n1zt,lzt,zw(1,1,inzee),zt(1,1,2*j-1))
         end if
       end do

       ! transform along x axis
       ! input: G2,R1,R3,(Rp3)
       lot=ncache/(4*n1)

       do j=1,n2eff,lot
         ma=j
         mb=min(j+(lot-1),n2eff)
         n1dfft=mb-ma+1

         i=1
         call fftstp(lzt,n1dfft,n1zt,lot,n1,zt(1,j,1),zw(1,1,1), &
&                    trig1,after1(i),now1(i),before1(i),-1)

         inzee=1
         do i=2,ic1
           call fftstp(lot,n1dfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
&                     trig1,after1(i),now1(i),before1(i),-1)
           inzee=3-inzee
         end do
         ! output: G2,G1,R3,(Rp3)

         ! input:  G2,G1,R3,Gp2,(Rp3)
         ! output: G1,G2,R3,Gp2,(Rp3)
         ! write(std_out,*) 'J2st,Jp2st',J2st,Jp2st
         if (nproc_fft == 1) then
           call unmpiswitch(j3,n1dfft,Jp2st,J2st,lot,n1,nd2proc,nd3proc,nproc_fft,option,zw(1,1,inzee),zmpi2)
         else
           call unmpiswitch(j3,n1dfft,Jp2st,J2st,lot,n1,nd2proc,nd3proc,nproc_fft,option,zw(1,1,inzee),zmpi1)
         end if
       end do

     end if
   end do ! j3

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Gp2,(Rp3)
   ! output: G1,G2,R3,Rp3,(Gp2)
   if (nproc_fft>1) then
     call timab(544,1,tsec)
     call xmpi_alltoall(zmpi1,2*n1*nd2proc*nd3proc, &
&                       zmpi2,2*n1*nd2proc*nd3proc,comm_fft,ierr)
     call timab(544,2,tsec)
   end if

   ! transform along z axis
   ! input: G1,G2,R3,(Gp2)
   lot=ncache/(4*n3)

   do j2=1,nd2proc
     if (me_fft*(nd2proc)+j2 <= n2eff) then
       do i1=1,n1,lot
         ma=i1
         mb=min(i1+(lot-1),n1)
         n1dfft=mb-ma+1

         ! input:  G1,G2,R3,(Gp2)
         ! output: G1,R3,G2,(Gp2)
         call unscramble(i1,j2,lot,n1dfft,n1,n3,nd2proc,nd3,zmpi2,zw(1,1,1))

         inzee=1
         do i=1,ic3
           call fftstp(lot,n1dfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
&            trig3,after3(i),now3(i),before3(i),-1)
           inzee=3-inzee
         end do

         call unfill(nd1,nd3,lot,n1dfft,n3,zw(1,1,inzee),zf(1,i1,1,j2,idat))
         ! output: G1,G3,G2,(Gp2)
       end do
     end if
   end do

 end do ! idat

 ABI_DEALLOCATE(trig1)
 ABI_DEALLOCATE(after1)
 ABI_DEALLOCATE(now1)
 ABI_DEALLOCATE(before1)
 ABI_DEALLOCATE(trig2)
 ABI_DEALLOCATE(after2)
 ABI_DEALLOCATE(now2)
 ABI_DEALLOCATE(before2)
 ABI_DEALLOCATE(trig3)
 ABI_DEALLOCATE(after3)
 ABI_DEALLOCATE(now3)
 ABI_DEALLOCATE(before3)
 ABI_DEALLOCATE(zmpi2)
 ABI_DEALLOCATE(zw)
 ABI_DEALLOCATE(zt)
 if (nproc_fft>1)  then
   ABI_DEALLOCATE(zmpi1)
 end if

end subroutine sg2002_forw
!!***

!----------------------------------------------------------------------

!!****f* m_sg2002/sg2002_mpiback_wf
!! NAME
!!  sg2002_mpiback_wf
!!
!! FUNCTION
!!   Does multiple 3-dim backward FFTs from Fourier into real space
!!   Adopt standard convention that isign=1 for backward transform
!!
!!   CALCULATES THE DISCRETE FOURIER TRANSFORM ZF(I1,I2,I3)=
!!
!!   S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) ZF(j1,j3,j2)
!!
!!   in parallel using MPI/OpenMP.
!!
!! INPUTS:
!!    icplexwf=1 if wavefunction is real, 2 if complex
!!    ndat=Number of wavefunctions to transform.
!!    n1,n2,n3: logical dimension of the transform. As transform lengths
!!              most products of the prime factors 2,3,5 are allowed.
!!              The detailed table with allowed transform lengths can be found in subroutine CTRIG
!!    nd1,nd2,nd3: Leading Dimension of ZR
!!    nd3proc=((nd3-1)/nproc_fft)+1 maximal number of big box 3rd dim slices for one proc
!!    max1 is positive or zero; m1 >=max1+1
!!      i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!!      then, if m1 > max1+1, one has min1=max1-m1+1 and
!!      i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!!    max2 and max3 have a similar definition of range
!!    m1,m2,m3=Size of the box enclosing the G-sphere.
!!    md1,md2,md3: Dimension of ZF given on the **small** FFT box.
!!    md2proc=((md2-1)/nproc_fft)+1 maximal number of small box 2nd dim slices for one proc
!!    nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!    comm_fft=MPI communicator for the FFT.
!!    ZF: input array (note the switch of i2 and i3)
!!          real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!          imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!
!! OUTPUTS
!!    ZR: output array
!!          ZR(1,i1,i2,i3,idat)=real(R(i1,i2,i3,idat))
!!          ZR(2,i1,i2,i3,idat)=imag(R(i1,i2,i3,idat))
!!        i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!!
!! NOTES
!!   The maximum number of processors that can reasonably be used is max(n2/2,n3/2)
!!
!!   It is very important to find the optimal
!!   value of NCACHE. NCACHE determines the size of the work array ZW, that
!!   has to fit into cache. It has therefore to be chosen to equal roughly
!!   half the size of the physical cache in units of real*8 numbers.
!!   The optimal value of ncache can easily be determined by numerical
!!   experimentation. A too large value of ncache leads to a dramatic
!!   and sudden decrease of performance, a too small value to a to a
!!   slow and less dramatic decrease of performance. If NCACHE is set
!!   to a value so small, that not even a single one dimensional transform
!!   can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!
!! SOURCE

subroutine sg2002_mpiback_wf(icplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc,&
&  max1,max2,max3,m1,m2,m3,md1,md2proc,md3,zf,zr,comm_fft)

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: icplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc
 integer,intent(in) :: max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft
 real(dp),intent(in) :: zf(2,md1,md3,md2proc,ndat)
 real(dp),intent(out) :: zr(2,nd1,nd2,nd3proc,ndat)

!Local variables-------------------------------
 integer :: i,j,i1,i2,ic1,ic2,ic3,idat,ierr,inzee,includelast
 integer :: ioption,j2,j3,j2st,jp2st,jeff,lot,lzt,m1zt,ma,mb,n1dfft,nnd3
 integer :: m2eff,ncache,n1eff,n1half,nproc_fft,me_fft
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: zt(:,:,:) ! work arrays for transpositions
 real(dp),allocatable :: zmpi1(:,:,:,:,:),zmpi2(:,:,:,:,:)  ! work arrays for MPI
 real(dp),allocatable :: zw(:,:,:) ! cache work array
! FFT work arrays
 real(dp),allocatable :: trig1(:,:),trig2(:,:),trig3(:,:)
 integer,allocatable :: after1(:),now1(:),before1(:),after2(:)
 integer,allocatable :: now2(:),before2(:),after3(:),now3(:),before3(:)
 real(dp) :: tsec(2)

! *************************************************************************

 ! call timab(541,1,tsec)
 ! FIXME must provide a default value but which one?
 ! ioption = 0
 ioption = 1
 !if (paral_kgb==1) ioption=1

 nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)

 ! Find cache size that gives optimal performance on machine
 ncache=4*max(n1,n2,n3,1024)
 if (ncache/(4*max(n1,n2,n3))<1) then
   write(msg,"(5a)") &
&    'ncache has to be enlarged to be able to hold at',ch10, &
&    'least one 1-d FFT of each size even though this will',ch10,&
&    'reduce the performance for shorter transform lengths'
    ABI_ERROR(msg)
 end if

 ! Effective m1 and m2 (complex-to-complex or real-to-complex)
 n1eff=n1; m2eff=m2; m1zt=n1
 if (icplexwf==1) then
   n1eff=(n1+1)/2; m2eff=m2/2+1; m1zt=2*(n1/2+1)
 end if

 lzt=m2eff
 if (mod(m2eff,2)==0) lzt=lzt+1
 if (mod(m2eff,4)==0) lzt=lzt+1

 ! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ABI_ALLOCATE(trig1,(2,n1))
 ABI_ALLOCATE(after1,(mdata))
 ABI_ALLOCATE(now1,(mdata))
 ABI_ALLOCATE(before1,(mdata))
 ABI_ALLOCATE(trig2,(2,n2))
 ABI_ALLOCATE(after2,(mdata))
 ABI_ALLOCATE(now2,(mdata))
 ABI_ALLOCATE(before2,(mdata))
 ABI_ALLOCATE(trig3,(2,n3))
 ABI_ALLOCATE(after3,(mdata))
 ABI_ALLOCATE(now3,(mdata))
 ABI_ALLOCATE(before3,(mdata))

 ! Allocate cache work array and work arrays for MPI transpositions.
 ABI_ALLOCATE(zw,(2,ncache/4,2))
 ABI_ALLOCATE(zt,(2,lzt,m1zt))
 ABI_ALLOCATE(zmpi2,(2,md1,md2proc,nnd3,ndat))
 if (nproc_fft>1)  then
   ABI_ALLOCATE(zmpi1,(2,md1,md2proc,nnd3,ndat))
 end if

 ! Compute twiddle coefficients.
 call ctrig(n3,trig3,after3,before3,now3,1,ic3)
 call ctrig(n1,trig1,after1,before1,now1,1,ic1)
 call ctrig(n2,trig2,after2,before2,now2,1,ic2)

!DEBUG
! write(std_out,'(2a,3i4)' )itoa(me_fft),': sg2002_mpiback_wf,zf n1,n2,n3',n1,n2,n3
! write(std_out,'(2a,3i4)' )itoa(me_fft),': nd1,nd2,nd3proc',nd1,nd2,nd3proc
! write(std_out,'(2a,3i4)' )itoa(me_fft),': m1,m2,m3',m1,m2,m3
! write(std_out,'(2a,3i4)' )itoa(me_fft),': max1,max2,max3',max1,max2,max3
! write(std_out,'(2a,3i4)' )itoa(me_fft),': md1,md2proc,md3',md1,md2proc,md3
! write(std_out,'(2a,3i4)' )itoa(me_fft),'n1eff,m2eff,m1zt',n1eff,m2eff,m1zt
!ENDDEBUG

 do idat=1,ndat

    ! transform along z axis
    ! input: G1,G3,G2,(Gp2)
    lot=ncache/(4*n3)

    zw(:,:,:)=zero
    zt(:,:,:)=zero

    ! Loop over the y planes treated by this node and trasform n1ddft G_z lines.
    do j2=1,md2proc

      ! if (me_fft*md2proc+j2<=m2eff) then !a faire plus tard

      do i1=1,m1,lot
        ma=i1
        mb=min(i1+(lot-1),m1)
        n1dfft=mb-ma+1

        ! zero-pad n1dfft G_z lines
        ! input:  G1,G3,G2,(Gp2)
        ! output: G1,R3,G2,(Gp2)
        call fill_cent(md1,md3,lot,n1dfft,max3,m3,n3,zf(1,i1,1,j2,idat),zw(1,1,1))

        ! Transform along z.
        inzee=1
        do i=1,ic3
          call fftstp(lot,n1dfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
&                     trig3,after3(i),now3(i),before3(i),1)
          inzee=3-inzee
        end do

        ! Local rotation.
        ! input:  G1,R3,G2,(Gp2)
        ! output: G1,G2,R3,(Gp2)
        call scramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zw(1,1,inzee),zmpi2(:,:,:,:,idat))
      end do
      !
    end do ! j2

    ! Interprocessor data transposition
    ! input:  G1,G2,R3,Rp3,(Gp2)
    ! output: G1,G2,R3,Gp2,(Rp3)
    if (nproc_fft>1) then
      call timab(543,1,tsec)
      call xmpi_alltoall(zmpi2(:,:,:,:,idat),2*md1*md2proc*nd3proc, &
&                        zmpi1(:,:,:,:,idat),2*md1*md2proc*nd3proc,comm_fft,ierr)
      call timab(543,2,tsec)
    end if

    ! Loop over the z treated by this node.
    do j3=1,nd3proc
      !j3glob = j3 + me_fft*nd3proc
      if (me_fft*nd3proc+j3 <= n3) then
        Jp2st=1; J2st=1

        lot=ncache/(4*n1)

        ! Loop over G_y in the small box.
        do j=1,m2eff,lot
          ma=j
          mb=min(j+(lot-1),m2eff)
          n1dfft=mb-ma+1

          ! Zero-pad input.
          ! input:  G1,G2,R3,JG2,(Rp3)
          ! output: G2,G1,R3,JG2,(Rp3)
          if (nproc_fft==1) then
            call mpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot,max1,md1,m1,n1,&
&             md2proc,nd3proc,nproc_fft,ioption,zmpi2(:,:,:,:,idat),zw(1,1,1),max2,m2,n2)
          else
            call mpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot,max1,md1,m1,n1,&
&             md2proc,nd3proc,nproc_fft,ioption,zmpi1(:,:,:,:,idat),zw(1,1,1),max2,m2,n2)
          end if

          ! Transform along x
          ! input:  G2,G1,R3,(Rp3)
          ! output: G2,R1,R3,(Rp3)
          inzee=1
          do i=1,ic1-1
            call fftstp(lot,n1dfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
&                       trig1,after1(i),now1(i),before1(i),1)
            inzee=3-inzee
          end do

          i=ic1
          call fftstp(lot,n1dfft,n1,lzt,m1zt,zw(1,1,inzee),zt(1,j,1), &
&                     trig1,after1(i),now1(i),before1(i),1)
        end do

        ! Transform along y axis (take into account c2c or c2r case).
        ! Must loop over the full box.
        lot=ncache/(4*n2)

        do j=1,n1eff,lot
          ma=j
          mb=min(j+(lot-1),n1eff)
          n1dfft=mb-ma+1
          includelast=1
          if (icplexwf==1) then
            jeff=2*j-1
            if (mb==n1eff .and. n1eff*2/=n1) includelast=0
          end if

          ! Zero-pad the input.
          ! input:  G2,R1,R3,(Rp3)
          ! output: R1,G2,R3,(Rp3)
          if (icplexwf==2) then
            call switch_cent(n1dfft,max2,m2,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
          else
            call switchreal_cent(includelast,n1dfft,max2,n2,lot,m1zt,lzt,zt(1,1,jeff),zw(1,1,1))
          end if

          ! input:  R1,G2,R3,(Rp3)
          ! output: R1,R2,R3,(Rp3)
          inzee=1
          do i=1,ic2-1
            call fftstp(lot,n1dfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
&                        trig2,after2(i),now2(i),before2(i),1)
            inzee=3-inzee
          end do

          i=ic2

        call fftstp(lot,n1dfft,n2,nd1,nd2,zw(1,1,inzee),zr(1,j,1,j3,idat), &
&                     trig2,after2(i),now2(i),before2(i),1)


        end do

        ! Treat real wavefunctions.
        if (icplexwf==1) then
          n1half=n1/2
          ! If odd
          if (n1half*2/=n1) then
            do i2=1,n2
              zr(1,n1,i2,j3,idat)=zr(1,n1eff,i2,j3,idat)
              zr(2,n1,i2,j3,idat)=zero
            end do
          end if
          do i2=1,n2
            do i1=n1half,1,-1
              zr(1,2*i1-1,i2,j3,idat)=zr(1,i1,i2,j3,idat)
              zr(1,2*i1  ,i2,j3,idat)=zr(2,i1,i2,j3,idat)
              zr(2,2*i1-1,i2,j3,idat)=zero
              zr(2,2*i1  ,i2,j3,idat)=zero
            end do
          end do
        end if

      end if

   end do ! j3
 end do ! idat

 ABI_DEALLOCATE(trig1)
 ABI_DEALLOCATE(after1)
 ABI_DEALLOCATE(now1)
 ABI_DEALLOCATE(before1)
 ABI_DEALLOCATE(trig2)
 ABI_DEALLOCATE(after2)
 ABI_DEALLOCATE(now2)
 ABI_DEALLOCATE(before2)
 ABI_DEALLOCATE(trig3)
 ABI_DEALLOCATE(after3)
 ABI_DEALLOCATE(now3)
 ABI_DEALLOCATE(before3)
 ABI_DEALLOCATE(zmpi2)
 ABI_DEALLOCATE(zw)
 ABI_DEALLOCATE(zt)
 if (nproc_fft>1)  then
   ABI_DEALLOCATE(zmpi1)
 end if

 !call timab(541,2,tsec)

end subroutine sg2002_mpiback_wf
!!***

!----------------------------------------------------------------------

!!****f* m_sg2002/sg2002_mpiforw_wf
!! NAME
!!  sg2002_mpiforw_wf
!!
!! FUNCTION
!!   Does multiple 3-dim backward FFTs from real into Fourier space
!!   Adopt standard convention that isign=-1 for forward transform
!!   CALCULATES THE DISCRETE FOURIERTRANSFORM
!!
!!   ZF(I1,I3,I2)=S_(j1,j2,j3) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) ZR(j1,j2,j3)
!!
!!   in parallel using MPI/OpenMP.
!!
!! INPUT:
!!   ZR: input array
!!        ZR(1,i1,i2,i3,idat)=real(R(i1,i2,i3,idat))
!!        ZR(2,i1,i2,i3,idat)=imag(R(i1,i2,i3,idat))
!!        i1=1,n1 , i2=1,n2 , i3=1,n3 , idat=1,ndat
!!   NOTE that ZR is changed by the routine
!!
!!   n1,n2,n3: logical dimension of the transform. As transform lengths
!!             most products of the prime factors 2,3,5 are allowed.
!!             The detailed table with allowed transform lengths can
!!             be found in subroutine CTRIG
!!   nd1,nd2,nd3: Dimension of ZR
!!   nd3proc=((nd3-1)/nproc_fft)+1  maximal number of big box 3rd dim slices for one proc
!!
!! OUTPUT:
!!   ZF: output array (note the switch of i2 and i3)
!!        real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!        imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!   max1 is positive or zero ; m1 >=max1+1
!!     i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!!     then, if m1 > max1+1, one has min1=max1-m1+1 and
!!     i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!!     i2 and i3 have a similar definition of range
!!   idat=1,ndat
!!   md1,md2,md3: Dimension of ZF
!!   md2proc=((md2-1)/nproc_fft)+1  maximal number of small box 2nd dim slices for one proc
!!   nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!   me_fft: [0:nproc-1] rank of the processor in the FFT communicator.
!!   comm_fft=MPI communicator for parallel FFT.
!!
!! NOTES
!!  The maximum number of processors that can reasonably be used is max(n2/2,n3/2)
!!
!!  It is very important to find the optimal
!!  value of NCACHE. NCACHE determines the size of the work array ZW, that
!!  has to fit into cache. It has therefore to be chosen to equal roughly
!!   half the size of the physical cache in units of real*8 numbers.
!!  The optimal value of ncache can easily be determined by numerical
!!  experimentation. A too large value of ncache leads to a dramatic
!!  and sudden decrease of performance, a too small value to a to a
!!  slow and less dramatic decrease of performance. If NCACHE is set
!!  to a value so small, that not even a single one dimensional transform
!!  can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!
!! SOURCE

subroutine sg2002_mpiforw_wf(icplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc,&
&        max1,max2,max3,m1,m2,m3,md1,md2proc,md3,zr,zf,comm_fft)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icplexwf,ndat,n1,n2,n3,nd1,nd2,nd3proc
 integer,intent(in) :: max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft
!arrays
 real(dp),intent(inout) :: zr(2,nd1,nd2,nd3proc,ndat)
 real(dp),intent(out) :: zf(2,md1,md3,md2proc,ndat)

!Local variables-------------------------------
!scalars
 integer :: i,j,i1,i2,i3,ic1,ic2,ic3,idat,ierr,inzee,nproc_fft,me_fft
 integer :: ioption,j2,j3,j2st,jp2st,lot,lzt,m1zt,ma,mb,n1dfft,nnd3
 integer :: m2eff,ncache,n1eff,n1half,i1inv,i2inv,i3inv
 character(len=500) :: msg
!arrays
 real(dp), allocatable :: zt(:,:,:) ! work arrays for transpositions
 real(dp), allocatable :: zmpi1(:,:,:,:,:),zmpi2(:,:,:,:,:) ! work arrays for MPI
 real(dp), allocatable :: zw(:,:,:) ! cache work array
! FFT work arrays
 real(dp), allocatable :: trig1(:,:),trig2(:,:),trig3(:,:)
 integer, allocatable :: after1(:),now1(:),before1(:),after2(:),now2(:),before2(:),after3(:),now3(:),before3(:)
 real(dp) :: tsec(2)

! *************************************************************************

 ! call timab(542,1,tsec)

 ! FIXME must provide a default value but which one?
 !ioption = 0
 ioption = 1
 !if (paral_kgb==1) ioption=1

 nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)

 ! find cache size that gives optimal performance on machine
 ncache=4*max(n1,n2,n3,1024)
 if (ncache/(4*max(n1,n2,n3))<1) then
   write(msg,'(5a)') &
&    'ncache has to be enlarged to be able to hold at',ch10, &
&    'least one 1-d FFT of each size even though this will',ch10,&
&    'reduce the performance for shorter transform lengths'
   ABI_ERROR(msg)
 end if

 ! Effective m1 and m2 (complex-to-complex or real-to-complex)
 n1eff=n1; m2eff=m2; m1zt=n1
 if (icplexwf==1) then
   n1eff=(n1+1)/2; m2eff=m2/2+1; m1zt=2*(n1/2+1)
 end if

 lzt=m2eff
 if (mod(m2eff,2)==0) lzt=lzt+1
 if (mod(m2eff,4)==0) lzt=lzt+1

 ! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ABI_ALLOCATE(trig1,(2,n1))
 ABI_ALLOCATE(after1,(mdata))
 ABI_ALLOCATE(now1,(mdata))
 ABI_ALLOCATE(before1,(mdata))
 ABI_ALLOCATE(trig2,(2,n2))
 ABI_ALLOCATE(after2,(mdata))
 ABI_ALLOCATE(now2,(mdata))
 ABI_ALLOCATE(before2,(mdata))
 ABI_ALLOCATE(trig3,(2,n3))
 ABI_ALLOCATE(after3,(mdata))
 ABI_ALLOCATE(now3,(mdata))
 ABI_ALLOCATE(before3,(mdata))
 ABI_ALLOCATE(zw,(2,ncache/4,2))
 ABI_ALLOCATE(zt,(2,lzt,m1zt))
 ABI_ALLOCATE(zmpi2,(2,md1,md2proc,nnd3,ndat))
 if (nproc_fft>1)  then
   ABI_ALLOCATE(zmpi1,(2,md1,md2proc,nnd3,ndat))
 end if

 call ctrig(n2,trig2,after2,before2,now2,-1,ic2)
 call ctrig(n1,trig1,after1,before1,now1,-1,ic1)
 call ctrig(n3,trig3,after3,before3,now3,-1,ic3)

!DEBUG
! write(std_out,'(2a,3i4)' )itoa(me_fft),'sg2002_mpiforw_wf, enter', i1,i2,i3,zr,n1,n2,n3',n1,n2,n3
! write(std_out,'(2a,3i4)' )itoa(me_fft),'nd1,nd2,nd3proc',nd1,nd2,nd3proc
! write(std_out,'(2a,3i4)' )itoa(me_fft),'m1,m2,m3',m1,m2,m3
! write(std_out,'(2a,3i4)' )itoa(me_fft),'max1,max2,max3',max1,max2,max3
! write(std_out,'(2a,3i4)' )itoa(me_fft),'md1,md2proc,md3',md1,md2proc,md3
! write(std_out,'(2a,3i4)' )itoa(me_fft),'n1eff,m2eff,m1zt',n1eff,m2eff,m1zt
!ENDDEBUG

  do idat=1,ndat
    ! Loop over the z-planes treated by this node
    do j3=1,nd3proc

       if (me_fft*nd3proc+j3 <= n3) then
         Jp2st=1
         J2st=1

         ! Treat real wavefunctions.
         if (icplexwf==1) then
           n1half=n1/2
           do i2=1,n2
             do i1=1,n1half
               zr(1,i1,i2,j3,idat)=zr(1,2*i1-1,i2,j3,idat)
               zr(2,i1,i2,j3,idat)=zr(1,2*i1  ,i2,j3,idat)
             end do
           end do
           ! If odd
           if(n1half*2/=n1)then
             do i2=1,n2
               zr(1,n1eff,i2,j3,idat)=zr(1,n1,i2,j3,idat)
               zr(2,n1eff,i2,j3,idat)=zero
             end do
           end if
         end if

         ! transform along y axis
         ! input: R1,R2,R3,(Rp3)
         ! input: R1,G2,R3,(Rp3)
         lot=ncache/(4*n2)

         do j=1,n1eff,lot
           ma=j
           mb=min(j+(lot-1),n1eff)
           n1dfft=mb-ma+1
           i=1
           call fftstp(nd1,n1dfft,nd2,lot,n2,zr(1,j,1,j3,idat),zw(1,1,1), &
&                      trig2,after2(i),now2(i),before2(i),-1)

           inzee=1
           do i=2,ic2
             call fftstp(lot,n1dfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
&                         trig2,after2(i),now2(i),before2(i),-1)
              inzee=3-inzee
           end do

           ! input:  R1,G2,R3,(Rp3)
           ! output: G2,R1,R3,(Rp3)
           if(icplexwf==2)then
             call unswitch_cent(n1dfft,max2,m2,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,j))
           else
             call unswitchreal_cent(n1dfft,max2,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,2*j-1))
           end if
         end do

         ! transform along x axis
         ! input: G2,R1,R3,(Rp3)
         lot=ncache/(4*n1)

         do j=1,m2eff,lot
           ma=j
           mb=min(j+(lot-1),m2eff)
           n1dfft=mb-ma+1
           i=1
           call fftstp(lzt,n1dfft,m1zt,lot,n1,zt(1,j,1),zw(1,1,1), &
&                      trig1,after1(i),now1(i),before1(i),-1)

           inzee=1
           do i=2,ic1
             call fftstp(lot,n1dfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
&                        trig1,after1(i),now1(i),before1(i),-1)
             inzee=3-inzee
           end do
           ! output: G2,G1,R3,(Rp3)

           ! input:  G2,G1,R3,Gp2,(Rp3)
           ! output: G1,G2,R3,Gp2,(Rp3)
           if (nproc_fft==1) then
             call unmpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot,max1,md1,m1,n1,&
&              md2proc,nd3proc,nproc_fft,ioption,zw(1,1,inzee),zmpi2(:,:,:,:,idat))
           else
             call unmpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot,max1,md1,m1,n1,&
&              md2proc,nd3proc,nproc_fft,ioption,zw(1,1,inzee),zmpi1(:,:,:,:,idat))
           end if
         end do

        end if
     end do ! j3

     ! Interprocessor data transposition
     ! input:  G1,G2,R3,Gp2,(Rp3)
     ! output: G1,G2,R3,Rp3,(Gp2)
     if (nproc_fft>1) then
        call timab(544,1,tsec)
        call xmpi_alltoall(zmpi1(:,:,:,:,idat),2*md1*md2proc*nd3proc, &
&                          zmpi2(:,:,:,:,idat),2*md1*md2proc*nd3proc,comm_fft,ierr)

        call timab(544,2,tsec)
     end if

     ! transform along z axis
     ! input: G1,G2,R3,(Gp2)
     lot=ncache/(4*n3)

     do j2=1,md2proc
       if (me_fft*md2proc+j2 <= m2eff) then
         ! write(std_out,*)' forwf_wf : before unscramble, j2,md2proc,me_fft,m2=',j2,md2proc,me_fft,m2
         do i1=1,m1,lot
           ma=i1
           mb=min(i1+(lot-1),m1)
           n1dfft=mb-ma+1

           ! input:  G1,G2,R3,(Gp2)
           ! output: G1,R3,G2,(Gp2)
           call unscramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zmpi2(:,:,:,:,idat),zw(1,1,1))

           inzee=1
           do i=1,ic3
             call fftstp(lot,n1dfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
&              trig3,after3(i),now3(i),before3(i),-1)
             inzee=3-inzee
           end do

           call unfill_cent(md1,md3,lot,n1dfft,max3,m3,n3,zw(1,1,inzee),zf(1,i1,1,j2,idat))
           ! output: G1,G3,G2,(Gp2)
         end do
       end if
     end do

     if (icplexwf==1) then
       ! Complete missing values with complex conjugate
       ! Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.
       do i3=1,m3
         i3inv=m3+2-i3
         if(i3==1)i3inv=1

         if (m2eff>1) then
           do i2=2,m2eff
             i2inv=m2+2-i2
             zf(1,1,i3inv,i2inv,idat)= zf(1,1,i3,i2,idat)
             zf(2,1,i3inv,i2inv,idat)=-zf(2,1,i3,i2,idat)
             do i1=2,m1
               i1inv=m1+2-i1
               zf(1,i1inv,i3inv,i2inv,idat)= zf(1,i1,i3,i2,idat)
               zf(2,i1inv,i3inv,i2inv,idat)=-zf(2,i1,i3,i2,idat)
             end do
           end do
         end if
       end do
     end if

 end do ! idat

 ABI_DEALLOCATE(trig1)
 ABI_DEALLOCATE(after1)
 ABI_DEALLOCATE(now1)
 ABI_DEALLOCATE(before1)
 ABI_DEALLOCATE(trig2)
 ABI_DEALLOCATE(after2)
 ABI_DEALLOCATE(now2)
 ABI_DEALLOCATE(before2)
 ABI_DEALLOCATE(trig3)
 ABI_DEALLOCATE(after3)
 ABI_DEALLOCATE(now3)
 ABI_DEALLOCATE(before3)
 ABI_DEALLOCATE(zmpi2)
 ABI_DEALLOCATE(zw)
 ABI_DEALLOCATE(zt)
 if (nproc_fft>1)  then
   ABI_DEALLOCATE(zmpi1)
 end if

 !call timab(542,2,tsec)

end subroutine sg2002_mpiforw_wf
!!***

!----------------------------------------------------------------------

!!****f* m_sg2002/sg2002_mpifourdp
!! NAME
!! sg2002_mpifourdp
!!
!! FUNCTION
!! Conduct Fourier transform of REAL or COMPLEX function f(r)=fofr defined on
!! fft grid in real space, to create complex f(G)=fofg defined on full fft grid
!! in reciprocal space, in full storage mode, or the reverse operation.
!! For the reverse operation, the final data is divided by nfftot.
!! REAL case when cplex=1, COMPLEX case when cplex=2
!! Usually used for density and potentials.
!!
!! INPUTS
!! cplex=1 if fofr is real, 2 if fofr is complex
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! ndat=Numbre of FFT transforms
!! isign=sign of Fourier transform exponent: current convention uses
!!    +1 for transforming from G to r
!!    -1 for transforming from r to G.
!! fftn2_distrib(2),ffti2_local(2)
!! fftn3_distrib(3),ffti3_local(3)
!! comm_fft=MPI communicator
!!
!! SIDE EFFECTS
!! Input/Output
!! fofg(2,nfft)=f(G), complex.
!! fofr(cplex*nfft)=input function f(r) (real or complex)
!!
!! TODO
!!  Write simplified API for sequential version.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!
!! SOURCE

subroutine sg2002_mpifourdp(cplex,nfft,ngfft,ndat,isign,&
&  fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,isign,nfft,ndat,comm_fft
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: fftn2_distrib(ngfft(2)),ffti2_local(ngfft(2))
 integer,intent(in) :: fftn3_distrib(ngfft(3)),ffti3_local(ngfft(3))
 real(dp),intent(inout) :: fofg(2,nfft*ndat),fofr(cplex*nfft*ndat)

!Local variables-------------------------------
!scalars
 integer :: n1,n2,n3,n4,n5,n6,nd2proc,nd3proc,nproc_fft,me_fft
!arrays
 real(dp),allocatable :: workf(:,:,:,:,:),workr(:,:,:,:,:)

! *************************************************************************

 ! Note the only c2c is supported in parallel.
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6)
 me_fft=ngfft(11); nproc_fft=ngfft(10)

 nd2proc=((n2-1)/nproc_fft) +1
 nd3proc=((n6-1)/nproc_fft) +1
 ABI_ALLOCATE(workr,(2,n4,n5,nd3proc,ndat))
 ABI_ALLOCATE(workf,(2,n4,n6,nd2proc,ndat))

 ! Complex to Complex
 select case (isign)
 case (1)
   ! G --> R
   call mpifft_fg2dbox(nfft,ndat,fofg,n1,n2,n3,n4,nd2proc,n6,fftn2_distrib,ffti2_local,me_fft,workf)

   call sg2002_back(2,ndat,n1,n2,n3,n4,n5,n6,n4,nd2proc,nd3proc,2,workf,workr,comm_fft)

   call mpifft_dbox2fr(n1,n2,n3,n4,n5,nd3proc,ndat,fftn3_distrib,ffti3_local,me_fft,workr,cplex,nfft,fofr)

 case (-1)
   ! R --> G
   call mpifft_fr2dbox(cplex,nfft,ndat,fofr,n1,n2,n3,n4,n5,nd3proc,fftn3_distrib,ffti3_local,me_fft,workr)

   call sg2002_forw(2,ndat,n1,n2,n3,n4,n5,n6,n4,nd2proc,nd3proc,2,workr,workf,comm_fft)

   ! Transfer FFT output to the original fft box.
   call mpifft_dbox2fg(n1,n2,n3,n4,nd2proc,n6,ndat,fftn2_distrib,ffti2_local,me_fft,workf,nfft,fofg)

 case default
   ABI_BUG("Wrong isign")
 end select

 ABI_DEALLOCATE(workr)
 ABI_DEALLOCATE(workf)

end subroutine sg2002_mpifourdp
!!***

!----------------------------------------------------------------------

!!****f* m_sg2002/sg2002_applypot
!! NAME
!!  sg2002_applypot
!!
!! FUNCTION
!! Applies the local real space potential to multiple wavefunctions in Fourier space
!!
!! INPUTS
!!   ZF: Wavefunction (input/output) (note the switch of i2 and i3)
!!        real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!        imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!   max1 is positive or zero ; m1 >=max1+1
!!   i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!!   then, if m1 > max1+1, one has min1=max1-m1+1 and
!!   i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!!   i2 and i3 have a similar definition of range
!!   idat=1,ndat
!!   md1,md2,md3: Dimension of ZF (input as well as output), distributed on different procs
!!   md2proc=((md2-1)/nproc_fft)+1  maximal number of small box 2nd dim slices for one proc
!!
!!   POT: Potential
!!        POT(cplex*i1,i2,i3)
!!        cplex=1 or 2 ,  i1=1,n1 , i2=1,n2 , i3=1,n3
!!   nd1,nd2,nd3: dimension of pot
!!   comm_fft: MPI communicator
!!   nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!   me_fft: [0:nproc_fft-1] number of processor as returned by MPI_COMM_RANK
!!    n1,n2,n3: logical dimension of the transform. As transform lengths
!!              most products of the prime factors 2,3,5 are allowed.
!!             The detailed table with allowed transform lengths can
!!             be found in subroutine CTRIG
!!
!! NOTES:
!!   PERFORMANCE CONSIDERATIONS:
!!   The maximum number of processors that can reasonably be used is max(n2/2,n3/2)
!!
!!   It is very important to find the optimal
!!   value of NCACHE. NCACHE determines the size of the work array ZW, that
!!   has to fit into cache. It has therefore to be chosen to equal roughly
!!    half the size of the physical cache in units of real*8 numbers.
!!   The optimal value of ncache can easily be determined by numerical
!!   experimentation. A too large value of ncache leads to a dramatic
!!   and sudden decrease of performance, a too small value to a to a
!!   slow and less dramatic decrease of performance. If NCACHE is set
!!   to a value so small, that not even a single one dimensional transform
!!   can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!
!! SOURCE


subroutine sg2002_applypot(icplexwf,cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd3proc,&
&  max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,&
&  max1o,max2o,max3o,m1o,m2o,m3o,comm_fft,nproc_fft,me_fft,pot,zf)

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: icplexwf,cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd3proc
 integer,intent(in) :: max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3
 integer,intent(in) :: max1o,max2o,max3o,m1o,m2o,m3o,comm_fft,nproc_fft,me_fft
 real(dp),intent(in) :: pot(cplex*nd1,nd2,nd3)
 real(dp),intent(inout) :: zf(2,md1,md3,md2proc,ndat)

!Local variables-------------------------------
!scalars
 integer,parameter :: unused0=0
 integer :: i,j,i1,i2,i3,ic1,ic2,ic3,idat,ierr,inzee,j3glob
 integer :: ioption,j2,j3,lot,lzt,m1zt,ma,mb,n1dfft,nnd3
 integer :: m2eff,ncache,n1eff,i1inv,i2inv,i3inv,jeff,includelast,j2stb
 integer :: jx,j2stf,Jp2stb,Jp2stf,m2ieff,m2oeff
!arrays
 real(dp) :: tsec(2)
 real(dp), allocatable :: zt(:,:,:) ! work arrays for transpositions
 real(dp), allocatable :: zmpi1(:,:,:,:),zmpi2(:,:,:,:) ! work arrays for MPI
 real(dp), allocatable :: zw(:,:,:) ! cache work array
! FFT work arrays
 real(dp), allocatable, dimension(:,:) :: btrig1,btrig2,btrig3
 real(dp), allocatable, dimension(:,:) :: ftrig1,ftrig2,ftrig3
 integer, allocatable, dimension(:) :: after1,now1,before1,after2,now2,before2,after3,now3,before3

! *************************************************************************

 !ioption=0 ! This was in the old version.
 ioption=1 ! This one is needed to be compatible with paral_kgb

 ncache=4*max(n1,n2,n3,1024)
 if (ncache/(4*max(n1,n2,n3)) < 1) then
   write(std_out,*) &
&    'ncache has to be enlarged to be able to hold at', &
&    'least one 1-d FFT of each size even though this will', &
&    'reduce the performance for shorter transform lengths'
   ABI_ERROR("Aborting now")
 end if

 ! Effective m1 and m2 (complex-to-complex or real-to-complex)
 n1eff=n1; m2ieff=m2i; m2oeff=m2o; m1zt=n1
 if (icplexwf==1) then
   n1eff=(n1+1)/2; m2ieff=m2i/2+1; m2oeff=m2o/2+1; m1zt=2*(n1/2+1)
 end if

 m2eff=max(m2ieff,m2oeff)
 lzt=m2eff
 if (mod(m2eff,2) == 0) lzt=lzt+1
 if (mod(m2eff,4) == 0) lzt=lzt+1

 ! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ABI_ALLOCATE(btrig1,(2,n1))
 ABI_ALLOCATE(ftrig1,(2,n1))
 ABI_ALLOCATE(after1,(mdata))
 ABI_ALLOCATE(now1,(mdata))
 ABI_ALLOCATE(before1,(mdata))
 ABI_ALLOCATE(btrig2,(2,n2))
 ABI_ALLOCATE(ftrig2,(2,n2))
 ABI_ALLOCATE(after2,(mdata))
 ABI_ALLOCATE(now2,(mdata))
 ABI_ALLOCATE(before2,(mdata))
 ABI_ALLOCATE(btrig3,(2,n3))
 ABI_ALLOCATE(ftrig3,(2,n3))
 ABI_ALLOCATE(after3,(mdata))
 ABI_ALLOCATE(now3,(mdata))
 ABI_ALLOCATE(before3,(mdata))

 ABI_ALLOCATE(zw,(2,ncache/4,2))
 ABI_ALLOCATE(zt,(2,lzt,m1zt))
 ABI_ALLOCATE(zmpi2,(2,md1,md2proc,nnd3))
 if (nproc_fft > 1)  then
   ABI_ALLOCATE(zmpi1,(2,md1,md2proc,nnd3))
 end if

 call ctrig(n3,btrig3,after3,before3,now3,1,ic3)
 call ctrig(n1,btrig1,after1,before1,now1,1,ic1)
 call ctrig(n2,btrig2,after2,before2,now2,1,ic2)

 do j=1,n1
   ftrig1(1,j)= btrig1(1,j)
   ftrig1(2,j)=-btrig1(2,j)
 end do
 do j=1,n2
   ftrig2(1,j)= btrig2(1,j)
   ftrig2(2,j)=-btrig2(2,j)
 end do
 do j=1,n3
   ftrig3(1,j)= btrig3(1,j)
   ftrig3(2,j)=-btrig3(2,j)
 end do

 do idat=1,ndat
   !
   ! transform along z axis
   ! input: G1,G3,G2,(Gp2)
   lot=ncache/(4*n3)
   do j2=1,md2proc
     if (me_fft*md2proc+j2 <= m2ieff) then
       do i1=1,m1i,lot
         ma=i1
         mb=min(i1+(lot-1),m1i)
         n1dfft=mb-ma+1

         ! zero-pad n1dfft G_z lines
         ! input: G1,G3,G2,(Gp2)
         call fill_cent(md1,md3,lot,n1dfft,max3i,m3i,n3,zf(1,i1,1,j2,idat),zw(1,1,1))

         inzee=1
         do i=1,ic3
           call fftstp(lot,n1dfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
&                      btrig3,after3(i),now3(i),before3(i),1)
           inzee=3-inzee
         end do

         ! Local rotation.
         ! input:  G1,R3,G2,(Gp2)
         ! output: G1,G2,R3,(Gp2)
         call scramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zw(1,1,inzee),zmpi2)
       end do
     end if
   end do

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Rp3,(Gp2)
   ! output: G1,G2,R3,Gp2,(Rp3)
   if (nproc_fft > 1) then
      call timab(543,1,tsec)
      call xmpi_alltoall(zmpi2,2*md1*md2proc*nd3proc,&
&                        zmpi1,2*md1*md2proc*nd3proc,comm_fft,ierr)
      call timab(543,2,tsec)
   end if

   do j3=1,nd3proc
     j3glob = j3 + me_fft*nd3proc

     if (me_fft*nd3proc+j3 <= n3) then
       Jp2stb=1; J2stb=1
       Jp2stf=1; J2stf=1

       ! transform along x axis
       lot=ncache/(4*n1)

       do j=1,m2ieff,lot
         ma=j
         mb=min(j+(lot-1),m2ieff)
         n1dfft=mb-ma+1

         ! Zero-pad input.
         ! input:  G1,G2,R3,G2,(Rp3)
         ! output: G2,G1,R3,G2,(Rp3)
         if (nproc_fft == 1) then
           call mpiswitch_cent(j3,n1dfft,Jp2stb,J2stb,lot,max1i,md1,m1i,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zmpi2,zw(1,1,1), unused0, unused0, unused0)
         else
           call mpiswitch_cent(j3,n1dfft,Jp2stb,J2stb,lot,max1i,md1,m1i,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zmpi1,zw(1,1,1), unused0, unused0, unused0)
         end if

         ! Transform along x
         ! input:  G2,G1,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         inzee=1
         do i=1,ic1-1
           call fftstp(lot,n1dfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
&                     btrig1,after1(i),now1(i),before1(i),1)
           inzee=3-inzee
         end do

         i=ic1
         call fftstp(lot,n1dfft,n1,lzt,m1zt,zw(1,1,inzee),zt(1,j,1), &
&                    btrig1,after1(i),now1(i),before1(i),1)
       end do

       ! Transform along y axis (take into account c2c or c2r case).
       ! Must loop over the full box.
       lot=ncache/(4*n2)

       if (icplexwf==1) then
         if(mod(lot,2).ne.0)lot=lot-1 ! needed to introduce jeff
       end if

       do j=1,n1eff,lot
         ma=j
         mb=min(j+(lot-1),n1eff)
         n1dfft=mb-ma+1
         jeff=j
         includelast=1

         if (icplexwf==1) then
           jeff=2*j-1
           includelast=1
           if (mb==n1eff .and. n1eff*2/=n1) includelast=0
         end if

         ! Zero-pad the input.
         !  input: G2,R1,R3,(Rp3)
         ! output: R1,G2,R3,(Rp3)
         if (icplexwf==2) then
           call switch_cent(n1dfft,max2i,m2i,n2,lot,n1,lzt,zt(1,1,jeff),zw(1,1,1))
         else
           call switchreal_cent(includelast,n1dfft,max2i,n2,lot,m1zt,lzt,zt(1,1,jeff),zw(1,1,1))
         end if

         ! input:  R1,G2,R3,(Rp3)
         ! output: R1,R2,R3,(Rp3)
         inzee=1
         do i=1,ic2
           call fftstp(lot,n1dfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
&                       btrig2,after2(i),now2(i),before2(i),1)
            inzee=3-inzee
         end do
         ! output: R1,R2,R3,(Rp3)

         ! Multiply with potential in real space
         jx=cplex*(jeff-1)+1
         call multpot(icplexwf,cplex,includelast,nd1,nd2,n2,lot,n1dfft,pot(jx,1,j3glob),zw(1,1,inzee))

         ! TRANSFORM BACK IN FOURIER SPACE
         ! transform along y axis
         ! input: R1,R2,R3,(Rp3)
         do i=1,ic2
           call fftstp(lot,n1dfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
&                       ftrig2,after2(i),now2(i),before2(i),-1)
           inzee=3-inzee
         end do

         !  input: R1,G2,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         if (icplexwf==2) then
           call unswitch_cent(n1dfft,max2o,m2o,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,jeff))
         else
           call unswitchreal_cent(n1dfft,max2o,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,jeff))
         end if

       end do ! j

       ! transform along x axis
       ! input:  R2,R1,R3,(Rp3)
       ! output: R2,G1,R3,(Rp3)
       lot=ncache/(4*n1)

       do j=1,m2oeff,lot
         ma=j
         mb=min(j+(lot-1),m2oeff)
         n1dfft=mb-ma+1
         i=1
         call fftstp(lzt,n1dfft,m1zt,lot,n1,zt(1,j,1),zw(1,1,1), &
&                    ftrig1,after1(i),now1(i),before1(i),-1)

         inzee=1
         do i=2,ic1
           call fftstp(lot,n1dfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
&                       ftrig1,after1(i),now1(i),before1(i),-1)
           inzee=3-inzee
         end do

         ! input:  G2,G1,R3,Gp2,(Rp3)
         ! output: G1,G2,R3,Gp2,(Rp3)
         if (nproc_fft == 1) then
           call unmpiswitch_cent(j3,n1dfft,Jp2stf,J2stf,lot,max1o,md1,m1o,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zw(1,1,inzee),zmpi2)
         else
           call unmpiswitch_cent(j3,n1dfft,Jp2stf,J2stf,lot,max1o,md1,m1o,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zw(1,1,inzee),zmpi1)
         end if
       end do ! j
     end if
   end do

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Gp2,(Rp3)
   ! output: G1,G2,R3,Rp3,(Gp2)
   if (nproc_fft > 1) then
     call timab(544,1,tsec)
     call xmpi_alltoall(zmpi1,2*md1*md2proc*nd3proc, &
&                       zmpi2,2*md1*md2proc*nd3proc,comm_fft,ierr)
     call timab(544,2,tsec)
   end if

   ! transform along z axis
   ! input: G1,G2,R3,(Gp2)
   lot=ncache/(4*n3)
   do j2=1,md2proc
     if (me_fft*md2proc+j2 <= m2oeff) then
       do i1=1,m1o,lot
         ma=i1
         mb=min(i1+(lot-1),m1o)
         n1dfft=mb-ma+1

         ! input:  G1,G2,R3,(Gp2)
         ! output: G1,R3,G2,(Gp2)
         call unscramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zmpi2,zw(1,1,1))

         inzee=1
         do i=1,ic3
           call fftstp(lot,n1dfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
&           ftrig3,after3(i),now3(i),before3(i),-1)
           inzee=3-inzee
         end do

         call unfill_cent(md1,md3,lot,n1dfft,max3o,m3o,n3,zw(1,1,inzee),zf(1,i1,1,j2,idat))
         ! output: G1,G3,G2,(Gp2)
       end do
     end if
   end do

   ! Complete missing values with complex conjugate
   ! Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.
   if (icplexwf==1) then
     do i3=1,m3o
       i3inv=m3o+2-i3
       if (i3==1) i3inv=1
       if (m2oeff>1)then
         do i2=2,m2oeff
           i2inv=m2o+2-i2
           zf(1,1,i3inv,i2inv,idat)= zf(1,1,i3,i2,idat)
           zf(2,1,i3inv,i2inv,idat)=-zf(2,1,i3,i2,idat)
           do i1=2,m1o
             i1inv=m1o+2-i1
             zf(1,i1inv,i3inv,i2inv,idat)= zf(1,i1,i3,i2,idat)
             zf(2,i1inv,i3inv,i2inv,idat)=-zf(2,i1,i3,i2,idat)
           end do
         end do
       end if
     end do
   end if

 end do ! idat

 ABI_DEALLOCATE(btrig1)
 ABI_DEALLOCATE(ftrig1)
 ABI_DEALLOCATE(after1)
 ABI_DEALLOCATE(now1)
 ABI_DEALLOCATE(before1)
 ABI_DEALLOCATE(btrig2)
 ABI_DEALLOCATE(ftrig2)
 ABI_DEALLOCATE(after2)
 ABI_DEALLOCATE(now2)
 ABI_DEALLOCATE(before2)
 ABI_DEALLOCATE(btrig3)
 ABI_DEALLOCATE(ftrig3)
 ABI_DEALLOCATE(after3)
 ABI_DEALLOCATE(now3)
 ABI_DEALLOCATE(before3)

 ABI_DEALLOCATE(zmpi2)
 ABI_DEALLOCATE(zw)
 ABI_DEALLOCATE(zt)
 if (nproc_fft > 1)  then
   ABI_DEALLOCATE(zmpi1)
 end if

end subroutine sg2002_applypot
!!***

!----------------------------------------------------------------------

!!****f* m_sg2002/sg2002_applypot_many
!! NAME
!!  sg2002_applypot_many
!!
!! FUNCTION
!! Applies the local real space potential to multiple wavefunctions in Fourier space
!!
!! INPUTS
!!   ZF: Wavefunction (input/output) (note the switch of i2 and i3)
!!        real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!        imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!   max1 is positive or zero ; m1 >=max1+1
!!   i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!!   then, if m1 > max1+1, one has min1=max1-m1+1 and
!!   i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!!   i2 and i3 have a similar definition of range
!!   idat=1,ndat
!!   md1,md2,md3: Dimension of ZF (input as well as output), distributed on different procs
!!   md2proc=((md2-1)/nproc_fft)+1  maximal number of small box 2nd dim slices for one proc
!!
!!   POT: Potential
!!        POT(cplex*i1,i2,i3)
!!        cplex=1 or 2 ,  i1=1,n1 , i2=1,n2 , i3=1,n3
!!   nd1,nd2,nd3: dimension of pot
!!   comm_fft: MPI communicator
!!   nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!   me_fft: [0:nproc_fft-1] number of processor as returned by MPI_COMM_RANK
!!    n1,n2,n3: logical dimension of the transform. As transform lengths
!!              most products of the prime factors 2,3,5 are allowed.
!!             The detailed table with allowed transform lengths can
!!             be found in subroutine CTRIG
!!
!! NOTES:
!!   PERFORMANCE CONSIDERATIONS:
!!   The maximum number of processors that can reasonably be used is max(n2/2,n3/2)
!!
!!   It is very important to find the optimal
!!   value of NCACHE. NCACHE determines the size of the work array ZW, that
!!   has to fit into cache. It has therefore to be chosen to equal roughly
!!    half the size of the physical cache in units of real*8 numbers.
!!   The optimal value of ncache can easily be determined by numerical
!!   experimentation. A too large value of ncache leads to a dramatic
!!   and sudden decrease of performance, a too small value to a to a
!!   slow and less dramatic decrease of performance. If NCACHE is set
!!   to a value so small, that not even a single one dimensional transform
!!   can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!
!! SOURCE


subroutine sg2002_applypot_many(icplexwf,cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd3proc,&
&  max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,&
&  max1o,max2o,max3o,m1o,m2o,m3o,comm_fft,nproc_fft,me_fft,pot,zf)

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: icplexwf,cplex,ndat,n1,n2,n3,nd1,nd2,nd3,nd3proc
 integer,intent(in) :: max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3
 integer,intent(in) :: max1o,max2o,max3o,m1o,m2o,m3o,comm_fft,nproc_fft,me_fft
 real(dp),intent(in) :: pot(cplex*nd1,nd2,nd3)
 real(dp),intent(inout) :: zf(2,md1,md3,md2proc,ndat)

!Local variables-------------------------------
!scalars
 integer,parameter :: unused0=0
 integer :: i,j,i1,i2,i3,ic1,ic2,ic3,idat,ierr,inzee,j3glob
 integer :: ioption,j2,j3,lot,lzt,m1zt,ma,mb,n1dfft,nnd3
 integer :: m2eff,ncache,n1eff,i1inv,i2inv,i3inv,jeff,includelast,j2stb
 integer :: jx,j2stf,Jp2stb,Jp2stf,m2ieff,m2oeff
!arrays
 integer :: requests(ndat)
 real(dp) :: tsec(2)
 real(dp), allocatable :: zt(:,:,:) ! work arrays for transpositions
 real(dp), allocatable :: zmpi1(:,:,:,:,:),zmpi2(:,:,:,:,:) ! work arrays for MPI
 real(dp), allocatable :: zw(:,:,:) ! cache work array
! FFT work arrays
 real(dp), allocatable, dimension(:,:) :: btrig1,btrig2,btrig3
 real(dp), allocatable, dimension(:,:) :: ftrig1,ftrig2,ftrig3
 integer, allocatable, dimension(:) :: after1,now1,before1,after2,now2,before2,after3,now3,before3

! *************************************************************************

 !ioption=0 ! This was in the old version.
 ioption=1 ! This one is needed to be compatible with paral_kgb

 ! call timab(541,1,tsec)
 ncache=4*max(n1,n2,n3,1024)
 if (ncache/(4*max(n1,n2,n3)) < 1) then
   write(std_out,*) &
&    'ncache has to be enlarged to be able to hold at', &
&    'least one 1-d FFT of each size even though this will', &
&    'reduce the performance for shorter transform lengths'
   ABI_ERROR("Aborting now")
 end if

 ! Effective m1 and m2 (complex-to-complex or real-to-complex)
 n1eff=n1; m2ieff=m2i; m2oeff=m2o; m1zt=n1
 if (icplexwf==1) then
   n1eff=(n1+1)/2; m2ieff=m2i/2+1; m2oeff=m2o/2+1; m1zt=2*(n1/2+1)
 end if

 m2eff=max(m2ieff,m2oeff)
 lzt=m2eff
 if (mod(m2eff,2) == 0) lzt=lzt+1
 if (mod(m2eff,4) == 0) lzt=lzt+1

 ! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ABI_ALLOCATE(btrig1,(2,n1))
 ABI_ALLOCATE(ftrig1,(2,n1))
 ABI_ALLOCATE(after1,(mdata))
 ABI_ALLOCATE(now1,(mdata))
 ABI_ALLOCATE(before1,(mdata))
 ABI_ALLOCATE(btrig2,(2,n2))
 ABI_ALLOCATE(ftrig2,(2,n2))
 ABI_ALLOCATE(after2,(mdata))
 ABI_ALLOCATE(now2,(mdata))
 ABI_ALLOCATE(before2,(mdata))
 ABI_ALLOCATE(btrig3,(2,n3))
 ABI_ALLOCATE(ftrig3,(2,n3))
 ABI_ALLOCATE(after3,(mdata))
 ABI_ALLOCATE(now3,(mdata))
 ABI_ALLOCATE(before3,(mdata))

 ABI_ALLOCATE(zw,(2,ncache/4,2))
 ABI_ALLOCATE(zt,(2,lzt,m1zt))
 ABI_ALLOCATE(zmpi2,(2,md1,md2proc,nnd3,ndat))
 if (nproc_fft > 1)  then
   ABI_ALLOCATE(zmpi1,(2,md1,md2proc,nnd3,ndat))
 end if

 call ctrig(n3,btrig3,after3,before3,now3,1,ic3)
 call ctrig(n1,btrig1,after1,before1,now1,1,ic1)
 call ctrig(n2,btrig2,after2,before2,now2,1,ic2)

 do j=1,n1
   ftrig1(1,j)= btrig1(1,j)
   ftrig1(2,j)=-btrig1(2,j)
 end do
 do j=1,n2
   ftrig2(1,j)= btrig2(1,j)
   ftrig2(2,j)=-btrig2(2,j)
 end do
 do j=1,n3
   ftrig3(1,j)= btrig3(1,j)
   ftrig3(2,j)=-btrig3(2,j)
 end do

 ! Here we take advantage of non-blocking IALLTOALL:
 ! Perform the first step of MPI-FFT for ndat wavefunctions.
 do idat=1,ndat

   !
   ! transform along z axis
   ! input: G1,G3,G2,(Gp2)
   lot=ncache/(4*n3)
   do j2=1,md2proc
     if (me_fft*md2proc+j2 <= m2ieff) then
       do i1=1,m1i,lot
         ma=i1
         mb=min(i1+(lot-1),m1i)
         n1dfft=mb-ma+1

         ! zero-pad n1dfft G_z lines
         ! input: G1,G3,G2,(Gp2)
         call fill_cent(md1,md3,lot,n1dfft,max3i,m3i,n3,zf(1,i1,1,j2,idat),zw(1,1,1))

         inzee=1
         do i=1,ic3
           call fftstp(lot,n1dfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
&                      btrig3,after3(i),now3(i),before3(i),1)
           inzee=3-inzee
         end do

         ! Local rotation.
         ! input:  G1,R3,G2,(Gp2)
         ! output: G1,G2,R3,(Gp2)
         call scramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zw(1,1,inzee),zmpi2(:,:,:,:,idat))
       end do
     end if
   end do

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Rp3,(Gp2)
   ! output: G1,G2,R3,Gp2,(Rp3)
   if (nproc_fft > 1) then
      call timab(543,1,tsec)
      call xmpi_ialltoall(zmpi2(:,:,:,:,idat),2*md1*md2proc*nd3proc,&
&                         zmpi1(:,:,:,:,idat),2*md1*md2proc*nd3proc,comm_fft,requests(idat))
      call timab(543,2,tsec)
   end if
 end do ! idat

 ! The second step of MPI-FFT
 do idat=1,ndat
    ! Make sure communication is completed.
    if (nproc_fft>1) call xmpi_wait(requests(idat),ierr)

   do j3=1,nd3proc
     j3glob = j3 + me_fft*nd3proc

     if (me_fft*nd3proc+j3 <= n3) then
       Jp2stb=1; J2stb=1
       Jp2stf=1; J2stf=1

       ! transform along x axis
       lot=ncache/(4*n1)

       do j=1,m2ieff,lot
         ma=j
         mb=min(j+(lot-1),m2ieff)
         n1dfft=mb-ma+1

         ! Zero-pad input.
         ! input:  G1,G2,R3,G2,(Rp3)
         ! output: G2,G1,R3,G2,(Rp3)
         if (nproc_fft == 1) then
           call mpiswitch_cent(j3,n1dfft,Jp2stb,J2stb,lot,max1i,md1,m1i,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zmpi2(:,:,:,:,idat),zw(1,1,1), unused0, unused0, unused0)
         else
           call mpiswitch_cent(j3,n1dfft,Jp2stb,J2stb,lot,max1i,md1,m1i,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zmpi1(:,:,:,:,idat),zw(1,1,1), unused0, unused0, unused0)
         end if

         ! Transform along x
         ! input:  G2,G1,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         inzee=1
         do i=1,ic1-1
           call fftstp(lot,n1dfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
&                     btrig1,after1(i),now1(i),before1(i),1)
           inzee=3-inzee
         end do

         i=ic1
         call fftstp(lot,n1dfft,n1,lzt,m1zt,zw(1,1,inzee),zt(1,j,1), &
&                    btrig1,after1(i),now1(i),before1(i),1)
       end do

       ! Transform along y axis (take into account c2c or c2r case).
       ! Must loop over the full box.
       lot=ncache/(4*n2)

       if (icplexwf==1) then
         if(mod(lot,2).ne.0)lot=lot-1 ! needed to introduce jeff
       end if

       do j=1,n1eff,lot
         ma=j
         mb=min(j+(lot-1),n1eff)
         n1dfft=mb-ma+1
         jeff=j
         includelast=1

         if (icplexwf==1) then
           jeff=2*j-1
           includelast=1
           if (mb==n1eff .and. n1eff*2/=n1) includelast=0
         end if

         ! Zero-pad the input.
         !  input: G2,R1,R3,(Rp3)
         ! output: R1,G2,R3,(Rp3)
         if (icplexwf==2) then
           call switch_cent(n1dfft,max2i,m2i,n2,lot,n1,lzt,zt(1,1,jeff),zw(1,1,1))
         else
           call switchreal_cent(includelast,n1dfft,max2i,n2,lot,m1zt,lzt,zt(1,1,jeff),zw(1,1,1))
         end if

         ! input:  R1,G2,R3,(Rp3)
         ! output: R1,R2,R3,(Rp3)
         inzee=1
         do i=1,ic2
           call fftstp(lot,n1dfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
&                       btrig2,after2(i),now2(i),before2(i),1)
            inzee=3-inzee
         end do
         ! output: R1,R2,R3,(Rp3)

         ! Multiply with potential in real space
         jx=cplex*(jeff-1)+1
         call multpot(icplexwf,cplex,includelast,nd1,nd2,n2,lot,n1dfft,pot(jx,1,j3glob),zw(1,1,inzee))

         ! TRANSFORM BACK IN FOURIER SPACE
         ! transform along y axis
         ! input: R1,R2,R3,(Rp3)
         do i=1,ic2
           call fftstp(lot,n1dfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
&                       ftrig2,after2(i),now2(i),before2(i),-1)
           inzee=3-inzee
         end do

         !  input: R1,G2,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         if (icplexwf==2) then
           call unswitch_cent(n1dfft,max2o,m2o,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,jeff))
         else
           call unswitchreal_cent(n1dfft,max2o,n2,lot,n1,lzt,zw(1,1,inzee),zt(1,1,jeff))
         end if

       end do ! j

       ! transform along x axis
       ! input:  R2,R1,R3,(Rp3)
       ! output: R2,G1,R3,(Rp3)
       lot=ncache/(4*n1)

       do j=1,m2oeff,lot
         ma=j
         mb=min(j+(lot-1),m2oeff)
         n1dfft=mb-ma+1
         i=1
         call fftstp(lzt,n1dfft,m1zt,lot,n1,zt(1,j,1),zw(1,1,1), &
&                    ftrig1,after1(i),now1(i),before1(i),-1)

         inzee=1
         do i=2,ic1
           call fftstp(lot,n1dfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
&                       ftrig1,after1(i),now1(i),before1(i),-1)
           inzee=3-inzee
         end do

         ! input:  G2,G1,R3,Gp2,(Rp3)
         ! output: G1,G2,R3,Gp2,(Rp3)
         if (nproc_fft == 1) then
           call unmpiswitch_cent(j3,n1dfft,Jp2stf,J2stf,lot,max1o,md1,m1o,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zw(1,1,inzee),zmpi2(:,:,:,:,idat))
         else
           call unmpiswitch_cent(j3,n1dfft,Jp2stf,J2stf,lot,max1o,md1,m1o,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zw(1,1,inzee),zmpi1(:,:,:,:,idat))
         end if
       end do ! j
     end if
   end do

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Gp2,(Rp3)
   ! output: G1,G2,R3,Rp3,(Gp2)
   if (nproc_fft > 1) then
     call timab(544,1,tsec)
     call xmpi_ialltoall(zmpi1(:,:,:,:,idat),2*md1*md2proc*nd3proc, &
&                        zmpi2(:,:,:,:,idat),2*md1*md2proc*nd3proc,comm_fft,requests(idat))
     call timab(544,2,tsec)
   end if
 end do ! idat

 do idat=1,ndat
   if (nproc_fft>1) call xmpi_wait(requests(idat),ierr)

   ! transform along z axis
   ! input: G1,G2,R3,(Gp2)
   lot=ncache/(4*n3)
   do j2=1,md2proc
     if (me_fft*md2proc+j2 <= m2oeff) then
       do i1=1,m1o,lot
         ma=i1
         mb=min(i1+(lot-1),m1o)
         n1dfft=mb-ma+1

         ! input:  G1,G2,R3,(Gp2)
         ! output: G1,R3,G2,(Gp2)
         call unscramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zmpi2(:,:,:,:,idat),zw(1,1,1))

         inzee=1
         do i=1,ic3
           call fftstp(lot,n1dfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
&           ftrig3,after3(i),now3(i),before3(i),-1)
           inzee=3-inzee
         end do

         call unfill_cent(md1,md3,lot,n1dfft,max3o,m3o,n3,zw(1,1,inzee),zf(1,i1,1,j2,idat))
         ! output: G1,G3,G2,(Gp2)
       end do
     end if
   end do

   ! Complete missing values with complex conjugate
   ! Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.
   if (icplexwf==1) then
     do i3=1,m3o
       i3inv=m3o+2-i3
       if (i3==1) i3inv=1
       if (m2oeff>1)then
         do i2=2,m2oeff
           i2inv=m2o+2-i2
           zf(1,1,i3inv,i2inv,idat)= zf(1,1,i3,i2,idat)
           zf(2,1,i3inv,i2inv,idat)=-zf(2,1,i3,i2,idat)
           do i1=2,m1o
             i1inv=m1o+2-i1
             zf(1,i1inv,i3inv,i2inv,idat)= zf(1,i1,i3,i2,idat)
             zf(2,i1inv,i3inv,i2inv,idat)=-zf(2,i1,i3,i2,idat)
           end do
         end do
       end if
     end do
   end if

 end do ! idat

 ABI_DEALLOCATE(btrig1)
 ABI_DEALLOCATE(ftrig1)
 ABI_DEALLOCATE(after1)
 ABI_DEALLOCATE(now1)
 ABI_DEALLOCATE(before1)
 ABI_DEALLOCATE(btrig2)
 ABI_DEALLOCATE(ftrig2)
 ABI_DEALLOCATE(after2)
 ABI_DEALLOCATE(now2)
 ABI_DEALLOCATE(before2)
 ABI_DEALLOCATE(btrig3)
 ABI_DEALLOCATE(ftrig3)
 ABI_DEALLOCATE(after3)
 ABI_DEALLOCATE(now3)
 ABI_DEALLOCATE(before3)

 ABI_DEALLOCATE(zmpi2)
 ABI_DEALLOCATE(zw)
 ABI_DEALLOCATE(zt)
 if (nproc_fft > 1)  then
   ABI_DEALLOCATE(zmpi1)
 end if

end subroutine sg2002_applypot_many
!!***

!----------------------------------------------------------------------

!!****f* m_sg2002/sg2002_accrho
!! NAME
!! sg2002_accrho
!!
!! FUNCTION
!! Accumulates the real space density rho from the ndat wavefunctions zf
!! by transforming zf into real space and adding all the amplitudes squared
!!
!! INPUTS:
!!   ZF: input array (note the switch of i2 and i3)
!!         real(F(i1,i3,i2,idat))=ZF(1,i1,i3,i2,idat)
!!         imag(F(i1,i3,i2,idat))=ZF(2,i1,i3,i2,idat)
!!   max1 is positive or zero ; m1 >=max1+1
!!   i1= 1... max1+1 corresponds to positive and zero wavevectors 0 ... max1
!!   then, if m1 > max1+1, one has min1=max1-m1+1 and
!!   i1= max1+2 ... m1 corresponds to negative wavevectors min1 ... -1
!!   i2 and i3 have a similar definition of range
!!   idat=1,ndat
!!   md1,md2,md3: Dimension of ZF
!!   md2proc=((md2-1)/nproc_fft)+1 ! maximal number of small box 2nd dim slices for one proc
!!   weight(ndat)= weight for the density accumulation
!!
!! OUTPUTS:
!!    RHOoutput(i1,i2,i3) = RHOinput(i1,i2,i3) + sum on idat of (Re(FFT(ZF))**2 *weight_r + weight_i*Im(FFT(ZF))**2
!!        i1=1,n1 , i2=1,n2 , i3=1,n3
!!   comm_fft: MPI communicator
!!   nproc_fft: number of processors used as returned by MPI_COMM_SIZE
!!   me_fft: [0:nproc_fft-1] number of processor as returned by MPI_COMM_RANK
!!    n1,n2,n3: logical dimension of the transform. As transform lengths
!!              most products of the prime factors 2,3,5 are allowed.
!!             The detailed table with allowed transform lengths can
!!             be found in subroutine CTRIG
!!    nd1,nd2,nd3: Dimension of RHO
!!   nd3proc=((nd3-1)/nproc_fft)+1 ! maximal number of big box 3rd dim slices for one proc
!!
!! NOTES:
!!   PERFORMANCE CONSIDERATIONS:
!!   The maximum number of processors that can reasonably be used is max(n2/2,n3/2)
!!
!!   It is very important to find the optimal
!!   value of NCACHE. NCACHE determines the size of the work array ZW, that
!!   has to fit into cache. It has therefore to be chosen to equal roughly
!!    half the size of the physical cache in units of real*8 numbers.
!!   The optimal value of ncache can easily be determined by numerical
!!   experimentation. A too large value of ncache leads to a dramatic
!!   and sudden decrease of performance, a too small value to a to a
!!   slow and less dramatic decrease of performance. If NCACHE is set
!!   to a value so small, that not even a single one dimensional transform
!!   can be done in the workarray zw, the program stops with an error message.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!
!! SOURCE

subroutine sg2002_accrho(icplexwf,ndat,n1,n2,n3,nd1,nd2,nd3,nd3proc,&
&  max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft,nproc_fft,me_fft,zf,rho,weight_r,weight_i)

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: icplexwf,ndat,n1,n2,n3,nd1,nd2,nd3,nd3proc
 integer,intent(in) :: max1,max2,max3,m1,m2,m3,md1,md2proc,md3,comm_fft,nproc_fft,me_fft
 real(dp),intent(in) :: zf(2,md1,md3,md2proc,ndat)
 real(dp),intent(in) :: weight_r(ndat), weight_i(ndat)
 real(dp),intent(inout) :: rho(nd1,nd2,nd3)

!Local variables-------------------------------
!scalars
 integer,parameter :: unused0=0
 integer :: i,j,i1,ic1,ic2,ic3,idat,ierr,inzee,j3glob
 integer :: ioption,j2,j3,j2st,jp2st,lot,lzt,m1zt,ma,mb,n1dfft,nnd3
 integer :: m2eff,ncache,n1eff,jeff,includelast
!arrays
 real(dp), allocatable :: zmpi1(:,:,:,:),zmpi2(:,:,:,:) ! work arrays for MPI
 real(dp), allocatable :: zt(:,:,:)  ! work arrays for transpositions
 real(dp), allocatable :: zw(:,:,:) ! cache work array
 real(dp) :: tsec(2)
! FFT work arrays
 real(dp), allocatable, dimension(:,:) :: trig1,trig2,trig3
 integer, allocatable, dimension(:) :: after1,now1,before1, after2,now2,before2,after3,now3,before3

! *************************************************************************

 !ioption=0 ! This was in the old version.
 ioption=1 ! This one is needed to be compatible with paral_kgb

 !nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)

! find cache size that gives optimal performance on machine
 ncache=4*max(n1,n2,n3,1024)
 if (ncache/(4*max(n1,n2,n3)) < 1) then
    write(std_out,*) &
&     'ncache has to be enlarged to be able to hold at', &
&     'least one 1-d FFT of each size even though this will', &
&     'reduce the performance for shorter transform lengths'
    ABI_ERROR("Aborting now")
 end if

!Effective m1 and m2 (complex-to-complex or real-to-complex)
 n1eff=n1; m2eff=m2 ; m1zt=n1
 if (icplexwf==1) then
   n1eff=(n1+1)/2; m2eff=m2/2+1; m1zt=2*(n1/2+1)
 end if

 lzt=m2eff
 if (mod(m2eff,2) == 0) lzt=lzt+1
 if (mod(m2eff,4) == 0) lzt=lzt+1

 ! maximal number of big box 3rd dim slices for all procs
 nnd3=nd3proc*nproc_fft

 ABI_ALLOCATE(trig1,(2,n1))
 ABI_ALLOCATE(after1,(mdata))
 ABI_ALLOCATE(now1,(mdata))
 ABI_ALLOCATE(before1,(mdata))
 ABI_ALLOCATE(trig2,(2,n2))
 ABI_ALLOCATE(after2,(mdata))
 ABI_ALLOCATE(now2,(mdata))
 ABI_ALLOCATE(before2,(mdata))
 ABI_ALLOCATE(trig3,(2,n3))
 ABI_ALLOCATE(after3,(mdata))
 ABI_ALLOCATE(now3,(mdata))
 ABI_ALLOCATE(before3,(mdata))

 ABI_ALLOCATE(zw,(2,ncache/4,2))
 ABI_ALLOCATE(zt,(2,lzt,m1zt))
 ABI_ALLOCATE(zmpi2,(2,md1,md2proc,nnd3))
 if (nproc_fft > 1)  then
   ABI_ALLOCATE(zmpi1,(2,md1,md2proc,nnd3))
 end if

 call ctrig(n3,trig3,after3,before3,now3,1,ic3)
 call ctrig(n1,trig1,after1,before1,now1,1,ic1)
 call ctrig(n2,trig2,after2,before2,now2,1,ic2)

 do idat=1,ndat
   ! transform along z axis
   ! input: I1,I3,J2,(Jp2)
   lot=ncache/(4*n3)

   ! Loop over the y planes treated by this node and trasform n1ddft G_z lines.
   do j2=1,md2proc
     if (me_fft*md2proc+j2 <= m2eff) then ! MG REMOVED TO BE COSISTENT WITH BACK_WF
       do i1=1,m1,lot
         ma=i1
         mb=min(i1+(lot-1),m1)
         n1dfft=mb-ma+1

         ! zero-pad n1dfft G_z lines
         !  input: G1,G3,G2,(Gp2)
         ! output: G1,R3,G2,(Gp2)
         call fill_cent(md1,md3,lot,n1dfft,max3,m3,n3,zf(1,i1,1,j2,idat),zw(1,1,1))

         ! Transform along z.
         inzee=1
         do i=1,ic3
           call fftstp(lot,n1dfft,n3,lot,n3,zw(1,1,inzee),zw(1,1,3-inzee), &
&                       trig3,after3(i),now3(i),before3(i),1)
           inzee=3-inzee
         end do

         ! Local rotation.
         ! input:  G1,R3,G2,(Gp2)
         ! output: G1,G2,R3,(Gp2)
         call scramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zw(1,1,inzee),zmpi2)
       end do
     end if
   end do

   ! Interprocessor data transposition
   ! input:  G1,G2,R3,Rp3,(Gp2)
   ! output: G1,G2,R3,Gp2,(Rp3)
   if (nproc_fft > 1) then
     call timab(543,1,tsec)
     call xmpi_alltoall(zmpi2,2*md1*md2proc*nd3proc, &
&                       zmpi1,2*md1*md2proc*nd3proc,comm_fft,ierr)
     call timab(543,2,tsec)
   end if

   ! Loop over the z treated by this node.
   do j3=1,nd3proc
     j3glob = j3 + me_fft*nd3proc
     !ABI_CHECK(j3glob <= n3, "j3glob")

     if (me_fft*nd3proc+j3 <= n3) then
       Jp2st=1; J2st=1

       lot=ncache/(4*n1)

       ! Loop over G_y in the small box.
       do j=1,m2eff,lot
         ma=j
         mb=min(j+(lot-1),m2eff)
         n1dfft=mb-ma+1

         ! Zero-pad input.
         ! input:  G1,G2,R3,JG2,(Rp3)
         ! output: G2,G1,R3,JG2,(Rp3)

         if (nproc_fft == 1) then
          call mpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot,max1,md1,m1,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zmpi2,zw(1,1,1),unused0, unused0,unused0)
         else
          call mpiswitch_cent(j3,n1dfft,Jp2st,J2st,lot,max1,md1,m1,n1,&
&           md2proc,nd3proc,nproc_fft,ioption,zmpi1,zw(1,1,1), unused0,unused0,unused0)
         end if

         ! Transform along x
         ! input:  G2,G1,R3,(Rp3)
         ! output: G2,R1,R3,(Rp3)
         inzee=1
         do i=1,ic1-1
           call fftstp(lot,n1dfft,n1,lot,n1,zw(1,1,inzee),zw(1,1,3-inzee), &
&                    trig1,after1(i),now1(i),before1(i),1)
           inzee=3-inzee
         end do

         i=ic1
         call fftstp(lot,n1dfft,n1,lzt,m1zt,zw(1,1,inzee),zt(1,j,1), &
&                    trig1,after1(i),now1(i),before1(i),1)
       end do

       ! Transform along y axis (take into account c2c or c2r case).
       ! Must loop over the full box.
       lot=ncache/(4*n2)
       if (icplexwf==1) then
         if (mod(lot,2).ne.0) lot=lot-1 ! needed to introduce jeff
       end if

       do j=1,n1eff,lot
         ma=j
         mb=min(j+(lot-1),n1eff)
         n1dfft=mb-ma+1
         jeff=j
         includelast=1

         if (icplexwf==1) then
           jeff=2*j-1
           includelast=1
           if (mb==n1eff .and. n1eff*2/=n1) includelast=0
         end if

         ! Zero-pad the input.
         ! input:  G2,R1,R3,(Rp3)
         ! output: R1,G2,R3,(Rp3)
         if (icplexwf==2) then
           call switch_cent(n1dfft,max2,m2,n2,lot,n1,lzt,zt(1,1,j),zw(1,1,1))
         else
           call switchreal_cent(includelast,n1dfft,max2,n2,lot,m1zt,lzt,zt(1,1,jeff),zw(1,1,1))
         end if

         inzee=1
         do i=1,ic2
           call fftstp(lot,n1dfft,n2,lot,n2,zw(1,1,inzee),zw(1,1,3-inzee), &
&                       trig2,after2(i),now2(i),before2(i),1)
           inzee=3-inzee
         end do

         ! Accumulate
         call addrho(icplexwf,includelast,nd1,nd2,n2,lot,n1dfft,&
&          zw(1,1,inzee),rho(jeff,1,j3glob),weight_r(idat),weight_i(idat))
       end do
       ! output: i1,i2,j3,(jp3)

      end if
    end do ! j3
 end do ! idat

 ABI_DEALLOCATE(trig1)
 ABI_DEALLOCATE(after1)
 ABI_DEALLOCATE(now1)
 ABI_DEALLOCATE(before1)
 ABI_DEALLOCATE(trig2)
 ABI_DEALLOCATE(after2)
 ABI_DEALLOCATE(now2)
 ABI_DEALLOCATE(before2)
 ABI_DEALLOCATE(trig3)
 ABI_DEALLOCATE(after3)
 ABI_DEALLOCATE(now3)
 ABI_DEALLOCATE(before3)

 ABI_DEALLOCATE(zmpi2)
 ABI_DEALLOCATE(zw)
 ABI_DEALLOCATE(zt)
 if (nproc_fft > 1)  then
   ABI_DEALLOCATE(zmpi1)
 end if

end subroutine sg2002_accrho
!!***

!!****f* m_sg2002/ctrig
!! NAME
!!  ctrig
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_sg2002
!!
!! CHILDREN
!!
!! SOURCE

subroutine ctrig(n,trig,after,before,now,isign,ic)

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: n,isign
 integer,intent(inout) :: ic
 integer,intent(inout) :: after(mdata),before(mdata),now(mdata)
 real(dp),intent(inout) :: trig(2,n)

!Local variables-------------------------------
!scalars
 integer :: i,itt,j,nh
 real(dp) :: angle,trigc,trigs

! *************************************************************************

 do i=1,ndata
   if (n.eq.ifftdata(1,i)) then
     ic=0
     do j=1,(mdata-1)
       itt=ifftdata(1+j,i)
       if (itt.gt.1) then
         ic=ic+1
         now(j)=ifftdata(1+j,i)
       else
         goto 1000
       end if
     end do
     goto 1000
   end if
 end do

 write(std_out,*) 'VALUE OF',n,'NOT ALLOWED FOR FFT, ALLOWED VALUES ARE:'
37 format(15(i5))
 write(std_out,37) (ifftdata(1,i),i=1,ndata)
 ABI_ERROR("Aborting now")

1000 continue
 after(1)=1
 before(ic)=1
 do i=2,ic
   after(i)=after(i-1)*now(i-1)
   before(ic-i+1)=before(ic-i+2)*now(ic-i+2)
 end do

 angle=isign*two_pi/n
 if (mod(n,2).eq.0) then
   nh=n/2
   trig(1,1)=one
   trig(2,1)=zero
   trig(1,nh+1)=-one
   trig(2,nh+1)=zero
   do i=1,nh-1
     trigc=cos(i*angle)
     trigs=sin(i*angle)
     trig(1,i+1)=trigc
     trig(2,i+1)=trigs
     trig(1,n-i+1)=trigc
     trig(2,n-i+1)=-trigs
   end do
 else
   nh=(n-1)/2
   trig(1,1)=one
   trig(2,1)=zero
   do i=1,nh
     trigc=cos(i*angle)
     trigs=sin(i*angle)
     trig(1,i+1)=trigc
     trig(2,i+1)=trigs
     trig(1,n-i+1)=trigc
     trig(2,n-i+1)=-trigs
   end do
 end if

end subroutine ctrig
!!***

!!****f* m_sg2002/fftstp
!! NAME
!!  fftstp
!!
!! FUNCTION
!!
!! INPUTS
!!   mm
!!   n1dfft
!!   m
!!   nn
!!   n
!!   zin
!!   trig
!!   after
!!   now
!!   before
!!   isign
!!
!! OUTPUT
!!   zout
!!
!! PARENTS
!!      m_sg2002
!!
!! CHILDREN
!!
!! SOURCE

subroutine fftstp(mm,n1dfft,m,nn,n,zin,zout,trig,after,now,before,isign)

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: after,before,mm,n1dfft,m,nn,n,now,isign
 real(dp),intent(in) :: trig(2,n),zin(2,mm,m)
 real(dp),intent(inout) :: zout(2,nn,n)

!Local variables-------------------------------
 integer :: atn,atb,ia,ias,ib,itrig,itt,j,nin1,nin2,nin3,nin4,nin5,nin6,nin7,nin8
 integer :: nout1,nout2,nout3,nout4,nout5,nout6,nout7,nout8
 real(dp) :: am,ap,bm,bp,ci3,ci4,ci5,ci6,ci7,ci8,cm,cos2,cos4,cp,cr2,cr3,cr4,cr5,cr6,cr7,cr8
 real(dp) :: dm,bb,ci2,dpp,r,r2,r25,r3,r34,r4,r5,r6,r7,r8,rt2i,s,r1,s1,s2,s3,s25,s34,s4,s5,s6,s7,s8
 real(dp) :: sin2,ui1,ui2,ui3,ur1,ur2,ur3,sin4,vi1,vi2,vi3,vr1,vr2,vr3

! *************************************************************************
        atn=after*now
        atb=after*before

!         sqrt(.5d0)
        rt2i=half_sqrt2
        if (now.eq.2) then
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib=1,before
          nin1=nin1+after
          nin2=nin1+atb
          nout1=nout1+atn
          nout2=nout1+after
            do j=1,n1dfft
            r1=zin(1,j,nin1)
            s1=zin(2,j,nin1)
            r2=zin(1,j,nin2)
            s2=zin(2,j,nin2)
            zout(1,j,nout1)= r2 + r1
            zout(2,j,nout1)= s2 + s1
            zout(1,j,nout2)= r1 - r2
            zout(2,j,nout2)= s1 - s2
          enddo
        enddo
        do 2000,ia=2,after
        ias=ia-1
        if (2*ias.eq.after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                          nin1=nin1+after
                          nin2=nin1+atb
                          nout1=nout1+atn
                          nout2=nout1+after
                          do j=1,n1dfft
                            r1=zin(1,j,nin1)
                            s1=zin(2,j,nin1)
                            r2=zin(2,j,nin2)
                            s2=zin(1,j,nin2)
                            zout(1,j,nout1)= r1 - r2
                            zout(2,j,nout1)= s2 + s1
                            zout(1,j,nout2)= r2 + r1
                            zout(2,j,nout2)= s1 - s2
                          enddo
                        enddo
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                          nin1=nin1+after
                          nin2=nin1+atb
                          nout1=nout1+atn
                          nout2=nout1+after
                            do j=1,n1dfft
                            r1=zin(1,j,nin1)
                            s1=zin(2,j,nin1)
                            r2=zin(2,j,nin2)
                            s2=zin(1,j,nin2)
                            zout(1,j,nout1)= r2 + r1
                            zout(2,j,nout1)= s1 - s2
                            zout(1,j,nout2)= r1 - r2
                            zout(2,j,nout2)= s2 + s1
                          enddo
                        enddo
                end if
        else if (4*ias.eq.after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do j=1,n1dfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r - s)*rt2i
                        s2=(r + s)*rt2i
                        zout(1,j,nout1)= r2 + r1
                        zout(2,j,nout1)= s2 + s1
                        zout(1,j,nout2)= r1 - r2
                        zout(2,j,nout2)= s1 - s2
                        enddo
                        enddo
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do j=1,n1dfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(s - r)*rt2i
                        zout(1,j,nout1)= r2 + r1
                        zout(2,j,nout1)= s2 + s1
                        zout(1,j,nout2)= r1 - r2
                        zout(2,j,nout2)= s1 - s2
                        enddo
                        enddo
                end if
        else if (4*ias.eq.3*after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do j=1,n1dfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(r - s)*rt2i
                        zout(1,j,nout1)= r1 - r2
                        zout(2,j,nout1)= s2 + s1
                        zout(1,j,nout2)= r2 + r1
                        zout(2,j,nout2)= s1 - s2
                        enddo
                        enddo
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        do j=1,n1dfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(s - r)*rt2i
                        s2=(r + s)*rt2i
                        zout(1,j,nout1)= r2 + r1
                        zout(2,j,nout1)= s1 - s2
                        zout(1,j,nout2)= r1 - r2
                        zout(2,j,nout2)= s2 + s1
                        enddo
                        enddo
                end if
        else
                itrig=ias*before+1
                cr2=trig(1,itrig)
                ci2=trig(2,itrig)
                nin1=ia-after
                nout1=ia-atn
                do ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nout1=nout1+atn
                nout2=nout1+after
                do j=1,n1dfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=r*cr2 - s*ci2
                s2=r*ci2 + s*cr2
                zout(1,j,nout1)= r2 + r1
                zout(2,j,nout1)= s2 + s1
                zout(1,j,nout2)= r1 - r2
                zout(2,j,nout2)= s1 - s2
                enddo
                enddo
        end if
2000        continue
        else if (now.eq.4) then
        if (isign.eq.1) then
                ia=1
                nin1=ia-after
                nout1=ia-atn
                do ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                do j=1,n1dfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(1,j,nin2)
                s2=zin(2,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r4=zin(1,j,nin4)
                s4=zin(2,j,nin4)
                r=r1 + r3
                s=r2 + r4
                zout(1,j,nout1) = r + s
                zout(1,j,nout3) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,j,nout2) = r - s
                zout(1,j,nout4) = r + s
                r=s1 + s3
                s=s2 + s4
                zout(2,j,nout1) = r + s
                zout(2,j,nout3) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,j,nout2) = r + s
                zout(2,j,nout4) = r - s
                enddo
                enddo
                do 4000,ia=2,after
                ias=ia-1
                if (2*ias.eq.after) then
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do j=1,n1dfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r-s)*rt2i
                        s2=(r+s)*rt2i
                        r3=zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(r + s)*rt2i
                        s4=(r - s)*rt2i
                        r=r1 - r3
                        s=r2 - r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 + r3
                        s=s2 - s4
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout4) = r + s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 - s3
                        s=r2 + r4
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout4) = r - s
                        enddo
                        enddo
                else
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do j=1,n1dfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout4) = r + s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout4) = r - s
                        enddo
                        enddo
                end if
4000                continue
        else
                ia=1
                nin1=ia-after
                nout1=ia-atn
                do ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                do j=1,n1dfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(1,j,nin2)
                s2=zin(2,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r4=zin(1,j,nin4)
                s4=zin(2,j,nin4)
                r=r1 + r3
                s=r2 + r4
                zout(1,j,nout1) = r + s
                zout(1,j,nout3) = r - s
                r=r1 - r3
                s=s2 - s4
                zout(1,j,nout2) = r + s
                zout(1,j,nout4) = r - s
                r=s1 + s3
                s=s2 + s4
                zout(2,j,nout1) = r + s
                zout(2,j,nout3) = r - s
                r=s1 - s3
                s=r2 - r4
                zout(2,j,nout2) = r - s
                zout(2,j,nout4) = r + s
                enddo
                enddo
                do 4100,ia=2,after
                ias=ia-1
                if (2*ias.eq.after) then
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do j=1,n1dfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(s - r)*rt2i
                        r3=zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(s - r)*rt2i
                        s4=(r + s)*rt2i
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 - r3
                        s=s2 + s4
                        zout(1,j,nout2) = r + s
                        zout(1,j,nout4) = r - s
                        r=s1 - s3
                        s=s2 - s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 + s3
                        s=r2 - r4
                        zout(2,j,nout2) = r - s
                        zout(2,j,nout4) = r + s
                        enddo
                        enddo
                else
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        do j=1,n1dfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=r1 + r3
                        s=r2 + r4
                        zout(1,j,nout1) = r + s
                        zout(1,j,nout3) = r - s
                        r=r1 - r3
                        s=s2 - s4
                        zout(1,j,nout2) = r + s
                        zout(1,j,nout4) = r - s
                        r=s1 + s3
                        s=s2 + s4
                        zout(2,j,nout1) = r + s
                        zout(2,j,nout3) = r - s
                        r=s1 - s3
                        s=r2 - r4
                        zout(2,j,nout2) = r - s
                        zout(2,j,nout4) = r + s
                        enddo
                        enddo
                end if
4100                continue
        end if
        else if (now.eq.8) then
        if (isign.eq.-1) then
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        do j=1,n1dfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(1,j,nin2)
                        s2=zin(2,j,nin2)
                        r3=zin(1,j,nin3)
                        s3=zin(2,j,nin3)
                        r4=zin(1,j,nin4)
                        s4=zin(2,j,nin4)
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r6=zin(1,j,nin6)
                        s6=zin(2,j,nin6)
                        r7=zin(1,j,nin7)
                        s7=zin(2,j,nin7)
                        r8=zin(1,j,nin8)
                        s8=zin(2,j,nin8)
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dpp=r + s
                        dm=r - s
                        zout(1,j,nout1) = ap + bp
                        zout(2,j,nout1) = cp + dpp
                        zout(1,j,nout5) = ap - bp
                        zout(2,j,nout5) = cp - dpp
                        zout(1,j,nout3) = am + dm
                        zout(2,j,nout3) = cm - bm
                        zout(1,j,nout7) = am - dm
                        zout(2,j,nout7) = cm + bm
                        r=r1 - r5
                        s=s3 - s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r3 - r7
                        bp=r + s
                        bm=r - s
                        r=s4 - s8
                        s=r2 - r6
                        cp=r + s
                        cm=r - s
                        r=s2 - s6
                        s=r4 - r8
                        dpp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( dm - cp)*rt2i
                        cp= ( cm + dpp)*rt2i
                        dpp = ( cm - dpp)*rt2i
                        zout(1,j,nout2) = ap + r
                        zout(2,j,nout2) = bm + s
                        zout(1,j,nout6) = ap - r
                        zout(2,j,nout6) = bm - s
                        zout(1,j,nout4) = am + cp
                        zout(2,j,nout4) = bp + dpp
                        zout(1,j,nout8) = am - cp
                        zout(2,j,nout8) = bp - dpp
                        enddo
                        enddo
                do 8000,ia=2,after
                ias=ia-1
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        itrig=itrig+itt
                        cr5=trig(1,itrig)
                        ci5=trig(2,itrig)
                        itrig=itrig+itt
                        cr6=trig(1,itrig)
                        ci6=trig(2,itrig)
                        itrig=itrig+itt
                        cr7=trig(1,itrig)
                        ci7=trig(2,itrig)
                        itrig=itrig+itt
                        cr8=trig(1,itrig)
                        ci8=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        do j=1,n1dfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=zin(1,j,nin5)
                        s=zin(2,j,nin5)
                        r5=r*cr5 - s*ci5
                        s5=r*ci5 + s*cr5
                        r=zin(1,j,nin6)
                        s=zin(2,j,nin6)
                        r6=r*cr6 - s*ci6
                        s6=r*ci6 + s*cr6
                        r=zin(1,j,nin7)
                        s=zin(2,j,nin7)
                        r7=r*cr7 - s*ci7
                        s7=r*ci7 + s*cr7
                        r=zin(1,j,nin8)
                        s=zin(2,j,nin8)
                        r8=r*cr8 - s*ci8
                        s8=r*ci8 + s*cr8
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dpp=r + s
                        dm=r - s
                        zout(1,j,nout1) = ap + bp
                        zout(2,j,nout1) = cp + dpp
                        zout(1,j,nout5) = ap - bp
                        zout(2,j,nout5) = cp - dpp
                        zout(1,j,nout3) = am + dm
                        zout(2,j,nout3) = cm - bm
                        zout(1,j,nout7) = am - dm
                        zout(2,j,nout7) = cm + bm
                        r=r1 - r5
                        s=s3 - s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r3 - r7
                        bp=r + s
                        bm=r - s
                        r=s4 - s8
                        s=r2 - r6
                        cp=r + s
                        cm=r - s
                        r=s2 - s6
                        s=r4 - r8
                        dpp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( dm - cp)*rt2i
                        cp= ( cm + dpp)*rt2i
                        dpp = ( cm - dpp)*rt2i
                        zout(1,j,nout2) = ap + r
                        zout(2,j,nout2) = bm + s
                        zout(1,j,nout6) = ap - r
                        zout(2,j,nout6) = bm - s
                        zout(1,j,nout4) = am + cp
                        zout(2,j,nout4) = bp + dpp
                        zout(1,j,nout8) = am - cp
                        zout(2,j,nout8) = bp - dpp
                        enddo
                        enddo
8000                continue

        else
                ia=1
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        do j=1,n1dfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r2=zin(1,j,nin2)
                        s2=zin(2,j,nin2)
                        r3=zin(1,j,nin3)
                        s3=zin(2,j,nin3)
                        r4=zin(1,j,nin4)
                        s4=zin(2,j,nin4)
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r6=zin(1,j,nin6)
                        s6=zin(2,j,nin6)
                        r7=zin(1,j,nin7)
                        s7=zin(2,j,nin7)
                        r8=zin(1,j,nin8)
                        s8=zin(2,j,nin8)
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dpp=r + s
                        dm=r - s
                        zout(1,j,nout1) = ap + bp
                        zout(2,j,nout1) = cp + dpp
                        zout(1,j,nout5) = ap - bp
                        zout(2,j,nout5) = cp - dpp
                        zout(1,j,nout3) = am - dm
                        zout(2,j,nout3) = cm + bm
                        zout(1,j,nout7) = am + dm
                        zout(2,j,nout7) = cm - bm
                        r= r1 - r5
                        s=-s3 + s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r7 - r3
                        bp=r + s
                        bm=r - s
                        r=-s4 + s8
                        s= r2 - r6
                        cp=r + s
                        cm=r - s
                        r=-s2 + s6
                        s= r4 - r8
                        dpp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( cp - dm)*rt2i
                        cp= ( cm + dpp)*rt2i
                        dpp= ( dpp - cm)*rt2i
                        zout(1,j,nout2) = ap + r
                        zout(2,j,nout2) = bm + s
                        zout(1,j,nout6) = ap - r
                        zout(2,j,nout6) = bm - s
                        zout(1,j,nout4) = am + cp
                        zout(2,j,nout4) = bp + dpp
                        zout(1,j,nout8) = am - cp
                        zout(2,j,nout8) = bp - dpp
                        enddo
                        enddo

                do 8001,ia=2,after
                ias=ia-1
                        itt=ias*before
                        itrig=itt+1
                        cr2=trig(1,itrig)
                        ci2=trig(2,itrig)
                        itrig=itrig+itt
                        cr3=trig(1,itrig)
                        ci3=trig(2,itrig)
                        itrig=itrig+itt
                        cr4=trig(1,itrig)
                        ci4=trig(2,itrig)
                        itrig=itrig+itt
                        cr5=trig(1,itrig)
                        ci5=trig(2,itrig)
                        itrig=itrig+itt
                        cr6=trig(1,itrig)
                        ci6=trig(2,itrig)
                        itrig=itrig+itt
                        cr7=trig(1,itrig)
                        ci7=trig(2,itrig)
                        itrig=itrig+itt
                        cr8=trig(1,itrig)
                        ci8=trig(2,itrig)
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nin6=nin5+atb
                        nin7=nin6+atb
                        nin8=nin7+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        nout6=nout5+after
                        nout7=nout6+after
                        nout8=nout7+after
                        do j=1,n1dfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=r*cr2 - s*ci2
                        s2=r*ci2 + s*cr2
                        r=zin(1,j,nin3)
                        s=zin(2,j,nin3)
                        r3=r*cr3 - s*ci3
                        s3=r*ci3 + s*cr3
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=r*cr4 - s*ci4
                        s4=r*ci4 + s*cr4
                        r=zin(1,j,nin5)
                        s=zin(2,j,nin5)
                        r5=r*cr5 - s*ci5
                        s5=r*ci5 + s*cr5
                        r=zin(1,j,nin6)
                        s=zin(2,j,nin6)
                        r6=r*cr6 - s*ci6
                        s6=r*ci6 + s*cr6
                        r=zin(1,j,nin7)
                        s=zin(2,j,nin7)
                        r7=r*cr7 - s*ci7
                        s7=r*ci7 + s*cr7
                        r=zin(1,j,nin8)
                        s=zin(2,j,nin8)
                        r8=r*cr8 - s*ci8
                        s8=r*ci8 + s*cr8
                        r=r1 + r5
                        s=r3 + r7
                        ap=r + s
                        am=r - s
                        r=r2 + r6
                        s=r4 + r8
                        bp=r + s
                        bm=r - s
                        r=s1 + s5
                        s=s3 + s7
                        cp=r + s
                        cm=r - s
                        r=s2 + s6
                        s=s4 + s8
                        dpp=r + s
                        dm=r - s
                        zout(1,j,nout1) = ap + bp
                        zout(2,j,nout1) = cp + dpp
                        zout(1,j,nout5) = ap - bp
                        zout(2,j,nout5) = cp - dpp
                        zout(1,j,nout3) = am - dm
                        zout(2,j,nout3) = cm + bm
                        zout(1,j,nout7) = am + dm
                        zout(2,j,nout7) = cm - bm
                        r= r1 - r5
                        s=-s3 + s7
                        ap=r + s
                        am=r - s
                        r=s1 - s5
                        s=r7 - r3
                        bp=r + s
                        bm=r - s
                        r=-s4 + s8
                        s= r2 - r6
                        cp=r + s
                        cm=r - s
                        r=-s2 + s6
                        s= r4 - r8
                        dpp=r + s
                        dm=r - s
                        r = ( cp + dm)*rt2i
                        s = ( cp - dm)*rt2i
                        cp= ( cm + dpp)*rt2i
                        dpp= ( dpp - cm)*rt2i
                        zout(1,j,nout2) = ap + r
                        zout(2,j,nout2) = bm + s
                        zout(1,j,nout6) = ap - r
                        zout(2,j,nout6) = bm - s
                        zout(1,j,nout4) = am + cp
                        zout(2,j,nout4) = bp + dpp
                        zout(1,j,nout8) = am - cp
                        zout(2,j,nout8) = bp - dpp
                        enddo
                        enddo
8001                continue

        end if
        else if (now.eq.3) then
!         .5d0*sqrt(3.d0)
        bb=isign*0.8660254037844387d0
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        do j=1,n1dfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r2=zin(1,j,nin2)
        s2=zin(2,j,nin2)
        r3=zin(1,j,nin3)
        s3=zin(2,j,nin3)
        r=r2 + r3
        s=s2 + s3
        zout(1,j,nout1) = r + r1
        zout(2,j,nout1) = s + s1
        r1=r1 - .5d0*r
        s1=s1 - .5d0*s
        r2=bb*(r2-r3)
        s2=bb*(s2-s3)
        zout(1,j,nout2) = r1 - s2
        zout(2,j,nout2) = s1 + r2
        zout(1,j,nout3) = r1 + s2
        zout(2,j,nout3) = s1 - r2
        enddo
        enddo
        do 3000,ia=2,after
        ias=ia-1
        if (4*ias.eq.3*after) then
        if (isign.eq.1) then
                nin1=ia-after
                nout1=ia-atn
                do ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do j=1,n1dfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(2,j,nin2)
                s2=zin(1,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r=r3 + r2
                s=s2 - s3
                zout(1,j,nout1) = r1 - r
                zout(2,j,nout1) = s + s1
                r1=r1 + .5d0*r
                s1=s1 - .5d0*s
                r2=bb*(r2-r3)
                s2=bb*(s2+s3)
                zout(1,j,nout2) = r1 - s2
                zout(2,j,nout2) = s1 - r2
                zout(1,j,nout3) = r1 + s2
                zout(2,j,nout3) = s1 + r2
                enddo
                enddo
        else
                nin1=ia-after
                nout1=ia-atn
                do ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do j=1,n1dfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r2=zin(2,j,nin2)
                s2=zin(1,j,nin2)
                r3=zin(1,j,nin3)
                s3=zin(2,j,nin3)
                r=r2 - r3
                s=s2 + s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s1 - s
                r1=r1 - .5d0*r
                s1=s1 + .5d0*s
                r2=bb*(r2+r3)
                s2=bb*(s2-s3)
                zout(1,j,nout2) = r1 + s2
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 - s2
                zout(2,j,nout3) = s1 - r2
                enddo
                enddo
        end if
        else if (8*ias.eq.3*after) then
        if (isign.eq.1) then
                nin1=ia-after
                nout1=ia-atn
                do ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do j=1,n1dfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r - s)*rt2i
                s2=(r + s)*rt2i
                r3=zin(2,j,nin3)
                s3=zin(1,j,nin3)
                r=r2 - r3
                s=s2 + s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s + s1
                r1=r1 - .5d0*r
                s1=s1 - .5d0*s
                r2=bb*(r2+r3)
                s2=bb*(s2-s3)
                zout(1,j,nout2) = r1 - s2
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 + s2
                zout(2,j,nout3) = s1 - r2
                enddo
                enddo
        else
                nin1=ia-after
                nout1=ia-atn
                do ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                do j=1,n1dfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=(r + s)*rt2i
                s2=(s - r)*rt2i
                r3=zin(2,j,nin3)
                s3=zin(1,j,nin3)
                r=r2 + r3
                s=s2 - s3
                zout(1,j,nout1) = r + r1
                zout(2,j,nout1) = s + s1
                r1=r1 - .5d0*r
                s1=s1 - .5d0*s
                r2=bb*(r2-r3)
                s2=bb*(s2+s3)
                zout(1,j,nout2) = r1 - s2
                zout(2,j,nout2) = s1 + r2
                zout(1,j,nout3) = r1 + s2
                zout(2,j,nout3) = s1 - r2
                enddo
                enddo
        end if
        else
        itt=ias*before
        itrig=itt+1
        cr2=trig(1,itrig)
        ci2=trig(2,itrig)
        itrig=itrig+itt
        cr3=trig(1,itrig)
        ci3=trig(2,itrig)
        nin1=ia-after
        nout1=ia-atn
        do ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        do j=1,n1dfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r=zin(1,j,nin2)
        s=zin(2,j,nin2)
        r2=r*cr2 - s*ci2
        s2=r*ci2 + s*cr2
        r=zin(1,j,nin3)
        s=zin(2,j,nin3)
        r3=r*cr3 - s*ci3
        s3=r*ci3 + s*cr3
        r=r2 + r3
        s=s2 + s3
        zout(1,j,nout1) = r + r1
        zout(2,j,nout1) = s + s1
        r1=r1 - .5d0*r
        s1=s1 - .5d0*s
        r2=bb*(r2-r3)
        s2=bb*(s2-s3)
        zout(1,j,nout2) = r1 - s2
        zout(2,j,nout2) = s1 + r2
        zout(1,j,nout3) = r1 + s2
        zout(2,j,nout3) = s1 - r2
        enddo
        enddo
        end if
3000        continue
        else if (now==5) then
!         cos(2.d0*pi/5.d0)
        cos2=0.3090169943749474d0
!         cos(4.d0*pi/5.d0)
        cos4=-0.8090169943749474d0
!        sin(2.d0*pi/5.d0)
        sin2=isign*0.9510565162951536d0
!         sin(4.d0*pi/5.d0)
        sin4=isign*0.5877852522924731d0
        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        do j=1,n1dfft
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        r2=zin(1,j,nin2)
        s2=zin(2,j,nin2)
        r3=zin(1,j,nin3)
        s3=zin(2,j,nin3)
        r4=zin(1,j,nin4)
        s4=zin(2,j,nin4)
        r5=zin(1,j,nin5)
        s5=zin(2,j,nin5)
        r25 = r2 + r5
        r34 = r3 + r4
        s25 = s2 - s5
        s34 = s3 - s4
        zout(1,j,nout1) = r1 + r25 + r34
        r = r1 + cos2*r25 + cos4*r34
        s = sin2*s25 + sin4*s34
        zout(1,j,nout2) = r - s
        zout(1,j,nout5) = r + s
        r = r1 + cos4*r25 + cos2*r34
        s = sin4*s25 - sin2*s34
        zout(1,j,nout3) = r - s
        zout(1,j,nout4) = r + s
        r25 = r2 - r5
        r34 = r3 - r4
        s25 = s2 + s5
        s34 = s3 + s4
        zout(2,j,nout1) = s1 + s25 + s34
        r = s1 + cos2*s25 + cos4*s34
        s = sin2*r25 + sin4*r34
        zout(2,j,nout2) = r + s
        zout(2,j,nout5) = r - s
        r = s1 + cos4*s25 + cos2*s34
        s = sin4*r25 - sin2*r34
        zout(2,j,nout3) = r + s
        zout(2,j,nout4) = r - s
        enddo
        enddo
        do 5000,ia=2,after
        ias=ia-1
        if (8*ias.eq.5*after) then
                if (isign.eq.1) then
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        do j=1,n1dfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r - s)*rt2i
                        s2=(r + s)*rt2i
                        r3=zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(r + s)*rt2i
                        s4=(r - s)*rt2i
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r25 = r2 - r5
                        r34 = r3 + r4
                        s25 = s2 + s5
                        s34 = s3 - s4
                        zout(1,j,nout1) = r1 + r25 - r34
                        r = r1 + cos2*r25 - cos4*r34
                        s = sin2*s25 + sin4*s34
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout5) = r + s
                        r = r1 + cos4*r25 - cos2*r34
                        s = sin4*s25 - sin2*s34
                        zout(1,j,nout3) = r - s
                        zout(1,j,nout4) = r + s
                        r25 = r2 + r5
                        r34 = r4 - r3
                        s25 = s2 - s5
                        s34 = s3 + s4
                        zout(2,j,nout1) = s1 + s25 + s34
                        r = s1 + cos2*s25 + cos4*s34
                        s = sin2*r25 + sin4*r34
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout5) = r - s
                        r = s1 + cos4*s25 + cos2*s34
                        s = sin4*r25 - sin2*r34
                        zout(2,j,nout3) = r + s
                        zout(2,j,nout4) = r - s
                        enddo
                        enddo
                else
                        nin1=ia-after
                        nout1=ia-atn
                        do ib=1,before
                        nin1=nin1+after
                        nin2=nin1+atb
                        nin3=nin2+atb
                        nin4=nin3+atb
                        nin5=nin4+atb
                        nout1=nout1+atn
                        nout2=nout1+after
                        nout3=nout2+after
                        nout4=nout3+after
                        nout5=nout4+after
                        do j=1,n1dfft
                        r1=zin(1,j,nin1)
                        s1=zin(2,j,nin1)
                        r=zin(1,j,nin2)
                        s=zin(2,j,nin2)
                        r2=(r + s)*rt2i
                        s2=(s - r)*rt2i
                        r3=zin(2,j,nin3)
                        s3=zin(1,j,nin3)
                        r=zin(1,j,nin4)
                        s=zin(2,j,nin4)
                        r4=(s - r)*rt2i
                        s4=(r + s)*rt2i
                        r5=zin(1,j,nin5)
                        s5=zin(2,j,nin5)
                        r25 = r2 - r5
                        r34 = r3 + r4
                        s25 = s2 + s5
                        s34 = s4 - s3
                        zout(1,j,nout1) = r1 + r25 + r34
                        r = r1 + cos2*r25 + cos4*r34
                        s = sin2*s25 + sin4*s34
                        zout(1,j,nout2) = r - s
                        zout(1,j,nout5) = r + s
                        r = r1 + cos4*r25 + cos2*r34
                        s = sin4*s25 - sin2*s34
                        zout(1,j,nout3) = r - s
                        zout(1,j,nout4) = r + s
                        r25 = r2 + r5
                        r34 = r3 - r4
                        s25 = s2 - s5
                        s34 = s3 + s4
                        zout(2,j,nout1) = s1 + s25 - s34
                        r = s1 + cos2*s25 - cos4*s34
                        s = sin2*r25 + sin4*r34
                        zout(2,j,nout2) = r + s
                        zout(2,j,nout5) = r - s
                        r = s1 + cos4*s25 - cos2*s34
                        s = sin4*r25 - sin2*r34
                        zout(2,j,nout3) = r + s
                        zout(2,j,nout4) = r - s
                        enddo
                        enddo
                end if
        else
                ias=ia-1
                itt=ias*before
                itrig=itt+1
                cr2=trig(1,itrig)
                ci2=trig(2,itrig)
                itrig=itrig+itt
                cr3=trig(1,itrig)
                ci3=trig(2,itrig)
                itrig=itrig+itt
                cr4=trig(1,itrig)
                ci4=trig(2,itrig)
                itrig=itrig+itt
                cr5=trig(1,itrig)
                ci5=trig(2,itrig)
                nin1=ia-after
                nout1=ia-atn
                do ib=1,before
                nin1=nin1+after
                nin2=nin1+atb
                nin3=nin2+atb
                nin4=nin3+atb
                nin5=nin4+atb
                nout1=nout1+atn
                nout2=nout1+after
                nout3=nout2+after
                nout4=nout3+after
                nout5=nout4+after
                do j=1,n1dfft
                r1=zin(1,j,nin1)
                s1=zin(2,j,nin1)
                r=zin(1,j,nin2)
                s=zin(2,j,nin2)
                r2=r*cr2 - s*ci2
                s2=r*ci2 + s*cr2
                r=zin(1,j,nin3)
                s=zin(2,j,nin3)
                r3=r*cr3 - s*ci3
                s3=r*ci3 + s*cr3
                r=zin(1,j,nin4)
                s=zin(2,j,nin4)
                r4=r*cr4 - s*ci4
                s4=r*ci4 + s*cr4
                r=zin(1,j,nin5)
                s=zin(2,j,nin5)
                r5=r*cr5 - s*ci5
                s5=r*ci5 + s*cr5
                r25 = r2 + r5
                r34 = r3 + r4
                s25 = s2 - s5
                s34 = s3 - s4
                zout(1,j,nout1) = r1 + r25 + r34
                r = r1 + cos2*r25 + cos4*r34
                s = sin2*s25 + sin4*s34
                zout(1,j,nout2) = r - s
                zout(1,j,nout5) = r + s
                r = r1 + cos4*r25 + cos2*r34
                s = sin4*s25 - sin2*s34
                zout(1,j,nout3) = r - s
                zout(1,j,nout4) = r + s
                r25 = r2 - r5
                r34 = r3 - r4
                s25 = s2 + s5
                s34 = s3 + s4
                zout(2,j,nout1) = s1 + s25 + s34
                r = s1 + cos2*s25 + cos4*s34
                s = sin2*r25 + sin4*r34
                zout(2,j,nout2) = r + s
                zout(2,j,nout5) = r - s
                r = s1 + cos4*s25 + cos2*s34
                s = sin4*r25 - sin2*r34
                zout(2,j,nout3) = r + s
                zout(2,j,nout4) = r - s
                enddo
                enddo
        end if
5000        continue
       else if (now.eq.6) then
!         .5d0*sqrt(3.d0)
        bb=isign*0.8660254037844387d0

        ia=1
        nin1=ia-after
        nout1=ia-atn
        do ib=1,before
        nin1=nin1+after
        nin2=nin1+atb
        nin3=nin2+atb
        nin4=nin3+atb
        nin5=nin4+atb
        nin6=nin5+atb
        nout1=nout1+atn
        nout2=nout1+after
        nout3=nout2+after
        nout4=nout3+after
        nout5=nout4+after
        nout6=nout5+after
        do j=1,n1dfft
        r2=zin(1,j,nin3)
        s2=zin(2,j,nin3)
        r3=zin(1,j,nin5)
        s3=zin(2,j,nin5)
        r=r2 + r3
        s=s2 + s3
        r1=zin(1,j,nin1)
        s1=zin(2,j,nin1)
        ur1 = r + r1
        ui1 = s + s1
        r1=r1 - .5d0*r
        s1=s1 - .5d0*s
        r=r2-r3
        s=s2-s3
        ur2 = r1 - s*bb
        ui2 = s1 + r*bb
        ur3 = r1 + s*bb
        ui3 = s1 - r*bb

        r2=zin(1,j,nin6)
        s2=zin(2,j,nin6)
        r3=zin(1,j,nin2)
        s3=zin(2,j,nin2)
        r=r2 + r3
        s=s2 + s3
        r1=zin(1,j,nin4)
        s1=zin(2,j,nin4)
        vr1 = r + r1
        vi1 = s + s1
        r1=r1 - .5d0*r
        s1=s1 - .5d0*s
        r=r2-r3
        s=s2-s3
        vr2 = r1 - s*bb
        vi2 = s1 + r*bb
        vr3 = r1 + s*bb
        vi3 = s1 - r*bb

        zout(1,j,nout1)=ur1+vr1
        zout(2,j,nout1)=ui1+vi1
        zout(1,j,nout5)=ur2+vr2
        zout(2,j,nout5)=ui2+vi2
        zout(1,j,nout3)=ur3+vr3
        zout(2,j,nout3)=ui3+vi3
        zout(1,j,nout4)=ur1-vr1
        zout(2,j,nout4)=ui1-vi1
        zout(1,j,nout2)=ur2-vr2
        zout(2,j,nout2)=ui2-vi2
        zout(1,j,nout6)=ur3-vr3
        zout(2,j,nout6)=ui3-vi3
        enddo
        enddo

        else
          ABI_ERROR('error fftstp')
        end if

end subroutine fftstp
!!***

end module m_sg2002
!!***
