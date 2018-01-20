!{\src2tex{textfont=tt}}
!!****f* ABINIT/ccfft
!! NAME
!! ccfft
!!
!! FUNCTION
!! Carry out complex-to-complex Fourier transforms between real
!! and reciprocal (G) space. Library of such routines.
!! Include machine-dependent F90 routines used with fftalg=200.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2017 ABINIT group (PT, XG, FF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  fftalga=govern the choice of the fft routine to be used
!!    if 1: SGoedecker routine
!!    if 2: Machine dependent routine, depending on the precompilation options
!!    if 3: FFTW library routine
!!    if 4: new SGoedecker routine, version 2002
!!          Warning : the second and third dimensions of the Fourier space
!!          array are switched, compared to the usual case
!!  fftcache=size of the cache (kB)
!!  isign= Integer specifying which sign to be used for the transformation.
!!         must be either +1 or -1.
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  n1,n2,n3=Actual integer dimensions (see ngfft) for the 3D sequence.
!!           Physical dimension of the transform.
!!  n4,n5,n6=Leading dimensions. Generally, n6 is not different to n3.
!!  ndat=number of FFT to do in //
!!  option= 1 if call from fourwf, 2 if call from other routine
!!  work1(2,n4*n5*n6)=Array to be transformed.
!!
!! OUTPUT
!!  inplace = 0 if result is in work2 ; =1 if result is in work1 (machine dependent)
!!  normalized=0 if the backward (isign=-1) FFT is not normalized, so has to be normalized outside of ccfft
!!            =1 otherwise
!!  work2(2,n4*n5*n6)=transformed array in case inplace=0.
!!
!! SIDE EFFECTS
!!  work1(2,n4*n5*n6)=at input, array to be transformed
!!                    at output, transformed array (in case inplace=1)
!!
!! NOTES
!! precompilation definitions :
!!   -D(machine_list) :  (case fftalga=200)
!!      choice of machine-dependent FFT library, if permitted
!!   -DHAVE_FFT_FFTW2   : (case fftalga=300) activate the FFTW lib
!!   -Dnolib  : (case fftalga=200) call SGoedecker routine,
!!      instead of machine-dependent one
!!
!! More about fftalga=200
!! Library routines for the following platforms have been implemented :
!!  Compaq/DEC
!!  HP          (in place FFT)
!!  SGI         (in place FFT)
!!  NEC         (in place FFT)
!! For all the other platforms, or if the CPP directive nolib is
!! activated, one uses the fft routine from S. Goedecker.
!!
!! PARENTS
!!      fourdp,fourwf
!!
!! CHILDREN
!!      sg2002_back,sg2002_forw,sg_fft_cc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ccfft(ngfft,isign,n1,n2,n3,n4,n5,n6,ndat,option,work1,work2,comm_fft)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_sgfft,     only : sg_fft_cc
 use m_sg2002,    only : sg2002_back, sg2002_forw

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ccfft'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isign,n1,n2,n3,n4,n5,n6,ndat,option,comm_fft
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(inout) :: work1(2,n4*n5*n6*ndat)
 real(dp),intent(inout) :: work2(2,n4*n5*n6*ndat) !vz_i

!Local variables ------------------------------
!scalars
 integer,parameter :: cplex2=2
 integer :: fftalg,fftalga,fftalgb,fftalgc,fftcache
 integer :: nd2proc,nd3proc,nproc_fft
 character(len=500) :: message

!*************************************************************************

 nproc_fft=ngfft(10)

 fftcache=ngfft(8)
 fftalg  =ngfft(7)
 fftalga =fftalg/100; fftalgb=mod(fftalg,100)/10; fftalgc=mod(fftalg,10)

 if(fftalga==2)then
   MSG_ERROR("Machine dependent FFTs are not supported anymore")

 else if(fftalga==3)then
   MSG_ERROR("Old interface with FFTW2 is not supported anymore")

 else if(fftalga<1 .or. fftalga>4)then
   write(message, '(a,a,a,i5,a,a)' )&
&   'The allowed values of fftalg(A) are 1, 2, 3, and 4 .',ch10,&
&   'The actual value of fftalg(A) is',fftalga,ch10,&
&   'Action : check the value of fftalg in your input file.'
   MSG_ERROR(message)
 end if

 ! This routine will be removed ASAP.
 ! Do not add new FFT libraries without previous discussion with Matteo Giantomassi
 ! inplace==1 or normalize==2 are not supported anymore in the caller (fourwf, fourdp)
 !inplace=0; normalized=0

 if (fftalga/=4) then
   ! Call Stefan Goedecker FFT
   call sg_fft_cc(fftcache,n1,n2,n3,n4,n5,n6,ndat,isign,work1,work2)

 else if (fftalga==4) then
   ! Call new version of Stefan Goedecker FFT
   nd2proc=((n2-1)/nproc_fft) +1
   nd3proc=((n6-1)/nproc_fft) +1

   if (isign==1) then 
     ! Fourier to real space (backward)
     call sg2002_back(cplex2,ndat,n1,n2,n3,n4,n5,n6,n4,nd2proc,nd3proc,option,work1,work2,comm_fft)
   else 
     ! isign=-1, real space to Fourier (forward)
     call sg2002_forw(cplex2,ndat,n1,n2,n3,n4,n5,n6,n4,nd2proc,nd3proc,option,work1,work2,comm_fft)
   end if
 end if

end subroutine ccfft
!!***
