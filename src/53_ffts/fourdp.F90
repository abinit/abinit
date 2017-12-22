!{\src2tex{textfont=tt}}
!!****f* ABINIT/fourdp
!! NAME
!! fourdp
!!
!! FUNCTION
!! Conduct Fourier transform of REAL or COMPLEX function f(r)=fofr defined on
!! fft grid in real space, to create complex f(G)=fofg defined on full fft grid 
!! in reciprocal space, in full storage mode, or the reverse operation.
!! For the reverse operation, the final data is divided by nfftot.
!! REAL case when cplex=1, COMPLEX case when cplex=2
!! Usually used for density and potentials.
!!
!! There are two different possibilities :
!!  fftalgb=0 means using the complex-to-complex FFT routine,
!!   irrespective of the value of cplex
!!  fftalgb=1 means using a real-to-complex FFT or a complex-to-complex FFT,
!!   depending on the value of cplex.
!!  The only real-to-complex FFT available is from SGoedecker library.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cplex=1 if fofr is real, 2 if fofr is complex
!! isign=sign of Fourier transform exponent: current convention uses
!!  +1 for transforming from G to r 
!!  -1 for transforming from r to G.
!! mpi_enreg=information about MPI parallelization
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! paral_kgb=Flag related to the kpoint-band-fft parallelism
!! tim_fourdp=timing code of the calling routine (can be set to 0 if not attributed)
!!
!! TODO
!!  Remove paral_kgb
!!
!! SIDE EFFECTS
!! Input/Output
!! fofg(2,nfft)=f(G), complex.
!! fofr(cplex*nfft)=input function f(r) (real or complex)
!!
!! PARENTS
!!      atm2fft,bethe_salpeter,calc_smeared_density,dfpt_atm2fft,dfpt_dyfro
!!      dfpt_eltfrxc,dfpt_looppert,dfpt_newvtr,dfpt_scfcv,dfpt_vlocal
!!      dfptnl_loop,dieltcel,energy,fock_getghc,forces,fourdp_6d,fresidrsp
!!      green_kernel,gstate,hartre,hartrestr,initro,jellium,laplacian,m_dvdb
!!      m_electronpositron,m_epjdos,m_fft_prof,m_hidecudarec,m_kxc,m_ppmodel
!!      m_screening,make_efg_el,mklocl_realspace,mklocl_recipspace,moddiel
!!      moddiel_csrb,mrgscr,newrho,newvtr,nonlinear,nres2vres,odamix,pawmknhat
!!      pawmknhat_psipsi,pawmkrho,posdoppler,prcref,prcref_PMA,recursion
!!      recursion_nl,redgr,respfn,rotate_rho,scfcv,screening,setup_positron
!!      sigma,stress,symrhg,tddft,transgrid,vlocalstr,xcden,xcpot
!!
!! CHILDREN
!!      ccfft,dfti_seqfourdp,fftw3_mpifourdp,fftw3_seqfourdp,fourdp_mpi
!!      ptabs_fourdp,sg2002_back,sg2002_forw,sg2002_mpifourdp,sg_fft_rc,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fourdp(cplex,fofg,fofr,isign,mpi_enreg,nfft,ngfft,paral_kgb,tim_fourdp)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_fftcore
 
 use m_mpinfo,      only : ptabs_fourdp
 use m_sgfft,       only : sg_fft_rc
 use m_sg2002,      only : sg2002_mpifourdp, sg2002_back, sg2002_forw
 use m_fftw3,       only : fftw3_seqfourdp, fftw3_mpifourdp
 use m_dfti,        only : dfti_seqfourdp
 use m_fft,         only : fourdp_mpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fourdp'
 use interfaces_18_timing
 use interfaces_53_ffts, except_this_one => fourdp
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,isign,nfft,paral_kgb,tim_fourdp
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(inout) :: fofg(2,nfft),fofr(cplex*nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1=1
 integer :: fftalg,fftalga,fftalgb,fftcache,i1,i2,i3,base
 integer :: n1,n1half1,n1halfm,n2,n2half1,n3,n4
 integer :: n4half1,n5,n5half1,n6 !nd2proc,nd3proc,i3_local,i2_local,
 integer :: comm_fft,nproc_fft,me_fft
 real(dp) :: xnorm
 character(len=500) :: message
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: work1(:,:,:,:),work2(:,:,:,:)
 real(dp),allocatable :: workf(:,:,:,:),workr(:,:,:,:)

! *************************************************************************

 ABI_UNUSED(paral_kgb)
 
 ! Keep track of timing
 call timab(260+tim_fourdp,1,tsec)

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6)
 me_fft=ngfft(11); nproc_fft=ngfft(10)
 comm_fft = mpi_enreg%comm_fft
 !write(std_out,*)"fourdp, nx,ny,nz,nfft =",n1,n2,n3,nfft

 fftcache=ngfft(8)
 fftalg  =ngfft(7); fftalga =fftalg/100; fftalgb =mod(fftalg,100)/10

 xnorm=one/dble(n1*n2*n3)
 !write(std_out,*)' fourdp :me_fft',me_fft,'nproc_fft',nproc_fft,'nfft',nfft

 if (fftalgb/=0 .and. fftalgb/=1) then
   write(message, '(a,i4,a,a,a,a,a)' )&
&   'The input algorithm number fftalg=',fftalg,' is not allowed.',ch10,&
&   'The second digit (fftalg(B)) must be 0 or 1.',ch10,&
&   'Action : change fftalg in your input file.'
   MSG_BUG(message)
 end if

 if (fftalgb==1 .and. ALL(fftalga/=(/1,3,4,5/)) )then
   write(message,'(a,i4,5a)')&
&   'The input algorithm number fftalg=',fftalg,' is not allowed.',ch10,&
&   'When fftalg(B) is 1, the allowed values for fftalg(A) are 1 and 4.',ch10,&
&   'Action: change fftalg in your input file.'
   MSG_BUG(message)
 end if

 if (n4<n1.or.n5<n2.or.n6<n3) then
   write(message,'(a,3i8,a,3i8)')'  Each of n4,n5,n6=',n4,n5,n6,'must be >= n1, n2, n3 =',n1,n2,n3
   MSG_BUG(message)
 end if

 ! Get the distrib associated with this fft_grid => for i2 and i3 planes
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 ! Branch immediately depending on nproc_fft 
 if (nproc_fft > 1) then
   call fourdp_mpi(cplex,nfft,ngfft,ndat1,isign,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)
   goto 100
 end if

 if (fftalga == FFT_FFTW3) then
   ! Call sequential or MPI FFTW3 version.
   if (nproc_fft == 1) then
     !call wrtout(std_out,"FFTW3 SEQFOURDP","COLL")
     call fftw3_seqfourdp(cplex,n1,n2,n3,n1,n2,n3,ndat1,isign,fofg,fofr)
   else
     !call wrtout(std_out,"FFTW3 MPIFOURDP","COLL")
     call fftw3_mpifourdp(cplex,nfft,ngfft,ndat1,isign,&
&     fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)
   end if
   ! Accumulate timing and return
   call timab(260+tim_fourdp,2,tsec); return
 end if

 if (fftalga==FFT_DFTI) then 
   ! Call sequential or MPI MKL.
   if (nproc_fft == 1) then
     call dfti_seqfourdp(cplex,n1,n2,n3,n1,n2,n3,ndat1,isign,fofg,fofr)
   else
     MSG_ERROR("MPI fourdp with MKL cluster DFT not implemented")
   end if
   ! Accumulate timing and return
   call timab(260+tim_fourdp,2,tsec); return
 end if

 ! Here, deal  with the new SG FFT, complex-to-complex case
 if (fftalga==FFT_SG2002 .and. (fftalgb==0 .or. cplex==2)) then
   call sg2002_mpifourdp(cplex,nfft,ngfft,ndat1,isign,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)
   !call sg2002_seqfourdp(cplex,nfft,ngfft,ndat1,isign,fftn2_fofg,fofr)
 end if

 ! Here, deal with the new SG FFT, with real-to-complex
 if (fftalga==FFT_SG2002 .and. fftalgb==1 .and. cplex==1) then
   ABI_CHECK(nproc_fft == 1,"fftalg 41x does not support nproc_fft > 1")

   n1half1=n1/2+1; n1halfm=(n1+1)/2
   n2half1=n2/2+1
   ! n4half1 or n5half1 are the odd integers >= n1half1 or n2half1
   n4half1=(n1half1/2)*2+1
   n5half1=(n2half1/2)*2+1
   ABI_ALLOCATE(workr,(2,n4half1,n5,n6))
   ABI_ALLOCATE(workf,(2,n4,n6,n5half1))

   if (isign==1) then
     do i3=1,n3
       do i2=1,n2half1
         base=n1*(i2-1+n2*(i3-1))
         do i1=1,n1
           workf(1,i1,i3,i2)=fofg(1,i1+base)
           workf(2,i1,i3,i2)=fofg(2,i1+base)
         end do
       end do
     end do

     !nd2proc=((n2-1)/nproc_fft) +1
     !nd3proc=((n6-1)/nproc_fft) +1

     ! change the call? n5half1 et n6 ?
     call sg2002_back(cplex,ndat1,n1,n2,n3,n4,n5,n6,n4half1,n5half1,n6,2,workf,workr,comm_fft)

     do i3=1,n3
       do i2=1,n2
         base=n1*(i2-1+n2*(i3-1))
         do i1=1,n1half1-1
           ! copy data
           fofr(2*i1-1+base)=workr(1,i1,i2,i3)
           fofr(2*i1  +base)=workr(2,i1,i2,i3)
         end do
         ! If n1 odd, must add last data
         if((2*n1half1-2)/=n1)then
           fofr(n1+base)=workr(1,n1half1,i2,i3)
         end if
       end do
     end do

   else if (isign==-1) then
     do i3=1,n3
       do i2=1,n2
         base=n1*(i2-1+n2*(i3-1))
         do i1=1,n1half1-1
           workr(1,i1,i2,i3)=fofr(2*i1-1+base)
           workr(2,i1,i2,i3)=fofr(2*i1  +base)
         end do
         ! If n1 odd, must add last data
         if((2*n1half1-2)/=n1)then
           workr(1,n1half1,i2,i3)=fofr(n1+base)
           workr(2,n1half1,i2,i3)=zero
         end if
       end do
     end do

     call sg2002_forw(cplex,ndat1,n1,n2,n3,n4,n5,n6,n4half1,n5half1,n6,2,workr,workf,comm_fft)

     ! Transfer fft output to the original fft box
     do i3=1,n3
       do i2=1,n2half1
         base=n1*(i2-1+n2*(i3-1))
         do i1=1,n1
           fofg(1,i1+base)=workf(1,i1,i3,i2)*xnorm
           fofg(2,i1+base)=workf(2,i1,i3,i2)*xnorm
         end do
       end do

       ! Complete missing values with complex conjugate
       ! Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.
       if(n2half1>2)then
         do i2=2,n2+1-n2half1
           base=n1*((n2+2-i2)-1)
           if(i3/=1)base=base+n1*n2*((n3+2-i3)-1)
           fofg(1,1+base)= workf(1,1,i3,i2)*xnorm
           fofg(2,1+base)=-workf(2,1,i3,i2)*xnorm
           do i1=2,n1
             fofg(1,n1+2-i1+base)= workf(1,i1,i3,i2)*xnorm
             fofg(2,n1+2-i1+base)=-workf(2,i1,i3,i2)*xnorm
           end do
         end do
       end if
     end do

   end if ! isign
   ABI_DEALLOCATE(workr)
   ABI_DEALLOCATE(workf)
 end if

 ! Here, one calls the complex-to-complex FFT subroutine
 if( (fftalgb==0 .or. cplex==2) .and. fftalga/=4 )then

   ABI_ALLOCATE(work1,(2,n4,n5,n6))
   ABI_ALLOCATE(work2,(2,n4,n5,n6))

   if (isign==1) then

     ! Transfer fofg to the expanded fft box
!$OMP PARALLEL DO PRIVATE(base)
     do i3=1,n3
       do i2=1,n2
         base=n1*(i2-1+n2*(i3-1))
         do i1=1,n1
           work1(1,i1,i2,i3)=fofg(1,i1+base)
           work1(2,i1,i2,i3)=fofg(2,i1+base)
         end do
       end do
     end do

     ! Call Stefan Goedecker C2C FFT
     !call sg_fft_cc(fftcache,n1,n2,n3,n4,n5,n6,ndat1,isign,work1,work2)
     call ccfft(ngfft,isign,n1,n2,n3,n4,n5,n6,ndat1,2,work1,work2,comm_fft)

     ! Take data from expanded box and put it in the original box.
     if (cplex==1) then
       ! REAL case
!$OMP PARALLEL DO PRIVATE(base)
       do i3=1,n3
         do i2=1,n2
           base=n1*(i2-1+n2*(i3-1))
           do i1=1,n1
             fofr(i1+base)=work2(1,i1,i2,i3)
           end do
         end do
       end do

     else
       ! COMPLEX case
!$OMP PARALLEL DO PRIVATE(base)
       do i3=1,n3
         do i2=1,n2
           base=2*n1*(i2-1+n2*(i3-1))
           do i1=1,n1
             fofr(2*i1-1+base)=work2(1,i1,i2,i3)
             fofr(2*i1  +base)=work2(2,i1,i2,i3)
           end do
         end do
       end do
     end if

   else if (isign==-1) then

     ! Insert fofr into the augmented fft box
     if (cplex==1) then 
       ! REAL case
!$OMP PARALLEL DO PRIVATE(base) 
       do i3=1,n3
         do i2=1,n2
           base=n1*(i2-1+n2*(i3-1))
           do i1=1,n1
             ! copy data
             work1(1,i1,i2,i3)=fofr(i1+base)
             work1(2,i1,i2,i3)=zero
           end do
         end do
       end do
     else
       ! COMPLEX case
!$OMP PARALLEL DO PRIVATE(base)
       do i3=1,n3
         do i2=1,n2
           base=2*n1*(i2-1+n2*(i3-1))
           do i1=1,n1
             ! copy data
             work1(1,i1,i2,i3)=fofr(2*i1-1+base)
             work1(2,i1,i2,i3)=fofr(2*i1  +base)
           end do
         end do
       end do
     end if ! cplex

     ! Call Stefan Goedecker C2C FFT
     !call sg_fft_cc(fftcache,n1,n2,n3,n4,n5,n6,ndat1,isign,work1,work2)
     call ccfft(ngfft,isign,n1,n2,n3,n4,n5,n6,ndat1,2,work1,work2,comm_fft)

     ! Transfer fft output to the original fft box
!$OMP PARALLEL DO PRIVATE(base)
     do i3=1,n3
       do i2=1,n2
         base=n1*(i2-1+n2*(i3-1))
         do i1=1,n1
           fofg(1,i1+base)=work2(1,i1,i2,i3)*xnorm
           fofg(2,i1+base)=work2(2,i1,i2,i3)*xnorm
         end do
       end do
     end do

   end if ! isign

   ABI_DEALLOCATE(work1)
   ABI_DEALLOCATE(work2)
 end if ! End simple algorithm

 ! Here sophisticated algorithm based on S. Goedecker routines, only for the REAL case.
 ! Take advantage of the fact that fofr is real, and that fofg has corresponding symmetry properties.
 if( (fftalgb==1 .and. cplex==1) .and. fftalga/=4 )then
   ABI_CHECK(nproc_fft==1,"nproc > 1 not supported")
   ABI_CHECK(ndat1==1,"ndat > 1 not supported")
   call sg_fft_rc(cplex,fofg,fofr,isign,nfft,ngfft)
 end if

 100 call timab(260+tim_fourdp,2,tsec)

end subroutine fourdp
!!***
