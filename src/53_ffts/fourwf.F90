!{\src2tex{textfont=tt}}
!!****f* ABINIT/fourwf
!! NAME
!! fourwf
!!
!! FUNCTION
!! Carry out composite Fourier transforms between real and reciprocal (G) space.
!! Wavefunctions, contained in a sphere in reciprocal space,
!! can be FFT to real space. They can also be FFT from real space
!! to a sphere. Also, the density maybe accumulated, and a local
!! potential can be applied.
!!
!! The different options are :
!! - option=0 --> reciprocal to real space and output the result.
!! - option=1 --> reciprocal to real space and accumulate the density.
!! - option=2 --> reciprocal to real space, apply the local potential to the wavefunction
!!                in real space and produce the result in reciprocal space.
!! - option=3 --> real space to reciprocal space.
!!                NOTE that in this case, fftalg=1x1 MUST be used. This may be changed in the future.
!!
!! The different sections of this routine corresponds to different
!! algorithms, used independently of each others :
!!(read first the description of the fftalg input variable in abinit_help)
!! - fftalg=xx0 : use simple complex-to-complex routines, without zero padding
!!     (rather simple, so can be used to understand how fourwf.f works);
!! - fftalg=1x1 : use S Goedecker routines, with zero padding
!!     (7/12 savings in execution time);
!! - fftalg=1x2 : call even more sophisticated coding also based on S Goedecker routines
!!
!! This routine contains many parts that differ only
!! by small details, in order to treat each case with the better speed.
!! Also for better speed, it uses no F90 construct, except the allocate command
!! and for zeroing arrays.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, FF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cplex= if 1 , denpot is real, if 2 , denpot is complex
!!    (cplex=2 only allowed for option=2, and istwf_k=1)
!!    not relevant if option=0 or option=3, so cplex=0 can be used to minimize memory
!! fofgin(2,npwin)=holds input wavefunction in G vector basis sphere.
!!                 (intent(in) but the routine sphere can modify it for another iflag)
!! gboundin(2*mgfft+8,2)=sphere boundary info for reciprocal to real space
!! gboundout(2*mgfft+8,2)=sphere boundary info for real to reciprocal space
!! istwf_k=option parameter that describes the storage of wfs
!! kg_kin(3,npwin)=reduced planewave coordinates, input
!! kg_kout(3,npwout)=reduced planewave coordinates, output
!! mgfft=maximum size of 1D FFTs
!! mpi_enreg=information about MPI parallelization
!! ndat=number of FFT to do in //
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! npwin=number of elements in fofgin array (for option 0, 1 and 2)
!! npwout=number of elements in fofgout array (for option 2 and 3)
!! n4,n5,n6=ngfft(4),ngfft(5),ngfft(6), dimensions of fofr.
!! option= if 0: do direct FFT
!!         if 1: do direct FFT, then sum the density
!!         if 2: do direct FFT, multiply by the potential, then do reverse FFT
!!         if 3: do reverse FFT only
!! paral_kgb=Flag related to the kpoint-band-fft parallelism
!! tim_fourwf=timing code of the calling routine (can be set to 0 if not attributed)
!! weight_r=weight to be used for the accumulation of the density in real space
!!         (needed only when option=1)
!! weight_i=weight to be used for the accumulation of the density in real space
!!         (needed only when option=1 and (fftalg=4 and fftalgc/=0))
!! fofginb(2,npwin)=holds second input wavefunction in G vector basis sphere.
!!                 (intent(in) but the routine sphere can modify it for another iflag)
!!                 (for non diagonal occupation)
!! use_ndo = use non diagonal occupations.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! for option==0, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                fofr(2,n4,n5,n6*ndat) contains the output Fourier Transform of fofgin;
!!                no use of denpot, fofgout and npwout.
!! for option==1, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                denpot(cplex*n4,n5,n6) contains the input density at input,
!!                and the updated density at output (accumulated);
!!                no use of fofgout and npwout.
!! for option==2, fofgin(2,npwin*ndat)=holds input wavefunction in G sphere;
!!                denpot(cplex*n4,n5,n6) contains the input local potential;
!!                fofgout(2,npwout*ndat) contains the output function;
!! for option==3, fofr(2,n4,n5,n6*ndat) contains the input real space wavefunction;
!!                fofgout(2,npwout*ndat) contains its output Fourier transform;
!!                no use of fofgin and npwin.
!!
!! TODO
!!  Remove paral_kgb, we are already passing mpi_enreg
!!
!! NOTES
!!   DO NOT CHANGE THE API OF THIS FUNCTION.
!!   If you need a specialized routine for the FFT of the wavefunctions, create
!!   a wrapper that uses fourwf to accomplish your task. This routine, indeed,
!!   has already too many parameters and each change in the API requires a careful
!!   modification of the different wrappers used for specialized FFTs such as FFTW3 and MKL-DFTI
!!
!! PARENTS
!!      dfpt_accrho,dfpt_mkrho,dfptnl_resp,fock_getghc,getgh1c,getghc
!!      gwls_hamiltonian,m_cut3d,m_epjdos,m_fft_prof,m_fock,mkrho,mlwfovlp
!!      pawmkaewf,pawsushat,posdoppler,prep_fourwf,spin_current,susk,suskmm
!!      tddft,vtowfk
!!
!! CHILDREN
!!      ccfft,cg_addtorho,cg_box2gsph,dcopy,dfti_seqfourwf,fftw3_seqfourwf
!!      fourwf_mpi,gpu_fourwf,ptabs_fourwf,sg_fftpad,sg_fftrisc,sg_fftrisc_2
!!      sphere,sphere_fft,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
&  kg_kin,kg_kout,mgfft,mpi_enreg,ndat,ngfft,npwin,npwout,n4,n5,n6,option,&
&  paral_kgb,tim_fourwf,weight_r,weight_i, &
&  use_gpu_cuda,use_ndo,fofginb) ! Optional arguments

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_cgtools

 use m_mpinfo,    only : ptabs_fourwf
 use m_fftcore,   only : sphere_fft, sphere
 use m_sgfft,     only : sg_fftpad, sg_fftrisc, sg_fftrisc_2
 use m_dfti,      only : dfti_seqfourwf
 use m_fftw3,     only : fftw3_seqfourwf
 use m_fft,       only : fourwf_mpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fourwf'
 use interfaces_18_timing
 use interfaces_53_ffts, except_this_one => fourwf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,mgfft,n4,n5,n6,ndat,npwin,npwout,option,paral_kgb
 integer,intent(in) :: tim_fourwf
 integer,intent(in),optional :: use_gpu_cuda,use_ndo
 real(dp),intent(in) :: weight_r,weight_i
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout),ngfft(18)
 real(dp),intent(inout) :: denpot(cplex*n4,n5,n6),fofgin(2,npwin*ndat)
 real(dp),intent(inout),optional :: fofginb(:,:) ! (2,npwin*ndat)
 real(dp),intent(inout) :: fofr(2,n4,n5,n6*ndat)
 real(dp),intent(out) :: fofgout(2,npwout*ndat)

!Local variables-------------------------------
!scalars
 integer :: fftalg,fftalga,fftalgc,fftcache,i1,i2,i2_local,i3,i3_local,i3_glob,idat,ier
 integer :: iflag,ig,comm_fft,me_g0,me_fft,n1,n2,n3,nd2proc,nd3proc
 integer :: nfftot,nproc_fft,option_ccfft 
 real(dp) :: fim,fre,xnorm
 character(len=500) :: message
 logical ::  luse_gpu_cuda,luse_ndo 
!arrays
 integer,parameter :: shiftg0(3)=0
 integer,parameter :: symmE(3,3)=reshape([1,0,0,0,1,0,0,0,1],[3,3])
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: work1(:,:,:,:),work2(:,:,:,:),work3(:,:,:,:)
 real(dp),allocatable :: work4(:,:,:,:),work_sum(:,:,:,:) 

! *************************************************************************

 ! Accumulate timing
 call timab(840+tim_fourwf,1,tsec)

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3); nfftot=n1*n2*n3
 fftcache=ngfft(8)
 fftalg=ngfft(7); fftalga=fftalg/100; fftalgc=mod(fftalg,10)
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)

 comm_fft = mpi_enreg%comm_fft; me_g0 = mpi_enreg%me_g0

 !if (ndat/=1) then
 !  write(std_out,*)fftalg
 !  MSG_ERROR("Really? I thought nobody uses ndat > 1")
 !end if 

 !if (weight_r /= weight_i) then
 !  write(std_out,*)fftalg
 !  MSG_ERROR("Really? I thought nobody uses weight_r != weight_i")
 !end if 

 !if (option == 0 .and. fftalgc == 0) then
 !  MSG_ERROR("Option 0 is buggy when fftalgc ==0 is used!")
 !end if 

!Cuda version of fourwf
 luse_gpu_cuda=PRESENT(use_gpu_cuda)
 if (luse_gpu_cuda) luse_gpu_cuda=(luse_gpu_cuda.and.(use_gpu_cuda==1))

 if(luse_gpu_cuda) then
#if defined HAVE_GPU_CUDA
   call gpu_fourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
&   kg_kin,kg_kout,mgfft,mpi_enreg,ndat,ngfft,npwin,npwout,n4,n5,n6,option,&
&   paral_kgb,tim_fourwf,weight_r,weight_i) !,&
!  &  use_ndo,fofginb)
#endif
   call timab(840+tim_fourwf,2,tsec); return
 end if

 if ((fftalgc<0 .or. fftalgc>2)) then
   write(message, '(a,i4,a,a,a,a,a)' )&
&   'The input algorithm number fftalg=',fftalg,' is not allowed.',ch10,&
&   'The third digit, fftalg(C), must be 0, 1, or 2',ch10,&
&   'Action: change fftalg in your input file.'
   MSG_ERROR(message)
 end if

 if (fftalgc/=0 .and. ALL(fftalga/=(/1,3,4,5/)) ) then
   write(message, '(a,i4,5a)' )&
&   'The input algorithm number fftalg=',fftalg,' is not allowed.',ch10,&
&   'The first digit must be 1,3,4 when the last digit is not 0.',ch10,&
&   'Action : change fftalg in your input file.'
   MSG_ERROR(message)
 end if

 if (option<0 .or. option>3)then
   write(message, '(a,i4,a,a,a)' )&
&   'The option number',option,' is not allowed.',ch10,&
&   'Only option=0, 1, 2 or 3 are allowed presently.'
   MSG_ERROR(message)
 end if

 if (option==1 .and. cplex/=1) then
   write(message, '(a,a,a,i4,a)' )&
&   'With the option number 1, cplex must be 1,',ch10,&
&   'but it is cplex=',cplex,'.'
   MSG_ERROR(message)
 end if

 if (option==2 .and. (cplex/=1 .and. cplex/=2)) then
   write(message, '(a,a,a,i4,a)' )&
&   'With the option number 2, cplex must be 1 or 2,',ch10,&
&   'but it is cplex=',cplex,'.'
   MSG_ERROR(message)
 end if

 ! DMFT uses its own FFT algorithm (that should be wrapped in a different routine!)
 luse_ndo=.false.
 if (present(use_ndo).and.present(fofginb)) then
   if(use_ndo==1) then 
     luse_ndo=.true.
     if((size(fofginb,2)==0)) then
       write(message, '(a,a,a,i4,i5)' )&
&       'fofginb has a dimension equal to zero and use_ndo==1',ch10,&
&       'Action : check dimension of fofginb',size(fofginb,2),use_ndo
       MSG_ERROR(message)
     end if
   end if
 end if

 if (luse_ndo) then 
   if (.not.(fftalgc==2 .and. option/=3)) then
     MSG_ERROR("luse_ndo but not .not.(fftalgc==2 .and. option/=3)")
   end if
   ABI_CHECK(nproc_fft==1, "DMFT with nproc_fft != 1")
   ABI_CHECK(ndat == 1, "use_ndo and ndat != 1 not coded")

   call sg_fftrisc_2(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,&
&   istwf_k,kg_kin,kg_kout,&
&   mgfft,ngfft,npwin,npwout,n4,n5,n6,option,weight_r,weight_2=weight_i,luse_ndo=luse_ndo,fofgin_p=fofginb)
   goto 100
 end if 

 ! Get the distrib associated with this fft_grid => for i2 and i3 planes
 call ptabs_fourwf(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 ! Branch immediately depending on nproc_fft 
 if (nproc_fft > 1 .and. fftalg /= 412) then
   call fourwf_mpi(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,&
&   istwf_k,kg_kin,kg_kout,me_g0,mgfft,ngfft,mpi_enreg%distribfft,n1,n2,n3,npwin,npwout,&
&   n4,n5,n6,ndat,option,weight_r,weight_i,comm_fft)
   goto 100
 end if

 select case (fftalga)

 case (FFT_FFTW3)
   if (luse_ndo) MSG_ERROR("luse_ndo not supported by FFTW3")
   if (nproc_fft == 1) then
!      call wrtout(std_out,"FFTW3_SEQFOURWF","COLL")
     call fftw3_seqfourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
&     kg_kin,kg_kout,mgfft,ndat,ngfft,npwin,npwout,n4,n5,n6,option,weight_r,weight_i)
   else
     MSG_ERROR("Not coded")
   end if

 case (FFT_DFTI) 
   if (luse_ndo) MSG_ERROR("luse_ndo not supported by DFTI")
   if (nproc_fft == 1) then
!     call wrtout(std_out,"DFTI_SEQFOURWF","COLL")
     call dfti_seqfourwf(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,istwf_k,&
&     kg_kin,kg_kout,mgfft,ndat,ngfft,npwin,npwout,n4,n5,n6,option,weight_r,weight_i)
   else
     MSG_ERROR("Not coded")
   end if

 case default 
   ! TODO: Clean this code!

   ! Here, use routines that make forwards FFT separately of backwards FFT,
   ! in particular, usual 3DFFT library routines, called in ccfft.
   if (fftalgc==0 .or. (fftalgc==1 .and. fftalga/=4) .or. &
&   (fftalgc==2 .and. fftalga/=4 .and. option==3) )then

     ABI_ALLOCATE(work1,(2,n4,n5,n6*ndat))

     if (option/=3)then 
       ! Insert fofgin into the fft box (array fofr)

       if (fftalga/=4)then
         iflag=1
         call sphere(fofgin,ndat,npwin,fofr,n1,n2,n3,n4,n5,n6,kg_kin,istwf_k,iflag,me_g0,shiftg0,symmE,one)

       else if (fftalga==4 .and. fftalgc==0) then
         ! Note the switch of n5 and n6, as they are only
         ! needed to dimension work2 inside "sphere"
         ABI_ALLOCATE(work2,(2,n4,n6,n5*ndat))

         iflag=2
         nd2proc=((n2-1)/nproc_fft) +1
         nd3proc=((n6-1)/nproc_fft) +1
         ABI_ALLOCATE(work3,(2,n4,n6,nd2proc*ndat))
         ABI_ALLOCATE(work4,(2,n4,n5,nd3proc*ndat))

         if (istwf_k == 1 .and. paral_kgb==1) then
           ! sphere dont need a big array
           work3=zero
           call sphere_fft(fofgin,ndat,npwin,work3,n1,n2,n3,n4,n6,kg_kin,&
&           mpi_enreg%distribfft%tab_fftwf2_local,nd2proc)
         else
           ! sphere needs a big array and communications
           if (nproc_fft == 1 .and. ndat == 1 .and. istwf_k == 1) then
             ! dimensions of tab work3 and work2 are identical no need to use work2
             work3=zero
             call sphere(fofgin,ndat,npwin,work3,n1,n2,n3,n4,n6,nd2proc,&
&             kg_kin,istwf_k,iflag,me_g0,shiftg0,symmE,one)
           else
             work2=zero
             call sphere(fofgin,ndat,npwin,work2,n1,n2,n3,n4,n6,n5,&
&             kg_kin,istwf_k,iflag,me_g0,shiftg0,symmE,one)

             if (paral_kgb==1 .and. istwf_k > 1) then
               ! Collect G-vectors on each node
               work3=zero
               ABI_ALLOCATE(work_sum,(2,n4,n6,n5*ndat))
               call timab(48,1,tsec)
               call xmpi_sum(work2,work_sum,2*n4*n6*n5*ndat,comm_fft,ier)
               call timab(48,2,tsec)

               ! Extract my list of G-vectors needed for MPI-FFT.
               do idat=1,ndat
                 do i2=1,n2
                   if( fftn2_distrib(i2) == me_fft) then
                     i2_local = ffti2_local(i2) + nd2proc*(idat-1)
                     do i3=1,n3
                       do i1=1,n1
                         work3(1,i1,i3,i2_local)=work_sum(1,i1,i3,i2+n5*(idat-1))
                         work3(2,i1,i3,i2_local)=work_sum(2,i1,i3,i2+n5*(idat-1))
                       end do
                     end do
                   end if
                 end do
               end do
               ABI_DEALLOCATE(work_sum)
             end if

             if (paral_kgb/=1) then
               do idat=1,ndat
                 do i2=1,n2
                   do i3=1,n3
                     do i1=1,n1
                       work3(1,i1,i3,i2+nd2proc*(idat-1))=work2(1,i1,i3,i2+n5*(idat-1))
                       work3(2,i1,i3,i2+nd2proc*(idat-1))=work2(2,i1,i3,i2+n5*(idat-1))
                     end do
                   end do
                 end do
               end do
             end if
           end if
         end if
         if (paral_kgb==1) then
           option_ccfft=1
         else
           option_ccfft=2
         end if
       end if

       ! Fourier transform fofr (reciprocal to real space)
       ! The output might be in work1 or fofr, depending on inplace
       if (fftalgc==0) then
         if (fftalga/=4) then 
           ! Call usual 3DFFT library routines 
           call ccfft(ngfft,+1,n1,n2,n3,n4,n5,n6,ndat,2,fofr,work1,comm_fft)
         else 
           ! SG simplest complex-to-complex routine
           call ccfft(ngfft,+1,n1,n2,n3,n4,n5,n6,ndat,option_ccfft,work3,work4,comm_fft)
           ABI_DEALLOCATE(work2)
           ABI_DEALLOCATE(work3)
         end if
       else
         ! Call SG routine, with zero padding
         call sg_fftpad(fftcache,mgfft,n1,n2,n3,n4,n5,n6,ndat,gboundin,+1,fofr,work1)
       end if
     end if ! option/=3

     ! Note that if option==0 everything is alright already, the output is available in fofr.
     ! MG: TODO: Rewrite this mess in human-readable form!
     if (option==0) then 
       if (fftalgc==0) then
         if (fftalga/=4) then 
           call DCOPY(2*n4*n5*n6*ndat,work1,1,fofr,1)
         else 
           call DCOPY(2*n4*n5*n6*ndat,work4,1,fofr,1)
         end if
       else
         ! Results are copied to fofr.
         call DCOPY(2*n4*n5*n6*ndat,work1,1,fofr,1)
       end if
     end if 

     if (option==1) then 
       ! Accumulate density
       if ((fftalgc==0) .and. (fftalga==4)) then
         do idat=1,ndat
           do i3=1,n3
             if( me_fft == fftn3_distrib(i3) ) then
               i3_local = ffti3_local(i3) + nd3proc*(idat-1)
               do i2=1,n2
                 do i1=1,n1
                   denpot(i1,i2,i3)=denpot(i1,i2,i3)+&
&                   weight_r*work4(1,i1,i2,i3_local)**2+&
&                   weight_i*work4(2,i1,i2,i3_local)**2
                 end do
               end do
             end if
           end do
         end do ! idat
       else
         call cg_addtorho(n1,n2,n3,n4,n5,n6,ndat,weight_r,weight_i,work1,denpot)
       end if
     end if ! option==1

     if (option==2) then 
       ! Apply local potential
       if (cplex==1) then

         if ((fftalgc==0) .and. (fftalga==4)) then
!$OMP PARALLEL DO PRIVATE(i3_local,i3_glob) 
           do idat=1,ndat
             do i3=1,n3
               if( me_fft == fftn3_distrib(i3) ) then
                 i3_local = ffti3_local(i3) + nd3proc*(idat-1)
                 i3_glob = i3+n3*(idat-1)
                 do i2=1,n2
                   do i1=1,n1
                     fofr(1,i1,i2,i3_glob)= denpot(i1,i2,i3)*work4(1,i1,i2,i3_local)
                     fofr(2,i1,i2,i3_glob)= denpot(i1,i2,i3)*work4(2,i1,i2,i3_local)
                   end do
                 end do
               end if
             end do
           end do
         end if
         if ((fftalgc/=0) .or. (fftalga/=4)) then
!$OMP PARALLEL DO PRIVATE(i3_glob)
           do idat=1,ndat
             do i3=1,n3
               if( me_fft == fftn3_distrib(i3) ) then
                 i3_glob = i3+n3*(idat-1)
                 do i2=1,n2
                   do i1=1,n1
                     fofr(1,i1,i2,i3_glob)=denpot(i1,i2,i3)*work1(1,i1,i2,i3+n3*(idat-1))
                     fofr(2,i1,i2,i3_glob)=denpot(i1,i2,i3)*work1(2,i1,i2,i3+n3*(idat-1))
                   end do
                 end do
               end if
             end do
           end do
         end if

       else if (cplex==2) then
         if ((fftalgc==0) .and. (fftalga==4)) then
!$OMP PARALLEL DO PRIVATE(fre,fim,i3_local,i3_glob) 
           do idat=1,ndat
             do i3=1,n3
               if( me_fft == fftn3_distrib(i3) ) then
                 i3_local = ffti3_local(i3) + nd3proc*(idat-1)
                 i3_glob = i3+n3*(idat-1)
                 do i2=1,n2
                   do i1=1,n1
                     fre=work4(1,i1,i2,i3_local)
                     fim=work4(2,i1,i2,i3_local)
                     fofr(1,i1,i2,i3_glob)=denpot(2*i1-1,i2,i3)*fre -denpot(2*i1,i2,i3)*fim
                     fofr(2,i1,i2,i3_glob)=denpot(2*i1-1,i2,i3)*fim +denpot(2*i1,i2,i3)*fre
                   end do
                 end do
               end if
             end do
           end do
         end if

         if ((fftalgc/=0) .or. (fftalga/=4)) then
!$OMP PARALLEL DO PRIVATE(fre,fim,i3_glob)
           do idat=1,ndat
             do i3=1,n3
               if( me_fft == fftn3_distrib(i3) ) then
                 i3_glob = i3+n3*(idat-1)
                 do i2=1,n2
                   do i1=1,n1
                     fre=work1(1,i1,i2,i3+n3*(idat-1))
                     fim=work1(2,i1,i2,i3+n3*(idat-1))
                     fofr(1,i1,i2,i3_glob)=denpot(2*i1-1,i2,i3)*fre -denpot(2*i1,i2,i3)*fim
                     fofr(2,i1,i2,i3_glob)=denpot(2*i1-1,i2,i3)*fim +denpot(2*i1,i2,i3)*fre
                   end do
                 end do
               end if
             end do
           end do
         end if
       end if ! cplex=2

     end if ! option=2

     ! The data for option==2 or option==3 is now in fofr.
     if (option==2 .or. option==3) then

       if (fftalgc==0) then 
         ! Call usual 3DFFT library routines or SG simplest complex-to-complex routine
         if (fftalga==FFT_SG2002) then
           ABI_DEALLOCATE(work1)
           ABI_ALLOCATE(work1,(2,n4,n6,n5*ndat))
         end if

         if (option==3 .or. fftalga/=4) then
           call ccfft(ngfft,-1,n1,n2,n3,n4,n5,n6,ndat,2,fofr,work1,comm_fft)
         else
           ! creation of small arrays
           ! nd3proc=((n5-1)/nproc_fft) +1
           nd3proc=((n6-1)/nproc_fft) +1
           nd2proc=((n2-1)/nproc_fft) +1
           ABI_ALLOCATE(work3,(2,n4,n5,nd3proc*ndat))
           ABI_ALLOCATE(work2,(2,n4,n6,nd2proc*ndat))

           if (paral_kgb==1) then

             if (cplex==1) then
               do idat=1,ndat
                 do i3=1,n3
                   if( me_fft == fftn3_distrib(i3) ) then
                     i3_local = ffti3_local(i3) + nd3proc*(idat-1)
                     do i2=1,n2
                       do i1=1,n1
                         work3(1,i1,i2,i3_local)=denpot(i1,i2,i3)*work4(1,i1,i2,i3_local)
                         work3(2,i1,i2,i3_local)=denpot(i1,i2,i3)*work4(2,i1,i2,i3_local)
                       end do
                     end do
                   end if
                 end do
               end do
             else
               do idat=1,ndat
                 do i3=1,n3
                   if( me_fft == fftn3_distrib(i3) ) then
                     i3_local = ffti3_local(i3) + nd3proc*(idat-1)
                     do i2=1,n2
                       do i1=1,n1
                         fre=work4(1,i1,i2,i3_local)
                         fim=work4(2,i1,i2,i3_local)
                         work3(1,i1,i2,i3_local) = denpot(2*i1-1,i2,i3)*fre-denpot(2*i1,i2,i3)*fim
                         work3(2,i1,i2,i3_local) = denpot(2*i1-1,i2,i3)*fim+denpot(2*i1,i2,i3)*fre
                       end do
                     end do
                   end if
                 end do
               end do
             end if
             option_ccfft=1

           else 
             if (nproc_fft /=1 .or. ndat /= 1 ) then
               do idat=1,ndat
                 do i3=1,n3
                   do i2=1,n2
                     do i1=1,n1
                       work3(1,i1,i2,i3+nd3proc*(idat-1))=fofr(1,i1,i2,i3+n3*(idat-1))
                       work3(2,i1,i2,i3+nd3proc*(idat-1))=fofr(2,i1,i2,i3+n3*(idat-1))
                     end do
                   end do
                 end do
               end do
               option_ccfft=2
             end if
           end if

           if (paral_kgb==1) then
             call ccfft(ngfft,-1,n1,n2,n3,n4,n5,n6,ndat,option_ccfft,work3,work2,comm_fft)
           else
             if (nproc_fft /=1 .or. ndat /= 1 ) then
               call ccfft(ngfft,-1,n1,n2,n3,n4,n5,n6,ndat,option_ccfft,work3,work2,comm_fft)
             else
               call ccfft(ngfft,-1,n1,n2,n3,n4,n5,n6,ndat,option_ccfft,fofr,work1,comm_fft)
             end if
           end if

           ! load of work1
           if ((paral_kgb==1) .and.  ( istwf_k > 1 )) work1(:,:,:,:)=zero

           if (paral_kgb==1) then
             if ( istwf_k > 1 ) then
               do idat=1,ndat
                 do i2=1,n2
                   if( me_fft == fftn2_distrib(i2) ) then
                     i2_local = ffti2_local(i2) + nd2proc*(idat-1)
                     do i3=1,n3
                       do i1=1,n1
                         work1(1,i1,i3,i2+n5*(idat-1))= work2(1,i1,i3,i2_local)
                         work1(2,i1,i3,i2+n5*(idat-1))= work2(2,i1,i3,i2_local)
                       end do
                     end do
                   end if
                 end do
               end do
             end if

           else
             if (nproc_fft /=1 .or. ndat /= 1 ) then
               do idat=1,ndat
                 do i2=1,n2
                   do i3=1,n3
                     do i1=1,n1
                       work1(1,i1,i3,i2+n5*(idat-1))=work2(1,i1,i3,i2+nd2proc*(idat-1))
                       work1(2,i1,i3,i2+n5*(idat-1))=work2(2,i1,i3,i2+nd2proc*(idat-1))
                     end do
                   end do
                 end do
               end do
             end if
           end if
           ABI_DEALLOCATE(work3)
           if ((paral_kgb==1) .and.  ( istwf_k > 1 )) then
             call timab(48,1,tsec)
             call xmpi_sum(work1,comm_fft,ier)
             call timab(48,2,tsec)
           end if
         end if

       else
         ! Call SG routine, with zero padding
         call sg_fftpad(fftcache,mgfft,n1,n2,n3,n4,n5,n6,ndat,gboundout,-1,fofr,work1)
       end if

       xnorm = one/dble(nfftot)

       if (fftalga/=4) then
         call cg_box2gsph(n1,n2,n3,n4,n5,n6,ndat,npwout,kg_kout,work1,fofgout, rscal=xnorm)
       else 
         ! if fftalga==4
         if ((paral_kgb==1) .and. ( istwf_k == 1 )) then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,i2_local)
           do idat=1,ndat
             do ig=1,npwout
               i1=kg_kout(1,ig); if(i1<0)i1=i1+n1 ; i1=i1+1
               i2=kg_kout(2,ig); if(i2<0)i2=i2+n2 ; i2=i2+1
               i3=kg_kout(3,ig); if(i3<0)i3=i3+n3 ; i3=i3+1
               i2_local = ffti2_local(i2) + nd2proc*(idat-1)
               fofgout(1,ig+npwout*(idat-1))= work2(1,i1,i3,i2_local)*xnorm
               fofgout(2,ig+npwout*(idat-1))= work2(2,i1,i3,i2_local)*xnorm
             end do
           end do
           ABI_DEALLOCATE(work2)
         else
!$OMP PARALLEL DO PRIVATE(i1,i2,i3) 
           do idat=1,ndat
             do ig=1,npwout
               i1=kg_kout(1,ig); if(i1<0)i1=i1+n1; i1=i1+1
               i2=kg_kout(2,ig); if(i2<0)i2=i2+n2; i2=i2+1
               i3=kg_kout(3,ig); if(i3<0)i3=i3+n3; i3=i3+1
               fofgout(1,ig+npwout*(idat-1))=work1(1,i1,i3,i2+n5*(idat-1))*xnorm
               fofgout(2,ig+npwout*(idat-1))=work1(2,i1,i3,i2+n5*(idat-1))*xnorm
             end do
           end do
         end if
       end if ! fftalga
     end if ! if option==2 or 3

     if (allocated(work1))  then
       ABI_DEALLOCATE(work1)
     end if
   end if

   ! Here, call more sophisticated specialized 3-dimensional fft
   ! (zero padding as well as maximize cache reuse) based on S Goedecker routines.
   ! Specially tuned for cache architectures.
   if (fftalga==FFT_SG .and. fftalgc==2 .and. option/=3) then
     call sg_fftrisc(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,&
&     istwf_k,kg_kin,kg_kout,mgfft,ndat,ngfft,npwin,npwout,n4,n5,n6,option,weight_r,weight_i)
   end if

   ! Here, call new FFT from S Goedecker, also sophisticated specialized 3-dimensional fft
   ! (zero padding as well as maximize cache reuse)
   if (fftalga==FFT_SG2002 .and. fftalgc/=0) then
     ! The args are not the same as fourwf, but might be
     call fourwf_mpi(cplex,denpot,fofgin,fofgout,fofr,gboundin,gboundout,&
&     istwf_k,kg_kin,kg_kout,me_g0,mgfft,ngfft,mpi_enreg%distribfft,n1,n2,n3,npwin,npwout,&
&     n4,n5,n6,ndat,option,weight_r,weight_i,comm_fft)
   end if

   if (allocated(work4))  then
     ABI_DEALLOCATE(work4)
   end if
   if (allocated(work2))  then
     ABI_DEALLOCATE(work2)
   end if

 end select

! Accumulate timing
 100 continue
 call timab(840+tim_fourwf,2,tsec)

end subroutine fourwf
!!***
