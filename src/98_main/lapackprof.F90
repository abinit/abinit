!!****p* ABINIT/lapackprof
!! NAME
!! lapackprof
!!
!! FUNCTION
!!  Utility for profiling Linear Algebra libraries used by Abinit.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  (main program)
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,abinit_doctor,cg_zaxpy,cg_zcopy,cg_zgemm
!!      cg_zgemv,cwtime,destroy_mpi_enreg,herald,init_mpi_enreg,lower,projbd
!!      pw_orthon,random_number,sqmat_itranspose,test_xginv,wrtout,xgerc,xheevx
!!      xhpev,xmpi_bcast,xmpi_end,xmpi_init,xomp_set_num_threads,xomp_show_info
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program lapackprof

 use defs_basis
 use m_build_info
 use m_abicore
 use m_xmpi
 use m_xomp
 use m_errors
 use m_hide_blas
 use m_cgtools
 use m_hide_lapack

 use defs_abitypes,   only : MPI_type
 use m_fstrings,      only : lower, itoa, strcat, sjoin
 use m_specialmsg,    only : specialmsg_getcount, herald
 use m_argparse,      only : get_arg, get_arg_list
 use m_time,          only : cwtime
 use m_io_tools,      only : prompt
 use m_numeric_tools, only : arth
 use m_mpinfo,        only : init_mpi_enreg, destroy_mpi_enreg

 implicit none

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0, nspinor1 = 1
 integer :: comm, npw, ierr, my_rank, ii, isz, jj, it, step, icall, nfound
 integer :: istwfk = 1, useoverlap, mcg, mgsc, band, g0, idx, ortalgo, abimem_level, prtvol
 real(dp) ::  ctime, wtime, gflops, abimem_limit_mb
 logical :: do_check
 character(len=24) :: skinds
 character(len=500) :: header, routname, command, arg, msg
 type(MPI_type) :: MPI_enreg
!arrays
 integer,allocatable :: sizes(:)
 real(dp) :: alpha(2), beta(2) ,dot(2)
 real(dp),allocatable :: cg(:,:), gsc(:,:), ortho_check(:,:,:)
 real(dp),allocatable :: cg1(:,:), cg2(:,:), cg3(:,:), ene(:), direc(:,:), scprod(:,:)
 complex(dpc),allocatable :: zvec(:), zmat(:,:), wmat(:,:), zpmat(:), evec(:,:)
! complex(spc),allocatable :: vec(:), mat(:,:)
 type(latime_t) :: Tres
 integer :: ncalls = 50, nband = 200, nsizes = 20, nthreads = 1
 logical :: debug_mode = .FALSE.
 integer :: size_arth(2) = [1000, 2000]

! *************************************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 call xmpi_init()
 comm  = xmpi_world; my_rank = xmpi_comm_rank(comm)

 ! Initialize memory profiling if it is activated
 ! if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
 ! note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
 ABI_CHECK(get_arg("abimem-level", abimem_level, msg, default=0) == 0, msg)
 ABI_CHECK(get_arg("abimem-limit-mb", abimem_limit_mb, msg, default=20.0_dp) == 0, msg)
#ifdef HAVE_MEM_PROFILING
 call abimem_init(abimem_level, limit_mb=abimem_limit_mb)
#endif

 call herald("LAPACKPROF", abinit_version, std_out)
 ABI_CHECK(get_arg("prtvol", prtvol, msg, default=0) == 0, msg)

 call get_command_argument(1, command)

 ! Command line options.
 do ii=2,command_argument_count()
   call get_command_argument(ii, arg)
   if (arg == "-v" .or. arg == "--version") then
     write(std_out,"(a)") trim(abinit_version); goto 100

   else if (arg == "-h" .or. arg == "--help") then
     ! Document the options.
     write(std_out,*)"-v, --version              Show version number and exit."
     write(std_out,*)"-h, --help                 Show this help and exit."
     write(std_out,*)" "
     write(std_out,*)"=== Options for developers ==="
     write(std_out,*)" "
     write(std_out,*)"test_v1complete FILE [--symv1scf 1] [--potfile foo.nc]"
     goto 100
   end if
 end do

 call init_mpi_enreg(mpi_enreg)

 if (my_rank == master) then
   write(std_out,'(a)')" Tool for profiling and testing Linear Algebra routines used in ABINIT."
 end if

 call xomp_set_num_threads(nthreads)
 call xomp_show_info(std_out)

 call xmpi_bcast(size_arth, master, comm, ierr)
 call xmpi_bcast(nsizes, master, comm, ierr)

 ABI_MALLOC(sizes, (nsizes))
 sizes = arth(size_arth(1), size_arth(2), nsizes)

 select case (command)
 case ("projbd")
   istwfk = 1
   useoverlap = 0
   routname = "projbd"
   header = strcat("BEGIN_BENCHMARK: ",routname)
   write(std_out,'(a)')TRIM(header)

   do isz=1,nsizes
     npw = sizes(isz)
     mcg  = npw * nband
     mgsc = mcg * useoverlap
     ABI_MALLOC(cg, (2, mcg))
     call random_number(cg)

     ABI_MALLOC(gsc, (2, mgsc))
     ABI_MALLOC(direc, (2, npw*nspinor1))
     direc = zero
     ABI_MALLOC(scprod, (2,nband))

     call cwtime(ctime,wtime,gflops,"start")

     call projbd(cg,direc,0,0,0,istwfk,mcg,mgsc,nband,npw,nspinor1,gsc,scprod,0,0,useoverlap,&
                 mpi_enreg%me_g0,mpi_enreg%comm_fft)

     call cwtime(ctime,wtime,gflops,"stop")
     write(std_out,'(2(a,i0),2(a,f9.6))')" npw = ",npw,", nband = ",nband,", cpu_time = ",ctime,", wall_time = ",wtime

     ABI_FREE(scprod)
     ABI_FREE(direc)
     ABI_FREE(cg)
     ABI_FREE(gsc)
   end do

   header = strcat("END_BENCHMARK: ",routname)
   write(std_out,'(a)')TRIM(header)

 case ("ortho")
   istwfk = 1
   useoverlap = 0

   do ortalgo=1,4,1
     routname = strcat(" pw_orthon, ortalgo = ",itoa(ortalgo))
     header = strcat("BEGIN_BENCHMARK: ",routname)
     write(std_out,'(a)')TRIM(header)

     do isz=1,nsizes
       npw = sizes(isz)
       mcg  = npw * nband
       mgsc = mcg * useoverlap
       ABI_MALLOC(cg,(2,mcg))
       call random_number(cg)
       ABI_MALLOC(gsc,(2,mgsc))

       if (istwfk/=1) then
         do band=1,nband
           g0 = 1 + (band-1)*npw
           cg(2,g0) = zero
         end do
       end if

       call cwtime(ctime, wtime, gflops, "start")

       call pw_orthon(0,0,istwfk,mcg,mgsc,npw,nband,ortalgo,gsc,useoverlap,cg,&
                      mpi_enreg%me_g0,mpi_enreg%comm_bandspinorfft)

       call cwtime(ctime,wtime,gflops,"stop")

       write(std_out,'(2(a,i0),2(a,f9.6))')" npw = ",npw,", nband = ",nband,", cpu_time = ",ctime,", wall_time = ",wtime

       if (debug_mode) then
         ABI_MALLOC(ortho_check,(2, nband, nband))
         if (istwfk/=1) then
           do band=1,nband
             g0 = 1 + (band-1)*npw
             cg(:,g0) = half * cg(:,g0)
           end do
         end if

         call cg_zgemm("C","N",npw,nband,nband,cg,cg,ortho_check)

         if (istwfk/=1) ortho_check = two * ortho_check

         do band=1,nband
           ortho_check(1,band,band) = ortho_check(1,band,band) - one
         end do

         write(std_out,*)"MAX ABS error",MAXVAL( ABS(RESHAPE(ortho_check,(/2*nband*nband/)) ))
         ABI_FREE(ortho_check)
       end if

       ABI_FREE(cg)
       ABI_FREE(gsc)
     end do

     header = strcat("END_BENCHMARK: ",routname)
     write(std_out,'(a)')TRIM(header)
   end do ! ortalgo

 case ("copy")
   do step=1,2
     if (step == 1) routname = " f90_zcopy"
     if (step == 2) routname = " cg_zcopy"
     header = strcat("BEGIN_BENCHMARK: ",routname)
     write(std_out,'(a)')TRIM(header)

     do isz=1,nsizes
       npw  = sizes(isz)
       ABI_MALLOC(cg1, (2, npw))
       ABI_MALLOC(cg2, (2, npw))
       cg1 = zero

       call cwtime(ctime, wtime, gflops, "start")
       if (step==1) then
         do ii=1,ncalls
           do jj=1,npw
             cg2(:,jj) = cg1(:,jj)
           end do
         end do
       else
         do ii=1,ncalls
           call cg_zcopy(npw, cg1, cg2)
         end do
       end if
       call cwtime(ctime, wtime, gflops, "stop")
       write(std_out,'(a,i0,2(a,f9.6))')" npw = ",npw,", cpu_time = ",ctime,", wall_time = ",wtime

       ABI_FREE(cg1)
       ABI_FREE(cg2)
     end do

     header = strcat("END_BENCHMARK: ",routname)
     write(std_out,'(a)')TRIM(header)
   end do

 case ("zdotc")
   do step=1,2
     if (step == 1) routname = " f90_zdotc"
     if (step == 2) routname = " cg_zdotc"
     header = strcat("BEGIN_BENCHMARK: ",routname)
     write(std_out,'(a)')TRIM(header)

     do isz=1,nsizes
       npw  = sizes(isz)
       ABI_MALLOC(cg1,(2, npw))
       ABI_MALLOC(cg2,(2, npw))
       call random_number(cg1)
       call random_number(cg2)

       call cwtime(ctime, wtime, gflops, "start")
       if (step == 1) then
         do jj=1,ncalls
           dot = zero
!$OMP PARALLEL DO REDUCTION(+:dot)
           do ii=1,npw
             dot(1) = dot(1) + cg1(1,ii)*cg2(1,ii) + cg1(2,ii)*cg2(2,ii)
             dot(2) = dot(2) + cg1(1,ii)*cg2(2,ii) - cg1(2,ii)*cg2(1,ii)
           end do
         end do
       else
         do ii=1,ncalls
           dot = cg_zdotc(npw,cg1,cg2)
         end do
       end if
       call cwtime(ctime, wtime, gflops, "stop")
       write(std_out,'(a,i0,2(a,f9.6))')" npw = ",npw,", cpu_time = ",ctime,", wall_time = ",wtime
       ABI_FREE(cg1)
       ABI_FREE(cg2)
     end do

     header = strcat("END_BENCHMARK: ",routname)
     write(std_out,'(a)')TRIM(header)
   end do

 case ("axpy")
   ! alpha = (/one,zero/)
   alpha = [one, two]
   do step=1,2
     if (step == 1) routname = " f90_axpy"
     if (step == 2) routname = " cg_axpy"
     header = strcat("BEGIN_BENCHMARK: ",routname)
     write(std_out,'(a)')TRIM(header)

     do isz=1,nsizes
       npw = sizes(isz)
       ABI_MALLOC(cg1, (2, npw))
       ABI_MALLOC(cg2, (2, npw))
       cg2 = zero

       do jj=1,npw
         cg1(:,jj) = jj
       end do

       call cwtime(ctime, wtime, gflops, "start")
       if (step == 1) then
         jj = 0
         do icall=1,ncalls
           !call random_number(cg1)
           call random_number(cg2(:,1:1))
           !cg2 = zero
           do ii=1,npw
             jj = jj+1
             cg2(1,ii) = alpha(1)*cg1(1,ii) - alpha(2)*cg1(2,ii) + cg2(1,ii)
             cg2(2,ii) = alpha(1)*cg1(2,ii) + alpha(2)*cg1(1,ii) + cg2(2,ii)
           end do
         end do
       else
         do icall=1,ncalls
           !call random_number(cg1)
           !call random_number(cg2)
           !cg2 = zero
           call random_number(cg2(:,1:1))
           call cg_zaxpy(npw, alpha, cg1, cg2)
         end do
       end if
       call cwtime(ctime, wtime, gflops, "stop")
       write(std_out,'(a,i0,2(a,f9.6))')" npw = ",npw,", cpu_time = ",ctime,", wall_time = ",wtime
       ABI_FREE(cg1)
       ABI_FREE(cg2)
     end do

     header = strcat("END_BENCHMARK: ",routname)
     write(std_out,'(a)')TRIM(header)
   end do

 case ("zgemv")
   alpha = [one, two]
   beta = [zero, zero]
   do step=1,2
     if (step == 1) routname = " f90_zgemv"
     if (step == 2) routname = " cg_zgemv"
     header = strcat("BEGIN_BENCHMARK: ",routname)
     write(std_out,'(a)')TRIM(header)

     do isz=1,nsizes
       npw = sizes(isz)
       ABI_MALLOC(cg1, (2, npw*npw))
       ABI_MALLOC(cg2, (2, npw))
       ABI_MALLOC(cg3, (2, npw))

       do jj=1,npw*npw
         cg1(:,jj) = jj
       end do
       cg2 = one
       cg3 = zero

       call cwtime(ctime, wtime, gflops, "start")
       if (step == 1) then
         !do icall=1,ncalls
         !do jj=1,npw
         !ar=scprod(1,iband);ai=scprod(2,iband)
         !do ipw=1,npw_sp
         !cg_re=cg(1,index1+ipw)
         !cg_im=cg(2,index1+ipw)
         !direc(1,ipw)=direc(1,ipw)-ar*cg_re+ai*cg_im
         !direc(2,ipw)=direc(2,ipw)-ar*cg_im-ai*cg_re
         !end do
         !end do
         !end do
         !!cg3(1,ii) = alpha(1)*cg1(1,ii) - alpha(2)*cg1(2,ii) + cg2(1,ii)
         !!cg3(2,ii) = alpha(1)*cg1(2,ii) + alpha(2)*cg1(1,ii) + cg2(2,ii)
         !end do
         !end do
       else
         do icall=1,ncalls
           call cg_zgemv("N", npw, npw, cg1, cg2, cg3)
         end do
       end if
       call cwtime(ctime, wtime, gflops, "stop")
       write(std_out,'(a,i0,2(a,f9.6))')" npw = ",npw,", cpu_time = ",ctime,", wall_time = ",wtime

       ABI_FREE(cg1)
       ABI_FREE(cg2)
       ABI_FREE(cg3)
     end do

     header = strcat("END_BENCHMARK: ",routname)
     write(std_out,'(a)')TRIM(header)
   end do

 case ("transpose")

   do step=1,2
     if (step == 1) routname = " f90_transpose"
     if (step == 2) routname = " mkl_transponse"
     header = strcat("BEGIN_BENCHMARK: ",routname)
     write(std_out,'(a)')TRIM(header)

     do isz=1,nsizes
       npw = sizes(isz)
       ABI_MALLOC(zmat, (npw, npw))
       ABI_MALLOC(wmat, (npw, npw))
       do jj=1,npw
         do ii=1,npw
           zmat(ii,jj) = DCMPLX(ii,jj)
         end do
       end do

       call cwtime(ctime, wtime, gflops, "start")
       if (step==1) then
         do icall=1,ncalls
           !wmat = TRANSPOSE(zmat)
           zmat = TRANSPOSE(zmat)
         end do
       else
         do icall=1,ncalls
           !call sqmat_otranspose(npw,zmat,wmat)
           call sqmat_itranspose(npw, zmat)
         end do
       end if
       call cwtime(ctime, wtime, gflops, "stop")
       write(std_out,'(a,i0,2(a,f9.6))')"npw = ",npw,", cpu_time ",ctime,", wall_time ",wtime

       ABI_FREE(zmat)
       ABI_FREE(wmat)
     end do

     header = strcat("END_BENCHMARK: ",routname)
     write(std_out,'(a)')TRIM(header)
   end do

  case ("zgerc")

   do step=1,2
     if (step == 1) routname = " f90_zgerc"
     if (step == 2) routname = " zgerc"
     header = strcat("BEGIN_BENCHMARK: ",routname)
     write(std_out,'(a)')TRIM(header)

     do isz=1,nsizes
       npw  = sizes(isz)
       ABI_MALLOC(zvec, (npw))
       ABI_MALLOC(zmat, (npw, npw))
       zvec = cone; zmat = czero

       call cwtime(ctime, wtime, gflops, "start")
       if (step==1) then ! Home made zgerc
         do it=1,ncalls
           zmat = czero
!$OMP PARALLEL DO
           do jj=1,npw
             do ii=1,npw
               zmat(ii,jj) = zmat(ii,jj) + CONJG(zvec(ii)) * zvec(jj)
             end do
           end do
         end do
       else
         do jj=1,ncalls
           zmat = czero
           call XGERC(npw,npw,(1._dp,0._dp),zvec,1,zvec,1,zmat,npw)
         end do
       end if
       call cwtime(ctime, wtime, gflops, "stop")
       write(std_out,'(a,i0,2(a,f9.6))')" npw = ",npw,", cpu_time ",ctime,", wall_time ",wtime

       ABI_FREE(zvec)
       ABI_FREE(zmat)
     end do

     header = strcat("END_BENCHMARK: ",routname)
     write(std_out,'(a)')TRIM(header)
   end do

   !do isz=1,nsizes
   !npw  = sizes(isz)
   !ABI_MALLOC(vec,(npw))
   !ABI_MALLOC(mat,(npw,npw))
   !vec = cone; mat = czero

   !call cwtime(ctime,wtime,gflops,"start")

   !do jj=1,ncalls
   !call XGERC(npw,npw,(1._sp,0._sp),vec,1,vec,1,mat,npw)
   !end do

   !call cwtime(ctime,wtime,gflops,"stop")

   !write(std_out,'(a,i0,2f9.3)')" CGERG size, cpu_time, wall_time, max_abserr ",npw,ctime,wtime

   !ABI_FREE(vec)
   !ABI_FREE(mat)
   !end do

 case ("xginv")
   do_check = .FALSE.
   !do_check = .TRUE.
   do ii=1,nsizes
     npw  = sizes(ii)
     call test_xginv(npw, skinds, do_check, Tres, comm)

     if (my_rank==master) then
       write(std_out,'(a,i0,3f9.3)')&
        " routine = xginv, size, cpu_time, wall_time, max_abserr ",Tres%msize,Tres%ctime,Tres%wtime,Tres%max_abserr
     end if
   end do

 case ("xhpev_vs_xheev")

   do step=1,2
     do isz=1,nsizes
       npw = sizes(isz)
       ABI_REMALLOC(ene, (npw))
       ABI_REMALLOC(evec, (npw, npw))
       ABI_REMALLOC(zpmat, (npw*(npw+1)/2))
       zpmat = czero
       idx = 0
       do jj=1,npw
         do ii=1,jj
           idx = idx + 1
           zpmat(idx) = cone
         end do
       end do

       if (step == 1) then
         call cwtime(ctime,wtime,gflops,"start")
         call xhpev("V", "U", npw, zpmat, ene, evec, npw)
       else
         ABI_MALLOC(zmat, (npw, npw))
         do jj=1,npw
           do ii=1,jj
             idx = ii + jj*(jj-1)/2
             zmat(ii,jj) = zpmat(idx)
           end do
         end do
         call cwtime(ctime, wtime, gflops, "start")
         !call xheev("V","U",npw,zmat,ene)
         call xheevx("V","A","U",npw,zmat,zero,zero,1,1,zero,nfound,ene,evec,npw)
         ABI_FREE(zmat)
       end if
       call cwtime(ctime, wtime, gflops, "stop")
       write(std_out,'(a,i0,2(a,f9.6))')"EIG PROBLEM: size = ",npw,", cpu_time ",ctime,", wall_time ",wtime
     end do

     ABI_FREE(ene)
     ABI_FREE(evec)
     ABI_FREE(zpmat)
   end do

 case default
   MSG_ERROR(sjoin("Wrong command:", command))
 end select

 ABI_FREE(sizes)
 call destroy_mpi_enreg(MPI_enreg)

 call wrtout(std_out,ch10//" Analysis completed.")
 call abinit_doctor("__lapackprof")

100 call xmpi_end()

 end program lapackprof
!!***
