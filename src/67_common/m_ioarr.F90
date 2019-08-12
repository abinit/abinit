!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ioarr
!! NAME
!! m_ioarr
!!
!! FUNCTION
!!  This module provides routines to read/write arrays given on the FFT mesh (densities, potentials ...).
!!  The code supports both Fortran files as well as netcdf files in a transparent way.
!!  The appropriate IO layer is selected from files extensions: netcdf primitives are used if the
!!  file ends with `.nc`. If all the other cases we read/write files in Fortran format.
!!  MPI-IO primitives are used when the FFT arrays are MPI distributed.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, MVer, MT, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_ioarr

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_wffile
 use m_errors
 use m_nctk
 use m_crystal
 use m_ebands
 use m_hdr
 use m_pawrhoij
#ifdef HAVE_MPI2
 use mpi
#endif
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_abitypes,   only : mpi_type, dataset_type
 use defs_datatypes,  only : ebands_t
 use defs_wvltypes,   only : wvl_denspot_type
 use m_time,          only : cwtime, cwtime_report
 use m_io_tools,      only : iomode_from_fname, iomode2str, open_file, get_unit
 use m_fstrings,      only : sjoin, itoa, endswith
 use m_numeric_tools, only : interpolate_denpot
 use m_geometry,      only : metric
 use m_mpinfo,        only : destroy_mpi_enreg, ptabs_fourdp, initmpi_seq
 use m_distribfft,    only : init_distribfft_seq
 use m_fourier_interpol,only : fourier_interpol

 implicit none

#ifdef HAVE_MPI1
 include 'mpif.h'
#endif

 private

 public :: ioarr                     ! Read or write rho(r) or v(r), either ground-state or response-functions.
 public :: fftdatar_write            ! Write an array in real space. IO library is automatically selected
                                     ! from the file extension and the number of FFT processors:
 public :: fftdatar_write_from_hdr   ! Write an array in real-space to file plus crystal_t and ebands_t
 public :: read_rhor                 ! Read rhor from DEN file.
 public :: fort_denpot_skip          ! Skip the header and the DEN/POT records (Fortran format)

 private :: denpot_spin_convert      ! Convert a density/potential from a spin representation to another

CONTAINS  !====================================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_ioarr/ioarr
!!
!! NAME
!! ioarr
!!
!! FUNCTION
!! Read or write rho(r) or v(r), either ground-state or response-functions.
!! If ground-state, these arrays are real, if response-functions, these arrays are complex.
!! (in general, an array stored in unformatted form on a real space fft grid).
!! rdwr=1 to read, 2 to write
!!
!! This subroutine should be called only by one processor in the writing mode
!!
!! INPUTS
!! (some may be output)
!! accessfil=
!!    0 for FORTRAN_IO
!!    3 for ETSF_IO
!!    4 for MPI_IO
!! cplex=1 for real array, 2 for complex
!! nfft=Number of FFT points treated by this node.
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! dtset <type(dataset_type)>=all input variables for this dataset
!! fform=integer specification for data type:
!!   2 for wf; 52 for density; 102 for potential
!!   old format (prior to ABINITv2.0): 1, 51 and 101.
!! fildata=file name
!! hdr <type(hdr_type)>=the header of wf, den and pot files
!!  if rdwr=1 , used to compare with the hdr of the read disk file
!!  if rdwr=2 , used as the header of the written disk file
!! mpi_enreg=information about MPI parallelization
!! rdwr=choice parameter, see above
!! rdwrpaw=1 only if rhoij PAW quantities have to be read (if rdwr=1)
!! [single_proc]=True if only ONE MPI process is calling this routine. This usually happens when
!!   master calls ioarr to read data that is then broadcasted in the caller. Default: False.
!!   Note that singleproc is not compatible with FFT parallelism because nfft is assumed to be
!!   the total number of points in the FFT mesh.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! arr(cplex*nfft,nspden)=array on real space grid, returned for rdwr=1, input for rdwr=2
!! etotal=total energy (Ha), returned for rdwr=1
!! === if rdwrpaw/=0 ===
!!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!
!! PARENTS
!!      gstate,outscfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine ioarr(accessfil,arr,dtset,etotal,fform,fildata,hdr,mpi_enreg, &
&                ngfft,cplex,nfft,pawrhoij,rdwr,rdwrpaw,wvl_den,single_proc)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accessfil,cplex,nfft,rdwr,rdwrpaw
 integer,intent(inout) :: fform
 real(dp),intent(inout) :: etotal
 character(len=*),intent(in) :: fildata
 logical,optional,intent(in) :: single_proc
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(wvl_denspot_type),optional, intent(in) :: wvl_den
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(inout),target :: arr(cplex*nfft,dtset%nspden)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
#ifdef HAVE_NETCDF
 integer :: ncid,ncerr
 character(len=fnlen) :: file_etsf
#endif
#ifdef HAVE_BIGDFT
 integer :: i,i1,i2,i3,ia,ind,n1,n2,n3
 integer :: zindex,zstart,zstop
#endif
!scalars
 integer,parameter :: master=0
 logical,parameter :: ALLOW_FFTINTERP=.True.
 logical :: need_fftinterp,icheck_fft,qeq0
 integer :: in_unt,out_unt,nfftot_in,nfftot_out,nspden,ncplxfft
 integer :: iomode,fform_dum,iarr,ierr,ispden,me,me_fft,comm_fft
 integer :: comm_cell,usewvl,unt
 integer :: restart,restartpaw,spaceComm,spaceComm_io
 real(dp) :: cputime,walltime,gflops
 character(len=500) :: message,errmsg
 character(len=fnlen) :: my_fildata
 character(len=nctk_slen) :: varname
 type(hdr_type) :: hdr0
 type(wffile_type) :: wff
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: ngfft_in(18),ngfft_out(18)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:),fftn3_distrib(:),ffti3_local(:)
 real(dp), ABI_CONTIGUOUS pointer :: arr_file(:,:),my_density(:,:)
 real(dp),allocatable :: rhor_file(:,:),rhog_in(:,:),rhor_out(:,:),rhog_out(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 ncplxfft = cplex*nfft

 restartpaw=0
 my_fildata = fildata
 nspden = dtset%nspden; usewvl = dtset%usewvl

!Check validity of arguments--only rho(r) (51,52) and V(r) (101,102) are presently supported
 if ( (fform-1)/2 /=25 .and. (fform-1)/2 /=50 ) then
   write(message,'(a,i0,a)')' Input fform= ',fform,' not allowed.'
   MSG_BUG(message)
 end if

!Print input fform
 if ( (fform-1)/2==25 .and. rdwr==1) then
   message = ' ioarr: reading density data '
 else if ( (fform-1)/2==25 .and. rdwr==2) then
   message = ' ioarr: writing density data'
 else if ( (fform-1)/2==50 .and. rdwr==1) then
   message = ' ioarr: reading potential data'
 else if ( (fform-1)/2==50 .and. rdwr==2) then
   message = ' ioarr: writing potential data'
 end if
 call wrtout(std_out,message,'COLL')

 call wrtout(std_out, 'ioarr: file name is: '//TRIM(fildata),'COLL')

#ifdef HAVE_NETCDF
 if (accessfil == IO_MODE_ETSF) then ! Initialize filename in case of ETSF file.
   file_etsf = nctk_ncify(fildata)
   call wrtout(std_out,sjoin('file name for ETSF access: ', file_etsf),'COLL')
 end if
#endif

!Some definitions for MPI-IO access
 spaceComm = mpi_enreg%comm_cell
 comm_cell = mpi_enreg%comm_cell
 comm_fft = mpi_enreg%comm_fft

 if (accessfil == 4) then
   iomode=IO_MODE_MPI
   if (rdwr==1) then
     spaceComm=mpi_enreg%comm_cell
   else
     spaceComm=mpi_enreg%comm_fft
   end if
   me=xmpi_comm_rank(spaceComm)
   if (mpi_enreg%nproc_fft>1) then
     me_fft=mpi_enreg%me_fft
     spaceComm_io=mpi_enreg%comm_fft
   else
     me_fft=0
     spaceComm_io=xmpi_comm_self
   end if
 end if
 if (usewvl==1) then
   spaceComm=mpi_enreg%comm_cell
   me=xmpi_comm_rank(spaceComm)
 end if

 ! Change communicators and ranks if we are calling ioarr with one single processor.
 if (present(single_proc)) then
   if (single_proc) then
     spaceComm = xmpi_comm_self
     spaceComm_io = xmpi_comm_self
     ABI_CHECK(mpi_enreg%nproc_fft == 1, "single_proc cannot be used when nproc_fft > 1")
     comm_cell = xmpi_comm_self
     comm_fft = xmpi_comm_self
     me = 0
   end if
 end if

!=======================================
!Handle input from disk file
!=======================================

 call cwtime(cputime, walltime, gflops, "start")

 if (rdwr==1) then
   if (accessfil == 0 .or. accessfil == 4) then

     ! Here master checks if the input rho(r) is given on a FFT mesh that quals
     ! the one used in the run. If not, we perform a Fourier interpolation, we write the
     ! interpolated rho(r) to a temporary file and we use this file to restart.
     if (ALLOW_FFTINTERP .and. usewvl==0) then
       need_fftinterp = .False.; icheck_fft = .True.
       ! only master checks the FFT mesh if MPI-IO. All processors read ngfft if Fortran-IO
       ! Note that, when Fortran-IO is used, we don't know if the routine is called
       ! by a single processor or by all procs in comm_cell hence we cannot broadcast my_fildata
       ! inside spaceComm as done if accessfil == 4
       if (accessfil == 4) icheck_fft = (xmpi_comm_rank(spaceComm)==master)

       if (icheck_fft) then
         if (open_file(fildata,message,newunit=in_unt,form='unformatted',status='old') /= 0) then
           MSG_ERROR(message)
         end if

         call hdr_io(fform_dum,hdr0,rdwr,in_unt)
         need_fftinterp = (ANY(hdr%ngfft/=hdr0%ngfft) )
         qeq0=(hdr%qptn(1)**2+hdr%qptn(2)**2+hdr%qptn(3)**2<1.d-14)
         ! FIXME: SHould handle double-grid if PAW
         nfftot_in = product(hdr0%ngfft(1:3))
         nfftot_out = product(hdr%ngfft(1:3))

         if (need_fftinterp) then
           write(message, "(2a,2(a,3(i0,1x)))")&
&           "Will perform Fourier interpolation since in and out ngfft differ",ch10,&
&           "ngfft in file: ",hdr0%ngfft,", expected ngfft: ",hdr%ngfft
           MSG_WARNING(message)

           ! Read rho(r) from file, interpolate it, write data and change fildata
           ABI_MALLOC(rhor_file, (cplex*nfftot_in, hdr0%nspden))
           ABI_MALLOC(rhog_in, (2, nfftot_in))
           ABI_MALLOC(rhor_out, (cplex*nfftot_out, hdr0%nspden))
           ABI_MALLOC(rhog_out, (2, nfftot_out))

           do ispden=1,hdr0%nspden
             read(in_unt, err=10, iomsg=errmsg) (rhor_file(iarr,ispden), iarr=1,cplex*nfftot_in)
           end do

           ngfft_in = dtset%ngfft; ngfft_out = dtset%ngfft
           ngfft_in(1:3) = hdr0%ngfft(1:3); ngfft_out(1:3) = hdr%ngfft(1:3)
           ngfft_in(4:6) = hdr0%ngfft(1:3); ngfft_out(4:6) = hdr%ngfft(1:3)
           ngfft_in(9:18) = 0; ngfft_out(9:18) = 0
           ngfft_in(10) = 1; ngfft_out(10) = 1

           call initmpi_seq(MPI_enreg_seq)
           ! Which one is coarse? Note that this part is not very robust and can fail!
           if (ngfft_in(2) * ngfft_in(3) < ngfft_out(2) * ngfft_out(3)) then
             call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfft_in(2),ngfft_in(3),'all')
             call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',ngfft_out(2),ngfft_out(3),'all')
           else
             call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',ngfft_in(2),ngfft_in(3),'all')
             call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfft_out(2),ngfft_out(3),'all')
           end if

           call fourier_interpol(cplex,hdr0%nspden,0,0,nfftot_in,ngfft_in,nfftot_out,ngfft_out,&
            MPI_enreg_seq,rhor_file,rhor_out,rhog_in,rhog_out)

           call destroy_mpi_enreg(MPI_enreg_seq)

           ! MG Hack: Change fildata so that we will use this file to read the correct rho(r)
           ! FIXME: This should be done in a cleaner way!
           my_fildata = trim(fildata)//"__fftinterp_rhor__"
           if (my_fildata == fildata) my_fildata = "__fftinterp_rhor__"
           if (open_file(my_fildata,message,newunit=out_unt,form='unformatted',status='unknown') /= 0) then
             MSG_ERROR(message)
           end if
           call hdr_io(fform_dum,hdr,2,out_unt)
           do ispden=1,hdr0%nspden
             write(out_unt, err=10, iomsg=errmsg) (rhor_out(iarr,ispden),iarr=1,cplex*nfftot_out)
           end do
           close(out_unt)

           ABI_FREE(rhor_file)
           ABI_FREE(rhog_in)
           ABI_FREE(rhor_out)
           ABI_FREE(rhog_out)
         end if ! need_fftinterp

         call hdr_free(hdr0)
         close(in_unt, err=10, iomsg=errmsg)
       end if ! master
       if (accessfil == 4) call xmpi_bcast(my_fildata,master,spaceComm,ierr)
     end if

     if (accessfil == 4) then
       unt = get_unit()
       call WffOpen(iomode,spaceComm,my_fildata,ierr,wff,0,me,unt,spaceComm_io)
       call hdr_io(fform_dum,hdr0,rdwr,wff)
       ! Compare the internal header and the header from the file
       call hdr_check(fform,fform_dum,hdr,hdr0,'COLL',restart,restartpaw)

     else
       if (open_file(my_fildata, message, newunit=unt, form="unformatted", status="old", action="read") /= 0) then
         MSG_ERROR(message)
       end if
       ! Initialize hdr0, thanks to reading of unwff1
       call hdr_io(fform_dum,hdr0,rdwr,unt)
       ! Compare the internal header and the header from the file
       call hdr_check(fform,fform_dum,hdr,hdr0,'COLL',restart,restartpaw)
     end if
     etotal=hdr0%etot

     ! NOTE: should check that restart is possible !!
     !call ptabs_fourdp(mpi_enreg,ngfft(2),ngfft(3),fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

     ! If nspden[file] /= nspden, need a temporary array
     if (hdr0%nspden/=nspden) then
       ABI_MALLOC(arr_file,(cplex*nfft,hdr0%nspden))
     else
       arr_file => arr
     end if

     ! Read data
     do ispden=1,hdr0%nspden
       if(accessfil == 4) then
         call xderiveRRecInit(wff,ierr)
         call xderiveRead(wff,arr(1:ncplxfft,ispden),ncplxfft,spaceComm_io,ierr)
         call xderiveRRecEnd(wff,ierr)

         !call xmpio_read_dp(mpi_fh,offset,sc_mode,ncount,buf,fmarker,mpierr,advance)
         !do idat=1,ndat
         !  do i3=1,n3
         !    if( fftn3_distrib(i3) == me_fft) then
         !      i3_local = ffti3_local(i3)
         !      i3_ldat = i3_local + (idat - 1) * nd3proc
         !      do i2=1,n2
         !        frbase=n1*(i2-1+n2*(i3_local-1)) + (idat - 1) * nfft
         !        do i1=1,n1
         !          fofr(i1+frbase)=workr(1,i1,i2,i3_ldat)
         !        end do
         !      end do
         !    end if
         !  end do
         !end do

       else
         read(unt, err=10, iomsg=errmsg) (arr_file(iarr,ispden),iarr=1,ncplxfft)
       end if
     end do

     if (accessfil == 4) then
       call wffclose(wff,ierr)
     else
       close (unit=unt, err=10, iomsg=errmsg)
     end if

#ifdef HAVE_NETCDF
   else if (accessfil == 3) then

     ! Read the header and broadcast it in comm_cell
     ! FIXME: Use xmpi_comm_self for the time-being because, in loper, ioarr
     ! is called by me==0
     call hdr_read_from_fname(hdr0, file_etsf, fform_dum, comm_cell)
     !call hdr_read_from_fname(hdr0, file_etsf, fform_dum, xmpi_comm_self)
     ABI_CHECK(fform_dum/=0, "hdr_read_from_fname returned fform 0")

     ! Compare the internal header and the header from the file
     call hdr_check(fform, fform_dum, hdr, hdr0, 'COLL', restart, restartpaw)

     ! If nspden[file] /= nspden, need a temporary array
     if (hdr0%nspden/=nspden) then
       ABI_MALLOC(arr_file,(cplex*nfft,hdr0%nspden))
     else
       arr_file => arr
     end if

     if (usewvl == 1) then
       ! Read the array
       if (fform==52) then ! density
         varname = "density"
       else if (fform==102) then ! all potential forms!!!!
         varname = "exchange_correlation_potential"
       end if

       ! Open the file
       NCF_CHECK(nctk_open_read(ncid, file_etsf, xmpi_comm_self))
       NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, varname), arr_file))
       NCF_CHECK(nf90_close(ncid))
     else
       ! Get MPI-FFT tables from input ngfft
       call ptabs_fourdp(mpi_enreg,ngfft(2),ngfft(3),fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

       ! Get the name of the netcdf variable from the ABINIT extension and read data.
       varname = varname_from_fname(file_etsf)
       ncerr = nctk_read_datar(file_etsf,varname,ngfft,cplex,nfft,hdr0%nspden,comm_fft,fftn3_distrib,ffti3_local,arr)
       NCF_CHECK(ncerr)
     end if
#endif

   else
     write(message,'(a,i0,a)')'Bad value for accessfil', accessfil, ' on read '
     MSG_BUG(message)
   end if

   call wrtout(std_out,sjoin("data read from disk file: ", fildata),'COLL')

   etotal=hdr0%etot

   ! Possibly need to convert the potential/density spin components
   if (hdr0%nspden/=nspden) then
     call denpot_spin_convert(arr_file,hdr0%nspden,arr,nspden,fform)
     ABI_FREE(arr_file)
   end if

!  Eventually copy (or distribute) PAW data
   if (rdwrpaw==1.and.restartpaw/=0) then
     if (size(hdr0%pawrhoij) /= size(pawrhoij)) then
       call pawrhoij_copy(hdr0%pawrhoij,pawrhoij,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab, &
                          keep_nspden=.true.)
     else
       call pawrhoij_copy(hdr0%pawrhoij,pawrhoij,keep_nspden=.true.)
     end if
   end if

   if (accessfil == 0 .or. accessfil == 3 .or. accessfil == 4) call hdr_free(hdr0)

 ! =======================================
 ! Set up for writing data
 ! =======================================
 else if (rdwr==2) then

!  In the wavelet case (isolated boundary counditions), the
!  arr array has a buffer that we need to remove.
   if (usewvl == 1) then
#ifdef HAVE_BIGDFT
     zindex = wvl_den%denspot%dpbox%nscatterarr(me, 3)
     if (wvl_den%denspot%rhod%geocode == 'F') then
       n1 = (wvl_den%denspot%dpbox%ndims(1) - 31) / 2
       n2 = (wvl_den%denspot%dpbox%ndims(2) - 31) / 2
       n3 = (wvl_den%denspot%dpbox%ndims(3) - 31) / 2
       zstart = max(15 - zindex, 0)
       zstop  = wvl_den%denspot%dpbox%nscatterarr(me, 2) + &
&       wvl_den%denspot%dpbox%nscatterarr(me, 4) - &
&       max(zindex + wvl_den%denspot%dpbox%nscatterarr(me, 2) &
&       - 2 * n3 - 15, 0)
     else
       MSG_ERROR('ioarr: WVL not implemented yet.')
     end if
     if (zstop - zstart + 1 > 0) then
!      Our slab contains (zstop - zstart + 1) elements
       ABI_ALLOCATE(my_density,((n1*2)*(n2*2)*(zstop-zstart),nspden))
!      We copy the data except the buffer to my_density
       ind = 0

       do i3 = zstart, zstop - 1, 1
         ia = (i3 - 1) * dtset%ngfft(1) * dtset%ngfft(2)
         do i2 = 0, 2 * n2 - 1, 1
           i = ia + (i2 + 14) * dtset%ngfft(1) + 14
           do i1 = 0, 2 * n1 - 1, 1
             i   = i + 1
             ind = ind + 1
             my_density(ind, :) = arr(i, :)
           end do
         end do
       end do
     else
       nullify(my_density)
     end if
#else
     BIGDFT_NOTENABLED_ERROR()
     if(.false. .and. present(wvl_den))then
       write(std_out,*)' One should not be here'
     endif
#endif
   end if

   ! Make sure ngfft agrees with hdr%ngfft.
   if (usewvl == 0) then
     if (any(ngfft(:3) /= hdr%ngfft(:3))) then
       write(message,"(2(a,3(1x,i0)))")"input ngfft: ",ngfft(:3),"differs from  hdr%ngfft: ",hdr%ngfft(:3)
       MSG_ERROR(message)
     end if
   end if

   if (accessfil == 0 .or. accessfil == 4) then
     if(accessfil == 4) then
       unt = get_unit()
       call WffOpen(iomode,spaceComm,fildata,ierr,wff,0,me,unt)
       call hdr_io(fform,hdr,rdwr,wff)
     else
       if (open_file(fildata, message, newunit=unt, form='unformatted', status='unknown', action="write") /= 0) then
         MSG_ERROR(message)
       end if

       ! Write header
       call hdr_io(fform,hdr,rdwr,unt)
     end if

     ! Write actual data
     do ispden=1,nspden
       if(accessfil == 4) then
         call xderiveWRecInit(wff,ierr,me_fft)
         call xderiveWrite(wff,arr(1:ncplxfft,ispden),ncplxfft,spaceComm_io,ierr)
         call xderiveWRecEnd(wff,ierr,me_fft)
       else
         if (usewvl == 0) then
           write(unt, err=10, iomsg=errmsg) (arr(iarr,ispden),iarr=1,ncplxfft)
         else
           write(unt, err=10, iomsg=errmsg) (my_density(iarr,ispden),iarr=1,size(my_density, 1))
         end if
       end if
     end do

     if(accessfil == 4) then
       call WffClose(wff,ierr)
     else
       close(unt, err=10, iomsg=errmsg)
     end if

#ifdef HAVE_NETCDF
   else if ( accessfil == 3 ) then

     ! Master in comm_fft creates the file and writes the header.
     if (xmpi_comm_rank(comm_fft) == 0) then
       call hdr_write_to_fname(hdr, file_etsf, fform)
     end if
     call xmpi_barrier(comm_fft)

     ! Write the array
     if (usewvl == 0) then
       ! Get MPI-FFT tables from input ngfft
       call ptabs_fourdp(mpi_enreg,ngfft(2),ngfft(3),fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

       varname = varname_from_fname(file_etsf)
       ncerr = nctk_write_datar(varname,file_etsf,ngfft,cplex,nfft,nspden,comm_fft,fftn3_distrib,ffti3_local,arr)
       NCF_CHECK(ncerr)
     else
       NCF_CHECK(nctk_open_modify(ncid, file_etsf, xmpi_comm_self))

       if (fform==52) then ! density
         varname = "density"
         if (usewvl == 0) then
           NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, varname), arr))
         else
           NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, varname), my_density))
         end if
       else if (fform==102) then ! all potential forms!!!!
         varname = "exchange_correlation_potential"
         NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, varname), arr))
       end if

       NCF_CHECK(nf90_close(ncid))
     end if
#endif

   else
     write(message,'(a,i0,a)')'Bad value for accessfil', accessfil, ' on write '
     MSG_ERROR(message)
   end if

   if (usewvl == 1 .and. associated(my_density)) then
     ABI_DEALLOCATE(my_density)
   end if

   call wrtout(std_out,sjoin(' Data written to disk file:', fildata),'COLL')

 else
   write(message,'(a,i0,a)')'Called with rdwr = ',rdwr,' not allowed.'
   MSG_BUG(message)
 end if

 call cwtime_report(" IO operation", cputime, walltime, gflops)

 DBG_EXIT("COLL")

 return

 ! Handle Fortran IO error
10 continue
 MSG_ERROR(errmsg)

end subroutine ioarr
!!***

!----------------------------------------------------------------------

!!****f* m_ioarr/fftdatar_write
!! NAME
!! fftdatar_write
!!
!! FUNCTION
!! Write an array in real space on the FFT box to file.
!! The array can be real or complex depending on cplex
!! IO library is automatically selected from the file extension and the number of FFT processors:
!!
!!   1) If path ends with ".nc", the netcdf library is used else Fortran format.
!!
!!   2) If nproc_fft > 1, parallel IO is used (if available)
!!
!! INPUTS
!! varname=Name of the variable to write (used if ETSF-IO).
!! path=File name
!! iomode=
!! hdr <type(hdr_type)>=the header of wf, den and pot files
!! crystal<crystal_t>= data type gathering info on symmetries and unit cell (used if etsf_io)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! cplex=1 for real array, 2 for complex
!! nfft=Number of FFT points treated by this node.
!! nspden=Number of spin-density components.
!! datar(cplex*nfft,nspden)=array on the real space FFT grid.
!! mpi_enreg=information about MPI parallelization
!! [ebands]<ebands_t>=data type with energies and occupations (used if etsf_io)
!!
!! OUTPUT
!!  Only writing
!!
!! NOTES
!!   The string passed to fftdatar_write (first argument) gives the name used to store the data in the netcdf file
!!   The function varname_from_fname defined in the module m_hdr.F90 gives the mapping between the Abinit
!!   file extension and the netcdf name e.g. foo_VHXC.nc --> vxc
!!   This function is used in cut3d so that we can immediately select the data to analyze without having
!!   to prompt the user
!!   Remember to update varname_from_fname if you add a new file or if you change the name of the variable.
!!
!!   fform i.e. the integer specification for data type is automatically initialized from varname.
!!
!! PARENTS
!!      cut3d,m_ioarr,outscfcv,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine fftdatar_write(varname,path,iomode,hdr,crystal,ngfft,cplex,nfft,nspden,datar,mpi_enreg,ebands)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,cplex,nfft,nspden
 character(len=*),intent(in) :: varname,path
 type(hdr_type),intent(inout) :: hdr
 type(crystal_t),intent(in) :: crystal
 type(ebands_t),optional,intent(in) :: ebands
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(inout) :: datar(cplex*nfft,nspden)
 !type(pawrhoij_type),optional,intent(inout) :: pawrhoij_all(hdr%usepaw*crystal%natom)

!Local variables-------------------------------
!!scalars
 integer,parameter :: master=0
 integer :: n1,n2,n3,comm_fft,nproc_fft,me_fft,iarr,ierr,ispden,unt,mpierr,fform
 integer :: i3_glob,my_iomode
 integer(kind=XMPI_OFFSET_KIND) :: hdr_offset,my_offset,nfft_tot
#ifdef HAVE_NETCDF
 integer :: ncid,ncerr
 character(len=fnlen) :: file_etsf
#endif
 real(dp) :: cputime,walltime,gflops
 character(len=500) :: msg,errmsg
 type(abifile_t) :: abifile
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:),fftn3_distrib(:),ffti3_local(:)
 integer(XMPI_OFFSET_KIND) :: bsize_frecord(nspden)

! *************************************************************************

 abifile = abifile_from_varname(varname)
 if (abifile%fform == 0) then
    MSG_ERROR(sjoin("Cannot find any abifile object associated to varname:", varname))
 end if
 ! Get fform from abifile. TODO: check file extension
 fform = abifile%fform

 comm_fft = mpi_enreg%comm_fft; nproc_fft = xmpi_comm_size(comm_fft); me_fft = mpi_enreg%me_fft
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); nfft_tot = n1*n2*n3

 ! Select iomode
 ! Use Fortran IO if nproc_fft 1, in principle this is not needed because the
 ! MPI-IO code should produce binary files that are readable with Fortran-IO
 ! but it seems that NAG uses its own binary format
 my_iomode = iomode
 if (my_iomode /= IO_MODE_ETSF .and. nproc_fft == 1) my_iomode = IO_MODE_FORTRAN
 if (nproc_fft > 1 .and. my_iomode == IO_MODE_FORTRAN) my_iomode = IO_MODE_MPI

 call wrtout(std_out, sjoin(ch10, "  fftdatar_write: About to write data to:", path, "with iomode:",iomode2str(my_iomode)))
 call cwtime(cputime, walltime, gflops, "start")

 ! Get MPI-FFT tables from input ngfft
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 select case (my_iomode)
 case (IO_MODE_FORTRAN)
   ABI_CHECK(nproc_fft == 1, "MPI-IO must be enabled when FFT parallelism is used")
   if (open_file(path, msg, newunit=unt, form='unformatted', status='unknown', action="write") /= 0) then
     MSG_ERROR(msg)
   end if
   call hdr_fort_write(hdr,unt,fform,ierr)
   ABI_CHECK(ierr==0,"ierr !=0")
   do ispden=1,nspden
     write(unt, err=10, iomsg=errmsg) (datar(iarr,ispden), iarr=1,cplex * nfft)
   end do
   close(unt, err=10, iomsg=errmsg)

   ! Write PAW rhoij
   !call pawrhoij_io(hdr%pawrhoij,unit,hdr%nsppol,hdr%nspinor,hdr%nspden,hdr%lmn_size,hdr%typat,headform,"Write")
   !call pawrhoij_io(rhoij_ptr,ncid,hdr%nsppol,hdr%nspinor,hdr%nspden,hdr%lmn_size,hdr%typat,&
   !                 HDR_LATEST_HEADFORM,"Write",form="netcdf")

#ifdef HAVE_MPI_IO
 case (IO_MODE_MPI)
   ! Find the first z-plane treated by this node.
   ! WARNING: Here I assume that the z-planes in real space
   ! are distributed in contiguous blocks (as usually done in MPI-FFT)
   do i3_glob=1,n3
     if (me_fft == fftn3_distrib(i3_glob)) exit
   end do
   ABI_CHECK(i3_glob /= n3 +1, "This processor does not have z-planes!")

   ! Master writes the header.
   if (me_fft == master) call hdr_write_to_fname(hdr,path,fform)
   call xmpi_barrier(comm_fft) ! TODO: Non-blocking barrier.

   call MPI_FILE_OPEN(comm_fft, path, MPI_MODE_RDWR, xmpio_info, unt, mpierr)
   ABI_CHECK_MPI(mpierr,"MPI_FILE_OPEN")

   ! Skip the header and get the offset of the header
   call hdr_mpio_skip(unt,fform,hdr_offset)
   !write(std_out,*)"i3_glob, nfft, hdr_offset,",i3_glob,nfft,hdr_offset,fftn3_distrib == me_fft

   ! Each proc writes a contiguous slice of the nspden records.
   ! my_offset is the position inside the Fortran record.
   do ispden=1,nspden
     my_offset = hdr_offset + xmpio_bsize_frm + ((ispden - 1) * 2 * xmpio_bsize_frm) + &
     ((i3_glob-1) * cplex * n1 * n2 * xmpi_bsize_dp)  + ((ispden-1) * cplex * nfft_tot * xmpi_bsize_dp)
     call MPI_FILE_WRITE_AT_ALL(unt,my_offset,datar(:,ispden),cplex*nfft,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,mpierr)
     ABI_CHECK_MPI(mpierr,"MPI_FILE_WRITE_AT_ALL")
   end do

   ! master writes the fortran record markers.
   if (me_fft == master) then
     bsize_frecord = cplex * nfft_tot * xmpi_bsize_dp
#if 1
     my_offset = hdr_offset
     do ispden=1,nspden
       call xmpio_write_frm(unt,my_offset,xmpio_single,bsize_frecord(ispden),mpierr)
       ABI_CHECK_MPI(mpierr,"xmpio_write_frm")
     end do
#else
     ! TODO: Understand why this code does not work!
     call xmpio_write_frmarkers(unt,hdr_offset,xmpio_single,nspden,bsize_frecord,ierr)
     ABI_CHECK(ierr==0, "xmpio_write_frmarkers")
#endif
   end if

   call MPI_FILE_CLOSE(unt,mpierr)
   ABI_CHECK_MPI(mpierr,"FILE_CLOSE!")

   ! Add full pawrhoij datastructure at the end of the file.
   !if (present(pawrhoij_all) .and. me_fft == master .and. hdr%usepaw == 1) then
   !  if (open_file(path, msg, newunit=unt, form='unformatted', status='old', action="write", access="append") /= 0) then
   !    MSG_ERROR(msg)
   !  end if
   !  call pawrhoij_io(pawrhoij_all,un,hdr%nsppol,hdr%nspinor,hdr%nspden,hdr%lmn_size,hdr%typat,hdr%headform,"Write")
   !  close(unt)
   !end if
#endif

#ifdef HAVE_NETCDF
 case (IO_MODE_ETSF)
   file_etsf = nctk_ncify(path)

   ! Write datar.
   ncerr = nctk_write_datar(varname,file_etsf,ngfft,cplex,nfft,nspden, &
   comm_fft,fftn3_distrib,ffti3_local,datar,action="create")
   NCF_CHECK(ncerr)
   call xmpi_barrier(comm_fft)

   ! Master writes the header.
   if (xmpi_comm_rank(comm_fft) == master) then
     NCF_CHECK(nctk_open_modify(ncid, file_etsf, xmpi_comm_self))
     NCF_CHECK(hdr_ncwrite(hdr, ncid, fform, nc_define=.True.))
     ! Add information on the crystalline structure.
     NCF_CHECK(crystal%ncwrite(ncid))
     if (present(ebands)) then
       NCF_CHECK(ebands_ncwrite(ebands, ncid))
     end if

     ! Add full pawrhoij datastructure.
     !if (present(pawrhoij_all) .and. me_fft == master .and. hdr%usepaw == 1) then
     !  call pawrhoij_io(pawrhoij_all,ncid,hdr%nsppol,hdr%nspinor,hdr%nspden,hdr%lmn_size,hdr%typat,hdr%headform,"Write", form="netcdf")
     !end if

     NCF_CHECK(nf90_close(ncid))
   end if
#endif

 case default
   MSG_ERROR(sjoin("Wrong iomode:",itoa(my_iomode)))
 end select

 call cwtime_report(" IO operation", cputime, walltime, gflops)

 return

 ! Handle Fortran IO error
10 continue
 MSG_ERROR(errmsg)

end subroutine fftdatar_write
!!***

!----------------------------------------------------------------------

!!****f* m_ioarr/fftdatar_write_from_hdr
!! NAME
!! fftdatar_write_from_hdr
!!
!! FUNCTION
!! Write an array in real space on the FFT box to file.
!! crystal and ebands are constructed from the Abinit header.
!!
!! TODO
!! This routine will be removed when crystal_t and ebands_t will become standard objects
!! available in the GS/DFPT part.
!!
!! INPUTS
!! [eigen](mband*hdr%nkpt*hdr%nsppol)=GS eigenvalues
!! See fftdatar_write for the meaning of the other variables.
!!
!! OUTPUT
!!
!! PARENTS
!!      dfpt_scfcv,scfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine fftdatar_write_from_hdr(varname,path,iomode,hdr,ngfft,cplex,nfft,nspden,datar,mpi_enreg,eigen)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iomode,cplex,nfft,nspden
 character(len=*),intent(in) :: varname,path
 type(hdr_type),intent(inout) :: hdr
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(inout) :: datar(cplex*nfft,nspden)
 real(dp),optional,intent(in) :: eigen(:)

!Local variables-------------------------------
!!scalars
 integer :: timrev,mband
 type(crystal_t) :: crystal
 type(ebands_t) :: ebands
!arrays
 real(dp),allocatable :: ene3d(:,:,:)

! *************************************************************************

 timrev = 2; if (any(hdr%kptopt == [3, 4])) timrev = 1
 crystal = hdr_get_crystal(hdr, timrev)

 if (present(eigen)) then
     mband = maxval(hdr%nband)
     ABI_CHECK(size(eigen) ==  mband * hdr%nkpt * hdr%nsppol, "Wrong size(eigen)")
     ABI_MALLOC(ene3d, (mband, hdr%nkpt, hdr%nsppol))
     call unpack_eneocc(hdr%nkpt, hdr%nsppol, mband, hdr%nband, eigen, ene3d)
     ebands = ebands_from_hdr(hdr, mband, ene3d)
     ABI_FREE(ene3d)

    call fftdatar_write(varname,path,iomode,hdr,crystal,ngfft,cplex,nfft,nspden,datar,mpi_enreg,ebands=ebands)
    call ebands_free(ebands)
 else
    call fftdatar_write(varname,path,iomode,hdr,crystal,ngfft,cplex,nfft,nspden,datar,mpi_enreg)
 end if

 call crystal%free()

end subroutine fftdatar_write_from_hdr
!!***

!----------------------------------------------------------------------

!!****f* m_ioarr/read_rhor
!! NAME
!! read_rhor
!!
!! FUNCTION
!!  Read the DEN file with name fname reporting the density on the real FFT mesh
!!  specified through the input variable ngfft. If the FFT mesh asked in input and that found
!!  on file differ, the routine performs a FFT interpolation and renormalize the density so that it
!!  integrates to the correct number of electrons. The interpolation is done only for NC.
!!  For PAW, this is not possible because one should include the onsite contribution so this task
!!  is delegated to the caller.
!!
!! INPUTS
!! fname=Name of the density file
!! cplex=1 if density is real, 2 if complex e.g. DFPT density.
!! nspden=Number of spin density components.
!! nfft=Number of FFT points (treated by this processor)
!! ngfft(18)=Info on the FFT mesh.
!! pawread= 1 if pawrhoij should be read from file, 0 otherwise. Meaningful only if usepaw==1.
!! mpi_enreg<MPI_type>=Information about MPI parallelization
!! comm=MPI communicator. See notes
!! [check_hdr] <type(hdr_type)>=Optional. Used to compare with the hdr read from disk file
!!   The routine will abort if restart cannot be performed.
!! [allow_interp]=If True, the density read from file will be interpolated if the mesh differs from the one
!!   expected by the caller. This option is usually used in **self-consistent** calculations.
!!   If False (default), the code stops if the two meshes are different.
!!
!! OUTPUT
!! orhor(cplex*nfft,nspden)=The density on the real space mesh.
!! ohdr=Abinit header read from file.
!! pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data. only
!!   if pawread==1. The arrays is supposed to be already allocated in the caller and its
!!   size must be consistent with the MPI communicator comm.
!!
!! NOTES
!!   if xmpi_comm_size(comm) == 1, nfft shall be equal to nfftot, and len(pawrhoij) == natom
!!   This means that one can call this routine with
!!
!!     if (xmpi_comm_rank(comm) == 0) call read_rhor(...., comm=xmpi_comm_self)
!!
!!   to get the full array and pawrhoij(natom) on the master node.
!!
!!   if xmpi_comm_size(comm) > 1, nfft represents the number of FFT points treated by this processor,
!!   and pawrhoij is dimensioned with my_natom
!!   All the processors inside comm and comm_atom should call this routine.
!!
!! PARENTS
!!      dfpt_looppert,dfptnl_loop,gstate,mrgscr,nonlinear,respfn,setup_positron
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine read_rhor(fname, cplex, nspden, nfft, ngfft, pawread, mpi_enreg, orhor, ohdr, pawrhoij, comm, &
  check_hdr, allow_interp) ! Optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,pawread,comm
 character(len=*),intent(in) :: fname
 type(MPI_type),intent(in) :: mpi_enreg
 type(hdr_type),intent(out) :: ohdr
 type(hdr_type),optional,intent(in) :: check_hdr
 logical,optional,intent(in) :: allow_interp
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(out) :: orhor(cplex*nfft,nspden)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: unt,fform,iomode,my_rank,mybase,globase,cplex_file
!integer :: optin,optout
 integer :: ispden,ifft,nfftot_file,nprocs,ierr,i1,i2,i3,i3_local,n1,n2,n3
 integer,parameter :: fform_den=52
 integer :: restart, restartpaw
#ifdef HAVE_NETCDF
 integer :: ncerr
#endif
 real(dp) :: ratio,ucvol
 real(dp) :: cputime,walltime,gflops
 logical :: need_interp,have_mpifft,allow_interp__
 character(len=500) :: msg,errmsg
 character(len=fnlen) :: my_fname
 character(len=nctk_slen) :: varname
 !type(mpi_type) :: mpi_enreg_seq
!arrays
!integer :: ngfft_file(18)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:),fftn3_distrib(:),ffti3_local(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
!real(dp) :: rhogdum(1)
 real(dp),allocatable :: rhor_file(:,:),rhor_tmp(:,:)
 type(pawrhoij_type),allocatable :: pawrhoij_file(:)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); have_mpifft = (nfft /= product(ngfft(1:3)))
 allow_interp__ = .False.; if (present(allow_interp)) allow_interp__ = allow_interp

 call wrtout(std_out, sjoin(" About to read data(r) from:", fname), do_flush=.True.)
 call cwtime(cputime, walltime, gflops, "start")

 ! Master node opens the file, read the header and the FFT data
 ! This approach facilitates the interpolation of the density if in_ngfft(1:3) /= file_ngfft(1:3)
 if (my_rank == master) then
   my_fname = fname
   if (nctk_try_fort_or_ncfile(my_fname, msg) /= 0 ) then
     MSG_ERROR(msg)
   end if

   iomode = iomode_from_fname(my_fname)
   select case (iomode)

   case (IO_MODE_FORTRAN, IO_MODE_MPI)
     if (open_file(my_fname, msg, newunit=unt, form='unformatted', status='old', action="read") /= 0) then
       MSG_ERROR(msg)
     end if

     call hdr_fort_read(ohdr, unt, fform)

     ! Check important dimensions.
     ABI_CHECK(fform /= 0, sjoin("fform == 0 while reading:", my_fname))
     if (fform /= fform_den) then
       write(msg, "(2a, 2(a, i0))")' File: ',trim(my_fname),' is not a density file. fform: ',fform,", expecting: ", fform_den
       MSG_WARNING(msg)
     end if
     cplex_file = 1
     if (ohdr%pertcase /= 0) then
       cplex_file = 2; if (ohdr%qptn(1)**2 + ohdr%qptn(2)**2 + ohdr%qptn(3)**2 <1.d-14) cplex_file= 1
     end if
     ABI_CHECK(cplex_file == cplex, "cplex_file != cplex")

     ! Read FFT array (full box)
     nfftot_file = product(ohdr%ngfft(:3))
     ABI_MALLOC(rhor_file, (cplex*nfftot_file, ohdr%nspden))
     do ispden=1,ohdr%nspden
       read(unt, err=10, iomsg=errmsg) (rhor_file(ifft,ispden), ifft=1,cplex*nfftot_file)
     end do
     close(unt)

#ifdef HAVE_NETCDF
   case (IO_MODE_ETSF)
     NCF_CHECK(nctk_open_read(unt, my_fname, xmpi_comm_self))
     call hdr_ncread(ohdr, unt, fform)

     ! Check important dimensions.
     ABI_CHECK(fform /= 0, sjoin("fform == 0 while reading:", my_fname))
     if (fform /= fform_den) then
       write(msg, "(2a, 2(a, i0))")' File: ',trim(my_fname),' is not a density file: fform= ',fform,", expecting:", fform_den
       MSG_WARNING(msg)
     end if

     cplex_file = 1
     if (ohdr%pertcase /= 0) then
       cplex_file = 2; if (ohdr%qptn(1)**2 + ohdr%qptn(2)**2 + ohdr%qptn(3)**2 <1.d-14) cplex_file= 1
     end if
     ABI_CHECK(cplex_file == cplex, "cplex_file != cplex")

     ! Read FFT array (full box)
     nfftot_file = product(ohdr%ngfft(:3))
     ABI_MALLOC(rhor_file, (cplex*nfftot_file, ohdr%nspden))

     varname = varname_from_fname(my_fname)
     ncerr= nf90_get_var(unt, nctk_idname(unt, varname), rhor_file, &
                        count=[cplex, ohdr%ngfft(1), ohdr%ngfft(2), ohdr%ngfft(3), ohdr%nspden])
     NCF_CHECK(ncerr)
     NCF_CHECK(nf90_close(unt))
#endif
   case default
     MSG_ERROR(sjoin("Wrong iomode:", itoa(iomode)))
   end select

   need_interp = any(ohdr%ngfft(1:3) /= ngfft(1:3))
   if (need_interp .and. allow_interp__) then
     MSG_COMMENT("Real space meshes (DEN file, input rhor) are different. Will interpolate rhor(r).")

#if 0
     ABI_CHECK(cplex == 1, "cplex != 1 not coded!")
     ngfft_file(1:3) = ohdr%ngfft(1:3)
     ngfft_file(4) = 2*(ngfft_file(1)/2)+1 ! 4:18 are used in fourdp
     ngfft_file(5) = 2*(ngfft_file(2)/2)+1
     ngfft_file(6) = ngfft_file(3)
     ngfft_file(7:18) = ngfft(7:18)
     optin  = 0 ! Input is taken from rhor
     optout = 0 ! Output is only in real space

     ! Fake MPI_type for the sequential part.
     !call initmpi_seq(MPI_enreg_seq)
     !call init_distribfft_seq(MPI_enreg_seq%distribfft, 'c', ngfftc(2), ngfftc(3), 'all')
     !call init_distribfft_seq(MPI_enreg_seq%distribfft, 'f', ngfftf(2), ngfftf(3), 'all')

     call fourier_interpol(cplex,ohdr%nspden,optin,optout,nfftot_file,ngfft_file,nfft,ngfft,&
       mpi_enreg,rhor_file,orhor,rhogdum,rhogdum)

     !call destroy_mpi_enreg(MPI_enreg_seq)
#endif

     ABI_MALLOC(rhor_tmp, (cplex*product(ngfft(1:3)), ohdr%nspden))
     call interpolate_denpot(cplex, ohdr%ngfft(1:3), ohdr%nspden, rhor_file, ngfft(1:3), rhor_tmp)

     ohdr%ngfft(1:3) = ngfft(1:3)
     nfftot_file = product(ohdr%ngfft(:3))
     ABI_FREE(rhor_file)
     ABI_MALLOC(rhor_file, (cplex*nfftot_file, ohdr%nspden))
     rhor_file = rhor_tmp
     ABI_FREE(rhor_tmp)

     ! Renormalize charge to avoid errors due to the interpolation.
     ! Do this only for NC since for PAW we should add the onsite contribution.
     ! This is left to the caller.
     !if (ohdr%usepaw == 0) then
     if (ohdr%usepaw == 0 .and. fform == fform_den) then
       call metric(gmet, gprimd, -1, rmet, ohdr%rprimd, ucvol)
       ratio = ohdr%nelect / (sum(rhor_file(:,1))*ucvol/ product(ngfft(1:3)))
       rhor_file = rhor_file * ratio
       write(msg,'(a,f8.2,a,f8.4)')' Expected nelect: ',ohdr%nelect,' renormalization ratio: ',ratio
       call wrtout(std_out,msg,'COLL')
     end if
   end if ! need_interp

   ! Read PAW Rhoij
   if (ohdr%usepaw == 1) then
     ABI_DT_MALLOC(pawrhoij_file, (ohdr%natom))
     call pawrhoij_nullify(pawrhoij_file)
     call pawrhoij_alloc(pawrhoij_file, ohdr%pawrhoij(1)%cplex_rhoij, ohdr%pawrhoij(1)%nspden, ohdr%pawrhoij(1)%nspinor, &
         ohdr%pawrhoij(1)%nsppol, ohdr%typat, lmnsize=ohdr%lmn_size, qphase=ohdr%pawrhoij(1)%qphase)
     call pawrhoij_copy(ohdr%pawrhoij, pawrhoij_file)
   end if

  ! Check that restart is possible !
  ! This check must be done here because we may have changed hdr% if need_interp
  if (present(check_hdr)) then
    ! FIXME: Temporary hack: fform_den to make hdr_check happy!
    call hdr_check(fform_den, fform_den, check_hdr, ohdr, "COLL", restart, restartpaw)
    !call hdr_check(fform_den, fform, check_hdr, ohdr, "COLL", restart, restartpaw)
  end if

 end if ! master

 if (nprocs == 1) then
   if (ohdr%nspden == nspden) then
     orhor = rhor_file
   else
     call denpot_spin_convert(rhor_file,ohdr%nspden,orhor,nspden,fform)
   end if
   if (pawread == 1) call pawrhoij_copy(pawrhoij_file, pawrhoij, keep_nspden=.true.)

 else
   call hdr_bcast(ohdr, master, my_rank, comm)
   call xmpi_bcast(fform, master, comm, ierr)

   ! Eventually copy (or distribute) PAW data
   if (ohdr%usepaw == 1 .and. pawread == 1) then
     if (my_rank /= master) then
       ABI_DT_MALLOC(pawrhoij_file, (ohdr%natom))
       call pawrhoij_nullify(pawrhoij_file)
       call pawrhoij_alloc(pawrhoij_file, ohdr%pawrhoij(1)%cplex_rhoij, ohdr%pawrhoij(1)%nspden, ohdr%pawrhoij(1)%nspinor, &
&           ohdr%pawrhoij(1)%nsppol, ohdr%typat, lmnsize=ohdr%lmn_size, qphase=ohdr%pawrhoij(1)%qphase)
     end if
     if (size(ohdr%pawrhoij) /= size(pawrhoij)) then
       call pawrhoij_copy(ohdr%pawrhoij,pawrhoij,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab, &
&                         keep_nspden=.true.)
     else
       call pawrhoij_copy(ohdr%pawrhoij,pawrhoij, keep_nspden=.true.)
     end if
   end if

   if (my_rank /= master) then
     nfftot_file = product(ohdr%ngfft(1:3))
     ABI_MALLOC(rhor_file, (cplex*nfftot_file, ohdr%nspden))
   end if
   call xmpi_bcast(rhor_file, master, comm,ierr)

   if (have_mpifft) then
     ! Extract slice treated by this MPI-FFT process.
     call ptabs_fourdp(mpi_enreg, ngfft(2), ngfft(3), fftn2_distrib, ffti2_local, fftn3_distrib, ffti3_local)
     if (ohdr%nspden==nspden) then
       do ispden=1,nspden
         do i3=1,n3
           if (fftn3_distrib(i3) /= mpi_enreg%me_fft) cycle
           i3_local = ffti3_local(i3)
           do i2=1,n2
             mybase = cplex * (n1 * (i2-1 + n2*(i3_local-1)))
             globase = cplex * (n1 * (i2-1 + n2*(i3-1)))
             do i1=1,n1*cplex
               orhor(i1+mybase,ispden) = rhor_file(i1+globase,ispden)
             end do
           end do
         end do
       end do
     else
       do i3=1,n3
         if (fftn3_distrib(i3) /= mpi_enreg%me_fft) cycle
         i3_local = ffti3_local(i3)
         do i2=1,n2
           mybase  = 1 + cplex * (n1 * (i2-1 + n2*(i3_local-1)))
           globase = 1 + cplex * (n1 * (i2-1 + n2*(i3-1)))
           call denpot_spin_convert(rhor_file,ohdr%nspden,orhor,nspden,fform,&
&                  istart_in=globase,istart_out=mybase,nelem=n1*cplex)
         end do
       end do
     end if
   else
     if (ohdr%nspden==nspden) then
       orhor = rhor_file
     else
       call denpot_spin_convert(rhor_file,ohdr%nspden,orhor,nspden,fform)
     end if
   end if
 end if ! nprocs > 1

 ABI_FREE(rhor_file)

 if (allocated(pawrhoij_file)) then
   call pawrhoij_free(pawrhoij_file)
   ABI_DT_FREE(pawrhoij_file)
 end if

! Non-collinear magnetism: avoid zero magnetization, because it produces numerical instabilities
! Add a small real to the magnetization
 if (nspden==4) orhor(:,4)=orhor(:,4)+tol14
 if (ohdr%usepaw==1.and.size(pawrhoij)>0) then
   if (pawrhoij(1)%nspden==4) then
     do i1=1,size(pawrhoij)
       pawrhoij(i1)%rhoijp(:,4)=pawrhoij(i1)%rhoijp(:,4)+tol10
     end do
   end if
 end if

 call cwtime_report(" read_rhor", cputime, walltime, gflops)
 return

 ! Handle Fortran IO error
10 continue
 MSG_ERROR(errmsg)

end subroutine read_rhor
!!***

!----------------------------------------------------------------------

!!****f* m_ioarr/fort_denpot_skip
!! NAME
!!  fort_denpot_skip
!!
!! FUNCTION
!!  Skip the header and the DEN/POT records. Mainly used to append data to a pre-existent file.
!!  Return exit code.
!!
!! INPUTS
!!  unit=Fortran unit number (already opened in the caller).
!!  msg=Error message if ierr /= 0
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function fort_denpot_skip(unit, msg) result(ierr)

!Arguments ------------------------------------
 integer,intent(in) :: unit
 character(len=*),intent(out) :: msg

!Local variables-------------------------------
 integer :: ii,fform,nspden
 type(hdr_type) :: hdr

! *********************************************************************

 ierr = 1
 call hdr_fort_read(hdr, unit, fform)
 if (fform == 0) then
    msg = "hdr_fort_read returned fform == 0"; return
 end if

 nspden = hdr%nspden
 call hdr_free(hdr)

 ! Skip the records with v1.
 do ii=1,nspden
   read(unit, iostat=ierr, iomsg=msg)
   if (ierr /= 0) return
 end do

 ierr = 0

end function fort_denpot_skip
!!***

!----------------------------------------------------------------------

!!****f* m_ioarr/denpot_spin_convert
!! NAME
!!  denpot_spin_convert
!!
!! FUNCTION
!!  Convert a density/potential from a spin representation to another
!!
!! INPUTS
!!  denpot_in(:,nspden_in)=input density//potential
!!  nspden_in=number of spin-component of the input density/potential
!!  fform=file format (density or potential)
!!  [istart_in]= --optional-- starting index in the denpot_in array; default is 1
!!  [istart_out]= --optional-- starting index in the denpot_out array; default is 1
!!  [nelem]= --optional-- number of elements to copy from denpot_in to denpot_out; default is all
!!
!! OUTPUT
!!  denpot_out(:,nspden_out)=output density//potential
!!  nspden_out=number of spin-component of the output density/potential
!!
!! NOTES
!!  More explicitely:
!!    We copy denpot_in(istar_in+1:istart_in+nelem,:)
!!       into denpot_out(istart_out+1:istart_out+nelem,:)
!!
!! PARENTS
!!      m_ioarr
!!
!! CHILDREN
!!
!! SOURCE

subroutine denpot_spin_convert(denpot_in,nspden_in,denpot_out,nspden_out,fform,&
&                              istart_in,istart_out,nelem) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nspden_in,nspden_out,fform
 integer,intent(in),optional :: istart_in,istart_out,nelem
!arrays
 real(dp),intent(in) :: denpot_in(:,:)
 real(dp),intent(out) :: denpot_out(:,:)

!Local variables-------------------------------
 integer :: iend_in,iend_out,ispden,my_istart_in,my_istart_out,my_nelem
 character(len=500) :: msg

! *********************************************************************

!Optional arguments
 my_istart_in=1;if (present(istart_in)) my_istart_in=istart_in
 my_istart_out=1;if (present(istart_out)) my_istart_out=istart_out
 iend_in=size(denpot_in,1) ; iend_out=size(denpot_out,1)
 my_nelem=min(iend_in-my_istart_in+1,iend_out-my_istart_out+1)
 if (present(nelem)) my_nelem=nelem

!Checks
 if (size(denpot_in,2)/=nspden_in) then
   msg='size(denpot_in,2)/=nspden_in!'
   MSG_BUG(msg)
 end if
 if (size(denpot_out,2)/=nspden_out) then
   msg='size(denpot_out,2)/=nspden_out!'
   MSG_BUG(msg)
 end if
 if (my_istart_in+my_nelem-1>size(denpot_in,1)) then
   msg='istart_in+nelem>size(denpot_in,1)!'
   MSG_BUG(msg)
 end if
 if (my_istart_out+my_nelem-1>size(denpot_out,1)) then
   msg='istart_out+nelem>size(denpot_out,1)!'
   MSG_BUG(msg)
 end if

!Simple copy if the number of spin-components is unchanged...
 if (nspden_in==nspden_out) then
   do ispden=1,nspden_in
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,ispden)= &
&      denpot_in(my_istart_in:my_istart_in+my_nelem-1,ispden)
   end do
   return
 end if

!...otherwise, we need to convert.
 if ((fform-1)/2==25) then

!  First case: DENSITY

   if      (nspden_in==1.and.nspden_out==2) then
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,1)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,2)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)*half
   else if (nspden_in==1.and.nspden_out==4) then
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,1)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,2)=zero
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,3)=zero
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,4)=zero
   else if (nspden_in==2.and.nspden_out==1) then
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,1)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)
   else if (nspden_in==2.and.nspden_out==4) then
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,1)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,2)=zero
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,3)=zero
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,4)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,2)*two &
&                                                        -denpot_in(my_istart_in:my_istart_in+my_nelem,1)
   else if (nspden_in==4.and.nspden_out==1) then
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,1)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)
   else if (nspden_in==4.and.nspden_out==2) then
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,1)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,2)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)*half &
&                                                        +denpot_in(my_istart_in:my_istart_in+my_nelem-1,4)*half
   end if

 else

!  Second case: POTENTIAL

   if      (nspden_in==1.and.nspden_out==2) then
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,1)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,2)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)
   else if (nspden_in==1.and.nspden_out==4) then
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,1)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,2)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,3)=zero
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,4)=zero
   else if (nspden_in==2.and.nspden_out==1) then
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,1)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)*half &
&                                                        +denpot_in(my_istart_in:my_istart_in+my_nelem-1,2)*half
   else if (nspden_in==2.and.nspden_out==4) then
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,1)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,2)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,2)
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,3)=zero
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,4)=zero
   else if (nspden_in==4.and.nspden_out==1) then
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,1)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)*half &
&                                                        +denpot_in(my_istart_in:my_istart_in+my_nelem-1,2)*half
   else if (nspden_in==4.and.nspden_out==2) then
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,1)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,1)
     denpot_out(my_istart_out:my_istart_out+my_nelem-1,2)=denpot_in(my_istart_in:my_istart_in+my_nelem-1,2)
   end if

 end if

end subroutine denpot_spin_convert
!!***

end module m_ioarr
!!***
