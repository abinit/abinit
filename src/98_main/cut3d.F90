!!****p* ABINIT/cut3d
!! NAME
!! cut3d
!!
!! FUNCTION
!! Main routine for the analysis of the density and potential files,
!! as well as other files with the ABINIT header.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2019 ABINIT group (GMR, RC, LSI, XG, NCJ, JFB, MCote, LPizzagalli)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main program)
!!
!! OUTPUT
!!  (main program)
!!
!! NOTES
!! natom = number of atoms in the unit cell
!! nr1,nr2,nr3 = grid size (nr1 x nr2 x nr3 = filrho dimension)
!! ntypat = number of atom types
!! ucvol = unit cell volume (> 0)
!! filrho = name of the density file (binary or netcdf)
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,abinit_doctor,crystal_free,crystal_from_hdr
!!      cut3d_hirsh,cut3d_lineint,cut3d_planeint,cut3d_pointint,cut3d_rrho
!!      cut3d_volumeint,cut3d_wffile,destroy_mpi_enreg,fftdatar_write
!!      flush_unit,hdr_echo,hdr_free,hdr_read_from_fname,herald
!!      init_distribfft_seq,initmpi_seq,metric,ngfft_seq,timein,wrtout,xmpi_end
!!      xmpi_init,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program cut3d

 use defs_basis
 use m_errors
 use m_build_info
 use m_xmpi
 use m_nctk
 use m_abicore
#ifdef HAVE_NETCDF
 use netcdf
#endif
#if defined FC_NAG
 use f90_unix_proc
#endif
 use m_hdr
 use m_cut3d
 use m_crystal

 use defs_abitypes,     only : MPI_type
 use m_specialmsg,      only : specialmsg_getcount, herald
 use m_fstrings,        only : endswith, sjoin, itoa
 use m_time,            only : timein
 use m_geometry,        only : xred2xcart, metric
 use m_mpinfo,          only : destroy_mpi_enreg, initmpi_seq
 use m_fftcore,         only : ngfft_seq
 use m_distribfft,      only : init_distribfft_seq
 use m_ioarr,           only : fftdatar_write
 use m_io_tools,        only : flush_unit, file_exists, open_file, is_open, get_unit, read_string
 implicit none

!Local variables-------------------------------
 character(len=1) :: outputchar,blank=' '
!scalars
 integer,parameter :: mfiles=10,exchn2n3d0=0
 integer :: fform0,gridshift1,gridshift2,gridshift3,i1,i2,i3
 integer :: iatom,ifiles,ii,ii1,ii2,ii3,index,iprompt,ir1,ir2,ir3,ispden,cplex
 integer :: itask,jfiles,natom,nfiles,nr1,nr2,unt,comm,iomode,nprocs,my_rank
 integer :: nr3,nr1_stored,nr2_stored,nr3_stored,nrws,nspden,nspden_stored,ntypat,timrev,nfft
 real(dp) :: dotdenpot,maxmz,normz,sumdenpot,ucvol,xm,xnow,xp,ym,ynow,yp,zm,znow,zp,tcpui,twalli
 character(len=24) :: codename
 character(len=fnlen) :: filnam,filrho,filrho_tmp
 character(len=nctk_slen) :: varname
 type(hdr_type) :: hdr
 type(abifile_t) :: abifile
 type(MPI_type) :: mpi_enreg
 type(crystal_t) :: cryst
!arrays
 integer, allocatable :: isdenpot(:)
 integer :: ngfft(18)
 real(dp) :: rprimd(3,3),shift_tau(3),tsec(2)
 real(dp) :: xcart2(3),gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: grid(:,:,:),grid_full(:,:,:,:),grid_full_stored(:,:,:,:,:),gridtt(:,:,:)
 real(dp),allocatable :: tau2(:,:),xcart(:,:),xred(:,:),rhomacu(:,:),gridmz(:,:,:),gridmy(:,:,:),gridmx(:,:,:)
 character(len=fnlen),allocatable :: filrho_stored(:)
 character(len=500) :: message

!******************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

!Initialize MPI
 call xmpi_init()
 comm = xmpi_world
 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 ABI_CHECK(nprocs == 1, "cut3d not programmed for parallel execution")

!Initialize memory profiling if it is activated
!if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 call timein(tcpui,twalli)

!Default for sequential use
!Other values of mpi_enreg are dataset dependent, and should NOT be initialized inside cut3d.F90.
 call initmpi_seq(mpi_enreg)

 codename='CUT3D '//repeat(' ',18)
 call herald(codename,abinit_version,std_out)

!BIG LOOP on files
 ABI_MALLOC(isdenpot,(mfiles))
 isdenpot=0
 ABI_MALLOC(filrho_stored,(mfiles))
 iomode = IO_MODE_FORTRAN

 do ifiles=1,mfiles

!  Get name of density file
   write(std_out,*)
   write(std_out,*) ' What is the name of the 3D function (density, potential or wavef) file ?'
   if (read_string(filrho, unit=std_in) /= 0) then
     MSG_ERROR("Fatal error!")
   end if
   filrho_tmp=adjustl(filrho)
   do ii=1,len_trim(filrho_tmp)
     if(filrho_tmp(ii:ii)==blank)then
       filrho=trim(filrho_tmp(1:ii-1))
       exit
     end if
   end do
   write(std_out,*) ' => Your 3D function file is: ',trim(filrho)
   write(std_out,*)
   ! Checking the existence of data file
   if (nctk_try_fort_or_ncfile(filrho, message) /= 0) then
     MSG_ERROR(message)
   end if

!  Treat the different cases: formatted or unformatted
   iomode = IO_MODE_FORTRAN; if (endswith(filrho, ".nc")) iomode = IO_MODE_ETSF
   if (iomode == IO_MODE_FORTRAN) then
     write(std_out,"(a)") '- Your file contains unformatted binary header + 3D data'
   else
     write(std_out,"(a)") '- Your file contains ETSF data'
   end if

   ! Read the header and extract dimensions.
   write(std_out,*)
   call hdr_read_from_fname(hdr, filrho, fform0, comm)
   ABI_CHECK(fform0 /= 0, "hdr_read returned fform = 0")
   abifile = abifile_from_fform(fform0)
   ABI_CHECK(abifile%fform /= 0, "Cannot detect abifile from fform")

!  Echo part of the header
   call hdr%echo(fform0, 4)

   nr1=hdr%ngfft(1); nr2=hdr%ngfft(2); nr3=hdr%ngfft(3)
   natom=hdr%natom
   nspden=hdr%nspden
   ntypat=hdr%ntypat
   rprimd(:,:)=hdr%rprimd(:,:)

!  Need to know natom in order to allocate xcart
   ABI_MALLOC(xcart,(3,natom))
   ABI_MALLOC(xred,(3,natom))
   xred(:,:)=hdr%xred(:,:)
   call xred2xcart(natom,rprimd,xcart,xred)

   ispden=0
   if (abifile%class == "density" .or. abifile%class == "potential") then
     if(nspden/=1)then
       write(std_out,'(a)' )' '
       write(std_out,'(a)' )' * This file contains more than one spin component,'
       write(std_out,'(a,i3,a)' )'  (indeed, nspden=',nspden,' )'
       write(std_out,'(a)' )'  Some of the tasks that you will define later will concern all spin components.'
       write(std_out,'(a)' )'  Others tasks might require you to have chosen among the following:'
     end if
     if(nspden==2)then
       write(std_out,'(a)' )'   ispden= 0 ==> Total density'
       write(std_out,'(a)' )'   ispden= 1 ==> spin-up density'
       write(std_out,'(a)' )'   ispden= 2 ==> spin-down density'
       write(std_out,'(a)' )'   ispden= 3 ==> spin-polarization (or magnetization) density'
       write(std_out,'(a)' )'                 spin up - spin down difference.'
     end if
     if(nspden==4)then
       write(std_out,'(a)' )'   ispden= 0 ==> Total density'
       write(std_out,'(a)' )'   ispden= 1 ==> magnetization in the x direction'
       write(std_out,'(a)' )'   ispden= 2 ==> magnetization in the y direction'
       write(std_out,'(a)' )'   ispden= 3 ==> magnetization in the z direction'
       write(std_out,'(a)' )'   ispden= 4 might be used to plot the magnetization (3D) in the XCrysDen format,'
     end if
     if(nspden/=1)then
       write(std_out,*)'  Please define ispden:'
       read(std_in,*)ispden
       write(std_out,'(a,i3)' )' You entered ispden=',ispden
     end if
   end if

   write(std_out,*)
   write(std_out,*) '==========================================================='
   write(std_out,*)

!  Echo the value of different input parameters
   write(std_out,*)'ECHO important input variables ...'
   write(std_out,*)
   write(std_out,*) ' Dimensional primitive vectors (ABINIT equivalent: rprimd):'
   write(std_out,'(3es16.6)' ) rprimd(1:3,1)
   write(std_out,'(3es16.6)' ) rprimd(1:3,2)
   write(std_out,'(3es16.6)' ) rprimd(1:3,3)

!  Compute ucvol and test the non-collinearity of rprimd vectors.
   call metric(gmet,gprimd,dev_null,rmet,rprimd,ucvol)

   write(std_out,'(a,3i5)' ) '  Grid density (ABINIT equivalent: ngfft): ',nr1,nr2,nr3
   write(std_out,*) ' Number of atoms       :',natom
   write(std_out,*) ' Number of atomic types:',ntypat

   write(std_out,*)
   write(std_out,*) '  #    Atomic positions (cartesian coordinates - Bohr)'
   do iatom=1,natom
     write(std_out,'(i4,3es16.6)' )iatom,xcart(1:3,iatom)
   end do
   write(std_out,*)

!  ------------------------------------------------------------------------
!  Branching: either WF file, or DEN/POT file.

   if (abifile%class == "wf_planewave") then
     write(std_out,*)' This file is a WF file. '
     isdenpot(ifiles)=0
     iprompt = 0 ! this needs to be initialized, as it is used after the loop on files...

     call cut3d_wffile(filrho,hdr%ecut_eff,exchn2n3d0,hdr%istwfk,hdr%kptns,natom,hdr%nband,hdr%nkpt,hdr%npwarr,&
&     nr1,nr2,nr3,hdr%nspinor,hdr%nsppol,ntypat,rprimd,xcart,hdr%typat,hdr%znucltypat)
     call hdr%free()

!    -------------------------------------------------------------------------
! This is a DEN/POT file
   else if (abifile%class == "density" .or. abifile%class == "potential") then

!    This should become a subroutine
     write(std_out,*)' This file is a Density or Potential file '
     isdenpot(ifiles)=1

!    Read the function on the 3D grid
     ABI_MALLOC(grid,(nr1,nr2,nr3))
     ABI_MALLOC(grid_full,(nr1,nr2,nr3,nspden))
     ABI_MALLOC(gridtt,(nr1,nr2,nr3))
     ABI_MALLOC(gridmx,(nr1,nr2,nr3))
     ABI_MALLOC(gridmy,(nr1,nr2,nr3))
     ABI_MALLOC(gridmz,(nr1,nr2,nr3))

     varname = varname_from_fname(filrho)
     if (iomode == IO_MODE_ETSF) then
       call wrtout(std_out, sjoin("- Reading netcdf variable: ", varname))
     end if

     call cut3d_rrho(filrho,varname,iomode,grid_full,nr1,nr2,nr3,nspden)

!    Do not forget that the first sub-array of a density file is the total density,
!    while the first sub-array of a potential file is the spin-up potential
     if (abifile%class == "density") then

!      gridtt= grid --> Total density or potential.
!      gridmx= grid --> spin-Up density, or magnetization density in X direction.
!      gridmy= grid --> spin-Down density, or magnetization density in Y direction.
!      gridmz= grid --> spin-polarization density (Magnetization),
!      or magnetization density in Z direction.
       gridtt(:,:,:)=grid_full(:,:,:,1)
       if(nspden==2)then
         gridmx = grid_full(:,:,:,2)
         gridmy = grid_full(:,:,:,1)-grid_full(:,:,:,2)
         gridmz = -grid_full(:,:,:,1)+two*grid_full(:,:,:,2)
       else if(nspden==4)then
         gridmx = grid_full(:,:,:,2)
         gridmy = grid_full(:,:,:,3)
         gridmz = grid_full(:,:,:,4)
       end if

       if(nspden==1)then
         grid = grid_full(:,:,:,1)
       else
         if(ispden==0)then
           grid = gridtt
         else if(ispden==1)then
           grid = gridmx
         else if(ispden==2)then
           grid = gridmy
         else if(ispden==3)then
           grid = gridmz
!          if(ispden==0)then
!          grid(:,:,:)=grid_full(:,:,:,1)
!          else if(ispden==1)then
!          grid(:,:,:)=grid_full(:,:,:,2)
!          else if(ispden==2)then
!          grid(:,:,:)=grid_full(:,:,:,1)-grid_full(:,:,:,2)
!          else if(ispden==-1)then
!          grid(:,:,:)=-grid_full(:,:,:,1)+two*grid_full(:,:,:,2)
         else if(ispden==4)then
           write(std_out,*) ' '
         else
           MSG_ERROR(sjoin('bad ispden value = ',itoa(ispden)))
         end if
       end if

     else if (abifile%class == "potential") then   ! Potential case
       if(ispden==0)then
         grid(:,:,:)=grid_full(:,:,:,1)
       else if(ispden==1 .or. ispden==2)then
         grid(:,:,:)=grid_full(:,:,:,ispden)
       else
         MSG_ERROR(sjoin('bad ispden value = ',itoa(ispden)))
       end if
       gridtt = grid
     end if

     write(std_out,*)
     write(std_out,*) ' 3D function was read. Ready for further treatment.'
     write(std_out,*)
     write(std_out,*) '==========================================================='
     write(std_out,*)

!    ------------------------------------------------------------------------

!    At this moment all the input is done
!    The code knows the geometry of the system,
!    and the data file (electron density, potential, etc).
!    It will further calculate the electron density by interpolation in
!    a point, along a line or in a plane.

     do
       do
         write(std_out,*) ' What is your choice ? Type:'
         write(std_out,*) '  0 => exit'
         write(std_out,*) '  1 => point  (interpolation of data for a single point)'
         write(std_out,*) '  2 => line   (interpolation of data along a line)'
         write(std_out,*) '  3 => plane  (interpolation of data in a plane)'
         write(std_out,*) '  4 => volume (interpolation of data in a volume)'
         write(std_out,*) '  5 => 3D formatted data (output the bare 3D data - one column)'
         write(std_out,*) '  6 => 3D indexed data (bare 3D data, preceeded by 3D index)'
         write(std_out,*) '  7 => 3D Molekel formatted data '
         write(std_out,*) '  8 => 3D data with coordinates (tecplot ASCII format)'
         write(std_out,*) '  9 => output .xsf file for XCrysDen'
         write(std_out,*) ' 11 => compute atomic charge using the Hirshfeld method'
         write(std_out,*) ' 14 => Gaussian/cube wavefunction module'
         write(std_out,*) ' 15 => Write data to netcdf file'
         read(std_in,*) itask
         write(std_out,'(a,a,i2,a)' ) ch10,' Your choice is ',itask,ch10

         if ((5 <= itask .and. itask <= 9) .or. any(itask == [14, 15]) )then
           write(std_out,*) ch10,'  Enter the name of an output file:'
           if (read_string(filnam, unit=std_in) /= 0) then
             MSG_ERROR("Fatal error!")
           end if
           write(std_out,*) '  The name of your file is: ',trim(filnam)
         end if

         select case(itask)

         case(1) ! point calculation
           call cut3d_pointint(gridtt,gridmx,gridmy,gridmz,nr1,nr2,nr3,nspden,rprimd)
           exit

         case(2) ! line calculation
           call cut3d_lineint(gridtt,gridmx,gridmy,gridmz,nr1,nr2,nr3,nspden,rprimd)
           exit

         case(3) ! plane calculation
           call cut3d_planeint(gridtt,gridmx,gridmy,gridmz,natom,nr1,nr2,nr3,nspden,rprimd,xcart)
           exit

         case(4) ! volume calculation
           write(std_out,*) ' Enter volume calculation'
           call cut3d_volumeint(gridtt,gridmx,gridmy,gridmz,natom,nr1,nr2,nr3,nspden,rprimd,xcart)
           exit

         case(5)
           ! Rewrite the data on a formatted file, just in one (or four) column(s)
           if (open_file(filnam,message, newunit=unt, status='unknown', action="write") /= 0) then
             MSG_ERROR(message)
           end if

           if(nspden==1)then
             do i3=1,nr3
               do i2=1,nr2
                 do i1=1,nr1
                   write(unt,'(4(es22.12))') grid(i1,i2,i3)
                 end do
               end do
             end do
           else
             do i3=1,nr3
               do i2=1,nr2
                 do i1=1,nr1
                   write(unt,'(4(es22.12))') gridtt(i1,i2,i3), gridmx(i1,i2,i3), gridmy(i1,i2,i3), gridmz(i1,i2,i3)
                 end do
               end do
             end do
           end if
           close(unt)
           exit

         case(6)
!            Rewrite the data on a formatted file, 3D index + density
           if (open_file(filnam,message, newunit=unt, status='unknown', action="write") /= 0) then
             MSG_ERROR(message)
           end if

           if(nspden==1)then
             write(unt,*)'   i1    i2    i3      data '
             do i3=1,nr3
               do i2=1,nr2
                 do i1=1,nr1
                   write(unt,'(3i6,4(es24.14))') i1,i2,i3,grid(i1,i2,i3)
                 end do
               end do
             end do
           else
             if(nspden==2)then
               write(unt,*)'   i1    i2    i3     non-spin-polarized spin up  spin down  difference  '
             else if(nspden==4)then
               write(unt,*)'   i1    i2    i3     non-spin-polarized   x       y      z   '
             end if
             do i3=1,nr3
               do i2=1,nr2
                 do i1=1,nr1
                   write(unt,'(3i6,4(es24.14))') i1,i2,i3,gridtt(i1,i2,i3),gridmx(i1,i2,i3),gridmy(i1,i2,i3),gridmz(i1,i2,i3)
                 end do
               end do
             end do
           end if ! nspden
           close(unt)
           exit

         case(7)
           if (open_file(filnam,message, newunit=unt, form='unformatted', action="write") /= 0) then
             MSG_ERROR(message)
           end if

           xm=0 ; xp=rprimd(1,1)*Bohr_Ang
           ym=0 ; yp=rprimd(2,2)*Bohr_Ang
           zm=0 ; zp=rprimd(3,3)*Bohr_Ang
           write(std_out,'(/,a,/)' )' Extremas (x,y,z) of the cube in which the molecule is placed, in Angstroms'
           write(std_out,'(5x,6f10.5)' ) xm,xp,ym,yp,zm,zp
           write(std_out,'(/,a,2x,3i5)' )' Number of points per side: ',nr1,nr2,nr3
           write(std_out,'(/,a,2x,i10,//)' )' Total number of points:', nr1*nr2*nr3
           write(unt) xm,xp,ym,yp,zm,zp,nr1,nr2,nr3
           ABI_MALLOC(rhomacu,(nr1,nr2))
           do i3=1,nr3
             do i2=1,nr2
               do i1=1,nr1
                 rhomacu(i1,i2)=grid(i1,i2,i3)
               end do
             end do
             write(unt) rhomacu(:,:)
           end do
           close(unt)
           exit

         case (8)
           if (open_file(filnam, message, newunit=unt, form='formatted', action="write") /= 0) then
             MSG_ERROR(message)
           end if

           write(std_out,'(/,a,/)' )' Extremas (x,y,z) of the cube in which the molecule is placed, in Angstroms'
           write(std_out,'(5x,6f10.5)' ) xm,xp,ym,yp,zm,zp
           write(std_out,'(/,a,2x,3i5)' )' Number of points per side: ',nr1,nr2,nr3
           write(std_out,'(/,a,2x,i10,//)' )' Total number of points:', nr1*nr2*nr3
           write(unt,'(a)') 'TITLE = "  " '
           write(unt,'(a)') 'VARIABLES = "X"  "Y"  "Z" (all three in Angstrom)  "DENSITY or POTENTIAL" (atomic units) '
           write(unt,'(3(a,i6),a)') 'ZONE I=',nr1, ' J=', nr2, ' K=', nr3, ' F=POINT'
           do i3=1,nr3
             do i2=1,nr2
               do i1=1,nr1
                 xnow = rprimd(1,1)*(i1-1)/nr1 + rprimd(1,2)*(i2-1)/nr2 + rprimd(1,3)*(i3-1)/nr3
                 ynow = rprimd(2,1)*(i1-1)/nr1 + rprimd(2,2)*(i2-1)/nr2 + rprimd(2,3)*(i3-1)/nr3
                 znow = rprimd(3,1)*(i1-1)/nr1 + rprimd(3,2)*(i2-1)/nr2 + rprimd(3,3)*(i3-1)/nr3
                 write(unt,'(4es22.15)') Bohr_Ang*xnow, Bohr_Ang*ynow, Bohr_Ang*znow, grid (i1,i2,i3)
               end do
             end do
           end do
           close(unt)
           exit

         case (9)
           if (open_file(filnam, message, newunit=unt, form='formatted', action="write") /= 0) then
             MSG_ERROR(message)
           end if
           xm=0 ; xp=rprimd(1,1)*Bohr_Ang
           ym=0 ; yp=rprimd(2,2)*Bohr_Ang
           zm=0 ; zp=rprimd(3,3)*Bohr_Ang
           write(std_out,'(/,a,/)' )' Extremas (x,y,z) of the cube in which the molecule is placed, in Angstroms'
           write(std_out,'(5x,6f10.5)' ) xm,xp,ym,yp,zm,zp
           write(std_out,'(/,a,2x,3i5)' )' Number of points per side: ',nr1+1,nr2+1,nr3+1
           write(std_out,'(/,a,2x,i10,//)' )' Total number of points:', (nr1+1)*(nr2+1)*(nr3+1)
           write(std_out,*) '  znucl = ', hdr%znucltypat, ' type = ', hdr%typat, ' ntypat = ', ntypat

           gridshift1 = 0
           gridshift2 = 0
           gridshift3 = 0
           write(std_out,*) 'Do you want to shift the grid along the x,y or z axis (y/n)?'
           write(std_out,*)
           shift_tau(:) = zero
           read(std_in,"(a)") outputchar
           if (outputchar == 'y' .or. outputchar == 'Y') then
             write(std_out,*) 'Give the three shifts (x,y,z < ',nr1,nr2,nr3,'):'
             write(std_out,*)
             read(std_in,*) gridshift1, gridshift2, gridshift3
             shift_tau(:) = gridshift1*rprimd(:,1)/(nr1+1) + gridshift2*rprimd(:,2)/(nr2+1) + gridshift3*rprimd(:,3)/(nr3+1)
           end if
!
!            Generate translated coordinates to match density shift
!
           ABI_MALLOC(tau2,(3,natom))
           do iatom = 1,natom
             tau2(:,iatom) = xcart(:,iatom) - shift_tau(:)
           end do
!            ################################################################### (LD)
!            Option only available for "xcrysden" format as documented at the beginning
           if (ispden==4) then
!              It is necessary to know previously how many atoms will be used.
!              in order to plot the necessary magnetization arrows only.
             write(std_out,*)'Is it possible to decrease the number of arrows in order to improve the'
             write(std_out,*)'visualization in the screen, and decrease the size of the xcrysden output file.'
             write(std_out,*)'How many arrows would you like to skip? 0 = take all. 1 = skip every other point...'
             read (std_in,*) nrws
             nrws=nrws+1
             index=natom
             maxmz=0.0
             do i1=1,nr1,nrws
               do i2=1,nr2,nrws
                 do i3=1,nr3,nrws
                   normz=gridmx(i1,i2,i3)**2+gridmy(i1,i2,i3)**2+gridmz(i1,i2,i3)**2
                   if(normz > maxmz) maxmz=normz
                 end do
               end do
             end do
             if(abs(maxmz)<tol10)then
               MSG_ERROR('At least, one of the components must differ from zero.')
             end if
             do i1=1,nr1,nrws
               do i2=1,nr2,nrws
                 do i3=1,nr3,nrws
                   normz=gridmx(i1,i2,i3)**2+gridmy(i1,i2,i3)**2+gridmz(i1,i2,i3)**2
                   if(0.1*maxmz <= normz) index=index+1
                 end do
               end do
             end do

             write(unt,'(1X,A)') 'CRYSTAL'
             write(unt,'(1X,A)') 'PRIMVEC'
             do i1 = 1,3
               write(unt,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(i2,i1), i2=1,3)
             end do
             write(unt,'(1X,A)') 'PRIMCOORD'
             write(unt,*) index, '1'

             ! write out atom types and positions
             do iatom = 1,natom
               write(unt,'(i9,3(3X,ES17.10))') nint(hdr%znucltypat(hdr%typat(iatom))),Bohr_Ang*tau2(1:3,iatom)
             end do

             ! write out magnetization vectors.
             ! xcrysden consider these as X (dummy) atoms.
             do i1=1,nr1,nrws
               do i2=1,nr2,nrws
                 do i3=1,nr3,nrws
                   normz=gridmx(i1,i2,i3)**2+gridmy(i1,i2,i3)**2+gridmz(i1,i2,i3)**2
                   if(0.1*maxmz <= normz) then
                     xcart2 = matmul (rprimd, (/(i1-one)/nr1, (i2-one)/nr2, (i3-one)/nr3/))
                     write(unt,'(A,1X,6(ES17.10,2X))')'X',&
                     Bohr_Ang*(xcart2(1)-shift_tau(1)),&
                     Bohr_Ang*(xcart2(2)-shift_tau(2)),&
                     Bohr_Ang*(xcart2(3)-shift_tau(3)),&
                     gridmx(i1,i2,i3),&
                     gridmy(i1,i2,i3),&
                     gridmz(i1,i2,i3)
                   end if
                 end do
               end do
             end do
           else
!              ################################################################### (LD)
!
!              normal case: output density or potential (scalar field)
             write(unt,'(1X,A)')  'DIM-GROUP'
             write(unt,*) '3  1'
             write(unt,'(1X,A)') 'PRIMVEC'
             do i1 = 1,3
               write(unt,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(i2,i1), i2=1,3)
             end do
             write(unt,'(1X,A)') 'PRIMCOORD'
             write(unt,*) natom, ' 1'
             do iatom = 1,natom
               write(unt,'(i9,3(3X,ES17.10))') nint(hdr%znucltypat(hdr%typat(iatom))),Bohr_Ang*tau2(1:3,iatom)
             end do
             write(unt,'(1X,A)') 'ATOMS'
             do iatom = 1,natom
               write(unt,'(i9,3(3X,ES17.10))') nint(hdr%znucltypat(hdr%typat(iatom))),Bohr_Ang*tau2(1:3,iatom)
             end do
!              write(31,'(1X,A)') 'FRAMES'
             write(unt,'(1X,A)') 'BEGIN_BLOCK_DATAGRID3D'
             write(unt,*) 'datagrids'
             write(unt,'(1X,A)') 'DATAGRID_3D_DENSITY'
             write(unt,*) nr1+1,nr2+1,nr3+1
             write(unt,*) '0.0 0.0 0.0 '
             do i1 = 1,3
               write(unt,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(i2,i1), i2=1,3)
             end do

             index = 0
             do ir3=gridshift3+1,nr3+1
               ii3=mod(ir3-1,nr3) + 1
               do ir2=gridshift2+1,nr2+1
                 ii2=mod(ir2-1,nr2) + 1
                 do ir1=gridshift1+1,nr1+1
                   ii1=mod(ir1-1,nr1) + 1
                   write(unt,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                   index = index+1
                   if (mod (index,6) == 0) write (unt,*)
                 end do
                 do ir1=1,gridshift1
                   ii1=mod(ir1-1,nr1) + 1
                   write(unt,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                   index = index+1
                   if (mod (index,6) == 0) write (unt,*)
                 end do
               end do
               do ir2=1,gridshift2
                 ii2=mod(ir2-1,nr2) + 1
                 do ir1=gridshift1+1,nr1+1
                   ii1=mod(ir1-1,nr1) + 1
                   write(unt,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                   index = index+1
                   if (mod (index,6) == 0) write (unt,*)
                 end do
                 do ir1=1,gridshift1
                   ii1=mod(ir1-1,nr1) + 1
                   write(unt,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                   index = index+1
                   if (mod (index,6) == 0) write (unt,*)
                 end do
               end do
             end do
             do ir3=1,gridshift3
               ii3=mod(ir3-1,nr3) + 1
               do ir2=gridshift2+1,nr2+1
                 ii2=mod(ir2-1,nr2) + 1
                 do ir1=gridshift1+1,nr1+1
                   ii1=mod(ir1-1,nr1) + 1
                   write(unt,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                   index = index+1
                   if (mod (index,6) == 0) write (unt,*)
                 end do
                 do ir1=1,gridshift1
                   ii1=mod(ir1-1,nr1) + 1
                   write(unt,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                   index = index+1
                   if (mod (index,6) == 0) write (unt,*)
                 end do
               end do
               do ir2=1,gridshift2
                 ii2=mod(ir2-1,nr2) + 1
                 do ir1=gridshift1+1,nr1+1
                   ii1=mod(ir1-1,nr1) + 1
                   write(unt,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                   index = index+1
                   if (mod (index,6) == 0) write (unt,*)
                 end do
                 do ir1=1,gridshift1
                   ii1=mod(ir1-1,nr1) + 1
                   write(unt,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                   index = index+1
                   if (mod (index,6) == 0) write (unt,*)
                 end do
               end do
             end do
             write (unt,*)
             write(unt,'(1X,A)') 'END_DATAGRID_3D'
             write(unt,'(1X,A)') 'END_BLOCK_DATAGRID3D'

           end if

           close(unt)
           exit

         case(11)
           call cut3d_hirsh(grid,natom,nr1,nr2,nr3,ntypat,rprimd,xcart,hdr%typat,hdr%zionpsp,hdr%znucltypat)
           exit

         case(14) ! CUBE file format from GAUSSIAN
           write(std_out,*)
           write(std_out,*) 'Output a cube file of 3D volumetric data'
           write(std_out,*)

!            EXAMPLE FROM THE WEB
!            CPMD CUBE FILE.
!            OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z
!            3    0.000000    0.000000    0.000000
!            40    0.283459    0.000000    0.000000
!            40    0.000000    0.283459    0.000000
!            40    0.000000    0.000000    0.283459
!            8    0.000000    5.570575    5.669178    5.593517
!            1    0.000000    5.562867    5.669178    7.428055
!            1    0.000000    7.340606    5.669178    5.111259
!            -0.25568E-04  0.59213E-05  0.81068E-05  0.10868E-04  0.11313E-04  0.35999E-05

           if (open_file(filnam,message,newunit=unt, status='unknown', form='formatted', action="write") /= 0) then
             MSG_ERROR(message)
           end if

           !%% call print_fofr_cube(nr1,nr2,n3,nr1,nr2,nr3,fofr,rprimd,natom,znucl_atom,xcart,unit=unt)
           write(unt,'(a)') 'ABINIT generated cube file'
           write(unt,'(a)') 'from cut3d tool'

           write(unt,'(i9,3(1x,f12.6))') natom,0.,0.,0.
           write(unt,'(i9,3(1x,f12.6))') nr1,(rprimd(ir2,1)/nr1, ir2=1,3)
           write(unt,'(i9,3(1x,f12.6))') nr2,(rprimd(ir2,2)/nr2, ir2=1,3)
           write(unt,'(i9,3(1x,f12.6))') nr3,(rprimd(ir2,3)/nr3, ir2=1,3)

           do iatom=1,natom
             write(unt,'(i9,4(3X,ES17.10))') nint(hdr%znucltypat(hdr%typat(iatom))),0.d0, &
&             xcart(1,iatom),xcart(2,iatom),xcart(3,iatom)
           end do

!            C ordering of the indexes
           do i1=1,nr1
             do i2=1,nr2
               do i3=1,nr3
                 write(unt,'(6(f12.6,2x))') grid(i1,i2,i3)
               end do
             end do
           end do

           close(unt)
           exit

         case (15)
           ! Write netcdf file.
           timrev = 2; if (any(hdr%kptopt == [3, 4])) timrev = 1
           cryst = hdr%get_crystal(timrev)
           call ngfft_seq(ngfft, [nr1, nr2, nr3])
           ngfft(4:6) = ngfft(1:3)
           nfft = product(ngfft(1:3))
           cplex = 1
           call init_distribfft_seq(mpi_enreg%distribfft, 'c', ngfft(2), ngfft(3), 'all')
           call init_distribfft_seq(mpi_enreg%distribfft, 'f', ngfft(2), ngfft(3), 'all')

           call fftdatar_write(varname,filnam,IO_MODE_ETSF,hdr,cryst,ngfft,cplex,nfft,nspden,grid_full,mpi_enreg)
           call cryst%free()

         case(0)
           write(std_out,*)' Exit requested by user'
           exit

         case default
           MSG_ERROR(sjoin("Wrong task:", itoa(itask)))
         end select
       end do

       write(std_out,*) ' Task ',itask,' has been done !'
       write(std_out,*)
       write(std_out,'(a)') ' More analysis of the 3D file ? ( 0=no ; 1=default=yes ; 2= treat another file - restricted usage)'
       read(std_in,*) iprompt
       if(iprompt/=1) then
         call hdr%free()
         exit
       else
         cycle
       end if
     end do

   else
     MSG_ERROR(sjoin("Don't know how to handle file class ", abifile%class))
   end if ! WF file or DEN/POT file

!  A maximum number of files had been previously specified, but set the actual number of files
!  to 1 if one does not read at least one other.
   if(ifiles==1)then
     nfiles=1
     if(iprompt==2)nfiles=mfiles

!    A data structure for storing the important information should be created ...
!    Here, one supposes that the files are compatible ...
     if(isdenpot(ifiles)==1)then
       ABI_MALLOC(grid_full_stored,(nr1,nr2,nr3,nspden,nfiles))
       nr1_stored=nr1
       nr2_stored=nr2
       nr3_stored=nr3
       nspden_stored=nspden
     else if(isdenpot(ifiles)/=1 .and. iprompt==2)then
       MSG_ERROR("in case of storage mode, the first file must be a density/potential file.")
     end if
   end if

   if(isdenpot(ifiles)==1) grid_full_stored(:,:,:,:,ifiles)=grid_full(:,:,:,:)
   if(isdenpot(ifiles)==1) filrho_stored(ifiles)=filrho

   if(allocated(xcart)) then
     ABI_FREE(xcart)
   end if
   if(allocated(xred)) then
     ABI_FREE(xred)
   end if
   if(allocated(grid)) then
     ABI_FREE(grid)
   end if
   if(allocated(grid_full)) then
     ABI_FREE(grid_full)
   end if
   if(allocated(gridtt)) then
     ABI_FREE(gridtt)
   end if
   if(allocated(gridmx)) then
     ABI_FREE(gridmx)
   end if
   if(allocated(gridmy)) then
     ABI_FREE(gridmy)
   end if
   if(allocated(gridmz)) then
     ABI_FREE(gridmz)
   end if
   if(allocated(rhomacu)) then
     ABI_FREE(rhomacu)
   end if
   if(allocated(tau2)) then
     ABI_FREE(tau2)
   end if

   if(iprompt/=2) then
     exit
   end if

 end do ! End big loop on files

!Will provide different information on the density and potential files
 do ifiles=1,nfiles
   if(isdenpot(ifiles)==1)then
     write(std_out,*)
     write(std_out,*) ' Provide some global information about the density and/or potential file(s)'
     exit
   end if
 end do
 do ifiles=1,nfiles
   if(isdenpot(ifiles)==1)then
     write(std_out,*)
     write(std_out, '(a,i5,3a)' ) '-  File number ',ifiles,', with name "',trim(filrho_stored(ifiles)),'"'
     write(std_out, '(a,i12,a,es14.6)' ) '  Number of grid points =',nr1*nr2*nr3,' ; Volume of real space cell (Bohr^3)=',ucvol
     do ispden=1,nspden
       sumdenpot=sum(grid_full_stored(:,:,:,ispden,ifiles))
       write(std_out, '(a,i5,3a)' ) '   Spin-component number ',ispden
       write(std_out, '(a,3es16.6)' ) '      Sum of values, mean, mean times cell volume=',&
&       sumdenpot,sumdenpot/real(nr1*nr2*nr3),sumdenpot*ucvol/real(nr1*nr2*nr3)
     end do
   end if
 end do

 if(nspden==1)then
!  At present, only nspden=1 is correctly implemented, due to specificities of the treatment of the spin-density
   do ifiles=1,nfiles
     if(isdenpot(ifiles)==1)then
       write(std_out,*)
       write(std_out,'(a)') ' Provide some global joint information about the stored density and potential file(s)'
       exit
     end if
   end do
   do ifiles=1,nfiles
     if(isdenpot(ifiles)==1)then
       do jfiles=ifiles,nfiles
         if(isdenpot(jfiles)==1)then
           write(std_out,*)
           write(std_out, '(a,2i5)' )'  File numbers: ',ifiles,jfiles
           do ispden=1,nspden
             dotdenpot=zero
             do ir1=1,nr1
               do ir2=1,nr2
                 do ir3=1,nr3
                   dotdenpot=dotdenpot+grid_full_stored(ir1,ir2,ir3,ispden,ifiles)*grid_full_stored(ir1,ir2,ir3,ispden,jfiles)
                 end do
               end do
             end do
             write(std_out, '(a,i5,3a)' ) '   Spin-component number ',ispden
             write(std_out, '(a,3es16.6)' ) '      Dot product of values, mean, mean times cell volume=',&
!            write(std_out, '(a,3es20.10)' ) '      Dot product of values, mean, mean times cell volume=',&
&             dotdenpot,dotdenpot/real(nr1*nr2*nr3),dotdenpot*ucvol/real(nr1*nr2*nr3)
           end do
         end if
       end do
     end if
   end do
 end if

 ABI_FREE(filrho_stored)

 if(allocated(grid_full_stored)) then
   ABI_FREE(grid_full_stored)
 end if

 if(allocated(isdenpot)) then
   ABI_FREE(isdenpot)
 end if

 call timein(tsec(1),tsec(2))
 tsec(1)=tsec(1)-tcpui
 tsec(2)=tsec(2)-twalli

 write(std_out, '(3a,f13.1,a,f13.1)' )'-',ch10,'- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)

 write(std_out,*)
 write(std_out,*) ' Thank you for using me'
 write(std_out,*)

 call flush_unit(std_out)
 call destroy_mpi_enreg(mpi_enreg)
 call abinit_doctor("__cut3d")
 call xmpi_end()

 end program cut3d
!!***
