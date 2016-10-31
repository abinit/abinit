!{\src2tex{textfont=tt}}
!!****p* ABINIT/optic
!! NAME
!! optic
!!
!! FUNCTION
!! Driver routine to call linopt and nlinopt, which calculate
!! the linear and non-linear optical responses in the RPA.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (SSharma,MVer,VRecoules,YG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!!  bantot
!!  doccde(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy.
!!  domega=frequency range
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree).
!!  eigen11(2*mband*mband*nkpt_rbz*nsppol)=first-order eigenvalues (hartree)
!!  in reciprocal direction 100
!!  eigen12(2*mband*mband*nkpt_rbz*nsppol)=first-order eigenvalues (hartree)
!!  in reciprocal direction 010
!!  eigen13(2*mband*mband*nkpt_rbz*nsppol)=first-order eigenvalues (hartree)
!!  in reciprocal direction 001
!!  ecut=kinetic energy planewave cutoff (hartree).
!!  entropy= entropy associated with the smearing (adimensional)
!!  fermie= fermi energy (Hartree)
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{2}$).
!!  gmet_inv(3,3)=inverse of reciprocal space metric ($\textrm{bohr}^{2}$).
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  nomega=number of frequency for conductivity computation
!!  mband=maximum number of bands.
!!  natom = number of atoms in the unit cell.
!!  nband(nkpt*nsppol)=number of bands at each RF k point for each spin.
!!  nelect=number of electrons per unit cell
!!  nkpt=number of k points in the IBZ for this perturbation
!!  ngfft(3)=integer fft box dimensions.
!!  nspinor=number of spinorial components of the wavefunctions.
!!  nsppol=1 for unpolarized, 2 for spin-polarized.
!!  ntypat = number of atom types.
!!  occ(mband*nkpt*nsppol)=occupation number for each band and k.
!!  occopt==option for occupancies
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).
!!  rprimd(3,3)=real space primitive translations.
!!  of primitive translations.
!!  broadening=smearing width (or temperature) in Hartree
!!  ucvol=unit cell volume in ($\textrm{bohr}^{3}$).
!!  maxomega=frequency windows for computations of sigma
!!  wtk(nkpt)=weight assigned to each k point.
!!  znucl(natom)=atomic number of atoms
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,abinit_doctor,ebands_free,ebands_init
!!      ebands_update_occ,flush_unit,hdr_bcast,hdr_copy,hdr_free,herald
!!      int2char4,linelop,linopt,mati3inv,matr3inv,metric,nctk_fort_or_ncfile
!!      nlinopt,nonlinopt,pmat2cart,pmat_renorm,sym2cart,timein,wfk_close
!!      wfk_open_read,wfk_read_eigk,wrtout,xmpi_bcast,xmpi_end,xmpi_init
!!      xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program optic

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_xmpi
 use m_xomp
 use m_profiling_abi
 use m_build_info
 use m_optic_tools
 use m_wfk
 use m_nctk
 use m_hdr
 use m_ebands
 use m_eprenorms
 use m_crystal

 use m_time ,     only : asctime
 use m_io_tools,  only : flush_unit, open_file, file_exists, num_opened_units, show_units
 use m_specialmsg,only : specialmsg_getcount
 use m_fstrings,  only : int2char4, itoa, sjoin, strcat

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'optic'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
 integer,parameter :: formeig0=0,formeig1=1,tim_rwwf=0,master=0
 integer :: bantot,bdtot0_index,bdtot_index!,dosdeltae
 integer :: headform,ierr,ii,ikpt,isym
 integer :: isppol,mband,nomega,natom,nband1
 integer :: nsym
 integer :: nkpt,nspinor,nsppol,ntypat
 integer :: occopt,nks_per_proc,work_size,lin1,lin2,nlin1,nlin2,nlin3
 integer :: linel1,linel2,linel3,nonlin1,nonlin2,nonlin3
 integer :: num_lin_comp=1,num_nonlin_comp=0,num_linel_comp=0,num_nonlin2_comp=0
 integer :: autoparal=0,max_ncpus=0
 integer :: nonlin_comp(27) = 0
 integer :: linel_comp(27) = 0
 integer :: nonlin2_comp(27) = 0
 integer :: lin_comp(9) = [11, 22 ,33, 12, 13, 21, 23, 31, 32]
 integer,allocatable :: nband(:)
 integer,allocatable :: symrel(:,:,:)
 integer,allocatable :: symrec(:,:,:)
 real(dp),allocatable :: symcart(:,:,:)
 real(dp) :: domega,ecut,fermie!,maxocc,entropy
 real(dp):: eff
 logical :: do_antiresonant, do_temperature
 integer,allocatable :: istwfk(:), npwarr(:)
 real(dp) :: nelect
 real(dp) :: broadening,ucvol,maxomega,scissor,tolerance,tphysel
 real(dp) :: tcpu,tcpui,twall,twalli
 real(dp) :: tsec(2)
 real(dp) :: gmet(3,3),gmet_inv(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),gprimd_trans(3,3)
 real(dp),allocatable :: kpt(:,:)
 real(dp),allocatable :: cond_kg(:),cond_nd(:),doccde(:)
 real(dp),allocatable :: eig0tmp(:), eigen0(:)
 complex(dpc),allocatable :: eigen0_c(:,:) ! states,temp
 real(dp),allocatable :: eigen11(:),eigen12(:),eigtmp(:)
 real(dp),allocatable :: eigen13(:),occ(:),wtk(:)
 complex(dpc),allocatable :: pmat(:,:,:,:,:)
 character(len=fnlen) :: filnam,wfkfile,ddkfile_1,ddkfile_2,ddkfile_3,filnam_out, epfile
!  for the moment this is imposed by the format in linopt.f and nlinopt.f
 character(len=256) :: fn_radix,tmp_radix
 character(len=10) :: s1,s2,s3
 character(len=24) :: codename
 character(len=24) :: start_datetime
 character(len=500) :: msg
 type(hdr_type) :: hdr
 type(ebands_t) :: BSt, EPBSt
 integer :: finunt
 namelist /FILES/ ddkfile_1, ddkfile_2, ddkfile_3, wfkfile
 namelist /PARAMETERS/ broadening, domega, maxomega, scissor, tolerance, do_antiresonant, do_temperature, &
                       autoparal, max_ncpus
 namelist /COMPUTATIONS/ num_lin_comp, lin_comp, num_nonlin_comp, nonlin_comp, &
&        num_linel_comp, linel_comp, num_nonlin2_comp, nonlin2_comp
 namelist /TEMPERATURE/ epfile
 !character(len=fnlen) :: test
 integer :: iomode
 integer :: comm,nproc,my_rank
 type(wfk_t) :: wfk0,wfk1,wfk2,wfk3

 !_EP_NC reading !
 type(eprenorms_t) :: Epren
 character(len=fnlen) :: ep_nc_fname
 integer :: ep_ntemp
 logical :: do_ep_renorm
 integer :: nband_tmp
 integer :: itemp
 character(len=10) :: stemp
 type(crystal_t) :: Cryst
 logical :: remove_inv

! *********************************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 call xmpi_init()

 comm = xmpi_world
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

!Initialize memory profiling if it is activated
!if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 call timein(tcpui,twall)

 call timein(tcpui,twalli)
 start_datetime = asctime()

 if (my_rank == master) then
   codename='OPTIC '//repeat(' ',18)
   call herald(codename,abinit_version,std_out)

   !Read data file name
   write(std_out,'(a)')' Please, give the name of the data file ...'
   read(5, '(a)')filnam
   write(std_out,'(a)')' The name of the data file is :',filnam
   write(std_out,'(a)')' Please, give the name of the output file ...'
   read(5, '(a)')filnam_out
   write(std_out,'(a)')' The name of the output file is :',filnam_out
   write(std_out,'(a)')' Please, give the root name for the (non)linear optical data output file ...'
   read(5, '(a)')fn_radix
   write(std_out,'(a)')' The root name of the output files is :',trim(fn_radix)

   !Read data file
   if (open_file(filnam,msg,newunit=finunt,form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if

   write(msg,'(3a)') "From version 7.11.4, optic uses namelists as input.",ch10,&
&   "See e.g. ~/tests/tutorespfn/Input/toptic_2.in"
   MSG_COMMENT(msg)

   ! Setup some default values:
   broadening = 1e-3_dp ! Ha
   domega = 1e-3_dp ! Ha
   maxomega = 1.0_dp ! Ha
   scissor = 0.0_dp ! no scissor by default
   tolerance = 1e-3_dp ! Ha
   do_antiresonant = .TRUE. ! do use antiresonant approximation (only resonant transitions in the calculation)
   do_temperature = .FALSE.

   ! Read input file
   read(finunt,nml=FILES)
   read(finunt,nml=PARAMETERS)
   read(finunt,nml=COMPUTATIONS)

   if(do_temperature) then
     read(finunt, nml=TEMPERATURE)
   end if

   close(finunt)

   ! Validate input
   if (num_nonlin_comp > 0) then
     if (all(nonlin_comp(1:num_nonlin_comp) == 0)) then
       MSG_ERROR("nonlin_comp must be specified when num_nonlin_comp > 0")
     end if
   end if

   if (num_linel_comp > 0) then
     if (all(linel_comp(1:num_linel_comp) == 0)) then
       MSG_ERROR("linel_comp must be specified when num_linel_comp > 0")
     end if
   end if

   if (num_nonlin2_comp > 0) then
     if (all(nonlin2_comp(1:num_nonlin2_comp) == 0)) then
       MSG_ERROR("nonlin2_comp must be specified when num_nonlin2_comp > 0")
     end if
   end if

   ! Open the Wavefunction files
   ! Note: Cannot use MPI-IO here because of prtwf=3.
   ! If prtwf==3, the DDK file does not contain the wavefunctions but 
   ! this info is not reported in the header and the offsets in wfk_compute_offsets 
   ! are always computed assuming the presence of the cg

   ! TODO: one should perform basic consistency tests for the GS WFK and the DDK files, e.g.
   ! k-points and their order, spins, number of bands could differ in the four files.
   ! Note indeed that we are assuming the same numer of bands in all the files.
   call nctk_fort_or_ncfile(wfkfile, iomode, msg)
   if (len_trim(msg) /= 0) MSG_ERROR(msg)
   if (iomode == IO_MODE_MPI) iomode = IO_MODE_FORTRAN
   call wfk_open_read(wfk0,wfkfile,formeig0,iomode,10,xmpi_comm_self)

   call nctk_fort_or_ncfile(ddkfile_1, iomode, msg)
   if (len_trim(msg) /= 0) MSG_ERROR(msg)
   if (iomode == IO_MODE_MPI) iomode = IO_MODE_FORTRAN
   call wfk_open_read(wfk1,ddkfile_1,formeig1,iomode,11,xmpi_comm_self)

   call nctk_fort_or_ncfile(ddkfile_2, iomode, msg)
   if (len_trim(msg) /= 0) MSG_ERROR(msg)
   if (iomode == IO_MODE_MPI) iomode = IO_MODE_FORTRAN
   call wfk_open_read(wfk2,ddkfile_2,formeig1,iomode,12,xmpi_comm_self)

   call nctk_fort_or_ncfile(ddkfile_3, iomode, msg)
   if (len_trim(msg) /= 0) MSG_ERROR(msg)
   if (iomode == IO_MODE_MPI) iomode = IO_MODE_FORTRAN
   call wfk_open_read(wfk3,ddkfile_3,formeig1,iomode,13,xmpi_comm_self)

   if (wfk_compare(wfk1, wfk2) /= 0) then
     MSG_ERROR("ddfile_1 and ddkfile_2 are not consistent. see above messages")
   end if
   if (wfk_compare(wfk1, wfk3) /= 0) then
     MSG_ERROR("ddfile_1 and ddkfile_3 are not consistent. see above messages")
   end if

   !  Read the header from the gs file
   call hdr_copy(wfk0%hdr, hdr)

   ep_nc_fname = 'test_EP.nc'
   if(do_temperature) then
     ep_nc_fname = epfile
   end if
   do_ep_renorm = file_exists(ep_nc_fname)
   if(do_ep_renorm) then
     call eprenorms_from_epnc(Epren,ep_nc_fname)
     ep_ntemp = Epren%ntemp
   else if(do_temperature) then
     MSG_ERROR("You have asked for temperature but the epfile is not present !")
   else
     ep_ntemp = 1
   end if


   ! autoparal section
   if (autoparal /= 0 .and. max_ncpus /= 0) then
     write(std_out,'(a)')"--- !Autoparal"
     write(std_out,"(a)")'#Autoparal section for Optic runs.'

     write(std_out,"(a)")   "info:"
     write(std_out,"(a,i0)")"    autoparal: ",autoparal
     write(std_out,"(a,i0)")"    max_ncpus: ",max_ncpus
     write(std_out,"(a,i0)")"    nspinor: ",hdr%nspinor
     write(std_out,"(a,i0)")"    nsppol: ",hdr%nsppol
     write(std_out,"(a,i0)")"    nkpt: ",hdr%nkpt
     write(std_out,"(a,i0)")"    mband: ",maxval(hdr%nband)

     work_size = hdr%nkpt !* hdr%nsppol

     ! List of configurations.
     write(std_out,"(a)")"configurations:"
     do ii=1,max_ncpus
       if (ii > work_size) cycle
       nks_per_proc = work_size / ii
       nks_per_proc = nks_per_proc + MOD(work_size, ii)
       eff = (one * work_size) / (ii * nks_per_proc)
       write(std_out,"(a,i0)")"    - tot_ncpus: ",ii !* omp_ncpus
       write(std_out,"(a,i0)")"      mpi_ncpus: ",ii
       write(std_out,"(a,i0)")"      omp_ncpus: ",1
       write(std_out,"(a,f12.9)")"      efficiency: ",eff
       !write(,"(a,f12.2)")"      mem_per_cpu: ",mempercpu_mb 
     end do

     write(std_out,'(a)')"..."
     MSG_ERROR_NODUMP("aborting now")
     !call xmpi_abort(msg="Aborting now")
   end if 

 end if

 call hdr_bcast(hdr,master,my_rank,comm)
 !TODO put parameters in datastructure
 call xmpi_bcast(broadening,master,comm,ierr)
 call xmpi_bcast(domega,master,comm,ierr)
 call xmpi_bcast(maxomega,master,comm,ierr)
 call xmpi_bcast(scissor,master,comm,ierr)
 call xmpi_bcast(tolerance,master,comm,ierr)
 call xmpi_bcast(num_lin_comp,master,comm,ierr)
 call xmpi_bcast(lin_comp,master,comm,ierr)
 call xmpi_bcast(num_nonlin_comp,master,comm,ierr)
 call xmpi_bcast(nonlin_comp,master,comm,ierr)
 call xmpi_bcast(num_linel_comp,master,comm,ierr)
 call xmpi_bcast(linel_comp,master,comm,ierr)
 call xmpi_bcast(num_nonlin2_comp,master,comm,ierr)
 call xmpi_bcast(nonlin2_comp,master,comm,ierr)
 call xmpi_bcast(do_antiresonant,master,comm,ierr)

 call xmpi_bcast(do_ep_renorm,master,comm,ierr)
 call xmpi_bcast(ep_ntemp,master,comm,ierr)

 if(do_ep_renorm) then
   call eprenorms_bcast(Epren, master, comm)
 end if

!Extract info from the header
 headform=hdr%headform
 bantot=hdr%bantot
 ecut=hdr%ecut_eff
 natom=hdr%natom
 nkpt=hdr%nkpt
 nspinor=hdr%nspinor
 nsppol=hdr%nsppol
 ntypat=hdr%ntypat
 occopt=hdr%occopt
 rprimd(:,:)=hdr%rprimd(:,:)
 ABI_ALLOCATE(nband,(nkpt*nsppol))
 ABI_ALLOCATE(occ,(bantot))
 !fermie=hdr%fermie
 ! YG Fermi energy contained in the header of a NSCF computation is always 0 !!
 occ(1:bantot)=hdr%occ(1:bantot)
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)

 nsym=hdr%nsym
 ABI_ALLOCATE(symrel,(3,3,nsym))
 ABI_ALLOCATE(symrec,(3,3,nsym))
 symrel(:,:,:) = hdr%symrel(:,:,:)
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
 end do

 ABI_ALLOCATE(kpt,(3,nkpt))
 kpt(:,:) = hdr%kptns(:,:)

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(nband(:))
 do ii=1,nkpt
   if (nband(ii) /= mband) then
     MSG_ERROR("nband must be constant across kpts")
   end if
 end do

 ! Initializes crystal
 remove_inv = .false.
 if(hdr%nspden == 4) remove_inv = .true.
 !space_group = 0
 call crystal_init(hdr%amu, Cryst, 0, hdr%natom, hdr%npsp, hdr%ntypat, &
&  hdr%nsym, rprimd, hdr%typat, hdr%xred, hdr%zionpsp, hdr%znuclpsp, 1, &
&  (hdr%nspden==2 .and. hdr%nsppol==1),remove_inv, hdr%title,&
&  symrel, hdr%tnons, hdr%symafm)

 if(my_rank == master) then
   write(std_out,*)
   write(std_out,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1:3,1)
   write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,2)
   write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,3)
   write(std_out,'(a,i8)')       ' natom             =',natom
   write(std_out,'(a,2i8)')      ' nkpt,mband        =',nkpt,mband
   write(std_out,'(a, f10.5,a)' ) ' ecut              =',ecut,' Ha'
 end if

!Read the eigenvalues of ground-state and ddk files
 ABI_ALLOCATE(eigen0,(mband*nkpt*nsppol))
 ABI_ALLOCATE(eigen0_c,(mband*nkpt*nsppol,ep_ntemp))
 ABI_ALLOCATE(eigen11,(2*mband*mband*nkpt*nsppol))
 ABI_ALLOCATE(eigen12,(2*mband*mband*nkpt*nsppol))
 ABI_ALLOCATE(eigen13,(2*mband*mband*nkpt*nsppol))

 if(my_rank == master) then

   ABI_ALLOCATE(eigtmp,(2*mband*mband))
   ABI_ALLOCATE(eig0tmp,(mband))
   bdtot0_index=0 ; bdtot_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       nband1=nband(ikpt+(isppol-1)*nkpt)
       eigtmp = zero
       eig0tmp = zero
       
       call wfk_read_eigk(wfk0,ikpt,isppol,xmpio_single,eig0tmp)
       eigen0(1+bdtot0_index:nband1+bdtot0_index)=eig0tmp(1:nband1)
       call wfk_read_eigk(wfk1,ikpt,isppol,xmpio_single,eigtmp)
       eigen11(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)
       call wfk_read_eigk(wfk2,ikpt,isppol,xmpio_single,eigtmp)
       eigen12(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)
       call wfk_read_eigk(wfk3,ikpt,isppol,xmpio_single,eigtmp)
       eigen13(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)

       !ABI_CHECK(wfk1%nband(ikpt,isppol) == nband1, "ddk1 nband1")
       !ABI_CHECK(wfk2%nband(ikpt,isppol) == nband1, "ddk2 nband1")
       !ABI_CHECK(wfk3%nband(ikpt,isppol) == nband1, "ddk3 nband1")

       bdtot0_index=bdtot0_index+nband1
       bdtot_index=bdtot_index+2*nband1**2

       !write(std_out,*)"eig0tmp:"
       !write(std_out,*)eig0tmp
       !write(std_out,*)"eigtmp",ikpt,isppol
       !do ii=1,2*nband1**2
       !  write(std_out,*)ii,eigtmp(ii)
       !end do

     end do
   end do

   call wfk_close(wfk0)
   call wfk_close(wfk1)
   call wfk_close(wfk2)
   call wfk_close(wfk3)
   
   ABI_DEALLOCATE(eigtmp)
   ABI_DEALLOCATE(eig0tmp)

 end if

 call xmpi_bcast(eigen0,master,comm,ierr)
 call xmpi_bcast(eigen11,master,comm,ierr)
 call xmpi_bcast(eigen12,master,comm,ierr)
 call xmpi_bcast(eigen13,master,comm,ierr)


!---------------------------------------------------------------------------------
!gmet inversion
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 call matr3inv(gmet,gmet_inv)

!---------------------------------------------------------------------------------
!derivative of occupation wrt the energy.
 ABI_ALLOCATE(wtk,(nkpt))
 wtk = hdr%wtk

 ABI_ALLOCATE(doccde,(mband*nkpt*nsppol))

 !Recompute fermie from header
 !WARNING no garantie that it works for other materials than insulators
 nelect = hdr%nelect
 tphysel = zero
 ABI_ALLOCATE(istwfk,(nkpt))
 ABI_ALLOCATE(npwarr,(nkpt))
 istwfk = hdr%istwfk
 npwarr = hdr%npwarr

 call ebands_init(bantot, BSt, nelect, doccde, eigen0, istwfk, kpt, &
& nband, nkpt, npwarr, nsppol, nspinor, tphysel, broadening, occopt, occ, wtk, &
& hdr%charge, hdr%kptopt, hdr%kptrlatt_orig, hdr%nshiftk_orig, hdr%shiftk_orig, &
& hdr%kptrlatt, hdr%nshiftk, hdr%shiftk)

 !YG : should we use broadening for ebands_init
 call ebands_update_occ(BSt, -99.99d0)
 fermie = BSt%fermie
 ABI_DEALLOCATE(istwfk)
 ABI_DEALLOCATE(npwarr)

!---------------------------------------------------------------------------------
!size of the frequency range
 nomega=int((maxomega+domega*0.001_dp)/domega)
 maxomega = dble(nomega)*domega
 ABI_ALLOCATE(cond_nd,(nomega))
 ABI_ALLOCATE(cond_kg,(nomega))

 if(my_rank == master) then
   write(std_out,'(a,f10.5,a,f10.5,a)' )' fermie            =',fermie,' Ha',fermie*Ha_eV,' eV'
   write(std_out,'(a,f10.5,a)')' Scissor shift     =', scissor, ' Ha'
   write(std_out,'(a,f10.5,a)')' Tolerance on closeness to singularities     =', tolerance, ' Ha'
   write(std_out,'(a,f10.5,a)')' Smearing factor      =', broadening, ' Ha'
   if(do_antiresonant) then
     write(std_out,'(a)') ' Will use the antiresonant approximation '
   else
     write(std_out,'(a)') ' Will not use the antiresonant approximation (only available for nonlin2 and linel components!) '
   end if
   write(std_out,'(a)') ' linear coeffs to be calculated : '
   write(std_out,'(9i3)') lin_comp(1:num_lin_comp)
   write(std_out,'(a)') ' non-linear coeffs to be calculated : '
   write(std_out,'(27i4)') nonlin_comp(1:num_nonlin_comp)
   write(std_out,'(a)') ' electronic part of electro-optic coeffs to be calculated :'
   write(std_out,'(27i4)') linel_comp(1:num_linel_comp)
   write(std_out,'(a)') ' non-linear coeffs (V2) to be calculated :'
   write(std_out,'(27i4)') nonlin2_comp(1:num_nonlin2_comp)
 end if
 
 ABI_ALLOCATE(symcart,(3,3,nsym))
 !YG : we need to transpose gprimd since matrinv give the transpose of the inverse !
 gprimd_trans = transpose(gprimd)
 call sym2cart(gprimd_trans,nsym,rprimd,symrel,symcart)

 ABI_ALLOCATE(pmat,(mband,mband,nkpt,3,nsppol))
 call wrtout(std_out," optic : Call pmat2cart","COLL")

 call pmat2cart(eigen11,eigen12,eigen13,mband,nkpt,nsppol,pmat,rprimd)

 call pmat_renorm(fermie, eigen0, mband, nkpt, nsppol, pmat, scissor)

 !TODO implement renormalization in ebands
 if(my_rank == master) then

   if(do_ep_renorm) then
     bdtot0_index=0 ; bdtot_index=0
     do isppol=1,nsppol
       do ikpt=1,nkpt
         call flush(std_out)
         nband1 = nband(ikpt+(isppol-1)*nkpt)
         nband_tmp=MIN(nband1,Epren%mband)

         !FIXME change check
         if (ANY(ABS(eigen0(1+bdtot0_index:MIN(10,nband_tmp)+bdtot0_index) - &
&             Epren%eigens(1:MIN(10,nband_tmp),ikpt,isppol)) > tol3)) then
           MSG_ERROR("Error in eigen !")
         end if
         if (ANY(ABS(occ(1+bdtot0_index:MIN(10,nband_tmp)+bdtot0_index) - &
&             Epren%occs(1:MIN(10,nband_tmp),ikpt,isppol)) > tol3)) then
           MSG_ERROR("Error in occ!")
         end if

         ! Upgrade energies
         do itemp = 1,Epren%ntemp
           eigen0_c(1+bdtot0_index:nband_tmp+bdtot0_index,itemp) = eigen0(1+bdtot0_index:nband_tmp+bdtot0_index) + Epren%renorms(1,1:nband_tmp,ikpt,isppol,itemp) + (0.0_dp,1.0_dp)*Epren%lifetimes(1,1:nband_tmp,ikpt,isppol,itemp)
         end do

         bdtot0_index=bdtot0_index+nband1
         bdtot_index=bdtot_index+2*nband1**2
       end do
     end do
   else
     ! Only 1 temperature (no temperature !)!
     eigen0_c(:,1) = eigen0(:)

   end if


 end if
 call xmpi_bcast(eigen0_c,master,comm,ierr)

!IN CALLED ROUTINE
!call linopt(nspin,,nkpt,wkpt,nsymcrys,symcrys,nstval,occv,evalv,efermi,pmat,v1,v2,nmesh,de,scissor,brod)
!
!v1,v2=desired component of the dielectric function(integer) 1=x,2=y,3=z
!nmesh=desired number of energy mesh points(integer)
!de=desired step in energy(real); nmesh*de=maximum energy
!scissor=scissors shift in Ha(real)
!brod=broadening in Ha(real)
!
 call wrtout(std_out," optic : Call linopt","COLL")

 do itemp=1,ep_ntemp
   call ebands_copy(BSt, EPBst)
   if(do_ep_renorm) then
     call renorm_bst(Epren, EPBst, Cryst, itemp, do_lifetime=.True.,do_check=.True.)
   end if
   do ii=1,num_lin_comp
     lin1 = int(lin_comp(ii)/10.0_dp)
     lin2 = mod(lin_comp(ii),10)
     write(msg,*) ' linopt ', lin1,lin2
     call wrtout(std_out,msg,"COLL")
     call int2char4(lin1,s1)
     call int2char4(lin2,s2)
     call int2char4(itemp,stemp)
     ABI_CHECK((s1(1:1)/='#'),'Bug: string length too short!')
     ABI_CHECK((s2(1:1)/='#'),'Bug: string length too short!')
     ABI_CHECK((stemp(1:1)/='#'),'Bug: string length too short!')
     if(do_ep_renorm) then
       tmp_radix = trim(fn_radix)//"_"//trim(s1)//"_"//trim(s2)//"_T"//trim(stemp)
     else
       tmp_radix = trim(fn_radix)//"_"//trim(s1)//"_"//trim(s2)
     end if
     call linopt(nsppol,ucvol,nkpt,wtk,nsym,symcart,mband,BSt,EPBst,fermie,pmat, &
     lin1,lin2,nomega,domega,scissor,broadening,tmp_radix,comm)
   end do
   call ebands_free(EPBst)
 end do
 
 if(do_ep_renorm) then
   call eprenorms_free(Epren)
 end if
 
 call wrtout(std_out," optic : Call nlinopt","COLL")
 
 do ii=1,num_nonlin_comp
   nlin1 = int( nonlin_comp(ii)/100.0_dp)
   nlin2 = int((nonlin_comp(ii)-nlin1*100.0_dp)/10.0_dp)
   nlin3 = mod( nonlin_comp(ii),10)
   write(msg,*) ' nlinopt ', nlin1,nlin2,nlin3
   call wrtout(std_out,msg,"COLL")
   call int2char4(nlin1,s1)
   call int2char4(nlin2,s2)
   call int2char4(nlin3,s3)
   ABI_CHECK((s1(1:1)/='#'),'Bug: string length too short!')
   ABI_CHECK((s2(1:1)/='#'),'Bug: string length too short!')
   ABI_CHECK((s3(1:1)/='#'),'Bug: string length too short!')
   tmp_radix = trim(fn_radix)//"_"//trim(s1)//"_"//trim(s2)//"_"//trim(s3)
   call nlinopt(nsppol,ucvol,nkpt,wtk,nsym,symcart,mband,eigen0,fermie,pmat,&
&   nlin1,nlin2,nlin3,nomega,domega,scissor,broadening,tolerance,tmp_radix,comm)
 end do

 call wrtout(std_out," optic : Call linelop","COLL")

 do ii=1,num_linel_comp
   linel1 = int( linel_comp(ii)/100.0_dp)
   linel2 = int((linel_comp(ii)-linel1*100.0_dp)/10.0_dp)
   linel3 = mod( linel_comp(ii),10)
   write(msg,*) ' linelop ',linel1,linel2,linel3
   call wrtout(std_out,msg,"COLL")
   call int2char4(linel1,s1)
   call int2char4(linel2,s2)
   call int2char4(linel3,s3)
   tmp_radix = trim(fn_radix)//"_"//trim(s1)//"_"//trim(s2)//"_"//trim(s3)
   call linelop(nsppol,ucvol,nkpt,wtk,nsym,symcart,mband,eigen0,occ,fermie,pmat,&
&   linel1,linel2,linel3,nomega,domega,scissor,broadening,tolerance,tmp_radix,do_antiresonant,comm)
 end do

 call wrtout(std_out," optic : Call nonlinopt","COLL")
 
 do ii=1,num_nonlin2_comp
   nonlin1 = int( nonlin2_comp(ii)/100.0_dp)
   nonlin2 = int((nonlin2_comp(ii)-nonlin1*100.0_dp)/10.0_dp)
   nonlin3 = mod( nonlin2_comp(ii),10)
   write(msg,*) ' nonlinopt ',nonlin1,nonlin2,nonlin3
   call wrtout(std_out,msg,"COLL")
   call int2char4(nonlin1,s1)
   call int2char4(nonlin2,s2)
   call int2char4(nonlin3,s3)
   tmp_radix = trim(fn_radix)//"_"//trim(s1)//"_"//trim(s2)//"_"//trim(s3)
   call nonlinopt(nsppol,ucvol,nkpt,wtk,nsym,symcart,mband,eigen0,occ,fermie,pmat,&
&   nonlin1,nonlin2,nonlin3,nomega,domega,scissor,broadening,tolerance,tmp_radix,do_antiresonant,comm)
 end do

 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(eigen11)
 ABI_DEALLOCATE(eigen12)
 ABI_DEALLOCATE(eigen13)
 ABI_DEALLOCATE(eigen0)
 ABI_DEALLOCATE(eigen0_c)
 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(wtk)
 ABI_DEALLOCATE(cond_nd)
 ABI_DEALLOCATE(cond_kg)

 ABI_DEALLOCATE(kpt)
 ABI_DEALLOCATE(symrel)
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(symcart)
 ABI_DEALLOCATE(pmat)

 call ebands_free(BSt)
 call hdr_free(hdr)
 call crystal_free(Cryst)

 call timein(tcpu,twall)

 tsec(1)=tcpu-tcpui
 tsec(2)=twall-twalli

 if(my_rank == master) then
   write(std_out,'(a,80a,a,a,a)' )ch10,('=',ii=1,80),ch10,ch10,' Calculation completed.'

   write(std_out, '(a,a,a,f13.1,a,f13.1)' ) &
&   '-',ch10,'- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
 end if

 call xmpi_sum(tsec,comm,ierr)

 if(my_rank == master) then
   ! Write YAML document with the final summary.
   ! we use this doc to test whether the calculation is completed.
   write(std_out,"(a)")"--- !FinalSummary"
   write(std_out,"(a)")"program: optic"
   write(std_out,"(2a)")"version: ",trim(abinit_version)
   write(std_out,"(2a)")"start_datetime: ",start_datetime
   write(std_out,"(2a)")"end_datetime: ",asctime()
   write(std_out,"(a,f13.1)")"overall_cpu_time: ",tsec(1)
   write(std_out,"(a,f13.1)")"overall_wall_time: ",tsec(2)
   write(std_out,"(a,i0)")"mpi_procs: ",xmpi_comm_size(xmpi_world)
   write(std_out,"(a,i0)")"omp_threads: ",xomp_get_num_threads(open_parallel=.True.)
   !write(std_out,"(a,i0)")"num_warnings: ",nwarning
   !write(std_out,"(a,i0)")"num_comments: ",ncomment
   write(std_out,"(a)")"..."
   call flush_unit(std_out)
 end if

!Write information on file about the memory before ending mpi module, if memory profiling is enabled
 call abinit_doctor(filnam)

 call xmpi_end()

 end program optic
!!***
