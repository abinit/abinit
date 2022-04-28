!!****m* ABINIT/m_conducti
!! NAME
!!  m_conducti
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula.
!!
!! COPYRIGHT
!!  Copyright (C) 2002-2022 ABINIT group (VRecoules, PGhosh, SMazevet, SM, SVinko, NBrouwer)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_conducti

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_wfk
 use m_hdr
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_abitypes,  only : MPI_type
 use m_io_tools,     only : open_file, close_unit, get_unit
 use m_fstrings,     only : sjoin
 use m_symtk,        only : matr3inv
 use m_hide_lapack,  only : jacobi
 use m_occ,          only : getnel
 use m_geometry,     only : metric
 use m_splines,      only : intrpl,splint,spline
 use m_mpinfo,       only : distrb2,init_mpi_enreg,destroy_mpi_enreg,proc_distrb_cycle,initmpi_band

 implicit none

 private
!!***

 public :: conducti_paw
 public :: conducti_paw_core
 public :: conducti_nc

!I/O parameters
!Set to true to use netcdf-MPIIO when available
 logical,parameter :: use_netcdf_mpiio=.true.

!!***

contains
!!***

!!****f* m_conducti/conducti_paw
!! NAME
!! conducti_paw
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula for PAW formalism
!!
!! INPUTS
!!  filnam=generic name for input data
!!  filnam_out=generic name for output data
!!  [varocc]=if true, read arbitrary occupations from a file
!!
!! OUTPUT
!!  Only printing
!!
!! NOTES
!!  bantot
!!  doccde(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy.
!!  dom=frequency range
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree).
!!  eigen11(2,nkpt,mband,mband,nsppol)=first-order eigenvalues (hartree)
!!  in direction x
!!  eigen12(2,nkpt,mband,mband,nsppol)=first-order eigenvalues (hartree)
!!  in direction y
!!  eigen13(2,nkpt,mband,mband,nsppol)=first-order eigenvalues (hartree)
!!  in direction z
!!  ecut=kinetic energy planewave cutoff (hartree).
!!  fermie= fermi energy (Hartree)
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{2}$).
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1).
!!  kin11= Onsager kinetic coeficient=optical conductivity
!!  kin12= Onsager kinetic coeficient
!!  kin21= Onsager kinetic coeficient
!!  kin22= Onsager kinetic coeficient
!!  Kth=thermal conductivity
!!  mom=number of frequency for conductivity computation
!!  mband=maximum number of bands.
!!  natom = number of atoms in the unit cell.
!!  nband(nkpt*nsppol)=number of bands at each RF k point for each spin.
!!  nkpt=number of k points in the IBZ for this perturbation
!!  ngfft(3)=integer fft box dimensions.
!!  nspinor=number of spinorial components of the wavefunctions.
!!  nsppol=1 for unpolarized, 2 for spin-polarized.
!!  ntypat = number of atom types.
!!  occ(mband*nkpt*nsppol)=occupation number for each band and k.
!!  occopt==option for occupancies
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).sigx(mom,nphicor))
!!  rprimd(3,3)=real space primitive translations.
!!  of primitive translations.
!!  Sth=thermopower
!!  tsmear=smearing width (or temperature) in Hartree
!!  ucvol=unit cell volume in ($\textrm{bohr}^{3}$).
!!  wind=frequency windows for computations of sigma
!!  wtk(nkpt)=weight assigned to each k point.
!!  znucl(natom)=atomic number of atoms
!!  np_sum=noziere-pines sumrule
!!
!! PARENTS
!!      conducti
!!
!! CHILDREN
!!      close_unit,getnel,hdr_fort_read,hdr_ncread,hdr%bcast,hdr%free
!!      metric,msig,nctk_fort_or_ncfile,nctk_open_read,nctk_set_collective
!!      nf90_close,nf90_get_var,open_file,spline,splint,xmpi_bcast,xmpi_sum
!!
!! SOURCE

 subroutine conducti_paw(filnam,filnam_out,varocc)

!Arguments -----------------------------------
!scalars
 character(len=fnlen) :: filnam,filnam_out
 integer, optional :: varocc

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: bsize,bd_stride,dimid,iomode,bantot,bdtot_index,ncid,varid,nb_per_proc,etiq
 integer :: comm,fform1,headform,iband,ijband,ierr,ikpt,master_band
 integer :: iom,isppol,jband,l1,l2,mband,me,mpierr,mom
 integer :: natom,nband_k,nkpt,nproc,nspinor,nsppol,ntypat
 integer :: occopt,iunt,opt_unt,occ_unt,iocc,my_iband
 integer :: lij_unt,sig_unt,kth_unt,ocond_unt,occunit,occnpt
 logical :: nc_unlimited,mykpt,myband,iomode_estf_mpiio,read_half_dipoles
 real(dp) :: dirac,del,deltae,deltae_min,deltae_min_tmp,dhdk2_g,diff_eig,diff_occ
 real(dp) :: dosdeltae,ecut,entropy,fact_diag,fermie,fermih,kin_fact,maxocc
 real(dp) :: np_sum,np_sum_k1,np_sum_k2,omin,omax,dom,oml,sig,socc,socc_k
 real(dp) :: Tatm,tphysel,tsmear,ucvol,eig_in_max,eig_in_min
 character(len=fnlen) :: filnam1,filnam_gen,occfile
 character(len=500) :: msg
 type(hdr_type) :: hdr
 type(MPI_type) :: mpi_enreg
!arrays
 integer :: nc_count_5(5),nc_count_6(6),nc_start_5(5),nc_start_6(6),nc_stride_5(5),nc_stride_6(6)
 integer,allocatable :: nband(:)
 real(dp) :: dhdk2_r(3,3),gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: cond_nd(:,:,:),cond_nd_k(:,:,:)
 real(dp),allocatable :: doccde(:),doccde_k(:),eig0_k(:),eigen0(:),eig0nc(:,:,:)
 real(dp),allocatable :: occ(:),occ_k(:),wtk(:),oml1(:),occ_in(:),eig_in(:),occ_tmp(:),ypp(:)
 real(dp),allocatable :: kin11(:,:),kin12(:),kin21(:),kin22(:)
 real(dp),allocatable :: kin11_k(:),kin12_k(:),kin21_k(:),kin22_k(:),Kth(:),Stp(:)
 real(dp),allocatable :: psinablapsi(:,:,:),sig_abs(:)
 
! *********************************************************************************
 
! ---------------------------------------------------------------------------------
! Read input data

!Global MPI communicator
 comm = xmpi_world
 nproc = xmpi_comm_size(comm)
 me = xmpi_comm_rank(comm)

!Read input parameters file
 if (me==master) then
   if (open_file(filnam,msg,newunit=iunt,form='formatted',action="read",status="old")/=0) then
     ABI_ERROR(msg)
   end if
   rewind(iunt)
   read(iunt,*)
   read(iunt,'(a)') filnam_gen ! Generic name for the files
   filnam1=trim(filnam_gen)//'_OPT'
!  Read frequency range
   read(iunt,*) dom,omin,omax,mom
!  In case of varocc read filename of occupation datafile
   if (present(varocc)) then
     if (me==master) then
       write(std_out,*) 'Warning, this undocumented feature is highly experimental'
       write(std_out,*) 'and of limited physical validity, proceed with extreme caution!!!'
     end if
     read(iunt,*) occfile !filename of the non-eq distribution function
   end if
   close(iunt)
 end if

!Send data to all procs
 call xmpi_bcast(dom,master,comm,mpierr)
 call xmpi_bcast(omin,master,comm,mpierr)
 call xmpi_bcast(omax,master,comm,mpierr)
 call xmpi_bcast(mom,master,comm,mpierr)

! ---------------------------------------------------------------------------------
! Read OPT file

!Check for FORTRAN/.nc OPT file and set iomode to IO_MODE_FORTRAN_MASTER/IO_MODE_ETSF
 if (me==master) then
   call nctk_fort_or_ncfile(filnam1,iomode,msg)
   if (iomode/=IO_MODE_ETSF) iomode=IO_MODE_FORTRAN_MASTER
 end if
 call xmpi_bcast(filnam1,master,comm,mpierr)
 call xmpi_bcast(iomode,master,comm,mpierr)

 !Open OPT file and read HEADER
 if (me==master) then
   if (iomode==IO_MODE_ETSF) then
     NCF_CHECK(nctk_open_read(ncid,filnam1,xmpi_comm_self))
     call hdr_ncread(hdr,ncid,fform1)
   else
     if (open_file(filnam1,msg,newunit=opt_unt,form="unformatted",status="old")/=0) then
       ABI_ERROR(msg)
     end if
     call hdr_fort_read(hdr,opt_unt,fform1,rewind=.true.)
   end if 
   ABI_CHECK(fform1/=0,sjoin("Error while reading ",filnam1))
   ABI_CHECK(fform1==610.or.fform1==620,"Conducti requires an OPT file with fform=610 or 620!")
 end if
 call hdr%bcast(master,me,comm)
 call xmpi_bcast(fform1,master,comm,mpierr)

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
 fermie=hdr%fermie
 tsmear=hdr%tsmear
 ABI_MALLOC(nband,(nkpt*nsppol))
 ABI_MALLOC(occ,(bantot))
 ABI_MALLOC(wtk,(nkpt))
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)
 occ(1:bantot)=hdr%occ(1:bantot)
 wtk(1:nkpt)=hdr%wtk(1:nkpt)
 mband=maxval(nband(:))  ! Get mband, as the maximum value of nband(nkpt)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol) ! Get metrics of simulation cell

!Read eigenvalues
 ABI_MALLOC(eigen0,(mband*nkpt*nsppol))
 if (me==master) then
   if (iomode==IO_MODE_ETSF) then
     varid=nctk_idname(ncid,"eigenvalues")
     ABI_MALLOC(eig0nc,(mband,nkpt,nsppol))
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_get_var(ncid,varid,eig0nc))
#endif
     eigen0 = reshape(eig0nc,[mband*nkpt*nsppol])
     ABI_FREE(eig0nc)
     !Close file here because the rest will possibly be read with collective I/O
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_close(ncid))
#endif
   else
     read(opt_unt)(eigen0(iband),iband=1,mband*nkpt*nsppol)
   end if
 end if
 call xmpi_bcast(eigen0,master,comm,mpierr)

!In case of varocc, use arbitrary occupations
!Read occfile and overwrite occ by interpolating the data in the input file
 if (present(varocc)) then
   if (me==master) then
     ABI_MALLOC(occ_tmp,(bantot))
     if (open_file(occfile,msg,newunit=occ_unt,form='formatted')/=0) then
       ABI_ERROR(msg)
     end if
!    Get units used in occfile (1=Ha, 2= eV)
     read(occ_unt,*) occunit
     read(occ_unt,*) occnpt
     ABI_MALLOC(occ_in,(occnpt))
     ABI_MALLOC(eig_in,(occnpt))
     do iocc=1,occnpt
       read(occ_unt,*) eig_in(iocc),occ_in(iocc)
       if (occunit==2) then !Convert units from eV to Hartree
         eig_in(iocc)=eig_in(iocc)/Ha_eV
         occ_in(iocc)=occ_in(iocc)
       end if
     end do
     close(occ_unt)
!    Interpolation
     ABI_MALLOC(ypp,(occnpt))
     call spline(eig_in,occ_in,occnpt,zero,zero,ypp)
!    Interpolate neccessary values
     eig_in_max=maxval(eig_in)
     eig_in_min=minval(eig_in)
     do iocc=1,bantot
!      Check for Extrapolation and set them to physically sound values
       if (eigen0(iocc)<eig_in_min) then
         occ_tmp(iocc)=one
       else if (eigen0(iocc)>eig_in_max) then
         occ_tmp(iocc)=zero
       else
         call splint(occnpt,eig_in,occ_in,ypp,1,eigen0(iocc),occ_tmp(iocc),ierr)
       end if
     end do
     occ=occ_tmp
!    Clean up
     ABI_FREE(ypp)
     ABI_FREE(occ_tmp)
     ABI_FREE(occ_in)
     ABI_FREE(eig_in)
   end if
   call xmpi_bcast(occ,master,comm,mpierr)
   occopt=2
 end if ! varocc?

!---------------------------------------------------------------------------------
! Prepare kpt/band parallelization

 call init_mpi_enreg(mpi_enreg)
 mpi_enreg%comm_kpt=comm
 mpi_enreg%me_kpt=me
 mpi_enreg%nproc_spkpt=nproc
 mpi_enreg%paralbd=1
 ABI_MALLOC(mpi_enreg%proc_distrb,(nkpt,mband,nsppol))
 ABI_MALLOC(mpi_enreg%my_kpttab,(nkpt))
 call distrb2(mband,nb_per_proc,nband,nkpt,nproc,nsppol,mpi_enreg)
 call initmpi_band(nkpt,mpi_enreg,nband,nkpt,nsppol)

!---------------------------------------------------------------------------------
!Print some data

 Tatm=tsmear*Ha_K
 if (me==master) then
   write(std_out,*)
   write(std_out,'(a)' )' Input data:'
   write(std_out,'(a,i8,3f10.5,a)')' npts,omin,omax,width      =',mom,omin,omax,dom,' Ha'
   write(std_out,*)
   write(std_out,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1:3,1)
   write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,2)
   write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,3)
   write(std_out,'(a,i8)')       ' natom             =',natom
   write(std_out,'(a,3i8)')      ' nkpt,mband,nsppol        =',nkpt,mband,nsppol
   write(std_out, '(a, f10.5,a)' ) ' ecut              =',ecut,' Ha'
   write(std_out,'(a,f10.5,a,f10.5,a)' )' fermie            =',fermie,' Ha',fermie*Ha_eV,' eV'
   write(std_out,'(a,f12.5,a,f12.5,a)') ' Temp              =',tsmear,' Ha ',Tatm,' Kelvin'
 end if

! ---------------------------------------------------------------------------------
! Compute derivative of occupations wrt the energy

 ABI_MALLOC(doccde,(mband*nkpt*nsppol))
 if (occopt==1) then
   if (me==master) then
     write(std_out,'(a,i4)')  ' occopt            =',occopt
   end if
   doccde=zero
 else
   tphysel=zero
   maxocc=two/(nsppol*nspinor)
 end if
 if (occopt==2) then
   tsmear=one
   entropy=one
   if (me==master) then
     write(std_out,'(a)') ' occopt=2 fixed occupation case'
   end if
 else
   call getnel(doccde,dosdeltae,eigen0,entropy,fermie,fermih,maxocc,mband,nband,&
&   socc,nkpt,nsppol,occ,occopt,1,tphysel,tsmear,12,wtk)
   entropy=tsmear*entropy
   if (me==master) then
     write(std_out, '(a,es22.12)') ' tsmear*entropy (Ha)           =',entropy
   end if
   entropy=entropy/socc
   if (me==master) then
     write(std_out, '(a,es22.12)') ' tsmear*entropy (Ha/e)         =',entropy
   end if
 endif

!---------------------------------------------------------------------------------
! Determine the frequency range and allocate frequency-dependent arrays

 del=(omax-omin)/(mom-1)
 ABI_MALLOC(oml1,(mom))
 do iom=1,mom
   oml1(iom)=omin+dble(iom-1)*del
 end do

 ABI_MALLOC(kin11,(mom,nsppol))
 ABI_MALLOC(kin12,(mom))
 ABI_MALLOC(kin21,(mom))
 ABI_MALLOC(kin22,(mom))
 ABI_MALLOC(cond_nd,(3,3,mom))
 ABI_MALLOC(sig_abs,(mom))
 ABI_MALLOC(Kth,(mom))
 ABI_MALLOC(Stp,(mom))
 kin11   = zero
 kin12   = zero
 kin21   = zero
 kin22   = zero
 cond_nd = zero
 sig_abs = zero
 Kth     = zero
 Stp     = zero

!---------------------------------------------------------------------------------
!Prepare valence-valence dipoles reading

 iomode_estf_mpiio=(iomode==IO_MODE_ETSF.and.nctk_has_mpiio.and.use_netcdf_mpiio)

 !In case of netCDF access to OPT file, prepare collective I/O
 if (iomode == IO_MODE_ETSF) then
   if (iomode_estf_mpiio) then
     NCF_CHECK(nctk_open_read(ncid,filnam1,comm))
     varid=nctk_idname(ncid,"dipole_valence_valence")
     if (nproc>1) then
       NCF_CHECK(nctk_set_collective(ncid,varid))
     end if
#ifdef HAVE_NETCDF
     nc_unlimited=(nf90_inq_dimid(ncid,"unlimited_bands",dimid)==NF90_NOERR)
     read_half_dipoles=(nf90_inq_dimid(ncid,"max_number_of_state_pairs",dimid)==NF90_NOERR)
#endif
   else
     if (me==master) then
       NCF_CHECK(nctk_open_read(ncid,filnam1,xmpi_comm_self))
       varid=nctk_idname(ncid,"dipole_valence_valence")
       !if (nctk_has_mpiio.and.(.not.use_netcdf_mpiio)) then
       !  NCF_CHECK(nctk_set_collective(ncid,varid))
       !end if
#ifdef HAVE_NETCDF
       nc_unlimited=(nf90_inq_dimid(ncid,"unlimited_bands",dimid)==NF90_NOERR)
       read_half_dipoles=(nf90_inq_dimid(ncid,"max_number_of_state_pairs",dimid)==NF90_NOERR)
#endif
     end if
     call xmpi_bcast(nc_unlimited,master,comm,ierr)
     call xmpi_bcast(read_half_dipoles,master,comm,ierr)
   end if
   if (nc_unlimited.and.read_half_dipoles) then
     msg="The OPT file has a wrong format!"
     ABI_BUG(msg)
   end if
 else
   read_half_dipoles=(fform1==620)
 end if

 if (iomode_estf_mpiio) then
   !If MPI-IO, store only ib elements for each jb
   ABI_MALLOC(psinablapsi,(2,3,mband))
 else
   !If not, store all pairs (or half)
   if (read_half_dipoles) then ! only ib>=jb
     ABI_MALLOC(psinablapsi,(2,3,(mband*(mband+1))/2))
   else
     ABI_MALLOC(psinablapsi,(2,3,mband*mband))
   end if
 end if

!---------------------------------------------------------------------------------
! Compute conductivity

 ABI_MALLOC(kin11_k,(mom))
 ABI_MALLOC(kin12_k,(mom))
 ABI_MALLOC(kin21_k,(mom))
 ABI_MALLOC(kin22_k,(mom))
 ABI_MALLOC(cond_nd_k,(3,3,mom))

 np_sum  = zero
 socc    = zero
 deltae  = zero
 deltae_min = 1.d99

!LOOP OVER SPINS/K
 bdtot_index = 0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     etiq=ikpt+(isppol-1)*nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     mykpt=.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me))
     master_band=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol))
     
!    In case of non MPI-IO, has to read all (n,m) dipoles for this k-point
!      Master node reads and send to relevant processor
     if (.not.iomode_estf_mpiio.and.me==master) then
       if (iomode==IO_MODE_ETSF) then
         psinablapsi=zero
#ifdef HAVE_NETCDF
         if (nc_unlimited) then
           nc_start_6=[1,1,1,ikpt,isppol,1] ; nc_count_6=[2,3,mband,1,1,mband] ; nc_stride_6=[1,1,1,1,1,1] 
           NCF_CHECK(nf90_get_var(ncid,varid,psinablapsi,start=nc_start_6,stride=nc_stride_6,count=nc_count_6))
         else if (.not.read_half_dipoles) then
           nc_start_6=[1,1,1,1,ikpt,isppol] ; nc_count_6=[2,3,mband,mband,1,1] ; nc_stride_6=[1,1,1,1,1,1] 
           NCF_CHECK(nf90_get_var(ncid,varid,psinablapsi,start=nc_start_6,stride=nc_stride_6,count=nc_count_6))
         else
           nc_start_5=[1,1,1,ikpt,isppol] ; nc_count_5=[2,3,(mband*(mband+1))/2,1,1] ; nc_stride_5=[1,1,1,1,1] 
           NCF_CHECK(nf90_get_var(ncid,varid,psinablapsi,start=nc_start_5,stride=nc_stride_5,count=nc_count_5))
         end if
#endif
       else
         psinablapsi=zero
         bsize=nband_k**2;if (read_half_dipoles) bsize=(nband_k*(nband_k+1))/2
         read(opt_unt)(psinablapsi(1:2,1,ijband),ijband=1,bsize)
         read(opt_unt)(psinablapsi(1:2,2,ijband),ijband=1,bsize)
         read(opt_unt)(psinablapsi(1:2,3,ijband),ijband=1,bsize)           
       end if
       if (.not.mykpt) then
         call xmpi_exch(psinablapsi,etiq,master,psinablapsi,master_band,comm,ierr)
       end if
     end if

!    Select k-points for current proc
     if (mykpt) then

       ABI_MALLOC(eig0_k,(nband_k))
       ABI_MALLOC(occ_k,(nband_k))
       ABI_MALLOC(doccde_k,(nband_k))

       cond_nd_k = zero
       kin11_k   = zero
       kin12_k   = zero
       kin21_k   = zero
       kin22_k   = zero
       np_sum_k1 = zero
       np_sum_k2 = zero
       socc_k    = zero

!      k-dependent data
       eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
       occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
       doccde_k(:)=doccde(1+bdtot_index:nband_k+bdtot_index)

!      In case of non MPI-IO, receive all (n,m) dipoles from master proc
!        Then broadcast them to all band processors
       if (.not.iomode_estf_mpiio) then
         if (me/=master.and.me==master_band) then
           call xmpi_exch(psinablapsi,etiq,master,psinablapsi,me,comm,ierr)
         end if
         call xmpi_bcast(psinablapsi,master,mpi_enreg%comm_band,mpierr)
       end if  

!      LOOP OVER BANDS n
       do iband=1,nband_k

         !If MPI-IO, store only ib elements for each jb
         !If not, store all (ib,jb) pairs
         my_iband=merge(1,iband,iomode_estf_mpiio)

!        Select bands for current proc
         myband=(mpi_enreg%proc_distrb(ikpt,iband,isppol)==me)
         if (myband) then
       
!          In case of MPI-IO, read valence-valence dipoles for band n
           if (iomode_estf_mpiio) then
#ifdef HAVE_NETCDF
             if (nc_unlimited) then
               nc_start_6=[1,1,iband,ikpt,isppol,1] ; nc_count_6=[2,3,1,1,1,mband] ; nc_stride_6=[1,1,1,1,1,1] 
               NCF_CHECK(nf90_get_var(ncid,varid,psinablapsi,start=nc_start_6,stride=nc_stride_6,count=nc_count_6))
             else if (.not.read_half_dipoles) then
               nc_start_6=[1,1,1,iband,ikpt,isppol] ; nc_count_6=[2,3,mband,1,1,1] ; nc_stride_6=[1,1,1,1,1,1] 
               NCF_CHECK(nf90_get_var(ncid,varid,psinablapsi,start=nc_start_6,stride=nc_stride_6,count=nc_count_6))
             else
               nc_start_5=[1,1,(iband*(iband-1))/2+1,ikpt,isppol] ; nc_count_5=[2,3,iband,1,1] ; nc_stride_5=[1,1,1,1,1] 
               NCF_CHECK(nf90_get_var(ncid,varid,psinablapsi,start=nc_start_5,stride=nc_stride_5,count=nc_count_5))
             end if
#endif
           end if

!          LOOP OVER BANDS m
           do jband=1,iband-1
             diff_occ = occ_k(iband)-occ_k(jband)
             diff_eig = eig0_k(iband)-eig0_k(jband)
             if (dabs(diff_occ)>=tol8) then

               dhdk2_r = zero
               dhdk2_g = zero

               if (read_half_dipoles) then
                 ijband=(my_iband*(my_iband-1))/2+jband
               else
                 !psinablapsi size is mband for netCDF I/O, nband_k for Fortran I/O
                 bd_stride=merge(mband,nband_k,iomode==IO_MODE_ETSF)
                 ijband=(my_iband-1)*bd_stride+jband
               end if
                 
               do l2=1,3
                 do l1=1,3
                   dhdk2_r(l1,l2)=dhdk2_r(l1,l2)+(&
&                    psinablapsi(1,l1,ijband)*psinablapsi(1,l2,ijband)&
&                   +psinablapsi(2,l1,ijband)*psinablapsi(2,l2,ijband))
                 end do
               end do
               do l1=1,3
                 dhdk2_g=dhdk2_g &
&                  +(psinablapsi(1,l1,ijband)*psinablapsi(1,l1,ijband) &
&                   +psinablapsi(2,l1,ijband)*psinablapsi(2,l1,ijband))
               end do


               !Minimal validity limit
               deltae_min_tmp=dabs(diff_eig)
               if ((deltae_min_tmp>=tol5).and.(deltae_min_tmp<=deltae_min)) deltae_min=deltae_min_tmp

               !Conductivity for each omega - Apply KG formula
               kin_fact=(eig0_k(iband)+eig0_k(jband))*half-(fermie+entropy)
               do iom=1,mom
                 oml=oml1(iom)
                 dirac=(dexp(-((diff_eig+oml)/dom)**2)-dexp(-((diff_eig-oml)/dom)**2))/oml ! Take into account (n,m) and (m,n)
                 sig=dhdk2_g*diff_occ*dirac
                 kin11_k(iom)=kin11_k(iom)+sig
                 kin12_k(iom)=kin12_k(iom)-sig*kin_fact
                 kin21_k(iom)=kin21_k(iom)-sig*kin_fact
                 kin22_k(iom)=kin22_k(iom)+sig*kin_fact**2
                 do l2=1,3
                   do l1=1,3
                     cond_nd_k(l1,l2,iom)=cond_nd_k(l1,l2,iom)+dhdk2_r(l1,l2)*diff_occ*dirac
                   end do
                 end do
               end do

               !Evaluate sumrule
               fact_diag=merge(one,two,iband==jband)
               if (dabs(diff_eig)>=tol10) then
                 np_sum_k1=np_sum_k1 - fact_diag*dhdk2_g*diff_occ/diff_eig
               else
                 np_sum_k2=np_sum_k2 - fact_diag*doccde_k(iband)*dhdk2_g
               end if

             end if ! diff_occ>tol8
           end do !jband

           socc_k=socc_k+occ_k(iband)

         end if ! my band?
       end do ! iband

!      Accumulate k-point contribution
       do iom=1,mom
         kin11(iom,isppol)=kin11(iom,isppol)+wtk(ikpt)*kin11_k(iom)
         kin12(iom)=kin12(iom)+wtk(ikpt)*kin12_k(iom)
         kin21(iom)=kin21(iom)+wtk(ikpt)*kin21_k(iom)
         kin22(iom)=kin22(iom)+wtk(ikpt)*kin22_k(iom)
         cond_nd(:,:,iom)=cond_nd(:,:,iom)+wtk(ikpt)*cond_nd_k(:,:,iom)
       end do
       np_sum=np_sum+wtk(ikpt)*(np_sum_k1+np_sum_k2)
       socc=socc+wtk(ikpt)*socc_k

!      Validity limit
       deltae=deltae+(eig0_k(nband_k)-fermie)

       ABI_FREE(eig0_k)
       ABI_FREE(occ_k)
       ABI_FREE(doccde_k)

!    End loop over kpt/spin
     end if ! My kpt?
     bdtot_index=bdtot_index+nband_k
   end do ! ikpt
 end do ! isppol

!Accumulate kpt/band contributions over processors
 call xmpi_sum(kin11,comm,mpierr)
 call xmpi_sum(kin12,comm,mpierr)
 call xmpi_sum(kin21,comm,mpierr)
 call xmpi_sum(kin22,comm,mpierr)
 call xmpi_sum(np_sum,comm,mpierr)
 call xmpi_sum(socc,comm,mpierr)
 call xmpi_sum(deltae,comm,mpierr)
 call xmpi_min(deltae_min,comm,mpierr)
 deltae=deltae/mpi_enreg%nproc_band
 
 ABI_FREE(psinablapsi)

!---------------------------------------------------------------------------------
! Output results

!Print file headers (only master node)
 if (me==master) then

!  Standard output
   write(std_out,'(a,3f10.5)')' sumrule           =',np_sum/socc/three/dble(nsppol),socc
   write(std_out,'(a,f10.5,a,f10.5,a)')&
&   ' Emax-Efermi       =',deltae/dble(nkpt*nsppol),' Ha', &
&                          deltae/dble(nkpt*nsppol)*Ha_eV,' eV'
   write(std_out,'(a,f10.5,a,f10.5,a)')&
&   ' DeltaE min        =',deltae_min,' Ha',deltae_min*Ha_eV,' eV'

!  _Lij file
   if (open_file(trim(filnam_out)//'_Lij',msg, newunit=lij_unt, form='formatted', action="write") /= 0) then
     ABI_ERROR(msg)
   end if
   write(lij_unt,'(a)')' # omega(ua) L12 L22 L22'

!  _sig file
   if (open_file(trim(filnam_out)//'_sig', msg, newunit=sig_unt, form='formatted', action="write") /= 0) then
     ABI_ERROR(msg)
   end if
   if (nsppol==1) then
     write(sig_unt,'(a)')' # omega(ua) hbar*omega(eV)    cond(ua)             cond(ohm.cm)-1'
   else
     write(sig_unt,'(2a)')' # omega(ua) hbar*omega(eV)      cond(ua)            cond(ohm.cm)-1',&
&     '      cond(ohm.cm)-1 UP      cond(ohm.cm)-1 DN'
   end if

!  _Kth file
   if (open_file(trim(filnam_out)//'_Kth', msg, newunit=kth_unt, form='formatted', action="write") /=0) then
     ABI_ERROR(msg)
   end if
   write(kth_unt,'(a)')&
&   " #Thermal conductivity following B. Holst et al Phys. Rev. B 83 (2011) 235120"
   write(kth_unt,'(a)')&
&   ' #omega(ua) hbar*omega(eV)  thermal cond(ua)   Kth(W/m/K)   thermopower(ua)   Stp(microohm/K)'

!  Output file
   if (open_file(trim(filnam_out)//'.out', msg, newunit=ocond_unt, form='formatted', action="write") /= 0) then
     ABI_ERROR(msg)
   end if
   write(ocond_unt,'(a)' )'#Conducti output file:'
   write(ocond_unt,'(a)' )'#Contains all results produced by conducti utility'
   write(ocond_unt,'(a)' )'#  '
   write(ocond_unt,'(a,i8,3f10.5,a)')'# npts,omin,omax,width     =' ,mom,omin,omax,dom,' Ha'
   write(ocond_unt,'(a,3f10.5,a)' )'# rprimd(bohr)      =',rprimd(1:3,1)
   write(ocond_unt,'(a,3f10.5,a)' )'#                    ',rprimd(1:3,2)
   write(ocond_unt,'(a,3f10.5,a)' )'#                    ',rprimd(1:3,3)
   write(ocond_unt,'(a,i8)' )      '# natom             =',natom
   write(ocond_unt,'(a,3i8)' )     '# nkpt,mband,nsppol        =',nkpt,mband,nsppol
   write(ocond_unt,'(a, f10.5,a)' )'# ecut                  =',ecut,' Ha'
   write(ocond_unt,'(a,f10.5,a,f10.5,a)' )'# fermie            =',fermie,' Ha ',fermie*Ha_eV,' eV'

   write(ocond_unt,'(a,f12.5,a,f12.5,a)') '# Temp              =',tsmear,' Ha ',Tatm,' Kelvin'
   write(ocond_unt,'(a,f15.5,2x,f15.5)' )'# sumrule           =',np_sum/socc/three/dble(nsppol),socc
   write(ocond_unt,'(a,f10.5,a,f10.5,a)' )&
&   '# Emax-Efermi       =',deltae/dble(nkpt*nsppol),' Ha',deltae/dble(nkpt*nsppol)*Ha_eV,' eV'
   write(ocond_unt,'(a)' )'# '
   write(ocond_unt,'(a)')'# omega(ua)       cond(ua)             thermal cond(ua)       thermopower(ua)'

 end if ! me==master?

!Compute (and print) thermal conductivity and thermopower
 do iom=1,mom
   oml=oml1(iom)

   do isppol=1,nsppol
     kin11(iom,isppol)=kin11(iom,isppol)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
     if (dabs(kin11(iom,isppol))<tol19) kin11(iom,isppol)=zero
     sig_abs(iom)=sig_abs(iom)+kin11(iom,isppol)
   end do

   kin21(iom)=kin21(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   kin12(iom)=kin12(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   kin22(iom)=kin22(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)

   Kth(iom)=kin22(iom)
   Stp(iom)=zero
   if (sig_abs(iom)/=zero)  then
     Kth(iom)=Kth(iom)-(kin12(iom)*kin21(iom)/sig_abs(iom))
     Stp(iom)=kin12(iom)/(sig_abs(iom)*tsmear)
   end if

   if (dabs(Kth(iom))<tol19) Kth(iom)=zero
   if (dabs(Stp(iom))<tol19) Stp(iom)=zero
   if (abs(kin12(iom))<1.d-80) kin12=zero
   if (abs(kin21(iom))<1.d-80) kin21=zero
   if (abs(kin22(iom))<1.d-80) kin22=zero

   if (me==master) then
     write(lij_unt,'(f12.5,4es22.12)')oml,kin12(iom),kin22(iom),kin22(iom)/Tatm*3.4057d9
     if (nsppol==1) then
       write(sig_unt,'(2f12.5,2es22.12)') oml,oml*Ha_eV,sig_abs(iom),sig_abs(iom)*Ohmcm
     else
       write(sig_unt,'(2f12.5,4es22.12)') oml,oml*Ha_eV,sig_abs(iom),sig_abs(iom)*Ohmcm, &
&        kin11(iom,1)*Ohmcm,kin11(iom,2)*Ohmcm
     end if
     write(kth_unt,'(2f12.5,4es22.12)') oml,oml*Ha_eV,Kth(iom),Kth(iom)*3.4057d9/Tatm, &
&     Stp(iom),Stp(iom)*3.6753d-2
     write(ocond_unt,'(1f12.5,3es22.12)') oml,sig_abs(iom),Kth(iom),Stp(iom)
   end if

 end do

!Compute the imaginary part of the conductivity (principal value)
!  +derived optical properties.
 if (me==master) then
   call msig(sig_abs,mom,oml1,filnam_out)
 end if
 
!---------------------------------------------------------------------------------
! End

!Close all files
 if (me==master) then
   write(std_out,'(2a)')ch10,'OUTPUT'
   write(std_out,'(a)')trim(filnam_out)//'_Lij : Onsager kinetic coefficients'
   write(std_out,'(a)')trim(filnam_out)//'_sig : Optical conductivity'
   write(std_out,'(a)')trim(filnam_out)//'_Kth : Thermal conductivity and thermopower'
   write(std_out,'(a)')trim(filnam_out)//'_eps : Dielectric function'
   write(std_out,'(a)')trim(filnam_out)//'_abs : n, k, reflectivity, absorption'
   close(lij_unt)
   close(sig_unt)
   close(kth_unt)
   close(ocond_unt)
 end if
 if (iomode == IO_MODE_ETSF) then
   if (iomode_estf_mpiio.or.me==master) then
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_close(ncid))
#endif
   end if
 else if (me==master) then
   ierr=close_unit(opt_unt,msg)
   ABI_CHECK(ierr==0,sjoin("Error while closing ",filnam1))
 end if

!Release memory space
 ABI_FREE(kin11)
 ABI_FREE(kin22)
 ABI_FREE(kin12)
 ABI_FREE(kin21)
 ABI_FREE(kin11_k)
 ABI_FREE(kin22_k)
 ABI_FREE(kin12_k)
 ABI_FREE(kin21_k)
 ABI_FREE(Stp)
 ABI_FREE(Kth)
 ABI_FREE(cond_nd)
 ABI_FREE(cond_nd_k)
 ABI_FREE(sig_abs)
 ABI_FREE(eigen0)
 ABI_FREE(nband)
 ABI_FREE(oml1)
 ABI_FREE(occ)
 ABI_FREE(doccde)
 ABI_FREE(wtk)
 call hdr%free()
 call destroy_mpi_enreg(mpi_enreg)

end subroutine conducti_paw
!!***

!----------------------------------------------------------------------

!!****f* m_conducti/conducti_paw_core
!! NAME
!! conducti_paw_core
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula for PAW formalism
!!
!! INPUTS
!!  filnam=generic name for input data
!!  filnam_out=generic name for output data
!!  [with_absorption]=optiona lflag to activate the computation of absorption (_sigX file) (default=TRUE)
!!  [with_emissivity]=optiona lflag to activate the computation of emissivity (_emisX file) (default=FALSE)
!!
!! OUTPUT
!!  Only printing
!!
!! NOTES
!!  bantot
!!  dom=frequency range
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree).
!!  ecut=kinetic energy planewave cutoff (hartree).
!!  fermie= fermi energy (Hartree)
!!  mom=number of frequency for conductivity computation
!!  mband=maximum number of bands.
!!  natom = number of atoms in the unit cell.
!!  nband(nkpt*nsppol)=number of bands at each RF k point for each spin.
!!  nkpt=number of k points in the IBZ for this perturbation
!!  ngfft(3)=integer fft box dimensions.
!!  nspinor=number of spinorial components of the wavefunctions.
!!  nsppol=1 for unpolarized, 2 for spin-polarized.
!!  ntypat = number of atom types.
!!  occ(mband*nkpt*nsppol)=occupation number for each band and k.
!!  occopt==option for occupancies
!!  psinablapsi2(2,3,mband,nphicor,natom)Matrix elements = <Phi_core|Nabla|Phi_i>
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).sigx(mom,nphicor))
!!  rprimd(3,3)=real space primitive translations.
!!  of primitive translations.
!!  ucvol=unit cell volume in ($\textrm{bohr}^{3}$).
!!  wind=frequency window for computations of sigma
!!  wtk(nkpt)=weight assigned to each k point.
!!
!! PARENTS
!!      conducti
!!
!! CHILDREN
!!      close_unit,getnel,hdr_fort_read,hdr_ncread,hdr%bcast,hdr%free
!!      metric,msig,nctk_fort_or_ncfile,nctk_get_dim,nctk_open_read,nctk_set_collective
!!      nf90_close,nf90_get_var,open_file,spline,splint,xmpi_bcast,xmpi_sum
!!
!! SOURCE

 subroutine conducti_paw_core(filnam,filnam_out,with_absorption,with_emissivity)

!Arguments -----------------------------------
!scalars
 character(len=fnlen) :: filnam,filnam_out
 logical,intent(in),optional :: with_absorption,with_emissivity

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iomode,atnbr,bantot,bdtot_index,comm,my_iband
 integer :: fform2,headform,iatom,iband,icor,ierr,ikpt
 integer :: iom,isppol,l1,mband,me,mom,mpierr,j2,etiq
 integer :: natom,nband_k,nkpt,nphicor,nproc,nspinor,nsppol,ntypat
 integer :: occopt,iunt,opt2_unt,ncid,varid,master_band,nb_per_proc
 integer :: sigx1_unt,sigx1_s_unt,sigx1_tot_unt,ems_unt,ems_s_unt,ems_tot_unt
 logical :: iomode_estf_mpiio,myband,mykpt,need_absorption,need_emissivity
 real(dp) :: del_sig,del_emis,deltae,diff_occ,ecut,fermie
 real(dp) :: omin,omax,omin_sig,omax_sig,omin_emis,omax_emis
 real(dp) :: oml,dom,dom_ctr,dom_max,dom_tan1,dom_tan2
 real(dp) :: Tatm,tsmear,ucvol
 character(len=fnlen) :: filnam2,filnam_gen
 character(len=500) :: msg
 type(hdr_type) :: hdr
 type(MPI_type) :: mpi_enreg
!arrays
 integer :: nc_count(7),nc_start(7),nc_stride(7)
 integer,allocatable :: nband(:),ncor(:),lcor(:),kappacor(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: dom_var1(:,:),dhdk2_g(:)
 real(dp),allocatable :: eig0_k(:),eigen0(:),eig0nc(:,:,:)
 real(dp),allocatable :: energy_cor(:),edge(:)
 real(dp),allocatable :: occ(:),occ_k(:),wtk(:)
 real(dp),allocatable :: oml_edge(:,:),oml_emis(:)
 real(dp),allocatable :: psinablapsi2(:,:,:,:,:)
 real(dp),allocatable :: sigx1(:,:,:,:),sigx1_av(:,:,:),sigx1_k(:,:,:)
 real(dp),allocatable :: sum_spin_sigx1(:,:,:),sum_spin_sigx1_av(:,:)
 real(dp),allocatable :: emisx(:,:,:,:),emisx_av(:,:,:),emisx_k(:,:,:)
 real(dp),allocatable :: sum_spin_emisx(:,:,:),sum_spin_emisx_av(:,:)

! *********************************************************************************

!optional flags
 need_absorption=.true. ;if (present(with_absorption)) need_absorption=with_absorption
 need_emissivity=.false.;if (present(with_emissivity)) need_emissivity=with_emissivity
 if ((.not.need_absorption).and.(.not.need_emissivity)) return

! ---------------------------------------------------------------------------------
! Read input data

!Global MPI communicators
 comm = xmpi_world
 nproc = xmpi_comm_size(comm)
 me = xmpi_comm_rank(comm)

!Read input parameters file
 if (me==master) then
   if (open_file(filnam,msg,newunit=iunt,form='formatted',action="read",status="old")/=0) then
     ABI_ERROR(msg)
   end if
   rewind(iunt)
   read(iunt,*)
   read(iunt,'(a)')filnam_gen ! Generic name for the files
   filnam2=trim(filnam_gen)//'_OPT2'
!  Read frequency range
   if (need_absorption) then
     read(iunt,err=11,end=11,fmt=*) dom,omin,omax,mom,atnbr,dom_max,dom_ctr
     goto 12
11   backspace(iunt) ; read(iunt,*) dom,omin,omax,mom,atnbr ; dom_max=zero ; dom_ctr=zero
12   continue
else if (need_emissivity) then
     read(iunt,*) dom,omin,omax,mom,atnbr
     dom_max=zero;dom_ctr=zero
   end if
   close(iunt)
   if (abs(dom_max)>tol10.and.dom_max<dom) then
     msg = 'dom_max must be higher than dom!'
     ABI_ERROR(msg)
   end if
 end if
   
!Send data to all procs
 call xmpi_bcast(dom,master,comm,mpierr)
 call xmpi_bcast(omin,master,comm,mpierr)
 call xmpi_bcast(omax,master,comm,mpierr)
 call xmpi_bcast(mom,master,comm,mpierr)
 call xmpi_bcast(atnbr,master,comm,mpierr)
 call xmpi_bcast(dom_max,master,comm,mpierr)
 call xmpi_bcast(dom_ctr,master,comm,mpierr)

! ---------------------------------------------------------------------------------
! Read OPT2 file

!Check for FORTRAN/.nc OPT2 file and set iomode to IO_MODE_FORTRAN_MASTER/IO_MODE_ETSF
 if (me==master) then
   call nctk_fort_or_ncfile(filnam2,iomode,msg)
   if (iomode/=IO_MODE_ETSF) iomode=IO_MODE_FORTRAN_MASTER
 end if
 call xmpi_bcast(filnam2,master,comm,mpierr)
 call xmpi_bcast(iomode,master,comm,mpierr)

!Open OPT2 file and read HEADER
 if (me==master) then
   if (iomode==IO_MODE_ETSF) then
     NCF_CHECK(nctk_open_read(ncid,filnam2,xmpi_comm_self))
     call hdr_ncread(hdr,ncid,fform2)
   else
     if (open_file(filnam2,msg,newunit=opt2_unt,form="unformatted",status="old")/=0) then
       ABI_ERROR(msg)
     end if
     call hdr_fort_read(hdr,opt2_unt,fform2,rewind=.true.)
   end if 
   ABI_CHECK(fform2/=0,sjoin("Error while reading ",filnam2))
   ABI_CHECK(fform2==611.or.fform2==612.or.fform2==613,"OPT2 file format should be fform=611/612/613!")
 end if
 call hdr%bcast(master,me,comm)
 call xmpi_bcast(fform2,master,comm,mpierr)

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
 fermie=hdr%fermie
 tsmear=hdr%tsmear
 ABI_MALLOC(nband,(nkpt*nsppol))
 ABI_MALLOC(occ,(bantot))
 ABI_MALLOC(wtk,(nkpt))
 occ(1:bantot)=hdr%occ(1:bantot)
 wtk(1:nkpt)=hdr%wtk(1:nkpt)
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)
 mband=maxval(nband(:))  ! Get mband, as the maximum value of nband(nkpt)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol) ! Get metrics of simulation cell

!Read eigenvalues
 ABI_MALLOC(eigen0,(mband*nkpt*nsppol))
 if (me==master) then
   if (iomode==IO_MODE_ETSF) then
     varid=nctk_idname(ncid,"eigenvalues")
     ABI_MALLOC(eig0nc,(mband,nkpt,nsppol))
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_get_var(ncid,varid,eig0nc))
#endif
     eigen0 = reshape(eig0nc,[mband*nkpt*nsppol])
     ABI_FREE(eig0nc)
   else
     read(opt2_unt)(eigen0(iband),iband=1,mband*nkpt*nsppol)
   end if
 end if
 call xmpi_bcast(eigen0,master,comm,mpierr)

!Read core states
 if (me==master) then
   if (iomode==IO_MODE_ETSF) then
#ifdef HAVE_NETCDF
     NCF_CHECK(nctk_get_dim(ncid,"number_of_core_states",nphicor))
#endif
   else
     read(opt2_unt) nphicor
   end if
  end if
 call xmpi_bcast(nphicor,master,comm,mpierr)

 ABI_MALLOC(ncor,(nphicor))
 ABI_MALLOC(lcor,(nphicor))
 ABI_MALLOC(kappacor,(nphicor))
 ABI_MALLOC(energy_cor,(nphicor))
 ABI_MALLOC(edge,(nphicor))

 if (me==master) then
   if (iomode==IO_MODE_ETSF) then
#ifdef HAVE_NETCDF
     varid=nctk_idname(ncid,"n_quantum_number_core")
     NCF_CHECK(nf90_get_var(ncid,varid,ncor))
     varid=nctk_idname(ncid,"l_quantum_number_core")
     NCF_CHECK(nf90_get_var(ncid,varid,lcor))
     varid=nctk_idname(ncid,"kappa_core")
     NCF_CHECK(nf90_get_var(ncid,varid,kappacor))
     varid=nctk_idname(ncid,"eigenvalues_core")
     NCF_CHECK(nf90_get_var(ncid,varid,energy_cor))
!Close here netcdf file here because the rest has to be read with collective I/O
     NCF_CHECK(nf90_close(ncid))
#endif
   else
     do icor=1,nphicor
       read(unit=opt2_unt,err=23,end=23) ncor(icor),lcor(icor),kappacor(icor),energy_cor(icor)
       goto 24
23     backspace(opt2_unt) ; read(unit=opt2_unt) ncor(icor),lcor(icor),energy_cor(icor)
       kappacor(icor)=0
24     continue
     end do
   end if
 end if ! master
 call xmpi_bcast(nphicor,master,comm,mpierr)
 call xmpi_bcast(ncor,master,comm,mpierr)
 call xmpi_bcast(lcor,master,comm,mpierr)
 call xmpi_bcast(kappacor,master,comm,mpierr)
 call xmpi_bcast(energy_cor,master,comm,mpierr)
 edge(1:nphicor)=fermie-energy_cor(1:nphicor)

!---------------------------------------------------------------------------------
! Prepare kpt/band parallelization

 call init_mpi_enreg(mpi_enreg)
 mpi_enreg%comm_kpt=comm
 mpi_enreg%me_kpt=me
 mpi_enreg%nproc_spkpt=nproc
 mpi_enreg%paralbd=1
 ABI_MALLOC(mpi_enreg%proc_distrb,(nkpt,mband,nsppol))
 ABI_MALLOC(mpi_enreg%my_kpttab,(nkpt))
 call distrb2(mband,nb_per_proc,nband,nkpt,nproc,nsppol,mpi_enreg)
 call initmpi_band(nkpt,mpi_enreg,nband,nkpt,nsppol)

!---------------------------------------------------------------------------------
!Print some data

 Tatm=tsmear*Ha_K
 if (me==master) then
   write(std_out,*)
   write(std_out,'(a)')'--------------------------------------------'
   write(std_out,'(a,i4)') 'selected atom for X ray emission',atnbr
   write(std_out,'(a)')'--------------------------------------------'
   if (need_absorption) then
     if (abs(dom_max)>tol10) then
       write(std_out,'(a)')'************************ ARCTAN SMEARING'
       write(std_out,'(a,i8,3f10.5,a)')' npts,omin,omax,width      =',mom,omin,omax,dom,' Ha'
       write(std_out,'(a,2f10.5,a)')' dom_max,center       =',dom_max,dom_ctr,' Ha'
     else
       write(std_out,'(a)')'************************ FIXED SMEARING'
       write(std_out,'(a,i8,3f10.5,a)')' npts,omin,omax,width      =',mom,omin,omax,dom,' Ha'
     endif
   end if
   write(std_out,*)
   write(std_out,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1,1:3)
   write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(2,1:3)
   write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(3,1:3)
   write(std_out,'(a,i8)')       ' natom             =',natom
   write(std_out,'(a,3i8)')      ' nkpt,mband,nsppol        =',nkpt,mband,nsppol
   write(std_out, '(a, f10.5,a)' ) ' ecut              =',ecut,' Ha'
   write(std_out,'(a,f10.5,a,f10.5,a)' )' fermie            =',fermie,' Ha',fermie*Ha_eV,' eV'
   write(std_out,'(a,f12.5,a,f12.5,a)') ' Temp              =',tsmear,' Ha ',Tatm,' Kelvin'
   write(std_out,*)
   write(std_out,*)
   write(std_out,'(a)')'--------------------------------------------'
   write(std_out,'(a,i4)') ' Number of core orbitals nc=',nphicor
   do icor=1,nphicor
     if (kappacor(icor)==0) then
       write(std_out,'(a,2i4,2f10.5)') ' n, l, Energy(Ha), Edge(Ha): ', &
         ncor(icor),lcor(icor),energy_cor(icor),edge(icor)
     else
       if (kappacor(icor)>0) then
         j2=2*lcor(icor)-1
         write(std_out,'(a,i4,i4,a,i4,2f10.5)') ' n, j, l, Energy(Ha), Edge(Ha): ',&
&          ncor(icor),j2,' / 2',lcor(icor),energy_cor(icor),edge(icor)
       else
         if (kappacor(icor)<-1) then
           j2=2*lcor(icor)+1
           write(std_out,'(a,i4,i4,a,i4,2f10.5)') ' n, j, l, Energy(Ha), Edge(Ha): ',&
&            ncor(icor),j2,'/2',lcor(icor),energy_cor(icor),edge(icor)
         else
           write(std_out,'(a,i4,a,i4,2f10.5)') ' n, j, l, Energy(Ha), Edge(Ha): ',&
&            ncor(icor),'   1/2',lcor(icor),energy_cor(icor),edge(icor)
         end if
       end if
     end if
   end do
   write(std_out,'(a)')'--------------------------------------------'
 end if ! master

!---------------------------------------------------------------------------------
! Determine the frequency range and allocate frequency-dependent arrays

 if (need_absorption) then
   ABI_MALLOC(oml_edge,(nphicor,mom))
   ABI_MALLOC(dom_var1,(nphicor,mom))
   omax_sig=omax ; omin_sig=omin
   del_sig=(omax_sig-omin_sig)/(mom-1)
   do iom=1,mom
     do icor=1,nphicor
       oml_edge(icor,iom)=-energy_cor(icor)+ dble(iom-1)*del_sig + omin_sig - one
       if ((oml_edge(icor,iom)<=edge(icor)).or.(abs(dom_max)<=tol10)) then
         dom_var1(icor,iom)= dom
       else
         dom_tan1= (oml_edge(icor,iom)-edge(icor))/dom_ctr
         dom_tan2=dom_tan1-one/dom_tan1**2
         dom_var1(icor,iom)=dom+dom_max*(half+piinv*datan(dom_tan2))
       endif
     enddo
   enddo
   ABI_MALLOC(sigx1,(nphicor,mom,natom,nsppol))
   sigx1=zero
 end if
 
 if (need_emissivity) then
   ABI_MALLOC(oml_emis,(mom))
   omin_emis=minval(eigen0)
   omax_emis=maxval(eigen0)
   del_emis=(omax_emis-omin_emis)/(mom-1)
   do iom=1,mom
     oml_emis(iom)=omin_emis+dble(iom-1)*del_emis
   end do
   ABI_MALLOC(emisx,(nphicor,mom,natom,nsppol))
   emisx=zero
 end if

!---------------------------------------------------------------------------------
! Prepare core-valence dipoles reading

 iomode_estf_mpiio=(iomode==IO_MODE_ETSF.and.nctk_has_mpiio.and.use_netcdf_mpiio)

 !In case of netCDF access to OPT file, prepare collective I/O
 if (iomode == IO_MODE_ETSF) then
   if (iomode_estf_mpiio) then
     NCF_CHECK(nctk_open_read(ncid,filnam2,comm))
     varid=nctk_idname(ncid,"dipole_core_valence")
     if (nproc>1) then
       NCF_CHECK(nctk_set_collective(ncid,varid))
     end if
   else if (me==master) then
     NCF_CHECK(nctk_open_read(ncid,filnam2,xmpi_comm_self))
     varid=nctk_idname(ncid,"dipole_core_valence")
     !if (nctk_has_mpiio.and.(.not.use_netcdf_mpiio)) then
     !  NCF_CHECK(nctk_set_collective(ncid,varid))
     !end if
   end if
 end if

 if (iomode_estf_mpiio) then
   !If MPI-IO, store only elements for one band
   ABI_MALLOC(psinablapsi2,(2,3,nphicor,natom,1))
 else
   !If not, store the elements for all bands
   ABI_MALLOC(psinablapsi2,(2,3,nphicor,natom,mband))
 end if

!---------------------------------------------------------------------------------
! Compute X absorption coefficient and/or X emissivity

 if (need_absorption) then
   ABI_MALLOC(sigx1_k,(nphicor,mom,natom))
 end if
 if (need_emissivity) then
   ABI_MALLOC(emisx_k,(nphicor,mom,natom))
 end if
 ABI_MALLOC(dhdk2_g,(nphicor))

 deltae  = zero

!LOOP OVER SPINS/K
 bdtot_index = 0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     etiq=ikpt+(isppol-1)*nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     mykpt=.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me))
     master_band=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol))

!    In case of non MPI-IO, has to read all (n,m) dipoles for this k-point
!      Master node reads and send to relevant processor
     if (.not.iomode_estf_mpiio.and.me==master) then
       if (iomode==IO_MODE_ETSF) then
         nc_start=[1,1,1,1,1,ikpt,isppol];nc_stride=[1,1,1,1,1,1,1] 
         nc_count=[2,3,nphicor,natom,mband,1,1]
#ifdef HAVE_NETCDF
         NCF_CHECK(nf90_get_var(ncid,varid,psinablapsi2,start=nc_start,stride=nc_stride,count=nc_count))
#endif
       else
         psinablapsi2=zero
         if (fform2==612) then ! New OPT2 file format
           read(opt2_unt) (((psinablapsi2(1:2,1,icor,iatom,iband),icor=1,nphicor),iatom=1,natom),iband=1,nband_k)
           read(opt2_unt) (((psinablapsi2(1:2,2,icor,iatom,iband),icor=1,nphicor),iatom=1,natom),iband=1,nband_k)
           read(opt2_unt) (((psinablapsi2(1:2,3,icor,iatom,iband),icor=1,nphicor),iatom=1,natom),iband=1,nband_k)
         else if (fform2==613) then ! Large OPT2 file format
           do iband=1,nband_k
             read(opt2_unt) ((psinablapsi2(1:2,1,icor,iatom,iband),icor=1,nphicor),iatom=1,natom)
             read(opt2_unt) ((psinablapsi2(1:2,2,icor,iatom,iband),icor=1,nphicor),iatom=1,natom)
             read(opt2_unt) ((psinablapsi2(1:2,3,icor,iatom,iband),icor=1,nphicor),iatom=1,natom)
           end do
         else
           !The old writing was not efficient (indexes order is bad)
           do iatom=1,natom
             read(opt2_unt) ((psinablapsi2(1:2,1,icor,iatom,iband),iband=1,nband_k),icor=1,nphicor)
             read(opt2_unt) ((psinablapsi2(1:2,2,icor,iatom,iband),iband=1,nband_k),icor=1,nphicor)
             read(opt2_unt) ((psinablapsi2(1:2,3,icor,iatom,iband),iband=1,nband_k),icor=1,nphicor)
           end do
         end if
       end if
       if (.not.mykpt) then
         call xmpi_exch(psinablapsi2,etiq,master,psinablapsi2,master_band,comm,ierr)
       end if
     end if

!    Select k-points for current proc
     if (mykpt) then

       ABI_MALLOC(eig0_k,(nband_k))
       ABI_MALLOC(occ_k,(nband_k))
       
       if (need_absorption) sigx1_k=zero
       if (need_emissivity) emisx_k=zero

!      k-dependent data
       eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
       occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)

!      In case of non MPI-IO, receive all (n,m) dipoles from master proc
!        Then broadcast them to all band processors
       if (.not.iomode_estf_mpiio) then
         if (me/=master.and.me==master_band) then
           call xmpi_exch(psinablapsi2,etiq,master,psinablapsi2,me,comm,ierr)
         end if
         call xmpi_bcast(psinablapsi2,master,mpi_enreg%comm_band,mpierr)
       end if  

!      LOOP OVER BANDS
       do iband=1,nband_k

         !If MPI-IO, store only ib elements for each iband
         !If not, store all elements
         my_iband=merge(1,iband,iomode_estf_mpiio)
         dhdk2_g   = zero

!        Select bands for current proc
         myband=(mpi_enreg%proc_distrb(ikpt,iband,isppol)==me)
         if (myband) then

           diff_occ = (two/dble(nsppol*nspinor))-occ_k(iband)

!          In case of MPI-IO, read core-valence dipoles for band n
           if (iomode_estf_mpiio) then
             nc_start=[1,1,1,1,iband,ikpt,isppol];nc_stride=[1,1,1,1,1,1,1] 
             nc_count=[2,3,nphicor,natom,1,1,1]
#ifdef HAVE_NETCDF
             NCF_CHECK(nf90_get_var(ncid,varid,psinablapsi2,start=nc_start,stride=nc_stride,count=nc_count))
#endif
           end if

!          LOOP OVER ATOMS
           do iatom=1,natom

             dhdk2_g = zero
             do icor=1,nphicor
               do l1=1,3
                 dhdk2_g(icor)=dhdk2_g(icor) &
&                 +(psinablapsi2(1,l1,icor,iatom,my_iband)*psinablapsi2(1,l1,icor,iatom,my_iband) &
&                  +psinablapsi2(2,l1,icor,iatom,my_iband)*psinablapsi2(2,l1,icor,iatom,my_iband))
               end do
             end do

!            Absoption for each omega
             if (need_absorption) then
               do iom=1,mom
                 do icor=1,nphicor
                   !Lorentzian var arctan
                   sigx1_k(icor,iom,iatom)=sigx1_k(icor,iom,iatom) &
&                     +dhdk2_g(icor)*diff_occ*oml_edge(icor,iom) &
&                     *(dom_var1(icor,iom)/((-energy_cor(icor)+eig0_k(iband) &
&                      -oml_edge(icor,iom))**2+(dom_var1(icor,iom)/two)**2))
                 end do
               end do
             end if

!            Emissivity for each omega
             if (need_emissivity) then
               do iom=1,mom
                 do icor=1,nphicor
                   oml=-energy_cor(icor)-oml_emis(iom)
                   emisx_k(icor,iom,iatom)=emisx_k(icor,iom,iatom) &
&                     +dhdk2_g(icor)*occ_k(iband)/oml &
&                     *dexp(-((-energy_cor(icor)+eig0_k(iband)-oml)/dom)**2)
                 end do
               end do
             end if

           end do ! iatom

         end if ! my band?
       end do ! iband

!      Accumulate k-point contribution
       if (need_absorption) then
         sigx1(1:nphicor,1:mom,1:natom,isppol)=sigx1(1:nphicor,1:mom,1:natom,isppol) &
&                                              +wtk(ikpt)*sigx1_k(1:nphicor,1:mom,1:natom)   
       end if
       if (need_emissivity) then
         emisx(1:nphicor,1:mom,1:natom,isppol)=emisx(1:nphicor,1:mom,1:natom,isppol) &
&                                              +wtk(ikpt)*emisx_k(1:nphicor,1:mom,1:natom)   
       end if

!      Validity limit
       deltae=deltae+eig0_k(nband_k)

       ABI_FREE(eig0_k)
       ABI_FREE(occ_k)

!    End loop over kpt/spin
     end if ! My kpt?
     bdtot_index=bdtot_index+nband_k
   end do ! ikpt
 end do ! isppol

!Accumulate kpt/band contributions over processors
 if (need_absorption) then
   call xmpi_sum(sigx1,comm,mpierr)
 end if
 if (need_emissivity) then
   call xmpi_sum(emisx,comm,mpierr)
 end if
 call xmpi_sum(deltae,comm,mpierr)
 deltae=deltae/mpi_enreg%nproc_band

!Release some memory
 ABI_FREE(dhdk2_g)
 ABI_FREE(psinablapsi2)
 if (need_absorption) then
   ABI_FREE(sigx1_k)
 end if
 if (need_emissivity) then
   ABI_FREE(emisx_k)
 end if

 !Close core-valence dipoles file
 if (iomode == IO_MODE_ETSF) then
   if (iomode_estf_mpiio.or.me==master) then
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_close(ncid))
#endif
   end if
 else if (me==master) then
   ierr=close_unit(opt2_unt,msg)
   ABI_CHECK(ierr==0,sjoin("Error while closing ",filnam2))
 end if

!---------------------------------------------------------------------------------
! Post-processing

!Absorption: post-processing
 if (need_absorption) then

!  Absorption: apply scaling factor
   do isppol=1,nsppol
     do iatom=1,natom
       do iom=1,mom
         do icor=1,nphicor
           sigx1(icor,iom,iatom,isppol)=sigx1(icor,iom,iatom,isppol)*two_pi*dble(natom)/ucvol/two/(dabs(energy_cor(icor)))
         end do
       end do
     end do
   end do

!  Absorption: average over atoms
   ABI_MALLOC(sigx1_av,(nphicor,mom,nsppol))
   sigx1_av=zero
   do isppol=1,nsppol
     do iatom=1,natom
       do iom=1,mom
         do icor=1,nphicor
           sigx1_av(icor,iom,isppol)=sigx1_av(icor,iom,isppol)+sigx1(icor,iom,iatom,isppol)
         end do
       end do
     end do
   end do

!  Absorption: spin treatment
   if(nsppol==2) then
     ABI_MALLOC(sum_spin_sigx1,(nphicor,mom,natom))
     ABI_MALLOC(sum_spin_sigx1_av,(nphicor,mom))
     sum_spin_sigx1=zero ; sum_spin_sigx1_av=zero
     do isppol=1,nsppol
       do iatom=1,natom
         do iom=1,mom
           do icor=1,nphicor
             sum_spin_sigx1(icor,iom,iatom)=sum_spin_sigx1(icor,iom,iatom) &
&                                          +sigx1(icor,iom,iatom,isppol)
           end do
         end do
       end do
     end do
     do isppol=1,nsppol
       do icor=1,nphicor
         do iom=1,mom
           sum_spin_sigx1_av(icor,iom)=sum_spin_sigx1_av(icor,iom)+sigx1_av(icor,iom,isppol)
         end do
       end do
     end do
   endif
 end if
   
!Emissivity: post-processing
 if (need_emissivity) then

!  Emissivity: filter low values
   do isppol=1,nsppol
     do iatom=1,natom
       do iom=1,mom
         do icor=1,nphicor
           if (emisx(icor,iom,iatom,isppol)<=tol16) emisx(icor,iom,iatom,isppol)=zero
         end do
       end do
     end do
   end do

!  Emissivity: apply scaling factor
   emisx=emisx*two_pi*third*dble(natom)/(dom*ucvol)*half/dsqrt(pi)

!  Emissivity: average over atoms
   ABI_MALLOC(emisx_av,(nphicor,mom,nsppol))
   emisx_av=zero
   do isppol=1,nsppol
     do iatom=1,natom
       do iom=1,mom
         do icor=1,nphicor
           emisx_av(icor,iom,isppol)=emisx_av(icor,iom,isppol)+emisx(icor,iom,iatom,isppol)
         end do
       end do
     end do
   end do
  emisx_av=emisx_av/dble(natom)

!  Emissivity: spin treatment
   if(nsppol==2) then
     ABI_MALLOC(sum_spin_emisx,(nphicor,mom,natom))
     ABI_MALLOC(sum_spin_emisx_av,(nphicor,mom))
     sum_spin_emisx=zero ; sum_spin_emisx_av=zero
     do isppol=1,nsppol
       do iatom=1,natom
         do iom=1,mom
           do icor=1,nphicor
             sum_spin_emisx(icor,iom,iatom)=sum_spin_emisx(icor,iom,iatom) &
&                                          +emisx(icor,iom,iatom,isppol)
           end do
         end do
       end do
     end do
     do isppol=1,nsppol
       do icor=1,nphicor
         do iom=1,mom
           sum_spin_emisx_av(icor,iom)=sum_spin_emisx_av(icor,iom)+emisx_av(icor,iom,isppol)
         end do
       end do
     end do
   endif
 end if
  
!---------------------------------------------------------------------------------
! Output results
   
!Only master node outputs results in files (only master node)
 if (me==master) then

!  Standard output
   if (need_absorption) then
     write(std_out,*) 'Absorption: valence state orbital energies: omin,omax',omin_sig,omax_sig
   end if
   if (need_emissivity) then
     write(std_out,*) 'Emissivity: valence state orbital energies: omin,omax',omin_emis,omax_emis
   end if
   if (need_absorption) then
     write(std_out,'(a,f10.5,a,f10.5,a)')&
&   ' Emax       =',deltae/dble(nkpt*nsppol),' Ha',deltae/dble(nkpt*nsppol)*Ha_eV,' eV'
   end if

!  _sigX file
   if (need_absorption) then
     if (open_file(trim(filnam_out)//'_sigX',msg,newunit=sigx1_unt,form='formatted',action="write")/=0) then
       ABI_ERROR(msg)
     end if
     if (abs(dom_max)>tol10) then
       write(sigx1_unt,'(a)')'#***************************************************** ARCTAN SMEARING ********'
       write(sigx1_unt,'(a,i8,3f10.5,a)')'# npts,omin,omax,width      =',mom,omin,omax,dom,' Ha'
       write(sigx1_unt,'(a,2f10.5,a)')'# dom_max,center        =',dom_max,dom_ctr,' Ha'
     else
       write(sigx1_unt,'(a)')'#***************************************************** FIXED SMEARING ********'
       write(sigx1_unt,'(a,i8,3f10.5,a)')'# npts,omin,omax,width      =',mom,omin,omax,dom,' Ha'
     endif
     write(sigx1_unt,'(a,3f10.5,a)' )'# rprimd(bohr)      =',rprimd(1:3,1)
     write(sigx1_unt,'(a,3f10.5,a)' )'#                    ',rprimd(1:3,2)
     write(sigx1_unt,'(a,3f10.5,a)' )'#                    ',rprimd(1:3,3)
     write(sigx1_unt,'(a,i8)' )      '# natom             =',natom
     write(sigx1_unt,'(a,3i8)' )     '# nkpt,mband,nsppol        =',nkpt,mband,nsppol
     write(sigx1_unt,'(a, f10.5,a)' )'# ecut                  =',ecut,' Ha'
     write(sigx1_unt,'(a,f10.5,a,f10.5,a)' )'# fermie            =',fermie,' Ha ',fermie*Ha_eV,' eV'
     write(sigx1_unt,'(a,f12.5,a,f12.5,a)') '# Temp              =',tsmear,' Ha ',Tatm,' Kelvin'
     write(sigx1_unt,'(a)')'#----------------------------------------------------------------------------'
     write(sigx1_unt,'(a,i4)') '# Number of core orbitals nc=',nphicor
     do icor=1,nphicor
       if (kappacor(icor)==0) then
         write(sigx1_unt,'(a,2i4,2f10.5)') '# n, l, Energy(Ha), Edge(Ha): ',ncor(icor),lcor(icor),energy_cor(icor),edge(icor)
       else
         if (kappacor(icor)>0) then
           j2=2*lcor(icor)-1
           write(sigx1_unt,'(a,i4,i4,a,i4,2f10.5)') '# n, j, l, Energy(Ha), Edge(Ha): ', &
&            ncor(icor),j2,' / 2',lcor(icor),energy_cor(icor),edge(icor)
         else
           if (kappacor(icor)<-1) then
             j2=2*lcor(icor)+1
             write(sigx1_unt,'(a,i4,i4,a,i4,2f10.5)') '# n, j, l, Energy(Ha), Edge(Ha): ', &
&              ncor(icor),j2,' / 2',lcor(icor),energy_cor(icor),edge(icor)
           else
             write(sigx1_unt,'(a,i4,a,i4,2f10.5)') '# n, j, l, Energy(Ha), Edge(Ha): ', &
&              ncor(icor),'   1 / 2',lcor(icor),energy_cor(icor),edge(icor)
           end if
         end if
       end if
     end do
     write(sigx1_unt,'(a)')'#----------------------------------------------------------------------------'
     write(sigx1_unt,'(a,f10.5,a,f10.5,a)')&
&     '# Emax       =',deltae/dble(nkpt*nsppol),' Ha',deltae/dble(nkpt*nsppol)*Ha_eV,' eV'
     write(sigx1_unt,'(a)')'#----------------------------------------------------------------------------'
     do iom=1,mom
       write(sigx1_unt,'(100(1x,e15.8))') &
&      (oml_edge(icor,iom),sigx1_av(icor,iom,1)/dble(natom),sigx1(icor,iom,atnbr,1),icor=1,nphicor)
     end do
     close(sigx1_unt)
   end if

!    _s_sigX and tot_sigX files
   if (need_absorption.and.nsppol==2) then
     if (open_file(trim(filnam_out)//'_s_sigX',msg,newunit=sigx1_s_unt,form='formatted',action="write")/=0) then
       ABI_ERROR(msg)
     end if
     if (open_file(trim(filnam_out)//'_tot_sigX',msg,newunit=sigx1_tot_unt,form='formatted',action="write")/=0) then
       ABI_ERROR(msg)
     end if
     do iom=1,mom
       write(sigx1_s_unt,'(100(1x,e15.8))') &
&        (oml_edge(icor,iom),sigx1_av(icor,iom,nsppol)/dble(natom),sigx1(icor,iom,atnbr,nsppol),icor=1,nphicor)
       write(sigx1_tot_unt,'(100(1x,e15.8))') &
&          (oml_edge(icor,iom),sum_spin_sigx1_av(icor,iom)/dble(natom),sum_spin_sigx1(icor,iom,atnbr),icor=1,nphicor)
     end do
     close(sigx1_s_unt)
     close(sigx1_tot_unt)
   end if

!  _emisX file
   if (need_emissivity) then
     if (open_file(trim(filnam_out)//'_emisX',msg,newunit=ems_unt,form='formatted',action="write")/=0) then
       ABI_ERROR(msg)
     end if
     write(ems_unt,*) '# conducti: Xray emission spectrum, all in atomic units by default '
     write(ems_unt,*) '# One block of 3 columns per core wavefunction'
     write(ems_unt,*) '# energy, sigx_av, sigx, etc... '
     do iom=1,mom
       write(ems_unt,'(3(3(1x,e15.8),2x))') &
&       ((-energy_cor(icor)+oml_emis(iom)),emisx_av(icor,iom,1),emisx(icor,iom,atnbr,1),icor=1,nphicor)
     end do
    close(ems_unt)
  end if
    
!    _s_emisX and tot_emisX files
   if (need_emissivity.and.nsppol==2) then
     if (open_file(trim(filnam_out)//'_s_emisX',msg,newunit=ems_s_unt,form='formatted',action="write")/=0) then
       ABI_ERROR(msg)
     end if
     if (open_file(trim(filnam_out)//'_tot_emisX',msg,newunit=ems_tot_unt,form='formatted',action="write")/=0) then
       ABI_ERROR(msg)
     end if
     do iom=1,mom
       write(ems_s_unt,'(3(3(1x,e15.8),2x))') &
&        ((-energy_cor(icor)+oml_emis(iom)),emisx_av(icor,iom,nsppol)/dble(natom),emisx(icor,iom,atnbr,nsppol),icor=1,nphicor)
       write(ems_tot_unt,'(3(3(1x,e15.8),2x))') &
&       ((-energy_cor(icor)+oml_emis(iom)),sum_spin_emisx_av(icor,iom)/dble(natom),sum_spin_emisx(icor,iom,atnbr),icor=1,nphicor)
     end do
     close(ems_s_unt)
     close(ems_tot_unt)
   end if
   
 endif ! master node

!---------------------------------------------------------------------------------
! End

!Release memory space
 if (need_absorption) then
   ABI_FREE(sigx1)
   ABI_FREE(sigx1_av)
   if (nsppol==2) then
     ABI_FREE(sum_spin_sigx1)
     ABI_FREE(sum_spin_sigx1_av)
   end if
   ABI_FREE(dom_var1)
   ABI_FREE(oml_edge)
 end if
 if (need_emissivity) then
   ABI_FREE(emisx)
   ABI_FREE(emisx_av)
   if (nsppol==2) then
     ABI_FREE(sum_spin_emisx)
     ABI_FREE(sum_spin_emisx_av)
   end if
   ABI_FREE(oml_emis)
 end if
 ABI_FREE(ncor)
 ABI_FREE(lcor)
 ABI_FREE(kappacor)
 ABI_FREE(energy_cor)
 ABI_FREE(edge)
 ABI_FREE(eigen0)
 ABI_FREE(nband)
 ABI_FREE(occ)
 ABI_FREE(wtk)
 call hdr%free()
 call destroy_mpi_enreg(mpi_enreg)

end subroutine conducti_paw_core
!!***

!----------------------------------------------------------------------

!!****f* m_conducti/conducti_nc
!! NAME
!! conducti_nc
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula.
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
!!  dom=frequency range
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
!!  kin11= Onsager kinetic coeficient=optical conductivity
!!  kin12= Onsager kinetic coeficient
!!  kin21= Onsager kinetic coeficient
!!  kin22= Onsager kinetic coeficient
!!  Kth=thermal conductivity
!!  mom=number of frequency for conductivity computation
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
!!  Sth=thermopower
!!  tsmear=smearing width (or temperature) in Hartree
!!  ucvol=unit cell volume in ($\textrm{bohr}^{3}$).
!!  wind=frequency windows for computations of sigma
!!  wtk(nkpt)=weight assigned to each k point.
!!  znucl(natom)=atomic number of atoms
!!  np_sum=noziere-pines sumrule
!!  cond_kg(mom)=kubo-greenwood conductivity
!!
!! PARENTS
!!      conducti
!!
!! CHILDREN
!!      ddk%close,ddk%compare,ddk%read_eigk,getnel,gswfk%read_eigk,gswfk%close
!!      jacobi,msig,nctk_fort_or_ncfile,openfile,wfk_open_read
!!
!! SOURCE


subroutine conducti_nc(filnam,filnam_out)

!Arguments -----------------------------------
!scalars
 character(len=fnlen) :: filnam,filnam_out

!Local variables-------------------------------
!scalars
 integer,parameter :: formeig0=0,formeig1=1
 integer :: bantot,bd2tot_index,bdtot0_index,bdtot_index
 integer :: headform,iband,ii,jj,ikpt,iunt
 integer :: index_1,iom,isppol,jband,l1,l2,mband,mom,natom,nband1
 integer :: nrot,iomode
 integer :: nband_k,nkpt,nlign,nrest,nspinor,nsppol,ntypat
 integer :: occopt,comm
 integer :: tens_unt,lij_unt,sig_unt,kth_unt,ocond_unt
 real(dp) :: deltae,dosdeltae,diff_occ,dom,ecut,entropy,fermie,maxocc
 real(dp) :: nelect,np_sum,np_sum_k1,np_sum_k2,omin,oml,socc,socc_k,sig
 real(dp) :: tphysel,tsmear,ucvol,wind,Tatm
 character(len=fnlen) :: filnam0,filnam1,filnam2,filnam3
 character(len=500) :: msg
 type(hdr_type) :: hdr
 type(wfk_t) :: gswfk,ddk1,ddk2,ddk3
!arrays
 integer,allocatable :: nband(:)
 real(dp) :: gmet(3,3),gmet_inv(3,3),gprimd(3,3),gprimd_inv(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: cond_kg(:,:,:),cond_kg_cart(:,:,:),cond_nd(:,:,:),dhdk2_r(:,:,:,:),dhdk2_g(:,:)
 real(dp),allocatable :: doccde(:),doccde_k(:),cond_kg_xx(:),cond_kg_yy(:),cond_kg_zz(:),trace(:)
 real(dp),allocatable :: eig0_k(:),eig0tmp(:),eig1_k(:,:),eigen0(:),eigen11(:)
 real(dp),allocatable :: eigen12(:),eigtmp(:)
 real(dp),allocatable :: eigen13(:),occ(:),occ_k(:),wtk(:),cond_tot(:),oml1(:)
 real(dp),allocatable :: kin11(:),kin12(:),kin21(:),kin22(:)
 real(dp),allocatable :: kin11_k(:),kin12_k(:),kin21_k(:),kin22_k(:),Kth(:),Stp(:)
 real(dp) :: cond_kg_w(3,3),z(3,3)
 real(dp) :: eig_cond(3)

! *********************************************************************************

!Read data file
 if (open_file(filnam,msg,newunit=iunt,form='formatted',status="old")/=0) then
   ABI_ERROR(msg)
 end if

 rewind(iunt)
 read(iunt,*)
 read(iunt,'(a)')filnam1       ! first ddk file
 read(iunt,'(a)')filnam2       ! second ddk file
 read(iunt,'(a)')filnam3       ! third ddk file
 read(iunt,'(a)')filnam0       ! ground-state data

!Open the GS Wavefunction file and the 3 DDK files.
! TODO: one should perform basic consistency tests for the GS WFK and the DDK files, e.g.
! k-points and their order, spins, number of bands could differ in the four files.
! Note indeed that we are assuming the same numer of bands in all the files.
 comm = xmpi_comm_self
 call nctk_fort_or_ncfile(filnam0, iomode, msg)
 if (len_trim(msg) /= 0) ABI_ERROR(msg)
 call wfk_open_read(gswfk,filnam0, formeig0, iomode, get_unit(), comm)

 call nctk_fort_or_ncfile(filnam1, iomode, msg)
 if (len_trim(msg) /= 0) ABI_ERROR(msg)
 call wfk_open_read(ddk1,filnam1, formeig1, iomode, get_unit(), comm, hdr_out=hdr)

 call nctk_fort_or_ncfile(filnam2, iomode, msg)
 if (len_trim(msg) /= 0) ABI_ERROR(msg)
 call wfk_open_read(ddk2,filnam2, formeig1, iomode, get_unit(), comm)

 call nctk_fort_or_ncfile(filnam3, iomode, msg)
 if (len_trim(msg) /= 0) ABI_ERROR(msg)
 call wfk_open_read(ddk3,filnam3, formeig1, iomode, get_unit(), comm)

 if (ddk1%compare(ddk2) /= 0) then
   ABI_ERROR("ddk1 and ddk2 are not consistent. see above messages")
 end if
 if (ddk1%compare(ddk3) /= 0) then
   ABI_ERROR("ddk1 and ddk3 are not consistent. see above messages")
 end if

!Extract params from the header of the first ddk file (might have been the GS file ?)

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
 ABI_MALLOC(nband,(nkpt*nsppol))
 ABI_MALLOC(occ,(bantot))
 fermie=hdr%fermie
 occ(1:bantot)=hdr%occ(1:bantot)
 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(nband(:))

 write(std_out,*)
 write(std_out,'(a,3f10.5,a)' )' rprimd(bohr)      =',rprimd(1:3,1)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,2)
 write(std_out,'(a,3f10.5,a)' )'                    ',rprimd(1:3,3)
 write(std_out,'(a,i8)')       ' natom             =',natom
 write(std_out,'(a,2i8)')      ' nkpt,mband        =',nkpt,mband
 write(std_out,'(a, f10.5,a)' ) ' ecut              =',ecut,' Ha'
 write(std_out,'(a,f10.5,a,f10.5,a)' )' fermie            =',fermie,' Ha',fermie*Ha_eV,' eV'

!Prepare the reading of ddk Wff files
 ABI_MALLOC(eigtmp,(2*mband*mband))
 ABI_MALLOC(eig0tmp,(mband))

!Read the eigenvalues of ground-state and ddk files
 ABI_MALLOC(eigen0,(mband*nkpt*nsppol))
 ABI_MALLOC(eigen11,(2*mband*mband*nkpt*nsppol))
 ABI_MALLOC(eigen12,(2*mband*mband*nkpt*nsppol))
 ABI_MALLOC(eigen13,(2*mband*mband*nkpt*nsppol))
 bdtot0_index=0 ; bdtot_index=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband1=nband(ikpt+(isppol-1)*nkpt)
     call gswfk%read_eigk(ikpt,isppol,xmpio_single,eig0tmp)
     eigen0(1+bdtot0_index:nband1+bdtot0_index)=eig0tmp(1:nband1)

     call ddk1%read_eigk(ikpt,isppol,xmpio_single,eigtmp)
     eigen11(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)

     call ddk2%read_eigk(ikpt,isppol,xmpio_single,eigtmp)
     eigen12(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)

     call ddk3%read_eigk(ikpt,isppol,xmpio_single,eigtmp)
     eigen13(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)

     bdtot0_index=bdtot0_index+nband1
     bdtot_index=bdtot_index+2*nband1**2
   end do
 end do

!Close files
 call gswfk%close()
 call ddk1%close()
 call ddk2%close()
 call ddk3%close()

 ABI_FREE(eigtmp)
 ABI_FREE(eig0tmp)

!---------------------------------------------------------------------------------
!Gmet inversion
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 call matr3inv(gmet,gmet_inv)
 call matr3inv(gprimd,gprimd_inv)

!---------------------------------------------------------------------------------
!Derivative of occupation wrt the energy.

 ABI_MALLOC(doccde,(mband*nkpt*nsppol))
 ABI_MALLOC(wtk,(nkpt))

 read(iunt,*)tsmear
 Tatm=tsmear*Ha_K
 write(std_out,'(a,f12.5,a,f12.5,a)') ' Temp              =',tsmear,' Ha ',Tatm,' Kelvin'
!
 nlign=nkpt/6
 nrest=nkpt-6*nlign
 index_1=0
 do ii=1,nlign
   read(iunt,*)wtk(1+index_1:6+index_1)
   index_1=index_1+6
 end do
 if (nrest/=0) then
   read(iunt,*)wtk(6*nlign+1:nkpt)
 end if
!
 if (occopt==1) then
   write(std_out,'(a,i4)')  ' occopt            =',occopt
   doccde=zero
 else
   tphysel=zero
   maxocc=two/(nsppol*nspinor)
   dosdeltae=zero
!  CP: using 1 and nband(0) as dummy value, because function
!  not implemented for occopt==9; adding fermih=fermie in the list of arguments as well
   call getnel(doccde,dosdeltae,eigen0,entropy,fermie,fermie,maxocc,mband,nband,&
&   nelect,nkpt,nsppol,occ,occopt,1,tphysel,tsmear,11,wtk,1,nband(1))
 end if

!---------------------------------------------------------------------------------
!Size of the frequency range

 read(iunt,*)dom,wind
 close(iunt)
 mom=int(wind/dom)
 ABI_MALLOC(oml1,(mom))
 do iom=1,mom
   oml1(iom)=tol10*1000._dp+dble(iom)*dom
 end do

 ABI_MALLOC(cond_nd,(mom,3,3))
 ABI_MALLOC(cond_kg,(mom,3,3))
 ABI_MALLOC(cond_kg_cart,(mom,3,3))
 ABI_MALLOC(cond_kg_xx,(mom))
 ABI_MALLOC(cond_kg_yy,(mom))
 ABI_MALLOC(trace,(mom))
 ABI_MALLOC(cond_kg_zz,(mom))
 ABI_MALLOC(cond_tot,(mom))
 ABI_MALLOC(kin11,(mom))
 ABI_MALLOC(kin12,(mom))
 ABI_MALLOC(kin21,(mom))
 ABI_MALLOC(kin22,(mom))
 ABI_MALLOC(kin11_k,(mom))
 ABI_MALLOC(kin12_k,(mom))
 ABI_MALLOC(kin21_k,(mom))
 ABI_MALLOC(kin22_k,(mom))
 ABI_MALLOC(Kth,(mom))
 ABI_MALLOC(Stp,(mom))
 write(std_out,'(a,i8,2f10.5,a)')' mom,wind,dom      =',mom,wind,dom,' Ha'

!---------------------------------------------------------------------------------

 kin11   = zero
 kin12   = zero
 kin21   = zero
 kin22   = zero
 np_sum  = zero
 socc    = zero
 cond_kg = zero

!LOOP OVER SPINS
 do isppol=1,nsppol
!
   bdtot_index = 0
   bd2tot_index = 0

   deltae  = zero
!
!  BIG FAT k POINT LOOP
!
   do ikpt=1,nkpt

     nband_k=nband(ikpt+(isppol-1)*nkpt)

     ABI_MALLOC(eig0_k,(nband_k))
     ABI_MALLOC(eig1_k,(2*nband_k**2,3))
     ABI_MALLOC(occ_k,(nband_k))
     ABI_MALLOC(doccde_k,(nband_k))
     ABI_MALLOC(dhdk2_r,(3,3,nband_k,nband_k))
     ABI_MALLOC(dhdk2_g,(nband_k,nband_k))

     cond_nd   = zero
     kin11_k   = zero
     kin12_k   = zero
     kin21_k   = zero
     kin22_k   = zero
     np_sum_k1 = zero
     np_sum_k2 = zero
     socc_k    = zero
     dhdk2_r   = zero
     dhdk2_g   = zero

!    eigenvalue for k-point
     eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
!    first derivative eigenvalues for k-point
     eig1_k(:,1)=eigen11(1+bd2tot_index:2*nband_k**2+bd2tot_index)
     eig1_k(:,2)=eigen12(1+bd2tot_index:2*nband_k**2+bd2tot_index)
     eig1_k(:,3)=eigen13(1+bd2tot_index:2*nband_k**2+bd2tot_index)
!    occupation numbers for k-point
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
!    derivative of occupation number for k-point
     doccde_k(:)=doccde(1+bdtot_index:nband_k+bdtot_index)

!    LOOP OVER BAND
     do iband=1,nband_k
       do jband=1,nband_k
!
!        TODO : replace with BLAS calls
         do l1=1,3
           do l2=1,3
             do ii=1,3
               do jj=1,3
                 dhdk2_r(l1,l2,iband,jband)=dhdk2_r(l1,l2,iband,jband)+(rprimd(l1,ii)&
&                 *eig1_k(2*iband-1+(jband-1)*2*nband_k,ii)*&
&                 rprimd(l2,jj)*eig1_k(2*iband-1+(jband-1)*2*nband_k,jj)&
&                 +rprimd(l1,ii)*eig1_k(2*iband  +(jband-1)*2*nband_k,ii)*&
&                 rprimd(l2,jj)*eig1_k(2*iband+(jband-1)*2*nband_k,jj))
               end do
             end do
           end do
         end do

         do l1=1,3
           do l2=1,3
             dhdk2_r(l1,l2,iband,jband)=dhdk2_r(l1,l2,iband,jband)/two_pi/two_pi
           end do
         end do

!        TODO: replace with BLAS calls
         do l1=1,3
           do l2=1,3
             dhdk2_g(iband,jband)=dhdk2_g(iband,jband)+gmet_inv(l1,l2)*( &
&             eig1_k(2*iband-1+(jband-1)*2*nband_k,l1)*&
&             eig1_k(2*iband-1+(jband-1)*2*nband_k,l2) &
&            +eig1_k(2*iband  +(jband-1)*2*nband_k,l1)*&
&             eig1_k(2*iband  +(jband-1)*2*nband_k,l2))
           end do
         end do
         dhdk2_g(iband,jband)=dhdk2_g(iband,jband)/two_pi/two_pi

         diff_occ = occ_k(iband)-occ_k(jband)
!        if (dabs(diff_occ)>=tol8) then

!        Conductivity for each omega
         omin = zero
         do iom=1,mom
           oml=oml1(iom)
           if (jband>iband) then
             sig= dhdk2_g(iband,jband)&
&             *(diff_occ)/oml*(dexp(-((eig0_k(jband)-eig0_k(iband)-oml)/dom)**2)&
&             -dexp(-((eig0_k(iband)-eig0_k(jband)-oml)/dom)**2))
             kin11_k(iom)=kin11_k(iom)+sig
             kin12_k(iom)=kin12_k(iom)-sig*(eig0_k(jband)-fermie)
             kin21_k(iom)=kin21_k(iom)-sig*(eig0_k(iband)-fermie)
             kin22_k(iom)=kin22_k(iom) + &
&             sig*(eig0_k(iband)-fermie)*(eig0_k(jband)-fermie)
           end if
           do l1=1,3
             do l2=1,3
               cond_nd(iom,l1,l2)=cond_nd(iom,l1,l2) +dhdk2_r(l1,l2,iband,jband)&
&               *(diff_occ)/oml*dexp(-((eig0_k(jband)-eig0_k(iband)-oml)/dom)**2)
             end do
           end do

         end do

!        Sumrule start
         if (dabs(eig0_k(iband)-eig0_k(jband))>=tol10) then
           np_sum_k1=np_sum_k1 -dhdk2_g(iband,jband)&
&           *(diff_occ)/(eig0_k(iband)-eig0_k(jband))
         else
           np_sum_k2=np_sum_k2 - doccde_k(iband)*dhdk2_g(iband,jband)
         end if


!        end loop over band
       end do
       socc_k=socc_k+occ_k(iband)
     end do

     do iom=1,mom
       kin11(iom)=kin11(iom)+wtk(ikpt)*kin11_k(iom)
       kin12(iom)=kin12(iom)+wtk(ikpt)*kin12_k(iom)
       kin21(iom)=kin21(iom)+wtk(ikpt)*kin21_k(iom)
       kin22(iom)=kin22(iom)+wtk(ikpt)*kin22_k(iom)
       do l1=1,3
         do l2=1,3
           cond_kg(iom,l1,l2)=cond_kg(iom,l1,l2)+wtk(ikpt)*cond_nd(iom,l1,l2)
         end do
       end do
     end do

     np_sum=np_sum + wtk(ikpt)*(np_sum_k1+np_sum_k2)
     socc=socc+wtk(ikpt)*socc_k

!    Validity limit
     deltae=deltae+(eig0_k(nband_k)-fermie)

     bd2tot_index=bd2tot_index+2*nband_k**2
     bdtot_index=bdtot_index+nband_k
     ABI_FREE(eig0_k)
     ABI_FREE(eig1_k)
     ABI_FREE(occ_k)
     ABI_FREE(doccde_k)
     ABI_FREE(dhdk2_r)
     ABI_FREE(dhdk2_g)
!    End loop over k
   end do

   write(std_out,'(a,3f10.5)')' sumrule           =',np_sum/socc/three,socc
   write(std_out,'(a,f10.5,a,f10.5,a)')&
&   ' Emax-Efermi       =',deltae/dble(nkpt),' Ha',deltae/dble(nkpt)*Ha_eV,' eV'

!End loop over spins
 end do

 cond_kg=cond_kg*two_pi*third/(dom*ucvol)*half/dsqrt(pi)


!Check that new output file does NOT exist
!Keep this line : prevent silly (compiler ?) bug on HP 8000
 write(std_out,*)' conducti : call isfile '
!
 if (open_file(trim(filnam_out)//'_tens',msg,newunit=tens_unt,form='formatted',action="write")/=0) then
   ABI_ERROR(msg)
 end if
 if (open_file(trim(filnam_out)//'_Lij',msg,newunit=lij_unt,form='formatted',action="write")/=0) then
   ABI_ERROR(msg)
 end if
 write(lij_unt,'(a)')' # omega(ua) L12 L21 L22 L22'

 if (open_file(trim(filnam_out)//'_sig',msg,newunit=sig_unt,form='formatted',action="write")/=0) then
   ABI_ERROR(msg)
 end if
 write(sig_unt,'(a)')' # omega(ua) hbar*omega(eV)    cond(ua)             cond(ohm.cm)-1'

 if (open_file(trim(filnam_out)//'_Kth',msg,newunit=kth_unt,form='formatted',action="write")/=0) then
   ABI_ERROR(msg)
 end if
 write(kth_unt,'(a)')&
& ' #omega(ua) hbar*omega(eV)  thermal cond(ua)   Kth(W/m/K)   thermopower(ua)   Stp(microohm/K)'

 if (open_file(trim(filnam_out)//'.out',msg,newunit=ocond_unt,form='formatted',action="write")/=0) then
   ABI_ERROR(msg)
 end if
 write(ocond_unt,'(a)' )' Conducti output file:'
 write(ocond_unt,'(a)' )' Contains all results produced by conducti utility'
 write(ocond_unt,'(a)' )' '
 write(ocond_unt,'(a)')' # omega(ua)       cond(ua)             thermal cond(ua)       thermopower(ua)'

!Keep this line : prevent silly (compiler ?) bug on HP 8000
 write(std_out,*)' conducti : after call isfile '

!Compute thermal conductivity and thermopower
 do iom=1,mom
   oml=oml1(iom)
   kin11(iom)=kin11(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   kin21(iom)=kin21(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   kin12(iom)=kin12(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   kin22(iom)=kin22(iom)*two_pi*third/(dom*ucvol)*half/dsqrt(pi)
   if (dabs(kin11(iom))<10.0d-20) kin11(iom)=zero
   Kth(iom)=kin22(iom)
   Stp(iom)=zero
   if(kin11(iom)/=zero)  then
     Kth(iom)=Kth(iom)-(kin12(iom)*kin21(iom)/kin11(iom))
     Stp(iom)=kin12(iom)/(kin11(iom)*Tatm)
   end if
   if (dabs(Kth(iom))<10.0d-20) Kth(iom)=zero
   if (dabs(Stp(iom))<10.0d-20) Stp(iom)=zero
   if (dabs(kin12(iom))<10.0d-20) kin12(iom)=zero
   if (dabs(kin21(iom))<10.0d-20) kin21(iom)=zero
   if (dabs(kin22(iom))<10.0d-20) kin22(iom)=zero

   write(lij_unt,'(f12.5,4es22.12)')oml,kin12(iom),kin21(iom),kin22(iom),kin22(iom)/Tatm*3.4057d9
   write(sig_unt,'(2f12.5,2es22.12)') oml,oml*Ha_eV,kin11(iom),kin11(iom)*Ohmcm
   write(kth_unt,'(2f12.5,4es22.12)') oml,oml*Ha_eV,Kth(iom),Kth(iom)*3.4057d9/Tatm,Stp(iom),Stp(iom)*3.6753d-2
   write(ocond_unt,'(1f12.5,3es22.12)') oml,kin11(iom),Kth(iom),Stp(iom)
 end do

 write(tens_unt,'(a)' )' Conductivity file '
 write(tens_unt,'(a)' )' ----------------- '
 write(tens_unt,'(a)' )' Contain first the full conductivity tensor, for the desired set of energies,'
 write(tens_unt,'(a)' )' then, the three principal values, for the desired set of energies'
 write(tens_unt,'(a)' )' (note that eigenvalues are not directly associated with xx,yy,zz)'
 write(tens_unt,'(a)' )' '

 write(ocond_unt,'(a)' )' '
 write(ocond_unt,'(a)' )' full conductivity tensor, for the desired set of energies'
 write(ocond_unt,'(a)' )' then, the three principal values, for the desired set of energies:'

 do iom=1,mom
   oml=oml1(iom)*Ha_eV
   write(tens_unt, '(a,es16.6,a)' ) ' energy (in eV) =',oml,', conductivity tensor (in Ohm.cm-1) follows :'
   write(ocond_unt, '(a,es16.6,a)' ) ' energy (in eV) =',oml,', conductivity tensor (in Ohm.cm-1) follows :'
   do l1=1,3
     write(tens_unt,"(3f25.15)") (cond_kg(iom,l1,l2)*Ohmcm,l2=1,3)
     write(ocond_unt,"(3f25.15)") (cond_kg(iom,l1,l2)*Ohmcm,l2=1,3)
   end do
 end do

!Diagonalizing the conductivity matrix for sigma_xx,sigma_yy,sigma_zz
 cond_kg_xx=0d0
 cond_kg_yy=0d0
 cond_kg_zz=0d0
!trace=0d0 ! Used for checking with the original version of the code
 do iom=1,mom
   oml=oml1(iom)*Ha_eV
   cond_kg_w=0d0
   do l1=1,3
     do l2=1,3
       cond_kg_w(l1,l2)=cond_kg(iom,l1,l2)
     end do
   end do
   call jacobi(cond_kg_w,3,3,eig_cond,z,nrot)

!  When the value is too small, set it to zero before printing
   if(abs(eig_cond(1))<tol10)eig_cond(1)=zero
   if(abs(eig_cond(2))<tol10)eig_cond(2)=zero
   if(abs(eig_cond(3))<tol10)eig_cond(3)=zero

   cond_kg_xx(iom)=eig_cond(1)
   cond_kg_yy(iom)=eig_cond(2)
   cond_kg_zz(iom)=eig_cond(3)
!  trace(iom)=cond_kg_xx(iom)+cond_kg_yy(iom)+cond_kg_zz(iom)
 end do

!DEBUG Keep this line : prevent silly (compiler ?) bug on HP 8000
!write(std_out,*)' conducti : after open '
!ENDDEBUG

 write(tens_unt,'(a,a)')ch10,' Now, print principal values of the conductivity tensor.'
 write(tens_unt,'(a)')' '
 write(tens_unt,'(a)')' #omega(ua)   cond_1(ua)     cond_2(ua) cond_3(ua)  cond_tot(ua)'

 write(ocond_unt,'(a)')' '
 write(ocond_unt,'(a,a)')ch10,' Now, print principal values of the conductivity tensor.'
 write(ocond_unt,'(a)')' '
 write(ocond_unt,'(a)')' #omega(ua)   cond_1(ua)     cond_2(ua) cond_3(ua)  cond_tot(ua)'


 do iom=1,mom
   cond_tot(iom)=cond_kg_xx(iom)+cond_kg_yy(iom)+cond_kg_zz(iom)
   write(tens_unt,'(f12.5,4es22.12)')oml1(iom),cond_kg_xx(iom),cond_kg_yy(iom),cond_kg_zz(iom),cond_tot(iom)
   write(ocond_unt,'(f12.5,4es22.12)')oml1(iom),cond_kg_xx(iom),cond_kg_yy(iom),cond_kg_zz(iom),cond_tot(iom)
 end do

 write(tens_unt,*)
 write(tens_unt,'(a)')' #hbar*omega(eV)    cond_1(ohm.cm)-1    cond_2(ohm.cm)-1    cond_3(ohm.cm)-1    cond_t(ohm.cm)-1'
 write(ocond_unt,*)
 write(ocond_unt,'(a)')' #hbar*omega(eV)    cond_1(ohm.cm)-1    cond_2(ohm.cm)-1    cond_3(ohm.cm)-1    cond_t(ohm.cm)-1'

 do iom=1,mom
   oml=oml1(iom)*Ha_eV
   cond_tot(iom)=cond_tot(iom)*Ohmcm
   cond_kg_xx(iom)=cond_kg_xx(iom)*Ohmcm
   cond_kg_yy(iom)=cond_kg_yy(iom)*Ohmcm
   cond_kg_zz(iom)=cond_kg_zz(iom)*Ohmcm
   write(tens_unt,'(f12.5,4es22.12)')oml,cond_kg_xx(iom),cond_kg_yy(iom),cond_kg_zz(iom),cond_tot(iom)
   write(ocond_unt,'(f12.5,4es22.12)')oml,cond_kg_xx(iom),cond_kg_yy(iom),cond_kg_zz(iom),cond_tot(iom)
 end do

!Calculate the imaginary part of the conductivity (principal value)
!+derived optical properties.
 call msig(kin11,mom,oml1,filnam_out)

 close(tens_unt)
 close(lij_unt)
 close(sig_unt)
 close(kth_unt)
 close(ocond_unt)

 ABI_FREE(nband)
 ABI_FREE(oml1)
 ABI_FREE(occ)
 ABI_FREE(eigen11)
 ABI_FREE(eigen12)
 ABI_FREE(eigen13)
 ABI_FREE(eigen0)
 ABI_FREE(doccde)
 ABI_FREE(wtk)
 ABI_FREE(cond_nd)
 ABI_FREE(cond_kg)
 ABI_FREE(cond_kg_cart)
 ABI_FREE(cond_kg_xx)
 ABI_FREE(cond_kg_yy)
 ABI_FREE(trace)
 ABI_FREE(cond_kg_zz)
 ABI_FREE(cond_tot)
 ABI_FREE(kin11)
 ABI_FREE(kin22)
 ABI_FREE(kin12)
 ABI_FREE(kin21)
 ABI_FREE(kin11_k)
 ABI_FREE(kin22_k)
 ABI_FREE(kin12_k)
 ABI_FREE(kin21_k)
 ABI_FREE(Stp)
 ABI_FREE(Kth)
 call hdr%free()

 end subroutine conducti_nc
!!***

!----------------------------------------------------------------------

!!****f* m_conducti/msig
!! NAME
!! msig
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula for PAW formalism
!!
!! INPUTS
!!  fcti(npti)=  conductivity, as calculated in conducti
!!  npti= number of points to calculate conductivity
!!  xi(npti)= energies where the conductivity is calculated
!!
!! OUTPUT
!!   no output, only files
!!
!! NOTES
!!     this program calculates the imaginary part of the conductivity (principal value)
!!     +derived optical properties.
!!     the calculation is performed on the same grid as the initial input
!!     to calculate the principal value, a trapezoidale integration +taylor expansion to
!!     third order is used (W.J. Thomson computer in physics vol 12 p94 1998)
!!    two input files are needed inppv.dat (parameters) and sigma.dat (energy,sigma_1)
!!     two output files ppsigma.dat (energy,sigma_1,sigma_2,epsilon_1,epsilon_2)
!!                      abs.dat     (energy,nomega,komega,romega,absomega)
!!     march 2002 s.mazevet
!!
!! PARENTS
!!      m_conducti
!!
!! CHILDREN
!!      intrpl,open_file
!!
!! SOURCE

subroutine msig(fcti,npti,xi,filnam_out_sig)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: npti
!arrays
 real(dp),intent(in) :: fcti(npti),xi(npti)
 character(len=fnlen),intent(in) :: filnam_out_sig

!Local variables-------------------------------
!scalars
 integer :: npt = 10000
 integer :: ii,ip,npt1,npt2,eps_unt,abs_unt
 real(dp),parameter :: del=0.001_dp,ohmtosec=9.d11
 real(dp) :: dx,dx1,dx2,eps1,eps2,idel,komega,pole,refl,sigma2,xsum
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: abso(:),fct(:),fct1(:),fct2(:),fct3(:),fct4(:),fct5(:)
 real(dp),allocatable :: fctii(:),fp(:),fpp(:),fppp(:),nomega(:),ppsig(:)
 real(dp),allocatable :: x1(:),x2(:)

! *********************************************************************************

 if (npti > 12000) then
   msg = "Sorry - the interpolator INTRPL is hard coded for maximum 12000 points." // &
&        ch10 // " Reduce the conducti input npti, or implement a better interpolator!"
   ABI_ERROR(msg)
 end if

 write(std_out,'(2a)')ch10,'Calculate the principal value and related optical properties'
 write(std_out,'(a)')'following W.J. Thomson computer in physics vol 12 p94 1998 for '
 write(std_out,'(a)')'the principal value.'
 write(std_out,'(a)')'OPTIONS'
 write(std_out,'(a)')'use default number of integration pts: npt=10000'
 write(std_out,'(a)')'Use default value for delta interval: del=1e-3'

 if (open_file(trim(filnam_out_sig)//'_eps',msg,newunit=eps_unt,status='replace',action="write")/=0) then
   ABI_ERROR(msg)
 end if
 write(eps_unt,'(a)')'#energy (eV),sigma_1(Ohm-1cm-1),sigma_2(Ohm-1cm-1),epsilon_1,epsilon_2'

 if (open_file(trim(filnam_out_sig)//'_abs',msg,newunit=abs_unt,status='replace',action="write")/=0) then
   ABI_ERROR(msg)
 end if
 write(abs_unt,'(a)')'#energy(eV),nomega,komega,refl.,abso.(cm-1)'

 ABI_MALLOC(fct,(npt))
 ABI_MALLOC(fct2,(npt))
 ABI_MALLOC(fct3,(npt))
 ABI_MALLOC(fct4,(npt))
 ABI_MALLOC(fct5,(npt))
 ABI_MALLOC(fp,(npt))
 ABI_MALLOC(fpp,(npt))
 ABI_MALLOC(fppp,(npt))
 ABI_MALLOC(x1,(npt))
 ABI_MALLOC(x2,(npt))
 ABI_MALLOC(fct1,(npt))
 ABI_MALLOC(ppsig,(npt))
 ABI_MALLOC(fctii,(npt))
 ABI_MALLOC(abso,(npt))
 ABI_MALLOC(nomega,(npt))

!loop on the initial energy grid
 do ip=1,npti

!  adjust the interval before and after the pole to reflect range/npt interval
   xsum=zero
   dx=(xi(npti)-xi(1))/dble(npt-1)
   pole=xi(ip)
   npt1=int((pole-del)/dx)
   dx1=zero
   if(npt1/=1) dx1=(pole-del)/(npt1-1)
   npt2=int((xi(npti)-pole-del)/dx)
   dx2=(xi(npti)-pole-del)/(npt2-1)

!  for the moment skip the pp calculation when the pole if too close to the end of the range
   if (npt1<=1.or.npt2<=1) then

     xsum=zero
     ppsig(ip)=zero

   else

!    define the fct for which the pp calculation is needed using xi^2-pole^2 factorization
     fctii(1:npti) = zero
     fctii(1:npti)=fcti(1:npti)*pole/(xi(1:npti)+pole)

!    define the grid on each side of the pole x1 before x2 after
     do ii=1,npt1
       x1(ii)=dx1*dble(ii-1)
     end do
     do ii=1,npt2
       x2(ii)=pole+del+dx2*dble(ii-1)
     end do

!    interpolate the initial fct fii on the new grids x1 and x2 (cubic spline)
!    write(std_out,*) npti,npt1

!    MJV 6/12/2008:
!    For each use of fctii should ensure that npt1 npt2 etc... are less than
!    npt=len(fctii)
!    TODO: move to spline/splint routines with no memory limitation
     call intrpl(npti,xi,fctii,npt1,x1,fct4,fct1,fct5,1)
     call intrpl(npti,xi,fctii,npt2,x2,fct3,fct2,fct5,1)

!    calculate the two integrals from 0-->pole-lamda and pole+lamda--> end range
!    trapezoidal integration
     do ii=1,npt1
       fct1(ii)=fct4(ii)/(x1(ii)-pole)
     end do
     do ii=1,npt2
       fct2(ii)=fct3(ii)/(x2(ii)-pole)
     end do

     do ii=2,npt1-1
       xsum=xsum+fct1(ii)*dx1
     end do
     do ii=2,npt2-1
       xsum=xsum+fct2(ii)*dx2
     end do
     xsum=xsum+half*(fct1(1)+fct1(npt1))*dx1+half*(fct2(1)+fct2(npt2))*dx2

!    Calculate the first and third derivative at the pole and add the taylor expansion
!    TODO: move to spline/splint routines with no memory limitation
     call intrpl(npti,xi,fctii,npti,xi,fct3,fct4,fct5,1)
     call intrpl(npti,xi,fct4,1,(/pole/),fp,fpp,fppp,1)

     idel=two*fp(1)*(del)+fppp(1)*(del**3)/nine
     xsum=xsum+idel

   end if

!  Calculate the derivated optical quantities and output the value
   sigma2=(-two/pi)*xsum
   eps1=one-(four_pi*sigma2/(pole))
   eps2=four*fcti(ip)*pi/(pole)

!  A special treatment of the case where eps2 is very small compared to eps1 is needed
   if(eps2**2 > eps1**2 * tol12)then
     nomega(ip)=sqrt(half*(eps1 + sqrt(eps1**2 + eps2**2)))
     komega=sqrt(half*(-eps1 + sqrt(eps1**2 + eps2**2)))
     abso(ip)=four_pi*fcti(ip)*ohmtosec*Ohmcm/nomega(ip)/(Sp_Lt_SI*100._dp)
   else if(eps1>zero)then
     nomega(ip)=sqrt(half*(eps1 + sqrt(eps1**2 + eps2**2)))
     komega=half*abs(eps2/sqrt(eps1))
     abso(ip)=four_pi*fcti(ip)*ohmtosec*Ohmcm/nomega(ip)/(Sp_Lt_SI*100._dp)
   else if(eps1<zero)then
     nomega(ip)=half*abs(eps2/sqrt(-eps1))
     komega=sqrt(half*(-eps1 + sqrt(eps1**2 + eps2**2)))
     abso(ip)=two*sqrt(-eps1)*pole*ohmtosec*Ohmcm/(Sp_Lt_SI*100._dp)
   end if

   refl=((one-nomega(ip))**2 + komega**2)/ &
&   ((one+nomega(ip))**2 + komega**2)

   write(eps_unt,'(5e18.10)') Ha_eV*pole,fcti(ip)*Ohmcm,sigma2*Ohmcm,eps1,eps2
   write(abs_unt,'(5e18.10)') Ha_eV*pole,nomega(ip),komega,refl,abso(ip)

 end do

 close(eps_unt)
 close(abs_unt)

 ABI_FREE(fct)
 ABI_FREE(x1)
 ABI_FREE(x2)
 ABI_FREE(fct2)
 ABI_FREE(fct3)
 ABI_FREE(fct4)
 ABI_FREE(fct5)
 ABI_FREE(fp)
 ABI_FREE(fpp)
 ABI_FREE(fppp)
 ABI_FREE(fct1)
 ABI_FREE(ppsig)
 ABI_FREE(fctii)
 ABI_FREE(abso)
 ABI_FREE(nomega)

end subroutine msig
!!***

end module m_conducti
!!***
