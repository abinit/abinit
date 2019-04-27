!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_iowf
!! NAME
!! m_iowf
!!
!! FUNCTION
!! Procedures for the IO of the WFK file.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, AR, MB, MVer, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_iowf

 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_abicore
 use m_errors
 use m_xmpi
 use m_wffile
 use m_abi_etsf
 use m_nctk
 use m_wfk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr
 use m_ebands

 use m_time,           only : cwtime, cwtime_report, timab
 use m_io_tools,       only : get_unit, flush_unit, iomode2str
 use m_fstrings,       only : endswith, sjoin
 use m_numeric_tools,  only : mask2blocks
 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_cgtools,        only : cg_zcopy
 use m_crystal,        only : crystal_t
 use m_rwwf,           only : rwwf
 use m_mpinfo,         only : proc_distrb_cycle
 use m_vkbr,           only : calc_vkb
 use m_wvl_rwwf,       only : wvl_write

 implicit none

 private

 public :: outwf

!!***

CONTAINS  !====================================================================================================
!!***

!!****f* m_iowf/outwf
!! NAME
!! outwf
!!
!! FUNCTION
!! Conduct output of a "wave-functions" file.
!!  - Compute the maximal residual
!!  - Then open a permanent file wff2 for final output of wf data
!!  - Create a new header for the file.
!!  - Write wave-functions (and energies)
!!
!! INPUTS
!!  cg(2,mcg)=wavefunction array (storage if nkpt>1)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  eigen( (2*mband)**response *mband*nkpt*nsppol)= eigenvalues (hartree) for all bands at each k point
!!  filnam= character string giving the root to form the name of the
!!   output WFK or WFQ file if response==0, otherwise it is the filename.
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kptns(3,nkpt)=k points in terms of recip primitive translations
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mkmem=Number of k-points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum number of plane waves
!!  natom=number of atoms in unit cell
!!  nband=number of bands
!!  nkpt=number of k points
!!  npwarr(nkpt)=number of planewaves in basis and on boundary for each k
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nstep=desired number of electron iteration steps
!!  occ(mband*nkpt*nsppol)=occupations for all bands at each k point
!!  resid(mband*nkpt*nsppol)=squared residuals for each band and k point
!!   where resid(n,k)=|<C(n,k)|(H-e(n,k))|C(n,k)>|^2 for the ground state
!!  response: if == 0, GS wavefunctions , if == 1, RF wavefunctions
!!  unwff2=unit for output of wavefunction
!!  wfs <type(wvl_projector_type)>=wavefunctions information for wavelets.
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! * The name of the file wff2 might be the same as that of the file wff1.
!!
!! PARENTS
!!      berryphase_new,dfpt_looppert,gstate
!!
!! CHILDREN
!!      xmpi_sum_master
!!
!! SOURCE

subroutine outwf(cg,dtset,psps,eigen,filnam,hdr,kg,kptns,mband,mcg,mkmem,&
&                mpi_enreg,mpw,natom,nband,nkpt,npwarr,&
&                nsppol,occ,resid,response,unwff2,&
&                wfs,wvl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mkmem,mpw,natom,nkpt,nsppol,response,unwff2
!integer,intent(in) :: nstep
 character(len=*),intent(in) :: filnam
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(hdr_type), intent(inout) :: hdr
 type(wvl_wf_type),intent(in) :: wfs
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 integer, intent(in) :: kg(3,mpw*mkmem),nband(nkpt*nsppol),npwarr(nkpt)
 real(dp), intent(inout) :: cg(2,mcg)
 real(dp), intent(in) :: eigen((2*mband)**response*mband*nkpt*nsppol),kptns(3,nkpt)
 real(dp), intent(in) :: occ(mband*nkpt*nsppol),resid(mband*nkpt*nsppol)

!Local variables-------------------------------
 integer,parameter :: nkpt_max=50
 integer :: iomode,action,band_index,fform,formeig,iband,ibdkpt,icg,iat,iproj
 integer :: ierr,ii,ikg,ikpt,spin,master,mcg_disk,me,me0,my_nspinor
 integer :: nband_k,nkpt_eff,nmaster,npw_k,option,rdwr,sender,source !npwtot_k,
 integer :: spaceComm,spaceComm_io,spacecomsender,spaceWorld,sread,sskip,tim_rwwf,xfdim2
#ifdef HAVE_MPI
 integer :: ipwnbd
#endif
 real(dp) :: residk,residm,resims,cpu,wall,gflops
 logical :: ihave_data,iwrite,iam_master,done
 character(len=500) :: msg
 type(wffile_type) :: wff2
 character(len=fnlen) :: path
!arrays
 integer,allocatable :: kg_disk(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cg_disk(:,:),eig_k(:),occ_k(:)
#ifdef HAVE_NETCDF
 integer :: ncid, ncerr, kg_varid, mpw_disk, npwk_disk, timrev
 real(dp),allocatable :: vkb(:,:,:),vkbd(:,:,:),vkbsign(:,:)
 type(crystal_t) :: crystal
#endif

! *************************************************************************
!For readability of the source file, define a "me" variable also in the sequential case

 DBG_ENTER("COLL")

 xfdim2 = natom+4
!Init mpi_comm
 spaceWorld= mpi_enreg%comm_cell
 spaceComm=spaceWorld
 spaceComm_io=xmpi_comm_self

 if (mpi_enreg%paral_kgb==1 ) spaceComm_io= mpi_enreg%comm_bandspinorfft
 if (mpi_enreg%paral_kgb==1 ) spaceComm= mpi_enreg%comm_cell

!Paral_kgb=1 and Fortran-I/O is not supported (only for testing purpose)
 if (mpi_enreg%paral_kgb==1.and.dtset%iomode==IO_MODE_FORTRAN) then
   spaceWorld=mpi_enreg%comm_kpt
   write(msg,'(7a)') &
&   'WF file is written using standard Fortran I/O',ch10,&
&   'and Kpt-band-FFT parallelization is active !',ch10,&
&   'This is only allowed for testing purposes.',ch10,&
&   'The produced WF file will be incomplete and not useable.'
   MSG_WARNING(msg)
 end if

!If parallel HF calculation
 if (mpi_enreg%paral_hf==1 ) spaceComm_io= mpi_enreg%comm_hf
 if (mpi_enreg%paral_hf==1 ) spaceComm= mpi_enreg%comm_cell

!Paral_hf=1 and Fortran-I/O is not supported (copy from paral_kgb... not tested)
 if (mpi_enreg%paral_hf==1.and.dtset%iomode==IO_MODE_FORTRAN) then
   spaceWorld=mpi_enreg%comm_kpt
   write(msg,'(7a)') &
&   'WF file is written using standard Fortran I/O',ch10,&
&   'and HF parallelization is active !',ch10,&
&   'This is only allowed for testing purposes.',ch10,&
&   'The produced WF file will be incomplete and not useable.'
   MSG_WARNING(msg)
 end if

 ! check consistency between dimensions and input hdr.
 !ABI_CHECK(mband == maxval(hdr%nband), "hdr:mband")
 !ABI_CHECK(nkpt == hdr%nkpt, "hdr:nkpt")
 !ABI_CHECK(nsppol == hdr%nsppol, "hdr:nsppol")
 !ABI_CHECK(all(hdr%npwarr == npwarr), "hdr:npwarr")
 !ABI_CHECK(all(hdr%nband == nband), "hdr:nband")
 !ABI_CHECK(maxval(hdr%npwarr) == mpw, "hdr:nband")

!Init me
 me=mpi_enreg%me_kpt
 me0=me
!Define master
 master=0

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 tim_rwwf =0
 source = master
 sread = master
 iam_master=(master==me)
 iwrite=iam_master
 sender=-1

!Compute mean square and maximum residual over all bands and k points and spins
!(disregard k point weights and occupation numbers here)
 band_index=sum(nband(1:nkpt*nsppol))
 resims=sum(resid(1:band_index))/dble(band_index)

!Find largest residual over bands, k points, and spins, except for nbdbuf highest bands
!Already AVAILABLE in hdr ?!
 ibdkpt=1
 residm=zero
 do spin=1,nsppol
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(spin-1)*nkpt)
     nband_k=max(1,nband_k-dtset%nbdbuf)
     residm=max(residm,maxval(resid(ibdkpt:ibdkpt+nband_k-1)))
     ibdkpt=ibdkpt+nband_k
   end do
 end do

 write(msg,'(a,2p,e12.4,a,e12.4)')' Mean square residual over all n,k,spin= ',resims,'; max=',residm
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 band_index=0
 nkpt_eff=nkpt
 if( (dtset%prtvol==0 .or. dtset%prtvol==1) .and. nkpt_eff>nkpt_max ) nkpt_eff=nkpt_max

!Loop over spin again
 do spin=1,nsppol
!  Give (squared) residuals for all bands at each k
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(spin-1)*nkpt)
!    Will not print all residuals when prtvol=0 or 1
     if(ikpt<=nkpt_eff)then
!      Find largest residual over all bands for given k point
       residk=maxval(resid(1+band_index:nband_k+band_index))
       write(msg,'(1x,3f8.4,3x,i2,1p,e13.5,a)')kptns(1:3,ikpt),spin,residk,' kpt; spin; max resid(k); each band:'
       if(dtset%prtvol>=2)then
         call wrtout(ab_out,msg,'COLL')
       endif
       call wrtout(std_out,msg,'COLL')
       do ii=0,(nband_k-1)/8
         write(msg,'(1x,1p,8e9.2)')(resid(iband+band_index),iband=1+ii*8,min(nband_k,8+ii*8))
         if(dtset%prtvol>=2)then
           call wrtout(ab_out,msg,'COLL')
         endif
         call wrtout(std_out,msg,'COLL')
       end do
     else if(ikpt==nkpt_eff+1)then
       write(msg,'(2a)')' outwf : prtvol=0 or 1, do not print more k-points.',ch10
       if(dtset%prtvol>=2)then
         call wrtout(ab_out,msg,'COLL')
       endif
       call wrtout(std_out,msg,'COLL')
     end if
     band_index=band_index+nband_k
   end do
 end do

!Will write the wavefunction file only when nstep>0
!MT 07 2015: writing reactivated when nstep=0
!if (nstep>0 .and. dtset%prtwf/=0) then
 if (dtset%prtwf/=0) then

   ! Only the master write the file, except if MPI I/O, but the
   ! full wff dataset should be provided to WffOpen in this case
   iomode=IO_MODE_FORTRAN_MASTER
   if (dtset%iomode==IO_MODE_MPI)  iomode = IO_MODE_MPI
   if (dtset%iomode==IO_MODE_ETSF) iomode = IO_MODE_ETSF

#ifdef HAVE_NETCDF
   if (dtset%iomode == IO_MODE_ETSF .and. dtset%usewvl == 0) then
     call cg_ncwrite(filnam,hdr,dtset,response,mpw,mband,nband,nkpt,nsppol,&
       dtset%nspinor,mcg,mkmem,eigen,occ,cg,npwarr,kg,mpi_enreg,done)
#ifdef HAVE_NETCDF_DEFAULT
     ABI_CHECK(done, "cg_ncwrite must handle the output of the WFK file.")
#endif

     ! Write KB form factors. Only master works. G-vectors are read from file to avoid
     ! having to deal with paral_kgb distribution.
     if (me == master .and. dtset%prtkbff == 1 .and. dtset%iomode == IO_MODE_ETSF .and. dtset%usepaw == 0) then
       ABI_CHECK(done, "cg_ncwrite was not able to generate WFK.nc in parallel. Perphase hdf5 is not working")
       path = nctk_ncify(filnam)
       call wrtout(std_out, sjoin("Writing KB form factors to:", path))
       NCF_CHECK(nctk_open_modify(ncid, path, xmpi_comm_self))
       NCF_CHECK(nf90_inq_varid(ncid, "reduced_coordinates_of_plane_waves", kg_varid))
       mpw_disk = maxval(hdr%npwarr)

       ! Dimensions needed by client code to allocate memory when reading.
       ncerr = nctk_def_dims(ncid, [ &
         nctkdim_t("mproj", psps%mproj), &
         nctkdim_t("mpsang", psps%mpsang), &
         nctkdim_t("mpssoang", psps%mpssoang), &
         nctkdim_t("lnmax", psps%lnmax), &
         nctkdim_t("lmnmax", psps%lmnmax) &
       ])
       NCF_CHECK(ncerr)

       ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "mpspso"])
       NCF_CHECK(ncerr)

       ! Write indlmn table (needed to access vkb arrays)
       ncerr = nctk_def_arrays(ncid, [ &
         nctkarr_t("indlmn", "int", "six, lmnmax, number_of_atom_species"), &
         nctkarr_t("vkbsign", "dp", "lnmax, number_of_atom_species"), &
         nctkarr_t("vkb", "dp", "max_number_of_coefficients, lnmax, number_of_atom_species, number_of_kpoints"), &
         nctkarr_t("vkbd", "dp", "max_number_of_coefficients, lnmax, number_of_atom_species, number_of_kpoints") &
       ], defmode=.True.)
       NCF_CHECK(ncerr)

       ! Switch to write mode.
       NCF_CHECK(nctk_set_datamode(ncid))
       NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "indlmn"), psps%indlmn))

       ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
         "mpspso"], &
         [psps%mpspso])
       NCF_CHECK(ncerr)

       ! Calculate KB form factors and derivatives.
       ! The arrays are allocated with lnmax to support pseudos with more than projector.
       ! Note that lnmax takes into account lloc hence arrays are in packed form and one should be
       ! accessed with the indices provided by psps%indlmn.
       ABI_MALLOC(vkbsign, (psps%lnmax, psps%ntypat))
       ABI_MALLOC(vkb, (mpw_disk, psps%lnmax, psps%ntypat))
       ABI_MALLOC(vkbd, (mpw_disk, psps%lnmax, psps%ntypat))
       ABI_MALLOC(kg_disk, (3, mpw_disk))

       timrev = 2 ! FIXME: Use abinit convention for timrev
       crystal = hdr_get_crystal(hdr, timrev)

       ! For each k-point: read full G-vector list from file, compute KB data and write to file.
       do ikpt=1,nkpt
         npwk_disk = hdr%npwarr(ikpt)
         NCF_CHECK(nf90_get_var(ncid, kg_varid, kg_disk, start=[1, 1, ikpt], count=[3, npwk_disk, 1]))
         vkb = zero; vkbd = zero
         call calc_vkb(crystal, psps, kptns(:, ikpt), npwk_disk, mpw_disk, kg_disk, vkbsign, vkb, vkbd)

         if (ikpt == 1) then
           ! This for the automatic tests.
           NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vkbsign"), vkbsign))
           if(dtset%prtvol>=2)then
             write(msg,'(a)') 'prtkbff: writing first and last G-components of the KB form factors'
             call wrtout(ab_out,msg,'COLL')
             do iat=1,psps%ntypat
               write(msg,'(a10,i5)') 'atype ',iat
               call wrtout(ab_out,msg,'COLL')
               do iproj=1,psps%lnmax
                 write(msg,'(a10,i5,a,a10,e12.4,a,2(a10,2e12.4,a))') &
                        'projector ', iproj,ch10, &
                        'vkbsign   ', vkbsign(iproj,iat), ch10, &
                        'vkb       ', vkb(1,iproj,iat),  vkb(npwk_disk,iproj,iat), ch10, &
                        'vkbd      ', vkbd(1,iproj,iat), vkbd(npwk_disk,iproj,iat), ''
                 call wrtout(ab_out,msg,'COLL')
               end do
             end do
           end if
         end if

         NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vkb"), vkb, start=[1, 1, 1, ikpt]))
         NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vkbd"), vkbd, start=[1, 1, 1, ikpt]))
       end do
       NCF_CHECK(nf90_close(ncid))

       ABI_FREE(kg_disk)
       ABI_FREE(vkbsign)
       ABI_FREE(vkb)
       ABI_FREE(vkbd)
       call crystal%free()
     end if

     if (done) return
     ! If cg_ncwrite cannot handle the IO because HDF5 + MPI-IO support is missing, we fallback to Fortran + MPI-IO.
     msg = "Could not produce a netcdf file in parallel (MPI-IO support is missing). Will fallback to MPI-IO with Fortran"
     MSG_WARNING(msg)
     iomode=IO_MODE_MPI
   end if
#endif

   call cwtime(cpu, wall, gflops, "start")
   call wrtout(std_out, sjoin(ch10,'  outwf: writing wavefunctions to:', trim(filnam), "with iomode:", iomode2str(iomode)))

   ! Create an ETSF file for the wavefunctions
   if (iomode == IO_MODE_ETSF) then
     ABI_CHECK(xmpi_comm_size(spaceComm) == 1, "Legacy etsf-io code does not support nprocs > 1")
#ifdef HAVE_ETSF_IO
     call abi_etsf_init(dtset, filnam, 2, .true., hdr%lmn_size, psps, wfs)
     !crystal = hdr_get_crystal(hdr, 2)
     !NCF_CHECK(crystal%ncwrite_path(nctk_ncify(filnam)))
     !call crystal%free()
     !ncerr = ebands_ncwrite_path(gs_ebands, filname, ncid)
     !NCF_CHECK(ncerr)
#else
     MSG_ERROR("ETSF_IO is not activated")
     ABI_UNUSED(psps%ntypat)
#endif
   end if

   call WffOpen(iomode,spaceComm,filnam,ierr,wff2,master,me0,unwff2,spaceComm_io)
   ! Conduct wavefunction output to wff2

   ABI_ALLOCATE(kg_disk,(3,mpw))

   mcg_disk=mpw*my_nspinor*mband
   formeig=0; if (response==1) formeig=1

   ABI_ALLOCATE(eig_k,( (2*mband)**formeig * mband))
   ABI_ALLOCATE(occ_k,(mband))

#ifdef HAVE_MPI
   call xmpi_barrier(spaceComm)
!  Compute mband and mpw
   ABI_MALLOC_OR_DIE(cg_disk,(2,mcg_disk), ierr)
#endif

   band_index=0
   icg=0
   if(mpi_enreg%paralbd==0) tim_rwwf=6
   if(mpi_enreg%paralbd==1) tim_rwwf=12

!  Write header info for new wf file
   rdwr=2
   if (dtset%usewvl==0) then
     fform=2
   else
     fform = 200 ! Use 200 as radical for naming file format used by wavelets.
   end if

   if (wff2%iomode < 2) then
     call hdr_io(fform,hdr,rdwr,wff2)
     call WffKg(wff2,1)
   else if (wff2%iomode==IO_MODE_ETSF .and. iam_master) then
#ifdef HAVE_NETCDF
     NCF_CHECK(hdr_ncwrite(hdr, wff2%unwff, fform, nc_define=.True.))
#endif
   end if

   do spin=1,nsppol
     ikg=0

     do ikpt=1,nkpt
       nband_k=nband(ikpt+(spin-1)*nkpt)
       npw_k=npwarr(ikpt)

#ifdef HAVE_MPI
       if (dtset%usewvl == 0) then
         call xmpi_barrier(spaceWorld)

!        Must transfer the wavefunctions to the master processor
!        Separate sections for paralbd=1 or other values ; might be merged
         if(mpi_enreg%paralbd==0)then
           nmaster=0
           source=minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,spin))
           ihave_data=.false.
           if(source==me)ihave_data=.true.
           action=0
!          I am the master node, and I have the data in cg or cg_disk
           if((iam_master).and.(ihave_data))action=1
!          I am not the master, and I have the data => send to master
           if((.not.iam_master).and.(ihave_data))action=2
!          I am the master, and I receive the data
           if((iam_master).and.(.not.ihave_data))action=3

!          I have the data in cg or cg_disk ( MPI_IO case)
           if (iomode==IO_MODE_MPI) then
             action = 0
             sender=-1
             iwrite=.false.
             if (ihave_data)then
               action=1
               iwrite=.true.
               sender=me
             end if
           end if

!          I am the master node, and I have the data in cg or cg_disk
!          I have the data in cg or cg_disk ( MPI_IO case)
           if(action==1)then
!            Copy from kg to kg_disk
             kg_disk(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
!            Copy from cg to cg_disk
             do ipwnbd=1,nband_k*npw_k*my_nspinor
               cg_disk(1,ipwnbd)=cg(1,ipwnbd+icg)
               cg_disk(2,ipwnbd)=cg(2,ipwnbd+icg)
             end do
           end if

!          I am not the master, and I have the data => send to master
!          I am the master, and I receive the data
           if ( action==2.or.action==3) then
             !write(std_out,*)npw_k,nband_k
             call timab(48,1,tsec)
             if(action==2)then
               call xmpi_exch(kg(:,1+ikg:npw_k+ikg),3*npw_k,source,kg_disk,nmaster,spaceWorld,ierr)
               call xmpi_exch(cg(:,icg+1:icg+nband_k*npw_k*my_nspinor),2*nband_k*npw_k*my_nspinor, &
&               source,cg_disk,nmaster,spaceWorld,ierr)
             else
               call xmpi_exch(kg_disk,3*npw_k,source,kg_disk,nmaster,spaceWorld,ierr)
               call xmpi_exch(cg_disk,2*nband_k*npw_k*my_nspinor,source,cg_disk,nmaster,spaceWorld,ierr)
             end if
             call timab(48,2,tsec)
           end if


         else if(mpi_enreg%paralbd==1)then
           nmaster=0
#ifdef HAVE_MPI_IO
           sender=IO_MODE_FORTRAN_MASTER
           if( iomode==IO_MODE_MPI) then
             nmaster=mpi_enreg%proc_distrb(ikpt,1,spin)
             sender=nmaster
           end if
#endif

!          Note the loop over bands
           do iband=1,nband_k

!            The message passing related to kg is counted as one band
             action=0

!            I am the master node, and I have the data in cg or cg_disk
             if( mpi_enreg%proc_distrb(ikpt,iband,spin)==nmaster .and. me==nmaster) then
               action=1
!              I am not the master, and I have the data => send to master
             elseif( mpi_enreg%proc_distrb(ikpt,iband,spin)==me .and. me/=nmaster ) then
               action = 2
!              I am the master, and I receive the data
             elseif( mpi_enreg%proc_distrb(ikpt,iband,spin)/=me .and. me==nmaster ) then
               action=3
             end if

             if(action==1) then
!              I am the master node, and I have the data in cg or cg_disk
!              Copy from kg to kg_disk
               if(iband==1)kg_disk(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
!              Copy from cg to cg_disk
               do ipwnbd=1,npw_k*my_nspinor
                 cg_disk(1,(iband-1)*npw_k*my_nspinor+ipwnbd) = cg(1,(iband-1)*npw_k*my_nspinor+ipwnbd+icg)
                 cg_disk(2,(iband-1)*npw_k*my_nspinor+ipwnbd) = cg(2,(iband-1)*npw_k*my_nspinor+ipwnbd+icg)
               end do
             end if  ! action=1

             if ( action==2.or.action==3) then
!              action=2 :  I am not the master, and I have the data => send to master
!              action=3 :  I am the master, and I receive the data
               call timab(48,1,tsec)
               if ( iband == 1 ) then
                 if (action==2) then
                   call xmpi_exch(kg(:,1+ikg:npw_k+ikg),3*npw_k,mpi_enreg%proc_distrb(ikpt,iband,spin), &
&                   kg_disk,nmaster,spaceWorld,ierr)
                 else
                   call xmpi_exch(kg_disk,3*npw_k,mpi_enreg%proc_distrb(ikpt,iband,spin),  &
&                   kg_disk,nmaster,spaceWorld,ierr)
                 end if
               end if       ! iband =1
               ipwnbd=(iband-1)*npw_k*my_nspinor
               if (action==2) then
                 call xmpi_exch( cg(:,ipwnbd+icg+1:ipwnbd+icg+npw_k*my_nspinor),2*npw_k*my_nspinor &
&                 ,mpi_enreg%proc_distrb(ikpt,iband,spin)                    &
&                 ,cg_disk(:,ipwnbd+1:ipwnbd+npw_k*my_nspinor),nmaster,spaceWorld,ierr)
               else
                 call xmpi_exch( cg_disk(:,ipwnbd+1:ipwnbd+npw_k*my_nspinor),2*npw_k*my_nspinor    &
&                 ,mpi_enreg%proc_distrb(ikpt,iband,spin)                    &
&                 ,cg_disk(:,ipwnbd+1:ipwnbd+npw_k*my_nspinor),nmaster,spaceWorld,ierr)
               end if

               call timab(48,2,tsec)
             end if        ! action=2 or action=3

             if(iomode==IO_MODE_MPI) then
!              I have the data in cg or cg_disk
               iwrite=.false.
               if (nmaster == me) iwrite=.true.
             end if

           end do ! End of loop over bands
         end if ! End of paralbd=1
       end if
#endif

!      Only the master will write to disk the final output wf file.
!      in MPI_IO case only iwrite will write to disk the final output wf file.
       if(iwrite) then
!        write(std_out,*) 'outwf : I am master and will write wf file'
         if(formeig==0)then
           eig_k(1:nband_k)=eigen(1+band_index:nband_k+band_index)
           occ_k(1:nband_k)=occ(1+band_index:nband_k+band_index)
         else
           eig_k(1:2*nband_k*nband_k)=eigen(1+band_index:2*nband_k*nband_k+band_index)
         end if
         option=2
         if(dtset%prtwf==3)option=5
!        if (dtset%prtwf == 2 .and. mkmem/=0) option=4

         if (dtset%usewvl == 0) then
#ifdef HAVE_MPI
           call rwwf(cg_disk,eig_k,formeig,0,0,ikpt,spin,kg_disk,mband,mcg_disk,mpi_enreg, &
&           nband_k, nband_k,npw_k,my_nspinor,occ_k,option,1,tim_rwwf,wff2)

#else
           kg_disk(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
           call rwwf(cg,eig_k,formeig,0,icg,ikpt,spin,kg_disk,mband,mcg,mpi_enreg,nband_k, &
&           nband_k, npw_k,my_nspinor,occ_k,option,1,tim_rwwf,wff2)
#endif
         else
           call wvl_write(dtset,eigen,mpi_enreg,option,hdr%rprimd,wff2,wfs,wvl,hdr%xred)
         end if
       end if

!      The wavefunctions for the present k point and spin are written
       if(response==0)band_index=band_index+nband_k
       if(response==1)band_index=band_index+2*nband_k*nband_k

       sskip=1
#ifdef HAVE_MPI
       if (dtset%usewvl == 0) then
         sskip=0
         if(.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,spin,me)))sskip=1
       end if
#endif
       if(sskip==1)then
         icg=icg+npw_k*my_nspinor*nband_k
         ikg=ikg+npw_k
       end if


#ifdef HAVE_MPI_IO
       spacecomsender=spaceComm
       if (mpi_enreg%paral_kgb==1) spacecomsender =mpi_enreg%comm_kpt
       if (mpi_enreg%paral_hf==1) spacecomsender =mpi_enreg%comm_kpt
       call WffOffset(wff2,sender,spacecomsender,ierr)
#endif

     end do ! ikpt
   end do ! spin

   ABI_DEALLOCATE(kg_disk)
#ifdef HAVE_MPI
   ABI_DEALLOCATE(cg_disk)
#endif

   ABI_DEALLOCATE(eig_k)
   ABI_DEALLOCATE(occ_k)

!  Close the wavefunction file (and do NOT delete it !)
   !if (wff2%iomode /= IO_MODE_NETCDF) then
   call WffClose(wff2,ierr)
   !end if

   call cwtime_report(" WFK output", cpu, wall, gflops)
 end if ! End condition of nstep>0

 ! Block here because we might need to read the WFK file in the caller.
 call xmpi_barrier(mpi_enreg%comm_cell)

 DBG_EXIT("COLL")

end subroutine outwf
!!***

!----------------------------------------------------------------------

!!****f* m_iowf/cg_ncwrite
!! NAME
!! cg_ncwrite
!!
!! FUNCTION
!! Conduct output of a "wave-functions" file with netcdf
!!
!! INPUTS
!!  fname= character string giving the root to form the name of the
!!   output WFK or WFQ file if response==0, otherwise it is the filename.
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  response: if == 0, GS wavefunctions , if == 1, RF wavefunctions
!!  mpw=maximum number of plane waves
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mkmem=maximum number of k-points treated by this node
!!  nkpt=number of k points
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  eigen((2*mband)**response *mband*nkpt*nsppol)= eigenvalues (hartree) for all bands at each k point
!!  occ(mband*nkpt*nsppol)=occupations for all bands at each k point
!!  cg(2,mcg)=wavefunction array (storage if nkpt>1)
!!  npwarr(nkpt)=number of planewaves in basis and on boundary for each k
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mpi_enreg=information about MPI parallelization
!!
!! OUTPUT
!!  done=.True if cg_ncwrite can handle the output of the WFK file in parallel.
!!
!! PARENTS
!!      m_iowf
!!
!! CHILDREN
!!      xmpi_sum_master
!!
!! SOURCE

#ifdef HAVE_NETCDF

subroutine cg_ncwrite(fname,hdr,dtset,response,mpw,mband,nband,nkpt,nsppol,nspinor,mcg,&
                      mkmem,eigen,occ,cg,npwarr,kg,mpi_enreg,done)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: response,mband,mcg,mkmem,mpw,nkpt,nsppol,nspinor
 character(len=*),intent(in) :: fname
 logical,intent(out) :: done
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(in) :: mpi_enreg
 type(hdr_type),intent(in) :: hdr
!arrays
 integer, intent(in) :: nband(nkpt*nsppol),kg(3,mpw*mkmem),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: eigen((2*mband)**response*mband*nkpt*nsppol),occ(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,fform2=2
 integer :: ii,iomode,icg,iband,ikg,ikpt,spin,me_cell,me_kpt,me_band,me_spinor,my_nspinor,nband_k,npw_k
 integer :: comm_cell,comm_fft,comm_bandfft,formeig
 integer :: cnt,min_cnt,max_cnt,ierr,action,source,ncid,ncerr,cg_varid,kg_varid !,eig_varid,
 integer :: timrev,paral_kgb,npwtot_k !,start_pwblock !,start_cgblock !count_pwblock,
 integer :: ipw,ispinor_index,npwso,npwsotot,npwtot,nspinortot,ikpt_this_proc,ispinor
 integer :: bandpp,nproc_band,nproc_fft,nproc_spinor,me_fft,nproc_cell,nwrites
 integer :: comm_mpiio,nranks,bstart,bcount !nbdblock,nblocks,
 !integer :: band_blocksize,band
 real(dp) :: cpu,wall,gflops
 logical :: ihave_data,iam_master,single_writer,same_layout,use_collective
 character(len=500) :: msg
 character(len=fnlen) :: path
 type(wfk_t) :: wfk
 type(crystal_t) :: crystal
 type(ebands_t) :: gs_ebands
!arrays
 integer,allocatable :: kg_k(:,:),iter2kscgkg(:,:),ind_cg_mpi_to_seq(:),rank_has_cg(:),ranks_io(:)!,gblock(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: eigen3d(:,:,:),occ3d(:,:,:),cg_k(:,:),my_cgblock(:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")
 done = .False.

 path = nctk_ncify(fname)
 call wrtout(std_out, sjoin("in cg_ncwrite with path:", path), "COLL")

 ! communicators and ranks
 comm_cell = mpi_enreg%comm_cell; me_cell = mpi_enreg%me_cell; nproc_cell = mpi_enreg%nproc_cell
 comm_fft = mpi_enreg%comm_fft; me_fft  = mpi_enreg%me_fft; nproc_fft  = mpi_enreg%nproc_fft
 comm_bandfft = mpi_enreg%comm_bandfft
 me_kpt = mpi_enreg%me_kpt; me_band = mpi_enreg%me_band; me_spinor = mpi_enreg%me_spinor
 iam_master = (me_kpt == master)

 paral_kgb = dtset%paral_kgb
 nproc_band = mpi_enreg%nproc_band
 bandpp     = mpi_enreg%bandpp
 nproc_spinor = mpi_enreg%nproc_spinor

 ! FIXME
 my_nspinor = max(1, nspinor/nproc_spinor)
 if (nspinor == 2 .and. my_nspinor == 1) then
   MSG_ERROR("Spinor parallelization not coded yet")
 end if

 ! TODO: Be careful with response == 1 in parallel because the distribution of the cg
 ! can be **VERY** different from cg if nprocs > nkpt * nsppol
 !ABI_CHECK(response==0, "response == 1 not coded")
 if (size(hdr%nband) == size(nband)) then
   ABI_CHECK(all(Hdr%nband == nband),"nband")
 else
   MSG_ERROR("hdr%nband and nband have different size!")
 end if

 if (xmpi_comm_size(comm_cell) == 1) then
   ABI_CHECK(all(npwarr == hdr%npwarr), "npwarr != hdr%npwarr")
 end if

 ! FIXME: Use abinit convention for timrev
 timrev = 2
 crystal = hdr_get_crystal(hdr, timrev)

 ! TODO
 ! Be careful with response == 1.
 ! gs_ebands contains the GS eigenvalues and occupation and will be written if this is a
 ! GS wfk. If we have a DFPT file, eigen stores the GKK matrix element, in this case
 ! we don't write gs_ebands but we define new variables in the netcdf file to store the GKK
 ABI_MALLOC(occ3d, (mband,nkpt,nsppol))
 call unpack_eneocc(nkpt,nsppol,mband,nband,occ,occ3d)

 if (response == 0) then
   formeig = 0
   ABI_MALLOC(eigen3d, (mband,nkpt,nsppol))
   call unpack_eneocc(nkpt,nsppol,mband,nband,eigen,eigen3d)
   gs_ebands = ebands_from_hdr(hdr, mband, eigen3d); gs_ebands%occ = occ3d
   ABI_FREE(eigen3d)
 else
   formeig = 1
 end if

 ! same_layout is set to True if the interal representation of the cgs
 ! is compatible with the representation on file.
 iomode = IO_MODE_ETSF
 if (response == 0) then
    same_layout = (paral_kgb == 0 .or. (paral_kgb == 1 .and. all([nproc_fft, nproc_band, nproc_spinor] == 1)))
 else
   ! For the time being, these cases are not implemented in the DFPT part.
   ABI_CHECK(nproc_fft==1, "nproc_fft != 1 not coded")
   ABI_CHECK(nproc_band==1, "nproc_band != 1 not coded")
   ABI_CHECK(nproc_spinor==1, "nproc_spinor != 1 not coded")

   ! Note: It would be possible to use collective IO if the cg1 are block-distributed
   same_layout = .True.
   spin_loop: do spin=1,nsppol
     do ikpt=1,nkpt
       nband_k = nband(ikpt + (spin-1)*nkpt)
       if (.not. proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,spin,me_kpt)) then
         if (any(mpi_enreg%proc_distrb(ikpt,:nband_k,spin) /= me_kpt)) then
           same_layout = .False.; exit spin_loop
         end if
       end if
     end do
   end do spin_loop
   call xmpi_land(same_layout, comm_cell)
 end if

 call cwtime(cpu, wall, gflops, "start")

 if (same_layout) then
   single_writer = .True.
   if (nctk_has_mpiio) single_writer = .False.
   !single_writer = .True.

   write(msg,'(5a,l1)')&
     "same layout Writing WFK file: ",trim(path),", with iomode: ",trim(iomode2str(iomode)),", single writer: ",single_writer
   call wrtout(std_out,msg,'PERS', do_flush=.True.)

   if (.not. single_writer) then

#ifdef HAVE_NETCDF
     ! master opens the file and write the metadata.
     if (xmpi_comm_rank(comm_cell) == master) then
       ncerr = nf90_einval
#ifdef HAVE_NETCDF_MPI
       ncerr = nf90_create(path, cmode=ior(ior(nf90_netcdf4, nf90_mpiio), nf90_write), &
         comm=xmpi_comm_self, info=xmpio_info, ncid=ncid)
#endif
       NCF_CHECK_MSG(ncerr, sjoin("create_par: ", path))

       call wfk_ncdef_dims_vars(ncid, hdr, fform2, write_hdr=.True.)
       NCF_CHECK(crystal%ncwrite(ncid))

       if (response == 0) then
         NCF_CHECK(ebands_ncwrite(gs_ebands, ncid))
       else
         call ncwrite_eigen1_occ(ncid, nband, mband, nkpt, nsppol, eigen, occ3d)
       end if

       NCF_CHECK(nf90_close(ncid))
     end if

     ! Compute the table for collective IO.
     ABI_CALLOC(iter2kscgkg, (4, nkpt*nsppol))
     cnt = 0; icg = 0
     do spin=1,nsppol
       ikg = 0
       do ikpt=1,nkpt
         nband_k = nband(ikpt + (spin-1)*nkpt)
         npw_k = npwarr(ikpt)
         if (.not. proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,spin,me_kpt)) then
           ! FIXME v2[01], v3[6], with 4 procs
           ! v6[7] and v7[68], v7[69] fail but due to a extra line with nprocs
           ! v7[96] with np=4 seems to be more serious but it crashes also in trunk.
           ! v2[88] is slow likely due to outscfv
           ABI_CHECK(all(mpi_enreg%proc_distrb(ikpt,:nband_k,spin) == me_kpt), "bands are distributed")
           cnt = cnt + 1
           iter2kscgkg(:,cnt) = [ikpt, spin, icg, ikg]
           icg = icg + npw_k*my_nspinor*nband_k
           ikg = ikg + npw_k
         end if
       end do
     end do
     if (cnt == 0) then
       write(std_out,*)"cnt == 0 for me_cell, me, me_kpt",mpi_enreg%me_cell, mpi_enreg%me, mpi_enreg%me_kpt
       !ABI_CHECK(cnt > 0, "This processor does not have wavefunctions!")
     end if

     call xmpi_min(cnt, min_cnt, comm_cell, ierr)
     call xmpi_max(cnt, max_cnt, comm_cell, ierr)

     ! Handle idle procs, i.e. processors that do not have wavefunctions
     ! This happens if paral_kgb == 0 and nprocs > nkpt * nsppol (Abinit does not stop anymore!)
     comm_mpiio = comm_cell

     if (min_cnt <= 0) then
       MSG_COMMENT("Will create subcommunicator to exclude idle processors from MPI-IO collective calls")
       ABI_CHECK(paral_kgb == 0, "paral_kgb == 1 with idle processors should never happen")

       ! Prepare the call to xmpi_subcomm that will replace comm_mpiio.
       ABI_CALLOC(rank_has_cg, (0:nproc_cell-1))
       do spin=1,nsppol
         do ikpt=1,nkpt
           nband_k = nband(ikpt + (spin-1)*nkpt)
           if (.not. proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,spin,me_kpt)) then
              rank_has_cg(me_kpt) = 1
              exit
            end if
         end do
       end do

       call xmpi_sum(rank_has_cg,comm_cell,ierr)
       nranks = count(rank_has_cg == 1)
       ABI_MALLOC(ranks_io, (nranks))
       cnt = 0
       do ii=0,nproc_cell-1
         if (rank_has_cg(ii) == 1) then
           cnt = cnt + 1
           ranks_io(cnt) = ii
         end if
       end do
       !write(std_out,*)"nranks, ranks_io:", nranks, ranks_io
       comm_mpiio = xmpi_subcomm(comm_cell, nranks, ranks_io)
       if (.not. rank_has_cg(me_kpt) == 1) then
         comm_mpiio = xmpi_comm_null
         ABI_CHECK(rank_has_cg(me_kpt) == 0, "rank_has_cg must be 0 or 1")
       end if
       ABI_FREE(ranks_io)
       ABI_FREE(rank_has_cg)
     end if

     ! Open the file in parallel inside comm_mpiio.
     call xmpi_barrier(comm_cell)
     if (comm_mpiio == xmpi_comm_null) goto 100

     ncerr = nf90_einval
#ifdef HAVE_NETCDF_MPI
     ncerr = nf90_open(path, mode=ior(ior(nf90_netcdf4, nf90_mpiio), nf90_write),&
                       comm=comm_mpiio, info=xmpio_info, ncid=ncid)
#endif
     NCF_CHECK_MSG(ncerr, sjoin("open_par: ", path))

     ! Use individual IO (default) for the G-vectors [3, mpw, nkpt]
     NCF_CHECK(nf90_inq_varid(ncid, "reduced_coordinates_of_plane_waves", kg_varid))
     spin = 1; ikg=0
     do ikpt=1,nkpt
       npw_k = npwarr(ikpt)
       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,1,spin,me_kpt)) cycle
       ncerr = nf90_put_var(ncid, kg_varid, kg(:, 1+ikg:), start=[1,1,ikpt], count=[3,npw_k,1])
       NCF_CHECK_MSG(ncerr, "put_kg_k")
       ikg = ikg + npw_k
     end do

     NCF_CHECK(nf90_inq_varid(ncid, "coefficients_of_wavefunctions", cg_varid))

     use_collective = (response == 0) ! or (response == 1 .and. nproc_cell == 1)

     if (use_collective) then
       call wrtout(std_out,"Using collective IO for the CGs","COLL")
       ! Use collective IO for the CGs
       ncerr = nf90_einval
#ifdef HAVE_NETCDF_MPI
       ncerr = nf90_var_par_access(ncid, cg_varid, nf90_collective)
#endif
       NCF_CHECK(ncerr)

       do cnt=1,max_cnt
          ikpt = iter2kscgkg(1,cnt)
          spin = iter2kscgkg(2,cnt)
          icg  = iter2kscgkg(3,cnt)
          ikg  = iter2kscgkg(4,cnt)

          ! The array on file has shape [cplex, mpw, nspinor, mband, nkpt, nsppol]
          if (ikpt /= 0) then
            nband_k = nband(ikpt + (spin-1)*nkpt)
            npw_k   = npwarr(ikpt)

            !write(std_out,*)"About to write ikpt, spin, npw_k, icg: ",ikpt, spin, npw_k, icg
            ncerr = nf90_put_var(ncid, cg_varid, cg(:, 1+icg:), start=[1,1,1,1,ikpt,spin], &
               count=[2,npw_k,nspinor,nband_k,1,1])
            NCF_CHECK_MSG(ncerr, "writing cg")
          else
            ! This happens when nkpt * nsppol // nprocs != 0
            ! Note that we are using collective MPI-IO hence all processors must call put_var
            ! Here we re-write the ug(0) of the first (k-point, spin) treated by the node.
            ikpt = iter2kscgkg(1,1)
            spin = iter2kscgkg(2,1)
            ncerr = nf90_put_var(ncid, cg_varid, cg, start=[1,1,1,1,ikpt,spin], count=[1,1,1,1,1,1])
            NCF_CHECK_MSG(ncerr, "re-writing cg")
          end if
       end do

     else
       call wrtout(std_out,"Using individual IO for the CGs","COLL")
       ! Individual IO of the CGs (for debugging purposes)
       icg = 0
       do spin=1,nsppol
         do ikpt=1,nkpt
           nband_k = nband(ikpt + (spin-1)*nkpt)
           npw_k = npwarr(ikpt)
           if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,spin,me_kpt)) cycle
           ncerr = nf90_put_var(ncid, cg_varid, cg(:, 1+icg:), start=[1,1,1,1,ikpt,spin], count=[2,npw_k,nspinor,nband_k,1,1])
           NCF_CHECK_MSG(ncerr, "put_var")
           icg = icg+npw_k*my_nspinor*nband_k
         end do
       end do
     end if

     NCF_CHECK(nf90_close(ncid))

100  call xmpi_barrier(comm_cell)
     ABI_FREE(iter2kscgkg)

     done = .True.
     call cwtime_report(" collective ncwrite", cpu, wall, gflops)
#endif
! HAVE_NETCDF
   else ! single_writer
     if (nproc_cell > 1) then
       MSG_WARNING("Slow version without MPI-IO support. Processors send data to master...")
     else
       call wrtout(std_out, "Using netcdf library without MPI-IO support")
     end if

     ABI_MALLOC(kg_k,(3,mpw))
     ABI_MALLOC_OR_DIE(cg_k,(2,mpw*my_nspinor*mband), ierr)

     if (iam_master) then
       call wfk_open_write(wfk,hdr,path,formeig,iomode,get_unit(),xmpi_comm_self,write_hdr=.True.)

       NCF_CHECK(crystal%ncwrite(wfk%fh))
       !write(std_out,*)"after crystal_ncwrite"

       ! Write eigenvalues and occupations (these arrays are not MPI-distributed)
       if (response == 0) then
         NCF_CHECK(ebands_ncwrite(gs_ebands, wfk%fh))
       else
         call ncwrite_eigen1_occ(wfk%fh, nband, mband, nkpt, nsppol, eigen, occ3d)
       end if
     end if

     icg = 0
     do spin=1,nsppol
       ikg = 0
       do ikpt=1,nkpt
         nband_k = nband(ikpt + (spin-1)*nkpt)
         npw_k   = npwarr(ikpt)

         call xmpi_barrier(comm_cell)

         ! Transfer the wavefunctions and the g-vectors to the master processor
         source = minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,spin))
         ihave_data = (source==me_kpt)

         action=0
         if (iam_master .and. ihave_data)    action=1 ! I am the master node, and I have the data in cg
         if (.not.iam_master.and.ihave_data) action=2 ! I am not the master, and I have the data => send to master
         if (iam_master.and..not.ihave_data) action=3 ! I am the master, and I receive the data

         if (action==1) then ! Copy from kg and cg
           kg_k(:,1:npw_k) = kg(:,ikg+1:ikg+npw_k)
           call cg_zcopy(npw_k*my_nspinor*nband_k, cg(1,icg+1), cg_k)
         end if

         ! Exchange data
         if (action==2.or.action==3) then
           call timab(48,1,tsec)
           if (action==2) then
             call xmpi_exch(kg(:,1+ikg:npw_k+ikg),3*npw_k,source,kg_k,master,comm_cell,ierr)
             call xmpi_exch(cg(:,icg+1:icg+nband_k*npw_k*my_nspinor),2*nband_k*npw_k*my_nspinor,&
&             source,cg_k,master,comm_cell,ierr)
           else
             call xmpi_exch(kg_k,3*npw_k,source,kg_k,master,comm_cell,ierr)
             call xmpi_exch(cg_k,2*nband_k*npw_k*my_nspinor,source,cg_k,master,comm_cell,ierr)
           end if
           call timab(48,2,tsec)
         end if

         ! Master writes this block of bands.
         if (iam_master) then
           if (response == 0) then
             call wfk_write_band_block(wfk,[1,nband_k],ikpt,spin,xmpio_single,kg_k=kg_k,cg_k=cg_k,&
               eig_k=gs_ebands%eig(:,ikpt,spin),occ_k=gs_ebands%occ(:,ikpt,spin))
             !write(std_out,*)"cg_k",cg_k(:,1:2)
           else
             call wfk_write_band_block(wfk,[1,nband_k],ikpt,spin,xmpio_single,kg_k=kg_k,cg_k=cg_k)
           end if
         end if

         if (.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,spin,me_kpt))) then
           icg = icg+npw_k*my_nspinor*nband_k
           ikg = ikg+npw_k
         end if

       end do !ikpt
     end do !spin

     ABI_FREE(kg_k)
     ABI_FREE(cg_k)

     if (iam_master) call wfk_close(wfk)
     call xmpi_barrier(comm_cell)

     done = .True.
     call cwtime_report(" individual ncwrite", cpu, wall, gflops)
   end if

 else ! not same_layout

   if (nctk_has_mpiio) then
#ifdef HAVE_NETCDF
     call wrtout(std_out, &
       sjoin("scattered data. writing WFK file",trim(path),", with iomode: ",iomode2str(iomode)), 'PERS', do_flush=.True.)

     ! master write the metadata.
     if (xmpi_comm_rank(comm_cell) == master) then
       ncerr = nf90_einval
#ifdef HAVE_NETCDF_MPI
       ncerr = nf90_create(path, cmode=ior(ior(nf90_netcdf4, nf90_mpiio), nf90_write), &
         comm=xmpi_comm_self, info=xmpio_info, ncid=ncid)
#endif
       NCF_CHECK_MSG(ncerr, sjoin("create_par:", path))

       call wfk_ncdef_dims_vars(ncid, hdr, fform2, write_hdr=.True.)
       NCF_CHECK(crystal%ncwrite(ncid))

       ! Write eigenvalues and occupations (these arrays are not MPI-distributed)
       if (response == 0) then
         NCF_CHECK(ebands_ncwrite(gs_ebands, ncid))
       else
         call ncwrite_eigen1_occ(ncid, nband, mband, nkpt, nsppol, eigen, occ3d)
       end if

       NCF_CHECK(nf90_close(ncid))
     end if

     ! Reopen the file inside comm_cell
     call xmpi_barrier(comm_cell)
     ncerr = nf90_einval
#ifdef HAVE_NETCDF_MPI
     ncerr = nf90_open(path, mode=ior(ior(nf90_netcdf4, nf90_mpiio), nf90_write), &
       comm=comm_cell, info=xmpio_info, ncid=ncid)
#endif
     NCF_CHECK_MSG(ncerr, sjoin("create_par:", path))

     ! Get var ids
     NCF_CHECK(nf90_inq_varid(ncid, "reduced_coordinates_of_plane_waves", kg_varid))
     NCF_CHECK(nf90_inq_varid(ncid, "coefficients_of_wavefunctions", cg_varid))

     ! Write the G-vectors
     ikg = 0
     do ikpt=1,nkpt
       npw_k = npwarr(ikpt)
       npwtot_k = hdr%npwarr(ikpt)
       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,1,1,me_kpt)) cycle
       if (me_spinor /= 0) cycle
       !write(std_out,*)"In G-vector loop ",ikpt,", with me_cell me_kpt me_band, me_spinor ",me_cell,me_kpt,me_band,me_spinor
       !write(std_out,*)" with ik, npw_k, npwtot_k: ",ikpt, npw_k, npwtot_k

       ABI_MALLOC(ind_cg_mpi_to_seq, (npw_k))
       if (allocated(mpi_enreg%my_kgtab)) then
         ikpt_this_proc = mpi_enreg%my_kpttab(ikpt)
         ind_cg_mpi_to_seq = mpi_enreg%my_kgtab(1:npw_k,ikpt_this_proc)
       else
         ABI_CHECK(nproc_fft==1, "nproc_fft !=1 and my_kgtab not allocated")
         ind_cg_mpi_to_seq(1:npw_k) = [(ipw, ipw=1,npw_k)]
       end if

       ABI_CALLOC(kg_k, (3, npwtot_k))
       do ipw=1,npw_k
         kg_k(:, ind_cg_mpi_to_seq(ipw)) = kg(:,ikg+ipw)
       end do
       call xmpi_sum_master(kg_k,master,comm_bandfft,ierr)
       if (xmpi_comm_rank(comm_bandfft) == master) then
         ncerr = nf90_put_var(ncid, kg_varid, kg_k, start=[1,1,ikpt], count=[3,npwtot_k,1])
         NCF_CHECK_MSG(ncerr, "putting kg_k")
       end if
       ABI_FREE(kg_k)

       ! gblock contains block-distributed G-vectors inside comm_fft
       !call kg2seqblocks(npwtot_k,npw_k,kg(:,ikg+1:),ind_cg_mpi_to_seq,comm_fft,start_pwblock,count_pwblock,gblock)
       !write(std_out,*)"gblock(:,2)",gblock(:,2)
       !ncerr = nf90_put_var(ncid, kg_varid, gblock, start=[1,start_pwblock,ikpt], count=[3,count_pwblock,1])
       !NCF_CHECK_MSG(ncerr, "putting kg_k")
       !ABI_FREE(gblock)

       ABI_FREE(ind_cg_mpi_to_seq)

       ikg = ikg+npw_k
     end do

     ! Write wavefunctions
     if (response == 1) then
       ! The cg1 is allocated with size mcg1=mpw1*dtset%nspinor*dtset%mband*mk1mem_rbz*dtset%nsppol (see dfpt_looppert)
       ! hence bands are not MPI distributed, this is the reason why we have to use the loop over bands
       ! and the check on nwrites
       icg = 0
       do spin=1,nsppol
         do ikpt=1,nkpt
           nband_k = nband(ikpt + (spin-1)*nkpt)
           npw_k = npwarr(ikpt)
           npwtot_k = hdr%npwarr(ikpt)
           nwrites = 0
           do iband=1,nband_k
              if (mpi_enreg%proc_distrb(ikpt,iband,spin)/=me_kpt) cycle
              ! The coefficients_of_wavefunctions on file have shape [cplex, mpw, nspinor, mband, nkpt, nsppol]
              nwrites = nwrites + 1
              ii = 1 + (iband-1)*npw_k*my_nspinor + icg
              ncerr = nf90_put_var(ncid, cg_varid, cg(1:,ii:), start=[1,1,1,iband,ikpt,spin], count=[2,npwtot_k,nspinor,1,1,1])
              NCF_CHECK_MSG(ncerr, "put_var cg")
           end do ! iband
           if (nwrites /= 0) icg = icg + npw_k*my_nspinor*nband_k
         end do !ikpt
       end do !spin

     else
       icg = 0
       do spin=1,nsppol
         do ikpt=1,nkpt
           nband_k = nband(ikpt + (spin-1)*nkpt)
           npw_k = npwarr(ikpt)
           npwtot_k = hdr%npwarr(ikpt)

           if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,spin,me_kpt)) cycle
           !write(std_out,*)"In u(g)-vector loop ",ikpt,", with me_cell me_kpt me_band, me_spinor ",me_cell,me_kpt,me_band,me_spinor
           !write(std_out,*)"nband_k, npw_k, npwtot_k: ",nband_k, npw_k, npwtot_k

           ! Part taken from writewf
           ! Note that ind_cg_mpi_to_seq is wrong in nspinor > 1
           npwtot=npw_k; npwso=npw_k*nspinor
           npwsotot=npwso
           nspinortot = min(2,(1+mpi_enreg%paral_spinor)*nspinor)

           ABI_MALLOC(ind_cg_mpi_to_seq, (npwso))
           if (allocated(mpi_enreg%my_kgtab)) then
             ikpt_this_proc=mpi_enreg%my_kpttab(ikpt)
             do ispinor=1,nspinor
               ispinor_index=ispinor
               if (nproc_spinor > 1) ispinor_index = mpi_enreg%me_spinor + 1
               ind_cg_mpi_to_seq(1+npw_k*(ispinor-1):npw_k*ispinor)=npwtot*(ispinor_index-1) &
&               + mpi_enreg%my_kgtab(1:npw_k,ikpt_this_proc)
             end do
           else
             ABI_CHECK(nproc_fft==1, "nproc_fft !=1 and my_kgtab not allocated")
             ind_cg_mpi_to_seq(1:npwso) = [(ipw, ipw=1,npwso)]
           end if

           ! TODO: Blocking + collective IO in comm_cell
           !band_blocksize = nband_k/nbdblock
           !ibandmin = 1
           !step=min(ii,MAXBAND, nband_disk)
           !do iblock=1,nband_disk/step+1
           !   ibandmax = min(ibandmin+step-1, nband_disk)
           !   nband_block = ibandmax-ibandmin+1
           !gbase = icg + npw_k*my_nspinor*nband_k

           ! Each MPI proc write bcount bands with npwtot_k G-vectors starting from bstart.
           call cg2seqblocks(npwtot_k,npw_k,nband_k,cg(:,icg+1:),ind_cg_mpi_to_seq,comm_bandfft,bstart,bcount,my_cgblock)

           ! The coefficients_of_wavefunctions on file have shape [cplex, mpw, nspinor, mband, nkpt, nsppol]
           ncerr = nf90_put_var(ncid, cg_varid, my_cgblock, start=[1,1,1,bstart,ikpt,spin], count=[2,npwtot_k,nspinor,bcount,1,1])
           NCF_CHECK_MSG(ncerr, "put_var cg")

           ABI_FREE(my_cgblock)
           ABI_FREE(ind_cg_mpi_to_seq)

           icg = icg+npw_k*my_nspinor*nband_k
         end do !ikpt
       end do !spin
     end if

     NCF_CHECK(nf90_close(ncid))
     done = .True.

     call cwtime_report(" scattered ncwrite", cpu, wall, gflops)
#endif
! HAVE_NETCDF
   end if !nctk_has_mpiio
 end if

 call crystal%free()
 if (response == 0) call ebands_free(gs_ebands)

 ABI_FREE(occ3d)

 DBG_EXIT("COLL")

end subroutine cg_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_iowf/ncwrite_eigen1_occ
!! NAME
!! ncwrite_eigen1_occ
!!
!! FUNCTION
!!  Write the first order DFPT eigenvalues and the occupations.
!!
!! INPUTS
!!  ncid=Netcdf file handler.
!!  nband(nkpt*nsppol)=Number of bands.
!!  mband=maximum number of bands
!!  nkpt=Total number of k points
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  eigen((2*mband**2*nkpt*nsppol)= eigenvalues (hartree) for all bands at each k point
!!  occ3d(mband*nkpt*nsppol)=occupations for all bands at each k point
!!    Note that occ3d differes from the occ arrays used in the rest of the code.
!!    occ3d is an array with constant stride `mband` whereas abinit (for reasons that are not clear to me)
!!    packs the occupations in a 1d vector with a k-dependent separation (nband_k).
!!    occ and occ3d differ only if nband is k-dependent but you should never assume this, hence
!!    remember to *convert* occ into occ3d before calling this routine.
!!
!! PARENTS
!!      m_iowf
!!
!! CHILDREN
!!      xmpi_sum_master
!!
!! SOURCE

subroutine ncwrite_eigen1_occ(ncid, nband, mband, nkpt, nsppol, eigen, occ3d)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,mband,nkpt,nsppol
!arrays
 integer, intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: eigen(2*mband**2*nkpt*nsppol),occ3d(mband,nkpt,nsppol)

!Local variables-------------------------------
!scalars
 integer :: idx,spin,ikpt,nband_k,ib2,ib1
 integer :: ncerr,occ_varid,h1mat_varid
!arrays
 real(dp),allocatable :: h1mat(:,:,:,:,:)

! *************************************************************************

 ! Declare h1 array with abinit conventions
 ! Cannot use a 3D array since eigen1 are packed and one could have different number of bands
 ! per k-points. (this is not an official etsf-io variable)!
 ! elements in eigen are packed in the first positions.
 ! Remember that eigen are not MPI-distributed so no communication is needed
 idx=1
 ABI_CALLOC(h1mat, (2,mband,mband,nkpt,nsppol))
 do spin=1,nsppol
   do ikpt=1,nkpt
     nband_k = nband(ikpt + (spin-1)*nkpt)
     do ib2=1,nband_k
       do ib1=1,nband_k
          h1mat(:,ib1,ib2,ikpt,spin) = eigen(idx:idx+1)
          idx=idx+2
        end do
     end do
   end do
 end do

 NCF_CHECK(nctk_set_defmode(ncid))

 ncerr = nctk_def_arrays(ncid, nctkarr_t('h1_matrix_elements', "dp", &
&"complex, max_number_of_states, max_number_of_states, number_of_kpoints, number_of_spins"))
 NCF_CHECK(ncerr)

 NCF_CHECK(nf90_inq_varid(ncid, "occupations", occ_varid))
 NCF_CHECK(nf90_inq_varid(ncid, "h1_matrix_elements", h1mat_varid))

 ! Write data
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK_MSG(nf90_put_var(ncid, occ_varid, occ3d), "putting occ3d")
 NCF_CHECK_MSG(nf90_put_var(ncid, h1mat_varid, h1mat), "putting h1mat")

 ABI_FREE(h1mat)

end subroutine ncwrite_eigen1_occ
!!***

#endif
! HAVE_NETCDF

!----------------------------------------------------------------------

!!****f* m_iowf/kg2seqblocks
!! NAME
!! kg2seqblocks
!!
!! FUNCTION
!!
!! INPUTS
!!  npwtot_k
!!  npw_k=
!!  kg(3,npw_k)=reduced planewave coordinates.
!!  gmpi2seq
!!  comm_fft=FFT communicator
!!
!! OUTPUT
!!  start_pwblock
!!  count_pwblock
!!  gblock(:,:)
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_sum_master
!!
!! SOURCE

subroutine kg2seqblocks(npwtot_k,npw_k,kg_k,gmpi2seq,comm_fft,start_pwblock,count_pwblock,gblock)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwtot_k,npw_k,comm_fft
 integer,intent(out) :: start_pwblock,count_pwblock
!arrays
 integer,intent(in) :: kg_k(3,npw_k),gmpi2seq(npw_k)
 integer,allocatable,intent(out) :: gblock(:,:)

!Local variables-------------------------------
!scalars
 integer :: me_fft,nproc_fft,rank,ig,igseq,maxnpw,ierr
!arrays
 integer,allocatable :: igstart_rank(:),gbuf(:,:)

! *************************************************************************

 me_fft = xmpi_comm_rank(comm_fft); nproc_fft = xmpi_comm_size(comm_fft)

 ! Handle sequential case.
 if (nproc_fft == 1) then
   start_pwblock = 1; count_pwblock = npw_k
   ABI_MALLOC(gblock, (3, npw_k))
   gblock(:,:) = kg_k
   return ! DOH
 end if

 ABI_MALLOC(igstart_rank, (0:nproc_fft))
 igstart_rank = [(1 + (npwtot_k/nproc_fft)*rank, rank=0,nproc_fft-1), 1 + npwtot_k]

 ! Get max dimension for workspace array
 ! Cannot use npwtot_k / nproc_fft because G-vectors are not equally distributed.
 call xmpi_max(npw_k, maxnpw, comm_fft, ierr)
 ABI_MALLOC(gbuf, (3, maxnpw))

 do rank=0,nproc_fft-1
   start_pwblock = igstart_rank(rank)
   count_pwblock = igstart_rank(rank+1) - igstart_rank(rank)

   gbuf = 0
   do ig=1,npw_k
     igseq = gmpi2seq(ig)
     if (igseq >= igstart_rank(rank) .and. igseq < igstart_rank(rank+1)) then
       igseq = igseq - start_pwblock + 1
       !ABI_CHECK(igseq <= maxnpw, "boom")
       gbuf(:,igseq) = kg_k(:,ig)
     end if
   end do
   call xmpi_sum_master(gbuf,rank,comm_fft,ierr)

   if (me_fft == rank) then
     ABI_MALLOC_OR_DIE(gblock, (3, count_pwblock), ierr)
     gblock = gbuf(:, :count_pwblock)
   end if
 end do

 start_pwblock = igstart_rank(me_fft)
 count_pwblock = igstart_rank(me_fft+1) - igstart_rank(me_fft)

 ABI_FREE(gbuf)
 ABI_FREE(igstart_rank)

end subroutine kg2seqblocks
!!***

!----------------------------------------------------------------------

!!****f* m_iowf/cg2seqblocks
!! NAME
!! cg2seqblocks
!!
!! FUNCTION
!!
!! INPUTS
!!  npwtot_k
!!  npw_k=
!!  cg_k(2,npw_k*nband_k)
!!  gmpi2seq
!!  comm_fft=FFT communicator
!!
!! OUTPUT
!!  bstart
!!  bcount
!!
!! PARENTS
!!      m_iowf
!!
!! CHILDREN
!!      xmpi_sum_master
!!
!! SOURCE

subroutine cg2seqblocks(npwtot_k,npw_k,nband,cg_k,gmpi2seq,comm_bandfft,bstart,bcount,my_cgblock)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwtot_k,npw_k,nband,comm_bandfft
 integer,intent(out) :: bstart,bcount
!arrays
 integer,intent(in) :: gmpi2seq(npw_k)
 real(dp),intent(in) :: cg_k(2,npw_k,nband)
 real(dp),allocatable,intent(out) :: my_cgblock(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: me,nprocs,rank,ig,ib,igseq,ierr,nbb,nbmax,band
!arrays
 integer,allocatable :: bstart_rank(:)
 real(dp),allocatable :: cgbuf(:,:,:)

! *************************************************************************

 me = xmpi_comm_rank(comm_bandfft); nprocs = xmpi_comm_size(comm_bandfft)

 ! Handle sequential case.
 if (nprocs == 1) then
   bstart = 1; bcount = nband
   ABI_MALLOC_OR_DIE(my_cgblock, (2, npwtot_k, nband), ierr)
   my_cgblock = cg_k
   return ! DOH
 end if

 ABI_MALLOC(bstart_rank, (0:nprocs))
 bstart_rank = [(1 + (nband/nprocs)*rank, rank=0,nprocs-1), 1 + nband]

 ! Allocate MPI buffer (same size on each MPI proc)
 nbmax = 0
 do rank=0,nprocs-1
   nbmax = max(nbmax, bstart_rank(rank+1) - bstart_rank(rank))
 end do
 ABI_MALLOC_OR_DIE(cgbuf, (2, npwtot_k, nbmax), ierr)

 nbb = bstart_rank(me+1) - bstart_rank(me)
 ABI_MALLOC_OR_DIE(my_cgblock, (2, npwtot_k, nbb), ierr)

 ! TODO: This should be replaced by gatherv but premature optimization....
 do rank=0,nprocs-1
   bstart = bstart_rank(rank)
   bcount = bstart_rank(rank+1) - bstart_rank(rank)

   do band=bstart, bstart+bcount-1
     ib = band - bstart + 1
     cgbuf(:,:,ib) = zero
     do ig=1,npw_k
       igseq = gmpi2seq(ig)
       cgbuf(:,igseq,ib) = cg_k(:,ig,band)
     end do
   end do ! band

   call xmpi_sum_master(cgbuf,rank,comm_bandfft,ierr)
   if (me == rank) my_cgblock = cgbuf(:,:,:bcount)
 end do ! rank

 bstart = bstart_rank(me)
 bcount = bstart_rank(me+1) - bstart_rank(me)

 ABI_FREE(cgbuf)
 ABI_FREE(bstart_rank)

end subroutine cg2seqblocks
!!***

!----------------------------------------------------------------------

END MODULE m_iowf
!!***
