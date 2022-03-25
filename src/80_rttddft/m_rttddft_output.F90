!!****m* ABINIT/m_rttddft_output
!! NAME
!!  m_rttddft_ouptut
!!
!! FUNCTION
!!  Manages most output of RT-TDDFT runs
!!
!! COPYRIGHT
!!  Copyright (C) 2021-2022 ABINIT group (FB, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!  m_rttddft_driver
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rttddft_output

 use defs_basis
 use defs_abitypes,   only: MPI_type
 use defs_datatypes,  only: pseudopotential_type, ebands_t

 use m_common,        only: prteigrs
 use m_crystal,       only: crystal_init, crystal_t
 use m_dtfil,         only: datafiles_type
 use m_dtset,         only: dataset_type
 use m_ebands,        only: ebands_init, ebands_free
 use m_epjdos,        only: dos_calcnwrite, partial_dos_fractions, &
                          & partial_dos_fractions_paw,epjdos_t,    &
                          & epjdos_new, prtfatbands, fatbands_ncwrite
 use m_errors,        only: msg_hndl, assert, netcdf_check
 use m_ioarr,         only: fftdatar_write
 use m_io_tools,      only: open_file
 use m_iowf,          only: outwf
 use m_mpinfo,        only: iwrite_fftdatar
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_nctk,          only: nctk_open_create
 use m_paral_atom,    only: get_my_atmtab, free_my_atmtab
 use m_rttddft_tdks,  only: tdks_type
 use m_specialmsg,    only: wrtout
 use m_xmpi,          only: xmpi_comm_rank
   
 implicit none

 private
!!***

 public :: rttddft_output
!!***

contains 

!!****f* m_rttddft_output/rttddft_output
!!
!! NAME
!!  rttddft_output
!!
!! FUNCTION
!!  Main output subroutine
!!
!! INPUTS
!!  dtfil <type datafiles_type> = infos about file names, file unit numbers
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_driver/rttddft
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_output(dtfil, dtset, istep, mpi_enreg, psps, tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                    intent(in)    :: istep
 type(datafiles_type),       intent(inout) :: dtfil
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 type(tdks_type),            intent(inout) :: tdks
 !arrays
 
 !Local variables-------------------------------
 integer :: me
 !scalars
 character            :: tmp
 character(len=500)   :: msg
 integer              :: stat
 !arrays

! *************************************************************************

 if (istep == tdks%first_step) then
    !Open energy file and write header if needed
   if (open_file(tdks%fname_tdener, msg, newunit=tdks%tdener_unit, status='unknown', form='formatted') /= 0) then
      write(msg,'(a,a)') 'Error while trying to open file ', tdks%fname_tdener
      ABI_ERROR(msg)
   end if 
   if (mpi_enreg%me == 0) then
      if (dtset%td_restart>0) then
         do
            read(tdks%tdener_unit,*,iostat=stat) tmp
            if (stat /= 0) exit
         end do
         backspace(tdks%tdener_unit)
      else
         write(msg,'(a)') "# RT-TDDFT -- Energy file. All quantities are in Hartree atomic units."
         call wrtout(tdks%tdener_unit,msg)
         write(msg,'(a)') "# step  time  E_total  E_kinetic  E_hartree  E_xc  E_ewald  &
                         & E_corepsp  E_localpsp  E_nonlocalpsp  E_paw  E_entropy  E_vdw"
         call wrtout(tdks%tdener_unit,msg)
      end if
   end if
 end if

!!FB: This is most probably not needed
!!Update header, with evolving variables
!call tdks%hdr%update(tdks%bantot,tdks%etot,tdks%energies%e_fermie,tdks%energies%e_fermih, &
!                   & tdks%hdr%residm,tdks%rprimd,tdks%occ,tdks%pawrhoij,                  &
!                   & tdks%xred,dtset%amu_orig,comm_atom=mpi_enreg%comm_atom,              &
!                   & mpi_atmtab=mpi_enreg%my_atmtab)

 !Write in output file
 write(msg,'(a,a,f14.6,a)') ch10, 'Total energy = ', tdks%etot,' Ha'
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 write(msg,'(a,f14.6,a)') 'Integrated density (ie. total nb of electrons) = ', &
                          & sum(tdks%rhor(:,1))*tdks%ucvol/tdks%nfftf, ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 !Writes in energy file
 write(msg,'(i8,f10.5,11(f14.8,1X))') istep-1, (istep-1)*tdks%dt, tdks%etot, tdks%energies%e_kinetic,                 &
                                    & tdks%energies%e_hartree, tdks%energies%e_xc, tdks%energies%e_ewald,             &
                                    & tdks%energies%e_corepsp, tdks%energies%e_localpsp, tdks%energies%e_nlpsp_vfock, &
                                    & tdks%energies%e_paw, tdks%energies%e_entropy, tdks%energies%e_vdw_dftd
 call wrtout(tdks%tdener_unit,msg)

!!Special case of last step
!!prints last WF for restart and close TD_ENER file
!if (istep == tdks%first_step+tdks%ntime-1) then
!   if (mod(istep,dtset%td_prtstr) == 0) then 
!      call prt_den(dtfil,dtset,istep,mpi_enreg,psps,tdks)
!      call prt_eig(dtfil,dtset,istep,mpi_enreg,psps,tdks)
!   end if
!   call prt_wfk(dtfil,dtset,istep,mpi_enreg,psps,tdks)
!   call prt_restart(dtfil,istep,mpi_enreg,tdks)
!   close(tdks%tdener_unit)
!else if (mod(istep,dtset%td_prtstr) == 0) then
!   call prt_den(dtfil,dtset,istep,mpi_enreg,psps,tdks)
!   if (dtset%prtwf > 0) then
!      call prt_wfk(dtfil,dtset,istep,mpi_enreg,psps,tdks)
!      call prt_restart(dtfil,istep,mpi_enreg,tdks)
!   end if
!end if

 if (mod(istep,dtset%td_prtstr) == 0) then
    call prt_den(dtfil,dtset,istep,mpi_enreg,psps,tdks)
    call prt_eig(dtfil,dtset,istep,mpi_enreg,tdks)
    if (dtset%prtwf > 0) then
       call prt_wfk(dtfil,dtset,istep,mpi_enreg,psps,tdks)
       call prt_restart(dtfil,istep,mpi_enreg,tdks)
    end if
 end if

 ! Special case of last step
 me = xmpi_comm_rank(mpi_enreg%comm_world)
 if (istep == tdks%first_step+tdks%ntime-1) then
    if (mod(istep,dtset%td_prtstr) /= 0 .or. dtset%prtwf <= 0) then 
      call prt_wfk(dtfil,dtset,istep,mpi_enreg,psps,tdks,force_write=.TRUE.)
      call prt_restart(dtfil,istep,mpi_enreg,tdks)
    end if
    close(tdks%tdener_unit)
 end if

end subroutine rttddft_output

!!****f* m_rttddft_output/prt_eig
!!
!! NAME
!!  prt_eig
!!
!! FUNCTION
!!  Outputs the eigenvalues
!!
!! INPUTS
!!  dtfil <type datafiles_type> = infos about file names, file unit numbers
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine prt_eig(dtfil, dtset, istep, mpi_enreg, tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                    intent(in)    :: istep
 type(datafiles_type),       intent(inout) :: dtfil
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(tdks_type),            intent(inout) :: tdks
 !arrays
 
 !Local variables-------------------------------
 !scalars
 integer,parameter     :: enunit=1, option=2
 integer               :: me
 integer               :: spacecomm
 real(dp)              :: vxcavg_dum
 character(len=fnlen)  :: fname
 character(len=24)     :: step_nb
 !arrays
 real(dp)              :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)

! *************************************************************************

 spacecomm = mpi_enreg%comm_cell
 me = xmpi_comm_rank(spacecomm)

 write(step_nb,*) istep
 fname = trim(dtfil%filnam_ds(4))//'_EIG_'//trim(adjustl(step_nb))
 resid = zero
 vxcavg_dum=zero

 if(me==0)then
   call prteigrs(tdks%eigen,enunit,tdks%energies%e_fermie,tdks%energies%e_fermih,  & 
               & fname,ab_out,dtset%iscf,dtset%kptns,dtset%kptopt,dtset%mband,     &
               & dtset%nband,dtset%nbdbuf,dtset%nkpt,0,dtset%nsppol,tdks%occ,      &
               & dtset%occopt, option,dtset%prteig,dtset%prtvol,resid,dtset%tolwfr, &
               & vxcavg_dum,dtset%wtk)
 end if

end subroutine prt_eig

!!****f* m_rttddft_output/prt_den
!!
!! NAME
!!  prt_den
!!
!! FUNCTION
!!  Computes and outputs the electronic density and/or DOS
!!
!! INPUTS
!!  dtfil <type datafiles_type> = infos about file names, file unit numbers
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine prt_den(dtfil, dtset, istep, mpi_enreg, psps, tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                    intent(in)    :: istep
 type(datafiles_type),       intent(inout) :: dtfil
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 type(tdks_type),            intent(inout) :: tdks
 !arrays
 
 !Local variables-------------------------------
 !scalars
 integer,parameter     :: master=0, cplex1=1
 integer               :: bantot
 integer               :: collect
 integer               :: iatom
 integer               :: spacecomm
 integer               :: my_comm_atom, my_natom
 integer               :: me
 integer               :: natom
!#ifdef HAVE_NETCDF
! integer               :: ncid
!#endif
 integer               :: timrev
 character(len=fnlen)  :: fname
 character(len=24)     :: step_nb
 logical               :: paral_atom
 logical               :: remove_inv
 logical               :: my_atmtab_allocated
 type(crystal_t)       :: crystal
 type(epjdos_t)        :: dos
 type(ebands_t)        :: ebands
 !arrays
 integer, pointer      :: my_atmtab(:)
 real(dp), allocatable :: doccde(:)

! *************************************************************************

 spacecomm = mpi_enreg%comm_cell
 me = xmpi_comm_rank(spacecomm)

 natom = dtset%natom
 my_natom = mpi_enreg%my_natom
 paral_atom=(my_natom/=natom)
 my_comm_atom = mpi_enreg%comm_atom
 nullify(my_atmtab)
 if (paral_atom) then
   call get_my_atmtab(mpi_enreg%comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)
 else
   ABI_MALLOC(my_atmtab, (natom))
   my_atmtab = (/ (iatom, iatom=1, natom) /)
   my_atmtab_allocated = .true.
 end if

 !FB: Maybe this should be moved out of that subroutine if needed in other outputs than densities?
 remove_inv=.false.
 timrev = 2; if (any(dtset%kptopt == [3, 4])) timrev= 1
 call crystal_init(dtset%amu_orig(:,1),crystal,dtset%spgroup,natom,dtset%npsp,psps%ntypat, &
   dtset%nsym,tdks%rprimd,dtset%typat,tdks%xred,dtset%ziontypat,dtset%znucl,timrev,&
   dtset%nspden==2.and.dtset%nsppol==1,remove_inv,tdks%hdr%title,&
   dtset%symrel,dtset%tnons,dtset%symafm)
 !Electron band energies.
 bantot= dtset%mband*dtset%nkpt*dtset%nsppol
 ABI_CALLOC(doccde, (bantot))
 call ebands_init(bantot,ebands,dtset%nelect,dtset%ne_qFD,dtset%nh_qFD,dtset%ivalence,         &
   doccde,tdks%eigen,dtset%istwfk,dtset%kptns,dtset%nband,dtset%nkpt,tdks%npwarr,dtset%nsppol, &
   dtset%nspinor,dtset%tphysel,dtset%tsmear,dtset%occopt,tdks%occ,dtset%wtk,&
   dtset%cellcharge(1),dtset%kptopt,dtset%kptrlatt_orig,dtset%nshiftk_orig,dtset%shiftk_orig, &
   dtset%kptrlatt,dtset%nshiftk,dtset%shiftk)
 ABI_FREE(doccde)

 write(step_nb,*) istep

 !** Outputs the density
 !Warnings :
 !- core charge is excluded from the charge density;
 !- the potential is the INPUT vtrial.
 if (iwrite_fftdatar(mpi_enreg) .and. dtset%prtden/=0) then
     fname = trim(dtfil%filnam_ds(4))//'_DEN_'//trim(adjustl(step_nb))
     call fftdatar_write("density",fname,dtset%iomode,tdks%hdr,crystal,tdks%pawfgr%ngfft, &
                       & cplex1,tdks%pawfgr%nfft,dtset%nspden,tdks%rhor,mpi_enreg,ebands=ebands)
 end if

 !** Generate DOS using the tetrahedron method or using Gaussians
 if (dtset%prtdos>=2.or.dtset%pawfatbnd>0) then
   dos = epjdos_new(dtset, psps, tdks%pawtab)

   if (dos%partial_dos_flag>=1 .or. dos%fatbands_flag==1)then
      ! Generate fractions for partial DOSs if needed partial_dos 1,2,3,4  give different decompositions
      collect = 1 !; if (psps%usepaw==1 .and. dos%partial_dos_flag /= 2) collect = 0
      if ((psps%usepaw==0.or.dtset%pawprtdos/=2) .and. dos%partial_dos_flag>=1) then
         call partial_dos_fractions(dos,crystal,dtset,tdks%eigen,tdks%occ,tdks%npwarr,tdks%kg,tdks%cg,tdks%mcg,collect,mpi_enreg)
      end if

      if (psps%usepaw==1 .and. dos%partial_dos_flag /= 2) then
         ! TODO: update partial_dos_fractions_paw for extra atoms - no PAW contribution normally, but check bounds and so on.
         call partial_dos_fractions_paw(dos,tdks%cprj,tdks%dimcprj,dtset,tdks%mcprj,dtset%mkmem,mpi_enreg,tdks%pawrad,tdks%pawtab)
      end if
   else
      dos%fractions(:,:,:,1)=one
   end if

   !Here, print out fatbands for the k-points given in file appended _FATBANDS
   if (me == master .and. dtset%pawfatbnd>0 .and. dos%fatbands_flag==1) then
      fname = trim(dtfil%filnam_ds(4))//'_FATBANDS_'//trim(adjustl(step_nb))
      call prtfatbands(dos,dtset,ebands,fname,dtset%pawfatbnd,tdks%pawtab)
   end if

   !Here, computation and output of DOS and partial DOS  _DOS
   if (dos%fatbands_flag == 0 .and. dos%prtdos /= 4) then
      fname = trim(dtfil%filnam_ds(4))//'_DOS_'//trim(adjustl(step_nb))
      call dos_calcnwrite(dos,dtset,crystal,ebands,fname,spacecomm)
   end if
 end if

!#ifdef HAVE_NETCDF
!   ! Write netcdf file with dos% results.
!   if (me == master) then
!     fname = trim(dtfil%filnam_ds(4))//'_FATBANDS_'//trim(adjustl(step_nb))//'.nc'
!     NCF_CHECK(nctk_open_create(ncid, fname, xmpi_comm_self))
!     call fatbands_ncwrite(dos, crystal, ebands, tdks%hdr, dtset, psps, tdks%pawtab, ncid)
!     NCF_CHECK(nf90_close(ncid))
!   end if
!#endif

 call dos%free()
 call crystal%free()
 call ebands_free(ebands)
 
end subroutine prt_den

!!****f* m_rttddft_output/prt_wfk
!!
!! NAME
!!  prt_wfk
!!
!! FUNCTION
!!  Outputs wavefunctions in WFK file
!!
!! INPUTS
!!  dtfil <type datafiles_type> = infos about file names, file unit numbers
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine prt_wfk(dtfil, dtset, istep, mpi_enreg, psps, tdks, force_write)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                    intent(in)    :: istep
 type(datafiles_type),       intent(inout) :: dtfil
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 type(tdks_type),            intent(inout) :: tdks
 logical,  optional,         intent(in)    :: force_write
 !arrays

 !Local variables-------------------------------
 !scalars
 integer,parameter     :: response=0
 character(len=fnlen)  :: fname
 character(len=24)     :: step_nb
 logical               :: lforce_write = .FALSE.
 !arrays

! *************************************************************************

 write(step_nb,*) istep
 !Use initial eigenvalues to ensure that we get the same occupation upon restart
 fname = trim(dtfil%filnam_ds(4))//'_WFK_'//trim(adjustl(step_nb))

 if (present(force_write)) then
    if (force_write) lforce_write = .TRUE.
 end if

 if (lforce_write) then
   call outwf(tdks%cg,dtset,psps,tdks%eigen0,fname,tdks%hdr,tdks%kg,dtset%kptns, &
             & dtset%mband,tdks%mcg,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,  &
             & dtset%nband,dtset%nkpt,tdks%npwarr,dtset%nsppol,tdks%occ,response, &
             & dtfil%unwff2,tdks%wvl%wfs,tdks%wvl%descr, force_write=.TRUE.)
 else
   call outwf(tdks%cg,dtset,psps,tdks%eigen0,fname,tdks%hdr,tdks%kg,dtset%kptns, &
             & dtset%mband,tdks%mcg,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,  &
             & dtset%nband,dtset%nkpt,tdks%npwarr,dtset%nsppol,tdks%occ,response, &
             & dtfil%unwff2,tdks%wvl%wfs,tdks%wvl%descr)
 end if

end subroutine prt_wfk

!!****f* m_rttddft_output/prt_restart
!!
!! NAME
!!  prt_restart
!!
!! FUNCTION
!!  Print restart file
!!
!! INPUTS
!!  dtfil <type datafiles_type> = infos about file names, file unit numbers
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  tdks <type(tdks_type)> = the tdks object to initialize
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine prt_restart(dtfil, istep, mpi_enreg, tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                    intent(in)    :: istep
 type(datafiles_type),       intent(inout) :: dtfil
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(tdks_type),            intent(inout) :: tdks
 !arrays
 
 !Local variables-------------------------------
 !scalars
 character(len=500)    :: msg
 character(len=fnlen)  :: fname
 character(len=24 )    :: step_nb
 !arrays

! *************************************************************************

 if (mpi_enreg%me == 0) then 
   write(step_nb,*) istep
   rewind(tdks%tdrestart_unit)
   write(msg,'(a)') step_nb
   call wrtout(tdks%tdrestart_unit,msg)
   write(msg,'(a)') tdks%fname_tdener
   call wrtout(tdks%tdrestart_unit,msg)
   fname = trim(dtfil%filnam_ds(4))//'_WFK_'//trim(adjustl(step_nb))
   write(msg,'(a)') fname
   call wrtout(tdks%tdrestart_unit,msg)
 end if

end subroutine prt_restart

end module m_rttddft_output
!!***
