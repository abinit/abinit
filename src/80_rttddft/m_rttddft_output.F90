!!****m* ABINIT/m_rttddft_output
!! NAME
!!  m_rttddft_ouptut
!!
!! FUNCTION
!!  Manages most output of RT-TDDFT runs
!!
!! COPYRIGHT
!!  Copyright (C) 2021 ABINIT group (FB, MT)
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

 use m_crystal,       only: crystal_init, crystal_t
 use m_dtfil,         only: datafiles_type
 use m_dtset,         only: dataset_type
 use m_ebands,        only: ebands_init, ebands_free
 use m_epjdos,        only: dos_calcnwrite, partial_dos_fractions, &
                          & partial_dos_fractions_paw,epjdos_t,    &
                          & epjdos_new, prtfatbands, fatbands_ncwrite
 use m_errors,        only: msg_hndl, assert, netcdf_check
 use m_io_tools,      only: open_file
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_nctk,          only: nctk_open_create
 use m_paral_atom,    only: get_my_atmtab, free_my_atmtab
 use m_rttddft_types, only: tdks_type
 use m_specialmsg,    only: wrtout
 use m_xmpi,          only: xmpi_comm_rank, xmpi_comm_self
   
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
 !scalars
 character(len=500)   :: msg
 !arrays

! *************************************************************************

 if (istep == 0) then
   !Open energy file and write header
   if (open_file(dtfil%fnameabo_td_ener, msg, newunit=tdks%tdener_unit, status='unknown', form='formatted') /= 0) then
      write(msg,'(a,a)') 'Error while trying to open file ', dtfil%fnameabo_td_ener
      ABI_ERROR(msg)
   end if 
   write(msg,'(a)') "# RT-TDDFT -- Energy file.  All quantities are in atomic units"
   call wrtout(tdks%tdener_unit,msg)
   write(msg,'(a)') "#   step    time      E_total      E_kinetic    E_hartree      E_xc        E_ewald   &
                   & E_corepsp   E_localpsp   E_nonlocalpsp   E_paw       E_entropy     E_vdw"
   call wrtout(tdks%tdener_unit,msg)

   write(msg,'(a,f14.6,a)') 'Integrated density (ie. total nb of electrons) = ', &
                           & SUM(tdks%rhor(:,1))*tdks%ucvol/tdks%nfftf, ch10
   call wrtout(ab_out,msg)
   if (do_write_log) call wrtout(std_out,msg)
 end if

 if (istep /= 0) then
   !Write in output file
   write(msg,'(a,f14.6,a)') 'Total energy = ', tdks%etot,' Ha'
   call wrtout(ab_out,msg)
   if (do_write_log) call wrtout(std_out,msg)

   write(msg,'(a,f14.6,a)') 'Integrated density (ie. total nb of electrons) = ', &
                            & SUM(tdks%rhor(:,1))*tdks%ucvol/tdks%nfftf, ch10
   call wrtout(ab_out,msg)
   if (do_write_log) call wrtout(std_out,msg)

   !Write in energy file
   write(msg,'(i8,f10.5,11(f12.6,X))') istep, istep*tdks%dt, tdks%etot, tdks%energies%e_kinetic,                       &
                                     & tdks%energies%e_hartree, tdks%energies%e_xc, tdks%energies%e_ewald,             &
                                     & tdks%energies%e_corepsp, tdks%energies%e_localpsp, tdks%energies%e_nlpsp_vfock, & 
                                     & tdks%energies%e_paw, tdks%energies%e_entropy, tdks%energies%e_vdw_dftd
   call wrtout(tdks%tdener_unit,msg)
 end if

 !Update header, with evolving variables
 call tdks%hdr%update(tdks%bantot,tdks%etot,tdks%energies%e_fermie,tdks%energies%e_fermih, &
                    & tdks%hdr%residm,tdks%rprimd,tdks%occ,tdks%pawrhoij,                  &
                    & tdks%xred,dtset%amu_orig,comm_atom=mpi_enreg%comm_atom,              &
                    & mpi_atmtab=mpi_enreg%my_atmtab)

 if (MOD(istep,dtset%td_prtstr) == 0) then
   if (dtset%prtdos>=2.or.dtset%pawfatbnd>0) then
       CALL prt_dos(dtfil,dtset,istep,mpi_enreg,psps,tdks)
   end if
 end if
 !Close files at the end
 if (istep == tdks%ntime) then
   close(tdks%tdener_unit)
 end if

end subroutine rttddft_output

!!****f* m_rttddft_output/prt_dos
!!
!! NAME
!!  prt_dos
!!
!! FUNCTION
!!  Compute and outputs the electronic DOS
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
subroutine prt_dos(dtfil, dtset, istep, mpi_enreg, psps, tdks)

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
 integer,parameter     :: master=0
 integer               :: bantot
 integer               :: collect
 integer               :: iatom
 integer               :: spacecomm
 integer               :: my_comm_atom, my_natom
 integer               :: me
 integer               :: natom
#ifdef HAVE_NETCDF
 integer               :: ncid
#endif
 integer               :: timrev
 character(len=500)    :: msg
 character(len=fnlen ) :: fname
 character(len=24 )    :: step_nb
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

 remove_inv=.false.
 if (dtset%nspden==4 .and. dtset%usedmft==1) remove_inv=.true. ! MG: why this?

 !FB: Maybe this should be moved out of that subroutine if needed in other outputs than DOS
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

 !Generate DOS using the tetrahedron method or using Gaussians
 WRITE(step_nb,*) istep 
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
 
end subroutine prt_dos

end module m_rttddft_output
!!***
