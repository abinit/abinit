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

 use m_dtfil,         only: datafiles_type
 use m_dtset,         only: dataset_type
 use m_errors,        only: msg_hndl, assert
 use m_io_tools,      only: open_file
 use m_rttddft_types, only: tdks_type
 use m_specialmsg,    only: wrtout
   
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
subroutine rttddft_output(dtfil, dtset, istep, mpi_enreg, tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,              intent(in)    :: istep
 type(datafiles_type), intent(inout) :: dtfil
 type(dataset_type),   intent(inout) :: dtset
 type(MPI_type),       intent(inout) :: mpi_enreg
 type(tdks_type),      intent(inout) :: tdks
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
   write(msg,'(a)') "# RT-TDDFT -- Energy file  All quantities are in atomic units"
   call wrtout(tdks%tdener_unit,msg)
   write(msg,'(a)') "#   step    time      E_total    E_kinetic    E_hartree     E_xc     E_ewald   &
                   & E_localpsp    E_corepsp   E_entropy    E_vdw      E_paw"
   call wrtout(tdks%tdener_unit,msg)

   write(msg,'(a,f14.6,a)') 'Integrated density (ie. total nb of electrons) = ', &
                           & SUM(tdks%rhor(:,1))*tdks%ucvol/tdks%nfftf, ch10
   call wrtout(ab_out,msg)
   if (do_write_log) call wrtout(std_out,msg)
 end if

 if (istep /= 0) then
   !Update header, with evolving variables
   call tdks%hdr%update(tdks%bantot,tdks%etot,tdks%energies%e_fermie,tdks%energies%e_fermih, &
                       & tdks%hdr%residm,tdks%rprimd,tdks%occ,tdks%pawrhoij,                 &
                       & tdks%xred,dtset%amu_orig,comm_atom=mpi_enreg%comm_atom,             &
                       & mpi_atmtab=mpi_enreg%my_atmtab)

   !Write in output file
   write(msg,'(a,f14.6,a)') 'Total energy = ', tdks%etot,' Ha'
   call wrtout(ab_out,msg)
   if (do_write_log) call wrtout(std_out,msg)

   write(msg,'(a,f14.6,a)') 'Integrated density (ie. total nb of electrons) = ', &
                            & SUM(tdks%rhor(:,1))*tdks%ucvol/tdks%nfftf, ch10
   call wrtout(ab_out,msg)
   if (do_write_log) call wrtout(std_out,msg)

   !Write in energy file
   write(msg,'(i8,f10.5,10f12.6)') istep, istep*tdks%dt, tdks%etot, tdks%energies%e_kinetic,                   &
                                 & tdks%energies%e_hartree, tdks%energies%e_xc, tdks%energies%e_ewald,         &
                                 & tdks%energies%e_localpsp, tdks%energies%e_corepsp, tdks%energies%e_entropy, &
                                 & tdks%energies%e_vdw_dftd, tdks%energies%e_paw
   call wrtout(tdks%tdener_unit,msg)
 end if

 !Close files at the end
 if (istep == tdks%ntime) then
   close(tdks%tdener_unit)
 end if

end subroutine rttddft_output

end module m_rttddft_output
!!***
