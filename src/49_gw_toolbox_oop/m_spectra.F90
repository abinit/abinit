!!****m* ABINIT/m_spectra
!! NAME
!! m_spectra
!!
!! FUNCTION
!! This module contains the definition of the specta_type data type
!! used to store results related to optical spectra with or without
!! nonlocal field effects as well as the electron energy loss function
!! for a given q. These quantities are obtained from the dielectric
!! matrix as calculated in the GW part of ABINIT (screening.F90)
!!
!! COPYRIGHT
!! Copyright (C) 2008-2019 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_spectra

 use defs_basis
 use m_errors
 use m_abicore

 use m_io_tools,      only : open_file
 use m_fstrings,      only : strcat

 implicit none

 private
!!***

!!****t* m_spectra/spectra_t
!! NAME
!! spectra_t
!!
!! FUNCTION
!!  Object used to store optical spectra with or without non-local field effects.
!!
!! SOURCE

 type,public :: spectra_t

 !scalars
  integer :: nomega
  ! number of frequencies

  integer :: nqpts
  ! number of q-points

 !arrays
  real(dp),allocatable :: omega(:)
  ! omega(nomega)
  ! Real frequency mesh for optical spectra.

  real(dp),allocatable :: qpts(:,:)
  ! qpts(3,nqpoints)
  ! The list of q-points used for the spectra

  real(dp),allocatable :: eelf(:,:)
  ! eelf(nomega,nqpoints)
  ! contains the Electron Energy Loss Function i.e. -\Im{ e^{-1}_{G1=0,G2=0}(q-->0,nomega)}

  complex(dpc),allocatable :: emacro_lf(:,:)
  ! emacro_lf(nomega,nqpoints)
  ! contains 1/e^{-1}_{G1=0,G2=0}(q-->0,nomega) (with Local field effects)

  complex(dpc),allocatable :: emacro_nlf(:,:)
  ! emacro_nlf(nomega,nqpoints)
  ! contains e_{G1=0,G2=0}(q-->0,nomega) (without Local field effects)

 contains

   procedure :: free => spectra_free
   ! Free memory.

   procedure :: write => spectra_write
   ! Write results on file.

   procedure :: repr => spectra_repr
   ! Return info on Macroscopic diel. constant in form of a string.

 end type spectra_t
!!***

 public :: spectra_init       ! Creation method.

 integer,public,parameter :: W_EM_LF  = 1
 integer,public,parameter :: W_EM_NLF = 2
 integer,public,parameter :: W_EELF   = 4

CONTAINS  !========================================================================================
!!***

!!****f* m_spectra/spectra_init
!! NAME
!!  spectra_init
!!
!! FUNCTION
!!  Initialize the object.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine spectra_init(Spectra,nomega,omega,nqpts,qpts)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,nqpts
!arrays
 real(dp),intent(in) :: omega(nomega),qpts(3,nqpts)
 type(spectra_t),intent(out) :: Spectra

! *********************************************************************

 Spectra%nomega = nomega
 Spectra%nqpts  = nqpts

 ABI_MALLOC(Spectra%qpts,(3,nqpts))
 Spectra%qpts = qpts

 ABI_MALLOC(Spectra%omega,(nomega))
 Spectra%omega = omega

 ABI_CALLOC(Spectra%emacro_lf,(nomega,nqpts))
 ABI_CALLOC(Spectra%emacro_nlf,(nomega,nqpts))
 ABI_CALLOC(Spectra%eelf,(nomega,nqpts))

end subroutine spectra_init
!!***

!----------------------------------------------------------------------

!!****f* m_spectra/spectra_free
!! NAME
!!  spectra_free
!!
!! FUNCTION
!!  Deallocate all associated pointers defined in the structure.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_screening,screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine spectra_free(Spectra)

!Arguments ------------------------------------
 class(spectra_t),intent(inout) :: spectra

! *********************************************************************

 ABI_SFREE(Spectra%omega)
 ABI_SFREE(Spectra%qpts)
 ABI_SFREE(Spectra%emacro_lf)
 ABI_SFREE(Spectra%emacro_nlf)
 ABI_SFREE(Spectra%eelf)

end subroutine spectra_free
!!***

!----------------------------------------------------------------------

!!****f* m_spectra/spectra_write
!! NAME
!!  spectra_write
!!
!! FUNCTION
!!  Write the optical spectra stored in the object on an external formatted file.
!!
!! INPUTS
!!  Spectra=The Object containing the spectra
!!  write_bits=Positive integer defining the quantities to be written (bit representation is used)
!!  fname=Name of the file to be written.
!!
!! OUTPUT
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine spectra_write(Spectra,write_bits,fname)

!Arguments ------------------------------------
!scalars
 class(spectra_t),intent(in) :: Spectra
 integer,intent(in) :: write_bits
 character(len=*),intent(in) :: fname

!Local variables-------------------------------
!scalars
 integer :: unt,io,iqpt
 real(dp) :: mino,maxo
 character(len=100) :: fmt
 character(len=500) :: msg

! *********************************************************************

 if (write_bits<0) RETURN

 mino = MINVAL(Spectra%omega)
 maxo = MAXVAL(Spectra%omega)

 if (open_file(fname,msg,newunit=unt) /= 0) then
   MSG_ERROR(msg)
 end if

 !write(unt,'(a,i5,2(a,f9.1),a)')'# nomega : ',Spectra%nomega,' from ',mino*Ha_eV,' up to ',maxo*Ha_eV,' [eV] '

 if ( IAND(write_bits,W_EM_NLF ) == W_EM_NLF ) then
   write(unt,'(a)')'#'
   write(unt,'(a)')'# Macroscopic Dielectric Function without local fields'
   call dump_qlist()
   write(unt,'(a)')'# Omega [eV]    Re epsilon_M       IM eps_M '
   write(fmt,'(a,i3,a)') '(1x,f7.3,7x,',Spectra%nqpts,'(2(e12.4),2x))'
   do io=1,Spectra%nomega
     write(unt,fmt) Spectra%omega(io)*Ha_eV, ( Spectra%emacro_nlf(io,iqpt), iqpt=1,Spectra%nqpts )
   end do
 end if
 !
 if ( IAND(write_bits,W_EM_LF ) == W_EM_LF ) then
   write(unt,'(a)')'#'
   write(unt,'(a)')'# Macroscopic Dielectric Function with local fields included'
   call dump_qlist()
   write(unt,'(a)')'# Omega [eV]    Re epsilon_M       Im eps_M '
   write(fmt,'(a,i3,a)') '(1x,f7.3,7x,',Spectra%nqpts,'(2(e12.4),2x))'
   do io=1,Spectra%nomega
     write(unt,fmt) Spectra%omega(io)*Ha_eV, ( Spectra%emacro_lf(io,iqpt), iqpt=1,Spectra%nqpts )
   end do
 end if
 !
 if ( IAND(write_bits,W_EELF) == W_EELF) then
   write(unt,'(a)')'#'
   write(unt,'(a)')'# Electron Energy Loss Function -Im(1/epsilon_M)'
   call dump_qlist()
   write(unt,'(a)')'# Omega [eV]    -Im(1/epsilon_M)'
   write(fmt,'(a,i3,a)') '(1x,f7.3,7x,',Spectra%nqpts,'(e12.4,2x))'
   do io=1,Spectra%nomega
     write(unt,fmt) Spectra%omega(io)*Ha_eV, ( Spectra%eelf(io,iqpt), iqpt=1,Spectra%nqpts ) ! -AIMAG(chi0(1,1,io))
   end do
 end if

 close(unt)

CONTAINS
!!***

!!****f* m_spectra/dump_Qlist
!! NAME
!!  dump_Qlist
!!
!! FUNCTION
!!  Helper function used to write the list of q-points used for the long-wavelength limit.
!!
!! INPUTS
!!  Spectra=The Object containing the spectra
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_spectra
!!
!! CHILDREN
!!
!! SOURCE

 subroutine dump_Qlist()

  integer :: iqpt
  write(unt,'(a,i3)')'# Q-point list, No. ',Spectra%nqpts
  do iqpt=1,Spectra%nqpts
    write(unt,'(a,i3,a,3f9.6,a)')'# ',iqpt,')  [',Spectra%qpts(:,iqpt),'] r.l.u. '
  end do
 end subroutine dump_Qlist

end subroutine spectra_write
!!***

!----------------------------------------------------------------------

!!****f* m_spectra/spectra_repr
!! NAME
!!  spectra_repr
!!
!! FUNCTION
!!  Returns a string reporting info on the calculated dielectric constant.
!!
!! INPUTS
!!  Spectra=The Object containing the spectra
!!
!! OUTPUT
!!
!! PARENTS
!!      m_screening,screening
!!
!! CHILDREN
!!
!! SOURCE

subroutine spectra_repr(Spectra,str)

!Arguments ------------------------------------
!scalars
 class(spectra_t),intent(in) :: Spectra
 character(len=*),intent(out) :: str

!Local variables-------------------------------
!scalars
 integer :: iqpt
 real(dp) :: epsilon0,epsilon0_nlf
 character(len=500) :: msg

! *********************************************************************

 !istatic = -1
 !do io = Spectra%nomega
 ! if (ABS(REAL(Spectra%omega(io)))<1.e-3.and.ABS(AIMAG(Spectra%omega(io)))<1.e-3) then
 !  istatic = io
 !  EXIT
 ! end if
 !end do

 str = ""
 do iqpt=1,Spectra%nqpts
   epsilon0    = REAL(Spectra%emacro_lf (1,iqpt))
   epsilon0_nlf= REAL(Spectra%emacro_nlf(1,iqpt))
   write(msg,'(a,3f9.6,a)')' For q-point: ',Spectra%qpts(:,iqpt),ch10
   str = strcat(str,msg)
   write(msg,'(1x,a,f8.4,a)')' dielectric constant = ',epsilon0,ch10
   str = strcat(str,msg)
   write(msg,'(1x,a,f8.4,a)')' dielectric constant without local fields = ',epsilon0_nlf,ch10
   str = strcat(str,msg)
 end do

end subroutine spectra_repr

!----------------------------------------------------------------------

END MODULE m_spectra
!!***
