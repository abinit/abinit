!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_eprenorms
!! NAME
!! m_eprenorms
!!
!! FUNCTION
!! This module contains datatypes to compute the renormalization of electronic states due to
!! eph coupling and temperature effects
!!
!! NOTES
!! This code is still under development and the API will change in the next versions.
!! Contact gmatteo
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (YG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_eprenorms

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_errors
 use m_xmpi
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_nctk

 use m_crystal,  only : crystal_t
 use m_kpts,     only : listkk

 implicit none

 private
!!***

!!****t* m_eprenorms/eprenorms_t
!! NAME
!! eprenorms_t
!!
!! FUNCTION
!! Datatype gathering data for electron-phonon renormalization of the band structure
!!
!! SOURCE

 type,public :: eprenorms_t

  !scalars
  integer :: nkpt
  ! Number of kpoints

  integer :: nsppol
  ! Number of spin channels

  integer :: mband
  ! Maximum number of bands

  integer :: ntemp
  ! Number of temperatures

  !arrays
  real(dp), allocatable :: kpts(:,:)
  ! kpt(3,nkpt)
  ! Kpoints

  real(dp), allocatable :: temps(:)
  ! temps(ntemp)
  ! Temperatures

  real(dp), allocatable :: eigens(:,:,:)
  ! eigens(mband,nkpt,nsppol)
  ! Kohn-Sham eigenvalues

  real(dp), allocatable :: occs(:,:,:)
  ! occ(mband,nkpt,nsppol)
  ! Occupation numbers

  real(dp), allocatable :: renorms(:,:,:,:,:)
  ! renorms(2,mband,nkpt,nsppol,ntemp)
  ! Renormalization of the eigenvalues for each temperature

  real(dp), allocatable :: linewidth(:,:,:,:,:)
  ! linewidth(2,mband,nkpt,nsppol,ntemp)
  ! Electron-phonon induced linewidth of the eigens

 end type eprenorms_t

 public :: eprenorms_init
 public :: eprenorms_free
 public :: eprenorms_from_epnc
 public :: eprenorms_bcast

 public :: renorm_bst
!!***

CONTAINS  !============================================================================
!!***

!!****f* m_eprenorms/eprenorms_init
!! NAME
!! eprenorms_init
!!
!! FUNCTION
!!  Initializes an eprenorms_t datatype
!!
!! INPUTS
!!
!! OUTPUT
!!  Epren<eprenorms_t>=Datatype gathering electron-phonon renormalizations
!!
!! PARENTS
!!      m_eprenorms
!!
!! CHILDREN
!!      listkk
!!
!! SOURCE

subroutine eprenorms_init(Epren,nkpt,nsppol,mband,ntemp)

!Arugments -----------------------------------
!scalars
 integer,intent(in) :: nkpt, nsppol, mband, ntemp
 type(eprenorms_t) :: Epren
!arrays

!*************************************************************************

 DBG_ENTER("COLL")

 Epren%nkpt = nkpt
 Epren%nsppol = nsppol
 Epren%mband = mband
 Epren%ntemp = ntemp

 ABI_MALLOC(Epren%temps,(Epren%ntemp))
 ABI_MALLOC(Epren%kpts,(3,Epren%nkpt))
 ABI_MALLOC(Epren%eigens,(Epren%mband,Epren%nkpt,Epren%nsppol))
 ABI_MALLOC(Epren%occs,(Epren%mband,Epren%nkpt,Epren%nsppol))
 ABI_MALLOC(Epren%renorms,(2,Epren%mband,Epren%nkpt,Epren%nsppol,Epren%ntemp))
 ABI_MALLOC(Epren%linewidth,(2,Epren%mband,Epren%nkpt,Epren%nsppol,Epren%ntemp))

 DBG_EXIT("COLL")

end subroutine eprenorms_init
!!***

!---------------------------------------------------------------------

!!****f* m_eprenorms/eprenorms_free
!! NAME
!! eprenorms_free
!!
!! FUNCTION
!! Deallocate all memory associated with eprenorms
!!
!! INPUTS
!! Epren<eprenorms_t>=The datatype to be freed
!!
!! PARENTS
!!      bethe_salpeter,optic
!!
!! CHILDREN
!!      listkk
!!
!! SOURCE

subroutine eprenorms_free(Epren)

!Arguments -----------------------------------
!scalars
 type(eprenorms_t),intent(inout) :: Epren

!*********************************************************************

 ABI_SFREE(Epren%temps)
 ABI_SFREE(Epren%kpts)
 ABI_SFREE(Epren%eigens)
 ABI_SFREE(Epren%occs)
 ABI_SFREE(Epren%renorms)
 ABI_SFREE(Epren%linewidth)

end subroutine eprenorms_free
!!***

!---------------------------------------------------------------------

!!****f* m_eprenorms/eprenorms_from_epnc
!! NAME
!! eprenorms_from_epnc
!!
!! FUNCTION
!! Allocates and initializes the datatype from a _EP.nc file
!!
!! INPUTS
!! filename = name of the file to be read
!!
!! SIDE EFFECTS
!! Epren<eprenorms_t> = fields are initialized and filled with data from filename
!!
!! PARENTS
!!      optic,setup_bse
!!
!! CHILDREN
!!      listkk
!!
!! SOURCE

subroutine eprenorms_from_epnc(Epren,filename)

!Arguments -----------------------------------
!scalars
 character(len=fnlen),intent(in) :: filename
 type(eprenorms_t),intent(inout) :: Epren

!Local variables------------------------------
 integer :: nkpt, mband, nsppol, ntemp
#ifdef HAVE_NETCDF
 integer :: ncid
#endif

! ************************************************************************

#ifdef HAVE_NETCDF
 NCF_CHECK(nctk_open_read(ncid, filename, xmpi_comm_self))
 NCF_CHECK(nctk_set_datamode(ncid))

 NCF_CHECK(nctk_get_dim(ncid, "number_of_kpoints", nkpt))
 NCF_CHECK(nctk_get_dim(ncid, "number_of_spins",nsppol))
 NCF_CHECK(nctk_get_dim(ncid, "max_number_of_states",mband))
 NCF_CHECK(nctk_get_dim(ncid, "number_of_temperature",ntemp))

 call eprenorms_init(Epren, nkpt, nsppol, mband, ntemp)

 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid,"reduced_coordinates_of_kpoints"), Epren%kpts))
 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid,"temperature"), Epren%temps))
 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid,"eigenvalues"), Epren%eigens))
 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid,"occupations"), Epren%occs))
 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid,"zero_point_motion"), Epren%renorms))
 ! TODO: This should be changed. What is stored is a linewidth, not a lifetime,
 ! we postone the change so as to not break compatibility
 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid,"lifetime"), Epren%linewidth))

#endif

end subroutine eprenorms_from_epnc
!!***

!---------------------------------------------------------------------

!!****f* m_eprenorms/eprenorms_bcast
!! NAME
!! eprenorms_bcast
!!
!! FUNCTION
!!  MPI broadcast all the content of eprenorms_t
!!
!! INPUTS
!!  master = Rank of master
!!  comm = MPI communicator
!!
!! SIDE EFFECTS
!!  Epren<eprenorms_t> = Data broadcasted on every node from master
!!
!! PARENTS
!!      optic,setup_bse
!!
!! CHILDREN
!!      listkk
!!
!! SOURCE

subroutine eprenorms_bcast(Epren,master,comm)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: master, comm
 type(eprenorms_t),intent(inout) :: Epren

!Local variables------------------------------
!scalars
 integer :: ierr

! ************************************************************************

 if (xmpi_comm_size(comm) == 1) return

 call xmpi_bcast(Epren%nkpt, master, comm, ierr)
 call xmpi_bcast(Epren%nsppol, master, comm, ierr)
 call xmpi_bcast(Epren%mband, master, comm, ierr)
 call xmpi_bcast(Epren%ntemp, master, comm, ierr)

 if (xmpi_comm_rank(comm) /= master) then
  call eprenorms_init(Epren, Epren%nkpt, Epren%nsppol, Epren%mband, Epren%ntemp)
 end if

 call xmpi_bcast(Epren%kpts, master, comm, ierr)
 call xmpi_bcast(Epren%temps, master, comm, ierr)
 call xmpi_bcast(Epren%eigens, master, comm, ierr)
 call xmpi_bcast(Epren%occs, master, comm, ierr)
 call xmpi_bcast(Epren%renorms, master, comm, ierr)
 call xmpi_bcast(Epren%linewidth, master, comm, ierr)

end subroutine eprenorms_bcast
!!***

!---------------------------------------------------------------------

!!****f* m_eprenorms/renorm_bst
!! NAME
!! renorm_bst
!!
!! FUNCTION
!!  Renormalize the band structure Bst from data contained Epren
!!
!! INPUTS
!!  Epren<eprenorms_t> = datatype containing the elphon renormalization
!!  itemp = index of the temperature you want to use
!!  do_lifetime = .true. if we want to use imaginary eigenvalues (lifetime field)
!!
!! SIDE EFFECTS
!!  Bst<bands_t> : eigens are changed according to epren
!!                 linewidth is allocated and filled with data if do_linewidth
!!
!! PARENTS
!!      m_exc_spectra,m_haydock,optic
!!
!! CHILDREN
!!      listkk
!!
!! SOURCE

subroutine renorm_bst(Epren,Bst,Cryst,itemp,do_lifetime,do_check)

!Arguments -----------------------------------
!scalars
 integer :: itemp
 logical,intent(in) :: do_lifetime
 logical,optional,intent(in) :: do_check
 type(eprenorms_t),intent(in) :: Epren
 type(ebands_t),intent(inout) :: Bst
 type(crystal_t),intent(in) :: Cryst

!Local variables------------------------------
!scalars
 integer :: isppol,ikpt,comm
 integer :: nband1, nband_tmp
 integer :: timrev, sppoldbl
 integer :: ik_eph
 real(dp) :: dksqmax
 logical :: check
!arrays
 integer,allocatable :: bs2eph(:,:)

! ************************************************************************

 ABI_CHECK(Bst%nsppol == Epren%nsppol, "Nsppol should be the same")

 comm = xmpi_comm_self

 if(do_lifetime) then
   ABI_MALLOC(Bst%linewidth,(1,Bst%mband,Bst%nkpt,Bst%nsppol))
 end if

 check = .TRUE.
 if(present(do_check)) then
   check = do_check
 end if

 sppoldbl = 1 !; if (any(Cryst%symafm == -1) .and. Epren%nsppol == 1) nsppoldbl=2
 ABI_MALLOC(bs2eph, (BSt%nkpt*sppoldbl, 6))
 timrev = 1
 call listkk(dksqmax, Cryst%gmet, bs2eph, Epren%kpts, BSt%kptns, Epren%nkpt, Bst%nkpt, Cryst%nsym, &
&   sppoldbl, Cryst%symafm, Cryst%symrel, timrev, comm, use_symrec=.False.)

 do isppol=1,Bst%nsppol
   do ikpt=1,Bst%nkpt
     nband1 = Bst%nband(ikpt+(isppol-1)*Bst%nkpt)
     nband_tmp=MIN(nband1,Epren%mband)

     ik_eph = bs2eph(ikpt,1)

     !FIXME change check
     if (check) then
       if (ANY(ABS(Bst%eig(1:MIN(10,nband_tmp),ikpt,isppol) - Epren%eigens(1:MIN(10,nband_tmp),ik_eph,isppol)) > tol3)) then
         !write(stdout,*) "eig : ",BSt%eig(1:MIN(10,nband_tmp),ikpt,isppol)
         !write(stdout,*) "eigens : ",Epren%eigens(1:MIN(10,nband_tmp),ikpt,isppol)
         MSG_ERROR("Error in eigenvalues, check the _EP.nc file with respect to your input file !")
       end if
       if (ANY(ABS(Bst%occ(1:MIN(10,nband_tmp),ikpt,isppol) - Epren%occs(1:MIN(10,nband_tmp),ik_eph,isppol)) > tol3)) then
         MSG_ERROR("Error in occupations, check the _EP.nc file with respect to your input file !")
       end if
     end if

     ! Upgrade energies
     Bst%eig(1:nband_tmp,ikpt,isppol) = BSt%eig(1:nband_tmp,ikpt,isppol) + Epren%renorms(1,1:nband_tmp,ik_eph,isppol,itemp)

     if (do_lifetime) then
       Bst%linewidth(1,1:nband_tmp,ikpt,isppol) = Epren%linewidth(1,1:nband_tmp,ik_eph,isppol,itemp)
     end if
   end do
 end do

 ABI_FREE(bs2eph)

end subroutine renorm_bst
!!***

!---------------------------------------------------------------------

end module m_eprenorms
!!***
