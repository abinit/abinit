!!****m* ABINIT/m_efmas
!! NAME
!! m_efmas
!!
!! FUNCTION
!! This module contains datatypes for efmas functionalities.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (JLJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_efmas

 use defs_basis
 use m_errors
 use m_abicore
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_efmas_defs
 use m_nctk
 use m_cgtools
 use m_dtset

 use defs_abitypes,         only : MPI_type
 use m_gaussian_quadrature, only : cgqf
 use m_io_tools, only : get_unit

 implicit none

 private

!public procedures.
 public :: efmasval_free
 public :: efmasval_free_array
 public :: efmasdeg_free
 public :: efmasdeg_free_array
 public :: efmas_ncread
 public :: check_degeneracies
 public :: print_tr_efmas
 public :: print_efmas
 public :: efmas_main
 public :: efmas_analysis

!private procedures.
 private :: MATMUL_ ! Workaround to make tests pass on ubu/buda slaves
 interface MATMUL_
  module procedure MATMUL_DP
  module procedure MATMUL_DPC
 end interface MATMUL_

!!***

CONTAINS
!===========================================================

!!****f* m_efmas/efmasval_free
!! NAME
!! efmasval_free
!!
!! FUNCTION
!! This routine deallocates an efmasval_type.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_efmas
!!
!! CHILDREN
!!
!! SOURCE

 subroutine efmasval_free(efmasval)

   !Arguments ------------------------------------
   type(efmasval_type),intent(inout) :: efmasval

   ! *********************************************************************

   if(allocated(efmasval%ch2c)) then
     ABI_FREE(efmasval%ch2c)
   end if
   if(allocated(efmasval%eig2_diag)) then
     ABI_FREE(efmasval%eig2_diag)
   end if

 end subroutine efmasval_free
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/efmasval_free_array
!! NAME
!! efmasval_free_array
!!
!! FUNCTION
!! This routine deallocates an efmasval_type or, optionally, an array of efmasval_type.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!
!! SOURCE

 subroutine efmasval_free_array(efmasval)

   !Arguments ------------------------------------
   type(efmasval_type),allocatable,intent(inout) :: efmasval(:,:)

   !!!Local variables-------------------------------
   integer :: i,j,n(2)

   ! *********************************************************************

   !XG20180810: please do not remove. Otherwise, I get an error on my Mac.
   write(std_out,*)' efmasval_free_array : enter '

   if(allocated(efmasval)) then
     n=shape(efmasval)
     do i=1,n(1)
       do j=1,n(2)
         call efmasval_free(efmasval(i,j))
       end do
     end do
     ABI_DATATYPE_DEALLOCATE(efmasval)
   end if

 end subroutine efmasval_free_array
!!***

!----------------------------------------------------------------------

!!****f* m_efmas/efmasdeg_free
!! NAME
!! efmasdeg_free
!!
!! FUNCTION
!! This routine deallocates an efmasdeg_type.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_efmas
!!
!! CHILDREN
!!
!! SOURCE

 subroutine efmasdeg_free(efmasdeg)

   !Arguments ------------------------------------
   type(efmasdeg_type),intent(inout) :: efmasdeg

   ! *********************************************************************
   if(allocated(efmasdeg%degs_bounds)) then
     ABI_FREE(efmasdeg%degs_bounds)
   end if
   if(allocated(efmasdeg%ideg)) then
     ABI_FREE(efmasdeg%ideg)
   end if
 end subroutine efmasdeg_free
!!***

!----------------------------------------------------------------------

!!****f* m_efmas/efmasdeg_free_array
!! NAME
!! efmasdeg_free_array
!!
!! FUNCTION
!! This routine deallocates an efmasdeg_type or, optionally, an array of efmasdeg_type.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!
!! SOURCE

 subroutine efmasdeg_free_array(efmasdeg)

   !Arguments ------------------------------------
   type(efmasdeg_type),allocatable,intent(inout) :: efmasdeg(:)

   !!!Local variables-------------------------------
   integer :: i,n

   ! *********************************************************************
   if(allocated(efmasdeg)) then
     n=size(efmasdeg)
     do i=1,n
       call efmasdeg_free(efmasdeg(i))
     end do
     ABI_DATATYPE_DEALLOCATE(efmasdeg)
   end if
 end subroutine efmasdeg_free_array
!!***

!----------------------------------------------------------------------

!!****f* m_efmas/check_degeneracies
!! NAME
!! check_degeneracies
!!
!! FUNCTION
!! This routine check for 0th order band degeneracies at given k-point.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      d2frnl
!!
!! CHILDREN
!!
!! SOURCE

 subroutine check_degeneracies(efmasdeg,bands,nband,eigen,deg_tol)

   !Arguments ------------------------------------
   type(efmasdeg_type),intent(out) :: efmasdeg
   integer,intent(in) :: bands(2),nband
   real(dp),intent(in) :: eigen(nband)
   real(dp),intent(in),optional :: deg_tol

   !!!Local variables-------------------------------
   integer :: deg_dim,iband, ideg
   integer, allocatable :: degs_bounds(:,:)
   real(dp) :: tol
   real(dp) :: eigen_tmp(nband)
   logical :: treated

   ! *********************************************************************

   tol=tol5; if(present(deg_tol)) tol=deg_tol

   !!! Determine sets of degenerate states in eigen0, i.e., at 0th order.
   efmasdeg%ndegs=1
   efmasdeg%nband=nband
   ABI_MALLOC(degs_bounds,(2,nband))
   ABI_MALLOC(efmasdeg%ideg, (nband))
   degs_bounds=0; degs_bounds(1,1)=1
   efmasdeg%ideg=0; efmasdeg%ideg(1)=1

   eigen_tmp(:) = eigen(:)

   do iband=2,nband
     if (ABS(eigen_tmp(iband)-eigen_tmp(iband-1))>tol) then
       degs_bounds(2,efmasdeg%ndegs) = iband-1
       efmasdeg%ndegs=efmasdeg%ndegs+1
       degs_bounds(1,efmasdeg%ndegs) = iband
     end if
     efmasdeg%ideg(iband) = efmasdeg%ndegs
   end do
   degs_bounds(2,efmasdeg%ndegs)=nband
   ABI_MALLOC(efmasdeg%degs_bounds,(2,efmasdeg%ndegs))
   efmasdeg%degs_bounds(1:2,1:efmasdeg%ndegs) = degs_bounds(1:2,1:efmasdeg%ndegs)
   ABI_FREE(degs_bounds)

   !!! Determine if treated bands are part of a degeneracy at 0th order.
   efmasdeg%deg_range=0
   deg_dim=0
   treated=.false.
   write(std_out,'(a,i6)') 'Number of sets of bands for this k-point:',efmasdeg%ndegs
   write(std_out,'(a)') 'Set index; range of bands included in the set; is the set degenerate?(T/F); &
&                        is the set treated by EFMAS?(T/F):'
   do ideg=1,efmasdeg%ndegs
     deg_dim = efmasdeg%degs_bounds(2,ideg) - efmasdeg%degs_bounds(1,ideg) + 1
     !If there is some level in the set that is inside the interval defined by bands(1:2), treat such set
     !The band range might be larger than the nband interval: it includes it, and also include degenerate states
     if(efmasdeg%degs_bounds(1,ideg)<=bands(2) .and. efmasdeg%degs_bounds(2,ideg)>=bands(1)) then
       treated = .true.
       if(efmasdeg%degs_bounds(1,ideg)<=bands(1)) then
         efmasdeg%deg_range(1) = ideg
       end if
       if(efmasdeg%degs_bounds(2,ideg)>=bands(2)) then
         efmasdeg%deg_range(2) = ideg
       end if
     end if
     write(std_out,'(2i6,a,i6,2l4)') ideg, efmasdeg%degs_bounds(1,ideg), ' -', efmasdeg%degs_bounds(2,ideg), &
&                                    (deg_dim>1), treated
   end do

!   write(std_out,*)'ndegs=',          efmasdeg%ndegs
!   write(std_out,*)'degs_bounds=',    efmasdeg%degs_bounds
!   write(std_out,*)'ideg=',           efmasdeg%ideg
!   write(std_out,*)'deg_range=',      efmasdeg%deg_range

  !!This first attempt WORKS, but only if the symmetries are enabled, see line 1578 of dfpt_looppert.F90.
  !use m_crystal,          only : crystal_t, crystal_init, crystal_free, crystal_print
  !use m_esymm
  !integer :: timrev
  !character(len=132),allocatable :: title(:)
  !type(crystal_t) :: Cryst
  !type(esymm_t) :: Bsym

  !timrev = 1
  !if(dtset%istwfk(1)/=1) timrev=2
  !ABI_ALLOCATE(title,(dtset%ntypat))
  !title(:) = "Bloup"
  !call crystal_init(Cryst,dtset%spgroup,dtset%natom,dtset%npsp,dtset%ntypat,dtset%nsym,dtset%rprimd_orig(:,:,1),&
  !&                 dtset%typat,dtset%xred_orig(:,:,1),dtset%ziontypat,dtset%znucl,timrev,.false.,.false.,title,&
  !&                 dtset%symrel,dtset%tnons,dtset%symafm)
  !call crystal_print(Cryst)
  !ABI_DEALLOCATE(title)
  !call esymm_init(Bsym,kpt_rbz(:,ikpt),Cryst,.false.,nspinor,1,mband,tol5,eigen0,dtset%tolsym)
  !write(std_out,*) 'DEBUG : Bsym. ndegs=',Bsym%ndegs
  !do iband=1,Bsym%ndegs
  !  write(std_out,*) Bsym%degs_bounds(:,iband)
  !end do

  !call crystal_free(Cryst)
  !call esymm_free(Bsym)

 end subroutine check_degeneracies
!!***

!----------------------------------------------------------------------

!!****f* m_efmas/print_efmas
!! NAME
!! print_efmas
!!
!! FUNCTION
!! This routine prints the information needed to compute rapidly the band effective masses,
!! namely, the generalized second-order k-derivatives of the eigenenergies,
!! see Eq.(66) of Laflamme2016.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine print_efmas(efmasdeg,efmasval,kpt,ncid)

!Arguments ------------------------------------
!scalars
 integer,            intent(in) :: ncid
!arrays
 real(dp),            intent(in) :: kpt(:,:)
 type(efmasdeg_type), intent(in) :: efmasdeg(:)
 type(efmasval_type), intent(in) :: efmasval(:,:)

!Local variables-------------------------------
 integer :: deg_dim,eig2_diag_arr_dim
 integer :: iband,ideg,ideg_tot,ieig,ikpt
 integer :: jband,mband,ndegs_tot,nkpt,nkptdeg,nkptval
 integer, allocatable :: nband_arr(:)
integer, allocatable :: ndegs_arr(:)
 integer, allocatable :: degs_range_arr(:,:)
 integer, allocatable :: ideg_arr(:,:)
 integer, allocatable :: degs_bounds_arr(:,:)
 real(dp), allocatable :: ch2c_arr(:,:,:,:)
 real(dp), allocatable :: eig2_diag_arr(:,:,:,:)
 character(len=500) :: msg
#ifdef HAVE_NETCDF
 integer :: ncerr
#endif
!----------------------------------------------------------------------

!XG20180519 Here, suppose that dtset%nkpt=nkpt_rbz (as done by Jonathan).
!To be reexamined/corrected at the time of parallelization.

#ifdef HAVE_NETCDF
 nkptdeg=size(efmasdeg,1)
 nkptval=size(efmasval,2)
 if(nkptdeg/=nkptval) then
   write(msg,'(a,i8,a,i8,a)') ' nkptdeg and nkptval =',nkptdeg,' and ',nkptval,' differ, which is inconsistent.'
   MSG_ERROR(msg)
 end if
 nkpt=nkptdeg
 if(nkpt/=size(kpt,2)) then
   write(msg,'(a,i8,a,i8,a)') ' nkptdeg and nkpt =',nkptdeg,' and ',nkpt,' differ, which is inconsistent.'
   MSG_ERROR(msg)
 end if

 mband=size(efmasval,1)

!Total number of (degenerate) sets over all k points
 ndegs_tot=sum(efmasdeg%ndegs)
!Total number of generalized second-order k-derivatives
 eig2_diag_arr_dim=zero
 do ikpt=1,nkpt
   do ideg=efmasdeg(ikpt)%deg_range(1),efmasdeg(ikpt)%deg_range(2)
     deg_dim = efmasdeg(ikpt)%degs_bounds(2,ideg) - efmasdeg(ikpt)%degs_bounds(1,ideg) + 1
     eig2_diag_arr_dim = eig2_diag_arr_dim + deg_dim**2
   enddo
 enddo

!Allocate the arrays to be nc-written
 ABI_MALLOC(nband_arr, (nkpt) )
 ABI_MALLOC(ndegs_arr, (nkpt) )
 ABI_MALLOC(degs_range_arr, (2,nkpt) )
 ABI_MALLOC(ideg_arr, (mband,nkpt) )
 ABI_MALLOC(degs_bounds_arr, (2,ndegs_tot) )
 ABI_MALLOC(ch2c_arr, (2,3,3,eig2_diag_arr_dim) )
 ABI_MALLOC(eig2_diag_arr, (2,3,3,eig2_diag_arr_dim) )

!Prepare the arrays to be nc-written
 ideg_tot=1
 ieig=1
 do ikpt=1,nkpt
   nband_arr(ikpt)=efmasdeg(ikpt)%nband
   ndegs_arr(ikpt)=efmasdeg(ikpt)%ndegs
   degs_range_arr(:,ikpt)=efmasdeg(ikpt)%deg_range(:)
   ideg_arr(:,ikpt)=0
   ideg_arr(1:efmasdeg(ikpt)%nband,ikpt)=efmasdeg(ikpt)%ideg(:)
   do ideg=1,efmasdeg(ikpt)%ndegs
     degs_bounds_arr(:,ideg_tot)=efmasdeg(ikpt)%degs_bounds(:,ideg)
     ideg_tot=ideg_tot+1
   enddo
   do ideg=efmasdeg(ikpt)%deg_range(1),efmasdeg(ikpt)%deg_range(2)
     deg_dim = efmasdeg(ikpt)%degs_bounds(2,ideg) - efmasdeg(ikpt)%degs_bounds(1,ideg) + 1
     do jband=1,deg_dim
       do iband=1,deg_dim
         ch2c_arr(1,:,:,ieig+iband-1)=real(efmasval(ideg,ikpt)%ch2c(:,:,iband,jband))
         ch2c_arr(2,:,:,ieig+iband-1)=aimag(efmasval(ideg,ikpt)%ch2c(:,:,iband,jband))
         eig2_diag_arr(1,:,:,ieig+iband-1)=real(efmasval(ideg,ikpt)%eig2_diag(:,:,iband,jband))
         eig2_diag_arr(2,:,:,ieig+iband-1)=aimag(efmasval(ideg,ikpt)%eig2_diag(:,:,iband,jband))
       enddo
       ieig=ieig+deg_dim
     enddo
   enddo
 enddo

!Define dimensions
 ncerr=nctk_def_dims(ncid, [ &
&  nctkdim_t("number_of_reduced_dimensions", 3), &
&  nctkdim_t("real_or_complex", 2), &
&  nctkdim_t("number_of_kpoints", nkpt), &
&  nctkdim_t("max_number_of_states", mband), &
&  nctkdim_t("total_number_of_degenerate_sets", ndegs_tot), &
&  nctkdim_t("eig2_diag_arr_dim", eig2_diag_arr_dim)&
&  ], defmode=.True.)
 NCF_CHECK(ncerr)

 ncerr = nctk_def_arrays(ncid, [ &
& nctkarr_t("reduced_coordinates_of_kpoints", "dp", "number_of_reduced_dimensions, number_of_kpoints"), &
& nctkarr_t("number_of_states", "int", "number_of_kpoints"), &
& nctkarr_t("number_of_degenerate_sets", "int", "number_of_kpoints"), &
& nctkarr_t("degs_range_arr", "int", "two, number_of_kpoints"), &
& nctkarr_t("ideg_arr", "int", "max_number_of_states, number_of_kpoints"), &
& nctkarr_t("degs_bounds_arr", "int", "two, total_number_of_degenerate_sets"), &
& nctkarr_t("ch2c_arr", "dp", "real_or_complex, number_of_reduced_dimensions, number_of_reduced_dimensions, eig2_diag_arr_dim"),  &
& nctkarr_t("eig2_diag_arr","dp","real_or_complex, number_of_reduced_dimensions, number_of_reduced_dimensions, eig2_diag_arr_dim")&
  ])
 NCF_CHECK(ncerr)

 ! Write data.
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "reduced_coordinates_of_kpoints"), kpt))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "number_of_states"), nband_arr))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "number_of_degenerate_sets"), ndegs_arr))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "degs_range_arr"),            degs_range_arr))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ideg_arr"),                  ideg_arr))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "degs_bounds_arr"),           degs_bounds_arr))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ch2c_arr"),                  ch2c_arr))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eig2_diag_arr"),             eig2_diag_arr))

!Deallocate the arrays
 ABI_FREE(nband_arr)
 ABI_FREE(ndegs_arr)
 ABI_FREE(degs_range_arr)
 ABI_FREE(ideg_arr)
 ABI_FREE(degs_bounds_arr)
 ABI_FREE(ch2c_arr)
 ABI_FREE(eig2_diag_arr)
#endif

end subroutine print_efmas
!!***

!----------------------------------------------------------------------

!!****f* m_efmas/efmas_ncread
!! NAME
!! efmas_ncread
!!
!! FUNCTION
!! This routine reads from an EFMAS NetCDF file the information needed
!! to compute rapidly the band effective masses,
!! namely, the generalized second-order k-derivatives of the eigenenergies,
!! see Eq.(66) of Laflamme2016.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine efmas_ncread(efmasdeg,efmasval,kpt,ncid)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
!arrays
 real(dp), allocatable,intent(out) :: kpt(:,:)
 type(efmasdeg_type), allocatable, intent(out) :: efmasdeg(:)
 type(efmasval_type), allocatable, intent(out) :: efmasval(:,:)

!Local variables-------------------------------
 integer :: deg_dim,eig2_diag_arr_dim
 integer :: iband,ideg,ideg_tot,ieig,ikpt
 integer :: jband,mband,nband,ndegs,ndegs_tot,nkpt
 integer, allocatable :: nband_arr(:)
 integer, allocatable :: ndegs_arr(:)
 integer, allocatable :: degs_range_arr(:,:)
 integer, allocatable :: ideg_arr(:,:)
 integer, allocatable :: degs_bounds_arr(:,:)
 real(dp), allocatable :: ch2c_arr(:,:,:,:)
 real(dp), allocatable :: eig2_diag_arr(:,:,:,:)
!----------------------------------------------------------------------

#ifdef HAVE_NETCDF
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nctk_get_dim(ncid, "number_of_kpoints", nkpt))
 NCF_CHECK(nctk_get_dim(ncid, "max_number_of_states", mband))
 NCF_CHECK(nctk_get_dim(ncid, "total_number_of_degenerate_sets", ndegs_tot))
 NCF_CHECK(nctk_get_dim(ncid, "eig2_diag_arr_dim", eig2_diag_arr_dim))

!Allocate the arrays to be read from NetCDF file
 ABI_MALLOC(kpt, (3,nkpt) )
 ABI_MALLOC(nband_arr, (nkpt) )
 ABI_MALLOC(ndegs_arr, (nkpt) )
 ABI_MALLOC(degs_range_arr, (2,nkpt) )
 ABI_MALLOC(ideg_arr, (mband,nkpt) )
 ABI_MALLOC(degs_bounds_arr, (2,ndegs_tot) )
 ABI_MALLOC(ch2c_arr, (2,3,3,eig2_diag_arr_dim) )
 ABI_MALLOC(eig2_diag_arr, (2,3,3,eig2_diag_arr_dim) )

!Read from NetCDF file
 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "reduced_coordinates_of_kpoints"), kpt))
 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "number_of_states"), nband_arr))
 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "number_of_degenerate_sets"), ndegs_arr))
 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "degs_range_arr"),            degs_range_arr))
 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "ideg_arr"),                  ideg_arr))
 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "degs_bounds_arr"),           degs_bounds_arr))
 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "ch2c_arr"),                  ch2c_arr))
 NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "eig2_diag_arr"),             eig2_diag_arr))

!Prepare the efmas* datastructures
 ABI_DT_MALLOC(efmasdeg,(nkpt))
 ABI_DT_MALLOC(efmasval,(mband,nkpt))

 ideg_tot=1
 ieig=1
 do ikpt=1,nkpt
   efmasdeg(ikpt)%deg_range(:)=degs_range_arr(:,ikpt)
   nband=nband_arr(ikpt)
   efmasdeg(ikpt)%nband=nband
   ABI_MALLOC(efmasdeg(ikpt)%ideg, (nband))
   efmasdeg(ikpt)%ideg=ideg_arr(1:nband,ikpt)
   ndegs=ndegs_arr(ikpt)
   efmasdeg(ikpt)%ndegs=ndegs
   ABI_MALLOC(efmasdeg(ikpt)%degs_bounds,(2,nband))
   do ideg=1,ndegs
     efmasdeg(ikpt)%degs_bounds(:,ideg)=degs_bounds_arr(:,ideg_tot)
     ideg_tot=ideg_tot+1
     if( efmasdeg(ikpt)%deg_range(1) <= ideg .and. ideg <= efmasdeg(ikpt)%deg_range(2) ) then
       deg_dim=efmasdeg(ikpt)%degs_bounds(2,ideg) - efmasdeg(ikpt)%degs_bounds(1,ideg) + 1
       ABI_MALLOC(efmasval(ideg,ikpt)%ch2c,(3,3,deg_dim,deg_dim))
       ABI_MALLOC(efmasval(ideg,ikpt)%eig2_diag,(3,3,deg_dim,deg_dim))
       efmasval(ideg,ikpt)%ch2c=zero
       efmasval(ideg,ikpt)%eig2_diag=zero
       do jband=1,deg_dim
         do iband=1,deg_dim
           efmasval(ideg,ikpt)%ch2c(:,:,iband,jband)=&
&           dcmplx(ch2c_arr(1,:,:,ieig+iband-1),ch2c_arr(2,:,:,ieig+iband-1))
           efmasval(ideg,ikpt)%eig2_diag(:,:,iband,jband)=&
&           dcmplx(eig2_diag_arr(1,:,:,ieig+iband-1),eig2_diag_arr(2,:,:,ieig+iband-1))
         enddo
         ieig=ieig+deg_dim
       enddo
     else
       ABI_MALLOC(efmasval(ideg,ikpt)%ch2c,(0,0,0,0))
       ABI_MALLOC(efmasval(ideg,ikpt)%eig2_diag,(0,0,0,0))
     end if
   end do
 enddo

!Deallocate the arrays
 ABI_FREE(nband_arr)
 ABI_FREE(ndegs_arr)
 ABI_FREE(degs_range_arr)
 ABI_FREE(ideg_arr)
 ABI_FREE(degs_bounds_arr)
 ABI_FREE(ch2c_arr)
 ABI_FREE(eig2_diag_arr)
#endif

 end subroutine efmas_ncread

!!***

!----------------------------------------------------------------------

!!****f* m_efmas/print_tr_efmas
!! NAME
!! print_tr_efmas
!!
!! FUNCTION
!! This routine prints the transport equivalent effective mass and others info
!! for a degenerate set of bands
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_efmas
!!
!! CHILDREN
!!
!! SOURCE

 subroutine print_tr_efmas(io_unit,kpt,band,deg_dim,mdim,ndirs,dirs,m_cart,rprimd,efmas_tensor,ntheta, &
&                       m_avg,m_avg_frohlich,saddle_warn,efmas_eigval,efmas_eigvec,transport_tensor_scale)

   !Arguments ------------------------------------
   integer, intent(in) :: io_unit, band, deg_dim, mdim, ndirs
   real(dp), intent(in) :: m_cart(ndirs,deg_dim), kpt(3), dirs(3,ndirs), rprimd(3,3), efmas_tensor(mdim,mdim,deg_dim)
   integer, intent(in) :: ntheta
   real(dp), intent(in) :: m_avg(deg_dim),m_avg_frohlich(deg_dim)
   logical, intent(in) :: saddle_warn(deg_dim)
   real(dp), intent(in), optional :: efmas_eigval(mdim,deg_dim)
   real(dp), intent(in), optional :: efmas_eigvec(mdim,mdim,deg_dim)
   real(dp), intent(in), optional :: transport_tensor_scale(deg_dim)

   !Local variables ------------------------------
   logical :: extras
   integer :: iband, adir
   character(len=22) :: format_eigvec
   character(len=500) :: msg, tmpstr
   real(dp) :: vec(3)

   if(deg_dim>1) then
     extras = present(efmas_eigval) .and. present(efmas_eigvec)
     if(mdim==3 .and. .not. extras) then
       write(msg,'(a,l1,a,i1,a)') 'Subroutine print_tr_efmas called with degenerate=',deg_dim>1,&
&            ' and mdim=',mdim,', but missing required arguments for this case.'
       MSG_ERROR(msg)
     end if
     if(mdim==2 .and. .not. (extras .or. present(transport_tensor_scale))) then
       write(msg,'(a,l1,a,i1,a)') 'Subroutine print_tr_efmas called with degenerate=',deg_dim>1,&
&            ' and mdim=',mdim,', but missing required arguments for this case.'
       MSG_ERROR(msg)
     end if
   else
     extras = present(efmas_eigval) .and. present(efmas_eigvec)
     if(mdim>1 .and. .not. extras) then
       write(msg,'(a,l1,a,i1,a)') 'Subroutine print_tr_efmas called with degenerate=',deg_dim>1,&
&            ' and mdim=',mdim,', but missing required arguments for this case.'
       MSG_ERROR(msg)
     end if
   end if

   if(deg_dim>1) then
     write(io_unit,'(2a)') ch10,' COMMENTS: '
     write(io_unit,'(a,3(f6.3,a),i5,a,i5)') ' - At k-point (',kpt(1),',',kpt(2),',',kpt(3),'), bands ',band,' through ',&
&          band+deg_dim-1
     if(mdim>1) then
       write(io_unit,'(a)') '   are DEGENERATE (effective mass tensor is therefore not defined).'
       if(mdim==3) then
         write(io_unit,'(a)') '   See Section IIIB Eqs. (67)-(70) and Appendix E of PRB 93 205147 (2016).' ! [[cite:Laflamme2016]]
         write(io_unit,'(a,i7,a)') &
&          ' - Angular average effective mass for Frohlich model is to be averaged over degenerate bands. See later.'
       elseif(mdim==2) then
         write(io_unit,'(a)') ' - Also, 2D requested (perpendicular to Z axis).'
         write(io_unit,'(a)') '   See Section IIIB and Appendix F, Eqs. (F12)-(F14) of PRB 93 205147 (2016).' ! [[cite:Laflamme2016]]
       end if
       write(io_unit,'(a,i7,a)') ' - Associated theta integrals calculated with ntheta=',ntheta,' points.'
     else
       write(io_unit,'(a)') '   are DEGENERATE.'
       write(io_unit,'(a)') ' - Also, 1D requested (parallel to X axis).'
     end if
   end if

   if(ANY(saddle_warn)) then
     write(msg,'(2a)') ch10,'Band(s)'
     do iband=1,deg_dim
       if(saddle_warn(iband)) then
         write(tmpstr,'(i5)') band+iband-1
         msg = TRIM(msg)//' '//TRIM(tmpstr)//','
       end if
     end do
     write(tmpstr,'(6a)') ch10,'are not band extrema, but saddle points;',ch10, &
&                       'the transport equivalent formalism breaks down in these conditions.',ch10, &
&                       'The associated tensor(s) will therefore not be printed.'
     MSG_WARNING_UNIT(TRIM(msg)//TRIM(tmpstr), io_unit)
   end if

   if(deg_dim>1 .and. mdim>1) then
     write(msg,'(a)') ' Transport equivalent effective mass tensor'
   else
     write(msg,'(a)') ' Effective mass tensor'
   end if

   if(mdim>1) then
     write(format_eigvec,'(a,i1,a)') '(i3,',mdim,'f14.10,a,3f14.10)'
   end if

   do iband=1,deg_dim
     write(io_unit,'(2a,3(f6.3,a),i5)') ch10,' K-point (',kpt(1),',',kpt(2),',',kpt(3),') | band = ',band+iband-1
     write(io_unit,'(a)') trim(msg)//':'
     if(.not. saddle_warn(iband)) then
       do adir=1,mdim
         write(io_unit,'(3f26.10)') efmas_tensor(adir,:,iband)
       end do
       if(present(efmas_eigval))then
         write(io_unit,'(a)') trim(msg)//' eigenvalues:'
         write(io_unit,'(3f26.10)') efmas_eigval(:,iband)
       endif
     else
       write(io_unit,'(a)') '     *** SADDLE POINT: TRANSPORT EQV. EFF. MASS NOT DEFINED (see WARNING above) ***'
     end if

     if(mdim>1) then
       if(mdim==2 .and. deg_dim>1) then
         write(io_unit,'(a,f26.10)') 'Scaling of transport tensor (Eq. (FXX)) = ',transport_tensor_scale(iband)
       end if
       if(.not. saddle_warn(iband)) then
         if(io_unit == std_out) then
           write(io_unit,'(a)') trim(msg)//' eigenvectors in cartesian / reduced coord.:'
           do adir=1,mdim
             if( count( abs(efmas_eigval(adir,iband)-efmas_eigval(:,iband))<tol4 ) > 1 ) then
               write(io_unit,'(i3,a)') adir, ' Eigenvalue degenerate => eigenvector undefined'
             else
               vec=zero; vec(1:mdim)=efmas_eigvec(adir,:,iband)
               vec=matmul(transpose(rprimd)/two_pi,vec); vec=vec/sqrt(sum(vec**2))
               write(io_unit,format_eigvec) adir, efmas_eigvec(adir,:,iband), ' / ', vec
             end if
           end do
         end if
       end if
     end if

     if(mdim==3)then
       !An exactly zero average effective masse is artificial (or from a saddle point with symmetry). Does not print.
       if(abs(m_avg(iband))>tol8)then
         write(io_unit,'(a,f14.10)') &
&          ' Angular average effective mass 1/(<1/m>)= ',m_avg(iband)
       endif
       if(abs(m_avg_frohlich(iband))>tol8)then
         write(io_unit,'(a,f14.10)') &
&          ' Angular average effective mass for Frohlich model (<m**0.5>)**2= ',m_avg_frohlich(iband)
       endif
     endif

     write(io_unit,'(a)') ' Effective masses along directions: (cart. coord. / red. coord. -> eff. mass)'
     do adir=1,ndirs
       vec=dirs(:,adir)
       vec=matmul(transpose(rprimd)/two_pi,vec); vec=vec/sqrt(sum(vec**2))
       write(io_unit,'(i5,a,3f10.6,a,3f10.6,a,f14.10)') adir,': ', dirs(:,adir), ' / ', vec, ' -> ', m_cart(adir,iband)
     end do
   end do

   if(deg_dim>1 .and. mdim==3) then
     write(io_unit,'(2a)') ch10,&
&     ' Angular average effective mass for Frohlich model, averaged over degenerate bands.'
     write(io_unit,'(a,es16.6)') &
&     ' Value of     (<<m**0.5>>)**2 = ',(sum(abs(m_avg_frohlich(1:deg_dim))**0.5)/deg_dim)**2
     write(io_unit,'(a,es16.6,a)') &
&     ' Absolute Value of <<m**0.5>> = ', sum(abs(m_avg_frohlich(1:deg_dim))**0.5)/deg_dim,ch10
   endif

 end subroutine print_tr_efmas
!!***

!----------------------------------------------------------------------

!!****f* m_efmas/efmas_main
!! NAME
!! efmas_main
!!
!! FUNCTION
!! This routine calculates the generalized second-order k-derivative, Eq.66 of Laflamme2016,
!! in reduced coordinates.
!!
!! INPUTS
!!  cg(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%nsppol*nkpt_rbz)=pw coefficients of GS wavefunctions at k.
!!  cg1_pert(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%nsppo*nkpt_rbz,3,mpert) = first-order wf in G
!!            space for each perturbation. The wavefunction is orthogonal to the
!!            active space.
!!  dim_eig2rf = 1 if cg1_pert, gh0c1_pert and gh1c_pert are allocated.
!!               0 otherwise.
!!  dtset = dataset structure containing the input variable of the calculation.
!!  efmasdeg(nkpt_rbz) <type(efmasdeg_type)>= information about the band degeneracy at each k point
!!  eigen0(nkpt_rbz*dtset%mband*dtset%nsppol) = 0-order eigenvalues at all K-points:
!!            <k,n'|H(0)|k,n'> (hartree).
!!  eigen1(nkpt_rbz*2*dtset%nsppol*dtset%mband**2,3,mpert) = matrix of first-order:
!!            <k+Q,n'|H(1)|k,n> (hartree) (calculated in dfpt_cgwf).
!!  gh0c1_pert(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%nsppol*nkpt_rbz,3,mpert) = matrix containing the
!!            vector:  <G|H(0)|psi(1)>, for each perturbation.
!!  gh1c_pert(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%nsppol*nkpt_rbz,3,mpert)) = matrix containing the
!!            vector:  <G|H(1)|n,k>, for each perturbation. The wavefunction is
!!            orthogonal to the active space.
!!  istwfk_pert(nkpt_rbz,3,mpert) = integer for choice of storage of wavefunction at
!!            each k point for each perturbation.
!!  mpert = maximum number of perturbations.
!!  mpi_enreg = informations about MPI parallelization.
!!  nkpt_rbz = number of k-points for each perturbation.
!!  npwarr(nkpt_rbz,mpert) = array of numbers of plane waves for each k-point
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  efmasval(mband,nkpt_rbz) <type(efmasdeg_type)>= generalized 2nd-order k-derivatives of eigenvalues
!!    efmasval(:,:)%ch2c INPUT : frozen wavefunction H2 contribution to generalized 2nd order k-derivatives of eigenenergy
!!    efmasval(:,:)%eig2_diag OUTPUT : generalized 2nd order k-derivatives of eigenenergy
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!      cgqf,dgemm,dgetrf,dgetri,dotprod_g,dsyev,print_tr_efmas,zgemm,zgetrf
!!      zgetri,zheev
!!
!! SOURCE

 subroutine efmas_main(cg,cg1_pert,dim_eig2rf,dtset,efmasdeg,efmasval,eigen0,&
&   eigen1,gh0c1_pert,gh1c_pert,istwfk_pert,mpert,mpi_enreg,nkpt_rbz,npwarr,rprimd)

 !Arguments ------------------------------------
 !scalars
  integer,            intent(in)    :: dim_eig2rf,mpert,nkpt_rbz
  type(dataset_type), intent(in)    :: dtset
  type(MPI_type),     intent(in) :: mpi_enreg
 !arrays
  integer,  intent(in) :: istwfk_pert(nkpt_rbz,3,mpert)
  integer,  intent(in) :: npwarr(nkpt_rbz,mpert)
  real(dp), intent(in) :: cg1_pert(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%nsppol*nkpt_rbz*dim_eig2rf,3,mpert)
  real(dp), intent(in) :: gh0c1_pert(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%nsppol*nkpt_rbz*dim_eig2rf,3,mpert)
  real(dp), intent(in) :: gh1c_pert(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%nsppol*nkpt_rbz*dim_eig2rf,3,mpert)
  real(dp), intent(in) :: eigen0(nkpt_rbz*dtset%mband*dtset%nsppol)
  real(dp), intent(in) :: eigen1(nkpt_rbz*2*dtset%nsppol*dtset%mband**2,3,mpert)
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(in) :: cg(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%nsppol*nkpt_rbz)
  type(efmasdeg_type), allocatable,intent(in) :: efmasdeg(:)
  type(efmasval_type),  allocatable,intent(inout) :: efmasval(:,:)

 !Local variables-------------------------------
  logical :: degenerate
  logical :: debug
  integer :: ipert
  integer :: isppol
  integer :: icg2   !TODOM : Reactivate the sections for icg2 / allow choice of k-point other than the first in the list.
  integer :: npw_k
  integer :: nband_k
  integer :: nspinor
  integer :: ideg
  integer :: ikpt
  integer :: istwf_k
  integer :: master,me,spaceworld
  integer :: band2tot_index
  integer :: bandtot_index
  integer :: iband, jband, kband
  integer :: adir,bdir
  integer :: deg_dim
  integer :: degl
  character(len=500) :: msg
  real(dp) :: deltae
  real(dp) :: dot2i,dot2r,dot3i,dot3r,doti,dotr
  real(dp), allocatable :: cg0(:,:)
  real(dp), allocatable :: cg1_pert2(:,:),cg1_pert1(:,:)
  real(dp), allocatable :: gh1c_pert2(:,:),gh1c_pert1(:,:),gh0c1_pert1(:,:)
  complex(dpc) :: eig2_part(3,3)
  complex(dpc) :: eig2_ch2c(3,3)
  complex(dpc) :: eig2_paral(3,3)
  complex(dpc) :: eig2_gauge_change(3,3)
  complex(dpc) :: eig1a, eig1b, g_ch
  complex(dpc), allocatable :: eigen1_deg(:,:)
  complex(dpc), allocatable :: eig2_diag(:,:,:,:)
  complex(dpc), allocatable :: eig2_diag_cart(:,:,:,:)

! *********************************************************************

  debug = .false. ! Prints additional info to std_out

! Init parallelism
  master =0
  spaceworld=mpi_enreg%comm_cell
  me=mpi_enreg%me_kpt

  write(msg,'(4a)') ch10,&
&  ' CALCULATION OF EFFECTIVE MASSES',ch10,&
&  ' NOTE : Additional infos (eff. mass eigenvalues, eigenvectors and, if degenerate, average mass) are available in stdout.'
  call wrtout(std_out,msg,'COLL')
  call wrtout(ab_out,msg,'COLL')

  if(dtset%nsppol/=1)then
    write(msg,'(a,i3,a)') 'nsppol=',dtset%nsppol,' is not yet treated in m_efmas.'
    MSG_ERROR(msg)
  end if
  if(dtset%nspden/=1)then
    write(msg,'(a,i3,a)') 'nspden=',dtset%nspden,' is not yet treated in m_efmas.'
    MSG_ERROR(msg)
  end if
  if(dtset%efmas_deg==0) then
    write(msg,'(a)') 'efmas_deg==0 is for debugging; the results for degenerate bands will be garbage.'
    MSG_WARNING(msg)
    MSG_WARNING_UNIT(msg, ab_out)
  end if

  ipert = dtset%natom+1
  isppol = 1

  icg2 = 0
  band2tot_index=0
  bandtot_index=0

  !XG20180519 : in the original coding by Jonathan, there is a lack of care about using dtset%nkpt or nkpt_rbz ...
  !Not important in the sequential case (?!) but likely problematic in the parallel case.
  do ikpt=1,dtset%nkpt
    npw_k = npwarr(ikpt,ipert)
    nband_k = dtset%nband(ikpt)
    nspinor = dtset%nspinor

    ABI_ALLOCATE(cg1_pert2,(2,npw_k*nspinor))
    ABI_ALLOCATE(cg1_pert1,(2,npw_k*nspinor))
    ABI_ALLOCATE(gh1c_pert2,(2,npw_k*nspinor))
    ABI_ALLOCATE(gh1c_pert1,(2,npw_k*nspinor))
    ABI_ALLOCATE(gh0c1_pert1,(2,npw_k*nspinor))
    ABI_ALLOCATE(cg0,(2,npw_k*nspinor))

    do ideg=efmasdeg(ikpt)%deg_range(1),efmasdeg(ikpt)%deg_range(2)
      deg_dim    = efmasdeg(ikpt)%degs_bounds(2,ideg) - efmasdeg(ikpt)%degs_bounds(1,ideg) + 1
      degenerate = (deg_dim>1) .and. (dtset%efmas_deg/=0)
      degl       = efmasdeg(ikpt)%degs_bounds(1,ideg)-1

      ABI_ALLOCATE(eigen1_deg,(deg_dim,deg_dim))
      !!! If treated band degenerate at 0th order, check that we are at extrema.
      if(degenerate) then
        do adir=1,3
          do iband=1,deg_dim
            do jband=1,deg_dim
              eigen1_deg(iband,jband) = cmplx(eigen1(2*(jband+degl)-1+(iband+degl-1)*2*nband_k,adir,ipert),&
&              eigen1(2*(jband+degl)+(iband+degl-1)*2*nband_k,adir,ipert),dpc)
            end do
          end do
          if (.not.(ALL(ABS(eigen1_deg)<tol5))) then
            write(msg,'(a,a)') ' Effective masses calculations require given k-point(s) to be band extrema for given bands, ',&
&            'but gradient of band(s) was found to be nonzero.'
            MSG_ERROR(msg)
          end if
        end do !adir=1,3
      end if !degenerate(1)
      ABI_DEALLOCATE(eigen1_deg)

      ABI_ALLOCATE(eig2_diag,(3,3,deg_dim,deg_dim))
      ABI_ALLOCATE(eig2_diag_cart,(3,3,deg_dim,deg_dim))
      eig2_diag = zero

      do iband=1,deg_dim
        write(std_out,*)"  In the set (possibly degenerate) compute band ",iband  ! This line here to avoid weird
        cg0(:,:) = cg(:,1+(degl+iband-1)*npw_k*nspinor+icg2:(degl+iband)*npw_k*nspinor+icg2)
        do jband=1,deg_dim
          eig2_part = zero
          eig2_ch2c = zero
          eig2_paral = zero
          eig2_gauge_change = zero
          do adir=1,3
            istwf_k = istwfk_pert(ikpt,adir,ipert)
            do bdir=1,3

              ! Calculate the gauge change (to be subtracted to go from parallel transport to diagonal gauge).
              do kband=1,nband_k
                !!! Equivalent to the gauge change in eig2stern.F90, but works also for other choices than the parallel gauge.
                eig1a = cmplx( eigen1(2*kband-1+(degl+iband-1)*2*nband_k+band2tot_index,adir,ipert), &
&                -eigen1(2*kband+(degl+iband-1)*2*nband_k+band2tot_index,adir,ipert), kind=dpc )
                eig1b = cmplx( eigen1(2*kband-1+(degl+jband-1)*2*nband_k+band2tot_index,bdir,ipert), &
&                eigen1(2*kband+(degl+jband-1)*2*nband_k+band2tot_index,bdir,ipert), kind=dpc )
                g_ch = eig1a*eig1b
                eig1a = cmplx( eigen1(2*kband-1+(degl+iband-1)*2*nband_k+band2tot_index,bdir,ipert), &
&                -eigen1(2*kband+(degl+iband-1)*2*nband_k+band2tot_index,bdir,ipert), kind=dpc )
                eig1b = cmplx( eigen1(2*kband-1+(degl+jband-1)*2*nband_k+band2tot_index,adir,ipert), &
&                eigen1(2*kband+(degl+jband-1)*2*nband_k+band2tot_index,adir,ipert), kind=dpc )
                g_ch = g_ch + eig1a*eig1b

                deltae = eigen0(kband+bandtot_index) - eigen0((degl+iband)+bandtot_index)
                if( kband<=degl.or.kband>degl+deg_dim) then
                  g_ch = g_ch/deltae
                else
                  g_ch = zero
                end if
                eig2_gauge_change(adir,bdir) = eig2_gauge_change(adir,bdir) + g_ch
              end do !kband

              cg1_pert2(:,:)   = cg1_pert(:,1+(degl+jband-1)*npw_k*nspinor+icg2:(degl+jband)*npw_k*nspinor+icg2,bdir,ipert)
              cg1_pert1(:,:)   = cg1_pert(:,1+(degl+iband-1)*npw_k*nspinor+icg2:(degl+iband)*npw_k*nspinor+icg2,adir,ipert)
              gh1c_pert1(:,:)  = gh1c_pert(:,1+(degl+iband-1)*npw_k*nspinor+icg2:(degl+iband)*npw_k*nspinor+icg2,adir,ipert)
              gh1c_pert2(:,:)  = gh1c_pert(:,1+(degl+jband-1)*npw_k*nspinor+icg2:(degl+jband)*npw_k*nspinor+icg2,bdir,ipert)
              gh0c1_pert1(:,:) = gh0c1_pert(:,1+(degl+iband-1)*npw_k*nspinor+icg2:(degl+iband)*npw_k*nspinor+icg2,adir,ipert)

              ! The first two dotprod corresponds to:  <Psi(1)|H(1)|Psi(0)> + cc.
              ! They are calculated using wavefunctions <Psi(1)| that are orthogonal to the active space.
              dotr=zero ; doti=zero
              call dotprod_g(dotr,doti,istwf_k,npw_k*nspinor,2,cg1_pert1,gh1c_pert2,mpi_enreg%me_g0,&
&              mpi_enreg%comm_spinorfft)
              dot2r=zero ; dot2i=zero
              call dotprod_g(dot2r,dot2i,istwf_k,npw_k*nspinor,2,gh1c_pert1,cg1_pert2,mpi_enreg%me_g0,&
&              mpi_enreg%comm_spinorfft)

              ! This dotprod corresponds to : <Psi(1)|H(0)- E(0)|Psi(1)>
              ! It is calculated using wavefunctions that are orthogonal to the active space.
              dot3r=zero ; dot3i=zero
              call dotprod_g(dot3r,dot3i,istwf_k,npw_k*nspinor,2,gh0c1_pert1,cg1_pert2,mpi_enreg%me_g0,&
&              mpi_enreg%comm_spinorfft)

              eig2_part(adir,bdir) = cmplx(dotr+dot2r+dot3r,doti+dot2i+dot3i,kind=dpc)
              !eig2_part(adir,bdir) = cmplx(dotr+dot2r,doti+dot2i,kind=dpc)  !DEBUG
              !eig2_part(adir,bdir) = cmplx(dotr,doti,kind=dpc)              !DEBUG

              eig2_ch2c(adir,bdir) = efmasval(ideg,ikpt)%ch2c(adir,bdir,iband,jband)

            end do !bdir
          end do  !adir

          do adir=1,3
            do bdir=1,3
              eig2_paral(adir,bdir) = eig2_part(adir,bdir) + eig2_part(bdir,adir) + eig2_ch2c(adir,bdir)
            end do
          end do

          eig2_diag(:,:,iband,jband) = eig2_paral - eig2_gauge_change
          efmasval(ideg,ikpt)%eig2_diag(:,:,iband,jband)=eig2_diag(:,:,iband,jband)

          !!! Decomposition of the hessian in into its different contributions.
          if(debug) then

            eig2_diag_cart(:,:,iband,jband) = matmul(matmul(rprimd,eig2_diag(:,:,iband,jband)),transpose(rprimd))/two_pi**2
            eig2_paral                 = matmul(matmul(rprimd,eig2_paral),                transpose(rprimd))/two_pi**2
            eig2_gauge_change          = matmul(matmul(rprimd,eig2_gauge_change),         transpose(rprimd))/two_pi**2
            eig2_ch2c                  = matmul(matmul(rprimd,eig2_ch2c),                 transpose(rprimd))/two_pi**2
            eig2_part                  = matmul(matmul(rprimd,eig2_part),                 transpose(rprimd))/two_pi**2

            write(std_out,'(a)') 'Hessian of eigenvalues               = H. in parallel gauge - Gauge transformation'
            do adir=1,3
              write(std_out,'(3f12.8,2(a,3f12.8))')&
&               real(eig2_diag_cart(adir,:,iband,jband),dp),' |',real(eig2_paral(adir,:),dp),' |',real(eig2_gauge_change(adir,:),dp)
            end do
            write(std_out,'(a)') 'H. in parallel gauge  = Second der. of H     + First derivatives    + First derivatives^T'
            do adir=1,3
              write(std_out,'(3f12.8,2(a,3f12.8))') real(eig2_paral(adir,:),dp),' |',real(eig2_ch2c(adir,:),dp),' |', &
&                                                   real(eig2_part(adir,:),dp)
            end do
          end if !debug

          if(.not. degenerate .and. iband==jband) then
            write(std_out,'(a,3f20.16)') 'Gradient of eigenvalues = ',&
&            matmul(rprimd,eigen1(2*(degl+iband)-1+(degl+iband-1)*2*nband_k+band2tot_index,:,ipert))/two_pi
          end if !.not.degenerate

        end do !jband
      end do !iband

      ABI_DEALLOCATE(eig2_diag)
      ABI_DEALLOCATE(eig2_diag_cart)

    end do !ideg

    ABI_DEALLOCATE(cg1_pert2)
    ABI_DEALLOCATE(cg1_pert1)
    ABI_DEALLOCATE(gh1c_pert2)
    ABI_DEALLOCATE(gh1c_pert1)
    ABI_DEALLOCATE(gh0c1_pert1)
    ABI_DEALLOCATE(cg0)

    icg2=icg2+npw_k*dtset%nspinor*nband_k
    bandtot_index=bandtot_index+nband_k
    band2tot_index=band2tot_index+2*nband_k**2
  end do ! ikpt

 end subroutine efmas_main
!!***

!----------------------------------------------------------------------

!!****f* m_efmas/efmas_analysis
!! NAME
!! efmas_analysis
!!
!! FUNCTION
!! This routine analyzes the generalized second-order k-derivatives of eigenenergies,
!! and compute the effective mass tensor
!! (inverse of hessian of eigenvalues with respect to the wavevector)
!! in cartesian coordinates along different directions in k-space, or also the transport equivalent effective mass.
!!
!! INPUTS
!!  dtset = dataset structure containing the input variable of the calculation.
!!  efmasdeg(nkpt_rbz) <type(efmasdeg_type)>= information about the band degeneracy at each k point
!!  efmasval(mband,nkpt_rbz) <type(efmasdeg_type)>= double tensor datastructure
!!    efmasval(:,:)%eig2_diag band curvature double tensor
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points.
!!  mpert = maximum number of perturbations.
!!  mpi_enreg = informations about MPI parallelization.
!!  nkpt_rbz = number of k-points for each perturbation.
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!
!! OUTPUT
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!      cgqf,dgemm,dgetrf,dgetri,dotprod_g,dsyev,print_tr_efmas,zgemm,zgetrf
!!      zgetri,zheev
!!
!! SOURCE

 subroutine efmas_analysis(dtset,efmasdeg,efmasval,kpt_rbz,mpi_enreg,nkpt_rbz,rprimd)

 !Arguments ------------------------------------
 !scalars
  integer,            intent(in)    :: nkpt_rbz
  type(dataset_type), intent(in)    :: dtset
  type(MPI_type),     intent(in) :: mpi_enreg
 !arrays
  real(dp), intent(in) :: rprimd(3,3)
  real(dp), intent(in) :: kpt_rbz(3,nkpt_rbz)
  type(efmasdeg_type), intent(in) :: efmasdeg(:)
  type(efmasval_type), intent(in) :: efmasval(:,:)

 !Local variables-------------------------------
  logical :: degenerate
  logical :: debug
  logical :: print_fsph
  logical, allocatable :: saddle_warn(:), start_eigf3d_pos(:)
  integer :: info
  integer :: isppol
  integer :: ideg,jdeg
  integer :: ikpt
  integer :: master,me,spaceworld
  integer :: iband, jband
  integer :: adir,bdir
  integer :: deg_dim
  integer :: degl
  integer :: lwork
  integer :: itheta, iphi
  integer :: ntheta, nphi
  integer :: mdim
  integer :: cdirs, ndirs
  integer :: ipiv(3)
  integer :: io_unit
  character(len=500) :: msg, filename
  real(dp) :: cosph,costh,sinph,sinth
  real(dp) :: f3d_scal
  real(dp) :: weight
  real(dp) :: gprimd(3,3)
  real(dp), allocatable :: unit_r(:), dr_dth(:), dr_dph(:)
  real(dp), allocatable :: eigenval(:), rwork(:)
  real(dp), allocatable :: eigf3d(:)
  real(dp), allocatable :: m_avg(:), m_avg_frohlich(:),m_cart(:,:)
  real(dp), allocatable :: deigf3d_dth(:), deigf3d_dph(:)
  real(dp), allocatable :: unit_speed(:,:), transport_tensor(:,:,:)
  real(dp), allocatable :: cart_rotation(:,:), transport_tensor_eig(:)
  real(dp), allocatable :: transport_eqv_m(:,:,:), transport_eqv_eigval(:,:), transport_eqv_eigvec(:,:,:)
  real(dp), allocatable :: transport_tensor_scale(:)
  real(dp), allocatable :: gq_points_th(:),gq_points_costh(:),gq_points_sinth(:),gq_weights_th(:)
  real(dp), allocatable :: gq_points_ph(:),gq_points_cosph(:),gq_points_sinph(:),gq_weights_ph(:)
  real(dp), allocatable :: dirs(:,:)
  real(dp),allocatable :: prodr(:,:)
  !real(dp), allocatable :: f3dfd(:,:,:)
  complex(dpc) :: matr2d(2,2)
  complex(dpc), allocatable :: eigenvec(:,:), work(:)
  complex(dpc), allocatable :: eig2_diag_cart(:,:,:,:)
  complex(dpc), allocatable :: f3d(:,:), df3d_dth(:,:), df3d_dph(:,:)
  complex(dpc), allocatable :: unitary_tr(:,:), eff_mass(:,:)
  complex(dpc),allocatable :: prodc(:,:)

! *********************************************************************

  debug = .false. ! Prints additional info to std_out
  print_fsph = .false. ! Open a file and print the angle dependent curvature f(\theta,\phi)
                       ! for each band & kpts treated; 1 file per degenerate ensemble of bands.
                       ! Angles are those used in the numerical integration.

! Init parallelism
  master =0
  spaceworld=mpi_enreg%comm_cell
  me=mpi_enreg%me_kpt

  isppol = 1

  mdim = dtset%efmas_dim

!HERE ALLOCATE

  gprimd = rprimd
  call dgetrf(mdim,mdim,gprimd,mdim,ipiv,info)
  ABI_ALLOCATE(rwork,(3))
  call dgetri(mdim,gprimd,mdim,ipiv,rwork,3,info)
  ABI_DEALLOCATE(rwork)
  gprimd = two_pi*transpose(gprimd)

  cdirs = dtset%efmas_calc_dirs
  ndirs = mdim
  if(cdirs/=0) ndirs = dtset%efmas_n_dirs
  ABI_ALLOCATE(dirs,(3,ndirs))
  if(cdirs==0) then
    dirs = zero
    do adir=1,ndirs
      dirs(adir,adir)=1.0_dp
    end do
  elseif(cdirs==1) then
    dirs(:,:) = dtset%efmas_dirs(:,1:ndirs)
    do adir=1,ndirs
      dirs(:,adir) = dirs(:,adir)/sqrt(sum(dirs(:,adir)**2))
    end do
  elseif(cdirs==2) then
    dirs(:,:) = matmul(gprimd,dtset%efmas_dirs(:,1:ndirs))
    do adir=1,ndirs
      dirs(:,adir) = dirs(:,adir)/sqrt(sum(dirs(:,adir)**2))
    end do
  elseif(cdirs==3) then
    dirs(1,:) = sin(dtset%efmas_dirs(1,1:ndirs)*pi/180)*cos(dtset%efmas_dirs(2,1:ndirs)*pi/180)
    dirs(2,:) = sin(dtset%efmas_dirs(1,1:ndirs)*pi/180)*sin(dtset%efmas_dirs(2,1:ndirs)*pi/180)
    dirs(3,:) = cos(dtset%efmas_dirs(1,1:ndirs)*pi/180)
  end if

  !!! Initialization of integrals for the degenerate case.
  ntheta   = dtset%efmas_ntheta
  nphi     = 2*ntheta
  ABI_ALLOCATE(gq_points_th,(ntheta))
  ABI_ALLOCATE(gq_points_costh,(ntheta))
  ABI_ALLOCATE(gq_points_sinth,(ntheta))
  ABI_ALLOCATE(gq_weights_th,(ntheta))
  ABI_ALLOCATE(gq_points_ph,(nphi))
  ABI_ALLOCATE(gq_points_cosph,(nphi))
  ABI_ALLOCATE(gq_points_sinph,(nphi))
  ABI_ALLOCATE(gq_weights_ph,(nphi))
  call cgqf(ntheta,1,zero,zero,zero,pi,gq_points_th,gq_weights_th)
  !XG180501 : TODO : There is no need to make a Gauss-Legendre integral for the phi variable,
  !since the function to be integrated is periodic...
  call cgqf(nphi,1,zero,zero,zero,2*pi,gq_points_ph,gq_weights_ph)
  do itheta=1,ntheta
    gq_points_costh(itheta)=cos(gq_points_th(itheta))
    gq_points_sinth(itheta)=sin(gq_points_th(itheta))
  enddo
  do iphi=1,nphi
    gq_points_cosph(iphi)=cos(gq_points_ph(iphi))
    gq_points_sinph(iphi)=sin(gq_points_ph(iphi))
  enddo

  ABI_ALLOCATE(eff_mass,(mdim,mdim))

!XG20180519 : incoherent, efmasdeg is dimensioned at nkpt_rbz, and not at dtset%nkpt ...
  do ikpt=1,dtset%nkpt
    do ideg=efmasdeg(ikpt)%deg_range(1),efmasdeg(ikpt)%deg_range(2)

     deg_dim    = efmasdeg(ikpt)%degs_bounds(2,ideg) - efmasdeg(ikpt)%degs_bounds(1,ideg) + 1
     degenerate = (deg_dim>1) .and. (dtset%efmas_deg/=0)
     degl       = efmasdeg(ikpt)%degs_bounds(1,ideg)-1

     !!! Allocations
     ABI_ALLOCATE(eigenvec,(deg_dim,deg_dim))
     ABI_ALLOCATE(eigenval,(deg_dim))

     ABI_ALLOCATE(eig2_diag_cart,(3,3,deg_dim,deg_dim))

     do iband=1,deg_dim
        write(std_out,*)"  Compute band ",iband  ! This line here to avoid weird
        do jband=1,deg_dim

          eff_mass=zero

          eig2_diag_cart(:,:,iband,jband)=efmasval(ideg,ikpt)%eig2_diag(:,:,iband,jband)
          eig2_diag_cart(:,:,iband,jband) = matmul(matmul(rprimd,eig2_diag_cart(:,:,iband,jband)),transpose(rprimd))/two_pi**2


          if(.not. degenerate .and. iband==jband) then

            !Compute effective mass tensor from second derivative matrix. Simple inversion.
            eff_mass(:,:) = eig2_diag_cart(1:mdim,1:mdim,iband,jband)
            call zgetrf(mdim,mdim,eff_mass(1:mdim,1:mdim),mdim,ipiv,info)
            ABI_ALLOCATE(work,(3))
            call zgetri(mdim,eff_mass(1:mdim,1:mdim),mdim,ipiv,work,3,info)
            ABI_DEALLOCATE(work)

            !DIAGONALIZATION
            ABI_ALLOCATE(transport_eqv_eigvec,(mdim,mdim,deg_dim))
            transport_eqv_eigvec=zero
            ABI_ALLOCATE(transport_eqv_eigval,(mdim,deg_dim))
            transport_eqv_eigval=zero
            transport_eqv_eigvec(:,:,iband) = real(eff_mass(1:mdim,1:mdim),dp)
            lwork=-1
            ABI_ALLOCATE(rwork,(1))
            call dsyev('V','U',mdim,transport_eqv_eigvec(:,:,iband),mdim,transport_eqv_eigval(:,iband),rwork,lwork,info)
            lwork = max(1, 3*mdim-1) ! lwork >= max(1, 3*mdim-1)
            ABI_DEALLOCATE(rwork)

            ABI_ALLOCATE(rwork,(lwork))
            rwork=zero
            call dsyev('V','U',mdim,transport_eqv_eigvec(:,:,iband),mdim,transport_eqv_eigval(:,iband),rwork,lwork,info)
            ABI_DEALLOCATE(rwork)
            transport_eqv_eigvec(:,:,iband) = transpose(transport_eqv_eigvec(:,:,iband)) !So that lines contain eigenvectors.

            !Frohlich average effective mass
            ABI_ALLOCATE(m_avg,(1))
            ABI_ALLOCATE(m_avg_frohlich,(1))
            ABI_ALLOCATE(saddle_warn,(1))
            ABI_ALLOCATE(unit_r,(mdim))
            ABI_ALLOCATE(start_eigf3d_pos,(1))

            m_avg=zero
            m_avg_frohlich=zero
            saddle_warn=.false.

            if(mdim==3)then
              !One has to perform the integral over the sphere
              do itheta=1,ntheta
                costh=gq_points_costh(itheta) ; sinth=gq_points_sinth(itheta)
                do iphi=1,nphi
                  cosph=gq_points_cosph(iphi) ; sinph=gq_points_sinph(iphi)
                  weight=gq_weights_th(itheta)*gq_weights_ph(iphi)

                  unit_r(1)=sinth*cosph
                  unit_r(2)=sinth*sinph
                  unit_r(3)=costh

                  f3d_scal=dot_product(unit_r(:),matmul(real(eig2_diag_cart(:,:,iband,jband),dp),unit_r(:)))
                  m_avg = m_avg + weight*sinth*f3d_scal
                  m_avg_frohlich = m_avg_frohlich + weight*sinth/(abs(f3d_scal)**half)

                  if(itheta==1 .and. iphi==1) start_eigf3d_pos = f3d_scal > 0
                  if(start_eigf3d_pos(1) .neqv. (f3d_scal>0)) then
                    saddle_warn(1)=.true.
                  end if
                enddo
              enddo
              m_avg = quarter/pi*m_avg
              m_avg = one/m_avg
              m_avg_frohlich = quarter/pi*m_avg_frohlich
              m_avg_frohlich = m_avg_frohlich**2
              m_avg_frohlich(1) = DSIGN(m_avg_frohlich(1),m_avg(1))

            endif ! mdim==3

            !EFMAS_DIRS
            ABI_ALLOCATE(m_cart,(ndirs,deg_dim))
            m_cart=zero
            do adir=1,ndirs
              m_cart(adir,1)=1.0_dp/dot_product(dirs(:,adir),matmul(real(eig2_diag_cart(:,:,iband,jband),dp),dirs(:,adir)))
            end do

            !PRINTING RESULTS
            call print_tr_efmas(std_out,kpt_rbz(:,ikpt),degl+iband,1,mdim,ndirs,dirs,m_cart,rprimd,real(eff_mass,dp), &
&             ntheta,m_avg,m_avg_frohlich,saddle_warn,&
&             transport_eqv_eigval(:,iband:iband),transport_eqv_eigvec(:,:,iband:iband))
            call print_tr_efmas(ab_out, kpt_rbz(:,ikpt),degl+iband,1,mdim,ndirs,dirs,m_cart,rprimd,real(eff_mass,dp), &
&             ntheta,m_avg,m_avg_frohlich,saddle_warn,&
&             transport_eqv_eigval(:,iband:iband),transport_eqv_eigvec(:,:,iband:iband))
            ABI_DEALLOCATE(m_cart)
            ABI_DEALLOCATE(transport_eqv_eigvec)
            ABI_DEALLOCATE(transport_eqv_eigval)
            ABI_DEALLOCATE(m_avg)
            ABI_DEALLOCATE(m_avg_frohlich)
            ABI_DEALLOCATE(unit_r)
            ABI_DEALLOCATE(saddle_warn)
            ABI_DEALLOCATE(start_eigf3d_pos)

          end if !.not.degenerate
        end do !jband
      end do !iband

      !!! EQV_MASS
      if(degenerate .and. mdim==3) then
        ABI_ALLOCATE(unit_r,(mdim))
        ABI_ALLOCATE(dr_dth,(mdim))
        ABI_ALLOCATE(dr_dph,(mdim))
        ABI_ALLOCATE(f3d,(deg_dim,deg_dim))
        ABI_ALLOCATE(df3d_dth,(deg_dim,deg_dim))
        ABI_ALLOCATE(df3d_dph,(deg_dim,deg_dim))
        ABI_ALLOCATE(unitary_tr,(deg_dim,deg_dim))
        ABI_ALLOCATE(eigf3d,(deg_dim))
        ABI_ALLOCATE(saddle_warn,(deg_dim))
        ABI_ALLOCATE(start_eigf3d_pos,(deg_dim))
        ABI_ALLOCATE(m_avg,(deg_dim))
        ABI_ALLOCATE(m_avg_frohlich,(deg_dim))
        ABI_ALLOCATE(m_cart,(ndirs,deg_dim))
        ABI_ALLOCATE(deigf3d_dth,(deg_dim))
        ABI_ALLOCATE(deigf3d_dph,(deg_dim))
        ABI_ALLOCATE(unit_speed,(mdim,deg_dim))
        ABI_ALLOCATE(transport_tensor,(mdim,mdim,deg_dim))
        ABI_ALLOCATE(transport_tensor_eig,(mdim))
        ABI_ALLOCATE(transport_eqv_m,(mdim,mdim,deg_dim))
        ABI_ALLOCATE(transport_eqv_eigval,(mdim,deg_dim))
        ABI_ALLOCATE(transport_eqv_eigvec,(mdim,mdim,deg_dim))
        ABI_ALLOCATE(prodc,(deg_dim,deg_dim))
        ABI_ALLOCATE(prodr,(mdim,mdim))
        !ABI_ALLOCATE(f3dfd,(2,nphi,deg_dim))
        unit_r=zero
        dr_dth=zero
        dr_dph=zero
        f3d=zero
        df3d_dth=zero
        df3d_dph=zero
        unitary_tr=zero
        eigf3d=zero
        saddle_warn=.false.
        start_eigf3d_pos=.true.
        m_avg=zero
        m_avg_frohlich=zero
        m_cart=zero
        deigf3d_dth=zero
        deigf3d_dph=zero
        unit_speed=zero
        transport_tensor=zero
        transport_tensor_eig=zero
        transport_eqv_m=zero
        transport_eqv_eigval=zero
        transport_eqv_eigvec=zero

        !Hack to print f(theta,phi) & weights to a file
        if(print_fsph) then
          write(msg,*) degl+1
          filename='f_band_'//TRIM(ADJUSTL(msg))//'-'
          write(msg,*) degl+deg_dim
          filename=TRIM(filename)//TRIM(ADJUSTL(msg))//'.dat'
          io_unit = get_unit()
          open(io_unit,file=TRIM(filename),status='replace')
          write(io_unit,*) 'ntheta=',ntheta,', nphi=',nphi
          write(io_unit,*) 'itheta, iphi, weight, f_n(theta,phi)'
        end if

        do itheta=1,ntheta
          costh=gq_points_costh(itheta) ; sinth=gq_points_sinth(itheta)
          do iphi=1,nphi
            cosph=gq_points_cosph(iphi) ; sinph=gq_points_sinph(iphi)
            weight=gq_weights_th(itheta)*gq_weights_ph(iphi)

            unit_r(1)=sinth*cosph
            unit_r(2)=sinth*sinph
            unit_r(3)=costh

            dr_dth(1)=costh*cosph
            dr_dth(2)=costh*sinph
            dr_dth(3)=-sinth

            dr_dph(1)=-sinph  !sin(theta)*
            dr_dph(2)=cosph   !cos(theta)*
            dr_dph(3)=zero

            do iband=1,deg_dim
              do jband=1,deg_dim
                f3d(iband,jband)=DOT_PRODUCT(unit_r,MATMUL(eig2_diag_cart(:,:,iband,jband),unit_r))
                df3d_dth(iband,jband)=DOT_PRODUCT(dr_dth,MATMUL(eig2_diag_cart(:,:,iband,jband),unit_r))+&
&                DOT_PRODUCT(unit_r,MATMUL(eig2_diag_cart(:,:,iband,jband),dr_dth))
                df3d_dph(iband,jband)=DOT_PRODUCT(dr_dph,MATMUL(eig2_diag_cart(:,:,iband,jband),unit_r))+&
&                DOT_PRODUCT(unit_r,MATMUL(eig2_diag_cart(:,:,iband,jband),dr_dph))
              end do
            end do
            !DIAGONALIZATION
            eigenvec = f3d        !IN
            lwork=-1
            ABI_ALLOCATE(work,(1))
            ABI_ALLOCATE(rwork,(3*deg_dim-2))
            call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
            lwork=int(work(1))
            ABI_DEALLOCATE(work)
            eigenval = zero
            ABI_ALLOCATE(work,(lwork))
            work=zero; rwork=zero
            call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
            ABI_DEALLOCATE(rwork)
            ABI_DEALLOCATE(work)
            unitary_tr = eigenvec !OUT
            eigf3d = eigenval     !OUT
            if(itheta==1 .and. iphi==1) start_eigf3d_pos = eigf3d > 0
            do iband=1,deg_dim
              if(start_eigf3d_pos(iband) .neqv. (eigf3d(iband)>0)) then
                saddle_warn(iband)=.true.
              end if
            end do

            !Hack to print f(theta,phi)
            if(print_fsph) write(io_unit,*) gq_points_th(itheta), gq_points_ph(iphi), weight, eigf3d(:)

            !!DEBUG-Mech.
            !!A=-4.20449; B=0.378191; C=5.309  !Mech's fit
            !A=-4.62503023; B=0.68699088; C=5.20516873 !My fit
            !R = sqrt(B**2 + C**2*sin(theta)**2*(cos(theta)**2 + sin(theta)**2*sin(phi)**2*cos(phi)**2))
            !eigf3d(1) = A - R
            !eigf3d(2) = A + R

            !!angular FD
            !f3dfd(2,iphi,:)=eigf3d(:)

            m_avg = m_avg + weight*sinth*eigf3d
            m_avg_frohlich = m_avg_frohlich + weight*sinth/(abs(eigf3d))**half

            prodc=MATMUL_(f3d,unitary_tr,deg_dim,deg_dim) ; f3d=MATMUL_(unitary_tr,prodc,deg_dim,deg_dim,transa='c')
            !f3d = MATMUL(CONJG(TRANSPOSE(unitary_tr)),MATMUL(f3d,unitary_tr))
            do iband=1,deg_dim
              eigf3d(iband) = real(f3d(iband,iband),dp)
            end do
            prodc=MATMUL_(df3d_dth,unitary_tr,deg_dim,deg_dim) ; df3d_dth=MATMUL_(unitary_tr,prodc,deg_dim,deg_dim,transa='c')
            !df3d_dth=MATMUL(CONJG(TRANSPOSE(unitary_tr)),MATMUL(df3d_dth,unitary_tr))
            do iband=1,deg_dim
              deigf3d_dth(iband) = real(df3d_dth(iband,iband),dp)
            end do
            prodc=MATMUL_(df3d_dph,unitary_tr,deg_dim,deg_dim) ; df3d_dph=MATMUL_(unitary_tr,prodc,deg_dim,deg_dim,transa='c')
            !df3d_dph = MATMUL(CONJG(TRANSPOSE(unitary_tr)),MATMUL(df3d_dph,unitary_tr))
            do iband=1,deg_dim
              deigf3d_dph(iband) = real(df3d_dph(iband,iband),dp)
            end do

            !!DEBUG-Mech.
            !eigf3d(1) = A - R
            !eigf3d(2) = A + R
            !deigf3d_dth(1) = -1./2./R*C**2*(2.*sin(theta)*cos(theta)*(cos(theta)**2 + sin(theta)**2*sin(phi)**2*cos(phi)**2) + 2.*sin(theta)**3*cos(theta)*(sin(phi)**2*cos(phi)**2 - 1))
            !deigf3d_dth(2) =  1./2./R*C**2*(2.*sin(theta)*cos(theta)*(cos(theta)**2 + sin(theta)**2*sin(phi)**2*cos(phi)**2) + 2.*sin(theta)**3*cos(theta)*(sin(phi)**2*cos(phi)**2 - 1))
            !deigf3d_dph(1) = -1./2./R*C**2*sin(theta)**3*(2.*sin(phi)*cos(phi)**3 - 2.*sin(phi)**3*cos(phi))
            !deigf3d_dph(2) =  1./2./R*C**2*sin(theta)**3*(2.*sin(phi)*cos(phi)**3 - 2.*sin(phi)**3*cos(phi))

            !!angular FD
            !if(iphi/=1 .and. itheta/=1) then
            !  deigf3d_dph(:) = (f3dfd(2,iphi,:)-f3dfd(2,iphi-1,:))/two_pi*nphi/sin(theta)
            !else
            !  deigf3d_dph(:) = zero
            !end if
            !if(itheta/=1) then
            !  deigf3d_dth(:) = (f3dfd(2,iphi,:)-f3dfd(1,iphi,:))/pi*ntheta
            !else
            !  deigf3d_dth(:) = zero
            !end if

            unit_speed(1,:) = 2._dp*sinth*cosph*eigf3d + costh*cosph*deigf3d_dth - sinph*deigf3d_dph      !/sin(theta)
            unit_speed(2,:) = 2._dp*sinth*sinph*eigf3d + costh*sinph*deigf3d_dth + cosph*deigf3d_dph      !/sin(theta)
            unit_speed(3,:) = 2._dp*costh*eigf3d - sinth*deigf3d_dth

            do jdeg=1,deg_dim
              do bdir=1,mdim
                do adir=1,mdim
                  transport_tensor(adir,bdir,jdeg) = transport_tensor(adir,bdir,jdeg) + &
&                weight*sinth*unit_speed(adir,jdeg)*unit_speed(bdir,jdeg)/(ABS(eigf3d(jdeg))**2.5_dp)
                end do
              end do
            end do
          end do !iphi
          !!angular FD
          !f3dfd(1,:,:) = f3dfd(2,:,:)
        end do !itheta

        !Hack to print f(theta,phi)
        if(print_fsph) close(io_unit)

        m_avg = quarter/pi*m_avg
        m_avg = one/m_avg

        m_avg_frohlich = quarter/pi*m_avg_frohlich
        m_avg_frohlich = m_avg_frohlich**2

        transport_tensor = 1.0_dp/2.0_dp*transport_tensor

        !Effective masses along directions.
        do adir=1,ndirs
          do iband=1,deg_dim
            do jband=1,deg_dim
              f3d(iband,jband) = dot_product(dirs(:,adir),matmul(eig2_diag_cart(:,:,iband,jband),dirs(:,adir)))
            end do
          end do
          !f3d(:,:) = eig2_diag_cart(adir,adir,:,:)
          eigenvec = f3d        !IN
          lwork=-1
          ABI_ALLOCATE(work,(1))
          ABI_ALLOCATE(rwork,(3*deg_dim-2))
          call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
          lwork=int(work(1))
          ABI_DEALLOCATE(work)
          eigenval = zero
          ABI_ALLOCATE(work,(lwork))
          work=zero; rwork=zero
          call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
          ABI_DEALLOCATE(rwork)
          ABI_DEALLOCATE(work)
          unitary_tr = eigenvec !OUT
          eigf3d = eigenval     !OUT
          m_cart(adir,:)=1._dp/eigf3d(:)
        end do

        do iband=1,deg_dim
          !DIAGONALIZATION
          transport_eqv_eigvec(:,:,iband) = transport_tensor(:,:,iband)
          lwork=-1
          ABI_ALLOCATE(rwork,(1))
          call dsyev('V','U',mdim,transport_eqv_eigvec(:,:,iband),mdim,transport_tensor_eig,rwork,lwork,info)
          lwork=int(rwork(1))
          ABI_DEALLOCATE(rwork)
          transport_tensor_eig = zero
          ABI_ALLOCATE(rwork,(lwork))
          rwork=zero
          call dsyev('V','U',mdim,transport_eqv_eigvec(:,:,iband),mdim,transport_tensor_eig,rwork,lwork,info)
          ABI_DEALLOCATE(rwork)
          transport_eqv_eigvec(:,:,iband) = transpose(transport_eqv_eigvec(:,:,iband)) !So that lines contain eigenvectors.

          prodr=MATMUL_(transport_tensor(:,:,iband),transport_eqv_eigvec(:,:,iband),mdim,mdim,transb='t')
          transport_tensor(:,:,iband)=MATMUL_(transport_eqv_eigvec(:,:,iband),prodr,mdim,mdim)
          !transport_tensor(:,:,iband) = MATMUL(transport_eqv_eigvec(:,:,iband), &
          !                              MATMUL(transport_tensor(:,:,iband),TRANSPOSE(transport_eqv_eigvec(:,:,iband))))

          transport_eqv_eigval(1,iband) = transport_tensor_eig(2)*transport_tensor_eig(3)*(3._dp/8._dp/pi)**2
          transport_eqv_eigval(2,iband) = transport_tensor_eig(3)*transport_tensor_eig(1)*(3._dp/8._dp/pi)**2
          transport_eqv_eigval(3,iband) = transport_tensor_eig(1)*transport_tensor_eig(2)*(3._dp/8._dp/pi)**2
          !The transport tensor loses the sign of the effective mass, this restores it.
          transport_eqv_eigval(:,iband) = DSIGN(transport_eqv_eigval(:,iband),m_avg(iband))
          transport_eqv_m(1,1,iband) = transport_eqv_eigval(1,iband)
          transport_eqv_m(2,2,iband) = transport_eqv_eigval(2,iband)
          transport_eqv_m(3,3,iband) = transport_eqv_eigval(3,iband)

          m_avg_frohlich(iband) = DSIGN(m_avg_frohlich(iband),m_avg(iband))

          prodr=MATMUL_(transport_eqv_m(:,:,iband),transport_eqv_eigvec(:,:,iband),mdim,mdim)
          transport_eqv_m(:,:,iband)=MATMUL_(transport_eqv_eigvec(:,:,iband),prodr,mdim,mdim,transa='t')
          !transport_eqv_m(:,:,iband) = MATMUL(TRANSPOSE(transport_eqv_eigvec(:,:,iband)), &
          !                             MATMUL(transport_eqv_m(:,:,iband),transport_eqv_eigvec(:,:,iband)))

        end do

        call print_tr_efmas(std_out,kpt_rbz(:,ikpt),degl+1,deg_dim,mdim,ndirs,dirs,m_cart,rprimd,transport_eqv_m, &
&                        ntheta,m_avg,m_avg_frohlich,saddle_warn,transport_eqv_eigval,transport_eqv_eigvec)
        call print_tr_efmas(ab_out, kpt_rbz(:,ikpt),degl+1,deg_dim,mdim,ndirs,dirs,m_cart,rprimd,transport_eqv_m, &
&                        ntheta,m_avg,m_avg_frohlich,saddle_warn,transport_eqv_eigval,transport_eqv_eigvec)

        ABI_DEALLOCATE(unit_r)
        ABI_DEALLOCATE(dr_dth)
        ABI_DEALLOCATE(dr_dph)
        ABI_DEALLOCATE(f3d)
        ABI_DEALLOCATE(df3d_dth)
        ABI_DEALLOCATE(df3d_dph)
        ABI_DEALLOCATE(unitary_tr)
        ABI_DEALLOCATE(eigf3d)
        ABI_DEALLOCATE(saddle_warn)
        ABI_DEALLOCATE(start_eigf3d_pos)
        ABI_DEALLOCATE(m_avg)
        ABI_DEALLOCATE(m_avg_frohlich)
        ABI_DEALLOCATE(m_cart)
        ABI_DEALLOCATE(deigf3d_dth)
        ABI_DEALLOCATE(deigf3d_dph)
        ABI_DEALLOCATE(unit_speed)
        ABI_DEALLOCATE(transport_tensor)
        ABI_DEALLOCATE(transport_tensor_eig)
        ABI_DEALLOCATE(transport_eqv_m)
        ABI_DEALLOCATE(transport_eqv_eigval)
        ABI_DEALLOCATE(transport_eqv_eigvec)
        ABI_DEALLOCATE(prodc)
        ABI_DEALLOCATE(prodr)
        !ABI_DEALLOCATE(f3dfd)

      elseif (degenerate .and. mdim==2) then

        ABI_ALLOCATE(unit_r,(mdim))
        ABI_ALLOCATE(dr_dph,(mdim))
        ABI_ALLOCATE(f3d,(deg_dim,deg_dim))
        ABI_ALLOCATE(df3d_dph,(deg_dim,deg_dim))
        ABI_ALLOCATE(unitary_tr,(deg_dim,deg_dim))
        ABI_ALLOCATE(eigf3d,(deg_dim))
        ABI_ALLOCATE(saddle_warn,(deg_dim))
        ABI_ALLOCATE(start_eigf3d_pos,(deg_dim))
        ABI_ALLOCATE(m_avg,(deg_dim))
        ABI_ALLOCATE(m_avg_frohlich,(deg_dim))
        ABI_ALLOCATE(m_cart,(ndirs,deg_dim))
        ABI_ALLOCATE(deigf3d_dph,(deg_dim))
        ABI_ALLOCATE(unit_speed,(mdim,deg_dim))
        ABI_ALLOCATE(transport_tensor,(mdim,mdim,deg_dim))
        ABI_ALLOCATE(cart_rotation,(mdim,mdim))
        ABI_ALLOCATE(transport_tensor_eig,(mdim))
        ABI_ALLOCATE(transport_eqv_m,(mdim,mdim,deg_dim))
        ABI_ALLOCATE(transport_eqv_eigval,(mdim,deg_dim))
        ABI_ALLOCATE(transport_eqv_eigvec,(mdim,mdim,deg_dim))
        ABI_ALLOCATE(transport_tensor_scale,(deg_dim))
        ABI_ALLOCATE(prodc,(deg_dim,deg_dim))
        ABI_ALLOCATE(prodr,(mdim,mdim))
        unit_r=zero
        dr_dph=zero
        f3d=zero
        df3d_dph=zero
        unitary_tr=zero
        eigf3d=zero
        saddle_warn=.false.
        start_eigf3d_pos=.true.
        m_avg=zero
        m_avg_frohlich=zero
        m_cart=zero
        deigf3d_dph=zero
        unit_speed=zero
        transport_tensor=zero
        cart_rotation=zero
        transport_tensor_eig=zero
        transport_eqv_m=zero
        transport_eqv_eigval=zero
        transport_eqv_eigvec=zero
        transport_tensor_scale=zero

        do iphi=1,nphi
          cosph=gq_points_cosph(iphi) ; sinph=gq_points_sinph(iphi)
          weight=gq_weights_ph(iphi)

          unit_r(1)=cosph
          unit_r(2)=sinph

          dr_dph(1)=-sinph
          dr_dph(2)=cosph

          do iband=1,deg_dim
            do jband=1,deg_dim
              matr2d = eig2_diag_cart(1:mdim,1:mdim,iband,jband)
              f3d(iband,jband)=DOT_PRODUCT(unit_r,MATMUL(matr2d,unit_r))
              df3d_dph(iband,jband)=DOT_PRODUCT(dr_dph,MATMUL(matr2d,unit_r))+&
&              DOT_PRODUCT(unit_r,MATMUL(matr2d,dr_dph))
            end do
          end do

          !DIAGONALIZATION
          eigenvec = f3d        !IN
          lwork=-1
          ABI_ALLOCATE(work,(1))
          ABI_ALLOCATE(rwork,(3*deg_dim-2))
          call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
          lwork=int(work(1))
          ABI_DEALLOCATE(work)
          eigenval = zero
          ABI_ALLOCATE(work,(lwork))
          work=zero; rwork=zero
          call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
          ABI_DEALLOCATE(rwork)
          ABI_DEALLOCATE(work)
          unitary_tr = eigenvec !OUT
          eigf3d = eigenval     !OUT
          if(iphi==1) start_eigf3d_pos = eigf3d > 0
          do iband=1,deg_dim
            if(start_eigf3d_pos(iband) .neqv. (eigf3d(iband)>0)) then
              saddle_warn(iband)=.true.
            end if
          end do

          m_avg = m_avg + weight*eigf3d
          m_avg_frohlich = m_avg_frohlich + weight/(abs(eigf3d))**half

          prodc=MATMUL_(f3d,unitary_tr,deg_dim,deg_dim) ; f3d=MATMUL_(unitary_tr,prodc,deg_dim,deg_dim,transa='c')
          !f3d = MATMUL(CONJG(TRANSPOSE(unitary_tr)),MATMUL(f3d,unitary_tr))
          do iband=1,deg_dim
            eigf3d(iband) = real(f3d(iband,iband),dp)
          end do
          prodc=MATMUL_(df3d_dph,unitary_tr,deg_dim,deg_dim) ; df3d_dph=MATMUL_(unitary_tr,prodc,deg_dim,deg_dim,transa='c')
          !df3d_dph = MATMUL(CONJG(TRANSPOSE(unitary_tr)),MATMUL(df3d_dph,unitary_tr))
          do iband=1,deg_dim
            deigf3d_dph(iband) = real(df3d_dph(iband,iband),dp)
          end do

          unit_speed(1,:) = 2._dp*cosph*eigf3d - sinph*deigf3d_dph
          unit_speed(2,:) = 2._dp*sinph*eigf3d + cosph*deigf3d_dph

          do jdeg=1,deg_dim
            do bdir=1,mdim
              do adir=1,mdim
                transport_tensor(adir,bdir,jdeg) = transport_tensor(adir,bdir,jdeg) + &
&                weight*unit_speed(adir,jdeg)*unit_speed(bdir,jdeg)/(ABS(eigf3d(jdeg))**2)
              end do
            end do
          end do

        end do !iphi

        !!!DEBUG

        m_avg = half/pi*m_avg
        m_avg = one/m_avg

        m_avg_frohlich = half/pi*m_avg_frohlich
        m_avg_frohlich = m_avg_frohlich**2

        transport_tensor = 1.0_dp/2.0_dp*transport_tensor

        !Effective masses along directions.
        do adir=1,ndirs
          do iband=1,deg_dim
            do jband=1,deg_dim
              f3d(iband,jband) = dot_product(dirs(:,adir),matmul(eig2_diag_cart(:,:,iband,jband),dirs(:,adir)))
            end do
          end do
          !f3d(:,:) = eig2_diag_cart(adir,adir,:,:)
          eigenvec = f3d        !IN
          lwork=-1
          ABI_ALLOCATE(work,(1))
          ABI_ALLOCATE(rwork,(3*deg_dim-2))
          call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
          lwork=int(work(1))
          ABI_DEALLOCATE(work)
          eigenval = zero
          ABI_ALLOCATE(work,(lwork))
          work=zero; rwork=zero
          call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
          ABI_DEALLOCATE(rwork)
          ABI_DEALLOCATE(work)
          unitary_tr = eigenvec !OUT
          eigf3d = eigenval     !OUT
          m_cart(adir,:)=1._dp/eigf3d(:)
        end do

        do iband=1,deg_dim
            !DIAGONALIZATION
          cart_rotation = transport_tensor(:,:,iband)
          lwork=-1
          ABI_ALLOCATE(rwork,(1))
          call dsyev('V','U',mdim,cart_rotation,mdim,transport_tensor_eig,rwork,lwork,info)
          lwork=int(rwork(1))
          ABI_DEALLOCATE(rwork)
          transport_tensor_eig = zero
          ABI_ALLOCATE(rwork,(lwork))
          rwork=zero
          call dsyev('V','U',mdim,cart_rotation,mdim,transport_tensor_eig,rwork,lwork,info)
          ABI_DEALLOCATE(rwork)
          transport_eqv_eigvec(:,:,iband) = transpose(cart_rotation(:,:)) !So that lines contain eigenvectors, not columns.

          prodr=MATMUL_(transport_tensor(:,:,iband),cart_rotation,mdim,mdim)
          transport_tensor(:,:,iband)=MATMUL_(cart_rotation,prodr,mdim,mdim,transa='t')
          !transport_tensor(:,:,iband) = MATMUL(TRANSPOSE(cart_rotation),MATMUL(transport_tensor(:,:,iband),cart_rotation))

          transport_eqv_eigval(1,iband) = 0.5*m_avg(iband)*(1.0 + transport_tensor_eig(2)/transport_tensor_eig(1))
          transport_eqv_eigval(2,iband) = transport_eqv_eigval(1,iband)*transport_tensor_eig(1)/transport_tensor_eig(2)
          !The transport tensor loses the sign of the effective mass, this restores it.
          transport_eqv_eigval(:,iband) = SIGN(transport_eqv_eigval(:,iband),m_avg(iband))
          transport_eqv_m(1,1,iband) = transport_eqv_eigval(1,iband)
          transport_eqv_m(2,2,iband) = transport_eqv_eigval(2,iband)
          transport_tensor_scale(iband) = sqrt(transport_tensor_eig(1)*transport_tensor_eig(2))/two_pi

          m_avg_frohlich(iband) = SIGN(m_avg_frohlich(iband),m_avg(iband))

          prodr=MATMUL_(transport_eqv_m(:,:,iband),cart_rotation,mdim,mdim,transb='t')
          transport_eqv_m(:,:,iband)=MATMUL_(cart_rotation,prodr,mdim,mdim)
          !transport_eqv_m(:,:,iband) = MATMUL(cart_rotation,MATMUL(transport_eqv_m(:,:,iband),TRANSPOSE(cart_rotation)))

        end do

        call print_tr_efmas(std_out,kpt_rbz(:,ikpt),degl+1,deg_dim,mdim,ndirs,dirs,m_cart,rprimd,transport_eqv_m, &
&                        ntheta,m_avg,m_avg_frohlich,saddle_warn,transport_eqv_eigval,transport_eqv_eigvec,transport_tensor_scale)
        call print_tr_efmas(ab_out, kpt_rbz(:,ikpt),degl+1,deg_dim,mdim,ndirs,dirs,m_cart,rprimd,transport_eqv_m, &
&                        ntheta,m_avg,m_avg_frohlich,saddle_warn,transport_eqv_eigval,transport_eqv_eigvec,transport_tensor_scale)

        ABI_DEALLOCATE(unit_r)
        ABI_DEALLOCATE(dr_dph)
        ABI_DEALLOCATE(f3d)
        ABI_DEALLOCATE(df3d_dph)
        ABI_DEALLOCATE(unitary_tr)
        ABI_DEALLOCATE(eigf3d)
        ABI_DEALLOCATE(saddle_warn)
        ABI_DEALLOCATE(start_eigf3d_pos)
        ABI_DEALLOCATE(m_avg)
        ABI_DEALLOCATE(m_avg_frohlich)
        ABI_DEALLOCATE(m_cart)
        ABI_DEALLOCATE(deigf3d_dph)
        ABI_DEALLOCATE(unit_speed)
        ABI_DEALLOCATE(transport_tensor)
        ABI_DEALLOCATE(cart_rotation)
        ABI_DEALLOCATE(transport_tensor_eig)
        ABI_DEALLOCATE(transport_eqv_m)
        ABI_DEALLOCATE(transport_eqv_eigval)
        ABI_DEALLOCATE(transport_eqv_eigvec)
        ABI_DEALLOCATE(transport_tensor_scale)
        ABI_DEALLOCATE(prodc)
        ABI_DEALLOCATE(prodr)

      elseif (degenerate .and. mdim==1) then

        ABI_ALLOCATE(f3d,(deg_dim,deg_dim))
        ABI_ALLOCATE(unitary_tr,(deg_dim,deg_dim))
        ABI_ALLOCATE(eigf3d,(deg_dim))
        ABI_ALLOCATE(m_cart,(ndirs,deg_dim))
        ABI_ALLOCATE(transport_eqv_m,(mdim,mdim,deg_dim))
        ABI_ALLOCATE(m_avg,(deg_dim))
        ABI_ALLOCATE(m_avg_frohlich,(deg_dim))
        ABI_ALLOCATE(saddle_warn,(deg_dim))


        f3d=zero
        unitary_tr=zero
        eigf3d=zero
        m_cart=zero
        transport_eqv_m=zero
        m_avg=zero
        m_avg_frohlich=zero
        saddle_warn=.false.

        f3d(:,:) = eig2_diag_cart(1,1,:,:)

        !DIAGONALIZATION
        eigenvec = f3d        !IN
        lwork=-1
        ABI_ALLOCATE(work,(1))
        ABI_ALLOCATE(rwork,(3*deg_dim-2))
        call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
        lwork=int(work(1))
        ABI_DEALLOCATE(work)
        eigenval = zero
        ABI_ALLOCATE(work,(lwork))
        work=zero; rwork=zero
        call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
        ABI_DEALLOCATE(rwork)
        ABI_DEALLOCATE(work)
        unitary_tr = eigenvec !OUT
        eigf3d = eigenval     !OUT

        transport_eqv_m(1,1,:)=1._dp/eigf3d(:)

        !Effective masses along directions.
        do adir=1,ndirs
          do iband=1,deg_dim
            do jband=1,deg_dim
              f3d(iband,jband) = dot_product(dirs(:,adir),matmul(eig2_diag_cart(:,:,iband,jband),dirs(:,adir)))
            end do
          end do
          eigenvec = f3d        !IN
          lwork=-1
          ABI_ALLOCATE(work,(1))
          ABI_ALLOCATE(rwork,(3*deg_dim-2))
          call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
          lwork=int(work(1))
          ABI_DEALLOCATE(work)
          eigenval = zero
          ABI_ALLOCATE(work,(lwork))
          work=zero; rwork=zero
          call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
          ABI_DEALLOCATE(rwork)
          ABI_DEALLOCATE(work)
          eigf3d = eigenval     !OUT
          m_cart(adir,:)=1._dp/eigf3d(:)
        end do

        call print_tr_efmas(std_out,kpt_rbz(:,ikpt),degl+1,deg_dim,mdim,ndirs,dirs,m_cart,rprimd,transport_eqv_m,&
&          ntheta,m_avg,m_avg_frohlich,saddle_warn)
        call print_tr_efmas(ab_out, kpt_rbz(:,ikpt),degl+1,deg_dim,mdim,ndirs,dirs,m_cart,rprimd,transport_eqv_m,&
&          ntheta,m_avg,m_avg_frohlich,saddle_warn)

        ABI_DEALLOCATE(f3d)
        ABI_DEALLOCATE(unitary_tr)
        ABI_DEALLOCATE(eigf3d)
        ABI_DEALLOCATE(m_cart)
        ABI_DEALLOCATE(transport_eqv_m)
        ABI_DEALLOCATE(m_avg)
        ABI_DEALLOCATE(m_avg_frohlich)
        ABI_DEALLOCATE(saddle_warn)


      end if !(degenerate)

      ABI_DEALLOCATE(eig2_diag_cart)

      ABI_DEALLOCATE(eigenval)

      ABI_DEALLOCATE(eigenvec)
    end do !ideg

  end do !ikpt

  ABI_DEALLOCATE(eff_mass)
  ABI_DEALLOCATE(dirs)
  ABI_DEALLOCATE(gq_points_th)
  ABI_DEALLOCATE(gq_points_costh)
  ABI_DEALLOCATE(gq_points_sinth)
  ABI_DEALLOCATE(gq_weights_th)
  ABI_DEALLOCATE(gq_points_ph)
  ABI_DEALLOCATE(gq_points_cosph)
  ABI_DEALLOCATE(gq_points_sinph)
  ABI_DEALLOCATE(gq_weights_ph)

  write(std_out,'(3a)') ch10,' END OF EFFECTIVE MASSES SECTION',ch10
  write(ab_out, '(3a)') ch10,' END OF EFFECTIVE MASSES SECTION',ch10

 end subroutine efmas_analysis
!!***

!----------------------------------------------------------------------

!!****f* m_efmas/MATMUL_DP
!! NAME
!! MATMUL_DP
!!
!! FUNCTION
!! Mimic MATMUL Fortran intrinsic function with BLAS3: C=A.B
!! This is a temporary workaround to make tests pass on intel/mkl architectures
!! Real version
!!
!! INPUTS
!!  aa(:,:),bb(:,:)= input matrices
!!  mm,nn= sizes of output matrix
!!  [transa,transb]= equivalent to transa, transb args of gemm ('n','t','c')
!!                   if not present, default is 'n'.
!!
!! OUTPUT
!!  MATMUL_DP(mm,nn)= output matrix A.B
!!
!! SOURCE

function MATMUL_DP(aa,bb,mm,nn,transa,transb)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mm,nn
 character(len=1),optional,intent(in) :: transa,transb
!arrays
 real(dp),intent(in) :: aa(:,:),bb(:,:)
 real(dp) :: MATMUL_DP(mm,nn)

!Local variables-------------------------------
 integer :: kk,lda,ldb
 character(len=1) :: transa_,transb_

! *************************************************************************

 transa_='n';if (present(transa)) transa_=transa
 transb_='n';if (present(transb)) transb_=transb

 lda=size(aa,1) ; ldb=size(bb,1)

 if (transa_=='n') then
   kk=size(aa,2)
   if (size(aa,1)/=mm) then
     MSG_BUG('Error in sizes!')
   end if
 else
   kk=size(aa,1)
   if (size(aa,2)/=mm) then
     MSG_BUG('Error in sizes!')
   end if
 end if

 if (transb_=='n') then
   if (size(bb,1)/=kk.or.size(bb,2)/=nn) then
     MSG_BUG('Error in sizes!')
   end if
 else
   if (size(bb,1)/=nn.or.size(bb,2)/=kk) then
     MSG_BUG('Error in sizes!')
   end if
 end if

 call DGEMM(transa_,transb_,mm,nn,kk,one,aa,lda,bb,ldb,zero,MATMUL_DP,mm)

end function MATMUL_DP
!!***

!----------------------------------------------------------------------

!!****f* m_efmas/MATMUL_DPC
!! NAME
!! MATMUL_DPC
!!
!! FUNCTION
!! Mimic MATMUL Fortran intrinsic function with BLAS3
!! This is a temporary workaround to make tests pass on intel/mkl architectures
!! Complex version
!!
!! INPUTS
!!  aa(:,:),bb(:,:)= input matrices
!!  mm,nn= sizes of output matrix
!!  [transa,transb]= equivalent to transa, transb args of gemm ('n','t','c')
!!                   if not present, default is 'n'.
!!
!! OUTPUT
!!  MATMUL_DPC(:,:)= output matrix A.B
!!
!! SOURCE

function MATMUL_DPC(aa,bb,mm,nn,transa,transb)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mm,nn
 character(len=1),optional,intent(in) :: transa,transb
!arrays
 complex(dpc),intent(in) :: aa(:,:),bb(:,:)
 complex(dpc) :: MATMUL_DPC(mm,nn)

!Local variables-------------------------------
 integer :: kk,lda,ldb
 character(len=1) :: transa_,transb_

! *************************************************************************

 transa_='n';if (present(transa)) transa_=transa
 transb_='n';if (present(transb)) transb_=transb

 lda=size(aa,1) ; ldb=size(bb,1)

 if (transa_=='n') then
   kk=size(aa,2)
   if (size(aa,1)/=mm) then
     MSG_BUG('Error in sizes!')
   end if
 else
   kk=size(aa,1)
   if (size(aa,2)/=mm) then
     MSG_BUG('Error in sizes!')
   end if
 end if

 if (transb_=='n') then
   if (size(bb,1)/=kk.or.size(bb,2)/=nn) then
     MSG_BUG('Error in sizes!')
   end if
 else
   if (size(bb,1)/=nn.or.size(bb,2)/=kk) then
     MSG_BUG('Error in sizes!')
   end if
 end if

 call ZGEMM(transa_,transb_,mm,nn,kk,cone,aa,lda,bb,ldb,czero,MATMUL_DPC,mm)

end function MATMUL_DPC
!!***

end module m_efmas
!!***
