!!****m* ABINIT/m_pawang
!! NAME
!!  m_pawang
!!
!! FUNCTION
!!  This module contains the definition of the pawang_type structured datatype,
!!  as well as related functions and methods.
!!  pawang_type variables define ANGular mesh discretization of PAW augmentation
!!  regions and related data.
!!
!! COPYRIGHT
!! Copyright (C) 2013-2024 ABINIT group (MT, FJ, BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

MODULE m_pawang

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 use m_paw_sphharm, only : initylmr, ylm_angular_mesh, mat_mlms2jmj, mat_slm2ylm, &
&                          lsylm, realgaunt, nablarealgaunt

 implicit none

 private

!Public procedures.
 public :: pawang_init       ! Constructor
 public :: pawang_free       ! Free memory

!!***

!----------------------------------------------------------------------


!!****t* m_pawang/pawang_type
!! NAME
!! pawang_type
!!
!! FUNCTION
!! For PAW: ANGular mesh discretization of PAW augmentation regions and related data.
!!
!! SOURCE

 type,public :: pawang_type

!Integer scalars

  integer :: angl_size=0
   ! Dimension of paw angular mesh
   ! angl_size=ntheta*nphi

  integer :: l_max=-1
   ! Maximum value of angular momentum l+1

  integer :: l_size_max=-1
   ! Maximum value of angular momentum +1
   ! leading to non-zero Gaunt coefficients
   ! l_size_max = 2*l_max-1

  integer :: ngnt=0
   ! Number of non zero Gaunt coefficients

  integer :: nnablagnt=0
   ! Number of non zero Gaunt coefficient derivatives

  integer :: ntheta, nphi
   ! Dimensions of paw angular mesh

  integer :: nsym
   ! Number of symmetry elements in space group

  integer :: nabgnt_option=-1
   ! Option for nablarealgaunt coefficients:
   ! nabgnt_option==0, nablarealgaunt coeffs are not computed (and not allocated)
   ! nabgnt_option==1, nablarealgaunt coeffs are computed up to l_max

  integer :: gnt_option=-1
   ! Option for Gaunt coefficients:
   !  gnt_option==0, Gaunt coeffs are not computed (and not allocated)
   !  gnt_option==1, Gaunt coeffs are computed up to 2*l_size_max-1
   !  gnt_option==2, Gaunt coeffs are computed up to l_size_max

  integer :: use_ls_ylm=0
   ! Flag: use_ls_ylm=1 if pawang%ls_ylm is allocated

  integer :: ylm_size=0
   ! Size of ylmr/ylmrgr arrays

!Integer arrays

  integer, allocatable :: gntselect(:,:)
   ! gntselect(l_size_max**2,l_max**2*(l_max**2+1)/2)
   ! Selection rules for Gaunt coefficients stored as (LM,ij) where ij is in packed form.
   ! (if gntselect>0, Gaunt coeff. is non-zero)

  integer, allocatable :: nablagntselect(:,:,:)
  ! nablagntselect(l_size_max**2,l_max**2,l_max**2)
  ! Selection rules for nablaGaunt coefficients
  ! (if nablagntselect>0, nablGaunt coeff. is non-zero)

!Real (real(dp)) arrays

  real(dp), allocatable :: anginit(:,:)
   ! anginit(3,angl_size)
   ! For each point of the angular mesh, gives the coordinates
   ! of the corresponding point on an unitary sphere
   ! Not used in present version (5.3)

  real(dp), allocatable :: angwgth(:)
   ! angwgth(angl_size)
   ! For each point of the angular mesh, gives the weight
   ! of the corresponding point on an unitary sphere

  real(dp), allocatable :: ls_ylm(:,:,:)
   ! ls_ylm(2,l_max**2*(l_max**2+1)/2,2)
   ! LS operator in the real spherical harmonics basis
   ! ls_ylm(ilm1m2,ispin)= <sigma, y_lm1| LS |y_lm2, sigma_prime>

  real(dp), allocatable :: nablarealgnt(:)
   ! nablarealgnt(2,nnablagnt)
   ! Non zero real nablaGaunt coefficients

  real(dp), allocatable :: realgnt(:)
   ! realgnt(ngnt)
   ! Non zero real Gaunt coefficients

  real(dp), allocatable :: ylmr(:,:)
   ! ylmr(ylm_size,angl_size)
   ! Real Ylm calculated in real space

  real(dp), allocatable :: ylmrgr(:,:,:)
   ! ylmrgr(1:3,ylm_size,angl_size)
   ! First gradients of real Ylm calculated in real space (cart. coordinates)

  real(dp), allocatable :: zarot(:,:,:,:)
   !  zarot(l_size_max,l_size_max,l_max,nsym)
   !  Coeffs of the transformation of real spherical
   !  harmonics under the symmetry operations symrec.

 end type pawang_type
!!***

CONTAINS

!===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_pawang/pawang_init
!! NAME
!! pawang_init
!!
!! FUNCTION
!!  Initialize a pawang datastructure
!!
!! INPUTS
!!  gnt_option=flag activated if pawang%gntselect and pawang%realgnt have to be allocated
!!             also determine the size of these pointers
!!  nabgnt_option=flag activated if pawang%nablagntselect and pawang%nablarealgnt have to be allocated
!!  lmax=maximum value of angular momentum l
!!  nphi,ntheta=dimensions of paw angular mesh
!!  nsym=number of symetries
!!  ngrad2_ylm=order of spherical harmonics gradients to be stored (0, 1 or 2)
!!  use_angular_grid=flag activated if angular grid data have to be allocated
!!                   (pawang%angwgth, pawang%anginit)
!!  use_ylm=flag activated if spherical harmonics have to be allocated and computed
!!          (pawang%ylmr, pawang%ylmrgr)
!!  use_ls_ylm=flag activated if LS operator has to be allocated and computed
!!          (pawang%ls_ylm)
!!
!! OUTPUT
!!  Pawang <type(pawang_type)>=ANGular mesh discretization and related data
!!
!! SOURCE

subroutine pawang_init(Pawang,gnt_option,nabgnt_option,lmax,nphi,ntheta,nsym,ngrad2_ylm,&
&                      use_angular_grid,use_ylm,use_ls_ylm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: gnt_option,nabgnt_option,lmax,nphi,nsym,ntheta,ngrad2_ylm
 integer,intent(in) :: use_angular_grid,use_ylm,use_ls_ylm
 type(Pawang_type),intent(inout) :: Pawang

!Local variables-------------------------------
!scalars
 integer :: ll,sz1,sz2
!arrays
 real(dp),allocatable :: rgnt_tmp(:)
 real(dp),allocatable :: nablargnt_tmp(:)

! *************************************************************************

 !@Pawang_type

 Pawang%l_max=lmax+1
 Pawang%l_size_max=2*Pawang%l_max-1
 Pawang%nsym=nsym

 if (use_angular_grid==1) then
   Pawang%nphi=nphi
   Pawang%ntheta=ntheta
   call ylm_angular_mesh(Pawang%ntheta,Pawang%nphi,Pawang%angl_size,Pawang%anginit,Pawang%angwgth)
 else
   Pawang%nphi=0
   Pawang%ntheta=0
   Pawang%angl_size=0
 end if

 if (use_ylm>0.and.Pawang%angl_size>0) then
   ll=Pawang%l_size_max+1
   Pawang%ylm_size=ll**2
   LIBPAW_ALLOCATE(Pawang%ylmr,(Pawang%ylm_size,Pawang%angl_size))
   if (ngrad2_ylm==2) then
     LIBPAW_ALLOCATE(pawang%ylmrgr,(9,Pawang%ylm_size,Pawang%angl_size))
     call initylmr(ll,0,pawang%angl_size,pawang%angwgth,3,pawang%anginit,pawang%ylmr,&
&                  ylmr_gr=pawang%ylmrgr)
   else if (ngrad2_ylm==1) then
       LIBPAW_ALLOCATE(Pawang%ylmrgr,(3,Pawang%ylm_size,Pawang%angl_size))
       call initylmr(ll,0,pawang%angl_size,pawang%angwgth,2,pawang%anginit,pawang%ylmr,&
&                    ylmr_gr=pawang%ylmrgr)
   else
     call initylmr(ll,0,pawang%angl_size,pawang%angwgth,1,pawang%anginit,pawang%ylmr)
   end if
 else
   Pawang%ylm_size=0
 end if

 Pawang%gnt_option=gnt_option
 if (Pawang%gnt_option==1.or.Pawang%gnt_option==2) then
   if (Pawang%gnt_option==1) then
     sz1=(Pawang%l_size_max)**2
     sz2=((Pawang%l_max**2)*(Pawang%l_max**2+1))/2
     LIBPAW_ALLOCATE(rgnt_tmp,(sz1*sz2))
     LIBPAW_ALLOCATE(pawang%gntselect,(sz1,sz2))
     call realgaunt(Pawang%l_max,Pawang%ngnt,Pawang%gntselect,rgnt_tmp)
   else if (Pawang%gnt_option==2) then
     sz1=(2*Pawang%l_size_max-1)**2
     sz2=((Pawang%l_size_max)**2*(Pawang%l_size_max**2+1))/2
     LIBPAW_ALLOCATE(rgnt_tmp,(sz1*sz2))
     LIBPAW_ALLOCATE(pawang%gntselect,(sz1,sz2))
     call realgaunt(Pawang%l_size_max,Pawang%ngnt,Pawang%gntselect,rgnt_tmp)
   end if
   if (allocated(pawang%realgnt))  then
     LIBPAW_DEALLOCATE(pawang%realgnt)
   end if
   LIBPAW_ALLOCATE(Pawang%realgnt,(Pawang%ngnt))
   Pawang%realgnt(1:Pawang%ngnt)=rgnt_tmp(1:Pawang%ngnt)
   LIBPAW_DEALLOCATE(rgnt_tmp)
 end if

 Pawang%nabgnt_option=nabgnt_option
 if (Pawang%nabgnt_option==1) then
   sz1=(Pawang%l_size_max)**2
   sz2=(Pawang%l_max)**2
   LIBPAW_ALLOCATE(nablargnt_tmp,(sz1*sz2*sz2))
   LIBPAW_ALLOCATE(pawang%nablagntselect,(sz1,sz2,sz2))
   call nablarealgaunt(pawang%l_size_max,pawang%l_max, &
&                      pawang%nnablagnt,pawang%nablagntselect,nablargnt_tmp)
   if (allocated(pawang%nablarealgnt)) then
     LIBPAW_DEALLOCATE(pawang%nablarealgnt)
   end if
   LIBPAW_ALLOCATE(pawang%nablarealgnt,(pawang%nnablagnt))
   Pawang%nablarealgnt(1:Pawang%nnablagnt)=nablargnt_tmp(1:Pawang%nnablagnt)
   LIBPAW_DEALLOCATE(nablargnt_tmp)
 end if

 Pawang%use_ls_ylm=use_ls_ylm
 if (use_ls_ylm>0) then
   LIBPAW_ALLOCATE(pawang%ls_ylm,(2,Pawang%l_max**2*(Pawang%l_max**2+1)/2,2))
   call lsylm(pawang%ls_ylm,lmax)
 end if

 if (nsym>0) then
   if (.not.allocated(Pawang%zarot)) then
     LIBPAW_ALLOCATE(Pawang%zarot,(Pawang%l_size_max,Pawang%l_size_max,Pawang%l_max,nsym))
   end if
 end if

end subroutine pawang_init
!!***

!----------------------------------------------------------------------

!!****f* m_pawang/pawang_free
!! NAME
!! pawang_free
!!
!! FUNCTION
!!  Free all dynamic memory and reset all flags stored in a pawang datastructure
!!
!! SIDE EFFECTS
!!  Pawang <type(pawang_type)>=ANGular mesh discretization and related data
!!
!! SOURCE

subroutine pawang_free(Pawang)

!Arguments ------------------------------------
!scalars
 type(Pawang_type),intent(inout) :: Pawang

! *************************************************************************

 !@Pawang_type
 if (allocated(pawang%angwgth))    then
   LIBPAW_DEALLOCATE(pawang%angwgth)
 end if
 if (allocated(pawang%anginit))    then
   LIBPAW_DEALLOCATE(pawang%anginit)
 end if
 if (allocated(pawang%zarot))      then
   LIBPAW_DEALLOCATE(pawang%zarot)
 end if
 if (allocated(pawang%gntselect))  then
   LIBPAW_DEALLOCATE(pawang%gntselect)
 end if
  if (allocated(pawang%nablagntselect))  then
   LIBPAW_DEALLOCATE(pawang%nablagntselect)
 end if
 if (allocated(pawang%realgnt))    then
   LIBPAW_DEALLOCATE(pawang%realgnt)
 end if
 if (allocated(pawang%nablarealgnt))    then
   LIBPAW_DEALLOCATE(pawang%nablarealgnt)
 end if
 if (allocated(pawang%ylmr))       then
   LIBPAW_DEALLOCATE(pawang%ylmr)
 end if
 if (allocated(pawang%ylmrgr))     then
   LIBPAW_DEALLOCATE(pawang%ylmrgr)
 end if
 if (allocated(pawang%ls_ylm))     then
   LIBPAW_DEALLOCATE(pawang%ls_ylm)
 end if

 pawang%angl_size =0
 pawang%ylm_size  =0
 pawang%use_ls_ylm=0
 pawang%l_max=-1
 pawang%l_size_max=-1
 pawang%gnt_option=-1
 pawang%nabgnt_option=-1
 pawang%ngnt=0

end subroutine pawang_free
!!***

!----------------------------------------------------------------------

END MODULE m_pawang
!!***
