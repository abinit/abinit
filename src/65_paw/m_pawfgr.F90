!!****m* ABINIT/m_pawfgr
!! NAME
!!  m_pawfgr
!!
!! FUNCTION
!!  This module contains the definition of the pawfgr_type structured datatype,
!!  as well as related functions and methods.
!!  pawfgr_type variables define Fine rectangular GRid parameters and related data.
!!
!! COPYRIGHT
!! Copyright (C) 2013-2020 ABINIT group (MT, FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!  * Routines tagged with "@type_name" are strongly connected to the definition of the data type.
!!    Strongly connected means that the proper functioning of the implementation relies on the
!!    assumption that the tagged procedure is consistent with the type declaration.
!!    Every time a developer changes the structure "type_name" adding new entries, he/she has to make sure
!!    that all the strongly connected routines are changed accordingly to accommodate the modification of the data type
!!    Typical examples of strongly connected routines are creation, destruction or reset methods.
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

MODULE m_pawfgr

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_dtset

 use m_kg,       only : getcut

 implicit none

 private

!public procedures.
 public :: pawfgr_init
 public :: pawfgr_destroy
 public :: pawfgr_nullify
 public :: indgrid
!!***

!----------------------------------------------------------------------


!!****t* m_pawfgr/pawfgr_type
!! NAME
!! pawfgr_type
!!
!! FUNCTION
!! For PAW, Fine GRid parameters and related data
!!
!! SOURCE

 type,public :: pawfgr_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

!Integer scalars

  integer :: mgfft, nfft
   ! Values of mffft and nfft for the fine grid:
   !   mgfft= max(ngfft(i)) [max. size of 1D FFT grid]
   !   nfft=ngfft1*ngfft2*ngfft3 [number of pts in the FFT box]

  integer :: mgfftc, nfftc
   ! Values of mffft and nfft for the COARSE grid:
   !   mgfftc= max(ngfftc(i)) [max. size of 1D FFT grid]
   !   nfftc=ngfftc1*ngfftc2*ngfftc3 [number of pts in the FFT box]

  integer :: usefinegrid
   ! Flag: =1 if a double-grid is used to convert spherical data
   !       to Fourier grid. =0 otherwise

!Integer arrays

  ! MGTODO: Replace with allocatable
  integer, pointer :: coatofin(:)
   ! coatofin(nfftc)
   ! Index of the points of the coarse grid on the fine grid

  integer, pointer :: fintocoa(:)
   ! fintocoa(nfft)
   ! Index of the points of the fine grid on the coarse grid
   !  (=0 if the point of the fine grid does not belong to the coarse grid)

  integer :: ngfft(18)
   ! ngfft(1:18)=integer array with FFT box dimensions and other
   ! information on FFTs, for the fine rectangular grid

  integer :: ngfftc(18)
   ! ngfft(1:18)=integer array with FFT box dimensions and other
   ! information on FFTs, for the COARSE rectangular grid

 end type pawfgr_type
!!***

CONTAINS

!===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_pawfgr/pawfgr_init
!! NAME
!! pawfgr_init
!!
!! FUNCTION
!!  Initialize a pawfgr_type datatype, reporting also info on the mesh
!!  according to the method used (norm-conserving PSP or PAW)
!!
!! INPUTS
!!  k0(3)=input k vector for k+G sphere
!!  Dtset <type(dataset_type)>=all input variables for this dataset
!!   %dilatmx
!!   %usepaw
!!   %usewvl
!!   %natom
!!   %ngfft
!!   %ngfftdg
!!   %nfft
!!   %mgfft
!!   %mgfftdg
!!   %dilatmx
!!   %pawecutdg
!!   %ecut
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!
!! OUTPUT
!!  ecut_eff=effective energy cutoff (hartree) for coarse planewave basis sphere
!!  ecutdg_eff=effective energy cutoff (hartree) for dense planewave basis sphere
!!  gsqcutc_eff=(PAW) Fourier cutoff on G^2 for "large sphere" of radius double for the coarse FFT grid
!   gsqcutf_eff=Fourier cutoff on G^2 for "large sphere" of radius double for the dense FFT grid
!!  nfftf=(effective) number of FFT grid points (for this proc), for dense FFT mesh
!!  mgfftf=maximum size of 1D FFTs, for dense FFT mesh
!!  ngfftc(18),ngfftf(18)=contain all needed information about 3D FFT, for coarse and dense FFT mesh, resp.
!!                        see ~abinit/doc/variables/vargs.htm#ngfft
!!  Pawfgr<pawfgr_type>=For PAW, Fine rectangular GRid parameters and related data
!!
!! PARENTS
!!      bethe_salpeter,eph,gstate,respfn,screening,sigma,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgr_init(Pawfgr,Dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
&                      gsqcutc_eff,gsqcutf_eff,gmet,k0) ! optional

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: nfftf,mgfftf
 real(dp),intent(out) :: ecut_eff,ecutdg_eff
 real(dp),intent(out),optional :: gsqcutf_eff,gsqcutc_eff
 type(dataset_type),intent(in) :: Dtset
 type(Pawfgr_type),intent(out) :: Pawfgr
!arrays
 real(dp),intent(in),optional :: gmet(3,3)
 integer,intent(out) :: ngfftc(18),ngfftf(18)
 real(dp),intent(in),optional :: k0(3)

!Local variables-------------------------------
 integer :: ii,nfftc_tot,nfftf_tot
 real(dp) :: boxcut,boxcutc
 character(len=500) :: msg

!************************************************************************

 DBG_ENTER("COLL")

 !@Pawfgr_type

 if ((present(gsqcutc_eff).or.present(gsqcutf_eff)).and.&
&    ((.not.present(gmet)).or.(.not.present(k0)))) then
   msg='To compute gsqcut[c,f]_eff, both k0 and gmet must be present as argument !'
   MSG_BUG(msg)
 end if

 ngfftc(:)=Dtset%ngfft(:)

 SELECT CASE (Dtset%usepaw)

 CASE (0)
  ! === Norm-conserving pseudopotentials ===
  nfftf=Dtset%nfft ; mgfftf=Dtset%mgfft ; ngfftf(:)=Dtset%ngfft(:)
  Pawfgr%usefinegrid=0
  ABI_ALLOCATE(Pawfgr%coatofin,(0))
  ABI_ALLOCATE(Pawfgr%fintocoa,(0))
  ecut_eff  =Dtset%ecut*Dtset%dilatmx**2
  ecutdg_eff=ecut_eff

 CASE (1)
  ! == PAW calculation ===
  if (any(Dtset%ngfftdg(1:3)/=Dtset%ngfft(1:3)) .and. Dtset%usewvl==0) then
! if (Dtset%pawecutdg>=1.0000001_dp*Dtset%ecut .and. Dtset%usewvl==0) then
   ! * Use fine FFT grid generated according to pawecutdg.
   nfftf=Dtset%nfftdg ; mgfftf=Dtset%mgfftdg ; ngfftf(:)=Dtset%ngfftdg(:)
   nfftc_tot =ngfftc(1)*ngfftc(2)*ngfftc(3)
   nfftf_tot =ngfftf(1)*ngfftf(2)*ngfftf(3)
   Pawfgr%usefinegrid=1
   ABI_ALLOCATE(Pawfgr%coatofin,(nfftc_tot))
   ABI_ALLOCATE(Pawfgr%fintocoa,(nfftf_tot))
   call indgrid(Pawfgr%coatofin,Pawfgr%fintocoa,nfftc_tot,nfftf_tot,ngfftc,ngfftf)
  else
   ! * Do not use fine FFT mesh. Simple transfer that can be done in parallel with only local info.
   nfftf=Dtset%nfft ; mgfftf=Dtset%mgfft ; ngfftf(:)=Dtset%ngfft(:)
   Pawfgr%usefinegrid=0
   ABI_ALLOCATE(Pawfgr%coatofin,(Dtset%nfft))
   ABI_ALLOCATE(Pawfgr%fintocoa,(Dtset%nfft))
   do ii=1,Dtset%nfft
    Pawfgr%coatofin(ii)=ii ; Pawfgr%fintocoa(ii)=ii
   end do
  end if
  ecutdg_eff=Dtset%pawecutdg*Dtset%dilatmx**2
  ecut_eff  =Dtset%ecut*Dtset%dilatmx**2

 CASE DEFAULT
  write(msg,'(a,i4)')' Wrong value of usepaw: ',Dtset%usepaw
  MSG_BUG(msg)
 END SELECT

! == Store useful dimensions in Pawfgr ===
 Pawfgr%nfftc=Dtset%nfft ; Pawfgr%mgfftc=Dtset%mgfft ; Pawfgr%ngfftc(:)=Dtset%ngfft(:)
 Pawfgr%nfft=nfftf       ; Pawfgr%mgfft=mgfftf       ; Pawfgr%ngfft (:)=ngfftf(:)

 !
 ! === Get boxcut for given gmet, ngfft, and ecut (center at k0) ===
 !     boxcut=ratio of basis sphere diameter to fft box side
 boxcut=-one
 if (Dtset%usepaw==1) then
   if (present(gsqcutc_eff)) then
     write(msg,'(2a)')ch10,' Coarse grid specifications (used for wave-functions):'
     call wrtout(std_out,msg,'COLL')
     call getcut(boxcutc,ecut_eff,gmet,gsqcutc_eff,Dtset%iboxcut,std_out,k0,ngfftc)
   end if
   if (present(gsqcutf_eff)) then
     write(msg,'(2a)')ch10,' Fine grid specifications (used for densities):'
     call wrtout(std_out,msg,'COLL')
     call getcut(boxcut,ecutdg_eff,gmet,gsqcutf_eff,Dtset%iboxcut,std_out,k0,ngfftf)
   end if
 else if (present(gsqcutc_eff)) then
   call getcut(boxcut,ecut_eff,gmet,gsqcutc_eff,Dtset%iboxcut,std_out,k0,ngfftc)
   gsqcutf_eff=gsqcutc_eff
 end if
 !
 ! === Check that boxcut>=2 if intxc=1; otherwise intxc must be set=0 ===
 if (boxcut>=zero .and. boxcut<two .and. Dtset%intxc==1) then
   write(msg,'(a,es12.4,5a)')&
&   ' boxcut=',boxcut,' is < 2.0  => intxc must be 0;',ch10,&
&   ' Need larger ngfft to use intxc=1.',ch10,&
&   ' Action : you could increase ngfft, or decrease ecut, or put intxc=0.'
   MSG_ERROR(msg)
 end if

 DBG_EXIT("COLL")

end subroutine pawfgr_init
!!***

!----------------------------------------------------------------------

!!****f* m_pawfgr/pawfgr_destroy
!! NAME
!!  pawfgr_destroy
!!
!! FUNCTION
!!  Deallocate pointers and nullify flags in a pawfgr structure
!!
!! SIDE EFFECTS
!!  pawfgr<type(pawfgr_type)>= Fine GRid parameters and related data
!!
!! PARENTS
!!      bethe_salpeter,eph,fourier_interpol,gstate,m_rec,respfn,screening,sigma
!!      wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgr_destroy(Pawfgr)

 implicit none

!Arguments ------------------------------------
!arrays
 type(Pawfgr_type),intent(inout) :: Pawfgr

!Local variables-------------------------------

! *************************************************************************

 DBG_ENTER("COLL")

!@Pawfgr_type

 if (associated(Pawfgr%coatofin))  then
   ABI_DEALLOCATE(Pawfgr%coatofin)
 end if
 if (associated(Pawfgr%fintocoa))  then
   ABI_DEALLOCATE(Pawfgr%fintocoa)
 end if

 Pawfgr%usefinegrid=0

 DBG_EXIT("COLL")

end subroutine pawfgr_destroy
!!***

!----------------------------------------------------------------------

!!****f* m_paw_ij/pawfgr_nullify
!! NAME
!!  pawfgr_nullify
!!
!! FUNCTION
!!  Nullify pointers and flags in a pawfgr structure
!!
!! SIDE EFFECTS
!!  Pawfgr<type(pawfgr_type)>=Fine GRid parameters and related data. Nullified in output
!!
!! PARENTS
!!      m_rec
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawfgr_nullify(Pawfgr)

 implicit none

!Arguments ------------------------------------
!arrays
 type(Pawfgr_type),intent(inout) :: Pawfgr

!Local variables-------------------------------

! *************************************************************************

 !@Pawfgr_type

  nullify(Pawfgr%coatofin)
  nullify(Pawfgr%fintocoa)

  Pawfgr%usefinegrid=0

end subroutine pawfgr_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_pawfgr/indgrid
!!
!! NAME
!! indgrid
!!
!! FUNCTION
!! Calculate the correspondance between the coarse grid and the fine grid
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! nfftc=total number of FFt grid=n1*n2*n3 for the coarse grid
!! nfftf=total number of FFt grid=n1*n2*n3 for the fine grid
!! ngfftc(18)=contain all needed information about 3D FFT, for the coarse grid,
!!        see ~abinit/doc/variables/vargs.htm#ngfft
!! ngfftf(18)=contain all needed information about 3D FFT, for the fine grid,
!!        see ~abinit/doc/variables/vargs.htm#ngfft
!!
!! OUTPUT
!! coatofin(nfftc)= index of the points of the coarse grid on the fine grid
!! fintocoa(nfftf)=index of the points of the fine grid on the
!!   coarse grid (=0 if the point of the fine grid does not belong to
!!   the coarse grid).
!!
!! PARENTS
!!      fourier_interpol,m_pawfgr,m_rec
!!
!! CHILDREN
!!
!! SOURCE

subroutine indgrid(coatofin,fintocoa,nfftc,nfftf,ngfftc,ngfftf)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftc,nfftf
!arrays
 integer,intent(in) :: ngfftc(18),ngfftf(18)
 integer,intent(out) :: coatofin(nfftc),fintocoa(nfftf)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,if1,if2,if3,ii,ing,n1c,n1f,n2c,n2f,n3c,n3f,narg1,narg2
!arrays
 integer :: id(3)
 integer,allocatable :: gc(:,:),gf(:,:)
 character(len=500) :: msg

! *************************************************************************
!

 DBG_ENTER("COLL")

 n1c=ngfftc(1);n2c=ngfftc(2);n3c=ngfftc(3)
 n1f=ngfftf(1);n2f=ngfftf(2);n3f=ngfftf(3)

 ABI_ALLOCATE(gc,(3,max(n1c,n2c,n3c)))
 do ii=1,3
   id(ii)=ngfftc(ii)/2+2
   do ing=1,ngfftc(ii)
     gc(ii,ing)=ing-(ing/id(ii))*ngfftc(ii)-1
   end do
 end do

 ABI_ALLOCATE(gf,(3,max(n1f,n2f,n3f)))
 do ii=1,3
   id(ii)=ngfftf(ii)/2+2
   do ing=1,ngfftf(ii)
     gf(ii,ing)=ing-(ing/id(ii))*ngfftf(ii)-1
   end do
 end do

 coatofin=0;fintocoa=0
 do i1=1,n1c
   do if1=1,n1f
     if(gc(1,i1)==gf(1,if1)) then
       do i2=1,n2c
         do if2=1,n2f
           if(gc(2,i2)==gf(2,if2)) then
             do i3=1,n3c
               do if3=1,n3f
                 if(gc(3,i3)==gf(3,if3)) then
                   narg1=i1+n1c*(i2-1+n2c*(i3-1))
                   narg2=if1+n1f*(if2-1+n2f*(if3-1))
                   coatofin(narg1)=narg2
                   fintocoa(narg2)=narg1
                   exit
                 end if
               end do
             end do
             exit ! To avoid N_fine * N_coarse scaling
           end if
         end do
       end do
       exit ! To avoid N_fine * N_coarse scaling
     end if
   end do
 end do

!Check coatofin to make sure there are no zeros!
 do ii=1,ubound(coatofin,1)
   if (coatofin(ii)==0) then
     msg = 'A zero was found in coatofin. Check that the fine FFT mesh is finer in each dimension than the coarse FFT mesh.'
     MSG_ERROR(msg)
   end if
 end do

 ABI_DEALLOCATE(gf)
 ABI_DEALLOCATE(gc)

 DBG_EXIT("COLL")

end subroutine indgrid
!!***

!----------------------------------------------------------------------

END MODULE m_pawfgr
!!***
