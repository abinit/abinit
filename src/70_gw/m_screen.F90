!!****m* ABINIT/m_screen
!! NAME
!!  m_screen
!!
!! FUNCTION
!!  Screening object used in the BSE code.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2025 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_screen

 use defs_basis
 use m_xmpi
 use m_hide_blas
 use m_errors
 use m_splines
 use m_abicore
 use m_kxc
 use m_screening
 use m_nctk
 use m_sort
 use m_yaml

 use m_gwdefs,         only : GW_TOLQ0, czero_gw
 use m_fstrings,       only : firstchar, endswith, strcat, itoa, sjoin
 use m_numeric_tools,  only : print_arr
 use m_geometry,       only : normv
 use m_crystal,        only : crystal_t
 use m_bz_mesh,        only : kmesh_t
 use m_gsphere,        only : gsphere_t
 use m_vcoul,          only : vcoul_t
 use m_io_screening,   only : read_screening, hscr_t, ncname_from_id, em1_ncname
 use m_ppmodel,        only : ppmodel_t, PPM_NONE, PPM_NOTAB

 implicit none

 private

 ! Flags defining the content of the %mat buffer in the fgg_t type.
 integer,public,parameter :: MAT_NOTYPE         = 0
 integer,public,parameter :: MAT_CHI0           = 1
 integer,public,parameter :: MAT_CHI            = 2
 integer,public,parameter :: MAT_EPSILON        = 3
 integer,public,parameter :: MAT_INV_EPSILON    = 4
 integer,public,parameter :: MAT_INV_EPSILON_M1 = 5
 integer,public,parameter :: MAT_W              = 6
 integer,public,parameter :: MAT_W_M1           = 7
 !
 ! Family vertex.
 integer,public,parameter :: VTX_FAMILY_NONE  = 0   ! No vertex correction.
 integer,public,parameter :: VTX_FAMILY_TDDFT = 1   ! TDDFT-based vertex.
 integer,public,parameter :: VTX_FAMILY_ADA   = 2   ! ADA vertex.
 !
 ! Test charge or test particle.
 integer,public,parameter :: VTX_TEST_CHARGE   = 0
 integer,public,parameter :: VTX_TEST_PARTICLE = 1
 !
 ! Named constants for the frequency mesh.
 integer,private,parameter :: WMESH_LINEAR    = 1
 integer,private,parameter :: WMESH_GAUSS_LEG = 2
 integer,private,parameter :: WMESH_TAN_GRID  = 3
 !
 ! Method used for the frequency integration.
 integer,public,parameter ::  WINT_NONE    = 0
 integer,public,parameter ::  WINT_PPMODEL = 1
 integer,public,parameter ::  WINT_CONTOUR = 2
 integer,public,parameter ::  WINT_AC      = 3
 !
 ! Parameters used for the model dielectric function.
 integer,public,parameter :: MDL_NONE      = 0
 integer,public,parameter :: MDL_BECHSTEDT = 1
 !
 ! Flags giving the status of the local buffers defined in fgg_t.
 integer,private,parameter :: MAT_NODATA    = 0
 integer,private,parameter :: MAT_ALLOCATED = 1
 integer,private,parameter :: MAT_STORED    = 2
 !
 ! Flags giving the status of the local buffers defined in fgg_t.
 integer,private,parameter :: FGG_QBZ_ISPOINTER  =1 ! Fgg_qbz is used to store the address in memory.
 integer,private,parameter :: FGG_QBZ_ISALLOCATED=2 ! Fgg_qbz is used as an allocable array.
!!***

!----------------------------------------------------------------------

!!****t* m_screen/screen_info_t
!! NAME
!! screen_info_t
!!
!! FUNCTION
!!  Container storing the parameters used to initialize a screen_t datatype or to
!!  calculate a new SCR file from the SUSC file containing the independent-particle
!!  polarizability.
!!
!! NOTES
!!  The list of parameters passed to screening_init is copied in W%Info.
!!  At present there is no need to provide a copy method since the structure does
!!  not contain pointers but such a method must be defined and used if
!!  dynamic entities are added to the datatype.
!!
!! SOURCE

type,public :: screen_info_t

  integer :: mat_type = MAT_NOTYPE
  ! Matrix identifier. See MAT_* flags.

  integer :: vtx_family = VTX_FAMILY_NONE
  ! Vertex correction family.

  integer :: invalid_freq = 0
  ! Sets the procedure to follow when a ppm frequency is invalid (negative or imaginary),
  ! see input variable gw_invalid_freq

  integer :: ixc = 0
  ! XC functional used for the TDDFT-based vertex.

  integer :: use_ada = 0
  ! >0 if ADA vertex is used.

  integer :: use_mdf = MDL_NONE
  ! >0 if model dielectric function is used.

  integer :: use_ppm = PPM_NONE
  ! >0 if ppmodel is used.

  integer :: vtx_test = VTX_TEST_CHARGE
  ! test charge or test particle.

  integer :: wint_method = WINT_NONE
  ! Defines the frequency integration technique. See WIN_ flags.
  ! NOTE that this flag can be changed at run time. For example
  ! one can switch from the CD to the PPm if the ppmodel parameters are in memory

  real(dp) :: ada_kappa = 2.1_dp
  ! Inverse smearing length used for ADA.

  real(dp) :: eps_inf = 12.0_dp
  ! Dielectric constant used for the model dielectric function.

  real(dp) :: drude_plsmf = zero
  ! Drude plasma frequency used for PPmodel 1.

  contains

    procedure :: print => screen_info_print

end type screen_info_t
!!***

!----------------------------------------------------------------------

!!****t* m_screen/fgg_t
!! NAME
!! fgg_t
!!
!! FUNCTION
!!  Object used to store F(G,G')(w) for a given q-point.
!!
!! SOURCE

 type,public :: fgg_t

  integer :: nomega
  ! Number of frequencies.

  integer :: npw
  ! Number of G vectors.

  integer :: nqlwl
  ! Number of points for the treatment of the long wave-length limit.

  integer :: has_mat = MAT_NODATA
  ! Flag giving the status of mat.

  complex(gwpc),allocatable :: mat(:,:,:)
  ! (npw, npw, nomega)
  ! The component of the two-point function $F_{G,G',w}$ for a given q.

  !complex(dpc),allocatable :: head(:,:,:)
  ! head(3,3,nomega)

  !complex(dpc),allocatable :: lwing(:,:,:)
  ! lwing(3,npwe,nomega)
  ! Lower wings

  !complex(dpc),allocatable :: uwing(:,:,:)
  ! uwing(3,npwe,nomega)
  ! Upper wings.

 contains

   procedure :: init => fgg_init   ! Creation method.
   !procedure :: free => fgg_free
 end type fgg_t

 public :: fgg_free   ! Free memory.
!!***

 interface fgg_free
   module procedure fgg_free_0D
   module procedure fgg_free_1D
 end interface fgg_free

!----------------------------------------------------------------------

!!****t* m_screen/screen_t
!! NAME
!! screen_t
!!
!! FUNCTION
!!  Object used to store the screening matrix in reciprocal space.
!!
!! SOURCE

 type,public :: screen_t

  ! scalars
  integer :: iomode               ! Flag defining the IO mode.
  integer :: debug_level=0        ! Internal Flag defining the debug level.
  integer :: mqmem                ! =0 for out-of-core solution, =nqibz if entire matrix is stored in memory.
  integer :: nI,nJ                ! Number of components (rows,columns) in chi|eps^-1. (1,1) if collinear.
  integer :: nqibz                ! Number of q-points in the IBZ used.
  integer :: nqlwl                ! Number of points used for the treatment of the long wave-length limit.
  integer :: nomega               ! Total Number of frequencies used.
  integer :: nomega_i             ! Number of purely imaginary frequencies used.
  integer :: nomega_r             ! Number of real frequencies used.
  integer :: npw                  ! Number of G vectors.
  integer :: prtvol               ! Verbosity level.
  integer :: has_ppmodel          ! 1 if PPmodel tables are stored.
  integer :: has_fgg              ! 1 if Fgg tables are stored.
  integer :: nfftf_tot
  integer :: nspden

  ! arrays
  integer :: ngfftf(18)          ! Info on the FFT mesh used for ae_rhor (used for the model dielectric function)

  real(dp),allocatable :: ae_rhor(:,:)
  ! ae_rhor(nfft,nspden)
  ! Density in real space used to construct the TDDFT kernel or the model dielectric function.
  ! NOTE that ae_rhor is given on the dense mesh as it contains the PAW onsite contribution.

  character(len=fnlen) :: fname = ABI_NOFILE  ! Name of the file used for the out-of-core solution.

  real(dp),allocatable :: qibz(:,:)
  ! (3,nqibz)
  ! q-points in reduced coordinates

  !type(kmesh_t) :: Qmesh
  ! Info the q-point sampling.

  real(dp),allocatable :: qlwl(:,:)
  ! (3,nqlwl)
  ! q-points used for the long wave-length limit treatment.

  complex(dpc),allocatable :: omega(:)
  ! (nomega)
  ! List of frequencies. Real frequencies are packed first.

  integer,allocatable :: gvec(:,:)
  ! (3,npw)
  ! G-vectors used to describe the two-point function (r.l.u.).

  !type(gsphere_t) :: Gsphere
  ! Info on the G-sphere. Note that the basis set does not depend on q.
  ! The G-vectors are ordered in shells to facilitate the application of symmetry operations.
  ! See m_gsphere.F90.

  logical,allocatable :: keep_qibz(:)
   ! (nqibz)
   ! Storage strategy: keep or not keep Em1(q) in memory.

  type(fgg_t),pointer :: Fgg(:)
  ! (nqibz)
  ! F_{G,G'}(q,w) for q in the IBZ.

  integer :: fgg_qbz_stat = FGG_QBZ_ISPOINTER
  ! Status of Fgg_qbz

  integer :: fgg_qbz_idx = 0
  ! The index of the q-point in BZ pointed by Fgg_qbz. Used for debugging purpose.

  type(fgg_t),pointer :: Fgg_qbz  => null()
  ! Buffer used for storing F_GG' at the point q_bz in the BZ
  ! If q_bz is in the IBZ, Fgg_qbz *points* to Fgg(iq_ibz)
  ! If q_bz is not in the IBZ, Fgg_qbz is *allocated* and used to store the symmetrized matrix.

  type(ppmodel_t) :: ppm
  ! Structure storing the plasmon-pole parameters.

  type(screen_info_t) :: Info
  ! Parameters used to construct the screening.

contains

  procedure :: nullify => screen_nullify
    ! Nullify all pointers before use.

  procedure :: init => screen_init
    ! Creation method.

  procedure :: print => screen_print
    ! Print info on object

  procedure :: free => screen_free
   ! Free dynamic memory

  procedure :: rotate_iqbz => screen_rotate_iqbz
    ! Prepare the object for applying W_qbz.

  procedure :: w0gemv => screen_w0gemv
    ! Matrix vector multiplication \sum_{G'} F_{G,G') |u(G')>.

  procedure :: calc_sigc => screen_calc_sigc
    ! Compute the frequency convolution.

  procedure :: ihave_fgg => screen_ihave_fgg
    ! Inquire the processor whether it has a particular F_{GG')(q) and with which status.

end type screen_t
!!***

contains
!----------------------------------------------------------------------

!!****f* m_screen/screen_info_print
!! NAME
!!  screen_info_print
!!
!! FUNCTION
!!  Printout object
!!
!! INPUTS
!!  [header]=String to be printed as header for additional info.
!!  [unit]=Unit number for output
!!  [prtvol]=Verbosity level
!!  [mode_paral]=Either "COLL" or "PERS"
!!
!! SOURCE

subroutine screen_info_print(W_info, header, unit, mode_paral, prtvol)

!Arguments ------------------------------------
!scalars
 class(screen_info_t),intent(in) :: W_info
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header

!Local variables-------------------------------
 integer :: my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg
! *********************************************************************

 !@screen_info_t
 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Info on the screen_info_t% object ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

!integer
 write(msg,'(a,i3)')" mat_type    ",W_info%mat_type
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,i3)')" vtx_family  ",W_info%vtx_family
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,i3)')" invalid_freq",W_info%invalid_freq
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,i3)')" ixc         ",W_info%ixc
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,i3)')" use_ada     ",W_info%use_ada
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,i3)')" use_mdf     ",W_info%use_mdf
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,i3)')" use_ppm     ",W_info%use_ppm
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,i3)')" vtx_test    ",W_info%vtx_test
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,i3)')" wint_method ",W_info%wint_method
 call wrtout(my_unt,msg,my_mode)

!real
 write(msg,'(a,f8.3)')" ada_kappa   ",W_info%ada_kappa
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,f8.3)')" eps_inf     ",W_info%eps_inf
 call wrtout(my_unt,msg,my_mode)
 write(msg,'(a,f8.3)')" drude_plsmf ",W_info%drude_plsmf
 call wrtout(my_unt,msg,my_mode)

end subroutine screen_info_print
!!***

!----------------------------------------------------------------------

!!****f* m_screen/fgg_free_0D
!! NAME
!! fgg_free_0D
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! SOURCE

subroutine fgg_free_0D(Fgg)

!Arguments ------------------------------------
 type(fgg_t),intent(inout) :: Fgg
! *************************************************************************

 ABI_SFREE(Fgg%mat)
 Fgg%has_mat = MAT_NODATA

end subroutine fgg_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_screen/fgg_free_1D
!! NAME
!! fgg_free_1D
!!
!! FUNCTION
!! Deallocate all the pointers in the structure that result to be associated.
!!
!! INPUT
!!  [keep_qibz(:)]=Optional logical mask used to select the q-points that are deallocated.
!!
!! SOURCE

subroutine fgg_free_1D(Fgg, keep_qibz)

!Arguments ------------------------------------
!scalars
 type(fgg_t),intent(inout) :: Fgg(:)
 logical,optional,intent(in) :: keep_qibz(:)

!Local variables ------------------------------
!scalars
 integer :: iq_ibz
 logical :: keep_it
! *************************************************************************

 do iq_ibz=LBOUND(Fgg,DIM=1),UBOUND(Fgg,DIM=1)
   keep_it = .FALSE.; if (PRESENT(keep_qibz)) keep_it = keep_qibz(iq_ibz)
   if (.not. keep_it) call fgg_free_0D(Fgg(iq_ibz))
 end do

end subroutine fgg_free_1D
!!***

!----------------------------------------------------------------------

!!****f* m_screen/fgg_init
!! NAME
!! fgg_init
!!
!! FUNCTION
!! Initialize the structure allocating the memory and initializing the internal variables.
!!
!! INPUT
!!  npw
!!  nqlwl
!!
!! SOURCE

subroutine fgg_init(Fgg, npw, nomega, nqlwl)

!Arguments ------------------------------------
!scalars
 class(fgg_t),intent(inout) :: Fgg
 integer,intent(in) :: npw, nqlwl, nomega

!Local variables ------------------------------
 integer :: ierr

! *************************************************************************

 Fgg%nomega = nomega; Fgg%npw = npw; Fgg%nqlwl = nqlwl

 if (npw > 0 .and. nomega > 0) then
   ABI_MALLOC_OR_DIE(Fgg%mat, (npw, npw, nomega), ierr)
   Fgg%has_mat = MAT_ALLOCATED
 end if

end subroutine fgg_init
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_fgg_qbz_set
!! NAME
!!  screen_fgg_qbz_set
!!
!! FUNCTION
!!  Helper function used to perform the setup W%Fgg_qbz setting also the internal
!!  flag that defines its status.
!!
!! INPUTS
!!  iq_bz=Index of the q-point in the BZ.
!!  nqlwl=Number of wings wanted.
!!  how= "Pointer" is a true pointer is wanted.
!!       "Allocated" if memory has to be allocated.
!!
!! NOTES
!!  iq_bz and nqlwl are not used if how="Pointer".
!!
!! SOURCE

subroutine screen_fgg_qbz_set(screen, iq_bz, nqlwl, how)

!Arguments ------------------------------------
!scalars
 class(screen_t),intent(inout) :: screen
 integer,intent(in) :: iq_bz,nqlwl
 character(len=*),intent(in) :: how

!Local variables ------------------------------
!scalars
 !character(len=500) :: msg

!************************************************************************

 screen%fgg_qbz_idx = iq_bz ! Save the index of the q-point in the BZ.

 if (firstchar(how, (/"P"/)) ) then
   ! We want a pointer.
   select case (screen%fgg_qbz_stat)
   case (FGG_QBZ_ISALLOCATED)
     call fgg_free_0D(screen%Fgg_qbz)
     ABI_FREE(screen%Fgg_qbz)
     nullify(screen%Fgg_qbz)
     screen%fgg_qbz_stat = FGG_QBZ_ISPOINTER

   case (FGG_QBZ_ISPOINTER)
     ! Set it to null().
     nullify(screen%Fgg_qbz)

   case default
     ABI_ERROR(sjoin("Wrong status:", itoa(screen%fgg_qbz_stat)))
   end select

 else if (firstchar(how, (/"A"/)) ) then
   ! We want an allocatable array.

   select case (screen%fgg_qbz_stat)
   case (FGG_QBZ_ISPOINTER)
     ! Allocate memory
     nullify(screen%Fgg_qbz)
     ABI_MALLOC(screen%Fgg_qbz,)

     call screen%Fgg_qbz%init(screen%npw, screen%nomega, nqlwl)
     screen%fgg_qbz_stat = FGG_QBZ_ISALLOCATED

   case (FGG_QBZ_ISALLOCATED)
     screen%Fgg_qbz%has_mat = MAT_ALLOCATED  ! STORED --> ALLOCATED

   case default
     ABI_ERROR(sjoin("Wrong status:", itoa(screen%fgg_qbz_stat)))
   end select

 else
   ABI_BUG(sjoin("Wrong how:", how))
 end if

end subroutine screen_fgg_qbz_set
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_ihave_fgg
!! NAME
!!  screen_ihave_fgg
!!
!! FUNCTION
!!  Inquire the processor whether it has a particular F_{GG')(q) and with which status.
!!
!! INPUTS
!!  iq_ibz=k-point index
!!  [how]=string defining which status is checked. By default the function returns
!!     .TRUE. if the wave is either MAT_ALLOCATED or MAT_STORED.
!!     Possible mutually exclusive values: "Allocated", "Stored".
!!     Only the first character is checked (no case-sensitive)
!!
!! NOTES
!!   A zero index can be used to inquire the status of the full set of q-points.
!!
!! SOURCE

logical pure function screen_ihave_fgg(screen, iq_ibz, how)

!Arguments ------------------------------------
!scalars
 class(screen_t),intent(in) :: screen
 integer,intent(in) :: iq_ibz
 character(len=*),optional,intent(in) :: how

!Local variables ------------------------------
 integer :: ii, check(2)

!************************************************************************

 check = [MAT_ALLOCATED, MAT_STORED]
 if (PRESENT(how)) then
   if (firstchar(how, (/"A","a"/))) check = [MAT_ALLOCATED, MAT_ALLOCATED]
   if (firstchar(how, (/"S","s"/))) check = [MAT_STORED, MAT_STORED]
 end if

 if (iq_ibz > 0) then
   screen_ihave_fgg = (screen%Fgg(iq_ibz)%has_mat == check(1) .or.&
                       screen%Fgg(iq_ibz)%has_mat == check(2) )
 else
   ! check the status of the full set of q-tables.
   screen_ihave_fgg=.TRUE.
   do ii=1,screen%nqibz
     screen_ihave_fgg = screen_ihave_fgg .and. &
                     (screen%Fgg(ii)%has_mat == check(1) .or.&
                      screen%Fgg(ii)%has_mat == check(2) )
   end do
 end if

end function screen_ihave_fgg
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_nullify
!! NAME
!! screen_nullify
!!
!! FUNCTION
!! Initialize the pointers to null()
!!
!! SOURCE

subroutine screen_nullify(screen)

!Arguments ------------------------------------
 class(screen_t),intent(inout) :: screen

! *************************************************************************

 nullify(screen%Fgg_qbz); screen%fgg_qbz_stat=FGG_QBZ_ISPOINTER ! Needed since the initial status is undefined.
 nullify(screen%Fgg)

end subroutine screen_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_free
!! NAME
!! screen_free
!!
!! FUNCTION
!! Free the memory allocate in the datatype.
!!
!! SOURCE

subroutine screen_free(screen)

!Arguments ------------------------------------
 class(screen_t),intent(inout) :: screen

! *************************************************************************

 ! integer
 ABI_SFREE(screen%gvec)

 !real
 ABI_SFREE(screen%ae_rhor)
 ABI_SFREE(screen%qibz)
 ABI_SFREE(screen%qlwl)

 !complex
 ABI_SFREE(screen%omega)

 ! logical
 ABI_SFREE(screen%keep_qibz)

 ! types
 ! Here be careful with dangling pointers.
 ! First Fgg_qbz that might point to one of the %Fgg then %Fgg.
 select case (screen%fgg_qbz_stat)

 case (FGG_QBZ_ISALLOCATED)
   call fgg_free_0D(screen%Fgg_qbz)
   ABI_FREE(screen%Fgg_qbz)
   nullify(screen%Fgg_qbz)
   screen%fgg_qbz_stat=FGG_QBZ_ISPOINTER

 case (FGG_QBZ_ISPOINTER)
   nullify(screen%Fgg_qbz)
   screen%fgg_qbz_stat=FGG_QBZ_ISPOINTER

 case default
   continue
 end select

 ! Free the Fgg matrices.
 if (associated(screen%Fgg)) then
   call fgg_free(screen%Fgg)
   ABI_FREE(screen%Fgg)
 end if

 ! Free the plasmon pole tables.
 call screen%ppm%free()

end subroutine screen_free
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_print
!! NAME
!! screen_print
!!
!! FUNCTION
!! Print info on the object.
!!
!! SOURCE

subroutine screen_print(screen, units, header)

!Arguments ------------------------------------
 class(screen_t),intent(in) :: screen
 integer,intent(in) :: units(:)
 character(len=*),optional,intent(in) :: header

!Local variables-------------------------------
 character(len=500) :: msg
 type(yamldoc_t) :: ydoc

! *************************************************************************

 msg = ' ==== Info on the screen_t object ==== '; if (present(header)) msg=' ==== '//trim(adjustl(header))//' ==== '
 call wrtout(units, msg)

 ydoc = yamldoc_open('screen_params') !, width=11, real_fmt='(3f8.3)')
 call ydoc%add_int("nomega_r", screen%nomega_r)
 call ydoc%add_int("nomega_i", screen%nomega_i)
 !call ydoc%add_real("drude_plsmf", ppm%drude_plsmf)
 !call ydoc%add_int1d("has_qibz", ppm%has_qibz)

 call ydoc%write_units_and_free(units)

end subroutine screen_print
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_init
!! NAME
!!  screen_init
!!
!! FUNCTION
!!  Initialize basic dimensions and other important arrays
!!  starting from a file containing either epsilon^{-1} (_SCR) or chi0 (_SUSC).
!!
!! INPUTS
!!  W_Info<screen_info_t>=The list of parameters used to construct the screen function.
!!  Cryst<crystal_t>=Info on the unit cell.
!!  Qmesh<kmesh_t>=Info on the Q-mesh.
!!  Gsph<gsphere_t>=Info on the plane-wave basis set used for the two-point function.
!!  Vcp<vcoul_t>=Structure gathering data on the Coulombian interaction.
!!  ifname=The name of the external file used to read the matrix.
!!  id_required=Identifier used to specify the type of two-point function that is wanted.
!!  iomode=Option defining the file format of the external file.
!!  mqmem=0 for out-of-core solution, /=0 if entire matrix has to be stored in memory.
!!  npw_asked=Number of G-vector to be used in the calculation, if <=0 use Max allowed number.
!!  ngfftf(18)=Info on the (fine) mesh used for the density.
!!  nfftf_tot=Total number of point in the FFT mesh for ae_rhor
!!  nsppol=Numer of independent spin polarizations.
!!  nspden=Number of spin density components in ae_rhor
!!  ae_rhor(nfftf_tot,nspden)
!!  prtvol=Verbosity level.
!!  comm=MPI communicator.
!!
!! SOURCE

subroutine screen_init(screen, W_Info, Cryst, Qmesh, Gsph, Vcp, ifname, mqmem, npw_asked, &
                       iomode, ngfftf, nfftf_tot, nsppol, nspden, ae_rhor, prtvol, comm)

!Arguments ------------------------------------
!scalars
 class(screen_t),intent(out) :: screen
 integer,intent(in) :: mqmem,iomode,npw_asked,comm,prtvol,nsppol, nfftf_tot,nspden
 character(len=fnlen),intent(in) :: ifname
 type(crystal_t),intent(in) :: Cryst
 type(gsphere_t),intent(in) :: Gsph
 type(vcoul_t),intent(in) :: Vcp
 type(kmesh_t),intent(in) :: Qmesh
 type(screen_info_t),intent(in) :: W_Info
!arrays
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(in) :: ae_rhor(nfftf_tot,nspden)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: option_test,approx_type,ixc_required,id_required !,nkxc
 integer :: fform,my_rank,mat_type_read, nqibz,nomega,iq_ibz,npw,nqlwl
 integer :: nI,nJ,iq_bz,mdf_type,ppmodel,ierr, iw,qsort,ii
 real(dp) :: eps_inf,drude_plsmf
 logical :: free_Fgg,found,from_file,is_qeq0,remove_dgg !only_one_kpt,
 character(len=500) :: msg
 character(len=fnlen) :: sus_fname,scr_fname
 character(len=nctk_slen) :: varname
 type(hscr_t) :: Hscr
!arrays
 integer :: units(2), g0(3), iperm(Qmesh%nibz)
 real(dp) :: wt_list(Qmesh%nibz)
 !complex(gwpc),ABI_CONTIGUOUS pointer :: em1_ggw(:,:,:)

! *********************************************************************

 DBG_ENTER("COLL")
 ABI_UNUSED(nsppol)

 my_rank = xmpi_comm_rank(comm); units = [std_out, ab_out]
 call screen%nullify()

 ! Initialize basic parameters
 screen%Info = W_info
 call screen%Info%print(header="W info", unit=std_out)

 id_required  = W_Info%mat_type
 approx_type  = W_Info%vtx_family
 option_test  = W_Info%vtx_test
 ixc_required = W_Info%ixc
 varname = ncname_from_id(id_required)

 if (all(id_required /= [MAT_INV_EPSILON])) then
   ABI_ERROR(sjoin("id_required:", itoa(id_required), " not available"))
 end if

 ! This part must be rationalized.
 remove_dgg = (id_required == MAT_W_M1)

 if (screen%Info%use_mdf == MDL_NONE) screen%fname = ifname
 screen%nI = 1; screen%nJ = 1

 ! The q-point sampling is initialized from qmesh.
 screen%nqibz = Qmesh%nibz
 ABI_MALLOC(screen%qibz, (3, screen%nqibz))
 screen%qibz= Qmesh%ibz

 screen%mqmem = mqmem; if (screen%mqmem /= 0) screen%mqmem = screen%nqibz !; screen%mqmem = 0

 ABI_MALLOC(screen%keep_qibz, (screen%nqibz))
 screen%keep_qibz = .TRUE.; if (screen%mqmem == 0) screen%keep_qibz = .False.

 if (screen%mqmem /= 0 .and. screen%mqmem < screen%nqibz) then
   ! Keep in memory the most representative q-points.
   screen%keep_qibz = .FALSE.
   wt_list = Qmesh%wt; iperm = (/(ii,ii=1,Qmesh%nibz)/)
   call sort_dp(Qmesh%nibz, wt_list, iperm, tol12)
   do qsort=Qmesh%nibz,Qmesh%nibz-mqmem+1,1
     iq_ibz = iperm(qsort)
     screen%keep_qibz(iq_ibz) = .TRUE.
   end do
 end if

 screen%fgg_qbz_idx = 0
 screen%iomode = iomode
 screen%prtvol = prtvol
 screen%has_ppmodel = 0; if (screen%Info%use_ppm /= PPM_NONE) screen%has_ppmodel = 1

 ! Copy the AE density for the model dielectric function or for the vertex corrections.
 screen%nspden     = nspden
 screen%ngfftf     = ngfftf
 screen%nfftf_tot  = nfftf_tot

 ABI_MALLOC(screen%ae_rhor, (nfftf_tot, nspden))
 screen%ae_rhor = ae_rhor

 free_Fgg = .FALSE.

 screen%has_fgg = 0; if (ANY(screen%Info%wint_method == [WINT_CONTOUR, WINT_AC])) screen%has_fgg = 1

 if (screen%has_fgg > 0 .and. screen%has_ppmodel > 0) then
   ABI_WARNING("Both PPmodel tables and F_(GG')(q,w) are stored in memory")
 end if

 ! Default values used if external file is not read.
 nqlwl = 0; nomega = 1
 screen%npw = npw_asked

 ! Model dielectric function does not require any external file.
 from_file = (screen%Info%use_mdf == MDL_NONE)

 if (from_file) then
   ! Open file and check its content.
   if (endswith(screen%fname, ".nc")) screen%iomode = IO_MODE_ETSF
   call hscr%from_file(screen%fname, fform, comm)
   ! Echo of the header
   if (my_rank == master .and. screen%prtvol > 0) call hscr%print()

   mat_type_read = Hscr%id
   nqlwl         = Hscr%nqlwl
   nomega        = Hscr%nomega
 end if

 screen%nqlwl  = nqlwl
 screen%nomega = nomega

 ABI_MALLOC(screen%qlwl, (3, screen%nqlwl))
 ABI_MALLOC(screen%omega, (screen%nomega))

 if (from_file) then
   screen%qlwl  = Hscr%qlwl
   screen%omega = Hscr%omega

   ! G-vectors.
   screen%npw = Hscr%npwe
   if (npw_asked > 0) then
     if (npw_asked > Hscr%npwe) then
       write(msg,'(a,i0,2a,i0)') &
        'The number of G-vectors saved on file is less than the value required: ',npw_asked,ch10,&
        'Calculation will proceed with the Max available npw: ',Hscr%npwe
       ABI_WARNING(msg)
     else
       screen%npw = npw_asked ! Redefine the no. of G"s for W.
       write(msg,'(a,i0,2a,i0)')&
        'The number of G-vectors saved on file is larger than the value required: ',npw_asked,ch10,&
        'Calculation will proceed with npw: ',screen%npw
       ABI_COMMENT(msg)
     end if
   end if

   ! Here consistency check on G-vectors and q-points.
   if (ANY(Hscr%gvec(:,1:screen%npw) /= Gsph%gvec(:,1:screen%npw))) then
     !write(std_out) W%gvec, Gsph%gvec
     ABI_ERROR("Hscr%gvec /= Gsph%gvec(1:W%npw)")
   end if
   ABI_CHECK(Hscr%nqibz == Qmesh%nibz, "Mismatch in the number of q-points in the IBZ")
   ierr = 0
   do iq_ibz=1,Hscr%nqibz
     if (ANY(ABS(Qmesh%ibz(:,iq_ibz) - Hscr%qibz(:,iq_ibz)) > tol6) ) then
       ierr = ierr + 1
       write(std_out,'(i0,2(3f7.3,1x))')iq_ibz, Qmesh%ibz(:,iq_ibz), Hscr%qibz(:,iq_ibz)
     end if
   end do
   ABI_CHECK(ierr == 0, "Wrong ordering in q-point list, Aborting now")
 end if

 ABI_MALLOC(screen%gvec, (3, screen%npw))
 screen%gvec = Gsph%gvec(:,1:screen%npw)

 ! Frequency mesh.
 screen%nomega_r = 1; screen%nomega_i = 0
 if (screen%nomega == 2 ) then
   screen%nomega_r = 1; screen%nomega_i = 1
 else
   ! Real frequencies are packed in the first locations.
   screen%nomega_r = 1
   do iw=1,screen%nomega
     if (DBLE(screen%omega(iw))>0.001*Ha_eV) screen%nomega_r=iw
   end do
   screen%nomega_i = screen%nomega - screen%nomega_r
 end if

 ! ------------------------------ Initialization completed --------------------------------
 !
 ! Just to keep the code below more readable.
 npw    = screen%npw
 nqibz  = screen%nqibz
 nomega = screen%nomega
 nI     = screen%ni
 nJ     = screen%nj
 ABI_MALLOC(screen%Fgg, (nqibz))

 if (from_file) then

   ! Read ab-initio em1 from file.
   select case (mat_type_read)
   case (MAT_INV_EPSILON)
     call wrtout(std_out, strcat("Em1 will be initialized from SCR file: ", screen%fname))

   case (MAT_CHI0)
     ! Should Write new SCR file.
     ABI_ERROR("Not coded yet")
     sus_fname = screen%fname; scr_fname="TESTING_SUS2SCR"

     screen%fname = scr_fname  ! Change the name of the file associated to W.

   case default
     write(msg,'(a,i0)')" Unsupported conversion from mat_type ",mat_type_read
     ABI_ERROR(msg)
   end select

   ! Begin reading.
   do iq_ibz=1,nqibz
     if (.not. screen%keep_qibz(iq_ibz)) then
       !call wrtout(std_out, strcat("Skipping iq_ibz: ",itoa(iq_ibz)))
       CYCLE
     end if

     nqlwl = 0; is_qeq0 = (normv(screen%qibz(:,iq_ibz),Cryst%gmet,'G') < GW_TOLQ0)
     if (is_qeq0) nqlwl=screen%nqlwl

     ! Allocate F_{GG'}(w)
     call screen%Fgg(iq_ibz)%init(npw, nomega, nqlwl)

     ! Read data from file (use MPI-IO if possible)
     if (screen%iomode /= IO_MODE_ETSF .and. xmpi_mpiio == 1) then
       !call wrtout(std_out, "read_screening with MPI_IO")
       call read_screening(varname, screen%fname, npw, 1, nomega, screen%Fgg(iq_ibz)%mat, IO_MODE_MPI, comm, iqiA=iq_ibz)
     else
       call read_screening(varname, screen%fname, npw, 1, nomega, screen%Fgg(iq_ibz)%mat, screen%iomode, comm, iqiA=iq_ibz)
     end if

     ! W contains Em1 and is ready to use.
     screen%Fgg(iq_ibz)%has_mat = MAT_STORED
   end do

 else

   ! Model dielectric function. Only epsm-1 is supported here.
   call wrtout(std_out," Calculating model dielectric function... ")
   ABI_CHECK(screen%nomega == 1, "Cannot use nomega > 1 in model dielectric function")

   do iq_ibz=1,nqibz
     if (.not.screen%keep_qibz(iq_ibz)) CYCLE

     ! The wings are not used here.
     nqlwl=0; is_qeq0= (normv(screen%qibz(:,iq_ibz),Cryst%gmet,'G')<GW_TOLQ0)

     ! Calculate the model. Note that mdielf awaits an index in the BZ.
     found = qmesh%has_bz_item(Qmesh%ibz(:,iq_ibz),iq_bz,g0)
     if (.not.found .or. any(g0 /= 0)) then
       ABI_ERROR("Problem in retrieving ibz point")
     end if

     ! Allocate F_{GG'}(w).
     call screen%Fgg(iq_ibz)%init(npw, nomega, nqlwl)

     eps_inf  =  screen%Info%eps_inf
     mdf_type =  screen%Info%use_mdf
     !em1_ggw  => screen%Fgg(iq_ibz)%mat

     ! Construct W TODO check the new implementation.
     call screen_mdielf(iq_bz,npw,nomega,mdf_type,eps_inf,Cryst,Qmesh,Vcp,Gsph,&
                        nspden,nfftf_tot,ngfftf,ae_rhor,"EM1",screen%Fgg(iq_ibz)%mat,comm)

     screen%Fgg(iq_ibz)%has_mat = MAT_STORED

     if (screen%prtvol > 0) then
       do iw=1,nomega
         write(msg,'(a,i3,a,i4,a)')'  Model symmetrical e^{-1} (q=',iq_ibz,', omega=',iw,', G,G'')'
         call wrtout(std_out,msg)
         call print_arr([std_out], screen%Fgg(iq_ibz)%mat(:,:,iw))
       end do
     end if
   end do ! iq_ibz
 end if

 ! Init plasmon-pole parameters from em1.
 if (screen%has_ppmodel > 0) then
   call wrtout(std_out, " Calling ppm_init ...")
   ppmodel = screen%Info%use_ppm; drude_plsmf = screen%Info%drude_plsmf
   call screen%ppm%init(screen%mqmem, screen%nqibz, screen%npw, ppmodel, drude_plsmf, screen%Info%invalid_freq)
   !call screen%ppm%print(units)

   do iq_ibz=1,nqibz
     if (screen%ihave_fgg(iq_ibz, how="Stored")) then
       !call wrtout(std_out, sjoin(" Calling ppm%new_setup for iq_ibz:", itoa(iq_ibz)))
       call screen%ppm%new_setup(iq_ibz, Cryst, Qmesh, npw, nomega, screen%omega, &
                                 screen%Fgg(iq_ibz)%mat, nfftf_tot, Gsph%gvec, ngfftf, screen%ae_rhor(:,1))
     end if
   end do

 end if
 !stop

 ! Deallocate Fgg if the matrices are not needed anymore.
 if (free_Fgg) then
   call screen_fgg_qbz_set(screen, 0, 0, "Pointer") ! Avoid dangling pointer.
   call fgg_free(screen%Fgg, keep_qibz=screen%keep_qibz)
 end if

 if (from_file) call Hscr%free()

 DBG_EXIT("COLL")

end subroutine screen_init
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_rotate_iqbz
!! NAME
!!  screen_rotate_iqbz
!!
!! FUNCTION
!!  Modify the status of the object so that the symmetrized component F(q_bz)_GG' is calculated
!!  (if needed) and is made available in the internal buffer. This routine must be called before
!!  performing any operation that involves the symmetrized component of the two-point function.
!!
!! INPUTS
!!  iq_bz=Index of the q-point in the BZ where F(q_bz)_GG' is wanted.
!!  Cryst=Crystal structure.
!!  Gsph=The G-sphere
!!  Qmesh=Structure defining the q-mesh used for sample the BZ.
!!
!! SIDE EFFECTS
!!   screen%ppm
!!   screen%Fgg_qbz
!!
!! NOTES
!!  In the present implementation, we are not considering a possible umklapp vector G0 in the
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required
!!  to reconstruct the BZ.
!!
!! SOURCE

subroutine screen_rotate_iqbz(screen, iq_bz, Cryst, Gsph, Qmesh, Vcp)

!Arguments ------------------------------------
!scalars
 class(screen_t),intent(inout) :: screen
 integer,intent(in) :: iq_bz
 type(crystal_t),intent(in) :: Cryst
 type(gsphere_t),intent(in) :: Gsph
 type(kmesh_t),intent(in) :: Qmesh
 type(vcoul_t),intent(in) :: Vcp

!Local variables-------------------------------
!scalars
 integer,parameter :: nqlwl0=0
 integer :: iq_ibz,isym_q,itim_q,npw,nomega,nqibz,mdf_type
 real(dp) :: eps_inf
 logical :: q_isirred
 !character(len=500) :: msg
!arrays
 real(dp) :: qbz(3)

! *********************************************************************

 DBG_ENTER("COLL")

 npw = screen%npw; nqibz = screen%nqibz; nomega = screen%nomega
 call qmesh%get_bz_item(iq_bz, qbz, iq_ibz, isym_q, itim_q, isirred=q_isirred)

 ! ========================================================
 ! ==== Branching for in-core or out-of-core solutions ====
 ! ========================================================
 if (screen%ihave_fgg(iq_ibz, how="Stored")) then

   if (q_isirred) then
     ! Symmetrization is not needed. Target the data in memory.
     call screen_fgg_qbz_set(screen, iq_bz, nqlwl0, "Pointer")
     screen%Fgg_qbz => screen%Fgg(iq_ibz)
   else
     ! Allocate space. ! TODO Wings are not symmetrized but oh well
     call screen_fgg_qbz_set(screen, iq_bz, nqlwl0, "Allocate")  ! Dimensions should not be changed.

     ! Out-of-place symmetrization.
     !em1_qibz => screen%Fgg(iq_ibz)%mat; em1_qbz  => screen%Fgg_qbz%mat
     call em1_symmetrize_op(iq_bz, npw, nomega, Gsph, Qmesh, screen%Fgg(iq_ibz)%mat, screen%Fgg_qbz%mat)
   end if

   if (screen%has_ppmodel > 0) then
     ! Symmetrize the ppmodel tables: em1_qibz => W%Fgg(iq_ibz)%mat
     call screen%ppm%rotate_iqbz(iq_bz, Cryst, Qmesh, Gsph, npw, nomega, screen%omega, screen%Fgg(iq_ibz)%mat, &
                                 screen%nfftf_tot, screen%ngfftf, screen%ae_rhor(:,1))

     !call screen%ppm_get_qbz(Gsph, Qmesh, iq_bz, botsq_qbz, otq_qbz, eig_qbz)
   end if

 else if (screen%ihave_fgg(iq_ibz, how="Allocated")) then
   ABI_ERROR("Fgg_iqibz is allocated but not initialized!")

 else
   ! Out of core branch

   if (screen%fgg_qbz_idx /= iq_bz) then
     ! Must compute em1_qbz here. em1_qbz => W%Fgg_qbz%mat
     ! Allocate the BZ buffer.
     call screen_fgg_qbz_set(screen, iq_bz, nqlwl0, "Allocate")

     if (screen%Info%use_mdf /= MDL_NONE) then
       ! Compute the model-dielectric function at qbz on-the fly and in sequential
       !call wrtout(std_out,"Will compute MDF on the fly")
       call screen_mdielf(iq_bz,npw,nomega,screen%Info%use_mdf,screen%Info%eps_inf,Cryst,Qmesh,Vcp,Gsph,&
                          screen%nspden,screen%nfftf_tot,screen%ngfftf,screen%ae_rhor,"EM1", &
                          screen%Fgg_qbz%mat,xmpi_comm_self)

     else
       ! Read W(q_ibz) and symmetrize it (do this only if we don't have the correct q_bz in memory).
       call wrtout(std_out,sjoin("Out of core with file: ",screen%fname))
       call read_screening(em1_ncname, screen%fname, npw, 1, nomega, screen%Fgg_qbz%mat, &
                           screen%iomode, xmpi_comm_self, iqiA=iq_ibz)

       ! In-place symmetrization to get the q-point in the BZ.
       if (.not. q_isirred) then
         call em1_symmetrize_ip(iq_bz, npw, nomega, Gsph, Qmesh, screen%Fgg_qbz%mat)
       end if
     end if

     screen%Fgg_qbz%has_mat = MAT_STORED
   end if

   ABI_CHECK(screen%Fgg_qbz%has_mat == MAT_STORED, "Wrong has_mat")

   ! Ppmodel calculations with ppm tables in memory.
   ! TODO treat the case in which IBZ tables are stored in memory.
   if (screen%has_ppmodel > 0) then
     ABI_ERROR("Not implemented error")
     ! Symmetrize the ppmodel using em1_qibz.
     call screen%ppm%rotate_iqbz(iq_bz, Cryst, Qmesh, Gsph, npw, nomega, screen%omega, &
                                 screen%Fgg_qbz%mat, screen%nfftf_tot, screen%ngfftf, screen%ae_rhor(:,1))
   end if
 end if

 ! Calculate model dielectric function for this q-point in the BZ.
 eps_inf = screen%Info%eps_inf; mdf_type = screen%Info%use_mdf

 ! Model dielectric function. Only epsm-1 is supported here.
 !call wrtout(std_out," Calculating model dielectric function... ")
 !ABI_CHECK(W%nomega==1,"Cannot use nomega>1 in model dielectric function")

 ! screen%Fgg_qbz%mat

 !%  call screen_mdielf(iq_bz,npw,nomega,mdf_type,eps_inf,Cryst,Qmesh,Vcp,Gsph,&
 !% &  screen%nspden,screen%nfftf_tot,screen%ngfftf,screen%ae_rhor,"EM1",em1_qbz,xmpi_comm_self)

 ! Store the index of the q-point in the BZ for checking purpose.
 screen%fgg_qbz_idx = iq_bz

 DBG_EXIT("COLL")

end subroutine screen_rotate_iqbz
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_w0gemv
!! NAME
!! screen_w0gemv
!!
!! FUNCTION
!!  Perform the matrix multiplication W x vector in reciprocal space.
!!
!! INPUTS
!!  in_npw=Number of G vectors in in_ket
!!  nspinor=Number of spinorial components.
!!  in_ket(in_npw)= |\phi> in reciprocal space.
!!  trans= On entry, TRANS specifies the operation to be performed as follows:
!!  TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!!  TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!!  TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y.
!!
!! OUTPUT
!!   out_ket(in_npw)= W |\phi\> in reciprocal space.
!!   ZGEMV  performs one of the matrix-vector operations
!!   *
!!   *     y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
!!   *
!!   *     y := alpha*A**H*x + beta*y,
!!   *
!!   *  where alpha and beta are scalars, x and y are vectors and A is an m by n matrix.
!!
!! SOURCE

subroutine screen_w0gemv(screen, trans, in_npw, nspinor, only_diago, alpha, beta, in_ket, out_ket)

!Arguments ------------------------------------
!scalars
 class(screen_t),intent(in) :: screen
 integer,intent(in) :: in_npw,nspinor
 complex(gwpc),intent(in) :: alpha,beta
 logical,intent(in) :: only_diago
 character(len=*),intent(in) ::  trans
!arrays
 complex(gwpc),intent(in) :: in_ket(in_npw*nspinor)
 complex(gwpc),intent(out) :: out_ket(in_npw*nspinor)

!Local variables-------------------------------
!scalars
 integer :: ig,lda
!arrays
 complex(gwpc),ABI_CONTIGUOUS pointer :: em1_qbz(:,:)

! *************************************************************************

 lda = screen%npw; em1_qbz => screen%Fgg_qbz%mat(:,:,1)

 if (.not.only_diago) then
   call xgemv(trans,in_npw,in_npw,alpha,em1_qbz,lda,in_ket,1,beta,out_ket,1)

 else
   if (beta /= czero_gw) then
     if (firstchar(trans, (/"C"/))) then
       do ig=1,in_npw
         out_ket(ig) = alpha * CONJG(em1_qbz(ig,ig)) * in_ket(ig) + beta * out_ket(ig)
       end do
     else if (firstchar(trans, (/"N","T"/))) then
       do ig=1,in_npw
         out_ket(ig) = alpha * em1_qbz(ig,ig) * in_ket(ig) + beta * out_ket(ig)
       end do
     else
       ABI_ERROR(sjoin("Wrong trans:", trans))
     end if

   else
     ! beta == 0
     if (firstchar(trans, (/"C"/)) ) then
       do ig=1,in_npw
         out_ket(ig) = alpha * CONJG(em1_qbz(ig,ig)) * in_ket(ig)
       end do
     else if (firstchar(trans, (/"N","T"/))) then
       do ig=1,in_npw
         out_ket(ig) = alpha * em1_qbz(ig,ig) * in_ket(ig)
       end do
     else
       ABI_ERROR(sjoin("Wrong trans:", trans))
     end if
   end if
 end if

end subroutine screen_w0gemv
!!***

!----------------------------------------------------------------------

!!****f* m_screen/screen_calc_sigc
!! NAME
!! screen_calc_sigc
!!
!! FUNCTION
!!
!! INPUTS
!!  npw_c=Number of G vectors in in_ket
!!  nspinor=Number of spinorial components.
!!  in_ket(npw_c)= |\phi> in reciprocal space.
!!  trans= On entry, TRANS specifies the operation to be performed as follows:
!!      TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!!      TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!!      TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!!
!! OUTPUT
!!
!! SOURCE

subroutine screen_calc_sigc(screen, trans, nomega, omegame0i, theta_mu_minus_e0i, zcut, &
                            nspinor, npw_x, npw_c, rhotwgp, out_ket, sigcme)

!Arguments ------------------------------------
!scalars
 class(screen_t),intent(in) :: screen
 character(len=*),intent(in) ::  trans
 integer,intent(in) :: nomega, nspinor, npw_x, npw_c
 real(dp),intent(in) :: theta_mu_minus_e0i, zcut
!arrays
 real(dp),intent(in) :: omegame0i(nomega)
 complex(gwpc),intent(in) :: rhotwgp(npw_x*nspinor)
 complex(gwpc),intent(inout) :: out_ket(npw_c*nspinor, nomega)
 complex(gwpc),intent(out) :: sigcme(nomega)

!Local variables-------------------------------
 complex(gwpc),allocatable :: botsq_conjg_transp(:,:), otq_transp(:,:)

! *************************************************************************

 out_ket = czero_gw

 select case (screen%info%wint_method)

 case (WINT_PPMODEL)
   ABI_CHECK_IGE(screen%has_ppmodel, 0, "has_ppmodel should be > 0")

   select case (trans)
   case ("N")
     associate (botsq => screen%ppm%bigomegatwsq_qbz_vals, &
                otq   => screen%ppm%omegatw_qbz_vals, &
                eig   => screen%ppm%eigpot_qbz_vals)

       call screen%ppm%calc_sigc(nspinor, npw_c, nomega, rhotwgp, botsq, otq, &
                                 omegame0i, zcut, theta_mu_minus_e0i, eig, npw_x, out_ket, sigcme)

     end associate

   case ("T")
     associate (ppm => screen%ppm, &
                botsq => screen%ppm%bigomegatwsq_qbz_vals, &
                otq   => screen%ppm%omegatw_qbz_vals, &
                eig   => screen%ppm%eigpot_qbz_vals)

     ABI_MALLOC(botsq_conjg_transp,(PPm%dm2_botsq,npw_c))
     botsq_conjg_transp=TRANSPOSE(botsq) ! Keep these two lines separated, otherwise gfortran messes up
     !botsq_conjg_transp=CONJG(botsq_conjg_transp)
     ABI_MALLOC(otq_transp,(ppm%dm2_otq, ppm%npwc))
     otq_transp=TRANSPOSE(otq)

     call screen%ppm%calc_sigc(nspinor, npw_c, nomega, rhotwgp, botsq_conjg_transp, otq_transp, &
                               omegame0i, zcut, theta_mu_minus_e0i, eig, npw_x, out_ket, sigcme)

     ABI_FREE(botsq_conjg_transp)
     ABI_FREE(otq_transp)
     end associate

   case default
    ABI_ERROR(sjoin("Invalid trans:", trans))
   end select

 case default
   ABI_ERROR(sjoin("Unsupported wint_method:", itoa(screen%info%wint_method)))
 end select

end subroutine screen_calc_sigc
!!***

!----------------------------------------------------------------------


!!****f* m_screen/em1_symmetrize_ip
!! NAME
!!  em1_symmetrize_ip
!!
!! FUNCTION
!!  Symmetrizes the two-point function in G-space. Symmetrization is done
!!  inplace thorugh an auxiliary work array of dimension (npw_c,npw_c)
!!
!! INPUTS
!!  nomega=All frequencies from 1 up to nomega are symmetrized.
!!  npw_c=Number of G vectors in the symmetrized matrix.
!!  Gsph<gsphere_t>=data related to the G-sphere
!!  Qmesh<kmesh_t>=Structure defining the q-mesh used for Er.
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required.
!!
!! SIDE EFFECTS
!!  epsm1(npw_c,npw_c,nomega)
!!   input:  filled with the matrix at the q-point that has to be symmetrized.
!!   output: symmetrised matrix.
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required
!!  to reconstruct the BZ.
!!
!!  * Remember the symmetry properties of E
!!    If q_bz=Sq_ibz+G0:
!!
!!    $ E_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau} E_{G1,G2)}(q)
!!
!!    The invariance under exchange of the real space position E(1,2) = E(2,1) leads to:
!!    $ E_{-G2,-G1}(-q) = E_{G1,G2)
!!
!! SOURCE

subroutine em1_symmetrize_ip(iq_bz, npw_c, nomega, Gsph, Qmesh, epsm1)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,nomega,npw_c
 type(gsphere_t),intent(in) :: Gsph
 type(kmesh_t),intent(in) :: Qmesh
!arrays
 complex(gwpc),intent(inout) :: epsm1(npw_c,npw_c,nomega)

!Local variables-------------------------------
!scalars
 integer :: iw,g1,g2,isg1,isg2,iq_ibz,itim_q,isym_q,ierr
 logical :: q_isirred
 complex(gwpc) :: phmsg1t,phmsg2t_star
 !character(len=500) :: msg
!arrays
 real(dp) :: qbz(3)
 complex(gwpc),allocatable :: work(:,:)

! *********************************************************************

 ! Get iq_ibz, and symmetries from iq_ibz.
 call qmesh%get_BZ_item(iq_bz,qbz,iq_ibz,isym_q,itim_q,isirred=q_isirred)

 if (q_isirred) RETURN ! Nothing to do

 !write(msg,'(a,f8.2,a)')" out of memory in work , requiring ",npw_c**2*gwpc*b2Mb," Mb"
 ABI_MALLOC_OR_DIE(work, (npw_c,npw_c), ierr)

!$OMP PARALLEL DO PRIVATE(isg2,isg1,phmsg1t,phmsg2t_star,work) IF (nomega > 1)
 do iw=1,nomega
   do g2=1,npw_c
     isg2 = Gsph%rottb(g2,itim_q,isym_q)
     phmsg2t_star = CONJG(Gsph%phmSGt(g2,isym_q))
     do g1=1,npw_c
       isg1 = Gsph%rottb(g1,itim_q,isym_q)
       phmsg1t = Gsph%phmSGt(g1,isym_q)
       work(isg1,isg2) = epsm1(g1,g2,iw) * phmsg1t * phmsg2t_star
     end do
   end do
   epsm1(:,:,iw) = work
 end do

 ABI_FREE(work)

 ! Account for time-reversal
 if (itim_q==2) then
!$OMP PARALLEL DO IF (nomega > 1)
   do iw=1,nomega
     call sqmat_itranspose(npw_c, epsm1(:,:,iw))
   end do
 end if

end subroutine em1_symmetrize_ip
!!***

!----------------------------------------------------------------------

!!****f* m_screen/em1_symmetrize_op
!! NAME
!!  em1_symmetrize_op
!!
!! FUNCTION
!!  Symmetrizes the two-point function in G-space. Symmetrization is done outofplace.
!!
!! INPUTS
!!  nomega=All frequencies from 1 up to nomega are symmetrized.
!!  npw_c=Number of G vectors in the symmetrized matrix.
!!  Gsph<gsphere_t>=data related to the G-sphere
!!  Qmesh<kmesh_t>=Structure defining the q-mesh used for Er.
!!  iq_bz=Index of the q-point in the BZ where epsilon^-1 is required.
!!  in_epsm1(npw_c,npw_c,nomega)
!!
!! OUTPUT
!!  out_epsm1(npw_c,npw_c,nomega)
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0 in the
!!  expression Sq = q+G0. Treating this case would require some changes in the G-sphere
!!  since we have to consider G-G0. The code however stops in sigma if a nonzero G0 is required
!!  to reconstruct the BZ.
!!
!!  * Remember the symmetry properties of E
!!    If q_bz=Sq_ibz+G0:
!!
!!    $ E_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau} E_{G1,G2)}(q)
!!
!!    The invariance under exchange of the real space position E(1,2) = E(2,1) leads to:
!!    $ E_{-G2,-G1}(-q) = E_{G1,G2)
!!
!! SOURCE

subroutine em1_symmetrize_op(iq_bz, npw_c, nomega, Gsph, Qmesh, in_epsm1, out_epsm1)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz,nomega,npw_c
 type(gsphere_t),target,intent(in) :: Gsph
 type(kmesh_t),intent(in) :: Qmesh
!arrays
 complex(gwpc),intent(in) :: in_epsm1(npw_c,npw_c,nomega)
 complex(gwpc),intent(out) :: out_epsm1(npw_c,npw_c,nomega)

!Local variables-------------------------------
!scalars
 integer :: iw,g1,g2,isg1,isg2,iq_ibz,itim_q,isym_q
 logical :: q_isirred
 complex(gwpc) :: phmsg1t,phmsg2t_star
!arrays
 real(dp) :: qbz(3)

! *********************************************************************

 ! Get iq_ibz, and symmetries from iq_ibz.
 call qmesh%get_BZ_item(iq_bz,qbz,iq_ibz,isym_q,itim_q,isirred=q_isirred)

 if (q_isirred) then
   out_epsm1 = in_epsm1; return
 end if

! grottb is a 1-1 mapping.
!$OMP PARALLEL DO PRIVATE(isg1,isg2,phmsg1t,phmsg2t_star) COLLAPSE(2) IF (nomega > 1)
 do iw=1,nomega
   do g2=1,npw_c
     isg2=Gsph%rottb(g2,itim_q,isym_q)
     phmsg2t_star = CONJG(Gsph%phmSGt(g2,isym_q))
     do g1=1,npw_c
       isg1=Gsph%rottb(g1,itim_q,isym_q)
       phmsg1t = Gsph%phmSGt(g1,isym_q)
       out_epsm1(isg1,isg2,iw) = in_epsm1(g1,g2,iw) * phmsg1t * phmsg2t_star
     end do
   end do
 end do

 ! Account for time-reversal
 if (itim_q==2) then
!$OMP PARALLEL DO IF (nomega > 1)
   do iw=1,nomega
     call sqmat_itranspose(npw_c, out_epsm1(:,:,iw))
   end do
 end if

end subroutine em1_symmetrize_op
!!***

!----------------------------------------------------------------------

end module m_screen
!!***
