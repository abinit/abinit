!!****m* ABINIT/m_libpaw_libxc_funcs
!! NAME
!!  m_libpaw_libxc_funcs
!!
!! FUNCTION
!!  Module containing interfaces to the LibXC library, for exchange
!!  correlation potentials and energies.
!!
!! COPYRIGHT
!! Copyright (C) 2015-2020 ABINIT group (MO, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!!
!! NOTES
!!  This file comes directly from m_libpaw_libxc.F90 module delivered with ABINIT.
!!  It defines a structured datatype (libpaw_libxc_type) and associated methods
!!  to initialize/finalize it and get properties from it.
!!  * It uses by default a global variable (paw_xc_global, libpaw_libxc_type) which has
!!    to be initialized/finalized with the libpaw_libxc_init and libpaw_libxc_end methods.
!!  * It is also possible to define a local (private) variable of type libpaw_libxc_type.
!!    For that, the different methods have to be called with an extra optional
!!    argument (called xc_funcs in this example):
!!    !!!!! call libpaw_libxc_init(ixc,nspden,xc_funcs)
!!    !!!!! call libpaw_libxc_end(xc_funcs)
!!
!! SOURCE

!Need iso C bindings provided by the compiler
#define LIBPAW_ISO_C_BINDING 1

#include "libpaw.h"

module m_libpaw_libxc_funcs

 USE_DEFS
 USE_MSG_HANDLING
 USE_MEMORY_PROFILING

!ISO C bindings are mandatory
#ifdef LIBPAW_ISO_C_BINDING
 use iso_c_binding
#endif

 implicit none
 private

!Public functions
 public :: libpaw_libxc_check              ! Check if the code has been compiled with libXC
 public :: libpaw_libxc_init               ! Initialize the desired XC functional, from libXC
 public :: libpaw_libxc_end                ! End usage of libXC functional
 public :: libpaw_libxc_fullname           ! Return full name of the XC functional
 public :: libpaw_libxc_getrefs            ! Get references of a XC functional
 public :: libpaw_libxc_getid              ! Return identifer of a XC functional, from its name
 public :: libpaw_libxc_family_from_id     ! Retrieve family of a XC functional, from its id
 public :: libpaw_libxc_ixc                ! The value of ixc used to initialize the XC functionals
 public :: libpaw_libxc_getvxc             ! Return XC potential and energy, from input density
 public :: libpaw_libxc_isgga              ! Return TRUE if the XC functional is GGA or meta-GGA
 public :: libpaw_libxc_ismgga             ! Return TRUE if the XC functional is meta-GGA
 public :: libpaw_libxc_needs_laplacian    ! Return TRUE if functional uses LAPLACIAN
 public :: libpaw_libxc_is_hybrid          ! Return TRUE if the XC functional is hybrid
 public :: libpaw_libxc_is_hybrid_from_id  ! Return TRUE if a XC functional is hybrid, from its id
 public :: libpaw_libxc_has_kxc            ! Return TRUE if Kxc (3rd der) is available for the XC functional
 public :: libpaw_libxc_has_k3xc           ! Return TRUE if K3xc (4th der) is available for the XC functional
 public :: libpaw_libxc_nspin              ! The number of spin components for the XC functionals
 public :: libpaw_libxc_get_hybridparams   ! Retrieve parameter(s) of a hybrid functional
 public :: libpaw_libxc_set_hybridparams   ! Change parameter(s) of a hybrid functional
 public :: libpaw_libxc_gga_from_hybrid    ! Return the id of the XC-GGA used for the hybrid

!Private functions
 private :: libpaw_libxc_constants_load    ! Load libXC constants from C headers
 private :: libpaw_libxc_compute_tb09      ! Compute c parameter for Tran-Blaha 2009 functional
#ifdef LIBPAW_ISO_C_BINDING
 private :: char_f_to_c                    ! Convert a string from Fortran to C
 private :: char_c_to_f                    ! Convert a string from C to Fortran
#endif

!Public constants (use libpaw_libxc_constants_load to init them)
 integer,public,save :: LIBPAW_XC_FAMILY_UNKNOWN       = -1
 integer,public,save :: LIBPAW_XC_FAMILY_LDA           =  1
 integer,public,save :: LIBPAW_XC_FAMILY_GGA           =  2
 integer,public,save :: LIBPAW_XC_FAMILY_MGGA          =  4
 integer,public,save :: LIBPAW_XC_FAMILY_LCA           =  8
 integer,public,save :: LIBPAW_XC_FAMILY_OEP           = 16
 integer,public,save :: LIBPAW_XC_FAMILY_HYB_GGA       = 32
 integer,public,save :: LIBPAW_XC_FAMILY_HYB_MGGA      = 64
 integer,public,save :: LIBPAW_XC_FLAGS_HAVE_EXC       =  1
 integer,public,save :: LIBPAW_XC_FLAGS_HAVE_VXC       =  2
 integer,public,save :: LIBPAW_XC_FLAGS_HAVE_FXC       =  4
 integer,public,save :: LIBPAW_XC_FLAGS_HAVE_KXC       =  8
 integer,public,save :: LIBPAW_XC_FLAGS_HAVE_LXC       = 16
 integer,public,save :: LIBPAW_XC_FLAGS_NEEDS_LAPLACIAN= 32768
 integer,public,save :: LIBPAW_XC_EXCHANGE             =  0
 integer,public,save :: LIBPAW_XC_CORRELATION          =  1
 integer,public,save :: LIBPAW_XC_EXCHANGE_CORRELATION =  2
 integer,public,save :: LIBPAW_XC_KINETIC              =  3
 integer,public,save :: LIBPAW_XC_HYB_NONE             =  0
 integer,public,save :: LIBPAW_XC_HYB_FOCK             =  1
 integer,public,save :: LIBPAW_XC_HYB_PT2              =  2
 integer,public,save :: LIBPAW_XC_HYB_ERF_SR           =  4
 integer,public,save :: LIBPAW_XC_HYB_YUKAWA_SR        =  8
 integer,public,save :: LIBPAW_XC_HYB_GAUSSIAN_SR      = 16
 integer,public,save :: LIBPAW_XC_HYB_SEMILOCAL        =  0
 integer,public,save :: LIBPAW_XC_HYB_HYBRID           =  1
 integer,public,save :: LIBPAW_XC_HYB_CAM              =  2
 integer,public,save :: LIBPAW_XC_HYB_CAMY             =  3
 integer,public,save :: LIBPAW_XC_HYB_CAMG             =  4
 integer,public,save :: LIBPAW_XC_HYB_DOUBLE_HYBRID    =  5
 integer,public,save :: LIBPAW_XC_HYB_MIXTURE          = 32768
 integer,public,save :: LIBPAW_XC_SINGLE_PRECISION     =  0
 logical,private,save :: libpaw_xc_constants_initialized=.false.

!XC functional public type
 type,public :: libpaw_libxc_type
   integer  :: id              ! identifier
   integer  :: family          ! LDA, GGA, etc.
   integer  :: kind            ! EXCHANGE, CORRELATION, etc.
   integer  :: nspin           ! # of spin components
   integer  :: abi_ixc         ! Abinit IXC id for this functional
   logical  :: has_exc         ! TRUE is exc is available for the functional
   logical  :: has_vxc         ! TRUE is vxc is available for the functional
   logical  :: has_fxc         ! TRUE is fxc is available for the functional
   logical  :: has_kxc         ! TRUE is kxc is available for the functional
   logical  :: needs_laplacian ! TRUE is functional needs laplacian of density
   logical  :: is_hybrid       ! TRUE is functional is a hybrid functional
   real(dp) :: hyb_mixing      ! Hybrid functional: mixing factor of Fock contribution (default=0)
   real(dp) :: hyb_mixing_sr   ! Hybrid functional: mixing factor of SR Fock contribution (default=0)
   real(dp) :: hyb_range       ! Range (for separation) for a hybrid functional (default=0)
   real(dp) :: xc_tb09_c       ! Special TB09 functional parameter
#ifdef LIBPAW_ISO_C_BINDING
   type(C_PTR),pointer :: conf => null() ! C pointer to the functional itself
#endif
 end type libpaw_libxc_type

!----------------------------------------------------------------------

!Private global XC functional
 type(libpaw_libxc_type),target,save :: paw_xc_global(2)

!----------------------------------------------------------------------

!Interfaces for C bindings
#ifdef LIBPAW_ISO_C_BINDING
 interface
   integer(C_INT) function xc_func_init(xc_func,functional,nspin) bind(C,name="xc_func_init")
     use iso_c_binding, only : C_INT,C_PTR
     integer(C_INT),value :: functional,nspin
     type(C_PTR) :: xc_func
   end function xc_func_init
 end interface
!
 interface
   subroutine xc_func_end(xc_func) bind(C,name="xc_func_end")
     use iso_c_binding, only : C_PTR
     type(C_PTR) :: xc_func
   end subroutine xc_func_end
 end interface
!
 interface
   integer(C_INT) function xc_functional_get_number(name) &
&                          bind(C,name="xc_functional_get_number")
     use iso_c_binding, only : C_INT,C_PTR
     type(C_PTR),value :: name
   end function xc_functional_get_number
 end interface
!
 interface
   type(C_PTR) function xc_functional_get_name(number) &
&                       bind(C,name="xc_functional_get_name")
     use iso_c_binding, only : C_INT,C_PTR
     integer(C_INT),value :: number
   end function xc_functional_get_name
 end interface
!
 interface
   integer(C_INT) function xc_family_from_id(id,family,number) &
&                          bind(C,name="xc_family_from_id")
     use iso_c_binding, only : C_INT,C_PTR
     integer(C_INT),value :: id
     type(C_PTR),value :: family,number
   end function xc_family_from_id
 end interface
!
 interface
   subroutine xc_hyb_cam_coef(xc_func,omega,alpha,beta) &
&             bind(C,name="xc_hyb_cam_coef")
     use iso_c_binding, only : C_DOUBLE,C_PTR
     real(C_DOUBLE) :: omega,alpha,beta
     type(C_PTR) :: xc_func
   end subroutine xc_hyb_cam_coef
 end interface
!
 interface
   subroutine libpaw_xc_get_lda(xc_func,np,rho,zk,vrho,v2rho2,v3rho3) &
&             bind(C,name="libpaw_xc_get_lda")
     use iso_c_binding, only : C_INT,C_PTR
     integer(C_INT),value :: np
     type(C_PTR),value :: rho,zk,vrho,v2rho2,v3rho3
     type(C_PTR) :: xc_func
   end subroutine libpaw_xc_get_lda
 end interface
!
 interface
   subroutine libpaw_xc_get_gga(xc_func,np,rho,sigma,zk,vrho,vsigma, &
&             v2rho2,v2rhosigma,v2sigma2,v3rho3,v3rho2sigma,v3rhosigma2,v3sigma3) &
&             bind(C,name="libpaw_xc_get_gga")
     use iso_c_binding, only : C_INT,C_PTR
     integer(C_INT),value :: np
     type(C_PTR),value :: rho,sigma,zk,vrho,vsigma,v2rho2,v2rhosigma,v2sigma2, &
&                         v3rho3,v3rho2sigma,v3rhosigma2,v3sigma3
     type(C_PTR) :: xc_func
   end subroutine libpaw_xc_get_gga
 end interface
!
 interface
   subroutine libpaw_xc_get_mgga(xc_func,np,rho,sigma,lapl,tau,zk,vrho,vsigma,vlapl,vtau, &
&             v2rho2,v2rhosigma,v2rholapl,v2rhotau,v2sigma2,v2sigmalapl, &
&             v2sigmatau,v2lapl2,v2lapltau,v2tau2) &
&             bind(C,name="libpaw_xc_get_mgga")
     use iso_c_binding, only : C_INT,C_PTR
     integer(C_INT),value :: np
     type(C_PTR),value :: rho,sigma,lapl,tau,zk,vrho,vsigma,vlapl,vtau, &
&                         v2rho2,v2sigma2,v2lapl2,v2tau2,v2rhosigma,v2rholapl,v2rhotau, &
&                         v2sigmalapl,v2sigmatau,v2lapltau
     type(C_PTR) :: xc_func
   end subroutine libpaw_xc_get_mgga
 end interface
!
 interface
   subroutine libpaw_xc_func_set_params(xc_func,params,n_params) bind(C)
     use iso_c_binding, only : C_INT,C_DOUBLE,C_PTR
     integer(C_INT),value :: n_params
     real(C_DOUBLE) :: params(*)
     type(C_PTR) :: xc_func
   end subroutine libpaw_xc_func_set_params
 end interface
!
 interface
   subroutine libpaw_xc_func_set_density_threshold(xc_func,dens_threshold) bind(C)
     use iso_c_binding, only : C_DOUBLE,C_PTR
     real(C_DOUBLE) :: dens_threshold
     type(C_PTR) :: xc_func
   end subroutine libpaw_xc_func_set_density_threshold
 end interface
!
 interface
   integer(C_INT) function libpaw_xc_func_is_hybrid_from_id(func_id) bind(C)
     use iso_c_binding, only : C_INT
     integer(C_INT),value :: func_id
   end function libpaw_xc_func_is_hybrid_from_id
 end interface
!
 interface
   subroutine libpaw_xc_get_singleprecision_constant(xc_cst_singleprecision) &
&             bind(C,name="libpaw_xc_get_singleprecision_constant")
     use iso_c_binding, only : C_INT
     integer(C_INT) :: xc_cst_singleprecision
   end subroutine libpaw_xc_get_singleprecision_constant
 end interface
!
 interface
   subroutine libpaw_xc_get_family_constants(xc_cst_unknown,xc_cst_lda,xc_cst_gga, &
&             xc_cst_mgga,xc_cst_lca,xc_cst_oep,xc_cst_hyb_gga,xc_cst_hyb_mgga) &
&             bind(C,name="libpaw_xc_get_family_constants")
     use iso_c_binding, only : C_INT
     integer(C_INT) :: xc_cst_unknown,xc_cst_lda,xc_cst_gga,xc_cst_mgga, &
&                      xc_cst_lca,xc_cst_oep,xc_cst_hyb_gga,xc_cst_hyb_mgga
   end subroutine libpaw_xc_get_family_constants
 end interface
!
 interface
   subroutine libpaw_xc_get_flags_constants(xc_cst_flags_have_exc,xc_cst_flags_have_vxc, &
              xc_cst_flags_have_fxc,xc_cst_flags_have_kxc,xc_cst_flags_have_lxc, &
&             xc_cst_flags_needs_lapl) &
&             bind(C,name="libpaw_xc_get_flags_constants")
     use iso_c_binding, only : C_INT
     integer(C_INT) :: xc_cst_flags_have_exc,xc_cst_flags_have_vxc,xc_cst_flags_have_fxc, &
&                      xc_cst_flags_have_kxc,xc_cst_flags_have_lxc,xc_cst_flags_needs_lapl
   end subroutine libpaw_xc_get_flags_constants
 end interface
!
 interface
   subroutine libpaw_xc_get_kind_constants(xc_cst_exchange,xc_cst_correlation, &
&             xc_cst_exchange_correlation,xc_cst_kinetic) &
&             bind(C,name="libpaw_xc_get_kind_constants")
     use iso_c_binding, only : C_INT
     integer(C_INT) :: xc_cst_exchange,xc_cst_correlation, &
&                      xc_cst_exchange_correlation,xc_cst_kinetic
   end subroutine libpaw_xc_get_kind_constants
 end interface
!
 interface
   subroutine libpaw_xc_get_hybrid_constants(xc_cst_hyb_none, &
              xc_cst_hyb_fock,xc_cst_hyb_pt2,xc_cst_hyb_erf_sr,xc_cst_hyb_yukawa_sr, &
              xc_cst_hyb_gaussian_sr,xc_cst_hyb_semilocal, xc_cst_hyb_hybrid,xc_cst_hyb_cam, &
              xc_cst_hyb_camy,xc_cst_hyb_camg,xc_cst_hyb_double_hybrid, &
              xc_cst_hyb_mixture) bind(C,name="libpaw_xc_get_hybrid_constants")
     use iso_c_binding, only : C_INT
     integer(C_INT) :: xc_cst_hyb_none, xc_cst_hyb_fock,xc_cst_hyb_pt2, xc_cst_hyb_erf_sr, &
                       xc_cst_hyb_yukawa_sr,xc_cst_hyb_gaussian_sr,xc_cst_hyb_semilocal, &
                       xc_cst_hyb_hybrid,xc_cst_hyb_cam,xc_cst_hyb_camy,xc_cst_hyb_camg, &
                       xc_cst_hyb_double_hybrid,xc_cst_hyb_mixture
   end subroutine libpaw_xc_get_hybrid_constants
 end interface
!
 interface
   type(C_PTR) function libpaw_xc_func_type_malloc() &
&                       bind(C,name="libpaw_xc_func_type_malloc")
     use iso_c_binding, only : C_PTR
   end function libpaw_xc_func_type_malloc
 end interface
!
 interface
   subroutine libpaw_xc_func_type_free(xc_func) &
&             bind(C,name="libpaw_xc_func_type_free")
     use iso_c_binding, only : C_PTR
     type(C_PTR) :: xc_func
   end subroutine libpaw_xc_func_type_free
 end interface
!
 interface
   type(C_PTR) function libpaw_xc_get_info_name(xc_func) &
&                       bind(C,name="libpaw_xc_get_info_name")
     use iso_c_binding, only : C_PTR
     type(C_PTR) :: xc_func
   end function libpaw_xc_get_info_name
 end interface
!
 interface
   type(C_PTR) function libpaw_xc_get_info_refs(xc_func,iref) &
&                       bind(C,name="libpaw_xc_get_info_refs")
     use iso_c_binding, only : C_INT,C_PTR
     type(C_PTR) :: xc_func
     integer(C_INT) :: iref
   end function libpaw_xc_get_info_refs
 end interface
!
 interface
   integer(C_INT) function libpaw_xc_get_info_flags(xc_func) &
&                          bind(C,name="libpaw_xc_get_info_flags")
     use iso_c_binding, only : C_INT,C_PTR
     type(C_PTR) :: xc_func
   end function libpaw_xc_get_info_flags
 end interface
!
 interface
   integer(C_INT) function libpaw_xc_get_info_kind(xc_func) &
&                          bind(C,name="libpaw_xc_get_info_kind")
     use iso_c_binding, only : C_INT,C_PTR
     type(C_PTR) :: xc_func
   end function libpaw_xc_get_info_kind
 end interface
#endif

contains
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_constants_load
!! NAME
!!  libpaw_libxc_constants_load
!!
!! FUNCTION
!!  Load libXC constants from C headers
!!
!! PARENTS
!!      m_libpaw_libxc
!!
!! CHILDREN
!!
!! SOURCE

 subroutine libpaw_libxc_constants_load()

!Local variables-------------------------------
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 integer(C_INT) :: i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13
#endif

! *************************************************************************

#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
  call libpaw_xc_get_singleprecision_constant(i1)
  LIBPAW_XC_SINGLE_PRECISION     = int(i1)
  call libpaw_xc_get_family_constants(i1,i2,i3,i4,i5,i6,i7,i8)
  LIBPAW_XC_FAMILY_UNKNOWN       = int(i1)
  LIBPAW_XC_FAMILY_LDA           = int(i2)
  LIBPAW_XC_FAMILY_GGA           = int(i3)
  LIBPAW_XC_FAMILY_MGGA          = int(i4)
  LIBPAW_XC_FAMILY_LCA           = int(i5)
  LIBPAW_XC_FAMILY_OEP           = int(i6)
  LIBPAW_XC_FAMILY_HYB_GGA       = int(i7)
  LIBPAW_XC_FAMILY_HYB_MGGA      = int(i8)
  call libpaw_xc_get_flags_constants(i1,i2,i3,i4,i5,i6)
  LIBPAW_XC_FLAGS_HAVE_EXC       = int(i1)
  LIBPAW_XC_FLAGS_HAVE_VXC       = int(i2)
  LIBPAW_XC_FLAGS_HAVE_FXC       = int(i3)
  LIBPAW_XC_FLAGS_HAVE_KXC       = int(i4)
  LIBPAW_XC_FLAGS_HAVE_LXC       = int(i5)
  LIBPAW_XC_FLAGS_NEEDS_LAPLACIAN= int(i6)
  call libpaw_xc_get_kind_constants(i1,i2,i3,i4)
  LIBPAW_XC_EXCHANGE             = int(i1)
  LIBPAW_XC_CORRELATION          = int(i2)
  LIBPAW_XC_EXCHANGE_CORRELATION = int(i3)
  LIBPAW_XC_KINETIC              = int(i4)
  call libpaw_xc_get_hybrid_constants(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13)
  LIBPAW_XC_HYB_NONE             = int(i1)
  LIBPAW_XC_HYB_FOCK             = int(i2)
  LIBPAW_XC_HYB_PT2              = int(i3)
  LIBPAW_XC_HYB_ERF_SR           = int(i4)
  LIBPAW_XC_HYB_YUKAWA_SR        = int(i5)
  LIBPAW_XC_HYB_GAUSSIAN_SR      = int(i6)
  LIBPAW_XC_HYB_SEMILOCAL        = int(i7)
  LIBPAW_XC_HYB_HYBRID           = int(i8)
  LIBPAW_XC_HYB_CAM              = int(i9)
  LIBPAW_XC_HYB_CAMY             = int(i10)
  LIBPAW_XC_HYB_CAMG             = int(i11)
  LIBPAW_XC_HYB_DOUBLE_HYBRID    = int(i12)
  LIBPAW_XC_HYB_MIXTURE          = int(i13)
  libpaw_xc_constants_initialized=.true.
#endif

 end subroutine libpaw_libxc_constants_load
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_check
!! NAME
!!  libpaw_libxc_check
!!
!! FUNCTION
!!  Check if the code has been compiled with an usable version of libXC
!!
!! INPUTS
!! [stop_if_error]=optional flag; if TRUE the code stops if libXC is not correctly used
!!
!! SOURCE

 function libpaw_libxc_check(stop_if_error)

!Arguments ------------------------------------
 logical :: libpaw_libxc_check
 logical,intent(in),optional :: stop_if_error
!Local variables-------------------------------
 character(len=100) :: msg

! *************************************************************************

#if defined LIBPAW_HAVE_LIBXC
#if defined LIBPAW_ISO_C_BINDING
 if (.not.libpaw_xc_constants_initialized) call libpaw_libxc_constants_load()
 if (LIBPAW_XC_SINGLE_PRECISION==1) then
   libpaw_libxc_check=.false.
   msg='LibXC should be compiled with double precision!'
 end if
#else
 libpaw_libxc_check=.false.
 msg='LibXC cannot be used without ISO_C_BINDING support by the Fortran compiler!'
#endif
#else
 libpaw_libxc_check=.false.
 msg='LibPAW was not compiled with LibXC support.'
#endif

 if (present(stop_if_error)) then
   if (stop_if_error.and.trim(msg)/="") then
     LIBPAW_ERROR(msg)
   end if
 end if

 end function libpaw_libxc_check
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_init
!! NAME
!!  libpaw_libxc_init
!!
!! FUNCTION
!!  Initialize the desired XC functional, from LibXC.
!!  * Call the LibXC initializer
!!  * Fill preliminary fields in module structures.
!!
!! INPUTS
!! ixc=XC code for Abinit
!! nspden=number of spin-density components
!! [xc_tb09_c]=special argument for the Tran-Blaha 2009 functional
!!
!! SIDE EFFECTS
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional argument
!!                     XC functionals to initialize
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine libpaw_libxc_init(ixc,nspden,xc_functionals,&
&                             xc_tb09_c) ! optional argument


!Arguments ------------------------------------
 integer, intent(in) :: nspden
 integer, intent(in) :: ixc
 real(dp),intent(in),optional :: xc_tb09_c
 type(libpaw_libxc_type),intent(inout),optional,target :: xc_functionals(2)
!Local variables-------------------------------
 integer :: ii,nspden_eff
 character(len=500) :: msg
 type(libpaw_libxc_type),pointer :: xc_func
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 integer :: flags
 integer(C_INT) :: func_id_c,iref_c,npar_c,nspin_c,success_c
 real(C_DOUBLE) :: alpha_c,beta_c,omega_c,param_c(1)
 character(kind=C_CHAR,len=1),pointer :: strg_c
 type(C_PTR) :: func_ptr_c
#endif

! *************************************************************************

!Check libXC
 if (.not.libpaw_libxc_check(stop_if_error=.true.)) return
 if (.not.libpaw_xc_constants_initialized) call libpaw_libxc_constants_load()

 nspden_eff=min(nspden,2)

!Select XC functional(s) identifiers
 if (present(xc_functionals)) then
   xc_functionals(1)%id = -ixc/1000
   xc_functionals(2)%id = -ixc + (ixc/1000)*1000
 else
   paw_xc_global(1)%id = -ixc/1000
   paw_xc_global(2)%id = -ixc + (ixc/1000)*1000
 end if

 do ii = 1,2

!  Select XC functional
   if (present(xc_functionals)) then
     xc_func => xc_functionals(ii)
   else
     xc_func => paw_xc_global(ii)
   end if

   xc_func%abi_ixc=ixc !Save abinit value for reference

   xc_func%family=LIBPAW_XC_FAMILY_UNKNOWN
   xc_func%kind=-1
   xc_func%nspin=nspden_eff
   xc_func%has_exc=.false.
   xc_func%has_vxc=.false.
   xc_func%has_fxc=.false.
   xc_func%has_kxc=.false.
   xc_func%needs_laplacian=.false.
   xc_func%is_hybrid=.false.
   xc_func%hyb_mixing=zero
   xc_func%hyb_mixing_sr=zero
   xc_func%hyb_range=zero
   xc_func%xc_tb09_c=99.99_dp

   if (xc_func%id<=0) cycle

!  Get XC functional family
   xc_func%family=libpaw_libxc_family_from_id(xc_func%id)
   if (xc_func%family/=LIBPAW_XC_FAMILY_LDA .and. &
&      xc_func%family/=LIBPAW_XC_FAMILY_GGA .and. &
&      xc_func%family/=LIBPAW_XC_FAMILY_MGGA.and. &
&      xc_func%family/=LIBPAW_XC_FAMILY_HYB_GGA) then
     write(msg, '(a,i8,2a,i8,6a)' )&
&      'Invalid IXC = ',ixc,ch10,&
&      'The LibXC functional family ',xc_func%family,&
&      'is currently unsupported by LibPAW',ch10,&
&      '(-1 means the family is unknown to the LibXC itself)',ch10,&
&      'Please consult the LibXC documentation',ch10
     LIBPAW_ERROR(msg)
   end if

#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING

!  Allocate functional
   func_ptr_c=libpaw_xc_func_type_malloc()
   call c_f_pointer(func_ptr_c,xc_func%conf)

!  Initialize functional
   func_id_c=int(xc_func%id,kind=C_INT)
   nspin_c=int(nspden_eff,kind=C_INT)
   success_c=xc_func_init(xc_func%conf,func_id_c,nspin_c)
   if (success_c/=0) then
     msg='Error in libXC functional initialization!'
     LIBPAW_ERROR(msg)
   end if

!  Special treatment for LDA_C_XALPHA functional
   if (xc_func%id==libpaw_libxc_getid('XC_LDA_C_XALPHA')) then
     param_c(1)=real(zero,kind=C_DOUBLE);npar_c=int(1,kind=C_INT)
     call libpaw_xc_func_set_params(xc_func%conf,param_c,npar_c)
   end if

!  Special treatment for XC_MGGA_X_TB09  functional
   if (xc_func%id==libpaw_libxc_getid('XC_MGGA_X_TB09')) then
     if (.not.present(xc_tb09_c)) then
       msg='xc_tb09_c argument is mandatory for TB09 functional!'
       LIBPAW_BUG(msg)
     end if
     xc_func%xc_tb09_c=xc_tb09_c
   end if

!  Get functional kind
   xc_func%kind=int(libpaw_xc_get_info_kind(xc_func%conf))

!  Get functional flags
   flags=int(libpaw_xc_get_info_flags(xc_func%conf))
   xc_func%has_exc=(iand(flags,LIBPAW_XC_FLAGS_HAVE_EXC)>0)
   xc_func%has_vxc=(iand(flags,LIBPAW_XC_FLAGS_HAVE_VXC)>0)
   xc_func%has_fxc=(iand(flags,LIBPAW_XC_FLAGS_HAVE_FXC)>0)
   xc_func%has_kxc=(iand(flags,LIBPAW_XC_FLAGS_HAVE_KXC)>0)

!  Retrieve parameters for metaGGA functionals
   if (xc_func%family==LIBPAW_XC_FAMILY_MGGA.or. &
&      xc_func%family==LIBPAW_XC_FAMILY_HYB_MGGA) then
     xc_func%needs_laplacian=(iand(flags,LIBPAW_XC_FLAGS_NEEDS_LAPLACIAN)>0)
   end if

!  Retrieve parameters for hybrid functionals
   xc_func%is_hybrid=(libpaw_xc_func_is_hybrid_from_id(xc_func%id)==1)
   if (xc_func%is_hybrid) then
     call xc_hyb_cam_coef(xc_func%conf,omega_c,alpha_c,beta_c)
     xc_func%hyb_mixing=real(alpha_c,kind=dp)
     xc_func%hyb_mixing_sr=real(beta_c,kind=dp)
     xc_func%hyb_range=real(omega_c,kind=dp)
   end if

!  Dump functional information
   call c_f_pointer(libpaw_xc_get_info_name(xc_func%conf),strg_c)
   call char_c_to_f(strg_c,msg)
   call wrtout(std_out,msg,'COLL')
   iref_c=0
   do while (iref_c>=0)
     call c_f_pointer(libpaw_xc_get_info_refs(xc_func%conf,iref_c),strg_c)
     if (associated(strg_c)) then
       call char_c_to_f(strg_c,msg)
       call wrtout(std_out,msg,'COLL')
       iref_c=iref_c+1
     else
       iref_c=-1
     end if
   end do

#else
   if (.False.) write(std_out,*)xc_tb09_c
#endif

 end do

end subroutine libpaw_libxc_init
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_end
!! NAME
!!  libpaw_libxc_end
!!
!! FUNCTION
!!  End usage of LibXC functional. Call LibXC end function,
!!  and deallocate module contents.
!!
!! SIDE EFFECTS
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional argument
!!                     XC functionals to initialize
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine libpaw_libxc_end(xc_functionals)

!Arguments ------------------------------------
 type(libpaw_libxc_type),intent(inout),optional,target :: xc_functionals(2)
!Local variables-------------------------------
 integer :: ii
 type(libpaw_libxc_type),pointer :: xc_func

! *************************************************************************

 do ii = 1,2

!  Select XC functional
   if (present(xc_functionals)) then
     xc_func => xc_functionals(ii)
   else
     xc_func => paw_xc_global(ii)
   end if

   if (xc_func%id <= 0) cycle
   xc_func%id=-1
   xc_func%family=-1
   xc_func%kind=-1
   xc_func%nspin=1
   xc_func%abi_ixc=huge(0)
   xc_func%has_exc=.false.
   xc_func%has_vxc=.false.
   xc_func%has_fxc=.false.
   xc_func%has_kxc=.false.
   xc_func%needs_laplacian=.false.
   xc_func%is_hybrid=.false.
   xc_func%hyb_mixing_sr=zero
   xc_func%hyb_range=zero
   xc_func%xc_tb09_c=99.99_dp
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
   if (associated(xc_func%conf)) then
     call xc_func_end(xc_func%conf)
     call libpaw_xc_func_type_free(c_loc(xc_func%conf))
   end if
#endif

 end do

 end subroutine libpaw_libxc_end
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_fullname
!! NAME
!!  libpaw_libxc_fullname
!!
!! FUNCTION
!!  Return full name of the XC functional
!!
!! INPUTS
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional argument
!!                     XC functionals to initialize
!!
!! SOURCE

 function libpaw_libxc_fullname(xc_functionals)

!Arguments ------------------------------------
 character(len=100) :: libpaw_libxc_fullname
 type(libpaw_libxc_type),intent(in),optional,target :: xc_functionals(2)
!Local variables-------------------------------
 integer :: nxc
 type(libpaw_libxc_type),pointer :: xc_funcs(:)
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 character(len=100) :: xcname
 character(kind=C_CHAR,len=1),pointer :: strg_c
#endif

! *************************************************************************

 libpaw_libxc_fullname='No XC functional'

 if (present(xc_functionals)) then
   xc_funcs => xc_functionals
 else
   xc_funcs => paw_xc_global
 end if

 nxc=size(xc_funcs)
 if (nxc<1) return

#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 if (nxc<2) then
   if (xc_funcs(1)%id /= 0) then
     call c_f_pointer(xc_functional_get_name(xc_funcs(1)%id),strg_c)
     call char_c_to_f(strg_c,libpaw_libxc_fullname)
   end if
 else if (xc_funcs(1)%id <= 0) then
   if (xc_funcs(2)%id /= 0) then
     call c_f_pointer(xc_functional_get_name(xc_funcs(2)%id),strg_c)
     call char_c_to_f(strg_c,libpaw_libxc_fullname)
   end if
 else if (xc_funcs(2)%id <= 0) then
   if (xc_funcs(1)%id /= 0) then
     call c_f_pointer(xc_functional_get_name(xc_funcs(1)%id),strg_c)
     call char_c_to_f(strg_c,libpaw_libxc_fullname)
   end if
 else
   call c_f_pointer(xc_functional_get_name(xc_funcs(1)%id),strg_c)
   call char_c_to_f(strg_c,libpaw_libxc_fullname)
   call c_f_pointer(xc_functional_get_name(xc_funcs(2)%id),strg_c)
   call char_c_to_f(strg_c,xcname)
   libpaw_libxc_fullname=trim(libpaw_libxc_fullname)//'+'//trim(xcname)
 end if
 libpaw_libxc_fullname=trim(libpaw_libxc_fullname)
#endif

end function libpaw_libxc_fullname
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_getrefs
!! NAME
!!  libpaw_libxc_getrefs
!!
!! FUNCTION
!!  Return the reference(s) of a XC functional
!!
!! INPUTS
!!  xc_functional=<type(libpaw_libxc_type)>, handle for XC functional
!!
!! OUTPUT
!!  xcrefs(:)= references(s) of the functional
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine libpaw_libxc_getrefs(xcrefs,xc_functional)

!Arguments ------------------------------------
 character(len=*),intent(out) :: xcrefs(:)
 type(libpaw_libxc_type),intent(in) :: xc_functional

!Local variables-------------------------------
#if defined LIBPAW_HAVE_LIBXC && defined HAVE_FC_ISO_C_BINDING
 integer(C_INT) :: iref_c
 character(kind=C_CHAR,len=1),pointer :: strg_c
#endif

! *************************************************************************

 xcrefs(:)=''

#if defined LIBPAW_HAVE_LIBXC && defined HAVE_FC_ISO_C_BINDING
 iref_c=0
 do while (iref_c>=0.and.iref_c<size(xcrefs))
   call c_f_pointer(libpaw_xc_get_info_refs(xc_functional%conf,iref_c),strg_c)
   if (associated(strg_c)) then
     call char_c_to_f(strg_c,xcrefs(iref_c+1))
     iref_c=iref_c+1
   else
     iref_c=-1
   end if
 end do
#else
 if (.False.) write(std_out,*)xc_functional%id
#endif

end subroutine libpaw_libxc_getrefs
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_family_from_id
!! NAME
!!  libpaw_libxc_family_from_id
!!
!! FUNCTION
!!  Return family of a XC functional from its id
!!
!! INPUTS
!!  xcid= id of a LibXC functional
!!
!! SOURCE

 function libpaw_libxc_family_from_id(xcid)

!Arguments ------------------------------------
 integer :: libpaw_libxc_family_from_id
 integer,intent(in) :: xcid
!Local variables-------------------------------
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 integer(C_INT) :: xcid_c
#endif

! *************************************************************************

#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 xcid_c=int(xcid,kind=C_INT)
 libpaw_libxc_family_from_id=int(xc_family_from_id(xcid_c,C_NULL_PTR,C_NULL_PTR))
#else
 libpaw_libxc_family_from_id=-1
 if (.false.) write(std_out,*) xcid
#endif

end function libpaw_libxc_family_from_id
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_getid
!! NAME
!!  libpaw_libxc_getid
!!
!! FUNCTION
!!  Return identifer of a XC functional from its name
!!  Return -1 if undefined
!!
!! INPUTS
!!  xcname= string containing the name of a XC functional
!!
!! SOURCE

 function libpaw_libxc_getid(xcname)

!Arguments ------------------------------------
 integer :: libpaw_libxc_getid
 character(len=*),intent(in) :: xcname
!Local variables-------------------------------
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 character(len=256) :: str
 character(kind=C_CHAR,len=1),target :: name_c(len_trim(xcname)+1)
 character(kind=C_CHAR,len=1),target :: name_c_xc(len_trim(xcname)-2)
 type(C_PTR) :: name_c_ptr
#endif

! *************************************************************************

#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 str=trim(xcname)
 if (xcname(1:3)=="XC_".or.xcname(1:3)=="xc_") then
   str=xcname(4:);name_c_xc=char_f_to_c(str)
   name_c_ptr=c_loc(name_c_xc)
 else
   name_c=char_f_to_c(str)
   name_c_ptr=c_loc(name_c)
 end if
 libpaw_libxc_getid=int(xc_functional_get_number(name_c_ptr))
#else
 libpaw_libxc_getid=-1
 if (.false.) write(std_out,*) xcname
#endif

end function libpaw_libxc_getid
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_ixc
!! NAME
!!  libpaw_libxc_ixc
!!
!! FUNCTION
!!  Return the value of ixc used to initialize the XC structure
!!
!! INPUTS
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional argument
!!                     XC functionals to initialize
!!
!! SOURCE

 function libpaw_libxc_ixc(xc_functionals)

!Arguments ------------------------------------
 integer :: libpaw_libxc_ixc
 type(libpaw_libxc_type),intent(in),optional :: xc_functionals(2)

! *************************************************************************

 if (present(xc_functionals)) then
   libpaw_libxc_ixc=xc_functionals(1)%abi_ixc
 else
   libpaw_libxc_ixc=paw_xc_global(1)%abi_ixc
 end if

end function libpaw_libxc_ixc
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_isgga
!! NAME
!!  libpaw_libxc_isgga
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  is a GGA or not
!!
!! INPUTS
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional argument
!!                     XC functionals to initialize
!!
!! SOURCE

 function libpaw_libxc_isgga(xc_functionals)

!Arguments ------------------------------------
 logical :: libpaw_libxc_isgga
 type(libpaw_libxc_type),intent(in),optional :: xc_functionals(2)

! *************************************************************************

 libpaw_libxc_isgga = .false.
 if (.not.libpaw_xc_constants_initialized) call libpaw_libxc_constants_load()

 if (present(xc_functionals)) then
   libpaw_libxc_isgga=(any(xc_functionals%family==LIBPAW_XC_FAMILY_GGA) .or. &
&                      any(xc_functionals%family==LIBPAW_XC_FAMILY_HYB_GGA))
 else
   libpaw_libxc_isgga=(any(paw_xc_global%family==LIBPAW_XC_FAMILY_GGA) .or. &
&                      any(paw_xc_global%family==LIBPAW_XC_FAMILY_HYB_GGA))
 end if

end function libpaw_libxc_isgga
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_ismgga
!! NAME
!!  libpaw_libxc_ismgga
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  is a Meta-GGA or not
!!
!! INPUTS
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional argument
!!                     XC functionals to initialize
!!
!! SOURCE

function libpaw_libxc_ismgga(xc_functionals)

!Arguments ------------------------------------
 logical :: libpaw_libxc_ismgga
 type(libpaw_libxc_type),intent(in),optional :: xc_functionals(2)

! *************************************************************************

 libpaw_libxc_ismgga = .false.
 if (.not.libpaw_xc_constants_initialized) call libpaw_libxc_constants_load()

 if (present(xc_functionals)) then
   libpaw_libxc_ismgga=(any(xc_functionals%family==LIBPAW_XC_FAMILY_MGGA) .or. &
&                       any(xc_functionals%family==LIBPAW_XC_FAMILY_HYB_MGGA))
 else
   libpaw_libxc_ismgga=(any(paw_xc_global%family==LIBPAW_XC_FAMILY_MGGA) .or. &
&                       any(paw_xc_global%family==LIBPAW_XC_FAMILY_HYB_MGGA))
 end if

end function libpaw_libxc_ismgga
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_needs_laplacian
!! NAME
!!  libpaw_libxc_needs_laplacian
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  needs the laplacian of the density or not
!!
!! INPUTS
!! [xc_functionals(2)]=<type(libxc_functional_type)>, optional argument
!!                     Handle for XC functionals
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function libpaw_libxc_needs_laplacian(xc_functionals)

!Arguments ------------------------------------
 implicit none
 logical :: libpaw_libxc_needs_laplacian
 type(libpaw_libxc_type),intent(in),optional :: xc_functionals(2)

! *************************************************************************

 libpaw_libxc_needs_laplacian = .false.

 if (present(xc_functionals)) then
   libpaw_libxc_needs_laplacian=(any(xc_functionals%needs_laplacian))
 else
   libpaw_libxc_needs_laplacian=(any(paw_xc_global%needs_laplacian))
 end if

 end function libpaw_libxc_needs_laplacian
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_is_hybrid
!! NAME
!!  libpaw_libxc_is_hybrid
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  is hybrid or not
!!
!! INPUTS
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional argument
!!                     XC functionals to initialize
!!
!! SOURCE

 function libpaw_libxc_is_hybrid(xc_functionals)

!Arguments ------------------------------------
 logical :: libpaw_libxc_is_hybrid
 type(libpaw_libxc_type),intent(in),optional :: xc_functionals(2)

! *************************************************************************

 libpaw_libxc_is_hybrid = .false.

 if (present(xc_functionals)) then
   libpaw_libxc_is_hybrid=(any(xc_functionals%is_hybrid))
 else
   libpaw_libxc_is_hybrid=(any(paw_xc_global%is_hybrid))
 end if

end function libpaw_libxc_is_hybrid
!!***

!----------------------------------------------------------------------

!!****f*  m_libpaw_libxc_funcs/libpaw_libxc_is_hybrid_from_id
!! NAME
!!  libpaw_libxc_is_hybrid_from_id
!!
!! FUNCTION
!!  Test function to identify whether a functional is hybrid or not, from its id
!!
!! INPUTS
!!  xcid= id of a LibXC functional
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function libpaw_libxc_is_hybrid_from_id(xcid)

!Arguments ------------------------------------
 logical :: libpaw_libxc_is_hybrid_from_id
 integer,intent(in) :: xcid
!Local variables-------------------------------
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 integer(C_INT) :: xcid_c
#endif

! *************************************************************************

#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 xcid_c=int(xcid,kind=C_INT)
 libpaw_libxc_is_hybrid_from_id =(libpaw_xc_func_is_hybrid_from_id(xcid_c)==1)
#else
 libpaw_libxc_is_hybrid_from_id = .false.
 if (.false.) write(std_out,*) xcid
#endif

end function libpaw_libxc_is_hybrid_from_id
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_has_kxc
!! NAME
!!  libpaw_libxc_has_kxc
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  provides Kxc or not (fxc in the libXC convention)
!!
!! INPUTS
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional argument
!!                     XC functionals
!!
!! SOURCE

function libpaw_libxc_has_kxc(xc_functionals)

!Arguments ------------------------------------
 logical :: libpaw_libxc_has_kxc
 type(libpaw_libxc_type),intent(in),optional,target :: xc_functionals(2)
!Local variables-------------------------------
 integer :: ii

! *************************************************************************

 libpaw_libxc_has_kxc=.true.

 do ii=1,2
   if (present(xc_functionals)) then
     if (.not.xc_functionals(ii)%has_fxc) libpaw_libxc_has_kxc=.false.
   else
     if (.not.paw_xc_global(ii)%has_fxc) libpaw_libxc_has_kxc=.false.
   end if
 end do

end function libpaw_libxc_has_kxc
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_has_k3xc
!! NAME
!!  libpaw_libxc_has_k3xc
!!
!! FUNCTION
!!  Test function to identify whether the presently used functional
!!  provides K3xc or not (kxc in the libXC convention)
!!
!! INPUTS
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional argument
!!                     Handle for XC functionals
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function libpaw_libxc_has_k3xc(xc_functionals)

!Arguments ------------------------------------
 logical :: libpaw_libxc_has_k3xc
 type(libpaw_libxc_type),intent(in),optional,target :: xc_functionals(2)
!Local variables-------------------------------
 integer :: ii

! *************************************************************************

 libpaw_libxc_has_k3xc=.true.

 do ii=1,2
   if (present(xc_functionals)) then
     if (.not.xc_functionals(ii)%has_kxc) libpaw_libxc_has_k3xc=.false.
   else
     if (.not.paw_xc_global(ii)%has_kxc) libpaw_libxc_has_k3xc=.false.
   end if
 end do

end function libpaw_libxc_has_k3xc
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_nspin
!! NAME
!!  libpaw_libxc_nspin
!!
!! FUNCTION
!!  Returns the number of spin components for the XC functionals
!!
!! INPUTS
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional argument
!!                     XC functionals to initialize
!!
!! SOURCE

function libpaw_libxc_nspin(xc_functionals)

!Arguments ------------------------------------
 integer :: libpaw_libxc_nspin
 type(libpaw_libxc_type),intent(in),optional :: xc_functionals(2)

! *************************************************************************

 libpaw_libxc_nspin = 1

 if (present(xc_functionals)) then
   if (any(xc_functionals%nspin==2)) libpaw_libxc_nspin=2
 else
   if (any(paw_xc_global%nspin==2)) libpaw_libxc_nspin=2
 end if

end function libpaw_libxc_nspin
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_getvxc
!! NAME
!!  libpaw_libxc_getvxc
!!
!! FUNCTION
!!  Return XC potential and energy, from input density (gradient etc...)
!!
!! INPUTS
!! ndvxc=size of dvxc
!! nd2vxc=size of d2vxc
!! npts=number of of points for the density
!! nspden=number of spin-density components
!! order=requested order of derivation
!! rho(npts,nspden)=electronic density
!! [grho2(npts,nspden)]=squared gradient of the density
!! [lrho(npts,nspden)]=laplacian of the density
!! [tau(npts,nspden)]= kinetic energy density
!!
!! OUTPUT
!! exc(npts)=XC energy density
!! vxc(npts,nspden)=derivative of the energy density wrt to the density
!! [vxclrho(npts,nspden)]=derivative of the energy density wrt to the density laplacian
!! [vxctau(npts,nspden)]=derivative of the energy density wrt to the kinetic energy density
!! [dvxc(npts,ndvxc)]=2nd derivative of the energy density wrt to the density
!! [vxcgr(npts,3)]=2nd derivative of the energy density wrt to the gradient
!!                 2nd derivative of the energy density wrt to the density and the gradient
!! [d2vxc(npts,nd2vxc)]=3rd derivative of the energy density wrt to the density
!!
!! SIDE EFFECTS
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional argument
!!                     XC functionals to initialize
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 subroutine libpaw_libxc_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho,exc,vxc,&
&           grho2,vxcgr,lrho,vxclrho,tau,vxctau,dvxc,d2vxc,xc_functionals) ! Optional arguments

!Arguments ------------------------------------
 integer, intent(in) :: ndvxc,nd2vxc,npts,nspden,order
 real(dp),intent(in)  :: rho(npts,nspden)
 real(dp),intent(out) :: vxc(npts,nspden),exc(npts)
 real(dp),intent(in),optional :: grho2(npts,2*min(nspden,2)-1)
 real(dp),intent(out),optional :: vxcgr(npts,3)
 real(dp),intent(in),optional :: lrho(npts,nspden)
 real(dp),intent(out),optional :: vxclrho(npts,nspden)
 real(dp),intent(in),optional :: tau(npts,nspden)
 real(dp),intent(out),optional :: vxctau(npts,nspden)
 real(dp),intent(out),optional :: dvxc(npts,ndvxc)
 real(dp),intent(out),optional :: d2vxc(npts,nd2vxc)
 type(libpaw_libxc_type),intent(inout),optional,target :: xc_functionals(2)
!Local variables -------------------------------
!scalars
 integer  :: ii,ipts
 logical :: is_gga,is_mgga,needs_laplacian
 character(len=500) :: msg
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 type(C_PTR) :: rho_c,sigma_c,lrho_c,tau_c
#endif
!arrays
 real(dp),target :: rhotmp(nspden),sigma(3),exctmp,vxctmp(nspden),vsigma(3)
 real(dp),target :: v2rho2(3),v2rhosigma(6),v2sigma2(6)
 real(dp),target :: v3rho3(4),v3rho2sigma(9),v3rhosigma2(12),v3sigma3(10)
 real(dp),target :: lrhotmp(nspden),tautmp(nspden),vlrho(nspden),vtau(nspden)
 type(libpaw_libxc_type),pointer :: xc_funcs(:)
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 type(C_PTR) :: exc_c(2),vxc_c(2),vsigma_c(2),vlrho_c(2),vtau_c(2)
 type(C_PTR) :: v2rho2_c(2),v2rhosigma_c(2),v2sigma2_c(2)
 type(C_PTR) :: v3rho3_c(2),v3rho2sigma_c(2),v3rhosigma2_c(2),v3sigma3_c(2)
#endif

! *************************************************************************

 if (.not.libpaw_xc_constants_initialized) call libpaw_libxc_constants_load()

!Select XC functional(s)
 if (present(xc_functionals)) then
   xc_funcs => xc_functionals
 else
   xc_funcs => paw_xc_global
 end if

 is_gga =libpaw_libxc_isgga (xc_funcs)
 is_mgga=libpaw_libxc_ismgga(xc_funcs)
 needs_laplacian=(libpaw_libxc_needs_laplacian(xc_funcs).and.present(lrho))

 if (is_gga.and.(.not.present(grho2))) then
   msg='GGA needs gradient of density!'
   LIBPAW_BUG(msg)
 end if
 if (is_mgga) then
   if (present(vxctau).and.(.not.present(tau))) then
     msg='meta-GGA needs tau!'
     LIBPAW_BUG(msg)
   end if
   if (needs_laplacian) then
     if (present(vxclrho).and.(.not.present(lrho))) then
       msg='meta-GGA needs lrho!'
       LIBPAW_BUG(msg)
     end if
   end if
 endif

!Inititalize all output arrays to zero
 exc=zero ; vxc=zero
 if (present(dvxc)) dvxc=zero
 if (present(d2vxc)) d2vxc=zero
 if (is_gga.or.is_mgga.and.present(vxcgr)) vxcgr=zero
 if (is_mgga.and.present(vxclrho)) vxclrho=zero
 if (is_mgga.and.present(vxctau)) vxctau=zero

!Determine which XC outputs can be computed
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 do ii = 1,2
   if (xc_funcs(ii)%has_exc) then
     exc_c(ii)=c_loc(exctmp)
   else
     exc_c(ii)=C_NULL_PTR
   end if
   if (xc_funcs(ii)%has_vxc) then
     vxc_c(ii)=c_loc(vxctmp)
     vsigma_c(ii)=c_loc(vsigma)
     vtau_c(ii)=c_loc(vtau)
     vlrho_c(ii)=c_loc(vlrho)
   else
     vxc_c(ii)=C_NULL_PTR
     vsigma_c(ii)=c_NULL_PTR
     vtau_c(ii)=C_NULL_PTR
     vlrho_c(ii)=C_NULL_PTR
   end if
   if ((xc_funcs(ii)%has_fxc).and.(abs(order)>1)) then
     v2rho2_c(ii)=c_loc(v2rho2)
     v2sigma2_c(ii)=c_loc(v2sigma2)
     v2rhosigma_c(ii)=c_loc(v2rhosigma)
   else
     v2rho2_c(ii)=C_NULL_PTR
     v2sigma2_c(ii)=C_NULL_PTR
     v2rhosigma_c(ii)=C_NULL_PTR
   end if
   if ((xc_funcs(ii)%has_kxc).and.(abs(order)>2)) then
     v3rho3_c(ii)=c_loc(v3rho3)
     v3sigma3_c(ii)=c_loc(v3sigma3)
     v3rho2sigma_c(ii)=c_loc(v3rho2sigma)
     v3rhosigma2_c(ii)=c_loc(v3rhosigma2)
   else
     v3rho3_c(ii)=C_NULL_PTR
     v3sigma3_c(ii)=C_NULL_PTR
     v3rho2sigma_c(ii)=C_NULL_PTR
     v3rhosigma2_c(ii)=C_NULL_PTR
   end if
 end do
#endif

!Initialize temporary arrays
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 rhotmp=zero ; rho_c=c_loc(rhotmp)
 if (is_gga.or.is_mgga) then
   sigma=zero ; sigma_c=c_loc(sigma)
 end if
 if (is_mgga) then
   tautmp=zero ; tau_c=c_loc(tautmp)
   lrhotmp=zero ; lrho_c=c_loc(lrhotmp)
 end if
#endif

!Some mGGA functionals require a special treatment
 if (is_mgga) then
   !TB09 functional requires the c parameter to be set
   call libpaw_libxc_compute_tb09(npts,nspden,rho,grho2,xc_funcs)
 end if

!Loop over points
 do ipts=1,npts

!  Convert the quantities provided to the ones needed by libxc
   if (nspden == 1) then
     ! rho_up is passed in the spin-unpolarized case, while the libxc
     ! expects the total density
     rhotmp(1:nspden) = two*rho(ipts,1:nspden)
   else
     rhotmp(1:nspden) = rho(ipts,1:nspden)
   end if
   if (is_gga.or.is_mgga) then
     if (nspden==1) then
       ! |grho_up|^2 is passed while Libxc needs |grho_tot|^2
       sigma(1) = four*grho2(ipts,1)
     else
       ! |grho_up|^2, |grho_dn|^2, and |grho_tot|^2 are passed
       ! while Libxc needs |grho_up|^2, grho_up.grho_dn, and |grho_dn|^2
       sigma(1) = grho2(ipts,1)
       sigma(2) = (grho2(ipts,3) - grho2(ipts,1) - grho2(ipts,2))/two
       sigma(3) = grho2(ipts,2)
     end if
   end if
   if (is_mgga) then
     if (nspden==1) then
       tautmp(1:nspden) = two*tau(ipts,1:nspden)
       if (needs_laplacian) lrhotmp(1:nspden) = two*lrho(ipts,1:nspden)
     else
       tautmp(1:nspden) = tau(ipts,1:nspden)
       if (needs_laplacian) lrhotmp(1:nspden) = lrho(ipts,1:nspden)
     end if
   end if

!  Loop over functionals
   do ii = 1,2
     if (xc_funcs(ii)%id<=0) cycle

!    Get the energy and the potential (and possibly the other derivatives)
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
     exctmp=zero ; vxctmp=zero
!    ===== LDA =====
     if (xc_funcs(ii)%family==LIBPAW_XC_FAMILY_LDA) then
       exctmp=zero ; vxctmp=zero ; v2rho2=zero ; v3rho3=zero
       call libpaw_xc_get_lda(xc_funcs(ii)%conf,1,rho_c, &
&                  exc_c(ii),vxc_c(ii),v2rho2_c(ii),v3rho3_c(ii))
!    ===== GGA =====
     else if (xc_funcs(ii)%family==LIBPAW_XC_FAMILY_GGA.or. &
&             xc_funcs(ii)%family==LIBPAW_XC_FAMILY_HYB_GGA) then
       exctmp=zero ; vxctmp=zero ; vsigma=zero
       v2rho2=zero ; v2sigma2=zero ; v2rhosigma=zero
       v3rho3=zero ; v3rho2sigma=zero ; v3rhosigma2=zero ; v3sigma3=zero
       call libpaw_xc_get_gga(xc_funcs(ii)%conf,1,rho_c,sigma_c, &
&                  exc_c(ii),vxc_c(ii),vsigma_c(ii), &
&                  v2rho2_c(ii),v2rhosigma_c(ii),v2sigma2_c(ii), &
&                  v3rho3_c(ii),v3rho2sigma_c(ii),v3rhosigma2_c(ii),v3sigma3_c(ii))
!    ===== mGGA =====
     else if (xc_funcs(ii)%family==LIBPAW_XC_FAMILY_MGGA.or. &
&             xc_funcs(ii)%family==LIBPAW_XC_FAMILY_HYB_GGA) then
       exctmp=zero ; vxctmp=zero ; vsigma=zero ; vlrho=zero ; vtau=zero
       call libpaw_xc_get_mgga(xc_funcs(ii)%conf,1,rho_c,sigma_c,lrho_c,tau_c, &
&                  exc_c(ii),vxc_c(ii),vsigma_c(ii),vlrho_c(ii),vtau_c(ii), &
&                  C_NULL_PTR,C_NULL_PTR,C_NULL_PTR,C_NULL_PTR,C_NULL_PTR, &
&                  C_NULL_PTR,C_NULL_PTR,C_NULL_PTR,C_NULL_PTR,C_NULL_PTR)
     end if
#endif

     exc(ipts) = exc(ipts) + exctmp
     vxc(ipts,1:nspden) = vxc(ipts,1:nspden) + vxctmp(1:nspden)

!    Deal with fxc and kxc
     if (abs(order)>1) then
!      ----- LDA -----
       if (xc_funcs(ii)%family==LIBPAW_XC_FAMILY_LDA) then
         if (nspden==1) then
           if(order>=2) then
             dvxc(ipts,1)=dvxc(ipts,1)+v2rho2(1)
             if(order>2) then
               d2vxc(ipts,1)=d2vxc(ipts,1)+v3rho3(1)
             endif
           else if (order==-2) then
             dvxc(ipts,1)=dvxc(ipts,1)+v2rho2(1)
             dvxc(ipts,2)=dvxc(ipts,2)+v2rho2(1)
           endif
         else
           dvxc(ipts,1)=dvxc(ipts,1)+v2rho2(1)
           dvxc(ipts,2)=dvxc(ipts,2)+v2rho2(2)
           dvxc(ipts,3)=dvxc(ipts,3)+v2rho2(3)
           if(abs(order)>2) then
             d2vxc(ipts,1)=d2vxc(ipts,1)+v3rho3(1)
             d2vxc(ipts,2)=d2vxc(ipts,2)+v3rho3(2)
             d2vxc(ipts,3)=d2vxc(ipts,3)+v3rho3(3)
             d2vxc(ipts,4)=d2vxc(ipts,4)+v3rho3(4)
           endif
         endif
!      ----- GGA -----
       else if (xc_funcs(ii)%family==LIBPAW_XC_FAMILY_GGA.or. &
&               xc_funcs(ii)%family==LIBPAW_XC_FAMILY_HYB_GGA) then
         if (xc_funcs(ii)%kind==LIBPAW_XC_EXCHANGE) then
           if (nspden==1) then
             dvxc(ipts,1)=v2rho2(1)*two
             dvxc(ipts,2)=dvxc(ipts,1)
             dvxc(ipts,3)=two*two*vsigma(1)
             dvxc(ipts,4)=dvxc(ipts,3)
             dvxc(ipts,5)=four*two*v2rhosigma(1)
             dvxc(ipts,6)=dvxc(ipts,5)
             dvxc(ipts,7)=two*four*four*v2sigma2(1)
             dvxc(ipts,8)=dvxc(ipts,7)
           else
             dvxc(ipts,1)=v2rho2(1)
             dvxc(ipts,2)=v2rho2(3)
             dvxc(ipts,3)=two*vsigma(1)
             dvxc(ipts,4)=two*vsigma(3)
             dvxc(ipts,5)=two*v2rhosigma(1)
             dvxc(ipts,6)=two*v2rhosigma(6)
             dvxc(ipts,7)=four*v2sigma2(1)
             dvxc(ipts,8)=four*v2sigma2(6)
           end if
         else if (xc_funcs(ii)%kind==LIBPAW_XC_CORRELATION) then
           if (nspden==1) then
             dvxc(ipts,9)=v2rho2(1)
             dvxc(ipts,10)=dvxc(ipts,9)
             dvxc(ipts,11)=dvxc(ipts,9)
             dvxc(ipts,12)=two*vsigma(1)
             dvxc(ipts,13)=two*v2rhosigma(1)
             dvxc(ipts,14)=dvxc(ipts,13)
             dvxc(ipts,15)=four*v2sigma2(1)
           else
             dvxc(ipts,9)=v2rho2(1)
             dvxc(ipts,10)=v2rho2(2)
             dvxc(ipts,11)=v2rho2(3)
             dvxc(ipts,12)=two*vsigma(1)
             dvxc(ipts,13)=two*v2rhosigma(1)
             dvxc(ipts,14)=two*v2rhosigma(6)
             dvxc(ipts,15)=four*v2sigma2(1)
           end if
         end if
       end if
     end if

!    Convert the quantities returned by Libxc
     if ((is_gga.or.is_mgga).and.present(vxcgr)) then
       if (nspden==1) then
         vxcgr(ipts,3) = vxcgr(ipts,3) + vsigma(1)*two
       else
         vxcgr(ipts,1) = vxcgr(ipts,1) + two*vsigma(1) - vsigma(2)
         vxcgr(ipts,2) = vxcgr(ipts,2) + two*vsigma(3) - vsigma(2)
         vxcgr(ipts,3) = vxcgr(ipts,3) + vsigma(2)
       end if
     end if
     if (is_mgga.and.present(vxctau)) then
       vxctau(ipts,1:nspden)  = vxctau(ipts,1:nspden)  + vtau(1:nspden)
     end if
     if (is_mgga.and.needs_laplacian.and.present(vxclrho)) then
       vxclrho(ipts,1:nspden) = vxclrho(ipts,1:nspden) + vlrho(1:nspden)
     end if

   end do ! ii
 end do   ! ipts

end subroutine libpaw_libxc_getvxc
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_get_hybridparams
!! NAME
!!  libpaw_libxc_get_hybridparams
!!
!! FUNCTION
!!  Returns the parameters of an hybrid functional (mixing coefficient(s) and range separation)
!!
!! INPUTS
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional argument
!!                     XC functionals to initialize
!!
!! OUTPUT
!!  [hyb_mixing]  = mixing factor of Fock contribution
!!  [hyb_mixing_sr]= mixing factor of short-range Fock contribution
!!  [hyb_range]    = Range (for separation)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine libpaw_libxc_get_hybridparams(hyb_mixing,hyb_mixing_sr,hyb_range,xc_functionals)

!Arguments ------------------------------------
 real(dp),intent(out),optional :: hyb_mixing,hyb_mixing_sr,hyb_range
 type(libpaw_libxc_type),intent(in),optional,target :: xc_functionals(2)
!Local variables -------------------------------
 integer :: ii
 character(len=500) :: msg
 type(libpaw_libxc_type),pointer :: xc_func

! *************************************************************************

 if (present(hyb_mixing   )) hyb_mixing   =zero
 if (present(hyb_mixing_sr)) hyb_mixing_sr=zero
 if (present(hyb_range    )) hyb_range    =zero

 do ii = 1, 2

!  Select XC functional
   if (present(xc_functionals)) then
     xc_func => xc_functionals(ii)
   else
     xc_func => paw_xc_global(ii)
   end if

!  Mixing coefficient for the Fock contribution
   if (present(hyb_mixing)) then
     if (abs(xc_func%hyb_mixing) > tol8) then
       if (abs(hyb_mixing) <= tol8) then
         hyb_mixing=xc_func%hyb_mixing
       else
         msg='Invalid XC functional: contains 2 hybrid exchange functionals!'
         LIBPAW_ERROR(msg)
       end if
     end if
   end if

!  Mixing coefficient for the short-range Fock contribution
   if (present(hyb_mixing_sr)) then
     if (abs(xc_func%hyb_mixing_sr) > tol8) then
       if (abs(hyb_mixing_sr) <= tol8) then
         hyb_mixing_sr=xc_func%hyb_mixing_sr
       else
         msg='Invalid XC functional: contains 2 hybrid exchange functionals!'
         LIBPAW_ERROR(msg)
       end if
     end if
   end if

!  Range separation
   if (present(hyb_range)) then
     if (abs(xc_func%hyb_range) > tol8) then
       if (abs(hyb_range) <= tol8) then
         hyb_range=xc_func%hyb_range
       else
         msg='Invalid XC functional: contains 2 hybrid exchange functionals!'
         LIBPAW_ERROR(msg)
       end if
     end if
   end if

 end do

end subroutine libpaw_libxc_get_hybridparams
!!***

!----------------------------------------------------------------------

!!****f* libpaw_libxc_funcs/libpaw_libxc_set_hybridparams
!! NAME
!!  libpaw_libxc_set_hybridparams
!!
!! FUNCTION
!!  Set the parameters of an hybrid functional (mixing coefficient(s) and range separation)
!!
!! INPUTS
!! [hyb_mixing]       = mixing factor of Fock contribution
!! [hyb_mixing_sr]    = mixing factor of short-range Fock contribution
!! [hyb_range]        = Range (for separation)
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional argument
!!                     XC functionals to initialize
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine libpaw_libxc_set_hybridparams(hyb_mixing,hyb_mixing_sr,hyb_range,xc_functionals)

!Arguments ------------------------------------
 real(dp),intent(in),optional :: hyb_mixing,hyb_mixing_sr,hyb_range
 type(libpaw_libxc_type),intent(in),optional,target :: xc_functionals(2)
!Local variables -------------------------------
 integer :: ii,id_pbe0,id_hse03,id_hse06
 logical :: is_pbe0,is_hse
 character(len=500) :: msg
 type(libpaw_libxc_type),pointer :: xc_func
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 integer(C_INT) :: npar_c
 real(C_DOUBLE) :: alpha_c,beta_c,omega_c,param_c(3)
#endif

! *************************************************************************

 is_pbe0=.false.
 is_hse =.false.
 id_pbe0=libpaw_libxc_getid('HYB_GGA_XC_PBEH')
 id_hse03=libpaw_libxc_getid('HYB_GGA_XC_HSE03')
 id_hse06=libpaw_libxc_getid('HYB_GGA_XC_HSE06')

 do ii = 1, 2

!  Select XC functional
   if (present(xc_functionals)) then
     xc_func => xc_functionals(ii)
   else
     xc_func => paw_xc_global(ii)
   end if

!  Doesnt work with all hybrid functionals
   if (is_pbe0.or.is_hse) then
     msg='Invalid XC functional: contains 2 hybrid exchange functionals!'
     LIBPAW_ERROR(msg)
   end if
   is_pbe0=(xc_func%id==id_pbe0)
   is_hse=((xc_func%id==id_hse03).or.(xc_func%id==id_hse06))
   if ((.not.is_pbe0).and.(.not.is_hse)) cycle

#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
!  First retrieve current values of parameters
   call xc_hyb_cam_coef(xc_func%conf,omega_c,alpha_c,beta_c)

!  New values for parameters
   if (present(hyb_mixing)) alpha_c=real(hyb_mixing,kind=C_DOUBLE)
   if (present(hyb_mixing_sr)) beta_c=real(hyb_mixing_sr,kind=C_DOUBLE)
   if (present(hyb_range)) omega_c=real(hyb_range,kind=C_DOUBLE)

!  PBE0: set parameters
   if (is_pbe0) then
     npar_c=int(1,kind=C_INT) ; param_c(1)=alpha_c
     call libpaw_xc_func_set_params(xc_func%conf,param_c,npar_c)
   end if

!  HSE: set parameters
   if (is_hse) then
     npar_c=int(3,kind=C_INT)
     param_c(1)=beta_c;param_c(2:3)=omega_c
     call libpaw_xc_func_set_params(xc_func%conf,param_c,npar_c)
   end if

#else
!  This is to avoid unused arguments
   if(.false. .and. present(hyb_mixing) .and. present(hyb_mixing_sr) .and. present(hyb_range))then
     msg='One should not be here'
   endif
#endif

 end do

 if ((.not.is_pbe0).and.(.not.is_hse)) then
   msg='Invalid XC functional: not able to change parameters for this functional!'
   LIBPAW_WARNING(msg)
 end if

end subroutine libpaw_libxc_set_hybridparams
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_gga_from_hybrid
!! NAME
!!  libpaw_libxc_gga_from_hybrid
!!
!! FUNCTION
!!  Returns a logical flag: TRUE if one can deduce, from the id of a hybrid functional,
!!  the id(s) of the GGA functional on which it is based.
!!  Optionally returns the id of the GGA functional on which the hybrid functional is based
!!  (2 integers defining the GGA X and C functionals).
!!  - If an id is provided as input argument, it is used as input id;
!!  - If not, the input id is taken from the optional xc_functionals datastructure;
!!  - If no input argument is given, the input id is taken from the global paw_xc_global datastructure.
!!
!! INPUTS
!! [hybrid_id]=<type(libpaw_libxc_type)>, optional : id of an input hybrid functional
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional : XC functionals from which
!!                     the id(s) can to be used
!!
!! OUTPUT
!! [gga_id(2)]=array that contains the GGA libXC id(s)
!! libpaw_libxc_gga_from_hybrid=.true. if the GGA has been found from the input id
!!
!! SOURCE

function libpaw_libxc_gga_from_hybrid(gga_id,hybrid_id,xc_functionals)

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: hybrid_id
 logical :: libpaw_libxc_gga_from_hybrid
!arrays
 integer,intent(out),optional :: gga_id(2)
 type(libpaw_libxc_type),intent(inout),optional,target :: xc_functionals(2)
!Local variables -------------------------------
!scalars
 integer :: ii
 logical :: is_hybrid
 character(len=100) :: c_name,x_name,msg
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 character(len=100) :: xc_name
 character(kind=C_CHAR,len=1),pointer :: strg_c
#endif
!arrays
 integer :: trial_id(2)

! *************************************************************************

 libpaw_libxc_gga_from_hybrid=.false.

 is_hybrid=.false.
 if (present(hybrid_id)) then
   trial_id(1)=hybrid_id
   trial_id(2)=0
   is_hybrid=libpaw_libxc_is_hybrid_from_id(trial_id(1))
 else if (present(xc_functionals)) then
   trial_id(1)=xc_functionals(1)%id
   trial_id(2)=xc_functionals(2)%id
   is_hybrid=libpaw_libxc_is_hybrid(xc_functionals)
 else
   trial_id(1)=paw_xc_global(1)%id
   trial_id(2)=paw_xc_global(2)%id
   is_hybrid=libpaw_libxc_is_hybrid(paw_xc_global)
 end if

 c_name="unknown" ; x_name="unknown"

 do ii = 1, 2

   if ((trial_id(ii)<=0).or.(.not.is_hybrid)) cycle

   if (libpaw_libxc_gga_from_hybrid) then
     msg='Invalid XC functional: contains 2 hybrid functionals!'
     LIBPAW_ERROR(msg)
   end if

#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING

   call c_f_pointer(xc_functional_get_name(trial_id(ii)),strg_c)
   call char_c_to_f(strg_c,xc_name)

!  AVAILABLE FUNCTIONALS

!  ===== PBE0 =====
   if (xc_name=="hyb_gga_xc_pbeh" .or. &
&      xc_name=="hyb_gga_xc_pbe0_13") then
     c_name="GGA_C_PBE"
     x_name="GGA_X_PBE"
     libpaw_libxc_gga_from_hybrid=.true.

!  ===== HSE =====
   else if (xc_name=="hyb_gga_xc_hse03" .or. &
&           xc_name=="hyb_gga_xc_hse06" ) then
     c_name="GGA_C_PBE"
     x_name="GGA_X_PBE"
     libpaw_libxc_gga_from_hybrid=.true.
   end if

#endif

 enddo ! ii

 if (present(gga_id)) then
   if (libpaw_libxc_gga_from_hybrid) then
     gga_id(1)=libpaw_libxc_getid(c_name)
     gga_id(2)=libpaw_libxc_getid(x_name)
   else
     gga_id(:)=-1
   end if
 end if

 if (.not.libpaw_libxc_gga_from_hybrid) then
   msg='Unable to find a GGA functional for this hybrid!'
   LIBPAW_ERROR(msg)
 end if

end function libpaw_libxc_gga_from_hybrid
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/libpaw_libxc_compute_tb09
!! NAME
!!  libpaw_libxc_compute_tb09
!!
!! FUNCTION
!!  Compute c parameter for Tran-Blaha 2009 functional and set it
!!
!! INPUTS
!! npts=number of of points for the density
!! nspden=number of spin-density components
!! rho(npts,nspden)=electronic density
!! grho2(npts,nspden)=squared gradient of the density
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! [xc_functionals(2)]=<type(libpaw_libxc_type)>, optional argument
!!                     XC functionals to initialize
!!
!! PARENTS
!!      m_libpaw_libxc
!!
!! CHILDREN
!!
!! SOURCE

 subroutine libpaw_libxc_compute_tb09(npts,nspden,rho,grho2,xc_functionals)

!Arguments ------------------------------------
 integer, intent(in) :: npts,nspden
 real(dp),intent(in)  :: rho(npts,nspden),grho2(npts,2*min(nspden,2)-1)
 type(libpaw_libxc_type),intent(inout),optional,target :: xc_functionals(2)
!Local variables -------------------------------
!scalars
 integer  :: ii,ipts
 logical :: fixed_c_tb09,is_mgga_tb09
 real(dp) :: cc
 character(len=500) :: msg
!arrays
 type(libpaw_libxc_type),pointer :: xc_funcs(:)
 real(dp),allocatable :: gnon(:)
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 integer(C_INT) :: npar_c=int(2,kind=C_INT)
 real(C_DOUBLE) :: param_c(2)
#endif

! *************************************************************************

 if (.not.libpaw_xc_constants_initialized) call libpaw_libxc_constants_load()

!Select XC functional(s)
 if (present(xc_functionals)) then
   xc_funcs => xc_functionals
 else
   xc_funcs => paw_xc_global
 end if

 is_mgga_tb09=(any(xc_funcs%id==libpaw_libxc_getid('XC_MGGA_X_TB09')))
 fixed_c_tb09=(any(abs(xc_funcs%xc_tb09_c-99.99_dp)>tol12))

 if (is_mgga_tb09) then

!  C is fixed by the user
   if (fixed_c_tb09) then
     cc=zero
     do ii=1,2
       if (abs(xc_funcs(ii)%xc_tb09_c-99.99_dp)>tol12) cc=xc_funcs(ii)%xc_tb09_c
     end do
     write(msg,'(2a,f9.6)' ) ch10,&
&    'In the mGGA functional TB09, c is fixed by the user and is equal to ',cc
     call wrtout(std_out,msg,'COLL')
!  C is computed
   else
     LIBPAW_ALLOCATE(gnon,(npts))
     do ipts=1,npts
       if (sum(rho(ipts,:))<=1e-7_dp) then
         gnon(ipts)=zero
       else
         if (nspden==1) then
           gnon(ipts)=sqrt(grho2(ipts,1))/rho(ipts,1)
         else
           gnon(ipts)=sqrt(grho2(ipts,3))/sum(rho(ipts,:))
         end if
       end if
     end do
     cc= -0.012_dp + 1.023_dp*sqrt(sum(gnon)/npts)
     LIBPAW_DEALLOCATE(gnon)
     write(msg,'(2a,f9.6)' ) ch10,'In the mGGA functional TB09, c = ',cc
     call wrtout(std_out,msg,'COLL')
   end if

!  Set c in XC data structure
   do ii=1,2
     if (xc_funcs(ii)%id==libpaw_libxc_getid('XC_MGGA_X_TB09')) then
#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
     param_c(1)=real(cc,kind=C_DOUBLE) ; param_c(2)=real(0._dp,kind=C_DOUBLE)
     call libpaw_xc_func_set_params(xc_funcs(ii)%conf,param_c,npar_c)
#endif
     end if
   end do
 end if

end subroutine libpaw_libxc_compute_tb09
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/char_f_to_c
!! NAME
!!  char_f_to_c
!!
!! FUNCTION
!! Helper function to convert a Fortran string to a C string
!! Based on a routine by Joseph M. Krahn
!!
!! INPUTS
!!  f_string=Fortran string
!!
!! OUTPUT
!!  c_string=C string
!!
!! SOURCE

#if defined LIBPAW_ISO_C_BINDING
function char_f_to_c(f_string) result(c_string)
!Arguments ------------------------------------
 character(len=*),intent(in) :: f_string
 character(kind=C_CHAR,len=1) :: c_string(len_trim(f_string)+1)
!Local variables -------------------------------
 integer :: ii,strlen
!! *************************************************************************
 strlen=len_trim(f_string)
 forall(ii=1:strlen)
   c_string(ii)=f_string(ii:ii)
 end forall
 c_string(strlen+1)=C_NULL_CHAR
 end function char_f_to_c
#endif
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_libxc_funcs/char_c_to_f
!! NAME
!!  char_c_to_f
!!
!! FUNCTION
!! Helper function to convert a C string to a Fortran string
!! Based on a routine by Joseph M. Krahn
!!
!! INPUTS
!!  c_string=C string
!!
!! OUTPUT
!!  f_string=Fortran string
!!
!! PARENTS
!!      m_libpaw_libxc
!!
!! CHILDREN
!!
!! SOURCE

#if defined LIBPAW_ISO_C_BINDING
subroutine char_c_to_f(c_string,f_string)
!Arguments ------------------------------------
 character(kind=C_CHAR,len=1),intent(in) :: c_string(*)
 character(len=*),intent(out) :: f_string
!Local variables -------------------------------
 integer :: ii
!! *************************************************************************
 ii=1
 do while(c_string(ii)/=C_NULL_CHAR.and.ii<=len(f_string))
   f_string(ii:ii)=c_string(ii) ; ii=ii+1
 end do
 if (ii<len(f_string)) f_string(ii:)=' '
 end subroutine char_c_to_f
#endif
!!***

!----------------------------------------------------------------------

end module m_libpaw_libxc_funcs
!!***

!======================================================================
!======================================================================

!!****m* ABINIT/m_libpaw_libxc
!! NAME
!!  m_libpaw_libxc
!!
!! FUNCTION
!!  Module used to interface libPAW with host code.
!!  At present, two cases are implemented:
!!   - Use of ABINIT m_libxc_functional module
!!   - Use of embedded m_libpaw_libxc_funcs module
!!
!! COPYRIGHT
!! Copyright (C) 2014-2020 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

module m_libpaw_libxc

#if defined HAVE_LIBPAW_ABINIT
 use libxc_functionals

#else
 use m_libpaw_libxc_funcs, only : &
& libxc_functionals_check            => libpaw_libxc_check, &
& libxc_functionals_init             => libpaw_libxc_init, &
& libxc_functionals_end              => libpaw_libxc_end, &
& libxc_functionals_fullname         => libpaw_libxc_fullname, &
& libxc_functionals_getid            => libpaw_libxc_getid, &
& libxc_functionals_family_from_id   => libpaw_libxc_family_from_id, &
& libxc_functionals_ixc              => libpaw_libxc_ixc, &
& libxc_functionals_getvxc           => libpaw_libxc_getvxc, &
& libxc_functionals_isgga            => libpaw_libxc_isgga, &
& libxc_functionals_ismgga           => libpaw_libxc_ismgga, &
& libxc_functionals_needs_laplacian  => libpaw_libxc_needs_laplacian, &
& libxc_functionals_is_hybrid        => libpaw_libxc_is_hybrid, &
& libxc_functionals_has_kxc          => libpaw_libxc_has_kxc, &
& libxc_functionals_nspin            => libpaw_libxc_nspin, &
& libxc_functionals_get_hybridparams => libpaw_libxc_get_hybridparams, &
& libxc_functionals_set_hybridparams => libpaw_libxc_set_hybridparams, &
& libxc_functionals_gga_from_hybrid  => libpaw_libxc_gga_from_hybrid
#endif

 implicit none

end module m_libpaw_libxc
!!***
