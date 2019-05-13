"""
Units Tests for Parser
"""
from __future__ import print_function, division, unicode_literals, absolute_import

from unittest import TestCase

from .parser import FortranKissParser
from .project import AbinitProject

class TestTrickyCode(TestCase):

    def test_tricky_module(self):
        """Parsing tricky Fortran module similar to 39_libpaw/m_libpaw_libxc.F90"""
        s = """\
#include "libpaw.h"

module m_libpaw_libxc_funcs

 USE_DEFS
 USE_MSG_HANDLING
 USE_MEMORY_PROFILING
#ifdef LIBPAW_ISO_C_BINDING
 use iso_c_binding
#endif

 implicit none
 private

!Public functions
 public :: libpaw_libxc_check              ! Check if the code has been compiled with libXC

!Private functions
 private :: libpaw_libxc_constants_load    ! Load libXC constants from C headers
 private :: char_c_to_f                    ! Convert a string from C to Fortran

!Public constants (use libpaw_libxc_constants_load to init them)
 integer,public,save :: LIBPAW_XC_FAMILY_UNKNOWN       = -1
 logical,private,save :: libpaw_xc_constants_initialized=.false.

!XC functional public type
 type,public :: libpaw_libxc_type
   integer  :: id              ! identifier
   logical  :: has_exc         ! TRUE is exc is available for the functional
   real(dp) :: hyb_mixing      ! Hybrid functional: mixing factor of Fock contribution (default=0)
#ifdef LIBPAW_ISO_C_BINDING
   type(C_PTR),pointer :: conf => null() ! C pointer to the functional itself
#endif
 end type libpaw_libxc_type

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
   subroutine xc_mgga(xc_func,np,rho,sigma,lapl,tau,zk,vrho,vsigma,vlapl,vtau, &
&             v2rho2,v2sigma2,v2lapl2,v2tau2,v2rhosigma,v2rholapl,v2rhotau, &
&             v2sigmalapl,v2sigmatau,v2lapltau) &
&             bind(C,name="xc_mgga")
     use iso_c_binding, only : C_INT,C_PTR
     integer(C_INT),value :: np
     type(C_PTR),value :: rho,sigma,lapl,tau,zk,vrho,vsigma,vlapl,vtau, &
&                         v2rho2,v2sigma2,v2lapl2,v2tau2,v2rhosigma,v2rholapl,v2rhotau, &
&                         v2sigmalapl,v2sigmatau,v2lapltau
     type(C_PTR) :: xc_func
   end subroutine xc_mgga
 end interface
#endif

contains

 subroutine libpaw_libxc_constants_load()

#if defined LIBPAW_HAVE_LIBXC && defined LIBPAW_ISO_C_BINDING
 integer(C_INT) :: i1,i2,i3,i4,i5,i6,i7,i8
#endif

 end subroutine libpaw_libxc_constants_load

 subroutine libpaw_libxc_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho,exc,vxc,&
&           grho2,vxcgr,lrho,vxclrho,tau,vxctau,dvxc,d2vxc,xc_tb09_c,xc_functionals) ! Optional arguments

 implicit none

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
 real(dp),intent(in),optional :: xc_tb09_c
 type(libpaw_libxc_type),intent(inout),optional,target :: xc_functionals(2)

! *************************************************************************

end subroutine libpaw_libxc_getvxc

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

#if defined LIBPAW_ISO_C_BINDING
subroutine char_c_to_f(c_string,f_string)

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

module m_libpaw_libxc

#if defined HAVE_LIBPAW_ABINIT
 use libxc_functionals
#else
 use m_libpaw_libxc_funcs, only : &
& libxc_functionals_check            => libpaw_libxc_check, &
& libxc_functionals_gga_from_hybrid  => libpaw_libxc_gga_from_hybrid
#endif

 implicit none

end module m_libpaw_libxc
"""
        p = FortranKissParser(macros=AbinitProject.MACROS)
        p.parse_string(s)
        assert len(p.modules) == 2
        mod1, mod2 = p.modules
        assert mod1.name == "m_libpaw_libxc_funcs" and mod2.name == "m_libpaw_libxc"
        assert "defs_basis" in mod1.uses
        assert len(mod1.types) == 1
        libxc_t = mod1.types[0]
        assert libxc_t.name == "libpaw_libxc_type"
        assert len(mod1.interfaces) == 5
