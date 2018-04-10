"""
Units Tests for Parser
"""
from __future__ import print_function, division, unicode_literals, absolute_import

from unittest import TestCase

from .parser import FortranKissParser


class TestFortranKissParser(TestCase):

    def test_regex(self):
        """Test regular expressions."""
        p = FortranKissParser

        # Note: Don't use strings with inlined comments e.g.
        #
        #   "program foo !this is an inlined comment"
        #
        # because the parser will strip them and operate stripped strings.

        # Program
        m = p.RE_PROG_START.match("program foo")
        assert m and m.group("name") == "foo"
        assert p.RE_PROG_START.match("program")
        m = p.RE_PROG_END.match("end program foo")
        assert m and m.group("name") == "foo"
        m = p.RE_PROG_END.match("end")
        assert m and m.group("name") is None

        # Module
        m = p.RE_MOD_START.match("module   foo")
        assert m and m.group("name") == "foo"
        #m = p.RE_MOD_START.match("module   foo !hello")
        #assert m and m.group("name") == "foo"
        m = p.RE_MOD_END.match("end   module   foo")
        assert m and m.group("name") == "foo"
        m = p.RE_MOD_END.match("end")
        assert m and m.group("name") is None

        # Subroutine
        m = p.RE_SUB_START.match("subroutine  foo")
        assert m and m.group("name") == "foo" and not m.group("prefix")
        m = p.RE_SUB_START.match("subroutine foo(a, b)")
        assert m and m.group("name") == "foo" and m.group("prefix") == ""
        m = p.RE_SUB_START.match("pure  subroutine foo(a, b)")
        assert m and m.group("name") == "foo" and m.group("prefix").strip() == "pure"
        m = p.RE_SUB_END.match("end  subroutine foo")
        assert m and m.group("name") == "foo"
        m = p.RE_SUB_END.match("end")
        assert m and m.group("name") is None

        # Functions: these are tricky so test them carefully.
        m = p.RE_FUNC_START.match("function  foo")
        assert m and m.group("name") == "foo" and not m.group("prefix")
        m = p.RE_FUNC_START.match("integer function foo(a, b)")
        assert m and m.group("name") == "foo" #and not m.group("prefix")
        m = p.RE_FUNC_START.match("integer(kind = i8b) function foo(a, b)")
        assert m and m.group("name") == "foo" #and not m.group("prefix")
        m = p.RE_FUNC_START.match("real(dp) function foo(a, b)")
        assert m and m.group("name") == "foo" #and not m.group("prefix")
        m = p.RE_FUNC_START.match("pure complex(dpc) function foo(a)")
        assert m and m.group("name") == "foo" #and not m.group("prefix")
        m = p.RE_FUNC_START.match("pure logical function foo(a, b)")
        assert m and m.group("name") == "foo"
        #assert m.group("prefix") == "foo"
        m = p.RE_FUNC_START.match("type(foo_t) function foo(a, b) result(res)")
        assert m and m.group("name") == "foo"

        m = p.RE_FUNC_END.match("end  function foo")
        assert m and m.group("name") == "foo"
        m = p.RE_FUNC_END.match("end")
        assert m and m.group("name") is None

        # Calls to subroutines
        m = p.RE_SUBCALL.match("call  foo")
        assert m and m.group("name") == "foo"
        m = p.RE_SUBCALL.match("call  foo ( a, b)")
        assert m and m.group("name") == "foo"
        m = p.RE_SUBCALL.match("if (a == 1) call foo(a=2)")
        assert m and m.group("name") == "foo"
        #m = p.RE_SUBCALL.match("if (cond) call obj%foo(a, b)")
        #assert m and m.group("name") == "obj%foo"
        #m = p.RE_SUBCALL.match("if (.True.) call obj%foo(a=1, b=3)")
        #assert m and m.group("name") == "obj%foo"

        # Interface
        m = p.RE_INTERFACE_START.match("abstract interface foo")
        assert m #and m.group("name") == "foo"
        m = p.RE_INTERFACE_END.match("end interface")
        assert m and not m.group("name")
        m = p.RE_INTERFACE_END.match("end  interface  foo")
        assert m and m.group("name").strip() == "foo"

    def test_simple_program(self):
        """Parsing Fortran program written following good programming standards."""
        s = """
#include "abi_common.h"

! A Very simple example
! of Fortran program written following good programming standards

program hello_world
    use defs_basis
    implicit none
    call foo()

contains ! contained procedure

    subroutine foo()
    end subroutine foo

    function bar()
    end function bar
end program
"""
        p = FortranKissParser()
        p.parse_lines(s.splitlines())
        assert len(p.programs) == 1
        assert not p.modules and not p.subroutines
        assert "abi_common.h" in p.includes
        assert "defs_basis" in p.uses

        prog = p.programs[0]
        assert prog.name == "hello_world" and prog.ancestor is None
        assert prog.is_program and not prog.is_subroutine
        assert len(prog.preamble) == 2
        assert prog.preamble[0] == "! a very simple example"
        assert "foo" in prog.children
        assert prog.to_string(verbose=2)

        assert len(prog.contains) == 2
        sub = prog.contains[0]
        assert sub.is_subroutine and sub.name == "foo" and sub.ancestor is prog
        assert sub.to_string(verbose=2)

        func = prog.contains[1]
        assert func.is_function and func.name == "bar" and func.ancestor is prog
        assert func.to_string(verbose=2)

    def test_simple_module(self):
        """Parsing Fortran module written following good programming standards."""
        s = """
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

! A Very simple example
! of Fortran module written following good programming standards

module m_crystal
    use defs_basis, only : dp
    use m_fstrings,  only : int2char10

    implicit none
    private

    type,public :: crystal_t
        integer :: nsym
        real(dp),allocatable :: xred(:,:)
    end type

    public :: crystal_free
    public :: idx_spatial_inversion

    interface operator (==)
      module procedure coeffs_compare
    end interface operator (==)

CONTAINS  !=============================
!!***

subroutine crystal_free(cryst)
    type(crystal_t),intent(inout) :: cryst
    if (allocated(cryst%zion)) then
      ABI_FREE(cryst%zion)
    end if
end subroutine crystal_free

pure function idx_spatial_inversion(cryst) result(inv_idx)
    integer :: inv_idx
    type(crystal_t),intent(in) :: cryst
    if (cryst%nsym > 1)
end function idx_spatial_inversion

logical function private_function(cryst) result(ans)
    type(crystal_t),intent(in) :: cryst
    ans = (cryst%nsym > 1)
end function idx_spatial_inversion

end module m_crystal
"""
        p = FortranKissParser()
        p.parse_lines(s.splitlines())
        assert len(p.modules) == 1
        assert not p.programs and not p.subroutines
        assert "abi_common.h" in p.includes
        assert "defs_basis" in p.uses
        assert "m_fstrings" in p.uses

        mod = p.modules[0]
        assert mod.name == "m_crystal" and mod.ancestor is None
        assert mod.is_module and not mod.is_subroutine
        assert len(mod.preamble) == 2
        assert mod.preamble[0] == "! a very simple example"
        assert not mod.children
        assert mod.to_string(verbose=2)

        assert len(mod.contains) == 3
        sub = mod.contains[0]
        assert sub.is_subroutine and sub.name == "crystal_free" and sub.ancestor is mod

        func = mod.contains[1]
        assert func.is_function and func.name == "idx_spatial_inversion" and func.ancestor is mod

    def test_tricky_module(self)
        """Parsing tricky Fortran module similar to 42_libpaw/m_libpaw_libxc.F90"""
        s = """
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
""""
