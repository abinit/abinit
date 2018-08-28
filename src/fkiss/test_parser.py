"""
Units Tests for Parser
"""
from __future__ import print_function, division, unicode_literals, absolute_import

from unittest import TestCase

from .parser import FortranKissParser
#from .project import AbinitProject


class TestFortranKissParser(TestCase):

    def test_simple_program(self):
        """Parsing Fortran program written following good programming standards."""
        s = """
#include "abi_common.h"

! A Very simple example
! of Fortran program written following good programming standards

program hello_world!my first test
    use defs_basis
    implicit none
    call foo()
    if (.True.) call foobar()

contains! contained procedure

  subroutine foo()!hello
  end subroutine foo

function bar()
    end function bar
end program
"""
        p = FortranKissParser(verbose=2)
        p.parse_string(s)
        assert len(p.programs) == 1
        assert not p.modules
        assert not p.subroutines
        assert "abi_common.h" in p.all_includes
        assert "defs_basis" in p.all_uses

        prog = p.programs[0]
        assert prog.name == "hello_world" and prog.ancestor is None
        assert prog.is_program and not prog.is_subroutine
        assert prog.to_string(verbose=2)
        preamble = prog.preamble.splitlines()
        assert len(preamble) == 2
        assert preamble[0] == "! A Very simple example"
        assert "foo" in prog.children
        assert "foobara" in prog.children

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
        ! nsym docstring
        real(dp),allocatable :: xred(:,:)
        ! xred
        ! docstring

    end type crystal_t

    public :: crystal_free
    public :: idx_spatial_inversion

    interface get_unit
      module procedure get_free_unit
      module procedure get_unit_from_fname

    end interface

    interface operator (==)
      module procedure coeffs_compare
    end interface operator (==)

CONTAINS  !=============================
!!***

subroutine crystal_free(cryst)
    type(crystal_t),intent(inout) :: cryst

    ! interfaces inside routines are not parsed
    interface!but the parser should not crash here.
          function integrand(x)
          integer, parameter :: dp=kind(1.0d0)
          real(dp) integrand
          real(dp),intent(in) :: x
          end function integrand
end interface

    if (allocated(cryst%zion)) then
      ABI_FREE(cryst%zion)
    end if
end subroutine crystal_free

pure function idx_spatial_inversion(cryst) result(inv_idx)
    integer :: inv_idx
    type(crystal_t),intent(in) :: cryst
    if (cryst%nsym > 1)
end function idx_spatial_inversion

logical function private_bool_function(cryst)
    logical :: private_bool_function
    type(crystal_t),intent(in) :: cryst
    private_bool_function = (cryst%nsym > 1)
end function idx_spatial_inversion

end module m_crystal
"""
        p = FortranKissParser()
        p.parse_string(s)
        assert len(p.modules) == 1
        assert not p.programs and not p.subroutines
        assert "abi_common.h" in p.all_includes
        assert "defs_basis" in p.all_uses
        assert "m_fstrings" in p.all_uses

        mod = p.modules[0]
        assert mod.name == "m_crystal" and mod.ancestor is None
        assert mod.is_module and not mod.is_subroutine
        preamble = mod.preamble.splitlines()
        assert len(preamble) == 2
        assert preamble[0] == "! A Very simple example"
        assert not mod.children
        assert mod.to_string(verbose=2)

        assert len(mod.contains) == 3
        sub = mod.contains[0]
        assert sub.is_subroutine and sub.name == "crystal_free" and sub.ancestor is mod

        func = mod.contains[1]
        assert func.is_function and func.name == "idx_spatial_inversion" and func.ancestor is mod
        # TODO
        #assert func.is_pure
        #assert func.result.name == "inv_idx"
        #assert func.is_public
        #assert func.result.ftype == "integer" and func.result.kind is None

        func = mod.contains[2]
        assert func.is_function and func.name == "private_bool_function" and func.ancestor is mod
        #assert not func.is_pure
        #assert func.result.name == func.name
        #assert not func.is_public and func.is_provate
        #assert func.result.ftype == "logical" and func.result.kind is None

        assert len(mod.interfaces) == 2
        assert mod.interfaces[0].name == "get_unit"
        # TODO
        #assert mod.interfaces[1].name == "(==)"

        # Test datatype
        assert len(mod.types) == 1
        crystal_t = mod.types[0]
        assert crystal_t.name == "crystal_t"
        assert crystal_t.ancestor.is_module  and crystal_t.ancestor.name == "m_crystal"
        str(crystal_t); repr(crystal_t)
        assert crystal_t.to_string(verbose=2)

        crystal_t.analyze()
        nsym, xred = crystal_t.variables["nsym"], crystal_t.variables["xred"]
        str(nsym); repr(nsym)
        assert nsym.ftype == "integer" and nsym.is_scalar and not nsym.is_pointer
        assert nsym.doc == "! nsym docstring"
        # TODO
        #assert xred.ftype == "real" and xred.kind == "dp"
        assert not xred.is_scalar and xred.is_array and xred.shape == "(:,:)"
        assert xred.doc.splitlines() == ["! xred", "! docstring"]

    def test_datatype(self):
        """Parsing Fortran datatype declaration"""
        s = """\
module foo_module!my first module

type, public :: foo_t

  integer ( i8b ) , public :: nsym
  ! number of symmetries

  real (kind = 4), private , allocatable:: xred(:, :)

  real (kind = 4), dimension(:,:) , allocatable:: xcart

  type(foo_t), pointer :: foo
  ! Another foo type

  double complex,allocatable::cvals(:)!inlined

  character (len=fnlen), allocatable , private :: strings(:)

  real(dp), allocatable :: ucvol, ecut  ! doc for ucvol and ecut

end type foo_t!foo_t

end module foo_module
"""
        p = FortranKissParser()
        p.parse_string(s)
        dt = p.modules[0].types[0]
        assert dt.name == "foo_t"
        dt.analyze()

        # Test nsym
        nsym = dt.variables["nsym"]
        assert nsym.ftype == "integer" and nsym.kind == "i8b"
        assert nsym.doc == "! number of symmetries"

        # Test xred && xcart
        # TODO: xcart
        for vname in ["xred",]: # "xcart"]:
            var = dt.variables[vname]
            assert var.ftype == "real" and var.kind == "4" and not var.doc
            assert var.is_allocatable and not var.is_scalar and var.is_array

        # Test foo
        foo = dt.variables["foo"]
        assert foo.is_pointer
        assert foo.ftype == "type" and foo.kind == "foo_t"
        assert foo.is_scalar
        #assert foo.is_pointer_set_to_null

        # Test cvals
        cvals = dt.variables["cvals"]
        assert cvals.ftype == "doublecomplex" and cvals.kind is None

        # Test strings
        strings = dt.variables["strings"]
        assert strings.ftype == "character"
        #assert strings.strlen == "fnlen"
        #assert len(strings.shape) == 1

        ucvol, ecut = dt.variables["ucvol"], dt.variables["ecut"]
        assert ucvol.is_scalar and ecut.is_scalar
        #assert ucvol.doc == "! doc for ucvol and ecut"
        assert ucvol.doc == ecut.doc
