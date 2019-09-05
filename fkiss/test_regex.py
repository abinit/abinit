"""
Units Tests for regular expressions.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

from unittest import TestCase

from .regex import HasRegex


class TestRegext(TestCase):

    def test_regex(self):
        """Test regular expressions."""
        p = HasRegex

        # A quoted expression.
        assert p.RE_QUOTED.match("'hello'")
        assert p.RE_QUOTED.match('""')

        m = p.RE_F90COMMENT.match("  ! hello")
        assert m and m.group("value") == " hello"

        assert p.RE_OMP_SENTINEL.match("!$OMP PARALLEL DO")

        # Program
        m = p.RE_PROG_START.match("program foo")
        assert m and m.group("name") == "foo"
        m = p.RE_PROC_END.match("end program foo")
        assert m and m.group("proc_type") == "program" and m.group("name") == "foo"

        # Module
        m = p.RE_MOD_START.match("module   foo")
        assert m and m.group("name") == "foo"
        m = p.RE_MOD_START.match("module   foo!hello")
        assert m and m.group("name") == "foo"
        m = p.RE_MOD_END.match("end   module   foo")
        assert m and m.group("name") == "foo"

        # Subroutine
        m = p.RE_SUB_START.match("subroutine  foo")
        assert m and m.group("name") == "foo" and not m.group("prefix")
        m = p.RE_SUB_START.match("subroutine foo(a, b)")
        assert m and m.group("name") == "foo" and m.group("prefix") == ""
        m = p.RE_SUB_START.match("pure  subroutine foo(a, b)")
        assert m and m.group("name") == "foo" and m.group("prefix").strip() == "pure"
        #m = p.RE_SUB_END.match("end  subroutine foo!hello")
        #assert m and m.group("name") == "foo"
        m = p.RE_PROC_END.match("end  subroutine foo!hello")
        assert m and m.group("proc_type") == "subroutine" and m.group("name") == "foo"

        # Functions: these are tricky so test them carefully.
        m = p.RE_FUNC_START.match("function  foo")
        assert m and m.group("name") == "foo" and not m.group("prefix")
        m = p.RE_FUNC_START.match("integer function foo(a, b)")
        assert m and m.group("name") == "foo" and m.group("prefix").strip() == "integer"
        m = p.RE_FUNC_START.match("integer(kind = i8b) function foo(a, b)")
        assert m and m.group("name") == "foo" and m.group("prefix").strip() == "integer(kind = i8b)"
        m = p.RE_FUNC_START.match("real(dp) function foo(a, b)")
        assert m and m.group("name") == "foo" and m.group("prefix").strip() == "real(dp)"
        m = p.RE_FUNC_START.match("pure complex(dpc) function foo(a)")
        assert m and m.group("name") == "foo" and m.group("prefix").strip() == "pure complex(dpc)"
        m = p.RE_FUNC_START.match("pure logical function foo(a, b)")
        assert m and m.group("name") == "foo" and m.group("prefix").strip() == "pure logical"
        m = p.RE_FUNC_START.match("type ( foo_t ) function foo(a, b) result(res)")
        assert m and m.group("name") == "foo" and m.group("prefix").strip() == "type ( foo_t )"

        m = p.RE_FUNC_END.match("end  function foo  !bar")
        assert m and m.group("name") == "foo"
        #m = p.RE_FUNC_END.match("end")
        #assert m and m.group("name") is None

        # Calls to subroutines
        m = p.RE_SUBCALL.match("call  foo")
        assert m and m.group("name") == "foo"
        m = p.RE_SUBCALL.match("call  foo ( a, b)")
        assert m and m.group("name") == "foo"
        m = p.RE_SUBCALL.match("if (a == 1) call foo(a=2)")
        assert m and m.group("name") == "foo"
        #m = p.RE_SUBCALL.match("call obj%foo(a=1, b=3)")
        #assert m and m.group("name") == "obj%foo" and m.group("method") == "foo"
        #m = p.RE_SUBCALL.match("call obj % foo")
        #assert m and m.group("name") == "obj%foo" and m.group("method") == "foo"
        #m = p.RE_SUBCALL.match("if (a == b) call obj%foo(a, b)")
        #assert m and m.group("name") == "obj%foo" and m.group("method") == "foo"

        # Interface
        m = p.RE_INTERFACE_END.match("end interface")
        assert m and not m.group("name")
        m = p.RE_INTERFACE_START.match("abstract interface foo")
        assert m and m.group("name") == "foo"
        m = p.RE_INTERFACE_END.match("end  interface  foo")
        assert m and m.group("name") == "foo"

        # Intrinsic type declarations.
        m = p.RE_CHARACTER_DEC.match("character (len= 10)")
        assert m and m.group("len") == "10"
        m = p.RE_CHARACTER_DEC.match("character(len=fnlen)")
        assert m and m.group("len") == "fnlen"
        m = p.RE_CHARACTER_DEC.match("character(len=*),intent(in)")
        assert m and m.group("len") == "*"

        m = p.RE_INTENT.search("integer, intent( inout )")
        assert m and m.group("value") == "inout"

        assert not p.RE_TYPECLASS_DEC.match("type, public :: foo_t")
        m = p.RE_TYPECLASS_DEC.match("type( foo_t) :: foo")
        assert m and m.group("ftype") == "type" and m.group("name") == "foo_t"
        m = p.RE_TYPECLASS_DEC.match("type(pawrhoij_type ),allocatable :: pawrhoij1_i1pert(:)")
        assert m and m.group("ftype") == "type" and m.group("name") == "pawrhoij_type"
        m = p.RE_TYPECLASS_DEC.match("class ( Circle ), intent(in) :: this")
        assert m and m.group("ftype") == "class" and m.group("name") == "Circle"

        # integer
        m = p.RE_NUMBOOL_DEC.match("integer :: timrev")
        assert m and m.group("ftype") == "integer" and m.group("kind") is None
        m = p.RE_NUMBOOL_DEC.match("integer ( i8b ) :: timrev")
        assert m and m.group("ftype") == "integer" and m.group("kind") == "i8b"
        m = p.RE_NUMBOOL_DEC.match("integer(kind=4) :: timrev")
        assert m and m.group("ftype") == "integer" and m.group("kind") == "4"
        m = p.RE_NUMBOOL_DEC.match("integer ( kind = i4b ) :: timrev")
        assert m and m.group("ftype") == "integer" and m.group("kind") == "i4b"
        # real
        m = p.RE_NUMBOOL_DEC.match("double precision :: ucvol")
        assert m and m.group("ftype") == "double precision" and m.group("kind") is None
        m = p.RE_NUMBOOL_DEC.match("double complex :: j")
        assert m and m.group("ftype") == "double complex" and m.group("kind") is None
        m = p.RE_NUMBOOL_DEC.match("real(dp) :: ucvol")
        assert m and m.group("ftype") == "real" and m.group("kind") == "dp"
        # boolean
        m = p.RE_NUMBOOL_DEC.match("logical :: use_antiferro")
        assert m and m.group("ftype") == "logical" and m.group("kind") is None

        # datatype declaration.
        m = p.RE_TYPE_START.match("type foo_t")
        assert m and m.group("name") == "foo_t" and not m.group("attribs").strip()
        m = p.RE_TYPE_START.match("type, public :: foo_t")
        assert m and m.group("name") == "foo_t" and "public" in m.group("attribs")
        m = p.RE_TYPE_START.match("type, public, foobar :: foo_t")
        assert m and m.group("name") == "foo_t" and "foobar" in m.group("attribs")
        assert not p.RE_TYPE_START.match("type(foo_t) :: foo")

        m = p.RE_PUB_OR_PRIVATE.match("public ")
        assert m and m.group("name") == "public"
        m = p.RE_PUB_OR_PRIVATE.match("private ! everything private ")
        assert m and m.group("name") == "private"
        assert not p.RE_PUB_OR_PRIVATE.match("private = 1 ! this is not private ")

        assert p.RE_CONTAINS.match("contains  ")
        assert p.RE_CONTAINS.match("contains!foo ")
        assert not p.RE_CONTAINS.match("contains=1")

        # Test continuation lines
        m = p.RE_CONTLINE_START.match("subroutine (a, &")
        assert m and m.group("prefix") == "subroutine (a, " and not m.group("postfix")
        m = p.RE_CONTLINE_START.match("subroutine (a, & !&hello")
        assert m and m.group("prefix") == "subroutine (a, " and m.group("postfix") == " !&hello"

        assert not p.RE_CONTLINE_START.match("& call foo(a)")
        m = p.RE_CONTLINE_START.match("call foo('me & you')")
        assert m and not m.group("postfix").strip().startswith("!")
        # This will fool the script
        #m = p.RE_CONTLINE_START.match("call foo('me &! you')")
        assert not p.RE_CONTLINE_START.match("! call &")

        m = p.RE_CONTLINE_NEXT.match("& call &")
        assert m and m.group("value") == " call "
        assert p.RE_CONTLINE_NEXT.match("& call &!debug")
        assert p.RE_CONTLINE_NEXT.match("& end ")
