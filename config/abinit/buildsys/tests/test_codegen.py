#!/usr/bin/python3

import re
import unittest

from abinit.buildsys.codegen import AbinitTemplate

# Generic config template for basic tests
my_template_text = """\
AC_LANG_PUSH([@%language%@])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
  [@%includes%@],
  [[
    @%source%@
  ]])],
  [@%action_ok%@],
  [@%action_fail%@])
AC_LANG_POP([@%language%@])
"""
my_patterns = ["action_fail", "action_ok", "includes", "language", "source"]

class TestTemplate(unittest.TestCase):

    def test_template_sets_correct_patterns(self):

        tpl = AbinitTemplate(my_template_text)
        pat = tpl.get_patterns()
        chk = [item for item in my_patterns if not item in pat] + \
            [item for item in pat if not item in my_patterns]

        self.assertTrue( (len(pat) == len(my_patterns)) and (len(chk) == 0) )


    def test_template_accepts_any_delimiters(self):

        tpl = AbinitTemplate("01pattern1234\n01other_pattern234",
            start_tag="01", end_tag="234")
        pat = tpl.get_patterns()

        self.assertTrue( (len(pat) == 2) and (pat[0] == "other_pattern") and \
            (pat[1] == "pattern1") )


    def test_template_rejects_simple_text(self):

        with self.assertRaises(ValueError):
            AbinitTemplate("Simple\ntext.")


    def test_template_get_missing_accepts_dict(self):

        tpl = AbinitTemplate(my_template_text)
        miss = tpl.get_missing({"language":"Fortran"})

        self.assertTrue( len(miss) == 4 )


    def test_template_returns_correct_missing(self):

        tpl = AbinitTemplate(my_template_text)
        miss = tpl.get_missing(
            ["action_fail", "action_ok", "includes", "source"])

        self.assertTrue( (len(miss) == 1) and (miss[0] == "language") )


    def test_template_get_undefined_accepts_dict(self):

        tpl = AbinitTemplate(my_template_text)
        undef = tpl.get_undefined({"language":"Fortran"})

        self.assertTrue( undef is None )


    def test_template_returns_correct_undefined(self):

        tpl = AbinitTemplate(my_template_text)
        undef = tpl.get_undefined(["language", "something"])

        self.assertTrue( len(undef) == 1 )


    def test_template_substitute_stops_on_missing(self):

        tpl = AbinitTemplate(my_template_text)
        rep = {
            "action_fail":"result='no'",
            "action_ok":"result='yes'",
            "includes":"",
            "source":"call some_routine(some_param)"}

        with self.assertRaises(ValueError):
            tpl.substitute(rep, err_miss=True, err_undef=False)


    def test_template_substitute_stops_on_undefined(self):

        tpl = AbinitTemplate(my_template_text)
        rep = {
            "action_fail":"result='no'",
            "action_ok":"result='yes'",
            "includes":"",
            "source":"call some_routine(some_param)",
            "something":"nothing"}

        with self.assertRaises(ValueError):
            tpl.substitute(rep, err_miss=False, err_undef=True)


    def test_template_substitute_passes_when_asked(self):

        tpl = AbinitTemplate(my_template_text)
        rep = {
            "action_fail":"result='no'",
            "action_ok":"result='yes'",
            "includes":"",
            "source":"call some_routine(some_param)",
            "something":"nothing"}
        txt = tpl.substitute(rep, err_miss=False, err_undef=False)

        self.assertTrue( re.search("@%language%@", txt, flags=re.MULTILINE) )

