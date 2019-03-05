#!/usr/bin/python3

import unittest

from abinit.buildsys.config import AbinitConfigParser

# Generic config template for basic tests
my_config_text = """\
[DEFAULT]
var0 = value0_default

[section1]
var1 = value1
VAR1u = value1U

[Section2]
var2 = value2
Var2m = value2M

[SECTION3]
var3 = value3
var0 = value0_section3
"""

class TestConfigParser(unittest.TestCase):


    def test_parser_preserves_section_case(self):

        cnf = AbinitConfigParser(my_config_text)
        titles = cnf.sections()
        assert ( ("section1" in titles) and ("Section2" in titles) and \
                 ("SECTION3" in titles) )


    def test_parser_preserves_varname_case(self):

        cnf = AbinitConfigParser(my_config_text)
        names = [item for item, blob in cnf.items("section1")] + \
            [item for item, blob in cnf.items("Section2")]
        assert ( ("VAR1u" in names) and ("Var2m" in names) )

