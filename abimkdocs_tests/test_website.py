# coding: utf-8
"""Tests abiref.bib file."""
from __future__ import division, print_function, unicode_literals, absolute_import

from .tools import patch_syspath, AbimkdocsTest
patch_syspath()

import os
from abimkdocs.website import Website


class WebsiteTest(AbimkdocsTest):

    def test_website(self):
        dirpath = os.path.join(os.path.dirname(__file__), "..", "..", "doc")
        website = Website.build("./doc", deploy=False, verbose=10)
        assert website is website.get()
        assert not website.warnings
        #website.find_unreferenced_mds()

        # Test (Abinit) wikilink syntax: [namespace:name#fragment|text]
        # URLs are supposed to be relative to page_rpath.
        def element(token, page_rpath="/tutorial/base1.md"):
            return website.get_wikilink(token, page_rpath)

        e = element("https://www.abinit.org|Abinit website")
        assert e.get("href") == "https://www.abinit.org" and e.text == "Abinit website"
        e = element("#internal_link|text")
        assert e.get("href") == "#internal_link" and e.text == "text"
        e = element("lesson:base1|base1")
        assert e.get("href") == "base1" and e.text == "base1"
        e = element("tutorial:base1|base1")
        assert e.get("href") == "base1" and e.text == "base1"
        # This requires howto_topic
        #e = element("topic_SelfEnergy|self-energy")
        #assert e.get("href") == "..//topics/SelfEnergy" and e.text == "self-energy"
        e = element("help:abinit|See Abinit help")
        assert e.get("href") == "../guide/abinit" and e.text == "See Abinit help"
        e = element("dipdip@anaddb|See anaddb")
        assert e.get("href") == "../variables/anaddb#dipdip" and e.text == "See anaddb"
        e = element("dipdip@anaddb")
        assert e.get("href") == "../variables/anaddb#dipdip" and e.text == "dipdip@anaddb"
        e = element("ecut")
        assert e.get("href") == "../variables/basic#ecut" and e.text == "ecut"

        e = element("cite:Gonze2009")
        assert e.get("href") == "../theory/bibliography#gonze2009" and e.text == "[Gonze2009]"

        e = element("~abinit/tests/v1/Input/t01.in|t01.in")
        assert e.get("href") == "../tests/v1/Input/t01.in" and e.text == "t01.in"
        e = element("~abinit/tests/Psps_for_tests/6c.lda.atompaw|6c.paw")
        assert e.get("href") == "../tests/Psps_for_tests/6c.lda.atompaw" and e.text == "6c.paw"

        e = element("ENERGY")
        #assert e.get("href") == "/guide/abinit#32-more-about-abinit-input-variables" and e.text == "ENERGY"
        e = element("AUTO_FROM_PSP")
        assert e.get("href") == "../variables/external_parameters#auto_from_psp" and e.text == "AUTO_FROM_PSP"

        # Wikilinks with namespace.
        e = element("anaddb:asr")
        assert e.get("href") == "../variables/anaddb#asr" and e.text == "asr"
        e = element("asr@anaddb")
        assert e.get("href") == "../variables/anaddb#asr" and e.text == "asr@anaddb"
        e = element("lesson:wannier90|w90")
        assert e.get("href") == "wannier90" and e.text == "w90"
        e = element("tutorial:wannier90|w90")
        assert e.get("href") == "wannier90" and e.text == "w90"
        e = element("help:abinit|Abinit help")
        assert e.get("href") == "../guide/abinit" and e.text == "Abinit help"
        # TODO howto_topic
        #e = element("topic:BSE|BSE topic")
        #assert e.get("href") == "../topics/bse" and e.text == "BSE topic"
        e = element("cite:Amadon2008|Read this")
        assert e.get("href") == "../theory/bibliography#amadon2008" and e.text == "Read this"
        e = element("theory:mbt|GW Notes")
        assert e.get("href") == "../theory/mbt" and e.text == "GW Notes"
        e = element("varset:allvars|All vars")
        assert e.get("href") == "../variables" and e.text == "All vars"
        e = element("varset:bse|BSE varset")
        assert e.get("href") == "../variables/bse" and e.text == "BSE varset"

        e = element("test:libxc_41")
        assert e.get("href") == "../tests/libxc/Input/t41.in" and e.text == "libxc[41]"
        e = element("src:94_scfcv/scfcv.F90")
        assert e.get("href") == "https://github.com/abinit/abinit/blob/master/src/94_scfcv/scfcv.F90"
        e = element("ac:abiref_gnu_9.2_debug.ac")
        assert e.get("href") == "../abichecks/buildsys/Refs/abiref_gnu_9.2_debug.ac"
        e = element("pdf:howto_chebfi.pdf|chebfi")
        #assert e.get("href") == "/build/config-examples/abiref_gnu_9.2_debug.ac" and e.text == "chebfi"
        #e = element("[gitsha:f74dba1ed8346ca586dc95fd10fe4b8ced108d5e]")
        assert not website.warnings
