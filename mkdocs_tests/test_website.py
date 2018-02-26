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
        #websitefind_unreferenced_mds(self):

        # Test (Abinit) wikilink syntax: [namespace:name#fragment|text]
        # URLs are supposed to be relative to page_rpath.
        def element(token, page_rpath="/mkdocs-tutorial/base1.md"):
            return website.get_wikilink(token, page_rpath)

        e = element("https://www.abinit.org|Abinit website")
        assert e.get("href") == "https://www.abinit.org" and e.text == "Abinit website"
        e = element("#internal_link|text")
        assert e.get("href") == "#internal_link" and e.text == "text"
        e = element("lesson_base1|base1")
        assert e.get("href") == "../tutorials/base1" and e.text == "base1"
        # This requires howto_topic
        #e = element("topic_SelfEnergy|self-energy")
        #assert e.get("href") == "..//mkdocs-topics/SelfEnergy" and e.text == "self-energy"
        e = element("help_abinit|See Abinit help")
        assert e.get("href") == "../mkdocs-user-guide/abinit" and e.text == "See Abinit help"
        e = element("dipdip@anaddb|See anaddb")
        assert e.get("href") == "../mkdocs-variables/anaddb#dipdip" and e.text == "See anaddb"
        e = element("dipdip@anaddb")
        assert e.get("href") == "../mkdocs-variables/anaddb#dipdip" and e.text == "dipdip@anaddb"
        e = element("ecut")
        assert e.get("href") == "../mkdocs-variables/basic#ecut" and e.text == "ecut"

        e = element("Gonze2009")
        assert e.get("href") == "../mkdocs-theory/bibliography#gonze2009" and e.text == "[Gonze2009]"

        e = element("~/abinit/tests/v1/Input/t01.in|t01.in")
        # FIXME
        #assert e.get("href") == "../tests/v1/input/t01.in" and e.text == "t01.in"
        e = element("~abinit/tests/Psps_for_tests/6c.lda.atompaw|6c.paw")
        #assert e.get("href") == "../tests/psps_for_tests/6c.lda.atompaw" and e.text == "6c.paw"

        e = element("ENERGY")
        #assert e.get("href") == "/mkdocs-user-guide/abinit#32-more-about-abinit-input-variables" and e.text == "ENERGY"
        e = element("AUTO_FROM_PSP")
        assert e.get("href") == "/mkdocs-variables/external_parameters#auto_from_psp" and e.text == "AUTO_FROM_PSP"

        # Wikilinks with namespace.
        e = element("anaddb:asr")
        assert e.get("href") == "/mkdocs-variables/anaddb#%asr" and e.text == "asr@anaddb"
        e = element("lesson:wannier90|w90")
        assert e.get("href") == "/mkdocs-tutorials/wannier90" and e.text == "w90"
        e = element("help:abinit|Abinit help")
        assert e.get("href") == "/mkdocs-user-guide/help_abinit" and e.text == "Abinit help"
        e = element("topic:BSE|BSE topic")
        assert e.get("href") == "/mkdocs-topics/bse" and e.text == "BSE topic"
        # TODO bib --> cite
        e = element("bib:Amadon2008|Read this")
        assert e.get("href") == "/mkdocs-theory/bibliography#amadon2008" and e.text == "Read this"
        e = element("theorydoc:mbpt|GW Notes")
        assert e.get("href") == "/mkdocs-theory/mbpt" and e.text == "GW Notes"
        e = element("varset:allvars|All vars")
        assert e.get("href") == "/mkdocs-variables/index" and e.text == "All vars"
        e = element("varset:BSE|BSE set")
        assert e.get("href") == "/mkdocs-variables/bse" and e.text == "BSE set"

        e = element("test:libxc_41")
        assert e.get("href") == "/tests/libxc/Input/t41.in" and e.text == "libxc[41]"

        e = element("src:94_scfcv/scfcv.F90")
        assert e.get("href") == "https://github.com/abinit/abinit/blob/master/src/94_scfcv/scfcv.F90"

        e = element("ac:abiref_gnu_5.3_debug.ac")
        assert e.get("href") == "/build/config-examples/abiref_gnu_5.3_debug.ac"

        e = element("pdf:howto_chebfi.pdf|chebfi")
        #assert e.get("href") == "/build/config-examples/abiref_gnu_5.3_debug.ac" and e.text == "chebfi"

        #e = element("[gitsha:f74dba1ed8346ca586dc95fd10fe4b8ced108d5e]")

        assert not website.warnings
