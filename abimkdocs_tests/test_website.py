"""Tests abiref.bib file."""

from pathlib import Path

from .tools import AbimkdocsTest, patch_syspath

patch_syspath()

from abimkdocs.website import Website, find_pdf_basename_collisions


class WebsiteTest(AbimkdocsTest):

    def test_website(self):
        repository_root = Path(__file__).resolve().parents[1]
        website = Website.build(str(repository_root / "doc"), deploy=False, verbose=10)
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
        assert e.get("href") == "" and e.text == "base1"
        e = element("tutorial:base1|base1")
        assert e.get("href") == "" and e.text == "base1"
        # This requires howto_topic
        #e = element("topic_SelfEnergy|self-energy")
        #assert e.get("href") == "..//topics/SelfEnergy" and e.text == "self-energy"
        e = element("help:abinit|See Abinit help")
        assert e.get("href") == "../guide/abinit.md" and e.text == "See Abinit help"
        e = element("dipdip@anaddb|See anaddb")
        assert e.get("href") == "../variables/anaddb.md#dipdip" and e.text == "See anaddb"
        e = element("dipdip@anaddb")
        assert e.get("href") == "../variables/anaddb.md#dipdip" and e.text == "dipdip@anaddb"
        e = element("ecut")
        assert e.get("href") == "../variables/basic.md#ecut" and e.text == "ecut"

        e = element("cite:Gonze2009")
        assert e.get("href") == "../theory/bibliography.md#gonze2009" and e.text == "[Gonze2009]"

        e = element("~abinit/tests/v1/Input/t01.in|t01.in")
        assert e.get("href") == "../tests/v1/Input/t01.in" and e.text == "t01.in"
        e = element("~abinit/tests/Pspdir/6c.lda.atompaw|6c.paw")
        assert e.get("href") == "../tests/Pspdir/6c.lda.atompaw" and e.text == "6c.paw"

        e = element("ENERGY")
        assert e.get("href") == "../guide/abinit.md#parameters" and e.text == "ENERGY"
        e = element("AUTO_FROM_PSP")
        assert e.get("href") == "../variables/external_parameters.md#auto_from_psp" and e.text == "AUTO_FROM_PSP"

        # Wikilinks with namespace.
        e = element("anaddb:asr")
        assert e.get("href") == "../variables/anaddb.md#asr" and e.text == "asr"
        e = element("asr@anaddb")
        assert e.get("href") == "../variables/anaddb.md#asr" and e.text == "asr@anaddb"
        e = element("lesson:wannier90|w90")
        assert e.get("href") == "wannier90.md" and e.text == "w90"
        e = element("tutorial:wannier90|w90")
        assert e.get("href") == "wannier90.md" and e.text == "w90"
        e = element("help:abinit|Abinit help")
        assert e.get("href") == "../guide/abinit.md" and e.text == "Abinit help"
        # TODO howto_topic
        #e = element("topic:BSE|BSE topic")
        #assert e.get("href") == "../topics/bse" and e.text == "BSE topic"
        e = element("cite:Amadon2008|Read this")
        assert e.get("href") == "../theory/bibliography.md#amadon2008" and e.text == "Read this"
        e = element("theory:mbt|GW Notes")
        assert e.get("href") == "../theory/mbt.md" and e.text == "GW Notes"
        e = element("varset:allvars|All vars")
        assert e.get("href") == "../variables/index.md" and e.text == "All vars"
        e = element("varset:bse|BSE varset")
        assert e.get("href") == "../variables/bse.md" and e.text == "BSE varset"

        e = element("test:libxc_41")
        assert e.get("href") == "../tests/libxc/Input/t41.abi" and e.text == "libxc[41]"
        e = element("src:94_scfcv/scfcv.F90")
        assert e.get("href") == "https://github.com/abinit/abinit/blob/master/src/94_scfcv/scfcv.F90"
        # FIXME: buildsys refs are not needed anymore (YP)
        #e = element("ac:abiref_gnu_9.2_debug.ac")
        #assert e.get("href") == "../abichecks/buildsys/Refs/abiref_gnu_9.2_debug.ac"
        # doc/topics/documents/howto_chebfi.pdf is now the only copy of this
        # file (doc/theory/howto_chebfi.pdf, a stale duplicate, was removed:
        # see find_pdf_basename_collisions()/Website.pdfs for why having two
        # PDFs share a basename is a real, silently-machine-dependent bug,
        # not just untidiness).
        e = element("pdf:howto_chebfi.pdf|chebfi")
        assert e.get("href") == "../topics/documents/howto_chebfi.pdf" and e.text == "chebfi"
        #e = element("[gitsha:f74dba1ed8346ca586dc95fd10fe4b8ced108d5e]")
        assert not website.warnings


class FindPdfBasenameCollisionsTest(AbimkdocsTest):
    """Unit tests for find_pdf_basename_collisions(), isolated from the
    (expensive) full Website.build() -- no need to construct a real site to
    exercise this pure function.
    """

    def test_no_collision_when_every_basename_is_unique(self):
        pairs = [("a.pdf", "/doc/x/a.pdf"), ("b.pdf", "/doc/y/b.pdf")]
        assert find_pdf_basename_collisions(pairs) == []

    def test_detects_a_collision_and_names_every_conflicting_path(self):
        pairs = [
            ("howto_chebfi.pdf", "/doc/theory/howto_chebfi.pdf"),
            ("howto_chebfi.pdf", "/doc/topics/documents/howto_chebfi.pdf"),
            ("unrelated.pdf", "/doc/x/unrelated.pdf"),
        ]
        messages = find_pdf_basename_collisions(pairs)

        assert len(messages) == 1
        msg = messages[0]
        assert "howto_chebfi.pdf" in msg
        assert "/doc/theory/howto_chebfi.pdf" in msg
        assert "/doc/topics/documents/howto_chebfi.pdf" in msg
        assert "unrelated.pdf" not in msg

    def test_names_the_alphabetically_last_path_as_the_winner(self):
        """The message must name whichever path OrderedDict(sorted(pairs))
        actually keeps -- the alphabetically-last one -- not just list the
        conflicting paths without saying which one wins.
        """
        pairs = [
            ("dup.pdf", "/doc/theory/dup.pdf"),
            ("dup.pdf", "/doc/topics/documents/dup.pdf"),
        ]
        messages = find_pdf_basename_collisions(pairs)

        assert len(messages) == 1
        assert "`/doc/topics/documents/dup.pdf`" in messages[0]

    def test_handles_more_than_two_colliding_files(self):
        pairs = [
            ("dup.pdf", "/doc/a/dup.pdf"),
            ("dup.pdf", "/doc/b/dup.pdf"),
            ("dup.pdf", "/doc/c/dup.pdf"),
        ]
        messages = find_pdf_basename_collisions(pairs)

        assert len(messages) == 1
        assert "Found 3 PDF files" in messages[0]
        for path in ("/doc/a/dup.pdf", "/doc/b/dup.pdf", "/doc/c/dup.pdf"):
            assert path in messages[0]
