# coding: utf-8
"""Tests abiref.bib file."""
from __future__ import division, print_function, unicode_literals, absolute_import

from .tools import patch_syspath, AbimkdocsTest
patch_syspath()

import os

from abimkdocs.website import MyEntry


class BibrefTest(AbimkdocsTest):

    def test_bibref(self):
        bibpath = os.path.join(os.path.dirname(__file__), "..", "doc", "abiref.bib")
        assert os.path.isfile(bibpath)

        # Get bibtex references and cast to MyEntry instance.
        from pybtex.database import parse_file
        bib_data = parse_file(bibpath, bib_format="bibtex")
        for entry in bib_data.entries.values():
            entry.__class__ = MyEntry

        # Mandatory fields
        type2fields = dict(
            eprint=("journal", "archivePrefix", "eprint", "year"),
            article=("journal", "volume", "pages", "year"),
            book=("publisher", "year", "isbn"),
            phdthesis=("school", "year"),
            misc=("year"),
        )
        type2fields["incollection"] = type2fields["book"]
        type2fields["mastersthesis"] = type2fields["phdthesis"]

        def validate_entry(entry):
            print("Testing bibtex key `%s` of type `%s`" % (entry.key, entry.type))
            fields = entry.fields
            assert "title" in fields and fields["title"]
            assert "author" in entry.persons and len(entry.persons["author"]) > 0
            for f in type2fields[entry.type]:
                assert f in fields
                #assert fields[f]

            if entry.type in ("article"):
                assert "url" in fields or "doi" in fields
                if "url" in fields: assert fields["url"]
                if "doi" in fields: assert fields["doi"]

        for key, entry in bib_data.entries.items():
            # TODO validate_entry(entry)
            assert entry.to_abimarkdown()
            assert entry.to_html()
            assert entry.to_bibtex()
            assert entry.get_bibtex_btn_modal(link=False)
