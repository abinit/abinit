"""Tests abiref.bib file."""

from .tools import AbimkdocsTest, patch_syspath

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
            eprint=("eprint", "year"),
            article=("journal", "year"),
            book=("publisher", "year"),
            phdthesis=("school", "year"),
            misc=("year",),
            incollection=("year",),
            inproceedings=("year",),
            mastersthesis=("school", "year"),
        )

        def validate_entry(entry):
            print("Testing bibtex key `%s` of type `%s`" % (entry.key, entry.type))
            fields = entry.fields
            assert fields.get("title")
            if entry.type not in ("book", "misc"):
                assert any(entry.persons.get(role) for role in ("author", "editor"))
            for f in type2fields[entry.type]:
                assert f in fields
                #assert fields[f]

            if entry.type == "article":
                # Older bibliography entries may provide neither URL nor DOI.
                if "url" in fields: assert fields["url"]
                if "doi" in fields: assert fields["doi"]

        for key, entry in bib_data.entries.items():
            # Validate mandatory metadata before exercising the output renderers.
            validate_entry(entry)
            assert entry.to_abimarkdown()
            assert entry.to_html()
            assert entry.to_bibtex()
            #assert entry.get_bibtex_btn_modal(link=False)
