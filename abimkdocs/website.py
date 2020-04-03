# coding: utf-8
"""
Classes and functions used to generate the (static) website with the Abinit documentation
from markdown files and the mkdocs static website generator.

For the different between Absolute, Relative, and Root-relative URLs see:

    <http://ifyoucodeittheywill.com/2009/03/absolute-relative-and-root-relative-urls/>
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import io
import time
import re
import shutil
import uuid
import pickle
import yaml
import markdown
import datetime

from collections import OrderedDict, defaultdict
from itertools import groupby
from pprint import pprint
from pybtex.database import parse_file, Entry, BibliographyData
from markdown.util import etree
from pygments import highlight
from pygments.lexers import BashLexer, BibTeXLexer
from pygments.formatters import HtmlFormatter
from doc.tests.pymods.termcolor import cprint
from .variables import lazy_property, Variable,  ABI_TOPICS, ABI_RELEVANCES


def my_unicode(s):
    """Convert string to unicode (needed for py2.7 DOH!)"""
    return unicode(s) if sys.version_info[0] <= 2 else str(s)


def escape(text, tag=None, cls=None):
    """Escape HTML entities in ``text`` string. Enclose new text in ``tag`` if tag with class ``cls``."""
    import html
    text = html.escape(text, quote=True)
    if tag:
        text = '<{tag} class="{cls}">\n{text}\n</{tag}>\n'.format(tag=tag, text=text, cls=cls if cls else "")
    return text


def gen_id(n=1, pre="uuid-"):
    """
    Generate ``n`` universally unique identifiers prepended with ``pre`` string.
    Return string if n == 1 or list of strings if n > 1
    """
    # The HTML4 spec says:
    # ID and NAME tokens must begin with a letter ([A-Za-z]) and may be followed by any number of letters,
    # digits ([0-9]), hyphens ("-"), underscores ("_"), colons (":"), and periods (".").
    if n == 1:
        return pre + str(uuid.uuid4())
    elif n > 1:
        return [pre + str(uuid.uuid4()) for i in range(n)]
    else:
        raise ValueError("n must be > 0 but got %s" % str(n))


def splitall(path):
    """Return list with the components of a ``path``."""
    allparts = []
    while True:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts


def sort_and_groupby(items, key, reverse=False):
    """Sort items using ``key`` function and invoke groupby to group items."""
    return groupby(sorted(items, key=key, reverse=reverse), key=key)


class MyEntry(Entry):
    """
    Extends pybtex Entry with useful methods for generating HTML output.
    See https://bitbucket.org/pybtex-devs/pybtex/
    """
    @lazy_property
    def authors(self):
        """String with authors. Empty if authors are not provided."""
        try:
            #return ", ".join(my_unicode(p) for p in self.persons["author"])
            return ", ".join(my_unicode(p).partition(',')[2] + " " +
                             my_unicode(p).partition(',')[0] for p in self.persons["author"])
        except KeyError:
            return ""

    def to_abimarkdown(self, bibtex_ui="button"):
        """
        Return markdown string with bibliographic entry. Can use Abinit markdown extensions

        Args:
            bibtex_ui: If not None a modal window with the bibtex entry is added.
                Possible values in [None, "link", "button"].
        """
        fields = self.fields
        # Remove {} from (Latex) title.
        # TODO: title must be present
        title = "*%s*" % fields["title"].replace("{", "").replace("}", "")
        authors = self.authors

        # FIXME: enforce format at the level of the unit tests
        if self.type == "article":
            s = '{}  \n{}  \n'.format(authors, title)
            if "eprint" in fields:
                s += "{} **{}**, {} ({})".format(fields["journal"], fields.get("archivePrefix", ""),
                        fields["eprint"], fields["year"])
            else:
                s += "{} **{}**, {} ({})".format(fields["journal"], fields["volume"],
                        fields.get("pages", ""), fields["year"])

        elif self.type in ("book", "inproceedings", "incollection"):
            # FIXME Better treatment for incollection
            #editors = ", ".join(str(e) for e in self.persons["editor"]])
            s = '{}  \n{}  \n'.format(authors, title)
            s += "{} ({})".format(fields["publisher"], fields["year"])
            if "isbn" in fields:
                s += "isbn: %s" % fields["isbn"]

        elif self.type in ("phdthesis", "mastersthesis"):
            s = '{}  \n{}  \n{} ({})'.format(authors, title, fields["school"], fields["year"])

        elif self.type in ("misc", "unpublished"):
            s = '{}  \n{} ({})'.format(authors, title, fields["year"])

        else:
            raise TypeError("Don't know how to convert type: `%s` into markdown string" % self.type)

        s += "  \n"
        if "url" in fields:
            s += 'URL: <a href="{url}" target="_blank">{url}</a><br>'.format(url=fields["url"])
            #s += 'URL: <{url}>  \n'.format(url=fields["url"])
        elif "doi" in fields:
            doi = fields["doi"]
            doi_root = "https://doi.org/"
            if not doi.startswith(doi_root): doi = doi_root + doi
            #s += 'DOI: <{doi}>  \n'.format(doi=doi)
            s += 'DOI: <a href="{doi}" target="_blank">{doi}</a><br>'.format(doi=doi)

        # Add modal window with bibtex button/link.
        if bibtex_ui is not None:
            assert bibtex_ui in ("link", "button")
            btn, modal = self.get_bibtex_btn_modal(link=bibtex_ui=="link")
            s += btn + modal

        return s

    def to_html(self):
        """Return string with entry in HTML format."""
        return markdown.markdown(self.to_abimarkdown())

    def to_bibtex(self):
        """Return the data as a unicode string in the given format."""
        return BibliographyData({self.key: self}).to_string("bibtex")

    def get_bibtex_btn_modal(self, link=False):
        """
        Build HTML string with bootstrap modal and link to open the modal.

        Args:
            link: True if a link instead of a button is wanted.

        Return: (link, modal)
        """
        # https://v4-alpha.getbootstrap.com/components/modal/#examples
        #text = escape(self.to_bibtex(), tag="pre")
        text = highlight(self.to_bibtex(), BibTeXLexer(), HtmlFormatter(cssclass="codehilite"))
        # Construct ids from self.key as they are unique.
        modal_id, modal_label_id = "modal-id-%s" % self.key, "modal-label-id-%s" % self.key

        if link:
            btn = """<a data-toggle="modal" href="#{modal_id}">bibtex</a>""".format(**locals())
        else:
            btn = """\
<button type="button" class="btn btn-primary btn-xsm btn-labeled small-text" data-toggle="modal" data-target="#{modal_id}">
  <span class="btn-label"><i class="fa fa-id-card" aria-hidden="true"></i></span>bibtex
</button>""".format(**locals())

        modal = """\
<div class="modal fade" id="{modal_id}" tabindex="-1" role="dialog" aria-labelledby="{modal_label_id}">
  <div class="modal-dialog modal-lg" role="document">
    <div class="modal-content">
      <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
        <h4 class="modal-title" id="{modal_label_id}">bibtex</h4>
      </div>
      <div class="modal-body">{text}</div>
    </div>
  </div>
</div>""".format(**locals())

        return btn, modal


_WEBSITE = None

class Website(object):
    """
    This object is a singleton. It stores all the information required to generated the HTML documentation
    (input variables, test suite, bibtex entries).
    It also provides methods such as `get_wikilink` that will be invoked by the python markdown parser
    to implement extensions to the standard markdown syntax.
    """
    # Regular expression for wikilinks.
    #WIKILINK_RE = r'\[\[([\w0-9_ -]+)\]\]'
    #WIKILINK_RE = r'\[\[([\w0-9_ -\./]+)\]\]'
    WIKILINK_RE = r'\[\[([^\[]+)\]\]'
    #WIKILINK_RE = r'(?![~`])\[\[([^\[]+)\]\]'

    @classmethod
    def build(cls, root, deploy, verbose):
        """
        Build Website object from directory ``root`` and cache it.
        Main entry point for client code.
        """
        global _WEBSITE
        if _WEBSITE is not None:
            raise RuntimeError("website has been already costructed")
        _WEBSITE = cls(root, deploy, verbose=verbose)
        return _WEBSITE

    @classmethod
    def get(cls):
        """Return Website instance. Assume object already initialized with build_website."""
        global _WEBSITE
        if _WEBSITE is None:
            raise RuntimeError("website must be constructuted by calling `Website.build`")
        return _WEBSITE

    def __init__(self, root, deploy, verbose=0):
        start = time.time()
        self.root = os.path.abspath(root)
        self.deploy = bool(deploy)
        self.verbose = verbose
        self.md_generated = []
        self.ignored_paths = []
        self.warnings = []

        # Read mkdocs configuration file.
        # TODO: Should read Abinit version from a centralized file.
        with io.open(os.path.join(self.root, "..", "mkdocs.yml"), "rt", encoding="utf-8") as fh:
            if hasattr(yaml, "FullLoader"):
                self.mkdocs_config = yaml.load(fh, Loader=yaml.FullLoader)
            else:
                self.mkdocs_config = yaml.load(fh)

        # Build parser to convert Markdown to HTML.
        # The parser must support the same extensions as those used by mkdocs
        # so we initialize it from the options specified in mkdocs.yml
        # This implies that all extensions requirining the website singlecto must post-pone the import.
        #   * extensions: A list of extensions, which can either
        #       be strings or objects.  See the docstring on Markdown.
        #   * configs: A dictionary mapping module names to config options

        extensions, extension_configs = [], {}
        for item in self.mkdocs_config["markdown_extensions"]:
            #print(item, type(item))
            if isinstance(item, dict):
                assert len(item) == 1 and len(item.values()) == 1
                modname = str(list(item.keys())[0])
                extensions.append(modname)
                v = list(item.values())[0]
                #print(v)
                if v is not None:
                    extension_configs[modname] =  v
            else:
                extensions.append(str(item))
        self.markdown = markdown.Markdown(extensions=extensions, extension_configs=extension_configs)

        # Build database with all input variables indexed by code name.
        from .variables import get_codevars
        self.codevars = get_codevars()

        # Get bibtex references and cast to MyEntry instance.
        bib_path = os.path.join(self.root, "abiref.bib")
        self.bib_data = parse_file(bib_path, bib_format="bibtex")
        for entry in self.bib_data.entries.values():
            entry.__class__ = MyEntry

        # Read code statistics from file and produce json file.
        self.abinit_stats = AbinitStats(os.path.join(self.root, "maintainers", "statistics.txt"))
        self.abinit_stats.json_dump(os.path.join(self.root, "developers", "statistics.json"))

        # Build flat list of tests.
        from doc import tests as tmod
        self.abinit_tests = tmod.abitests
        tests = tmod.abitests.select_tests(suite_args=[], regenerate=True, flat_list=True)

        # Construct dictionary rpath --> test. Use OrderedDict to have deterministic behaviour.
        self.rpath2test = {}
        for t in tests:
            key = os.path.join(*splitall(t.inp_fname)[-4:])
            self.rpath2test[key] = t
        self.rpath2test = OrderedDict([(k, self.rpath2test[k]) for k in sorted(self.rpath2test.keys())])
        #print(self.rpath2test.keys())

        # Find variables used in tests.
        for var in self.codevars.iter_allvars():
            assert not hasattr(var, "tests")
            var.tests = []
            var.tests_info = {}

        def test_get_varnames(test, varnames):
            # TODO: This should become a method of BaseTest and must be improved.
            # See new method of BaseTest...
            with io.open(test.inp_fname, "rt", encoding="utf-8") as fh:
                s = fh.read()
            vused = [v for v in varnames if v in s]
            return vused

        codes_without_vars = set()
        white_list = set([
            "atompaw", "cut3d",
            #"multibinit",
            "mrgddb", "mrggkk", "band2eps",
            #"ujdet",
            "fold2Bloch", "fftprof",
            #"conducti",
            "mrgscr",
            #"macroave",
            "mrgdv", "a-TDEP"
        ])

        for test in tests:
            vd = self.codevars.get(test.executable, None)
            # Not all codes have variables documented in the database e.g. multibinit
            if vd is None:
                if test.executable not in codes_without_vars:
                    codes_without_vars.add(test.executable)
                    if test.executable not in white_list:
                        cprint("WARNING: Cannot find variables associated to code: `%s`" % test.executable, "yellow")
                continue
            for vname in test_get_varnames(test, list(vd.keys())):
                var = vd[vname]
                var.tests.append(test)

        # Pre-compute vars.tests and their frequency.
        for var in self.codevars.iter_allvars():
            var.tests_info["num_all_tests"] = len([t for t in tests if t.executable == var.executable])
            var.tests_info["num_all_tutorial_tests"] = len([t for t in tests
                if t.executable == var.executable and t.suite_name.startswith("tuto")])
            var.tests_info["num_tests_in_tutorial"] = len([t for t in var.tests
                if t.executable == var.executable and t.suite_name.startswith("tuto")])

        # Find pdf files and sort them by basename.
        self.pdfs = OrderedDict(sorted([t for t in self.walk_filepath() if t[0].endswith(".pdf")],
                                key=lambda t: t[0]))

        cprint("Initial website generation completed in %.2f [s]" % (time.time() - start), "green")

        self.do_not_edit_comment = """\
<!--
This file is automatically generated by mksite.py. All changes will be lost.
Change the input yaml files or the python code
-->
"""

    def walk_filepath(self):
        """
        Iterate over the files stored in the doc directory. Return (filename, path).
        Files in site and ~abinit/doc/tests are excluded.
        """
        excludes = [os.path.join(self.root, f) for f in ("site", os.path.join("doc", "tests"))]
        for root, dirs, files in os.walk(self.root, topdown=True):
            if any(root.startswith(e) for e in excludes):
                print("Excluding root", root)
                dirs[:] = []
                continue
            #print(root)
            for f in files:
                if f.startswith("_"): continue
                #if f == "README.md": continue
                yield f, os.path.join(root, f)

    def warn(self, msg):
        """Print warning message to terminal and save it for future reference."""
        msg = "WARNING: %s" % msg
        self.warnings.append(msg)
        cprint(msg, color="yellow")

    def convert_markdown(self, source):
        """"
        Convert markdown string `source` to serialized HTML.
        """
        self.markdown.reset()
        return my_unicode(self.markdown.convert(source))

    def new_mdfile(self, dirname, mdname, meta=None, with_comment=True):
        """
        Create new markdown file with name `mdname` in directory `dirname`.
        `meta` is an optional dictionary with meta-variables added to the front matter.

        Return: File object.

        .. warning::

            Unicode characters in meta are not supported (annoying portability issue with py2.7)
        """
        dirpath = os.path.join(self.root, dirname)
        if not os.path.isdir(dirpath): os.mkdir(dirpath)
        path = os.path.join(dirpath, mdname)
        assert path not in self.md_generated
        self.md_generated.append(path)
        if self.verbose: print("Generating markdown file: `%s`" % path)

        mdf = io.open(path, "wt", encoding="utf-8")
        if meta is not None:
            # Must convert to ASCII to avoid !!python/unicode tags in YAML doc
            # (mkdocs does not use yaml to parse the front matter).
            if sys.version_info[0] <= 2:
                meta = {k.encode("ascii", errors="strict"): meta[k] for k in meta}
            s = yaml.dump(meta, indent=4, default_flow_style=False).strip().replace(" !!python/unicode", "")
            mdf.write("---\n%s\n---\n" % s)

        if with_comment:
            mdf.write(self.do_not_edit_comment)
        mdf.rpath = "/" + os.path.relpath(path, self.root)

        return mdf

    def generate_mdindex(self, dirname):
        """
        Generate the index.md file from the meta data section given in
        the md files stored in directory `dirname`.
        """
        workdir = os.path.join(self.root, dirname)
        pages = [MarkdownPage(os.path.join(workdir, fname), self) for fname in os.listdir(workdir)
                if fname.endswith(".md") and fname != "index.md"]

        #index_md = ["# Index of files in %s.\n\n" % dirname]
        for page in sorted(pages, key=lambda p: p.basename):
            try:
                desc = page.meta["description"]
            except KeyError:
                raise KeyError(
                    "Markdown page `%s` does not have `description` key in front matter.\n"
                    "This is required to generate index.md automatically in python." % page.path)
            desc = desc[0]
            #index_md.append("*  [%s](%s): %s" % (page.basename.replace(".md", ""), page.url, desc))

        #with self.new_mdfile(dirname, "index.md") as mdf:
        #    mdf.write("\n".join(index_md))

    def copy_readme_files(self):
        """
        Copy README_*.md files from ~abint to ~abinit/doc and *git ignore* them.
        Files must be included in mkdocs.yml in the `Installation` section.
        """
        top = os.path.abspath(os.path.join(self.root, ".."))
        for f in os.listdir(top):
            if f.startswith("README_") and f.endswith(".md"):
                src = os.path.join(top, f)
                dest = os.path.join(self.root, f)
                shutil.copy(src, dest)
                self.ignored_paths.append(dest)

    def generate_page_with_ac_examples(self):
        """Generate markdown pages with all ac exaples found in config-examples."""
        dirpath = os.path.join(self.root, "build", "config-examples")
        md_lines = []
        app = md_lines.append
        app("""
# Autoconf examples

This page gathers the autoconf files used by the buildbot testfarm

""")
        for f in os.listdir(dirpath):
            path = os.path.join(dirpath, f)
            if os.path.isdir(path) or path.endswith(".swp"): continue
            app("## %s  " %  f)
            with io.open(path, "rt", encoding="utf-8") as fh:
                # Remove all comments except for options that are specified.
                #print(path)
                ac_lines = []
                inblock = False
                for l in reversed(fh.readlines()):
                    l = l.strip()
                    if not l:
                        inblock = False
                        continue
                    if l and l[0].isalpha():
                        inblock = True
                        ac_lines.append(l + "\n")
                    elif inblock:
                        ac_lines.append(l)

                app("\n\n```shell")
                md_lines.extend(reversed(ac_lines))
                app("```\n")

        # Write MD file.
        with self.new_mdfile("developers", "autoconf_examples.md") as mdf:
            mdf.write("\n".join(md_lines))

    def generate_markdown_files(self):
        """Generate markdown files using the data stored in the bibtex file, the abivars file ..."""
        start = time.time()

        self.copy_readme_files()
        self.generate_page_with_ac_examples()

        # Write index.md with the description of the input variables.
        meta = {"description": "Complete list of Abinit input variables"}
        with self.new_mdfile("variables", "index.md", meta=meta) as mdf:
            mdf.write("\n\n# Input variables \n\n")
            for code, vd in self.codevars.items():
                mdf.write("## %s variables   \n\n" % code)
                mdf.write(vd.get_vartabs_html(self, mdf.rpath))

            # This for the table of variables implemented by Jordan
            mdf.write(self.build_varsearch_html(mdf.rpath))

        # Build markdown page with external parameters.
        with self.new_mdfile("variables", "external_parameters.md") as mdf:
            mdf.write("""\
This document lists and provides the description of the name (keywords) of external parameters
that are not input variables, but that are used in the documentation of other variables,
typically compilation parameters, available libraries, or number of processors.
You can change these parameters at compile or run time usually.

""")
            for pname, info in self.codevars.all_external_params.items():
                mdf.write("## %s  \n%s  \n\n" % (pname, info))

        # Build markdown pages for the different sets of variables.
        for code, vd in self.codevars.items():
            cprint("Generating markdown files with input variables of code: `%s`..." % vd.executable, "green")
            for varset in vd.my_varset_list:
                var_list = [v for v in vd.values() if v.varset == varset]
                meta = {"description": "%s input variables" % varset}
                with self.new_mdfile("variables", varset + ".md", meta=meta) as mdf:
                    mdf.write("""\
# {varset} input variables

This document lists and provides the description of the name (keywords) of the
{varset} input variables to be used in the input file for the {executable} executable.

""".format(varset=varset, executable=vd.executable))

                    for i, var in enumerate(var_list):
                        mdf.write(var.to_abimarkdown(with_hr=False))

        # Add plotly figures.
        # TODO: Replace it with dot
        if False and self.deploy:
            with self.new_mdfile("variables", "connections.md", meta={"plotly": True}) as mdf:
                mdf.write("# Dependency graphs  \n")
                mdf.write("""
These graphs show the dependencies of the input variables towards each other.
The colormap gives the number of input variables connected to the node.

""")
                for code, vd in self.codevars.items():
                    for varset in vd.my_varset_list:
                        mdf.write("## %s, varset: %s  \n\n" % (code, varset))

        # Write Markdown page with statistics.
        with self.new_mdfile("variables", "varset_stats.md") as mdf:
            mdf.write("""
# Input variables, statistics

This document lists the input variables for ABINIT and three post-processors of ABINIT,
in order of number of occurrence in the input files provided with the package.

""")
            for code, vd in self.codevars.items():
                num_tests = len([test for test in self.rpath2test.values() if test.executable == code])
                mdf.write("\n\n## %s \n\n" % code)
                mdf.write("%d tests\n\n" % num_tests)
                # TODO The number of tests is smaller than ecut! Count Tutorial
                items = sorted([(len(v.tests), v) for v in vd.values()], key=lambda t: t[0], reverse=True)
                # https://www.w3schools.com/bootstrap/bootstrap_list_groups.asp
                lines = ['<ul class="list-group">']
                for count, group in groupby(items, key=lambda t: t[0]):
                    vlist = [item[1] for item in sorted(group, key=lambda t: t[1].name)]
                    rpath = os.path.join(mdf.rpath.replace(".md", ""), "index.md")
                    s = ", ".join(v.internal_link(self, rpath) for v in vlist)
                    # Set color depending on coverage.
                    ratio = 100 * count / num_tests
                    if ratio > 40:
                        cls = "list-group-item-success"
                    elif ratio > 2:
                        cls = "list-group-item-warning"
                    else:
                        cls = "list-group-item-danger"
                    lines.append('<li class="list-group-item %s"> %s <span class="badge"> %d </span></li>' % (cls, s, count))
                mdf.write("\n".join(lines) + "</ul>")

        # Topics
        cprint("Generating Markdown files with topics ...", "green")
        self.all_topics = ABI_TOPICS
        self.all_relevances = ABI_RELEVANCES

        dirpath = os.path.join(self.root, "topics")
        all_mdfiles = [f for f in os.listdir(dirpath) if f.endswith(".md") and f.startswith("_")]
        for topic in self.all_topics:
            mdname = "_" + topic + ".md"
            try:
                all_mdfiles.remove(mdname)
            except ValueError:
                cprint("Cannot find `%s` in all_mdfiles" % mdname, "yellow")

        if all_mdfiles:
            raise RuntimeError("Found md files in topics not listed in `variables_code.py` modules\n%s" % (
                str(all_mdfiles)))

        # datastructures needed for topics index.md
        index_md = ["# Alphabetical list of topics\n"]
        self.howto_topic = {}
        for topic in self.all_topics:
            # Read description from md file.
            with io.open(os.path.join(dirpath, "_" + topic + ".md"), "rt", encoding="utf-8") as fh:
                for line in fh:
                    if "description:" in line:
                        self.howto_topic[topic] = line.replace("description:", "").strip()
                        break
                else:
                    raise RuntimeError("Cannot find `description:` in topic file: `%s`" % topic)

            # Find list of variables associated to this topic
            # Order and group vlist by relevances and write list with links.
            # TODO: Can we have multiple relevances with the same topic?
            related_variables = "No variable associated to this topic."
            vlist = [var for var in self.codevars.iter_allvars() if topic in var.topic2relevances]
            if vlist:
                lines = []
                def sort_relevances(t):
                    # TODO: Add rank to ABI_RELEVANCES
                    try:
                        return {"compulsory": 0, "basic": 1, "useful": 2, "expert": 3, "internal": 4,
                                "prpot": 5, "prfermi": 6, "prden": 7, "prgeo": 8, "prdos": 9, "prgs": 10,
                                "prngs": 11, "prmisc": 12}[t[0]]
                    except KeyError:
                        raise KeyError("Cannot find relevance `%s` in dict. Add it to sort_relevances with the proper rank."
                                % str(t))

                # Build list of (relevance, variable) tuple then sort and group by relevance.
                items = [(v.topic2relevances[topic][0], v) for v in vlist]
                for num, group in sort_and_groupby(items, key=lambda t: sort_relevances(t)):
                    # Alphabetical order inside group.
                    group = list(sorted(group, key=lambda t: t[1].name))
                    relevance = group[0][0]
                    lines.append("*%s:*\n" % relevance)
                    lines.extend("- %s  %s" % (v.wikilink, v.mnemonics) for (_, v) in group)
                    lines.append(" ")
                related_variables = "\n".join(lines)

            # Find tests associated to this `topic`
            # Group tests by `suite_name` and write markdown list with links.
            items = [(rpath, test) for (rpath, test) in self.rpath2test.items() if topic in test.topics]
            selected_input_files = "No input file associated to this topic."
            if items:
                lines = []
                for suite_name, group in sort_and_groupby(items, key=lambda t: t[1].suite_name):
                    lines.append("*%s:*\n" % suite_name)
                    lines.extend("- [[%s]]" % rpath for (rpath, test) in group)
                    lines.append(" ")
                selected_input_files = "\n".join(lines)

            # Read template, interpolate and write md file included in mkdocs.yml.
            with io.open(os.path.join(self.root, "topics", "_" + topic + ".md"), "rt", encoding="utf-8") as fh:
                template = fh.read()
                template = template.replace("is the source file for this topics. Can be edited."," file has been generated automatically from the corresponding _* source file. DO NOT EDIT. Edit the source file instead.")
                template = template.replace("{{ related_variables }}", related_variables)
                template = template.replace("{{ selected_input_files }}", selected_input_files)

            with self.new_mdfile("topics", topic + ".md", with_comment=False) as mdf:
                mdf.write(template)

        # Now write topics index.md (sorted by first character)
        for firstchar, group in sort_and_groupby(self.all_topics, key=lambda t: t[0].upper()):
            index_md.append("## %s" % firstchar)
            index_md.extend("- [[topic:%s|%s]]: %s" % (topic, topic, self.howto_topic[topic]) for topic in group)

        meta = {"description": "List of Abinit topics"}
        with self.new_mdfile("topics", "index.md", meta=meta) as mdf:
            mdf.write("\n".join(index_md))

        # Build page with full list of tests grouped by `suite_name`.
        cprint("Generating Markdown file with tests ...", "green")
        meta = {"description": "List of Abinit tests"}
        items = [(rpath, test) for (rpath, test) in self.rpath2test.items()]

        with self.new_mdfile("developers", "testsuite.md", meta=meta) as mdf:
            for suite_name, group in sort_and_groupby(items, key=lambda t: t[1].suite_name):
                group = list(group)
                mdf.write('## %s  \n\n' % suite_name)
                for i, (rpath, test) in enumerate(group):
                    mdf.write('### [[%s]]   \n\n' % rpath)
                    mdf.write(my_unicode(test.description))
                    mdf.write("\n\n")
                    mdf.write("Executable: %s   \n" % test.executable)
                    if test.keywords:
                        mdf.write("Keywords(s): %s   \n" % ", ".join(k for k in sorted(test.keywords)))
                    if test.topics:
                        mdf.write("Topic(s): %s  \n" % ", ".join("[[topic:%s]]" % t for t in test.topics))
                    if test.authors and "Unknown" not in test.authors:
                        mdf.write("Author(s): %s  \n" % ", ".join(a for a in sorted(test.authors)))
                    if i != len(group) - 1:
                        mdf.write("\n\n* * *\n\n")

        # All markdown files have been generated. Now scan all md files to find all wikilinks,
        # in particular the bibliographic references needed to generate backlinks.
        self.analyze_pages()

        # Now generate page with bibliography.
        cprint("Generating Markdown file with bibliographic entries ...", "green")
        citation2pages = defaultdict(list)
        for page in self.md_pages:
            for citation in page.citations:
                citation2pages[citation].append(page)

        meta = {"description": "Bibliographical references mentioned in the ABINIT documentation"}
        with self.new_mdfile("theory", "bibliography.md", meta=meta) as mdf:
            lines = []
            lines.append("""\
# Bibliography

This document lists all the bibliographical references mentioned in the ABINIT documentation,
with link(s) to the Web pages where such references are mentioned, as well as to the bibtex formatted reference.
The bibtex file is available [here](../abiref.bib).

""")
            for name in sorted(self.bib_data.entries.keys()):
                entry = self.bib_data.entries[name]
                lines.append("\n\n## **%s** \n\n" % name)
                try:
                    lines.append(entry.to_abimarkdown())
                except Exception as exc:
                    raise ValueError("Exception while trying to convert bibtex entry `%s`\n%s\n" % (name, str(exc)))
                if citation2pages[name]:
                    lines.append("Referred to in: %s" % ", ".join('[{url}]({url})'.format(url=url)
                        for url in sorted([page.url for page in citation2pages[name]])))

            mdf.write("\n".join(lines))

        meta = {"description": "List of PDF files provided by the Abinit documentation"}
        with self.new_mdfile("theory", "documents.md", meta=meta) as mdf:
            mdf.write("# PDF files  \n")
            for fname, path in self.pdfs.items():
                mdf.write("## %s  \n" % fname)
                rpdf = "/" + os.path.relpath(path, self.root)
                src = os.path.relpath(rpdf, mdf.rpath)
                html = '<embed src="{src}" type="application/pdf" width="100%" height="480px">\n\n'.format(src=src)
                mdf.write(html)

        #for dirname in ["theory"]:
        #    self.generate_mdindex(dirname)
        #topic2pages = defaultdict(list)
        #for page in self.md_pages:
        #    for topic in page.topics:
        #        topic2pages[topic].append(page)

        with open(os.path.join(self.root, ".gitignore"), "wt") as fh:
            fh.write("# The following md files have been copied from ~abinit and should be `git ignored`\n")
            for p in self.ignored_paths:
                fh.write(os.path.relpath(p, self.root) + "\n")

            fh.write("# The following md files have been automatically generated and should be `git ignored`\n")
            for p in self.md_generated:
                fh.write(os.path.relpath(p, self.root) + "\n")

        cprint("Markdown files generation completed in %.2f [s]" % (time.time() - start), "green")

    def analyze_pages(self):
        """
        Analyze all markdown pages, find wiklinks in pages required to generate backlinks in docs.
        """
        cprint("Analyzing markdown pages ...", "green")
        start = time.time()

        #ignored = set(["doc/developers/markdown.md"])

        self.md_pages, self.html_pages = [], []
        for f, path in self.walk_filepath():
            if f.startswith("_"): continue
            #if os.path.relpath(path, self.root) in ignored: continue
            #if f == "README.md": continue
            if f.endswith(".md"):
                self.md_pages.append(MarkdownPage(path, self))
            elif f.endswith(".html"):
                self.html_pages.append(HtmlPage(path, self))

        self.find_unreferenced_mds()
        cprint("Completed in %.2f [s]" % (time.time() - start), "green")

    def find_unreferenced_mds(self):
        """
        Extract all md pages listed in mkdocs.yml and compare them with the md files
        in docs directory. Issue a warning if the two sets are not equal.
        """
        def find_mds(obj):
            """Return list of md files reported in mkdocs.yml"""
            md_files = []
            if isinstance(obj, list):
                for item in obj:
                    md_files.extend(find_mds(item))
            elif isinstance(obj, dict):
                for key, value in obj.items():
                    md_files.extend(find_mds(value))
            elif hasattr(obj, "endswith"):
                # Assume string
                assert obj.endswith(".md")
                md_files.append(obj)
            else:
                raise TypeError("Don't know how to handle type %s\n%s" % (type(obj), str(obj)))

            return md_files

        pages_in_toolbar = []
        for entry in self.mkdocs_config["pages"]:
            pages_in_toolbar.extend(find_mds(entry))
        #for p in pages_in_toolbar: print(p)

        # Find elements in `pages_on_disk` not in `pages_in_toolbar`
        pages_in_toolbar = set(pages_in_toolbar)
        pages_on_disk = set(p.relpath for p in self.md_pages)
        diff = pages_on_disk.difference(pages_in_toolbar)
        if diff:
            self.warn("Found markdown files on disk not included in mkdocs.yml:\n%s" % "\n".join(diff))
        diff = pages_in_toolbar.difference(pages_on_disk)
        if diff:
            self.warn("Found markdown files in mkdocsyml not present in directories:\n%s" % "\n".join(diff))

    def slugify(self, value):
        """
        Slugify a string, to make it URL friendly. Use same convention as TOC extensions of python markdown.
        """
        from markdown.extensions.toc import slugify
        return slugify(value, separator="-")

    def preprocess_mdlines(self, lines):
        """Preprocess markdown lines."""
        lines = self._preprocess_aliases(lines)
        lines = self._preprocess_include(lines)
        lines = self._preprocess_macros(lines)
        return lines

    def _preprocess_macros(self, lines):
        """Preprocess markdown lines and replace [TUTORIAL_README] string."""

        tutorial_readme = """

!!! note

    Supposing you made your own install of ABINIT, the input files to run the examples
    are in the *~abinit/tests/* directory where *~abinit* is the absolute path of the abinit top-level directory.
    If you have NOT made your own install, ask your system administrator where to find the package, especially the executable and test files.

    To execute the tutorials, create a working directory (`Work*`) and
    copy there the input files and the *files* file of the lesson. This will be explicitly mentioned in the first lessons,
    that will tell you more about the *files* file (see also [[help:abinit#intro|section 1.1]]).
    The *files* file ending with *_x* (e.g. *tbase1_x.files*) **must be edited** every time you start to use a new input file.

    Most of the tutorials do not rely on parallelism (except specific [[tutorial:basepar|tutorials on parallelism]]).
    However you can run most of the tutorial examples in parallel, see the [[topic:parallelism|topic on parallelism]].

    In case you work on your own PC or workstation, to make things easier, we suggest you define some handy environment variables by
    executing the following lines in the terminal:

    ```bash
    export ABI_HOME=Replace_with_the_absolute_path_to_the_abinit_top_level_dir
    export PATH=$ABI_HOME/src/98_main/:$PATH
    export ABI_TESTS=$ABI_HOME/tests/
    export ABI_PSPDIR=$ABI_TESTS/Psps_for_tests/  # Pseudopotentials used in examples.
    ```

    Examples in this tutorial use these shell variables: copy and paste
    the code snippets into the terminal (**remember to set ABI_HOME first!**).
    The 'export PATH' line adds the directory containing the executables to your [PATH](http://www.linfo.org/path_env_var.html)
    so that you can invoke the code by simply typing *abinit* in the terminal instead of providing the absolute path.

"""

        new_lines = []
        for line in lines:
            if "[TUTORIAL_README]" in line:
                new_lines.extend(tutorial_readme.splitlines())
            else:
                new_lines.append(line)

        return new_lines

    def _preprocess_aliases(self, lines):
        """
        Handle aliases.
        |token| will be replaced by value by the Markdown preprocessor
        NB: white spaces in token are not allowed, `token` is ignored
        """
        def repl(matchobj):
            key = matchobj.group("key")
            if self.verbose: print("Found possible alias:", key)
            if key == "today":
                return datetime.date.today().strftime("%B %d, %Y")

            value = self.mkdocs_config["extra"]["abimkdocs_aliases"].get(key)
            if value is not None:
                if self.verbose: print("Returning", value)
                return " " + value + " "
            else:
                if self.verbose: print("Returning full match:", matchobj.group(0))
                return matchobj.group(0)

        alias_syntax = re.compile(r"[^`\$]\|(?P<key>\w+)\|")
        #alias_syntax = re.compile(r"(?!`+)\|(?P<key>\w+)\|")
        return [re.sub(alias_syntax, repl, line) for line in lines]

    def _preprocess_include(self, lines):
        """Handle {action ...} syntax."""
        inc_syntax = re.compile(r'^\{%\s*(.+?)\s*%\}')
        new_lines = []
        for line in lines:
            m = inc_syntax.search(line)
            if not m:
                new_lines.append(line)
            else:
                args = m.group(1).split()
                action = args.pop(0)
                if self.verbose: print("Triggering action:", action, "with args:", str(args))

                # Dispatch according to action.
                if action == "modal":
                    if len(args) > 1:
                        new_lines.extend(self.modal_with_tabs(args).splitlines())
                    else:
                        new_lines.extend(self.modal_from_filename(args[0]).splitlines())
                elif action == "dialog":
                    if len(args) > 1:
                        new_lines.extend(self.dialogs_from_filenames(args).splitlines())
                    else:
                        new_lines.extend(self.dialog_from_filename(args[0]).splitlines())
                elif action == "include":
                    with io.open(args[0], "rt", encoding="utf-8") as f:
                        new_lines.extend([l.rstrip() for l in f])
                else:
                    raise ValueError("Don't know how to handle action: `%s` in token: `%s`" % (action, m.group(1)))

        return new_lines

    @staticmethod
    def parse_wikilink_token(token):
        """
        Parse wikilink token of the form `namespace:name#fragment|text||args`
        where namespace, fragment and text are optional

        Return: (namespace, name, fragment, text)
            Individual entries are set to None if non present in token.
        """
        #args = ""
        #if "||" in token:
        #    token, args = token.split("||")

        text = None
        if "|" in token:
            token, text = token.split("|")
            text = text.strip()

        fragment = None
        if "#" in token:
            token, fragment = token.split("#")
            fragment = fragment.strip()

        namespace = None
        if ":" in token:
            namespace, name = token.split(":")
            namespace, name = namespace.strip(), name.strip()
        else:
            name = token.strip()
            if not name: name = None

        return namespace, name, fragment, text

    def get_wikilink(self, token, page_rpath):
        """
        Invoked by the wikilink extension to implement the wikilink syntax: [namespace:name#fragment|text]

        Args:
            token: The string enclosed between square brackets.
            page_rpath: The root-relative path of the markdown file (needed to generate relative links).

        Return:
            :class:`etree.Element` instance representing the HTML anchor. classes are automatically
                addeded to the link so that we can style them with CSS.
        """
        token = token.strip()
        if not token:
            self.warn("Empty wikilink in %s" % page_rpath)
            return ""

        #if token.startswith("~~") and token.endswith("~~"):
        #    token = token[2:-2]
        #    try:
        #        a = self.get_wikilink(token, page_rpath)
        #        return a.text
        #    except:
        #        return token

        html_classes = ["wikilink"]
        target = ""
        a = etree.Element("a")

        if any(token.startswith(prefix) for prefix in ("www.", "http:", "https:", "ftp:", "file:")):
            # Handle [[www.google.com|text]]
            url, a.text = token, token
            if "|" in token: url, a.text = token.split("|")
            a.set('href', url)
            a.set('target', "_blank")
            return a

        # [[namespace:name#fragment|text]]
        try:
            namespace, name, fragment, a.text = self.parse_wikilink_token(token)
        except ValueError:
            raise ValueError("Cannot parse wikilink token `%s`" % token)

        if namespace is not None and name is None:
            raise ValueError("Wrong wikilink token: `%s` in `%s`.\nnamespace is not None and name is None" %
                    (token, page_rpath))

        # Treat different cases and define `url` and `text`
        # Note that url is a root-relative URL that will be converted to relative URL at the end.
        if namespace is None:
            if name is None:
                # Handle [[#internal_link|text]]
                assert fragment is not None
                url = ""
                if a.text is None: a.text = fragment
            else:

                if "@" in name:
                    # Handle [[dipdip@anaddb|text]]
                    vname, code = name.split("@")
                    var = self.codevars[code][vname]
                    url = "/variables/%s#%s" % (var.varset, var.name)
                    if a.text is None: a.text = name
                    html_classes.append("codevar-wikilink")

                elif name in self.codevars["abinit"]:
                    # Handle link to Abinit variable e.g. [[ecut|text]]
                    var = self.codevars["abinit"][name]
                    url = "/variables/%s#%s" % (var.varset, var.name)
                    html_classes.append("codevar-wikilink")
                    if a.text is None:
                        a.text = var.name if not var.is_internal else "%%%s" % var.name

                elif name.startswith("tests/") or name.startswith("~abinit/tests/"):
                    assert fragment is None
                    if a.text is None: a.text = name
                    if "Psps_for_tests" in name:
                        # Handle [[~abinit/tests/Psps_for_tests/6c.lda.atompaw]]
                        nm = name.replace("~abinit/", "")
                        url = "/" + nm
                    else:
                        # Handle [[tests/tutorial/Refs/tbase1_2.out|text]]
                        #if not text.startswith("~abinit/"): text = "~abinit/" + text
                        nm = name.replace("~abinit/", "")
                        url = "/" + nm

                        # Add popover with test description if input file.
                        if nm in self.rpath2test:
                            test = self.rpath2test[nm]
                            content = test.description # + "\n\n" + ", ".join(test.authors)
                            add_popover(a, content=content)

                    target = "_blank"
                    html_classes.append("abifile-wikilink")

                elif name in self.codevars.all_characteristics:
                    # handle [[ENERGY]] by building internal link to abinit user guide
                    url = "/guide/abinit#parameters"
                    if a.text is None: a.text = name

                elif name in self.codevars.all_external_params:
                    # handle [[AUTO_FROM_PSP]] by building link with popover
                    content = ("This is an external parameter\n"
                               "typically compilation parameters, available libraries, or number of processors.\n"
                               "You can change these parameters at compile or runtime usually.\n")
                    url = "/variables/external_parameters#%s" % self.slugify(name)
                    if a.text is None: a.text = name
                    add_popover(a, title=self.codevars.all_external_params[name], content=content)

                else:
                    self.warn("Don't know how to handle wikilink token `%s` in `%s`" % (token, page_rpath))
                    url, a.text = "FAKE_URL", "FAKE_URL"

        else:
            # namespace is defined
            if namespace in self.codevars:
                # Handle [[anaddb:asr|text]] or [[abinit:ecut|text]]
                assert fragment is None
                var = self.codevars[namespace][name]
                url = "/variables/%s#%s" % (var.varset, var.name)
                html_classes.append("codevar-wikilink")
                if a.text is None:
                    a.text = var.name if not var.is_internal else "%%%s" % var.name

            elif namespace == "lesson" or namespace == "tutorial" :
                # Handle [[tutorial:wannier90|text]]
                if name == "index":
                    url = "/tutorial/"
                    if a.text is None: a.text = "tutorial home page"
                else:
                    url = "/tutorial/%s" % name
                    if a.text is None: a.text = "%s %s" % (name, namespace)
                html_classes.append("lesson-wikilink")

            elif namespace == "help" or namespace == "guide" :
                # Handle [[help:optic|text] NB: [[help:codename]] is echoed "codename help file"
                if name == "index":
                    url = "/guide/"
                    if a.text is None: a.text = "user-guide home page"
                else:
                    url = "/guide/%s" % name
                    if a.text is None: a.text = "%s help file" % name
                html_classes.append("user-guide-wikilink")

            elif namespace == "about" :
                # Handle [[about:release-notes|text] NB: [[about:file]] is echoed "file"
                if name == "index":
                    url = "/about/"
                    if a.text is None: a.text = "no index for about at present"
                else:
                    url = "/about/%s" % name
                    if a.text is None: a.text = "%s" % name
                html_classes.append("about-wikilink")

            elif namespace == "developers" :
                # Handle [[developers:psp8_info|text] NB: [[developers:filename]] is echoed "filename developer doc"
                if name == "index":
                    url = "/developers/"
                    if a.text is None: a.text = "no index for developers at present"
                else:
                    url = "/developers/%s" % name
                    if a.text is None: a.text = "%s developer doc" % name
                html_classes.append("developers-wikilink")

            elif namespace == "topic":
                # Handle [[topic:BSE|text]]
                html_classes.append("topic-wikilink")
                if name == "index":
                    url = "/topics/"
                    if a.text is None: a.text = "Topics index"
                else:
                    url = "/topics/%s" % name
                    if a.text is None: a.text = "%s_%s" % (namespace, name)
                    add_popover(a, content=self.howto_topic[name])

            elif namespace == "cite":
                # Handle [[bib:biblio|bibliography]]
                if name == "biblio":
                    url = "/theory/bibliography/"
                    if a.text is None: a.text = "bibliography"
                else:
                    # Handle [[bib:Amadon2008]]
                    try:
                        ref = self.bib_data.entries[name]
                        url = "/theory/bibliography#%s" % self.slugify(name)
                        content = ref.fields["title"].replace("{", "").replace("}", "") #+ "\n\n" + ref.authors
                        add_popover(a, content=content)
                        if a.text is None: a.text = "[%s]" % name
                        html_classes.append("citation-wikilink")
                    except Exception as exc:
                        self.warn("Exception `%s:%s`\nwhile treating wikilink token: `%s` in `%s`" %
                                (exc.__class__, str(exc), token, page_rpath))
                        url, a.text = "FAKE_URL", "FAKE_URL"

            elif namespace == "theory":
                # Handle [[theorydoc:mbpt|text]]
                url = "/theory/%s" % name
                html_classes.append("theory-wikilink")
                if a.text is None: a.text = name

            elif namespace == "varset":
                # Handle [[varset:BSE|text]]
                assert fragment is None
                if name == "allvars":
                    url = "/variables/"
                else:
                    url = "/variables/%s" % name
                if a.text is None: a.text = "%s varset" % name

            elif namespace == "test":
                # Handle [[test:libxc_41]] (syntax for suite) [[test:gspw_01]] (syntax for subsuite)
                tokens = name.split("_")
                prefix, tnum = "_".join(tokens[:-1]), tokens[-1]
                if prefix in self.abinit_tests.all_subsuite_names:
                    # [[test:gspw_01]]  --> Need to get the name of suite from subsuite.
                    suite_name = self.abinit_tests.suite_of_subsuite(prefix).name
                    url = "/tests/%s/Input/t%s.in" % (suite_name, name)
                else:
                    # [[test:libxc_41]]
                    url = "/tests/%s/Input/t%s.in" % (prefix, tnum)

                if a.text is None: a.text = "%s[%s]" % (prefix, tnum)
                test = self.rpath2test[url[1:]]
                content = test.description # + "\n\n" + ", ".join(test.authors)
                add_popover(a, content=content)
                target = "_blank"
                html_classes.append("abifile-wikilink")

            elif namespace == "src":
                # Handle [[src:94_scfcv/scfcv.F90]]
                url = "https://github.com/abinit/abinit/blob/master/src/%s" % name
                if a.text is None: a.text = name
                target = "_blank"
                html_classes.append("abifile-wikilink")

            elif namespace == "ac":
                # Handle [[ac:abiref_gnu_9.2_debug.ac]]
                # The following is incorrect: files in /build/config-examples are generated when makemake is issued.
                # url = "/build/config-examples/%s" % name
                # By contrast, the following is a permanent reference
                # FIXME: buildsys refs are not generated anymore (YP)
                #url = "/abichecks/buildsys/Refs/%s" % name
                #if a.text is None: a.text = name
                #target = "_blank"
                #html_classes.append("abifile-wikilink")
                url = "/build/config-template.ac9"
                pass

            elif namespace == "pdf":
                # Handle [[pdf:howto_chebfi.pdf]] or [[pdf:howto_chebfi]]
                if not name.endswith(".pdf"): name += ".pdf"
                try:
                    path = self.pdfs[name]
                    url = "/" + os.path.relpath(path, self.root)
                except KeyError:
                    self.warn("Cannot find pdf file `%s` specified in wikilink `%s` in `%s`" % (name, token, page_rpath))
                    url, a.text = "FAKE_URL", "FAKE_URL"

                if a.text is None: a.text = name
                target = "_blank"
                html_classes.append("abifile-wikilink")

            elif namespace == "gitsha":
                # Handle [gitsha:f74dba1ed8346ca586dc95fd10fe4b8ced108d5e]
                url = "https://github.com/abinit/abinit/commit/%s" % name
                if a.text is None: a.text = name[:7]
                target = "_blank"
                html_classes.append("abigit-wikilink")

            # TODO? Issue
            #Fix issue https://github.com/abinit/abinit/issues/1
            else:
                self.warn("Don't know how to handle wikilink token `%s` in `%s`" % (token, page_rpath))
                url, a.text = "FAKE_URL", "FAKE_URL"

        a.set("class", " ".join(html_classes))
        if fragment is not None: url = "%s#%s" % (url, fragment)

        if sys.version_info[0] <= 2:
            from urlparse import urlparse
        else:
            from urllib.parse import urlparse

        o = urlparse(url)
        if o.scheme:
            a.set('href', url)
            return a

        # From root-relative url to relative url.
        end = ""
        if "#" in url:
            url, end = url.split("#")

        if not url:
            # Handle `#internal_link`
            url = "#" + end
        else:
            if not page_rpath.startswith("/"): page_rpath = "/" + page_rpath
            page_rpath = os.path.dirname(page_rpath.replace(".md", ""))
            url = os.path.relpath(url, page_rpath)
            if end: url = "%s#%s" % (url, end)

        if self.verbose: print("token", token, "page_rpath", page_rpath, "url", url)
        a.set('href', url.strip())
        if target: a.set('target', target)
        return a

    def build_varsearch_html(self, page_rpath):
        """Build single dictionary mapping varname --> var. Add @code if not abinit."""
        allvars = {}
        for code, vd in self.codevars.items():
            allvars.update({v.abivarname: v for v in vd.values()})

        tabs = "\n".join("""\
<a class="TabLetterLink" href="#{cap_char}" onClick="openLetter(event,'{cap_char}')" id="click{cap_char}">{cap_char}</a>""".format(cap_char=cap_char) for cap_char in sorted(set([k[0].upper() for k in allvars])))

        html_vars = ""
        for char, group in sort_and_groupby(list(allvars.items()), key=lambda t: t[0][0].upper()):
            lis = "\n".join("<li>{link}</li>".format(
                link=var.internal_link(self, page_rpath, label=var.abivarname, cls="small-grey-link")) for _, var in group)

        #for char, group in sort_and_groupby(allvars, key=lambda t: t[0][0].upper()):
        #    group = list(group)
        #    lis = []
        #    for i, (abivarname, var) in enumerate(group):
        #        if (i % 4) == 0 and i != 0: lis.append('</div>')
        #        if (i % 4) == 0 and i != len(group) - 1 : lis.append('<div class="row">')
        #        lis.append("""<li class="{col_cls}">{link}</li>""".format(
        #            col_cls="col-md-3",
        #            link=var.internal_link(self, page_rpath, label=abivarname, cls="")))
        #    if lis[-1] != '</div>': lis.append('</div>')
        #    lis = "\n".join(lis)

            html_vars += """
<li><ul id="{char}" class="TabContentLetter">
<li class="HeaderLetter">{char}</li> {lis} </ul></li>""".format(char=char, lis=lis)

        # NB: <form> is needed in order not to trigger the f/s keydown event registered by mkdocs-material.
        search_form = """
<div class="md-container">
  <div class="input-group custom-search-form">
    <form>
      <input type="text" class="form-control" id="InputSearch" onkeyup="searchInput()"
	onClick="searchInput()" placeholder="Search">
    </form>
    <span class="input-group-btn">
      <button class="btn btn-primary" type="submit" onClick="searchInput()">
        <span class="glyphicon glyphicon-search"></span>
      </button>
    </span>
  </div>
</div>

<script> $(function() {defaultClick(true);}); </script>
"""
        return """

## All variables

See aim, anaddb or optic for the subset of input variables for the executables AIM(Bader), ANADDB and OPTIC.
Such input variables are specifically labelled @aim, @anaddb, or @optic in the input variable database.
Enter any string to search in the database. Clicking without any request will give all variables.

{search_form}

<div class="TabsLetter">
{tabs}
</div>

<ul id="Letters">
{html_vars}
</ul>""".format(**locals())

    def dialogs_from_filenames(self, paths):
        buttons, dialogs = [], []
        for path in paths:
            btn, dialog = self.dialog_from_filename(path, ret_btn_dialog=True)
            buttons.append(btn)
            dialogs.append(dialog)

        button_group = '<div class="text-center"><div class="btn-group-vertical">\n%s\n</div></div>' % "\n".join(buttons)
        return button_group + "\n".join(dialogs)

    def dialog_from_filename(self, path, title=None, ret_btn_dialog=False):
        """Build customized jquery dialog to show the content of filepath `path`."""
        title = path if title is None else title
        with io.open(os.path.join(self.root, path), "rt", encoding="utf-8") as fh:
            if path.endswith(".in"):
                text = highlight(fh.read(), BashLexer(), HtmlFormatter(cssclass="codehilite small-text"))
            else:
                text = escape(fh.read(), tag="pre", cls="small-text")

        btn_id, dialog_id = gen_id(n=2)
        button = """\
<button type="button" id="{btn_id}" class="btn btn-default btn-labeled">
  <span class="btn-label"><i class="fa fa-window-restore" aria-hidden="true"></i></span>View {path}
</button>""".format(**locals())

        dialog = """
<div id="{dialog_id}" class="my_dialog" title="{title}" hidden><div>{text}</div></div>

<script> $(function() {{ abidocs_jqueryui_dialog("#{dialog_id}", "#{btn_id}") }}); </script>
""".format(**locals())

        if not ret_btn_dialog:
            button = '<div class="text-center">%s</div>' % button
            return button + dialog
        else:
            return button, dialog

    def modal_from_filename(self, path, title=None):
        """Return HTML string with bootstrap modal and content taken from file `path`."""
        # Based on https://v4-alpha.getbootstrap.com/components/modal/#examples
        # See also https://stackoverflow.com/questions/14971766/load-content-with-ajax-in-bootstrap-modal
        title = path if title is None else title
        with io.open(os.path.join(self.root, path), "rt", encoding="utf-8") as fh:
            text = escape(fh.read(), tag="pre", cls="small-text")

        return """\
<div class="text-center"> <!-- Button trigger modal -->
  <button type="button" class="btn btn-primary btn-labeled" data-toggle="modal" data-target="#{modal_id}">
    <span class="btn-label"><i class="glyphicon glyphicon-modal-window" aria-hidden="true"></i></span>View {path}
  </button>
</div>

<!-- Modal -->
<div class="modal fade" id="{modal_id}" tabindex="-1" role="dialog" aria-labelledby="{modal_label_id}">
  <div class="modal-dialog modal-lg" role="document">
    <div class="modal-content">
      <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
        <h4 class="modal-title" id="{modal_label_id}">{title}</h4>
      </div>
      <div class="modal-body">{text}</div>
    </div>
  </div>
</div>""".format(modal_id=gen_id(), modal_label_id=gen_id(), **locals())

    def modal_with_tabs(self, paths, title=None):
        # Based on http://jsfiddle.net/n__o/19rhfnqm/
        title = title if title else ""
        apaths = [os.path.join(self.root, p) for p in paths]
        button_label = "View " + ", ".join(paths)

        text_list = []
        for p in apaths:
            with io.open(p, "rt", encoding="utf-8") as fh:
                text_list.append(escape(fh.read(), tag="pre", cls="small-text"))
        tab_ids = gen_id(n=len(apaths))
        #print("paths", paths, "\ntab_ids", tab_ids)

        s = """\
<div class="text-center"> <!-- Button trigger modal -->
  <button type="button" class="btn btn-primary btn-labeled" data-toggle="modal" data-target="#{modal_id}">
    <span class="btn-label"><i class="glyphicon glyphicon-modal-window" aria-hidden="true"></i></span>{button_label}
  </button>
</div>

<!-- Modal -->
<div class="modal fade" id="{modal_id}" tabindex="-1" role="dialog" aria-labelledby="{modal_label_id}" aria-hidden="true">
  <div class="modal-dialog modal-lg" role="document">
    <div class="modal-content">
      <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
        <h4 class="modal-title" id="{modal_label_id}">{title}</h4>
      </div>
      <div class="modal-body">
        <div role="tabpanel">
          <!-- Nav tabs -->
          <ul class="nav nav-tabs" role="tablist">""".format(modal_id=gen_id(), modal_label_id=gen_id(), **locals())

        for i, (path, tid) in enumerate(zip(paths, tab_ids)):
            s += """\
          <li role="presentation" class="{li_class}">
          <a href="{href}" aria-controls="uploadTab" role="tab" data-toggle="tab">{path}</a>
          </li> """.format(li_class="active" if i == 0 else " ", href="#%s" % tid, path=path)

        s +=  """\
          </ul>
          <!-- Tab panes -->
          <div class="tab-content">"""

        for i, (text, tid) in enumerate(zip(text_list, tab_ids)):
            s += """<div role="tabpanel" class="tab-pane {active}" id="{tid}">{text}</div>""".format(
                    active="active" if i == 0 else " ", tid=tid, text=text)

        s += 6 * "</div>"

        return s


class Page(object):

    def __init__(self, path, website):
        self.path = os.path.abspath(path)
        self.website = website
        self.citations = set()
        self.topics = set()

    @property
    def basename(self):
        return os.path.basename(self.path)

    @property
    def relpath(self):
        return os.path.relpath(self.path, self.website.root)

    @property
    def relurl(self):
        return os.path.relpath(self.path, self.website.root).replace(".md", "")

    @property
    def url(self):
        return ("/" + self.relpath).replace(".md", "")


def add_popover(element, title=None, content=None, html=False):
    """
    Helper function to add popover to an anchor element.
    """
    # NB: Unfortunately, cannot subclass etree.Element in py2.7.
    def tos(s):
        return s if html else escape(s)
    element.set("data-toggle", "popover")
    if title: element.set("title", tos(title))
    element.set("data-placement", "auto right")
    element.set("data-trigger", "hover focus")
    if content: element.set("data-content", tos(content))
    if html: element.set("data-html", "true")


def a2s(element, cls=None):
    """Convert element tree element into HTML string."""
    cls = element.get("class") if cls is None else cls
    return '<a href="%s" class="%s">%s</a>' % (element.get("href"), cls, element.text)


class MarkdownPage(Page):

    def __init__(self, path, website):
        super(MarkdownPage, self).__init__(path, website)
        self.meta = {}
        with io.open(self.path, "rt", encoding="utf-8") as fh:
           string = fh.read()
        lines = string.split("\n")
        #""" Parse Meta-Data and store in Markdown.Meta. """
        # https://github.com/Python-Markdown/markdown/blob/master/markdown/extensions/meta.py
        #self.meta = self._get_meta(string.split("\n"))

        # Note: this logic is able to detect backlinks only if wikilinks syntax is used.
        for m in re.finditer(website.WIKILINK_RE, string):
            token = m.group(1).strip()
            try:
                link = website.get_wikilink(token, self.url)
            except Exception as exc:
                cprint("Exception while trying to handle wikilink `%s` in `%s`" % (token, self.path))
                raise

            if hasattr(link, "get"):
                link_class = link.get("class", "")
                if "citation-wikilink" in link_class:
                    self.citations.add(token)  # TODO Should be name
                elif "topic-wikilink" in link_class:
                    self.topics.add(token) # TODO: Should be name


class HtmlPage(Page):
    def __init__(self, path, website):
        super(HtmlPage, self).__init__(path, website)


class AbinitStats(object):
    """
    This object parses the data stored in statistics.txt and produces the JSON document
    used by plotly to plot the results on the web-site.
    """
    def __init__(self, path):
        self.path = os.path.abspath(path)
        self.parse()

    def parse(self):
        """
        =====================================================================
        Version  Date       Size         Number        Number        Number
                 released   tar.gz       of F90 files  of F90 lines  of tests
                            (10e6Bytes)
        =====================================================================
        4.3      2004 Feb   14.5          726          252602        432
        """
        keys = ("versions", "dates", "targz_sizes", "num_f90files", "num_f90lines", "num_tests")
        self.data = OrderedDict([(k, []) for k in keys])

        with io.open(self.path, "rt", encoding="utf-8") as fh:
            indata = False
            for line in fh:
                indata = indata or line.startswith("4.3")
                if indata:
                    tokens = line.split()
                    values = [tokens[0], "-".join([tokens[1], tokens[2]])] + tokens[3:]
                    for k, v in zip(keys, values):
                        if k not in ("versions", "dates"): v = float(v)
                        self.data[k].append(v)

    def update(self):
        """
        Size of tar.gz      : ls -l *tar.gz    (on shiva)
        Number of F90 files : ls src/*/*.F90 | wc
        Number of F90 lines : cat src/*/*.F90 | wc
        Number of tests     : ls tests/*/Input/t*in | wc
        """
        root = os.path.join(os.path.dirnname(self.path), "..", "..")
        src_dir = os.path.join(root, "src")
        if not os.path.isdir(src_dir):
            raise RuntimeError("Cannot find Abinit src directory. Someone moved statistics.txt file!")
        # Use same shell-based approach to compute stats to be consistent with previous data.
        from subprocess import check_output
        num_f90files = int(check_output(["ls -l %s/*/*.F90 | wc" % src_dir], shell=True))
        num_f90lines = int(check_output(["cat %s/*/*.F90 | wc" % src_dir], shell=True))
        num_tests = int(check_output(["ls %s/tests/*/Input/t*in | wc" % src_dir], shell=True))
        self.parse()

    def json_dump(self, path):
        """Write data in JSON format to file `path`."""
        import json
        with io.open(path, "wt", encoding="utf-8") as fh:
            fh.write(my_unicode(json.dumps(self.data, ensure_ascii=False)))


class HTMLValidator(object):
    """
    This object checks HTML validity by sending requests to <https://validator.w3.org/>

    Arg:
        verbose: Verbosity level.
    """
    def __init__(self, verbose):
        self.verbose = bool(verbose)

    def validate_website(self, dirpath):
        """
        Validate all html pages inside directory `dirpath`. Return exit status.
        """
        print("Validating website in directory:", dirpath)
        retcode = 0
        for top, dirs, files in os.walk(dirpath):
            for f in files:
                if not (f.endswith(".html") or f.endswith(".htm")): continue
                retcode += self.validate_htmlpage(os.path.join(top, f))

        return retcode

    def validate_htmlpage(self, path):
        """Validate html page. Return exit status."""
        # https://bitbucket.org/nmb10/py_w3c
        # import HTML validator and create validator instance
        import urllib
        from py_w3c.validators.html.validator import HTMLValidator
        vld = HTMLValidator()
        num_err, num_ignored, num_warn = 0, 0, 0

        # Ignore error messages containing the following substrings.
        exclude_substrings = [
            "element is obsolete. Use CSS",
            "An “img” element must have an “alt” attribute",
            "query: “|” is not allowed.",
            "Element “div” not allowed as child of element “label” in this",
        ]

        # Also, ignore messages of the form:
        #   'Attribute “autocorrect” not allowed on element “input” at this '
        element2attrs = {
            "input": ["autocorrect", "autocapitalize"],
        }
        for element, attrs in element2attrs.items():
            exclude_substrings.extend('Attribute “%s” not allowed on element “%s”' % (attr, element) for attr in attrs)

        try:
            vld.validate_file(path)
        except urllib.error.URLError as exc:
            cprint("Exception while validating %s.\n%s\nWill try again after 2 sec...\n" % (path, str(exc)), "magenta")
            import time
            time.sleep(2)
            vld.validate_file(path)

        # errors and warnings are list of dicts.
        num_warn += len(vld.warnings)
        if self.verbose:
            for warn in vld.warnings:
                warn["File"] = path
                pprint(warn, indent=4)
                print(80 * "=")

        for err in vld.errors:
            if any(s in err["message"] for s in exclude_substrings):
                num_ignored += 1
                continue
            err["File"] = path
            pprint(err, indent=4)
            print(80 * "=")
            num_err += 1

        cprint("Errors %s, Ignored Errors %s, Warnings: %s in file: %s" % (num_err, num_ignored, num_warn, path),
               color="red" if num_err else "green")

        return num_err
