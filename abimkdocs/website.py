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
import uuid
import pickle
import yaml
import markdown

from collections import OrderedDict, defaultdict
from itertools import groupby
from pprint import pprint
from html2text import html2text
from pybtex.database import parse_file, Entry, BibliographyData
from markdown.util import etree
from pygments import highlight
from pygments.lexers import BashLexer, BibTeXLexer
from pygments.formatters import HtmlFormatter
from doc.tests.pymods.termcolor import cprint
from .variables import Variable


ABINIT_REPO = "/Users/gmatteo/git_repos/abinit/"
if not os.path.exists(ABINIT_REPO):
    raise ValueError("ABINIT_REPO: %s does not exist\n. Please change the global variable in the python module." %
            ABINIT_REPO)


class lazy_property(object):
    """
    lazy_property descriptor

    Used as a decorator to create lazy attributes. Lazy attributes
    are evaluated on first use.
    """

    def __init__(self, func):
        self.__func = func
        from functools import wraps
        wraps(self.__func)(self)

    def __get__(self, inst, inst_cls):
        if inst is None:
            return self

        if not hasattr(inst, '__dict__'):
            raise AttributeError("'%s' object has no attribute '__dict__'"
                                 % (inst_cls.__name__,))

        name = self.__name__
        if name.startswith('__') and not name.endswith('__'):
            name = '_%s%s' % (inst_cls.__name__, name)

        value = self.__func(inst)
        inst.__dict__[name] = value
        return value

    @classmethod
    def invalidate(cls, inst, name):
        """Invalidate a lazy attribute.

        This obviously violates the lazy contract. A subclass of lazy
        may however have a contract where invalidation is appropriate.
        """
        inst_cls = inst.__class__

        if not hasattr(inst, '__dict__'):
            raise AttributeError("'%s' object has no attribute '__dict__'"
                                 % (inst_cls.__name__,))

        if name.startswith('__') and not name.endswith('__'):
            name = '_%s%s' % (inst_cls.__name__, name)

        if not isinstance(getattr(inst_cls, name), cls):
            raise AttributeError("'%s.%s' is not a %s attribute"
                                 % (inst_cls.__name__, name, cls.__name__))

        if name in inst.__dict__:
            del inst.__dict__[name]


def my_unicode(s):
    """Convert string to unicode (needed for py2.7 DOH!)"""
    return unicode(s) if sys.version_info[0] <= 2 else str(s)


def escape(text, tag=None, cls=None):
    """Escape HTML entities in ``text`` string. Enclose new text in ``tag`` if tag with class ``cls``."""
    import cgi
    text = cgi.escape(text, quote=True)
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
            return ", ".join(my_unicode(p) for p in self.persons["author"])
        except KeyError:
            return ""

    def to_markdown(self, bibtex_ui="button"):
        """
        Return markdown string with bibliographic entry.

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
            #s += 'URL: <a href="{url}" target="_blank">{url}</a><br>'.format(url=fields["url"])
            s += 'URL: <{url}>  \n'.format(url=fields["url"])
        elif "doi" in fields:
            doi = fields["doi"]
            doi_root = "https://doi.org/"
            if not doi.startswith(doi_root): doi = doi_root + doi
            s += 'DOI: <{doi}>  \n'.format(doi=doi)
            #s += 'DOI: <a href="{doi}" target="_blank">{doi}</a><br>'.format(doi=doi)

        # Add modal window with bibtex button/link.
        if bibtex_ui is not None:
            assert bibtex_ui in ("link", "button")
            btn, modal = self.get_bibtex_btn_modal(link=bibtex_ui=="link")
            s += btn + modal
        return s

    def to_html(self):
        """Return string with entry in HTML format."""
        return markdown.markdown(self.to_markdown())

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
        self.warnings = []

        # Read mkdocs configuration file.
        # TODO: Should read Abinit version from a centralized file.
        with io.open(os.path.join(self.root, "..", "mkdocs.yml"), "rt", encoding="utf-8") as fh:
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
        from .variables import get_variables_code
        self.variables_code = get_variables_code()

        # Get bibtex references and cast to MyEntry instance.
        # FIXME
        bib_path = os.path.join(self.root, "mkdocs-abiref.bib")
        #bib_path = os.path.join(self.root, "biblio/origin_files/abiref.bib")
        self.bib_data = parse_file(bib_path, bib_format="bibtex")
        for entry in self.bib_data.entries.values():
            entry.__class__ = MyEntry

        # FIXME
        # Read code statistics from file and produce json file.
        #self.abinit_stats = AbinitStats(os.path.join(self.root, "statistics.txt"))
        #self.abinit_stats.json_dump(os.path.join(self.root, "data", "statistics.json"))

        # Build AbinitTestSuite object.
        from doc import tests as tmod
        tests = []
        for t in tmod.abitests.select_tests(suite_args=[], regenerate=True):
            # DO NOT use isinstance to check if ChainOfTests but rely on duck typing.
            # See https://stackoverflow.com/questions/9006740/isinstance-and-type-equivelence-failure-due-to-import-mechanism-python-djan
            #if isinstance(t, ChainOfTests):
            if hasattr(t, "tests"):
                tests.extend(t.tests)
            else:
                tests.append(t)

        # Construct dictionary rpath --> test. Use OrderedDict to have deterministic behaviour.
        self.rpath2test = {}
        for t in tests:
            key = os.path.join(*splitall(t.inp_fname)[-4:])
            self.rpath2test[key] = t
        self.rpath2test = OrderedDict([(k, self.rpath2test[k]) for k in sorted(self.rpath2test.keys())])
        #print(self.rpath2test.keys())

        # Find variables used in tests.
        for var in self.variables_code.iter_allvars():
            assert not hasattr(var, "tests")
            var.tests = []
            var.tests_info = {}

        def test_get_varnames(test, varnames):
            # TODO: This should become a method of BaseTest and must be improved.
            with io.open(test.inp_fname, "rt", encoding="utf-8") as fh:
                s = fh.read()
            vused = [v for v in varnames if v in s]
            return vused

        codes_without_vars = set()
        white_list = set([
            "atompaw",
            "cut3d",
            #"multibinit",
            "mrgddb",
            "mrggkk",
            "band2eps",
            #"ujdet",
            "fold2Bloch",
            "fftprof",
            #"conducti",
            "mrgscr",
            #"macroave",
            "mrgdv",
        ])

        for test in tests:
            vd = self.variables_code.get(test.executable, None)
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
        for var in self.variables_code.iter_allvars():
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

    def new_mdfile(self, dirname, mdname, meta=None):
        """
        Create new markdown file with name `mdname` in directory `dirname`.
        `meta` is an optional dictionary with meta-variables added to the front matter.

        Return: File object.

        .. warning::

            Unicode characters in meta are not supported (annoying portability issue with py2.7)
        """
        path = os.path.join(self.root, dirname, mdname)
        assert path not in self.md_generated
        self.md_generated.append(path)
        if self.verbose: print("Generating markdown file: `%s`" % path)
        rpath = "/" + os.path.relpath(path, self.root)
        if meta is None: meta = {}
        assert "rpath" not in meta
        #meta["rpath"] = rpath
        # Must convert to ASCII to avoid !!python/unicode tags in YAML doc
        # (mkdocs does not use pyaml to parse the front matter).
        if sys.version_info[0] <= 2:
            meta = {k.encode("ascii", errors="strict"): meta[k] for k in meta}
        s = yaml.dump(meta, indent=4, default_flow_style=False).strip().replace(" !!python/unicode", "")
        mdf = io.open(path, "wt", encoding="utf-8")
        mdf.write("---\n%s\n---\n" % s)
        mdf.write(self.do_not_edit_comment)
        mdf.rpath = rpath
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

    def generate_markdown_files(self):
        """Generate markdown files using the data stored in the bibtex file, the abivars file ..."""
        start = time.time()

        # Convert test description from markdown to HTML
        #for t in self.rpath2test.values():
        #    t.description = self.convert_markdown(t.description)

        # Write index.md with the description of the input variables.
        meta = {"description": "Complete list of Abinit input variables"}
        with self.new_mdfile("mkdocs-variables", "index.md", meta=meta) as mdf:
            for code, vd in self.variables_code.items():
                mdf.write("## %s variables   \n\n" % code)
                mdf.write(vd.get_vartabs_html(self, mdf.rpath))
                #mdf.write(2*"\n" + "* * *\n")
            #mdf.write(2*"\n" + "* * *\n")

            # This for the table of variables implemented by Jordan
            mdf.write(self.build_varsearch_html(mdf.rpath))

        # Build markdown page with external parameters.
        with self.new_mdfile("mkdocs-variables", "external_parameters.md") as mdf:
            mdf.write("""\
This document lists and provides the description of the name (keywords) of external parameters
that are not input variables, but that are used in the documentation of other variables,
typically compilation parameters, available libraries, or number of processors.
You can change these parameters at compile or run time usually.

""")
            for pname, info in self.variables_code.external_params.items():
                mdf.write("## %s  \n%s  \n\n" % (pname, info))

        # Build markdown pages for the different sets of variables.
        for code, vd in self.variables_code.items():
            cprint("Generating markdown files with input variables of code: `%s`..." % vd.executable, "green")
            for varset in vd.all_varset:
                var_list = [v for v in vd.values() if v.varset == varset]
                meta = {"description": "%s input variables" % varset}
                with self.new_mdfile("mkdocs-variables", varset + ".md", meta=meta) as mdf:
                    mdf.write("""\
# {varset} input variables

This document lists and provides the description of the name (keywords) of the
{varset} input variables to be used in the input file for the {executable} executable.

""".format(varset=varset, executable=vd.executable))
                    for i, var in enumerate(var_list):
                        mdf.write(var.to_markdown(with_hr=False))

        # Add plotly figures.
        # TODO: Replace it with dot
        if self.deploy:
            with self.new_mdfile("mkdocs-variables", "connections.md", meta={"plotly": True}) as mdf:
                mdf.write("# Dependency graphs  \n")
                mdf.write("""
These graphs show the dependencies of the input variables towards each other.
The colormap gives the number of input variables connected to the node.

""")
                for code, vd in self.variables_code.items():
                    for varset in vd.all_varset:
                        mdf.write("## %s, varset: %s  \n\n" % (code, varset))
                        mdf.write(vd.get_plotly_networkx(varset=varset, include_plotlyjs=False))

        # Write Markdown page with statistics.
        with self.new_mdfile("mkdocs-variables", "varset_stats.md") as mdf:
            mdf.write("""
# Input variables, statistics

This document lists the input variables for ABINIT and three post-processors of ABINIT,
in order of number of occurrence in the input files provided with the package.

""")
            for code, vd in self.variables_code.items():
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

        cprint("Generating Markdown files with topics ...", "green")
        repo_path = os.path.join(ABINIT_REPO, "doc/topics/origin_files/")
        with io.open(os.path.join(repo_path, "list_of_topics.yml"), "rt", encoding="utf-8") as fh:
            self.all_topics = sorted(yaml.load(fh), key=lambda t: t[0].upper())
        with io.open(os.path.join(repo_path, "list_relevances.yml"), "rt", encoding="utf-8") as fh:
            # tribe_name --> description
            self.all_tribes = OrderedDict(yaml.load(fh))

        # datastructures needed for topics index.md
        index_md = ["# Alphabetical list of topics\n"]
        self.howto_topic = {}
        for topic in self.all_topics:
            # Read data from yaml file and generate markdown string.
            with io.open(os.path.join(repo_path, "topic_" + topic + ".yml"), "rt", encoding="utf-8") as fh:
                top = yaml.load(fh)[0]
                title = html2text(top.title)
                introduction = html2text(top.introduction)
                howto = html2text(top.howto).strip().lstrip()
                self.howto_topic[topic] = "How to " + howto if not howto.startswith("to ") else "How " + howto
                tutorials = top.tutorials.strip()

            # Find list of variables associated to this topic
            # Order and group vlist by tribes and write list with links.
            # TODO: Can we have multiple tribes with the same topic?
            related_variables = "No variable associated to this topic."
            vlist = [var for var in self.variables_code.iter_allvars() if topic in var.topic_tribes]
            if vlist:
                lines = []
                def sort_tribes(t):
                    try:
                        return {"basic": 0, "compulsory": 1, "expert": 2, "useful": 3, "internal": 4,
                                "prpot": 5, "prfermi": 6, "prden": 7, "prgeo": 8, "prdos": 9, "prgs": 10,
                                "prngs": 11, "prmisc": 12}[t[0]]
                    except KeyError:
                        raise KeyError("Cannot find tribe `%s` in dict. Add it to sort_tribes with the proper rank."
                                % str(t))

                items = sorted([(v.topic_tribes[topic][0], v) for v in vlist], key=lambda t: sort_tribes(t))
                for tribe, group in sort_and_groupby(items, key=lambda t: t[0]):
                    lines.append("*%s:*\n" % tribe)
                    lines.extend("- %s  %s" % (v.mdlink, v.mnemonics) for (_, v) in group)
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

            # Build markdown text and write md file.
            text = """
This page gives hints on how to {howto} with the ABINIT package.

## Introduction

{introduction}

## Related Input Variables

{related_variables}

## Selected Input Files

{selected_input_files}

""".format(**locals())

            if tutorials:
                text += """\
## Tutorials

{tutorials}""".format(tutorials=html2text(tutorials))

            meta = {"authors": top.authors, "description": "%s Abinit topic" % topic}
            with self.new_mdfile("mkdocs-topics", topic + ".md", meta=meta) as mdf:
            #with self.new_mdfile("topics", topic + ".md", meta=meta) as mdf:
                mdf.write(text)

        # Now write topics index.md (sorted by first character)
        for firstchar, group in sort_and_groupby(self.all_topics, key=lambda t: t[0].upper()):
            index_md.append("## %s" % firstchar)
            index_md.extend("- [[topic:%s|%s]]: %s" % (topic, topic, self.howto_topic[topic]) for topic in group)

        meta = {"description": "List of Abinit topics"}
        with self.new_mdfile("mkdocs-topics", "index.md", meta=meta) as mdf:
            mdf.write("\n".join(index_md))

        # Build page with full list of tests grouped by `suite_name`.
        cprint("Generating Markdown file with tests ...", "green")
        meta = {"description": "List of Abinit tests"}
        items = [(rpath, test) for (rpath, test) in self.rpath2test.items()]

        with self.new_mdfile("mkdocs-developers", "testsuite.md", meta=meta) as mdf:
            for suite_name, group in sort_and_groupby(items, key=lambda t: t[1].suite_name):
                group = list(group)
                mdf.write('## %s  \n\n' % suite_name)
                for i, (rpath, test) in enumerate(group):
                    mdf.write('### [[%s]]   \n\n' % rpath)
                    #mdf.write('### <a href="{rpath}"> {rpath} </a> \n\n'.format(rpath=rpath))
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
        with self.new_mdfile("mkdocs-theory", "bibliography.md", meta=meta) as mdf:
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
                    lines.append(entry.to_markdown())
                except Exception as exc:
                    raise ValueError("Exception while trying to convert bibtex entry `%s`\n%s\n" % (name, str(exc)))
                if citation2pages[name]:
                    lines.append("Referred to in: %s" % ", ".join('[{url}]({url})'.format(url=url)
                        for url in sorted([page.url for page in citation2pages[name]])))

            mdf.write("\n".join(lines))

        # Write acknowledgments page.
        repo_path = os.path.join(ABINIT_REPO, "doc/biblio/origin_files/")
        with io.open(os.path.join(repo_path, "bibfiles.yml"), "rt", encoding="utf-8") as fh:
            for comp in yaml.load(fh):
                if comp.name == "acknow": break
            else:
                raise RuntimeError("Cannot find `acknow` section in components")

        meta = {"description": "Suggested acknowledgments and references"}
        with self.new_mdfile("mkdocs-theory", "acknowledgments.md", meta=meta) as mdf:
            mdf.write("# Acknowledgments  \n")
            mdf.write(html2text(comp.purpose))
            mdf.write(html2text(comp.introduction))

        meta = {"description": "List of PDF files provided by the Abinit documentation"}
        with self.new_mdfile("mkdocs-theory", "documents.md", meta=meta) as mdf:
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
            self.warn("Found markdown files on disk not included in mkdocs.py:\n%s" % "\n".join(diff))
        diff = pages_in_toolbar.difference(pages_on_disk)
        if diff:
            self.warn("Found markdown files in mkdocs.py not present in directories:\n%s" % "\n".join(diff))

    def slugify(self, value):
        """
        Slugify a string, to make it URL friendly. Use same convention as TOC extensions of python markdown.
        """
        from markdown.extensions.toc import slugify
        return slugify(value, separator="-")

    def preprocess_mdlines(self, lines):
        """Preprocess markdown lines."""
        INC_SYNTAX = re.compile(r'^\{%\s*(.+?)\s*%\}')
        new_lines = []
        #print(lines)
        for line in lines:
            m = INC_SYNTAX.search(line)
            if not m:
                new_lines.append(line)
            else:
                args = m.group(1).split()
                action = args.pop(0)
                if self.verbose: print("Triggering action:", action, "with args:", str(args))

                if action == "editor":
                    if len(args) > 1:
                        new_lines.extend(self.editor_tabs(args, title=None).splitlines())
                    else:
                        new_lines.extend(self.editor_panel(args[0], title=None).splitlines())
                elif action == "modal":
                    if len(args) > 1:
                        new_lines.extend(self.modal_with_tabs(args).splitlines())
                    else:
                        new_lines.extend(self.modal_from_filename(args[0]).splitlines())
                elif action == "dialog":
                    if len(args) > 1:
                        new_lines.extend(self.dialogs_from_filenames(args).splitlines())
                    else:
                        new_lines.extend(self.dialog_from_filename(args[0]).splitlines())
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
        Involked by the wikilink extension to implement the wikilink syntax: [namespace:name#fragment|text]

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
                if name.startswith("lesson_"):
                    # Handle [[lesson_gw1|text]]
                    url = "/tutorials/%s" % name.replace("lesson_" , " ", 1).strip()
                    if a.text is None: a.text = name
                    html_classes.append("lesson-wikilink")

                elif name.startswith("topic_"):
                    # Handle [[topic_SelfEnergy|text]]
                    name = name.replace("topic_" , " ", 1).strip()
                    url = "/mkdocs-topics/%s" % name
                    if a.text is None: a.text = "%s topic" % name
                    html_classes.append("topic-wikilink")
                    add_popover(a, content=self.howto_topic[name])

                elif name.startswith("help_"):
                    # Handle [[help_abinit|text]]
                    code = name.replace("help_" , " ", 1).strip()
                    url = "/mkdocs-user-guide/%s" % code
                    if a.text is None: a.text = "%s help file" % code
                    html_classes.append("user-guide-wikilink")

                elif "@" in name:
                    # Handle [[dipdip@anaddb|text]]
                    vname, code = name.split("@")
                    var = self.variables_code[code][vname]
                    url = "/mkdocs-variables/%s#%s" % (var.varset, var.name)
                    if a.text is None: a.text = name
                    html_classes.append("codevar-wikilink")

                elif name in self.variables_code["abinit"]:
                    # Handle link to Abinit variable e.g. [[ecut|text]]
                    var = self.variables_code["abinit"][name]
                    url = "/mkdocs-variables/%s#%s" % (var.varset, var.name)
                    html_classes.append("codevar-wikilink")
                    if a.text is None:
                        a.text = var.name if not var.is_internal else "%%%s" % var.name

                elif name in self.bib_data.entries:
                    # Handle citation
                    ref = self.bib_data.entries[name]
                    url = "/mkdocs-theory/bibliography#%s" % self.slugify(name)
                    content = ref.fields.get("title", "Unknown")
                    if content == "Unknown":
                        self.warn("Entry for %s does not provide title" % name)
                    add_popover(a, content=content) #+ "\n\n" + ref.authors
                    if a.text is None: a.text = "[%s]" % name
                    html_classes.append("citation-wikilink")

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

                elif name in self.variables_code.characteristics:
                    # handle [[ENERGY]] by building internal link to abinit user guide
                    url = "/mkdocs-user-guide/abinit/#32-more-about-abinit-input-variables"
                    if a.text is None: a.text = name

                elif name in self.variables_code.external_params:
                    # handle [[AUTO_FROM_PSP]] by building link with popover
                    content = ("This is an external parameter\n"
                               "typically compilation parameters, available libraries, or number of processors.\n"
                               "You can change these parameters at compile or runtime usually.\n")
                    url = "/mkdocs-variables/external_parameters#%s" % self.slugify(name)
                    if a.text is None: a.text = name
                    add_popover(a, title=self.variables_code.external_params[name], content=content)

                else:
                    self.warn("Don't know how to handle wikilink token `%s` in `%s`" % (token, page_rpath))
                    url, a.text = "FAKE_URL", "FAKE_URL"

        else:
            # namespace is defined
            if namespace in self.variables_code:
                # Handle [[anaddb:asr|text]] or [[abinit:ecut|text]]
                assert fragment is None
                var = self.variables_code[namespace][name]
                url = "/mkdocs-variables/%s#%s" % (var.varset, var.name)
                html_classes.append("codevar-wikilink")
                if a.text is None:
                    a.text = var.name if not var.is_internal else "%%%s" % var.name

            elif namespace == "lesson":
                # Handle [[lesson:wannier90|text]]
                url = "/mkdocs-tutorials/%s" % name
                if a.text is None: a.text = "%s %s" % (name, namespace)
                html_classes.append("lesson-wikilink")

            elif namespace == "help":
                # Handle [[help:optic|text] NB: [[help_codename]] is echoed "codename help file"
                url = "/mkdocs-user-guide/%s" % name
                if a.text is None: a.text = "%s help file" % name
                html_classes.append("user-guide-wikilink")

            elif namespace == "topic":
                # Handle [[topic:BSE|text]]
                url = "/mkdocs-topics/%s" % name
                html_classes.append("topic-wikilink")
                if a.text is None: a.text = "%s_%s" % (namespace, name)
                add_popover(a, content=self.howto_topic[name])

            elif namespace in ("bib", "cite"):
                if namespace == "bib":
                    self.warn("%s in %s is deprecated" % (token, page_rpath))
                # Handle [[bib:biblio|bibliography]]
                if name == "biblio":
                    url = "/mkdocs-theory/bibliography/"
                    if a.text is None: a.text = "bibliography"
                else:
                    # Handle [[bib:Amadon2008]]
                    # TODO bib --> cite
                    try:
                        ref = self.bib_data.entries[name]
                        url = "/mkdocs-theory/bibliography#%s" % self.slugify(name)
                        add_popover(a, content=ref.fields["title"]) #+ "\n\n" + ref.authors
                        if a.text is None: a.text = "[%s]" % name
                        html_classes.append("citation-wikilink")
                    except Exception as exc:
                        self.warn("Exception `%s:%s`\nwhile treating wikilink token: `%s` in `%s`" %
                                (exc.__class__, str(exc), token, page_rpath))
                        url, a.text = "FAKE_URL", "FAKE_URL"

            elif namespace in ("theorydoc", "theory"):
                if namespace == "theorydoc":
                    self.warn("%s in %s is deprecated" % (token, page_rpath))
                # TODO theorydoc --> theory
                # Handle [[theorydoc:mbpt|text]]
                url = "/mkdocs-theory/%s" % name
                html_classes.append("theory-wikilink")
                if a.text is None: a.text = name

            elif namespace == "varset":
                # Handle [[varset:BSE|text]]
                assert fragment is None
                if name == "allvars":
                    url = "/mkdocs-variables/index"
                else:
                    url = "/mkdocs-variables/%s" % name
                if a.text is None: a.text = "%s varset" % name

            elif namespace == "test":
                # Handle [[test:libxc_41]]
                # TODO: Treat subsuite
                tokens = name.split("_")
                suite_name, tnum = "_".join(tokens[:-1]), tokens[-1]
                url = "/tests/%s/Input/t%s.in" % (suite_name, tnum)
                if a.text is None: a.text = "%s[%s]" % (suite_name, tnum)
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
                # Handle [[ac:abiref_gnu_5.3_debug.ac]]
                url = "/build/config-examples/%s" % name
                if a.text is None: a.text = name
                target = "_blank"
                html_classes.append("abifile-wikilink")

            elif namespace == "pdf":
                # Handle [[pdf:howto_chebfi.pdf]]
                try:
                    path = self.pdfs[name]
                    url = "/" + os.path.relpath(path, self.root)
                except KeyError:
                    self.warn("Don't know how to handle wikilink token `%s` in `%s`" % (token, page_rpath))
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
        a.set('href', url)
        if target: a.set('target', target)
        return a

    def build_varsearch_html(self, page_rpath):
        # Build single dictionary mapping varname --> var. Add @code if not abinit.
        allvars = {}
        for code, vd in self.variables_code.items():
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
</button>.""".format(**locals())

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

    def editor_panel(self, path, title=None):
        title = path if title is None else str(title)
        path = os.path.join(self.root, path)
        with io.open(path, "rt", encoding="utf-8") as fh:
            text = escape(fh.read(), tag="pre", cls="small-text")

        return """\
<div class="md-container">
  <div class="panel panel-default">
    <div class="panel-heading">{title}</div>
    <div class="panel-body"><div class="editor" hidden id="{editor_id}">{text}</div></div>
</div></div>""".format(editor_id=gen_id(), **locals())

    def editor_tabs(self, paths, title=None):
        title = "EditorTabs" if title is None else str(title)
        apaths = [os.path.join(self.root, p) for p in paths]

        text_list = []
        for path in apaths:
            with io.open(path, "rt", encoding="utf-8") as fh:
                text_list.append(escape(fh.read(), tag="pre", cls="small-text"))
        tab_ids = gen_id(n=len(text_list))
        editor_ids = gen_id(n=len(text_list))

        # https://codepen.io/wizly/pen/BlKxo?editors=1000
        s = """\
<div class="md-container">
  <div>{title}</div>
    <div>
      <!-- Nav tabs -->
      <ul class="nav nav-pills nav-justified">""".format(title=title)

        for i, (path, tid) in enumerate(zip(paths, tab_ids)):
            s += """<li class="{li_class}"><a href="{href}" data-toggle="pill">{path}</a></li>""".format(
                li_class="active" if i == 0 else " ", href="#%s" % tid, path=path)

        s +=  """\
        </ul>
        <!-- Tab panes -->
        <div class="tab-content clearfix">"""

        for i, (text, tid, editor_id) in enumerate(zip(text_list, tab_ids, editor_ids)):
            s += """\
<div class="tab-pane {active}" id="{tid}">
<div id="{editor_id}" class="editor" hidden>{text}</div></div>
""".format(active="fade in active" if i == 0 else "fade", tid=tid, editor_id=editor_id, text=text)

        return s + 3 * "</div>"


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
        #self.meta = self._get_meta(string.split("\n"))

        # Note: this logic is able to detect backlinks only if wikilinks syntax is used.
        for m in re.finditer(website.WIKILINK_RE, string):
            token = m.group(1).strip()
            try:
                link = website.get_wikilink(token, self.url)
            except Exception as exc:
                cprint("Exception while trying to handle wikilink `%s` in `%s`" % (token, self.path))
                raise
                print(exc)
                continue

            if hasattr(link, "get"):
                link_class = link.get("class", "")
                if "citation-wikilink" in link_class:
                    self.citations.add(token)  # TODO Should be name
                elif "topic-wikilink" in link_class:
                    self.topics.add(token) # TODO: Should be name

        """
        # Add rpath to meta (useful to give the origin of errors in markdown extensions)
        if len(lines) > 1 and lines[0].startswith("---"):
            for i, l in enumerate(lines[1:]):
                if l.startswith("---"):
                    i += 1
                    break
            else:
                raise RuntimeError("Cannot find second `---` marker in %s" % path)

            # Cannot used OrderedDict because markdown parser does not understand !!python/object
            # If py2, convert strings to ascii to avoid !!unicode in meta!
            d = dict(**yaml.load("\n".join(lines[1:i])))

            rpath = "/" + os.path.relpath(path, website.root)
            if "rpath" not in d:
                raise RuntimeError("rpath front matter entry missing in %s" % self.path)
            if d["rpath"] != rpath:
                raise RuntimeError("Wrong rpath in %s.\nExpecting `%s` but got `%s`" % (self.path, rpath, d["rpath"]))

            # This to add rpat automatically to md pages. WARNING: Requires py3k.
            if False and "rpath" not in d or d["rpath"] != rpath and sys.version_info[0] >= 3:
                d["rpath"] = rpath
                del lines[1:i]
                lines.insert(1, yaml.dump(d, indent=4, default_flow_style=False).strip())
                with io.open(self.path, "wt", encoding="utf-8") as fh:
                    fh.write("\n".join(lines))

        #print(self)
        """

    def _get_meta(self, lines):
        """ Parse Meta-Data and store in Markdown.Meta. """
        # https://github.com/Python-Markdown/markdown/blob/master/markdown/extensions/meta.py
        from markdown.extensions.meta import BEGIN_RE, END_RE, META_RE, META_MORE_RE
        meta = {}
        key = None
        if lines and BEGIN_RE.match(lines[0]):
            lines.pop(0)
        while lines:
            line = lines.pop(0)
            m1 = META_RE.match(line)
            if line.strip() == '' or END_RE.match(line):
                break  # blank line or end of YAML header - done
            if m1:
                key = m1.group('key').lower().strip()
                value = m1.group('value').strip()
                try:
                    meta[key].append(value)
                except KeyError:
                    meta[key] = [value]
            else:
                m2 = META_MORE_RE.match(line)
                if m2 and key:
                    # Add another line to existing key
                    meta[key].append(m2.group('value').strip())
                else:
                    lines.insert(0, line)
                    break  # no meta data - done
        return meta


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


def build_modal_with_ajax():
    # https://stackoverflow.com/questions/19663555/bootstrap-3-how-to-load-content-in-modal-body-via-ajax#answer-27718674
    modal = """

<!-- Default bootstrap modal example -->
<div class="modal fade" id="myModal" tabindex="-1" role="dialog" aria-labelledby="myModalLabel" aria-hidden="true">
  <div class="modal-dialog">
    <div class="modal-content">
      <div class="modal-header">
        <button type="button" class="close" data-dismiss="modal" aria-label="Close">
          <span aria-hidden="true">&times;</span>
        </button>
        <h4 class="modal-title" id="myModalLabel">Modal title</h4>
      </div>
      <div class="modal-body">
        ...
      </div>
    </div>
  </div>
</div>

<script>
// Fill modal with content from link href
$(function() {
    $("#myModal").on("show.bs.modal", function(e) {
        var link = $(e.relatedTarget);
        $(this).find(".modal-body").load(link.attr("href"));
    });
});
</script>

"""
    return modal
    #mdf.write("""### <a href="{rpath}" data-toggle="modal" data-target="#myModal" data-remote="false"> {rpath} </a>  \n\n""".format(rpath=rpath))


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