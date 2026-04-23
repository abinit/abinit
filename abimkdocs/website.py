"""
Classes and functions used to generate the (static) website with the Abinit documentation
from markdown files and the mkdocs static website generator.

For the different between Absolute, Relative, and Root-relative URLs see:

    <http://ifyoucodeittheywill.com/2009/03/absolute-relative-and-root-relative-urls/>
"""

from __future__ import annotations

import datetime
import os
import re
import shutil
import sys
import time
import uuid
from typing import Any, Iterator, Iterable, Callable, TypeVar, IO

#from markdown.util import etree
import xml.etree.ElementTree as etree
from collections import OrderedDict, defaultdict
from itertools import groupby
from pprint import pprint

import markdown
import yaml
from pybtex.database import BibliographyData, Entry, parse_file
from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import BashLexer, BibTeXLexer, PythonLexer

from doc.tests.pymods.termcolor import cprint

from .variables import ABI_RELEVANCES, ABI_TOPICS, lazy_property


def my_unicode(s: Any) -> str:
    """
    Convert a string or object to a unicode string.

    Args:
        s: The object to convert.

    Returns:
        str: The unicode string representation.
    """
    return unicode(s) if sys.version_info[0] <= 2 else str(s)


def escape(text: str, tag: str | None = None, cls: str | None = None) -> str:
    """
    Escape HTML entities in the given text and optionally enclose it in an HTML tag.

    Args:
        text: The string to escape.
        tag: Optional HTML tag name to wrap the text.
        cls: Optional CSS class to apply to the tag.

    Returns:
        str: The escaped (and possibly wrapped) string.
    """
    import html
    text = html.escape(text, quote=True)
    if tag:
        text = '<{tag} class="{cls}">\n{text}\n</{tag}>\n'.format(tag=tag, text=text, cls=cls or "")
    return text


def gen_id(n: int = 1, pre: str = "uuid-") -> str | list[str]:
    """
    Generate universally unique identifiers (UUIDs).

    Args:
        n: Number of IDs to generate.
        pre: Prefix to prepend to each ID.

    Returns:
        str or list: A single ID string if n=1, otherwise a list of ID strings.

    Raises:
        ValueError: If n <= 0.
    """
    # The HTML4 spec says:
    # ID and NAME tokens must begin with a letter ([A-Z a-z]) and may be followed by any number of letters,
    # digits ([0-9]), hyphens ("-"), underscores ("_"), colons (":"), and periods (".").
    if n == 1:
        return pre + str(uuid.uuid4())
    if n > 1:
        return [pre + str(uuid.uuid4()) for i in range(n)]
    raise ValueError("n must be > 0 but got %s" % str(n))


def splitall(path: str) -> list[str]:
    """
    Split a file path into all its component parts.

    Args:
        path: The path string to split.

    Returns:
        list: List of path components.
    """
    allparts = []
    while True:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        if parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        path = parts[0]
        allparts.insert(0, parts[1])
    return allparts


T = TypeVar("T")
def sort_and_groupby(items: Iterable[T], key: Callable[[T], Any], reverse: bool = False) -> Iterator[tuple[Any, Iterator[T]]]:
    """
    Sort a list and group its elements based on a key function.

    Args:
        items: Iterable of items to group.
        key: Function used for sorting and grouping.
        reverse: If True, sort in reverse order.

    Returns:
        itertools.groupby: An iterator of grouped elements.
    """
    return groupby(sorted(items, key=key, reverse=reverse), key=key)


class MyEntry(Entry):
    """
    Extends pybtex Entry with useful methods for generating HTML output.
    See https://bitbucket.org/pybtex-devs/pybtex/
    """
    @lazy_property
    def authors(self) -> str:
        """
        String containing the list of authors formatted for display.
        Returns an empty string if authors are not provided.
        """
        try:
            #return ", ".join(my_unicode(p) for p in self.persons["author"])
            return ", ".join(my_unicode(p).partition(",")[2] + " " +
                             my_unicode(p).partition(",")[0] for p in self.persons["author"])
        except KeyError:
            return ""

    def to_abimarkdown(self, bibtex_ui: str | None = "button") -> str:
        """
        Format the bibliographic entry as a markdown string with ABINIT extensions.

        Args:
            bibtex_ui: UI element for the bibtex source. Possible values are:
                - None: No UI element is added.
                - "link": A simple link to open a modal window with bibtex content.
                - "button": A button to open the modal window.

        Returns:
            str: The markdown representation of the entry.

        Raises:
            TypeError: If the entry type is unknown.
        """
        fields = self.fields
        # Remove {} from (Latex) title.
        # TODO: title must be present
        title = "*%s*" % fields["title"].replace("{", "").replace("}", "")
        authors = self.authors

        # FIXME: enforce format at the level of the unit tests
        if self.type == "article":
            s = f"{authors}  \n{title}  \n"
            if "eprint" in fields:
                s += "{} **{}**, {} ({})".format(fields["journal"], fields.get("archivePrefix", ""),
                        fields["eprint"], fields["year"])
            else:
                s += "{} **{}**, {} ({})".format(fields["journal"], fields["volume"],
                        fields.get("pages", ""), fields["year"])

        elif self.type in ("book", "inproceedings", "incollection"):
            # FIXME Better treatment for incollection
            #editors = ", ".join(str(e) for e in self.persons["editor"]])
            s = f"{authors}  \n{title}  \n"
            s += "{} ({})".format(fields["publisher"], fields["year"])
            if "isbn" in fields:
                s += "isbn: %s" % fields["isbn"]

        elif self.type in ("phdthesis", "mastersthesis"):
            s = "{}  \n{}  \n{} ({})".format(authors, title, fields["school"], fields["year"])

        elif self.type in ("misc", "unpublished"):
            s = "{}  \n{} ({})".format(authors, title, fields["year"])

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
            s += f'DOI: <a href="{doi}" target="_blank">{doi}</a><br>'

        if bibtex_ui is not None:
            # Add modal window with bibtex button/link.
            # Use https://github.com/kylefox/jquery-modal
            bibtex = highlight(self.to_bibtex(), BibTeXLexer(), HtmlFormatter(cssclass="codehilite"))
            modal_id = "modal-id-%s" % self.key
            s = f"""
{s}
<a href="#{modal_id}" rel="modal:open">bibtex</a> <!-- Link to open the modal -->
<div id="{modal_id}" class="modal" style="height: initial; max-height: 90%; max-width: 90%; width: initial;">
{bibtex}
<a href="#" rel="modal:close">Close</a>
</div>
"""
        return s

    def to_html(self) -> str:
        """
        Convert the entry to an HTML string.

        Returns:
            str: The HTML representation of the bibliographic entry.
        """
        return markdown.markdown(self.to_abimarkdown())

    def to_bibtex(self) -> str:
        """
        Convert the entry to a BibTeX formatted string.

        Returns:
            str: The BibTeX string.
        """
        return BibliographyData({self.key: self}).to_string("bibtex")

_WEBSITE = None


class Website:
    """
    Singleton object that stores all information required to generate the ABINIT HTML documentation.

    This class manages input variables, test suites, bibliographic entries, and
    provides methods to convert Markdown to HTML with ABINIT-specific extensions.

    Attributes:
        root: Absolute path to the root directory of the documentation.
        deploy: Boolean indicating if the site is being built for deployment.
        verbose: Verbosity level for debugging and logging.
        md_generated: List of paths to markdown files generated during the build.
        ignored_paths: List of paths to be ignored by the build system.
        warnings: List of warning messages generated during the build.
        mkdocs_config: Dictionary containing the parsed mkdocs.yml configuration.
        markdown: `markdown.Markdown` instance configured with ABINIT extensions.
        codevars: Database of ABINIT input variables.
        bib_data: Bibliographic data parsed from BibTeX files.
        abinit_stats: Statistics about the ABINIT source code.
        abinit_tests: Database of ABINIT tests.
        rpath2test: Mapping from root-relative paths to `Test` objects.
        pdfs: Ordered dictionary mapping file basenames to paths of PDF documents.
    """
    # Regular expression for wikilinks.
    #WIKILINK_RE = r'\[\[([\w0-9_ -]+)\]\]'
    #WIKILINK_RE = r'\[\[([\w0-9_ -\./]+)\]\]'
    WIKILINK_RE = r"\[\[([^\[]+)\]\]"
    #WIKILINK_RE = r'(?![~`])\[\[([^\[]+)\]\]'

    @classmethod
    def build(cls, root: str, deploy: bool, verbose: int) -> Website:
        """
        Build the Website singleton from the given root directory.

        Args:
            root: Root directory of the documentation.
            deploy: Whether to build for deployment.
            verbose: Verbosity level.

        Returns:
            Website: The instantiated Website singleton.

        Raises:
            RuntimeError: If the Website singleton has already been constructed.
        """
        global _WEBSITE
        if _WEBSITE is not None:
            raise RuntimeError("website has been already costructed")
        _WEBSITE = cls(root, deploy, verbose=verbose)
        return _WEBSITE

    @classmethod
    def get(cls) -> Website:
        """
        Get the Website singleton instance.

        Returns:
            Website: The Website instance.

        Raises:
            RuntimeError: If the Website singleton has not been constructed yet.
        """
        global _WEBSITE
        if _WEBSITE is None:
            raise RuntimeError("website must be constructed by calling `Website.build`")
        return _WEBSITE

    def __init__(self, root: str, deploy: bool, verbose: int = 0):
        """
        Initialize the Website instance.

        Args:
            root: Path to the documentation root directory.
            deploy: If True, the site is being built for deployment.
            verbose: Verbosity level.
        """
        start = time.time()
        self.root = os.path.abspath(root)
        self.deploy = bool(deploy)
        self.verbose = verbose
        self.md_generated = []
        self.ignored_paths = []
        self.warnings = []

        # Read mkdocs configuration file.
        # TODO: Should read Abinit version from a centralized file.
        with open(os.path.join(self.root, "..", "mkdocs.yml"), encoding="utf-8") as fh:
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
            with open(test.inp_fname, encoding="utf-8") as fh:
                s = fh.read()
            vused = [v for v in varnames if v in s]
            return vused

        codes_without_vars = set()
        white_list = set([
            "atompaw", "cut3d",
            #"multibinit",
            "mrgddb", "mrggkk", "band2eps",
            #"lruj",
            "fold2Bloch", "fftprof",
            #"conducti",
            "mrgscr",
            #"macroave",
            "mrgdv", "atdep"
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

    def walk_filepath(self) -> Iterator[tuple[str, str]]:
        """
        Iterate over the files stored in the documentation directory.

        Yields:
            tuple: (filename, absolute_path) for each file found.
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

    def warn(self, msg: str) -> None:
        """Print warning message to terminal and save it for future reference."""
        msg = "WARNING: %s" % msg
        self.warnings.append(msg)
        cprint(msg, color="yellow")

    def convert_markdown(self, source: str) -> str:
        """
        Convert a markdown string to a serialized HTML string.

        Args:
            source: The input markdown string.

        Returns:
            str: The converted HTML string.
        """
        self.markdown.reset()
        return my_unicode(self.markdown.convert(source))

    def new_mdfile(self, dirname: str, mdname: str, meta: dict | None = None, with_comment: bool = True, hide_navigation: bool = False, hide_toc: bool = False) -> IO[str]:
        """
        Create a new markdown file and register it in the build system.

        Args:
            dirname: Directory name relative to the root.
            mdname: Name of the markdown file.
            meta: Optional dictionary for YAML front matter.
            with_comment: If True, add a "DO NOT EDIT" comment.
            hide_navigation: If True, hide the navigation sidebar.
            hide_toc: If True, hide the table of contents.

        Returns:
            file: The opened file object for writing.
        """
        dirpath = os.path.join(self.root, dirname)
        if not os.path.isdir(dirpath): os.mkdir(dirpath)
        path = os.path.join(dirpath, mdname)
        assert path not in self.md_generated
        self.md_generated.append(path)
        if self.verbose: print("Generating markdown file: `%s`" % path)

        mdf = open(path, "w", encoding="utf-8")

        if hide_navigation or hide_toc:
            # https://squidfunk.github.io/mkdocs-material/setup/setting-up-navigation/#hide-the-sidebars
            # hide:
            #   - navigation # Hide navigation
            #   - toc        # Hide table of contents
            if meta is None: meta = {}
            assert "hide" not in meta
            meta["hide"] = []
            if hide_navigation: meta["hide"].append("navigation")
            if hide_toc: meta["hide"].append("top")

        if meta is not None:
            # Must convert to ASCII to avoid !!python/unicode tags in YAML doc
            # (mkdocs does not use yaml to parse the front matter).
            if sys.version_info[0] <= 2:
                meta = {k.encode("ascii", errors="strict"): meta[k] for k in meta}
            s = yaml.dump(meta, indent=4, default_flow_style=False).strip().replace(" !!python/unicode", "")
            mdf.write("---\n%s\n---\n" % s)

        if with_comment: mdf.write(self.do_not_edit_comment)

        mdf.rpath = "/" + os.path.relpath(path, self.root)

        return mdf

    def generate_mdindex(self, dirname: str) -> None:
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

    def copy_install_files(self) -> None:
        """
        Copy INSTALL_*.md files from ~abinit to ~abinit/doc and *git ignore* them.
        Files must be included in mkdocs.yml in the `Installation` section.
        """
        top = os.path.abspath(os.path.join(self.root, ".."))
        for f in os.listdir(top):
            if f.startswith("INSTALL_") and f.endswith(".md"):
                src = os.path.join(top, f)
                dest = os.path.join(self.root, f)
                shutil.copy(src, dest)
                self.ignored_paths.append(dest)

    def generate_page_with_ac_examples(self) -> None:
        """Generate markdown pages with all ac examples found in config-examples."""
        dirpath = os.path.join(self.root, "build", "config-examples")
        md_lines = []
        app = md_lines.append
        app("""
# Autoconf examples

This page gathers the autoconf files used by the buildbot testfarm. The different
bots are described in Abinit web site [server matrix](https://github.com/abinit/abinit_web/blob/main/docs/servers.md)
and [builder matrix](https://github.com/abinit/abinit_web/blob/main/docs/builder.md).

""")
        for f in os.listdir(dirpath):
            path = os.path.join(dirpath, f)
            if os.path.isdir(path) or path.endswith(".swp") or path.endswith(".ac"): continue
            app("## %s  " %  f)
            with open(path, encoding="utf-8") as fh:
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

    def generate_markdown_files(self) -> None:
        """
        Main orchestration method to generate all dynamic markdown files for the site.
        Includes variables, tutorial tests, bibliography, and statistics.
        """
        start = time.time()

        self.copy_install_files()
        self.generate_page_with_ac_examples()

        # Write index.md with the description of the input variables.
        meta = {"description": "Complete list of Abinit input variables"}
        with self.new_mdfile("variables", "index.md", meta=meta, hide_toc=True) as mdf:
            mdf.write("\n\n# Input variables \n\n")
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
                    mdf.write(f"""\
# {varset} input variables

This document lists and provides the description of the name (keywords) of the
{varset} input variables to be used in the input file for the {vd.executable} executable.

""")  #.format(varset=varset, executable=vd.executable))

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
            with open(os.path.join(dirpath, "_" + topic + ".md"), encoding="utf-8") as fh:
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
                    group = sorted(group, key=lambda t: t[1].name)
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
            with open(os.path.join(self.root, "topics", "_" + topic + ".md"), encoding="utf-8") as fh:
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
                mdf.write("## %s  \n\n" % suite_name)
                for i, (rpath, test) in enumerate(group):
                    mdf.write("### [[%s]]   \n\n" % rpath)
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
        # TODO: Should profile this part, I believe that most of the time in mkdocs in spent to convert
        # this huge md file to html.
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
The full bibtex file is available [here](../abiref.bib).

""")
            for name in sorted(self.bib_data.entries.keys()):
                entry = self.bib_data.entries[name]
                lines.append("\n\n## **%s** \n\n" % name)
                try:
                    lines.append(entry.to_abimarkdown())
                except Exception as exc:
                    raise ValueError("Exception while trying to convert bibtex entry `%s`\n%s\n" % (name, str(exc)))
                if citation2pages[name]:
                    lines.append("Referred to in: %s" % ", ".join(f"[{url}]({url})"
                        for url in sorted([page.url for page in citation2pages[name]])))

            mdf.write("\n".join(lines))

        meta = {"description": "List of PDF files provided by the Abinit documentation"}
        with self.new_mdfile("theory", "documents.md", meta=meta) as mdf:
            mdf.write("# PDF files  \n")
            for fname, path in self.pdfs.items():
                mdf.write("## %s  \n" % fname)
                rpdf = "/" + os.path.relpath(path, self.root)
                src = os.path.relpath(rpdf, mdf.rpath)
                html = f'<embed src="{src}" type="application/pdf" width="100%" height="480px">\n\n'
                mdf.write(html)

        #for dirname in ["theory"]:
        #    self.generate_mdindex(dirname)
        #topic2pages = defaultdict(list)
        #for page in self.md_pages:
        #    for topic in page.topics:
        #        topic2pages[topic].append(page)

        with open(os.path.join(self.root, ".gitignore"), "w") as fh:
            fh.write("# The following md files have been copied from ~abinit and should be `git ignored`\n")
            fh.writelines(os.path.relpath(p, self.root) + "\n" for p in self.ignored_paths)

            fh.write("# The following md files have been automatically generated and should be `git ignored`\n")
            fh.writelines(os.path.relpath(p, self.root) + "\n" for p in self.md_generated)

        cprint("Markdown files generation completed in %.2f [s]" % (time.time() - start), "green")

    def analyze_pages(self) -> None:
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
                try:
                    self.md_pages.append(MarkdownPage(path, self))
                except Exception as exc:
                    cprint("Exception while parsing markdown page: %s" % path, color="red")
                    raise exc

            elif f.endswith(".html"):
                self.html_pages.append(HtmlPage(path, self))

        self.find_unreferenced_mds()
        cprint("Completed in %.2f [s]" % (time.time() - start), "green")

    def find_unreferenced_mds(self) -> None:
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
        entries = self.mkdocs_config.get("nav")
        #entries = self.mkdocs_config.get("pages")
        #if entries is None: entries = self.mkdocs_config.get("nav") # Old mkdocs syntax
        for entry in entries:
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

    def slugify(self, value: str) -> str:
        """
        Slugify a string, to make it URL friendly. Use same convention as TOC extensions of python markdown.
        """
        from markdown.extensions.toc import slugify
        return slugify(value, separator="-")

    def preprocess_mdlines(self, lines: list[str]) -> list[str]:
        """
        Preprocess the markdown lines before parsing.

        Args:
            lines: List of markdown lines as strings.

        Returns:
            list: Preprocessed markdown lines.
        """
        lines = self._preprocess_aliases(lines)
        lines = self._preprocess_include(lines)
        lines = self._preprocess_macros(lines)
        return lines

    def _preprocess_macros(self, lines: list[str]) -> list[str]:
        """
        Handle macro substitutions in markdown lines (e.g., [TUTORIAL_README]).

        Args:
            lines: List of markdown lines.

        Returns:
            list: Markdown lines with macros expanded.
        """
        tutorial_readme = """

!!! note

    Supposing you made your own installation of ABINIT, the input files to run the examples
    are in the *~abinit/tests/* directory where *~abinit* is the **absolute path** of the abinit top-level directory.
    If you have NOT made your own install, ask your system administrator where to find the package,
    especially the executable and test files.

    In case you work on your own PC or workstation, to make things easier, we suggest you define
    some handy environment variables by executing the following lines in the terminal:

    ```bash
    export ABI_HOME=Replace_with_absolute_path_to_abinit_top_level_dir # Change this line
    export PATH=$ABI_HOME/src/98_main/:$PATH      # Do not change this line: path to executable
    export ABI_TESTS=$ABI_HOME/tests/             # Do not change this line: path to tests dir
    export ABI_PSPDIR=$ABI_TESTS/Pspdir/  # Do not change this line: path to pseudos dir
    ```

    Examples in this tutorial use these shell variables: copy and paste
    the code snippets into the terminal (**remember to set ABI_HOME first!**) or, alternatively,
    source the `set_abienv.sh` script located in the *~abinit* directory:

    ```sh
    source ~abinit/set_abienv.sh
    ```

    The 'export PATH' line adds the directory containing the executables to your [PATH](http://www.linfo.org/path_env_var.html)
    so that you can invoke the code by simply typing *abinit* in the terminal instead of providing the absolute path.

    To execute the tutorials, create a working directory (`Work*`) and
    copy there the input files of the lesson.

    Most of the tutorials do not rely on parallelism (except specific [[tutorial:basepar|tutorials on parallelism]]).
    However you can run most of the tutorial examples in parallel with MPI, see the [[topic:parallelism|topic on parallelism]].
"""

        hdiago_readme = """

!!! important

    In this lesson, we rely on the iterative KS eigensolvers to compute empty states.
    However, when a large number of unoccupied states is required, a direct diagonalization
    of the KS Hamiltonian is generally more efficient.
    This can be done by setting:

        [[optdriver]] = 6
        [[gwr_task]] = "HDIAGO"

    and, when available, enabling the ELPA library for optimal performance
    For additional information, please consult the [gwr_intro](/tutorial/gwr_intro) page

"""
        new_lines = []
        for line in lines:
            if "[TUTORIAL_README]" in line:
                new_lines.extend(tutorial_readme.splitlines())
            elif "[HDIAGO_README]" in line:
                new_lines.extend(hdiago_readme.splitlines())
            else:
                new_lines.append(line)

        return new_lines


    def _preprocess_aliases(self, lines: list[str]) -> list[str]:
        """
        Handle aliases of the form |token|.

        Args:
            lines: List of markdown lines.

        Returns:
            list: Markdown lines with aliases replaced.
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
            if self.verbose: print("Returning full match:", matchobj.group(0))
            return matchobj.group(0)

        alias_syntax = re.compile(r"[^`\$]\|(?P<key>\w+)\|")
        #alias_syntax = re.compile(r"(?!`+)\|(?P<key>\w+)\|")
        return [re.sub(alias_syntax, repl, line) for line in lines]

    def _preprocess_include(self, lines: list[str]) -> list[str]:
        """
        Handle custom ABINIT inclusion syntax: {% action ... %}.

        Args:
            lines: List of markdown lines.

        Returns:
            list: Markdown lines with inclusions processed.

        Raises:
            ValueError: If an unknown action is encountered.
        """
        inc_syntax = re.compile(r"^\{%\s*(.+?)\s*%\}")
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
                if action == "dialog":
                    if len(args) > 1:
                        new_lines.extend(self.dialogs_from_filenames(args).splitlines())
                    else:
                        new_lines.extend(self.dialog_from_filename(args[0]).splitlines())
                elif action == "include":
                    with open(args[0], encoding="utf-8") as f:
                        new_lines.extend([l.rstrip() for l in f])
                else:
                    raise ValueError("Don't know how to handle action: `%s` in token: `%s`" % (action, m.group(1)))

        return new_lines

    @staticmethod
    def parse_wikilink_token(token: str) -> tuple[str | None, str | None, str | None, str | None]:
        """
        Parse a wikilink token of the form `namespace:name#fragment|text`.

        Args:
            token: The string enclosed between square brackets.

        Returns:
            tuple: (namespace, name, fragment, text). Individual entries are None if not present.
        """
        #args = ""
        #if "||" in token:
        #    token, args = token.split("||")

        if token.startswith(":") and token.endswith(":"):
            # Handle special cases with POSIX regex e.g. [[:digit:]]
            return None, None, None, token

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

    def get_wikilink(self, token: str, page_rpath: str) -> etree.Element | str:
        """
        Resolve a wikilink token and return an HTML anchor element.

        Supported syntaxes:
        - [[www.google.com|text]] (external)
        - [[ecut]] (Abinit variable)
        - [[tutorial:wannier90]] (namespace:name)
        - [[#fragment|text]] (internal fragment)

        Args:
            token: The raw string within [[...]].
            page_rpath: The root-relative path of the current page.

        Returns:
            etree.Element or str: The HTML anchor element or an empty string on warning.
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
            a.set("href", url)
            a.set("target", "_blank")
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
                #assert fragment is not None
                url = ""
                if a.text is None: a.text = fragment
            else:

                if "@" in name:
                    # Handle [[dipdip@anaddb|text]]
                    vname, code = name.split("@")
                    var = self.codevars[code][vname.lower()]
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
                    if "Pspdir" in name:
                        # Handle [[~abinit/tests/Pspdir/6c.lda.atompaw]]
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
                    add_popover(a, content=content, title=self.codevars.all_external_params[name])

                else:
                    self.warn("Don't know how to handle wikilink token `%s` in `%s`" % (token, page_rpath))
                    url, a.text = "FAKE_URL", "FAKE_URL"

        # namespace is defined
        elif namespace in self.codevars:
            # Handle [[anaddb:asr|text]] or [[abinit:ecut|text]]
            assert fragment is None
            var = self.codevars[namespace][name.lower()]
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
            # Handle [[cite:biblio]]
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
            # Handle [[test:libxc_41]] (syntax for suite) [[test:gswvl_01]] (syntax for subsuite)
            tokens = name.split("_")
            prefix, tnum = "_".join(tokens[:-1]), tokens[-1]
            if prefix in self.abinit_tests.all_subsuite_names:
                # [[test:gspw_01]]  --> Need to get the name of suite from subsuite.
                suite_name = self.abinit_tests.suite_of_subsuite(prefix).name
                url = "/tests/%s/Input/t%s.abi" % (suite_name, name)
            else:
                # [[test:libxc_41]]
                url = "/tests/%s/Input/t%s.abi" % (prefix, tnum)

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
            a.set("href", url)
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
            # Hacking previous implementation to make it work with new mkdocs (?)
            #url = os.path.relpath(url, page_rpath)
            if end: url = "%s#%s" % (url, end)
            #print("url", url)

        if self.verbose: print("token", token, "page_rpath", page_rpath, "url", url)

        a.set("href", url.strip())
        if target: a.set("target", target)
        return a

    def build_varsearch_html(self, page_rpath: str) -> str:
        """
        Build the HTML for the interactive variable search page.

        Args:
            page_rpath: The root-relative path of the current page.

        Returns:
            str: The generated HTML fragment.
        """
        # Build single dictionary mapping varname --> var. Add @code if not abinit.
        allvars = {}
        for code, vd in self.codevars.items():
            allvars.update({v.abivarname: v for v in vd.values()})

        tabs = "\n".join(f"""\
<a class="TabLetterLink" href="#{cap_char}" onClick="openLetter(event,'{cap_char}')" id="click{cap_char}">{cap_char}</a>""" for cap_char in sorted(set([k[0].upper() for k in allvars])))

        html_vars = ""
        for char, group in sort_and_groupby(list(allvars.items()), key=lambda t: t[0][0].upper()):
            lis = "\n".join("<li>{link}</li>".format(
                link=var.internal_link(self, page_rpath, label=var.abivarname, cls="small-grey-link")) for _, var in sorted(group))

            html_vars += f"""
<li><ul id="{char}" class="TabContentLetter">
<li class="HeaderLetter">{char}</li> {lis} </ul></li>""" #.format(char=char, lis=lis)

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
        return f"""

Enter any string to search in the database.
Enter empty string to show all variables for all executables.
Variables specific to `executable` are denoted with the syntax `varname@executable` e.g. `a2fsmear@anaddb`.
Enter `@anaddb` in the search bar to show only the variables of `anaddb`.

{search_form}

<div class="TabsLetter"> {tabs} </div>

<ul id="Letters">
{html_vars}
</ul>"""

    def dialogs_from_filenames(self, paths: list[str]) -> str:
        """Build customized jquery dialog to show the content of a list of filepaths."""
        buttons, dialogs = [], []
        for path in paths:
            btn, dialog = self.dialog_from_filename(path, ret_btn_dialog=True)
            buttons.append(btn)
            dialogs.append(dialog)

        button_group = '<div class="text-center"><div class="btn-group-vertical">\n%s\n</div></div>' % "\n".join(buttons)
        return button_group + "\n".join(dialogs)

    def dialog_from_filename(self, path: str, title: str | None = None, ret_btn_dialog: bool = False) -> str | tuple[str, str]:
        """
        Create a jQuery UI dialog to display the content of a file.

        Args:
            path: Relative path to the file.
            title: Optional title for the dialog.
            ret_btn_dialog: If True, return a (button, dialog) tuple instead of a concatenated string.

        Returns:
            str or tuple: The generated HTML (and possibly the button separately).
        """
        abs_path = os.path.join(self.root, path)

        title = path if title is None else title
        with open(os.path.join(self.root, path), encoding="utf-8") as fh:
            if path.endswith(".abi") or path.endswith(".in"):
                text = highlight(fh.read(), BashLexer(), HtmlFormatter(cssclass="codehilite small-text"))
            elif path.endswith(".py"):
                text = highlight(fh.read(), PythonLexer(), HtmlFormatter(cssclass="codehilite small-text"))
            else:
                text = escape(fh.read(), tag="pre", cls="small-text")

        btn_id, dialog_id = gen_id(n=2)

        dialog = f"""\
<div id="{dialog_id}" class="my_dialog" title="{title}" hidden><div>{text}</div></div>

<script> $(function() {{ abidocs_jqueryui_dialog("#{dialog_id}", "#{btn_id}") }}); </script>
"""

        button = f'<button id="{btn_id}" class="md-button md-button--secondary"> View {path} </button>'

        if not ret_btn_dialog:
            button = '<div class="text-center">%s</div>' % button
            return button + dialog
        return button, dialog


class Page:
    """
    Base class representing a documentation page.

    Attributes:
        path: Absolute path to the source file.
        website: Reference to the `Website` singleton.
        citations: Set of bibliographic keys cited in this page.
        topics: Set of topic names associated with this page.
    """

    def __init__(self, path: str, website: Website):
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


def add_popover(element: etree.Element, content: str | None = None, title: str | None = None, html: bool = False) -> None:
    """
    Attach a Tippy.js popover to an etree Element.

    Args:
        element: The etree Element instance (usually an anchor).
        content: The text content for the popover.
        title: Optional title for the popover.
        html: If True, treat content as raw HTML (otherwise escape it).
    """
    # NB: Unfortunately, cannot subclass etree.Element in py2.7.
    def tos(s):
        return s if html else escape(s)

    element.set("data-toggle", "popover")
    if title: element.set("title", tos(title))
    if content:
        element.set("data-tippy-content", tos(content))


def a2s(element: etree.Element, cls: str | None = None) -> str:
    """Convert element tree element into HTML string."""
    cls = element.get("class") if cls is None else cls
    return '<a href="%s" class="%s">%s</a>' % (element.get("href"), cls, element.text)


class MarkdownPage(Page):
    """Representation of a page generated from a markdown file."""

    def __init__(self, path: str, website: Website):
        super().__init__(path, website)
        self.meta = {}
        with open(self.path, encoding="utf-8") as fh:
           string = fh.read()

        #lines = string.split("\n")
        #""" Parse Meta-Data and store in Markdown.Meta. """
        # https://github.com/Python-Markdown/markdown/blob/master/markdown/extensions/meta.py
        #self.meta = self._get_meta(string.split("\n"))

        # Note: this logic is able to detect backlinks only if wikilinks syntax is used.
        for m in re.finditer(website.WIKILINK_RE, string):
            token = m.group(1).strip()
            try:
                link = website.get_wikilink(token, self.url)
            except Exception:
                cprint("Exception while trying to handle wikilink `%s` in `%s`" % (token, self.path))
                raise

            if hasattr(link, "get"):
                link_class = link.get("class", "")
                if "citation-wikilink" in link_class:
                    self.citations.add(token)  # TODO Should be name
                elif "topic-wikilink" in link_class:
                    self.topics.add(token) # TODO: Should be name


class HtmlPage(Page):
    """Representation of a page generated from an HTML file."""
    def __init__(self, path: str, website: Website):
        super().__init__(path, website)


class AbinitStats:
    """
    Parser for ABINIT source code statistics (lines of code, number of files, etc.).
    Exports data to JSON format for visualization.
    """
    def __init__(self, path: str):
        self.path = os.path.abspath(path)
        self.parse()

    def parse(self) -> None:
        """
        Parse the content of the statistics file.
        """
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

        with open(self.path, encoding="utf-8") as fh:
            indata = False
            for line in fh:
                indata = indata or line.startswith("4.3")
                if indata:
                    tokens = line.split()
                    values = [tokens[0], "-".join([tokens[1], tokens[2]])] + tokens[3:]
                    for k, v in zip(keys, values):
                        if k not in ("versions", "dates"): v = float(v)
                        self.data[k].append(v)

    def update(self) -> None:
        """
        Recompute statistics by scanning the source directory.
        Calculates number of F90 files, lines of code, and number of tests.
        """
        root = os.path.join(os.path.dirnname(self.path), "..", "..")
        src_dir = os.path.join(root, "src")
        if not os.path.isdir(src_dir):
            raise RuntimeError("Cannot find Abinit src directory. Someone moved statistics.txt file!")
        # Use same shell-based approach to compute stats to be consistent with previous data.
        from subprocess import check_output
        num_f90files = int(check_output(["ls -l %s/*/*.F90 | wc" % src_dir], shell=True))
        num_f90lines = int(check_output(["cat %s/*/*.F90 | wc" % src_dir], shell=True))
        num_tests = int(check_output(["ls %s/tests/*/Input/t*abi | wc" % src_dir], shell=True))
        self.parse()

    def json_dump(self, path: str) -> None:
        """
        Export the collected statistics to a JSON file.

        Args:
            path: Path where the JSON file will be written.
        """
        import json
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(my_unicode(json.dumps(self.data, ensure_ascii=False)))


class HTMLValidator:
    """
    Validates HTML pages using the W3C validator service.

    Attributes:
        verbose: Verbosity level for validation reports.
    """
    def __init__(self, verbose: int | bool):
        self.verbose = bool(verbose)

    def validate_website(self, dirpath: str) -> int:
        """
        Validate all HTML files within a directory and its subdirectories.

        Args:
            dirpath: Path to the directory containing HTML files.

        Returns:
            int: The total number of validation errors found.
        """
        print("Validating website in directory:", dirpath)
        retcode = 0
        for top, dirs, files in os.walk(dirpath):
            for f in files:
                if not (f.endswith(".html") or f.endswith(".htm")): continue
                retcode += self.validate_htmlpage(os.path.join(top, f))

        return retcode

    def validate_htmlpage(self, path: str) -> int:
        """
        Validate a single HTML file.

        Args:
            path: Path to the HTML file.

        Returns:
            int: Number of validation errors found in the file.
        """
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
            exclude_substrings.extend("Attribute “%s” not allowed on element “%s”" % (attr, element) for attr in attrs)

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
