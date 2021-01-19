---
authors: MG, XG
---

This page describes the details of the documentation system of Abinit and how to contribute to it. 

Most of the documentation is written in [Markdown](https://en.wikipedia.org/wiki/Markdown)
a lightweight markup language with plain text 
[formatting syntax](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet).
The documentation includes the User Guide, the Abinit tutorial, the topics, the release notes
as well as the pages with the [input variables](../variables/) and the [bibliographic references](../theory/bibliography.md)
that are generated *automatically* in python from the information reported in 
`~abinit/mkdocs/variables_abinit.py` (and similar files in the same directory for other main executables) and the bibtex 
entries given in the `~abinit/doc/abiref.bib` file.

The website is automatically generated with [MkDocs](http://www.mkdocs.org/)
a static site generator geared towards project documentation.
MkDocs employs [Python-Markdown](https://pypi.python.org/pypi/Markdown) to parse the Markdown documentation
and use a single [YAML](http://en.wikipedia.org/wiki/YAML) configuration file (`mkdocs.yml`) 
defining the organization of the pages on the website.
Navigation bars, header and footer are generated *automatically* by the framework
using the [jinja template engine](http://jinja.pocoo.org/).
Previous versions of the documentation can be consulted using the drop down menu 
at the top-right corner. (==not yet operational==)

MkDocs includes a couple built-in themes as well as various third party themes,
all of which can easily be customized with extra CSS or JavaScript or overridden from the theme directory. 
The Abinit website uses [Mkdocs-Material](http://squidfunk.github.io/mkdocs-material/), a theme
built using Google's [Material Design](https://www.google.com/design/spec/material-design) guidelines.
We also use [fontawesome icons](https://fontawesome.com/) and
[Bootstrap](http://getbootstrap.com/) a popular HTML, CSS, and Javascript framework 
for developing responsive, mobile first projects on the web 
(shrink the browser window to see how the menu and the navigation bars react).

Note that the majority of the Abinit developers do not need to know how to use these technologies
since they will mainly interact with markdown files (plain text files that can be easily modified in the editor)
while Mkdocs will handle the HTML/CSS/Javascript part.

In addition to the basic markdown syntax, the Abinit documentation supports extensions and shortcuts
to ease the writing of hyperlinks and the inclusion of bibliographic citations.
A detailed description of *our markdown dialect* is given in [our markdown page](markdown).
Also [MathJax](https://www.mathjax.org/) for equations in LaTeX is activated, 
and the (few) specificities of its usage in the Abinit docs are explained [in this section](markdown.md#mathjax).

As a net result, Abinit developers can write nice-looking documentation and release notes without having to use 
HTML explicitly while working in an environment that is well-integrated with the Abinit ecosystem 
(the yaml database of input variables, the test suite, bibtex citations).
Adding new content is straightforward: write a new page in Markdown, add the new entry to `mkdocs.yml` 
and finally regenerate the website with MkDocs.


## Getting started

Make sure you are in the top ABINIT directory.
Install the python packages required to build the website with:

```sh
pip install -r requirements.txt --user
```

!!! note
    Python 3.6 is strongly recommended although the code works with python2.7 as well.
    The entire documentation supports Unicode so feel free to use unicode symbols in the docs.

MkDocs comes with a built-in dev-server that lets you preview your documentation as you work on it. 
First issue */*/makemake (actually config/scripts/makemake). Then start *our customized* server 
by running the `./mksite.py serve` command:

```sh
./mksite.py serve

Regenerating database...
Saving database to /Users/gmatteo/git_repos/abidocs/doc/tests/test_suite.cpkl
Initial website generation completed in 9.17 [s]
Generating markdown files with input variables of code: `abinit`...
...
...
INFO    -  Building documentation...
INFO    -  Cleaning site directory
[I 170826 03:37:05 server:283] Serving on http://127.0.0.1:8000
[I 170826 03:37:05 handlers:60] Start watching changes
[I 170826 03:37:05 handlers:62] Start detecting changes
```

Open up `http://127.0.0.1:8000/` in your browser, and you'll see the default home page being displayed.
Note that the generation of the website takes 1-2 minutes but this is a price that must be paid only once.
The web server, indeed, reloads automatically the source files that are modified by the user
so that one can easily change the documentation and inspect the changes in the corresponding HTML files.

!!! tip
    Use `./mksite.py serve --dirtyreload` to enable the live reloading in the development server, 
    but only re-build files that have changed. 
    This option is designed for site development purposes and is **much faster** than the default live reloading.

!!! warning
    The server re-builds automatically the pages generated from changed `.md` files, 
    but not the ones from changed `~abinit/doc/abiref.bib`
    neither from changed `~abinit/abimkdocs/\*.py` . This means that the upgrade of the description of an input variable 
    or a bibtex reference is done by closing the server and reissue the adequate `./mksite.py` command.
    Also, the case of the `.md` files in the `~abinit/doc/topics` directory is similar, as the `.md` source files, 
    prepended with an underscore, must be preprocessed by `./mksite.py` to deliver the `.md` files, without underscore, 
    that are live reloaded.

`./mksite serve` builds the website in a temporary directory. If you need to inspect the HTML files produced 
by the script, use:

    ./mksite.py build

The HTML pages will be available in the `site` directory.
It's also possible to validate all the HTML documents in `site` by using:

    ./mksite.py validate

This command will connect to the [W3C validator service](https://validator.w3.org/) and possible 
errors are printed to the terminal.
To validate a given list of HTML files use:

    ./mksite.py validate site/index.html site/variables/index.html

At present (v8.7.7), many html files are not compliant with the strict html syntax, so this procedure is ==not yet operational==.

Note that the HTML files are produced in a temporary directory, thus they **are not under revision control**.
The real source is represented by the `.md` files and the other `.yml` files. These are the files that can be 
changed by the developers and are therefore under revision control.
The Markdown pages generated by `./mksite.py` are automatically listed in `~abinit/doc/.gitignore` 
and are thus ignored by git.

The `~abinit/doc/mksite.py` script generates the website by converting markdown files into HTML.
The script:

* Starts by creating python objects using the information reported in 
    - the python files in abimkdocs with the input variables,
    - the `~abinit/doc/abiref.bib` for the list of Bibliographic references,
    - the input files contained in `~abinit/tests/*/Input`. 
* Performs initial consistency checks.
* Generate the markdown files for variables, citations, etc.  
* Invoke `mkdocs` to parse the markdown files declared in `mkdocs.yml`
* Expands special strings, of the form <span style="background-color: #E0E0E0;font-size:90%;"> &#91; [namespace:name#section|text] &#93; </span> to create HTML links.
* Creates the needed HMTL files 

The expansion of special strings is documented in the [links section](markdown.md#links). 
It can be used in all the YAML files mentioned below. 
For equations/formulas, [Mathjax](http://docs.mathjax.org/en/latest/mathjax.html) is activated, and allows
to process and visualize LaTeX formulas, see also [this section](markdown.md#MathJax) for further details.


## Writing docs

The markdown files are stored inside the `doc` directory according to the following structure:

```console
├── doc
│   ├── about
│   ├── css
│   ├── developers
│   ├── extra_javascript
│   ├── images
│   ├── variables
│   ├── tests
│   ├── theory
│   ├── topics
│   ├── tutorial
│   └── guide
```

* about: Files with release notes, license
* *css*: Extra CSS files used by the website
* developers: Documentation for developers (documentation howtos, git, coding rules...)
* *extra_javascript*: Extra javascript code used by the website
* *images*: logos and favicon
* variables: files with input variables (automatically generated).
* *tests*: symbolic links to the `~abinit/tests` directory.
* theory: files with theoretical notes
* topics: files with Abinit topics
* tutorial: official Abinit tutorials
* guide: help files for main executables

The directory in *italic* are mainly used to build the website and are not visible outside.
The other directories contain markdown files, each directory is associated to an 
entry in the website menu (see `pages` in `mkdocs.yml`).
The [pages configuration](http://www.mkdocs.org/user-guide/writing-your-docs/) in `mkdocs.yml` 
defines which pages are built by MkDocs and how they appear in the documentation navigation. 

Each directory contains an `index.md` file that is supposed to be a "general" page with an overview 
of the topics treated in that directory.
Some of these `index.md` files are automatically generated by python (e.g. `variables/index.md`)
but others such as `tutorial/index.md` are not. 
So make sure that new documentation pages are properly mentioned and linked in the corresponding 
`index.md` file when a new page is added.

Images and additional material (e.g. scripts) associated to a markdown page are stored 
in the corresponding "assets" directory whose name is constructed from the base name of the markdown page.
For instance the figures used in `tutorial/bse.md` are stored in `tutorial/bse_assets`

<!--
!!! note
    If `True`, use `<page_name>/index.html` style files with hyperlinks to
    the directory. If `False`, use `<page_name>.html` style file with hyperlinks to the file.
    True generates nicer URLs, but False is useful if browsing the output on a filesystem.
-->


### Front matter

Front matter is the first section of the markdown file and must take the form of valid YAML 
document enclosed between triple-dashed lines. Here is a basic example:

```yaml
---
title: Documenting Code Like a Hacker
authors: MG
---
```

Between these triple-dashed lines, you can set predefined variables (see below for a reference) 
or even create custom ones of your own. 
These variables will then be available to the framework.
For instance, the list of authors is reported in the HTML page footer while `title` is added to 
to the HTML meta section.


## Documentation Guide lines

* Each paragraph name should be short enough to fit nicely in the menu, but also long enough to stand 
  on its own to a reasonable extent. 
  The titles set here are used in the navigation menu and the page title that displays in the browser tab. 

* Each page should start with a paragraph that explains what will be covered.
  The first heading in a page should be Heading1 (`#` in Markdown). 
  All others should be in H2  (`##`) and H3 (`###`), only where necessary. 
  If you find yourself wanting to use H4, consider if it's truly necessary.

* Don't use terms like "previous page", etc. because we may add or re-arrange pages in the future. 
  Instead, use a hyperlink to the chapter.
  Also avoid sentences like *If you follow the tutorial, you should go back to the tutorial window now.*

* Number the paragraphs only if really needed: the links in the navigation bar are not readable and besides
  the number will appear in the [permalink](/developers/markdown.md#permalinks).
  This means that you may need to change several links if you decide to add a new section to the page later on.
  Users may want to share or bookmarks links to the Abinit documentation so broken links should be avoided 
  as much as possible.

* The fact the [wikilink syntax](markdown.md#wikilinks) facilitates the inclusion of hyperlinks does not 
  mean that we have to add links *everywhere*. 
  This is especially true in the documentation of the input variables in which it does not make sense to
  put links to the same variable we are describing or the same link over and over again in the same paragraph.

* Text in all uppercase is significantly more difficult to read than lower and mixed case text.
  Writing in all caps is like shouting so use all caps sparingly.


## How to add/modify an input variable 

The yaml database has been replaced by python modules.
The variables are now declared in `~abinit/abimkdocs/variables_CODENAME.py`.

This file consists of a list of dictionaries, each dictionary
contains the declaration of a single variable and the associated documentation in markdown format.
Wikilinks, latex and markdown extensions can be used inside `text`.
Adding a new variable is easy. Edit the python module and add a new item at the end of the list. 
A template is provided.

Remember that `\` is an escaping character in python so the interpreter may 
raise an Exception if you start to add Latex equations in the documentation e.g.

```python
Error: '\alpha' is an unrecognized escape in character string starting ""^\alpha"
```

The solution is simple, declare the string as `raw string` by prepending `r` e.g.:

```python
    text=r"""
The [[spmeth]] input variable defines the method used to calculate the
irreducible polarizability $\chi^{(0)}_{KS}$."""
```

Note that input variables for the executables other than the main abinit (e.g. anaddb, aim, optic) are 
denoted `input_variable_name@executable`, e.g. `dipdip@anaddb`
(this allows to waive the ambiguity with the dipdip input variable used in the main abinit).

After having edited the python modules you **must rerun** `./mksite serve` to see the changes.

!!! important

    Use ```pytest abimkdocs_tests/test_variables.py``` to validate your changes
    before rebuilding the documentation.

    Well, at present (v8.7.7) this script detect too many problems for this procedure to be useful. So this is (==not yet operational==).


## How to add a bibliographic reference

Bibliographic references must be in bibtex format and should provide enough information so that the python code
can generate appropriate links in the website.
The central bibliography database is presently located in `~abinit/doc/abiref.bib`.

For published work with a DOI, we strongly recommend *avoiding* a *cut&paste* from your own bibtex file
to the central bibliography database. 
Indeed, there are units tests to enforce the presence of particular entries in the bibtex document and
your bibtex may not fulfill these requirements.

Providing bibtex data from the publisher site is a better method.
If you know the DOI of the article, it is also possible to use [BetterBib](https://github.com/nschloe/betterbib)
to fetch data from [Crossref](http://www.crossref.org/) and produce the bibtex entry.
BetterBib is available from the Python Package Index, so simply type:

    pip install betterbib

and then use `doi2bibtex` from the command line:

```text
betterbib-doi2bibtex 10.1103/PhysRevLett.96.066402

@article{bibtex,
  author = {Amadon, B. and Biermann, S. and Georges, A. and Aryasetiawan, F.},
  doi = {10.1103/physrevlett.96.066402},
  issn = {0031-9007, 1079-7114},
  journal = {Physical Review Letters},
  month = feb,
  number = {6},
  publisher = {American Physical Society (APS)},
  source = {Crossref},
  title = {{The α−γ Transition} of Cerium Is Entropy Driven},
  url = {http://dx.doi.org/10.1103/physrevlett.96.066402},
  volume = {96},
  year = {2006}
}
```

Add the entry to the bibtex file and use the `FirstAuthorYear` convention for the key
(make sure it's not a duplicated entry).
Note that the bibtex ID must be of the form "FirstauthornameYEAR", e.g. "Amadon2008"
(start with an uppercase letter, then lower case, then four-digit year).
Possibly, a letter might be added in case of ambiguity: e.g. there exists also `Amadon2008a`
Then, build the HTML pages using `./mksite.py serve`.

Run the tests with:

    pytest abimkdocs_tests/test_bibtex.py

with pytest to validate your changes.

In order to refer to a bibliography entry, use the [Wikilink syntax](markdown#wikilinks) with the "cite" namespace.

## Topics

The topic files are written in Markdown and can be found in ~abinit/doc/topics.
The source files start with an underscore e.g. `_AbiPy.md`.
These are **template files** containing the text and two variables:

```
## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}
```

that will be filled by `./mksite.py` by inspecting the database of variables and the tests of the test suite..
A new Markdown file **without underscore** will be generated and included in `mkdocs.yml`.

!!! important

    Developers are supposed to edit the version with the underscore and provide enough
    information in the declaration of the variable and in the `TEST_INFO` section
    so that `./mksite.py` can fill the template.
    Remember to restart `./mksite.py` to see the changes.

## How to a add a new document

In order to add a new tutorial, create a new Markdown file in doc/tutorial and 
register it in `mkdocs.yml` 
Then, build the HTML using `./mksite.py serve` and start to enjoy the Markdown syntax.

The organization of help files and theory documents is very similar to the one for the other tutorials.

### Topics and relevances

Since the beginning of the ABINIT HTML documentation, every input variable 
has been required to belong to a **varset** (set of variables, e.g. `varbas`, `varfil`).
However, starting in Summer 2017, we require every input variable to be also mentioned in at least one of the
documentation **topics** and, for such topic, to be characterized by a **relevance**.

The allowed list of relevancies (a generic list, irrespective of the topic) is declared in
`~abinit/abimkdocs/variables.py`. 
Standard names are:

- *compulsory* (when such input variable **must** be present in the input file when the "feature" of the topic is activated)
- *basic* (when such input variable is usually explicitly specified in the standard usage, although the default might be adequate)
- *useful* (when the default value is used most of the time)
- *expert* (when only expert users should use other values than the default)

Other relevance names have been allowed for specific topics, in which such a classification 
(compulsory/basic/useful/expert) is not a relevant one.

In order to specify the (possibly several) combinations of topic+relevance to which an input variable is attached,
the field "topics" is used inside the `~abinit/abimkdocs/variables_abinit.py` file
(and similar files in the same directory for the other executables).

Some examples:

* for dmatpawu: *DFT+U_useful*
* for mdwall: *MolecularDynamics_expert*
* for gwpara: *parallelism_useful, GW_basic*

The latter is a case where one input variable is associated to two topics, with a different relevance
for topic "parallelism" and topic "GW".


## Release Notes

Release notes are written in Markdown so it is possible to use the [wikilink syntax](markdown.md#wikilinks)
to insert links to new tests, new autoconf files and even links to pull-requests and issues that will redirect
the reader to the Abinit repository on github.
For example, the following markdown text

```md
B.1   
Implementation of algorithms to interpolate the electronic band structure.
See the new input variables [[einterp]], [[nkpath]], and [[prtebands]], 
and the new tests [[test:v8_04]], [[test:libxc_41]].
Added in [[gitsha:f74dba1ed8346ca586dc95fd10fe4b8ced108d5e]]

B.2
Added subsuite syntax [[test:gspw_01]]

C.2  
New versions of Fortran compilers have been integrated in the test farm:

- intel 16.0
- gnu 6.1 and 6.2
- IBM xlf compiler 14.1
- NAG 7.0

Corresponding examples are available in [[ac:abiref_nag_7.0_openmpi.ac]]
```

produces a nice report with links to the features available in the new version:


B.1   
Implementation of algorithms to interpolate the electronic band structure.
See the new input variables [[einterp]], [[nkpath]], and [[prtebands]], 
and the new tests [[test:v8_04]], [[test:libxc_41]].
Added in [[gitsha:f74dba1ed8346ca586dc95fd10fe4b8ced108d5e]].

B.2
Added subsuite syntax [[test:gspw_01]]

C.2  
New versions of Fortran compilers have been integrated in the test farm:

- intel 16.0
- gnu 6.1 and 6.2
- IBM xlf compiler 14.1
- NAG 7.0

Corresponding examples are available in [[ac:abiref_nag_7.0_openmpi.ac]].

!!! important
    We are already using Markdown on gitlab to document our merge requests.
    This means that we can easily integrate all this gitlab documentation
    with the release notes published on the website.


## Variable object

It is the type that contains the other fields.  

abivarname 
: The name of the variable. Note that the name for input variables 
  of the executables anaddb, aim and optic is always finished with @anaddb, @aim or @optic.

characteristics 
: Possibly, a specific characteristics of the input variable. 
  To be chosen among the names in `~abinit/doc/input_variables/origin_files/characteristics.yml`.

commentdefault
: Possibly, some comment about a default value.

commentdims
: Possibly, some comment about the dimension of an array.

defaultval
: Must be an integer or real value, possibly specified using the types presented below (e.g. !multiplevalue)

dimensions
: Either scalar or a list of dimensions, using YML syntax.

excludes
: Possible excluded values

mnemonics
: A longer description of the variable role, in a few words

requires
: The input variable is relevant only if this condition is fulfilled

text
: Free text describing the input variable

topics
: A string, specified in [topics_and_relevances](#topics-and-relevances)

varset
: a unique "set of variables" to which the variable belong. 
  To be chosen among the names in `~abinit/doc/input_variables/origin_files/varsets.yml`.

vartype
: to be chosen among integer, real or string
  If there is no information of a type for a specific variable, its value must be "null".


### MultipleValue object

This is the equivalent to the X * Y syntax in the Abinit parser.

<code>
  X * Y
</code>

will become

```yaml
  !multiplevalue
    number : X
    value : Y
```

If X is null, it means that you want to do *Y (all Y)


### Range object

```yaml
  !range
     start: 1
     stop: N
```

As a default value, it means that the default value is 1, 2, ... N


### ValueWithConditions object

This type allows to specify conditions on values:

```yaml
!valuewithconditions
    defaultval: -[[diemix]]
    '70 < [[iprcel]] and [[iprcel]] < 80': '[[diemix]]'
    '[[iscf]]<10': '[[diemix]]'
    '[[iprcel]]==0': '[[diemix]]'
```

defaultval is the default value if no condition is fulfilled.
As condition, please use strings with the most basic expressions, 
containing <, < =, >, >=, ==, !=, +, -, *, /, etc to allow for further simple parsing !

As a convention, we use "pythonic" way for expressions, so you can use "or", "and" and "in" 
also as <span style="background-color: #E0E0E0;font-size:90%;"> &#91; [varname] &#93; in [1,2,5]</span> for example ...


### ValueWithUnit object

This type allows to specify values with units:

```yaml
!valuewithunit
    units: eV
    value: 100.0
```

means "100 eV".


### Constraints between variables

In the YML file (and via the GUI), there are some constraints between variables that have been introduced.
You can specify "requires: CONDITION" and "excludes: CONDITION" in the YML file 
(or fill the fields requires and excludes in the GUI).

If a varname has "requires: CONDITION", it means that the variable is only relevant when CONDITION is fulfilled.
If a varname has as "excludes: CONDITION", it means that the specification of the variable in the input file forbids 
the CONDITION to be fulfilled.

Pay attention to strings. If it is recognized as string directly, you don't need ticks (' ').
Otherwise, you need to put ticks. 
For example, if you want to use a link as a value, use a link shortcut like <span style="background-color: #E0E0E0;font-size:90%;"> &#91; [abivarname] &#93; </span>. 
See the doc about link shortcuts at [links shortcuts](markdown.md#links). 
