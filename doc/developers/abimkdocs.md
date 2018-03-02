---
authors: MG, XG
---

This page describes the details of the documentation system of Abinit and how to contribute to it. 
A *Proof of concept* website is available at <https://gmatteo.github.io/test_abidocs>.

Most of the documentation is written in [Markdown](https://en.wikipedia.org/wiki/Markdown)
a lightweight markup language with plain text 
[formatting syntax](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet).
The documentation includes the User Guide, the Abinit lessons, the topics, the release notes
as well as the pages with the [input variables](../variables/) and the [bibliographic references](../theory/bibliography.md)
that are generated *automatically* in python from the information reported in `abinit_vars.yml` and the bibtex 
entries given in the `abiref.bib` file.

The website is automatically generated with [MkDocs](http://www.mkdocs.org/)
a static site generator geared towards project documentation.
MkDocs employs [Python-Markdown](https://pythonhosted.org/Markdown) to parse the Markdown documentation
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
We also use [fontawesome icons](http://fontawesome.io/) and
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
and the (few) specificities of its usage in the Abinit docs are explained [in this section](markdown.md#MathJax).

As a net result, Abinit developers can write nice-looking documentation and release notes without having to use 
HTML explicitly while working in an environment that is well-integrated with the Abinit ecosystem 
(the yaml database of input variables, the test suite, bibtex citations).
Adding new content is straightforward: write a new page in Markdown, add the new entry to `mkdocs.yml` 
and finally regenerate the website with MkDocs.


## Getting started

Install the python packages required to build the website with:

```sh
cd ~abinit/docs
pip install -r requirements.txt
```

!!! note
    Python 3.6 is strongly recommended although the code works with python2.7 as well.
    The entire documentation supports Unicode so feel free to use unicode symbols in the docs.

MkDocs comes with a built-in dev-server that lets you preview your documentation as you work on it. 
Make sure you are in `~abinit/docs`, and then start *our customized* server 
by running the `mksite.py serve` command:

```sh
cd ~abinit/docs
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
    Use `mksite.py serve --dirtyreload` to enable the live reloading in the development server, 
    but only re-build files that have changed. 
    This option is designed for site development purposes and is **much faster** than the default live reloading.


`mksite serve` builds the website in a temporary directory. If you need to inspect the HTML files produced 
by the script, use:

    mksite.py build

The HTML pages will be available in the `site` directory.
It's also possible to validate all the HTML documents in `site` by using:

    mksite.py validate

This command will connect to the [W3C validator service](https://validator.w3.org/) and possible 
errors are printed to the terminal.
To validate a given list of HTML files use:

    mksite.py validate site/index.html site/variables/index.html


Note that the HTML files are produced in a temporary directory, thus they **are not under revision control**.
The real source is represented by the `.md` files and the other `.yml` files, these are the files that can be 
changed by the developers and are therefore under revision control).
The Markdown pages generated by `mksite.py` are automatically listed in `~abinit/doc/.gitignore` 
and are thus ignored by git.

The `~abinit/doc/mksite.py` script generates the website by converting markdown files into HTML.
The script:

* Starts by creating python objects using the information reported in 
    - the `abivars.yml file` for the input variables,
    - the `abiref.bib` for the list of Bibliographic references,
    - a set of YAML files, 
    - the input files contained in `tests/*/Input`. 
* Performs initial consistency checks.
* Generate the markdown files for variables, citations, etc.  
* Invoke `mkdocs` to parse the markdown files declared in `mkdocs.yml`
* Expands special strings, of the form `[[namespace:name#section|text]]` to create HTML links.
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
│   ├── data
│   ├── developers
│   ├── extra_javascript
│   ├── images
│   ├── variables
│   ├── tests
│   ├── theory
│   ├── topics
│   ├── tutorials
│   └── user-guide
```

* about: Files with release notes, license
* *css*: Extra CSS files used by the website
* *data*: 
* developers: Documentation for developers (documentation howtos, git, coding rules...)
* *extra_javascript*: Extra javascript code used by the website
* *images*: logos and favicon
* variables: files with input variables (automatically generated).
* *tests*: symbolic links to the `~abinit/tests` directory.
* theory: files with theoretical notes
* topics: files with Abinit topics
* tutorials: official Abinit lessons
* user-guide: help files for main executables

The directory in *italic* are mainly used to build the website and are not visible outside.
The other directories contain markdown files, each directory is associated to an 
entry in the website menu (see `pages` in `mkdocs.yml`).
The [pages configuration](http://www.mkdocs.org/user-guide/writing-your-docs/) in `mkdocs.yml` 
defines which pages are built by MkDocs and how they appear in the documentation navigation. 

Each directory contains an `index.md` file that is supposed to be a "general" page with an overview 
of the topics treated in that directory.
Some of these `index.md` files are automatically generated by python (e.g. `variables/index.md`)
but others such as `lessons/index.md` are not. 
So make sure that new documentation pages are properly mentioned and linked in the corresponding 
`index.md` file when a new page is added.

Images and additional material (e.g. scripts) associated to a markdown page are stored 
in the corresponding "assets" directory whose name is constructed from the base name of the markdown page.
For instance the figures used in `lesson/bse.md` are stored in `lesson/bse_assets`

Note that some files have been renamed when compared to the previous website so
`lessons/lesson_bse.html` becomes `lesson/bse.md`.
This modification is necessary to have nicer URLs, it allows power users to hack the URL more easily.
and reduce the amount of typing required to insert link.

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
title: Blogging Like a Hacker
authors: MG
rpath: /doc/developers/abimkdocs.md
---
```

Between these triple-dashed lines, you can set predefined variables (see below for a reference) 
or even create custom ones of your own. 
These variables will then be available to the framework.
For instance, the list of authors is reported in the HTML page footer while `title` is added to 
to the HTML meta section.
`rpath` gives the location of the markdown page with respect to `~abinit/docs` and it's the only 
required entry in the front matter.


## Documentation Guide lines

==To be discussed==

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

* The fact the [wikilink syntax](markdown.md#wiki-links) facilitates the inclusion of hyperlinks does not 
  mean that we have to add links *everywhere*. 
  This is especially true in the documentation of the input variables in which it does not make sense to
  put links to the same variable we are describing or the same link over and over again in the same paragraph.

* Text in all uppercase is significantly more difficult to read than lower and mixed case text.
  Writing in all caps is like shouting so use all caps sparingly.


## How to add/modify?


### Input variables

Edit the file `abinit_vars.yml` with the `Abivars.jar` GUI:

    java -jar Abivars.jar
  
A window should open, and you can fill the requested information.
To add a new variable, click on `Edit` (upper right) then click on `Add new variable`.

Note that input variables for the executables other than the main abinit (e.g. anaddb, aim, optic) are 
denoted `input_variable_name@executable`, e.g. `dipdip@anaddb`
(this allows to waive the ambiguity with the dipdip input variable used in the main abinit).

After having edited the info related to one input variable 
(see the [file format specifications](abinit_varsyml-format))
you must still **save** the entire file (click on "File" (upper right) then "Save"). 
Just editing one input variable and saving the changes will only go to some internal variable, 
and this is NOT saving your modifications in the YAML storage file.
Then, build the HTML pages using `mksite.py serve`.


### Bibliographic reference

Edit the file `~abinit/doc/abiref.bib` with your preferred editor (it's a standard bibtex file).
Note that the bibtex ID must be of the form "FirstauthornameYEAR", e.g. "Amadon2008" 
(start with an uppercase letter, then lower case, then four-digit year). 
Possibly, a letter might be added in case of ambiguity: e.g. there exists also `Amadon2008a`
Then, build the HTML pages using `mksite.py serve`.

If you know the DOI of the article, it's possible to use [BetterBib](https://github.com/nschloe/betterbib)
to fetch data from [Crossref](http://www.crossref.org/) and produce the bibtex entry.
BetterBib is available from the Python Package Index, so simply type:

    pip install betterbib

and then use `doi2bibtex` from the command line:

```text
doi2bibtex 10.1103/PhysRevLett.96.066402

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


### Topics

The topic HTML files are assembled by `mksite.py` from different sources.

The high-level information is contained in `~abinit/doc/topics/origin_files/topic_NAME.yml`,
where NAME is the name of the topic. 
The first section ("introduction") is also found in this file,
as well as the information on lessons of the tutorial that are relevant for this topic. 
The "text" content of the "introduction" section is in plain HTML.

At variance, the other sections of the `topic_NAME.html` are created from other sources. 
The list of input variables that are relevant to this topics is assembled from the information
given for these input variables, see [Input variables: how_to_add_modify](#how-to-addmodify)
as well as [Topics and tribes](topics-and-tribes), while the list of relevant input files
is assembled from the information in each of these input files (add a line "#%% topics = " in the last section of the input file.). 
The bibliography list is assembled from the references cited in the "introduction" section, 
that use [Shortcuts for Web links](markdown.md#links).

Note the file `default_topic.yml`, whose components are used in case they are not specified explicitly 
in the more specific file `topic_NAME.yml`.
The list of topics is found in the file `~abinit/doc/topics/origin_files/list_of_topics.yml`.

Thus, if you want to a modify "topic" Web page, the "introduction" and "tutorial" sections as well as 
the name and some other high level info can be modified by editing `~abinit/doc/topics/origin_files/topic_NAME.yml`, 
while you have to edit the relevant input variables and input files to modify the next sections. 
You can modify the bibliography section only through a modification of the "introduction" and "tutorial" sections. 

To add a new topic, add the name in `list_of_topics.yml`, and then create a corresponding `topic_NAME.yml` file.
The different components are used by the script generate_doc.py as follows:

* `name`: must be the name of the topics, also used in the name of file (`topic_name.yml`)
* `keyword`: will be inserted in the HTML header to create the Web name of the HTML page, if the default header is used
* `authors`: will be inserted in the copyright, if the default copyright is used
* `howto`: will be inserted in the subtitle of the Web page "How to ... ?" if the default subtitle is ised
* `header`: the HTML header (usually, take the default, that uses the component "keyword")
* `title`: the title that will be echoed in the HTML page
* `subtitle`: the subtitle (usually, take the default, that uses the component "howto")
* `copyright`: the copyright (usually, take the default, that uses the component "authors")
* `links`: list of links to other pages (usually, take the default)
* `introduction`: see above
* `tutorials`: see above
* `end`: final tags, take the default

Then, build the HTML using `mksite.py serve`.


### Lessons

The major part of each lesson HTML file comes from `~abinit/doc/topics/origin_files/lesson_NAME.yml`,
although selected high-level information (name, keyword, author and subtitle) 
is contained in `~abinit/doc/topics/origin_files/lessons.yml`. 
Note the last section of the latter file, that gives a default value for each component of a lesson file 
that would not have been specified in either `lesson_NAME.yml` or the specific section of `lessons.yml`.

The content of `~abinit/doc/topics/origin_files/lesson_NAME.yml` is either:

* an "intro" section and a "body" section,
* or (this is preferred), an "intro" section and a list of sections, each having a "title" and a "body".

In the latter case, a table of content is automatically generated by generate_doc.py .
The latter structure can be seen in `lesson_dmft.yml`.

The "text" content of these section is in plain HTML.
Note that the indentation is important in YAML. 
The "text" lines in `lesson_NAME.yml` must be indented by at least two blanks.

In order to add a new lesson, introduce a new section in `lessons.yml`, and create a new 
`~abinit/doc/topics/origin_files/lesson_NAME.yml`.

Then, build the HTML using `mksite.py`.


### Help files

The structuration for help files is very similar to the one for the lessons of the tutorial.
The major part of comes from `~abinit/doc/users/origin_files/help_NAME.yml`,
although selected high-level information (name, keyword, author and subtitle) 
is contained in `~abinit/doc/topics/origin_files/helps.yml`.

Do not forget to build the HTML using generate_doc.py.


### Theory documents

The structuration for theory documents is very similar to the one for the lessons of the tutorial.
The major part of comes from `~abinit/doc/users/origin_files/theorydoc_NAME.yml`,
although selected high-level information (name, keyword, author and subtitle) 
is contained in `~abinit/doc/topics/origin_files/theorydocs.yml`.

Do not forget to build the HTML using generate_doc.py.


### Topics and tribes

Since the beginning of the ABINIT HTML documentation, every input variable 
has been required to belong to a **varset** (set of variables, e.g. `varbas`, `varfil`).
However, starting in Summer 2017, we require every input variable to be also mentioned in at least one of the
documentation "topics" and, for such topic, to be characterized by a "tribe".

The allowed list of tribes (a generic list, irrespective of the topic) is contained in
`~abinit/doc/topics/origin_files/list_tribes.yml`. 
Standard names are:

- `compulsory` (when such input variable **must** be present in the input file when the "feature" of the topic is activated)
- `basic` (when such input variable is usually explicitly specified in the standard usage, although the default might be adequate)
- `useful` (when the default value is used most of the time)
- `expert` (when only expert users should use other values than the default)

Other tribe names have been allowed for specific topics, in which such a classification 
(compulsory/basic/useful/expert) is not a relevant one.

In order to specify the (possibly several) combinations of topic+tribe to which an input variable is attached,
the field "topics" is used inside the `~abinit/doc/input_variables/generated_doc/abinit_vars.yml` file
(and can be filled thanks to the use of the Abivars.jar GUI).

Some examples:

* for dmatpawu: "DFT+U_useful"
* for mdwall: "MolecularDynamics_expert"
* for gwpara: "parallelism_useful, GW_basic"

The latter is a case where one input variable is associated to two topics, with a different tribe 
for topic "parallelism" and topic "GW".


### Release Notes

Release notes are written in Markdown so it's possible to use the [wikilink syntax](markdown.md#wiki-links)
to insert links to new tests, new autoconf files and even links to pull-requests and issues that will redirect
the reader to the Abinit repository on github.
For example, the following markdown text

```md
B.1   
Implementation of algorithms to interpolate the electronic band structure.
See the new input variables [[einterp]], [[nkpath]], and [[prtebands]], 
and the new tests [[test:v8_04]], [[test:libxc_41]].
Added in [[gitsha:f74dba1ed8346ca586dc95fd10fe4b8ced108d5e]]

C.2  
New versions of Fortran compilers have been integrated in the test farm:

- intel 16.0
- gnu 6.1 and 6.2
- IBM xlf compiler 14.1
- NAG 5.3

Corresponding examples are available in [[ac:abiref_gnu_5.3_debug.ac]]
```

produces a nice report with links to the features available in the new version:

B.1   
Implementation of algorithms to interpolate the electronic band structure.
See the new input variables [[einterp]], [[nkpath]], and [[prtebands]], 
and the new tests [[test:v8_04]], [[test:libxc_41]].
Added in [[gitsha:f74dba1ed8346ca586dc95fd10fe4b8ced108d5e]]

C.2  
New versions of Fortran compilers have been integrated in the test farm:

- intel 16.0
- gnu 6.1 and 6.2
- IBM xlf compiler 14.1
- NAG 5.3

Corresponding examples are available in [[ac:abiref_gnu_5.3_debug.ac]]

!!! important
    We are already using Markdown on gitlab to document our merge requests.
    This means that we can easily integrate all this gitlab documentation
    with the release notes published on the website.


## abinit_vars.yml format

As values in the `~abinit/doc/input_variables/origin_files/abinit_vars.yml` file, 
you can specify numbers, string, arrays, following the [YAML specifications](http://en.wikipedia.org/wiki/YAML).
Several "types" are defined to allow sufficient flexibility in the specifications, as follows.


### !variable object

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
: A string, specified in [topics_and_tribes](#topics-and-tribes)

varset
: a unique "set of variables" to which the variable belong. 
  To be chosen among the names in `~abinit/doc/input_variables/origin_files/varsets.yml`.

vartype
: to be chosen among integer, real or string
  If there is no information of a type for a specific variable, its value must be "null".


### !multiplevalue object

This is the equivalent to the X * Y syntax in fortran.

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


### !range object

```yaml
  !range
     start: 1
     stop: N
```

As a default value, it means that the default value is 1, 2, ... N


### !valuewithconditions object

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
also as `[[varname]] in [1,2,5]` for example ...


### !valuewithunit object

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
For example, if you want to use a link as a value, use a link shortcut like `[[abivarname]]`. 
See the doc about link shortcuts at [links shortcuts](markdown.md#links). 
