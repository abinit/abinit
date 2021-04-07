## Getting started with mkdocs

!!! warning
    
    The code supports py2.7 and python3.6 but py3k is **strongly** suggested
    especially when building the final version before deploying 
    as py3k has native support for unicode.

Install the python packages required to build the static website with:

    pip install -r requirements.txt

then install the mkdocs plugin with:

    cd abimkdocs_plugin
    python setup.py install

If you use conda, you may want to create a new environment based on python3.6 with:

    conda create -n abinit-abimkdocs-2 python=3.6
    source activate abinit-abimkdocs

and then install the packages with pip (see above commands).

MkDocs comes with a built-in dev-server that lets you preview your documentation as you work on it. 
Make sure you are in `~abinit`, and then start *our customized* server 
by running the `mksite.py` serve command with the `--dirtyreload` option:

```console
$ ./mksite.py serve --dirtyreload
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

Use:

    $ ./mksite.py --help

to get the list of commands and:

    $ ./mksite.py COMMAND --help

to get the documentation of `COMMAND`.

## Notes

topics are still in Yaml format
wikilink syntax: <span style="background-color: #E0E0E0;font-size:90%;"> &#91; [names|text] &#93; </span> will become <span style="background-color: #E0E0E0;font-size:90%;"> &#91; [text|namespace] &#93; </span>.

## Troubleshooting

If you get error message about ASCII being used as default encoding on your machine add:

    export LC_ALL=en_US.UTF-8  
    export LANG=en_US.UTF-8

to your ~/.bashrc and source it

## How to add new variables

The yaml database has been replaced by python modules.
The variables are now declared in `~abinit/abimkdocs/variables_CODENAME.py`.

This file consists of a list of dictionaries, each dictionary
contains the declaration of a single variable and the associated documentation in markdown format.
Wikilinks, latex and markdown extensions can be used inside `text`.
This is, for example, the declaration of the `accuracy` variable in python:

```python
Variable(
    abivarname="accuracy",
    varset="basic",
    vartype="integer",
    topics=['Planewaves_basic', 'SCFControl_basic'],
    dimensions="scalar",
    defaultval=0,
    mnemonics="ACCURACY",
    text="""
Allows to tune the accuracy of a calculation by setting automatically the
variables according to the following table:

accuracy         | 1         | 2          | 3            | 4            | 5         | 6
---              |---        |---         |---           |---           |---        |---
[[ecut]]         | E_min     | E_med      | E_med        | E_max        | E_max     | E_max
[[pawecutdg]]    | ecut      | ecut       | 1.2 * ecut   | 1.5 * ecut   | 2 * ecut  | 2 * ecut
[[fband]]        | 0.5       | 0.5        | 0.5          | 0.5          | 0.75      | 0.75
...

For a parallel calculation, [[timopt]] is enforced to be 0.
...
""",
),
```

Adding a new variable is easy. Edit the python module and add a new item at the end of the list. 
A template is provided ...


## How to add a new bibtex entry

Citations must be in bibtex format and provide enough information so that the python code
can generate appropriated links in the website.
For published work with a DOI, we strongly recommend *avoiding* a `cut&paste` from your bibtex files
(there are units tests to enforce the presence of particular entries in the bibtex document and
your bibtex may not fullfill these requirements).

A much better solution is to use BetterBib and the DOI of the article to fetch data 
from Crossref and produce the bibtex entry. 
BetterBib is available from the Python Package Index, so simply type:

    pip install betterbib

and then use doi2bibtex from the command line:

    doi2bibtex 10.1103/PhysRevLett.96.066402

Add the entry to the bibtex file and use the `FirstAuthorYear` convention for the key 
(make sure it's not a duplicated entry).
Note that the bibtex ID must be of the form "FirstauthornameYEAR", e.g. "Amadon2008" 
(start with an uppercase letter, then lower case, then four-digit year). 
Possibly, a letter might be added in case of ambiguity: e.g. there exists also `Amadon2008a`
Then, build the HTML pages using `mksite.py serve`.

Run the tests in `./tests/test_bibtex.py` with pytest (see next section) to validate your changes.

## Running the unit tests

Unit tests are located in the ./tests directory. 
To execute the tests, install `pytest` with:

    pip install pytest

and then:

    pytest -v ./tests/

Use 

    pytest -v ./tests/test_variables.py

to execute a particular module.

## Checking links with linkchecker

Build the website with:

    ./mksite.py build

then use

    linkchecker site/index.html > links.err

!!! important

    For the time being, linkchecker python2.7
