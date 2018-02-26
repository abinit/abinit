## Getting started

Install the python packages required to build the static website with:

    $ cd ~abinit/docs
    $ pip install -r requirements.txt

then install the mkdocs plugin with:

    $ cd mkdocs_plugins
    $ python setup.py install

MkDocs comes with a built-in dev-server that lets you preview your documentation as you work on it. 
Make sure you are in `~abinit/docs`, and then start *our customized* server 
by running the `mksite.py` serve command:

```console
$ cd ~abinit/docs
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

## How to add new variables

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
