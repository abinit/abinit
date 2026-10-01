# ABINIT documentation support code

The `abimkdocs` package provides ABINIT-specific data models and Markdown extensions used to build the documentation website.
Its two central responsibilities are maintaining the input-variable database and preparing the generated Markdown consumed by MkDocs.

For installation, website-building commands, Markdown authoring, and the complete documentation workflow, see the [ABINIT documentation guide](../doc/developers/abimkdocs.md).
The [ABINIT Markdown guide](../doc/developers/markdown.md) describes the supported wikilinks, citations, equations, includes, and other extensions.

## Directory map

| Path | Purpose |
|---|---|
| `variables.py` | Variable data model, validation, formatting, and the database loader. |
| `variables_CODE.py` | Declarative input-variable database for the executable named `CODE`. |
| `website.py` | Website model, generated-page orchestration, wikilink resolution, and consistency checks. |
| `preprocessor.py`, `treepreprocessor.py` | Markdown preprocessing used by the MkDocs integration. |
| `wikilinks.py` | Python-Markdown extension for ABINIT wikilinks. |
| `mdx_figcaption.py` | Python-Markdown extension for figure captions. |

`mksite.py` in the repository root is the command-line frontend.
It constructs a `Website`, asks it to generate the dynamic Markdown files, and then invokes MkDocs.

## Declaring input variables

Variables for each executable are declared in `variables_CODE.py`.
Each module has the following basic structure:

```python
executable = "CODE"

try:
    from abimkdocs.variables import MultipleValue, Range, ValueWithUnit
except ImportError:
    pass

ValueWithConditions = dict
Variable = dict

variables = [
    Variable(
        # One declaration per variable.
    ),
]
```

The declarations intentionally use dictionary-compatible constructors.
When `variables.py` loads a module, `InputVariables.from_pyfile` converts each dictionary into a validated `Variable` object.
Keeping the declaration modules data-oriented also makes them easier to consume from other Python projects.

Add new declarations to the `variables` list and keep the style of the surrounding file.
The more detailed contributor instructions are in [How to add or modify an input variable](../doc/developers/abimkdocs.md#how-to-addmodify-an-input-variable).

### Executable names and variable sets

The spelling of `abivarname` identifies the program that owns a variable:

| Program | `abivarname` | `varset` |
|---|---|---|
| Main `abinit` program | `ecut` | An ABINIT variable group such as `basic`, `dev`, or `eph`. |
| Other program | `asr@anaddb` | The executable name, here `anaddb`. |
| Other program | `dipdip_prt@multibinit` | The executable name, here `multibinit`. |

An unsuffixed name is always interpreted as an `abinit` variable.
Variables belonging to another executable must therefore use the `name@executable` form.
For suffixed variables, the current implementation asserts that the suffix and `varset` are identical.

This convention prevents ambiguity when two executables use the same keyword.
It also determines generated links: ABINIT variables are grouped into pages such as `variables/basic.md`, while variables for another executable are grouped on a page such as `variables/anaddb.md`.

The module-level `executable` value must agree with the filename and with the suffixes used by its declarations.

### Required fields

Every `Variable` requires these fields:

| Field | Meaning |
|---|---|
| `abivarname` | Keyword, qualified with `@executable` when it does not belong to the main program. |
| `varset` | ABINIT variable group, or the executable name for a non-ABINIT variable. |
| `vartype` | One of `integer`, `real`, `string`, or `integer or string`. |
| `topics` | One or more `Topic_Relevance` strings. |
| `dimensions` | `"scalar"` or a description of the array dimensions. |
| `added_in_version` | Version in which the variable appeared, or the established historical marker such as `before_v9`. |
| `text` | Main documentation in ABINIT Markdown. |

Most public variables should also define `defaultval` and `mnemonics`, even though the constructor does not currently require them.

The topic must be listed in `ABI_TOPICS` and the relevance must be a key in `ABI_RELEVANCES`, both defined in `variables.py`.
Typical declarations look like `GroundState_basic`, `GW_expert`, or `LatticeModel_useful`.
Validation rejects unknown variable types, topics, relevances, and characteristics.

### Optional fields

| Field | Use |
|---|---|
| `characteristics` | Special properties written as wikilinks, for example `[[ENERGY]]` or `[[INTERNAL_ONLY]]`. |
| `requires` | Condition under which the variable is relevant. |
| `excludes` | Variables or conditions incompatible with this variable. |
| `commentdefault` | Additional explanation of the default value. |
| `commentdims` | Additional explanation of the dimensions. |
| `alternative_name` | Alias or historical name. |

`[[INTERNAL_ONLY]]` excludes implementation variables from public usage and test-coverage reporting.
Use it only for variables that users should not set directly.

### Dimensions and defaults

Use `dimensions="scalar"` for scalar variables.
Use a list for arrays, with wikilinks when a dimension depends on another input variable:

```python
dimensions=[3, "[[natom]]"]
```

The helper classes in `variables.py` describe defaults and dimensions that cannot be represented by a single literal:

| Helper | Meaning |
|---|---|
| `ValueWithUnit(value=..., units=...)` | A value carrying an explicit physical unit. |
| `ValueWithConditions({...})` | A value or dimension selected by conditions; include a `defaultval` entry. |
| `MultipleValue(number=..., value=...)` | A repeated value, optionally with a variable repetition count. |
| `Range(start=..., stop=...)` | A range of allowed or generated values. |

Refer to existing declarations and their tests before introducing a new representation.
Use wikilinks such as `[[natom]]` in dimensions and conditions so dependency information and links can be derived from the declaration.

### Writing the description

The `text` field accepts the same ABINIT Markdown dialect as ordinary documentation pages.
Use wikilinks for variables, topics, bibliography entries, tests, source files, and other supported namespaces.
See the [Markdown guide](../doc/developers/markdown.md) for the complete syntax.

Use a raw triple-quoted string when the documentation contains LaTeX backslashes:

```python
text=r"""
The tolerance applies to $\alpha$ and depends on [[ecut]].
""",
```

Avoid repeating a link to the variable currently being documented.
`Variable.to_abimarkdown` converts self-references into emphasis to reduce redundant links.

### Complete example

```python
Variable(
    abivarname="example_tol",
    varset="dev",
    vartype="real",
    topics=["Development_expert"],
    dimensions="scalar",
    defaultval=ValueWithUnit(value=1.0e-8, units="Ha"),
    mnemonics="EXAMPLE TOLerance",
    requires="[[optdriver]] == 0",
    added_in_version="10.6.0",
    text=r"""
Defines the example convergence tolerance.
It is relevant when [[optdriver]] is zero.
""",
),
```

For a variable owned by `anaddb`, use `abivarname="example_tol@anaddb"`, `varset="anaddb"`, and place it in `variables_anaddb.py`.

## How the variable database is loaded

`get_codevars()` lazily constructs and caches a `VarDatabase`.
`VarDatabase.from_pyfiles()` scans this directory for every `variables_*.py` module and loads it through `InputVariables.from_pyfile()`.

The resulting structure has two levels:

```text
VarDatabase
  executable name -> InputVariables
                       normalized variable name -> Variable
```

Variable names are normalized to lowercase and stripped of the `@executable` suffix before they become dictionary keys.
The original qualified spelling remains available as `Variable.abivarname`, and `Variable.executable` is derived from that spelling.

`Variable.validate()` checks individual declarations.
The tests additionally check duplicate names, links between dimensions and variables, agreement with variables recognized by the Fortran code, and usage in test inputs.

## How `website.py` works

`website.py` is the integration layer between the repository data and MkDocs.
It does not replace MkDocs; it prepares content and implements ABINIT-specific semantics before MkDocs renders the final HTML.

### Construction

`Website.build(root, deploy, verbose)` creates the process-wide `Website` singleton.
Call `Website.get()` only after this construction step.
For the normal command-line workflow, `root` is the top-level `doc` directory.

During initialization, `Website`:

1. Reads the ABINIT version from `.current_version`.
2. Reads the version-controlled `mkdocs.yml.in` template and substitutes the version.
3. Configures a Python-Markdown parser with the extensions declared by MkDocs.
4. Loads the complete input-variable database.
5. Parses `doc/abiref.bib` and wraps its entries for ABINIT Markdown rendering.
6. Loads source statistics and writes the generated statistics JSON file.
7. Loads the ABINIT test suite and indexes tests by repository-relative input path.
8. Associates documented variables with tests for the corresponding executable.
9. Indexes PDFs used by documentation pages.

The executable attached to each test is important.
Only variables from the matching `InputVariables` collection are considered when collecting usage statistics.
This is another reason to qualify non-ABINIT variable names correctly.

### Generated Markdown

`Website.generate_markdown_files()` creates the dynamic Markdown inputs consumed by MkDocs.
These include:

- the input-variable index and search interface;
- one variable page per ABINIT variable set and per auxiliary executable;
- the external-parameter page;
- topic pages and their variable/test backlinks;
- the test-suite page;
- the bibliography and citation backlinks;
- source statistics and selected installation or configuration pages.

Generated files contain a warning that they must not be edited directly.
`Website.new_mdfile()` records generated paths, writes front matter, and updates `doc/.gitignore` through the generation workflow.
Make changes in the declaration modules, templates, bibliography, tests, or source Markdown instead.

### Wikilinks and preprocessing

Before rendering, `Website` expands ABINIT-specific aliases, includes, and macros.
`get_wikilink()` resolves tokens to the correct relative URL for variables, topics, tutorials, tests, citations, source files, PDFs, and other namespaces.

Common variable forms are:

```text
[[ecut]]
[[abinit:ecut]]
[[asr@anaddb]]
[[anaddb:asr]]
```

The current page path is part of resolution, allowing `Website` to emit relative URLs that MkDocs can validate.
The supported syntax and examples are documented in the [wikilink section of the Markdown guide](../doc/developers/markdown.md#wikilinks).

After generation, `analyze_pages()` parses Markdown and HTML pages, collects citations and other references, and compares pages on disk with the navigation declared in `mkdocs.yml.in`.
Problems are accumulated in `Website.warnings` and are exercised by the website unit tests.

## Development workflow

Install the documentation dependencies and preview the website as described in [Getting started](../doc/developers/abimkdocs.md#getting-started).
After changing a variable declaration or Python support code, restart `mksite.py`; the live MkDocs reload cannot reconstruct the in-memory variable database by itself.

Run the focused tests first:

```sh
pytest abimkdocs_tests/test_variables.py
pytest abimkdocs_tests/test_website.py
```

Then build or serve the documentation:

```sh
./mksite.py build
./mksite.py serve
```

The generated HTML belongs under `site/` and must not be committed.
The dynamic Markdown written under `doc/` is also generated and ignored; edit its source instead.

When adding this README or another file under `abimkdocs`, run:

```sh
./config/scripts/makemake
```

`abimkdocs` is registered as documentation data in `config/specs/buildsys.conf`.
`makemake` refreshes `config/dist/auto-abimkdocs.lst`, which is how `make dist` learns that the file belongs in the source archive.
