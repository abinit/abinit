---
name: abinit-documentation-check
description: Validate ABINIT documentation, variables, citations, wikilinks, generated Markdown, navigation, and the MkDocs site. Use after changes under doc/, abimkdocs/, abinit_theme/, abimkdocs_plugin/, mkdocs.yml.in, or mksite.py; do not use for unrelated Python tests.
---

# Check ABINIT documentation

Work from the top-level ABINIT repository.
When editing files under `doc/`, read `doc/AGENTS.md` completely and follow its writing,
linking, template, and asset conventions.

Use the requested Conda environment for every Python and MkDocs command:

```sh
source /Users/giantomassi/miniconda3/etc/profile.d/conda.sh
conda activate env3.14
```

Run activation and the dependent command in the same shell invocation.

## Choose checks by change type

Start with the smallest relevant check, then broaden validation when shared machinery or navigation
changed.

| Changed content | Minimum check |
|---|---|
| Input-variable definitions or rendering | `pytest abimkdocs_tests/test_variables.py -q` |
| `abimkdocs/abiref.bib` or citation handling | `pytest abimkdocs_tests/test_bibtex.py -q` |
| Wikilinks, Markdown extensions, `Website`, generated pages, or MkDocs integration | `pytest abimkdocs_tests/test_website.py -q` |
| Shared code under `abimkdocs/` or `abimkdocs_plugin/` | `pytest abimkdocs_tests -q` |
| Markdown links, navigation, templates, theme, JavaScript, CSS, or deployment configuration | Full strict site build |

For a broad documentation change, run the complete `abimkdocs_tests` suite before building the site.

## Fast checks

Run from the repository root:

```sh
pytest abimkdocs_tests -q
```

These tests validate variable metadata, BibTeX entries, wikilink behavior, and the Python website
machinery.
They do not prove that every Markdown page is present in navigation or that MkDocs accepts every link.

When a test fails, report the failing table or assertion and the responsible source path.
Do not silence a consistency error by expanding an exception list unless the exception represents a
confirmed design rule.

## Full site build

Use the ABINIT wrapper rather than invoking MkDocs directly:

```sh
python mksite.py build --strict
```

The wrapper:

1. Reads the ABINIT version information.
2. Generates top-level `mkdocs.yml` from `mkdocs.yml.in`.
3. Constructs the `Website` database.
4. Generates derived Markdown pages and expands underscore-prefixed topic templates.
5. Invokes MkDocs to build `site/`.
6. Adds custom `Website` warnings to its exit status.

Treat every warning from the custom website layer or MkDocs strict mode as a validation failure.
Resolve the original link, metadata, navigation, or generation problem rather than disabling strict
checking.

### Mutation warning

`mksite.py` writes `mkdocs.yml` before it processes the requested command, including `--help`.
A build also writes generated Markdown and the `site/` directory.
Inspect `git status` before and after running it, preserve unrelated user changes, and do not delete
pre-existing generated output unless the user authorizes cleanup.

Use `mkdocs.yml.in` as the version-controlled configuration source.
Do not make a lasting fix only in generated `mkdocs.yml`.

## HTML validation

After a successful site build, the optional HTML validator can inspect the entire site:

```sh
python mksite.py validate
```

It can also inspect selected pages:

```sh
python mksite.py validate site/index.html site/variables/index.html
```

This path requires the `py_w3c` validator dependency, which is not part of the core `mksite`
dependency group.
If it is unavailable, report that exact limitation; do not install it without authorization.

## Review rules

- Register every new documentation page in the `nav` section of `mkdocs.yml.in`.
- Edit underscore-prefixed topic templates such as `doc/topics/_BSE.md`, not only their generated
  counterparts.
- Use relative Markdown links with explicit `.md` targets and `index.md` for directory indexes.
- Prefer ABINIT wikilinks for variables, citations, topics, tests, tutorials, source files, and PDFs.
- Put page-specific figures in the page's dedicated assets directory.
- Use raw Python strings for variable descriptions containing LaTeX backslashes.
- Confirm whether an apparent missing page is generated before creating a duplicate static file.
- Preserve existing comments and generated-file warnings.

## Troubleshooting

### A page is reported as unreferenced

Check whether the page should be registered in `mkdocs.yml.in`, is generated from a template, or is an
intentional non-page asset.
Do not add arbitrary navigation entries merely to eliminate the warning.

### A relative link is unrecognized

Resolve it from the Markdown file containing the link.
Add the `.md` suffix for a page or `index.md` for a directory target.
Use MkDocs' suggested target only after confirming it represents the intended page.

### A wikilink cannot be resolved

Identify its namespace and check the corresponding database:

- Variables: `abimkdocs/variables_*.py`
- Citations: `abimkdocs/abiref.bib`
- Tests: registered files under `tests/`
- Topics and tutorials: `mkdocs.yml.in` and documentation sources
- Source files and PDFs: the referenced repository path

Fix the authoritative data rather than hardcoding generated HTML.

### Tests pass but the site build fails

The unit tests cover representative behavior, while the full build processes the complete navigation
and Markdown corpus.
Use the first MkDocs error or warning as the primary diagnostic; later messages may be consequences of
the same missing page, extension, or malformed metadata.

## Completion report

State which checks ran, their pass/fail counts, whether strict MkDocs completed, and whether optional
HTML validation was available.
List unresolved warnings with their source pages.
Mention generated or modified files left in the working tree, especially `mkdocs.yml`, generated topic
pages, and `site/`.
