# Documentation Guidelines for AI Agents

This document provides instructions for automated agents, AI assistants, and bots that contribute to or modify the Abinit documentation. Please adhere to these guidelines when generating, formatting, or updating content.

## General Approach

The Abinit documentation is built using [MkDocs](https://www.mkdocs.org/) with the [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/) theme. All documentation is written in Markdown and must be compatible with the extensions configured in `mkdocs.yml.in`.

- **Markdown-first**: Documentation is written in plain text Markdown.
- **Directory Structure**: Markdown files are located inside the `doc/` directory (e.g., `doc/tutorial`, `doc/topics`, `doc/theory`, etc.). All new pages must be registered in the `nav` section of `mkdocs.yml.in`.
- **Previewing**: The site is generated using a customized python script. Running `./mksite.py serve` starts a local server at `http://127.0.0.1:8000/`.

## Writing and Formatting Rules

### 1. File Structure and Front Matter
- Include YAML front matter at the top of new `.md` files (e.g., `title`, `authors`).
- Paragraph names should be concise for navigation.
- Start pages with an introductory paragraph, using one `#` (H1) heading for the title, followed by `##` (H2) and `###` (H3) as needed.

### 2. Variables and Python Documentation
- Documentation for Abinit input variables is stored in python files (e.g., `..//abimkdocs/variables_abinit.py`) as lists of dictionaries.
- **Important**: When writing variable descriptions in these python files, always use Python raw strings (`r"""..."""`) to prevent escape sequence errors when using LaTeX (e.g., `\alpha`).

### 3. Links and Cross-Referencing
- **Root-Relative URLs**: Use relative URLs for internal markdown links (e.g., `[MBPT document](../theory/mbt.md)`) instead of root-relative paths, as this is the recommended practice to avoid link breakage.
- **Wikilinks**: Abinit uses a custom wikilink extension. Use this syntax extensively:
  - Variables: `[[ecut]]` for Abinit variables or `[[dipdip@anaddb]]` if dipdip is an anaddb variable
  - Citations: `[[cite:Amadon2008]]`
  - Topics: `[[topic:BSE]]`
  - Tests: `[[test:libxc_41]]`
  - Tutorials: `[[tutorial:gw1]]`
  - Source files: `[[src:94_scfcv/m_scfcv.F90]]`

### 4. Wikilink Internals and Link Validation
- **Python Generation**: The custom wikilink extension is backed by `abimkdocs/website.py` (specifically the `get_wikilink` function). It dynamically parses tokens and generates compliant HTML links. It uses an internal toggle (`self.use_relative_urls = True`) to enforce standard relative URLs.
- **Directory vs File Links**: When generating or hardcoding markdown links, links pointing to directories *must* include `index.md` (e.g., `../tutorial/index.md`), and links pointing to specific files *must* include the `.md` extension (e.g., `../theory/mbt.md`). Failure to include `.md` will cause MkDocs to reject the link as unrecognized.
- **Testing Links**: Always run `python3 ./mksite.py build` (or `./mksite.py build`) to trigger MkDocs' strict link validator. Scan the output for `INFO` or `WARNING` messages (e.g., "contains an unrecognized relative link"). MkDocs often provides the correct resolution string ("Did you mean '...'?"), which should be used to fix the target links directly in the Python macros or static Markdown files.

### 5. Topic Templates
- Files in `doc/topics/` that start with an underscore (e.g., `_AbiPy.md`) are template files. They are automatically processed by `../mksite.py` to fill in `{{ related_variables }}` and `{{ selected_input_files }}`.
Agents must modify the template file with the underscore, not the auto-generated one.

### 6. Math and Equations
- LaTeX is fully supported via MathJax.
- Use `$...$` for inline equations and `$$...$$` or `\begin{equation}...\end{equation}` for block equations.
- In LaTeX, standard macros are available (e.g., `\rr`, `\GG`, `\kk`, `\qq`, `\kq`). See `../abinit_theme/main.html` for details.

### 7. Extensions
- **Admonitions**: Emphasize important text using MkDocs Material admonitions (`!!! note`, `!!! warning`, `!!! tip`, `!!! danger`).
- **Collapsible Blocks**: Use the Details extension for long content (`??? note "Title"` or `???+ note "Title"` for an initially open block).
- **Figures**: Place images in a dedicated `assets` directory named after the markdown file (e.g., `doc/tutorial/bse_assets/` for `doc/tutorial/bse.md`).
- **Videos & PDFs**: Use `[[pdf:filename]]` to link to internal PDFs.

By adhering to these rules, you will ensure that the documentation remains consistent, correctly processed by the custom Abinit MkDocs pipeline, and seamlessly integrated with the rest of the project.
