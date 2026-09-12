# ABINIT's MkDocs theme overrides

This directory contains the Jinja templates that customize the Material for MkDocs theme used by the ABINIT documentation website.
It is the presentation layer of a larger documentation pipeline implemented by `mksite.py`, the `abimkdocs` Python package, and MkDocs.

For instructions on installing the documentation dependencies, writing pages, and building or serving the website, see the [ABINIT documentation guide](../doc/developers/abimkdocs.md).
The [ABINIT Markdown guide](../doc/developers/markdown.md) documents the custom link and authoring syntax.
The implementation of variable declarations and `website.py` is summarized in the [`abimkdocs` README](../abimkdocs/README.md).

## Role in the documentation stack

The website is built in three main stages:

```text
Repository sources
  Markdown, variable declarations, tests, bibliography, version
                              |
                              v
                       mksite.py
  generates mkdocs.yml and asks Website to generate Markdown
                              |
                              v
                abimkdocs + MkDocs
  preprocess aliases/includes, resolve wikilinks, convert Markdown
                              |
                              v
          Material theme + abinit_theme overrides
  render navigation, page content, footer, and ABINIT UI additions
                              |
                              v
                         site/*.html
```

`abinit_theme` does not generate variable pages, inspect tests, or resolve wikilinks.
Those operations belong to `abimkdocs.website.Website` and the ABINIT Markdown extensions.
The theme receives the resulting page content and metadata from MkDocs and controls how they are presented in HTML.

## How MkDocs selects the theme

The authoritative configuration is the top-level `mkdocs.yml.in` template.
Its theme section contains:

```yaml
theme:
    name: material
    custom_dir: abinit_theme
```

`name: material` selects Material for MkDocs as the base theme.
`custom_dir: abinit_theme` tells MkDocs to look in this directory before falling back to the installed Material templates.
A file here therefore overrides the Material template with the same relative path, while every template not present here continues to come from Material.

Do not copy the entire Material theme into this directory.
Keep overrides small so that ABINIT can continue to benefit from upstream theme fixes and so that differences remain reviewable.

`mksite.py` reads `mkdocs.yml.in`, replaces the `ABINIT_VERSION` placeholder, and writes the generated `mkdocs.yml` used by the MkDocs command.
Edit `mkdocs.yml.in`, not the generated `mkdocs.yml`, when changing theme configuration.

## Files in this directory

| Path | Responsibility |
|---|---|
| `main.html` | Extends Material's `base.html`, adds external libraries and conditional assets, and appends the floating navigation control. |
| `partials/footer.html` | Replaces Material's footer to show ABINIT version and page authors while preserving navigation and social links. |

Both files use the Jinja context supplied by MkDocs and Material, including `page`, `page.meta`, `nav`, `config`, and `mkdocs_version`.
Consult the version of Material pinned by the project before using a new block or context variable because upstream template interfaces can change.

## `main.html`

`main.html` extends Material's `base.html` and overrides two blocks.

### Library block

The `libs` block calls `{{ super() }}` first, preserving the scripts and styles loaded by Material.
It then adds libraries required by ABINIT-specific HTML generated elsewhere in the documentation infrastructure:

- jQuery and jQuery UI for dialogs and older interactive components;
- jQuery Modal;
- Popper and Tippy.js for popovers;
- Font Awesome icons used by the floating action button;
- Pace for page-loading feedback;
- MathJax configuration and ABINIT mathematical macros;
- Inconsolata for input-variable presentation.

Some assets are enabled through page front matter:

```yaml
---
plotly: true
light_gallery: true
---
```

When `page.meta.plotly` is true, the template loads Plotly.
When `page.meta.light_gallery` is true, it loads the LightGallery assets.
Keep expensive page-specific libraries behind metadata checks instead of loading them globally.

The template deliberately avoids Bootstrap because its CSS and JavaScript conventions conflict with Material for MkDocs.
Before adding a frontend dependency, check its selectors and global objects against both Material and the existing libraries.

### Content block

The `content` block calls `{{ super() }}` to render the normal Material page and then appends a floating navigation button.
The button provides a return-to-top action and a compact menu derived from the top-level MkDocs navigation.

The corresponding behavior is implemented in `doc/extra_javascript/abidocs.js`, and its styling is primarily in `doc/css/extra.css`.
Changing the generated HTML structure may therefore require coordinated template, JavaScript, and CSS changes.

## Footer override

`partials/footer.html` replaces Material's footer partial.
It retains previous/next-page navigation and the Material social-links partial.
It adds ABINIT-specific information from the MkDocs context:

- `config.extra.version`, populated from the `ABINIT_VERSION` substitution performed by `mksite.py`;
- `config.copyright`, defined in `mkdocs.yml.in`;
- `page.meta.authors`, read from a page's YAML front matter.

A source page can provide authors with:

```yaml
---
authors: MG, XG
---
```

The footer obtains this metadata only after MkDocs and the ABINIT plugin have processed the page.
The theme should display metadata, not attempt to parse source Markdown itself.

## Interaction with `mksite.py`

The normal entry point is:

```sh
./mksite.py build
./mksite.py serve
```

For `build`, `serve`, and `gh-deploy`, `mksite.py` performs these operations before delegating to the standard MkDocs command-line interface:

1. Detects the ABINIT version.
2. Generates `mkdocs.yml` from `mkdocs.yml.in`.
3. Constructs the `abimkdocs.website.Website` singleton with `doc` as its root.
4. Generates variable pages, topic pages, bibliography entries, test documentation, and other dynamic Markdown.
5. Invokes MkDocs, which loads the plugin, Markdown extensions, Material theme, and these template overrides.

The generated `mkdocs.yml` provides the connection to this directory through `theme.custom_dir`.
It also declares the CSS, JavaScript, plugins, and Markdown extensions on which the templates and generated markup depend.

## Interaction with the `abimkdocs` package

The `abimkdocs` package supplies the semantic layer used before theme rendering.
In particular:

- `Website.generate_markdown_files()` creates variable indexes, variable pages, topic pages, bibliography pages, and test-suite pages.
- `abimkdocs.preprocessor` expands aliases, includes, and ABINIT macros.
- `abimkdocs.wikilinks` resolves constructs such as `[[ecut]]`, `[[anaddb:asr]]`, and `[[cite:Gonze2009]]`.
- `abimkdocs.mdx_figcaption` provides figure-caption support.
- `Website.get_wikilink()` emits page-relative links used by the Markdown extension.

These extensions are listed in `mkdocs.yml.in` under `markdown_extensions`.
They access the `Website` singleton created by `mksite.py`, so invoking MkDocs directly can bypass required initialization and generated content.
Use `mksite.py` as the normal website entry point.

## Interaction with the `abimkdocs` MkDocs plugin

The separate `abimkdocs_plugin` package registers an MkDocs plugin named `abimkdocs`.
It is enabled by this entry in `mkdocs.yml.in`:

```yaml
plugins:
    - abimkdocs
    - search
```

During MkDocs' `on_page_markdown` event, the plugin reconstructs the page's YAML front matter and adds `rpath`, the page path relative to the documentation root.
The reconstructed metadata is then visible to the ABINIT Markdown processing stage.
Later, MkDocs exposes page metadata to Jinja as `page.meta`, allowing the theme to read fields such as `authors`, `plotly`, and `light_gallery`.

The flow for metadata is therefore:

```text
Markdown front matter
       |
       v
MkDocs Page.meta
       |
       +-- abimkdocs plugin adds rpath for Markdown processing
       |
       `-- Jinja page.meta controls theme rendering
```

Do not move semantic Markdown transformations into Jinja templates.
Do not generate visual page structure in the plugin when it belongs in the theme.
Keeping these responsibilities separate makes both Markdown processing and HTML rendering testable.

## CSS and JavaScript

The custom static assets do not live in `abinit_theme`.
They are stored below the MkDocs `docs_dir` and registered in `mkdocs.yml.in`:

```yaml
extra_css:
    - css/extra.css
    - css/gallery.css
    - css/my_codehilite.css

extra_javascript:
    - extra_javascript/abidocs.js
```

This arrangement lets MkDocs copy and serve the assets as documentation content while the Jinja templates concentrate on page structure.

`abidocs.js` initializes Tippy popovers, controls the floating button, implements the input-variable search tabs, and opens dialogs generated by `website.py`.
The HTML produced by `Website.build_varsearch_html()` and `Website.dialog_from_filename()` uses IDs and classes expected by this JavaScript and by `extra.css`.

When changing an interactive component, trace all three parts:

| Layer | Typical location |
|---|---|
| Generated markup | `abimkdocs/website.py` or source Markdown. |
| Behavior | `doc/extra_javascript/abidocs.js`. |
| Presentation | `doc/css/extra.css` and, when necessary, `abinit_theme/main.html`. |

## Making theme changes

Use these guidelines when editing the templates:

1. Override the smallest possible Material block or partial.
2. Call `super()` when extending a block whose upstream contents must remain.
3. Use MkDocs or Material URL filters for internal resources and navigation links.
4. Guard page-specific resources with front-matter flags.
5. Keep content generation and link resolution in `abimkdocs`, not in templates.
6. Check associated CSS and JavaScript whenever IDs, classes, or element structure change.
7. Verify that the template still works with pages that have no optional metadata.

The base theme is pinned in `pyproject.toml` and `requirements.txt`.
Review the Material release notes and upstream template changes before upgrading it, especially because the footer is a full partial override.

## Testing and previewing

Activate the documentation Python environment and run the focused infrastructure tests:

```sh
pytest abimkdocs_tests/test_website.py
```

Build the complete website with strict diagnostics when possible:

```sh
./mksite.py build --strict
```

For visual development, use:

```sh
./mksite.py serve
```

Inspect at least:

- an ordinary documentation page;
- a generated variable page;
- a page with authors in its front matter;
- a page using Plotly or LightGallery metadata if those paths changed;
- previous/next links, social links, the floating menu, popovers, and dialogs;
- narrow and wide viewport layouts.

Template correctness cannot be established fully by Python unit tests.
A local site build and browser inspection are required for changes that affect markup, styling, or behavior.

## Distribution through the Autotools build system

`abinit_theme` is registered as documentation data in `config/specs/buildsys.conf`.
After adding, renaming, or removing a file in this directory, run:

```sh
./config/scripts/makemake
```

The generator refreshes `config/dist/auto-abinit_theme.lst`.
That generated list is consumed by the top-level Automake files so `make dist` includes the custom theme in ABINIT source archives.

For a distribution-sensitive change, verify with:

```sh
make dist
```

The generated site under `site/` is not source material and must not be committed.
