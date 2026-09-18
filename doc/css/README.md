# Documentation stylesheets

This directory contains the stylesheets loaded by the MkDocs site through `extra_css` in `mkdocs.yml.in`.

- `extra.css` contains project-wide theme overrides, typography, tables, admonitions, and interactive widget styles.
- `gallery.css` contains the layout and hover effects used by image galleries.
- `my_codehilite.css` contains syntax-highlighting styles for code blocks.

## Maintenance

Keep general site and theme changes in `extra.css`.
Use the more specialized files only for gallery layout or syntax highlighting.
When changing selectors, check their use in the documentation templates and JavaScript before removing or renaming them.

Run `./mksite.py build` from the repository root to validate stylesheet changes with the complete documentation site.
