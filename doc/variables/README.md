Abinit vars in YAML (2014)
==========================

Y. Gillet, F. Abreu Araujo and X. Gonze

The documentation about the input variables is contained in the YAML file `initial_files/abinit_vars.yml`.
See the README.md there.

## Generating the new html

In the ~abinit/doc directory, run:

    python generate_doc.py

The new html files will be updated in ~abinit/doc/input_variables/generated_files
Do not modify files in the latter directory. It is useless since they are automatically generated.
