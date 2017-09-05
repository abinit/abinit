Abinit topics
=============

F. Jollet and X. Gonze

The information to generate the ABINIT topics is contained in the YAML files in ~abinit/doc/topics/initial_files

## Generating the new html

In the ~abinit/doc directory, run:

    python generate_doc.py

The new html files will be updated in ~abinit/doc/topics/generated_files .
Do not modify files in the latter directory. It is useless since they are automatically generated.
