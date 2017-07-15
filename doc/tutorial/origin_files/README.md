Abinit tutorial
===============

X. Gonze

The information to generate the lessons of the ABINIT tutorial is contained in the YAML files in ~abinit/doc/tutorial/initial_files .
To add a new lesson, you must modify the "lessons.yml" file, adding some structuring information on the lesson,
and add a new html file that gives the intro and the content for the new lesson. 
The defaults for each lesson are contained in lessons.yml, see the last part of this file,
and are possibly superceded by the specified sections in the new YAML file.

## Generating the new html

In the ~abinit/doc directory, run:

    python generate_doc.py

The new html files will be updated in ~abinit/doc/tutorial/generated_files .
Do not modify files in the latter directory. It is useless since they are automatically generated.
