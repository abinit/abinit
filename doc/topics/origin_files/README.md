Abinit topics
=============

F. Jollet and X. Gonze

The information to generate the ABINIT topics is contained in the YAML files in ~abinit/doc/topics/initial_files
To add a new topic, you must modify the "list_of_topics.yml" file, and add a new YAML file that gives the 
information for the new topics. The defaults for each section are contained in default_topics.yml
and are possibly superceded by the specified sections in the new YAML file.

## Generating the new html

In the ~abinit/doc directory, run:

    python generate_doc.py

The new html files will be updated in ~abinit/doc/topics/generated_files .
Do not modify files in the latter directory. It is useless since they are automatically generated.
