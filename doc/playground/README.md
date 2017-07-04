Automatic generation of documentation for ABINIT : input variables and topics.
==============================================================================

Y. Gillet, F. Abreu Araujo, F. Jollet and X. Gonze (2014-2017)

The documentation about the input variables is contained in the directory yml_files.
In particular, the YAML file `abinit_vars.yml` (for the input variables),
and the different files topic*.yml .

To edit the file abinit_vars.yml, it is recommended to use the Abivars Java Previewer
with the command (issued in the directory yml_files):

    java -jar Abivars.jar

It permits an easy graphical view of the different variables and links between them.

This application opens the file `abinit_vars.yml` by default.
When you are done, you can overwrite `abinit_vars.yml` file from the `File -> Save` menu.
You can also open any other file from the `File -> Open` menu.

## Generating the new html files (the different files that describe the input variables, 
## as well as the topic* files

Simply run:

    python abi_yml2html.py

The new html files will be updated in `html_automatically_generated`
