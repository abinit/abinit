Abinit vars in YAML (2014)
==========================

Y. Gillet, F. Abreu Araujo and X. Gonze

The documentation about the input variables is contained in the YAML file `abinit_vars.yml`.
To edit the file abinit_vars.yml, it is recommended to use the Abivars Java Previewer
with the command:

    java -jar Abivars.jar

It permits an easy graphical view of the different variables and links between them.

This application opens the file `abinit_vars.yml` by default.
When you are done, you can overwrite `abinit_vars.yml` file from the `File -> Save` menu.
You can also open any other file from the `File -> Open` menu.

## Generating the new html

Simply run:

    python abi_yml2html.py

The new html files will be updated in `html_automatically_generated`
