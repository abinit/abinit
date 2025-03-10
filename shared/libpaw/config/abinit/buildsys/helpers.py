import os

import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


class LibraryBuilder(object):


    def __init__(self, cfg_path):

        with open(yaml_path, "r") as yaml_file:
            self.specs = load(yaml_file, Loader=Loader)

        if ( not "entity" in self.specs.keys() ):
            raise KeyError("Missing entity description in {}".format(cfg_path))
        if ( self.specs["entity"] != "srcdir" ):
            raise ValueError("Invalid entity description in {}".format(cfg_path))


    def get_automake(self):
