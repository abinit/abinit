
try:
    from yaml import CDumper as Dumper
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader


class LibraryBuilder:


    def __init__(self, cfg_path):

        with open(yaml_path) as yaml_file:
            self.specs = load(yaml_file, Loader=Loader)

        if ( "entity" not in self.specs.keys() ):
            raise KeyError(f"Missing entity description in {cfg_path}")
        if ( self.specs["entity"] != "srcdir" ):
            raise ValueError(f"Invalid entity description in {cfg_path}")


    #def get_automake(self):
