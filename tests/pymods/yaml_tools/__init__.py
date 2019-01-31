'''
    This package gather all tools used by the Abinit test suite for manipulating YAML formated data.
'''
from __future__ import print_function, division, unicode_literals
try:
    import yaml

    yaml_parse = yaml.load
    yaml_print = yaml.dump

    is_available = True

except ImportError:
    is_available = False
