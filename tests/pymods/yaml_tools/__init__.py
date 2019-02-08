'''
    This package gather all tools used by the Abinit test suite for manipulating YAML formated data.
'''
from __future__ import print_function, division, unicode_literals
try:
    import yaml

    is_available = True

except ImportError:
    is_available = False

if is_available:

    class CorruptedDocument(object):
        def __init__(self, error):
            self.error = error

        def why(self):
            return self.error.message

    def yaml_parse(*args, **kwargs):
        from . import structures
        try:
            return yaml.load(*args, **kwargs)
        except Exception as e:
            return CorruptedDocument(e)

    yaml_print = yaml.dump
