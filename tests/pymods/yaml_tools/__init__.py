'''
    This package gather all tools used by the Abinit test suite for
    manipulating YAML formated data.
'''
from __future__ import print_function, division, unicode_literals
try:
    import yaml
    import numpy  # numpy is also required

    is_available = True

except ImportError:
    is_available = False

if is_available:

    class CorruptedDocument(object):
        def __init__(self, error, context=''):
            self.error = error
            self.context = context

        def __repr__(self):
            return "<Corrupted document {}: {}>".format(
                self.error.__class__.__name__, str(self.error))

    def yaml_parse(content, *args, **kwargs):
        from . import structures
        try:
            return yaml.load(content, *args, **kwargs)
        except yaml.YAMLError as e:
            return CorruptedDocument(e, context=content)

    yaml_print = yaml.dump
