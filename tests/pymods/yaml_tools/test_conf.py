from __future__ import print_function, division, unicode_literals
import os
from yaml import YAMLError
from .conf_parser import conf_parser
from . import yaml_parse
from .errors import ConfigContextError


def get_default_conf():
    '''
        Load and parse default test config file.
    '''
    dir_path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(dir_path, 'default_test.yaml')
    with open(file_path) as f:
        try:
            return yaml_parse(f.read())
        except YAMLError:
            return {}


class TestConf:
    '''
        Interface to access parameters and constraints defined by the
        configuration file by following the traversal of the data tree.
    '''
    def __init__(self, src):
        self.known_params = conf_parser.parameters.copy()
        self.param_stack = []
        self.constraints_stack = []

        self.current_path = []

        self.tree = conf_parser.make_tree(get_default_conf())
        try:
            conf = yaml_parse(src)
        except YAMLError:
            conf = {}

        new_tree = conf_parser.make_tree(conf)
        self.tree.update(new_tree)

    @classmethod
    def from_file(cls, filename):
        '''
            Create a new instance of TestConf from a configuration file.
        '''
        with open(filename) as f:
            return cls(f.read())

    @property
    def path(self):
        return tuple(self.current_path)

    def get_top_level_constraints(self):
        '''
            Return a list of the constraints defined at the tol level
            of configuration
        '''
        return self.tree.get_new_constraints_at(())

    def get_constraint_for(self, obj):
        '''
            Return a list of the constraints in the current scope that
            apply to obj.
        '''
        constraints = []
        exclude = set()

        cursor = len(self.param_stack) - 1
        top = cursor

        while cursor >= 0:
            for name, cons in self.constraints_stack[cursor].items():
                if cursor < top and not cons.inherited:
                    # pass if we are no longer on
                    # the same level as the caller
                    continue
                elif name in exclude:
                    continue
                else:
                    if cons.apply_to(obj):
                        exclude.update(cons.exclude)
                        constraints.append(cons)
            cursor -= 1

        return constraints

    def get_param(self, name):
        '''
            Return the value of the asked parameter as defined in
            the nearest scope or its default value (depending on
            wether or not it can be inherited from another scope
            and wether or not it effectvely has been defined)
        '''
        default = self.known_params[name]['default']
        cursor = len(self.param_stack) - 1

        # browse scope from deeper to the top until param is
        # define
        while cursor >= 0:
            if name in self.param_stack[cursor]:
                return self.param_stack[cursor][name]
            elif not self.known_params[name]['inherited']:
                return default
            else:
                cursor -= 1

        return default

    def go_down(self, child):
        '''
            Go deeper in the tree.
        '''
        # Append the new level to the path
        self.current_path.append(child)

        # Append the newly defined parameters and constraints to the scopes
        self.param_stack.append(
            self.tree.get_new_params_at(self.current_path))
        self.constraints_stack.append(
            self.tree.get_new_constraints_at(self.current_path))

        self.will_enter = True
        return self

    def __enter__(self):
        '''
            Act as a context manager.
        '''
        if not self.will_enter:  # Should always use go_down when using with
            raise ConfigContextError(self.current_path)
        self.will_enter = False
        return self

    def __exit__(self, type, value, traceback):
        '''
            Automatically go back when leaving with block.
        '''
        self.go_up()
        # FIXME Is there error cases (type, value, traceback)
        # where we should stop error from propating ?

    def go_up(self):
        '''
            Go back to a higher level of the tree.
        '''
        try:
            self.current_path.pop()
            self.param_stack.pop()
            self.constraints_stack.pop()
        except IndexError:
            pass
