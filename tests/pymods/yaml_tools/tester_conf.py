from __future__ import print_function, division, unicode_literals
import os
from yaml import YAMLError
from .conf_parser import conf_parser
from . import yaml_parse
from .errors import ConfigContextError
from .abinit_iterators import ITERATORS, IterStateFilter


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


def state_hash(d):
    st = []
    for it in ITERATORS:
        if it in d:
            st.append(it + str(d[it]))
    return hash(''.join(st))


class TesterConf:
    '''
        Interface to access parameters and constraints defined by the
        configuration file by following the traversal of the data tree.
    '''
    def __init__(self, src=None):
        self.known_params = conf_parser.parameters.copy()
        self.param_stack = []
        self.constraints_stack = []

        self.current_path = []

        # defaut conf is not supposed to use filters
        self.tree = conf_parser.make_trees(get_default_conf())[0]['__default']

        self.current_filter = None
        if src is not None:
            try:
                conf = yaml_parse(src)
            except YAMLError:
                conf = {}
            self.trees, self.filters = conf_parser.make_trees(conf)
            self.tree.update(self.trees['__default'])
            self.trees['__default'] = self.tree.copy()
        else:
            self.trees = {}
            self.filters = {}

        self.__tree_cache = {}

    @classmethod
    def from_file(cls, filename):
        '''
            Create a new instance of TesterConf from a configuration file.
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

    def get_constraints_for(self, obj):
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
            and wether or not it effectively has been defined)
        '''
        default = self.known_params[name]['__default']
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

    def use_filter(self, state):
        '''
            Start using filtered configurations if available.
        '''
        if state_hash(state) in self.__tree_cache:
            self.tree = self.__tree_cache[state_hash(state)]
        else:
            # Order filters from the most general to the most specific
            filters = sorted(
                [name for name, filt in self.filters if filt.match(state)],
                key=IterStateFilter.key
            )

            # Apply filtered trees, filters may be []
            for name in filters:
                self.tree.update(self.trees[name])

            self.__tree_cache[state_hash(state)] = self.tree

        self.will_enter = True
        return self

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

    def __enter__(self):
        '''
            Act as a context manager.
        '''
        # Should always use go_down or apply_filter when using 'with' block
        if not self.will_enter:
            raise ConfigContextError(self.current_path)
        self.will_enter = False
        return self

    def __exit__(self, type, value, traceback):
        '''
            Automatically go back when leaving with block.
        '''
        if not self.current_path:
            # already on top level, their is only filtered config that can
            # be cleaned
            self.tree = self.trees['__default'].copy()
        else:
            self.go_up()
