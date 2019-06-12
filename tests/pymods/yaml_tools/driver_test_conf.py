from __future__ import print_function, division, unicode_literals
import os
from copy import copy
from yaml import YAMLError
from .conf_parser import conf_parser
from . import yaml_parse
from .errors import ConfigContextError
from .abinit_iterators import ITERATORS


THIS_PATH = os.path.dirname(os.path.realpath(__file__))
DEFAULT_CONF_PATH = os.path.join(THIS_PATH, 'default_test.yaml')


def get_default_conf(filename):
    '''
        Load and parse default test config file.
    '''
    with open(filename) as f:
        try:
            return yaml_parse(f.read()) or {}
        except YAMLError:
            return {}


class DriverTestConf:
    '''
        Interface to access parameters and constraints defined by the
        configuration file by following the traversal of the data tree.
    '''
    default_conf = DEFAULT_CONF_PATH

    def __init__(self, src=None, metadata={}):
        self.known_params = conf_parser.parameters.copy()
        self.param_stack = []
        self.constraints_stack = []

        self.current_path = []
        self.current_state = {}

        self._infos = []

        # defaut conf is not supposed to use filters
        self.tree = conf_parser.make_trees(
            get_default_conf(self.default_conf),
            {'file name': 'default file'}
        )[0]['__default']

        self.current_filter = None
        if src is not None:
            try:
                conf = yaml_parse(src)
            except YAMLError as e:
                conf = {}
                self.warning('An error occured while parsing source:\n'
                             '{}: {}'.format(type(e).__name__, str(e)))
            self.trees, self.filters = conf_parser.make_trees(conf, metadata)
            self.tree.update(self.trees['__default'])
        else:
            self.info('No source have been provided apart from default'
                      ' config.')
            self.trees = {}
            self.filters = {}

        self.debug = False

        self.trees['__default'] = self.tree.copy()
        self._tree_cache = {}

    @classmethod
    def from_file(cls, filename):
        '''
            Create a new instance of DriverTestConf from a configuration file.
        '''
        with open(filename) as f:
            return cls(f.read(), {'file name': filename})

    def extra_info(self):
        return ['# ' + inf for inf in self._infos]

    def info(self, msg):
        self._infos.append('[INFO] ' + msg)

    def warning(self, msg):
        self._infos.append('[WARNING] ' + msg)

    @property
    def path(self):
        return tuple(self.current_path)

    def get_top_level_constraints(self):
        '''
            Return a list of the constraints defined at the tol level
            of configuration
        '''
        return self.tree.get_new_constraints_at(())

    def get_top_level_params(self):
        '''
            Return a dict of the parameters defined at the tol level
            of configuration
        '''
        return self.tree.get_new_params_at(())

    def get_constraints_for(self, obj):
        '''
            Return a list of the constraints in the current scope that apply
            to obj. If obj is None, return all available constraints in the
            scope.
        '''
        constraints = []
        exclude = set()
        already_defined = set()

        cursor = len(self.param_stack) - 1
        top = cursor  # top of the stack, bottom of the hierarchy

        def look_in(cons_dict, caller_lvl):
            for name, cons in cons_dict.items():
                if name in exclude or name in already_defined:
                    # if the constraint have been either already
                    # overridden or excluded only apply its exlusion
                    exclude.update(cons.exclude)
                elif (caller_lvl or cons.inherited) \
                        and (obj is None or cons.apply_to(obj)):
                    exclude.update(cons.exclude)
                    constraints.append(cons)
                    already_defined.add(name)

        while cursor >= 0:  # loop from bottom to top
            look_in(self.constraints_stack[cursor], cursor == top)
            cursor -= 1  # next level

        # finish with top level
        look_in(self.get_top_level_constraints(), False)

        return constraints

    def get_param(self, name):
        '''
            Return the value of the asked parameter as defined in
            the nearest scope or its default value (depending on
            wether or not it can be inherited from another scope
            and wether or not it effectively has been defined)
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

        top_params = self.get_top_level_params()
        if name in top_params:
            return top_params[name]
        else:
            return default

    def use_filter(self, state):
        '''
            Start using filtered configurations if available.
        '''
        def state_hash(d):
            st = []
            for it in ITERATORS:
                if it in d:
                    st.append(it + str(d[it]))
            return hash(''.join(st))

        self.current_state = state
        if state_hash(state) in self._tree_cache:
            self.tree = self._tree_cache[state_hash(state)]
        else:
            # Order filters from the most general to the most specific
            filters = sorted(
                ((filt, name) for name, filt in self.filters.items()
                 if filt.match(state)),
                reverse=True  # sort from the more specific/restrictive to the
                              # more general. In case of equality, sort using
                              # the reversed lexicographic order of the name
            )

            # Apply filtered trees, filters may be []
            for filt, name in filters:
                self.tree.update(self.trees[name])

            self._tree_cache[state_hash(state)] = self.tree

        # Rebuild stacks with the new tree
        self.rebuild_stacks()

        self.will_enter = True
        return self

    def clean_filter(self):
        '''
            Restore default filter state
        '''
        self.current_state = {}
        self.tree = self.trees['__default'].copy()

    def rebuild_stacks(self):
        '''
            Rebuild parameters and constraints stacks.
        '''
        path = copy(self.current_path)
        self.current_path = []
        self.param_stack = []
        self.constraints_stack = []
        for sp in path:
            self.go_down(sp)

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
        if self.current_path:  # if not already at top
            self.current_path.pop()
            self.param_stack.pop()
            self.constraints_stack.pop()

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
            self.clean_filter()
        else:
            self.go_up()
