from __future__ import print_function, division, unicode_literals
from inspect import isclass
from copy import deepcopy
from numpy import ndarray
from .errors import (UnknownParamError, ValueTypeError, InvalidNodeError,
                     AlreadySetKeyError, IllegalFilterNameError)
from .abinit_iterators import IterStateFilter


def empty_tree():
    return {
        'spec': {},
        'constraints': {},
        'parameters': {}
    }


def make_apply_to(type_):
    '''
        Return a function that take in argument the constraint and
        an object from the data tree an return True it the constraints apply to
        the object.
    '''
    if type_ == 'number':
        def apply_to(self, obj):
            return isinstance(obj, (int, float, complex))
    elif type_ == 'real' or type_ == float:
        def apply_to(self, obj):
            return isinstance(obj, float)

    elif type_ == 'integer' or type_ == int:
        def apply_to(self, obj):
            return isinstance(obj, int)

    elif type_ == 'complex' or type_ == complex:
        def apply_to(self, obj):
            return isinstance(obj, complex)

    elif type_ == 'Array' or type_ == ndarray:
        def apply_to(self, obj):
            return isinstance(obj, ndarray)

    elif type_ == 'this':
        def apply_to(self, obj):
            return True

    elif isclass(type_):
        def apply_to(self, obj):
            return isinstance(obj, type_)

    else:  # dummy
        def apply_to(self, obj):
            return False

    return apply_to


class Constraint(object):
    '''
        Represent a constraint to be applied to some piece of data.
    '''
    def __init__(self, name, test, val_type, inherited, use_params, exclude,
                 apply_to):
        self.name = name
        self.test = test
        self.type = val_type
        self.inherited = inherited
        self.use_params = use_params
        self.exclude = exclude
        self.value = None

        self.__apply_to = make_apply_to(apply_to)

    def check(self, ref, tested, conf):
        '''
            Return True if the constraint is verified.
        '''
        params = [conf.get_param(p) for p in self.use_params]
        return self.test(self.value, ref, tested, *params)

    def apply_to(self, obj):
        '''
            __apply_to is not a method so we have to pass self explicitly
        '''
        return self.__apply_to(self, obj)

    def copy(self):
        '''
            Create a copy of self.
        '''
        cp = Constraint(self.name, self.test, self.type, self.inherited,
                        self.use_params, self.exclude, 'dummy')
        cp.__apply_to = self.__apply_to
        return cp

    def with_value(self, val):
        '''
            Create a copy of self with the value attribute set.
        '''
        if not isinstance(val, self.type):
            raise ValueTypeError(self.name, self.type, val)
        cp = self.copy()
        cp.value = val
        return cp


class ConfTree(object):
    '''
        Configuration tree accessor. The actual tree is a dictionary.
    '''
    def __init__(self, dict_tree):
        self.dict = dict_tree

    def add(self, key, d):
        if key in self.dict:
            raise AlreadySetKeyError(key)
        self.dict[key] = d

    def copy(self):
        return ConfTree(deepcopy(self.dict))

    def update(self, tree):
        '''
            Update self with values found in tree.
        '''
        self.__update(self.dict, tree.dict)

    def __update(self, old_d, new_d):
        '''
            Recursively update the content of old_d with new_d.
        '''
        for key, cons in new_d['constraints'].items():
            old_d['constraints'][key] = cons

        for param, value in new_d['parameters'].items():
            old_d['parameters'][param] = value

        for spec, spec_d in new_d['spec'].items():
            if spec not in old_d['spec']:
                old_d['spec'][spec] = empty_tree()
            self.__update(old_d['spec'][spec], spec_d)

    def get_spec_at(self, path):
        '''
            Get specializations defined at a given node in the tree.
            Return an empty dictionary if the path does not exists.
        '''
        d = self.dict
        for spec in path:
            if spec in d['spec']:
                d = d['spec'][spec]
            else:
                return {}
        return d['spec']

    def get_new_params_at(self, path):
        '''
            Get params defined at a given node in the tree.
            Return an empty dictionary if the path does not exists.
        '''
        d = self.dict
        for spec in path:
            if spec in d['spec']:
                d = d['spec'][spec]
            else:
                return {}
        return d['parameters']

    def get_new_constraints_at(self, path):
        '''
            Get constraints defined at a given node in the tree.
            Return an empty dictionary if the path does not exists.
        '''
        d = self.dict
        for spec in path:
            if spec in d['spec']:
                d = d['spec'][spec]
            else:
                return {}
        return d['constraints']


class ConfParser(object):
    '''
        Test configuration loader and parser. It take output from yaml parser
        and build the actual configuration tree.
    '''
    def __init__(self):
        self.parameters = {}
        self.constraints = {}

    def import_parser(self, conf_parser):
        '''
            Import constraints and parameters from another conf_parser.
            In case of name conflict the local values are overridden by the new
            values.
        '''
        self.parameters.update(conf_parser.parameters)
        self.constraints.update(conf_parser.constraints)

    def parameter(self, token, default=None, value_type=float, inherited=True):
        '''
            Register a parameter to be recognised while parsing config.
        '''
        self.parameters[token] = {
            'type': value_type,
            'inherited': inherited,
            'default': default,
        }

    def constraint(self, name=None, value_type=float, inherited=True,
                   apply_to='number', use_params=[], exclude=set()):
        '''
            Register a constraints to be recognised while parsing config.
        '''
        def register(fun):
            if name is None:
                nname = fun.__name__
            else:
                nname = name
            for param in use_params:
                if param not in self.parameters:
                    raise UnknownParamError(nname, param)
            self.constraints[nname] = Constraint(
                nname, fun, value_type, inherited, use_params,
                set(exclude), apply_to
            )
            return fun
        return register

    def __make_tree(self, parsed_src):
        '''
            Recursively build the configuration tree
        '''
        tree = empty_tree()
        if parsed_src:
            for key, val in parsed_src.items():
                if key in self.constraints:
                    tree['constraints'][key] = \
                            self.constraints[key].with_value(val)
                elif key in self.parameters:
                    if not isinstance(val, self.parameters[key]['type']):
                        raise ValueTypeError(key, self.parameters[key]['type'],
                                             val)
                    tree['parameters'][key] = val
                elif isinstance(val, dict):
                    tree['spec'][key] = self.__make_tree(val)
                else:
                    raise InvalidNodeError(key, val)
        return tree

    def make_trees(self, parsed_src):
        '''
            Create a ConfTree instance from the yaml parser output.
        '''
        filters = {}
        trees = {}

        if 'filters' in parsed_src:
            # Get filters definitions
            for name, filt in parsed_src['filters'].items():
                if name == '__default' or name in self.constraints \
                   or name in self.parameters:
                    raise IllegalFilterNameError(name)
                filters[name] = IterStateFilter(filt)

                # Parse each filtered tree then remove it from the source tree
                if name in parsed_src:
                    trees[name] = ConfTree(self.__make_tree(parsed_src[name]))
                    del parsed_src[name]
                else:
                    pass  # Should we raise an error ?

            # Remove the filters field from the dict
            del parsed_src['filters']

        # Parse the fields remaining as the default tree
        trees['__default'] = ConfTree(self.__make_tree(parsed_src))

        return trees, filters
