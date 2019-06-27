from __future__ import print_function, division, unicode_literals
from inspect import isclass
from copy import deepcopy
from .errors import (UnknownParamError, ValueTypeError, InvalidNodeError,
                     IllegalFilterNameError)
from .abinit_iterators import IterStateFilter
from .tricks import cstm_isinstance
from .common import Undef, normalize_attr, string, BaseArray, FailDetail
from warnings import warn


def make_apply_to(type_):
    '''
        Return a function that takes in argument the constraint and
        an object from the data tree and returns True it the constraints
        apply to the object.
    '''
    if type_ == 'number':
        def apply_to(self, obj):
            return isinstance(obj, (int, float, complex))

    elif type_ == 'real':
        def apply_to(self, obj):
            return isinstance(obj, float)

    elif type_ == 'integer':
        def apply_to(self, obj):
            return isinstance(obj, int)

    elif type_ == 'complex':
        def apply_to(self, obj):
            return isinstance(obj, (float, complex))

    elif type_ == 'Array':
        def apply_to(self, obj):
            return getattr(obj, '_is_base_array', False)

    elif type_ == 'this':
        def apply_to(self, obj):
            return True

    elif isclass(type_):
        def apply_to(self, obj):
            return cstm_isinstance(obj, type_)

    else:  # dummy
        def apply_to(self, obj):
            return False

    return apply_to


class Constraint(object):
    '''
        Represent a constraint to be applied to some piece of data.
    '''
    def __init__(self, name, test, val_type, inherited, use_params, exclude,
                 apply_to, handle_undef, value=None, metadata={}):
        self.name = name
        self.test = test
        self.type = val_type
        self.inherited = inherited
        self.use_params = use_params
        self.exclude = exclude
        self.value = value
        self.metadata = metadata
        self.handle_undef = handle_undef

        self._apply_to_type = apply_to
        self._apply_to = make_apply_to(apply_to)

    def __repr__(self):
        return 'Constraint({})'.format(', '.join((
            str(self.name),
            str(self.test),
            str(self.type),
            str(self.inherited),
            str(self.use_params),
            str(self.exclude),
            str(self.value),
            str(self.metadata),
        )))

    def check(self, ref, tested, conf):
        '''
            Return True if the constraint is verified.
        '''
        # apply to floats at least
        if getattr(ref, '_not_available', False):
            return FailDetail(
                'This constraint was to be applied to a document that is not'
                ' available (check warnings).'
            )
        if self.handle_undef:
            if isinstance(ref, (float, complex)) and self._apply_to(self, 1.0):
                if conf.get_param('allow_undef'):
                    if Undef.is_undef(ref):
                        return True
                elif Undef.is_undef(ref) or Undef.is_undef(tested):
                    return FailDetail('undef value have been found.')
            elif (getattr(ref, '_is_base_array', False)
                  and self._apply_to(self, BaseArray((0,)))):
                if conf.get_param('allow_undef'):
                    if ref._has_undef:
                        return True
                elif ref._has_undef or tested._has_undef:
                    return FailDetail('undef value have been found.')

        params = [conf.get_param(p) for p in self.use_params]
        return self.test(self.value, ref, tested, *params)

    def apply_to(self, obj):
        '''
            _apply_to is not a method so we have to pass self explicitly
        '''
        return self._apply_to(self, obj)

    def copy(self):
        '''
            Create a copy of self.
        '''
        cp = Constraint(self.name, self.test, self.type, self.inherited,
                        self.use_params, self.exclude, 'copying',
                        self.handle_undef)
        cp._apply_to = self._apply_to
        return cp

    def with_value(self, val, metadata={}):
        '''
            Create a copy of self with the value attribute set.
        '''
        if not isinstance(val, self.type):
            raise ValueTypeError(self.name, self.type, val)
        cp = self.copy()
        cp.value = val
        cp.metadata = metadata
        return cp

    def __eq__(self, other):
        return (
            isinstance(other, Constraint)
            and self.name == other.name
            and self.type == other.type
            and self.inherited == other.inherited
            and self.use_params == other.use_params
            and self.exclude == other.exclude
            and self.value == other.value
            and self.handle_undef == other.handle_undef
        )

    def __ne__(self, other):
        return not (self == other)


class SpecKey(object):
    def __init__(self, name, hardreset=False):
        self.name = normalize_attr(name)
        self.hardreset = hardreset

    @classmethod
    def parse(cls, name):
        hardr = False
        if isinstance(name, int):
            name = string(name)
        elif name.endswith('!'):
            hardr = True
            name = name[:-1]

        return cls(name, hardreset=hardr)

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return isinstance(other, SpecKey) and self.name == other.name

    def __ne__(self, other):
        return not isinstance(other, SpecKey) or self.name != other.name

    def __repr__(self):
        return '"' + self.name + ('!' if self.hardreset else '') + '"'


class ConfTree(object):
    '''
        Configuration tree wrapper. Give access to constraints and parameters
        defined in the nodes.
    '''
    def __init__(self, dict_tree):
        self.dict = dict_tree

    @staticmethod
    def _empty_tree():
        return {
            'spec': {},
            'constraints': {},
            'parameters': {}
        }

    @classmethod
    def make_tree(cls, src, parser):
        '''
            Create a new instance tree instance from a valid dictionary
        '''
        params, cons, ctx = parser.parameters, parser.constraints, parser.ctx()

        def mk(src):
            '''
                Recursively build the configuration tree
            '''
            tree = cls._empty_tree()
            if src:
                for key, val in src.items():
                    if key in cons:  # add constraint
                        tree['constraints'][key] = cons[key].with_value(val,
                                                                        ctx)
                    elif key in params:  # add parameter
                        if not cstm_isinstance(val, params[key]['type']):
                            raise ValueTypeError(key, params[key]['type'], val)

                        tree['parameters'][key] = val

                    elif isinstance(val, dict):  # add specialization
                        tree['spec'][SpecKey.parse(key)] = mk(val)

                    else:  # unknown key
                        raise InvalidNodeError(key, val)
            return tree

        return cls(mk(src))

    def copy(self):
        return ConfTree(deepcopy(self.dict))

    def update(self, tree):
        '''
            Update self with values found in tree.
        '''
        def up(old_d, new_d):
            '''
                Recursively update the content of old_d with new_d.
            '''
            for key, cons in new_d['constraints'].items():
                old_d['constraints'][key] = cons

            for param, value in new_d['parameters'].items():
                old_d['parameters'][param] = value

            for spec, spec_d in new_d['spec'].items():
                if spec not in old_d['spec']:
                    old_d['spec'][spec] = self._empty_tree()
                if spec.hardreset:  # simply override
                    old_d['spec'][spec] = deepcopy(spec_d)
                else:  # recursively override individual items
                    up(old_d['spec'][spec], spec_d)

        up(self.dict, tree.dict)

    def get_spec_at(self, path):
        '''
            Get specializations defined at a given node in the tree.
            Return an empty dictionary if the path does not exist.
        '''
        d = self.dict
        for spec in path:
            spec = SpecKey.parse(spec)
            if spec in d['spec']:
                d = d['spec'][spec]
            else:
                return {}
        return [repr(sp) for sp in d['spec']]

    def get_new_params_at(self, path):
        '''
            Get params defined at a given node in the tree.
            Return an empty dictionary if the path does not exist.
        '''
        d = self.dict
        for spec in path:
            spec = SpecKey.parse(spec)
            if spec in d['spec']:
                d = d['spec'][spec]
            else:
                return {}
        return d['parameters']

    def get_new_constraints_at(self, path):
        '''
            Get constraints defined at a given node in the tree.
            Return an empty dictionary if the path does not exist.
        '''
        d = self.dict
        for spec in path:
            spec = SpecKey.parse(spec)
            if spec in d['spec']:
                d = d['spec'][spec]
            else:
                return {}
        return d['constraints']

    def __repr__(self):
        return 'ConfTree({})'.format(self.dict)

    def __eq__(self, other):
        return isinstance(other, ConfTree) and self.dict == other.dict

    def __ne__(self, other):
        return not (self == other)


class ConfParser(object):
    '''
        Test configuration loader and parser. It takes output from yaml parser
        and build the actual configuration tree.
    '''
    def __init__(self):
        self.parameters = {
            'allow_undef': {
                'type': bool,
                'inherited': True,
                'default': True
            }
        }
        self.constraints = {}
        self.metadata = {}

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
                   apply_to='number', use_params=[], exclude=set(),
                   handle_undef=True):
        '''
            Register a constraints to be recognised while parsing config.
        '''
        def register(fun):
            if name is None:
                name_ = fun.__name__
            else:
                name_ = name

            if apply_to == 'this':
                inherited_ = False
            else:
                inherited_ = inherited

            for param in use_params:
                if param not in self.parameters:
                    raise UnknownParamError(name_, param)

            self.constraints[name_] = Constraint(
                name_, fun, value_type, inherited_, use_params,
                set(exclude), apply_to, handle_undef
            )
            return fun
        return register

    def ctx(self):
        return self.metadata.copy()

    def make_trees(self, parsed_src, metadata={}):
        '''
            Create a ConfTree instance from the yaml parser output.
        '''
        assert isinstance(parsed_src, dict), ('parsed_src have to be derivated'
                                              ' from a dictionary but it is'
                                              ' a {}'.format(type(parsed_src)))
        filters = {}
        trees = {}
        self.metadata = metadata.copy()

        if 'filters' in parsed_src:
            # Get filters definitions
            for name, filt in parsed_src['filters'].items():
                if name == '__default' or name in self.constraints \
                   or name in self.parameters:
                    raise IllegalFilterNameError(name)
                filters[name] = IterStateFilter(filt)

                # Parse each filtered tree and remove it from the source tree
                if name in parsed_src:
                    self.metadata['tree'] = name
                    trees[name] = ConfTree.make_tree(parsed_src.pop(name),
                                                     self)
                else:
                    # Should we raise an error ?
                    warn('In YAML config {} filter is defined but not used.'
                         .format(name))

            # Remove the filters field from the dict
            del parsed_src['filters']

        # Parse the fields remaining as the default tree
        self.metadata['tree'] = 'default tree'
        trees['__default'] = ConfTree.make_tree(parsed_src, self)

        return trees, filters
