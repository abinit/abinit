'''
The tokens recognised by the configuration parser are declared here.
There are two kind of tokens:
- the parameters are just variables with a default value. They are used by constraints.
- the constraints represent an actual check of the data. They consist in a
function returning a boolean. The docstring of the function is used as the
description for the testconf_explorer. There is no point in defining
value_type, inherited, apply_to or use_params in the docstring but you
should explain when the test fails and when it succeeds.
'''
from __future__ import print_function, division, unicode_literals

from numpy import ndarray
from numpy.linalg import norm
from .meta_conf_parser import ConfParser
from .errors import MissingCallbackError

conf_parser = ConfParser()

# Parameters
# arguments are token, default=None, value_type=float, inherited=True
conf_parser.parameter('tol_eq', default=1e-8, inherited=True)


# Constraints
# default arguments for constraints are:
# value_type=float, inherited=True, apply_to='number' use_params=[], exclude={}
# handle_undef = True
@conf_parser.constraint(exclude={'tol', 'ceil', 'ignore'})
def tol_rel(tol, ref, tested):
    '''
        Valid if the relative difference between the values is below the
        given tolerance.
    '''
    if abs(ref) + abs(tested) == 0.0:
        return True
    return abs(ref - tested) / (abs(ref) + abs(tested)) < tol


@conf_parser.constraint(exclude={'tol', 'ceil', 'ignore'})
def tol_abs(tol, ref, tested):
    '''
        Valid if the absolute difference between the values is below the
        given tolerance.
    '''
    return abs(ref - tested) < tol


@conf_parser.constraint(apply_to='Array', inherited=True)
def tol_vec(tol, ref, tested):
    '''
        Valid if the cartesian norm of the vector (ref - tested) is below
        the given tolerance.
    '''
    return norm(ref - tested) < tol


@conf_parser.constraint(exclude={'tol', 'tol_abs', 'tol_rel', 'ignore'})
def ceil(ceil_val, ref, tested):
    '''
        Valid if the absolute value of the tested file tested is below the
        given tolerance.
    '''
    return abs(tested) < ceil_val


@conf_parser.constraint(value_type=bool, exclude={'ceil', 'tol', 'tol_rel',
                                                  'tol_abs'})
def ignore(yes, ref, tested):
    '''
        Override numbers tests and always return the same result.
    '''
    return yes


@conf_parser.constraint(value_type=str, inherited=False, apply_to='this',
                        use_params=['tol_eq'])
def equation(eq, ref, tested, tol_eq):
    '''
        If the given expression returns a number, its absolute value is compared
        to the optional tol_eq parameter. If the expression returns a vector,
        it is the euclidian norm that is used.
    '''
    res = eval(eq, {}, {'this': tested, 'ref': ref})
    if isinstance(res, ndarray):
        return norm(res) < tol_eq
    else:
        return abs(res) < tol_eq


@conf_parser.constraint(value_type=list, inherited=False, apply_to='this',
                        use_params=['tol_eq'])
def equations(eqs, ref, tested, tol_eq):
    '''
        See equation. Same thing with a list of equations.
    '''
    for eq in eqs:
        res = eval(eq, {}, {'this': tested, 'ref': ref})
        if isinstance(res, ndarray):
            if norm(res) >= tol_eq:
                return False
        else:
            if abs(res) >= tol_eq:
                return False


@conf_parser.constraint(value_type=dict, apply_to='this')
def callback(locs, ref, tested):
    '''
        Call a method of the reference data with the tested data as first parameter
        and the other parameters as keyword arguments. Return the result of the call.
    '''
    method = locs.pop('method')
    if hasattr(ref, method):
        return getattr(ref, method)(tested, **locs)
    else:
        raise MissingCallbackError(ref, method)
