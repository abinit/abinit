'''
    Define basic structures without particular requirements.
'''
from __future__ import print_function, division, unicode_literals
from ..register_tag import yaml_auto_map, yaml_scalar, yaml_not_available_tag
from ..common import FailDetail
from .pandas_commons import has_pandas


@yaml_auto_map
class Etot(object):
    '''
    Component if total energy.
    '''
    __yaml_tag = 'ETOT'

    not_components = {
        'Etotal',
        'comment',
        'Band energy',
        'Total energy(eV)'
    }

    def __init__(self, comment='no comment'):
        self.comment = comment

    @classmethod
    def from_map(cls, map):
        new = super(Etot, cls).from_map(map)
        new.components = {
            name: value for name, value in new.__dict__.items()
            if name not in cls.not_components
        }
        return new


@yaml_auto_map
class EtotDC(object):
    '''
    Components of total energy in Double Counting.
    '''
    not_components = {
        'Etotal (DC)',
        'comment',
        'Band energy',
        '-kT*entropy',
        'Total DC energy(eV)'
    }


@yaml_auto_map
class ResultsGS(object):
    '''
    Miscellaneous results from ground state computations.
    '''


if has_pandas:
    from .pandas_commons import Table

    @yaml_scalar
    class EtotIters(Table):
        is_dict_like = False  # prevent tester from browsing columns

        def last_iter(self, other, **opts):
            '''
                Expects opts to be a dictionary with keys being column names and
                values being 'ceil': ceiling_tol_value or 'tol': tolerance_value.
                The checks are only performed on the last values of each columns.
                An additional optional key of opts is 'tol_iter' giving a tolerance
                for the variation of number of iterations. The default value is 5.
            '''
            tol_iter = opts.get('tol_iter', 5)

            def chk_tol(a, b, tol):
                return abs(a - b) / (abs(a) + abs(b)) < tol

            def chk_ceil(b, ceil):
                return abs(b) < ceil

            o_n, s_n = other.shape[0] - 1, self.shape[0] - 1
            for key in self:
                oserie, sserie = other[key], self[key]
                # index -1 does not work on series

                if key in opts:  # for each column look for a constraint
                    if 'ceil' in opts[key]:
                        ceil = opts[key]['ceil']
                        if not chk_ceil(oserie[o_n], ceil):
                            msg = ('Last item of {} column does not match the'
                                   ' ceil {}: value is {}.')
                            return FailDetail(
                                msg.format(key, ceil, oserie[o_n])
                            )
                    if 'tol' in opts[key]:
                        tol = opts[key]['tol']
                        if not chk_tol(sserie[s_n], oserie[o_n], tol):
                            msg = ('Last item of {} column does not match the'
                                   ' tolerance {}: difference is {}.')
                            return FailDetail(
                                msg.format(key, tol, sserie[s_n] - oserie[o_n])
                            )
            if abs(s_n - o_n) > tol_iter:
                return FailDetail('Difference between number of iteration'
                                  ' is above tol_iter')
            return True

else:
    yaml_not_available_tag('EtotIters', 'Pandas module is not available')
