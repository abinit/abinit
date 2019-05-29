'''
    Define basic structures without particular requirements.
'''
from __future__ import print_function, division, unicode_literals
from ..register_tag import yaml_auto_map, yaml_scalar, yaml_not_available_tag
from ..common import FailDetail
from .pandas_commons import has_pandas


@yaml_auto_map
class Etot(object):
    __yaml_tag = 'ETOT'

    def __init__(self, label='nothing', comment='no comment'):
        self.label = label
        self.comment = comment

    @classmethod
    def from_map(cls, map):
        new = super(Etot, cls).from_map(map)
        new.components = {
            name: value for name, value in new.__dict__.items()
            if name not in [
                'Etotal',
                'label',
                'comment',
                'Band energy',
                'Total energy(eV)'
            ]
        }
        return new


if has_pandas:
    from .pandas_commons import Table

    @yaml_scalar
    class EtotIters(Table):
        is_dict_like = False  # prevent tester from browsing columns

        residues = {
            'deltaE(h)',
            'residm',
            'vres2'
        }

        def last_iter(self, other, **opts):
            tol = opts.get('tol', 1.0e-10)
            ceil = opts.get('ceil', 1.0e-10)
            tol_iter = opts.get('tol_iter', 5)

            def chk_tol(a, b):
                return abs(a - b) / (abs(a) + abs(b)) < tol

            def chk_ceil(b):
                return abs(b) < ceil

            o_n, s_n = other.shape[0] - 1, self.shape[0] - 1
            for key in self:
                oserie, sserie = other[key], self[key]
                # index -1 does not work on series

                if key in self.residues:
                    if not chk_ceil(oserie[o_n]):
                        msg = ('Last item of {} column does not match the'
                               ' ceil {}.')
                        return FailDetail(msg.format(key, ceil))
                else:
                    if not chk_tol(sserie[s_n], oserie[o_n]):
                        msg = ('Last item of {} column does not match the'
                               ' tolerance {}.')
                        return FailDetail(msg.format(key, tol))
            if abs(s_n - o_n) > tol_iter:
                return FailDetail('Difference between number of iteration'
                                  ' is above tol_iter')
            return True
else:
    yaml_not_available_tag('EtotIters', 'Pandas module is not available')
