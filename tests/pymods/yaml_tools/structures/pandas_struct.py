from ..common import FailDetail
from pandas import read_csv, DataFrame
from pandas.compat import StringIO
from ..register_tag import yaml_scalar


@yaml_scalar
class Table(DataFrame):
    # imply that the class implement a complete dict-like interface
    is_dict_like = True
    table_sep = r'\s+'

    @classmethod
    def from_scalar(cls, scal):
        return cls(read_csv(StringIO(scal), sep=cls.table_sep))

    def to_scalar(self):
        return self.to_string()


@yaml_scalar
class GwSigmaData(Table):
    pass


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
