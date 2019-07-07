'''
    Structures dedicated to e-e self-energy e.g. GW computations
'''
from ..register_tag import yaml_scalar, yaml_auto_map, yaml_not_available_tag
from .pandas_commons import has_pandas


@yaml_auto_map
class SelfEnergy_ee(object):
    """Results for e-e self-energy for a single k-point and spin."""


if has_pandas:
    from .pandas_commons import Table

    @yaml_scalar
    class SigmaeeData(Table):
        """Store Sigma_nk matrix elements"""
else:
    yaml_not_available_tag('SigmaeeData', 'Pandas module is not available')
