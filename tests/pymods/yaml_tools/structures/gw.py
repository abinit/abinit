'''
    Structures dedicated to GW computation
'''
from ..register_tag import yaml_scalar, yaml_auto_map, yaml_not_available_tag
from .pandas_commons import has_pandas


@yaml_auto_map
class GwSigma(object):
    pass


if has_pandas:
    from .pandas_commons import Table

    @yaml_scalar
    class GwSigmaData(Table):
        pass
else:
    yaml_not_available_tag('GwSigmaData', 'Pandas module is not available')
