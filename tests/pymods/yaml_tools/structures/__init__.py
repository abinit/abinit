'''
    Import all available structures defined in this package.
    the possible operations on the extracted data.
'''
from __future__ import print_function, division, unicode_literals
from .. import has_pandas
from ..register_tag import yaml_not_available_tag

from .common_struct import *
from .numpy_struct import *

if has_pandas:
    from .pandas_struct import *
else:
    yaml_not_available_tag('Table', 'Pandas module is not available')
    yaml_not_available_tag('GwSigmaData', 'Pandas module is not available')
    yaml_not_available_tag('EtotIters', 'Pandas module is not available')
