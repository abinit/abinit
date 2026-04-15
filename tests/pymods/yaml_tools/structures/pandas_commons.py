"""Base classes for documents containing tabular data represented with pandas DataFrames."""
from .. import has_pandas
from ..register_tag import yaml_not_available_tag, yaml_scalar

if has_pandas:
    from pandas import DataFrame, read_csv
    try:
        from pandas.compat import StringIO
    except ImportError:
        from io import StringIO

    @yaml_scalar
    class Table(DataFrame):
        # assume the class implements a complete dict-like interface
        is_dict_like = True
        table_sep = r"\s+"

        @classmethod
        def from_scalar(cls, scal):
            return cls(read_csv(StringIO(scal), sep=cls.table_sep))

        def to_scalar(self):
            return self.to_string()
else:
    yaml_not_available_tag("Table", "Pandas module is not available")
