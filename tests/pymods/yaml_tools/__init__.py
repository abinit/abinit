"""
This package gathers all tools used by the Abinit test suite for manipulating YAML formatted data.
"""
from __future__ import annotations

import warnings
from typing import Any

from .errors import NoYAMLSupportError, TagMismatchError, UntaggedDocumentError

try:
    import numpy  # numpy is also required
    import yaml
    is_available = True

except ImportError:
    warnings.warn("\nCannot import numpy or yaml package.\nUse `pip install numpy pyyaml --user`"
                  "\nto install the packages in user mode.")
    is_available = False

try:
    import pandas
    has_pandas = True
except ImportError:
    has_pandas = False
    warnings.warn("\nCannot import pandas package. Use `pip install pandas --user`"
                  "\nto install the package in user mode.")


if is_available:
    # use the Yaml C binding (faster) if possible
    if hasattr(yaml, "CSafeLoader"):
        Loader = yaml.CSafeLoader
    else:
        warnings.warn("The libyaml binding is not available, tests will take"
                      " more time. Using python3 may solve the problem. If it"
                      " doesn't, you may have to install libyaml yourself.")
        Loader = yaml.SafeLoader

    from .common import get_yaml_tag, string

    def yaml_parse(content: str, *args: Any, **kwargs: Any) -> Any:
        from . import structures
        return yaml.load(content, *args, Loader=Loader, **kwargs)

    yaml_print = yaml.dump


class Document:
    """
    A document with all its metadata extracted from the original file.
    """
    def __init__(self, iterators: dict[str, Any], start: int, lines: list[str], tag: str | None = None) -> None:
        """
        Initialize the Document object.

        Args:
            iterators (dict): State of the iterators for this document.
            start (int): Starting line number in the original file.
            lines (list): List of lines belonging to the document.
            tag (str, optional): YAML tag of the document.
        """
        self.iterators = iterators
        self.start = start
        self.end = -1
        self.lines = lines
        self._tag = tag
        self._obj: Any = None
        self._corrupted = False
        self._id: str | None = None

    def _parse(self) -> None:
        """
        Parse lines, set `obj` property.
        Raise an error if the document is untagged.
        """
        if is_available:
            content = "\n".join(self.lines)
            try:
                self._obj = yaml_parse(content)
            except yaml.YAMLError as e:
                print("Exception in Document._parse()\ncontent:\n", content, "\nException:\n", e)
                self._obj = e
                self._corrupted = True
                self._tag = "Corrupted"

            # use type in instead of isinstance because inheritance is fine
            if type(self._obj) in {dict, list, tuple, string}:
                raise UntaggedDocumentError(self.start)
            tag = get_yaml_tag(type(self._obj))
            if self._tag is not None and tag != self._tag:
                self._corrupted = True
                self._obj = TagMismatchError(self.start, tag, self._tag)
            else:
                self._tag = tag

            # MG: Get iterators at this level.
            #self.iterators = self._obj["iterator_state"]
        else:
            raise NoYAMLSupportError("Try to access YAML document but YAML is"
                                     " not available in this environment.")

    @property
    def id(self) -> str:
        """
        Produce a string id that should be unique.
        """
        # MG: FIXME: Well this is not unique. I don't think a document should have an id!
        if self._id is None:
            state = []
            for key, val in self.iterators.items():
                state.append(f"{key}={val}")

            self._id = ",".join(state) + " " + self.tag
        return self._id

    @property
    def obj(self) -> Any:
        """
        The python object constructed by Pyyaml.
        """
        if self._obj is None: self._parse()
        return self._obj

    @property
    def tag(self) -> str | None:
        """
        The document tag.
        """
        if self._tag is None: self._parse()
        return self._tag

    @property
    def corrupted(self) -> bool:
        """
        True if Yaml document is corrupted.
        """
        if self._obj is None: self._parse()
        return self._corrupted
