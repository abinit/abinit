from __future__ import annotations

"""
Implement the steps to extract data from an Abinit output file.
Extract lines associated with their "meta character" (that makes sense in
fldiff), and valid YAML documents associated with their iteration context.

This module provides the DataExtractor class, which processes ABINIT output files
to identify significant lines (marked with meta-characters for fldiff) and
extract structured YAML documents associated with their iteration context.
"""
import re

from .yaml_tools import Document
from .yaml_tools import is_available as has_yaml
from .yaml_tools.abinit_iterators import ITERATOR_RANKS
from .yaml_tools.errors import DuplicateDocumentError, NoIteratorDefinedError

# Tag is only recognised if it is a valid a word ([A-Za-z0-9_]+)
# It won't recognise serialized tags for example
doc_start_re = re.compile(r"---(?: !(\w+))?\n?$")
doc_end_re = re.compile(r"\.\.\.\n?$")


class DataExtractor:
    """
    Setup extraction of formatted documents and significant lines.

    The DataExtractor scans lines and identifies whether they belong to a YAML
    document or if they are "significant" lines that should be compared using
    the fldiff algorithm. It tracks the state of ABINIT iterators to provide
    context for extracted YAML documents.

    Attributes:
        IGNORE_LINES_STARTING_WITH (list[str]): List of prefixes for lines that should be ignored entirely.
    """

    IGNORE_LINES_STARTING_WITH = [
        "MPI startup(): Warning: I_MPI_PMI_LIBRARY",
        "MPI startup(): PMI server not found.",
    ]

    def __init__(self, use_yaml: bool, ignore: bool = True, ignoreP: bool = True, xml_mode: bool = False):
        """
        Initialize the DataExtractor.

        Args:
            use_yaml (bool): True to use YAML mode.
            ignore (bool): Whether to ignore certain lines.
            ignoreP (bool): Whether to ignore 'P' meta-character lines.
            xml_mode (bool): True if parsing XML documents.
        """
        self.use_yaml = use_yaml and has_yaml
        # do not use fldiff on data that have explicitly been written for YAML use
        self.use_fl_for_yaml = not use_yaml
        self.ignore = ignore
        self.ignoreP = ignoreP
        self.iterators_state = {}
        self.xml_mode = xml_mode
        self.corrupted_docs = []
        self.abinit_messages = []

    def _get_metachar(self, line: str) -> str:
        """
        Return a meta character which gives the behaviour of the line.

        Args:
            line (str): The line to analyze.

        Returns:
            str: The meta-character.
        """
        if not line or line.isspace():  # blank line
            c = "-"
        elif line[0].isspace():
            c = " "
            # dirty fix for compatibility
            # I think xml should not be compared with the basic algorithm
            if self.xml_mode and "timeInfo" in line:
                c = "."
        else:
            c = line[0]
            if c == ",":
                if self.ignore:
                    c = "-"
                else:
                    c = "+"
            elif c == "P":
                if self.ignoreP:
                    c = "-"
                else:
                    c = "+"
        return c

    def ignore_line(self, line: str) -> bool:
        """
        Check if the line should be ignored.

        Args:
            line (str): The line to check.

        Returns:
            bool: True if the line should be ignored.
        """
        if (any(line.startswith(l) for l in self.IGNORE_LINES_STARTING_WITH)): return True
        return False

    def extract(self, src_lines: list[str]) -> tuple[list[tuple[int, str, str]], dict[str, Document], list[tuple[int, str]]]:
        """
        Extract formatted documents and significant lines from src_lines.

        Args:
            src_lines (list): List of strings to process.

        Returns:
            tuple: (significant_lines, documents, ignored_lines) where:
                significant_lines (list): List of (index, metachar, content) tuples.
                documents (dict): Dictionary mapping document IDs to Document objects.
                ignored_lines (list): List of (index, content) tuples for lines not belonging to a doc.
        """
        # Reset internal state to allow several extractions with the same instance.
        self.iterators_state = {}
        self.corrupted_docs = []
        lines, docs, ignored = [], {}, []

        current_doc = None
        for i, line in enumerate(src_lines):

            if self.ignore_line(line): continue

            # TODO
            # Ignore Yaml documents matching e.g. `--- !tagname # fldiff_ignore

            if current_doc is not None:
                # accumulate source lines
                current_doc.lines.append(line)

                if line.startswith("...") and doc_end_re.match(line):
                    # reached the end of the doc
                    if self.use_yaml:
                        current_doc.end = i

                        if getattr(current_doc.obj, "_is_iter_start", False):
                            # special case of IterStart
                            curr_it = current_doc.obj.iterator

                            # Update current iterators state
                            # list freeze the key list to allow deleting in the loop
                            for iterator in list(self.iterators_state):
                                if ITERATOR_RANKS[curr_it] < ITERATOR_RANKS[iterator]:
                                    del self.iterators_state[iterator]
                            self.iterators_state[curr_it] = current_doc.obj.iteration

                        elif current_doc.corrupted:
                            # Signal corruption but ignore the document
                            self.corrupted_docs.append(current_doc)

                        elif getattr(current_doc.obj, "_is_abinit_message", False):
                            # Special case of Warning, Error etc.. store it for later use
                            self.abinit_messages.append(current_doc)

                        elif current_doc.obj is not None:
                            if not current_doc.iterators:
                                # This is not normal!
                                raise NoIteratorDefinedError(current_doc)

                            if current_doc.id in docs:
                                raise DuplicateDocumentError(line, current_doc.id)

                            docs[current_doc.id] = current_doc

                    elif self.use_fl_for_yaml:
                         # let fldiff compare lines if YAML test is disabled
                        lines.extend((current_doc.start + i, " ", " " + line) for i, line in enumerate(current_doc.lines))

                    # go back to normal mode
                    current_doc = None

            elif self._get_metachar(line) == "-":
                # starting a yaml doc
                if line.startswith("---") and doc_start_re.match(line):
                    tag = doc_start_re.match(line).group(1)
                    #iterators_state =
                    current_doc = Document(self.iterators_state.copy(), i, [line], tag=tag)
                else:
                    ignored.append((i, line))
            else:
                # significant line not in a doc
                lines.append((i, self._get_metachar(line), line))

        return lines, docs, ignored
