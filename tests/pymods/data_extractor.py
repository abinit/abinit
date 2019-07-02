'''
    Implement the steps to extract data from an Abinit output file.
    Extract lines associated with their "meta character" (that make sense in
    fldiff), and valid YAML documents associated with there iteration context.
'''
from __future__ import print_function, division, unicode_literals
import re
from .yaml_tools import Document
from .yaml_tools.abinit_iterators import ITERATOR_RANKS
from .yaml_tools.errors import NoIteratorDefinedError, DuplicateDocumentError

# Tag is only recognised if it is a valid a word ([A-Za-z0-9_]+)
# It won't recognise serialized tags for example
doc_start_re = re.compile(r'---( !!?\w+)?\n?$')
doc_end_re = re.compile(r'\.\.\.\n?$')


class DataExtractor(object):
    '''
        Setup extraction of formatted documents and significant lines.
    '''

    def __init__(self, use_yaml, ignore=True, ignoreP=True, xml_mode=False):
        self.use_yaml = use_yaml
        self.ignore = ignore
        self.ignoreP = ignoreP
        self.iterators_state = {}
        self.xml_mode = xml_mode
        self.corrupted_docs = []
        self.abinit_messages = []

    def _get_metachar(self, line):
        '''
            Return a meta character which gives the behaviour of the line
            independently from options.
        '''
        if not line or line.isspace():  # blank line
            c = '-'
        elif line[0].isspace():
            c = ' '
            # dirty fix for compatibility
            # I think xml should not be compared with the basic algorithm
            if self.xml_mode and 'timeInfo' in line:
                c = '.'
        else:
            c = line[0]
            if c == ',':
                if self.ignore:
                    c = '-'
                else:
                    c = '+'
            elif c == 'P':
                if self.ignoreP:
                    c = '-'
                else:
                    c = '+'
        return c

    def extract(self, src_lines):
        '''
            Extract formatted documents and significant lines from list of strings.
        '''
        # Reset internal state to allow several extractions with the same instance
        self.iterators_state = {}
        self.corrupted_docs = []
        lines, docs, ignored = [], {}, []

        current_doc = None
        for i, line in enumerate(src_lines):
            if current_doc is not None:
                # accumulate source lines
                current_doc.lines.append(line)
                if line.startswith('...') and doc_end_re.match(line):
                    # reached the end of the doc
                    if self.use_yaml:
                        current_doc.end = i

                        if getattr(current_doc.obj, '_is_iter_start', False):
                            # special case of IterStart
                            curr_it = current_doc.obj.iterator

                            # Update current iterators state
                            # list freeze the key list to allow deleting in the
                            # loop
                            for iterator in list(self.iterators_state):
                                if ITERATOR_RANKS[curr_it] \
                                   < ITERATOR_RANKS[iterator]:
                                    del self.iterators_state[iterator]
                            self.iterators_state[curr_it] = \
                                current_doc.obj.iteration

                        elif current_doc.corrupted:
                            # Signal corruption but ignore the document
                            self.corrupted_docs.append(current_doc)

                        elif getattr(current_doc.obj, '_is_abinit_message', False):
                            # Special case of Warning, Error etc..
                            # store it for later use
                            self.abinit_messages.append(current_doc)

                        elif current_doc.obj is not None:
                            if not current_doc.iterators:
                                # This is not normal !
                                raise NoIteratorDefinedError(current_doc)

                            if current_doc.id in docs:
                                raise DuplicateDocumentError(current_doc.id)

                            docs[current_doc.id] = current_doc
                    current_doc = None  # go back to normal mode

            elif self._get_metachar(line) == '-':
                # starting a yaml doc
                if line.startswith('---') and doc_start_re.match(line):
                    current_doc = Document(self.iterators_state.copy(),
                                           i, [line])
                else:
                    ignored.append((i, line))
            else:
                # significant line not in a doc
                lines.append((i, self._get_metachar(line), line))

        return lines, docs, ignored
