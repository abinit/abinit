'''
    Implement the steps to extract data from an Abinit output file.
'''

from __future__ import print_function, division, unicode_literals
import re
from .yaml_tools import is_available as has_yaml
from .yaml_tools.abinit_iterators import ITERATOR_RANKS
from .yaml_tools.errors import NoIteratorDefinedError
if has_yaml:
    from .yaml_tools import yaml_parse
else:
    def yaml_parse(x):
        pass

doc_start_re = re.compile(r'---( !.*)?\n?$')
doc_end_re = re.compile(r'\.\.\.')


def parse_doc(doc):
    if doc['type'] == 'yaml':
        obj = yaml_parse(''.join(doc['lines']))
        doc['obj'] = obj
    return doc


class DataExtractor:
    '''
        Setup extraction of formated documents and significant lines.
    '''

    def __init__(self, ignore=True, ignoreP=True, xml_mode=False):
        self.ignore = ignore
        self.ignoreP = ignoreP
        self.iterators_state = {}
        self.xml_mode = xml_mode
        self.has_corrupted_doc = False
        self.abinit_messages = []

    def __get_metachar(self, line):
        '''
            Return a metacharacter wich give the behaviour of the line
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
            Extract formated documents and significant lines for the source.
        '''
        # Reset those states to allow several extract with the same instance
        self.iterators_state = {}
        self.has_corrupted_doc = False
        lines, docs, ignored = [], [], []

        current_doc = None
        for i, line in enumerate(src_lines):
            if current_doc is not None:
                # accumulate source lines
                current_doc['lines'].append(line)
                if doc_end_re.match(line):  # reached the end of the doc
                    if not has_yaml:
                        # ignore the document
                        pass
                    else:
                        current_doc['end'] = i
                        # parse source
                        parse_doc(current_doc)

                        if hasattr(current_doc['obj'], '_is_iter_start'):
                            # special case of IterStart
                            curr_it = current_doc['obj'].iterator

                            # Update current iterators state
                            for iterator in self.iterators_state:
                                if ITERATOR_RANKS[curr_it] \
                                   < ITERATOR_RANKS[iterator]:
                                    del self.iterators_state[iterator]
                            self.iterators_state[curr_it] = \
                                current_doc['obj'].iteration

                        elif hasattr(current_doc['obj'], '_is_corrupted_doc'):
                            # Signal corruption but ignore the document
                            self.has_corrupted_doc = True

                        elif hasattr(current_doc['obj'], '_is_abinit_message'):
                            # Special case of Warning, Error etc..
                            # store it for later use
                            self.abinit_messages.append(current_doc)

                        elif current_doc['obj'] is not None:
                            if not current_doc['iterators']:
                                # This is not normal !
                                raise NoIteratorDefinedError(current_doc)
                            docs.append(current_doc)
                    current_doc = None  # go back to normal mode

            else:
                if self.__get_metachar(line) == '-':
                    if doc_start_re.match(line):  # starting a yaml doc
                        current_doc = {
                            'type': 'yaml',
                            # save iterations states
                            'iterators': self.iterators_state.copy(),
                            'start': i,
                            'end': -1,
                            'lines': [line],
                            'obj': None
                        }
                    else:
                        ignored.append((i, line))
                else:  # significant line not in a doc
                    lines.append((i, self.__get_metachar(line), line))

        return lines, docs, ignored
