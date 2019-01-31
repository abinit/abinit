'''
    Implement the steps to extract data from an Abinit output file.
'''

from __future__ import print_function, division, unicode_literals
import re
from .yaml_tools import is_available as has_yaml
from .yaml_tools.abinit_metadata import ITERATOR_RANKS
if has_yaml:
    from .yaml_tools.structures import IterStart
    from .yaml_tools import yaml_parse
else:
    IterStart = object

    def yaml_parse(x):
        pass

doc_start_re = re.compile(r'--- (!.*)?')
doc_end_re = re.compile(r'...')


def parse_doc(doc):
    if has_yaml and doc['type'] == 'yaml':
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

    def __get_metachar(self, line):
        '''
            Return a metacharacter wich give the behaviour of the line independently from options.
        '''
        if not line or line[0].isspace():
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
        lines, documents, ignored = [], [], []

        docmode = False
        current_document = None
        for i, line in enumerate(src_lines):
            if docmode:
                current_document['lines'].append(line)  # accumulate source lines
                if doc_end_re.match(line):  # reached the end of the document
                    current_document['end'] = i
                    # parse source
                    parse_doc(current_document)

                    # special case of IterStart
                    if isinstance(current_document['obj'], IterStart):
                        for iterator in self.iterators_state:
                            if ITERATOR_RANKS[current_document['obj'].iterator] < ITERATOR_RANKS[iterator]:
                                del self.iterators_state[iterator]

                        self.iterators_state[current_document['obj'].iterator] = current_document['obj'].iteration

                    documents.append(current_document)
            else:
                if self.__get_metachar(line) == '-':
                    if doc_start_re.match(line):  # starting a yaml document
                        current_document = {
                            'type': 'yaml',
                            'iterators': self.iterators_state.copy(),  # save iterations states
                            'start': i,
                            'end': -1,
                            'lines': [line],
                            'obj': None
                        }
                    else:
                        ignored.append((i, line))
                else:  # significant line not in a document
                    lines.append((i, self.__get_metachar(line), line))

        return lines, documents, ignored
