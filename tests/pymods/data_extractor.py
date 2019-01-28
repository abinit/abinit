'''
Implement the steps to extract data from an Abinit output file.
'''

from __future__ import print_function, division, unicode_literals
import re
from yaml_tools import yaml_parse, yaml_object


docstart_re = re.compile(r'--- (!.*)?')
docend_re = re.compile(r'...')

iterator_ranks = {  # associate an iterator with its deepness in the global computation
    key: i for i, key in enumerate([
        'idtset',
        'itimimage',
        'iimage',
        'itime',
        'istep'
    ])
}


# this decorator makes IterStart a suitable object for the !IterStart tag and
# defines defaults values for attributes corresponding to the keys of iterator_ranks
@yaml_object(**{key: -1 for key in iterator_ranks})
class IterStart(object):
    @property  # define a getter for self.iterator
    def iterator(self):
        for key in iterator_ranks:
            if self.__getattr__(key) != -1:
                return key
        # I probably should raise an error.

    @property
    def iteration(self):  # define a getter for self.iteration
        for key in iterator_ranks:
            if self.__getattr__(key) != -1:
                return self.__getattr__(key)
        # I probably should raise an error.


def parse_doc(doc):
    if doc['type'] == 'yaml':
        obj = yaml_parse(''.join(doc['lines']))
        doc['obj'] = obj
    return doc


class DataExtractor:

    def __init__(self, ignore=True, ignoreP=True):
        self.ignore = ignore
        self.ignoreP = ignoreP
        self.iterators_state = {}

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
            Split lines into two groups: ignored and analysed
        '''
        lines, documents, ignored = [], [], []

        docmode = False
        current_document = None
        for i, line in enumerate(src_lines):
            if docmode:
                current_document['lines'].append(line)
                if docend_re.match(line):
                    current_document['end'] = i
                    parse_doc(current_document)
                    if isinstance(current_document['obj'], IterStart):
                        for iterator in self.iterators_state:
                            if iterator_ranks[current_document['obj'].iterator] < iterator_ranks[iterator]:
                                del self.iterators_state[iterator]
                        self.iterators_state[current_document['obj'].iterator] = current_document['obj'].iteration
                    documents.append(current_document)
            else:
                if self.__get_metachar(line) == '-':
                    if docstart_re.match(line):
                        current_document = {
                            'type': 'yaml',
                            'iterations': self.iterator_ranks.copy(),
                            'start': i,
                            'lines': [line],
                            'obj': None
                        }
                    else:
                        ignored.append((i, line))
                else:
                    lines.append((i, line))

        return lines, documents, ignored
