from __future__ import print_function, division, unicode_literals
from .errors import EmptySetError, NotOrderedOverlappingSetError
ITERATORS = [
    'dtset',
    'timimage',
    'image',
    'time',
    'step'
]

# associate an iterator with its deepness in the global computation
ITERATOR_RANKS = {key: i for i, key in enumerate(ITERATORS)}


class IntSet(object):
    '''
        Represent a subset of the natural integers.
    '''
    def __init__(self, obj):
        if isinstance(obj, int):
            self._type = 'singleton'
            self.value = obj

            def test(v):
                if isinstance(v, IntSet):
                    if v._type == 'singleton':
                        return v.value == self.value
                    else:
                        return False
                else:
                    return v == obj
        elif isinstance(obj, list):
            self.values = frozenset(obj)

            def test(v):
                if isinstance(v, IntSet):
                    if v._type == 'singleton':
                        return v.value in self
                    elif v._type == 'finite':
                        for val in v.values:
                            if val not in self:
                                return False
                        return True
                    elif v._type == 'bounded':
                        return set(range(v.min, v.max+1)).issubset(self.values)
                    else:
                        return False
                else:
                    return v in self.values
            self._type = 'finite'

        elif isinstance(obj, dict):
            fr = obj.get('from', 1)
            to = obj.get('to', -1)
            if to == -1:
                self._type = 'half-bounded'
                self.min = fr

                def test(v):
                    if isinstance(v, IntSet):
                        if v._type == 'singleton':
                            return v.value in self
                        elif v._type == 'finite':
                            for val in v.values:
                                if val not in self:
                                    return False
                            return True
                        elif v._type == 'half-bounded' or v._type == 'bounded':
                            return v.min >= self.min
                        else:
                            return False
                    else:
                        return v >= fr
            elif to <= fr:
                raise EmptySetError(obj)
            else:
                self._test = 'bounded'
                self.min = fr
                self.max = to

                def test(v):
                    if isinstance(v, IntSet):
                        if v._type == 'singleton':
                            return v.value in self
                        elif v._type == 'finite':
                            for val in v.values:
                                if val not in self:
                                    return False
                            return True
                        elif v._type == 'bounded':
                            return v.min >= self.min and v.max <= self.max
                        else:
                            return False
                    else:
                        return v >= fr and v <= to

        elif obj == 'all':
            def test(v):
                return True
            self._type = 'natural'

        self.__test = test

    def __contains__(self, v):
        return self.__test(v)

    def __repr__(self):
        if self._type == 'singleton':
            return 'IntSet({})'.format(self.value)
        elif self._type == 'finite':
            return 'IntSet({})'.format(', '.join(str(i) for i in self.values))
        elif self._type == 'bounded':
            return 'IntSet({{"from": {}, "to": {} }})'.format(self.min,
                                                              self.max)
        elif self._type == 'half-bounded':
            return 'IntSet({{"from": {} }})'.format(self.min)
        else:
            return 'IntSet("all")'


class IterStateFilter(object):
    '''
        Represent a set of conditions on the iterator state of a document.
    '''
    def __init__(self, d):
        self.filters = {}
        for it in ITERATORS:
            if it in d:
                self.filters[it] = IntSet(d[it])

    def match(self, state):
        for f, s in self.filters.items():
            if state[f] not in s:
                return False
        return True

    def include(self, filt):
        for it in ITERATORS:
            if it in self.filters:
                if it not in filt.filters:
                    # this means that filt.filters['filt'] = all
                    return False
                else:
                    if filt.filters[it] in self.filters[it]:
                        continue
                    else:
                        return False
        return True

    def __repr__(self):
        return ('IterStateFilter({'
                + ', '.join('"{}": {}'.format(n, s)
                            for n, s in self.filters.items())
                + '})')

    @staticmethod
    def cmp(filt1, filt2):
        '''
            Return 1 or -1 if their is a relation of order between the two
            members. Else raise an error.
        '''
        if filt1.include(filt2):
            return -1
        elif filt2.include(filt1):
            return 1
        else:
            raise NotOrderedOverlappingSetError(filt1, filt2)
