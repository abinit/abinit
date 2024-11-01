"""
Define classes and contants to represent the state of iteration of a document
as well as the operations possible on this state. This is used in filter
applications by the configuration handler.
"""
from __future__ import print_function, division, unicode_literals

from .errors import EmptySetError, NotOrderedOverlappingSetError

ITERATORS = [  # order matters
    'dtset',
    'timimage',
    'image',
    'time',
    'step',
]

# associate an iterator with its deepness in the global computation
ITERATOR_RANKS = {key: i for i, key in enumerate(ITERATORS)}


class IntSet(object):
    """
    Represent a subset of the natural integers.
    """
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
            self._type = 'finite'
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
                        return set(range(v.min, v.max + 1)).issubset(self.values)
                    else:
                        return False
                else:
                    return v in self.values

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
                self._type = 'bounded'
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
            self._type = 'natural'

            def test(v):
                return True

        else:
            raise TypeError('Unknown input for IntSet: {}'.format(obj))

        self._test = test

    def __contains__(self, v):
        return self._test(v)

    def __eq__(self, other):
        return isinstance(other, IntSet) and self in other and other in self

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        if self._type == 'singleton':
            return 'IntSet({})'.format(self.value)
        elif self._type == 'finite':
            return 'IntSet([{}])'.format(
                ', '.join(str(i) for i in self.values)
            )
        elif self._type == 'bounded':
            return 'IntSet({{"from": {}, "to": {} }})'.format(self.min,
                                                              self.max)
        elif self._type == 'half-bounded':
            return 'IntSet({{"from": {} }})'.format(self.min)
        else:
            return 'IntSet("all")'


class IterStateFilter(object):
    """
    Represent a set of conditions on the iterator state of a document.
    Alternatively it can be seen as a cartesian product of subsets of
    the natural integers. The implicit subset for each component is N*
    ({1, 2, 3, 4...}).
    For example IterStateFilter({'dtset': 4, 'image': {1, 5}}) is
    {4} x N* x {1, 2, 3, 4, 5} x N* x N*
    """
    def __init__(self, d):
        self.filters = {}
        for it in ITERATORS:
            if it in d:
                self.filters[it] = IntSet(d[it])

    def match(self, state):
        """
        Does a given state match this filter?
        Is a given tuple in this set?
        """
        for it, int_set in self.filters.items():
            if it in state and state[it] not in int_set:
                return False
        return True

    def include(self, filt):
        """
        Return True if filt is included (see the set interpretation
        in class docstring) in self, False otherwise.
        """
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

    def cmp(self, other):
        """
        Return 1 or -1 if their is a relation of order between the two
        members else raise an error.
        """
        assert isinstance(other, IterStateFilter), (
            "IterStateFilter cannot be compared with {}".format(other)
        )
        if self.include(other):
            if other.include(self):  # A c B & B c A => A = B
                return 0
            else:
                return -1
        elif other.include(self):
            return 1
        else:
            raise NotOrderedOverlappingSetError(self, other)

    def __eq__(self, other):
        if not isinstance(other, IterStateFilter):
            return False
        for it in ITERATORS:
            if it in self.filters:
                if it in other.filters:
                    if self.filters[it] != other.filters[it]:
                        return False
                else:
                    return False
            else:
                if it in other.filters:
                    return False

        return True

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self.cmp(other) == 1

    def __le__(self, other):
        return self.cmp(other) >= 0

    def __gt__(self, other):
        return self.cmp(other) == -1

    def __ge__(self, other):
        return self.cmp(other) <= 0
