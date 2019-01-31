from __future__ import print_function, division, unicode_literals
ITERATORS = [
    'idtset',
    'itimimage',
    'iimage',
    'itime',
    'istep'
]

# associate an iterator with its deepness in the global computation
ITERATOR_RANKS = {key: i for i, key in enumerate(ITERATORS)}
