import pytest
from .errors import (EmptySetError, NotOrderedOverlappingSetError)
from .abinit_iterators import IterStateFilter, iter_state_cmp


class TestStateFilter(object):
    def test_empty(self):
        with pytest.raises(EmptySetError):
            IterStateFilter({'dtset': {'from': 5, 'to': 2}})

    def test_match_singleton(self):
        f1 = IterStateFilter({'dtset': 8})
        assert not f1.match({'dtset': 4})
        assert not f1.match({'dtset': 9})
        assert f1.match({'dtset': 8})

    def test_match_finite(self):
        f1 = IterStateFilter({'dtset': [1, 2, 3, 8]})
        assert not f1.match({'dtset': 4})
        assert not f1.match({'dtset': 9})
        assert f1.match({'dtset': 1})
        assert f1.match({'dtset': 3})
        assert f1.match({'dtset': 8})

    def test_match_bounded(self):
        f1 = IterStateFilter({'dtset': {'from': 1, 'to': 5}})
        assert f1.match({'dtset': 1})
        assert f1.match({'dtset': 4, 'image': 7})
        assert f1.match({'dtset': 5})
        assert not f1.match({'dtset': 6})

    def test_match_half_bounded(self):
        f1 = IterStateFilter({'dtset': {'from': 5}})
        assert not f1.match({'dtset': 4})
        assert f1.match({'dtset': 5})
        assert f1.match({'dtset': 50000})

    def test_include(self):
        f1 = IterStateFilter({
            'dtset': {'from': 5, 'to': 8},
            'image': [1, 2, 3, 5],
        })

        f2 = IterStateFilter({
            'dtset': [7, 6],
            'image': [1, 3, 5],
        })

        f3 = IterStateFilter({
            'dtset': [7, 6],
            'image': [1, 2, 5],
        })

        assert not f2.include(f1)
        assert f1.include(f2)

        assert f1.include(f3)
        assert not f3.include(f1)

        assert not f2.include(f3)
        assert not f3.include(f2)

    def test_cmp(self):
        f1 = IterStateFilter({
            'dtset': {'from': 5, 'to': 8},
            'image': [1, 2, 3, 5],
        })

        f2 = IterStateFilter({
            'dtset': [7, 6],
            'image': [1, 3, 5],
        })

        f3 = IterStateFilter({
            'dtset': [7, 6],
            'image': [1, 2, 5],
        })

        assert iter_state_cmp(f1, f2) == -1
        assert iter_state_cmp(f2, f1) == 1
        assert iter_state_cmp(f1, f3) == -1
        assert iter_state_cmp(f3, f1) == 1
        with pytest.raises(NotOrderedOverlappingSetError):
            iter_state_cmp(f2, f3)

    def test_sort(self):
        f1 = IterStateFilter({
            'dtset': {'from': 5, 'to': 8},
            'image': [1, 2, 3, 5],
        })

        f2 = IterStateFilter({
            'dtset': [7, 6],
            'image': [1, 3, 5],
        })

        f3 = IterStateFilter({
            'dtset': [7, 6],
            'image': [1, 2, 5],
        })

        f4 = IterStateFilter({
            'dtset': [7, 6],
            'image': [1, 3, 5],
            'time': {'from': 2, 'to': 5}
        })

        assert sorted([f4, f1, f2], key=IterStateFilter.key) == [f1, f2, f4]
        with pytest.raises(NotOrderedOverlappingSetError):
            sorted([f4, f3, f1, f2], key=IterStateFilter.key)
