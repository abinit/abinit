from __future__ import print_function, division, unicode_literals
import pytest
from .errors import (EmptySetError, NotOrderedOverlappingSetError)
from .abinit_iterators import IterStateFilter, iter_state_cmp
from .conf_parser import conf_parser
from .meta_conf_parser import ConfTree, ConfParser


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


class TestMetaConfParser(object):
    src1 = {
        'doc1': {
            'sp1': {},
            'sp2': {
                'subsp': {}
            }
        },
        'doc2': {
            'sp1': {
                'subsp': {
                    'subsubsp': {}
                }
            },
            'sp2': {
                'subsp': {}
            }

        }
    }

    src2 = {
        'doc1': {
            'tol_abs': 3.0,
            'sp1': {
                'tol_abs': 1.0,
                'tol_rel': 1.0
            },
            'sp2': {
                'ceil': 5.0,
                'subsp': {
                    'ceil': 7.0
                }
            }
        },
    }

    src3 = {
        'doc1': {
            'tol_abs': 5.0,
            'sp2': {
                'ceil': 8.0,
                'subsp2': {
                    'tol_abs': 12.0
                }
            }
        },
        'doc2': {
            'tol_abs': 1.1,
            'sp3': {
                'tol_rel': 1.5
            }
        }
    }

    src23 = {  # src2.update(src3)
        'doc1': {
            'tol_abs': 5.0,
            'sp1': {
                'tol_abs': 1.0,
                'tol_rel': 1.0
            },
            'sp2': {
                'ceil': 8.0,
                'subsp': {
                    'ceil': 7.0
                },
                'subsp2': {
                    'tol_abs': 12.0
                }
            }
        },
        'doc2': {
            'tol_abs': 1.1,
            'sp3': {
                'tol_rel': 1.5
            }
        }
    }

    src4 = {
        'filters': {
            'f1': {
                'dtset': 5,
                'image': [2, 5, 7]
            },
            'f2': {
                'dtset': {'from': 7},
                'image': {'to': 8}
            }
        },
        'f1': {},
        'f2': {}
    }

    def test_make_tree_empty(self):
        cp = ConfParser()
        trees, _ = cp.make_trees({})
        assert trees['__default'] == ConfTree({'spec': {}, 'constraints': {},
                                              'parameters': {}})

    def test_make_tree_specs(self):
        cp = ConfParser()
        trees, _ = cp.make_trees(self.src1)
        assert trees['__default'] == ConfTree({
            'spec': {
                'doc1': {
                    'spec': {
                        'sp1': {
                            'spec': {},
                            'constraints': {},
                            'parameters': {}
                        },
                        'sp2': {
                            'spec': {
                                'subsp': {
                                    'spec': {},
                                    'constraints': {},
                                    'parameters': {}
                                }
                            },
                            'constraints': {},
                            'parameters': {}
                        },
                    },
                    'constraints': {},
                    'parameters': {}
                },
                'doc2': {
                    'spec': {
                        'sp1': {
                            'spec': {
                                'subsp': {
                                    'spec': {
                                        'subsubsp': {
                                            'spec': {},
                                            'constraints': {},
                                            'parameters': {}
                                        }
                                    },
                                    'constraints': {},
                                    'parameters': {}
                                }
                            },
                            'constraints': {},
                            'parameters': {}
                        },
                        'sp2': {
                            'spec': {
                                'subsp': {
                                    'spec': {},
                                    'constraints': {},
                                    'parameters': {}
                                }
                            },
                            'constraints': {},
                            'parameters': {}
                        },
                    },
                    'constraints': {},
                    'parameters': {}
                }
            },
            'constraints': {},
            'parameters': {}
        })

    def test_make_tree_constraints(self):
        cp = ConfParser()

        @cp.constraint()
        def tol_abs(tol, ref, test):
            pass

        @cp.constraint()
        def tol_rel(tol, ref, test):
            pass

        @cp.constraint()
        def ceil(tol, ref, test):
            pass

        trees, _ = cp.make_trees(self.src2)
        ref = ConfTree({
            'spec': {
                'doc1': {
                    'spec': {
                        'sp1': {
                            'spec': {},
                            'constraints': {
                                'tol_abs':
                                    cp.constraints['tol_abs'].with_value(1.0),
                                'tol_rel':
                                    cp.constraints['tol_rel'].with_value(1.0),
                            },
                            'parameters': {}
                        },
                        'sp2': {
                            'spec': {
                                'subsp': {
                                    'spec': {},
                                    'constraints': {
                                        'ceil': (cp.constraints['ceil']
                                                 .with_value(7.0))
                                    },
                                    'parameters': {}
                                }
                            },
                            'constraints': {
                                'ceil': cp.constraints['ceil'].with_value(5.0)
                            },
                            'parameters': {}
                        },
                    },
                    'constraints': {
                        'tol_abs': cp.constraints['tol_abs'].with_value(3.0)
                    },
                    'parameters': {}
                }
            },
            'constraints': {},
            'parameters': {}
        })
        assert trees['__default'] == ref

    def test_make_tree_filters(self):
        cp = ConfParser()

        trees, filters = cp.make_trees(self.src4)

        assert filters == {
            'f1': IterStateFilter({'dtset': 5, 'image': [2, 5, 7]}),
            'f2': IterStateFilter({'dtset': {'from': 7},
                                   'image': {'from': 1, 'to': 8}})
        }

        assert trees == {
            '__default': ConfTree({
                'spec': {},
                'constraints': {},
                'parameters': {}
            }),
            'f1': ConfTree({
                'spec': {},
                'constraints': {},
                'parameters': {}
            }),
            'f2': ConfTree({
                'spec': {},
                'constraints': {},
                'parameters': {}
            })
        }

    def test_update_tree(self):
        cp = ConfParser()

        @cp.constraint()
        def tol_abs(tol, ref, test):
            pass

        @cp.constraint()
        def tol_rel(tol, ref, test):
            pass

        @cp.constraint()
        def ceil(tol, ref, test):
            pass

        ref = cp.make_trees(self.src23)[0]['__default']

        test = cp.make_trees(self.src2)[0]['__default']
        test.update(cp.make_trees(self.src3)[0]['__default'])

        assert test == ref


class TestConfParser(object):
    def test_conf_parser(self):
        assert conf_parser
