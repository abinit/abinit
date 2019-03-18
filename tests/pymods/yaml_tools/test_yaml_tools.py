from __future__ import print_function, division, unicode_literals
import pytest
from .errors import (EmptySetError, NotOrderedOverlappingSetError)
from .abinit_iterators import IterStateFilter
from .meta_conf_parser import ConfTree, ConfParser, SpecKey
from .driver_test_conf import DriverTestConf


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

        assert f1 > f2
        assert f2 < f1
        assert f1 > f3
        assert f3 < f1
        with pytest.raises(NotOrderedOverlappingSetError):
            f2 < f3

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

        assert sorted([f4, f1, f2], reverse=True) == [f1, f2, f4]
        with pytest.raises(NotOrderedOverlappingSetError):
            sorted([f4, f3, f1, f2])


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

    tree1 = {
        'spec': {
            SpecKey('doc1'): {
                'spec': {
                    SpecKey('sp1'): {
                        'spec': {},
                        'constraints': {},
                        'parameters': {}
                    },
                    SpecKey('sp2'): {
                        'spec': {
                            SpecKey('subsp'): {
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
            SpecKey('doc2'): {
                'spec': {
                    SpecKey('sp1'): {
                        'spec': {
                            SpecKey('subsp'): {
                                'spec': {
                                    SpecKey('subsubsp'): {
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
                    SpecKey('sp2'): {
                        'spec': {
                            SpecKey('subsp'): {
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

    src5 = {
        'doc1': {
            'sp1!': {
                'tol_abs': 2.22
            },
            'sp2!': {
                'tol_abs': 1.11
            }
        }
    }

    src25 = {
        'doc1': {
            'tol_abs': 3.0,
            'sp1': {
                'tol_abs': 2.22
            },
            'sp2': {
                'tol_abs': 1.11
            }
        },
    }

    def test_make_tree_empty(self):
        cp = ConfParser()
        trees, _ = cp.make_trees({})
        assert trees['__default'] == ConfTree({'spec': {}, 'constraints': {},
                                              'parameters': {}})

    def test_make_tree_specs(self):
        cp = ConfParser()
        trees, _ = cp.make_trees(self.src1)
        assert trees['__default'] == ConfTree(self.tree1)

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
                SpecKey('doc1'): {
                    'spec': {
                        SpecKey('sp1'): {
                            'spec': {},
                            'constraints': {
                                'tol_abs':
                                    cp.constraints['tol_abs'].with_value(1.0),
                                'tol_rel':
                                    cp.constraints['tol_rel'].with_value(1.0),
                            },
                            'parameters': {}
                        },
                        SpecKey('sp2'): {
                            'spec': {
                                SpecKey('subsp'): {
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

        ref1 = IterStateFilter({'dtset': 5, 'image': [2, 5, 7]})
        ref2 = IterStateFilter({'dtset': {'from': 7},
                                'image': {'from': 1, 'to': 8}})

        assert filters == {'f1': ref1, 'f2': ref2}
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

    def test_update_tree_hardreset(self):
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

        ref = cp.make_trees(self.src25)[0]['__default']

        test = cp.make_trees(self.src2)[0]['__default']
        test.update(cp.make_trees(self.src5)[0]['__default'])

        assert test == ref


class TestDriverTestConf(object):
    src1 = '''\
tol_abs: 1.2e-7
sp1:
    sp3:
        tol_abs: 1.8e-9

dt1:
    tol_abs: 1.3e-6
    tol_rel: 1.4e-8
    sp1:
        sp2 with spaces:
            ceil: 1.0e-8

im1:
    sp1:
        ceil: 1.5e-8

filters:
    dt1:
        dtset: 1
        image: {from: 1, to: 8}
    im1:
        dtset: 1
        image: 1
'''

    src2 = '''\
tol_abs: 1.2e-7
sp1:
    sp3:
        tol_abs: 1.8e-9

dt1:
    tol_abs: 1.3e-6
    tol_rel: 1.4e-8
    sp1:
        sp2:
            ceil: 1.0e-8

im1:
    sp1!:
        ceil: 1.5e-8

filters:
    dt1:
        dtset: 1
        image: {from: 1, to: 8}
    im1:
        dtset: 1
        image: 1
'''

    def test_get_constraints(self):
        DriverTestConf.default_conf = '/dev/null'
        driver = DriverTestConf(src=self.src1)

        with driver.use_filter({'dtset': 1, 'image': 4}):
            with driver.go_down('sp1'):
                with driver.go_down('sp2 with spaces'):
                    constraints = driver.get_constraints_for(1.0)
                    assert len(constraints) == 1
                    assert constraints[0].name == 'ceil'
                    assert constraints[0].value == 1.0e-8
                with driver.go_down('sp3'):
                    constraints = driver.get_constraints_for(1.0)
                    assert len(constraints) == 2
                    for c in constraints:
                        assert c.name in ('tol_abs', 'tol_rel')
                        if c.name == 'tol_abs':
                            assert c.value == 1.8e-9
                        elif c.name == 'tol_rel':
                            assert c.value == 1.4e-8

        with driver.use_filter({'dtset': 1, 'image': 1}):
            with driver.go_down('sp1'):
                with driver.go_down('sp2 with spaces'):
                    constraints = driver.get_constraints_for(1.0)
                    assert len(constraints) == 1
                    assert constraints[0].name == 'ceil'
                    assert constraints[0].value == 1.0e-8
                with driver.go_down('sp3'):
                    constraints = driver.get_constraints_for(1.0)
                    print(constraints)
                    assert len(constraints) == 1
                    assert constraints[0].name == 'tol_abs'
                    assert constraints[0].value == 1.8e-9
                with driver.go_down('sp4'):
                    constraints = driver.get_constraints_for(1.0)
                    assert len(constraints) == 1
                    assert constraints[0].name == 'ceil'
                    assert constraints[0].value == 1.5e-8

    def test_get_constraints_hardreset(self):
        DriverTestConf.default_conf = '/dev/null'
        driver = DriverTestConf(src=self.src2)

        with driver.use_filter({'dtset': 1, 'image': 1}):
            with driver.go_down('sp1'):
                with driver.go_down('sp2'):
                    constraints = driver.get_constraints_for(1.0)
                    assert len(constraints) == 1
                    assert constraints[0].name == 'ceil'
                    assert constraints[0].value == 1.5e-8
                with driver.go_down('sp3'):
                    constraints = driver.get_constraints_for(1.0)
                    print(constraints)
                    assert len(constraints) == 1
                    assert constraints[0].name == 'ceil'
                    assert constraints[0].value == 1.5e-8
                with driver.go_down('sp4'):
                    constraints = driver.get_constraints_for(1.0)
                    assert len(constraints) == 1
                    assert constraints[0].name == 'ceil'
                    assert constraints[0].value == 1.5e-8


class TestTester(object):
    pass
