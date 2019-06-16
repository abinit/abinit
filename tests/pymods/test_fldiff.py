from __future__ import print_function, division, unicode_literals
import pytest
from .fldiff import (
    Differ,
    LineCountDifference,
    # MetaCharDifference,
    FloatDifference,
    TextDifference,
    ForcedDifference,
)
from .data_extractor import DataExtractor
from .yaml_tools.errors import NoIteratorDefinedError


class TestDiffer:
    lines1 = [
        ' Here are some regular lines of numbers : 0.4546\n',
        ' 5.8787 44.537e+056\n',
        ' .7856 5.0\n',
        '- This lines should be ignored\n',
        '+ This line appears in output but is not counted as errorneous\n',
        ': This line is always compared as characters 0.457896321545\n',
        '. This one too but it will never be counted a errorneous 2.4354364\n',
        '% This line have a fixed tolerance of 1.01e-2: 2.043643 5.5473684\n',
        'P P and , lines are handled according to parameters as + or - lines\n'
    ]

    lines2 = [  # numerical differences
        ' Here are some regular lines of numbers : 0.4546\n',
        ' 5.8707 44.535e+056\n',
        ' .7856 8.0\n',
        '- This lines should be ignored\n',
        '+ This line appears in output but is not counted as errorneous\n',
        ': This line is always compared as characters 0.457896321545\n',
        '. This one too but it will never be counted a errorneous 2.4354364\n',
        '% This line have a fixed tolerance of 1.01e-2: 2.043644 5.4\n',
        'P P and , lines are handled according to parameters as + or - lines\n'
    ]

    lines3 = [  # text differences
        ' Here are some regular lines of floats : 0.4546\n',
        ' 5.8787 44.537e+056\n',
        ' .7856 5.0\n',
        '- Should this line be ignored ?\n',
        '+ This line appears in output but is not counted as errorneous\n',
        ': This line is always compared as characters 0.457896321546\n',
        '. This one too but it will never be counted a errorneous 2.4354360\n',
        '% This line have a fixed tolerance of 1.01e-2: 2.043643 5.5473684\n',
        'P , and P lines are handled according to parameters as + or - lines\n'
    ]

    lines4 = [  # line number differences but not significant
        '- We can append as much ignored lines as we want\n',
        'P It should not change the result\n',
        ' Here are some regular lines of numbers : 0.4546\n',
        ' 5.8787 44.537e+056\n',
        ' .7856 5.0\n',
        '- This lines should be ignored\n',
        '+ This line appears in output but is not counted as errorneous\n',
        '- Ignored lines are everywhere\n',
        ': This line is always compared as characters 0.457896321545\n',
        '. This one too but it will never be counted a errorneous 2.4354364\n',
        '% This line have a fixed tolerance of 1.01e-2: 2.043643 5.5473684\n',
        'P P and , lines are handled according to parameters as + or - lines\n'
    ]

    lines5 = [  # line number differences significant
        ' Here are some regular lines of numbers : 0.4546\n',
        ' 78.73687 98.5763\n',
        ' 5.8787 44.537e+056\n',
        ' .7856 5.0\n',
        '- This lines should be ignored\n',
        '+ This line appears in output but is not counted as errorneous\n',
        ': This line is always compared as characters 0.457896321545\n',
        '. This one too but it will never be counted a errorneous 2.4354364\n',
        '% This line have a fixed tolerance of 1.01e-2: 2.043643 5.5473684\n',
        'P P and , lines are handled according to parameters as + or - lines\n'
    ]

    def test_default(self):
        diff = Differ()
        assert diff.options['tolerance_abs'] == 1.01e-10
        assert diff.options['tolerance_rel'] == 1.01e-10
        assert diff.options['ignore']
        assert diff.options['ignoreP']

    def test_diff_lines_same(self):
        diff = Differ()
        differences = diff._diff_lines(self.lines1, self.lines1)[0]
        assert len(differences) == 1
        assert isinstance(differences[0], ForcedDifference)

    def test_diff_lines_float(self):
        diff = Differ()
        differences = diff._diff_lines(self.lines1, self.lines2)[0]
        print(differences)
        assert len(differences) == 5
        d3 = differences.pop(3)

        assert isinstance(d3, ForcedDifference)
        assert all(isinstance(d, FloatDifference) for d in differences)

    def test_diff_lines_text(self):
        diff = Differ()
        differences = diff._diff_lines(self.lines1, self.lines3)[0]
        assert len(differences) == 4

        d2 = differences.pop(1)
        assert isinstance(d2, ForcedDifference)

        assert all(isinstance(d, TextDifference) for d in differences)

        d4 = differences.pop(2)
        assert d4.silent

        assert not any(d.silent for d in differences)

    def test_diff_lines_number_not_significant(self):
        diff = Differ()
        differences = diff._diff_lines(self.lines1, self.lines4)[0]
        assert len(differences) == 1
        assert isinstance(differences[0], ForcedDifference)

    def test_diff_lines_number_significant(self):
        diff = Differ()
        differences = diff._diff_lines(self.lines1, self.lines5)[0]
        assert len(differences) == 1
        assert isinstance(differences[0], LineCountDifference)
        assert differences[0].more == 'file 2'

    def test_diff_lines_float_format(self):
        diff = Differ()
        differences = diff._diff_lines(
            [' .0007  564.5e-3  7000.0\n'],
            [' 7.0e-4  5.645D-1  7.0f3\n']
        )[0]
        assert len(differences) == 0

    def test_diff_ignore_blanks(self):
        diff = Differ()
        differences = diff._diff_lines(
            [
                ' One normal line\n',
                '.One messy dot   \t line\n',
                ':A colon line.\n',
                ' And a last line\n'
            ],
            [
                ' One  \tnormal line\n\r',
                '.One messy\tdot line  \n',
                ':A colon \t line.\n\r',
                ' And a last line'
            ]
        )[0]
        assert len(differences) == 0


class TestResult:
    '''
        Result only exists to reproduce the historical fldiff.pl so
        the fact that all tests pass is enough.
        It may be removed in the future.
    '''
    pass


class TestDataExtractor:
    def test_default(self):
        dext = DataExtractor(True)
        assert dext.ignore
        assert dext.ignoreP
        assert not dext.xml_mode

    def test_get_metachar(self):
        dext = DataExtractor(True)
        assert dext._get_metachar('-truc') == '-'
        assert dext._get_metachar('+truc') == '+'
        assert dext._get_metachar(' truc') == ' '
        assert dext._get_metachar('.truc') == '.'

        # ignore blank and empty lines
        assert dext._get_metachar('  \t\n') == '-'
        assert dext._get_metachar('') == '-'

        assert dext._get_metachar('Ptruc') == '-'
        assert dext._get_metachar(',truc') == '-'

        dext = DataExtractor(True, ignore=False)
        assert dext._get_metachar(',truc') == '+'

        dext = DataExtractor(True, ignoreP=False)
        assert dext._get_metachar('Ptruc') == '+'

    def test_extract_ignore_minus_meta(self):
        dext = DataExtractor(True)
        lines = [
            '- first',
            'P second',
            '   \t \n',
            ''
        ]

        linesres, _, ignored = dext.extract(lines)
        assert linesres == []
        assert ignored == [(i, line) for i, line in enumerate(lines)]

    def test_extract_keep_all_non_minus(self):
        dext = DataExtractor(True, ignore=False, ignoreP=False)
        lines = [
            '+ first',
            'P second',
            ', third',
            '  fourth',
            '. fifth',
        ]

        linesres, _, ignored = dext.extract(lines)
        assert linesres == [
            (0, '+', '+ first'),
            (1, '+', 'P second'),
            (2, '+', ', third'),
            (3, ' ', '  fourth'),
            (4, '.', '. fifth'),
        ]
        assert ignored == []

    def test_extract_require_iterstart(self):

        dext = DataExtractor(True)
        dext.iterators_state = {'dtset': 1}
        lines = '''\
---
label: a doc
a field: 58
another: 78
a list of strings:
- "a string"
- "two strings"
- "..."
...'''
        lines = [line + '\n' for line in lines.split('\n')]
        with pytest.raises(NoIteratorDefinedError):
            _, documents, _ = dext.extract(lines)

    def test_extract_find_yaml_doc(self):
        dext = DataExtractor(True)
        dext.iterators_state = {'dtset': 1}
        lines = '''\
--- !IterStart
dtset: 1
...

 Garbage here

---
label: a doc
a field: 58
another: 78
a list of strings:
- "a string"
- "two strings"
- "..."
...'''

        lines = [line + '\n' for line in lines.split('\n')]
        _, documents, _ = dext.extract(lines)

        # IterStart documents should not be in the document list
        assert len(documents) == 1
        assert documents['dtset=1 a doc'].iterators == {'dtset': 1}
        assert documents['dtset=1 a doc'].start == 6
        assert documents['dtset=1 a doc'].end == 14
        assert documents['dtset=1 a doc'].lines == lines[6:]
        assert documents['dtset=1 a doc'].obj == {
            'label': 'a doc',
            'a field': 58,
            'another': 78,
            'a list of strings': ['a string', 'two strings', '...']
        }
