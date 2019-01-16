'''
Test suite for the fldiff module.
Usage:
    python -m test_fldiff
'''
import unittest
from fldiff import Differ, Result
import fldiff


class TestFlDiffResult(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        unittest.TestCase.__init__(self, *args, **kwargs)

        self.r1 = Result([])
        self.r1.ndiff_lines = 13
        self.r1.max_abs_err = 1.5e-5
        self.r1.max_rel_err = 2.1e-6
        self.r1.details = 'msg'

        self.r2 = Result([])
        self.r2.ndiff_lines = 13
        self.r2.max_abs_err = 1.54364345e-12
        self.r2.max_abs_ln = 19
        self.r2.max_rel_err = 1.13288473436874
        self.r2.max_rel_ln = 23
        self.r2.details = 'msg'

    def test_passed_within_tols(self):
        self.r1.fatal_error = True
        self.r1.success = False
        isok, status, msg = self.r1.passed_within_tols(1000, 1000, 1000)
        self.assertFalse(isok)
        self.assertEqual(status, 'failed')

    def test_passed_within_tols_success(self):
        self.r1.fatal_error = False
        self.r1.success = True
        isok, status, msg = self.r1.passed_within_tols(0., 0., 0.)
        self.assertTrue(isok)
        self.assertEqual(status, 'succeeded')

    def test_passed_within_tols_passed(self):
        self.r1.fatal_error = False
        self.r1.success = False
        isok, status, msg = self.r1.passed_within_tols(15, 2e-5, 2.3e-6)
        self.assertTrue(isok)
        self.assertEqual(status, 'passed')

    def test_passed_within_tols_fail_lncount(self):
        self.r1.fatal_error = False
        self.r1.success = False
        isok, status, msg = self.r1.passed_within_tols(10, 1, 1)
        self.assertFalse(isok)
        self.assertEqual(status, 'failed')

    def test_passed_within_tols_fail_rel_err(self):
        self.r1.fatal_error = False
        self.r1.success = False
        isok, status, msg = self.r1.passed_within_tols(100, 1, 1.0e-15)
        self.assertFalse(isok)
        self.assertEqual(status, 'failed')

    def test_passed_within_tols_fail_abs_err(self):
        self.r1.fatal_error = False
        self.r1.success = False
        isok, status, msg = self.r1.passed_within_tols(100, 1.0e-15, 1)
        self.assertFalse(isok)
        self.assertEqual(status, 'failed')

    def test_get_summary_fail(self):
        self.r2.fatal_error = True
        self.assertEqual(self.r2.get_summary(), 'Summary: fldiff fatal error.\n')

    def test_get_summary_success(self):
        self.r2.fatal_error = False
        self.r2.success = True
        self.assertEqual(self.r2.get_summary(), 'Summary: no significant difference has been found')

    def test_get_summary_pass_with_label(self):
        self.r2.success = False
        self.r2.label = 'xyz'

        self.assertEqual(
            self.r2.get_summary(),
            'Summary xyz: different lines=13, max abs_diff=1.544e-12 (l.19), max rel_diff=1.133e+00 (l.23).'
        )

    def test_get_summary_pass_without_label(self):
        self.r2.success = False
        self.r2.label = None

        self.assertEqual(
            self.r2.get_summary(),
            'Summary: different lines=13, max abs_diff=1.544e-12 (l.19), max rel_diff=1.133e+00 (l.23).'
        )

    def test_analyse_no_difference(self):
        r1 = Result([])
        self.assertEqual(r1.max_abs_err, 0.0)
        self.assertEqual(r1.max_rel_err, 0.0)
        self.assertEqual(r1.ndiff_lines, 0)
        self.assertTrue(r1.success)
        self.assertFalse(r1.fatal_error)

    def test_analyse_forced(self):
        r1 = Result([fldiff.ForcedDifference(1, 1, 'l1\n', 'l2\n')])
        self.assertEqual(r1.max_abs_err, 0.0)
        self.assertEqual(r1.max_rel_err, 0.0)
        self.assertEqual(r1.ndiff_lines, 0)
        self.assertTrue(r1.success)
        self.assertFalse(r1.fatal_error)

    def test_analyse_fatal(self):
        r1 = Result([fldiff.LineCountDifference('f1', 'f2')])
        self.assertFalse(r1.success)
        self.assertTrue(r1.fatal_error)

        r1 = Result([fldiff.MetaCharDifference(1, 1, ' ', '%')])
        self.assertFalse(r1.success)
        self.assertTrue(r1.fatal_error)

    def test_analyse_float(self):
        r1 = Result([fldiff.FloatDifference(0, 1, 'l1\n', 'l2\n', 1.13e-2, 2.23e-3)])
        self.assertFalse(r1.success)
        self.assertFalse(r1.fatal_error)
        self.assertEqual(r1.max_abs_err, 1.13e-2)
        self.assertEqual(r1.max_abs_ln, 1)
        self.assertEqual(r1.max_rel_err, 2.23e-3)
        self.assertEqual(r1.max_rel_ln, 1)
        self.assertEqual(r1.ndiff_lines, 1)

    def test_analyse_text(self):
        r1 = Result([fldiff.TextDifference(1, 2, 'l1\n', 'l2\n')])
        self.assertFalse(r1.success)
        self.assertFalse(r1.fatal_error)
        self.assertEqual(r1.max_abs_err, 0.0)
        self.assertEqual(r1.max_rel_err, 0.0)
        self.assertEqual(r1.ndiff_lines, 1)

    def test_analyse_text_silent(self):
        r1 = Result([fldiff.TextDifference(1, 2, 'l1\n', 'l2\n', silent=True)])
        self.assertTrue(r1.success)
        self.assertFalse(r1.fatal_error)
        self.assertEqual(r1.max_abs_err, 0.0)
        self.assertEqual(r1.max_rel_err, 0.0)
        self.assertEqual(r1.ndiff_lines, 0)

    def test_analyse_two_diff_one_line(self):
        r1 = Result([
            fldiff.TextDifference(1, 2, 'l1\n', 'l2\n'),
            fldiff.FloatDifference(1, 2, 'l1\n', 'l2\n', 0.0, 0.0),

        ])
        self.assertEqual(r1.ndiff_lines, 1)


class TestDiffer(unittest.TestCase):
    lines1 = [
        ' Here are some regular lines of numbers : 0.4546\n',
        ' 5.8787 44.537e+056\n',
        ' .7856 5.0\n',
        '- This lines should be ignored\n',
        '+ This line should always appear in the output but is not counted as errorneous\n',
        ': This line is always compared as characters 0.457896321545\n',
        '. This one too but it will never be counted a errorneous 2.3735435434354364\n',
        '% This line have a fixed tolerance of 1.01e-2: 2.043643 5.5473684\n',
        'P P and , lines are handled according to parameters as + or - lines\n'
    ]

    lines2 = [  # numerical differences
        ' Here are some regular lines of numbers : 0.4546\n',
        ' 5.8707 44.535e+056\n',
        ' .7856 8.0\n',
        '- This lines should be ignored\n',
        '+ This line should always appear in the output but is not counted as errorneous\n',
        ': This line is always compared as characters 0.457896321545\n',
        '. This one too but it will never be counted a errorneous 2.3735435434354364\n',
        '% This line have a fixed tolerance of 1.01e-2: 2.043644 5.5\n',
        'P P and , lines are handled according to parameters as + or - lines\n'
    ]

    lines3 = [  # text differences
        ' Here are some regular lines of floats : 0.4546\n',
        ' 5.8787 44.537e+056\n',
        ' .7856 5.0\n',
        '- Should this line be ignored ?\n',
        '+ This line should always appear in the output but is not counted as errorneous\n',
        ': This line is always compared as characters 0.457896321546\n',
        '. This one too but it will never be counted a errorneous 2.3735435434354360\n',
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
        '+ This line should always appear in the output but is not counted as errorneous\n',
        '- Ignored lines are everywhere\n',
        ': This line is always compared as characters 0.457896321545\n',
        '. This one too but it will never be counted a errorneous 2.3735435434354364\n',
        '% This line have a fixed tolerance of 1.01e-2: 2.043643 5.5473684\n',
        'P P and , lines are handled according to parameters as + or - lines\n'
    ]

    lines5 = [  # line number differences significant
        ' Here are some regular lines of numbers : 0.4546\n',
        ' 78.73687 98.5763\n',
        ' 5.8787 44.537e+056\n',
        ' .7856 5.0\n',
        '- This lines should be ignored\n',
        '+ This line should always appear in the output but is not counted as errorneous\n',
        ': This line is always compared as characters 0.457896321545\n',
        '. This one too but it will never be counted a errorneous 2.3735435434354364\n',
        '% This line have a fixed tolerance of 1.01e-2: 2.043643 5.5473684\n',
        'P P and , lines are handled according to parameters as + or - lines\n'
    ]

    def test_default_parameters(self):
        diff = Differ()
        self.assertEqual(diff.options['tolerance_abs'], 1.01e-10)
        self.assertEqual(diff.options['tolerance_rel'], 1.01e-10)
        self.assertTrue(diff.options['ignore'])
        self.assertTrue(diff.options['ignoreP'])

    def test_get_metachar(self):
        class Dummy(Differ):
            def __init__(self, ignore=True, ignoreP=True):
                self.xml_mode = False
                self.options = {
                    'ignore': ignore,
                    'ignoreP': ignoreP,
                }
        self.assertEqual(Differ._Differ__get_metachar(Dummy(), '- a line\n'), '-')
        self.assertEqual(Differ._Differ__get_metachar(Dummy(), '+ a line\n'), '+')
        self.assertEqual(Differ._Differ__get_metachar(Dummy(), '. a line\n'), '.')

        self.assertEqual(Differ._Differ__get_metachar(Dummy(), ', a line\n'), '-')
        self.assertEqual(Differ._Differ__get_metachar(Dummy(False, False), ', a line\n'), '+')

        self.assertEqual(Differ._Differ__get_metachar(Dummy(), 'P a line\n'), '-')
        self.assertEqual(Differ._Differ__get_metachar(Dummy(False, False), 'P a line\n'), '+')

        self.assertEqual(Differ._Differ__get_metachar(Dummy(), '\n'), ' ')
        self.assertEqual(Differ._Differ__get_metachar(Dummy(), ''), ' ')

    def test_clean(self):
        class Dummy(Differ):
            def __init__(self, ignore=True, ignoreP=True):
                self.xml_mode = False
                self.options = {
                    'ignore': ignore,
                    'ignoreP': ignoreP,
                }
        dummy = Dummy()
        self.assertEqual(dummy._Differ__clean([
            '- an ignored line\n',
            ' a valid line\n',
            '% another valid line\n',
            'P another ignored line\n'
        ]), ([
            (1, ' a valid line\n'),
            (2, '% another valid line\n'),
        ], [
            (0, '- an ignored line\n'),
            (3, 'P another ignored line\n')
        ]))

    def test_diff_lines_same(self):
        diff = Differ()
        differences = diff._Differ__diff_lines(self.lines1, self.lines1)
        self.assertEqual(len(differences), 1)
        self.assertIsInstance(differences[0], fldiff.ForcedDifference)

    def test_diff_lines_float(self):
        diff = Differ()
        differences = diff._Differ__diff_lines(self.lines1, self.lines2)
        self.assertEqual(len(differences), 5)
        d3 = differences.pop(3)

        self.assertIsInstance(d3, fldiff.ForcedDifference)
        self.assertTrue(all(isinstance(d, fldiff.FloatDifference) for d in differences))

    def test_diff_lines_text(self):
        diff = Differ()
        differences = diff._Differ__diff_lines(self.lines1, self.lines3)
        self.assertEqual(len(differences), 4)

        d2 = differences.pop(1)
        self.assertIsInstance(d2, fldiff.ForcedDifference)

        self.assertTrue(all(isinstance(d, fldiff.TextDifference) for d in differences))

        d4 = differences.pop(2)
        self.assertTrue(d4.silent)

        self.assertFalse(any(d.silent for d in differences))

    def test_diff_lines_number_not_significant(self):
        diff = Differ()
        differences = diff._Differ__diff_lines(self.lines1, self.lines4)
        self.assertEqual(len(differences), 1)
        self.assertIsInstance(differences[0], fldiff.ForcedDifference)

    def test_diff_lines_number_significant(self):
        diff = Differ()
        differences = diff._Differ__diff_lines(self.lines1, self.lines5)
        self.assertEqual(len(differences), 1)
        self.assertIsInstance(differences[0], fldiff.LineCountDifference)
        self.assertEqual(differences[0].more, 'file 2')

    def test_diff_lines_float_format(self):
        diff = Differ()
        differences = diff._Differ__diff_lines(
            [' .0007  564.5e-3  7000.0\n'],
            [' 7.0e-4  5.645D-1  7.0f3\n']
        )
        self.assertEqual(len(differences), 0)

    def test_diff_ignore_blanks(self):
        diff = Differ()
        differences = diff._Differ__diff_lines(
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
        )
        self.assertEqual(len(differences), 0)


if __name__ == '__main__':
    unittest.main()
