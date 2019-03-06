'''
Object: compare 2 output files from ABINIT line by line with arithmetic
comparisons of floating point substrings

The first character of each line in both files indicates the mode of
comparison; these first characters MUST be identical in corresponding lines.
By default, a floating point comparison with a tolerance of 1.01d-10 is done on
all floating point strings and a character comparison is done on all other
strings, disregarding multiple spaces.  For compatibility reasons with
historical fldiff.pl a float is must contains a dot and at least one digit
after it to be recognized as float.  In order for two numbers to be declared
different, BOTH the absolute difference and the relative difference (difference
divided by the sum of absolute values) must be bigger than the tolerance.

Some special characters at the beginning of lines require a different handling:
-	mark lines as same regardless to their content (i. e. ignore lines)
(can be be present in the 2 files or not, but must begin with -) _ mark lines
as same regardless to their content (must be present in the 2 files, but can
begin with _ in only one of them) +	mark lines as different regardless to
their content ,	handle as + if ignore option is False and as - else P	handle
as + if ignoreP option is False and as - else %	floating point comparisons are
done with a tolerance of 1.01e-2 ;	floating point comparisons are done
irrespective of signs :	ignore floating point numbers and do a characters
comparison .	do a characters comparison, but do not count this line in the
Summary

Both files should have the same number of non - starting lines.

The ignore options affects the treatment of the ',' special character in the
first column (see above)

The ignoreP options affects the treatment of the 'P' special character in the
first column (see above)

The label option, if specified, is appended at the end of the summary

the tolerance option set the tolerance for comparision of floats, the default
is 1.01e-10.  This modifications do not apply to the tolerance determined by
the '%',and '.' first-column special signs.
'''

from __future__ import print_function, division, unicode_literals
import re
from math import floor

from .data_extractor import DataExtractor
from .yaml_tools import is_available as has_yaml

if has_yaml:
    from .yaml_tools.driver_test_conf import DriverTestConf as YDriverConf
    from .yaml_tools.tester import Tester as YTester
else:
    from .yaml_tools import Mock

    class YDriverConf(Mock):
        pass

    class YTester(Mock):
        pass

# Match floats. Minimal float is .0 for historical reasons.
# In consequence integers will be compared as strings
float_re = re.compile(r'([+-]?[0-9]*\.[0-9]+(?:[eEdDfF][+-]?[0-9]+)?)')


def norm_spaces(s):
    '''
        Normalize all blanks ( \\n\\r\\t).
    '''
    # the join/split technic remove all blanks and put one space between
    # non-blanks words
    return ' '.join(s.split())


def relative_truncate(f, n):
    '''
        >>> rel_truncate(1.8367387367, 2)
        1.83
        >>> rel_truncate(1.8367387367e-5, 5)
        1.83673e-05
        >>> rel_truncate(1.8367387367e+7, 4)
        18367000.0
    '''
    ten_n = 10.0**n
    if f == 0.0:
        return 0.0
    elif abs(f) >= 1:
        ten_p = 10.0
        while abs(f) > ten_p:
            ten_p *= 10.0
        fact = 10 * ten_n / ten_p
        return floor(f * fact) / fact
    else:
        ten_p = 0.1
        p = -1
        while abs(f) < ten_p:
            ten_p *= 0.1
            p -= 1
        fact = ten_n / ten_p
        return floor(f * fact) / fact


class ConstDict(object):
    '''
        Represent an immutable dict.
    '''
    def __init__(self, d=None):
        if d:
            self.__d = dict(d)
        else:
            self.__d = {}

    def __getitem__(self, key):
        return self.__d[key]

    def get_dict(self):
        return self.__d.copy()


default_options = ConstDict({
    'ignore': True,
    'ignoreP': True,
    'tolerance_abs': 1.01e-10,
    'tolerance_rel': 1.01e-10,
    'label': None,
    'use_fl': True,
    'use_yaml': True
})


class LineDifference(object):
    '''
        Base class for representing a difference.
    '''
    def __init__(self, p1, p2, l1, l2):
        self.lines = (p1 + 1, p2 + 1)
        if l1 == '' or l1[-1] not in '\n\r':
            l1 += '\n'
        self.content = (l1, l2)

    def __eq__(self, other):
        '''
            Implement the == test.
        '''
        return self.lines == other.lines and self.content == other.content

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        '''
            Default representation of difference is inspired by gnu diff tool.
        '''
        return (
            '{}\n'.format(*self.lines)
            + '< ' + self.content[0]
            + '> ' + self.content[1]
        )


class LineCountDifference(LineDifference):
    '''
        Represent a difference between line counts.
    '''
    def __init__(self, more, less, line_count=(0, 0)):
        '''
            more: the name of the file with more lines
            less: the name of the file with less lines
        '''
        LineDifference.__init__(self, 0, 0, '', '')
        self.more = more
        self.less = less
        self.line_count = line_count

    def __repr__(self):
        if self.line_count != (0, 0):
            return ('{} have more significant lines than {} ({} > {}).\n'
                    .format(self.more, self.less, *self.line_count))

        return '{} have more significant lines than {}.\n'.format(self.more,
                                                                  self.less)


class MetaCharDifference(LineDifference):
    '''
        Represent a difference between themetacharacters of lines.
    '''
    def __init__(self, p1, p2, m1, m2):
        LineDifference.__init__(self, p1, p2, '', '')
        self.metas = (m1, m2)

    def __repr__(self):
        return ('At line {} (in file 1), line {} (in file 2), the leading'
                ' characters where differents: {} and {}.\n').format(
                    *(self.lines + self.metas))


class FloatDifference(LineDifference):
    '''
        Represent a difference between floating point values.
    '''
    def __init__(self, p1, p2, line1, line2, abs_err, rel_err):
        LineDifference.__init__(self, p1, p2, line1, line2)
        self.abs_err = abs_err
        self.rel_err = rel_err


class TextDifference(LineDifference):
    '''
        Represent a difference between text parts of a lines.
    '''
    def __init__(self, p1, p2, line1, line2, silent=False):
        LineDifference.__init__(self, p1, p2, line1, line2)
        self.silent = silent


class ForcedDifference(LineDifference):
    '''
        A difference is arbitrarly declared.
    '''
    pass


class Result(object):
    '''
        Analyse and summarize the set of differences found by a diff.
    '''
    def __init__(self, fl_diff, yaml_diff, extra_info=[], label=None):
        '''
            differences is expected to be a list of Difference instances
        '''
        self.fl_diff = fl_diff
        self.yaml_diff = yaml_diff
        self.extra_info = extra_info
        self.fatal_error = False
        self.yaml_error = False
        self.success = True
        self.max_abs_err = 0.0
        self.max_rel_err = 0.0
        self.max_abs_ln = 0
        self.max_rel_ln = 0
        self.ndiff_lines = 0
        self.label = label

        self.details = self.__analyse()

    def __analyse(self):
        '''
            Analyse a difference list and extract summary informations and
            details.  Sumary informations are
            - self.max_abs_err: maximum absolute difference
            - self.max_rel_err: maximu relative difference
            - self.max_abs_ln: line number where the maximum absolute
              difference is reached for the first time
            - self.max_rel_ln: line number where the maximum relative
              difference is reached for the first time
            - self.ndiff_lines: number of lines flagged as different (excluding
              "silent" differences: line starting with '.' '+' and depending of
              Diff options ',' and 'P')
        '''
        details = []
        error_lines = set()


        if self.yaml_diff:
            details.append('# Start YAML based comparison report\n')
            self.yaml_error = True
        for diff in self.yaml_diff:
            self.success = False
            details.append(repr(diff) + '\n\n')

        if self.fl_diff:
            details.append('# Start legacy fldiff comparision report\n')
        for diff in self.fl_diff:
            if isinstance(diff, LineCountDifference) \
               or isinstance(diff, MetaCharDifference):
                self.fatal_error = True
                self.success = False
                details = str(diff)

            elif isinstance(diff, ForcedDifference):
                pass  # Silent differences: not counted as different line

            elif isinstance(diff, FloatDifference):
                if diff.lines[0] not in error_lines:
                    self.ndiff_lines += 1
                if diff.abs_err > self.max_abs_err:
                    self.max_abs_err = diff.abs_err
                    self.max_abs_ln = diff.lines[0]
                if diff.rel_err > self.max_rel_err:
                    self.max_rel_err = diff.rel_err
                    self.max_rel_ln = diff.lines[0]
                self.success = False

            elif isinstance(diff, TextDifference):
                if not diff.silent:
                    if diff.lines[0] not in error_lines:
                        self.ndiff_lines += 1
                    self.success = False

            else:  # any other Difference
                assert isinstance(diff, LineDifference), \
                        'Unknown type of Difference.'
                if diff.lines[0] not in error_lines:
                    self.ndiff_lines += 1
                self.success = False

            if self.fatal_error:
                break  # stop analysing in case of fatal error

            if diff.lines[0] not in error_lines:
                details.append(str(diff))
                error_lines.add(diff.lines[0])

        return details

    def get_summary(self):
        '''
            Return a textual summary of the diff.
        '''
        if self.yaml_error:
            summary = 'yaml_test errors.'
        elif self.fatal_error:
            summary = 'fldiff fatal error.'
        elif self.success:
            summary = 'no significant difference has been found.'
        else:
            summary = ('different lines={}, max abs_diff={:.3e} (l.{}),'
                       ' max rel_diff={:.3e} (l.{}).').format(
                self.ndiff_lines,
                self.max_abs_err,
                self.max_abs_ln,
                self.max_rel_err,
                self.max_rel_ln
            )
        if self.label is not None:
            summary = 'Summary ' + self.label + ': ' + summary
        else:
            summary = 'Summary: ' + summary

        return summary

    def dump_details(self, file=None):
        '''
            Either return a string describing all detected differences
            or write it into the given file (expected to be a writable stream).
        '''
        if file is None:
            return ('\n'.join(self.extra_info) + '\n' + ''.join(self.details)
                    + self.get_summary())
        else:
            file.write('\n'.join(self.extra_info) + '\n')
            file.writelines(self.details)
            file.write(self.get_summary() + '\n')
            return None

    def passed_within_tols(self, tolnlines, tolabs, tolrel):
        '''
            Check the result of the diff against the given tolerances.
        '''
        if self.yaml_error:
            status = 'failed'
            msg = 'yaml_test errors. First is:\n{}\n'.format(self.yaml_diff[0])
        elif self.fatal_error:
            status = 'failed'
            msg = 'fldiff fatal error:\n' + self.details
        elif self.success:
            status = 'succeeded'
            msg = 'succeeded'
        else:
            # truncate to prevent fldiff from printing 1.000 < 1.000
            # compatibility fix, this may be removed later
            abs_error = relative_truncate(self.max_abs_err, 3)
            rel_error = relative_truncate(self.max_rel_err, 3)
            ndiff_lines = self.ndiff_lines
            status = 'failed'
            fact = 1.0

            locs = locals()
            if ndiff_lines > tolnlines:
                msg = 'failed: erroneous lines {ndiff_lines} > {tolnlines}'
            elif abs_error > tolabs * fact and rel_error < tolrel:
                msg = 'failed: absolute error {abs_error:.4} > {tolabs}'
            elif rel_error > tolrel * fact and abs_error < tolabs:
                msg = 'failed: relative error {rel_error:.4} > {tolrel}'
            elif abs_error > tolabs * fact and rel_error > tolrel * fact:
                msg = ('failed: absolute error {abs_error:.4} > {tolabs},'
                       ' relative error {rel_error:.4} > {tolrel}')
            else:
                status = 'passed'
                msg = ('passed: absolute error {abs_error:.4} < {tolabs},'
                       ' relative error {rel_error:.4} < {tolrel}')

            msg = msg.format(**locs)

        isok = status in ('succeeded', 'passed')
        return isok, status, msg


class Differ(object):
    def __init__(self, yaml_test=None, **options):
        '''
            Init a differ with some parameters.
            Known parameters are:
                - ignore: bool (default True)
                - ignoreP: bool (default True)
                - tolerance: float (tolerance for both relative and absolute
                  difference)
                - tolerance_abs: float (default 1.01e-10)
                - tolerance_rel: float (default 1.01e-10)
                - label: str (default None)
                - use_yaml: bool (default True)
                - use_fl: bool (default True)
        '''
        self.xml_mode = False  # this is the first dirty fix.
        self.options = default_options.get_dict()
        self.options.update(options)
        if 'tolerance' in options:
            self.options['tolerance_abs'] = options['tolerance']
            self.options['tolerance_rel'] = options['tolerance']

        self.use_fl = self.options['use_fl']
        self.yaml_test = YDriverConf()
        self.use_yaml = False
        if has_yaml:
            self.use_yaml = self.options['use_yaml']
            if self.use_yaml:
                if yaml_test and 'file' in yaml_test and yaml_test['file']:
                    self.yaml_test = YDriverConf.from_file(yaml_test['file'])
                elif yaml_test and 'yaml' in yaml_test and yaml_test['yaml']:
                    self.yaml_test = YDriverConf(yaml_test['yaml'])
                else:
                    self.yaml_test = YDriverConf()

    def diff(self, file1, file2):
        '''
            Compute the diff of file 1 (reference) and file 2 (out)
            and return a Result instance.
        '''
        if file1.endswith('.xml'):
            self.xml_mode = True

        with open(file1, 'rt') as f:
            lines1 = f.readlines()

        with open(file2, 'rt') as f:
            lines2 = f.readlines()

        return Result(*self.diff_lines(lines1, lines2),
                      extra_info=(self.yaml_test.extra_info() if self.use_yaml
                                  else ()),
                      label=self.options['label'])

    def diff_lines(self, lines1, lines2):
        dext = DataExtractor(xml_mode=self.xml_mode,
                             ignore=self.options['ignore'],
                             ignoreP=self.options['ignoreP'])
        lines1, documents1, _ = dext.extract(lines1)
        lines2, documents2, _ = dext.extract(lines2)

        if self.use_fl:
            lines_differences = self.__diff_lines(lines1, lines2)
        else:
            lines_differences = []

        if self.use_yaml:
            doc_differences, _ = self.__test_doc(documents1, documents2)
        else:
            doc_differences = []

        return lines_differences, doc_differences

    def __test_doc(self, docs1, docs2):
        '''
            Compare docs2 to docs1 and apply tests on docs2.
        '''
        return YTester(docs1, docs2, self.yaml_test).run()

    def __diff_lines(self, lines1, lines2):
        '''
            Compute the effective comparision between two set of lines.
            LineCountDifference and MetaCharDifference are both fatal so
            they are returned alone if encountered.
        '''
        differences = []
        if len(lines1) > len(lines2):
            return [LineCountDifference('file 1', 'file 2',
                                        (len(lines1), len(lines2)))]
        elif len(lines1) < len(lines2):
            return [LineCountDifference('file 2', 'file 1',
                                        (len(lines2), len(lines1)))]
        else:
            for (i1, meta1, line1), (i2, meta2, line2) in zip(lines1, lines2):
                if meta1 != meta2:
                    if meta1 != '_' and meta2 != '_':
                        return [MetaCharDifference(i1, i2, meta1, meta2)]
                else:
                    if meta1 == '_':  # ignore these lines
                        pass

                    elif meta1 == '+':  # these lines are arbitrarly different
                        differences.append(ForcedDifference(
                            i1, i2, line1, line2
                        ))

                    elif meta1 == ':':  # do a character comparison
                        if norm_spaces(line1) != norm_spaces(line2):
                            differences.append(TextDifference(
                                i1, i2, line1, line2
                            ))

                    elif meta1 == '.':  # do a silent character comparison
                        if norm_spaces(line1) != norm_spaces(line2):
                            differences.append(TextDifference(
                                i1, i2, line1, line2, silent=True
                            ))

                    else:  # compare numerical values
                        is_float = False
                        splitted1 = float_re.split(line1)
                        splitted2 = float_re.split(line2)

                        # not the same number of floats on the line
                        if len(splitted1) != len(splitted2):
                            differences.append(TextDifference(i1, i2,
                                                              line1, line2))
                        else:
                            for v1, v2 in zip(splitted1, splitted2):
                                if is_float:
                                    tol = self.options['tolerance_abs']
                                    tolrel = self.options['tolerance_rel']
                                    f1 = float(v1.lower().replace('d', 'e')
                                               .replace('f', 'e'))
                                    f2 = float(v2.lower().replace('d', 'e')
                                               .replace('f', 'e'))

                                    if meta1 == ';':  # compare absolute values
                                        f1, f2 = abs(f1), abs(f2)
                                    elif meta1 == '%':  # force tolerance
                                        tol = 1.01e-2

                                    diff = abs(f1 - f2)
                                    if abs(f1 + f2) == 0.0:
                                        if diff > 0:
                                            diffrel = 1
                                        else:
                                            diffrel = 0.0
                                    else:
                                        diffrel = diff / (abs(f1) + abs(f2))

                                    if diff > tol and diffrel > tolrel:
                                        differences.append(
                                            FloatDifference(
                                                i1, i2, line1, line2,
                                                diff, diffrel
                                            )
                                        )

                                else:
                                    if norm_spaces(v1) != norm_spaces(v2):
                                        differences.append(TextDifference(
                                            i1, i2, line1, line2)
                                        )
                                # re.split always altern between separator and
                                # match (here floats)
                                is_float = not is_float
        return differences
