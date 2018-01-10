#!/usr/bin/env python
from __future__ import print_function, division, absolute_import #, unicode_literals

import sys
import abc
import difflib
import tempfile
import functools

from collections import OrderedDict, Counter, namedtuple
from pprint import pprint

# Object: compare 2 output files from ABINIT line by line with arithmetic
# comparisons of floating point substrings
#
# Usage: fldiff [-context] [ -ignore | -include ] [ -ignoreP | -includeP ] [ -easy | -medium | -ridiculous ] file1 file2 [label]
#
# The first character of each line in both files indicates the mode of comparison;
# these first characters MUST be identical in corresponding lines.
# By default, a floating point comparison with a tolerance of 1.01d-10 is done
# on all floating point strings and a character comparison is done on all other
# strings, disregarding multiple spaces.
# In order for two numbers to be declared different, BOTH the absolute difference
# and the relative difference (difference divided by the sum of absolute values)
# must be bigger than the tolerance.
#
# Some special characters at the beginning of lines require a different handling:
# -	mark lines as same regardless to their content (i. e. ignore lines)
#     (can be be present in the 2 files or not, but must begin with -)
# _ mark lines as same regardless to their content
#     (must be present in the 2 files, but can begin with _ in only one of them)
# +	mark lines as different regardless to their content
# ,	handle as + with -include option and as - with -ignore option
# P	handle as + with -includeP option and as - with -ignoreP option
# %	floating point comparisons are done with a tolerance of 1.01e-2
# ;	floating point comparisons are done irrespective of signs
# :	ignore floating point numbers and do a characters comparison
# .	do a characters comparison, but do not count this line in the Summary
#
# Both files should have the same number of non - starting lines.

# With -context option, save character strings for context and print it
# with line number when floating difference is found.
#
# The -ignore and -include options affects the treatment of the ',' 
# special character in the first column (see above) 
#
# The -ignoreP and -includeP options affects the treatment of the 'P'
# special character in the first column (see above)
#
# If -ridiculous   is specified, the default tolerance is set to 1.01e-2 
# If -easy   is specified, the default tolerance is set to 1.01e-5 
# If -medium is specified, the default tolerance is set to 1.01e-8 
# These modifications do not apply to the tolerance determined by the
# '%',and '.' first-column special signs
#
# If "label" is specified, it will appear in the last, summary line.

# Helper functions.

#def iflat(iterables):
#    """
#    Iterator over all elements of a nested iterable. It's recursive!
#
#    >>> list(iflat([[0], [1,2, [3,4]]]))
#    [0, 1, 2, 3, 4]
#    """
#    for item in iterables:
#        if not hasattr(item, "__iter__"):
#            yield item
#        else:
#            # iterable object.
#            for it in iflat(item):
#                yield it

def dict2str(d):
    return ", ".join("%s=%s" % (k, v) for k, v in d.items())

def isnumber(string):
    """
    Returns True if string is a number
    """
    try:
        float(string)
        return True
    except ValueError:
        return False


def isint(string):
    """
    True if string represents an integer.

    >>> isint("1")
    True
    >>> isint("3.")
    False
    """
    try:
        return int(string) == float(string) and "." not in string
    except ValueError:
        return False


def num_sigdigits(string):
    """
    Returns the number of significant digits in a string representing 
    a floating point number in the form 1.234 
    (exponential format is not supported here)

    >>> num_sigdigits("3. ")
    0
    >>> num_sigdigits("-3.27")
    2
    """
    string = string.rstrip()
    idx = string.find(".")
    if idx == -1 or idx == len(string): return 0
    digits = string[idx+1:]
    assert all(c.isdigit() for c in digits)
    return len(digits)


class AttrDict(dict):
    """Allows to access dict keys as obj.foo in addition to the traditional way obj['foo']"""
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


class DiffOptions(AttrDict):
    """
    Dictionary with the options used to compare two files.
    """
    # Default values
    _DEFAULTS = dict(
        tiny=1e-12,    # Absolute difference
        tolnlines=0,   # Tolerance on the number of different lines.
        tolabs=1e-5,   # Tolerance of the absolute error.
        tolrel=1e-2,   # Tolerance of the relative error.
        verbose=0,     # Verbosity level.
    )

    def __init__(self, *args, **kwargs):
        # Ignore keys that are not in _DEFAULTS
        kwargs = {k: v for k,v in kwargs.items() if k in self._DEFAULTS}
        super(DiffOptions, self).__init__(*args, **kwargs)

        # Add defaults if not specified.
        for k, v in self._DEFAULTS.items():
            if k not in self: self[k] = v


class Token(object):
    """
    Base class for tokens.

    A Token provides a diff method.

    Attributes::
        lineno:
        colpos:
        string:
        value:
        dtype:
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, string):
        self.string = string

    @classmethod
    def from_string(cls, s):
        if isnumber(s):
            if isint(s): 
                return IntegerToken(s)
            else:
                return FloatToken(s)
        else:
            return StringToken(s)

    def __str__(self):
        return self.string

    def diff_class(self, other):
        """
        Test if two tokens belong to the same class. 
        Returns `Diff` object that evalues to True if the class differs.
        """
        cls1, cls2 = self.__class__ , other.__class__
        msg = "Different class: %s != %s" % (cls1, cls2) if cls1 != cls2 else ""
        diff_codes = "e" if msg else ""
        return Diff(self, other, diff_codes, msg)

    @property
    def is_numeric(self):
        """True if self represents a number."""
        return isinstance(self, NumericToken)

    @abc.abstractproperty
    def value(self):
        """The representation of the token"""

    @abc.abstractmethod
    def diff(self, other, diff_opts):
        """
        Compare self with other diff using the options given in diff_opts.
        Returns `Diff` instance.
        """


class StringToken(Token):
    """A token that represents a string."""
    dtype = "s"

    def diff(self, other, diff_opts):
        diff = self.diff_class(other)
        if diff: return diff

        # character comparison (disregarding multiple spaces is not needed)
        msg = "Different strings: '%s' != '%s'" % (self.value, other.value) if self.value != other.value else ""
        diff_codes = "s" if msg else ""
        return Diff(self, other, diff_codes, msg)

    @property
    def value(self):
        return self.string


class NumericToken(Token):
    """Base class for numeric tokens."""
    __metaclass__ = abc.ABCMeta

    #@abc.abstractproperty
    #def mant_exp(self):
                          
    #@abc.abstractproperty
    #def eps(self):
    #    """
    #    By definition, epsilon is the smallest number such as 1 + eps != 1, so
    #    there should be exactly one ULP between 1 and 1 + eps
    #    """

    def diff(self, other, diff_opts):
        diff = self.diff_class(other)
        if diff: return diff

        abs_err = abs(self.value - other.value)

        if abs(self.value) > diff_opts.tiny:
            rel_err = abs_err / self.value
        else:
            rel_err = 0.0
            # TODO
            # http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
            # Even now our function isn't perfect. In general this function will behave poorly for numbers around zero. 
            # The positive number closest to zero and the negative number closest to zero are extremely close to each other,
            # yet this function will correctly calculate that they have a huge relative error of 2.0. 
            # If you want to count numbers near zero but of opposite sign as being equal then you need to add a maxAbsoluteError check also. 
            # The function would then return true if either the absoluteError or the relativeError were smaller than the maximums passed in. 
            # A typical value for this backup maxAbsoluteError would be very small - FLT_MAX or less, depending on whether the platform supports subnormals.
            # Slightly better AlmostEqual function - still not recommended
            #bool AlmostEqualRelativeOrAbsolute(float A, float B, float maxRelativeError, float maxAbsoluteError)
            #{
            #    if (fabs(A - B) < maxAbsoluteError)
            #        return true;
            #    float relativeError;
            #    if (fabs(B) > fabs(A))
            #        relativeError = fabs((A - B) / B);
            #    else
            #        relativeError = fabs((A - B) / A);
            #    if (relativeError <= maxRelativeError)
            #        return true;
            #    return false;
            #}

        lines = []
        app = lines.append

        diff_codes = ""
        if abs_err > diff_opts.tolabs:
             app("abs_err: %s > %s" % (abs_err, diff_opts.tolabs))
             diff_codes += "a"

        if rel_err > diff_opts.tolrel:
            app("rel_err: %s > %s" % (rel_err, diff_opts.tolrel))
            diff_codes += "r"

        msg = "; ".join(lines) if lines else ""
        return Diff(self, other, diff_codes, msg, abs_err=abs_err, rel_err=rel_err)


class IntegerToken(NumericToken):
    """A token that represents an integer."""
    dtype = "i"

    @property
    def value(self):
        return int(self.string)

    @property
    def eps(self):
        return 0

    #@property
    #def mant_exp(self):
    #    return self.value, 0


class FloatToken(NumericToken):
    """A token that represents a float."""
    dtype = "f"

    @property
    def value(self):
        return float(self.string)

    @property
    def has_eform(self):
        """True if the number is in scientific notation."""
        return "E" in self.string

    @property
    def eps(self):
        if self.has_eform:
            raise NotImplementedError
        else:
            nsdigs = num_sigdigits(self.string)
            return 10**-nsdigs

    #@property
    #def mant_exp(self):


@functools.total_ordering
class Diff(object):
    """
    Stores the results of a token comparison.

    The object evaluates to True in boolean context if the difference is significative.
    Diffs can be ordered 

    A diff has a code that is a string that can contain the following characaters:
        "s" for string difference
        "a" for absolute difference
        "r" for relative difference
        "e" if an error occured (diff of objects of different class
    """
    POSSIBLE_CODES = ["e", "r", "a", "s", ]

    def __init__(self, token1, token2, diff_codes, msg, **kwargs):
        self.token1, self.token2 = token1, token2
        self.codes = diff_codes
        self.msg = msg
        self.extra = kwargs

        # None if string
        self.abs_err = kwargs.get("abs_err", None)
        self.rel_err = kwargs.get("rel_err", None)
        if self.abs_err is not None: assert  self.rel_err is not None

    @classmethod
    def make_counter(cls):
        return OrderedDict([(k, 0) for k in cls.POSSIBLE_CODES])

    def __str__(self):
        if self.is_numeric:
            return str([str(self.token1), str(self.token2), self.msg])
        else:
            return self.msg

    def __bool__(self):
        return bool(self.msg)

    __nonzero__ = __bool__

    # We want to order diffs so that we can output the most important one that occurs on a line
    # Numerical differences are more important that string difference
    # When comparing two numerical differences, the relative error is more important than the absolute one. 
    def __eq__(self, other):
        left, right = (self.rel_err, self.abs_err), (other.rel_err, other.abs_err)
        return left == right

    def __gt__(self, other):
        """Python2 supports 2 > None but this syntax raises TypeError in py3k."""
        left, right = (self.rel_err, self.abs_err), (other.rel_err, other.abs_err)

        # Compute absolute values
        try:
            left = tuple(map(abs, left))
        except TypeError:
            # (None, None) > x is always False
            return False
        
        try:
            right = tuple(map(abs, right))
        except TypeError:
            # (num1, num2) > (None, None) is always True
            return False

        # Now we can safely compare two tuples of absolute values
        # TODO: Fix problem with strings
        return left > right

    @property
    def is_numeric(self):
        """True if we have a diff of numbers."""
        ans = self.abs_err is not None
        if ans: assert self.rel_err is not None
        return ans
    


#import re
#re.split('([;, ])', 'foo,bar spam;neggs')

import string as __string
#__punctuation = '!"#$%&\'()*,/:;<=>?@[\\]^_`{|}~'
__punctuation = '!"#%\'()*,/:;<=>@[\\]^{|}'
__replace_punctuation = __string.maketrans(__punctuation, ' ' * len(__punctuation))


def parse_numstrings(line):
    """
    Split a lines in a list of numeric and string tokens.
    """
    # Replace punctuation characters with blank spaces (except for ".")
    # Preserve the first column as it has a special meaning.
    line = line[0] + line[1:].translate(__replace_punctuation)
    return [Token.from_string(s) for s in line.split()]


def ignore_numbers(line):
    """Split a string in a list of string tokens, ignore numbers"""
    return [StringToken(s) for s in line.split()]


class Tokenizer(object):
    """
    The `Tokenizer` receives two strings, looks at the character in the first column (they must 
    be equal) and split the lines in a list of tokens depending on the algorithm specified by 
    the character in the first column
    """
    _HANDLERS = {
    #  "-":                         # -	mark lines as same regardless to their content (i. e. ignore lines)
    #                               #     (can be be present in the 2 files or not, but must begin with -)
    #  "_":                         # _ mark lines as same regardless to their content
                                    #     (must be present in the 2 files, but can begin with _ in only one of them)
    #  "+":                         # +	mark lines as different regardless to their content
      "+": ignore_numbers,          # +	mark lines as different regardless to their content
    #  ",":                         # ,	handle as + with -include option and as - with -ignore option
    #  "P":                         # P	handle as + with -includeP option and as - with -ignoreP option
    #  "%":                         # %	floating point comparisons are done with a tolerance of 1.01e-2
    #  ";":                         # ;	floating point comparisons are done irrespective of signs
      ":":  ignore_numbers,         # :	ignore floating point numbers and do a characters comparison
      ".":  ignore_numbers,         # .	do a characters comparison, but do not count this line in the Summary
    }

    def __init__(self, diff_opts):
        self.diff_opts = diff_opts

    def make_diff(self, line1, opos1, line2, opos2):
        """
        Return sorted(diff_list), exceptions
        """
        exceptions = []
        eapp = exceptions.append

        ch1, ch2 = line1[0], line2[0]
        if ch1 != ch2: 
            msg = "Different chars in first column"
            eapp(DiffLineException(msg))

        # Select the function that will tokenize the lines.
        func = self._HANDLERS.get(ch1, parse_numstrings)

        tokens1 = func(line1)
        tokens2 = func(line2)

        # If we don't have the same number of tokens in line
        # we restrict the comparison to the minimum set of tokens
        # and we register the exception.
        if len(tokens1) != len(tokens2):
            msg ="Different number of tokens"
            eapp(DifferentNumberOfTokens(msg))

        diff_list = []

        for tok1, tok2 in zip(tokens1, tokens2):
            diff = tok1.diff(tok2, self.diff_opts)
            if not diff: continue
            diff_list.append(diff)

        # Return **sorted** list of diffs and exceptions
        return sorted(diff_list, reverse=True), exceptions


class DiffException(Exception):
    """Base exception class."""


class DiffLineException(DiffException):
    """Exceptions triggered by line comparison."""


#class DiffTokenException(DiffException):
#    """Exceptions triggered by token comparison."""


class DiffNumberOfTokens(DiffException):
    """Line compasison cannot be done..."""


class FileLine(namedtuple("FileLine", "opos, string")):
    """
        opos:
            Original position in the file
        string:
            String with the line content
    """


class Matcher(object):
    def __init__(self, file1, file2):
        # Open files, raise immediately if files are not found.
        try:
            with open(file1, "r") as fh1: all_lines1 = fh1.readlines()
        except IOError:
            raise IOError("Cannot open the reference file %s. This should not happen!" % file1)
                                                                                                           
        try:
            with open(file2, "r") as fh2: all_lines2 = fh2.readlines()
        except IOError:
            raise IOError("Cannot open the output file %s. Likely a problem with the abinit run!" % file1)
        
        self.op2line1, self.ignored1, done1 = self._polish(all_lines1)
        self.op2line2, self.ignored2, done2 = self._polish(all_lines2)

        # Handle non-critical exceptions
        self.exceptions = []

        # Make sure that both files are not empty
        if not len(self.op2line1) or not len(self.op2line2):
            msg = ""
            if not len(self.op2line1): msg += "Output file %s is empty!\n" % file1
            if not len(self.op2line2): msg += "Output file %s is empty!\n" % file2
            self.exceptions.append(DiffException(msg))
                                                                                  
        # Here we check if the two files have a different number of lines.
        if len(self.op2line1) != len(self.op2line2):
            msg = "Files have different number of lines! Will try to do my best!"
            self.exceptions.append(DiffException(msg))
            # TODO Use diff algorith to match the strings

        if not done1 or not done2:
            msg = ""
            if not done1: msg += "Output file %s is not completed\n" % file1
            if not done1: msg += "Output file %s is not completed\n" % file2
            self.exceptions.append(DiffException(msg))

    def show(self):
        for fline1, fline2 in self:
            print("Ref@%d: " % fline1.opos + fline1.string, end="")
            print("New@%d: " % fline2.opos + fline2.string, end="")
            print("")

    def diffalgo(self):
        # Each line of a Differ delta begins with a two-letter code:
        # '- '	line unique to sequence 1
        # '+ '	line unique to sequence 2
        # '  '	line common to both sequences
        # '? '	line not present in either input sequence
        # See https://docs.python.org/2.7/library/difflib.html
        deltas = list(difflib.Differ.compare(text1, text2))

        class MyDiffer(difflib.Differ):
            def my_compare(self, a, b):
                """
                Compare two sequences of lines; generate the resulting delta.

                Each sequence must contain individual single-line strings ending with
                newlines. Such sequences can be obtained from the `readlines()` method
                of file-like objects.  The delta generated also consists of newline-
                terminated strings, ready to be printed as-is via the writeline()
                method of a file-like object.

                Example:

                >>> print ''.join(Differ().compare('one\ntwo\nthree\n'.splitlines(1),
                ...                                'ore\ntree\nemu\n'.splitlines(1))),
                - one
                ?  ^
                + ore
                ?  ^
                - two
                - three
                ?  -
                + tree
                + emu
                """
                cruncher = difflib.SequenceMatcher(self.linejunk, a, b)
                for tag, alo, ahi, blo, bhi in cruncher.get_opcodes():
                    if tag == 'replace':
                        g = self._fancy_replace(a, alo, ahi, b, blo, bhi)
                    elif tag == 'delete':
                        g = self._dump('-', a, alo, ahi)
                    elif tag == 'insert':
                        g = self._dump('+', b, blo, bhi)
                    elif tag == 'equal':
                        g = self._dump(' ', a, alo, ahi)
                    else:
                        raise ValueError, 'unknown tag %r' % (tag,)

                    for line in g:
                        yield line

        # Create list of tuples (data, line) with the lines that are not common
        #data = []
        #for l in delta:
        #    code, l = l[:2], l[2:]
        #    if code == "  ": continue
        #    data.append(code, l)

        ## Create list of lines.
        #left, right = [], []
        #for code, line in reversed(code):
        #    if code == "? ":

        ## Reverse
        #left, right = list(reversed(left)), list(reversed(left))

    @staticmethod
    def _polish(lines):
        """
        Returns (op2line, ignores, completed)
        
            op2line: OrderedDict mapping opos --> line.
                     where opos is the line number in the **original** file.
            ignored: List of lines that are ignored
            completed: True if the output file is completed
        """
        op2line, ignored = OrderedDict(), OrderedDict()
        for opos, line in enumerate(lines):
            # Ignore lines that start with "-" or empty lines.
            #if line.startswith("-"): continue
            if line.startswith("-") or not line.strip(): 
                ignored[opos] = line
            else:
                op2line[opos] = line

        # TODO Check that the calculation is completed.
        completed = True

        return op2line, ignored, completed

    def __iter__(self):
        for (op1, s1), (op2, s2) in zip(self.op2line1.items(), self.op2line2.items()):
            yield FileLine(opos=op1, string=s1), FileLine(opos=op2, string=s2)

    def get_line1(self, opos):
        return self.op2line1[opos]

    def get_line2(self, opos):
        return self.op2line2[opos]


class FileDiffer(object):
    def __init__(self, diff_opts):
        self.diff_opts = diff_opts

    def analyze(self, file1, file2):
        """Main entry point for client code."""
        self.file1, self.file2 = file1, file2

        self.diffs_map = OrderedDict()
        self.exceptions = []

        try:
            self.matcher = matcher = Matcher(file1, file2)
            #matcher.show()
        except IOError:
            raise 

        # Add exceptions encountered by matcher.
        self.exceptions.extend(matcher.exceptions)

        # Number of different lines.
        nlines = 0
        tokenizer = Tokenizer(self.diff_opts)

        for fline1, fline2 in matcher:
            line1, opos1 = fline1.string, fline1.opos
            line2, opos2 = fline2.string, fline1.opos

            diffs, excs = tokenizer.make_diff(line1, opos1, line2, opos2)
            
            if excs:
                self.exceptions.extend(excs)

            if diffs: 
                # Store the line in an OrderedDict indexed by (opos1, opos2)
                key = (fline1.opos, fline2.opos)
                self.diffs_map[key] = diffs
                nlines += 1

        return 1 if self.exceptions or nlines else 0

    def nlines_totnumdiffs(self):
        """Return the number of lines with significant differences and the total number of diffs"""
        nlines, tot_numdiffs = 0, 0

        count = Diff.make_counter()

        for (pos1, pos2), diff_list in self.diffs_map.items():
            if any(diff_list): nlines += 1
            sig_diffs = [d for d in diff_list if d]
            tot_numdiffs += len(sig_diffs)
            for diff in sig_diffs:
                for c in diff.codes: count[c] += 1

        return nlines, tot_numdiffs, count

    def get_max_absdiff(self):
        return self._get_amax("abs_err")

    def get_max_reldiff(self):
        return self._get_amax("rel_err")

    def _get_amax(self, key):
        # Build flat list of significant numerical diffs (i.e. diffs that evaluate to True)
        # Use standard decorate and sort algorithm to get the position of the max.
        items = []
        for i, (_, diff_list) in enumerate(self.diffs_map.items()):
            items.extend([(d, i) for d in diff_list if d and d.is_numeric])
                                                                                            
        # Sort numerical diffs according to |d.key|
        abserr_sort = sorted(items, key=lambda t: abs(getattr(t[0], key)), reverse=True)
        #pprint([str(s) for s in abserr_sort])

        if len(abserr_sort): 
            # Return diff, index
            return abserr_sort[0] 
        else:
            # Return None if list is empty. 
            # This happens when there's no significant difference. 
            return None, None

    def write_results(self, out=sys.stdout, verbose=0):
        """
        Print a summary of the comparison.
        """
        summary = ["< %s\n> %s\n" % (self.file1, self.file2)]
        app = summary.append

        app("Number of exceptions: %d" % len(self.exceptions))
        for i, exc in enumerate(self.exceptions):
            app("Exception %d): %s" % (i, str(exc)))

        nlines, tot_numdiffs, count = self.nlines_totnumdiffs()
        app("Number of different lines: %d" % nlines)
        app("Total number of significant differences: %d" % tot_numdiffs)
        app("Diff counters: %s" % dict2str(count))

        # Find the max absolute error and the max relative error.
        max_absdiff, idx_amax = self.get_max_absdiff()
        max_reldiff, idx_rmax = self.get_max_reldiff()

        # Print info on the max absolute error and the max relative error.
        #if max_absdiff is not None:
        #if max_reldiff is not None:
        app("Max absolute difference at %s: %s" % (idx_amax, max_absdiff))
        app("Max relative difference at %s: %s" % (idx_rmax, max_reldiff))

        summary = "\n".join(summary)

        out.write(summary)
        html = HtmlDiffFile()
        html.make_summary(summary)

        for i, ((pos1, pos2), diff_list) in enumerate(self.diffs_map.items()):
            # No printout if there's not significant diff for this line.
            if not any(diff_list): continue

            # the header has the format: diff_index, (pos1, pos2), diff counter
            dcount = dict2str(Counter("".join(d.codes for d in diff_list)))
            header = "%d, (%d, %d), %s" % (i, pos1, pos2, dcount)

            # Write lines.
            s1 = "< " + self.matcher.get_line1(pos1)
            s2 = "> " + self.matcher.get_line2(pos2)

            out.write(header + "\n")
            out.write(s1)
            out.write(s2)

            if not options.verbose:
                # Print the most important difference
                print(diff_list[0])
            else:
                print(str(diff_list))

            html.add_lines_and_diff(i, header, s1, s2, diff_list, verbose=verbose)

        out.write("\n")

        #html.write(os.path.abspath("foo.html"))
        html.browse()


def browse_localpath(filepath):
    """
    Open filepath in the default browser.

    .. warning: This code is not portable since we should pass a url.
    """
    import webbrowser
    if not filepath.startswith("file://"): filepath = "file://" + filepath

    try:
        webbrowser.open(filepath)
    except webbrowser.Error as exc:
        # Warn the user and ignore the exception.
        warnings.warn(str(exc))


class HtmlDiffFile(object):
    """
    This object produces a HTML file with a summary of the comparison.
    """
    def __init__(self):
        # For html see also http://bonrouge.com/~togglit
        js = """
<script type="text/javascript">
function toggle(obj) {
          var obj=document.getElementById(obj);
          if (obj.style.display == "block") obj.style.display = "none";
          else obj.style.display = "block";
}
</script>"""
        self.text = [js]
        self.make_help()

    def make_help(self):
        help_id = "-1"
        s = """
        <a href="javascript: void(0);" onClick="toggle(%(help_id)s)"> 
        Help <br>
        </a>
        <div id="%(help_id)s" style="display:none;">
        <hr>
        The header has the format: diff index, (pos1, pos2), counter <br>

        where:<br>

          index is the sequential index of the significant difference<br>
          (pos1, pos2) gives the position in the two files<br>
          counter counts the occurrence of each type of differences per line<br>

        Codes:<br>
        "s" for string <br>
        "a" for absolute error <br>
        "r" for relative error <br>
        "e" for exceptions <br>
        <hr>
        </div>
        """ % locals()
                                                                           
        self.text.append(s)

    def make_summary(self, summary):
        summary = self._str2html(summary)
        s = """
            %(summary)s <br> 
        """ % locals()

        #<a name="max_abserr"></a> 
        #Click <a href="#max_abserr">here</a> to read chapter 4. 

        self.text.append(s)

    def add_lines_and_diff(self, doc_id, header, s1, s2, diff_list, verbose=0):
        if not verbose:
            info = str(diff_list[0])
        else:
            info = str(diff_list)

        s = """
        <a href="javascript: void(0);" onClick="toggle(%(doc_id)s)"> 
        %(header)s <br>
        </a>
        %(s1)s <br> 
        %(s2)s <br>
        <div id="%(doc_id)s" style="display:none;"><hr>%(info)s<hr> </div>
        """ % locals()

        self.text.append(s)

    @staticmethod
    def _str2html(string, end="<br>"):
        """Returns a HTML string."""
        lines = string.splitlines()
        return "<br>".join(lines) + end

    def write(self, filepath):
        with open(filepath, "w") as fh:
            fh.write(" ".join(self.text))

    def browse(self):
        # Create temporary file, write data, open file in the browser 
        _, tmp_path = tempfile.mkstemp(text=True)

        self.write(tmp_path)
        browse_localpath(tmp_path)

# Unit tests
from unittest import TestCase

class TestTokens(TestCase):

    def test_base(self):
        """test compare methods"""
        atrue, afalse = self.assertTrue, self.assertFalse

        diff_opts = DiffOptions()

        s = Token.from_string("hello")
        afalse(s.is_numeric)

        i = Token.from_string("1")
        atrue(i.is_numeric and i.dtype=="i" and i.eps==0)
        atrue(s.diff(i, diff_opts))

        f = Token.from_string("-3.25")
        atrue(f.is_numeric and f.dtype=="f" and not f.has_eform and f.eps == 0.01)

        e = Token.from_string("-3.25E+00")
        atrue(e.is_numeric and f.dtype=="f" and e.has_eform)
        #atrue(e.eps == 0.01)

        # f and e are equal numbers
        afalse(e.diff(f, diff_opts))

        # Test the ordering of diffs
        diff2str = lambda s1, s2 : Token.from_string(s1).diff(Token.from_string(s2), diff_opts)

        d0 = diff2str("hello", "helli")
        d1 = diff2str("hello", "world")
        atrue(d0 < d1)

        d0 = Token.from_string("hello").diff(Token.from_string("world"))
        d1 = Token.from_string("1").diff(Token.from_string("10"))
        #diffs = [d0, d1]

        #sort_diffs = sorted(diffs)
        #atrue(diffs == sort_diffs)
        #assert 0


class TestTokenizer(TestCase):

    def test_parse_numstrings(self):
        aequal = self.assertEqual

        mapstr = lambda tokens : map(str, tokens)
        dtypes = lambda tokens : [t.dtype for t in tokens]

        tokens = parse_numstrings("P  mgfft = 16   mkmem=2 mpssoang=   3 mpw =188")
        aequal(mapstr(tokens), ['P', 'mgfft', '16', 'mkmem', '2', 'mpssoang', '3', 'mpw', '188'])
        aequal(dtypes(tokens), ['s', 's', 'i', 's', 'i', 's', 'i', 's', 'i'])

        tokens = parse_numstrings(" DATASET    2 : space group Fd -3 m (#227); Bravais cF (face-center cubic)")
        aequal(mapstr(tokens), ['DATASET', '2', 'space', 'group', 'Fd', '-3', 'm', '227', 'Bravais', 'cF', 'face-center', 'cubic'])
        aequal(dtypes(tokens), ['s', 'i', 's', 's', 's', 'i', 's', 'i', 's', 's', 's', 's'])

        tokens = parse_numstrings("  sigma(1 1)=  9.50962026E-05  sigma(3 2)=  0.00000000E+00")
        aequal(mapstr(tokens), ['sigma', '1', '1', '9.50962026E-05', 'sigma', '3', '2', '0.00000000E+00'])
        aequal(dtypes(tokens), ['s', 'i', 'i', 'f', 's', 'i', 'i', 'f'])

        tokens = parse_numstrings(". (10% accuracy) 1, (-1.E+01, 3.41) (kssform=1), [-5.2,2.2]")
        aequal(mapstr(tokens), ['.', '10', 'accuracy', '1', '-1.E+01', '3.41', 'kssform', '1', '-5.2', '2.2'])
        aequal(dtypes(tokens), ['s', 'i', 's', 'i', 'f', 'f', 's', 'i', 'f', 'f'])

        tokens = parse_numstrings(" 17:31:21 EDT 1994 12x  {-0.25, 1.); {-0.3, getscr/=0")
        aequal(mapstr(tokens), ['17', '31', '21', 'EDT', '1994', '12x', '-0.25', '1.', '-0.3', 'getscr', '0'])
        aequal(dtypes(tokens), ['i', 'i', 'i', 's', 'i', 's',  'f', 'f', 'f', 's', 'i'])


def main():
    def str_examples():
        examples = (
          "\n"
          "Usage example:\n\n" 
          "fldiff.py file1 file2    ==> Compare file1 with file2\n"
          "                                                     \n"
        )
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg:  sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    import argparse 
    parser = argparse.ArgumentParser(epilog=str_examples(),formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
                         help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('files', nargs=2, help='Files to be compared')

    # Parse the command line.
    try:
        options = parser.parse_args()
    except Exception:
        show_examples_and_exit(error_code=1)

    # Extract diff_opts from command line options.
    diff_opts = DiffOptions(vars(options))
    differ = FileDiffer(diff_opts)

    file1, file2 = options.files
    retcode = differ.analyze(file1, file2)
    if retcode != 0: differ.write_results()

    return retcode 


if __name__ == "__main__":
    sys.exit(main())
