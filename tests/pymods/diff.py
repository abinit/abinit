#!/usr/bin/env python
""" Command line interface to difflib.py providing diffs in four formats:

* ndiff:    lists every line and highlights interline changes.
* context:  highlights clusters of changes in a before/after format.
* unified:  highlights clusters of changes in an inline format.
* html:     generates side by side comparison with change highlights.
"""
from __future__ import print_function, division, absolute_import, unicode_literals

import sys
import os
import time
import difflib
import argparse


def abinit_junk(line):
    return (
        line.startswith('-')
        or line.startswith('+')
        or line.startswith(',')
        or line.startswith('P')
        or line.startswith('%')
        or line.startswith(';')
        or line.startswith('.')
        or line.startswith('_')
        or line.isspace()
    )


def main():
    # Configure the option parser
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', action='store_true', default=False,
                        help='Produce a context format diff (default)')

    parser.add_argument('-u', action='store_true', default=False,
                        help='Produce a unified format diff')

    hlp = 'Produce HTML side by side diff (can use -c and -l in conjunction)'
    parser.add_argument('-m', action='store_true', default=False, help=hlp)

    hlp = 'Produce HTML table of side by side diff (can use -c and -l in conjunction)'
    parser.add_argument('-t', action='store_true', default=False, help=hlp)

    parser.add_argument('-n', action='store_true', default=False,
                        help='Produce a ndiff format diff')
    parser.add_argument('-l', '--lines', type=int, default=3,
                        help='Set number of context lines (default 3)')

    parser.add_argument('-f', '--file', type=str, default='',
                        help='Write diff to file FILE. stdout is used if not specified', metavar='FILE')

    parser.add_argument('-j', '--abinit-junk', action='store_true', default=False,
                        help='Use builtin Abinit output specific heuristic instead of builtin one.')

    parser.add_argument('fromfile', help='Reference file')
    parser.add_argument('tofile', help='Compared file')

    options = parser.parse_args()

    n = options.lines
    fromfile, tofile = options.fromfile, options.tofile

    # we're passing these as arguments to the diff function
    fromdate = time.ctime(os.stat(fromfile).st_mtime)
    todate = time.ctime(os.stat(tofile).st_mtime)
    fromlines = open(fromfile, 'U').readlines()
    tolines = open(tofile, 'U').readlines()

    if options.abinit_junk:
        junk = abinit_junk
    else:
        junk = None

    if options.u:
        diff = difflib.unified_diff(fromlines, tolines, fromfile, tofile,
                                    fromdate, todate, n=n)
    elif options.n:
        diff = difflib.ndiff(fromlines, tolines, linejunk=junk)

    elif options.m:
        diff = difflib.HtmlDiff(linejunk=junk).make_file(
            fromlines, tolines, fromfile, tofile, context=options.c,
            numlines=n
        )
    elif options.t:
        diff = difflib.HtmlDiff(linejunk=junk).make_table(
            fromlines, tolines, fromfile, tofile, context=options.c,
            numlines=n
        )
    else:
        diff = difflib.context_diff(fromlines, tolines, fromfile, tofile,
                                    fromdate, todate, n=n)

    # writelines because diff is a generator
    if options.file:
        fh = open(options.file, 'w')
        fh.writelines(diff)
        fh.close()
    else:
        sys.stdout.writelines(diff)

    sys.exit(0)


if __name__ == '__main__':
    main()
