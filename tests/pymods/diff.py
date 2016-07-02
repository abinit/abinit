#!/usr/bin/env python
""" Command line interface to difflib.py providing diffs in four formats:

* ndiff:    lists every line and highlights interline changes.
* context:  highlights clusters of changes in a before/after format.
* unified:  highlights clusters of changes in an inline format.
* html:     generates side by side comparison with change highlights.
"""
from __future__ import print_function, division, absolute_import #, unicode_literals

import sys, os, time, difflib, optparse

def main():
     # Configure the option parser
    usage = "usage: %prog [options] fromfile tofile"
    parser = optparse.OptionParser(usage)

    parser.add_option("-c", action="store_true", default=False,
                      help='Produce a context format diff (default)')

    parser.add_option("-u", action="store_true", default=False,
                      help='Produce a unified format diff')

    hlp = 'Produce HTML side by side diff (can use -c and -l in conjunction)'
    parser.add_option("-m", action="store_true", default=False, help=hlp)

    hlp = 'Produce HTML table of side by side diff (can use -c and -l in conjunction)'
    parser.add_option("-t", action="store_true", default=False, help=hlp)

    parser.add_option("-n", action="store_true", default=False,
                      help='Produce a ndiff format diff')
    parser.add_option("-l", "--lines", type="int", default=3,
                      help='Set number of context lines (default 3)')

    parser.add_option("-f", "--file", type="string", default="",
                      help="Write diff to file FILE. stdout is used if not specified", metavar="FILE")

    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    if len(args) != 2:
        parser.error("need to specify both a fromfile and tofile")

    n = options.lines
    fromfile, tofile = args # as specified in the usage string

    # we're passing these as arguments to the diff function
    fromdate = time.ctime(os.stat(fromfile).st_mtime)
    todate = time.ctime(os.stat(tofile).st_mtime)
    fromlines = open(fromfile, 'U').readlines()
    tolines = open(tofile, 'U').readlines()

    if options.u:
        diff = difflib.unified_diff(fromlines, tolines, fromfile, tofile,
                                    fromdate, todate, n=n)
    elif options.n:
        diff = difflib.ndiff(fromlines, tolines)

    elif options.m:
        diff = difflib.HtmlDiff().make_file(fromlines, tolines, fromfile,
                                            tofile, context=options.c,
                                            numlines=n)
    elif options.t:
        diff = difflib.HtmlDiff().make_table(fromlines, tolines, fromfile,
                                             tofile, context=options.c,
                                             numlines=n)
    else:
        diff = difflib.context_diff(fromlines, tolines, fromfile, tofile,
                                    fromdate, todate, n=n)

    # writelines because diff is a generator
    if options.file:
        fh = open(options.file, "w")
        fh.writelines(diff)
        fh.close()
    else:
        sys.stdout.writelines(diff)

    sys.exit(0)

if __name__ == '__main__':
    main()
