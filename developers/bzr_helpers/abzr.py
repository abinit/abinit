#!/usr/bin/env python
"""
This script automates the editing of a long list of source files, for example
when we have several conflicts after a bzr merge or when we want to edit 
several files containing a particular token.
"""
from __future__ import print_function

import sys
import os
import argparse 
import fnmatch
import subprocess 

def is_string(s):
    """True if s behaves like a string (duck typing test)."""
    try:
        dummy = s + " "
        return True

    except TypeError:
        return False

def list_strings(arg):
    """
    Always return a list of strings, given a string or list of strings as
    input.

    :Examples:

    >>> list_strings('A single string')
    ['A single string']

    >>> list_strings(['A single string in a list'])
    ['A single string in a list']

    >>> list_strings(['A','list','of','strings'])
    ['A', 'list', 'of', 'strings']
    """
    if is_string(arg):
        return [arg]
    else:
        return arg


class Editor(object):
    def __init__(self, editor=None):
        if editor is None:
            self.editor = os.getenv("EDITOR", "vi")
        else:
            self.editor = str(editor)

    def edit_file(self, fname):
        from subprocess import call
        retcode = call([self.editor, fname])
        if retcode != 0:
            import warnings
            warnings.warn("Error while trying to edit file: %s" % fname)
        return retcode

    def edit_files(self, fnames, ask_for_exit=True):
        exit_status = 0
        for idx, fname in enumerate(fnames):
            exit_status = self.edit_file(fname)
            if ask_for_exit and idx != len(fnames)-1 and self.user_wants_to_exit():
                break
        return exit_status

    @staticmethod
    def user_wants_to_exit():
        try:
            answer = raw_input("Do you want to continue [Y/n]")
        except EOFError:
            #print "exc"
            return True

        return answer.lower().strip() in ["n", "no"]


class WildCard(object):
    """
    This object provides an easy-to-use interface for
    filename matching with shell patterns (fnmatch).

    .. example:

    >>> w = WildCard("*.nc|*.pdf")
    >>> w.filter(["foo.nc", "bar.pdf", "hello.txt"])
    ['foo.nc', 'bar.pdf']

    >>> w.filter("foo.nc")
    ['foo.nc']
    """
    def __init__(self, wildcard, sep="|"):
        """
        Initializes a WildCard.

        Args:
            wildcard (str): String of tokens separated by sep. Each token
                represents a pattern.
            sep (str): Separator for shell patterns.
        """
        self.pats = ["*"]
        if wildcard:
            self.pats = wildcard.split(sep)

    def __str__(self):
        return "<%s, patterns = %s>" % (self.__class__.__name__, self.pats)

    def filter(self, names):
        """
        Returns a list with the names matching the pattern.
        """
        names = list_strings(names)

        fnames = []
        for f in names:
            for pat in self.pats:
                if fnmatch.fnmatch(f, pat):
                    fnames.append(f)

        return fnames

    def match(self, name):
        """
        Returns True if name matches one of the patterns.
        """
        for pat in self.pats:
            if fnmatch.fnmatch(name, pat):
                return True

        return False


def edit_conflicts():
    """
    Call bzr status, extract the paths of the files with conflicts 
    and edit them with the editor specified by $EDITOR 
    """
    # Run "bzr status" and get the output
    output = subprocess.check_output(["bzr", "status"])

    # Find the files with conflicts.
    sentinel = "  Text conflict in"
    conflicts = []
    for line in output.splitlines():
        if line.startswith(sentinel):
            conflicts.append(line[len(sentinel):].strip())

    # Call the editor to solve the conflicts.
    if conflicts:
        return Editor().edit_files(conflicts)
    else:
        print("No conflicts found in local bzr repository!")
        return 0


def grep_and_edit(top, token, wildcard=None):
    """
    Find all the files located within top that contains token
    and match the shell pattern specified by wildcard.
    Then call $EDITOR to edit the files.
    """
    if wildcard is None: wildcard = WildCard("*.F90|*.finc")

    files = []
    for dirpath, dirnames, filenames in os.walk(top):
        filenames = wildcard.filter(filenames)
        if not filenames: continue

        for f in filenames:
            path = os.path.join(dirpath, f)

            try:
                output = subprocess.check_output(["grep", str(token), path])
            except subprocess.CalledProcessError as exc:
                # Grep returns 1 if token is not found. 2 if real error.
                if exc.returncode == 1:  
                    output = exc.output 
                else:
                    raise

            if output: files.append(path)

    # Call the editor to edit the files.
    print("Will edit %s files" % len(files))
    #print(files)
    return Editor().edit_files(files)


def main():
    def str_examples():
        examples = (
          "\n"
          "Usage example:\n\n" 
          "abzr.py solve                ==> Use $EDITOR to edit the file(s) with conflicts.\n"
          "abzr.py grep 'call wrtout'   ==> Use $EDITOR to edit the file(s) located within .  that containing the token 'call wrtout'.\n"
          "abzr.py grep wrtout 98_main  ==> Use $EDITOR to edit the file(s) located in 98_main containing the token wrtout.\n"
        )
        return examples

    def show_examples_and_exit(err_msg=None, error_code=1):
        """Display the usage of the script."""
        sys.stderr.write(str_examples())
        if err_msg: 
            sys.stderr.write("Fatal Error\n" + err_msg + "\n")
        sys.exit(error_code)

    parser = argparse.ArgumentParser(epilog=str_examples(),formatter_class=argparse.RawDescriptionHelpFormatter)

    # Create the parsers for the sub-commands.
    subparsers = parser.add_subparsers(dest='command', help='sub-command help', description="Valid subcommands")

    # Subparser for single command.
    p_solve = subparsers.add_parser('solve', help="Use $EDITOR to edit the file(s) with conflicts.") 

    p_grep = subparsers.add_parser('grep', help="Use $EDITOR to edit the file(s) containing the given token") 
    p_grep.add_argument("grep_args", nargs="+", help="token [directories].")

    # Parse the command line.
    try:
        options = parser.parse_args()
    except Exception:
        show_examples_and_exit(error_code=1)

    if options.command == "solve":
        return edit_conflicts()

    elif options.command == "grep":
        # Extract token and tops
        if len(options.grep_args) > 1:
            token, tops = options.grep_args[0], options.grep_args[1:]
        else:
            token, tops = options.grep_args[0], ["."]

        for top in tops:
            retcode = grep_and_edit(top, token)
            if retcode: return retcode

if __name__ == "__main__":
    sys.exit(main())
