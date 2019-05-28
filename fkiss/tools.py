from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os


# Helper functions (coming from AbiPy)
class lazy_property(object):
    """
    lazy_property descriptor

    Used as a decorator to create lazy attributes.
    Lazy attributes are evaluated on first use.
    """

    def __init__(self, func):
        self.__func = func
        from functools import wraps
        wraps(self.__func)(self)

    def __get__(self, inst, inst_cls):
        if inst is None:
            return self

        if not hasattr(inst, '__dict__'):
            raise AttributeError("'%s' object has no attribute '__dict__'"
                                 % (inst_cls.__name__,))

        name = self.__name__
        if name.startswith('__') and not name.endswith('__'):
            name = '_%s%s' % (inst_cls.__name__, name)

        value = self.__func(inst)
        inst.__dict__[name] = value
        return value

    @classmethod
    def invalidate(cls, inst, name):
        """Invalidate a lazy attribute.

        This obviously violates the lazy contract. A subclass of lazy
        may however have a contract where invalidation is appropriate.
        """
        inst_cls = inst.__class__

        if not hasattr(inst, '__dict__'):
            raise AttributeError("'%s' object has no attribute '__dict__'"
                                 % (inst_cls.__name__,))

        if name.startswith('__') and not name.endswith('__'):
            name = '_%s%s' % (inst_cls.__name__, name)

        if not isinstance(getattr(inst_cls, name), cls):
            raise AttributeError("'%s.%s' is not a %s attribute"
                                 % (inst_cls.__name__, name, cls.__name__))

        if name in inst.__dict__:
            del inst.__dict__[name]


class Editor(object):
    """Python interface to text editors."""
    def __init__(self, editor=None):
        if editor is None:
            self.editor = os.getenv("EDITOR", "vi")
        else:
            self.editor = str(editor)

    def edit_file(self, fname, lineno=None):
        from subprocess import call
        if lineno is None:
            retcode = call([self.editor, fname])
        else:
            # FIXME This works only for vi
            retcode = call([self.editor, fname, "+%s" % str(lineno)])

        if retcode != 0:
            warnings.warn("Error while trying to edit file: %s" % fname)

        return retcode

    def edit_files(self, fnames, ask_for_exit=True):
        """
        Edit a list of files, if assk_for_exit is True, we ask
        whether the user wants to exit from the cycle at each iteration.
        """
        exit_status = 0
        for idx, fname in enumerate(fnames):
            exit_status = self.edit_file(fname)
            if ask_for_exit and idx != len(fnames)-1 and user_wants_to_exit():
                break

        return exit_status


def user_wants_to_exit():
    """Interactive problem, return False if user entered `n` or `no`."""
    try:
        answer = prompt("Do you want to continue [Y/n]")
    except EOFError:
        return True

    return answer.lower().strip() in ["n", "no"]


def prompt(question):
    if sys.version_info >= (3, 0):
        my_input = input
    else:
        # python 2.x.
        my_input = raw_input

    return my_input(question)


def pprint_table(table, out=sys.stdout, rstrip=False):
    """
    Prints out a table of data, padded for alignment
    Each row must have the same number of columns.

    Args:
        out:
            Output stream (file-like object)
        table:
            The table to print. A list of lists.
        rstrip:
            if true, trailing withespaces are removed from the entries.
    """
    def max_width_col(table, col_idx):
        """Get the maximum width of the given column index"""
        return max([len(row[col_idx]) for row in table])

    if rstrip:
        for row_idx, row in enumerate(table):
            table[row_idx] = [c.rstrip() for c in row]

    col_paddings = []
    ncols = len(table[0])
    for i in range(ncols):
        col_paddings.append(max_width_col(table, i))

    for row in table:
        # left col
        out.write( row[0].ljust(col_paddings[0] + 1) )
        # rest of the cols
        for i in range(1, len(row)):
            col = row[i].rjust(col_paddings[i] + 2)
            out.write(col)
        out.write("\n")
