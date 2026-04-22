"""Tools extracted by AbiPy."""
# TODO: Rationalize modules, merged with pymods

import abc
import os
import sys
import tempfile

from .termcolor import cprint


# Helper functions (coming from AbiPy)
class lazy_property:
    """
    Descriptor for implementing lazy-evaluated properties.

    The decorated method is evaluated only once on first access, after which
    the result is cached on the instance.
    """

    def __init__(self, func):
        self.__func = func
        from functools import wraps
        wraps(self.__func)(self)

    def __get__(self, inst, inst_cls):
        if inst is None:
            return self

        if not hasattr(inst, "__dict__"):
            raise AttributeError("'%s' object has no attribute '__dict__'"
                                 % (inst_cls.__name__,))

        name = self.__name__
        if name.startswith("__") and not name.endswith("__"):
            name = "_%s%s" % (inst_cls.__name__, name)

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

        if not hasattr(inst, "__dict__"):
            raise AttributeError("'%s' object has no attribute '__dict__'"
                                 % (inst_cls.__name__,))

        if name.startswith("__") and not name.endswith("__"):
            name = "_%s%s" % (inst_cls.__name__, name)

        if not isinstance(getattr(inst_cls, name), cls):
            raise AttributeError("'%s.%s' is not a %s attribute"
                                 % (inst_cls.__name__, name, cls.__name__))

        if name in inst.__dict__:
            del inst.__dict__[name]


class Editor:
    """
    Python interface to system text editors.

    Args:
        editor: The name of the editor executable (e.g., "vi", "emacs").
            If None, defaults to the $EDITOR environment variable.
    """
    def __init__(self, editor=None):
        if editor is None:
            self.editor = os.getenv("EDITOR", "vi")
        else:
            self.editor = str(editor)

    def edit_file(self, fname, lineno=None):
        """
        Open a file in the editor, optionally jumping to a specific line.

        Args:
            fname: Path to the file to edit.
            lineno: Optional line number to jump to.

        Returns:
            int: The exit status of the editor process.
        """
        from subprocess import call
        if lineno is None:
            retcode = call([self.editor, fname])
        else:
            # FIXME This works only for vi
            retcode = call([self.editor, fname, "+%s" % str(lineno)])

        if retcode != 0:
            import warnings
            warnings.warn("Error while trying to edit file: %s" % fname)

        return retcode

    def edit_files(self, fnames, ask_for_exit=True):
        """
        Iteratively edit a list of files.

        Args:
            fnames: List of file paths to edit.
            ask_for_exit: If True, prompt the user after each file to see
                if they want to stop the cycle.

        Returns:
            int: The exit status of the last editor invocation.
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
    Print a 2D table of data with aligned columns.

    Args:
        table: A list of lists (rows and columns).
        out: File-like object to write the table to.
        rstrip: If True, remove trailing whitespace from each entry.
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


def print_dataframe(frame, title=None, precision=6, sortby=None, file=sys.stdout, display=None):
    """
    Print a pandas DataFrame in a formatted table.

    Args:
        frame: The pandas DataFrame to print.
        title: Optional title string to print above the table.
        precision: Floating point output precision.
        sortby: Name or list of names to sort the DataFrame by.
        file: File-like object to write to (defaults to sys.stdout).
            If "string", returns the formatted table as a string.
        display: Optional rich display format (e.g., "html") for Jupyter.

    Returns:
        str or None: The formatted string if file="string", else None.
    """
    return_string = file == "string"
    if return_string:
        from io import StringIO
        file = StringIO()

    if title is not None: print(title, file=file)
    if sortby is not None and sortby in frame:
        frame = frame.sort_values(sortby, inplace=False)

    import pandas as pd
    with pd.option_context("display.max_rows", len(frame),
                           "display.max_columns", len(list(frame.keys())),
                           "display.precision", precision,
                           ):
        if display is None:
            print(frame, file=file)
            print(" ", file=file)
            if return_string: return file.getvalue()
        else:
            from IPython.core.display import HTML
            output = getattr(frame, "_repr_%s_" % display)()
            return HTML(output)


# Here metaclass is needed for abc but then I should import six --> leave it at it is without metaclass
#@six.add_metaclass(abc.ABCMeta)
class NotebookWriter: #metaclass=abc.ABCMeta):
    """
    Mixin class for objects capable of generating Jupyter notebooks.

    Subclasses must provide a concrete implementation of `write_notebook`.
    """
    def make_and_open_notebook(self, nbpath=None, foreground=False):  # pragma: no cover
        """
        Generate an jupyter_ notebook and open it in the browser.

        Args:
            nbpath: If nbpath is None, a temporary file is created.
            foreground: By default, jupyter is executed in background and stdout, stderr are redirected
            to devnull. Use foreground to run the process in foreground

        Return:
            system exit code.

        Raise:
            `RuntimeError` if jupyter_ is not in $PATH
        """
        nbpath = self.write_notebook(nbpath=nbpath)

        if which("jupyter") is None:
            raise RuntimeError("Cannot find jupyter in $PATH. Install it with `conda install jupyter or `pip install jupyter`")

        # Use jupyter-lab instead of classic notebook if possible.
        has_jupyterlab = which("jupyter-lab") is not None
        #has_jupyterlab = True
        appname = "jupyter-lab" if has_jupyterlab else "jupyter notebook"

        if foreground:
            return os.system("%s %s" % (appname, nbpath))
        fd, tmpname = tempfile.mkstemp(text=True)
        print(tmpname)
        cmd = "%s %s" % (appname, nbpath)
        print("Executing:", cmd)
        print("stdout and stderr redirected to %s" % tmpname)
        import subprocess
        process = subprocess.Popen(cmd.split(), shell=False, stdout=fd, stderr=fd)
        cprint("pid: %s" % str(process.pid), "yellow")
        return 0

    @staticmethod
    def get_nbformat_nbv():
        """Return nbformat module, notebook version module"""
        import nbformat
        nbv = nbformat.v4
        return nbformat, nbv

    def get_nbformat_nbv_nb(self, title=None):
        """
        Return ``nbformat`` module, notebook version module
        and new notebook with title and import section
        """
        nbformat, nbv = self.get_nbformat_nbv()
        nb = nbv.new_notebook()

        if title is not None:
            nb.cells.append(nbv.new_markdown_cell("## %s" % title))

        nb.cells.extend([
            nbv.new_code_cell("""\
#%matplotlib notebook
""")
        ])

        return nbformat, nbv, nb

    @abc.abstractmethod
    def write_notebook(self, nbpath=None):
        """
        Write a jupyter notebook to nbpath. If nbpath is None, a temporary file is created.
        Return path to the notebook. A typical template:

        .. code-block:: python

            # Preable.
            nbformat, nbv, nb = self.get_nbformat_nbv_nb(title=None)

            #####################
            # Put your code here
            nb.cells.extend([
                nbv.new_markdown_cell("# This is a markdown cell"),
                nbv.new_code_cell("a = 1"),
            ])
            #####################

            # Call _write_nb_nbpath
            return self._write_nb_nbpath(nb, nbpath)
        """

    @staticmethod
    def _write_nb_nbpath(nb, nbpath):
        """
        This method must be called at the end of ``write_notebook``.
        nb is the jupyter notebook and nbpath the argument passed to ``write_notebook``.
        """
        import os
        import tempfile
        if nbpath is None:
            _, nbpath = tempfile.mkstemp(prefix="abinb_", suffix=".ipynb", dir=os.getcwd(), text=True)

        # Write notebook
        import nbformat
        with open(nbpath, "w", encoding="utf8") as fh:
            nbformat.write(nb, fh)
            return nbpath

    @classmethod
    def pickle_load(cls, filepath):
        """
        Loads the object from a pickle file.
        """
        with open(filepath, "rb") as fh:
            new = pickle.load(fh)
            #assert cls is new.__class__
            return new

    def pickle_dump(self, filepath=None):
        """
        Save the status of the object in pickle format.
        If filepath is None, a temporary file is created.

        Return:
            name of the pickle file.
        """
        if filepath is None:
            _, filepath = tempfile.mkstemp(suffix=".pickle")

        with open(filepath, "wb") as fh:
            pickle.dump(self, fh)
            return filepath

    @abc.abstractmethod
    def yield_figs(self, **kwargs):  # pragma: no cover
        """
        This function *generates* a predefined list of matplotlib figures with minimal input from the user.
        Used in abiview.py to get a quick look at the results.
        """

    def expose(self, slide_mode=False, slide_timeout=None, **kwargs):
        """
        Shows a predefined list of matplotlib figures with minimal input from the user.
        """
        from abipy.tools.plotting import MplExpose
        with MplExpose(slide_mode=slide_mode, slide_timeout=slide_mode, verbose=1) as e:
            e(self.yield_figs(**kwargs))


def which(cmd):
    """
    Returns full path to a executable.

    Args:
        cmd (str): Executable command to search for.

    Returns:
        (str) Full path to command. None if it is not found.

    Example::

        full_path_to_python = which("python")
    """

    def is_exe(fp):
        return os.path.isfile(fp) and os.access(fp, os.X_OK)

    fpath, fname = os.path.split(cmd)
    if fpath:
        if is_exe(cmd):
            return cmd
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, cmd)
            if is_exe(exe_file):
                return exe_file
    return None
