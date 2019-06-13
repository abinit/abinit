from __future__ import print_function, division, absolute_import #, unicode_literals

import os
import sys
import shutil
import tempfile
import warnings

from subprocess import Popen, PIPE, call


__version__ = "0.1"
__author__ = "Matteo Giantomassi"

__all__ = [
    "RestrictedShell",
    "StringColorizer",
    "Editor",
]

# Helper functions


def patch(fromfile, tofile):
    """
    Use the unix tools diff and patch to patch tofile.
    Returns 0 if success
    """
    # Temporary patch file.
    _, tmp = tempfile.mkstemp(suffix='.patch')

    # An exit status of 0 means no differences were found, 1 means some
    # differences were found, and 2 means trouble.
    #diff_cmd = "diff -c %s %s > %s" % (fromfile, tofile, tmp)
    diff_cmd = "diff -c %s %s > %s" % (tofile, fromfile, tmp)
    print(diff_cmd)

    retcode = os.system(diff_cmd)

    #if retcode >= 2:
    #    warnings.warn("%s returned %s, won't patch!" % (diff_cmd, retcode) )
    #    return retcode

    # Keep a backup copy of tofile.
    bkp = tofile + ".orig"
    shutil.copy(tofile, bkp)

    #print(open(tmp, "r").readlines())
    #patch_cmd = "patch -p1 -i %s -o %s" % (tmp, tofile)
    #patch_cmd = "patch -p0 < %s" % tmp
    patch_cmd = "cp %s %s" % (fromfile, tofile)
    #print(patch_cmd)

    retcode = os.system(patch_cmd)
    if retcode != 0:
        warnings.warn("%s returned %s, reverting to original file" % (patch_cmd, retcode))
        shutil.move(bkp, tofile)
        return retcode

    try:
        os.remove(tmp)
    except IOError:
        pass

    return retcode


def unzip(gz_fname, dest=None):
    """decompress a gz file."""
    import gzip

    if not gz_fname.endswith(".gz"):
        raise ValueError("%s should end with .gz" % gz_fname)

    try:
        gz_fh = gzip.open(gz_fname, 'rb')
        file_content = gz_fh.read()
    finally:
        gz_fh.close() # Cannot use try, except, finally in python2-4

    try:
        if dest is None: dest = gz_fname[:-3]
        out_fh = open(dest, "wb")
        out_fh.write(file_content)
    finally:
        out_fh.close()


def touch(fname, times=None):
    """Emulate unix touch."""
    import os
    with open(fname, 'a'):
        os.utime(fname, times)


def tail_file(fname, n, aslist=False):
    """Emulate unix tail. Assumes a unix-like system."""
    args = ["tail", "-n " + str(n), fname]

    if sys.version_info >= (3, 0):
        p = Popen(args, shell=False, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    else:
        p = Popen(args, shell=False, stdout=PIPE, stderr=PIPE)

    ret_code = p.wait()

    if ret_code != 0:
        raise RuntimeError("return_code = %s, cmd = %s" %(ret_code, " ".join(args) ))

    if aslist:
        return p.stdout.readlines()
    else:
        return p.stdout.read()


def which(program):
    """
    python version of the unix tool which locate a program file in the user's path
    Return:
        None if program cannot be found.
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def tonumber(s):
    """Convert string to number, raise ValueError if s cannot be converted."""
    # Duck test. Much more readable than the ugly strfltrem routine in fldiff.pl
    try:
        stnum = s.upper().replace("D","E")  # D-01 is not recognized by python: Replace it with E.
        stnum = strip_punct(stnum)          # Remove punctuation chars.
        return float(stnum)                 # Try to convert.
    except ValueError:
        raise
    except:
        raise RuntimeError("Don't know how to handle string: " + s)


def nums_and_text(line):
    """split line into (numbers, text)."""
    tokens = line.split()
    text = ""
    numbers = []
    for tok in tokens:
        try:
            numbers.append( tonumber(tok) )
        except ValueError:
            text += " " + tok
    return numbers, text


class RShellError(Exception):
    """Exceptions raised by RestrictedShell"""


class RestrictedShell(object):
    """
    This object executes a restricted set of shell commands.
    It's main goal is to provide a restricted access to the
    computing environment as we want to avoid executing
    arbitrary code passed through the TEST_INFO sections.

    At present, it supports rm, cp, mv and touch
    """
    _key2command = {
        # key (function,   nargs)
        #"rm": (shutil.rmtree, 2),
        "cp": (shutil.copy,   2),
        "mv": (shutil.move,   2),
        "touch":  (touch,     1),
    }

    Error = RShellError

    def __init__(self, inp_dir, workdir, psps_dir):
        """Helper function executing simple commands passed via a string."""

        self.exceptions = []

        self.prefix2dir = {
            "i": os.path.abspath(inp_dir),
            "w": os.path.abspath(workdir),
            "p": os.path.abspath(psps_dir),
        }

    def empty_exceptions(self):
        self.exceptions = []

    def execute(self, string):
        """
        Don't raise exceptions since python threads get stuck.
        Exceptions are stored in self.exceptions
        """
        #print("executing %s" % string)
        _key2command = RestrictedShell._key2command

        tokens = string.split()
        try:
            key, args  = tokens[0], tokens[1:]
            pres, key = key.split("_") # Assume command in the form pre_cmd
            cmd   = _key2command[key][0]
            expected_nargs = _key2command[key][1]
            #print pres, key, cmd, expected_nargs
            nargs = len(args)
            if nargs != expected_nargs:
                err_msg = " Too many arguments, cmd = %s, args = %s " % (cmd, args)
                self.exceptions.append(self.Error(err_msg))
                return
        except:
            err_msg = "Not able to interpret the string: %s " % string
            self.exceptions.append(self.Error(err_msg))
            return

        new_args = []
        for pref, arg in zip(pres, args):
            new_args.append(os.path.join(self.prefix2dir[pref], arg))
        #print "new_args: ",new_args

        try:
            if nargs == 1:
                # Touch
                assert pres == "w"
                return cmd(new_args[0])

            elif nargs == 2:
                # Copy or Move
                src, dest = new_args[0], new_args[1]

                # If src does not exist, look for the netcdf version.
                # and change the extension of dest accordingly.
                # This trick is needed so that we can run the same
                # test both with Fortran-IO as well as with Netcdf
                # without having to use two different TEST_INFO sections..
                if key in ("cp", "mv"):
                    if not os.path.exists(src):
                        print("File %s does not exist, will try netcdf version" % src)
                        if os.path.exists(src + "-etsf.nc"):
                            src += "-etsf.nc"
                            dest += "-etsf.nc"
                        elif os.path.exists(src + ".nc"):
                            src += ".nc"
                            dest += ".nc"
                        else:
                            exc = self.Error("Cannot find neither %s nor netcdf versions in workdir" % src)
                            return self.exceptions.append(exc)
                        print("Fortran file not found, will do:", key, src, dest)

                # Execute command
                return cmd(src, dest)

            else:
                raise NotImplementedError("nargs = %s is too large" % nargs)

        except:
            import sys
            err_msg = "Executing: " + string + "\n" + str(sys.exc_info()[1])
            self.exceptions.append(self.Error(err_msg))


def stream_has_colours(stream):
    """
    True if stream supports colours. Python cookbook, #475186
    """
    if not hasattr(stream, "isatty"):
        return False

    if not stream.isatty():
        return False # auto color only on TTYs
    try:
        import curses
        curses.setupterm()
        return curses.tigetnum("colors") > 2
    except:
        # guess false in case of error
        return False


class StringColorizer(object):
    colours = {
        "default": "",
        "blue":  "\x1b[01;34m",
        "cyan": "\x1b[01;36m",
        "green": "\x1b[01;32m",
        "red": "\x1b[01;31m",
        # lighting colours.
        #"lred":    "\x1b[01;05;37;41m"
        }

    def __init__(self, stream):
        self.has_colours = stream_has_colours(stream)

    def __call__(self, string, colour):
        if self.has_colours:
            code = self.colours.get(colour, "")
            if code:
                return code + string + "\x1b[00m"
            else:
                return string
        else:
            return string


def prompt(question):
    if sys.version_info >= (3, 0):
        my_input = input
    else:
        # python 2.x.
        my_input = raw_input

    return my_input(question)


def user_wants_to_exit():
    """Interactive problem, return False if user entered `n` or `no`."""
    try:
        answer = prompt("Do you want to continue [Y/n]")
    except EOFError:
        return True

    return answer.lower().strip() in ["n", "no"]


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


class PatcherError(Exception):
    """Base class for errors raised by patcher"""


class Patcher(object):
    """Python interface to differt/patcher tools."""
    Error = PatcherError

    # Interactive programs (I known) with support for patches.
    interactive_patchers = [
        "vimdiff",
        "kdiff3",
        "tkdiff",
    ]

    # The revered unix tool for automatic patches.
    auto_patchers = ["patch"]

    known_patchers = interactive_patchers + auto_patchers

    def __init__(self, patcher=None):
        """
        Args:
            patcher:
                string with the name of the utility we want to use
                for generating/applying patchers
                If patcher is None, we use the applications specified
                in the env variables PATHCHER or vimdiff if $PATHCHER is
                not defined
        """
        if patcher is None:
            self.patcher = os.getenv("PATCHER", "vimdiff")

        else:
            if patcher not in Patcher.known_patchers:
                raise ValueError("%s is not supported" % patcher)
            self.patcher = str(patcher)

        if which(self.patcher) is None:
            raise ValueError("Cannot find executable %s in $PATH" % self.patcher)

    @property
    def is_interactive(self):
        """True if this is an interactive patcher."""
        return self.patcher in Patcher.interactive_patchers

    def patch(self, fromfile, tofile):
        """Patch a file."""
        if self.patcher == "patch":
            try:
                return patch(fromfile, tofile)
            except Exception as exc:
                raise self.Error("%s: trying to patch  %s %s:\n%s" % (self.patcher, fromfile, tofile, str(exc)))

        else:
            try:
                return call([self.patcher, fromfile, tofile])
            except Exception as exc:
                raise self.Error("%s: trying to patch  %s %s:\n%s" % (self.patcher, fromfile, tofile, str(exc)))

    def patch_files(self, fromfiles, tofiles):
        """Patch a list of files."""
        assert len(fromfiles) == len(tofiles)
        exit_status = 0

        nfiles = len(fromfiles)
        if not nfiles:
            print("Nothing to patch")
            return exit_status

        if not self.is_interactive:
            ans = prompt("Will patch %d files. Continue [y/n]" % nfiles)
            if ans.lower() not in ["y", "yes"]:
                print("Exit requested by user. Nothing is changed.")
                return 0

        for idx, (fro, to) in enumerate(zip(fromfiles, tofiles)):

            if self.is_interactive:
                #print("from %s, to %s" % (fro, to))
                exit_status = self.patch(fro, to)
                if idx != (nfiles-1) and user_wants_to_exit():
                    break

            else:
                # Automatic patch.
                # 1) Keep a backup copy of to.
                shutil.copy(to, to + ".orig")
                # 2) Apply the patch
                exit_status = self.patch(fro, to)

                if exit_status != 0:
                    msg = "Error while trying to patch %s with %s, interrupting" % (to, fro)
                    break

        return exit_status


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


def ascii_wasp():
    return \
r"""
           | )/ )
        \\ |//,' __
        (")(_)-"()))=-
           (\\

         " ,  ,
            ", ,
               ""     _---.    ..;%%%;, .
                 "" .",  ,  .==% %%%%%%% ' .
                   "", %%%   =%% %%%%%%;  ; ;-_
                   %; %%%%%  .;%;%%%"%p ---; _  '-_
                   %; %%%%% __;%%;p/; O        --_ "-,_
                    q; %%% /v \;%p ;%%%%%;--__    "'-__'-._
                    //\\" // \  % ;%%%%%%%;',/%\_  __  "'-_'\_
                    \  / //   \/   ;%% %; %;/\%%%%;;;;\    "- _\
                       ,"             %;  %%;  %%;;'  ';%       -\-_
                  -=\="             __%    %%;_ |;;    %%%\          \
                                  _/ _=      \==_;;,_ %%%; % -_      /
                                 / /-          =%- ;%%%%; %%;  "--__/
                                //=             ==%-%%;  %; %
                                /             _=_-  d  ;%; ;%;  :F_P:
                                \            =,-"    d%%; ;%%;
                                            //        %  ;%%;
                                           //          d%%%"
                                            \           %%
                                                        V
"""


def ascii_scream():
    return r"""
---;;;;;;;-----'''''''''``'  --- `'  .,,ccc$$hcccccc,.  `' ,;;!!!'``,;;!!'
;;;;,,.,;-------''''''' ,;;!!-    .zJ$$$$$$$$$$$$$$$$$$$c,. `' ,;;!!!!' ,;
  ```'    -;;;!'''''-  `.,..   .zJ$$$$$$$$$$$$$$$$$$$$$$$$$$c, `!!'' ,;!!'
!!-  ' `,;;;;;;;;;;'''''```' ,c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$c,  ;!!'' ,;
,;;;!!!!!!!!''``.,;;;;!'`'  z$$$$$$$$??''''''.,,.`"?$$$$$$$$$$$  ``,;;!!!
;;..       --''```_..,;;!  J$$$$$$??,zcd$$$$$$$$$$$$$$$$$$$$$$$$h  ``'``'
```'''   ,;;''``.,.,;;,  ,$$$$$$F,z$$$$$$$$$$$$$$$$$$$c,`""?$$$$$h
!!!!;;;;,   --`!'''''''  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$h.`"$$$$h .
`'''``.,;;;!;;;--;;   zF,$$$$$$$$$$?????$$$$$$$$$$$$$?????$$r ;?$$$ $.
!;.,..,.````.,;;;;  ,$P'J"$$$$$$P" .,c,,.J$$$$$$$$$"',cc,_`?h.`$$$$ $L
'``````'    .,..  ,$$". $ $$$$P",c$$$$$$$$$$$$$$$$',$$$$$$$$$$ $$$$ $$c,
!!!!!!!!!!!!!'''  J$',$ $.`$$P c$$$$$$$$$$$$$$$$$$,$$$$$$$$$$$ $$$$ $$$$C
   ``            J$ ,$P $$ ?$',$$$$???$$$$$$$$$$$$$$$??'''?$$$ <$$$ $$$$$
c           ;,  z$F,$$  `$$ $ ?$"      "$$$.?$$$ $$$P c??c, ?$.<$$',$$$$$F
$$h.  -!>   ('  $" $F ,F ?$ $ F ,="?$$c,`$$F $$"z$$',$' ,$$P $h.`$ ?$$$$$r
$$$$$hc,. ``'  J$ $P J$ . $$F L ",,J$$$F <$hc$$ "$L,`??????,J$$$.` z$$$$$
$$$$$$$$$$c,'' ?F,$',$F.: $$ c$c,,,,,c,,J$$$$$$$ ?$$$c,,,c$$$$$$F. $$$$$$
`"$$$$$$$$$$$c, $$',$$ :: $$$$$$$$F"',$$$$$$$$$$h ?$$$L;;$$$??$$$$ $$$$$$
   "?$$$$$$$$$$ $$$$$$ : .`F"$$$$$$$$$$$$'''"?'''h $$$$$$$"$,J$$$$ $$$$$'
      "?$$$$$$$ $$$$$$.`.` h `$$$$$$$$$$$cccc$$c,zJ$$$$$P' $$$$$P',$$$$P
$.       `""?$$ $$$$$$$  ` "$c "?$$$$$$$$$$$$??$$$$$$$$" ,J$$$P",J$$$$P
..           `" ?$$$$$$h    ?$$c.`?$$$$$$$$$' . <$$$$$' ,$$$"  ,$$$$$"
!!>. .          `$$$$$$$h  . "$$$c,"$$$$$$$' `' `$$$P  ,$$$' ,c$$$$$'   ;!
```<!!!>         `$$$$$$$c     "$$$c`?$$$$$  : : $$$  ,$$P' z$$$$$$'   ;!!
$hc ```'  ;       `$$$$$$$.      ?$$c ?$$$$ .: : $$$  $$F ,J$$$$$$'   ;!!
.,..      '        `$$$$$$$       "$$h`$$$$ .' ' $$$ ,$$ ,J$$$$$$'    !!!
????P               `$$$$$$L       $$$ $$$F :.: J$$P J$F J$$$$$P     ;!!
-=<                  ?$$."$$       `$$ ?$$' `' z$$$F $P  $$$$$$'     !!'
cc                   `$$$c`?        ?$.`$$hc, cd$$F ,$'  $$$$$$     ;!!
                      $$$$c         `$$c$$$$$$$$$",c$'   $$$$$$     `!!
                      $$$$$          `?$$$$$$$$$$$$P'    $$$$$$> ..
                      $$$$$            `"?$$$$$$$P"      $$$$$$L $$c,
          !!         <$$$$$            zc,`'''',         <$$$$$$.`$$$$cc,
          !!         J$$$$P            `$$$$$$$' !'       $$$$$$L `$$$$$$h
         ;,          $$$$$L          `! J$$$$$',!!        $$$$$$$  `$$$$$$
          '         <$$$$$.           ! $$$$$$ !!         ?$$$$$$   `$$$$$
                   ,$$$$$$$c          `,`???? ;'         c,?$$$$'    `?$$$
                   $$$$$$$??           `!;;;;!     .     `h."?$P      `$$$
                  ,$$$$$$$h.            `'''      `'      `$$$P        `?$
                   $$$$$$$$h                      `!'      `"'           `
                  `$$$$$$$$F          !;     !    ;,
                   `$$$$$$$'         `!!>         `!
c,        ;,        `?$$$$P           !!>             .
"""

def ascii_abinit():
    return r"""
                        -////-                                                        .`
                       `:////:`                                                   :+ymd-
                      `:::///:`                                             `.:oyhdNdo.
                     .::-:////.                                        `.:osddddddho-
   `.:+-`          `-::` -////:`                                `-:/ysddddhhdddo:.
:shhdNmh-         `-/:`  .:///:.                  `.`  ``-.o+yhdhddhhhdddyy+:.
mo++ydhs:-``     `::.     .////-     `.`  .. `:-/:////dhhhhhhhhhmhhso/.-```
hds++++oyyydshhoooos++++++oo+ooo++++s+/ssdmyhhhhyo///yymyhsoy/:..`      `::.
 :+yhhhyyyssooooooosssssssyyyyyyyyyyyo/hyhhsyyhy+oo:-:... `..``    `.`  `/:``
       -:--/ooooosyyyyyyys+//+oo++/:-://:::/:::`  -/. `/::::/:::.  -/:`:://::-
             -/:.             .:://:.:/:.`  `./:. -/. `//-`  `.:/. -/:  `/:.
           `-/:`                `    :/:      ./: -/. `/-      ./- -/:  `/:.
          .:/-`                      -/:.   `./:. -/. `/-      `/- -/:  `/:.
         .::-`                        .:::::::-.  -/. `/-      `/- -/:   -/::-
       `-/:.                            ``..``    `.` `..      `.`.:o/-.``.://:-----::::-````.
       -:-`                                          ``.-.::::/++/+MMMd+oooo++///+ssssooyysysyss++:.
        `                                ` ..-.://+/+/+/::/-.--```-ydh+`                ```.-/+oyyys
                                  ``.-:://+/+::-..```              ` `                       .:/+o/:
                            ``.:://///--.``                                                   `.`
                        `.-////-..`
                      .:::-.`
"""


from functools import wraps


class lazy_property(object):
    """
    lazy_property descriptor

    Used as a decorator to create lazy attributes. Lazy attributes
    are evaluated on first use.
    """

    def __init__(self, func):
        self.__func = func
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



if __name__ == "__main__":
    print(ascii_abinit())
