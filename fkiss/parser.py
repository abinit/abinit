"""
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import io
import re
import traceback

from textwrap import TextWrapper
from pprint import pformat
from collections import OrderedDict, defaultdict, deque
from .termcolor import cprint
from .tools import lazy_property
from .regex import HasRegex


class Node(object):

    def __repr__(self):
        if self.ancestor is not None:
            return "<%s: %s.%s>" % (self.__class__.__name__, self.ancestor.name, self.name)
        else:
            return "<%s: %s>" % (self.__class__.__name__, self.name)

    def __str__(self):
        # FIXME: This is ABC
        return self.to_string()

    @lazy_property
    def is_procedure(self):
        """Subroutine/function/module/program."""
        return isinstance(self, Procedure)

    @staticmethod
    def terminal_highlight(s, bg="dark"):
        try:
            from pygments import highlight
        except ImportError:
            return s

        from pygments.lexers import FortranLexer
        from pygments.formatters import TerminalFormatter
        return highlight(s, FortranLexer(), TerminalFormatter(bg=bg))


class FortranVariable(Node):
    """
    Fortran Variable.
    type-spec[ [, att] ... :: ] v[/c-list/][, v[/c-list/]] ...
    """
    #def __init__(self,name,vartype,parent,attribs=[],intent="",
    #             optional=False,permission="public",kind=None,
    #             strlen=None,proto=None,doc=[],points=False,initial_value=None):

    def __init__(self, name, ancestor, ftype, shape, kind=None, strlen=None,
                 attribs=None, initial_value=None, doc=None):

        self.name, self.ancestor = name.strip(), ancestor
        self.attribs = () if attribs is None else tuple(filter(None, (a.strip() for a in attribs)))
        self.ftype, self.shape = ftype, shape
        # TODO: dimension
        #print(self.attribs)
        #if "dimension" in attribs:
        #    self.shape = attribs["dimension"].replace(" ", "").replace("dimension", "")
        #    print("shape", self.shape)

        self.initial_value = initial_value.strip() if initial_value else None
        self.doc = "" if doc is None else doc

        if ftype == "character":
            if kind is not None:
                raise ValueError("ftype: %s, kind: %s, strlen: %s" % (ftype, kind, strlen))
            self.strlen = strlen
        else:
            if strlen is not None:
                raise ValueError("ftype: %s, kind: %s, strlen: %s" % (ftype, kind, strlen))
            self.kind = kind

    def to_string(self, verbose=0):
        return "foo"

    @lazy_property
    def is_scalar(self):
        return not bool(self.shape)

    @lazy_property
    def is_array(self):
        return bool(self.shape)

    @lazy_property
    def is_allocatable(self):
        return "allocatable" in self.attribs

    @lazy_property
    def intent(self):
        for i, a in enumerate(self.attribs):
            if a.startswith("intent"):
                return self.attribs[i].replace("intent", "").replace("(", "").replace(")", "").strip()
        return None

    @lazy_property
    def is_pointer(self):
        return "pointer" in self.attribs

    #@lazy_property
    #def is_pointer_set_to_null(self):
    #    if not self.is_pointer: return False
    #    return self.initial_value == "null()"


class Datatype(Node, HasRegex):

    def __init__(self, name, ancestor, preamble, attribs, lines):
        self.name, self.ancestor = name, ancestor
        self.attribs = attribs
        self.preamble = "\n".join(preamble) if preamble else ""
        self.lines = lines
        self._analyzed = False

    def to_string(self, verbose=0):
        s = "\n".join(self.lines)
        return s

    def check_abirules(self, verbose=0):
        retcode = 0
        return retcode

    def analyze(self, verbose=0):
        if self._analyzed: return
        self._analyzed = True
        if verbose: print(self.terminal_highlight("\n".join(self.lines)))
        self.variables = OrderedDict()
        doc = []
        names = None
        visibility = "public"

        # TODO: Handle contains
        for i, line in enumerate(self.lines):
            #line = line.strip()
            if not line or line.startswith("#"): continue
            if i == 0 and not line.lower().startswith("type"): raise ValueError(line)
            if i == len(self.lines) - 1 and not line.lower().startswith("end"): raise ValueError(line)
            if i in (0, len(self.lines) -1): continue

            if line in ("private", "public"):
                visibility = line
                continue

            #fvars = self.parse_variables(line)

            if line.startswith("!"):
                doc.append(line)
            else:
                # New variable found.
                if names is not None:
                    for i, name in enumerate(names):
                        if verbose > 1: print("Creating var:", name)
                        self.variables[name] = FortranVariable(name, self, ftype, shapes[i], kind=kind, strlen=None,
                                                  attribs=attribs, initial_value=initial_values[i], doc="n".join(doc))
                        doc = []

                # Assume: integer, attr1, attr2 :: varname(...)
                line = line.replace(" ", "")
                if verbose: print(line)
                #print(line)
                # Handle inlined comment.
                i = line.find("!")
                if i != -1:
                    doc.append(line[i:])
                    line = line[:i]

                if "::" not in line: raise ValueError(line)
                pre, post = line.split("::")
                toks = pre.split(",")
                ftype = toks[0]
                attribs = [] if len(toks) == 1 else toks[1:]

                # Extract ftype and kind
                # TODO
                kind, strlen = None, None
                m = self.RE_CHARACTER_DEC.match(ftype)
                if m:
                    ftype = "character"
                    strlen, kind = m.group("len"), None

                m = self.RE_TYPECLASS_DEC.match(ftype)
                if m:
                    ftype, kind, strlen = m.group("ftype"), m.group("name"), None

                m = self.RE_NUMBOOL_DEC.match(ftype)
                if m:
                    ftype, kind, strlen = m.group("ftype"), m.group("kind"), None

                # TODO: a(1,2), b, c(3, 4)
                if ")" in post:
                    vlist = post.split("),")
                    for i, v in enumerate(vlist):
                        if "(" in v and not v.endswith(")"):
                            vlist[i] = v + ")"
                else:
                    vlist = post.split(",")

                # Extract default values from vlist (e.g. `a = zero`)
                initial_values = [None] * len(vlist)
                for i, v in enumerate(vlist):
                    if "=" in v:
                        #print(v)
                        v, default = v.split("=", 1)
                        vlist[i] = v
                        initial_values[i] = default

                # Extract shapes from vlist (None if scalar)
                names = [None] * len(vlist)
                shapes = [None] * len(vlist)
                for iv, v in enumerate(vlist):
                    names[iv] = v
                    j = v.find("(")
                    if j != -1:
                        if v[-1] != ")": raise ValueError(v)
                        names[iv] = v[:j]
                        shapes[iv] = v[j:]

                #if verbose:
                #    print(f"ftype={ftype}, attribs={attribs}, vlist={vlist}, names={names}, shapes={shapes}, initial_values={initial_values}\n".format(locals()))

        for i, name in enumerate(names):
            var = FortranVariable(name, self, ftype, shapes[i], kind=kind, strlen=strlen,
                                  attribs=attribs, initial_value=initial_values[i], doc="\n".join(doc))
            self.variables[name] = var


class Interface(Node):

    def __init__(self, name, ancestor, lines):
        self.name, self.ancestor = name, ancestor
        self.lines = lines
        self._analyzed = False

    def to_string(self, verbose=0):
        s = "\n".join(self.lines)
        return s

    def analyze(self, verbose=0):
        if self._analyzed: return
        self._analyzed = True


class Procedure(Node):
    """
    Base class for programs/routines/functions/modules.

    contains: List of contained `Procedure`
    local_uses: List of strings with the name of the modules used (explicit) by this procedure
    includes: List of strings with the name of the included file.
    parents: List of `Procedure` calling this one
    children: List of string with the name of the subroutines called by this procedure.
    """

    def __init__(self, name, ancestor, preamble, line="", prefix=None, arg_names=None, path="<UnknownFile>"):
        self.name = name.strip()
        self.ancestor = ancestor
        self.preamble = "\n".join(preamble) if preamble else ""
        self.line = line
        self.attribs = () if not prefix else tuple(s.strip() for s in prefix.split(","))
        # Init dictionary arg_name --> None (None will replaced by FortranVariable) afterwards.
        self.args = OrderedDict()
        if arg_names is not None:
            for aname in arg_names:
                self.args[aname] = None
        self.path = path
        self.basename = os.path.basename(self.path)
        # This trick is needed for F90.in files that will be post-processed by the build system
        if self.basename.endswith(".F90.in"): self.basename = self.basename[:-3]

        self.num_f90lines, self.num_doclines, self.num_omp_statements = 0, 0, 0
        self.contains, self.local_uses, self.includes = [], [], []
        self.parents, self.children = [], []
        self.types = []
        self.interfaces = []

        # TODO
        # Initialize visibility.
        # The real value will be set by analyzing the module.
        self.visibility = "public"
        #self.has_implicit_none = False

    @lazy_property
    def is_program(self):
        return isinstance(self, Program)

    @lazy_property
    def is_module(self):
        return isinstance(self, Module)

    @lazy_property
    def is_subroutine(self):
        return isinstance(self, Subroutine)

    @lazy_property
    def is_function(self):
        return isinstance(self, Function)

    #@lazy_property
    #def visibility(self):
    #    # visibility of modules is assumed to be initialized.
    #    return self._visibility
    #    if self._visibility is not None: return self._visibility
    #    if self.is_program: return True

    #    if self.ancestor is not None:
    #        if self.ancestor.is_subroutine or self.ancestor.is_function or self.ancestor.is_program:
    #            self._visibility = False
    #            return self._visibility

    #    ancestor = self.ancestor
    #    while ancestor is not None:
    #        if ancestor.is_module:
    #            self._visibility = ancestor.visibility
    #            return self._visibility
    #        ancestor = self.ancestor

    #    # Procedure outside module
    #    #print(repr(self))
    #    #raise RuntimeError("you should not be here!")
    #    self._visibility = True
    #    return self._visibility

    @lazy_property
    def dirpath(self):
        """Absolute path of the directory in which the procedure is located."""
        return None if self.path is None else os.path.dirname(self.path)

    @lazy_property
    def dirname(self):
        """name of the directory in which the procedure is located."""
        return None if self.path is None else os.path.basename(os.path.dirname(self.path))

    @lazy_property
    def dirlevel(self):
        if self.dirname is None:
            return -1
        else:
            # 72_response --> 72
            try:
                return int(self.dirname.split("_")[0])
            except Exception as exc:
                cprint("Cannot extract dirlevel from dirname: `%s`" % self.dirname)
                raise exc

    @property
    def is_public(self):
        return self.visibility == "public"

    @property
    def is_private(self):
        return self.visibility == "private"

    #@lazy_property
    #def global_uses(self)
    #    """String with all the modules used by this procedure (locals + globals)"""

    @property
    def is_contained(self):
        if self.is_module: return False
        return self.ancestor is not None

    @lazy_property
    def public_procedures(self):
        """List of public procedures."""
        if self.is_program:
            return [self]
        elif self.is_module:
            return [self] + [p for p in self.contains if p.is_public]
        elif self.is_subroutine or self.is_function:
            return [self] if self.is_public else []
        raise TypeError("Don't know how to find public entities of type: %s" % type(self))

    def stree(self, level=0):
        lines = [level * "\t" + repr(self)]; app = lines.append
        level += 1
        for p in self.contains:
            app(p.stree(level=level))

        return "\n".join(lines)

    def to_string(self, verbose=0, width=90):
        """
        String representation with verbosity level `verbose`.
        Text is wrapped at `width` columns.
        """
        w = TextWrapper(initial_indent="\t", subsequent_indent="\t", width=width)
        lines = []; app = lines.append

        app("%s: %s\n" % (self.__class__.__name__.upper(), self.name))
        app("Directory: %s" % os.path.basename(self.dirname))

        if self.ancestor is not None:
            app("ANCESTOR:\n\t%s (%s)" % (self.ancestor.name, self.ancestor.__class__.__name__))
        if self.uses:
            app("USES:\n%s\n" % w.fill(", ".join(self.uses)))
            diff = sorted(set(self.local_uses) - set(self.uses))
            if diff:
                app("LOCAL USES:\n%s\n" % w.fill(", ".join(diff)))
        if self.includes:
            app("INCLUDES:\n%s\n" % w.fill(", ".join(self.includes)))
        if self.contains:
            app("CONTAINS:\n%s\n" % w.fill(", ".join(c.name for c in self.contains)))

        if self.types:
            app("DATATYPES:\n%s\n" % w.fill(", ".join(d.name for d in self.types)))
        if self.interfaces:
            app("INTERFACES:\n%s\n" % w.fill(", ".join(i.name for i in self.interfaces)))

        app("PARENTS:\n%s\n" % w.fill(", ".join(sorted(p.name for p in self.parents))))
        #if verbose:
        # Add directory of parents
        dirnames = sorted(set(os.path.basename(p.dirname) for p in self.parents))
        app("PARENT_DIRS:\n%s\n" % w.fill(", ".join(dirnames)))

        app("CHILDREN:\n%s\n" % w.fill(", ".join(sorted(c for c in self.children))))

        if verbose > 1:
            app("")
            app("number of Fortran lines:%s" % self.num_f90lines)
            app("number of doc lines: %s" % self.num_doclines)
            app("number of OpenMP statements: %s" % self.num_omp_statements)
            # Add directory of children
            #dirnames = sorted(set(os.path.basename(p.dirname) for p in self.children))
            #app("CHILDREN_DIRS:\n%s\n" % w.fill(", ".join(dirnames)))
            app("PREAMBLE:\n%s" % self.preamble)

        return "\n".join(lines)

    @lazy_property
    def uses(self):
        """
        List of strings with the modules used by this procedure.
        The list includes the modules used explicitly inside the procedure as well as
        the modules imported at the module level.
        """
        uses = self.local_uses[:]
        if self.is_module or self.is_program:
            for p in self.contains:
                uses.extend(p.local_uses)
                # TODO: should be recursive
                for c in p.contains:
                    uses.extend(c.local_uses)

        elif self.is_subroutine or self.is_function:
            ancestor = self.ancestor
            while ancestor is not None:
                uses.extend(ancestor.local_uses)
                ancestor = ancestor.ancestor

        else:
            raise TypeError("Don't know how to handle %s" % type(self))

        return sorted(set(uses))


class Program(Procedure):
    """Fortran program."""
    proc_type = "program"

    def check_abirules(self, verbose=0):
        retcode = 0
        return retcode


class Function(Procedure):
    """Fortran function."""
    proc_type = "function"

    def check_abirules(self, verbose=0):
        retcode = 0
        return retcode


class Subroutine(Procedure):
    """Fortran subroutine."""
    proc_type = "subroutine"

    def check_abirules(self, verbose=0):
        retcode = 0
        if not self.preamble:
            cprint("Empty preamble in %s" % repr(self), "red")
            retcode += 1

        try:
            header = RobodocHeader.from_string(self.preamble)
        except Exception:
            cprint("Wrong Robodoc header in %s" % repr(self), "red")
            cprint(traceback.format_exc(), "red")
            retcode += 1
        #retcode += header.check_abirules(self)

        return retcode


class Module(Procedure):
    """Fortran module."""
    proc_type = "module"

    def __init__(self, name, ancestor, preamble, path=None):
        super(Module, self).__init__(name, ancestor, preamble, path=path)
        self.default_visibility = True
        #self.variables = OrderedDict()
        #self.public_procedure_names = []
        #self.private_procedure_names = []
        #self.usedby_mods = []

    def to_string(self, verbose=0, width=90):
        lines = []; app = lines.append
        app(super(Module, self).to_string(verbose=verbose, width=width))
        #w = TextWrapper(initial_indent="\t", subsequent_indent="\t", width=width)
        return "\n".join(lines)

    def check_abirules(self, verbose=0):
        retcode = 0

        if not self.preamble:
            cprint("Empty preamble in %s" % repr(self), "red")
            retcode += 1

        for proc in self.contains:
            retcode += proc.check_abirules(verbose=verbose)
        for dtype in self.types:
            retcode += dtype.check_abirules(verbose=verbose)

        return retcode


class FortranKissParser(HasRegex):
    """
    Parse fortran code.
    """

    def __init__(self, macros=None, verbose=0, strict=False):
        self.verbose = verbose
        self.strict = strict
        self.macros = {} if macros is None else macros

    def parse_file(self, path, include_files=False):
        """
        Parse Fortran file in `path`. Include external files if `include_files`.
        """
        with io.open(path, "rt", encoding="utf8") as fh:
            # Include Fortran files?
            if include_files:
                lines = []
                for line in fh:
                    l =  line.strip().replace("'", "").replace('"', "")
                    if l.startswith("#include") and (l.endswith(".finc") or l.endswith(".F90")):
                        basename = l.split()[-1]
                        with io.open(os.path.join(os.path.dirname(path), basename), "rt", encoding="utf8") as incfh:
                            lines.extend(il for il in incfh)
                    else:
                        lines.append(line)
                string = "\n".join(lines)
            else:
                string = fh.read()

            if string and string[0] == "#":
                raise ValueError("Found `#` as first character in file `%s`\n" % path +
                                 "Please avoid it because it has a special meaning and it confuses ctags")

            return self.parse_string(string, path=path)

    def preproc_string(self, string):
        # Preprocess string to facilitate further analysis.
        # Use approach similar to the one used in Ford:
        #
        #     1. remove trailing white-space.
        #     2. ignore blank lines.
        #     3. combine line continuations into single string.
        #       (inlined comments and comment lines inside continuation section are discarded)
        #     4. split lines along semicolons.

        # Replace macros. Needed e.g. to treat USE_DEFS macros in libpaw and tetralib.
        for macro, value in self.macros.items():
            string = re.sub(macro, value, string)

        # 1) Build deque because we are gonna pop a lot!
        lines = deque(l.strip() for l in string.splitlines())

        # List storing pre-processed lines.
        new_lines = []
        napp = new_lines.append

        while lines:
            line = lines.popleft()
            #if not line: continue
            icomm = line.find("!")
            if icomm == 0:
                napp(line)
                continue

            m = self.RE_CONTLINE_START.match(line)
            if not m:
                # 2) "standard" line without continuation.
                ismcol = line.find(";")
                if ismcol == -1:
                    napp(line)
                else:
                    # Be careful when splitting: `integer :: foo ! hello; word`
                    comment = None
                    if icomm != 0:
                        toks = self.quote_split('!', line, strip=False)
                        line = toks[0]
                        if len(toks) > 1: comment = "! " + "".join(toks[1:])
                    new_lines.extend(self.quote_split(';', line, strip=True))
                    if comment:
                        napp(comment)
                continue

            # 3) Handle continuation line
            # http://fortranwiki.org/fortran/show/Continuation+lines
            # NB: I'm not gonna consider `;` in line because it's costly.
            cont_line = m.group("prefix")
            #comment = m.group("postfix")
            while lines:
                line = lines.popleft()
                # Ignore comments inside continuation lines
                if line.startswith("!"): continue

                # `& foo` triggers immediate exit.
                if line.startswith("&") and line.count("&") == 1:
                    sv = self.trim_comment(line.replace("&", ""))
                    cont_line += sv
                    if sv:
                        napp(cont_line)
                        break

                m = self.RE_CONTLINE_NEXT.match(line)
                if m:
                    # `& call &"` OR `& end` OR ` continue &`
                    value = m.group("value")
                    cont_line += value
                    sv = self.trim_comment(value)
                    #print("value", value, "\nsv:", sv)
                    # This is not clear now but it works
                    if sv and sv[0] != "&" and sv[-1] == "&":
                        napp(cont_line)
                        break
                else:
                    # First line without & e.g. `end`
                    # Add it to cont_line and exit.
                    cont_line += line
                    napp(cont_line)
                    break

        if self.verbose >= 3:
            print("WILL OPERATE ON LINES\n", "\n".join(new_lines))

        return deque(new_lines)

    def parse_string(self, string, path=None):
        if path is None: path = "<UnknownPath>"
        self.path = path

        self.lines = self.preproc_string(string)
        self.warnings = []

        self.num_doclines, self.num_f90lines, self.num_omp_statements = 0, 0, 0
        self.preamble, self.stack = [], []
        self.all_includes, self.all_uses  = [], []
        self.ancestor = None

        # Invokations of Fortran functions are difficult to handle
        # without inspecting local variables so we only handle explicit calls to routines.
        # in principle I may re-read the source and use regex for val = foo() where foo is
        # one of the functions in the project but it's gonna be costly.
        while self.lines:
            line = self.lines.popleft()
            if not line: continue
            if self.handle_comment(line): continue
            # Convert to lower case here so that we don't have to deal with case.
            line = line.lower()
            if self.handle_use_statement(line): continue
            if self.handle_cpp_line(line): continue
            if self.handle_contains(line): continue
            if self.handle_call(line): continue
            # subroutine or function declaration.
            if self.handle_procedure(line): continue
            if self.consume_module_header(line): continue
            # Handle `end module`
            m = self.RE_MOD_END.match(line)
            if m:
                self.close_stack_entry(line, end_proc_type="module", end_name=m.group("name"))
                continue

            #print("Ignored line:", line)
            self.num_f90lines += 1

        self.all_includes = sorted(set(self.all_includes))
        self.all_uses = sorted(set(self.all_uses))

        # Extract data from stack.
        # TODO: Support visibility But I need to parse the first portion of the header.
        # to handle e.g. public :: foo
        self.programs, self.modules, self.subroutines, self.functions = [], [], [], []

        while self.stack:
            p, status = self.stack.pop(0)
            if status != "closed":
                self.warn("Unclosed %s with status %s, ancestor: %s" % (repr(p), status, repr(p.ancestor)))

            # Sort entries here.
            p.local_uses = sorted(p.local_uses)
            p.children = sorted(set(p.children))

            if p.ancestor is not None:
                #print("Adding %s to ancestor %s" % (repr(p), repr(p.ancestor)))
                p.ancestor.contains.append(p)
            else:
                if p.is_module: self.modules.append(p)
                elif p.is_program: self.programs.append(p)
                # Here only if p is subroutine or function outside module.
                elif p.is_subroutine: self.subroutines.append(p)
                elif p.is_function: self.functions.append(p)
                else: raise ValueError("Don't know how to handle type `%s`" % type(p))

        return self

    @staticmethod
    def trim_comment(line):
        i = line.find("!")
        if i != -1: line = line[:i]
        return line.strip()

    def warn(self, msg):
        cprint(msg, color="yellow")
        if not self.strict:
            self.warnings.append(msg)
        else:
            raise RuntimeError(msg)

    def handle_contains(self, line):
        m = self.RE_CONTAINS.match(line)
        if not m: return False
        self.ancestor = self.stack[-1][0]
        self.preamble = []
        if self.verbose > 1: print("Setting ancestor to:", repr(self.ancestor))
        return True

    def handle_cpp_line(self, line):
        # Handle include statement (CPP or Fortran version).
        m = self.RE_INCLUDE.match(line)
        if m:
            path = re.sub(r"'|\"", "", m.group("path"))
            if self.stack:
                self.stack[-1][0].includes.append(path)
            else:
                self.all_includes.append(path)
            return True

        return True if line[0] == "#" else False

    def handle_comment(self, line):
        # Count number of comments and code line
        # Inlined comments are not counted (also because I don't like them)
        m = self.RE_F90COMMENT.match(line)
        if not m:
            #self.num_f90lines += 1
            return False

        # Count (non-emtpy) comments and OMP (preamble is included in num_doclines)
        omp = self.RE_OMP_SENTINEL.match(line)
        if omp:
            self.num_omp_statements += 1

        if not omp and m.group("value").strip():
            self.num_doclines += 1

        # Robodoc (deactivated for the time being)
        if False:
            m = RobodocHeader.RE_HEADER_START.match(line)
            if m:
                robo_lines = [line]
                while self.lines:
                    line = self.lines.popleft()
                    if not line.startswith("!!"):
                        self.lines.appendleft(line)
                        try:
                            header = RobodocHeader.from_lines(robo_lines)
                            #break
                        except Exception as exc:
                            cprint("Wrong Robodoc header in %s" % self.path, "red")
                            cprint(str(exc), "red")
                        finally:
                            break
                    else:
                        robo_lines.append(line)
                else:
                    raise ValueError("You should not be here.")

        # TODO Must handle comments after contains here!
        #if not self.stack or (self.stack and self.stack[-1][1] != "open"):
        if self.preamble is not None:
            # Ignore stupid robodoc marker.
            if line != "!!***":
                self.preamble.append(line)

        return True

    def handle_use_statement(self, line):
        # Find use statements and the corresponding module
        if not line.startswith("use "): return False
        smod = line.split()[1].split(",")[0].lower()
        # Remove comment at the end of the line if present.
        i = smod.find("!")
        if i != -1: smod = smod[:i]
        if self.verbose > 1: print("Found used module:", smod)
        self.stack[-1][0].local_uses.append(smod)
        self.all_uses.append(smod)
        return True

    def handle_call(self, line):
        # At this level subname is a string that will be replaced by a Procedure object afterwards
        # TODO: should handle `call obj%foo()` syntax.
        m = self.RE_SUBCALL.match(line)
        if not m: return False
        subname = m.group("name")
        assert subname
        if self.verbose: print("Adding %s to children of %s" % (subname, repr(self.stack[-1][0])))
        self.stack[-1][0].children.append(subname)
        return True

    def consume_module_header(self, line):
        m = self.RE_MOD_START.match(line)
        if not m: return False
        name = m.group("name")
        if self.verbose > 1: print("Entering module:", name)
        module = Module(name, self.ancestor, self.preamble, path=self.path)
        self.add_node_to_stack(module)

        while self.lines:
            line = self.lines.popleft()
            if not line: continue
            if self.handle_comment(line): continue
            # Convert to lower case here so that we don't have to deal with case.
            line = line.lower()
            if self.handle_use_statement(line): continue

            m = self.RE_PUB_OR_PRIVATE.match(line)
            if m:
                module.default_visibility = m.group("name")
                continue
            # TODO: Procedure declaration statement.
            #proc_names
            #module.public_procedure_names.extend(proc_names)
            #module.private_procedure_names.extend(proc_names)

            # Interface declaration.
            if self.consume_interface(line): continue

            # Datatype declaration.
            if self.consume_datatype(line): continue

            # Exit here
            if self.handle_contains(line): return True

            # or here if the module does not have *contains*
            m = self.RE_MOD_END.match(line)
            if m:
                self.close_stack_entry(line, end_proc_type="module", end_name=m.group("name"))
                return True

        else:
            raise ValueError("Cannot find `contains` in %s" % self.path)

    def consume_interface(self, line):
        m = self.RE_INTERFACE_START.match(line)
        if not m: return False
        #assert self.stack[-1][1] == "open"
        buflines = [line]
        name = m.group("name")
        if self.verbose > 1: print("begin interface", name, "in line:", line)
        while self.lines:
            line = self.lines.popleft()
            buflines.append(line)
            end_match = self.RE_INTERFACE_END.match(line)
            if end_match:
                if self.verbose > 1: print("end interface", line)
                # Add interface to the last item on the stack.
                # NB Don't enforce name `end interface [name]`
                self.stack[-1][0].interfaces.append(Interface(name, self.ancestor, buflines))
                return True
        else:
            raise ValueError("Cannot find `end interface %s` in %s" % (name, self.path))

    def consume_datatype(self, line):
        m = self.RE_TYPE_START.match(line)
        if not m: return False
        buflines = [line]
        name = m.group("name")
        if self.verbose > 1: print("begin datatype", name, "in line:", line)
        # Extract attributes and put them in a tuple.
        attribs = m.group("attribs").lower().replace(":", "")
        attribs = tuple(filter(None, (s.strip() for s in attribs.split(",")))) if attribs else ()

        while self.lines:
            line = self.lines.popleft()
            buflines.append(line)
            if self.handle_comment(line): continue
            line = line.lower()
            end_match = self.RE_TYPE_END.match(line)
            if end_match:
                end_name = end_match.group("name")
                if self.verbose > 1: print("end datatype %s in %s" % (end_name, line))
                if name != end_name:
                    self.warn("Cannot find `end type %s` in %s" % (name, self.path))
                    #raise ValueError("Cannot find `end type %s` in %s" % (name, self.path))
                # Add type to the last item on the stack.
                dtype = Datatype(name, self.ancestor, self.preamble, attribs, buflines)
                #try:
                #    dtype.analyze()
                #except:
                #    pass
                self.stack[-1][0].types.append(dtype)
                return True
        else:
            raise ValueError("Cannot find `end type %s` in %s" % (name, self.path))

    def handle_procedure(self, line):
        if not self.RE_SEARCH_PROC.search(line): return False
        # Find if (subroutine|function|program) and select regex for end tag.
        proc_type = "subroutine"
        m = self.RE_SUB_START.match(line)
        re_end = self.RE_SUB_END
        if not m:
            m = self.RE_FUNC_START.match(line)
            if m:
                proc_type, re_end = "function", self.RE_FUNC_END
        if not m:
            m = self.RE_PROG_START.match(line)
            if m:
                proc_type, re_end = "program", self.RE_PROG_END
        if not m:
            return False

        # Extract procedure name.
        name = m.group("name")
        if not name:
            raise ValueError("Cannot find %s name in line `%s`" % (proc_type, line))

        if self.verbose > 1:
            print("Found `%s %s`" % (proc_type, name), "at line:\n\t", line)
            print("Ancestor is set to", repr(self.ancestor))

        # Extract procedure arguments from prototype string.
        arg_names = []
        if proc_type != "program":
            margs = self.RE_SUB_ARGS.search(line)
            # This for `subroutine start_exx` Is it so difficult to add ()?
            empty_args = not margs and line.strip().endswith(name)
            if not margs and not empty_args:
                self.warn("Cannot extract arguments from line:\n%s\n" % line)
            else:
                arg_names = margs.group("args") if not empty_args else []
                if arg_names: arg_names = [a.strip() for a in arg_names.split(",")]

        if proc_type == "subroutine":
            prefix = m.group("prefix")
            new_node = Subroutine(name, self.ancestor, self.preamble,
                line=line, prefix=prefix, arg_names=arg_names, path=self.path)
        elif proc_type == "function":
            prefix = m.group("prefix")
            new_node = Function(name, self.ancestor, self.preamble,
                line=line, prefix=prefix, arg_names=arg_names, path=self.path)
        elif proc_type == "program":
            new_node = Program(name, self.ancestor, self.preamble, arg_names=arg_names, path=self.path)
        else:
            raise ValueError(proc_type)

        self.add_node_to_stack(new_node)

        has_contains = False
        while self.lines:
            line = self.lines.popleft()
            if not line: continue
            if self.handle_comment(line): continue
            # Convert to lower case here so that we don't have to deal with case.
            line = line.lower()
            if self.handle_call(line): continue
            if self.handle_args(line): continue
            if self.handle_use_statement(line): continue
            if self.handle_cpp_line(line): continue
            if self.consume_interface(line): continue
            cont = self.handle_contains(line)
            if cont: has_contains = True
            if has_contains:
                # Search for contained procedures (recursively).
                if self.handle_procedure(line): continue

            # Exit if `end proc_type`
            m = self.RE_PROC_END.match(line)
            if m:
                end_name = m.group("name")
                end_proc_type =m.group("proc_type")
                # Warn if `end [proc_type [name]]`
                if end_proc_type != proc_type:
                    self.warn("Cannot find `end %s %s` in %s" % (proc_type, name, self.path))
                if end_name != name:
                    self.warn("Cannot find `end %s %s` in %s" % (proc_type, name, self.path))

                self.close_stack_entry(line, end_proc_type, end_name)
                break

            self.num_f90lines += 1
        else:
            raise ValueError("Cannot find `end %s `%s`" % (proc_type, name))

        return True

    def handle_args(self, line):
        m = self.RE_INTENT.search(line)
        if not m: return False
        return True

        fvars = self.parse_variables(line)
        #if fvars:
        #   docs = []
        #   while True:
        #       line = self.lines.popleft()
        #       if not line.startswith("!"):
        #          for v.name in fvars:
        #               v.doc = "\n".join(docs)
        #          self.lines.appendleft(line)
        #          break
        #        docs.append(line)
        #   continue

        # Set procedure arguments.
        proc = self.stack[-1][0]
        for var in fvars:
            if var.name not in proc.args:
                #raise RuntimeError("Cannot find `%s` in `%s`\nin file: %s" % (var.name, list(proc.args.keys()), self.path))
                self.warn("Cannot find `%s` in `%s`\nin file: %s" % (var.name, list(proc.args.keys()), self.path))
            else:
                proc.args[var.name] = var

        return True

    def add_node_to_stack(self, node):
        self.ancestor = node
        self.stack.append([node, "open"])
        self.preamble = None
        # TODO: Recheck this part (open, end, accumulate?)
        self.num_f90lines, self.num_doclines, self.num_omp_statements = 0, 0, 0

    def close_stack_entry(self, line, end_proc_type, end_name):
        if end_name:
            # Close the last entry in the stack with name == end_name.
            for item in reversed(self.stack):
                if item[0].name == end_name:
                    node = item[0]
                    item[1] = "closed"
                    break
            else:
                raise RuntimeError("Cannot find end_name `%s` in stack:\n%s\npath: %s\nLast line:%s" % (
                    end_name, pformat([s[0].name for s in self.stack]), self.path, line))
        else:
            # Close the last entry in the stack with end_proc_type.
            if end_proc_type is not None:
                self.warn("Found `end %s` without name in %s:%s" % (end_proc_type, self.path, line))
                for item in reversed(self.stack):
                    if item[0].proc_type == end_proc_type:
                        node = item[0]
                        item[1] = "closed"
                        break
                else:
                    raise RuntimeError("Cannot find end_proc_type `%s` in stack:\n%s\nLast line:%s" % (
                        end_proc_type, pformat([s[0].proc_type for s in self.stack]), line))
            else:
                # This is the best I can do without any info.
                self.warn("Found plain `end` without procedure_type and name in %s:%s" % (self.path, line))
                self.stack[-1][1] = "closed"
                node = stack[-1][0]

        if self.verbose > 1: print("Closing", repr(node))
        if self.ancestor is not None and self.ancestor.name == end_name:
            self.ancestor = self.ancestor.ancestor

        # Set attributes of last node.
        node.num_f90lines = self.num_f90lines
        node.num_doclines = self.num_doclines
        node.num_omp_statements = self.num_omp_statements

        self.num_f90lines, self.num_doclines, self.num_omp_statements = 0, 0, 0

        if self.verbose and self.preamble:
            self.warn("%s, preamble:\n%s" % (self.path, "\n".join(self.preamble)))

        self.preamble = []

    def parse_variables(self, line):
        # Remove inlined comment from line (if any)
        icomm, doc = line.find("!"), ""
        if icomm != -1:
            line, doc = line[:icomm], line[icomm:]

        # ftype [attribs] :: var_list
        if "::" not in line:
            # TODO
            raise ValueError(line)
            toks = line.split()
            ftype, post = toks[0], toks[1:]
        else:
            pre, post = line.split("::")
            toks = pre.split(",")
            ftype = toks[0]
            attribs = [] if len(toks) == 1 else toks[1:]

        # Extract ftype and kind
        # TODO
        kind, strlen = None, None
        m = self.RE_CHARACTER_DEC.match(ftype)
        if m:
            ftype = "character"
            strlen, kind = m.group("len"), None
        if not m:
            m = self.RE_TYPECLASS_DEC.match(ftype)
            if m:
                ftype, kind, strlen = m.group("ftype"), m.group("name"), None
        if not m:
            m = self.RE_NUMBOOL_DEC.match(ftype)
            if m:
                ftype, kind, strlen = m.group("ftype"), m.group("kind"), None
        if not m:
            #raise ValueError("Cannot find Fortran type in line: %s. file: %s" % (line, self.path))
            self.warn("Cannot find Fortran type in line: %s. file: %s" % (line, self.path))
            return []

        # TODO: a(1,2), b, c(3, 4)
        if ")" in post:
            vlist = post.split("),")
            for i, v in enumerate(vlist):
                if "(" in v and not v.endswith(")"):
                    vlist[i] = v + ")"
        else:
            vlist = post.split(",")

        # Extract default values from vlist (e.g. `a = zero`)
        initial_values = [None] * len(vlist)
        for i, v in enumerate(vlist):
            if "=" in v:
                v, default = v.split("=", 1)
                vlist[i] = v
                initial_values[i] = default.lower().strip()

        # Strip and convert to lower case.
        vlist = [v.strip().lower() for v in vlist]
        if kind is not None: kind = kind.strip().lower()

        # Extract shapes from vlist (None if scalar)
        names = [None] * len(vlist)
        shapes = [None] * len(vlist)
        for iv, v in enumerate(vlist):
            names[iv] = v
            j = v.find("(")
            if j != -1:
                if v[-1] != ")": raise ValueError(v)
                names[iv] = v[:j]
                shapes[iv] = v[j:]

        #if verbose:
        #    print(f"ftype={ftype}, attribs={attribs}, vlist={vlist}, names={names}, shapes={shapes}, initial_values={initial_values}\n".format(locals()))

        fvars = []
        for i, name in enumerate(names):
            var = FortranVariable(name, self, ftype, shapes[i], kind=kind, strlen=strlen,
                                  attribs=attribs, initial_value=initial_values[i], doc=doc)
            fvars.append(var)

        return fvars

    @staticmethod
    def quote_split(sep, string, strip=False):
        """
        Splits the strings into pieces divided by sep, when sep in not inside quotes.
        Copied from https://github.com/Fortran-FOSS-Programmers/ford/blob/master/ford/utils.py
        """
        if len(sep) != 1: raise Exception("Separation string must be one character long")
        retlist = []
        squote = False
        dquote = False
        left = 0
        i = 0
        while i < len(string):
            if string[i] == '"' and not dquote:
                if not squote:
                    squote = True
                elif (i+1) < len(string) and string[i+1] == '"':
                    i += 1
                else:
                    squote = False
            elif string[i] == "'" and not squote:
                if not dquote:
                    dquote = True
                elif (i+1) < len(string) and string[i+1] == "'":
                    i += 1
                else:
                    dquote = False
            elif string[i] == sep and not dquote and not squote:
                retlist.append(string[left:i])
                left = i + 1
            i += 1

        retlist.append(string[left:])

        return [s.strip() for s in retlist] if strip else retlist


class RobodocHeader(OrderedDict):
    # See config/robodoc/robodoc-html.rc

    #ALL_KEYS = [
    #    "NAME", "COPYRIGHT", "FUNCTION",
    #    "INPUTS", "OUTPUT", "OUTPUTS", "SIDE EFFECTS",
    #    "NOTES", "TODO", "PARENTS", "CHILDREN", "SOURCE",
    #]

    ALL_KEYS = [
    "SOURCE",
    "NAME",
    "COPYRIGHT",
    "SYNOPSIS",
    "USAGE",
    "FUNCTION",
    "DESCRIPTION",
    "PURPOSE",
    "AUTHOR",
    "CREATION DATE",
    "MODIFICATION HISTORY",
    "HISTORY",
    "INPUTS",
    "ARGUMENTS",
    "OPTIONS",
    "PARAMETERS",
    "SWITCHES",
    "OUTPUT",
    "OUTPUTS",    # This is not in the official list
    "SIDE EFFECTS",
    "RESULT",
    "RETURN VALUE",
    "EXAMPLE",
    "NOTES",
    "DIAGNOSTICS",
    "WARNINGS",
    "ERRORS",
    "BUGS",
    "TODO",
    "IDEAS",
    "PORTABILITY",
    "SEE ALSO",
    "METHODS",
    "NEW METHODS",
    "ATTRIBUTES",
    "NEW ATTRIBUTES",
    "TAGS",
    "COMMANDS",
    "DERIVED FROM",
    "DERIVED BY",
    "USES",
    "CHILDREN",
    "USED BY",
    "PARENTS",
    ]

    header_markers = "!!****"
    remark_markers = "!!"
    end_markers = "!!***"

    #headertypes:
    #    p "Programs" robo_programs
    #    m "Modules"  robo_modules
    #    d "Directories" robo_directories

    # Detect robodoc header (!****)
    RE_HEADER_START = re.compile(r"^!!\*{4}(?P<ptype>[a-z])\*\s+(?P<name>.+?)$", re.MULTILINE)

    @classmethod
    def from_string(cls, s):
        return cls.from_lines(s.splitlines())

    @classmethod
    def from_lines(cls, lines):
        #m = cls.RE_HEADER_START.search(s)
        m = cls.RE_HEADER_START.match(lines[0])
        if not m:
            raise ValueError("Cannot find robodoc header in string:\n%s" % lines[0])

        new = cls()
        # Save ptype and name
        new.ptype, new.name = m.group("ptype"), m.group("name")

        #s = s[m.end(2) + 1:]
        #lines = s.splitlines()
        active_key = None
        for line in lines[1:]:
            #print(line)
            if not line.startswith("!!"): break
            s = line[2:]
            if active_key is None and not s: continue
            k = s.strip()
            if k in cls.ALL_KEYS:
                if new[k]:
                    raise ValueError("Key %s already present in robodoc header.\n%s" % (k, "\n".join(lines)))
                new[k] = []
                active_key = k
                continue

            if active_key is None:
                raise ValueError("active_key cannot be None.\n%s" % "\n".join(lines))

            new[active_key].append(line)

        for key, lines in new.items():
            new[key] = "\n".join(lines)

        return new

    def __init__(self, *args, **kwargs):
        super(RobodocHeader, self).__init__(*args, **kwargs)
        for k in self.ALL_KEYS:
            self[k] = []

    def __str__(self):
        return self.to_string(self)

    def to_string(self, verbose=0):
        lines = []
        app = lines.append
        for key, value in self.items():
            if not value: continue
            app("!! %s" % key)
            app(self[key])
            app("!!")

        return "\n".join(lines)
