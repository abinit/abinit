"""
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import re

from textwrap import TextWrapper
from pprint import pformat
from collections import OrderedDict, defaultdict
from .tools import lazy_property


class Procedure(object):
    """
    Base class

    contains: List of contained `Procedure`
    local_uses: List of strings with the name of the modules used (explicit) by this procedure
    includes: List of strings with the name of the included file.
    parents: List of `Procedure` calling this one
    children: List of `Procedure` called by this one
    """

    def __init__(self, name, ancestor, preamble, path=None):
        self.name = name.strip()
        self.ancestor = ancestor
        self.preamble = preamble
        self.path = path

        self.num_f90lines, self.num_doclines = 0, 0
        self.contains, self.local_uses, self.includes = [], [], []
        self.parents, self.children = [], []
        #self.visibility = "public"

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

    @lazy_property
    def dirname(self):
        return None if self.path is None else os.path.dirname(self.path)

    def __repr__(self):
        #anc_name = self.ancestor.name if self.ancestor is not None else "None"
        #return "<%s: %s, ancestor %s>" % (self.__class__.__name__, self.name, anc_name)
        return "<%s: %s>" % (self.__class__.__name__, self.name)

    @lazy_property
    def public_procedures(self):
        if self.is_program:
            return [self]
        elif self.is_module:
            return [self] + [p for p in self.contains] # if p.is_public
        elif self.is_subroutine or self.is_function:
            return [self] # if self.is_public else []
        raise TypeError("Don't know how to find public entities of type: %s" % type(self))

    def stree(self, level=0):
        lines = [level * "\t" + repr(self)]; app = lines.append
        level += 1
        #if self.is_module:
        for p in self.contains:
            app(p.stree(level=level))

        return "\n".join(lines)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0, width=90):
        w = TextWrapper(initial_indent="\t", subsequent_indent="\t", width=width)
        lines = []; app = lines.append

        app("%s: %s" % (self.__class__.__name__, self.name)),
        if self.ancestor is not None:
            app("Ancestor: %s (%s)" % (self.ancestor.name, self.ancestor.ftype))
        app("")
        if self.uses:
            app("USES:\n%s\n" % w.fill(", ".join(self.uses)))
            diff = sorted(set(self.local_uses) - set(self.uses))
            if diff:
                app("LOCAL USES:\n%s\n" % w.fill(", ".join(diff)))
        if self.includes:
            app("INCLUDES:\n%s\n" % w.fill(", ".join(self.includes)))
        if self.contains:
            app("CONTAINS:\n%s\n" % w.fill(", ".join(c.name for c in self.contains)))

        app("PARENTS:\n%s\n" % w.fill(", ".join(p.name for p in self.parents)))
        app("CHILDREN:\n%s\n" % w.fill(", ".join(c for c in self.children)))

        return "\n".join(lines)

    #@property
    #def is_public(self):
    #    return self.visibility == "public"

    #@property
    #def is_private(self):
    #    return self.visibility == "private"

    #@lazy_property
    #def global_uses(self)
    #    """String with all the modules used by thie procedure (locals + globals)"""

    @lazy_property
    def uses(self):
        """
        List of strings with all the modules used by this procedure (locals + globals)
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

    #def analyze(self):
    #    if self.is_module:
    #        for p in self.contains:
    #            p.visibility = self.visibility
    #    else:
    #        for p in self.contains:
    #            p.visibility = "private"


class Program(Procedure):
    """Fortran program."""
    ftype = "program"


class Function(Procedure):
    """Fortran function."""
    ftype = "function"


class Subroutine(Procedure):
    """Fortran subroutine."""
    ftype = "subroutine"


class Module(Procedure):
    """Fortran module."""
    ftype = "module"

    def to_string(self, verbose=0, width=90):
        lines = []; app = lines.append
        app(super(Module, self).to_string(verbose=verbose, width=width))
        #w = TextWrapper(initial_indent="\t", subsequent_indent="\t", width=width)

        return "\n".join(lines)

# This one is problematic. The firs one is stricter but it does not play well with the other code!!
# I use the old one used by abilint
#RE_FUNC_END = re.compile('^[ \t]*end[ \t]*(function\n)',re.I).match

#Detect call (not after an if or ;)
#re_call = re.compile('(^[ \t]*call[ ]*)(\w+)', re.MULTILINE+re.I)

# Taken from https://github.com/cmacmackin/ford/blob/master/ford/sourceform.py
#VAR_TYPE_STRING = "^integer|real|double\s*precision|character|complex|logical|type|class|procedure|enumerator"
#VARKIND_RE = re.compile("\((.*)\)|\*\s*(\d+|\(.*\))")
#KIND_RE = re.compile("kind\s*=\s*",re.I)
#LEN_RE = re.compile("len\s*=\s*",re.I)
#ATTRIBSPLIT_RE = re.compile(",\s*(\w.*?)::\s*(.*)\s*")
#ATTRIBSPLIT2_RE = re.compile("\s*(::)?\s*(.*)\s*")
#ASSIGN_RE = re.compile("(\w+\s*(?:\([^=]*\)))\s*=(?!>)(?:\s*([^\s]+))?")
#POINT_RE = re.compile("(\w+\s*(?:\([^=>]*\)))\s*=>(?:\s*([^\s]+))?")
#EXTENDS_RE = re.compile("extends\s*\(\s*([^()\s]+)\s*\)")
#DOUBLE_PREC_RE = re.compile("double\s+precision",re.I)
#QUOTES_RE = re.compile("\"([^\"]|\"\")*\"|'([^']|'')*'",re.I)
#PARA_CAPTURE_RE = re.compile("<p>.*?</p>",re.I|re.DOTALL)
#COMMA_RE = re.compile(",(?!\s)")
#NBSP_RE = re.compile(" (?= )|(?<= ) ")
#DIM_RE = re.compile("^\w+\s*(\(.*\))\s*$")

#ATTRIB_RE = re.compile("^(asynchronous|allocatable|bind\s*\(.*\)|data|dimension|external|intent\s*\(\s*\w+\s*\)|optional|parameter|pointer|private|protected|public|save|target|value|volatile)(?:\s+|\s*::\s*)((/|\(|\w).*?)\s*$",re.I)
#END_RE = re.compile("^end\s*(?:(module|submodule|subroutine|function|procedure|program|type|interface|enum|block\sdata|block|associate)(?:\s+(\w.*))?)?$",re.I)
#BLOCK_RE = re.compile("^(\w+\s*:)?\s*block\s*$",re.I)
#BLOCK_DATA_RE = re.compile('^block\s*data\s*(\w+)?\s*$',re.I)
#ASSOCIATE_RE = re.compile("^(\w+\s*:)?\s*associate\s*\((.+)\)\s*$",re.I)
#ENUM_RE = re.compile("^enum\s*,\s*bind\s*\(.*\)\s*$",re.I)
#MODPROC_RE = re.compile("^(module\s+)?procedure\s*(?:::|\s)\s*(\w.*)$",re.I)
#MODULE_RE = re.compile("^module(?:\s+(\w+))?$",re.I)
#SUBMODULE_RE = re.compile("^submodule\s*\(\s*(\w+)\s*(?::\s*(\w+))?\s*\)\s*(?:::|\s)\s*(\w+)$",re.I)
#PROGRAM_RE = re.compile("^program(?:\s+(\w+))?$",re.I)
#SUBROUTINE_RE = re.compile("^\s*(?:(.+?)\s+)?subroutine\s+(\w+)\s*(\([^()]*\))?(?:\s*bind\s*\(\s*(.*)\s*\))?$",re.I)
#TYPE_RE = re.compile("^type(?:\s+|\s*(,.*)?::\s*)((?!(?:is\s*\())\w+)\s*(\([^()]*\))?\s*$",re.I)
#RE_MOD_END = re.compile(r'end(\s*module\s*\w*|)\Z', re.I)
##~ ABS_INTERFACE_RE = re.compile("^abstract\s+interface(?:\s+(\S.+))?$",re.I)
#BOUNDPROC_RE = re.compile("^(generic|procedure)\s*(\([^()]*\))?\s*(.*)\s*::\s*(\w.*)",re.I)
#COMMON_RE = re.compile("^common(?:\s*/\s*(\w+)\s*/\s*|\s+)(\w+.*)",re.I)
#COMMON_SPLIT_RE = re.compile("\s*(/\s*\w+\s*/)\s*",re.I)
#FINAL_RE = re.compile("^final\s*::\s*(\w.*)",re.I)
#USE_RE = re.compile("^use(?:\s*(?:,\s*(?:non_)?intrinsic\s*)?::\s*|\s+)(\w+)\s*($|,.*)",re.I)
#ARITH_GOTO_RE = re.compile("go\s*to\s*\([0-9,\s]+\)",re.I)
#CALL_RE = re.compile("(?:^|[^a-zA-Z0-9_% ]\s*)(\w+)(?=\s*\(\s*(?:.*?)\s*\))",re.I)


class FortranKissParser(object):
    """
    Parse fortran code.
    """
    # PROGRAM [name]
    # END [PROGRAM [name]]
    RE_PROG_START = re.compile(r"program\s*(?P<name>\w*)\Z", re.I)
    RE_PROG_END = re.compile(r"end(\s*program\s*(?P<name>\w*)|)\Z", re.I)

    # MODULE <name>
    # END [MODULE [name]]
    RE_MOD_START= re.compile(r"module\s+(?P<name>\w+)\Z", re.I)
    RE_MOD_END = re.compile(r"end(\s*module\s*(?P<name>\w*)|)\Z", re.I)

    # [<prefix>] <SUBROUTINE> <name> [(<args>)] [<suffix>]
    # END [SUBROUTINE [name]]
    RE_SUB_START = re.compile(r'(?P<prefix>(recursive|pure|elemental|\s)*)subroutine\s*(?P<name>\w+)', re.I)
    RE_SUB_END = re.compile(r'end(\s*subroutine\s*(?P<name>\w*)|)\Z', re.I)

    # [<prefix>] FUNCTION <name> ([<dummy-arg-list>]) [<suffix>]
    # END [FUNCTION [name]]
    RE_FUNC_START = re.compile(r"""
(?P<prefix>
(
recursive | pure | elemental |
logical | integer | integer(\s*\(.+\)\s*)? |
double\s+precision | real(\s*\(.+\))? |
complex(\s*\(.+\))? |
character\s*\(len=\w+\s*\) |
type\s*\(\w+\)
)
\s+
)*
\s*function\s+(?P<name>\w+)\s*""",
re.I | re.VERBOSE)
    #RE_FUNC_START = re.compile('^[ \t]*(([^!\'"\n]*?)function)',re.I)
    #RE_FUNC_START = re.compile("^(?:(.+?)\s+)?function\s+(\w+)\s*(\([^()]*\))?(?=(?:.*result\s*\(\s*(\w+)\s*\))?)(?=(?:.*bind\s*\(\s*(.*)\s*\))?).*$", re.I)
    RE_FUNC_END = re.compile(r"end(\s*function\s*(?P<name>\w*)|)\Z", re.I)

    RE_INTERFACE_START = re.compile("(abstract\s+)?interface(?:\s+(\S.+))?$", re.I)
    RE_INTERFACE_END = re.compile("end\s+interface(?:\s+(\S.+))?$", re.I)
    #RE_INTERFACE_END = re.compile("end\s+interface(?P<name>\s+\w*|)\Z", re.I)

    # [if ()] call <name> [([<dummy-arg-list>])]
    #RE_SUBCALL = re.compile("^(?:if\s*\(.*\)\s*)?call\s+(?P<name>\w+)\s*(?:\(\s*(.*?)\s*\))?\Z", re.I)
    RE_SUBCALL = re.compile("^(?:if\s*\(.*\)\s*)?call\s+(?P<name>\w+)", re.I)

    def __init__(self, verbose=0):
        self.verbose = verbose

    #@staticmethod
    #def rstrip_comment(s):
    #    for i in range(len(s), -1, -1)
    #        if s[i] == "!"

    def parse_file(self, path):
        with open(path, "rt") as fh:
            return self.parse_lines(fh.readlines(), path=path)

    #def parse_string(self, s):
    #    return self.parse_lines(s.splitlines())

    def parse_lines(self, lines, path=None):
        lines = [l.strip().lower() for l in lines]

        num_doclines, num_f90lines = 0, 0
        preamble, stack, includes  = [], [], []
        ancestor = None
        inblock = None

        while lines:
            line = lines.pop(0)
            if not line: continue
            #print(line)

            # Count number of comments and code line
            # Inlined comments are not counted (also because I don't like them)
            if line.startswith("!"):
                if not stack or (stack and stack[-1][1] != "open"):
                    preamble.append(line)
                if line.replace("!", ""): num_doclines += 1
                continue
            else:
                num_f90lines += 1

            # Interface declaration.
            if self.RE_INTERFACE_START.match(line):
                if self.verbose: print("begin interface", line)
                inblock = "interface"
                continue
                #self.consume_interface_block()

            if inblock == "interface":
                if self.verbose: print("in interface", line)
                if self.RE_INTERFACE_END.match(line):
                    inblock = None
                continue

            # Handle include statement (CPP or Fortran version).
            if line.startswith("#include ") or line.startswith("include "):
                what = line.split()[1].replace("'", "").replace('"', "")
                if stack:
                    stack[-1][0].includes.append(what)
                else:
                    includes.append(what)
                continue

            # Find use statements and the corresponding module
            # TODO in principle one could have use A; use B
            if line.startswith("use "):
                smod = line.split()[1].split(",")[0].lower()
                if smod.startswith("interfaces_"): continue
                # Remove comment at the end of the line if present.
                i = smod.find("!")
                if i != -1: smod = smod[:i]
                if self.verbose > 1: print("Found used module:", smod)
                # TODO
                #if stack:
                stack[-1][0].local_uses.append(smod)
                #uses.append(smod)
                #else:
                continue

            if self.simple_match(line, "contains"):
                assert stack
                ancestor = stack[-1][0]
                if self.verbose > 1: print("Setting ancestor to:", repr(ancestor))
                continue

            # Extract name from (module | program).
            m = self.RE_PROG_START.match(line)
            if m:
                name = m.group("name")
                if self.verbose > 1: print("Entering program:", name)
                if not name:
                    raise ValueError("Cannot find program name in line `%s`" % line)
                stack.append([Program(name, ancestor, preamble, path=path), "open"])
                preamble = []
                continue

            m = self.RE_MOD_START.match(line)
            if m:
                name = m.group("name")
                if self.verbose: print("Entering module:", name)
                # Ignore module procedure
                #if name == "procedure": continue
                # TODO
                #assert ancestor is None
                stack.append([Module(name, ancestor, preamble, path=path), "open"])
                preamble = []
                continue

            # Subroutine declaration.
            m = self.RE_SUB_START.match(line)
            #if not m: m = self.RE_FUNC_START.match(line)
            if m:
                subname = m.group("name")
                if not subname:
                    raise ValueError("Cannot find procedure name in line `%s`" % line)
                if self.verbose: print("Found subroutine:", subname, "in line:\n\t", line)
                stack.append([Subroutine(subname, ancestor, preamble, path=path), "open"])
                preamble = []
                continue

            # Function declaration
            m = self.RE_FUNC_START.match(line)
            if m:
                func_name = m.group("name")
                if not func_name:
                    raise ValueError("Cannot find procedure name in line:\n\t%s" % line)
                stack.append([Function(func_name, ancestor, preamble, path=path), "open"])
                preamble = []
                continue

            # Invokations of Fortran functions are difficult to handle
            # without inspecting locals so we only handle explicit calls to routines.
            m = self.RE_SUBCALL.match(line)
            if m:
                subname = m.group("name")
                assert subname
                stack[-1][0].children.append(subname)
                continue

            #isend, ftype, name = self.maybe_end_procedure(line):
            #if isend:
            if line.startswith("end "):
                tokens = line.split()
                end_ftype = None if len(tokens) == 1 else tokens[1]
                if end_ftype not in {"program", "function", "subroutine", "module", None}: continue
                try:
                    end_name = tokens[2]
                except IndexError:
                    # TODO
                    #print("index error in:", line)
                    end_name = None

                #if end_ftype == "module": print("got end", end_name)
                #print(line, end_name)
                #last_name = stack[-1][0].name
                #if last_name != end_name:
                #    print("WARNING: end_name:", end_name, " != last_name:", last_name)
                #    print("line", line, "end_ftype", end_ftype)

                #self.close_stack_entry(end_name, end_ftype)
                if end_name is not None:
                    # Close the last entry in the stack with name == end_name.
                    for item in reversed(stack):
                        if item[0].name == end_name:
                            item[1] = "closed"
                            break
                    else:
                        raise RuntimeError("Cannot find end_name `%s` in stack:\n%s" % (
                            end_name, pformat([s[0].name for s in stack])))
                else:
                    # Close the last entry in the stack with end_ftype.
                    if end_ftype is not None:
                        for item in reversed(stack):
                            if item[0].ftype == end_ftype:
                                item[1] = "closed"
                                break
                        else:
                            raise RuntimeError("Cannot find end_ftype `%s` in stack:\n%s" % (
                                end_ftype, pformat([s[0].ftype for s in stack])))

                    else:
                        stack[-1][1] = "closed"

                #print("Closing procedure", stack[-1][0].name)
                if ancestor is not None and ancestor.name == end_name:
                    ancestor = ancestor.ancestor

        self.includes = sorted(set(includes))
        #self.uses = sorted(set(uses))

        # Extract data from stack.
        self.programs, self.modules, self.subroutines, self.functions = [], [], [], []
        while stack:
            p, status = stack.pop(0)
            if status != "closed":
                print("WARNING: unclosed", repr(p), status, "ancestor", repr(p.ancestor))

            # Sort entries here.
            p.local_uses = sorted(p.local_uses)
            p.children = sorted(set(p.children))

            if p.ancestor is not None:
                #print("Adding %s to ancestor %s" % (repr(p), repr(p.ancestor)))
                p.ancestor.contains.append(p)
                #p.ancestor.contains[p.name] = p
            else:
                if p.is_module: self.modules.append(p)
                elif p.is_subroutine: self.subroutines.append(p)
                elif p.is_function: self.functions.append(p)
                elif p.is_program: self.programs.append(p)
                else: raise ValueError("Don't know how to handle %s" % type(p))

        return self

    @staticmethod
    def simple_match(s, token):
        i = s.find("!")
        if i != -1: s = s[:i]
        return s.strip() == token

    #@staticmethod
    #def maybe_end_procedure(s):
    #    if not line.startswith("end"):
    #        return False, False, False
    #    tokens = line.split()
    #    what = tokens[1]
    #    assert what in ("program", "function", "subroutine", "module"):
    #    try:
    #        fname = tokens[2]
    #    except IndexError:
    #        fname = None

    #    return isend, ftype, name
