"""
Regular expressions for Fortran code.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import re

# This one is problematic. The firs one is stricter but it does not play well with the other code!!
# I use the old one used by abilint
#RE_FUNC_END = re.compile('^[ \t]*end[ \t]*(function\n)',re.I).match

#Detect call (not after an if or ;)
#re_call = re.compile('(^[ \t]*call[ ]*)(\w+)', re.MULTILINE+re.I)

# Taken from https://github.com/cmacmackin/ford/blob/master/ford/sourceform.py
#VAR_TYPE_STRING = "integer|real|double\s*precision|character|complex|logical|type|class|procedure|enumerator"
#VAR_TYPE_RE = re.compile(r"(?P<ftype>(integer|real|double\s*precision|character|complex|logical|type|class|procedure|enumerator))", re.I)
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

class HasRegex(object):
    """
    Mixin class providing regular expressions used to analyze Fortran code.
    """

    # https://software.intel.com/en-us/fortran-compiler-18.0-developer-guide-and-reference-interface

    # PROGRAM [name]
    # END [PROGRAM [name]]
    RE_PROG_START = re.compile(r"^program\s*(?P<name>\w*)\Z", re.I)
    RE_PROG_END = re.compile(r"^end(\s*program\s*(?P<name>\w*)|)\Z", re.I)

    # MODULE <name>
    # END [MODULE [name]]
    RE_MOD_START= re.compile(r"^module\s+(?P<name>\w+)\Z", re.I)
    RE_MOD_END = re.compile(r"^end(\s*module\s*(?P<name>\w*)|)\Z", re.I)

    # [<prefix>] <SUBROUTINE> <name> [(<args>)] [<suffix>]
    # END [SUBROUTINE [name]]
    RE_SUB_START = re.compile(r'(?P<prefix>(recursive|pure|elemental|\s)*)subroutine\s*(?P<name>\w+)', re.I)
    RE_SUB_END = re.compile(r'^end(\s*subroutine\s*(?P<name>\w*)|)\Z', re.I)

    # [<prefix>] FUNCTION <name> ([<dummy-arg-list>]) [<suffix>]
    # END [FUNCTION [name]]
    RE_FUNC_START = re.compile(r"""
(?P<prefix>
(
recursive (.*) | pure (.*) | elemental (.*) |
logical | integer | integer(\s*\(.+\)\s*)? |
double\s+precision | real(\s*\(.+\))? | complex(\s*\(.+\))? |
character\s*\(\s*len=\w+\s*\) | type\s*\(\s*\w+\s*\)
)
\s+
)*
\s*function\s+(?P<name>\w+)\s*""",
re.I | re.VERBOSE)
    #RE_FUNC_START = re.compile('^[ \t]*(([^!\'"\n]*?)function)',re.I)
    #RE_FUNC_START = re.compile("^(?:(.+?)\s+)?function\s+(\w+)\s*(\([^()]*\))?(?=(?:.*result\s*\(\s*(\w+)\s*\))?)(?=(?:.*bind\s*\(\s*(.*)\s*\))?).*$", re.I)
    RE_FUNC_END = re.compile(r"^end(\s*function\s*(?P<name>\w*)|)\Z", re.I)

    # INTERFACE [generic-spec]
    # [interface-body]...
    # [[MODULE]PROCEDURE name-list]...
    # END INTERFACE [generic-spec]
    #RE_INTERFACE_START = re.compile("^(abstract\s+)?interface(?:\s+(\S.+))?$", re.I)
    #RE_INTERFACE_END = re.compile("^end\s+interface(?:\s+(\S.+))?$", re.I)

    RE_INTERFACE_START = re.compile(r"(abstract\s+)?interface\s*(?P<name>\w*)", re.I)
    RE_INTERFACE_END = re.compile(r"end\s+interface\s*(?P<name>\w*)", re.I)

    # [if ()] call <name> [([<dummy-arg-list>])]
    #RE_SUBCALL = re.compile("^(?:if\s*\(.*\)\s*)?call\s+(?P<name>\w+)\s*(?:\(\s*(.*?)\s*\))?\Z", re.I)
    RE_SUBCALL = re.compile("^(?:if\s*\(.*\)\s*)?call\s+(?P<name>\w+)", re.I)

    # TYPE [[,attr-list] :: ] name [(type-param-name-list)]
    #    [type-param-def-stmts]
    #    [PRIVATE statement or SEQUENCE statement]. . .
    #    [component-definition]. . .
    #    [type-bound-procedure-part]
    # END TYPE [name]
    #RE_TYPE_START = re.compile(r'type\s*(|.*::)\s*(?P<name>\w+)', re.I)
    #RE_TYPE_END = re.compile(r'^end\s+type', re.I)
    RE_TYPE_START = re.compile(r'^type(?:\s+|\s*(,.*)?::\s*)(?P<name>\w+)\Z', re.I)
    RE_TYPE_END = re.compile(r'^end(\s*type\s*(?P<name>\w*)|)\Z', re.I)

    RE_PUB_OR_PRIVATE = re.compile(r"^(?P<name>public|private)\s*(\!+\s*\w*|\Z)", re.I)
    #RE_CONTAINS = re.compile(r"^(contains)\s*(\!+\s*\w*|\Z)", re.I)

    #re_assumed_shape = re.compile("\(:.*\)", re.I)
    #Ampersand (continuation line)
    #re_amp = re.compile('[ ]*&[ ]*')

    #Detection of character as style f77 + f95
    #re_character_f77 = re.compile('[ ]*character[ ]*([*]?[(][^)]+[)]|[*][0-9]+|)', re.I)
    #Detect a group after =
    #re_equal = re.compile('[ ]*=.*')

    # For character(len=xxx) or character(len=*)
    RE_CHARACTER_DEC = re.compile(r'character\s*\(\s*len\s*=\s*(?P<len>(\*|\w+))\)', re.I)

    #For complex(kind=dp) or complex(kind(dp))
    #re_complex = re.compile('complex[(](kind[=(])?(?P<kind>[^)]+)[)]+', re.I)
    #For real(kind=dp) or real(kind(dp))
    #re_real = re.compile('real[(](kind[=(])?(?P<kind>[^)]+)[)]+', re.I)
    # For type(xxx)
    RE_TYPECLASS_DEC = re.compile('(?P<ftype>(type|class))\s*\(\s*(?P<name>\w+)\s*\)', re.I)

    #Continuation
    #re_continuation = re.compile("&[ \t]*(![^\n]*)?\n")

    #Include command
    #re_include = re.compile('^[ ]*include.*?\n',re.MULTILINE+re.IGNORECASE)

    #Detect robodoc header (!***)
    #re_robodoc_header = re.compile("!!\*\*\*\*[a-z]\*")
