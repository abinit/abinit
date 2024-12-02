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
    Many regexs use `^` so we assume source lines have been already stripped.
    """
    # https://software.intel.com/en-us/fortran-compiler-18.0-developer-guide-and-reference-interface

    # A quoted expression.
    RE_QUOTED = re.compile(r"(\".*?\")|('.*?')")

    RE_F90COMMENT = re.compile(r"\s*!(?P<value>.*)")
    RE_OMP_SENTINEL = re.compile(r"\s*!\$OMP", re.I)

    RE_SEARCH_PROC = re.compile(r"(subroutine|function|program)", re.I)

    # We don't allow plain `end`.
    # end subroutine name is mandatory
    # Perhaps this is not standard-compliant
    # but it makes the implementation more robust and the source code more readable.
    RE_PROC_END = re.compile(r"^end\s*(?P<proc_type>subroutine|function|program)\s+(?P<name>\w*)", re.I)

    # PROGRAM [name] && END [PROGRAM [name]]
    # NB: we enforce PROGRAM and NAME
    RE_PROG_START = re.compile(r"^program\s+(?P<name>\w+)", re.I)
    RE_PROG_END = re.compile(r"^end\s+program\s+(?P<name>\w+)", re.I)

    # MODULE <name> && END [MODULE [name]]
    # NB: we enforce MODULE and NAME in END
    RE_MOD_START= re.compile(r"^module\s+(?P<name>\w+)", re.I)
    RE_MOD_END = re.compile(r"^end\s+module\s+(?P<name>\w+)", re.I)

    #[ MODULE ] PROCEDURE <procedure-name-list>
    #re.compile(r"(module\s*|)procedure\s(?P<namelist>\w+)", re.I)

    # [<prefix>] <SUBROUTINE> <name> [(<args>)] [<suffix>]
    # END [SUBROUTINE [name]]
    # NB: we enforce name and name
    RE_SUB_START = re.compile(r'(?P<prefix>(recursive|pure|elemental|\s)*)subroutine\s*(?P<name>\w+)', re.I)
    RE_SUB_END = re.compile(r'^end\s+subroutine\s+(?P<name>\w+)', re.I)
    RE_SUB_ARGS = re.compile(r"\((?P<args>.*?)\)")

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

    # Enforce function and name
    RE_FUNC_END = re.compile(r"^end\s+function\s*(?P<name>\w+)", re.I)
    #result_re = re.compile(r'result\s*\((.*?)\)', re.I)

    # INTERFACE [generic-spec]
    # [interface-body]...
    # [[MODULE] PROCEDURE name-list]...
    # END INTERFACE [generic-spec]
    RE_INTERFACE_START = re.compile(r"(abstract\s+)?interface\s*(?P<name>\w*)", re.I)
    RE_INTERFACE_END = re.compile(r"^end\s+interface\s*(?P<name>\w*)", re.I)

    # [if ()] call <name> [([<dummy-arg-list>])]
    #RE_SUBCALL = re.compile("^(?:if\s*\(.*\)\s*)?call\s+(?P<name>\w+)\s*(?:\(\s*(.*?)\s*\))?\Z", re.I)
    RE_SUBCALL = re.compile("^(?:if\s*\(.*\)\s*)?call\s+(?P<name>\w+)", re.I)

    # TYPE [[,attr-list] :: ] name [(type-param-name-list)]
    #    [type-param-def-stmts]
    #    [PRIVATE statement or SEQUENCE statement]. . .
    #    [component-definition]. . .
    #    [type-bound-procedure-part]
    # END TYPE [name]
    # NB: Enforcing name in END
    #RE_TYPE_START = re.compile(r'^type(?:\s+|\s*(,.*)?::\s*)(?P<name>\w+)', re.I)
    RE_TYPE_START = re.compile(r'^type(?P<attribs>(?:\s+|\s*(,.*)?::\s*))(?P<name>\w+)', re.I)
    RE_TYPE_END = re.compile(r'^end\s+type\s+(?P<name>\w+)', re.I)

    RE_PUB_OR_PRIVATE = re.compile(r"^(?P<name>public|private)\s*(\!+\s*\w*|\Z)", re.I)

    #public = re.compile('(^public$)|(^public\s*(\w+)\s*$)|(^public\s*::\s*(\w+)(\s*,\s*\w+)*$)', re.I)
    #private = re.compile('(^private$)|(^private\s*(\w+)\s*$)|(^private\s*::\s*(\w+)(\s*,\s*\w+)*$)', re.I)

    RE_CONTAINS = re.compile(r"^(contains)\s*(\!+\s*\w*|\Z)", re.I)
    #RE_PUB_OR_PRIVATE_MODPROC = re.compile((r"^(?P<name>public|private)\s*(\!+\s*\w*|\Z)", re.I)

    # Continuation line:
    # To be improved with Lookahead/Lookbehind
    RE_CONTLINE_START = re.compile(r"^(?P<prefix>[^&\!]+)&(?P<postfix>\s*(\!+.*)?)", re.I)

    RE_CONTLINE_NEXT = re.compile(r"^&?(?P<value>.*?)&\s*(\!+.*)?", re.I)

    #Detection of character as style f77 + f95
    #re_character_f77 = re.compile('[ ]*character[ ]*([*]?[(][^)]+[)]|[*][0-9]+|)', re.I)
    #Detect a group after =
    #re_equal = re.compile('[ ]*=.*')

    # For character(len=xxx) or character(len=*)
    RE_CHARACTER_DEC = re.compile(r'^character\s*\(\s*len\s*=\s*(?P<len>(\*|\w+))\)', re.I)

    # To search for intent.
    RE_INTENT = re.compile(r",\s+intent\s*\(\s*(?P<value>in|out|inout)\s*\)", re.I)

    # For type(xxx)
    RE_TYPECLASS_DEC = re.compile('^(?P<ftype>(type|class))\s*\(\s*(?P<name>\w+)\s*\)', re.I)

    RE_NUMBOOL_DEC = re.compile(r"""
^(?P<ftype>(integer|real|double\s*precision|complex|double\s*complex|logical))\s*
(\(\s*(kind\s*=\s*)?(?P<kind>\w*)\s*\))?
""",
re.I | re.VERBOSE)

    # Include statement (CPP and F version).
    RE_INCLUDE = re.compile("#?include\s+(?P<path>.+)", re.I)

    #re_cpp_begin= re.compile("^[ ]*#[ \t]*if[ \t]*defined[ \t]*"+my_active_macro)
    #re_cpp_end  = re.compile("^[ ]*#[ \t]*endif")
    #re_define   = re.compile("^[ ]*#[ \t]*define[ \t]+")

    #RE_ASSUMED_SHAPE = re.compile("\(\s+:[,:\s]?\)", re.I)
    #re_result = re.compile('result[(](?P<result>[^)]+)[)]')

    #Use statement
    #re_use = re.compile('^[ \t]*use[ ]*(?P<name>\w+)', re.MULTILINE+re.IGNORECASE)
    #re_use_prefix = re.compile('^[ \t]*use[ ]*'+prefix+'.*?\n', re.MULTILINE+re.IGNORECASE)

    #RE_POINTS_TO = re.compile("\s*=>\s*", re.I)
    #MODPROC_RE = re.compile("^(module\s+)?procedure\s*(?:::|\s)\s*(\w.*)$", re.I)
    #FINAL_RE = re.compile("^final\s*::\s*(\w.*)", re.I)
    #VARIABLE_STRING = "^(integer|real|double\s*precision|character|complex|logical|type(?!\s+is)|class(?!\s+is|\s+default)|procedure|enumerator{})\s*((?:\(|\s\w|[:,*]).*)$"


FORTRAN_INTRINSICS = {
    'abort','abs','abstract','access','achar','acos','acosh','adjustl',
    'adjustr','aimag','aint','alarm','all','allocatable','allocate',
    'allocated','and','anint','any','asin','asinh','assign','associate',
    'associated','asynchronous','atan','atan2','atanh','atomic_add',
    'atomic_and','atomic_cas','atomic_define','atomic_fetch_add',
    'atomic_fetch_and','atomic_fetch_or','atomic_fetch_xor','atomic_or',
    'atomic_ref','atomic_xor','backtrace','backspace','bessel_j0',
    'bessel_j1','bessel_jn','bessel_y0','bessel_y1','bessel_yn','bge',
    'bgt','bind','bit_size','ble','block','block data','blt','btest',
    'c_associated','c_f_pointer','c_f_procpointer','c_funloc','c_loc',
    'c_sizeof','cabs','call','case','case default','cdabs','ceiling',
    'char','character','chdir','chmod','class','close','cmplx',
    'codimension','co_broadcast','co_max','co_min','co_reduce','co_sum',
    'command_argument_count','common','compiler_options',
    'compiler_version','complex','concurrent','conjg','contains',
    'contiguous','continue','cos','cosh','count','cpu_time','critical',
    'cshift','cycle','data','ctime','dabs','date_and_time','dble',
    'dcmplx','deallocate','deferred','digits','dim','dimension','do',
    'while','dlog','dlog10','dmax1','dmin1',
    'dot_product','double precision','dprod','dreal','dshiftl','dshiftr',
    'dsqrt','dtime','elemental','else','else if','elseif','elsewhere',
    'end','end associate','end block','end block data','end critical',
    'end do','end enum','end forall','end function','end if',
    'end interface','end module','end program','end select',
    'end submodule','end subroutine','end type','end where','endfile',
    'endif','entry','enum','enumerator','eoshift','epsilon',
    'equivalence','erf','erfc','erfc_scaled','etime','error stop',
    'execute_command_line','exit','exp','exponent','extends',
    'extends_type_of','external','fget','fgetc','final','findloc',
    'fdate','floor','flush','fnum','forall','format','fput','fputc',
    'fraction','function','free','fseek','fstat','ftell','gamma',
    'generic','gerror','getarg','get_command','get_command_argument',
    'getcwd','getenv','get_environment_variable','go to','goto','getgid',
    'getlog','getpid','getuid','gmtime','hostnm','huge','hypot','iabs',
    'iachar','iall','iand','iany','iargc','ibclr','ibits','ibset','ichar',
    'idate','ieee_class','ieee_copy_sign','ieee_get_flag',
    'ieee_get_halting_mode','ieee_get_rounding_mode','ieee_get_status',
    'ieee_get_underflow_mode','ieee_is_finite','ieee_is_nan',
    'ieee_is_negative','ieee_is_normal','ieee_logb','ieee_next_after',
    'ieee_rem','ieee_rint','ieee_scalb','ieee_selected_real_kind',
    'ieee_set_flag','ieee_set_halting_mode','ieee_set_rounding_mode',
    'ieee_set_status','ieee_support_datatype','ieee_support_denormal',
    'ieee_support_divide','ieee_support_flag','ieee_support_halting',
    'ieee_support_inf','ieee_support_io','ieee_support_nan',
    'ieee_support_rounding','ieee_support_sqrt','ieee_support_standard',
    'ieee_support_underflow_control','ieee_unordered','ieee_value',
    'ieor','ierrno','if','imag','image_index','implicit',
    'implicit none','import','include','index','inquire','int','integer',
    'intent','interface','intrinsic','int2','int8','ior','iparity',
    'irand','is','is_contiguous','is_iostat_end','is_iostat_eor',
    'isatty','ishft','ishftc','isnan','itime','kill','kind','lbound',
    'lcobound','leadz','len','len_trim','lge','lgt','link','lle','llt',
    'lock','lnblnk','loc','log','log_gamma','log10','logical','long',
    'lshift','lstat','ltime','malloc','maskl','maskr','matmul','max',
    'max0','maxexponent','maxloc','maxval','mclock','mclock8','merge',
    'merge_bits','min','min0','minexponent','minloc','minval','mod',
    'module','module procedure','modulo','move_alloc','mvbits','namelist',
    'nearest','new_line','nint','non_overridable','none','nopass','norm2',
    'not','null','nullify','num_images','only','open','or','operator',
    'optional','pack','parameter','parity','pass','pause','pointer',
    'perror','popcnt','poppar','precision','present','print','private',
    'procedure','product','program','protected','public','pure','radix',
    'ran','rand','random_number','random_seed','range','rank','read',
    'real','recursive','rename','repeat','reshape','result','return',
    'rewind','rewrite','rrspacing','rshift','same_type_as','save',
    'scale','scan','secnds','second','select','select case','select type',
    'selected_char_kind','selected_int_kind','selected_real_kind',
    'sequence','set_exponent','shape','shifta','shiftl','shiftr','sign',
    'signal','sin','sinh','size','sizeof','sleep','spacing','spread',
    'sqrt','srand','stat','stop','storage_size','submodule','subroutine',
    'sum','sync all','sync images','sync memory','symlnk','system',
    'system_clock','tan','tanh','target','then','this_image','time',
    'time8','tiny','trailz','transfer','transpose','trim','ttynam',
    'type','type_as','ubound','ucobound','umask','unlock','unlink',
    'unpack','use','value','verify','volatile','wait','where','while',
    'write','xor','zabs',
}
