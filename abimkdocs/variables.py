# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

import sys
import os
import json

from collections import OrderedDict, defaultdict
from itertools import groupby

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


def is_string(s):
    """True if s behaves like a string (duck typing test)."""
    try:
        s + " "
        return True
    except TypeError:
        return False


def list_strings(arg):
    """
    Always return a list of strings, given a string or list of strings as input.

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

def splitall(path):
    """Return list with all components of a path."""
    allparts = []
    while True:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path: # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts


# Unit names supported in Abinit input.
ABI_UNITS = [
    'au',
    'Angstr',
    'Angstrom',
    'Angstroms',
    'Bohr',
    'Bohrs',
    'eV',
    'Ha',
    'Hartree',
    'Hartrees',
    'K',
    'Ry',
    'Rydberg',
    'Rydbergs',
    'T',
    'Tesla',
    'Second',
    'S',
    'Sec',
]

# Operators supported by parser
ABI_OPS = ['sqrt', 'end', '*', '/']


# List of strings with possible character of variables.
# This is the reference set that will checked against the input
# given by the developer in the variables_CODENAME modules.
ABI_CHARACTERISTICS = [
    "DEVELOP",
    "EVOLVING",
    "ENERGY",
    "INPUT_ONLY",
    "INTERNAL_ONLY",
    "LENGTH",
    "MAGNETIC_FIELD",
    "NO_MULTI",
]

# external parametersare not input variables,
# but are used in the documentation of other variables.
ABI_EXTERNAL_PARAMS = OrderedDict([
    ("AUTO_FROM_PSP", "Means that the value is read from the PSP file"),
    ("CUDA", "True if CUDA is enabled (compilation)"),
    ("ETSF_IO", "True if NetCDF is enabled (compilation)"),
    ("FFTW3", "True if FFTW3 is enabled (compilation)"),
    ("MPI_IO", "True if MPI_IO is enabled (compilation)"),
    ("NPROC", "Number of processors used for Abinit"),
    ("PARALLEL", "True if the code is compiled with MPI"),
    ("SEQUENTIAL", "True if the code is compiled without MPI"),
])

# List of topics
ABI_TOPICS = [
    "Abipy",
    "APPA",
    "Artificial",
    "AtomManipulator",
    "AtomTypes",
    "Bader",
    "Band2eps",
    "Berry",
    "BandOcc",
    "BoundingProcess",
    "BSE",
    "ConstrainedDFT",
    "ConstrainedPol",
    "Control",
    "Coulomb",
    "CrossingBarriers",
    "CRPA",
    "crystal",
    "DFT+U",
    "DeltaSCF",
    "DensityPotential",
    "Dev",
    "DFPT",
    "DMFT",
    "DynamicsMultibinit",
    "EffectiveMass",
    "EFG",
    "Elastic",
    "ElPhonInt",
    "ElPhonTransport",
    "ElecDOS",
    "ElecBandStructure",
    "FileFormats",
    "FitProcess",
    "ForcesStresses",
    "FrequencyMeshMBPT",
    "GeoConstraints",
    "GeoOpt",
    "Git",
    "GSintroduction",
    "GW",
    "GWls",
    "Hybrids",
    "k-points",
    "LatticeModel",
    "LDAminushalf",
    "Longwave" ,
    "LOTF",
    "MagField",
    "MagMom",
    "MolecularDynamics",
    "Macroave",
    "multidtset",
    "nonlinear",
    "Optic",
    "Output",
    "parallelism",
    "PAW",
    "PIMD",
    "Planewaves",
    "Phonons",
    "PhononBands",
    "PhononWidth",
    "PortabilityNonRegression",
    "positron",
    "printing",
    "PseudosPAW",
    "q-points",
    "RandStopPow",
    "Recursion",
    "RPACorrEn",
    "SCFControl",
    "SCFAlgorithms",
    "SelfEnergy",
    "SmartSymm",
    "spinpolarisation",
    "SpinDynamicsMultibinit",
    "STM",
    "Susceptibility",
    "TDDFT",
    "Tdep",
    "TDepES",
    "Temperature",
    "TransPath",
    "TuningSpeed",
    "Unfolding",
    "UnitCell",
    "vdw",
    "Verification",
    "Wannier",
    "Wavelets",
    "xc",
]

# Relevance associated to the topic
ABI_RELEVANCES = OrderedDict([
    ("compulsory", 'Compulsory input variables'),
    ("basic", 'Basic input variables'),
    ("useful", 'Useful input variables'),
    ("internal", 'Relevant internal variables'),
    ("prpot", 'Printing input variables for potentials'),
    ("prfermi", 'Printing input variables for fermi level or surfaces'),
    ("prden", 'Printing input variables for density, eigenenergies, k-points and wavefunctions'),
    ("prgeo", 'Printing input variables for geometry'),
    ("prdos", "Printing DOS-related input variables"),
    ("prgs", 'Printing other ground-state input variables'),
    ("prngs", 'Printing non-ground-state input variables'),
    ("prmisc", 'Printing miscellaneous files'),
    ("expert",  'Input variables for experts'),
])


class Variable(object):
    """
    This object gathers information about a single variable. name, associated topics, description etc
    It's constructed from the variables_CODENANE.py modules but client code usually
    interact with variables via the :class:`VarDatabase` dictionary.
    """
    def __init__(self,
                 abivarname=None,
                 varset=None,
                 vartype=None,
                 topics=None,
                 dimensions=None,
                 defaultval=None,
                 mnemonics=None,
                 characteristics=None,
                 excludes=None,
                 requires=None,
                 commentdefault=None,
                 commentdims=None,
                 added_in_version=None,
                 alternative_name=None,
                 text=None,
                ):
        """
        Args:
            abivarname (str): Name of the variable (including @code if not abinit e.g asr@anaddb).
                Required
            varset (str): The group this variable belongs to (could be code if code has no group).
                Required
            vartype (str): The type of the variable. Required
            topics (list): List of strings with topics. Required
            dimensions: List of strings with dimensions or "scalar". Required.
            defaultval: Default value. None if no default is provided. Other possibilities are ...
                Either constant number, formula or another variable
            mnemonics (str): Mnemonic string (required).
            characteristics (list): List of characteristics or None
            excludes (str): String with variables that are excluded if this variable is given.
            requires (str): String with variables that are required.
            commentdefault=None,
            commentdims=None,
            added_in_version (str): String with the Abinit version in which this variable was added.
            alternative_name: alias name (used if a new variable with a different name was introduced, in place
                of of an old variable that is still supported.
            text: markdown string with documentation. Required.
        """
        self.abivarname = abivarname
        self.varset = varset
        self.vartype = vartype
        self.topics = topics
        self.dimensions = dimensions
        self.defaultval = defaultval
        self.mnemonics = mnemonics
        self.characteristics = characteristics
        self.excludes = excludes
        self.requires = requires
        self.commentdefault = commentdefault
        self.commentdims = commentdims
        self.added_in_version = added_in_version
        self.alternative_name = alternative_name
        self.text = my_unicode(text)

        errors = []
        for a in ("abivarname", "varset", "vartype", "topics", "dimensions", "added_in_version", "text"):
            if getattr(self, a) is None:
                errors.append("attribute %s is mandatory" % a)
        if errors:
            raise ValueError("Errors in %s:\n%s" % (self.abivarname, "\n".join(errors)))

    @lazy_property
    def name(self):
        """Name of the variable without the executable name."""
        return self.abivarname if "@" not in self.abivarname else self.abivarname.split("@")[0]

    @lazy_property
    def executable(self):
        """string with the name of the code associated to this variable."""
        if "@" in self.abivarname:
            code = self.abivarname.split("@")[1]
            assert code == self.varset
        else:
            code = "abinit"
        return code

    @lazy_property
    def website_url(self):
        """
        The absolute URL associated to this variable on the Abinit website.
        """
        # This is gonna be the official API on the server
        #docs.abinit.org/vardocs/CODENAME/VARNAME?version=8.6.2
        #return "https://docs.abinit.org/vardocs/%s/%s" % (self.executable, self.name)

        # For the time being, we have to use:
        # variables/eph/#asr
        # variables/anaddb#asr
        if self.executable == "abinit":
            return "https://docs.abinit.org/variables/%s#%s" % (self.varset, self.name)
        else:
            return "https://docs.abinit.org/variables/%s#%s" % (self.executable, self.name)

    @lazy_property
    def topic2relevances(self):
        """topic --> list of relevances"""
        assert self.topics is not None
        od = OrderedDict()
        for tok in self.topics:
            topic, relevance = [s.strip() for s in tok.split("_")]
            if topic not in od: od[topic] = []
            od[topic].append(relevance)
        return od

    @lazy_property
    def is_internal(self):
        """True if this is an internal variable."""
        return self.characteristics is not None and '[[INTERNAL_ONLY]]' in self.characteristics

    @lazy_property
    def wikilink(self):
        """Abinit wikilink."""
        return "[[%s:%s]]" % (self.executable, self.name)

    def __repr__(self):
        """Variable name + mnemonics"""
        return self.abivarname + "  <" + str(self.mnemonics) + ">"

    def to_string(self, verbose=0):
        """String representation with verbosity level `verbose`."""
        return "Variable " + str(self.abivarname) + " (default = " + str(self.defaultval) + ")"

    def __str__(self):
        return self.to_string()

    def __hash__(self):
        # abivarname is unique
        return hash(self.abivarname)

    def __eq__(self, other):
        if other is None: return False
        return self.abivarname == other.abivarname

    def __ne__(self, other):
        return not (self == other)

    @lazy_property
    def info(self):
        """String with extra info on the variable."""
        attrs = [
            "vartype", "characteristics",  "mnemonics", "dimensions", "defaultval",
            "abivarname", "commentdefault", "commentdims", "varset",
            "requires", "excludes",
            "added_in_version", "alternative_name",
            ]

        def astr(obj):
            return str(obj).replace("[[", "").replace("]]", "")

        d = {k: astr(getattr(self, k)) for k in attrs if getattr(self, k) is not None}
        return json.dumps(d, indent=4, sort_keys=True)

    def _repr_html_(self):
        """Integration with jupyter notebooks."""
        try:
            import markdown
        except ImportError:
            markdown = None

        if markdown is None:
            html = "<h2>Default value:</h2>" + my_unicode(self.defaultval) + "<br/><h2>Description</h2>" + self.text
            return html.replace("[[", "<b>").replace("]]", "</b>")
        else:
            md = self.text.replace("[[", "<b>").replace("]]", "</b>")
            return markdown.markdown("""
## Default value:
{defaultval}

## Description:
{md}
""".format(defaultval=my_unicode(self.defaultval), md=my_unicode(md)))

    def browse(self):
        """Open variable documentation in browser."""
        import webbrowser
        return webbrowser.open(self.website_url)

    @lazy_property
    def isarray(self):
        """True if this variable is an array."""
        return not (is_string(self.dimensions) and self.dimensions == "scalar")

    def depends_on_dimension(self, dimname):
        """
        True if variable is an array whose shape depends on dimension name `dimname`.

        Args: dimname: String of :class:`Variable` object.
        """
        if not self.isarray: return False
        if isinstance(dimname, Variable): dimname = dimname.name
        # This test is not very robust and can fail.
        # Assume no space between `[` and name (there should be a test for this...)
        key = "[[%s]]" % dimname
        for d in self.dimensions:
            if key in str(d): return True
        return False

    def html_link(self, label=None):
        """String with the URL of the web page."""
        label = self.name if label is None else label
        return '<a href="%s" target="_blank">%s</a>' % (self.website_url, label)

    def get_parent_names(self):
        """
        Return set of strings with the name of the parents
        i.e. the variables that are connected to this variable
        (either because they are present in dimensions on in requires).
        """
        #if hasattr(self, ...
        import re
        parent_names = []
        WIKILINK_RE = r'\[\[([\w0-9_ -]+)\]\]'
        # TODO
        #  parent = self[parent]
        # KeyError: "'nzchempot'
        #WIKILINK_RE = r'\[\[([^\[]+)\]\]'
        if isinstance(self.dimensions, (list, tuple)):
            for dim in self.dimensions:
                dim = str(dim)
                m = re.match(WIKILINK_RE, dim)
                if m:
                    parent_names.append(m.group(1))

        if self.requires is not None:
            parent_names.extend([m.group(1) for m in re.finditer(WIKILINK_RE, self.requires) if m])

        # Convert to set and remove possibile self-reference.
        parent_names = set(parent_names)
        parent_names.discard(self.name)
        return parent_names

    def internal_link(self, website, page_rpath, label=None, cls=None):
        """String with the website internal URL."""
        token = "%s:%s" % (self.executable, self.name)
        a = website.get_wikilink(token, page_rpath)
        cls = a.get("class") if cls is None else cls
        return '<a href="%s" class="%s">%s</a>' % (a.get("href"), cls, a.text if label is None else label)

    @staticmethod
    def format_dimensions(dimensions):
        """Pretty print dimensions."""
        if dimensions is None:
            s = ''
        elif dimensions == "scalar":
            s = 'scalar'
        else:
            #s = str(dimensions)
            if isinstance(dimensions, (list, tuple)):
                s = '('
                for dim in dimensions:
                    s += str(dim) + ','
                s = s[:-1]
                s += ')'
            else:
                s = str(dimensions)

        return s

    def to_abimarkdown(self, with_hr=True):
        """
        Return markdown string. Can use Abinit markdown extensions.
        """
        lines = []; app = lines.append

        app("## **%s** \n\n" % self.name)
        app("*Mnemonics:* %s  " % str(self.mnemonics))
        if self.characteristics:
            app("*Characteristics:* %s  " % ", ".join(self.characteristics))
        if self.topic2relevances:
            app("*Mentioned in topic(s):* %s  " % ", ".join("[[topic:%s]]" % k for k in self.topic2relevances))
        app("*Variable type:* %s  " % str(self.vartype))
        if self.dimensions:
            app("*Dimensions:* %s  " % self.format_dimensions(self.dimensions))
        if self.commentdims:
            app("*Commentdims:* %s  " % self.commentdims)
        app("*Default value:* %s  " % self.defaultval)
        if self.commentdefault:
            app("*Comment:* %s  " % self.commentdefault)
        if self.requires:
            app("*Only relevant if:* %s  " % str(self.requires))
        if self.excludes:
            app("*The use of this variable forbids the use of:* %s  " % self.excludes)
        app("*Added in version:* %s  " % self.added_in_version)

        # Add links to tests.
        if hasattr(self, "tests") and not self.is_internal:
            # Constitutes an usage report e.g.
            # Rarely used, in abinit tests [8/888], in tuto abinit tests [2/136].
            assert hasattr(self, "tests_info")
            tests_info = self.tests_info
            ratio_all = len(self.tests) / tests_info["num_all_tests"]
            frequency = "Rarely used"
            if ratio_all > 0.5:
                frequency = "Very frequently used"
            elif ratio_all > 0.01:
                frequency = "Moderately used"

            info = "%s, [%d/%d] in all %s tests, [%d/%d] in %s tutorials" % (
                frequency, len(self.tests), tests_info["num_all_tests"], self.executable,
                tests_info["num_tests_in_tutorial"], tests_info["num_all_tutorial_tests"], self.executable)

            # Use https://facelessuser.github.io/pymdown-extensions/extensions/details/
            # Truncate list of tests if we have more that `max_ntests` entries.
            count, max_ntests = 0, 20
            app('\n??? note "Test list (click to open). %s"' % info)
            tlist = sorted(self.tests, key=lambda t: t.suite_name)
            d = {}
            for suite_name, tests_in_suite in groupby(tlist, key=lambda t: t.suite_name):
                ipaths = [os.path.join(*splitall(t.inp_fname)[-4:]) for t in tests_in_suite]
                count += len(ipaths)
                d[suite_name] = ipaths

            for suite_name, ipaths in d.items():
                if count > max_ntests: ipaths = ipaths[:min(3, len(ipaths))]
                s = "- " + suite_name + ":  " + ", ".join("[[%s|%s]]" % (p, os.path.basename(p)) for p in ipaths)
                if count > max_ntests: s += " ..."
                app("    " + s)
            app("\n\n")

        # Add text with description.
        app(2 * "\n")
        # Replace all occurrences of [[name]] with **name** to reduce number of html links in docs
        new_text = self.text.replace("[[%s]]" % self.name, " **%s** " % self.name)
        app(new_text)
        if with_hr: app("* * *" + 2*"\n")

        return "\n".join(lines)

    def validate(self):
        """Validate variable. Raises ValueError if not valid."""
        errors = []
        eapp = errors.append

        try:
            svar = str(self)
        except Exception as exc:
            svar = "Unknown"
            eapp(str(exc))

        if self.abivarname is None:
            eapp("Variable `%s` has no name" % svar)

        if self.vartype is None:
            eapp("Variable `%s` has no vartype" % svar)
        elif not self.vartype in ("integer", "real", "string"):
            eapp("%s must have vartype in ['integer', 'real', 'string'].")

        if self.topics is None:
            eapp("%s does not have at least one topic and the associated relevance" % svar)

        for topic, relevances in self.topic2relevances.items():
            if topic not in ABI_TOPICS:
                eapp("%s delivers topic `%s` that does not belong to the allowed list" % (sname, topic))
            for relevance in relevances:
                if relevance not in ABI_RELEVANCES:
                    eapp("%s delivers relevance `%s` that does not belong to the allowed list" % (sname, relevance))

	# Compare the characteristics of this variable with the refs to detect possible typos.
        if self.characteristics is not None:
            if not isinstance(self.characteristics, list):
                eapp("The field characteristics of %s is not a list" % svar)
            else:
                for cat in self.characteristics:
                    if cat.replace("[[", "").replace("]]", "") not in ABI_CHARACTERISTICS:
                        eapp("The characteristics %s of %s is not valid" % (cat, svar))

        if self.dimensions is None:
            eapp("%s does not have a dimension. If it is a *scalar*, it must be declared so." % svar)
        else:
            if self.dimensions != "scalar":
                if not isinstance(self.dimensions, (list, ValueWithConditions)):
                    eapp('The dimensions field of %s is not a list neither a valuewithconditions' % svar)

        if self.varset is None:
            eapp('`%s` does not have a varset' % svar)
        #else:
        #    if not isinstance(self.varset, str) or self.varset not in ref_varset:
        #        print('The field varset of %s should be one of the valid varsets' % str(self))

        if len(self.name) > 25:
            eapp("Lenght of `%s` is longer than 25 characters." % self.name)

        if errors:
            raise ValueError("\n".join(errors))


class ValueWithUnit(object):
    """
    This type allows to specify values with units:
    """
    def __init__(self, value=None, units=None):
        self.value = value
        self.units = units

    def __str__(self):
        return str(self.value) + " " + str(self.units)

    def __repr__(self):
        return str(self)


class Range(object):
    """
    Specifies a range (start:stop:step)
    """
    start = None
    stop = None

    def __init__(self, start=None, stop=None):
        self.start = start
        self.stop = stop

    def isin(self, value):
        """True if value is in range."""
        isin = True
        if self.start is not None:
            isin = isin and (self.start <= self.value)
        if stop is not None:
            isin = isin and self.stop > self.value
        return str(self)

    def __repr__(self):
        # Add whitespace after `[` or before `]` to avoid [[[ and ]]] patterns
        # that enter into conflict with wikiling syntax [[...]]
        if self.start is not None and self.stop is not None:
            return "[ " + str(self.start) + " ... " + str(self.stop) + " ]"
        if self.start is not None:
            return "[ " + str(self.start) + "; ->"
        if self.stop is not None:
            return "<-;" + str(self.stop) + " ]"
        else:
            return None


class ValueWithConditions(dict):
    """
    Used for variables whose value depends on a list of conditions.

    .. example:

        ValueWithConditions({'[[paral_kgb]]==1': '6', 'defaultval': 2}),

        Means that the variable is set to 6 if paral_kgb == 1 else 2
    """
    def __repr__(self):
        s = ''
        for key in self:
            if key != 'defaultval':
                s += str(self[key]) + ' if ' + str(key) + ',\n'
        s += str(self["defaultval"]) + ' otherwise.\n'
        return s

    def __str__(self):
        return self.__repr__()


class MultipleValue(object):
    """
    Used for variables that can assume multiple values.
    This is the equivalent to the X * Y syntax in the Abinit parser.
    If X is null, it means that you want to do *Y (all Y)
    """
    def __init__(self, number=None, value=None):
        self.number = number
        self.value = value

    def __repr__(self):
        if self.number is None:
            return "*" + str(self.value)
        else:
            return str(self.number) + " * " + str(self.value)


def my_unicode(s):
    """Convert string to unicode (needed for py2.7 DOH!)"""
    return unicode(s) if sys.version_info[0] <= 2 else str(s)

##############
# Public API #
##############

_VARS = None

def get_codevars():
    """
    Return the database of variables indexed by code name and cache it.
    Main entry point for client code.
    """
    global _VARS
    if _VARS is None: _VARS = VarDatabase.from_pyfiles()
    return _VARS


class VarDatabase(OrderedDict):
    """
    This object stores the full set of input variables for all the Abinit executables.
    in a dictionary mapping the name of the code to a subdictionary of variables.
    """

    all_characteristics = ABI_CHARACTERISTICS
    all_external_params = ABI_EXTERNAL_PARAMS

    @classmethod
    def from_pyfiles(cls, dirpath=None):
        """
        Initialize the object from python modules inside dirpath.
        If dirpath is None, the directory of the present module is used.
        """
        if dirpath is None:
            dirpath = os.path.dirname(os.path.abspath(__file__))
        pyfiles = [os.path.join(dirpath, f) for f in os.listdir(dirpath) if
                   f.startswith("variables_") and f.endswith(".py")]
        new = cls()
        for pyf in pyfiles:
            vd = InputVariables.from_pyfile(pyf)
            new[vd.executable] = vd

        return new

    def iter_allvars(self):
        """Iterate over all variables. Flat view."""
        for vd in self.values():
            for var in vd.values():
                yield var

    def get_version_endpoints(self):
        """
        API used by the webser to serve the documentation of a variable given codename, varname, [version]:

            docs.abinit.org/vardocs/abinit/asr?version=8.6.2

        # asr@anaddb at /variables/anaddb#asr
        # asr@abinit at /variables/eph#asr
        # asr@abinit at /variables/abinit/eph#asr
        """
        code_urls = {}
        for codename, vard in self.items():
            code_urls[codename] = d = {}
            for vname, var in var.items():
                # This is the internal convention used to build the mkdocs site.
                d[vname] = "/variables/%s/%s#%s" % (codename, var.varset, var.name)
	# TODO: version and change mkdocs.yml
        return version, code_urls

    def update_json_endpoints(self, json_path, indent=4):
        """
        Update the json file with the mapping varname --> relative url
        used by the webserve to implement the `vardocs` API.
        """
        with open(json_path, "rt") as fh:
            oldd = json.load(fh)

        new_version, newd = self.get_version_endpoints()
        assert new_version not in oldd
        oldd[new_version] = newd
        with open(json_path, "wt") as fh:
            json.dump(oldd, fh, indent=indent)

    def _write_pymods(self, dirpath="."):
        """
        Internal method used to regenerate the python modules.
        """
        dirpath = os.path.abspath(dirpath)
        from pprint import pformat
        def nones2arg(obj, must_be_string=False):
            if obj is None:
                if must_be_string:
                    raise TypeError("obj must be string.")
                return None
            elif isinstance(obj, str):
                s = str(obj).rstrip()
                if "\n" in s: return '"""%s"""' % s
                if "'" in s: return '"%s"' % s
                if '"' in s: return "'%s'" % s
                return '"%s"' % s
            else:
                raise TypeError("%s: %s" % (type(obj), str(obj)))

        def topics2arg(obj):
            if isinstance(obj, str):
                if "," in obj:
                    obj = [s.strip() for s in obj.split(",")]
                else:
                    obj = [obj]
            if isinstance(obj, (list, tuple)): return pformat(obj)

            raise TypeError("%s: %s" % (type(obj), str(obj)))

        def dimensions2arg(obj):
            if isinstance(obj, str) and obj == "scalar": return '"scalar"'
            if isinstance(obj, (ValueWithUnit, MultipleValue, Range, ValueWithConditions)):
                return "%s(%s)" % (obj.__class__.__name__, pformat(obj.__dict__))
            if isinstance(obj, (list, tuple)): return pformat(obj)

            raise TypeError("%s, %s" % (type(obj), str(obj)))

        def defaultval2arg(obj):
            if obj is None: return obj
            if isinstance(obj, (ValueWithUnit, MultipleValue, Range, ValueWithConditions)):
                return "%s(%s)" % (obj.__class__.__name__, pformat(obj.__dict__))
            if isinstance(obj, (list, tuple)): return pformat(obj)
            if isinstance(obj, str): return '"%s"' % str(obj)
            if isinstance(obj, (int, float)): return obj

            raise TypeError("%s, %s" % (type(obj), str(obj)))

        for code in self:
            varsd = self[code]

            lines = ["""\
from __future__ import print_function, division, unicode_literals, absolute_import

from abimkdocs.variables import ValueWithUnit, MultipleValue, Range
ValueWithConditions = dict

Variable=dict\nvariables = ["""
]
            for name in sorted(varsd.keys()):
                var = varsd[name]
                text = '"""\n' + var.text.rstrip() + '\n"""'
                s = """\
Variable(
    abivarname={abivarname},
    varset={varset},
    vartype={vartype},
    topics={topics},
    dimensions={dimensions},
    defaultval={defaultval},
    mnemonics={mnemonics},
    characteristics={characteristics},
    excludes={excludes},
    requires={requires},
    commentdefault={commentdefault},
    commentdims={commentdims},
    added_in_version=None,
    alternative_name=None,
    text={text},
),
""".format(vartype='"%s"' % var.vartype,
          characteristics=None if var.characteristics is None else pformat(var.characteristics),
          mnemonics=nones2arg(var.mnemonics, must_be_string=True),
          requires=nones2arg(var.requires),
          excludes=nones2arg(var.excludes),
          dimensions=dimensions2arg(var.dimensions),
          varset='"%s"' % var.varset,
          abivarname='"%s"' % var.abivarname,
          commentdefault=nones2arg(var.commentdefault),
          topics=topics2arg(var.topics),
          commentdims=nones2arg(var.commentdims),
          defaultval=defaultval2arg(var.defaultval),
	  added_in_version=var.added_in_version,
	  alternative_name=var.alternative_name,
          text=text,
          )

                lines.append(s)
                #print(s)

            lines.append("]")
            # Write file
            with open(os.path.join(dirpath, "variables_%s.py" % code), "wt") as fh:
                fh.write("\n".join(lines))
                fh.write("\n")


class InputVariables(OrderedDict):
    """
    Dictionary storing the variables used by one executable.

    .. attributes:

	executable: Name of executable e.g. anaddb
    """
    @classmethod
    def from_pyfile(cls, filepath):
        """Initialize the object from python file."""
        import imp
        module = imp.load_source(filepath, filepath)
        #from importlib.machinery import SourceFileLoader
        #module = SourceFileLoader(filepath, filepath).load_module()
        vlist = [Variable(**d) for d in module.variables]
        new = cls()
        new.executable = module.executable
        for v in sorted(vlist, key=lambda v: v.name):
            new[v.name] = v
        return new

    @lazy_property
    def my_varset_list(self):
        """Set with the all the varset strings found in the database."""
        return sorted(set(v.varset for v in self.values()))

    @lazy_property
    def name2varset(self):
        """Dictionary mapping the name of the variable to the varset section."""
        d = {}
        for name, var in self.items():
            d[name] = var.varset
        return d

    @lazy_property
    def my_characteristics(self):
        """Set with all characteristics found in the database. NB [] are removed from the string."""
        allchars = []
        for var in self.values():
            if var.characteristics is not None:
                allchars.extend([c.replace("[", "").replace("]", "") for c in var.characteristics])
        return set(allchars)

    def get_all_vnames(self, with_internal=False):
        """
        Return set with all the variable names including possible aliases.
        """
        doc_vnames = []
        for name, var in self.items():
            if not with_internal and var.is_internal: continue
            doc_vnames.append(name)
            if var.alternative_name is not None:
                doc_vnames.append(var.alternative_name)
        return set(doc_vnames)

    def groupby_first_letter(self):
        """Return ordered dict mapping first_char --> list of variables."""
        keys = sorted(self.keys(), key=lambda n: n[0].upper())
        od = OrderedDict()
        for char, group in groupby(keys, key=lambda n: n[0].upper()):
            od[char] = [self[name] for name in group]
        return od

    def get_vartabs_html(self, website, page_rpath):
        """Return HTML string with all the variabes in tabular format."""
        ch2vars = self.groupby_first_letter()
        ch2vars["All"] = self.values()
        # http://getbootstrap.com/javascript/#tabs
        html = """\
<div>
<!-- Nav tabs -->
<ul class="nav nav-pills" role="tablist">\n"""

        idname = self.executable + "-tabs"
        for i, char in enumerate(ch2vars):
            id_char = "#%s-%s" % (idname, char)
            if i == 0:
                html += """\n
<li role="presentation" class="active"><a href="%s" role="tab" data-toggle="tab">%s</a></li>\n""" % (id_char, char)
            else:
                html += """\
<li role="presentation"><a href="%s" role="tab" data-toggle="tab">%s</a></li>\n""" % (id_char, char)
        html += """\
</ul>
<!-- Tab panes -->
<div class="tab-content">
        """
        for i, (char, vlist) in enumerate(ch2vars.items()):
            id_char = "%s-%s" % (idname, char)
            p = " ".join(v.internal_link(website, page_rpath, cls="small-grey-link") for v in vlist)
            if i == 0:
                html += '<div role="tabpanel" class="tab-pane active" id="%s">\n%s\n</div>\n' % (id_char, p)
            else:
                html += '<div role="tabpanel" class="tab-pane" id="%s">\n%s\n</div>\n' % (id_char, p)

        return html + "</div></div>"

    def group_by_varset(self, names):
        """
        Group a list of variable in sections.

        Args:
	    names: string or list of strings with ABINIT variable names.

        Return:
            Ordered dict mapping section_name to the list of variable names belonging to the section.
            The dict uses the same ordering as those in `self.sections`
        """
        d = defaultdict(list)

        for name in list_strings(names):
            try:
                sec = self.name2varset[name]
                d[sec].append(name)
            except KeyError as exc:
                msg = ("`%s` is not a registered variable of code `%s`.\nPerhaps you are using an old " +
                       "version of the database with a more recent Abinit?") % (name, self.executable)
                raise KeyError(msg)

        return OrderedDict([(sec, d[sec]) for sec in self.my_varset_list if d[sec]])

    def apropos(self, varname):
        """Return the list of :class:`Variable` objects that are related` to the given varname"""
        var_list = []
        for v in self.values():
            if (v.text and varname in v.text or
               (v.dimensions is not None and varname in str(v.dimensions)) or
               (v.requires is not None and varname in v.requires) or
               (v.excludes is not None and varname in v.excludes)):
                var_list.append(v)

        return var_list

    def vars_with_varset(self, sections):
        """
        List of :class:`Variable` associated to the given sections.
        sections can be a string or a list of strings.
        """
        sections = set(list_strings(sections))
        varlist = []
        for v in self.values():
            if v.varset in sections:
                varlist.append(v)

        return varlist

    def vars_with_char(self, chars):
        """
        Return list of :class:`Variable` with the specified characteristic.
        chars can be a string or a list of strings.
        """
        chars = ["[[" + c + "]]" for c in list_strings(chars)]
        varlist = []
        for v in self.values():
            if v.characteristics is None: continue
            if any(c in v.characteristics for c in chars):
                varlist.append(v)

        return varlist

    def get_graphviz_varname(self, varname, engine="automatic", graph_attr=None, node_attr=None, edge_attr=None):
        """
        Generate task graph in the DOT language (only parents and children of this task).

        Args:
            varname: Name of the variable.
            engine: ['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage']
            graph_attr: Mapping of (attribute, value) pairs for the graph.
            node_attr: Mapping of (attribute, value) pairs set for all nodes.
            edge_attr: Mapping of (attribute, value) pairs set for all edges.

        Returns: graphviz.Digraph <https://graphviz.readthedocs.io/en/stable/api.html#digraph>
        """
        var = self[varname]

        # https://www.graphviz.org/doc/info/
        from graphviz import Digraph
        graph = Digraph("task", engine="dot" if engine == "automatic" else engine)
        #graph.attr(label=repr(var))
        #graph.node_attr.update(color='lightblue2', style='filled')
        #cluster_kwargs = dict(rankdir="LR", pagedir="BL", style="rounded", bgcolor="azure2")

        # These are the default attrs for graphviz
        default_graph_attr = {
            'rankdir': 'LR',
            #'size': "8.0, 12.0",
        }
        if graph_attr is None: graph_attr = default_graph_attr

        default_node_attr = {
            #'shape': 'box',
            #'fontsize': 10,
            #'height': 0.25,
            #'fontname': '"Vera Sans, DejaVu Sans, Liberation Sans, '
            #            'Arial, Helvetica, sans"',
            #'style': '"setlinewidth(0.5)"',
        }
        if node_attr is None: node_attr = default_node_attr

        default_edge_attr = {
            #'arrowsize': '0.5',
            #'style': '"setlinewidth(0.5)"',
        }
        if edge_attr is None: edge_attr = default_edge_attr

        # Add input attributes.
        graph.graph_attr.update(**graph_attr)
        graph.node_attr.update(**node_attr)
        graph.edge_attr.update(**edge_attr)

        def node_kwargs(var):
            return dict(
                shape="box",
                fontsize="10",
                height="0.25",
                #color=var.color_hex,
                label=str(var),
                URL=var.website_url,
                target="_top",
                tooltip=str(var.mnemonics),
            )

        edge_kwargs = dict(arrowType="vee", style="solid")

        graph.node(var.name, **node_kwargs(var))
        for parent in var.get_parent_names():
            parent = self[parent]
            graph.node(parent.name, **node_kwargs(parent))
            graph.edge(parent.name, var.name, **edge_kwargs) #, label=edge_label, color=self.color_hex

        with_children = True
        if with_children: # > threshold
            # Connect task to children.
            for oname, ovar in self.items():
                if oname == varname: continue
                if varname not in ovar.get_parent_names(): continue
                graph.node(ovar.name, **node_kwargs(ovar))
                graph.edge(var.name, ovar.name, **edge_kwargs) #, label=edge_label, color=self.color_hex

        return graph

    def get_graphviz(self, varset=None, vartype=None, engine="automatic", graph_attr=None, node_attr=None, edge_attr=None):
        """
        Generate graph in the DOT language (only parents and children of this task).

        Args:
            varset: Select variables with this `varset`. Include all if None
	    vartype: Select variables with this `vartype`. Include all
            engine: ['dot', 'neato', 'twopi', 'circo', 'fdp', 'sfdp', 'patchwork', 'osage']
            graph_attr: Mapping of (attribute, value) pairs for the graph.
            node_attr: Mapping of (attribute, value) pairs set for all nodes.
            edge_attr: Mapping of (attribute, value) pairs set for all edges.

        Returns: graphviz.Digraph <https://graphviz.readthedocs.io/en/stable/api.html#digraph>
        """
        # https://www.graphviz.org/doc/info/
        from graphviz import Digraph
        graph = Digraph("task", engine="dot" if engine == "automatic" else engine)
        #graph.attr(label=repr(var))
        #graph.node_attr.update(color='lightblue2', style='filled')
        #cluster_kwargs = dict(rankdir="LR", pagedir="BL", style="rounded", bgcolor="azure2")

        # These are the default attrs for graphviz
        default_graph_attr = {
            'rankdir': 'LR',
            #'size': "8.0, 12.0",
        }
        if graph_attr is None: graph_attr = default_graph_attr

        default_node_attr = {
            #'shape': 'box',
            #'fontsize': 10,
            #'height': 0.25,
            #'fontname': '"Vera Sans, DejaVu Sans, Liberation Sans, '
            #            'Arial, Helvetica, sans"',
            #'style': '"setlinewidth(0.5)"',
        }
        if node_attr is None: node_attr = default_node_attr

        default_edge_attr = {
            #'arrowsize': '0.5',
            #'style': '"setlinewidth(0.5)"',
        }
        if edge_attr is None: edge_attr = default_edge_attr

        # Add input attributes.
        graph.graph_attr.update(**graph_attr)
        graph.node_attr.update(**node_attr)
        graph.edge_attr.update(**edge_attr)

        def node_kwargs(var):
            return dict(
                shape="box",
                fontsize="10",
                height="0.25",
                #color=var.color_hex,
                label=str(var),
                URL=var.website_url,
                target="_top",
                tooltip=str(var.mnemonics),
            )

        edge_kwargs = dict(arrowType="vee", style="solid")
        with_children = False

        for name, var in self.items():
            if vartype is not None and var.vartype != vartype: continue
            if varset is not None and var.varset != varset: continue

            graph.node(var.name, **node_kwargs(var))
            for parent in var.get_parent_names():
                parent = self[parent]
                graph.node(parent.name, **node_kwargs(parent))
                graph.edge(parent.name, var.name, **edge_kwargs) #, label=edge_label, color=self.color_hex

            if with_children: # > threshold
                # Connect task to children.
                for oname, ovar in self.items():
                    if oname == varname: continue
                    if varname not in ovar.get_parent_names(): continue
                    graph.node(ovar.name, **node_kwargs(ovar))
                    graph.edge(var.name, ovar.name, **edge_kwargs) #, label=edge_label, color=self.color_hex

        return graph
