import json
import os
import sys
from collections import OrderedDict, defaultdict
from itertools import groupby

# Helper functions (coming from AbiPy)


class lazy_property:
    """
    Descriptor for implementing lazy attributes.

    Lazy attributes are evaluated only on the first access, and the result
    is cached in the object's `__dict__`.
    """

    def __init__(self, func):
        self.__func = func
        from functools import wraps

        wraps(self.__func)(self)

    def __get__(self, inst, inst_cls):
        if inst is None:
            return self

        if not hasattr(inst, "__dict__"):
            raise AttributeError("'%s' object has no attribute '__dict__'" % (inst_cls.__name__,))

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
            raise AttributeError("'%s' object has no attribute '__dict__'" % (inst_cls.__name__,))

        if name.startswith("__") and not name.endswith("__"):
            name = "_%s%s" % (inst_cls.__name__, name)

        if not isinstance(getattr(inst_cls, name), cls):
            raise AttributeError(
                "'%s.%s' is not a %s attribute" % (inst_cls.__name__, name, cls.__name__)
            )

        if name in inst.__dict__:
            del inst.__dict__[name]


def is_string(s):
    """
    Determine if an object behaves like a string.

    Args:
        s: The object to test.

    Returns:
        bool: True if it behaves like a string, False otherwise.
    """
    try:
        s + " "
        return True
    except TypeError:
        return False


def list_strings(arg):
    """
    Ensure the output is a list of strings.

    If the input is already a list of strings, it is returned as is.
    If the input is a single string, it is wrapped in a list.

    Args:
        arg: A string or a list of strings.

    Returns:
        list: A list of strings.

    Examples:
        >>> list_strings("A single string")
        ['A single string']
        >>> list_strings(["A", "list"])
        ['A', 'list']
    """
    if is_string(arg):
        return [arg]
    return arg


def splitall(path):
    """
    Split a file path into all its component parts.

    Args:
        path: The path string to split.

    Returns:
        list: List of path components.
    """
    allparts = []
    while True:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        if parts[1] == path:  # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        path = parts[0]
        allparts.insert(0, parts[1])
    return allparts


# Unit names supported in Abinit input.
ABI_UNITS = [
    "au",
    "Angstr",
    "Angstrom",
    "Angstroms",
    "Bohr",
    "Bohrs",
    "eV",
    "Ha",
    "Hartree",
    "Hartrees",
    "K",
    "Ry",
    "Rydberg",
    "Rydbergs",
    "T",
    "Tesla",
    "Second",
    "S",
    "Sec",
]

# Operators supported by parser
ABI_OPS = ["sqrt", "end", "*", "/"]


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
ABI_EXTERNAL_PARAMS = OrderedDict(
    [
        ("AUTO_FROM_PSP", "Means that the value is read from the PSP file"),
        ("CUDA", "True if ABINIT has been compiled using Nvidia CUDA (compilation for Nvidia GPU)"),
        ("ETSF_IO", "True if NetCDF is enabled (compilation)"),
        ("FFTW3", "True if FFTW3 is enabled (compilation)"),
        (
            "GPU",
            "True if ABINIT has been compiled using one of the GPU implementations (CUDA, OPENMP_OFFLOAD, KOKKOS)",
        ),
        (
            "KOKKOS",
            "True if ABINIT has been compiled using KOKKOS performance library (compilation for GPU accelerators)",
        ),
        ("MPI_IO", "True if MPI_IO is enabled (compilation)"),
        ("NPROC", "Number of processors used for Abinit"),
        ("NVTX", "True if ABINIT has been linked to the NVIDIA® Tools Extension SDK (NVTX)"),
        (
            "OPENMP",
            "True if ABINIT has been compiled using OPENMP multithreading (compilation for multicore processors)",
        ),
        (
            "OPENMP_OFFLOAD",
            "True if ABINIT has been compiled using OPENMP_OFFLOAD (openMP v5+) (compilation for GPU accelerators)",
        ),
        ("PARALLEL", "True if the code is compiled with MPI"),
        ("ROCTX", "True if ABINIT has been linked to the AMD ROCm Tools Extension SDK (ROCTX)"),
        ("SEQUENTIAL", "True if the code is compiled without MPI"),
    ]
)

# List of topics
# The topics should be declared both in this file and in mkdocs.yml.in
ABI_TOPICS = [
    "Abipy",
    "APPA",
    "Artificial",
    "aTDEP",
    "AtomCentered",
    "AtomManipulator",
    "AtomTypes",
    "Bader",
    "Band2eps",
    "Berry",
    "BandOcc",
    "BoundingProcess",
    "BSE",
    "ConstrainedDFPT",
    "ConstrainedDFT",
    "ConstrainedPol",
    "Control",
    "Coulomb",
    "CrossingBarriers",
    "CalcUJ",
    "crystal",
    "DFT+U",
    "DeltaSCF",
    "DensityPotential",
    "Dev",
    "DFPT",
    "DMFT",
    "DmftTriqsCthyb",
    "DynamicsMultibinit",
    "EffectiveMass",
    "EFG",
    "Elastic",
    "ElPhonInt",
    "ElPhonTransport",
    "ElecDOS",
    "ElecBandStructure",
    "ExtFPMD",
    "FileFormats",
    "FitProcess",
    "ForcesStresses",
    "FrequencyMeshMBPT",
    "GeoConstraints",
    "GeoOpt",
    "Git",
    "GSintroduction",
    "GW",
    "GWR",
    "GWls",
    "Hybrids",
    "Install",
    "k-points",
    "LatticeModel",
    "LatticeWannier",
    "LWFModel",
    "LDAminushalf",
    "longwave",
    "LOTF",
    "MagField",
    "MagMom",
    "MolecularDynamics",
    "Macroave",
    "multidtset",
    "NMR",
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
    "Polaron",
    "PortabilityNonRegression",
    "positron",
    "printing",
    "PseudosPAW",
    "q-points",
    "RandStopPow",
    "Recursion",
    "RPACorrEn",
    "RTTDDFT",
    "SCFControl",
    "SCFAlgorithms",
    "SelfEnergy",
    "SmartSymm",
    "spinpolarisation",
    "SpinDynamicsMultibinit",
    "STM",
    "Susceptibility",
    "TDDFT",
    "TDepES",
    "Temperature",
    "TransPath",
    "TuningSpeedMem",
    "Unfolding",
    "UnitCell",
    "vdw",
    "Verification",
    "Wannier",
    "Wavelets",
    "xc",
]

# Relevance associated to the topic
ABI_RELEVANCES = OrderedDict(
    [
        ("compulsory", "Compulsory input variables"),
        ("basic", "Basic input variables"),
        ("useful", "Useful input variables"),
        ("internal", "Relevant internal variables"),
        ("prpot", "Printing input variables for potentials"),
        ("prfermi", "Printing input variables for fermi level or surfaces"),
        (
            "prden",
            "Printing input variables for density, eigenenergies, k-points and wavefunctions",
        ),
        ("prgeo", "Printing input variables for geometry"),
        ("prdos", "Printing DOS-related input variables"),
        ("prgs", "Printing other ground-state input variables"),
        ("prngs", "Printing non-ground-state input variables"),
        ("prmisc", "Printing miscellaneous files"),
        ("expert", "Input variables for experts"),
    ]
)


class Variable:
    """
    Gather information about an ABINIT input variable.

    Attributes:
        abivarname: Fully qualified variable name (e.g., asr@anaddb).
        varset: The variable set/group name.
        vartype: The type of the variable (integer, real, string).
        topics: List of associated topics and relevances.
        dimensions: Dimensions description or "scalar".
        defaultval: Default value of the variable.
        mnemonics: Short mnemonic description.
        characteristics: List of variable characteristics (e.g., ENERGY).
        excludes: List of variables that cannot be used with this one.
        requires: List of variables required by this one.
        commentdefault: Optional comment about the default value.
        commentdims: Optional comment about the dimensions.
        added_in_version: ABINIT version when the variable was introduced.
        alternative_name: Alias or old name.
        text: Markdown string containing the main documentation.
    """

    def __init__(
        self,
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
        Initialize a Variable object.

        Args:
            abivarname: Variable name, including @code suffix if applicable.
            varset: The group name defining the variable's category.
            vartype: String specifying the data type.
            topics: List of topic strings in 'TopicName_Relevance' format.
            dimensions: List of strings for dimensions or "scalar".
            defaultval: Default value or formula.
            mnemonics: Brief mnemonic description.
            characteristics: List of property flags.
            excludes: Comma-separated variable names that conflict with this one.
            requires: Comma-separated variable names required by this one.
            commentdefault: Documentation comment for the default value.
            commentdims: Documentation comment for the dimensions.
            added_in_version: ABINIT version identifier.
            alternative_name: Alias name for backward compatibility.
            text: Markdown documentation content.

        Raises:
            ValueError: If mandatory attributes are missing.
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
        for a in (
            "abivarname",
            "varset",
            "vartype",
            "topics",
            "dimensions",
            "added_in_version",
            "text",
        ):
            if getattr(self, a) is None:
                errors.append("attribute %s is mandatory" % a)
        if errors:
            raise ValueError("Errors in %s:\n%s" % (self.abivarname, "\n".join(errors)))

    @lazy_property
    def name(self):
        """
        Generate the normalized variable name.

        Returns:
            str: Lowercase variable name without the executable suffix.
        """
        return (
            self.abivarname.lower()
            if "@" not in self.abivarname
            else self.abivarname.split("@")[0].lower()
        )

    @lazy_property
    def executable(self):
        """
        Identify the executable associated with the variable.

        Returns:
            str: Code name (e.g., 'abinit', 'anaddb').
        """
        if "@" in self.abivarname:
            code = self.abivarname.split("@")[1]
            assert code == self.varset
        else:
            code = "abinit"
        return code

    @lazy_property
    def website_url(self):
        """
        Construct the documentation URL for the variable.

        Returns:
            str: Absolute URL to the official documentation page.
        """
        # This is gonna be the official API on the server
        # docs.abinit.org/vardocs/CODENAME/VARNAME?version=8.6.2
        # return "https://docs.abinit.org/vardocs/%s/%s" % (self.executable, self.name)

        # For the time being, we have to use:
        # variables/eph/#asr
        # variables/anaddb#asr
        if self.executable == "abinit":
            return "https://docs.abinit.org/variables/%s#%s" % (self.varset, self.name)
        return "https://docs.abinit.org/variables/%s#%s" % (self.executable, self.name)

    @lazy_property
    def topic2relevances(self):
        """
        Map topics to their associated list of relevances.

        Returns:
            OrderedDict: Mapping of topic names to lists of relevance strings.
        """
        assert self.topics is not None
        od = OrderedDict()
        for tok in self.topics:
            topic, relevance = [s.strip() for s in tok.split("_")]
            if topic not in od:
                od[topic] = []
            od[topic].append(relevance)
        return od

    @lazy_property
    def is_internal(self):
        """
        Check if the variable is for internal use only.

        Returns:
            bool: True if identified as INTERNAL_ONLY.
        """
        return self.characteristics is not None and "[[INTERNAL_ONLY]]" in self.characteristics

    @lazy_property
    def wikilink(self):
        """
        Generate the ABINIT wikilink syntax for the variable.

        Returns:
            str: Wikilink string.
        """
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
        if other is None:
            return False
        return self.abivarname == other.abivarname

    def __ne__(self, other):
        return not (self == other)

    @lazy_property
    def info(self):
        """
        Produce a JSON-formatted string containing variable metadata.

        Returns:
            str: JSON representation of selected attributes.
        """
        attrs = [
            "vartype",
            "characteristics",
            "mnemonics",
            "dimensions",
            "defaultval",
            "abivarname",
            "commentdefault",
            "commentdims",
            "varset",
            "requires",
            "excludes",
            "added_in_version",
            "alternative_name",
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
            html = (
                "<h2>Default value:</h2>"
                + my_unicode(self.defaultval)
                + "<br/><h2>Description</h2>"
                + self.text
            )
            return html.replace("[[", "<b>").replace("]]", "</b>")
        md = self.text.replace("[[", "<b>").replace("]]", "</b>")
        return markdown.markdown(f"""
## Default value:
{my_unicode(self.defaultval)}

## Description:
{my_unicode(md)}
""")

    def browse(self):
        """
        Open the variable's documentation in the default web browser.

        Returns:
            bool: True if the browser was successfully opened.
        """
        import webbrowser

        return webbrowser.open(self.website_url)

    @lazy_property
    def isarray(self):
        """
        Determine if the variable is an array.

        Returns:
            bool: True if dimensions are not 'scalar'.
        """
        return not (is_string(self.dimensions) and self.dimensions == "scalar")

    def depends_on_dimension(self, dimname):
        """
        Check if the variable's shape depends on a specific dimension.

        Args:
            dimname: Name of the dimension variable or the Variable object itself.

        Returns:
            bool: True if the dimension name is found in the dimensions list.
        """
        if not self.isarray:
            return False
        if isinstance(dimname, Variable):
            dimname = dimname.name
        # This test is not very robust and can fail.
        # Assume no space between `[` and name (there should be a test for this...)
        key = "[[%s]]" % dimname
        for d in self.dimensions:
            if key in str(d):
                return True
        return False

    def html_link(self, label=None):
        """
        Generate an HTML anchor tag for the variable's documentation.

        Args:
            label: Optional link text. Defaults to the variable name.

        Returns:
            str: HTML anchor tag.
        """
        label = self.name if label is None else label
        return '<a href="%s" target="_blank">%s</a>' % (self.website_url, label)

    def get_parent_names(self):
        """
        Identify variables that directly influence this one.

        This includes variables mentioned in the dimensions or the 'requires' list.

        Returns:
            set: Set of parent variable names.
        """
        # if hasattr(self, ...
        import re

        parent_names = []
        WIKILINK_RE = r"\[\[([\w0-9_ -]+)\]\]"
        # TODO
        #  parent = self[parent]
        # KeyError: "'nzchempot'
        # WIKILINK_RE = r'\[\[([^\[]+)\]\]'
        if isinstance(self.dimensions, (list, tuple)):
            for dim in self.dimensions:
                dim = str(dim)
                m = re.match(WIKILINK_RE, dim)
                if m:
                    parent_names.append(m.group(1))

        if self.requires is not None:
            parent_names.extend([m.group(1) for m in re.finditer(WIKILINK_RE, self.requires) if m])

        # Convert to set and remove possible self-reference.
        parent_names = set(parent_names)
        parent_names.discard(self.name)
        return parent_names

    def internal_link(self, website, page_rpath, label=None, cls=None):
        """
        Generate an HTML link for use within the ABINIT website.

        Args:
            website: The Website instance.
            page_rpath: The root-relative path of the current page.
            label: Optional link text.
            cls: Optional CSS class for the anchor tag.

        Returns:
            str: HTML anchor tag with relative URL.
        """
        token = "%s:%s" % (self.executable, self.name)
        a = website.get_wikilink(token, page_rpath, is_html=True)
        cls = a.get("class") if cls is None else cls
        return '<a href="%s" class="%s">%s</a>' % (
            a.get("href"),
            cls,
            a.text if label is None else label,
        )

    @staticmethod
    def format_dimensions(dimensions):
        """Pretty print dimensions."""
        if dimensions is None:
            s = ""
        elif dimensions == "scalar":
            s = "scalar"
        # s = str(dimensions)
        elif isinstance(dimensions, (list, tuple)):
            s = "("
            for dim in dimensions:
                s += str(dim) + ","
            s = s[:-1]
            s += ")"
        else:
            s = str(dimensions)

        return s

    def to_abimarkdown(self, with_hr=True):
        """
        Convert the variable metadata and description to a Markdown string.

        Supports ABINIT-specific Markdown extensions.

        Args:
            with_hr: If True, append a horizontal rule at the end.

        Returns:
            str: Formatted Markdown string.
        """
        lines = []
        app = lines.append

        app("## **%s** \n\n" % self.name)
        app("*Mnemonics:* %s  " % str(self.mnemonics))
        if self.characteristics:
            app("*Characteristics:* %s  " % ", ".join(self.characteristics))
        if self.topic2relevances:
            app(
                "*Mentioned in topic(s):* %s  "
                % ", ".join("[[topic:%s]]" % k for k in self.topic2relevances)
            )
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
                frequency,
                len(self.tests),
                tests_info["num_all_tests"],
                self.executable,
                tests_info["num_tests_in_tutorial"],
                tests_info["num_all_tutorial_tests"],
                self.executable,
            )

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
                if count > max_ntests:
                    ipaths = ipaths[: min(3, len(ipaths))]
                s = (
                    "- "
                    + suite_name
                    + ":  "
                    + ", ".join("[[%s|%s]]" % (p, os.path.basename(p)) for p in ipaths)
                )
                if count > max_ntests:
                    s += " ..."
                app("    " + s)
            app("\n\n")

        # Add text with description.
        app(2 * "\n")
        # Replace all occurrences of [[name]] with **name** to reduce number of html links in docs
        new_text = self.text.replace("[[%s]]" % self.name, " **%s** " % self.name)
        app(new_text)
        if with_hr:
            app("* * *" + 2 * "\n")

        return "\n".join(lines)

    def validate(self):
        """
        Validate the variable's attributes for consistency and correctness.

        Checks mandatory fields, allowed topics, relevances, and formatting.

        Raises:
            ValueError: If validation errors are found.
        """
        errors = []
        eapp = errors.append

        try:
            svar = str(self)
        except Exception as exc:
            svar = "Unknown"
            eapp(str(exc))

        if self.abivarname is None:
            eapp("Variable `%s` has no name" % svar)

        allowed_vartypes = {"integer", "real", "string", "integer or string"}

        if self.vartype is None:
            eapp("Variable `%s` has no vartype" % svar)
        elif self.vartype not in allowed_vartypes:
            eapp(f"{svar} must be in {allowed_vartypes} while it is {self.vartype}")

        if self.topics is None:
            eapp("%s does not have at least one topic and the associated relevance" % svar)

        for topic, relevances in self.topic2relevances.items():
            if topic not in ABI_TOPICS:
                eapp(
                    "%s delivers topic `%s` that does not belong to the allowed list"
                    % (sname, topic)
                )
            for relevance in relevances:
                if relevance not in ABI_RELEVANCES:
                    eapp(
                        "%s delivers relevance `%s` that does not belong to the allowed list"
                        % (sname, relevance)
                    )

        # Compare the characteristics of this variable with the refs to detect possible typos.
        if self.characteristics is not None:
            if not isinstance(self.characteristics, list):
                eapp("The field characteristics of %s is not a list" % svar)
            else:
                for cat in self.characteristics:
                    if cat.replace("[[", "").replace("]]", "") not in ABI_CHARACTERISTICS:
                        eapp("The characteristics %s of %s is not valid" % (cat, svar))

        if self.dimensions is None:
            eapp(
                "%s does not have a dimension. If it is a *scalar*, it must be declared so." % svar
            )
        #elif self.dimensions != "scalar":
        #    if not isinstance(self.dimensions, (list, ValueWithConditions)):
        #        eapp(
        #            "The dimensions field of %s is not a list neither a valuewithconditions" % svar
        #        )

        if self.varset is None:
            eapp("`%s` does not have a varset" % svar)
        # else:
        #    if not isinstance(self.varset, str) or self.varset not in ref_varset:
        #        print('The field varset of %s should be one of the valid varsets' % str(self))

        #if len(self.name) > 25:
        #    eapp("Length of `%s` is longer than 25 characters." % self.name)

        if errors:
            raise ValueError("\n".join(errors))


class ValueWithUnit:
    """
    Representation of a value associated with a specific unit.

    Attributes:
        value: The numerical or formula value.
        units: String specifying the unit (e.g., 'Ha', 'Bohr').
    """

    def __init__(self, value=None, units=None):
        self.value = value
        self.units = units

    def __str__(self):
        return str(self.value) + " " + str(self.units)

    def __repr__(self):
        return str(self)


class Range:
    """
    Representation of a numerical range [start, stop].

    Attributes:
        start: The lower bound of the range.
        stop: The upper bound of the range.
    """

    start = None
    stop = None

    def __init__(self, start=None, stop=None):
        """
        Initialize a Range.

        Args:
            start: Lower bound.
            stop: Upper bound.
        """
        self.start = start
        self.stop = stop

    def isin(self, value):
        """
        Check if a value falls within the range.

        Args:
            value: The numerical value to check.

        Returns:
            bool: True if inside the range, False otherwise.
        """
        isin = True
        if self.start is not None:
            isin = isin and (self.start <= self.value)
        if self.stop is not None:
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
        return None


class ValueWithConditions(dict):
    """
    Dictionary-like object for values that depend on specific conditions.

    Example:
        `ValueWithConditions({'[[paral_kgb]]==1': '6', 'defaultval': 2})`
        means the value is 6 if `paral_kgb` is 1, otherwise it is 2.
    """

    def __repr__(self):
        s = ""
        for key in self:
            if key != "defaultval":
                s += str(self[key]) + " if " + str(key) + ",\n"
        s += str(self["defaultval"]) + " otherwise.\n"
        return s

    def __str__(self):
        return self.__repr__()


class MultipleValue:
    """
    Representation of repeating values (equivalent to `X * Y` in ABINIT input).

    Attributes:
        number: Number of repetitions.
        value: The value to repeat.
    """

    def __init__(self, number=None, value=None):
        self.number = number
        self.value = value

    def __repr__(self):
        if self.number is None:
            return "*" + str(self.value)
        return str(self.number) + " * " + str(self.value)


def my_unicode(s):
    """
    Convert a string or object to a unicode string.

    Args:
        s: The object to convert.

    Returns:
        str: The unicode string representation.
    """
    return unicode(s) if sys.version_info[0] <= 2 else str(s)


##############
# Public API #
##############


_VARS = None


def get_codevars():
    """
    Get the global variable database, initializing and caching it if necessary.

    Returns:
        VarDatabase: The cached variable database instance.
    """
    global _VARS
    if _VARS is None:
        _VARS = VarDatabase.from_pyfiles()
    return _VARS


class VarDatabase(OrderedDict):
    """
    Registry for all input variables across all ABINIT executables.

    Items are indexed by executable name (e.g., 'abinit', 'anaddb').
    """

    all_characteristics = ABI_CHARACTERISTICS
    all_external_params = ABI_EXTERNAL_PARAMS

    @classmethod
    def from_pyfiles(cls, dirpath=None):
        """
        Initialize the database by scanning Python modules in a directory.

        Modules must name `variables_CODE.py`.

        Args:
            dirpath: Directory path to scan. Defaults to the directory of this module.

        Returns:
            VarDatabase: Loaded database instance.
        """
        if dirpath is None:
            dirpath = os.path.dirname(os.path.abspath(__file__))
        pyfiles = [
            os.path.join(dirpath, f)
            for f in os.listdir(dirpath)
            if f.startswith("variables_") and f.endswith(".py")
        ]
        new = cls()
        for pyf in pyfiles:
            vd = InputVariables.from_pyfile(pyf)
            new[vd.executable] = vd

        return new

    def iter_allvars(self):
        """
        Iterate over every variable in the database across all executables.

        Yields:
            Variable: Each variable instance in the database.
        """
        for vd in self.values():
            for var in vd.values():
                yield var

    def get_version_endpoints(self):
        """
        Generate mapping of variable names to their relative URLs for the web server.

        Returns:
            tuple: (version_string, code_urls_dict)
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
        Update a JSON file mapping variable names to relative URLs.

        Args:
            json_path: Path to the JSON file to update.
            indent: JSON indentation level.

        Raises:
            AssertionError: If the version already exists in the file.
        """
        with open(json_path) as fh:
            oldd = json.load(fh)

        new_version, newd = self.get_version_endpoints()
        assert new_version not in oldd
        oldd[new_version] = newd
        with open(json_path, "w") as fh:
            json.dump(oldd, fh, indent=indent)

    def _write_pymods(self, dirpath="."):
        """
        Regenerate the Python modules defining the variables.

        Internal method for database maintenance.

        Args:
            dirpath: Target directory for the generated .py files.
        """
        dirpath = os.path.abspath(dirpath)
        from pprint import pformat

        def nones2arg(obj, must_be_string=False):
            if obj is None:
                if must_be_string:
                    raise TypeError("obj must be string.")
                return None
            if isinstance(obj, str):
                s = str(obj).rstrip()
                if "\n" in s:
                    return '"""%s"""' % s
                if "'" in s:
                    return '"%s"' % s
                if '"' in s:
                    return "'%s'" % s
                return '"%s"' % s
            raise TypeError("%s: %s" % (type(obj), str(obj)))

        def topics2arg(obj):
            if isinstance(obj, str):
                if "," in obj:
                    obj = [s.strip() for s in obj.split(",")]
                else:
                    obj = [obj]
            if isinstance(obj, (list, tuple)):
                return pformat(obj)

            raise TypeError("%s: %s" % (type(obj), str(obj)))

        def dimensions2arg(obj):
            if isinstance(obj, str) and obj == "scalar":
                return '"scalar"'
            if isinstance(obj, (ValueWithUnit, MultipleValue, Range, ValueWithConditions)):
                return "%s(%s)" % (obj.__class__.__name__, pformat(obj.__dict__))
            if isinstance(obj, (list, tuple)):
                return pformat(obj)

            raise TypeError("%s, %s" % (type(obj), str(obj)))

        def defaultval2arg(obj):
            if obj is None:
                return obj
            if isinstance(obj, (ValueWithUnit, MultipleValue, Range, ValueWithConditions)):
                return "%s(%s)" % (obj.__class__.__name__, pformat(obj.__dict__))
            if isinstance(obj, (list, tuple)):
                return pformat(obj)
            if isinstance(obj, str):
                return '"%s"' % str(obj)
            if isinstance(obj, (int, float)):
                return obj

            raise TypeError("%s, %s" % (type(obj), str(obj)))

        for code in self:
            varsd = self[code]

            lines = [
                """\
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
""".format(
                    vartype='"%s"' % var.vartype,
                    characteristics=None
                    if var.characteristics is None
                    else pformat(var.characteristics),
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
                    text=text,
                )

                lines.append(s)
                # print(s)

            lines.append("]")
            # Write file
            with open(os.path.join(dirpath, "variables_%s.py" % code), "w") as fh:
                fh.write("\n".join(lines))
                fh.write("\n")


class InputVariables(OrderedDict):
    """
    Collection of input variables for a specific executable.

    Attributes:
        executable: Name of the associated executable (e.g., 'anaddb').
    """

    @classmethod
    def from_pyfile(cls, filepath):
        """
        Load variables for an executable from a formatted Python module.

        Args:
            filepath: Path to the Python file.

        Returns:
            InputVariables: Loaded instance.
        """
        try:
            import imp

            module = imp.load_source(filepath, filepath)
        except ModuleNotFoundError:
            from importlib.machinery import SourceFileLoader

            module = SourceFileLoader(filepath, filepath).load_module()
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
        Retrieve all variable names associated with the executable.

        Args:
            with_internal: If True, include internal variables.

        Returns:
            set: Set of variable names and their alternative names.
        """
        doc_vnames = []
        for name, var in self.items():
            if not with_internal and var.is_internal:
                continue
            doc_vnames.append(name)
            if var.alternative_name is not None:
                doc_vnames.append(var.alternative_name)
        return set(doc_vnames)

    def groupby_first_letter(self):
        """
        Group variables by their first letter.

        Returns:
            OrderedDict: Mapping of uppercase letters to lists of Variable objects.
        """
        keys = sorted(self.keys(), key=lambda n: n[0].upper())
        od = OrderedDict()
        for char, group in groupby(keys, key=lambda n: n[0].upper()):
            od[char] = [self[name] for name in group]
        return od

    def group_by_varset(self, names):
        """
        Group a list of variable names by their associated varset identifier.

        Args:
            names: A single variable name string or a list of names.

        Returns:
            OrderedDict: Mapping of varset strings to lists of variable names.

        Raises:
            KeyError: If a variable name is not found in the database.
        """
        d = defaultdict(list)

        for name in list_strings(names):
            try:
                sec = self.name2varset[name]
                d[sec].append(name)
            except KeyError:
                msg = (
                    "`%s` is not a registered variable of code `%s`.\nPerhaps you are using an old "
                    "version of the database with a more recent Abinit?"
                ) % (name, self.executable)
                raise KeyError(msg)

        return OrderedDict([(sec, d[sec]) for sec in self.my_varset_list if d[sec]])

    def apropos(self, varname):
        """
        Find all variables that refer to or are related to a specific name.

        Args:
            varname: The substring or name to search for across text, dimensions,
                and dependencies.

        Returns:
            list: List of matching Variable objects.
        """
        var_list = []
        for v in self.values():
            if (
                (v.text and varname in v.text)
                or (v.dimensions is not None and varname in str(v.dimensions))
                or (v.requires is not None and varname in v.requires)
                or (v.excludes is not None and varname in v.excludes)
            ):
                var_list.append(v)

        return var_list

    def vars_with_varset(self, sections):
        """
        Filter variables that belong to specific varset sections.

        Args:
            sections: A single section name or a list of section names.

        Returns:
            list: List of matching Variable objects.
        """
        sections = set(list_strings(sections))
        varlist = []
        for v in self.values():
            if v.varset in sections:
                varlist.append(v)

        return varlist

    def vars_with_char(self, chars):
        """
        Filter variables by their associated characteristics.

        Args:
            chars: A single characteristic name or a list of names.

        Returns:
            list: List of matching Variable objects.
        """
        chars = ["[[" + c + "]]" for c in list_strings(chars)]
        varlist = []
        for v in self.values():
            if v.characteristics is None:
                continue
            if any(c in v.characteristics for c in chars):
                varlist.append(v)

        return varlist

    def get_graphviz_varname(
        self, varname, engine="automatic", graph_attr=None, node_attr=None, edge_attr=None
    ):
        """
        Generate a Graphviz Digraph showing the dependencies of a specific variable.

        Includes both parents (variables it depends on) and children (variables
        that depend on it).

        Args:
            varname: Name of the target variable.
            engine: Graphviz layout engine.
            graph_attr: Custom graph attributes.
            node_attr: Custom default node attributes.
            edge_attr: Custom default edge attributes.

        Returns:
            graphviz.Digraph: The generated dependency graph.
        """
        var = self[varname]

        # https://www.graphviz.org/doc/info/
        from graphviz import Digraph

        graph = Digraph("task", engine="dot" if engine == "automatic" else engine)
        # graph.attr(label=repr(var))
        # graph.node_attr.update(color='lightblue2', style='filled')
        # cluster_kwargs = dict(rankdir="LR", pagedir="BL", style="rounded", bgcolor="azure2")

        # These are the default attrs for graphviz
        default_graph_attr = {
            "rankdir": "LR",
            # 'size': "8.0, 12.0",
        }
        if graph_attr is None:
            graph_attr = default_graph_attr

        default_node_attr = {
            # 'shape': 'box',
            # 'fontsize': 10,
            # 'height': 0.25,
            # 'fontname': '"Vera Sans, DejaVu Sans, Liberation Sans, '
            #            'Arial, Helvetica, sans"',
            # 'style': '"setlinewidth(0.5)"',
        }
        if node_attr is None:
            node_attr = default_node_attr

        default_edge_attr = {
            # 'arrowsize': '0.5',
            # 'style': '"setlinewidth(0.5)"',
        }
        if edge_attr is None:
            edge_attr = default_edge_attr

        # Add input attributes.
        graph.graph_attr.update(**graph_attr)
        graph.node_attr.update(**node_attr)
        graph.edge_attr.update(**edge_attr)

        def node_kwargs(var):
            return dict(
                shape="box",
                fontsize="10",
                height="0.25",
                # color=var.color_hex,
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
            # , label=edge_label, color=self.color_hex
            graph.edge(parent.name, var.name, **edge_kwargs)

        with_children = True
        if with_children:  # > threshold
            # Connect task to children.
            for oname, ovar in self.items():
                if oname == varname:
                    continue
                if varname not in ovar.get_parent_names():
                    continue
                graph.node(ovar.name, **node_kwargs(ovar))
                # , label=edge_label, color=self.color_hex
                graph.edge(var.name, ovar.name, **edge_kwargs)

        return graph

    def get_graphviz(
        self,
        varset=None,
        vartype=None,
        engine="automatic",
        graph_attr=None,
        node_attr=None,
        edge_attr=None,
    ):
        """
        Generate a Graphviz Digraph for multiple variables based on filtering.

        Args:
            varset: Filter variables by this varset.
            vartype: Filter variables by this type.
            engine: Graphviz layout engine.
            graph_attr: Custom graph attributes.
            node_attr: Custom default node attributes.
            edge_attr: Custom default edge attributes.

        Returns:
            graphviz.Digraph: The generated dependency graph.
        """
        # https://www.graphviz.org/doc/info/
        from graphviz import Digraph

        graph = Digraph("task", engine="dot" if engine == "automatic" else engine)
        # graph.attr(label=repr(var))
        # graph.node_attr.update(color='lightblue2', style='filled')
        # cluster_kwargs = dict(rankdir="LR", pagedir="BL", style="rounded", bgcolor="azure2")

        # These are the default attrs for graphviz
        default_graph_attr = {
            "rankdir": "LR",
            # 'size': "8.0, 12.0",
        }
        if graph_attr is None:
            graph_attr = default_graph_attr

        default_node_attr = {
            # 'shape': 'box',
            # 'fontsize': 10,
            # 'height': 0.25,
            # 'fontname': '"Vera Sans, DejaVu Sans, Liberation Sans, '
            #            'Arial, Helvetica, sans"',
            # 'style': '"setlinewidth(0.5)"',
        }
        if node_attr is None:
            node_attr = default_node_attr

        default_edge_attr = {
            # 'arrowsize': '0.5',
            # 'style': '"setlinewidth(0.5)"',
        }
        if edge_attr is None:
            edge_attr = default_edge_attr

        # Add input attributes.
        graph.graph_attr.update(**graph_attr)
        graph.node_attr.update(**node_attr)
        graph.edge_attr.update(**edge_attr)

        def node_kwargs(var):
            return dict(
                shape="box",
                fontsize="10",
                height="0.25",
                # color=var.color_hex,
                label=str(var),
                URL=var.website_url,
                target="_top",
                tooltip=str(var.mnemonics),
            )

        edge_kwargs = dict(arrowType="vee", style="solid")
        with_children = False

        for name, var in self.items():
            if vartype is not None and var.vartype != vartype:
                continue
            if varset is not None and var.varset != varset:
                continue

            graph.node(var.name, **node_kwargs(var))
            for parent in var.get_parent_names():
                parent = self[parent]
                graph.node(parent.name, **node_kwargs(parent))
                # , label=edge_label, color=self.color_hex
                graph.edge(parent.name, var.name, **edge_kwargs)

            if with_children:  # > threshold
                # Connect task to children.
                for oname, ovar in self.items():
                    if oname == name:
                        continue
                    if name not in ovar.get_parent_names():
                        continue
                    graph.node(ovar.name, **node_kwargs(ovar))
                    # , label=edge_label, color=self.color_hex
                    graph.edge(var.name, ovar.name, **edge_kwargs)

        return graph
