# coding: utf-8
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import io

from html2text import html2text

try:
    import yaml
except ImportError:
    raise ImportError("pyyaml package is not installed. Install it with `pip install pyyaml`")

from itertools import groupby


def format_dimensions(dimensions):

  if dimensions is None:
    s = ''
  elif dimensions == "scalar":
    s = 'scalar'
  else:
    #s = str(dimensions)
    if isinstance(dimensions,list):
      s = '('
      for dim in dimensions:
        s += str(dim) + ','

      s = s[:-1]
      s += ')'
    else:
      s = str(dimensions)

  return s

class literal(str): pass

def literal_unicode_representer(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')

yaml.add_representer(literal, literal_unicode_representer)

def valuewithunit_representer(dumper, data):
    return dumper.represent_mapping('!valuewithunit', data.__dict__)

####################################################################################################

# Classes

####################################################################################################

class Variable(yaml.YAMLObject):
    vartype = ''                   # String containing the type
    characteristics = None          # String containing the characteristics
    mnemonics = None               # String containing the mnemonics
    dimensions = None              # Array containing either int, formula or another variable  TODO: shape?
    defaultval = None              # Either constant number, formula or another variable
    text = None                    # Description (str)
    abivarname = None              # Name of the variable (str)   code@name
    commentdefault = None
    commentdims = None
    varset = None
    range = None        # TODO This is not used and should be removed
    requires = None
    excludes = None
    executables = None  # TODO This is not used and should be removed
    topics = None

    yaml_tag = u'!variable'

    def attrs(self):
        return ['vartype', 'characteristics', 'mnemonics', 'dimensions', 'defaultval', 'text',
                'abivarname', 'varset', 'executables', 'topics']

    def __init__(self, vartype=None, characteristics=None,
                 mnemonics=None, dimensions=None, defaultval=None,
                 text=None, abivarname=None, executables=None, varset=None, range=None, excludes=None, requires=None,
                 commentdefault=None, commentdims=None, topics=None):
        self.vartype = vartype
        self.characteristics = characteristics
        self.mnemonics = mnemonics
        self.dimensions = dimensions
        self.defaultval = defaultval
        self.text = literal(text)
        self.abivarname = abivarname
        self.varset = varset
        self.executables = executables
        self.commentdefault = commentdefault
        self.commentdims = commentdims
        self.range = range
        self.topics = topics
        self.excludes = excludes
        self.requires = requires

    @property
    def name(self):
        """Name of the variable without the executable name."""
        return self.abivarname if "@" not in self.abivarname else self.abivarname.split("@")[0]

    @property
    def executable(self):
        """str with the name of the code associated to this variable."""
        if "@" in self.abivarname:
            code = self.abivarname.split("@")[1]
            assert code == self.varset
        else:
            code = "abinit"
        return code

    @property
    def absolute_url(self):
        #docs.abinit.org/vardocs/CODENAME/VARNAME?version=8.6.2
        return "https://www.docs.abinit.org/vardoc/%s/%s" % (self.executable, self.name)

    @property
    def topic_tribes(self):
        """topic --> list of tribes"""
        if hasattr(self, "_topic_tribes"):
            return self._topic_tribes
        else:
            assert self.topics is not None
            od = OrderedDict()
            #for tok in self.topics.split(","):
            for tok in self.topics:
                topic, tribe = [s.strip() for s in tok.split("_")]
                if topic not in od: od[topic] = []
                od[topic].append(tribe)
            self._topic_tribes = od
            return od

    @property
    def is_internal(self):
        """True if this is an internal variable."""
        return self.characteristics is not None and '[[INTERNAL_ONLY]]' in self.characteristics

    @property
    def mdlink(self):
        """Abinit wikilink."""
        return "[[%s:%s]]" % (self.executable, self.name)

    def __str__(self):
        return self.to_string()

    # MG new code
    #def __repr__(self):
    #    return "Variable " + str(self.abivarname) + " (default = " + str(self.defaultval) + ")"

    def to_string(self, verbose=0):
        return "Variable " + str(self.abivarname) + " (default = " + str(self.defaultval) + ")"

    def __hash__(self):
        # abivarname is unique
        return hash(self.abivarname)

    def __eq__(self, other):
        if other is None: return False
        return self.abivarname == other.abivarname

    def __ne__(self, other):
        return not (self == other)

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

    def to_markdown(self, with_hr=True):
        lines = []; app = lines.append

        app("## **%s** \n\n" % self.name)
        app("*Mnemonics:* %s  " % str(self.mnemonics))
        if self.characteristics:
            app("*Characteristics:* %s  " % ", ".join(self.characteristics))
        if self.topic_tribes:
            app("*Mentioned in topic(s):* %s  " % ", ".join("[[topic:%s]]" % k for k in self.topic_tribes))
        app("*Variable type:* %s  " % str(self.vartype))
        if self.dimensions:
           app("*Dimensions:* %s  " % format_dimensions(self.dimensions))
        if self.commentdims:
            app("*Commentdims:* %s  " % self.commentdims)
        app("*Default value:* %s  " % self.defaultval)
        if self.commentdefault:
            app("*Comment:* %s  " % self.commentdefault)
        if self.requires:
            app("*Only relevant if:* %s  " % str(self.requires))
        if self.excludes:
            app("*The use of this variable forbids the use of:* %s  " % self.excludes)

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
        if self.text is not None:
            #md_text = html2text(self.text)
            md_text = self.text
            app(2 * "\n")
            app(md_text)
        else:
            raise ValueError("Variable `%s` does not have description text!" % self.name)

        if with_hr:
            app("* * *" + 2*"\n")

        return "\n".join(lines)

    def validate(self):
        """Validate variable. Raises ValueError if not valid."""
        errors = []
        eapp = errors.append

        svar = "None"
        try:
            svar = str(self)
        except Exception as exc:
            eapp(str(exc))

        if self.abivarname is None:
            eapp("Variable `%s` has no name" % svar)

        if self.vartype is None:
            eapp("Variable `%s` has no vartype" % svar)
        elif not self.vartype in ("integer", "real", "string"):
            eapp("%s must have vartype in ['integer', 'real', 'string'].")

        if self.topics is None:
            eapp("%s does not have at least one topic and the associated relevance" % svar)

        """
        topics_name_relevance = var.topics.split(',')
        for topic_name_relevance in topics_name_relevance:
              name_relevance = topic_name_relevance.split('_')
              if not name_relevance[0].strip() in topics:
                    print('FAIL: ', abivarname, ' delivers topicname_relevance ',name_relevance,
                          ' with topicname ',name_relevance[0].strip(),' that does not belong to the allowed list')
                    retcode += 1
              if not name_relevance[1].strip() in relevance_names:
                    print('FAIL: ', abivarname, ' delivers topicname_relevance ',name_relevance,
                          ' with relevance ',name_relevance[1].strip(),' that does not belong to the allowed list')
                    retcode += 1
        """

        if self.characteristics is not None:
            if not isinstance(self.characteristics, list):
                eapp("The field characteristics of %s is not a list" % svar)
            #else:
            #    for cat in self.characteristics:
            #        if cat.replace("[[", "").replace("]]", "") not in characteristics:
            #            eapp(The characteristics %s  of %s is not valid" % (cat, self))

        if self.dimensions is None:
            eapp("%s does not have a dimension. If it is a *scalar*, it must be declared so." % svar)
        else:
            if self.dimensions != "scalar":
                if not isinstance(self.dimensions, (list, ValueWithConditions)):
                    eapp('The dimensions field of %s is not a list neither a valuewithconditions' % svar)

        if self.varset is None:
            eapp('`%s` does not have a varset' % svar)
        #else:
        #    if not isinstance(self.varset, str) or self.varset not in varset_names:
        #        print('The field varset of %s should be one of the valid varsets' % str(self))

        if len(self.name) > 20:
            eapp("Lenght of `%s` is longer than 20 characters." % svar)

        if errors:
            raise ValueError("\n".join(errors))


####################################################################################################

class Components(yaml.YAMLObject):
    name = None  # String containing section name
    keyword = '' # String containing the short description of the topics, to be echoed in the title of the section file.
    authors = '' # String containing the list of authors.
    howto  = ''  # For the 'topics' Should complete the sentence beginning with "How to"
    header = ''  # Header of the file, possibly the 'default' one
    title  = ''  # Title  of the file, possibly the 'default' one
    subtitle  = ''  # Subtitle  of the file, possibly the 'default' one
    purpose   = ''  # Purpose  of the file, possibly the 'default' one
    advice    = ''  # Advice  of the file, possibly the 'default' one
    copyright = ''  # Copyright of the file, possibly the 'default' one
    introduction = ''  # Introduction
    links     = ''  # Links of the file, possibly the 'default' one
    menu      = ''  # Menu of the file, possibly the 'default' one
    tofcontent_header      = ''  # Header of the table of content of the file, possibly the 'default' one
    tutorials    = '' # List of relevant tutorials
    examples     = '' # Relevant examples
    end       = ''  # End of the file, possibly the 'default' one

    yaml_tag = u'!components'

    def attrs(self):
        return ['name', 'keyword', 'authors', 'howto', 'header', 'title', 'subtitle', 'purpose', 'advice',
                'copyright', 'introduction', 'links', 'menu', 'tofcontent_header', 'tutorials', 'examples', 'end']

    #Note that the default values are actually not initialized here, but in the data file, in order to ease the maintenance.
    def __init__(self, name=None, keyword=None, authors=None, howto=None, header=None, title=None, subtitle=None, purpose=None, advice=None,
                 copyright=None, links=None, menu=None, tofcontent_header=None, tutorials=None, examples=None, end=None):
        self.name = name
        self.keyword = keyword
        self.authors = authors
        self.howto = howto
        self.header = header
        self.title  = title
        self.subtitle = subtitle
        self.purpose  = purpose
        self.advice   = advice
        self.copyright= copyright
        self.introduction= introduction
        self.links    = links
        self.menu     = menu
        self.tofcontent_header = tofcontent_header
        self.tutorials = tutorials
        self.examples = examples
        self.end      = end

####################################################################################################

#class ValueWithUnit(yaml.YAMLObject):
class ValueWithUnit(object):
    value = None
    units = None
    yaml_tag = u'!valuewithunit'

    def __init__(self, value=None, units=None):
        self.value = value
        self.units = units

    def __str__(self):
        return str(self.value) + " " + str(self.units)

    def __repr__(self):
        return str(self)

####################################################################################################

class Range(yaml.YAMLObject):
    """
    Specifies a range (start:stop:step)
    """
    start = None
    stop = None

    yaml_tag = u'!range'

    def __init__(self, start=None, stop=None):
        self.start = start
        self.stop = stop

    def isin(self, value):
        isin = True
        if self.start is not None:
            isin = isin and (self.start <= self.value)
        if stop is not None:
            isin = isin and self.stop > self.value
        return str(self)

    def __repr__(self):
        # Add whitespace after `[` or before `]` to avoid [[[ and ]]] pattersn
        # that enter into conflict with wikiling syntax [[...]]
        if self.start is not None and self.stop is not None:
            return "[ " + str(self.start) + " ... " + str(self.stop) + " ]"
        if self.start is not None:
            return "[ " + str(self.start) + "; ->"
        if self.stop is not None:
            return "<-;" + str(self.stop) + " ]"
        else:
            return None

####################################################################################################

#class ValueWithConditions(yaml.YAMLObject):
class ValueWithConditions(dict):
    """
    Used for variables whose value depends on a list of conditions.

    .. example:

        ValueWithConditions({'[[paral_kgb]]==1': '6', 'defaultval': 2}),

        Means that the variable is set to 6 if paral_kgb == 1 else 2
    """
    yaml_tag = u'!valuewithconditions'

    def __repr__(self):
        s = ''
        for key in self:
            if key != 'defaultval':
                s += str(self[key]) + ' if ' + str(key) + ',\n'
        s += str(self["defaultval"]) + ' otherwise.\n'
        return s

    def __str__(self):
        return self.__repr__()

####################################################################################################

#class MultipleValue(yaml.YAMLObject):
class MultipleValue(dict):
    """
    Used for variables that can assume multiple values.

    abivarname="istwfk",
    defaultval=MultipleValue({'number': None, 'value': 0}),
    """

    number = None
    value = None

    yaml_tag = u'!multiplevalue'

    def __init__(self, number=None, value=None):
        self.number = number
        self.value = value

    def __repr__(self):
        if self.number is None:
            return "*" + str(self.value)
        else:
            return str(self.number) + "*" + str(self.value)


# My Code
from collections import OrderedDict

_VARS = None

def get_variables_code(use_pymods=True):
    """
    Return the database of variable and cache it. Main entry point for client code.
    """
    global _VARS
    if _VARS is None:
        if use_pymods:
            _VARS = VarDatabase.from_pyfiles()
        else:
            yaml_path = os.path.join(os.path.dirname(__file__), "..", "doc", "variables", "origin_files", "abinit_vars.yml")
            _VARS = VarDatabase.from_file(yaml_path)

    return _VARS


class VarDatabase(OrderedDict):
    """
    This object stores the full set of input variables for all the Abinit executables.
    in a dictionary mapping the name of the code to a subdictionary of variables.
    """
    @classmethod
    def from_pyfiles(cls, dirpath=None):
        """Initialize the object from python modules."""
        if dirpath is None:
            dirpath = os.path.dirname(os.path.abspath(__file__))
        pyfiles = [os.path.join(dirpath, f) for f in os.listdir(dirpath) if
                   f.startswith("variables_") and f.endswith(".py")]
        new = cls()
        #for exname in sorted(set(v.executable for v in vlist)):
        for pyf in pyfiles:
            print("pyf", pyf)
            vd = InputVariables.from_pyfile(pyf)
            #items = [(v.name, v) for v in vlist if v.executable == exname]
            #vd = InputVariables(sorted(items, key=lambda t: t[0]))
            #vd.executable = exname
            #vd.all_varset = sorted(set(v.varset for v in vd.values()))
            new[vd.executable] =  vd

        # Read list of strings with possible character of variables.
        yaml_path = "/Users/gmatteo/git_repos/abinit/doc/variables/origin_files/"
        with io.open(os.path.join(os.path.dirname(yaml_path), "characteristics.yml"), "rt", encoding="utf-8") as f:
            new.characteristics = yaml.load(f)

        # Read list of `external_params` i.e. external parameters that are not input variables,
        # but that are used in the documentation of other variables
        # then convert to dict {name --> description}
        with io.open(os.path.join(os.path.dirname(yaml_path), "list_externalvars.yml"), "rt", encoding="utf-8") as f:
            d = {k: v for k, v in yaml.load(f)}
            new.external_params = OrderedDict([(k, d[k]) for k in sorted(d.keys())])

        return new

    @classmethod
    def from_file(cls, yaml_path):
        with io.open(yaml_path, 'rt', encoding="utf-8") as f:
            vlist = yaml.load(f)

        new = cls()
        # Read list of strings with possible character of variables.
        with io.open(os.path.join(os.path.dirname(yaml_path), "characteristics.yml"), "rt", encoding="utf-8") as f:
            new.characteristics = yaml.load(f)

        # Read list of `external_params` i.e. external parameters that are not input variables,
        # but that are used in the documentation of other variables
        # then convert to dict {name --> description}
        with io.open(os.path.join(os.path.dirname(yaml_path), "list_externalvars.yml"), "rt", encoding="utf-8") as f:
            d = {k: v for k, v in yaml.load(f)}
            new.external_params = OrderedDict([(k, d[k]) for k in sorted(d.keys())])

        for exname in sorted(set(v.executable for v in vlist)):
            items = [(v.name, v) for v in vlist if v.executable == exname]
            vd = InputVariables(sorted(items, key=lambda t: t[0]))
            vd.executable = exname
            vd.all_varset = sorted(set(v.varset for v in vd.values()))
            new[exname] =  vd

        return new

    def iter_allvars(self):
        """Iterate over all variables. Flat view."""
        for vd in self.values():
            for var in vd.values():
                yield var

    #def get_version_endpoints(self):
    #    """
    #    docs.abinit.org/vardocs/abinit/asr?version=8.6.2
    #    """
    #    code_urls = {}
    #    for codename, vard in self.items()
    #        code_urls[codename] = d = {}
    #        for vname, var in var.items():
    #            d[vname] = var.url
    #    return version, code_urls

    #def update_json_endpoints(self, json_path, indent=4):
    #    import json
    #    with open(json_path, "rt") as fh:
    #        oldd = json.load(fh)
    #    new_version, newd = self.get_version_endpoints()
    #    assert new_version not in oldd
    #    oldd[new_version] = newd
    #    with open(json_path, "wt") as fh:
    #        json.dump(oldd, fh, indent=indent)

    def write_pymods(self, dirpath="."):
        """
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

        from html2text import html2text
        #for code in ["abinit", ]:
        for code in self:
            varsd = self[code]

            lines = ["""\
from __future__ import print_function, division, unicode_literals, absolute_import

# We don't install with setup.py hence we have to add this directory to sys.path
#import sys
#import os
#sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from abimkdocs.variables import ValueWithUnit, MultipleValue, Range
ValueWithConditions = dict

Variable=dict\nvariables = ["""
]
            for name in sorted(varsd.keys()):
                var = varsd[name]
                text = html2text(var.text)
                text = '"""\n' + text.rstrip() + '\n"""'
                assert var.range is None
                assert var.executables is None
                #executables=None,
                #range=None,
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
    commentdefault={commentdefault},
    commentdims={commentdims},
    excludes={excludes},
    requires={requires},
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
          text=text,
          defaultval=defaultval2arg(var.defaultval),
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
    """
    @classmethod
    def from_pyfile(cls, filepath):
        import imp
        module = imp.load_source(filepath, filepath)
        new = cls()
        new.executable = module.executable
        vlist = [Variable(**d) for d in module.variables]
        vlist = sorted(vlist, key=lambda v: v.name)
        for v in vlist:
            new[v.name] = v
        new.all_varset = sorted(set(v.varset for v in new.values()))
        return new

    def groupby_first_letter(self):
        keys = sorted(self.keys(), key=lambda n: n[0].upper())
        od = OrderedDict()
        for char, group in groupby(keys, key=lambda n: n[0].upper()):
            od[char] = [self[name] for name in group]
        return od

    def get_vartabs_html(self, website, page_rpath):
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

    def get_plotly_networkx(self, varset="all", layout_type="spring", include_plotlyjs=False):
        # https://plot.ly/python/network-graphs/
        # Build the graph
        import networkx as nx
        g, edge_labels = nx.Graph(), {}
        for i, (name, var) in enumerate(self.items()):
            if varset != "all" and var.varset != varset: continue
            parents = var.get_parent_names()
            if not parents: continue
            g.add_node(var, name=name)
            for parent in parents:
                #print(parent, "is parent of ", name)
                parent = self[parent]
                g.add_edge(parent, var)

                # TODO: Add getters! What about locked nodes!
                #i = [dep.node for dep in child.deps].index(task)
                #edge_labels[(task, child)] = " ".join(child.deps[i].exts)

        # Get positions for all nodes using layout_type.
        # e.g. pos = nx.spring_layout(g, iterations=80)
        pos = getattr(nx, layout_type + "_layout")(g) #, scale=100000, iterations=30)
        #pos = nx.spring_layout(g, k=2, iterations=100)

        # Add edges as disconnected lines in a single trace and nodes as a scatter trace
        import plotly
        import plotly.graph_objs as go
        from textwrap import fill
        edge_trace = go.Scatter(
            x=[],
            y=[],
            line=go.Line(width=2.0, color='#888', dash="dot"),
            hoverinfo='none',
            mode='lines',
        )

        for edge in g.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_trace['x'] += [x0, x1, None]
            edge_trace['y'] += [y0, y1, None]

        node_trace = go.Scatter(
            x=[],
            y=[],
            text=[v.name for v in g],
            textposition='bottom',
            mode='markers+text',
            hoverinfo='text',
            #hovertext=[fill(html2text(v.text), width=90) for v in g.nodes()],
            hovertext=[v.mnemonics for v in g.nodes()],
            marker=go.Marker(
                showscale=True,
                # colorscale options
                # 'Greys' | 'Greens' | 'Bluered' | 'Hot' | 'Picnic' | 'Portland' |
                # Jet' | 'RdBu' | 'Blackbody' | 'Earth' | 'Electric' | 'YIOrRd' | 'YIGnBu'
                colorscale='YIGnBu',
                reversescale=True,
                color=[],
                size=20,
                colorbar=dict(
                    thickness=15,
                    title='Node Connections',
                    xanchor='left',
                    titleside='right'
                ),
                line=dict(width=2)))

        for node in g.nodes():
            x, y = pos[node]
            node_trace['x'].append(x)
            node_trace['y'].append(y)

        # Color node points by the number of connections.
        for node, adjacencies in enumerate(g.adjacency_list()):
            node_trace['marker']['color'].append(len(adjacencies))
            node_info = '# of connections: '+str(len(adjacencies))
            #node_trace['text'].append(node_info)

        fig = go.Figure(data=go.Data([edge_trace, node_trace]),
                     layout=go.Layout(
                        title='<br>Network graph for varset %s' % varset,
                        titlefont=dict(size=16),
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20,l=5,r=5,t=40),
                        annotations=[ dict(
                            text="Network for varset %s" % varset,
                            showarrow=False,
                            xref="paper", yref="paper",
                            x=0.005, y=-0.002 ) ],
                        xaxis=go.XAxis(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=go.YAxis(showgrid=False, zeroline=False, showticklabels=False)))

        #py.iplot(fig, filename='networkx')
        #plotly.offline.plot(fig)
        s = plotly.offline.plot(fig, include_plotlyjs=include_plotlyjs, output_type='div')
        #print(s)

        return s

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
                URL=var.absolute_url,
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

    def get_graphviz(self, vartype=None, engine="automatic", graph_attr=None, node_attr=None, edge_attr=None):
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
        #if vartype
        #var = self[varname]

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
                URL=var.absolute_url,
                target="_top",
                tooltip=str(var.mnemonics),
            )

        edge_kwargs = dict(arrowType="vee", style="solid")
        with_children = False

        for name, var in self.items():
            if vartype is not None and var.vartype != vartype: continue

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


def splitall(path):
    """Return list with all components of a path."""
    import os, sys
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