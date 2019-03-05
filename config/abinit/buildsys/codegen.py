#!/usr/bin/python3

import os
import re
import sys

                    ########################################

#
# Generic template handler
#

class AbinitTemplate(object):

    def __init__(self, template_text, start_tag="@%", end_tag="%@"):
        """
            Creates an Abinit Template from a text string or file, with the
            specified pattern delimiters. Patterns are expected to contain
            lower-case letters, digits and underscores only. For instance,
            "action1_ok" is a valid pattern, while "Action1-ok" --- which
            contains a capital letter and a dash --- is not.

            This constructor will raise a ValueError if:
                * the delimiters are improperly specified;
                * no pattern is found in the template.

            :param template_text: text string or text file containing a
                                  template with patterns
            :param delimiters: list of length 2 containing an opening and
                               a closing delimiter for patterns.
        """

        # Read configuration from file or text string
        re_newline = re.compile(r"\n", flags=re.MULTILINE+re.DOTALL)
        if ( (not re_newline.search(template_text)) and \
             os.access(template_text, os.R_OK) ):
            with open(template_text, "r") as tpl_file:
              self.template = tpl_file.read()
        else:
            self.template = template_text

        # Set delimiters
        self.start_tag = start_tag
        self.end_tag = end_tag

        # Prepare regular expressions
        self.re_patterns = re.compile("%s([a-z0-9_]+)%s" % \
            (self.start_tag, self.end_tag), flags=re.MULTILINE)

        # Extract patterns within the template
        patterns = self.re_patterns.findall(self.template)
        if ( patterns ):
            self.patterns = list(set(patterns))
            self.patterns.sort()
        else:
            raise ValueError("no pattern found in template")


    def get_missing(self, patterns):
        """
            Returns the elements missing in the provided patterns to fully
            substitute the patterns of the template.

            :param patterns: a list or dictionary containing patterns
            :return: a list, or None if all elements are found
        """

        # Accept both lists and dicts
        if ( isinstance(patterns, dict) ):
          tmp_patterns = patterns.keys()
        else:
          tmp_patterns = patterns

        missing = [item for item in self.patterns if not item in tmp_patterns]

        if ( len(missing) > 0 ):
            return missing
        else:
            return None


    def get_patterns(self):
        """
            Returns the patterns defined in the template.

            :return: a list of strings
        """

        return self.patterns


    def get_undefined(self, patterns):
        """
            Returns the elements present in the provided patterns but not
            defined in the template.

            :param patterns: a list or dictionary containing patterns
            :return: a list, or None if all elements are found
        """

        # Accept both lists and dicts
        if ( isinstance(patterns, dict) ):
          tmp_patterns = patterns.keys()
        else:
          tmp_patterns = patterns

        mismatch = [item for item in tmp_patterns if not item in self.patterns]

        if ( len(mismatch) > 0 ):
            return mismatch
        else:
            return None


    def substitute(self, replacements, err_miss=False, err_undef=True):
        """
            Returns a text string resulting from the substitution of the
            template patterns with the provided replacements. By default,
            an exception will be raised if there is a mismatch between
            the template and the replacements, but this behavior can be
            tuned.

            :param replacements: a dictionary for pattern substitution
            :param err_miss: whether to throw an exception if the
                             replacements are incomplete
            :param err_undef: whether to throw an exception if the
                              replacements contain patterns undefined
                              in the template
            :return: the substituted text
        """

        # Check whether replacements perfectly match template patterns
        if ( err_miss ):
            pat_miss = self.get_missing(replacements)
            if ( pat_miss ):
                raise ValueError(
                    "incomplete replacements provided:\n    %s" % pat_miss)
        if ( err_undef ):
            pat_undef = self.get_undefined(replacements)
            if ( pat_undef ):
                raise ValueError(
                    "provided replacements contain undefined patterns:\n    %s" % \
                    pat_undef)

        # Substitute all provided patterns
        ret = self.template
        for (key, val) in replacements.items():
            ret = re.sub("%s%s%s" % (self.start_tag, key, self.end_tag),
                val, ret)

        return ret

