#!/usr/bin/env python
"""Extract the names of the entities defined in a fortran datatype."""
from __future__ import print_function, division

import sys
import re


def remove_text_in_parenthesis(s):
    return re.sub(r'\([^)]*\)', '', s)


class FortypeParserError(Exception):
    """Exceptions raised by FortypeParser."""


class FortypeParser(object):
    """
    Extracts the name of the variables defined in a Fortran datatype
    Assumes plain F90 format without CPP macros.
    """
    Error = FortypeParserError

    # Name of the type ("type xxxx", "type :: xxxx" or "type, yyy :: xxxx"
    re_name_type = re.compile('^[ \t]*type[ ]*(|.*::)[ ]*(?P<name>\w+)', re.IGNORECASE)

    # Detect "end type"
    re_end_type = re.compile('^[ \t]*end[ ]+type[ ]*(?P<name>\w+)', re.IGNORECASE)

    def __init__(self, filepath, ftype_name):
        """
        Args:
            filepath:
                Fortran file with the declaraion of the datatype
            ftype_name:
                Name of the datatype
        """
        # Open the file and find the section with the type declaration.
        indec, self.declaration = 0, []
        with open(filepath) as fh:
            for line in fh:
                line = line.lower().strip()
                if not line or line.startswith("!"):
                    continue

                m = self.re_name_type.match(line)
                if m and m.groups()[1] == ftype_name:
                    #print(m, m.groups())
                    indec = 1

                if indec:
                    self.declaration.append(line)
                    m = self.re_end_type.match(line)
                    if m:
                        #print(m.groups())
                        if m.groups()[0] != ftype_name:
                            raise self.Error("Cannot find end to type declaration!")
                        break

        if not self.declaration:
            raise self.Error("Cannot find declaration of Fortran type: %s" % ftype_name)

        # Here we assume abinit coding rules e.g. real :: var
        # and we use "::" to find the section with the variables.
        self.varnames, badlines = [], []
        for dec in self.declaration[1:-1]:
            s = dec.find("::")
            if s == -1:
                badlines.append(dec)
            else:
                # Strip comments.
                dec = dec[s+2:]
                e = dec.find("!")
                if e != -1:
                    dec = dec[:e]

                # Dec is the string containing our variables.
                # Remove text inside (..) to avoid problems with arrays and split using ,
                # print("dec: ", dec)
                tokens = remove_text_in_parenthesis(dec)
                tokens = [s.strip() for s in tokens.split(",")]

                self.varnames.extend(tokens)

        if badlines:
            print(badlines)
            err_msg = "You are not following abinit coding rules"
            raise ValueError(err_msg)

    def __str__(self):
        return "\n".join(self.varnames)

    def json_dump(self, stream=sys.stdout):
        import json
        json.dump(self.varnames, stream, indent=4)


def main():
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('filepath', help="Fortran file with the declaration of the datatype")
    parser.add_argument('datatype', help="Name of the datatype to be analyzed")

    options = parser.parse_args()

    p = FortypeParser(options.filepath, options.datatype)
    #import pprint
    #pprint.pprint(p.varnames)
    p.json_dump()

    return 0


if __name__ == "__main__":
    sys.exit(main())
