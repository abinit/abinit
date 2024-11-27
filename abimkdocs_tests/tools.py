from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import os

from unittest import TestCase

_PATCH_DONE = False

def patch_syspath():
    global _PATCH_DONE
    if _PATCH_DONE: return
    # We don't install with setup.py hence we have to add the directory [...]/abinit/tests to $PYTHONPATH
    pack_dir = os.path.dirname(os.path.abspath(__file__))
    pack_dir = os.path.join(pack_dir, "..")
    sys.path.insert(0, pack_dir)

    # This needed to import doc.tests
    sys.path.insert(0, os.path.join(pack_dir, "doc"))
    _PATCH_DONE = True

class AbimkdocsTest(TestCase):

    @staticmethod
    def get_abinit_varnames_from_f90():
        # construct list of input keywords that appear in chkvars.F90
        home_dir = os.path.join(os.path.dirname(__file__) , "..")
        path = os.path.join(home_dir, "src/44_abitypes_defs/m_dtset.F90")

        in_block = False
        words = []
        with open(path, "rt") as fh:
            for line in fh:
                if line.find("admitted variable names") > 0: in_block = True
                if line.find("Extra token") > 0: in_block = False
                if in_block and line.find("list_var") > 0:
                    line_words = (line.split("'")[1]).split()
                    for i in range(len(line_words)):
                        words.append(line_words[i])

        if not words:
            print("Found empty list of words in %s " % path)
            print("Perhaps someone changed the format of the file?")
            print("Please modify the code in " + __file__)
            raise RuntimeError("")

        return set(words)

    @staticmethod
    def get_anaddb_varnames_from_f90():
        # Scan the source and search for the calls to intagm. Parse the arguments
        # and extract the name of the variable. The prototype of intagm is:
        #    call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'brav',tread,'INT')
        import re
        re_call = re.compile(r'\s*call\s+intagm\((.+)\)\w*', re.I)

        # construct list of key words appearing in anaddb input
        home_dir = os.path.join(os.path.dirname(__file__) , "..")
        path = os.path.join(home_dir, "src/77_ddb/m_anaddb_dataset.F90")

        words = []
        with open(path, "rt") as fh:
            for line in fh:
                m = re_call.match(line)
                if m:
                  tokens = m.group(1).split(",")
                  assert len(tokens) == 9
                  words.append(tokens[-3].replace("'","").replace('"',""))

            if not words:
                print("Found empty list of words in file %s" % path)
                print("Perhaps someone changed the format of the file?")
                print("Please modify the code in " + __file__)
                raise RuntimeError()

        return set(words)
