from __future__ import division, print_function, unicode_literals, absolute_import

import sys
import os

from unittest import TestCase

class AbimkdocsTest(TestCase):
    """hello"""

_PATCH_DONE = False

def patch_syspath():
    global _PATCH_DONE
    if _PATCH_DONE: return
    # We don't install with setup.py hence we have to add the directory [...]/abinit/tests to $PYTHONPATH
    pack_dir = os.path.dirname(os.path.abspath(__file__))
    #pack_dir = os.path.dirname(os.path.abspath(filepath))
    pack_dir = os.path.join(pack_dir, "..")
    sys.path.insert(0, pack_dir)

    # This needed to import doc.tests
    sys.path.insert(0, os.path.join(pack_dir, "doc"))
    _PATCH_DONE = True
