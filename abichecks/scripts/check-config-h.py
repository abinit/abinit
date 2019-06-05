#!/usr/bin/env python
from __future__ import unicode_literals, division, print_function, absolute_import

import re
import os
import sys

from abirules_tools import find_src_dirs

IGNORED_DIRS = ["libpaw"]

# Init
re_srcfile = re.compile("\.(F|F90)$")
re_config  = re.compile(\
  "#if defined HAVE_CONFIG_H\n#include .config\.h.\n#endif\n",re.MULTILINE)


def main():
  for top in find_src_dirs():
      assert os.path.exists(top)
      retval = 0
      for root,dirs,files in os.walk(top):
        # Sort dirs
        dirs.sort()
        for item in files:
          # Skip some dirs
          ignored = 0
          for dd in IGNORED_DIRS:
            if root.find(dd) != -1: ignored=1

            if ignored == 0:
              if re_srcfile.search(item):
                path = os.path.join(root, item)
                with open(path, "rt") as fh:
                  src_data = fh.read()
                src_count = len(re.findall(re_config,src_data))

                if src_count == 0:
                  sys.stderr.write("%s/%s: missing include of config.h\n" % (root,item))
                  retval = 1
                elif src_count > 1:
                  sys.stdout.write("%s/%s: %d includes of config.h\n" % (root,item,src_count))

  return retval

if __name__ == "__main__":
  sys.exit(main())
