#!/usr/bin/python
#
# Copyright (C) 2011 Douglas Schilling Landgraf <dougsland@redhat.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

from snack import *

screen = SnackScreen()

ret = EntryWindow(screen, 'Title', 'My super agenda', 
		['name', 'lastname', 'age'])


screen.finish()

status = ret[0]
values = ret[1]

# OK or Cancel
print "Pressed %s" % (status)

# Print every single item
for item in values:
	print item

