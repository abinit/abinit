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

lbox = Listbox(height = 40, width = 90, returnExit = 1)
lbox.append("Fedora", 1)
lbox.append("Red Hat Enterprise Linux", 2)
lbox.append("Ubuntu", 3)
lbox.append("Slackware", 4)
lbox.append("CentOS", 5)

grid = GridForm(screen, "Select your favorite ditro", 1, 1)
grid.add(lbox, 0, 0)

result = grid.runOnce()

screen.finish()

#print "listbox:", lbox.current()
if lbox.current() == 1:
	print "Selected Fedora!"
elif lbox.current() == 2:
	print "Selected Red Hat Enterprise Linux!"
elif lbox.current() == 3:
	print "Selected Ubuntu!"
elif lbox.current() == 4:
	print "Selected Slackware!"
elif lbox.current() == 5:
	print "Selected CentOS!"
