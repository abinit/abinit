#! /usr/bin/python

#    Copyright (C) 2010-2020 ABINIT group
#
#    Written by Matthieu Verstraete in python (compatible v1.9).
#    This is free software, and you are welcome to redistribute it
#    under certain conditions (GNU General Public License,
#    see ~abinit/COPYING or http://www.gnu.org/copyleft/gpl.txt).
#
#    ABINIT is a project of the Universite Catholique de Louvain,
#    Corning Inc. and other collaborators, see ~abinit/doc/developers/contributors.txt.
#    Please read ~abinit/doc/biblio/generated_files/bib_acknow.html for suggested
#    acknowledgments of the ABINIT effort.
#
#    For more information, see https://www.abinit.org .
#
#  This script is to be used with the PHON code (or equivalent)
#  to calculate frozen phonon frequencies, free energies, etc...
#  It takes a DISP file of atomic displacements and an SPOSCAR file
#  with a with a supercell structure (VASP/PHON formats), and creates
#  the necessary lines for an abinit input file to calculate the
#  forces in displaced configurations.
#  See http://chianti.geol.ucl.ac.uk/~dario or
#  D. Alfe, Computer Physics Communications 180,2622-2633 (2009)
# 
#  NOTE: the symmetries in the present (1.28 8/2010) version of PHON
#  are not functioning properly in some cases. It is your own
#  responsibility to check it, and has nothing to do with ABINIT.
#
#  How to use:
#
#  1) run abinit for the relaxed structure with prtposcar 1.
#  this creates a reference XXX_POSCAR file, which you should rename
#  POSCAR.
#  2) run phon with LSUPER=.true.  or (e.g.)  phonopy -d --dim="2 2 3"
#     to create the SPOSCAR and DISP files
#  3) run this script (phondisp2abi.py).
#  4) copy script output to the abinit input file (removing duplicate
#  input variables etc...)
#  5) run abinit for each of the given datasets (and prtposcar 1 still)
#  6) concatenate the resulting XXX_FORCES files into one FORCES file
#  You also need to include the header lines for each displacement,
#  which are given by phondisp2abi.py in comments for each dataset
#  7) run phon again to get the desired phonons and properties.
#
#

import re
import string
import numpy
import numpy.linalg

#
#  convert PHON DISP and SPOSCAR files into ABINIT datasets with appropriately displaced atoms
#

fp_disp = open('DISP')
lines_disp = fp_disp.readlines()
fp_sposcar = open('SPOSCAR')
lines_sposcar = fp_sposcar.readlines()

# make unit cell input line
rprimd = numpy.zeros((3,3))
for idir  in range(3):
	line = lines_sposcar[2+idir]
	tokens = map(float,string.split(line))
	rprimd[0][idir] = tokens[0]
	rprimd[1][idir] = tokens[1]
	rprimd[2][idir] = tokens[2]


# get equilibirum positions
equilxred=[]
for line in lines_sposcar[7:]:
	equilxred.append(numpy.array(map(float,string.split(line))))

# output unit cell input line
print "# Add this to the abinit input file to do the PHON displacements"
print "#   given in the DISP file, with respect to the supercell in SPOSCAR"
print "#"
print "# Remember the POSCAR files have sorted the atomic types so the positions"
print "#   and displacements are now type ordered (fix typat, spinat, etc!)"
print "#"
print "ndtset ", len(lines_disp)
print "# supercell lattice vectors "
print "acell 1 1 1 Angstr"
print "rprim"
print " %24.14f %24.14f %24.14f" % (rprimd[0][0], rprimd[1][0], rprimd[2][0]) 
print " %24.14f %24.14f %24.14f" % (rprimd[0][1], rprimd[1][1], rprimd[2][1]) 
print " %24.14f %24.14f %24.14f" % (rprimd[0][2], rprimd[1][2], rprimd[2][2]) 


idtset=1
# for each displacement,
for line in lines_disp:
	tokens = string.split(line)
# get displacement in reduced coordinates
	iatom = int(tokens[1])
	dispred = numpy.array(map(float,tokens[2:5]))
# add displacement to correct atom
	xred = list(equilxred)
	xred[iatom-1] = xred[iatom-1] + dispred
	
# output xred for this dataset
	print "# add the following line, without the #, to the FORCES file for this dtset, when concatenating"
	print "# %d %24.14f %24.14f %24.14f" % (iatom, dispred[0], dispred[1], dispred[2])
	print "xred%d" % (idtset,)
	for xred_1at in xred:
		print " %24.14f %24.14f %24.14f" % (xred_1at[0], xred_1at[1], xred_1at[2])
# increment dataset counter
	idtset=idtset+1	
