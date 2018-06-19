#! /usr/bin/python

#
#    Copyright (C) 2003-2018 ABINIT group
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
#  This module makes a file into a list of tokens, conserving
#   or eliminating comments of the type "#  bla bla   \n"
#  Probably redundant with other, more powerful python modules, but hey.
#


def tokenize_file(infilename):
	"""Take input file name and return list of tokens, eliminating comments"""
#  reg expss
	import re
	import sys
	import fileinput
	import math

#
#  open and parse input 
#

	try:
		fp = open (infilename, 'r')
	except IOError:
		print "Error opening file"
		raise

	lines = fp.readlines ()

#
#  put all tokens into tokens and remove comments
#
	tokens = []
	for line in lines:
		tmp = re.split ('[ \t\n]*',line)
#		print "tmp = ", tmp
		for tok in tmp:
			if (tok != ''):
				if (re.compile('[#!][.]*').match(tok)):
					break
				tokens.append(tok)
#	print "tokens = ", tokens

	fp.close()

	return tokens

def tokenize_file_keep_comments(infilename):
	"""Take input file name and return list of tokens, keep comments"""
#  reg expss
	import re
	import sys
	import fileinput
	import math

#
#  open and parse input 
#

	try:
		fp = open (infilename, 'r')
	except IOError:
		print "Error opening file"
		raise

	lines = fp.readlines ()

#
#  put all tokens into tokens and remove comments
#
	tokens = []
	for line in lines:
		tmp = re.split ('[ \t\n]*',line)
#		print "tmp = ", tmp
		for tok in tmp:
			if (tok != ''):
				tokens.append(tok)
#	print "tokens = ", tokens

	fp.close()

	return tokens
