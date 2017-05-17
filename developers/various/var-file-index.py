# -*- coding: us-ascii -*-
#------------------------------------------------------
#Build an index for all input variables
#inside the given input files of the ABINIT distribution.
#Date: 22/09/2004
#
#------------------------------------------------------
# Copyright (C) 2004-2017 ABINIT group (td)
# This file is distributed under the terms of the
# GNU General Public License, see ~abinit/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
#------------------------------------------------------

#Help message
def usage(error=0):
    print """
var_file_index [--help]
    Build the file 'Infos/varfileindex.html' which links
    the input files and the input variables.
"""
    sys.exit(error)

#Some import modules
import getopt
import os
import re
import sys
import time

version = time.strftime("%Y-%m-%d",time.gmtime())

#Detect options.
try:
    optlist, args = getopt.getopt(sys.argv[1:],'hu',["help"])
except getopt.error:
    print 'Error in arguments'
    usage(error=1)
    

for opt,a in optlist:
    if opt == '-h' or opt == '--help':
        usage(error=0)

if len(args) != 0:
    sys.stderr.write('Please no arguments\n')
    usage(error=1)

#Check if the directory is "Infos"
directory = os.path.basename(os.getcwd())

#if directory != "Infos":
#    sys.stderr.write("Works only in the 'Infos' directory.\n")
#    sys.exit(1)
    
#Source directory
SRC = ".."

#Regular expressions used
re_comment = re.compile('[!#].*')
re_keywords = re.compile('[a-zA-Z]+')

#Class for ABINIT
class project:
    "Class for ABINIT"
    def __init__(self,dir,pattern_dir='.*',pattern_file='.*',nofile=[],nokeyword=[]):
        "Initialisation"
        self.ROOT = dir
        self.keywords = dict()
        self.nofile = nofile
        self.nokeyword = nokeyword
        self.ndirs = 0
        self.nfiles = 0
        self.add(self.ROOT,pattern_dir,pattern_file)
    
    def add(self,dir,pattern_dir='.*',pattern_file='.*'):
        "Add keywords from files or directories"
        re_dir = re.compile(pattern_dir)
        re_file = re.compile(pattern_file)
        files = os.listdir(dir)
        for file in files:
            dd = "%s/%s" % (dir,file)
            if os.path.isdir(dd):
                if re_dir.match(file):
                    print "[%s:" % dd,
                    self.ndirs = self.ndirs + 1
                    #Add all files in this sub-directory
                    self.add(dd,pattern_dir,pattern_file)
                    print "]"
                else:
                    #We ignore it.
                    continue
            elif os.path.isfile(dd):
                if re_file.match(file) and file not in self.nofile:
                    self.add_file(dir,file)
    
    def add_file(self,dir,file):
        "Add the keywords of a file"
        self.nfiles = self.nfiles + 1
        print "(%d-%s)" % (self.nfiles,file),
        fd = open("%s/%s" % (dir,file))
        pathname = [file,dir]
        for line in fd:
            line = re_comment.sub('',line)
            list_keywords = re_keywords.findall(line)
            for keyword in list_keywords:
                keyword = keyword.lower()
                #Check if the keyword is toot small (D,E) or undesirable.
                if len(keyword) <= 1 or keyword in self.nokeyword:
                    continue
                if keyword not in self.keywords:
                    self.keywords[keyword] = []
                if pathname not in self.keywords[keyword]:
                    self.keywords[keyword].append(pathname)
        fd.close()

#List of files which we do not process
nofile = [ "ab.in", "t41.in", "t43.in", "t44.in", "t45.in", "t53.in", "t61.in", "t78.cut.in" ]
#List of keywords which we do not consider
nokeyword = [ "angstrom", "au", "ry", "hartree", "bohr" ]

#List of directories Test_...
abinit = project(SRC,pattern_dir="Test_",pattern_file=".*[.]in$",nofile=nofile,nokeyword=nokeyword)

print "Build the index of input variables"
print "Number of directories: %d" % abinit.ndirs
print "Number of files: %d" % abinit.nfiles

#Header of the html page
html_head = '''<?xml version="1.0"?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<head>
   <title>ABINIT: list of test input files.</title>
   <style type="text/css" media="all">@import "./formabinit.css";</style>
   <meta http-equiv="Content-Type" content="text/html; charset=us-ascii"/>

    <meta name="keywords" content="ab initio, DFT" />
    <meta name="description" content="List of test input files for ABINIT" />
    <meta name="revision" content="%s" />
</head>

<body>

<hr/>

<h1>Index of test input files using each input variable</h1>

<h3>
This document lists all the input files of the ABINIT Test directories 
according to the different input variables that are used. 
</h3>
<h4>Note that no distinction
is made between the different programs included in the ABINIT package :
input files (and variables) of the main ABINIT program
are mixed with those of MRGDDB, ANADDB, OPTIC, CONDUCTI, AIM, ...
</h3>

<h2>Alphabetical list of input variables:</h2>

''' % version

html_letter = ' '*3+'<p><b>%s.</b>\n'

html_tail_letter = ' '*3+'</p>\n'

html_keyword = ' '*6+'<a href="#%s">%s</a>&nbsp;\n'

html_second = '''
<hr/>

<h2>List of test input files for each input variable:</h2>

'''

#html to write the input variable.
html_head_var = ' '*3+'<h3><a name="%s">%s</a></h3>\n'+' '*3+'<p>\n'
html_tail_var = ' '*3+'</p>\n'

#html to write the file name
html_file = ' '*6+'<a href="%s">%s</a>&nbsp;\n'

#End of the html file.
html_tail = '''

</body>
</html>
'''

#Build the html page
allkeys = abinit.keywords.keys()
allkeys.sort()
fd = open("varfileindex.html","w")
fd.write(html_head)

#Build the alphabetic list
letter = ""
for keyword in allkeys:
    if letter == "":
        letter = keyword[0]
        fd.write(html_letter % letter.upper())
    elif letter != keyword[0]:
        letter = keyword[0]
        fd.write(html_tail_letter)
        fd.write(html_letter % letter.upper())
    fd.write(html_keyword % (keyword,keyword))
fd.write(html_tail_letter)

#Build all input file references for each input variables.
fd.write(html_second)
for keyword in allkeys:
    fd.write(html_head_var % (keyword,keyword))
    allfiles = abinit.keywords[keyword]
    allfiles.sort()
    for file in allfiles:
        pathname = "%s/%s" % (file[1],file[0])
        fd.write(html_file % (pathname,file[0]))
    fd.write(html_tail_var)

fd.write(html_tail)
fd.close()

sys.stderr.write("The file 'varfileindex.html' is ready.\n")

sys.exit(0)
