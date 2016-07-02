#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import commands

class Input:
    
    # This class can read input file for APPA
    # This is useful for run the no windozs version of APPA

#-----------------------------#
#-------CONSTRUCTOR-----------#
#-----------------------------#

  def __init__(self, pnamefile):
    self.input_file = str(pnamefile) 



#-----------------------------#
#-----------METHODS-----------#
#-----------------------------#

  def read(self):
    
    try:
      import re
      re_lib = True
    except:
      re_lib = False

    if re_lib :
      return self.expression_regular()
      


  def expression_regular(self):
    import re
    
    input = {}
    file = open(self.input_file)

    for line in file:
      #Remove the comment line
      m = re.search('^\s*#.*', line)
      if m:
        continue

      m = re.search('^\s*$', line)
      if m:
        continue

      m = re.search('^\s*file\s+(\S+)', line)
      if m:
        file = m.group(1)
        input[m.group(1)] = {}
        continue
      
      m = re.search('^\s*stepmin\s+(\S+)', line)
      if m:
        input[file]["stepmin"] = m.group(1)
        continue
      
      m = re.search('^\s*stepmax\s+(\S+)', line)
      if m:
        input[file]["stepmax"] = m.group(1)
        continue

      m = re.search('^\s*potential_energy\s+(\S*)\s*', line)
      if m:
        if m.group(1):
          input[file]['potential_energy'] = m.group(1)
        else:
          input[file]['potential_energy'] = 3
        continue

      m = re.search('^\s*total_energy\s+(\S*)\s*', line)
      if m:
        if m.group(1):
          input[file]['total_energy'] = m.group(1)
        else:
          input[file]['total_energy'] = 3         
        continue

      m = re.search('^\s*kinetic_energy\s+(\S*)\s*', line)
      if m:
        if m.group(1):
          input[file]['kinetic_energy'] = m.group(1)
        else:
          input[file]['kinetic_energy'] = 3         
        continue

      m = re.search('^\s*pressure\s+(\S*)\s*', line)
      if m:
        if m.group(1):
          input[file]['pressure'] = m.group(1)
        else:
          input[file]['pressure'] = 3         
        continue

      
      m = re.search('^\s*temperature\s+(\S*)\s*', line)
      if m:
        if m.group(1):
          input[file]['temperatures'] = m.group(1)
        else:
          input[file]['temperatures'] = 3 
        continue

      m = re.search('^\s*vacf\s+(\S*)\s*', line)
      if m:
        if m.group(1):
          input[file]['vacf'] = m.group(1)
        else:
          input[file]['vacf'] = 3 
        continue

      m = re.search('^\s*vdos\s+(\S*)\s*', line)
      if m:
        if m.group(1):
          input[file]['vdos'] = m.group(1)
        else:
          input[file]['vdos'] = 3 
        continue

      m = re.search('^\s*vdosres\s+(\S*)\s*', line)
      if m:
        if m.group(1):
          input[file]['vdosres'] = m.group(1)
        else:
          input[file]['vdosres'] = 8 
        continue


      m = re.search('^\s*stress\s+(\S*)\s*', line)
      if m:
        if m.group(1):
          input[file]['stress'] = m.group(1)
        else:
          input[file]['stress'] = 3 
        continue

      m = re.search('^\s*rdf\s+(\S*)\s*', line)
      if m:
        if m.group(1):
          input[file]['rdf'] = m.group(1)
        else:
          input[file]['rdf'] = 8 
        continue

      
      print "!!!Warning "+ line.replace('\n','')+ " is not a keyword"

    return input


#temp = commands.getoutput(' grep \'*\' '+self.input_file+' | awk \'{print $1,$2}\'  ')
    #print temp
