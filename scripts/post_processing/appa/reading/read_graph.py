#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import commands

class ReadGraph:
    
    # This class can read graph in a file
    # This is useful to add graphics

#-----------------------------#
#-------CONSTRUCTOR-----------#
#-----------------------------#

  def __init__(self, pnamefile):
    self.input_file = str(pnamefile) 
    self.number_column = 0
    file = open(self.input_file)

#-----------------------------#
#-----------METHODS-----------#
#-----------------------------#

  def read(self,c1,c2):
      file = open(self.input_file)
      x = []
      y = []
      idx1 = int(c1) - 1
      idx2 = int(c2) - 1
      for line in file:
          temp =  line.split()
          try:
            x.append(float(temp[idx1]))
          except:
            pass
          try:
            y.append(float(temp[idx2]))
          except:
            pass
      
      if len(x)==len(y):
        return [x,y]
      else:
        return [0,0]
  
  def getFile(self):
      file = open(self.input_file)
      return file.read()

  def getNbColumn(self):
      file = open(self.input_file)
      for line in file:
          test =  line.split()
          nb   =  len(test)
          if nb != self.number_column :
              self.number_column = nb
      file.close()
      return self.number_column 
