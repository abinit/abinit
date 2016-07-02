#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import sys,os
import string
from PyQt4 import Qt,QtGui,QtCore




#---------------------------------------------------#
#---------------------------------------------------#
#----------------BUTTON CLASS-----------------------#
#---------------------------------------------------#
#---------------------------------------------------#
class Button(QtGui.QPushButton):

   #-------------Constructor--------------#
   def __init__(self,pname,pnumber,parent = None):
       QtGui.QPushButton.__init__(self,pname, parent)
       self.par = parent
       self.nb = pnumber
	
   #------------Methods--------------------#	
   def clic(self):
       self.emit(QtCore.SIGNAL("change(int)"),self.nb)
       
       
   def getNumber(self):
       return self.nb
