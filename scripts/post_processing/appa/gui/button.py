#!/usr/bin/env python

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

from PyQt4 import QtCore, QtGui


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
