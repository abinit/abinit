#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2014
"""

import sys,os,time,commands
import string, math

#GUI 
import gui.graph as Graph
import gui.conv  as Conv

#Utility
import utility.writeHIST as Write
import utility.analysis as Analysis

try:
    from PyQt4 import Qt,QtGui,QtCore
except:
    pass;

from numpy import *

#----------------------------------------------------------------#
#--------------------------NETCDF WRITER-------------------------#
#----------------------------------------------------------------#
class winNetcdf(QtGui.QWidget):

    PTOE = Analysis.PeriodicTableElement()

    def __init__(self, file, parent = None,name =''):

        self.file = file
        self.name = name
        self.ni   = 1
        self.nf   = 2
        self.initUI(parent)
        self.raise_()

    def initUI(self, parent):

        #-----------------Creation of the windows----------------------------#
        QtGui.QWidget.__init__(self, parent)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle(self.name + ' Save Netcdf')
        self.setFixedSize(600, 450)
        self.center()
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)

        self.lbl1 = QtGui.QLabel(" Name :", self)
        self.lbl1.setFixedWidth(95)

        self.lname = QtGui.QLineEdit()
        self.lname.setFixedWidth(200)
        
        self.pbClose = QtGui.QPushButton("Close")
        self.pbClose.setFixedSize(70,20)
        self.connect(self.pbClose,QtCore.SIGNAL("clicked()"),QtCore.SLOT('close()'))

        self.pbSave = QtGui.QPushButton("Save")
        self.pbSave.setFixedSize(70,20)
        self.connect(self.pbSave,QtCore.SIGNAL("clicked()"),self.save)

        self.layout.addWidget(self.lbl1    , 1, 0, 1, 1, QtCore.Qt.AlignRight)
        self.layout.addWidget(self.lname   , 1, 1, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.pbClose , 7, 0, 1, 2, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.pbSave  , 7, 1, 1, 4, QtCore.Qt.AlignCenter)

        self.show()
        #------------------------------------------------------------------------#

    def save(self):
        Write.writeHIST(self.file,self.name,self.ni,self.nf)

    def close(self):
        del self.graphMSD
        del self
    
    def closeEvent(self, event):
        try:
            del self.graphMSD
        except:
            pass
        try:
            del self
        except:
            pass     
        
    def center(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size =  self.geometry()
        self.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)
