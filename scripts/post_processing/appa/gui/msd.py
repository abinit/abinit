#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import sys,os,time,commands
import string, math

#GUI 
import gui.graph as Graph
import gui.conv  as Conv

#Utility
import utility.write as Write
import utility.analysis as Analysis

try:
    from PyQt4 import Qt,QtGui,QtCore
except:
    pass;

from numpy import *

#----------------------------------------------------------------#
#---------------WINDOWS-MEAN SQUARED DEPLACEMENT-----------------#
#----------------------------------------------------------------#
class winMSD(QtGui.QWidget):

    PTOE = Analysis.PeriodicTableElement()

    def __init__(self, file, parent = None,name =''):

        self.file = file
        self.name = name
        self.initUI(parent)
        self.displayGraph()
        self.raise_()

    def initUI(self, parent):

        #-----------------Creation of the windows----------------------------#
        QtGui.QWidget.__init__(self, parent)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle(self.name + ' MSD option')
        self.setFixedSize(200, 150)
        self.center()
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)

        self.lbl1 = QtGui.QLabel(" Atom type 1 :", self)
        self.lbl1.setFixedWidth(95)

        self.CBox1 = QtGui.QComboBox()
        self.CBox1.setFixedWidth(70)
        
        for i in range(len(self.file.getZnucl())):
            self.CBox1.addItem(str(self.PTOE.getName(self.file.getZnucl()[i])))
        self.connect(self.CBox1,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.displayGraph)        


        self.pbClose = QtGui.QPushButton("close")
        self.pbClose.setFixedSize(70,20)
        self.connect(self.pbClose,QtCore.SIGNAL("clicked()"),QtCore.SLOT('close()'))


        self.layout.addWidget(self.lbl1    , 1, 0, 1, 1, QtCore.Qt.AlignRight)
        self.layout.addWidget(self.CBox1   , 1, 1, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.pbClose , 7, 0, 1, 2, QtCore.Qt.AlignCenter)

        self.show()
        #------------------------------------------------------------------------#


    def displayGraph(self):

        atom = self.CBox1.currentIndex() + 1
        
        self.MeanSquaredDeplacement = Analysis.MSD(self.file,atom)

        x = self.MeanSquaredDeplacement.getX()
        y = self.MeanSquaredDeplacement.getMSD()
         
        try:
            self.graphMSD.update(x,y,'step', "Mean squared deplacement",name = self.name)
            self.graphMSD.addPlot(x,linspace(1,1,len(x)))
            self.graphMSD.show()
        except:
            self.graphMSD = Graph.graphic(x,y,'step', "Mean squared deplacement", average=False,name = self.name)
            self.connect(self.graphMSD, QtCore.SIGNAL("myCustomizedSignal()"), self.close)
            self.graphMSD.show()


    def update(self,pfile):
        self.file = pfile
        atom = self.CBox1.currentIndex() + 1
        try:
            self.MeanSquaredDeplacement = Analysis.MSD(self.file,atom)
            x = self.MeanSquaredDeplacement.getX()
            y = self.MeanSquaredDeplacement.getMSD()
            self.graphMSD.update(x,y,'step', "Mean squared deplacement",name = self.name)
            self.graphMSD.addPlot(x,linspace(1,1,len(x)))
        except:
            pass

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
