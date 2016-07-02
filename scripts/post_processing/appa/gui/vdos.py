#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import sys,os
import string, math

#gui
import gui.graph as Graph

#Utility
import utility.analysis as Analysis
try:
    from PyQt4 import Qt,QtGui,QtCore
except:
    pass
from numpy import sqrt,zeros,conjugate,arange,linspace,exp,log,sin


class winVDOS(QtGui.QWidget):

    def __init__(self, pVACF,  pDtion, parent = None,name =''):

        self.VACF = pVACF
        self.name = name
        self.dtion = pDtion
        self.initUI(parent)
        self.displayGraph()
        self.raise_()

    def initUI(self, parent):

        #-----------------Creation of the windows----------------------------#
        QtGui.QWidget.__init__(self, parent)
        self.setWindowTitle(self.name + ' option')
        self.setFixedSize(200, 150)
        self.center()
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)

        self.lbl = QtGui.QLabel("Gaussian Width", self)
        self.lbl.setFixedWidth(140)


        self.sbres = QtGui.QSpinBox()
        self.sbres.setMaximum(100)
        self.sbres.setMinimum(1)
        self.sbres.setValue(8)
        self.sbres.setFixedSize(70,20)
        self.connect(self.sbres,QtCore.SIGNAL('valueChanged(int )'),self.displayGraph)


        self.pbok = QtGui.QPushButton("close")
        self.pbok.setFixedSize(70,20)
        self.connect(self.pbok,QtCore.SIGNAL("clicked()"),QtCore.SLOT('close()'))

        self.checkbox =QtGui.QCheckBox("Show discret spectrum")
        self.connect(self.checkbox, QtCore.SIGNAL('clicked()'), self.showSprectrum)

        self.layout.addWidget(self.lbl   , 1, 0, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.sbres , 1, 1, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.checkbox , 2, 0, 1, 2, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.pbok  , 3, 0, 1, 2, QtCore.Qt.AlignCenter)

        self.show()
        #------------------------------------------------------------------------#


    def displayGraph(self):
        res = self.sbres.value()
        pdos = Analysis.DOS(self.VACF,res,self.dtion)

        y = pdos.getDOS()
        x = pdos.getFrequencies()
        
        try:
            self.GraphDOS.update(x,y,'E (meV)', "Phonons DOS (1/meV)",name = self.name, adjust=True)
        except:
            self.GraphDOS = Graph.graphic(x,y,'E (meV)', "Phonons DOS (1/meV)", average=False, adjust=True,name = self.name)
            self.connect(self.GraphDOS, QtCore.SIGNAL("myCustomizedSignal()"), self.close)
            self.GraphDOS.show()

        if self.checkbox.isChecked():
            self.addSprectum()


    def showSprectrum(self):

        if self.checkbox.isChecked():
            self.addSprectum()

        else:
            self.displayGraph()

    def addSprectum(self):
        res = self.sbres.value()
        pdos = Analysis.DOS(self.VACF,res,self.dtion)

        y = pdos.getSprectumDOS()
        x = pdos.getFrequencies()

        self.GraphDOS.addPlot(x,y,bar = True)

    def center(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size =  self.geometry()
        self.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)


    def close(self):
        del self.GraphDOS
        del self
    
    def closeEvent(self, event):
        try:
            del self.GraphDOS
        except:
            pass
        try:
            del self
        except:
            pass     







