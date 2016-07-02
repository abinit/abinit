#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import sys,os,time,commands
import string, math, re
import threading

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
#---------------------WINDOWS-TOPOLOGY---------------------------#
#----------------------------------------------------------------#
class winRDF(QtGui.QWidget):

    PTOE = Analysis.PeriodicTableElement()

    def __init__(self, file, parent = None,name =''):

        self.file = file
        self.name = name
        self.initUI(parent)
        self.raise_()

        #-----------------Creation of the windows----------------------------#
    def initUI(self, parent):

        QtGui.QWidget.__init__(self, parent)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle(self.name + ' RDF option')

        typat = self.file.getTypat()
        ntypat = max(typat)
        nat = 0

        for i in enumerate(typat):
            nat = nat+1

        self.center()
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)

        self.mode = QtGui.QLabel("Calculation :", self)
        self.mode.setFixedWidth(85)

        self.sbm = QtGui.QComboBox()
        self.sbm.setFixedWidth(70)

        self.sbm.addItem(str('RDF'))
        self.sbm.addItem(str('ADF'))
        self.sbm.addItem(str('NDF'))
        self.sbm.addItem(str('Proba'))

        if ntypat == 1:

            self.lbl1 = QtGui.QLabel(" Atom type :", self)
            self.lbl1.setFixedWidth(80)
            
            self.CBox1 = QtGui.QComboBox()
            self.CBox1.setFixedWidth(70)
            
            for i in range(len(self.file.getZnucl())):
                
                self.CBox1.addItem(str(self.PTOE.getName(self.file.getZnucl()[i])))
                
        elif ntypat >= 2:
            
            self.lbl1 = QtGui.QLabel(" Atom type 1 :", self)
            self.lbl1.setFixedWidth(85)
            
            self.lbl2 = QtGui.QLabel(" Atom type 2 :", self)
            self.lbl2.setFixedWidth(85)

            self.CBox1 = QtGui.QComboBox()
            self.CBox1.setFixedWidth(70)
            self.CBox2 = QtGui.QComboBox()
            self.CBox2.setFixedWidth(70)
            
            self.CBox1.addItem(str('all'))
            self.CBox2.addItem(str('all'))
            
            for i in range(len(self.file.getZnucl())):

                self.CBox1.addItem(str(self.PTOE.getName(self.file.getZnucl()[i])))
                self.CBox2.addItem(str(self.PTOE.getName(self.file.getZnucl()[i])))



        self.lbl3 = QtGui.QLabel("delta R :", self)
        self.lbl3.setFixedWidth(55)

        self.sbdr = QtGui.QDoubleSpinBox()
        self.sbdr.setSingleStep(0.01)
        self.sbdr.setValue(0.1)
        self.sbdr.setMinimum(0.01)
        self.sbdr.setMaximum(2)
        self.sbdr.setFixedWidth(70)

        self.lbl4 = QtGui.QLabel("Box wide :", self)
        self.lbl4.setFixedWidth(68)

        self.sbrm = QtGui.QDoubleSpinBox()
        self.sbrm.setSingleStep(0.5)
        self.sbrm.setMinimum(-0.5)
        self.sbrm.setMaximum(2.0)
        self.sbrm.setValue(0)
        self.sbrm.setFixedWidth(70)

        self.lbl5 = QtGui.QLabel("Step :", self)
        self.lbl5.setFixedWidth(36)
        self.sbStep = QtGui.QSpinBox()
        self.sbStep.setMinimum(1)
        self.sbStep.setMaximum(self.file.getNbTime())
        self.sbStep.setFixedWidth(70)
        self.sbStep.setValue(5)

        self.lbl8 = QtGui.QLabel("Deconvolution :", self)
        self.lbl8.setFixedWidth(105)
        self.sbdec = QtGui.QCheckBox()

        self.pbClose = QtGui.QPushButton("close")
        self.pbClose.setFixedSize(70,20)
        self.connect(self.pbClose,QtCore.SIGNAL("clicked()"),QtCore.SLOT('close()'))

        self.pbDraw = QtGui.QPushButton("g(r)")
        self.pbDraw.setFixedSize(70,20)
        self.connect(self.pbDraw,QtCore.SIGNAL("clicked()"),self.Graphics)

        self.changeBox()

    def changeBox(self):

        mode = self.sbm.currentText()

        typat = self.file.getTypat()
        ntypat = max(typat)
        nat = 0

        if ntypat == 1:

            s1 = 0
            
        elif ntypat >= 2:
        
            s1 = 50

        if mode == 'RDF':

            self.setFixedSize(200, 250+s1)

            if ntypat == 1:

                self.layout.addWidget(self.mode    , 1, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbm     , 1, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl1    , 2, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.CBox1   , 2, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl3    , 3, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbdr    , 3, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl4    , 4, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbrm    , 4, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl5    , 5, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbStep  , 5, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl8    , 6, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbdec   , 6, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.pbClose , 8, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.pbDraw  , 8, 1, 1, 1, QtCore.Qt.AlignCenter)
            
            else:
            
                self.layout.addWidget(self.mode    , 1, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbm     , 1, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl1    , 2, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.CBox1   , 2, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl2    , 3, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.CBox2   , 3, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl3    , 4, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbdr    , 4, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl4    , 5, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbrm    , 5, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl5    , 6, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbStep  , 6, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl8    , 7, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbdec   , 7, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.pbClose , 9, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.pbDraw  , 9, 1, 1, 1, QtCore.Qt.AlignCenter)

        elif mode == 'ADF':

            self.setFixedSize(200, 200+s1)

            if ntypat == 1:

                self.layout.addWidget(self.mode    , 1, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbm     , 1, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl1    , 2, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.CBox1   , 2, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl6    , 3, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbdt    , 3, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl7    , 4, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbn     , 4, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl5    , 5, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbStep  , 5, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.pbClose , 7, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.pbDraw  , 7, 1, 1, 1, QtCore.Qt.AlignCenter)
            
            else:
            
                self.layout.addWidget(self.mode    , 1, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbm     , 1, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl1    , 2, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.CBox1   , 2, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl2    , 3, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.CBox2   , 3, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl6    , 4, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbdt    , 4, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl7    , 5, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbn     , 5, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl5    , 6, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbStep  , 6, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.pbClose , 8, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.pbDraw  , 8, 1, 1, 1, QtCore.Qt.AlignCenter)


        elif mode == 'NDF':

            self.setFixedSize(200, 200+s1)

            if ntypat == 1:

                self.layout.addWidget(self.mode    , 1, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbm     , 1, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl1    , 2, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.CBox1   , 2, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl3    , 3, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbdr    , 3, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl7    , 4, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbn     , 4, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl5    , 5, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbStep  , 5, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.pbClose , 7, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.pbDraw  , 7, 1, 1, 1, QtCore.Qt.AlignCenter)

            else:

                self.layout.addWidget(self.mode    , 1, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbm     , 1, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl1    , 2, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.CBox1   , 2, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl2    , 3, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.CBox2   , 3, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl3    , 4, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbdr    , 4, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl7    , 5, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbn     , 5, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.lbl5    , 6, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.sbStep  , 6, 1, 1, 1, QtCore.Qt.AlignCenter)
                self.layout.addWidget(self.pbClose , 8, 0, 1, 1, QtCore.Qt.AlignRight)
                self.layout.addWidget(self.pbDraw  , 8, 1, 1, 1, QtCore.Qt.AlignCenter)

        elif mode == 'Proba':

            self.setFixedSize(200, 170)

            self.layout.addWidget(self.mode    , 1, 0, 1, 1, QtCore.Qt.AlignRight)
            self.layout.addWidget(self.sbm     , 1, 1, 1, 1, QtCore.Qt.AlignCenter)
            self.layout.addWidget(self.lbl1    , 2, 0, 1, 1, QtCore.Qt.AlignRight)
            self.layout.addWidget(self.CBox1   , 2, 1, 1, 1, QtCore.Qt.AlignCenter)
            self.layout.addWidget(self.lbl7    , 3, 0, 1, 1, QtCore.Qt.AlignRight)
            self.layout.addWidget(self.sbn     , 3, 1, 1, 1, QtCore.Qt.AlignCenter)
            self.layout.addWidget(self.pbClose , 5, 0, 1, 1, QtCore.Qt.AlignRight)
            self.layout.addWidget(self.pbDraw  , 5, 1, 1, 1, QtCore.Qt.AlignCenter)

        self.show()

	self.connect(self.sbm,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.changeMode)

    def changeMode(self):

        mode = self.sbm.currentText()

        # erase text on the window

        try:
            self.layout.removeWidget(self.lbl1)
        except:
            pass
        try:
            self.layout.removeWidget(self.lbl2)
        except:
            pass
        try:
            self.layout.removeWidget(self.lbl3)
        except:
            pass
        try:
            self.layout.removeWidget(self.lbl4)
        except:
            pass
        try:
            self.layout.removeWidget(self.lbl5)
        except:
            pass
        try:
            self.layout.removeWidget(self.lbl6)
        except:
            pass
        try:
            self.layout.removeWidget(self.lbl7)
        except:
            pass
        try:
            self.layout.removeWidget(self.lbl8)
        except:
            pass
        try:
            self.lbl1.setParent(None)
        except:
            pass
        try:
            self.lbl2.setParent(None)
        except:
            pass
        try:
            self.lbl3.setParent(None)
        except:
            pass
        try:
            self.lbl4.setParent(None)
        except:
            pass
        try:
            self.lbl5.setParent(None)
        except:
            pass
        try:
            self.lbl6.setParent(None)
        except:
            pass
        try:
            self.lbl7.setParent(None)
        except:
            pass
        try:
            self.lbl8.setParent(None)
        except:
            pass

        try:
            self.layout.removeWidget(self.CBox1)
        except:
            pass
        try:
            self.layout.removeWidget(self.CBox2)
        except:
            pass            
        try:
            self.layout.removeWidget(sbdr)
        except:
            pass
        try:
            self.layout.removeWidget(sbrm)
        except:
            pass
        try:
            self.layout.removeWidget(sbstep)
        except:
            pass
        try:
            self.layout.removeWidget(sbdt)
        except:
            pass
        try:
            self.layout.removeWidget(sbn)
        except:
            pass
        try:
            self.layout.removeWidget(sdec)
        except:
            pass
        try:
            self.layout.removeWidget(pbClose)
        except:
            pass
        try:
            self.layout.removeWidget(pbDraw)
        except:
            pass
        try:
            self.CBox1.setParent(None)
        except:
            pass
        try:
            self.CBox2.setParent(None)
        except:
            pass
        try:
            self.sbdr.setParent(None)
        except:
            pass
        try:
            self.sbdm.setParent(None)
        except:
            pass
        try:
            self.sbstep.setParent(None)
        except:
            pass
        try:
            self.sbdt.setParent(None)
        except:
            pass
        try:
            self.sbn.setParent(None)
        except:
            pass
        try:
            self.sbdec.setParent(None)
        except:
            pass
        try:
            self.pbClose.setParent(None)
        except:
            pass
        try:
            self.pbDraw.setParent(None)
        except:
            pass


        typat = self.file.getTypat()
        ntypat = max(typat)
        nat = 0

        for i in enumerate(typat):
            nat = nat+1

        #initialisation of texts

        if mode == 'RDF':

            if ntypat == 1:

                self.lbl1 = QtGui.QLabel(" Atom type :", self)
                self.lbl1.setFixedWidth(80)
            
                self.CBox1 = QtGui.QComboBox()
                self.CBox1.setFixedWidth(70)
            
                for i in range(len(self.file.getZnucl())):
                
                    self.CBox1.addItem(str(self.PTOE.getName(self.file.getZnucl()[i])))

            elif ntypat >= 2:

                self.lbl1 = QtGui.QLabel(" Atom type 1 :", self)
                self.lbl1.setFixedWidth(85)
                
                self.lbl2 = QtGui.QLabel(" Atom type 2 :", self)
                self.lbl2.setFixedWidth(85)

                self.CBox1 = QtGui.QComboBox()
                self.CBox1.setFixedWidth(70)
                self.CBox2 = QtGui.QComboBox()
                self.CBox2.setFixedWidth(70)
            
                self.CBox1.addItem(str('all'))
                self.CBox2.addItem(str('all'))
            
                for i in range(len(self.file.getZnucl())):
                    
                    self.CBox1.addItem(str(self.PTOE.getName(self.file.getZnucl()[i])))
                    self.CBox2.addItem(str(self.PTOE.getName(self.file.getZnucl()[i])))

            self.lbl3 = QtGui.QLabel("delta R :", self)
            self.lbl3.setFixedWidth(55)
            self.sbdr = QtGui.QDoubleSpinBox()
            self.sbdr.setSingleStep(0.01)
            self.sbdr.setValue(0.1)
            self.sbdr.setMinimum(0.01)
            self.sbdr.setMaximum(2)
            self.sbdr.setFixedWidth(70)

            self.lbl4 = QtGui.QLabel("Box wide :", self)
            self.lbl4.setFixedWidth(68)
            self.sbrm = QtGui.QDoubleSpinBox()
            self.sbrm.setSingleStep(0.5)
            self.sbrm.setMinimum(-0.5)
            self.sbrm.setMaximum(2.0)
            self.sbrm.setValue(0)
            self.sbrm.setFixedWidth(70)

            self.lbl5 = QtGui.QLabel("Step :", self)
            self.lbl5.setFixedWidth(36)
            self.sbStep = QtGui.QSpinBox()
            self.sbStep.setMinimum(1)
            self.sbStep.setMaximum(self.file.getNbTime())
            self.sbStep.setFixedWidth(70)
            self.sbStep.setValue(5)

            self.lbl8 = QtGui.QLabel("Deconvolution :", self)
            self.lbl8.setFixedWidth(105)
            self.sbdec = QtGui.QCheckBox()

        elif mode == 'ADF':

            if ntypat == 1:

                self.lbl1 = QtGui.QLabel(" Atom type :", self)
                self.lbl1.setFixedWidth(80)
            
                self.CBox1 = QtGui.QComboBox()
                self.CBox1.setFixedWidth(70)
            
                for i in range(len(self.file.getZnucl())):
                
                    self.CBox1.addItem(str(self.PTOE.getName(self.file.getZnucl()[i])))

            elif ntypat >= 2:

                self.lbl1 = QtGui.QLabel(" Atom type 1 :", self)
                self.lbl1.setFixedWidth(85)
                
                self.lbl2 = QtGui.QLabel(" Atom type 2 :", self)
                self.lbl2.setFixedWidth(85)

                self.CBox1 = QtGui.QComboBox()
                self.CBox1.setFixedWidth(70)
                self.CBox2 = QtGui.QComboBox()
                self.CBox2.setFixedWidth(70)
            
                self.CBox1.addItem(str('all'))
                self.CBox2.addItem(str('all'))
            
                for i in range(len(self.file.getZnucl())):

                    self.CBox1.addItem(str(self.PTOE.getName(self.file.getZnucl()[i])))
                    self.CBox2.addItem(str(self.PTOE.getName(self.file.getZnucl()[i])))

            self.lbl5 = QtGui.QLabel("Step :", self)
            self.lbl5.setFixedWidth(36)
            self.sbStep = QtGui.QSpinBox()
            self.sbStep.setMinimum(1)
            self.sbStep.setMaximum(self.file.getNbTime())
            self.sbStep.setFixedWidth(70)
            self.sbStep.setValue(5)

            self.lbl6 = QtGui.QLabel("delta theta :", self)
            self.lbl6.setFixedWidth(80)
            self.sbdt = QtGui.QDoubleSpinBox()
            self.sbdt.setSingleStep(0.2)
            self.sbdt.setValue(1)
            self.sbdt.setMinimum(0.2)
            self.sbdt.setMaximum(5)
            self.sbdt.setFixedWidth(70)

            self.lbl7 = QtGui.QLabel("cutoff :", self)
            self.lbl7.setFixedWidth(60)
            self.sbn = QtGui.QDoubleSpinBox()
            self.sbn.setSingleStep(0.1)
            self.sbn.setValue(4)
            self.sbn.setMinimum(1)
            self.sbn.setMaximum(20)
            self.sbn.setFixedWidth(70)

        elif mode == 'NDF':

            if ntypat == 1:

                self.lbl1 = QtGui.QLabel(" Atom type :", self)
                self.lbl1.setFixedWidth(80)
            
                self.CBox1 = QtGui.QComboBox()
                self.CBox1.setFixedWidth(70)
            
                for i in range(len(self.file.getZnucl())):
                
                    self.CBox1.addItem(str(self.PTOE.getName(self.file.getZnucl()[i])))

            elif ntypat >= 2:

                self.lbl1 = QtGui.QLabel(" Atom type 1 :", self)
                self.lbl1.setFixedWidth(85)
            
                self.lbl2 = QtGui.QLabel(" Atom type 2 :", self)
                self.lbl2.setFixedWidth(85)

                self.CBox1 = QtGui.QComboBox()
                self.CBox1.setFixedWidth(70)
                self.CBox2 = QtGui.QComboBox()
                self.CBox2.setFixedWidth(70)
                
                self.CBox1.addItem(str('all'))
                self.CBox2.addItem(str('all'))
                
                for i in range(len(self.file.getZnucl())):

                    self.CBox1.addItem(str(self.PTOE.getName(self.file.getZnucl()[i])))
                    self.CBox2.addItem(str(self.PTOE.getName(self.file.getZnucl()[i])))

            self.lbl3 = QtGui.QLabel("delta R :", self)
            self.lbl3.setFixedWidth(55)
            self.sbdr = QtGui.QDoubleSpinBox()
            self.sbdr.setSingleStep(0.01)
            self.sbdr.setValue(0.05)
            self.sbdr.setMinimum(0.01)
            self.sbdr.setMaximum(1)
            self.sbdr.setFixedWidth(70)

            self.lbl5 = QtGui.QLabel("Step :", self)
            self.lbl5.setFixedWidth(36)
            self.sbStep = QtGui.QSpinBox()
            self.sbStep.setMinimum(1)
            self.sbStep.setMaximum(self.file.getNbTime())
            self.sbStep.setFixedWidth(70)
            self.sbStep.setValue(5)

            self.lbl7 = QtGui.QLabel("neighbors :", self)
            self.lbl7.setFixedWidth(71)
            self.sbn = QtGui.QSpinBox()
            self.sbn.setSingleStep(1)
            self.sbn.setValue(1)
            self.sbn.setMinimum(1)
            self.sbn.setMaximum(nat-1)
            self.sbn.setFixedWidth(70)

        elif mode == 'Proba':

            self.lbl7 = QtGui.QLabel("neighbors :", self)
            self.lbl7.setFixedWidth(71)
            self.sbn = QtGui.QSpinBox()
            self.sbn.setSingleStep(1)
            self.sbn.setValue(1)
            self.sbn.setMinimum(1)
            self.sbn.setMaximum(nat-1)
            self.sbn.setFixedWidth(70)

        self.pbClose = QtGui.QPushButton("close")
        self.pbClose.setFixedSize(70,20)
        self.connect(self.pbClose,QtCore.SIGNAL("clicked()"),self.close)

        self.pbDraw = QtGui.QPushButton("g(r)")
        self.pbDraw.setFixedSize(70,20)
        self.connect(self.pbDraw,QtCore.SIGNAL("clicked()"),self.Graphics)

        self.changeBox()

        #------------------------------------------------------------------------#

    def Graphics(self):
        
        mode1 = self.sbm.currentText()

        if mode1 == 'RDF':

            typat = self.file.getTypat()
            ntypat = max(typat)

            if ntypat == 1:
                
                atom1 = 1
                atom2 = 1

            elif ntypat >= 2:

                atom1 = self.CBox1.currentIndex()
                atom2 = self.CBox2.currentIndex()
                
            box  = self.sbrm.value()
            dr    = self.sbdr.value()
            step  = self.sbStep.value()
                

            if self.sbdec.isChecked():
                
                self.RadialDistrib = Analysis.RDF(self.file,1,atom1,atom2,box,dr,step)
                
                data = self.RadialDistrib.getDATA()

                x = self.RadialDistrib.getR()

                for i in range(1,len(typat)-1):

                    self.Deconvolution = Analysis.DEC(self.file,i,atom1,atom2,data,box,dr,step)
                
                    y = self.Deconvolution.getNEI()

                    try:
                        self.GraphDEC.addGraph(x,y,False)
                    except:
                        self.GraphDEC = Graph.graphic(x,y,'R (Bohr)', "Radial Distribution", average=False,name = self.name)
                    
                self.GraphDEC.show()

            else:

                self.RadialDistrib = Analysis.RDF(self.file,0,atom1,atom2,box,dr,step)
                
                x = self.RadialDistrib.getR()
                
                y = self.RadialDistrib.getRDF()
        
                try:
                    self.GraphRDF.update(x,y,'R (Bohr)', "Radial Distribution",name = self.name)
                    self.GraphRDF.addPlot(x,linspace(1,1,len(x)))
                except:
                    self.GraphRDF = Graph.graphic(x,y,'R (Bohr)', "Radial Distribution", average=False,name = self.name)
                    self.GraphRDF.addPlot(x,linspace(1,1,len(x)))
            #self.connect(self.GraphRDF, QtCore.SIGNAL("myCustomizedSignal()"), self.close)
                self.GraphRDF.show()

                def close(self):
                    self.hide()
                    try:
                        del self.GraphRDF
                    except:
                        pass

                y = self.RadialDistrib.getINT()
            
                try:
                    self.GraphINT.update(x,y,'R (Bohr)', "Integral of the Radial Distribution",name = self.name)
                except:
                    self.GraphINT = Graph.graphic(x,y,'R (Bohr)', "Integral of the Radial Distribution", average=False,name = self.name)
            #self.connect(self.GraphRDF, QtCore.SIGNAL("myCustomizedSignal()"), self.close)
                self.GraphINT.show()
    
                def close(self):
                    self.hide()
                    try:
                        del self.GraphINT
                    except:
                        pass

        elif mode1 == 'ADF':

            neib = self.sbn.value()

            typat = self.file.getTypat()
            ntypat = max(typat)
        
            if ntypat == 1:

                atom1 = 1
                atom2 = 1

            elif ntypat >= 2:

                atom1 = self.CBox1.currentIndex()
                atom2 = self.CBox2.currentIndex()

            dtheta = self.sbdt.value()
            step  = self.sbStep.value()

            self.AngularDistrib = Analysis.ADF(neib,self.file,atom1,atom2,dtheta,step)

            x = self.AngularDistrib.getTheta()

            y = self.AngularDistrib.getADF()
        
            try:
                self.GraphADF.update(x,y,'Theta (degres)', "Angular Distribution",name = self.name)
            except:
                self.GraphADF = Graph.graphic(x,y,'Theta (degres)', "Angular Distribution", average=False,name = self.name)
            #self.connect(self.GraphRDF, QtCore.SIGNAL("myCustomizedSignal()"), self.close)
            self.GraphADF.show()

            def close(self):
                self.hide()
                try:
                    del self.GraphADF                    
                except:
                    pass


        elif mode1 == 'NDF':


            neib = self.sbn.value()

            typat = self.file.getTypat()
            ntypat = max(typat)
        
            if ntypat == 1:

                atom1 = 1
                atom2 = 1

            elif ntypat >= 2:

                atom1 = self.CBox1.currentIndex()
                atom2 = self.CBox2.currentIndex()

            dr = self.sbdr.value()
            step  = self.sbStep.value()

            self.NeighborDistrib = Analysis.NDF(neib,self.file,atom1,atom2,dr,step)

            x = self.NeighborDistrib.getR()

            y = self.NeighborDistrib.getNDF()
        
            try:
                self.GraphNDF.update(x,y,'R (Bohr)', "Neighbor distance",name = self.name)
            except:
                self.GraphNDF = Graph.graphic(x,y,'R (Bohr)', "Neighbor distance", average=False,name = self.name)
            #self.connect(self.GraphRDF, QtCore.SIGNAL("myCustomizedSignal()"), self.close)
            self.GraphNDF.show()

            def close(self):
                self.hide()
                try:
                    del self.GraphNDF                    
                except:
                    pass


        elif mode1 == 'Proba':


            neib = self.sbn.value()

            typat = self.file.getTypat()
            ntypat = max(typat)
        
            if ntypat == 1:

                atom1 = 1

            elif ntypat >= 2:

                atom1 = self.CBox1.currentIndex()

            self.Probability = Analysis.Proba(neib,self.file,atom1)

            x = self.Probability.getN()

            y = self.Probability.getProba()
        
            try:
                self.GraphProba.update(x,y,'Neighbor', "Probability",name = self.name)
            except:
                self.GraphProba = Graph.graphic(x,y,'Neighbor', "Probability", average=False,name = self.name)
            #self.connect(self.GraphRDF, QtCore.SIGNAL("myCustomizedSignal()"), self.close)
            self.GraphProba.show()



    def center(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size =  self.geometry()
        self.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)


    def close(self):
        del self.graphPosition
        del self
    
    def closeEvent(self, event):
        try:
            del self.graphPosition
        except:
            pass
        try:
            del self
        except:
            pass     

