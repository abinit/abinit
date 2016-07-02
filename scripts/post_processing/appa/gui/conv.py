#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import sys,os
import string, math
from PyQt4 import Qt,QtGui,QtCore
from numpy import pi

class Conversion(QtGui.QWidget):
	
    def __init__(self, parent = None):
	    	
        self.initUI(parent)

        #Unit by default (from the output files):
        self.default      = {'Pressure':['GPa',1],'Energy':['Ha',1],'Temperature':['K',0],'Volume':["Bohr^3",1],\
                                 'Distance':["Bohr",1],'Angle':["Degree",1]}
	
	      #Unit Change:
        self.units_change = {'Pressure':['GPa',1],'Energy':["Ha",1],'Temperature':['K',0],'Volume':["Bohr^3",1],\
                                 'Distance':["Bohr",1],'Angle':["Degree",1]}
	
        self.setBox(self.default)	    
	    
    
    def initUI(self, parent):
	    
	#-----------------Creation of the windows----------------------------#
        QtGui.QWidget.__init__(self, parent)
        self.setWindowTitle('Units')
        self.setFixedSize(260, 250)
        self.center()
	
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
	
	
        self.pbapply = QtGui.QPushButton("apply")
        self.pbapply.setFixedSize(70,20)
        self.pbapply.setEnabled(False)
        self.connect(self.pbapply,QtCore.SIGNAL("clicked()"),self.apply)
        
        self.pbcancel = QtGui.QPushButton("cancel")
        self.pbcancel.setFixedSize(70,20)
        self.connect(self.pbcancel,QtCore.SIGNAL("clicked()"),self.cancel)
        
        self.pbok = QtGui.QPushButton("ok")
        self.pbok.setFixedSize(70,20)
        self.connect(self.pbok,QtCore.SIGNAL("clicked()"),self.ok)
        
        self.lble     = QtGui.QLabel('energy :')
        self.lble.setFixedSize(100,20)
        
        self.lblDistance  = QtGui.QLabel('Distance :')
        self.lblDistance.setFixedSize(100,20)
        
        self.lblvol   = QtGui.QLabel('Volume :')
        self.lblvol.setFixedSize(100,20)
	
        self.lblpress = QtGui.QLabel('Pressure :')
        self.lblpress.setFixedSize(100,20)
	
        self.lbltemp  = QtGui.QLabel('Temperature :')
        self.lbltemp.setFixedSize(100,20)

        self.lblangle  = QtGui.QLabel('Angle :')
        self.lblangle.setFixedSize(100,20)

	
        self.CBoxe = QtGui.QComboBox()
        self.CBoxe.addItem("Ha")
        self.CBoxe.addItem("eV")
        self.connect(self.CBoxe,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.boxChange)
        self.CBoxe.setFixedSize(100,20)
	
        self.CBoxvol = QtGui.QComboBox()
        self.CBoxvol.addItem("Bohr^3")
        self.CBoxvol.addItem("Angstrom^3")
        self.connect(self.CBoxvol,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.boxChange)
        self.CBoxvol.setFixedSize(100,20)
	
	
        self.CBoxDistance = QtGui.QComboBox()
        self.CBoxDistance.addItem("Bohr")
        self.CBoxDistance.addItem("Angstrom")
        self.connect(self.CBoxDistance,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.boxChange)
        self.CBoxDistance.setFixedSize(100,20)

        self.CBoxpress = QtGui.QComboBox()
        self.CBoxpress.addItem("Ha/bohr^3")
        self.CBoxpress.addItem("GPa")
        self.connect(self.CBoxpress,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.boxChange)
        self.CBoxpress.setFixedSize(100,20)	
        
        self.CBoxtemp = QtGui.QComboBox()
        self.CBoxtemp.addItem("K")
        self.CBoxtemp.addItem("C")
        self.connect(self.CBoxtemp,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.boxChange)		
        self.CBoxtemp.setFixedSize(100,20)	

        self.CBoxAngle = QtGui.QComboBox()
        self.CBoxAngle.addItem("Degree")
        self.CBoxAngle.addItem("Radian")
        self.connect(self.CBoxAngle,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.boxChange)		
        self.CBoxAngle.setFixedSize(100,20)	

        
        
        self.layout.addWidget(self.lble     , 1, 0, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.CBoxe    , 1, 1, 1, 2, QtCore.Qt.AlignCenter)
        
        self.layout.addWidget(self.lblvol   , 2, 0, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.CBoxvol  , 2, 1, 1, 2, QtCore.Qt.AlignCenter)

        self.layout.addWidget(self.lblDistance  , 3, 0, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.CBoxDistance , 3, 1, 1, 2, QtCore.Qt.AlignCenter)
	
        self.layout.addWidget(self.lblpress , 4, 0, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.CBoxpress, 4, 1, 1, 2, QtCore.Qt.AlignCenter)
        
        self.layout.addWidget(self.lbltemp  , 5, 0, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.CBoxtemp , 5, 1, 1, 2, QtCore.Qt.AlignCenter)
        
        self.layout.addWidget(self.lblangle  , 6, 0, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.CBoxAngle , 6, 1, 1, 2, QtCore.Qt.AlignCenter)
        
        self.layout.addWidget(self.pbcancel , 8, 0, 1, 1, QtCore.Qt.AlignCenter)	
        self.layout.addWidget(self.pbapply  , 8, 1, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.pbok     , 8, 2, 1, 1, QtCore.Qt.AlignCenter)
	#--------------------------------------------------------------------------#
    
    
    def boxChange(self,unit):
        self.pbapply.setEnabled(True)
		
    def setBox(self,units):
        self.CBoxe.setCurrentIndex(self.CBoxe.findText(units['Energy'][0]))	
        self.CBoxvol.setCurrentIndex(self.CBoxvol.findText(units['Volume'][0]))
        self.CBoxpress.setCurrentIndex(self.CBoxpress.findText(units['Pressure'][0]))
        self.CBoxtemp.setCurrentIndex(self.CBoxtemp.findText(units['Temperature'][0]))
        self.CBoxDistance.setCurrentIndex(self.CBoxDistance.findText(units['Distance'][0]))
        self.CBoxAngle.setCurrentIndex(self.CBoxAngle.findText(units['Angle'][0]))

    def ok(self):
        if self.pbapply.isEnabled():
            self.changeUnits()
            self.emit(QtCore.SIGNAL("myCustomizedSignal(PyQt_PyObject)"),self.units_change)
            self.hide()
        else:
            self.hide()
    def apply(self):
        self.changeUnits()
        self.emit(QtCore.SIGNAL("myCustomizedSignal(PyQt_PyObject)"),self.units_change)	
        self.pbapply.setEnabled(False)
	
	
    def cancel(self):
        self.hide()

    def changeUnits(self):	    
	      #-----------------Conversion value (relative to the self.default)--------------:
        self.conv = {'eV':27.2113834,'C':273.15,"Ha/bohr^3":1/29421.033e0,"Angstrom^3":0.5291772085936**3,\
                         "Angstrom":0.5291772085936,"Radian":pi/180}
	
        self.units_change['Energy'][0]      = str(self.CBoxe.currentText())
        self.units_change['Volume'][0]      = str(self.CBoxvol.currentText())
        self.units_change['Pressure'][0]    = str(self.CBoxpress.currentText())
        self.units_change['Temperature'][0] = str(self.CBoxtemp.currentText())
        self.units_change['Distance'][0]    = str(self.CBoxDistance.currentText())
        self.units_change['Angle'][0]       = str(self.CBoxAngle.currentText())
	
        #---------------Change the energy value--------------------#
        if self.units_change['Energy'][0] != self.default['Energy'][0]:
            self.units_change['Energy'][1] = self.conv[self.units_change['Energy'][0]]	
        else :
            self.units_change['Energy'][1] = self.default['Energy'][1]	
        #-----------------------------------------------------------#	
	    
	      #---------------Change the temperature value-----------------#   
        if self.units_change['Temperature'][0] != self.default['Temperature'][0]:
            self.units_change['Temperature'][1] = self.conv[self.units_change['Temperature'][0]]	
        else :
            self.units_change['Temperature'][1] = self.default['Temperature'][1]		
 	      #-----------------------------------------------------------#
	
	      #---------------Change the Pressure value-----------------#   
        if self.units_change['Pressure'][0] != self.default['Pressure'][0]:
            self.units_change['Pressure'][1] = self.conv[self.units_change['Pressure'][0]]	
        else :
            self.units_change['Pressure'][1] = self.default['Pressure'][1]		
        #-----------------------------------------------------------#
	
	      #---------------Change the Volume value-----------------#   
        if self.units_change['Volume'][0] != self.default['Volume'][0]:
            self.units_change['Volume'][1] = self.conv[self.units_change['Volume'][0]]	
        else :
            self.units_change['Volume'][1] = self.default['Volume'][1]		
	      #-----------------------------------------------------------#
	
	      #---------------Change the Distancess value-----------------#   
        if self.units_change['Distance'][0] != self.default['Distance'][0]:
            self.units_change['Distance'][1] = self.conv[self.units_change['Distance'][0]]	
        else :
            self.units_change['Distance'][1] = self.default['Distance'][1]		
	      #-----------------------------------------------------------#
	
	      #---------------Change the Angles value-----------------#   
        if self.units_change['Angle'][0] != self.default['Angle'][0]:
            self.units_change['Angle'][1] = self.conv[self.units_change['Angle'][0]]	
        else :
            self.units_change['Angle'][1] = self.default['Angle'][1]		
	      #-----------------------------------------------------------#
	
    def showUnits(self):
        self.setBox(self.units_change)    
        self.pbapply.setEnabled(False)	
        self.show()
	
    def getUnits(self):
        return self.units_change
    
    def center(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size =  self.geometry()
        self.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)
        
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
