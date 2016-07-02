#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import sys,os
import string
try:
    from PyQt4 import Qt,QtGui,QtCore
except:
    pass

import utility.global_variable as var

class About(QtGui.QWidget):

    def __init__(self, parent = None):
        self.name = "About"
        self.initUI(parent)
      	self.raise_()

    def initUI(self, parent):
        QtGui.QWidget.__init__(self, parent)
        self.setWindowTitle(self.name)
        self.setFixedSize(450,350)
        self.center()
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)

        self.lbl = QtGui.QLabel("   APPA", self)
        self.lbl.setFixedWidth(140)
        self.lbl.setFont(QtGui.QFont("calibri", 25))

        self.lblversion = QtGui.QLabel(var.version, self)
        self.lblversion.setFont(QtGui.QFont("calibri", 10))

        self.tbro = QtGui.QTextBrowser()
        self.tbro.setFixedHeight(250)
        self.tbro.setFixedWidth(430)
        self.tbro.setText(' 05/20/2014 -- Version 1.0.8'+'\n'
                          +'-Position display improved (type of atoms)\n'
                          +'-Save simulation in NetCDf Format\n\n'
                          +'01/17/2014 -- Version 1.0.7'+'\n'\
                          +'-periodic boundary conditions are better taken into account\n'
                          +'-RDF can now be calculated further the cell limit\n'\
                          +'-the RDF is corrected\n'\
                          +'-add integral of the RDF\n'\
                          +'-add angular distribution function (ADF)\n'\
                          +'-add neibhor distribution function (NDF)\n'\
                          +'-add deconvolution of the RDF\n'\
                          +'-add survival probability (BETA)\n'\
                          +'-add graphic options for axis view : XY XZ or YZ\n'\
                          +'-add graphic options for the type of view : standard, PBC or average\n'\
                          +'-graphic now can show two kind of atoms\n\n'\
                          +'04/29/2013 -- Version 1.0.6 '+'\n'\
                          +'-able to read PIMD output files\n'\
                          +'-adding graphic of acell\n'\
                          +'-adding graphic of angles\n'\
                          +'-adding graphic of volume\n'\
                          +'-able to read NPT MD (G(r), MSD, VDOS, VACF not available)\n'\
                          +'-mean squared displacement calculation\n'\
                          +'\n03/18/2013 -- Version 1.0.5'+'\n'\
                          +'-no windows version can read input file\n'\
                          +'-add fotran calculation for g(r)\n'\
                          +'-add fotran calculation for average\n'\
                          +'-add fotran calculation for standard deviation\n'\
                          +'\n02/12/2013 -- Version 1.0.4'+'\n'\
                          +'-Corecting bug relating to the g(r)\n'\
                          +'\n02/12/2013 -- Version 1.0.3'+'\n'\
                          +'-Reading file use threading\n'\
                          +'-Change the loading bar when reading file\n'\
                          +'-improvement calculation of kinetic energy\n'\
                          +'-improvement reading of the output file\n'\
                          +'-Correction of the bug of the vibrational density of states\n'\
                          +'\n01/21/2013 -- Version 1.0.2'+'\n'\
                          +'-Removing question when closing windows\n'\
                          +'-Change the reading of classic output file\n'\
                          +'-Correcting bug relating to the G(r)\n'\
                          +'-Correcting bug relating to reading of the stress\n'\
                          +'\n01/12/2013 -- Version 1.0.1'+'\n'\
                          +'- bugs correction\n'\
                          +'- \"ok\" bouton become \"update\" button\n'\
                          +'- Adding loading bar for RDF calculation\n'\
                          +'- Correcting several bugs with the loading bar for the no windows version\n'\
                          +'- Adding the About windows and menu\n'\
                          +'- Adding visualization of the particles in 2D via \"Positions\" button\n'\
                          +'- Adding distance units (bohr and angstrom)\n'\
                          +'\n01/07/2013 -- Version 1.0.0 First Release')

        self.layout.addWidget(self.lbl    , 1, 0, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.lblversion, 2, 0, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.tbro   , 3, 0, 10, 1, QtCore.Qt.AlignCenter)
        self.show()



    def center(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size =  self.geometry()
        self.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)


