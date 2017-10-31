#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
authors: Martin Alexandre, Vincent Stutzmann
last edited: May 2014
"""

#GUI 
import gui.graph as Graph
import gui.conv  as Conv

#Utility
import utility.write as Write
import utility.positions as Atom
import utility.analysis as Analysis

try:
    from PyQt4 import Qt,QtGui,QtCore
except:
    pass;

from numpy import *

#------------------------------------------------------#
#-------------------WINDOWS POSITIONS------------------#
#------------------------------------------------------#

class WinPOS(QtGui.QWidget):

    def __init__(self, file,punits, parent = None,name =''):

        self.file = file
        self.name = name
        self.units = punits #Dictionary with the units ans the conversion
        self.x = []
        self.y = []
        self.legend_x = "" 
        self.legend_y = ""
        self.typat = self.file.getTypat()
        self.ntypat = max(self.typat)
        self.PTOE = Analysis.PeriodicTableElement()
        self.initUI(parent) #Creation of postion windows
        self.changeMode()   #Generation of graphic of position
        self.raise_()

    def initUI(self, parent):

        #-----------------Creation of the windows----------------------------#
        QtGui.QWidget.__init__(self, parent)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle('Position options')
        self.setFixedSize(200, 150)

        self.center()
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)

        self.mode1 = QtGui.QLabel("Position :", self)
        self.mode1.setFixedWidth(65)

        self.sbm1 = QtGui.QComboBox()
        self.sbm1.setFixedWidth(95)

        self.sbm1.addItem(str('Standard'))
        self.sbm1.addItem(str('PBC (cubic only)'))
        self.sbm1.addItem(str('Average'))
        self.connect(self.sbm1,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.changeMode)        

        self.mode2 = QtGui.QLabel("View :", self)
        self.mode2.setFixedWidth(45)

        self.sbm2 = QtGui.QComboBox()
        self.sbm2.setFixedWidth(70)

        self.sbm2.addItem(str('XY'))
        self.sbm2.addItem(str('XZ'))
        self.sbm2.addItem(str('YZ'))
        self.connect(self.sbm2,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.changeMode)

        self.lbl3 = QtGui.QLabel("Reduced :", self)
        self.lbl3.setFixedWidth(65)
        self.sbrel = QtGui.QCheckBox()
        self.connect(self.sbrel,QtCore.SIGNAL('clicked()'),self.changeMode)

        self.pbClose = QtGui.QPushButton("close")
        self.pbClose.setFixedSize(70,20)
        self.connect(self.pbClose,QtCore.SIGNAL("clicked()"),QtCore.SLOT('close()')) 

        self.layout.addWidget(self.mode1    , 1, 0, 1, 1, QtCore.Qt.AlignRight)
        self.layout.addWidget(self.sbm1     , 1, 1, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.mode2    , 2, 0, 1, 1, QtCore.Qt.AlignRight)
        self.layout.addWidget(self.sbm2     , 2, 1, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.lbl3     , 3, 0, 1, 1, QtCore.Qt.AlignRight)
        self.layout.addWidget(self.sbrel    , 3, 1, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.pbClose  , 7, 0, 1, 2, QtCore.Qt.AlignCenter)

        self.show()
        #------------------------------------------------------------------------#

    def changeMode(self):

        mode1 = self.sbm1.currentText()
        mode2 = self.sbm2.currentText()
        
        self.marker_size = 2 #Set maker_size to 2 (except average)
        self.markerscale = 5 #Set the relative size of legend marker vs original

        if self.sbrel.isChecked():
            reduced = True
            unit    = 1  # No Conversion 
            legend  = "" 
        else:
            reduced = False
            unit   = self.units['Distance'][1]
            legend =  " ("+str(self.units['Distance'][0])+")"
        option = 0

        if mode1 == 'Average':
            option = 1
            self.marker_size = 5 #set marker_size to 5 (better for visualization)
            self.markerscale = 2
        elif mode1 == 'PBC (cubic only)':
            option = 2 

        for itype in range(1,self.ntypat+1):
            
            indexAtom = array([i for i,x in enumerate(self.typat) if x == itype],dtype=int)

            self.pos_conv = Atom.convert_position(self.file,indexAtom,option=option,reduced=reduced)
            
            if mode2 == 'XY':
                self.x = self.pos_conv.getX()
                self.y = self.pos_conv.getY()
                self.legend_x = "x" + legend
                self.legend_y = "y" + legend

            elif mode2 == 'XZ':
                self.x = self.pos_conv.getX()
                self.y = self.pos_conv.getZ()
                self.legend_x = "x" + legend
                self.legend_y = "z" + legend

            else :
                self.x = self.pos_conv.getY()
                self.y = self.pos_conv.getZ()
                self.legend_x = "y" + legend
                self.legend_y = "z "+ legend
                
            if itype == 1:
                try: 
                    self.graphPosition.updatePos(self.x*unit,self.y*unit,\
                                                     self.legend_x,\
                                                     self.legend_y,\
                                                     marker_size=self.marker_size, name = self.name)
                except:
                    self.graphPosition = Graph.graphic(self.x*unit,self.y*unit,\
                                                           self.legend_x,\
                                                           self.legend_y,\
                                                           average=False,point = True,marker='.',\
                                                           marker_size=self.marker_size,name = self.name)
            else:
                self.graphPosition.addPlot(self.x*unit,self.y*unit,\
                                               point=True,marker_size=self.marker_size)

        self.legend = []
        for i in range(self.ntypat):
            self.legend.append(str(self.PTOE.getName(self.file.getZnucl()[i])))
            
        self.graphPosition.addLegend(self.legend,markerscale=self.markerscale)
        self.graphPosition.show()

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
