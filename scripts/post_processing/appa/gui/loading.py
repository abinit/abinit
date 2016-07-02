#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import os, sys, time
import string, math, re

try:
    from PyQt4 import Qt,QtGui,QtCore
except:
    pass


class loading(QtGui.QWidget):

#-----------------------------#
#-------CONSTRUCTOR-----------#
#-----------------------------#

    def __init__(self, parent = None, message = "please wait"):
            self.message = message
            self.initUI(parent)           
            self.raise_()
            
    def initUI(self, parent):

	#-----------------Creation of the windows----------------------------#
        QtGui.QWidget.__init__(self, parent)
        self.setWindowTitle('Please wait')
        self.setFixedSize(160, 60)
        self.center()

        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)

        self.lblload = QtGui.QLabel(self.message)
        
        self.progressbar = QtGui.QProgressBar(self)
        self.progressbar.setMinimum(0)
        self.progressbar.setMaximum(0)
        self.progressbar.setFixedSize(140,20)

        self.layout.addWidget(self.lblload    , 1, 0, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.progressbar, 2, 0, 1, 1, QtCore.Qt.AlignCenter)
        
        self.show()

    def cancel(self):
        self.thread.stop()
        return

    def update(self,pmessage):
        self.lblload.setText(pmessage)
        self.raise_()

    def center(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size =  self.geometry()
        self.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)

    def __del__(self):
        pass;
 

    def closeEvent(self, event):
        self.raise_()
        event.ignore()
