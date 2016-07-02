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


#utility
import utility.global_variable as var

#utility
import reading.read_graph as ReadGraph

class OpenGraph(QtGui.QDialog):

    def __init__(self, parent = None):
        self.name = "import graphics (BETA)"
        QtGui.QDialog.__init__(self,parent)
        self.initUI(parent)
      	self.raise_()

    def initUI(self, parent):
        QtGui.QWidget.__init__(self, parent)
        self.setWindowTitle(self.name)
        self.setFixedSize(600,400)
        self.center()
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)

        self.lbl2 = QtGui.QLabel("select file :", self)
        self.lbl2.setFixedSize(70,36)
        self.tefile = QtGui.QTextEdit()
        self.tefile.setReadOnly(True)
        self.tefile.setFixedSize(400,36)
 
        self.browse =  QtGui.QPushButton('&Browse', self)
        self.browse.setFixedSize(100,36)
        self.connect(self.browse ,QtCore.SIGNAL("clicked()"),self.openFile)	

        
        self.tbro = QtGui.QTextBrowser()
        self.tbro.setFixedHeight(300)
        self.tbro.setFixedWidth(530)
        self.tbro.setLineWrapMode(QtGui.QTextEdit.NoWrap)

        self.lbl3 = QtGui.QLabel("x :", self)
        self.lbl3.setFixedSize(20,36)

        self.CBox1 = QtGui.QComboBox()
        self.CBox1.addItem("")
        self.CBox1.setFixedSize(70,36)
        self.connect(self.CBox1,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.changeData)

        self.lbl4 = QtGui.QLabel("y :", self)
        self.lbl4.setFixedSize(20,36)

        self.CBox2= QtGui.QComboBox()
        self.CBox2.addItem("")
        self.CBox2.setFixedSize(70,36)
        self.connect(self.CBox2,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.changeData)


        self.imp =  QtGui.QPushButton('&import', self)
        self.imp.setFixedSize(100,36)
        self.connect(self.imp ,QtCore.SIGNAL("clicked()"),self.close)	


        self.layout.addWidget(self.lbl2   , 0, 0, 1, 1, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.tefile , 0, 1, 1, 4,QtCore.Qt.AlignLeft)
        self.layout.addWidget(self.browse , 0 ,5, 1, 1,QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.tbro   , 1, 0, 10, 6, QtCore.Qt.AlignCenter)
        self.layout.addWidget(self.lbl3   ,12, 0, 1, 1, QtCore.Qt.AlignRight)
        self.layout.addWidget(self.CBox1  ,12, 1, 1, 1, QtCore.Qt.AlignLeft)
        self.layout.addWidget(self.lbl4   ,12, 2, 1, 1, QtCore.Qt.AlignRight)
        self.layout.addWidget(self.CBox2  ,12, 3, 1, 1, QtCore.Qt.AlignLeft)
        self.layout.addWidget(self.imp    ,12, 5, 1, 1, QtCore.Qt.AlignLeft)
        self.show()

    def openFile(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file', var.path(), "FILE (*)")
        #Save the current path
        pathFile=str(fname)
        var.global_path = pathFile
        self.tefile.setText(pathFile)
        self.file = ReadGraph.ReadGraph(pathFile)
        self.tbro.setText(self.file.getFile())
        
        #add item in the comboxbox:
        nb = self.file.getNbColumn()
        if nb != 0:
            self.CBox1.clear()
            self.CBox2.clear()
            self.CBox1.addItem("")
            self.CBox2.addItem("")

            for i in range(nb):
                self.CBox1.addItem(str(i+1))
                self.CBox2.addItem(str(i+1))
        if nb == 2:
            self.CBox1.setCurrentIndex(1)
            self.CBox2.setCurrentIndex(2)
                
    def center(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size =  self.geometry()
        self.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)


    def changeData(self):
            if (self.CBox1.currentText()!='' and  self.CBox2.currentText()!=''):
                self.data = self.file.read(self.CBox1.currentText(),self.CBox2.currentText())
                
                
    def importGraph(self):
        return  self.data

        
