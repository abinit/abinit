#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import sys,os,commands,threading
import string, math, re

#GUI 
import gui.graph as Graph
import gui.md as MD
import gui.gs  as GS
import gui.elasticConstant as ElasticConstant
import gui.conv as Conv
import gui.about as About
import gui.loading as load

#Reading
import reading.read as Read

#Utility
import utility.write as Write
import utility.thread as thread
import utility.global_variable as var

from PyQt4 import Qt,QtGui,QtCore
import numpy as np


#---------------------------------------------------#
#---------------------------------------------------#
#----------------MAIN FRAME-------------------------#
#---------------------------------------------------#
#---------------------------------------------------#

class MainFrame(QtGui.QMainWindow):


    def __init__(self):
        super(MainFrame, self).__init__()
        
        self.conv = Conv.Conversion()
        self.connect(self.conv, QtCore.SIGNAL("myCustomizedSignal(PyQt_PyObject)"), self.changeUnits)
        self.units  = self.conv.getUnits()

        self.initUI()


    def initUI(self):
        #----------MainFrame parameters----------#
        self.setWindowTitle("APPA "+var.version)
        self.setFixedSize(700, 550)
        self.center()
        #-------------------------------------#

        #---------------Creation of menubar----------------------------------#
        self.open = QtGui.QAction( '&Open', self)
        self.open.setShortcut('Ctrl+O')
        self.open.setStatusTip('Open File')
        self.connect(self.open, QtCore.SIGNAL('triggered()'), self.showDialog)

        self.close = QtGui.QAction('&Exit', self)
        self.close.setShortcut('Ctrl+Q')
        self.close.setStatusTip('Exit application')
        self.connect(self.close, QtCore.SIGNAL('triggered()'), QtCore.SLOT('close()'))

        self.save = QtGui.QAction('&Save', self)
        self.save.setShortcut('Ctrl+S')
        self.save.setStatusTip('Save simulation data')
        self.connect(self.save, QtCore.SIGNAL('triggered()'), self.showSave)

        self.export = QtGui.QAction('E&xport (.xyz)', self)
        self.export.setShortcut('Ctrl+X')
        self.export.setStatusTip('Export data to XYZ file')
        self.connect(self.export, QtCore.SIGNAL('triggered()'), self.showExport)

        self.menubar = self.menuBar()
        self.fileMenu1 = self.menubar.addMenu('&File')
        self.fileMenu1.addAction(self.open)
        self.fileMenu1.addAction(self.save)
        self.fileMenu1.addAction(self.export)
        self.fileMenu1.addAction(self.close)

        self.ec = QtGui.QAction( '&Elastics constants', self)
        self.ec.setShortcut('Ctrl+E')
        self.ec.setStatusTip('Calculation of Elastics constants')
        self.connect(self.ec, QtCore.SIGNAL('triggered()'), self.showElastics)

        self.fileMenu2 = self.menubar.addMenu('&Calculation')
        self.fileMenu2.addAction(self.ec)

        self.unit = QtGui.QAction( '&Units', self)
        self.unit.setShortcut('Ctrl+U')
        self.unit.setStatusTip('change physical units')
        self.connect(self.unit, QtCore.SIGNAL('triggered()'), self.showConv)

        self.fileMenu3 = self.menubar.addMenu('&Option')
        self.fileMenu3.addAction(self.unit)

        self.about = QtGui.QAction( '&About', self)
        self.about.setShortcut('Ctrl+A')
        self.about.setStatusTip('About software')
        self.connect(self.about, QtCore.SIGNAL('triggered()'), self.showAbout)

        self.fileMenu3 = self.menubar.addMenu('&APPA')
        self.fileMenu3.addAction(self.about)
        #---------------------------------------------------------------------#



        #----------------Creation of statusBar--------------------#
        self.setStatusBar(QtGui.QStatusBar())
        #---------------------------------------------------------#


        #-------Creation of CentralWidget-----------------------------------------#
        self.widget = QtGui.QWidget()
        self.widget_layout = QtGui.QGridLayout()
        self.widget.setLayout(self.widget_layout)


        self.box1 = QtGui.QGroupBox()
        self.box1layout = QtGui.QGridLayout()
        self.box1.setLayout(self.box1layout)
        
        self.lbltitle = QtGui.QLabel("Abinit Post-Process Application")
        self.lbltitle.setFont(QtGui.QFont("calibri", 25))
        self.lbltitle.setFixedWidth(520);
        self.box1layout.addWidget(self.lbltitle,1,0)

        self.tab = QtGui.QTabWidget()
        self.tab.setTabsClosable (True)
        self.connect(self.tab,QtCore.SIGNAL('tabCloseRequested (int)'),self.closeTab)
        self.tab.setTabPosition(1)

        #----------Try to open the last .nc and .HIST files------------#
        MD_file = Read.MolecularDynamicFile("")
        if MD_file.isGoodFile():
            self.page1 = MD.Netcdf_MD(MD_file,self.units)
            self.tab.addTab(self.page1,MD_file.getNameFile())

        #----------Try to open the last Ground State file(BETA)--------#
        #GS_file = Read.outputFile("")
        #if GS_file.isGoodFile():
        #    self.page2 = GS.Ouput_GS(GS_file,self.units)
        #    self.tab.addTab(self.page2,str(GS_file.getNameFile()))

        #Connection of Signal (for the threading):
        self.connect(self,QtCore.SIGNAL("Reading(PyQt_PyObject)"),self.add)


        self.widget_layout.addWidget(self.box1,1,0,1,2)
        self.widget_layout.addWidget(self.tab,2,0,5,2)

        self.setCentralWidget(self.widget)
        #------------------------------------------------------------------------#
        self.show()
        
        if self.tab.count() == 0:
            self.showDialog()
        

#----------------------------------Methods---------------------------------------------#
    def showDialog(self):
        path = QtGui.QFileDialog.getOpenFileName(self, 'Open file', var.path(), "FILE (*_HIST *_OUT.nc *.out* *HIST.nc)")
        pathFile=str(path)
        var.global_path = pathFile
        del path
        
        if pathFile !="":
            if pathFile.find(' ') != -1 :
                #Sometimes the space caracter in the pathfile have to be replace by '\ '
                if os.path.exists(pathFile.replace(' ','\ ')):
                    pathFile = pathFile.replace(' ','\ ')

            #Check the existence of the both Netcdf Files before reading:
            if ( pathFile.find('_OUT.nc') != -1 and os.path.exists(pathFile.replace('_OUT.nc','_HIST')) )\
                    or ( pathFile.find('_HIST') != -1 and os.path.exists(pathFile.replace('_OUT.nc','_OUT.nc')) ):

                            self.read = threading.Thread(target=self.read, args=(pathFile,))
                            self.read.setDaemon(True)
                            self.read.start()      
                            del self.read

                            self.progressbar = load.loading(message="Reading output")
                            self.progressbar.show()                        
                            self.progressbar.raise_()
                            return
                        
            #Read the ASCII FILE:
            elif (pathFile.find('.out') != -1):
                    
                self.read = threading.Thread(target=self.read, args=(pathFile,))
                self.read.setDaemon(True)
                self.read.start()
                del self.read
                
                self.progressbar = load.loading(message="Reading output")
                self.progressbar.show()                        
                self.progressbar.raise_()
                return

                #------------BETA------------#
                # TEMP = Read.outputFile(pathFile)
                # if TEMP.isGoodFile():
                #     self.GS_page = GS.Ouput_GS(TEMP,self.units)
                #     self.tab.addTab(self.GS_page,str(TEMP.getNameFile()))
                #     self.tab.setCurrentIndex(self.tab.indexOf(self.GS_page))
                #     return

            else:
                self.showError("This file can't be read by APPA")

            #The file is not good for appa
            if (self.goodFile == False):
                self.showError("This file can't be read by APPA")
                try:
                    del self.progressbar
                except:
                    pass;


    def read(self,pathFile):
        TEMP = Read.MolecularDynamicFile(pathFile)
        if TEMP.isGoodFile():
            self.emit(QtCore.SIGNAL("Reading(PyQt_PyObject)"), TEMP)
        try:
           del self.progressbar
        except:
            pass;

    def add(self,pfile):
        MD_page = MD.Netcdf_MD(pfile,self.units)
        self.tab.addTab(MD_page,str(pfile.getNameFile()))
        self.tab.setCurrentIndex(self.tab.indexOf(MD_page))        
        return
        
    def showElastics(self):
        self.page3 = ElasticConstant.Elastic()
        self.tab.addTab(self.page3,"Elastic constant")
        self.tab.setCurrentIndex(self.tab.indexOf(self.page3))

    def showAbout(self):
        self.aboutPage = About.About()
        self.aboutPage.raise_()

    def showSave(self):
        fname = QtGui.QFileDialog.getSaveFileName(self,"Save Graphics",os.getcwd(), "FILE")
        if (fname !=""):
            print 'test'+fname
            try:
               Write.SaveFile(fname).saveData(self.tab.currentWidget().getData())
            except:
               self.showError("This file is not correct or no file open")


    def showExport(self):
        try:
            fname = QtGui.QFileDialog.getSaveFileName(self,"Export data",os.getcwd(), "XYZ file (*.xyz)")
            if (fname !=""):
                if 'xyz' in fname.split('.'):
                    pass
                else:
                    fname += '.xyz'
                pos   = (self.tab.currentWidget().getFile()).getXCart()   * 0.5291772085936 # Angstrom
                acell = (self.tab.currentWidget().getFile()).getAcell() * 0.5291772085936 # Angstrom
                typat = (self.tab.currentWidget().getFile()).getTypat()
                znucl = (self.tab.currentWidget().getFile()).getZnucl()
                Write.SaveFile(fname).xyzFormat(pos,acell,typat,znucl)
        except:
            self.showError("This file is not molecular dynamics file")


    def showError(self,perror):
        QtGui.QMessageBox.critical(self,"Warning",perror)



    def showConv(self):
        self.conv.showUnits()
        self.conv.raise_()




    def changeUnits(self,punits):
        self.units = punits
        for i in range(self.tab.count()):
            self.tab.widget(i).updateUnits(self.units)




    def closeTab(self,index):
        reply = QtGui.QMessageBox.question(self, 'Warning',
            "Are you sure you want to close?", QtGui.QMessageBox.Yes |
            QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            try:
                self.tab.currentWidget().restart()
                self.tab.currentWidget().closeGraphic()
            except:
                pass
            self.tab.removeTab(index)


    def closeEvent(self, event):
        sys.exit(0)
        event.accept()



    def center(self):
        screen = QtGui.QDesktopWidget().screenGeometry()
        size =  self.geometry()
        self.move((screen.width()-size.width())/2, (screen.height()-size.height())/2)

#----------------------------------------------------------------------------------------------------#
