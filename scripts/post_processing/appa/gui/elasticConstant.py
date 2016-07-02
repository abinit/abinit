#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os,re

#gui
import button as Button

#Reading
import reading.read as Read

#Utility
import utility.analysis as Analysis
import utility.global_variable as var


from PyQt4 import Qt,QtGui,QtCore
from numpy import linspace,zeros,dot

#---------------------------------------------------#
#---------------------------------------------------#
#----------------ELASTIC CLASS----------------------#
#---------------------------------------------------#
#---------------------------------------------------#
            
class Elastic(QtGui.QStackedWidget):

    def __init__(self, parent = None):
        QtGui.QStackedWidget.__init__(self, parent)
        self.files = {}
        self.te = {}
        self.ni = 200

        self.initUI()
      
	    
    def initUI(self):
        self.page0 = QtGui.QWidget()
        self.page0_layout = QtGui.QGridLayout()
        self.page0.setLayout(self.page0_layout)
        self.lbltitle = QtGui.QLabel("Abinit elastic constant  calculation :")
        self.lbltitle.setFont(QtGui.QFont("calibri", 25))
        self.lbltitle.setFixedWidth(570);
	    
        self.lbl1 = QtGui.QLabel("Number of files : ")
        self.lbl1.setFixedSize(100,20)
        self.sb1 = QtGui.QSpinBox()
        self.sb1.setFixedSize(50,20)
        self.sb1.setMaximum(6) 
        self.sb1.setMinimum(1)
        self.next1 = QtGui.QPushButton('&next', self)
        self.next1.setFixedSize(150,20)
        self.connect(self.next1 ,QtCore.SIGNAL("clicked()"),self.createPage1)
        
        self.page0_layout.addWidget(self.lbltitle,1,0,1,2,QtCore.Qt.AlignCenter)
        self.page0_layout.addWidget(self.lbl1,2,0,1,1,QtCore.Qt.AlignRight)
        self.page0_layout.addWidget(self.sb1,2,1,1,1,QtCore.Qt.AlignLeft)
        self.page0_layout.addWidget(self.next1,3,0,1,2,QtCore.Qt.AlignCenter)
			
			
        self.addWidget(self.page0)
        self.setCurrentIndex(self.indexOf(self.page0))

    def createPage1(self):
	    
        self.nbFiles = self.sb1.value()	
        try:
            self.removeWidget(self.page1)
        except:
            pass
		
        self.page1 = QtGui.QWidget()
        self.page1_layout = QtGui.QGridLayout()
        self.page1.setLayout(self.page1_layout)
        
        self.lbl2 = QtGui.QLabel("File without Deformation : ")
        self.lbl2.setFixedSize(160,20)
        self.teWD = QtGui.QTextEdit()
        self.teWD.setReadOnly(True)
        self.teWD.setFixedSize(330,36)
        self.te['WD'] = self.teWD
        self.browseWD = Button.Button('&Browse',0, self)
        self.browseWD.setFixedSize(100,20)
        self.connect(self.browseWD ,QtCore.SIGNAL("clicked()"),self.openWDFile)	
        self.page1_layout.addWidget(self.lbl2,0,0,1,1,QtCore.Qt.AlignCenter)
        self.page1_layout.addWidget(self.teWD,0,1,1,1,QtCore.Qt.AlignLeft)
        self.page1_layout.addWidget(self.browseWD,0,2,1,1,QtCore.Qt.AlignCenter)
        
        self.preview1 = QtGui.QPushButton('&Previous', self)
        self.connect(self.preview1 ,QtCore.SIGNAL("clicked()"),self.showPage_0)
        self.preview1.setFixedSize(130,20)
        self.next2 = QtGui.QPushButton('&next', self)
        self.next2.setFixedSize(130,20)
        self.connect(self.next2 ,QtCore.SIGNAL("clicked()"),self.verif_page1)
	
        for i in range(self.nbFiles-1):
            globals()['a%s' % i] = QtGui.QLabel("Deformation "+str(i+1)+ " : ")
            globals()['a%s' % i].setFixedSize(160,20)
            globals()['b%s' % i] = QtGui.QTextEdit()
            globals()['b%s' % i].setFixedSize(330,36)
            globals()['b%s' % i].setReadOnly(True)
            globals()['c%s' % i] = Button.Button('&Browse',i,parent=self)
            globals()['c%s' % i].setFixedSize(100,20)
            self.te[i] = globals()['b%s' % i]
            self.connect(globals()['c%s' % i] ,QtCore.SIGNAL("clicked()"),globals()['c%s' % i].clic)
            self.connect(globals()['c%s' % i] ,QtCore.SIGNAL("change(int)"),self.openFile)
            self.page1_layout.addWidget(globals()['a%s' % i],i+1,0,QtCore.Qt.AlignCenter)
            self.page1_layout.addWidget(globals()['b%s' % i],i+1,1,QtCore.Qt.AlignCenter)
            self.page1_layout.addWidget(globals()['c%s' % i],i+1,2,QtCore.Qt.AlignCenter)

        self.lblWD = QtGui.QLabel("Dataset without Deformation : ")
        self.lblWD.setFixedSize(190,20)
        self.sbWD = QtGui.QSpinBox()
        self.sbWD.setValue(-1)
        self.sbWD.setFixedSize(70,20)	    
        self.sbWD.setMinimum(1)
        self.lblWD.hide()
        self.sbWD.hide()
	
        self.page1_layout.addWidget(self.lblWD,self.nbFiles,0,1,1,QtCore.Qt.AlignRight)
        self.page1_layout.addWidget(self.sbWD,self.nbFiles,1,1,2,QtCore.Qt.AlignLeft)	
        self.page1_layout.addWidget(self.preview1,self.nbFiles+1,1,1,1,QtCore.Qt.AlignRight)
        self.page1_layout.addWidget(self.next2,self.nbFiles+1,2,1,1,QtCore.Qt.AlignLeft)
	
        self.addWidget(self.page1)	
        self.showPage_1()
	
    def createPage3(self):
        self.page2 = QtGui.QWidget()
        self.page2_layout = QtGui.QGridLayout()
        self.page2.setLayout(self.page2_layout)
        self.lblni = QtGui.QLabel('initial step : ')
        self.lblni.setFixedSize(70,20)
        self.sbni = QtGui.QSpinBox()	
        self.sbni.setMinimum(1)
        self.sbni.setMaximum(self.maxNi())
        self.sbni.setFixedSize(70,20)
        self.sbni.setValue(self.ni)	
        self.connect(self.sbni ,QtCore.SIGNAL("valueChanged(int)"),self.refresh_page3)
        self.preview2 = QtGui.QPushButton('&previous', self)
        self.preview2.setFixedSize(150,20)
        self.connect(self.preview2 ,QtCore.SIGNAL("clicked()"),self.showPage_1)


        self.elasticC = Analysis.elasticCalculation(self.files,self.ni,self.sbWD.value()-1).getElasticsConstants()

        if len(self.elasticC) == 0 :
            self.showPage_1()  
            return
	 
        i = 1
        j = 0
        for key in self.elasticC:
            globals()['a%s' % key] = QtGui.QLabel(str(key)+ "(GPa) : ")
            globals()['a%s' % key].setFixedSize(70,20)
            globals()['b%s' % key] = QtGui.QTextEdit()
            globals()['b%s' % key].setFixedSize(100,30)
            globals()['b%s' % key].setReadOnly(True)
            globals()['b%s' % key].setText(str(round(self.elasticC[key],2)))

            self.page2_layout.addWidget(globals()['a%s' % key],i,j,QtCore.Qt.AlignRight)
            self.page2_layout.addWidget(globals()['b%s' % key],i,j+1,QtCore.Qt.AlignLeft)

            if j == 4:
                j = 0
                i += 1
            else :
                j += 2

        if self.nbFiles > 1:
            self.page2_layout.addWidget(self.lblni,i+1,0,1,3,QtCore.Qt.AlignRight)
            self.page2_layout.addWidget(self.sbni,i+1,3,1,3,QtCore.Qt.AlignLeft)
            self.page2_layout.addWidget(self.preview2,i+2,0,1,6,QtCore.Qt.AlignCenter)
            
            self.addWidget(self.page2)
            self.setCurrentIndex(self.indexOf(self.page2))
	
	
    def showPage_1(self):
        self.setCurrentIndex(self.indexOf(self.page1))

    def showPage_0(self):
        self.setCurrentIndex(self.indexOf(self.page0))
    
    def verif_page1(self):
        if len(self.files) == self.nbFiles:
            if len(self.files) == 1 :
                val = self.sbWD.value()

                try:
                    if val > self.files['WD'].getNbDataset(): 
                        print 'error, dataset number'
                        return
                except :
                    print 'error, no dataset in this file'
                    return
                self.createPage3()
            else:
                self.createPage3()
                return
        else:
            print 'error, not enough files'
	    
    def type(self,fic):
        try:
            fic.getNbDataSet()
            return 'GS'
        except:
            return 'MD'
	    
	    
    def refresh_page3(self):	    
        try:
            self.ni = self.sbni.value()
            self.elasticC = Analysis.elasticCalculation(self.files,self.ni,self.sbWD.value()-1).getElasticsConstants()
            for key in self.elasticC:
                globals()['b%s' % key].setText(str(round(self.elasticC[key],2)))	    	
        except:
            pass
	
    def openWDFile(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file', var.path(), "FILE (*.nc *.HIST *.out)")
        #Save the current path
        pathFile=str(fname)
        var.global_path = pathFile

        if fname !="":
            fname = fname.replace(' ','\ ')
            if str(fname).find('_OUT.nc') != -1 :	
                TEMP = Read.MolecularDynamicFile(fname)
                if TEMP.isGoodFile():
                    self.files['WD']=TEMP
                    self.teWD.setText(str(fname))
                    self.lblWD.hide()
                    self.sbWD.hide()

            else:
                TEMP = Read.outputFile(fname)
                if TEMP.isGoodFile():
                    self.files['WD']=TEMP
                    self.teWD.setText(str(fname))		   
                    self.lblWD.show()
                    self.sbWD.show()


    def openFile(self,num):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file',var.path(), "FILE (*.nc .HIST *.out)")
        #Save the current path
        pathFile=str(fname)
        var.global_path = pathFile
        if fname !="":
            fname = fname.replace(' ','\ ')
            if str(fname).find('_OUT.nc') != -1 :
                TEMP = Read.MolecularDynamicFile(fname)
                if TEMP.isGoodFile():
                    self.files[num]=TEMP
                    self.te[num].setText(str(fname))		   
                    self.lblWD.hide()
                    self.sbWD.hide()
            else:
                TEMP = Read.outputFile(fname)
                if TEMP.isGoodFile():
                    self.files[num]=TEMP
                    self.te[num].setText(str(fname))	
                    self.lblWD.show()
                    self.sbWD.show()
		   		   	
    def maxNi(self):
        try:
            val = len(self.files['WD'].getStress())
            for i in range(len(self.files)-1):
                val2 = len(self.files[i].getStress())
                if (val > val2):
                    val = val2
		
                    return val
        except:
            return 1


    def goodFiles(self):
        for i in range(len(self.files)-1):
            if self.files[i].isGoodFile() == False:
                return False
        return True
	
    
   



    
