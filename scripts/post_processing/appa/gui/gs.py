#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import sys,os,re
from PyQt4 import Qt,QtGui,QtCore
from numpy import linspace

#GUI 
import gui.graph as Graph

#Reading
import reading.read as Read

#Utility
import utility.write as Write
import utility.analysis as Analysis

#---------------------------------------------------#
#---------------------------------------------------#
#----------------GS CLASS---------------------------#
#---------------------------------------------------#
#---------------------------------------------------#
            
class Ouput_GS(QtGui.QWidget):

    def __init__(self,pfile,punits, parent = None):


	self.file2  = pfile
        self.units  = punits #Dictionary with the units ans the conversion
	self.isShow = '' 
	self.initUI(parent)
	
    def initUI(self,parent):
	 
	QtGui.QWidget.__init__(self, parent)
	
    	self.layout = QtGui.QGridLayout()
	self.setLayout(self.layout)
	
	self.lblNbdataset = QtGui.QLabel("Number of dataset: "+str(self.file2.getNbDataset()), self)
	self.lblDate = QtGui.QLabel("Date : "+str(self.file2.getDate()), self)
	self.lblHour = QtGui.QLabel("Starting at : "+str(self.file2.getHour()), self)
	self.lblName = QtGui.QLabel("Name: "+str(self.file2.getNameFile()), self)
	self.lblData1 = QtGui.QLabel("Data 1: ", self)
	self.lblData2 = QtGui.QLabel("Data 2: ", self)
	
	self.lblNbdataset.setFixedSize(250,25)
	self.lblDate.setFixedSize(250,25)
	self.lblHour.setFixedSize(250,25)
	self.lblName.setFixedSize(250,25)
	self.lblData1.setFixedSize(100,25)
	self.lblData2.setFixedSize(100,25)
	
	self.CBox1 = QtGui.QComboBox()
	self.CBox1.addItem("")
	self.CBox1.addItem("DATASET")
	self.CBox1.addItem("ENERGY")
	self.CBox1.addItem("PRESSURE")
	self.CBox1.addItem("VOLUME")
	self.CBox1.addItem("ECUT")
	self.CBox1.addItem("A")
	self.CBox1.addItem("B")
	self.CBox1.addItem("C")
	self.connect(self.CBox1,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.changeData)
		
	self.CBox2 = QtGui.QComboBox()
	self.CBox2.addItem("")
	self.CBox2.addItem("DATASET")
	self.CBox2.addItem("ENERGY")
	self.CBox2.addItem("PRESSURE")
	self.CBox2.addItem("VOLUME")
	self.CBox2.addItem("ECUT")
	self.CBox2.addItem("A")
	self.CBox2.addItem("B")	
	self.CBox2.addItem("C")	
	self.connect(self.CBox2,QtCore.SIGNAL('currentIndexChanged(const QString&)'),self.changeData)
	
	self.checkbox1 = QtGui.QCheckBox("Show Graphic")
        self.connect(self.checkbox1,QtCore.SIGNAL('clicked()'),self.GraphOutput)
	self.checkbox1.setFixedSize(250,25)
	
	self.tbro = QtGui.QTextBrowser()
	
	self.birchCheckbox = QtGui.QCheckBox('&fit Energy / Volume curve', self)
        self.birchCheckbox.setStatusTip('Calculation of Bo by fitting Energy/Volume curve with Birch Murnaghan Equation ')
	self.birchCheckbox.setFixedSize(250,25)
	
	self.connect(self.birchCheckbox, QtCore.SIGNAL('clicked()'), self.showBirch)
	
	self.layout.addWidget(self.lblName,1,0,1,2)
	self.layout.addWidget(self.lblDate,2,0,1,2)
	self.layout.addWidget(self.lblHour,3,0,1,2)
	self.layout.addWidget(self.lblNbdataset,4,0,1,2)
	self.layout.addWidget(self.lblData1,5,0,QtCore.Qt.AlignLeft)
	self.layout.addWidget(self.lblData2,6,0,QtCore.Qt.AlignLeft)	
	self.layout.addWidget(self.CBox1,5,1,QtCore.Qt.AlignLeft)
	self.layout.addWidget(self.CBox2,6,1,QtCore.Qt.AlignLeft)
	self.layout.addWidget(self.checkbox1,7,0,1,2,QtCore.Qt.AlignLeft)
	self.layout.addWidget(self.birchCheckbox,8,0,1,2,QtCore.Qt.AlignLeft)
	self.layout.addWidget(self.tbro,1,2,10,4)
	self.birchCheckbox.setVisible(False)

    def isGoodFile(self):
	return self.file2.isGoodFile()
    	
    def setFile(self,pfile):
	self.file2 = pfile
		
    def restart(self):
	if(self.file2.isGoodFile()):
	    self.lblNbdataset.setText("Number of dataset: "+str(self.file2.getNbDataset()))
	    self.lblDate.setText("Date : "+str(self.file2.getDate()))
	    self.lblHour.setText("Starting at : "+str(self.file2.getHour()))
	    self.lblName.setText("Name: "+str(self.file2.getName()))
	    self.CBox1.setCurrentIndex(0)
	    self.CBox2.setCurrentIndex(0)
	    self.checkbox1.setChecked(False)
	    self.tbro.clear()
	    try:
 	        del self.Graph1
	    except:
	        pass 

    def toString(self,xlabel,x,ylabel,y):
	res = ""
	tab = ""
	if (len(x)!=0 and len(y)!=0):
   	    if(len(x) == len(y)):
		if (len(xlabel) == 0 and len(ylabel) == 0):
		    tab = ' '	    
	        elif (len(xlabel) >= 8 and len(str(x[0])) <= 8 ) :
		    tab = ' \t\t '
	            res+=xlabel+'\t'+ylabel+'\n'
		    
	        elif (len(xlabel) <= 8 and len(str(x[0])) >= 8 ) :
		    tab = ' \t '
	            res+=xlabel+'\t\t'+ylabel+'\n'
		    		    
	        elif (len(xlabel) > 8 and len(str(x[0])) > 8):
		    tab = '\t'
		    res+=xlabel+' \t '+ylabel+'\n'		    
					
		else:
		    tab = '\t'
		    res+=xlabel+'\t'+ylabel+'\n'
		    		    			
		for i in range(len(x)):
		    res+=str(x[i])+tab+str(y[i])+'\n'
		    
	        return res
		
	    else:
	        self.showError("Datas haven't the same size") 
		return ""
	else:
	    if (len(x)==0):
		xlabel = ylabel
		x = y
	    res+=xlabel+'\n'
	    for i in range(len(x)):
		res+=str(x[i])+'\n'
	    return res

    def GraphOutput(self):
	#The purpose of this fonction is to create or update
	#the graph of the output file
	if self.checkbox1.isChecked():
		
	    boxData = self.getBoxData()    
	    xlbl = boxData[0]
	    ylbl = boxData[1]
	    x = boxData[2]
	    y = boxData[3]
	    if (len(x)!=len(y) or (len(x)==0 or len(y) ==0)):
	        return 
	    
	    try:
		self.graph.update(x,y,xlbl,ylbl,name = str(self.file2.getNameFile()) + ' ' + ylbl)
	    except:
		self.graph=Graph.graphic(x,y,xlbl,ylbl,average=False,name = str(self.file2.getNameFile()) + ' ' + ylbl)
		self.connect(self.graph, QtCore.SIGNAL("myCustomizedSignal()"), self.closeGraph)
		self.graph.show()
	    if self.birchCheckbox.isChecked():
	        self.showBirch()  
	else:
	    self.closeGraph()

    def changeData(self):
	    
	boxData = self.getBoxData()
	xlabel = boxData[0]
	ylabel = boxData[1]
	data1 = boxData[2]
	data2 = boxData[3]
	self.tbro.setText(self.toString(xlabel,data1,ylabel,data2))
	self.isShow = 'GS'

		
	if (str(xlabel).find("ENERGY") != -1 and str(ylabel).find("VOLUME")  != -1): 
	    self.birchCheckbox.setVisible(True)  
	elif ( str(ylabel).find("ENERGY") != -1 and str(xlabel).find("VOLUME")  != -1):
	    self.birchCheckbox.setVisible(True)          
	else :
	    self.birchCheckbox.setVisible(False)
	    self.birchCheckbox.setCheckState(QtCore.Qt.Unchecked)

 	if (len(data1)==len(data2) and (len(data1)!=0 and len(data2) !=0)):
	    self.GraphOutput()   
    
    def getBoxData(self):
	xlabel = self.CBox1.currentText()
	ylabel = self.CBox2.currentText()
	data1 = []
	data2 = []
	
	if xlabel == '':
	    pass
	elif xlabel == "DATASET":
	    Nbdata = self.file2.getNbDataset()
	    data1 = linspace(1,Nbdata,Nbdata)
	elif xlabel == "ENERGY":
	    xlabel+=" ("+str(self.units['Energy'][0])+")"
	    data1 = self.file2.getE_Tot() * self.units['Energy'][1]
	elif xlabel == "PRESSURE":
	    xlabel+=" ("+str(self.units['Pressure'][0])+")"
	    data1 = self.file2.getPress() * self.units['Pressure'][1]
	elif xlabel == "VOLUME":
	    xlabel+=" ("+str(self.units['Volume'][0])+")"
	    data1 = self.file2.getVol() * self.units['Volume'][1]
	elif xlabel == "ECUT":
	    xlabel+=" ("+str(self.units['Energy'][0])+")"
	    data1 = self.file2.getEcut() * self.units['Energy'][1]
	elif xlabel == "A":
	    xlabel+=" ("+str(self.units['Lengh'][0])+")"
	    data1 = self.file2.getA()
	elif xlabel == "B":
	    xlabel+=" ("+str(self.units['Lengh'][0])+")"
	    data1 = self.file2.getB()
	elif xlabel == "C":
	    xlabel+=" ("+str(self.units['Lengh'][0])+")"
	    data1 = self.file2.getC()		
	
	if ylabel == '':
	    pass				
	if ylabel == "DATASET":
	    Nbdata = self.file2.getNbDataset()
	    data2 = linspace(1,Nbdata,Nbdata)
	elif ylabel == "ENERGY":
	    ylabel+=" ("+str(self.units['Energy'][0])+")"
	    data2 = self.file2.getE_Tot() * self.units['Energy'][1]
	elif ylabel == "PRESSURE":
	    ylabel+=" ("+str(self.units['Pressure'][0])+")"
	    data2 = self.file2.getPress() * self.units['Pressure'][1]
	elif ylabel == "VOLUME":
	    ylabel+=" ("+str(self.units['Volume'][0])+")"
	    data2 = self.file2.getVol() *  self.units['Volume'][1] 
	elif ylabel == "ECUT":
	    ylabel+=" ("+str(self.units['Energy'][0])+")"
	    data2 = self.file2.getEcut() * self.units['Energy'][1]
	elif ylabel == "A":
	    ylabel+=" ("+str(self.units['Lengh'][0])+")"
	    data2 = self.file2.getA() 
	elif ylabel == "B":
	    ylabel+=" ("+str(self.units['Lengh'][0])+")"
	    data2 = self.file2.getB() 
	elif ylabel == "C":
	    ylabel+=" ("+str(self.units['Lengh'][0])+")"
	    data2 = self.file2.getC()     
		
	return [xlabel,ylabel,data1,data2]
	    
    def showBirch(self):
 
	vol = self.file2.getVol()
	
	Energie = self.file2.getE_Tot()
        
	if (len(vol)!=len(Energie) or (len(vol)==0 or len(Energie) ==0)):
	        return 
	
	data = self.toString('',vol,'',Energie)
	
	try :
	    self.value
	except:
	    self.value = Analysis.birch(data,len(vol)).getValue()
	    pass
	
	V0 = self.value['V0']  
	B0 = self.value['B0'] / 29421.033e0
	E0 = self.value['E0']
	dB0 = self.value['dB0']
	RMS = self.value['RMS']
	
	x = linspace( min(vol), max(vol), 1000) 

	y = E0 + ( (B0 * x) / dB0 ) * (  ((V0/x)**(dB0)) / (dB0-1) + 1  ) - ( (B0*V0)/(dB0-1) )
	
	xlabel ="VOLUME ("+str(self.units['Volume'][0])+")"
	
	ylabel ="ENERGY ("+str(self.units['Energy'][0])+")"
	
	x *=  self.units['Volume'][1]
	
	y *=  self.units['Energy'][1]

	#Update the textBox:
	self.res  = 'Results of the fit : \n\n'	
	
	self.res += 'B0   = '  +str(self.value['B0'] * self.units['Pressure'][1]) +' ' +str(self.units['Pressure'][0]) + '\n'
	self.res += 'B0\'  = ' +str(dB0)+'\n'
	self.res += 'E0   = '  +str(E0 * self.units['Energy'][1])+' ' +str(self.units['Energy'][0]) + '\n'
	self.res += 'V0   = '  +str(V0 * self.units['Volume'][1])+'  '+str(self.units['Volume'][0])+'\n'
	self.res += 'RMS = '   +str(RMS* self.units['Energy'][1])+' ' +str(self.units['Energy'][0]) + '\n'
	     
	self.tbro.setText(self.res)
    	self.isShow = 'BIRCH'

	#update the graph:
	if self.birchCheckbox.isChecked():
	    boxData = self.getBoxData()
	    xlabel  = boxData[0]
	    ylabel  = boxData[1]
	    data1   = boxData[2]
	    data2   = boxData[3]
	    try :
		self.graph.setPlot(data1,data2,xlabel,ylabel,point = True)
						
		if xlabel == "ENERGY ("+str(self.units['Energy'][0])+")" :	            
		    self.graph.addPlot(y,x)
		else:
		    self.graph.addPlot(x,y)
		    
		self.graph.addLegend(["Simulation Data","Fit"])
	    except:
		pass
	else:
	    self.changeData()
	        
	
    def closeGraph(self):
	try:
	    self.checkbox1.setCheckState(QtCore.Qt.Unchecked)
	    del self.graph
	except:
	    pass
		    
    def closeBirchGraph(self):
	try:
	    self.birch = False
	    del self.graphBirch
	except:
	    pass
		    
    def updateUnits(self, punits):
        self.units =  punits
	if self.isShow == 'GS' :
	    self.showBirch()
	    self.changeData()

	elif self.isShow == 'BIRCH' :
	    self.changeData()
            self.showBirch()
    	    
	    
    def showError(self,perror):
        QtGui.QMessageBox.critical(self,"Warning",perror)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
