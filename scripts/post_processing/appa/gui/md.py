#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: january 2014
"""

import sys,os,re
import threading

#GUI 
import gui.graph  as Graph
import gui.vdos   as VDOS
import gui.rdf    as RDF
import gui.msd    as MSD
import gui.netcdf as NETCDF
import gui.atom_pos   as POS
#Reading
import reading.read as Read

#Utility
import utility.analysis as Analysis
import utility.write    as Write
import utility.writeHIST as WriteHIST
import utility.global_variable as Var

#Fortran code
import fortran.math as Math

from PyQt4 import Qt,QtGui,QtCore
from numpy import linspace,hstack

#---------------------------------------------------#
#---------------------------------------------------#
#----------------MD CLASS---------------------------#
#---------------------------------------------------#
#---------------------------------------------------#
            
class Netcdf_MD(QtGui.QWidget):


#----------------------Constructor---------------------------#
    #This class show the panel of the Molacular dynamics
    #Parameters :
    #                   - pfile  => objet of NetcdfFile
    #                   - punits => array with units
    #You can see the result of your simulation, you can choose the initial and final step.
    #You also can display severals graphics.

    def __init__(self, pfile, punits, parent = None):

        self.file1 = pfile
        self.ni = self.file1.getNi()#Departure of the dataset
        self.nf = self.file1.getNf()#End of the dataset
        self.units = punits #Dictionary with the units ans the conversion
        
        QtGui.QWidget.__init__(self, parent)
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
                
        self.box2 = QtGui.QGroupBox("Simulation data:")
        self.box2layout = QtGui.QGridLayout()
        self.box2.setLayout(self.box2layout)        

        E_tot = self.units['Energy'][1] * Math.average(self.file1.getE_Tot())
        deviation = self.file1.getStandardDeviation(self.units['Energy'][1] *self.file1.getE_Tot(), E_tot)      
        strETOT = self.displayData(E_tot,deviation)

        vol = self.units['Volume'][1] * Math.average(self.file1.getVol())
        deviation = self.file1.getStandardDeviation(self.units['Volume'][1] *self.file1.getVol(), vol)        
        strVol = self.displayData(vol,deviation)

        Temp = self.file1.getTemp()        
        ATemp = Math.average(Temp) - self.units['Temperature'][1]        
        deviation = Math.standard_deviation( Temp - self.units['Temperature'][1] , ATemp)
        strTemp = self.displayData(ATemp,deviation)

        Press = Math.average( self.file1.getPress() ) * self.units['Pressure'][1]
        deviation = Math.standard_deviation(self.units['Pressure'][1] * self.file1.getPress(), Press)
        strPress = self.displayData(Press,deviation)

        self.lbl1 = QtGui.QLabel("Number of cell(s)  : "+str(self.file1.getNImage()), self)
        self.lbl2 = QtGui.QLabel("Number of atoms  : "+str(self.file1.getNatom()), self)
        self.lbl3 = QtGui.QLabel("Number of steps  : "+str(self.file1.getNbTime()), self)
        self.lbl4 = QtGui.QLabel("Total Energy ("+str(self.units['Energy'][0])+"): "+strETOT, self)                
        self.lbl5 = QtGui.QLabel("Volume ("+str(self.units['Volume'][0])+"): "+strVol, self)
        self.lbl6 = QtGui.QLabel("Temperature ("+str(self.units['Temperature'][0])+")  : "+strTemp, self)
        self.lbl7 = QtGui.QLabel("Pressure ("+str(self.units['Pressure'][0])+")    : "+strPress, self)        

        self.lbl1.setFixedWidth(300)
        self.lbl2.setFixedWidth(300)
        self.lbl3.setFixedWidth(300)
        self.lbl4.setFixedWidth(300)
        self.lbl5.setFixedWidth(300)
        self.lbl6.setFixedWidth(300)
        self.lbl7.setFixedWidth(300)
                
        self.box2layout.addWidget(self.lbl1,1,0)
        self.box2layout.addWidget(self.lbl2,2,0)
        self.box2layout.addWidget(self.lbl3,3,0)
        self.box2layout.addWidget(self.lbl4,4,0)
        self.box2layout.addWidget(self.lbl5,5,0)
        self.box2layout.addWidget(self.lbl6,6,0)
        self.box2layout.addWidget(self.lbl7,7,0)

        self.box3 = QtGui.QGroupBox("Option")
        self.box3layout = QtGui.QGridLayout()
        self.box3.setLayout(self.box3layout)
        
        self.lbl3_1 = QtGui.QLabel('Choose the initial and final step\nfor the molecular dynamics :')
        self.lbl3_1.setFixedSize(300,50)
        self.lblni = QtGui.QLabel('initial step')
        self.lblnf = QtGui.QLabel('final step')
        self.sbni = QtGui.QSpinBox()
        self.sbnf = QtGui.QSpinBox()
        self.sbni.setMaximum(self.file1.getNbTime())
        self.sbnf.setMaximum(self.file1.getNbTime()) 
        self.sbni.setMinimum(1)
        self.sbnf.setMinimum(1) 
        self.sbni.setFixedSize(150,20)
        self.sbnf.setFixedSize(150,20)
        self.sbni.setValue(self.ni)
        self.sbnf.setValue(self.nf)
        self.connect(self.sbni,QtCore.SIGNAL('keyPressed'),self.verifStep)
        self.connect(self.sbnf,QtCore.SIGNAL('keyPressed'),self.verifStep)

        self.pbok = QtGui.QPushButton("update")
        self.pbok.setFixedSize(100,20)
        self.connect(self.pbok,QtCore.SIGNAL("clicked()"),self.updateStep)
        
        self.box3layout.addWidget(self.lbl3_1, 1, 0, 1, 2,QtCore.Qt.AlignCenter)
        self.box3layout.addWidget(self.lblni, 2, 0, 1, 1,QtCore.Qt.AlignCenter)
        self.box3layout.addWidget(self.sbni, 2, 1, 1, 1,QtCore.Qt.AlignCenter)
        self.box3layout.addWidget(self.lblnf, 3, 0, 1, 1,QtCore.Qt.AlignCenter)
        self.box3layout.addWidget(self.sbnf, 3, 1, 1 , 1,QtCore.Qt.AlignCenter)
        self.box3layout.addWidget(self.pbok, 4, 0,1,2,QtCore.Qt.AlignCenter)
        
        
        self.box4 = QtGui.QGroupBox("Graphics")
        self.box4layout = QtGui.QGridLayout()
        self.box4.setLayout(self.box4layout)

        self.potentialEnergy = QtGui.QPushButton('Potential Energy', self)
        self.potentialEnergy.setStatusTip('Show graphic of  Potential energy')
        
        self.connect(self.potentialEnergy, QtCore.SIGNAL('clicked()'), self.showPotentialEnergy)

        self.totalEnergy = QtGui.QPushButton('Total Energy', self)
        self.totalEnergy.setStatusTip('Show graphic of total energy')
        self.connect(self.totalEnergy, QtCore.SIGNAL('clicked()'), self.showTotalEnergy)

        self.kineticEnergy = QtGui.QPushButton('Kinetic Energy', self)
        self.kineticEnergy.setStatusTip('Show graphic of kinetic energy')
        self.connect(self.kineticEnergy, QtCore.SIGNAL('clicked()'), self.showKineticEnergy)
        
        self.temperature = QtGui.QPushButton('Temperature', self)
        self.temperature.setStatusTip('Show graphic of temperature')
        self.connect(self.temperature, QtCore.SIGNAL('clicked()'), self.showTemperature)
        
        self.pressure = QtGui.QPushButton('Pressure', self)
        self.pressure.setStatusTip('Show graphic of pressure')
        self.connect(self.pressure, QtCore.SIGNAL('clicked()'), self.showPressure)

        self.stress = QtGui.QPushButton('Stress', self)
        self.stress.setStatusTip('Show graphic of stress')
        self.connect(self.stress, QtCore.SIGNAL('clicked()'), self.showStress)

        self.volume = QtGui.QPushButton('Volume', self)
        self.volume.setStatusTip('Show graphic of volume')
        self.connect(self.volume, QtCore.SIGNAL('clicked()'), self.showVolume)


        self.acell = QtGui.QPushButton('Cell', self)
        self.acell.setStatusTip('Show graphic of cell')
        self.connect(self.acell, QtCore.SIGNAL('clicked()'), self.showAcell)


        self.angles = QtGui.QPushButton('Angles', self)
        self.angles.setStatusTip('Show graphic of Angles')
        self.connect(self.angles, QtCore.SIGNAL('clicked()'), self.showAngles)

        
        self.VAF = QtGui.QPushButton('VACF', self)
        self.VAF.setStatusTip('Show graphic of Velocity autocorrelation function')
#        self.connect(self.VAF, QtCore.SIGNAL('clicked()'), self.showVAF)        
        
        self.DOS = QtGui.QPushButton('VDOS', self)
        self.DOS.setStatusTip('Show graphic of Vibrational Density Of States')
#        self.connect(self.DOS, QtCore.SIGNAL('clicked()'), self.showDOS)
                        
        self.RDF = QtGui.QPushButton('Topology', self)
        self.RDF.setStatusTip('Show graphic of Radial Distribution Function')
#        self.connect(self.RDF, QtCore.SIGNAL('clicked()'), self.showRDF)        
                
        self.position = QtGui.QPushButton('Positions', self)
        self.position.setStatusTip('Show positions of all the particules for all the step in 2D ')
#        self.connect(self.position, QtCore.SIGNAL('clicked()'), self.showPosition)

        self.MSD = QtGui.QPushButton('MSD(beta)', self)
        self.MSD.setStatusTip('Show graphic of Means Squared Displacement(beta)')
#        self.connect(self.MSD, QtCore.SIGNAL('clicked()'), self.showMSD)        

        self.netcdf = QtGui.QPushButton('save netcdf', self)
        self.netcdf.setStatusTip('Save molecular dynamics in necdf format')
#        self.connect(self.netcdf, QtCore.SIGNAL('clicked()'), self.showNetcdf)

        self.xyz = QtGui.QPushButton('save .xyz', self)
        self.xyz.setStatusTip('Save position in .xyz format')
        self.connect(self.xyz, QtCore.SIGNAL('clicked()'), self.showExport)
        
        self.box4layout.addWidget(self.potentialEnergy,1,0)
        self.box4layout.addWidget(self.totalEnergy,1,1)
        self.box4layout.addWidget(self.kineticEnergy ,1,2)
        self.box4layout.addWidget(self.temperature,1,3)
        self.box4layout.addWidget(self.pressure,1,4)
        self.box4layout.addWidget(self.stress,1,5)
        self.box4layout.addWidget(self.volume,2,0)
        self.box4layout.addWidget(self.acell,2,1)
        self.box4layout.addWidget(self.angles,2,2)
        self.box4layout.addWidget(self.VAF,3,0)
        self.box4layout.addWidget(self.DOS,3,1)
        self.box4layout.addWidget(self.RDF,3,2)
        self.box4layout.addWidget(self.MSD,3,3)
        self.box4layout.addWidget(self.position,2,3)
        self.box4layout.addWidget(self.netcdf,2,5)
        self.box4layout.addWidget(self.xyz,3,5)
        self.layout.addWidget(self.box2,1,0,5,1)
        self.layout.addWidget(self.box3,1,1,5,1)
        self.layout.addWidget(self.box4,6,0,1,2)
    
#-------------------------------------------------------------------------#        

#----------------------------------Method---------------------------------#        
    def setFile(self,pfile):
        self.file1 = pfile

    def getFile(self):
        return self.file1

    def closeEvent(self, event):
        reply = QtGui.QMessageBox.question(self, 'Message',
            "Are you sure you want to quit?", QtGui.QMessageBox.Yes | 
            QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            sys.exit(0)
            event.accept()
        else:
            event.ignore()

    def showError(self,perror):
        QtGui.QMessageBox.critical(self,"Warning",perror)

    def showMessage(self,pmessage):
        QtGui.QMessageBox.information(self,"Information",pmessage)
        
        
    def showPotentialEnergy(self):
        if self.file1.isGoodFile():
            x = linspace(self.ni,self.nf-1,self.nf-self.ni)# Temporarily !!!
            epot = self.file1.getE_pot() * self.units['Energy'][1]
            self.file1.getNameFile()
            self.GraphPotentialEnergy = Graph.graphic(x,epot,'Step',"Potential Energy ("+str(self.units['Energy'][0])+")",\
                                                          name = self.file1.getNameFile()+" Potential Energy")
            self.GraphPotentialEnergy.show()
        else:
            self.showError("Your NetCDF file is not correct")

    def showTotalEnergy(self):
        if self.file1.isGoodFile():
            x = linspace(self.ni,self.nf-1,self.nf-self.ni)# Temporarily !!!
            eTot = self.file1.getE_Tot() * self.units['Energy'][1]
            self.GraphTotalEnergy = Graph.graphic(x,eTot,'Step',"Total Energy("+str(self.units['Energy'][0])+")",\
                                                      name = self.file1.getNameFile()+" Total Energy")
            self.GraphTotalEnergy.show()
        else:
            self.showError("Your NetCDF file is not correct")

    def showKineticEnergy(self):
        if self.file1.isGoodFile():
            x = linspace(self.ni,self.nf-1,self.nf-self.ni)# Temporarily !!!
            eKin = self.file1.getE_kin()  * self.units['Energy'][1]
            self.GraphKineticEnergy = Graph.graphic(x,eKin,'Step',"Kinetic Energy ("+str(self.units['Energy'][0])+")",\
                                                        name = self.file1.getNameFile() +" kinetic Energy")
            self.GraphKineticEnergy.show()
        else:
            self.showError("Your NetCDF file is not correct") 

    def showTemperature(self):
        if self.file1.isGoodFile():
            x = linspace(self.ni,self.nf-1,self.nf-self.ni)# Temporarily !!!
            temp = self.file1.getTemp() - self.units['Temperature'][1]
            self.GraphTemperature = Graph.graphic(x,temp,'Step',"Temperature ("+str(self.units['Temperature'][0])+")",\
                                                      name = self.file1.getNameFile() +" Temperature")
            self.GraphTemperature.show()
        else:
            self.showError("Your NetCDF file is not correct")

    def showPressure(self):
        if self.file1.isGoodFile():
            x = linspace(self.ni,self.nf-1,self.nf-self.ni)# Temporarily !!!
            pressure = self.file1.getPress() * self.units['Pressure'][1]
            self.GraphPressure = Graph.graphic(x,pressure,'Step',"Pressure ("+str(self.units['Pressure'][0])+")",\
                                                   name = self.file1.getNameFile() +" Pressure")
            self.GraphPressure.show()
        else:
            self.showError("Your NetCDF file is not correct")
            
    def showStress(self):        
        if self.file1.isGoodFile():
            x = linspace(self.ni,self.nf-1,self.nf-self.ni)# Temporarily !!!
            st =self.file1.getStress()

            self.GraphStress= Graph.graphic(x,st * self.units['Pressure'][1],'Step', "Stress ("+str(self.units['Pressure'][0])+")",\
                                                average=False,name = self.file1.getNameFile() +" Stress")
            self.GraphStress.addLegend([r'$\sigma_1$',r'$\sigma_2$',r'$\sigma_3$',r'$\sigma_4$',r'$\sigma_5$',r'$\sigma_6$'])
            self.GraphStress.show()
                        
    def showVolume(self):
        if self.file1.isGoodFile():
            x = linspace(self.ni,self.nf-1,self.nf-self.ni)
            vol = self.file1.getVol() * self.units['Volume'][1]
            self.GraphVolume = Graph.graphic(x,vol,'Step',"Volume ("+str(self.units['Volume'][0])+")",\
                                                   name = self.file1.getNameFile() +" Volume")
            self.GraphVolume.show()
        else:
            self.showError("Your NetCDF file is not correct")


    def showAcell(self):
        if self.file1.isGoodFile():
            x = linspace(self.ni,self.nf-1,self.nf-self.ni)
            acell = self.file1.getAcell() * self.units['Distance'][1]
            self.GraphAcell = Graph.graphic(x,acell,'Step',"Acell ("+str(self.units['Distance'][0])+")",\
                                                   average=False,name = self.file1.getNameFile() +" Distance")
            self.GraphAcell.addLegend([r'$a$',r'$b$',r'$c$'])
            self.GraphAcell.show()
        else:
            self.showError("Your NetCDF file is not correct")

                    
    def showVAF(self):
        if self.file1.isGoodFile():            
            try:
                x = linspace(self.ni,self.nf,len(self.vacf))
                self.GraphVAF = Graph.graphic(x,self.vacf,'step',"VACF",average=False,name = self.file1.getNameFile() +" VACF")
                self.GraphVAF.show()
            except:
                self.vacf = Analysis.Correlation(self.file1.getVel()).getCorrelationFunction(normalize = True)
                x = linspace(self.ni,self.nf,len(self.vacf))
                self.GraphVAF = Graph.graphic(x,self.vacf,'step',"VACF",average=False,name = self.file1.getNameFile() +" VACF")
                self.GraphVAF.show()
        else:
            self.showError("Your NetCDF file is not correct")        

    def showAngles(self):
        if self.file1.isGoodFile():
            x = linspace(self.ni,self.nf-1,self.nf-self.ni)
            angles = self.file1.getAngles() * self.units['Angle'][1]
            self.GraphAngles = Graph.graphic(x,angles,'Step',"Angles ("+str(self.units['Angle'][0])+")",\
                                                   average=False,name = self.file1.getNameFile() +" Angles")
            self.GraphAngles.addLegend([r'$\alpha$',r'$\beta$',r'$\gamma$'])
            self.GraphAngles.show()
        else:
            self.showError("Your NetCDF file is not correct")

    def showDOS(self):
        if self.file1.isGoodFile():
            try:
                self.wdos = VDOS.winVDOS(self.vacf,name = self.file1.getNameFile() +" VDOS")
            except:
                self.vacf = Analysis.Correlation(self.file1.getVel()).getCorrelationFunction(normalize = True)
                self.wdos = VDOS.winVDOS(self.vacf,self.file1.getDtion(),name = self.file1.getNameFile() +" VDOS")
        else:
            self.showError("Your NetCDF file is not correct")
            
    def showRDF(self):
        if self.file1.isGoodFile():
            try:
                self.rdf.show()
                self.rdf.raise_()
            except:                
                self.rdf = RDF.winRDF(self.file1,name = self.file1.getNameFile() +" Radial pair distribution")
        else:
            self.showError("Your NetCDF file is not correct")        

    def showMSD(self):
        if self.file1.isGoodFile():
            try:
                self.msd.show()
                self.msd.raise_()
            except:                
                self.msd = MSD.winMSD(self.file1,name = self.file1.getNameFile() +" Mean squared displacement")
        else:
            self.showError("Your NetCDF file is not correct")        

            
    def showPosition(self):
        if self.file1.isGoodFile():
            try:
                self.pos.show()
                self.pos.raise_()
            except:                
                self.pos = POS.WinPOS(self.file1,self.units,name = self.file1.getNameFile() +"Atomic positions")
        else:
            self.showError("Your NetCDF file is not correct")        

    def showExport(self):
        try:
            fname = QtGui.QFileDialog.getSaveFileName(self,"Export data in xyz format",Var.path(), "XYZ file (*.xyz)")
            if (fname !=""):
                if 'xyz' in fname.split('.'):
                    pass
                else:
                    fname += '.xyz'
                pos   = self.file1.getXCart()   * 0.5291772085936 # Angstrom
                acell = self.file1.getAcell() * 0.5291772085936 # Angstrom
                typat = self.file1.getTypat()
                znucl = self.file1.getZnucl()
                Write.SaveFile(fname).xyzFormat(pos,acell,typat,znucl)
                self.showMessage("file save in "+ str(fname))
        except:
            self.showError("Unable to save in xyz format ")

    def showNetcdf(self):
        if self.file1.isGoodFile():
            try:
                fname = QtGui.QFileDialog.getSaveFileName(self,"Export data in xyz format",Var.path(), "HIST file (*_HIST)")
                if (fname !=""):
                    if 'HIST' in fname.split('_'):
                        pass
                    else:
                        fname += '_HIST'
                    WriteHIST.writeHIST(self.file1,fname,self.ni,self.nf)
                    self.showMessage("file save in "+ str(fname))
            except:
                self.showError("Unable to save your simulation ")
             
    def verifStep(self):
        if self.file1.isGoodFile():
            ni = self.sbni.value()
            nf = self.sbnf.value()
            if(ni>nf):
                self.showError("initial step must be inferior to final step")
                self.sbni.setValue(self.ni)
                self.sbnf.setValue(self.nf)
                return False
            if(nf==ni):
                self.showError("final step must be superior to initial step")
                self.sbni.setValue(self.ni)
                self.sbnf.setValue(self.nf)
                return False
            if((ni == self.ni)and(nf == self.nf)):
                return False
            return True
        else:
            self.showError("Your NetCDF file is not correct")
            return False
        

    def closeGraphic(self):
        try:
            del self.GraphPressure
        except:
            pass 
        try:
            del self.GraphPotentialEnergy
        except:
            pass 
        try:
            del self.GraphTotalEnergy
        except:
            pass 
        try:
            del self.GraphKineticEnergy
        except:
            pass 
        try:
            del self.GraphTemperature
        except:
            pass         
        try:
            del self.GraphStress
        except:
            pass         
        try:
            del self.GraphVAL
        except:
            pass 
        try:
            del self.GraphVol
        except:
            pass
        try:
            del self.GraphAcell
        except:
            pass
        try:
            del self.GraphAngles
        except:
            pass

    def updateUnits(self,punits):
           
        self.units = punits
        
        #-----Change the label and the value:------#        
        self.updateLabel()
        
        #---------Change the Graph units:----------#
        self.updateGraph( units = True )


    def updateStep(self):
  
        if self.verifStep():
                
            self.file1.setNi(self.sbni.value())
            self.file1.setNf(self.sbnf.value())
            self.ni = self.sbni.value()
            self.nf = self.sbnf.value()

            x = linspace(self.ni,self.nf-1,self.nf-self.ni)# Temporarily !!!
            
            #------Change the label :--------#        
            self.updateLabel()
            
            #------Change the Graph :--------#
            self.updateGraph()                    
        


    def updateLabel(self):
        E_tot = self.units['Energy'][1] * Math.average(self.file1.getE_Tot())
        deviation = self.file1.getStandardDeviation(self.units['Energy'][1] *self.file1.getE_Tot(), E_tot)      
        strETOT = self.displayData(E_tot,deviation)

        vol = self.units['Volume'][1] * Math.average(self.file1.getVol())
        deviation = self.file1.getStandardDeviation(self.units['Volume'][1] *self.file1.getVol(), vol)        
        strVol = self.displayData(vol,deviation)

        Temp = self.file1.getTemp()        
        ATemp = Math.average(Temp) - self.units['Temperature'][1]        
        deviation = Math.standard_deviation( Temp - self.units['Temperature'][1] , ATemp)
        strTemp = self.displayData(ATemp,deviation)

        Press = Math.average( self.file1.getPress() ) * self.units['Pressure'][1]
        deviation = Math.standard_deviation(self.units['Pressure'][1] * self.file1.getPress(), Press)
        strPress = self.displayData(Press,deviation)

        self.lbl4.setText("Total Energy ("+str(self.units['Energy'][0])+"): "+strETOT)
        self.lbl5.setText("Volume ("+str(self.units['Volume'][0])+"): "+strVol)
        self.lbl6.setText("Temperature ("+str(self.units['Temperature'][0])+")  : "+strTemp)
        self.lbl7.setText("Pressure ("+str(self.units['Pressure'][0])+")    : "+strPress)

        
    def updateGraph(self,units = False):

        x = linspace(self.ni,self.nf-1,self.nf-self.ni) # Temporarily !!!
        
        try:
            self.GraphPressure
            pressure = self.file1.getPress() * self.units['Pressure'][1]
            self.GraphPressure.update(x,pressure,'Step',"Pressure ("+str(self.units['Pressure'][0])+")",\
                                          name = self.file1.getNameFile()+" Pressure")
        except:
            pass 

        try:
            self.GraphPotentialEnergy
            epot = self.file1.getE_pot() * self.units['Energy'][1]
            self.GraphPotentialEnergy.update(x,epot,'Step',"Potential Energy ("+str(self.units['Energy'][0])+")",\
                                                 name = self.file1.getNameFile()+" Potential Energy")
        except:
            pass 
        
        try:
            self.GraphTotalEnergy
            eTot = self.file1.getE_Tot() * self.units['Energy'][1]
            self.GraphTotalEnergy.update(x,eTot,'Step',"Total Energy("+str(self.units['Energy'][0])+")",\
                                             name = self.file1.getNameFile()+" Total Energy")
        except:
            pass 

        try:
            self.GraphKineticEnergy
            eKin = self.file1.getE_kin()  * self.units['Energy'][1]
            self.GraphKineticEnergy.update(x,eKin,'Step',"Kinetic Energy ("+str(self.units['Energy'][0])+")",\
                                               name = self.file1.getNameFile()+" Kinetic Energy")
        except:
            pass 

        try:
            self.GraphTemperature
            temp = self.file1.getTemp() - self.units['Temperature'][1]
            self.GraphTemperature.update(x,temp,'Step',"Temperature ("+str(self.units['Temperature'][0])+")",\
                                             name = self.file1.getNameFile()+" Temperature")
        except:
            pass

        try:
            self.GraphVolume
            vol = self.file1.getVol()  * self.units['Volume'][1]
            self.GraphVolume.update(x,vol,'Step',"Volume ("+str(self.units['Volume'][0])+")",\
                                        name = self.file1.getNameFile()+" Volume")
        except:
            pass
        try:
            self.GraphStress
            st   = self.file1.getStress()  * self.units['Pressure'][1]
            self.GraphStress.update(x,st,'Step',"Stress ("+str(self.units['Pressure'][0])+")",name = self.file1.getNameFile()+" Stress")
            self.GraphStress.addLegend([r'$\sigma_1$',r'$\sigma_2$',r'$\sigma_3$',r'$\sigma_4$',r'$\sigma_5$',r'$\sigma_6$'])
        except:
            pass

        try:        
            self.GraphAcell
            acell = self.file1.getAcell()  * self.units['Distance'][1]
            self.GraphAcell.update(x,acell,'Step',"Acell ("+str(self.units['Distance'][0])+")",name = self.file1.getNameFile()+" Acell")
            self.GraphAcell.addLegend([r'$a$',r'$b$',r'$c$'])
        except:
            pass

        try:        
            self.GraphAngles
            angles = self.file1.getAngles()  * self.units['Angle'][1]
            self.GraphAngles.update(x,angles,'Step',"Angles ("+str(self.units['Angle'][0])+")",name = self.file1.getNameFile()+" Acell")
            self.GraphAngles.addLegend([r'$\alpha$',r'$\beta$',r'$\gamma$'])
        except:
            pass

        try:
#            if(units):
#                self.pos.units = self.units
#                self.pos.updateGraph()
#            else:
            self.pos.file = self.file1
            self.changeMode()
        except:
            pass 
                            
        try:
            self.msd.update(self.file1)
        except:
            pass

        try:
            self.GraphVAF
            #del self.vacf
            if (units == False):
                #The velocity autocorelation will not be update if the units change
                self.vacf = Analysis.Correlation(self.file1.getVel()).getCorrelationFunction(normalize = True)
                x = linspace(self.ni,self.nf-1,len(self.vacf))# Temporarily !!!
                self.GraphVAF.update(x,self.vacf,'step',"VAF",name = self.file1.getNameFile()+" VAF")        
        except:
            pass 
        

    def getData(self):
        data = ""
        data += self.lbl1.text()  + " " + self.lblnbatom.text()  +"\n"
        data += self.lbl2.text()  + " " + self.lblNbTime.text()  +"\n"
        data += self.lbl3.text()  + " " + self.lblE_tot.text()   +"\n"
        data += self.lbl4.text()  + " " + self.lblVol.text()     +"\n"
        data += self.lbl5.text()  + " " + self.lblTemp.text()    +"\n"
        data += self.lbl6.text()  + " " + self.lblPress.text()   +"\n"
        data += self.lblni.text() + " " + str(self.sbni.value()) +"\n"
        data += self.lblnf.text() + " " + str(self.sbnf.value()) +"\n"
        return data
                
    def keyPressEvent(self, event):
        if ((event.key() == 16777220 or event.key() == 16777221)):
            self.updateStep()
            #self.closeGraphic()
            event.accept()
            
    def displayData(self,data,deviation):
        return str("%.4e" %data) + ' +/- ' + str("%.1e" %deviation)
        # if(abs(data) < 1e-03 and deviation < 1e-03):
        #     return str("%.3g" %data) + ' +/- ' + str("%.2g" %deviation)
        # if(abs(data) < 1e-03 and deviation > 1e-03):
        #     return str("%.4g" %data) + ' +/- ' + str("%.2g" %deviation)
        # if(abs(data) > 1e+03 and deviation < 1e-03 ):
        #     return str("%.4g" %data) + ' +/- ' + str("%.2g" %deviation)
        # if(abs(data) > 1e+03 and deviation > 1e-03 ):
        #     return str("%.4g" %data) + ' +/- ' + str(round(deviation,2))
        # if(abs(data) > 1e-03 and deviation < 1e-03 ):
        #     return str(round(data,2)) + ' +/- ' + str("%.2g" %deviation)
        # if(abs(data) > 1e-03 and deviation > 1e-03 ):
        #     return str(round(data,2)) + ' +/- ' + str(round(deviation,2))
        # return ''	
