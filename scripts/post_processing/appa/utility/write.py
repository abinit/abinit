#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import sys,os,commands
import string, math

#Utility
import utility.analysis as Analysis

from numpy import ndarray

class SaveFile:

    namefiles = ""

    def __init__(self, pnamefile):
        self.namefile = pnamefile
    
    def saveGraph(self,x,y,xlbl,ylbl):
        file_out = open( self.namefile ,"w")
        #if y have more than 1 quantity (stress for example)
        #if True:
        try:
            line ="#"+xlbl
            for i in range(len(y[0])):
                if(len(x)==len(y[:,i])):
                    line+="\t"+ylbl+str(i)
            line+='\n'
            file_out.write(line)
            for t in range (len(x)):
                line = str(x[t])
                for i in range(len(y[0])):
                    line +=" \t "+str(y[t][i])
                line += "\n"
                file_out.write(line)
                                                    
                                                        
        #if y have  1 quantity (total energy for example)
        #else:
        except:                        
            if(len(x)==len(y)):
                line ="#"+xlbl+"\t"+ylbl+'\n'
                file_out.write(line)
                for i in range (len(x)):
                    line = str(x[i])+" \t "+ str(y[i]) + "\n"
                    file_out.write(line)

                                    
        file_out.close()

    def saveData(self,data):
        file_out = open( self.namefile ,"w")
        if data != "":
            file_out.write(data) 
        file_out.close()
            
            
    def xyzFormat(self,pos,acell,typat,znucl):  
        PTOE = Analysis.PeriodicTableElement()        
        file_out = open( self.namefile ,"w")
        for t in range(len(pos)):
            file_out.write(str(len(pos[0]))+'\n')
            file_out.write('acell ' + str(acell[t][0])+' '+str(acell[t][1])+' '+str(acell[t][2])+' '+'\n')
            for i in range(len(pos[0])):
                typeofatom = PTOE.getName(znucl[typat[i] - 1 ])
                file_out.write(str(typeofatom) +' '+ str(pos[t][i][0]) +' '+ str(pos[t][i][1]) +' '+ str(pos[t][i][2]) + '\n')
