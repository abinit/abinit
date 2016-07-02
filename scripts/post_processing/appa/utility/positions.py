#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
authors: Martin Alexandre, Vincent Stutzmann
last edited: May 2014
"""

from numpy import array,reshape,ones

#Fortran Code
import fortran.positions as posit

#----------------------------------------------------------------#
#-------------------CONVERSION POSITIONS-------------------------#
#----------------------------------------------------------------#

class convert_position:

    def __init__(self,pfile,indexAtom,option=0,reduced=False):

        # Option:
        #       0 return the positions of each atom (Default)  
        #       1 return the average positions 
        #       2 return the positions of each atom using periodic boundary condistions 
        # Reduced :
        #       True  return the result in reduced coordinates
        #       False return the result in cartesian coordinates

        self.file = pfile
        natom  = len(indexAtom)    # Number particule
        
        if (reduced):
            self.pos = self.file.getXRed()   # position of atom (cartesian)

        else:
            self.pos = self.file.getXCart()  # position of atom (reduced)

        nbtime = len(self.pos)        # number of step
        acell = self.file.getAcell()  # acell

        if(reduced):
            a = ones(nbtime)          # Set a,b and c to 1
            b = ones(nbtime)          # Only for cubic (must be generalize)
            c = ones(nbtime)
        else:
            a = acell[:,0]
            b = acell[:,1]
            c = acell[:,2]
    

        self.result = [[],[],[]]        

        if option==0:
            self.result = posit.atom_std(nbtime,self.pos,indexAtom+1)
        
        if option==1:
            self.result = posit.atom_ave(nbtime,self.pos,indexAtom+1)

        if option==2:
            self.result = posit.atom_sym(nbtime,self.pos,a,b,c,indexAtom+1)


    def getX(self):
        # Return array with x postion of 1 type of particule
        # This array is reshape by type x[ntime]
        return self.result[0]

    def getY(self):
        # Return array with y postion of 1 type of particule
        # This array is reshape by type y[ntime]
        return self.result[1]

    def getZ(self):
        # Return array with z postion of 1 type of particule
        # This array is reshape by type z[ntime]
        return self.result[2]
