#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import sys,os,commands,time
import string, math
from numpy import array,sqrt,zeros,conjugate,arange,linspace,exp,log,sin,dot,degrees,arccos
from numpy.fft import rfft,fft, ifft
from numpy.linalg import inv

#Fortran Code
import fortran.math as math
import fortran.topology as topo

#----------------------------------------------------------------#
#--------------(AUTO)CORRELATION FUNCTION------------------------#
#----------------------------------------------------------------#
class Correlation:

   #-------------Constructor--------------#
    def __init__(self, serie1,serie2=None):

        # This class aim is to calculate the (auto)correlation function
        # the parameter is :
        #    - pSerie witch is the serie  ([time,atom,[x,y,z]] )

        self.nbtime = len(serie1)               # Number of steps in the simulation
        self.nbatom = len(serie1[0])               # Number of atoms in the simulation
        self.CorrelationFunction = zeros(2*self.nbtime,complex) # (Auto)correlation function array
        self.calculation(serie1,serie2)

   #------------Methods-----------------------#
    def calculation(self, pserie1, pserie2):

        serie1 = pserie1

        if pserie2 == None:
            serie2 = pserie1
        else:
            serie2 = pserie2

        for at in range(self.nbatom):

            for j in range(3):

                self.CorrelationFunction +=  ifft(conjugate(fft(serie1[:,at,j],2*self.nbtime,0))*fft(serie2[:,at,j],2*self.nbtime,0))

        #Normalisation by the therm 1/ ((nbtime - t) * 3 * N). N is the number of atoms, nbtime the number of step and t the step.
        self.CorrelationFunction = self.CorrelationFunction.real[:self.nbtime] / ((self.nbtime-arange(self.nbtime)) * 3 * self.nbatom)


    def getCorrelationFunction(self,normalize = False):
        if len(self.CorrelationFunction ) >= 0 :
            if normalize:
                return self.CorrelationFunction/self.CorrelationFunction[0]
            else:
                return self.CorrelationFunction
        else:
            return 0



#----------------------------------------------------------------#
#------------------DENSITY OF STATES (PHONONS)-------------------#
#----------------------------------------------------------------#
class DOS:
    def __init__(self, pdata, pres, pdtion):

        # This class aim is to calculate the density of phonon states (DOS),as being the Fourier transforms of the VACF
        # the parameters are  :
        #    - pdata witch is the array (1D) of the VACF
        #    - pres witch is the resolution
        #    - pdtion witch is ion time steps in atomic units of time
        self.atu = 2.41884*10**-5          # one atomic time unit (in ps)
        self.dtion = pdtion                   # ion time steps in atomic units of time
        self.VACF = pdata                # Velocity autocorelation Function (array 1D)
        self.dt = self.dtion * self.atu        # time step (in ps)
        self.nbtime = len(pdata)        # Number of steps
        self.times = self.dt * arange(self.nbtime)# 1D array with the times of the simulation
        self.fact = 0.6582054674822119  # conversion factor ps => 1/meV


        
        #Calculation of sigma for the gaussian function
        self.resolution = pres
        self.sigma = 1.0/(1.5192669*self.resolution/(2.0*sqrt(2.0*log(2.0))))



    def getDOS(self):
        #Return array (1D) phonons density of states (1/mev)
        DOS = self.windowGaussian(self.VACF, self.times, self.sigma)
        return fft(DOS,len(DOS),0).real[:len(self.VACF)] * self.dt * self.fact


    def getFrequencies(self):
        #Return array (1D) with the frenquencies (mev)
        frequencies = 4.1356 * arange(self.nbtime)/(2 * self.nbtime*self.dt)
        return frequencies


    def windowGaussian(self, inputSeries, x, s = 50, x0 = 0):
        # Return the input serie multiplying by the gaussian function
        # The size of the output array is double that the input serie,
        # because the final array is periodized function


        gauss= exp(-0.5*( (x - x0)**2 / (2*s**2)))


        # Creation of the final array:
        res = zeros(2*len(inputSeries))

        # Multiplying the input serie by the gaussian:
        win = gauss*inputSeries

        # Periodic function is creating :
        res[:len(inputSeries)] = win
        res[len(inputSeries):] = win[-1:-len(inputSeries)-1:-1]


        #old algorithm:
        #res = zeros((2*len(inputSeries) - 2,))
        #win = gauss*inputSeries
        #res[:len(inputSeries)] = win
        #res[len(inputSeries):] = win[-2:0:-1]

        return res


    def getSprectumDOS(self):
        return  rfft(self.VACF,2*len(self.VACF)-1,0).real * 2 * self.dt * self.fact




#----------------------------------------------------------------#
#-------------------Periodic table of the Element----------------#
#----------------------------------------------------------------#


class PeriodicTableElement :

    # This class returns the atomic number of an element.
    # Or the element corresponding to the atomic number.

    def __init__(self):
        self.table1 = { 'Ru' : 44, 'Re' : 75, 'Rf' : 104, 'Rg' : 111, 'Ra' : 88, 'Rb' : 37, 'Rn' : 86, 'Rh' : 45, 'Be' : 4, 'Ba' : 56,\
                        'Bh' : 107, 'Bi' : 83, 'Bk' : 97, 'Br' : 35, 'Ho' : 67, 'H' : 1, 'P' : 15, 'Os' : 76, 'Es' : 99, 'Hg' : 80,\
                        'Ge' : 32, 'Gd' : 64, 'Ga' : 31, 'He' : 2, 'Pr' : 59, 'Pt' : 78, 'Pu' : 94, 'Mg' : 12, 'Pb' : 82, 'Pa' : 91,\
                        'Pd' : 46, 'Xe' : 54, 'Po' : 84, 'Pm' : 61, 'Uut' : 113, 'Uuq' : 114, 'Uup' : 115, 'Uus' : 117, 'Uuo' : 118,\
                        'Uuh' : 116, 'Hf' : 72, 'K' : 19, 'Uub' : 112, 'Md' : 101, 'C' : 6, 'Mo' : 42, 'Mn' : 25, 'O' : 8, 'Mt' : 109,\
                        'S' : 16, 'W' : 74, 'Zn' : 30, 'Eu' : 63, 'Zr' : 40, 'Er' : 68, 'Ni' : 28, 'No' : 102, 'Na' : 11, 'Nb' : 41,\
                        'Nd' : 60, 'Ne' : 10, 'Np' : 93, 'Fr' : 87, 'Fe' : 26, 'Fm' : 100, 'B' : 5, 'F' : 9, 'Sr' : 38, 'N' : 7, 'Kr' : 36,\
                        'Si' : 14, 'Sn' : 50, 'Sm' : 62, 'V' : 23, 'Sc' : 21, 'Sb' : 51, 'Sg' : 106, 'Se' : 34, 'Co' : 27, 'Cm' : 96,\
                        'Cl' : 17, 'Ca' : 20, 'Cf' : 98, 'Ce' : 58, 'Cd' : 48, 'Lu' : 71, 'Cs' : 55, 'Cr' : 24, 'Cu' : 29, 'La' : 57,\
                        'Li' : 3, 'Tl' : 81, 'Tm' : 69, 'Lr' : 103, 'Th' : 90, 'Ti' : 22, 'Te' : 52, 'Tb' : 65, 'Tc' : 43, 'Ta' : 73,\
                        'Yb' : 70, 'Db' : 105, 'Dy' : 66, 'Ds' : 110, 'I' : 53, 'In' : 49, 'U' : 92, 'Y' : 39, 'Ac' : 89, 'Ag' : 47,\
                        'Ir' : 77, 'Am' : 95, 'Al' : 13, 'As' : 33, 'Ar' : 18, 'Au' : 79, 'At' : 85, 'Hs' : 108}

        self.table2 = { 44 : 'Ru', 75 : 'Re', 104 : 'Rf', 111 : 'Rg', 88 : 'Ra', 37 : 'Rb', 86 : 'Rn', 45 : 'Rh', 4 : 'Be', 56 : 'Ba',\
                        107 : 'Bh', 83 : 'Bi', 97 : 'Bk', 35 : 'Br', 67 : 'Ho', 1 : 'H', 15 : 'P', 76 : 'Os', 99 : 'Es', 80 : 'Hg',\
                        32 : 'Ge', 64 : 'Gd', 31 : 'Ga', 2 : 'He', 59 : 'Pr', 78 : 'Pt', 94 : 'Pu', 12 : 'Mg', 82 : 'Pb', 91 : 'Pa',\
                        46 : 'Pd', 54 : 'Xe', 84 : 'Po', 61 : 'Pm', 113 : 'Uut', 114 : 'Uuq', 115 : 'Uup', 117 : 'Uus', 118 : 'Uuo',\
                        116 : 'Uuh', 72 : 'Hf', 19 : 'K', 112 : 'Uub', 101 : 'Md', 6 : 'C', 42 : 'Mo', 25 : 'Mn', 8 : 'O', 109 : 'Mt',\
                        16 : 'S', 74 : 'W', 30 : 'Zn', 63 : 'Eu', 40 : 'Zr', 68 : 'Er', 28 : 'Ni', 102 : 'No', 11 : 'Na', 41 : 'Nb',\
                        60 : 'Nd', 10 : 'Ne', 93 : 'Np', 87 : 'Fr', 26 : 'Fe', 100 : 'Fm', 5 : 'B', 9 : 'F', 38 : 'Sr', 7 : 'N',\
                        36 : 'Kr', 14 : 'Si', 50 : 'Sn', 62 : 'Sm', 23 : 'V', 21 : 'Sc', 51 : 'Sb', 106 : 'Sg', 34 : 'Se', 27 : 'Co',\
                        96 : 'Cm', 17 : 'Cl', 20 : 'Ca', 98 : 'Cf', 58 : 'Ce', 48 : 'Cd', 71 : 'Lu', 55 : 'Cs', 24 : 'Cr', 29 : 'Cu',\
                        57 : 'La', 3 : 'Li', 81 : 'Tl', 69 : 'Tm', 103 : 'Lr', 90 : 'Th', 22 : 'Ti', 52 : 'Te', 65 : 'Tb', 43 : 'Tc',\
                        73 : 'Ta', 70 : 'Yb', 105 : 'Db', 66 : 'Dy', 110 : 'Ds', 53 : 'I', 49 : 'In', 92 : 'U', 39 : 'Y', 89 : 'Ac',\
                        47 : 'Ag', 77 : 'Ir', 95 : 'Am', 13 : 'Al', 33 : 'As', 18 : 'Ar', 79 : 'Au', 85 : 'At', 108 : 'Hs'}

        self.table3 = {   1:    1.007940,  2:    4.002602,  3:    6.941000,  4:    9.012182,  5:   10.811000,  6:   12.011000,  7:   14.006740,\
                          8:   15.999400,  9:   18.998404, 10:   20.179701, 11:   22.989767, 12:   24.305000, 13:   26.981539, 14:   28.085501,\
                          15:   30.973763, 16:   32.066002, 17:   35.452702, 18:   39.948002, 19:   39.098301, 20:   40.077999, 21:   44.955910,\
                          22:   47.880001, 23:   50.941502, 24:   51.996101, 25:   54.938049, 26:   55.847000, 27:   58.933201, 28:   58.689999, 29:   63.546001, 30:   65.389999,\
                          31:   69.723000, 32:   72.610001, 33:   74.921593, 34:   78.959999, 35:   79.903999, 36:   83.800003, 37:   85.467796, 38:   87.620003, 39:   88.905853,\
                          40:   91.223999, 41:   92.906380, 42:   95.940002, 43:   98.906197, 44:  101.070000, 45:  102.905502, 46:  106.419998, 47:  107.868202, 48:  112.411003,\
                          49:  114.820000, 50:  118.709999, 51:  121.752998, 52:  127.599998, 53:  126.904472, 54:  131.289993, 55:  132.905426, 56:  137.326996, 57:  138.905502,\
                          58:  140.115005, 59:  140.907654, 60:  144.240005, 61:  147.910004, 62:  150.360001, 63:  151.964996, 64:  157.250000, 65:  158.925339, 66:  162.500000, \
                          67:  164.930313, 68:  167.259995, 69:  168.934204, 70:  173.039993, 71:  174.966995, 72:  178.490005, 73:  180.947906, 74:  183.850006, 75:  186.207001, \
                          76:  190.199997, 77:  192.220001, 78:  195.080002, 79:  196.966537, 80:  200.589996, 81:  204.383301, 82:  207.199997, 83:  208.980377, 84:  209.000000, \
                          85:  210.000000, 86:  222.000000, 87:  223.000000, 88:  226.025406, 89:  230.000000, 90:  232.038101, 91:  231.035904, 92:  238.028900, 93:  237.048203, \
                          94:  242.000000, 95:  243.000000, 96:  247.000000, 97:  247.000000, 98:  249.000000, 99:  254.000000,100:  253.000000,101:  256.000000,102:  254.000000,\
                          103:  257.000000}

    def getZnucl(self,name):
        try :
            return self.table1[name]
        except:
            return 0

    def getName(self,znucl):
        try :
            return self.table2[znucl]
        except:
            return 'Unknown Element'

    def getMass(self,znucl):
        try :
            return self.table3[znucl]
        except:
            return 'Unknown Element'


#----------------------------------------------------------------#
#-------------------RADIAL DISTRIBUTION FUNCTION-----------------#
#----------------------------------------------------------------#
class RDF:

    def __init__(self, pfile, mode, atom1, atom2, box, pdr, pstep, No_Windows=False):

        # This class aim is to calculate the Radial Distrubution Function
        # the parameters are :
        #    - mode = 0 for normal RDF, 1 for its deconvolution
        #    - pfile  is the output file (ascii/netcdf)
        #    - atom1/2 give the number of the atom in typat
        #    - box : give the number of box wide for the calculculation of g(r)
        #    - pdr   : give the radial precision of the calculation
        #    - pstep : give the self.increment for the step loop

        file = pfile

        pos = file.getXCart()          # position of atom (cartesian)

        inc = pstep                    # incrementation

        acell = file.getAcell()        # acell
        a = acell[:,0]
        b = acell[:,1]
        c = acell[:,2]
        a_max = max(acell[:,0])
        b_max = max(acell[:,1])
        c_max = max(acell[:,2])

        rprim = zeros((3,3))           # primitives vectors of the cell
        rprim[0] = file.getRPrim()[0,0]
        rprim[1] = file.getRPrim()[0,1]
        rprim[2] = file.getRPrim()[0,2]

        rho = 1/file.getVol()        # density

        if box >= 0:

            f = (1/2.+box)           # factor of maximum radius
            
        else: # box < 0
            
            f = 3**0.5/2.

        rmax = f*min(a_max,b_max,c_max)

        deltaR = pdr                      # delta r  (Bohr)

        maxbin = int((rmax/deltaR))       # number of iteration

        typat = file.getTypat()           # Type of particule

        if atom1 == 0:
            indexAtom1 = array([i for i,x in enumerate(typat) if x == 1 and 2],dtype=int) + 1
        else:
            indexAtom1 = array([i for i,x in enumerate(typat) if x == atom1],dtype=int) + 1 #get list of the index of typat1 (+1 for fortran)

            
        if atom2 == 0:
            indexAtom2 = array([i for i,x in enumerate(typat) if x == 1 and 2 ],dtype=int) + 1
        else:
            indexAtom2 = array([i for i,x in enumerate(typat) if x == atom2],dtype=int) + 1 #get list of the index of typat2 (+1 for fortran)
                
        nbtime =  len(pos)                    # final step

        self.r = arange(maxbin)*deltaR

        if mode == 0:

            self.g = topo.pair_distribution(nbtime,pos,rprim,f,inc,a,b,c,deltaR,rho,indexAtom1,indexAtom2,self.r)

            self.r += deltaR/2.

        else: #mode == 1

            nba = len(pos[0])                   # number of atom

            it = (int(nbtime/inc)+1)*nba*(nba-1)

            self.m = topo.m_distribution(0,nbtime,pos,inc,it,a,b,c,indexAtom1,indexAtom2)

    def getR(self):
        return self.r

    def getRDF(self):
        return self.g[0]

    def getINT(self):
        return self.g[1]

    def getDATA(self):
        return self.m[2]



class DEC:

    def __init__(self, pfile, n, atom1, atom2, data, box, pdr, step, No_Windows=False):

        # This class aim is to calculate the Deconvolution of the RDF
        # the parameters are :
        #    - pfile  is the output file (ascii/netcdf)
        #    - n = is the neibhor number
        #    - atom1/2 give the number of the atom in typat
        #    - box : give the number of box wide for the calculculation of g(r)
        #    - pdr   : give the radial precision of the calculation
        #    - pstep : give the self.increment for the step loop

        file = pfile

        pos = file.getXCart()            # position of atom (cartesian)

        nei = n

        inc = step                    #incrementation

        acell = file.getAcell()        #acell
        a_max = max(acell[:,0])
        b_max = max(acell[:,1])
        c_max = max(acell[:,2])

        if box >= 0:

            f = (1/2.+ box)           # factor of maximum radius

        else: # box < 0
            
            f = 3**0.5/2.

        rmax = f*min(a_max,b_max,c_max)

        deltaR = pdr                           # delta r  (Bohr)

        maxbin = int((rmax/deltaR))            # number of iteration

        typat = file.getTypat()              # Type of particule

        if atom1 == 0:
            indexAtom1 = array([i for i,x in enumerate(typat) if x == 1 and 2],dtype=int) + 1
        else:
            indexAtom1 = array([i for i,x in enumerate(typat) if x == atom1],dtype=int) + 1 #get list of the index of typat1 (+1 for fortran)

            
        if atom2 == 0:
            indexAtom2 = array([i for i,x in enumerate(typat) if x == 1 and 2 ],dtype=int) + 1
        else:
            indexAtom2 = array([i for i,x in enumerate(typat) if x == atom2],dtype=int) + 1 #get list of the index of typat2 (+1 for fortran)
                
        nbtime =  len(pos)                    # final step

        self.r = arange(maxbin)*deltaR

        self.n = topo.n_distribution(nei,nbtime,data,inc,deltaR,indexAtom1,indexAtom2,0,self.r)        

    def getNEI(self):
        return self.n

#-----------------------------------------------------------------#
#-------------------ANGULAR DISTRIBUTION FUNCTION-----------------#
#-----------------------------------------------------------------#

class ADF:

    def __init__(self, nei, pfile, atom1, atom2, pdt, pstep, No_Windows=False):

        # This class aim is to calculate the Angular Distrubution Function
        # the parameters are :
        #    - nei gives the number of neibours we consider
        #    - pfile  is the output file (ascii/netcdf)
        #    - atom1/2 give the number of the atom in typat
        #    - pdt   : give the dtheta for the calculation
        #    - pstep : give the self.increment for the step loop

        file = pfile

        pos = file.getXCart()            # position of atom (cartesian)

        inc = pstep                    #incrementation

        acell = file.getAcell()        #acell
        a = acell[:,0]
        b = acell[:,1]
        c = acell[:,2]

        rprim = zeros((3,3))           # primitives vectors of the cell
        rprim[0] = file.getRPrim()[0,0]
        rprim[1] = file.getRPrim()[0,1]
        rprim[2] = file.getRPrim()[0,2]

        deltaTheta = pdt                           # delta Theta  (degres)

        maxbin = int((180./deltaTheta)) + 1            # number of iteration

        typat = file.getTypat()              # Type of particule

        if atom1 == 0:
            indexAtom1 = array([i for i,x in enumerate(typat) if x == 1 and 2],dtype=int) + 1
        else:
            indexAtom1 = array([i for i,x in enumerate(typat) if x == atom1],dtype=int) + 1 #get list of the index of typat1 (+1 for fortran)


        if atom2 == 0:
            indexAtom2 = array([i for i,x in enumerate(typat) if x == 1 and 2 ],dtype=int) + 1
        else:
            indexAtom2 = array([i for i,x in enumerate(typat) if x == atom2],dtype=int) + 1 #get list of the index of typat2 (+1 for fortran)            

        nbtime =  len(pos)                    # final step

        self.t = arange(maxbin)*deltaTheta

        self.a = topo.angular_distribution(nei,nbtime,pos,rprim,inc,a,b,c,deltaTheta,indexAtom1,indexAtom2,self.t)

    def getTheta(self):
        return self.t

    def getADF(self):
        return self.a


#------------------------------------------------------------------#
#-------------------NEIGHBOR DISTRIBUTION FUNCTION-----------------#
#------------------------------------------------------------------#

class NDF:

    def __init__(self, nei, pfile, atom1, atom2, pdr, pstep, No_Windows=False):

        # This class aim is to calculate the Neighbor Distrubution Function
        # the parameters are :
        #    - nei gives the neith neibours we consider.
        #    - pfile  is the output file (ascii/netcdf)
        #    - atom1 give the number of the atom in typat
        #    - pdt   : give the dtheta for the calculation
        #    - pstep : give the self.increment for the step loop

        file = pfile

        pos = file.getXCart()            # position of atom (cartesian)

        inc = pstep                    #incrementation

        acell = file.getAcell()        #acell
        a = acell[:,0]
        b = acell[:,1]
        c = acell[:,2]

        rprim = zeros((3,3))           # primitives vectors of the cell
        rprim[0] = file.getRPrim()[0,0]
        rprim[1] = file.getRPrim()[0,1]
        rprim[2] = file.getRPrim()[0,2]

        deltaR = pdr

        typat = file.getTypat()              # Type of particule

        if atom1 == 0:
            indexAtom1 = array([i for i,x in enumerate(typat) if x == 1 and 2],dtype=int) + 1
        else:
            indexAtom1 = array([i for i,x in enumerate(typat) if x == atom1],dtype=int) + 1 #get list of the index of typat1 (+1 for fortran)


        if atom2 == 0:
            indexAtom2 = array([i for i,x in enumerate(typat) if x == 1 and 2 ],dtype=int) + 1
        else:
            indexAtom2 = array([i for i,x in enumerate(typat) if x == atom2],dtype=int) + 1 #get list of the index of typat2 (+1 for fortran)

        nbtime =  len(pos)                    # final step

        nba = len(pos[0])                     # number of atom

        it = (int(nbtime/inc)+1)*nba*(nba-1)

        self.m = topo.m_distribution(nei,nbtime,pos,rprim,inc,it,a,b,c,indexAtom1,indexAtom2)

        maxi = self.m[0]

        mini = self.m[1]

        data = self.m[2]

        self.r = []

        for i in range(int(mini*1000),int(maxi*1000+1),int(deltaR*1000)):

            i = i/1000.

            self.r.append(i) 

        self.n = topo.n_distribution(nei,nbtime,data,inc,deltaR,indexAtom1,indexAtom2,mini,self.r)

    def getR(self):
        return self.r

    def getNDF(self):
        return self.n


#--------------------------------------------------------#
#-------------------Probability function-----------------#
#--------------------------------------------------------#

class Proba:

    def __init__(self, nei, pfile, atom1, No_Windows=False):

        # This class aim is to calculate the Pobability that one neighbor stay in the same neighborhood of an atom between the first step and the last step considered
        # the parameters are :
        #    - nei give the number of neighbors considered
        #    - pfile  is the output file (ascii/netcdf)
        #    - atom1 give the number of the atom in typat

        file = pfile

        pos = file.getXCart()            # position of atom (cartesian)

        acell = file.getAcell()        #acell
        a = acell[:,0]
        b = acell[:,1]
        c = acell[:,2]

        rprim = zeros((3,3))           # primitives vectors of the cell
        rprim[0] = file.getRPrim()[0,0]
        rprim[1] = file.getRPrim()[0,1]
        rprim[2] = file.getRPrim()[0,2]

        typat = file.getTypat()              # Type of particule

        if atom1 == 0:
            indexAtom1 = array([i for i,x in enumerate(typat) if x == 1 and 2],dtype=int) + 1
        else:
            indexAtom1 = array([i for i,x in enumerate(typat) if x == atom1],dtype=int) + 1 #get list of the index of typat1 (+1 for fortran)
            
        nbtime =  len(pos)                    # final step

        self.n = arange(1,nei+1)

        self.p = topo.probability(nei,nbtime,pos,rprim,a,b,c,indexAtom1,indexAtom1)

    def getN(self):
        return self.n

    def getProba(self):
        return self.p




#-----------------------------------------------------------------#
#-------------------MEANS SQUARE DISPLACEMENT---------------------#
#-----------------------------------------------------------------#
class MSD:

    def __init__(self, pfile, atom1):

        # This class aim is to calculate the Mean square displacement
        # the parameters are :
        #    - pfile  is the output file (ascii/netcdf)
        #    - atom1 give the number of the atom in typat

        file = pfile

        pos = file.getXCart()          # position of atom (cartesian)

        typat = file.getTypat()        # Type of particule
        
        indexAtom1 = array([i for i,x in enumerate(typat) if x == atom1],dtype=int) + 1 #get list of the index of typat1 (+1 for fortran)

        self.mds =  math.mean_square_displacement(pos,indexAtom1)
        self.x = arange(len(self.mds))

    def getX(self):
        return self.x

    def getMSD(self):
        return self.mds






#---------------------------------------------------#
#---------ELASTIC CALCULATION CLASS-----------------#
#---------------------------------------------------#
class elasticCalculation:

   #-------------Constructor--------------#
    def __init__(self,pfiles,pni,pWDData):
        # This class aim is to calculate the Elastic constant
        # the parameter is :
        #        -pfiles : witch is dictionnay array with the path of the files

        self.files   = pfiles         #Dictionnary array
        self.elastic = {}         #Dictionnary array with the elastics constants
        self.ni = pni                 #Departure of the step
        self.datasetWD = pWDData #Dataset without deformation

        self.struct      = structure(self.files).getStructure()

        self.deformation = deformation(self.files,self.ni,self.datasetWD).getDeformation()

        if len(self.deformation) != 0 :
            if (self.struct == 'cubic'):
                self.elastic['C11'] = 0
                self.elastic['C12'] = 0
                self.elastic['C44'] = 0

                self.M = 0 # M = C11 - C12  using in the tetragonal deformotation
                self.T = 0 # T = C11 + 2C12 using in the isotropic deformotation

                #print self.deformation
                #for key in self.deformation:
                        #print key

                sigma_0 = self.deformation['WD'][2]

                for i in range(len(self.deformation)-1):

                    sigma_i = self.deformation[i][2]
                    delta = self.deformation[i][1]

                    if self.deformation[i][0] == 'single':

                        C11 = ( sigma_i[2] - sigma_0[2] ) / delta
                        C12 = ( sigma_i[1] - sigma_0[1] ) / delta
                        C44 = ( sigma_i[5] - sigma_0[5] ) / (2*delta)

                        self.elastic['C11'] += C11
                        self.elastic['C12'] += C12
                        self.elastic['C44'] += C44

                    if self.deformation[i][0] == 'orto':
                        C44 = ( sigma_i[5] - sigma_0[5] ) / (2*delta)

                        self.elastic['C44'] += C44

                    if self.deformation[i][0] == 'tetra':
                        self.M += ( sigma_i[1]-  sigma_0[1] ) / delta

                    if self.deformation[i][0] == 'iso':

                        self.T += ( sigma_i[1] - sigma_0[1] ) / delta

                     #----------ADD OTHER DEFORMATION HERE--------#

                     #--------------------------------------------#
                try:
                    self.M /= self.count(self.deformation,'tetra')
                    self.T /= self.count(self.deformation,'iso')
                except:
                    pass

                nbOfSingleDef = self.count(self.deformation,'single')


                if (self.M == 0 and self.T != 0 or self.M != 0 and self.T == 0):
                    print 'Deformation error 1'
                    self.elastic = {}
                    return

                elif (self.M != 0 and self.T != 0):

                    self.elastic['C11'] +=   ( 2 * self.M + self.T  ) / 3
                    self.elastic['C12'] +=  -( self.M - self.T  ) / 3

                    nbC11 = nbOfSingleDef +1

                else :
                    nbC11 = nbOfSingleDef



                if nbC11  != 0 :

                    self.elastic['C11'] /= nbC11
                    self.elastic['C12'] /= nbC11

                self.elastic['C44'] /= (nbOfSingleDef + self.count(self.deformation,'orto')  )


                return


                #----------ADD OTHER STRUCTURE HERE--------#
                if (self.struct == 'hexagonal'):

                    return
                #------------------------------------------#
        else:
            self.elastic = {}
            print 'Deformation error 2'
            return
    #------------Methods--------------------#


    def getElasticsConstants(self):
        return self.elastic


    def count(self,dic,pkey):
        i = 0

        for key in dic:
            if  dic[key][0] == pkey :
                i += 1
        return i

#---------------------------------------------------#
#---------------------------------------------------#
#----------------STRUCTURE CLASS--------------------#
#--------------------NOT YET IMPLEMENTED------------#
#---------------------------------------------------#
class structure:

    #-------------Constructor--------------#
    def __init__(self, pfiles):
        # This class aim is to determinate the structure of the system
        # the parameters are :
        #        -self.files
        self.files       = pfiles

    def getStructure(self):
        return "cubic"

#---------------------------------------------------#
#---------------------------------------------------#
#----------------DEFORMATION CLASS------------------#
#---------------------------------------------------#
#---------------------------------------------------#
class deformation:

    #-------------Constructor--------------#
    def __init__(self, pfiles,pni,pWDData):
        # This class aim is to determinate the deformations used for the calculation of elastic constant
        # the parameters are :
        #        -self.files
        #        -self.ni

        self.files       = pfiles
        self.ni          = pni
        self.datasetWD   = pWDData #Dataset without deformation
        self.deformation = {}
        self.goodDeformation = True

        self.calculation()

    #------------Methods--------------------#
    def calculation(self):
        self.rprim0 = []
        self.sigma = {}
        self.deformation = {}

        if self.files['WD'].getIonMove()[0] in [1,6,7,8,9,12] :

            self.files['WD'].setNi(self.ni)
            sigma = self.files['WD'].getStress()
            sigma = sum(sigma[:,:]) / len(sigma)
            self.deformation['WD'] = ['WD', 0 ,sigma]

            self.rprim_0 = self.files['WD'].getRPrim()[0,:,:]

            for i in range(len(self.files)-1):
                self.files[i].setNi(self.ni)
                sigma = self.files[i].getStress()
                sigma = sum(sigma[:,:]) / len(sigma)

                rprim = self.files[i].getRPrim()[0,:,:]
                d = dot(inv(self.rprim_0),rprim)
                self.deform(i,d,sigma)
            return

        if self.files['WD'].getIonMove()[0] == 0 :
            if 1 == 1:
                sigma = self.files['WD'].getStress()[self.datasetWD]
                self.deformation['WD'] = ['WD', 0 ,sigma]
                self.rprim_0 = self.files['WD'].getRPrim()[self.datasetWD]


                self.nbdataset = self.files['WD'].getNbDataset()
                j = 0

                for i in range(self.nbdataset) :

                    if i != self.datasetWD :
                        sigma = self.files['WD'].getStress()[i]
                        rprim = self.files['WD'].getRPrim()[i]
                        d = dot(inv(self.rprim_0),rprim)
                        self.deform(j,d,sigma)
                        j += 1

                for i in range(len(self.files)-1):
                    self.nbdataset = self.files[i].getNbDataset()
                    for j in range(self.nbdataset) :

                        sigma = self.files[i].getStress()[j]
                        rprim = self.files[i].getRPrim()[j]
                        d = dot(inv(self.rprim_0),rprim)
                        self.deform(i,d,sigma)
                return
        else :
            self.goodDeformation = False

    def deform(self,i,d,sigma):

        if( round(d[0][1],6) == round(d[1][0],6) and round(d[0][0],6) == round(d[1][1],6) and round(d[2][2],6) == (1 + round(d[0][1],6)) ):
            self.deformation[i] = ['single',round(d[0][1],5),sigma]
            return


        if( round(d[0][0],6) == round(d[1][1],6) and round(d[2][2],4) ==  round( 1.0 / ( d[0][0]**2 ),4 ) ):\

            if ( abs(d[0][1]) <= 1e-10 and abs(d[1][0]) <= 1e-10 and abs(d[2][0]) <= 1e-10 and  d[2][1]<= 1e-10 and abs(d[1][0]) <= 1e-10 ):
                self.deformation[i] = ['tetra',d[0][0]-1,sigma]
                return


        if( round(d[0][1],6) == round(d[1][0],6) and round(d[2][2],4) == round( 1.0 / (1 - d[0][1]**2 ),4 ) ):
            if (d[0][0] == d[1][1] and abs(d[0][2]) <= 1e-10 and abs(d[1][2]) <= 1e-10):
                self.deformation[i] = ['orto',round(d[0][1],4),sigma]
                return

        if( round(d[0][0],6) == round(d[1][1],6) and round(d[1][1],6) == round(d[2][2],6) and  round(d[0][0],6) != 0 ):
            if ( abs(d[0][1]) <= 1e-10 and  abs(d[0][2]) <= 1e-10 and  abs(d[0][1]) <= 1e-10 and  abs(d[1][0]) <= 1e-10):
                if (abs(d[1][0]) <= 1e-10 and  abs(d[1][2])<= 1e-10 and  abs(d[1][0]) <= 1e-10 ):
                    self.deformation[i] = ['iso',d[0][0]-1,sigma]
            return

        #----------ADD OTHER DEFORMATION HERE------#

        #------------------------------------------#
        self.goodDeformation = False

    #--------Accessor-----------#
    def getDeformation(self):
        if (self.goodDeformation):
            return self.deformation
        else:
            return []



#----------------------------------------------------------------#
#-------------------CALCULATION OF BULK MODULUS------------------#
#----------------------------------------------------------------#

class birch:


    def __init__(self,data,number):
        #This class uses FORTRAN code for the fit
        self.value = {}

        self.namefiles ="fit.inp"

        self.workDir = os.getcwd()

        self.homeDir = commands.getoutput('echo ~/')


        path = os.path.dirname(__file__) + "/" # Get the path of the source Code

        os.system('cp -r '+path + 'BIRCH '+ self.homeDir) # Copy the BIRCH DIR to the home DIR

        os.chdir(self.homeDir + 'BIRCH')  # Move to the BIRCH DIR
        print  os.getcwd()
        file_out = open( self.namefiles ,"w") # Creation of the inputFIle for the FIT
        file_out.write('  2      ! Nb of columns\n')
        file_out.write('  1      ! Column of X-axis\n')
        file_out.write('  2      ! Column of Y-axis\n')
        file_out.write('  1      ! What is X: 1 for Volume, 2 for Acell\n')
        file_out.write('  1      ! 1 for V=a^3/4, 2 for V=a^3/2\n')
        file_out.write('  2      ! Unit for L: 1 for angstrom, 2 for u.a.\n')
        file_out.write('  3      ! Unit for E: 1 for Ryd, 2 for Ev, 3 for Ha, 4 for ?\n')
        file_out.write('  1      ! Mupliplier for V (V will be multiplied by this factor)\n')
        file_out.write('  1      ! Mupliplier for E (E will be divided by this factor)\n')
        file_out.write('  '+str(number)+'     ! Nb of input points\n')
        file_out.write(data)
        file_out.write('  3     ! Order of v^(-2/3) polynomial\n')
        file_out.write(' -2     ! ?\n')
        file_out.write('  0     ! ?\n')
        file_out.write('  0     ! ?\n')
        file_out.write('  0     ! ?\n')
        file_out.close()

        benf = commands.getoutput('./benf ')  # execution of the FIT program
        if benf == '':
                print 'fit program computes successfully!'
        else:
            commands.getoutput('ifort -O -o benf benf.f')# compilation of the fit program
            benf = commands.getoutput('./benf ')
            if benf =='':
                print 'fit program computes successfully!'
            else:
                print 'unable to launch the fit program'


        ## Get the value of the fit
        self.value['E0']  = float( (commands.getoutput('grep -P \'Eo = .* Rydbergs\' fit.res').split())[2]) / 2.0
        self.value['B0']  = float( (commands.getoutput('grep -P \'Ko = .* GPa\' fit.res').split())[2] )
        self.value['dB0'] = float( (commands.getoutput('grep -P \"Ko\'=.*\" fit.res').split())[1] )
        self.value['V0']  = float( (commands.getoutput('grep -P \'Vo =.* Bohr.*\' fit.res').split())[2])
        self.value['RMS'] = float( (commands.getoutput('grep -P \"RMS.*Ry\" fit.res').split())[6]) / 2.0

        os.chdir(self.homeDir) # Move to the HOME DIR

        os.system('rm -r BIRCH') # Remove the BIRCH DIR

        os.chdir(self.workDir) # Move to the original path

    def getValue(self):
        return self.value


