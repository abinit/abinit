#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: May 2013
"""

import sys,os,commands
import string, math, shutil
from scipy.io import netcdf
import numpy as np

class writeHIST:

    def __init__(self,pfile,pname,pni,pnf):
        self.namefile = pname
        self.file     = pfile
        self.ni       = pni
        self.nf       = pnf
        self.write()

    def write(self):
        
#       Copy _OUT.nc file
        name_out_nc = str(self.file.getNameFile(fullpath=True)).replace('_HIST','_OUT.nc')
        self.new_out_nc = str(self.namefile).replace('_HIST','_OUT.nc')
        shutil.copyfile(name_out_nc,self.new_out_nc)

#       Creation of new HIST file
        f = netcdf.netcdf_file(self.namefile, 'w')

#       dimensions:
        f.createDimension('natom',self.file.getNatom()) 
        f.createDimension('xyz',3) 
        f.createDimension('time',self.file.getNbStep())
        f.createDimension('six', 6) 

#       variables:
        xcart = f.createVariable('xcart','d',('time','natom','xyz'))
        xcart[:,:,:] = self.file.xcart[self.ni-1:self.nf,0]
        xcart.units = "bohr"
        xcart.mnemonics = "vectors (X) of atom positions in CARTesian coordinates" ;

        xred = f.createVariable('xred','d',('time','natom','xyz'))
        xred[:,:,:] = self.file.xred[self.ni-1:self.nf,0]
        xred.units = "dimensionless" ;
        xred.mnemonics = "vectors (X) of atom positions in REDuced coordinates" ;

        fcart = f.createVariable('fcart','d',('time','natom','xyz'))
        fcart[:,:,:] = self.file.fcart[self.ni-1:self.nf,0]
        fcart.units = "Ha/bohr" ;
        fcart.mnemonics = "atom Forces in CARTesian coordinates" ;

        fred = f.createVariable('fred','d',('time','natom','xyz'))
        fred[:,:,:] = self.file.fred[self.ni-1:self.nf,0]
        fred.units = "dimensionless" ;
        fred.mnemonics = "atom Forces in REDuced coordinates" ;

        vel = f.createVariable('vel','d',('time','natom','xyz'))
        vel[:,:,:] = self.file.vel[self.ni-1:self.nf,0]
        vel.units = "bohr*Ha/hbar" ;
        vel.mnemonics = "VELocity" ;

        acell = f.createVariable('acell','d',('time','xyz'))
        acell[:,:] = self.file.acell[self.ni-1:self.nf,0]
        acell.units = "bohr" ;
        acell.mnemonics = "CELL lattice vector scaling" ;

        rprimd = f.createVariable('rprimd','d',('time','xyz','xyz'))
        rprimd[:,:,:] = self.file.rprimd[self.ni-1:self.nf,0]
        rprimd.units = "bohr" ;
        rprimd.mnemonics = "Real space PRIMitive translations, Dimensional" ;

        etotal = f.createVariable('etotal','d',('time',))
        etotal[:] = self.file.E_pot[self.ni-1:self.nf]
        etotal.units = "Ha" ;
        etotal.mnemonics = "TOTAL Energy" ;

        ekin = f.createVariable('ekin','d',('time',))
        ekin[:] = self.file.E_kin_ion[self.ni-1:self.nf]
        ekin.units = "Ha" ;
        ekin.mnemonics = "Energy KINetic ionic" ;

        mdtime = f.createVariable('mdtime','d',('time',))
        mdtime[:] = self.file.mdtime[self.ni-1:self.nf]
        mdtime.units = "hbar/Ha" ;
        mdtime.mnemonics = "Molecular Dynamics TIME" ;

        strten = f.createVariable('strten','d',('time','six'))
        strten[:,:] = self.file.stress[self.ni-1:self.nf,0]
        strten.units = "Ha/bohr^3" ;
        strten.mnemonics = "STRess tensor" ;

        f.close()
