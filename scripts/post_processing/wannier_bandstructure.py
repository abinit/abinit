#!/usr/bin/python
#=================================================================#
#  Script to plot the bandstructure from an abinit bandstructure  #
#  _EIG.nc netcdf file or from a wannier bandstructure, or from   #
#  an _EIG.nc file+GW file+ bandstructure _EIG.nc file            #
#=================================================================#

#########
#IMPORTS#
#########

import numpy as N
import matplotlib.pyplot as P
import netCDF4 as nc
import sys
import os
import argparse
import time

#############
##VARIABLES##
#############

class VariableContainer:pass

#Constants
csts = VariableContainer()

csts.hartree2ev = N.float(27.211396132)
csts.ev2hartree = N.float(1/csts.hartree2ev)
csts.sqrtpi = N.float(N.sqrt(N.pi))
csts.invsqrtpi = N.float(1/csts.sqrtpi)
csts.TOLKPTS = N.float(0.00001)

###########
##CLASSES##
###########

class PolynomialFit(object):
    def __init__(self):
        self.degree = 2

class EigenvalueContainer(object):
    nsppol = None
    nkpt = None
    mband = None
    eigenvalues = None
    units = None
    wtk = None
    filename = None
    filefullpath = None
    bd_indices = None
    eigenvalue_type = None
    kpoints = None
    #kpoint_sampling_type: can be Monkhorst-Pack or Bandstructure
    KPT_W90_TOL = N.float(1.0e-6)
    KPT_DFT_TOL = N.float(1.0e-8)
    kpoint_sampling_type = 'Monkhorst-Pack'
    inputgvectors = None
    gvectors = None
    special_kpoints = None
    special_kpoints_names = None
    special_kpoints_indices = None
    kpoint_path_values = None
    kpoint_reduced_path_values = None
    kpoint_path_length = None
    #reduced_norm = None
    norm_paths = None
    norm_reduced_paths = None
    def __init__(self,directory=None,filename=None):
        if filename == None:return
        if directory == None:directory='.'
        self.filename = filename
        self.filefullpath = '%s/%s' %(directory,filename)
        self.file_open(self.filefullpath)
    def set_kpoint_sampling_type(self,kpoint_sampling_type):
        if kpoint_sampling_type != 'Monkhorst-Pack' and kpoint_sampling_type != 'Bandstructure':
            print 'ERROR: kpoint_sampling_type "%s" does not exists' %kpoint_sampling_type
            print '       it should be "Monkhorst-Pack" or "Bandstructure" ... exit'
            sys.exit()
        self.kpoint_sampling_type = kpoint_sampling_type
    def correct_kpt(self,kpoint,tolerance=N.float(1.0e-6)):
        kpt_correct = N.array(kpoint,N.float)
        changed = False
        for ii in range(3):
            if N.allclose(kpoint[ii],N.float(1.0/3.0),atol=tolerance):
                kpt_correct[ii] = N.float(1.0/3.0)
                changed = True
            elif N.allclose(kpoint[ii],N.float(1.0/6.0),atol=tolerance):
                kpt_correct[ii] = N.float(1.0/6.0)
                changed = True
            elif N.allclose(kpoint[ii],N.float(-1.0/6.0),atol=tolerance):
                kpt_correct[ii] = N.float(-1.0/6.0)
                changed = True
            elif N.allclose(kpoint[ii],N.float(-1.0/3.0),atol=tolerance):
                kpt_correct[ii] = N.float(-1.0/3.0)
                changed = True
        if changed:
            print 'COMMENT: kpoint %15.12f %15.12f %15.12f has been changed to %15.12f %15.12f %15.12f' %(kpoint[0],kpoint[1],kpoint[2],kpt_correct[0],kpt_correct[1],kpt_correct[2])
        return kpt_correct
    def find_special_kpoints(self,gvectors=None):
        if self.kpoint_sampling_type != 'Bandstructure':
            print 'ERROR: special kpoints are usefull only for bandstructures ... returning find_special_kpoints'
            return
        if self.eigenvalue_type == 'W90':
            correct_kpt_tolerance = N.float(1.0e-4)
            KPT_TOL = self.KPT_W90_TOL
        elif self.eigenvalue_type == 'DFT':
            correct_kpt_tolerance = N.float(1.0e-6)
            KPT_TOL = self.KPT_DFT_TOL
        else:
            print 'ERROR: eigenvalue_type is "%s" while it should be "W90" or "DFT" ... returning find_special_kpoints' %self.eigenvalue_type
            return
        if gvectors == None:
            self.inputgvectors = False
            self.gvectors = N.identity(3,N.float)
        else:
            if N.shape(gvectors) != (3, 3):
                print 'ERROR: wrong gvectors ... exiting now'
                sys.exit()
            self.inputgvectors = True
            self.gvectors = gvectors
        full_kpoints = N.zeros((self.nkpt,3),N.float)
        for ikpt in range(self.nkpt):
            full_kpoints[ikpt,:] = self.kpoints[ikpt,0]*self.gvectors[0,:]+self.kpoints[ikpt,1]*self.gvectors[1,:]+self.kpoints[ikpt,2]*self.gvectors[2,:]
        delta_kpt = full_kpoints[1,:]-full_kpoints[0,:]
        self.special_kpoints_indices = list()
        self.special_kpoints = list()
        self.special_kpoints_indices.append(0)
        self.special_kpoints.append(self.correct_kpt(self.kpoints[0,:],tolerance=correct_kpt_tolerance))
        for ikpt in range(1,self.nkpt-1):
            thisdelta = full_kpoints[ikpt+1,:]-full_kpoints[ikpt,:]
            if not N.allclose(thisdelta,delta_kpt,atol=KPT_TOL):
                delta_kpt = thisdelta
                self.special_kpoints_indices.append(ikpt)
                self.special_kpoints.append(self.correct_kpt(self.kpoints[ikpt,:],tolerance=correct_kpt_tolerance))
        self.special_kpoints_indices.append(N.shape(self.kpoints)[0]-1)
        self.special_kpoints.append(self.correct_kpt(self.kpoints[-1,:],tolerance=correct_kpt_tolerance))
        print 'Special Kpoints : '
        print ' {0:d} : {1[0]: 8.8f} {1[1]: 8.8f} {1[2]: 8.8f}'.format(1,self.kpoints[0,:])
        self.norm_paths = N.zeros((N.shape(self.special_kpoints_indices)[0]-1),N.float)
        self.norm_reduced_paths = N.zeros((N.shape(self.special_kpoints_indices)[0]-1),N.float)
        for ispkpt in range(1,N.shape(self.special_kpoints_indices)[0]):
            self.norm_paths[ispkpt-1] = N.linalg.norm(full_kpoints[self.special_kpoints_indices[ispkpt]]-full_kpoints[self.special_kpoints_indices[ispkpt-1]])
            self.norm_reduced_paths[ispkpt-1] = N.linalg.norm(self.special_kpoints[ispkpt]-self.special_kpoints[ispkpt-1])
            print '     {2:d}-{3:d} path length : {0: 8.8f} | reduced path length : {1: 8.8f}'.\
                  format(self.norm_paths[ispkpt-1],self.norm_reduced_paths[ispkpt-1],ispkpt,ispkpt+1)
            print ' {0:d} : {1[0]: 8.8f} {1[1]: 8.8f} {1[2]: 8.8f}'.format(ispkpt+1,self.kpoints[self.special_kpoints_indices[ispkpt],:])
        self.kpoint_path_length = N.sum(self.norm_paths)
        self.kpoint_reduced_path_length = N.sum(self.norm_reduced_paths)
        self.normalized_kpoint_path_norm = self.norm_paths/self.kpoint_path_length
        self.normalized_kpoint_reduced_path_norm = self.norm_reduced_paths/self.kpoint_reduced_path_length
        
        kptredpathval = list()
        kptpathval = list()
        kptredpathval.append(N.float(0.0))
        kptpathval.append(N.float(0.0))
        curlen = N.float(0.0)
        redcurlen = N.float(0.0)
        for ispkpt in range(1,N.shape(self.special_kpoints_indices)[0]):
            kptredpathval.extend(N.linspace(redcurlen,redcurlen+self.norm_reduced_paths[ispkpt-1],self.special_kpoints_indices[ispkpt]-self.special_kpoints_indices[ispkpt-1]+1)[1:])
            kptpathval.extend(N.linspace(curlen,curlen+self.norm_paths[ispkpt-1],self.special_kpoints_indices[ispkpt]-self.special_kpoints_indices[ispkpt-1]+1)[1:])
            redcurlen = redcurlen + self.norm_reduced_paths[ispkpt-1]
            curlen = curlen + self.norm_paths[ispkpt-1]
        self.kpoint_path_values = N.array(kptpathval,N.float)
        self.kpoint_reduced_path_values = N.array(kptredpathval,N.float)
        self.normalized_kpoint_path_values = self.kpoint_path_values/self.kpoint_path_length
        self.normalized_kpoint_reduced_path_values = self.kpoint_reduced_path_values/self.kpoint_reduced_path_length
        self.special_kpoints = N.array(self.special_kpoints,N.float)
    def file_open(self,filefullpath):
        if filefullpath[-3:] == '_GW':
            self.gw_file_open(filefullpath)
        elif filefullpath[-7:] == '_EIG.nc':
            self.nc_eig_open(filefullpath)
        elif filefullpath[-4:] == '.dat':
            self.wannier_bs_file_open(filefullpath)
    def has_eigenvalue(self,nsppol,isppol,kpoint,iband):
        if self.nsppol != nsppol:
            return False
        for ikpt in range(self.nkpt):
            if N.absolute(self.kpoints[ikpt,0]-kpoint[0]) < csts.TOLKPTS and \
               N.absolute(self.kpoints[ikpt,1]-kpoint[1]) < csts.TOLKPTS and \
               N.absolute(self.kpoints[ikpt,2]-kpoint[2]) < csts.TOLKPTS:
                if iband >= self.bd_indices[isppol,ikpt,0]-1 and iband < self.bd_indices[isppol,ikpt,1]:
                    return True
                return False
        return False
    def get_eigenvalue(self,nsppol,isppol,kpoint,iband):
        for ikpt in range(self.nkpt):
            if N.absolute(self.kpoints[ikpt,0]-kpoint[0]) < csts.TOLKPTS and \
               N.absolute(self.kpoints[ikpt,1]-kpoint[1]) < csts.TOLKPTS and \
               N.absolute(self.kpoints[ikpt,2]-kpoint[2]) < csts.TOLKPTS:
                return self.eigenvalues[isppol,ikpt,iband]
    def wannier_bs_file_open(self,filefullpath):
        if not (os.path.isfile(filefullpath)):
            print 'ERROR : file "%s" does not exists' %filefullpath
            print '... exiting now ...'
            sys.exit()
        print 'WARNING: no spin polarization reading yet for Wannier90 bandstructure files!'
        self.eigenvalue_type = 'W90'
        self.nsppol = None
        self.nkpt = None
        self.mband = None
        self.eigenvalues = None
        self.units = None
        self.filefullpath = filefullpath
        reader = open(self.filefullpath,'r')
        filedata = reader.readlines()
        reader.close()
        for iline in range(len(filedata)):
            if filedata[iline].strip() == '':
                self.nkpt = iline
                break
        self.mband = N.int(len(filedata)/self.nkpt)
        self.nsppol = 1
        self.eigenvalues = N.zeros([self.nsppol,self.nkpt,self.mband],N.float)
        self.kpoints = N.zeros([self.nkpt,3],N.float)
        iline = 0
        kpt_file = '%s.kpt' %filefullpath[:-4]
        if os.path.isfile(kpt_file):
            reader = open(kpt_file,'r')
            kptdata = reader.readlines()
            reader.close()
            if N.int(kptdata[0]) != self.nkpt:
                print 'ERROR : the number of kpoints in file "%s" is not the same as in "%s" ... exit' %(self.filefullpath,kpt_file)
                sys.exit()
            for ikpt in range(self.nkpt):
                linesplit = kptdata[ikpt+1].split()
                self.kpoints[ikpt,0] = N.float(linesplit[0])
                self.kpoints[ikpt,1] = N.float(linesplit[1])
                self.kpoints[ikpt,2] = N.float(linesplit[2])
        else:
            for ikpt in range(self.nkpt):
                self.kpoints[ikpt,0] = N.float(filedata[ikpt].split()[0])
        for iband in range(self.mband):
            for ikpt in range(self.nkpt):
                self.eigenvalues[0,ikpt,iband] = N.float(filedata[iline].split()[1])
                iline = iline+1
            iline = iline+1
        self.eigenvalues = self.eigenvalues*csts.ev2hartree
        self.units = 'Hartree'
    def gw_file_open(self,filefullpath):
        if not (os.path.isfile(filefullpath)):
            print 'ERROR : file "%s" does not exists' %filefullpath
            print '... exiting now ...'
            sys.exit()
        self.eigenvalue_type = 'GW'
        self.nsppol = None
        self.nkpt = None
        self.mband = None
        self.eigenvalues = None
        self.units = None
        self.filefullpath = filefullpath
        reader = open(self.filefullpath,'r')
        filedata = reader.readlines()
        reader.close()
        self.nkpt = N.int(filedata[0].split()[0])
        self.kpoints = N.ones([self.nkpt,3],N.float)
        self.nsppol = N.int(filedata[0].split()[1])
        self.bd_indices = N.zeros((self.nsppol,self.nkpt,2),N.int)
        icur = 1
        nbd_kpt = N.zeros([self.nsppol,self.nkpt],N.int)
        for isppol in range(self.nsppol):
            for ikpt in range(self.nkpt):
                self.kpoints[ikpt,:] = N.array(filedata[icur].split()[:],N.float)
                icur = icur + 1
                nbd_kpt[isppol,ikpt] = N.int(filedata[icur])
                self.bd_indices[isppol,ikpt,0] = N.int(filedata[icur+1].split()[0])
                self.bd_indices[isppol,ikpt,1] = N.int(filedata[icur+nbd_kpt[isppol,ikpt]].split()[0])
                icur = icur + nbd_kpt[isppol,ikpt] + 1
        self.mband = N.max(self.bd_indices[:,:,1])
        self.eigenvalues = N.zeros([self.nsppol,self.nkpt,self.mband],N.float)
        self.eigenvalues[:,:,:] = N.nan
        ii = 3
        for isppol in range(self.nsppol):
            for ikpt in range(self.nkpt):
                for iband in range(self.bd_indices[isppol,ikpt,0]-1,self.bd_indices[isppol,ikpt,1]):
                    self.eigenvalues[isppol,ikpt,iband] = N.float(filedata[ii].split()[1])
                    ii = ii + 1
                ii = ii + 2
        self.eigenvalues = csts.ev2hartree*self.eigenvalues
        self.units = 'Hartree'
    def pfit_gw_file_write(self,polyfitlist,directory=None,filename=None,bdgw=None,energy_pivots=None,gwec=None):
        if filename == None:return
        if directory == None:directory='.'
        filefullpath = '%s/%s' %(directory,filename)
        if (os.path.isfile(filefullpath)):
            user_input = raw_input('WARNING : file "%s" exists, do you want to overwrite it ? (y/n)' %filefullpath)
            if not (user_input == 'y' or user_input == 'Y'):
                return
        writer = open(filefullpath,'w')
        writer.write('%12s%12s\n' %(self.nkpt,self.nsppol))
        if gwec == None:
            for ikpt in range(self.nkpt):
                for isppol in range(self.nsppol):
                    writer.write('%10.6f%10.6f%10.6f\n' %(self.kpoints[ikpt,0],self.kpoints[ikpt,1],self.kpoints[ikpt,2]))
                    writer.write('%4i\n' %(bdgw[1]-bdgw[0]+1))
                    for iband in range(bdgw[0]-1,bdgw[1]):
                        delta = N.polyval(polyfitlist[-1],csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband])
                        for ipivot in range(len(energy_pivots)):
                            if csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband] <= energy_pivots[ipivot]:
                                delta = N.polyval(polyfitlist[ipivot],csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband])
                                break
                        writer.write('%6i%9.4f%9.4f%9.4f\n' %(iband+1,csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband]+delta,delta,0.0))
        else:
            for ikpt in range(self.nkpt):
                for isppol in range(self.nsppol):
                    writer.write('%10.6f%10.6f%10.6f\n' %(self.kpoints[ikpt,0],self.kpoints[ikpt,1],self.kpoints[ikpt,2]))
                    writer.write('%4i\n' %(bdgw[1]-bdgw[0]+1))
                    for iband in range(bdgw[0]-1,bdgw[1]):
                        if gwec.has_eigenvalue(self.nsppol,isppol,self.kpoints[ikpt],iband):
                            gw_eig = gwec.get_eigenvalue(self.nsppol,isppol,self.kpoints[ikpt],iband)
                            writer.write('%6i%9.4f%9.4f%9.4f\n' %(iband+1,csts.hartree2ev*gw_eig,csts.hartree2ev*(gw_eig-self.eigenvalues[isppol,ikpt,iband]),0.0))
                        else:
                            delta = N.polyval(polyfitlist[-1],csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband])
                            for ipivot in range(len(energy_pivots)):
                                if csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband] <= energy_pivots[ipivot]:
                                    delta = N.polyval(polyfitlist[ipivot],csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband])
                                    break
                            writer.write('%6i%9.4f%9.4f%9.4f\n' %(iband+1,csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband]+delta,delta,0.0))
        writer.close()
    def pfit_dft_to_gw_bs_write(self,polyfitlist,directory=None,filename=None,bdgw=None,energy_pivots=None,gwec=None):
        if filename == None:return
        if directory == None:directory='.'
        filefullpath = '%s/%s' %(directory,filename)
        if (os.path.isfile(filefullpath)):
            user_input = raw_input('WARNING : file "%s" exists, do you want to overwrite it ? (y/n)' %filefullpath)
            if not (user_input == 'y' or user_input == 'Y'):
                return
        writer = open(filefullpath,'w')
        if gwec == None:
            for ikpt in range(self.nkpt):
                writer.write('%s' %ikpt)
                for isppol in range(self.nsppol):
                    for iband in range(bdgw[0]-1,bdgw[1]):
                        delta = N.polyval(polyfitlist[-1],csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband])
                        for ipivot in range(len(energy_pivots)):
                            if csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband] <= energy_pivots[ipivot]:
                                delta = N.polyval(polyfitlist[ipivot],csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband])
                                break
                        writer.write(' %s' %(csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband]+delta))
                writer.write('\n')
        else:
            print 'NOT SUPPORTED YET'
            sys.exit()
        writer.close()
    def nc_eig_open(self,filefullpath):
        if not (os.path.isfile(filefullpath)):
            print 'ERROR : file "%s" does not exists' %filefullpath
            print '... exiting now ...'
            sys.exit()
        ncdata = nc.Dataset(filefullpath)
        self.eigenvalue_type = 'DFT'
        self.nsppol = None
        self.nkpt = None
        self.mband = None
        self.eigenvalues = None
        self.units = None
        self.filefullpath = filefullpath
        for dimname,dimobj in ncdata.dimensions.iteritems():
            if dimname == 'nsppol':self.nsppol = N.int(len(dimobj))
            if dimname == 'nkpt':self.nkpt = N.int(len(dimobj))
            if dimname == 'mband':self.mband = N.int(len(dimobj))
        for varname in ncdata.variables:
            if varname == 'Eigenvalues':
                varobj = ncdata.variables[varname]
                varshape = N.shape(varobj[:])
                self.units = None
                for attrname in varobj.ncattrs():
                    if attrname == 'units':
                        self.units = varobj.getncattr(attrname)
                if self.units == None:
                    print 'WARNING : units are not specified'
                    print '... assuming "Hartree" units ...'
                    self.units = 'Hartree'
                elif self.units != 'Hartree':
                    print 'ERROR : units are unknown : "%s"' %self.units
                    print '... exiting now ...'
                    sys.exit()
                self.eigenvalues = N.reshape(N.array(varobj,N.float),varshape)
                self.nsppol = varshape[0]
                self.nkpt = varshape[1]
                self.kpoints = -1*N.ones((self.nkpt,3),N.float)
                self.mband = varshape[2]
                self.bd_indices = N.zeros((self.nsppol,self.nkpt,2),N.int)
                self.bd_indices[:,:,0] = 1
                self.bd_indices[:,:,1] = self.mband
                break
        for varname in ncdata.variables:
            if varname == 'Kptns':
                varobj = ncdata.variables[varname]
                varshape = N.shape(varobj[:])
                self.kpoints = N.reshape(N.array(varobj,N.float),varshape)
    def write_bandstructure_to_file(self,filename,option_kpts='bohrm1_units'):
        #if option_kpts is set to 'normalized', the path of the bandstructure will be normalized to 1 (and special k-points correctly chosen)
        if self.kpoint_sampling_type != 'Bandstructure':
            print 'ERROR: kpoint_sampling_type is not "Bandstructure" ... returning from write_bandstructure_to_file'
            return
        if self.nsppol > 1:
            print 'ERROR: number of spins is more than 1, this is not yet coded ... returning from write_bandstructure_to_file'
            return
        writer = open(filename,'w')
        writer.write('# BANDSTRUCTURE FILE FROM DAVID\'S SCRIPT\n')
        writer.write('# nsppol = %s\n' %self.nsppol)
        writer.write('# nband = %s\n' %self.mband)
        writer.write('# eigenvalue_type = %s\n' %self.eigenvalue_type)
        if self.inputgvectors:
            writer.write('# inputgvectors = 1 (%s)\n' %self.inputgvectors)
        else:
            writer.write('# inputgvectors = 0 (%s)\n' %self.inputgvectors)
        writer.write('# gvectors(1) = %20.17f %20.17f %20.17f \n' %(self.gvectors[0,0],self.gvectors[0,1],self.gvectors[0,2]))
        writer.write('# gvectors(2) = %20.17f %20.17f %20.17f \n' %(self.gvectors[1,0],self.gvectors[1,1],self.gvectors[1,2]))
        writer.write('# gvectors(3) = %20.17f %20.17f %20.17f \n' %(self.gvectors[2,0],self.gvectors[2,1],self.gvectors[2,2]))
        writer.write('# special_kpoints_number = %s\n' %(len(self.special_kpoints_indices)))
        writer.write('# list of special kpoints : (given in reduced coordinates, value_path is in Bohr^-1, value_red_path has its total path normalized to 1)\n')
        for ii in range(len(self.special_kpoints_indices)):
            ispkpt = self.special_kpoints_indices[ii]
            spkpt = self.special_kpoints[ii]
            writer.write('#    special_kpt_index %5s : %20.17f %20.17f %20.17f (value_path = %20.17f | value_red_path = %20.17f)\n' %(ispkpt,spkpt[0],spkpt[1],spkpt[2],self.kpoint_path_values[ispkpt],self.kpoint_reduced_path_values[ispkpt]))
        writer.write('# special_kpoints_names :\n')
        for ii in range(len(self.special_kpoints_indices)):
            ispkpt = self.special_kpoints_indices[ii]
            spkpt = self.special_kpoints[ii]
            writer.write('#    special_kpt_name %3s : "%s" : %20.17f %20.17f %20.17f\n' %(ii+1,self.special_kpoints_names[ii],spkpt[0],spkpt[1],spkpt[2]))
        writer.write('# kpoint_path_length = %20.17f \n' %(self.kpoint_path_length))
        writer.write('# kpoint_path_number = %s \n' %(self.nkpt))
        if self.inputgvectors:
            writer.write('# kpoint_path_units = %s\n' %(option_kpts))
        else:
            writer.write('# kpoint_path_units =  %s (!!! CONSIDERING UNITARY GVECTORS MATRIX !!!)\n' %(option_kpts))
        writer.write('#BEGIN\n')
        if option_kpts == 'bohrm1_units':
            values_path = self.kpoint_path_values
        elif option_kpts == 'reduced':
            values_path = self.kpoint_reduced_path_values
        elif option_kpts == 'bohrm1_units_normalized':
            values_path = self.normalized_kpoint_path_values
        elif option_kpts == 'reduced_normalized':
            values_path = self.normalized_kpoint_reduced_path_values
        else:
            print 'ERROR: wrong option_kpts ... exit'
            writer.write('... CANCELLED (wrong option_kpts)')
            writer.close()
            sys.exit()
        for isppol in range(self.nsppol):
            writer.write('#isppol %s\n' %isppol)
            for iband in range(self.mband):
                writer.write('#iband %5s (band number %s)\n' %(iband,iband+1))
                for ikpt in range(self.nkpt):
                    writer.write('%20.17f %20.17f\n' %(values_path[ikpt],self.eigenvalues[isppol,ikpt,iband]))
                writer.write('\n')
        writer.write('#END\n')
        writer.write('\n#KPT_LIST\n')
        for ikpt in range(self.nkpt):
            writer.write('# %6d : %20.17f %20.17f %20.17f\n' %(ikpt,self.kpoints[ikpt,0],self.kpoints[ikpt,1],self.kpoints[ikpt,2]))
        writer.close()
    def read_bandstructure_from_file(self,filename):
        reader = open(filename,'r')
        bs_data = reader.readlines()
        reader.close()
        self.gvectors = N.identity(3,N.float)
        self.kpoint_sampling_type = 'Bandstructure'
        self.special_kpoints_indices = list()
        self.special_kpoints = list()
        for ii in range(len(bs_data)):
            if bs_data[ii] == '#BEGIN\n':
                ibegin = ii
                break
            elif bs_data[ii][:10] == '# nsppol =':
                self.nsppol = N.int(bs_data[ii][10:])
            elif bs_data[ii][:9] == '# nband =':
                self.mband = N.int(bs_data[ii][9:])
            elif bs_data[ii][:19] == '# eigenvalue_type =':
                self.eigenvalue_type = bs_data[ii][19:].strip()
            elif bs_data[ii][:17] == '# inputgvectors =':
                tt = N.int(bs_data[ii][18])
                if tt == 1:
                    self.inputgvectors = True
                elif tt == 0:
                    self.inputgvectors = False
                else:
                    print 'ERROR: reading inputgvectors ... exit'
                    sys.exit()
            elif bs_data[ii][:15] == '# gvectors(1) =':
                sp = bs_data[ii][15:].split()
                self.gvectors[0,0] = N.float(sp[0])
                self.gvectors[0,1] = N.float(sp[1])
                self.gvectors[0,2] = N.float(sp[2])
            elif bs_data[ii][:15] == '# gvectors(2) =':
                sp = bs_data[ii][15:].split()
                self.gvectors[1,0] = N.float(sp[0])
                self.gvectors[1,1] = N.float(sp[1])
                self.gvectors[1,2] = N.float(sp[2])
            elif bs_data[ii][:15] == '# gvectors(3) =':
                sp = bs_data[ii][15:].split()
                self.gvectors[2,0] = N.float(sp[0])
                self.gvectors[2,1] = N.float(sp[1])
                self.gvectors[2,2] = N.float(sp[2])
            elif bs_data[ii][:26] == '# special_kpoints_number =':
                special_kpoints_number = N.int(bs_data[ii][26:])
                self.special_kpoints_names = ['']*special_kpoints_number
            elif bs_data[ii][:22] == '#    special_kpt_index':
                sp = bs_data[ii][22:].split()
                self.special_kpoints_indices.append(N.int(sp[0]))
                self.special_kpoints.append(N.array([sp[2],sp[3],sp[4]]))
            elif bs_data[ii][:21] == '#    special_kpt_name':
                sp = bs_data[ii][21:].split()
                ispkpt = N.int(sp[0])-1
                self.special_kpoints_names[ispkpt] = sp[2][1:-1]
            elif bs_data[ii][:22] == '# kpoint_path_length =':
                self.kpoint_path_length = N.float(bs_data[ii][22:])
            elif bs_data[ii][:22] == '# kpoint_path_number =':
                self.nkpt = N.int(bs_data[ii][22:])
            elif bs_data[ii][:21] == '# kpoint_path_units =':
                kpoint_path_units = bs_data[ii][21:].strip()
        self.special_kpoints_indices = N.array(self.special_kpoints_indices,N.int)
        self.special_kpoints = N.array(self.special_kpoints,N.float)
        if len(self.special_kpoints_indices) != special_kpoints_number or len(self.special_kpoints) != special_kpoints_number:
            print 'ERROR: reading the special kpoints ... exit'
            sys.exit()
        self.kpoint_path_values = N.zeros([self.nkpt],N.float)
        self.kpoint_reduced_path_values = N.zeros([self.nkpt],N.float)
        if kpoint_path_units == 'bohrm1_units':
            jj = 0
            for ii in range(ibegin+1,len(bs_data)):
                if bs_data[ii][:7] == '#isppol' or bs_data[ii][:6] == '#iband':continue
                if bs_data[ii] == '\n':
                    break
                self.kpoint_path_values[jj] = N.float(bs_data[ii].split()[0])
                jj = jj + 1
            if jj != self.nkpt:
                print 'ERROR: reading bandstructure file ... exit'
                sys.exit()
            self.normalized_kpoint_path_values = self.kpoint_path_values/self.kpoint_path_length
        if kpoint_path_units == 'bohrm1_units_normalized':
            jj = 0
            for ii in range(ibegin+1,len(bs_data)):
                if bs_data[ii][:7] == '#isppol' or bs_data[ii][:6] == '#iband':continue
                if bs_data[ii] == '\n':
                    break
                self.normalized_kpoint_path_values[jj] = N.float(bs_data[ii].split()[0])
                jj = jj + 1
            if jj != self.nkpt:
                print 'ERROR: reading bandstructure file ... exit'
                sys.exit()
            self.kpoint_path_values = self.normalized_kpoint_path_values*self.kpoint_path_length
        elif kpoint_path_units == 'reduced_normalized':
            jj = 0
            for ii in range(ibegin+1,len(bs_data)):
                if bs_data[ii][:7] == '#isppol' or bs_data[ii][:6] == '#iband':continue
                if bs_data[ii] == '\n':
                    break
                self.normalized_kpoint_reduced_path_values[jj] = N.float(bs_data[ii].split()[0])
                jj = jj + 1
            if jj != self.nkpt:
                print 'ERROR: reading bandstructure file ... exit'
                sys.exit()
            self.kpoint_reduced_path_values = self.normalized_kpoint_reduced_path_values/self.kpoint_reduced_path_length
        elif kpoint_path_units == 'reduced':
            jj = 0
            for ii in range(ibegin+1,len(bs_data)):
                if bs_data[ii][:7] == '#isppol' or bs_data[ii][:6] == '#iband':continue
                if bs_data[ii] == '\n':
                    break
                self.kpoint_reduced_path_values[jj] = N.float(bs_data[ii].split()[0])
                jj = jj + 1
            if jj != self.nkpt:
                print 'ERROR: reading bandstructure file ... exit'
                sys.exit()
            self.normalized_kpoint_reduced_path_values = self.kpoint_reduced_path_values/self.kpoint_reduced_path_length
        self.eigenvalues = N.zeros([self.nsppol,self.nkpt,self.mband],N.float)
        check_nband = 0
        for ii in range(ibegin+1,len(bs_data)):
            if bs_data[ii][:7] == '#isppol':
                isppol = N.int(bs_data[ii][7:])
            elif bs_data[ii][:6] == '#iband':
                iband = N.int(bs_data[ii][6:].split()[0])
                ikpt = 0
            elif bs_data[ii][:4] == '#END':
                break
            elif bs_data[ii] == '\n':
                check_nband = check_nband + 1
            else:
                self.eigenvalues[isppol,ikpt,iband] = N.float(bs_data[ii].split()[1])
                ikpt = ikpt + 1

def check_gw_vs_dft_parameters(dftec,gwec):
    if gwec.eigenvalue_type != 'GW' or dftec.eigenvalue_type != 'DFT':
        print 'ERROR: eigenvalue files do not contain GW and DFT eigenvalues ... exiting now'
        sys.exit()
    if dftec.nsppol != gwec.nsppol or dftec.nkpt != gwec.nkpt:
        print 'ERROR: the number of spins/kpoints is not the same in the GW and DFT files used to make the interpolation ... exiting now'
        sys.exit()
    for ikpt in range(dftec.nkpt):
        if N.absolute(dftec.kpoints[ikpt,0]-gwec.kpoints[ikpt,0]) > csts.TOLKPTS or \
           N.absolute(dftec.kpoints[ikpt,1]-gwec.kpoints[ikpt,1]) > csts.TOLKPTS or \
           N.absolute(dftec.kpoints[ikpt,2]-gwec.kpoints[ikpt,2]) > csts.TOLKPTS:
            print 'ERROR: the kpoints are not the same in the GW and DFT files used to make the interpolation ... exiting now'
            sys.exit()

def plot_gw_vs_dft_eig(dftec,gwec,vbm_index,energy_pivots=None,polyfit_degrees=None):
    if gwec.eigenvalue_type != 'GW' or dftec.eigenvalue_type != 'DFT':
        print 'ERROR: eigenvalue containers do not contain GW and DFT eigenvalues ... exiting now'
        sys.exit()
    if dftec.nsppol != gwec.nsppol or dftec.nkpt != gwec.nkpt:
        print 'ERROR: the number of spins/kpoints is not the same in the GW and DFT containers ... exiting now'
        sys.exit()
    valdftarray = N.array([],N.float)
    conddftarray = N.array([],N.float)
    valgwarray = N.array([],N.float)
    condgwarray = N.array([],N.float)
    for isppol in range(dftec.nsppol):
        for ikpt in range(dftec.nkpt):
            ibdmin = N.max([dftec.bd_indices[isppol,ikpt,0],gwec.bd_indices[isppol,ikpt,0]])-1
            ibdmax = N.min([dftec.bd_indices[isppol,ikpt,1],gwec.bd_indices[isppol,ikpt,1]])-1
            valdftarray = N.append(valdftarray,csts.hartree2ev*dftec.eigenvalues[isppol,ikpt,ibdmin:vbm_index])
            valgwarray = N.append(valgwarray,csts.hartree2ev*gwec.eigenvalues[isppol,ikpt,ibdmin:vbm_index])
            conddftarray = N.append(conddftarray,csts.hartree2ev*dftec.eigenvalues[isppol,ikpt,vbm_index:ibdmax+1])
            condgwarray = N.append(condgwarray,csts.hartree2ev*gwec.eigenvalues[isppol,ikpt,vbm_index:ibdmax+1])
    if energy_pivots == None:
        if plot_figures == 1:
            P.figure(1)
            P.hold(True)
            P.grid(True)
            P.plot(valdftarray,valgwarray,'bx')
            P.plot(conddftarray,condgwarray,'rx')
            P.xlabel('DFT eigenvalues (in eV)')
            P.ylabel('GW eigenvalues (in eV)')
            P.figure(2)
            P.hold(True)
            P.grid(True)
            P.plot(valdftarray,valgwarray-valdftarray,'bx')
            P.plot(conddftarray,condgwarray-conddftarray,'rx')
            P.xlabel('DFT eigenvalues (in eV)')
            P.ylabel('GW correction to the DFT eigenvalues (in eV)')
            P.show()
            return
    polyfitlist = list()
    if len(polyfit_degrees) == 1:
        print 'ERROR: making a fit with only one interval is not allowed ... exiting now'
        sys.exit()
    dftarray = N.append(valdftarray,conddftarray)
    gwarray = N.append(valgwarray,condgwarray)
    dftarray_list = list()
    gwarray_list = list()
    for iinterval in range(len(polyfit_degrees)):
        tmpdftarray = N.array([],N.float)
        tmpgwarray = N.array([],N.float)
        if iinterval == 0:
            emin = None
            emax = energy_pivots[0]
            for ii in range(len(dftarray)):
                if dftarray[ii] <= emax:
                    tmpdftarray = N.append(tmpdftarray,[dftarray[ii]])
                    tmpgwarray = N.append(tmpgwarray,[gwarray[ii]])
        elif iinterval == len(polyfit_degrees)-1:
            emin = energy_pivots[-1]
            emax = None
            for ii in range(len(dftarray)):
                if dftarray[ii] >= emin:
                    tmpdftarray = N.append(tmpdftarray,[dftarray[ii]])
                    tmpgwarray = N.append(tmpgwarray,[gwarray[ii]])
        else:
            emin = energy_pivots[iinterval-1]
            emax = energy_pivots[iinterval]
            for ii in range(len(dftarray)):
                if dftarray[ii] >= emin and dftarray[ii] <= emax:
                    tmpdftarray = N.append(tmpdftarray,[dftarray[ii]])
                    tmpgwarray = N.append(tmpgwarray,[gwarray[ii]])
        dftarray_list.append(tmpdftarray)
        gwarray_list.append(tmpgwarray)
        pfit = N.polyfit(tmpdftarray,tmpgwarray-tmpdftarray,polyfit_degrees[iinterval])
        polyfitlist.append(pfit)
    if plot_figures == 1:
        linspace_npoints = 200
        valpoly_x = N.linspace(N.min(valdftarray),N.max(valdftarray),linspace_npoints)
        condpoly_x = N.linspace(N.min(conddftarray),N.max(conddftarray),linspace_npoints)
        P.figure(3)
        P.hold(True)
        P.grid(True)
        P.plot(valdftarray,valgwarray-valdftarray,'bx')
        P.plot(conddftarray,condgwarray-conddftarray,'rx')
        [x_min,x_max] = P.xlim()
        for iinterval in range(len(polyfit_degrees)):
            if iinterval == 0:
                tmppoly_x = N.linspace(x_min,energy_pivots[iinterval],linspace_npoints)
            elif iinterval == len(polyfit_degrees)-1:
                tmppoly_x = N.linspace(energy_pivots[iinterval-1],x_max,linspace_npoints)
            else:
                tmppoly_x = N.linspace(energy_pivots[iinterval-1],energy_pivots[iinterval],linspace_npoints)
            P.plot(tmppoly_x,N.polyval(polyfitlist[iinterval],tmppoly_x),'k')
        for ipivot in range(len(energy_pivots)):
            en = energy_pivots[ipivot]
            P.plot([en,en],[N.polyval(polyfitlist[ipivot],en),N.polyval(polyfitlist[ipivot+1],en)],'k-.')
        P.xlabel('DFT eigenvalues (in eV)')
        P.ylabel('GW correction to the DFT eigenvalues (in eV)')
        P.figure(4)
        P.hold(True)
        P.grid(True)
        for iinterval in range(len(polyfit_degrees)):
            P.plot(dftarray_list[iinterval],gwarray_list[iinterval]-dftarray_list[iinterval]-N.polyval(polyfitlist[iinterval],dftarray_list[iinterval]),'bx')
        [x_min,x_max] = P.xlim()
        P.plot([x_min,x_max],[0,0],'k-')
        P.xlabel('DFT eigenvalues (in eV)')
        P.ylabel('Error in the fit (in eV)')
        P.show()
    return polyfitlist

def compare_bandstructures(ec_ref,ec_test):
    nspkpt_ref = len(ec_ref.special_kpoints)
    nspkpt_test = len(ec_test.special_kpoints)
    if nspkpt_ref != nspkpt_test:
        print 'ERROR: The number of special kpoints is different in the two files ... exit'
        sys.exit()
    eig_type_ref = ec_ref.eigenvalue_type
    eig_type_test = ec_test.eigenvalue_type
    print eig_type_ref,eig_type_test
    if eig_type_ref == 'DFT' and eig_type_test == 'W90':
        TOL_KPTS = N.float(1.0e-4)
    else:
        TOL_KPTS = N.float(1.0e-6)
    print TOL_KPTS
    for ispkpt in range(nspkpt_ref):
        print 'difference between the two :',ec_ref.special_kpoints[ispkpt,:]-ec_test.special_kpoints[ispkpt,:]
        if not N.allclose(ec_ref.special_kpoints[ispkpt,:],ec_test.special_kpoints[ispkpt,:],atol=TOL_KPTS):
            print 'ERROR: The kpoints are not the same :'
            print '       Kpt #%s ' %ispkpt
            print '         Reference => %20.17f %20.17f %20.17f' %(ec_ref.special_kpoints[ispkpt,0],ec_ref.special_kpoints[ispkpt,1],ec_ref.special_kpoints[ispkpt,2])
            print '         Compared  => %20.17f %20.17f %20.17f' %(ec_test.special_kpoints[ispkpt,0],ec_test.special_kpoints[ispkpt,1],ec_test.special_kpoints[ispkpt,2])
            print '       ... exit'
            sys.exit()
    mband_comparison = N.min([ec_ref.mband,ec_test.mband])
    if mband_comparison < ec_ref.mband:
        print 'Number of bands in the test bandstructure is lower than the number of bands in the reference (%s)' %ec_ref.mband
        print '  => Comparison will proceed with %s bands' %ec_test.mband
    elif mband_comparison < ec_test.mband:
        print 'Number of bands in the reference bandstructure is lower than the number of bands in the test bandstructure (%s)' %ec_test.mband
        print '  => Comparison will only proceed with %s bands of the test bandstructure' %ec_ref.mband
    else:
        print 'Number of bands in the reference and test bandstructure is the same'
        print '  => Comparison will proceed with %s bands' %mband_comparison
#    eig_test_ref_path = ec_ref.eigenvalues[:,:,:mband_comparison]
    rmsd_per_band = N.zeros([ec_ref.nsppol,mband_comparison],N.float)
    nrmsd_per_band = N.zeros([ec_ref.nsppol,mband_comparison],N.float)
    mae_per_band = N.zeros([ec_ref.nsppol,mband_comparison],N.float)
    for isppol in range(ec_ref.nsppol):
        for iband in range(mband_comparison):
            interp = N.interp(ec_ref.normalized_kpoint_path_values,ec_test.normalized_kpoint_path_values,ec_test.eigenvalues[isppol,:,iband])
            rmsd_per_band[isppol,iband] = N.sqrt(N.sum((csts.hartree2ev*interp-csts.hartree2ev*ec_ref.eigenvalues[isppol,:,iband])**2)/ec_ref.nkpt)
            mae_per_band[isppol,iband] = N.sum(N.abs(csts.hartree2ev*interp-csts.hartree2ev*ec_ref.eigenvalues[isppol,:,iband]))/ec_ref.nkpt
    P.figure(1)
    P.plot(mae_per_band[0,:])
    P.figure(2)
    P.plot(rmsd_per_band[0,:])
    P.show()

def get_gvectors():
    if os.path.isfile('.gvectors.bsinfo'):
        print 'File ".gvectors.bsinfo found with the following gvectors information :"'
        try:
            gvectors_reader = open('.gvectors.bsinfo','r')
            gvectors_data = gvectors_reader.readlines()
            gvectors_reader.close()
            trial_gvectors = N.identity(3,N.float)
            trial_gvectors[0,0] = N.float(gvectors_data[0].split()[0])
            trial_gvectors[0,1] = N.float(gvectors_data[0].split()[1])
            trial_gvectors[0,2] = N.float(gvectors_data[0].split()[2])
            trial_gvectors[1,0] = N.float(gvectors_data[1].split()[0])
            trial_gvectors[1,1] = N.float(gvectors_data[1].split()[1])
            trial_gvectors[1,2] = N.float(gvectors_data[1].split()[2])
            trial_gvectors[2,0] = N.float(gvectors_data[2].split()[0])
            trial_gvectors[2,1] = N.float(gvectors_data[2].split()[1])
            trial_gvectors[2,2] = N.float(gvectors_data[2].split()[2])
            print ' gvectors(1) = [ %20.17f %20.17f %20.17f ]' %(trial_gvectors[0,0],trial_gvectors[0,1],trial_gvectors[0,2])
            print ' gvectors(2) = [ %20.17f %20.17f %20.17f ]' %(trial_gvectors[1,0],trial_gvectors[1,1],trial_gvectors[1,2])
            print ' gvectors(3) = [ %20.17f %20.17f %20.17f ]' %(trial_gvectors[2,0],trial_gvectors[2,1],trial_gvectors[2,2])
        except:
            print 'ERROR: file ".gvectors.bsinfo" might be corrupted (empty or not formatted correctly ...)'
            print '       you should remove the file and start again or check the file ... exit'
            sys.exit()
        test = raw_input('Press <ENTER> to use these gvectors (any other character to enter manually other gvectors)\n')
        if test == '':
            gvectors = trial_gvectors
        else:
            gvectors = N.identity(3,N.float)
            test = raw_input('Enter G1 (example : "0.153 0 0") : \n')
            gvectors[0,0] = N.float(test.split()[0])
            gvectors[0,1] = N.float(test.split()[1])
            gvectors[0,2] = N.float(test.split()[2])
            test = raw_input('Enter G2 (example : "0.042 1.023 0") : \n')
            gvectors[1,0] = N.float(test.split()[0])
            gvectors[1,1] = N.float(test.split()[1])
            gvectors[1,2] = N.float(test.split()[2])
            test = raw_input('Enter G3 (example : "0 0 1.432") : \n')
            gvectors[2,0] = N.float(test.split()[0])
            gvectors[2,1] = N.float(test.split()[1])
            gvectors[2,2] = N.float(test.split()[2])
            test = raw_input('Do you want to overwrite the gvectors contained in the file ".gvectors.bsinfo" ? (<ENTER> for yes, anything else for no)\n')
            if test == '':
                print 'Writing gvectors to file ".gvectors.bsinfo" ...'
                gvectors_writer = open('.gvectors.bsinfo','w')
                gvectors_writer.write('%20.17f %20.17f %20.17f\n' %(trial_gvectors[0,0],trial_gvectors[0,1],trial_gvectors[0,2]))
                gvectors_writer.write('%20.17f %20.17f %20.17f\n' %(trial_gvectors[1,0],trial_gvectors[1,1],trial_gvectors[1,2]))
                gvectors_writer.write('%20.17f %20.17f %20.17f\n' %(trial_gvectors[2,0],trial_gvectors[2,1],trial_gvectors[2,2]))
                gvectors_writer.close()
                print '... done'
    else: 
        test = raw_input('Do you want to enter the the reciprocal space primitive vectors (y/n)\n')
        if test == 'y':
            gvectors = N.identity(3,N.float)
            test = raw_input('Enter G1 (example : "0.153 0 0") : ')
            gvectors[0,0] = N.float(test.split()[0])
            gvectors[0,1] = N.float(test.split()[1])
            gvectors[0,2] = N.float(test.split()[2])
            test = raw_input('Enter G2 (example : "0.042 1.023 0") : ')
            gvectors[1,0] = N.float(test.split()[0])
            gvectors[1,1] = N.float(test.split()[1])
            gvectors[1,2] = N.float(test.split()[2])
            test = raw_input('Enter G3 (example : "0 0 1.432") : ')
            gvectors[2,0] = N.float(test.split()[0])
            gvectors[2,1] = N.float(test.split()[1])
            gvectors[2,2] = N.float(test.split()[2])
            test = raw_input('Do you want to write the gvectors to file ".gvectors.bsinfo" ? (<ENTER> for yes, anything else for no)\n')
            if test == '':
                print 'Writing gvectors to file ".gvectors.bsinfo" ...'
                gvectors_writer = open('.gvectors.bsinfo','w')
                gvectors_writer.write('%20.17f %20.17f %20.17f\n' %(gvectors[0,0],gvectors[0,1],gvectors[0,2]))
                gvectors_writer.write('%20.17f %20.17f %20.17f\n' %(gvectors[1,0],gvectors[1,1],gvectors[1,2]))
                gvectors_writer.write('%20.17f %20.17f %20.17f\n' %(gvectors[2,0],gvectors[2,1],gvectors[2,2]))
                gvectors_writer.close()
                print '... done'
        else:
            gvectors = None
    return gvectors

# Parse the command line options

parser = argparse.ArgumentParser(description='Tool for analyzing wannier90 bandstructures')
parser.add_argument('files',help='files to be opened',nargs='*')
args = parser.parse_args()
args_dict = vars(args)

if args_dict['files']:
    print 'will open the file'
else:
    print 'ERROR: you should provide some bandstructure file ! exiting now ...'
    sys.exit()
w90_file_list = args_dict['files']
ec_w90_list = list()
if len(w90_file_list) == 1:
    gvectors = get_gvectors()
    w90_file = w90_file_list[0]
    ec_w90 = EigenvalueContainer(directory='.',filename=w90_file)
    ec_w90.set_kpoint_sampling_type('Bandstructure')
    ec_w90.find_special_kpoints(gvectors)
    ec_w90_list.append(ec_w90)
elif len(w90_file_list) > 1:
    print 'There is more than 1 file. What do you want to do ?'
    print ' 1 => merge all the files with different states (k-point path should be the same !) (default)'
    test = raw_input(' 0 => quit\n... ')
    if test == '0':
        print '... exiting now'
        sys.exit()
    elif test == '1' or test == '':
        option_multiple = 1
        print 'Script will try to merge the files with the different states'
        gvectors = get_gvectors()
    
        for iw90 in range(len(w90_file_list)):
            w90_file = w90_file_list[iw90]
            ec_w90 = EigenvalueContainer(directory='.',filename=w90_file)
            ec_w90.set_kpoint_sampling_type('Bandstructure')
            ec_w90.find_special_kpoints(gvectors)
            ec_w90_list.append(ec_w90)
    if option_multiple == 1:
        nband_tot = 0
        for iw90 in range(len(ec_w90_list)):
            nband_tot = nband_tot + ec_w90_list[iw90].mband
        ec_w90 = EigenvalueContainer()
        ec_w90.mband = nband_tot
        ec_w90.nkpt = ec_w90_list[0].nkpt
        ec_w90.units = 'Hartree'
        ec_w90.kpoints = ec_w90_list[0].kpoints
        ec_w90.eigenvalue_type = 'W90'
        ec_w90.nsppol = 1
        isppol = 0
        ec_w90.eigenvalues = N.zeros((ec_w90.nsppol,ec_w90.nkpt,ec_w90.mband))
        #ec_w90.reduced_norm=ec_w90_list[0].reduced_norm
        ec_w90.kpoint_path_values=ec_w90_list[0].kpoint_path_values
        ec_w90.kpoint_reduced_path_values=ec_w90_list[0].kpoint_reduced_path_values
        ec_w90.normalized_kpoint_path_values = ec_w90_list[0].normalized_kpoint_path_values
        ec_w90.normalized_kpoint_reduced_path_values = ec_w90_list[0].normalized_kpoint_reduced_path_values
        ec_w90.inputgvectors = ec_w90_list[0].inputgvectors
        ec_w90.gvectors = ec_w90_list[0].gvectors
        ec_w90.special_kpoints = ec_w90_list[0].special_kpoints
        ec_w90.special_kpoints_indices = ec_w90_list[0].special_kpoints_indices
        ec_w90.kpoint_path_length = ec_w90_list[0].kpoint_path_length
        ec_w90.norm_paths = ec_w90_list[0].norm_paths
        ec_w90.norm_reduced_paths = ec_w90_list[0].norm_reduced_paths
        ec_w90.kpoint_sampling_type = ec_w90_list[0].kpoint_sampling_type

        iband = 0
        for iw90 in range(len(ec_w90_list)):
            for ibd_thisw90 in range(ec_w90_list[iw90].mband):
                for ikpt in range(ec_w90.nkpt):
                    ec_w90.eigenvalues[isppol,ikpt,iband] = ec_w90_list[iw90].eigenvalues[isppol,ikpt,ibd_thisw90]
                iband = iband+1

print 'Number of bands in the file : %s' %(N.shape(ec_w90.eigenvalues)[2])
test = raw_input('Enter the number of bands to be plotted (<ENTER> : %s) : \n' %(N.shape(ec_w90.eigenvalues)[2]))
if test == '':
    nbd_plot = N.shape(ec_w90.eigenvalues)[2]
else:
    nbd_plot = N.int(test)
if nbd_plot > N.shape(ec_w90.eigenvalues)[2]:
    print 'ERROR: the number of bands to be plotted is larger than the number available ... exit'
    sys.exit()

ec_w90.special_kpoints_names = ['']*len(ec_w90.special_kpoints_indices)
for ii in range(len(ec_w90.special_kpoints_indices)):
    ec_w90.special_kpoints_names[ii] = 'k%s' %(ii+1)
print 'List of special kpoints :'
for ii in range(len(ec_w90.special_kpoints_indices)):
    spkpt = ec_w90.kpoints[ec_w90.special_kpoints_indices[ii]]
    print ' Kpoint %s : %s %s %s' %(ii+1,spkpt[0],spkpt[1],spkpt[2])
print 'Enter the name of the %s special k-points :' %(len(ec_w90.special_kpoints_indices))
test = raw_input('')
if len(test.split()) == len(ec_w90.special_kpoints_indices):
    for ii in range(len(ec_w90.special_kpoints_indices)):
        ec_w90.special_kpoints_names[ii] = test.split()[ii]

#ec_dft = EigenvalueContainer(filename='quartz_bso_EIG.nc')
#ec_dft.set_kpoint_sampling_type('Bandstructure')
#ec_dft.find_special_kpoints(gvectors)

#compare_bandstructures(ec_dft,ec_w90)

test = raw_input('Enter base name for bandstructure file : \n')
ec_w90.write_bandstructure_to_file('%s.bandstructure' %test)


P.figure(1,figsize=(3.464,5))
P.hold('on')
P.grid('on')
P.xticks(N.take(ec_w90.kpoint_reduced_path_values,N.array(ec_w90.special_kpoints_indices,N.int)),ec_w90.special_kpoints_names)
for iband in range(nbd_plot):
    P.plot(ec_w90.kpoint_reduced_path_values,ec_w90.eigenvalues[0,:,iband]*csts.hartree2ev,'k-',linewidth=2)
#ec_w90.read_bandstructure_from_file('mytest.bandstructure')
#ec_w90.write_bandstructure_to_file('mytest2.bandstructure')
#ec_w90_2 = EigenvalueContainer()
#ec_w90_2.read_bandstructure_from_file('mytest.bandstructure')
#ec_w90_2.write_bandstructure_to_file('mytest2.bandstructure')
#for iband in range(nbd_plot):
#    P.plot(ec_dft.kpoint_reduced_path_values,ec_dft.eigenvalues[0,:,iband]*csts.hartree2ev,'r--',linewidth=2)
P.show()
