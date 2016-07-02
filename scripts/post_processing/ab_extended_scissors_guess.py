#!/usr/bin/python
#====================================================================#
#  Script to get the eigenvalues from an abinit _EIG.nc netcdf file  #
#====================================================================#

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
csts.fig_width = 8
csts.fig_height = 6
csts.markersize = 10
csts.markeredgewidth = 2


###########
##METHODS##
###########

#Get the coefficients of the line going through 2 points of xi,yi and xj,yj. By default, begin and end points
#If the indices are given (2 integers), the line must go through these two points
def line_ab(x,y,indices=None):
    if len(N.shape(x)) != 1 or len(N.shape(y)) != 1:
        print 'ERROR: the array x and/or y is not 1D ... exit'
        sys.exit()
    if len(x) != len(y):
        print 'ERROR: x and y arrays have different lengths ... exit'
        sys.exit()
    if indices == None:
        indices = N.array([0,len(x)-1])
    else:
        if indices[0] < 0 or indices[1] >= len(x) or indices[0] == indices[1]:
            print 'ERROR: indices (0 <= indices[0]=%s < indices[1]=%s < len(x)=%s) are wrong ... exit' %(indices[0],indices[1],len(x))
            sys.exit()
    a = (y[indices[0]]-y[indices[1]])/(x[indices[0]]-x[indices[1]])
    b = (y[indices[0]]*x[indices[1]]-x[indices[0]]*y[indices[1]])/(x[indices[1]]-x[indices[0]])
    return N.array([a,b],N.float)

#Get the coefficients of the polynom of degree 2 going through 2 points xi,yi and xj,yj and using a least squares procedure for the other points
#If the indices are given (2 integers), the polynom must go through these two points
def poly2d_ab(x,y,indices=None):
    if len(N.shape(x)) != 1 or len(N.shape(y)) != 1:
        print 'ERROR: the array x and/or y is not 1D ... exit'
        sys.exit()
    if len(x) != len(y):
        print 'ERROR: x and y arrays have different lengths ... exit'
        sys.exit()
    if indices == None:
        indices = N.array([0,len(x)-1])
    else:
        if indices[0] < 0 or indices[1] >= len(x) or indices[0] == indices[1]:
            print 'ERROR: indices (0 <= indices[0]=%s < indices[1]=%s < len(x)=%s) are wrong ... exit' %(indices[0],indices[1],len(x))
            sys.exit()
    x1 = x[indices[0]]
    x2 = x[indices[1]]
    y1 = y[indices[0]]
    y2 = y[indices[1]]
    x3 = (x1+x2)/2
    newy = y - y1*(x-x2)*(x-x3)/((x1-x2)*(x1-x3)) - y2*(x-x1)*(x-x3)/((x2-x1)*(x2-x3))
    A = N.vstack([4*(x*x-(x1+x2)*x+x1*x2)/(2*x1*x2-x2*x2-x1*x1)]).T
    y3 = N.linalg.lstsq(A,newy)[0]
    pp = N.polyfit(N.array([x1,x2,x3]),N.array([y1,y2,y3]),2)
    return pp

#Get the coefficients of the polynom of degree "degree" going through 2 points xi,yi and xj,yj and using a least squares procedure for the other points
#If the indices are given (2 integers), the polynom must go through these two points
def polynd_ab(x,y,degree,indices=None):
    if degree < 1:
        print 'ERROR: cannot find a polynomial going through two points with this degree (%s) ... exit' %degree
        sys.exit()
    if len(N.shape(x)) != 1 or len(N.shape(y)) != 1:
        print 'ERROR: the array x and/or y is not 1D ... exit'
        sys.exit()
    if len(x) != len(y):
        print 'ERROR: x and y arrays have different lengths ... exit'
        sys.exit()
    if indices == None:
        indices = N.array([0,len(x)-1])
    else:
        if indices[0] < 0 or indices[1] >= len(x) or indices[0] == indices[1]:
            print 'ERROR: indices (0 <= indices[0]=%s < indices[1]=%s < len(x)=%s) are wrong ... exit' %(indices[0],indices[1],len(x))
            sys.exit()
    if degree == 1:
        pp = N.polyfit(N.array([x[indices[0]],x[indices[1]]]),N.array([y[indices[0]],y[indices[1]]]),degree)
        return pp
    x1 = x[indices[0]]
    x2 = x[indices[1]]
    y1 = y[indices[0]]
    y2 = y[indices[1]]
    xm = N.linspace(N.min(x),N.max(x),degree,endpoint=False)[1:]
    prod1 = (x-x2)/(x1-x2)*y1
    prod2 = (x-x1)/(x2-x1)*y2
    coeff_list = list()
    for ii in range(len(xm)):
        prod1 = prod1*(x-xm[ii])/(x1-xm[ii])
        prod2 = prod2*(x-xm[ii])/(x2-xm[ii])
        prod_ii = (x-x2)*(x-x1)/((xm[ii]-x2)*(xm[ii]-x1))
        for jj in range(len(xm)):
            if ii != jj:
                prod_ii = prod_ii*(x-xm[jj])/(xm[ii]-xm[jj])
        coeff_list.append(prod_ii)
    p1 = prod1 + prod2
    newy = y - p1
    A = N.vstack(N.array(coeff_list)).T
    ym = N.linalg.lstsq(A,newy)[0]
    xx = N.array([x1])
    yy = N.array([y1])
    for ii in range(len(xm)):
        xx = N.append(xx,[xm[ii]])
        yy = N.append(yy,[ym[ii]])
    xx = N.append(xx,[x2])
    yy = N.append(yy,[y2])
    pp = N.polyfit(xx,yy,degree)
    return pp

#Get the coefficients of the polynom of degree "degree" going through 1 point xi,yi and using a least squares procedure for the other points
#If the indice is given (1 integer), the polynom must go through this point, otherwise, it takes the first point in the list by default
def polynd_a(x,y,degree,indices=None):
    if degree < 1:
        print 'ERROR: cannot find a polynomial going through one points with this degree (%s) ... exit' %degree
        sys.exit()
    if len(N.shape(x)) != 1 or len(N.shape(y)) != 1:
        print 'ERROR: the array x and/or y is not 1D ... exit'
        sys.exit()
    if len(x) != len(y):
        print 'ERROR: x and y arrays have different lengths ... exit'
        sys.exit()
    if indices == None:
        indices = N.array([0])
    elif len(indices) != 1:
        print 'ERROR: there should be only one index of a point through which the polynom has to go through ... exit'
        sys.exit()
    else:
        if indices[0] < 0 or indices[0] >= len(x):
            print 'ERROR: index (0 <= indices[0]=%s < len(x)=%s) is wrong ... exit' %(indices[0],len(x))
            sys.exit()
    x1 = x[indices[0]]
    y1 = y[indices[0]]
    xm = N.linspace(N.min(x),N.max(x),degree+1,endpoint=True)[1:]
    prod1 = y1
    coeff_list = list()
    for ii in range(len(xm)):
        prod1 = prod1*(x-xm[ii])/(x1-xm[ii])
        prod_ii = (x-x1)/(xm[ii]-x1)
        for jj in range(len(xm)):
            if ii != jj:
                prod_ii = prod_ii*(x-xm[jj])/(xm[ii]-xm[jj])
        coeff_list.append(prod_ii)
    newy = y - prod1
    A = N.vstack(N.array(coeff_list)).T
    ym = N.linalg.lstsq(A,newy)[0]
    xx = N.array([x1])
    yy = N.array([y1])
    for ii in range(len(xm)):
        xx = N.append(xx,[xm[ii]])
        yy = N.append(yy,[ym[ii]])
    pp = N.polyfit(xx,yy,degree)
    return pp

#Given the polyfit_list and energy_pivots, finds the 2nd degree that goes to a zero slope (where the value will be separated by delta_energy)
#starting from a given energy (derivative is the same at this point). Returns the polynom and the values of the vertex of the polynom (maximum or minimum)
def smoothend(energy_pivots,polyfit_list,energy,delta_energy_ev=None):
    method = 2
    if delta_energy_ev == None:
        if method == 1:
            delta_energy_ev = 0.05
        elif method == 2:
            delta_energy_ev = 1.00
    if method == 1:
        xi = energy
        if xi < energy_pivots[0]:
            print 'Error: energy should be larger than the first energy pivot ...'
            sys.exit()
        if xi > energy_pivots[-1]:
            ii = len(energy_pivots)
        else:
            ii = N.argwhere(energy_pivots>xi)[0]
        fpi = N.polyval(N.polyder(polyfit_list[ii]),xi)
        fi = N.polyval(polyfit_list[ii],xi)
        if fpi == 0:
            print 'TODO, the derivative is 0 ... easy but lazy :-)'
            return
        elif fpi > 0:
            fv = fi + delta_energy_ev
        elif fpi < 0:
            fv = fi - delta_energy_ev
        aa = fpi**2/(4*(fi-fv))
        bb = fpi - 2*aa*xi
        cc = fi - aa*xi**2 - bb*xi
        xv = -bb/(2*aa)
        new_energy_pivots = N.zeros(ii+2,N.float)
        for jj in N.arange(ii):
            new_energy_pivots[jj] = energy_pivots[jj]
        new_energy_pivots[-2] = energy
        new_energy_pivots[-1] = xv
        new_polyfit_list = list()
        for jj in N.arange(ii+1):
            new_polyfit_list.append(polyfit_list[jj])
        new_polyfit_list.append([aa,bb,cc])
        new_polyfit_list.append([fv])
        return new_energy_pivots,new_polyfit_list
    if method == 2:
        xi = energy
        if xi < energy_pivots[0]:
            print 'Error: energy should be larger than the first energy pivot ...'
            sys.exit()
        if xi > energy_pivots[-1]:
            ii = len(energy_pivots)
        else:
            ii = N.argwhere(energy_pivots>xi)[0]
        fpi = N.polyval(N.polyder(polyfit_list[ii]),xi)
        fi = N.polyval(polyfit_list[ii],xi)
        if fpi == 0:
            new_energy_pivots = N.zeros(ii+1,N.float)
            for jj in N.arange(ii):
                new_energy_pivots[jj] = energy_pivots[jj]
            new_energy_pivots[-1] = energy
            new_polyfit_list = list()
            for jj in N.arange(ii+1):
                new_polyfit_list.append(polyfit_list[jj])
            new_polyfit_list.append([fi])
            return new_energy_pivots,new_polyfit_list
        else:
            xv = xi + delta_energy_ev
            bb = fpi/(1.0-xi/xv)
            aa = -bb/(2.0*xv)
            cc = fi - aa*xi**2 - bb*xi
            pp = [aa,bb,cc]
            fv = N.polyval(pp,xv)
            new_energy_pivots = N.zeros(ii+2,N.float)
            for jj in N.arange(ii):
                new_energy_pivots[jj] = energy_pivots[jj]
            new_energy_pivots[-2] = xi
            new_energy_pivots[-1] = xv
            new_polyfit_list = list()
            for jj in N.arange(ii+1):
                new_polyfit_list.append(polyfit_list[jj])
            new_polyfit_list.append(pp)
            new_polyfit_list.append([fv])
            if N.abs(fv-fi) > 0.05:
                print 'WARNING : the last energy pivot is more than 0.05 eV from the constant correction'
            else:
                print 'COMMENT : smoothing the end of the graph starting at energy {0: 8.8f} eV'.format(xi)
            print '  => constant correction for higher states : fv = {0: 8.8f} eV'.format(fv)
            print '  => last energy pivot :                     fi = {0: 8.8f} eV'.format(fi)
            print '  => fv - fi = {0: 8.8f} eV'.format(fv-fi)
            return new_energy_pivots,new_polyfit_list

def write_polyfit(filename,energypivots,polyfitlist,energypivots_2=None,polyfitlist_2=None):
    writer = open(filename,'w')
    if energypivots_2 == None:
        writer.write('nsppol 1\n')
    else:
        print 'write_polyfit not yet implemented for nsppol = 2 ... returning'
        writer.write('NOT IMPLEMENTED\n')
        writer.close()
        return
    writer.write('%s\n' %len(polyfitlist))
    for ie in range(len(energypivots)):
        writer.write('%s  ' %energypivots[ie])
    writer.write('\n')
    for ip in range(len(polyfitlist)):
        pfit = polyfitlist[ip]
        for ic in range(len(pfit)):
            writer.write('%s  ' %pfit[ic])
        writer.write('\n')
    writer.close()

def read_polyfit(filename):
    reader = open(filename,'r')
    data = reader.readlines()
    if data[0][:8] != 'nsppol 1':
        print data[0]
        print data[0][:8]
        print 'read_polyfit not yet implemented for nsppol != 1 ... returning'
        reader.close()
        return
    npfit = N.int(data[1])
    energypivots = N.zeros(npfit-1,N.float)
    polyfitlist = list()
    for ie, ep in enumerate(data[2].split()):
        energypivots[ie] = N.float(ep)
    for ip in range(npfit):
        sp = data[3+ip].split()
        tmp = N.zeros(len(sp),N.float)
        for ic, cc in enumerate(sp):
            tmp[ic] = N.float(cc)
        polyfitlist.append(tmp)
    return energypivots,polyfitlist

###########
##CLASSES##
###########

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
    GROUP_BANDS_BS_TOL_EV = N.float(0.01)
    GROUP_BANDS_BS_TOL = GROUP_BANDS_BS_TOL_EV*csts.ev2hartree
    GROUP_BANDS_TOL_EV = N.float(0.2)
    GROUP_BANDS_TOL = GROUP_BANDS_TOL_EV*csts.ev2hartree
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
    def file_open(self,filefullpath):
        if filefullpath[-3:] == '_GW':
            self.gw_file_open(filefullpath)
        else:
            self.nc_eig_open(filefullpath)
    def set_kpoint_sampling_type(self,kpoint_sampling_type):
        if kpoint_sampling_type != 'Monkhorst-Pack' and kpoint_sampling_type != 'Bandstructure':
            print 'ERROR: kpoint_sampling_type "%s" does not exists' %kpoint_sampling_type
            print '       it should be "Monkhorst-Pack" or "Bandstructure" ... exit'
            sys.exit()
        self.kpoint_sampling_type = kpoint_sampling_type
    def find_band_groups(self,bandstructure_file=None,tolerance_ev=None,spinchoice=None):
        if self.nsppol > 1:
            print 'WARNING: find_band_groups is carefully checked only for nsppol = 1'
        if spinchoice == None:
            print 'COMMENT: find_band_groups handles spins up and down on equal footing'
            spinchoice = 'common'
        elif spinchoice == 'common':
            print 'COMMENT: find_band_groups handles spins up and down on equal footing'
        elif spinchoice == 'separate':
            print 'COMMENT: find_band_groups handles spins up and down as 2 different band structures'
        if bandstructure_file != None:
            ec_bs = EigenvalueContainer(filename=bandstructure_file)
            eigenvalues = ec_bs.eigenvalues
            nkpt = ec_bs.nkpt
            nband = ec_bs.mband
            nsppol = ec_bs.nsppol
        else:
            eigenvalues = self.eigenvalues
            nkpt = self.nkpt
            nband = self.mband
            nsppol = self.nsppol
        if tolerance_ev == None:
            if bandstructure_file == None:
                tolerance = self.GROUP_BANDS_TOL
            else:
                tolerance = self.GROUP_BANDS_BS_TOL
        else:
            tolerance = tolerance_ev*csts.ev2hartree
        if spinchoice == 'common':
            energy_pivots_list = list()
            band = eigenvalues[:,:,0]
            for iband in range(1,nband):
                if N.min(eigenvalues[:,:,iband]) - N.max(band) > tolerance:
                    energy_pivots_list.append((N.min(eigenvalues[:,:,iband]) + N.max(band))/2)
                band = eigenvalues[:,:,iband]
            return N.array(energy_pivots_list)
        elif spinchoice == 'separate':
            energy_pivots_list_up = list()
            energy_pivots_list_down = list()
            bandup = eigenvalues[0,:,0]
            banddown = eigenvalues[1,:,0]
            for iband in range(1,nband):
                if N.min(eigenvalues[0,:,iband]) - N.max(bandup) > tolerance:
                    energy_pivots_list_up.append((N.min(eigenvalues[0,:,iband]) + N.max(bandup))/2)
                bandup = eigenvalues[0,:,iband]
                if N.min(eigenvalues[1,:,iband]) - N.max(banddown) > tolerance:
                    energy_pivots_list_down.append((N.min(eigenvalues[1,:,iband]) + N.max(banddown))/2)
                banddown = eigenvalues[1,:,iband]
            return N.array(energy_pivots_list_up),N.array(energy_pivots_list_down)
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
        for ikpt in range(self.nkpt):
            for isppol in range(self.nsppol):
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
        for ikpt in range(self.nkpt):
            for isppol in range(self.nsppol):
                for iband in range(self.bd_indices[isppol,ikpt,0]-1,self.bd_indices[isppol,ikpt,1]):
                    self.eigenvalues[isppol,ikpt,iband] = N.float(filedata[ii].split()[1])
                    ii = ii + 1
                ii = ii + 2
        self.eigenvalues = csts.ev2hartree*self.eigenvalues
        self.units = 'Hartree'
    def pfit_gw_eigenvalues_ha(self,polyfitlist_up,energy_pivots_up=None,nband=None,polyfitlist_down=None,energy_pivots_down=None,ecgw=None):
        if polyfitlist_down == None and energy_pivots_down != None:
            print 'ERROR: list of polyfits and energy pivots are not coherent ... exit'
            sys.exit()
        if polyfitlist_down != None and energy_pivots_down == None:
            print 'ERROR: list of polyfits and energy pivots are not coherent ... exit'
            sys.exit()
        if nband == None:
            mband = N.shape(self.eigenvalues)[2]
        else:
            mband = nband
        pfit_eigenvalues = csts.hartree2ev*N.array(self.eigenvalues)
        if polyfitlist_down == None:
            for ikpt in range(self.nkpt):
                for isppol in range(self.nsppol):
                    if ecgw == None:
                        ibdmin = 0
                        ibdmax = mband
                    else:
                        ibdmin = ecgw.bd_indices[isppol,ikpt,0]-1
                        ibdmax = ecgw.bd_indices[isppol,ikpt,1]
                    for iband in range(ibdmin,ibdmax):
                        delta = N.polyval(polyfitlist_up[-1],csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband])
                        for ipivot in range(len(energy_pivots_up)):
                            if csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband] <= energy_pivots_up[ipivot]:
                                delta = N.polyval(polyfitlist_up[ipivot],csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband])
                                break
                        pfit_eigenvalues[isppol,ikpt,iband] = self.eigenvalues[isppol,ikpt,iband]*csts.hartree2ev + delta
            return pfit_eigenvalues*csts.ev2hartree
        else:
            for ikpt in range(self.nkpt):
                isppol = 0
                if ecgw == None:
                    ibdmin = 0
                    ibdmax = mband
                else:
                    print ecgw.bd_indices
                    ibdmin = ecgw.bd_indices[isppol,0,0]-1
                    ibdmax = ecgw.bd_indices[isppol,0,1]
                for iband in range(ibdmin,ibdmax):
                    delta = N.polyval(polyfitlist_up[-1],csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband])
                    for ipivot in range(len(energy_pivots_up)):
                        if polyfitlist_up[ipivot] != None:
                            if csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband] <= energy_pivots_up[ipivot]:
                                delta = N.polyval(polyfitlist_up[ipivot],csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband])
                                break
                    pfit_eigenvalues[isppol,ikpt,iband] = self.eigenvalues[isppol,ikpt,iband]*csts.hartree2ev + delta
                isppol = 1
                if ecgw == None:
                    ibdmin = 0
                    ibdmax = mband
                else:
                    ibdmin = ecgw.bd_indices[isppol,0,0]-1
                    ibdmax = ecgw.bd_indices[isppol,0,1]
                for iband in range(ibdmin,ibdmax):
                    delta = N.polyval(polyfitlist_down[-1],csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband])
                    for ipivot in range(len(energy_pivots_down)):
                        if polyfitlist_down[ipivot] != None:
                            if csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband] <= energy_pivots_down[ipivot]:
                                delta = N.polyval(polyfitlist_down[ipivot],csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband])
                                break
                    pfit_eigenvalues[isppol,ikpt,iband] = self.eigenvalues[isppol,ikpt,iband]*csts.hartree2ev + delta
            return pfit_eigenvalues*csts.ev2hartree
    def pfit_gw_eigenvalues(self,polyfitlist,energy_pivots=None,nband=None):
        if nband == None:
            mband = N.shape(self.eigenvalues)[2]
        else:
            mband = nband
        pfit_eigenvalues = N.zeros((self.nsppol,self.nkpt,mband))
        for ikpt in range(self.nkpt):
            for isppol in range(self.nsppol):
                for iband in range(mband):
                    delta = N.polyval(polyfitlist[-1],csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband])
                    for ipivot in range(len(energy_pivots)):
                        if csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband] <= energy_pivots[ipivot]:
                            delta = N.polyval(polyfitlist[ipivot],csts.hartree2ev*self.eigenvalues[isppol,ikpt,iband])
                            break
                    pfit_eigenvalues[isppol,ikpt,iband] = self.eigenvalues[isppol,ikpt,iband]*csts.hartree2ev + delta
        return pfit_eigenvalues
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
            print 'ERROR: number of spins is more than 1, this is not fully tested ... use with care !'
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
#    def write_bandstructure_to_file(self,filename,option_kpts='bohrm1_units'):
#        #if option_kpts is set to 'normalized', the path of the bandstructure will be normalized to 1 (and special k-points correctly chosen)
#        if self.kpoint_sampling_type != 'Bandstructure':
#            print 'ERROR: kpoint_sampling_type is not "Bandstructure" ... returning from write_bandstructure_to_file'
#            return
#        if self.nsppol > 1:
#            print 'ERROR: number of spins is more than 1, this is not yet coded ... returning from write_bandstructure_to_file'
#            return
#        writer = open(filename,'w')
#        writer.write('# BANDSTRUCTURE FILE FROM DAVID\'S SCRIPT\n')
#        writer.write('# nsppol = %s\n' %self.nsppol)
#        writer.write('# nband = %s\n' %self.mband)
#        writer.write('# eigenvalue_type = %s\n' %self.eigenvalue_type)
#        if self.inputgvectors:
#            writer.write('# inputgvectors = 1 (%s)\n' %self.inputgvectors)
#        else:
#            writer.write('# inputgvectors = 0 (%s)\n' %self.inputgvectors)
#        writer.write('# gvectors(1) = %20.17f %20.17f %20.17f \n' %(self.gvectors[0,0],self.gvectors[0,1],self.gvectors[0,2]))
#        writer.write('# gvectors(2) = %20.17f %20.17f %20.17f \n' %(self.gvectors[1,0],self.gvectors[1,1],self.gvectors[1,2]))
#        writer.write('# gvectors(3) = %20.17f %20.17f %20.17f \n' %(self.gvectors[2,0],self.gvectors[2,1],self.gvectors[2,2]))
#        writer.write('# special_kpoints_number = %s\n' %(len(self.special_kpoints_indices)))
#        writer.write('# list of special kpoints : (given in reduced coordinates, value_path is in Bohr^-1, value_red_path has its total path normalized to 1)\n')
#        for ii in range(len(self.special_kpoints_indices)):
#            ispkpt = self.special_kpoints_indices[ii]
#            spkpt = self.special_kpoints[ii]
#            writer.write('#    special_kpt_index %5s : %20.17f %20.17f %20.17f (value_path = %20.17f | value_red_path = %20.17f)\n' %(ispkpt,spkpt[0],spkpt[1],spkpt[2],self.kpoint_path_values[ispkpt],self.kpoint_reduced_path_values[ispkpt]))
#        writer.write('# special_kpoints_names :\n')
#        for ii in range(len(self.special_kpoints_indices)):
#            ispkpt = self.special_kpoints_indices[ii]
#            spkpt = self.special_kpoints[ii]
#            writer.write('#    special_kpt_name %3s : "%s" : %20.17f %20.17f %20.17f\n' %(ii+1,self.special_kpoints_names[ii],spkpt[0],spkpt[1],spkpt[2]))
#        writer.write('# kpoint_path_length = %20.17f \n' %(self.kpoint_path_length))
#        writer.write('# kpoint_path_number = %s \n' %(self.nkpt))
#        if self.inputgvectors:
#            writer.write('# kpoint_path_units = %s\n' %(option_kpts))
#        else:
#            writer.write('# kpoint_path_units =  %s (!!! CONSIDERING UNITARY GVECTORS MATRIX !!!)\n' %(option_kpts))
#        writer.write('#BEGIN\n')
#        if option_kpts == 'bohrm1_units':
#            values_path = self.kpoint_path_values
#        elif option_kpts == 'reduced':
#            values_path = self.kpoint_reduced_path_values
#        elif option_kpts == 'bohrm1_units_normalized':
#            values_path = self.normalized_kpoint_path_values
#        elif option_kpts == 'reduced_normalized':
#            values_path = self.normalized_kpoint_reduced_path_values
#        else:
#            print 'ERROR: wrong option_kpts ... exit'
#            writer.write('... CANCELLED (wrong option_kpts)')
#            writer.close()
#            sys.exit()
#        for isppol in range(self.nsppol):
#            writer.write('#isppol %s\n' %isppol)
#            for iband in range(self.mband):
#                writer.write('#iband %5s (band number %s)\n' %(iband,iband+1))
#                for ikpt in range(self.nkpt):
#                    writer.write('%20.17f %20.17f\n' %(values_path[ikpt],self.eigenvalues[isppol,ikpt,iband]))
#                writer.write('\n')
#        writer.write('#END\n')
#        writer.write('\n#KPT_LIST\n')
#        for ikpt in range(self.nkpt):
#            writer.write('# %6d : %20.17f %20.17f %20.17f\n' %(ikpt,self.kpoints[ikpt,0],self.kpoints[ikpt,1],self.kpoints[ikpt,2]))
#        writer.close()
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

def classify_eigenvalues(eigarray,energy_pivots,eig2array=None):
    eigarray_list = list()
    if eig2array != None:
        eig2array_list = list()
        for iinterval in range(len(energy_pivots)+1):
            print iinterval,' : '
            tmpeigarray = N.array([],N.float)
            tmpeig2array = N.array([],N.float)
            if iinterval == 0:
                emin = None
                emax = energy_pivots[0]
                print emin,emax
                for ii in range(len(eigarray)):
                    if eigarray[ii] <= emax:
                        tmpeigarray = N.append(tmpeigarray,[eigarray[ii]])
                        tmpeig2array = N.append(tmpeig2array,[eig2array[ii]])
            elif iinterval == len(energy_pivots):
                emin = energy_pivots[-1]
                emax = None
                print emin,emax
                for ii in range(len(eigarray)):
                    if eigarray[ii] >= emin:
                        tmpeigarray = N.append(tmpeigarray,[eigarray[ii]])
                        tmpeig2array = N.append(tmpeig2array,[eig2array[ii]])
            else:
                emin = energy_pivots[iinterval-1]
                emax = energy_pivots[iinterval]
                print emin,emax
                for ii in range(len(eigarray)):
                    if eigarray[ii] >= emin and eigarray[ii] <= emax:
                        tmpeigarray = N.append(tmpeigarray,[eigarray[ii]])
                        tmpeig2array = N.append(tmpeig2array,[eig2array[ii]])
            eigarray_list.append(tmpeigarray)
            eig2array_list.append(tmpeig2array)
        return eigarray_list,eig2array_list
    else:
        for iinterval in range(len(energy_pivots)+1):
            tmpeigarray = N.array([],N.float)
            if iinterval == 0:
                emin = None
                emax = energy_pivots[0]
                for ii in range(len(eigarray)):
                    if eigarray[ii] <= emax:
                        tmpeigarray = N.append(tmpeigarray,[eigarray[ii]])
            elif iinterval == len(energy_pivots):
                emin = energy_pivots[-1]
                emax = None
                for ii in range(len(eigarray)):
                    if eigarray[ii] >= emin:
                        tmpeigarray = N.append(tmpeigarray,[eigarray[ii]])
            else:
                emin = energy_pivots[iinterval-1]
                emax = energy_pivots[iinterval]
                for ii in range(len(eigarray)):
                    if eigarray[ii] >= emin and eigarray[ii] <= emax:
                        tmpeigarray = N.append(tmpeigarray,[eigarray[ii]])
            eigarray_list.append(tmpeigarray)
        return eigarray_list

def plot_gw_vs_dft_eig(dftec,gwec,vbm_index,energy_pivots_up=None,energy_pivots_down=None,polyfit_degrees_up=None,polyfit_degrees_down=None,limitpoints=None,spinchoice=None,smooth_end=True,smooth_energy=None,smooth_delta_energy=None):
    DELTA_ENERGY_END = 2.0
    if gwec.eigenvalue_type != 'GW' or dftec.eigenvalue_type != 'DFT':
        print 'ERROR: eigenvalue containers do not contain GW and DFT eigenvalues ... exiting now'
        sys.exit()
    if dftec.nsppol != gwec.nsppol or dftec.nkpt != gwec.nkpt:
        print 'ERROR: the number of spins/kpoints is not the same in the GW and DFT containers ... exiting now'
        sys.exit()
    if dftec.nsppol == 1:
        spinchoice = 'common'
    valdftarray = N.array([],N.float)
    conddftarray = N.array([],N.float)
    valgwarray = N.array([],N.float)
    condgwarray = N.array([],N.float)
    if dftec.nsppol == 2 and spinchoice == 'separate':
        upvaldftarray = N.array([],N.float)
        upconddftarray = N.array([],N.float)
        upvalgwarray = N.array([],N.float)
        upcondgwarray = N.array([],N.float)
        downvaldftarray = N.array([],N.float)
        downconddftarray = N.array([],N.float)
        downvalgwarray = N.array([],N.float)
        downcondgwarray = N.array([],N.float)
    if spinchoice == None or spinchoice == 'common':
        for ikpt in range(dftec.nkpt):
            for isppol in range(dftec.nsppol):
                ibdmin = N.max([dftec.bd_indices[isppol,ikpt,0],gwec.bd_indices[isppol,ikpt,0]])-1
                ibdmax = N.min([dftec.bd_indices[isppol,ikpt,1],gwec.bd_indices[isppol,ikpt,1]])-1
                valdftarray = N.append(valdftarray,csts.hartree2ev*dftec.eigenvalues[isppol,ikpt,ibdmin:vbm_index])
                valgwarray = N.append(valgwarray,csts.hartree2ev*gwec.eigenvalues[isppol,ikpt,ibdmin:vbm_index])
                conddftarray = N.append(conddftarray,csts.hartree2ev*dftec.eigenvalues[isppol,ikpt,vbm_index:ibdmax+1])
                condgwarray = N.append(condgwarray,csts.hartree2ev*gwec.eigenvalues[isppol,ikpt,vbm_index:ibdmax+1])
    elif spinchoice == 'separate':
        for ikpt in range(dftec.nkpt):
            isppol = 0
            ibdmin = N.max([dftec.bd_indices[isppol,ikpt,0],gwec.bd_indices[isppol,ikpt,0]])-1
            ibdmax = N.min([dftec.bd_indices[isppol,ikpt,1],gwec.bd_indices[isppol,ikpt,1]])-1
            upvaldftarray = N.append(upvaldftarray,csts.hartree2ev*dftec.eigenvalues[isppol,ikpt,ibdmin:vbm_index])
            upvalgwarray = N.append(upvalgwarray,csts.hartree2ev*gwec.eigenvalues[isppol,ikpt,ibdmin:vbm_index])
            upconddftarray = N.append(upconddftarray,csts.hartree2ev*dftec.eigenvalues[isppol,ikpt,vbm_index:ibdmax+1])
            upcondgwarray = N.append(upcondgwarray,csts.hartree2ev*gwec.eigenvalues[isppol,ikpt,vbm_index:ibdmax+1])
            isppol = 1
            ibdmin = N.max([dftec.bd_indices[isppol,ikpt,0],gwec.bd_indices[isppol,ikpt,0]])-1
            ibdmax = N.min([dftec.bd_indices[isppol,ikpt,1],gwec.bd_indices[isppol,ikpt,1]])-1
            downvaldftarray = N.append(downvaldftarray,csts.hartree2ev*dftec.eigenvalues[isppol,ikpt,ibdmin:vbm_index])
            downvalgwarray = N.append(downvalgwarray,csts.hartree2ev*gwec.eigenvalues[isppol,ikpt,ibdmin:vbm_index])
            downconddftarray = N.append(downconddftarray,csts.hartree2ev*dftec.eigenvalues[isppol,ikpt,vbm_index:ibdmax+1])
            downcondgwarray = N.append(downcondgwarray,csts.hartree2ev*gwec.eigenvalues[isppol,ikpt,vbm_index:ibdmax+1])
    if energy_pivots_up == None:
        if plot_figures == 1:
            if dftec.nsppol == 2:
                P.figure(1,figsize=(csts.fig_width,csts.fig_height))
                P.hold(True)
                P.grid(True)
                P.plot(upvaldftarray,upvalgwarray,'bx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
                P.plot(upconddftarray,upcondgwarray,'rx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
                P.xlabel('DFT eigenvalues - spin UP (in eV)')
                P.ylabel('GW eigenvalues (in eV)')
                P.figure(2,figsize=(csts.fig_width,csts.fig_height))
                P.hold(True)
                P.grid(True)
                P.plot(downvaldftarray,downvalgwarray,'bo',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
                P.plot(downconddftarray,downcondgwarray,'ro',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
                P.xlabel('DFT eigenvalues - spin UP (in eV)')
                P.ylabel('GW eigenvalues (in eV)')
            else:
                P.figure(1,figsize=(csts.fig_width,csts.fig_height))
                P.hold(True)
                P.grid(True)
                P.plot(valdftarray,valgwarray,'bx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
                P.plot(conddftarray,condgwarray,'rx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
                P.xlabel('DFT eigenvalues (in eV)')
                P.ylabel('GW eigenvalues (in eV)')
            if dftec.nsppol == 2:
                P.figure(3,figsize=(csts.fig_width,csts.fig_height))
                P.hold(True)
                P.grid(True)
                P.plot(upvaldftarray,upvalgwarray-upvaldftarray,'bx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
                P.plot(upconddftarray,upcondgwarray-upconddftarray,'rx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
                P.xlabel('DFT eigenvalues - spin UP (in eV)')
                P.ylabel('GW correction to the DFT eigenvalues (in eV)')
                P.figure(4,figsize=(csts.fig_width,csts.fig_height))
                P.hold(True)
                P.grid(True)
                P.plot(downvaldftarray,downvalgwarray-downvaldftarray,'bo',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
                P.plot(downconddftarray,downcondgwarray-downconddftarray,'ro',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
                P.xlabel('DFT eigenvalues - spin DOWN (in eV)')
                P.ylabel('GW correction to the DFT eigenvalues (in eV)')
            else:
                P.figure(2,figsize=(csts.fig_width,csts.fig_height))
                P.hold(True)
                P.grid(True)
                P.plot(valdftarray,valgwarray-valdftarray,'bx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
                P.plot(conddftarray,condgwarray-conddftarray,'rx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
                P.xlabel('DFT eigenvalues(in eV)')
                P.ylabel('GW correction to the DFT eigenvalues (in eV)')
            P.show()
            return
    if spinchoice == None or spinchoice == 'common':
        polyfitlist = list()
        if len(polyfit_degrees_up) == 1:
            print 'ERROR: making a fit with only one interval is not allowed ... exiting now'
            sys.exit()
        dftarray = N.append(valdftarray,conddftarray)
        gwarray = N.append(valgwarray,condgwarray)
        dftarray_list,gwarray_list = classify_eigenvalues(dftarray,energy_pivots_up,gwarray)
        for iinterval in range(len(polyfit_degrees_up)):
            tmpdftarray = dftarray_list[iinterval]
            tmpgwarray = gwarray_list[iinterval]
            if len(tmpdftarray) > 0:
                if limitpoints == 'least-squares' or (polyfit_degrees_up[iinterval] <= 0 and iinterval != len(polyfit_degrees_up)-1):
                    pfit = N.polyfit(tmpdftarray,tmpgwarray-tmpdftarray,N.abs(polyfit_degrees_up[iinterval]))
                elif limitpoints == 'least-squares_last-fixed':
                    if iinterval == len(polyfit_degrees_up)-1:
                        idftmin = N.argmin(tmpdftarray)
                        idftmax = N.argmax(tmpdftarray)
                        igwmin = N.argmin(tmpgwarray)
                        igwmax = N.argmax(tmpgwarray)
                        if idftmin == igwmin:
                            myimin = idftmin
                        else:
                            print 'COMMENT: the minimum for DFT and GW are not the same for band group #%s' %(iinterval+1)
                            print '         => the gw minimum is taken'
                            myimin = igwmin
                        pfit = polynd_a(tmpdftarray,tmpgwarray-tmpdftarray,polyfit_degrees_up[iinterval],indices=[myimin])
                    else:
                        pfit = N.polyfit(tmpdftarray,tmpgwarray-tmpdftarray,N.abs(polyfit_degrees_up[iinterval]))
                elif limitpoints == 'endpoints-fixed' or limitpoints == 'endpoints-fixed_last-flat':
                    idftmin = N.argmin(tmpdftarray)
                    idftmax = N.argmax(tmpdftarray)
                    igwmin = N.argmin(tmpgwarray)
                    igwmax = N.argmax(tmpgwarray)
                    if idftmin == igwmin:
                        myimin = idftmin
                    else:
                        print 'COMMENT: the minimum for DFT and GW are not the same for band group #%s' %(iinterval+1)
                        print '         => the gw minimum is taken'
                        myimin = igwmin
                    if iinterval == len(polyfit_degrees_up)-1:
                        if limitpoints == 'endpoints-fixed':
                            pfit = polynd_a(tmpdftarray,tmpgwarray-tmpdftarray,polyfit_degrees_up[iinterval],indices=[myimin])
                        elif limitpoints == 'endpoints-fixed_last-flat':
                            pfit = [N.polyval(polyfitlist[-1],energy_pivots_up[-1])]
                    else:
                        if idftmax == igwmax:
                            myimax = idftmax
                        else:
                            print 'COMMENT: the maximum for DFT and GW are not the same for band group #%s' %(iinterval+1)
                            print '         => the gw maximum is taken'
                            myimax = igwmax
                        pfit = polynd_ab(tmpdftarray,tmpgwarray-tmpdftarray,polyfit_degrees_up[iinterval],indices=[myimin,myimax])
            else:
                pfit = None
            polyfitlist.append(pfit)
        if smooth_end:
            if smooth_energy == None:
                smoothenergy = N.max(dftarray)
            else:
                smoothenergy = smooth_energy
            smoothdeltaenergy = None
            if smooth_delta_energy != None:
                smoothdeltaenergy = smooth_delta_energy
            oldpolyfitlist = list(polyfitlist)
            oldenergypivots = N.array(energy_pivots_up)
            energy_pivots_up,polyfitlist = smoothend(energy_pivots_up,polyfitlist,smoothenergy,delta_energy_ev=smoothdeltaenergy)
            dftarray_list,gwarray_list = classify_eigenvalues(dftarray,energy_pivots_up,gwarray)
        if plot_figures == 1:
            linspace_npoints = 200
            valpoly_x = N.linspace(N.min(valdftarray),N.max(valdftarray),linspace_npoints)
            condpoly_x = N.linspace(N.min(conddftarray),N.max(conddftarray),linspace_npoints)
            P.figure(3,figsize=(csts.fig_width,csts.fig_height))
            P.hold(True)
            P.grid(True)
            P.plot(valdftarray,valgwarray-valdftarray,'bx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            P.plot(conddftarray,condgwarray-conddftarray,'rx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            [x_min,x_max] = P.xlim()
            [y_min,y_max] = P.ylim()
            if smooth_end:
                x_max = energy_pivots_up[-1]+DELTA_ENERGY_END
            for iinterval in range(len(polyfitlist)):
                if iinterval == 0:
                    tmppoly_x = N.linspace(x_min,energy_pivots_up[iinterval],linspace_npoints)
                elif iinterval == len(polyfitlist)-1:
                    tmppoly_x = N.linspace(energy_pivots_up[iinterval-1],x_max,linspace_npoints)
                else:
                    tmppoly_x = N.linspace(energy_pivots_up[iinterval-1],energy_pivots_up[iinterval],linspace_npoints)
                if polyfitlist[iinterval] != None:
                    P.plot(tmppoly_x,N.polyval(polyfitlist[iinterval],tmppoly_x),'k',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            for ipivot in range(len(energy_pivots_up)):
                en = energy_pivots_up[ipivot]
                if polyfitlist[ipivot] != None and polyfitlist[ipivot+1] != None:
                    P.plot([en,en],[N.polyval(polyfitlist[ipivot],en),N.polyval(polyfitlist[ipivot+1],en)],'k-.',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            P.xlabel('DFT eigenvalues (in eV)')
            P.ylabel('GW correction to the DFT eigenvalues (in eV)')
            P.ylim([y_min,y_max])
            P.figure(4,figsize=(csts.fig_width,csts.fig_height))
            P.hold(True)
            P.grid(True)
            for iinterval in range(len(polyfitlist)):
                if polyfitlist[iinterval] != None:
                    P.plot(dftarray_list[iinterval],gwarray_list[iinterval]-dftarray_list[iinterval]-N.polyval(polyfitlist[iinterval],dftarray_list[iinterval]),'bx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            [x_min,x_max] = P.xlim()
            P.plot([x_min,x_max],[0,0],'k-')
            P.xlabel('DFT eigenvalues (in eV)')
            P.ylabel('Error in the fit (in eV)')
            P.show()
        return energy_pivots_up,polyfitlist
    elif spinchoice == 'separate':
        polyfitlist_up = list()
        polyfitlist_down = list()
        if len(polyfit_degrees_up) == 1 or len(polyfit_degrees_down) == 1:
            print 'ERROR: making a fit with only one interval is not allowed ... exiting now'
            sys.exit()
        updftarray = N.append(upvaldftarray,upconddftarray)
        upgwarray = N.append(upvalgwarray,upcondgwarray)
        downdftarray = N.append(downvaldftarray,downconddftarray)
        downgwarray = N.append(downvalgwarray,downcondgwarray)
        updftarray_list,upgwarray_list = classify_eigenvalues(updftarray,energy_pivots_up,upgwarray)
        downdftarray_list,downgwarray_list = classify_eigenvalues(downdftarray,energy_pivots_down,downgwarray)
        for iinterval in range(len(polyfit_degrees_up)):
            tmpdftarray = updftarray_list[iinterval]
            tmpgwarray = upgwarray_list[iinterval]
            if len(tmpdftarray) > 0:
                if limitpoints == 'least-squares' or polyfit_degrees_up[iinterval] <= 0:
                    pfit = N.polyfit(tmpdftarray,tmpgwarray-tmpdftarray,N.abs(polyfit_degrees_up[iinterval]))
                elif limitpoints == 'endpoints-fixed':
                    idftmin = N.argmin(tmpdftarray)
                    idftmax = N.argmax(tmpdftarray)
                    igwmin = N.argmin(tmpgwarray)
                    igwmax = N.argmax(tmpgwarray)
                    if idftmin == igwmin:
                        myimin = idftmin
                    else:
                        print 'COMMENT: the minimum for DFT and GW are not the same for band group #%s' %(iinterval+1)
                        print '         => the gw minimum is taken'
                        myimin = igwmin
                    if iinterval == len(polyfit_degrees_up)-1:
                        pfit = polynd_a(tmpdftarray,tmpgwarray-tmpdftarray,polyfit_degrees_up[iinterval],indices=[myimin])
                    else:
                        if idftmax == igwmax:
                            myimax = idftmax
                        else:
                            print 'COMMENT: the maximum for DFT and GW are not the same for band group #%s' %(iinterval+1)
                            print '         => the gw maximum is taken'
                            myimax = igwmax
                        pfit = polynd_ab(tmpdftarray,tmpgwarray-tmpdftarray,polyfit_degrees_up[iinterval],indices=[myimin,myimax])
            else:
                pfit = None
            polyfitlist_up.append(pfit)
        if smooth_end:
            if smooth_energy == None:
                smoothenergy = N.max(dftarray)
            else:
                smoothenergy = smooth_energy
            smoothdeltaenergy = None
            if smooth_delta_energy != None:
                smoothdeltaenergy = smooth_delta_energy
            oldpolyfitlist_up = list(polyfitlist_up)
            oldenergypivots_up = N.array(energy_pivots_up)
            energy_pivots_up,polyfitlist_up = smoothend(energy_pivots_up,polyfitlist_up,smoothenergy,delta_energy_ev=smoothdeltaenergy)
            updftarray_list,upgwarray_list = classify_eigenvalues(updftarray,energy_pivots_up,upgwarray)
        for iinterval in range(len(polyfit_degrees_down)):
            tmpdftarray = downdftarray_list[iinterval]
            tmpgwarray = downgwarray_list[iinterval]
            if len(tmpdftarray) > 0:
                if limitpoints == 'least-squares' or polyfit_degrees_down[iinterval] <= 0:
                    pfit = N.polyfit(tmpdftarray,tmpgwarray-tmpdftarray,N.abs(polyfit_degrees_down[iinterval]))
                elif limitpoints == 'endpoints-fixed':
                    idftmin = N.argmin(tmpdftarray)
                    idftmax = N.argmax(tmpdftarray)
                    igwmin = N.argmin(tmpgwarray)
                    igwmax = N.argmax(tmpgwarray)
                    if idftmin == igwmin:
                        myimin = idftmin
                    else:
                        print 'COMMENT: the minimum for DFT and GW are not the same for band group #%s' %(iinterval+1)
                        print '         => the gw minimum is taken'
                        myimin = igwmin
                    if iinterval == len(polyfit_degrees_down)-1:
                        pfit = polynd_a(tmpdftarray,tmpgwarray-tmpdftarray,polyfit_degrees_down[iinterval],indices=[myimin])
                    else:
                        if idftmax == igwmax:
                            myimax = idftmax
                        else:
                            print 'COMMENT: the maximum for DFT and GW are not the same for band group #%s' %(iinterval+1)
                            print '         => the gw maximum is taken'
                            myimax = igwmax
                        pfit = polynd_ab(tmpdftarray,tmpgwarray-tmpdftarray,polyfit_degrees_down[iinterval],indices=[myimin,myimax])
            else:
                pfit = None
            polyfitlist_down.append(pfit)
        if smooth_end:
            if smooth_energy == None:
                smoothenergy = N.max(dftarray)
            else:
                smoothenergy = smooth_energy
            smoothdeltaenergy = None
            if smooth_delta_energy != None:
                smoothdeltaenergy = smooth_delta_energy
            oldpolyfitlist_down = list(polyfitlist_down)
            oldenergypivots_down = N.array(energy_pivots_down)
            energy_pivots_down,polyfitlist_down = smoothend(energy_pivots_down,polyfitlist_down,smoothenergy,delta_energy_ev=smoothdeltaenergy)
            downdftarray_list,downgwarray_list = classify_eigenvalues(downdftarray,energy_pivots_down,downgwarray)
        if plot_figures == 1:
            linspace_npoints = 200
            upvalpoly_x = N.linspace(N.min(upvaldftarray),N.max(upvaldftarray),linspace_npoints)
            upcondpoly_x = N.linspace(N.min(upconddftarray),N.max(upconddftarray),linspace_npoints)
            downvalpoly_x = N.linspace(N.min(downvaldftarray),N.max(downvaldftarray),linspace_npoints)
            downcondpoly_x = N.linspace(N.min(downconddftarray),N.max(downconddftarray),linspace_npoints)
            P.figure(3,figsize=(csts.fig_width,csts.fig_height))
            P.hold(True)
            P.grid(True)
            P.plot(upvaldftarray,upvalgwarray-upvaldftarray,'bx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            P.plot(upconddftarray,upcondgwarray-upconddftarray,'rx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            [x_min,x_max] = P.xlim()
            [y_min,y_max] = P.ylim()
            #for iinterval in range(len(polyfit_degrees_up)):
            for iinterval in range(len(polyfitlist_up)):
                if iinterval == 0:
                    tmppoly_x = N.linspace(x_min,energy_pivots_up[iinterval],linspace_npoints)
                elif iinterval == len(polyfitlist_up)-1:
                    tmppoly_x = N.linspace(energy_pivots_up[iinterval-1],x_max,linspace_npoints)
                else:
                    tmppoly_x = N.linspace(energy_pivots_up[iinterval-1],energy_pivots_up[iinterval],linspace_npoints)
                if polyfitlist_up[iinterval] != None:
                    P.plot(tmppoly_x,N.polyval(polyfitlist_up[iinterval],tmppoly_x),'k',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            for ipivot in range(len(energy_pivots_up)):
                en = energy_pivots_up[ipivot]
                if polyfitlist_up[ipivot] != None and polyfitlist_up[ipivot+1] != None:
                    P.plot([en,en],[N.polyval(polyfitlist_up[ipivot],en),N.polyval(polyfitlist_up[ipivot+1],en)],'k-.',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            P.xlabel('DFT eigenvalues (in eV) - spin UP')
            P.ylabel('GW correction to the DFT eigenvalues (in eV)')
            P.ylim([y_min,y_max])
            P.figure(4,figsize=(csts.fig_width,csts.fig_height))
            P.hold(True)
            P.grid(True)
            P.plot(downvaldftarray,downvalgwarray-downvaldftarray,'bo',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            P.plot(downconddftarray,downcondgwarray-downconddftarray,'ro',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            [x_min,x_max] = P.xlim()
            [y_min,y_max] = P.ylim()
            for iinterval in range(len(polyfitlist_down)):
                if iinterval == 0:
                    tmppoly_x = N.linspace(x_min,energy_pivots_down[iinterval],linspace_npoints)
                elif iinterval == len(polyfitlist_down)-1:
                    tmppoly_x = N.linspace(energy_pivots_down[iinterval-1],x_max,linspace_npoints)
                else:
                    tmppoly_x = N.linspace(energy_pivots_down[iinterval-1],energy_pivots_down[iinterval],linspace_npoints)
                if polyfitlist_down[iinterval] != None:
                    P.plot(tmppoly_x,N.polyval(polyfitlist_down[iinterval],tmppoly_x),'k',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            for ipivot in range(len(energy_pivots_down)):
                en = energy_pivots_down[ipivot]
                if polyfitlist_down[ipivot] != None and polyfitlist_down[ipivot+1] != None:
                    P.plot([en,en],[N.polyval(polyfitlist_down[ipivot],en),N.polyval(polyfitlist_down[ipivot+1],en)],'k-.',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            P.xlabel('DFT eigenvalues (in eV) - spin DOWN')
            P.ylabel('GW correction to the DFT eigenvalues (in eV)')
            P.ylim([y_min,y_max])
            P.figure(5,figsize=(csts.fig_width,csts.fig_height))
            P.hold(True)
            P.grid(True)
            #for iinterval in range(len(polyfit_degrees_up)):
            for iinterval in range(len(polyfitlist_up)):
                if polyfitlist_up[iinterval] != None:
                    P.plot(updftarray_list[iinterval],upgwarray_list[iinterval]-updftarray_list[iinterval]-N.polyval(polyfitlist_up[iinterval],updftarray_list[iinterval]),'bx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            [x_min,x_max] = P.xlim()
            P.plot([x_min,x_max],[0,0],'k-')
            P.xlabel('DFT eigenvalues (in eV) - spin UP')
            P.ylabel('Error in the fit (in eV)')
            P.figure(6,figsize=(csts.fig_width,csts.fig_height))
            P.hold(True)
            P.grid(True)
            #for iinterval in range(len(polyfit_degrees_down)):
            for iinterval in range(len(polyfitlist_down)):
                if polyfitlist_down[iinterval] != None:
                    P.plot(downdftarray_list[iinterval],downgwarray_list[iinterval]-downdftarray_list[iinterval]-N.polyval(polyfitlist_down[iinterval],downdftarray_list[iinterval]),'bx',markersize=csts.markersize,markeredgewidth=csts.markeredgewidth)
            [x_min,x_max] = P.xlim()
            P.plot([x_min,x_max],[0,0],'k-')
            P.xlabel('DFT eigenvalues (in eV) - spin DOWN')
            P.ylabel('Error in the fit (in eV)')
            P.show()
        return energy_pivots_up,energy_pivots_down,polyfitlist_up,polyfitlist_down

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

parser = argparse.ArgumentParser(description='Tool for eigenvalue analysis')
parser.add_argument('-g','--graphical',help='use the graphical user interface',action='store_true')
parser.add_argument('-c','--command_line',help='use the command line interface',action='store_true')
parser.add_argument('files',help='files to be opened',nargs=2)
args = parser.parse_args()
args_dict = vars(args)
if args_dict['command_line'] and args_dict['graphical']:
    raise StandardError('Use either "-g/--graphical" or "-c/--command_line"')
elif args_dict['command_line']:
    use_gui = False
else:
    use_gui = False

if not use_gui:
    if args_dict['files']:
        if len(args_dict['files']) != 2:
            print 'ERROR: you should provide EIG.nc and _GW files ! exiting now ...'
            sys.exit()
        file_1 = args_dict['files'][0]
        file_2 = args_dict['files'][1]
        if file_1[-6:] == 'EIG.nc':
            eig_file = file_1
            if file_2[-3:] == '_GW':
                gw_file = file_2
            else:
                print 'ERROR: you should provide 1 _GW file with your EIG.nc file ! exiting now ...'
                sys.exit()
        elif file_1[-3:] == '_GW':
            gw_file = file_1
            if file_2[-6:] == 'EIG.nc':
                eig_file = file_2
            else:
                print 'ERROR: you should provide 1 EIG.nc file with your _GW file ! exiting now ...'
                sys.exit()
        else:
            print 'ERROR: you should provide 1 EIG.nc and 1 _GW files ! exiting now ...'
            sys.exit()
    else:
        print 'ERROR: you should provide EIG.nc and _GW files ! exiting now ...'
        sys.exit()

ec_dft = EigenvalueContainer(directory='.',filename=eig_file)
ec_gw = EigenvalueContainer(directory='.',filename=gw_file)
check_gw_vs_dft_parameters(ec_dft,ec_gw)
user_input = raw_input('Do you want to plot the figures ? (y/n)\n')
if user_input == 'y' or user_input == 'Y':
    plot_figures = 1
else:
    plot_figures = 0
user_input = raw_input('Enter the index of the valence band maximum :\n')
vbm_index = N.int(user_input)
user_input = raw_input('Do you want the script to automatically find groups of bands (y/n) ?\n')
if user_input == 'y':
    user_input = raw_input('Enter the name of the bandstructure file used to find groups of bands\n(<ENTER> for finding groups of bands on the regular grid -- file "%s" ... not recommended)\n' %eig_file)
    if user_input == '':
        if ec_dft.nsppol > 1:
            energy_pivots_up_ha,energy_pivots_down_ha = ec_dft.find_band_groups(spinchoice='separate')
            energy_pivots_up = csts.hartree2ev*energy_pivots_up_ha
            energy_pivots_down = csts.hartree2ev*energy_pivots_down_ha
        else:
            energy_pivots = csts.hartree2ev*ec_dft.find_band_groups()
    else:
        if ec_dft.nsppol > 1:
            energy_pivots_up_ha,energy_pivots_down_ha = ec_dft.find_band_groups(bandstructure_file=user_input,spinchoice='separate')
            energy_pivots_up = csts.hartree2ev*energy_pivots_up_ha
            energy_pivots_down = csts.hartree2ev*energy_pivots_down_ha
        else:
            energy_pivots = csts.hartree2ev*ec_dft.find_band_groups(bandstructure_file=user_input)
    if ec_dft.nsppol > 1:
        nfittingintervals_up = len(energy_pivots_up)+1
        nfittingintervals_down = len(energy_pivots_down)+1
    else:
        nfittingintervals = len(energy_pivots)+1
else:
    if plot_figures == 1:
        plot_gw_vs_dft_eig(ec_dft,ec_gw,vbm_index,spinchoice='separate')
    if ec_dft.nsppol == 1:
        user_input = raw_input('How many fitting intervals do you want ? (default is 2 : valence/conduction => press <ENTER>)\n')
        if user_input == '':
            nfittingintervals = 2
            energy_pivots = N.zeros(nfittingintervals-1,N.float)
            energy_pivots[0] = csts.hartree2ev*(N.min(ec_dft.eigenvalues[:,:,vbm_index])+N.max(ec_dft.eigenvalues[:,:,vbm_index-1]))/2
        else:
            nfittingintervals = N.int(user_input)
            energy_pivots = N.zeros(nfittingintervals-1,N.float)
            user_input = raw_input('Enter the %s energy "pivots" that splits the dft eigenvalues in %s fitting intervals (in eV) :\n' %(nfittingintervals-1,nfittingintervals))
            energy_pivots = N.array(user_input.split(),N.float)
            if len(energy_pivots) != nfittingintervals-1:
                print 'ERROR: you asked %s fitting intervals and provided %s energy "pivots".' %(nfittingintervals,len(energy_pivots))
                print '       you should provide %s energy "pivots" ... exiting now' %(nfittingintervals-1)
                sys.exit()
            for ienergy in range(1,len(energy_pivots)):
                if energy_pivots[ienergy] <= energy_pivots[ienergy-1]:
                    print 'ERROR: the energy pivots have to be entered increasingly'
                    print '       you should provide energy "pivots" with increasing energies ... exiting now'
                    sys.exit()
    elif ec_dft.nsppol == 2:
        user_input = raw_input('How many fitting intervals do you want for spin up ? (default is 2 : valence/conduction => press <ENTER>)\n')
        if user_input == '':
            nfittingintervals_up = 2
            energy_pivots_up = N.zeros(nfittingintervals_up-1,N.float)
            energy_pivots_up[0] = csts.hartree2ev*(N.min(ec_dft.eigenvalues[0,:,vbm_index])+N.max(ec_dft.eigenvalues[0,:,vbm_index-1]))/2
        else:
            nfittingintervals_up = N.int(user_input)
            energy_pivots_up = N.zeros(nfittingintervals_up-1,N.float)
            user_input = raw_input('Enter the %s energy "pivots" that splits the dft eigenvalues (spin up) in %s fitting intervals (in eV) :\n' %(nfittingintervals_up-1,nfittingintervals_up))
            energy_pivots_up = N.array(user_input.split(),N.float)
            if len(energy_pivots_up) != nfittingintervals_up-1:
                print 'ERROR: you asked %s fitting intervals and provided %s energy "pivots".' %(nfittingintervals_up,len(energy_pivots_up))
                print '       you should provide %s energy "pivots" ... exiting now' %(nfittingintervals_up-1)
                sys.exit()
            for ienergy in range(1,len(energy_pivots_up)):
                if energy_pivots_up[ienergy] <= energy_pivots_up[ienergy-1]:
                    print 'ERROR: the energy pivots have to be entered increasingly'
                    print '       you should provide energy "pivots" with increasing energies ... exiting now'
                    sys.exit()
        user_input = raw_input('How many fitting intervals do you want for spin down ? (default is 2 : valence/conduction => press <ENTER>)\n')
        if user_input == '':
            nfittingintervals_down = 2
            energy_pivots_down = N.zeros(nfittingintervals_down-1,N.float)
            energy_pivots_down[0] = csts.hartree2ev*(N.min(ec_dft.eigenvalues[0,:,vbm_index])+N.max(ec_dft.eigenvalues[0,:,vbm_index-1]))/2
        else:
            nfittingintervals_down = N.int(user_input)
            energy_pivots_down = N.zeros(nfittingintervals_down-1,N.float)
            user_input = raw_input('Enter the %s energy "pivots" that splits the dft eigenvalues (spin down) in %s fitting intervals (in eV) :\n' %(nfittingintervals_down-1,nfittingintervals_down))
            energy_pivots_down = N.array(user_input.split(),N.float)
            if len(energy_pivots_down) != nfittingintervals_down-1:
                print 'ERROR: you asked %s fitting intervals and provided %s energy "pivots".' %(nfittingintervals_down,len(energy_pivots_down))
                print '       you should provide %s energy "pivots" ... exiting now' %(nfittingintervals_down-1)
                sys.exit()
            for ienergy in range(1,len(energy_pivots_down)):
                if energy_pivots_down[ienergy] <= energy_pivots_down[ienergy-1]:
                    print 'ERROR: the energy pivots have to be entered increasingly'
                    print '       you should provide energy "pivots" with increasing energies ... exiting now'
                    sys.exit()
if ec_dft.nsppol > 1:
    print 'Script will use the following energy pivots for the interpolation (spin up)'
    print energy_pivots_up
    user_input = raw_input('Enter the degree of polynomials used to fit the GW corrections (spin up) \nfor each interval (%s values, default is 3rd order polynomials with "fixed points" for each group of bands => press <ENTER>)\n or enter "options" to enter specific options' %nfittingintervals_up)
    if user_input == '':
        polyfit_degrees_up = 3*N.ones(nfittingintervals_up,N.int)
        option_limit_points = 'endpoints-fixed'
    elif user_input == 'options':
        print 'ERROR: this option is not yet coded ... exit'
        sys.exit()
    else:
        polyfit_degrees_up = N.array(user_input.split(),N.int)
        option_limit_points = 'endpoints-fixed'
    print 'Script will use the following energy pivots for the interpolation (spin down)'
    print energy_pivots_down
    user_input = raw_input('Enter the degree of polynomials used to fit the GW corrections (spin down) \nfor each interval (%s values, default is 3rd order polynomials with "fixed points" for each group of bands => press <ENTER>)\n or enter "options" to enter specific options' %nfittingintervals_down)
    if user_input == '':
        polyfit_degrees_down = 3*N.ones(nfittingintervals_down,N.int)
        option_limit_points = 'endpoints-fixed'
    elif user_input == 'options':
        print 'ERROR: this option is not yet coded ... exit'
        sys.exit()
    else:
        polyfit_degrees_down = N.array(user_input.split(),N.int)
        option_limit_points = 'endpoints-fixed'
    new_energy_pivots_up,new_energy_pivots_down,polyfit_list_up,polyfit_list_down = plot_gw_vs_dft_eig(ec_dft,ec_gw,vbm_index,energy_pivots_up=energy_pivots_up,energy_pivots_down=energy_pivots_down,polyfit_degrees_up=polyfit_degrees_up,polyfit_degrees_down=polyfit_degrees_down,limitpoints=option_limit_points,spinchoice='separate')
else:
    print 'Script will use the following energy pivots for the interpolation (same for all spins)'
    print energy_pivots
    user_input = raw_input('Enter the degree of polynomials used to fit the GW corrections \nfor each interval (%s values, default is 3rd order polynomials with "fixed points" for each group of bands => press <ENTER>)\n or enter "options" to enter specific options' %nfittingintervals)
    if user_input == '':
        polyfit_degrees = 3*N.ones(nfittingintervals,N.int)
        option_limit_points = 'endpoints-fixed'
    elif user_input == 'options':
        print 'ERROR: this option is not yet coded ... exit'
        sys.exit()
    else:
        if user_input.split()[-1] == 'x':
            tmp = user_input.split()
            tmp[-1] = '0'
            polyfit_degrees = N.array(tmp,N.int)
            option_limit_points = 'endpoints-fixed_last-flat'
        else:
            polyfit_degrees = N.array(user_input.split(),N.int)
            option_limit_points = 'endpoints-fixed'
    user_input = raw_input('Enter specific options for the end of the polyfit ? (y/n) [<ENTER> to continue without entering specific options]')
    if user_input == 'y':
        user_input = raw_input('Enter the end smooth energy (to be documented ...) : ')
        smoothenergy = N.float(user_input)
        user_input = raw_input('Enter the end smooth delta energy (to be documented ...) [<ENTER> for default]: ')
        if user_input == '':
            smoothdeltaenergy = None
        else:
            smoothdeltaenergy = N.float(user_input)
        new_energypivots,polyfit_list = plot_gw_vs_dft_eig(ec_dft,ec_gw,vbm_index,energy_pivots_up=energy_pivots,polyfit_degrees_up=polyfit_degrees,limitpoints=option_limit_points,smooth_end=True,smooth_energy=smoothenergy,smooth_delta_energy=smoothdeltaenergy)
    else:
        new_energypivots,polyfit_list = plot_gw_vs_dft_eig(ec_dft,ec_gw,vbm_index,energy_pivots_up=energy_pivots,polyfit_degrees_up=polyfit_degrees,limitpoints=option_limit_points)

if ec_dft.nsppol > 1:
    print polyfit_list_up
    print polyfit_list_down
else:
    print polyfit_list
    write_polyfit('mytest.pfitlist',new_energypivots,polyfit_list)

gw_interpolate = False
user_input = raw_input('Do you want to make an interpolated _GW file ? (y/n)\n')
if user_input == 'y' or user_input == 'Y':
    gw_interpolate = True
if gw_interpolate:
    nc_eig_file = raw_input('Enter the name of the EIG.nc file you want to extrapolate to GW :\n')
    new_ec_dft = EigenvalueContainer(directory='.',filename=nc_eig_file)
    user_input = raw_input('For which "bdgw"\'s do you want the interpolation (indices of \nthe smallest \
                            valence and largest conduction bands) ? bdgw(1)<=vbm_index<bdgw(2)<=%s \
                            (usually "1 something")\n' %N.min(new_ec_dft.bd_indices[:,:,1]))
    if user_input == '':
        bdgw_interpolated = N.array([1,N.min(new_ec_dft.bd_indices[:,:,1])])
    else:
        bdgw_interpolated = N.array(user_input.split(),N.int)
    
    filename = '%s_polyfit_GW' %(nc_eig_file)
    new_ec_dft.pfit_gw_file_write(polyfit_list,filename=filename,bdgw=bdgw_interpolated,energy_pivots=new_energypivots,gwec=ec_gw)
user_input = raw_input('Do you want to make an interpolated bandstructure file ? (y/n)\n')
if user_input == 'y' or user_input == 'Y':
    nc_eig_file = raw_input('Enter the name of the bandstructure EIG.nc file you want to extrapolate to GW :\n')
    new_ec_dft = EigenvalueContainer(directory='.',filename=nc_eig_file)

    gvectors = get_gvectors()

    if ec_dft.nsppol > 1:
        gw_eigenvalues = new_ec_dft.pfit_gw_eigenvalues_ha(polyfit_list_up,energy_pivots_up=energy_pivots_up,polyfitlist_down=polyfit_list_down,energy_pivots_down=energy_pivots_down,ecgw=ec_gw)
        new_ec_dft.eigenvalues = gw_eigenvalues
    else:
        #gw_eigenvalues = new_ec_dft.pfit_gw_eigenvalues_ha(polyfit_list,energy_pivots_up=energy_pivots,ecgw=ec_gw)
        gw_eigenvalues = new_ec_dft.pfit_gw_eigenvalues_ha(polyfit_list,energy_pivots_up=new_energypivots,ecgw=None)
        new_ec_dft.eigenvalues = gw_eigenvalues
    new_ec_dft.set_kpoint_sampling_type('Bandstructure')
    new_ec_dft.find_special_kpoints(gvectors)
    
    print 'Number of bands in the file : %s' %(N.shape(new_ec_dft.eigenvalues)[2])
    test = raw_input('Enter the number of bands to be plotted (<ENTER> : %s) : \n' %(N.shape(new_ec_dft.eigenvalues)[2]))
    if test == '':
        nbd_plot = N.shape(new_ec_dft.eigenvalues)[2]
    else:
        nbd_plot = N.int(test)
    if nbd_plot > N.shape(new_ec_dft.eigenvalues)[2]:
        print 'ERROR: the number of bands to be plotted is larger than the number available ... exit'
        sys.exit()
    
    new_ec_dft.special_kpoints_names = ['']*len(new_ec_dft.special_kpoints_indices)
    for ii in range(len(new_ec_dft.special_kpoints_indices)):
        new_ec_dft.special_kpoints_names[ii] = 'k%s' %(ii+1)
    print 'List of special kpoints :'
    for ii in range(len(new_ec_dft.special_kpoints_indices)):
        spkpt = new_ec_dft.kpoints[new_ec_dft.special_kpoints_indices[ii]]
        print ' Kpoint %s : %s %s %s' %(ii+1,spkpt[0],spkpt[1],spkpt[2])
    print 'Enter the name of the %s special k-points :' %(len(new_ec_dft.special_kpoints_indices))
    test = raw_input('')
    if len(test.split()) == len(new_ec_dft.special_kpoints_indices):
        for ii in range(len(new_ec_dft.special_kpoints_indices)):
            new_ec_dft.special_kpoints_names[ii] = test.split()[ii]
    
    test = raw_input('Enter base name for bandstructure file : \n')
    new_ec_dft.write_bandstructure_to_file('%s.bandstructure' %test)
    
    
    P.figure(1,figsize=(3.464,5))
    P.hold('on')
    P.grid('on')
    P.xticks(N.take(new_ec_dft.kpoint_reduced_path_values,N.array(new_ec_dft.special_kpoints_indices,N.int)),new_ec_dft.special_kpoints_names)
    for iband in range(nbd_plot):
        if new_ec_dft.nsppol == 1:
            P.plot(new_ec_dft.kpoint_reduced_path_values,new_ec_dft.eigenvalues[0,:,iband]*csts.hartree2ev,'k-',linewidth=2)
        else:
            P.plot(new_ec_dft.kpoint_reduced_path_values,new_ec_dft.eigenvalues[0,:,iband]*csts.hartree2ev,'k-',linewidth=2)
            P.plot(new_ec_dft.kpoint_reduced_path_values,new_ec_dft.eigenvalues[1,:,iband]*csts.hartree2ev,'r-',linewidth=2)
    P.show()
