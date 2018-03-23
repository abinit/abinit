#!/usr/bin/env python
# Author: Samuel Ponc\'e
# Date: 26/06/2013
# Script to create an corrected electronic bandstructure with lifetime broadening

try:
  from rf_final import system
except ImportError:
  import warnings
  warnings.warn("The system module is missing!")
  raise
try:
  import numpy as N
except ImportError:
  import warnings
  warnings.warn("The numpy module is missing!")
  raise
from numpy import zeros
try:
  import netCDF4 as nc
except ImportError:
  import warnings
  warnings.warn("The netCDF4 module is missing!")
  raise
import matplotlib.pyplot as P
from numpy.linalg import inv
from scipy.interpolate import spline


#############
# Constants #
#############
class VariableContainer:pass
csts = VariableContainer()
csts.hartree2ev = N.float(27.211396132)
csts.ev2hartree = N.float(1/csts.hartree2ev)
csts.sqrtpi = N.float(N.sqrt(N.pi))
csts.invsqrtpi = N.float(1/csts.sqrtpi)
csts.TOLKPTS = N.float(0.00001)

def angle_kpt(vector1,vector2):
  dot_prod = N.dot(vector1,vector2.T)
  arg = dot_prod/(N.linalg.norm(vector1)*N.linalg.norm(vector2))
  if arg > 1: arg=1
  if arg < -1: arg=-1
  theta = N.arccos(arg)/(N.pi*180)
  return theta

class find_special_kpt:
  def __init__(self,kpt):
    nkpt = len(kpt[:,0])
    self.special_kpoints = list()
    self.position_special = []
    angle = 0
    ref = 0
    for ikpt in N.arange(nkpt):
      if ikpt == 0:
        self.special_kpoints.append(kpt[ikpt,:])
        self.position_special.append(ikpt)
      elif ikpt == 1:
        vector2 = kpt[ikpt,:]-kpt[ikpt-1,:]
      elif ikpt == nkpt-1:
        self.special_kpoints.append(kpt[ikpt,:])
        self.position_special.append(ikpt)
      else:
        vector1 = kpt[ikpt,:]-kpt[ikpt-1,:]
        if N.allclose(kpt[ikpt-1,:],[3.0/8,3.0/8,3.0/4],1e-8):
          self.special_kpoints.append(kpt[ikpt-1,:])
          self.position_special.append(ikpt-1)
        else:
          if angle_kpt(vector1,vector2) < csts.TOLKPTS:
            pass
          else:
            self.special_kpoints.append(kpt[ikpt-1,:])
            self.position_special.append(ikpt-1)
        vector2 = vector1
     

# Interaction with the user
print """
  _____  _       _          ______ _____        ____   _____ 
 |  __ \| |     | |        |  ____|  __ \      |  _ \ / ____|
 | |__) | | ___ | |_ ______| |__  | |__) |_____| |_) | (___  
 |  ___/| |/ _ \| __|______|  __| |  ___/______|  _ < \___ \ 
 | |    | | (_) | |_       | |____| |          | |_) |____) |
 |_|    |_|\___/ \__|      |______|_|          |____/|_____/
"""
print """ 
This script allows you to plot an electronic bandstructure from en _EP.nc 
file. If you computed the file including lifetime you can also plot 
"fat-band".
"""
# Read the input file
user_input = raw_input('Enter name of the _EP.nc file\n')
EP_file = user_input.strip()
EP = system(directory='.',filename=EP_file)

gprimd = inv(EP.rprimd)
full_kpt = N.matrix(EP.kpt)*gprimd
#print full_kpt

special = find_special_kpt(EP.kpt)

print "The special k-points are:"
for ikpt in N.arange(len(special.special_kpoints)):
  print str(special.special_kpoints[ikpt])
user_input = raw_input('Enter the name of the '+str(len(special.special_kpoints))+' special k-points \n')

special_name = user_input.split()
for ii in N.arange(len(special_name)):
  special_name[ii] = "$\mathbf{"+str(special_name[ii])+"}$"

user_input = raw_input('Enter energy limit for the plot in eV: e.g. -5 10\n')
user_tmp = user_input.split()
if len(user_tmp) != 2:
  raise Exception("You should provide only 2 numbers")
else: # Append and TRIM the input string with STRIP
  lower = N.float(user_tmp[0])
  upper = N.float(user_tmp[1])

if EP.ntemp > 0:
  print "Enter the temperature at which you want to do the Bandstructure plot (in K)"
  user_input = raw_input('The possible temperature are:\n'+str(EP.temp[:])+'\n')
  if len(user_input.split()) != 1:
    raise Exception("You should provide only 1 number")
  else: # Append and TRIM the input string with STRIP
    temp = N.float(user_input)
    for itemp in N.arange(EP.ntemp):
      if N.allclose(temp,EP.temp[itemp],csts.TOLKPTS):
        temp_index = itemp
        break  

# Compute the first conduction band
for iband in N.arange(EP.nband):
  if (EP.occ[0,0,iband]==0):
    cond = iband
    print "First Cond band ",cond
    break

# Definie the Fermi level as the highest band before 0 occupations
fermi = -1000000.0
for ikpt in N.arange(EP.nkpt):
  for iband in N.arange(cond):
    if (EP.eigenvalues[0,ikpt,iband] > fermi):
      fermi = EP.eigenvalues[0,ikpt,iband]


xspan = N.arange(0,EP.nkpt,1)

xfine = N.arange(0,(EP.nkpt-1)+0.2,0.2)

P.figure(1,figsize=(5.196,7.5))
P.rc('text',usetex = True)
P.hold('on')
P.grid('on')
yprops = dict(rotation=0, horizontalalignment='right',verticalalignment='center',x=-0.01)
ax = P.gca()
#ax.set_xticks(special.position_special,special_name)

P.xticks(special.position_special,special_name,fontsize=16)
P.yticks(fontsize=16)
ylabel = ax.set_ylabel('Energy [eV]',fontsize=16,**yprops)
ax.yaxis.set_label_coords(0.17, 1.05)
ax.set_ylim([lower,upper])

# Plot Eigenenergies
for iband in N.arange(EP.nband):
  if iband < cond:
    eigen =(EP.eigenvalues[0,:,iband]-fermi)*csts.hartree2ev
    eig_fine = spline(xspan,eigen,xfine)
#    P.plot(xfine,eig_fine,linewidth=2,color='#0099FF')    
    P.plot(xfine,eig_fine,linewidth=2,color='k')    
  else:
    eig_fine = spline(xspan,(EP.eigenvalues[0,:,iband]-fermi)*csts.hartree2ev,xfine)
    P.plot(xfine,eig_fine,linewidth=2,color='k')
#    P.plot(xfine,eig_fine,linewidth=2,color='#FF6600')

# Plot renormalization
if EP.ntemp > 0: 
  for iband in N.arange(EP.nband):
    if iband < cond:
      #renormalization = real part of zpr
      renorm = (EP.eigenvalues[0,:,iband]-fermi + EP.zpm[temp_index,0,:,iband,0])*csts.hartree2ev
      bandwith = EP.zpm[temp_index,0,:,iband,1]*csts.hartree2ev
      renorm_fine = spline(xspan,renorm,xfine)
      bandwith_fine = spline(xspan,bandwith,xfine)
      P.fill_between(xfine,renorm_fine+bandwith_fine/2, renorm_fine-bandwith_fine/2, alpha=.3,color='b')
      P.plot(xfine,renorm_fine,color='b',linestyle='--',linewidth=2)
    else:
      #renormalization = real part of zpr
      renorm = (EP.eigenvalues[0,:,iband]-fermi + EP.zpm[temp_index,0,:,iband,0])*csts.hartree2ev
      bandwith = EP.zpm[temp_index,0,:,iband,1]*csts.hartree2ev
      renorm_fine = spline(xspan,renorm,xfine)
      bandwith_fine = spline(xspan,bandwith,xfine)
      P.fill_between(xfine,renorm_fine+bandwith_fine/2, renorm_fine-bandwith_fine/2, alpha=.3, color='r')
      P.plot(xfine,renorm_fine,color='r',linestyle='--',linewidth=2)
  bbox_props = dict(boxstyle="square", fc="w", ec="0.5", alpha=1.0)
  ax.text(EP.nkpt,upper, "Temperature: "+str(temp)+" K", ha="right", va="top", size=16,
        bbox=bbox_props)

P.show()










