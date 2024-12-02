#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
author: Nicholas Pike
Email : Nicholas.pike@smn.uio.no

Purpose: Use the extracted elastic tensor to calculate the Debye temperature 
         and other thermodynamic properties that can be derived from the elastic
         properties.
         
Reference: Phys. Rev. B 95, 155206 (2017) and (more importantly) 
           J. Phys. Chem. Solids 24, 909 (1963)
           
"""
#import useful modules
import os
import sys
import linecache
import numpy as np

#user defined variables
amu_to_kg    = 1.66054E-27
ang_to_meter = 1E-10
h            = 6.626070040E-34
kb           = 1.38064852E-23
Na           = 6.0221409E23
kbar_to_GPa  = 0.01
m_to_cm      = 100.0
kgm3_to_cm3  = 1000.0


def find_file(filetype):
    """ 
    Author: Nicholas Pike
    Email : Nicholas.Pike@smn.uio.no
    
    Purpose: Determine if the files from your calculation is found in 
             the directory. It does this by looking for a specific set of data 
             in the file
    """
    _filename = ''
    
    if filetype == 'OUTCAR':
        #look for the file in the directory by first listing the files
        files = [f for f in os.listdir('.') if os.path.isfile(f) and not f.endswith(".png") and not f.endswith('.py') and not f.endswith('~') and not f.endswith('.nc') and not f == 'derived_elastic' ]
        
        #now loop through the files and search for keywords
        found_vasp    = False
        found_elastic = False
        
        for f in files:
            with open(f,'r') as d:
                for i, line in enumerate(d):
                    if 'vasp' in line:
                        found_vasp = True
                        _filename = f
                    elif 'ELASTIC MODULI  (kBar)' in line:
                        found_elastic = True
                        
        if found_elastic == False and found_vasp == True:
            print('ERROR: While a VASP output file was found, it did not contain the Elastic tensor.'
                  '\n Check calculation before trying again.')
            sys.exit()
        elif found_elastic == False and found_vasp == False:
            print('ERROR: No VASP OUTCAR file was found in the directory.'
                  '\n Check the directory before trying again')
            sys.exit()
    
    elif filetype == 'ANADDB':  
        #look for the file in the directory by first listing the files
        files = [f for f in os.listdir('.') if os.path.isfile(f) and not f.endswith(".png") and not f.endswith('.py') and not f.endswith('~') and not f.endswith('.nc') and not f == 'derived_elastic' ]
        #now loop through the files and search for keywords
        found_anaddb  = False
        found_elastic = False
        found_ddb     = False
        file_output   = ''
        file_ddb      = ''
        
        for f in files:
            with open(f,'r') as d:
                for i, line in enumerate(d):
                    if ' ANADDB comes with ABSOLUTELY NO WARRANTY.' in line:
                        found_anaddb = True
                        file_output = f
                    elif 'Elastic Tensor (relaxed ion) (unit:10^2GP)' in line:
                        found_elastic = True
                    elif ' **** DERIVATIVE DATABASE ****' in line:
                        found_ddb  = True
                        file_ddb = f
                        
        if found_elastic == False and found_anaddb == True:
            print('ERROR: While a ANADDB output file was found, it did not contain the Elastic tensor.'
                  '\n Check calculation before trying again.')
            sys.exit()
        elif found_elastic == False and found_anaddb == False:
            print('ERROR: No ANADDB output file was found in the directory.'
                  '\n Check the directory before trying again')
            sys.exit()   
        elif found_ddb == False:
            print('ERROR: No DDB file was found in the director.'
                  '\n Check the directory before trying again.')
            sys.exit()
        
        _filename = [file_output,file_ddb]      

    return _filename

def gather_from_OUTCAR(filename):
    """
    Author: Nicholas Pike
    Email : Nicholas.pike@smn.uio.no
    
    Purpose: Gather data from the output of a VASP calculation for the elastic tensor
    
    OUTPUT: Returns the system name, elastic tensor, ion mass, volume in (m^3), 
            and the nultiplicity of each atom in the unit cell
    
    """
    system_name = ''
    elasten     = np.empty(shape=(6,6))
    ionmass     = []
    vol         = 0.0
    #open file and extract data
    with open(filename,'r') as d:
       for i, line in enumerate(d):
          if 'POSCAR =' in line:
              data = linecache.getline(filename,i+1).split()
              system_name = data[2:] 
          elif 'SYMMETRIZED ELASTIC MODULI (kBar)' in line:
              for j in range(6):
                  data = linecache.getline(filename,i+j+4).split()
                  elasten[j][:] = data[1:]
          elif 'Mass of Ions in am' in line:
              data = linecache.getline(filename,i+2).split()
              for j in range(2,len(data)):
                  ionmass = np.append(float(data[j]),ionmass)
          elif 'volume of cell :' in line:
              data = linecache.getline(filename,i+1).split()
              vol  = float(data[4]) 
          
    atommult = []
    for i in range(len(system_name)):
        namelist = list(system_name[i])
        number   = namelist[len(namelist)-1]
        atommult = np.append(float(number),atommult)
         
    data_array = [system_name,elasten,ionmass,vol*ang_to_meter**3,atommult]
    
    return data_array

def gather_from_ANADDB(fileout,fileddb):
    """
    Author: Nicholas Pike
    Email : Nicholas.pike@smn.uio.no
    
    Purpose: Gather data from the output of an ANADDB calculation for the elastic tensor
    
    OUTPUT: Returns the system name, elastic tensor, ion mass, volume in (m^3), 
            and the nultiplicity of each atom in the unit cell
    
    """
    elasten     = np.empty(shape=(6,6))
    ionmass     = []
    vol         = 0.0
    ntypat      = 0.0
    typat       = []
    
    #open file and extract data
    with open(fileout,'r') as d:
       for i, line in enumerate(d):
          if 'Elastic Tensor (relaxed ion) (unit:10^2GP):' in line:
              for j in range(6):
                  data = linecache.getline(fileout,i+j+4).split()
                  elasten[j][:] = data[:]
          elif ' Unit cell volume ucvol= ' in line:
              data = linecache.getline(fileout,i+1).split()
              vol  = float(data[4]) 
                       
    with open(fileddb,'r') as d:
        for i, line in enumerate(d):
            if 'amu' in line:
               data = linecache.getline(fileddb,i+1).split()
               ionmass = data[1:]
            elif 'ntypat' in line:
                data = linecache.getline(fileddb,i+1).split()
                ntypat = data[1]   
            elif 'typat' in line:
                data = linecache.getline(fileddb,i+1).split()
                typat = data[1:]                   
    
    #replace D with E
    for i in range(len(ionmass)):
        ionmass[i]=float(ionmass[i].replace('D','E'))
              
    atommult = []
    u,numtyat = np.unique(typat, return_counts= True)
    for i in range(int(ntypat)):
        atommult = np.append(numtyat[i],atommult)
             
    data_array = ['',elasten*1000.0,ionmass,vol*ang_to_meter**3,atommult]
    return data_array

def sound_velocities(data):
    """
    Author: Nicholas Pike
    Email : Nicholas.pike@smn.uio.no
    
    Purpose: Calculate the longitudional, shear and average sound velocities
    """
    #declare variables
    vl      = 0
    vs      = 0
    va      = 0
    density = 0
    B       = 0
    BR      = 0
    G       = 0
    GR      = 0
    YXX     = 0
    YYY     = 0
    YZZ     = 0
    nuxy    = 0
    nuyz    = 0
    nuxz    = 0
    
    #calculate compliance tensor
    s = np.linalg.inv(data[1])
    
    #calculate density
    masstot = 0
    for i in range(len(data[2])):
        masstot +=data[2][i]*data[4][i] #data[2] is the mass data[4] is the multiplicity
    density = masstot*amu_to_kg/data[3]/kgm3_to_cm3 #data[3] is the volume in g/cm3
    
    #Calculate bulk modulus and shear modulus
    #note: data[1] is the elastic tensor
    B  = 1.0/9.0*((data[1][0][0]+data[1][1][1]+data[1][2][2])+2.0*(data[1][0][1]+data[1][0][2]+data[1][1][2]))
    BR = 1.0/((s[0,0]+s[1,1]+s[2,2])+2.0*(s[0,1]+s[1,2]+s[0,2]))
    
    G  = 1.0/15.0*((data[1][0][0]+data[1][1][1]+data[1][2][2])-(data[1][0][1]+data[1][0][2]+data[1][1][2])+3.0*(data[1][3][3]+data[1][4][4]+data[1][5][5]))
    GR = 15.0/(4.0*(s[0,0]+s[1,1]+s[3,3])-4.0*(s[0,1]+s[1,2]-s[0,2])+3.0*(s[3,3]+s[4,4]+s[5,5])) 
    
    #calculate Youngs modulus
    YXX = 1.0/s[0,0]
    YYY = 1.0/s[1,1]
    YZZ = 1.0/s[2,2]
    
    #calculate vl
    if B + 4.0/3.0*G <= 0:
        vl = 0
    else:  
    	vl  = np.sqrt((B+4.0/3.0*G)*kbar_to_GPa/(density/kgm3_to_cm3))*m_to_cm # convert kbar to GPa and then g/cm3 to kg /m3
    
    #calculate vs
    if G <= 0:
	G = 0
    else:
        vs = np.sqrt(G*kbar_to_GPa/(density/kgm3_to_cm3))*m_to_cm # convert kbar to GPa and then g/cm3 to kg /m3
    
    #calculate va
    va = (1.0/3.0*((1.0/vl**3)+(2.0/vs**3)))**(-1.0/3.0)
    
    #calculate poisson ratio
    nuxy = -s[0,1]*YXX
    nuyz = -s[1,2]*YYY
    nuxz = -s[0,2]*YZZ
    
    #print results 
    f1 = open('derived_elastic','a')
    f1.write('Elastic Tensor\n')
    f1.write('Direction XX\t  YY\t   ZZ\t    XY\t     YZ\t      ZX\n' )
    f1.write('XX \t %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n'%(data[1][0][0],data[1][0][1],data[1][0][2],data[1][0][3],data[1][0][4],data[1][0][5] ))
    f1.write('YY \t %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n'%(data[1][1][0],data[1][1][1],data[1][1][2],data[1][1][3],data[1][1][4],data[1][1][5] ))
    f1.write('ZZ \t %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n'%(data[1][2][0],data[1][2][1],data[1][2][2],data[1][2][3],data[1][2][4],data[1][2][5] ))
    f1.write('XY \t %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n'%(data[1][3][0],data[1][3][1],data[1][3][2],data[1][3][3],data[1][3][4],data[1][3][5] ))
    f1.write('YZ \t %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n'%(data[1][4][0],data[1][4][1],data[1][4][2],data[1][4][3],data[1][4][4],data[1][4][5] ))
    f1.write('ZX \t %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n'%(data[1][5][0],data[1][5][1],data[1][5][2],data[1][5][3],data[1][5][4],data[1][5][5] ))
    f1.write('\n')
    f1.write('Compliance Tensor\n')
    f1.write('Direction XX\t  YY\t   ZZ\t    XY\t     YZ\t      ZX\n' )
    f1.write('XX \t %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n'%(s[0][0],s[0][1],s[0][2],s[0][3],s[0][4],s[0][5] ))
    f1.write('YY \t %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n'%(s[1][0],s[1][1],s[1][2],s[1][3],s[1][4],s[1][5] ))
    f1.write('ZZ \t %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n'%(s[2][0],s[2][1],s[2][2],s[2][3],s[2][4],s[2][5] ))
    f1.write('XY \t %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n'%(s[3][0],s[3][1],s[3][2],s[3][3],s[3][4],s[3][5] ))
    f1.write('YZ \t %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n'%(s[4][0],s[4][1],s[4][2],s[4][3],s[4][4],s[4][5] ))
    f1.write('ZX \t %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n'%(s[5][0],s[5][1],s[5][2],s[5][3],s[5][4],s[5][5] ))
    f1.write('\n')
    f1.write('Voigt Approximation:\n')
    f1.write(' Bulk Modulus (kbar):    %4.2f\n'%B)
    f1.write(' Shear Modulus (kbar):   %4.2f\n'%G)
    f1.write(' Long Velocity (cm/s):   %4.2f\n'%vl)
    f1.write(' Shear Velocity (cm/s):  %4.2f\n'%vs)
    f1.write(' Avg. Velocity (cm/s):   %4.2f\n'%va)
    f1.write('\n')
    f1.write('Reuss Approximation:\n')
    f1.write(' Bulk Modulus (kbar):    %4.2f\n'%BR)
    f1.write(' Shear Modulus (kbar):   %4.2f\n'%GR)
    f1.write('\n')
    f1.write("Young's Modulus:\n")
    f1.write(' Y_xx (kbar):           %4.2f\n' %YXX)
    f1.write(' Y_yy (kbar):           %4.2f\n' %YYY)
    f1.write(' Y_zz (kbar):           %4.2f\n' %YZZ)
    f1.write('\n')
    f1.write("Poisson's Ratio (from elastic tensor)\n")
    f1.write(' nu_xy:                 %0.2f\n' %nuxy)
    f1.write(' nu_yz:                 %0.2f\n' %nuyz)
    f1.write(' nu_xz:                 %0.2f\n' %nuxz)
    f1.write('\n')
    f1.write('Additional information:\n')
    f1.write(' Density (g/cm**3):      %4.4f\n'%density)
    f1.write('\n')
    f1.close()
    
    return vl, vs, va

def Debye_temperature(va,data):
    """
    Author: Nicholas Pike
    Email : Nicholas.pike@smn.uio.no
    
    Purpose: Calculate the Debye Temperature
    """
    #declare variables
    thetaD        = 0
    natom         = 0
    density       = 0
    masstot       = 0
    molarmass     = 0
    atom_molecule = 0

    #calculate density
    for i in range(len(data[2])):
        masstot +=data[2][i]*data[4][i] #data[2] is the mass data[4] is the multiplicity
    density = masstot*amu_to_kg/data[3]/kgm3_to_cm3 #data[3] is the volume in g/cm3
    
    #number of unit cells 
    uc_count = 0
    for i in range(len(data[4])):
        uc_count = find_gcd(uc_count,data[4][i])
    
    #calculate molar mass
    for i in range(len(data[2])):
        molarmass += data[2][i]
       
    #calculate number of atoms 
    for i in range(len(data[4])):
        natom +=data[4][i]   
    
    #calculate number of atoms per molecule
    atom_molecule = natom/uc_count

    #Calculate Debye temperature
    thetaD = h/kb*((3.0*atom_molecule*Na*density)/(4.0*np.pi*molarmass))**(1.0/3.0)*va*natom**(-1.0/3.0)*m_to_cm
    
    #print results 
    f1 = open('derived_elastic','a')
    f1.write(' Debye Temp (K):         %4.2f\n'%thetaD)
    f1.write('\n')
    f1.close()
    
    return thetaD

#Hack to get gcd of an array
def find_gcd(x, y):
    """
    Determines the gcd for an array of values.
    """
    while(y):
        x, y = y, x % y
    return x

def gruneisen(vs,vl):
    """
    Author: Nicholas Pike
    Email : Nicholas.pike@smn.uio.no
    
    Purpose: Calculate the gruneisen parameter using poissons ratio
    """
    #declare variables
    poisson = 0
    grun    = 0
    
    #calculate poissons ratio
    poisson = (1.0-2.0*(vs/vl)**2)/(2.0-2.0*(vs/vl)**2)
    
    #calculate gruneisen parameter
    grun = 3.0/2.0*((1.0+poisson)/(2.0-3.0*poisson))
    
    #print results
    f1 = open('derived_elastic','a')
    f1.write(' Poisson Ratio (from sound) %0.4f\n'%poisson)
    f1.write(' Gruneisen Param:           %4.4f\n'%grun)    
    f1.write('\n')
    f1.close()
    
    return grun

def thermal_cond(data,thetaD,grun):
    """    
    Author: Nicholas Pike
    Email : Nicholas.pike@smn.uio.no
    
    Purpose: Calculate the thermal conductivity using Slack 
             J. Phys. Chem. Solids 34, 321 (1973)
    
    """
    #declare variables
    kappa   = 0.0
    A       = 0
    massav  = 0
    natom   = 0
    masstot = 0
    temp    = 300.0
        
    #calculate A
    A = 2.43E-8/(1.0-(0.514/grun)+(0.228/grun**2))
    
    #calculate the number of atoms
    for i in range(len(data[4])):
        natom +=data[4][i]
    
    #Calculate average atomic mass
    for i in range(len(data[2])):
        masstot += data[2][i]*data[4][i] #data[2] is the mass data[4] is the multiplicity
    massav = masstot/natom
    
    #calculate kappa
    kappa = A*massav*(data[3]/ang_to_meter**3)**(1.0/3.0)*thetaD**3/(grun**2*temp)
    
    #print results
    f1 = open('derived_elastic','a')
    f1.write(' Kappa @300K (W/cmK):    %2.4f\n'%kappa)
    f1.close() 
    
    return kappa

"""
Start main program
"""
#starts main program for windows machines... has no effect for other machine types
if __name__ == '__main__':
    __author__     = 'Nicholas Pike'
    __copyright__  = 'none'
    __credits__    = 'none'
    __license__    = 'none'
    __version__    = '0.0'
    __maintainer__ = 'Nicholas Pike'
    __email__      = 'Nicholas.pike@smn.uio.no'
    __status__     = 'experimental'
    __date__       = 'April 2018'
    
#determine name of input file which will be read
if len(sys.argv) < 2: 
    print('Use python Elastic_debye.py --help to view the help menu and \nto learn how to run the program.')
    sys.exit()
if sys.argv[1].startswith('--'):
    if sys.argv[1] == '--help':
        print('--help\t\t Prints this help menu.\n') 
        print('--usage\t\t Allows the user to determine how to use this program.')
        print('--VASP\t\t Determines the derived properties using vasp output.')
        print('--HT\t\t Returns thee Debye Temperature using vasp output data.')
        print('--ANADDB\t Determines the derived properties using anaddb output.')
        print('\n')
        sys.exit()
    elif sys.argv[1] == '--usage':
        print('--usage\t To use this program use the following in the command line.\n python Elastic_debye.py --KEY')
        sys.exit()
    elif sys.argv[1] == '--VASP':
        print('Calculation of the derived elastic properties will begin after\n'
          ' reading the output file.')
        #find outcar file
        filename = find_file('OUTCAR')
        
        #gather data from OUTCAR file
        data     = gather_from_OUTCAR(filename)
        
        #Calculate sound velocities
        vl,vs,va = sound_velocities(data)
        
        #calculate Debye Temperature
        thetaD = Debye_temperature(va,data)
        
        #calculate Gruneisen parameter
        grun = gruneisen(vs,vl)
        
        #calculate thermal conductivity
        kappa = thermal_cond(data,thetaD,grun)
        
        #print references
        f1 = open('derived_elastic','a')
        f1.write('\nReferences:\n')
        f1.write('Derived Elastic properties from:\n')
        f1.write('Phys. Rev. B 95, 155206 (2017)\n')
        f1.write('Pike et al. Unpublished (2018)\n\n')
        f1.write('Slack Model for Thermal Conductivity from:\n')
        f1.write('Slack J. Phys. Chem. Solids 34, 321 (1973)\n')
        f1.close() 
        
        #calculation complete
        print('Calculation complete.  Thank you for using this program.')
        
    elif sys.argv[1] == '--HT':
        #find outcar file
        filename = find_file('OUTCAR')
        
        #gather data from OUTCAR file
        data     = gather_from_OUTCAR(filename)
        
        #Calculate sound velocities
        vl,vs,va = sound_velocities(data)
        
        #calculate Debye Temperature
        thetaD = Debye_temperature(va,data)
        print(thetaD) #prints debye temperature for high-throughput calculations
        
        #calculate Gruneisen parameter
        grun = gruneisen(vs,vl)
        
        #calculate thermal conductivity
        kappa = thermal_cond(data,thetaD,grun)
        
        #print references
        f1 = open('derived_elastic','a')
        f1.write('\nReferences:\n')
        f1.write('Derived Elastic properties from:\n')
        f1.write('Phys. Rev. B 95, 155206 (2017)\n')
        f1.write('Pike et al. Unpublished (2018)\n\n')
        f1.write('Slack Model for Thermal Conductivity from:\n')
        f1.write('Slack J. Phys. Chem. Solids 34, 321 (1973)\n')
        f1.close() 
        
    elif sys.argv[1] == '--ANADDB':
        print('Calculation of the derived elastic properties will begin after\n'
          ' reading the output file.')
        #find anaddb output file
        filename = find_file('ANADDB')
        
        #gather data from OUTCAR file
        data     = gather_from_ANADDB(filename[0],filename[1])
        
        #Calculate sound velocities
        vl,vs,va = sound_velocities(data)
        
        #calculate Debye Temperature
        thetaD = Debye_temperature(va,data)
        
        #calculate Gruneisen parameter
        grun = gruneisen(vs,vl)
        
        #calculate thermal conductivity
        kappa = thermal_cond(data,thetaD,grun)
        
        #print references
        f1 = open('derived_elastic','a')
        f1.write('\nReferences:\n')
        f1.write('Derived Elastic properties from:\n')
        f1.write('Phys. Rev. B 95, 155206 (2017)\n')
        f1.write('Pike et al. Phys. Rev. Mat. 2, 063608 (2018)\n\n')
        f1.write('Slack Model for Thermal Conductivity from:\n')
        f1.write('Slack J. Phys. Chem. Solids 34, 321 (1973)\n')
        f1.close() 
        
        #calculation complete
        print('Calculation complete. Thank you for using this program.')
    else:
        print('ERROR: This tag is unknown.')
else:
    print('Please run program with the --help tag')
    sys.exit()
    
    
    
