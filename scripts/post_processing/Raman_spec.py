# -*- coding: utf-8 -*-
"""
Created: Nov 12, 2016
Author : Nicholas Pike
Email  : Nicholas.pike@ulg.ac.be

Purpose: To calculate the Raman spectrum, at a user defined orientation and angle, 
         by reading in data from an abinit calculation.  This program will read
         the anaddb output file "anaddb.out" or from a user specified file, 
         extract the needed data, and output the result for plotting with your 
         favorite program. 
         
To run:  To run this script simply execute the following code
         python Raman_spec.py "name of input file"
         
Version: 0.0 - Initial building of program
         0.1 - Additional for angle dependent calculation
         0.2 - Additional user input required and number of output files reduced
         0.3 - Moved all user prompted input to an input file

Output: Program will output a text file containing the raman spectrum vs frequency
        and an outfile which outlines what happens in the calculation.

Bugs:   If you find a bug with this program or wish to see a feature added to 
        this script please let me know at Nicholas.pike@ulg.ac.be
        
        Happy Hunting!
"""
#begin main program

"""
Start of definations and other useful information
"""

def GET_UNIT(string):
    if string == 'Ha':
        return 0
    if string == 'Hz':
        return 1
    return 2

def READ_INPUT(user_filein):
    """
    Author: Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
    
    Purpose: Reads the input file and determines if the user's input is acceptable.
    
    Output:  array of input information used by program
    """
    #declare array of  values
    vararray = [0,0,[0,2],[0,2],0,0,False,2,False,1000,[-1.0,2],[-1.0,2]]
    
    #check if file exists
    check = CHECK_FILE(user_filein)  
    
    if check == False:
        print('')
        print('The input file was not found in the directory. \n Please correct this.')
        print('\n Remember the input file should be formated as follows:\n\n '\
                 'filename "name of file"\n outname "name of outfile"\n temp "temperature in Kelvin"\n'
                 ' frequency "frequency in cm^-1"\n spread "spread of lorentz in cm^-1"\n '\
                 'calctype "type of calculation 0- abort, 1- powder, 2-ij polarization, 3- angle"\n')
        sys.exit()
    else:
        for line in open(user_filein):
            li=line.strip('\n')        #Removes any newline command
            if not li.startswith("#"): #Checks if the line does not start with a # 
                                       #character
                l = li.split()
                if len(l)>0:
                  if l[0] == 'filename':          # name of the anaddb output file
                      vararray[0] = str(l[1])
                  elif l[0] == 'temp':            # temperature
                      vararray[1] = float(l[1])
                  elif l[0] == 'laser_freq':       # laser frequency
                      vararray[2][0] = float(l[1])
                      if len(l)>2:
                          vararray[2][1] = GET_UNIT(l[2])
                  elif l[0] == 'spread':          # spread of lorentz
                      vararray[3][0] = float(l[1])
                      if len(l)>2:
                          vararray[3][1] = GET_UNIT(l[2])
                  elif l[0] == 'calctype':        # calculation type
                      vararray[4] = int(l[1])
                  elif l[0] == 'outname':         # output file name
                      vararray[5] = str(l[1])
                  elif l[0] == 'relative_intensity':        # relative intensities or not
                      vararray[6] = True
                  elif l[0] == 'freq_unit':       # unit for output frequencies
                      vararray[7] = GET_UNIT(l[1])
                  elif l[0] == 'keep_file':
                      vararray[8] = True
                  elif l[0] == 'n_freq':
                      vararray[9] = int(l[1])
                  elif l[0] == 'min_freq':
                      vararray[10][0] = float(l[1])
                      if len(l)>2:
                          vararray[10][1] = GET_UNIT(l[2])
                  elif l[0] == 'max_freq':
                      vararray[11][0] = float(l[1])
                      if len(l)>2:
                          vararray[11][1] = GET_UNIT(l[2])
        #set output file name
        global outname
        #check for output file
        outname = vararray[5]
        delete  = vararray[8]
        outname = CHECK_REPEAT(outname,delete)  

        vararray[5] = outname
        printout('Input file located and read in.\n')
                    
    #Now check that the user put a valid name in for the anaddb output file
    # Sanity check...
    if vararray[0] == '' : # no file name given
        print('WARNING: No input file name given. Please try again.')
        sys.exit()
    elif vararray[0].endswith('.nc'):
        print('WARNING: This program uses the output file and not a NetCDF file.')
        sys.exit()
    else:
        check = CHECK_FILE(vararray[0])
    
    if not check:
        printout('')
        printout('The anaddb output file was not found in the directory.\n Please correct this.')
        sys.exit()
    else:
        printout('The anaddb output file was found in the directory.')


    return vararray

def CHECK_FILE(filename): 
    """
    Author: Nicholas Pike
    Email : Nicholas.pike@ulg.ac.be
        
    Purpose: Checks if the file "filename" exists in the current directory and outputs
    a boolean value indicating whether or not the file is found.
    """
    if os.path.isfile(filename):
        logic  = True
    else:
        logic  = False
        
    return logic
    
def PRINT_HEADER():
    """
    Author: Nicholas Pike
    Email : Nicholas.pike@ulg.ac.be
        
    Purpose: Prints a header to the output file and terminal when the user
    starts to run the program. 
    """
       
    #Header for output file
    printout('')
    printout('++++++++++++++++++++++++ Version %s ++++++++++++++++++++++++\n'%__version__)
    printout('')
    printout('This program generates a Raman Spectrum after an anaddb '\
             'calculation.\n')
    printout('In this version, the program reads an the output file of an '\
             'anaddb\n calculation, finds the information it needs, and outputs '\
             'that\n information to a file.')
    printout('')
    printout('Author: Nicholas Pike')
    printout('Email:  Nicholas.pike@ulg.ac.be')
    printout('')
           
    return 
    
def printout(to_output):
    """
    Author: Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
        
    Purpose: This defintion should print to an output file, known as the output,
    in which all comments, warnings, and results are printed too. 
    """
    #print to file
    f1= open(outname,'a')
    f1.write(to_output+'\n')
    f1.close()
        
    return
    
def printoutfile(output,outfile):
    """
    Author: Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
    
    Purpose: Print raman spectrum as a function of frequency to the output file.
    """
    f1= open(outfile,'a')
    f1.write(output+'\n')
    f1.close()
    return
    
def CHECK_REPEAT(filename,delete):
    """
    Author: Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
        
    In order to create a new output file a file the name must not be used. If 
    this file already exists then we rename the file until it no longer exists.
    """
    c = 1
    out = CHECK_FILE(filename)
    newfilename = filename #If the file does not exist
    if out and delete:
        os.remove(filename)
    else:
        while out:
            if out == True:
                string = filename+ str(c)
                c += 1
                newfilename = string
                out = CHECK_FILE(newfilename)
            else:
                newfilename = filename
                break
        
    return newfilename

def UNIT_TO_HA(val,unit):
    if not (unit in [0,1,2]):
        return -1
    if unit==0:# val in Ha
        return val
    if unit==1:# val in Hz
        return val*Hz_to_Ha
    if unit==2:# val in cm1
        return val*cm1_to_hartree

def HA_TO_UNIT(val,unit):
    if not (unit in [0,1,2]):
        return -1
    if unit==0:# val in Ha
        return val
    if unit==1:# val in Hz
        return val/Hz_to_Ha
    if unit==2:# val in cm1
        return val/cm1_to_hartree

def CALL_RAMAN_MENU(output,keywords,vararray):
    """
    Author:Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
    
    Purpose: To call and activate the raman spectrum part of this calculation.
    """
    T          = vararray[1]
    [laser,laser_unit] = vararray[2]
    [width,width_unit] = vararray[3]
    relative   = vararray[6]
    freq_unit  = vararray[7]
    delete     = vararray[8]
    freqsteps  = vararray[9]
    [min_freq,min_unit] = vararray[10]
    [max_freq,max_unit] = vararray[11]

    # Compute spectra or not
    with_spreading = True
    if width<=0:
        with_spreading = False
    else:
        width = UNIT_TO_HA(width,width_unit)

    #Convert frequencies in Ha
    laser = UNIT_TO_HA(laser,laser_unit)

    #Convert min/max frequencies in freq_unit (after a conversion in Ha)
    min_freq = UNIT_TO_HA(min_freq,min_unit)
    min_freq = HA_TO_UNIT(min_freq,freq_unit)
    max_freq = UNIT_TO_HA(max_freq,max_unit)
    max_freq = HA_TO_UNIT(max_freq,freq_unit)
    #Get some data
    modeenergy = output[1]
    modedata   = output[2] 
    qdirs      = output[4]
    n_modes    = len(modeenergy[0])

    #remainder of printing goes to log file only
    printout('')
    printout('Anaddb created with Abinit Version: %s' %keywords[0])
    printout('')

    # Print header of output file
    if not with_spreading:
        outfile_base = outname+'_intensity'
    else:
        outfile_base = outname+'_spec'

    # Start loop on q directions
    for idir,qdir in enumerate(qdirs):

        # LO modes or not?
        LO = False
        if np.sqrt(qdir[0]**2+qdir[1]**2+qdir[2]**2)>10**(-10):
            LO = True

        printout('************************************************************************')
        if LO:
          printout('*** LO modes, in the q direction (cartesian coordinates) : %s %s %s '%(qdir[0],qdir[1],qdir[2]))
          outfile = outfile_base + '_%04.0f_%04.0f_%04.0f'%(qdir[0]*1000,qdir[1]*1000,qdir[2]*1000)
        else:
          printout('*** TO modes')
          outfile = outfile_base
        printout('************************************************************************')
        printout('')

        #Remove existing file
        out = CHECK_FILE(outfile)
        if out and delete:
            os.remove(outfile)
       
        #Modify the input data silently.
        menergy = GET_MODEENERGY(modeenergy[idir])
        rarray  = GET_RAMANMATRIX(modedata[idir])

        # Compute all intensities
        printout('Calculating the raman intensities.')

        I_powder,I_powder_rel = RAMAN_INTENSITIES(menergy,rarray,laser,T,'POWDER')
        I_xx,I_xx_rel         = RAMAN_INTENSITIES(menergy,rarray,laser,T,'XX')
        I_xy,I_xy_rel         = RAMAN_INTENSITIES(menergy,rarray,laser,T,'XY')
        I_xz,I_xz_rel         = RAMAN_INTENSITIES(menergy,rarray,laser,T,'XZ')
        I_yy,I_yy_rel         = RAMAN_INTENSITIES(menergy,rarray,laser,T,'YY')
        I_yz,I_yz_rel         = RAMAN_INTENSITIES(menergy,rarray,laser,T,'YZ')
        I_zz,I_zz_rel         = RAMAN_INTENSITIES(menergy,rarray,laser,T,'ZZ')

        # Print header in the output file for this q direction
        if LO:
            printoutfile('#*** LO modes : %s %s %s '%(qdir[0],qdir[1],qdir[2]),outfile)
        else:
            printoutfile('#*** TO modes',outfile)
        
        string_freq = '  freq (Ha)'
        if freq_unit == 1:
            string_freq = '  freq (Hz)'
        if freq_unit == 2:
            string_freq = ' freq (cm1)'

        if not relative:
            printoutfile('#'+string_freq+' abs Ipowder   abs I(XX)   abs I(XY)   abs I(XZ)   abs I(YY)   abs I(YZ)   abs I(ZZ)',outfile)
        else:
            printoutfile('#'+string_freq+' rel Ipowder   rel I(XX)   rel I(XY)   rel I(XZ)   rel I(YY)   rel I(YZ)   rel I(ZZ)',outfile)

        if with_spreading:
            minf = np.inf # infinity
            maxf = 0.0
            # If required by the user : set the minimum and/or maximum values
            if min_freq >= 0.0:
                minf = min_freq
            if max_freq >= 0.0:
                maxf = max_freq

        max_rarray = max(abs(np.ravel(rarray)))
        # Loop on modes
        for j in range(n_modes):
            rmatj = rarray[j]
            max_rmatj = max(abs(np.ravel(rmatj)))
            if max_rmatj < 10**(-3) * max_rarray:
                printout(' -- Mode %i: %10.2f (cm-1) : Raman tensor is negligible'%(j+1,menergy[j][2]))
            if max_rmatj < 10**(-3) * max_rarray:
                printout(' -- Mode %i: %10.2f (cm-1) : Raman tensor is negligible'%(j+1,menergy[j][2]))
            else:
                printout(' -- Mode %i: %10.2f (cm-1)'%(j+1,menergy[j][2]))
                printout('    -- Raman tensor :')
                printout('    %17.9f %17.9f %17.9f' %(rmatj[0][0],rmatj[0][1],rmatj[0][2]))
                printout('    %17.9f %17.9f %17.9f' %(rmatj[1][0],rmatj[1][1],rmatj[1][2]))
                printout('    %17.9f %17.9f %17.9f' %(rmatj[2][0],rmatj[2][1],rmatj[2][2]))
                printout('    -- Raman intensities :')
                printout('        | powder     | (XX)      | (XY)      | (XZ)      | (YY)      | (YZ)      | (ZZ)')
                l_I = [I_powder[j],I_xx[j],I_xy[j],I_xz[j],I_yy[j],I_yz[j],I_zz[j]]
                l_I_rel = [I_powder_rel[j],I_xx_rel[j],I_xy_rel[j],I_xz_rel[j],I_yy_rel[j],I_yz_rel[j],I_zz_rel[j]]
                printout('    abs | %10.5e|%10.5e|%10.5e|%10.5e|%10.5e|%10.5e|%10.5e'%(l_I[0],
                         l_I[1],l_I[2],l_I[3],l_I[4],l_I[5],l_I[6]))
                printout('    rel | %11.2f|%11.5f|%11.5f|%11.5f|%11.5f|%11.5f|%11.5f'%(l_I_rel[0],
                         l_I_rel[1],l_I_rel[2],l_I_rel[3],l_I_rel[4],l_I_rel[5],l_I_rel[6]))
                printout('')

                # if width <=0 : write information to the output file with a simple format
                if not with_spreading:
#                    printout('Printing results to an output file named %s' %outfile)
                    if relative:
                        printoutfile(' %11.2f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f' %(menergy[j][freq_unit],
                    l_I_rel[0],l_I_rel[1],l_I_rel[2],l_I_rel[3],l_I_rel[4],l_I_rel[5],l_I_rel[6]),outfile)
                    else:
                        printoutfile(' %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e' %(menergy[j][freq_unit],
                        l_I[0],l_I[1],l_I[2],l_I[3],l_I[4],l_I[5],l_I[6]),outfile)
                # For accumulation of intensities
                else:
                    freq = menergy[j][freq_unit]
                    if min_freq < 0.0 and minf >= freq:
                        minf = freq*0.95
                    if max_freq < 0.0 and maxf <= freq:
                        maxf = freq*1.05

        #Compute Raman Spectra
        if with_spreading:
            maxval = np.zeros(7)
            ramanplot = np.zeros([8,freqsteps])
            # For each freqstep : accumulate all modes and store maximum values
            for i in range(freqsteps):
                freq = minf+i*(maxf-minf)/freqsteps
                freq_Ha = UNIT_TO_HA(freq,freq_unit)
                intstep = np.zeros(7)
                for j in range(n_modes):
                    rmatj = rarray[j]
                    norm_rmatj = sum(abs(np.ravel(rmatj)))
                    if norm_rmatj > 10**(-10):
                        l_I = [I_powder[j],I_xx[j],I_xy[j],I_xz[j],I_yy[j],I_yz[j],I_zz[j]]
                        spread = width/((freq_Ha-menergy[j][0])**2+width**2)
                        for iint in range(7):
                            intstep[iint] += spread*l_I[iint]
                for iint in range(7):
                    if maxval[iint] < intstep[iint]:
                        maxval[iint] = intstep[iint]
                    ramanplot[iint][i] = intstep[iint]
            # For each freqstep : write the output intensity
            for i in range(freqsteps):
                freq = minf+i*(maxf-minf)/freqsteps # Frequency in cm1
                if relative:
                    for iint in range(7):
                        if maxval[iint] > 0.0:
                            ramanplot[iint][i] = ramanplot[iint][i] / maxval[iint]
                printoutfile('%10.5e %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e' %(freq,
                    ramanplot[0][i],ramanplot[1][i],ramanplot[2][i],ramanplot[3][i],
                    ramanplot[4][i],ramanplot[5][i],ramanplot[6][i]),outfile)

    return
        
def LOAD_ANADDB(infile):
    """
    Author: Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
        
    Loads abinit data by reading the output file and storing necessary information
    in the correct arrays.  This is done by reading the file twice. First to 
    find keywords, and the second time to extract dat
    
    """
    keywords =['', #Abinit version
               False, #dieflag - frequency dependent dielectric constant
               False, #nfreq   - number of frequency steps
               False, #nlflag  - Raman tensor and nonlinear optical tensor
               False, #elaflag - elastic tensor flag
               False, #piezoflag - Piezoelectric tensor
               ]
    #First, read the output file and determine if any of the keywords are present
    with open(infile,'r') as f:
        for num,line in enumerate(f,1):
            if line.startswith( '.Version'):
                l = line.strip('\n').split()
                keywords[0] = str(l[1])
            elif line.startswith( '     dieflag'):
                l = line.strip('\n').split()
                dieflag = int(l[len(l)-1])
                if dieflag > 0 and dieflag != 2:
                    keywords[1] = True
                else:
                    keywords[1] = False
            elif line.startswith( '       nfreq'):
                l = line.strip('\n').split()
                nfreq = int(l[len(l)-1])
                if nfreq > 10:
                    keywords[2] = nfreq
                else:
                    keywords[2] = 0
                    keywords[1] = False
            elif line.startswith( '      nlflag'):
                l = line.strip('\n').split()
                if int(l[len(l)-1])> 0:
                    keywords[3] = True
                else: 
                    keywords[3] = False
    
                #Now that keywords are read in and stored. We can read the file again and
                #figure out what information needs to be stored
    
    
    #declare storage arrays    
    outinfo= [0,0,0,0,0]
    qdirs = []
    modeenergy = []
    modedata = []
    diedata = []
    
    if keywords[1]: #dielectric tensor as a function of frequency
        #Do something
        printout('Starting extraction of the dielectric Tensor as a function of frequency.')
        diedata=np.zeros(shape=(keywords[2],7))
        with open(infile,'r') as f:
            for num,line in enumerate(f,1):
                if line.startswith(' Frequency(Hartree)    Dielectric constant                Reflectivity'):
                    c=2
                    for c in range(2,keywords[2]+2):
                        data = linecache.getline(infile,num+c)
                        l=data.strip('\n').split()                            
                        diedata[c-2][0]= l[0]
                        diedata[c-2][1]= l[1]
                        diedata[c-2][2]= l[2]
                        diedata[c-2][3]= l[3]
                        diedata[c-2][4]= l[4]
                        diedata[c-2][5]= l[5]
                        diedata[c-2][6]= l[6]
                            
        outinfo[3] = diedata
        printout('Finished extraction of the dielectric Tensor as a function of frequency.')

    if keywords[3]: #Raman tensor and modes
        printout('Starting extraction of the raman tensor.')
        with open(infile,'r') as f:
            iqdir = 0
            for num,line in enumerate(f,1):
                if ' Raman susceptibilit' in line:
                    modeenergy.append([])
                    modedata.append([])
                    iqdir += 1
                    if 'transverse' in line : # TO mode
                        qdir = np.array([0.0,0.0,0.0])
                    if 'non-analyticity' in line : # LO mode
                        data = linecache.getline(infile,num+1)
                        l = data.strip('\n').split()
                        qdir = np.zeros(3)
                        for ii in range(3):
                          qdir[ii] = float(l[ii+3])
                    qdirs.append(qdir)
                if line.startswith(' Mod') and 'cm-1)' in line:
                    l = line.strip('\n').split('(')
                    onemode_energy = float(l[1].split()[0])
                    onemode_raman = np.zeros([3,3])
                    for ii in range(3):  
                        data = linecache.getline(infile,num+ii+1)
                        l = data.strip('\n').split()
                        for jj in range(3):
                          onemode_raman[ii][jj] =l[jj+1]
                    modeenergy[iqdir-1].append(onemode_energy)
                    modedata[iqdir-1].append(onemode_raman)

        outinfo[1] = modeenergy
        outinfo[2] = modedata
        outinfo[4] = qdirs

        printout('Finished extraction of the raman tensor.')
        printout('')
            
    return keywords,outinfo
    
def GET_MODEENERGY(modeenergy):
    """
    Author: Nicholas Pike
    Email: Nicholas.Pike@ulg.ac.be
    
    Purpose: Convert the Mode energies from cm-1 to Ha for later calculations
    """
    nmode   = int(len(modeenergy))
    energy = np.zeros(shape=(nmode,3))
    
    for j in range(nmode):
        energy[j][0] = float(modeenergy[j])*cm1_to_hartree 
        energy[j][1] = float(modeenergy[j])*cm1_to_hz
        energy[j][2] = float(modeenergy[j])
        
    return energy

def GET_RAMANMATRIX(modedata):
    """
    Author: Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
    
    Purpose: Convert the long list of matrix elements to 3x3 arrays for each 
    mode index
    """
    nmode = len(modedata)
    for imode in range(nmode):
        for ii in range(3):
            for jj in range(3):
                value = modedata[imode][ii][jj]
                if abs(value)<10**(-10):
                  modedata[imode][ii][jj] = 0.0
    return modedata

def GET_NORM(m):
  norm = 0
  for ii in range(3):
    for jj in range(3):
      norm += m[ii][jj]**2
  return norm

def GET_G012(m):
    m_I = np.zeros([3,3],dtype=np.float64)
    m_A = np.zeros([3,3],dtype=np.float64)
    m_S = np.zeros([3,3],dtype=np.float64)
    trace_m = m[0][0] + m[1][1] + m[2][2]
    for ii in range(3):
      m_I[ii][ii] =  1.0/3.0*trace_m
      m_S[ii][ii] = -1.0/3.0*trace_m
      for jj in range(3):
        m_A[ii][jj]  = (m[ii][jj] - m[jj][ii])/2.0
        m_S[ii][jj] += (m[ii][jj] + m[jj][ii])/2.0
    m_sum = m_I + m_S + m_A
    if sum(abs(np.ravel(m_sum-m))) > 1e-14:
      print("ERROR in GET_G012 : m_sum is not m!")
      sys.exit()
    G0 = GET_NORM(m_I)
    G1 = GET_NORM(m_A)
    G2 = GET_NORM(m_S)
#    G0 = 1.0/3.0*(m[0][0]+m[1][1]+m[2][2])**2
#    G1 = 1.0/2.0*((m[0][1]-m[1][0])**2+(m[0][2]-m[2][0])**2+(m[1][2]-m[2][1])**2)
#    G2 = 1.0/2.0*((m[0][1]+m[1][0])**2+(m[0][2]+m[2][0])**2+(m[1][2]+m[2][1])**2)
#    G2+= 1.0/3.0*((m[0][0]-m[1][1])**2+(m[1][1]-m[2][2])**2+(m[2][2]-m[0][0])**2)
    return G0,G1,G2

def RAMAN_INTENSITIES(menergy,rarray,laser,T,option):
    """
    Author: Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
    
    Purpose: To compute the intensities of each modes
    """
    nmode = len(menergy)
    Intensities     = np.zeros(nmode)
    Intensities_rel = np.zeros(nmode)
#    if option == 'POWDER':
#      printout(" hbar/kT = {:12.5e}".format(hplank/(kb*T)))
    for j in range(nmode):
        if menergy[j][2] > 1: # in cm1
            mbose = 1.0/(np.exp(hplank*menergy[j][1]/(kb*T))-1.0)
            pref = 1.0/(2.0*menergy[j][0])*(menergy[j][0] - laser)**4/clight**4*(mbose+1.0)
            if option == 'POWDER':
#                printout("mode {:d} : {:12.5e} mbose = {:8.5f}".format(j+1,menergy[j][1],mbose+1.0))
#                for ii in range(3):
#                  printout(" {:12.5e} {:12.5e} {:12.5e} ".format(rarray[j][ii][0],rarray[j][ii][1],rarray[j][ii][2]))
                G0,G1,G2 = GET_G012(rarray[j])
#                printout("G1 = {:12.5e} G1 = {:12.5e} G2 = {:12.5e}".format(G0,G1,G2))
                pref  = (2*np.pi)*pref
                Ipar  = pref*(10.0*G0+4.0*G2)
                Iperp = pref*( 5.0*G1+3.0*G2)
                Intensities[j] = Ipar + Iperp
            else:
                Inpol = []
                Outpol = []
                if  option == 'XX':
                    Inpol    = np.array([1,0,0])
                    Outpol   = np.array([1,0,0])
                elif option == 'XY':
                    Inpol    = np.array([1,0,0])
                    Outpol   = np.array([0,1,0]) 
                elif option == 'XZ':
                    Inpol    = np.array([1,0,0])
                    Outpol   = np.array([0,0,1]) 
                elif option == 'YY':
                    Inpol    = np.array([0,1,0])
                    Outpol   = np.array([0,1,0])
                elif option == 'YZ':
                    Inpol    = np.array([0,1,0])
                    Outpol   = np.array([0,0,1]) 
                elif option == 'ZZ':
                    Inpol    = np.array([0,0,1])
                    Outpol   = np.array([0,0,1])
                rj1 = [rarray[j][0][0],rarray[j][0][1],rarray[j][0][2]]
                rj2 = [rarray[j][1][0],rarray[j][1][1],rarray[j][1][2]]
                rj3 = [rarray[j][2][0],rarray[j][2][1],rarray[j][2][2]]
                modearray = np.array([rj1,rj2,rj3])
                alpha_ij = np.dot(Outpol,np.dot(modearray,Inpol))
                Intensities[j] = alpha_ij**2*pref
      
    max_int = max(Intensities)
    if max_int > 0.0:
        Intensities_rel = Intensities/max_int
                        
    return Intensities,Intensities_rel
    
def DETER_MENU(var_array):
    """
    Author: Nicholas Pike
    Email : Nicholas.pike@ulg.ac.be
    
    Determine what calculation to ask the user about depending on what anaddb
    calculation was done.  This is done by reading the file and looking for 
    specific flags in the output file
    
    Args: infile - name of input file (assumed to be anaddb.out)
    """
    infile = var_array[0]
    
    keywords,output = LOAD_ANADDB(infile)
           
    if keywords[1]:
        #run module to print Dielectric tensor as a function of frequency
        outfile = outname+'_dielectric_freq'
        printoutfile('# Freq    dielectric      Reflectivity',outfile)
        printoutfile('#           x y z x y z  ',outfile)
        for j in range(int(keywords[2])):
            o3j=output[3][j]
            printoutfile('%e %e %e %e %e %e %e' %(o3j[0],o3j[1],o3j[2],o3j[3],o3j[4],o3j[5],o3j[6]),outfile)
    if keywords[3]:
        #Choose to identify the modes characteristics or make a Raman Spectrum  
        printout('Entering Raman Spectrum Calculation.')
        printout('')
        CALL_RAMAN_MENU(output,keywords,var_array) #Calls the menu for the Raman calculation
              
    printout('Thank you for using this program.')
    
    return 
    
    
"""
Execution Code below this line
"""
   
#Program runs via this if statement (needed for windows computers):
if __name__ == '__main__':
    #declare important information
    __author__     = 'Nicholas Pike'
    __copyright__  = 'none'
    __credits__    = 'none'
    __license__    = 'none'
    __version__    = '0.3'
    __maintainer__ = 'Nicholas Pike'
    __email__      = 'Nicholas.pike@ulg.ac.be'
    __status__     = 'production'
    __date__       = 'November 2016'
    
    #Import useful python programs
    import os 
    import sys
    import numpy as np
    import linecache
    systemversion = sys.version_info
    
    #Declare global variables
    global outname
    global cm1_to_hartree 
    global cm1_to_hz
    global Hz_to_Ha
    global hplank
    global kb
    global T
    global clight
    global width
    
    cm1_to_hartree = 4.55633E-6         # conversion factor between cm-1 and Hartree
    cm1_to_hz      = 2.99793E10         # conversion factor between cm-1 and Hz
    Hz_to_Ha       = 1.519828500716E-16 # conversion factor between Hz and Hartree 
    kb             = 8.6173324E-5       # Boltzmann constant in eV/K
    T              = 0                  # Default Temperature (user input variable)
    hplank         = 6.58211928E-16     # h in eVs
    clight         = 137.0359997566     # Speed of light in atomic units
    width          = 0.0                # spread of lorentian (user input variable)

    if any('SPYDER' in name for name in os.environ):
        user_inputfile = 'input_raman'
    else:
        #Determines what the program is to do if it is run from command line
        user_inputfile = sys.argv[1] #input should be python program_name tfile 
          
    #Name of default input file
    vararray = READ_INPUT(user_inputfile) #The program will look for the specified input file
    
       
    #Print header file and start the calculation
    PRINT_HEADER()
    
    #Read anaddb.out file and determine what type of calculation to run
    DETER_MENU(vararray)
    
    
#Ends program
