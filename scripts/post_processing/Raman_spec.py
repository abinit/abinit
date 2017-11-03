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
def READ_INPUT(user_filein,outname):
    """
    Author: Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
    
    Purpose: Reads the input file and determines if the user's input is acceptable.
    
    Output:  array of input information used by program
    """
    #declare array of  values
    vararray = [0,0,0,0,0,0]
    
    #check if file exists
    check = CHECK_FILE(user_filein)
    
    #check for output file
    outname = CHECK_REPEAT(outname)  
    
    vararray[5] = outname  #saved here just in case
    
    if check == False:
        printout('')
        printout('The input file was not found in the directory. \n Please correct this.')
        printout('\n Remember the input file should be formated as follows:\n\n '\
                 'file name "name of file"\n temp "temperature in Kelvin"\n frequency '\
                 '"frequency in cm^-1"\n spread "spread of lorentz in cm^-1"\n '\
                 'calctype "type of calculation 0- abort, 1- powder, 2-ij polarization, 3- angle"\n')
        sys.exit()
    else:
        printout('Input file located.  Reading input file.')
        for line in open(user_filein):
            li=line.strip('\n')        #Removes any newline command
            if not li.startswith("#"): #Checks if the line does not start with a # 
                                       #character
                l = li.split(' ')
                if l[0] == 'filename':          # name of the anaddb output file
                    vararray[0] = str(l[1])
                elif l[0] == 'temp':            # temperature
                    vararray[1] = float(l[1])
                elif l[0] == 'frequency':       # laser frequency
                    vararray[2] = float(l[1])
                elif l[0] == 'spread':          # spread of lorentz
                    vararray[3] = float(l[1])
                elif l[0] == '':                # passes over spaces
                    fill = 0
                elif l[0] == 'calctype':        # calculation type
                    vararray[4] = int(l[1])
                    
    #Now check that the user put a valid name in for the anaddb output file
    check = CHECK_FILE(vararray[0])
    if check == False:
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
    
def CHECK_REPEAT(filename):
    """
    Author: Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
        
    In order to create a new output file a file the name must not be used. If 
    this file already exists then we rename the file until it no longer exists.
    """
    c = 1
    out = CHECK_FILE(filename)
    newfilename = filename #If the file does not exist 
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


    
def CALL_RAMAN_MENU(output,keywords,vararray):
    """
    Author:Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
    
    Purpose: To call and activate the raman spectrum part of this calculation.
    """
    T       = vararray[1]
    input2  = vararray[2]
    width   = vararray[3]
    choice1 = vararray[4]
    
    #Decide what to do when the user makes a choice
    ramanplot = 0
    if choice1 ==0:
        printout('Calculation Aborted.')
        sys.exit()

    elif choice1 == 1: 
        printout('Calculating the raman spectrum for a powder sample.')
        
        #remainder of printing goes to log file only
        printout('')
        printout('Anaddb created with Abinit Version: %s' %keywords[0])
        printout('')
        printout('Phonon Mode Energies (cm-1).\n')
        printout('\n')
        for j in range(int(len(output[1])/2)):
            printout('Mode %i: %s\n'%(j+1,output[1][j]))
    
        printout('')
        printout('Raman intensity Matrix (linearlized):')
        printout('')        
        for j in range(int(len(output[1])/2)):
            printout('Mode %i:\n' %(j+1))        
            printout('%s %s %s %s %s %s %s %s %s\n' %(output[2][j][0],output[2][j][1],output[2][j][2],output[2][j][3],output[2][j][4],output[2][j][5],output[2][j][6],output[2][j][7],output[2][j][8]))
        printout('')
        #Modify the input data silently.
        menergy = GET_MODEENERGY(output)
        rarray  = GET_RAMANMATRIX(output)
        ramanplot = RAMAN_POWDER(menergy,rarray,input2,width,T)
        
        #After the calculation completes we print the results to a file
        outfile = 'raman_powder.out'
        printout('Printing results to an output file named %s' %outfile)
        printout('')
        
        printoutfile('#Freq:        Raman Intensity',outfile)
        for j in range(ramanplot.shape[0]):
            printoutfile('%e %e ' %(ramanplot[j][0],ramanplot[j][1]),outfile)
            
        printout('Printing to output file complete.')
        printout('')
                
    elif choice1 == 2:
        printout('Calculating the raman spectrum for ij polarization.')
                       
        printout('Anaddb created with Abinit Version: %s' %keywords[0])
        printout('')
        printout('Phonon Mode Energies (cm-1).')
        printout('')
        for j in range(int(len(output[1])/2)):
            printout('Mode %i: %s\n'%(j+1,output[1][j]))
    
        printout('')
        printout('Raman intensity Matrix (linearlized):')
        printout('')        
        for j in range(int(len(output[1])/2)):
            printout('Mode %i:' %(j+1))        
            printout('%s %s %s %s %s %s %s %s %s' %(output[2][j][0],output[2][j][1],output[2][j][2],output[2][j][3],output[2][j][4],output[2][j][5],output[2][j][6],output[2][j][7],output[2][j][8]))
        printout('')
        #Modify the input data silently.
        menergy = GET_MODEENERGY(output)
        rarray  = GET_RAMANMATRIX(output)
        ramanxx = RAMAN_POLAR(menergy,rarray,input2,'XX',width,T)
        ramanxy = RAMAN_POLAR(menergy,rarray,input2,'XY',width,T)
        ramanxz = RAMAN_POLAR(menergy,rarray,input2,'XZ',width,T)
        ramanyy = RAMAN_POLAR(menergy,rarray,input2,'YY',width,T)
        ramanyz = RAMAN_POLAR(menergy,rarray,input2,'YZ',width,T)
        ramanzz = RAMAN_POLAR(menergy,rarray,input2,'ZZ',width,T)
        
        #After the calculation completes we print the results to a file
        outfile = 'raman_ij.out'
        printout('Printing results to an output file named %s' %outfile)
        printout('')
        
        printoutfile('#Freq:        Raman Intensity ( XX, XY, XZ, YY, YZ, ZZ):',outfile)
        for j in range(ramanxx.shape[0]):
            printoutfile('%e %e %e %e %e %e %e ' %(ramanxx[j][0],ramanxx[j][1],
                                                   ramanxy[j][1],ramanxz[j][1],
                                                   ramanyy[j][1],ramanyz[j][1],
                                                   ramanzz[j][1]),outfile)
            
        printout('Printing to output file complete.')
        printout('')
        
    elif choice1 == 3:
        printout('Calculating the raman spectrum as a function of the angle theta.')
                
        printout('Anaddb created with Abinit Version: %s' %keywords[0])
        printout('')
        printout('Phonon Mode Energies (cm-1).')
        printout('')
        for j in range(int(len(output[1])/2)):
            printout('Mode %i: %s\n'%(j+1,output[1][j]))
    
        printout('')
        printout('Raman intensity Matrix (linearlized):')
        printout('')        
        for j in range(int(len(output[1])/2)):
            printout('Mode %i:\n' %(j+1))        
            printout('%s %s %s %s %s %s %s %s %s\n' %(output[2][j][0],output[2][j][1],output[2][j][2],output[2][j][3],output[2][j][4],output[2][j][5],output[2][j][6],output[2][j][7],output[2][j][8]))
        printout('')
        #Modify the input data silently.
        menergy = GET_MODEENERGY(output)
        rarray  = GET_RAMANMATRIX(output)
        ramanplot = RAMAN_POLAR(menergy,rarray,input2,'THETA',width,T)
        
        #After the calculation completes we print the results to a file
        outfile = 'raman_theta.out'
        printout('Printing results to an output file named %s' %outfile)
        printout('')
        
        printoutfile('#Theta:        Raman Intensity',outfile)
        freq = ''
        for j in range(ramanplot.shape[0]):
            freq = '%e '%ramanplot[j][0]
            string = ''
            for k in range(1,ramanplot.shape[1]):
                string += '%e ' %ramanplot[j][k]                
            printoutfile(freq + string,outfile )
            
        printout('Printing to output file complete.')
        printout('')
                        
    elif choice1 >8:
        printout('Your choice of %s, is not accepted.\n Please try again.')
        sys.exit()
        
    return 
        
def LOAD_ANADDB(infile):
    """
    Author: Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
        
    Loads abinit data by reading the output file and storing necessary information
    in the correct arrays.  This is done by reading the file twice. First to 
    find keywords, and the second time to extract dat
    
    """
    keywords = np.array(['', #Abinit version
                         False, #dieflag - frequency dependent dielectric constant
                         False, #nfreq   - number of frequency steps
                         False, #nlflag  - Raman tensor and nonlinear optical tensor
                         False, #elaflag - elastic tensor flag
                         False, #piezoflag - Piezoelectric tensor
                         ])
    #First, read the output file and determine if any of the keywords are present
    with open(infile,'r') as f:
        for num,line in enumerate(f,1):
            if line.startswith( '.Version'):
                l = line.strip('\n').split(' ')
                keywords[0] = str(l[1])
            elif line.startswith( '     dieflag'):
                l = line.strip('\n').split(' ')
                if int(l[len(l)-1]) > 0:
                    keywords[1] = True 
                else:
                    keywords[1] = False
            elif line.startswith( '       nfreq'):
                l = line.strip('\n').split(' ')
                if int(l[len(l)-1]) > 10:
                    keywords[2] = int(l[len(l)-1])
                else:
                    keywords[2] = 0
                    keywords[1] = False
            elif line.startswith( '      nlflag'):
                l = line.strip('\n').split(' ')
                if int(l[len(l)-1])> 0:
                    keywords[3] = True
                else: 
                    keywords[3] = False
            else:
                keywords[0]=keywords[0]
    
                #Now that keywords are read in and stored. We can read the file again and
                #figure out what information needs to be stored
    
    
    #declare storage arrays    
    repeat = 0
    outinfo= [0,0,0,0]
    modeenergy = []
    modedata = []
    diedata = []
    
    #Now look through keywords and extract desired information
    for i in range(len(keywords)): 
        if keywords[i] and i == 1 and keywords[i] != 'False': #dielectric tensor as a function of frequency
            #Do something
            printout('Starting extraction of the dielectric Tensor as a function of frequency.')
            diedata=np.zeros(shape=(int(keywords[2]),7))
            with open(infile,'r') as f:
                for num,line in enumerate(f,1):
                    if line.startswith(' Frequency(Hartree)    Dielectric constant                Reflectivity'):
                        c=2
                        for c in range(2,int(keywords[2])+2):
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

        elif keywords[i] and i == 3: #Raman tensor and modes
            printout('Starting extraction of the raman tensor.')
            with open(infile,'r') as f:
                for num,line in enumerate(f,1):
                    if line.startswith(' Phonon frequencies in cm-1    :'):
                        c= 1 #starts couting at 1
                        readstop = True
                        while readstop:
                            data = linecache.getline(infile,num+c)
                            l = data.strip('\n').split()
                            for j in range(1,len(l)):
                                modeenergy = np.append(modeenergy,l[j])
                            if len(l) <6:
                                readstop = False
                            else:
                                c+=1  
            outinfo[1] = modeenergy
                            
            modedata=np.zeros((int(len(modeenergy)/2),9))
            with open(infile,'r') as g:
                d = 1
                for num,line in enumerate(g,1):
                    c = 1 #starts couting at 1
                    if line.startswith(' Mod ') and repeat <=len(modeenergy)/2-1 :
                        readstop = True
                        while readstop:
                            data = linecache.getline(infile,num+c)
                            l = data.strip('\n').split()
                            modedata[d-1][3*(c-1)+0] =l[1]
                            modedata[d-1][3*(c-1)+1] =l[2]
                            modedata[d-1][3*(c-1)+2] =l[3]
                            if c >=3:
                                d+=1
                                readstop = False
                            else:
                                c+=1
                        repeat += 1                            
                
            outinfo[2] = modedata
            printout('Finished extraction of the raman tensor.')
            printout('')
            
    return keywords,outinfo
    
def GET_MODEENERGY(outinfo):
    """
    Author: Nicholas Pike
    Email: Nicholas.Pike@ulg.ac.be
    
    Purpose: Convert the Mode energies from cm-1 to Ha for later calculations
    """
    menergy = outinfo[1]
    energy = np.zeros(shape = (int(len(menergy)/2),2))
    
    for j in range(int(len(menergy)/2)):
        energy[j][0] = float(menergy[j])*cm1_to_hartree 
        energy[j][1] = float(menergy[j])*cm1_to_hz
        
    return energy

def GET_RAMANMATRIX(outinfo):
    """
    Author: Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
    
    Purpose: Convert the long list of matrix elements to 3x3 arrays for each 
    mode index
    """
    array = np.zeros(shape= (int(len(outinfo[1])/2),3,3))
    for x in range(int(len(outinfo[1])/2)):
        array[x][0][0] = outinfo[2][x][0]
        array[x][0][1] = outinfo[2][x][1]
        array[x][0][2] = outinfo[2][x][2]
        array[x][1][0] = outinfo[2][x][3]
        array[x][1][1] = outinfo[2][x][4]
        array[x][1][2] = outinfo[2][x][5]
        array[x][2][0] = outinfo[2][x][6]
        array[x][2][1] = outinfo[2][x][7]  
        array[x][2][2] = outinfo[2][x][8]
        
    return array
    
def RAMAN_POWDER(menergy,rarray,laser,width,T):
    """
    Author: Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
    
    Purpose: To calculate the Raman Spectrum of a power sample over a frequency 
    range a little larger than the frequency range of the phonon modes.
    """
    freqsteps = 1000
    ramanplot = np.zeros(shape = (freqsteps,2))
    #Compute Bose Factor
    mbose = np.zeros(len(menergy))
    for i in range(len(menergy)):
        if menergy[i][1] == 0.0:
            mbose[i] = 0.0
        else:
            mbose[i] = 1.0/(np.exp(hplank*menergy[i][1]/(kb*T))-1.0)
    
    #convert Laser frequency to Hartree.
    las = laser*cm1_to_hartree
    
    #Find range on frequencies with 10% buffer on each side
    minf = 0.0 
    maxf = 1.0
    for i in range(len(menergy)):
        if minf >= menergy[i][1]:
            minf = menergy[i][1]*0.95
        elif maxf <= menergy[i][1]:
            maxf = menergy[i][1]*1.05
    
    #Calculate Raman Spectrum
    for i in range(freqsteps):
        freq = minf+i*(maxf-minf)/freqsteps # Frequency in Hz
        ramanplot[i][0] = freq*Hz_to_Ha
        intstep = 0.0
        #counts modes
        for j in range(len(menergy)):
            spread = width/((freq*Hz_to_Ha-menergy[j][0])**2+width**2)
            G0 = 1.0/3.0*(rarray[j][0][0]+rarray[j][1][1]+rarray[j][2][2])**2
            G1 = 1.0/2.0*((rarray[j][0][1]-rarray[j][1][2])**2+(rarray[j][1][2]-rarray[j][2][0])**2+(rarray[j][2][0]-rarray[j][0][1])**2)
            G2 = 1.0/2.0*((rarray[j][0][1]+rarray[j][1][2])**2 +(rarray[j][1][2]+rarray[j][2][0])**2+(rarray[j][2][0]+rarray[j][0][1])**2) + 1.0/3.0*((rarray[j][0][0]-rarray[j][1][1])**2+(rarray[j][1][1]-rarray[j][2][2])**2+(rarray[j][2][2]-rarray[j][0][0])**2)
            Ipar  = (np.pi/menergy[j][0])*(menergy[j][0] - las)**4/clight**4*(mbose[j]+1.0)*(10.0*G0+4.0*G2)    
            Iperp = (np.pi/menergy[j][0])*(menergy[j][0] - las)**4/clight**4*(mbose[j]+1.0)*(5.0*G1+3.0*G2)
            Gterm = Ipar +Iperp
            if menergy[j][0]==0.0:
                intstep += 0.0
            else:
                intstep += spread*Gterm   
        ramanplot[i][1] = intstep
        
    return ramanplot
    
def RAMAN_POLAR(menergy,rarray,laser,option,width,T):
    """
    Author: Nicholas Pike
    Email: Nicholas.pike@ulg.ac.be
    
    Purpose: Calculate the Raman Spectrum for light polarized along 1 and 
    measured along a polarization along 2.
    """
    freqsteps = 1000
    #Compute Bose Factor
    mbose = np.zeros(len(menergy))
    for i in range(len(menergy)):
        if menergy[i][1] == 0.0:
            mbose[i] = 0.0
        else:
            mbose[i] = 1.0/(np.exp(hplank*menergy[i][1]/(kb*T))-1.0)
    
    #convert Laser frequency to Hartree.
    las = laser*cm1_to_hartree
    
    #Find range on frequencies with 10% buffer on each side
    minf = 0.0 
    maxf = 1.0
    for i in range(len(menergy)):
        if minf >= menergy[i][1]:
            minf = menergy[i][1]*0.95
        elif maxf <= menergy[i][1]:
            maxf = menergy[i][1]*1.05
    Inpol = []
    Outpol = []
    #Determine polarization directions from "option"
    if  option == 'XX':
        Inpol    = np.array([[1],[0],[0]])
        Outpol   = np.array([1,0,0])
        ramanplot = np.zeros(shape = (freqsteps,2))

        #Calculate Raman Spectrum
        for i in range(freqsteps):
            freq = minf+i*(maxf-minf)/freqsteps # Frequency in Hz
            ramanplot[i][0] = freq*Hz_to_Ha
            intstep = 0.0
            #counts modes
            for j in range(len(menergy)):
                spread = width/((freq*Hz_to_Ha-menergy[j][0])**2+width**2)
                modearray = np.array([[rarray[j][0][0],rarray[j][0][1],rarray[j][0][2]],[rarray[j][1][0],rarray[j][1][1],rarray[j][1][2]],[rarray[j][2][0],rarray[j][2][1],rarray[j][2][2]]])
                Gterm = np.dot(Outpol,np.dot(modearray,Inpol))
                if menergy[j][0]==0.0:
                    intstep += 0.0
                else:
                    intstep += spread*Gterm[0]**2*(1.0/(2.0*menergy[j][0]))*(menergy[j][0] - las)**4/clight**4*(mbose[j]+1.0) 
                ramanplot[i][1] = intstep
        
      
    elif option == 'XY':
        Inpol    = np.array([[1],[0],[0]])
        Outpol   = np.array([0,1,0]) 
        ramanplot = np.zeros(shape = (freqsteps,2))

        #Calculate Raman Spectrum
        for i in range(freqsteps):
            freq = minf+i*(maxf-minf)/freqsteps # Frequency in Hz
            ramanplot[i][0] = freq*Hz_to_Ha
            intstep = 0.0
            #counts modes
            for j in range(len(menergy)):
                spread = width/((freq*Hz_to_Ha-menergy[j][0])**2+width**2)
                modearray = np.array([[rarray[j][0][0],rarray[j][0][1],rarray[j][0][2]],[rarray[j][1][0],rarray[j][1][1],rarray[j][1][2]],[rarray[j][2][0],rarray[j][2][1],rarray[j][2][2]]])
                Gterm = np.dot(Outpol,np.dot(modearray,Inpol))
                if menergy[j][0]==0.0:
                    intstep += 0.0
                else:
                    intstep += spread*Gterm[0]**2*(1.0/(2.0*menergy[j][0]))*(menergy[j][0] - las)**4/clight**4*(mbose[j]+1.0) 
                ramanplot[i][1] = intstep
                
    elif option == 'XZ':
        Inpol    = np.array([[1],[0],[0]])
        Outpol   = np.array([0,0,1]) 
        ramanplot = np.zeros(shape = (freqsteps,2))

        #Calculate Raman Spectrum
        for i in range(freqsteps):
            freq = minf+i*(maxf-minf)/freqsteps # Frequency in Hz
            ramanplot[i][0] = freq*Hz_to_Ha
            intstep = 0.0
            #counts modes
            for j in range(len(menergy)):
                spread = width/((freq*Hz_to_Ha-menergy[j][0])**2+width**2)
                modearray = np.array([[rarray[j][0][0],rarray[j][0][1],rarray[j][0][2]],[rarray[j][1][0],rarray[j][1][1],rarray[j][1][2]],[rarray[j][2][0],rarray[j][2][1],rarray[j][2][2]]])
                Gterm = np.dot(Outpol,np.dot(modearray,Inpol))
                if menergy[j][0]==0.0:
                    intstep += 0.0
                else:
                    intstep += spread*Gterm[0]**2*(1.0/(2.0*menergy[j][0]))*(menergy[j][0] - las)**4/clight**4*(mbose[j]+1.0) 
                ramanplot[i][1] = intstep
                
    elif option == 'YY':
        Inpol    = np.array([[0],[1],[0]])
        Outpol   = np.array([0,1,0]) 
        ramanplot = np.zeros(shape = (freqsteps,2))

        #Calculate Raman Spectrum
        for i in range(freqsteps):
            freq = minf+i*(maxf-minf)/freqsteps # Frequency in Hz
            ramanplot[i][0] = freq*Hz_to_Ha
            intstep = 0.0
            #counts modes
            for j in range(len(menergy)):
                spread = width/((freq*Hz_to_Ha-menergy[j][0])**2+width**2)
                modearray = np.array([[rarray[j][0][0],rarray[j][0][1],rarray[j][0][2]],[rarray[j][1][0],rarray[j][1][1],rarray[j][1][2]],[rarray[j][2][0],rarray[j][2][1],rarray[j][2][2]]])
                Gterm = np.dot(Outpol,np.dot(modearray,Inpol))
                if menergy[j][0]==0.0:
                    intstep += 0.0
                else:
                    intstep += spread*Gterm[0]**2*(1.0/(2.0*menergy[j][0]))*(menergy[j][0] - las)**4/clight**4*(mbose[j]+1.0) 
                ramanplot[i][1] = intstep
                
    elif option == 'YZ':
        Inpol    = np.array([[0],[1],[0]])
        Outpol   = np.array([0,0,1]) 
        ramanplot = np.zeros(shape = (freqsteps,2))

        #Calculate Raman Spectrum
        for i in range(freqsteps):
            freq = minf+i*(maxf-minf)/freqsteps # Frequency in Hz
            ramanplot[i][0] = freq*Hz_to_Ha
            intstep = 0.0
            #counts modes
            for j in range(len(menergy)):
                spread = width/((freq*Hz_to_Ha-menergy[j][0])**2+width**2)
                modearray = np.array([[rarray[j][0][0],rarray[j][0][1],rarray[j][0][2]],[rarray[j][1][0],rarray[j][1][1],rarray[j][1][2]],[rarray[j][2][0],rarray[j][2][1],rarray[j][2][2]]])
                Gterm = np.dot(Outpol,np.dot(modearray,Inpol))
                if menergy[j][0]==0.0:
                    intstep += 0.0
                else:
                    intstep += spread*Gterm[0]**2*(1.0/(2.0*menergy[j][0]))*(menergy[j][0] - las)**4/clight**4*(mbose[j]+1.0) 
                ramanplot[i][1] = intstep
                
    elif option == 'ZZ':
        Inpol    = np.array([[0],[0],[1]])
        Outpol   = np.array([0,0,1]) 
        ramanplot = np.zeros(shape = (freqsteps,2))

        #Calculate Raman Spectrum
        for i in range(freqsteps):
            freq = minf+i*(maxf-minf)/freqsteps # Frequency in Hz
            ramanplot[i][0] = freq*Hz_to_Ha
            intstep = 0.0
            #counts modes
            for j in range(len(menergy)):
                spread = width/((freq*Hz_to_Ha-menergy[j][0])**2+width**2)
                modearray = np.array([[rarray[j][0][0],rarray[j][0][1],rarray[j][0][2]],[rarray[j][1][0],rarray[j][1][1],rarray[j][1][2]],[rarray[j][2][0],rarray[j][2][1],rarray[j][2][2]]])
                Gterm = np.dot(Outpol,np.dot(modearray,Inpol))
                if menergy[j][0]==0.0:
                    intstep += 0.0
                else:
                    intstep += spread*Gterm[0]**2*(1.0/(2.0*menergy[j][0]))*(menergy[j][0] - las)**4/clight**4*(mbose[j]+1.0) 
                ramanplot[i][1] = intstep
                
    elif option == 'THETA':
        #Raman Spectrum is calculated twice, once for x and once for y initial
        #polarization with the outgoing polarization at an angle theta
        thetasteps = 720
        ramanplot = np.zeros(shape = (thetasteps,len(menergy)+1))
        #Calculate Raman Spectrum
        for i in range(thetasteps):
            theta = 0.0+i*2.0*np.pi/thetasteps
            ramanplot[i][0] = theta
            intstep = 0.0
            Inpol1 = np.array([[1],[0],[0]])
            Inpol2 = np.array([[0],[1],[0]])
            Outpol = np.array([np.cos(theta),np.sin(theta),0])
            for j in range(len(menergy)):
                modearray = -np.array([[rarray[j][0][0],rarray[j][0][1],rarray[j][0][2]],[rarray[j][1][0],rarray[j][1][1],rarray[j][1][2]],[rarray[j][2][0],rarray[j][2][1],rarray[j][2][2]]])
                Gterm1 = np.dot(Outpol,np.dot(modearray,Inpol1))
                Gterm2 = np.dot(Outpol,np.dot(modearray,Inpol2))
                if menergy[j][0]==0.0:
                    intstep = 0.0
                else:
                    intstep = (Gterm1[0]**2+Gterm2[0]**2)*(1.0/(2.0*menergy[j][0]))*(menergy[j][0] - las)**4/clight**4*(mbose[j]+1.0)
                ramanplot[i][j+1] = intstep

    else:
        print('How did this happen?') #Should be impossible to get here
    
    return ramanplot
    
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
           
    for i in range(len(keywords)):
        if keywords[i]== 'True' and i == 1:
            #run module to print Dielectric tensor as a function of frequency
            outfile = 'dielectric.freq'
            printoutfile('# Freq    dielectric      Reflectivity',outfile)
            printoutfile('#           x y z x y z  ',outfile)
            for j in range(int(keywords[2])):
                printoutfile('%e %e %e %e %e %e %e' %(output[3][j][0],output[3][j][1],output[3][j][2],output[3][j][3],output[3][j][4],output[3][j][5],output[3][j][6]),outfile)
        elif keywords[i]== 'True' and i == 3:
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
    
    outname        = 'Ramanspec.out'
    cm1_to_hartree = 4.55633E-6         # conversion factor between cm-1 and Hartree
    cm1_to_hz      = 2.99793E10         # conversion factor between cm-1 and Hz
    Hz_to_Ha       = 1.519828500716E-16 # conversion factor between Hz and Hartree 
    kb             = 8.6173324E-5       # Boltzmann constant in eV/K
    T              = 0                  # Default Temperature (user input variable)
    hplank         = 6.58211928E-15     # h in eVs
    clight         = 137.0359997566     # Speed of light in atomic units
    width          = 0.0                # spread of lorentian (user input variable)

    if any('SPYDER' in name for name in os.environ):
        user_inputfile = 'input_raman'
    else:
        #Determines what the program is to do if it is run from command line
        user_inputfile = sys.argv[1] #input should be python program_name tfile 
          
    #Name of default input file 
    vararray = READ_INPUT(user_inputfile,outname) #The program will look for the specified input file
    
    
    #Print header file and start the calculation
    PRINT_HEADER()
    
    #Read anaddb.out file and determine what type of calculation to run
    DETER_MENU(vararray)
    
    
#Ends program
