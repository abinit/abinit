#=================================
# AbinitStructureViewer.py
version = 'beta'
#=================================
# written by Benjamin Tardif
# benjamin.tardif@umontreal.ca

# last modified : october 25 2006 by Guillaume Dumont
# added support for default rprim
#=================================
headline = '\n==========================\n AbinitStructureViewer.py\n version %s\n==========================' %version


#=====================================================================================================================================================================
#IMPORTS
import os
import sys
from Numeric import *
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#VARIABLES
class VariableContainer:pass

#control variables
ctrl = VariableContainer()
ctrl.arglist = list()                  # list of strings containing the input arguments and keywords in the command line
ctrl.arglist = sys.argv[1:]
ctrl.validkeywords = list()            # list of strings containing the valid keywords who can be used in the command line
ctrl.validkeywords = ['-debug','-ldebug','-setup']
ctrl.debugmode = bool()                # True if debug mode is activated (by adding keyword "-debug" in the command line)
ctrl.ldebugmode = bool()               # True if ldebug mode is activated (by adding keyword "-ldebug" in the command line)
ctrl.autolaunch = str()                # 'yes': the autolaunch mode is activated (can be set using the -setup keyword)
                                       # 'no': otherwise
ctrl.launchcommand = str()             # the command used to launch the jmol application
ctrl.defaultlaunchcommand = str()      # the default launch command used to launch the jmol application
ctrl.defaultlaunchcommand = 'java -jar /Applications/jmol-10.00/Jmol.jar'
ctrl.changelaunchcommand = str()       # 'yes': the user want to change the launch command (while in setup mode)
                                       # 'no': the user do not want to change the launch command (while in setup mode)
ctrl.filename = str()                  # name of the file treated
ctrl.filetype = str()                  # 'out': the file correspond to a .out file
                                       # 'log': otherwise
ctrl.filedata = list()                 # list of strings each containing one line of the file
ctrl.relaxationtype = int()            # 0: no relaxation no relaxation (ionmov == 0 and optcell == 0)
                                       # 1: atoms relaxation (ionmov != 0 and optcell == 0)
                                       # 2: atoms and cell relaxation (ionmov != 0 and optcell != 0)
ctrl.bohrtoangst = float()             # conversion factor between bohrs and angstroms (angstrom/bohr)
ctrl.bohrtoangst = float(0.5291772108)
ctrl.completed = bool()                # True if the calculation was completed in the file
ctrl.askreplicate = str()              # 3 integers separated by spaces indicating the cell replication
ctrl.replicatevalidator = bool()       # True if the format of askreplicate is valid
ctrl.replicate = zeros(3,Int)          # integer array of size 1x3 containing the number of times the cell is replicated in each direction
ctrl.natomreplicated = int()           # total number of atoms displayed, including the replicated ones 
ctrl.xyzfilename = str()               # name of the file created
ctrl.flagacell = bool()                # True when acell is assigned
ctrl_flagrprim = bool()                # True when rprim is assigned
ctrl.flagxangst = bool()               # True when xangst is assigned
ctrl.periodictable = {}                # dictionary mapping the atomic number to the atomic symbol
ctrl.periodictable = {
      1:'H ',                                                                                                                                                  2:'He',
      3:'Li',  4:'Be',                                                                                            5:'B ',  6:'C ',  7:'N ',  8:'O ',  9:'F ', 10:'Ne',
     11:'Na', 12:'Mg',                                                                                           13:'Al', 14:'Si', 15:'P ', 16:'S ', 17:'Cl', 18:'Ar',
     19:'K ', 20:'Ca', 21:'Sc', 22:'Ti', 23:'V ', 24:'Cr', 25:'Mn', 26:'Fe', 27:'Co', 28:'Ni', 29:'Cu', 30:'Zn', 31:'Ga', 32:'Ge', 33:'As', 34:'Se', 35:'Br', 36:'Kr',
     37:'Rb', 38:'Sr', 39:'Y ', 40:'Zr', 41:'Nb', 42:'Mo', 43:'Tc', 44:'Ru', 45:'Rh', 46:'Pd', 47:'Ag', 48:'Cd', 49:'In', 50:'Sn', 51:'Sb', 52:'Te', 53:'I ', 54:'Xe',
     55:'Cs', 56:'Ba',          72:'Hf', 73:'Ta', 74:'W ', 75:'Re', 76:'Os', 77:'Ir', 78:'Pt', 79:'Au', 80:'Hg', 81:'Tl', 82:'Pb', 83:'Bi', 84:'Po', 85:'At', 86:'Rn',
     87:'Fr', 88:'Ra',         104:'Ku',105:'Ha',106:'Unh',107:'Uns',108:'Uno',109:'Une',

                                57:'La', 58:'Ce', 59:'Pr', 60:'Nd', 61:'Pm', 62:'Sm', 63:'Eu', 64:'Gd', 65:'Tb', 66:'Dy', 67:'Ho', 68:'Er', 69:'Tm', 70:'Yb', 71:'Lu',
                                89:'Ac', 90:'Th', 91:'Pa', 92:'U ', 93:'Np', 94:'Pu', 95:'Am', 96:'Cm', 97:'Bk', 98:'Cf', 99:'Es',100:'Fm',101:'Md',102:'No',103:'Lr'}

#data variables
data = VariableContainer()
data.natom = int()                     # number of atoms in the unit cell
data.ntypat = int()                    # number of types of atoms
data.typat = list()                    # list of integers corresponding to the type of each atom
data.ionmov = int()                    # indicates the algorithm used to optimize the atomic positions (see Abinit variable "ionmov")
data.znucl = list()                    # list in integers corresponding to the atomic number of each atom
data.optcell = int()                   # indicates the algorithm used to optimize the cell geometry (see Abinit variable "optcell")
data.acell = list()                    # list of array(3,Float) corresponding to the scale cell (in bohrs) (see Abinit variable "acell")
data.rprim = list()                    # list of array((3,3),Float) corresponding the the real space primitive translations (see Abinit variable "rprim")
data.a1 = list()                       # list of array(3,Float) corresponding to the 1st translation vector (in angstroms)
data.a2 = list()                       # list of array(3,Float) corresponding to the 2nd translation vector (in angstroms)
data.a3 = list()                       # list of array(3,Float) corresponding to the 3rd translation vector (in angstroms)
data.xangst = list()                   # list of array((data.natom,3),Float) corresponding to the atomic positions (in angstroms)
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#METHODS
def detectfile(filename,directory): # type(filename) = type(directory) = string
    # method detectfile returns True if the filename is found in the specified directory 
    if filename in os.listdir(directory):
        return True
    else:
        return False

def clean(list): # type(list) = list of strings
    # method clean removes character strings '\n' and '\r' and empty lines from a list
    # (the list is usually obtained with the ".readlines()" method)
    L = len(list)
    for i in range(L):
        list[L-1-i] = list[L-1-i].replace('\n','')
        list[L-1-i] = list[L-1-i].replace('\r','')
        if list[L-1-i].split() == []:
            list.pop(L-1-i)

def rmdoubleentries(list): # type(list) = any list of paired repeated elements
    # method rmdoubleentries removes paired repeated entries of a list and keeps only one copy of each
    # example : rmdoubleentries([1,1,2,2,3,3,4,4,5,5]) = [1,2,3,4,5]
    # example : rmdoubleentries([1,1,2,2,3,3,4,4,5])   = [1,2,3,4,5]
    for i in range(int(len(list)/2)):
        list.pop(-2-i)
#=====================================================================================================================================================================


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#MAIN
print headline


#=====================================================================================================================================================================
#AUTOLAUNCH FILE

if detectfile('.AbinitStructureViewer_autolaunch',sys.path[0]) == False:
    #autolaunch file not found, create a default one
    writer = open(sys.path[0]+'/.AbinitStructureViewer_autolaunch','w')
    writer.write('this file is used by the program AbinitStructureViewer.py (version %s)\n' %version)
    writer.write('(this file is not essential and can be deleted if needed)\n\n')
    writer.write('jmol launch command :\n%s\n\n' %ctrl.defaultlaunchcommand)
    writer.write('automatically launch jmol each time a .xyz file is created ?\nno')
    writer.close()
else:
    #get checkversion, launchcommand, autolaunch
    reader = open(sys.path[0]+'/.AbinitStructureViewer_autolaunch','r')
    launchfile = reader.readlines()
    reader.close()
    clean(launchfile)
    checkversion = launchfile[0].split()[-1].split(')')[0]
    launchcommand = launchfile[3]
    autolaunch = launchfile[5]
    if checkversion != version:
        #versions are different, update the autolaunch file
        writer = open(sys.path[0]+'/.AbinitStructureViewer_autolaunch','w')
        writer.write('this file is used by the program AbinitStructureViewer.py (version %s)\n' %version)
        writer.write('(this file is not essential and can be deleted if needed)\n\n')
        writer.write('jmol launch command :\n%s\n\n' %launchcommand)
        writer.write('automatically launch jmol each time a .xyz file is created ?\n%s' %autolaunch)
        writer.close()
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#COMMAND LINE

#abort if a keyword is not valid
for arg in ctrl.arglist:
    if arg[0] == '-':
        # a keyword is found
        if arg not in ctrl.validkeywords:
            # the keyword is not valid
            print '\n- operation aborted -\n\n%s is not a valid keyword' %arg
            print '\nvalid keywords are : -debug, -ldebug, -setup\n'
            sys.exit()

#abort if a keyword is repeated
for keyword in ctrl.validkeywords:
    if ctrl.arglist.count(keyword) > 1:
        print '\n- operation aborted -\n\nkeyword %s is repeated %s times\n' %(keyword,ctrl.arglist.count(keyword))
        sys.exit()

#keyword -setup
if '-setup' in ctrl.arglist:
    # user activated the autolaunch setup
    ctrl.arglist.pop(ctrl.arglist.index('-setup'))

    # get launchcommand
    reader = open(sys.path[0]+'/.AbinitStructureViewer_autolaunch','r')
    launchfile = reader.readlines()
    reader.close()
    clean(launchfile)
    ctrl.launchcommand = launchfile[3]

    # change launchcommand
    print '\ncurrent jmol launch command is :\n%s\n' %ctrl.launchcommand
    while ctrl.changelaunchcommand not in ['yes','no']:    
        ctrl.changelaunchcommand = raw_input('do you wish to change it (yes ; no) ? ')

    if ctrl.changelaunchcommand == 'yes':
        ctrl.launchcommand = raw_input('\nenter the new jmol launch command :\n')

    # change autolaunch
    ctrl.autolaunch = raw_input('\nautomatically launch jmol each time a .xyz file is created (yes ; no) ? ')
    while ctrl.autolaunch not in ['yes','no']:
        ctrl.autolaunch = raw_input('automatically launch jmol each time a .xyz file is created (yes ; no) ? ')

    # overwrite autolaunch file
    writer = open(sys.path[0]+'/.AbinitStructureViewer_autolaunch','w')
    writer.write('this file is used by the program AbinitStructureViewer.py (version %s)\n' %version)
    writer.write('(this file is not essential and can be deleted if needed)\n\n')
    writer.write('jmol launch command :\n%s\n\n' %ctrl.launchcommand)
    writer.write('automatically launch jmol each time a .xyz file is created ?\n%s' %ctrl.autolaunch)
    writer.close()

    print '\n- modifications done -\n'
    sys.exit()

#keyword -debug
if '-debug' in ctrl.arglist:
    # user activated the debug mode    
    ctrl.debugmode=True
    ctrl.arglist.pop(ctrl.arglist.index('-debug'))

#keyword -ldebug
if '-ldebug' in ctrl.arglist:
    # user activated the ldebug mode    
    ctrl.debugmode=True
    ctrl.ldebugmode=True
    ctrl.arglist.pop(ctrl.arglist.index('-ldebug'))

#(put additionnal keywords here)

#get filename
if len(ctrl.arglist) == 0:
    # user entered no filename in the command line
    ctrl.filename = raw_input('\nEnter the filename : \n')
elif len(ctrl.arglist) == 1:
    # user entered the filename in the command line
    ctrl.filename = ctrl.arglist[0]
elif len(ctrl.arglist) > 1:
    # user entered too much arguments in the command line
    print '\n- too much arguments entered in the command line -\n'
    sys.exit()

#abort if the file does not exists
if detectfile(ctrl.filename,'.') == False:
    print '\n"%s" - file not found -\n' %ctrl.filename
    sys.exit()

#activate debugmode
if ctrl.debugmode==True:
    print '\n- DEBUG MODE -'
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#READ THE FILE

#read file and acquire data
reader = open(ctrl.filename,"r")
ctrl.filedata = reader.readlines()
reader.close()
clean(ctrl.filedata)

#compute filetype
if ctrl.filename.split('.')[-1][:3] == 'out':
    ctrl.filetype = 'out'
else:
    ctrl.filetype = 'log'
if ctrl.debugmode==True:print 'filetype : "%s" will be treated as a < %s > file' %(ctrl.filename,ctrl.filetype)
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#EXTRACT DATA FROM THE FILE

#get natom
for line in ctrl.filedata:
    if line.split()[0] == 'natom':
        data.natom = int(line.split()[1])
if ctrl.debugmode==True:print 'natom = %s' %data.natom

#get ntypat
for line in ctrl.filedata:
    if line.split()[0] == 'ntypat':
        data.ntypat = int(line.split()[1])
if ctrl.debugmode==True:print 'ntypat = %s' %data.ntypat

#get typat
for i in range(len(ctrl.filedata)):
    if ctrl.filedata[i].split()[0] == 'typat':
        k = i
        while len(data.typat) < int(data.natom):
            for j in range(len(ctrl.filedata[k].split())):
                if ctrl.filedata[k].split()[j] != 'typat':
                    data.typat.append(int(ctrl.filedata[k].split()[j]))
            k = k+1
if ctrl.debugmode==True:print 'typat = %s' %data.typat

#get znucl
for i in range(len(ctrl.filedata)):
    if ctrl.filedata[i].split()[0] == 'znucl':
        k = i
        while len(data.znucl) < int(data.ntypat):
            for j in range(len(ctrl.filedata[k].split())):
                if ctrl.filedata[k].split()[j] != 'znucl':
                    data.znucl.append(int(float(ctrl.filedata[k].split()[j])))
            k = k+1
if ctrl.debugmode==True:print 'znucl = %s' %data.znucl

#abort if znucl is not present in the file
if data.znucl == []:
    if ctrl.filetype == 'out':
        print '\n- operation aborted -\n"znucl" not found in the given file\n'
    elif ctrl.filetype == 'log':
        print '\n- operation aborted -\n"znucl" not found in the given file\nmaybe "%s" is not a valid log file\n' %ctrl.filename
    sys.exit()

#get ionmov
for line in ctrl.filedata:
    if line.split()[0] == 'ionmov':
        data.ionmov = int(line.split()[1])
if ctrl.debugmode==True:print 'ionmov = %s' %data.ionmov

#get optcell
for line in ctrl.filedata:
    if line.split()[0] == 'optcell':
        data.optcell = int(line.split()[1])
if ctrl.debugmode==True:print 'optcell = %s' %data.optcell

#compute relaxationtype
if data.ionmov == 0 and data.optcell == 0:
    ctrl.relaxationtype = 0
    if ctrl.debugmode==True:print 'relaxationtype : no relaxation (ionmov == 0 and optcell == 0)'

if data.ionmov != 0 and data.optcell == 0:
    ctrl.relaxationtype = 1
    if ctrl.debugmode==True:print 'relaxationtype : atoms relaxation (ionmov != 0 and optcell == 0)'

if data.ionmov != 0 and data.optcell != 0:
    ctrl.relaxationtype = 2
    if ctrl.debugmode==True:print 'relaxationtype : atoms and cell relaxation (ionmov != 0 and optcell != 0)'

#get acell, rprim, xangst
#----------relaxationtype=0----------
if ctrl.relaxationtype == 0:
    for line in ctrl.filedata:
        if ctrl.flagacell==False:
            if line.split()[0] == 'acell':
                data.acell.append(ctrl.bohrtoangst*array([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])],Float))
                ctrl.flagacell=True
    
    for i in range(len(ctrl.filedata)):
        if ctrl_flagrprim==False:
            if ctrl.filedata[i].split()[0] == 'rprim':
                data.rprim.append(array([\
                    [float(ctrl.filedata[i  ].split()[1]),float(ctrl.filedata[i  ].split()[2]),float(ctrl.filedata[i  ].split()[3])],\
                    [float(ctrl.filedata[i+1].split()[0]),float(ctrl.filedata[i+1].split()[1]),float(ctrl.filedata[i+1].split()[2])],\
                    [float(ctrl.filedata[i+2].split()[0]),float(ctrl.filedata[i+2].split()[1]),float(ctrl.filedata[i+2].split()[2])]],Float))
                ctrl_flagrprim=True
    
    if ctrl_flagrprim==False:
        data.rprim = [ array([[1.0 , 0.0, 0.0], [0.0 , 1.0, 0.0], [0.0 , 0.0, 1.0]]) ]*len(data.acell)
	
    for k in range(len(ctrl.filedata)):
        if ctrl.flagxangst==False:
            if ctrl.filedata[k].split()[0] == 'xcart':
                data.xangst.append(zeros((data.natom,3),Float))
                data.xangst[-1][0][0] = ctrl.bohrtoangst*float(ctrl.filedata[k].split()[1])
                data.xangst[-1][0][1] = ctrl.bohrtoangst*float(ctrl.filedata[k].split()[2])
                data.xangst[-1][0][2] = ctrl.bohrtoangst*float(ctrl.filedata[k].split()[3])
                for i in range(data.natom-1):
                    for j in range(3):
                        data.xangst[-1][i+1][j] = ctrl.bohrtoangst*float(ctrl.filedata[k+1+i].split()[j])
                ctrl.flagxangst=True

    if data.xangst == []:
            print '\n- operation aborted -\n"xcart" not found in the given file'
            sys.exit()

#----------relaxationtype=1----------
if ctrl.relaxationtype == 1:
    for line in ctrl.filedata:
        if ctrl.flagacell==False:
            if line.split()[0] == 'acell':
                data.acell.append(ctrl.bohrtoangst*array([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])],Float))
                ctrl.flagacell=True
                
    for i in range(len(ctrl.filedata)):
        if ctrl_flagrprim==False:
            if ctrl.filedata[i].split()[0] == 'rprim':
                data.rprim.append(array([\
                    [float(ctrl.filedata[i  ].split()[1]),float(ctrl.filedata[i  ].split()[2]),float(ctrl.filedata[i  ].split()[3])],\
                    [float(ctrl.filedata[i+1].split()[0]),float(ctrl.filedata[i+1].split()[1]),float(ctrl.filedata[i+1].split()[2])],\
                    [float(ctrl.filedata[i+2].split()[0]),float(ctrl.filedata[i+2].split()[1]),float(ctrl.filedata[i+2].split()[2])]],Float))
                ctrl_flagrprim=True
		
    if ctrl_flagrprim==False:
        data.rprim = [ array([[1.0 , 0.0, 0.0], [0.0 , 1.0, 0.0], [0.0 , 0.0, 1.0]]) ]*len(data.acell)
	
			      
    for k in range(len(ctrl.filedata)):
        if ctrl.filedata[k] == ' Cartesian coordinates (bohr)':
            data.acell.append(data.acell[-1])
            data.rprim.append(data.rprim[-1])
            data.xangst.append(zeros((data.natom,3),Float))
            for i in range(data.natom):
                for j in range(3):
                    data.xangst[-1][i][j] = ctrl.bohrtoangst*float(ctrl.filedata[k+1+i].split()[j])

    if data.xangst == []:
        print '\n- operation aborted -\n"Cartesian coordinates" not found in the given file'
        sys.exit()

    data.acell.pop(0)
    data.rprim.pop(0)

    if ctrl.filetype == 'log':
        rmdoubleentries(data.acell)
        rmdoubleentries(data.rprim)
        rmdoubleentries(data.xangst)

#----------relaxationtype=2----------
if ctrl.relaxationtype == 2:
    for line in ctrl.filedata:
        if line.split()[0] == 'acell=':
            data.acell.append(ctrl.bohrtoangst*array([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])],Float))

    for i in range(len(ctrl.filedata)):
        if ctrl.filedata[i].split()[0] == 'rprim=':
            data.rprim.append(array([\
                [float(ctrl.filedata[i  ].split()[1]),float(ctrl.filedata[i  ].split()[2]),float(ctrl.filedata[i  ].split()[3])],\
                [float(ctrl.filedata[i+1].split()[0]),float(ctrl.filedata[i+1].split()[1]),float(ctrl.filedata[i+1].split()[2])],\
                [float(ctrl.filedata[i+2].split()[0]),float(ctrl.filedata[i+2].split()[1]),float(ctrl.filedata[i+2].split()[2])]],Float))

    if ctrl_flagrprim==False:
        data.rprim = [ array([[1.0 , 0.0, 0.0], [0.0 , 1.0, 0.0], [0.0 , 0.0, 1.0]]) ]*len(data.acell)
	
    for k in range(len(ctrl.filedata)):
        if ctrl.filedata[k] == ' Cartesian coordinates (bohr)':
            data.xangst.append(zeros((data.natom,3),Float))
            for i in range(data.natom):
                for j in range(3):
                    data.xangst[-1][i][j] = ctrl.bohrtoangst*float(ctrl.filedata[k+1+i].split()[j])

    if data.xangst == []:
        print '\n- operation aborted -\n"Cartesian coordinates" not found in the given file'
        sys.exit()

    while len(data.acell) != len(data.xangst):
        data.acell.pop(-1)
    while len(data.rprim) != len(data.xangst):
        data.rprim.pop(-1)

    if ctrl.filetype == 'log':
        rmdoubleentries(data.acell)
        rmdoubleentries(data.rprim)
        rmdoubleentries(data.xangst)
#------------------------------------

if ctrl.debugmode==True and ctrl.ldebugmode==False:
    print 'acell  = [%s element(s)]' %len(data.acell)
    print 'rprim  = [%s element(s)] for a more detailed version, type "-ldebug" instead of "-debug" in the command line' %len(data.rprim)
    print 'xangst = [%s element(s)]' %len(data.xangst)

if ctrl.debugmode==True and ctrl.ldebugmode==True:
    print 'acell  = [%s element(s)]' %len(data.acell)
    print 'rprim  = [%s element(s)]' %len(data.rprim)
    print 'xangst = [%s element(s)]' %len(data.xangst)

    L = len(data.acell)
    for i in range(L):
        print '\nacell %s of %s :\n%s' %(i+1,L,data.acell[i])
    L = len(data.rprim)
    for i in range(L):
        print '\nrprim %s of %s :\n%s' %(i+1,L,data.rprim[i])
    L = len(data.xangst)
    for i in range(L):
        print '\nxangst %s of %s :\n%s' %(i+1,L,data.xangst[i])


#compute primitive vectors (in angstroms)
for i in range(len(data.acell)):
    data.a1.append(data.acell[i][0]*data.rprim[i][0])
    data.a2.append(data.acell[i][1]*data.rprim[i][1])
    data.a3.append(data.acell[i][2]*data.rprim[i][2])

#display number of configurations found
print '\n%s configuration(s) found' %len(data.acell)
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#WARNING IF CALCULATION NOT COMPLETED

for line in ctrl.filedata:
    if line == ' Calculation completed.':
        ctrl.completed=True

if ctrl.completed == False:
    print '- WARNING : calculation not completed in the given file -' 
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#REPLICATION PARAMETERS

#user entered replication parameters
while ctrl.replicatevalidator == False:
    ctrl.askreplicate = raw_input('\nEnter the number of times you wish to replicate the primitive cell\n(If you only want the primitive cell, enter : 1 1 1)\n')
    if len(ctrl.askreplicate.split()) != 3:
        ctrl.replicatevalidator = False
    elif ctrl.askreplicate.split()[0].isdigit() + ctrl.askreplicate.split()[1].isdigit() + ctrl.askreplicate.split()[2].isdigit() != 3:
        ctrl.replicatevalidator = False
    elif int(ctrl.askreplicate.split()[0]) == 0 or int(ctrl.askreplicate.split()[1]) == 0 or int(ctrl.askreplicate.split()[2]) == 0:
        ctrl.replicatevalidator = False
    else:
        ctrl.replicatevalidator=True
    if ctrl.replicatevalidator == False:
        print '- invalid entry - (enter 3 non zero integers separated by a space)'

#computed replication parameters
ctrl.replicate = array([int(ctrl.askreplicate.split()[0]),int(ctrl.askreplicate.split()[1]),int(ctrl.askreplicate.split()[2])],Int)
ctrl.natomreplicated = int(data.natom*ctrl.replicate[0]*ctrl.replicate[1]*ctrl.replicate[2])
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#WRITE THE .XYZ FILE

#compute the .xyz filename
if ctrl.filetype == 'out':
    ctrl.xyzfilename = '%s_%s.xyz' %(ctrl.filename.split('.')[0],ctrl.filename.split('.')[1])
else:
    ctrl.xyzfilename = '%s.xyz' %ctrl.filename

#write the file
writer = open(ctrl.xyzfilename,"w")

configuration = 0
for i in range(len(data.xangst)):
    if ctrl.relaxationtype == 0:
        writer.write('%s\nangstrom\n' %ctrl.natomreplicated)
    else:
        writer.write('%s\nCONFIGURATION %s\n' %(ctrl.natomreplicated,configuration))
        configuration+=1
    for j in range(data.natom):
        for a in range(ctrl.replicate[0]):
            for b in range(ctrl.replicate[1]):
                for c in range(ctrl.replicate[2]):
                    atompos = data.xangst[i][j]+a*data.a1[i]+b*data.a2[i]+c*data.a3[i]
                    atomposx = '%e' %atompos[0]
                    atomposy = '%e' %atompos[1]
                    atomposz = '%e' %atompos[2]
                    if atompos[0]>=0:atomposx = ' %e' %atompos[0]
                    if atompos[1]>=0:atomposy = ' %e' %atompos[1]
                    if atompos[2]>=0:atomposz = ' %e' %atompos[2]
                    writer.write('%s   %s   %s   %s\n' %(ctrl.periodictable[data.znucl[data.typat[j]-1]],atomposx,atomposy,atomposz))
writer.close()

print '\n"%s" file created successfully\n' %ctrl.xyzfilename
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#AUTOLAUNCH

#get autolaunch and launchcommand
reader = open(sys.path[0]+'/.AbinitStructureViewer_autolaunch','r')
launchfile = reader.readlines()
reader.close()
clean(launchfile)
ctrl.autolaunch = launchfile[5]
ctrl.launchcommand = launchfile[3]

#autolaunch
if ctrl.autolaunch == 'yes':
    print 'launching jmol using command :\n%s %s &' %(ctrl.launchcommand,ctrl.xyzfilename)
    os.system('%s %s &' %(ctrl.launchcommand,ctrl.xyzfilename))
#=====================================================================================================================================================================


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
