#===================================
program = 'AbinitBandStructureMaker.py'
version = '1.3'
#===================================
# last modified : november 16 2010
# written by Benjamin Tardif
# benjamin.tardif@umontreal.ca
# Modified by Paul Boulanger 
# 1.3 : Converted to numpy
#     : corrected a bug related to termination of the segments (line 114)
#     : defined a second angletol to distinguish between specialkpt detection and segment construction
#===================================


#=====================================================================================================================================================================
#IMPORTS
import os
import re
import sys
import time
from numpy import *
#from Numeric import *
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#VARIABLES
class VariableContainer:pass

#default variables
default = VariableContainer()
default.setupfilename = '.%s_setup' %program[:-3]
default.launchcommand = 'xmgrace'
default.autolaunch = 'no'
default.energyshift = 'yes'
default.valbandcolor = 'blue'
default.conbandcolor = 'red'
default.bandlinewidth = 1
default.fermilinewidth = 1
default.separatorlinewidth = 1
default.emptyspacewidth = 10

#control variables
ctrl = VariableContainer()
ctrl.arglist = sys.argv[1:]              # list of strings containing all the input arguments in the command line
# list of valid keywords which can be used in the command line
ctrl.validkeywords = ['-setautolaunch','-setenergyshift','-setlinecolor','-setlinewidth','-setspacewidth','-setup','-setdefault','-debug'] 
ctrl.debugmode = False                   # True if debug mode is activated (by adding keyword "-debug" in the command line)
ctrl.launchcommand = str()            # string containing the xmgrace launch command
ctrl.autolaunch = str()               # 'yes'or 'no', indicating if xmgrace will be automatically launched each time a .agr file is created
ctrl.energyshift = str()              # 'yes' or 'no', indicating if the energies will be shifted to bring the fermi energy to zero
ctrl.valbandcolor = str()
ctrl.conbandcolor = str()
ctrl.bandlinewidth = int()
ctrl.fermilinewidth = int()
ctrl.separatorlinewidth = int()
ctrl.emptyspacewidth = int()


# dictionary maping color name with color number in xmgrace 
ctrl.xmgracecolor = {
'white'    : 0,
'black'    : 1,
'red'      : 2,
'green'    : 3,
'blue'     : 4,
'yellow'   : 5,
'brown'    : 6,
'grey'     : 7,
'violet'   : 8,
'cyan'     : 9,
'magenta'  :10,
'orange'   :11,
'indigo'   :12,
'maroon'   :13,
'turquoise':14,
'green4'   :15}

ctrl.filename = str()                 # name of the file entered (*.out or *.dbs)
ctrl.filetype = str()                 # 'out' or 'dbs' according to the type of the file entered
ctrl.filedata = list()                # list of strings each containing one line of the file entered

ctrl.dbsfilename = str()              # name of the file produced (*.dbs)
ctrl.agrfilename = str()              # name of the file produced (*.agr)

ctrl.angletol = 1                # maximum angle between 2 k-points under which they will be considered being in the same direction
ctrl.angletol2 = 0.1
ctrl.bandstructurescheme = str() # scheme of the type '{lA}-{lB}-{lC}' corresponding to the band structure to be plotted
                                 # each {} corresponds to the name given to a special k-point
                                 # {} 1st character : 'l' for 'letter' or 's' for 'symbol'
                                 # {} 2nd character : a letter which combined with the 1st character will give a caption associated with this special k-point 
                                 # '-' between {} indicates that a band structure must be plotted beteen the 2 corresponding k-points
                                 # ' ' (empty space) between {} indicates that an empty space must be inserted in the band structure
ctrl.dicospecialkpt = dict()     # dictionary mapping , keys=[vectors], values={lx}
ctrl.segmentcaptionlist = list() # list of lists of string, each list having the form
                                 # [{lA}, {lB}] for a segment joining 2 points
                                 # ['empty space'] for a empty space in the band structure
ctrl.spacepercent = 10             # total percentage (in % units) of the graph to be occupied by empty spaces, if any
ctrl.segmentcartlength = array([]) # array containing the cartesian length (in bohrs) of each segment of the band structure
ctrl.segmentrellength = array([])  # array containing the relative length (dimensionless) of each segment of the band structure

ctrl.dicoxkpt = {}        #dictionary linking x  to kpt (array)
ctrl.captiontick = list()
ctrl.X = array([])
ctrl.Y = list()           # list of array(nband,Float)
ctrl.nvalenceband = int()
ctrl.bandgap = float()
ctrl.hartree_to_eV = float(27.2113845) #eV/hartree

ctrl.ndataset = int() # number of datasets found in the .out file
ctrl.datasetlocation = list() # list containing the starting and ending line index of each dataset

ctrl.useddataset = list() # list containing the number of used dataset
ctrl.databasekey = list() # list of list containing the informations describing the calculation parameters used to generate the k-point database
ctrl.alphabet={
     1:'A', 2:'B', 3:'C', 4:'D', 5:'E', 6:'F', 7:'G', 8:'H', 9:'I',10:'J',11:'K',12:'L',13:'M',
    14:'N',15:'O',16:'P',17:'Q',18:'R',19:'S',20:'T',21:'U',22:'V',23:'W',24:'X',25:'Y',26:'Z'}

#data variables
data = VariableContainer()
data.nband = int()
data.kpt = list() # list of array(3,Float)
data.energy = list() # list of array(data.nband,Float)
data.specialkpt = list() # list of array(3,Float)
data.G = list() # list of array(3,Float) reciprocal vectors
data.fermienergy = float() #fermi energy, in eV 
data.typat = list() # list containing the type of each atom
data.units = str() # 'eV' or 'hartree'

#graph variables
graph = VariableContainer()
graph.title = str() # title of the graph
graph.worldymin = float()
graph.worldymax = float()

#feedback variables
feedback = VariableContainer()
feedback.feedback = False
feedback.filename = '.%s_feedback' %program[:-3]
feedback.ntimeused = int()
feedback.email = 'benjamin.tardif@umontreal.ca'
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#METHODS : general

def header(program,version):
    # type(program) = string
    # type(version) = string
    # returns a header to be printed in the shell each time the program is launched
    L = len(program)+len(version)+9+2
    line = L*'='
    header = '\n%s\n %s version %s\n%s' %(line,program,version,line)
    return header

def detectfile(filename,path='.'):
    # type(filename) = string
    # type(path) = string
    # returns True if the given file is found in the specified path
    if filename in os.listdir(path):
        return True
    else:
        return False

def floatable(x):
    # type(x) = string, int or float
    # returns True if given x can be converted to a float
    try:
        float(x)
        return True
    except:
        return False
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#METHODS : list manipulation

def clean(list):
    # type(list) = list of strings (usually obtained with the ".readlines()" method)
    # removes "\n" and "\r" and empty lines from given list
    L = len(list)
    for i in range(L):
        list[L-1-i] = list[L-1-i].replace('\n','')
        list[L-1-i] = list[L-1-i].replace('\r','')
        if list[L-1-i].split() == []:
            list.pop(L-1-i)

def clean2(list):
    # type(list) = list of strings (usually obtained with the ".readlines()" method)
    # removes "\n" and "\r" from given list and replaces empty lines by "#"
    L = len(list)
    for i in range(L):
        list[L-1-i] = list[L-1-i].replace('\n','')
        list[L-1-i] = list[L-1-i].replace('\r','')
        if list[L-1-i].split() == []:
            list[L-1-i] = "#"

def rmrepetitions(list,pairedlist=None):
    # type(list) = any list whith all elements of the same length
    # removes repeated entries in the list, keeping only the first occurence
    # (if a paired list is specified, data removed from list will be removes from pairedlist too)
    #     example : rmrepetition([1,2,2,3,4,4,4,8,7,6,7]) = [1,2,3,4,8,7,6]
    L = len(list)
    try:s = len(list[0])
    except:s = 1
    i = 0
    while i < len(list)-1:
        j = i+1
        while j < len(list):
            if sum(list[j] == list[i]) == s:
                list.pop(j)
                if pairedlist:pairedlist.pop(j)
                j-=1
            j+=1
        i+=1

def rmsuccessiverepetitions(list,pairedlist=None):
    # type(list) = any list whith all elements of the same length
    # removes repeated successives entries in the list, keeping only the first occurence.
    # (if a paired list is specified, data removed from list will be removes from pairedlist too)
    #     example : rmrepetition([1,2,2,3,4,4,4,1,2,3,3]) = [1,2,3,4,1,2,3]
    L = len(list)
    try:s = len(list[0])
    except:s = 1
    i = 0
    while i < len(list)-1:
        j = i+1
        if sum(list[j] == list[i]) == s:
            list.pop(j)
            if pairedlist:pairedlist.pop(j)
            i-=1
        i+=1
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#METHODS : vector operations
def norm(vector):
    # type(vector) = array(3,Float)
    # returns the norm of a vector
    x = vector[0]
    y = vector[1]
    z = vector[2]
    norm = (x**2+y**2+z**2)**0.5
    return norm

def angle(vector1,vector2):
    # type(vector1) = array(3,Float)
    # type(vector2) = array(3,Float)
    # returns the angle (in degree) between the two vectors
    arg = dot(vector1,vector2)/norm(vector1)/norm(vector2)
    if arg >  1:arg= 1
    if arg < -1:arg=-1
    theta = (arccos(arg))/pi*180
    return theta

def kpt_red_to_cart(kptred,primvectors):
    # type(kptred) = array(3,Float) representing the reduced coordinates of a k-point
    # type(primvectors) = a list of 3 array(3,Float) each representing a primitive vector in cartesian coordinates
    # returns an array(3,Float) containing the coordinates of the given k-point in cartesian coordinates
    kptcart = kptred[0]*primvectors[0] + kptred[1]*primvectors[1] + kptred[2]*primvectors[2]
    return kptcart
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#METHODS : setup system

def writesetupfile(setupfilename=default.setupfilename,path=sys.path[0],\
          launchcommand=default.launchcommand,autolaunch=default.autolaunch,\
          energyshift=default.energyshift,\
          valbandcolor=default.valbandcolor,conbandcolor=default.conbandcolor,\
          bandlinewidth=default.bandlinewidth,fermilinewidth=default.fermilinewidth,separatorlinewidth=default.separatorlinewidth,\
          emptyspacewidth=default.emptyspacewidth):

    writer = open(sys.path[0]+'/'+setupfilename,'w')

    writer.write('----------------------------------------------------------------------------')
    writer.write('\n this file is used by the program %s (version %s)' %(program,version))
    writer.write('\n (this file is not essential and can be deleted if needed)')
    writer.write('\n----------------------------------------------------------------------------')

    writer.write('\n\n============================================================================')
    writer.write('\n-setautolaunch')
    writer.write('\n\nXMGRACE LAUNCH COMMAND:\n%s' %launchcommand)
    writer.write('\n\nXMGRACE AUTOLAUNCH:\n%s' %autolaunch)
    writer.write('\n============================================================================')

    writer.write('\n\n============================================================================')
    writer.write('\n-setenergyshift')
    writer.write('\n\nSHIFT FERMI ENERGY TO ZERO:\n%s' %energyshift)
    writer.write('\n============================================================================')

    writer.write('\n\n============================================================================')
    writer.write('\n-setlinecolor')
    writer.write('\n\nVALENCE BANDS COLOR:\n%s' %valbandcolor)
    writer.write('\n\nCONDUCTION BANDS COLOR:\n%s' %conbandcolor)
    writer.write('\n============================================================================')

    writer.write('\n\n============================================================================')
    writer.write('\n-setlinewidth')
    writer.write('\n\nBAND LINES WIDTH:\n%s' %bandlinewidth)
    writer.write('\n\nFERMI ENERGY LINE WIDTH:\n%s' %fermilinewidth)
    writer.write('\n\nSEPARATOR LINES WIDTH:\n%s' %separatorlinewidth)
    writer.write('\n============================================================================')

    writer.write('\n\n============================================================================')
    writer.write('\n-setspacewidth')
    writer.write('\n\nEMPTY SPACE(S) WIDTH PERCENTAGE:\n%s' %emptyspacewidth)
    writer.write('\n============================================================================')

    writer.close()

def setupfilecompatibility(oldsetupfile,newsetupfile):
    reader = open(sys.path[0]+'/'+default.setupfilename,'r')
    setup1 = reader.readlines()
    reader.close()
    reader = open(sys.path[0]+'/'+default.setupfilename+'2','r')
    setup2 = reader.readlines()
    reader.close()
    if len(setup1) == len(setup2):
        i = 1 #skip the first three lines of the file
        while i < len(setup1)-1:
            i+=1
            if setup1[i] != setup2[i]: return False
            if ':' in setup1[i]: i+=1
        return True
    else:
        return False
        
def getsettings():
    # --> ctrl.launchcommand,ctrl.autolaunch,ctrl.energyshift,ctrl.valbandcolor,ctrl.conbandcolor,ctrl.bandlinewidth,ctrl.fermilinewidth,ctrl.emptyspacewidth
    reader = open(sys.path[0]+'/'+default.setupfilename,'r')
    setupfile = reader.readlines()
    reader.close()
    clean2(setupfile)
    for i in range(len(setupfile)):
        if setupfile[i] == 'XMGRACE LAUNCH COMMAND:':
            ctrl.launchcommand = setupfile[i+1]
        elif setupfile[i] == 'XMGRACE AUTOLAUNCH:':
            ctrl.autolaunch = setupfile[i+1]
        elif setupfile[i] == 'SHIFT FERMI ENERGY TO ZERO:':
            ctrl.energyshift = setupfile[i+1]
        elif setupfile[i] == 'VALENCE BANDS COLOR:':
            ctrl.valbandcolor = setupfile[i+1]
        elif setupfile[i] == 'CONDUCTION BANDS COLOR:':
            ctrl.conbandcolor = setupfile[i+1]
        elif setupfile[i] == 'BAND LINES WIDTH:':
            ctrl.bandlinewidth = setupfile[i+1]
        elif setupfile[i] == 'FERMI ENERGY LINE WIDTH:':
            ctrl.fermilinewidth = setupfile[i+1]
        elif setupfile[i] == 'SEPARATOR LINES WIDTH:':
            ctrl.separatorlinewidth = setupfile[i+1]
        elif setupfile[i] == 'EMPTY SPACE(S) WIDTH PERCENTAGE:':
            ctrl.emptyspacewidth = setupfile[i+1]
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#METHODS : feedback system

def feedbackbugged(jobtype='?'):
    if detectfile(feedback.filename,sys.path[0]) == False:
        #feedback file does not exists, create a default one
        writer = open(sys.path[0]+'/'+feedback.filename,'w')
        writer.write('times used / date used / version used / job type / job status / {crash reason}\n')
        writer.write('\n%s\t%s\t%s\t%s\t%s' %(feedback.ntimeused+1,time.ctime(),version,jobtype,'BUGGED'))
        writer.close()
    else:
        #feedback file already exist, update it
        reader = open(sys.path[0]+'/'+feedback.filename,'r')
        filedata = reader.readlines()
        reader.close()
        writer = open(sys.path[0]+'/'+feedback.filename,'w')
        if jobtype =='?':
            #this is the first time the method is being called
            feedback.ntimeused = int(filedata[-1].split()[0])
            writer.writelines(filedata)
            writer.write('\n%s\t%s\t%s\t%s\t%s' %(feedback.ntimeused+1,time.ctime(),version,jobtype,'BUGGED'))
            writer.close()
        else:
            #this is not the first time the method is being called
            filedata.pop(-1)
            writer.writelines(filedata)
            writer.write('%s\t%s\t%s\t%s\t%s' %(feedback.ntimeused+1,time.ctime(),version,jobtype,'BUGGED'))            

def feedbackcrashed(jobtype,reason):
    reader = open(sys.path[0]+'/'+feedback.filename,'r')
    filedata = reader.readlines()
    reader.close()
    filedata.pop(-1)
    writer = open(sys.path[0]+'/'+feedback.filename,'w')
    writer.writelines(filedata)
    writer.write('%s\t%s\t%s\t%s\t%s\t\t--> %s' %(feedback.ntimeused+1,time.ctime(),version,jobtype,'CRASHED',reason))
    writer.close()
    try:
        os.system('mail -s "%s feedback #%s" %s < %s' %(program,feedback.ntimeused+1,feedback.email,sys.path[0]+'/'+feedback.filename))
    except:
        pass

def feedbackcompleted(jobtype):
    reader = open(sys.path[0]+'/'+feedback.filename,'r')
    filedata = reader.readlines()
    reader.close()
    filedata.pop(-1)
    writer = open(sys.path[0]+'/'+feedback.filename,'w')
    writer.writelines(filedata)
    writer.write('%s\t%s\t%s\t%s\t%s' %(feedback.ntimeused+1,time.ctime(),version,jobtype,'COMPLETED'))
    writer.close()
    try:
        os.system('mail -s "%s feedback #%s" %s < %s' %(program,feedback.ntimeused+1,feedback.email,sys.path[0]+'/'+feedback.filename))
    except:
        pass
#=====================================================================================================================================================================


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#MAIN
print header(program,version)
if feedback.feedback==True:feedbackbugged()

#=====================================================================================================================================================================
#SETUP FILE

if detectfile(default.setupfilename,sys.path[0]) == False:
    #setup file not found, create a default one
    writesetupfile()
else:
    #a setup file already exists
    reader = open(sys.path[0]+'/'+default.setupfilename,'r')
    setupfile = reader.readlines()
    reader.close()
    clean2(setupfile)
    checkversion = setupfile[1].split()[-1].split(')')[0]
    #update the setup file if the checkversion is different from the version    
    if checkversion != version:
        print '\n- WARNING -\nnew version detected\n%s was upgraded from version %s to version %s' %(program,checkversion,version)
        writesetupfile(default.setupfilename+'2')
        if setupfilecompatibility(default.setupfilename,default.setupfilename+'2')==True:
            getsettings()
            writesetupfile(setupfilename=default.setupfilename,path=sys.path[0],\
                           launchcommand=ctrl.launchcommand,autolaunch=ctrl.autolaunch,\
                           energyshift=ctrl.energyshift,\
                           valbandcolor=ctrl.valbandcolor,conbandcolor=ctrl.conbandcolor,\
                           bandlinewidth=ctrl.bandlinewidth,fermilinewidth=ctrl.fermilinewidth,separatorlinewidth=ctrl.separatorlinewidth,\
                           emptyspacewidth=ctrl.emptyspacewidth)
        else:
            print '\n- WARNING -\nsetup file system has changed since your last version\nall settings restored to default\nyou may have to reset some of your previous settings'
            writesetupfile()
        os.system('rm -f %s/%s' %(sys.path[0],default.setupfilename+'2'))
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#COMMAND LINE

#abort if a keyword is not valid
for arg in ctrl.arglist:
    if arg[0] == '-':
        #a keyword is found
        if arg not in ctrl.validkeywords:
            #the keyword is not valid
            print '\n- ERROR -\n%s is not a valid keyword' %arg
            validkeywords = str()
            for keyword in ctrl.validkeywords:
                validkeywords+=keyword+', '
            validkeywords = validkeywords[:-2]
            print '\nvalid keywords are :\n%s\n' %validkeywords
            if feedback.feedback==True:feedbackcrashed('?','invalid keyword')
            sys.exit()

#abort if a keyword is repeated
for keyword in ctrl.validkeywords:
    if ctrl.arglist.count(keyword) > 1:
        print '\n- ERROR -\nkeyword %s is repeated %s times\n' %(keyword,ctrl.arglist.count(keyword))
        if feedback.feedback==True:feedbackcrashed('?','repeated keyword')
        sys.exit()

#get keywords
setautolaunch = False
setenergyshift = False
setlinecolor = False
setlinewidth = False
setspacewidth = False
setup= False
setdefault = False

if '-setautolaunch' in ctrl.arglist:
    ctrl.arglist.pop(ctrl.arglist.index('-setautolaunch'))
    setautolaunch = True

if '-setenergyshift' in ctrl.arglist:
    ctrl.arglist.pop(ctrl.arglist.index('-setenergyshift'))
    setenergyshift = True

if '-setlinecolor' in ctrl.arglist:
    ctrl.arglist.pop(ctrl.arglist.index('-setlinecolor'))
    setlinecolor = True

if '-setlinewidth' in ctrl.arglist:
    ctrl.arglist.pop(ctrl.arglist.index('-setlinewidth'))
    setlinewidth = True

if '-setspacewidth' in ctrl.arglist:
    ctrl.arglist.pop(ctrl.arglist.index('-setspacewidth'))
    setspacewidth = True

if '-setup' in ctrl.arglist:
    ctrl.arglist.pop(ctrl.arglist.index('-setup'))
    setautolaunch = True
    setenergyshift = True
    setlinecolor = True
    setlinewidth = True
    setspacewidth = True

if '-setdefault' in ctrl.arglist:
    ctrl.arglist.pop(ctrl.arglist.index('-setdefault'))
    setdefault = True

if '-debug' in ctrl.arglist:
    ctrl.arglist.pop(ctrl.arglist.index('-debug'))
    ctrl.debugmode = True

#(put additionnal keywords here)

#SETUP MODE
if setdefault==True:
    #create a default setup file
    print '\n--> starting SETUP MODE'
    if feedback.feedback==True:feedbackbugged('SETUP')
    writesetupfile()
    print '\ndefault setup restored'
    print '\n--> leaving SETUP MODE\n'
    if feedback.feedback == True:feedbackcompleted('SETUP')
    sys.exit()

getsettings()
if setautolaunch+setenergyshift+setlinecolor+setlinewidth+setspacewidth+setup!=0:
    print '\n--> starting SETUP MODE'
    if feedback.feedback==True:feedbackbugged('SETUP')

    if setautolaunch==True:
        #change launchcommand --> ctrl.launchcommand
        print '\ncurrent xmgrace launch command is :\n%s\n' %ctrl.launchcommand
        answer = str()
        while answer not in ['yes','no']:    
            answer = raw_input('do you wish to change it (yes ; no) ? ')
        if answer == 'yes':
            ctrl.launchcommand = raw_input('\nenter the new xmgrace launch command :\n')

        #change autolaunch --> ctrl.autolaunch
        ctrl.autolaunch = raw_input('\nautomatically launch xmgrace each time a .agr file is created (yes ; no) ? ')
        while ctrl.autolaunch not in ['yes','no']:
            ctrl.autolaunch = raw_input('automatically launch xmgrace each time a .agr file is created (yes ; no) ? ')

    if setenergyshift==True:
        #change energy shift --> ctrl.energyshift
        ctrl.energyshift = raw_input('\nshift energy eigeivalues to bring the fermi energy to zero (yes ; no) ? ')
        while ctrl.energyshift not in ['yes','no']:
            ctrl.energyshift = raw_input('shift energy eigeivalues to bring the fermi energy to zero (yes ; no) ? ')

    if setlinecolor==True:
        #change valence bands color --> ctrl.valbandcolor
        ctrl.valbandcolor = raw_input('\nChoose the color of the valence bands : ')
        while ctrl.valbandcolor not in ctrl.xmgracecolor.keys():
            colors = str()
            for color in ctrl.xmgracecolor.keys():
                colors += '%s, ' %color
            colors = colors[:-2]
            print '\n- invalid entry -\npossible answers are :\n%s' %colors
            ctrl.valbandcolor = raw_input('\nChoose the color of the valence bands : ')

        #change conduction bands color --> ctrl.conbandcolor
        ctrl.conbandcolor = raw_input('\nChoose the color of the conduction bands : ')
        while ctrl.conbandcolor not in ctrl.xmgracecolor.keys():
            colors = str()
            for color in ctrl.xmgracecolor.keys():
                colors += '%s, ' %color
            colors = colors[:-2]
            print '\n- invalid entry -\npossible answers are :\n%s' %colors
            ctrl.conbandcolor = raw_input('\nChoose the color of the conduction bands : ')

    if setlinewidth==True:
        #change band lines width --> ctrl.bandlinewidth
        ctrl.bandlinewidth = raw_input('\nChoose the width of the band lines : ')
        while floatable(ctrl.bandlinewidth) == False:
            ctrl.bandlinewidth = raw_input('Choose the width of the band lines : ')

        #change fermi energy line width --> ctrl.fermilinewidth
        ctrl.fermilinewidth = raw_input('\nChoose the width of the fermi energy line : ')
        while floatable(ctrl.fermilinewidth) == False:
            ctrl.fermilinewidth = raw_input('Choose the width of the fermi energy line : ')

        #change separator lines width --> ctrl.separatorlinewidth
        ctrl.separatorlinewidth = raw_input('\nChoose the width of the separator lines : ')
        while floatable(ctrl.separatorlinewidth) == False:
            ctrl.separatorlinewidth = raw_input('Choose the width of the separator lines : ')

    if setspacewidth==True:
        #change empty space(s) width --> ctrl.emptyspacewidth
        ctrl.emptyspacewidth = raw_input('\nChoose the total width (in percentage) of the empty space(s) on the graph, if any : ')
        while floatable(ctrl.emptyspacewidth) == False:
            ctrl.emptyspacewidth = raw_input('Choose the total width (in percentage) of the empty space(s) on the graph, if any : ')

    #overwrite setup file
    writesetupfile(setupfilename=default.setupfilename,path=sys.path[0],\
                   launchcommand=ctrl.launchcommand,autolaunch=ctrl.autolaunch,\
                   energyshift=ctrl.energyshift,\
                   valbandcolor=ctrl.valbandcolor,conbandcolor=ctrl.conbandcolor,\
                   bandlinewidth=ctrl.bandlinewidth,fermilinewidth=ctrl.fermilinewidth,separatorlinewidth=ctrl.separatorlinewidth,\
                   emptyspacewidth=ctrl.emptyspacewidth)

    print '\n--> leaving SETUP MODE\n'
    if feedback.feedback==True:feedbackcompleted('SETUP')
    sys.exit()

#get the filename --> ctrl.filename
if len(ctrl.arglist) == 0:
    #user entered no filename in the command line
    ctrl.filename = raw_input('\nEnter the filename : \n')
elif len(ctrl.arglist) == 1:
    #user entered the filename in the command line
    ctrl.filename = ctrl.arglist[0]
elif len(ctrl.arglist) > 1:
    #user entered too much arguments in the command line
    print '\n- ERROR -\ntoo many arguments entered in the command line\n'
    if feedback.feedback==True:feedbackcrashed('?','too many arguments entered')
    sys.exit()

#compute file type --> ctrl.filetype
#abort if the file type is not valid
if ctrl.filename.split('.')[-1][:3] == 'out':
    ctrl.filetype = 'out'
elif ctrl.filename.split('.')[-1][:3] == 'dbs':
    ctrl.filetype = 'dbs'
else:
    print '\n- ERROR -\ninvalid file type (must be .out or .dbs)\n'
    if feedback.feedback==True:feedbackcrashed('?','invalid filetype')
    sys.exit()

#abort if the file does not exists
if detectfile(ctrl.filename,'.') == False:
    print '\n- ERROR -\n"%s" file not found\n' %ctrl.filename
    if feedback.feedback==True:feedbackcrashed('?','file not found')
    sys.exit()

#activate debugmode, if needed
if ctrl.debugmode==True:
    print '\n--> DEBUG MODE'
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#READ THE FILE

#read file and acquire data --> ctrl.filedata
if ctrl.debugmode==True:print '\n--> file "%s" will be treated as a < %s > file' %(ctrl.filename,ctrl.filetype)
reader = open(ctrl.filename,"r")
ctrl.filedata = reader.readlines()
reader.close()
if ctrl.debugmode==True:print '\n--> file read successfully\n    %s line(s) read' %len(ctrl.filedata)
clean2(ctrl.filedata)
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#EXTRACT DATA FROM THE FILE

if ctrl.filetype == 'out':
    if feedback.feedback==True:feedbackbugged('OUT')

    #warning if the calculation is not completed
    calculationcompleted = False
    for line in ctrl.filedata:
        if line == ' Calculation completed.':
            calculationcompleted = True
    if calculationcompleted == False:
        print '\n- WARNING -\ncalculation not completed'

    #compute number of datasets --> ctrl.ndataset
    #compute first and last line of each dataset --> ctrl.datasetlocation
    for i in range(len(ctrl.filedata)):
        if ctrl.filedata[i][:10] == '== DATASET':
            ctrl.ndataset += 1
            if ctrl.datasetlocation != list():
                ctrl.datasetlocation[-1][4] = i
            ctrl.datasetlocation.append(['DATASET %s' %ctrl.ndataset,'first line =',i+1,'last line =','not found'])
        if ctrl.filedata[i][:17] == '== END DATASET(S)':
            ctrl.datasetlocation[-1][4] = i
    if ctrl.debugmode==True:
        print '\n--> dataset locations computed'
        for line in ctrl.datasetlocation:
            print '    %s' %line

    #compute list of datasets to use --> ctrl.useddataset
    validanswer = False
    allowedchars = ['-',',','0','1','2','3','4','5','6','7','8','9']
    print '\n%s dataset(s) detected' %ctrl.ndataset
    if ctrl.datasetlocation[-1][-1] == 'not found':
        print 'the last dataset is not completed and will be ignored'
        ctrl.ndataset -= 1

    if ctrl.ndataset == 0:
        print '\n- ERROR -\nno completed dataset available\n'
        if feedback.feedback==True:feedbackcrashed('OUT','no completed dataset')
        sys.exit()
    elif ctrl.ndataset == 1:
        ctrl.useddataset = [1]
    elif ctrl.ndataset > 1:
        while validanswer == False:
            answer = raw_input('\nWhich dataset(s) do you want to use (1 to %s) ? ' %ctrl.ndataset)

            ctrl.useddataset = list()
            validanswer = True
            
            #removes empty spaces from answer
            answersplit = answer.split()
            answer = str()
            for splitted in answersplit:
                answer+=splitted

            #compute ctrl.useddataset
            try:
                S = answer.split(',')
                for i in range(len(S)):
                    if '-' in S[i]:
                        a = int(S[i].split('-')[0])
                        b = int(S[i].split('-')[1])
                        ctrl.useddataset += range(a,b+1)
                    else:
                        ctrl.useddataset += [int(S[i])]
                rmrepetitions(ctrl.useddataset)
                ctrl.useddataset = sort(ctrl.useddataset)
            except:
                validanswer = False

            #verify validity
            for number in ctrl.useddataset:
                if number < 1 or number > ctrl.ndataset:
                    validanswer = False

            #show format instructions to user if invalid entry
            if validanswer == False:
                    print '\n- invalid entry -'
                    print 'use commas to separate different datasets'
                    print 'you can use minus signs to specify a group of successive datasets'
                    print 'for example, if you want to use the datasets 1, 3, 4, 5, 6 and 8, type : 1,3-6,8'

    if ctrl.debugmode==True:print '\n--> list of used datasets computed\n    %s' %ctrl.useddataset

    #get type of each atom --> data.typat
    #(assuming only one occurence of "typat" is present in the .out file)
    try:
        flag_typat = False
        k=0
        for i in range(len(ctrl.filedata)):
            if flag_typat==False and ctrl.filedata[i].split()[0] == 'typat':
                flag_typat = True
                k = i

        data.typat = ctrl.filedata[k].split()[1:]
        while ctrl.filedata[k+1].split()[0].isdigit()==True:
            k+=1
            for j in range(len(ctrl.filedata[k].split())):
                data.typat.append(ctrl.filedata[k].split()[j])

        for i in range(len(data.typat)):
            data.typat[i] = int(data.typat[i])
    except:
        data.typat = '?'

    if ctrl.debugmode==True:print '\n--> typat found\n    %s' %data.typat

    #compute number of valence bands --> ctrl.nvalenceband
    #(assuming only one occurence of "- pspini" is present for each atom type)
    try:
        nion = list()
        for i in range(len(ctrl.filedata)):
            if ctrl.filedata[i][:9] == '- pspini:':
                nion.append(float(ctrl.filedata[i+3].split()[2]))

        for i in range(len(data.typat)):
            ctrl.nvalenceband += int(nion[data.typat[i]-1])

        ctrl.nvalenceband = ctrl.nvalenceband/2.0
        if ctrl.nvalenceband%1 == 0:
            # ctrl.nvalenceband is an integer
            ctrl.nvalenceband = int(ctrl.nvalenceband)
        else:
            # ctrl.nvalenceband is not an integer
            ctrl.nvalenceband = int(ctrl.nvalenceband) + 1
    except:
        ctrl.nvalenceband = '?'

    if ctrl.debugmode==True:print '\n--> number of valence bands computed\n    %s' %ctrl.nvalenceband

    #get fermi energy --> data.fermienergy
    #(assuming only one occurence of "Fermi energy" is present in the .out file)
    try:
        for i in range(len(ctrl.filedata)):
            if ctrl.filedata[i].split()[0] == 'Fermi':
                data.fermienergy = float(ctrl.filedata[i].split()[-5])*hartree_to_eV
    except:
        pass

    if data.fermienergy == float(0):
        data.fermienergy = 'automatic'
        
    if ctrl.debugmode==True:print '\n--> fermi energy found\n    %s' %data.fermienergy

#------------------------------
    parser_template = '{}\d?\s?=?\s*([\d.E+]+)'
    #compute k-points and energy eigenvalues for each dataset
    starter = ctrl.filedata[:ctrl.datasetlocation[0][2]-1]
    for d in range(len(ctrl.useddataset)):
        n = ctrl.useddataset[d] #number of the dataset
        dataset = ctrl.filedata[ctrl.datasetlocation[n-1][2]-1:ctrl.datasetlocation[n-1][4]]

        #compute the dataset key --> datasetkey
        flag_nband = False
        datasetkey = list()

        #ecut
        datasetkey.append(['ecut:'])
        for i in range(len(starter)):
            if starter[i].split()[0] == 'ecut%s' %n or starter[i].split()[0] == 'ecut':
                value = float(re.search(parser_template.format('ecut'), starter[i]).group(1))
                datasetkey[-1].append(value)
        if len(datasetkey[-1]) == 1:
            datasetkey[-1].append('notfound')

        #natom
        datasetkey.append(['natom:'])
        for i in range(len(starter)):
            if starter[i].split()[0] == 'natom%s' %n or starter[i].split()[0] == 'natom':
                value = float(re.search(parser_template.format('natom'), starter[i]).group(1))
                datasetkey[-1].append(value)
        if len(datasetkey[-1]) == 1:
            datasetkey[-1].append(float(1)) #default

        #nband
        datasetkey.append(['nband:'])
        for i in range(len(starter)):
            if starter[i].split()[0] == 'nband%s' %n or starter[i].split()[0] == 'nband':
                value = float(re.search(parser_template.format('nband'), starter[i]).group(1))
                datasetkey[-1].append(value)
        if len(datasetkey[-1]) == 1:
            datasetkey[-1].append(float(1)) #default

        #occopt
        datasetkey.append(['occopt:'])
        for i in range(len(starter)):
            if starter[i].split()[0] == 'occopt%s' %n or starter[i].split()[0] == 'occopt':
                value = float(re.search(parser_template.format('occopt'), starter[i]).group(1))
                datasetkey[-1].append(value)
        if len(datasetkey[-1]) == 1:
            datasetkey[-1].append(float(1)) #default

        #set fermi energy to "automatic" if occopt is non metallic 
        if datasetkey[-1][-1] in [0,1,2]:
            data.fermienergy = 'automatic'

        #toldfe
        datasetkey.append(['toldfe:'])
        for i in range(len(starter)):
            if starter[i].split()[0] == 'toldfe%s' %n or starter[i].split()[0] == 'toldfe':
                datasetkey[-1].append(float(starter[i].split()[1]))
        if len(datasetkey[-1]) == 1:
            datasetkey[-1].append(float(0)) #default

        #toldff
        datasetkey.append(['toldff:'])
        for i in range(len(starter)):
            if starter[i].split()[0] == 'toldff%s' %n or starter[i].split()[0] == 'toldff':
                datasetkey[-1].append(float(starter[i].split()[1]))
        if len(datasetkey[-1]) == 1:
            datasetkey[-1].append(float(0)) #default

        #tolvrs
        datasetkey.append(['tolvrs:'])
        for i in range(len(starter)):
            if starter[i].split()[0] == 'tolvrs%s' %n or starter[i].split()[0] == 'tolvrs':
                datasetkey[-1].append(float(starter[i].split()[1]))
        if len(datasetkey[-1]) == 1:
            datasetkey[-1].append(float(0)) #default

        #tolwfr
        datasetkey.append(['tolwfr:'])
        for i in range(len(starter)):
            if starter[i].split()[0] == 'tolwfr%s' %n or starter[i].split()[0] == 'tolwfr':
                datasetkey[-1].append(float(starter[i].split()[1]))
        if len(datasetkey[-1]) == 1:
            datasetkey[-1].append(float(0)) #default

        #typat
        datasetkey.append(['typat:'])
        for i in range(len(starter)):
            if starter[i].split()[0] == 'typat%s' %n or starter[i].split()[0] == 'typat':
                k = i
                temp = list()
                temp = starter[k].split()[1:]
                while starter[k+1].split()[0].isdigit()==True:
                    k+=1
                    temp += starter[k].split()
                for j in range(len(temp)):
                    temp[j] = int(temp[j])
                datasetkey[-1]+=temp
        if len(datasetkey[-1]) == 1:
            datasetkey[-1].append(float(1)) #default

        #reciprocalvectors
        datasetkey.append(['reciprocalvectors:'])
        for i in range(len(dataset)):
            if dataset[i].split()[0] == 'R(1)=':
                for j in range(3):
                    linesplit = dataset[i+j].split()
                    datasetkey[-1].append(float(linesplit[-3]))
                    datasetkey[-1].append(float(linesplit[-2]))
                    datasetkey[-1].append(float(linesplit[-1]))
        if len(datasetkey[-1]) == 1:
            datasetkey[-1].append('notfound') #default

        #reduced coordinates
        datasetkey.append(['reducedcoordinates:'])
        for i in range(len(dataset)):
            if dataset[i][:20] == ' reduced coordinates':
                k = i
                while len(dataset[k+1].split()) == 3:
                    datasetkey[-1].append(float(dataset[k+1].split()[0]))
                    datasetkey[-1].append(float(dataset[k+1].split()[1]))
                    datasetkey[-1].append(float(dataset[k+1].split()[2]))
                    k+=1

        #verify the dataset key
        if d == 0:
            #compute database key --> ctrl.databasekey
            ctrl.databasekey = datasetkey
            refdataset = n
        else:
            if datasetkey != ctrl.databasekey:
                print '\n- ERROR -\nDATASET %s is not compatible with DATASET %s' %(n,refdataset)

                #given reason
                for i in range(len(datasetkey)):
                    if datasetkey[i] != ctrl.databasekey[i]:
                        print '%s are different' %datasetkey[i][0][:-1]
                print ''
                if feedback.feedback==True:feedbackcrashed('OUT','datasetkey not compatible')
                sys.exit()
            else:
                pass

        #get eigenvalue energy units
        for line in dataset:
            if line.split()[0] == 'Eigenvalues':
                data.units = line.replace('(','').replace(')','').split()[1]
        
        #get k-points --> data.kpt
        #get energy eigenvalues --> data.energy
        kptlist = list()
        for i in range(len(dataset)):
            if dataset[i].split()[0][:4] == 'kpt#':
                linesplit = dataset[i].split()
                kptlist.append(array([float(linesplit[-5]),float(linesplit[-4]),float(linesplit[-3])],dtype=float))
                k=i+1
                energylist = list()
                while dataset[k].split()[0].replace('-','').replace('.','').replace('e','').replace('E','').isdigit():
                    linesplit = dataset[k].split()
                    for j in range(len(linesplit)):
                        energylist.append(float(linesplit[j]))
                    k+=1
                if data.units == 'hartree':
                    data.energy.append(array(energylist)*ctrl.hartree_to_eV)
                else:
                    data.energy.append(array(energylist))
        data.kpt += kptlist
        if ctrl.debugmode==True:
            print '\n--> k-points found for DATASET %s\n    {%s element(s)}' %(n,len(kptlist))
        if ctrl.debugmode==True:
            print '\n--> energy eigenvalues found for DATASET %s\n    {%s element(s)} for each k-point' %(n,len(data.energy[0]))

#-------------------------------------



    #compute special k-points --> data.specialkpt
    rmsuccessiverepetitions(data.kpt,data.energy)
    for i in range(len(data.kpt)):
        if i == 0:
            data.specialkpt.append(data.kpt[i])
        elif i == 1:
            vector2 = data.kpt[i]-data.kpt[i-1]
        elif i == len(data.kpt)-1:
            data.specialkpt.append(data.kpt[i])
        else:
            vector1 = data.kpt[i] - data.kpt[i-1]
            
            if angle(vector1,vector2) < ctrl.angletol:
                pass
            else:
                data.specialkpt.append(data.kpt[i-1])
            vector2 = vector1
    if ctrl.debugmode==True:
        print '\n--> special k-points computed\n    %s element(s)' %len(data.specialkpt)
        for i in range(len(data.specialkpt)):
            print '    %s' %data.specialkpt[i]

    #compute band structure scheme --> ctrl.bandstructurescheme 
    L = 0
    dico = dict()
    for i in range(len(data.specialkpt)):
        k = str(data.specialkpt[i])
        if not k in dico.keys():
            L+=1
            dico[k] = '{l%s}' %ctrl.alphabet[L]
        ctrl.bandstructurescheme += '-%s' %dico[k]
    ctrl.bandstructurescheme = ctrl.bandstructurescheme[1:]
    if ctrl.debugmode==True:print '\n--> band structure scheme computed\n    %s' %ctrl.bandstructurescheme

if ctrl.filetype == 'dbs':
    if feedback.feedback==True:feedbackbugged('DBS')

    #get graph title --> graph.title
    for i in range (len(ctrl.filedata)):
        if ctrl.filedata[i].split()[0] == 'GRAPH':
            graph.title = ctrl.filedata[i+1]
            while graph.title[0] == ' ':graph.title = graph.title[1:]
            while graph.title[-1] == ' ':graph.title = graph.title[:-1]

    #get reciprocal vectors --> data.G
    for i in range(len(ctrl.filedata)):
        if ctrl.filedata[i].split()[0] == 'reciprocalvectors:':
            linesplit = ctrl.filedata[i].split()
            data.G.append(array([float(linesplit[1]),float(linesplit[2]),float(linesplit[3])],dtype=float))
            data.G.append(array([float(linesplit[4]),float(linesplit[5]),float(linesplit[6])],dtype=float))
            data.G.append(array([float(linesplit[7]),float(linesplit[8]),float(linesplit[9])],dtype=float))
    if ctrl.debugmode==True:
        print '\n--> reciprocal vectors found'
        for i in range(3):
            print '    G(%s)= %s' %(i+1,data.G[i])

    #get number of valence bands
    for i in range(len(ctrl.filedata)):
        if ctrl.filedata[i].split()[0] == 'NUMBER':
            ctrl.nvalenceband = int(ctrl.filedata[i+1].split()[0])

    #get fermi energy
    for i in range(len(ctrl.filedata)):
        if ctrl.filedata[i].split()[0] == 'FERMI':
            data.fermienergy = ctrl.filedata[i+1].split()[0]
    if ctrl.debugmode==True:print '\n--> fermi energy found\n    %s' %data.fermienergy

    #get special k-points --> ctrl.dicospecialkpt {caption:array}
    for i in range(len(ctrl.filedata)):
        linesplit = ctrl.filedata[i].split()
        if '}=' in linesplit[0] and linesplit[0][0] == '{':
            kptcaption = linesplit[0].split('=')[0]
            ctrl.dicospecialkpt[kptcaption] = array([float(linesplit[1]),float(linesplit[2]),float(linesplit[3])],dtype=float)
    if ctrl.debugmode==True:
        print '\n--> special k-points found\n    %s element(s)' %len(ctrl.dicospecialkpt)
        for i in range(len(ctrl.dicospecialkpt)):
            print '    %s : %s' %(ctrl.dicospecialkpt.keys()[i],ctrl.dicospecialkpt.values()[i])

    #get band structure scheme --> ctrl.bandstructurescheme
    for i in range(len(ctrl.filedata)):
        if ctrl.filedata[i].split()[0][0] == '{' and ctrl.filedata[i].split()[-1][-1] == '}':
            ctrl.bandstructurescheme = ctrl.filedata[i]
            while ctrl.bandstructurescheme[0]  != '{':ctrl.bandstructurescheme = ctrl.bandstructurescheme[1:]
            while ctrl.bandstructurescheme[-1] != '}':ctrl.bandstructurescheme = ctrl.bandstructurescheme[:-1]
    if ctrl.debugmode==True:print '\n--> band structure scheme found\n    %s' %ctrl.bandstructurescheme

    #get k-points --> data.kpt
    #get energy eigenvalues --> data.energy
    for i in range(len(ctrl.filedata)):
        if ctrl.filedata[i].split()[0] == 'kpt':
            linesplit = ctrl.filedata[i+1].split()
            data.kpt.append(array([float(linesplit[0]),float(linesplit[1]),float(linesplit[2])],dtype=float))
            energieslist = list()
            for energy in ctrl.filedata[i+2].split():
                energieslist.append(float(energy))
            data.energy.append(array(energieslist))
    if ctrl.debugmode==True:print '\n--> k-points found\n    {%s element(s)}' %len(data.kpt)
    if ctrl.debugmode==True:print '\n--> energy eigenvalues found\n    {%s element(s)} for each k-point' %len(data.energy[0])

    #compute segment caption list --> ctrl.segmentcaptionlist
    for captiongroup in ctrl.bandstructurescheme.split():
        captions = captiongroup.split('-')
        for i in range(len(captions)-1):
            ctrl.segmentcaptionlist.append([captions[i],captions[i+1]])
        ctrl.segmentcaptionlist.append(['empty space'])
    ctrl.segmentcaptionlist.pop(-1)
    if ctrl.debugmode==True:
        print '\n--> segment caption list computed\n    %s element(s)' %len(ctrl.segmentcaptionlist)
        for i in range(len(ctrl.segmentcaptionlist)):
            print '    %s' %ctrl.segmentcaptionlist[i]

    #compute segment cartesian length --> ctrl.segmentcartlength
    nvac = 0
    nseg = 0
    totallen = 0
    segmentcartlength = list()
    for caption in ctrl.segmentcaptionlist:
        if caption[0] == 'empty space':
            nvac+=1
            segmentcartlength.append('empty space')
        else:
            nseg+=1
            ki = kpt_red_to_cart(ctrl.dicospecialkpt[caption[0]],data.G)
            kf = kpt_red_to_cart(ctrl.dicospecialkpt[caption[1]],data.G)
            segmentcartlength.append(norm(kf-ki))
            totallen += segmentcartlength[-1]
    
    if nvac != 0:
        spacelen = (float(ctrl.spacepercent)/100)*totallen/nvac/(1-float(ctrl.spacepercent)/100)
        for i in range(len(segmentcartlength)):
            if segmentcartlength[i] == 'empty space':
                segmentcartlength[i]=spacelen

    ctrl.segmentcartlength = array(segmentcartlength)
    if ctrl.debugmode==True:
        print '\n--> segment cartesian length computed\n    %s element(s)' %len(ctrl.segmentcartlength)
        for i in range(len(ctrl.segmentcartlength)):
            print '    %s' %ctrl.segmentcartlength[i]

    #compute segment relative length --> ctrl.segmentrellength
    totallen = sum(ctrl.segmentcartlength)
    segmentrellength = list()
    for length in ctrl.segmentcartlength:
        segmentrellength.append(length/totallen)
    ctrl.segmentrellength = array(segmentrellength)
    if ctrl.debugmode==True:
        print '\n--> segment relative length computed\n    %s element(s)' %len(ctrl.segmentrellength)
        for i in range(len(ctrl.segmentrellength)):
            print '    %s' %ctrl.segmentrellength[i]
    
    #compute positions of xticks --> ctrl.xtick 
    xtick = list()
    for i in range(len(ctrl.segmentrellength)):
        xtick.append(sum(ctrl.segmentrellength[:i+1]))
    xtick.insert(0,float(0))
    ctrl.xtick = array(xtick)
    if ctrl.debugmode==True:
        print '\n--> positions of xticks computed\n    %s element(s)' %len(ctrl.xtick)
        for i in range(len(ctrl.xtick)):
            print '    %s' %ctrl.xtick[i]


    #compute captions of xticks --> ctrl.captiontick
    ctrl.captiontick = ctrl.bandstructurescheme.replace('-',' ').split()
    for i in range(len(ctrl.captiontick)):
        if ctrl.captiontick[i][1] == 'l':
            ctrl.captiontick[i] = '"%s"' %ctrl.captiontick[i][2]
        elif ctrl.captiontick[i][1] == 's':
            ctrl.captiontick[i] = '"\\f{Symbol}%s"' %ctrl.captiontick[i][2]


    if ctrl.debugmode==True:
        print '\n--> captions of xticks computed\n    %s elements(s)' %len(ctrl.captiontick)
        for i in range(len(ctrl.captiontick)):
            print '    %s' %ctrl.captiontick[i]
    

    #compute dictionary mapping x coordinate on graph with k-points --> ctrl.dicoxkpt {x,kpt} (type(x) = float, type(kpt)=array(3,Float))
    for i in range(len(ctrl.segmentcaptionlist)):
        caption = ctrl.segmentcaptionlist[i]
        if caption == ['empty space']:
            pass
        else:
            correct = 0
            ki = ctrl.dicospecialkpt[caption[0]]
            kf = ctrl.dicospecialkpt[caption[1]]
            for k in range(len(data.kpt)):
                kpt=data.kpt[k]
                goodkpt=False
                if list(kpt)==list(ki):
                    goodkpt = True
                    xfrac = 0
                elif list(kpt)==list(kf):
                    goodkpt = True
                    xfrac=1
                elif angle(kpt-ki,kf-ki) < ctrl.angletol2 and dot(kf-kpt,kf-ki) > 0.0:
                    goodkpt = True
                    xfrac=dot(kpt-ki,kf-ki)/norm(kf-ki)/norm(kf-ki)
                if goodkpt == True:
                    correct+=1
                    ctrl.dicoxkpt[xfrac*ctrl.segmentrellength[i]+ctrl.xtick[i]] = kpt

    #compute abcissa array --> ctrl.X
    ctrl.X = sort(ctrl.dicoxkpt.keys())
    
    #compute ordinate arrays --> ctrl.Y
    xsort = sort(ctrl.dicoxkpt.keys())
    for i in range(len(xsort)):
        x = xsort[i]
        k = ctrl.dicoxkpt[x]
        for j in range(len(data.kpt)):
            if list(data.kpt[j]) == list(k):
                index = j
        ctrl.Y.append(data.energy[index])
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#WRITE THE FILE

if ctrl.filetype == 'out':

    #compute the .dbs filename --> ctrl.dbsfilename
    ctrl.dbsfilename = '%s.dbs' %ctrl.filename

    #open the writer
    writer = open(ctrl.dbsfilename,"w")

    #write the default graph title
    writer.write('GRAPH TITLE:\nBand Structure from %s\n' %ctrl.filename)

    #write number of valence bands
    writer.write('\nNUMBER OF VALENCE BANDS:\n%s\n' %ctrl.nvalenceband)

    #write the fermi energy
    writer.write('\nFERMI ENERGY (eV):\n%s\n' %data.fermienergy)

    #write the special kpts
    rmrepetitions(data.specialkpt)
    writer.write('\nSPECIAL K-POINTS (reduced coord):\n')
    for i in range(len(data.specialkpt)):
        k = data.specialkpt[i]
        kx = '%1.4f' %k[0]
        ky = '%1.4f' %k[1]
        kz = '%1.4f' %k[2]
        if k[0]>=0:kx = ' %1.4f' %k[0]
        if k[1]>=0:ky = ' %1.4f' %k[1]
        if k[2]>=0:kz = ' %1.4f' %k[2]
        if i <= 25:
            writer.write('{l%s}= %s %s %s\n' %(ctrl.alphabet[i+1],kx,ky,kz))
        else:
            writer.write('{l%s}= %s %s %s\n' %('x',kx,ky,kz))

    #write the band structure scheme
    writer.write('\nBAND STRUCTURE SCHEME:\n')
    writer.write('%s\n' %ctrl.bandstructurescheme)

    #write the kpts and energies
    rmrepetitions(data.kpt,data.energy)
    writer.write('\n\n\nDATABASE:\n')
    for i in range(len(data.kpt)):
        k = data.kpt[i]
        kx = '%1.4f' %k[0]
        ky = '%1.4f' %k[1]
        kz = '%1.4f' %k[2]
        if k[0]>=0:kx = ' %1.4f' %k[0]
        if k[1]>=0:ky = ' %1.4f' %k[1]
        if k[2]>=0:kz = ' %1.4f' %k[2]
        writer.write('kpt %s\n%s %s %s\n' %(i+1,kx,ky,kz))
        for e in data.energy[i]:
            if e >=0:writer.write(' %e ' %e)
            else:writer.write('%e ' %e)
        writer.write('\n')

    #write the database key
    writer.write('\nDATABASE KEY:\n')
    for i in range(len(ctrl.databasekey)):
        for j in range(len(ctrl.databasekey[i])):
            writer.write('%s ' %ctrl.databasekey[i][j])
        writer.write('\n')

    #close the writer
    writer.close()

    print '\n"%s" file created successfully\n' %ctrl.dbsfilename
    if feedback.feedback==True:feedbackcompleted('OUT')


if ctrl.filetype == 'dbs':
    
    #compute the .agr filename --> ctrl.agrfilename
    ctrl.agrfilename = '%s.agr' %ctrl.filename[:-4]

    #compute fermi energy value in eV --> ctrl.fermienergy
    if data.fermienergy == 'automatic':
        maxlist = list()
        for i in range(len(ctrl.Y)):
            maxlist.append(ctrl.Y[i][ctrl.nvalenceband-1])
        data.fermienergy = max(maxlist)
    else:
        data.fermienergy = float(data.fermienergy)

    #compute the energy gap --> ctrl.bandgap
    maxHOMOlist = list()
    minLUMOlist = list()
    for i in range(len(ctrl.Y)):
        maxHOMOlist.append(ctrl.Y[i][ctrl.nvalenceband-1])
        minLUMOlist.append(ctrl.Y[i][ctrl.nvalenceband])
    maxHOMO = max(maxHOMOlist)
    minLUMO = min(minLUMOlist)
    ctrl.bandgap = minLUMO - maxHOMO

    #compute gap style --> ctrl.gapstyle
    ctrl.gapstyle = '(indirect)' #default
    maxHOMOindexlist = list()
    minLUMOindexlist = list()

    i = 0
    while maxHOMO in maxHOMOlist:
        maxHOMOindexlist.append(maxHOMOlist.index(maxHOMO)+i)
        maxHOMOlist.pop(maxHOMOlist.index(maxHOMO))
        i += 1

    i = 0
    while minLUMO in minLUMOlist:
        minLUMOindexlist.append(minLUMOlist.index(minLUMO)+i)
        minLUMOlist.pop(minLUMOlist.index(minLUMO))
        i += 1

    for M in maxHOMOindexlist:
        if M in minLUMOindexlist:
            ctrl.gapstyle = '(direct)'

    #shift energies to bring the fermi energy to zero, if wanted
    if ctrl.energyshift == 'yes':
        for i in range(len(ctrl.Y)):
            ctrl.Y[i] = ctrl.Y[i] - data.fermienergy
        data.fermienergy = 0.0

    #compute plot.worldymin --> graph.worldymin
    minlist = list()
    for array in ctrl.Y:
        minlist.append(min(array))
    graph.worldymin = min(minlist)
    
    #compute plot.worldymax --> graph.worldymax
    maxlist = list()
    for array in ctrl.Y:
        maxlist.append(max(array))
    graph.worldymax = max(maxlist)

    #adjust worldymin et worldymax
    width = graph.worldymax - graph.worldymin
    graph.worldymin = graph.worldymin - 0.05*width
    graph.worldymax = graph.worldymax + 0.05*width

    #open the writer
    writer = open(ctrl.agrfilename,"w")

    #write the file
    writer.write('# file produced using %s version %s\n' %(program,version))
    writer.write('\n')
    writer.write('# version of xmgrace:\n')
    writer.write('@   version 50114\n')
    writer.write('\n')

    writer.write('# graph title:\n')
    writer.write('@   title "%s"\n' %graph.title)
    writer.write('\n')

    writer.write('# graph range:\n')
    writer.write('@   world xmin 0\n')
    writer.write('@   world xmax 1\n')
    writer.write('@   world ymin %s\n' %graph.worldymin)
    writer.write('@   world ymax %s\n' %graph.worldymax)
    writer.write('\n')

    writer.write('# X axis properties:\n')
    writer.write('@   xaxis tick major size 0.0\n') #height of x tick lines
    writer.write('@   xaxis tick spec type both\n') #???
    writer.write('@   xaxis tick spec 15\n') #???
    for i in range(len(ctrl.xtick)):
        writer.write('@   xaxis tick major %s, %s\n' %(i,ctrl.xtick[i]))
        writer.write('@   xaxis ticklabel %s, %s\n' %(i,ctrl.captiontick[i]))
    if ctrl.bandgap > 0:
         writer.write('@   xaxis label "E\sG %s\N = %s eV"\n' %(ctrl.gapstyle,ctrl.bandgap))
    writer.write('\n')

    writer.write('# Y axis properties:\n')
    writer.write('@   yaxis label "Energy (eV)"\n')
    writer.write('@   yaxis tick major 5\n')
    writer.write('@   yaxis tick minor 1\n')
    writer.write('@   yaxis tick place normal\n')
    writer.write('\n')

    writer.write('# alternate Y axis properties:\n')
    writer.write('@   altyaxis on\n')
    writer.write('@   altyaxis ticklabel on\n')
    writer.write('@   altyaxis ticklabel place opposite\n')
    writer.write('@   altyaxis ticklabel type spec\n')
    writer.write('@   altyaxis tick type spec\n')
    writer.write('@   altyaxis tick spec 1\n')
    writer.write('@   altyaxis tick major 0, %s\n' %data.fermienergy)
    writer.write('@   altyaxis ticklabel 0, "\\f{Symbol}e\\f{}\sF\N"\n') #epsilon fermi symbol
    writer.write('\n')

    writer.write('# frame properties:\n')
    writer.write('@   frame linewidth %s\n' %ctrl.separatorlinewidth)
    writer.write('\n')
    
    s = 0

    writer.write('# plot of energy bands:\n')
    for i in range(len(ctrl.Y[0])):
        if i+1 > ctrl.nvalenceband:
            color = ctrl.xmgracecolor[ctrl.conbandcolor] 
        else:
            color = ctrl.xmgracecolor[ctrl.valbandcolor]
        writer.write('@   s%s line linewidth %s\n' %(s,ctrl.bandlinewidth))
        writer.write('@   s%s line color %s\n' %(s,color))
        s+=1
        writer.write('@   TYPE xy\n')
        for j in range(len(ctrl.X)):
            writer.write('      %s \t %s\n' %(ctrl.X[j],ctrl.Y[j][i]))
        writer.write('    &\n')
    writer.write('\n')

    writer.write('# plot of fermi energy line:\n')
    writer.write('@   s%s linewidth %s\n' %(s,ctrl.fermilinewidth))
    writer.write('@   s%s linestyle 2\n' %s)
    writer.write('@   s%s line color 1\n' %s)
    s+=1
    writer.write('@   TYPE xy\n')
    writer.write('      %s \t %s\n' %(0,data.fermienergy))
    writer.write('      %s \t %s\n' %(1,data.fermienergy))
    writer.write('    &\n')
    writer.write('\n')

    writer.write('# plot of empty spaces:\n')
    for i in range(len(ctrl.segmentcaptionlist)):
        if ctrl.segmentcaptionlist[i] == ['empty space']:
            writer.write('@   s%s linewidth %s\n' %(s,ctrl.bandlinewidth))
            writer.write('@   s%s line color 0\n' %s)
            s+=1
            xi = ctrl.xtick[i]
            xf = ctrl.xtick[i+1]
            xsort = list(sort(ctrl.dicoxkpt.keys()))
            index = xsort.index(xi)
            writer.write('@   TYPE xy\n')
            writer.write('      %s \t %s\n' %(xi,data.fermienergy))
            writer.write('      %s \t %s\n' %(xf,data.fermienergy))
            writer.write('    &\n')
            for j in range(len(ctrl.Y[0])):
                writer.write('@   s%s linewidth %s\n' %(s,ctrl.bandlinewidth))
                writer.write('@   s%s line color 0\n' %s)
                s+=1
                writer.write('@   TYPE xy\n')
                writer.write('      %s \t %s\n' %(xi,ctrl.Y[index][j]))
                writer.write('      %s \t %s\n' %(xf,ctrl.Y[index+1][j]))
            writer.write('    &\n')
    writer.write('\n')

    writer.write('# plot of vertical separators:\n')
    for i in range(len(ctrl.xtick)-2):
        writer.write('@   s%s linewidth %s\n' %(s,ctrl.separatorlinewidth))
        writer.write('@   s%s line color 1\n' %s)
        s+=1
        writer.write('@   TYPE xy\n')
        writer.write('      %s \t %s\n' %(ctrl.xtick[i+1],graph.worldymin))
        writer.write('      %s \t %s\n' %(ctrl.xtick[i+1],graph.worldymax))
        writer.write('    &\n')
    writer.write('\n')
        
    #close the writer
    writer.close()

    print '\n"%s" file created successfully\n' %ctrl.agrfilename
    if feedback.feedback==True:feedbackcompleted('DBS')
    
#=====================================================================================================================================================================


#=====================================================================================================================================================================
#AUTOLAUNCH

if ctrl.filetype == 'dbs' and ctrl.autolaunch == 'yes':
    print 'launching xmgrace using command :\n> %s %s &\n' %(ctrl.launchcommand,ctrl.agrfilename)
    os.system('%s %s &' %(ctrl.launchcommand,ctrl.agrfilename))
#=====================================================================================================================================================================


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
