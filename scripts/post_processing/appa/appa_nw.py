#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Abinit Post Process Application
author: Martin Alexandre
last edited: April 2013
"""
import sys,readline,os,commands,subprocess

try:
    import fortran.math
except:
    print 'To use APPA, you first need to compile a FORTRAN file:'
    print '   cd ~appa/fortran'
    print '   f2py -c math.f90 -m math'
    sys.exit()

#Reading
import reading.read as Read
import reading.input as Input
#Utility
import utility.write as Write
import utility.analysis as Analysis
import utility.complet as Completer

try:
    import utility.canvas as Canvas
    canv = True
except:
    canv = False

#Fortran code
import fortran.math as Math

from numpy import linspace,zeros,where

#---------------------------------------------#
#  __  __      _   _               _          #
# |  \/  | ___| |_| |__   ___   __| |___   _  #
# | |\/| |/ _ \ __| '_ \ / _ \ / _` / __| (_) #
# | |  | |  __/ |_| | | | (_) | (_| \__ \  _  #
# |_|  |_|\___|\__|_| |_|\___/ \__,_|___/ (_) #
#					                                    #
#---------------------------------------------#     
                                      			      

readline.parse_and_bind("tab: complete")

#Units array (in progress):
units  = {'Pressure':['GPa',1],'Energy':["Ha",1],'Temperature':['K',0],'Volume':["Bohr^3",1],'Lengh':["Bohr",1],'Distance':["Angstrom",0.5291772085936]}


def displayData(data,deviation):
    return str("%.10e" %data) + ' +/- ' + str("%.1e" %deviation)

def openFile(choose = 'please choose file (help):'):
    while True :
        
        words = commands.getoutput('ls').split()
        completer = Completer.Completer(words)
        readline.set_completer(completer.complete)
        
        sub = raw_input('Command : ')
        
        if sub == 'ls':
            subprocess.call(sub, shell=True)
            
        elif sub == 'pwd':
            print commands.getoutput(sub)
            
        elif sub == 'cd ..':
            os.chdir('..')
            
        elif sub == 'cd':
            os.chdir(commands.getoutput('echo ~/'))
            
        elif sub.find('cd') == 0:
            dir_files = sub.replace('cd','').split()
            dir_file = str(dir_files[0])

            if len(dir_files) > 1 :
                for i in range(1,len(dir_files)):
                    dir_file += ' '+str(dir_files[i])		    	        
                
            try:
                os.chdir(dir_file)
            except:
                print 'No such file or directory : '+str(dir_file)
                
        elif (sub.find('.nc')!=-1 or sub.find('HIST')!=-1 or sub.find('out') != -1 ) :
            MD_file = Read.MolecularDynamicFile(str(os.getcwd())+'/'+str(sub))
            if MD_file.getNameFile() != '':
                print '\nloading successful!'
                print 'File :' + str(MD_file.getNameFile())+'\n\n'
                return MD_file   
                
                
        elif sub == 'clear':
            clear()
        
        
        elif sub == 'quit':
            break;
        elif (sub == 'help') or (sub == 'help'):
            print 'available commands:'
            print '    ls '
            print '    cd '
            print '    pwd '
            print '    clear '
            print '    quit'
        
        else:
            print 'Unknown command'

def saveFile(type_data,name = ""):
    fname = name
    if fname=="":
        fname = raw_input("Choose name file ("+str(type_data)+") : ")
        
    if str(type_data) in fname.split('.'):
        pass
    else:
        fname += '.'+str(type_data)
        fname = str(os.getcwd()) +'/' + fname 
    return fname



def showData(file,units):

    print "------------Simulation datas:------------"
                
    print 'Number of atoms : ' +  str(file.getNatom())        
    
    print 'Step: ' + str(file.getNbTime())
    print "initial step :  " + str(file.getNi())
    print "final   step :  " + str(file.getNf())                
                               
    E_tot = units['Energy'][1] *  Math.average(file.getE_Tot())
    deviation = file.getStandardDeviation(units['Energy'][1] * file.getE_Tot(), E_tot)        
    print 'Total Energy ('+ str(units['Energy'][0]) +"): " + displayData(E_tot,deviation)
        
        
    Vol = units['Volume'][1] * Math.average(file.getVol())
    deviation = file.getStandardDeviation( file.getVol() * units['Volume'][1] , Vol)
    print 'Volume ('+str(units['Volume'][0])+")  : " + displayData(Vol,deviation)
    
    Temp = file.getTemp()        
    ATemp =  Math.average(Temp) - units['Temperature'][1]        
    deviation = file.getStandardDeviation( Temp - units['Temperature'][1] , ATemp)
    print 'Temperature ('+str(units['Temperature'][0])+")  : " + displayData(ATemp,deviation)
        
                
                       
    Press =  Math.average(file.getPress() ) * units['Pressure'][1]
    deviation = file.getStandardDeviation(units['Pressure'][1] * file.getPress(), Press)
    print 'Pressure ('+str(units['Pressure'][0])+")   : "  + displayData(Press,deviation)    

    Acell =  Math.average(file.getAcell()[:,0]) * units['Distance'][1]
    deviation = file.getStandardDeviation(units['Distance'][1] * file.getAcell()[:,0], Acell)
    print 'a ('+str(units['Distance'][0])+")     : "  + displayData(Acell,deviation)    

    Acell =  Math.average(file.getAcell()[:,1]) * units['Distance'][1]
    deviation = file.getStandardDeviation(units['Distance'][1] * file.getAcell()[:,1], Acell)
    print 'b ('+str(units['Distance'][0])+")     : "  + displayData(Acell,deviation)    

    Acell =  Math.average(file.getAcell()[:,2]) * units['Distance'][1]
    deviation = file.getStandardDeviation(units['Distance'][1] * file.getAcell()[:,2], Acell)
    print 'c ('+str(units['Distance'][0])+")     : "  + displayData(Acell,deviation)    
    
    print "-----------------------------------------"
        
def logo():
    print"      ___                                     ___                              "
    print"     /  /\         _____        ___          /__/\        ___           ___    "
    print"    /  /::\       /  /::\      /  /\         \  \:\      /  /\         /  /\   "
    print"   /  /:/\:\     /  /:/\:\    /  /:/          \  \:\    /  /:/        /  /:/   "
    print"  /  /:/~/::\   /  /:/~/::\  /__/::\      _____\__\:\  /__/::\       /  /:/    "
    print" /__/:/ /:/\:\ /__/:/ /:/\:| \__\/\:\__  /__/::::::::\ \__\/\:\__   /  /::\    "
    print" \  \:\/:/__\/ \  \:\/:/~/:/    \  \:\/\ \  \:\~~\~~\/    \  \:\/\ /__/:/\:\   "
    print"  \  \::/       \  \::/ /:/      \__\::/  \  \:\  ~~~      \__\::/ \__\/  \:\  "
    print"   \  \:\        \  \:\/:/       /__/:/    \  \:\          /__/:/       \  \:\ "
    print"    \  \:\        \  \::/        \__\/      \  \:\         \__\/         \__\/ "
    print"     \__\/         \__\/                     \__\/                             "
    print""

def clear():
    os.system('/usr/bin/clear')
    logo()    
    try:
        global MD_file
        global units
        print 'File :' + str(MD_file.getNameFile())+'\n'
        if MD_file.isGoodFile():
            showData(MD_file,units)
    except:
        print "Loading File please wait..."
        #Try to read the last output:
        MD_file = Read.MolecularDynamicFile("")
        if MD_file.isGoodFile():
            print 'loading successful!'
            print 'File :' + str(MD_file.getNameFile())+'\n\n'
        else:
         #the last output does'nt exist        
            print 'Unable to load file,',  
            MD_file = openFile()
            if MD_file == None:
                sys.exit(0)
              
    
#---------------------------------#       
#  __  __    _    ___ _   _       #
# |  \/  |  / \  |_ _| \ | |  _   #
# | |\/| | / _ \  | ||  \| | (_)  #
# | |  | |/ ___ \ | || |\  |  _   #
# |_|  |_/_/   \_\___|_| \_| (_)  #
#                                 #
#-------------------------------- #                               
#Check the argument file (for input .appa)
if len(sys.argv)>1:
    logo()
    A = Input.Input(sys.argv[1])
    input_data =  A.read()
    if len(input_data)!=0:
        
        input_file = True
        size_input = len(input_data)
        i=0
    else:
        input_file = False
        size_input = 0 
        size_calculation=0
        i=-1
        j=-1
else:
    input_file = False
    clear()  
    ni = 0
    nf = 0
    showData(MD_file,units)
    size_input = 0
    size_calculation=0
    i=-1
    j=-1
    
#Main loop:
while i<size_input :

    choice = 0
    
    if input_file == True: 
        input_name=input_data.keys()[i]
        if os.path.isfile(input_name):
            pass
        else:
            print "\n!!!WARNING "+input_name+" doesnt exist"
            i+=1
            continue
       
        MD_file = Read.MolecularDynamicFile(str(input_name))
        if (MD_file.getNameFile() != '' and MD_file.isGoodFile()):
            print '\nFile :' + str(input_name)
            print 'loading successful!\n'
            #Set the inial step
            if("stepmin" in input_data[input_name]):
                try:
                    ni = int(input_data[input_name]["stepmin"])
                except:
                    ni = 1
                if ni >=1 and ni < MD_file.getNbTime():
                    MD_file.setNi(ni)
                else:
                    ni = 1
                    MD_file.setNi(ni)
                print "Set initial step to "+str(ni)
            #set the final step
            if("stepmax" in input_data[input_name]):
                try:
                    nf = int(input_data[input_name]["stepmax"])
                except:
                    nf = MD_file.getNbTime()
               
                if nf <= MD_file.getNbTime() and ni < nf :
                    MD_file.setNf(nf)
                else:
                    nf = MD_file.getNbTime()
                    MD_file.setNf(nf)
                print "Set final step to "+str(nf)     

            j=0
            size_calculation=len(input_data[input_name])
            choice = 4
            showData(MD_file,units)
    else:
        print " "
        print "your options are:"
        print " "
        print "1  => Open new File"
        print "2  => Change initial/final step"
        print "3  => Save simulation (XYZ format)"
        print "4  => Calculation of quantities"
        print "5  => Elastic constant calculation"
        print "6  => Clear"
        print "7  => Help"
        print "8  => Quit"
        print " "
        
        try:
            choice = input("Choose your option: ")
        except:
            pass
    
    if choice == 1:
        temp = openFile()
        if temp != None:
            MD_file = temp
            clear()
            
            
    if choice == 2:
        while True:
            if input_file == True:
                try:
                    MD_file.setNi(int(ni))
                    break;
                except:
                    pass
            else:
                try:
                    ni = input("Choose initial step : ")
                except :
                    pass; 
                if ni >=1 and ni < MD_file.getNbTime() : 
                    MD_file.setNi(int(ni))
                break;
        while 2:
            try :
                nf = input("Choose final step (ni<nf or 0=max) : ")
            except :
                pass;
            if nf==0:
                nf = MD_file.getNbTime()
                MD_file.setNf(int(MD_file.getNbTime()))
                break;
            else:
                if nf <= MD_file.getNbTime() and ni < nf : 
                    MD_file.setNf(int(nf))
                    break;            
        MD_file.setNf(int(nf))
        clear()
        
    elif choice == 3:
        try:    
            fname = saveFile('xyz')
            pos   = MD_file.getPos()
            acell = MD_file.getAcell()
            typat = MD_file.getTypat()
            znucl = MD_file.getZnucl()
            print "Saving File, please wait..."
            Write.SaveFile(fname).xyzFormat(pos,acell,typat,znucl)
            print  fname 
            print 'File save sucessful!'
        except:
            print 'Unable to save file'
            
            
    elif choice == 4:
        while j < size_calculation:
            choiceQuantitie = 0
            name_quantitie=""
            if input_file == False:
                print " "
                print "your options are:"
                print " "
                print "1  => Total Energy" 
                print "2  => Potential Energy"
                print "3  => Kinetic Energy" 
                print "4  => Temperature" 
                print "5  => Pressure" 
                print "6  => Stress" 
                print "7  => Acell" 
                print "8  => VACF" 
                print "9  => VDOS"  
                print "10 => RDF"
                print "11 => Help"
                print "12 => Return"
                print " "
        
                try:
                    choiceQuantitie = input("Choose your option: ")
                except :
                    pass;    
                
                if choiceQuantitie in [1,2,3,4,5,6,7,8,9]:
                    print " "
                    print "Choose format:"
                    print "1  => ASCII" 
                    if canv==True:
                        print "2  => PDF"
                        print "3  => PDF+ASCII"
                    else:
                        print "No canvas library to plot PDF"
                    print "4  => Return"

                    while True:        
                        try:
                            choiceformat = input("Choose your format: ")
                            if choiceformat == 4 :
                                choiceQuantitie = 0
                            break
                            if (canv==False and (choiceformat==3 or choiceformat==2)):
                                choiceQuantitie = 0
                                break
                        except :
                            pass; 
            else:
                calculation = input_data[input_name].keys()[j]
                if   calculation == "potential_energy":
                    choiceQuantitie=1
                elif calculation == "total_energy":
                    choiceQuantitie=2
                elif calculation == "kinetic_energy":
                    choiceQuantitie=3
                elif calculation == "temperatures":
                    choiceQuantitie=4
                elif calculation == "pressure":
                    choiceQuantitie=5
                elif calculation == "stress":
                    choiceQuantitie=6
                elif calculation == "acell":
                    choiceQuantitie=7
                elif calculation == "vacf":
                    choiceQuantitie=8
                elif calculation == "vdos":
                    choiceQuantitie=9
                    res = 0
                    if "vdosres" in input_data[input_name].keys():
                        res = input_data[input_name]["vdosres"]
                        
                        try:
                            res = int(res)
                            if (res <=0 or res > 100):
                                res = 8
                        except:
                            res = 8
                        
                    else:
                        res = 8
                elif calculation == "rdf":
                    pass

                if choiceQuantitie in [1,2,3,4,5,6,7,8]:
                    print calculation + " calculation"
                    name_quantitie=input_name+"_"+calculation

                choiceformat = input_data[input_name][calculation]
                try:
                    choiceformat = int(choiceformat)
                    if choiceformat in [1,2,3]:
                        pass
                    else:
                        choiceformat = 3
                except:
                    choiceformat = 3


            if choiceQuantitie == 1:
                
                x = linspace(MD_file.getNi(),MD_file.getNf()-1,\
                             MD_file.getNf()-MD_file.getNi())# Temporarily !!!
                y = MD_file.getE_Tot() * units['Energy'][1]                
                
                if choiceformat == 1 or choiceformat == 3:
                    fname = saveFile('dat',name_quantitie)
                    Write.SaveFile(fname).saveGraph(x,y,'Step',"Total Energy("+str(units['Energy'][0])+")")           
                    print  fname 
                    print 'File save sucessful!'
                if choiceformat == 2 or choiceformat == 3:
                    fname = saveFile('pdf',name_quantitie)
                    pdf = Canvas.Canvas(width=10, height=8, dpi=100,x=x,y=y,pxlbl='Step',pylbl="Total Energy("+str(units['Energy'][0])+")")
                    pdf.print_figure(fname)           
                    print  fname 
                    print 'File save sucessful!'
                
                
                
            elif choiceQuantitie == 2:
                x = linspace(MD_file.getNi(),MD_file.getNf()-1,MD_file.getNf()-MD_file.getNi())# Temporarily !!!
                y = MD_file.getE_pot() * units['Energy'][1]
                
                if choiceformat == 1 or choiceformat == 3:
                    fname = saveFile('dat',name_quantitie)        
                    Write.SaveFile(fname).saveGraph(x,y,'Step',"Potential Energy("+str(units['Energy'][0])+")")           
                    print  fname 
                    print 'File save sucessful!'
                if choiceformat == 2 or choiceformat == 3:
                    fname = saveFile('pdf',name_quantitie)
                    pdf = Canvas.Canvas(width=10, height=8, dpi=100,x=x,y=y,pxlbl='Step',pylbl="Potential Energy("+str(units['Energy'][0])+")")
                    pdf.print_figure(fname)
                    print  fname 
                    print 'File save sucessful!'
                
                
                
            elif choiceQuantitie == 3:
                    
                x = linspace(MD_file.getNi(),MD_file.getNf()-1,MD_file.getNf()-MD_file.getNi())# Temporarily !!!
                y = MD_file.getE_kin() * units['Energy'][1]                
                
                if choiceformat == 1 or choiceformat == 3:        
                    fname = saveFile('dat',name_quantitie)        
                    Write.SaveFile(fname).saveGraph(x,y,'Step',"Kinetic Energy("+str(units['Energy'][0])+")")
                    print  fname 
                    print 'File save sucessful!'
                if choiceformat == 2 or choiceformat == 3:
                    fname = saveFile('pdf',name_quantitie)
                    pdf = Canvas.Canvas(width=10, height=8, dpi=100,x=x,y=y,pxlbl='Step',pylbl="Kinetic Energy("+str(units['Energy'][0])+")")
                    pdf.print_figure(fname)
                    print  fname 
                    print 'File save sucessful!'
                
                
                
            elif choiceQuantitie == 4:
                        
                x = linspace(MD_file.getNi(),MD_file.getNf()-1,MD_file.getNf()-MD_file.getNi())# Temporarily !!!
                y = MD_file.getTemp() - units['Temperature'][1]
                
                if choiceformat == 1 or choiceformat == 3:        
                    fname = saveFile('dat',name_quantitie)        
                    Write.SaveFile(fname).saveGraph(x,y,'Step',"Temperature ("+str(units['Temperature'][0])+")")
                    print  fname 
                    print 'File save sucessful!'
                if choiceformat == 2 or choiceformat == 3:
                    fname = saveFile('pdf',name_quantitie)
                    pdf = Canvas.Canvas(width=10, height=8, dpi=100,x=x,y=y,pxlbl='Step',pylbl="Temperature ("+str(units['Temperature'][0])+")")
                    pdf.print_figure(fname)
                    print  fname 
                    print 'File save sucessful!'
                
            elif choiceQuantitie == 5:
                x = linspace(MD_file.getNi(),MD_file.getNf()-1,MD_file.getNf()-MD_file.getNi())# Temporarily !!!
                y = MD_file.getPress() * units['Pressure'][1]
                    
                if choiceformat == 1 or choiceformat == 3:        
                    fname = saveFile('dat',name_quantitie)        
                    Write.SaveFile(fname).saveGraph(x,y,'Step',"Pressure ("+str(units['Pressure'][0])+")")
                    print  fname 
                    print 'File save sucessful!'
                if choiceformat == 2 or choiceformat == 3:
                    fname = saveFile('pdf',name_quantitie)
                    pdf = Canvas.Canvas(width=10, height=8, dpi=100,x=x,y=y,pxlbl='Step',pylbl="Pressure ("+str(units['Pressure'][0])+")")
                    pdf.print_figure(fname)
                    print  fname 
                    print 'File save sucessful!'        
                        
                        
            elif choiceQuantitie == 6:
                x = linspace(MD_file.getNi(),MD_file.getNf()-1,MD_file.getNf()-MD_file.getNi())# Temporarily !!!
                y = MD_file.getStress()  * units['Pressure'][1]
                        
                if choiceformat == 1 or choiceformat == 3:        
                    fname = saveFile('dat',name_quantitie)        
                    Write.SaveFile(fname).saveGraph(x,y,'Step',"("+str(units['Pressure'][0])+") s" )
                    print  fname 
                    print 'File save sucessful!'
                if choiceformat == 2 or choiceformat == 3:
                    fname = saveFile('pdf',name_quantitie)
                    pdf = Canvas.Canvas(width=10, height=8, dpi=100,x=x,y=y,pxlbl='Step',pylbl="Stress ("+str(units['Pressure'][0])+")")
                    pdf.addLegend([r'$\sigma_1$',r'$\sigma_2$',r'$\sigma_3$',r'$\sigma_4$',r'$\sigma_5$',r'$\sigma_6$'])
                    pdf.print_figure(fname)
                    print  fname 
                    print 'File save sucessful!'

            elif choiceQuantitie == 7:
                x = linspace(MD_file.getNi(),MD_file.getNf()-1,MD_file.getNf()-MD_file.getNi())
                y = MD_file.getAcell() * units['Distance'][1]
                if choiceformat == 1 or choiceformat == 3:
                    fname = saveFile('dat',name_quantitie)        
                    Write.SaveFile(fname).saveGraph(x,y,'Step',"acell")
                    print  fname 
                    print 'File save sucessful!'
                if choiceformat == 2 or choiceformat == 3:
                    fname = saveFile('pdf',name_quantitie)
                    pdf = Canvas.Canvas(width=10, height=8, dpi=100,x=x,y=y,pxlbl='Step',pylbl="Acell"+" ("+str(units['Distance'][0])+")")
                    pdf.addLegend([r'$a$',r'$b$',r'$c$'],markerscale=1)
                    pdf.print_figure(fname)        
                    print  fname 
                    print 'File save sucessful!'
                    
                        
            elif choiceQuantitie == 8:
                y = Analysis.Correlation(MD_file.getVel()).getCorrelationFunction(normalize = True)
                x = linspace(MD_file.getNi(),MD_file.getNf()-1,MD_file.getNf()-MD_file.getNi())# Temporarily !!!
                    
                if choiceformat == 1 or choiceformat == 3:        
                    fname = saveFile('dat',name_quantitie)        
                    Write.SaveFile(fname).saveGraph(x,y,'Step',"VAF")
                    print  fname 
                    print 'File save sucessful!'
                if choiceformat == 2 or choiceformat == 3:
                    fname = saveFile('pdf',name_quantitie)
                    pdf = Canvas.Canvas(width=10, height=8, dpi=100,x=x,y=y,pxlbl='Step',pylbl="VAF")
                    pdf.print_figure(fname)        
                    print  fname 
                    print 'File save sucessful!'
                
                
            elif choiceQuantitie == 9:
                if input_file == False:
                    res = 0
                    while res<1 or res>100:
                        try:
                            res = input("Choose gaussian width (1-100) : ")
                        except:
                            pass

                vacf = Analysis.Correlation(MD_file.getVel()).getCorrelationFunction(normalize = True)
                pdos = Analysis.DOS(vacf,res,MD_file.getDtion())

                y = pdos.getDOS()
                x = pdos.getFrequencies()
                
                if choiceformat == 1 or choiceformat == 3:        
                    fname = saveFile('dat',name_quantitie)        
                    Write.SaveFile(fname).saveGraph(x,y,'E (meV)',"Phonons DOS (nm^2/ps)")
                    print  fname 
                    print 'File save sucessful!'
                if choiceformat == 2 or choiceformat == 3:
                    fname = saveFile('pdf',name_quantitie)
                    pdf = Canvas.Canvas(width=10, height=8, dpi=100,x=x,y=y,pxlbl='E (meV)',pylbl= "Phonons DOS (nm^2/ps)",adjust=True)
                    pdf.print_figure(str(fname))        
                    print  fname 
                    print 'File save sucessful!'            
            

            elif choiceQuantitie == 10:
                    
                PTOE = Analysis.PeriodicTableElement()
                nameAtom = MD_file.getAtomName()
                
                if len(nameAtom) == 1:
                    atom0 = 1
                    atom1 = 1
                else:
                    for i in range(2):
                        stringChoice = 'Choose the atom' + str(i+1)+' ( '
                        for k in range(len(nameAtom)):
                            stringChoice += str(nameAtom[k])
                            if k != len(nameAtom)-1:
                                stringChoice += ' or '
                            else:
                                stringChoice+= ' ):'
                        while True:
                            globals()['atom%s' % i] = raw_input(stringChoice)
                            if globals()['atom%s' % i] in nameAtom:
                                globals()['atom%s' % i] = nameAtom.index(globals()['atom%s' % i])+1
                                
                        
                while True:
                    try:
                        prmax = input("Choose maximum radius ("+str(MD_file.getRPrim()[0][0]*0.5)+") : ")
                        if prmax <= MD_file.getRPrim()[0][0]:
                            break
                    except:
                        pass;                        

                while True:
                    try:
                        pdr = input("Choose radius incrementation (0.1) : ")
                        if pdr > 0:
                            break
                    except:
                            pass;        

                while True:
                    try:
                        pstep = input("Choose time step incrementation (15) : ")
                        if pstep > 0:
                            break
                    except:
                        pass;        

                print 'please wait....'  

                try:
                    rdf = Analysis.RDF(MD_file,atom0,atom1, prmax, pdr, pstep,No_Windows=True)
                except:
                    print 'Unable to calculate le radial density function'
                    break

                x = rdf.getR()
                y = rdf.getRDF()

                if choiceformat == 1 or choiceformat == 3:        
                    fname = saveFile('dat',name_quantitie)        
                    Write.SaveFile(fname).saveGraph(x,y,'R (bohr)',"RDF")
                    print  fname 
                    print 'File save sucessful!'
                if choiceformat == 2 or choiceformat == 3:
                    fname = saveFile('pdf',name_quantitie)
                    pdf = Canvas.Canvas(width=10, height=8, dpi=100,x=x,y=y,pxlbl='R (bohr)',pylbl= "RDF")
                    pdf.print_figure(str(fname))        
                    print  fname 
                    print 'File save sucessful!'                            
                break;


            elif choiceQuantitie == 12:
              break;

            if (input_file==True):
                j+=1
          
          
    elif choice == 5:
        files = {} 
        while True:
            try:
                nbfile = input("Choose the number of file : ")
                if nbfile > 0:
                    break
            except:
                pass;
        while True:
            try:
                MD_file = openFile(choose = 'Choose the file without deformation (help):')
                if MD_file.isGoodFile():
                    files['WD']=MD_file
                    break
            except:
                pass;
        for i in range(nbfile-1):
            while True:
                try:
                    MD_file = openFile(choose = 'Choose the file with deformation (help):')
                    if MD_file.isGoodFile():
                        files[i]=MD_file
                        break
                except:
                    pass;
        while True:
            try:
                step = input("Choose initial step : ")
                if step >=1 and step < files['WD'].getNbTime() : 
                    break;
            except:
                pass;
        elasticC = Analysis.elasticCalculation(files,step,1).getElasticsConstants()
        print '\nThe elastic constants are:'        
        print elasticC
         
    elif choice == 6:
        clear()        
        
    elif choice == 7:
        clear();        
        
    elif choice == 8:
        break;

    if(input_file == True):
        i+=1

print "Thank you for using Abinit Post Process, see you soon!"
