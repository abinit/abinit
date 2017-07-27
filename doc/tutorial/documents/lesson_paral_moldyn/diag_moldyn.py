#!/usr/bin/env  python
import os, sys, time
import string, math
from scipy.io  import netcdf

#-----------------------------------------------------------------|
def moy (file_in) :
#-----------------------------------------------------------------|
	FILE_IN = open(file_in, 'r')
	neff= len(FILE_IN.readlines())
	FILE_IN.close()
	# calcul des moyennes et des fluctuations
	FILE_IN = open(file_in, 'r')
	sum = 0.
	sum2 = 0.
	for line in FILE_IN.readlines() :
		y = string.atof(line)
		sum = sum + y
		sum2 = sum2 + y * y
	amoy = sum / float(neff)
	fluc = math.sqrt(abs(sum2 / float(neff) - amoy * amoy))
	prec = 1.
	if abs(amoy) > prec :
			print '%6s = %8.1f +/ %8.4f ' % (file_in, amoy, fluc)
	else :
			print '%6s = %14.5e +/ %14.5e ' % (file_in, amoy, fluc)
	FILE_IN.close()

#-----------------------------------------------------------------|
def lireNetcdf() :
#-----------------------------------------------------------------|
#Recherche les derniers fichiers *_HIST/*_OUT.nc cree
#Recherche le dernier fichier *_MOLDYN.nc cree
	lesfichiers = os.listdir(os.getcwd())
	f_fdr ='ok' ; fiche='' ; OUT_list = [] ; MOLDYN_list = []
	for fichier in lesfichiers :
		stats = os.stat(fichier)
		fic_tuple = time.localtime(stats[8]), fichier
		if fichier.find('_OUT.nc') == len(fichier)-7 and len(fichier)>7 :
			OUT_list.append(fic_tuple)
		if fichier.find('_MOLDYN.nc') == len(fichier)-10 and len(fichier)>10 :
			MOLDYN_list.append(fic_tuple)
	ok_OUT=-1;ok_MOLDYN=-1
	if len(OUT_list) > 0 :
		OUT_list.sort() ; OUT_list.reverse()
		fic_HIST=OUT_list[0][1].replace('_OUT.nc','_HIST.nc')
		if os.path.exists(fic_HIST) :
			fiche = fic_HIST
			ok_OUT=0
	elif ok_OUT==-1 and len(MOLDYN_list) > 0 :
		MOLDYN_list.sort() ; MOLDYN_list.reverse()
		fiche = MOLDYN_list[0][1]
		ok_MOLDYN=0

	if ok_OUT == -1 and ok_MOLDYN == -1 :
		print "Vous n'avez pas de fichier netcdf pour diagnostics !"
		sys.exit(1)

	ConvNetCdf_Ascii(fiche)
	return f_fdr


def ConvNetCdf_Ascii(file) :
	#Lit les fichiers mis au format netcdf et les transforme en ASCII
        #Deux cas: fichier _MOLDYN.nc (old format) ou fichier _HIST/_OUT.nc (new format)
	
        Ha_eV = 27.2113834e0   # Facteur de conversion 1 Hartree = 27,2113834 eV
	kb_eVK = 8.617342e-5   # Cste de Boltzman en eV/K  = (1.380658e-23)/(1.60217733e-19) #eV/K
	amu_emass = 1.660538782e-27/9.10938215e-31 # 1 atomic mass unit, in electronic mass
	fact = 29421.033e0     # Facteur de conversion  hartree/bohr^3 en GPa

        TypeFichier=0
	if file.find('_MOLDYN.nc') == len(file)-10 :
		TypeFichier=1
		fichier1=file
		fichier2=''
	elif file.find('_HIST.nc') == len(file)-8 :
		TypeFichier=2
		fichier1=file
		fichier2=file.replace('_HIST.nc','_OUT.nc')
	elif file.find('_OUT.nc') == len(file)-7 :
		TypeFichier=3
		fichier1=file.replace('_OUT.nc','_HIST.nc')
		fichier2=file
	if TypeFichier ==0 :
		print "Bug (formats des fichiers) !"
		sys.exit(1)

	#Cas 1 ====================================================================
	if TypeFichier == 1 :
		ncfile = netcdf.netcdf_file(fichier1, 'r')
		Nbatomes= ncfile.dimensions['NbAtoms']
		DimTensor= ncfile.dimensions['DimTensor']
		Vol = ncfile.variables['Cell_Volume'].data[0]
		#Potential, kinetic and total energies
		E_pot= Ha_eV*ncfile.variables['E_pot'].data
		Nbtimes= E_pot.shape[0]
		EcrireFichierTemporaire('EPOT',[Nbtimes], E_pot)
		E_kin= Ha_eV*ncfile.variables['E_kin'].data
		EcrireFichierTemporaire('ECIN',[Nbtimes], E_kin)
		ETOT= E_kin + E_pot
		EcrireFichierTemporaire('ETOT',[Nbtimes], ETOT)
		#Temperature
		TEMPER = 2.* E_kin/(3. * kb_eVK * Nbatomes)
		EcrireFichierTemporaire('TEMPER',[Nbtimes], TEMPER)
		#Stress stensor
		Stress=  ncfile.variables['Stress'].data
		EcrireFichierTemporaire('Stress',[Nbtimes,DimTensor], Stress)
		#Pressure in GPa : -1/3 trace(stress tensor)  + (Nbatomes/V)*kb_eVK*T
		Pression = zeros(Nbtimes, float64)
		for i in range ( 0,  Nbtimes) :
			Pionique= (Nbatomes/Vol)* kb_eVK * TEMPER[i] /Ha_eV #Hartree/bohr^3
			Pression[i]= -(Stress[i][0]+ Stress[i][1] + Stress[i][2])/3. + Pionique
			Pression[i]= fact * Pression[i] #GPa
		EcrireFichierTemporaire('PRESS',[Nbtimes], Pression)
		ncfile.close()
	#Cas 2 or 3 ===============================================================
	if TypeFichier == 2 or TypeFichier ==3 :
		ncfile1 = netcdf.netcdf_file(fichier1, 'r')
		ncfile2 = netcdf.netcdf_file(fichier2, 'r')
		#Simulation cell data
		natom=ncfile1.dimensions['natom']
		strten_size=ncfile1.dimensions['six']
		typat=ncfile2.variables['typat'].data
		amu=ncfile2.variables['amu'].data
		mass=[]
		for typ in typat:
			mass.extend([amu_emass*amu[typ-1]])
		rprimd=ncfile1.variables['rprimd'].data[0]
		Vol=rprimd[0][0]*(rprimd[1][1]*rprimd[2][2]-rprimd[1][2]*rprimd[2][1]) \
		 +  rprimd[0][1]*(rprimd[1][2]*rprimd[2][0]-rprimd[1][0]*rprimd[2][2]) \
		 +  rprimd[0][2]*(rprimd[1][0]*rprimd[2][1]-rprimd[1][1]*rprimd[2][0])
		#Potential energy
		E_pot = Ha_eV*ncfile1.variables['etotal'].data
		Nbtimes =E_pot.shape[0] -1  # Temporarily !!!
		EcrireFichierTemporaire('EPOT',[Nbtimes], E_pot)
		#Kinetic energy from velocities and masses
		vel = ncfile1.variables['vel'].data
		E_kin=[]
		for itim in range(Nbtimes):
			jtim=itim+1 # Temporarily
			Ek=0.
			for iat in range(natom):
				Ek=Ek+mass[iat]*(vel[jtim][iat][0]**2+vel[jtim][iat][1]**2+vel[jtim][iat][2]**2)
			E_kin.extend([0.5*Ha_eV*Ek])
		EcrireFichierTemporaire('ECIN',[Nbtimes], E_kin)
		#Total energy
		ETOT = E_kin[0:Nbtimes] + E_pot[0:Nbtimes]
		EcrireFichierTemporaire('ETOT',[Nbtimes], ETOT)
		#Temperature
		TEMPER=[]
		for itim in range(Nbtimes):
			TEMPER.extend([2.*E_kin[itim]/(3.*kb_eVK*natom)])
		EcrireFichierTemporaire('TEMPER',[Nbtimes], TEMPER)
		#Stress stensor
		Stress=ncfile1.variables['strten'].data
		EcrireFichierTemporaire('Stress',[Nbtimes, strten_size], Stress)
		#Pressure
		Pression=[]
		for itim in range(Nbtimes):
			Pionique = (natom/Vol)*kb_eVK*TEMPER[itim]/Ha_eV
			PP=-(Stress[itim][0]+Stress[itim][1]+Stress[itim][2])/3. + Pionique
 			Pression.extend([fact*PP])
		EcrireFichierTemporaire('PRESS',[Nbtimes], Pression)
		#Positions
		xred = ncfile1.variables['xred'].data
        file_xred = open('XRED',"w")
        for itim in range(Nbtimes):
            line="# TIME STEP %d \n" % (itim+1)
            file_xred.write(line)
            for ia in range(natom):
                line='%.10e  %.10e  %.10e \n' % (xred[itim][ia][0],xred[itim][ia][1],xred[itim][ia][2])
                file_xred.write(line)
        file_xred.close
     	xcart = ncfile1.variables['xcart'].data
        file_xcart = open('XCART',"w")
        for itim in range(Nbtimes):
            line="# TIME STEP %d \n" % (itim+1)
            file_xcart.write(line)
            for ia in range(natom):
                line='%.10e  %.10e  %.10e \n' % (xcart[itim][ia][0],xcart[itim][ia][1],xcart[itim][ia][2])
                file_xcart.write(line)
        file_xcart.close

        ncfile1.close()
        ncfile2.close() 

def EcrireFichierTemporaire (nom_fic , dimensions,  data) :
    file_out = open( nom_fic ,"w")
    if (len (dimensions) == 1 ) :
	for i in range (0, dimensions[0]):
	    line = "%.10e\n"   % data[i]
	    file_out.write(line)
    if (len (dimensions) == 2 ) :
	for i in range (0, dimensions[0]):
	    line = ' '
	    for j in range (0, dimensions[1]):
		lineC = "%.10e "   % data[i][j]
		line = line + lineC
	    line = line + "\n"
	    file_out.write(line)
    file_out.close()

#-----------------------------------------------------------------|
#-----------------------------------------------------------------|
#----------------   Program main     -----------------------------|
#-----------------------------------------------------------------|
#-----------------------------------------------------------------|
nbArg = len(sys.argv) - 1
localDir=os.getcwd()
FDR = lireNetcdf()
FILE_IN = open('ETOT', 'r')
neff= len(FILE_IN.readlines())
FILE_IN.close()
moy('ETOT')
moy('TEMPER')
moy('PRESS')
