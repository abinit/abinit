#!/usr/bin/python

import sys
import string

if len(sys.argv)<2 or len(sys.argv)>11:
    print ('Usage: python band_comp_abinit2abinit.py [abinit _EIG file1] ...')
    print ('          ... [abinit _EIG file2] [align #kpt #band] ...')
    print ('          ... [nbdbuf #] [Fermi #abinit_file_1 #abinit_file_2] [eV]')
    print ('          ... (align, eV, nbdbuf and Fermi sections optional)')
    sys.exit()

# Check third argument
eVconv = ''
str_nbdbuf = ''
nbdbuf = 0
align_values = 0
align_ikpt = 0
align_iband = 0
align_values_fermi = 0
align_abinit1_fermi = 0.0
align_abinit2_fermi = 0.0
iarg = 3
while iarg < len(sys.argv):
    if str(sys.argv[iarg])=='eV':
        eVconv = str(sys.argv[iarg]) # Assume Ha values and convert to eV
        print ('# Values assumed to be in Ha and converted to eV')
    if str(sys.argv[iarg])=='align':
        align_values = 1
        align_ikpt = int(sys.argv[iarg+1])
        align_iband = int(sys.argv[iarg+2])
        print ('# Values aligned at kpt:','%4i'%align_ikpt,\
	      ' and band:','%4i'%align_iband)
        iarg = iarg + 2
    if (str(sys.argv[iarg])=='Fermi' or \
            str(sys.argv[iarg])=='fermi'): # Align Fermi energies
        align_values_fermi = 1
        align_abinit1_fermi = float(sys.argv[iarg+1])
        align_abinit2_fermi = float(sys.argv[iarg+2])
        iarg = iarg + 2
    if str(sys.argv[iarg])=='nbdbuf':
        nbdbuf = int(sys.argv[iarg+1])
        print ('# nbdbuf set, last:','%4i'%nbdbuf,' bands will be ignored')
        iarg = iarg + 1
    iarg = iarg + 1

# Parse abinit file 1
input_file1_name = str(sys.argv[1]) # name of first input file (first argument)
input_file1_r = open(input_file1_name,'r')  # open it as read file       
abinit_file1_data = input_file1_r.readlines()
input_file1_r.close()
# Read in k-point data as a list of numbers, by k-point
abinit_band1_data = []
k_point_list = []
nspinpol = 1
# Check if we have a spin-polarised calc.
print (abinit_file1_data[0])
if abinit_file1_data[0].find('SPIN') > -1:
    print ('# System is spin-polarised')
    nspinpol = 2
print
for iline in range(2,len(abinit_file1_data)):
    # Skip line if it is a k-point spec.
    if abinit_file1_data[iline].find('kpt') > -1:
        continue
    # Accumulate values
#    new_values = map(float,string.split(abinit_file1_data[iline]))
    new_values = map(float,abinit_file1_data[iline].split())
    k_point_list.extend(new_values)
    # If we are on last line, finish appending
    if iline == len(abinit_file1_data)-1:
        abinit_band1_data.append(k_point_list)
        break
    # If the next line is a k-point spec., append.
    if abinit_file1_data[iline+1].find('kpt') > -1:
        abinit_band1_data.append(k_point_list)
        k_point_list = []
        continue

nkpt1 = len(abinit_band1_data)//nspinpol
nbands1 = len(abinit_band1_data[0])

# Parse abinit file 2
input_file2_name = str(sys.argv[2]) # name of second input file (second argument)
input_file2_r = open(input_file2_name,'r')  # open it as read file       
abinit_file2_data = input_file2_r.readlines()
input_file2_r.close()
# Read in k-point data as a list of numbers, by k-point
abinit_band2_data = []
k_point_list = []
for iline in range(2,len(abinit_file2_data)):
    # Skip line if it is a k-point spec.
    if abinit_file2_data[iline].find('kpt') > -1:
        continue
    # Accumulate values
    new_values = map(float,abinit_file2_data[iline].split())
    k_point_list.extend(new_values)
    # If we are on last line, finish appending
    if iline == len(abinit_file2_data)-1:
       abinit_band2_data.append(k_point_list)
       break
    # If the next line is a k-point spec., append.
    if abinit_file2_data[iline+1].find('kpt') > -1:
       abinit_band2_data.append(k_point_list)
       k_point_list = []
       continue

nkpt2 = len(abinit_band2_data)//nspinpol
nbands2 = len(abinit_band2_data[0])

# Check that the file contains the same number of
# nkpt and nbands
if nkpt2!=nkpt1 or nbands2!=nbands1:
   print (' ERROR: number of k-points or bands not the same!')
   sys.exit()


shift1 = 0.0
shift2 = 0.0
# If there is alignment at a certain k-pooint
if align_values:
    if align_ikpt>nkpt1:    
        print ('ERROR: index of kpt for alignment is larger than nkpt in data!')
    if align_iband>(nbands1-nbdbuf):    
        print ('ERROR: index of band for alignment is larger than nbands in data!')
    shift1 = abinit_band1_data[align_ikpt][align_iband]
    shift2 = abinit_band2_data[align_ikpt][align_iband]

# If Fermi alignment is done, print info
if align_values_fermi:
    print ('# Abinit file 1 Fermi energy (Ha): ',\
          '%18.9E'%align_abinit1_fermi)
    print ('# Abinit file 2 Fermi energy (Ha): ',\
          '%18.9E'%align_abinit2_fermi)
    shift1 = align_abinit1_fermi
    shift2 = align_abinit2_fermi


# Now we begin the reading and comparison of the data
Ha_to_eV = 27.21138386
avg_diff = []
input_file2_current_line='\n'
previous_kpt_val = -1.0
diff = []
max_diff = []
min_diff = []

for sppol in range(nspinpol):
    in_min_diff = 1000000000.0
    in_max_diff = -100000000.0
    in_avg_diff = 0.0
    nvals = 0
    for band in range((nbands1-nbdbuf)):
        bdiff = []
        for kpt in range(nkpt1):
#            print sppol,band,kpt
            # Calculate difference,average,max,min
            abinit_band1_data[kpt+nkpt1*sppol][band] = abinit_band1_data[kpt+nkpt1*sppol][band] - shift1
            abinit_band2_data[kpt+nkpt1*sppol][band] = abinit_band2_data[kpt+nkpt1*sppol][band] - shift2
            if eVconv=='eV':
                abinit_band1_data[kpt+nkpt1*sppol][band] = abinit_band1_data[kpt+nkpt1*sppol][band]*Ha_to_eV
                abinit_band2_data[kpt+nkpt1*sppol][band] = abinit_band2_data[kpt+nkpt1*sppol][band]*Ha_to_eV
                in_diff = abinit_band2_data[kpt+nkpt1*sppol][band] - abinit_band1_data[kpt+nkpt1*sppol][band]
            else:
                in_diff = abinit_band2_data[kpt+nkpt1*sppol][band] - abinit_band1_data[kpt+nkpt1*sppol][band]
            in_avg_diff = in_avg_diff + abs(in_diff)
            nvals = nvals + 1
            if in_diff<in_min_diff:
                in_min_diff = in_diff
            if in_diff>in_max_diff:
                in_max_diff = in_diff
            bdiff.append(in_diff)
#            print bdiff
        diff.append(bdiff)
    in_avg_diff = in_avg_diff/float(nvals)
    avg_diff.append(in_avg_diff)
    min_diff.append(in_min_diff)
    max_diff.append(in_max_diff)


#print nbands1,nkpt1,avg_diff,min_diff,max_diff
#print len(diff[0])
#print len(diff)
#print len(abinit_band2_data[0])
#print len(abinit_band2_data)

# Start the output
if nspinpol == 1:
    print ('#   k-point       abinit1 val          abinit2 val             diff        ')
    for band in range((nbands1-nbdbuf)):
        for kpt in range(nkpt1):
            if abs(abinit_band1_data[kpt][band])<0.1 or abs(abinit_band2_data[kpt][band])<0.1:
                print ('%7.2f'%float(kpt),'%18.9E'%abinit_band1_data[kpt][band],\
                      '%18.9E'%abinit_band2_data[kpt][band],'%12.3E'%diff[band][kpt])
            else:
                print ('%7.2f'%float(kpt),'%14.9f'%abinit_band1_data[kpt][band],\
	              '%18.9f'%abinit_band2_data[kpt][band],'%16.3E'%diff[band][kpt])
        print ('     ')
    print ('')
    print ('#')
    print ('#        nvals:','%5i'%nvals)
    if eVconv=='eV':
        print ('# average diff:','%12.6F'%avg_diff[0],' eV')
        print ('# minimum diff:','%12.6F'%min_diff[0],' eV')
        print ('# maximum diff:','%12.6F'%max_diff[0],' eV')
    else:
        print ('# average diff:','%12.6F'%avg_diff[0],' Ha')
        print ('# minimum diff:','%12.6F'%min_diff[0],' Ha')
        print ('# maximum diff:','%12.6F'%max_diff[0],' Ha')

else:
    print ('#   k-point       abinit1 val (spin up)  abinit2 val (spin up)    diff (spin up)   abinit1 val (spin down)  abinit2 val (spin down)    diff (spin down)    ')
    for band in range((nbands1-nbdbuf)):
        for kpt in range(nkpt1):
            if abs(abinit_band1_data[kpt][band])<0.1 or abs(abinit_band2_data[kpt][band])<0.1:
                print ('%7.2f'%float(kpt),'%18.9E'%abinit_band1_data[kpt][band],\
                      '%18.9E'%abinit_band2_data[kpt][band],'%12.3E'%diff[band][kpt],'   ',\
                      '%18.9E'%abinit_band1_data[kpt+nkpt1][band],\
                      '%18.9E'%abinit_band2_data[kpt+nkpt1][band],'%12.3E'%diff[band+nbands1][kpt])
            else:
                print ('%7.2f'%float(kpt),'%14.9f'%abinit_band1_data[kpt][band],\
	              '%18.9f'%abinit_band2_data[kpt][band],'%16.3E'%diff[band][kpt],'   ',\
                      '%14.9f'%abinit_band1_data[kpt+nkpt1][band],\
                      '%18.9f'%abinit_band2_data[kpt+nkpt1][band],'%16.3E'%diff[band+nbands1][kpt])
        print ('     ')
    print ('')
    print ('#')
    print ('#        nvals:','%5i'%nvals)
    if eVconv=='eV':
        print ('#   average diff (spin up):','%12.6F'%avg_diff[0],' eV')
        print ('#   minimum diff (spin up):','%12.6F'%min_diff[0],' eV')
        print ('#   maximum diff (spin up):','%12.6F'%max_diff[0],' eV')
        print ('#')
        print ('# average diff (spin down):','%12.6F'%avg_diff[1],' eV')
        print ('# minimum diff (spin down):','%12.6F'%min_diff[1],' eV')
        print ('# maximum diff (spin down):','%12.6F'%max_diff[1],' eV')
    else:
        print ('#   average diff (spin up):','%12.6F'%avg_diff[0],' Ha')
        print ('#   minimum diff (spin up):','%12.6F'%min_diff[0],' Ha')
        print ('#   maximum diff (spin up):','%12.6F'%max_diff[0],' Ha')
        print ('#')
        print ('# average diff (spin down):','%12.6F'%avg_diff[1],' Ha')
        print ('# minimum diff (spin down):','%12.6F'%min_diff[1],' Ha')
        print ('# maximum diff (spin down):','%12.6F'%max_diff[1],' Ha')

print ('#')
print ('# NOTE: Abinit values are read in fixed format with five decimal')
print ('#       places. For low values, four or three decimal figures')
print ('#       may be the highest precision you can get.')
