#!/usr/bin/python

import sys, string 

if len(sys.argv)<2 or len(sys.argv)>11:
    print ('Usage: python comp_band_abinit2elk.py [abinit _EIG file] ...')
    print ('          ... [elk BAND.OUT file] [align #kpt #band] ...')
    print ('          ... [nbdbuf #] [Fermi #abinit #elk]')
    print ('              (align, eV, nbdbuf and Fermi sections optional)')
    sys.exit()


# Check arguments
eVconv = ''
str_nbdbuf = ''
nbdbuf = 0
align_values = 0
align_ikpt = 0
align_iband = 0
align_values_fermi = 0
align_abinit_fermi = 0.0
align_elk_fermi = 0.0
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
        align_abinit_fermi = float(sys.argv[iarg+1])
        align_elk_fermi = float(sys.argv[iarg+2])
        iarg = iarg + 2
    if str(sys.argv[iarg])=='nbdbuf':
        nbdbuf = int(sys.argv[iarg+1])
        print ('# nbdbuf set, last:','%4i'%nbdbuf,' bands will be ignored')
        iarg = iarg + 1
    iarg = iarg + 1

# Parse abinit file
input_file1_name = str(sys.argv[1]) # name of first input file (first argument)
input_file1_r = open(input_file1_name,'r')  # open it as read file       
abinit_file_data = input_file1_r.readlines()
input_file1_r.close()
# Read in k-point data as a list of numbers, by k-point
abinit_band_data = []
k_point_list = []
for iline in range(2,len(abinit_file_data)):
    # Skip line if it is a k-point spec.
    if abinit_file_data[iline].find('kpt#') > -1:
        continue
    # Accumulate values
    new_values = map(float,abinit_file_data[iline].split())
    k_point_list.extend(new_values)
    # If we are on last line, finish appending
    if iline == len(abinit_file_data)-1:
        abinit_band_data.append(k_point_list)
        break
    # If the next line is a k-point spec., append.
    if abinit_file_data[iline+1].find('kpt#') > -1:
        abinit_band_data.append(k_point_list)
        k_point_list = []
        continue

#print 'nkpt:',len(abinit_band_data)       
#print 'nbnds:',len(abinit_band_data[0])

nkpt = len(abinit_band_data)
nbands = len(abinit_band_data[0])

# Initialise variables for elk file
input_file2_name = str(sys.argv[2])  # name of second input file(second argument)
input_file2_r = open(input_file2_name, "r")  # open it as read file

# If alignment is done, parse elk file and find alignment value
if align_values:
    elk_file_data = input_file2_r.readlines()
    band = 0
    kpt = -1
    previous_kpt_val = -1.0
    for iline in range(0,len(elk_file_data)-1):
        kpt = kpt + 1
        if elk_file_data[iline].split() == []: # skip blank lines
           kpt = -1
           band = band + 1
           continue
        val_list = list(map(float,elk_file_data[iline].split()))
        elk_align_value = val_list[1]
        if previous_kpt_val==elk_align_value: # skip repeated kpts
            kpt = kpt - 1
            band = band - 1
            continue
        if kpt+1==align_ikpt and band+1==align_iband:
            break
    abinit_align_value = abinit_band_data[align_ikpt-1][align_iband-1]
    print ('# Abinit alignment value: ',\
          '%18.9E'%abinit_align_value)
    print ('#    Elk alignment value: ',\
          '%18.9E'%elk_align_value)
    input_file2_r.close()
    input_file2_r = open(input_file2_name, "r")

# If Fermi alignment is done, print info
if align_values_fermi:
    print ('# Abinit Fermi energy (Ha): ',\
          '%18.9E'%align_abinit_fermi)
    print ('#    Elk Fermi energy (Ha): ',\
          '%18.9E'%align_elk_fermi)

print ('#   k-point       abinit val          elk val             diff        %diff')

# Now we begin the reading and comparison of the data
Ha_to_eV = 27.21138386
avg_diff = 0.0
avg_percentage = 0.0
nvals = 0
min_diff = 1000000000.0
max_diff = -100000000.0
input_file2_current_line='\n'
band = 0
kpt = -1
previous_kpt_val = -1.0
reset_done = 1
while input_file2_current_line.endswith('\n'):
    kpt = kpt + 1
    #print kpt
    input_file2_current_line = input_file2_r.readline()
    if not input_file2_current_line.endswith('\n'): # stop at last line
        break
    if input_file2_current_line.split() == []: # skip blank lines
        kpt = -1
        band = band + 1
        if band == nbands-nbdbuf:
          break
        print ('     ')
        continue
    # Parse the records of the lines
    numbers_list2 = []
    numbers2 = map(float, input_file2_current_line.split())
    numbers_list2.extend(numbers2)
    # Cycle if k-point value is repeated in elk
    if numbers_list2[0] == previous_kpt_val:
        kpt = kpt - 1
        continue
    previous_kpt_val = numbers_list2[0]
    # Calculate difference,average,max,min
    if align_values:
        abinit_band_data[kpt][band] = abinit_band_data[kpt][band] \
	              - abinit_align_value 
        numbers_list2[1] = numbers_list2[1] - elk_align_value
    if align_values_fermi:
        abinit_band_data[kpt][band] = abinit_band_data[kpt][band] \
	              - align_abinit_fermi 
        numbers_list2[1] = numbers_list2[1] - align_elk_fermi
    if eVconv=='eV':
        abinit_band_data[kpt][band] = abinit_band_data[kpt][band]*Ha_to_eV
        numbers_list2[1] = numbers_list2[1]*Ha_to_eV
    diff = numbers_list2[1] - abinit_band_data[kpt][band]
    if numbers_list2[1]!=0.0:
        percentage = (abinit_band_data[kpt][band]/numbers_list2[1] - 1.0)*100.0
    else:
        percentage = 0.0
    avg_diff = avg_diff + abs(diff)
    avg_percentage = avg_percentage + abs(percentage)
    nvals = nvals + 1
    if diff<min_diff:
        min_diff = diff
        min_percentage = percentage
    if diff>max_diff:
        max_diff = diff
        max_percentage = percentage
    if align_values and band+1==align_iband and reset_done: # Save current averages for the
        occ_avg_diff = avg_diff             # occupied states and reset
        occ_avg_percentage = avg_percentage
        occ_min_diff = min_diff
        occ_min_percentage = min_percentage
        occ_max_diff = max_diff
        occ_max_percentage = max_percentage
        occ_nvals = nvals
        nvals = 0
        avg_diff = 0.0
        avg_percentage = 0.0
        min_diff = 1000000000.0
        max_diff = -100000000.0
        reset_done = 0

    # Print the output data
    if abs(abinit_band_data[kpt][band])<0.1 or abs(numbers_list2[1])<0.1:
        print ('%14.9f'%numbers_list2[0],'%18.9E'%abinit_band_data[kpt][band],\
              '%18.9E'%numbers_list2[1],'%12.3E'%diff,' (','%7.3f'%percentage,'%)')
    else:
        print ('%14.9f'%numbers_list2[0],'%14.9f'%abinit_band_data[kpt][band],\
	      '%18.9f'%numbers_list2[1],'%16.3E'%diff,' (','%7.3f'%percentage,'%)')

input_file2_r.close

avg_diff = avg_diff/float(nvals)
avg_percentage = avg_percentage/float(nvals)

if eVconv=='eV':
  eVconv=' eV'
else:
  eVconv=' Ha'

if align_values:
    occ_avg_diff = occ_avg_diff/float(occ_nvals)
    occ_avg_percentage = occ_avg_percentage/float(occ_nvals)
    print ('#')
    print ('# AVERAGES FOR OCCUPIED STATES:') 
    print ('#        nvals:','%5i'%occ_nvals) 
    print ('# average diff:','%12.6F'%occ_avg_diff,eVconv)
    print ('# minimum diff:','%12.6F'%occ_min_diff,eVconv)
    print ('# maximum diff:','%12.6F'%occ_max_diff,eVconv)
    print ('#')
    print ('# AVERAGES FOR UNOCCUPIED STATES:')
    print ('#        nvals:','%5i'%nvals)
    print ('# average diff:','%12.6F'%avg_diff,eVconv)
    print ('# minimum diff:','%12.6F'%min_diff,eVconv)
    print ('# maximum diff:','%12.6F'%max_diff,eVconv)
    print ('#')
    print ('# NOTE: Abinit values are read in fixed format with five decimal')
    print ('#       places. For low values, ~E-04 or ~E-03 may be the highest')
    print ('#       precision you can get.')
else:
    print ('#        nvals:','%5i'%nvals)
    print ('# average diff:','%12.6F'%avg_diff,eVconv)
    print ('# minimum diff:','%12.6F'%min_diff,eVconv)
    print ('# maximum diff:','%12.6F'%max_diff,eVconv)
    print ('#')
    print ('# NOTE: Abinit values are read in fixed format with five decimal')
    print ('#       places. For low values, ~E-04 or ~E-03 may be the highest')
    print ('#       precision you can get.')

