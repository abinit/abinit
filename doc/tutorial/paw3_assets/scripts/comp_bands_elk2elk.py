#!/usr/bin/python

import sys, string 

if len(sys.argv)<3 or len(sys.argv)>4:
    print 'Usage: python band_comp_elk_elk.py [elk BAND.OUT file 1] ...'
    print '          ... [elk BAND.OUT file 2] [eV](optional)'
    sys.exit()

# Initialise variables
input_file1_name = str(sys.argv[1]) # name of first input file (first argument)
input_file1_r = open(input_file1_name,'r')  # open it as read file       
input_file2_name = str(sys.argv[2])  # name of second input file(second argument)
input_file2_r = open(input_file2_name, "r")  # open it as read file
eVconv = ''
if len(sys.argv)==4:
    eVconv = str(sys.argv[3]) # Assume Ha values and convert to eV
if eVconv!='':
    if eVconv!='eV':
        sys.exit('If third argument is provided, it must be: eV ')
    else:
        print '# Values assumed to be in Ha and converted to eV'

print '#   k-point       band val 1     band val 2       diff        %diff'

# Now we begin the reading and comparison of the data
Ha_to_eV = 27.21138386
avg_diff = 0.0
avg_percentage = 0.0
nvals = 0
min_diff = 1000000000.0
max_diff = -100000000.0
input_file1_current_line='\n'
while input_file1_current_line.endswith('\n'):
    input_file1_current_line = input_file1_r.readline()
    input_file2_current_line = input_file2_r.readline()
    if not input_file1_current_line.endswith('\n'): # stop at last line
        break
    if string.split(input_file1_current_line) == []: # skip blank lines
        print '     '
        continue
    # Parse the records of the lines
    numbers_list1 = []
    numbers1 = map(float, string.split(input_file1_current_line))
    numbers_list1.extend(numbers1)
    numbers_list2 = []
    numbers2 = map(float, string.split(input_file2_current_line))
    numbers_list2.extend(numbers2)
    if numbers_list1[0]!=numbers_list2[0]:
        sys.exit('ERROR: k-points differ!')
    # Calculate difference,average,max,min
    if eVconv=='eV':
        numbers_list1[1] = numbers_list1[1]*Ha_to_eV
	numbers_list2[1] = numbers_list2[1]*Ha_to_eV
        diff = numbers_list2[1] - numbers_list1[1]
    else:
        diff = numbers_list2[1] - numbers_list1[1]
    if numbers_list2[1]<=0.0:
        percentage = (numbers_list1[1]/numbers_list2[1] - 1.0)*100.0
    else:
        percentage = (numbers_list1[1]/numbers_list2[1] - 1.0)*100.0
    avg_diff = avg_diff + abs(diff)
    avg_percentage = avg_percentage + abs(percentage)
    nvals = nvals + 1
    if diff<min_diff:
        min_diff = diff
	min_percentage = percentage
    if diff>max_diff:
        max_diff = diff
	max_percentage = percentage
    if abs(numbers_list1[1])<0.1 or abs(numbers_list2[1])<0.1:
        print '%14.9f'%numbers_list1[0],'%18.9E'%numbers_list1[1],\
              '%18.9E'%numbers_list2[1],'%12.3E'%diff,' (','%7.3f'%percentage,'%)'
    else:
        print '%14.9f'%numbers_list1[0],'%14.9f'%numbers_list1[1],\
	      '%18.9f'%numbers_list2[1],'%16.3E'%diff,' (','%7.3f'%percentage,'%)'

avg_diff = avg_diff/float(nvals)
avg_percentage = avg_percentage/float(nvals)

print '#        nvals:','%5i'%nvals
print '# average diff:','%12.3E'%avg_diff,' (','%7.3f'%avg_percentage,'%)'
print '# minimum diff:','%12.3E'%min_diff,' (','%7.3f'%min_percentage,'%)'
print '# maximum diff:','%12.3E'%max_diff,' (','%7.3f'%max_percentage,'%)'
