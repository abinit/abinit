#=================================
# xmlTagger.py
version = '1.0'
#=================================
# last modified : january 17 2006
# written by Benjamin Tardif
# benjamin.tardif@umontreal.ca
#=================================
header = '\n#==============\n# xmlTagger.py\n# version %s\n#==============' %version

#====================================================================================================
#IMPORTS
import os
import sys
#====================================================================================================


#====================================================================================================
#METHODS
def detectfile(filename,path): # type(filename) = type(path) = string
    # method detectfile returns True if the specified file is found in the specified path
    return filename in os.listdir(path)

def clean(list): # type(list) = list of strings
    # method clean removes character strings '\n' and '\r' and empty lines from a string list
    # (the string list is usually obtained with the ".readlines()" method)
    L = len(list)
    for i in range(L):
        list[L-1-i] = list[L-1-i].replace('\n','')
        list[L-1-i] = list[L-1-i].replace('\r','')
        if list[L-1-i].split() == []:
            list.pop(L-1-i)
#====================================================================================================


#----------------------------------------------------------------------------------------------------
#MAIN
print header


#====================================================================================================
#COMMAND LINE

#get xmlfilename
if len(sys.argv) > 2:
    # user entered too many arguments in the command line
    print '\n- ERROR -\ntoo many arguments in the command line'
    sys.exit()
elif len(sys.argv) == 2:
    # user entered the xmlfilename in the command line
    xmlfilename = sys.argv[1]
else:
    # user entered no xmlfilename in the command line
    xmlfilename = raw_input('\nEnter the name of the xml file to tag :\n')

#abort if file not found
if detectfile(xmlfilename,'.') == False:
    print '\n- ERROR -\nfile not found\n'
    sys.exit()

#abort if the file is not a xml file
if xmlfilename[-4:] != '.xml':
    print '\n- ERROR -\nyou must enter a xml file (*.xml)\n'
    sys.exit()

#abort if the file is already a tagged xml file
if xmlfilename[-8:] == '_tag.xml':
    print '\n- ERROR -\nthis file is already tagged\n'
    sys.exit()
#====================================================================================================


#====================================================================================================
#READ AND TREAT THE FILE

#read the file
reader = open(xmlfilename,'r')
filedata = reader.readlines()
reader.close()
clean(filedata)

#for each line, remove all characters before '<' and after '>'
for i in range(len(filedata)):
    while filedata[i][0] != '<':
        filedata[i] = filedata[i][1:]
    while filedata[i][-1] != '>':
        filedata[i] = filedata[i][:-1]

#compute len_max (number of digits of the number of the last line of the xml file)
len_max = len(str(len(filedata)))

#compute tagxmlfilename (name of the tagged xml file)
tagxmlfilename = xmlfilename[:-4]+'_tag.xml'
#====================================================================================================


#====================================================================================================
#WRITE THE TAGGED XML FILE

writer = open(tagxmlfilename,'w')

tag=0
for line in filedata:
    if line.split()[0][1] == '/':
        # </Element>
        tag-=0
        len_tag = len(str(tag))
        writer.write((len_max+7)*' '+'%s\n' %line)
    elif line.split()[-1][-2] == '/':
        # <Element/>
        tag+=1
        len_tag = len(str(tag))
        writer.write((len_max-len_tag)*' '+'<!--%i-->'%tag+line[:-2]+" tag='%i'/>\n"%tag)
    else:
        # <Element>
        tag+=1
        len_tag = len(str(tag))
        writer.write((len_max-len_tag)*' '+'<!--%i-->'%tag+line[:-1]+" tag='%i'>\n"%tag)

writer.close()

print '\n"%s" file created successfully\n' %tagxmlfilename
#====================================================================================================


#----------------------------------------------------------------------------------------------------
