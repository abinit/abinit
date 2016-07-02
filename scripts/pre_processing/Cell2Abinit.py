#! /usr/bin/python
####################################################################
#   MaterialsStudio Cell File To Abinit Config file Converter
#                           zhoubo
#                       2006.9.26
####################################################################

import sys
import re

if __name__=="__main__":
    if  len(sys.argv) >2:
        useage()
        sys.exit()
    try:
        fp=open(sys.argv[1])
    except:
        print "Cannot open Cell file !"
        sys.exit()

    data=fp.readlines()

    #lattice paramenters
    a=b=c=0.0
    #primitive vectors
    a_=b_=c_=[]
    # atoms
    atoms=[]
    atompos=[]
    ntypat=""
    natom=0
    typat=""

    n=0
    while n<len(data):
        if data[n].find("%BLOCK LATTICE_CART")>-1:
            t1=re.split("\s+",data[n+1].strip())
            a_=map(float,t1)
            for x in a_:
                if x > a: a=x
            for x in xrange(len(a_)):
                a_[x]=a_[x]/a
            t1=re.split("\s+",data[n+2].strip())
            b_=map(float,t1)
            for x in b_:
                if x > b: b=x
            for x in xrange(len(b_)):
                b_[x]=b_[x]/b

            t1=re.split("\s+",data[n+3].strip())
            c_=map(float,t1)
            for x in c_:
                if x > c: c=x
            for x in xrange(len(c_)):
                c_[x]=c_[x]/c

            n+=4
            
        if data[n].find("%ENDBLOCK LATTICE_CART") > -1:
            n+=1
            continue

        if data[n].find("%BLOCK POSITIONS_FRAC") > -1:
            n+=1
            while data[n].find( "%ENDBLOCK POSITIONS_FRAC")==-1:
                t1=re.split("\s+",data[n].strip())
                if t1[0] not in atoms:
                    atoms.append(t1[0])
                    ntypat+=str(len(atoms))+" "
                    natom+=1
                    typat+=str(len(atoms))+" "
                    atompos.append(map(float,t1[1:]))
                else:
                    natom+=1
                    typat+=str(atoms.index(t1[0])+1)+" "
                    atompos.append(map(float,t1[1:]))
                n+=1
        if data[n].find("%ENDBLOCK POSITIONS_FRAC") > -1:
            break

        n+=1


    # Prepare the output
    print "####################################################"
    print "# Atoms Stucture created by CellToAbinit Converter !"
    print "####################################################"
    print ""
    print "# Definition of the unit cell"
    print "#The length of the primitive vectors "
    print "#                    1 Bohr=0.5291772108 Angstroms "
    print "# acell %f %f %f" % (round(a,6),round(b,6),round(c,6)),"Angstroms"
    print "acell %f %f %f" % (round(a/0.5291772108,6),round(b/0.5291772108,6),round(c/0.5291772108,6))
    print "rprim "
    print "    %8f  %8f  %8f"%(a_[0],a_[1],a_[2])
    print "    %8f  %8f  %8f"%(b_[0],b_[1],b_[2])
    print "    %8f  %8f  %8f"%(c_[0],c_[1],c_[2])
    print ""
    print "#Definition of the atom types"
    print "#      %d kind of atoms" %( natom)
    print "#      ",atoms    
    print "ntypat "+ntypat
    print "znucl  "+"Needed"
    print ""
    print "#Definition of tha atoms"
    print "natom ",natom
    print "typat ",typat
    print "xred"
    for x in atompos:
        print "   %8s %8s %8s"%(round(x[0],6),round(x[1],6),round(x[2],6))
    print ""
    print "##### End of CellToAbinit Conveter ! #####"
