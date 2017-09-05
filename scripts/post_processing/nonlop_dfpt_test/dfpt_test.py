#! /usr/bin/env python
# -*- coding: iso-8859-15 -*-
#set -x

# ===============================================================
# = Tool for ABINIT nonlop routine testing                      =
# =                                                             =
# = Uses an instrumented version of forstrnps routine;          =
# = To activate this instrumented version, uncomment all        =
# = 'TESTDFPT' sections in forstrnps.F90 file.                 =
# ===============================================================

import os,sys

# TO BE CUSTOMIZED BY USER
ABINIT_EXE='../../../build/atlas-fftw3/src/98_main/abinit'

#Read argument(s)
n_arg=len(sys.argv)-1
if n_arg < 1 or sys.argv[1] == "--help"  or sys.argv[1] == "-help":
  print >> sys.stderr, 'Syntax: '+ sys.argv[0]+' --choice choice [--test_case test_case]'
  print >> sys.stderr, '         [--signs signs] [--signsdfpt signsdfpt]'
  print >> sys.stderr, '         [--idir1 idir1] [--idir2 idir2] [--iatom iatom]'
  print >> sys.stderr, '         [--iband iband] --enl [enl]'
  print >> sys.stderr, '   choice= 2  : d/d_atm (force)'
  print >> sys.stderr, '           3  : d/d_+strain (stress)'
  print >> sys.stderr, '           5  : d/d_k (ddk)'
  print >> sys.stderr, '           51 : d(right)/d_k, partial dk derivative'
  print >> sys.stderr, '           54k: d2/d_atm.d_k, finite difference on k (effective charge)'
  print >> sys.stderr, '           54a: d2/d_atm.d_k, finite difference on atm (effective charge)'
  print >> sys.stderr, '           55k: d2/d_strain.d_k, finite difference on k (piezo)'
  print >> sys.stderr, '           55s: d2/d_strain.d_k, finite difference on strain (piezo)'
  print >> sys.stderr, '           8  : d2/d_k.d_k, full dk derivative (effective mass)'
  print >> sys.stderr, '           81 : d/d_k[.d(right)_k], partial dk derivative (d_k.d_field)'
  print >> sys.stderr, '   test_case = which test case to use (see input_template directory), optional (default=TEST_CAO-1) '
  print >> sys.stderr, '   signs = option for nonlop calculation (1:<Psi|Vnl|Psi> or 2:<g|Vnl|Psi>), optional (default=1) '
  print >> sys.stderr, '   signsdfpt= option for DFPT nonlop calc. (1:<Psi|Vnl^(1)|Psi> or 2:<g|Vnl^(1)|Psi>), optional (default=signs) '
  print >> sys.stderr, '   idir1 = direction of 1st perturbation, optional (default=1)'
  print >> sys.stderr, '   idir2 = direction of 2nd perturbation, optional (default=1) '
  print >> sys.stderr, '   iatom = index of perturbed atom (direction idir1), optional (default=1) '
  print >> sys.stderr, '   iband = compute only band iband, optional (default=all bands)'
  print >> sys.stderr, '   enl   = option for NL operator (sij: use Sij, dij: use Dij), optional (default=sij)'
  sys.exit()

test_type='$';test_case='TEST_CAO-1'
dir1=0;dir2=0;atom=0;band=0;signs=1;signsdfpt=-1;enl='sij'
for ii in range(n_arg+1)[1:]:
  arg=sys.argv[ii]
  if arg=='-c' or arg=='-choice' or arg=='--choice':
    test_type=sys.argv[ii+1]
  if arg=='-t' or arg=='-test_case' or arg=='--test_case':
    test_case=sys.argv[ii+1]
  if arg=='-i1' or arg=='-idir1' or arg=='--idir1':
    dir1=int(sys.argv[ii+1])
  if arg=='-i2' or arg=='-idir2' or arg=='--idir2':
    dir2=int(sys.argv[ii+1])
  if arg=='-ia' or arg=='-iatom' or arg=='--iatom':
    atom=int(sys.argv[ii+1])
  if arg=='-ib' or arg=='-iband' or arg=='--iband':
    band=int(sys.argv[ii+1])
  if arg=='-s' or arg=='-signs' or arg=='--signs':
    signs=int(sys.argv[ii+1])
  if arg=='-s' or arg=='-signsdfpt' or arg=='--signsdfpt':
    signsdfpt=int(sys.argv[ii+1])
  if arg=='-e' or arg=='-enl' or arg=='--enl':
    enl=sys.argv[ii+1]
if signsdfpt==-1:signsdfpt=signs

if test_type!='2'   and test_type!='3'   and \
   test_type!='5'   and test_type!='51'  and \
   test_type!='54k' and test_type!='54a' and \
   test_type!='55k' and test_type!='55s' and \
   test_type!='8'   and test_type!='81':
  print >> sys.stderr, 'Error: wrong value for choice!'
  sys.exit()
if dir2<0 or dir2>6 or dir1<0 or dir1>6:
  print >> sys.stderr, 'Error: wrong values for dir1/dir2!'
  sys.exit()
if atom<0:
  print >> sys.stderr, 'Error: wrong value for iatom!'
  sys.exit()
if band<0:
  print >> sys.stderr, 'Error: wrong value for iband!'
  sys.exit()
if signs!=1 and signs!=2:
  print >> sys.stderr, 'Error: wrong value for signs!'
  sys.exit()
if signsdfpt!=1 and signsdfpt!=2:
  print >> sys.stderr, 'Error: wrong value for signsdfpt!'
  sys.exit()
if enl!='sij' and enl!='dij':
  print >> sys.stderr, 'Error: wrong value for enl!'
  sys.exit()

if (signsdfpt==2 and (test_type=='55k' or test_type=='55s')):
  print >> sys.stderr, 'Error: signsdfpt=%s not allowed with choice=%s!' %(signsdfpt,test_type)
  sys.exit()


#Name of input dir
GENERIC_INPUT_DIR=test_case
if test_type  =='2':       # Force
  INPUT_DIR=GENERIC_INPUT_DIR+'_ATM'
elif test_type=='3':       # Stress
  INPUT_DIR=GENERIC_INPUT_DIR+'_STR'
elif test_type=='5':       # ddk
  INPUT_DIR=GENERIC_INPUT_DIR+'_K'
elif test_type=='51':      # d(right)dk
  INPUT_DIR=GENERIC_INPUT_DIR+'_K'
elif test_type=='54k':     # Effective charge
  INPUT_DIR=GENERIC_INPUT_DIR+'_K'
elif test_type=='54a':     # Effective charge
  INPUT_DIR=GENERIC_INPUT_DIR+'_ATM'
elif test_type=='55k':     # Piezo
  INPUT_DIR=GENERIC_INPUT_DIR+'_K'
elif test_type=='55s':     # Piezo
  INPUT_DIR=GENERIC_INPUT_DIR+'_STR'
elif test_type=='8':       # D2dk
  INPUT_DIR=GENERIC_INPUT_DIR+'_K'
elif test_type=='81':      # D2dk partial
  INPUT_DIR=GENERIC_INPUT_DIR+'_K'

#natom,nband, delta:
#These values should be consistent with the ABINIT input file
#Eventually read the values from data file
natom=2 ; nband=8 ; DELTA=0.0001 ; pseudos=['Ca.LDA_PW-JTH.xml']
ff=open('./input_template/'+INPUT_DIR+'/data','r')
fflines=ff.readlines()
ff.close()
for lgn in fflines:
 if lgn.find('natom=')!=-1: natom=int(lgn.split()[1])
 if lgn.find('nband=')!=-1: nband=int(lgn.split()[1])
 if lgn.find('delta=')!=-1: DELTA=float(lgn.split()[1])
 if lgn.find('pseudo=')!=-1: pseudos=(lgn.split()[1:])

if atom>natom:
  print >> sys.stderr, 'Error: iatom>natom!'
  sys.exit()
if band>nband:
  print >> sys.stderr, 'Error: iband>nband!'
  sys.exit()

#Select nonlop input arguments
dir1_list=[1];dir2_list=[1];atom_list=[0]
if test_type=='2':         # Force
  choice='choice1'
  choicedfpt='choicedfpt2'
  dir1_list=[1,2,3]
  dir2_list=[0]
  atom_list=range(natom+1)[1:]
  df_conv_factor=1.
elif test_type=='3':       # Stress
  choice='choice1'
  choicedfpt='choicedfpt3'
  dir1_list=[1,2,3,4,5,6]
  dir2_list=[0]
  atom_list=[0]
  df_conv_factor=1.
if test_type=='5':         # ddk
  choice='choice1'
  choicedfpt='choicedfpt5'
  dir1_list=[1,2,3]
  dir2_list=[0]
  atom_list=[0]
  df_conv_factor=1.
if test_type=='51':        # d(right)dk
  choice='choice1'
  choicedfpt='choicedfpt51'
  dir1_list=[1,2,3]
  dir2_list=[0]
  atom_list=[0]
  df_conv_factor=2.
elif test_type=='54k':     # Effective charge
  choice='choice2'
  choicedfpt='choicedfpt54'
  dir1_list=[1,2,3]
  dir2_list=[1,2,3]
  atom_list=range(natom+1)[1:]
  df_conv_factor=2.
elif test_type=='54a':     # Effective charge
  choice='choice51'
  choicedfpt='choicedfpt54'
  dir1_list=[1,2,3]
  dir2_list=[1,2,3]
  atom_list=range(natom+1)[1:]
  df_conv_factor=1.
elif test_type=='55k':     # Piezo
  choice='choice3'
  choicedfpt='choicedfpt55'
  dir1_list=[1,2,3,4,5,6]
  dir2_list=[1,2,3]
  atom_list=[0]
  df_conv_factor=2.
elif test_type=='55s':     # Piezo
  choice='choice5'
  choicedfpt='choicedfpt55'
  dir1_list=[1,2,3,4,5,6]
  dir2_list=[1,2,3]
  atom_list=[0]
  df_conv_factor=2.
elif test_type=='8':       # D2dk
  choice='choice5'
  choicedfpt='choicedfpt8'
  dir1_list=[1,2,3]
  dir2_list=[1,2,3]
  atom_list=[0]
  df_conv_factor=1.
elif test_type=='81':      # D2dk partial
  choice='choice51'
  choicedfpt='choicedfpt81'
  dir1_list=[1,2,3]
  dir2_list=[1,2,3]
  atom_list=[0]
  df_conv_factor=1.
#  choice='choice5'
#  df_conv_factor=2.

#If requested, overwrite default values for pert. dirs
if dir1>0: dir1_list=[dir1]
if dir2>0: dir2_list=[dir2]
if atom>0: atom_list=[atom]

#Print title
print("===========================================")
print("NONLOP TEST, CHOICE=%s, SIGNS=%d/%d, ENL=%s" % (test_type,signs,signsdfpt,enl))
sys.stdout.flush()

ab=[[1,1],[2,2],[3,3],[3,2],[3,1],[2,1]]
ba=[1,6,5,6,2,4,5,4,3] # from xy to Voigt

#Loop on perturbations
for iatom in atom_list:   # Atom (optional)
  for idir1 in dir1_list:   # Perturbation 1
   for idir2 in dir2_list:  # Perturbation 2 (optional)

      print("===========================================")

      if test_type=='2':         # Force
        idir       =0
        idirdfpt   =idir1
        input_index=3*(iatom-1)+idir1
        print ("atm=%d, atm_dir=%d" %(iatom,idir1))
      elif test_type=='3':       # Stress
        idir       =0
        idirdfpt   =idir1
        input_index=idir1
        alpha=ab[idir1-1][0];beta =ab[idir1-1][1]
        print ("str1_dir=%d, str2_dir=%d" %(alpha,beta))
      elif test_type=='5':       # ddk
        idir       =0
        idirdfpt   =idir1
        input_index=idir1
        print ("k_dir=%d" %(idir1))
      elif test_type=='51':      # d(right)dk
        idir       =0
        idirdfpt   =idir1
        input_index=idir1
        print ("k_dir=%d" %(idir1))
      elif test_type=='54k':     # Effective charge
        idir       =idir1
        idirdfpt   =3*(idir1-1)+idir2
        input_index=idir2
        print ("atm=%d, atm_dir=%d, k_dir=%d" %(iatom,idir1,idir2))
      elif test_type=='54a':     # Effective charge
        idir       =idir2
        idirdfpt   =3*(idir1-1)+idir2
        input_index=3*(iatom-1)+idir1
        print ("atm=%d, atm_dir=%d, k_dir=%d" %(iatom,idir1,idir2))
      elif test_type=='55k':     # Piezo
        idir       =idir1
        idirdfpt   =3*(idir1-1)+idir2
        input_index=idir2
        alpha=ab[idir1-1][0];beta =ab[idir1-1][1]
        print ("str1_dir=%d, str2_dir=%d, k_dir=%d" %(alpha,beta,idir2))
      elif test_type=='55s':     # Effective charge
        idir       =idir2
        idirdfpt   =3*(idir1-1)+idir2
        input_index=idir1
        alpha=ab[idir1-1][0];beta =ab[idir1-1][1]
        print ("str1_dir=%d, str2_dir=%d, k_dir=%d" %(alpha,beta,idir2))
      elif test_type=='8':       # D2dk
        idir       =idir2
        if signsdfpt==1:
          idirdfpt =ba[3*(idir1-1)+idir2-1]
        if signsdfpt==2:
          idirdfpt =3*(idir1-1)+idir2
        input_index=idir1
        print ("k1_dir=%d, k2_dir=%d" %(idir1,idir2))
      elif test_type=='81':      # D2dk partial
        idir       =idir2
        idirdfpt   =3*(idir1-1)+idir2
        input_index=idir1
        print ("k1_dir=%d, k2_dir=%d" %(idir1,idir2))
      sys.stdout.flush()

#     Create temporary files
      os.system('rm -rf config')
      os.system('mkdir config')
      os.system('touch config/'+choice)
      os.system('touch config/'+choicedfpt)
      os.system('cp -rf input_template/'+INPUT_DIR+'/inputDFk'+str(input_index)+'.in config/inputDF.in')
      os.system('touch config/idir'+str(idir))
      os.system('touch config/idirdfpt'+str(idirdfpt))
      os.system('touch config/signs'+str(signs))
      os.system('touch config/signsdfpt'+str(signsdfpt))
      if iatom>0:
        os.system('touch config/iatom'+str(iatom))
      if band>0:
        os.system('touch config/iband'+str(band))
      if enl=='sij':
         os.system('touch config/sij')
      if enl=='dij':
         os.system('touch config/dij')

#     Add "useria" flag in input file (to activate nonlop_test routine)
      if 'useria 112233' not in open('config/inputDF.in').read():
        with open('config/inputDF.in', "a") as ff:
          ff.write('\nuseria 112233  ! Activate nonlop_test routine\n')

#     Create files filesdf
      ff=open('./exec/filesdf','w')
      ff.write('config/inputDF.in\n')
      ff.write('exec/outputDF.out\n')
      ff.write('temp/nonlop_test_inp\n')
      ff.write('temp/nonlop_test_out\n')
      ff.write('temp/nonlop_test_tmp\n')
      for pseudo in pseudos:
        ff.write('pseudo/'+pseudo+'\n')
      ff.close()

#     Run ABINIT
      os.system('rm -rf ./temp/* ./exec/outputDF.out* ./exec/logdf')
      os.system(ABINIT_EXE+'< ./exec/filesdf > ./exec/logdf')
      os.system('rm -rf malloc.prc')

#     Extract relevant lines from ABINIT log
      os.system('grep TESTDFPT ./exec/logdf > ./exec/res')
      ff=open('./exec/res','r')
      fflines=ff.readlines()
      ff.close()

#     Set some data to locate results in log
      dtset_shift=nband+1
      if band>0:dtset_shift=2
      first_line=1+dtset_shift
      df1_shift =first_line
      df2_shift =first_line+  dtset_shift
      dfpt_shift=first_line+2*dtset_shift

#     Print result for each band
      band_list=range(nband)
      if band>0:band_list=[band-1]
      for iband in band_list:
        band_indx=band_list.index(iband)

        df1 =float(fflines[df1_shift +band_indx].split()[2])
        df2 =float(fflines[df2_shift +band_indx].split()[2])
        dfpt=float(fflines[dfpt_shift+band_indx].split()[2])

        df=float((df2-df1)/(2.*DELTA)/df_conv_factor)

        if abs(df)>5.e-8 or abs(dfpt)>5.e-8:
          diff=abs(df-dfpt)/abs(df)*100.
        else:
          diff=0.
        if diff>0.001:
          print ("  band=%d: diff=%15.12g perc. !!!!!" %(iband+1,diff))
        else:
          print ("  band=%d: diff=%15.12g perc." %(iband+1,diff))
      sys.stdout.flush()
