#! /usr/bin/python

#
#    Copyright (C) 2003-2020 ABINIT group
#
#    Written by Matthieu Verstraete in python (compatible v1.9).
#    This is free software, and you are welcome to redistribute it
#    under certain conditions (GNU General Public License,
#    see ~abinit/COPYING or http://www.gnu.org/copyleft/gpl.txt).
#
#    ABINIT is a project of the Universite Catholique de Louvain,
#    Corning Inc. and other collaborators, see ~abinit/doc/developers/contributors.txt.
#    Please read ~abinit/doc/biblio/generated_files/bib_acknow.html for suggested
#    acknowledgments of the ABINIT effort.
#
#    For more information, see https://www.abinit.org .
#
#    Input is files from ABINIT: ground state or phonon runs, and from anaddb
#    for the Gamma point at least.
#
#    TODO:
#    * Does not yet understand non-orthogonal basis vectors.
#    * LO/TO splitting is not taken into account yet for the calculation
#      of the characters in anaddb so the characters are for the TO modes
#      only. The splitting can be identified afterwards.
#    * Add the primes and plus/minus to the irrep symbols according to the
#      2-axis and mirror symmetry character of each irrep.
#    * Irreps are numbered according to order returned by Bilbao server,
#      which may not be canonical order.
#    * If the symmetries are one day output to the anaddb files, only
#      one input file will be needed (ie parse it for syms, instead of 
#      "abifilename").
#


#  reg expss
import re
import string
import sys
#import fileinput
import math
import xml.dom.minidom
from tokenize_file import *
import urllib


elemabbrev = ['XX','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Un']
ha2cmm1 = 219474.63

def get_sym_elem(chardom):
	"""Retrieve symmetry elements from XML file"""
	isym=0
	sym_elem = []
	sym_tnons = []
	sym_elem_nodes = chardom.getElementsByTagName('genpos')
#	print " number of sym elements ", len(sym_elem_nodes)
	for sym_elem_node in sym_elem_nodes:
		#print sym_elem_node
		sym_elem.append([])
		sym_tnons.append([])
		rows = sym_elem_node.getElementsByTagName('row')
		for row in rows:
			try:
				data = row.childNodes[0].data
				#print data
				tmp = re.split ('[ \n\t\[\]]*',data)
				tmp.pop()
				sym_tnons[isym].append(tmp.pop())
				tmp.pop(0)
				#print tmp
				sym_elem[isym].extend(map(string.atoi,tmp))
			except AttributeError:
				print "value error in get_sym_elem" 
				raise
				
		isym = isym+1
	
#	print sym_elem
	return sym_elem,sym_tnons

def get_irreps_chars(chardom,nsym):
	"""Retrieve characters for sym elements from XML file"""
	irreps_node = chardom.getElementsByTagName('irreps')
	if (len(irreps_node) > 1):
		print "Error : more than one irreps field in XML"
		raise ValueError
	try :
		nirreps = irreps_node[0].getAttribute('num')
	except ValueError:
		print "no nirreps field in irreps tag"
		raise
#	print " number of irreps ", nirreps
	repr_nodes = irreps_node[0].getElementsByTagName('repr')
	irrep_character = []
	irrep_name = []
	iirrep=0
	for repr_node in repr_nodes:
		irrep_character.append([])
		#print repr_node
		irrep_name.append(repr_node.getAttribute('name'))
		dim = string.atoi(repr_node.getAttribute('dim'))
		for mat in repr_node.getElementsByTagName('mat'):
			idim=0
			character = 0.0
			for row in mat.getElementsByTagName('row'):
				try:
					data = row.childNodes[0].data
					#print data
					tmp = re.split ('[ \n\t\[\]\(\)\,]*',data)
					tmp.pop()
					tmp.pop(0)
					#print tmp
                                        char = float(tmp[2*idim])*math.cos(math.pi*float(tmp[2*idim+1])/180.0)
					#print char
					character = character + char
				except AttributeError:
					print "value error in get_irreps_chars" 
					raise
				
				idim=idim+1
                        irrep_character[iirrep].append(int(character))
		iirrep=iirrep+1
	
	print irrep_name,irrep_character
	return irrep_name,irrep_character

def get_abinit_syms(tokens):
	"""Retrieve symmetry elements from abinit file tokenized previously"""
	
	ind_nsym = tokens.index('nsym')
	nsym = int (tokens[ind_nsym+2])
#	print "nsym = ", nsym
	ind_symrel = tokens.index('symrel')
	symrel = []
	ii = 0
	for isym in range(nsym):
		symrel.append(map(int,tokens[ind_symrel+ii+1:ind_symrel+ii+10]))
		ii = ii+9
	tnons = []
	try:
		ind_tnons = tokens.index('tnons')
		ii = 0
		for itnons in range(nsym):
			tnons.append(map(float,tokens[ind_tnons+ii+1:ind_tnons+ii+3]))
			ii = ii+3
	except ValueError:
		for itnons in range(nsym):
			tnons.append((0.0,0.0,0.0))
		
#	print "symrel = ", symrel
	return symrel

def get_abinit_group(tokens):
	"""Retrieve symmetry group from abinit output file"""
	groupnum = 0
	try:
		ind_spgroup = tokens.index('spgroup')
		groupnum = int(tokens[ind_spgroup+1])
	except ValueError:
		ind_spgroup = tokens.index('space')
		for ii in range(10):
			if tokens[ind_spgroup+1] == 'group':
				while not re.compile('\(\#[0-9]*\)').match(tokens[ind_spgroup]):
					ind_spgroup = ind_spgroup+1
				groupnum = int(re.sub('[\(\)\#;]*','',tokens[ind_spgroup]))
				break
			else:
				tokens.pop(ind_spgroup)
#	print groupnum
	return groupnum
#		if not re.compile('degenerate').match(tokens[ind1]):

			
def get_abinit_chars(tokens):
	"""Retrieve characters and frequencies for modes from anaddb output file"""

#  find the correct occurence of Phonon which is followed by the gamma point freqs
	for ii in range(tokens.count('Phonon')):
		ind1 = tokens.index('Phonon')
		if tokens[ind1+1] != 'wavevector' or tokens[ind1+2] != '(reduced' :
			tokens.pop(ind1)
		else:
#  if gamma point
			if (float(tokens[ind1+5]) == 0.0 and \
				float(tokens[ind1+6]) == 0.0 and \
				float(tokens[ind1+7]) == 0.0 ) :
				tokens.pop(ind1)
			break
	abi_freq = []
	ind1=ind1+12
#	print tokens[ind1],tokens[ind1+1]
#	while re.compile('[-]?[0-9].[0-9]*[Ee]?[-]?[0-9]*').match(tokens[ind1]):
	while re.compile('[-]?[0-9]+.[0-9]*[Ee]?[-+]?[0-9]*').match(tokens[ind1]):
		abi_freq.append(float(tokens[ind1]))
		#print abi_freq
		ind1=ind1+1

	for ii in range(tokens.count('Analysis')):
		ind1 = tokens.index('Analysis')
		if tokens[ind1+2] != 'degeneracies' or tokens[ind1+4] != 'characters':
			tokens.pop(ind1)
		else:
			break
	while not re.compile('Symmetry').match(tokens[ind1]):
		ind1=ind1+1
	imode=0
	abi_characters = []
	while not re.compile('=[=]*').match(tokens[ind1]):
		abi_characters.append([])
		while not re.compile('[0-9]').match(tokens[ind1]):
			ind1=ind1+1
		imode_file = int(tokens[ind1])-1
		if (imode_file != imode):
			print "Probably error : modes don't match ", imode, imode_file
		ind1=ind1+1
		#print "tokens now ", tokens[ind1:ind1+20]
		#print "abi_characters ", abi_characters
		degen = 1
		if not re.compile('degenerate').match(tokens[ind1]):
			while re.compile('[-]?[0-9]+\.[0-9]+').match(tokens[ind1]):
				charlist = int(float(tokens[ind1]))
				abi_characters[imode].append(charlist)
				ind1 = ind1+1
		else:
#  Find degeneracy of present mode
			while not  re.compile('[0-9]+').match(tokens[ind1]):
				ind1=ind1+1
			if re.compile('to').match(tokens[ind1+1]):
				ind1=ind1+2
			jmode = int(tokens[ind1])
			degen = jmode-imode
			
			ind1=ind1+1
			while not re.compile('[-]?[0-9]+\.[0-9]+').match(tokens[ind1]):
				ind1=ind1+1
			while re.compile('[-]?[0-9]+\.[0-9]+').match(tokens[ind1]):
				charlist = int(float(tokens[ind1]))
				abi_characters[imode].append(charlist)
				ind1 = ind1+1
			
			
# Duplicate characters for degenerate modes
		for ii in range(degen-1):
			abi_characters.append(abi_characters[imode])
		imode=imode+degen

	#print "abi_characters = ", abi_characters,abi_freq
	return abi_characters,abi_freq

	
def get_sym_corresp(abinit_sym_elem,bilbao_sym_elem):
	"""determine correspondence between Bilbao and abinit symmetry lists"""
	abi_to_bilbao = []
	nsym = len(abinit_sym_elem)
	print " nsyms = ", nsym, len(bilbao_sym_elem)
	print "abinit_sym_elem = "
	for symelem in abinit_sym_elem:
		print symelem
	print "bilbao_sym_elem = "
	for symelem in bilbao_sym_elem:
		print symelem
	if nsym != len(bilbao_sym_elem):
		print "Error: different number of syms in abinit/bilbao"
		raise ValueError
	for iabi in range(nsym):
		for ibil in range(nsym):
			found = 1
			for ielem in range(9):
				if abinit_sym_elem[iabi] != bilbao_sym_elem[ibil]:
					found=0
					break
			if found == 1:
				abi_to_bilbao.append(ibil)
				break
	if len(abi_to_bilbao) != nsym:
		print "Error : not all correspondences found"
		print "abi_to_bilbao = "
		print abi_to_bilbao
		raise ValueError
#	print abi_to_bilbao
	return abi_to_bilbao

def get_bilbao_symbols(bilbao_sym_elem,bilbao_irreps_chars):
	"""Deteremine irreps names from symmetry considerations"""
	
	iident = find_ident(bilbao_sym_elem)
	iinv = -1
	try:
		iinv = find_inv(bilbao_sym_elem)
		has_inv = 1
	except ValueError:
		print "Inv not found"
		has_inv = 0
	irot,rot_order = find_princ_rotation(bilbao_sym_elem,has_inv,iinv)
#	print " irot,rot_order = ",irot,rot_order
	if has_inv==1:
		print "ident,inv,irot,rot_order = ", iident,iinv,irot,rot_order
	elif has_inv==0:
		print "ident,irot,rot_order = ", iident,irot,rot_order

	ia = 1
	ib = 1
	ie = 1
	it = 1

	bilbao_symbols=[]
	for iirrep in range(len(bilbao_irreps_chars)):
		if bilbao_irreps_chars[iirrep][iident] == 1:
			if bilbao_irreps_chars[iirrep][irot] == 1:
				irrepletter = 'A'
				irrepnum = ia
				ia = ia+1
			elif bilbao_irreps_chars[iirrep][irot] == -1:
				irrepletter = 'B'
				irrepnum = ib
				ib = ib+1
			else:
				print "couldn't classify the irrep into A or B"
				raise ValueError
		elif bilbao_irreps_chars[iirrep][iident] == 2:
			irrepletter = 'E'
			irrepnum = ie
			ie = ie+1
		elif bilbao_irreps_chars[iirrep][iident] == 3:
			irrepletter = 'T'
			irrepnum = it
			it = it+1
		else:
			print "couldn't classify the irrep into A,B,E, or T"
			raise ValueError
		
		geradeletter = ' '
		if has_inv==1:
			if bilbao_irreps_chars[iirrep][iinv] > 0:
				geradeletter = 'g'
			elif bilbao_irreps_chars[iirrep][iinv] < 0:
				geradeletter = 'u'
			irrepnum = int((irrepnum+1.0)/2.0)
			
		print "%1c%1d%1c" % (irrepletter,irrepnum,geradeletter)
		bilbao_symbols.append("%1c%1d%1c" % (irrepletter,irrepnum,geradeletter))

	return bilbao_symbols

def find_ident(sym_elem):
	"""find the index of the identity operation in the symmetries"""
	for isym in range(len(sym_elem)):
		if sym_elem[isym] == [1,0,0,0,1,0,0,0,1]:
			return isym
	print "Error : ident not found"
	raise ValueError

def find_inv(sym_elem):
	"""find the index of the inverse operation in the symmetries, if it exists"""
	for isym in range(len(sym_elem)):
		if sym_elem[isym] == [-1,0,0,0,-1,0,0,0,-1]:
			return isym
#	print "Error : inv not found"
	raise ValueError

def find_princ_rotation(sym_elem,has_inv,iinv):
	"""Find the principal rotation of the group, as the highest order member which doesn't give -I"""
	morder = 0
	msym = 0
	for isym in range(len(sym_elem)):
		is_Sop = 0
		#if has_inv==1:
		#	if 
		# suppose that order is at most 6 for elements of sym
		tmpsym = [1,0,0,0,1,0,0,0,1]
		for iorder in range(1,7):
#			print "iorder= ", iorder
			tmpsym = sym_elem_prod(tmpsym,sym_elem[isym])
			if tmpsym == [-1,0,0,0,-1,0,0,0,-1]:
				is_Sop = 1
				break
			if tmpsym == [1,0,0,0,1,0,0,0,1]:
				break
		# if we found an S operation (axis+Inverse) dont count it
		#   should have a corresponding C operation
		#   the danger is that order(S) = 2*order(C) and we want the C
		#   to deterine the principal rotation
		if is_Sop == 1:
			break
		if iorder > morder:
			msym = isym
			morder = iorder
	return msym,morder

def sym_elem_prod (s1,s2):
	""" calculate (non commuting) product of s2*s1"""
	sym = []
	for ii in range(3):
		for jj in range(3):
			ss = 0
			for kk in range(3):
				ss = ss + s2[3*ii+kk]*s1[3*kk+jj]
			sym.append(ss)
	return sym

def print_abi_irreps (abi_to_bilbao,abinit_chars,abinit_freq,irreps_names,irreps_chars):
	"""Print out irreps for each abinit mode"""
	
	order = len(irreps_chars[0])
#	print "abinit_chars = ", abinit_chars
	for ichar in range(len(abinit_chars)):
		freqcmm1 = abinit_freq[ichar]*ha2cmm1
		ichar_eff=ichar+1
		print "Abinit mode # ", "%2d = " % ichar_eff, " %15.2f cm-1 " % freqcmm1,
		for iirrep in range(len(irreps_chars)):
			sum = 0
			for isym in range(len(irreps_chars[iirrep])):
				sum = sum + irreps_chars[iirrep][abi_to_bilbao[isym]]*abinit_chars[ichar][isym]
			if sum != 0:
				print int(sum/order), "x", irreps_names[iirrep], '(Irrep num. ', iirrep, ')',
                                # The following should never happen
                                if (sum % order != 0):
                                        print '(Sum ', sum, ', order ', order, ')',
		print
	

if __name__=='__main__':
#	for ielem in range(len(elemabbrev)):
#		print ielem, elemabbrev[ielem]

#
	if len(sys.argv) < 2:
		print "Error : need two filename arguments : abinit output file and anaddb output file"
		raise ValueError
	else:
#		charfilename = sys.argv[1]
#		groupnum = int(sys.argv[1])
		abifilename = sys.argv[1]
		anaddbfilename = sys.argv[2]

#
#  get group and symmetry operations from abinit output file
#
	abitokens=tokenize_file(abifilename)
	#print abitokens
	
	groupnum = get_abinit_group(abitokens)
	print groupnum
	
	abinit_sym_elem = get_abinit_syms(abitokens)
	
#
#  get characters of different modes from anaddb output file. Get the information
#    only for the Gamma point.
#
	anaddbtokens=tokenize_file_keep_comments(anaddbfilename)
	#print anaddbtokens

	abinit_chars,abinit_freq = get_abinit_chars(anaddbtokens)
	
#
#   open url and request xml input for characters and symops of 
#     point group "groupnum"
#
	print "http://www.cryst.ehu.es"+\
		"/cgi-bin/cryst/xml/nph-get_doc?p=irreps&g=%d&k=0,0,0,GM,p\n" % groupnum
	try:
		fhan = urllib.urlopen("http://www.cryst.ehu.es"+\
			"/cgi-bin/cryst/xml/nph-get_doc?p=irreps&g=%d&k=0,0,0,GM,p&gor=2\n" % groupnum)
	except IOError:
		print "Error opening http://www.cryst.ehu.es"
		raise
	
	bilbao_xml_lines = fhan.readlines()
	XML = ''
	for line in bilbao_xml_lines:
		XML = XML + string.strip(line)
#		print bilbao_xml_lines
	

#
#   parse the XML input obtained from www.cryst.ehu.es Bilbao server
#    use the minimalistic minidom. No dtd support in minidom, although 
#    one exists on the server.
#
	chardom = xml.dom.minidom.parseString(XML)
	#print chardom.documentElement.tagName
	if (chardom.documentElement.hasAttributes()):
		print chardom.documentElement.attributes
#	else:
#		print "No attributes here!"


#
#  extract symmetry elements from Bilbao xml,
#
	bilbao_sym_elem,bilbao_sym_tnons = get_sym_elem(chardom)

#
#  extract irreps from Bilbao xml,
#
	bilbao_irreps_names,bilbao_irreps_chars=get_irreps_chars(chardom,len(bilbao_sym_elem))

#
#  get correspondance between abinit and Bilbao symmetry elements
#
	abi_to_bilbao = get_sym_corresp(abinit_sym_elem,bilbao_sym_elem)

#
#  generate names of irreps A1u ... from the characters of each irrep for -I, main rotation axis...
#
	bilbao_irreps_names = get_bilbao_symbols(bilbao_sym_elem,bilbao_irreps_chars)

#
#  print out abinit frequencies and irreps for each phonon mode.
#
	print_abi_irreps (abi_to_bilbao,abinit_chars,abinit_freq,bilbao_irreps_names,bilbao_irreps_chars)
	

