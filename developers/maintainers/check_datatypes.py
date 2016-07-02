#!/usr/bin/env python
# ======================================================================
# == Python script checking if some Fortran modules are judiciously   ==
# == used in ABINIT src files.                                        ==
# ==                                                                  ==
# == At present, the following modules are tested:                    ==
# ==   defs_datatypes                                                 ==
# ==   defs_abitypes                                                  ==
# ==   m_pawrhoij                                                     ==
# ==   m_pawcprj                                                      ==
# ==   m_pawtab                                                       ==
# ==   m_pawrad                                                       ==
# ==   m_pawang                                                       ==
# ==   m_sphharm                                                      ==
# ==   m_paral_atom                                                   ==
# ==   m_paw_an                                                       ==                                                         ==
# ==   m_paw_ij                                                       ==
# ==   m_pawfgrtab                                                    ==
# ==   m_paw_numeric                                                  ==
# ==   m_paw_finegrid                                                 ==
# ==   m_pawdij                                                       ==
# ==                                                                  ==
# ==   Usage:                                                         ==
# ==     check_datatypes [--all] [--unused] [--suggest] [--ok]        ==
# ==                     [--only module_name] [--autofix]             ==
# ==     --unused  : list unused datatypes                            ==
# ==     --suggest : list suggested changes in use statements         ==
# ==     --ok      : list correctly used datatypes                    ==
# ==     --all     : list all (default)                               ==
# ==     --only    : treat only module named "module_name"            ==
# ==     --autofix : automatically repair errors in src files         ==
# ==     Only one option in (unused, suggest, ok, all) allowed !      ==
# ==                                                                  ==
# ==                              M. Torrent - 09/2012, rev. 06/2013  ==
# ======================================================================
from __future__ import print_function

import os
import re
import sys

# ---------------------------------------------------------------------------
MODULES_LIST = ["defs_datatypes","defs_abitypes","m_pawrhoij","m_pawcprj", \
                "m_pawtab","m_pawang","m_pawrad","m_sphharm","m_paral_atom", \
                "m_paw_an","m_paw_ij","m_pawfgrtab","m_paw_numeric","m_pawdij", \
                "m_paw_finegrid"]

#Content of module "defs_datatypes"
DEFS_DATATYPES_TYPES = [  # Note: use only lowercases
"bandstructure_type",
"bcp_type",
"macro_uj_type",
"pseudopotential_gth_type",
"pseudopotential_type",
"pspheader_type"]
DEFS_DATATYPES_ROUTINES = []
DEFS_DATATYPES_FUNCTIONS= []

#Content of module "defs_abitypes"
DEFS_ABITYPES_TYPES = [  # Note: use only lowercases
"aim_dataset_type",
"anaddb_dataset_type",
"dataset_type",
"mpi_type",
"datafiles_type",
"hdr_type",
"ab_dimensions",
"macro_uj_type"]
DEFS_ABITYPES_ROUTINES = []
DEFS_ABITYPES_FUNCTIONS= []

#Content of module "m_pawrhoij"
M_PAWRHOIJ_TYPES = [  # Note: use only lowercases
"pawrhoij_type"]
M_PAWRHOIJ_ROUTINES = [  # Note: use only lowercases
"pawrhoij_alloc",
"pawrhoij_destroy",
"pawrhoij_nullify",
"pawrhoij_copy",
"pawrhoij_gather",
"pawrhoij_bcast",
"pawrhoij_redistribute",
"pawrhoij_io",
"pawrhoij_unpack",
"pawrhoij_init_unpacked",
"pawrhoij_destroy_unpacked",
"pawrhoij_mpisum_unpacked",
"symrhoij"]
M_PAWRHOIJ_FUNCTIONS = []

#Content of module "m_pawcprj"
M_PAWCPRJ_TYPES = [  # Note: use only lowercases
"pawcprj_type"]
M_PAWCPRJ_ROUTINES = [  # Note: use only lowercases
"pawcprj_alloc",
"pawcprj_destroy",
"pawcprj_nullify",
"pawcprj_set_zero",
"pawcprj_copy",
"pawcprj_axpby",
"pawcprj_zaxpby",
"pawcprj_lincom",
"pawcprj_output",
"pawcprj_diskinit_r",
"pawcprj_diskinit_w",
"pawcprj_diskskip",
"pawcprj_get",
"pawcprj_put",
"pawcprj_reorder",
"pawcprj_mpi_allgather",
"pawcprj_bcast",
"pawcprj_transpose",
"pawcprj_gather_spin",
"pawcprj_mpi_exch",
"pawcprj_mpi_send",
"pawcprj_mpi_recv",
"pawcprj_getdim"]
M_PAWCPRJ_FUNCTIONS = [  # Note: use only lowercases
"paw_overlap"]

#Content of module "m_pawtab"
M_PAWTAB_TYPES = [  # Note: use only lowercases
"wvlpaw_rholoc_type",
"wvlpaw_type",
"pawtab_type"]
M_PAWTAB_ROUTINES = [   # Note: use only lowercases
"pawtab_nullify",
"pawtab_destroy",
"pawtab_set_flags",
"pawtab_print",
"pawtab_bcast",
"wvlpaw_allocate",
"wvlpaw_destroy",
"wvlpaw_nullify",
"wvlpaw_rholoc_destroy",
"wvlpaw_rholoc_nullify"]
M_PAWTAB_FUNCTIONS= []  # Note: use only lowercases

#Content of module "m_pawang"
M_PAWANG_TYPES = [  # Note: use only lowercases
"pawang_type"]
M_PAWANG_ROUTINES = [   # Note: use only lowercases
"pawang_init",
"pawang_nullify",
"pawang_destroy",
"pawang_lsylm",
"initang",
"realgaunt",
"gauleg",
"mat_mlms2jmj",
"mat_slm2ylm"]
M_PAWANG_FUNCTIONS= [  # Note: use only lowercases
"gaunt",
"rfactorial",
"perms"]

#Content of module "m_pawrad"
M_PAWRAD_TYPES = [  # Note: use only lowercases
"pawrad_type"]
M_PAWRAD_ROUTINES = [   # Note: use only lowercases
"pawrad_init",
"pawrad_nullify",
"pawrad_destroy",
"pawrad_print",
"pawrad_isame",
"pawrad_copy",
"pawrad_deducer0",
"pawrad_bcast",
"simp_gen",
"nderiv_gen",
"nderiv_lin",
"bound_deriv",
"poisson",
"calc_slatradl"]
M_PAWRAD_FUNCTIONS= [  # Note: use only lowercases
"pawrad_ifromr"]

#Content of module "m_sphharm"
M_SPHHARM_TYPES = []  # Note: use only lowercases
M_SPHHARM_ROUTINES = [   # Note: use only lowercases
"ylmcd",
"ylm_cmplx",
"initylmr",
"ys",
"lxyz",
"slxyzs",
"plm_coeff",
"plm_d2theta",
"pl_deriv",
"mat_mlms2jmj",
"mat_slm2ylm"]
M_SPHHARM_FUNCTIONS= [  # Note: use only lowercases
"ylmc",
"ass_leg_pol",
"plm_dphi",
"plm_dtheta"]

#Content of module "m_paral_atom"
M_PARAL_ATOM_TYPES = []  # Note: use only lowercases
M_PARAL_ATOM_ROUTINES = [  # Note: use only lowercases
"get_my_natom",
"get_my_atmtab",
"free_my_atmtab"]  
M_PARAL_ATOM_FUNCTIONS = []  # Note: use only lowercases

#Content of module "m_paw_an"
M_PAW_AN_TYPES = [  # Note: use only lowercases
"paw_an_type"]
M_PAW_AN_ROUTINES = [  # Note: use only lowercases
"paw_an_init",
"paw_an_destroy",
"paw_an_nullify",
"paw_an_copy",
"paw_an_print",
"paw_an_gather",
"paw_an_redistribute",
"paw_an_reset_flags"]
M_PAW_AN_FUNCTIONS = []  # Note: use only lowercases

#Content of module "m_paw_ij"
M_PAW_IJ_TYPES = [  # Note: use only lowercases
"paw_ij_type"]
M_PAW_IJ_ROUTINES = [  # Note: use only lowercases
"paw_ij_init",
"paw_ij_destroy",
"paw_ij_nullify",
"paw_ij_copy",
"paw_ij_print",
"paw_ij_gather",
"paw_ij_redistribute",
"paw_ij_reset_flags"]
M_PAW_IJ_FUNCTIONS = []  # Note: use only lowercases

#Content of module "m_pawfgrtab"
M_PAWFGRTAB_TYPES = [  # Note: use only lowercases
"pawfgrtab_type"]
M_PAWFGRTAB_ROUTINES = [  # Note: use only lowercases
"pawfgrtab_init",
"pawfgrtab_destroy",
"pawfgrtab_nullify",
"pawfgrtab_copy",
"pawfgrtab_print",
"pawfgrtab_gather",
"pawfgrtab_redistribute"]
M_PAWFGRTAB_FUNCTIONS = []  # Note: use only lowercases

#Content of module "m_paw_numeric"
M_PAW_NUMERIC_TYPES = []  # Note: use only lowercases
M_PAW_NUMERIC_ROUTINES = [  # Note: use only lowercases
"paw_spline",
"paw_splint",
"paw_smooth",
"paw_sort_dp",
"jbessel",
"solvbes",
"jbessel_4spline"]
M_PAW_NUMERIC_FUNCTIONS = []  # Note: use only lowercases

#Content of module "m_pawdij"
M_PAWDIJ_TYPES = []  # Note: use only lowercases
M_PAWDIJ_ROUTINES = [  # Note: use only lowercases
"pawdij",
"pawdijhartree",
"pawdijxc",
"pawdijxcm",
"pawdijhat",
"pawdiju",
"pawdijso",
"pawdijexxc",
"pawdijfr",
"pawpupot",
"pawxpot",
"symdij",
"symdij_all"]
M_PAWDIJ_FUNCTIONS = []  # Note: use only lowercases

#Content of module "m_paw_finegrid"
M_PAW_FINEGRID_TYPES = []  # Note: use only lowercases
M_PAW_FINEGRID_ROUTINES = [  # Note: use only lowercases
"pawgylm",
"pawexpiqr"]
M_PAW_FINEGRID_FUNCTIONS = []  # Note: use only lowercases

#Global lists
TYPES_LIST    = [DEFS_DATATYPES_TYPES,DEFS_ABITYPES_TYPES, \
                 M_PAWRHOIJ_TYPES,M_PAWCPRJ_TYPES, \
                 M_PAWTAB_TYPES,M_PAWANG_TYPES,M_PAWRAD_TYPES, \
                 M_SPHHARM_TYPES,M_PARAL_ATOM_TYPES, \
                 M_PAW_AN_TYPES,M_PAW_IJ_TYPES,M_PAWFGRTAB_TYPES, \
                 M_PAW_NUMERIC_TYPES,M_PAWDIJ_TYPES,M_PAW_FINEGRID_TYPES]
ROUTINES_LIST = [DEFS_DATATYPES_ROUTINES,DEFS_ABITYPES_ROUTINES, \
                 M_PAWRHOIJ_ROUTINES,M_PAWCPRJ_ROUTINES, \
                 M_PAWTAB_ROUTINES,M_PAWANG_ROUTINES,M_PAWRAD_ROUTINES, \
                 M_SPHHARM_ROUTINES,M_PARAL_ATOM_ROUTINES, \
                 M_PAW_AN_ROUTINES,M_PAW_IJ_ROUTINES,M_PAWFGRTAB_ROUTINES, \
                 M_PAW_NUMERIC_ROUTINES,M_PAWDIJ_ROUTINES,M_PAW_FINEGRID_ROUTINES]
FUNCTIONS_LIST= [DEFS_DATATYPES_FUNCTIONS,DEFS_ABITYPES_FUNCTIONS, \
                 M_PAWRHOIJ_FUNCTIONS,M_PAWCPRJ_FUNCTIONS, \
                 M_PAWTAB_FUNCTIONS,M_PAWANG_FUNCTIONS,M_PAWRAD_FUNCTIONS, \
                 M_SPHHARM_FUNCTIONS,M_PARAL_ATOM_FUNCTIONS, \
                 M_PAW_AN_FUNCTIONS,M_PAW_IJ_FUNCTIONS,M_PAWFGRTAB_FUNCTIONS, \
                 M_PAW_NUMERIC_FUNCTIONS,M_PAWDIJ_FUNCTIONS,M_PAW_FINEGRID_FUNCTIONS]

#Lines beginning with the following statement will not be considered
COMMENT = "!"

#Files beginning with the following statement will not be considered
IGNORES_FILES = ["interfaces_"]

# ---------------------------------------------------------------------------

if (len(sys.argv)==1) or ("--help" in sys.argv):
  print ("Usage:")
  print ("  check_datatypes [--all] [--unused] [--suggest] [--ok]")
  print ("                  [--only module_name] [--autofix]")
  print ("  --unused  : list unused datatypes")
  print ("  --suggest : list suggested changes in use statements")
  print ("  --ok      : list correctly used datatypes")
  print ("  --all     : list all (default)")
  print ("  --only    : treat only module named module_name")
  print ("  --autofix : automatically repair errors in src files")
  print ("  Only one option in (unused, suggest, ok, all) allowed !")
  exit()

noselection=((not "--unused" in sys.argv) and (not "--suggest" in sys.argv) and \
             (not "--ok" in sys.argv))
do_unused =(("--all" in sys.argv) or ("--unused"  in sys.argv) or (noselection))
do_suggest=(("--all" in sys.argv) or ("--suggest" in sys.argv) or (noselection))
do_ok     =(("--all" in sys.argv) or ("--ok"      in sys.argv) or (noselection))
autofix   =("--autofix" in sys.argv)
if "--only" in sys.argv:
  if len(sys.argv)>sys.argv.index("--only"):
    modules_restricted_list=[sys.argv[sys.argv.index("--only")+1]]
  else:
    modules_restricted_list=[]
else:
  modules_restricted_list=MODULES_LIST

# ---------------------------------------------------------------------------

print()
print( '---------------------------------------------------------------------')
print( ' Looking for "use module" statements in ABINIT Fortran source files. ')
print( '  - check if they are useful                                         ')
print( '  - list missing "only" statements                                   ')
print( '---------------------------------------------------------------------\n')

re_srcfile = re.compile("\.([Ff]|[Ff]90)$")
file_total_count=0

#Loop over files in src folder
for (root, dirs, files) in os.walk("./src"):
  for src in files:
    if (re_srcfile.search(src)):
      file_total_count+=1
      filename=os.path.join(root,src)
      with open(filename, "r") as fh:
        src_data = fh.readlines()

#     Ignore some files
      ignore=False
      for item in IGNORES_FILES:
        ignore=(src.find(item)==0)
      if not ignore:

#       Loop over modules
        for module in MODULES_LIST:
          if module in modules_restricted_list:
            mod_index=MODULES_LIST.index(module)
            module_found=False

#           Loop over lines in file (searching use of module)
            for line in src_data:
              if line.find(COMMENT) != 0:
                line_lower=re.sub(" ","",line.lower()) # Lower case, no space

#               Check use of module
                if line_lower.find("use"+module) != -1 and line_lower.find("use"+module+"_") == -1:
                  module_found=True
                  line_use_module=line_lower
                  line_index=src_data.index(line)
                  for i in range(3):
                    if line_use_module.find("&")==len(line_use_module)-2:
                      for line2 in src_data[1+line_index:]:
                        line2_lower=re.sub(" ","",line2.lower())
                        if line2_lower.find("&")==0:
                          line_use_module=line_use_module[:line_use_module.find("&")-1]+line2_lower[1:]
                          line_index=src_data.index(line2)
                          break
                  break

#           Continue if module has been found
            if module_found:
              types_list_mod=[];routines_list_mod=[];functions_list_mod=[]
              sugg_only_string="";missing_count=0;module_unused=False
              module_stat_found=False;contains_stat_found=False
              subroutine_count=0;function_count=0

#             Loop over lines in file (again)
              for line in src_data:
                if line.find(COMMENT) != 0:
                  line_lower=re.sub(" ","",line.lower()) # Lower case, no space

#                 Looking for datatypes in line
                  for type in TYPES_LIST[mod_index]:
                    if line_lower.find("type("+type+")") != -1:
                      if not type in types_list_mod: types_list_mod.append(type)

#                 Looking for routines in line
                  for routine in ROUTINES_LIST[mod_index]:
                    if line_lower.find("call"+routine+"(") != -1:
                      if not routine in routines_list_mod: routines_list_mod.append(routine)

#                 Looking for functions in line
                  for function in FUNCTIONS_LIST[mod_index]:
                    if line_lower.find(function+"(") != -1:
                      if not function in functions_list_mod: functions_list_mod.append(function)

#                 Looking for multiple "subroutines"
                  if not module_stat_found: module_stat_found=(line_lower.find("module") != -1)
                  if not contains_stat_found: contains_stat_found=(line_lower.find("contains") != -1)
                  if line_lower.find("endsubroutine")!=-1:subroutine_count+=1
                  if line_lower.find("endfunction")!=-1:function_count+=1

              multiple_sub=((not module_stat_found and not contains_stat_found) and \
                            (subroutine_count>1 or function_count>1))
              list_all=types_list_mod+routines_list_mod+functions_list_mod

#             First case: module not used - print a warning
              if (len(list_all))==0:
                module_unused=True
                if do_unused and not autofix:
                  print ("File %s: module %s present but not used !" % (filename,module))

#             Second case: module used - suggest "only" statement
              else:
                len_cur=0
                for item in list_all:
                  if line_use_module.find(item) == -1: missing_count+=1
                  if sugg_only_string != "":
                    sugg_only_string+=", ";len_cur=len_cur+1
                  if len_cur>75:
                    if list_all.index(item)!=-1:
                      sugg_only_string+=" &\n&"+' '.ljust(13+len(module))
                      len_cur=0
                  sugg_only_string+=item;len_cur=len_cur+len(item)
                if missing_count==0:
                  if do_ok:
                    print ("File %s: module %s correctly used !" % (filename,module))
                else:
                  if do_suggest and not autofix:
                    print ("File %s: module %s, suggested use:" % (filename,module))
                    print ("  * use "+module+", only : "+sugg_only_string)
                    if multiple_sub:
                      print ("  * WARNING: several subroutines/functions in file !")

#             === AUTOFIX ===
              if autofix:
                if (module_unused and do_unused) or (missing_count>0 and do_suggest):
                  print ("==> FIXING file %s" %(filename))
                  if module_unused and do_unused:
                    print ('    ELIMINATING module %s !' % (module))
                  if missing_count>0 and do_suggest:
                    print ('    REWRITING "use module" line for %s !' % (module))
                    if multiple_sub:
                      print ('    WARNING: several subroutines/functions in file (check it MANUALLY) !')
#                 Open a temporary file for writing
                  ErrorEncountered=False
                  filenametmp=filename+'.tmp'
                  try: filout=open(filenametmp,'w')
                  except:
                    print ('File %s, error: could not open tmp file !' % (filename))
                    ErrorEncountered=True
                  if not ErrorEncountered:
#                   Loop over lines in the file
                    for line in src_data:
                      line_lower=re.sub(" ","",line.lower()) # Lower case, no space
                      do_write_line=1
                      if line_lower.find("use"+module) != -1 and line.find(COMMENT) != 0:
                        if module_unused and do_unused: do_write_line=0
                        if missing_count>0 and do_suggest: do_write_line=2
                      if do_write_line>0:
                        newline=line
                        if do_write_line==2: newline=" use "+module+", only : "+sugg_only_string+"\n"
                        try: filout.write(newline)
                        except:
                          print ('File %s, error: could not write into tmp file !' % (filename))
                          ErrorEncountered=True;break
                  filout.close()
#                 Replace current file by temporary file
                  if not ErrorEncountered:
                    try: 
                        os.system('mv -f '+filenametmp+' '+filename)
                    except: 
                        print ('File %s, error: could not move tmp file !' % (filename))

#Final printing
print ("--------------------------------------")
print ("Done !, %s files have been explored." % (file_total_count))
