#*****************************************************************************
#                                                                            *
#       Author : J.Jansen                                                    *
#       Version : 2.4                                                        *
#       Date : 14 February 2007                                              *
#       Purpose : Abinit make file for OpenVMS                               *
#                                                                            *
#*****************************************************************************
#

MAINS=aim,anaddb,cut3d,mrgddb,lwf,conducti

# Archives for Src* directories
AR_DRIVE=[.src.95_drive]lib95_drive.olb
AR_CUT3D=[.src.83_cut3d]lib83_cut3d.olb
AR_SEQ=[.src.79_seqpar_mpi]lib18abinis.olb
AR_PAR=[.src.79_seqpar_mpi]lib18abinip.olb
AR_DDB=[.src.77_ddb]lib77_ddb.olb
AR_LWF=[.src.77_lwf]lib77_lwf.olb
AR_SUSCEP=[.src.77_suscep]lib77_suscep.olb
AR_RESPONSE=[.src.72_response]lib72_response.olb
AR_COMMON=[.src.67_common]lib67_common.olb
AR_IOWFDENPOT=[.src.62_iowfdenpot]lib62_iowfdenpot.olb
AR_WFS=[.src.14wfs]lib14wfs.olb
AR_WVL_WFS=[.src.62_wvl_wfs]lib62_wvl_wfs.olb
AR_RECURSION=[.src.67_recursion]lib67_recursion.olb
AR_RSPRC=[.src.68_rsprc]lib68_rsprc.olb
AR_GW=[.src.68_gw]lib68_gw.olb
AR_IOVARS=[.src.57_iovars]lib57_iovars.olb
AR_IONETCDF=[.src.59_ionetcdf]lib59_ionetcdf.olb
AR_PAW=[.src.13paw]lib13paw.olb
AR_RECIPSPACE=[.src.56_recipspace]lib56_recipspace.olb
AR_XC=[.src.56_xc]lib56_xc.olb
AR_XML=[.src.47_xml]lib47_xml.olb
AR_BADER=[.src.62_bader]lib62_bader.olb
AR_ABIUTIL=[.src.55_abiutil]lib55_abiutil.olb
AR_NONLOCAL=[.src.13nonlocal]lib13nonlocal.olb
AR_FFTS=[.src.53_ffts]lib53_ffts.olb
AR_PSP=[.src.13psp]lib13psp.olb
AR_GEOMETRY=[.src.42_geometry]lib42_geometry.olb
AR_GEOMOPTIM=[.src.45_geomoptim]lib45_geomoptim.olb
AR_PARSER=[.src.42_parser]lib42_parser.olb
AR_SPACEPAR=[.src.54_spacepar]lib54_spacepar.olb
AR_MGMPI=[.src.51_manage_mpi]lib51_manage_mpi.olb
AR_CONTRACT=[.src.32_contract]lib32_contract.olb
AR_CG=[.src.62_cg_noabirule]lib62_cg_noabirule.olb
AR_HIDEMPI=[.src.12_hide_mpi]lib12_hide_mpi.olb
AR_UTIL=[.src.32_util]lib32_util.olb
AR_HIDEWRITE=[.src.14_hidewrite]lib14_hidewrite.olb
AR_HIDELEAVE=[.src.16_hideleave]lib16_hideleave.olb
AR_TIMING=[.src.18_timing]lib18_timing.olb
AR_DEFS=[.src.10_defs]libdefs.olb
AR_OCCEIG=[.src.62_occeig]lib62_occeig.olb
AR_NLSTRAIN=[.src.42_nlstrain]lib42_nlstrain.olb
AR_IOMPI=[.src.56_io_mpi]lib56_io_mpi.olb
AR_POISSON=[.src.62_poisson]lib62_poisson.olb

AR_NUMERIC=[.lib.numeric]libnumeric.olb
AR_NUMERICF90=[.lib.numericf90]libnumericf90.olb
AR_FFTNEW=[.src.52_fft_mpi_noabirule]lib52_fft_mpi_noabirule.olb

allseq : version,crea_descrip_mms.exe,vms_prepare_input.exe,prepare_input
	$(MMS) numeric,abinit,$(MAINS)

crea_descrip_mms.exe : crea_descrip_mms.obj
	link crea_descrip_mms

crea_descrip_mms.obj : [.vms]crea_descrip_mms.f90
	cc/define=(VMS)/comment=as_is/prep=[.vms]crea_descrip_mms.f_\
	[.vms]crea_descrip_mms.f90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	[.vms]crea_descrip_mms.f_
	delete [.vms]crea_descrip_mms.f_;*

vms_prepare_input.exe : vms_prepare_input.obj
	link vms_prepare_input

vms_prepare_input.obj : [.vms]vms_prepare_input.f90
	cc/define=(VMS)/comment=as_is/prep=[.vms]vms_prepare_input.f_\
	[.vms]vms_prepare_input.f90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	[.vms]vms_prepare_input.f_
	delete [.vms]vms_prepare_input.f_;*

numeric :
	$(MMS) [.lib.numeric]descrip.mms
	set default [.lib.numeric]
	$(MMS)
	set default [--]

[.src.52_fft_mpi_noabirule]lib52_fft_mpi_noabirule.olb :
	$(MMS) [.src.52_fft_mpi_noabirule]descrip.mms
	set default [.src.52_fft_mpi_noabirule]
	$(MMS)
	set default [--]

[.lib.numericf90]libnumericf90.olb :
	$(MMS) [.lib.numericf90]descrip.mms
	set default [.lib.numericf90]
	$(MMS)
	set default [--]

DEP_ABINIT=$(AR_DRIVE) $(AR_SEQ) $(AR_SUSCEP) \
 $(AR_RESPONSE) $(AR_COMMON) $(AR_IOWFDENPOT) $(AR_WFS) $(AR_WVL_WFS) \
 $(AR_GW) $(AR_RECURSION) $(AR_RSPRC)  \
 $(AR_IOVARS) $(AR_IOMPI) $(AR_IONETCDF) $(AR_PAW) $(AR_RECIPSPACE) $(AR_XC) \
 $(AR_XML) $(AR_NONLOCAL) $(AR_FFTS) $(AR_PSP) $(AR_GEOMETRY) $(AR_PARSER) \
 $(AR_SPACEPAR) \
 $(AR_UTIL) $(AR_CONTRACT) $(AR_MGMPI) $(AR_BASIS) $(AR_DEFS) $(AR_FFTNEW)\
 $(AR_NUMERIC) $(AR_NUMERICF90) $(AR_GEOMOPTIM) $(AR_OCCEIG) $(AR_NLSTRAIN)\
 $(AR_POISSON) $(AR_CG) $(AR_HIDEMPI)

abinit : abinit.exe
	write sys$output " abinit.exe has been made "
	write sys$output " "

abinit.exe : [.src.98_main]abinit.obj $(DEP_ABINIT)
	write sys$output " abinit.exe will be made "
	write sys$output " "
	open/write optf abinit.opt
	write optf "$(AR_DRIVE)/lib"
	write optf "$(AR_SEQ)/lib"
	write optf "$(AR_SUSCEP)/lib"
	write optf "$(AR_RESPONSE)/lib"
	write optf "$(AR_GEOMOPTIM)/lib"
	write optf "$(AR_COMMON)/lib"
	write optf "$(AR_GW)/lib"
	write optf "$(AR_RECURSION)/lib"
	write optf "$(AR_RSPRC)/lib"
	write optf "$(AR_IOWFDENPOT)/lib"
	write optf "$(AR_WFS)/lib"
	write optf "$(AR_WVL_WFS)/lib"
	write optf "$(AR_OCCEIG)/lib"
	write optf "$(AR_IOVARS)/lib"
	write optf "$(AR_IOMPI)/lib"
	write optf "$(AR_IONETCDF)/lib"
	write optf "$(AR_PAW)/lib"
	write optf "$(AR_RECIPSPACE)/lib"
	write optf "$(AR_XC)/lib"
	write optf "$(AR_XML)/lib"
	write optf "$(AR_NONLOCAL)/lib"
	write optf "$(AR_FFTS)/lib"
	write optf "$(AR_PSP)/lib"
	write optf "$(AR_NLSTRAIN)/lib"
	write optf "$(AR_POISSON)/lib"
	write optf "$(AR_GEOMETRY)/lib"
	write optf "$(AR_PARSER)/lib"
	write optf "$(AR_SPACEPAR)/lib"
	write optf "$(AR_UTIL)/lib"
	write optf "$(AR_CONTRACT)/lib"
	write optf "$(AR_MGMPI)/lib"
	write optf "$(AR_FFTNEW)/lib"
	write optf "$(AR_CG)/lib"
	write optf "$(AR_HIDEMPI)/lib"
	write optf "$(AR_BASIS)/lib"
	write optf "$(AR_DEFS)/lib"
	write optf "$(AR_NUMERIC)/lib"
	write optf "$(AR_NUMERICF90)/lib"
	write optf "sys$library:libnetcdf/lib"
	write optf "sys$library:libfox/lib"
	close optf
	link/exec=abinit.exe [.src.98_main]abinit.obj,[]abinit.opt/opt

[.src.98_main]abinis.obj : [.src.98_main]abinit.F90 [.src.10_defs]libdefs.olb
	set default [.src.98_main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=abinit.f_ abinit.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=abinit.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.10_defs]) abinit.f_
	delete abinit.f_;*
	set default [--]

DEP_AIM= $(AR_BADER) $(AR_GEOMETRY) $(AR_IOWFDENPOT) $(AR_UTIL) $(AR_BASIS)\
	$(AR_DEFS) $(AR_MGMPI) $(AR_NUMERIC) $(AR_PARSER) $(AR_NUMERICF90)\
	$(AR_UTIL) $(AR_IOMPI) $(AR_UTIL)

DEP_AIM_LIB=$(AR_BADER)/lib,$(AR_GEOMETRY)/lib,$(AR_IOWFDENPOT)/lib,\
	$(AR_PARSER)/lib,$(AR_UTIL)/lib,$(AR_MGMPI)/lib,$(AR_BASIS)/lib,\
	$(AR_DEFS)/lib,$(AR_NUMERIC)/lib,$(AR_NUMERICF90)/lib,$(AR_UTIL)/lib,\
	$(AR_IOMPI)/lib,$(AR_UTIL)/lib

aim : aim.exe
	write sys$output " aim.exe has been made "
	write sys$output " "

aim.exe : [.src.98_main]aim.obj $(DEP_AIM)
	write sys$output " aim.exe will be made "
	write sys$output " "
	link/exec=aim.exe [.src.98_main]aim.obj,$(DEP_AIM_LIB)

[.src.98_main]aim.obj : [.src.98_main]aim.F90 [.src.10_defs]libdefs.olb
	set default [.src.98_main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=aim.f_ aim.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=aim.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.10_defs]) aim.f_
	delete aim.f_;*
	set default [--]

DEP_ANADDB=$(AR_DDB) $(AR_RESPONSE) $(AR_WFS) $(AR_COMMON) $(AR_IOWFDENPOT)\
 $(AR_RECIPSPACE) $(AR_GEOMETRY) $(AR_NONLOCAL) $(AR_PARSER) $(AR_SPACEPAR) \
 $(AR_UTIL) $(AR_MGMPI) $(AR_BASIS) $(AR_DEFS) $(AR_NUMERIC) $(AR_NUMERICF90) \
 $(AR_UTIL) $(AR_OCCEIG) $(AR_HIDEMPI) $(AR_IOMPI) $(AR_WVL_WFS)

DEP_ANADDB_LIB=$(AR_DDB)/lib,$(AR_RESPONSE)/lib,$(AR_WFS)/lib,\
	$(AR_WVL_WFS)/lib,\
	$(AR_COMMON)/lib,$(AR_IOWFDENPOT)/lib,$(AR_OCCEIG)/lib,\
	$(AR_RECIPSPACE)/lib,$(AR_IOMPI)/lib,$(AR_HIDEMPI)/lib,\
	$(AR_GEOMETRY)/lib,$(AR_NONLOCAL)/lib,$(AR_PARSER)/lib,\
	$(AR_SPACEPAR)/lib,$(AR_UTIL)/lib,$(AR_MGMPI)/lib,$(AR_BASIS)/lib,\
	$(AR_DEFS)/lib,$(AR_NUMERIC)/lib,$(AR_NUMERICF90)/lib,$(AR_UTIL)/lib

anaddb : anaddb.exe
	write sys$output " anaddb.exe has been made "
	write sys$output " "

anaddb.exe : [.src.98_main]anaddb.obj $(DEP_ANADDB)
	write sys$output " anaddb.exe will be made "
	write sys$output " "
	link/exec=anaddb.exe [.src.98_main]anaddb.obj,$(DEP_ANADDB_LIB),\
	sys$library:libnetcdf/lib

[.src.98_main]anaddb.obj : [.src.98_main]anaddb.F90 [.src.10_defs]libdefs.olb
	set default [.src.98_main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=anaddb.f_ anaddb.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=anaddb.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.10_defs]) anaddb.f_
	delete anaddb.f_;*
	set default [--]

DEP_CUT3D=$(AR_CUT3D) $(AR_SEQ) $(AR_COMMON) $(AR_WFS) $(AR_IOWFDENPOT) \
	$(AR_GEOMETRY) $(AR_FFTS) $(AR_SPACEPAR) $(AR_RECIPSPACE) $(AR_UTIL) \
	$(AR_CONTRACT) $(AR_MGMPI) $(AR_BASIS) $(AR_DEFS) $(AR_FFTNEW)\
	$(AR_NUMERIC) $(AR_NONLOCAL) $(AR_NUMERICF90) $(AR_UTIL)\
	$(AR_OCCEIG) $(AR_PARSER) $(AR_HIDEMPI) $(AR_IOMPI) $(AR_WVL_WFS)

DEP_CUT3D_LIB=$(AR_CUT3D)/lib,$(AR_SEQ)/lib,$(AR_COMMON)/lib,$(AR_WFS)/lib,\
	$(AR_WVL_WFS)/lib,\
	$(AR_IOWFDENPOT)/lib,$(AR_OCCEIG)/lib,$(AR_GEOMETRY)/lib,\
	$(AR_FFTS)/lib,\
	$(AR_SPACEPAR)/lib,$(AR_RECIPSPACE)/lib,$(AR_UTIL)/lib,\
	$(AR_PARSER)/lib,$(AR_IOMPI)/lib,$(AR_HIDEMPI)/lib,\
	$(AR_NONLOCAL)/lib,$(AR_CONTRACT)/lib,$(AR_MGMPI)/lib,$(AR_BASIS)/lib,\
	$(AR_DEFS)/lib,$(AR_FFTNEW)/lib,$(AR_NUMERIC)/lib,\
	$(AR_NUMERICF90)/lib,$(AR_UTIL)/lib,sys$library:libnetcdf/lib

cut3d : cut3d.exe
	write sys$output " cut3d.exe has been made "
	write sys$output " "

cut3d.exe : [.src.98_main]cut3d.obj $(DEP_CUT3D)
	write sys$output " cut3d.exe will be made "
	write sys$output " "
	link/exec=cut3d.exe [.src.98_main]cut3d.obj,$(DEP_CUT3D_LIB)

[.src.98_main]cut3d.obj : [.src.98_main]cut3d.F90 [.src.10_defs]libdefs.olb
	set default [.src.98_main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=cut3d.f_ cut3d.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=cut3d.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.10_defs]) cut3d.f_
	delete cut3d.f_;*
	set default [--]

DEP_MRGDDB=$(AR_DDB) $(AR_RESPONSE) $(AR_UTIL) $(AR_MGMPI) $(AR_BASIS)\
	$(AR_DEFS) $(AR_NUMERIC)

DEP_MRGDDB_LIB=$(AR_DDB)/lib,$(AR_RESPONSE)/lib,$(AR_UTIL)/lib,\
	$(AR_MGMPI)/lib,$(AR_BASIS)/lib,$(AR_DEFS)/lib,$(AR_NUMERIC)/lib

mrgddb : mrgddb.exe
	write sys$output " mrgddb.exe has been made "
	write sys$output " "

mrgddb.exe : [.src.98_main]mrgddb.obj $(DEP_MRGDDB)
	write sys$output " mrgddb.exe will be made "
	write sys$output " "
	link/exec=mrgddb.exe [.src.98_main]mrgddb.obj,$(DEP_MRGDDB_LIB)

[.src.98_main]mrgddb.obj : [.src.98_main]mrgddb.F90 [.src.10_defs]libdefs.olb
	set default [.src.98_main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=mrgddb.f_ mrgddb.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=mrgddb.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.10_defs]) mrgddb.f_
	delete mrgddb.f_;*
	set default [--]

DEP_LWF=$(AR_LWF) $(AR_PARSER) $(AR_UTIL) $(AR_MGMPI) $(AR_BASIS) $(AR_DEFS)\
	$(AR_NUMERIC)

DEP_LWF_LIB=$(AR_LWF)/lib,$(AR_PARSER)/lib,$(AR_UTIL)/lib,$(AR_MGMPI)/lib,\
	$(AR_BASIS)/lib,$(AR_DEFS)/lib,$(AR_NUMERIC)/lib

lwf : lwf.exe
	write sys$output " lwf.exe has been made "
	write sys$output " "

lwf.exe : [.src.98_main]lwf.obj $(DEP_LWF)
	write sys$output " lwf.exe will be made "
	write sys$output " "
	link/exec=lwf.exe [.src.98_main]lwf.obj,$(DEP_LWF_LIB)

[.src.98_main]lwf.obj : [.src.98_main]lwf.F90 [.src.10_defs]libdefs.olb
	set default [.src.98_main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=lwf.f_ lwf.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=lwf.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.10_defs]) lwf.f_
	delete lwf.f_;*
	set default [--]

DEP_CONDUCTI=$(AR_COMMON) $(AR_IOWFDENPOT) $(AR_GEOMETRY) $(AR_UTIL)\
	$(AR_MGMPI) $(AR_BASIS) $(AR_DEFS) $(AR_NUMERIC) $(AR_NUMERICF90)\
	$(AR_UTIL) $(AR_OCCEIG) $(AR_IOMPI) $(AR_WVL_WFS)

DEP_CONDUCTI_LIB=$(AR_COMMON)/lib,$(AR_IOWFDENPOT)/lib,$(AR_OCCEIG)/lib,\
	$(AR_GEOMETRY)/lib,$(AR_IOMPI)/lib,$(AR_WVL_WFS)/lib,\
	$(AR_UTIL)/lib,$(AR_MGMPI)/lib,$(AR_BASIS)/lib,$(AR_DEFS)/lib,\
	$(AR_NUMERIC)/lib,$(AR_NUMERICF90)/lib,$(AR_UTIL)/lib,\
	sys$library:libnetcdf/lib

conducti : conducti.exe
	write sys$output " conducti.exe has been made "
	write sys$output " "

conducti.exe : [.src.98_main]conducti.obj $(DEP_CONDUCTI)
	write sys$output " conducti.exe will be made "
	write sys$output " "
	link/exec=conducti.exe [.src.98_main]conducti.obj,$(DEP_CONDUCTI_LIB)

[.src.98_main]conducti.obj : [.src.98_main]conducti.F90 [.src.10_defs]libdefs.olb
	set default [.src.98_main]
	cc/define=(VMS,HAVE_NETCDF,HAVE_XMLF90)/comment=as_is/prep=conducti.f_ conducti.F90
	f90/source=free/convert=big_endian/name=lowercase\
	/opti=(tune=host,nopipe)/arch=host/debug/object=conducti.obj\
	/check=(noarg,bo,form,fp,noout,over,powe,under)\
	/include=([--],[--.src.10_defs]) conducti.f_
	delete conducti.f_;*
	set default [--]

[.src.10_defs]libdefs.olb :
	$(MMS) [.src.10_defs]descrip.mms
	set default [.src.10_defs]
	$(MMS) m_pseudo_types.obj
	$(MMS) defs_basis.obj
	$(MMS) defs_datatypes.obj
	$(MMS) interfaces_54_spacepar.obj
	$(MMS) interfaces_56_recipspace.obj
	$(MMS) interfaces_32_util.obj
	$(MMS) interfaces_55_abiutil.obj
	$(MMS) interfaces_62_cg_noabirule.obj
	$(MMS) interfaces_56_xc.obj
	$(MMS) interfaces_53_ffts.obj
	$(MMS)
	set default [--]

[.src.95_drive]lib95_drive.olb :
	$(MMS) [.src.95_drive]descrip.mms
	set default [.src.95_drive]
	$(MMS)
	set default [--]

[.src.62_occeig]lib62_occeig.olb :
	$(MMS) [.src.62_occeig]descrip.mms
	set default [.src.62_occeig]
	$(MMS)
	set default [--]

[.src.42_nlstrain]lib42_nlstrain.olb :
	$(MMS) [.src.42_nlstrain]descrip.mms
	set default [.src.42_nlstrain]
	$(MMS)
	set default [--]

[.src.62_poisson]lib62_poisson.olb :
	$(MMS) [.src.62_poisson]descrip.mms
	set default [.src.62_poisson]
	$(MMS)
	set default [--]

[.src.83_cut3d]lib83_cut3d.olb :
	$(MMS) [.src.83_cut3d]descrip.mms
	set default [.src.83_cut3d]
	$(MMS)
	set default [--]

[.src.79_seqpar_mpi]lib18abinis.olb :
	$(MMS) [.src.79_seqpar_mpi]descrip.mms
	set default [.src.79_seqpar_mpi]
	$(MMS)
	copy lib79_seqpar_mpi.olb lib18abinis.olb
	set default [--]

[.src.77_suscep]lib77_suscep.olb :
	$(MMS) [.src.77_suscep]descrip.mms
	set default [.src.77_suscep]
	$(MMS)
	set default [--]

[.src.77_ddb]lib77_ddb.olb :
	$(MMS) [.src.77_ddb]descrip.mms
	set default [.src.77_ddb]
	$(MMS)
	set default [--]

[.src.77_lwf]lib77_lwf.olb :
	$(MMS) [.src.77_lwf]descrip.mms
	set default [.src.77_lwf]
	$(MMS)
	set default [--]

[.src.72_response]lib72_response.olb :
	$(MMS) [.src.72_response]descrip.mms
	set default [.src.72_response]
	$(MMS)
	set default [--]

[.src.67_common]lib67_common.olb :
	$(MMS) [.src.67_common]descrip.mms
	set default [.src.67_common]
	$(MMS)
	set default [--]

[.src.62_iowfdenpot]lib62_iowfdenpot.olb :
	$(MMS) [.src.62_iowfdenpot]descrip.mms
	set default [.src.62_iowfdenpot]
	$(MMS)
	set default [--]

[.src.14wfs]lib14wfs.olb :
	$(MMS) [.src.14wfs]descrip.mms
	set default [.src.14wfs]
	$(MMS)
	set default [--]

[.src.62_wvl_wfs]lib62_wvl_wfs.olb :
	$(MMS) [.src.62_wvl_wfs]descrip.mms
	set default [.src.62_wvl_wfs]
	$(MMS)
	set default [--]

[.src.68_gw]lib68_gw.olb :
	$(MMS) [.src.68_gw]descrip.mms
	set default [.src.68_gw]
	$(MMS)
	set default [--]

[.src.67_recursion]lib67_recursion.olb :
	$(MMS) [.src.67_recursion]descrip.mms
	set default [.src.67_recursion]
	$(MMS)
	set default [--]

[.src.68_rsprc]lib68_rsprc.olb :
	$(MMS) [.src.68_rsprc]descrip.mms
	set default [.src.68_rsprc]
	$(MMS)
	set default [--]

[.src.62_cg_noabirule]lib62_cg_noabirule.olb :
	$(MMS) [.src.62_cg_noabirule]descrip.mms
	set default [.src.62_cg_noabirule]
	$(MMS)
	set default [--]

[.src.12_hide_mpi]lib12_hide_mpi.olb :
	$(MMS) [.src.12_hide_mpi]descrip.mms
	set default [.src.12_hide_mpi]
	$(MMS)
	set default [--]

[.src.56_io_mpi]lib56_io_mpi.olb :
        $(MMS) [.src.56_io_mpi]descrip.mms
        set default [.src.56_io_mpi]
        $(MMS)
        set default [--]

[.src.57_iovars]lib57_iovars.olb :
	$(MMS) [.src.57_iovars]descrip.mms
	set default [.src.57_iovars]
	$(MMS)
	set default [--]

[.src.59_ionetcdf]lib59_ionetcdf.olb :
	$(MMS) [.src.59_ionetcdf]descrip.mms
	set default [.src.59_ionetcdf]
	$(MMS)
	set default [--]

[.src.13paw]lib13paw.olb :
	$(MMS) [.src.13paw]descrip.mms
	set default [.src.13paw]
	$(MMS)
	set default [--]

[.src.56_recipspace]lib56_recipspace.olb :
	$(MMS) [.src.56_recipspace]descrip.mms
	set default [.src.56_recipspace]
	$(MMS)
	set default [--]

[.src.56_xc]lib56_xc.olb :
	$(MMS) [.src.56_xc]descrip.mms
	set default [.src.56_xc]
	$(MMS)
	set default [--]

[.src.47_xml]lib47_xml.olb :
	$(MMS) [.src.47_xml]descrip.mms
	set default [.src.47_xml]
	$(MMS)
	set default [--]

[.src.13nonlocal]lib13nonlocal.olb :
	$(MMS) [.src.13nonlocal]descrip.mms
	set default [.src.13nonlocal]
	$(MMS)
	set default [--]

[.src.53_ffts]lib53_ffts.olb :
	$(MMS) [.src.53_ffts]descrip.mms
	set default [.src.53_ffts]
	$(MMS)
	set default [--]

[.src.13psp]lib13psp.olb :
	$(MMS) [.src.13psp]descrip.mms
	set default [.src.13psp]
	$(MMS) smoothvlocal.obj
	$(MMS)
	set default [--]

[.src.42_geometry]lib42_geometry.olb :
	$(MMS) [.src.42_geometry]descrip.mms
	set default [.src.42_geometry]
	$(MMS)
	set default [--]

[.src.45_geomoptim]lib45_geomoptim.olb :
	$(MMS) [.src.45_geomoptim]descrip.mms
	set default [.src.45_geomoptim]
	$(MMS)
	set default [--]

[.src.42_parser]lib42_parser.olb :
	$(MMS) [.src.42_parser]descrip.mms
	set default [.src.42_parser]
	$(MMS)
	set default [--]

[.src.54_spacepar]lib54_spacepar.olb :
	$(MMS) [.src.54_spacepar]descrip.mms
	set default [.src.54_spacepar]
	$(MMS)
	set default [--]

[.src.62_bader]lib62_bader.olb :
	$(MMS) [.src.62_bader]descrip.mms
	set default [.src.62_bader]
	$(MMS)
	set default [--]

[.src.55_abiutil]lib55_abiutil.olb :
        $(MMS) [.src.55_abiutil]descrip.mms
        set default [.src.55_abiutil]
        $(MMS)
        set default [--]

[.src.32_contract]lib32_contract.olb :
	$(MMS) [.src.32_contract]descrip.mms
	set default [.src.32_contract]
	$(MMS)
	set default [--]

[.src.51_manage_mpi]lib51_manage_mpi.olb :
	$(MMS) [.src.51_manage_mpi]descrip.mms
	set default [.src.51_manage_mpi]
	$(MMS)
	set default [--]

[.src.32_util]lib32_util.olb :
	$(MMS) [.src.32_util]descrip.mms
	set default [.src.32_util]
	$(MMS)
	set default [--]

[.src.14_hidewrite]lib14_hidewrite.olb :
	$(MMS) [.src.14_hidewrite]descrip.mms
	set default [.src.14_hidewrite]
	$(MMS)
	set default [--]

[.src.16_hideleave]lib16_hideleave.olb :
        $(MMS) [.src.16_hideleave]descrip.mms
        set default [.src.16_hideleave]
        $(MMS)
        set default [--]

[.src.18_timing]lib18_timing.olb :
        $(MMS) [.src.18_timing]descrip.mms
        set default [.src.18_timing]
        $(MMS)
        set default [--]

[.src.52_fft_mpi_noabirule]descrip.mms : crea_descrip_mms.exe
	set def [.src.52_fft_mpi_noabirule]
	run [--]crea_descrip_mms
	set def [--]

[.lib.numeric]descrip.mms : crea_descrip_mms.exe
	set def [.lib.numeric]
	run [--]crea_descrip_mms
	set def [--]

[.lib.numericf90]descrip.mms : crea_descrip_mms.exe
	set def [.lib.numericf90]
	run [--]crea_descrip_mms
	set def [--]

[.src.14_hidewrite]descrip.mms : crea_descrip_mms.exe
	set def [.src.14_hidewrite]
	run [--]crea_descrip_mms
	set def [--]

[.src.16_hideleave]descrip.mms : crea_descrip_mms.exe
        set def [.src.16_hideleave]
        run [--]crea_descrip_mms
        set def [--]

[.src.18_timing]descrip.mms : crea_descrip_mms.exe
        set def [.src.18_timing]
        run [--]crea_descrip_mms
        set def [--]

[.src.51_manage_mpi]descrip.mms : crea_descrip_mms.exe
	set def [.src.51_manage_mpi]
	run [--]crea_descrip_mms
	set def [--]

[.src.32_contract]descrip.mms : crea_descrip_mms.exe
	set def [.src.32_contract]
	run [--]crea_descrip_mms
	set def [--]

[.src.32_util]descrip.mms : crea_descrip_mms.exe
	set def [.src.32_util]
	run [--]crea_descrip_mms
	set def [--]

[.src.55_abiutil]descrip.mms : crea_descrip_mms.exe
        set def [.src.55_abiutil]
        run [--]crea_descrip_mms
        set def [--]

[.src.62_bader]descrip.mms : crea_descrip_mms.exe
	set def [.src.62_bader]
	run [--]crea_descrip_mms
	set def [--]

[.src.53_ffts]descrip.mms : crea_descrip_mms.exe
	set def [.src.53_ffts]
	run [--]crea_descrip_mms
	set def [--]

[.src.42_geometry]descrip.mms : crea_descrip_mms.exe
	set def [.src.42_geometry]
	run [--]crea_descrip_mms
	set def [--]

[.src.13nonlocal]descrip.mms : crea_descrip_mms.exe
	set def [.src.13nonlocal]
	run [--]crea_descrip_mms
	set def [--]

[.src.42_parser]descrip.mms : crea_descrip_mms.exe
	set def [.src.42_parser]
	run [--]crea_descrip_mms
	set def [--]

[.src.13psp]descrip.mms : crea_descrip_mms.exe
	set def [.src.13psp]
	run [--]crea_descrip_mms
	set def [--]

[.src.54_spacepar]descrip.mms : crea_descrip_mms.exe
	set def [.src.54_spacepar]
	run [--]crea_descrip_mms
	set def [--]

[.src.68_gw]descrip.mms : crea_descrip_mms.exe
	set def [.src.68_gw]
	run [--]crea_descrip_mms
	set def [--]

[.src.67_recursion]descrip.mms : crea_descrip_mms.exe
	set def [.src.67_recursion]
	run [--]crea_descrip_mms
	set def [--]

[.src.68_rsprc]descrip.mms : crea_descrip_mms.exe
	set def [.src.68_rsprc]
	run [--]crea_descrip_mms
	set def [--]

[.src.62_cg_noabirule]descrip.mms : crea_descrip_mms.exe
	set def [.src.62_cg_noabirule]
	run [--]crea_descrip_mms
	set def [--]

[.src.12_hide_mpi]descrip.mms : crea_descrip_mms.exe
	set def [.src.12_hide_mpi]
	run [--]crea_descrip_mms
	set def [--]

[.src.56_io_mpi]descrip.mms : crea_descrip_mms.exe
        set def [.src.56_io_mpi]
        run [--]crea_descrip_mms
        set def [--]

[.src.57_iovars]descrip.mms : crea_descrip_mms.exe
	set def [.src.57_iovars]
	run [--]crea_descrip_mms
	set def [--]

[.src.59_ionetcdf]descrip.mms : crea_descrip_mms.exe
	set def [.src.59_ionetcdf]
	run [--]crea_descrip_mms
	set def [--]

[.src.13paw]descrip.mms : crea_descrip_mms.exe
	set def [.src.13paw]
	run [--]crea_descrip_mms
	set def [--]

[.src.56_recipspace]descrip.mms : crea_descrip_mms.exe
	set def [.src.56_recipspace]
	run [--]crea_descrip_mms
	set def [--]

[.src.56_xc]descrip.mms : crea_descrip_mms.exe
	set def [.src.56_xc]
	run [--]crea_descrip_mms
	set def [--]

[.src.47_xml]descrip.mms : crea_descrip_mms.exe
	set def [.src.47_xml]
	run [--]crea_descrip_mms
	set def [--]

[.src.62_iowfdenpot]descrip.mms : crea_descrip_mms.exe
	set def [.src.62_iowfdenpot]
	run [--]crea_descrip_mms
	set def [--]

[.src.62_wvl_wfs]descrip.mms : crea_descrip_mms.exe
	set def [.src.62_wvl_wfs]
	run [--]crea_descrip_mms
	set def [--]

[.src.67_common]descrip.mms : crea_descrip_mms.exe
	set def [.src.67_common]
	run [--]crea_descrip_mms
	set def [--]

[.src.72_response]descrip.mms : crea_descrip_mms.exe
	set def [.src.72_response]
	run [--]crea_descrip_mms
	set def [--]

[.src.77_ddb]descrip.mms : crea_descrip_mms.exe
	set def [.src.77_ddb]
	run [--]crea_descrip_mms
	set def [--]

[.src.77_lwf]descrip.mms : crea_descrip_mms.exe
	set def [.src.77_lwf]
	run [--]crea_descrip_mms
	set def [--]

[.src.77_suscep]descrip.mms : crea_descrip_mms.exe
	set def [.src.77_suscep]
	run [--]crea_descrip_mms
	set def [--]

[.src.79_seqpar_mpi]descrip.mms : crea_descrip_mms.exe
	set def [.src.79_seqpar_mpi]
	run [--]crea_descrip_mms
	set def [--]

[.src.83_cut3d]descrip.mms : crea_descrip_mms.exe
	set def [.src.83_cut3d]
	run [--]crea_descrip_mms
	set def [--]

[.src.95_drive]descrip.mms : crea_descrip_mms.exe
	set def [.src.95_drive]
	run [--]crea_descrip_mms
	set def [--]

[.src.42_nlstrain]descrip.mms : crea_descrip_mms.exe
	set def [.src.42_nlstrain]
	run [--]crea_descrip_mms
	set def [--]

[.src.62_poisson]descrip.mms : crea_descrip_mms.exe
	set def [.src.62_poisson]
	run [--]crea_descrip_mms
	set def [--]

[.src.45_geomoptim]descrip.mms : crea_descrip_mms.exe
	set def [.src.45_geomoptim]
	run [--]crea_descrip_mms
	set def [--]

[.src.62_occeig]descrip.mms : crea_descrip_mms.exe
	set def [.src.62_occeig]
	run [--]crea_descrip_mms
	set def [--]

[.src.10_defs]descrip.mms : crea_descrip_mms.exe
	set def [.src.10_defs]
	run [--]crea_descrip_mms
	set def [--]

prepare_input :
	set def [.tests.tutorial.Input]
	mcr [---]vms_prepare_input t1x.files
	mcr [---]vms_prepare_input t2x.files
	mcr [---]vms_prepare_input t3x.files
	mcr [---]vms_prepare_input t4x.files
	set def [---]

# Build-in tests
tests : test1 test2 test3 test4 test5 test6
	write sys$output "Performed test 1-6"

test12 : test1 test2 
	write sys$output "Performed test 1-2"

test123 : test1 test2 test3 
	write sys$output "Performed test 1-3"

test124 : test1 test2 test4 
	write sys$output "Performed test 1,2,4"

test125 : test1 test2 test5 
	write sys$output "Performed test 1,2,5"

test1 :
	@ @[.tests.Scripts]run-basic-tests.com [.tests.built-in] 1

test2 :
	@ @[.tests.Scripts]run-basic-tests.com [.tests.built-in] 2

test3 :
	@ @[.tests.Scripts]run-basic-tests.com [.tests.built-in] 3

test4 :
	@ @[.tests.Scripts]run-basic-tests.com [.tests.built-in] 4

test5 :
	@ @[.tests.Scripts]run-basic-tests.com [.tests.built-in] 5

test6 :
	@ @[.tests.Scripts]run-basic-tests.com [.tests.built-in] 6
