!!****m* ABINIT/m_rdmft
!! NAME
!!  m_rdmft
!!
!! FUNCTION
!! Calculate the RDMFT energy.
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2020 ABINIT group (MG, GMR, VO, LR, RWG, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rdmft

 use defs_basis
 use m_gwdefs
 use defs_wvltypes
 use m_dtset
 use m_xmpi
 use m_xomp
 use m_errors
 use m_abicore
 use m_ab7_mixing
 use m_nctk
 use m_kxc
 use m_distribfft
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr
 use m_wfd
 use m_dtfil
 use m_crystal

 use defs_datatypes,  only : pseudopotential_type, ebands_t
 use defs_abitypes,   only : MPI_type
 use m_time,          only : timab
 use m_numeric_tools, only : imax_loc
 use m_fstrings,      only : strcat, sjoin, itoa, basename, ktoa, ltoa
 use m_hide_blas,     only : xdotc
 use m_io_tools,      only : open_file, file_exists, iomode_from_fname
 use m_mpinfo,        only : destroy_mpi_enreg, initmpi_seq
 use m_geometry,      only : normv, mkrdim, metric
 use m_fftcore,       only : print_ngfft
 use m_fft_mesh,      only : get_gftt, setmesh
 use m_fft,           only : fourdp
 use m_ioarr,         only : fftdatar_write, read_rhor
 use m_ebands,        only : ebands_update_occ, ebands_copy, ebands_report_gap, get_valence_idx, get_bandenergy, &
                             ebands_free, ebands_init, ebands_ncwrite, ebands_interpolate_kpath, get_eneocc_vect, &
                             enclose_degbands, get_gaps, gaps_t
 use m_energies,      only : energies_type, energies_init
 use m_bz_mesh,       only : kmesh_t, kmesh_free, littlegroup_t, littlegroup_init, littlegroup_free, &
                             kmesh_init, has_BZ_item, isamek, get_ng0sh, kmesh_print, &
                             get_bz_item, has_IBZ_item, find_qmesh
 use m_gsphere,       only : gsphere_t, gsph_init, gsph_free, merge_and_sort_kg, gsph_extend, setshells
 use m_kg,            only : getph
 use m_xcdata,        only : get_xclevel
 use m_vcoul,         only : vcoul_t, vcoul_init, vcoul_free
 use m_qparticles,    only : wrqps, rdqps, rdgw, show_QP, updt_m_ks_to_qp
 use m_sigma,         only : sigma_init, sigma_free, sigma_ncwrite, sigma_t, sigma_get_exene, &
                             write_sigma_header, write_sigma_results
 use m_esymm,         only : esymm_t, esymm_free, esymm_failed
 use m_melemts,       only : melflags_reset, melements_print, melements_free, melflags_t, melements_t, melements_zero
 use m_pawang,        only : pawang_type
 use m_pawrad,        only : pawrad_type
 use m_pawtab,        only : pawtab_type, pawtab_print, pawtab_get_lsize
 use m_paw_an,        only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
 use m_paw_ij,        only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify, paw_ij_print
 use m_pawfgrtab,     only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free, pawfgrtab_print
 use m_pawrhoij,      only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, &
                             pawrhoij_inquire_dim, pawrhoij_symrhoij, pawrhoij_unpack
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free, paw_overlap
 use m_pawdij,        only : pawdij, symdij_all
 use m_pawfgr,        only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_paw_pwaves_lmn,only : paw_pwaves_lmn_t, paw_pwaves_lmn_init, paw_pwaves_lmn_free
 use m_pawpwij,       only : pawpwff_t, pawpwff_init, pawpwff_free
 use m_paw_slater,    only : paw_mkdijexc_core, paw_dijhf
 use m_paw_dmft,      only : paw_dmft_type
 use m_paw_sphharm,   only : setsym_ylm
 use m_paw_mkrho,     only : denfgr
 use m_paw_nhat,      only : nhatgrid,pawmknhat
 use m_paw_tools,     only : chkpawovlp,pawprt
 use m_paw_denpot,    only : pawdenpot
 use m_paw_init,      only : pawinit,paw_gencond
 use m_classify_bands,only : classify_bands
 use m_wfk,           only : wfk_read_eigenvalues
 use m_io_kss,        only : make_gvec_kss
 use m_vhxc_me,       only : calc_vhxc_me
 use m_sigx,          only : calc_sigx_me
 use m_setvtr,        only : setvtr
 use m_mkrho,         only : prtrhomxmn
 use m_pspini,        only : pspini
 !use m_paw_correlations,only : pawpuxinit
 use m_spacepar,      only : hartre

 implicit none

 private
!!***

 public :: rdmft
!!***

contains
!!***

!!****f* m_rmdft/rdmft
!! NAME
!! rdfmt
!!
!! FUNCTION
!! Calculate the matrix elements of the self-energy operator.
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! codvsn=code version
!! Dtfil<type(datafiles_type)>=variables related to files
!! Dtset<type(dataset_type)>=all input variables for this dataset
!! Pawang<type(pawang_type)>=paw angular mesh and related data
!! Pawrad(ntypat*usepaw)<type(pawrad_type)>=paw radial mesh and related data
!! Pawtab(ntypat*usepaw)<type(pawtab_type)>=paw tabulated starting data
!! Psps<type(pseudopotential_type)>=variables related to pseudopotentials
!!   Before entering the first time in sigma, a significant part of Psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,ntypat,n1xccc,usepaw,useylm,
!!   and the arrays dimensioned to npsp. All the remaining components of Psps are to be initialized in
!!   the call to pspini. The next time the code enters screening, Psps might be identical to the
!!   one of the previous Dtset, in which case, no reinitialisation is scheduled in pspini.F90.
!! rprim(3,3)=dimensionless real space primitive translations
!!
!! OUTPUT
!!  Output is written on the main abinit output file. Some results are stored in external files
!!
!! PARENTS
!!      driver
!!
!! NOTES
!!
!! ON THE USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut) for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ... are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg) for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ... Total density, potentials, ... are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used. It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf) are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! CHILDREN
!!      calc_sigc_me,calc_sigx_me,calc_ucrpa,calc_vhxc_me,chkpawovlp
!!      classify_bands,cohsex_me,denfgr,destroy_mpi_enreg
!!      ebands_copy,ebands_free,ebands_interpolate_kpath,ebands_report_gap
!!      ebands_update_occ,em1results_free,energies_init,esymm_free
!!      fftdatar_write,fourdp,get_gftt,getph,gsph_free,hdr_free
!!      init_distribfft_seq,initmpi_seq,kmesh_free,kxc_ada,kxc_driver
!!      littlegroup_free,littlegroup_init,melements_free,melements_print
!!      melements_zero,melflags_reset,metric,mkdump_er,mkrdim,nhatgrid
!!      paw_an_free,paw_an_init,paw_an_nullify,paw_check_symcprj,paw_dijhf
!!      paw_gencond,paw_ij_free,paw_ij_init,paw_ij_nullify,paw_ij_print
!!      paw_mkdijexc_core,paw_pwaves_lmn_free,paw_pwaves_lmn_init,paw_qpscgw
!!      pawcprj_alloc,pawcprj_free,pawdenpot,pawdij,pawfgr_destroy,pawfgr_init
!!      pawfgrtab_free,pawfgrtab_init,pawfgrtab_print,pawinit,pawmknhat,pawprt
!!      pawpuxinit,pawpwff_free,pawpwff_init,pawrhoij_alloc,pawrhoij_copy
!!      pawrhoij_free,pawtab_get_lsize,pawtab_print,ppm_free,ppm_init
!!      prep_calc_ucrpa,print_ngfft,prtrhomxmn,pspini,rdgw,rdqps,read_rhor
!!      setsym_ylm,setup_ppmodel,setup_sigma,setvtr,show_qp,sigma_bksmask
!!      sigma_free,sigma_init,sigma_tables,sigparams_free,solve_dyson,symdij
!!      symdij_all,test_charge,timab,updt_m_ks_to_qp,vcoul_free
!!      wfd_change_ngfft,wfd_copy,wfd_distribute_bands,wfd_free,wfd_get_cprj
!!      wfd_init,wfd_mkrho,wfd_print,wfd_read_wfk,wfd_reset_ur_cprj,wfd_rotate
!!      wfd_test_ortho,write_sigma_header,write_sigma_results,wrqps,wrtout
!!      xmpi_barrier,xmpi_bcast,xmpi_sum
!!
!! SOURCE

subroutine rdmft(acell,codvsn,Dtfil,Dtset,Pawang,Pawrad,Pawtab,Psps,rprim,converged)
!Arguments ------------------------------------
!scalars
 logical,intent(out) :: converged
 character(len=8),intent(in) :: codvsn
 type(Datafiles_type),intent(in) :: Dtfil
 type(Dataset_type),intent(inout) :: Dtset
 type(Pawang_type),intent(inout) :: Pawang
 type(Pseudopotential_type),intent(inout) :: Psps
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3)
 type(Pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Psps%usepaw)
 type(Pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Psps%usepaw)
!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp5=5,master=0,cplex1=1
 integer :: my_rank,comm,nprocs,ierr,usexcnhat,mgfftf,nfftf,nfftf_tot,psp_gencond
 real(dp) :: ucvol,ecore,ecut_eff,ecutdg_eff,gsqcutf_eff,gsqcutc_eff
 character(len=500) :: msg
 character(len=fnlen) :: wfk_fname
 type(Energies_type) :: DM_energies
 type(MPI_type) :: MPI_enreg_seq
 type(pawfgr_type) :: Pawfgr

!arrays
 integer :: gwc_ngfft(18),ngfftc(18),ngfftf(18),gwx_ngfft(18)
 real(dp),parameter ::  k0(3)=zero
 real(dp) :: gmet(3,3),gprimd(3,3),rprimd(3,3),rmet(3,3)

!************************************************************************
 write(msg,'(7a)')&
 ' RDMFT: Calculation of the GS energy ',ch10,ch10,&
 ' Incorporated in ABINIT by M. Rodriguez-Mayorga.'
 call wrtout([std_out, ab_out], msg)

#if defined HAVE_GW_DPC
 if (gwpc/=8) then
   write(msg,'(6a)')ch10,&
    'Number of bytes for double precision complex /=8 ',ch10,&
    'Cannot continue due to kind mismatch in BLAS library ',ch10,&
    'Some BLAS interfaces are not generated by abilint '
   MSG_ERROR(msg)
 end if
 write(msg,'(a,i2,a)')'.Using double precision arithmetic ; dmpc = ',gwpc,ch10
#else
 write(msg,'(a,i2,a)')'.Using single precision arithmetic ; dmpc = ',gwpc,ch10
#endif
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 !=== Initialize MPI variables, and parallelization level ===
 comm = xmpi_world; my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 converged = .FALSE.

 if (my_rank == master) then
   wfk_fname = dtfil%fnamewffk
   if (nctk_try_fort_or_ncfile(wfk_fname, msg) /= 0) then
     MSG_ERROR(msg)
   end if
 end if
 call xmpi_bcast(wfk_fname, master, comm, ierr)

! === Some variables need to be initialized/nullify at start ===
 call energies_init(DM_energies)
 usexcnhat=0
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,ab_out,rmet,rprimd,ucvol)
 !
 ! === Define FFT grid(s) sizes ===
 ! * Be careful! This mesh is only used for densities, potentials and the matrix elements of v_Hxc. It is NOT the
 ! (usually coarser) GW FFT mesh employed for the oscillator matrix elements that is defined in setmesh.F90.
 ! See also NOTES in the comments at the beginning of this file.
 ! NOTE: This mesh is defined in invars2m using ecutwfn, in GW Dtset%ecut is forced to be equal to Dtset%ecutwfn.

 call pawfgr_init(Pawfgr,Dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
& gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=gmet,k0=k0)

 ! Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfftc(2),ngfftc(3),'all')
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',ngfftf(2),ngfftf(3),'all')

 call print_ngfft(ngfftf,header='Dense FFT mesh used for densities and potentials')
 nfftf_tot=PRODUCT(ngfftf(1:3))

 ! Open and read pseudopotential files ===
 call pspini(Dtset,Dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,Pawrad,Pawtab,Psps,rprimd,comm_mpi=comm)




end subroutine rdmft
!!***

end module m_rdmft
!!***
