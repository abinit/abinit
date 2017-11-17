!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_sigma
!! NAME
!! setup_sigma
!!
!! FUNCTION
!!  Initialize the data type containing parameters for a sigma calculation.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2017 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! wfk_fname=Name of the WFK file.
!! Dtset<type(dataset_type)>=all input variables for this dataset
!! Dtfil<type(datafiles_type)>=variables related to files
!! rprim(3,3)=dimensionless real space primitive translations
!! ngfft(18)=information on the (fine) FFT grid used for the density.
!! Psps <Pseudopotential_type)>=Info on pseudopotential, only for consistency check of the WFK file
!!
!! OUTPUT
!! Sigp<sigparams_t>=Parameters governing the self-energy calculation.
!! Kmesh <kmesh_t>=Structure describing the k-point sampling.
!! Qmesh <kmesh_t>=Structure describing the q-point sampling.
!! Cryst<crystal_t>=Info on unit cell and symmetries.
!! Gsph_Max<gsphere_t>=Info on the G-sphere
!! Gsph_c<gsphere_t>=Info on the G-sphere for W and Sigma_c
!! Gsph_x<gsphere_t>=Info on the G-sphere for and Sigma_x
!! Hdr_wfk<hdr_type>=The header of the WFK file
!! Hdr_out<hdr_type>=The header to be used for the results of sigma calculations.
!! Vcp<vcoul_t>= Datatype gathering information on the coulombian interaction and the cutoff technique.
!! Er<Epsilonm1_results>=Datatype storing data used to construct the screening (partially Initialized in OUTPUT)
!! KS_BSt<ebands_t>=The KS energies and occupation factors.
!! gwc_ngfft(18), gwx_ngfft(18)= FFT meshes for the oscillator strengths used for the correlated and the
!!   exchange part of the self-energy, respectively.
!! comm=MPI communicator.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setup_sigma(codvsn,wfk_fname,acell,rprim,ngfftf,Dtset,Dtfil,Psps,Pawtab,&
& gwx_ngfft,gwc_ngfft,Hdr_wfk,Hdr_out,Cryst,Kmesh,Qmesh,KS_BSt,Gsph_Max,Gsph_x,Gsph_c,Vcp,Er,Sigp,comm)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_nctk
 use m_hdr

 use m_gwdefs,        only : GW_Q0_DEFAULT, SIG_GW_AC, sigparams_t, sigma_is_herm, sigma_needs_w
 use m_io_tools,      only : file_exists
 use m_fstrings,      only : basename, sjoin, ktoa, ltoa
 use m_crystal,       only : crystal_print, idx_spatial_inversion, crystal_t
 use m_crystal_io,    only : crystal_from_hdr
 use m_bz_mesh,       only : kmesh_t, kmesh_init, has_BZ_item, isamek, get_ng0sh, kmesh_print,&
&                            get_bz_item, has_IBZ_item, find_qmesh
 use m_ebands,        only : ebands_init, enclose_degbands, get_valence_idx, ebands_update_occ, ebands_report_gap, &
&                            get_gaps, gaps_free, gaps_t, gaps_print
 use m_vcoul,         only : vcoul_t, vcoul_init, vcoul_free
 use m_fft_mesh,      only : setmesh
 use m_gsphere,       only : gsphere_t, gsph_init, merge_and_sort_kg, gsph_extend, setshells
 use m_screening,     only : init_er_from_file, epsilonm1_results
 use m_pawtab,        only : pawtab_type
 use m_pawrhoij,      only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free
 use m_io_kss,        only : make_gvec_kss
 use m_wfk,           only : wfk_read_eigenvalues
 use m_xcdata,        only : get_xclevel

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_sigma'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_56_io_mpi
 use interfaces_70_gw, except_this_one => setup_sigma
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 character(len=6),intent(in) :: codvsn
 character(len=*),intent(in) :: wfk_fname
 type(Datafiles_type),intent(in) :: Dtfil
 type(Dataset_type),intent(inout) :: Dtset
 type(Pseudopotential_type),intent(in) :: Psps
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Dtset%usepaw)
 type(sigparams_t),intent(out) :: Sigp
 type(Epsilonm1_results),intent(out) :: Er
 type(ebands_t),intent(out) :: KS_BSt
 type(kmesh_t),intent(out) :: Kmesh,Qmesh
 type(crystal_t),intent(out) :: Cryst
 type(gsphere_t),intent(out) :: Gsph_Max,Gsph_x,Gsph_c
 type(Hdr_type),intent(out) :: Hdr_wfk,Hdr_out
 type(vcoul_t),intent(out) :: Vcp
!arrays
 integer,intent(in) :: ngfftf(18)
 integer,intent(out) :: gwc_ngfft(18),gwx_ngfft(18)
 real(dp),intent(in) :: acell(3),rprim(3,3)

!Local variables-------------------------------
!scalars
 integer,parameter :: pertcase0=0,master=0
 integer :: bantot,enforce_sym,gwcalctyp,ib,ibtot,icutcoul_eff,ii,ikcalc,ikibz,io,isppol,itypat,jj,method
 integer :: mod10,mqmem,mband,ng_kss,nsheps,ikcalc2bz,ierr,gap_err,ng
 integer :: gwc_nfftot,gwx_nfftot,nqlwl,test_npwkss,my_rank,nprocs,ik,nk_found,ifo,timrev,usefock_ixc
 integer :: iqbz,isym,iq_ibz,itim,ic,pinv,ig1,ng_sigx,spin,gw_qprange,ivcoul_init,nvcoul_init,xclevel_ixc
 real(dp),parameter :: OMEGASIMIN=0.01d0,tol_enediff=0.001_dp*eV_Ha
 real(dp) :: domegas,domegasi,ucvol,rcut
 logical,parameter :: linear_imag_mesh=.TRUE.
 logical :: ltest,remove_inv,changed,found
 character(len=500) :: msg
 character(len=fnlen) :: fname,fcore,string
 type(wvl_internal_type) :: wvl
 type(gaps_t) :: gaps
!arrays
 integer :: ng0sh_opt(3),G0(3),q_umklp(3),kpos(6)
 integer,allocatable :: npwarr(:),val_indeces(:,:)
 integer,pointer :: gvec_kss(:,:),gsphere_sigx_p(:,:)
 integer,pointer :: test_gvec_kss(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),sq(3),q_bz(3),gamma_point(3,1)
 real(dp),pointer :: energies_p(:,:,:)
 real(dp),allocatable :: doccde(:),eigen(:),occfact(:),qlwl(:,:)
 type(Pawrhoij_type),allocatable :: Pawrhoij(:)
 type(vcoul_t) :: Vcp_ks

! *************************************************************************

 DBG_ENTER('COLL')

 ! Check for calculations that are not implemented
 ltest=ALL(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol)==Dtset%nband(1))
 ABI_CHECK(ltest,'Dtset%nband(:) must be constant')

 my_rank = xmpi_comm_rank(comm); nprocs  = xmpi_comm_size(comm)

 ! Basic parameters
 Sigp%ppmodel    = Dtset%ppmodel
 Sigp%gwcalctyp  = Dtset%gwcalctyp
 Sigp%nbnds      = Dtset%nband(1)
 Sigp%symsigma   = Dtset%symsigma
 Sigp%zcut       = Dtset%zcut
 Sigp%mbpt_sciss = Dtset%mbpt_sciss

 timrev=  2 ! This information is not reported in the header
            ! 1 => do not use time-reversal symmetry
            ! 2 => take advantage of time-reversal symmetry
 if (any(dtset%kptopt == [3, 4])) timrev = 1

 ! === For HF, SEX or COHSEX use Hybertsen-Louie PPM (only $\omega=0$) ===
 ! * Use fake screening for HF.
 ! FIXME Why, we should not redefine Sigp%ppmodel
 gwcalctyp=Sigp%gwcalctyp
 mod10 =MOD(Sigp%gwcalctyp,10)
 if (mod10==5.or.mod10==6.or.mod10==7) Sigp%ppmodel=2
 if (mod10<5.and.MOD(Sigp%gwcalctyp,1)/=1) then ! * One shot GW (PPM or contour deformation).
   if (Dtset%nomegasrd==1) then ! avoid division by zero!
     Sigp%nomegasrd  =1
     Sigp%maxomega4sd=zero
     Sigp%deltae     =zero
   else
     Sigp%nomegasrd   = Dtset%nomegasrd
     Sigp%maxomega4sd = Dtset%omegasrdmax
     Sigp%deltae     = (2*Sigp%maxomega4sd)/(Sigp%nomegasrd-1)
   endif
 else
   ! For AC no need to evaluate derivative by finite differences.
   Sigp%nomegasrd  =1
   Sigp%maxomega4sd=zero
   Sigp%deltae     =zero
 end if

 ! For analytic continuation define the number of imaginary frequencies for Sigma
 ! Tests show than more than 12 freqs in the Pade approximant worsen the results!
 Sigp%nomegasi=0

 if (mod10==1) then
   Sigp%nomegasi  =Dtset%nomegasi
   Sigp%omegasimax=Dtset%omegasimax
   Sigp%omegasimin=OMEGASIMIN
   write(msg,'(4a,i3,2(2a,f8.3),a)')ch10,&
&    ' Parameters for analytic continuation : ',ch10,&
&    '  number of imaginary frequencies for sigma =  ',Sigp%nomegasi,ch10,&
&    '  min frequency for sigma on imag axis [eV] =  ',Sigp%omegasimin*Ha_eV,ch10,&
&    '  max frequency for sigma on imag axis [eV] =  ',Sigp%omegasimax*Ha_eV,ch10
   call wrtout(std_out,msg,'COLL')

   !TODO this should not be done here but in init_sigma_t
   ABI_MALLOC(Sigp%omegasi,(Sigp%nomegasi))

   if (linear_imag_mesh) then  ! * Linear mesh along the imaginary axis.
     domegasi=Sigp%omegasimax/(Sigp%nomegasi-1)
     do io=1,Sigp%nomegasi
       Sigp%omegasi(io)=CMPLX(zero,(io-1)*domegasi)
     end do
   else ! * Logarithmic mesh along the imaginary axis.
     MSG_ERROR("AC + log mesh not implemented")
     !domegasi=(Sigp%omegasimax/Sigp%omegasimin)**(one/(Sigp%nomegasi-1))
     !Sigp%omegasi(1)=czero; ldi=domegasi
     !do io=2,Sigp%nomegasi
     ! omega(io)=CMPLX(zero,ldi*Sigp%omegasimin)
     ! Sigp%omegasi(io)=ldi*domegasi
     !end do
   end if

   write(msg,'(4a)')ch10,&
&    ' setup_sigma : calculating Sigma(iw)',&
&    ' at imaginary frequencies [eV] (Fermi Level set to 0) ',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   do io=1,Sigp%nomegasi
     write(msg,'(2(f10.3,2x))')Sigp%omegasi(io)*Ha_eV
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end do

   ltest=(Sigp%omegasimax>0.1d-4.and.Sigp%nomegasi>0)
   ABI_CHECK(ltest,'Wrong value of omegasimax or nomegasi')
   if (Sigp%gwcalctyp/=1) then ! only one shot GW is allowed for AC.
     MSG_ERROR("SC-GW with analytic continuation is not coded")
   end if
 end if

 if (Sigp%symsigma/=0.and.gwcalctyp>=20) then
   MSG_WARNING("SC-GW with symmetries is still under development. Use at your own risk!")
 end if

 ! Setup parameters for Spectral function.
 if (Dtset%gw_customnfreqsp/=0) then
   Sigp%nomegasr = Dtset%gw_customnfreqsp
   MSG_WARNING('Custom grid for spectral function specified. Assuming experienced user.')
   if (Dtset%gw_customnfreqsp/=0) then
     Dtset%nfreqsp = Dtset%gw_customnfreqsp
     MSG_WARNING('nfreqsp has been set to the same number as gw_customnfreqsp')
   end if
 else
   Sigp%nomegasr  =Dtset%nfreqsp
   Sigp%minomega_r=Dtset%freqspmin
   Sigp%maxomega_r=Dtset%freqspmax
 end if

 if (Sigp%nomegasr>0) then
   if (Dtset%gw_customnfreqsp==0) then
     ! Check
     if (Sigp%minomega_r >= Sigp%maxomega_r) then
       MSG_ERROR('freqspmin must be smaller than freqspmax!')
     end if
     if(Sigp%nomegasr==1) then
      domegas=0.d0
     else
      domegas=(Sigp%maxomega_r-Sigp%minomega_r)/(Sigp%nomegasr-1)
     endif
     !TODO this should be moved to Sr% and done in init_sigma_t
     ABI_MALLOC(Sigp%omega_r,(Sigp%nomegasr))
     do io=1,Sigp%nomegasr
       Sigp%omega_r(io) = CMPLX(Sigp%minomega_r + domegas*(io-1),zero)
     end do
     write(msg,'(4a,i8,3(2a,f8.3),a)')ch10,&
&      ' Parameters for the calculation of the spectral function : ',ch10,&
&      '  Number of points    = ',Sigp%nomegasr,ch10,&
&      '  Min frequency  [eV] = ',Sigp%minomega_r*Ha_eV,ch10,&
&      '  Max frequency  [eV] = ',Sigp%maxomega_r*Ha_eV,ch10,&
&      '  Frequency step [eV] = ',domegas*Ha_eV,ch10
     call wrtout(std_out,msg,'COLL')
   else
     Sigp%minomega_r = MINVAL(Dtset%gw_freqsp(:))
     Sigp%maxomega_r = MAXVAL(Dtset%gw_freqsp(:))
     !TODO this should be moved to Sr% and done in init_sigma_t
     ABI_MALLOC(Sigp%omega_r,(Sigp%nomegasr))
     do io=1,Sigp%nomegasr
       Sigp%omega_r(io) = CMPLX(Dtset%gw_freqsp(io),zero)
     end do
     write(msg,'(4a,i8,2(2a,f8.3),3a)')ch10,&
&      ' Parameters for the calculation of the spectral function : ',ch10,&
&      '  Number of points    = ',Sigp%nomegasr,ch10,&
&      '  Min frequency  [eV] = ',Sigp%minomega_r*Ha_eV,ch10,&
&      '  Max frequency  [eV] = ',Sigp%maxomega_r*Ha_eV,ch10,&
&      '  A custom set of frequencies is used! See the input file for values.',ch10
     call wrtout(std_out,msg,'COLL')
   end if
 else
   !In indefo all these quantities are set to zero
   !Sigp%nomegasr=1
   !allocate(Sigp%omega_r(Sigp%nomegasr))
   !Sigp%omega_r(1)=0
 end if

 ! Dimensional primitive translations rprimd (from input), gprimd, metrics and unit cell volume
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 Sigp%npwwfn=Dtset%npwwfn
 Sigp%npwx  =Dtset%npwsigx

 ! Read parameters of the WFK, verifify them and retrieve all G-vectors.
 call wfk_read_eigenvalues(wfk_fname,energies_p,Hdr_wfk,comm)
 mband = MAXVAL(Hdr_wfk%nband)

 remove_inv = .FALSE.
 call hdr_vs_dtset(Hdr_wfk,Dtset)

 test_npwkss = 0
 call make_gvec_kss(Dtset%nkpt,Dtset%kptns,Hdr_wfk%ecut_eff,Dtset%symmorphi,Dtset%nsym,Dtset%symrel,Dtset%tnons,&
&  gprimd,Dtset%prtvol,test_npwkss,test_gvec_kss,ierr)
 ABI_CHECK(ierr==0,"Fatal error in make_gvec_kss")

 ABI_MALLOC(gvec_kss,(3,test_npwkss))
 gvec_kss = test_gvec_kss
 ng_kss = test_npwkss

 ng = MIN(SIZE(gvec_kss,DIM=2),SIZE(test_gvec_kss,DIM=2))
 ierr = 0
 do ig1=1,ng
   if (ANY(gvec_kss(:,ig1)/=test_gvec_kss(:,ig1))) then
     ierr=ierr+1
     write(std_out,*)" gvec_kss ",ig1,"/",ng,gvec_kss(:,ig1),test_gvec_kss(:,ig1)
   end if
 end do
 ABI_CHECK(ierr == 0, "Mismatch between gvec_kss and test_gvec_kss")

 ABI_FREE(test_gvec_kss)

 ! Get important dimensions from the WFK header
 Sigp%nsppol =Hdr_wfk%nsppol
 Sigp%nspinor=Hdr_wfk%nspinor
 Sigp%nsig_ab=Hdr_wfk%nspinor**2  ! TODO Is it useful calculating only diagonal terms?

 if (Sigp%nbnds>mband) then
   Sigp%nbnds    =mband
   Dtset%nband(:)=mband
   Dtset%mband   =MAXVAL(Dtset%nband)
   write(msg,'(3a,i4,a)')&
&    'Number of bands found less then required',ch10,&
&    'calculation will proceed with nbnds = ',mband,ch10
   MSG_WARNING(msg)
 end if

 ! Check input
 if (Sigp%ppmodel==3.or.Sigp%ppmodel==4) then
   if (gwcalctyp>=10) then
     write(msg,'(a,i3,a)')' The ppmodel chosen and gwcalctyp ',Dtset%gwcalctyp,' are not compatible. '
     MSG_ERROR(msg)
   end if
   if (Sigp%nspinor==2) then
     write(msg,'(a,i3,a)')' The ppmodel chosen and nspinor ',Sigp%nspinor,' are not compatible. '
     MSG_ERROR(msg)
   end if
 end if

 ! Create crystal_t data type
 call crystal_from_hdr(Cryst,Hdr_wfk,timrev,remove_inv)
 call crystal_print(Cryst)

 if (Sigp%npwwfn>ng_kss) then ! cannot use more G"s for the wfs than those stored on file
   Sigp%npwwfn  =ng_kss
   Dtset%npwwfn =ng_kss
   write(msg,'(2a,(a,i8,a))')&
&    'Number of G-vectors for WFS found in the KSS file is less than required',ch10,&
&    'calculation will proceed with npwwfn  = ',Sigp%npwwfn,ch10
   MSG_WARNING(msg)
 end if

 if (Sigp%npwx>ng_kss) then
   ! Have to recalcuate the (large) sphere for Sigma_x.
   pinv=1; if (remove_inv.and.Cryst%timrev==2) pinv=-1
   gamma_point(:,1) = (/zero,zero,zero/); nullify(gsphere_sigx_p)

   call merge_and_sort_kg(1,gamma_point,Dtset%ecutsigx,Cryst%nsym,pinv,Cryst%symrel,&
&   Cryst%gprimd,gsphere_sigx_p,Dtset%prtvol)

   ng_sigx=SIZE(gsphere_sigx_p,DIM=2)
   Sigp%npwx     = ng_sigx
   Dtset%npwsigx = ng_sigx

   write(msg,'(2a,(a,i8,a))')&
&    'Number of G-vectors for Sigma_x found in the KSS file is less than required',ch10,&
&    'calculation will proceed with npwsigx = ',Sigp%npwx,ch10
   MSG_WARNING(msg)

   ltest = (Sigp%npwx >= ng_kss)
   ABI_CHECK(ltest,"Sigp%npwx<ng_kss!")

   ! * Fill gvec_kss with larger sphere.
   ABI_FREE(gvec_kss)
   ABI_MALLOC(gvec_kss,(3,Sigp%npwx))
   gvec_kss = gsphere_sigx_p
   ABI_FREE(gsphere_sigx_p)
 end if

 ! Set up of the k-points and tables in the whole BZ ===
 ! TODO Recheck symmorphy and inversion
 call kmesh_init(Kmesh,Cryst,Hdr_wfk%nkpt,Hdr_wfk%kptns,Dtset%kptopt,wrap_1zone=.FALSE.)
 !call kmesh_init(Kmesh,Cryst,Hdr_wfk%nkpt,Hdr_wfk%kptns,Dtset%kptopt,wrap_1zone=.TRUE.)

 ! Some required information are not filled up inside kmesh_init
 ! So doing it here, even though it is not clean
 Kmesh%kptrlatt(:,:) =Dtset%kptrlatt(:,:)
 Kmesh%nshift        =Dtset%nshiftk
 ABI_ALLOCATE(Kmesh%shift,(3,Kmesh%nshift))
 Kmesh%shift(:,:)    =Dtset%shiftk(:,1:Dtset%nshiftk)

 call kmesh_print(Kmesh,"K-mesh for the wavefunctions",std_out,Dtset%prtvol,"COLL")
 call kmesh_print(Kmesh,"K-mesh for the wavefunctions",ab_out, 0,           "COLL")

 ! === Initialize the band structure datatype ===
 ! * Copy WFK energies and occupations up to Sigp%nbnds==Dtset%nband(:)
 bantot=SUM(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol))
 ABI_MALLOC(doccde,(bantot))
 ABI_MALLOC(eigen,(bantot))
 ABI_MALLOC(occfact,(bantot))
 doccde(:)=zero; eigen(:)=zero; occfact(:)=zero

 jj=0; ibtot=0
 do isppol=1,Dtset%nsppol
   do ikibz=1,Dtset%nkpt
     do ib=1,Hdr_wfk%nband(ikibz+(isppol-1)*Dtset%nkpt)
       ibtot=ibtot+1
       if (ib<=Sigp%nbnds) then
         jj=jj+1
         occfact(jj)=Hdr_wfk%occ(ibtot)
         eigen  (jj)=energies_p(ib,ikibz,isppol)
       end if
     end do
   end do
 end do
 ABI_FREE(energies_p)

 ! Make sure that Dtset%wtk==Kmesh%wt due to the dirty treatment of
 ! symmetry operations in the old GW code (symmorphy and inversion)
 ltest=(ALL(ABS(Dtset%wtk(1:Kmesh%nibz)-Kmesh%wt(1:Kmesh%nibz))<tol6))
 ABI_CHECK(ltest,'Mismatch between Dtset%wtk and Kmesh%wt')

 ABI_MALLOC(npwarr,(Dtset%nkpt))
 npwarr(:)=Sigp%npwwfn

 call ebands_init(bantot,KS_BSt,Dtset%nelect,doccde,eigen,Dtset%istwfk,Kmesh%ibz,Dtset%nband,&
&  Kmesh%nibz,npwarr,Dtset%nsppol,Dtset%nspinor,Dtset%tphysel,Dtset%tsmear,Dtset%occopt,occfact,Kmesh%wt,&
&  dtset%charge, dtset%kptopt, dtset%kptrlatt_orig, dtset%nshiftk_orig, dtset%shiftk_orig,&
&  dtset%kptrlatt, dtset%nshiftk, dtset%shiftk)

 ABI_FREE(doccde)
 ABI_FREE(eigen)
 ABI_FREE(npwarr)

 ! Calculate KS occupation numbers and ks_vbk(nkibz,nsppol) ====
 ! ks_vbk gives the (valence|last Fermi band) index for each k and spin.
 ! spinmagntarget is passed to fermi.F90 to fix the problem with newocc in case of magnetic metals
 call ebands_update_occ(KS_BSt,Dtset%spinmagntarget,prtvol=0)

 gap_err = get_gaps(KS_BSt, gaps)
 call gaps_print(gaps, unit=std_out)
 call ebands_report_gap(KS_BSt, unit=std_out)

 ABI_MALLOC(val_indeces,(KS_BSt%nkpt,KS_BSt%nsppol))
 val_indeces = get_valence_idx(KS_BSt)

 ! Create Sigma header
 ! TODO Fix problems with symmorphy and k-points
 call hdr_init(KS_BSt,codvsn,Dtset,Hdr_out,Pawtab,pertcase0,Psps,wvl)

 ! Get Pawrhoij from the header of the WFK file
 ABI_DT_MALLOC(Pawrhoij,(Cryst%natom*Dtset%usepaw))
 if (Dtset%usepaw==1) then
   call pawrhoij_alloc(Pawrhoij,1,Dtset%nspden,Dtset%nspinor,Dtset%nsppol,Cryst%typat,pawtab=Pawtab)
   call pawrhoij_copy(Hdr_wfk%Pawrhoij,Pawrhoij)
 end if
 call hdr_update(hdr_out,bantot,1.0d20,1.0d20,1.0d20,Cryst%rprimd,occfact,Pawrhoij,Cryst%xred,dtset%amu_orig(:,1))

 ABI_FREE(occfact)
 call pawrhoij_free(Pawrhoij)
 ABI_DT_FREE(Pawrhoij)

 ! ===========================================================
 ! ==== Setup of k-points and bands for the GW corrections ====
 ! ===========================================================
 ! * maxbdgw and minbdgw are the Max and min band index for GW corrections over k-points.
 !   They are used to dimension the wavefunctions and to calculate the matrix elements.
 !
 if (dtset%nkptgw==0) then
   ! Use qp_range to select the interesting k-points and the corresponing bands.
   !
   !    0 --> Compute the QP corrections only for the fundamental and the optical gap.
   ! +num --> Compute the QP corrections for all the k-points in the irreducible zone and include `num`
   !           bands above and below the Fermi level.
   ! -num --> Compute the QP corrections for all the k-points in the irreducible zone.
   !          Include all occupied states and `num` empty states.

   call wrtout(std_out, "nkptgw == 0 ==> Automatic selection of k-points and bands for the corrections.")
   if (gap_err /=0 .and. dtset%gw_qprange==0) then
     msg = "Problem while computing the fundamental and optical gap (likely metal). Will replace gw_qprange=0 with gw_qprange=1"
     MSG_WARNING(msg)
     dtset%gw_qprange = 1
   end if
   gw_qprange = dtset%gw_qprange

   if (dtset%ucrpa>0) then
     dtset%nkptgw=Kmesh%nbz
     Sigp%nkptgw =dtset%nkptgw
     ABI_MALLOC(Sigp%kptgw,(3,Sigp%nkptgw))
     ABI_MALLOC(Sigp%minbnd,(Sigp%nkptgw,Sigp%nsppol))
     ABI_MALLOC(Sigp%maxbnd,(Sigp%nkptgw,Sigp%nsppol))
     Sigp%kptgw(:,:)=Kmesh%bz(:,:)
     Sigp%minbnd=1
     Sigp%maxbnd=Sigp%nbnds

   else if (gw_qprange/=0) then
     ! Include all the k-points in the IBZ.
     ! Note that kptgw == ebands%kptns so we can use a single ik index in the loop over k-points
     ! No need to map kptgw onto ebands%kptns.
     dtset%nkptgw=Kmesh%nibz
     Sigp%nkptgw =dtset%nkptgw
     ABI_MALLOC(Sigp%kptgw,(3,Sigp%nkptgw))
     ABI_MALLOC(Sigp%minbnd,(Sigp%nkptgw,Sigp%nsppol))
     ABI_MALLOC(Sigp%maxbnd,(Sigp%nkptgw,Sigp%nsppol))
     Sigp%kptgw(:,:)=Kmesh%ibz(:,:)
     Sigp%minbnd=1
     Sigp%maxbnd=Sigp%nbnds

     if (gw_qprange>0) then
       ! All k-points: Add buffer of bands above and below the Fermi level.
       do spin=1,Sigp%nsppol
         do ik=1,Sigp%nkptgw
           Sigp%minbnd(ik,spin) = MAX(val_indeces(ik,spin) - gw_qprange, 1)
           Sigp%maxbnd(ik,spin) = MIN(val_indeces(ik,spin) + gw_qprange + 1, Sigp%nbnds)
         end do
       end do

     else
       ! All k-points: include all occupied states and -gw_qprange empty states.
       Sigp%minbnd = 1
       do spin=1,Sigp%nsppol
         do ik=1,Sigp%nkptgw
           Sigp%maxbnd(ik,spin) = MIN(val_indeces(ik,spin) - gw_qprange, Sigp%nbnds)
         end do
       end do
     end if

   else
     ! gw_qprange is not specified in the input.
     ! Include the optical and the fundamental KS gap.
     ! The main problem here is that kptgw and nkptgw do not depend on the spin and therefore
     ! we have compute the union of the k-points where the fundamental and the optical gaps are located.
     !
     ! Find the list of `interesting` kpoints.
     ABI_CHECK(gap_err == 0, "gw_qprange 0 cannot be used because I cannot find the gap (gap_err !=0)")
     nk_found = 1; kpos(1) = gaps%fo_kpos(1,1)

     do spin=1,Sigp%nsppol
       do ifo=1,3
         ik = gaps%fo_kpos(ifo, spin)
         found = .FALSE.; jj = 0
         do while (.not. found .and. jj < nk_found)
           jj = jj + 1; found = (kpos(jj) == ik)
         end do
         if (.not. found) then
           nk_found = nk_found + 1; kpos(nk_found) = ik
         end if
       end do
     end do

     ! Now we can define the list of k-points and the bands range.
     dtset%nkptgw=nk_found
     Sigp%nkptgw =dtset%nkptgw

     ABI_MALLOC(Sigp%kptgw,(3,Sigp%nkptgw))
     ABI_MALLOC(Sigp%minbnd,(Sigp%nkptgw,Sigp%nsppol))
     ABI_MALLOC(Sigp%maxbnd,(Sigp%nkptgw,Sigp%nsppol))

     do ii=1,Sigp%nkptgw
       ik = kpos(ii)
       Sigp%kptgw(:,ii)=Kmesh%ibz(:,ik)
       do spin=1,Sigp%nsppol
         Sigp%minbnd(ii,spin) = val_indeces(ik,spin)
         Sigp%maxbnd(ii,spin) = val_indeces(ik,spin) + 1
       end do
     end do
   end if

 else
   ! Treat only the k-points and bands specified in the input file.
   Sigp%nkptgw=dtset%nkptgw
   ABI_MALLOC(Sigp%kptgw,(3,Sigp%nkptgw))
   ABI_MALLOC(Sigp%minbnd,(Sigp%nkptgw,Sigp%nsppol))
   ABI_MALLOC(Sigp%maxbnd,(Sigp%nkptgw,Sigp%nsppol))

   do spin=1,Sigp%nsppol
     Sigp%minbnd(:,spin)=dtset%bdgw(1,:,spin)
     Sigp%maxbnd(:,spin)=dtset%bdgw(2,:,spin)
   end do

   do ii=1,3
     do ikcalc=1,Sigp%nkptgw
       Sigp%kptgw(ii,ikcalc)=Dtset%kptgw(ii,ikcalc)
     end do
   end do

   do spin=1,Sigp%nsppol
     do ikcalc=1,Sigp%nkptgw
       if (Dtset%bdgw(2,ikcalc,spin)>Sigp%nbnds) then
         write(msg,'(a,2i0,2(a,i0),2a,i0)')&
&          "For (k,s) ",ikcalc,spin," bdgw= ",Dtset%bdgw(2,ikcalc,spin), " > nbnds=",Sigp%nbnds,ch10,&
&          "Calculation will continue with bdgw =",Sigp%nbnds
         MSG_COMMENT(msg)
         Dtset%bdgw(2,ikcalc,spin)=Sigp%nbnds
       end if
     end do
   end do

 end if

 ! Make sure that all the degenerate states are included.
 ! * We will have to average the GW corrections over degenerate states if symsigma=1 is used.
 ! * KS states belonging to the same irreducible representation should be included in the basis set used for SCGW.
 if (Sigp%symsigma/=0 .or. gwcalctyp>=10) then
   do isppol=1,Sigp%nsppol
     do ikcalc=1,Sigp%nkptgw

       if (has_IBZ_item(Kmesh,Sigp%kptgw(:,ikcalc),ikibz,G0)) then
         call enclose_degbands(KS_BSt,ikibz,isppol,Sigp%minbnd(ikcalc,isppol),Sigp%maxbnd(ikcalc,isppol),changed,tol_enediff)
         if (changed) then
           write(msg,'(2(a,i0),2a,2(1x,i0))')&
&            "Not all the degenerate states at ikcalc= ",ikcalc,", spin= ",isppol,ch10,&
&            "were included in the bdgw set. bdgw has been changed to: ",Sigp%minbnd(ikcalc,isppol),Sigp%maxbnd(ikcalc,isppol)
           MSG_COMMENT(msg)
         end if
       else
         MSG_ERROR(sjoin('k-point', ktoa(Sigp%kptgw(:,ikcalc)), 'not in IBZ'))
       end if

     end do
   end do
 end if

 !if (.not. associated(Dtset%bdgw)) then
 !  ABI_MALLOC(Dtset%bdgw, (2,Sigp%nkptgw,Sigp%nsppol))
 !end if
 !do spin=1,Sigp%nsppol
 !  Dtset%bdgw(1,:,spin) = Sigp%minbnd(:,spin)
 !  Dtset%bdgw(2,:,spin) = Sigp%maxbnd(:,spin)
 !end do

 Sigp%minbdgw=MINVAL(Sigp%minbnd)
 Sigp%maxbdgw=MAXVAL(Sigp%maxbnd)

 ABI_MALLOC(Sigp%kptgw2bz,(Sigp%nkptgw))
 !
 !=== Check if the k-points are in the BZ ===
 !FB TODO Honestly the code is not able to treat k-points, which are not in the IBZ.
 !This extension should require to change the code in different places.
 !Therefore, one should by now prevent the user from calculating sigma for a k-point not in the IBZ.

 do ikcalc=1,Sigp%nkptgw
   if (has_BZ_item(Kmesh,Sigp%kptgw(:,ikcalc),ikcalc2bz,G0)) then
     !found = has_IBZ_item(Kmesh,Sigp%kptgw(:,ikcalc),ikcalc2bz,G0)
     Sigp%kptgw2bz(ikcalc) = ikcalc2bz
   else
     MSG_ERROR(sjoin('k-point:', ktoa(Sigp%kptgw(:,ikcalc)), 'not in the kbz set'))
   end if
 end do

 ! Check if there are duplicated k-point in Sigp%
 do ii=1,Sigp%nkptgw
   do jj=ii+1,Sigp%nkptgw
     if (isamek(Sigp%kptgw(:,ii),Sigp%kptgw(:,jj),G0)) then
       write(msg,'(5a)')&
&        'kptgw contains duplicated k-points. This is not allowed since ',ch10,&
&        'the QP corrections for this k-point will be calculated more than once. ',ch10,&
&        'Check your input file. '
       MSG_ERROR(msg)
     end if
   end do
 end do
 !
 ! Warn the user if SCGW run and not all the k-points are included.
 if (gwcalctyp>=10 .and. Sigp%nkptgw/=Hdr_wfk%nkpt) then
   write(msg,'(3a,2(a,i0),2a)')ch10,&
&    " COMMENT: In a self-consistent GW run, the QP corrections should be calculated for all the k-points of the KSS file ",ch10,&
&    " but nkptgw= ",Sigp%nkptgw," and WFK nkpt= ",Hdr_wfk%nkpt,ch10,&
&    " Assuming expert user. Execution will continue. "
   call wrtout(ab_out,msg,"COLL")
 end if

 ! Setup of the table used in the case of SCGW on wavefunctions to reduce the number
 ! of elements <i,kgw,s|\Sigma|j,kgw,s> that have to be calculated. No use of symmetries, except for Hermiticity.
 call sigma_tables(Sigp,Kmesh)

 ! === Read external file and initialize basic dimension of Er% ===
 ! TODO use mqmem as input variable instead of gwmem

 ! === If required, use a matrix for $\Sigma_c$ which is smaller than that stored on file ===
 ! * By default the entire matrix is read and used,
 ! * Define consistently npweps and ecuteps for \Sigma_c according the input
 if (Dtset%npweps>0.or.Dtset%ecuteps>0) then
   ! This should not happen: the Dtset array should not be modified after having been initialized.
   if (Dtset%npweps>0) Dtset%ecuteps=zero
   nsheps=0
   call setshells(Dtset%ecuteps,Dtset%npweps,nsheps,Dtset%nsym,gmet,gprimd,Dtset%symrel,'eps',ucvol)
 end if

 mqmem=0; if (Dtset%gwmem/10==1) mqmem=1

 if (Dtset%getscr/=0.or.Dtset%irdscr/=0) then
   fname=Dtfil%fnameabi_scr
 else if (Dtset%getsuscep/=0.or.Dtset%irdsuscep/=0) then
   fname=Dtfil%fnameabi_sus
 else
   fname=Dtfil%fnameabi_scr
   !FIXME this has to be cleaned, in tgw2_3 Dtset%get* and Dtset%ird* are  not defined
   !MSG_ERROR("getsuscep or irdsuscep are not defined")
 end if
 !
 ! === Setup of q-mesh in the whole BZ ===
 ! * Stop if a nonzero umklapp is needed to reconstruct the BZ. In this case, indeed,
 !   epsilon^-1(Sq) should be symmetrized in csigme using a different expression (G-G_o is needed)
 !
 if (sigma_needs_w(Sigp)) then
   if (.not. file_exists(fname)) then
     fname = nctk_ncify(fname)
     MSG_COMMENT(sjoin("File not found. Will try netcdf file:", fname))
   end if

   call init_Er_from_file(Er,fname,mqmem,Dtset%npweps,comm)

   Sigp%npwc=Er%npwe
   if (Sigp%npwc>Sigp%npwx) then
     Sigp%npwc=Sigp%npwx
     MSG_COMMENT("Found npw_correlation > npw_exchange, Imposing npwc=npwx")
     ! There is a good reason for doing so, see csigme.F90 and the size of the arrays
     ! rhotwgp and rhotwgp: we need to define a max size and we opt for Sigp%npwx.
   end if
   Er%npwe=Sigp%npwc
   Dtset%npweps=Er%npwe
   call kmesh_init(Qmesh,Cryst,Er%nqibz,Er%qibz,Dtset%kptopt)

 else
   Er%npwe     =1
   Sigp%npwc   =1
   Dtset%npweps=1
   call find_qmesh(Qmesh,Cryst,Kmesh)
   ABI_MALLOC(Er%gvec,(3,1))
   Er%gvec(:,1) = (/0,0,0/)
 end if

 call kmesh_print(Qmesh,"Q-mesh for screening function",std_out,Dtset%prtvol,"COLL")
 call kmesh_print(Qmesh,"Q-mesh for screening function",ab_out ,0           ,"COLL")

 do iqbz=1,Qmesh%nbz
   call get_BZ_item(Qmesh,iqbz,q_bz,iq_ibz,isym,itim,umklp=q_umklp)

   if (ANY(q_umklp/=0)) then
     sq = (3-2*itim)*MATMUL(Cryst%symrec(:,:,isym),Qmesh%ibz(:,iq_ibz))
     write(std_out,*) sq,Qmesh%bz(:,iqbz)
     write(msg,'(a,3f6.3,a,3f6.3,2a,9i3,a,i2,2a)')&
&      'qpoint ',Qmesh%bz(:,iqbz),' is the symmetric of ',Qmesh%ibz(:,iq_ibz),ch10,&
&      'through operation ',Cryst%symrec(:,:,isym),' and itim ',itim,ch10,&
&      'however a non zero umklapp G_o vector is required and this is not yet allowed'
     MSG_ERROR(msg)
   end if
 end do
 !
 ! === Find optimal value for G-sphere enlargment due to oscillator matrix elements ===
 ! * Here I have to be sure that Qmesh%bz is always inside the BZ, not always true size bz is buggy
 ! * -one is used because we loop over all the possibile differences, unlike screening

 call get_ng0sh(Sigp%nkptgw,Sigp%kptgw,Kmesh%nbz,Kmesh%bz,Qmesh%nbz,Qmesh%bz,-one,ng0sh_opt)
 call wrtout(std_out, sjoin(' Optimal value for ng0sh ', ltoa(ng0sh_opt)), "COLL")
 Sigp%mG0=ng0sh_opt

 ! G-sphere for W and Sigma_c is initialized from the SCR file.
 call gsph_init(Gsph_c,Cryst,Er%npwe,gvec=Er%gvec)
 call gsph_init(Gsph_x,Cryst,Sigp%npwx,gvec=gvec_kss)

 ! === Make biggest G-sphere of Sigp%npwvec vectors ===
 Sigp%npwvec=MAX(Sigp%npwwfn,Sigp%npwx)
 call gsph_init(Gsph_Max,Cryst,Sigp%npwvec,gvec=gvec_kss)

!BEGINDEBUG
 ! Make sure that the two G-spheres are equivalent.
 ierr=0
 if (sigma_needs_w(Sigp)) then
   ng = MIN(SIZE(Gsph_c%gvec,DIM=2),SIZE(gvec_kss,DIM=2))
   do ig1=1,ng
     if (ANY(Gsph_c%gvec(:,ig1)/=gvec_kss(:,ig1))) then
       ierr=ierr+1
       write(std_out,*)" Gsph_c, gvec_kss ",ig1,"/",ng,Gsph_c%gvec(:,ig1),gvec_kss(:,ig1)
     end if
   end do
   ABI_CHECK(ierr==0,"Mismatch between Gsph_c and gvec_kss")
 end if
 ierr=0
 ng = MIN(SIZE(Gsph_x%gvec,DIM=2),SIZE(gvec_kss,DIM=2))
 do ig1=1,ng
   if (ANY(Gsph_x%gvec(:,ig1)/=gvec_kss(:,ig1))) then
     ierr=ierr+1
     write(std_out,*)" Gsph_x, gvec_kss ",ig1,"/",ng,Gsph_x%gvec(:,ig1),gvec_kss(:,ig1)
   end if
 end do
 ABI_CHECK(ierr==0,"Mismatch between Gsph_x and gvec_kss")
!ENDDEBUG

 ABI_FREE(gvec_kss)
 !
 ! === Get Fourier components of the Coulombian for all q-points in the IBZ ===
 ! * If required, use a cutoff in the interaction
 ! * Pcv%vc_sqrt contains Vc^{-1/2}
 ! * Setup also the analytical calculation of the q->0 component
 ! FIXME recheck ngfftf since I got different charge outside the cutoff region

 if (Dtset%gw_nqlwl==0) then
   nqlwl=1
   ABI_MALLOC(qlwl,(3,nqlwl))
   qlwl(:,1)= GW_Q0_DEFAULT
 else
   nqlwl=Dtset%gw_nqlwl
   ABI_MALLOC(qlwl,(3,nqlwl))
   qlwl(:,:)=Dtset%gw_qlwl(:,1:nqlwl)
 end if

!The Coulomb interaction used here might have two terms : 
!the first term generates the usual sigma self-energy, but possibly, one should subtract
!from it the Coulomb interaction already present in the Kohn-Sham basis, 
!if the usefock associated to ixc is one.
!The latter excludes (in the present implementation) mod(Dtset%gwcalctyp,10)==5
 call get_xclevel(Dtset%ixc,xclevel_ixc,usefock_ixc)
 if(usefock_ixc==1)then
   if(mod(Dtset%gwcalctyp,10)==5)then
     write(msg,'(4a,i3,a,i3,4a,i5)')ch10,&
&     ' The starting wavefunctions were obtained from self-consistent calculations in the planewave basis set',ch10,&
&     ' with ixc = ',Dtset%ixc,' associated with usefock =',usefock_ixc,ch10,&
&     ' In this case, the present implementation does not allow that the self-energy for sigma corresponds to',ch10,&
&     '  mod(gwcalctyp,10)==5, while your gwcalctyp= ',Dtset%gwcalctyp
     MSG_ERROR(msg)
   endif
 endif

 nvcoul_init=1
 if(usefock_ixc==1)nvcoul_init=2

 do ivcoul_init=1,nvcoul_init
   rcut = Dtset%rcut
   icutcoul_eff=Dtset%icutcoul
   Sigp%sigma_mixing=one
   if( ((mod(Dtset%gwcalctyp,10)==5).and.ivcoul_init==1) .or. ivcoul_init==2)then
     if(abs(Dtset%hyb_mixing)>tol8)then
!      Warning : the absolute value is needed, because of the singular way used to define the default for this input variable
       Sigp%sigma_mixing=abs(Dtset%hyb_mixing)
     else if(abs(Dtset%hyb_mixing_sr)>tol8)then
       Sigp%sigma_mixing=abs(Dtset%hyb_mixing_sr)
       icutcoul_eff=5
     endif
     if(abs(rcut)<tol6 .and. abs(Dtset%hyb_range_fock)>tol8)rcut=one/Dtset%hyb_range_fock
   endif

!#if 1
   if(ivcoul_init==1)then
     if (Gsph_x%ng > Gsph_c%ng) then
       call vcoul_init(Vcp,Gsph_x,Cryst,Qmesh,Kmesh,rcut,icutcoul_eff,Dtset%vcutgeo,&
&        Dtset%ecutsigx,Gsph_x%ng,nqlwl,qlwl,ngfftf,comm)
     else
       call vcoul_init(Vcp,Gsph_c,Cryst,Qmesh,Kmesh,rcut,icutcoul_eff,Dtset%vcutgeo,&
&        Dtset%ecutsigx,Gsph_c%ng,nqlwl,qlwl,ngfftf,comm)
     end if
   else
!    Use a temporary Vcp_ks to compute the Coulomb interaction already present in the Fock part of the Kohn-Sham Hamiltonian
     if (Gsph_x%ng > Gsph_c%ng) then
       call vcoul_init(Vcp_ks,Gsph_x,Cryst,Qmesh,Kmesh,rcut,icutcoul_eff,Dtset%vcutgeo,&
&        Dtset%ecutsigx,Gsph_x%ng,nqlwl,qlwl,ngfftf,comm)
     else
       call vcoul_init(Vcp_ks,Gsph_c,Cryst,Qmesh,Kmesh,rcut,icutcoul_eff,Dtset%vcutgeo,&
&        Dtset%ecutsigx,Gsph_c%ng,nqlwl,qlwl,ngfftf,comm)
     end if
!    Now compute the residual Coulomb interaction
     Vcp%vc_sqrt_resid=sqrt(Vcp%vc_sqrt**2-Vcp_ks%vc_sqrt**2)
     call vcoul_free(Vcp_ks)
   endif
!#else
!   call vcoul_init(Vcp,Gsph_Max,Cryst,Qmesh,Kmesh,rcut,icutcoul_eff,ivcoul_init,Dtset%vcutgeo,&
!&    Dtset%ecutsigx,Sigp%npwx,nqlwl,qlwl,ngfftf,comm)
!#endif

 enddo

#if 0
 ! Using the random q for the optical limit is one of the reasons
 ! why sigma breaks the initial energy degeneracies.
 Vcp%i_sz=zero
 Vcp%vc_sqrt(1,:)=czero
 Vcp%vcqlwl_sqrt(1,:)=czero
#endif

 ABI_FREE(qlwl)

 Sigp%ecuteps = Dtset%ecuteps
 Sigp%ecutwfn = Dtset%ecutwfn
 Sigp%ecutsigx = Dtset%ecutsigx

 ! === Setup of the FFT mesh for the oscilator strengths ===
 ! * gwc_ngfft(7:18)==Dtset%ngfft(7:18) which is initialized before entering screening.
 ! * Here we redefine gwc_ngfft(1:6) according to the following options :
 !
 ! method==0 --> FFT grid read from fft.in (debugging purpose)
 ! method==1 --> Normal FFT mesh
 ! method==2 --> Slightly augmented FFT grid to calculate exactly rho_tw_g (see setmesh.F90)
 ! method==3 --> Doubled FFT grid, same as the the FFT for the density,
 !
 ! enforce_sym==1 ==> Enforce a FFT mesh compatible with all the symmetry operation and FFT library
 ! enforce_sym==0 ==> Find the smallest FFT grid compatbile with the library, do not care about symmetries
 !
 gwc_ngfft(1:18)=Dtset%ngfft(1:18)
 gwx_ngfft(1:18)=Dtset%ngfft(1:18)

 method=2
 if (Dtset%fftgw==00 .or. Dtset%fftgw==01) method=0
 if (Dtset%fftgw==10 .or. Dtset%fftgw==11) method=1
 if (Dtset%fftgw==20 .or. Dtset%fftgw==21) method=2
 if (Dtset%fftgw==30 .or. Dtset%fftgw==31) method=3
 enforce_sym=MOD(Dtset%fftgw,10)

 ! FFT mesh for sigma_x.
 call setmesh(gmet,Gsph_Max%gvec,gwx_ngfft,Sigp%npwvec,Sigp%npwx,Sigp%npwwfn,&
&  gwx_nfftot,method,Sigp%mG0,Cryst,enforce_sym)

 ! FFT mesh for sigma_c.
 call setmesh(gmet,Gsph_Max%gvec,gwc_ngfft,Sigp%npwvec,Er%npwe,Sigp%npwwfn,&
&  gwc_nfftot,method,Sigp%mG0,Cryst,enforce_sym,unit=dev_null)

 !call new_setmesh(Cryst,ecut_osc,ecutwfn,nkpt,kpoints,method,Sigp%mG0,enforce_sym,gwx_ngfft,gwx_nfftot)
 !call new_setmesh(Cryst,ecut_osc,ecutwfn,nkpt,kpoints,method,Sigp%mG0,enforce_sym,gwc_ngfft,gwc_nfftot)

 ! ======================================================================
 ! ==== Check for presence of files with core orbitals, for PAW only ====
 ! ======================================================================
 Sigp%use_sigxcore=0
 if (Dtset%usepaw==1.and.Dtset%gw_sigxcore==1) then
   ii = 0
   do itypat=1,Cryst%ntypat
     string = Psps%filpsp(itypat)
     fcore = "CORE_"//TRIM(basename(string))
     ic = INDEX (TRIM(string), "/" , back=.TRUE.)
     if (ic>0 .and. ic<LEN_TRIM(string)) then
       ! string defines a path, prepend path to fcore
       fcore = Psps%filpsp(itypat)(1:ic)//TRIM(fcore)
     end if
     if (file_exists(fcore)) then
       ii = ii+1
     else
       MSG_WARNING(sjoin("HF decoupling is required but cannot find file:", fcore))
     end if
   end do

   Sigp%use_sigxcore=1
   if (ii/=Cryst%ntypat) then
     MSG_ERROR("Files with core orbitals not found")
   end if
 end if ! PAW+HF decoupling
 !
 ! ==============================
 ! ==== Extrapolar technique ====
 ! ==============================
 Sigp%gwcomp   = Dtset%gwcomp
 Sigp%gwencomp = Dtset%gwencomp

 if (Sigp%gwcomp==1) then
   write(msg,'(6a,e11.4,a)')ch10,&
&    'Using the extrapolar approximation to accelerate convergence',ch10,&
&    'with respect to the number of bands included',ch10,&
&    'with gwencomp: ',Sigp%gwencomp*Ha_eV,' [eV]'
   call wrtout(std_out,msg,'COLL')
 end if
 !
 ! ===================================
 ! ==== Final compatibility tests ====
 ! ===================================
 if (ANY(KS_BSt%istwfk/=1)) then
   MSG_WARNING('istwfk/=1 is still under development')
 end if

 ltest=(KS_BSt%mband==Sigp%nbnds.and.ALL(KS_BSt%nband==Sigp%nbnds))
 ABI_CHECK(ltest,'BUG in definition of KS_BSt%nband')

 ! FIXME
 if (Dtset%symsigma/=0 .and. Sigp%nomegasr/=0) then
   if (idx_spatial_inversion(Cryst) == 0) then
     write(msg,'(5a)')' setup_sigma : BUG :',ch10,&
&      'It is not possible to use symsigma/=0 to calculate the spectral function ',ch10,&
&      'when the system does not have the spatial inversion. Please use symsigma=0 '
     MSG_WARNING(msg)
   end if
 end if

 if (mod10==SIG_GW_AC) then
   if (Sigp%gwcalctyp/=1) MSG_ERROR("Self-consistency with AC not implemented")
   if (Sigp%gwcomp==1) MSG_ERROR("AC with extrapolar technique not implemented")
 end if

 call gaps_free(gaps)

 ABI_FREE(val_indeces)

 DBG_EXIT('COLL')

end subroutine setup_sigma
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/sigma_tables
!! NAME
!! sigma_tables
!!
!! FUNCTION
!!  Build symmetry tables used to speedup self-consistent GW calculations.
!!
!! INPUTS
!! Kmesh <kmesh_t>=Structure describing the k-point sampling.
!! [Bnd_sym(Kmesh%nibz,Sigp%nsppol)] <type(Bands_Symmetries)>
!!
!! SiDE EFFECTS
!! Sigp<sigparams_t>=This routine initializes the tables:
!!   %Sigcij_tab
!!   %Sigxij_tab
!!  that are used to select the matrix elements of the self-energy that have to be calculated.
!!
!! PARENTS
!!      setup_sigma,sigma
!!
!! CHILDREN
!!
!! SOURCE


subroutine sigma_tables(Sigp,Kmesh,Bnd_sym)

 use defs_basis
 use defs_datatypes
 use m_gwdefs
 use m_profiling_abi
 use m_errors

 use m_bz_mesh,     only : kmesh_t
 use m_esymm,       only : esymm_t, esymm_failed

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigma_tables'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(sigparams_t),intent(inout) :: Sigp
 type(kmesh_t),intent(in) :: Kmesh
!arrays
 type(esymm_t),optional,intent(in) :: Bnd_sym(Kmesh%nibz,Sigp%nsppol)

!Local variables-------------------------------
!scalars
 integer :: gwcalctyp,spin,ikcalc,ik_ibz,bmin,bmax,bcol,brow
 integer :: ii,idx_x,idx_c,irr_idx1,irr_idx2
!arrays
 integer,allocatable :: sigc_bidx(:),sigx_bidx(:)
 logical :: use_sym_at(Kmesh%nibz,Sigp%nsppol)

! *************************************************************************

 gwcalctyp=Sigp%gwcalctyp

 ! Recreate the Sig_ij tables taking advantage of the classification of the bands.
 if (allocated(Sigp%Sigxij_tab)) then
   call sigijtab_free(Sigp%Sigxij_tab)
   ABI_DT_FREE(Sigp%Sigxij_tab)
 end if
 if (allocated(Sigp%Sigcij_tab)) then
   call sigijtab_free(Sigp%Sigcij_tab)
   ABI_DT_FREE(Sigp%Sigcij_tab)
 end if

 ABI_DT_MALLOC(Sigp%Sigcij_tab,(Sigp%nkptgw,Sigp%nsppol))
 ABI_DT_MALLOC(Sigp%Sigxij_tab,(Sigp%nkptgw,Sigp%nsppol))

 use_sym_at=.FALSE.
 if (PRESENT(Bnd_sym)) then
   do spin=1,Sigp%nsppol
     do ikcalc=1,Sigp%nkptgw
      ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
      use_sym_at(ik_ibz,spin) = ( .not.esymm_failed(Bnd_sym(ik_ibz,spin)) )
     end do
   end do
 end if

 do spin=1,Sigp%nsppol
   do ikcalc=1,Sigp%nkptgw
     ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))

     if (use_sym_at(ik_ibz,spin)) then
       if (gwcalctyp<20) then
         MSG_ERROR("You should not be here!")
       end if

       bmin=Sigp%minbnd(ikcalc,spin); bmax=Sigp%maxbnd(ikcalc,spin)
       ABI_DT_MALLOC(Sigp%Sigxij_tab(ikcalc,spin)%col,(bmin:bmax))
       ABI_DT_MALLOC(Sigp%Sigcij_tab(ikcalc,spin)%col,(bmin:bmax))

       do bcol=bmin,bmax
         ABI_MALLOC(sigc_bidx,(bmax-bmin+1))
         ABI_MALLOC(sigx_bidx,(bmax-bmin+1))

         if (Bnd_sym(ik_ibz,spin)%err_status/=0) then   ! Band classification failed.
           sigc_bidx = (/(ii,ii=bmin,bmax)/)
           idx_c = bmax-bmin+1
           sigx_bidx = (/(ii,ii=bmin,bcol)/) ! Hermitian
           idx_x = bcol-bmin+1
         else
           irr_idx2 = Bnd_sym(ik_ibz,spin)%b2irrep(bcol)
           idx_c = 0
           do brow=bmin,bmax
             irr_idx1 = Bnd_sym(ik_ibz,spin)%b2irrep(brow)
             if (sigma_is_herm(Sigp).and.bcol<brow) CYCLE  ! Only the upper triangle for HF, SEX, or COHSEX.
             if (irr_idx1 == irr_idx2) then ! same character, add this row to the list.
               idx_c = idx_c +1
               sigc_bidx(idx_c) = brow
             end if
           end do
           idx_x = 0
           do brow=bmin,bcol
             irr_idx1 = Bnd_sym(ik_ibz,spin)%b2irrep(brow)
             if (bcol<brow) CYCLE  ! Sig_x is always Hermitian.
             if (irr_idx1 == irr_idx2) then ! same character, add this row to the list.
               idx_x = idx_x +1
               sigx_bidx(idx_x) = brow
             end if
           end do
         end if
         !
         ! Table for Sigma_x matrix elements taking into account symmetries of the bands.
         ABI_MALLOC(Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx,(idx_x))

         Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%size1= idx_x
         Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(:) = sigx_bidx(1:idx_x)
         !write(std_out,*)" Sigxij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol
         !write(std_out,*)" size: ",idx_x,(Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(ii),ii=1,idx_x)
         !
         ! Table for Sigma_c matrix elements taking into account symmetries of the bands.
         ABI_MALLOC(Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx,(idx_c))

         Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%size1= idx_c
         Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(:) = sigc_bidx(1:idx_c)
         !write(std_out,*)" Sigcij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol
         !write(std_out,*)" size: ",idx_c,(Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(ii), ii=1,idx_c)

         ABI_FREE(sigx_bidx)
         ABI_FREE(sigc_bidx)
       end do ! bcol

     else  ! Symmetries cannot be used for this (k,s).

       bmin=Sigp%minbnd(ikcalc,spin); bmax=Sigp%maxbnd(ikcalc,spin)
       ABI_DT_MALLOC(Sigp%Sigcij_tab(ikcalc,spin)%col,(bmin:bmax))
       ABI_DT_MALLOC(Sigp%Sigxij_tab(ikcalc,spin)%col,(bmin:bmax))

       if (gwcalctyp<20) then  ! QP wavefunctions == KS, therefore only diagonal elements are calculated.
         do bcol=bmin,bmax
           ABI_MALLOC(Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx,(1:1))
           Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%size1= 1
           Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(1) = bcol
           ABI_MALLOC(Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx,(1:1))
           Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%size1= 1
           Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(1) = bcol
         end do
       else
         ! Use QP wavefunctions, Sigma_ij matrix is sparse but we have to classify the states in sigma.
         ! The only thing we can do here is filling the entire matrix taking advantage of Hermiticity (if any).
         do bcol=bmin,bmax
           ABI_MALLOC(Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx,(bcol-bmin+1))
           Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%size1= bcol-bmin+1
           Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(:) = (/(ii,ii=bmin,bcol)/) ! Sigma_x is Hermitian.
           !write(std_out,*)"Sigxij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol,Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(:)

           ABI_MALLOC(sigc_bidx,(bmax-bmin+1))
           idx_c = 0
           do brow=bmin,bmax
             if (sigma_is_herm(Sigp).and.bcol<brow) CYCLE  ! Only the upper triangle of Sigc_ij is needed (SEX, COHSEX).
             idx_c = idx_c +1
             sigc_bidx(idx_c) = brow
           end do
           ABI_MALLOC(Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx,(idx_c))
           Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%size1= idx_c
           Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(:) = sigc_bidx(1:idx_c)
           ABI_FREE(sigc_bidx)
           !write(std_out,*)"Sigcij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol,Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(:)
         end do
       end if
     end if

   end do !ikcalc
 end do !spin

end subroutine sigma_tables
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/sigma_bksmask
!! NAME
!! sigma_bksmask
!!
!! FUNCTION
!!  Compute tables for the distribution and the storage of the wavefunctions in the SIGMA code.
!!
!! INPUTS
!! Dtset<type(dataset_type)>=all input variables for this dataset
!! Sigp<sigparams_t>=Parameters governing the self-energy calculation.
!! Kmesh <kmesh_t>=Structure describing the k-point sampling.
!! nprocs=Total number of MPI processors
!! my_rank=Rank of this this processor.
!!
!! OUTPUT
!! my_spins(:)=Pointer to NULL in input. In output: list of spins treated by this node.
!! bks_mask(Sigp%nbnds,Kmesh%nibz,Sigp%nsppol)=True if this node will treat this state.
!! keep_ur(Sigp%nbnds,Kmesh%nibz,Sigp%nsppol)=True if this node will store this state in real space.
!! ierr=Exit status.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigma_bksmask(Dtset,Sigp,Kmesh,my_rank,nprocs,my_spins,bks_mask,keep_ur,ierr)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_gwdefs
 use m_errors
 use m_profiling_abi

 use m_bz_mesh,       only : kmesh_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigma_bksmask'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_rank,nprocs
 integer,intent(out) :: ierr
 type(Dataset_type),intent(in) :: Dtset
 type(sigparams_t),intent(in) :: Sigp
 type(kmesh_t),intent(in) :: Kmesh
!arrays
 integer,allocatable,intent(out) :: my_spins(:)
 logical,intent(out) :: bks_mask(Sigp%nbnds,Kmesh%nibz,Sigp%nsppol)
 logical,intent(out) :: keep_ur(Sigp%nbnds,Kmesh%nibz,Sigp%nsppol)

!Local variables-------------------------------
!scalars
 integer :: my_nspins,my_maxb,my_minb,isp,spin,band,nsppol,rank_spin
 logical :: store_ur
!arrays
 integer :: tmp_spins(Sigp%nsppol),nprocs_spin(Sigp%nsppol)

! *************************************************************************

 ierr=0; nsppol=Sigp%nsppol

 ! List of spins for each node, number of processors per each spin
 ! and the MPI rank in the "spin" communicator.
 my_nspins=nsppol
 ABI_MALLOC(my_spins, (nsppol))
 my_spins= [(isp, isp=1,nsppol)]
 nprocs_spin = nprocs; rank_spin = my_rank

 if (nsppol==2 .and. nprocs>1) then
   ! Distribute spins (optimal distribution if nprocs is even)
   nprocs_spin(1) = nprocs/2
   nprocs_spin(2) = nprocs - nprocs/2
   my_nspins=1
   my_spins(1)=1
   if (my_rank+1>nprocs/2) then
     ! I will treate spin=2, compute shifted rank.
     my_spins(1)=2
     rank_spin = my_rank - nprocs/2
   end if
 end if

 store_ur = (MODULO(Dtset%gwmem,10)==1)
 keep_ur=.FALSE.; bks_mask=.FALSE.

 select case (Dtset%gwpara)
 case (1)
   ! Parallelization over transitions **without** memory distributions (Except for the spin).
   my_minb=1; my_maxb=Sigp%nbnds
   do isp=1,my_nspins
     spin = my_spins(isp)
     bks_mask(my_minb:my_maxb,:,spin)=.TRUE.
     if (store_ur) keep_ur(my_minb:my_maxb,:,spin)=.TRUE.
   end do

 case (2)
   ! Distribute bands and spin (alternating planes of bands)
   do isp=1,my_nspins
     spin = my_spins(isp)
     do band=1,Sigp%nbnds
       if (xmpi_distrib_with_replicas(band,Sigp%nbnds,rank_spin,nprocs_spin(spin))) then
       !if (MODULO(band-1,nprocs_spin(spin))==rank_spin) then
         bks_mask(band,:,spin)=.TRUE.
         if (store_ur) keep_ur(band,:,spin)=.TRUE.
       end if
     end do
   end do

#if 0
   ! Each node store the full set of occupied states to speed up Sigma_x.
   do isp=1,my_nspins
     spin = my_spins(isp)
     do ik_ibz=1,Kmesh%nibz
       ks_iv=ks_vbik(ik_ibz,spin) ! Valence index for this (k,s)
       bks_mask(1:ks_iv,:,spin)=.TRUE.
       if (store_ur) keep_ur(1:ks_iv,:,spin)=.TRUE.
     end do
   end do
#endif

 case default
   MSG_WARNING("Wrong value for gwpara")
   ierr = 1
 end select

 ! Return my_spins with correct size.
 tmp_spins(1:my_nspins) = my_spins(1:my_nspins)

 ABI_FREE(my_spins)
 ABI_MALLOC(my_spins, (my_nspins))
 my_spins = tmp_spins(1:my_nspins)

end subroutine sigma_bksmask
!!***
