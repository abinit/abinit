!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_screening
!! NAME
!! setup_screening
!!
!! FUNCTION
!!  Initialize the Ep% data type containing the parameters used during the screening calculation.
!!  as well as basic objects describing the BZ sampling .... TODO list to be completed
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2016 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! wfk_fname=Name of the input WFK file.
!! acell(3)=length scales of primitive translations (Bohr).
!! rprim(3,3)=dimensionless real space primitive translations.
!! ngfftf(18)=Contain all needed information about the 3D FFT for densities and potentials.
!! dtfil <type(datafiles_type)>=variables related to files
!!
!! OUTPUT
!! ngfft_gw(18)=Contain all needed information about the 3D FFT for the oscillator strengths.
!!  See ~abinit/doc/input_variables/vargs.htm#ngfft
!! Ltg_q(:)<littlegroup_t>,=
!! Ep<em1params_t>=Parameters for the screening calculation.
!!  Most part of it is Initialized and checked.
!! Hdr_wfk type(Hdr_type)=Header of the KSS file.
!! Cryst<crystal_t>=Definition of the unit cell and its symmetries.
!! Kmesh<kmesh_t>=Structure defining the k-point sampling (wavefunctions).
!! Qmesh<kmesh_t>=Structure defining the q-point sampling (screening)
!! Gsph_wfn<gsphere_t>=Structure defining the G-sphere for the wavefunctions (not k-dependent).
!! Gsph_epsG0<gsphere_t>=The G-sphere for the screening, enlarged to take into account for umklapps.
!! Psps <Pseudopotential_type)>=Info on pseudopotential, only for consistency check of the KSS file
!! Vcp <type vcoul_t> datatype gathering information on the coulombian cutoff technique
!! comm=MPI communicator.
!!
!! SIDE EFFECTS
!! Dtset<Dataset_type>=All input variables for this dataset.
!!  %ecutwfn, %npwwfn, 
!!  %ecuteps, %npweps
!!   might be redefined in setshells in order to close the shell.
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      xmpi_split_work2_i4b
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setup_screening(codvsn,acell,rprim,ngfftf,wfk_fname,dtfil,Dtset,Psps,Pawtab,&
& ngfft_gw,Hdr_wfk,Hdr_out,Cryst,Kmesh,Qmesh,KS_BSt,Ltg_q,Gsph_epsG0,Gsph_wfn,Vcp,Ep,comm)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr

 use m_gwdefs,        only : GW_TOLQ0, GW_TOLQ, GW_Q0_DEFAULT, czero_gw, em1params_t
 use m_fstrings,      only : strcat, sjoin
 use m_io_tools,      only : file_exists
 use m_geometry,      only : normv
 use m_crystal,       only : crystal_print, crystal_t
 use m_crystal_io,    only : crystal_from_hdr
 use m_bz_mesh,       only : kmesh_t, kmesh_init, get_ng0sh, kmesh_print, find_qmesh, get_BZ_item,&
&                            littlegroup_t, littlegroup_init, make_mesh, kmesh_free
 use m_ebands,        only : ebands_init
 use m_vcoul,         only : vcoul_t, vcoul_init
 use m_fft_mesh,      only : setmesh
 use m_gsphere,       only : gsphere_t, gsph_init, merge_and_sort_kg, gsph_free, setshells
 use m_pawtab,        only : pawtab_type
 use m_pawrhoij,      only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free
 use m_io_kss,        only : testkss, make_gvec_kss
 use m_wfk,           only : wfk_read_eigenvalues

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_screening'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_56_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 character(len=6),intent(in) :: codvsn
 character(len=fnlen),intent(in) :: wfk_fname
 type(Dataset_type),intent(inout) :: Dtset !INOUT is due to setshells
 type(datafiles_type),intent(in) :: dtfil
 type(Pseudopotential_type),intent(in) :: Psps
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Dtset%usepaw)
 type(em1params_t),intent(out) :: Ep
 type(Hdr_type),intent(out) :: Hdr_wfk,Hdr_out
 type(ebands_t),intent(out) :: KS_BSt
 type(kmesh_t),intent(out) :: Kmesh,Qmesh
 type(crystal_t),intent(out) :: Cryst
 type(gsphere_t),intent(out) :: Gsph_epsG0,Gsph_wfn
 type(vcoul_t),intent(out) :: Vcp
!arrays
 integer,intent(in) :: ngfftf(18)
 integer,intent(out) :: ngfft_gw(18)
 real(dp),intent(in) :: acell(3),rprim(3,3)
 type(littlegroup_t),pointer :: Ltg_q(:)

!Local variables-------------------------------
!scalars
 integer,parameter :: NOMEGAGAUSS=30,NOMEGAREAL=201
 integer,parameter :: pertcase0=0,master=0
 integer :: bantot,ib,ibtot,ikibz,iq,iqp,isppol,ig,ng,ierr
 integer :: jj,mod10,mband,ng_kss,iqbz,isym,iq_ibz,itim
 integer :: timrev,use_umklp,ncerr
 integer :: npwepG0,nshepspG0,method,enforce_sym,nfftgw_tot !,spin,band,ik_ibz,
 integer :: istart,iend,test_npwkss,my_rank,nprocs !ii
 real(dp),parameter :: OMEGAERMAX=100.0/Ha_eV
 real(dp) :: ecutepspG0,ucvol,domegareal
 logical :: remove_inv,ltest,found,is_static,has_q0
 character(len=500) :: msg
 type(wvl_internal_type) :: wvl
!arrays
 integer :: ng0sh_opt(3)
 integer,allocatable :: npwarr(:)
 integer,pointer :: gvec_kss(:,:)
 integer,pointer :: test_gvec_kss(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),qtmp(3),sq(3),qbz(3)
 real(dp),pointer :: energies_p(:,:,:)
 real(dp),allocatable :: doccde(:),eigen(:),occfact(:)
 type(Pawrhoij_type),allocatable :: Pawrhoij(:)

! *************************************************************************

 DBG_ENTER('COLL')

 ! === Check for calculations that are not implemented ===
 ltest=ALL(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol)==Dtset%nband(1))
 ABI_CHECK(ltest,'Dtset%nband(:) must be constant')

 my_rank = xmpi_comm_rank(comm)
 nprocs  = xmpi_comm_size(comm)

 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 ! === Set up basic parameters of the calculation ===
 Ep%gwcalctyp            =Dtset%gwcalctyp
 Ep%plasmon_pole_model   =.TRUE.
 Ep%analytic_continuation=.FALSE.
 Ep%contour_deformation  =.FALSE.

 mod10=MOD(Ep%gwcalctyp,10)
 if (mod10/=0.and.mod10/=8) Ep%plasmon_pole_model   =.FALSE.
 if (mod10==1)              Ep%analytic_continuation=.TRUE.
 if (mod10==2.or.mod10==9)  Ep%contour_deformation  =.TRUE.
 is_static=(mod10==5.or.mod10==6.or.mod10==7)

 Ep%nbnds  =Dtset%nband(1)
 Ep%symchi =Dtset%symchi
 Ep%inclvkb=Dtset%inclvkb; if (Dtset%usepaw/=0) Ep%inclvkb=0
 Ep%zcut   =Dtset%zcut

 write(msg,'(2a,i4,2a,f10.6,a)')ch10,&
&  ' GW calculation type              = ',Ep%gwcalctyp,ch10,&
&  ' zcut to avoid poles in chi0 [eV] = ',Ep%zcut*Ha_eV,ch10
 call wrtout(std_out,msg,'COLL')
 !
 Ep%awtr  =Dtset%awtr
 Ep%npwe  =Dtset%npweps
 Ep%npwwfn=Dtset%npwwfn
 Ep%npwvec=MAX(Ep%npwe,Ep%npwwfn)

 timrev = 2 ! This information is not reported in the header
            ! 1 --> do not use time-reversal symmetry
            ! 2 --> take advantage of time-reversal symmetry
 if (timrev==1.and.Dtset%awtr/=0) then
   msg = "awtr/=0 cannot be used when time-reversal symmetry doesn't hold"
   MSG_ERROR(msg)
 end if
 !
 ! === Read parameters from (WFK|KSS) and verifify them ===
 call wfk_read_eigenvalues(wfk_fname,energies_p,Hdr_wfk,comm)
 mband = MAXVAL(Hdr_wfk%nband)

 call hdr_vs_dtset(Hdr_wfk,Dtset)
 remove_inv=.FALSE.

 test_npwkss = 0
 call make_gvec_kss(Dtset%nkpt,Dtset%kptns,Hdr_wfk%ecut_eff,Dtset%symmorphi,Dtset%nsym,Dtset%symrel,Dtset%tnons,&
&   gprimd,Dtset%prtvol,test_npwkss,test_gvec_kss,ierr)
 ABI_CHECK(ierr==0,"Fatal error in make_gvec_kss")

 ABI_MALLOC(gvec_kss,(3,test_npwkss))
 gvec_kss = test_gvec_kss
 ng_kss = test_npwkss

 if (Ep%npwvec>ng_kss) then
   Ep%npwvec=ng_kss
   if (Ep%npwwfn> ng_kss) Ep%npwwfn=ng_kss
   if (Ep%npwe  > ng_kss) Ep%npwe  =ng_kss
   write(msg,'(3a,3(a,i6,a))')ch10,&
&    ' Number of G-vectors found less then required. Calculation will proceed with ',ch10,&
&    '  npwvec = ',Ep%npwvec,ch10,&
&    '  npweps = ',Ep%npwe  ,ch10,&
&    '  npwwfn = ',Ep%npwwfn,ch10
   MSG_WARNING(msg)
 end if

 ng = MIN(SIZE(gvec_kss,DIM=2),SIZE(test_gvec_kss,DIM=2))
 ierr = 0
 do ig=1,ng
   if (ANY(gvec_kss(:,ig)/=test_gvec_kss(:,ig))) then
     ierr=ierr+1
     write(std_out,*)" gvec_kss ",ig,"/",ng,gvec_kss(:,ig),test_gvec_kss(:,ig)
   end if
 end do
 ABI_CHECK(ierr==0,"Mismatch between gvec_kss and test_gvec_kss")
 ABI_FREE(test_gvec_kss)

 ! === Get important dimension from Hdr_wfk ===
 ! * Check also the consistency btw Hdr_wfk and Dtset.
 Ep%nsppol=Hdr_wfk%nsppol
 Ep%nkibz =Hdr_wfk%nkpt

 if (Ep%nbnds>mband) then
   Ep%nbnds=mband
   Dtset%nband(:)=mband
   write(msg,'(4a,i4,a)')ch10,&
&    ' Number of bands found less then required. ',ch10,&
&    ' Calculation will proceed with nbnds = ',mband,ch10
   MSG_WARNING(msg)
 end if

 call crystal_from_hdr(Cryst,Hdr_wfk,timrev,remove_inv)
 call crystal_print(Cryst,mode_paral='COLL')

 ! === Create basic data types for the calculation ===
 ! * Kmesh defines the k-point sampling for the wavefunctions.
 ! * Qmesh defines the q-point sampling for chi0, all possible differences k1-k2 reduced to the IBZ.
 ! TODO Kmesh%bz should be [-half,half[ but this modification will be painful!

 call kmesh_init(Kmesh,Cryst,Ep%nkibz,Hdr_wfk%kptns,Dtset%kptopt,wrap_1zone=.FALSE.)
 !call kmesh_init(Kmesh,Cryst,Ep%nkibz,Hdr_wfk%kptns,Dtset%kptopt,wrap_1zone=.TRUE.)
 ! Some required information are not filled up inside kmesh_init
 ! So doing it here, even though it is not clean
 Kmesh%kptrlatt(:,:) =Dtset%kptrlatt(:,:)
 Kmesh%nshift        =Dtset%nshiftk
 ABI_ALLOCATE(Kmesh%shift,(3,Kmesh%nshift))
 Kmesh%shift(:,:)    =Dtset%shiftk(:,1:Dtset%nshiftk)

 call kmesh_print(Kmesh,"K-mesh for the wavefunctions",std_out,Dtset%prtvol,"COLL")
 call kmesh_print(Kmesh,"K-mesh for the wavefunctions",ab_out, 0,           "COLL")
 !
 ! === Find Q-mesh, and do setup for long wavelength limit ===
 ! * Stop if a nonzero umklapp is needed to reconstruct the BZ. In this case, indeed,
 !   epsilon^-1(Sq) should be symmetrized in csigme using a different expression (G-G_o is needed)
 !
 call find_qmesh(Qmesh,Cryst,Kmesh)

 call kmesh_print(Qmesh,"Q-mesh for the screening function",std_out,Dtset%prtvol,"COLL")
 call kmesh_print(Qmesh,"Q-mesh for the screening function",ab_out ,0           ,"COLL")

 do iqbz=1,Qmesh%nbz
   call get_BZ_item(Qmesh,iqbz,qbz,iq_ibz,isym,itim)
   sq = (3-2*itim)*MATMUL(Cryst%symrec(:,:,isym),Qmesh%ibz(:,iq_ibz))
   if (ANY(ABS(qbz-sq )>1.0d-4)) then
     write(msg,'(a,3f6.3,a,3f6.3,2a,9i3,a,i2,2a)')&
&      ' qpoint ',qbz,' is the symmetric of ',Qmesh%ibz(:,iq_ibz),ch10,&
&      ' through operation ',Cryst%symrec(:,:,isym),' and itim ',itim,ch10,&
&      ' however a non zero umklapp G_o vector is required and this is not yet allowed'
     MSG_ERROR(msg)
   end if
 end do

 ! Write the list of qpoints for the screening in netcdf format and exit.
 ! This file is used by abipy to generate multiple input files.
 if (Dtset%nqptdm == -1) then
   if (my_rank==master) then
#ifdef HAVE_NETCDF
      ncerr = nctk_write_ibz(strcat(dtfil%filnam_ds(4), "_qptdms.nc"), qmesh%ibz, qmesh%wt)
      NCF_CHECK(ncerr)
#endif
   end if
   MSG_ERROR_NODUMP("Aborting now")
 end if

 if (Dtset%gw_nqlwl==0) then
   Ep%nqlwl=1
   ABI_MALLOC(Ep%qlwl,(3,Ep%nqlwl))
   Ep%qlwl(:,1)=GW_Q0_DEFAULT ! Use default shift 0.000010, 0.000020, 0.000030
 else
   Ep%nqlwl=Dtset%gw_nqlwl
   ABI_MALLOC(Ep%qlwl,(3,Ep%nqlwl))
   Ep%qlwl(:,:)=Dtset%gw_qlwl(:,1:Ep%nqlwl)
   ABI_CHECK(Ep%nqlwl==1,"nqlwl/=1 not coded")
 end if
 !write(std_out,*)" Using qlwl = ",Ep%qlwl

 ! === Find optimal value for G-sphere enlargment due to oscillator matrix elements ===
 call get_ng0sh(Kmesh%nbz,Kmesh%bz,Qmesh%nibz,Qmesh%ibz,Kmesh%nbz,Kmesh%bz,GW_TOLQ0,ng0sh_opt)

 write(msg,'(a,3i2)') ' Optimal value for ng0sh = ',ng0sh_opt(:)
 call wrtout(std_out,msg,"COLL")

 Ep%mG0(:)=ng0sh_opt(:)
 !Ep%mG0(:)=(/3,3,3/)

 ! === In case of symmetrization, find the little group of the q"s ===
 ! * For the long-wavelength limit we consider a small but finite q. However the oscillators are
 !  evaluated setting q==0. Thus it is possible to take advantage of symmetries also when q --> 0.
 ! * Here we calculate the enlargement of the G-sphere, npwepG0, needed to account for umklapps.
 ! TODO Switch on use_umklp, write all this stuff to ab_out

 Ep%npwepG0=Ep%npwe
 ABI_DT_MALLOC(Ltg_q,(Qmesh%nibz))

 do iq=1,Qmesh%nibz
   qtmp(:)=Qmesh%ibz(:,iq); if (normv(qtmp,gmet,'G')<GW_TOLQ0) qtmp(:)=zero; use_umklp=0
   call littlegroup_init(qtmp,Kmesh,Cryst,use_umklp,Ltg_q(iq),Ep%npwe,gvec=gvec_kss)
 end do

 ecutepspG0 = Dtset%ecuteps
 if (Ep%symchi/=0) then
   ecutepspG0=MAXVAL(Ltg_q(:)%max_kin_gmG0)+tol6; npwepG0=0; nshepspG0=0
   write(std_out,*)" Due to umklapp processes : ecutepspg0= ",ecutepspG0
   call setshells(ecutepspG0,npwepG0,nshepspG0,Cryst%nsym,gmet,gprimd,Cryst%symrel,'eps_pG0',Cryst%ucvol)
   Ep%npwepG0=npwepG0
 end if

 if (Ep%npwepG0>Ep%npwvec) then
   write(msg,'(3a,i5,a,i5)')&
&    ' npwepG0 > npwvec, decrease npweps or increase npwwfn. ',ch10,&
&    ' npwepG0 = ',Ep%npwepG0,' npwvec = ',Ep%npwvec
   MSG_ERROR(msg)
 end if

 ! === Create structure describing the G-sphere used for chi0/espilon and Wfns ===
 ! * The cutoff is >= ecuteps to allow for umklapp
#if 0
 call gsph_init(Gsph_wfn,Cryst,0,ecut=Dtset%ecutwfn)
 Dtset%npwwfn = Gsph_wfn%ng
 Ep%npwwfn = Gsph_wfn%ng
 ierr = 0
 do ig=1,MIN(Gsph_wfn%ng, ng_kss)
   if ( ANY(Gsph_wfn%gvec(:,ig) /= gvec_kss(:,ig)) ) then
     write(std_out,*)ig, Gsph_wfn%gvec(:,ig), gvec_kss(:,ig)
   end if
 end do
 ABI_CHECK(ierr==0,"Wrong gvec_wfn")
#else
 call gsph_init(Gsph_wfn,Cryst,Ep%npwvec,gvec=gvec_kss)
#endif

#if 0
 call gsph_init(Gsph_epsG0,Cryst,0,ecut=ecutepspG0)
 Ep%npwepG0 = Gsph_epsG0%ng
 ierr = 0
 do ig=1,MIN(Gsph_epsG0%ng, ng_kss)
   if ( ANY(Gsph_epsG0%gvec(:,ig) /= gvec_kss(:,ig)) ) then
     write(std_out,*)ig, Gsph_epsG0%gvec(:,ig), gvec_kss(:,ig)
   end if
 end do
 ABI_CHECK(ierr==0,"Wrong gvec_epsG0")
#else
 call gsph_init(Gsph_epsG0,Cryst,Ep%npwepG0,gvec=gvec_kss)
#endif
 !
 ! =======================================================================
 ! ==== Setup of the FFT mesh used for the oscillator matrix elements ====
 ! =======================================================================
 ! * ngfft_gw(7:18) is the same as Dtset%ngfft(7:18), initialized before entering setup_screening.
 !   Here we just redefine ngfft_gw(1:6) according to the following options:
 !
 !    method==0 ==> FFT grid read from __fft.in__ (only for debugging purpose)
 !    method==1 ==> normal FFT grid
 !    method==2 ==> slightly augmented FFT grid to calculate exactly rho_tw_g (see setmesh.F90)
 !    method==3 ==> doubled FFT grid, to treat exactly the convolution defining the density,
 !      Useful in sigma if ppmodel=[2,3,4] since rho(G-Gp) or to calculate matrix elements of v_Hxc.
 !
 !    enforce_sym==1 ==> enforce a direct space FFT mesh compatible with all symmetries operation
 !    enforce_sym==0 ==> Find the smallest FFT grid compatibile with the library, do not care about symmetries
 !
 ngfft_gw(1:18)=Dtset%ngfft(1:18); method=2
 if (Dtset%fftgw==00 .or. Dtset%fftgw==01) method=0
 if (Dtset%fftgw==10 .or. Dtset%fftgw==11) method=1
 if (Dtset%fftgw==20 .or. Dtset%fftgw==21) method=2
 if (Dtset%fftgw==30 .or. Dtset%fftgw==31) method=3
 enforce_sym=MOD(Dtset%fftgw,10)

 ! npwepG0 to account for umklapps.
 call setmesh(gmet,gvec_kss,ngfft_gw,Ep%npwvec,Ep%npwepG0,Ep%npwwfn,nfftgw_tot,method,Ep%mG0,Cryst,enforce_sym)
 !call new_setmesh(Cryst,ecut_osc,ecutwfn,nkpt,kpoints,method,Ep%mG0,enforce_sym,ngfft_gw,nfftgw_tot)

 ABI_FREE(gvec_kss)

 ! FIXME this wont work if nqptdm/=0
 call vcoul_init(Vcp,Gsph_epsG0,Cryst,Qmesh,Kmesh,Dtset%rcut,Dtset%icutcoul,Dtset%vcutgeo,Dtset%ecuteps,Ep%npwe,Ep%nqlwl,&
&  Ep%qlwl,ngfftf,comm)

#if 0
 ! Using the random q for the optical limit is one of the reasons
 ! why sigma breaks the initial energy degeneracies.
 Vcp%i_sz=zero
 Vcp%vc_sqrt(1,:)=czero
 Vcp%vcqlwl_sqrt(1,:)=czero
#endif

 ! Value of scissor energy
 Ep%mbpt_sciss=Dtset%mbpt_sciss

 ! === Define the frequency mesh for epsilon according to the method used ===
 Ep%nomegaei=1
 Ep%nomegaer=1; if (is_static) Ep%nomegaer=0
 Ep%nomegaec=0
 Ep%omegaermax=zero

 ! === For ppmodels 2,3,4, only omega=0 is needed ===
 if (Ep%plasmon_pole_model.and.Dtset%nfreqre==1.and.Dtset%nfreqim==0) then
   Ep%nomegaer=1; Ep%nomegaei=0
   write(msg,'(7a)')ch10,&
&    ' The inverse dielectric matrix will be calculated on zero frequency only',ch10,&
&    ' please note that the calculated epsilon^-1 cannot be used ',ch10,&
&    ' to calculate QP corrections using plasmonpole model 1',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if
 !
 ! === Max number of omega along the imaginary axis ===
 if (Ep%analytic_continuation.or.Ep%contour_deformation) then
   Ep%nomegaei=Dtset%nfreqim
   if (Dtset%gw_frqim_inzgrid==1) then
     write(msg,'(a)')' iomega = z/1-z transfom grid will be used for imaginary frequency grid'
     MSG_WARNING(msg)
   end if
   if (Dtset%cd_customnimfrqs/=0) then
     write(msg,'(a)')' Custom imaginary grid specified. Assuming experienced user.'
     MSG_WARNING(msg)
     Ep%nomegaei=Dtset%cd_customnimfrqs
   end if
   if (Ep%nomegaei==-1) then
     Ep%nomegaei=NOMEGAGAUSS
     write(msg,'(a,i5)')' Number of imaginary frequencies set to default= ',NOMEGAGAUSS
     MSG_WARNING(msg)
   end if
   if (Ep%nomegaei==0) then
     write(msg,'(a)')' nfreqim = 0 ! Assuming experienced user merging several frequency calculations.'
     MSG_WARNING(msg)
   end if
 end if

 ! === Range and total number of real frequencies ===
 Ep%omegaermin = zero
 if (Ep%contour_deformation) then
   Ep%nomegaer=Dtset%nfreqre; Ep%omegaermin=Dtset%freqremin; Ep%omegaermax=Dtset%freqremax
   if (Dtset%gw_frqre_tangrid==1) then
     Ep%omegaermax=Dtset%cd_max_freq
     write(msg,'(a)')' Tangent transfom grid will be used for real frequency grid'
     MSG_WARNING(msg)
   end if
   if (Dtset%gw_frqre_tangrid==1) then
     write(msg,'(a)')' Tangent transfom grid will be used for real frequency grid'
     MSG_WARNING(msg)
   end if
   if (Ep%nomegaer==-1) then 
     Ep%nomegaer=NOMEGAREAL
     write(msg,'(a,i5)')' Number of real frequencies set to default= ',NOMEGAREAL
     MSG_WARNING(msg)
   end if
   if (Ep%nomegaer==0) then
     write(msg,'(a)')'  nfreqre = 0 ! Assuming experienced user.'
     MSG_WARNING(msg)
   end if
   if (ABS(Ep%omegaermin)<TOL16) then
     Ep%omegaermin=zero
     write(msg,'(a,f8.4)')' Min real frequency set to default [Ha] = ',Ep%omegaermin
     MSG_WARNING(msg)
   end if
   if (Ep%omegaermin>Ep%omegaermax) then
     MSG_ERROR(' freqremin > freqremax !')
   end if
   if (Ep%omegaermax<TOL16) then
     Ep%omegaermax=OMEGAERMAX
     write(msg,'(a,f8.4)')' Max real frequency set to default [Ha] = ',OMEGAERMAX
     MSG_WARNING(msg)
   end if
   ! Check if a subset of the frequencies is to be used
   if (Dtset%cd_subset_freq(1)/=0) then
     istart = Dtset%cd_subset_freq(1)
     iend   = Dtset%cd_subset_freq(2)
     if (istart>iend.or.istart<0.or.iend<0) then
       MSG_ERROR(' check indices of cd_subset_freq!')
     end if
     write(msg,'(2(a,i0))')' Using cd_subset_freq to only do freq. from ',istart,' to ',iend
     MSG_WARNING(msg)
     ! Reset the numbers
     if (Dtset%gw_frqre_tangrid/=1) then ! Normal equidistant grid 
       Ep%nomegaer = iend-istart+1
       domegareal=(Ep%omegaermax-Ep%omegaermin)/(Ep%nomegaer-1)
       Ep%omegaermin = Ep%omegaermin+(istart-1)*domegareal
       Ep%omegaermax = Ep%omegaermin+(iend-1)*domegareal
     else
       Ep%nomegaer = iend-istart+1
     end if
   end if
 end if

 ! Check full grid calculations
 if (Dtset%cd_full_grid/=0) then
   MSG_WARNING(' FULL GRID IN COMPLEX PLANE CALCULATED.')
   MSG_WARNING(' YOU MIGHT NOT BE ABLE TO USE SCREENING FILES!')
   if (Dtset%cd_subset_freq(1)/=0) then
     MSG_ERROR(' cd_subset_freq cannot be used with cd_full_grid!')
   end if
   Ep%nomegaec = Ep%nomegaei*(Ep%nomegaer-1)
 end if

 Ep%nomega=Ep%nomegaer+Ep%nomegaei+Ep%nomegaec ! Total number of frequencies.

 ! ==== Setup of the spectral method ====
 Ep%spmeth  =Dtset%spmeth
 Ep%nomegasf=Dtset%nomegasf
 Ep%spsmear =Dtset%spbroad

 if (Ep%spmeth/=0) then
   write(msg,'(2a,i3,2a,i8)')ch10,&
&    ' setup_screening : using spectral method = ',Ep%spmeth,ch10,&
&    '  Number of frequencies for imaginary part = ',Ep%nomegasf
   call wrtout(std_out,msg,'COLL')
   if (Ep%spmeth==2) then
     write(msg,'(a,f8.5,a)')' gaussian broadening = ',Ep%spsmear*Ha_eV,' [eV]'
     call wrtout(std_out,msg,'COLL')
   end if
 end if

 Ep%nI=1; Ep%nJ=1
 if (Dtset%nspinor==2) then
   if (Dtset%usepaw==1.and.Dtset%pawspnorb>0) then
     Ep%nI=1; Ep%nJ=4
   end if
   ! For spin-spin interaction
   ! Ep%nI=4; Ep%nJ=4
   ABI_CHECK(Ep%npwepG0==Ep%npwe,"npwepG0 must be == npwe if spinor==2")
   ABI_CHECK(Ep%symchi==0,"symchi/=0 and nspinor=2 not available")
 end if
 !
 ! === Enable the calculations of chi0 on user-specified q-points ===
 Ep%nqibz=Qmesh%nibz
 ABI_MALLOC(Ep%qibz,(3,Ep%nqibz))
 Ep%qibz(:,:)=Qmesh%ibz(:,:)

 Ep%nqcalc=Ep%nqibz
 if (Dtset%nqptdm>0) Ep%nqcalc=Dtset%nqptdm

 ABI_MALLOC(Ep%qcalc,(3,Ep%nqcalc))
 if (Ep%nqcalc/=Ep%nqibz) then
   write(msg,'(6a)')ch10,&
&    ' Dielectric matrix will be calculated only for some ',ch10,&
&    ' selected q points provided by the user through the input variables ',ch10,&
&    ' nqptdm and qptdm'
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   ltest=(Ep%nqcalc<=Qmesh%nibz)
   msg = 'nqptdm should not exceed the number of q points in the IBZ'
   ABI_CHECK(ltest,msg)
   Ep%qcalc(:,:)=Dtset%qptdm(:,1:Ep%nqcalc)
   do iq=1,Ep%nqcalc ! * Check whether the q-points provided are correct.
     found=.FALSE.
     do iqp=1,Qmesh%nibz
       qtmp(:)=Ep%qcalc(:,iq)-Qmesh%ibz(:,iqp)
       found=(normv(qtmp,gmet,'G')<GW_TOLQ)
       if (found) EXIT
     end do
     msg = 'One or more points specified by Dtset%qptdm do not satisfy q=k1-k2'
     ABI_CHECK(found,msg)
   end do
 else
   Ep%qcalc(:,:)=Ep%qibz(:,:)
 end if

 ! To write the SCR header correctly, with heads and wings, we have
 ! to make sure that q==0, if present, is the first q-point in the list.
 !has_q0=(ANY(normv(Ep%qcalc(:,:),gmet,'G')<GW_TOLQ0)) !commented to avoid problems with sunstudio12
 has_q0=.FALSE.
 do iq=1,Ep%nqcalc
   if (normv(Ep%qcalc(:,iq),gmet,'G')<GW_TOLQ0) then
     has_q0=.TRUE.; EXIT
   end if
 end do

 if (has_q0.and.normv(Ep%qcalc(:,1),gmet,'G')>=GW_TOLQ0) then
   write(msg,'(5a)')&
&    'The list of q-points to be calculated contains the Gamma point, ',ch10,&
&    'however Gamma is not the first point in the list. ' ,ch10,&
&    'Please, change your input file accordingly. '
   MSG_ERROR(msg)
 end if

 ! === Initialize the band structure datatype ===
 ! * Copy KSS energies and occupations up to Ep%nbnds==Dtset%nband(:)
 ! TODO Recheck symmorphy and inversion

 bantot=SUM(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol))

 ABI_MALLOC(doccde,(bantot))
 ABI_MALLOC(eigen,(bantot))
 ABI_MALLOC(occfact,(bantot))
 doccde(:)=zero; eigen(:)=zero; occfact(:)=zero

 jj=0; ibtot=0
 do isppol=1,Dtset%nsppol
   do ikibz=1,Dtset%nkpt
     do ib=1,Hdr_wfk%nband(ikibz+Dtset%nkpt*(isppol-1))
       ibtot=ibtot+1
       if (ib<=Ep%nbnds) then
         jj=jj+1
         occfact(jj)=Hdr_wfk%occ(ibtot)
         eigen  (jj)=energies_p(ib,ikibz,isppol)
       end if
     end do
   end do
 end do
 ABI_FREE(energies_p)

 ! Make sure that Dtset%wtk==Kmesh%wt due to the dirty treatment of
 ! the symmetry operations in the old GW code (symmorphy and inversion)
 ltest=(ALL(ABS(Dtset%wtk(1:Kmesh%nibz)-Kmesh%wt(1:Kmesh%nibz))<tol6))
 ABI_CHECK(ltest,'Mismatch between Dtset%wtk and Kmesh%wt')

 ABI_MALLOC(npwarr,(Hdr_wfk%nkpt))
 npwarr(:)=Ep%npwwfn

 call ebands_init(bantot,KS_BSt,Dtset%nelect,doccde,eigen,Dtset%istwfk,Kmesh%ibz,Dtset%nband,&
& Kmesh%nibz,npwarr,Dtset%nsppol,Dtset%nspinor,Dtset%tphysel,Dtset%tsmear,Dtset%occopt,occfact,Kmesh%wt,&
& dtset%charge, dtset%kptopt, dtset%kptrlatt_orig, dtset%nshiftk_orig, dtset%shiftk_orig, &
& dtset%kptrlatt, dtset%nshiftk, dtset%shiftk)

 ! TODO modify outkss in order to calculate the eigenvalues also if NSCF calculation.
 ! this fails simply because in case of NSCF occ  are zero
 !ltest=(ALL(ABS(occfact-KS_BSt%occ)<1.d-2))
 !call assert(ltest,'difference in occfact')
 !write(std_out,*)MAXVAL(ABS(occfact(:)-KS_BSt%occ(:)))

 !TODO call ebands_update_occ here
 !$call ebands_update_occ(KS_BSt,spinmagntarget,stmbias,Dtset%prtvol)

 ABI_FREE(doccde)
 ABI_FREE(eigen)
 ABI_FREE(npwarr)

 ! === Initialize abinit header for the screening part ===
 call hdr_init(KS_BSt,codvsn,Dtset,Hdr_out,Pawtab,pertcase0,Psps,wvl)

 ! === Get Pawrhoij from the header of the KSS file ===
 ABI_DT_MALLOC(Pawrhoij,(Cryst%natom*Dtset%usepaw))
 if (Dtset%usepaw==1) then
   call pawrhoij_alloc(Pawrhoij,1,Dtset%nspden,Dtset%nspinor,Dtset%nsppol,Cryst%typat,pawtab=Pawtab)
   call pawrhoij_copy(Hdr_wfk%Pawrhoij,Pawrhoij)
 end if

 call hdr_update(Hdr_out,bantot,1.0d20,1.0d20,1.0d20,Cryst%rprimd,occfact,Pawrhoij,Cryst%xred,dtset%amu_orig(:,1))

 ABI_FREE(occfact)
 call pawrhoij_free(Pawrhoij)
 ABI_DT_FREE(Pawrhoij)
 !
 ! ==== Setup of extrapolar technique ====
 Ep%gwcomp   = Dtset%gwcomp
 Ep%gwencomp = Dtset%gwencomp

 if (Ep%gwcomp==1) then
   write(msg,'(a,f8.2,a)')' Using the completeness correction with gwencomp ',Ep%gwencomp*Ha_eV,' [eV] '
   call wrtout(std_out,msg,'COLL')
 end if

 ! === Final compatibility tests ===
 if (ANY(KS_BSt%istwfk/=1)) then
   MSG_WARNING('istwfk/=1 is still under development')
 end if

 ltest=(KS_BSt%mband==Ep%nbnds.and.ALL(KS_BSt%nband==Ep%nbnds))
 ABI_CHECK(ltest,'BUG in definition of KS_BSt%nband')

 if (Ep%gwcomp==1.and.Ep%spmeth>0) then
   MSG_ERROR("Hilbert transform and extrapolar method are not compatible")
 end if

 DBG_EXIT('COLL')

end subroutine setup_screening
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/chi0_bksmask
!! NAME
!! chi0_bksmask
!!
!! FUNCTION
!!  Compute tables for the distribution and the storage of the wavefunctions in the SCREENING code.
!!
!! INPUTS
!! Dtset<type(dataset_type)>=all input variables for this dataset
!! Ep<em1params_t>=Parameters for the screening calculation.
!! Kmesh <kmesh_t>=Structure describing the k-point sampling.
!! nbvw = Max. number of fully/partially occupied states over spin
!! nbcw = Max. number of unoccupied states considering the spin
!! nprocs=Total number of MPI processors
!! my_rank=Rank of this this processor.
!!
!! OUTPUT
!! bks_mask(Ep%nbnds,Kmesh%nibz,Sigp%nsppol)=True if this node will treat this state.
!! keep_ur(Ep%nbnds,Kmesh%nibz,Sigp%nsppol)=True if this node will store this state in real space.
!! ierr=Exit status.
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      xmpi_split_work2_i4b
!!
!! SOURCE

subroutine chi0_bksmask(Dtset,Ep,Kmesh,nbvw,nbcw,my_rank,nprocs,bks_mask,keep_ur,ierr)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_gwdefs
 use m_errors
 use m_xmpi

 use m_bz_mesh,       only : kmesh_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chi0_bksmask'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_rank,nprocs,nbvw,nbcw
 integer,intent(out) :: ierr
 type(Dataset_type),intent(in) :: Dtset
 type(em1params_t),intent(in) :: Ep
 type(kmesh_t),intent(in) :: Kmesh
!arrays
 logical,intent(out) :: bks_mask(Ep%nbnds,Kmesh%nibz,Dtset%nsppol)
 logical,intent(out) :: keep_ur(Ep%nbnds,Kmesh%nibz,Dtset%nsppol)

!Local variables-------------------------------
!scalars
 integer :: my_nspins,my_maxb,my_minb,isp,spin,distr_err,nsppol,band,rank_spin,ib
 character(len=500) :: msg
 logical :: store_ur
!arrays
 integer :: my_spins(Dtset%nsppol),nprocs_spin(Dtset%nsppol) 
 integer,allocatable :: istart(:),istop(:)

! *************************************************************************

 ierr=0; nsppol=Dtset%nsppol

 my_nspins=Dtset%nsppol; my_spins= [(isp,isp=1,nsppol)]

 ! List of spins for each node, number of processors per each spin 
 ! and the MPI rank in the "spin" communicator.
 nprocs_spin = nprocs; rank_spin = my_rank

 if (nsppol==2.and.nprocs>1) then 
   ! Distribute spins (optimal distribution if nprocs is even)
   nprocs_spin(1) = nprocs/2
   nprocs_spin(2) = nprocs - nprocs/2
   my_nspins=1; my_spins(1)=1
   if (my_rank+1>nprocs/2) then
     my_spins(1)=2
     rank_spin = my_rank - nprocs/2 
   end if
 end if

 store_ur = (MODULO(Dtset%gwmem,10)==1) 
 bks_mask=.FALSE.; keep_ur=.FALSE.

 select case (Dtset%gwpara)
 case (1)
   ! Parallelization over transitions **without** memory distributions (Except for the spin).
   my_minb=1; my_maxb=Ep%nbnds
   do isp=1,my_nspins
     spin = my_spins(isp)
     bks_mask(my_minb:my_maxb,:,spin) = .TRUE.
     if (store_ur) keep_ur(my_minb:my_maxb,:,spin)=.TRUE.
   end do

 case (2)
   ! Distribute bands and spin.
   do isp=1,my_nspins
     spin = my_spins(isp)

     if (nprocs_spin(spin) <= nbcw) then
       ! Distribute nbcw empty bands among nprocs_spin (block of bands without replicas).
       ! Bands are distributed in contiguous blocks because
       ! this distribution is well suited for the Hilber transform 
       ! since each node will allocate only a smaller frequency interval
       ! for the spectral function whose size scales with the number of MPI nodes.
       ! Note it is now meaningless to distinguish gwcomp=0 or 1 since the workload is well balanced later on
       ABI_MALLOC(istart,(nprocs_spin(spin)))
       ABI_MALLOC(istop,(nprocs_spin(spin)))

       call xmpi_split_work2_i4b(nbcw,nprocs_spin(spin),istart,istop,msg,distr_err)

       if (distr_err==2) then
         ! Too many processors.
         !MSG_WARNING(msg)
         ierr=1
       end if

       my_minb = nbvw + istart(rank_spin+1)
       my_maxb = nbvw + istop (rank_spin+1)

       ABI_FREE(istart)
       ABI_FREE(istop)

       if (my_maxb-my_minb+1<=0) then
         write(msg,'(3a,2(i0,a),2a)')&
&         'One or more processors has zero number of bands ',ch10,&
&         'my_minb= ',my_minb,' my_maxb= ',my_maxb,ch10,&
&         'This is a waste, decrease the number of processors.'
         MSG_ERROR(msg)
       end if

       bks_mask(my_minb:my_maxb,:,spin)=.TRUE.
       if (store_ur) keep_ur(my_minb:my_maxb,:,spin)=.TRUE.

     else
       ! New version (alternate bands with replicas if nprocs > nbcw)
       ! FIXME: Fix segmentation fault with Hilbert transform.
       do ib=1,nbcw
         if (xmpi_distrib_with_replicas(ib,nbcw,rank_spin,nprocs_spin(spin))) then
           band = ib + nbvw
           bks_mask(band,:,spin)=.TRUE.
           if (store_ur) keep_ur(band,:,spin)=.TRUE.
         end if
       end do
     end if

     ! This is needed to have all the occupied states on each node.
     bks_mask(1:nbvw,:,spin) = .TRUE.
     if (store_ur) keep_ur(1:nbvw,:,spin)=.TRUE.
   end do ! isp

 case default
   ierr = 1
   MSG_WARNING("Wrong value for gwpara")
 end select

end subroutine chi0_bksmask
!!***
