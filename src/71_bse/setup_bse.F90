!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_bse
!! NAME
!!  setup_bse
!!
!! FUNCTION
!!  This routine performs the initialization of basic objects and quantities used for Bethe-Salpeter calculations.
!!  In particular the excparam data type that defines the parameters of the calculation is completely
!!  initialized starting from the content of Dtset and the parameters read from the external WFK and SCR (SUSC) file.
!!
!! COPYRIGHT
!! Copyright (C) 1992-2009 EXC group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida)
!! Copyright (C) 2009-2016 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! codvsn=Code version
!! ngfft_gw(18)=Information about 3D FFT for density and potentials, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! acell(3)=Length scales of primitive translations (bohr)
!! rprim(3,3)=Dimensionless real space primitive translations.
!! Dtset<dataset_type>=All input variables for this dataset.
!!  Some of them might be redefined here TODO
!! Dtfil=filenames and unit numbers used in abinit.
!! Psps <pseudopotential_type>=variables related to pseudopotentials
!! Pawtab(Psps%ntypat*Dtset%usepaw)<pawtab_type>=PAW tabulated starting data
!!
!! OUTPUT
!! Cryst<crystal_t>=Info on the crystalline Structure.
!! Kmesh<kmesh_t>=Structure defining the k-sampling for the wavefunctions.
!! Qmesh<kmesh_t>=Structure defining the q-sampling for the symmetrized inverse dielectric matrix.
!! Gsph_x<gsphere_t=Data type gathering info on the G-sphere for wave functions and e^{-1},
!! KS_BSt<ebands_t>=The KS band structure (energies, occupancies, k-weights...)
!! Vcp<vcoul_t>=Structure gathering information on the Coulomb interaction in reciprocal space,
!!   including a possible cutoff in real space.
!! ngfft_osc(18)=Contain all needed information about the 3D FFT for the oscillator matrix elements.
!!   See ~abinit/doc/input_variables/vargs.htm#ngfft
!! Bsp<excparam>=Basic parameters defining the Bethe-Salpeter run. Completely initialed in output.
!! Hdr_wfk<Hdr_type>=The header of the WFK file.
!! Hdr_bse<Hdr_type>=Local header initialized from the parameters used for the Bethe-Salpeter calculation.
!! BS_files<excfiles>=Files used in the calculation.
!! w_file=File name used to construct W. Set to ABI_NOFILE if no external file is used.
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      apply_scissor,bsp_calctype2str,crystal_from_hdr,crystal_print
!!      ebands_copy,ebands_init,ebands_print,ebands_report_gap
!!      ebands_update_occ,find_qmesh,get_bz_item,get_ng0sh,gsph_extend
!!      gsph_init,hdr_init,hdr_update,hdr_vs_dtset,hscr_bcast,hscr_free
!!      hscr_from_file,hscr_print,init_transitions,kmesh_init,kmesh_print
!!      make_mesh,matrginv,metric,mkrdim,pawrhoij_alloc,pawrhoij_copy
!!      pawrhoij_free,print_bs_files,print_bs_parameters,print_gsphere
!!      print_ngfft,rdgw,setmesh,vcoul_init,wfk_read_eigenvalues,wrtout
!!      xmpi_bcast,xmpi_max,xmpi_split_work
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setup_bse(codvsn,acell,rprim,ngfftf,ngfft_osc,Dtset,Dtfil,BS_files,Psps,Pawtab,BSp,&
& Cryst,Kmesh,Qmesh,KS_BSt,QP_bst,Hdr_wfk,Gsph_x,Gsph_c,Vcp,Hdr_bse,w_fname,comm,Wvl)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_bs_defs
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_nctk
 use m_hdr

 use m_gwdefs,        only : GW_Q0_DEFAULT
 use m_fstrings,      only : toupper, sjoin
 use m_io_tools,      only : file_exists, open_file
 use m_geometry,      only : normv
 use m_crystal,       only : crystal_print, idx_spatial_inversion, crystal_t
 use m_crystal_io,    only : crystal_from_hdr
 use m_bz_mesh,       only : kmesh_t, kmesh_init, get_ng0sh, kmesh_print, get_BZ_item, find_qmesh, make_mesh
 use m_ebands,        only : ebands_init, ebands_print, ebands_copy, ebands_free, &
&                            ebands_update_occ, get_valence_idx, apply_scissor, ebands_report_gap
 use m_vcoul,         only : vcoul_t, vcoul_init
 use m_fftcore,       only : print_ngfft
 use m_fft_mesh,      only : setmesh
 use m_gsphere,       only : gsphere_t, gsph_init, print_gsphere, merge_and_sort_kg, gsph_extend
 use m_io_screening,  only : hscr_t, hscr_free, hscr_io, hscr_bcast, hscr_from_file, hscr_print
 use m_pawtab,        only : pawtab_type
 use m_pawrhoij,      only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free
 use m_qparticles,    only : rdgw
 use m_screen,        only : MDL_BECHSTEDT
 use m_wfk,           only : wfk_read_eigenvalues

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_bse'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_56_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 character(len=6),intent(in) :: codvsn
 character(len=fnlen),intent(out) :: w_fname
 type(dataset_type),intent(inout) :: Dtset
 type(datafiles_type),intent(in) :: Dtfil
 type(pseudopotential_type),intent(in) :: Psps
 type(excparam),intent(inout) :: Bsp
 type(hdr_type),intent(out) :: Hdr_wfk,Hdr_bse
 type(crystal_t),intent(out) :: Cryst
 type(kmesh_t),intent(out) :: Kmesh,Qmesh
 type(gsphere_t),intent(out) :: Gsph_x,Gsph_c
 type(ebands_t),intent(out) :: KS_BSt,QP_Bst
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Dtset%usepaw)
 type(vcoul_t),intent(out) :: Vcp
 type(excfiles),intent(out) :: BS_files
 type(wvl_internal_type), intent(in) :: Wvl
!arrays
 integer,intent(in) :: ngfftf(18)
 integer,intent(out) :: ngfft_osc(18)
 real(dp),intent(in) :: acell(3),rprim(3,3)

!Local variables ------------------------------
!scalars
 integer,parameter :: pertcase0=0,master=0
 integer(i8b) :: work_size,tot_nreh,neh_per_proc,il
 integer :: bantot,enforce_sym,ib,ibtot,ik_ibz,isppol,jj,method,iat,ount !ii,
 integer :: mband,io,nfftot_osc,spin,hexc_size,nqlwl,iq
 integer :: timrev,iq_bz,isym,iq_ibz,itim
 integer :: my_rank,nprocs,fform,npwe_file,ierr,my_k1, my_k2,my_nbks
 integer :: first_dig,second_dig,it
 real(dp) :: ucvol,qnorm
 real(dp):: eff,mempercpu_mb,wfsmem_mb,nonscal_mem,ug_mem,ur_mem,cprj_mem
 logical,parameter :: remove_inv=.FALSE.
 logical :: ltest,occ_from_dtset
 character(len=500) :: msg
 character(len=fnlen) :: gw_fname,test_file,wfk_fname
 type(hscr_t) :: Hscr
!arrays
 integer :: ng0sh_opt(3),val_idx(Dtset%nsppol)
 integer,allocatable :: npwarr(:),val_indeces(:,:),nlmn_atm(:)
 real(dp) :: qpt_bz(3),minmax_tene(2)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),sq(3)
 real(dp) :: qred2cart(3,3),qcart2red(3,3)
 real(dp),allocatable :: doccde(:),eigen(:),occfact(:),qlwl(:,:)
 real(dp),allocatable :: igwene(:,:,:)
 real(dp),pointer :: energies_p(:,:,:)
 complex(dpc),allocatable :: gw_energy(:,:,:)
 type(Pawrhoij_type),allocatable :: Pawrhoij(:)
!Interp@BSE
 !integer :: mode
 !integer :: kmult(3)
 !integer :: unt
 !integer :: rl_nb
 !logical :: interp_params_exists, prepare, sum_overlaps
 !namelist /interp_params/ mode,kmult,prepare,rl_nb,sum_overlaps
 !character(len=fnlen) :: tmp_fname

!************************************************************************

 DBG_ENTER("COLL")

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 ! === Check for calculations that are not implemented ===
 ltest=ALL(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol)==Dtset%nband(1))
 ABI_CHECK(ltest,'Dtset%nband must be constant')
 ABI_CHECK(Dtset%nspinor==1,"nspinor==2 not coded")

 ! === Dimensional primitive translations rprimd (from input), gprimd, metrics and unit cell volume ===
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 ! Read energies and header from the WFK file. 
 wfk_fname = dtfil%fnamewffk 
 if (.not. file_exists(wfk_fname)) then
   wfk_fname = nctk_ncify(wfk_fname)
   MSG_COMMENT(sjoin("File not found. Will try netcdf file: ", wfk_fname))
 end if

 call wfk_read_eigenvalues(wfk_fname,energies_p,Hdr_wfk,comm)
 mband = MAXVAL(Hdr_wfk%nband)

 call hdr_vs_dtset(Hdr_wfk,Dtset)

 ! === Create crystal_t data type ===
 !remove_inv= .FALSE. !(nsym_kss/=Hdr_wfk%nsym)
 timrev=  2 ! This information is not reported in the header
            ! 1 => do not use time-reversal symmetry
            ! 2 => take advantage of time-reversal symmetry

 call crystal_from_hdr(Cryst,Hdr_wfk,timrev,remove_inv)
 call crystal_print(Cryst)
 !
 ! Setup of the k-point list and symmetry tables in the  BZ -----------------------------------
 if (Dtset%chksymbreak==0) then
   MSG_WARNING("Calling make_mesh")
   call make_mesh(Kmesh,Cryst,Dtset%kptopt,Dtset%kptrlatt,Dtset%nshiftk,Dtset%shiftk,break_symmetry=.TRUE.)
   ! TODO
   !Check if kibz from KSS file corresponds to the one returned by make_mesh.
 else
   call kmesh_init(Kmesh,Cryst,Hdr_wfk%nkpt,Hdr_wfk%kptns,Dtset%kptopt)
 end if
 BSp%nkibz = Kmesh%nibz  !We might allow for a smaller number of points....

 call kmesh_print(Kmesh,"K-mesh for the wavefunctions",std_out,Dtset%prtvol,"COLL")
 call kmesh_print(Kmesh,"K-mesh for the wavefunctions",ab_out, 0,           "COLL")

 nqlwl = 0
 w_fname = ABI_NOFILE
 if (Dtset%getscr/=0.or.Dtset%irdscr/=0) then
   w_fname=Dtfil%fnameabi_scr
 else if (Dtset%getsuscep/=0.or.Dtset%irdsuscep/=0) then
   w_fname=Dtfil%fnameabi_sus
   MSG_ERROR("(get|ird)suscep not implemented")
 end if

 if (w_fname /= ABI_NOFILE) then ! Read dimensions from the external file

   if (.not. file_exists(w_fname)) then
     w_fname = nctk_ncify(w_fname)
     MSG_COMMENT(sjoin("File not found. Will try netcdf file: ", w_fname))
   end if

   if (my_rank==master) then 
     ! Master reads npw and nqlwl from SCR file.
     call wrtout(std_out,sjoin('Testing file: ', w_fname),"COLL")
     
     call hscr_from_file(hscr,w_fname,fform,xmpi_comm_self)
     ! Echo the header.
     if (Dtset%prtvol>0) call hscr_print(Hscr)

     npwe_file = Hscr%npwe ! Have to change %npweps if it was larger than dim on disk.
     nqlwl     = Hscr%nqlwl

     if (Dtset%npweps>npwe_file) then
       write(msg,'(2(a,i0),2a,i0)')&
&       "The number of G-vectors stored on file (",npwe_file,") is smaller than Dtset%npweps = ",Dtset%npweps,ch10,&
&       "Calculation will proceed with the maximum available set, npwe_file = ",npwe_file
       MSG_WARNING(msg)
       Dtset%npweps = npwe_file
     else  if (Dtset%npweps<npwe_file) then
       write(msg,'(2(a,i0),2a,i0)')&
&       "The number of G-vectors stored on file (",npwe_file,") is larger than Dtset%npweps = ",Dtset%npweps,ch10,&
&       "Calculation will proceed with Dtset%npweps = ",Dtset%npweps
       MSG_COMMENT(msg)
     end if
   end if

   call hscr_bcast(Hscr,master,my_rank,comm)
   call xmpi_bcast(Dtset%npweps,master,comm,ierr)
   call xmpi_bcast(nqlwl,master,comm,ierr)

   if (nqlwl>0) then
     ABI_MALLOC(qlwl,(3,nqlwl))
     qlwl = Hscr%qlwl
   end if
   !
   ! Init Qmesh from the SCR file.
   call kmesh_init(Qmesh,Cryst,Hscr%nqibz,Hscr%qibz,Dtset%kptopt)

   ! The G-sphere for W and Sigma_c is initialized from the gvectors found in the SCR file.
   call gsph_init(Gsph_c,Cryst,Dtset%npweps,gvec=Hscr%gvec)

   call hscr_free(Hscr)
 else
   ! Init Qmesh from the K-mesh reported in the WFK file.
   call find_qmesh(Qmesh,Cryst,Kmesh)

   ! The G-sphere for W and Sigma_c is initialized from ecutesp.
   call gsph_init(Gsph_c,Cryst,0,ecut=Dtset%ecuteps)
   Dtset%npweps = Gsph_c%ng 
 end if

 BSp%npweps = Dtset%npweps
 BSp%ecuteps = Dtset%ecuteps

 if (nqlwl==0) then
   nqlwl=1
   ABI_MALLOC(qlwl,(3,nqlwl))
   qlwl(:,nqlwl)= GW_Q0_DEFAULT
   write(msg,'(3a,i2,a,3f9.6)')&
&    "The Header of the screening file does not contain the list of q-point for the optical limit ",ch10,&
&    "Using nqlwl= ",nqlwl," and qlwl = ",qlwl(:,1)
   MSG_COMMENT(msg)
 end if
 write(std_out,*)"nqlwl and qlwl for Coulomb singularity and e^-1",nqlwl,qlwl

 ! === Setup of q-mesh in the whole BZ ===
 ! * Stop if a nonzero umklapp is needed to reconstruct the BZ. In this case, indeed,
 !   epsilon^-1(Sq) should be symmetrized in csigme using a different expression (G-G_o is needed)
 !
 call kmesh_print(Qmesh,"Q-mesh for the screening function",std_out,Dtset%prtvol,"COLL")
 call kmesh_print(Qmesh,"Q-mesh for the screening function",ab_out ,0           ,"COLL")

 do iq_bz=1,Qmesh%nbz
   call get_BZ_item(Qmesh,iq_bz,qpt_bz,iq_ibz,isym,itim)
   sq = (3-2*itim)*MATMUL(Cryst%symrec(:,:,isym),Qmesh%ibz(:,iq_ibz))
   if (ANY(ABS(Qmesh%bz(:,iq_bz)-sq )>1.0d-4)) then
     write(std_out,*) sq,Qmesh%bz(:,iq_bz)
     write(msg,'(a,3f6.3,a,3f6.3,2a,9i3,a,i2,2a)')&
&      'qpoint ',Qmesh%bz(:,iq_bz),' is the symmetric of ',Qmesh%ibz(:,iq_ibz),ch10,&
&      'through operation ',Cryst%symrec(:,:,isym),' and itim ',itim,ch10,&
&      'however a non zero umklapp G_o vector is required and this is not yet allowed'
     MSG_ERROR(msg)
   end if
 end do

 BSp%algorithm = Dtset%bs_algorithm
 BSp%nstates   = Dtset%bs_nstates
 Bsp%nsppol    = Dtset%nsppol
 Bsp%hayd_term = Dtset%bs_hayd_term
 !
 ! Define the algorithm for solving the BSE.
 if (BSp%algorithm == BSE_ALGO_HAYDOCK) then
   BSp%niter       = Dtset%bs_haydock_niter
   BSp%haydock_tol = Dtset%bs_haydock_tol

 else if (BSp%algorithm == BSE_ALGO_CG) then
   ! FIXME For the time being use an hardcoded value.
   ! TODO change name in Dtset%
   BSp%niter       = Dtset%nstep !100
   BSp%cg_tolwfr   = Dtset%tolwfr
   BSp%nline       = Dtset%nline
   BSp%nbdbuf      = Dtset%nbdbuf
   BSp%nstates     = Dtset%bs_nstates
   MSG_WARNING("Check CG setup")
 else
   !BSp%niter       = 0
   !BSp%tol_iter    = HUGE(one)
 end if
 !
 ! Shall we include Local field effects?
 SELECT CASE (Dtset%bs_exchange_term)
 CASE (0,1)
   BSp%exchange_term = Dtset%bs_exchange_term
 CASE DEFAULT
   write(msg,'(a,i0)')" Wrong bs_exchange_term: ",Dtset%bs_exchange_term
   MSG_ERROR(msg)
 END SELECT
 !
 ! Treatment of the off-diagonal coupling block.
 SELECT CASE (Dtset%bs_coupling)
 CASE (0)
   BSp%use_coupling = 0
   msg = 'RESONANT ONLY CALCULATION'
 CASE (1)
   BSp%use_coupling = 1
   msg = ' RESONANT+COUPLING CALCULATION '
 CASE DEFAULT
   write(msg,'(a,i0)')" Wrong bs_coupling: ",Dtset%bs_coupling
   MSG_ERROR(msg)
 END SELECT
 call wrtout(std_out,msg,"COLL")

 BSp%use_diagonal_Wgg = .FALSE.
 Bsp%use_coulomb_term = .TRUE.
 BSp%eps_inf=zero
 Bsp%mdlf_type=0

 first_dig =MOD(Dtset%bs_coulomb_term,10)
 second_dig=Dtset%bs_coulomb_term/10

 Bsp%wtype = second_dig
 SELECT CASE (second_dig)
 CASE (BSE_WTYPE_NONE)
   call wrtout(std_out,"Coulomb term won't be calculated","COLL")
   Bsp%use_coulomb_term = .FALSE.

 CASE (BSE_WTYPE_FROM_SCR)
   call wrtout(std_out,"W is read from an external SCR file","COLL")
   Bsp%use_coulomb_term = .TRUE.

 CASE (BSE_WTYPE_FROM_MDL)
   call wrtout(std_out,"W is approximated with the model dielectric function","COLL")
   Bsp%use_coulomb_term = .TRUE.
   BSp%mdlf_type = MDL_BECHSTEDT
   BSp%eps_inf = Dtset%mdf_epsinf
   ABI_CHECK(Bsp%eps_inf > zero, "mdf_epsinf <= 0")

 CASE DEFAULT
   write(msg,'(a,i0)')" Wrong second digit in bs_coulomb_term: ",Dtset%bs_coulomb_term
   MSG_ERROR(msg)
 END SELECT
 !
 ! Diagonal approximation or full matrix?
 BSp%use_diagonal_Wgg = .TRUE.
 if (Bsp%wtype /= BSE_WTYPE_NONE) then
   SELECT CASE (first_dig)
   CASE (0)
     call wrtout(std_out,"Using diagonal approximation W_GG","COLL")
     BSp%use_diagonal_Wgg = .TRUE.
   CASE (1)
     call wrtout(std_out,"Using full W_GG' matrix ","COLL")
     BSp%use_diagonal_Wgg = .FALSE.
   CASE DEFAULT
     write(msg,'(a,i0)')" Wrong first digit in bs_coulomb_term: ",Dtset%bs_coulomb_term
     MSG_ERROR(msg)
   END SELECT
 end if

 !TODO move the initialization of the parameters for the interpolation in setup_bse_interp
 
 BSp%use_interp = .FALSE.
 BSp%interp_mode = BSE_INTERP_YG
 BSp%interp_kmult(1:3) = 0
 BSp%prep_interp = .FALSE.
 BSp%sum_overlaps = .TRUE. ! Sum over the overlaps

 ! Printing ncham
 BSp%prt_ncham = .FALSE.

 ! Deactivate Interpolation Technique by default
! if (.FALSE.) then

 ! Reading parameters from the input file
 BSp%use_interp = (dtset%bs_interp_mode /= 0)
 BSp%prep_interp = (dtset%bs_interp_prep == 1)

 SELECT CASE (dtset%bs_interp_mode)
 CASE (0)
   ! No interpolation, do not print anything !
 CASE (1)
   call wrtout(std_out,"Using interpolation technique with energies and wavefunctions from dense WFK","COLL")
 CASE (2)
   call wrtout(std_out,"Interpolation technique with energies and wfn on dense WFK + treatment ABC of divergence","COLL")
 CASE (3)
   call wrtout(std_out,"Interpolation technique + divergence ABC along diagonal","COLL")
 CASE (4)
   call wrtout(std_out,"Using interpolation technique mode 1 with full computation of hamiltonian","COLL")
 CASE DEFAULT
   write(msg,'(a,i0)')" Wrong interpolation mode for bs_interp_mode: ",dtset%bs_interp_mode
   MSG_ERROR(msg)
 END SELECT

 ! Read from dtset
 if(BSp%use_interp) then
   BSp%interp_method = dtset%bs_interp_method
   BSp%rl_nb = dtset%bs_interp_rl_nb
   BSp%interp_m3_width = dtset%bs_interp_m3_width
   BSp%interp_kmult(1:3) = dtset%bs_interp_kmult(1:3)
   BSp%interp_mode = dtset%bs_interp_mode
 end if

 ! Dimensions and parameters of the calculation.
 ! TODO one should add npwx as well
 !BSp%npweps=Dtset%npweps
 !BSp%npwwfn=Dtset%npwwfn

 ABI_MALLOC(Bsp%lomo_spin, (Bsp%nsppol))
 ABI_MALLOC(Bsp%homo_spin, (Bsp%nsppol)) 
 ABI_MALLOC(Bsp%lumo_spin, (Bsp%nsppol)) 
 ABI_MALLOC(Bsp%humo_spin, (Bsp%nsppol)) 
 ABI_MALLOC(Bsp%nbndv_spin, (Bsp%nsppol)) 
 ABI_MALLOC(Bsp%nbndc_spin, (Bsp%nsppol))

 ! FIXME use bs_loband(nsppol)
 Bsp%lomo_spin = Dtset%bs_loband
 write(std_out,*)"bs_loband",Dtset%bs_loband
 !if (Bsp%nsppol == 2) Bsp%lomo_spin(2) = Dtset%bs_loband

 ! Check lomo correct only for unpolarized semiconductors
 !if (Dtset%nsppol == 1 .and. Bsp%lomo > Dtset%nelect/2) then 
 !  write(msg,'(a,i0,a,f8.3)') " Bsp%lomo = ",Bsp%lomo," cannot be greater than nelect/2 = ",Dtset%nelect/2
 !  MSG_ERROR(msg)
 !end if
 !
 ! ==============================================
 ! ==== Setup of the q for the optical limit ====
 ! ==============================================
 Bsp%inclvkb = Dtset%inclvkb

 qred2cart = two_pi*Cryst%gprimd
 qcart2red = qred2cart
 call matrginv(qcart2red,3,3)

 if (Dtset%gw_nqlwl==0) then
   BSp%nq = 6
   ABI_MALLOC(BSp%q,(3,BSp%nq))
   BSp%q(:,1) = (/one,zero,zero/)  ! (100)
   BSp%q(:,2) = (/zero,one,zero/)  ! (010)
   BSp%q(:,3) = (/zero,zero,one/)  ! (001)
   BSp%q(:,4) = MATMUL(qcart2red,(/one,zero,zero/)) ! (x)
   BSp%q(:,5) = MATMUL(qcart2red,(/zero,one,zero/)) ! (y)
   BSp%q(:,6) = MATMUL(qcart2red,(/zero,zero,one/)) ! (z)
 else
   BSp%nq = Dtset%gw_nqlwl
   ABI_MALLOC(BSp%q,(3,BSp%nq))
   BSp%q = Dtset%gw_qlwl
 end if

 do iq=1,BSp%nq ! normalization
   qnorm = normv(BSp%q(:,iq),Cryst%gmet,"G")
   BSp%q(:,iq) = BSp%q(:,iq)/qnorm
 end do
 !
 ! ======================================================
 ! === Define the flags defining the calculation type ===
 ! ======================================================
 Bsp%calc_type = Dtset%bs_calctype

 BSp%mbpt_sciss = zero ! Shall we use the scissors operator to open the gap?
 if (ABS(Dtset%mbpt_sciss)>tol6) BSp%mbpt_sciss = Dtset%mbpt_sciss

!now test input parameters from input and WFK file and assume some defaults
!
! TODO Add the possibility of using a randomly shifted k-mesh with nsym>1.
! so that densities and potentials are correctly symmetrized but
! the list of the k-point in the IBZ is not expanded.

 if (mband < Dtset%nband(1)) then
   write(msg,'(2(a,i0),3a,i0)')&
&    'WFK file contains only ', mband,' levels instead of ',Dtset%nband(1),' required;',ch10,&
&    'The calculation will be done with nbands= ',mband
   MSG_WARNING(msg)
   Dtset%nband(:) = mband
 end if

 BSp%nbnds = Dtset%nband(1) ! TODO Note the change in the meaning of input variables

 if (BSp%nbnds<=Dtset%nelect/2) then
   write(msg,'(2a,a,i0,a,f8.2)')&
&    'BSp%nbnds cannot be smaller than homo ',ch10,&
&    'while BSp%nbnds = ',BSp%nbnds,' and Dtset%nelect = ',Dtset%nelect
   MSG_ERROR(msg)
 end if

!TODO add new dim for exchange part and consider the possibility of having npwsigx > npwwfn (see setup_sigma).

 ! === Build enlarged G-sphere for the exchange part ===
 call gsph_extend(Gsph_c,Cryst,Dtset%ecutwfn,Gsph_x)
 call print_gsphere(Gsph_x,unit=std_out,prtvol=Dtset%prtvol)

 ! NPWVEC as the biggest between npweps and npwwfn. MG RECHECK this part.
 !BSp%npwwfn = Dtset%npwwfn
 Bsp%npwwfn = Gsph_x%ng  ! FIXME temporary hack
 BSp%npwvec=MAX(BSp%npwwfn,BSp%npweps)
 Bsp%ecutwfn = Dtset%ecutwfn

 ! Compute Coulomb term on the largest G-sphere.
 if (Gsph_x%ng > Gsph_c%ng ) then
   call vcoul_init(Vcp,Gsph_x,Cryst,Qmesh,Kmesh,Dtset%rcut,Dtset%icutcoul,Dtset%vcutgeo,Dtset%ecutsigx,Gsph_x%ng,&  
&    nqlwl,qlwl,ngfftf,comm)
 else
   call vcoul_init(Vcp,Gsph_c,Cryst,Qmesh,Kmesh,Dtset%rcut,Dtset%icutcoul,Dtset%vcutgeo,Dtset%ecutsigx,Gsph_c%ng,&  
&    nqlwl,qlwl,ngfftf,comm)
 end if

 ABI_FREE(qlwl)

 bantot=SUM(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol))
 ABI_MALLOC(doccde,(bantot))
 ABI_MALLOC(eigen,(bantot))
 ABI_MALLOC(occfact,(bantot))
 doccde=zero; eigen=zero; occfact=zero

 ! Get occupation from input if occopt == 2
 occ_from_dtset = (Dtset%occopt == 2) 

 jj=0; ibtot=0
 do isppol=1,Dtset%nsppol
   do ik_ibz=1,Dtset%nkpt
     do ib=1,Hdr_wfk%nband(ik_ibz+(isppol-1)*Dtset%nkpt)
       ibtot=ibtot+1
       if (ib<=BSP%nbnds) then
         jj=jj+1
         eigen  (jj)=energies_p(ib,ik_ibz,isppol)
         if (occ_from_dtset) then
           occfact(jj)=Dtset%occ_orig(ibtot)
         else
           occfact(jj)=Hdr_wfk%occ(ibtot)
         end if
       end if
     end do
   end do
 end do

 ABI_FREE(energies_p)
 !
 ! * Make sure that Dtset%wtk==Kmesh%wt due to the dirty treatment of
 !   symmetry operations in the old GW code (symmorphy and inversion)
 ltest=(ALL(ABS(Dtset%wtk(1:Kmesh%nibz)-Kmesh%wt(1:Kmesh%nibz))<tol6))
 ABI_CHECK(ltest,'Mismatch between Dtset%wtk and Kmesh%wt')

 ABI_MALLOC(npwarr,(Dtset%nkpt))
 npwarr=BSP%npwwfn

 call ebands_init(bantot,KS_BSt,Dtset%nelect,doccde,eigen,Dtset%istwfk,Kmesh%ibz,Dtset%nband,&
&  Kmesh%nibz,npwarr,Dtset%nsppol,Dtset%nspinor,Dtset%tphysel,Dtset%tsmear,Dtset%occopt,occfact,Kmesh%wt,&
&  dtset%charge, dtset%kptopt, dtset%kptrlatt_orig, dtset%nshiftk_orig, dtset%shiftk_orig, &
&  dtset%kptrlatt, dtset%nshiftk, dtset%shiftk)

 ABI_FREE(doccde)
 ABI_FREE(eigen)
 ABI_FREE(npwarr)

 !TODO Occupancies are zero if NSCF. One should calculate the occupancies from the energies when
 ! the occupation scheme for semiconductors is used.
 call ebands_update_occ(KS_BSt,Dtset%spinmagntarget,prtvol=Dtset%prtvol)

 call ebands_print(KS_BSt,"Band structure read from the WFK file",unit=std_out,prtvol=Dtset%prtvol)

 call ebands_report_gap(KS_BSt,header=" KS band structure",unit=std_out,mode_paral="COLL")

 ABI_MALLOC(val_indeces,(KS_BSt%nkpt,KS_BSt%nsppol))
 val_indeces = get_valence_idx(KS_BSt)

 do spin=1,KS_BSt%nsppol
   val_idx(spin) = val_indeces(1,spin)
   write(msg,'(a,i2,a,i0)')" For spin : ",spin," val_idx ",val_idx(spin)
   call wrtout(std_out,msg,"COLL")
   if ( ANY(val_indeces(1,spin) /= val_indeces(:,spin)) ) then
     MSG_ERROR("BSE code does not support metals")
   end if
 end do

 ABI_FREE(val_indeces)
 !
 ! === Create the BSE header ===
 call hdr_init(KS_BSt,codvsn,Dtset,Hdr_bse,Pawtab,pertcase0,Psps,wvl)

 ! === Get Pawrhoij from the header of the WFK file ===
 ABI_DT_MALLOC(Pawrhoij,(Cryst%natom*Dtset%usepaw))
 if (Dtset%usepaw==1) then
   call pawrhoij_alloc(Pawrhoij,1,Dtset%nspden,Dtset%nspinor,Dtset%nsppol,Cryst%typat,pawtab=Pawtab)
   call pawrhoij_copy(Hdr_wfk%Pawrhoij,Pawrhoij)
 end if

 call hdr_update(hdr_bse,bantot,1.0d20,1.0d20,1.0d20,Cryst%rprimd,occfact,Pawrhoij,Cryst%xred,dtset%amu_orig(:,1))

 ABI_FREE(occfact)

 if (Dtset%usepaw==1) call pawrhoij_free(Pawrhoij)
 ABI_DT_FREE(Pawrhoij)

 ! === Find optimal value for G-sphere enlargment due to oscillator matrix elements ===

 ! We will split k-points over processors
 call xmpi_split_work(Kmesh%nbz,comm,my_k1,my_k2,msg,ierr)
 if (ierr/=0) then
   MSG_WARNING(msg)
 end if

 ! If there is no work to do, just skip the computation
 if (my_k2-my_k1+1 <= 0) then
   ng0sh_opt(:)=(/zero,zero,zero/)
 else
   ! * Here I have to be sure that Qmesh%bz is always inside the BZ, not always true since bz is buggy
   ! * -one is used because we loop over all the possibile differences, unlike screening
   call get_ng0sh(my_k2-my_k1+1,Kmesh%bz(:,my_k1:my_k2),Kmesh%nbz,Kmesh%bz,&
&    Qmesh%nbz,Qmesh%bz,-one,ng0sh_opt)
 end if

 call xmpi_max(ng0sh_opt,BSp%mg0,comm,ierr)

 write(msg,'(a,3(i0,1x))') ' optimal value for ng0sh = ',BSp%mg0
 call wrtout(std_out,msg,"COLL")

 ! === Setup of the FFT mesh for the oscilator strengths ===
 ! * ngfft_osc(7:18)==Dtset%ngfft(7:18) which is initialized before entering screening.
 ! * Here we redefine ngfft_osc(1:6) according to the following options :
 !
 ! method==0 --> FFT grid read from fft.in (debugging purpose)
 ! method==1 --> Normal FFT mesh
 ! method==2 --> Slightly augmented FFT grid to calculate exactly rho_tw_g (see setmesh.F90)
 ! method==3 --> Doubled FFT grid, same as the the FFT for the density,
 !
 ! enforce_sym==1 ==> Enforce a FFT mesh compatible with all the symmetry operation and FFT library
 ! enforce_sym==0 ==> Find the smallest FFT grid compatbile with the library, do not care about symmetries
 !
 ngfft_osc(1:18)=Dtset%ngfft(1:18); method=2
 if (Dtset%fftgw==00 .or. Dtset%fftgw==01) method=0
 if (Dtset%fftgw==10 .or. Dtset%fftgw==11) method=1
 if (Dtset%fftgw==20 .or. Dtset%fftgw==21) method=2
 if (Dtset%fftgw==30 .or. Dtset%fftgw==31) method=3
 enforce_sym=MOD(Dtset%fftgw,10)

 call setmesh(gmet,Gsph_x%gvec,ngfft_osc,BSp%npwvec,BSp%npweps,BSp%npwwfn,nfftot_osc,method,BSp%mg0,Cryst,enforce_sym)
 nfftot_osc=PRODUCT(ngfft_osc(1:3))

 call print_ngfft(ngfft_osc,"FFT mesh for oscillator matrix elements",std_out,"COLL",prtvol=Dtset%prtvol)
 !
 ! BSp%homo gives the 
 !BSp%homo  = val_idx(1)
 ! highest occupied band for each spin
 BSp%homo_spin = val_idx

 ! TODO generalize the code to account for this unlikely case.
 !if (Dtset%nsppol==2) then
 !  ABI_CHECK(BSp%homo == val_idx(2),"Different valence indeces for spin up and down")
 !end if

 !BSp%lumo = BSp%homo + 1
 !BSp%humo = BSp%nbnds
 !BSp%nbndv = BSp%homo  - BSp%lomo + 1
 !BSp%nbndc = BSp%nbnds - BSp%homo

 BSp%lumo_spin = BSp%homo_spin + 1
 BSp%humo_spin = BSp%nbnds
 BSp%nbndv_spin = BSp%homo_spin  - BSp%lomo_spin + 1
 BSp%nbndc_spin = BSp%nbnds - BSp%homo_spin
 BSp%maxnbndv = MAXVAL(BSp%nbndv_spin(:))
 BSp%maxnbndc = MAXVAL(BSp%nbndc_spin(:))

 BSp%nkbz = Kmesh%nbz

 call ebands_copy(KS_BSt,QP_bst)
 ABI_MALLOC(igwene,(QP_bst%mband,QP_bst%nkpt,QP_bst%nsppol))
 igwene=zero

 call bsp_calctype2str(Bsp,msg)
 call wrtout(std_out,"Calculation type: "//TRIM(msg))

 SELECT CASE (Bsp%calc_type)
 CASE (BSE_HTYPE_RPA_KS)
   if (ABS(BSp%mbpt_sciss)>tol6) then
     write(msg,'(a,f8.2,a)')' Applying a scissors operator energy= ',BSp%mbpt_sciss*Ha_eV," [eV] on top of the KS energies."
     call wrtout(std_out,msg,"COLL")
     call apply_scissor(QP_BSt,BSp%mbpt_sciss)
   else
     write(msg,'(a,f8.2,a)')' Using KS energies since mbpt_sciss= ',BSp%mbpt_sciss*Ha_eV," [eV]."
     call wrtout(std_out,msg,"COLL")
   end if

 CASE (BSE_HTYPE_RPA_QPENE) ! Read _GW files with the corrections TODO here I should introduce variable getgw
   gw_fname=TRIM(Dtfil%filnam_ds(4))//'_GW'
   gw_fname="__in.gw__"
   if (.not.file_exists(gw_fname)) then
     msg = " File "//TRIM(gw_fname)//" not found. Aborting now"
     MSG_ERROR(msg)
   end if

   call rdgw(QP_Bst,gw_fname,igwene,extrapolate=.FALSE.) ! here gwenergy is real

   do isppol=1,Dtset%nsppol
     write(std_out,*) ' k       GW energies [eV]'
     do ik_ibz=1,Kmesh%nibz
       write(std_out,'(i3,7x,10f7.2/50(10x,10f7.2/))')ik_ibz,(QP_bst%eig(ib,ik_ibz,isppol)*Ha_eV,ib=1,BSp%nbnds)
     end do
     write(std_out,*) ' k       Im GW energies [eV]'
     do ik_ibz=1,Kmesh%nibz
       write(std_out,'(i3,7x,10f7.2/50(10x,10f7.2/))')ik_ibz,(igwene(ib,ik_ibz,isppol)*Ha_eV,ib=1,BSp%nbnds)
     end do
   end do
   !
   ! If required apply the scissors operator on top of the QP bands structure (!)
   if (ABS(BSp%mbpt_sciss)>tol6) then
     write(msg,'(a,f8.2,a)')' Applying a scissors operator ',BSp%mbpt_sciss*Ha_eV," [eV] on top of the QP energies!"
     MSG_COMMENT(msg)
     call apply_scissor(QP_BSt,BSp%mbpt_sciss)
   end if

 CASE (BSE_HTYPE_RPA_QP)
   MSG_ERROR("Not implemented error!")

 CASE DEFAULT
   write(msg,'(a,i0)')"Unknown value for Bsp%calc_type= ",Bsp%calc_type
   MSG_ERROR(msg)
 END SELECT

 call ebands_report_gap(QP_BSt,header=" QP band structure",unit=std_out,mode_paral="COLL")

 ! Transitions are ALWAYS ordered in c-v-k mode with k being the slowest index.
 ! FIXME: linewidths not coded.
 ABI_MALLOC(gw_energy,(BSp%nbnds,Kmesh%nibz,Dtset%nsppol))
 gw_energy = QP_BSt%eig

 BSp%have_complex_ene = ANY(igwene > tol16)

 ! Compute the number of resonant transitions, nreh, for the two spin channels and initialize BSp%Trans.
 ABI_MALLOC(Bsp%nreh,(Bsp%nsppol))

 ! Possible cutoff on the transitions.
 BSp%ircut = Dtset%bs_eh_cutoff(1)
 BSp%uvcut = Dtset%bs_eh_cutoff(2)

 call init_transitions(BSp%Trans,BSp%lomo_spin,BSp%humo_spin,BSp%ircut,Bsp%uvcut,BSp%nkbz,Bsp%nbnds,Bsp%nkibz,&
&  BSp%nsppol,Dtset%nspinor,gw_energy,QP_BSt%occ,Kmesh%tab,minmax_tene,Bsp%nreh)

 ! Setup of the frequency mesh for the absorption spectrum.
 ! If not specified, use the min-max resonant transition energy and make it 10% smaller|larger.

 !if (ABS(Dtset%bs_freq_mesh(1)) < tol6) then
 !   Dtset%bs_freq_mesh(1) = MAX(minmax_tene(1) - minmax_tene(1) * 0.1, zero)
 !end if

 if (ABS(Dtset%bs_freq_mesh(2)) < tol6) then
    Dtset%bs_freq_mesh(2) = minmax_tene(2) + minmax_tene(2) * 0.1
 end if
                                                                                        
 Bsp%omegai = Dtset%bs_freq_mesh(1)
 Bsp%omegae = Dtset%bs_freq_mesh(2)
 Bsp%domega = Dtset%bs_freq_mesh(3)
 BSp%broad  = Dtset%zcut
                                                                                        
 ! The frequency mesh (including the complex imaginary shift)
 BSp%nomega = (BSp%omegae - BSp%omegai)/BSp%domega + 1
 ABI_MALLOC(BSp%omega,(BSp%nomega))
 do io=1,BSp%nomega
   BSp%omega(io) = (BSp%omegai + (io-1)*BSp%domega)  + j_dpc*BSp%broad
 end do

 ABI_FREE(gw_energy)
 ABI_FREE(igwene)

 do spin=1,Bsp%nsppol
   write(msg,'(a,i2,a,i0)')" For spin: ",spin,' the number of resonant e-h transitions is: ',BSp%nreh(spin)
   call wrtout(std_out,msg,"COLL")
 end do

 if (ANY(Bsp%nreh/=Bsp%nreh(1))) then
   write(msg,'(a,2(i0,1x))')" BSE code with different number of transitions for the two spin channels: ",Bsp%nreh
   MSG_WARNING(msg)
 end if
 !
 ! Create transition table vcks2t
 Bsp%lomo_min = MINVAL(BSp%lomo_spin)
 Bsp%homo_max = MAXVAL(BSp%homo_spin)
 Bsp%lumo_min = MINVAL(BSp%lumo_spin)
 Bsp%humo_max = MAXVAL(BSp%humo_spin)

 ABI_MALLOC(Bsp%vcks2t,(BSp%lomo_min:BSp%homo_max,BSp%lumo_min:BSp%humo_max,BSp%nkbz,Dtset%nsppol))
 Bsp%vcks2t = 0

 do spin=1,BSp%nsppol
   do it=1,BSp%nreh(spin)
     BSp%vcks2t(BSp%Trans(it,spin)%v,BSp%Trans(it,spin)%c,BSp%Trans(it,spin)%k,spin) = it
   end do
 end do

 hexc_size = SUM(Bsp%nreh); if (Bsp%use_coupling>0) hexc_size=2*hexc_size
 if (Bsp%nstates<=0) then
   Bsp%nstates=hexc_size
 else
   if (Bsp%nstates>hexc_size) then
      Bsp%nstates=hexc_size
      write(msg,'(2(a,i0),2a)')&
&      "Since the total size of excitonic Hamiltonian ",hexc_size," is smaller than Bsp%nstates ",Bsp%nstates,ch10,&
&      "the number of excitonic states nstates has been modified"
     MSG_WARNING(msg)
   end if
 end if

 msg=' Fundamental parameters for the solution of the Bethe-Salpeter equation:'
 call print_bs_parameters(BSp,unit=std_out,header=msg,mode_paral="COLL",prtvol=Dtset%prtvol)
 call print_bs_parameters(BSp,unit=ab_out, header=msg,mode_paral="COLL")

 if (ANY (Cryst%symrec(:,:,1) /= RESHAPE ( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )) .or. &
&    ANY( ABS(Cryst%tnons(:,1)) > tol6) ) then
   write(msg,'(3a,9i2,2a,3f6.3,2a)')&
&    "The first symmetry operation should be the Identity with zero tnons while ",ch10,&
&    "symrec(:,:,1) = ",Cryst%symrec(:,:,1),ch10,&
&    "tnons(:,1)    = ",Cryst%tnons(:,1),ch10,&
&    "This is not allowed, sym_rhotwgq0 should be changed."
   MSG_ERROR(msg)
 end if
 !
 ! Prefix for generic output files.
 BS_files%out_basename = TRIM(Dtfil%filnam_ds(4))
 !
 ! Search for files to restart from.
 if (Dtset%gethaydock/=0 .or. Dtset%irdhaydock/=0) then
   BS_files%in_haydock_basename = TRIM(Dtfil%fnameabi_haydock)
 end if

 test_file = Dtfil%fnameabi_bsham_reso
 if (file_exists(test_file)) then
   BS_files%in_hreso = test_file
 else
   BS_files%out_hreso = TRIM(Dtfil%filnam_ds(4))//'_BSR'
 end if

 test_file = Dtfil%fnameabi_bsham_coup
 if (file_exists(test_file) ) then
   BS_files%in_hcoup = test_file
 else
   BS_files%out_hcoup = TRIM(Dtfil%filnam_ds(4))//'_BSC'
 end if
 !
 ! in_eig is the name of the input file with eigenvalues and eigenvectors
 ! constructed from getbseig or irdbseig. out_eig is the name of the output file
 ! produced by this dataset. in_eig_exists checks for the presence of the input file.
 !
 if (file_exists(Dtfil%fnameabi_bseig)) then
   BS_files%in_eig = Dtfil%fnameabi_bseig
 else
   BS_files%out_eig = TRIM(BS_files%out_basename)//"_BSEIG"
 end if

 call print_bs_files(BS_files,unit=std_out)

 !
 ! ==========================================================
 ! ==== Final check on the parameters of the calculation ====
 ! ==========================================================
 if ( Bsp%use_coupling>0 .and. ALL(Bsp%algorithm /= [BSE_ALGO_DDIAGO, BSE_ALGO_HAYDOCK]) ) then
   msg = "Resonant+Coupling is only available with the direct diagonalization or the haydock method."
   MSG_ERROR(msg)
 end if

 ! autoparal section
 if (dtset%max_ncpus /=0 .and. dtset%autoparal /=0 ) then
   ount = ab_out
   ! TODO:
   ! nsppol and calculation with coupling!

   ! Temporary table needed to estimate memory
   ABI_MALLOC(nlmn_atm,(Cryst%natom))
   if (Dtset%usepaw==1) then
     do iat=1,Cryst%natom
       nlmn_atm(iat)=Pawtab(Cryst%typat(iat))%lmn_size
     end do
   end if

   tot_nreh = SUM(BSp%nreh)
   work_size = tot_nreh * (tot_nreh + 1) / 2

   write(ount,'(a)')"--- !Autoparal"
   write(ount,"(a)")'#Autoparal section for Bethe-Salpeter runs.'

   write(ount,"(a)")   "info:"
   write(ount,"(a,i0)")"    autoparal: ",dtset%autoparal
   write(ount,"(a,i0)")"    max_ncpus: ",dtset%max_ncpus
   write(ount,"(a,i0)")"    nkibz: ",Bsp%nkibz
   write(ount,"(a,i0)")"    nkbz: ",Bsp%nkbz
   write(ount,"(a,i0)")"    nsppol: ",dtset%nsppol
   write(ount,"(a,i0)")"    nspinor: ",dtset%nspinor
   write(ount,"(a,i0)")"    lomo_min: ",Bsp%lomo_min
   write(ount,"(a,i0)")"    humo_max: ",Bsp%humo_max
   write(ount,"(a,i0)")"    tot_nreh: ",tot_nreh
   !write(ount,"(a,i0)")"    nbnds: ",Ep%nbnds

   ! Wavefunctions are not distributed. We read all the bands 
   ! from 1 up to Bsp%nbnds because we have to recompute rhor 
   ! but then we deallocate all the states that are not used for the construction of the e-h 
   ! before allocating the EXC hamiltonian. Hence we can safely use  (humo - lomo + 1) instead of Bsp%nbnds.
   !my_nbks = (Bsp%humo - Bsp%lomo +1) * Bsp%nkibz * Dtset%nsppol
 
   ! This one overestimates the memory but it seems to be safer.
   my_nbks = Bsp%nbnds * Dtset%nkpt * Dtset%nsppol

   ! Memory needed for Fourier components ug.
   ug_mem = two*gwpc*Dtset%nspinor*Bsp%npwwfn*my_nbks*b2Mb
                                                                               
   ! Memory needed for real space ur.
   ur_mem = zero
   if (MODULO(Dtset%gwmem,10)==1) then
     ur_mem = two*gwpc*Dtset%nspinor*nfftot_osc*my_nbks*b2Mb
   end if
                                                                               
   ! Memory needed for PAW projections Cprj 
   cprj_mem = zero
   if (Dtset%usepaw==1) cprj_mem = dp*Dtset%nspinor*SUM(nlmn_atm)*my_nbks*b2Mb
                                                                               
   wfsmem_mb = ug_mem + ur_mem + cprj_mem

   ! Non-scalable memory in Mb i.e. memory that is not distributed with MPI:  wavefunctions + W
   nonscal_mem = (wfsmem_mb + two*gwpc*BSp%npweps**2*b2Mb) * 1.1_dp

   ! List of configurations.
   write(ount,"(a)")"configurations:"
   do il=1,dtset%max_ncpus
     if (il > work_size) cycle
     neh_per_proc = work_size / il
     neh_per_proc = neh_per_proc + MOD(work_size, il)
     eff = (one * work_size) / (il * neh_per_proc)

     ! EXC matrix is distributed.
     mempercpu_mb = nonscal_mem + two * dpc * neh_per_proc * b2Mb

     write(ount,"(a,i0)")"    - tot_ncpus: ",il
     write(ount,"(a,i0)")"      mpi_ncpus: ",il
     !write(ount,"(a,i0)")"      omp_ncpus: ",omp_ncpus
     write(ount,"(a,f12.9)")"      efficiency: ",eff
     write(ount,"(a,f12.2)")"      mem_per_cpu: ",mempercpu_mb 
   end do

   write(ount,'(a)')"..."

   ABI_FREE(nlmn_atm)
   MSG_ERROR_NODUMP("aborting now")
 end if

 DBG_EXIT("COLL")

end subroutine setup_bse
!!***
