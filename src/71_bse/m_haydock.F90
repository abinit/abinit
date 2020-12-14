!!****m* ABINIT/m_haydock
!! NAME
!! m_haydock
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (M.Giantomassi, Y. Gillet, L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida)
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

MODULE m_haydock

 use defs_basis
 use m_abicore
 use m_bs_defs
 use m_xmpi
 use m_errors
 use m_nctk
 use m_haydock_io
 use m_linalg_interfaces
 use m_ebands
 use m_hdr
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_time,              only : timab
 use m_fstrings,          only : strcat, sjoin, itoa, int2char4
 use m_io_tools,          only : file_exists, open_file
 use defs_datatypes,      only : ebands_t, pseudopotential_type
 use m_geometry,          only : normv
 use m_hide_blas,         only : xdotc, xgemv
 use m_hide_lapack,       only : matrginv
 use m_numeric_tools,     only : print_arr, symmetrize, hermitianize, continued_fract, wrap2_pmhalf, iseven
 use m_kpts,              only : listkk
 use m_crystal,           only : crystal_t
 use m_bz_mesh,           only : kmesh_t, findqg0, get_bz_item
 use m_double_grid,       only : double_grid_t, get_kpt_from_indices_coarse, compute_corresp
 use m_paw_hr,            only : pawhur_t
 use m_wfd,               only : wfd_t
 use m_bse_io,            only : exc_write_optme
 use m_pawtab,            only : pawtab_type
 use m_vcoul,             only : vcoul_t
 use m_hexc,              only : hexc_init, hexc_interp_init, hexc_free, hexc_interp_free, &
&                                hexc_build_hinterp, hexc_matmul_tda, hexc_matmul_full, hexc_t, hexc_matmul_elphon, hexc_interp_t
 use m_exc_spectra,       only : exc_write_data, exc_eps_rpa, exc_write_tensor, mdfs_ncwrite
 use m_eprenorms,         only : eprenorms_t, renorm_bst
 use m_wfd_optic,         only : calc_optical_mels

 implicit none

 private
!!***

 public :: exc_haydock_driver     ! Driver for the Haydock method (main entry point for client code).

CONTAINS  !=======================================================================
!!***

!!****f* m_haydock/exc_haydock_driver
!! NAME
!! exc_haydock_driver
!!
!! FUNCTION
!!  Calculate the imaginary part of the macroscopic dielectric function with the Haydock recursive method.
!!
!! INPUTS
!! BSp<type(excparam)=The parameter for the Bethe-Salpeter run.
!! BS_files<excparam>=Files associated to the bethe_salpeter code.
!! Cryst<crystal_t>=Info on the crystalline structure.
!! Kmesh<type(kmesh_t)>=The list of k-points in the BZ, IBZ and symmetry tables.
!! Cryst<type(crystal_t)>=Info on the crystalline structure.
!! Hdr_bse
!! KS_BSt=The KS energies.
!! QP_BSt=The QP energies.
!! Wfd<wfd_t>=Wavefunction descriptor (input k-mesh)
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials.
!! Pawtab(Cryst%ntypat*usepaw)<pawtab_type>=PAW tabulated starting data.
!! Hur(Cryst%natom*usepaw)<type(pawhur_t)>=Only for PAW and DFT+U, quantities used to evaluate the commutator [H_u,r].
!!
!! OUTPUT
!!  The imaginary part of the macroscopic dielectric function is written on the external file _EXC_MDF
!!
!! PARENTS
!!      m_bethe_salpeter
!!
!! CHILDREN
!!
!! SOURCE

subroutine exc_haydock_driver(BSp,BS_files,Cryst,Kmesh,Hdr_bse,KS_BSt,QP_Bst,Wfd,Psps,Pawtab,Hur,Epren,&
& Kmesh_dense, KS_BSt_dense, QP_BSt_dense, Wfd_dense, Vcp_dense, grid)

!Arguments ------------------------------------
!scalars
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(kmesh_t),intent(in) :: Kmesh
 type(crystal_t),intent(in) :: Cryst
 type(Hdr_type),intent(in) :: Hdr_bse
 type(wfd_t),intent(inout) :: Wfd
 type(pseudopotential_type),intent(in) :: Psps
 type(ebands_t),intent(in) :: KS_BSt,QP_Bst
 type(double_grid_t),intent(in),optional :: grid
 type(kmesh_t),intent(in),optional :: Kmesh_dense
 type(wfd_t),intent(inout),optional :: Wfd_dense
 type(ebands_t),intent(in),optional :: KS_BSt_dense, QP_Bst_dense
 type(vcoul_t),intent(in),optional :: Vcp_dense
 type(eprenorms_t),intent(in) :: Epren
!arrays
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(pawhur_t),intent(in) :: Hur(Cryst%natom*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: io,my_rank,iq,itt,ierr
 integer :: hsize,comm,my_t1,my_t2,nsppol,nkets,nproc,ncid
 integer :: spin,spad,ik_bz,iv,ic,trans_idx,lomo_min,max_band
 real(dp) :: omegaev,rand_phi !,norm
 complex(dpc) :: ks_avg,gw_avg,exc_avg
 logical :: use_mpio,prtdos
 character(len=500) :: msg
 type(hexc_t) :: hexc
 type(hexc_interp_t) :: hexc_i
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: dos(:),dos_gw(:),dos_ks(:)
 complex(dpc),allocatable :: green(:,:)
 complex(dpc),allocatable :: opt_cvk(:,:,:,:,:),kets(:,:)
 complex(dpc),allocatable :: eps_rpanlf(:,:),eps_gwnlf(:,:)
 complex(dpc),allocatable :: tensor_cart(:,:),tensor_cart_rpanlf(:,:),tensor_cart_gwnlf(:,:)
 complex(dpc),allocatable :: tensor_red(:,:),tensor_red_rpanlf(:,:),tensor_red_gwnlf(:,:)

 !Temperature
 real(dp) :: dksqmax, en
 integer,allocatable :: bs2eph(:,:)
 integer :: sppoldbl, timrev
 logical :: do_ep_renorm, do_ep_lifetime
 integer :: ntemp
 character(len=4) :: ts
 character(len=fnlen) :: prefix
 character(len=fnlen) :: path
 complex(dpc),allocatable :: ep_renorms(:)
 integer :: ep_ik, ik, ireh, isppol
 integer :: itemp
 type(ebands_t) :: EPBSt, EP_QPBSt

!************************************************************************

 call timab(690,1,tsec) ! exc_haydock_driver
 call timab(691,1,tsec) ! exc_haydock_driver(read)

 if (BSp%have_complex_ene) then
   MSG_ERROR("Complex energies are not supported yet")
 end if

 my_rank = Wfd%my_rank
 comm    = Wfd%comm
 nsppol  = Wfd%nsppol
 nproc   = Wfd%nproc

 use_mpio=.FALSE.
#ifdef HAVE_MPI_IO
 use_mpio = (nproc > 1)
 !use_mpio = .TRUE.
#endif
 use_mpio=.FALSE.
 !use_mpio = .TRUE.

 ! Hsize refers to the size of the individual blocks (resonant and coupling).
 ! Thanks to the symmetry property of the starting vector, the Haydock method
 ! can be reformulated in terms of matrix-vector multiplication involving the
 ! blocks thus avoiding to allocation of the full matrix ( R   C )
 !                                                        -C* -R*)
 hsize=SUM(BSp%nreh)

 !YG2014
 call hexc_init(hexc, BSp, BS_files, Cryst, Kmesh, Wfd, KS_BSt, QP_BSt, comm)

 !YG2014
 if(BSp%use_interp) then
   call hexc_interp_init(hexc_i, hexc, BSp%interp_m3_width, BSp%interp_method,&
&     Kmesh_dense, Vcp_dense, grid, Wfd_dense, &
&     KS_BSt_dense, QP_BSt_dense, Psps, Pawtab)

 end if

 call timab(691,2,tsec) ! exc_haydock_driver(read)
 call timab(692,1,tsec) ! exc_haydock_driver(prep)
 !
 ! Prepare the starting vectors for the Lanczos chain.
 nkets=Bsp%nq

 prtdos=.FALSE. !prtdos=.TRUE.
 if (prtdos) then
   nkets=nkets+1
   if (Bsp%use_coupling>0) then
     MSG_ERROR("DOS with coupling not coded")
     nkets=nkets+1
   end if
 end if

 !YG2014
 ABI_MALLOC_OR_DIE(kets,(hexc%hsize,nkets), ierr)
 kets=czero
 !
 ! Prepare the kets for the macroscopic dielectric function.
 lomo_min=Bsp%lomo_min; max_band=Bsp%nbnds

 !YG2014
 ABI_MALLOC_OR_DIE(opt_cvk,(lomo_min:max_band,lomo_min:max_band,hexc%nbz,Wfd%nsppol,BSp%nq), ierr)

 do iq=1,Bsp%nq
 ! Note KS_BSt is used here to calculate the commutator.
   call calc_optical_mels(hexc%Wfd,hexc%Kmesh,hexc%KS_BSt,Cryst,Psps,Pawtab,Hur, &
&     BSp%inclvkb,BSp%lomo_spin,lomo_min,max_band,hexc%nbz,BSp%q(:,iq),opt_cvk(:,:,:,:,iq))

 ! Fill ket0 using the same ordering for the indeces as the one used for the excitonic Hamiltonian.
 ! Note that only the resonant part is used here.
   do spin=1,nsppol

     if(BSp%use_interp) then
       spad=(spin-1)*BSp%nreh_interp(spin)
     else
       spad=(spin-1)*BSp%nreh(spin)
     end if

     do ik_bz=1,hexc%nbz
       do iv=BSp%lomo_spin(spin),BSp%homo_spin(spin)
         do ic=BSp%lumo_spin(spin),BSp%nbnds

           if(BSp%use_interp) then
             trans_idx = BSp%vcks2t_interp(iv,ic,ik_bz,spin)
           else
             trans_idx = BSp%vcks2t(iv,ic,ik_bz,spin)
           end if

           if (trans_idx>0) kets(trans_idx+spad,iq)=opt_cvk(ic,iv,ik_bz,spin,iq)
         end do
       end do
     end do
   end do

 end do

 !
 ! ========================================================
 ! === Write the Optical Matrix Elements to NetCDF file ===
 ! ========================================================

 !if (.false.) then
 !  ome_fname='test_OME.nc'
 !  call exc_write_optme(ome_fname,minb,maxb,BSp%nkbz,Wfd%nsppol,BSp%nq,opt_cvk,ierr)
 !end if

 ! Free WFD descriptor, we don't need ur and ug anymore !
 ! We make space for interpolated hamiltonian
 call wfd%wave_free("All")
 if(BSp%use_interp) call wfd_dense%wave_free("All")

 ! Build interpolated hamiltonian
 if(BSp%use_interp) then
   if (any(BSp%interp_mode == [2,3,4])) then
     call hexc_build_hinterp(hexc, hexc_i)
   end if
 end if

 call timab(692,2,tsec) ! exc_haydock_driver(prep)
 call timab(693,1,tsec) ! exc_haydock_driver(wo lf) - that is, without local field

 do_ep_renorm = .FALSE.
 ntemp = 1
 do_ep_lifetime = .FALSE.

 if(BSp%do_ep_renorm) then
   if (BSp%nsppol == 2) then
     MSG_ERROR('Elphon renorm with nsppol == 2 not yet coded !')
   end if
   do_ep_renorm = .TRUE.
   ntemp = Epren%ntemp
   if(BSp%do_lifetime) then
     do_ep_lifetime = .TRUE.
   end if

   ! Force elphon linewidth
   do_ep_lifetime = .TRUE.

   ! Map points from BSE to elphon kpoints
   sppoldbl = 1 !; if (any(Cryst%symafm == -1) .and. Epren%nsppol == 1) nsppoldbl=2
   ABI_MALLOC(bs2eph, (Kmesh%nbz*sppoldbl, 6))
   ABI_MALLOC(ep_renorms, (hsize))
   timrev = 1
   call listkk(dksqmax, Cryst%gmet, bs2eph, Epren%kpts, Kmesh%bz, Epren%nkpt, Kmesh%nbz, Cryst%nsym, &
      sppoldbl, Cryst%symafm, Cryst%symrel, timrev, comm, use_symrec=.False.)
 end if

 call timab(693,2,tsec) ! exc_haydock_driver(wo lf    - that is, without local field
 call timab(694,1,tsec) ! exc_haydock_driver(apply

 prefix = ""
 do itemp = 1, ntemp

   call ebands_copy(hexc%KS_BSt, EPBSt)
   call ebands_copy(hexc%QP_BSt, EP_QPBSt)

   ! =================================================
   ! == Calculate elphon vector in transition space ==
   ! =================================================
   if (do_ep_renorm) then

     ! Will perform elphon renormalization for itemp

     call int2char4(itemp,ts)
     prefix = TRIM("_T") // ts

     ! No scissor with KSBST
     call renorm_bst(Epren, EPBSt, Cryst, itemp, do_lifetime=.TRUE.,do_check=.TRUE.)

     call renorm_bst(Epren, EP_QPBSt, Cryst, itemp, do_lifetime=.TRUE.,do_check=.FALSE.)

     do isppol = 1, BSp%nsppol
       do ireh = 1, BSp%nreh(isppol)
         ic = BSp%Trans(ireh,isppol)%c
         iv = BSp%Trans(ireh,isppol)%v
         ik = BSp%Trans(ireh,isppol)%k ! In the full bz
         en = BSp%Trans(ireh,isppol)%en

         ep_ik = bs2eph(ik,1)

         !TODO support multiple spins !
         if(ABS(en - (Epren%eigens(ic,ep_ik,isppol)-Epren%eigens(iv,ep_ik,isppol)+BSp%mbpt_sciss)) > tol3) then
           MSG_ERROR("Eigen from the transition does not correspond to the EP file !")
         end if
         ep_renorms(ireh) = (Epren%renorms(1,ic,ik,isppol,itemp) - Epren%renorms(1,iv,ik,isppol,itemp))

         ! Add linewith
         if(do_ep_lifetime) then
           ep_renorms(ireh) = ep_renorms(ireh) - j_dpc*(Epren%linewidth(1,ic,ik,isppol,itemp) +&
&                                                       Epren%linewidth(1,iv,ik,isppol,itemp))
         end if

       end do
     end do
   end if


   ! =======================================================
   ! === Make EPS RPA and GW without local-field effects ===
   ! =======================================================
   ABI_MALLOC(eps_rpanlf,(BSp%nomega,BSp%nq))
   ABI_MALLOC(dos_ks,(BSp%nomega))
   ABI_MALLOC(eps_gwnlf ,(BSp%nomega,BSp%nq))
   ABI_MALLOC(dos_gw,(BSp%nomega))

   call wrtout(std_out," Calculating RPA NLF and QP NLF epsilon","COLL")

   call exc_eps_rpa(BSp%nbnds,BSp%lomo_spin,BSp%lomo_min,BSp%homo_spin,hexc%Kmesh,EPBSt,BSp%nq,nsppol,&
&    opt_cvk,Cryst%ucvol,BSp%broad,BSp%nomega,BSp%omega,eps_rpanlf,dos_ks)

   call exc_eps_rpa(BSp%nbnds,BSp%lomo_spin,BSp%lomo_min,BSp%homo_spin,hexc%Kmesh,EP_QPBSt,BSp%nq,nsppol,&
&    opt_cvk,Cryst%ucvol,Bsp%broad,BSp%nomega,BSp%omega,eps_gwnlf,dos_gw)

   if (my_rank==master) then ! Only master works.
     !
     ! Master node writes final results on file.
     call exc_write_data(BSp,BS_files,"RPA_NLF_MDF",eps_rpanlf,prefix=prefix,dos=dos_ks)

     call exc_write_data(BSp,BS_files,"GW_NLF_MDF",eps_gwnlf,prefix=prefix,dos=dos_gw)

     ! Computing and writing tensor in files

     ! RPA_NLF
     ABI_MALLOC(tensor_cart_rpanlf,(BSp%nomega,6))
     ABI_MALLOC(tensor_red_rpanlf,(BSp%nomega,6))

     call wrtout(std_out," Calculating RPA NLF dielectric tensor","COLL")
     call haydock_mdf_to_tensor(BSp,Cryst,eps_rpanlf,tensor_cart_rpanlf, tensor_red_rpanlf, ierr)

     if(ierr == 0) then
        ! Writing tensor
        call exc_write_tensor(BSp,BS_files,"RPA_NLF_TSR_CART",tensor_cart_rpanlf)
        call exc_write_tensor(BSp,BS_files,"RPA_NLF_TSR_RED",tensor_red_rpanlf)
     else
        write(msg,'(3a)')&
&         'The RPA_NLF dielectric complex tensor cannot be computed',ch10,&
&         'There must be 6 different q-points in long wavelength limit (see gw_nqlwl)'
        MSG_COMMENT(msg)
     end if

     ABI_FREE(tensor_cart_rpanlf)
     ABI_FREE(tensor_red_rpanlf)

     ! GW_NLF
     ABI_MALLOC(tensor_cart_gwnlf,(BSp%nomega,6))
     ABI_MALLOC(tensor_red_gwnlf,(BSp%nomega,6))

     call wrtout(std_out," Calculating GW NLF dielectric tensor","COLL")

     call haydock_mdf_to_tensor(BSp,Cryst,eps_gwnlf,tensor_cart_gwnlf, tensor_red_gwnlf, ierr)

     if(ierr == 0) then
        ! Writing tensor
        call exc_write_tensor(BSp,BS_files,"GW_NLF_TSR_CART",tensor_cart_gwnlf)
        call exc_write_tensor(BSp,BS_files,"GW_NLF_TSR_RED",tensor_red_gwnlf)
     else
        write(msg,'(3a)')&
&         'The GW_NLF dielectric complex tensor cannot be computed',ch10,&
&         'There must be 6 different q-points in long wavelength limit (see gw_nqlwl)'
        MSG_COMMENT(msg)
     end if

     ABI_FREE(tensor_cart_gwnlf)
     ABI_FREE(tensor_red_gwnlf)

     !call wrtout(std_out," Checking Kramers Kronig on Excitonic Macroscopic Epsilon","COLL")
     !call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_exc(:,1))

     !call wrtout(std_out," Checking Kramers Kronig on RPA NLF Macroscopic Epsilon","COLL")
     !call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_rpanlf(:,1))

     !call wrtout(std_out," Checking Kramers Kronig on GW NLF Macroscopic Epsilon","COLL")
     !call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_gwnlf(:,1))

     !call wrtout(std_out," Checking f-sum rule on Excitonic Macroscopic Epsilon","COLL")

     !if (BSp%exchange_term>0) then
     !  MSG_COMMENT(' f-sum rule should be checked without LF')
     !end if
     !call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_exc(:,1)),drude_plsmf)

     !call wrtout(std_out," Checking f-sum rule on RPA NLF Macroscopic Epsilon","COLL")
     !call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_rpanlf(:,1)),drude_plsmf)

     !call wrtout(std_out," Checking f-sum rule on GW NLF Macroscopic Epsilon","COLL")
     !call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_gwnlf(:,1)),drude_plsmf)
   end if ! my_rank==master

   !call xmpi_barrier(comm)
   !
   ! The ket for the approximated DOS.
   if (prtdos) then
     MSG_WARNING("Calculating DOS with Haydock method")
     ABI_CHECK(BSp%use_coupling==0,"DOS with coupling not coded")
     iq = BSp%nq + 1
     if (my_rank==master) then
       !call random_seed()
       do itt=1,SUM(Bsp%nreh)
         call RANDOM_NUMBER(rand_phi)
         rand_phi = two_pi*rand_phi
         kets(itt,iq) = CMPLX( COS(rand_phi), SIN(rand_phi) )
       end do
       ! Normalize the vector.
       !norm = SQRT( DOT_PRODUCT(kets(:,iq), kets(:,iq)) )
       !kets(:,iq) = kets(:,iq)/norm
     end if
     call xmpi_bcast(kets(:,iq),master,comm,ierr)
   end if

   ABI_MALLOC(green,(BSp%nomega,nkets))

   if (BSp%use_coupling==0) then
     if(do_ep_renorm) then
       call haydock_bilanczos(BSp,BS_files,Cryst,Hdr_bse,hexc,hexc_i,hsize,hexc%my_t1,hexc%my_t2,nkets,kets,ep_renorms,green,comm)
     else
       !YG2014
       call haydock_herm(BSp,BS_files,Cryst,Hdr_bse,hexc%my_t1,hexc%my_t2,&
&         nkets,kets,green,hexc,hexc_i,comm)
     end if
   else
     if (BSp%use_interp) then
       MSG_ERROR("BSE Interpolation with coupling is not supported")
     else
       call haydock_psherm(BSp,BS_files,Cryst,Hdr_bse,hexc,hexc_i,hsize,my_t1,my_t2,nkets,kets,green,comm)
     end if
   end if

   !
   ! Add 1 to have the real part right.
   green = one + green

   if (my_rank==master) then ! Master writes the final results.
     !
     if (prtdos) then
       ABI_MALLOC(dos,(BSp%nomega))
       dos = -AIMAG(green(:,BSp%nq+1))
       call exc_write_data(BSp,BS_files,"EXC_MDF",green,prefix=prefix,dos=dos)
       ABI_FREE(dos)
     else
       call exc_write_data(BSp,BS_files,"EXC_MDF",green,prefix=prefix)
     end if
     !
     ! =========================
     ! === Write out Epsilon ===
     ! =========================

     ABI_MALLOC(tensor_cart,(BSp%nomega,6))
     ABI_MALLOC(tensor_red,(BSp%nomega,6))

     call wrtout(std_out," Calculating EXC dielectric tensor","COLL")
     call haydock_mdf_to_tensor(BSp,Cryst,green,tensor_cart,tensor_red,ierr)

     if (ierr == 0) then
         ! Writing tensor
         call exc_write_tensor(BSp,BS_files,"EXC_TSR_CART",tensor_cart)
         call exc_write_tensor(BSp,BS_files,"EXC_TSR_RED",tensor_red)
     else
         write(msg,'(3a)')&
&          'The EXC dielectric complex tensor cannot be computed',ch10,&
&          'There must be 6 different q-points in long wavelength limit (see gw_nqlwl)'
         MSG_COMMENT(msg)
     end if

     ABI_FREE(tensor_cart)
     ABI_FREE(tensor_red)
     !
     ! This part will be removed when fldiff will be able to compare two mdf files.
     write(ab_out,*)" "
     write(ab_out,*)"Macroscopic dielectric function:"
     write(ab_out,*)"omega [eV] <KS_RPA_nlf>  <GW_RPA_nlf>  <BSE> "
     do io=1,MIN(BSp%nomega,10)
       omegaev = REAL(BSp%omega(io))*Ha_eV
       ks_avg  = SUM( eps_rpanlf(io,:)) / Bsp%nq
       gw_avg  = SUM( eps_gwnlf (io,:)) / Bsp%nq
       exc_avg = SUM( green     (io,:)) / BSp%nq
       write(ab_out,'(7f9.4)')omegaev,ks_avg,gw_avg,exc_avg
     end do
     write(ab_out,*)" "

     ! Write MDF file with the final results.
     ! FIXME: It won't work if prtdos == True
#ifdef HAVE_NETCDF
     path = strcat(BS_files%out_basename,strcat(prefix,"_MDF.nc"))
     NCF_CHECK(nctk_open_create(ncid, path, xmpi_comm_self))
     NCF_CHECK(cryst%ncwrite(ncid))
     NCF_CHECK(ebands_ncwrite(QP_bst, ncid))
     call mdfs_ncwrite(ncid, Bsp, green, eps_rpanlf, eps_gwnlf)
     NCF_CHECK(nf90_close(ncid))
#else
     ABI_UNUSED(ncid)
#endif
   end if

   ABI_FREE(green)
   ABI_FREE(eps_rpanlf)
   ABI_FREE(eps_gwnlf)
   ABI_FREE(dos_ks)
   ABI_FREE(dos_gw)

   call ebands_free(EPBSt)
   call ebands_free(EP_QPBst)

 end do ! itemp loop

 ABI_FREE(opt_cvk)

 ABI_FREE(kets)

 call timab(694,2,tsec) ! exc_haydock_driver(apply
 call timab(695,1,tsec) ! exc_haydock_driver(end)

 !YG2014
 call hexc_free(hexc)
 call hexc_interp_free(hexc_i)

 if (do_ep_renorm) then
   ABI_FREE(ep_renorms)
   ABI_FREE(bs2eph)
 end if

 call timab(695,2,tsec) ! exc_haydock_driver(end)
 call timab(690,2,tsec) ! exc_haydock_driver

end subroutine exc_haydock_driver
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_herm
!! NAME
!! haydock_herm
!!
!! FUNCTION
!!  Reads the excitonic Hamiltonian from file and construct the Lanczos set of vectors
!!  by iterative matrix-vector multiplications.
!!
!! INPUTS
!! BSp<excparam>=Parameters for the Bethe-Salpeter calculation.
!! BS_files<excparam>=Files associated to the bethe_salpeter code.
!! Cryst<crystal_t>=Info on the crystalline structure.
!! Pawtab(Cryst%ntypat*usepaw)<pawtab_type>=PAW tabulated starting data.
!! hize=Size of the excitonic matrix.
!! my_t1,my_t2=First and last columns treated by this node.
!! nkets=Number of starting vectors for Haydock method.
!! kets(hsize,nkets)=The kets in the eh representation.
!! comm=MPI communicator.
!!
!! OUTPUT
!!  green(BSp%nomega,nkets)=
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine haydock_herm(BSp,BS_files,Cryst,Hdr_bse,my_t1,my_t2,&
& nkets,kets,green,hexc,hexc_i,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_t1,my_t2,nkets,comm
 type(crystal_t),intent(in) :: Cryst
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(Hdr_type),intent(in) :: Hdr_bse
 type(hexc_t),intent(inout) :: hexc
 type(hexc_interp_t),intent(inout) :: hexc_i
!arrays
 complex(dp),intent(out) :: green(BSp%nomega,nkets)
 complex(dpc),intent(in) :: kets(hexc%hsize,nkets)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: inn,nproc,my_rank,ierr
 integer :: niter_file,niter_max,niter_done,nsppol,iq,my_nt,term_type,n_all_omegas
 real(dp) :: norm,nfact
 logical :: can_restart,is_converged
 complex(dpc) :: factor
 character(len=500) :: msg
 character(len=fnlen),parameter :: tag_file="_HAYDR_SAVE"
 character(len=fnlen) :: restart_file,out_file
 type(haydock_type) :: haydock_file
!arrays
 real(dp),allocatable :: bb_file(:)
 real(dp),allocatable :: bb(:)
 complex(dpc),allocatable :: aa(:),phi_nm1(:),phi_n(:),hphi_n(:),hphi_nm1(:)
 complex(dpc),allocatable :: aa_file(:),phi_n_file(:),phi_nm1_file(:)
 complex(dpc),allocatable :: ket0(:),all_omegas(:),green_temp(:,:)
! complex(dpc),allocatable :: diag_dense(:)
 logical :: check(2)

!************************************************************************

 ABI_CHECK(Bsp%nsppol==1,"nsppol > 1 not implemented yet")

 nproc  = xmpi_comm_size(comm); my_rank= xmpi_comm_rank(comm)
 nsppol = Hdr_bse%nsppol

 if (BSp%use_interp) then
   MSG_COMMENT("No parallelization in Interpolation")
   my_nt = SUM(Bsp%nreh_interp)
 else
   my_nt = my_t2-my_t1+1
 end if

 ABI_CHECK(my_nt>0,"One of the processors has zero columns")

 write(msg,'(a,i0)')' Haydock algorithm with MAX number of iterations: ',BSp%niter
 call wrtout(std_out,msg,"COLL")
 !
 ! Select the terminator for the continued fraction.
 term_type=0; if (Bsp%hayd_term>0) term_type=1
 call wrtout(std_out,sjoin("Using terminator type: ",itoa(term_type)),"COLL")
 !
 ! Check for presence of the restart file.
 can_restart=.FALSE.
 if ( BS_files%in_haydock_basename /= BSE_NOFILE) then
   restart_file = TRIM(BS_files%in_haydock_basename)//TRIM(tag_file)
   if (file_exists(restart_file) ) then
     can_restart=.TRUE.
     msg = " Restarting Haydock calculation from file: "//TRIM(restart_file)
     call wrtout(std_out,msg,"COLL")
     call wrtout(ab_out,msg,"COLL")
   else
     can_restart=.FALSE.
     call wrtout(ab_out," WARNING: cannot find restart file: "//TRIM(restart_file),"COLL")
   end if
 end if
 ABI_CHECK(.not.can_restart,"restart not yet implemented")

 ! Open the file and write basic dimensions and info.
 if (my_rank==master) then
   out_file = TRIM(BS_files%out_basename)//TRIM(tag_file)
   call open_haydock(out_file,haydock_file)
   haydock_file%hsize = hexc%hsize
   haydock_file%use_coupling = Bsp%use_coupling
   haydock_file%op = BSE_HAYD_IMEPS
   haydock_file%nq = nkets
   haydock_file%broad = Bsp%broad
   call write_dim_haydock(haydock_file)
 end if

 !
 ! Calculate green(w) for the different starting points.
 green=czero
 do iq=1,nkets
   ABI_MALLOC(ket0,(hexc%hsize))
   ket0=kets(:,iq)

   !
   niter_file=0
   if (can_restart) then
     call haydock_restart(BSp,restart_file,BSE_HAYD_IMEPS,iq,hexc%hsize,&
&      niter_file,aa_file,bb_file,phi_nm1_file,phi_n_file,comm)
   end if
   !
   ! For n>1, we have:
   !  1) a_n = <n|H|n>
   !  2) b_n = || H|n> - a_n|n> -b_{n-1}|n-1> ||
   !  3) |n+1> = [H|n> -a_n|n> -b_{n-1}|n-1>]/b_n
   !
   ! The sequences starts with |1> normalized to 1 and b_0 =0, therefore:
   !  a_1 = <1|H|1>
   !  b_1 = || H|1> - a_1|1> ||
   !  |2> = [H|1> - a_1|1>]/b_1
   !
   ABI_MALLOC(hphi_n,(hexc%hsize))
   ABI_MALLOC(hphi_nm1,(hexc%hsize))
   ABI_MALLOC(phi_nm1,(my_nt))
   ABI_MALLOC(phi_n,(my_nt))

   niter_max = niter_file + Bsp%niter
   ABI_MALLOC(aa,(niter_max))
   ABI_MALLOC(bb,(niter_max))
   aa=czero; bb=zero

   if (niter_file==0) then       ! Calculation from scratch.
     phi_nm1=ket0(my_t1:my_t2)   ! Select the slice treated by this node.
     norm = DZNRM2(hexc%hsize,ket0,1) ! Normalization
     phi_nm1=phi_nm1/norm

     call hexc_matmul_tda(hexc,hexc_i,phi_nm1,hphi_n)

     aa(1)=xdotc(my_nt,phi_nm1,1,hphi_n(my_t1:),1)
     call xmpi_sum(aa(1:1),comm,ierr)

     phi_n = hphi_n(my_t1:my_t2) - aa(1)*phi_nm1

     bb(1) = xdotc(my_nt,phi_n,1,phi_n,1)
     call xmpi_sum(bb(1:1),comm,ierr)
     bb(1) = SQRT(bb(1))

     phi_n = phi_n/bb(1)
     niter_done=1

   else ! Use the previous a and b.
     niter_done=niter_file
     aa(1:niter_done) = aa_file
     bb(1:niter_done) = bb_file
     phi_nm1=phi_nm1_file(my_t1:my_t2)   ! Select the slice treated by this node.
     phi_n  =phi_n_file  (my_t1:my_t2)
   end if

   if (can_restart) then
     ABI_FREE(aa_file)
     ABI_FREE(bb_file)
     ABI_FREE(phi_nm1_file)
     ABI_FREE(phi_n_file)
   end if

   ! Multiplicative factor (k-point sampling and unit cell volume)
   ! TODO be careful with the spin here
   ! TODO four_pi comes from the coulomb term 1/|q| is already included in the
   ! oscillators hence the present approach wont work if a cutoff interaction is used.
   nfact = -four_pi/(Cryst%ucvol*hexc%nbz)
   if (nsppol==1) nfact=two*nfact

   factor = nfact*(DZNRM2(hexc%hsize,ket0,1)**2)

   ! Which quantity should be checked for convergence?
   check = (/.TRUE.,.TRUE./)
   if (ABS(Bsp%haydock_tol(2)-one)<tol6) check = (/.TRUE. ,.FALSE./)
   if (ABS(Bsp%haydock_tol(2)-two)<tol6) check = (/.FALSE.,.TRUE./)

   ! Create new frequencies "mirror" in negative range to add
   ! their contributions. Can be improved by computing only once
   ! zero frequency, but loosing clearness
   n_all_omegas = 2*BSp%nomega

   ABI_MALLOC(all_omegas,(n_all_omegas))
   ! Put all omegas with frequency > 0 in table
   all_omegas(BSp%nomega+1:n_all_omegas) = BSp%omega
   ! Put all omegas with frequency < 0
   ! Warning, the broadening must be kept positive
   all_omegas(1:BSp%nomega) = -DBLE(BSp%omega(BSp%nomega:1:-1)) + j_dpc*AIMAG(BSp%omega(BSp%nomega:1:-1))

   ABI_MALLOC(green_temp,(n_all_omegas,nkets))

   call haydock_herm_algo(niter_done,niter_max,n_all_omegas,all_omegas,BSp%haydock_tol(1),check,&
&    my_t1,my_t2,factor,term_type,aa,bb,phi_nm1,phi_n,&
&    green_temp(:,iq),inn,is_converged,&
&    hexc, hexc_i, comm)

   ! Computing result from two ranges of frequencies
   ! The real part is added, the imaginary part is substracted
   green(:,iq) = green_temp(BSp%nomega+1:n_all_omegas,iq)+CONJG(green_temp(BSp%nomega:1:-1,iq))

   ABI_FREE(all_omegas)
   ABI_FREE(green_temp)
   !
   ! Save the a"s and the b"s for possible restarting.
   ! 1) Info on the Q.
   ! 2) Number of iterations performed.
   ! 3) do iter=1,niter_performed
   !      aa(iter),bb(iter)
   !    end do
   ! 4) |n-1>
   !    |n>
   !
   hphi_nm1 = czero
   hphi_nm1(my_t1:my_t2) = phi_nm1
   call xmpi_sum_master(hphi_nm1,master,comm,ierr)

   hphi_n = czero
   hphi_n(my_t1:my_t2) = phi_n
   call xmpi_sum_master(hphi_n,master,comm,ierr)

   if (my_rank==master) then
     ! Write data for restarting
     call write_haydock(haydock_file, hexc%hsize, Bsp%q(:,iq), aa, bb, hphi_n, hphi_nm1, MIN(inn,niter_max), factor)
   end if

   ABI_FREE(hphi_n)
   ABI_FREE(hphi_nm1)
   ABI_FREE(phi_nm1)
   ABI_FREE(phi_n)
   ABI_FREE(aa)
   ABI_FREE(bb)
   ABI_FREE(ket0)
 end do ! iq

 if (my_rank==master) call close_haydock(haydock_file)

 call xmpi_barrier(comm)

end subroutine haydock_herm
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_herm_algo
!! NAME
!! haydock_herm_algo
!!
!! FUNCTION
!!
!! INPUTS
!!  niter_done=Number of iterations already performed (0 if the run starts from scratch).
!!  niter_max=Max number of iterations. Always > niter_done
!!  nomega=Number of Frequency points for the evaluation of the matrix element.
!!  omega(nomega)=Frequency set (imaginary part is already included).
!!  tol_iter=Tolerance used to stop the algorithm.
!!  check(2)=Logical flags to specify where both the real and the imaginary part of the
!!    matrix elements of the Green functions have to be checked for convergence.
!!  hsize=Size of the blocks.
!!  my_t1,my_t2=Indices of the first and last column stored treated by this done.
!!  term_type=0 if no terminator is used, 1 otherwise.
!!  hmat(hsize,my_t1:my_t2)=The columns of the block.
!!  factor
!!  ntrans = Number of transitions
!!  corresp = mapping between coarse points and neighbours
!!  overlaps = overlaps of wavefunctions between dense k-point coarse neighbours and bands
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  green(nomega)=Output matrix elements.
!!  inn=Last iteration performed.
!!  is_converged=.TRUE. of the algorithm converged.
!!
!! SIDE EFFECTS
!!  phi_nm1(my_t2-my_t1+1), phi_n(my_t2-my_t1+1)
!!    input: vectors used to initialize the iteration
!!    output: the vectors obtained in the last iteration
!!  aa(niter_max) and bb(niter_max)
!!    if niter_done>0: aa(1:niter_done), bb(1:niter_done) store the coefficients of the previous run.
!!    when the routine returns aa(1:inn) and bb(1:inn) contain the matrix elements of the tridiagonal form.
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine haydock_herm_algo(niter_done,niter_max,nomega,omega,tol_iter,check,&
& my_t1,my_t2,factor,term_type,aa,bb,phi_nm1,phi_n,&
& green,inn,is_converged,&
& hexc, hexc_i, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: niter_max,niter_done,nomega
 integer,intent(in) :: my_t1,my_t2,term_type
 integer,intent(in) :: comm
 integer,intent(out) :: inn
 logical,intent(out) :: is_converged
 real(dp),intent(in) :: tol_iter
 complex(dpc),intent(in) :: factor
 type(hexc_t),intent(in) :: hexc
 type(hexc_interp_t),intent(in) :: hexc_i
!arrays
 real(dp),intent(inout) :: bb(niter_max)
 complex(dpc),intent(out) :: green(nomega)
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(inout) :: aa(niter_max)
 complex(dpc),intent(inout) :: phi_nm1(my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phi_n  (my_t2-my_t1+1)
 logical,intent(in) :: check(2)

!Local variables ------------------------------
!scalars
 integer :: ierr,my_nt,niter_min,nconv
 character(len=500) :: msg
 logical,parameter :: force_real=.TRUE.
!arrays
 real(dp) :: abs_err(nomega,2) !,rel_err(nomega,2)
 complex(dpc),allocatable :: oldg(:),newg(:)
 complex(dpc),allocatable :: phi_np1(:),hphi_n(:),cfact(:)
 logical :: test(2)

!************************************************************************

 ! The sequences starts with |1> normalized to 1 and b_0 =0, therefore:
 !  a_1 = <1|H|1>
 !  b_1 = || H|1> - a_1|1> ||
 !  |2> = [H|1> - a_1|1>]/b_1
 !
 ! For n>1 we have
 !  1) a_n = <n|H|n>
 !  2) b_n = || H|n> - a_n|n> -b_{n-1}|n-1> ||
 !  3) |n+1> = [H|n> -a_n|n> -b_{n-1}|n-1>]/b_n
 !
 my_nt = my_t2-my_t1+1

 ABI_MALLOC_OR_DIE(hphi_n,(hexc%hsize), ierr)

 ABI_MALLOC(phi_np1,(my_nt))

 ABI_MALLOC(oldg,(nomega))
 oldg=czero
 ABI_MALLOC(newg,(nomega))
 newg=czero
 ABI_MALLOC(cfact,(nomega))
 cfact=czero

 nconv=0
 do inn=niter_done+1,niter_max

   !YG2014
   call hexc_matmul_tda(hexc,hexc_i,phi_n,hphi_n)

   aa(inn) = xdotc(my_nt,phi_n,1,hphi_n(my_t1:),1)
   call xmpi_sum(aa(inn:inn),comm,ierr)
   if (force_real) aa(inn) = DBLE(aa(inn)) ! Matrix is Hermitian.

   ! |n+1> = H|n> - A(n)|n> - B(n-1)|n-1>
   phi_np1 = hphi_n(my_t1:my_t2) - aa(inn)*phi_n - bb(inn-1)*phi_nm1

   bb(inn) = xdotc(my_nt,phi_np1,1,phi_np1,1)
   call xmpi_sum(bb(inn),comm,ierr)
   bb(inn) = SQRT(bb(inn))

   phi_np1 = phi_np1/bb(inn)

   phi_nm1 = phi_n
   phi_n   = phi_np1

   write(msg,'(a,i0,a,3es12.4)')' Iteration number ',inn,', b_i RE(a_i) IM(a_i) ',bb(inn),REAL(aa(inn)),AIMAG(aa(inn))
   call wrtout(std_out,msg,"COLL")

   call continued_fract(inn,term_type,aa,bb,nomega,omega,cfact)

   newg= factor*cfact
   !
   ! Avoid spurious convergence.
   niter_min=4; if (niter_done>1) niter_min=niter_done+1
   if (inn>niter_min) then
     test=.TRUE.
     abs_err(:,1) = ABS(DBLE (newg-oldg))
     abs_err(:,2) = ABS(AIMAG(newg-oldg))
     !
     if (tol_iter>zero) then
       ! Test on the L1 norm.
       if (check(1)) test(1) = SUM(abs_err(:,1)) < tol_iter*SUM(ABS(DBLE (newg)))
       if (check(2)) test(2) = SUM(abs_err(:,2)) < tol_iter*SUM(ABS(AIMAG(newg)))
     else
       ! Stringent test for each point.
       if (check(1)) test(1) = ALL( abs_err(:,1) < -tol_iter*ABS(DBLE (newg)))
       if (check(2)) test(2) = ALL( abs_err(:,2) < -tol_iter*ABS(AIMAG(newg)))
     end if
     !
     if (ALL(test)) then
       nconv = nconv+1
     else
       nconv = 0
     end if
     if (nconv==2) then
       if(inn<100)then
         write(msg,'(a,es10.2,a)')&
&          " >>> Haydock algorithm converged twice within haydock_tol= ",tol_iter," after less than 100 iterations."
       else
         write(msg,'(a,es10.2,a,i0,a)')&
&          " >>> Haydock algorithm converged twice within haydock_tol= ",tol_iter," after ",inn," iterations."
       endif
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       EXIT
     end if
   end if

   oldg = newg
 end do ! inn

 green = newg
 if (nconv/=2) then
   write(msg,'(a,es10.2,a,i0,a)')&
&    " WARNING: Haydock algorithm did not converge within ",tol_iter," after ",niter_max," iterations."
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

 is_converged = (nconv==2)

 ABI_FREE(oldg)
 ABI_FREE(newg)
 ABI_FREE(cfact)
 ABI_FREE(hphi_n)
 ABI_FREE(phi_np1)

end subroutine haydock_herm_algo
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_restart
!! NAME
!! haydock_restart
!!
!! FUNCTION
!! Restart the Haydock method from file reading the data produced in a previous run.
!!
!! INPUTS
!!  BSp<type(excparam)>=Parameters defining the Bethe-Salpeter calculation.
!!    omega(BSp%nomega)=Frequency mesh for the macroscopic dielectric function (broadening is already included).
!!  iq_search=The index of the q-point to be searched.
!!  hsize
!!  comm=MPI communicator.
!!  nsppol
!!  restart_file
!!
!! OUTPUT
!!  niter_file=Number of iterations already performed. 0 to signal that an error occurred during the reading
!!  bb_file(:)
!!  aa_file(:)
!!  phi_n_file(:)
!!  phi_nm1_file(:)
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine haydock_restart(BSp,restart_file,ftype,iq_search,hsize,niter_file,aa_file,bb_file,phi_nm1_file,phi_n_file,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,hsize,iq_search,ftype
 integer,intent(out) :: niter_file
 character(len=*),intent(in) :: restart_file
 type(excparam),intent(in) :: BSp
!arrays
 real(dp),allocatable,intent(out) :: bb_file(:)
 complex(dpc),allocatable,intent(out) :: aa_file(:),phi_n_file(:),phi_nm1_file(:)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: nproc,my_rank,ierr,op_file
 integer :: hsize_file,use_coupling_file
 complex(dpc) :: factor_file
 character(len=500) :: msg
 type(haydock_type) :: haydock_file

!************************************************************************

 nproc = xmpi_comm_size(comm); my_rank= xmpi_comm_rank(comm)

 if (my_rank==master) then
   call open_haydock(restart_file, haydock_file)

   call read_dim_haydock(haydock_file)

   if (haydock_file%op/=ftype) then
     write(msg,"(2(a,i0))")" Expecting restart file with filetype: ",ftype," but found ",op_file
     MSG_ERROR(msg)
   end if

   if (haydock_file%hsize/=hsize) then
     write(msg,"(2(a,i0))")&
&      " Rank of H_exc read from file: ",hsize_file," differs from the one used in this run: ",hsize
     MSG_ERROR(msg)
   end if

   if (haydock_file%use_coupling /= BSp%use_coupling) then
     write(msg,'(2(a,i0))')&
&      " use_coupling_file: ",use_coupling_file," differs from input file value: ",BSp%use_coupling
     MSG_ERROR(msg)
   end if

   call read_haydock(haydock_file, Bsp%q(:,iq_search), aa_file, bb_file, &
&                   phi_n_file, phi_nm1_file, niter_file, factor_file)

   if (niter_file == 0) then
     write(msg,"(a,3f8.4,3a)")&
&      " Could not find q-point: ",BSp%q(:,iq_search)," in file ",TRIM(restart_file),&
&      " Cannot restart Haydock iterations for this q-point"
     MSG_COMMENT(msg)
   else
     write(msg,'(a,i0)')" Number of iterations already performed: ",niter_file
     call wrtout(std_out,msg,"COLL")
     call wrtout(ab_out,msg,"COLL")

     if ( ABS(haydock_file%broad - BSp%broad) > tol6) then
       write(msg,'(2a,2(a,f8.4),a)')&
&        " Restart file has been produced with a different Lorentzian broadening: ",ch10,&
&        " broad_file: ",haydock_file%broad," input broadening: ",BSp%broad," Continuing anyway. "
       MSG_WARNING(msg)
     end if

     call close_haydock(haydock_file)
   end if
 end if
 !
 ! Master broadcasts the data.
 call xmpi_bcast(niter_file,master,comm,ierr)

 if (my_rank/=master) then
   ABI_MALLOC(aa_file,(niter_file))
   ABI_MALLOC(bb_file,(niter_file))
   ABI_MALLOC(phi_nm1_file,(hsize))
   ABI_MALLOC(phi_n_file,(hsize))
 end if

 call xmpi_bcast(aa_file,master,comm,ierr)
 call xmpi_bcast(bb_file,master,comm,ierr)
 call xmpi_bcast(phi_nm1_file,master,comm,ierr)
 call xmpi_bcast(phi_n_file,master,comm,ierr)

end subroutine haydock_restart
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_mdf_to_tensor
!! NAME
!! haydock_mdf_to_tensor
!!
!! FUNCTION
!! Transform macroscopic dielectric function from green function to each components of the tensor in red and cart coord.
!!
!! INPUTS
!!  BSp<type(excparam)>=Parameters defining the Bethe-Salpeter calculation.
!!    omega(BSp%nomega)=Frequency mesh for the macroscopic dielectric function (broadening is already included).
!!  Cryst=Parameters of the crystal
!!  eps(BSp%nomega,BSp%nq) = Macroscopic dielectric function to be written.
!!
!! OUTPUT
!!  tensor_cart(BSp%nomega,6) = dielectric tensor for each frequency, order (11,22,33,12,13,23) in cart. coord.
!!  tensor_red(BSp%nomega, 6) = idem in reduced coordinated
!!  ierr = 0 if the tensors have been successfully computed
!!      \= 0 if the system is ill-posed in terms of q-points (not enough or not independent q-points)
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine haydock_mdf_to_tensor(BSp,Cryst,eps,tensor_cart,tensor_red,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: ierr
 type(excparam),intent(in) :: BSp
 type(crystal_t),intent(in) :: Cryst
!arrays
 complex(dpc),intent(in) :: eps(BSp%nomega,BSp%nq)
 complex(dpc),intent(out) :: tensor_cart(BSp%nomega,6), tensor_red(BSp%nomega,6)

!Local variables ------------------------------
!scalars
 integer :: iq,info
 real(dp) :: normqcart, normqred
!arrays
 integer,allocatable :: ipiv(:)
 real(dp) :: qcart(3), qtmet(3)
 real(dp) :: qred2cart(3,3),qcart2red(3,3)
 complex(dpc) :: qqcart(BSp%nq,6), qqred(BSp%nq,6)
 complex(dpc) :: b(6,BSP%nomega)

!************************************************************************

 ! Error flag
 ierr = 0

 if(BSp%nq /= 6) then
    ierr = -1
    return
 end if

 ! Transformation matrices from reduced coordinates to cartesian coordinates
 qred2cart = two_pi*Cryst%gprimd
 qcart2red = qred2cart
 call matrginv(qcart2red,3,3)
 do iq = 1, 6

   ! Computing cartesian q-vector
   qcart = MATMUL(qred2cart, BSp%q(:,iq))

   ! Computing product 'metric - qred' to form quadratic form
   qtmet = (two_pi**2)*MATMUL(Cryst%gmet, BSp%q(:,iq))

   ! squared norms
   normqcart = qcart(1)**2+qcart(2)**2+qcart(3)**2
   normqred = (normv(BSp%q(:,iq),Cryst%gmet,"G"))**2

   ! Compute line 'iq' for matrix in cartesian coord
   qqcart(iq,1) = (qcart(1))**2
   qqcart(iq,2) = (qcart(2))**2
   qqcart(iq,3) = (qcart(3))**2
   qqcart(iq,4) = 2*(qcart(1)*qcart(2))
   qqcart(iq,5) = 2*(qcart(1)*qcart(3))
   qqcart(iq,6) = 2*(qcart(2)*qcart(3))

   ! Compute line 'iq' for matrix in reduced coord
   qqred(iq,1) = (qtmet(1))**2
   qqred(iq,2) = (qtmet(2))**2
   qqred(iq,3) = (qtmet(3))**2
   qqred(iq,4) = 2*(qtmet(1)*qtmet(2))
   qqred(iq,5) = 2*(qtmet(1)*qtmet(3))
   qqred(iq,6) = 2*(qtmet(2)*qtmet(3))

   ! Renormalize line
   qqcart(iq,:) = qqcart(iq,:)/normqcart
   qqred(iq,:) = qqred(iq,:)/normqred
 end do

 ABI_MALLOC(ipiv,(6))

 ! Solving linear system
 b = TRANSPOSE(eps)
 call ZGESV(6,BSp%nomega,qqcart,6,ipiv,b,6,info)
 tensor_cart = TRANSPOSE(b)

 if(info /= 0) then
   ! Skipping the rest of the routine
   ierr = info
   ABI_FREE(ipiv)
   return
 end if

 b = TRANSPOSE(eps)
 call ZGESV(6,BSp%nomega,qqred,6,ipiv,b,6,info)
 tensor_red = TRANSPOSE(b)

 if(info /= 0) ierr = info

 ABI_FREE(ipiv)

end subroutine haydock_mdf_to_tensor
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_psherm
!! NAME
!! haydock_psherm
!!
!! FUNCTION
!!  Reads the excitonic Hamiltonian from file and construct the Lanczos set of vectors
!!  by iterative matrix-vector multiplications.
!!
!! INPUTS
!!  BSp<type(excparam)>=Parameters defining the Bethe-Salpeter calculation.
!!    omega(BSp%nomega)=Frequency mesh for the macroscopic dielectric function (broadening is already included).
!! hize
!! my_t1,my_t2
!! hreso(hsize,my_t1:my_t2)
!! hcoup(hsize,my_t1:my_t2)
!! nkets
!! kets(hsize,nkets)
!! comm=MPI communicator.
!!
!! OUTPUT
!!  green(BSp%nomega)=The imaginary part of the macroscopic dielectric function.
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine haydock_psherm(BSp,BS_files,Cryst,Hdr_bse,hexc,hexc_i,hsize,my_t1,my_t2,nkets,kets,green,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: hsize,my_t1,my_t2,nkets,comm
 type(crystal_t),intent(in) :: Cryst
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(Hdr_type),intent(in) :: Hdr_bse
!arrays
 complex(dp),intent(out) :: green(BSp%nomega,BSp%nq)
 complex(dpc),intent(in) :: kets(hsize,nkets)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: inn,itt,out_unt,nproc,my_rank,ierr
 integer :: niter_file,niter_max,niter_done,nsppol,iq,my_nt,term_type
 real(dp) :: ket0_hbar_norm,nfact
 logical :: can_restart,is_converged
 complex(dpc) :: factor
 character(len=fnlen),parameter :: tag_file="_HAYDC_SAVE"
 character(len=500) :: msg
 character(len=fnlen) :: restart_file,out_file
 type(hexc_t),intent(in) :: hexc
 type(hexc_interp_t),intent(in) :: hexc_i
!arrays
 real(dp),allocatable :: bb_file(:)
 real(dp),allocatable :: bb(:)
 complex(dpc),allocatable :: aa(:),cc(:),phi_np1(:),phi_n(:),phi_nm1(:),cbuff(:)
 complex(dpc),allocatable :: aa_file(:),phi_n_file(:),phi_np1_file(:),cc_file(:)
 complex(dpc),allocatable :: ket0(:)
 logical :: check(2)

!************************************************************************

 MSG_WARNING("Haydock + coupling is still under development")

 if(BSp%use_interp) then
   MSG_ERROR("Coupling is not yet implemented with interpolation")
 end if

 nproc  = xmpi_comm_size(comm)
 my_rank= xmpi_comm_rank(comm)
 nsppol = Hdr_bse%nsppol

 my_nt = my_t2-my_t1+1
 ABI_CHECK(my_nt>0,"One of the processors has zero columns")

 ! Multiplicative factor (k-point sampling and unit cell volume)
 ! TODO be careful with the spin here
 ! TODO four_pi comes from the coulomb term 1/|q| is already included in the
 ! oscillators hence the present approach wont work if a cutoff interaction is used.
 nfact = four_pi/(Cryst%ucvol*BSp%nkbz)
 if (nsppol==1) nfact=two*nfact

 write(msg,'(a,i0)')' Haydock algorithm with MAX number of iterations: ',BSp%niter
 call wrtout(std_out,msg,"COLL")
 !
 ! Check for presence of the restart file.
 can_restart=.FALSE.

 if ( BS_files%in_haydock_basename /= BSE_NOFILE) then
   restart_file = strcat(BS_files%in_haydock_basename,tag_file)
   if (file_exists(restart_file) ) then
     can_restart=.TRUE.
     msg = strcat(" Restarting Haydock calculation from file: ",restart_file)
     call wrtout(std_out,msg,"COLL")
     call wrtout(ab_out,msg,"COLL")
     MSG_ERROR("Restart is not tested")
   else
     can_restart=.FALSE.
     MSG_WARNING(strcat("Cannot find restart file: ",restart_file))
   end if
 end if
 !
 ! Open the file and writes basic dimensions and info.
 if (my_rank==master) then
   out_file = TRIM(BS_files%out_basename)//TRIM(tag_file)
   if (open_file(out_file,msg,newunit=out_unt,form="unformatted") /= 0) then
     MSG_ERROR(msg)
   end if
   ! write header TODO: standardize this part.
   write(out_unt)hsize,Bsp%use_coupling,BSE_HAYD_IMEPS,nkets,Bsp%broad
 end if
 !
 ! Select the terminator for the continued fraction.
 term_type=0 !; if (Bsp%hayd_term>0) term_type=2
 call wrtout(std_out,sjoin("Using terminator type: ",itoa(term_type)),"COLL")
 !
 ! Calculate green(w) for the different starting kets.
 green=czero
 do iq=1,nkets
   ABI_MALLOC(ket0,(my_nt))
   ket0 = kets(my_t1:my_t2,iq)
   !
   niter_file=0

   if (can_restart) then
     call haydock_restart(BSp,restart_file,BSE_HAYD_IMEPS,iq,hsize,&
&      niter_file,aa_file,bb_file,phi_np1_file,phi_n_file,comm)
   end if
   !
   ABI_MALLOC(phi_nm1,(my_nt))
   ABI_MALLOC(phi_n,(my_nt))
   ABI_MALLOC(phi_np1,(my_nt))
   !
   ! TODO: Note the different convention used for the coefficients
   ! Should use the same convention in the Hermitian case.
   niter_max = niter_file + Bsp%niter
   ABI_MALLOC(aa,(niter_max))
   ABI_MALLOC(bb,(niter_max+1))
   ABI_MALLOC(cc,(niter_max+1))
   aa=czero; bb=czero; cc=czero

   if (niter_file==0) then ! Calculation from scratch.
     phi_n   = ket0
     call hexc_matmul_full(hexc, hexc_i, phi_n, phi_np1, -1)
     !phi_np1 = MATMUL(hreso,ket0) - MATMUL(hcoup,CONJG(ket0))
     ket0_hbar_norm = SQRT(two*DBLE(DOT_PRODUCT(phi_n,phi_np1)))
     phi_n   = phi_n  /ket0_hbar_norm
     phi_np1 = phi_np1/ket0_hbar_norm
     !ket0    = ket0/ket0_hbar_norm
     cc(1)=zero ! <P|F|P>
     !cc(1) =  DOT_PRODUCT(ket0,phi_np1)
     write(std_out,*)" cc(1), ket0_hbar_norm =",cc(1),ket0_hbar_norm

     phi_nm1 = czero
     niter_done=0  ! TODO Be careful here

   else ! Use the previously calculates a and b.
     niter_done=niter_file
     MSG_ERROR("Restart not coded")
     !aa(1:niter_done) = aa_file
     !bb(1:niter_done) = bb_file
     !phi_np1=phi_np1_file(my_t1:my_t2)   ! Select the slice treated by this node.
     !phi_n  =phi_n_file  (my_t1:my_t2)
   end if

   if (can_restart) then
     ABI_FREE(aa_file)
     ABI_FREE(bb_file)
     ABI_FREE(cc_file)
     ABI_FREE(phi_np1_file)
     ABI_FREE(phi_n_file)
   end if

   ! This factor gives the correct results
   factor = -nfact*ket0_hbar_norm / SQRT(two)

   ! Which quantity should be checked for convergence?
   check = (/.TRUE.,.TRUE./)
   if (ABS(Bsp%haydock_tol(2)-one)<tol6) check = (/.TRUE. ,.FALSE./)
   if (ABS(Bsp%haydock_tol(2)-two)<tol6) check = (/.FALSE.,.TRUE./)

   call haydock_psherm_optalgo(niter_done,niter_max,BSp%nomega,BSp%omega,BSp%haydock_tol(1),check,hexc,hexc_i,&
&    hsize,my_t1,my_t2,factor,term_type,aa,bb,cc,ket0,ket0_hbar_norm,phi_nm1,phi_n,phi_np1,&
&    green(:,iq),inn,is_converged,comm)
   !
   ! Save the a"s and the b"s for possible restarting.
   ! 1) Info on the Q.
   ! 2) Number of iterations performed.
   ! 3) do iter=1,niter_performed
   !      aa(iter),bb(iter)
   !    end do
   ! 4) |n-1>
   !    |n>
   !    |n+1>
   !
   if (my_rank==master) then ! Open the file and writes basic dimensions and info.
     write(out_unt)Bsp%q(:,iq)
     write(out_unt)MIN(inn,niter_max)  ! NB: if the previous loop completed inn=niter_max+1
     do itt=1,MIN(inn,niter_max)        !     if we exited then inn is not incremented by one.
       write(out_unt)itt,aa(itt),bb(itt)
     end do
   end if
   !
   ! cbuff is used as workspace to gather |n-1>, |n> and |n+1>.
   ABI_MALLOC(cbuff,(hsize))
   cbuff=czero; cbuff(my_t1:my_t2) = phi_nm1
   call xmpi_sum_master(cbuff,master,comm,ierr)
   if (my_rank==master) write(out_unt) cbuff ! |n-1>

   cbuff=czero; cbuff(my_t1:my_t2) = phi_n
   call xmpi_sum_master(cbuff,master,comm,ierr)
   if (my_rank==master) write(out_unt) cbuff ! |n>

   cbuff=czero; cbuff(my_t1:my_t2) = phi_np1
   call xmpi_sum_master(cbuff,master,comm,ierr)
   if (my_rank==master) write(out_unt) cbuff ! |n+1>

   ABI_FREE(phi_nm1)
   ABI_FREE(phi_n)
   ABI_FREE(phi_np1)
   ABI_FREE(cbuff)
   ABI_FREE(aa)
   ABI_FREE(bb)
   ABI_FREE(cc)
   ABI_FREE(ket0)
 end do ! iq

 if (my_rank==master) close(out_unt)

 call xmpi_barrier(comm)

end subroutine haydock_psherm
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_psherm_optalgo
!! NAME
!! haydock_psherm_optalgo
!!
!! FUNCTION
!!  Haydock algorithm for pseudo-hermitian matrix
!!
!! INPUTS
!!  niter_done=Number of iterations already performed (0 if the run starts from scratch).
!!  niter_tot=Max number of iterations. Always > niter_done
!!  nomega=Number of Frequency points for the evaluation of the matrix element.
!!  omega(nomega)=Frequency set (imaginary part is already included).
!!  tol_iter=Tollerance used to stop the the algorithm.
!!  check(2)=Logical flags to specify where both the real and the imaginary part of the
!!    matrix elements of the Green functions have to be checked for convergence.
!!  hsize=Size of the blocks.
!!  my_t1,my_t2=Indeces of the first and last column stored treated by this done.
!!  term_type=0 if no terminator is used, 1 otherwise.
!!  hreso(hsize,my_t1:my_t2)=The columns of the resonant block.
!!  hcoup(hsize,my_t1:my_t2)=The columns of the coupling block.
!!  factor
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  green(nomega)=Output matrix elements.
!!  inn=Last iteration performed.
!!  is_converged=.TRUE. of the algorithm converged.
!!
!! SIDE EFFECTS
!!  phi_nm1(my_t2-my_t1+1), phi_n(my_t2-my_t1+1)
!!    input: vectors used to initialize the iteration
!!    output: the vectors obtained in the last iteration
!!  aa(niter_tot) and bb(niter_tot+1)
!!    if niter_done>0: aa(1:niter_done), bb(1:niter_done) store the coefficients of the previous run.
!!    when the routine returns aa(1:inn) and bb(1:inn) contain the matrix elements of the tridiagonal form.
!!  cc(niter_tot+1)
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine haydock_psherm_optalgo(niter_done,niter_tot,nomega,omega,tol_iter,check,hexc,hexc_i,hsize,my_t1,my_t2,&
&  factor,term_type,aa,bb,cc,ket0,ket0_hbar_norm,phi_nm1,phi_n,phi_np1,green,inn,is_converged,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: niter_tot,niter_done,nomega,comm,hsize,my_t1,my_t2,term_type
 integer,intent(out) :: inn
 logical,intent(out) :: is_converged
 real(dp),intent(in) :: tol_iter,ket0_hbar_norm
 complex(dpc),intent(in) :: factor
 type(hexc_t),intent(in) :: hexc
 type(hexc_interp_t),intent(in) :: hexc_i
!arrays
 real(dp),intent(inout) :: bb(niter_tot+1)
 complex(dpc),intent(out) :: green(nomega)
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(inout) :: aa(niter_tot),cc(niter_tot+1)
 complex(dpc),intent(in) :: ket0(my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phi_nm1(my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phi_n  (my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phi_np1(my_t2-my_t1+1)
 logical,intent(in) :: check(2)

!Local variables ------------------------------
!scalars
 integer :: my_nt,niter_min,nconv,parity,ii,jj,tdim,ierr
 integer :: row_max,col_max,nlev
 character(len=500) :: msg
 real(dp) :: max_err,mean_err,mean_err2,std_dev,err
 logical :: keep_vectors=.TRUE.
!arrays
 real(dp) :: abs_err(nomega,2) !,ww_err(nomega,2)
 complex(dpc) :: gn0(nomega,niter_tot)
 complex(dpc),allocatable :: oldg(:),newg(:)
 complex(dpc),allocatable :: hphi_n(:),save_phi(:,:)
 complex(dpc),allocatable ::  alpha(:,:),beta(:,:),ovlp(:,:)
 complex(dpc),allocatable :: phi_test(:),phi_test2(:),g00(:)
 logical :: test(2)

!************************************************************************

 ABI_UNUSED(ket0_hbar_norm)

 my_nt = my_t2-my_t1+1

 ABI_MALLOC(oldg,(nomega))
 ABI_MALLOC(newg,(nomega))
 ABI_MALLOC(g00,(nomega))
 oldg=czero; newg=czero; g00=czero
 nconv=0

 keep_vectors = (keep_vectors.and.xmpi_comm_size(comm)==1)
 if (keep_vectors) then
   ABI_MALLOC_OR_DIE(save_phi,(my_t2-my_t1+1,niter_tot), ierr)
   save_phi=czero
 end if

 ABI_MALLOC(hphi_n,(hsize))

 do inn=niter_done+1,niter_tot
   !
   ! a(n) = <Vn+1|F|Vn+1> = <Vn|HFH|Vn>) = 0 by symmetry.
   aa(inn)=zero

   ! |n+1> = |n+1> - a(n)|Vn> - a(n)|n-1>
   phi_np1 = phi_np1 - bb(inn)*phi_nm1
   !
   ! |n-1> = |n>
   ! |n>   = |n+1>
   phi_nm1 = phi_n
   phi_n   = phi_np1
   !
   !|n+1> = H |n> using all eh components.
   parity = (-1)**(inn+1)
   call hexc_matmul_full(hexc, hexc_i, phi_n, phi_np1, parity)

   !phi_np1 = MATMUL(hreso,phi_n) + parity * MATMUL(hcoup,CONJG(phi_n))
   !call xmpi_sum(hphi_np1,comm,ierr)
   !
   ! B(n+1)= <n|F|n+1>^(1/2) = <n|FH|n>^(1/2))= (2*Re(<n|V+1>))^(1/2)
   ! by symmetry, where the dot_product is done in the resonant eh sub-space.
   !
   bb(inn+1)=SQRT(two*DBLE(DOT_PRODUCT(phi_n,phi_np1)))
   !bb(inn+1)=two*DBLE(DOT_PRODUCT(phi_n,phi_np1))
   !call xmpi_sum(bb(inn+1),comm,ierr)
   !bb(inn+1)=SQRT(bb(inn+1)
   !
   !|n+1> =|n+1>/B(n+1)
   phi_n   = phi_n  /bb(inn+1)
   phi_np1 = phi_np1/bb(inn+1)

   if (keep_vectors) save_phi(:,inn) = phi_n

   parity = (-1)**(inn+1)
   !if (parity==-1) then
   !  cc(inn+1)=czero
   !else
     cc(inn+1)=DOT_PRODUCT(ket0,phi_n) + parity * DOT_PRODUCT(phi_n,ket0)
   !end if
   !call xmpi_sum(cc(inn+1),comm,ierr)

   write(msg,'(a,i0,a,3es12.4)')' Iteration number ',inn,', b_i RE(c_i+1) IM(c_i+1) ',bb(inn),REAL(cc(inn+1)),AIMAG(cc(inn+1))
   call wrtout(std_out,msg,"COLL")

   call continued_fract(inn,term_type,aa,bb(2:),nomega,omega,g00)
   gn0(:,1) = g00

   if (.FALSE.) then
     gn0(:,2) = (one - omega(:)*g00(:))/bb(2)
     do ii=3,inn
       gn0(:,ii) = -(-bb(ii)*gn0(:,ii-2) -omega(:)*gn0(:,ii-1))/bb(ii+1)
     end do
   else
     do ii=2,inn
       nlev = inn-ii
       call continued_fract(nlev,term_type,aa,bb(ii+1:),nomega,omega,g00)
       gn0(:,ii) = +bb(ii+1) * g00 * gn0(:,ii-1)
     end do
   end if

   newg=czero
   do ii=1,inn
     newg(:) = newg + cc(ii)* gn0(:,ii)
   end do
   newg = factor*newg
   !
   ! Avoid spurious convergence.
   niter_min=4; if (niter_done>1) niter_min=niter_done+1
   if (inn>niter_min) then
     test=.TRUE.
     abs_err(:,1) = ABS(DBLE (newg-oldg))
     abs_err(:,2) = ABS(AIMAG(newg-oldg))
     !
     if (tol_iter>zero) then
       ! Test on the L1 norm.
       if (check(1)) test(1) = SUM(abs_err(:,1)) < tol_iter*SUM(ABS(DBLE (newg)))
       if (check(2)) test(2) = SUM(abs_err(:,2)) < tol_iter*SUM(ABS(AIMAG(newg)))
     else
       ! Stringent test for each point.
       if (check(1)) test(1) = ALL( abs_err(:,1) < -tol_iter*ABS(DBLE (newg)))
       if (check(2)) test(2) = ALL( abs_err(:,2) < -tol_iter*ABS(AIMAG(newg)))
     end if
     !
     if (ALL(test)) then
       nconv = nconv+1
     else
       nconv = 0
     end if
     if (nconv==2) then
       if(inn<100)then
         write(msg,'(a,es10.2,a)')&
&          " >>> Haydock algorithm converged twice within haydock_tol= ",tol_iter," after less than 100 iterations."
       else
         write(msg,'(a,es10.2,a,i0,a)')&
&          " >>> Haydock algorithm converged twice within haydock_tol= ",tol_iter," after ",inn," iterations."
       endif
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       EXIT
     end if
   end if
   !
   oldg = newg
 end do ! inn

 green = newg
 if (nconv/=2) then
   write(msg,'(a,es10.2,a,i0,a)')&
&    " WARNING: Haydock algorithm did not converge within ",tol_iter," after ",niter_tot," iterations."
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

 is_converged = (nconv==2)

 ABI_FREE(oldg)
 ABI_FREE(newg)
 ABI_FREE(g00)
 ABI_FREE(hphi_n)

 if (keep_vectors) then
   tdim = MIN(inn,niter_tot)
   ABI_MALLOC(ovlp,(tdim,tdim))

   ABI_MALLOC(phi_test,(hsize))
   ABI_MALLOC(phi_test2,(hsize))

   max_err=smallest_real; mean_err=zero; mean_err2=zero; row_max=-1
   do ii=1,tdim
     parity = (-1)**(ii+1)
     phi_test  = save_phi(:,ii)
     call hexc_matmul_full(hexc, hexc_i, phi_test, phi_test2, parity)
     !phi_test2 = MATMUL(hreso,phi_test) + parity * MATMUL(hcoup,CONJG(phi_test))
     ovlp(ii,ii) = DOT_PRODUCT(phi_test,phi_test2) + DOT_PRODUCT(phi_test2,phi_test)
     err = ABS(ovlp(ii,ii)-cone)
     mean_err  = mean_err + err
     mean_err2 = mean_err2 + err**2
     if (err > max_err) then
       max_err = err
       row_max = ii
     end if
   end do
   mean_err = mean_err/tdim
   std_dev = mean_err2/tdim -mean_err**2
   write(std_out,'(a,i0,1x,3es14.6)')&
&   " Error in normalization (ii, max_err,mean,std_dev): ",row_max,max_err,mean_err,std_dev

   ABI_FREE(phi_test)
   ABI_FREE(phi_test2)

   ABI_MALLOC(alpha,(hsize,tdim))

   ! Less efficient but for sake of simplicity with hexc_matmul
   ! TODO possibility to call hreso * phi, and hcoup * phi separately
   do ii=1,tdim
     parity = (-1)**(ii+1)
     call hexc_matmul_full(hexc, hexc_i, save_phi(:,ii), alpha(:,ii), parity)
   end do

   !alpha = MATMUL(hreso,save_phi(:,1:tdim))
   !
   !do ii=1,tdim
   !  parity = (-1)**(ii+1)
   !  alpha(:,ii) =  alpha(:,ii) + parity*MATMUL(hcoup,CONJG(save_phi(:,ii)))
   !end do

   ovlp = MATMUL(TRANSPOSE(CONJG(save_phi(:,1:tdim))),alpha)

   ABI_MALLOC(beta,(hsize,tdim))
   do ii=1,tdim
     parity = (-1)**(ii+1)
     beta(:,ii)  =  parity*save_phi(:,ii)
     alpha(:,ii) = -parity*alpha(:,ii)
   end do

   ovlp = ovlp - MATMUL(TRANSPOSE(CONJG(beta)),alpha)

   max_err=smallest_real; row_max=-1; col_max=-1
   mean_err=zero; mean_err2=zero
   do jj=1,tdim
     do ii=1,jj
       err = ABS(ovlp(ii,jj))
       if (ii==jj) err = ABS(err - one)
       mean_err  = mean_err + err
       mean_err2 = mean_err2 + err**2
       if (err > max_err) then
         max_err = err
         row_max=ii
         col_max=jj
       end if
     end do
   end do

   mean_err = mean_err/(tdim*(tdim+1)/2)
   std_dev = mean_err2/(tdim*(tdim+1)/2) - mean_err**2
   write(std_out,'(a,2(i0,1x),3es14.6)')&
&     " Error in Hbar-ortho (i,j), max_err, mean, std_dev ",row_max,col_max,max_err,mean_err,std_dev
   !call print_arr(ovlp,max_r=185,max_c=10,unit=std_out)

   ABI_FREE(alpha)
   ABI_FREE(beta)
   ABI_FREE(ovlp)
   ABI_FREE(save_phi)
 end if

end subroutine haydock_psherm_optalgo
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_bilanczos
!! NAME
!! haydock_bilanczos
!!
!! FUNCTION
!!  Reads the excitonic Hamiltonian from file and construct the Lanczos set of vectors
!!  by iterative matrix-vector multiplications for any general matrix.
!!
!! INPUTS
!!  BSp<type(excparam)>=Parameters defining the Bethe-Salpeter calculation.
!!    omega(BSp%nomega)=Frequency mesh for the macroscopic dielectric function (broadening is already included).
!! hize
!! my_t1,my_t2
!! hreso(hsize,my_t1:my_t2)
!! hcoup(hsize,my_t1:my_t2)
!! nkets
!! kets(hsize,nkets)
!! comm=MPI communicator.
!!
!! OUTPUT
!!  green(BSp%nomega)=The imaginary part of the macroscopic dielectric function.
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine haydock_bilanczos(BSp,BS_files,Cryst,Hdr_bse,hexc,hexc_i,hsize,my_t1,my_t2,nkets,kets,ep_renorms,green,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: hsize,my_t1,my_t2,nkets,comm
 type(crystal_t),intent(in) :: Cryst
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(Hdr_type),intent(in) :: Hdr_bse
!arrays
 complex(dp),intent(out) :: green(BSp%nomega,BSp%nq)
 complex(dpc),intent(in) :: kets(hsize,nkets)
 complex(dpc),intent(in) :: ep_renorms(hsize)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: inn,itt,out_unt,nproc,my_rank,ierr
 integer :: niter_file,niter_max,niter_done,nsppol,iq,my_nt,term_type,n_all_omegas
 real(dp) :: ket0_hbar_norm,nfact,norm
 logical :: can_restart,is_converged
 complex(dpc) :: factor
 character(len=fnlen),parameter :: tag_file="_HAYDC_SAVE"
 character(len=500) :: msg
 character(len=fnlen) :: restart_file,out_file
 type(hexc_t),intent(in) :: hexc
 type(hexc_interp_t),intent(in) :: hexc_i
!arrays
 complex(dpc),allocatable :: aa_file(:),bb_file(:),cc_file(:)
 complex(dpc),allocatable :: aa(:),bb(:),cc(:)
 complex(dpc),allocatable :: phi_np1(:),phi_n(:),phi_nm1(:)
 complex(dpc),allocatable :: phit_np1(:),phit_n(:),phit_nm1(:)
 complex(dpc),allocatable :: cbuff(:), phi_n_file(:),phi_np1_file(:)
 complex(dpc),allocatable :: ket0(:)
 complex(dpc),allocatable :: hphi_n(:), hphit_n(:)
 complex(dpc),allocatable :: all_omegas(:),green_temp(:,:)
 logical :: check(2)

!************************************************************************

 MSG_WARNING("Haydock with Bilanczos is still under development")

 if(BSp%use_interp) then
   MSG_ERROR("Bilanczos is not yet implemented with interpolation")
 end if

 nproc  = xmpi_comm_size(comm)
 my_rank= xmpi_comm_rank(comm)
 nsppol = Hdr_bse%nsppol

 my_nt = my_t2-my_t1+1
 ABI_CHECK(my_nt>0,"One of the processors has zero columns")

 ! Multiplicative factor (k-point sampling and unit cell volume)
 ! TODO be careful with the spin here
 ! TODO four_pi comes from the coulomb term 1/|q| is already included in the
 ! oscillators hence the present approach wont work if a cutoff interaction is used.
 nfact = four_pi/(Cryst%ucvol*BSp%nkbz)
 if (nsppol==1) nfact=two*nfact

 write(msg,'(a,i0)')' Bi-Lanczos algorithm with MAX number of iterations: ',BSp%niter
 call wrtout(std_out,msg,"COLL")
 !
 ! Check for presence of the restart file.
 can_restart=.FALSE.

 if ( BS_files%in_haydock_basename /= BSE_NOFILE) then
   restart_file = strcat(BS_files%in_haydock_basename,tag_file)
   if (file_exists(restart_file) ) then
     can_restart=.TRUE.
     msg = strcat(" Restarting Haydock calculation from file: ",restart_file)
     call wrtout(std_out,msg,"COLL")
     call wrtout(ab_out,msg,"COLL")
     MSG_ERROR("Restart is not implemented")
   else
     can_restart=.FALSE.
     MSG_WARNING(strcat("Cannot find restart file: ",restart_file))
   end if
 end if
 !
 ! Open the file and writes basic dimensions and info.
 if (my_rank==master) then
   out_file = TRIM(BS_files%out_basename)//TRIM(tag_file)
   if (open_file(out_file,msg,newunit=out_unt,form="unformatted") /= 0) then
     MSG_ERROR(msg)
   end if
   ! write header TODO: standardize this part.
   write(out_unt)hsize,Bsp%use_coupling,BSE_HAYD_IMEPS,nkets,Bsp%broad
 end if
 !
 ! Select the terminator for the continued fraction.
 term_type=0 !; if (Bsp%hayd_term>0) term_type=2
 call wrtout(std_out,sjoin("Using terminator type: ",itoa(term_type)),"COLL")
 !
 ! Calculate green(w) for the different starting kets.
 green=czero
 do iq=1,nkets
   ABI_MALLOC(ket0,(hexc%hsize))
   ket0 = kets(:,iq)
   !
   niter_file=0

   if (can_restart) then
!     call haydock_restart(BSp,restart_file,BSE_HAYD_IMEPS,iq,hsize,&
!&      niter_file,aa_file,bb_file,phi_np1_file,phi_n_file,comm)
   end if
   !
   ABI_MALLOC(phi_nm1,(my_nt))
   ABI_MALLOC(phi_n,(my_nt))
   ABI_MALLOC(phi_np1,(my_nt))
   ABI_MALLOC(phit_nm1,(my_nt))
   ABI_MALLOC(phit_n,(my_nt))
   ABI_MALLOC(phit_np1,(my_nt))
   ABI_MALLOC(hphi_n,(hexc%hsize))
   ABI_MALLOC(hphit_n,(hexc%hsize))
   !
   ! TODO: Note the different convention used for the coefficients
   ! Should use the same convention in the Hermitian case.
   niter_max = niter_file + Bsp%niter
   ABI_MALLOC(aa,(niter_max))
   ABI_MALLOC(bb,(niter_max))
   ABI_MALLOC(cc,(niter_max))
   aa=czero; bb=czero; cc=czero

   if (niter_file==0) then ! Calculation from scratch.
     phi_nm1 = ket0(my_t1:my_t2)
     phit_nm1 = ket0(my_t1:my_t2)
     norm = DZNRM2(hexc%hsize,ket0,1)
     phi_nm1=phi_nm1/norm
     phit_nm1=phit_nm1/norm

     call hexc_matmul_elphon(hexc,phi_nm1,hphi_n,'N',ep_renorms)
     call hexc_matmul_elphon(hexc,phit_nm1,hphit_n,'C',ep_renorms)

     aa(1)=xdotc(my_nt,phit_nm1,1,hphi_n(my_t1:),1)
     call xmpi_sum(aa(1:1),comm,ierr)

     phi_n = hphi_n(my_t1:my_t2) - aa(1)*phi_nm1
     phit_n = hphit_n(my_t1:my_t2) - CONJG(aa(1))*phit_nm1

     bb(1)=xdotc(my_nt,phi_n,1,phi_n,1)
     call xmpi_sum(bb(1:1),comm,ierr)
     bb(1) = SQRT(bb(1))

     cc(1)=xdotc(my_nt,phit_n,1,phi_n,1)
     call xmpi_sum(cc(1:1),comm,ierr)
     cc(1) = cc(1)/bb(1)

     phi_n   = phi_n  /bb(1)
     phit_n  = phit_n /CONJG(cc(1))
     niter_done=1  ! TODO Be careful here

   else ! Use the previously calculates a and b.
     niter_done=niter_file
     MSG_ERROR("Restart not coded")
     !aa(1:niter_done) = aa_file
     !bb(1:niter_done) = bb_file
     !phi_np1=phi_np1_file(my_t1:my_t2)   ! Select the slice treated by this node.
     !phi_n  =phi_n_file  (my_t1:my_t2)
   end if

   if (can_restart) then
     ABI_FREE(aa_file)
     ABI_FREE(bb_file)
     ABI_FREE(cc_file)
     ABI_FREE(phi_np1_file)
     ABI_FREE(phi_n_file)
   end if

   ! This factor gives the correct results
   factor = -nfact*(DZNRM2(hexc%hsize,ket0,1)**2)

   ! Which quantity should be checked for convergence?
   check = (/.TRUE.,.TRUE./)
   if (ABS(Bsp%haydock_tol(2)-one)<tol6) check = (/.TRUE. ,.FALSE./)
   if (ABS(Bsp%haydock_tol(2)-two)<tol6) check = (/.FALSE.,.TRUE./)
   ! Create new frequencies "mirror" in negative range to add
   ! their contributions. Can be improved by computing only once
   ! zero frequency, but loosing clearness
   n_all_omegas = 2*BSp%nomega

   ABI_MALLOC(all_omegas,(n_all_omegas))
   ! Put all omegas with frequency > 0 in table
   all_omegas(BSp%nomega+1:n_all_omegas) = BSp%omega
   ! Put all omegas with frequency < 0
   ! Warning, the broadening must be kept positive
   all_omegas(1:BSp%nomega) = -DBLE(BSp%omega(BSp%nomega:1:-1)) + j_dpc*AIMAG(BSp%omega(BSp%nomega:1:-1))

   ABI_MALLOC(green_temp,(n_all_omegas,nkets))



   call haydock_bilanczos_optalgo(niter_done,niter_max,n_all_omegas,all_omegas,BSp%haydock_tol(1),check,hexc,hexc_i,&
&    hsize,my_t1,my_t2,factor,term_type,ep_renorms,aa,bb,cc,ket0,ket0_hbar_norm,phi_nm1,phi_n,phi_np1,&
&    phit_nm1,phit_n,phit_np1,green_temp(:,iq),inn,is_converged,comm)


   ! Computing result from two ranges of frequencies
   ! The real part is added, the imaginary part is substracted
   green(:,iq) = green_temp(BSp%nomega+1:n_all_omegas,iq)+CONJG(green_temp(BSp%nomega:1:-1,iq))

   ABI_FREE(all_omegas)
   ABI_FREE(green_temp)


   !
   ! Save the a"s and the b"s for possible restarting.
   ! 1) Info on the Q.
   ! 2) Number of iterations performed.
   ! 3) do iter=1,niter_performed
   !      aa(iter),bb(iter)
   !    end do
   ! 4) |n-1>
   !    |n>
   !    |n+1>
   !
   if (my_rank==master) then ! Open the file and writes basic dimensions and info.
     write(out_unt)Bsp%q(:,iq)
     write(out_unt)MIN(inn,niter_max)  ! NB: if the previous loop completed inn=niter_max+1
     do itt=1,MIN(inn,niter_max)        !     if we exited then inn is not incremented by one.
       write(out_unt)itt,aa(itt),bb(itt)
     end do
   end if
   !
   ! cbuff is used as workspace to gather |n-1>, |n> and |n+1>.
   ABI_MALLOC(cbuff,(hsize))
   cbuff=czero; cbuff(my_t1:my_t2) = phi_nm1
   call xmpi_sum_master(cbuff,master,comm,ierr)
   if (my_rank==master) write(out_unt) cbuff ! |n-1>

   cbuff=czero; cbuff(my_t1:my_t2) = phi_n
   call xmpi_sum_master(cbuff,master,comm,ierr)
   if (my_rank==master) write(out_unt) cbuff ! |n>

   cbuff=czero; cbuff(my_t1:my_t2) = phi_np1
   call xmpi_sum_master(cbuff,master,comm,ierr)
   if (my_rank==master) write(out_unt) cbuff ! |n+1>

   ABI_FREE(phi_nm1)
   ABI_FREE(phi_n)
   ABI_FREE(phi_np1)
   ABI_FREE(phit_nm1)
   ABI_FREE(phit_n)
   ABI_FREE(phit_np1)
   ABI_FREE(hphi_n)
   ABI_FREE(hphit_n)
   ABI_FREE(cbuff)
   ABI_FREE(aa)
   ABI_FREE(bb)
   ABI_FREE(cc)
   ABI_FREE(ket0)
 end do ! iq

 if (my_rank==master) close(out_unt)

 call xmpi_barrier(comm)

end subroutine haydock_bilanczos
!!***

!----------------------------------------------------------------------

!!****f* m_haydock/haydock_bilanczos_optalgo
!! NAME
!! haydock_bilanczos_optalgo
!!
!! FUNCTION
!!  Haydock algorithm for general matrix
!!
!! INPUTS
!!  niter_done=Number of iterations already performed (0 if the run starts from scratch).
!!  niter_tot=Max number of iterations. Always > niter_done
!!  nomega=Number of Frequency points for the evaluation of the matrix element.
!!  omega(nomega)=Frequency set (imaginary part is already included).
!!  tol_iter=Tollerance used to stop the the algorithm.
!!  check(2)=Logical flags to specify where both the real and the imaginary part of the
!!    matrix elements of the Green functions have to be checked for convergence.
!!  hsize=Size of the blocks.
!!  my_t1,my_t2=Indeces of the first and last column stored treated by this done.
!!  term_type=0 if no terminator is used, 1 otherwise.
!!  hreso(hsize,my_t1:my_t2)=The columns of the resonant block.
!!  hcoup(hsize,my_t1:my_t2)=The columns of the coupling block.
!!  factor
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  green(nomega)=Output matrix elements.
!!  inn=Last iteration performed.
!!  is_converged=.TRUE. of the algorithm converged.
!!
!! SIDE EFFECTS
!!  phi_nm1(my_t2-my_t1+1), phi_n(my_t2-my_t1+1)
!!    input: vectors used to initialize the iteration
!!    output: the vectors obtained in the last iteration
!!  aa(niter_tot) and bb(niter_tot+1)
!!    if niter_done>0: aa(1:niter_done), bb(1:niter_done) store the coefficients of the previous run.
!!    when the routine returns aa(1:inn) and bb(1:inn) contain the matrix elements of the tridiagonal form.
!!  cc(niter_tot+1)
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine haydock_bilanczos_optalgo(niter_done,niter_tot,nomega,omega,tol_iter,check,hexc,hexc_i,hsize,my_t1,my_t2,&
&  factor,term_type,ep_renorms,aa,bb,cc,ket0,ket0_hbar_norm,phi_nm1,phi_n,phi_np1,phit_nm1,phit_n,phit_np1,&
&  green,inn,is_converged,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: niter_tot,niter_done,nomega,comm,hsize,my_t1,my_t2,term_type
 integer,intent(out) :: inn
 logical,intent(out) :: is_converged
 real(dp),intent(in) :: tol_iter,ket0_hbar_norm
 complex(dpc),intent(in) :: factor
 type(hexc_t),intent(in) :: hexc
 type(hexc_interp_t),intent(in) :: hexc_i
!arrays
 complex(dpc),intent(inout) :: bb(niter_tot+1)
 complex(dpc),intent(out) :: green(nomega)
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(inout) :: aa(niter_tot),cc(niter_tot+1)
 complex(dpc),intent(in) :: ket0(my_t2-my_t1+1)
 complex(dpc),intent(in) :: ep_renorms(hsize)
 complex(dpc),intent(inout) :: phi_nm1(my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phi_n  (my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phi_np1(my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phit_nm1(my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phit_n  (my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phit_np1(my_t2-my_t1+1)
 logical,intent(in) :: check(2)

!Local variables ------------------------------
!scalars
 integer :: my_nt,niter_min,nconv !,ierr
 character(len=500) :: msg
 logical :: keep_vectors=.TRUE.
!arrays
 real(dp) :: abs_err(nomega,2) !,ww_err(nomega,2)
 complex(dpc),allocatable :: oldg(:),newg(:)
 complex(dpc),allocatable :: hphi_np1(:),hphit_np1(:),save_phi(:,:),save_phit(:,:)
 complex(dpc),allocatable :: g00(:)
 logical :: test(2)
 integer :: ierr

!************************************************************************

 ABI_UNUSED(ket0_hbar_norm)
 ABI_UNUSED(ket0(1))
 ABI_UNUSED(hexc_i%hsize_dense)

 my_nt = my_t2-my_t1+1

 ABI_MALLOC(oldg,(nomega))
 ABI_MALLOC(newg,(nomega))
 ABI_MALLOC(g00,(nomega))
 oldg=czero; newg=czero; g00=czero
 nconv=0

 keep_vectors = (keep_vectors.and.xmpi_comm_size(comm)==1)
 if (keep_vectors) then
   ABI_MALLOC(save_phi,(my_t2-my_t1+1,niter_tot))
   ABI_MALLOC_OR_DIE(save_phit,(my_t2-my_t1+1,niter_tot),ierr)
   save_phi=czero
   save_phit=czero
 end if

 ABI_MALLOC_OR_DIE(hphi_np1,(hexc%hsize),ierr)
 ABI_MALLOC_OR_DIE(hphit_np1,(hexc%hsize),ierr)

 do inn=niter_done+1,niter_tot

   !|n+1> = H |n> using all eh components.
   call hexc_matmul_elphon(hexc, phi_n, hphi_np1, 'N', ep_renorms)
   call hexc_matmul_elphon(hexc, phit_n, hphit_np1, 'C', ep_renorms)

   ! a(n) = < phit_n | H  | phi_n >
   aa(inn)=xdotc(my_nt,phit_n,1,hphi_np1(my_t1:),1)
   call xmpi_sum(aa(inn),comm,ierr)

   ! |n+1> = |n+1> - a(n)|Vn> - c(n)|n-1>
   phi_np1 = hphi_np1(my_t1:my_t2) - aa(inn)*phi_n - cc(inn-1)*phi_nm1
   phit_np1 = hphit_np1(my_t1:my_t2) - CONJG(aa(inn))*phit_n - CONJG(bb(inn-1))*phit_nm1

   bb(inn) = xdotc(my_nt,phi_np1,1,phi_np1,1)
   call xmpi_sum(bb(inn),comm,ierr)
   bb(inn) = SQRT(bb(inn))

   cc(inn) = xdotc(my_nt,phit_np1,1,phi_np1,1)
   call xmpi_sum(cc(inn),comm,ierr)
   cc(inn) = cc(inn)/bb(inn)

   phi_np1 = phi_np1 / bb(inn)
   phit_np1 = phit_np1 / CONJG(cc(inn))

   !
   ! |n-1> = |n>
   ! |n>   = |n+1>
   phi_nm1 = phi_n
   phi_n   = phi_np1
   phit_nm1 = phit_n
   phit_n = phit_np1

   if (keep_vectors) then
     save_phi(:,inn) = phi_n
     save_phit(:,inn) = phit_n
   end if
   write(msg,'(a,i0,a,3es12.4)')' Iteration number ',inn,', b_i RE(c_i) IM(c_i) ',REAL(bb(inn)),REAL(cc(inn)),AIMAG(cc(inn))
   call wrtout(std_out,msg,"COLL")

   call continued_fract_general(inn,term_type,aa,bb,cc,nomega,omega,g00)
   newg = factor*g00
   !
   ! Avoid spurious convergence.
   niter_min=4; if (niter_done>1) niter_min=niter_done+1
   if (inn>niter_min) then
     test=.TRUE.
     abs_err(:,1) = ABS(DBLE (newg-oldg))
     abs_err(:,2) = ABS(AIMAG(newg-oldg))
     !
     if (tol_iter>zero) then
       ! Test on the L1 norm.
       if (check(1)) test(1) = SUM(abs_err(:,1)) < tol_iter*SUM(ABS(DBLE (newg)))
       if (check(2)) test(2) = SUM(abs_err(:,2)) < tol_iter*SUM(ABS(AIMAG(newg)))
     else
       ! Stringent test for each point.
       if (check(1)) test(1) = ALL( abs_err(:,1) < -tol_iter*ABS(DBLE (newg)))
       if (check(2)) test(2) = ALL( abs_err(:,2) < -tol_iter*ABS(AIMAG(newg)))
     end if
     !
     if (ALL(test)) then
       nconv = nconv+1
     else
       nconv = 0
     end if
     if (nconv==2) then
       if(inn<100)then
         write(msg,'(a,es10.2,a)')&
&          " >>> Haydock algorithm converged twice within haydock_tol= ",tol_iter," after less than 100 iterations."
       else
         write(msg,'(a,es10.2,a,i0,a)')&
&          " >>> Haydock algorithm converged twice within haydock_tol= ",tol_iter," after ",inn," iterations."
       endif
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       EXIT
     end if
   end if
   !
   oldg = newg
 end do ! inn

 green = newg
 if (nconv/=2) then
   write(msg,'(a,es10.2,a,i0,a)')&
&    " WARNING: Haydock algorithm did not converge within ",tol_iter," after ",niter_tot," iterations."
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

 is_converged = (nconv==2)

 ABI_FREE(oldg)
 ABI_FREE(newg)
 ABI_FREE(g00)
 ABI_FREE(hphi_np1)
 ABI_FREE(hphit_np1)

 if (keep_vectors) then
   ABI_FREE(save_phi)
   ABI_FREE(save_phit)
 end if

 !! if (keep_vectors) then
 !!   tdim = MIN(inn,niter_tot)
 !!   ABI_MALLOC(ovlp,(tdim,tdim))

 !!   ABI_MALLOC(phi_test,(hsize))
 !!   ABI_MALLOC(phi_test2,(hsize))

 !!   max_err=smallest_real; mean_err=zero; mean_err2=zero; row_max=-1
 !!   do ii=1,tdim
 !!     parity = (-1)**(ii+1)
 !!     phi_test  = save_phi(:,ii)
 !!     call hexc_matmul_full(hexc, hexc_i, phi_test, phi_test2, parity)
 !!     !phi_test2 = MATMUL(hreso,phi_test) + parity * MATMUL(hcoup,CONJG(phi_test))
 !!     ovlp(ii,ii) = DOT_PRODUCT(phi_test,phi_test2) + DOT_PRODUCT(phi_test2,phi_test)
 !!     err = ABS(ovlp(ii,ii)-cone)
 !!     mean_err  = mean_err + err
 !!     mean_err2 = mean_err2 + err**2
 !!     if (err > max_err) then
 !!       max_err = err
 !!       row_max = ii
 !!     end if
 !!   end do
 !!   mean_err = mean_err/tdim
 !!   std_dev = mean_err2/tdim -mean_err**2
 !!   write(std_out,'(a,i0,1x,3es14.6)')&
 !!&    " Error in normalization (ii, max_err,mean,std_dev): ",row_max,max_err,mean_err,std_dev

 !!   ABI_FREE(phi_test)
 !!   ABI_FREE(phi_test2)
 !!
 !!   ABI_MALLOC(alpha,(hsize,tdim))

 !!   ! Less efficient but for sake of simplicity with hexc_matmul
 !!   ! TODO possibility to call hreso * phi, and hcoup * phi separately
 !!   do ii=1,tdim
 !!     parity = (-1)**(ii+1)
 !!     call hexc_matmul_full(hexc, hexc_i, save_phi(:,ii), alpha(:,ii), parity)
 !!   end do

 !!   !alpha = MATMUL(hreso,save_phi(:,1:tdim))
 !!   !
 !!   !do ii=1,tdim
 !!   !  parity = (-1)**(ii+1)
 !!   !  alpha(:,ii) =  alpha(:,ii) + parity*MATMUL(hcoup,CONJG(save_phi(:,ii)))
 !!   !end do

 !!   ovlp = MATMUL(TRANSPOSE(CONJG(save_phi(:,1:tdim))),alpha)

 !!   ABI_MALLOC(beta,(hsize,tdim))
 !!   do ii=1,tdim
 !!     parity = (-1)**(ii+1)
 !!     beta(:,ii)  =  parity*save_phi(:,ii)
 !!     alpha(:,ii) = -parity*alpha(:,ii)
 !!   end do

 !!   ovlp = ovlp - MATMUL(TRANSPOSE(CONJG(beta)),alpha)

 !!   max_err=smallest_real; row_max=-1; col_max=-1
 !!   mean_err=zero; mean_err2=zero
 !!   do jj=1,tdim
 !!     do ii=1,jj
 !!       err = ABS(ovlp(ii,jj))
 !!       if (ii==jj) err = ABS(err - one)
 !!       mean_err  = mean_err + err
 !!       mean_err2 = mean_err2 + err**2
 !!       if (err > max_err) then
 !!         max_err = err
 !!         row_max=ii
 !!         col_max=jj
 !!       end if
 !!     end do
 !!   end do

 !!   mean_err = mean_err/(tdim*(tdim+1)/2)
 !!   std_dev = mean_err2/(tdim*(tdim+1)/2) - mean_err**2
 !!   write(std_out,'(a,2(i0,1x),3es14.6)')&
 !!      " Error in Hbar-ortho (i,j), max_err, mean, std_dev ",row_max,col_max,max_err,mean_err,std_dev
 !!   !call print_arr(ovlp,max_r=185,max_c=10,unit=std_out)

 !!   ABI_FREE(alpha)
 !!   ABI_FREE(beta)
 !!   ABI_FREE(ovlp)
 !!   ABI_FREE(save_phi)
 !! end if

end subroutine haydock_bilanczos_optalgo
!!***


!----------------------------------------------------------------------

!!****f* m_numeric_tools/continued_fract_general
!! NAME
!!  continued_fract
!!
!! FUNCTION
!!  This routine calculates the continued fraction:
!!
!!                        1
!! f(z) =  _______________________________
!!           z - a1 -        b1^2
!!                   _____________________
!!                     z - a2 -    b2^2
!!                             ___________
!!                                z -a3 -    ........
!!
!! INPUTS
!!  nlev=Number of "levels" in the continued fraction.
!!  term_type=Type of the terminator.
!!    0 --> No terminator.
!!   -1 --> Assume constant coefficients for a_i and b_i for i>nlev with a_inf = a(nlev) and b_inf = b(nleb)
!!    1 --> Same as above but a_inf and b_inf are obtained by averaging over the nlev values.
!!  aa(nlev)=Set of a_i coefficients.
!!  bb(nlev)=Set of b_i coefficients.
!!  nz=Number of points on the z-mesh.
!!  zpts(nz)=z-mesh.
!!
!! OUTPUT
!!  spectrum(nz)=Contains f(z) on the input mesh.
!!
!! PARENTS
!!      bsepostproc,m_haydock
!!
!! CHILDREN
!!
!! SOURCE

subroutine continued_fract_general(nlev,term_type,aa,bb,cc,nz,zpts,spectrum)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nlev,term_type,nz
!arrays
 complex(dpc),intent(in) :: bb(nlev)
 complex(dpc),intent(in) :: cc(nlev)
 complex(dpc),intent(in) :: aa(nlev)
 complex(dpc),intent(in) :: zpts(nz)
 complex(dpc),intent(out) :: spectrum(nz)

!Local variables ------------------------------
!scalars
 integer :: it
 complex(dpc) ::  bb_inf,bg,bu,swap
 complex(dpc) :: aa_inf
 character(len=500) :: msg
!arrays
 complex(dpc),allocatable :: div(:),den(:)

!************************************************************************

 ABI_MALLOC(div,(nz))
 ABI_MALLOC(den,(nz))

 select case (term_type)
 case (0) ! No terminator.
   div=czero
 case (-1,1)
   MSG_ERROR("Not yet implemented")
   if (term_type==-1) then
     bb_inf=bb(nlev)
     aa_inf=aa(nlev)
   else
     bb_inf=SUM(bb)/nlev
     aa_inf=SUM(aa)/nlev
   end if
   ! Be careful with the sign of the SQRT.
   div(:) = half*(bb(nlev)/(bb_inf))**2 * ( zpts-aa_inf - SQRT((zpts-aa_inf)**2 - four*bb_inf**2) )
 case (2)
   MSG_ERROR("Not yet implemented")
   div = zero
   if (nlev>4) then
     bg=zero; bu=zero
     do it=1,nlev,2
       if (it+2<nlev) bg = bg + bb(it+2)
       bu = bu + bb(it)
     end do
     bg = bg/(nlev/2+MOD(nlev,2))
     bu = bg/((nlev+1)/2)
     !if (iseven(nlev)) then
     if (.not.iseven(nlev)) then
       swap = bg
       bg = bu
       bu = bg
     end if
     !write(std_out,*)nlev,bg,bu
     !Here be careful with the sign of SQRT
     do it=1,nz
       div(it) = half/zpts(it) * (bb(nlev)/bu)**2 * &
&        ( (zpts(it)**2 +bu**2 -bg**2) - SQRT( (zpts(it)**2+bu**2-bg**2)**2 -four*(zpts(it)*bu)**2) )
     end do
   end if

 case default
   write(msg,'(a,i0)')" Wrong value for term_type : ",term_type
   MSG_ERROR(msg)
 end select

 do it=nlev,2,-1
   den(:) = zpts(:) - aa(it) - div(:)
   div(:) = (bb(it-1)*cc(it-1) )/ den(:)
 end do

 den = zpts(:) - aa(1) - div(:)
 div = one/den(:)

 spectrum = div
 ABI_FREE(div)
 ABI_FREE(den)

end subroutine continued_fract_general
!!***

!----------------------------------------------------------------------

end module m_haydock
!!***
