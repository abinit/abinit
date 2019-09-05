!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_exc_analyze
!! NAME
!!  m_exc_analyze
!!
!! FUNCTION
!!  Postprocessing tools for BSE calculations.
!!
!! COPYRIGHT
!!  Copyright (C) 1992-2009 EXC group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida)
!!  Copyright (C) 2009-2019 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
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

module m_exc_analyze

 use defs_basis
 use m_abicore
 use m_bs_defs
 use m_xmpi
 use m_errors

 use defs_datatypes,      only : pseudopotential_type, ebands_t
 use m_io_tools,          only : open_file
 use m_numeric_tools,     only : iseven, wrap2_zero_one
 use m_bz_mesh,           only : kmesh_t, get_BZ_item
 use m_crystal,           only : crystal_t
 use m_wfd,               only : wfd_t
 use m_bse_io,            only : exc_read_eigen
 use m_pptools,           only : printxsf
 use m_pawrad,            only : pawrad_type
 use m_pawtab,            only : pawtab_type,pawtab_get_lsize
 use m_pawfgrtab,         only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free, pawfgrtab_print
 use m_pawcprj,           only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_paw_pwaves_lmn,    only : paw_pwaves_lmn_t, paw_pwaves_lmn_init, paw_pwaves_lmn_free
 use m_paw_nhat,          only : nhatgrid

 implicit none

 private
!!***

 public :: exc_den
 !public :: exc_plot      Plots the excitonic wavefunction in real space.

contains
!!***

!!****f* m_exc_analyze/m_exc_analyze
!! NAME
!!  exc_plot
!!
!! FUNCTION
!!  Plots the excitonic wavefunction in real space.
!!
!! INPUTS
!! Bsp<excparam>=Container storing the input variables of the run.
!! BS_files<excfiles>=filenames used in the Bethe-Salpeter part.
!! Kmesh<kmesh_t>=Info on the k-point sampling for wave functions.
!! Cryst<crystal_t>=Structure defining the crystalline structure.
!! KS_Bst<ebands_t>
!! Pawtab(Cryst%ntypat*usepaw)<pawtab_type>=PAW tabulated starting data
!! Pawrad(ntypat*usepaw)<type(pawrad_type)>=paw radial mesh and related data.
!! Psps <pseudopotential_type>=variables related to pseudopotentials.
!! Wfd<wfd_t>=Handler for the wavefunctions.
!! ngfftf(18)=Info on the dense FFT mesh used for plotting the excitonic wavefunctions.
!! nrcell(3)=Number of cell replicas (the code will enforce odd number so that the e or the
!!  h can be centered in the bix box.
!! eh_rcoord(3)=Reduced coordinated of the (electron|hole)
!! which_fixed= 1 to fix the electron at eh_red, 2 for the hole.
!! spin_opt=1 for resolving the up-component, 2 for the down component, 3 if spin resolution is not wanted.
!! paw_add_onsite=.TRUE. if the onsite contribution is taken into account.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      exc_read_eigen,nhatgrid,paw_pwaves_lmn_free,paw_pwaves_lmn_init
!!      pawcprj_alloc,pawcprj_free,pawfgrtab_free,pawfgrtab_init
!!      pawfgrtab_print,pawtab_get_lsize,printxsf,wfd_change_ngfft,wfd_sym_ur
!!      wrap2_zero_one,wrtout,xmpi_bcast,xmpi_sum_master
!!
!! SOURCE

subroutine exc_plot(Bsp,Bs_files,Wfd,Kmesh,Cryst,Psps,Pawtab,Pawrad,paw_add_onsite,spin_opt,which_fixed,eh_rcoord,nrcell,ngfftf)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: which_fixed,spin_opt
 logical,intent(in) :: paw_add_onsite
 type(excparam),intent(in) :: Bsp
 type(excfiles),intent(in) :: BS_files
 type(kmesh_t),intent(in) :: Kmesh
 type(crystal_t),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
 type(wfd_t),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfftf(18),nrcell(3)
 real(dp),intent(in) :: eh_rcoord(3)
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: cplex1=1
 integer :: nsppol,usepaw,nspinor,comm,cond,val
 integer :: optcut,optgr0,optgr1,optgr2,optrad
 integer :: ik_bz,ierr,my_rank !ik_ibz, istwf_k, isym_k,itim_k, my_nbbp, npw_k,
 integer :: spin,spin_start,spin_stop,reh
 integer :: rt_idx,art_idx,ii,iatom,nr_tot
 integer :: irc,ir1,ir2,ir3,wp1,wp2,wp3,wp_idx,eh_fft_idx,eh_rr,rr
 integer :: hsize,xsf_unt,ncells,nvec
 integer :: sc_natom,master
 real(dp) :: ene_rt,k_dot_r12
 complex(dpc) :: res_coeff,antres_coeff,eikr12
 character(len=500) :: msg
 character(len=fnlen) :: eig_fname,out_fname
!arrays
 integer :: nrcl(3),bbox(3)
 !integer :: bbp_distrb(Wfd%mband,Wfd%mband),got(Wfd%nproc)
 integer,allocatable :: rcl2fft(:),sc_typat(:),l_size_atm(:),vec_idx(:)
 real(dp),parameter :: origin0(3)=(/zero,zero,zero/)
 real(dp) :: k_bz(3),eh_red(3),eh_shift(3),r12(3),sc_rprimd(3,3),scart_shift(3)
 real(dp),allocatable :: rclred(:,:),exc_phi2(:),sc_xcart(:,:) !,sc_znucl(:)
 complex(gwpc),allocatable :: exc_phi(:)
 complex(gwpc),allocatable :: ur_v(:),ur_c(:)
 complex(dpc),allocatable :: vec_list(:,:)
 !logical :: bbp_mask(Wfd%mband,Wfd%mband)
 type(pawcprj_type),allocatable :: Cp_v(:,:),Cp_c(:,:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)

!************************************************************************

 MSG_WARNING("Exc plot is still under development")

 call wrtout(std_out," Calculating excitonic wavefunctions in real space","COLL")
 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")
 ABI_CHECK(Wfd%usepaw==0,"PAW not coded")
 ABI_CHECK(Wfd%nproc==1,"nproc>1 not supported")

 ! If needed, prepare FFT tables to have u(r) on the ngfftf mesh.
 if ( ANY(ngfftf(1:3) /= Wfd%ngfft(1:3)) ) then
   call wfd%change_ngfft(Cryst,Psps,ngfftf)
 end if

 comm    = Wfd%comm
 my_rank = Wfd%my_rank
 master  = Wfd%master

 nsppol  = Wfd%nsppol
 nspinor = Wfd%nspinor
 usepaw  = Wfd%usepaw
 !
 ! TODO recheck this part.
 ! Prepare tables needed for splining the onsite term on the FFT box.
 if (Wfd%usepaw==1.and.paw_add_onsite) then
   MSG_WARNING("Testing the calculation of AE PAW wavefunctions.")
   ! Use a local pawfgrtab to make sure we use the correction in the paw spheres
   ! the usual pawfgrtab uses r_shape which may not be the same as r_paw.
   call pawtab_get_lsize(Pawtab,l_size_atm,Cryst%natom,Cryst%typat)
   ABI_DT_MALLOC(Pawfgrtab,(Cryst%natom))
   call pawfgrtab_init(Pawfgrtab,cplex1,l_size_atm,Wfd%nspden,Cryst%typat)
   ABI_FREE(l_size_atm)

   optcut=1                     ! use rpaw to construct Pawfgrtab.
   optgr0=0; optgr1=0; optgr2=0 ! dont need gY terms locally.
   optrad=1                     ! do store r-R.

   call nhatgrid(Cryst%atindx1,Cryst%gmet,Cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,Cryst%ntypat,&
&    optcut,optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)
   !Pawfgrtab is ready to use

   if (Wfd%pawprtvol>0) then
     call pawfgrtab_print(Pawfgrtab,natom=Cryst%natom,unit=std_out,&
&                         prtvol=Wfd%pawprtvol,mode_paral="COLL")
   end if

   ABI_DT_MALLOC(Paw_onsite,(Cryst%natom))
   call paw_pwaves_lmn_init(Paw_onsite,Cryst%natom,Cryst%natom,Cryst%ntypat,Cryst%rprimd,Cryst%xcart,&
&                           Pawtab,Pawrad,Pawfgrtab)
 end if
 !
 ! Number of cells in the big box. Odd numbers are needed to place the e or the h in the middle.
 do ii=1,3
   nrcl(ii) = nrcell(ii)
   if (iseven(nrcell(ii))) then
     nrcl(ii) = nrcell(ii)+1
     write(msg,'(2(a,i0))')" Enforcing odd number of cell replicas ",nrcell(ii)," --> ",nrcl(ii)
     MSG_WARNING(msg)
   end if
 end do

 ncells = PRODUCT(nrcl)
 nr_tot = Wfd%nfftot * ncells  ! Total number of points in the big box.
 bbox = nrcl*ngfftf(1:3)
 !
 ! rcl2fft: The image of the point in the small box.
 ! rcl2red: The reduced coordinates of the point in the big box in terms of rprimd.
 ABI_MALLOC(rcl2fft,(nr_tot))
 ABI_MALLOC(rclred,(3,nr_tot))

 irc = 0
 do ir3=0,bbox(3)-1 ! Loop over the points in the big box.
   do ir2=0,bbox(2)-1
     do ir1=0,bbox(1)-1
       irc = 1+irc
       !
       wp1=MODULO(ir1,ngfftf(1)) ! The FFT index of the point wrapped into the original cell.
       wp2=MODULO(ir2,ngfftf(2))
       wp3=MODULO(ir3,ngfftf(3))
       wp_idx = 1 + wp1 + wp2*ngfftf(1) + wp3*ngfftf(1)*ngfftf(2)
       rcl2fft(irc)  = wp_idx
       rclred(1,irc) = DBLE(ir1)/ngfftf(1) ! Reduced coordinates in terms of the origina cell.
       rclred(2,irc) = DBLE(ir2)/ngfftf(2)
       rclred(3,irc) = DBLE(ir3)/ngfftf(3)
     end do
   end do
 end do
 !
 ! Wrap the point to [0,1[ where 1 is not included (tol12)
 call wrap2_zero_one(eh_rcoord(1),eh_red(1),eh_shift(1))
 call wrap2_zero_one(eh_rcoord(2),eh_red(2),eh_shift(2))
 call wrap2_zero_one(eh_rcoord(3),eh_red(3),eh_shift(3))
 write(std_out,*)"Initial Position of (e|h) in reduced coordinates:",eh_red," shift: ",eh_shift
 !
 ! Initial position on the FFT grid (the closest one)
 eh_red(:)=NINT(eh_red(:)*ngfftf(1:3))
 !
 ! Not translate the (electron|hole) in the center of the big box.
 do ii=1,3
   if (nrcl(ii)/=1) eh_red(ii) = eh_red(ii) + (nrcl(ii)/2) * ngfftf(ii)
 end do
 write(std_out,*)"Position of (e|h) in the center of the big box in FFT coordinates:",eh_red(:)

 eh_rr = 1 + eh_red(1) + eh_red(2)*bbox(1) + eh_red(3)*bbox(1)*bbox(2)
 eh_fft_idx = rcl2fft(eh_rr)

 write(std_out,*)" Reduced coordinates of (e|h) in the center of the big box ",rclred(:,eh_rr)
 !
 ! Allocate wavefunctions.
 ABI_MALLOC(ur_v,(Wfd%nfft*nspinor))
 ABI_MALLOC(ur_c,(Wfd%nfft*nspinor))

 if (usepaw==1) then
   ABI_DT_MALLOC(Cp_v,(Wfd%natom,nspinor))
   call pawcprj_alloc(Cp_v,0,Wfd%nlmn_atm)
   ABI_DT_MALLOC(Cp_c,(Wfd%natom,nspinor))
   call pawcprj_alloc(Cp_c,0,Wfd%nlmn_atm)
 end if
 !
 ! Read excitonic eigenvector.
 hsize = SUM(BSp%nreh); if (BSp%use_coupling>0) hsize=2*hsize
 nvec=1

 ABI_MALLOC_OR_DIE(vec_list,(hsize,nvec), ierr)

 ABI_MALLOC(vec_idx,(nvec))
 vec_idx = (/(ii, ii=1,nvec)/)

 if (BS_files%in_eig /= BSE_NOFILE) then
   eig_fname = BS_files%in_eig
 else
   eig_fname = BS_files%out_eig
 end if

 if (my_rank==master) then
   call exc_read_eigen(eig_fname,hsize,nvec,vec_idx,vec_list,Bsp=Bsp)
 end if

 call xmpi_bcast(vec_list,master,comm,ierr)

 ABI_FREE(vec_idx)
 !
 ! Allocate the excitonic wavefunction on the big box.
 ABI_MALLOC(exc_phi,(nr_tot))
 exc_phi=czero
 !
 !got=0;
 !bbp_mask=.FALSE.; bbp_mask(minb:maxb,minb:maxb)=.TRUE.
 spin_start=1; spin_stop=Wfd%nsppol
 if (spin_opt==1) spin_stop =1
 if (spin_opt==2) spin_start=2

 ! TODO
 ! Distribute (b,b',k,s) states where k is in the *full* BZ.
 !call wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,my_band_list,got,bmask)

 do spin=spin_start,spin_stop
   do reh=1,Bsp%nreh(spin)
     ene_rt = Bsp%Trans(reh,spin)%en
     ik_bz  = Bsp%Trans(reh,spin)%k
     val    = Bsp%Trans(reh,spin)%v
     cond   = Bsp%Trans(reh,spin)%c
     k_bz   = Kmesh%bz(:,ik_bz)
     !
     ! The resonant component.
     rt_idx = reh + (spin-1)*Bsp%nreh(1)
     res_coeff = vec_list(rt_idx,1)
     !
     ! The anti-resonant component.
     antres_coeff = czero
     if (Bsp%use_coupling>0) then
       art_idx  = reh + (spin-1)*Bsp%nreh(1) + SUM(Bsp%nreh)
       antres_coeff = vec_list(art_idx,1)
     end if

     call wfd%sym_ur(Cryst,Kmesh,val ,ik_bz,spin,ur_v) ! TODO recheck this one.
     call wfd%sym_ur(Cryst,Kmesh,cond,ik_bz,spin,ur_c)

     !call wfd_paw_get_aeur(Wfd,band,ik_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur_ae,ur_ae_onsite,ur_ps_onsite)

     if (which_fixed==1) then ! electron
       do rr=1,nr_tot
         wp_idx = rcl2fft(rr)
         r12 = eh_red - rclred(:,rr)
         k_dot_r12 = two_pi * DOT_PRODUCT(k_bz,r12)
         eikr12 = DCMPLX(COS(k_dot_r12), SIN(k_dot_r12))
         exc_phi(rr) = exc_phi(rr) + eikr12 * ur_v(wp_idx) * CONJG(ur_c(wp_idx))
       end do
     else ! hole
       MSG_ERROR("Not coded")
     end if
     !
  end do
 end do

 ABI_FREE(vec_list)
 ABI_FREE(rcl2fft)
 ABI_FREE(rclred)
 ABI_FREE(ur_v)
 ABI_FREE(ur_c)
 !
 if (Wfd%usepaw==1) then
   call pawcprj_free(Cp_v)
   ABI_DT_FREE(Cp_v)
   call pawcprj_free(Cp_c)
   ABI_DT_FREE(Cp_c)
   if (paw_add_onsite) then
     call pawfgrtab_free(Pawfgrtab)
     ABI_DT_FREE(Pawfgrtab)
     call paw_pwaves_lmn_free(Paw_onsite)
     ABI_DT_FREE(Paw_onsite)
   end if
 end if
 !
 ! Collect results on master node.
 call xmpi_sum_master(exc_phi,master,comm,ierr)
 !
 ! The exciton density probability.
 ABI_MALLOC(exc_phi2,(nr_tot))

 do rr=1,nr_tot
   exc_phi2(rr) = ABS(exc_phi(rr))**2
 end do
 !
 ! Construct the supercell.
 sc_natom = ncells*Cryst%natom
 ABI_MALLOC(sc_xcart,(3,sc_natom))
 ABI_MALLOC(sc_typat,(sc_natom))

 ii=0
 do ir3=0,nrcl(3)-1  ! Loop over the replicaed cells.
   do ir2=0,nrcl(2)-1
     do ir1=0,nrcl(1)-1
      !
      sc_rprimd(:,1) = ir1 * Cryst%rprimd(:,1) ! Workspace.
      sc_rprimd(:,2) = ir2 * Cryst%rprimd(:,2)
      sc_rprimd(:,3) = ir3 * Cryst%rprimd(:,3)
      scart_shift = SUM(sc_rprimd,DIM=2)
      do iatom=1,Cryst%natom
        ii=ii+1
        sc_xcart(:,ii) = Cryst%xcart(:,iatom) + scart_shift
        sc_typat(ii)   = Cryst%typat(iatom)
      end do
      !
    end do
   end do
 end do
 !
 ! Lattice vectors of the big box.
 do ii=1,3
   sc_rprimd(:,ii) = nrcl(ii) * Cryst%rprimd(:,ii)
 end do
 !
 ! Master writes the data.
 if (my_rank==master) then
   out_fname = TRIM(BS_files%out_basename)//"_EXC_WF"
   if (open_file(out_fname,msg,newunit=xsf_unt,form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if

   call printxsf(bbox(1),bbox(2),bbox(3),exc_phi2,sc_rprimd,origin0,sc_natom,Cryst%ntypat,sc_typat,&
&    sc_xcart,Cryst%znucl,xsf_unt,0)

   close(xsf_unt)
 end if

 ABI_FREE(sc_xcart)
 ABI_FREE(sc_typat)

 ABI_FREE(exc_phi)
 ABI_FREE(exc_phi2)

end subroutine exc_plot
!!***

!!****f* m_exc_analyze/exc_den
!! NAME
!!  exc_den
!!
!! FUNCTION
!!  This routines calculates the electron-hole excited state density.
!!
!! INPUTS
!!  filbseig=Name of the file containing the excitonic eigenvectors and eigenvalues.
!!
!! OUTPUT
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      wfd_get_ur,wrtout
!!
!! SOURCE

subroutine exc_den(BSp,BS_files,ngfft,nfftot,Kmesh,ktabr,Wfd)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(kmesh_t),intent(in) :: Kmesh
 type(wfd_t),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: ktabr(nfftot,BSp%nkbz)

!Local variables ------------------------------
!scalars
 integer :: nfft1,nfft2,nfft3,spin,ierr
 integer :: it,itp,ik_ibz,ikp_bz,ik_bz,band,iv,ivp,ic,icp,ir,i1,i2,i3,iik
 integer :: state,nstates,min_idx,eig_unt,den_unt,sden_unt,hsize,homo
 real(dp) :: min_ene,cost
 character(len=500) :: msg
 character(len=fnlen) :: filbseig
!arrays
 real(dp) :: n0(nfftot),rho_eh(nfftot),nexc(nfftot)
 real(dp),allocatable :: exc_ene(:)
 complex(dpc) :: rhor_h(nfftot),rhor_e(nfftot)
 complex(gwpc),allocatable :: wfr(:,:,:), wfrk(:,:)
 complex(dpc),allocatable :: exc_ene_cplx(:),exc_vec(:)

!************************************************************************

 ABI_CHECK(Wfd%nsppol==1,"nsppol==2 not coded")

 if (BS_files%in_eig /= BSE_NOFILE) then
   filbseig = BS_files%in_eig
 else
   filbseig = BS_files%out_eig
 end if

 call wrtout(std_out,' Calculating electron-hole, excited state density',"COLL")
 call wrtout(std_out," Reading eigenstates from: "//TRIM(filbseig),"COLL")
 MSG_ERROR("Not tested")

 ! Prepare FFT tables to have u(r) on the ngfft_osc mesh.
 !if ( ANY(ngfftf(1:3) /= Wfd%ngfft(1:3)) ) then
 !  call wfd%change_ngfft(Cryst,Psps,ngfftf)
 !end if

 nfft1 = ngfft(1)
 nfft2 = ngfft(2)
 nfft3 = ngfft(3)
 ABI_CHECK(nfftot==PRODUCT(ngfft(1:3)),"Mismatch in FFT size")

!allocate and load wavefunctions in real space
 ABI_MALLOC_OR_DIE(wfr,(nfftot*Wfd%nspinor,BSp%nbnds,Wfd%nkibz), ierr)

 spin=1 ! SPIN support is missing.

 do ik_ibz=1,Wfd%nkibz
   do band=1,BSp%nbnds
     call wfd%get_ur(band,ik_ibz,spin,wfr(:,band,ik_ibz))
   end do
 end do

 if (open_file(filbseig,msg,newunit=eig_unt,form="unformatted", status="old", action="read") /= 0) then
   MSG_ERROR(msg)
 end if

 read(eig_unt) hsize,nstates

 if (nstates/=Bsp%nreh(1)) then
   write(std_out,*) 'not resonant calculation'
   close(eig_unt)
   RETURN
 end if

 ABI_MALLOC(exc_ene_cplx,(nstates))
 ABI_MALLOC(exc_ene,(nstates))

 read(eig_unt) exc_ene_cplx(1:nstates)

 ! Small imaginary part is always neglected.
 exc_ene(:) = exc_ene_cplx(:)
 ABI_FREE(exc_ene_cplx)
 !
 ! Find the lowest non-negative eigenvalue.
 min_ene = greatest_real
 do state=1,nstates
   if (exc_ene(state) < zero) cycle
   if (exc_ene(state) < min_ene) then
     min_ene = exc_ene(state)
     min_idx = state
   end if
 end do
 ABI_FREE(exc_ene)

 write(std_out,*)' Lowest eigenvalue ', min_idx, min_ene*Ha_eV
 !
 ! skip other eigenvectors
 do state=1,min_idx-1
   read(eig_unt)
 end do
 !
 ! read "lowest" eigenvector
 ABI_MALLOC(exc_vec,(hsize))
 read(eig_unt) exc_vec

 close(eig_unt)

 ABI_MALLOC_OR_DIE(wfrk,(nfftot,BSp%nbnds), ierr)

!calculate ground state density
 n0(:) = zero
 spin = 1
 homo = Bsp%homo_spin(spin)
 do ik_bz = 1, BSp%nkbz
   ik_ibz = Kmesh%tab(ik_bz)
   iik = (3-Kmesh%tabi(ik_bz))/2

   if (iik==1) then
     do ir=1,nfftot
       wfrk(ir,1:homo) = wfr(ktabr(ir,ik_bz),1:homo,ik_ibz)
     end do
   else
     do ir=1,nfftot
       wfrk(ir,1:homo) = conjg(wfr(ktabr(ir,ik_bz),1:homo,ik_ibz))
     end do
   end if

   do iv=1,homo
     n0(:) = n0(:) + conjg(wfrk(:,band)) * wfrk(:,band)
   end do
 end do
 !
 ! calculate electron and hole density
 rhor_h=czero
 rhor_e=czero
 ABI_CHECK(BSp%nsppol==1,"nsppol=2 not coded")

 do it=1,Bsp%nreh(1)
   ik_bz = Bsp%Trans(it,1)%k
   iv    = Bsp%Trans(it,1)%v
   ic    = Bsp%Trans(it,1)%c
   ik_ibz = Kmesh%tab(ik_bz)
   iik = (3-Kmesh%tabi(ik_bz))/2

   if (iik==1) then
     do ir = 1, nfftot
       wfrk(ir,:) = wfr(ktabr(ir,ik_bz),:,ik_ibz)
     end do
   else
     do ir = 1, nfftot
       wfrk(ir,:) = conjg(wfr(ktabr(ir,ik_bz),:,ik_ibz))
     end do
   end if
   !
   do itp=1,Bsp%nreh(1)
     ikp_bz = Bsp%Trans(itp,1)%k
     if (ikp_bz/=ik_bz) CYCLE
     icp = Bsp%Trans(itp,1)%c
     ivp = Bsp%Trans(itp,1)%v
     if(icp==ic) then
       rhor_h = rhor_h - conjg(exc_vec(it)) * exc_vec(itp) * wfrk(:,iv) * conjg(wfrk(:,ivp))
     end if
     if(ivp==iv) then
       rhor_e = rhor_e + conjg(exc_vec(it)) * exc_vec(itp) * conjg(wfrk(:,ic)) * wfrk(:,icp)
     end if
   end do
 end do

 ABI_FREE(exc_vec)
 !
 ! calculate excited state density minus ground state density
 ! n* - n0 = rhor_e + rhor_h
 rho_eh(:) = rhor_e + rhor_h
 !
 !calculate excited state density
 ! n* = n0 + rhor_e + rhor_h
 nexc(:) = n0(:) + rho_eh(:)

 if (open_file('out.den',msg,newunit=den_unt) /= 0) then
   MSG_ERROR(msg)
 end if

 ! here conversion to cartesian through a1, a2, a3
 cost = zero
 do i1=0,nfft1-1
   do i2=0,nfft2-1
     do i3=0,nfft3-1
       ir=1 + i1 + i2*nfft1 + i3*nfft1*nfft2
       write(den_unt,'(3i3,2x,5e11.4)')i1,i2,i3,n0(ir),nexc(ir),rho_eh(ir),real(rhor_e(ir)),real(rhor_h(ir))
       cost = cost + n0(ir)
     end do
   end do
 end do

 write(std_out,*) 'density normalization constant ', cost
 close(den_unt)

 if (open_file('out.sden',msg,newunit=sden_unt) /= 0) then
   MSG_ERROR(msg)
 end if

 ! we are looking for the plane between (100) (111)
 ! so it's the place where (111) [ (100) x v ] = 0 (mixed product)
 ! finally v2 - v3 = 0
 do i2=0, nfft2-1
   do i3=0, nfft3-1
     if(i2==i3) then
       write(std_out,*) i2, i3
       write(sden_unt,*) (rho_eh(1+i1+i2*nfft1+i3*nfft1*nfft2),i1=0,nfft1-1)
     end if
   end do
 end do

 close(sden_unt)

end subroutine exc_den
!!***

end module m_exc_analyze
!!***
