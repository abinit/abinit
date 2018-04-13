!{\src2tex{textfont=tt}}
!!****f* ABINIT/exc_build_ham
!! NAME
!!  exc_build_ham
!!
!! FUNCTION
!!  Calculate and write the excitonic Hamiltonian on an external binary file (Fortran file open
!!  in random mode) for subsequent treatment in the Bethe-Salpeter code.
!!
!! COPYRIGHT
!! Copyright (C) 1992-2009 EXC group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida)
!! Copyright (C) 2009-2018 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  BSp<excparam>=The parameters for the Bethe-Salpeter calculation.
!!  BS_files<excfiles>=File names internally used in the BS code.
!!  Cryst<crystal_t>=Info on the crystalline structure.
!!  Kmesh<kmesh_t>=The list of k-points in the BZ, IBZ and symmetry tables.
!!  Qmesh<kmesh_t>=The list of q-points for epsilon^{-1} and related symmetry tables.
!!  ktabr(nfftot_osc,BSp%nkbz)=The FFT index of $(R^{-1}(r-\tau))$ where R is symmetry needed to obtains
!!    the k-points from the irreducible image.  Used to symmetrize u_Sk where S = \transpose R^{-1}
!!  Gsph_x<gsphere_t>=Info on the G-sphere used to describe wavefunctions and W (the largest one is actually stored).
!!  Gsph_c<gsphere_t>=Info on the G-sphere used to describe the correlation part.
!!  Vcp<vcoul_t>=The Coulomb interaction in reciprocal space. A cutoff can be used
!!  W<screen_t>=Data type gathering info and data for W.
!!  nfftot_osc=Total Number of FFT points used for the oscillator matrix elements.
!!  ngfft_osc(18)=Info on the FFT algorithm used to calculate the oscillator matrix elements.
!!  Psps<Pseudopotential_type>=Variables related to pseudopotentials
!!  Pawtab(Psps%ntypat)<pawtab_type>=PAW tabulated starting data.
!!  Pawang<pawang_type>=PAW angular mesh and related data.
!!  Paw_pwff(Cryst%ntypat*Wfd%usepaw)<pawpwff_t>=Form factor used to calculate the onsite mat. elements of a plane wave.
!!  Wfd<wfd_t>=Handler for the wavefunctions.
!!
!! OUTPUT
!!  The excitonic Hamiltonian is saved on an external binary file (see below).
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      cwtime,get_bz_item,gsph_fft_tabs,paw_rho_tw_g,pawcprj_alloc
!!      pawcprj_free,pawpwij_free,pawpwij_init,rho_tw_g,timab,wfd_change_ngfft
!!      wfd_get_cprj,wfd_get_ur,wrtout,xmpi_distab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine exc_build_ham(BSp,BS_files,Cryst,Kmesh,Qmesh,ktabr,Gsph_x,Gsph_c,Vcp,&
& Wfd,W,Hdr_bse,nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_bs_defs
 use m_xmpi
 use m_errors

 use m_time,         only : timab
 use m_gwdefs,       only : czero_gw
 use m_crystal,      only : crystal_t
 use m_gsphere,      only : gsphere_t
 use m_vcoul,        only : vcoul_t
 use m_bz_mesh,      only : kmesh_t
 use m_pawpwij,      only : pawpwff_t
 use m_pawang,       only : pawang_type
 use m_pawtab,       only : pawtab_type
 use m_pawcprj,      only : pawcprj_type
 use m_wfd,          only : wfd_t
 use m_screen,       only : screen_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_build_ham'
 use interfaces_14_hidewrite
 use interfaces_71_bse, except_this_one => exc_build_ham
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot_osc
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(screen_t),intent(inout) :: W
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(crystal_t),intent(in) :: Cryst
 type(vcoul_t),intent(in) :: Vcp
 type(gsphere_t),intent(in) :: Gsph_x,Gsph_c
 type(Pseudopotential_type),intent(in) :: Psps
 type(Hdr_type),intent(inout) :: Hdr_bse
 type(pawang_type),intent(in) :: Pawang
 type(wfd_t),target,intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfft_osc(18)
 integer,intent(in) :: ktabr(nfftot_osc,Kmesh%nbz)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Wfd%usepaw)
 type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 logical :: do_resonant,do_coupling
 !character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 complex(gwpc),allocatable :: all_mgq0(:,:,:,:,:)

!************************************************************************

 call timab(670,1,tsec)

 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")
 ABI_CHECK(nfftot_osc==PRODUCT(ngfft_osc(1:3)),"mismatch in FFT size")

 if (BSp%have_complex_ene) then
   MSG_ERROR("Complex energies are not supported yet")
 end if

 ! Do we have to compute some block?
 do_resonant = (BS_files%in_hreso == BSE_NOFILE)
 do_coupling = (BS_files%in_hcoup == BSE_NOFILE)

 if (BSp%use_coupling == 0) then
   if (.not.do_resonant) then
     call wrtout(std_out,"Will skip the calculation of resonant block (will use BSR file)","COLL")
     goto 100
   end if
 else
   if (.not. do_resonant .and. .not. do_coupling) then
     call wrtout(std_out,"Will skip the calculation of both resonant and coupling block (will use BSR and BSC files)","COLL")
     goto 100
   end if
 end if

 ! Compute M_{k,q=0}^{b,b}(G) for all k-points in the IBZ and each pair b, b'
 ! used for the exchange part and part of the Coulomb term.
 call wrtout(std_out," Calculating all matrix elements for q=0 to save CPU time","COLL")

 call wfd_all_mgq0(Wfd,Cryst,Qmesh,Gsph_x,Vcp,Psps,Pawtab,Paw_pwff,&
&  Bsp%lomo_spin,Bsp%homo_spin,Bsp%humo_spin,nfftot_osc,ngfft_osc,Bsp%npweps,all_mgq0)

 ! ========================
 ! ==== Resonant Block ====
 ! ========================
 if (do_resonant) then
   call timab(672,1,tsec)
   call exc_build_block(BSp,Cryst,Kmesh,Qmesh,ktabr,Gsph_x,Gsph_c,Vcp,&
&    Wfd,W,Hdr_bse,nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff,all_mgq0,.TRUE.,BS_files%out_hreso)
   call timab(672,2,tsec)
 end if

 ! ========================
 ! ==== Coupling Block ====
 ! ========================
 if (do_coupling.and.BSp%use_coupling>0) then
   call timab(673,1,tsec)
   call exc_build_block(BSp,Cryst,Kmesh,Qmesh,ktabr,Gsph_x,Gsph_c,Vcp,&
&    Wfd,W,Hdr_bse,nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff,all_mgq0,.FALSE.,BS_files%out_hcoup)
   call timab(673,2,tsec)
 end if
 !
 ! * Free memory.
 ABI_FREE(all_mgq0)

100 call timab(670,2,tsec)

end subroutine exc_build_ham
!!***

!!****f* ABINIT/wfd_all_mgq0
!! NAME
!!  wfd_all_mgq0
!!
!! FUNCTION
!!
!! INPUTS
!!  Wfd<wfd_t>=Handler for the wavefunctions.
!!  Cryst<crystal_t>=Info on the crystalline structure.
!!  Qmesh<kmesh_t>=The list of q-points for epsilon^{-1} and related symmetry tables.
!!  Gsph_x<gsphere_t>=G-sphere with the G-vectors in mgq0.
!!  Vcp<vcoul_t>=The Coulomb interaction in reciprocal space. A cutoff can be used
!!  Psps<Pseudopotential_type>=Variables related to pseudopotentials
!!  Pawtab(Psps%ntypat)<pawtab_type>=PAW tabulated starting data.
!!  Paw_pwff(Cryst%ntypat*Wfd%usepaw)<pawpwff_t>=Form factor used to calculate the onsite mat. elements of a plane wave.
!!  lomo_spin(Wfd%nsppol)=Lowest occupied band for each spin
!!  homo_spin(Wfd%nsppol)=Highest occupied band for each spin
!!  humo_spin(Wfd%nsppol)=Highest unoccupied band for each spin
!!  nfftot_osc=Total Number of FFT points used for the oscillator matrix elements.
!!  ngfft_osc(18)=Info on the FFT algorithm used to calculate the oscillator matrix elements.
!!  npweps=Number of G-vectors in mgq0.
!!
!! OUTPUT
!!   mgq0(npweps,lomo_min:humo_max,lomo_min:humo_max,Wfd%nkibz,Wfd%nsppol)
!!     Allocated here and filled with the matrix elements on each node.
!!
!! PARENTS
!!      exc_build_ham
!!
!! CHILDREN
!!      cwtime,get_bz_item,gsph_fft_tabs,paw_rho_tw_g,pawcprj_alloc
!!      pawcprj_free,pawpwij_free,pawpwij_init,rho_tw_g,timab,wfd_change_ngfft
!!      wfd_get_cprj,wfd_get_ur,wrtout,xmpi_distab,xmpi_sum
!!
!! SOURCE

subroutine wfd_all_mgq0(Wfd,Cryst,Qmesh,Gsph_x,Vcp,&
& Psps,Pawtab,Paw_pwff,lomo_spin,homo_spin,humo_spin,nfftot_osc,ngfft_osc,npweps,mgq0)

 use defs_basis
 use m_profiling_abi
 use m_xmpi
 use m_errors

 use defs_datatypes,  only : pseudopotential_type
 use m_gwdefs,        only : czero_gw
 use m_time,          only : cwtime, timab
 use m_crystal,       only : crystal_t
 use m_gsphere,       only : gsphere_t, gsph_fft_tabs
 use m_vcoul,         only : vcoul_t
 use m_bz_mesh,       only : kmesh_t, get_BZ_item
 use m_pawpwij,       only : pawpwff_t, pawpwij_t, pawpwij_init, pawpwij_free, paw_rho_tw_g
 use m_pawtab,        only : pawtab_type
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_wfd,           only : wfd_t, wfd_get_ur, wfd_get_cprj, wfd_change_ngfft, wfd_ihave_ur, wfd_distribute_bbp
 use m_oscillators,   only : rho_tw_g

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_all_mgq0'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot_osc,npweps
 type(kmesh_t),intent(in) :: Qmesh
 type(crystal_t),intent(in) :: Cryst
 type(vcoul_t),intent(in) :: Vcp
 type(gsphere_t),intent(in) :: Gsph_x
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfd_t),target,intent(inout) :: Wfd
!arrays
 integer,intent(in) :: lomo_spin(Wfd%nsppol),homo_spin(Wfd%nsppol),humo_spin(Wfd%nsppol)
 integer,intent(in) :: ngfft_osc(18)
 complex(gwpc),allocatable,intent(out) :: mgq0(:,:,:,:,:)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: map2sphere1=1,dim_rtwg1=1,ndat1=1
 integer :: use_padfft,mgfft_osc,fftalga_osc,ii
 integer :: ik_ibz,itim_k,isym_k,iq_bz,iq_ibz,isym_q,itim_q,iqbz0
 integer :: ierr,iv,ic,spin,lomo_min,humo_max !,inv_ipw,ipw
 real(dp) :: cpu,wall,gflops !q0vol,fcc_const
 complex(dpc) :: ph_mkt
 character(len=500) :: msg
!arrays
 integer,allocatable :: igfftg0(:),task_distrib(:,:,:,:)
 integer,allocatable :: gbound(:,:),id_tab(:)
 real(dp) :: qbz(3),spinrot_k(4),tsec(2)
 complex(gwpc),allocatable :: rhotwg1(:)
 complex(gwpc),target,allocatable :: ur1(:),ur2(:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: ptr_ur1(:),ptr_ur2(:)
 type(pawcprj_type),allocatable :: Cp1(:,:),Cp2(:,:)
 type(pawpwij_t),allocatable :: Pwij_q0(:)

!************************************************************************

 call timab(671,1,tsec)

 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")
 ABI_CHECK(nfftot_osc==PRODUCT(ngfft_osc(1:3)),"mismatch in FFT size")

 lomo_min = MINVAL(lomo_spin); humo_max = MAXVAL(humo_spin)

 if ( ANY(ngfft_osc(1:3) /= Wfd%ngfft(1:3)) ) then
   call wfd_change_ngfft(Wfd,Cryst,Psps,ngfft_osc)
 end if

 mgfft_osc   = MAXVAL(ngfft_osc(1:3))
 fftalga_osc = ngfft_osc(7)/100 !; fftalgc_osc=MOD(ngfft_osc(7),10)

 ! (temporary) Table used for the wavefunction in the IBZ.
 ABI_MALLOC(id_tab, (Wfd%nfft))
 id_tab = (/(ii, ii=1,Wfd%nfft)/)

 ! Analytic integration of 4pi/q^2 over the volume element:
 ! $4pi/V \int_V d^3q 1/q^2 =4pi bz_geometric_factor V^(-2/3)$
 ! i_sz=4*pi*bz_geometry_factor*q0_vol**(-two_thirds) where q0_vol= V_BZ/N_k
 ! bz_geometry_factor: sphere=7.79, fcc=7.44, sc=6.188, bcc=6.946, wz=5.255
 ! (see gwa.pdf, appendix A.4)

 ! If q=0 and C=V then set up rho-twiddle(G=0) to reflect an
 ! analytic integration of q**-2 over the volume element:
 ! <q**-2> = 7.44 V**(-2/3)   (for fcc cell)

 ! q0vol = (8.0*pi**3) / (Cryst%ucvol*Kmesh%nbz)
 ! fcc_const = SQRT(7.44*q0vol**(-2.0/3.0))
 ! rtw = (6.0*pi**2/(Cryst%ucvol*Kmesh%nkbz))**(1./3.)
 ! Average of (q+q')**-2 integration for head of Coulomb matrix
 ! INTRTW(QL) = (2*pi*rtw + pi*(rtw**2/QL-QL)*LOG((QL+rtw)/(QL-rtw)))
 ! &              * (Cryst%ucvol*Kmesh%nbz)/(2*pi)**3. * QL*QL

 if (Wfd%usepaw==1) then
   ABI_DT_MALLOC(Cp1,(Wfd%natom,Wfd%nspinor))
   call pawcprj_alloc(Cp1,0,Wfd%nlmn_atm)
   ABI_DT_MALLOC(Cp2,(Wfd%natom,Wfd%nspinor))
   call pawcprj_alloc(Cp2,0,Wfd%nlmn_atm)
 end if

 ABI_MALLOC(ur1,(nfftot_osc*Wfd%nspinor))
 ABI_MALLOC(ur2,(nfftot_osc*Wfd%nspinor))

 ! Identify q==0
 iqbz0=0
 do iq_bz=1,Qmesh%nbz
   if (ALL(ABS(Qmesh%bz(:,iq_bz))<tol3)) iqbz0 = iq_bz
 end do
 ABI_CHECK(iqbz0/=0,"q=0 not found in q-point list!")

 ! * Get iq_ibz, and symmetries from iqbz0.
 call get_BZ_item(Qmesh,iqbz0,qbz,iq_ibz,isym_q,itim_q)

 if (Wfd%usepaw==1) then ! Prepare onsite contributions at q==0
   ABI_DT_MALLOC(Pwij_q0,(Cryst%ntypat))
   call pawpwij_init(Pwij_q0,npweps,Qmesh%bz(:,iqbz0),Gsph_x%gvec,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
 end if
 !
 ! Tables for the FFT of the oscillators.
 !  a) FFT index of the G sphere (only vertical transitions, unlike cchi0, no need to shift the sphere).
 !  b) gbound table for the zero-padded FFT performed in rhotwg.
 ABI_MALLOC(igfftg0,(Gsph_x%ng))
 ABI_MALLOC(gbound,(2*mgfft_osc+8,2))
 call gsph_fft_tabs(Gsph_x,(/0,0,0/),mgfft_osc,ngfft_osc,use_padfft,gbound,igfftg0)
 if ( ANY(fftalga_osc == (/2,4/)) ) use_padfft=0 ! Pad-FFT is not coded in rho_tw_g
 if (use_padfft==0) then
   ABI_FREE(gbound)
   ABI_MALLOC(gbound,(2*mgfft_osc+8,2*use_padfft))
 end if

 ABI_MALLOC(rhotwg1,(npweps))

 ABI_STAT_MALLOC(mgq0, (npweps,lomo_min:humo_max,lomo_min:humo_max,Wfd%nkibz,Wfd%nsppol), ierr)
 ABI_CHECK(ierr==0, "out-of-memory in mgq0")
 mgq0 = czero

 call cwtime(cpu,wall,gflops,"start")

 do spin=1,Wfd%nsppol
   ! Distribute the calculation of the matrix elements.
   ! processors have the entire set of wavefunctions hence we divide the workload
   ! without checking if the pair of states is available. Last dimension is fake.
   ABI_MALLOC(task_distrib,(lomo_spin(spin):humo_spin(spin),lomo_spin(spin):humo_spin(spin),Wfd%nkibz,1))
   call xmpi_distab(Wfd%nproc,task_distrib)

   ! loop over the k-points in IBZ
   do ik_ibz=1,Wfd%nkibz
     if ( ALL(task_distrib(:,:,ik_ibz,1)/= Wfd%my_rank) ) CYCLE

     ! Don't need to symmetrize the wavefunctions.
     itim_k=1; isym_k=1; ph_mkt=cone; spinrot_k=Cryst%spinrot(:,isym_k)

     do iv=lomo_spin(spin),humo_spin(spin) ! Loop over band V
       if ( ALL(task_distrib(:,iv,ik_ibz,1)/=Wfd%my_rank) ) CYCLE

       if (wfd_ihave_ur(Wfd,iv,ik_ibz,spin,how="Stored")) then
         ptr_ur1 =>  Wfd%Wave(iv,ik_ibz,spin)%ur
       else
         call wfd_get_ur(Wfd,iv,ik_ibz,spin,ur1)
         ptr_ur1 =>  ur1
       end if

       if (Wfd%usepaw==1) then
         call wfd_get_cprj(Wfd,iv,ik_ibz,spin,Cryst,Cp1,sorted=.FALSE.)
       end if

       ! Loop over band C
       do ic=lomo_spin(spin),humo_spin(spin)
         if ( task_distrib(ic,iv,ik_ibz,1)/=Wfd%my_rank ) CYCLE

         if (wfd_ihave_ur(Wfd,ic,ik_ibz,spin,how="Stored")) then
           ptr_ur2 =>  Wfd%Wave(ic,ik_ibz,spin)%ur
         else
           call wfd_get_ur(Wfd,ic,ik_ibz,spin,ur2)
           ptr_ur2 =>  ur2
         end if

         if (Wfd%usepaw==1) then
           call wfd_get_cprj(Wfd,ic,ik_ibz,spin,Cryst,Cp2,sorted=.FALSE.)
         end if

         call rho_tw_g(Wfd%nspinor,npweps,nfftot_osc,ndat1,ngfft_osc,map2sphere1,use_padfft,igfftg0,gbound,&
&          ptr_ur1,1,id_tab,ph_mkt,spinrot_k,&
&          ptr_ur2,1,id_tab,ph_mkt,spinrot_k,&
&          dim_rtwg1,rhotwg1)

         if (Wfd%usepaw==1) then ! Add PAW onsite contribution.
           call paw_rho_tw_g(npweps,dim_rtwg1,Wfd%nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,Gsph_x%gvec,&
&            Cp1,Cp2,Pwij_q0,rhotwg1)
         end if

         ! If q=0 treat Exchange and Coulomb-term independently
         if (iv <= homo_spin(spin) .and. ic <= homo_spin(spin) .or. &
&            iv >  homo_spin(spin) .and. ic >  homo_spin(spin)) then

           if (iv/=ic) then !COULOMB term: C/=V: ignore them
             rhotwg1(1) = czero_gw
           else
             ! If q=0 and C=V then set up rho-twiddle(G=0) to reflect an
             ! analytic integration of q**-2 over the volume element:
             ! <q**-2> = 7.44 V**(-2/3)   (for fcc cell)
             !rhotwg1(1) = fcc_const * qpg(1,iqbz0)
             rhotwg1(1) = SQRT(GWPC_CMPLX(Vcp%i_sz,zero)) / Vcp%vcqlwl_sqrt(1,1)
             !if (vcut) rhotwg1(1) = 1.0
           end if

         else
           ! At present this term is set to zero
           ! EXCHANGE term: limit value.
           ! Set up rho-twiddle(G=0) using small vector q instead of zero and k.p perturbation theory (see notes)
           rhotwg1(1) = czero_gw
         end if

         mgq0(:,iv,ic,ik_ibz,spin) = rhotwg1(:)
       end do !ic
     end do !iv
   end do !ik_ibz

   ABI_FREE(task_distrib)
 end do !spin

 ! TODO: One can speedup the calculation by computing the upper triangle of the
 ! matrix in (b,b') space and then take advantage of the symmetry property:
 !
 !   M_{k,0}{{bb'}(G)^* = M{k,0}{b'b'}(-G)

#if 0
 !!!! $OMP PARALLEL DO COLLAPSE(3) PRIVATE(inv_ipw)
 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz
     do iv=lomo_spin(spin),humo_spin(spin)
       do ic=1,iv-1
         do ipw=1,npweps
           inv_ipw = gsph_x%g2mg(ipw)
           mgq0(inv_ipw,ic,iv,ik_ibz,spin) = mgq0(ipw,iv,ic,ik_ibz,spin)
         end do
       end do
     end do
   end do
 end do
#endif
 !
 ! Gather matrix elements on each node.
 call xmpi_sum(mgq0,Wfd%comm,ierr)

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f9.6))')"cpu_time = ",cpu,", wall_time = ",wall
 call wrtout(std_out,msg,"PERS")

 ABI_FREE(rhotwg1)
 ABI_FREE(igfftg0)
 ABI_FREE(gbound)
 ABI_FREE(ur1)
 ABI_FREE(ur2)
 ABI_FREE(id_tab)

 if (Wfd%usepaw==1) then
   ! Deallocation for PAW.
   call pawpwij_free(Pwij_q0)
   ABI_DT_FREE(Pwij_q0)
   call pawcprj_free(Cp1)
   ABI_DT_FREE(Cp1)
   call pawcprj_free(Cp2)
   ABI_DT_FREE(Cp2)
 end if

 call timab(671,2,tsec)

end subroutine wfd_all_mgq0
!!***
