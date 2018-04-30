!{\src2tex{textfont=tt}}
!!****f* ABINIT/paw_qpscgw
!! NAME
!! paw_qpscgw
!!
!! FUNCTION
!!  This routine is called during QP self-consistent GW calculations. It calculates the new QP on-site quantities
!!  using the QP amplitudes read from the QPS file.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Wfd<wfd_t>=Datatype gathering data on QP amplitudes.
!!  nscf=Number of QPSCGW iterations done so far (read from the QPS file).
!!  nfftf=Number of points in the fine FFT grid.
!!  ngfft(18)=information on the fine FFT grid used for densities and potentials.
!!  Dtset<dataset_type>=All input variables for this dataset.
!!  Cryst<crystal_t>=Info on unit cell and symmetries.
!!  Kmesh<kmesh_t>=Structure describing the k-point sampling.
!!  Psps<Pseudopotential_type)>=Info on pseudopotential, only for consistency check of the KSS file
!!  Pawang<pawang_type>=PAW angular mesh and related data.
!!  Pawrad(ntypat*usepaw)<type(pawrad_type)>=paw radial mesh and related data
!!  Pawtab(ntypat*usepaw)<type(pawtab_type)>=paw tabulated starting data
!!  Pawfgrtab(natom)<Pawfgrtab_type>= For PAW, various arrays giving data related to fine grid for a given atom.
!!  prev_Pawrhoij(Cryst%natom))<Pawrhoij_type>=Previous QP rhoij used for mixing if nscf>0 and rhoqpmix/=one.
!!  MPI_enreg=Information about MPI parallelization.
!!  QP_BSt<ebands_t>=QP band structure.
!!  QP_energies<Energies_type>=Simple datastructure to gather all part of total energy.
!!  nhatgrdim= 0 if pawgrnhat array is not used ; 1 otherwise
!!
!! OUTPUT
!!  QP_pawrhoij(Cryst%natom))<Pawrhoij_type>=on-site densities calculated from the QP amplitudes.
!!  qp_nhat(nfftf,Dtset%nspden)=Compensation charge density calculated from the QP amplitudes.
!!  qp_nhatgr(nfftf,Dtset%nspden,3*nhatgrdim)=Derivatives of the QP nhat on fine rectangular grid (and derivatives).
!!  qp_compch_sph=QP compensation charge integral inside spheres computed over spherical meshes.
!!  qp_compch_fft=QP compensation charge inside spheres computed over fine fft grid.
!!  QP_paw_ij(Cryst%natom)<Paw_ij_type>=Non-local D_ij strengths of the QP Hamiltonian.
!!  QP_paw_an(Cryst%natom)<Paw_an_type>=Various arrays related to the Hamiltonian given 
!!   on ANgular mesh or ANgular moments.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      paw_an_init,paw_an_nullify,paw_ij_init,paw_ij_nullify,pawdenpot
!!      pawmknhat,pawrhoij_alloc,pawrhoij_unpack,symrhoij,wfd_pawrhoij,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine paw_qpscgw(Wfd,nscf,nfftf,ngfftf,Dtset,Cryst,Kmesh,Psps,QP_BSt,&
&  Pawang,Pawrad,Pawtab,Pawfgrtab,prev_Pawrhoij,&
&  QP_pawrhoij,QP_paw_ij,QP_paw_an,QP_energies,qp_nhat,nhatgrdim,qp_nhatgr,qp_compch_sph,qp_compch_fft)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_gwdefs,        only : sigparams_t
 use m_bz_mesh,       only : kmesh_t
 use m_crystal,       only : crystal_t
 use m_ebands,        only : get_eneocc_vect
 use m_energies,      only : energies_type
 use m_pawang,        only : pawang_type
 use m_pawrad,        only : pawrad_type
 use m_pawtab,        only : pawtab_type
 use m_paw_an,        only : paw_an_type, paw_an_init, paw_an_nullify
 use m_paw_ij,        only : paw_ij_type, paw_ij_init, paw_ij_nullify
 use m_pawfgrtab,     only : pawfgrtab_type
 use m_pawrhoij,      only : pawrhoij_type, pawrhoij_alloc, pawrhoij_unpack, pawrhoij_get_nspden, symrhoij
 use m_wfd,           only : wfd_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_qpscgw'
 use interfaces_14_hidewrite
 use interfaces_65_paw
 use interfaces_69_wfdesc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,nscf,nhatgrdim
 real(dp),intent(out) :: qp_compch_fft,qp_compch_sph
 type(kmesh_t),intent(in) :: Kmesh
 type(crystal_t),intent(in) :: Cryst
 type(Dataset_type),intent(in) :: Dtset
 type(Pseudopotential_type),intent(in) :: Psps
 type(Pawang_type),intent(in) :: Pawang
 type(ebands_t),intent(in) :: QP_BSt
 type(Energies_type),intent(inout) :: QP_energies
 type(wfd_t),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(out) :: qp_nhat(nfftf,Dtset%nspden)
 real(dp),intent(out) :: qp_nhatgr(nfftf,Dtset%nspden,3*nhatgrdim)
 type(Pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom)
 type(Pawrad_type),intent(in) :: Pawrad(Psps%ntypat*Psps%usepaw)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 type(Pawrhoij_type),intent(inout) :: prev_Pawrhoij(Cryst%natom)
 type(Pawrhoij_type),intent(out) :: QP_pawrhoij(Cryst%natom)
 type(Paw_ij_type),intent(out) :: QP_paw_ij(Cryst%natom)
 type(Paw_an_type),intent(inout) :: QP_paw_an(Cryst%natom)

!Local variables ------------------------------
!scalars
 integer :: choice,cplex,has_dijU,has_dijso,iat,ider,idir,ipert
 integer :: izero,nkxc1,nspden_rhoij,nzlmopt
 integer :: option,optrhoij,usexcnhat
 character(len=500) :: msg
!arrays
 real(dp) :: k0(3)

!************************************************************************

 ABI_UNUSED(Kmesh%nibz)
 !
 !  * 0 if Vloc in atomic data is Vbare    (Blochl s formulation)
 !  * 1 if Vloc in atomic data is VH(tnzc) (Kresse s formulation)
 usexcnhat=MAXVAL(Pawtab(:)%usexcnhat)
 !
 ! Calculate new rhoij_qp from updated Cprj_ibz, note use_rhoij_=1.
 nspden_rhoij=pawrhoij_get_nspden(Dtset%nspden,Dtset%nspinor,Dtset%pawspnorb)

 call pawrhoij_alloc(QP_pawrhoij,Dtset%pawcpxocc,nspden_rhoij,Dtset%nspinor,Dtset%nsppol,Cryst%typat,&
&                 pawtab=Pawtab,use_rhoij_=1,use_rhoijres=1)

 ! FIXME kptop should be passed via Kmesh, in GW time reversal is always assumed.
 call wfd_pawrhoij(Wfd,Cryst,QP_Bst,Dtset%kptopt,QP_pawrhoij,Dtset%pawprtvol)
 !
 ! * Symmetrize QP $\rho_{ij}$.
 choice=1; optrhoij=1; ipert=0
 call symrhoij(QP_pawrhoij,QP_pawrhoij,choice,Cryst%gprimd,Cryst%indsym,ipert,&
&              Cryst%natom,Cryst%nsym,Cryst%ntypat,optrhoij,Pawang,Dtset%pawprtvol,&
&              Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)

 ! ======================
 ! ==== Make QP nhat ====
 ! ======================
 cplex=1; ider=2*nhatgrdim; idir=0; ipert=0; izero=0; k0(:)=zero

 call pawmknhat(qp_compch_fft,cplex,ider,idir,ipert,izero,Cryst%gprimd,&
&  Cryst%natom,Cryst%natom,nfftf,ngfftf,nhatgrdim,Dtset%nspden,Cryst%ntypat,Pawang,&
&  Pawfgrtab,qp_nhatgr,qp_nhat,QP_Pawrhoij,QP_Pawrhoij,Pawtab,k0,Cryst%rprimd,&
&  Cryst%ucvol,Dtset%usewvl,Cryst%xred)

 ! Allocate quantities related to the PAW spheres for the QP Hamiltonian.
 ! TODO call paw_ij_init in scfcv and respfn, fix small issues
 has_dijso=Dtset%pawspnorb; has_dijU=Dtset%usepawu

 call paw_ij_nullify(QP_paw_ij)
 call paw_ij_init(QP_paw_ij,cplex,Dtset%nspinor,Dtset%nsppol,&
&  Dtset%nspden,Dtset%pawspnorb,Cryst%natom,Cryst%ntypat,Cryst%typat,Pawtab,&
&  has_dij=1,has_dijhartree=1,has_dijhat=1,has_dijxc=1,has_dijxc_hat=1,has_dijxc_val=1,&
&  has_dijso=has_dijso,has_dijU=has_dijU,has_exexch_pot=1,has_pawu_occ=1)

 call paw_an_nullify(QP_paw_an); nkxc1=0 ! No kernel
 call paw_an_init(QP_paw_an,Cryst%natom,Cryst%ntypat,nkxc1,Dtset%nspden,&
&     cplex,Dtset%pawxcdev,Cryst%typat,Pawang,Pawtab,has_vxc=1,has_vxcval=1)

 ! =====================================================
 ! ==== Optional mixing of the PAW onsite densities ====
 ! =====================================================
 if (nscf>0 .and. (ABS(Dtset%rhoqpmix-one)>tol12) ) then
   write(msg,'(a,f6.3)')' sigma: mixing on-site QP rho_ij densities using rhoqpmix= ',Dtset%rhoqpmix
   call wrtout(std_out,msg,'COLL')
   ! qp_rhor = prev_rhor + Dtset%rhoqpmix*(qp_rhor-prev_rhor)

   call pawrhoij_unpack(QP_Pawrhoij)   ! Unpack new QP %rhoijp
   call pawrhoij_unpack(prev_Pawrhoij) ! Unpack previous QP %rhoijp

   do iat=1,Cryst%natom
     QP_pawrhoij(iat)%rhoij_ = prev_Pawrhoij(iat)%rhoij_ &
&      + Dtset%rhoqpmix * (QP_pawrhoij(iat)%rhoij_ - prev_pawrhoij(iat)%rhoij_)

     prev_pawrhoij(iat)%use_rhoij_=0
     ABI_DEALLOCATE(prev_pawrhoij(iat)%rhoij_)
   end do
   !
   ! * Re-Symmetrize mixed QP $\rho_{ij}$.
   choice=1; optrhoij=1; ipert=0
   call symrhoij(QP_pawrhoij,QP_pawrhoij,choice,Cryst%gprimd,Cryst%indsym,ipert,&
&                Cryst%natom,Cryst%nsym,Cryst%ntypat,optrhoij,Pawang,Dtset%pawprtvol,&
&                Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)
 end if

 do iat=1,Cryst%natom
   QP_pawrhoij(iat)%use_rhoij_=0
   ABI_DEALLOCATE(QP_pawrhoij(iat)%rhoij_)
 end do
 !
 ! =================================================================================
 ! ==== Evaluate on-site energies, potentials, densities using (mixed) QP rhoij ====
 ! =================================================================================
 ! * Initialize also "lmselect" (index of non-zero LM-moments of densities).

 nzlmopt=-1; option=0; qp_compch_sph=greatest_real

 call pawdenpot(qp_compch_sph,QP_energies%e_paw,QP_energies%e_pawdc,&
&  ipert,Dtset%ixc,Cryst%natom,Cryst%natom,Dtset%nspden,&
&  Cryst%ntypat,Dtset%nucdipmom,nzlmopt,option,QP_paw_an,QP_paw_an,&
&  QP_paw_ij,Pawang,Dtset%pawprtvol,Pawrad,QP_pawrhoij,Dtset%pawspnorb,&
&  Pawtab,Dtset%pawxcdev,Dtset%spnorbscl,Dtset%xclevel,Dtset%xc_denpos,Cryst%ucvol,Psps%znuclpsp)

end subroutine paw_qpscgw
!!***
