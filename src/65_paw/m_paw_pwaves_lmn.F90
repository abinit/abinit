!!****m* ABINIT/m_paw_pwaves_lmn
!! NAME
!!  m_paw_pwaves_lmn
!!
!! FUNCTION
!!  This module contains the definition of the paw_pwaves_lmn structured datatype,
!!  as well as related functions and methods.
!!  paw_pwaves_lmn is used to store the 3D values of the all-electron and of
!!  the pseudized part of the PAW partial waves on the set of FFT points falling
!!  inside the spheres around each atom.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! * Routines tagged with "@type_name" are strongly connected to the definition of the data type.
!!   Strongly connected means that the proper functioning of the implementation relies on the
!!   assumption that the tagged procedure is consistent with the type declaration.
!!   Every time a developer changes the structure "type_name" adding new entries, he/she has to make sure
!!   that all the strongly connected routines are changed accordingly to accommodate the modification of the data type.
!!   Typical examples of strongly connected routines are creation, destruction or reset methods.
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

MODULE m_paw_pwaves_lmn

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_sort

 use m_numeric_tools, only : wrap2_zero_one
 use m_geometry,      only : xcart2xred
 use m_pawrad,        only : pawrad_type, nderiv_gen, pawrad_deducer0
 use m_pawtab,        only : pawtab_type
 use m_pawfgrtab,     only : pawfgrtab_type
 use m_paw_numeric,   only : paw_spline, paw_splint
 use m_paw_sphharm,   only : initylmr
 use m_paral_atom,    only : get_my_atmtab, free_my_atmtab

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_paw_pwaves_lmn/paw_pwaves_lmn_t
!! NAME
!! paw_pwaves_lmn_t
!!
!! FUNCTION
!!  Datatype used to store the 3D values of the all-electron and of the pseudized part of the
!!  PAW partial waves on the set of FFT points falling inside the spheres around each atom.
!!  The data is mainly used for plotting the true PAW wavefunctions in real space.
!!
!! SOURCE

 type,public ::  paw_pwaves_lmn_t

  integer :: nfgd

  integer :: lmn_size

  !$integer :: ngfft(18)

  integer,allocatable :: r0shift(:,:)
  ! r0shift(3,nfgd)

  !real(dp),allocatable:: phk_atm(:,:)
  ! phk_atmt(2,nfgd)

  real(dp),allocatable :: phi(:,:)
  ! phi (nfgd,lmn_size)
  ! \phi_{nlm}(ifgd) for each point of the FFT mesh located inside the PAW sphere (see pawfgrtab_type).

  real(dp),allocatable :: tphi(:,:)
  ! tphi (nfgd,lmn_size)
  ! \tphi_{nlm}(ifgd) for each point of the FFT mesh located inside the PAW sphere (see pawfgrtab_type).

  real(dp),allocatable :: phi_gr(:,:,:)
  ! phi_gr (3,nfgd,lmn_size)
  ! gradient, in cartesian coordinates, of \phi_{nlm}(ifgd) for each point of the FFT mesh
  ! located inside the PAW sphere (see pawfgrtab_type).

  real(dp),allocatable :: tphi_gr(:,:,:)
  ! tphi_gr (3,nfgd,lmn_size)
  ! gradient, in cartesian coordinates, of \tphi_{nlm}(ifgd) for each point of the FFT mesh
  ! located inside the PAW sphere (see pawfgrtab_type).

 end type paw_pwaves_lmn_t

 public :: paw_pwaves_lmn_init
 public :: paw_pwaves_lmn_free
!!***

CONTAINS !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_pwaves_lmn/paw_pwaves_lmn_init
!! NAME
!! paw_pwaves_lmn_init
!!
!! FUNCTION
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!
!! OUTPUT
!!
!! PARENTS
!!      classify_bands,exc_plot,m_wfd,pawmkaewf,screening,sigma,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_pwaves_lmn_init(Paw_onsite,my_natom,natom,ntypat,rprimd,xcart,Pawtab, &
&                              Pawrad,local_pawfgrtab,optgrad,&
&                              mpi_atmtab,comm_atom) ! optional arguments (parallelism)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,natom,ntypat
 integer,optional,intent(in) :: optgrad,comm_atom
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: xcart(3,natom),rprimd(3,3)
 type(pawtab_type),target,intent(in) :: Pawtab(ntypat)
 type(pawrad_type),intent(in) :: Pawrad(ntypat)
 type(pawfgrtab_type),intent(in) :: local_pawfgrtab(my_natom)
 type(paw_pwaves_lmn_t),intent(out) :: Paw_onsite(my_natom)

!Local variables-------------------------------
!scalars
 integer :: itypat,ln_size,lmn_size,mesh_size,inl,iatom,iatom1,my_comm_atom,my_optgrad
 integer :: nfgd,ifgd,ipsang,option_ylmr,normchoice,ii,jlmn,jl,jm,jlm,jln
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: phj,rR,tphj,ybcbeg,ybcend
!arrays
 integer, allocatable :: iperm(:)
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: yvals(4),red(3),shift(3)
 real(dp),allocatable :: ff(:),nrm(:),nrm_sort(:),phigrd(:,:),tphigrd(:,:),ylm_tmp(:,:),ylm(:,:),ylm_gr(:,:,:)
 real(dp),allocatable :: rsph_red(:,:),rsph_cart(:,:),phigrd_gr(:,:),tphigrd_gr(:,:),gg(:)
 type(paw_pwaves_lmn_t),allocatable :: Paw_lmn_spline(:)

! *************************************************************************

!@paw_pwaves_lmn_t

 if (my_natom==0) return
 my_optgrad = -1; if (present(optgrad)) my_optgrad = optgrad

 ABI_CHECK(all(local_pawfgrtab(:)%rfgd_allocated==1),"R vectors not allocated in pawfgrtab")

 ! Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 ! Prepare the spline. Calculate 2nd derivatives of partial waves for each atom type.
 ABI_MALLOC(Paw_lmn_spline,(ntypat))

 do itypat=1,ntypat
   ln_size  =Pawtab(itypat)%basis_size
   mesh_size=Pawtab(itypat)%mesh_size

   ABI_MALLOC(Paw_lmn_spline(itypat)%phi ,(mesh_size,ln_size))
   ABI_MALLOC(Paw_lmn_spline(itypat)%tphi,(mesh_size,ln_size))

   do inl=1,ln_size ! Calculate 2nd derivatives of %phi and %tphi for each ln component.
     ybcbeg=zero; ybcend=zero
     call paw_spline(Pawrad(itypat)%rad,Pawtab(itypat)%phi(:,inl), mesh_size,&
&                    ybcbeg,ybcend,Paw_lmn_spline(itypat)%phi(:,inl))

     ybcbeg=zero; ybcend=zero
     call paw_spline(Pawrad(itypat)%rad,Pawtab(itypat)%tphi(:,inl),mesh_size,&
&                    ybcbeg,ybcend,Paw_lmn_spline(itypat)%tphi(:,inl))
   end do
 end do
 !
 ! === spline for each atom ===
 ! * FFT points within PAW sphere depend on the atom site.
 do iatom1=1,my_natom
   iatom=iatom1;if (paral_atom) iatom=my_atmtab(iatom1)

   itypat   = local_pawfgrtab(iatom1)%itypat
   ln_size  = Pawtab(itypat)%basis_size
   lmn_size = Pawtab(itypat)%lmn_size
   mesh_size= Pawrad(itypat)%mesh_size
   nfgd     = local_pawfgrtab(iatom1)%nfgd ! no. of points in the fine grid for this PAW sphere
   indlmn   => Pawtab(itypat)%indlmn

   ! The points in the PAW sphere might belong to a different unit cell, in this case one has to
   ! reconstruct the contribution to the AE psi_k given by the lattice-symmetric atom in another sphere.
   ! Here I wrap rr back into the first unit cell keeping trace of the lattice vector L
   ! needed so that rr = rr_first_ucell + L
   ! The contribution to the AE u(r) in the first unit cell has to be multiplied by e^{-ikL}.

   Paw_onsite(iatom1)%nfgd     = nfgd
   Paw_onsite(iatom1)%lmn_size = lmn_size
   !
   ABI_MALLOC(Paw_onsite(iatom1)%r0shift,(3,nfgd))

   ABI_MALLOC(rsph_red,(3,nfgd))
   ABI_MALLOC(rsph_cart,(3,nfgd))
   do ifgd=1,nfgd
     rsph_cart(:,ifgd) = local_pawfgrtab(iatom1)%rfgd(:,ifgd) + xcart(:,iatom)
   end do
   call xcart2xred(nfgd,rprimd,rsph_cart,rsph_red) ! go to reduced coordinates.
   ABI_FREE(rsph_cart)

   do ifgd=1,nfgd
     call wrap2_zero_one(rsph_red(1,ifgd),red(1),shift(1)) ! rr = r_cell + shift
     call wrap2_zero_one(rsph_red(2,ifgd),red(2),shift(2))
     call wrap2_zero_one(rsph_red(3,ifgd),red(3),shift(3))
     Paw_onsite(iatom1)%r0shift(:,ifgd) = NINT(shift)
     !if (ANY( ABS(shift) > tol12)) then
       !MSG_WARNING("rmR_red is outside the first unit cell.")
       !write(ab_out,*)rsph_red(:,ifgd),shift
     !end if
   end do
   ABI_FREE(rsph_red)
   !
   ! * Obtain |r-R| on fine grid, note that rfgd is given in Cartesian coordinates.
   ABI_MALLOC(nrm,(nfgd))
   do ifgd=1,nfgd
     nrm(ifgd) = sqrt(dot_product(local_pawfgrtab(iatom1)%rfgd(:,ifgd),local_pawfgrtab(iatom1)%rfgd(:,ifgd)))
   end do
   !
   ! * Compute Ylm for each r-R vector.
   ipsang = 1 + (Pawtab(itypat)%l_size-1)/2 ! recall l_size=2*l_max-1 where l_max is shifted by 1.
   ABI_MALLOC(ylm_tmp,(ipsang**2,nfgd))
   normchoice = 1 ! Use computed norms of input vectors.
   if (my_optgrad==1) then
     option_ylmr=2 ! Compute Ylm(r-R) and its gradient
     ABI_MALLOC(ylm_gr,(3,ipsang**2,nfgd))
   else
     option_ylmr= 1 ! To compute Ylm(r-R).
     ABI_MALLOC(ylm_gr,(3,3,0))
   end if
   call initylmr(ipsang,normchoice,nfgd,nrm,option_ylmr,local_pawfgrtab(iatom1)%rfgd,ylm_tmp,ylm_gr)
   !
   !  Exchange dimensions for better memory access.
   ABI_MALLOC(ylm,(nfgd,ipsang**2))
   do ii=1,ipsang**2
     ylm(:,ii) = ylm_tmp(ii,:)
   end do
   ABI_FREE(ylm_tmp)
   !
   ! In order to do spline fits, the |r-R| data must be sorted
   ! Here we sort the nrm points, keeping track of which goes where
   ABI_MALLOC(nrm_sort,(nfgd))
   nrm_sort = nrm

   ABI_MALLOC(iperm,(nfgd))
   do ifgd=1,nfgd
     iperm(ifgd)=ifgd
   end do

   call sort_dp(nfgd,nrm_sort,iperm,tol8)

   ! Now make spline fits of phi and tphi  onto the fine grid around the atom
   ABI_MALLOC(phigrd,(nfgd,ln_size))
   ABI_MALLOC(tphigrd,(nfgd,ln_size))
   ABI_MALLOC(ff,(nfgd))
   if (my_optgrad==1) then
     ABI_MALLOC(phigrd_gr,(nfgd,ln_size))
     ABI_MALLOC(tphigrd_gr,(nfgd,ln_size))
     ABI_MALLOC(gg,(mesh_size))
   end if

   do inl=1,ln_size
     !
     ! * splint phi onto points and reorder indices.
     call paw_splint(mesh_size,Pawrad(itypat)%rad,Pawtab(itypat)%phi(:,inl),&
&                    Paw_lmn_spline(itypat)%phi(:,inl),nfgd,nrm_sort,ff)
     do ifgd=1,nfgd
       ii=iperm(ifgd)
       phigrd(ii,inl) = ff(ifgd)
     end do
     !
     ! * compute d phi/dr, interpolate onto points and reorder indices.
     if (my_optgrad==1) then
       ybcbeg=zero; ybcend=zero
       call nderiv_gen(gg,Pawtab(itypat)%phi(:,inl),Pawrad(itypat))
       call paw_spline(Pawrad(itypat)%rad,gg,mesh_size,ybcbeg,ybcend,&
&                      Paw_lmn_spline(itypat)%phi(:,inl))
       call paw_splint(mesh_size,Pawrad(itypat)%rad,Paw_lmn_spline(itypat)%phi(:,inl),&
&                      Paw_lmn_spline(itypat)%phi(:,inl),nfgd,nrm_sort,ff)
       do ifgd=1,nfgd
         ii=iperm(ifgd)
         phigrd_gr(ii,inl)  = ff(ifgd)
       end do
     end if
     !
     ! * compute d tphi/dr, interpolate onto points and reorder indices.
     call paw_splint(mesh_size,Pawrad(itypat)%rad,Pawtab(itypat)%tphi(:,inl),&
&                    Paw_lmn_spline(itypat)%tphi(:,inl),nfgd,nrm_sort,ff)
     do ifgd=1,nfgd
       ii=iperm(ifgd)
       tphigrd(ii,inl) = ff(ifgd)
     end do
     if (my_optgrad==1) then
       ybcbeg=zero; ybcend=zero
       call nderiv_gen(gg,Pawtab(itypat)%tphi(:,inl),Pawrad(itypat))
       call paw_spline(Pawrad(itypat)%rad,gg,mesh_size,ybcbeg,ybcend,&
&                      Paw_lmn_spline(itypat)%tphi(:,inl))
       call paw_splint(mesh_size,Pawrad(itypat)%rad,Paw_lmn_spline(itypat)%tphi(:,inl),&
&                      Paw_lmn_spline(itypat)%phi(:,inl),nfgd,nrm_sort,ff)
       do ifgd=1,nfgd
         ii=iperm(ifgd)
         tphigrd_gr(ii,inl)  = ff(ifgd)
       end do
     end if
   end do !inl

   ABI_FREE(ff)
   if (my_optgrad==1) then
     ABI_FREE(gg)
   end if
   !
   ! === Calculate AE and PS partial waves inside the sphere ===
   ! * recall that <r|phi>=u(r)*Slm(r^)/r, hence avoid division by zero except for s-waves.
   ABI_MALLOC(Paw_onsite(iatom1)%phi ,(nfgd,lmn_size))
   ABI_MALLOC(Paw_onsite(iatom1)%tphi,(nfgd,lmn_size))

   if (my_optgrad==1) then
     ABI_MALLOC(Paw_onsite(iatom1)%phi_gr ,(3,nfgd,lmn_size))
     ABI_MALLOC(Paw_onsite(iatom1)%tphi_gr,(3,nfgd,lmn_size))
   end if

   do jlmn=1,lmn_size
     jl  = indlmn(1,jlmn)
     jm  = indlmn(2,jlmn)
     jlm = indlmn(4,jlmn)
     jln = indlmn(5,jlmn)

     do ifgd=1,nfgd ! loop over fine grid points in current PAW sphere
       !if (nrm(ifgd)>tol16) then
       if (nrm(ifgd)>tol10) then ! tol10 to be consistent with initylmr.
         rR  = nrm(ifgd) ! value of |r-R|
         !write(ab_out,*) 'rR:',rR,' phigrd:',phigrd(ifgd,jln),' tphigrd:',tphigrd(ifgd,jln),' ylm:',ylm(ifgd,jlm)
         phj = phigrd (ifgd,jln)*ylm(ifgd,jlm)/rR
         tphj= tphigrd(ifgd,jln)*ylm(ifgd,jlm)/rR
         Paw_onsite(iatom1)%phi (ifgd,jlmn) =  phj
         Paw_onsite(iatom1)%tphi(ifgd,jlmn) = tphj

         if (my_optgrad==1) then
           Paw_onsite(iatom1)%phi_gr (1:3,ifgd,jlmn) = phigrd (ifgd,jln)*ylm_gr(1:3,jlm,ifgd) &
&            + phigrd_gr(ifgd,jln)*local_pawfgrtab(iatom1)%rfgd(1:3,ifgd)*ylm(ifgd,jlm)
           Paw_onsite(iatom1)%tphi_gr (1:3,ifgd,jlmn) = tphigrd (ifgd,jln)*ylm_gr(1:3,jlm,ifgd) &
&            + tphigrd_gr(ifgd,jln)*local_pawfgrtab(iatom1)%rfgd(1:3,ifgd)*ylm(ifgd,jlm)
         end if

       else ! Extrapolate if the point is at the origin
         yvals(1) = zero
         if (jl==0) then
           yvals(2:4) = Pawtab(itypat)%phi(2:4,jln)/Pawrad(itypat)%rad(2:4)
           call pawrad_deducer0(yvals,4,Pawrad(itypat))
         end if
         Paw_onsite(iatom1)%phi(ifgd,jlmn) = yvals(1) * ylm(ifgd,jlm)
         yvals(1) = zero
         if (jl==0) then
           yvals(2:4) = Pawtab(itypat)%tphi(2:4,jln)/Pawrad(itypat)%rad(2:4)
           call pawrad_deducer0(yvals,4,pawrad(itypat))
         end if
         Paw_onsite(iatom1)%tphi(ifgd,jlmn) = yvals(1) * ylm(ifgd,jlm)
         ! The gradient is expected to go to zero at the origin
         if (my_optgrad==1) then
           Paw_onsite(iatom1)%phi_gr (1:3,ifgd,jlmn) = zero
           Paw_onsite(iatom1)%tphi_gr (1:3,ifgd,jlmn) = zero
         end if
       end if

     end do !nfgd
   end do !jlmn

   ABI_FREE(nrm)
   ABI_FREE(nrm_sort)
   ABI_FREE(iperm)
   ABI_FREE(phigrd)
   ABI_FREE(tphigrd)
   ABI_FREE(ylm)
   ABI_FREE(ylm_gr)
   if (my_optgrad==1) then
     ABI_FREE(phigrd_gr)
     ABI_FREE(tphigrd_gr)
   end if
 end do !iatom
 !
 !* Free 2nd derivates used for spline.
 call paw_pwaves_lmn_free(Paw_lmn_spline)
 ABI_DATATYPE_DEALLOCATE(Paw_lmn_spline)

 ! Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine paw_pwaves_lmn_init
!!***

!----------------------------------------------------------------------

!!****f* m_paw_pwaves_lmn/paw_pwaves_lmn_free
!! NAME
!! paw_pwaves_lmn_free
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      classify_bands,exc_plot,m_paw_pwaves_lmn,m_wfd,pawmkaewf,screening
!!      sigma,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_pwaves_lmn_free(Paw_onsite)

 implicit none

!Arguments ------------------------------------
!arrays
 type(paw_pwaves_lmn_t),intent(inout) :: Paw_onsite(:)

!Local variables-------------------------------
!scalars
 integer :: iatom

! *************************************************************************

!@paw_pwaves_lmn_t

 do iatom=LBOUND(Paw_onsite,DIM=1),UBOUND(Paw_onsite,DIM=1)
   if (allocated(Paw_onsite(iatom)%phi)) then
     ABI_FREE(Paw_onsite(iatom)%phi)
   end if
   if (allocated(Paw_onsite(iatom)%tphi)) then
     ABI_FREE(Paw_onsite(iatom)%tphi)
   end if
   if (allocated(Paw_onsite(iatom)%r0shift)) then
     ABI_FREE(Paw_onsite(iatom)%r0shift)
   end if
   if (allocated(Paw_onsite(iatom)%phi_gr )) then
     ABI_FREE(Paw_onsite(iatom)%phi_gr)
   end if
   if (allocated(Paw_onsite(iatom)%tphi_gr)) then
     ABI_FREE(Paw_onsite(iatom)%tphi_gr)
   end if
 end do

end subroutine paw_pwaves_lmn_free
!!***

!----------------------------------------------------------------------

END MODULE m_paw_pwaves_lmn
!!***
