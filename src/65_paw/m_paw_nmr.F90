!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_nmr
!! NAME
!!  m_paw_nmr
!!
!! FUNCTION
!!  This module contains routines related to Nuclear Magnetic Resonance (NMR)
!!   observables (PAW approach).
!!
!! COPYRIGHT
!! Copyright (C) 2018-2019 ABINIT group (JWZ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_nmr


 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi

 use m_symtk,      only : matpointsym
 use m_pawang,     only : pawang_type
 use m_paw_sphharm, only : slxyzs
 use m_pawtab,     only : pawtab_type
 use m_pawrad,     only : pawrad_type,pawrad_deducer0,simp_gen
 use m_pawtab,     only : pawtab_type
 use m_paw_an,     only : paw_an_type
 use m_pawrhoij,   only : pawrhoij_type
 use m_paw_denpot, only : pawdensities
 use m_paral_atom, only : get_my_atmtab,free_my_atmtab

 implicit none

 private

!public procedures.
 public :: make_efg_onsite ! Compute the electric field gradient due to PAW on-site densities
 public :: make_fc_paw     ! Compute the PAW on-site contribution to the Fermi-contact
 public :: make_orbl_paw     ! Compute the PAW on-site contribution for orbital magnetism, 1/2<L>

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_nmr/make_orbl_paw
!! NAME
!! make_orbl_paw
!!
!! FUNCTION
!! Compute the onsite contribution to orbital magnetism, 1/2 <L>
!!
!! INPUTS
!!  idir=cartesian direction of interest
!!  natom=number of atoms in cell.
!!  ntypat=number of atom types
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat(ntypat)
!!
!! OUTPUT
!!  orbl=complex(dpc) 1/2<L_dir>

!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_orbl_paw(idir,natom,ntypat,orbl,pawrad,pawtab,typat)

 implicit none

!Arguments ------------------------------------
 !scalars
 integer,intent(in) :: idir,natom,ntypat
 complex(dpc),intent(out) :: orbl
 !arrays
 integer,intent(in) :: typat(ntypat)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),target,intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
 !scalars
 integer :: iatom,il,im,ilmn,itypat,jl,jm,jlmn,klmn,kln,mesh_size
 real(dp) :: intg
 complex(dpc) :: orbl_me
 !arrays
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp),allocatable :: ff(:)

! ************************************************************************

 DBG_ENTER("COLL")

 !loop over atoms in cell
 orbl = czero
 do iatom = 1, natom
   itypat=typat(iatom)
   indlmn => pawtab(itypat)%indlmn
   
   mesh_size=pawtab(itypat)%mesh_size
   ABI_ALLOCATE(ff,(mesh_size))

!    loop over basis elements for this atom
!    ----
     do jlmn=1,pawtab(itypat)%lmn_size
       jl=indlmn(1,jlmn)
       jm=indlmn(2,jlmn)
       do ilmn=1,pawtab(itypat)%lmn_size
         il=indlmn(1,ilmn)
         im=indlmn(2,ilmn)

         klmn=max(jlmn,ilmn)*(max(jlmn,ilmn)-1)/2 + min(jlmn,ilmn)
         kln = pawtab(itypat)%indklmn(2,klmn) ! need this for mesh selection below

         ! compute <L_dir>
         call slxyzs(jl,jm,idir,il,im,orbl_me)

         ! compute integral of phi_i*phi_j - tphi_i*tphi_j
         if (abs(orbl_me) > tol8) then
            ff(1:mesh_size)=pawtab(itypat)%phiphj(1:mesh_size,kln) - pawtab(itypat)%tphitphj(1:mesh_size,kln)
            call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
            call simp_gen(intg,ff,pawrad(itypat))
            
            orbl = orbl + half*orbl_me*intg

         end if ! end check that |L_dir| > 0, otherwise ignore term

      end do ! end loop over ilmn
   end do ! end loop over jlmn
      
   ABI_DEALLOCATE(ff)

 end do     ! Loop on atoms

 DBG_EXIT("COLL")

 end subroutine make_orbl_paw
!!***

!!****f* m_paw_nmr/make_efg_onsite
!! NAME
!! make_efg_onsite
!!
!! FUNCTION
!! Compute the electric field gradient due to onsite densities
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nsym=number of symmetries in space group
!!  ntypat=number of atom types
!!  paw_an(my_natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  rprimd(3,3), conversion from crystal coordinates to cartesian coordinates
!!  symrel(3,3,nsym)=symmetry operators in terms of action on primitive translations
!!  tnons(3,nsym) = nonsymmorphic translations
!!  xred(3,natom), location of atoms in crystal coordinates.
!!
!! OUTPUT
!!  efg(3,3,natom), the 3x3 efg tensor at each site due to nhat

!! NOTES
!! This routine computes the electric field gradient, specifically the components
!! $\partial^2 V/\partial x_\alpha \partial x_\beta$ of the potential generated by the valence
!! electrons, at each atomic site in the unit cell. Key references: Kresse and Joubert, ``From
!! ultrasoft pseudopotentials to the projector augmented wave method'', Phys. Rev. B. 59, 1758--1775 (1999) [[cite:Kresse1999]],
!! and Profeta, Mauri, and Pickard, ``Accurate first principles prediction of $^{17}$O NMR parameters in
!! SiO$_2$: Assignment of the zeolite ferrierite spectrum'', J. Am. Chem. Soc. 125, 541--548 (2003) [[cite:Profeta2003]]. See in particular
!! Eq. 11 and 12 of Profeta et al., but note that their sum over occupied states times 2 for occupation number is
!! replaced in the Kresse and Joubert formulation by the sum over $\rho_{ij}$ occupations for each basis element pair.
!!
!! PARENTS
!!      calc_efg
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,matpointsym,pawdensities,pawrad_deducer0
!!      simp_gen,xmpi_sum
!!
!! SOURCE

subroutine make_efg_onsite(efg,my_natom,natom,nsym,ntypat,paw_an,pawang,pawrhoij,pawrad,pawtab, &
&                          rprimd,symrel,tnons,xred,&
&                          mpi_atmtab,comm_atom) ! optional arguments (parallelism)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,natom,nsym,ntypat
 integer,optional,intent(in) :: comm_atom
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: rprimd(3,3),tnons(3,nsym),xred(3,natom)
 real(dp),intent(out) :: efg(3,3,natom)
 type(paw_an_type),intent(in) :: paw_an(my_natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: cplex,iatom,iatom_tot,ictr,ierr,imesh_size,ispden,itypat
 integer :: lm,lm_size,lnspden,mesh_size,my_comm_atom,nzlmopt,nspden
 integer :: opt_compch,opt_dens,opt_l,opt_print
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: c1,c2,c3,compch_sph,intg
!arrays
 integer,pointer :: my_atmtab(:)
 logical,allocatable :: lmselectin(:),lmselectout(:)
 real(dp),allocatable :: ff(:),nhat1(:,:,:),rho1(:,:,:),trho1(:,:,:)

! ************************************************************************

 DBG_ENTER("COLL")

 if (my_natom>0) then
   ABI_CHECK(pawrhoij(1)%qphase==1,'make_efg_onsite: not supposed to be called with qphqse=2!')
 end if

 efg(:,:,:) = zero

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
& my_natom_ref=my_natom)

!the following factors arise in expanding the angular dependence of the electric field gradient tensor in
!terms of real spherical harmonics. The real spherical harmonics are as in the routine initylmr.F90; see
!in particular also http://www.unioviedo.es/qcg/art/Theochem419-19-ov-BF97-rotation-matrices.pdf
 c1 = sqrt(16.0*pi/5.0)
 c2 = sqrt(4.0*pi/5.0)
 c3 = sqrt(12.0*pi/5.0)

!loop over atoms in cell
 do iatom = 1, my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
   itypat=pawrhoij(iatom)%itypat

   lm_size = paw_an(iatom)%lm_size
   if (lm_size < 5) cycle ! if lm_size < 5 then the on-site densities for this atom have no L=2 component
!  and therefore nothing to contribute to the on-site electric field gradient

   mesh_size=pawtab(itypat)%mesh_size
   ABI_ALLOCATE(ff,(mesh_size))

   cplex = pawrhoij(iatom)%qphase
   nspden = pawrhoij(iatom)%nspden
   ABI_ALLOCATE(lmselectin,(lm_size))
   ABI_ALLOCATE(lmselectout,(lm_size))
   lmselectin = .true. ! compute all moments of densities
   nzlmopt = -1
   opt_compch = 0
   compch_sph = zero
   opt_dens = 0 ! compute all densities
   opt_l = -1 ! all moments contribute
   opt_print = 0 ! do not print out moments

   ABI_ALLOCATE(nhat1,(cplex*mesh_size,lm_size,nspden))
   ABI_ALLOCATE(rho1,(cplex*mesh_size,lm_size,nspden))
   ABI_ALLOCATE(trho1,(cplex*mesh_size,lm_size,nspden))

!  loop over spin components
!  nspden = 1: just a single component
!  nspden = 2: two components, loop over them
!  nspden = 4: total is in component 1, only one of interest
   if ( nspden == 2 ) then
     lnspden = 2
   else
     lnspden = 1
   end if
   do ispden=1,lnspden

!    construct multipole expansion of on-site charge densities for this atom
     call pawdensities(compch_sph,cplex,iatom_tot,lmselectin,lmselectout,lm_size,&
&     nhat1,nspden,nzlmopt,opt_compch,opt_dens,opt_l,opt_print,&
&     pawang,0,pawrad(itypat),pawrhoij(iatom),pawtab(itypat),&
&     rho1,trho1)

     do lm = 5, 9 ! loop on L=2 components of multipole expansion

       if(.not. lmselectout(lm)) cycle ! skip moments that contributes zero

!      the following is r^2*(n1-tn1-nhat)/r^3 for this multipole moment
!      use only the real part of the density in case of cplex==2
       do imesh_size = 2, mesh_size
         ictr = cplex*(imesh_size - 1) + 1
         ff(imesh_size)=(rho1(ictr,lm,ispden)-trho1(ictr,lm,ispden)-&
&         nhat1(ictr,lm,ispden))/&
&         pawrad(itypat)%rad(imesh_size)
       end do
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))
       select case (lm)
       case (5) ! S_{2,-2}
         efg(1,2,iatom_tot) = efg(1,2,iatom_tot) - c3*intg ! xy case
       case (6) ! S_{2,-1}
         efg(2,3,iatom_tot) = efg(2,3,iatom_tot) - c3*intg ! yz case
       case (7) ! S_{2, 0}
         efg(1,1,iatom_tot) = efg(1,1,iatom_tot) + c2*intg ! xx case
         efg(2,2,iatom_tot) = efg(2,2,iatom_tot) + c2*intg ! yy case
         efg(3,3,iatom_tot) = efg(3,3,iatom_tot) - c1*intg ! zz case
       case (8) ! S_{2,+1}
         efg(1,3,iatom_tot) = efg(1,3,iatom_tot) - c3*intg ! xz case
       case (9) ! S_{2,+2}
         efg(1,1,iatom_tot) = efg(1,1,iatom_tot) - c3*intg ! xx case
         efg(2,2,iatom_tot) = efg(2,2,iatom_tot) + c3*intg ! yy case
       end select

     end do  ! end loop over LM components with L=2

   end do    ! Loop on spin components

!  Symmetrization of EFG
   efg(2,1,iatom_tot) = efg(1,2,iatom_tot)
   efg(3,1,iatom_tot) = efg(1,3,iatom_tot)
   efg(3,2,iatom_tot) = efg(2,3,iatom_tot)

   ABI_DEALLOCATE(lmselectin)
   ABI_DEALLOCATE(lmselectout)
   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(nhat1)
   ABI_DEALLOCATE(rho1)
   ABI_DEALLOCATE(trho1)

 end do     ! Loop on atoms

!Reduction in case of parallelisation over atoms
 if (paral_atom) then
   call xmpi_sum(efg,my_comm_atom,ierr)
 end if

! symmetrize tensor at each atomic site using point symmetry operations
 do iatom = 1, natom
   call matpointsym(iatom,efg(:,:,iatom),natom,nsym,rprimd,symrel,tnons,xred)
 end do

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

 end subroutine make_efg_onsite
!!***

!----------------------------------------------------------------------

!!****f* m_paw_nmr/make_fc_paw
!! NAME
!! make_fc_paw
!!
!! FUNCTION
!! Compute the Fermi-contact term due to the PAW cores
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nspden=number of spin ensity component
!!  ntypat=number of atom types
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  fc(nspden,natom)=the Fermi-contact interaction at each site due to PAW for each spin density
!!
!! NOTES
!! The Fermi contact interaction is the electron density evaluated exactly at the nuclear site.
!! For a nuclear site at R, we are thus computing the expectation value of $\delta^3(R)$, the
!! the three-dimensional delta function at vector position $R$. In terms of the radial variable only
!! the delta function is $\delta(r)/4\pi r^2$.  Because this observable is
!! absolutely confined within the PAW radius, only the response due to the AE PAW functions is
!! needed, the pseudo wavefunctions and pseudo PAW functions cancel each other out. We then
!! must compute the integral of $u_i/r times u_j/r \delta(R)d^3r$, for the $l=0$ angular momentum
!! states only. This is simplified with the use of L'H\^{o}spital's theorem to take the limit
!! as $r\rightarrow 0$, yielding $u_i'(r) u_j'(r)$. To compute the derivatives we just fit the
!! first 5 points of the $u$ functions to a line through the origin, using the least squares
!! procedure resulting from $\chi = sum_i (y_i - m*x_i)^2$ . This is more stable than
!! computing the derivative of the whole function and extrapolating it to zero.
!! See Zwanziger, J. Phys. Conden. Matt. 21, 15024-15036 (2009) [[cite:Zwanziger2009]].
!!
!! PARENTS
!!      calc_fc
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,xmpi_sum
!!
!! SOURCE

subroutine make_fc_paw(fc,my_natom,natom,nspden,ntypat,pawrhoij,pawrad,pawtab,&
&                      mpi_atmtab,comm_atom) ! optional arguments (parallelism)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,natom,nspden,ntypat
 integer,optional,intent(in) :: comm_atom
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(out) :: fc(nspden,natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
 type(pawtab_type),target,intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iatom,iatom_tot,ierr,irhoij,islope,ispden,itypat
 integer :: ilmn,il,iln,ilm,im,jl,jlm,jlmn,jln,jm,j0lmn,jrhoij
 integer :: klmn,kln,mesh_size,my_comm_atom,nslope
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: mi,mj,xi,xxsum,xysumi,xysumj,yi,yj
!arrays
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 integer,pointer :: my_atmtab(:)

! ************************************************************************

 DBG_ENTER("COLL")

 if (my_natom>0) then
   ABI_CHECK(pawrhoij(1)%qphase==1,'make_fc_paw: not supposed to be called with qphqse=2!')
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Initialization
 fc(:,:)=zero

!number of points to use in computing initial slopes of radial functions
 nslope = 5

!loop over atoms in cell
 do iatom = 1, my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
   itypat = pawrhoij(iatom)%itypat
   mesh_size=pawtab(itypat)%mesh_size
   indlmn => pawtab(itypat)%indlmn

!  loop over spin components
   do ispden=1,nspden

!    loop over basis elements for this atom
!    ----
     do jlmn=1,pawtab(itypat)%lmn_size
       jl=indlmn(1,jlmn)
       jm=indlmn(2,jlmn)
       jlm=indlmn(4,jlmn)
       jln=indlmn(5,jlmn)
       j0lmn=jlmn*(jlmn-1)/2
       do ilmn=1,jlmn
         il=indlmn(1,ilmn)
         im=indlmn(2,ilmn)
         iln=indlmn(5,ilmn)
         ilm=indlmn(4,ilmn)
         klmn=j0lmn+ilmn
         kln = pawtab(itypat)%indklmn(2,klmn) ! need this for mesh selection below

         if (jl==0 .and. il==0) then ! select only s-states

!          Loop over non-zero elements of rhoij
           jrhoij=1
           do irhoij=1,pawrhoij(iatom)%nrhoijsel
             if (klmn==pawrhoij(iatom)%rhoijselect(irhoij)) then ! rho_ij /= 0 for this klmn
               xxsum = 0 ! these three variables will be used to compute the slopes
               xysumi = 0
               xysumj = 0
               do islope=1, nslope
                 xi=0
                 if(pawrad(itypat)%mesh_type == 1) xi = (islope - 1)*pawrad(itypat)%rstep
                 if(pawrad(itypat)%mesh_type == 2) xi = pawrad(itypat)%rstep * &
&                 (exp(pawrad(itypat)%lstep * (islope - 1)) - 1)
                 if(pawrad(itypat)%mesh_type == 3) then
                   if (islope == 1) then
                     xi = 0
                   else
                     xi = pawrad(itypat)%rstep * exp(pawrad(itypat)%lstep*(islope-1))
                   end if
                 end if
                 if(pawrad(itypat)%mesh_type == 4) xi = &
&                 -pawrad(itypat)%rstep*log(1.0-(islope-1)/pawrad(itypat)%mesh_size)
                 yi = pawtab(itypat)%phi(islope,iln) ! function value for u_i
                 yj = pawtab(itypat)%phi(islope,jln) ! function value for u_j
                 xxsum =  xxsum + xi*xi
                 xysumi = xysumi + xi*yi
                 xysumj = xysumj + xi*yj
               end do
!              the slopes of the radial functions are obtained by minimizing
!              chi = sum(y_i - m*x_i)^2 (in other words, a linear least squares
!              fit constrained to go through the origin)
!              the result is m = sum(y_i*x_i)/sum(x_i*x_i)
               mi = xysumi/xxsum
               mj = xysumj/xxsum
!              accumulate the rho_ij contribution to the fermi contact for this spin density:
               fc(ispden,iatom_tot)=fc(ispden,iatom_tot)+&
&                 pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(jrhoij,ispden)*mi*mj/four_pi
             end if ! end selection on klmn for nonzero rho_ij
             jrhoij=jrhoij+pawrhoij(iatom)%cplex_rhoij
           end do ! end loop over nonzero rho_ij
         end if ! end l=l'=0 selection
       end do ! end loop over ilmn
     end do ! end loop over jlmn
   end do ! end loop over spin densities
 end do     ! Loop on atoms

!Reduction in case of parallelisation over atoms
 if (paral_atom) then
   call xmpi_sum(fc,my_comm_atom,ierr)
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

 end subroutine make_fc_paw
!!***

!----------------------------------------------------------------------

END MODULE m_paw_nmr
!!***
