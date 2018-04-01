!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawdfptenergy
!! NAME
!! pawdfptenergy
!!
!! FUNCTION
!! This routine compute the Hartree+XC PAW on-site contributions to a 1st-order or 2nd-order energy.
!!  These contributions are equal to:
!!    E_onsite=
!!       Int{ VHxc[n1_a^(1);nc^(1)].n1_b }
!!      -Int{ VHxc[tild_n1_a^(j1)+hat_n1_a^(j1);tild_n_c^(j1)].(tild_n1_b+n1_b)^(j2) }
!! Some typical uses:
!!  A-Contribution to non-stationary expression of the 2nd-order total energy:
!!    In that case, n1_a^(1)[r]=n1^(j1)[r] and n1_b[r]=delta_n1^(j2)[r]
!!    where j1 and j2 are two given perturbations,
!!    and delta_n1^(j)[r] is the 1s-order density only due to change of WF overlap.
!!    See PRB 78, 035105 (2008), Eq.(80)
!!    E_onsite=
!!       Int{ VHxc[n1^(j1);nc^(j1)].delta_n1^(j2) }
!!      -Int{ VHxc[tild_n1^(j1)+hat_n1^(j1);tild_n_c^(j1)].delta_(tild_n1+hat_n1)^(j2) }
!!  B-Contribution to first-order Fermi energy:
!!    In that case, n1_a^(1)[r]=n1^(j1)[r] and n1_b[r]=n1[r,EFermi]
!!    where j1 is the current perturbation, and n1[r,EFermi] is the density at Fermi level.
!!    E_onsite=
!!       Int{ VHxc[n1^(j1);nc^(j1)].n1[r,EFermi] }
!!      -Int{ VHxc[tild_n1^(j1)+hat_n1^(j1);tild_n_c^(j1)].(tild_n1+hat_n1)[r,EFermi] }
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ipert1,ipert2=indexes of perturbations (j1) and (j2)
!!                if ipert2<=0, we compute a first-order energy
!!                if ipert2> 0, we compute a second-order energy
!!  ixc= choice of exchange-correlation scheme
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  ntypat=number of types of atoms in unit cell.
!!  nzlmopt_a= For the n1_a density:
!!            if -1, compute all LM-moments of the density and use non-zero LM-moments
!!            if  0, compute all LM-moments of the density and use all LM-moments
!!            if +1, compute only non-zero LM-moments of the density (stored before)
!!  nzlmopt_b= For the n1_b density:
!!            if -1, compute all LM-moments of the density and use non-zero LM-moments
!!            if  0, compute all LM-moments of the density and use all LM-moments
!!            if +1, compute only non-zero LM-moments of the density (stored before)
!!  paw_an0(natom) <type(paw_an_type)>=paw arrays for 0th-order quantities given on angular mesh
!!  paw_an1(natom) <type(paw_an_type)>=paw arrays for 1st-order quantities given on angular mesh
!!                                     This corresponds to (j1) perturbation
!!  paw_ij1(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!                                     This corresponds to (j1) perturbation
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij_a(natom) <type(pawrhoij_type)>= paw rhoij 1st-order occupancies for the (j1) perturbation
!!  pawrhoij_b(natom) <type(pawrhoij_type)>=
!!    if ipert2> 0: paw rhoij 1st-order occupancies for the (j2) perturbation corrsponding to n1_b^(j2)[r]
!!    if ipert2<=0: paw rhoij occupancies corresponding to n1_b[r]
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  xclevel= XC functional level
!!
!! OUTPUT
!!  delta_energy(2)= real and imaginary parts of contributions to non-stationary expression for the
!!              second derivative of the total energy
!!
!! SIDE EFFECTS
!!    ==== if paw_an1(:)%has_vxc<2, compute 1st-order XC potentials
!!      paw_an1(natom)%vxc1(cplex_a*mesh_size,:,nspden) =AE 1st-order XC potential Vxc^(j1)
!!      paw_an1(natom)%vxct1(cplex_a*mesh_size,:,nspden)=PS 1st-order XC potential tVxc^(j1)
!!    ==== if paw_ij1(:)%has_dijhartree<2, compute 1st-order Dij_hartree
!!      paw_ij1(natom)%dijhartree(cplex_a*lmn2_size)=Hartree contribution to Dij^(j1)
!!
!! PARENTS
!!      dfpt_nstpaw,newfermie1
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawdensities,pawdijhartree,pawxc_dfpt
!!      pawxcm_dfpt,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawdfptenergy(delta_energy,ipert1,ipert2,ixc,my_natom,natom,ntypat,nzlmopt_a,nzlmopt_b,&
&                    paw_an0,paw_an1,paw_ij1,pawang,pawprtvol,pawrad,pawrhoij_a,pawrhoij_b,&
&                    pawtab,pawxcdev,xclevel, &
&                    mpi_atmtab,comm_atom) ! optional arguments (parallelism)


 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi, only : xmpi_comm_self,xmpi_sum

 use m_time,       only : timab
 use m_pawang,     only : pawang_type
 use m_pawrad,     only : pawrad_type
 use m_pawtab,     only : pawtab_type
 use m_paw_an,     only : paw_an_type
 use m_paw_ij,     only : paw_ij_type
 use m_pawrhoij,   only : pawrhoij_type
 use m_pawdij,     only : pawdijhartree
 use m_pawxc,      only : pawxc_dfpt, pawxcm_dfpt
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdfptenergy'
 use interfaces_18_timing
 use interfaces_65_paw, except_this_one => pawdfptenergy
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ipert1,ipert2,ixc,my_natom,natom,ntypat,nzlmopt_a,nzlmopt_b
 integer,intent(in) :: pawprtvol,pawxcdev,xclevel
 integer,optional,intent(in) :: comm_atom
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(out) :: delta_energy(2)
 type(paw_an_type),intent(in) :: paw_an0(my_natom)
 type(paw_an_type),intent(inout) :: paw_an1(my_natom)
 type(paw_ij_type),intent(inout) :: paw_ij1(my_natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij_a(my_natom),pawrhoij_b(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: cplex_a,cplex_b,cplex_dijh1,iatom,iatom_tot,ierr,irhoij,ispden,itypat,jrhoij
 integer :: klmn,lm_size_a,lm_size_b,mesh_size,my_comm_atom,nspden,nspdiag,opt_compch,optexc,optvxc
 integer :: usecore,usetcore,usexcnhat
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: compch,eexc,eexc_im
 character(len=500) :: msg
!arrays
 integer,pointer :: my_atmtab(:)
 logical,allocatable :: lmselect_a(:),lmselect_b(:),lmselect_tmp(:)
 real(dp) :: dij(2),delta_energy_h(2),delta_energy_xc(2),ro(2),tsec(2)
 real(dp),allocatable :: kxc_dum(:,:,:),nhat1(:,:,:),rho1(:,:,:),trho1(:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(567,1,tsec)

 if (.not.(ipert1==natom+1.or.ipert1==natom+10.or.ipert1==natom+11 &
& .or.ipert2==natom+1.or.ipert2==natom+10.or.ipert2==natom+11)) then
   if((abs(nzlmopt_a)/=1.and.nzlmopt_a/=0).or.(abs(nzlmopt_b)/=1.and.nzlmopt_b/=0)) then
     msg='invalid value for nzlmopt !'
     MSG_BUG(msg)
   end if
   if (my_natom>0) then
     if(paw_ij1(1)%has_dijhartree==0) then
       msg='dijhartree must be allocated !'
       MSG_BUG(msg)
     end if
     if(paw_an1(1)%has_vxc==0) then
       msg='vxc1 and vxct1 must be allocated !'
       MSG_BUG(msg)
     end if
     if(paw_an0(1)%has_kxc==0) then
       msg='kxc1 must be allocated !'
       MSG_BUG(msg)
     end if
     if ((ipert1<=natom.or.ipert1==natom+1.or.ipert1==natom+10.or.ipert1==natom+11).and.paw_an0(1)%has_kxc/=2) then
       msg='XC kernels for ground state must be in memory !'
       MSG_BUG(msg)
     end if
     if (paw_ij1(1)%cplex/=paw_an1(1)%cplex) then
       msg='paw_ij1()%cplex and paw_an1()%cplex must be equal !'
       MSG_BUG(msg)
     end if
     if (pawrhoij_a(1)%cplex<paw_an1(1)%cplex.or.pawrhoij_b(1)%cplex<paw_an1(1)%cplex) then
       msg='pawrhoij()%cplex must be >=paw_an1()%cplex  !'
       MSG_BUG(msg)
     end if
     if (pawrhoij_a(1)%nspden/=pawrhoij_b(1)%nspden) then
       msg='pawrhoij_a()%nspden must =pawrhoij_b()%nspden  !'
       MSG_BUG(msg)
     end if
   end if
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Init contribution to 1st-order (or 2nd-order) energy
 delta_energy(1:2)=zero

!For some perturbations, nothing else to do
 if (ipert1==natom+1.or.ipert1==natom+10.or.ipert1==natom+11 .or. &
& ipert2==natom+1.or.ipert2==natom+10.or.ipert2==natom+11) return

!Various inits
 opt_compch=0;optvxc=1;optexc=3
 usecore=0;usetcore=0  ! This is true for phonons and Efield pert.
 usexcnhat=maxval(pawtab(1:ntypat)%usexcnhat)
 delta_energy_xc(1:2)=zero;delta_energy_h(1:2)=zero
 dij(1:2)=zero;ro(1:2)=zero


!================ Loop on atomic sites =======================
 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

   itypat=pawrhoij_a(iatom)%itypat
   mesh_size=pawtab(itypat)%mesh_size
   nspden=pawrhoij_a(iatom)%nspden
   cplex_a=pawrhoij_a(iatom)%cplex
   cplex_b=pawrhoij_b(iatom)%cplex
   cplex_dijh1=paw_ij1(iatom)%cplex
   lm_size_a=paw_an1(iatom)%lm_size
   if (ipert2<=0) lm_size_b=paw_an0(iatom)%lm_size
   if (ipert2> 0) lm_size_b=paw_an1(iatom)%lm_size

!  If Vxc potentials are not in memory, compute them
   if (paw_an1(iatom)%has_vxc/=2) then
     ABI_ALLOCATE(rho1 ,(cplex_a*mesh_size,lm_size_a,nspden))
     ABI_ALLOCATE(trho1,(cplex_a*mesh_size,lm_size_a,nspden))
     ABI_ALLOCATE(nhat1,(cplex_a*mesh_size,lm_size_a,nspden*usexcnhat))
     ABI_ALLOCATE(lmselect_a,(lm_size_a))
     lmselect_a(:)=paw_an1(iatom)%lmselect(:)
     ABI_ALLOCATE(lmselect_tmp,(lm_size_a))
     lmselect_tmp(:)=.true.
     if (nzlmopt_a==1) lmselect_tmp(:)=lmselect_a(:)
!    Compute on-site 1st-order densities
     call pawdensities(compch,cplex_a,iatom_tot,lmselect_tmp,lmselect_a,&
&     lm_size_a,nhat1,nspden,nzlmopt_a,opt_compch,1-usexcnhat,-1,0,pawang,pawprtvol,&
&     pawrad(itypat),pawrhoij_a(iatom),pawtab(itypat),rho1,trho1)
     ABI_DEALLOCATE(lmselect_tmp)
!    Compute on-site 1st-order xc potentials
     if (pawxcdev/=0) then
       call pawxcm_dfpt(pawtab(itypat)%coredens,cplex_a,cplex_a,eexc,ixc,paw_an0(iatom)%kxc1,&
&       lm_size_a,lmselect_a,nhat1,paw_an0(iatom)%nkxc1,mesh_size,nspden,optvxc,&
&       pawang,pawrad(itypat),rho1,usecore,0,&
&       paw_an1(iatom)%vxc1,xclevel)
       call pawxcm_dfpt(pawtab(itypat)%tcoredens(:,1),&
&       cplex_a,cplex_a,eexc,ixc,paw_an0(iatom)%kxct1,&
&       lm_size_a,lmselect_a,nhat1,paw_an0(iatom)%nkxc1,mesh_size,nspden,optvxc,&
&       pawang,pawrad(itypat),trho1,usetcore,2*usexcnhat,&
&       paw_an1(iatom)%vxct1,xclevel)
     else
       call pawxc_dfpt(pawtab(itypat)%coredens,cplex_a,cplex_a,eexc,ixc,paw_an0(iatom)%kxc1,&
&       lm_size_a,lmselect_a,nhat1,paw_an0(iatom)%nkxc1,mesh_size,nspden,optvxc,&
&       pawang,pawrad(itypat),rho1,usecore,0,&
&       paw_an0(iatom)%vxc1,paw_an1(iatom)%vxc1,xclevel)
       call pawxc_dfpt(pawtab(itypat)%tcoredens(:,1),&
&       cplex_a,cplex_a,eexc,ixc,paw_an0(iatom)%kxct1,&
&       lm_size_a,lmselect_a,nhat1,paw_an0(iatom)%nkxc1,mesh_size,nspden,optvxc,&
&       pawang,pawrad(itypat),trho1,usetcore,2*usexcnhat,&
&       paw_an0(iatom)%vxct1,paw_an1(iatom)%vxct1,xclevel)
     end if

     paw_an1(iatom)%has_vxc=2
     ABI_DEALLOCATE(lmselect_a)
     ABI_DEALLOCATE(rho1)
     ABI_DEALLOCATE(trho1)
     ABI_DEALLOCATE(nhat1)
   end if ! has_vxc

!  If Dij_hartree are not in memory, compute them
   if (paw_ij1(iatom)%has_dijhartree/=2) then
     call pawdijhartree(cplex_dijh1,paw_ij1(iatom)%dijhartree,paw_ij1(iatom)%nspden,&
&     pawrhoij_a(iatom),pawtab(itypat))
     paw_ij1(iatom)%has_dijhartree=2
   end if

!  Compute contribution to 1st-order (or 2nd-order) energy from 1st-order XC potential
   ABI_ALLOCATE(rho1 ,(cplex_b*mesh_size,lm_size_b,nspden))
   ABI_ALLOCATE(trho1,(cplex_b*mesh_size,lm_size_b,nspden))
   ABI_ALLOCATE(nhat1,(cplex_b*mesh_size,lm_size_b,nspden*usexcnhat))
   ABI_ALLOCATE(lmselect_b,(lm_size_b))
   if (ipert2<=0) lmselect_b(:)=paw_an0(iatom)%lmselect(:)
   if (ipert2> 0) lmselect_b(:)=paw_an1(iatom)%lmselect(:)
   ABI_ALLOCATE(lmselect_tmp,(lm_size_b))
   lmselect_tmp(:)=.true.
   if (nzlmopt_b==1) lmselect_tmp(:)=lmselect_b(:)
!  Compute on-site 1st-order densities
   call pawdensities(compch,cplex_b,iatom_tot,lmselect_tmp,lmselect_b,&
&   lm_size_b,nhat1,nspden,nzlmopt_b,opt_compch,1-usexcnhat,-1,0,pawang,pawprtvol,&
&   pawrad(itypat),pawrhoij_b(iatom),pawtab(itypat),rho1,trho1)
   ABI_DEALLOCATE(lmselect_tmp)
!  Compute contributions to 1st-order (or 2nd-order) energy
   if (pawxcdev/=0) then
     ABI_ALLOCATE(kxc_dum,(mesh_size,pawang%angl_size,0))
     call pawxcm_dfpt(pawtab(itypat)%coredens,cplex_b,cplex_a,eexc,ixc,kxc_dum,&
&     lm_size_b,lmselect_b,nhat1,0,mesh_size,nspden,optexc,pawang,pawrad(itypat),&
&     rho1,usecore,0,paw_an1(iatom)%vxc1,xclevel,d2enxc_im=eexc_im)

     delta_energy_xc(1)=delta_energy_xc(1)+eexc
     delta_energy_xc(2)=delta_energy_xc(2)+eexc_im
     call pawxcm_dfpt(pawtab(itypat)%tcoredens(:,1),&
&     cplex_b,cplex_a,eexc,ixc,kxc_dum,&
&     lm_size_b,lmselect_b,nhat1,0,mesh_size,nspden,optexc,pawang,pawrad(itypat),&
&     trho1,usetcore,2*usexcnhat,paw_an1(iatom)%vxct1,xclevel,&
&     d2enxc_im=eexc_im)
     ABI_DEALLOCATE(kxc_dum)
     delta_energy_xc(1)=delta_energy_xc(1)-eexc
     delta_energy_xc(2)=delta_energy_xc(2)-eexc_im
   else
     ABI_ALLOCATE(kxc_dum,(mesh_size,lm_size_b,0))
     call pawxc_dfpt(pawtab(itypat)%coredens,cplex_b,cplex_a,eexc,ixc,kxc_dum,&
&     lm_size_b,lmselect_b,nhat1,0,mesh_size,nspden,optexc,pawang,pawrad(itypat),&
&     rho1,usecore,0,paw_an0(iatom)%vxc1,paw_an1(iatom)%vxc1,xclevel,d2enxc_im=eexc_im)
     delta_energy_xc(1)=delta_energy_xc(1)+eexc
     delta_energy_xc(2)=delta_energy_xc(2)+eexc_im
     call pawxc_dfpt(pawtab(itypat)%tcoredens(:,1),&
&     cplex_b,cplex_a,eexc,ixc,kxc_dum,&
&     lm_size_b,lmselect_b,nhat1,0,mesh_size,nspden,optexc,pawang,pawrad(itypat),&
&     trho1,usetcore,2*usexcnhat,paw_an0(iatom)%vxct1,paw_an1(iatom)%vxct1,xclevel,&
&     d2enxc_im=eexc_im)
     ABI_DEALLOCATE(kxc_dum)
     delta_energy_xc(1)=delta_energy_xc(1)-eexc
     delta_energy_xc(2)=delta_energy_xc(2)-eexc_im
   end if
   ABI_DEALLOCATE(lmselect_b)
   ABI_DEALLOCATE(rho1)
   ABI_DEALLOCATE(trho1)
   ABI_DEALLOCATE(nhat1)

!  Compute contribution to 1st-order(or 2nd-order) energy from 1st-order Hartree potential
   nspdiag=1;if (nspden==2) nspdiag=2
   do ispden=1,nspdiag
     if (cplex_dijh1==1) then
       jrhoij=1
       do irhoij=1,pawrhoij_b(iatom)%nrhoijsel
         klmn=pawrhoij_b(iatom)%rhoijselect(irhoij)
         dij(1)=paw_ij1(iatom)%dijhartree(klmn)
         ro(1)=pawrhoij_b(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
         delta_energy_h(1)=delta_energy_h(1)+ro(1)*dij(1)
         if (cplex_b==2) then
           ro(2)=pawrhoij_b(iatom)%rhoijp(jrhoij+1,ispden)*pawtab(itypat)%dltij(klmn)
           delta_energy_h(2)=delta_energy_h(2)+ro(2)*dij(1)
         end if
         jrhoij=jrhoij+cplex_b
       end do
     else ! cplex_dijh1==2
       jrhoij=1
       do irhoij=1,pawrhoij_b(iatom)%nrhoijsel
         klmn=pawrhoij_b(iatom)%rhoijselect(irhoij)
         dij(1:2)=paw_ij1(iatom)%dijhartree(2*klmn-1:2*klmn)
         ro(1)=pawrhoij_b(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
         delta_energy_h(1)=delta_energy_h(1)+ro(1)*dij(1)
         delta_energy_h(2)=delta_energy_h(2)-ro(1)*dij(2)
         if (cplex_b==2) then
           ro(2)=pawrhoij_b(iatom)%rhoijp(jrhoij+1,ispden)*pawtab(itypat)%dltij(klmn)
           delta_energy_h(1)=delta_energy_h(1)+ro(2)*dij(2)
           delta_energy_h(2)=delta_energy_h(2)+ro(2)*dij(1)
         end if
         jrhoij=jrhoij+cplex_b
       end do
     end if
   end do

!  ================ End loop oon atomic sites =======================
 end do

!Final building of 1st-order (or 2nd-order) energy
 delta_energy(1:2)=delta_energy_xc(1:2)+delta_energy_h(1:2)

!Reduction in case of parallelism
 if (paral_atom) then
   call xmpi_sum(delta_energy,my_comm_atom,ierr)
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 call timab(567,2,tsec)

 DBG_EXIT("COLL")

end subroutine pawdfptenergy
!!***
