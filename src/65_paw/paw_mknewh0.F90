!!****f* ABINIT/paw_mknewh0
!! NAME
!! paw_mknewh0
!!
!! FUNCTION
!! Calculates the new bare PAW Hamiltonian in the case of quasi-particle self-consistent GW calculations.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2017 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nspden=number of spin-density components
!!  nfftf=(effective) number of FFT grid points (for this proc) for the "fine" grid
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated
!!  pawprtvol=control print volume and debugging output for PAW
!!  Cryst<crystal_t>=Info on unit cell and its symmetries
!!  Pawtab(ntypat*usepaw)<type(pawtab_type)>=paw tabulated starting data
!!  Paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  Pawang<type(pawang_type)>=paw angular mesh and related data
!!  Pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  vxc(nfftf,nspden)=exchange-correlation potential
!!  vxc_val(nfftf,nspden)=valence only exchange-correlation potential
!!  vtrial(nfftf,nspden)=potential (Hartree+XC+loc)
!!
!! SIDE EFFECTS
!!  Paw_ij(natom*usepaw)<Paw_ij_type)>=paw arrays given on (i,j) channels
!!     At output: new value for Paw_ij()%dij
!!
!! PARENTS
!!      calc_vhxc_me
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawgylm,symdij,symdij_all,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine paw_mknewh0(my_natom,nsppol,nspden,nfftf,pawspnorb,pawprtvol,Cryst,&
&          Pawtab,Paw_an,Paw_ij,Pawang,Pawfgrtab,vxc,vxc_val,vtrial,&
&          mpi_atmtab,comm_atom) ! optional arguments (parallelism)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi, only : xmpi_comm_self

 use m_crystal,      only : crystal_t
 use m_pawang,       only : pawang_type
 use m_pawtab,       only : pawtab_type
 use m_paw_an,       only : paw_an_type
 use m_paw_ij,       only : paw_ij_type
 use m_pawfgrtab,    only : pawfgrtab_type
 use m_pawdij,       only : symdij, symdij_all
 use m_paw_finegrid, only : pawgylm
 use m_paral_atom,   only : get_my_atmtab, free_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_mknewh0'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,nsppol,nspden,nfftf,pawprtvol,pawspnorb
 integer,optional,intent(in) :: comm_atom
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: vxc(nfftf,nspden),vxc_val(nfftf,nspden),vtrial(nfftf,nspden)
 type(crystal_t),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(Pawtab_type),target,intent(in) :: Pawtab(Cryst%ntypat)
 type(Paw_an_type),intent(in) :: Paw_an(my_natom)
 type(Paw_ij_type),intent(inout) :: Paw_ij(my_natom)
 type(Pawfgrtab_type),intent(inout) :: Pawfgrtab(my_natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: ipert0=0
 integer :: iat,iat_tot,idij,cplex,ndij,option_dij
 integer :: itypat,lmn_size,j0lmn,jlmn,ilmn,klmn,klmn1,klm
 integer :: lmin,lmax,mm,isel,lm_size,lmn2_size,my_comm_atom,cplex_dij
 integer :: ils,ilslm,ic,lm0
 integer :: nsploop,is2fft
 real(dp) :: gylm,qijl
 logical :: ltest,my_atmtab_allocated,paral_atom
 character(len=500) :: msg
!arrays
 integer,ABI_CONTIGUOUS pointer :: indklmn_(:,:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: rdum(1),rdum2(1)
 real(dp),allocatable :: prod_hloc(:,:),prodhxc_core(:,:)
 real(dp),allocatable :: dijhl_hat(:,:),dijhmxc_val(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call wrtout(std_out,'Assembling PAW strengths for the bare Hamiltonian','COLL')

!== Set up parallelism over atoms ===
 paral_atom=(present(comm_atom).and.(my_natom/=Cryst%natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,Cryst%natom,my_natom_ref=my_natom)

 if (my_natom>0) then

!  === Test if required pointers in paw_ij are allocated ===
   ltest = (allocated(Paw_ij(1)%dijxc).and.allocated(Paw_ij(1)%dijxc_val) )
   ABI_CHECK(ltest,'dijxc or dijxc_val not calculated')

   ltest=(allocated(Paw_ij(1)%dijhat)) !.and.Paw_ij(1)%has_dijhat==2)
   ABI_CHECK(ltest,'dijhat not calculated')

   ltest=(allocated(Paw_ij(1)%dijhartree)) !.and.Paw_ij(1)%has_dijhartree==2)
   ABI_CHECK(ltest,'dijhartree not calculated')

   if (ANY(Pawtab(:)%usepawu>0)) then
     do iat=1,my_natom
       iat_tot=iat;if (paral_atom) iat_tot=my_atmtab(iat)
       itypat=Cryst%typat(iat_tot)
       if (Pawtab(itypat)%usepawu>0) then
         ltest=(allocated(Paw_ij(iat)%dijU) ) !.and.Paw_ij(iat)%has_dijU==2)
         write(msg,'(a,i3,a)')" For atom no. ",iat," %dijU(iat) has not been calculated."
         ABI_CHECK(ltest,msg)
       end if
     end do
   end if

   if (pawspnorb>0) then
     do iat=1,my_natom
       ltest=(allocated(Paw_ij(iat)%dijso) ) !.and.Paw_ij(iat)%has_dijso==2)
       write(msg,'(a,i3,a)')" For atom no. ",iat," %dijso(iat) has not been calculated."
       ABI_CHECK(ltest,msg)
     end do
   end if
 end if ! my_natom>0

!== Construct the new PAW H0 Hamiltonian ===
 do iat=1,my_natom
   iat_tot=iat;if (paral_atom) iat_tot=my_atmtab(iat)

   itypat    = Cryst%typat(iat_tot)
   lmn_size  = Pawtab(itypat)%lmn_size
   lmn2_size = Pawtab(itypat)%lmn2_size
   lm_size   = Paw_an(iat)%lm_size
   cplex     = Paw_ij(iat)%cplex
   cplex_dij = Paw_ij(iat)%cplex_dij
   ndij      = Paw_ij(iat)%ndij

   ABI_CHECK(cplex==1,'cplex/=1 not implemented')
   ABI_CHECK(cplex_dij==1,'cplex_dij/=1 not implemented')
!  
!  Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
   if (Pawfgrtab(iat)%gylm_allocated==0) then
     if (allocated(Pawfgrtab(iat)%gylm))  then
       ABI_DEALLOCATE(Pawfgrtab(iat)%gylm)
     end if
     ABI_ALLOCATE(Pawfgrtab(iat)%gylm,(Pawfgrtab(iat)%nfgd,lm_size))
     Pawfgrtab(iat)%gylm_allocated=2

     call pawgylm(Pawfgrtab(iat)%gylm,rdum,rdum2,lm_size,&
&     Pawfgrtab(iat)%nfgd,1,0,0,Pawtab(itypat),Pawfgrtab(iat)%rfgd)
   end if

!  === Calculate LM contribution to dijhmxc_val for this atom ===
!  * Dijxc contains also the Hat term on the FFT mesh while Dijxc_val does not
!  contain neither the hat term nor the LM sum of onsite terms (they should cancel each other)
!  FIXME change paw_dij,  otherwise I miss tnc in vxc
!  * prodhxc_core is used to assemble $\int g_l Ylm (vtrial - vxc_val[tn+nhat] dr$ on the FFT mesh ===
!  * The following quantities do not depend on ij
   ABI_ALLOCATE(prod_hloc   ,(lm_size,ndij))
   ABI_ALLOCATE(prodhxc_core,(lm_size,ndij))
   prod_hloc   =zero
   prodhxc_core=zero
   do idij=1,ndij
     do ilslm=1,lm_size
       do ic=1,Pawfgrtab(iat)%nfgd
         is2fft=Pawfgrtab(iat)%ifftsph(ic)
         gylm=Pawfgrtab(iat)%gylm(ic,ilslm)
         prod_hloc (ilslm,idij)=prod_hloc (ilslm,idij) + (vtrial(is2fft,idij)-vxc(is2fft,idij))*gylm
!        prodhxc_core(ilslm,idij)=prodhxc_core(ilslm,idij) + (vxc_val(is2fft,idij))*gylm
         prodhxc_core(ilslm,idij)=prodhxc_core(ilslm,idij) + (vtrial(is2fft,idij)-vxc_val(is2fft,idij))*gylm
       end do
     end do
   end do !idij

!  === Assembly the "Hat" contribution for this atom ====
   ABI_ALLOCATE(dijhl_hat  ,(cplex_dij*lmn2_size,ndij))
   ABI_ALLOCATE(dijhmxc_val,(cplex_dij*lmn2_size,ndij))
   dijhl_hat  =zero
   dijhmxc_val=zero
   indklmn_ => Pawtab(itypat)%indklmn(1:6,1:lmn2_size)

   do idij=1,ndij
     do klmn=1,lmn2_size
       klm =indklmn_(1,klmn)
       lmin=indklmn_(3,klmn)
       lmax=indklmn_(4,klmn)

!      === $\sum_lm q_ij^l prod* for each idij$ ===
       do ils=lmin,lmax,2
         lm0=ils**2+ils+1
         do mm=-ils,ils
           ilslm=lm0+mm
           isel=Pawang%gntselect(lm0+mm,klm)
           if (isel>0) then
             qijl=Pawtab(itypat)%qijl(ilslm,klmn)
             dijhl_hat  (klmn,idij)=dijhl_hat  (klmn,idij) +  prod_hloc (ilslm,idij)*qijl
             dijhmxc_val(klmn,idij)=dijhmxc_val(klmn,idij) +prodhxc_core(ilslm,idij)*qijl
           end if
         end do
       end do
     end do
   end do

   ABI_DEALLOCATE(prod_hloc)
   ABI_DEALLOCATE(prodhxc_core)

!  * Normalization factor due to integration on the FFT mesh
   dijhl_hat  = dijhl_hat  *Cryst%ucvol/DBLE(nfftf)
   dijhmxc_val= dijhmxc_val*Cryst%ucvol/DBLE(nfftf)

!  === Now assembly the bare Hamiltonian ===
!  * Loop over density components overwriting %dij
   nsploop=nsppol; if (Paw_ij(iat)%ndij==4) nsploop=4

   do idij=1,nsploop
     klmn1=1

     do jlmn=1,lmn_size
       j0lmn=jlmn*(jlmn-1)/2
       do ilmn=1,jlmn
         klmn=j0lmn+ilmn

!        The following gives back the input dij.
!        since dijxc contains the hat term done on the FFT mesh
         if (.FALSE.) then
           Paw_ij(iat)%dij(klmn,idij) =        &
&           Pawtab(itypat)%dij0    (klmn)      &
&           +Paw_ij(iat)%dijhartree(klmn)      &
&           +Paw_ij(iat)%dijxc     (klmn,idij) &
&           +dijhl_hat   (klmn,idij)

         else
!          === Make nonlocal part of h0 removing the valence contribution ===
!          Remeber that XC contains already the Hat contribution
           Paw_ij(iat)%dij(klmn,idij) =        &
&           Pawtab(itypat)%dij0      (klmn)    &
&           +Paw_ij(iat)%dijhartree(klmn)      &
&           +Paw_ij(iat)%dijxc     (klmn,idij) &  ! 2 lines to get the d1-dt1 XC core contribution + XC hat (core+val)
&          -Paw_ij(iat)%dijxc_val (klmn,idij) &  ! I suppose that the "hat" term on the FFT mesh in included in both.
&          +dijhmxc_val(klmn,idij)               ! Local + Hartree - XC val contribution to the "hat" term.

!          Add the U contribution to the 
!          if (.FALSE. .and. Pawtab(itypat)%usepawu>0) then
           if (.TRUE. .and. Pawtab(itypat)%usepawu>0) then
             Paw_ij(iat)%dij(klmn,idij) = Paw_ij(iat)%dij(klmn,idij) + Paw_ij(iat)%dijU(klmn,idij)
           end if
         end if
!        TODO dijso, dijU, vpawx?
!        Just to be consistent, update some values.
!$Paw_ij(iat)%dijhat(klmn,idij)=Paw_ij(iat)%dijhat(klmn,idij)-dijhmxc_val(klmn,idij)

       end do !ilmn
     end do !jlmn
   end do !idij

!  this is to be consistent?
!  deallocate(Paw_ij(iat)%dijvxc_val)
   ABI_DEALLOCATE(dijhl_hat)
   ABI_DEALLOCATE(dijhmxc_val)
 end do !iat

!=== Symmetrize total Dij ===
 option_dij=0 ! For total Dij.
#if 0
 if (paral_atom) then
   call symdij(Cryst%gprimd,Cryst%indsym,ipert0,my_natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,option_dij,&
&   Paw_ij,Pawang,pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec,&
&   comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
 else
   call symdij(Cryst%gprimd,,Cryst%indsym,ipert0,my_natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,option_dij,&
&   Paw_ij,Pawang,pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec)
 end if
#else
 if (paral_atom) then
   call symdij_all(Cryst%gprimd,Cryst%indsym,ipert0,my_natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,&
&   Paw_ij,Pawang,pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec,&
&   comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
 else
   call symdij_all(Cryst%gprimd,Cryst%indsym,ipert0,my_natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,&
&   Paw_ij,Pawang,pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec)
 end if
#endif

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine paw_mknewh0
!!***
