!{\src2tex{textfont=tt}}
!!****f* ABINIT/initrhoij
!! NAME
!! initrhoij
!!
!! FUNCTION
!! Initialize PAW rhoij occupancies (in packed storage)
!! from atomic ones
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cplex=1 if rhoij are REAL, 2 if they are complex
!!  lexexch(ntypat)=l on which local exact-exchange is applied for a given type of atom
!!  lpawu(ntypat)=l on which U is applied for a given type of atom (PAW+U)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of atom types
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated in PAW augmentation regions
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!                                     (containing initial rhoij)
!!  spinat(3,natom)=initial spin of each atom, in unit of hbar/2.
!!  typat(natom)=type of each atom
!!  === Optional arguments
!!    ngrhoij=number of gradients to be allocated (OPTIONAL, default=0)
!!    nlmnmix=number of rhoij elements to be mixed during SCF cycle (OPTIONAL, default=0)
!!    use_rhoij_=1 if pawrhoij(:)%rhoij_ has to be allocated (OPTIONAL, default=0)
!!    use_rhoijres=1 if pawrhoij(:)%rhoijres has to be allocated (OPTIONAL, default=0)

!!
!! OUTPUT
!!  pawrhoij(natom) <type(pawrhoij_type)>=rhoij quantities for each atom
!!                                        in packed storage
!!
!! PARENTS
!!      gstate,respfn,setup_positron
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawrhoij_alloc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initrhoij(cplex,lexexch,lpawu,my_natom,natom,&
&                    nspden,nspinor,nsppol,ntypat,pawrhoij,pawspnorb,pawtab,spinat,typat,&
&                    ngrhoij,nlmnmix,use_rhoij_,use_rhoijres,& ! optional arguments
&                    mpi_atmtab,comm_atom) ! optional arguments (parallelism)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi, only : xmpi_comm_self

 use m_pawtab,      only : pawtab_type
 use m_pawrhoij,    only : pawrhoij_type, pawrhoij_alloc, pawrhoij_get_nspden
 use m_paral_atom,  only : get_my_atmtab, free_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initrhoij'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,my_natom,natom,nspden,nspinor,nsppol,ntypat,pawspnorb
 integer,intent(in),optional :: comm_atom,ngrhoij,nlmnmix,use_rhoij_,use_rhoijres
 character(len=500) :: message
!arrays
 integer,intent(in) :: lexexch(ntypat),lpawu(ntypat)
 integer,intent(in) :: typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: spinat(3,natom)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!Arrays
!scalars
 integer :: iatom,iatom_rhoij,ilmn,ispden,itypat,j0lmn,jl,jlmn,jspden,klmn,klmn1,my_comm_atom
 integer :: ngrhoij0,nlmnmix0,nselect,nselect1,nspden_rhoij,use_rhoij_0,use_rhoijres0
 real(dp) :: ratio,ro,roshift,zratio,zz
 logical :: my_atmtab_allocated,paral_atom,test_exexch,test_pawu
!arrays
 integer,pointer :: my_atmtab(:)

!************************************************************************

 DBG_ENTER("COLL")

!PAW+U and local exact-exchange restriction
 do itypat=1,ntypat
   if (lpawu(itypat)/=lexexch(itypat).and. lpawu(itypat)/=-1.and.lexexch(itypat)/=-1) then
     message = ' lpawu must be equal to lexexch !'
     MSG_ERROR(message)
   end if
 end do

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 nspden_rhoij=pawrhoij_get_nspden(nspden,nspinor,pawspnorb)
 ratio=one;if (nspden_rhoij==2) ratio=half

 if (my_natom>0) then
   ngrhoij0=0;if (present(ngrhoij)) ngrhoij0=ngrhoij
   nlmnmix0=0;if (present(nlmnmix)) nlmnmix0=nlmnmix
   use_rhoij_0=0;if (present(use_rhoij_)) use_rhoij_0=use_rhoij_
   use_rhoijres0=0;if (present(use_rhoijres)) use_rhoijres0=use_rhoijres
   if (paral_atom) then
     call pawrhoij_alloc(pawrhoij,cplex,nspden_rhoij,nspinor,nsppol,typat,&
&     ngrhoij=ngrhoij0,nlmnmix=nlmnmix0,use_rhoij_=use_rhoij_0,use_rhoijres=use_rhoijres0,&
&     pawtab=pawtab,comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
   else
     call pawrhoij_alloc(pawrhoij,cplex,nspden_rhoij,nspinor,nsppol,typat,pawtab=pawtab,&
&     ngrhoij=ngrhoij0,nlmnmix=nlmnmix0,use_rhoij_=use_rhoij_0,use_rhoijres=use_rhoijres0)
   end if
 end if

 do iatom_rhoij=1,my_natom
   iatom=iatom_rhoij;if (paral_atom) iatom=my_atmtab(iatom_rhoij)
   itypat=typat(iatom)
   nselect=0

!  Determine Z (trace of rhoij0 or part of it)
   zz=zero
   do jlmn=1,pawtab(itypat)%lmn_size
     jl=pawtab(itypat)%indlmn(1,jlmn)
     j0lmn=jlmn*(jlmn-1)/2
     test_pawu=(lpawu(itypat)==-1.or.lpawu(itypat)==jl)
     test_exexch=(lexexch(itypat)==-1.or.lexexch(itypat)==jl)
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       if ((ilmn==jlmn).and.test_pawu.and.test_exexch) &
&       zz=zz+pawtab(itypat)%rhoij0(klmn)
     end do
   end do

!  Compute rhoij from tabulated value and magnetization
   do ispden=1,nspden_rhoij

     zratio=zero
     roshift=one
     ratio=one
     if (nspden_rhoij==2) then
       ratio=half
       if ((spinat(3,iatom)>zero.and.ispden==1).or.&
&       (spinat(3,iatom)<zero.and.ispden==2)) then
         if(abs(zz)>tol12)then
           zratio=two*abs(spinat(3,iatom))/zz
         else
           zratio=zero
         end if
       end if
     else if (nspden_rhoij==4.and.ispden>=2) then
       roshift=zero
       if(abs(zz)>tol12)then
         zratio=spinat(ispden-1,iatom)/zz
       else
         zratio=zero
       end if
     end if

     nselect=0;nselect1=1-cplex
     do jlmn=1,pawtab(itypat)%lmn_size
       jl=pawtab(itypat)%indlmn(1,jlmn)
       j0lmn=jlmn*(jlmn-1)/2
       test_pawu=(lpawu(itypat)==-1.or.lpawu(itypat)==jl)
       test_exexch=(lexexch(itypat)==-1.or.lexexch(itypat)==jl)
       do ilmn=1,jlmn
         klmn=j0lmn+ilmn
         ro=pawtab(itypat)%rhoij0(klmn)
         if ((ilmn==jlmn).and.test_pawu.and.test_exexch) then
           ro=ro*ratio*(roshift+zratio)
         else
           ro=ro*ratio*roshift
         end if

         klmn1=cplex*(klmn-1)+1
         if (abs(ro)>tol10) then
           pawrhoij(iatom_rhoij)%rhoijp(klmn1,ispden)=ro
         else
           pawrhoij(iatom_rhoij)%rhoijp(klmn1,ispden)=zero
         end if

         if (ispden==nspden_rhoij) then
           if (any(abs(pawrhoij(iatom_rhoij)%rhoijp(klmn1,:))>tol10)) then
             nselect=nselect+1;nselect1=nselect1+cplex
             pawrhoij(iatom_rhoij)%rhoijselect(nselect)=klmn
             do jspden=1,nspden_rhoij
               pawrhoij(iatom_rhoij)%rhoijp(nselect1,jspden)=pawrhoij(iatom_rhoij)%rhoijp(klmn1,jspden)
             end do
           end if
         end if

       end do
     end do

   end do
   pawrhoij(iatom_rhoij)%nrhoijsel=nselect

 end do ! iatom_rhoij

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine initrhoij
!!***
