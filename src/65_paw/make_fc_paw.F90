!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_fc_paw
!! NAME
!! make_fc_paw
!!
!! FUNCTION
!! Compute the Fermi-contact term due to the PAW cores
!!
!! COPYRIGHT
!! Copyright (C) 2005-2017 ABINIT group (JZ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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
!! See Zwanziger, J. Phys. Conden. Matt. 21, 15024-15036 (2009).
!!
!! PARENTS
!!      calc_fc
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine make_fc_paw(fc,my_natom,natom,nspden,ntypat,pawrhoij,pawrad,pawtab,&
&                      mpi_atmtab,comm_atom) ! optional arguments (parallelism)


 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi, only : xmpi_comm_self

 use m_paral_atom, only : get_my_atmtab, free_my_atmtab
 use m_pawrad,     only : pawrad_type
 use m_pawtab,     only : pawtab_type
 use m_pawrhoij,   only : pawrhoij_type
 use m_xmpi,       only : xmpi_sum

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_fc_paw'
!End of the abilint section

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
 integer :: ilmn,il,iln,ilm,im,jl,jlm,jlmn,jln,jm,j0lmn
 integer :: klmn,kln,mesh_size,my_comm_atom,nslope
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: mi,mj,xi,xxsum,xysumi,xysumj,yi,yj
!arrays
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 integer,pointer :: my_atmtab(:)

! ************************************************************************

 DBG_ENTER("COLL")

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
               if (pawrhoij(iatom)%cplex == 1) then
                 fc(ispden,iatom_tot)=fc(ispden,iatom_tot)+&
&                 pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(irhoij,ispden)*mi*mj/four_pi
               else
                 fc(ispden,iatom_tot)=fc(ispden,iatom_tot)+&
&                 pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(2*irhoij-1,ispden)*mi*mj/four_pi
               end if
             end if ! end selection on klmn for nonzero rho_ij
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
