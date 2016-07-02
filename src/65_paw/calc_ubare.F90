!{\src2tex{textfont=tt}} 
!!****f* ABINIT/calc_ubare
!! NAME
!! calc_ubare
!!
!! FUNCTION
!! Calculate the bare interaction on atomic orbitals
!!
!! COPYRIGHT
!! Copyright (C) 1999-2015 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  itypatcor = value of itypat for correlated species
!!  lpawu = angular momentum for correlated species
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data:
!!  pawang
!!     %lmax=Maximum value of angular momentum l+1
!!     %gntselect((2*l_max-1)**2,l_max**2,l_max**2)=
!!                     selection rules for Gaunt coefficients
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data:
!!     %mesh_size=Dimension of radial mesh
!!     %rad(mesh_size)=The coordinates of all the points of the radial mesh
!!     %radfact(mesh_size)=Factor used to compute radial integrals
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      pawpuxinit
!!
!! CHILDREN
!!      poisson,simp_gen,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine calc_ubare(itypatcor,lpawu,pawang,pawrad,pawtab,rmax)

 use defs_basis
 use m_profiling_abi
 use m_xmpi
 use m_errors

 use m_pawang,  only : pawang_type, pawang_init, pawang_free
 use m_pawrad,  only : pawrad_type, simp_gen, nderiv_gen, poisson, pawrad_ifromr
 use m_pawtab,  only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_ubare'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none
!Arguments ------------------------------------
 integer, intent(in)   :: itypatcor,lpawu
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),target,intent(in) :: pawtab
 real(dp), optional, intent(in) :: rmax
 
!Local variables ------------------------------
!scalars
 integer :: ilmn,ilmn1,iln,iln1,isel,isel1,itypat,jlmn,jlmn1,jln,jln1
 integer :: klm,klm1,klmn,klmn1,ll,lm0
 integer :: lmin,lmin1,lmax,lmn2_size,mesh_size,meshsz,mm
 real(dp) :: norm,r_for_intg,rg,rg1,ubare,uint,uint_tmp
 character(len=800) :: message
!arrays
 real(dp),allocatable :: ff(:),gg(:),phiphj(:),phiphj1(:)

!************************************************************************

 itypat=itypatcor
 write(message,'(11a,f12.4,2a,i7,2a,f12.4,2a,i7,2a,f12.4)') &
& ch10," =======================================================================",ch10, &
& "  == Calculation of diagonal bare Coulomb interaction on ATOMIC orbitals ",ch10, &
& "     (it is assumed that the wavefunction for the first reference ",ch10, &
& "             energy in PAW atomic data is an atomic eigenvalue)",ch10,ch10, &
& " Max value of the radius in atomic data file   =", pawrad%rmax ,ch10, &
& " Max value of the mesh   in atomic data file   =", pawrad%mesh_size,ch10, &
& " PAW radius is                                 =", pawtab%rpaw,ch10, &
& " PAW value of the mesh for integration is      =", pawrad%int_meshsz,ch10, &
& " Integral of atomic wavefunction until rpaw    =", pawtab%ph0phiint(1)
 if(.not.present(rmax)) then
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 mesh_size=pawrad%mesh_size

!  Definition of the mesh used for integration.
 if(present(rmax)) then
   if(rmax>pawrad%rmax)  then 
     write(message, '(a)' ) 'calc_ubare: the radius cannot be larger than the maximum radius of the mesh'
     MSG_ERROR(message)
   end if
   meshsz=pawrad_ifromr(pawrad,rmax)+5
   r_for_intg=rmax
 else
   meshsz=pawtab%partialwave_mesh_size
   r_for_intg=pawrad%rad(meshsz)  ! (we could use r_for_intg=-1)
 end if

 lmn2_size=pawtab%lmn2_size
 ABI_ALLOCATE(ff,(mesh_size))
 ABI_ALLOCATE(gg,(mesh_size))
 ABI_ALLOCATE(phiphj,(mesh_size))
 ABI_ALLOCATE(phiphj1,(mesh_size))
 do klmn=1,lmn2_size
   ilmn=pawtab%indklmn(7,klmn);jlmn=pawtab%indklmn(8,klmn)
   ! Select lpawu and first projectors il=jl=lpawu and first proj only
   if (( pawtab%indklmn(3,klmn)+pawtab%indklmn(4,klmn)==2*lpawu).and. &  
&   (-pawtab%indklmn(3,klmn)+pawtab%indklmn(4,klmn)==2*lpawu).and. &  
&   (pawtab%indlmn(3,ilmn)==1).and.(pawtab%indlmn(3,jlmn)==1) ) then
     klm=pawtab%indklmn(1,klmn);iln=pawtab%indlmn(5,ilmn);jln=pawtab%indlmn(5,jlmn)
     lmin=pawtab%indklmn(3,klmn);lmax=pawtab%indklmn(4,klmn)
     phiphj(1:meshsz)=pawtab%phi(1:meshsz,iln)*pawtab%phi(1:meshsz,jln)
     !write(6,*) "A",klmn,pawtab%klmntomn(1,klmn),pawtab%klmntomn(2,klmn),&
     !&pawtab%indklmn(7,klmn),pawtab%indklmn(8,klmn),pawtab%klmntomn(3,klmn),pawtab%klmntomn(4,klmn)
     do ll=lmin,lmin,2
       lm0=ll*ll+ll+1
       ff(1:meshsz)=phiphj(1:meshsz)
       call simp_gen(norm,ff,pawrad,r_for_intg=r_for_intg)
       call poisson(ff,ll,pawrad,gg)
       do klmn1=klmn,lmn2_size
         ilmn1=pawtab%indklmn(7,klmn);jlmn1=pawtab%indklmn(8,klmn)
         ! Select lpawu and first projectors il=jl=lpawu and first proj only
         if (( pawtab%indklmn(3,klmn1)+pawtab%indklmn(4,klmn1)==2*lpawu).and. &
&         (-pawtab%indklmn(3,klmn1)+pawtab%indklmn(4,klmn1)==2*lpawu).and. &
&         (pawtab%indlmn(3,ilmn1)==1).and.(pawtab%indlmn(3,jlmn1)==1) ) then
           !write(6,*) "A1",klmn1,pawtab%klmntomn(1,klmn1),pawtab%klmntomn(2,klmn1),&
           !&pawtab%indklmn(7,klmn1),pawtab%indklmn(8,klmn1),pawtab%klmntomn(3,klmn1),pawtab%klmntomn(4,klmn1)
           klm1=pawtab%indklmn(1,klmn1);iln1=pawtab%indlmn(5,ilmn1);jln1=pawtab%indlmn(5,jlmn1)
           phiphj1(1:meshsz)=pawtab%phi(1:meshsz,iln1)*pawtab%phi(1:meshsz,jln1)
           uint_tmp=zero
           if ((ll==lmin1)) then
             ff(1)=zero
             ff(2:meshsz)=phiphj1(2:meshsz)*gg(2:meshsz)*two/pawrad%rad(2:meshsz)
             call simp_gen(uint_tmp,ff,pawrad,r_for_intg=r_for_intg)
           end if
           uint=zero
           do mm=-ll,ll
             isel =pawang%gntselect(lm0+mm,klm)
             isel1=pawang%gntselect(lm0+mm,klm1)
             if (isel>0.and.isel1>0) then
               rg =pawang%realgnt(isel)
               rg1=pawang%realgnt(isel1)
               uint=uint+uint_tmp*rg*rg1*two_pi
             end if
           end do
           if((pawtab%indklmn(5,klmn)==pawtab%indklmn(6,klmn)).and.&
&           (pawtab%indklmn(5,klmn1)==pawtab%indklmn(6,klmn1)).and.&
&           (pawtab%indklmn(5,klmn)==pawtab%indklmn(5,klmn1))) then
             ubare=uint*Ha_eV
           end if
         end if
       end do
     end do
   end if
 end do
 ABI_DEALLOCATE(gg)
 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(phiphj)
 ABI_DEALLOCATE(phiphj1)

 write(message,'(a,3(a,f12.4,a),2a,f12.4,a)') ch10," For an atomic wfn truncated at rmax =",r_for_intg,ch10,&
& "     The norm of the wfn is                    =",norm,ch10,&
& "     The bare interaction (no renormalization) =",ubare," eV",ch10,&
& "     The bare interaction (for a renorm. wfn ) =",ubare/norm/norm," eV"
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 if(r_for_intg < 10_dp .and. .not.present(rmax)) then
   write(message,'(a,f6.2,4a)') '   ( WARNING: The radial mesh in the atomic data file is cut at',r_for_intg,ch10,&
&   '   Use XML atomic data files to compute the bare Coulomb interaction',ch10,&
&   '   on a true normalized atomic wavefunction )'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if
 if(present(rmax)) then
   write(message,'(2a)')  " =======================================================================",ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 end subroutine calc_ubare
!!***
