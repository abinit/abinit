!{\src2tex{textfont=tt}}
!!****f* ABINIT/getcprjb
!! NAME
!! getcprjb
!!
!! FUNCTION
!!  Compute <Proj_i,k|Cn,k+b> for one wave function |Cn,k+b> expressed in reciprocal space.
!!  The projector is at  kpoint k, while the wavefunction is at kpoint k+b. This construction
!!  is needed for orbital magnetization, in the computation of the <u_n,k1|H_k2k2|u_m,k3>
!!  matrix elements.
!!  |Proj_i> are non-local projectors (for each atom and each l,m,n)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  cwavef(2,npwb)=wavefunction at k+b
!!  dtorbmag=orbital magnetization data structure
!!  idir=direction of b vector
!!  ifor=1 or 2 for +b or -b
!!  ikptb=index of kptb=k+b vector
!!  mpsang=max angular momentum across all PAW sets
!!  natom=number of atoms in cell
!!  npwb=number of planewaves around k+b
!!  nspinor=number of spinors
!!  ntypat=number of types of atoms in unit cell
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat=typat(natom) list of atom types
!!  ylm_kb(npwb,mpsang*mpsang)=real spherical harmonics at each plane wave
!!
!! OUTPUTS
!!  cwaveprj(natom,nspinor) <type(pawcprj_type)>=projected input wave function <Proj_i,k|Cn,k+b> with all NL projectors
!!
!! TODO
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

subroutine getcprjb(cwavef,cwaveprj,dtorbmag,idir,ifor,ikptb,&
     &mpsang,natom,npwb,nspinor,ntypat,pawang,pawrad,pawtab,typat,ylm_kb)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_orbmag

 use m_pawang,           only : pawang_type
 use m_pawcprj,  only : pawcprj_type
 use m_pawrad,           only : pawrad_type, simp_gen
 use m_pawtab, only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getcprjb'
!End of the abilint section

 implicit none

!Arguments -------------------------------
 !scalars
 integer,intent(in) :: idir,ifor,ikptb,mpsang,natom,npwb,nspinor,ntypat
 type(orbmag_type), intent(inout) :: dtorbmag
 type(pawang_type),intent(in) :: pawang
 !arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: cwavef(2,npwb),ylm_kb(npwb,mpsang*mpsang)
 type(pawcprj_type),intent(inout) ::  cwaveprj(natom,nspinor)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)


!Local variables-------------------------------
 !scalars
 integer :: iatom,idx,il,im,ilm,iln,ilmn,ipw,isel,ispinor,itypat
 integer :: kllt,ll,llmm,ltmt,lt,mm,mt,mesh_size
 real(dp) :: intg,sxpi2
 complex(dpc) :: c1,c2,cpc,etb,etk
 !arrays
 real(dp),allocatable :: ff(:)
 ! the following is (i)^L mod 4.
 complex(dpc),dimension(0:3) :: iexpl(0:3)=(/cone,j_dpc,-cone,-j_dpc/)


! *********************************************************************

 DBG_ENTER('COLL')

 idx = 2*(idir-1)+ifor
 
 sxpi2 = four_pi*four_pi

 ! this restriction will be removed in the future
 ispinor = 1 

 
 do iatom = 1, natom
   
   itypat = typat(iatom)
   mesh_size = pawtab(itypat)%mesh_size

   ABI_ALLOCATE(ff,(mesh_size))

    !    here is exp(-i b.R) for current atom
   etb = dtorbmag%expibi(idx,iatom)

   do ilmn=1,pawtab(itypat)%lmn_size
     il=pawtab(itypat)%indlmn(1,ilmn)
     im=pawtab(itypat)%indlmn(2,ilmn)
     ilm=pawtab(itypat)%indlmn(4,ilmn)
     iln=pawtab(itypat)%indlmn(5,ilmn)

     cpc = czero

     do ll=0,pawtab(itypat)%l_size-1
       do lt=0,pawtab(itypat)%l_size-1
             ! require |ll-lt| <= il <= ll+lt
         if ((il .GT. (ll+lt)) .OR. (il .LT. abs(ll-lt))) cycle
         if ( mod((il+ll+lt),2) .NE. 0 ) cycle

         do ipw = 1, npwb

                ! ff(1:mesh_size) = pawrad(itypat)%rad(1:mesh_size)*&
                !      &pawtab(itypat)%tproj(1:mesh_size,iln)*&
                !      &dtorbmag%jb_bessel(idir,itypat,1:mesh_size,ll+1)*&
                !      &dtorbmag%jkg_bessel(itypat,1:mesh_size,ikptb,ipw,lt+1)
                ! call simp_gen(intg,ff,pawrad(itypat))
           
           etk = dtorbmag%phkgi(ikptb,ipw,iatom)

           if (dtorbmag%has_pjj_integral(itypat,iln,ll+1,lt+1,idir,ikptb)) then
             intg = dtorbmag%pjj_integral(itypat,iln,ll+1,lt+1,idir,ikptb,ipw)
           else
             intg = zero
             write(std_out,'(a)')'JWZ Debug: no pjj'
           end if

           c1 = sxpi2*etb*etk*iexpl(mod(il,4))*iexpl(mod(lt,4))*intg

           do mm = -ll, ll
             llmm = ll*ll+ll+mm+1

             do mt = -lt, lt

               if ( (mm+mt+im) .NE. 0) cycle
               
               ltmt = lt*lt+lt+mt+1

               if (llmm .LE. ltmt) then
                 kllt = ltmt*(ltmt-1)/2 + llmm
               else
                 kllt = llmm*(llmm-1)/2 + ltmt
               end if

                      ! select on non-zero Gaunt integral G^{il,im}_{ll,mm,lt,mt}
               isel = pawang%gntselect(ilm,kllt)
               if (isel > 0) then

                 c2 = pawang%realgnt(isel)*dtorbmag%ylmb(idx,llmm)*ylm_kb(ipw,ltmt)
                 cpc = cpc + c1*c2*cmplx(cwavef(1,ipw),cwavef(2,ipw))

               end if ! selection on non-zero Gaunt integrals

             end do ! end loop on mt

           end do ! end loop on mm

         end do ! end loop on npwb

       end do ! end loop on lt

     end do ! end loop on ll
     
     cwaveprj(iatom,ispinor)%cp(1,ilmn) = real(cpc)            
     cwaveprj(iatom,ispinor)%cp(2,ilmn) = aimag(cpc)

   end do ! end loop on ilmn
   
   ABI_DEALLOCATE(ff)

 end do ! end loop over atoms

 DBG_EXIT('COLL')

 end subroutine getcprjb
!!***
