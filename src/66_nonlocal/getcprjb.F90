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
!!  bvec(3)=offset vector appearing in k+b in terms of recip. translations
!!  gprimd(3,3)=reciprocal lattice vectors
!!  kbg(3,npwb)=plane wave coordinates around k+b
!!  kpointb(3)=k+b point in terms of recip. translation
!!  natom=number of atoms in cell
!!  npwb=number of planewaves around k+b
!!  ntypat=number of types of atoms in unit cell
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat=typat(natom) list of atom types
!!  xred(3,natom)=atomic locations, reduced coordinates
!!  ylm_kb(npwb,psps%mpsang*psps%mpsang)=real spherical harmonics at each plane wave
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

subroutine getcprjb(bvec,gprimd,kbg,kpointb,natom,npwb,ntypat,&
     & pawang,pawrad,pawtab,psps,typat,xred,ylm_kb)

 use defs_basis
 use defs_abitypes
 use defs_datatypes, only : pseudopotential_type
 use m_profiling_abi
 use m_errors

 use m_special_funcs, only : sbf8
 use m_pawang, only : pawang_type
 use m_pawrad, only : pawrad_type, simp_gen
 use m_pawtab, only : pawtab_type
 use m_paw_sphharm, only : initylmr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getcprjb'
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments -------------------------------
 !scalars
 integer,intent(in) :: natom,npwb,ntypat
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 !arrays
 integer,intent(in) :: kbg(3,npwb),typat(natom)
 real(dp),intent(in) :: bvec(3),gprimd(3,3),kpointb(3),xred(3,natom)
 real(dp),intent(in) :: ylm_kb(npwb,psps%mpsang*psps%mpsang*psps%useylm)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
 !scalars
 integer :: iatom,ipw,ir,itypat,mesh_size,ylmr_normchoice,ylmr_npts,ylmr_option
 real(dp) :: arg,bnorm,kpgnorm,phfac,sxpi2
 complex(dpc) :: etb
 !arrays
 real(dp) :: bb(3),bbn(3),bcart(3),kpgcart(3),kpgvec(3),ylmgr(1,1,0),ylmr_nrm(1)
 real(dp),allocatable :: calc_expibi(:,:),jb_bessel(:,:),jkg_bessel(:,:,:),ylmb(:),sb_out(:)
 complex(dpc),allocatable :: phkgi(:,:)
 ! the following is (i)^L mod 4.
 complex(dpc),dimension(0:3) :: il(0:3)=(/cone,j_dpc,-cone,-j_dpc/)


! *********************************************************************

 DBG_ENTER('COLL')

 sxpi2 = four_pi*four_pi

 ylmr_normchoice = 0 ! input to initylmr are normalized
 ylmr_npts = 1 ! only 1 point to compute in initylmr
 ylmr_nrm(1) = one ! weight of normed point for initylmr
 ylmr_option = 1 ! compute only ylm's in initylmr
 ABI_ALLOCATE(sb_out, (pawang%l_size_max))

 ! compute exp(-i b*I) at each atomic position I
 ABI_ALLOCATE(calc_expibi,(2,natom))
 call expibi(calc_expibi,bvec,natom,xred)

 ABI_ALLOCATE(phkgi,(npwb,natom))
 do ipw = 1, npwb
    kpgvec(:) = kbg(:,ipw)+kpointb(:)
    ! compute phase factors exp(i*(kb+Gb)*I) for each plane wave at each position
    do iatom = 1, natom
       phfac=two_pi*DOT_PRODUCT(kpgvec(:),xred(:,iatom))
       phkgi(ipw,iatom)=cmplx(cos(phfac),sin(phfac),kind=dpc)
    end do
 end do

 do iatom = 1, natom
    
    itypat = typat(iatom)
    mesh_size = pawtab(itypat)%mesh_size

    ABI_ALLOCATE(jb_bessel,(mesh_size,pawang%l_size_max))
    ABI_ALLOCATE(jkg_bessel,(mesh_size,npwb,pawang%l_size_max))
    ABI_ALLOCATE(ylmb,(pawang%l_size_max*pawang%l_size_max))

    !    here is exp(-i b.R) for current atom
    etb = cmplx(calc_expibi(1,iatom),calc_expibi(2,iatom))

    bb(:) = -bvec(:)

    !    reference bb to cartesian axes
    bcart(1:3)=MATMUL(gprimd(1:3,1:3),bb(1:3))

    !    bbn is b-hat (the unit vector in the b direction)
    bnorm=dsqrt(dot_product(bcart,bcart))
    bbn(:) = bcart(:)/bnorm

    !    as an argument to the bessel function, need 2pi*b*r = 1 so b is re-normed to two_pi
    bnorm = two_pi*bnorm
    do ir=1,mesh_size
       arg=bnorm*pawrad(itypat)%rad(ir)
       call sbf8(pawang%l_size_max,arg,sb_out) ! spherical bessel functions at each mesh point
       jb_bessel(ir,:) = sb_out
       do ipw=1,npwb
          kpgvec(:) = kbg(:,ipw)+kpointb(:)
          ! reference to cartesian axes
          kpgcart(1:3)=MATMUL(gprimd(1:3,1:3),kpgvec(1:3))
          kpgnorm=dsqrt(dot_product(kpgcart,kpgcart))
          kpgnorm=two_pi*kpgnorm
          arg=kpgnorm*pawrad(itypat)%rad(ir)
          call sbf8(pawang%l_size_max,arg,sb_out) ! spherical bessel functions at each mesh point
          jkg_bessel(ir,ipw,:)=sb_out
       end do
    end do ! end loop over mesh

    !    compute Y_LM(b) here
    call initylmr(pawang%l_size_max,ylmr_normchoice,ylmr_npts,ylmr_nrm,ylmr_option,bbn,ylmb(:),ylmgr)

    ABI_DEALLOCATE(jb_bessel)
    ABI_DEALLOCATE(jkg_bessel)
    ABI_DEALLOCATE(ylmb)
 end do ! end loop over atoms

 ABI_DEALLOCATE(sb_out)
 ABI_DEALLOCATE(calc_expibi)
 ABI_DEALLOCATE(phkgi)
      
 DBG_EXIT('COLL')

 end subroutine getcprjb
!!***
