!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_orbmag
!! NAME
!!  m_orbmag
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  used to handle orbital magnetization
!!
!! COPYRIGHT
!! Copyright (C) 2011-2017 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
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

module m_orbmag

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_pawcprj, only : pawcprj_type, pawcprj_free

 implicit none

 private
!!***


!!****t* defs_datatypes/orbmag_type
!! NAME
!! orbmag_type
!!
!! FUNCTION
!! variables used in orbital magnetism calculation
!!
!! SOURCE

 type, public :: orbmag_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer variables
  integer :: orbmag              ! value of orbmag input variable in use
  integer :: fmkmem              ! number of k-points in the FBZ per cpu
  integer :: fmkmem_max          ! max of fmkmem
  integer :: fnkpt               ! number of k-points in the FBZ
  integer :: lmax
  integer :: lmnmax
  integer :: lmn2max
  integer :: mkmem_max           ! max of mkmem
  integer :: natom               ! number of atoms in unit cell
  integer :: my_natom            ! number of atoms treated by current proc
  integer :: mband_occ           ! max number of occupied bands (over spin)
                                 ! this number must be the same for every k
  integer :: nspinor             ! nspinor input from data set
  integer :: nsym
  integer :: usepaw              ! 1 if a PAW calculation, 0 else

! Real(dp) scalars
  real(dp) :: sdeg               ! spin degeneracy: sdeg = 2 if nsppol = 1

  ! Real(dp) arrays
  real(dp) :: chern(2,3)           ! result of chern number calculation
  
  real(dp) :: dkvecs(3,3)        ! dkvec(:,idir) = vector between a k-point and its nearest neighbour along idir

  real(dp) :: orbmagvec(2,3)     ! result of orbital magnetization calculation

  ! Integer pointers
  integer, allocatable :: atom_indsym(:,:,:) ! atom_indsym(4,nsym,natom)
                                         ! this is data on how the symmetries map the atoms in the cell
                                         ! see symatm.F90 for full description
  integer, allocatable :: cgindex(:,:)    ! cgindex(nkpt,nsppol)
                                      ! for each k-point, stores the location
                                      ! of the WF in the cg array
  integer, allocatable :: cprjindex(:,:)  ! cprjindex(nkpt,nsppol)
                                      ! for each k-point, stores the location
                                      ! of the cprj in the cprj array (used only
                                      ! for PAW calculations)
  integer, allocatable :: fkgindex(:)     ! same as kgindex, but defined
                                      ! for the FBZ and intended to use
                                      ! with pwindf
  integer, allocatable :: ikpt_dk(:,:,:)  ! ikpt_dk(nkpt,2,3)
                                      ! ikpt_dp(ikpt,ii,idir) = index of the
                                      ! k-point at k+dk (ii=1) and k-dk (ii=2)
  integer, allocatable :: indkk_f2ibz(:,:)   ! indkk_f2ibz(1:dtorbmag%fnkpt,1:6)
                                         ! information needed to fold a
                                         ! k-point in the FBZ into the IBZ;
                                         ! the second index (1:6)
                                         ! is as described in listkk
  integer, allocatable :: i2fbz(:)           ! i2fbz(1:nkpt) gives index of IBZ
                                         ! k-points in the FBZ k-point list

  integer, allocatable :: kgindex(:)      ! kgind(nkpt)
                                      ! kgind(ikpt) = ikg

  integer, allocatable :: lmn_size(:)        ! lmn_size(ntypat)
  integer, allocatable :: lmn2_size(:)       ! lmn2_size(ntypat)

  integer, allocatable :: nband_occ(:)       ! nband_occ(nsppol) = actual number of occupied bands
                                             !  can be different for spin up and down!!!
! Real(dp) allocatables

  real(dp), allocatable :: fkptns(:,:)       ! fkptns(3,1:dtorbmag%fnkpt) k-points in FBZ

  real(dp),allocatable :: twdij0(:,:,:,:)    ! twdij0(2,24,lmn2max,natom) k1/k2/k3 twisted Dij0 terms
  
  real(dp), allocatable :: zarot(:,:,:,:)
   !  zarot(l_size_max,l_size_max,l_max,nsym)
   !  Coeffs of the transformation of real spherical
   !  harmonics under the symmetry operations. These are needed when the
   ! cprj's need to be computed in the full BZ, that is,
   ! in the PAW case with kptopt /= 3.

 end type orbmag_type

 ! Bound methods:
 public :: destroy_orbmag
 public :: pawtwdij_2b
!!***

contains

!!****f* m_orbmag/destroy_orbmag
!! NAME
!!
!! FUNCTION
!!   deallocate fields in orbmag structure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_orbmag(dtorbmag)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_orbmag'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(orbmag_type),intent(inout) :: dtorbmag

! ************************************************************************

! Integer pointers
  if(allocated(dtorbmag%atom_indsym))  then
    ABI_DEALLOCATE(dtorbmag%atom_indsym)
  end if
  if(allocated(dtorbmag%cgindex))  then
    ABI_DEALLOCATE(dtorbmag%cgindex)
  end if
  if(allocated(dtorbmag%cprjindex))  then
    ABI_DEALLOCATE(dtorbmag%cprjindex)
  end if
  if(allocated(dtorbmag%fkgindex))  then
    ABI_DEALLOCATE(dtorbmag%fkgindex)
  end if
  if(allocated(dtorbmag%ikpt_dk))  then
    ABI_DEALLOCATE(dtorbmag%ikpt_dk)
  end if
  if(allocated(dtorbmag%indkk_f2ibz))  then
    ABI_DEALLOCATE(dtorbmag%indkk_f2ibz)
  end if
  if(allocated(dtorbmag%i2fbz))  then
    ABI_DEALLOCATE(dtorbmag%i2fbz)
  end if
  if(allocated(dtorbmag%kgindex))  then
    ABI_DEALLOCATE(dtorbmag%kgindex)
  end if
  if(allocated(dtorbmag%lmn_size))  then
    ABI_DEALLOCATE(dtorbmag%lmn_size)
  end if
  if(allocated(dtorbmag%lmn2_size))  then
    ABI_DEALLOCATE(dtorbmag%lmn2_size)
  end if
  if(allocated(dtorbmag%nband_occ))  then
    ABI_DEALLOCATE(dtorbmag%nband_occ)
  end if
! Real(dp) pointers

  if(allocated(dtorbmag%fkptns))  then
    ABI_DEALLOCATE(dtorbmag%fkptns)
  end if
  if(allocated(dtorbmag%twdij0)) then
     ABI_DEALLOCATE(dtorbmag%twdij0)
  end if
  if(allocated(dtorbmag%zarot))  then
    ABI_DEALLOCATE(dtorbmag%zarot)
  end if

end subroutine destroy_orbmag
!!***

!----------------------------------------------------------------------

!!****f* m_orbmag/pawtwdij_2b
!! NAME
!! pawtwdij_2b
!!
!! FUNCTION
!! compute phase-twisted PAW on-site term 2b, that is,
!! VH[Znc] - VH[tZnc] in PAW spheres with a k point difference
!!
!! COPYRIGHT
!! Copyright (C) 2005-2018 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dkvecs(3) :: $\Delta k$ input vector
!!  expibi(2,my_natom,3) :: phase factors at each atomic site for given k offset
!!  gprimd(3,3)=dimensioned primitive translations of reciprocal lattice
!!  lmn2max :: lmnmax*(lmnmax+1)/2
!!  natom=number of atoms in unit cell
!!  ntypat=number of types of atoms in unit cell
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat=typat(natom) list of atom types
!!
!! OUTPUT
!!  calc_qijb(2,lmn2max,natom) :: PAW on-site overlaps of wavefunctions at neighboring
!!                                   k point
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      initylmr,sbf8,simp_gen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawtwdij_2b(calc_twdij,dtorbmag,gprimd,ntypat,pawang,pawrad,pawtab,typat,xred)

 use m_profiling_abi

 use defs_basis
 use m_errors
 use m_xmpi, only : xmpi_sum
 use m_pawang, only : pawang_type
 use m_pawrad, only : pawrad_type, simp_gen
 use m_pawtab, only : pawtab_type
 use m_paw_sphharm, only : initylmr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawtwdij_2b'
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments---------------------------
 !scalars
 type(orbmag_type),intent(in) :: dtorbmag
 integer,intent(in) :: ntypat
 type(pawang_type),intent(in) :: pawang
 real(dp),intent(out) :: calc_twdij(2,24,dtorbmag%lmn2max,dtorbmag%natom)
!arrays
 integer,intent(in) :: typat(dtorbmag%natom)
 real(dp),intent(in) :: gprimd(3,3),xred(3,dtorbmag%natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables---------------------------
!scalars
 integer :: adir,bdir,bfor,bsigma,eabg,epsabg,gdir,gfor,gsigma,iatom,ir,isel,itypat
 integer :: klm,kln,klmn,lbess,lbesslm,lmin,lmax,mbess,mesh_size,twindx
 integer :: ylmr_normchoice,ylmr_npts,ylmr_option
 real(dp) :: arg,bessg,bnorm,intg,rterm
 complex(dpc) :: cterm,etb,ifac
!arrays 
 real(dp) :: bb(3),bbn(3),bcart(3),dkb(3),dkg(3),dkvecs(3)
 real(dp) :: calc_expibi(2,dtorbmag%natom),ylmgr(1,1,0),ylmr_nrm(1)
 real(dp),allocatable :: ff(:),j_bessel(:,:),ylmb(:),sb_out(:)
! the following is (i)^L mod 4.
 complex(dpc),dimension(0:3) :: il(0:3)=(/cone,j_dpc,-cone,-j_dpc/)

 ! *************************************************************************

 calc_twdij(:,:,:,:) = zero

 ylmr_normchoice = 0 ! input to initylmr are normalized
 ylmr_npts = 1 ! only 1 point to compute in initylmr
 ylmr_nrm(1) = one ! weight of normed point for initylmr
 ylmr_option = 1 ! compute only ylm's in initylmr

 ABI_ALLOCATE(sb_out, (pawang%l_size_max))

 do adir = 1, 3

    do epsabg = 1, -1, -2

       if (epsabg .EQ. 1) then
          bdir = modulo(adir,3)+1
          gdir = modulo(adir+1,3)+1
          eabg = 1
       else
          bdir = modulo(adir+1,3)+1
          gdir = modulo(adir,3)+1
          eabg = 2
       end if

       do bfor = 1, 2
          if (bfor .EQ. 1) then
             bsigma = 1
          else
             bsigma = -1
          end if

          dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)

          do gfor = 1, 2
             if (gfor .EQ. 1) then
                gsigma = 1
             else
                gsigma = -1
             end if

             dkg(1:3) = gsigma*dtorbmag%dkvecs(1:3,gdir)

             ! <u_kb|H|u_kg>
             ! ket is ahead of bra
             ! based off <u_nk|u_nk+b> in Berrys phase code
             ! here based on B_\alpha<u_beta|H|u_gamma>
             dkvecs(:) = dkg(:) - dkb(:)

             call expibi(calc_expibi,dkvecs,dtorbmag%natom,xred)

             twindx = (adir-1)*8 + (eabg-1)*4 + (bfor-1)*2 + gfor

             do iatom = 1, dtorbmag%natom

                itypat = typat(iatom)
                mesh_size = pawtab(itypat)%mesh_size

                ABI_ALLOCATE(j_bessel,(mesh_size,pawang%l_size_max))
                ABI_ALLOCATE(ff,(mesh_size))
                ABI_ALLOCATE(ylmb,(pawang%l_size_max*pawang%l_size_max))

                !    here is exp(-i (k_ket - k_bra).R) for current atom
                etb = cmplx(calc_expibi(1,iatom),calc_expibi(2,iatom))

                !    note the definition used for the k-dependence of the PAW basis functions:
                !$|\phi_{i,k}\rangle = exp(-i k\cdot r)|\phi_i\rangle
                !    see Umari, Gonze, and Pasquarello, PRB 69,235102 Eq. 23. Thus the k-vector on the
                !    bra side enters as k, while on the ket side it enters as -k.
                bb(:) = -dkvecs(:)

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
                   j_bessel(ir,:) = sb_out
                end do ! end loop over mesh

                !    compute Y_LM(b) here
                call initylmr(pawang%l_size_max,ylmr_normchoice,ylmr_npts,ylmr_nrm,ylmr_option,bbn,ylmb(:),ylmgr)
     
                do klmn = 1, pawtab(itypat)%lmn2_size
                   klm =pawtab(itypat)%indklmn(1,klmn)
                   kln =pawtab(itypat)%indklmn(2,klmn)
                   lmin=pawtab(itypat)%indklmn(3,klmn)
                   lmax=pawtab(itypat)%indklmn(4,klmn)
                   do lbess = lmin, lmax, 2    ! only possible choices for L s.t. Gaunt integrals
                      !        will be non-zero
                      ifac = il(mod(lbess,4))
                      do mbess = -lbess, lbess
                         lbesslm = lbess*lbess+lbess+mbess+1
                         isel=pawang%gntselect(lbesslm,klm)
                         if (isel > 0) then
                            bessg = pawang%realgnt(isel)
                            ff(1:mesh_size)=(pawtab(itypat)%phiphj(1:mesh_size,kln)*&
                                 & pawtab(itypat)%VHnZC(1:mesh_size) - &
                                 & pawtab(itypat)%tphitphj(1:mesh_size,kln)*&
                                 & pawtab(itypat)%vhtnzc(1:mesh_size))*&
                                 & j_bessel(1:mesh_size,lbess+1)
                            call simp_gen(intg,ff,pawrad(itypat))
                            rterm = four_pi*bessg*intg*ylmb(lbesslm)
                            cterm = etb*ifac*rterm
                            calc_twdij(1,twindx,klmn,iatom) = &
                                 &               calc_twdij(1,twindx,klmn,iatom) + dreal(cterm)
                            calc_twdij(2,twindx,klmn,iatom) = &
                                 &               calc_twdij(2,twindx,klmn,iatom) + dimag(cterm)
             
                         end if ! end selection on non-zero Gaunt factors
                      end do ! end loop on mbess = -lbess, lbess
                   end do ! end loop on lmin-lmax bessel l values
                end do ! end loop on lmn2_size klmn basis pairs

                ABI_DEALLOCATE(j_bessel)
                ABI_DEALLOCATE(ff)
                ABI_DEALLOCATE(ylmb)
             end do ! end loop over atoms
          end do ! end loop over gfor
       end do ! end loop over bfor
    end do ! end loop over epsabg
 end do ! end loop over adir
 
 ABI_DEALLOCATE(sb_out)

 end subroutine pawtwdij_2b
!!***
 

end module m_orbmag
!!***
