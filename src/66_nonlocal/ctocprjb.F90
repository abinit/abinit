!{\src2tex{textfont=tt}}
!!****f* ABINIT/ctocprjb
!! NAME
!! ctocprjb
!!
!! FUNCTION
!! Compute <p_k+b|u_k> cprj's as needed by orbital magnetization,
!! at all k points and all bands
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
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

subroutine ctocprjb(atindx1,cg,cprj_kb_k,dtorbmag,dtset,gmet,gprimd,&
     & istwf_k,kg,mcg,mpi_enreg,nattyp,ncpgr,npwarr,pawtab,psps,rmet,rprimd,ucvol,xred)

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_xmpi
 use m_errors
 use m_profiling_abi
 use m_orbmag

 use m_kg, only : getph, ph1d3d
 use m_initylmg, only : initylmg
 use m_mkffnl,           only : mkffnl
 use m_pawtab,           only : pawtab_type
 use m_pawcprj,  only :  pawcprj_alloc, pawcprj_free, pawcprj_getdim, pawcprj_put, pawcprj_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ctocprjb'
 use interfaces_66_nonlocal, except_this_one => ctocprjb
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 !scalars
 integer,intent(in) :: istwf_k,mcg,ncpgr
 real(dp),intent(in) :: ucvol
 type(dataset_type),intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg
 type(orbmag_type), intent(inout) :: dtorbmag
 type(pseudopotential_type),intent(in) :: psps

 !arrays
 integer,intent(in) :: atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem)
 integer,intent(in) :: nattyp(dtset%ntypat),npwarr(dtset%nkpt)
 real(dp),intent(in) :: cg(2,mcg),gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),xred(3,dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
 type(pawcprj_type),intent(inout) :: cprj_kb_k(dtorbmag%fnkpt,3,2,dtset%natom,dtorbmag%nspinor*dtset%mband)

 !Locals------------------------------------
 !scalars
 integer :: bdir,bfor,bsigma,choice,cpopt,dimffnl,dimph1d,ia,iband,icg,ider,idir,ikg,ikpt
 integer :: n1,n2,n3,nband_k,nkpg,npw_k,optder
 real(dp) :: arg

 !arrays
 integer :: nband_dum(1),npwarr_dum(1)
 integer,allocatable :: dimlmn(:),kg_k(:,:)
 real(dp) :: dkb(3),kpoint(3),kpointb(3),kptns(3,1)
 real(dp),allocatable :: cwavef(:,:),ffnl(:,:,:,:),kpg_k_dummy(:,:)
 real(dp),allocatable :: ph1d(:,:),ph3d(:,:,:),phkxred(:,:)
 real(dp),allocatable :: ylm_k(:,:),ylm_k_gr(:,:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

 ! ***********************************************************************

 nband_k = dtorbmag%mband_occ

 ABI_ALLOCATE(dimlmn,(dtset%natom))
 call pawcprj_getdim(dimlmn,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'R')

 ABI_DATATYPE_ALLOCATE(cwaveprj,(dtset%natom,1))
 call pawcprj_alloc(cwaveprj,ncpgr,dimlmn)

 ABI_ALLOCATE(phkxred,(2,dtset%natom))

 n1=dtset%ngfft(1); n2=dtset%ngfft(2); n3=dtset%ngfft(3)
 dimph1d=dtset%natom*(2*(n1+n2+n3)+3)
 ABI_ALLOCATE(ph1d,(2,dimph1d))
 call getph(atindx1,dtset%natom,n1,n2,n3,ph1d,xred)

 do ikpt=1,dtorbmag%fnkpt

    kpoint(:)=dtorbmag%fkptns(:,ikpt)
    npw_k = npwarr(ikpt)
    ABI_ALLOCATE(cwavef,(2,npw_k))

    dimffnl=1
    ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,dtset%ntypat))

    icg = dtorbmag%cgindex(ikpt,dtset%nsppol)

    ikg = dtorbmag%fkgindex(ikpt)
    ABI_ALLOCATE(kg_k,(3,npw_k))
    kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
    nkpg = 0
    ABI_ALLOCATE(kpg_k_dummy,(npw_k,nkpg))

    ABI_ALLOCATE(ph3d,(2,npw_k,dtset%natom))

    ! data for initylmg call below
    optder=0
    nband_dum(1) = nband_k
    npwarr_dum(1) = npw_k
    ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang))
    ABI_ALLOCATE(ylm_k_gr,(npw_k,3+6*(optder/2),psps%mpsang*psps%mpsang))

    do bdir=1, 3
       do bfor=1, 2

          if (bfor .EQ. 1) then
             bsigma = 1
          else
             bsigma = -1
          end if

          dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)

          kpointb(:) = kpoint(:) + dkb(:)

          do ia=1, dtset%natom
             arg=two_pi*(kpointb(1)*xred(1,ia)+kpointb(2)*xred(2,ia)+kpointb(3)*xred(3,ia))
             phkxred(1,ia)=cos(arg);phkxred(2,ia)=sin(arg)
          end do

          call ph1d3d(1,dtset%natom,kg_k,dtset%natom,dtset%natom,npw_k,n1,n2,n3,phkxred,ph1d,ph3d)

          kptns(:,1) = kpointb(:)
          call initylmg(gprimd,kg_k,kptns,1,mpi_enreg,psps%mpsang,npw_k,&
               & nband_dum,1,npwarr_dum,dtset%nsppol,optder,rprimd,ylm_k,ylm_k_gr)

          !      Compute nonlocal form factors ffnl at all (k+b+G_k):
          ider=0;idir=0
          call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
               & gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k_dummy,kpointb,psps%lmnmax,&
               & psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
               & npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
               & psps%usepaw,psps%useylm,ylm_k,ylm_k_gr)

          choice=1;cpopt=0;idir=0
          do iband = 1, nband_k
             cwavef(1,1:npw_k) = cg(1,icg+(iband-1)*npw_k+1:icg+iband*npw_k)
             cwavef(2,1:npw_k) = cg(2,icg+(iband-1)*npw_k+1:icg+iband*npw_k)

             call getcprj(choice,cpopt,cwavef,cwaveprj,ffnl,&
                  & idir,psps%indlmn,istwf_k,kg_k,kpg_k_dummy,kpointb,psps%lmnmax,&
                  & dtset%mgfft,mpi_enreg,&
                  & dtset%natom,nattyp,dtset%ngfft,dtset%nloalg,npw_k,dtset%nspinor,dtset%ntypat,&
                  & phkxred,ph1d,ph3d,ucvol,psps%useylm)
             
             call pawcprj_put(atindx1,cwaveprj,cprj_kb_k(ikpt,bdir,bfor,:,:),dtset%natom,&
                  & iband,0,ikpt,0,1,nband_k,1,dtset%natom,1,nband_k,dimlmn,dtset%nspinor,dtset%nsppol,0)

          end do ! end loop over bands

       end do ! end loop over bfor
    end do ! end loop over bdir

    ABI_DEALLOCATE(kg_k)
    ABI_DEALLOCATE(ph3d)
    ABI_DEALLOCATE(kpg_k_dummy)
    ABI_DEALLOCATE(cwavef)
    ABI_DEALLOCATE(ffnl)
    ABI_DEALLOCATE(ylm_k)
    ABI_DEALLOCATE(ylm_k_gr)
    
 end do ! end loop over nkpt

 ABI_DEALLOCATE(dimlmn)
 call pawcprj_free(cwaveprj)
 ABI_DATATYPE_DEALLOCATE(cwaveprj)

 ABI_DEALLOCATE(phkxred)
 ABI_DEALLOCATE(ph1d)
  
end subroutine ctocprjb
!!***
