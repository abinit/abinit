!{\src2tex{textfont=tt}}
!!****f* ABINIT/chern_number
!! NAME
!! chern_number
!!
!! FUNCTION
!! This routine computes the Chern number based on input wavefunctions.
!! It is assumed that only completely filled bands are present.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT  group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!! cg(2,mcg)=planewave coefficients of wavefunctions
!! cprj(natom,mcprj*usecrpj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!! dtset <type(dataset_type)>=all input variables in this dataset
!! mband=maximum number of bands
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!! mkmem=number of k points treated by this node
!! mpi_enreg=information about MPI parallelization
!! mpw=maximum dimensioned size of npw
!! my_natom=number of atoms treated by current processor
!! natom=number of atoms in cell
!! nkpt=number of k points
!! npwarr(nkpt)=number of planewaves in basis at this k point
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! ntypat=number of types of atoms in unit cell
!! nkpt=number of k-points
!! pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!! typat(natom)=type integer for each atom in cell
!! usecprj=1 if cprj datastructure has been allocated
!! usepaw=1 if PAW calculation
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! dtorbmag <type(orbmag_type)> = variables related to orbital magnetization
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

subroutine chern_number(atindx1,cg,cprj,dtset,dtorbmag,&
     &            mcg,mcprj,mpi_enreg,npwarr,ntypat,pawtab,usecprj,usepaw)

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_xmpi
 use m_errors
 use m_orbmag
 use m_profiling_abi

 use m_pawtab,   only : pawtab_type
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_get, pawcprj_getdim

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chern_number'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 !scalars
 integer,intent(in) :: mcg,mcprj,ntypat,usecprj,usepaw
 type(dataset_type),intent(in) :: dtset
 type(MPI_type), intent(in) :: mpi_enreg
 type(orbmag_type), intent(inout) :: dtorbmag

 !arrays
 integer,intent(in) :: atindx1(dtset%natom),npwarr(dtset%nkpt)
 real(dp), intent(in) :: cg(2,mcg)
 type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj*usecprj)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

 !Local variables -------------------------
 !scalars
 integer :: bband,bbs,bra_end,bra_start,iatom,ikpt,ilmn,isppol,itypat
 integer :: kband,kbs,ket_end,ket_start,klmn,jlmn,ipw,ispinor
 integer :: my_nspinor,nband_k,ncpgr,npw_k,spnshft,spnipw
 real(dp) :: err_ovlp,mag_ovlp,max_err_ovlp,ovlp_i,ovlp_r,paw_i,paw_r,tot_r,tot_i
 complex(dpc) :: cpb,cpk,cterm

 !arrays
 integer,allocatable :: dimlmn(:),nattyp_dum(:)
 real(dp),allocatable :: bra(:,:),ket(:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:)

 ! ***********************************************************************
 my_nspinor=max(1,dtorbmag%nspinor/mpi_enreg%nproc_spinor)

 if (usepaw == 1) then ! cprj allocation
    ncpgr = cprj(1,1)%ncpgr
    ABI_ALLOCATE(dimlmn,(dtset%natom))
    call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,ntypat,dtset%typat,pawtab,'R')
    ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,dtorbmag%nspinor*dtset%mband))
    call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
 end if

 ! =======================================
 ! code to test orthonormality of cg_k
 ! =======================================
 
 ikpt = 3
 npw_k = npwarr(ikpt)
 isppol = 1
 nband_k = dtorbmag%mband_occ
 ABI_ALLOCATE(bra,(2,npw_k*my_nspinor))
 ABI_ALLOCATE(ket,(2,npw_k*my_nspinor))
 max_err_ovlp=0.0

 call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,dtorbmag%cprjindex(ikpt,isppol),ikpt,0,isppol,dtset%mband,&
&                 dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)
 do bband = 1, nband_k
    bra_start = dtorbmag%cgindex(ikpt,dtset%nsppol)+1+(bband-1)*npw_k*my_nspinor
    bra_end = bra_start + npw_k*my_nspinor - 1
    bra(1:2,1:npw_k*my_nspinor) = cg(1:2,bra_start:bra_end)
    do kband = 1, nband_k
       ket_start = dtorbmag%cgindex(ikpt,dtset%nsppol)+1+(kband-1)*npw_k*my_nspinor
       ket_end = ket_start + npw_k*my_nspinor - 1
       ket(1:2,1:npw_k*my_nspinor) = cg(1:2,ket_start:ket_end)

       tot_r = 0.0; tot_i = 0.0     
       do ispinor = 1, my_nspinor
          ovlp_r = 0.0; ovlp_i = 0.0
          spnshft = (ispinor-1)*npw_k
          do ipw = 1, npw_k
             spnipw = ipw + spnshft
             ovlp_r = ovlp_r + bra(1,spnipw)*ket(1,spnipw)+bra(2,spnipw)*ket(2,spnipw)
             ovlp_i = ovlp_i - bra(2,spnipw)*ket(1,spnipw)+bra(1,spnipw)*ket(2,spnipw)
          end do ! end loop over ipw
          paw_r = 0.0; paw_i = 0.0
          do iatom = 1, dtset%natom
             itypat = dtset%typat(iatom)
             do ilmn = 1, dtorbmag%lmn_size(itypat)
                do jlmn = 1, dtorbmag%lmn_size(itypat)
                   klmn=max(ilmn,jlmn)*(max(ilmn,jlmn)-1)/2 + min(ilmn,jlmn)
                   bbs = my_nspinor*(bband-1)+ispinor
                   kbs = my_nspinor*(kband-1)+ispinor
                   cpb=cmplx(cprj_k(iatom,bbs)%cp(1,ilmn),cprj_k(iatom,bbs)%cp(2,ilmn))
                   cpk=cmplx(cprj_k(iatom,kbs)%cp(1,jlmn),cprj_k(iatom,kbs)%cp(2,jlmn))
                   cterm = conjg(cpb)*pawtab(itypat)%sij(klmn)*cpk
                   paw_r = paw_r + real(cterm)
                   paw_i = paw_i + aimag(cterm)
                end do ! end loop over jlmn
             end do ! end loop over ilmn
          end do ! end loop over iatom
          tot_r = tot_r + ovlp_r + paw_r
          tot_i = tot_i + ovlp_i + paw_i
       end do ! end loop over ispinor

       write(std_out,'(a,2i4,2es16.8)')' JWZ Debug: chern_number bband kband ovlp : ',&
&                                        bband,kband,tot_r,tot_i
       mag_ovlp =  tot_r*tot_r + tot_i*tot_i
       if(bband==kband) then
          err_ovlp=abs(mag_ovlp-1.0)
       else
          err_ovlp=abs(mag_ovlp)
       end if
       max_err_ovlp=MAX(max_err_ovlp,err_ovlp)
    end do ! end loop over kband
 end do ! end loop over bband
 write(std_out,'(a,i4,es16.8)')' JWZ Debug: chern_number ikpt ovlp err : ',&
&                                ikpt,max_err_ovlp
 ABI_DEALLOCATE(bra)
 ABI_DEALLOCATE(ket)

! =========================================
! end code to test orthonormality of cg_k
! =========================================

 if (usepaw == 1) then
    ABI_DEALLOCATE(dimlmn)
    call pawcprj_free(cprj_k)
    ABI_DATATYPE_DEALLOCATE(cprj_k)
 end if

end subroutine chern_number
!!***
