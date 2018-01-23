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
!! gmet(3,3)=metric in reciprocal space
!! gprimd(3,3)=reciprocal space dimensional primitive translations
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!! mpi_enreg=information about MPI parallelization
!! npwarr(nkpt)=number of planewaves in basis at this k point
!! pawang <type(pawang_type)>=paw angular mesh and related data
!! pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!! pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!! pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!! pwind_alloc = first dimension of pwind
!! usecprj=1 if cprj datastructure has been allocated
!! usepaw=1 if PAW calculation
!! xred(3,natom) = location of atoms in unit cell
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

subroutine chern_number(atindx1,cg,cprj,dtset,dtorbmag,gmet,gprimd,&
     &            mcg,mcprj,mpi_enreg,npwarr,pawang,pawrad,pawtab,pwind,pwind_alloc,&
     &            usecprj,usepaw,xred)

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_xmpi
 use m_errors
 use m_orbmag
 use m_profiling_abi

 use m_pawang,           only : pawang_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_get, pawcprj_getdim

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chern_number'
 use interfaces_32_util
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 !scalars
 integer,intent(in) :: mcg,mcprj,pwind_alloc,usecprj,usepaw
 type(dataset_type),intent(in) :: dtset
 type(MPI_type), intent(in) :: mpi_enreg
 type(orbmag_type), intent(inout) :: dtorbmag
 type(pawang_type),intent(in) :: pawang

 !arrays
 integer,intent(in) :: atindx1(dtset%natom),npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
 real(dp), intent(in) :: cg(2,mcg),gmet(3,3),gprimd(3,3),xred(3,dtset%natom)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*usepaw)
 type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj*usecprj)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)

 !Local variables -------------------------
 !scalars
 integer :: adir,bdir,bpw,bsigma,epsabg,gdir,gpw,gsigma,icg,icgb,icgg
 integer :: ikg,ikpt,ikptb,ikptg,ipw,isppol,my_nspinor
 integer :: nn,nnp,nnpp,nband_k,ncpgr,npw_ka,npw_kb,npw_kg
 real(dp) :: bdelta,doti,dotr,gdelta,tmpi,tmpr
 complex(dpc) :: c_nn_nnp,c_nnp_nnpp,c_nnpp_nn,cpaw

 !arrays
 integer,allocatable :: dimlmn(:),nattyp_dum(:),pwind_kb(:),pwind_kbg(:),pwind_kg(:)
 real(dp) :: dkb(3),dkg(3),dkbg(3),mdkg(3)
 real(dp),allocatable :: kkb_paw(:,:,:),kbkg_paw(:,:,:),kgk_paw(:,:,:)
 real(dp),allocatable :: vecta(:,:),vectb(:,:),vectg(:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kb(:,:),cprj_kg(:,:)

 ! ***********************************************************************
 ! my_nspinor=max(1,dtorbmag%nspinor/mpi_enreg%nproc_spinor)

 ! TODO: generalize to nsppol > 1
 isppol = 1
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 
 nband_k = dtorbmag%mband_occ

 if (usepaw == 1) then ! cprj allocation
    ncpgr = cprj(1,1)%ncpgr
    ABI_ALLOCATE(dimlmn,(dtset%natom))
    call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,dtset%ntypat,dtset%typat,pawtab,'R')
    ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,dtorbmag%nspinor*dtset%mband))
    call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
    ABI_DATATYPE_ALLOCATE(cprj_kb,(dtset%natom,dtorbmag%nspinor*dtset%mband))
    call pawcprj_alloc(cprj_kb,ncpgr,dimlmn)
    ABI_DATATYPE_ALLOCATE(cprj_kg,(dtset%natom,dtorbmag%nspinor*dtset%mband))
    call pawcprj_alloc(cprj_kg,ncpgr,dimlmn)
 end if

 ABI_ALLOCATE(kkb_paw,(2,dtset%mband,dtset%mband))
 ABI_ALLOCATE(kbkg_paw,(2,dtset%mband,dtset%mband))
 ABI_ALLOCATE(kgk_paw,(2,dtset%mband,dtset%mband))
 ABI_ALLOCATE(vecta,(2,0:dtset%mpw*dtorbmag%nspinor))
 ABI_ALLOCATE(vectb,(2,0:dtset%mpw*dtorbmag%nspinor))
 ABI_ALLOCATE(vectg,(2,0:dtset%mpw*dtorbmag%nspinor))
 ABI_ALLOCATE(pwind_kb,(dtset%mpw))
 ABI_ALLOCATE(pwind_kbg,(dtset%mpw))
 ABI_ALLOCATE(pwind_kg,(dtset%mpw))
 vecta(:,0) = zero; vectb(:,0) = zero; vectg(:,0) = zero


 do adir = 1, 3

    dtorbmag%chern(:,adir) = zero

    do epsabg = 1, -1, -2

       if (epsabg .EQ. 1) then
          bdir = modulo(adir,3)+1
          gdir = modulo(adir+1,3)+1
       else
          bdir = modulo(adir+1,3)+1
          gdir = modulo(adir,3)+1
       end if

       ! loop over kpts, assuming for now kptopt 3, nsppol = 1, nspinor = 1, no parallelism, no symmorphic symmetry elements
       do ikpt = 1, dtorbmag%fnkpt
          icg = dtorbmag%cgindex(ikpt,isppol)
          ikg = dtorbmag%fkgindex(ikpt)
          npw_ka = npwarr(ikpt)
          call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,dtorbmag%cprjindex(ikpt,isppol),&
               & ikpt,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
               & my_nspinor,dtset%nsppol,0)
       
          do bsigma = -1, 1, 2
             ikptb = dtorbmag%ikpt_dk(ikpt,(bsigma+3)/2,bdir)
             icgb = dtorbmag%cgindex(ikptb,isppol)
             npw_kb = npwarr(ikptb)
             pwind_kb(1:npw_ka) = pwind(ikg+1:ikg+npw_ka,(bsigma+3)/2,bdir)
             dkb = bsigma*dtorbmag%dkvecs(1:3,bdir)
             bdelta = sqrt(DOT_PRODUCT(dkb,MATMUL(gmet,dkb)))
             call pawcprj_get(atindx1,cprj_kb,cprj,dtset%natom,1,dtorbmag%cprjindex(ikptb,isppol),&
                  & ikptb,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                  & my_nspinor,dtset%nsppol,0)

             call overlap_k1k2_paw(cprj_k,cprj_kb,dkb,gprimd,kkb_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
                  & dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

             do gsigma = -1, 1, 2
                ikptg = dtorbmag%ikpt_dk(ikpt,(gsigma+3)/2,gdir)
                icgg = dtorbmag%cgindex(ikptg,isppol)
                npw_kg = npwarr(ikptg)
                pwind_kg(1:npw_ka) = pwind(ikg+1:ikg+npw_ka,(gsigma+3)/2,gdir)
                dkg = gsigma*dtorbmag%dkvecs(1:3,gdir)
                gdelta = sqrt(DOT_PRODUCT(dkg,MATMUL(gmet,dkg)))
                call pawcprj_get(atindx1,cprj_kg,cprj,dtset%natom,1,dtorbmag%cprjindex(ikptg,isppol),&
                     & ikptg,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                     & my_nspinor,dtset%nsppol,0)
                
                dkbg = dkg - dkb
                call overlap_k1k2_paw(cprj_kb,cprj_kg,dkbg,gprimd,kbkg_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
                     & dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                mdkg = -dkg
                call overlap_k1k2_paw(cprj_kg,cprj_k,mdkg,gprimd,kgk_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
                     & dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                do nn = 1, dtorbmag%mband_occ
                   vecta(1:2,1:npw_ka) = cg(1:2,icg+(nn-1)*npw_ka+1:icg+nn*npw_ka)
                   do nnp = 1, dtorbmag%mband_occ
                      vectb(1:2,1:npw_kb) = cg(1:2,icgb+(nnp-1)*npw_kb+1:icg+nnp*npw_kb)

                      call overlap_g(doti,dotr,dtset%mpw,npw_ka,npw_kb,1,pwind_kb,vecta,vectb)
!                      dotr=zero;doti=zero
                      c_nn_nnp = cmplx(dotr+kkb_paw(1,nn,nnp),doti+kkb_paw(2,nn,nnp))
                   
                      do nnpp = 1, dtorbmag%mband_occ
                         vectg(1:2,1:npw_kg) = cg(1:2,icgg+(nnpp-1)*npw_kg+1:icg+nnpp*npw_kg)

!                         call overlap_g(doti,dotr,dtset%mpw,npw_kb,npw_kg,1,pwind_kb,vectb,vectg)
                         dotr=zero;doti=zero
                         do ipw = 1, npw_ka
                            bpw = pwind_kb(ipw)
                            gpw = pwind_kg(ipw)
                            dotr = dotr + vectb(1,bpw)*vectg(1,gpw)-vectb(2,bpw)*vectg(2,gpw)
                            doti = doti + vectb(1,bpw)*vectg(2,gpw)+vectb(2,bpw)*vectg(1,gpw)
                         end do
                         c_nnp_nnpp = cmplx(kbkg_paw(1,nnp,nnpp)+dotr,kbkg_paw(2,nnp,nnpp)+doti)

                         call overlap_g(doti,dotr,dtset%mpw,npw_ka,npw_kg,1,pwind_kg,vecta,vectg)
!                         dotr=zero;doti=zero
                         c_nnpp_nn = cmplx(kgk_paw(1,nnpp,nn)+dotr,kgk_paw(2,nnpp,nn)-doti)

                         cpaw = epsabg * bsigma * gsigma * c_nn_nnp * c_nnp_nnpp * c_nnpp_nn/(2.0*bdelta*2.0*gdelta)
                         dtorbmag%chern(1,adir) = dtorbmag%chern(1,adir) + real(cpaw)
                         dtorbmag%chern(2,adir) = dtorbmag%chern(2,adir) + aimag(cpaw)
                      end do ! end loop over nnpp bands
                   end do ! end loop over nnp bands
                end do ! end loop over nn bands
             end do ! end loop over gsigma
          end do ! end loop over bsigma
       end do ! end loop over ikpt
    end do ! end loop over epsabg
 end do ! end loop over adir

 do adir = 1, 3
    tmpr = -dtorbmag%chern(2,adir)/(two_pi*dtorbmag%fnkpt)
    tmpi = dtorbmag%chern(1,adir)/(two_pi*dtorbmag%fnkpt)
    dtorbmag%chern(1,adir) = tmpr
    dtorbmag%chern(2,adir) = tmpi
    write(std_out,'(a,i4,2es16.8)')' JWZ Debug: idir chern_number : ',&
         &   adir,dtorbmag%chern(1,adir),dtorbmag%chern(2,adir)
 end do
 
 
!  ! =======================================
!  ! code to test orthonormality of cg_k
!  ! =======================================
 
!  ikpt = 2
!  npw_k = npwarr(ikpt)
!  isppol = 1
!  nband_k = dtorbmag%mband_occ
!  ABI_ALLOCATE(bra,(2,npw_k*my_nspinor))
!  ABI_ALLOCATE(ket,(2,npw_k*my_nspinor))
!  max_err_ovlp=0.0

!  call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,dtorbmag%cprjindex(ikpt,isppol),ikpt,0,isppol,dtset%mband,&
! &                 dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)
!  do bband = 1, nband_k
!     bra_start = dtorbmag%cgindex(ikpt,dtset%nsppol)+1+(bband-1)*npw_k*my_nspinor
!     bra_end = bra_start + npw_k*my_nspinor - 1
!     bra(1:2,1:npw_k*my_nspinor) = cg(1:2,bra_start:bra_end)
!     do kband = 1, nband_k
!        ket_start = dtorbmag%cgindex(ikpt,dtset%nsppol)+1+(kband-1)*npw_k*my_nspinor
!        ket_end = ket_start + npw_k*my_nspinor - 1
!        ket(1:2,1:npw_k*my_nspinor) = cg(1:2,ket_start:ket_end)

!        tot_r = 0.0; tot_i = 0.0     
!        do ispinor = 1, my_nspinor
!           ovlp_r = 0.0; ovlp_i = 0.0
!           spnshft = (ispinor-1)*npw_k
!           do ipw = 1, npw_k
!              spnipw = ipw + spnshft
!              ovlp_r = ovlp_r + bra(1,spnipw)*ket(1,spnipw)+bra(2,spnipw)*ket(2,spnipw)
!              ovlp_i = ovlp_i - bra(2,spnipw)*ket(1,spnipw)+bra(1,spnipw)*ket(2,spnipw)
!           end do ! end loop over ipw
!           paw_r = 0.0; paw_i = 0.0
!           do iatom = 1, dtset%natom
!              itypat = dtset%typat(iatom)
!              do ilmn = 1, dtorbmag%lmn_size(itypat)
!                 do jlmn = 1, dtorbmag%lmn_size(itypat)
!                    klmn=max(ilmn,jlmn)*(max(ilmn,jlmn)-1)/2 + min(ilmn,jlmn)
!                    bbs = my_nspinor*(bband-1)+ispinor
!                    kbs = my_nspinor*(kband-1)+ispinor
!                    cpb=cmplx(cprj_k(iatom,bbs)%cp(1,ilmn),cprj_k(iatom,bbs)%cp(2,ilmn))
!                    cpk=cmplx(cprj_k(iatom,kbs)%cp(1,jlmn),cprj_k(iatom,kbs)%cp(2,jlmn))
!                    cterm = conjg(cpb)*pawtab(itypat)%sij(klmn)*cpk
!                    paw_r = paw_r + real(cterm)
!                    paw_i = paw_i + aimag(cterm)
!                 end do ! end loop over jlmn
!              end do ! end loop over ilmn
!           end do ! end loop over iatom
!           tot_r = tot_r + ovlp_r + paw_r
!           tot_i = tot_i + ovlp_i + paw_i
!        end do ! end loop over ispinor

!        write(std_out,'(a,2i4,2es16.8)')' JWZ Debug: chern_number bband kband ovlp : ',&
! &                                        bband,kband,tot_r,tot_i
!        mag_ovlp =  tot_r*tot_r + tot_i*tot_i
!        if(bband==kband) then
!           err_ovlp=abs(mag_ovlp-1.0)
!        else
!           err_ovlp=abs(mag_ovlp)
!        end if
!        max_err_ovlp=MAX(max_err_ovlp,err_ovlp)
!     end do ! end loop over kband
!  end do ! end loop over bband
!  write(std_out,'(a,i4,es16.8)')' JWZ Debug: chern_number ikpt ovlp err : ',&
! &                                ikpt,max_err_ovlp
!  ABI_DEALLOCATE(bra)
!  ABI_DEALLOCATE(ket)

! ! =========================================
! ! end code to test orthonormality of cg_k
! ! =========================================

 if (usepaw == 1) then
    ABI_DEALLOCATE(dimlmn)
    call pawcprj_free(cprj_k)
    ABI_DATATYPE_DEALLOCATE(cprj_k)
    call pawcprj_free(cprj_kb)
    ABI_DATATYPE_DEALLOCATE(cprj_kb)
    call pawcprj_free(cprj_kg)
    ABI_DATATYPE_DEALLOCATE(cprj_kg)
 end if

 ABI_DEALLOCATE(kkb_paw)
 ABI_DEALLOCATE(kbkg_paw)
 ABI_DEALLOCATE(kgk_paw)
 ABI_DEALLOCATE(vecta)
 ABI_DEALLOCATE(vectb)
 ABI_DEALLOCATE(vectg)
 ABI_DEALLOCATE(pwind_kb)
 ABI_DEALLOCATE(pwind_kbg)
 ABI_DEALLOCATE(pwind_kg)

end subroutine chern_number
!!***
