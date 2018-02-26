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
!! Copyright (C) 2003-2017 ABINIT  group (JZwanziger, MVeithen)
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
!! kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
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
!! symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!   reciprocal space primitive translations
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
!! See Ceresoli et al, PRB 74, 024408 (2006), and Gonze and Zwanziger, PRB 84
!! 064446 (2011). This routine computes the Chern number as
!! $C_\alpha = \frac{i}{2\pi}\int_{\mathrm{BZ}} dk \epsilon_{\alpha\beta\gamma}
!! \mathrm{Tr}[\rho_k \partial_\beta \rho_k (1 - \rho_k) \partial_gamma\rho_k] $
!! The derivative of the density operator is obtained from a discretized formula
!! $\partial_\beta \rho_k = \frac{1}{2\Delta}(\rho_{k+b} - \rho_{k-b})$ with
!! $\Delta = |b|$. When reduced to wavefunction overlaps the computation amounts to
!! multiple calls to smatrix.F90, exactly as in other Berry phase computations, with
!! the one additional complication of overlaps like $\langle u_{n,k+b}|u_{n',k+g}\rangle$.
!! At this stage mkpwind_k is invoked, which generalizes the code in initberry
!! and initorbmag necessary to index plane waves around different k points.
!! Direct questions and comments to J Zwanziger
!!
!! PARENTS
!!      afterscfloop
!!
!! CHILDREN
!!      mkpwind_k,overlap_k1k2_paw,pawcprj_alloc,pawcprj_free,pawcprj_get
!!      pawcprj_getdim,smatrix,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine chern_number(atindx1,cg,cprj,dtset,dtorbmag,gmet,gprimd,kg,&
&            mcg,mcprj,mpi_enreg,npwarr,pawang,pawrad,pawtab,pwind,pwind_alloc,&
&            symrec,usecprj,usepaw,xred)

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_xmpi
 use m_errors
 use m_profiling_abi

 use m_orbmag

 use m_fftcore, only : kpgsph
 use m_pawang,           only : pawang_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_copy, pawcprj_free,&
      & pawcprj_get, pawcprj_getdim, pawcprj_set_zero

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chern_number'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_56_recipspace
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcg,mcprj,pwind_alloc,usecprj,usepaw
 type(dataset_type),intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg
 type(orbmag_type), intent(inout) :: dtorbmag
 type(pawang_type),intent(in) :: pawang

!arrays
 integer,intent(in) :: atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem)
 integer,intent(in) :: npwarr(dtset%nkpt),pwind(pwind_alloc,2,3),symrec(3,3,dtset%nsym)
 real(dp), intent(in) :: cg(2,mcg),gmet(3,3),gprimd(3,3),xred(3,dtset%natom)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*usepaw)
 type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj*usecprj)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)

!Local variables -------------------------
!scalars
 integer :: adir,bdir,bfor,bsigma,ddkflag,epsabg,gdir,gfor,gsigma
 integer :: icg,icgb,icgg,icprj,icprjb,icprjg
 integer :: ikg,ikpt,ikptb,ikptg,isppol,itrs,job
 integer :: mcg1_k,my_nspinor,nband_k,ncpgr,nn,n1,n2,n3,npw_k,npw_kb,npw_kg,shiftbd
 real(dp) :: deltab,deltag
 complex(dpc) :: IA,IB,t1A,t2A,t3A,t1B,t2B,t3B,t4B,tprodA,tprodB
 character(len=500) :: message
 !arrays
 integer,allocatable :: dimlmn(:),nattyp_dum(:),pwind_kb(:),pwind_kg(:),pwind_bg(:),sflag_k(:)
 real(dp) :: cnum(2,3),dkb(3),dkg(3),dkbg(3),dtm_k(2)
 real(dp),allocatable :: cg1_k(:,:),kk_paw(:,:,:),pwnsfac_k(:,:),smat_all(:,:,:,:,:),smat_inv(:,:,:)
 real(dp),allocatable :: smat_kk(:,:,:)
 logical,allocatable :: has_smat(:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kb(:,:),cprj_kg(:,:)
 type(pawcprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:)

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
   if (dtset%kptopt /= 3) then
     ABI_DATATYPE_ALLOCATE(cprj_ikn,(dtset%natom,dtorbmag%nspinor*dtset%mband))
     ABI_DATATYPE_ALLOCATE(cprj_fkn,(dtset%natom,dtorbmag%nspinor*dtset%mband))
     call pawcprj_alloc(cprj_ikn,ncpgr,dimlmn)
     call pawcprj_alloc(cprj_fkn,ncpgr,dimlmn)
   end if
 else
   message = ' usepaw /= 1 but Chern number calculation requires PAW '
   MSG_ERROR(message)
 end if

 ABI_ALLOCATE(kk_paw,(2,dtset%mband,dtset%mband))
 ABI_ALLOCATE(pwind_kb,(dtset%mpw))
 ABI_ALLOCATE(pwind_kg,(dtset%mpw))
 ABI_ALLOCATE(pwind_bg,(dtset%mpw))
 ABI_ALLOCATE(pwnsfac_k,(4,dtset%mpw))
 pwnsfac_k(1,:) = 1.0_dp ! bra real
 pwnsfac_k(2,:) = 0.0_dp ! bra imag
 pwnsfac_k(3,:) = 1.0_dp ! ket real
 pwnsfac_k(4,:) = 0.0_dp ! ket imag

 mcg1_k = dtset%mpw*nband_k
 ABI_ALLOCATE(cg1_k,(2,mcg1_k))
 ABI_ALLOCATE(sflag_k,(nband_k))
 ABI_ALLOCATE(smat_inv,(2,nband_k,nband_k))
 ABI_ALLOCATE(smat_kk,(2,nband_k,nband_k))

 ddkflag = 1
 
!itrs = 0 means do not invoke time reversal symmetry in smatrix.F90
 itrs = 0
 
 job = 1
 shiftbd = 1

 ABI_ALLOCATE(has_smat,(dtorbmag%fnkpt,dtorbmag%fnkpt))
 ABI_ALLOCATE(smat_all,(2,nband_k,nband_k,dtorbmag%fnkpt,dtorbmag%fnkpt))
 has_smat(:,:)=.FALSE.
 
 do adir = 1, 3

   IA = czero
   IB = czero

   do epsabg = 1, -1, -2

     if (epsabg .EQ. 1) then
       bdir = modulo(adir,3)+1
       gdir = modulo(adir+1,3)+1
     else
       bdir = modulo(adir+1,3)+1
       gdir = modulo(adir,3)+1
     end if

       ! loop over kpts, assuming for now kptopt 3, nsppol = 1, nspinor = 1
       ! and no parallelism, no symmorphic symmetry elements
     do ikpt = 1, dtorbmag%fnkpt

       icprj = dtorbmag%cprjindex(ikpt,isppol)
       
       npw_k = npwarr(ikpt)
       icg = dtorbmag%cgindex(ikpt,dtset%nsppol)

       ikg = dtorbmag%fkgindex(ikpt)

       call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprj,ikpt,0,isppol,dtset%mband,&
&       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)

       do bfor = 1, 2
         if (bfor .EQ. 1) then
           bsigma = 1
         else
           bsigma = -1
         end if

         dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)
         deltab = sqrt(DOT_PRODUCT(dkb,MATMUL(gmet,dkb)))

         ikptb = dtorbmag%ikpt_dk(ikpt,bfor,bdir)
         icprjb = dtorbmag%cprjindex(ikptb,isppol)
         
         npw_kb = npwarr(ikptb)
         icgb = dtorbmag%cgindex(ikptb,dtset%nsppol)

         pwind_kb(1:npw_k) = pwind(ikg+1:ikg+npw_k,bfor,bdir)

         call pawcprj_get(atindx1,cprj_kb,cprj,dtset%natom,1,icprjb,&
&         ikptb,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
&         my_nspinor,dtset%nsppol,0)

         if (.NOT. has_smat(ikpt,ikptb)) then
           
           call overlap_k1k2_paw(cprj_k,cprj_kb,dkb,gprimd,kk_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
&           dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

           sflag_k=0
           call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icg,icgb,itrs,job,nband_k,&
&           mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kb,my_nspinor,&
&           pwind_kb,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,usepaw)

           do nn = 1, nband_k
             do n1 = 1, nband_k
               smat_all(1,nn,n1,ikpt,ikptb) =  smat_kk(1,nn,n1)
               smat_all(2,nn,n1,ikpt,ikptb) =  smat_kk(2,nn,n1)
               smat_all(1,n1,nn,ikptb,ikpt) =  smat_kk(1,nn,n1)
               smat_all(2,n1,nn,ikptb,ikpt) = -smat_kk(2,nn,n1)
             end do
           end do
           
           has_smat(ikpt,ikptb) = .TRUE.
           has_smat(ikptb,ikpt) = .TRUE.

         end if
         
         do gfor = 1, 2
           if (gfor .EQ. 1) then
             gsigma = 1
           else
             gsigma = -1
           end if

           dkg(1:3) = gsigma*dtorbmag%dkvecs(1:3,gdir)
           deltag = sqrt(DOT_PRODUCT(dkg,MATMUL(gmet,dkg)))

           ikptg = dtorbmag%ikpt_dk(ikpt,gfor,gdir)
           icprjg = dtorbmag%cprjindex(ikptg,isppol)
           
           npw_kg = npwarr(ikptg)
           icgg = dtorbmag%cgindex(ikptg,dtset%nsppol)

           pwind_kg(1:npw_k) = pwind(ikg+1:ikg+npw_k,gfor,gdir)

           call pawcprj_get(atindx1,cprj_kg,cprj,dtset%natom,1,icprjg,&
&           ikptg,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
&           my_nspinor,dtset%nsppol,0)

           if (.NOT. has_smat(ikpt,ikptg)) then
             
             call overlap_k1k2_paw(cprj_k,cprj_kg,dkg,gprimd,kk_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
&             dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

             sflag_k=0
             call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icg,icgg,itrs,job,nband_k,&
&             mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kg,my_nspinor,&
&             pwind_kg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,usepaw)

             do nn = 1, nband_k
               do n1 = 1, nband_k
                 smat_all(1,nn,n1,ikpt,ikptg) =  smat_kk(1,nn,n1)
                 smat_all(2,nn,n1,ikpt,ikptg) =  smat_kk(2,nn,n1)
                 smat_all(1,n1,nn,ikptg,ikpt) =  smat_kk(1,nn,n1)
                 smat_all(2,n1,nn,ikptg,ikpt) = -smat_kk(2,nn,n1)
               end do
             end do

             has_smat(ikpt,ikptg) = .TRUE.
             has_smat(ikptg,ikpt) = .TRUE.

           end if

           dkbg = dkg - dkb

           if (.NOT. has_smat(ikptb,ikptg)) then
             
             call overlap_k1k2_paw(cprj_kb,cprj_kg,dkbg,gprimd,kk_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
&             dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

             call mkpwind_k(dkbg,dtset,dtorbmag%fnkpt,dtorbmag%fkptns,gmet,dtorbmag%indkk_f2ibz,ikptb,ikptg,&
&             kg,dtorbmag%kgindex,mpi_enreg,npw_kb,pwind_bg,symrec)

             sflag_k=0
             call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icgb,icgg,itrs,job,nband_k,&
&             mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_kb,npw_kg,my_nspinor,&
&             pwind_bg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,usepaw)

             do nn = 1, nband_k
               do n1 = 1, nband_k
                 smat_all(1,nn,n1,ikptb,ikptg) =  smat_kk(1,nn,n1)
                 smat_all(2,nn,n1,ikptb,ikptg) =  smat_kk(2,nn,n1)
                 smat_all(1,n1,nn,ikptg,ikptb) =  smat_kk(1,nn,n1)
                 smat_all(2,n1,nn,ikptg,ikptb) = -smat_kk(2,nn,n1)
               end do
             end do

             has_smat(ikptb,ikptg) = .TRUE.
             has_smat(ikptg,ikptb) = .TRUE.

           end if

           do nn = 1, nband_k
             do n1 = 1, nband_k

               t1A = cmplx(smat_all(1,nn,n1,ikpt,ikptb),smat_all(2,nn,n1,ikpt,ikptb))
               t1B = t1A

               do n2 = 1, nband_k

                 t2A = cmplx(smat_all(1,n1,n2,ikptb,ikptg),smat_all(2,n1,n2,ikptb,ikptg))
                 t3A = conjg(cmplx(smat_all(1,nn,n2,ikpt,ikptg),smat_all(2,nn,n2,ikpt,ikptg)))

                 t2B = conjg(cmplx(smat_all(1,n2,n1,ikpt,ikptb),smat_all(2,n2,n1,ikpt,ikptb)))

                 do n3 = 1, nband_k

                   t3B = cmplx(smat_all(1,n2,n3,ikpt,ikptg),smat_all(2,n2,n3,ikpt,ikptg))
                   t4B=conjg(cmplx(smat_all(1,nn,n3,ikpt,ikptg),smat_all(2,nn,n3,ikpt,ikptg)))

                   tprodB = t1B*t2B*t3B*t4B
                   IB = IB - epsabg*bsigma*gsigma*tprodB/(2.0*deltab*2.0*deltag)
                 end do

                 tprodA = t1A*t2A*t3A
                 IA = IA + epsabg*bsigma*gsigma*tprodA/(2.0*deltab*2.0*deltag)

               end do
             end do
           end do
           
         end do
         

       end do
       
     end do ! end loop over fnkpt
   end do ! end loop over epsabg

   cnum(1,adir) = real(IA+IB)
   cnum(2,adir) = aimag(IA+IB)
 end do ! end loop over adir

 cnum(1,1:3) = MATMUL(gprimd,cnum(1,1:3))
 cnum(2,1:3) = MATMUL(gprimd,cnum(2,1:3))
 dtorbmag%chern(1,1:3) = -cnum(2,1:3)/(two_pi*dtorbmag%fnkpt)
 dtorbmag%chern(2,1:3) =  cnum(1,1:3)/(two_pi*dtorbmag%fnkpt)

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')
 
 write(message,'(a)')' Chern number C from orbital magnetization '
 call wrtout(ab_out,message,'COLL')
 write(message,'(a,a)')'----C is a real vector, given along Cartesian directions----',ch10
 call wrtout(ab_out,message,'COLL')

 do adir = 1, 3
   write(message,'(a,i4,a,2es16.8)')' C(',adir,') : real, imag ',&
&   dtorbmag%chern(1,adir),dtorbmag%chern(2,adir)
   call wrtout(ab_out,message,'COLL')
 end do

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

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
   if (dtset%kptopt /= 3) then
     call pawcprj_free(cprj_ikn)
     call pawcprj_free(cprj_fkn)
     ABI_DATATYPE_DEALLOCATE(cprj_ikn)
     ABI_DATATYPE_DEALLOCATE(cprj_fkn)
   end if
 end if

 ABI_DEALLOCATE(kk_paw)
 ABI_DEALLOCATE(cg1_k)
 ABI_DEALLOCATE(sflag_k)
 ABI_DEALLOCATE(smat_inv)
 ABI_DEALLOCATE(smat_kk)
 ABI_DEALLOCATE(pwind_kb)
 ABI_DEALLOCATE(pwind_kg)
 ABI_DEALLOCATE(pwind_bg)
 ABI_DEALLOCATE(pwnsfac_k)

 ABI_DEALLOCATE(has_smat)
 ABI_DEALLOCATE(smat_all)

end subroutine chern_number
!!***
