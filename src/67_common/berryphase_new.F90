!{\src2tex{textfont=tt}}
!!****f* ABINIT/berryphase_new
!! NAME
!! berryphase_new
!!
!! FUNCTION
!! This routine computes the Berry Phase polarization
!!  and the finite difference expression of the ddk.
!!  [see for example Na Sai et al., PRB 66, 104108 (2002)]
!!
!! COPYRIGHT
!! Copyright (C) 2003-2016 ABINIT  group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!! cg(2,mcg)=planewave coefficients of wavefunctions
!! cprj(natom,mcprj*usecrpj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!! dtfil <type(datafiles_type)>=variables related to files
!! dtset <type(dataset_type)>=all input variables in this dataset
!! psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! gprimd(3,3)=reciprocal space dimensional primitive translations
!! hdr <type(hdr_type)>=the header of wf, den and pot files
!! indlmn(6,lmnmax,ntypat)
!!   array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
!!                                  or i=lmn (if useylm=1)
!! kg(3,mpw*mkmem)=reduced planewave coordinates
!! lmnmax  If useylm=0, max number of (l,m,n) comp. over all type of psps (lnproj)
!!         If useylm=1, max number of (l,n)   comp. over all type of psps (lmnproj)
!! mband=maximum number of bands
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!! mkmem=number of k points treated by this node
!! mpi_enreg=informations about MPI parallelization
!! mpw=maximum dimensioned size of npw
!! my_natom=number of atoms treated by current processor
!! natom=number of atoms in cell
!! nkpt=number of k points
!! npwarr(nkpt)=number of planewaves in basis at this k point
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! ntypat=number of types of atoms in unit cell
!! nkpt=number of k-points
!! option = 1: compute Berryphase polarization
!!          2: compute finite difference expression of the ddk
!!          3: compute polarization & ddk
!! pawrhoij(natom*usepaw) <type(pawrhoij_type)> atomic occupancies
!! pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!! pwind_alloc = first dimension of pwind
!! pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! typat(natom)=type integer for each atom in cell
!! ucvol=unit cell volume in bohr**3.
!! unit_out= unit for output of the results (usually the .out file of ABINIT)
!!   The option unit_out = 0 is allowed. In this case, no information is written
!!   to the output file but only to the log file.
!! usecprj=1 if cprj datastructure has been allocated
!! usepaw= 1: use paw framework. 0:do not use paw.
!! xred(3,natom)=reduced atomic coordinates
!! zion(ntypat)=valence charge of each type of atom
!!
!! OUTPUT
!! ptot(3) = total polarization including correction for jumps
!! red_ptot(3) = total polarization including correction for jumps reduced units
!! pel(3) = reduced coordinates of the electronic polarization (a. u.)
!! pelev(3)= expectation value polarization term (PAW only) in cartesian coordinates (already contained in pel)
!! pion(3)= reduced coordinates of the ionic polarization (a. u.)
!!
!! SIDE EFFECTS
!! Input/Output
!! dtefield <type(efield_type)> = variables related to Berry phase
!!       and electric field calculations (see initberry.f).
!!       In case berryopt = 4, the overlap matrices computed
!!       in this routine are stored in dtefield%smat in order
!!       to be used in the electric field calculation.
!!
!! TODO
!!  - Use the analytical relation between the overlap matrices
!!    S(k,k+dk) and S(k+dk,k) to avoid to recompute them
!!    when ifor = 2.
!!  - change name of option input variable to something more explicit
!!
!! NOTES
!! - pel and pion do not take into account the factor 1/ucvol
!! - In case of a ddk calculation, the eigenvalues are not computed.
!! - The ddk computed by this routine should not be used to
!!   compute the electronic dielectric tensor.
!!
!! PARENTS
!!      elpolariz,update_e_field_vars
!!
!! CHILDREN
!!      appdig,dsdr_k_paw,outwf,pawcprj_alloc,pawcprj_copy,pawcprj_free
!!      pawcprj_get,pawcprj_getdim,pawcprj_mpi_allgather,pawcprj_mpi_recv
!!      pawcprj_mpi_send,pawcprj_put,pawcprj_symkn,pawpolev,polcart,rhophi
!!      smatrix,smatrix_k_paw,wrtout,xmpi_allgather,xmpi_bcast,xmpi_recv
!!      xmpi_send,xmpi_sum,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine berryphase_new(atindx1,cg,cprj,dtefield,dtfil,dtset,psps,&
&  gprimd,hdr,indlmn,kg,lmnmax,mband,mcg,mcprj,&
&  mkmem,mpi_enreg,mpw,my_natom,natom,npwarr,nsppol,ntypat,&
&  nkpt,option,pawrhoij,pawtab,pel,pelev,pion,ptot,red_ptot,pwind,&  !!REC
&  pwind_alloc,pwnsfac,&
&  rprimd,typat,ucvol,unit_out,usecprj,usepaw,xred,zion)

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use defs_wvltypes
 use m_xmpi
 use m_errors
 use m_efield
 use m_profiling_abi

 use m_numeric_tools, only : rhophi
 use m_io_tools, only : open_file
 use m_iowf,     only : outwf
 use m_pawtab,   only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_mpi_allgather, &
&                       pawcprj_put, pawcprj_copy, pawcprj_mpi_recv,  &
&                       pawcprj_mpi_send, pawcprj_free, pawcprj_getdim, pawcprj_symkn

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'berryphase_new'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: lmnmax,mband,mcg,mcprj,mkmem,mpw,my_natom,natom,nkpt
 integer, intent(in) :: nsppol,ntypat,option
 integer, intent(in) :: pwind_alloc,unit_out,usecprj,usepaw
 real(dp), intent(in) :: ucvol
 type(MPI_type), intent(inout) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(efield_type), intent(inout) :: dtefield
 type(hdr_type), intent(inout) :: hdr
!arrays
 integer, intent(in) :: atindx1(natom),indlmn(6,lmnmax,ntypat),kg(3,mpw*mkmem)
 integer, intent(in) :: npwarr(nkpt),pwind(pwind_alloc,2,3)
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: cg(2,mcg),gprimd(3,3)
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
 real(dp), intent(in) :: rprimd(3,3),zion(ntypat)
 real(dp), intent(inout) :: xred(3,natom)
 real(dp), intent(out) :: pel(3),pelev(3),pion(3)
 real(dp), intent(out) :: ptot(3),red_ptot(3) !!REC 
 type(pawrhoij_type), intent(in) :: pawrhoij(my_natom*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
 type(pawcprj_type),intent(in) ::  cprj(natom,mcprj*usecprj)

!Local variables -------------------------
 integer :: count,count1,ddkflag,dest,fdir,unt
 integer :: iatom,iband,icg,icg1,idir,idum,ikpt1i_sp
 integer :: ierr,ifor,ikg,ikpt,ikpt1,ikpt_loc!,ikpt2,ikpt2i,npw_k2, itrs
 integer :: icp1, icp2,icpgr_offset,iproc
! integer :: ii ! appears commented out below in a debug section
 integer :: inibz,ikpt1i
 integer :: isppol,istr,itypat,jband,jkpt,jkstr,job,jsppol
 integer :: maxbd,mcg1_k
 integer :: minbd,my_nspinor,nband_k,ncpgr,nfor,npw_k1,ntotcp,n2dim,nproc,pertcase
 integer :: response,shiftbd,source,spaceComm,tag
 integer :: jj,jstr,kk,ineigh_str
 integer :: istep,jstep,kpt_mark(dtefield%fnkpt),nkstr,nstr,iunmark,berrystep
 integer :: jkpt2, jkpt2i
 real(dp) :: det_mod,dkinv,dphase,dtm_real,dtm_imag,fac,gmod,phase0
 real(dp) :: pol,polbtot,polion,politot,poltot,rho
 logical :: calc_epaw3_force,calc_epaw3_stress,efield_flag
!!REC start
 integer :: jump
 real(dp),save ::pol0(3)
 logical, save :: first=.true.
 logical :: lexist
!!REC end
 real(dp) :: dphase_new,dphase_init
 character(len=fnlen) :: fiwf1o
 character(len=500) :: message
 type(wvl_wf_type) :: wfs
 type(wvl_internal_type) :: wvl
 integer,allocatable :: dimlmn(:),ikpt1_recv(:), sflag_k(:)!,pwind_k(:)
 integer,allocatable :: ikpt3(:), ikpt3i(:), sflag_k_mult(:,:),nattyp_dum(:),npw_k3(:)
 integer,allocatable :: idxkstr_mult(:,:), pwind_k_mult(:,:),itrs_mult(:)
 real(dp) :: det_average(2),dk(3),dtm_k(2),gpard(3),pel_cart(3),pion_cart(3)
 real(dp) :: polb(nsppol),ptot_cart(3),rel_string(2),xcart(3,natom)
 real(dp) :: delta_str(2),dist_,dstr(2)
 real(dp),allocatable :: buffer(:,:),buffer1(:),buffer2(:)
 real(dp),allocatable :: cg1(:,:),cg1_k(:,:),cgq(:,:)
 real(dp),allocatable :: det_string(:,:),dsdr(:,:,:,:,:),dsdr_sum(:,:,:),dsds_sum(:,:,:)
 real(dp),allocatable :: eig_dum(:),epawf3_str(:,:,:),epaws3_str(:,:,:)
 real(dp),allocatable :: occ_dum(:),polberry(:),resid(:),pwnsfac_k(:,:)
 real(dp),allocatable :: smat_inv(:,:,:),smat_k(:,:,:),smat_k_paw(:,:,:)
 real(dp),allocatable :: str_flag(:)
! real(dp),allocatable :: dist_str(:,:),det_string_test(:,:)
 real(dp),allocatable :: dtm_mult(:,:,:), coef(:,:), polb_mult(:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kb(:,:),cprj_buf(:,:),cprj_gat(:,:)
 type(pawcprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:)

! integer :: bband,bbs,bra_start,bra_end,ilmn,ipw,ispinor,jlmn,kband,kbs,ket_start,ket_end,klmn,npw_k
! integer :: nspinor,spnipw,spnshft
! real(dp) :: err_ovlp,mag_ovlp,max_err_ovlp, ovlp_r, ovlp_i, paw_r, paw_i
! real(dp) :: tot_r, tot_i
! real(dp),allocatable :: bra(:,:),ket(:,:)
! complex(dpc) :: cpb,cpk,cterm

!BEGIN TF_CHANGES
 integer :: me
!END TF_CHANGES

!no_abirules

! ***********************************************************************

!DEBUG
!write(std_out,*)' berryphase_new : enter'
!do ii=1,pwind_alloc
!write(std_out,*)ii,pwnsfac(:,ii)
!end do
!stop
!ENDDEBUG

!Init MPI
 spaceComm=mpi_enreg%comm_cell
 nproc=xmpi_comm_size(spaceComm)
 me=mpi_enreg%me_kpt

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

!allocate(pwind_k(mpw))
 ABI_ALLOCATE(pwnsfac_k,(4,mpw))
 ABI_ALLOCATE(sflag_k,(dtefield%nband_occ))
!pwind_k(:) = 0
 pwnsfac_k(1,:) = 1.0_dp ! bra real
 pwnsfac_k(2,:) = 0.0_dp ! bra imag
 pwnsfac_k(3,:) = 1.0_dp ! ket real
 pwnsfac_k(4,:) = 0.0_dp ! ket imag

 if (maxval(dtset%istwfk(:)) /= 1) then
   write(message, '(a,a,a)' )&
&   'This routine does not work yet with istwfk /= 1.',ch10,&
&   'This should have been tested previously ...'
   MSG_BUG(message)
 end if

 if (usepaw == 1 .and. usecprj /= 1) then
   message = ' PAW calculation but cprj datastructure has not been allocated !'
   MSG_BUG(message)
 end if

! useful flags for various efield possibilities 
 efield_flag = ( dtset%berryopt == 4 .or. dtset%berryopt == 6 .or. dtset%berryopt == 7 .or. &
& dtset%berryopt ==14 .or. dtset%berryopt ==16 .or. dtset%berryopt ==17 )
 calc_epaw3_force = ( efield_flag .and. dtset%optforces /= 0 .and. usepaw == 1 )
 calc_epaw3_stress = ( efield_flag .and. dtset%optstress /= 0  .and. usepaw == 1 )

 mcg1_k = mpw*mband
 shiftbd = 1
 if (option > 1) then
   ABI_ALLOCATE(cg1,(2,mcg))
   ABI_ALLOCATE(eig_dum,(2*mband*mband*nkpt*nsppol))
   ABI_ALLOCATE(occ_dum,(mband*nkpt*nsppol))
   eig_dum(:) = zero
   occ_dum(:) = dtefield%sdeg
 end if

!initialize variable tied to multiple step computation
 berrystep=dtset%berrystep
 ABI_ALLOCATE(ikpt3,(berrystep))
 ABI_ALLOCATE(ikpt3i,(berrystep))
 ABI_ALLOCATE(sflag_k_mult,(dtefield%nband_occ,berrystep))
 ABI_ALLOCATE(npw_k3,(berrystep))
 ABI_ALLOCATE(pwind_k_mult,(mpw,berrystep))
 ABI_ALLOCATE(itrs_mult,(berrystep))
 ABI_ALLOCATE(coef,(berrystep,berrystep))
 ABI_ALLOCATE(polb_mult,(nsppol,berrystep))
!coefficient for berryphase computation
 coef(:,:) = 0.0_dp
 do jstep = 1, berrystep
   coef(jstep,1) = 1.d0/real(jstep*jstep,dp)
   if(jstep/=1)coef(jstep,1)=coef(jstep,1)/real(1-jstep*jstep,dp)
 end do
 do istep = 2, berrystep
   do jstep = 1, berrystep
     coef(jstep, istep) = real(istep*istep,dp)*coef(jstep,istep-1)
     if(jstep /= istep)coef(jstep, istep)=coef(jstep,istep)/real(istep*istep-jstep*jstep,dp)
   end do
 end do
!the berryphase using the strings of steps dk, 2*dk, ..., istep*dk is :
!coef(1,istep)*berryphase(dk) + coef(2,istep)*berryphase(2*dk) + ... + coef(istep,istep)*berryphase(istep*dk)
!DEBUG
!write(std_out,*)'coef, sum coef'
!do istep=1,step
!write(std_out,*)coef(:,istep), sum(coef(1:istep,istep))
!end do
!ENDDEBUG

!allocate(dtm(2,dtefield%fnkpt*nsppol))
 ABI_ALLOCATE(dtm_mult,(2,dtefield%fnkpt*nsppol,berrystep))
 ABI_ALLOCATE(cg1_k,(2,mcg1_k))

 if (usepaw == 1) then ! cprj allocation
   ncpgr = cprj(1,1)%ncpgr
   if ( calc_epaw3_force ) then
     ABI_ALLOCATE(dsdr_sum,(natom,3,dtefield%fnkpt*nsppol))
     ABI_ALLOCATE(epawf3_str,(natom,3,3))
   end if 
   if ( calc_epaw3_stress ) then
     ABI_ALLOCATE(dsds_sum,(natom,6,dtefield%fnkpt*nsppol))
     ABI_ALLOCATE(epaws3_str,(natom,3,6))
   end if 
   ABI_ALLOCATE(dimlmn,(natom))
   call pawcprj_getdim(dimlmn,natom,nattyp_dum,ntypat,typat,pawtab,'R')
   ABI_DATATYPE_ALLOCATE(cprj_k,(natom,dtefield%nspinor*mband))
   ABI_DATATYPE_ALLOCATE(cprj_kb,(natom,dtefield%nspinor*mband))
   ABI_DATATYPE_ALLOCATE(cprj_gat,(natom,nproc*dtefield%nspinor*mband))
   call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
   call pawcprj_alloc(cprj_kb,ncpgr,dimlmn)
   call pawcprj_alloc(cprj_gat,ncpgr,dimlmn)
   if (dtset%kptopt /= 3) then
     ABI_DATATYPE_ALLOCATE(cprj_ikn,(natom,dtefield%nspinor*mband))
     ABI_DATATYPE_ALLOCATE(cprj_fkn,(natom,dtefield%nspinor*mband))
     call pawcprj_alloc(cprj_ikn,ncpgr,dimlmn)
     call pawcprj_alloc(cprj_fkn,ncpgr,dimlmn)
   end if

   n2dim = dtefield%nspinor*mband
   ntotcp = n2dim*SUM(dimlmn(:))
   if (nproc>1) then
     ABI_DATATYPE_ALLOCATE(cprj_buf,(natom,dtefield%nspinor*mband))
     call pawcprj_alloc(cprj_buf,ncpgr,dimlmn)
   end if

   if ( efield_flag ) then
     write(message,'(2a,i5,2a)')ch10,&
&     ' nkpt = ',nkpt,ch10,' copy cprj to dtefield%cprj '
     call wrtout(std_out,message,'COLL')

     do isppol = 1, nsppol

       ikpt_loc = 0
       ikpt1 = 0
       do while (ikpt_loc < mkmem)

         if (ikpt_loc < mkmem) ikpt1 = ikpt1 + 1
         if ((ikpt1 > nkpt).and.(ikpt_loc < mkmem)) exit
         nband_k = dtset%nband(ikpt1) 

         if ( (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt1,1,nband_k,isppol,me)) .and. &
&         (ikpt_loc <= mkmem) ) cycle

         ikpt_loc = ikpt_loc + 1

         ABI_ALLOCATE(ikpt1_recv,(nproc))
         call xmpi_allgather(ikpt1,ikpt1_recv,spaceComm,ierr)
         call pawcprj_get(atindx1,cprj_k,cprj,natom,1,(ikpt_loc-1)*nband_k*my_nspinor,ikpt1,0,isppol,mband,&
&         mkmem,natom,nband_k,nband_k,my_nspinor,nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         call pawcprj_mpi_allgather(cprj_k,cprj_gat,natom,n2dim,dimlmn,ncpgr,nproc,spaceComm,ierr,rank_ordered=.true.)
         do iproc = 1, nproc
           icp2=nband_k*(iproc-1)*my_nspinor
           call pawcprj_get(atindx1,cprj_k,cprj_gat,natom,1,icp2,ikpt1,0,isppol,mband,&
&           nproc,natom,nband_k,nband_k,my_nspinor,1,0,&
&           mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
           icp1 = nband_k*(ikpt1_recv(iproc)-1)*my_nspinor
           call pawcprj_put(atindx1,cprj_k,dtefield%cprj,natom,1,icp1,ikpt1,0,isppol,&
&           mband,dtefield%fnkpt,natom,nband_k,nband_k,dimlmn,my_nspinor,nsppol,0,&
&           mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         end do
         ABI_DEALLOCATE(ikpt1_recv)

       end do ! close loop over k-points

     end do ! end loop over nsppol

   end if ! end check on efield
 end if

!!=======================================
!! code to test orthonormality of cg_k
!!=======================================
!
!ikpt = 3
!npw_k = npwarr(ikpt)
!isppol = 1
!nband_k = dtefield%nband_occ
!ABI_ALLOCATE(bra,(2,npw_k*my_nspinor))
!ABI_ALLOCATE(ket,(2,npw_k*my_nspinor))
!max_err_ovlp=0.0
!call pawcprj_get(atindx1,cprj_k,cprj,natom,1,dtefield%cprjindex(ikpt,isppol),ikpt,0,isppol,mband,&
!&         mkmem,natom,nband_k,nband_k,my_nspinor,nsppol,0)
!do bband = 1, nband_k
!bra_start = dtefield%cgindex(ikpt,nsppol)+1+(bband-1)*npw_k*my_nspinor
!bra_end = bra_start + npw_k*my_nspinor - 1
!bra(1:2,1:npw_k*my_nspinor) = cg(1:2,bra_start:bra_end)
!do kband = 1, nband_k
!ket_start = dtefield%cgindex(ikpt,nsppol)+1+(kband-1)*npw_k*my_nspinor
!ket_end = ket_start + npw_k*my_nspinor - 1
!ket(1:2,1:npw_k*my_nspinor) = cg(1:2,ket_start:ket_end)
!
!tot_r = 0.0; tot_i = 0.0     
!do ispinor = 1, my_nspinor
!ovlp_r = 0.0; ovlp_i = 0.0
!spnshft = (ispinor-1)*npw_k
!do ipw = 1, npw_k
!spnipw = ipw + spnshft
!ovlp_r = ovlp_r + bra(1,spnipw)*ket(1,spnipw)+bra(2,spnipw)*ket(2,spnipw)
!ovlp_i = ovlp_i - bra(2,spnipw)*ket(1,spnipw)+bra(1,spnipw)*ket(2,spnipw)
!end do ! end loop over ipw
!paw_r = 0.0; paw_i = 0.0
!do iatom = 1, natom
!itypat = typat(iatom)
!do ilmn = 1, dtefield%lmn_size(itypat)
!do jlmn = 1, dtefield%lmn_size(itypat)
!klmn=max(ilmn,jlmn)*(max(ilmn,jlmn)-1)/2 + min(ilmn,jlmn)
!bbs = my_nspinor*(bband-1)+ispinor
!kbs = my_nspinor*(kband-1)+ispinor
!cpb=cmplx(cprj_k(iatom,bbs)%cp(1,ilmn),cprj_k(iatom,bbs)%cp(2,ilmn))
!cpk=cmplx(cprj_k(iatom,kbs)%cp(1,jlmn),cprj_k(iatom,kbs)%cp(2,jlmn))
!cterm = conjg(cpb)*pawtab(itypat)%sij(klmn)*cpk
!paw_r = paw_r + real(cterm)
!paw_i = paw_i + aimag(cterm)
!end do ! end loop over jlmn
!end do ! end loop over ilmn
!end do ! end loop over iatom
!tot_r = tot_r + ovlp_r + paw_r
!tot_i = tot_i + ovlp_i + paw_i
!end do ! end loop over ispinor
!
!!     write(std_out,'(a,2i4,2es16.8)')' JWZ Debug: berryphase_new bband kband ovlp : ',&
!!&           bband,kband,tot_r,tot_i
!mag_ovlp =  tot_r*tot_r + tot_i*tot_i
!if(bband==kband) then
!err_ovlp=abs(mag_ovlp-1.0)
!else
!err_ovlp=abs(mag_ovlp)
!end if
!max_err_ovlp=MAX(max_err_ovlp,err_ovlp)
!end do ! end loop over kband
!end do ! end loop over bband
!write(std_out,'(a,i4,es16.8)')' JWZ Debug: berrphase_new ikpt ovlp err : ',&
!&           ikpt,max_err_ovlp
!ABI_DEALLOCATE(bra)
!ABI_DEALLOCATE(ket)
!
!!=========================================
!! end code to test orthonormality of cg_k
!!=========================================

 pel(:) = zero ; pelev(:)=zero ; pion(:) = zero ; ptot(:)=zero ; red_ptot(:)=zero

 minbd = 1   ;  maxbd = dtefield%nband_occ

 if(calc_epaw3_force) dtefield%epawf3(:,:,:) = zero
 if(calc_epaw3_stress) dtefield%epaws3(:,:,:) = zero

 do idir = 1, 3

!  dtm(:,:) = zero
   dtm_mult(:,:,:) = zero
   if (calc_epaw3_force) dsdr_sum(:,:,:) = zero
   if (calc_epaw3_stress) dsds_sum(:,:,:) = zero

   if (dtset%rfdir(idir) == 1) then

     if (abs(dtefield%efield_dot(idir)) < tol12) dtefield%sflag(:,:,:,idir) = 0

!    Check whether the polarization or the ddk must be computed

!    nfor = 1 : to compute P, I only need the WF at k + dk
!    nfor = 2 : to compute the ddk, I need the WF at k + dk and k - dk
!    dkinv    : +-1/2dk

     if (option > 1) then

       ddkflag = 1
       nfor = 2
       job = 1
       cg1(:,:) = zero
       if (option == 3) job = 11

     else if (option == 1) then

       ddkflag = 0
       nfor = 1
       job = 10
! electric fields with PAW also needs S_inverse for forces and stresses
       if (calc_epaw3_force .or. calc_epaw3_stress) job = 11

     end if

     dk(:) = dtefield%dkvecs(:,idir)
     gpard(:) = dk(1)*gprimd(:,1) + dk(2)*gprimd(:,2) + dk(3)*gprimd(:,3)
     gmod = sqrt(dot_product(gpard,gpard))
     if (option > 1) dkinv = one/(two*dk(idir))

     write(message,'(a,a,a,3f9.5,a,a,3f9.5,a)')ch10,&
&     ' Computing the polarization (Berry phase) for reciprocal vector:',ch10,&
&     dk(:),' (in reduced coordinates)',ch10,&
&     gpard(1:3),' (in cartesian coordinates - atomic units)'
     call wrtout(std_out,message,'COLL')
     if (unit_out /= 0) then
       call wrtout(unit_out,message,'COLL')
     end if

     write(message,'(a,i5,a,a,i5)')&
&     ' Number of strings: ',dtefield%nstr(idir),ch10,&
&     ' Number of k points in string:', dtefield%nkstr(idir)
     call wrtout(std_out,message,'COLL')
     if (unit_out /= 0) then
       call wrtout(unit_out,message,'COLL')
     end if

     if ((option == 2).or.(option == 3)) then

       write(message,'(a,a,a,3f9.5,a,a,3f9.5,a)')ch10,&
&       ' Computing the ddk (Berry phase) for reciprocal vector:',ch10,&
&       dk(:),' (in reduced coordinates)',ch10,&
&       gpard(1:3),' (in cartesian coordinates - atomic units)'
       call wrtout(std_out,message,'COLL')
       if (unit_out /= 0) then
         call wrtout(unit_out,message,'COLL')
       end if

     end if

     do ifor = 1, nfor

       if (ifor == 2) then
         dk(:) = -1_dp*dk(:)
         job = 1   ! only the inverse of the overlap matrix is required
         dkinv = -1_dp*dkinv
       end if


!      Compute the determinant and/or the inverse of the overlap matrix
!      for each pair of k-points < u_nk | u_nk+dk >

       icg = 0 ; icg1 = 0
       ABI_ALLOCATE(smat_k,(2,dtefield%nband_occ,dtefield%nband_occ))
       ABI_ALLOCATE(smat_inv,(2,dtefield%nband_occ,dtefield%nband_occ))
       ABI_ALLOCATE(smat_k_paw,(2,usepaw*dtefield%nband_occ,usepaw*dtefield%nband_occ))
       if (calc_epaw3_force .or. calc_epaw3_stress) then ! dsdr needed for forces and stresses in electric field with PAW
         ABI_ALLOCATE(dsdr,(2,natom,ncpgr,usepaw*dtefield%nband_occ,usepaw*dtefield%nband_occ))
         dsdr = zero
       end if
       

!      Loop on the values of ikpt_loc and ikpt1 :
!      ikpt1 is incremented one by one, and number the k points in the FBZ
!      ikpt1i refer to the k point numbering in the IBZ
!      ikpt_loc differs from ikpt1 only in the parallel case, and gives
!      the index of the k point in the FBZ, in the set treated by the present processor
!      NOTE : in order to allow synchronisation, ikpt_loc contain information about
!      ikpt AND ISPPOL !
!      It means that the following loop is equivalent to a double loop :
!      do isppol = 1, nsppol
!      do ikpt1 =  1, dtefield%fmkmem
!      
       do ikpt_loc = 1, dtefield%fmkmem_max*nsppol

         ikpt1=mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc,1)
         isppol=mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc,2)

         if (ikpt1 > 0 .and. isppol > 0) then

           ikpt1i = dtefield%indkk_f2ibz(ikpt1,1)
           nband_k = dtset%nband(ikpt1i + (isppol-1)*dtset%nkpt)

!          DEBUG
!          Please keep this debugging feature
!          write(std_out,'(a,5i4)' )' berryphase_new : ikpt_loc,ikpt1,isppol,idir,ifor=',&
!          &                                  ikpt_loc,ikpt1,isppol,idir,ifor
!          ENDDEBUG

           inibz=0
           if (dtset%kptns(1,ikpt1i) == dtefield%fkptns(1,ikpt1) .and. &
&           dtset%kptns(2,ikpt1i) == dtefield%fkptns(2,ikpt1) .and. &
&           dtset%kptns(3,ikpt1i) == dtefield%fkptns(3,ikpt1)) inibz=1

           ikg = dtefield%fkgindex(ikpt1)
!          ikpt2 = dtefield%ikpt_dk(ikpt1,ifor,idir)
!          ikpt2i = dtefield%indkk_f2ibz(ikpt2,1)

!          ikpt3(istep) : index of kpt1 + istep*dk in the FBZ
!          ikpt3i(istep) : index of kpt1 + istep*dk in the IBZ
           ikpt3(1) = dtefield%ikpt_dk(ikpt1,ifor,idir)
           ikpt3i(1) = dtefield%indkk_f2ibz(ikpt3(1),1)
           do istep = 1, berrystep-1
             ikpt3(istep+1) = dtefield%ikpt_dk(ikpt3(istep),ifor,idir)
             ikpt3i(istep+1) = dtefield%indkk_f2ibz(ikpt3(istep+1),1)
           end do

!          itrs = 0
!          if (dtefield%indkk_f2ibz(ikpt1,6) == 1 ) itrs = itrs + 1
!          if (dtefield%indkk_f2ibz(ikpt2,6) == 1 ) itrs = itrs + 10

           itrs_mult(:)=0
           if (dtefield%indkk_f2ibz(ikpt1,6) == 1 ) itrs_mult(:) = itrs_mult(:) + 1
           do istep=1,berrystep
             if (dtefield%indkk_f2ibz(ikpt3(istep),6) == 1 ) itrs_mult(istep) = itrs_mult(istep) + 10
           end do

           npw_k1 = npwarr(ikpt1i)
!          npw_k2 = npwarr(ikpt2i)

           do istep = 1, berrystep
             npw_k3(istep)=npwarr(ikpt3i(istep))
           end do

!          ji: the loop is over the FBZ, but sflag and smat only apply to the IBZ
           if ( efield_flag .and. inibz == 1) then  !!HONG
             ikpt1i_sp=ikpt1i+(isppol-1)*dtset%nkpt
             smat_k(:,:,:) = dtefield%smat(:,:,:,ikpt1i_sp,ifor,idir)
           else
             smat_k(:,:,:) = zero
           end if

!          pwind_k(1:npw_k1) = pwind(ikg+1:ikg+npw_k1,ifor,idir)
           pwnsfac_k(1,1:npw_k1) = pwnsfac(1,ikg+1:ikg+npw_k1)
           pwnsfac_k(2,1:npw_k1) = pwnsfac(2,ikg+1:ikg+npw_k1)

!          the array needed to compute the overlap matrix between k and k+istep*dk (with multiple steps)
!          the 0-case (no corresponding pw in k and k+dk) could be handled better (k+2*dk could have a corresponding pw ?)
           pwind_k_mult(1:npw_k1,1)=pwind(ikg+1:ikg+npw_k1,ifor,idir)
           do istep = 1, berrystep-1
             do jj=1, npw_k1
               if(pwind_k_mult(jj,istep)/=0)then
                 pwind_k_mult(jj,istep+1) = pwind(dtefield%fkgindex(ikpt3(istep))+pwind_k_mult(jj,istep),ifor,idir)
               else
                 pwind_k_mult(jj,istep+1) = 0
               end if
             end do
           end do

!          DEBUG
!          write(std_out,*)' berryphase_new : dtset%berryopt,inibz,ikpt1i,isppol,dtset%nkpt,ifor,idir', &
!          &          dtset%berryopt,inibz,ikpt1i,isppol,dtset%nkpt,ifor,idir
!          write(std_out,'(a,4i4)' )' berryphase_new : sflag_k(:)=',sflag_k(:)
!          ENDDEBUG

           if ( efield_flag .and. inibz == 1) then  !!HONG
             ikpt1i_sp=ikpt1i+(isppol-1)*dtset%nkpt
             sflag_k(:) = dtefield%sflag(:,ikpt1i_sp,ifor,idir)
           else
             sflag_k(:) = 0
           end if

           if (usepaw == 1) then
             icp1=dtefield%cprjindex(ikpt1i,isppol)
             call pawcprj_get(atindx1,cprj_k,cprj,natom,1,icp1,ikpt1i,0,isppol,&
&             mband,mkmem,natom,dtefield%nband_occ,dtefield%nband_occ,&
&             my_nspinor,nsppol,0,mpicomm=mpi_enreg%comm_kpt,&
&             proc_distrb=mpi_enreg%proc_distrb)

             if ( ikpt1i /= ikpt1 ) then
               call pawcprj_copy(cprj_k,cprj_ikn)
               call pawcprj_symkn(cprj_fkn,cprj_ikn,dtefield%atom_indsym,dimlmn,-1,indlmn,&
&               dtefield%indkk_f2ibz(ikpt1,2),dtefield%indkk_f2ibz(ikpt1,6),&
&               dtefield%fkptns(:,dtefield%i2fbz(ikpt1i)),&
&               dtefield%lmax,dtefield%lmnmax,mband,natom,dtefield%nband_occ,my_nspinor,&
&               dtefield%nsym,ntypat,typat,dtefield%zarot)
               call pawcprj_copy(cprj_fkn,cprj_k)
             end if

           end if ! end if usepaw

!          DEBUG
!          write(std_out,'(a,4i4)' )' berryphase_new : sflag_k(:)=',sflag_k(:)
!          ENDDEBUG
           
!          DEBUG
!          write(std_out,'(a,7i4)')'me, idir,ifor, ikpt_loc, ikpt1, isppol = ',&
!          & me,idir,ifor,ikpt_loc,ikpt1,isppol
!          write(std_out,'(a,10i3)')'pwind_k(1:10) = ',pwind_k(1:10)
!          ENDDEBUG
           
           do istep=1,berrystep
             sflag_k_mult(:,istep) = sflag_k(:)
           end do
           
         end if ! end check that ikpt1 > 0 and isppol > 0

!        --------------------------------------------------------------------------------
!        Communication
!        --------------------------------------------------------------------------------

         do istep=1,berrystep

!          if(ikpt_loc <= nsppol*dtefield%fmkmem) then
           if (ikpt1 > 0 .and. isppol > 0) then ! I currently have a true kpt to use

             count = npw_k3(istep)*my_nspinor*nband_k
             ABI_ALLOCATE(cgq,(2,count))
             cgq = zero
             source = me
             if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt3i(istep),1,nband_k,isppol,me)) then
!              I need the datas from someone else
               source = mpi_enreg%proc_distrb(ikpt3i(istep),1,isppol)
             end if
           else
             source = -1 ! I do not have a kpt to use
           end if

           do dest = 0, nproc-1

             if ((dest==me) .and. (ikpt1>0) .and. (isppol>0)) then
!              I am destination and I have something to do
!              if (mpi_enreg%paral_compil_kpt == 1) write(std_out,*) &
!              &               'coucou 2, mpi_enreg%proc_distrb(ikpt3i(istep),1:nband_k,isppol) : ', &
!              &               mpi_enreg%proc_distrb(ikpt3i(istep),1:nband_k,isppol)
!              write(std_out,*)'ikpt3i(istep) ', ikpt3i(istep)
!              write(std_out,*)'nband_k ',nband_k
!              write(std_out,*)'isppol ', isppol
!              write(std_out,*)'mpi_enreg%proc_distrb',mpi_enreg%proc_distrb

               if (source == me) then
!                I am destination and source
!                DEBUG
!                write(std_out,*)'copying ... '
!                write(std_out,*)'me: ',me, 'ikpt3i(istep) ', ikpt3i(istep), 'isppol ', isppol
!                ENDDEBUG

!                pwnsfac
                 idum = dtefield%fkgindex(ikpt3(istep))
                 pwnsfac_k(3,1:npw_k3(istep)) = pwnsfac(1,idum+1:idum+npw_k3(istep))
                 pwnsfac_k(4,1:npw_k3(istep)) = pwnsfac(2,idum+1:idum+npw_k3(istep))

!                cgq (and cprj)
                 icg1 = dtefield%cgindex(ikpt3i(istep),isppol)

                 if (usepaw == 1) then
                   icp2=dtefield%cprjindex(ikpt3i(istep),isppol)
                   call pawcprj_get(atindx1,cprj_kb,cprj,natom,1,icp2,ikpt3i(istep),0,isppol,&
&                   mband,mkmem,natom,dtefield%nband_occ,dtefield%nband_occ,my_nspinor,&
&                   nsppol,0,mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
                 end if

                 cgq(:,1:count)  = cg(:,icg1+1:icg1+count)
!                if (usepaw == 1) then
!                  call pawcprj_copy(cprj_buf,cprj_kb)
!                end if

!                if ((source /= me)) then
               else
!                I am the destination but not the source -> receive
!                DEBUG
!                write(std_out,'(a)')'receiving ...'
!                write(std_out,'(a,i4,a,i4,a,i4,a,i4)')'me: ',me, 'source ', source,'ikpt3i(istep) ', ikpt3i(istep), 'isppol ', isppol
!                ENDDEBUG

!                receive pwnsfac
                 ABI_ALLOCATE(buffer,(2,npw_k3(istep)))
                 tag = ikpt3(istep) + (isppol - 1)*dtefield%fnkpt
                 call xmpi_recv(buffer,source,tag,spaceComm,ierr)
                 pwnsfac_k(3,1:npw_k3(istep)) = buffer(1,1:npw_k3(istep))
                 pwnsfac_k(4,1:npw_k3(istep)) = buffer(2,1:npw_k3(istep))
                 ABI_DEALLOCATE(buffer)

!                receive cgq (and cprj)
                 tag = ikpt3i(istep) + (isppol - 1)*nkpt
                 call xmpi_recv(cgq,source,tag,spaceComm,ierr)

                 if (usepaw == 1) then
                   call pawcprj_mpi_recv(natom,n2dim,dimlmn,ncpgr,cprj_kb,source,spaceComm,ierr)
                 end if

               end if

             else if (dest /= me) then 

!              jkpt is the kpt which is being treated by dest
!              jsppol is his isppol
               jkpt = mpi_enreg%kpt_loc2fbz_sp(dest, ikpt_loc,1)
               jsppol = mpi_enreg%kpt_loc2fbz_sp(dest, ikpt_loc,2)

               if (jkpt > 0 .and. jsppol > 0) then ! dest is treating a true kpt

                 jkpt2 = dtefield%ikpt_dk(jkpt,ifor,idir)
                 jkpt2i = dtefield%indkk_f2ibz(jkpt2,1)

!                check if I am his source
                 if((mpi_enreg%proc_distrb(jkpt2i,1,jsppol) == me))  then
!                  I know something about jkpt3i and I must send it
!                  DEBUG
!                  write(std_out,'(a)')'sending ...'
!                  write(std_out,'(a,i4,a,i4,a,i4,a,i4)')'dest: ',dest,' me: ',me,&
!                  &                          ' jkpt2i ',jkpt2i,' jsppol: ',jsppol
!                  ENDDEBUG

!                  pwnsfac
                   tag = jkpt2 + (jsppol - 1)*dtefield%fnkpt
                   count1 = npwarr(jkpt2i)
                   ABI_ALLOCATE(buffer,(2,count1))
                   idum = dtefield%fkgindex(jkpt2)
                   buffer(1,1:count1)  = pwnsfac(1,idum+1:idum+count1)
                   buffer(2,1:count1)  = pwnsfac(2,idum+1:idum+count1)
                   call xmpi_send(buffer,dest,tag,spaceComm,ierr)
                   ABI_DEALLOCATE(buffer)

!                  cgq (and cprj)
                   icg1 = dtefield%cgindex(jkpt2i,jsppol)

                   if (usepaw == 1) then
                     icp2=dtefield%cprjindex(jkpt2i,jsppol)
                     call pawcprj_get(atindx1,cprj_buf,cprj,natom,1,icp2,jkpt2i,0,jsppol,&
&                     mband,mkmem,natom,dtefield%nband_occ,dtefield%nband_occ,&
&                     my_nspinor,nsppol,0,mpicomm=mpi_enreg%comm_kpt,&
&                     proc_distrb=mpi_enreg%proc_distrb)
                   end if

                   tag = jkpt2i + (jsppol - 1)*nkpt
                   count1 = npwarr(jkpt2i)*my_nspinor*nband_k
                   ABI_ALLOCATE(buffer,(2,count1))
                   buffer(:,1:count1)  = cg(:,icg1+1:icg1+count1)
                   call xmpi_send(buffer,dest,tag,spaceComm,ierr)
                   ABI_DEALLOCATE(buffer)

                   if (usepaw == 1 ) then
                     call pawcprj_mpi_send(natom,n2dim,dimlmn,ncpgr,cprj_buf,dest,spaceComm,ierr)
                   end if

                 end if ! end check that I am his source

               end if ! end check that jkpt > 0 and jsppol > 0

             end if ! end if statements on dest == me or dest /= me

           end do  ! end loop over dest = 0, nproc - 1

           if (ikpt1 > 0 .and. isppol > 0) then ! if I am treating a kpt, compute the smatrix

             if (usepaw == 1) then
               if (ikpt3(istep) /= ikpt3i(istep)) then ! cprj_kb refers to ikpt3i(istep), must compute ikpt3(istep) value
                 call pawcprj_copy(cprj_kb,cprj_ikn)

                 call pawcprj_symkn(cprj_fkn,cprj_ikn,dtefield%atom_indsym,dimlmn,-1,indlmn,&
&                 dtefield%indkk_f2ibz(ikpt3(istep),2),dtefield%indkk_f2ibz(ikpt3(istep),6),&
&                 dtefield%fkptns(:,dtefield%i2fbz(ikpt3i(istep))),&
&                 dtefield%lmax,dtefield%lmnmax,mband,natom,&
&                 dtefield%nband_occ,my_nspinor,dtefield%nsym,ntypat,typat,&
&                 dtefield%zarot)
                 call pawcprj_copy(cprj_fkn,cprj_kb)
               end if
               call smatrix_k_paw(cprj_k,cprj_kb,dtefield,idir,ifor,mband,natom,smat_k_paw,typat)
!              write(std_out,'(a,5i4)')' JWZ berryphase_new : ikpt_loc,ikpt1,ikpt1i,ikpt2,ikpt2i ',ikpt_loc,ikpt1,ikpt1i,ikpt3(istep),ikpt3i(istep)
!              call smatrix_k0_paw(atindx1,cprj_k,cprj_k,dtefield,ikpt1i,idir,ifor,&
!              &                                  mband,mpi_enreg,natom,ntypat,pawtab,smat_k_paw,typat)
               if (calc_epaw3_force .or. calc_epaw3_stress) then
                 call dsdr_k_paw(cprj_k,cprj_kb,dsdr,dtefield,idir,ifor,mband,natom,ncpgr,typat)
               end if
             end if

             icg1 = 0
             icg = dtefield%cgindex(ikpt1i,isppol)
!            DEBUG
!            if(istep<=2)then
!            if(ikpt1==1)then
!            write(std_out,'(a,2i4,3e15.4)')'istep ikpt3, kpt, cgq', istep, ikpt3(istep), dtefield%fkptns(:,ikpt3(istep))
!            write(std_out,*) cgq
!            write(std_out,*)
!            end if
!            end if
!            ENDDEBUG
             call smatrix(cg,cgq,cg1_k,ddkflag,dtm_k,icg,icg1,itrs_mult(istep),job,maxbd,&
&             mcg,count,mcg1_k,minbd,&
&             mpw,dtefield%nband_occ,&
&             npw_k1,npw_k3(istep),my_nspinor,pwind_k_mult(:,istep),pwnsfac_k,sflag_k_mult(:,istep),&
&             shiftbd,smat_inv,smat_k,smat_k_paw,usepaw)
             
! in finite electric field case with paw must save additional F3 term in forces 
             if(calc_epaw3_force) then
! when ncpgr = 3, gradients are wrt to atom displacements
! but when ncpgr = 9, first 6 gradients are wrt strains, last three are displacements
               icpgr_offset = 0
               if (ncpgr == 9) icpgr_offset = 6
               do iatom = 1, natom
                 do fdir = 1, 3
                   dsdr_sum(iatom,fdir,ikpt1+(isppol-1)*dtefield%fnkpt) = zero
                   do iband = 1, dtefield%nband_occ
                     do jband = 1, dtefield%nband_occ
! collect Im{Trace{S^{-1}.dS/dR}} for this k point
                       dsdr_sum(iatom,fdir,ikpt1+(isppol-1)*dtefield%fnkpt) = &
&                       dsdr_sum(iatom,fdir,ikpt1+(isppol-1)*dtefield%fnkpt) + &
&                       smat_inv(2,iband,jband)*dsdr(1,iatom,icpgr_offset+fdir,jband,iband) + &
&                       smat_inv(1,iband,jband)*dsdr(2,iatom,icpgr_offset+fdir,jband,iband)
                     end do ! end sum over jband
                   end do ! end sum over iband
                 end do ! end sum over fdir
               end do ! end sum over iatom
             end if ! end check on calc_epaw3_force
! in finite electric field case with paw must save additional F3 term in stress
! note that when strains are present they are always saved before forces
! therefore no need for icpgr_offset in this case
             if(calc_epaw3_stress) then
               do iatom = 1, natom
                 do fdir = 1, 6
                   dsds_sum(iatom,fdir,ikpt1+(isppol-1)*dtefield%fnkpt) = zero
                   do iband = 1, dtefield%nband_occ
                     do jband = 1, dtefield%nband_occ
! collect Im{Trace{S^{-1}.dS/de}} for this k point
                       dsds_sum(iatom,fdir,ikpt1+(isppol-1)*dtefield%fnkpt) = &
&                       dsds_sum(iatom,fdir,ikpt1+(isppol-1)*dtefield%fnkpt) + &
&                       smat_inv(2,iband,jband)*dsdr(1,iatom,fdir,jband,iband) + &
&                       smat_inv(1,iband,jband)*dsdr(2,iatom,fdir,jband,iband)
                     end do ! end sum over jband
                   end do ! end sum over iband
                 end do ! end sum over fdir
               end do ! end sum over iatom
             end if ! end check on calc_epaw3_stress               

             if ((job == 10).or.(job == 11)) then

               if (sqrt(dtm_k(1)*dtm_k(1) + dtm_k(2)*dtm_k(2)) < tol12) then
                 write(message,'(a,i5,a,a,a)')&
&                 '  For k-point #',ikpt1,',',ch10,&
&                 '  the determinant of the overlap matrix is found to be 0.'
                 MSG_BUG(message)
               end if

               dtm_mult(1,ikpt1+(isppol-1)*dtefield%fnkpt,istep) = dtm_k(1)
               dtm_mult(2,ikpt1+(isppol-1)*dtefield%fnkpt,istep) = dtm_k(2)

             end if

             if ( efield_flag .and. inibz == 1 .and. istep == 1)  then  !!HONG
               ikpt1i_sp=ikpt1i+(isppol-1)*dtset%nkpt
               dtefield%smat(:,:,:,ikpt1i_sp,ifor,idir) = &
&               smat_k(:,:,:)
               dtefield%sflag(:,ikpt1i_sp,ifor,idir) = &
&               sflag_k_mult(:,1)
             end if

             if ((option>1.and.((job == 1).or.(job == 11))).and.(inibz == 1) .and. istep == 1) then
               cg1(:,icg + 1: icg + npw_k1*my_nspinor*nband_k) = &
               cg1(:,icg + 1:icg + npw_k1*my_nspinor*nband_k) + &
               dkinv*cg1_k(:,1:npw_k1*my_nspinor*nband_k)
             end if

             ABI_DEALLOCATE(cgq)

           end if ! end if ikpt1 > 0 and isppol > 0

         end do ! end loop over istep

!        if (ikpt_loc <= dtefield%fmkmem) sflag_k(:) = sflag_k_mult(:,1)
         if (ikpt1 > 0) sflag_k(:) = sflag_k_mult(:,1)

       end do ! close loop over ikpt_loc (k-points, isppol)

       ABI_DEALLOCATE(smat_inv)
       ABI_DEALLOCATE(smat_k)
       ABI_DEALLOCATE(smat_k_paw)
       if (calc_epaw3_force .or. calc_epaw3_stress) then
         ABI_DEALLOCATE(dsdr)
       end if

     end do   ! close loop over ifor

     if (nproc>1) then
       count = 2*dtefield%fnkpt*nsppol*berrystep
       ABI_ALLOCATE(buffer1,(count))
       ABI_ALLOCATE(buffer2,(count))
       buffer1(:) = reshape(dtm_mult(:,:,:),(/count/))
       call xmpi_sum(buffer1,buffer2,count,spaceComm,ierr)
       dtm_mult(:,:,:) = reshape(buffer2(:),(/2,dtefield%fnkpt*nsppol,berrystep/))
       ABI_DEALLOCATE(buffer1)
       ABI_DEALLOCATE(buffer2)
       if (calc_epaw3_force) then
         count = natom*3*dtefield%fnkpt*nsppol
         ABI_ALLOCATE(buffer1,(count))
         ABI_ALLOCATE(buffer2,(count))
         buffer1(:) = reshape(dsdr_sum(:,:,:),(/count/))
         call xmpi_sum(buffer1,buffer2,count,spaceComm,ierr)
         dsdr_sum(:,:,:) = reshape(buffer2(:),(/natom,3,dtefield%fnkpt*nsppol/))
         ABI_DEALLOCATE(buffer1)
         ABI_DEALLOCATE(buffer2)
       end if 
       if (calc_epaw3_stress) then
         count = natom*6*dtefield%fnkpt*nsppol
         ABI_ALLOCATE(buffer1,(count))
         ABI_ALLOCATE(buffer2,(count))
         buffer1(:) = reshape(dsds_sum(:,:,:),(/count/))
         call xmpi_sum(buffer1,buffer2,count,spaceComm,ierr)
         dsds_sum(:,:,:) = reshape(buffer2(:),(/natom,6,dtefield%fnkpt*nsppol/))
         ABI_DEALLOCATE(buffer1)
         ABI_DEALLOCATE(buffer2)
       end if 
     end if

!    DEBUG
!    write(std_out,*)
!    write(std_out,*)'istep = 1, nsppol =',nsppol
!    istep=1
!    isppol=1
!    do jkpt = 1, dtefield%fnkpt
!    write(std_out,'(a,i4,3e15.4,2e15.4)')'jkpt, kpt, dtm_mult(:,kpt,1)', jkpt, dtefield%fkptns(:,jkpt),  dtm_mult(:,jkpt+(isppol-1)*dtefield%fnkpt,istep)
!    end do
!    write(std_out,*)
!    write(std_out,*) "istep = 2"
!    if(berrystep>=2)then
!    istep=2
!    isppol=1
!    do jkpt = 1, dtefield%fnkpt
!    write(std_out,'(a,i4,3e15.4,2e15.4)')'jkpt, kpt, dtm_mult(:,kpt,2)', jkpt, dtefield%fkptns(:,jkpt),  dtm_mult(:,jkpt+(isppol-1)*dtefield%fnkpt,istep)
!    end do
!    end if
!    ENDDEBUG

!    ===========================================================================

!    Compute the Berry phase polarization

     if ((option == 1).or.(option == 3)) then

!      Compute the electronic Berry phase

       polb_mult(:,:)=zero
       do istep = 1,berrystep

         if(istep /= 1) then
!          construct the strings for a step of istep*dk
!          string length
           istr=1
           nkstr=1
           ikpt1=1
           do ikpt=1,dtefield%fnkpt
             do jstep = 1,istep
               ikpt1 = dtefield%ikpt_dk(ikpt1,1,idir)
             end do
             if (ikpt1 == 1) exit
             nkstr = nkstr + 1
           end do
!          Check that the string length is a divisor of nkpt
           if(mod(dtefield%fnkpt,nkstr) /= 0) then
             write(message,'(a,i5,a,i5,a,i7)')&
&             '  For istep = ', istep,&
&             '  The string length = ',nkstr,&
&             ', is not a divisor of fnkpt =',dtefield%fnkpt
             MSG_BUG(message)
           end if
           nstr = dtefield%fnkpt/nkstr

           write(message,'(a,i1,a,i2,a,i3,a,i6)')&
&           '  berryphase_new: for direction ',idir, ' and istep ', istep, ', nkstr = ',nkstr,&
&           ', nstr = ',nstr
           call wrtout(std_out,message,'COLL')
           call wrtout(ab_out,message,'COLL')

           ABI_ALLOCATE(idxkstr_mult,(nkstr,nstr))
           iunmark = 1
           kpt_mark(:)=0
           do istr=1,nstr
             do while(kpt_mark(iunmark) /= 0)
               iunmark = iunmark + 1
             end do
             idxkstr_mult(1,istr) = iunmark
             kpt_mark(iunmark)=1

             ikpt1 = idxkstr_mult(1,istr)
             do jkstr=2, nkstr
               do jstep = 1, istep
                 ikpt1 = dtefield%ikpt_dk(ikpt1,1,idir)
               end do
               idxkstr_mult(jkstr,istr) = ikpt1
               kpt_mark(ikpt1) = 1
             end do
           end do
         else
           nstr = dtefield%nstr(idir)
           nkstr = dtefield%nkstr(idir)
           ABI_ALLOCATE(idxkstr_mult,(nkstr,nstr))
           idxkstr_mult(:,:) = dtefield%idxkstr(1:nkstr,1:nstr,idir)
         end if
!        DEBUG
!        do istr=1,nstr
!        write(std_out,*)'string ', idxkstr_mult(:,istr)
!        end do
!        ENDBEBUG

         ABI_ALLOCATE(det_string,(2,nstr))
         ABI_ALLOCATE(polberry,(nstr))

         polbtot = zero
         do isppol = 1, nsppol

           det_string(1,:) = one ; det_string(2,:) = zero
           dtm_k(:) = one
           det_average(:) = zero


           do istr = 1, nstr

             if(calc_epaw3_force) epawf3_str(:,:,:) = zero
             if(calc_epaw3_stress) epaws3_str(:,:,:) = zero

             do jkstr = 1, nkstr

               ikpt=idxkstr_mult(jkstr,istr)

               dtm_real=dtm_mult(1,ikpt+(isppol-1)*dtefield%fnkpt,istep)
               dtm_imag=dtm_mult(2,ikpt+(isppol-1)*dtefield%fnkpt,istep)

               dtm_k(1) = det_string(1,istr)*dtm_real - &
&               det_string(2,istr)*dtm_imag
               dtm_k(2) = det_string(1,istr)*dtm_imag + &
&               det_string(2,istr)*dtm_real
               det_string(1:2,istr) = dtm_k(1:2)
!              DEBUG
!              write(std_out,'(a,i4,3e15.4,2e15.4)')'ikpt, kpt, dtm', ikpt, dtefield%fkptns(:,ikpt),  dtm_k
!              ENDDEBUG

               if(calc_epaw3_force) then
                 do iatom = 1, natom
                   do fdir = 1, 3
                     epawf3_str(iatom,idir,fdir) = epawf3_str(iatom,idir,fdir) + &
&                     dsdr_sum(iatom,fdir,ikpt+(isppol-1)*dtefield%fnkpt)
                   end do ! end loop over fdir
                 end do ! end loop over natom
               end if ! end check on calc_epaw3_force
               if(calc_epaw3_stress) then
                 do iatom = 1, natom
                   do fdir = 1, 6
                     epaws3_str(iatom,idir,fdir) = epaws3_str(iatom,idir,fdir) + &
&                     dsds_sum(iatom,fdir,ikpt+(isppol-1)*dtefield%fnkpt)
                   end do ! end loop over fdir
                 end do ! end loop over natom
               end if ! end check on calc_epaw3_stress

             end do

             if(calc_epaw3_force) then
               do iatom = 1, natom
                 do fdir = 1, 3
                   dtefield%epawf3(iatom,idir,fdir) = dtefield%epawf3(iatom,idir,fdir) + &
&                   epawf3_str(iatom,idir,fdir)
                 end do ! end loop over fdir
               end do ! end loop over natom
             end if ! end check on calc_epaw3_force
             if(calc_epaw3_stress) then
               do iatom = 1, natom
                 do fdir = 1, 6
                   dtefield%epaws3(iatom,idir,fdir) = dtefield%epaws3(iatom,idir,fdir) + &
&                   epaws3_str(iatom,idir,fdir)
                 end do ! end loop over fdir
               end do ! end loop over natom
             end if ! end check on calc_epaw3_stress

             det_average(:) = det_average(:) + &
&             det_string(:,istr)/dble(nstr)

           end do


!          correction to obtain a smouth logarithm of the determinant
           ABI_ALLOCATE(str_flag,(nstr))
!          DEBUG
!          since we don't have any case of non-nul Chern number,
!          we must change the det_string value "by brute force" if we want debug this
!          allocate(det_string_test(2,dtefield%nstr(idir)))
!          det_string_test(:,:)=det_string(:,:)
!          kk=0
!          det_string(1,1)=cos(2._dp*Pi*real(kk,dp)/four)
!          det_string(2,1)=sin(2._dp*Pi*real(kk,dp)/four)
!          jj=dtefield%str_neigh(1,1,idir)
!          ll=dtefield%str_neigh(2,1,idir)
!          do while (jj/=1)
!          kk=kk+1
!          det_string(1,jj)=cos(2._dp*Pi*real(kk,dp)/four)
!          det_string(2,jj)=sin(2._dp*Pi*real(kk,dp)/four)
!          det_string(1,ll)=cos(-2._dp*Pi*real(kk,dp)/four)
!          det_string(2,ll)=sin(-2._dp*Pi*real(kk,dp)/four)
!          jj=dtefield%str_neigh(1,jj,idir)
!          ll=dtefield%str_neigh(2,ll,idir)
!          enddo
!          ENDDEBUG
           if (istep==1) then
             do ineigh_str = 1,2
               str_flag(:)=0
               delta_str(:) = &
&               dtefield%coord_str(:,dtefield%str_neigh(ineigh_str,1,idir),idir) - dtefield%coord_str(:,1,idir)
               dstr(:)= delta_str(:) - nint(delta_str(:)) - real(dtefield%strg_neigh(ineigh_str,1,:,idir),dp)
               dist_=0._dp
               do kk = 1,2
                 do jj = 1,2
                   dist_ = dist_ + dstr(kk)*dtefield%gmet_str(kk,jj,idir)*dstr(jj)
                 end do
               end do
               dist_=sqrt(dist_)
               do istr = 1,dtefield%nstr(idir)
                 if(str_flag(istr)==0)then
!                  write(std_out,*)'new string'
                   str_flag(istr)=1
                   call rhophi(det_string(:,istr),dphase,rho)
!                  write(std_out,'(i4,e15.4,e15.4,e15.4)')istr, det_string(:,istr),dphase
                   dphase_init=dphase
                   jstr = dtefield%str_neigh(ineigh_str,istr,idir)
                   do while (istr/=jstr)
                     str_flag(jstr)=1
                     call rhophi(det_string(:,jstr),dphase_new,rho)
                     jj=nint((dphase_new-dphase)/(2._dp*Pi))
!                    DEBUG
!                    write(std_out,'(i4,e15.4,e15.4,e15.4,e15.4,i4)')jstr, det_string(:,jstr),dphase_new,dphase_new-dphase,jj
!                    ENDDEBUG
                     dphase_new=dphase_new-two*Pi*real(jj,dp)
                     if(jj/=0)then
                       write(message,'(6a)') ch10,&
&                       ' berryphase_new : WARNING -',ch10,&
&                       '  the berry phase has some huge variation in the space of strings of k-points',ch10,&
&                       '  ABINIT is trying to correct the berry phase, but it is highly experimental'
                       call wrtout(std_out,message,'PERS')
                     end if
!                    if(jj/=0)write(std_out,'(i4,e15.4,e15.4,e15.4,e15.4)')jstr, det_string(:,jstr),dphase_new,dphase_new-dphase
                     dphase=dphase_new
                     jstr=dtefield%str_neigh(ineigh_str,jstr,idir)
                   end do
!                  write(std_out,*)dphase_init, dphase, (dphase-dphase_init)/(2._dp*Pi),nint((dphase-dphase_init)/(2._dp*Pi))
                 end if
               end do
             end do
           end if
           ABI_DEALLOCATE(str_flag)
!          DEBUG
!          deallocate(dist_str)
!          det_string(:,:)=det_string_test(:,:)
!          deallocate(det_string_test)
!          ENDDEBUG

!          First berry phase that corresponds to det_average
!          phase0 = atan2(det_average(2),det_average(1))
           call rhophi(det_average,phase0,rho)
           det_mod = det_average(1)**2+det_average(2)**2

!          Then berry phase that corresponds to each string relative to the average
           do istr = 1, nstr

             rel_string(1) = (det_string(1,istr)*det_average(1) + &
             det_string(2,istr)*det_average(2))/det_mod
             rel_string(2) = (det_string(2,istr)*det_average(1) - &
             det_string(1,istr)*det_average(2))/det_mod
!            dphase = atan2(rel_string(2),rel_string(1))
             call rhophi(rel_string,dphase,rho)
             polberry(istr) = dtefield%sdeg*(phase0 + dphase)/two_pi
             polb_mult(isppol,istep) = polb_mult(isppol,istep) + polberry(istr)/(istep*dtefield%nstr(idir))
             polb(isppol) = zero
             do jstep=1, istep
               polb(isppol)=polb(isppol)+coef(jstep,istep)*polb_mult(isppol,jstep)
             end do

           end do

           if(berrystep>1)then
             write(message,'(9x,a,7x,e16.9,1x,a,i4,a,i4,a)')&
&             'total',polb_mult(isppol,istep),'(isppol=',isppol,', istep=',istep,')'!,ch10
             call wrtout(std_out,message,'COLL')

             write(message,'(3x,a,7x,e16.9,1x,a,i4,a,i4,a,a)')&
&             '+correction',polb(isppol),'(isppol=',isppol,', istep=1..',istep,')',ch10
             call wrtout(std_out,message,'COLL')

           else

             write(message,'(9x,a,7x,e16.9,1x,a,i4,a)')&
&             'total',polb_mult(isppol,istep),'(isppol=',isppol,')'!,ch10
             call wrtout(std_out,message,'COLL')
           end if

           polbtot = polbtot + polb(isppol)

         end do    ! isppol

!        Fold into interval [-1,1]
         polbtot = polbtot - 2_dp*nint(polbtot/2_dp)

         ABI_DEALLOCATE(det_string)
         ABI_DEALLOCATE(polberry)

!        ==========================================================================

!        Compute the ionic Berry phase

         call xred2xcart(natom,rprimd,xcart,xred)
         politot = zero
         write(message,'(a)')' Compute the ionic contributions'
         call wrtout(std_out,message,'COLL')

         write(message,'(a,2x,a,2x,a,15x,a)')ch10,&
&         'itom', 'itypat', 'polion'
         call wrtout(std_out,message,'COLL')

         do iatom = 1, natom
           itypat = typat(iatom)

!          The ionic phase can be computed much easier
           polion = zion(itypat)*xred(idir,iatom)

!          Fold into interval (-1,1)
           polion = polion - 2_dp*nint(polion/2_dp)
           politot = politot + polion
           write(message,'(2x,i2,5x,i2,10x,e16.9)') iatom,itypat,polion
           call wrtout(std_out,message,'COLL')
         end do

!        Fold into interval [-1,1] again
         politot = politot - 2_dp*nint(politot/2_dp)
         pion(idir) = politot

         write(message,'(9x,a,7x,es19.9)') 'total',politot
         call wrtout(std_out,message,'COLL')


!        ==========================================================================

!        Compute the total polarization

         poltot = politot + polbtot

         if (berrystep==1)then
           write(message,'(a,a)')ch10,&
&           ' Summary of the results'
           call wrtout(std_out,message,'COLL')
           if (unit_out /= 0) then
             call wrtout(unit_out,message,'COLL')
           end if
         else
           write(message,'(a,a,i4)')ch10,&
&           ' Summary of the results for istep =',istep
           call wrtout(std_out,message,'COLL')
           if (unit_out /= 0) then
             call wrtout(unit_out,message,'COLL')
           end if
         end if

         write(message,'(a,es19.9)')&
&         ' Electronic Berry phase ' ,polbtot
         call wrtout(std_out,message,'COLL')
         if (unit_out /= 0) then
           call wrtout(unit_out,message,'COLL')
         end if

         write(message,'(a,es19.9)') &
&         '            Ionic phase ', politot
         call wrtout(std_out,message,'COLL')
         if (unit_out /= 0) then
           call wrtout(unit_out,message,'COLL')
         end if

         write(message,'(a,es19.9)') &
&         '            Total phase ', poltot
         call wrtout(std_out,message,'COLL')
         if (unit_out /= 0) then
           call wrtout(unit_out,message,'COLL')
         end if

!        REC start
         if(abs(dtset%polcen(idir))>tol8)then
           poltot = poltot-dtset%polcen(idir)
           write(message,'(a,f15.10)') &
&           '    Translating Polarization by P0 for centrosymmetric cell: ',&
&           dtset%polcen(idir)
           call wrtout(std_out,message,'COLL')
           if (unit_out /= 0) then
             call wrtout(unit_out,message,'COLL')
           end if
         end if
!        REC end

         poltot = poltot - 2.0_dp*nint(poltot/2._dp)
         write(message,'(a,es19.9)') &
&         '    Remapping in [-1,1] ', poltot
         call wrtout(std_out,message,'COLL')
         if (unit_out /= 0) then
           call wrtout(unit_out,message,'COLL')
         end if

!        ! REC and HONG
!        =====================================================================================
!        Polarization branch control  (start)
!        -------------------------------------------------------------------------------------
!        berrysav == 0,  for non fixed D/d calculation, polarizaion is in [-1,1],done above
!        for fixed D/d calculation, choose polarization to minimize internal 
!        energy, or minimize |red_efiled|. (red_dfield=red_efiled+red_ptot)  
!        (d=e+p, as (26) of Stengel, Suppl.)
!        This is default value.
!        
!        berrysav == 1,  keep the polarization on the same branch, which saved in file POLSAVE
!        ======================================================================================

!        for fixed D/d calculation, choose polarization to minimize internal energy, or to minimize reduced electric field |red_efield|
         if(dtset%berrysav ==0 .and. (dtset%berryopt == 6 .or. dtset%berryopt == 7 .or. & 
&         dtset%berryopt == 16 .or. dtset%berryopt == 17))  then   

           jump=-nint(dtset%red_dfield(idir) - poltot)   ! red_efield = red_dfield - poltot

           if(jump /= 0)then
             write(message,'(a,i1,a,es19.9,a,i2)') &
&             ' P(',idir,') Shifted polarization branch to minimize red_efield &
&             k from ',poltot, ' by ',jump
             call wrtout(std_out,message,'COLL')
             if (unit_out /= 0) then
               call wrtout(unit_out,message,'COLL')
             end if
             poltot=poltot-jump
           end if
           pol0(idir)=poltot

         end if


!        keep the polarization on the same branch. 
         if (dtset%berrysav == 1) then 

!          use saved polarization to keep on same branch
           inquire(file='POLSAVE',exist=lexist)
           if(lexist)then
             if(idir==1)then
               if(mpi_enreg%me==0)then
                 if (open_file('POLSAVE',message,newunit=unt,status='OLD') /= 0) then
                   MSG_ERROR(message)
                 end if
                 read(unt,*)pol0
                 write(message,'(a,3f20.12)')'Reading old polarization:',pol0
                 call wrtout(std_out,message,'COLL')
                 if (unit_out /= 0) then
                   call wrtout(unit_out,message,'COLL')
                 end if
                 close(unt)
               end if
               call xmpi_bcast(pol0,0,spaceComm,ierr)
             end if
           else   
             pol0(idir)=poltot
           end if
           jump=nint(poltot-pol0(idir))
           if(jump /= 0)then
             write(message,'(a,i1,a,es19.9,a,i2)') &
&             ' P(',idir,') jumped to new branch. Shifting bac&
&             k from ',poltot, ' by ',jump
             call wrtout(std_out,message,'COLL')
             if (unit_out /= 0) then
               call wrtout(unit_out,message,'COLL')
             end if
             poltot=poltot-jump
           end if

           pol0(idir)=poltot

         end if

!        =====================================================================================
!        Polarization branch control  (end)
!        =====================================================================================


!        Transform the phase into a polarization
         fac = 1._dp/(gmod*dtefield%nkstr(idir))
!        !REC         fac = fac/ucvol
!        !REC         pol = fac*poltot
         red_ptot(idir)=poltot !!REC
         pol = fac*red_ptot(idir)/ucvol !!REC
         ptot(idir)=red_ptot(idir)/ucvol !!REC
         write(message,'(a,a,es19.9,a,a,a,es19.9,a,a)')ch10,&
&         '           Polarization ', pol,' (a.u. of charge)/bohr^2',ch10,&
&         '           Polarization ', pol*(e_Cb)/(Bohr_Ang*1d-10)**2,&
&         ' C/m^2',ch10
         call wrtout(std_out,message,'COLL')
         if (unit_out /= 0) then
           call wrtout(unit_out,message,'COLL')
         end if


         ABI_DEALLOCATE(idxkstr_mult)

       end do !istep

       pel(idir) = polbtot

     end if   ! option == 1 or option == 3

!    Write the ddk WF to a file

     if ((option == 2).or.(option == 3)) then

       pertcase = idir + 3*natom
       response = 1
       call appdig(pertcase,dtfil%fnameabo_1wf,fiwf1o)
       ABI_ALLOCATE(resid,(mband*nkpt*nsppol))
       resid(:) = zero

       call outwf(cg1,dtset,psps,eig_dum,fiwf1o,hdr,kg,dtset%kptns,&
&       mband,mcg,mkmem,mpi_enreg,mpw,natom,dtset%nband,&
&       nkpt,npwarr,nsppol,&
&       occ_dum,resid,response,dtfil%unwff2,wfs,wvl)

       ABI_DEALLOCATE(resid)
     end if  ! option == 2 or option == 3

   end if   ! rfdir(idir) == 1

 end do    ! Close loop over idir

!!REC start
 if (dtset%berrysav == 1) then
   if(mpi_enreg%me==0)then
     if (open_file('POLSAVE',message,newunit=unt,status='UNKNOWN') /= 0) then
       MSG_ERROR(message)
     end if 
     write(unt,'(3F20.12)') pol0
     close(unt)
   end if
   first=.false.
 end if
!!REC end

!Compute polarization in cartesian coordinates
 if ((dtset%rfdir(1) == 1).and.(dtset%rfdir(2) == 1).and.&
& (dtset%rfdir(3) == 1)) then

   if(usepaw.ne.1) then
     pelev=zero
   else
     call pawpolev(my_natom,natom,ntypat,pawrhoij,pawtab,pelev,&
&     comm_atom=mpi_enreg%comm_atom)
!    note that in the PAW case, the pelev contribution is already
!    implicitly included in the electronic polarization, from the
!    discretized derivative operator. In the NCPP case no such
!    terms exist anyway. Actually in the PAW formulation
!    such terms are included to all orders, unlike in USPP where only
!    zeroth and first-order terms are. In USPP the first-order term 
!    is pelev. Here we compute pelev separately only for reporting 
!    purposes in polcart, it is not added into pel or used in the the
!    PAW finite field code in make_grad_berry.F90
!    13 June 2012 J Zwanziger
   end if
   call polcart(red_ptot,pel,pel_cart,pelev,pion,pion_cart,3,&
&   ptot_cart,rprimd,ucvol,unit_out,usepaw)

 end if

!deallocate(pwind_k, dtm)
 ABI_DEALLOCATE(pwnsfac_k)
 ABI_DEALLOCATE(sflag_k)
 ABI_DEALLOCATE(cg1_k)
 if (option > 1)  then
   ABI_DEALLOCATE(cg1)
   ABI_DEALLOCATE(eig_dum)
   ABI_DEALLOCATE(occ_dum)
 end if

 if (usepaw == 1) then
   ABI_DEALLOCATE(dimlmn)
   call pawcprj_free(cprj_k)
   call pawcprj_free(cprj_kb)
   call pawcprj_free(cprj_gat)
   ABI_DATATYPE_DEALLOCATE(cprj_k)
   ABI_DATATYPE_DEALLOCATE(cprj_kb)
   ABI_DATATYPE_DEALLOCATE(cprj_gat)
   if (dtset%kptopt /= 3) then
     call pawcprj_free(cprj_ikn)
     call pawcprj_free(cprj_fkn)
     ABI_DATATYPE_DEALLOCATE(cprj_ikn)
     ABI_DATATYPE_DEALLOCATE(cprj_fkn)
   end if
   if (calc_epaw3_force) then
     ABI_DEALLOCATE(dsdr_sum)
     ABI_DEALLOCATE(epawf3_str)
   end if
   if (calc_epaw3_stress) then
     ABI_DEALLOCATE(dsds_sum)
     ABI_DEALLOCATE(epaws3_str)
   end if

   if (nproc>1) then
     call pawcprj_free(cprj_buf)
     ABI_DATATYPE_DEALLOCATE(cprj_buf)
   end if 

 end if

 ABI_DEALLOCATE(ikpt3)
 ABI_DEALLOCATE(ikpt3i)
 ABI_DEALLOCATE(sflag_k_mult)
 ABI_DEALLOCATE(npw_k3)
 ABI_DEALLOCATE(pwind_k_mult)
 ABI_DEALLOCATE(itrs_mult)
 ABI_DEALLOCATE(coef)
 ABI_DEALLOCATE(polb_mult)
 ABI_DEALLOCATE(dtm_mult)

!DEBUG
!write(std_out,*)'berryphase_new exit'
!END_DEBUG
end subroutine berryphase_new
!!***
