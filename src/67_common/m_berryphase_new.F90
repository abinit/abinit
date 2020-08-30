!!****m* ABINIT/m_berryphase_new
!! NAME
!!  m_berryphase_new
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2003-2020 ABINIT  group (MVeithen)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_berryphase_new

 use defs_basis
 use defs_wvltypes
 use m_efield
 use m_errors
 use m_abicore
 use m_xmpi
 use m_hdr
 use m_dtset
 use m_dtfil

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes,  only : MPI_type
 use m_berrytk,      only : smatrix, polcart
 use m_cgprj,        only : ctocprj
 use m_fftcore,      only : kpgsph
 use m_geometry,     only : xred2xcart, metric
 use m_io_tools,     only : open_file
 use m_iowf,         only : outwf, outresid
 use m_kg,           only : getph
 use m_kpts,         only : listkk, smpbz
 use m_mpinfo,       only : proc_distrb_cycle
 use m_numeric_tools,only : rhophi
 use m_pawang,       only : pawang_type
 use m_pawcprj,      only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_mpi_allgather, &
                            pawcprj_put, pawcprj_copy, pawcprj_mpi_recv,  &
                            pawcprj_mpi_send, pawcprj_free, pawcprj_getdim, pawcprj_symkn
 use m_paw_dfpt,     only : dsdr_k_paw
 use m_paw_efield,   only : pawpolev
 use m_pawrad, only : pawrad_type
 use m_pawrhoij,     only : pawrhoij_type
 use m_paw_sphharm, only : setsym_ylm
 use m_pawtab,       only : pawtab_type
 use m_paw_overlap,  only : expibi,qijb_kk,smatrix_k_paw
 use m_symtk,   only : symatm
 use m_time,    only : timab

 implicit none

 private
!!***

 public :: berryphase_new
 public :: prtefield
 public :: init_e_field_vars
 public :: initberry
 public :: update_e_field_vars
!!***

contains
!!***

!!****f* ABINIT/berryphase_new
!! NAME
!! berryphase_new
!!
!! FUNCTION
!! This routine computes the Berry Phase polarization
!!  and the finite difference expression of the ddk.
!!  See for example Na Sai et al., PRB 66, 104108 (2002) [[cite:Sai2002]]
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
!! mpi_enreg=information about MPI parallelization
!! mpw=maximum dimensioned size of npw
!! my_natom=number of atoms treated by current processor
!! natom=number of atoms in cell
!! nkpt=number of k points
!! npwarr(nkpt)=number of planewaves in basis at this k point
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! ntypat=number of types of atoms in unit cell
!! nkpt=number of k-points
!! calc_pol_ddk = 1: compute Berryphase polarization
!!                2: compute finite difference expression of the ddk
!!                3: compute polarization & ddk
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

subroutine berryphase_new(atindx1,cg,cprj,dtefield,dtfil,dtset,psps,&
&  gprimd,hdr,indlmn,kg,lmnmax,mband,mcg,mcprj,&
&  mkmem,mpi_enreg,mpw,my_natom,natom,npwarr,nsppol,ntypat,&
&  nkpt,calc_pol_ddk,pawrhoij,pawtab,pel,pelev,pion,ptot,red_ptot,pwind,&  !!REC
&  pwind_alloc,pwnsfac,&
&  rprimd,typat,ucvol,unit_out,usecprj,usepaw,xred,zion)

!Arguments ------------------------------------
 integer, intent(in) :: lmnmax,mband,mcg,mcprj,mkmem,mpw,my_natom,natom,nkpt
 integer, intent(in) :: nsppol,ntypat,calc_pol_ddk
 integer, intent(in) :: pwind_alloc,unit_out,usecprj,usepaw
 real(dp), intent(in) :: ucvol
 type(MPI_type), intent(in) :: mpi_enreg
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
 integer :: count,count1,dest,fdir,unt
 integer :: iatom,iband,icg,icg1,idir,idum,ikpt1i_sp
 integer :: ierr,ifor,ikg,ikpt,ikpt1,ikpt_loc!,ikpt2,ikpt2i,npw_k2, itrs
 integer :: icp1, icp2,icpgr_offset,iproc
! integer :: ii ! appears commented out below in a debug section
 integer :: inibz,ikpt1i
 integer :: isppol,istr,itypat,jband,jkpt,jkstr,jsppol
 integer :: det_inv_smat, det_smat, inv_smat
 integer :: maxbd,mcg1_k
 integer :: minbd,my_nspinor,nband_k,ncpgr,nfor,npw_k1,ntotcp,n2dim,nproc,pertcase
 integer :: response,shiftbd,source,spaceComm,tag
 integer :: jj,jstr,kk,ineigh_str
 integer :: istep,jstep,kpt_mark(dtefield%fnkpt),nkstr,nstr,iunmark,berrystep
 integer :: jkpt2, jkpt2i
 real(dp) :: det_mod,dkinv,dphase,dtm_real,dtm_imag,fac,gmod,phase0
 real(dp) :: pol,polbtot,polion,politot,poltot,rho
 logical :: calc_epaw3_force,calc_epaw3_stress,efield_flag
 integer :: polflag, ddkflag

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

 polflag = 1
 ddkflag = 1
 if (calc_pol_ddk == 1) then
   ddkflag = 0
 else if (calc_pol_ddk == 2) then
   polflag = 0
 end if

!allocate(pwind_k(mpw))
 ABI_ALLOCATE(pwnsfac_k,(4,mpw))
 ABI_ALLOCATE(sflag_k,(dtefield%mband_occ))
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
 if (ddkflag==1) then
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
 ABI_ALLOCATE(sflag_k_mult,(dtefield%mband_occ,berrystep))
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
         call pawcprj_mpi_allgather(cprj_k,cprj_gat,natom,n2dim,1,dimlmn,ncpgr,nproc,spaceComm,ierr,rank_ordered=.true.)
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
!nband_k = dtefield%mband_occ
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

 minbd = 1   ;  maxbd = dtefield%mband_occ

 if(calc_epaw3_force) dtefield%epawf3(:,:,:) = zero
 if(calc_epaw3_stress) dtefield%epaws3(:,:,:) = zero

 do idir = 1, 3

!  dtm(:,:) = zero
   dtm_mult(:,:,:) = zero
   if (calc_epaw3_force) dsdr_sum(:,:,:) = zero
   if (calc_epaw3_stress) dsds_sum(:,:,:) = zero

   if (dtset%rfdir(idir) /= 1) cycle

   if (abs(dtefield%efield_dot(idir)) < tol12) dtefield%sflag(:,:,:,idir) = 0

! calculate vector steps in k space
   dk(:) = dtefield%dkvecs(:,idir)
   gpard(:) = dk(1)*gprimd(:,1) + dk(2)*gprimd(:,2) + dk(3)*gprimd(:,3)
   gmod = sqrt(dot_product(gpard,gpard))

   write(message,'(a,a,a,3f9.5,a,a,3f9.5,a)')ch10,&
&   ' Computing the polarization (Berry phase) for reciprocal vector:',ch10,&
&   dk(:),' (in reduced coordinates)',ch10,&
&   gpard(1:3),' (in cartesian coordinates - atomic units)'
   call wrtout(std_out,message,'COLL')
   if (unit_out /= 0) then
     call wrtout(unit_out,message,'COLL')
   end if

   write(message,'(a,i5,a,a,i5)')&
&   ' Number of strings: ',dtefield%nstr(idir),ch10,&
&   ' Number of k points in string:', dtefield%nkstr(idir)
   call wrtout(std_out,message,'COLL')
   if (unit_out /= 0) then
     call wrtout(unit_out,message,'COLL')
   end if

!  Check whether the polarization or the ddk must be computed

!  nfor = 1 : to compute P, I only need the WF at k + dk
!  nfor = 2 : to compute the ddk, I need the WF at k + dk and k - dk
!  dkinv    : +-1/2dk


!  default for polarization
   nfor = 1
   if (ddkflag == 1) then
     nfor = 2
   end if

   if (ddkflag == 1) then

     cg1(:,:) = zero
     dkinv = one/(two*dk(idir))

     write(message,'(a,a,a,3f9.5,a,a,3f9.5,a)')ch10,&
&     ' Computing the ddk (Berry phase) for reciprocal vector:',ch10,&
&     dk(:),' (in reduced coordinates)',ch10,&
&     gpard(1:3),' (in cartesian coordinates - atomic units)'
     call wrtout(std_out,message,'COLL')
     if (unit_out /= 0) then
       call wrtout(unit_out,message,'COLL')
     end if

   end if

! From smatrix routine: det_inv_smat = type of calculation
!        1 : compute inverse of the overlap matrix
!       10 : compute determinant of the overlap matrix
!       11 : compute determinant and inverse of the overlap matrix
   inv_smat = 0
   det_smat = 0

!  for ddk need inverse matrix
   if (ddkflag == 1) then
     inv_smat = 1
   end if

!   if polarization is requested need smat determinant as well
   if (polflag == 1) then
     det_smat = 1
   end if

! electric fields with PAW also needs S_inverse for forces and stresses, even just for polarization
   if (calc_epaw3_force .or. calc_epaw3_stress) then
     inv_smat = 1
   end if

   det_inv_smat = 10*det_smat  + inv_smat

!--------------------------------------------------------------------
!  for each dk we require, calculate the smatrix, derivatives etc...
!--------------------------------------------------------------------
   do ifor = 1, nfor

     if (ifor == 2) then
       dk(:) = -1_dp*dk(:)
!      only the inverse of the overlap matrix is required on second pass, speeds things up a bit
       det_inv_smat = 1
       dkinv = -1_dp*dkinv
     end if


!    Compute the determinant and/or the inverse of the overlap matrix
!    for each pair of k-points < u_nk | u_nk+dk >

     icg = 0 ; icg1 = 0
     ABI_ALLOCATE(smat_k,(2,dtefield%mband_occ,dtefield%mband_occ))
     ABI_ALLOCATE(smat_inv,(2,dtefield%mband_occ,dtefield%mband_occ))
     ABI_ALLOCATE(smat_k_paw,(2,usepaw*dtefield%mband_occ,usepaw*dtefield%mband_occ))
     if (calc_epaw3_force .or. calc_epaw3_stress) then ! dsdr needed for forces and stresses in electric field with PAW
       ABI_ALLOCATE(dsdr,(2,natom,ncpgr,usepaw*dtefield%mband_occ,usepaw*dtefield%mband_occ))
       dsdr = zero
     end if


!    Loop on the values of ikpt_loc and ikpt1 :
!    ikpt1 is incremented one by one, and number the k points in the FBZ
!    ikpt1i refer to the k point numbering in the IBZ
!    ikpt_loc differs from ikpt1 only in the parallel case, and gives
!    the index of the k point in the FBZ, in the set treated by the present processor
!    NOTE : in order to allow synchronisation, ikpt_loc contain information about
!    ikpt AND ISPPOL !
!    It means that the following loop is equivalent to a double loop :
!    do isppol = 1, nsppol
!    do ikpt1 =  1, dtefield%fmkmem
!
     do ikpt_loc = 1, dtefield%fmkmem_max*nsppol

       ikpt1=mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc,1)
       isppol=mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc,2)

!      if this k and spin are for me do it
       if (ikpt1 > 0 .and. isppol > 0) then

         ikpt1i = dtefield%indkk_f2ibz(ikpt1,1)
         nband_k = dtset%nband(ikpt1i + (isppol-1)*dtset%nkpt)

!        DEBUG
!        Please keep this debugging feature
!        write(std_out,'(a,5i4)' )' berryphase_new : ikpt_loc,ikpt1,isppol,idir,ifor=',&
!        &                                  ikpt_loc,ikpt1,isppol,idir,ifor
!        ENDDEBUG

         inibz=0
         if (dtset%kptns(1,ikpt1i) == dtefield%fkptns(1,ikpt1) .and. &
&         dtset%kptns(2,ikpt1i) == dtefield%fkptns(2,ikpt1) .and. &
&         dtset%kptns(3,ikpt1i) == dtefield%fkptns(3,ikpt1)) inibz=1

         ikg = dtefield%fkgindex(ikpt1)
!        ikpt2 = dtefield%ikpt_dk(ikpt1,ifor,idir)
!        ikpt2i = dtefield%indkk_f2ibz(ikpt2,1)

!        ikpt3(istep) : index of kpt1 + istep*dk in the FBZ
!        ikpt3i(istep) : index of kpt1 + istep*dk in the IBZ
         ikpt3(1) = dtefield%ikpt_dk(ikpt1,ifor,idir)
         ikpt3i(1) = dtefield%indkk_f2ibz(ikpt3(1),1)
         do istep = 1, berrystep-1
           ikpt3(istep+1) = dtefield%ikpt_dk(ikpt3(istep),ifor,idir)
           ikpt3i(istep+1) = dtefield%indkk_f2ibz(ikpt3(istep+1),1)
         end do

!        itrs = 0
!        if (dtefield%indkk_f2ibz(ikpt1,6) == 1 ) itrs = itrs + 1
!        if (dtefield%indkk_f2ibz(ikpt2,6) == 1 ) itrs = itrs + 10

         itrs_mult(:)=0
         if (dtefield%indkk_f2ibz(ikpt1,6) == 1 ) itrs_mult(:) = itrs_mult(:) + 1
         do istep=1,berrystep
           if (dtefield%indkk_f2ibz(ikpt3(istep),6) == 1 ) itrs_mult(istep) = itrs_mult(istep) + 10
         end do

         npw_k1 = npwarr(ikpt1i)
!        npw_k2 = npwarr(ikpt2i)

         do istep = 1, berrystep
           npw_k3(istep)=npwarr(ikpt3i(istep))
         end do

!        ji: the loop is over the FBZ, but sflag and smat only apply to the IBZ
         if ( efield_flag .and. inibz == 1) then  !!HONG
           ikpt1i_sp=ikpt1i+(isppol-1)*dtset%nkpt
           smat_k(:,:,:) = dtefield%smat(:,:,:,ikpt1i_sp,ifor,idir)
         else
           smat_k(:,:,:) = zero
         end if

!        pwind_k(1:npw_k1) = pwind(ikg+1:ikg+npw_k1,ifor,idir)
         pwnsfac_k(1,1:npw_k1) = pwnsfac(1,ikg+1:ikg+npw_k1)
         pwnsfac_k(2,1:npw_k1) = pwnsfac(2,ikg+1:ikg+npw_k1)

!        the array needed to compute the overlap matrix between k and k+istep*dk (with multiple steps)
!        the 0-case (no corresponding pw in k and k+dk) could be handled better (k+2*dk could have a corresponding pw ?)
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

!        DEBUG
!        write(std_out,*)' berryphase_new : dtset%berryopt,inibz,ikpt1i,isppol,dtset%nkpt,ifor,idir', &
!        &          dtset%berryopt,inibz,ikpt1i,isppol,dtset%nkpt,ifor,idir
!        write(std_out,'(a,4i4)' )' berryphase_new : sflag_k(:)=',sflag_k(:)
!        ENDDEBUG

         if ( efield_flag .and. inibz == 1) then  !!HONG
           ikpt1i_sp=ikpt1i+(isppol-1)*dtset%nkpt
           sflag_k(:) = dtefield%sflag(:,ikpt1i_sp,ifor,idir)
         else
           sflag_k(:) = 0
         end if

         if (usepaw == 1) then
           icp1=dtefield%cprjindex(ikpt1i,isppol)
           call pawcprj_get(atindx1,cprj_k,cprj,natom,1,icp1,ikpt1i,0,isppol,&
&           mband,mkmem,natom,dtefield%mband_occ,dtefield%mband_occ,&
&           my_nspinor,nsppol,0,mpicomm=mpi_enreg%comm_kpt,&
&           proc_distrb=mpi_enreg%proc_distrb)

           if ( ikpt1i /= ikpt1 ) then
             call pawcprj_copy(cprj_k,cprj_ikn)
             call pawcprj_symkn(cprj_fkn,cprj_ikn,dtefield%atom_indsym,dimlmn,-1,indlmn,&
&             dtefield%indkk_f2ibz(ikpt1,2),dtefield%indkk_f2ibz(ikpt1,6),&
&             dtefield%fkptns(:,dtefield%i2fbz(ikpt1i)),&
&             dtefield%lmax,dtefield%lmnmax,mband,natom,dtefield%mband_occ,my_nspinor,&
&             dtefield%nsym,ntypat,typat,dtefield%zarot)
             call pawcprj_copy(cprj_fkn,cprj_k)
           end if

         end if ! end if usepaw

!        DEBUG
!        write(std_out,'(a,4i4)' )' berryphase_new : sflag_k(:)=',sflag_k(:)
!        ENDDEBUG

!        DEBUG
!        write(std_out,'(a,7i4)')'me, idir,ifor, ikpt_loc, ikpt1, isppol = ',&
!        & me,idir,ifor,ikpt_loc,ikpt1,isppol
!        write(std_out,'(a,10i3)')'pwind_k(1:10) = ',pwind_k(1:10)
!        ENDDEBUG

         do istep=1,berrystep
           sflag_k_mult(:,istep) = sflag_k(:)
         end do

       end if ! end check that ikpt1 > 0 and isppol > 0

!      --------------------------------------------------------------------------------
!      Communication
!      --------------------------------------------------------------------------------

       do istep=1,berrystep

!        if(ikpt_loc <= nsppol*dtefield%fmkmem) then
         if (ikpt1 > 0 .and. isppol > 0) then ! I currently have a true kpt to use

           count = npw_k3(istep)*my_nspinor*nband_k
           ABI_ALLOCATE(cgq,(2,count))
           cgq = zero
           source = me
           if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt3i(istep),1,nband_k,isppol,me)) then
!            I need the datas from someone else
             source = mpi_enreg%proc_distrb(ikpt3i(istep),1,isppol)
           end if
         else
           source = -1 ! I do not have a kpt to use
         end if

         do dest = 0, nproc-1

           if ((dest==me) .and. (ikpt1>0) .and. (isppol>0)) then
!            I am destination and I have something to do
!            if (mpi_enreg%paral_compil_kpt == 1) write(std_out,*) &
!            &               'coucou 2, mpi_enreg%proc_distrb(ikpt3i(istep),1:nband_k,isppol) : ', &
!            &               mpi_enreg%proc_distrb(ikpt3i(istep),1:nband_k,isppol)
!            write(std_out,*)'ikpt3i(istep) ', ikpt3i(istep)
!            write(std_out,*)'nband_k ',nband_k
!            write(std_out,*)'isppol ', isppol
!            write(std_out,*)'mpi_enreg%proc_distrb',mpi_enreg%proc_distrb

             if (source == me) then
!              I am destination and source
!              DEBUG
!              write(std_out,*)'copying ... '
!              write(std_out,*)'me: ',me, 'ikpt3i(istep) ', ikpt3i(istep), 'isppol ', isppol
!              ENDDEBUG

!              pwnsfac
               idum = dtefield%fkgindex(ikpt3(istep))
               pwnsfac_k(3,1:npw_k3(istep)) = pwnsfac(1,idum+1:idum+npw_k3(istep))
               pwnsfac_k(4,1:npw_k3(istep)) = pwnsfac(2,idum+1:idum+npw_k3(istep))

!              cgq (and cprj)
               icg1 = dtefield%cgindex(ikpt3i(istep),isppol)

               if (usepaw == 1) then
                 icp2=dtefield%cprjindex(ikpt3i(istep),isppol)
                 call pawcprj_get(atindx1,cprj_kb,cprj,natom,1,icp2,ikpt3i(istep),0,isppol,&
&                 mband,mkmem,natom,dtefield%mband_occ,dtefield%mband_occ,my_nspinor,&
&                 nsppol,0,mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
               end if

               cgq(:,1:count)  = cg(:,icg1+1:icg1+count)
!              if (usepaw == 1) then
!                call pawcprj_copy(cprj_buf,cprj_kb)
!              end if

!              if ((source /= me)) then
             else
!              I am the destination but not the source -> receive
!              DEBUG
!              write(std_out,'(a)')'receiving ...'
!              write(std_out,'(a,i4,a,i4,a,i4,a,i4)')'me: ',me, 'source ', source,'ikpt3i(istep) ', ikpt3i(istep), 'isppol ', isppol
!              ENDDEBUG

!              receive pwnsfac
               ABI_ALLOCATE(buffer,(2,npw_k3(istep)))
               tag = ikpt3(istep) + (isppol - 1)*dtefield%fnkpt
               call xmpi_recv(buffer,source,tag,spaceComm,ierr)
               pwnsfac_k(3,1:npw_k3(istep)) = buffer(1,1:npw_k3(istep))
               pwnsfac_k(4,1:npw_k3(istep)) = buffer(2,1:npw_k3(istep))
               ABI_DEALLOCATE(buffer)

!              receive cgq (and cprj)
               tag = ikpt3i(istep) + (isppol - 1)*nkpt
               call xmpi_recv(cgq,source,tag,spaceComm,ierr)

               if (usepaw == 1) then
                 call pawcprj_mpi_recv(natom,n2dim,dimlmn,ncpgr,cprj_kb,source,spaceComm,ierr)
               end if

             end if

           else if (dest /= me) then

!            jkpt is the kpt which is being treated by dest
!            jsppol is his isppol
             jkpt = mpi_enreg%kpt_loc2fbz_sp(dest, ikpt_loc,1)
             jsppol = mpi_enreg%kpt_loc2fbz_sp(dest, ikpt_loc,2)

             if (jkpt > 0 .and. jsppol > 0) then ! dest is treating a true kpt

               jkpt2 = dtefield%ikpt_dk(jkpt,ifor,idir)
               jkpt2i = dtefield%indkk_f2ibz(jkpt2,1)

!              check if I am his source
               if((mpi_enreg%proc_distrb(jkpt2i,1,jsppol) == me))  then
!                I know something about jkpt3i and I must send it
!                DEBUG
!                write(std_out,'(a)')'sending ...'
!                write(std_out,'(a,i4,a,i4,a,i4,a,i4)')'dest: ',dest,' me: ',me,&
!                &                          ' jkpt2i ',jkpt2i,' jsppol: ',jsppol
!                ENDDEBUG

!                pwnsfac
                 tag = jkpt2 + (jsppol - 1)*dtefield%fnkpt
                 count1 = npwarr(jkpt2i)
                 ABI_ALLOCATE(buffer,(2,count1))
                 idum = dtefield%fkgindex(jkpt2)
                 buffer(1,1:count1)  = pwnsfac(1,idum+1:idum+count1)
                 buffer(2,1:count1)  = pwnsfac(2,idum+1:idum+count1)
                 call xmpi_send(buffer,dest,tag,spaceComm,ierr)
                 ABI_DEALLOCATE(buffer)

!                cgq (and cprj)
                 icg1 = dtefield%cgindex(jkpt2i,jsppol)

                 if (usepaw == 1) then
                   icp2=dtefield%cprjindex(jkpt2i,jsppol)
                   call pawcprj_get(atindx1,cprj_buf,cprj,natom,1,icp2,jkpt2i,0,jsppol,&
&                   mband,mkmem,natom,dtefield%mband_occ,dtefield%mband_occ,&
&                   my_nspinor,nsppol,0,mpicomm=mpi_enreg%comm_kpt,&
&                   proc_distrb=mpi_enreg%proc_distrb)
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
&               dtefield%indkk_f2ibz(ikpt3(istep),2),dtefield%indkk_f2ibz(ikpt3(istep),6),&
&               dtefield%fkptns(:,dtefield%i2fbz(ikpt3i(istep))),&
&               dtefield%lmax,dtefield%lmnmax,mband,natom,&
&               dtefield%mband_occ,my_nspinor,dtefield%nsym,ntypat,typat,&
&               dtefield%zarot)
               call pawcprj_copy(cprj_fkn,cprj_kb)
             end if
             call smatrix_k_paw(cprj_k,cprj_kb,dtefield,idir,ifor,mband,natom,smat_k_paw,typat)
!            write(std_out,'(a,5i4)')' JWZ berryphase_new : ikpt_loc,ikpt1,ikpt1i,ikpt2,ikpt2i ',ikpt_loc,ikpt1,ikpt1i,ikpt3(istep),ikpt3i(istep)
!            call smatrix_k0_paw(atindx1,cprj_k,cprj_k,dtefield,ikpt1i,idir,ifor,&
!            &                                  mband,mpi_enreg,natom,ntypat,pawtab,smat_k_paw,typat)
             if (calc_epaw3_force .or. calc_epaw3_stress) then
               call dsdr_k_paw(cprj_k,cprj_kb,dsdr,dtefield,idir,ifor,mband,natom,ncpgr,typat)
             end if
           end if

           icg1 = 0
           icg = dtefield%cgindex(ikpt1i,isppol)
!          DEBUG
!          if(istep<=2)then
!          if(ikpt1==1)then
!          write(std_out,'(a,2i4,3e15.4)')'istep ikpt3, kpt, cgq', istep, ikpt3(istep), dtefield%fkptns(:,ikpt3(istep))
!          write(std_out,*) cgq
!          write(std_out,*)
!          end if
!          end if
!          ENDDEBUG
           call smatrix(cg,cgq,cg1_k,ddkflag,dtm_k,icg,icg1,itrs_mult(istep),det_inv_smat,maxbd,&
&           mcg,count,mcg1_k,minbd,&
&           mpw,dtefield%mband_occ,dtefield%nband_occ(isppol),&
&           npw_k1,npw_k3(istep),my_nspinor,pwind_k_mult(:,istep),pwnsfac_k,sflag_k_mult(:,istep),&
&           shiftbd,smat_inv,smat_k,smat_k_paw,usepaw)

! in finite electric field case with paw must save additional F3 term in forces
           if(calc_epaw3_force) then
! when ncpgr = 3, gradients are wrt to atom displacements
! but when ncpgr = 9, first 6 gradients are wrt strains, last three are displacements
             icpgr_offset = 0
             if (ncpgr == 9) icpgr_offset = 6
             do iatom = 1, natom
               do fdir = 1, 3
                 dsdr_sum(iatom,fdir,ikpt1+(isppol-1)*dtefield%fnkpt) = zero
                 do iband = 1, dtefield%nband_occ(isppol)
                   do jband = 1, dtefield%nband_occ(isppol)
! collect Im{Trace{S^{-1}.dS/dR}} for this k point
                     dsdr_sum(iatom,fdir,ikpt1+(isppol-1)*dtefield%fnkpt) = &
&                     dsdr_sum(iatom,fdir,ikpt1+(isppol-1)*dtefield%fnkpt) + &
&                     smat_inv(2,iband,jband)*dsdr(1,iatom,icpgr_offset+fdir,jband,iband) + &
&                     smat_inv(1,iband,jband)*dsdr(2,iatom,icpgr_offset+fdir,jband,iband)
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
                 do iband = 1, dtefield%nband_occ(isppol)
                   do jband = 1, dtefield%nband_occ(isppol)
! collect Im{Trace{S^{-1}.dS/de}} for this k point
                     dsds_sum(iatom,fdir,ikpt1+(isppol-1)*dtefield%fnkpt) = &
&                     dsds_sum(iatom,fdir,ikpt1+(isppol-1)*dtefield%fnkpt) + &
&                     smat_inv(2,iband,jband)*dsdr(1,iatom,fdir,jband,iband) + &
&                     smat_inv(1,iband,jband)*dsdr(2,iatom,fdir,jband,iband)
                   end do ! end sum over jband
                 end do ! end sum over iband
               end do ! end sum over fdir
             end do ! end sum over iatom
           end if ! end check on calc_epaw3_stress

           if ((det_inv_smat == 10).or.(det_inv_smat == 11)) then

             if (sqrt(dtm_k(1)*dtm_k(1) + dtm_k(2)*dtm_k(2)) < tol12) then
               ! EB: the MSG_BUG has been replaced here by what is done in 67_common/m_cgwf.F90
               ! This avoid the code to stop for phonons under E-field too.
               ! TODO: Since the same is done in m_cgwf.F90 and in m_berryphase_new.F90,
               ! rationalization should be done with one single module.
               write(message,'(a,i5,a,a,a)')&
&               '  For k-point #',ikpt1,',',ch10,&
&               '  the determinant of the overlap matrix is found to be 0. Fixing...'
                ! Try this:
                write(std_out,*)message,dtm_k(1:2)
                if(abs(dtm_k(1))<=1d-12)dtm_k(1)=1d-12
                if(abs(dtm_k(2))<=1d-12)dtm_k(2)=1d-12
                write(std_out,*)' Changing to:',dtm_k(1:2)
!               MSG_BUG(message)
             end if

             dtm_mult(1,ikpt1+(isppol-1)*dtefield%fnkpt,istep) = dtm_k(1)
             dtm_mult(2,ikpt1+(isppol-1)*dtefield%fnkpt,istep) = dtm_k(2)

           end if

           if ( efield_flag .and. inibz == 1 .and. istep == 1)  then  !!HONG
             ikpt1i_sp=ikpt1i+(isppol-1)*dtset%nkpt
             dtefield%smat(:,:,:,ikpt1i_sp,ifor,idir) = &
&             smat_k(:,:,:)
             dtefield%sflag(:,ikpt1i_sp,ifor,idir) = &
&             sflag_k_mult(:,1)
           end if

! for IBZ k-points and first step, add
           if ((ddkflag==1 .and.((det_inv_smat == 1).or.(det_inv_smat == 11))) .and. inibz == 1 .and. istep == 1) then
             cg1(:,icg + 1: icg + npw_k1*my_nspinor*nband_k) = &
             cg1(:,icg + 1:icg + npw_k1*my_nspinor*nband_k) + &
             dkinv*cg1_k(:,1:npw_k1*my_nspinor*nband_k)
           end if

           ABI_DEALLOCATE(cgq)

         end if ! end if ikpt1 > 0 and isppol > 0

       end do ! end loop over istep

!      if (ikpt_loc <= dtefield%fmkmem) sflag_k(:) = sflag_k_mult(:,1)
       if (ikpt1 > 0) sflag_k(:) = sflag_k_mult(:,1)

     end do ! close loop over ikpt_loc (k-points, isppol)

     ABI_DEALLOCATE(smat_inv)
     ABI_DEALLOCATE(smat_k)
     ABI_DEALLOCATE(smat_k_paw)
     if (calc_epaw3_force .or. calc_epaw3_stress) then
       ABI_DEALLOCATE(dsdr)
     end if

   end do   ! close loop over ifor

!  MPI communicate stuff between everyone
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
   end if ! if parallel

!  DEBUG
!  write(std_out,*)
!  write(std_out,*)'istep = 1, nsppol =',nsppol
!  istep=1
!  isppol=1
!  do jkpt = 1, dtefield%fnkpt
!  write(std_out,'(a,i4,3e15.4,2e15.4)')'jkpt, kpt, dtm_mult(:,kpt,1)', jkpt, dtefield%fkptns(:,jkpt),  dtm_mult(:,jkpt+(isppol-1)*dtefield%fnkpt,istep)
!  end do
!  write(std_out,*)
!  write(std_out,*) "istep = 2"
!  if(berrystep>=2)then
!  istep=2
!  isppol=1
!  do jkpt = 1, dtefield%fnkpt
!  write(std_out,'(a,i4,3e15.4,2e15.4)')'jkpt, kpt, dtm_mult(:,kpt,2)', jkpt, dtefield%fkptns(:,jkpt),  dtm_mult(:,jkpt+(isppol-1)*dtefield%fnkpt,istep)
!  end do
!  end if
!  ENDDEBUG

!  ===========================================================================
!  in DDK case everything has been calculated above from finite difference
!  Now write the ddk WF to a file
!  ===========================================================================

   if (ddkflag == 1) then

     pertcase = idir + 3*natom
     response = 1
     call appdig(pertcase,dtfil%fnameabo_1wf,fiwf1o)
     ABI_ALLOCATE(resid,(mband*nkpt*nsppol))
     resid(:) = zero

     call outresid(dtset,dtset%kptns,mband,&
&                dtset%nband,nkpt,&
&                nsppol,resid)

     call outwf(cg1,dtset,psps,eig_dum,fiwf1o,hdr,kg,dtset%kptns,&
&     mband,mcg,mkmem,mpi_enreg,mpw,natom,dtset%nband,&
&     nkpt,npwarr,nsppol,&
&     occ_dum,response,dtfil%unwff2,wfs,wvl)

     ABI_DEALLOCATE(resid)
   end if  ! ddkflag == 1
! end of ddk part for this idir


!  ===========================================================================
!  Compute the Berry phase polarization
!  ===========================================================================

   if (polflag == 1) then

!    Compute the electronic Berry phase

     polb_mult(:,:)=zero
     do istep = 1,berrystep

       if(berrystep==1) then
         write(message,'(a,a)')ch10,&
&         ' Compute the electronic contribution to polarization'
         call wrtout(std_out,message,'COLL')
       else
         write(message,'(a,a,i4,a)')ch10,&
&         ' Compute the electronic contribution to polarization for a step of istep=',&
&         istep,'*dk'
         call wrtout(std_out,message,'COLL')
       end if

       if(istep /= 1) then
!        construct the strings for a step of istep*dk
!        string length
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
!        Check that the string length is a divisor of nkpt
         if(mod(dtefield%fnkpt,nkstr) /= 0) then
           write(message,'(a,i5,a,i5,a,i7)')&
&           '  For istep = ', istep,&
&           '  The string length = ',nkstr,&
&           ', is not a divisor of fnkpt =',dtefield%fnkpt
           MSG_BUG(message)
         end if
         nstr = dtefield%fnkpt/nkstr

         write(message,'(a,i1,a,i2,a,i3,a,i6)')&
&         '  berryphase_new: for direction ',idir, ' and istep ', istep, ', nkstr = ',nkstr,&
&         ', nstr = ',nstr
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
!      DEBUG
!      do istr=1,nstr
!      write(std_out,*)'string ', idxkstr_mult(:,istr)
!      end do
!      ENDBEBUG

       ABI_ALLOCATE(det_string,(2,nstr))
       ABI_ALLOCATE(polberry,(nstr))
       write(message,'(a,10x,a,10x,a)')ch10,&
&       'istr','polberry(istr)'
       call wrtout(std_out,message,'COLL')

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
&             det_string(2,istr)*dtm_imag
             dtm_k(2) = det_string(1,istr)*dtm_imag + &
&             det_string(2,istr)*dtm_real
             det_string(1:2,istr) = dtm_k(1:2)
!            DEBUG
!            write(std_out,'(a,i4,3e15.4,2e15.4)')'ikpt, kpt, dtm', ikpt, dtefield%fkptns(:,ikpt),  dtm_k
!            ENDDEBUG

             if(calc_epaw3_force) then
               do iatom = 1, natom
                 do fdir = 1, 3
                   epawf3_str(iatom,idir,fdir) = epawf3_str(iatom,idir,fdir) + &
&                   dsdr_sum(iatom,fdir,ikpt+(isppol-1)*dtefield%fnkpt)
                 end do ! end loop over fdir
               end do ! end loop over natom
             end if ! end check on calc_epaw3_force
             if(calc_epaw3_stress) then
               do iatom = 1, natom
                 do fdir = 1, 6
                   epaws3_str(iatom,idir,fdir) = epaws3_str(iatom,idir,fdir) + &
&                   dsds_sum(iatom,fdir,ikpt+(isppol-1)*dtefield%fnkpt)
                 end do ! end loop over fdir
               end do ! end loop over natom
             end if ! end check on calc_epaw3_stress

           end do

           if(calc_epaw3_force) then
             do iatom = 1, natom
               do fdir = 1, 3
                 dtefield%epawf3(iatom,idir,fdir) = dtefield%epawf3(iatom,idir,fdir) + &
&                 epawf3_str(iatom,idir,fdir)
               end do ! end loop over fdir
             end do ! end loop over natom
           end if ! end check on calc_epaw3_force
           if(calc_epaw3_stress) then
             do iatom = 1, natom
               do fdir = 1, 6
                 dtefield%epaws3(iatom,idir,fdir) = dtefield%epaws3(iatom,idir,fdir) + &
&                 epaws3_str(iatom,idir,fdir)
               end do ! end loop over fdir
             end do ! end loop over natom
           end if ! end check on calc_epaw3_stress

           det_average(:) = det_average(:) + &
&           det_string(:,istr)/dble(nstr)

         end do


!        correction to obtain a smooth logarithm of the determinant
         ABI_ALLOCATE(str_flag,(nstr))
!        DEBUG
!        since we don't have any case of non-nul Chern number,
!        we must change the det_string value "by brute force" if we want debug this
!        allocate(det_string_test(2,dtefield%nstr(idir)))
!        det_string_test(:,:)=det_string(:,:)
!        kk=0
!        det_string(1,1)=cos(2._dp*Pi*real(kk,dp)/four)
!        det_string(2,1)=sin(2._dp*Pi*real(kk,dp)/four)
!        jj=dtefield%str_neigh(1,1,idir)
!        ll=dtefield%str_neigh(2,1,idir)
!        do while (jj/=1)
!        kk=kk+1
!        det_string(1,jj)=cos(2._dp*Pi*real(kk,dp)/four)
!        det_string(2,jj)=sin(2._dp*Pi*real(kk,dp)/four)
!        det_string(1,ll)=cos(-2._dp*Pi*real(kk,dp)/four)
!        det_string(2,ll)=sin(-2._dp*Pi*real(kk,dp)/four)
!        jj=dtefield%str_neigh(1,jj,idir)
!        ll=dtefield%str_neigh(2,ll,idir)
!        enddo
!        ENDDEBUG
         if (istep==1) then
           do ineigh_str = 1,2
             str_flag(:)=0
             delta_str(:) = &
&             dtefield%coord_str(:,dtefield%str_neigh(ineigh_str,1,idir),idir) - dtefield%coord_str(:,1,idir)
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
!                write(std_out,*)'new string'
                 str_flag(istr)=1
                 call rhophi(det_string(:,istr),dphase,rho)
!                write(std_out,'(i4,e15.4,e15.4,e15.4)')istr, det_string(:,istr),dphase
                 dphase_init=dphase
                 jstr = dtefield%str_neigh(ineigh_str,istr,idir)
                 do while (istr/=jstr)
                   str_flag(jstr)=1
                   call rhophi(det_string(:,jstr),dphase_new,rho)
                   jj=nint((dphase_new-dphase)/(2._dp*Pi))
!                  DEBUG
!                  write(std_out,'(i4,e15.4,e15.4,e15.4,e15.4,i4)')jstr, det_string(:,jstr),dphase_new,dphase_new-dphase,jj
!                  ENDDEBUG
                   dphase_new=dphase_new-two*Pi*real(jj,dp)
                   if(jj/=0)then
                     write(message,'(6a)') ch10,&
&                     ' berryphase_new : WARNING -',ch10,&
&                     '  the berry phase has some huge variation in the space of strings of k-points',ch10,&
&                     '  ABINIT is trying to correct the berry phase, but it is highly experimental'
                     call wrtout(std_out,message,'PERS')
                   end if
!                  if(jj/=0)write(std_out,'(i4,e15.4,e15.4,e15.4,e15.4)')jstr, det_string(:,jstr),dphase_new,dphase_new-dphase
                   dphase=dphase_new
                   jstr=dtefield%str_neigh(ineigh_str,jstr,idir)
                 end do
!                write(std_out,*)dphase_init, dphase, (dphase-dphase_init)/(2._dp*Pi),nint((dphase-dphase_init)/(2._dp*Pi))
               end if
             end do
           end do
         end if
         ABI_DEALLOCATE(str_flag)
!        DEBUG
!        deallocate(dist_str)
!        det_string(:,:)=det_string_test(:,:)
!        deallocate(det_string_test)
!        ENDDEBUG

!        First berry phase that corresponds to det_average
!        phase0 = atan2(det_average(2),det_average(1))
         call rhophi(det_average,phase0,rho)
         det_mod = det_average(1)**2+det_average(2)**2

!        Then berry phase that corresponds to each string relative to the average
         do istr = 1, nstr

           rel_string(1) = (det_string(1,istr)*det_average(1) + &
           det_string(2,istr)*det_average(2))/det_mod
           rel_string(2) = (det_string(2,istr)*det_average(1) - &
           det_string(1,istr)*det_average(2))/det_mod
!          dphase = atan2(rel_string(2),rel_string(1))
           call rhophi(rel_string,dphase,rho)
           polberry(istr) = dtefield%sdeg*(phase0 + dphase)/two_pi
           polb_mult(isppol,istep) = polb_mult(isppol,istep) + polberry(istr)/(istep*dtefield%nstr(idir))
           polb(isppol) = zero
           do jstep=1, istep
             polb(isppol)=polb(isppol)+coef(jstep,istep)*polb_mult(isppol,jstep)
           end do

           write(message,'(10x,i6,7x,e16.9)')istr,polberry(istr)
           call wrtout(std_out,message,'COLL')

         end do

         if(berrystep>1)then
           write(message,'(9x,a,7x,e16.9,1x,a,i4,a,i4,a)')&
&           'total',polb_mult(isppol,istep),'(isppol=',isppol,', istep=',istep,')'!,ch10
           call wrtout(std_out,message,'COLL')

           write(message,'(3x,a,7x,e16.9,1x,a,i4,a,i4,a,a)')&
&           '+correction',polb(isppol),'(isppol=',isppol,', istep=1..',istep,')',ch10
           call wrtout(std_out,message,'COLL')

         else

           write(message,'(9x,a,7x,e16.9,1x,a,i4,a)')&
&           'total',polb_mult(isppol,istep),'(isppol=',isppol,')'!,ch10
           call wrtout(std_out,message,'COLL')
         end if

         polbtot = polbtot + polb(isppol)

       end do    ! isppol

!      Fold into interval [-1,1]
       polbtot = polbtot - 2_dp*nint(polbtot/2_dp)

       ABI_DEALLOCATE(det_string)
       ABI_DEALLOCATE(polberry)

!      ==========================================================================

!      Compute the ionic Berry phase

       call xred2xcart(natom,rprimd,xcart,xred)
       politot = zero
       write(message,'(a)')' Compute the ionic contributions'
       call wrtout(std_out,message,'COLL')

       write(message,'(a,2x,a,2x,a,15x,a)')ch10,&
&       'itom', 'itypat', 'polion'
       call wrtout(std_out,message,'COLL')

       do iatom = 1, natom
         itypat = typat(iatom)

!        The ionic phase can be computed much easier
         polion = zion(itypat)*xred(idir,iatom)

!        Fold into interval (-1,1)
         polion = polion - 2_dp*nint(polion/2_dp)
         politot = politot + polion
         write(message,'(2x,i2,5x,i2,10x,e16.9)') iatom,itypat,polion
         call wrtout(std_out,message,'COLL')
       end do

!      Fold into interval [-1,1] again
       politot = politot - 2_dp*nint(politot/2_dp)
       pion(idir) = politot

       write(message,'(9x,a,7x,es19.9)') 'total',politot
       call wrtout(std_out,message,'COLL')


!      ==========================================================================

!      Compute the total polarization

       poltot = politot + polbtot

       if (berrystep==1)then
         write(message,'(a,a)')ch10,&
&         ' Summary of the results'
         call wrtout(std_out,message,'COLL')
         if (unit_out /= 0) then
           call wrtout(unit_out,message,'COLL')
         end if
       else
         write(message,'(a,a,i4)')ch10,&
&         ' Summary of the results for istep =',istep
         call wrtout(std_out,message,'COLL')
         if (unit_out /= 0) then
           call wrtout(unit_out,message,'COLL')
         end if
       end if

       write(message,'(a,es19.9)')&
&       ' Electronic Berry phase ' ,polbtot
       call wrtout(std_out,message,'COLL')
       if (unit_out /= 0) then
         call wrtout(unit_out,message,'COLL')
       end if

       write(message,'(a,es19.9)') &
&       '            Ionic phase ', politot
       call wrtout(std_out,message,'COLL')
       if (unit_out /= 0) then
         call wrtout(unit_out,message,'COLL')
       end if

       write(message,'(a,es19.9)') &
&       '            Total phase ', poltot
       call wrtout(std_out,message,'COLL')
       if (unit_out /= 0) then
         call wrtout(unit_out,message,'COLL')
       end if

!      REC start
       if(abs(dtset%polcen(idir))>tol8)then
         poltot = poltot-dtset%polcen(idir)
         write(message,'(a,f15.10)') &
&         '    Translating Polarization by P0 for centrosymmetric cell: ',&
&         dtset%polcen(idir)
         call wrtout(std_out,message,'COLL')
         if (unit_out /= 0) then
           call wrtout(unit_out,message,'COLL')
         end if
       end if
!      REC end

       poltot = poltot - 2.0_dp*nint(poltot/2._dp)
       write(message,'(a,es19.9)') &
&       '    Remapping in [-1,1] ', poltot
       call wrtout(std_out,message,'COLL')
       if (unit_out /= 0) then
         call wrtout(unit_out,message,'COLL')
       end if

!      ! REC and HONG
!      =====================================================================================
!      Polarization branch control  (start)
!      -------------------------------------------------------------------------------------
!      berrysav == 0,  for non fixed D/d calculation, polarizaion is in [-1,1],done above
!      for fixed D/d calculation, choose polarization to minimize internal
!      energy, or minimize |red_efiled|. (red_dfield=red_efiled+red_ptot)
!      (d=e+p, as (26) of Stengel, Suppl.) [[cite:Stengel2009]]
!      This is default value.
!
!      berrysav == 1,  keep the polarization on the same branch, which saved in file POLSAVE
!      ======================================================================================

!      for fixed D/d calculation, choose polarization to minimize internal energy, or to minimize reduced electric field |red_efield|
       if(dtset%berrysav ==0 .and. (dtset%berryopt == 6 .or. dtset%berryopt == 7 .or. &
&       dtset%berryopt == 16 .or. dtset%berryopt == 17))  then

         jump=-nint(dtset%red_dfield(idir) - poltot)   ! red_efield = red_dfield - poltot

         if(jump /= 0)then
           write(message,'(a,i1,a,es19.9,a,i2)') &
&           ' P(',idir,') Shifted polarization branch to minimize red_efield &
&           k from ',poltot, ' by ',jump
           call wrtout(std_out,message,'COLL')
           if (unit_out /= 0) then
             call wrtout(unit_out,message,'COLL')
           end if
           poltot=poltot-jump
         end if
         pol0(idir)=poltot

       end if


!      keep the polarization on the same branch.
       if (dtset%berrysav == 1) then

!        use saved polarization to keep on same branch
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
&           ' P(',idir,') jumped to new branch. Shifting bac&
&           k from ',poltot, ' by ',jump
           call wrtout(std_out,message,'COLL')
           if (unit_out /= 0) then
             call wrtout(unit_out,message,'COLL')
           end if
           poltot=poltot-jump
         end if

         pol0(idir)=poltot

       end if

!      =====================================================================================
!      Polarization branch control  (end)
!      =====================================================================================


!      Transform the phase into a polarization
       fac = 1._dp/(gmod*dtefield%nkstr(idir))
!      !REC         fac = fac/ucvol
!      !REC         pol = fac*poltot
       red_ptot(idir)=poltot !!REC
       pol = fac*red_ptot(idir)/ucvol !!REC
       ptot(idir)=red_ptot(idir)/ucvol !!REC
       write(message,'(a,a,es19.9,a,a,a,es19.9,a,a)')ch10,&
&       '           Polarization ', pol,' (a.u. of charge)/bohr^2',ch10,&
&       '           Polarization ', pol*(e_Cb)/(Bohr_Ang*1d-10)**2,&
&       ' C/m^2',ch10
       call wrtout(std_out,message,'COLL')
       if (unit_out /= 0) then
         call wrtout(unit_out,message,'COLL')
       end if


       ABI_DEALLOCATE(idxkstr_mult)

     end do !istep

     pel(idir) = polbtot

   end if   ! if calculate polarization polflag==1

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

!-------------------------------------------------
!   Compute polarization in cartesian coordinates
!-------------------------------------------------
 if (all(dtset%rfdir(:) == 1)) then

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
 if (ddkflag == 1)  then
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

!!****f* ABINIT/update_e_field_vars
!! NAME
!! update_e_field_vars
!!
!! FUNCTION
!! This routine updates E field variables
!!
!! COPYRIGHT
!! Copyright (C) 2003-2020 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! atindx(natom)=index table for atoms, inverse of atindx (see gstate.f)
!! atindx1(natom)=index table for atoms (see gstate.f)
!! cg(2,mcg)=planewave coefficients of wavefunctions
!! dimcprj(usepaw*natom)=lmn_size for each atom
!! dtfil <type(datafiles_type)>=variables related to files
!! gmet(3,3)=metric in reciprocal space
!! gprimd(3,3)=reciprocal space dimensional primitive translations
!! idir = determines directions for derivatives computed in ctocprj (0 for all)
!! kg(3,mpw*mkmem)=reduced planewave coordinates
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mkmem=number of k points treated by this node.
!! mpw=maximum dimensioned size of npw
!! my_natom=number of atoms treated by current processor
!! natom=number of atoms in cell
!! nattyp(ntypat)=number of atoms of each type
!! ngfft(18)=contain all needed information about 3D FFT, see ~ABINIT/Infos/vargs.htm#ngfft
!! nkpt=number of k-points
!! npwarr(nkpt)=number of planewaves in basis at this k point
!! ntypat=number of types of atoms in unit cell
!! pawrhoij(natom*usepaw) <type(pawrhoij_type)> atomic occupancies
!! pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!! pwind_alloc = first dimension of pwind
!! pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!! rmet(3,3)=metric in real space
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! scfcv_level= 0 if calling before scf loop, 1 if during
!! scfcv_quit=signals whether calling during scf quit (see scfcv.F90)
!! scfcv_step=istep value of loop counter from scfcv.F90
!! ucvol=unit cell volume in bohr**3.
!! unit_out= unit for output of the results (usually the .out file of ABINIT)
!!   The option unit_out = 0 is allowed. In this case, no information is written
!!   to the output file but only to the log file.
!! usepaw= 1: use paw framework. 0:do not use paw.
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for
!!     each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real
!!     spherical harmonics
!!
!! OUTPUT
!! efield_old_cart(3)=updating cartesian values of efield (used in berryopt
!!                    6,16,17)
!! pel_cg(3)=electronic polarization
!! pelev(3)=leading order PAW contribution in pel_cg (for reporting purposes
!!          only)
!! pion(3)=ionic part of polarization
!! ptot(3)=total polarization
!! red_efield2=updating efield used in berryopt 16,17
!! red_efield2_old=updating efield used in berryopt 16.17
!! red_ptot=updating efield used in berryopt 16.17
!!
!! SIDE EFFECTS
!! Input/Output
!! dtset <type(dataset_type)>=all input variables in this dataset
!! dtefield <type(efield_type)> = efield variables
!! hdr <type(hdr_type)>=the header of wf, den and pot files
!! mpi_enreg=information about MPI parallelization
!! ptot_cart(3)=total polarization in cartesian coordinates
!! xred(3,natom)=reduced atomic coordinates
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

subroutine update_e_field_vars(atindx,atindx1,cg,dimcprj,dtefield,dtfil,dtset,&
&  efield_old_cart,gmet,gprimd,hdr,idir,kg,mcg,&
&  mkmem,mpi_enreg,mpw,my_natom,natom,nattyp,ngfft,nkpt,npwarr,ntypat,&
&  pawrhoij,pawtab,pel_cg,pelev,pion,psps,ptot,ptot_cart,pwind,&
&  pwind_alloc,pwnsfac,red_efield2,red_efield2_old,red_ptot,rmet,rprimd,&
&  scfcv_level,scfcv_quit,scfcv_step,ucvol,unit_out,&
&  usepaw,xred,ylm,ylmgr)

  !Arguments ------------------------------------
  integer, intent(in) :: idir,mcg,mkmem,mpw,my_natom,natom,nkpt,ntypat
  integer, intent(in) :: pwind_alloc,scfcv_level,scfcv_quit,scfcv_step,unit_out,usepaw
  real(dp), intent(in) :: ucvol
  type(datafiles_type), intent(in) :: dtfil
  type(pseudopotential_type),intent(in) :: psps
  type(dataset_type), intent(inout) :: dtset
  type(efield_type), intent(inout) :: dtefield
  type(hdr_type), intent(inout) :: hdr
  type(MPI_type), intent(inout) :: mpi_enreg
  !arrays
  integer, intent(in) :: atindx(natom),atindx1(natom),dimcprj(usepaw*natom)
  integer, intent(in) :: kg(3,mpw*mkmem),nattyp(ntypat)
  integer, intent(in) :: ngfft(18),npwarr(nkpt),pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: cg(2,mcg),gmet(3,3),gprimd(3,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(in) :: rmet(3,3),rprimd(3,3)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(inout) :: ptot_cart(3),xred(3,natom),efield_old_cart(3) !vz_i
  real(dp), intent(out) :: pel_cg(3),pelev(3),pion(3) !vz_i
  real(dp), intent(inout) :: red_efield2(3),red_efield2_old(3) !vz_i
  real(dp), intent(out) :: ptot(3),red_ptot(3) !vz_i
  type(pawrhoij_type), intent(in) :: pawrhoij(my_natom*usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

  !Local variables -------------------------
  !scalars
  character(len=500) :: message
  integer :: ctocprj_choice,iatom,ii,iorder_cprj,mcprj,my_nspinor,ncpgr
  integer :: optberry,usecprj
  logical :: calc_epaw3_force, calc_epaw3_stress, efield
  !arrays
  real(dp) :: efield_test_cart(3),red_efield1(3)
  real(dp),allocatable :: ph1d(:,:)
  type(pawcprj_type),allocatable :: cprj(:,:)

  ! *************************************************************************

  efield = .false.

  if ( dtset%berryopt == 4 .or. &
       & dtset%berryopt == 6 .or. &
       & dtset%berryopt == 7 .or. &
       & dtset%berryopt ==14 .or. &
       & dtset%berryopt ==16 .or. &
       & dtset%berryopt ==17 ) efield = .true.
  calc_epaw3_force = ( efield .and. dtset%optforces /= 0 .and. usepaw == 1 )
  calc_epaw3_stress = ( efield .and. dtset%optstress /= 0  .and. usepaw == 1 )

  usecprj=1; if (psps%usepaw==0)  usecprj = 0
  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
  mcprj=my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol

  ncpgr = 0
  ctocprj_choice = 1 ! no derivs
  if ( efield .and. psps%usepaw == 1) then
     ABI_DATATYPE_ALLOCATE(cprj,(dtset%natom,mcprj))
     !  finite electric field may need gradients for forces, stress
     if (calc_epaw3_force .and. .not. calc_epaw3_stress) then
        ncpgr = 3; ctocprj_choice = 2 ! derivs w.r.t. position
     else if (.not. calc_epaw3_force .and. calc_epaw3_stress) then
        ncpgr = 6; ctocprj_choice = 3 ! derivs w.r.t strain
     else if (calc_epaw3_force .and. calc_epaw3_stress) then
        ncpgr = 9; ctocprj_choice = 23 ! derivs w.r.t. position and strain
     end if
     call pawcprj_alloc(cprj,ncpgr,dimcprj)
     iatom=0 ; iorder_cprj=1 ! retain ordering of input list
     !  all arguments to ctocprj are defined already except ph1d, do that here
     ABI_ALLOCATE(ph1d,(2,3*(2*dtset%mgfft+1)*dtset%natom))
     call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred)
     call ctocprj(atindx,cg,ctocprj_choice,cprj,gmet,gprimd,iatom,idir,iorder_cprj,&
          &   dtset%istwfk,kg,dtset%kptns,mcg,mcprj,dtset%mgfft,dtset%mkmem,&
          &   mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,nattyp,dtset%nband,&
          &   dtset%natom,ngfft,dtset%nkpt,dtset%nloalg,npwarr,dtset%nspinor,&
          &   dtset%nsppol,dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,&
          &   dtset%typat,ucvol,dtfil%unpaw,xred,ylm,ylmgr)
     ABI_DEALLOCATE(ph1d)
  else
     ABI_DATATYPE_ALLOCATE(cprj,(0,0))
  end if ! end update of cprj

  if ( efield ) then ! compute polarization and if necessary store cprj in efield
     optberry=1
     pel_cg(:) = zero;pelev=zero
     call berryphase_new(atindx1,cg,cprj,dtefield,dtfil,dtset,psps,gprimd,hdr,psps%indlmn,kg,&
          &   psps%lmnmax,dtset%mband,mcg,mcprj,dtset%mkmem,mpi_enreg,dtset%mpw,my_natom,&
          &   dtset%natom,npwarr,dtset%nsppol,psps%ntypat,dtset%nkpt,optberry,pawrhoij,pawtab,&
          &   pel_cg,pelev,pion,ptot,red_ptot,pwind,&
          &   pwind_alloc,pwnsfac,rprimd,dtset%typat,ucvol,&
          &   unit_out,usecprj,psps%usepaw,xred,psps%ziontypat)

     dtefield%red_ptot1(:)=red_ptot(:)

  end if ! end compute polarization and store cprj for efield

  if (efield .and. (scfcv_level == 0) ) then ! do this before scfcv loop

     efield_old_cart(:)=dtset%efield(:)   !!HONG

     !  save this value in order to print the final value of real electric field, comparing with the desired red_fieldbar
     dtefield%efield2(:)=dtset%efield(:)

     if ( dtset%berryopt ==16 .or. dtset%berryopt ==17) then   !!HONG
        do ii=1,3
           red_efield2(ii)=zero
           red_efield2_old(ii)  =(ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,ii))
        end do
     end if

     if (dtset%berryopt == 14 .and. scfcv_quit /=1) then
        !    ! Convert polarization to cartesian coords

        ptot_cart(:)=zero
        do ii = 1,3
           ptot_cart(ii)=rprimd(ii,1)*red_ptot(1) + rprimd(ii,2)*red_ptot(2) + &
                &       rprimd(ii,3)*red_ptot(3)
        end do
        ptot_cart(:)=ptot_cart(:)/ucvol

        do ii=1,3
           dtefield%efield_dot(ii) = dot_product(dtset%efield(:),rprimd(:,ii))
        end do

        !    !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
        write(message,'(a,a)')   ch10, 'scfcv: Constant reduced ebar-field:'

        call wrtout(std_out,message,'COLL')
        call prtefield(dtset,dtefield,std_out,rprimd)

        if(dtset%prtvol>=10)then
           call wrtout(ab_out,message,'COLL')
           call prtefield(dtset,dtefield,ab_out,rprimd)
        end if

        !    updating E field
        do ii =1,3   ! desired E field
           efield_test_cart(ii)=gprimd(ii,1)*dtset%red_efieldbar(1) + &
                &       gprimd(ii,2)*dtset%red_efieldbar(2)+gprimd(ii,3)*dtset%red_efieldbar(3)
        end do

        !    if not convergence well, need to add some code here to make sure efield_test_cart(:) not change much
        dtset%efield(:) = efield_test_cart(:)

     end if  ! berryopt ==14

  end if ! end efield .and. scfcv_level 0 tasks

!!!
!!! Various printing and update steps for the different efield options
!!!

  if (efield .and. (scfcv_level == 1) ) then ! do this each scf step

     if (dtset%prtvol >= 10)then
        write(message,'(6(a),3(e16.9,2x),a,a,3(e16.9,2x))')ch10,&
             &     ' scfcv: New value of the polarization:',ch10,&
             &     ' (reduced coordinates, a. u.)',ch10,&
             &     '     Electronic berry phase:       ', (pel_cg(ii), ii = 1, 3)
        call wrtout(ab_out,message,'COLL')
        call wrtout(std_out,message,'COLL')
        if(psps%usepaw==1) then
           write(message,'(a,3(e16.9,2x))')&
                &       '     ...includes PAW on-site term: ', (pelev(ii), ii = 1, 3)
           call wrtout(ab_out,message,'COLL')
           call wrtout(std_out,message,'COLL')
        end if
        write(message,'(a,3(e16.9,2x),a,a,3(e16.9,2x))')&
             &     '     Ionic:                        ', (pion(ii), ii = 1, 3), ch10, &
             &     '     Total:                        ', (red_ptot(ii), ii = 1, 3) !!REC
        call wrtout(ab_out,message,'COLL')
        call wrtout(std_out,message,'COLL')
     end if ! end prtvol >= 10 output

     ptot_cart(:)=zero
     do ii = 1,3
        ptot_cart(ii)=rprimd(ii,1)*red_ptot(1) + rprimd(ii,2)*red_ptot(2) + &
             &     rprimd(ii,3)*red_ptot(3)
     end do
     ptot_cart(:)=ptot_cart(:)/ucvol

     !  !===================================================================================================
     !  !                                       OUTPUT  for fixed E
     !  !===================================================================================================

     if (dtset%berryopt == 4) then

        !    !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
        write(message,'(a,a)')   ch10, 'scfcv: Constant unreduced E-field:'
        call wrtout(std_out,message,'COLL')
        call prtefield(dtset,dtefield,std_out,rprimd)
        if(dtset%prtvol>=10)then
           call wrtout(ab_out,message,'COLL')
           call prtefield(dtset,dtefield,ab_out,rprimd)
        end if
     end if ! end berryopt 4 output

     !  =====================================================================================
     !  !                                      fixed D calculation
     !  !====================================================================================
     if (dtset%berryopt == 6) then
        if (scfcv_step > 1) then

           !      ! update efield taking damping into account dfield is in cartesian in dtset structure (contains input value)
           !      ! same goes for efield - update the dtset%efield value
           efield_test_cart(:)=dtset%ddamp*(dtset%dfield(:)-4.0d0*pi*ptot_cart(:))+&
                &       (1.0d0-dtset%ddamp)*efield_old_cart(:)

           !      ! test whether change in efield in any direction exceed maxestep, if so, set the
           !      ! change to maxestep instead   ! need optimized !
           do ii = 1,3

             if (dabs(efield_test_cart(ii)-efield_old_cart(ii)) > dabs(dtset%maxestep)) then

               write(std_out,'(a,a,i5)') "JH - ","  E-field component:",ii
               write(std_out,'(a,es13.5,a,es13.5,a,es13.5,a,es13.5)') " E(n)=",efield_test_cart(ii), &
&               ",    E(n-1)=",efield_old_cart(ii), ",    E(n)-E(n-1)=", efield_test_cart(ii)-efield_old_cart(ii), &
&               ",    maxestep=",dtset%maxestep


               if (efield_test_cart(ii) > efield_old_cart(ii)) then
                 efield_test_cart(ii) = efield_old_cart(ii) + dabs(dtset%maxestep)
               else
                 efield_test_cart(ii) = efield_old_cart(ii) - dabs(dtset%maxestep)
               end if
             end if
           end do

           dtset%efield(:) = efield_test_cart(:)

           !      !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
           write(message,'(a,a)')   ch10, 'scfcv: Constant unreduced D-field  - updating E-field:'
           call wrtout(std_out,message,'COLL')
           call prtefield(dtset,dtefield,std_out,rprimd)
           if(dtset%prtvol>=10)then
              call wrtout(ab_out,message,'COLL')
              call prtefield(dtset,dtefield,ab_out,rprimd)
           end if

           !      ! need to update dtset%efield_dot(:) with new value
           dtefield%efield_dot(1) = dot_product(dtset%efield(:),rprimd(:,1))
           dtefield%efield_dot(2) = dot_product(dtset%efield(:),rprimd(:,2))
           dtefield%efield_dot(3) = dot_product(dtset%efield(:),rprimd(:,3))

        else

           write(message,'(a,a)')   ch10, 'scfcv: Constant unreduced D-field  - Pre E-field:'
           call wrtout(std_out,message,'COLL')
           call prtefield(dtset,dtefield,std_out,rprimd)
           if(dtset%prtvol>=10)then
              call wrtout(ab_out,message,'COLL')
              call prtefield(dtset,dtefield,ab_out,rprimd)
           end if

        end if  ! scfcv_step >1

        efield_old_cart(:)=dtset%efield(:)
     end if  ! berryopt ==6
     !  !===================================================================================================
     !  !                                      fixed reduced d calculation
     !  !===================================================================================================
     if (dtset%berryopt == 16) then

        if (scfcv_step > 1) then
           !      ! update efield taking damping into account reduced red_dfield
           !      red_efield2 is reduced electric field, defined by Eq.(25) of Nat. Phys. suppl. (2009) [[cite:Stengel2009]]

           red_efield2(:)=dtset%ddamp*(dtset%red_dfield(:)-red_ptot(:))+ (1.0d0-dtset%ddamp)*red_efield2_old(:)

           !      to calculate unreduced E
           efield_test_cart(:)=(4*pi/ucvol)*(rprimd(:,1)*red_efield2(1)+rprimd(:,2)*red_efield2(2)+rprimd(:,3)*red_efield2(3))

           !      ! test whether change in efield in any direction exceed maxestep, if so, set the
           !      ! change to maxestep instead   ! need optimized !
           do ii = 1,3

             if (dabs(efield_test_cart(ii)-efield_old_cart(ii)) > dabs(dtset%maxestep)) then

               write(std_out,'(a,a,i5)') "JH - ","  E-field component:",ii
               write(std_out,'(a,es13.5,a,es13.5,a,es13.5,a,es13.5)') " E(n)=",efield_test_cart(ii), &
&               ",    E(n-1)=",efield_old_cart(ii), ",    E(n)-E(n-1)=", efield_test_cart(ii)-efield_old_cart(ii), &
&               ",    maxestep=",dtset%maxestep

               if (efield_test_cart(ii) > efield_old_cart(ii)) then
                 efield_test_cart(ii) = efield_old_cart(ii) + dabs(dtset%maxestep)
               else
                 efield_test_cart(ii) = efield_old_cart(ii) - dabs(dtset%maxestep)
               end if
             end if
           end do

           dtset%efield(:) = efield_test_cart(:)

           !      !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
           write(message,'(a,a)')   ch10, 'scfcv: Constant reduced d-field  - updating E-field:'
           call wrtout(std_out,message,'COLL')
           call prtefield(dtset,dtefield,std_out,rprimd)
           if(dtset%prtvol>=10)then
              call wrtout(ab_out,message,'COLL')
              call prtefield(dtset,dtefield,ab_out,rprimd)
           end if

           !      ! need to update dtset%efield_dot(:) with new value
           !      ! This needs to be deleted  when efield_dot is deleted
           dtefield%efield_dot(1) = dot_product(dtset%efield(:),rprimd(:,1))
           dtefield%efield_dot(2) = dot_product(dtset%efield(:),rprimd(:,2))
           dtefield%efield_dot(3) = dot_product(dtset%efield(:),rprimd(:,3))

        else

           write(message,'(a,a)')   ch10, 'scfcv: Constant reduced d-field  - Pre E-field:'
           call wrtout(std_out,message,'COLL')
           call prtefield(dtset,dtefield,std_out,rprimd)
           if(dtset%prtvol>=10)then
              call wrtout(ab_out,message,'COLL')
              call prtefield(dtset,dtefield,ab_out,rprimd)
           end if

        end if  ! scfcv_step > 1

        efield_old_cart(:)=dtset%efield(:)
        red_efield2_old(:)=red_efield2(:)
     end if  ! berryopt ==16


     !  !===================================================================================================
     !  !                                      fixed reduced d and ebar calculation (mixed BC)
     !  !===================================================================================================
     if (dtset%berryopt == 17) then

        if (scfcv_step > 1) then
           !      ! update efield taking damping into account reduced red_dfield
           !      red_efield1 and red_efield2 is reduced electric field, defined by Eq.(25) of Nat. Phys. suppl. (2009) [[cite:Stengel1999]]
           !      red_efield1 for fixed ebar, red_efield2 for fixed d calculation

           !      save this value in order to print the final value of real electric field, comparing with the desired red_fieldbar
           dtefield%efield2(:)=dtset%efield(:)

           !      write(*,'(a,3i4)') "jfielddir=", (dtset%jfielddir(ii),ii=1,3)

           do ii=1,3
              if (dtset%jfielddir(ii) ==2 ) then    ! direction under fixed d
                 dtset%red_efieldbar(ii) = dot_product(dtset%efield(:),rprimd(:,ii)) !  update ebar which is not fixed
                 dtefield%efield_dot(ii) = dot_product(dtset%efield(:),rprimd(:,ii))
                 red_efield2(ii)=dtset%ddamp*(dtset%red_dfield(ii) - red_ptot(ii)) +  &
                      &           (1.0d0-dtset%ddamp)*red_efield2_old(ii)         ! d(ii) is fixed, update e(ii)  may need ddamping here

                 !          write(message,'(a,a,i5,a,i5)')   ch10, 'direction  ', ii,'   for fixed d, value is (2)  ', dtset%jfielddir(ii)
                 !          call wrtout(ab_out,message,'COLL')
                 !          call wrtout(std_out,message,'COLL')

              else if (dtset%jfielddir(ii) ==1 ) then   ! direction under fixed ebar
                 red_efield2(ii)= (ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,ii)) !  update e which is not fixed
                 dtset%red_dfield(ii)=red_ptot(ii) +  (ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,ii))  ! update d

                 !          write(message,'(a,a,i5,a,i5)')   ch10, 'direction  ', ii,'   for fixed ebar, value is (1)  ', dtset%jfielddir(ii)
                 !          call wrtout(ab_out,message,'COLL')
                 !          call wrtout(std_out,message,'COLL')

              end if
           end do

           do ii=1,3
              red_efield1(ii)  =(ucvol/(4*pi))*dot_product(dtset%red_efieldbar(:),gmet(:,ii))
           end do


           dtset%red_efield(:)=(red_efield1(:) + red_efield2(:))/2.0d0 ! average reduced efield,
           !      one is from fixed ebar part,
           !      the other is from fixed d part.
           !      This may need to be optimized !!

           write(message,'(a,a,a,a,3(es16.9,2x),a)')   ch10, 'Reduced efield from fixed ebar:', ch10, &
                &       '       e:  ', (red_efield1(ii),ii=1,3), ch10

           !      call wrtout(ab_out,message,'COLL')
           call wrtout(std_out,message,'COLL')

           write(message,'(a,a,a,a,3(es16.9,2x),a)')   ch10, 'Reduced efield from fixed d:', ch10, &
                &       '       e:  ', (red_efield2(ii),ii=1,3), ch10

           !      call wrtout(ab_out,message,'COLL')
           call wrtout(std_out,message,'COLL')

           write(message,'(a,a,a,a,3(es16.9,2x),a)')   ch10, 'Average reduced efield:', ch10, &
                &       '       e:  ', (dtset%red_efield(ii),ii=1,3), ch10

           !      call wrtout(ab_out,message,'COLL')
           call wrtout(std_out,message,'COLL')

           !      to calculate unreduced E
           do ii=1,3
              efield_test_cart(ii)  = (4*pi/ucvol)* dot_product(dtset%red_efield(:),rprimd(:,ii))
           end do

           !      ! test whether change in efield in any direction exceed maxestep, if so, set the
           !      ! change to maxestep instead   ! need optimized !
           do ii = 1,3
             if (dabs(efield_test_cart(ii)-efield_old_cart(ii)) > dabs(dtset%maxestep)) then

               write(std_out,'(a,a,i5)') "JH - ","  E-field component:",ii
               write(std_out,'(a,es13.5,a,es13.5,a,es13.5,a,es13.5)') " E(n)=",efield_test_cart(ii), &
&               ",    E(n-1)=",efield_old_cart(ii), ",    E(n)-E(n-1)=", efield_test_cart(ii)-efield_old_cart(ii), &
&               ",    maxestep=",dtset%maxestep

               if (efield_test_cart(ii) > efield_old_cart(ii)) then
                 efield_test_cart(ii) = efield_old_cart(ii) + dabs(dtset%maxestep)
               else
                 efield_test_cart(ii) = efield_old_cart(ii) - dabs(dtset%maxestep)
               end if
             end if
           end do

           dtset%efield(:) = efield_test_cart(:)

           !      !write the field parameters: D, E, P, d, e, p, dbar, ebar, pbar
           write(message,'(a,a)')   ch10, 'scfcv: Constant reduced ebar and d-field  - updating E-field:'
           call wrtout(std_out,message,'COLL')
           call prtefield(dtset,dtefield,std_out,rprimd)
           if(dtset%prtvol>=10)then
              call wrtout(ab_out,message,'COLL')
              call prtefield(dtset,dtefield,ab_out,rprimd)
           end if


           !      ! need to update dtset%efield_dot(:) with new value
           !      ! This needs to be deleted  when efield_dot is deleted
           dtefield%efield_dot(1) = dot_product(dtset%efield(:),rprimd(:,1))
           dtefield%efield_dot(2) = dot_product(dtset%efield(:),rprimd(:,2))
           dtefield%efield_dot(3) = dot_product(dtset%efield(:),rprimd(:,3))

        else

           write(message,'(a,a)')   ch10, 'scfcv: Constant reduced ebar and d-field  - Pre E-field:'
           call wrtout(std_out,message,'COLL')
           call prtefield(dtset,dtefield,std_out,rprimd)
           if(dtset%prtvol>=10)then
              call wrtout(ab_out,message,'COLL')
              call prtefield(dtset,dtefield,ab_out,rprimd)
           end if

        end if  ! scfcv_step > 1

        efield_old_cart(:)=dtset%efield(:)
        red_efield2_old(:)=red_efield2(:)

     end if  ! berryopt ==17

  end if ! end efield .and. scfcv_level 1 tasks

  !deallocate cprj
  if ( efield .and. psps%usepaw == 1) then
     call pawcprj_free(cprj)
  end if
  ABI_DATATYPE_DEALLOCATE(cprj)

end subroutine update_e_field_vars
!!***

!!****f* ABINIT/prtefield
!!
!! NAME
!! prtefield
!!
!! FUNCTION
!! Print components of electric field, displacement field and polarization in nice format
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, LBoeri, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryopt
!!   | efield
!!   | dfield
!!   | red_efield
!!   | red_efieldbar
!!   | red_dfield
!!  dtefield <type(efield_type)>
!!   | efield2
!!   | red_ptot1
!!  iunit = unit number to which the data is printed
!!  rprimd
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      gstate,update_e_field_vars
!!
!! CHILDREN
!!      metric,wrtout
!!
!! SOURCE

subroutine prtefield(dtset,dtefield,iunit,rprimd)

  !Arguments ------------------------------------
  integer :: iunit
  real(dp),intent(in) :: rprimd(3,3)
  type(efield_type),intent(in) :: dtefield
  type(dataset_type),intent(inout) :: dtset

  !Local variables-------------------------------
  ! Do not modify the length of this string
  !scalars
  integer :: idir,ii
  character(len=1500) :: message
  character(len=7)   :: flag_field(3)

  real(dp) ::    ucvol
  ! arrays
  real(dp) :: ptot_cart(3),gmet(3,3),gprimd(3,3),rmet(3,3),red_pbar(3),red_dbar(3),red_dfieldbar(3)
  real(dp) :: red_efieldbar_lc(3),red_efield_lc(3)


  ! *************************************************************************

  !DEBUG
  !write(iout, '(a)') ' prtefield : enter '
  !ENDDEBUG

  !write here

  call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)


  ptot_cart(:)=zero
  do idir = 1,3
     ptot_cart(idir)=rprimd(idir,1)*dtefield%red_ptot1(1) + rprimd(idir,2)*dtefield%red_ptot1(2) + &
          &   rprimd(idir,3)*dtefield%red_ptot1(3)
  end do
  ptot_cart(:)=ptot_cart(:)/ucvol

  if (dtset%berryopt == 4) then

     !  to calculate e Eq.(25)
     do idir=1,3
        dtset%red_efield(idir)  =(ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,idir))
     end do

     !  to calculate ebar Eq.(25)
     do idir=1,3
        dtset%red_efieldbar(idir)  =dot_product(dtset%efield(:),rprimd(:,idir))
     end do


     !  to calculate pbar
     do idir=1,3
        red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
     end do

     !MGNAG: This message is too long and causes
     ! Runtime Error: wrtout_cpp.f90, line 893: Buffer overflow on output
     ! with NAG in test seq_tsv6_125 where we write to std_out!
     ! I cannot change the RECLEN of std_out!

     write(message,'(a,a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ' (a. u.)', ch10,&
          &   '       E:  ', (dtset%efield(ii), ii=1,3), ch10, &
          &   '       P:  ', (ptot_cart(ii), ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10,&
          &   '    ebar:  ', (dtset%red_efieldbar(ii),ii=1,3), ch10, &   !!HONG need to change
          &  '    pbar:  ', (red_pbar(ii),ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10,&
          &   '       e:  ', (dtset%red_efield(ii),ii=1,3), ch10, &
          &   '       p:  ', (dtefield%red_ptot1(ii), ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a)') ch10,&
          &   ' (S.I.), that is V/m for E, and C/m^2 for P', ch10, &
          &   '-      E:  ', (dtset%efield(ii)*(Ha_J/(e_Cb*Bohr_Ang*1d-10)), ii=1,3), ch10, &    !(Ha_J/(e_Cb*Bohr_Ang*1d-10))= 5.14220652*1d+11
          &  '       P:  ', (ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2, ii=1,3),ch10
     call wrtout(iunit,message,'COLL')

  end if ! berryopt ==4


  if (dtset%berryopt == 6) then
     do idir=1,3
        dtset%red_efield(idir)  =(ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,idir))
     end do

     !  to calculate ebar   !! Need to be changed
     do idir=1,3
        dtset%red_efieldbar(idir)  = dot_product(dtset%efield(:),rprimd(:,idir))
     end do

     !  to calculate red_pbar
     do idir=1,3
        red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
     end do

     !  to calculate red_dbar
     do idir=1,3
        red_dbar(idir)  = dot_product(dtset%dfield(:),rprimd(:,idir))
     end do

     !  to calculate d
     do idir=1,3
        dtset%red_dfield(idir)  =(ucvol/(4*pi))*dot_product(dtset%dfield(:),gprimd(:,idir))
     end do

     write(message,'(a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ' (a. u.)', ch10,&
          &   '       e:  ', (dtset%red_efield(ii),ii=1,3), ch10, &
          &   '       p:  ', (dtefield%red_ptot1(ii), ii=1,3), ch10, &
          &   '       d:  ', (dtset%red_dfield(ii),ii = 1, 3), ch10, &
          &   ' e  +  p:  ', (dtset%red_efield(ii)+dtefield%red_ptot1(ii),ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10,&
          &   '    ebar:  ', (dtset%red_efieldbar(ii),ii=1,3), ch10, &   !!HONG need to change
          &  '    pbar:  ', (red_pbar(ii),ii=1,3), ch10, &
          &   '    dbar:  ', (red_dbar(ii),ii=1,3), ch10, &
          &   ' eba+pba:  ', (dtset%red_efieldbar(ii)+red_pbar(ii),ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
          &   '       E:  ', (dtset%efield(ii), ii=1,3), ch10, &
          &   '       P:  ', (ptot_cart(ii), ii=1,3), ch10, &
          &   '       D:  ', (dtset%dfield(ii),ii = 1, 3), ch10, &
          &   'E+4*pi*P:  ', (dtset%efield(ii)+4.0d0*pi*ptot_cart(ii),ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a)') ch10,&
&     ' (S.I.), that is V/m for E, and C/m^2 for P and D', ch10, &
&     '-      E:  ', (dtset%efield(ii)*(Ha_J/(e_Cb*Bohr_Ang*1d-10)), ii=1,3), ch10,& !(Ha_J/(e_Cb*Bohr_Ang*1d-10))= 5.14220652*1d+11
&     '       P:  ', (ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2, ii=1,3), ch10,&
&     '       D:  ', ((1.0d0/(4*pi))*dtset%dfield(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2,ii = 1, 3),ch10,&
&     'eps0*E+P:  ', (dtset%efield(ii)*eps0*(Ha_J/(e_Cb*Bohr_Ang*1d-10))+ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2,ii=1,3),ch10
        ! eps0*(Ha_J/(e_Cb*Bohr_Ang*1d-10))=8.854187817620*5.14220652*1d-1
     call wrtout(iunit,message,'COLL')

     !MGNAG Runtime Error: wrtout_cpp.f90, line 896: Buffer overflow on output

  end if  ! berryopt ==6


  if (dtset%berryopt == 14)  then

     do idir=1,3   ! ebar local
        red_efieldbar_lc(idir)=dot_product(dtefield%efield2(:),rprimd(:,idir))
     end do


     do idir=1,3
        dtset%red_efield(idir)  =(ucvol/(4*pi))*dot_product(dtefield%efield2(:),gprimd(:,idir))
     end do

     !  to calculate pbar
     do idir=1,3
        red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
     end do

     write(message,'(a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))')  ' (a. u.)', ch10,&
          &   '   ebar0:  ', (dtset%red_efieldbar(ii),ii=1,3), ch10, &
          &   '    ebar:  ', (red_efieldbar_lc(ii),ii=1,3), ch10, &
          &   '    pbar:  ', (red_pbar(ii),ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
          &   '       e:  ', (dtset%red_efield(ii),ii=1,3), ch10, &
          &   '       p:  ', (dtefield%red_ptot1(ii), ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
          &   '       E:  ', (dtefield%efield2(ii), ii=1,3), ch10, &
          &   '       P:  ', (ptot_cart(ii), ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a)') ch10, &
          &   ' (S.I.), that is V/m for E, and C/m^2 for P', ch10, &
          &   '-      E:  ', (dtefield%efield2(ii)*(Ha_J/(e_Cb*Bohr_Ang*1d-10)), ii=1,3), ch10, &    !(Ha_J/(e_Cb*Bohr_Ang*1d-10))= 5.14220652*1d+11
          &  '       P:  ', (ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2, ii=1,3),ch10
     call wrtout(iunit,message,'COLL')


  end if  ! berryopt ==14


  if (dtset%berryopt == 16) then

     !  to calculate e Eq.(25)
     do idir=1,3
        dtset%red_efield(idir)  =(ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,idir))
     end do

     !  to calculate ebar
     do idir=1,3
        dtset%red_efieldbar(idir)  = dot_product(dtset%efield(:),rprimd(:,idir))
     end do

     !  to calculate pbar
     do idir=1,3
        red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
     end do

     !  to calculate dbar
     do idir=1,3
        red_dfieldbar(idir)  = (4*pi/ucvol)*dot_product(dtset%red_dfield(:),rmet(:,idir))
     end do

     !  to calculate D
     do idir=1,3
        dtset%dfield(idir)  =(4*pi/ucvol)*dot_product(dtset%red_dfield(:),rprimd(:,idir))
     end do

     write(message,'(a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ' (a. u.)', ch10,&
          &   '       e:  ', (dtset%red_efield(ii),ii=1,3), ch10, &
          &   '       p:  ', (dtefield%red_ptot1(ii), ii=1,3), ch10, &
          &   '       d:  ', (dtset%red_dfield(ii),ii = 1, 3), ch10, &
          &   ' e  +  p:  ', (dtset%red_efield(ii)+dtefield%red_ptot1(ii),ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
          &   '    ebar:  ', (dtset%red_efieldbar(ii),ii=1,3), ch10, &
          &   '    pbar:  ', (red_pbar(ii),ii=1,3), ch10, &
          &   '    dbar:  ', (red_dfieldbar(ii),ii=1,3), ch10, &
          &   ' eba+pba:  ', (dtset%red_efieldbar(ii)+red_pbar(ii),ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
          &   '       E:  ', (dtset%efield(ii), ii=1,3), ch10, &
          &   '       P:  ', (ptot_cart(ii), ii=1,3), ch10, &
          &   '       D:  ', (dtset%dfield(ii),ii = 1, 3), ch10, &
          &   'E+4*pi*P:  ', (dtset%efield(ii)+4.0d0*pi*ptot_cart(ii),ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a)') ch10, &
&     ' (S.I.), that is V/m for E, and C/m^2 for P and D', ch10, &
&     '-      E:  ', (dtset%efield(ii)*(Ha_J/(e_Cb*Bohr_Ang*1d-10)), ii=1,3), ch10, &    !(Ha_J/(e_Cb*Bohr_Ang*1d-10))= 5.14220652*1d+11
&     '       P:  ', (ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2, ii=1,3), ch10, &
&     '       D:  ', ((1.0d0/(4*pi))*dtset%dfield(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2,ii = 1, 3), ch10, &
&     'eps0*E+P:  ', (dtset%efield(ii)*eps0*(Ha_J/(e_Cb*Bohr_Ang*1d-10))+ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2,ii=1,3),ch10
     call wrtout(iunit,message,'COLL')

  end if  ! berryopt ==16

  if (dtset%berryopt == 17) then

     do idir=1,3   ! ebar local
        red_efieldbar_lc(idir)=dot_product(dtefield%efield2(:),rprimd(:,idir))
     end do

     do idir=1,3
        dtset%red_efield(idir)  =(ucvol/(4*pi))*dot_product(dtefield%efield2(:),gprimd(:,idir))
     end do

     !  to calculate pbar
     do idir=1,3
        red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
     end do


     !  do idir=1,3
     !  if (dtset%rfdir(idir)==1) then   ! fixed ebar
     !  red_efieldbar_lc(idir)=dot_product(dtefield%efield2(:),rprimd(:,idir))  ! local efieldbar
     !  dtset%red_efield(idir)  =(ucvol/(4*pi))*dot_product(dtefield%efield2(:),gprimd(:,idir))
     !  red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
     !  dtset%dfield(idir)=dtefield%efield2(idir)+4*pi*ptot_cart(idir)
     !  dtset%red_dfield(idir)=dtset%red_efield+dtefield%red_ptot1(idir)
     !  dtset%red_dfieldbar(idir)=red_efieldbar_lc(idir)+red_pbar(idir)
     !  E_lc(idir)=dtefield%efield2(idir)
     !  e_lc(idir)=red_efieldbar_lc(idir)
     !  ebar_lc(idir)=dtset%red_efieldbar(idir)
     !  else if (dtset%rfdir(idir)==2) then ! fixed d
     !  dtset%red_efield(idir)  =(ucvol/(4*pi))*dot_product(dtefield%efield2(:),gprimd(:,idir))
     !  dtset%red_efieldbar(idir)  = dot_product(dtefield%efield2(:),rprimd(:,idir))
     !  red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
     !  red_dfieldbar(idir)  = (4*pi/ucvol)*dot_product(dtset%red_dfield(:),rmet(:,idir))
     !  dtset%dfield(idir)  =(4*pi/ucvol)*dot_product(dtset%red_dfield(:),rprimd(:,idir))
     !  E_lc(idir)=dtefield%efield2(idir)
     !  e_lc(idir)=dtset%red_efield(idir)
     !  ebar_lc(idir)=dtset%red_efieldbar(idir)
     !  end if
     !  enddo



     do idir=1,3
        red_efield_lc(idir)= (ucvol/(4*pi))*dot_product(dtefield%efield2(:),gprimd(:,idir))
        red_efieldbar_lc(idir)=dot_product(dtefield%efield2(:),rprimd(:,idir))  ! local efieldbar
        red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
        red_dfieldbar(idir)  = (4*pi/ucvol)*dot_product(dtset%red_dfield(:),rmet(:,idir))
        dtset%dfield(idir)  =(4*pi/ucvol)*dot_product(dtset%red_dfield(:),rprimd(:,idir))

     end do

     do idir=1,3
        if(dtset%jfielddir(idir)==1) then
           flag_field(idir)="E-field"
        else
           flag_field(idir)="D-field"
        end if
     end do

     write(message,'(a,a,a,6x,a,11x,a,11x,a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') &
          &   ' (a. u.)', ch10,&
          &   '           ', (flag_field(ii),ii=1,3),ch10, &
          &   '   ebar0:  ', (dtset%red_efieldbar(ii),ii=1,3), ch10, &
          &   '    ebar:  ', (red_efieldbar_lc(ii),ii=1,3), ch10, &
          &   '       d:  ', (dtset%red_dfield(ii),ii = 1, 3), ch10, &
          &   ' e  +  p:  ', (red_efield_lc(ii)+dtefield%red_ptot1(ii),ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
          &   '       e:  ', (red_efield_lc(ii),ii=1,3), ch10, &
          &   '       p:  ', (dtefield%red_ptot1(ii), ii=1,3), ch10, &
          &   '       d:  ', (dtset%red_dfield(ii),ii = 1, 3), ch10, &
          &   ' e  +  p:  ', (red_efield_lc(ii)+dtefield%red_ptot1(ii),ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
          &   '    ebar:  ', (red_efieldbar_lc(ii),ii=1,3), ch10, &
          &   '    pbar:  ', (red_pbar(ii),ii=1,3), ch10, &
          &   '    dbar:  ', (red_dfieldbar(ii),ii=1,3), ch10, &
          &   ' eba+pba:  ', (red_efieldbar_lc(ii)+red_pbar(ii),ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
          &   '       E:  ', (dtefield%efield2(ii), ii=1,3), ch10, &
          &   '       P:  ', (ptot_cart(ii), ii=1,3), ch10, &
          &   '       D:  ', (dtset%dfield(ii),ii = 1, 3), ch10, &
          &   'E+4*pi*P:  ', (dtset%efield(ii)+4.0d0*pi*ptot_cart(ii),ii=1,3)
     call wrtout(iunit,message,'COLL')

     write(message,'(a,a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a)')  ch10, &
&     ' (S.I.), that is V/m for E, and C/m^2 for P and D', ch10, &
&     '       E:  ', (dtefield%efield2(ii)*(Ha_J/(e_Cb*Bohr_Ang*1d-10)), ii=1,3), ch10,& !(Ha_J/(e_Cb*Bohr_Ang*1d-10))= 5.14220652*1d+11
&     '       P:  ', (ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2, ii=1,3), ch10, &
&     '       D:  ', ((1.0d0/(4*pi))*dtset%dfield(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2,ii = 1, 3), ch10, &
&     'eps0*E+P:  ', (dtefield%efield2(ii)*eps0*(Ha_J/(e_Cb*Bohr_Ang*1d-10))+ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2,ii=1,3),ch10
     call wrtout(iunit,message,'COLL')

  end if  ! berryopt ==17

end subroutine prtefield
!!***

!!****f* ABINIT/init_e_field_vars
!! NAME
!! init_e_field_vars
!!
!! FUNCTION
!! Initialization of variables and data structures used in polarization
!! calculations
!!
!! COPYRIGHT
!! Copyright (C) 2004-2020 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables in this dataset
!!  gmet(3,3) = reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3) = primitive translations in recip space
!!  kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!!  mpi_enreg=information about MPI parallelization
!!  npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!!  occ(mband*nkpt*nsppol) = occup number for each band at each k point
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3) = dimensional primitive vectors
!!  symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!    reciprocal space primitive translations
!!  xred(3,natom) = location of atoms in reduced units
!!
!! OUTPUT
!!  dtefield <type(efield_type)> :: initialized polarization variables
!!  pwind(pwind_alloc,2,3) = array used to compute the overlap matrix smat
!!                         between k-points k and k +- dk where dk is
!!                         parallel to the direction idir
!!  pwind_alloc = first dimension of pwind and pwnsfac
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!
!! SIDE EFFECTS
!!
!! TO DO
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      initberry
!!
!! SOURCE

subroutine init_e_field_vars(dtefield,dtset,gmet,gprimd,kg,&
     &              mpi_enreg,npwarr,occ,pawang,pawrad,pawtab,psps,&
     &              pwind,pwind_alloc,pwnsfac,rprimd,symrec,xred)

  !Arguments ------------------------------------
  !scalars
  integer,intent(out) :: pwind_alloc
  type(MPI_type),intent(inout) :: mpi_enreg
  type(dataset_type),intent(inout) :: dtset
  type(efield_type),intent(inout) :: dtefield !vz_i needs efield2
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  !arrays
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),npwarr(dtset%nkpt)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  integer,pointer :: pwind(:,:,:)
  real(dp),intent(in) :: gmet(3,3),gprimd(3,3),occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in) :: rprimd(3,3),xred(3,dtset%natom)
  real(dp),pointer :: pwnsfac(:,:)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

  !Local variables-------------------------------
  logical :: initfield
  !scalars

  ! *************************************************************************

  initfield = .false.

  !initialization
  dtefield%has_qijb = 0
  dtefield%has_epawf3 = 0
  dtefield%has_epaws3 = 0
  dtefield%has_expibi = 0
  dtefield%has_rij = 0
  dtefield%usecprj = 0
  dtefield%berryopt = 0

  if ((dtset%berryopt < 0).or.(dtset%berryopt == 4) .or. (dtset%berryopt == 6) .or.(dtset%berryopt == 7) .or. &
       & (dtset%berryopt == 14) .or.(dtset%berryopt == 16) .or.(dtset%berryopt == 17)) then
     nullify(pwind,pwnsfac)
     call initberry(dtefield,dtset,gmet,gprimd,kg,&
          &   dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,&
          &   dtset%natom,dtset%nkpt,npwarr,dtset%nsppol,&
          &   dtset%nsym,dtset%ntypat,occ,pawang,pawrad,pawtab,&
          &   psps,pwind,pwind_alloc,pwnsfac,rprimd,symrec,&
          &   dtset%typat,psps%usepaw,xred)
     initfield = .true.
  end if

  if (.not. initfield .and. dtset%orbmag == 0) then
     ! initorbmag.F90 also allocates pwind and pwnsfac
     pwind_alloc = 1
     ABI_ALLOCATE(pwind,(pwind_alloc,2,3))
     ABI_ALLOCATE(pwnsfac,(2,pwind_alloc))
     pwind(:,:,:)=0
     pwnsfac(:,:)=zero
  end if

end subroutine init_e_field_vars
!!***

!!****f* ABINIT/initberry
!! NAME
!! initberry
!!
!! FUNCTION
!! Initialization of Berryphase calculation of the polarization, the
!! ddk and the response of an insulator to a homogenous electric field.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2020 ABINIT group (MVeithen).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables in this dataset
!!  gmet(3,3) = reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3) = primitive translations in recip space
!!  kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!!  mband = maximum number of bands
!!  mkmem = maximum number of k-points in core memory
!!  mpw = maximum number of plane waves
!!  natom = number of atoms in unit cell
!!  nkpt = number of k points
!!  npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!!  nsppol = 1 for unpolarized, 2 for spin-polarized
!!  nsym = number of symmetry operations
!!  ntypat = number of types of atoms in unit cell
!!  occ(mband*nkpt*nsppol) = occup number for each band at each k point
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3) = dimensional primitive vectors
!!  symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!    reciprocal space primitive translations
!!  typat = typat(natom) list of atom types
!!  usepaw = flag for PAW (1 PAW, 0 NCPP)
!!  xred(3,natom) = location of atoms in reduced units
!!
!! OUTPUT
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations
!!  pwind(pwind_alloc,2,3) = array used to compute the overlap matrix smat
!!                         between k-points k and k +- dk where dk is
!!                         parallel to the direction idir
!!    jpw = pwind(ipw,ifor,idir)
!!      * ipw = index of plane wave vector G for a given k-point k
!!      * ifor = 1: k + dk
!!               2: k - dk
!!      * idir = direction of the polarization/ddk calculation [dk(idir)
!!               is the only non-zero element of dk(:)]
!!      * jpw = index of plane wave vector G (+dG) at k +- dk
!!              where dG is a shift of one reciprocal lattice vector
!!              (required to close the strings of k-points using the
!!               periodic gauge condition)
!!    In case a G-vector of the basis sphere of plane waves at k
!!    does not belong to the basis sphere of plane waves at k+dk, jpw = 0.
!!   pwind_alloc = first dimension of pwind and pwnsfac
!!   pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!
!! SIDE EFFECTS
!!  mpi_enreg = informations about MPI parallelization
!!    kptdstrb(nproc,nneighbour,fmkmem_max*nsppol) : Array required
!!      by berryphase_new.f for MPI // over k-points. Defined
!!      for k-points in the fBZ
!!      but for k-points in the iBZ. Used by vtorho.f
!!           nproc = number of cpus
!!           nneighbour = number of neighbours for each k-point (= 6)
!!
!! PARENTS
!!      init_e_field_vars
!!
!! CHILDREN
!!      expibi,kpgsph,listkk,metric,pawcprj_alloc,pawcprj_getdim,qijb_kk
!!      setsym_ylm,smpbz,symatm,timab,wrtout,xmpi_max,xmpi_sum
!!
!! SOURCE

subroutine initberry(dtefield,dtset,gmet,gprimd,kg,mband,&
     &              mkmem,mpi_enreg,mpw,natom,nkpt,npwarr,nsppol,&
     &              nsym,ntypat,occ,pawang,pawrad,pawtab,psps,&
     &              pwind,pwind_alloc,pwnsfac,&
     &              rprimd,symrec,typat,usepaw,xred)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: mband,mkmem,mpw,natom,nkpt,nsppol,nsym,ntypat,usepaw
  integer,intent(out) :: pwind_alloc
  type(MPI_type),intent(inout) :: mpi_enreg
  type(dataset_type),intent(inout) :: dtset
  type(efield_type),intent(inout) :: dtefield !vz_i
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  !arrays
  integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt)
  integer,intent(in) :: symrec(3,3,nsym),typat(natom)
  integer,pointer :: pwind(:,:,:)
  real(dp),intent(in) :: gmet(3,3),gprimd(3,3),occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
  real(dp),pointer :: pwnsfac(:,:)
  type(pawrad_type),intent(in) :: pawrad(ntypat)
  type(pawtab_type),intent(in) :: pawtab(ntypat)

  !Local variables-------------------------------
  !scalars
  integer :: exchn2n3d,flag,flag_kpt,fnkpt_computed,iband,icg,icprj
  integer :: idir,idum,idum1,ierr,ifor,ikg,ikg1,ikpt,ikpt1,ikpt1f
  integer :: ikpt1i,ikpt2,ikpt_loc,ikptf,ikpti,ikstr,index,ineigh,ipw,ipwnsfac
  integer :: isppol,istr,istwf_k,isym,isym1,itrs,itypat,iunmark,jpw,klmn,lmax,lmn2_size_max
  integer :: me,me_g0,mkmem_,my_nspinor,nband_k,mband_occ_k,ncpgr,nkstr,nproc,npw_k,npw_k1,spaceComm
  integer :: option, brav, mkpt, nkptlatt
  integer :: jstr,ii,jj,isign
  integer :: dk_flag, coord1, coord2
  integer :: mult
  real(dp) :: c1,ecut_eff,eg,eg_ev,rdum,diffk1,diffk2,diffk3
  real(dp) :: dist_, max_dist, last_dist, dist,kpt_shifted1,kpt_shifted2,kpt_shifted3
  real(dp) :: gprimdlc(3,3),rmetllc(3,3),gmetllc(3,3),ucvol_local
  ! gprimd(3,3) = inverse of rprimd
  ! rmetlcl(3,3)=real-space metric (same as rmet in metric.F90)
  ! gmetlcl(3,3)= same as gmet in metric.F90
  ! ucvol = volume of the unit cell in Bohr**3
  character(len=500) :: message
  logical :: calc_epaw3_force,calc_epaw3_stress,fieldflag
  !arrays
  integer :: dg(3),iadum(3),iadum1(3),neigh(6)
  integer,allocatable :: dimlmn(:),kg1_k(:,:),kpt_mark(:),nattyp_dum(:)
  real(dp) :: diffk(3),dk(3),dum33(3,3),eg_dir(3)
  real(dp) :: kpt1(3)
  real(dp) :: delta_str3(2), dstr(2),dk_str(2,2,3)
  real(dp) :: tsec(2)
  real(dp),allocatable :: calc_expibi(:,:),calc_qijb(:,:,:),spkpt(:,:)

  ! *************************************************************************

  DBG_ENTER("COLL")

  call timab(1001,1,tsec)
  call timab(1002,1,tsec)

  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  me=xmpi_comm_rank(spaceComm)

  !save the current value of berryopt
  dtefield%berryopt = dtset%berryopt
  !save the current value of nspinor
  dtefield%nspinor = dtset%nspinor

  !----------------------------------------------------------------------------
  !-------------------- Obtain k-point grid in the full BZ --------------------
  !----------------------------------------------------------------------------

  if(dtset%kptopt==1 .or. dtset%kptopt==2 .or. dtset%kptopt==4)then
     !  Compute the number of k points in the G-space unit cell
     nkptlatt=dtset%kptrlatt(1,1)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,3) &
          &   +dtset%kptrlatt(1,2)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,1) &
          &   +dtset%kptrlatt(1,3)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,2) &
          &   -dtset%kptrlatt(1,2)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,3) &
          &   -dtset%kptrlatt(1,3)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,1) &
          &   -dtset%kptrlatt(1,1)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,2)

     !  Call smpbz to obtain the list of k-point in the full BZ - without symmetry reduction
     option = 0
     brav = 1
     mkpt=nkptlatt*dtset%nshiftk
     ABI_ALLOCATE(spkpt,(3,mkpt))
     call smpbz(1,ab_out,dtset%kptrlatt,mkpt,fnkpt_computed,dtset%nshiftk,option,dtset%shiftk,spkpt)
     dtefield%fnkpt = fnkpt_computed
     ABI_ALLOCATE(dtefield%fkptns,(3,dtefield%fnkpt))
     dtefield%fkptns(:,:)=spkpt(:,1:dtefield%fnkpt)
     ABI_DEALLOCATE(spkpt)
  else if(dtset%kptopt==3.or.dtset%kptopt==0)then
     dtefield%fnkpt=nkpt
     ABI_ALLOCATE(dtefield%fkptns,(3,dtefield%fnkpt))
     dtefield%fkptns(1:3,1:dtefield%fnkpt)=dtset%kpt(1:3,1:dtefield%fnkpt)
     if(dtset%kptopt==0)then
        write(message,'(10a)') ch10,&
             &     ' initberry : WARNING -',ch10,&
             &     '  you have defined manually the k-point grid with kptopt = 0',ch10,&
             &     '  the berry phase calculation works only with a regular k-points grid,',ch10,&
             &     '  abinit doesn''t check if your grid is regular...'
        call wrtout(std_out,message,'PERS')
     end if
  end if

  !call listkk to get mapping from FBZ to IBZ
  rdum=1.0d-5  ! cutoff distance to decide when two k points match
  ABI_ALLOCATE(dtefield%indkk_f2ibz,(dtefield%fnkpt,6))

  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

  !ji: The following may need modification in the future
  !**** no spin-polarization doubling ; allow use of time reversal symmetry ****

  !Here is original call
  !
  !call listkk(rdum,gmet,dtefield%indkk_f2ibz,dtset%kptns,dtefield%fkptns,nkpt,&
  !& dtefield%fnkpt,dtset%nsym,1,dtset%symafm,dtset%symrel,1, spaceComm)

  call timab(1002,2,tsec)
  call timab(1003,1,tsec)

  call listkk(rdum,gmet,dtefield%indkk_f2ibz,dtset%kptns,dtefield%fkptns,nkpt,&
       & dtefield%fnkpt,dtset%nsym,1,dtset%symafm,symrec,1, spaceComm, use_symrec=.True.)

  call timab(1003,2,tsec)
  call timab(1004,1,tsec)

  !Construct i2fbz and f2ibz
  ABI_ALLOCATE(dtefield%i2fbz,(nkpt))
  idum=0
  do ikpt=1,dtefield%fnkpt
     if (dtefield%indkk_f2ibz(ikpt,2)==1 .and. &
          &   dtefield%indkk_f2ibz(ikpt,6) == 0 .and. &
          &   maxval(abs(dtefield%indkk_f2ibz(ikpt,3:5))) == 0 ) then
        dtefield%i2fbz(dtefield%indkk_f2ibz(ikpt,1))=ikpt
        idum=idum+1
     end if
  end do
  if (idum/=nkpt)then
     MSG_ERROR('Found wrong number of k-points in IBZ')
  end if

  !set flags for fields, forces, stresses
  fieldflag = ( (dtset%berryopt== 4) .or. (dtset%berryopt== 6) .or. (dtset%berryopt== 7)  &
       & .or. (dtset%berryopt==14) .or. (dtset%berryopt==16) .or. (dtset%berryopt==17) )
  ! following two flags activates computation of projector gradient contributions to force and
  ! stress in finite field PAW calculations
  calc_epaw3_force = (fieldflag .and. usepaw == 1 .and. dtset%optforces /= 0)
  calc_epaw3_stress = (fieldflag .and. usepaw == 1 .and. dtset%optstress /= 0)



  !----------------------------------------------------------------------------
  !------------- Allocate PAW space if necessary ------------------------------
  !----------------------------------------------------------------------------

  if (usepaw == 1) then

     dtefield%usepaw   = usepaw
     dtefield%natom    = natom
     dtefield%my_natom = mpi_enreg%my_natom

     ABI_ALLOCATE(dtefield%lmn_size,(ntypat))
     ABI_ALLOCATE(dtefield%lmn2_size,(ntypat))
     do itypat = 1, ntypat
        dtefield%lmn_size(itypat) = pawtab(itypat)%lmn_size
        dtefield%lmn2_size(itypat) = pawtab(itypat)%lmn2_size
     end do

     lmn2_size_max = psps%lmnmax*(psps%lmnmax+1)/2
     dtefield%lmn2max = lmn2_size_max

     ! expibi and qijb_kk are NOT parallelized over atoms
     ! this may change in the future (JZwanziger 18 March 2014)
     ABI_ALLOCATE(dtefield%qijb_kk,(2,lmn2_size_max,dtefield%natom,3))
     ABI_ALLOCATE(dtefield%expibi,(2,dtefield%natom,3))
     dtefield%has_expibi = 1
     dtefield%has_qijb = 1

     if ( fieldflag .and. dtefield%has_rij==0) then
        lmn2_size_max = psps%lmnmax*(psps%lmnmax+1)/2
        ABI_ALLOCATE(dtefield%rij,(lmn2_size_max,ntypat,3))
        dtefield%has_rij = 1
     end if

     ! additional F3-type force term for finite electric field with PAW. Same term
     ! might also apply for other displacement-type field calculations, but not sure yet
     ! JZwanziger 4 April 2014
     if ( calc_epaw3_force ) then
        ABI_ALLOCATE(dtefield%epawf3,(dtefield%natom,3,3))
        dtefield%has_epawf3 = 1
     end if
     if ( calc_epaw3_stress ) then
        ABI_ALLOCATE(dtefield%epaws3,(dtefield%natom,3,6))
        dtefield%has_epaws3 = 1
     end if

     ncpgr = 0
     if ( fieldflag .and. dtefield%usecprj == 0) then
        ABI_ALLOCATE(dimlmn,(natom))
        call pawcprj_getdim(dimlmn,natom,nattyp_dum,ntypat,typat,pawtab,'R')
        !    allocate space for cprj at kpts in BZ (IBZ or FBZ)
        ABI_DATATYPE_ALLOCATE(dtefield%cprj,(natom, mband*dtset%nspinor*dtset%nkpt*nsppol))
        !    write(std_out,*) "initberry alloc of cprj ", shape(dtefield%cprj)
        if (calc_epaw3_force .and. .not. calc_epaw3_stress) ncpgr = 3
        if (.not. calc_epaw3_force .and. calc_epaw3_stress) ncpgr = 6
        if (calc_epaw3_force .and. calc_epaw3_stress) ncpgr = 9
        call pawcprj_alloc(dtefield%cprj,ncpgr,dimlmn)
        dtefield%usecprj = 1
        ABI_DEALLOCATE(dimlmn)
     end if

     ABI_ALLOCATE(dtefield%cprjindex,(nkpt,nsppol))
     dtefield%cprjindex(:,:) = 0

     if (dtset%kptopt /= 3) then
        ABI_ALLOCATE(dtefield%atom_indsym,(4,nsym,natom))
        call symatm(dtefield%atom_indsym,natom,nsym,symrec,dtset%tnons,tol8,typat,xred)
        lmax = psps%mpsang - 1
        ABI_ALLOCATE(dtefield%zarot,(2*lmax+1,2*lmax+1,lmax+1,nsym))
        call setsym_ylm(gprimd,lmax,nsym,1,rprimd,symrec,dtefield%zarot)
        dtefield%nsym = nsym
        dtefield%lmax = lmax
        dtefield%lmnmax = psps%lmnmax
     end if

  end if

  !------------------------------------------------------------------------------
  !------------------- Compute variables related to MPI // ----------------------
  !------------------------------------------------------------------------------




  if (nproc==1) then
     dtefield%fmkmem = dtefield%fnkpt
     dtefield%fmkmem_max = dtefield%fnkpt
     dtefield%mkmem_max = nkpt
  else
     dtefield%fmkmem = 0
     do ikpt = 1, dtefield%fnkpt
        ikpti = dtefield%indkk_f2ibz(ikpt,1)
        nband_k = dtset%nband(ikpti)
        if (.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,-1,me))) &
             &     dtefield%fmkmem = dtefield%fmkmem + 1
     end do
     !  Maximum value of mkmem and fmkmem
     call xmpi_max(dtefield%fmkmem,dtefield%fmkmem_max,spaceComm,ierr)
     !  I have to use the dummy variable mkmem_ because
     !  mkmem is declared as intent(in) while the first
     !  argument of xmpi_max must be intent(inout)
     mkmem_ = mkmem
     call xmpi_max(mkmem_,dtefield%mkmem_max,spaceComm,ierr)
  end if

  ABI_ALLOCATE(mpi_enreg%kpt_loc2fbz_sp,(0:nproc-1,1:dtefield%fmkmem_max*nsppol, 1:2))
  ABI_ALLOCATE(mpi_enreg%kpt_loc2ibz_sp,(0:nproc-1,1:dtefield%mkmem_max*nsppol, 1:2))
  ABI_ALLOCATE(mpi_enreg%kptdstrb,(nproc,6,dtefield%fmkmem_max*nsppol*2))
  ABI_ALLOCATE(mpi_enreg%mkmem,(0:nproc-1))
  mpi_enreg%kpt_loc2fbz_sp(:,:,:) = 0
  mpi_enreg%kpt_loc2ibz_sp(:,:,:) = 0
  mpi_enreg%kptdstrb(:,:,:)       = 0
  mpi_enreg%mkmem(:)              = 0

  if (fieldflag) then
     ABI_ALLOCATE(dtefield%cgqindex,(3,6,nkpt*nsppol))
     ABI_ALLOCATE(dtefield%nneigh,(nkpt))
     dtefield%cgqindex(:,:,:) = 0 ; dtefield%nneigh(:) = 0
  end if

  pwind_alloc = mpw*dtefield%fmkmem_max
  ABI_ALLOCATE(pwind,(pwind_alloc,2,3))
  ABI_ALLOCATE(pwnsfac,(2,pwind_alloc))

  !------------------------------------------------------------------------------
  !---------------------- Compute efield_type variables -------------------------
  !------------------------------------------------------------------------------

  !Initialization of efield_type variables
  mult=dtset%useria+1
  dtefield%efield_dot(:) = zero
  dtefield%dkvecs(:,:) = zero
  dtefield%maxnstr = 0    ; dtefield%maxnkstr  = 0
  dtefield%nstr(:) = 0    ; dtefield%nkstr(:) = 0
  ABI_ALLOCATE(dtefield%ikpt_dk,(dtefield%fnkpt,2,3))
  ABI_ALLOCATE(dtefield%cgindex,(nkpt,nsppol))
  ABI_ALLOCATE(dtefield%kgindex,(nkpt))
  ABI_ALLOCATE(dtefield%fkgindex,(dtefield%fnkpt))
  dtefield%ikpt_dk(:,:,:) = 0
  dtefield%cgindex(:,:) = 0
  dtefield%mband_occ = 0
  ABI_ALLOCATE(dtefield%nband_occ,(nsppol))
  dtefield%kgindex(:) = 0
  dtefield%fkgindex(:) = 0

  if (fieldflag) then
     dtset%rfdir(1:3) = 1
  end if


  !Compute spin degeneracy
  if (nsppol == 1 .and. dtset%nspinor == 1) then
     dtefield%sdeg = two
  else if (nsppol == 2 .or. my_nspinor == 2) then
     dtefield%sdeg = one
  end if

  !Compute the number of occupied bands and check that
  !it is the same for each k-point

  index = 0
  do isppol = 1, nsppol
     dtefield%nband_occ(isppol) = 0
     do ikpt = 1, nkpt

        mband_occ_k = 0
        nband_k = dtset%nband(ikpt + (isppol - 1)*nkpt)

        do iband = 1, nband_k
           index = index + 1
           if (abs(occ(index) - dtefield%sdeg) < tol8) mband_occ_k = mband_occ_k + 1
        end do

        if (fieldflag) then
           if (nband_k /= mband_occ_k) then
              write(message,'(a,a,a)')&
                   &         '  In a finite electric field, nband must be equal ',ch10,&
                   &         '  to the number of valence bands.'
              MSG_ERROR(message)
           end if
        end if

        if (ikpt > 1) then
           if (dtefield%nband_occ(isppol) /= mband_occ_k) then
              message = "The number of valence bands is not the same for every k-point of present spin channel"
              MSG_ERROR(message)
           end if
        else
           dtefield%mband_occ         = max(dtefield%mband_occ, mband_occ_k)
           dtefield%nband_occ(isppol) = mband_occ_k
        end if

     end do                ! close loop over ikpt
  end do                ! close loop over isppol

  if (fieldflag) then
     ABI_ALLOCATE(dtefield%smat,(2,dtefield%mband_occ,dtefield%mband_occ,nkpt*nsppol,2,3))

     dtefield%smat(:,:,:,:,:,:) = zero
  end if

  ABI_ALLOCATE(dtefield%sflag,(dtefield%mband_occ,nkpt*nsppol,2,3))
  dtefield%sflag(:,:,:,:) = 0

  !Compute the location of each wavefunction

  icg = 0
  icprj = 0
  !ikg = 0
  do isppol = 1, nsppol
     do ikpt = 1, nkpt

        nband_k = dtset%nband(ikpt + (isppol-1)*nkpt)

        if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle

        dtefield%cgindex(ikpt,isppol) = icg
        npw_k = npwarr(ikpt)
        icg = icg + npw_k*dtefield%nspinor*nband_k

        if (usepaw == 1) then
           dtefield%cprjindex(ikpt,isppol) = icprj
           icprj = icprj + dtefield%nspinor*nband_k
        end if

     end do
  end do

  ikg = 0
  do ikpt = 1, nkpt
     if ((proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,1,me)).and.&
          &   (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,nsppol,me))) cycle

     npw_k = npwarr(ikpt)
     dtefield%kgindex(ikpt) = ikg
     ikg = ikg + npw_k
  end do

  !Need to use dtset%red_efieldbar in the whole code
  !Compute the reciprocal lattice coordinates of the electric field
  if (fieldflag) then

     call  metric(gmetllc,gprimdlc,-1,rmetllc,rprimd,ucvol_local)

     if (dtset%berryopt == 4 .or. dtset%berryopt == 6 .or. dtset%berryopt == 7) then

        do ii=1,3
           dtset%red_efieldbar(ii) = dot_product(dtset%efield(:),rprimd(:,ii))
           dtefield%efield_dot(ii) =  dtset%red_efieldbar(ii)
        end do

        !    dtefield%efield_dot(1) = dot_product(dtset%efield(:),rprimd(:,1))
        !    dtefield%efield_dot(2) = dot_product(dtset%efield(:),rprimd(:,2))
        !    dtefield%efield_dot(3) = dot_product(dtset%efield(:),rprimd(:,3))

        write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
             &     ' initberry: Reduced electric field (ebar)',ch10,&
             &     '  red_efieldbar(1:3) = ',dtset%red_efieldbar(1:3),ch10
        call wrtout(std_out,message,'COLL')

     end if

     if (dtset%berryopt == 6 .or. dtset%berryopt ==7 ) then

        do ii=1,3
           dtset%red_dfield(ii)= (dot_product(dtset%dfield(:),gprimdlc(:,ii)))*ucvol_local/(4.d0*pi)
        end do

        write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
             &     ' initberry: Reduced electric displacement field',ch10,&
             &     '  red_dfield(1:3) = ',dtset%red_dfield(1:3),ch10
        call wrtout(std_out,message,'COLL')

     end if


     if (  dtset%berryopt == 14 ) then
        !    transfer to unreduced electric field.
        do idir=1,3
           dtset%efield(idir)= dot_product(dtset%red_efieldbar(:),gprimdlc(:,idir))
           dtefield%efield_dot(idir) = dtset%red_efieldbar(idir)
           !      dtefield%efield2(idir)=dtset%red_efieldbar(idir)
        end do

        !    dtefield%efield_dot(1) = dtset%red_efieldbar(1)
        !    dtefield%efield_dot(2) = dtset%red_efieldbar(2)
        !    dtefield%efield_dot(3) = dtset%red_efieldbar(3)

        write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
             &     ' initberry: Unreduced electric field (a.u.)',ch10,&
             &     '  efield(1:3) = ',dtset%efield(1:3),ch10
        call wrtout(std_out,message,'COLL')

        write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
             &     ' initberry: Reduced electric field (ebar)',ch10,&
             &     '  red_efieldbar(1:3) = ',dtset%red_efieldbar(1:3),ch10
        call wrtout(std_out,message,'COLL')

     end if


     if ( dtset%berryopt == 16 ) then

        !    to calculate D
        do ii=1,3
           dtset%dfield(ii)  =(4*pi/ucvol_local)*dot_product(dtset%red_dfield(:),rprimd(:,ii))
        end do

        do idir=1,3
           dtset%efield(idir)= (4*pi/ucvol_local)*dot_product(dtset%red_efield(:),rprimd(:,idir))
        end do

        do idir=1,3
           dtset%red_efieldbar(idir)= (4*pi/ucvol_local)*dot_product(dtset%red_efield(:),rmetllc(:,idir))
           dtefield%efield_dot(idir) = dtset%red_efieldbar(idir)
        end do

        !    dtefield%efield_dot(1) = dtset%red_efieldbar(1)
        !    dtefield%efield_dot(2) = dtset%red_efieldbar(2)
        !    dtefield%efield_dot(3) = dtset%red_efieldbar(3)


        write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
             &     ' initberry: Unreduced electric displacement field (a.u.)',ch10,&
             &     '  dfield(1:3) = ',dtset%dfield(1:3),ch10
        call wrtout(std_out,message,'COLL')

        write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
             &     ' initberry: Unreduced electric field (a.u.)',ch10,&
             &     '  efield(1:3) = ',dtset%efield(1:3),ch10
        call wrtout(std_out,message,'COLL')

        write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
             &     ' initberry: Reduced electric field (ebar)',ch10,&
             &     '  red_efieldbar(1:3) = ',dtset%red_efieldbar(1:3),ch10
        call wrtout(std_out,message,'COLL')

     end if

     if ( dtset%berryopt ==17) then

        !    to calculate D

        do idir=1,3
           dtset%efield(idir)= dot_product(dtset%red_efieldbar(:),gprimdlc(:,idir))  ! from ebar
           dtset%dfield(idir)  =(4*pi/ucvol_local)*dot_product(dtset%red_dfield(:),rprimd(:,idir))
           !      dtset%red_efield(idir) = (ucvol_local/(4*pi))*dot_product(dtset%red_efieldbar(:),gmetllc(:,idir))
           dtefield%efield_dot(idir) = dtset%red_efieldbar(idir)
        end do

        write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
             &     ' initberry: Reduced electric field (ebar)',ch10,&
             &     '  red_efieldbar(1:3) = ',dtset%red_efieldbar(1:3),ch10
        call wrtout(std_out,message,'COLL')


        write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
             &     ' initberry: Unreduced electric field (a.u.)',ch10,&
             &     '  efield(1:3) = ',dtset%efield(1:3),ch10
        call wrtout(std_out,message,'COLL')

        write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
             &     ' initberry: Reduced electric displacement field (a.u.)',ch10,&
             &     '  red_dfield(1:3) = ',dtset%red_dfield(1:3),ch10
        call wrtout(std_out,message,'COLL')

        write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
             &     ' initberry: Unreduced electric displacement field (a.u.)',ch10,&
             &     '  dfield(1:3) = ',dtset%dfield(1:3),ch10
        call wrtout(std_out,message,'COLL')


     end if



  end if

  call timab(1004,2,tsec)

  !------------------------------------------------------------------------------
  !---------------------- Compute dk --------------------------------------------
  !------------------------------------------------------------------------------

  call timab(1005,1,tsec)

  do idir = 1, 3

     if (dtset%rfdir(idir) == 1) then

        !    Compute dk(:), the vector between a k-point and its nearest
        !    neighbour along the direction idir

        dk(:) = zero
        dk(idir) = 1._dp   ! 1 mean there is no other k-point un the direction idir
        do ikpt = 2, dtefield%fnkpt
           diffk(:) = abs(dtefield%fkptns(:,ikpt) - dtefield%fkptns(:,1))
           if ((diffk(1) < dk(1)+tol8).and.(diffk(2) < dk(2)+tol8).and.&
                &       (diffk(3) < dk(3)+tol8)) dk(:) = diffk(:)
        end do
        dtefield%dkvecs(:,idir) = dk(:)
        !    DEBUG
        !    write(std_out,*)' initberry : idir, dk', idir, dk
        !    ENDDEBUG

        !    For each k point, find k_prim such that k_prim= k + dk mod(G)
        !    where G is a vector of the reciprocal lattice

        do ikpt = 1, dtefield%fnkpt

           !      First k+dk, then k-dk
           do isign=-1,1,2
              kpt_shifted1=dtefield%fkptns(1,ikpt)- isign*dk(1)
              kpt_shifted2=dtefield%fkptns(2,ikpt)- isign*dk(2)
              kpt_shifted3=dtefield%fkptns(3,ikpt)- isign*dk(3)
              ! Note that this is still a order fnkpt**2 algorithm.
              ! It is possible to implement a order fnkpt algorithm, see listkk.F90.
              do ikpt1 = 1, dtefield%fnkpt
                 diffk1=dtefield%fkptns(1,ikpt1) - kpt_shifted1
                 if(abs(diffk1-nint(diffk1))>tol8)cycle
                 diffk2=dtefield%fkptns(2,ikpt1) - kpt_shifted2
                 if(abs(diffk2-nint(diffk2))>tol8)cycle
                 diffk3=dtefield%fkptns(3,ikpt1) - kpt_shifted3
                 if(abs(diffk3-nint(diffk3))>tol8)cycle
                 dtefield%ikpt_dk(ikpt,(isign+3)/2,idir) = ikpt1
                 exit
              end do   ! ikpt1
           end do     ! isign

           !      OLD CODING
           !      First: k + dk
           !      do ikpt1 = 1, dtefield%fnkpt
           !      diffk(:) = abs(dtefield%fkptns(:,ikpt1) - &
           !      &         dtefield%fkptns(:,ikpt) - dk(:))
           !      if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
           !      dtefield%ikpt_dk(ikpt,1,idir) = ikpt1
           !      exit
           !      end if
           !      end do

           !      Second: k - dk
           !      do ikpt1 = 1, dtefield%fnkpt
           !      diffk(:) = abs(dtefield%fkptns(:,ikpt1) - &
           !      &         dtefield%fkptns(:,ikpt) + dk(:))
           !      if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
           !      dtefield%ikpt_dk(ikpt,2,idir) = ikpt1
           !      exit
           !      end if
           !      end do

        end do     ! ikpt

        !    Find the string length, starting from k point 1
        !    (all strings must have the same number of points)

        nkstr = 1
        ikpt1 = 1
        do ikpt = 1, dtefield%fnkpt
           ikpt1 = dtefield%ikpt_dk(ikpt1,1,idir)
           if (ikpt1 == 1) exit
           nkstr = nkstr + 1
        end do

        !    Check that the string length is a divisor of nkpt
        if(mod(dtefield%fnkpt,nkstr) /= 0) then
           write(message,'(a,i5,a,i7)')&
                &       ' The string length = ',nkstr,&
                &       ', is not a divisor of fnkpt =',dtefield%fnkpt
           MSG_BUG(message)
        end if

        dtefield%nkstr(idir) = nkstr
        dtefield%nstr(idir)  = dtefield%fnkpt/nkstr

     end if      ! dtset%rfdir(idir) == 1

     write(message,'(a,i1,a,i3,a,i6)')&
          &   '  initberry: for direction ',idir,', nkstr = ',dtefield%nkstr(idir),&
          &   ', nstr = ',dtefield%nstr(idir)
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')

  end do     ! close loop over idir

  call timab(1005,2,tsec)
  call timab(1006,1,tsec)

  dtefield%maxnstr  = maxval(dtefield%nstr(:))
  dtefield%maxnkstr = maxval(dtefield%nkstr(:))
  ABI_ALLOCATE(dtefield%idxkstr,(dtefield%maxnkstr,dtefield%maxnstr,3))
  dtefield%idxkstr(:,:,:) = 0

  !for the geometry of the string space :
  ABI_ALLOCATE(dtefield%coord_str,(2,dtefield%maxnstr,3))
  ABI_ALLOCATE(dtefield%str_neigh,(-2:2,dtefield%maxnstr,3))
  ABI_ALLOCATE(dtefield%strg_neigh,(-2:2,dtefield%maxnstr,2,3))
  dtefield%coord_str(:,:,:) = 0.d0
  dtefield%str_neigh(:,:,:)=0
  dtefield%strg_neigh(:,:,:,:)=0
  dtefield%gmet_str(:,:,:)=0.d0

  !------------------------------------------------------------------------------
  !---------------------- Build the strings -------------------------------------
  !------------------------------------------------------------------------------

  ABI_ALLOCATE(kpt_mark,(dtefield%fnkpt))
  do idir = 1, 3

     if (dtset%rfdir(idir) == 1) then

        iunmark = 1
        kpt_mark(:) = 0
        do istr = 1, dtefield%nstr(idir)

           do while(kpt_mark(iunmark) /= 0)
              iunmark = iunmark + 1
           end do
           dtefield%idxkstr(1,istr,idir) = iunmark
           kpt_mark(iunmark) = 1
           do ikstr = 2, dtefield%nkstr(idir)
              ikpt1 = dtefield%idxkstr(ikstr-1,istr,idir)
              ikpt2 = dtefield%ikpt_dk(ikpt1,1,idir)
              dtefield%idxkstr(ikstr,istr,idir) = ikpt2
              kpt_mark(ikpt2) = 1
           end do

        end do    ! istr

        !    compute distance between strings
        !    compute the metric matrix of the strings space in the direction idir
        do ii = 1,3
           do jj = 1,3
              if (ii<idir.and.jj<idir) dtefield%gmet_str(ii  ,jj  ,idir) = &
                   &         gmet(ii,jj) - gmet(ii,idir)*gmet(jj,idir)/gmet(idir,idir)
              if (ii<idir.and.jj>idir) dtefield%gmet_str(ii  ,jj-1,idir) = &
                   &         gmet(ii,jj) - gmet(ii,idir)*gmet(jj,idir)/gmet(idir,idir)
              if (ii>idir.and.jj<idir) dtefield%gmet_str(ii-1,jj  ,idir) = &
                   &         gmet(ii,jj) - gmet(ii,idir)*gmet(jj,idir)/gmet(idir,idir)
              if (ii>idir.and.jj>idir) dtefield%gmet_str(ii-1,jj-1,idir) = &
                   &         gmet(ii,jj) - gmet(ii,idir)*gmet(jj,idir)/gmet(idir,idir)
           end do
        end do
        !    DEBUG
        !    write(std_out,*)'gmet'
        !    do ii=1,3
        !    write(std_out,*)gmet(ii,:)
        !    end do
        !    write(std_out,*)'gmet_str'
        !    do ii=1,2
        !    write(std_out,*)dtefield%gmet_str(ii,:,idir)
        !    end do
        !    ENDDEBUG
        do istr = 1, dtefield%nstr(idir)
           do ii = 1,3
              if (ii<idir) dtefield%coord_str(ii,istr,idir)=dtefield%fkptns(ii,dtefield%idxkstr(1,istr,idir))
              if (ii>idir) dtefield%coord_str(ii-1,istr,idir)=dtefield%fkptns(ii,dtefield%idxkstr(1,istr,idir))
           end do
        end do

        !    the following is very similar to getshell
        dist_ = 0._dp
        do ii = 1,2
           dist_ = dist_ + dtefield%gmet_str(ii,ii,idir)
        end do
        max_dist = 2._dp * dist_ * 2._dp

        dk_str(:,:,idir) = 0._dp
        last_dist = 0._dp
        !    ishell = 0
        !    dtefield%str_neigh(:,:,:) = 0
        dk_flag = 0
        do while (dk_flag /= 2)
           !      Advance shell counter
           !      ishell = ishell + 1

           !      Search the smallest distance between two strings
           dist = max_dist
           do istr = 1,dtefield%nstr(idir)
              delta_str3(:) = dtefield%coord_str(:,1,idir) - dtefield%coord_str(:,istr,idir)
              do coord1 = -1,1  !two loop to search also on the border of the BZ
                 do coord2 = -1,1
                    dist_ = 0._dp
                    dstr(:) = delta_str3(:) - nint(delta_str3(:))
                    dstr(1) = dstr(1) + real(coord1,dp)
                    dstr(2) = dstr(2) + real(coord2,dp)
                    do ii = 1,2
                       do jj = 1,2
                          dist_ = dist_ + dstr(ii)*dtefield%gmet_str(ii,jj,idir)*dstr(jj)
                       end do
                    end do
                    if ((dist_ < dist).and.(dist_ - last_dist > tol8)) then
                       dist = dist_
                    end if
                 end do
              end do
           end do

           last_dist = dist

           !      search the connecting vectors for that distance
           do istr = 1,dtefield%nstr(idir)
              delta_str3(:) = dtefield%coord_str(:,istr,idir) - dtefield%coord_str(:,1,idir)
              do coord1 = -1,1
                 do coord2 = -1,1
                    dist_ = 0._dp
                    dstr(:) = delta_str3(:) - nint(delta_str3(:))
                    dstr(1) = dstr(1) + real(coord1,dp)
                    dstr(2) = dstr(2) + real(coord2,dp)
                    do ii = 1,2
                       do jj = 1,2
                          dist_ = dist_ + dstr(ii)*dtefield%gmet_str(ii,jj,idir)*dstr(jj)
                       end do
                    end do
                    if (abs(dist_ - dist) < tol8) then
                       if (dk_flag == 0) then
                          dk_str(:,1,idir) = dstr(:)
                          dk_flag = 1
                          !                DEBUG
                          !                write(std_out,'(a,i4,2e15.4)')'1st connect', istr, dstr
                          !                ENDDEBUG
                       elseif (dk_str(1,1,idir)*dstr(2)-dk_str(2,1,idir)*dstr(1) > tol8) then
                          dk_str(:,2,idir) = dstr(:)
                          dk_flag = 2
                          !                DEBUG
                          !                write(std_out,'(a,i4,2e15.4)')'2nd connect', istr, dstr
                          !                ENDDEBUG
                          exit
                       end if
                    end if
                 end do
                 if (dk_flag == 2) exit
              end do
              if (dk_flag == 2) exit
           end do

        end do ! do while

        !    search the two neighbours for each string
        do istr = 1,dtefield%nstr(idir)
           dtefield%str_neigh(0,istr,idir) = istr
           dtefield%strg_neigh(0,istr,:,idir) = 0
           do jstr = 1,dtefield%nstr(idir)
              delta_str3(:) = dtefield%coord_str(:,jstr,idir) - dtefield%coord_str(:,istr,idir)
              do coord1 = -1,1
                 do coord2 = -1,1
                    dist_ = 0._dp
                    dstr(:) = delta_str3(:) - nint(delta_str3(:))
                    dstr(1) = dstr(1) + real(coord1,dp)
                    dstr(2) = dstr(2) + real(coord2,dp)
                    do ii = 1,2
                       if (sum(abs(dstr(:)-dk_str(:,ii,idir)))<tol8) then
                          dtefield%str_neigh(ii,istr,idir) = jstr
                          dtefield%strg_neigh(ii,istr,1,idir) = coord1
                          dtefield%strg_neigh(ii,istr,2,idir) = coord2
                       elseif (sum(abs(dstr(:)+dk_str(:,ii,idir)))<tol8) then
                          dtefield%str_neigh(-ii,istr,idir) = jstr
                          dtefield%strg_neigh(-ii,istr,1,idir) = coord1
                          dtefield%strg_neigh(-ii,istr,2,idir) = coord2
                       end if
                    end do
                 end do
              end do
           end do
        end do

        !    DEBUG
        !    write(std_out,'(a,e15.4,e15.4,e15.4,e15.4)')'dk_str',dk_str(1,1,idir),dk_str(2,1,idir),dk_str(1,2,idir),dk_str(2,2,idir)
        !    write(std_out,*)'istr, neigh1, strg(1,:), neigh2, strg(2,:),neigh-1, strg(-1,:), neigh-2, strg(-2,:)'
        !    do istr=1,dtefield%nstr(idir)
        !    write(std_out,'(13i4)')istr, &
        !    &       dtefield%str_neigh(1,istr,idir), dtefield%strg_neigh(1,istr,:,idir),&
        !    &       dtefield%str_neigh(2,istr,idir), dtefield%strg_neigh(2,istr,:,idir),&
        !    &       dtefield%str_neigh(-1,istr,idir), dtefield%strg_neigh(-1,istr,:,idir),&
        !    &       dtefield%str_neigh(-2,istr,idir), dtefield%strg_neigh(-2,istr,:,idir)
        !    end do
        !    ENDDEBUG


     end if         ! rfdir(idir) == 1

  end do           ! close loop over idir

  ABI_DEALLOCATE(kpt_mark)

  call timab(1006,2,tsec)
  call timab(1007,1,tsec)

  !------------------------------------------------------------------------------
  !------------ Compute PAW on-site terms if necessary --------------------------
  !------------------------------------------------------------------------------

  if (usepaw == 1 .and. dtefield%has_expibi == 1) then
     ABI_ALLOCATE(calc_expibi,(2,natom))
     do idir = 1, 3
        dk = dtefield%dkvecs(1:3,idir)
        calc_expibi = zero
        call expibi(calc_expibi,dk,natom,xred)
        dtefield%expibi(1:2,1:natom,idir) = calc_expibi
     end do
     !   call expibi(dtefield%expibi,dtefield%dkvecs,natom,xred)
     dtefield%has_expibi = 2
     ABI_DEALLOCATE(calc_expibi)
  end if

  if (usepaw == 1 .and. dtefield%has_qijb == 1) then
     ABI_ALLOCATE(calc_qijb,(2,dtefield%lmn2max,natom))

     do idir = 1, 3
        dk = dtefield%dkvecs(1:3,idir)
        calc_qijb = zero
        call qijb_kk(calc_qijb,dk,dtefield%expibi(1:2,1:natom,idir),&
             &     gprimd,dtefield%lmn2max,natom,ntypat,pawang,pawrad,pawtab,typat)
        dtefield%qijb_kk(1:2,1:dtefield%lmn2max,1:natom,idir) = calc_qijb
        !    call qijb_kk(dtefield%qijb_kk,dtefield%dkvecs,dtefield%expibi,&
        ! &   gprimd,dtefield%lmn2max,natom,ntypat,pawang,pawrad,pawtab,typat)
     end do
     dtefield%has_qijb = 2
     ABI_DEALLOCATE(calc_qijb)
  end if

  if (usepaw == 1 .and. dtefield%has_rij == 1) then
     c1=sqrt(four_pi/three)
     do itypat = 1, ntypat
        do klmn = 1, pawtab(itypat)%lmn2_size
           dtefield%rij(klmn,itypat,1) = c1*pawtab(itypat)%qijl(4,klmn) ! S_{1,1} ~ x
           dtefield%rij(klmn,itypat,2) = c1*pawtab(itypat)%qijl(2,klmn) ! S_{1,-1} ~ y
           dtefield%rij(klmn,itypat,3) = c1*pawtab(itypat)%qijl(3,klmn) ! S_{1,0} ~ z
        end do ! end loop over klmn
     end do ! end loop over itypat
     dtefield%has_rij = 2
  end if !

  call timab(1007,2,tsec)
  call timab(1008,1,tsec)

  !------------------------------------------------------------------------------
  !------------ Build the array pwind that is needed to compute the -------------
  !------------ overlap matrices at k +- dk                         -------------
  !------------------------------------------------------------------------------

  ecut_eff = dtset%ecut*(dtset%dilatmx)**2
  exchn2n3d = 0 ; istwf_k = 1 ; ikg1 = 0
  pwind(:,:,:) = 0
  pwnsfac(1,:) = 1.0_dp
  pwnsfac(2,:) = 0.0_dp
  ABI_ALLOCATE(kg1_k,(3,mpw))

  ipwnsfac = 0

  do idir = 1, 3

     if (dtset%rfdir(idir) == 1) then

        dk(:) = dtefield%dkvecs(:,idir)

        do ifor = 1, 2

           if (ifor == 2) dk(:) = -1._dp*dk(:)

           !      Build pwind and kgindex
           !      NOTE: The array kgindex is important for parallel execution.
           !      In case nsppol = 2, it may happent that a particular processor
           !      treats k-points at different spin polarizations.
           !      In this case, it is not possible to address the elements of
           !      pwind correctly without making use of the kgindex array.

           ikg = 0 ; ikpt_loc = 0 ; isppol = 1
           do ikpt = 1, dtefield%fnkpt

              ikpti = dtefield%indkk_f2ibz(ikpt,1)
              nband_k = dtset%nband(ikpti)
              ikpt1f = dtefield%ikpt_dk(ikpt,ifor,idir)
              ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

              if ((proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,1,me)).and.&
                   &         (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,nsppol,me))) cycle

              ikpt_loc = ikpt_loc + 1

              !        Build basis sphere of plane waves for the nearest neighbour of
              !        the k-point (important for MPI //)

              kg1_k(:,:) = 0
              kpt1(:) = dtset%kptns(:,ikpt1i)
              call kpgsph(ecut_eff,exchn2n3d,gmet,ikg1,ikpt,istwf_k,kg1_k,kpt1,&
                   &         1,mpi_enreg,mpw,npw_k1)
              me_g0=mpi_enreg%me_g0


              !        ji: fkgindex is defined here !
              dtefield%fkgindex(ikpt) = ikg

              !
              !        Deal with symmetry transformations
              !

              !        bra k-point k(b) and IBZ k-point kIBZ(b) related by
              !        k(b) = alpha(b) S(b)^t kIBZ(b) + G(b)
              !        where alpha(b), S(b) and G(b) are given by indkk_f2ibz
              !
              !        For the ket k-point:
              !        k(k) = alpha(k) S(k)^t kIBZ(k) + G(k) - GBZ(k)
              !        where GBZ(k) takes k(k) to the BZ
              !

              isym  = dtefield%indkk_f2ibz(ikpt,2)
              isym1 = dtefield%indkk_f2ibz(ikpt1f,2)

              !        Construct transformed G vector that enters the matching condition:
              !        alpha(k) S(k)^{t,-1} ( -G(b) - GBZ(k) + G(k) )

              dg(:) = -dtefield%indkk_f2ibz(ikpt,3:5) &
                   &         -nint(-dtefield%fkptns(:,ikpt) - dk(:) - tol10 + &
                   &         dtefield%fkptns(:,ikpt1f)) &
                   &         +dtefield%indkk_f2ibz(ikpt1f,3:5)

              !        old code
              !        iadum(:)=0
              !        do idum=1,3
              !        iadum(:)=iadum(:)+ symrec(:,idum,isym1)*dg(idum)
              !        end do

              !        new code
              iadum(:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),dg(:))

              dg(:) = iadum(:)

              if ( dtefield%indkk_f2ibz(ikpt1f,6) == 1 ) dg(:) = -dg(:)

              !        Construct S(k)^{t,-1} S(b)^{t}

              dum33(:,:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),symrec(:,:,isym))

              !        Construct alpha(k) alpha(b)

              if (dtefield%indkk_f2ibz(ikpt,6) == dtefield%indkk_f2ibz(ikpt1f,6)) then
                 itrs=0
              else
                 itrs=1
              end if


              npw_k  = npwarr(ikpti)
              !        npw_k1 = npwarr(ikpt1i)

              !        loop over bra G vectors
              do ipw = 1, npw_k

                 !          NOTE: the bra G vector is taken for the sym-related IBZ k point,
                 !          not for the FBZ k point
                 iadum(:) = kg(:,dtefield%kgindex(ikpti) + ipw)

                 !          Store non-symmorphic operation phase factor exp[i2\pi \alpha G \cdot t]

                 if ( ipwnsfac == 0 ) then
                    !            old code
                    rdum=0.0_dp
                    do idum=1,3
                       rdum=rdum+dble(iadum(idum))*dtset%tnons(idum,isym)
                    end do
                    rdum=two_pi*rdum
                    if ( dtefield%indkk_f2ibz(ikpt,6) == 1 ) rdum=-rdum
                    pwnsfac(1,ikg+ipw) = cos(rdum)
                    pwnsfac(2,ikg+ipw) = sin(rdum)
                    !
                    !            new code
                    !            rdum = DOT_PRODUCT(dble(iadum(:)),dtset%tnons(:,isym))
                    !            rdum= two_pi*rdum
                    !            if ( dtefield%indkk_f2ibz(ikpt,6) == 1 ) rdum=-rdum
                    !            pwnsfac(1,ikg+ipw) = cos(rdum)
                    !            pwnsfac(2,ikg+ipw) = sin(rdum)

                 end if

                 !          to determine r.l.v. matchings, we transformed the bra vector
                 !          Rotation
                 iadum1(:)=0
                 do idum1=1,3
                    iadum1(:)=iadum1(:)+dum33(:,idum1)*iadum(idum1)
                 end do
                 iadum(:)=iadum1(:)
                 !          Time reversal
                 if (itrs==1) iadum(:)=-iadum(:)
                 !          Translation
                 iadum(:) = iadum(:) + dg(:)

                 do jpw = 1, npw_k1
                    iadum1(1:3) = kg1_k(1:3,jpw)
                    if ( (iadum(1) == iadum1(1)).and. &
                         &             (iadum(2) == iadum1(2)).and. &
                         &             (iadum(3) == iadum1(3)) ) then
                       pwind(ikg + ipw,ifor,idir) = jpw
                       !              write(std_out,'(a,2x,3i4,2x,i4)') 'Found !:',iadum1(:),jpw
                       exit
                    end if
                 end do
              end do

              ikg  = ikg + npw_k

           end do    ! close loop over ikpt

           ipwnsfac = 1

        end do    ! close loop over ifor

     end if      ! rfdir(idir) == 1

  end do        ! close loop over idir


  call timab(1008,2,tsec)
  call timab(1009,1,tsec)

  !Build mpi_enreg%kptdstrb
  !array required to communicate the WFs between cpus in berryphase_new.f
  !(MPI // over k-points)
  if (nproc>1) then
     do idir = 1, 3
        if (dtset%rfdir(idir) == 1) then
           do ifor = 1, 2

              ikpt_loc = 0
              do isppol = 1, nsppol

                 do ikpt = 1, dtefield%fnkpt

                    ikpti = dtefield%indkk_f2ibz(ikpt,1)
                    nband_k = dtset%nband(ikpti)
                    ikpt1f = dtefield%ikpt_dk(ikpt,ifor,idir)
                    ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

                    if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,isppol,me)) cycle

                    ikpt_loc = ikpt_loc + 1
                    mpi_enreg%kptdstrb(me + 1,ifor+2*(idir-1),ikpt_loc) = &
                         &             ikpt1i + (isppol - 1)*nkpt

                    mpi_enreg%kptdstrb(me+1,ifor+2*(idir-1),&
                         &             ikpt_loc+dtefield%fmkmem_max*nsppol) = &
                         &             ikpt1f + (isppol - 1)*dtefield%fnkpt

                 end do   ! ikpt
              end do     ! isppol
           end do       ! ifor
        end if         ! dtset%rfdir(idir) == 1
     end do           ! idir
  end if             ! nproc>1

  !build mpi_enreg%kpt_loc2fbz_sp
  ikpt_loc = 0
  do isppol = 1, nsppol
     do ikpt = 1, dtefield%fnkpt

        ikpti = dtefield%indkk_f2ibz(ikpt,1)
        nband_k = dtset%nband(ikpti)

        if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,isppol,me)) cycle

        ikpt_loc = ikpt_loc + 1

        mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc, 1) = ikpt
        mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc, 2) = isppol

     end do
  end do


  !parallel case only :
  !build mpi_enreg%kpt_loc2ibz_sp, dtefield%cgqindex and dtefield%nneigh
  if ((fieldflag).and.(nproc>1)) then
     ikpt_loc = 0
     do isppol = 1, nsppol
        do ikpt = 1, nkpt

           ikptf = dtefield%i2fbz(ikpt)
           nband_k = dtset%nband(ikpti)

           neigh(:) = 0 ; icg = 0 ; ikg = 0 ; flag_kpt = 0; icprj = 0
           do idir=1, 3

              !        skip idir values for which efield_dot(idir) = 0
              if (abs(dtefield%efield_dot(idir)) < tol12) cycle

              do ifor = 1, 2

                 flag = 0

                 ikpt1f = dtefield%ikpt_dk(ikptf,ifor,idir)
                 ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

                 dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ikg
                 ikg = ikg + npwarr(ikpt1i)

                 !          check if this neighbour is also a previous neighbour
                 do ineigh = 1, (ifor+2*(idir-1))
                    if (neigh(ineigh) == ikpt1i) then
                       flag = 1
                       dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ineigh
                       dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
                            &               dtefield%cgqindex(2,ineigh,ikpt+(isppol-1)*nkpt)
                       exit
                    end if
                 end do
                 !          create the cgqindex of the neighbour if necessary
                 if (flag == 0) then
                    neigh(ifor+2*(idir-1)) = ikpt1i
                    dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
                         &             ifor+2*(idir-1)
                    dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = icg
                    if (isppol == 1) dtefield%nneigh(ikpt) = dtefield%nneigh(ikpt) + 1
                    icg = icg + npwarr(ikpt1i)*dtefield%nspinor*nband_k
                 end if
              end do !ifor
           end do !idir

           if (.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me))) then
              !        ikpt is one of my kpt_loc
              ikpt_loc = ikpt_loc + 1
              mpi_enreg%kpt_loc2ibz_sp(me, ikpt_loc, 1) = ikpt
              mpi_enreg%kpt_loc2ibz_sp(me, ikpt_loc, 2) = isppol
           end if

        end do !ikpt
     end do !isppol
  end if !nproc>1

  !should be temporary
  !unassigned mpi_enreg%kpt_loc2fbz_sp are empty ; inform other cpu (there are better ways...)
  mpi_enreg%mkmem(me) = mkmem
  !do ii=ikpt_loc+1,dtefield%fmkmem_max
  !mpi_enreg%kpt_loc2fbz_sp(me, ii, 1) = -1
  !end do


  !(same as mpi_enreg%kptdstrb but for k-points in the iBZ),
  !dtefield%cgqindex and dtefield%nneigh

  if ((fieldflag).and.(nproc>1)) then

     ikpt_loc = 1
     do isppol = 1, nsppol
        do ikpt = 1, nkpt

           nband_k = dtset%nband(ikpt)
           ikptf = dtefield%i2fbz(ikpt)

           neigh(:) = 0 ; icg = 0 ; ikg = 0 ; flag_kpt = 0; icprj = 0
           do idir = 1, 3

              !        Skip idir values for which efield_dot(idir) = 0
              if (abs(dtefield%efield_dot(idir)) < tol12 .and. (fieldflag)) cycle

              do ifor = 1, 2

                 ikpt1f = dtefield%ikpt_dk(ikptf,ifor,idir)
                 ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

                 !          dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ikg
                 ikg = ikg + npwarr(ikpt1i)

                 flag = 0
                 do ineigh = 1, (ifor+2*(idir-1))
                    if (neigh(ineigh) == ikpt1i) then
                       flag = 1
                       !              dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ineigh
                       !              dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
                       !              &               dtefield%cgqindex(2,ineigh,ikpt+(isppol-1)*nkpt)
                       exit
                    end if
                 end do
                 if (flag == 0) then
                    !            neigh(ifor+2*(idir-1)) = ikpt1i
                    !            dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
                    !            &             ifor+2*(idir-1)
                    !            dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = icg
                    !            if (isppol == 1) dtefield%nneigh(ikpt) = dtefield%nneigh(ikpt) + 1
                    !            icg = icg + npwarr(ikpt1i)*dtset%nspinor*nband_k
                 end if

                 if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle

                 flag_kpt = 1

                 !          MVeithen: the if condition allows to avoid that the same wavefunction
                 !          is send several times to a particular cpu

              end do    ! ifor
           end do    ! idir

           if (flag_kpt == 1) ikpt_loc = ikpt_loc + 1

        end do    ! ikpt
     end do    ! isppol

  end if   ! fieldflag

  call xmpi_sum(mpi_enreg%kptdstrb,spaceComm,ierr)
  call xmpi_sum(mpi_enreg%kpt_loc2fbz_sp,spaceComm,ierr)
  if (fieldflag) then
     call xmpi_sum(mpi_enreg%kpt_loc2ibz_sp,spaceComm,ierr)
     call xmpi_sum(mpi_enreg%mkmem,spaceComm,ierr)
  end if

  !------------------------------------------------------------------------------
  !------------------------ Estimate critical field -----------------------------
  !------------------------------------------------------------------------------

  !Compute the minimal value of the bandgap required to be below
  !the critical field as defined by the relation
  !| E_i*a_i | < E_g/n_i

  if (fieldflag) then

     do idir = 1, 3
        !    eg_dir(idir) = abs(dtefield%efield_dot(idir))*dtefield%nkstr(idir)
        eg_dir(idir) = abs(dtset%red_efieldbar(idir))*dtefield%nkstr(idir)
     end do


     eg = maxval(eg_dir)
     eg_ev = eg*Ha_eV

     if (dtset%optcell ==0 .and. (dtset%berryopt == 4 .or. dtset%berryopt == 14)) then
        write(message,'(a,a,a,a,a,a,a,a,f7.2,a,a)')ch10,&
             &     ' initberry: COMMENT - ',ch10,&
             &     '  As a rough estimate,',ch10,&
             &     '  to be below the critical field, the bandgap of your system',ch10,&
             &     '  should be larger than ',eg_ev,' eV.',ch10
        call wrtout(ab_out,message,'COLL')
        call wrtout(std_out,message,'COLL')

     else

        write(message,'(a,a,a,a,a,a,a)') ch10,&
             &     ' initberry: COMMENT - ',ch10,&
             &     '  The estimation of critical electric field should be checked after calculation.',ch10,&
             &     '  It is printed out just after total energy.' ,ch10

        call wrtout(ab_out,message,'COLL')
        call wrtout(std_out,message,'COLL')

     end if

  end if

  ABI_DEALLOCATE(kg1_k)

  call timab(1009,2,tsec)
  call timab(1001,2,tsec)

  DBG_EXIT("COLL")

end subroutine initberry
!!***


end module m_berryphase_new
!!***
