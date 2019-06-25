!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_mlwfovlp_qp
!! NAME
!!  m_mlwfovlp_qp
!!
!! FUNCTION
!!  Interpolate GW corrections with Wannier functions
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (DRH)
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

module m_mlwfovlp_qp

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wannier90
 use m_errors
 use m_abicore
 use m_xmpi
 use m_hdr

 use m_mpinfo,         only : destroy_mpi_enreg, initmpi_seq
 use m_pawtab,         only : pawtab_type
 use m_pawcprj,        only : pawcprj_type, paw_overlap, pawcprj_getdim, pawcprj_alloc, pawcprj_free
 use m_numeric_tools,  only : isordered
 use m_geometry,       only : metric
 use m_crystal,        only : crystal_t
 use m_kpts,           only : listkk
 use m_bz_mesh,        only : kmesh_t, kmesh_init, kmesh_free
 use m_ebands,         only : ebands_init, ebands_free
 use m_qparticles,     only : rdqps, rdgw
 use m_sort,           only : sort_dp

 implicit none

 private
!!***

 public :: mlwfovlp_qp
!!***

contains
!!***

!!****f* m_mlwfovlp_qp/mlwfovlp_qp
!! NAME
!! mlwfovlp_qp
!!
!! FUNCTION
!! Routine which computes replaces LDA wave functions and eigenvalues with
!! GW quasiparticle ones using previously computed qp wave functions in
!! LDA bloch function representation for Wannier code (www.wannier.org f90 version).
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  dtfil <type(datafiles_type)>=variables related to files
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mkmem =number of k points treated by this node.
!!  mpw=maximum dimensioned size of npw.
!!  natom=number of atoms in cell.
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  Hdr<Hdr_type>=The m_mlwfovlp_qp header.
!!  MPI_enreg=information about MPI parallelization
!!  Cprj_BZ(natom,nspinor*mband*mkmem*nsppol)= <p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=planewave coefficients of wavefunctions
!!   replaced by quasiparticle wavefunctions
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues replaced by qp eigenvalues(hartree)
!!
!! NOTES
!!  Number of bands for wannier calculation must be identical to number used
!!   for gw calculation.  Bands not wanted for wannier calculation must be
!!   excluded in exclude_band statement in wannier90.win file.
!!  Full plane-wave basis for LDA wavefunctions must be used in GW calculation,
!!   or inaccuracies may result.
!!  This is at best a beta version of this code, with little consistency
!!   checking, so the user must be very careful or the results may be invalid.
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      crystal_free,crystal_from_hdr,destroy_mpi_enreg,ebands_free,ebands_init
!!      initmpi_seq,kmesh_free,kmesh_init,listkk,metric,pawcprj_getdim,rdgw
!!      rdqps,sort_dp,update_cprj,wrtout,zgemm
!!
!! SOURCE

subroutine mlwfovlp_qp(cg,Cprj_BZ,dtset,dtfil,eigen,mband,mcg,mcprj,mkmem,mpw,natom,&
& nkpt,npwarr,nspden,nsppol,ntypat,Hdr,Pawtab,rprimd,MPI_enreg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mcprj,mkmem,mpw,nkpt,nspden,natom,ntypat
 integer,intent(in) :: nsppol
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(Hdr_type),intent(in) :: Hdr
 type(MPI_type),intent(in) :: MPI_enreg
 type(pawcprj_type),target,intent(inout) :: Cprj_BZ(natom,mcprj)
 type(Pawtab_type),intent(in) :: Pawtab(ntypat*Dtset%usepaw)
!arrays
 integer,intent(in) :: npwarr(nkpt)
 real(dp),intent(inout) :: cg(2,mcg)
 real(dp),intent(inout) :: eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer,parameter :: from_QPS_FILE=1,from_GW_FILE=2
 integer :: sppoldbl,timrev,bantot_ibz,ikibz,ikbz,dimrho
 integer :: iband,icg,icg_shift,ii,ipw,isppol,my_nspinor,nband_k,ord_iband
 integer :: nfftot,ikpt,irzkpt,npw_k,ikg
 integer :: nscf,nbsc,itimrev,band_index,nkibz,nkbz
 integer :: gw_timrev,input !,jb_idx,ib_idx,ijpack, jband,
 integer :: nprocs,ios
 real(dp) :: TOL_SORT=tol12
 real(dp) :: dksqmax,ucvol !ortho_err,
 logical :: ltest,qpenek_is_ordered,g0w0_exists
 character(len=500) :: msg
 character(len=fnlen) :: gw_fname
 type(ebands_t) :: QP_bst
 type(crystal_t)  :: Cryst
 type(kmesh_t) :: Kibz_mesh
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: indkk(nkpt,6),my_ngfft(18)
 integer,allocatable :: npwarr_ibz(:),nband_ibz(:),ibz2bz(:,:),istwfk_ibz(:)
 integer,allocatable :: dimlmn(:),iord(:),nattyp_dum(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3) !,paw_ovlp(2)
 real(dp),allocatable :: qp_rhor(:,:),sorted_qpene(:)
 real(dp),allocatable :: kibz(:,:),wtk_ibz(:)
 real(dp),allocatable :: doccde_ibz(:),occfact_ibz(:),eigen_ibz(:)
 real(dp),allocatable ::  igwene(:,:,:)
 complex(dpc),allocatable :: m_lda_to_qp(:,:,:,:),m_lda_to_qp_BZ(:,:,:,:) !,ortho(:)
 complex(dpc),allocatable :: m_tmp(:,:),cg_k(:,:),cg_qpk(:,:)
 type(Pawrhoij_type),allocatable :: prev_Pawrhoij(:)
 !type(pawcprj_type),pointer :: Cp1(:,:),Cp2(:,:)

!************************************************************************

 ABI_UNUSED(mkmem)

 DBG_ENTER("COLL")

 write(msg,'(17a)')ch10,&
&  ' mlwfovlp_qp: WARNING',ch10,&
&  '  The input *_WFK file of LDA wavefunctions to be  converted',ch10,&
&  '  to GW quasiparticle wavefunctions MUST have been written in',ch10,&
&  '  the run that produced the GW *_KSS file using kssform 3,',ch10,&
&  '  the ONLY value of kssform permitted for GW Wannier functions.',ch10,&
&  '  Otherwise, the *_QPS file needed here will be inconsistent,',ch10,&
&  '  and the output quasiparticle wavefunctions will be garbage.',ch10,&
&  '  No internal check that can verify this is presently possible.',ch10
 call wrtout(std_out,msg,'COLL')

 ! === Some features are not implemented yet ===
 ABI_CHECK(Dtset%nspinor==1,'nspinor==2 not implemented')
 ABI_CHECK(Dtset%nsppol==1,'nsppol==2 not implemented, check wannier90')
 ltest=ALL(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol)==Dtset%nband(1))
 ABI_CHECK(ltest,'nband(:) should be constant')
 !
 ! MPI initialization
 nprocs=MPI_enreg%nproc_cell

 if (nprocs/=1) then
   MSG_ERROR("mlwfovlp_qp not programmed for parallel execution")
 end if

 ! Compute reciprocal space metric gmet for unit cell of disk wf
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 ! Compute k points from gw irreducible set equivalent to full-zone wannier set
 sppoldbl=1 ; timrev=1 ; my_nspinor=max(1,Dtset%nspinor/MPI_enreg%nproc_spinor)
 call listkk(dksqmax,gmet,indkk,dtset%kptgw,dtset%kpt,dtset%nkptgw,nkpt,&
&  dtset%nsym,sppoldbl,dtset%symafm,dtset%symrel,timrev,xmpi_comm_self)

 if (dksqmax>tol8) then
   write(msg,'(5a)')&
&    'Set of GW irreducible-zone kptgw in input file is inconsistent',ch10,&
&    'with full-zone set being used for wannier90 setup.',ch10,&
&    'Action: correct input data'
   MSG_ERROR(msg)
 end if
 !
 ! === Initialize object defining the Band strucuture ===
 ! * Initialize with KS results using IBZ indexing.
 ! * After rdqps, QP_bst will contain the QP amplitudes.
 nkibz=Dtset%nkptgw
 ABI_MALLOC(kibz,(3,nkibz))
 ABI_MALLOC(wtk_ibz,(nkibz))
 kibz=Dtset%kptgw(:,1:Dtset%nkptgw)

 ! MG: This part is needed to get the IBZ weight that will be reported
 ! on ab_out thus we should be consistent. Ideally Cryst should be
 ! one of the basic abinit objects and it should be passed to this routine.

 gw_timrev=1; if (timrev==1) gw_timrev=2 !different conventions are used in GW and abinit!!

 cryst = hdr_get_crystal(Hdr, gw_timrev)
 call kmesh_init(Kibz_mesh,Cryst,nkibz,kibz,Dtset%kptopt)
 wtk_ibz=Kibz_mesh%wt
 call cryst%free()
 call kmesh_free(Kibz_mesh)

 ABI_MALLOC(ibz2bz,(nkibz,6))
 call listkk(dksqmax,gmet,ibz2bz,dtset%kpt,dtset%kptgw,nkpt,dtset%nkptgw,&
&  dtset%nsym,sppoldbl,dtset%symafm,dtset%symrel,timrev,xmpi_comm_self)

 ltest=ALL(ibz2bz(:,2)==1)
 ABI_CHECK(ltest,'Not able to found irreducible points in the BZ set!')

 if (dksqmax>tol8) then
    write(msg,'(5a)')&
&     'Set of GW irreducible-zone kptgw in input file is inconsistent',ch10,&
&     'with full-zone set being used for wannier90 setup.',ch10,&
&     'Action: correct input data'
    MSG_ERROR(msg)
 end if

 ABI_MALLOC(npwarr_ibz,(nkibz))
 ABI_MALLOC(istwfk_ibz,(nkibz))
 ABI_MALLOC(nband_ibz,(nkibz*nsppol))

 do isppol=1,nsppol
   do ikibz=1,nkibz
     ikbz=ibz2bz(ikibz+(sppoldbl-1)*(isppol-1)*nkibz,1)
     npwarr_ibz(ikibz)=      npwarr(ikbz)
     istwfk_ibz(ikibz)=Dtset%istwfk(ikbz)
     nband_ibz(ikibz+(isppol-1)*nkibz)=Dtset%nband(ikbz+(isppol-1)*nkpt)
   end do
 end do

 bantot_ibz=SUM(nband_ibz)
 ABI_MALLOC(doccde_ibz,(bantot_ibz))
 ABI_MALLOC(eigen_ibz,(bantot_ibz))
 ABI_MALLOC(occfact_ibz,(bantot_ibz))
 doccde_ibz(:)=zero ; eigen_ibz(:)=zero ; occfact_ibz(:)=zero

 band_index=0
 do isppol=1,nsppol
   do ikibz=1,nkibz
     ikbz=ibz2bz(ikibz+(sppoldbl-1)*(isppol-1)*nkibz,1)
     nband_k=nband_ibz(ikibz+(isppol-1)*nkibz)
     ii=SUM(Dtset%nband(1:ikbz+(isppol-1)*nkpt))-nband_k
     eigen_ibz(band_index+1:band_index+nband_k)=eigen(ii+1:ii+nband_k)
     band_index=band_index+nband_k
   end do
 end do

 call ebands_init(bantot_ibz,QP_bst,Dtset%nelect,doccde_ibz,eigen_ibz,istwfk_ibz,kibz,nband_ibz,&
&  nkibz,npwarr_ibz,nsppol,Dtset%nspinor,Dtset%tphysel,Dtset%tsmear,Dtset%occopt,occfact_ibz,wtk_ibz,&
&  dtset%charge,dtset%kptopt,dtset%kptrlatt_orig,dtset%nshiftk_orig,dtset%shiftk_orig,&
&  dtset%kptrlatt,dtset%nshiftk,dtset%shiftk)

 ABI_FREE(kibz)
 ABI_FREE(wtk_ibz)
 ABI_FREE(ibz2bz)
 ABI_FREE(npwarr_ibz)
 ABI_FREE(istwfk_ibz)
 ABI_FREE(nband_ibz)
 ABI_FREE(doccde_ibz)
 ABI_FREE(eigen_ibz)
 ABI_FREE(occfact_ibz)

 ! === Read in quasiparticle information ===
 ! * Initialize QP amplitudes with KS, QP_bst% presently contains KS energies.
 ! * If file not found return, everything has been already initialized with KS values
 !   Here qp_rhor is not needed thus dimrho=0
 ABI_MALLOC(m_lda_to_qp,(mband,mband,dtset%nkptgw,nsppol))
 m_lda_to_qp=czero
 do iband=1,mband
   m_lda_to_qp(iband,iband,:,:)=cone
 end do

 ! Fake MPI_type for rdqps
 call initmpi_seq(MPI_enreg_seq)

 my_ngfft=Dtset%ngfft; if (Dtset%usepaw==1.and.ALL(Dtset%ngfftdg(1:3)/=0)) my_ngfft=Dtset%ngfftdg
 nfftot=PRODUCT(my_ngfft(1:3)); dimrho=0

 ! Change gw_fname to read a GW file instead of the QPS file.
 ! TODO not so sure that wannier90 can handle G0W0 eigenvalues that are not ordered, though!
 gw_fname = "g0w0"
 g0w0_exists = .FALSE.
 inquire(file=gw_fname,iostat=ios,exist=g0w0_exists)
 if (ios/=0) then
   MSG_ERROR('File g0w0 exists but iostat returns nonzero value!')
 end if

 if (.not.g0w0_exists) then ! read QPS file (default behavior).
   input = from_QPS_FILE
   ABI_DT_MALLOC(prev_Pawrhoij,(Cryst%natom*Dtset%usepaw))
   ABI_MALLOC(qp_rhor,(nfftot,nspden*dimrho))

   call rdqps(QP_bst,Dtfil%fnameabi_qps,Dtset%usepaw,Dtset%nspden,dimrho,nscf,&
    nfftot,my_ngfft,ucvol,Cryst,Pawtab,MPI_enreg_seq,nbsc,m_lda_to_qp,qp_rhor,prev_Pawrhoij)

   ABI_FREE(qp_rhor)
   ABI_DT_FREE(prev_Pawrhoij)

 else
   ! Read GW file (m_lda_to_qp has been already set to 1, no extrapolation is performed)
   MSG_WARNING(' READING GW CORRECTIONS FROM FILE g0w0 !')
   input = from_GW_FILE
   ABI_MALLOC(igwene,(QP_bst%mband,QP_bst%nkpt,QP_bst%nsppol))
   call rdgw(QP_bst,gw_fname,igwene,extrapolate=.FALSE.)
   ABI_FREE(igwene)
 end if

 ! === Begin big loop over full-zone k points and spin (not implemented) ===
 ! * Wannier90 treats only a single spin, changes in wannier90 are needed
 ABI_MALLOC(cg_k,(mpw,mband))
 ABI_MALLOC(cg_qpk,(mpw,mband))
 ABI_MALLOC(m_tmp,(mband,mband))

 band_index=0 ; icg=0 ; ikg=0
 do isppol=1,nsppol
   do ikpt=1,nkpt

    irzkpt =indkk(ikpt+(sppoldbl-1)*(isppol-1)*nkpt,1)
    itimrev=indkk(ikpt+(sppoldbl-1)*(isppol-1)*nkpt,6)
    npw_k=npwarr(ikpt)
    nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)

    if (nband_k/=mband) then
      write(msg,'(a,i0,7a)')&
       'Number of bands for k point ',ikpt,' is inconsistent with number',ch10,&
       'specified for wannier90 calculation',ch10,&
       'Action: correct input so all band numbers are equal for GW',ch10,&
       'and wannier90 datasets.'
      MSG_ERROR(msg)
    end if

    ! Load KS states for this kbz and spin
    do iband=1,nband_k
      icg_shift=npw_k*my_nspinor*(iband-1)+icg
      do ipw=1,npw_k
        cg_k(ipw,iband)=DCMPLX(cg(1,ipw+icg_shift),cg(2,ipw+icg_shift))
      end do
    end do

    ! If time reversal is used for relating ikpt to irzkpt, then multiply by
    ! the complex conjugate of the lda-to-qp transformation matrix
    if (itimrev==0) then
      m_tmp(:,:)=m_lda_to_qp(:,:,irzkpt,isppol)
    else if (itimrev==1) then
      m_tmp(:,:)=conjg(m_lda_to_qp(:,:,irzkpt,isppol))
    else
      write(msg,'(2(a,i0))')'Invalid indkk(ikpt,6) ',itimrev,'from routine listkk for k-point ',ikpt
      MSG_BUG(msg)
    end if

    call ZGEMM('N','N',npw_k,mband,mband,cone,cg_k,mpw,m_tmp,mband,czero,cg_qpk,mpw)

    ! === Orthonormality test ===
    ! * nband >= maxval(bndgw) for this to pass, but may be less than nband used in GW.
    ! * Unfortunately, does not test WFK and QPS consistency.
    !allocate(ortho(nband_k*(nband_k+1)/2))
    !ortho=czero; ijpack=0
    !do jband=1,nband_k
    !  jb_idx=band_index+jband
    !  if (dtset%usepaw==1) Cp2 => Cprj_BZ(:,jband:jband+(my_nspinor-1))
    !  do iband=1,jband
    !    ib_idx=band_index+iband
    !    ijpack=ijpack+1
    !    ortho(ijpack)=sum(conjg(cg_qpk(1:npw_k,iband))*cg_qpk(1:npw_k,jband))
    !    if (dtset%usepaw==1) then
    !      Cp1 => Cprj_BZ(:,iband:iband+(my_nspinor-1))
    !      paw_ovlp = paw_overlap(Cp2,Cp1,Cryst%typat,Pawtab)
    !      ortho(ijpack) = ortho(ijpack) + CMPLX(paw_ovlp(1),paw_ovlp(2))
    !    end if
    !    if (jband==iband) ortho(ijpack)=ortho(ijpack)-cone
    !  end do
    !end do
    !ortho_err=maxval(abs(ortho))

    !write(std_out,*)' drh - mlwfovlp_qp: ikpt,ortho_err',ikpt,ortho_err
    !if (ortho_err>tol6) then
    !  write(msg, '(3a,i4,a,i6,a,1p,e8.1,3a)' )&
    !&    '  orthonormality error for quasiparticle wave functions.',ch10,&
    !&    '  spin=',isppol,'  k point=',ikpt,'  ortho_err=',ortho_err,' >1E-6',ch10,&
    !&    '  Action: Be sure input nband>=maxval(bndgw)'
    !  MSG_ERROR(msg)
    !end if
    !deallocate(ortho)

    ! Replace lda wave functions and eigenvalues with quasiparticle ones.
    qpenek_is_ordered = isordered(nband_k,QP_bst%eig(:,irzkpt,isppol),">",TOL_SORT)

    if (input==from_QPS_FILE .and. .not.qpenek_is_ordered) then
      write(msg,'(3a)')&
      " QP energies read from QPS file are not ordered, likely nband_k>nbdgw. ",ch10,&
      " Change nband in the input file so that it equals the number of GW states calculated"
      MSG_WARNING(msg)
    end if

    if ( .TRUE. ) then
      do iband=1,nband_k
        icg_shift=npw_k*my_nspinor*(iband-1)+icg
        eigen(iband+band_index)=QP_bst%eig(iband,irzkpt,isppol)
        do ipw=1,npw_k
          cg(1,ipw+icg_shift)= real(cg_qpk(ipw,iband))
          cg(2,ipw+icg_shift)=aimag(cg_qpk(ipw,iband))
        end do
      end do
    else
      ! FIXME There's a problem in twannier90 since nband_k > nbdgw and therefore we also read KS states from the QPS file!
      ! Automatic test has to be changed!
      write(msg,'(2a,3f8.4,3a)')ch10,&
&       "QP energies at k-point ",QP_bst%kptns(:,irzkpt)," are not sorted in ascending numerical order!",ch10,&
&       "Performing reordering of energies and wavefunctions to be written on the final WKF file."
      MSG_ERROR(msg)
      !write(std_out,*)"eig",(QP_bst%eig(ii,irzkpt,isppol),ii=1,nband_k)
      ABI_MALLOC(sorted_qpene,(nband_k))
      ABI_MALLOC(iord,(nband_k))
      sorted_qpene = QP_bst%eig(1:nband_k,irzkpt,isppol)
      iord = (/(ii, ii=1,nband_k)/)

      call sort_dp(nband_k,sorted_qpene,iord,TOL_SORT)
      do ii=1,nband_k
        write(std_out,*)"%eig, sorted_qpene, iord",QP_bst%eig(ii,irzkpt,isppol)*Ha_eV,sorted_qpene(ii)*Ha_eV,iord(ii)
      end do

      do iband=1,nband_k
        ord_iband = iord(iband)
        icg_shift=npw_k*my_nspinor*(iband-1)+icg
        !eigen(iband+band_index)=QP_bst%eig(iband,irzkpt,isppol)
        eigen(iband+band_index)=QP_bst%eig(ord_iband,irzkpt,isppol)
        do ipw=1,npw_k
          !cg(1,ipw+icg_shift)= real(cg_qpk(ipw,iband))
          !cg(2,ipw+icg_shift)=aimag(cg_qpk(ipw,iband))
          cg(1,ipw+icg_shift)= real(cg_qpk(ipw,ord_iband))
          cg(2,ipw+icg_shift)=aimag(cg_qpk(ipw,ord_iband))
        end do
      end do
      ABI_FREE(sorted_qpene)
      ABI_FREE(iord)
    end if

    band_index=band_index+nband_k
    icg=icg+npw_k*my_nspinor*nband_k
    ikg=ikg+npw_k
   end do !ikpt
 end do !isppol

 ABI_FREE(cg_k)
 ABI_FREE(cg_qpk)
 ABI_FREE(m_tmp)

 ! === If PAW, update projections in BZ ===
 ! * Since I am lazy and here I do not care about memory, I just reconstruct m_lda_to_qp in the BZ.
 ! * update_cprj will take care of updating the PAW projections to get <p_lmn|QP_{nks]>
 !   This allows some CPU saving, no need to call ctocprj.
 ! FIXME this part should be tested, automatic test to be provided
 if (Dtset%usepaw==1) then
   ABI_MALLOC(dimlmn,(natom))
   call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,ntypat,Dtset%typat,pawtab,'R')

   nkbz=nkpt
   ABI_MALLOC(m_lda_to_qp_BZ,(mband,mband,nkbz,nsppol))
   do isppol=1,nsppol
     do ikbz=1,nkbz
       ikibz  =indkk(ikibz+(sppoldbl-1)*(isppol-1)*nkbz,1)
       itimrev=indkk(ikibz+(sppoldbl-1)*(isppol-1)*nkbz,6)
       select case (itimrev)
       case (0)
         m_lda_to_qp_BZ(:,:,ikbz,isppol)=m_lda_to_qp(:,:,ikibz,isppol)
       case (1)
         m_lda_to_qp_BZ(:,:,ikbz,isppol)=CONJG(m_lda_to_qp(:,:,ikibz,isppol))
       case default
         write(msg,'(a,i3)')"Wrong itimrev= ",itimrev
         MSG_BUG(msg)
       end select
     end do
   end do

   call update_cprj(natom,nkbz,mband,nsppol,my_nspinor,m_lda_to_qp_BZ,dimlmn,Cprj_BZ)
   ABI_FREE(dimlmn)
   ABI_FREE(m_lda_to_qp_BZ)
 end if !PAW

 write(msg,'(6a)')ch10,&
  ' mlwfovlp_qp: Input KS wavefuctions have been converted',ch10,&
  '  to GW quasiparticle wavefunctions for maximally localized wannier',ch10,&
  '  function construction by wannier90.'
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 ABI_FREE(m_lda_to_qp)
 call ebands_free(QP_bst)
 call destroy_mpi_enreg(MPI_enreg_seq)

 DBG_EXIT("COLL")

end subroutine mlwfovlp_qp
!!***

!!****f* ABINIT/update_cprj
!! NAME
!! update_cprj
!!
!! FUNCTION
!!  Update the matrix elements of the PAW projectors in case of self-consistent GW.
!!
!! INPUTS
!!  dimlmn(natom)=number of (l,m,n) components for each atom (only for PAW)
!!  nkibz=number of k-points
!!  nsppol=number of spin
!!  nbnds=number of bands in the present GW calculation
!!  m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)= expansion of the QP amplitudes in terms of KS wavefunctions
!!  natom=number of atomd in unit cell
!!
!! OUTPUT
!!  Cprj_ibz(natom,nspinor*nkibz*nbnds*nsppol) <type(pawcprj_type)>=projected wave functions
!!   <Proj_i|Cnk> with all NL projectors. On exit, it contains the projections onto the
!!   QP amplitudes.
!!
!! TODO
!! To be moved to cprj_utils, although here we use complex variables.
!!
!! PARENTS
!!      mlwfovlp_qp
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free
!!
!! SOURCE
!!

subroutine update_cprj(natom,nkibz,nbnds,nsppol,nspinor,m_lda_to_qp,dimlmn,Cprj_ibz)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nbnds,nkibz,nsppol,nspinor
!arrays
 integer,intent(in) :: dimlmn(natom)
 complex(dpc),intent(in) :: m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)
 type(pawcprj_type),intent(inout) :: Cprj_ibz(natom,nspinor*nbnds*nkibz*nsppol)

!Local variables-------------------------------
!scalars
 integer :: iat,ib,ik,is,shift,indx_kibz,ilmn,nlmn,ispinor,ibsp,spad,ibdx
!arrays
 real(dp),allocatable :: re_p(:),im_p(:),vect(:,:),umat(:,:,:)
 type(pawcprj_type),allocatable :: Cprj_ks(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 ABI_DATATYPE_ALLOCATE(Cprj_ks,(natom,nspinor*nbnds))
 call pawcprj_alloc(Cprj_ks,0,dimlmn)

 ABI_ALLOCATE(re_p,(nbnds))
 ABI_ALLOCATE(im_p,(nbnds))
 ABI_ALLOCATE(vect,(2,nbnds))
 ABI_ALLOCATE(umat,(2,nbnds,nbnds))
 !
 ! $ \Psi^{QP}_{r,b} = \sum_n \Psi^{KS}_{r,n} M_{n,b} $
 !
 ! therefore the updated PAW projections are given by:
 !
 ! $ \<\tprj_j|\Psi^{QP}_a\> = sum_b M_{b,a} <\tprj_j|\Psi^{KS}_b\> $.
 !
 do is=1,nsppol
   do ik=1,nkibz

    shift=nspinor*nbnds*nkibz*(is-1)
    indx_kibz=nspinor*nbnds*(ik-1)+shift
    ibsp=0
    do ib=1,nbnds
      do ispinor=1,nspinor
        ibsp=ibsp+1
        do iat=1,natom
          Cprj_ks(iat,ibsp)%cp(:,:)=Cprj_ibz(iat,indx_kibz+ibsp)%cp(:,:)
        end do
      end do
    end do

    umat(1,:,:)=TRANSPOSE( REAL (m_lda_to_qp(:,:,ik,is)) )
    umat(2,:,:)=TRANSPOSE( AIMAG(m_lda_to_qp(:,:,ik,is)) )

    do iat=1,natom
      nlmn=dimlmn(iat)
      do ilmn=1,nlmn

        do ispinor=1,nspinor
           ! * Retrieve projections for this spinor component, at fixed atom and ilmn.
           spad=(ispinor-1)
           ibdx=0
           do ib=1,nbnds*nspinor,nspinor
            ibdx=ibdx+1
            vect(1,ibdx)=Cprj_ks(iat,ib+spad)%cp(1,ilmn)
            vect(2,ibdx)=Cprj_ks(iat,ib+spad)%cp(2,ilmn)
           end do

           re_p(:)= &
&            MATMUL(umat(1,:,:),vect(1,:)) &
&           -MATMUL(umat(2,:,:),vect(2,:))

           im_p(:)= &
&            MATMUL(umat(1,:,:),vect(2,:)) &
&           +MATMUL(umat(2,:,:),vect(1,:))

           ! === Save values ===
           ibdx=0
           do ib=1,nbnds*nspinor,nspinor
            ibdx=ibdx+1
            Cprj_ibz(iat,indx_kibz+spad+ib)%cp(1,ilmn)=re_p(ibdx)
            Cprj_ibz(iat,indx_kibz+spad+ib)%cp(2,ilmn)=im_p(ibdx)
           end do
        end do !ispinor

      end do !ilmn
    end do !iat

   end do !ik
 end do !is

 ABI_DEALLOCATE(re_p)
 ABI_DEALLOCATE(im_p)
 ABI_DEALLOCATE(vect)
 ABI_DEALLOCATE(umat)

 call pawcprj_free(Cprj_ks)
 ABI_DATATYPE_DEALLOCATE(Cprj_ks)

 DBG_EXIT("COLL")

end subroutine update_cprj
!!***

end module m_mlwfovlp_qp
!!***
