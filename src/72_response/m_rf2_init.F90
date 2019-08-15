!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_rf2_init
!! NAME
!!  m_rf2_init
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2015-2019 ABINIT group (LB,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! TODO
!!  Can be merged with m_rf2
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rf2_init

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_wfk
 use m_hamiltonian
 use m_cgtools
 use m_rf2
 use m_abicore
 use m_dtset
 use m_dtfil

 use m_time   , only : timab
 use m_pawcprj, only : pawcprj_type,pawcprj_alloc,pawcprj_copy,pawcprj_get,pawcprj_free,pawcprj_output
 use m_cgprj,   only : getcprj

 implicit none

 private
!!***

 public :: rf2_init
!!***

contains
!!***

!!****f* ABINIT/rf2_init
!!
!! NAME
!! rf2_init
!!
!! FUNCTION
!! Compute terms needed for the 2nd order Sternheimer equation.
!! All terms are stored in a rf2_t object.
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*nsppol)=planewave coefficients of wavefunctions at k
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!              projected with non-local projectors: cprj=<p_i|Cnk>
!!  rf2 : the object we want to initialize (see m_rf2.F90 for more information)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eig0_k(mband*nsppol)=GS eigenvalues at k (hartree)
!!  eig1_k(2*mband*mband*nsppol)=2nd-order eigenvalues at k,q (hartree)
!!  ffnl1=nonlocal form factors
!!  ffnl1_test=nonlocal form factors used for tests (i.e. when dtset%nonlinear_info>2)
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  ibg=shift to be applied on the location of data in the array cprj
!!  icg=shift to be applied on the location of data in the array cg
!!  idir=direction of the perturbation
!!  ikpt=number of the k-point
!!  ipert=type of the perturbation
!!  isppol=index of current spin component
!!  mkmem =number of k points trated by this node (GS data).
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  nband_k=number of bands at this k point for that spin polarization
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  rf_hamkq <type(rf_hamiltonian_type)>=all data for the 1st-order Hamiltonian at k,q
!!  rf_hamk_dir2 <type(rf_hamiltonian_type)>= (used only when ipert=natom+11, so q=0)
!!    same as rf_hamkq, but the direction of the perturbation is different
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  rocceig(nband_k,nband_k)= (occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n))
!!  ddk_f<wfk_t>=struct info for DDK file.
!!
!! OUTPUT
!!  rf2%RHS_Stern
!!
!! NOTES
!!
!! PARENTS
!!      dfpt_vtowfk
!!
!! CHILDREN
!!      cg_zaxpy,dotprod_g,getcprj,pawcprj_alloc,pawcprj_free,pawcprj_get
!!      rf2_accumulate_bands,rf2_apply_hamiltonian,rf2_getidirs,sqnorm_g
!!      wfk_read_bks,wrtout,xmpi_allgather,xmpi_barrier
!!
!! SOURCE

subroutine rf2_init(cg,cprj,rf2,dtset,dtfil,eig0_k,eig1_k,ffnl1,ffnl1_test,gs_hamkq,ibg,icg,idir,ikpt,ipert,isppol,mkmem,&
                     mpi_enreg,mpw,nband_k,nsppol,rf_hamkq,rf_hamk_dir2,occ_k,rocceig,ddk_f)

! *************************************************************************
!Arguments -------------------------------
!scalars
 integer,intent(in) :: ibg,icg,idir,ipert,isppol,ikpt
 integer,intent(in) :: mkmem,mpw,nband_k,nsppol
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout),target :: rf_hamkq,rf_hamk_dir2
 type(MPI_type),intent(in) :: mpi_enreg

!arrays
 real(dp),intent(in),target :: cg(2,mpw*gs_hamkq%nspinor*dtset%mband*mkmem*nsppol)
 real(dp),intent(in) :: eig0_k(dtset%mband)
 real(dp),intent(inout) :: eig1_k(2*dtset%mband**2) ! Here eig1_k contains 2nd order eigenvalues...
 real(dp),intent(in) :: ffnl1(:,:,:,:),ffnl1_test(:,:,:,:)
 real(dp),intent(in) :: occ_k(nband_k),rocceig(nband_k,nband_k)
 type(pawcprj_type),intent(in) :: cprj(gs_hamkq%natom,gs_hamkq%nspinor*dtset%mband*mkmem*nsppol*gs_hamkq%usecprj)
 type(rf2_t),intent(inout) :: rf2
 type(wfk_t),intent(inout) :: ddk_f(4)
!
!Local variables-------------------------------
!scalars
 integer,parameter :: iorder_cprj=0
 integer :: choice_cprj,cpopt_cprj,iband,icpgr_loc,idir1,idir2,idir_cprj,ierr
 integer :: indb,ipert1,ipert2,iproc,jband,kdir1
 integer :: me,my_nband,natom,ncpgr_loc,nproc_band,debug_mode
 integer :: size_cprj,size_wf,shift_band1,shift_band2,shift_cprj_band1,shift_cprj_dir1,shift_proc
 integer :: shift_dir1_lambda,shift_dir2_lambda,shift_dir1,shift_dir1_loc,shift_dir2,shift_jband_lambda
 logical :: has_cprj_jband,has_dudkprj
 real(dp) :: doti,dotr,dot2i,dot2r,invocc,tol_final,factor
 character(len=500) :: msg
!arrays
 integer :: file_index(2)
 real(dp) :: lambda_ij(2),tsec(2)
 real(dp),allocatable :: cg_jband(:,:,:),ddk_read(:,:),dudkdk(:,:),dudk_dir2(:,:)
 real(dp),allocatable :: eig1_read(:),gvnlx1(:,:),h_cwave(:,:),s_cwave(:,:),dsusdu_loc(:,:),dsusdu_gather(:,:)
 real(dp),allocatable,target :: dsusdu(:,:),dudk(:,:),eig1_k_stored(:)
 real(dp), ABI_CONTIGUOUS pointer :: cwave_dudk(:,:),cwave_i(:,:),cwave_j(:,:),eig1_k_jband(:)
 real(dp),pointer :: rhs_j(:,:)
 type(pawcprj_type),target :: cprj_empty(0,0)
 type(pawcprj_type),allocatable,target :: cprj_jband(:,:),dudkprj(:,:)
 type(pawcprj_type),pointer :: cprj_dudk(:,:),cprj_j(:,:)
 type(rf_hamiltonian_type),pointer :: rf_hamk_idir

! *********************************************************************

 DBG_ENTER("COLL")

 call timab(514,1,tsec)

!my mpi rank :
 me=mpi_enreg%me_kpt

 size_wf=gs_hamkq%npw_k*gs_hamkq%nspinor
 size_cprj=gs_hamkq%nspinor
 natom = gs_hamkq%natom
 debug_mode = 0
 if (dtset%nonlinear_info>2) debug_mode = 1 ! also active a lot of tests

!Define some attributes of the rf2 object
 rf2%nband_k = nband_k
 rf2%size_wf = size_wf
 rf2%size_cprj = size_cprj

 if(ipert<natom+10.or.ipert>natom+11) then
   write(msg,'(a)') 'ipert must be equal to natom+10 or natom+11 for rf2 calculations.'
   MSG_BUG(msg)
 end if

!Define perturbations and idirs
 rf2%iperts(1) = natom+1
 rf2%iperts(2) = natom+1
 if (ipert==natom+11)  rf2%iperts(2) = natom+2

 if (ipert==natom+10.and.idir<=3) then ! One perturbation, one direction
   rf2%ndir=1
   rf2%idirs(1)=idir ; rf2%idirs(2)=idir
 else ! Two perturbations or/and two directions
   rf2%ndir=2
   call rf2_getidirs(idir,idir1,idir2)
   rf2%idirs(1)=idir1
   rf2%idirs(2)=idir2
 end if

! **************************************************************************************************
! Get info from ddk files
! **************************************************************************************************

!Allocate work spaces
 ABI_ALLOCATE(eig1_read,(2*nband_k))
 ABI_ALLOCATE(ddk_read,(2,size_wf))
 eig1_read(:)=zero
 ddk_read(:,:)=zero

! "eig1_k_stored" contains dLambda_{nm}/dpert every bands n and m and ndir (=1 or 2) directions
! pert = k_dir (wavevector) or E_dir (electric field)
 ABI_ALLOCATE(eig1_k_stored,(2*rf2%ndir*nband_k**2))
 eig1_k_stored=zero

! "dudk" contains du/dpert1 for every bands and ndir (=1 or 2) directions
 ABI_MALLOC_OR_DIE(dudk,(2,rf2%ndir*nband_k*size_wf), ierr)
 dudk=zero
 has_dudkprj=.false.
 if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
   ABI_DATATYPE_ALLOCATE(dudkprj,(natom,rf2%ndir*nband_k*size_cprj))
   ncpgr_loc=1;if(ipert==natom+10.or.ipert==natom+11) ncpgr_loc=3
   call pawcprj_alloc(dudkprj,ncpgr_loc,gs_hamkq%dimcprj)
   choice_cprj=5 ; cpopt_cprj=0
   has_dudkprj=.true.
 else
   ABI_DATATYPE_ALLOCATE(dudkprj,(natom,0))
 end if

 if (debug_mode/=0) then
   write(msg,'(4(a,i2))') 'RF2_INIT : ipert-natom = ',ipert-natom,' , idir = ',idir,&
   ' , ikpt = ',ikpt,' , isppol = ',isppol
   call wrtout(std_out,msg,'COLL')
 end if

 file_index(1)=1 ! dir1
 file_index(2)=2 ! dir2
 if (ipert==natom+11) then ! see dfpt_looppert.F90
   file_index(1)=3 ! dir1
   file_index(2)=2 ! dir2
 end if

 do kdir1=1,rf2%ndir
   idir1=rf2%idirs(kdir1)
   ipert1=rf2%iperts(kdir1)
   do iband=1,nband_k
     call ddk_f(file_index(kdir1))%read_bks(iband,ikpt,isppol,xmpio_single,cg_bks=ddk_read,eig1_bks=eig1_read)
!    Copy ddk_read in "dudk"
     shift_band1=(iband-1)*size_wf
     shift_dir1=(kdir1-1)*nband_k*size_wf
     dudk(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)=ddk_read(:,:)
!    Copy eig1_read in "eig1_k_stored"
     shift_band1=(iband-1)*2*nband_k
     shift_dir1_lambda=2*(kdir1-1)*nband_k**2
     eig1_k_stored(1+shift_band1+shift_dir1_lambda:2*nband_k+shift_band1+shift_dir1_lambda)=eig1_read(:)
!    Get this dudk projected on NL projectors
     if (has_dudkprj.and.mpi_enreg%proc_distrb(ikpt,iband,isppol)==me) then
       shift_cprj_band1=(iband-1)*size_cprj
       shift_cprj_dir1=(kdir1-1)*nband_k*size_cprj
       cprj_dudk => dudkprj(:,1+shift_cprj_band1+shift_cprj_dir1: &
&       size_cprj+shift_cprj_band1+shift_cprj_dir1)
       idir_cprj=0;if (dudkprj(1,1)%ncpgr/=3) idir_cprj=idir1
       call getcprj(choice_cprj,cpopt_cprj,ddk_read,cprj_dudk,gs_hamkq%ffnl_k,idir_cprj,&
&       gs_hamkq%indlmn,gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kpg_k,gs_hamkq%kpt_k,&
&       gs_hamkq%lmnmax,gs_hamkq%mgfft,mpi_enreg,gs_hamkq%natom,gs_hamkq%nattyp,gs_hamkq%ngfft,&
&       gs_hamkq%nloalg,gs_hamkq%npw_k,gs_hamkq%nspinor,gs_hamkq%ntypat,gs_hamkq%phkxred,&
&       gs_hamkq%ph1d,gs_hamkq%ph3d_k,gs_hamkq%ucvol,gs_hamkq%useylm)
     end if
   end do
 end do

 ABI_ALLOCATE(dudkdk,(2,0))
 ABI_ALLOCATE(dudk_dir2,(2,0))

!Get dudkdk for ipert==natom+11
 if(ipert==natom+11) then
   ABI_DEALLOCATE(dudkdk)
   ABI_ALLOCATE(dudkdk,(2,nband_k*size_wf))
   if (idir>3) then
     ABI_DEALLOCATE(dudk_dir2)
     ABI_ALLOCATE(dudk_dir2,(2,nband_k*size_wf))
   end if
   do iband=1,nband_k
     call ddk_f(1)%read_bks(iband,ikpt,isppol,xmpio_single,cg_bks=ddk_read,eig1_bks=eig1_read)
     shift_band1=(iband-1)*size_wf
     dudkdk(:,1+shift_band1:size_wf+shift_band1)=ddk_read(:,:)
!    Check that < u^(0) | u^(2) > = - Re[< u^(1) | u^(1) >]
     if (debug_mode/=0 .and. gs_hamkq%usepaw==0) then
!      Compute < u^(0) | u^(2) >
       do jband=1,nband_k
         cwave_j => cg(:,1+shift_band1+icg:size_wf+shift_band1+icg)
         call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cwave_j,ddk_read,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
         if (idir<=3 .and. iband==jband .and. abs(occ_k(iband))>tol8) then
!          Compute < u^(1) | u^(1) > = Re[< u^(1) | u^(1) >]
           cwave_dudk => dudk(:,1+shift_band1:size_wf+shift_band1)
           call sqnorm_g(dot2r,gs_hamkq%istwf_k,size_wf,cwave_dudk,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
           dotr = dotr + dot2r
           dotr = sqrt(dotr**2+doti**2)
           if (dotr > tol7) then
             write(msg,'(a,i2,a,es22.13E3)') 'RF2 TEST dudkdk iband = ',iband,&
             ' : NOT PASSED. | < u^(0) | u^(2) > + Re[< u^(1) | u^(1) >] | = ',dotr
             call wrtout(std_out,msg)
             call wrtout(ab_out,msg)
!           else
!             write(msg,'(a,i2,a)') 'RF2 TEST dudkdk iband = ',iband,' : OK.'
!             call wrtout(std_out,msg)
           end if
         end if ! idir<=3
       end do ! jband
     end if ! debug_mode
!    Read ddk for idir2
     if (idir>3) then
       call ddk_f(4)%read_bks(iband,ikpt,isppol,xmpio_single,cg_bks=ddk_read,eig1_bks=eig1_read)
       dudk_dir2(:,1+shift_band1:size_wf+shift_band1)=ddk_read(:,:)
     end if
   end do !iband
 end if ! ipert=natom+11

 ABI_DEALLOCATE(ddk_read)
 ABI_DEALLOCATE(eig1_read)

! **************************************************************************************************
! COMPUTATION OF "dsusdu", A PART OF "A_mn" AND A PART OF "Lambda_mn" (see defs in m_rf2)
! **************************************************************************************************

!Allocate work spaces for one band
 ABI_ALLOCATE(h_cwave,(2,size_wf))
 ABI_ALLOCATE(s_cwave,(2,size_wf))
 ABI_ALLOCATE(gvnlx1,(2,size_wf))
 h_cwave(:,:) = zero
 s_cwave(:,:) = zero
 gvnlx1(:,:) = zero

! "dsusdu" contains dS/dpert_dir |u_band> + S|du_band/dpert1> for every bands and ndir (=1 or 2) directions
 ABI_MALLOC_OR_DIE(dsusdu,(2,rf2%ndir*nband_k*size_wf), ierr)
 dsusdu=zero

 ABI_ALLOCATE(rf2%amn,(2,nband_k**2))
 rf2%amn=zero

 ABI_ALLOCATE(rf2%lambda_mn,(2,nband_k**2))
 rf2%lambda_mn(:,:)=zero

!Allocate work spaces when debug_mode is activated
 has_cprj_jband=.false.
 if (debug_mode/=0) then ! Only for test purposes
   ABI_ALLOCATE(cg_jband,(2,size_wf*nband_k,2))
   cg_jband(:,:,1) = cg(:,1+icg:size_wf*nband_k+icg)
   if (ipert==natom+11) then ! Note the multiplication by "i"
     if (idir<=3) then
       cg_jband(1,:,2) = -dudk(2,1:size_wf*nband_k) ! for dir1
       cg_jband(2,:,2) =  dudk(1,1:size_wf*nband_k) ! for dir1
     else
       cg_jband(1,:,2) = -dudk_dir2(2,1:size_wf*nband_k) ! for dir2
       cg_jband(2,:,2) =  dudk_dir2(1,1:size_wf*nband_k) ! for dir2
     end if
   end if
   if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
     ABI_DATATYPE_ALLOCATE(cprj_jband,(natom,size_cprj*nband_k))
     has_cprj_jband=.true.
   else
     ABI_DATATYPE_ALLOCATE(cprj_jband,(natom,0))
   end if
 else
   ABI_ALLOCATE(cg_jband,(2,0,2))
   ABI_DATATYPE_ALLOCATE(cprj_jband,(natom,0))
 end if

 factor=one
 if(ipert==natom+10 .and. idir<=3) factor=two ! in order to not compute same terms twice

 do kdir1=1,rf2%ndir
!  First iteration (kdir1=1) :
!  pert1 = rf2%iperts(1) along rf2%idirs(1)
!  pert2 = rf2%iperts(2) along rf2%idirs(2)
!  Second iteration (kdir1=2) :
!  pert1 = rf2%iperts(2) along rf2%idirs(2)
!  pert2 = rf2%iperts(1) along rf2%idirs(1)
   idir1=rf2%idirs(kdir1)
   ipert1=rf2%iperts(kdir1)
   shift_dir1=(kdir1-1)*nband_k*size_wf
   shift_cprj_dir1=(kdir1-1)*nband_k*size_cprj
   shift_dir1_lambda=(kdir1-1)*2*nband_k**2
   if(ipert==natom+10 .and. idir<=3) then
     shift_dir2=0
     idir2 = idir1
     ipert2 = ipert1
     rf_hamk_idir => rf_hamkq
   else
     shift_dir2=(2-kdir1)*nband_k*size_wf
     idir2 = rf2%idirs(3-kdir1)
     ipert2 = rf2%iperts(3-kdir1)
     if (kdir1==1) rf_hamk_idir => rf_hamkq
     if (kdir1==2) rf_hamk_idir => rf_hamk_dir2
   end if

!  Load projected WF according to ipert1 and idir1
   cprj_j => cprj_empty ; cprj_dudk => cprj_empty
   if (has_cprj_jband) then
     call pawcprj_free(cprj_jband)
     ncpgr_loc= 3;if(ipert1==natom+1.or.ipert1==natom+2) ncpgr_loc=1
     icpgr_loc=-1;if(ipert1==natom+1.or.ipert1==natom+2) icpgr_loc=idir1
     call pawcprj_alloc(cprj_jband,ncpgr_loc,gs_hamkq%dimcprj)
     call pawcprj_get(gs_hamkq%atindx1,cprj_jband,cprj,natom,1,ibg,ikpt,iorder_cprj,&
&     isppol,dtset%mband,mkmem,natom,nband_k,nband_k,gs_hamkq%nspinor,nsppol,dtfil%unpaw,&
&     mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb,ncpgr=3,icpgr=icpgr_loc)
   end if

!  LOOP OVER BANDS
   do jband=1,nband_k ! = band n

!    Skip bands not treated by current proc
     if(mpi_enreg%proc_distrb(ikpt,jband,isppol)/=me) cycle

     shift_band1=(jband-1)*size_wf
     shift_cprj_band1=(jband-1)*size_cprj
     shift_jband_lambda=(jband-1)*2*nband_k

     if (abs(occ_k(jband))>tol8) then

!      Extract first order wavefunction and eigenvalues for jband
       eig1_k_jband => eig1_k_stored(1+shift_jband_lambda+shift_dir1_lambda:2*nband_k+shift_jband_lambda+shift_dir1_lambda)
       cwave_dudk => dudk(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)
       if (has_dudkprj) cprj_dudk => dudkprj(:,1+shift_cprj_band1+shift_cprj_dir1:size_cprj+shift_cprj_band1+shift_cprj_dir1)

!      Compute H^(0) | du/dpert1 > (in h_cwave) and S^(0) | du/dpert1 > (in s_cwave)
       call rf2_apply_hamiltonian(cg_jband,cprj_jband,cwave_dudk,cprj_dudk,h_cwave,s_cwave,&
&       eig0_k,eig1_k_jband,jband,gs_hamkq,gvnlx1,0,0,ikpt,isppol,mkmem,&
&       mpi_enreg,nband_k,nsppol,debug_mode,dtset%prtvol,rf_hamk_idir,size_cprj,size_wf)

       if (gs_hamkq%usepaw==0) s_cwave(:,:)=cwave_dudk(:,:) ! Store | du/dpert1 > in s_cwave

!      Copy infos in dsusdu
       dsusdu(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)=s_cwave(:,:)&
       +dsusdu(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)

       if (debug_mode/=0) then
         write(msg,'(2(a,i2))') 'RF2 TEST before accumulate_bands choice = 1 kdir1 = ',kdir1,' jband = ',jband
         call wrtout(std_out,msg)
       end if

!      For every occupied iband, we compute :
!      < du/dpert2(iband) | H^(0) | du/dpert1(jband) > and add it to lambda_mn
!      < du/dpert2(iband) | S^(0) | du/dpert1(jband) > and add it to amn
       do iband=1,rf2%nband_k  ! = band m
         if (abs(occ_k(iband))>tol8) then
           shift_band2=(iband-1)*size_wf
           cwave_dudk => dudk(:,1+shift_band2+shift_dir2:size_wf+shift_band2+shift_dir2)
           call rf2_accumulate_bands(rf2,1,gs_hamkq,mpi_enreg,iband,idir1,idir2,ipert1,ipert2,&
           jband,debug_mode,cwave_dudk,h_cwave,s_cwave)
         end if
       end do

!      Extract GS wavefunction for jband
       cwave_j => cg(:,1+shift_band1+icg:size_wf+shift_band1+icg)
       if(has_cprj_jband) cprj_j => cprj_jband(:,1+shift_cprj_band1:size_cprj+shift_cprj_band1)

       if (ipert1==natom+2) then
!        Extract ddk and multiply by i :
         if(idir<=3) then ! in this case : idir1=idir2
           gvnlx1(1,:) = -dudk(2,1+shift_band1:size_wf+shift_band1)
           gvnlx1(2,:) =  dudk(1,1+shift_band1:size_wf+shift_band1)
         else
           gvnlx1(1,:) = -dudk_dir2(2,1+shift_band1:size_wf+shift_band1)
           gvnlx1(2,:) =  dudk_dir2(1,1+shift_band1:size_wf+shift_band1)
         end if
       end if

!      Compute dH/dpert1 | u^(0) > (in h_cwave) and dS/dpert1 | u^(0) > (in s_cwave)
       call rf2_apply_hamiltonian(cg_jband,cprj_jband,cwave_j,cprj_j,h_cwave,s_cwave,&
&       eig0_k,eig1_k_jband,jband,gs_hamkq,gvnlx1,idir1,ipert1,ikpt,isppol,&
&       mkmem,mpi_enreg,nband_k,nsppol,debug_mode,dtset%prtvol,rf_hamk_idir,size_cprj,size_wf)

!      Copy infos in dsusdu
       if (gs_hamkq%usepaw==1) then
         dsusdu(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)=s_cwave(:,:)&
         +dsusdu(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)
       end if

       if (debug_mode/=0) then
         write(msg,'(2(a,i2))') 'RF2 TEST before accumulate_bands choice = 2 kdir1 = ',kdir1,' jband = ',jband
         call wrtout(std_out,msg)
       end if

!      For every occupied iband, we compute :
!      < du/dpert2(iband) | dH/dpert1 | u^(0)(jband) > and add it to lambda_mn
!      < du/dpert2(iband) | dS/dpert1 | u^(0)(jband) > and add it to amn
       do iband=1,rf2%nband_k  ! = band m
         if (abs(occ_k(iband))>tol8) then
           shift_band2=(iband-1)*size_wf
           cwave_dudk => dudk(:,1+shift_band2+shift_dir2:size_wf+shift_band2+shift_dir2)
           call rf2_accumulate_bands(rf2,2,gs_hamkq,mpi_enreg,iband,idir1,idir2,ipert1,ipert2,&
           jband,debug_mode,cwave_dudk,h_cwave,s_cwave)
         end if
       end do

     end if ! empty band test
   end do ! jband
 end do ! idir1

! Allgather dsusdu
 nproc_band = xmpi_comm_size(mpi_enreg%comm_band)
 if (nproc_band>1) then

   my_nband = nband_k/nproc_band;if (mod(nband_k,nproc_band)/=0) my_nband=my_nband+1
   ABI_ALLOCATE(dsusdu_loc,(2,size_wf*my_nband*rf2%ndir))
   ABI_ALLOCATE(dsusdu_gather,(2,size_wf*my_nband*rf2%ndir*nproc_band))
   dsusdu_loc(:,:) = zero
   dsusdu_gather(:,:) = zero

   do kdir1=1,rf2%ndir
     indb = 1
     shift_dir1=(kdir1-1)*size_wf*nband_k
     shift_dir1_loc=(kdir1-1)*size_wf*my_nband
     do jband=1,nband_k
!      Skip bands not treated by current proc
       if(mpi_enreg%proc_distrb(ikpt,jband,isppol)/=me) cycle

       shift_band1=(jband-1)*size_wf
       dsusdu_loc(:,indb+shift_dir1_loc:indb-1+size_wf+shift_dir1_loc) = &
       dsusdu(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)
       indb = indb + size_wf
     end do
   end do

   call xmpi_allgather(dsusdu_loc,2*size_wf*my_nband*rf2%ndir,dsusdu_gather,mpi_enreg%comm_band,ierr)

   do kdir1=1,rf2%ndir
     shift_dir1=(kdir1-1)*size_wf*nband_k
     shift_dir1_loc=(kdir1-1)*size_wf*my_nband
     do iproc=1,nproc_band
       shift_proc = (iproc-1)*size_wf*my_nband*rf2%ndir
       indb = 1
       do jband=1,my_nband
         iband = jband+(iproc-1)*my_nband
         if(iband<=nband_k) then
           shift_band1=(iband-1)*size_wf
           dsusdu(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1) = &
           dsusdu_gather(:,indb+shift_dir1_loc+shift_proc:indb-1+size_wf+shift_dir1_loc+shift_proc)
         end if
         indb = indb + size_wf
       end do
     end do
   end do
   ABI_DEALLOCATE(dsusdu_loc)
   ABI_DEALLOCATE(dsusdu_gather)

 end if

! **************************************************************************************************
! COMPUTATION OF "RHS_Stern", THE LAST PART OF "A_mn" AND A PART OF "Lambda_mn"
! **************************************************************************************************

 ABI_MALLOC_OR_DIE(rf2%RHS_Stern,(2,nband_k*size_wf), ierr)
 rf2%RHS_Stern(:,:)=zero

!Computation of terms containing H^(2)
 if (ipert/=natom+11 .or. gs_hamkq%usepaw==1) then ! Otherwise H^(2) = 0

!Load projected WF according to ipert and idir
   cprj_j => cprj_empty
   if (has_cprj_jband) then
     call pawcprj_free(cprj_jband)
     ncpgr_loc= 3;if(ipert==natom+1.or.ipert==natom+2) ncpgr_loc=1
     icpgr_loc=-1;if(ipert==natom+1.or.ipert==natom+2) icpgr_loc=idir
     call pawcprj_alloc(cprj_jband,ncpgr_loc,gs_hamkq%dimcprj)
     call pawcprj_get(gs_hamkq%atindx1,cprj_jband,cprj,natom,1,ibg,ikpt,iorder_cprj,&
&     isppol,dtset%mband,mkmem,natom,nband_k,nband_k,gs_hamkq%nspinor,nsppol,dtfil%unpaw,&
&     mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb,ncpgr=3,icpgr=icpgr_loc)
   end if

   if (ipert==natom+10) then
     rf_hamk_idir => rf_hamkq !     all info are in rf_hamkq
   else if (ipert==natom+11) then
     rf_hamk_idir => rf_hamk_dir2 ! all info are in rf_hamk_dir2
   end if

   do jband=1,nband_k

!    Skip bands not treated by current proc
     if(mpi_enreg%proc_distrb(ikpt,jband,isppol)/=me) cycle

     if (abs(occ_k(jband))>tol8) then
       shift_band1=(jband-1)*size_wf
       shift_cprj_band1=(jband-1)*size_cprj
       shift_jband_lambda=(jband-1)*2*nband_k

!      Extract GS wavefunction
       cwave_j => cg(:,1+shift_band1+icg:size_wf+shift_band1+icg)
       if(has_cprj_jband) cprj_j => cprj_jband(:,1+shift_cprj_band1:size_cprj+shift_cprj_band1)

!      Not used here, but a null pointer is not allowed... (called by rf2_apply_hamiltonian but not used)
       eig1_k_jband => eig1_k_stored(1+shift_jband_lambda:2*nband_k+shift_jband_lambda)

       if (ipert==natom+11) then
!        Extract ddk and multiply by i :
         if(idir<=3) then ! in this case : idir1=idir2
           gvnlx1(1,:) = -dudk(2,1+shift_band1:size_wf+shift_band1)
           gvnlx1(2,:) =  dudk(1,1+shift_band1:size_wf+shift_band1)
         else
           gvnlx1(1,:) = -dudk_dir2(2,1+shift_band1:size_wf+shift_band1)
           gvnlx1(2,:) =  dudk_dir2(1,1+shift_band1:size_wf+shift_band1)
         end if
       end if

!      Compute  : d^2H/(dpert1 dpert2)|u^(0)>  (in h_cwave)
!      and      : d^2S/(dpert1 dpert2)|u^(0)>  (in s_cwave)
       call rf2_apply_hamiltonian(cg_jband,cprj_jband,cwave_j,cprj_j,h_cwave,s_cwave,&
&       eig0_k,eig1_k_jband,jband,gs_hamkq,gvnlx1,idir,ipert,ikpt,isppol,&
&       mkmem,mpi_enreg,nband_k,nsppol,debug_mode,dtset%prtvol,rf_hamk_idir,size_cprj,size_wf,&
&       ffnl1=ffnl1,ffnl1_test=ffnl1_test)

       if (debug_mode/=0) then
         write(msg,'(a,i2)') 'RF2 TEST before accumulate_bands choice = 3 jband = ',jband
         call wrtout(std_out,msg)
       end if

!      For every occupied iband, we compute :
!      < u^(0)(iband) | d^2H/(dpert1 dpert2) | u^(0)(jband) > and add it to lambda_mn
!      < u^(0)(iband) | d^2S/(dpert1 dpert2) | u^(0)(jband) > and add it to amn
       do iband=1,rf2%nband_k  ! = band m
         if (abs(occ_k(iband))>tol8) then
           shift_band2=(iband-1)*size_wf
           cwave_i => cg(:,1+shift_band2+icg:size_wf+shift_band2+icg)
           if(ipert == natom+10) then
             ipert1 = natom+1
             ipert2 = natom+1
           else
             ipert1 = natom+1
             ipert2 = natom+2
           end if
           call rf2_getidirs(idir,idir1,idir2)
           call rf2_accumulate_bands(rf2,3,gs_hamkq,mpi_enreg,iband,idir1,idir2,ipert1,ipert2,&
           jband,debug_mode,cwave_i,h_cwave,s_cwave)
         end if
       end do

!      Add d^2H/(dk_dir1 dk_dir2)|u^(0)> to RHS_Stern :
       if (gs_hamkq%usepaw==1) h_cwave(:,:)=h_cwave(:,:)-eig0_k(jband)*s_cwave(:,:) ! if PAW : we add H^(2)-eps^(0) S^(2)
       rhs_j => rf2%RHS_Stern(:,1+shift_band1:size_wf+shift_band1)
       call cg_zaxpy(size_wf,(/one,zero/),h_cwave,rhs_j)

     end if ! empty band test
   end do ! jband
 end if ! H^(2) exists

!Computation of terms containing H^(1)
 do kdir1=1,rf2%ndir
!  First iteration (kdir1=1) :
!  pert1 = rf2%iperts(1) along rf2%idirs(1)
!  pert2 = rf2%iperts(2) along rf2%idirs(2)
!  Second iteration (kdir1=2) :
!  pert1 = rf2%iperts(2) along rf2%idirs(2)
!  pert2 = rf2%iperts(1) along rf2%idirs(1)
   shift_dir1=(kdir1-1)*nband_k*size_wf
   shift_cprj_dir1=(kdir1-1)*nband_k*size_cprj
   shift_dir1_lambda=(kdir1-1)*2*nband_k**2
   idir1=rf2%idirs(kdir1)
   ipert1=rf2%iperts(kdir1)
   if(ipert==natom+10 .and. idir<=3) then
     idir2=idir1
     ipert2=ipert1
     shift_dir2=0
     shift_dir2_lambda=0
     rf_hamk_idir => rf_hamkq
   else
     idir2=rf2%idirs(2-kdir1+1)
     ipert2=rf2%iperts(2-kdir1+1)
     shift_dir2=(2-kdir1)*nband_k*size_wf
     shift_dir2_lambda=(2-kdir1)*2*nband_k**2
     if (kdir1==1) rf_hamk_idir => rf_hamk_dir2 ! dir2
     if (kdir1==2) rf_hamk_idir => rf_hamkq ! dir1
   end if

!  Load projected WF according to ipert2 and idir2
   cprj_j => cprj_empty ;  ; cprj_dudk => cprj_empty
   if (has_cprj_jband) then
     call pawcprj_free(cprj_jband)
     ncpgr_loc= 3;if(ipert2==natom+1.or.ipert2==natom+2) ncpgr_loc=1
     icpgr_loc=-1;if(ipert2==natom+1.or.ipert2==natom+2) icpgr_loc=idir2
     call pawcprj_alloc(cprj_jband,ncpgr_loc,gs_hamkq%dimcprj)
     call pawcprj_get(gs_hamkq%atindx1,cprj_jband,cprj,natom,1,ibg,ikpt,iorder_cprj,&
&     isppol,dtset%mband,mkmem,natom,nband_k,nband_k,gs_hamkq%nspinor,nsppol,dtfil%unpaw,&
&     mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb,ncpgr=3,icpgr=icpgr_loc)
   end if

   do jband=1,nband_k

!    Skip bands not treated by current proc
     if(mpi_enreg%proc_distrb(ikpt,jband,isppol)/=me) cycle

     if (abs(occ_k(jband))>tol8) then
       shift_band1=(jband-1)*size_wf
       shift_cprj_band1=(jband-1)*size_cprj
       shift_jband_lambda=(jband-1)*2*nband_k

!      Extract first order wavefunction | du/dpert1 > and eigenvalues
       eig1_k_jband => eig1_k_stored(1+shift_jband_lambda+shift_dir2_lambda:2*nband_k+shift_jband_lambda+shift_dir2_lambda)
       cwave_dudk => dudk(:,1+shift_band1+shift_dir1:size_wf+shift_band1+shift_dir1)
       if (has_dudkprj) cprj_dudk => dudkprj(:,1+shift_cprj_band1+shift_cprj_dir1:size_cprj+shift_cprj_band1+shift_cprj_dir1)

       if (ipert2==natom+2) then
!        Extract dkdk and multiply by i :
         gvnlx1(1,:) = -dudkdk(2,1+shift_band1:size_wf+shift_band1)
         gvnlx1(2,:) =  dudkdk(1,1+shift_band1:size_wf+shift_band1)
       end if

!      Compute dH/dpert2 | du/dpert1 > (in h_cwave) and dS/dpert2 | du/dpert1 > (in s_cwave)
       call rf2_apply_hamiltonian(cg_jband,cprj_jband,cwave_dudk,cprj_dudk,h_cwave,s_cwave,&
&       eig0_k,eig1_k_jband,jband,gs_hamkq,gvnlx1,idir2,ipert2,ikpt,isppol,&
&       mkmem,mpi_enreg,nband_k,nsppol,debug_mode,dtset%prtvol,rf_hamk_idir,size_cprj,size_wf)

       if (debug_mode/=0) then
         write(msg,'(2(a,i2))') 'RF2 TEST before accumulate_bands choice = 4 kdir1 = ',kdir1,' jband = ',jband
         call wrtout(std_out,msg)
       end if

!      For every occupied iband, we compute :
!      < u^(0)(iband) | dH/dpert2 | du/dpert1(jband) > and add it to lambda_mn
!      < u^(0)(iband) | dS/dpert2 | du/dpert1(jband) > and add it to amn
       do iband=1,rf2%nband_k  ! = band m
         if (abs(occ_k(iband))>tol8) then
           shift_band2=(iband-1)*size_wf
           cwave_i => cg(:,1+shift_band2+icg:size_wf+shift_band2+icg)
           call rf2_accumulate_bands(rf2,4,gs_hamkq,mpi_enreg,iband,idir1,idir2,ipert1,ipert2,&
           jband,debug_mode,cwave_i,h_cwave,s_cwave)
         end if
       end do

!      Add dH/dpert2 | du/dpert1 > to RHS_Stern :
       if (gs_hamkq%usepaw==1) h_cwave(:,:)=h_cwave(:,:)-eig0_k(jband)*s_cwave(:,:) ! if PAW : we add H^(1)-eps^(0) S^(1)
       rhs_j => rf2%RHS_Stern(:,1+shift_band1:size_wf+shift_band1)
       call cg_zaxpy(size_wf,(/factor*one,zero/),h_cwave,rhs_j)

!      Compute : -factor * sum_iband ( dLambda/dpert1_{iband,jband} * dsusdu_{iband} )
       do iband=1,nband_k
         if (abs(occ_k(iband))>tol8) then ! if empty band, nothing to do

!          Extract lambda_ij(iband,jband) for dir1
           lambda_ij(1)=eig1_k_stored(2*iband-1+shift_jband_lambda+shift_dir1_lambda)
           lambda_ij(2)=eig1_k_stored(2*iband  +shift_jband_lambda+shift_dir1_lambda)

!          Extract dsusdu for iband and pert2 (in cwave_i)
           shift_band2=(iband-1)*size_wf
           cwave_i => dsusdu(:,1+shift_band2+shift_dir2:size_wf+shift_band2+shift_dir2)

!          Compute Lambda_{iband,jband} * dsusdu_{iband} and add it to RHS_Stern
           rhs_j => rf2%RHS_Stern(:,1+shift_band1:size_wf+shift_band1)
           call cg_zaxpy(size_wf,(/-factor*lambda_ij(1),-factor*lambda_ij(2)/),cwave_i,rhs_j) !do not forget the minus sign!

         end if ! empty iband test
       end do ! iband

     end if ! empty jband test
   end do ! jband

 end do ! kdir1

 ABI_DEALLOCATE(gvnlx1)
 ABI_DEALLOCATE(h_cwave)
 ABI_DEALLOCATE(s_cwave)
 ABI_DEALLOCATE(cg_jband)
 ABI_DEALLOCATE(dudk)
 ABI_DEALLOCATE(dudkdk)
 ABI_DEALLOCATE(dudk_dir2)
 ABI_DEALLOCATE(dsusdu)
 ABI_DEALLOCATE(eig1_k_stored)
 if (has_cprj_jband) call pawcprj_free(cprj_jband)
 ABI_DATATYPE_DEALLOCATE(cprj_jband)
 if (has_dudkprj) call pawcprj_free(dudkprj)
 ABI_DATATYPE_DEALLOCATE(dudkprj)

! Compute the part of 2nd order wavefunction that belongs to the space of empty bands
 do jband=1,nband_k

!  Skip bands not treated by current proc
   if(mpi_enreg%proc_distrb(ikpt,jband,isppol)/=me) cycle

   shift_band1=(jband-1)*size_wf
   if (abs(occ_k(jband))>tol8) then
     invocc = one/occ_k(jband)
     rhs_j => rf2%RHS_Stern(:,1+shift_band1:size_wf+shift_band1)
     do iband=1,nband_k
       if (iband /= jband) then
         if (debug_mode/=0) then
           if (abs(occ_k(iband) - occ_k(jband)) > tol12 .and. occ_k(iband) > tol8) then
             write(msg,'(a,i2,a,i2)') 'RF2 TEST ACTIVE SPACE : jband = ',jband,' iband = ',iband
             call wrtout(std_out,msg)
             call wrtout(ab_out,msg)
             write(msg,'(a)') 'ERROR : occ_k(iband) /= occ_k(jband) (and both are >0)'
             call wrtout(std_out,msg)
             call wrtout(ab_out,msg)
           end if
           if (abs(eig0_k(iband) - eig0_k(jband)) < tol8 ) then
             write(msg,'(a,i2,a,i2)') 'RF2 TEST ACTIVE SPACE : jband = ',jband,' iband = ',iband
             call wrtout(std_out,msg)
             write(msg,'(a,es22.13e3)') 'WARNING : DEGENERATE BANDS  Eig(jband) = Eig(jband) = ',eig0_k(jband)
             call wrtout(std_out,msg)
           end if
           if ( (eig0_k(iband) - eig0_k(jband) < -tol12) .and. (jband < iband) ) then
             write(msg,'(a,i2,a,i2)') 'RF2 TEST ACTIVE SPACE : jband = ',jband,' iband = ',iband
             call wrtout(std_out,msg)
             call wrtout(ab_out,msg)
             write(msg,'(a)') 'ERROR : Eig(jband) < Eig(iband) with jband < iband'
             call wrtout(std_out,msg)
             call wrtout(ab_out,msg)
             write(msg,'(a,es22.13e3)') 'Eig(jband) = ',eig0_k(jband)
             call wrtout(std_out,msg)
             call wrtout(ab_out,msg)
             write(msg,'(a,es22.13e3)') 'Eig(iband) = ',eig0_k(iband)
             call wrtout(std_out,msg)
             call wrtout(ab_out,msg)
           end if
         end if ! end tests
         if ( abs(occ_k(iband))<tol8 ) then ! for empty bands only
           shift_band2=(iband-1)*size_wf
           cwave_i => cg(:,1+shift_band2+icg:size_wf+shift_band2+icg)
           call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cwave_i,rhs_j,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
!          Store it in a_mn
!          /!\ There is a factor "-2" to simplify the use of amn in the following.
!          /!\ Occupied and empty bands will be treated in a same way.
           rf2%amn(:,iband+(jband-1)*nband_k)=-two*rocceig(iband,jband)*invocc*(/dotr,doti/)&
           +rf2%amn(:,iband+(jband-1)*nband_k)
         end if ! empty band test
       end if ! iband \= jband
     end do ! iband
   end if  ! empty band test
 end do ! jband

! **************************************************************************************************
!  COMPUTATION OF "dcwavef" AND "Lambda_mn" FROM "A_mn"
! **************************************************************************************************

 ABI_MALLOC_OR_DIE(rf2%dcwavef,(2,nband_k*size_wf), ierr)
 rf2%dcwavef=zero

 do jband=1,nband_k

!  Skip bands not treated by current proc
   if(mpi_enreg%proc_distrb(ikpt,jband,isppol)/=me) cycle

   shift_band1=(jband-1)*size_wf
   if (abs(occ_k(jband))>tol8) then
     do iband=1,nband_k
       shift_band2=(iband-1)*size_wf

!      Extract GS wavefunction for iband
       cwave_i => cg(:,1+shift_band2+icg:size_wf+shift_band2+icg)

       call cg_zaxpy(size_wf,-half*rf2%amn(:,iband+(jband-1)*nband_k), &
&       cwave_i,rf2%dcwavef(:,1+shift_band1))

       if (abs(occ_k(iband))>tol8 .and. abs(occ_k(jband))>tol8) then
         rf2%lambda_mn(:,iband+(jband-1)*nband_k) = rf2%lambda_mn(:,iband+(jband-1)*nband_k) &
         -half*(eig0_k(iband)+eig0_k(jband))*rf2%amn(:,iband+(jband-1)*nband_k)

         eig1_k(2*iband-1+(jband-1)*2*nband_k) = rf2%lambda_mn(1,iband+(jband-1)*nband_k)
         eig1_k(2*iband  +(jband-1)*2*nband_k) = rf2%lambda_mn(2,iband+(jband-1)*nband_k)

       end if ! empty band test
     end do ! iband
   end if ! empty band test
 end do ! jband

! For the following, "rf2%lambda_mn" and "rf2%RHS_Stern" must be computed for every bands
 call xmpi_barrier(mpi_enreg%comm_band)

! **************************************************************************************************
!  FINAL TEST
! **************************************************************************************************

 tol_final = tol6
 if (debug_mode/=0) then
   do jband=1,nband_k

!    Skip bands not treated by current proc
     if(mpi_enreg%proc_distrb(ikpt,jband,isppol)/=me) cycle

     if (abs(occ_k(jband))>tol8) then
!       write(msg,'(3(a,i2))') 'RF2 TEST FINAL : ipert=',ipert-natom,' idir=',idir,' jband=',jband
!       call wrtout(std_out,msg)
       shift_band1=(jband-1)*size_wf
       rhs_j => rf2%RHS_Stern(:,1+shift_band1:size_wf+shift_band1)
       cwave_j => cg(:,1+shift_band1+icg:size_wf+shift_band1+icg)
       call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cwave_j,rhs_j,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
       dot2r = dotr - rf2%lambda_mn(1,jband+(jband-1)*nband_k)
       dot2i = doti - rf2%lambda_mn(2,jband+(jband-1)*nband_k)
       dot2r = sqrt(dot2r**2+dot2i**2)
       if (dot2r > tol_final) then
         write(msg,'(a,i2,a,es22.13E3)') 'RF2 TEST FINAL iband = ',jband,' : NOT PASSED dotr = ',dotr
         call wrtout(std_out,msg)
         call wrtout(ab_out,msg)
         write(msg,'(2(a,es22.13E3))') ' < cwave_j | rhs_j > =',dotr,',',doti
         call wrtout(std_out,msg)
         call wrtout(ab_out,msg)
         write(msg,'(2(a,es22.13E3))') '           lambda_jj =',&
&         rf2%lambda_mn(1,jband+(jband-1)*nband_k),',',rf2%lambda_mn(2,jband+(jband-1)*nband_k)
         call wrtout(std_out,msg)
         call wrtout(ab_out,msg)
       else
         if (dot2r<tol9) dot2r = zero ! in order to hide the numerical noise
         write(msg,'(a,i2,a,es22.13E3,a,es7.1E2)') &
         'RF2 TEST FINAL iband = ',jband,' : OK. |test| = ',dot2r,' < ',tol_final
         call wrtout(std_out,msg)
       end if
     end if
   end do
 end if

! **************************************************************************************************
!  JOB FINISHED
! **************************************************************************************************

! Deallocations of arrays
 if (debug_mode==0) then
   ABI_DEALLOCATE(rf2%amn)
 end if

 call timab(514,2,tsec)

 DBG_EXIT("COLL")

end subroutine rf2_init
!!***

end module m_rf2_init
!!***
