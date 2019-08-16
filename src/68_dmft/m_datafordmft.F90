!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_datafordmft
!! NAME
!!  m_datafordmft
!!
!! FUNCTION
!! This module produces inputs for the DMFT calculation
!!
!! COPYRIGHT
!! Copyright (C) 2006-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

MODULE m_datafordmft

 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_abicore
 use m_errors
 use m_xmpi
 use m_dtset

 use defs_datatypes, only : pseudopotential_type
 use m_io_tools,  only : open_file
 use m_crystal, only : crystal_t
 use m_matlu, only: matlu_type,init_matlu,sym_matlu,copy_matlu,print_matlu,diff_matlu,destroy_matlu, checkdiag_matlu, &
                    gather_matlu,add_matlu
 use m_oper, only : init_oper,oper_type,identity_oper,loc_oper,destroy_oper,diff_oper, upfold_oper,copy_oper,prod_oper
 use m_pawang, only : pawang_type
 use m_pawtab, only : pawtab_type
 use m_paw_ij, only : paw_ij_type
 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_free
 use m_paw_dmft, only: paw_dmft_type
 use m_mpinfo,   only : proc_distrb_cycle
 use m_matrix, only : invsqrt_matrix

 implicit none

 private

 public :: datafordmft
 public :: compute_levels
 public :: psichi_renormalization
 public :: hybridization_asymptotic_coefficient
!!***

contains

!!****f* ABINIT/datafordmft
!! NAME
!! datafordmft
!!
!! FUNCTION
!!  Compute psichi (and print some data for check)
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!        -gprimd(3,3)=dimensional reciprocal space primitive translations
!!        -indsym(4,nsym,natom)=indirect indexing array for atom labels
!!        -symrec(3,3,nsym)=symmetry operations in reciprocal space
!!        - nsym= number of symetry operations
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                          and each |p_lmn> non-local projector
!!  dimcprj(natom) = dimension for cprj
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  fermie= Fermi energy
!!  lda_occup <type(oper_type)> = occupations in the correlated orbitals in LDA
!!  mband=maximum number of bands
!!  mkmem =number of k points treated by this node
!!  mpi_enreg=information about MPI parallelization
!!  nkpt=number of k points.
!!  my_nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol) = occupancies of KS states.
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  unpaw = file number for cprj
!!
!! OUTPUT
!!  paw_dmft%psichi(nsppol,nkpt,mband,my_nspinor,dtset%natom,(2*maxlpawu+1))): projections <Psi|chi>
!!  paw_dmft%eigen(paw_dmft%nsppol,paw_dmft%nkpt,paw_dmft%mband)
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv,vtorho
!!
!! CHILDREN
!!
!! SOURCE

subroutine datafordmft(cryst_struc,cprj,dimcprj,dtset,eigen,fermie,&
& lda_occup,mband,mband_cprj,mkmem,mpi_enreg,nkpt,my_nspinor,nsppol,occ,&
& paw_dmft,paw_ij,pawang,pawtab,psps,usecprj,unpaw,nbandkss)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mband_cprj,mkmem
 integer,intent(in) :: nkpt,my_nspinor,nsppol
 integer,intent(in) :: unpaw,usecprj
 integer, optional, intent(in) :: nbandkss
 real(dp),intent(in) :: fermie
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(oper_type), intent(inout) :: lda_occup !vz_i
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(crystal_t),intent(in) :: cryst_struc
!arrays
 integer, intent(in) :: dimcprj(cryst_struc%natom)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol)
 type(paw_ij_type),intent(in) :: paw_ij(:)
! type(pawcprj_type) :: cprj(cryst_struc%natom,my_nspinor*mband*mkmem*nsppol)
 type(pawcprj_type), intent(in) :: cprj(cryst_struc%natom,my_nspinor*mband_cprj*mkmem*nsppol*usecprj)
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
!Local variables-------------------------------
!scalars
 integer :: band_index,dimpsichi,facpara
 integer :: iat,iatom,ib,iband,ibandc,ibg,icat,icount_proj_ilmn,idijeff,ierr,ierrr,ikpt
 integer :: ilmn,im,im1,iorder_cprj,ispinor,ispinor1,isppol,itypat,ilmn1
 integer :: jj1,ldim,lmn_size,lpawu
 integer :: m1,m1_x2my2d,m1_t2g,m1_t2g_mod,maxnproju,me,natom,nband_k,nband_k_cprj
 integer :: nbandi,nbandf,nnn,nprocband,nsploop,option,opt_renorm,spaceComm,unt
 real(dp) :: ph0phiint_used
 character(len=500) :: message
!arrays
 real(dp) :: chinorm
 complex(dpc), allocatable :: buffer1(:)
 type(matlu_type), allocatable :: loc_occ_check(:)
 type(matlu_type), allocatable :: loc_norm_check(:)
 type(matlu_type), allocatable :: xocc_check(:)
 type(matlu_type), allocatable :: xnorm_check(:)
 type(matlu_type), allocatable :: matlu_temp(:)
 logical :: lprojchi,t2g,x2my2d
 integer,parameter :: spinor_idxs(2,4)=RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
 type(pawcprj_type),allocatable :: cwaveprj(:,:)
!************************************************************************

!DBG_ENTER("COLL")
!Fake test to keep fermie as argument. REMOVE IT AS SOON AS POSSIBLE ...
 if(fermie>huge(zero))chinorm=zero

 facpara=1 !mpi_enreg%nproc
 if(abs(dtset%pawprtvol)>=3) then
   write(message,*) " number of k-points used is nkpt=nkpt ", nkpt
   call wrtout(std_out,  message,'COLL')
   write(message,*) " warning: parallelised version        ", nkpt
   call wrtout(std_out,  message,'COLL')
   write(message,*) " weights k-points used is wtk=wtk"
   call wrtout(std_out,  message,'COLL')
 end if

 if(usecprj==0) then
   write(message,*) "  usecprj=0 : BUG in datafordmft",usecprj
   MSG_BUG(message)
 end if

 if(my_nspinor/=dtset%nspinor) then
   write(message,*) "  my_nspinor=/dtset%nspinor, datafordmft not working in this case",&
&   my_nspinor,dtset%nspinor
   MSG_ERROR(message)
 end if


!do ib=1,my_nspinor*mband_cprj*mkmem*nsppol*usecprj
!write(std_out,'(a,i6,3e16.7)') "cprj",ib,cprj(1,ib)%cp(1,19),cprj(1,ib)%cp(2,19),cprj(1,ib)%cp(1,19)**2+cprj(1,ib)%cp(2,19)**2
!enddo

!----------------------------------- MPI-------------------------------------

!Init parallelism
 spaceComm=mpi_enreg%comm_cell
 if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_kpt
 me=mpi_enreg%me_kpt
 iorder_cprj=0
 ABI_CHECK(dtset%mkmem/=0,"mkmem==0 not supported anymore!")
!todo_ab: extract cprj from file unpaw in the following..
!call abi_abort('COLL')

!----------------------------------- MPI-------------------------------------

 nbandi=paw_dmft%dmftbandi
 nbandf=paw_dmft%dmftbandf
 lprojchi=.false.
 lprojchi=.true.
 t2g=(paw_dmft%dmftqmc_t2g==1)
 x2my2d=(paw_dmft%dmftqmc_x2my2d==1)
 natom=cryst_struc%natom

!if(mpi_enreg%me==0) write(7886,*) "in datafordmft", mpi_enreg%me, mpi_enreg%nproc
!if(mpi_enreg%me==1) write(7887,*) "in datafordmft", mpi_enreg%me, mpi_enreg%nproc
!if(mpi_enreg%me==2) write(7888,*) "in datafordmft", mpi_enreg%me, mpi_enreg%nproc
 write(message,'(2a)') ch10,&
& '  == Prepare data for DMFT calculation  '
 call wrtout(std_out,message,'COLL')
 if(abs(dtset%pawprtvol)>=3) then
   write(message, '(a,a)' ) ch10,&
&   '---------------------------------------------------------------'
!  call wrtout(ab_out,message,'COLL')
!  call wrtout(std_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,a,a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&   '  Print useful data (as a check)',ch10,&
&   '  - Overlap of KS wfc with atomic orbital inside sphere',ch10,&
&   '  - Eigenvalues',ch10,&
&   '  - Weights of k-points',ch10,&
&   '  - Number of spins ',ch10,&
&   '  - Number of states'
!  call wrtout(ab_out,message,'COLL')
!  call wrtout(std_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,a)' ) ch10,&
&   '---------------------------------------------------------------'
 end if
 if(dtset%nstep==0.and.dtset%nbandkss==0) then
   message = 'nstep should be greater than 1'
   MSG_BUG(message)
 end if

!********************* Max Values for U terms.
!maxlpawu=0
 maxnproju=0
 do iatom=1,natom
   if(pawtab(dtset%typat(iatom))%lpawu.ne.-1 .and. pawtab(dtset%typat(iatom))%nproju.gt.maxnproju)&
&   maxnproju=pawtab(dtset%typat(iatom))%nproju
 end do
!*****************   in forlb.eig
 if(me.eq.0.and.abs(dtset%pawprtvol)>=3) then
   if (open_file('forlb.eig',message,newunit=unt,form='formatted',status='unknown') /= 0) then
     MSG_ERROR(message)
   end if
   rewind(unt)
   write(unt,*) "Number of bands,   spins, and  k-point; and spin-orbit flag"
   write(unt,*) mband,nsppol,nkpt,my_nspinor,nbandi,nbandf
   write(unt,*) " For each k-point, eigenvalues for each band"
   write(unt,*) (dtset%wtk(ikpt)*facpara,ikpt=1,nkpt)
   band_index=0
   do isppol=1,nsppol
     write(unt,*) " For spin"
     write(unt,*)  isppol
     do ikpt=1,nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
       ibandc=0
       write(unt,*) " For k-point"
       write(unt,*)  ikpt
       do iband=1,mband
         if(paw_dmft%band_in(iband)) then
           ibandc=ibandc+1
           write(unt, '(2i6,4x,f20.15)' ) ibandc,ikpt,eigen(iband+band_index)*2.d0
         end if
       end do
       band_index=band_index+nband_k
     end do
   end do
   close(unt)
 end if ! proc=me

!==   put eigen into eigen_lda
 band_index=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
     ibandc=0
     do iband=1,mband
       if(paw_dmft%band_in(iband)) then
         ibandc=ibandc+1
         paw_dmft%eigen_lda(isppol,ikpt,ibandc)=eigen(iband+band_index) ! in Ha
        ! paw_dmft%eigen_lda(isppol,ikpt,ibandc)=fermie
       end if
     end do
     band_index=band_index+nband_k
   end do
 end do

 if(abs(dtset%pawprtvol)>=3) then
   write(message, '(a,a)' ) ch10,&
&   '   datafordmft :  eigenvalues written'
   call wrtout(std_out,  message,'COLL')
 end if
!==========================================================================
!***************** Compute  <Psi|Chi>=\sum_{proja} <Psi|P_a><phi_a|Chi>
!==========================================================================
!write(std_out,*) "size(cprj,dim=1)",size(cprj,dim=1),size(cprj,dim=2),dtset%mband,dtset%mkmem,dtset%nkpt

!Allocate temporary cwaveprj storage
 ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,my_nspinor))
!write(std_out,*) "before alloc cprj"
!write(std_out,*) size(cwaveprj,dim=1),size(cwaveprj,dim=2),size(dimcprj,dim=1)

 call pawcprj_alloc(cwaveprj,0,dimcprj)
!write(std_out,*) "after alloc cprj"

 nprocband=(mband/mband_cprj)
 ibg=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
     nband_k_cprj=nband_k/nprocband

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle

!    write(2011,*) ikpt
!    ibsp=ibg
     ibandc=0
!    LOOP OVER BANDS
     ib=0
     do iband=1,nband_k

       if(paw_dmft%band_in(iband)) ibandc=ibandc+1

!      Parallelization: treat only some bands
       if (dtset%paral_kgb==1) then
         if (mod((iband-1)/mpi_enreg%bandpp,mpi_enreg%nproc_band)/=mpi_enreg%me_band) cycle
       else
         if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me) cycle
       end if
       ib=ib+1

       do ispinor=1,my_nspinor
!        ibsp =~ (ikpt-1)*nband*my_nspinor+iband
!        ibsp=ibsp+1
         icat=1
!        write(std_out,*) isppol,ikpt,iband,ispinor
         iat=0 ! to declare
         do itypat=1,dtset%ntypat
           lmn_size=pawtab(itypat)%lmn_size
!          write(std_out,*) isppol,ikpt,iband,ispinor
           do iatom=icat,icat+cryst_struc%nattyp(itypat)-1
             lpawu=pawtab(itypat)%lpawu
!            ----------   t2g case
             if(paw_dmft%dmftqmc_t2g==1.and.lpawu/=-1) then
               if(lpawu==2) then
!                lpawu==2 must be chosen in input and thus in
!                pawtab. On the contrary, paw_dmft now has
!                lpawu=1
                 m1_t2g=0 ! index for psichi which has a dimension 3
               else
                 write(message,'(a,a,i4,i4,2a)')  ch10,&
&                 '  Wrong use of dmftqmc_t2g',paw_dmft%dmftqmc_t2g,lpawu,ch10,&
&                 ' Action: desactivate qmftqmc_t2g or use lpawu=2'
                 MSG_ERROR(message)
               end if
             end if
!            ----------   t2g case
!            ----------   x2my2d case
             if(paw_dmft%dmftqmc_x2my2d==1.and.lpawu/=-1) then
               if(lpawu==2) then
!                lpawu==2 must be chosen in input and thus in
!                pawtab. On the contrary, paw_dmft now has
!                lpawu=1
                 m1_x2my2d=0 ! index for psichi which has a dimension 1
               else
                 write(message,'(a,a,i4,i4,2a)')  ch10,&
&                 '  Wrong use of dmftqmc_x2my2d',paw_dmft%dmftqmc_x2my2d,lpawu,ch10,&
&                 ' Action: desactivate dmftqmc_x2my2d or use lpawu=2'
                 MSG_ERROR(message)
               end if
             end if
!            ----------   x2my2d case
!            if(isppol==2) write(std_out,*) "ee",size(cprj(iatom,ibsp)%cp(:,:))

             iat=iat+1
             jj1=0
             if(lpawu.ne.-1) then

               call pawcprj_get(cryst_struc%atindx1,cwaveprj,cprj,natom,ib,ibg,ikpt,&
&               iorder_cprj,isppol,mband_cprj,mkmem,natom,1,nband_k_cprj,&
&               my_nspinor,nsppol,unpaw,mpicomm=mpi_enreg%comm_kpt,&
&               proc_distrb=mpi_enreg%proc_distrb)
!              write(std_out,*) "M-2",isppol,ikpt,iband,iatom,&
!              &             (cwaveprj(iatom,ispinor)%cp(1,13)**2+cwaveprj(iatom,ispinor)%cp(2,13)**2) ,ibsp

!              chinorm=(pawtab(itypat)%phiphjint(1))
               chinorm=1.d0
!              write(std_out,*) isppol,ikpt,iband,ispinor,iat
               do ilmn=1,lmn_size
!                write(std_out,*) ilmn
!                ------------ Select l=lpawu.
                 if (psps%indlmn(1,ilmn,itypat)==lpawu) then
!                  ------------ Check that the band is choosen (within nbandi and nbandf)
                   if(paw_dmft%band_in(iband)) then
!                    if(ilmn==13) then
!                    write(std_out,*) "M-2c",isppol,ikpt,iband,iatom
!                    write(std_out,*) "M-2b",isppol,ikpt,ibandc,&
!                    &             (cwaveprj(iatom,ispinor)%cp(1,13)**2+cwaveprj(iatom,ispinor)%cp(2,13)**2)
!                    endif
!                    write(std_out,*) "inside paw_dmft%band_in",iband
                     jj1=jj1+1
                     if(jj1>pawtab(itypat)%nproju*(2*lpawu+1)) then
                       write(message,'(a,a,a,a)')  ch10,&
&                       ' jj1 is not correct in datafordmft',ch10,&
&                       ' Action: CONTACT Abinit group'
                       MSG_BUG(message)
                     end if ! jj1
                     icount_proj_ilmn=psps%indlmn(3,ilmn,itypat)  ! n: nb of projector
!                    write(std_out,*) "icount_proj_ilmn",icount_proj_ilmn
                     m1=psps%indlmn(2,ilmn,itypat)+psps%indlmn(1,ilmn,itypat)+1
!                    ---- if lprochi=true do the sum over every projectors
!                    ---- if lprochi=false: . do the sum over only ground state projector
!                    ----                   . and take ph0phiint(icount_proj_ilmn)=1

!                    call pawcprj_get(cryst_struc%atindx1,cwaveprj,cprj,natom,&
!                    &             ib,ibg,ikpt,iorder_cprj,isppol,mband_cprj,mkmem,&
!                    &             natom,1,nband_k_cprj,my_nspinor,nsppol,unpaw,&
!                    &             mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)

                     if(lprojchi.or.icount_proj_ilmn==1) then
                       if(.not.lprojchi) ph0phiint_used=one
                       if(lprojchi) ph0phiint_used=pawtab(itypat)%ph0phiint(icount_proj_ilmn)
                       if(paw_dmft%dmftqmc_t2g==1) then ! t2g case

                         if(m1==1.or.m1==2.or.m1==4) then ! t2g case
                           m1_t2g=m1_t2g+1  ! t2g case1
                           m1_t2g_mod=mod(m1_t2g-1,3)+1
!                          write(std_out,*) "M0",isppol,ikpt,iband,ilmn,cprj(iatom,ibsp)%cp(1,ilmn)
!                          write(std_out,*) "M1",m1,m1_t2g,iband,ilmn,icount_proj_ilmn
                           paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1_t2g_mod)=&  ! t2g case
&                          paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1_t2g_mod)+&  ! t2g case
!                          &                           cmplx(cprj(iatom,ibsp)%cp(1,ilmn)*ph0phiint_used,&  ! t2g case
!                          &                           cprj(iatom,ibsp)%cp(2,ilmn)*ph0phiint_used,kind=dp)  ! t2g case
&                          cmplx(cwaveprj(iatom,ispinor)%cp(1,ilmn)*ph0phiint_used,&  ! t2g case
&                          cwaveprj(iatom,ispinor)%cp(2,ilmn)*ph0phiint_used,kind=dp)  ! t2g case
                         end if  !t2g case
                       else if(paw_dmft%dmftqmc_x2my2d==1) then ! x2my2d case
                         if(m1==5) then ! x2my2d case
                           m1_x2my2d=1  ! x2my2d case1
                           paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1_x2my2d)=&      ! x2my2d case
&                          paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1_x2my2d)+&      ! x2my2d case
&                          cmplx(cwaveprj(iatom,ispinor)%cp(1,ilmn)*ph0phiint_used,&       ! x2my2d case
&                          cwaveprj(iatom,ispinor)%cp(2,ilmn)*ph0phiint_used,kind=dp)      ! x2my2d case
                         end if  !x2my2d case
                       else
                         paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1)=&
&                         paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1)+&
!                        &                         cmplx(cprj(iatom,ibsp)%cp(1,ilmn)*ph0phiint_used,&
!                        &                         cprj(iatom,ibsp)%cp(2,ilmn)*ph0phiint_used,kind=dp)
&                         cmplx(cwaveprj(iatom,ispinor)%cp(1,ilmn)*ph0phiint_used,&
&                         cwaveprj(iatom,ispinor)%cp(2,ilmn)*ph0phiint_used,kind=dp)

!                        if(ibandc==3.and.iat==1.and.m1==1) then
!                        write(std_out,'(a,3i5)') "psichi integers",isppol,ikpt,ispinor
!                        write(std_out,'(a,2i5,2e16.7)') "psichi IB3 iAT1 IM1",ilmn,icount_proj_ilmn,&
!                        &             real(paw_dmft%psichi(isppol,ikpt,3,ispinor,1,1)), imag(paw_dmft%psichi(isppol,ikpt,3,ispinor,1,1))
!                        write(std_out,'(a,2i5,2e16.7)') "cwaveprj IB3 iAT1 IM1",ilmn,icount_proj_ilmn,cwaveprj(iatom,ispinor)%cp(1,ilmn) &
!                        &                    , cwaveprj(iatom,ispinor)%cp(2,ilmn)
!                        endif
                       end if
                     end if ! lprojchi=.true. (always)
                   end if ! nband belongs to paw_dmft%band_in
                 end if ! L=lpawu
               end do !ilmn : loop over L,M,N
             end if ! If lpawu/=-1
           end do ! iatom
           icat=icat+cryst_struc%nattyp(itypat)
         end do ! itypat
       end do ! ispinor
     end do !iband
     ibg=ibg+nband_k_cprj*my_nspinor
!    bdtot_index=bdtot_index+nband_k ! useless  (only for occ)
   end do !ikpt
 end do ! isppol
!do isppol=1,nsppol
!do ikpt=1,nkpt
!do ispinor=1,my_nspinor
!write(std_out,*) "psichi integers",isppol,ikpt,ispinor
!write(std_out,*) "psichi IB3 iAT1 IM1",&
!&             real(paw_dmft%psichi(isppol,ikpt,3,ispinor,1,1)), imag(paw_dmft%psichi(isppol,ikpt,3,ispinor,1,1))
!
!enddo
!enddo
!enddo
!call abi_abort('COLL')
 if(abs(dtset%pawprtvol)>=3) then
   write(message,*) "chinorm used here =",chinorm
   call wrtout(std_out,  message,'COLL')
 end if

!deallocate temporary cwaveprj/cprj storage
 call pawcprj_free(cwaveprj)
 ABI_DATATYPE_DEALLOCATE(cwaveprj)

!==========================================================================
!********************* Gather information for MPI before printing
!==========================================================================

 dimpsichi=2*nsppol*nkpt*mband*my_nspinor*natom*(2*paw_dmft%maxlpawu+1)
 ABI_ALLOCATE(buffer1,(dimpsichi))
 buffer1 = zero
 nnn=0
!write(176,*) "beg",psichi
 do isppol=1,nsppol
   do ikpt=1,nkpt
     do ibandc=1,paw_dmft%mbandc
       do ispinor=1,my_nspinor
         do iat=1,natom
           do m1=1,2*paw_dmft%maxlpawu+1
!            do m=1,2
             nnn=nnn+1
             buffer1(nnn)=paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1)
!            enddo
           end do
         end do
       end do
     end do
   end do
 end do
 call xmpi_barrier(spaceComm)
 call xmpi_sum(buffer1,spaceComm ,ierr)
 if (dtset%paral_kgb==1.and.nprocband>1) then
   call xmpi_sum(buffer1,mpi_enreg%comm_band,ierr) !Build sum over band processors
 end if
 call xmpi_barrier(spaceComm)
 nnn=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     do ibandc=1,paw_dmft%mbandc
       do ispinor=1,my_nspinor
         do iat=1,natom
           do m1=1,2*paw_dmft%maxlpawu+1
!            do m=1,2
             nnn=nnn+1
             paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1)=buffer1(nnn)
!            enddo
             ! if(ibandc==1) then
             !   write(6,*) "m1",m1
             !   write(6,*) "psichi datafordmft",paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1)
             ! endif
           end do
         end do
       end do
     end do
   end do
 end do
 ABI_DEALLOCATE(buffer1)
!write(177,*) "end",psichi

!do isppol=1,nsppol
!do ikpt=1,nkpt
!do ibandc=1,paw_dmft%mbandc
!do ispinor=1,my_nspinor
!write(std_out,*) "psichigather",isppol,ikpt,ibandc,&
!&             real(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,1,1))**2+&
!&             imag(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,1,1))**2
!
!enddo
!enddo
!enddo
!enddo

 call xmpi_barrier(spaceComm)
!if(mpi_enreg%me.eq.0) write(177,*) "end",psichi
!if(mpi_enreg%me.eq.1) write(178,*) "end",psichi
!if(mpi_enreg%me.eq.2) write(179,*) "end",psichi
 call xmpi_barrier(spaceComm)
!==========================================================================
!********* WRITE psichi in file for reference
!==========================================================================
 if(me.eq.0) then
   call psichi_print(dtset,cryst_struc%nattyp,cryst_struc%ntypat,nkpt,my_nspinor,&
&   nsppol,paw_dmft,pawtab,psps,t2g,x2my2d)
 end if ! proc=me

!********************* Check normalization  and occupations ***************
! Only if the complete BZ is sampled (ie paw_dmft%kspectralfunc=0)
!==========================================================================
 if(paw_dmft%dmft_kspectralfunc==0) then
   ABI_DATATYPE_ALLOCATE(xocc_check,(natom))
   ABI_DATATYPE_ALLOCATE(xnorm_check,(natom))
   call init_matlu(natom,my_nspinor,nsppol,paw_dmft%lpawu,xocc_check)
   call init_matlu(natom,my_nspinor,nsppol,paw_dmft%lpawu,xnorm_check)
   call psichi_check(dtset,cryst_struc%nattyp,nkpt,my_nspinor,&
&   nsppol,cryst_struc%ntypat,paw_dmft,pawtab,psps,xocc_check,xnorm_check)
!==========================================================================
!***************  write checks  *******************************************
!==========================================================================
   if(abs(dtset%pawprtvol)>=3) then
     write(message,*) "normalization computed"
     call wrtout(std_out,  message,'COLL')
   end if

   ABI_DATATYPE_ALLOCATE(loc_occ_check,(natom))
   ABI_DATATYPE_ALLOCATE(loc_norm_check,(natom))
   call init_matlu(natom,my_nspinor,nsppol,paw_dmft%lpawu,loc_occ_check)
   call init_matlu(natom,my_nspinor,nsppol,paw_dmft%lpawu,loc_norm_check)
   call copy_matlu(xocc_check,loc_occ_check,natom)
   call copy_matlu(xnorm_check,loc_norm_check,natom)

   write(message,'(2a,i4)')  ch10," == Check: Occupations and Norm from psichi are"
   call wrtout(std_out,  message,'COLL')

   if(paw_dmft%dmftcheck>=1) then
!    print occupations
     write(message,'(2a,i4)')  ch10,'  ------ Unsymetrised Occupation'
     call wrtout(std_out,  message,'COLL')

     call print_matlu(xocc_check,natom,dtset%pawprtvol)

!    print norms
     write(message,'(2a,i4)')  ch10,'  ------ Unsymetrised Norm'
     call wrtout(std_out,  message,'COLL')

     call print_matlu(xnorm_check,natom,dtset%pawprtvol)
   end if

!  symetrise and print occupations
   call sym_matlu(cryst_struc,loc_occ_check,pawang,paw_dmft)

   write(message,'(2a,i4)')  ch10,'  ------ Symetrised Occupation'
   call wrtout(std_out,  message,'COLL')

   call print_matlu(loc_occ_check,natom,dtset%pawprtvol)

!  symetrise and print norms
   call sym_matlu(cryst_struc,loc_norm_check,pawang,paw_dmft)

   write(message,'(2a,i4)')  ch10,'  ------ Symetrised Norm'
   call wrtout(std_out,  message,'COLL')

   call print_matlu(loc_norm_check,natom,dtset%pawprtvol)

!  deallocations
   do iatom=1,natom
     lda_occup%matlu(iatom)%mat=loc_occ_check(iatom)%mat
   end do

!  Tests density matrix LDA+U and density matrix computed here.
   if(paw_dmft%dmftcheck==2.or.(paw_dmft%dmftbandi==1)) then
     ABI_DATATYPE_ALLOCATE(matlu_temp,(natom))
     call init_matlu(natom,paw_dmft%nspinor,paw_dmft%nsppol,paw_dmft%lpawu,matlu_temp)
     do iatom=1,natom
       if(paw_dmft%lpawu(iatom).ne.-1) then
         ldim=2*paw_dmft%lpawu(iatom)+1
         nsploop=max(paw_dmft%nsppol,paw_dmft%nspinor**2)
         do idijeff=1,nsploop
           ispinor=0
           ispinor1=0
           if(nsploop==2) then
             isppol=spinor_idxs(1,idijeff)
             ispinor=1
             ispinor1=1
           else if(nsploop==4) then
             isppol=1
             ispinor=spinor_idxs(1,idijeff)
             ispinor1=spinor_idxs(2,idijeff)
           else if(nsploop==1) then
             isppol=1
             ispinor=1
             ispinor1=1
           else
             write(message,'(2a)') " BUG in datafordmft: nsploop should be equal to 2 or 4"
             call wrtout(std_out,message,'COLL')
           end if
           do im1 = 1 , ldim
             do im = 1 ,  ldim
               if(my_nspinor==2) matlu_temp(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=&
&               cmplx(paw_ij(iatom)%noccmmp(1,im,im1,idijeff),paw_ij(iatom)%noccmmp(2,im,im1,idijeff),kind=dp)
               if(my_nspinor==1) matlu_temp(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=&
&               cmplx(paw_ij(iatom)%noccmmp(1,im,im1,idijeff),zero,kind=dp)
             end do
           end do
         end do
       end if
     end do
     if(paw_dmft%dmftcheck==2) option=1
     if(paw_dmft%dmftcheck<=1) option=0
     call diff_matlu("LDA+U density matrix from INPUT wfk",&
&     "Direct calculation of density matrix with psichi from DIAGONALIZED wfk",&
&     matlu_temp,lda_occup%matlu,natom,option,tol3,ierrr) !tol1 tol2 tol3
     if(ierrr==-1) then
       write(message,'(10a)') ch10,&
&       '    -> These two quantities should agree if three conditions are fulfilled',ch10,&
&       '         -  input wavefunctions come from the same Hamiltonien (e.g LDA/GGA)',ch10,&
&       '         -  dmatpuopt is equal to 1',ch10,&
&       '         -  all valence states are in the valence',ch10,&
&       '    (for experts users: it is not compulsary that these conditions are fulfilled)'
       call wrtout(std_out,message,'COLL')
     end if
!    write(message,'(2a)') ch10,&
!    &   '  ***** => Calculations of density matrices with projections and in LDA+U are coherent****'
!    call wrtout(std_out,message,'COLL')

     call destroy_matlu(matlu_temp,natom)
     ABI_DATATYPE_DEALLOCATE(matlu_temp)
   else
     write(message,'(2a)') ch10,&
&     '  Warning: Consistency of density matrices computed from projection has not been checked: use dmftcheck>=2 '
     call wrtout(std_out,message,'COLL')
   end if

   call destroy_matlu(loc_norm_check,natom)
   ABI_DATATYPE_DEALLOCATE(loc_norm_check)
   call destroy_matlu(loc_occ_check,natom)
   ABI_DATATYPE_DEALLOCATE(loc_occ_check)

   call destroy_matlu(xocc_check,natom)
   call destroy_matlu(xnorm_check,natom)
   ABI_DATATYPE_DEALLOCATE(xocc_check)
   ABI_DATATYPE_DEALLOCATE(xnorm_check)
 endif

 if(present(nbandkss)) then
   if((me.eq.0.and.nbandkss/=0).or.(paw_dmft%dmft_kspectralfunc==1)) then
!     opt_renorm=1 ! if ucrpa==1, no need for individual orthonormalization
     opt_renorm=3
     if(dtset%ucrpa>=1.or.paw_dmft%dmft_kspectralfunc==1) opt_renorm=2
     call psichi_renormalization(cryst_struc,paw_dmft,pawang,opt=opt_renorm)
     call psichi_print(dtset,cryst_struc%nattyp,cryst_struc%ntypat,nkpt,my_nspinor,&
&     nsppol,paw_dmft,pawtab,psps,t2g,x2my2d)
   end if ! proc=me
 end if
!!***

 CONTAINS

!!****f* m_datafordmft/psichi_print
!! NAME
!!  psichi_print
!!
!! FUNCTION
!!  Print psichi for reference
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  maxnproju = maximum number of projector for LDA+U species
!!  nattyp(ntypat)= # atoms of each type
!!  mband= number of bands
!!  nkpt=number of k points
!!  my_nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  paw_dmft <type(paw_dmft)>=paw data for the self-consistency
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psichi(2,nsppol,nkpt,mband,my_nspinor,dtset%natom,(2*paw_dmft%maxlpawu+1))) projections for DMFT
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!!
!! SIDE EFFECTS
!!  print psichi in forlb.ovlp
!!
!! PARENTS
!!      datafordmft
!!
!! CHILDREN
!!
!! SOURCE

subroutine psichi_print(dtset,nattyp,ntypat,nkpt,my_nspinor,&
&nsppol,paw_dmft,pawtab,psps,t2g,x2my2d)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,my_nspinor,nsppol,ntypat
!arrays
 integer, intent(in) :: nattyp(ntypat)
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 logical t2g,x2my2d
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_dmft_type), intent(in) :: paw_dmft
!Local variables ------------------------------------
 integer :: ibg,isppol,ikpt,iband,ibandc,ispinor,icat,itypat,lmn_size
 integer :: iat,iatom,jj1,ilmn,m1,nband_k,unt
 integer :: m1_t2g,ll,m1_x2my2d
 real(dp) :: chinorm
 character(len=500) :: msg

! *********************************************************************
   ll=1

   if (open_file('forlb.ovlp',msg,newunit=unt,form='formatted',status='unknown') /= 0) then
     MSG_ERROR(msg)
   end if
   rewind(unt)

!  Header for calc_uCRPA.F90
   if  (COUNT(pawtab(:)%lpawu.NE.-1).EQ.1) then
     do  itypat=1,ntypat
       if(t2g) then
         if(pawtab(itypat)%lpawu.ne.-1) write(unt,*) "l= ",ll,itypat
       else if(x2my2d) then
         if(pawtab(itypat)%lpawu.ne.-1) write(unt,*) "l= ",ll-1,itypat
       else
         if(pawtab(itypat)%lpawu.ne.-1) write(unt,*) "l= ",pawtab(itypat)%lpawu,itypat
       end if
     end do
   else
     write(unt,*) "More than one correlated species"
   end if

   write(unt,*) "Bands ",paw_dmft%dmftbandi,paw_dmft%dmftbandf

   ibg=0
   do isppol=1,nsppol
     do ikpt=1,nkpt
!      rewind(1023)
       write(unt,'(a6,2x,i6)') "ikpt =",ikpt
       nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
       ibandc=0
       do iband=1,nband_k
         if(paw_dmft%band_in(iband)) then
           ibandc=ibandc+1
           write(unt,'(a8,2x,i6)') " iband =",iband
         end if
         do ispinor=1,my_nspinor
           icat=1
!          write(std_out,*) isppol,ikpt,iband,ispinor
           iat=0 ! to declare
           do itypat=1,dtset%ntypat
             lmn_size=pawtab(itypat)%lmn_size
!            write(std_out,*) isppol,ikpt,iband,ispinor
             do iatom=icat,icat+nattyp(itypat)-1
               iat=iat+1
               jj1=0
               if(pawtab(itypat)%lpawu.ne.-1) then
!                chinorm=(pawtab(itypat)%phiphjint(1))
                 chinorm=1.d0
!                write(std_out,*) isppol,ikpt,iband,ispinor,iat
                 m1_t2g=0
                 m1_x2my2d=0
                 do ilmn=1,lmn_size
!                  write(std_out,*) ilmn
!                  ------------ Select l=lpawu.  ---------------------------------------
                   if (psps%indlmn(1,ilmn,itypat)==pawtab(itypat)%lpawu.and.psps%indlmn(3,ilmn,itypat)==1) then
!                    ------------ Check that the band is choosen (within nbandi and nbandf)
                     if(paw_dmft%band_in(iband)) then
                       jj1=jj1+1
                       if(jj1>(2*pawtab(itypat)%lpawu+1)) then
                         write(message,'(a,a,i4,i5,i4)') ch10," Error 2 in datafordmft",jj1,pawtab(itypat)%lpawu
                         call wrtout(std_out,  message,'COLL')
                         MSG_ERROR("Aborting now")
                       end if ! jj1
!                      if(jj1>pawtab(dtset%typat(iatom))%nproju*(2*lpawu+1)) then
!                      write(message,'(a,a,a,a)')  ch10,&
!                      &                         'BUG: jj1 is not correct in datafordmft psichi_print',ch10,&
!                      &                         'Action: CONTACT Abinit group'
!                      call wrtout(std_out,  message,'COLL')
!                      call abi_abort('COLL')
!                      end if ! jj1
                       m1=psps%indlmn(2,ilmn,itypat)+psps%indlmn(1,ilmn,itypat)+1
!                      ----- Print only when the sum over projectors is done
!                      write(std_out,*) ilmn,m1
                       if(t2g) then
                         if(m1==1.or.m1==2.or.m1==4) then
                           m1_t2g=m1_t2g+1
                           write(unt,'(3i6,3x,2f23.15)') isppol, iat, m1,&
&                           real(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1_t2g))/chinorm,&
&                           aimag(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1_t2g))/chinorm
                         end if
                       else if(x2my2d) then
                         if(m1==5) then
                           m1_x2my2d=1
                           write(unt,'(3i6,3x,2f23.15)') isppol, iat, m1,&
&                           real(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1_x2my2d))/chinorm,&
&                           aimag(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1_x2my2d))/chinorm
                         end if
                       else
                         write(unt,'(3i6,3x,2f23.15)') isppol, iat, m1,&
&                         real(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1))/chinorm,&
&                         aimag(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1))/chinorm
                       end if
!                      if natom=1 then jj1 maximal value should be 2*lpawu+1
                     end if ! paw_dmft%band_in
                   end if
                 end do !ilmn
               end if ! lpawu.ne.-1
             end do ! iatom
             icat=icat+nattyp(itypat)
           end do ! itypat
         end do ! ispinor
       end do !iband
       ibg=ibg+nband_k*my_nspinor
     end do !ikpt
   end do ! isppol
!   write(unt,*) "Fermi level (in Ryd)="
!   write(unt,*) fermie*two
   close(unt)
 end subroutine psichi_print
!!***

!!****f* m_datafordmft/psichi_check
!! NAME
!!  psichi_check
!!
!! FUNCTION
!!  Check psichi: compute occupations
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  maxnproju = maximum number of projector for LDA+U species
!!  nattyp(ntypat)= # atoms of each type
!!  mband= number of bands
!!  nkpt=number of k points
!!  my_nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat= number of species
!!  paw_dmft <type(paw_dmft)>=paw data for the self-consistency
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!!  OUTPUTS:
!!  xocc_check(nsppol,my_nspinor,my_nspinor,natom,2*maxlpawu+1,2*maxlpawu+1): density matrix
!!  xnorm_check(nsppol,my_nspinor,my_nspinor,natom,2*maxlpawu+1,2*maxlpawu+1): matrix of norms
!!
!! SIDE EFFECTS
!!  check psichi: compute norm and occupations
!!
!! PARENTS
!!      datafordmft
!!
!! CHILDREN
!!
!! SOURCE

subroutine psichi_check(dtset,nattyp,nkpt,my_nspinor,&
& nsppol,ntypat,paw_dmft,pawtab,psps,xocc_check,xnorm_check)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,my_nspinor,nsppol,ntypat
!arrays
 integer, intent(in) :: nattyp(ntypat)
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(matlu_type), intent(inout) :: xocc_check(dtset%natom) !vz_i
 type(matlu_type), intent(inout) :: xnorm_check(dtset%natom) !vz_i
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
!Local variables ------------------------------------
 integer :: band_index,facnsppol,ibg,isppol,ikpt,iband,ibandc,ispinor,icat,itypat
 integer :: iat,iatom,ilmn,lmn_size,m,m1,nband_k
 complex(dpc) :: psichic,psichic1
 real(dp) :: chinorm

! *********************************************************************
   facnsppol=1
   if(my_nspinor==1.and.nsppol==1) then
     facnsppol=2
   end if

   ibg=0
   band_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
       ibandc=0
       do iband=1,nband_k
         if(paw_dmft%band_in(iband)) ibandc=ibandc+1
         do ispinor=1,my_nspinor
           icat=1
           iat=0
           do itypat=1,dtset%ntypat
             lmn_size=pawtab(itypat)%lmn_size
             do iatom=icat,icat+nattyp(itypat)-1
               iat=iat+1
!              ------------ Select correlated atoms
               if(paw_dmft%lpawu(iatom).ne.-1) then
                 chinorm=1.d0
                 do ilmn=1,lmn_size
!                  ------------ Select l=lpawu.
                   if (psps%indlmn(1,ilmn,itypat)==paw_dmft%lpawu(iatom).and.&
&                   psps%indlmn(3,ilmn,itypat)==1) then
                     do ilmn1=1,lmn_size
!                      ------------ Select l=lpawu and do not sum over projectors
!                      (this is already done in paw_dmft%psichi)
                       if (psps%indlmn(1,ilmn1,itypat)==paw_dmft%lpawu(iatom).and.&
&                       psps%indlmn(3,ilmn1,itypat)==1) then
!                        ------------ Check that the band is choosen (within nbandi and nbandf)
                         if(paw_dmft%band_in(iband)) then
                           m=psps%indlmn(2,ilmn,itypat)+psps%indlmn(1,ilmn,itypat)+1
                           m1=psps%indlmn(2,ilmn1,itypat)+psps%indlmn(1,ilmn,itypat)+1
                           if(psps%indlmn(3,ilmn,itypat)==1) then
                             do ispinor1=1,my_nspinor
                               psichic=paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m)
                               psichic1=paw_dmft%psichi(isppol,ikpt,ibandc,ispinor1,iat,m1)
!                              ------------ Compute occupation matrix
                               xocc_check(iatom)%mat(m,m1,isppol,ispinor,ispinor1)=&
&                               xocc_check(iatom)%mat(m,m1,isppol,ispinor,ispinor1)&
!                              &               +conjg(psichic)*psichic1*dtset%wtk(ikpt)*facpara*occ(iband+band_index)
&                               +conjg(psichic1)*psichic*dtset%wtk(ikpt)*facpara*occ(iband+band_index)/facnsppol
!                              ------------ Compute norm (should be equal to noccmmp
!                              (dmatpuopt=1) if all bands are taken into account)
                               xnorm_check(iatom)%mat(m,m1,isppol,ispinor,ispinor1)=&
&                               xnorm_check(iatom)%mat(m,m1,isppol,ispinor,ispinor1)&
!                              &               +conjg(psichic)*psichic1*dtset%wtk(ikpt)*facpara
&                               +conjg(psichic1)*psichic*dtset%wtk(ikpt)*facpara
                             end do ! ispinor1
                           end if
                         end if ! paw_dmft%band_in
                       end if
                     end do !ilmn1
                   end if
                 end do !ilmn
               end if ! lpawu.ne.-1
             end do ! iatom
             icat=icat+nattyp(itypat)
           end do ! itypat
         end do ! ispinor
       end do !iband
       band_index=band_index+nband_k
       ibg=ibg+nband_k*my_nspinor
     end do !ikpt
   end do ! isppol

 end subroutine psichi_check
!DBG_EXIT("COLL")
end subroutine datafordmft
!!***

!!****f* m_datafordmft/compute_levels
!! NAME
!! compute_levels
!!
!! FUNCTION
!! Compute levels for ctqmc
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      m_hubbard_one,qmc_prep_ctqmc
!!
!! CHILDREN
!!      checkdiag_matlu,loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

 subroutine compute_levels(cryst_struc,energy_level,hdc,pawang,paw_dmft,nondiag)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(pawang_type), intent(in) :: pawang
 type(oper_type), intent(in) :: hdc
 type(paw_dmft_type), intent(in)  :: paw_dmft
 type(oper_type),intent(inout)  :: energy_level !vz_i
 logical, optional, intent(out) :: nondiag

!Local variables ------------------------------
! scalars
 integer :: iatom,iband,ikpt,im1,isppol,ispinor
 integer :: lpawu,mbandc,natom,nspinor,nsppol,nkpt
 character(len=500) :: message
! arrays
!************************************************************************

 mbandc=paw_dmft%mbandc
 nkpt=paw_dmft%nkpt
 nsppol=paw_dmft%nsppol
 nspinor=paw_dmft%nspinor
 natom=paw_dmft%natom
 if(present(nondiag)) nondiag=.false.

!========================
!Get KS eigenvalues
!========================
 do iband=1,mbandc
   do ikpt=1,nkpt
     do isppol=1,nsppol
!      Take \epsilon_{nks}
!      ========================
       energy_level%ks(isppol,ikpt,iband,iband)=paw_dmft%eigen_lda(isppol,ikpt,iband)
     end do
   end do
 end do


!======================================================================
!Compute atomic levels from projection of \epsilon_{nks} and symetrize
!======================================================================
 call loc_oper(energy_level,paw_dmft,1)
! write(message,'(a,2x,a,f13.5)') ch10," == Print Energy levels before sym and only LDA"
! call wrtout(std_out,message,'COLL')
! call print_matlu(energy_level%matlu,natom,1)
 do iatom = 1 , natom
   lpawu=paw_dmft%lpawu(iatom)
   if(lpawu/=-1) then
     do isppol=1,nsppol
       do ispinor=1,nspinor
         do im1=1,2*lpawu+1
           energy_level%matlu(iatom)%mat(im1,im1,isppol,ispinor,ispinor)= &
&           energy_level%matlu(iatom)%mat(im1,im1,isppol,ispinor,ispinor)&
&           -hdc%matlu(iatom)%mat(im1,im1,isppol,ispinor,ispinor)-paw_dmft%fermie
         end do
       end do
     end do
!    write(std_out,*) "DC,fermie",hdc%matlu(iatom)%mat(1,1,1,1,1),paw_dmft%fermie
   end if
 end do ! natom
 call sym_matlu(cryst_struc,energy_level%matlu,pawang,paw_dmft)
 if(present(nondiag)) call checkdiag_matlu(energy_level%matlu,natom,tol7,nondiag)

 write(message,'(a,2x,a,f13.5)') ch10," == Print Energy levels for Fermi Level=",paw_dmft%fermie
 call wrtout(std_out,message,'COLL')
!call print_oper(energy_level,1,paw_dmft,1)
 call print_matlu(energy_level%matlu,natom,1)



 end subroutine compute_levels
!!***


!!****f* m_datafordmft/psichi_renormalization
!! NAME
!! psichi_renormalization
!!
!! FUNCTION
!! Renormalize psichi.
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>= crystal structure data.
!!  paw_dmft =  data for LDA+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!
!! OUTPUT
!!  paw_dmft%psichi(nsppol,nkpt,mband,nspinor,dtset%natom,(2*maxlpawu+1))): projections <Psi|chi> are orthonormalized.
!!
!! NOTES
!!
!! PARENTS
!!      m_datafordmft,m_dmft
!!
!! CHILDREN
!!      add_matlu,destroy_oper,gather_matlu,identity_oper,init_oper
!!      invsqrt_matrix,loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine psichi_renormalization(cryst_struc,paw_dmft,pawang,opt)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 type(pawang_type), intent(in) :: pawang
 integer, optional, intent(in) :: opt

!Local variables ------------------------------
!scalars
 integer :: iatom,ib,ikpt,im,ispinor,isppol,jkpt
 integer :: natom,mbandc,ndim,nkpt,nspinor,nsppol,option
 real(dp), pointer ::  temp_wtk(:) => null()
 real(dp) ::  pawprtvol
 type(oper_type) :: norm
 type(oper_type) :: oper_temp
 character(len=500) :: message
!arrays
! real(dp),allocatable :: e0pde(:,:,:),omegame0i(:)
!************************************************************************

 DBG_ENTER("COLL")

 option=3
 if(present(opt)) then
   if(opt==2) option=2
   if(opt==3) option=3
 end if
 pawprtvol=2

 nsppol  = paw_dmft%nsppol
 nkpt    = paw_dmft%nkpt
 mbandc  = paw_dmft%mbandc
 natom   = cryst_struc%natom
 nspinor = paw_dmft%nspinor


!== Normalize psichi
 if(option==1) then
!  ====================================
!  == simply renormalize psichi =======
!  ====================================
   write(message,'(2a)') ch10," Psichi are renormalized  "
   call wrtout(std_out,  message,'COLL')
   do isppol=1,nsppol
     do ikpt=1,nkpt
       do ib=1,mbandc
         do iatom=1,natom
           if(paw_dmft%lpawu(iatom).ne.-1) then
             ndim=2*paw_dmft%lpawu(iatom)+1
             do im=1,ndim
               do ispinor=1,nspinor
!                write(std_out,*) "psichi1",paw_dmft%psichi(isppol,ikpt,ib,ispinor,iatom,im)
                 paw_dmft%psichi(isppol,ikpt,ib,ispinor,iatom,im)=     &
&                 paw_dmft%psichi(isppol,ikpt,ib,ispinor,iatom,im)/    &
&                 sqrt(real(norm%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)))
               end do ! ispinor
             end do ! im
           end if
         end do ! iatom
       end do ! ib
     end do ! ikpt
   end do ! isppol
!  todo_ab introduce correct orthonormalization in the general case.

 else if(option==2) then ! option==2
!  ====================================
!  == renormalize k-point after k-point
!  ====================================
   ABI_ALLOCATE(temp_wtk,(1))

   write(message,'(2a,i5)') ch10,' Nb of K-point',nkpt
   call wrtout(std_out,message,'COLL')
   do jkpt=1,nkpt  ! jkpt
     write(message,'(2a,i5)') ch10,'  == Renormalization for K-point',jkpt
     call wrtout(std_out,message,'COLL')
     temp_wtk(1)=one

     call normalizepsichi(cryst_struc,1,paw_dmft,pawang,temp_wtk,jkpt)
   end do ! jkpt
   ABI_DEALLOCATE(temp_wtk)

 else if(option==3)  then  ! option==3
!  ====================================
!  == renormalize the sum over k-points
!  ====================================
!  allocate(temp_wtk(nkpt))
   temp_wtk=>paw_dmft%wtk
   write(message,'(6a)') ch10,'  ====================================== '&
&   ,ch10,'  == Renormalization for all K-points == '&
&   ,ch10,'  ======================================='
   call wrtout(std_out,message,'COLL')
   call normalizepsichi(cryst_struc,nkpt,paw_dmft,pawang,temp_wtk)
!  deallocate(temp_wtk)

 end if

!== Change back repr for norm


!===============================================
!==  Compute norm with new psichi
!===============================================
 call init_oper(paw_dmft,norm,nkpt=paw_dmft%nkpt,wtk=paw_dmft%wtk)
!== Build identity for norm%ks (option=1)
 call identity_oper(norm,1)

!== Deduce norm%matlu from norm%ks with new psichi
 if (paw_dmft%dmft_kspectralfunc==1) then
   call init_oper(paw_dmft,oper_temp)
   do jkpt=1,nkpt  ! jkpt
     call loc_oper(norm,paw_dmft,1,jkpt=jkpt)
     write(message,'(2a,i5)') &
&     ch10,"  == Check: Overlap with renormalized psichi for k-point",jkpt
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm%matlu,natom,prtopt=1)
!== Check that norm is now the identity
     call identity_oper(oper_temp,2)
     call diff_matlu('Overlap after renormalization','Identity',&
&     norm%matlu,oper_temp%matlu,norm%natom,1,tol6,zero_or_one=1)
   enddo
   call destroy_oper(oper_temp)
 else !dmft_kspectralfunc
   call loc_oper(norm,paw_dmft,1)


!== Print norm%matlu unsymetrized with new psichi
   if(pawprtvol>2) then
     write(message,'(4a,2a)') &
&     ch10,"  == Check: Overlap with renormalized psichi without symetrization is == "
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm%matlu,natom,prtopt=1)
   end if


!== Symetrise norm%matlu with new psichi
   call sym_matlu(cryst_struc,norm%matlu,pawang,paw_dmft)

!== Print norm%matlu symetrized with new psichi
   if(pawprtvol>2) then
     write(message,'(4a,2a)') &
&     ch10,"  == Check: Overlap with renormalized psichi and symetrization is =="
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm%matlu,natom,prtopt=1,opt_diag=-1)
   end if

!== Check that norm is now the identity
   call init_oper(paw_dmft,oper_temp)
   call identity_oper(oper_temp,2)
   call diff_oper('Overlap after renormalization','Identity',&
&   norm,oper_temp,1,tol6)
   call destroy_oper(oper_temp)

   call destroy_oper(norm)
 endif ! dmft_kspectralfunc

 paw_dmft%lpsichiortho=1

 DBG_EXIT("COLL")

 CONTAINS
!===========================================================
!!***

!!****f* psichi_renormalization/normalizepsichi
!! NAME
!!  normalizepsichi
!!
!! FUNCTION
!!  normalizepsichi psichi from norm
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  change psichi: normalizepsichi it
!!
!! PARENTS
!!      psichi_renormalization
!!
!! CHILDREN
!!      add_matlu,destroy_oper,gather_matlu,identity_oper,init_oper
!!      invsqrt_matrix,loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine normalizepsichi(cryst_struc,nkpt,paw_dmft,pawang,temp_wtk,jkpt)

!Arguments ------------------------------------
 integer,intent(in) :: nkpt
 integer,optional,intent(in) :: jkpt
 real(dp),pointer :: temp_wtk(:)
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 type(pawang_type), intent(in) :: pawang

!Local variables ------------------------------
!scalars
 integer :: diag,iatom,ib,ikpt1,im,im1,ispinor,ispinor1,isppol,isppol1,jc,jc1
 integer :: tndim
 integer :: natom,mbandc,ndim,nspinor,nsppol,iortho,natomcor
 integer :: itot,itot1,dimoverlap,iatomcor
 real(dp) :: pawprtvol
 type(oper_type) :: norm1,norm2,norm3
 character(len=500) :: message
 complex(dpc),allocatable :: wan(:,:,:),sqrtmatinv(:,:),wanall(:)
 type(coeff2c_type), allocatable :: overlap(:)
 complex(dpc), allocatable :: largeoverlap(:,:)
 complex(dpc), allocatable :: psichivect(:,:)
!arrays
! real(dp),allocatable :: e0pde(:,:,:),omegame0i(:)
!************************************************************************
   nsppol  = paw_dmft%nsppol
   mbandc  = paw_dmft%mbandc
   natom   = cryst_struc%natom
   nspinor = paw_dmft%nspinor
   pawprtvol=3
   diag=0
   natomcor=0
   dimoverlap=0
   do iatom=1,natom
     if(paw_dmft%lpawu(iatom).ne.-1) then
       natomcor=natomcor+1
       ndim=2*paw_dmft%lpawu(iatom)+1
       tndim=nspinor*ndim
       dimoverlap=dimoverlap+tndim
      ! write(6,*) "atom, dimoverlap",iatom,dimoverlap,natomcor
     end if
   end do

   if(nkpt/=1.and.present(jkpt)) then
     message = 'BUG in psichi_normalization'
     MSG_ERROR(message)
   end if

   iortho=1
  ! write(6,*) "nkpt, iortho",nkpt,iortho
   !if (natomcor>1) iortho=2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  First case: usual case (should in fact be used only for one atom and nkpt=1)
   if((.not.present(jkpt)).and.iortho==1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    *********************************************************************
     call init_oper(paw_dmft,norm1,nkpt=nkpt,wtk=temp_wtk)

!    == Build identity for norm1%ks (option=1)
     call identity_oper(norm1,1)

     if(nkpt==1.and.present(jkpt)) then
       call loc_oper(norm1,paw_dmft,1,jkpt=jkpt)
     end if
     if(.not.present(jkpt)) then
       call loc_oper(norm1,paw_dmft,1)
     end if
     if(nkpt>1) then
       call sym_matlu(cryst_struc,norm1%matlu,pawang,paw_dmft)
     end if

     if(pawprtvol>2) then
       write(message,'(2a)') ch10,'  - Print norm with current psichi '
       call wrtout(std_out,message,'COLL')
       call print_matlu(norm1%matlu,natom,prtopt=1,opt_exp=1)
     end if
!    ==-------------------------------------
!    == Start loop on atoms
     ABI_DATATYPE_ALLOCATE(overlap,(natom))
     do iatom=1,natom
       if(paw_dmft%lpawu(iatom).ne.-1) then
         ndim=2*paw_dmft%lpawu(iatom)+1
         tndim=nsppol*nspinor*ndim
         ABI_ALLOCATE(overlap(iatom)%value,(tndim,tndim))
         overlap(iatom)%value=czero
       end if
     end do
!    ==-------------------------------------

!    built large overlap matrix
     write(message,'(2a)') ch10,'  - Overlap (before orthonormalization) -'
     call wrtout(std_out,message,'COLL')
     call gather_matlu(norm1%matlu,overlap,cryst_struc%natom,option=1,prtopt=1)
     call destroy_oper(norm1)



     do iatom=1,natom
       if(paw_dmft%lpawu(iatom).ne.-1) then
         ndim=2*paw_dmft%lpawu(iatom)+1
         tndim=nsppol*nspinor*ndim
         ABI_ALLOCATE(sqrtmatinv,(tndim,tndim))

!        == Compute Inverse Square root of overlap : O^{-0.5}
         !do im=1,tndim
         !  do im1=1,tndim
         !    !write(message,'(a,1x,a,e21.14,a,e21.14,a)') "overlap", &
         !    !"(",real(overlap(1)%value(im,im1)),",",aimag(overlap(1)%value(im,im1)),")"
         !    write(6,*) "overlap",overlap(iatom)%value(im,im1)
         !  enddo
         !enddo
         !stop
         !call wrtout(std_out,message,'COLL')
         if(diag==0) then
           call invsqrt_matrix(overlap(iatom)%value,tndim)
           sqrtmatinv=overlap(iatom)%value
         else
           sqrtmatinv(:,:)=czero
           do ib=1,tndim
             sqrtmatinv(ib,ib)=cone/(sqrt(overlap(iatom)%value(ib,ib)))
           end do
         end if

!        == Apply O^{-0.5} on psichi
         ABI_ALLOCATE(wan,(nsppol,nspinor,ndim))
!        write(std_out,*) mbandc,nsppol,nspinor,ndim
!        write(std_out,*)  paw_dmft%psichi(1,1,1,1,1,1)
         do ikpt=1,nkpt
           do ib=1,mbandc
             if(present(jkpt)) then
               ikpt1=jkpt
             else
               ikpt1=ikpt
             end if
             jc=0
             wan=czero
             do isppol=1,nsppol
               do ispinor=1,nspinor
                 do im=1,ndim
!                  write(std_out,*) "psichi", paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)
                   jc=jc+1
                   jc1=0
                   do isppol1=1,nsppol
                     do ispinor1=1,nspinor
                       do im1=1,ndim
                         jc1=jc1+1
                         wan(isppol,ispinor,im)= wan(isppol,ispinor,im) &
&                         + paw_dmft%psichi(isppol1,ikpt1,ib,ispinor1,iatom,im1)*sqrtmatinv(jc,jc1)
                       end do ! ispinor1
                     end do ! isppol1
                   end do ! im1
                 end do ! im
               end do ! ispinor
             end do !  isppol
             do isppol=1,nsppol
               do ispinor=1,nspinor
                 do im=1,ndim
                   paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)=wan(isppol,ispinor,im)
!                  write(std_out,*) "psichi2", paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)
                 end do ! ispinor
               end do ! isppol
             end do ! im
           end do ! ib
         end do ! ikpt
         ABI_DEALLOCATE(wan)
         ABI_DEALLOCATE(sqrtmatinv)
!        write(std_out,*)  paw_dmft%psichi(1,1,1,1,1,1)

!        ==-------------------------------------
       end if ! lpawu.ne.-1
     end do ! iatom
!    == End loop on atoms
!    ==-------------------------------------
     do iatom=1,natom
       if(paw_dmft%lpawu(iatom).ne.-1) then
         ABI_DEALLOCATE(overlap(iatom)%value)
       end if
     end do
     ABI_DATATYPE_DEALLOCATE(overlap)

!    ======================================================================
!    == Check norm with new psichi.
!    ======================================================================

     call init_oper(paw_dmft,norm1,nkpt=nkpt,wtk=temp_wtk)

     call identity_oper(norm1,1)

     if(nkpt==1.and.present(jkpt)) then
       call loc_oper(norm1,paw_dmft,1,jkpt=jkpt)
     end if
     if(.not.present(jkpt)) then
       call loc_oper(norm1,paw_dmft,1)
     end if

     if (nkpt>1) then
       call sym_matlu(cryst_struc,norm1%matlu,pawang,paw_dmft)
     end if

     if(pawprtvol>2) then
       write(message,'(2a)') ch10,'  - Print norm with new psichi '
       call wrtout(std_out,message,'COLL')
       call print_matlu(norm1%matlu,natom,prtopt=1)
     end if

!    ======================================================================
!    == Check that norm-identity is zero
!    ======================================================================
     call init_oper(paw_dmft,norm2,nkpt=nkpt,wtk=temp_wtk)
     call init_oper(paw_dmft,norm3,nkpt=nkpt,wtk=temp_wtk)
     call identity_oper(norm2,2)
     call add_matlu(norm1%matlu,norm2%matlu,norm3%matlu,natom,-1)
     call destroy_oper(norm2)
     if(pawprtvol>2) then
       write(message,'(2a)') ch10,'  - Print norm with new psichi minus Identity '
       call wrtout(std_out,message,'COLL')
       call print_matlu(norm3%matlu,natom,prtopt=1,opt_exp=1)
     end if
     call destroy_oper(norm3)

     call destroy_oper(norm1)
!    call flush(std_out)           ! debug debug  debug   debug
!    MSG_ERROR("Stop for debugging")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  New implementation, several atoms, general case.
   else if(present(jkpt).or.iortho==2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ABI_ALLOCATE(largeoverlap,(dimoverlap,dimoverlap))
     ABI_ALLOCATE(psichivect,(mbandc,dimoverlap))
     ABI_ALLOCATE(sqrtmatinv,(dimoverlap,dimoverlap))
     ABI_ALLOCATE(wanall,(dimoverlap))


!    Big loop over isppol
     do isppol=1,nsppol

       do ib=1,mbandc
         itot=0
         do iatom=1,natom
           if(paw_dmft%lpawu(iatom).ne.-1) then
             ndim=2*paw_dmft%lpawu(iatom)+1
             do im=1,ndim
               do ispinor=1,nspinor
                 itot=itot+1
                 if(itot>dimoverlap) write(std_out,*) "itot>ndim",itot,ndim
                ! write(6,*) "ib,iatom,im,ispinor",ib,iatom,im,ispinor,jkpt
                 psichivect(ib,itot)= paw_dmft%psichi(isppol,jkpt,ib,ispinor,iatom,im)
               enddo ! ispinor
             enddo ! im
           endif
         enddo ! iatom
       enddo ! ib


!     calculation of overlap
       largeoverlap=czero
       do ib=1,mbandc
         do itot=1,dimoverlap
           do itot1=1,dimoverlap
              largeoverlap(itot,itot1)=largeoverlap(itot,itot1)+ &
&              psichivect(ib,itot)*conjg(psichivect(ib,itot1))
           enddo ! itot1
         enddo ! itot
       enddo ! ib

!     Math: orthogonalisation of overlap
       write(std_out,*)"jkpt=",jkpt
       do itot=1,dimoverlap
         write(std_out,'(100f7.3)') (largeoverlap(itot,itot1),itot1=1,dimoverlap)
       enddo
       call invsqrt_matrix(largeoverlap,dimoverlap)
       sqrtmatinv=largeoverlap
       write(std_out,*)"jkpt=",jkpt
       do itot=1,dimoverlap
         write(std_out,'(100f7.3)') (sqrtmatinv(itot,itot1),itot1=1,dimoverlap)
       enddo
       write(std_out,*)"jkpt=",jkpt
       do itot=1,dimoverlap
         write(std_out,'(100e9.3)') (sqrtmatinv(itot,itot1),itot1=1,dimoverlap)
       enddo


       do ib=1,mbandc
         wanall=czero
         do itot=1,dimoverlap
           do itot1=1,dimoverlap
              wanall(itot)= wanall(itot)+psichivect(ib,itot1)*sqrtmatinv(itot,itot1)
           enddo ! itot1
     !      write(std_out,'(3i3,2x,i3,2x,2e15.5,2x,2e15.5)') jkpt,isppol,ib,itot,psichivect(ib,itot),wanall(itot)
         enddo ! itot
         iatomcor=0
         do itot=1,dimoverlap
           psichivect(ib,itot)=wanall(itot)
         enddo
        ! do iatom=1,natom
        !   if(paw_dmft%lpawu(iatom).ne.-1) then
        !     ndim=2*paw_dmft%lpawu(iatom)+1
        !     iatomcor=iatomcor+1
        !     do im=1,ndim
        !       do ispinor=1,nspinor
        !         paw_dmft%psichi(isppol,jkpt,ib,ispinor,iatom,im)=wanall(iatomcor,isppol,ispinor,im)
        !       end do ! ispinor
        !     end do ! im
        !   endif
        ! enddo ! iatom
       enddo ! ib


!     calculation of overlap (check)
       largeoverlap=czero
       do ib=1,mbandc
         do itot=1,dimoverlap
           do itot1=1,dimoverlap
              largeoverlap(itot,itot1)=largeoverlap(itot,itot1)+ &
&              psichivect(ib,itot)*conjg(psichivect(ib,itot1))
           enddo ! itot1
         enddo ! itot
       enddo ! ib

       write(std_out,*)"jkpt=",jkpt
       do itot=1,dimoverlap
         write(std_out,'(100f7.3)') (largeoverlap(itot,itot1),itot1=1,dimoverlap)
       enddo


!      psichivect -> psichi
       do ib=1,mbandc
         itot=0
         do iatom=1,natom
           if(paw_dmft%lpawu(iatom).ne.-1) then
             ndim=2*paw_dmft%lpawu(iatom)+1
             iatomcor=iatomcor+1
             do im=1,ndim
               do ispinor=1,nspinor
                 itot=itot+1
                 paw_dmft%psichi(isppol,jkpt,ib,ispinor,iatom,im)=psichivect(ib,itot)
               end do ! ispinor
             end do ! im
           endif
         enddo ! iatom
       enddo ! ib


!   End big loop over isppol
     enddo !ispppol

     ABI_DEALLOCATE(psichivect)
     ABI_DEALLOCATE(sqrtmatinv)
     ABI_DEALLOCATE(wanall)
     ABI_DEALLOCATE(largeoverlap)

   endif

 end subroutine normalizepsichi

end subroutine psichi_renormalization
!!***


!!****f* m_datafordmft/hybridization_asymptotic_coefficient
!! NAME
!! hybridization_asymptotic_coefficient
!!
!! FUNCTION
!! Compute some components for the limit of hybridization
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  lda_occup
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawtab <type(pawtab)>
!!
!! OUTPUT
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!
!! NOTES
!!
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!      add_matlu,copy_oper,destroy_oper,init_oper,loc_oper,prod_oper,sym_matlu
!!
!! SOURCE

subroutine hybridization_asymptotic_coefficient(cryst_struc,paw_dmft,pawang,hybri_coeff)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(pawang_type), intent(in) :: pawang
 type(matlu_type), intent(inout) :: hybri_coeff(paw_dmft%natom)
!Local variables ------------------------------
 type(oper_type)  :: ham_a
 type(oper_type)  :: ham_b
 type(oper_type)  :: ham_squarelocal
 type(oper_type)  :: ham_squareks
 integer :: iband1,iband2,ikpt,isppol
!************************************************************************

! call init_oper(paw_dmft,self_minus_hdc)
 call init_oper(paw_dmft,ham_a)
 call init_oper(paw_dmft,ham_b)
 call init_oper(paw_dmft,ham_squareks)
 call init_oper(paw_dmft,ham_squarelocal)

! Create self_minus_hdc%matlu = Sigma-Hdc in local orbitals
! call add_matlu(self%oper(paw_dmft%dmft_nwlo)%matlu,self%hdc%matlu,&
!&             self_minus_hdc%matlu,cryst_struc%natom,-1)

!   write(message,'(a,2x,a)') ch10,        "  == self_minus_hdc (1)"
!   call wrtout(std_out,message,'COLL')
!   call print_matlu(self_minus_hdc%matlu,paw_dmft%natom,1,opt_exp=1)

!! Create self_minus_hdc%ks = Upfold Sigma-Hdc
! call upfold_oper(self_minus_hdc,paw_dmft,1)
! call loc_oper(self_minus_hdc,paw_dmft,1)

!   write(message,'(a,2x,a)') ch10,        "  == self_minus_hdc (2)"
!   call wrtout(std_out,message,'COLL')
!   call print_matlu(self_minus_hdc%matlu,paw_dmft%natom,1,opt_exp=1)

! Create ham_a%ks = H_ks  in KS basis
!----------------------------------------------------
 do iband1=1,paw_dmft%mbandc
   do iband2=1,paw_dmft%mbandc
     do ikpt=1,paw_dmft%nkpt
       do isppol=1,paw_dmft%nsppol
         if(iband1==iband2) then
           ham_a%ks(isppol,ikpt,iband1,iband2) = cmplx(paw_dmft%eigen_lda(isppol,ikpt,iband1),0.d0,kind=dp)
         else
           ham_a%ks(isppol,ikpt,iband1,iband2) = czero
         end if
       end do
     end do
   end do
 end do

! Create ham_a%matlu = H_ks in local orbitals
!---------------------------------------------
 call loc_oper(ham_a,paw_dmft,1)

! Symetrise the local quantity (energy levels)
!---------------------------------------------
 call sym_matlu(cryst_struc,ham_a%matlu,pawang,paw_dmft)

! Create ham_b%ks : Duplicate both ks and local part of ham_a into ham_b
!-----------------------------------------------------------------------
 call copy_oper(ham_a,ham_b)

! Compute ham_squareks%ks   : Make a product of the two KS version
!------------------------------------------------------------------
 call prod_oper(ham_a,ham_b,ham_squareks,1)

! Compute ham_squareks%matlu
!---------------------------
 call loc_oper(ham_squareks,paw_dmft,1)

! Symetrise ham_squareks%matlu
!------------------------------
 call sym_matlu(cryst_struc,ham_squareks%matlu,pawang,paw_dmft)

!   write(message,'(a,2x,a)') ch10,        "  == squareks"
!   call wrtout(std_out,message,'COLL')
!   call print_matlu(ham_squareks%matlu,paw_dmft%natom,1,opt_exp=1)


! Compute ham_squarelocal%matlu
!-------------------------------
 call prod_oper(ham_a,ham_b,ham_squarelocal,2)

! Compute the product in local orbitals
!--------------------------------------
 call sym_matlu(cryst_struc,ham_squarelocal%matlu,pawang,paw_dmft)

!   write(message,'(a,2x,a)') ch10,        "  == squarelocal"
!   call wrtout(std_out,message,'COLL')
!   call print_matlu(ham_squarelocal%matlu,paw_dmft%natom,1,opt_exp=1)

! The difference of ham_squareks and ham_squarelocal
! gives the coefficient that we are looking for ( such that F_ij(iw_n) = -C_ij/(iw_n) ).
!----------------------------------------------------------------------------------------
 call add_matlu(ham_squareks%matlu,ham_squarelocal%matlu,hybri_coeff,cryst_struc%natom,-1)

 !  write(message,'(a,2x,a)') ch10,        "  == Coeff C_ij before sym"
 !  call wrtout(std_out,message,'COLL')
 !  call print_matlu(hybri_coeff,paw_dmft%natom,1,opt_exp=1)

! Symetrise the local quantity
!------------------------------
 call sym_matlu(cryst_struc,hybri_coeff,pawang,paw_dmft)

 !  write(message,'(a,2x,a)') ch10,        "  == Coeff C_ij after sym"
 !  call wrtout(std_out,message,'COLL')
 !  call print_matlu(hybri_coeff,paw_dmft%natom,1,opt_exp=1)

 call destroy_oper(ham_squarelocal)
! call destroy_oper(self_minus_hdc)
 call destroy_oper(ham_a)
 call destroy_oper(ham_b)
 call destroy_oper(ham_squareks)


end subroutine hybridization_asymptotic_coefficient
!!***

END MODULE m_datafordmft
!!***

