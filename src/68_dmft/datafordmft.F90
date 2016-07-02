!{\src2tex{textfont=tt}}
!!****f* ABINIT/datafordmft
!! NAME
!! datafordmft
!!
!! FUNCTION
!!  Compute psichi (and print some data for check)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2016 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
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
!!  mpi_enreg=informations about MPI parallelization
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

subroutine datafordmft(cryst_struc,cprj,dimcprj,dtset,eigen,fermie,&
& lda_occup,mband,mband_cprj,mkmem,mpi_enreg,nkpt,my_nspinor,nsppol,occ,&
& paw_dmft,paw_ij,pawang,pawtab,psps,usecprj,unpaw,nbandkss)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_io_tools,  only : open_file
 use m_matlu, only: matlu_type,init_matlu,sym_matlu,copy_matlu,print_matlu,diff_matlu,destroy_matlu
 use m_crystal, only : crystal_t
 use m_oper, only : oper_type

 use m_pawang, only : pawang_type
 use m_pawtab, only : pawtab_type
 use m_paw_ij, only : paw_ij_type
 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_free
 use m_paw_dmft, only: paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'datafordmft'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_68_dmft, except_this_one => datafordmft
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mband_cprj,mkmem
 integer,intent(in) :: nkpt,my_nspinor,nsppol
 integer,intent(in) :: unpaw,usecprj
 integer, optional, intent(in) :: nbandkss
 real(dp),intent(in) :: fermie
 type(MPI_type),intent(inout) :: mpi_enreg
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
 integer :: m1,m1_t2g,m1_t2g_mod,maxnproju,me,natom,nband_k,nband_k_cprj
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
 logical :: lprojchi,t2g
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
!call leave_new('COLL')

!----------------------------------- MPI-------------------------------------

 nbandi=paw_dmft%dmftbandi
 nbandf=paw_dmft%dmftbandf
 lprojchi=.false.
 lprojchi=.true.
 t2g=(paw_dmft%dmftqmc_t2g==1)
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
             lpawu=pawtab(dtset%typat(iatom))%lpawu
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
&                 ' Action: desactivate qmftqmc_t2g or use lpawu=1'
                 MSG_ERROR(message)
               end if
             end if
!            ----------   t2g case
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
                     if(jj1>pawtab(dtset%typat(iatom))%nproju*(2*lpawu+1)) then
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
!call leave_new('COLL')
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
&   nsppol,paw_dmft,pawtab,psps,t2g)
 end if ! proc=me
!==========================================================================
!********************* Check normalization  and occupations ***************
!==========================================================================
 ABI_DATATYPE_ALLOCATE(xocc_check,(natom))
 ABI_DATATYPE_ALLOCATE(xnorm_check,(natom))
 call init_matlu(natom,my_nspinor,nsppol,paw_dmft%lpawu,xocc_check)
 call init_matlu(natom,my_nspinor,nsppol,paw_dmft%lpawu,xnorm_check)
 call psichi_check(dtset,cryst_struc%nattyp,nkpt,my_nspinor,&
& nsppol,cryst_struc%ntypat,paw_dmft,pawtab,psps,xocc_check,xnorm_check)
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
!  print occupations
   write(message,'(2a,i4)')  ch10,'  ------ Unsymetrised Occupation'
   call wrtout(std_out,  message,'COLL')

   call print_matlu(xocc_check,natom,dtset%pawprtvol)

!  print norms
   write(message,'(2a,i4)')  ch10,'  ------ Unsymetrised Norm'
   call wrtout(std_out,  message,'COLL')

   call print_matlu(xnorm_check,natom,dtset%pawprtvol)
 end if

!symetrise and print occupations
 call sym_matlu(cryst_struc,loc_occ_check,pawang)

 write(message,'(2a,i4)')  ch10,'  ------ Symetrised Occupation'
 call wrtout(std_out,  message,'COLL')

 call print_matlu(loc_occ_check,natom,dtset%pawprtvol)

!symetrise and print norms
 call sym_matlu(cryst_struc,loc_norm_check,pawang)

 write(message,'(2a,i4)')  ch10,'  ------ Symetrised Norm'
 call wrtout(std_out,  message,'COLL')

 call print_matlu(loc_norm_check,natom,dtset%pawprtvol)

!deallocations
 do iatom=1,natom
   lda_occup%matlu(iatom)%mat=loc_occ_check(iatom)%mat
 end do

!Tests density matrix LDA+U and density matrix computed here.
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
&             cmplx(paw_ij(iatom)%noccmmp(1,im,im1,idijeff),paw_ij(iatom)%noccmmp(2,im,im1,idijeff),kind=dp)
             if(my_nspinor==1) matlu_temp(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=&
&             cmplx(paw_ij(iatom)%noccmmp(1,im,im1,idijeff),zero,kind=dp)
           end do
         end do
       end do
     end if
   end do
   if(paw_dmft%dmftcheck==2) option=1
   if(paw_dmft%dmftcheck<=1) option=0
   call diff_matlu("LDA+U density matrix from INPUT wfk",&
&   "Direct calculation of density matrix with psichi from DIAGONALIZED wfk",&
&   matlu_temp,lda_occup%matlu,natom,option,tol3,ierrr) !tol1 tol2 tol3
   if(ierrr==-1) then
     write(message,'(10a)') ch10,&
&     '    -> These two quantities should agree if three conditions are fulfilled',ch10,&
&     '         -  input wavefunctions come from the same Hamiltonien (e.g LDA/GGA)',ch10,&
&     '         -  dmatpuopt is equal to 1',ch10,&
&     '         -  all valence states are in the valence',ch10,&
&     '    (for experts users: it is not compulsary that these conditions are fulfilled)'
     call wrtout(std_out,message,'COLL')
   end if
!  write(message,'(2a)') ch10,&
!  &   '  ***** => Calculations of density matrices with projections and in LDA+U are coherent****'
!  call wrtout(std_out,message,'COLL')

   call destroy_matlu(matlu_temp,natom)
   ABI_DATATYPE_DEALLOCATE(matlu_temp)
 else
   write(message,'(2a)') ch10,&
&   '  Warning: Consistency of density matrices computed from projection has not been checked: use dmftcheck>=2 '
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

 if(present(nbandkss)) then
   if(me.eq.0.and.nbandkss/=0) then
!     opt_renorm=1 ! if ucrpa==1, no need for individual orthonormalization
     opt_renorm=3
     if(dtset%ucrpa>=1) opt_renorm=2
     call psichi_renormalization(cryst_struc,paw_dmft,pawang,opt=opt_renorm)
     call psichi_print(dtset,cryst_struc%nattyp,cryst_struc%ntypat,nkpt,my_nspinor,&
&     nsppol,paw_dmft,pawtab,psps,t2g)
   end if ! proc=me
 end if

!*********************
!deallocate(psichi)
!DEBUG
!write(std_out,*)' datafordmft : exit'
!stop
!ENDDEBUG
 CONTAINS
!===========================================================
!!***

!!****f* datafordmft/psichi_print
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
&nsppol,paw_dmft,pawtab,psps,t2g)

 use m_profiling_abi
 use m_io_tools,  only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psichi_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,my_nspinor,nsppol,ntypat
!arrays
 integer, intent(in) :: nattyp(ntypat)
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 logical t2g
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_dmft_type), intent(in) :: paw_dmft
!Local variables ------------------------------------
 integer :: ibg,isppol,ikpt,iband,ibandc,ispinor,icat,itypat,lmn_size
 integer :: iat,iatom,jj1,ilmn,m1,nband_k,unt
 integer :: m1_t2g,ll
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
               if(pawtab(dtset%typat(iatom))%lpawu.ne.-1) then
!                chinorm=(pawtab(itypat)%phiphjint(1))
                 chinorm=1.d0
!                write(std_out,*) isppol,ikpt,iband,ispinor,iat
                 m1_t2g=0
                 do ilmn=1,lmn_size
!                  write(std_out,*) ilmn
!                  ------------ Select l=lpawu.  ---------------------------------------
                   if (psps%indlmn(1,ilmn,itypat)==pawtab(dtset%typat(iatom))%lpawu.and.psps%indlmn(3,ilmn,itypat)==1) then
!                    ------------ Check that the band is choosen (within nbandi and nbandf)
                     if(paw_dmft%band_in(iband)) then
                       jj1=jj1+1
                       if(jj1>(2*pawtab(dtset%typat(iatom))%lpawu+1)) then
                         write(message,'(a,a,i4,i5,i4)') ch10," Error 2 in datafordmft",jj1,pawtab(dtset%typat(iatom))%lpawu
                         call wrtout(std_out,  message,'COLL')
                         MSG_ERROR("Aborting now")
                       end if ! jj1
!                      if(jj1>pawtab(dtset%typat(iatom))%nproju*(2*lpawu+1)) then
!                      write(message,'(a,a,a,a)')  ch10,&
!                      &                         'BUG: jj1 is not correct in datafordmft psichi_print',ch10,&
!                      &                         'Action: CONTACT Abinit group'
!                      call wrtout(std_out,  message,'COLL')
!                      call leave_new('COLL')
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

!!****f* datafordmft/psichi_check
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

 use m_profiling_abi

 use m_matlu, only: matlu_type,init_matlu,sym_matlu

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psichi_check'
!End of the abilint section

 implicit none

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
