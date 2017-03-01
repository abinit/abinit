!{\src2tex{textfont=tt}}
!!****f* abinit/prep_nonlop
!! NAME
!! prep_nonlop
!!
!! FUNCTION
!! this routine prepares the data to the call of nonlop.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (FB,MT,GZ,MD,FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  choice: chooses possible output:
!!    choice=1 => a non-local energy contribution
!!          =2 => a gradient with respect to atomic position(s)
!!          =3 => a gradient with respect to strain(s)
!!          =23=> a gradient with respect to atm. pos. and strain(s)
!!          =4 => a 2nd derivative with respect to atomic pos.
!!          =24=> a gradient and 2nd derivative with respect to atomic pos.
!!          =5 => a gradient with respect to k wavevector
!!          =6 => 2nd derivatives with respect to strain and atm. pos.
!!          =7 => no operator, just projections
!!  blocksize= size of block for FFT
!!  cpopt=flag defining the status of cwaveprj=<Proj_i|Cnk> scalars (see below, side effects)
!!  cwavef(2,npw*my_nspinor*blocksize)=planewave coefficients of wavefunction.
!!  gvnlc=matrix elements <G|Vnonlocal|C>
!!  hamk <type(gs_hamiltonian_type)>=data defining the Hamiltonian at a given k (NL part involved here)
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2),
!!                        - k point direction in the case (choice=5,signs=2)
!!                        - strain component (1:6) in the case (choice=2,signs=2) or (choice=6,signs=1)
!!  lambdablock(blocksize)=factor to be used when computing (Vln-lambda.S) - only for paw_opt=2
!!  mpi_enreg=informations about mpi parallelization
!!  nnlout=dimension of enlout (when signs=1):
!!  ntypat=number of types of atoms in cell
!!  paw_opt= define the nonlocal operator concerned with
!!  signs= if 1, get contracted elements (energy, forces, stress, ...)
!!         if 2, applies the non-local operator to a function in reciprocal space
!!  tim_nonlop=timing code of the calling routine (can be set to 0 if not attributed)
!!
!! OUTPUT
!!  ==== if (signs==1) ====
!!  enlout_block(nnlout)=
!!    if paw_opt==0, 1 or 2: contribution of this block of states to the nl part of various properties
!!    if paw_opt==3:         contribution of this block of states to <c|S|c>  (where S=overlap when PAW)
!!  ==== if (signs==2) ====
!!    if paw_opt==0, 1, 2 or 4:
!!       gvnlc(2,my_nspinor*npw)=result of the aplication of the nl operator
!!                        or one of its derivative to the input vect.
!!    if paw_opt==3 or 4:
!!       gsc(2,my_nspinor*npw*(paw_opt/3))=result of the aplication of (I+S)
!!                        to the input vect. (where S=overlap when PAW)
!!
!! SIDE EFFECTS
!!  ==== ONLY IF useylm=1
!!  cwaveprj(natom,my_nspinor) <type(pawcprj_type)>=projected input wave function |c> on non-local projector
!!                                  =<p_lmn|c> and derivatives
!!                                  Treatment depends on cpopt parameter:
!!                     if cpopt=-1, <p_lmn|in> (and derivatives)
!!                                  have to be computed (and not saved)
!!                     if cpopt= 0, <p_lmn|in> have to be computed and saved
!!                                  derivatives are eventually computed but not saved
!!                     if cpopt= 1, <p_lmn|in> and first derivatives have to be computed and saved
!!                                  other derivatives are eventually computed but not saved
!!                     if cpopt= 2  <p_lmn|in> are already in memory;
!!                                  only derivatives are computed here and not saved
!! (if useylm=0, should have cpopt=-1)
!!
!! PARENTS
!!      energy,forstrnps,m_invovl,vtowfk
!!
!! CHILDREN
!!      dcopy,nonlop,prep_index_wavef_bandpp,prep_sort_wavef_spin,timab
!!      xmpi_allgather,xmpi_alltoallv
!!
!! NOTES
!!  cprj (as well as cg) is distributed over band processors.
!!  Only the mod((iband-1)/mpi_enreg%bandpp,mpi_enreg%nproc_band) projected WFs are stored on each proc.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prep_nonlop(choice,cpopt,cwaveprj,enlout_block,hamk,idir,lambdablock,&
&                      blocksize,mpi_enreg,nnlout,paw_opt,signs,gsc,&
&                      tim_nonlop,cwavef,gvnlc,already_transposed)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_bandfft_kpt, only : bandfft_kpt,bandfft_kpt_get_ikpt
 use m_hamiltonian, only : gs_hamiltonian_type
 use m_pawcprj, only : pawcprj_type
 use m_nonlop

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_nonlop'
 use interfaces_18_timing
 use interfaces_66_wfs, except_this_one => prep_nonlop
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: blocksize,choice,cpopt,idir,signs,nnlout,paw_opt
 logical,optional,intent(in) :: already_transposed
 real(dp),intent(in) :: lambdablock(blocksize)
 real(dp),intent(out) :: enlout_block(nnlout*blocksize),gvnlc(:,:),gsc(:,:)
 real(dp),intent(inout) :: cwavef(:,:)
 type(gs_hamiltonian_type),intent(in) :: hamk
 type(mpi_type),intent(inout) :: mpi_enreg
 type(pawcprj_type),intent(inout) :: cwaveprj(:,:)

!Local variables-------------------------------
!scalars
 integer :: bandpp,ier,ikpt_this_proc,my_nspinor,ndatarecv,nproc_band,npw,nspinortot
 integer :: old_me_g0,spaceComm=0,tim_nonlop
 logical :: do_transpose
 character(len=500) :: msg
!arrays
 integer,allocatable :: index_wavef_band(:)
 integer,  allocatable :: rdisplsloc(:),recvcountsloc(:),sdisplsloc(:),sendcountsloc(:)
 integer,ABI_CONTIGUOUS  pointer :: rdispls(:),recvcounts(:),sdispls(:),sendcounts(:)
 real(dp) :: lambda_nonlop(mpi_enreg%bandpp),tsec(2)
 real(dp), allocatable :: cwavef_alltoall2(:,:),gvnlc_alltoall2(:,:),gsc_alltoall2(:,:)
 real(dp), allocatable :: cwavef_alltoall1(:,:),gvnlc_alltoall1(:,:),gsc_alltoall1(:,:)
 real(dp),allocatable :: enlout(:)

! *************************************************************************

 DBG_ENTER('COLL')

 call timab(570,1,tsec)

 do_transpose = .true.
 if(present(already_transposed)) then
   if(already_transposed) do_transpose = .false.
 end if

 nproc_band = mpi_enreg%nproc_band
 bandpp     = mpi_enreg%bandpp
 spaceComm=mpi_enreg%comm_fft
 if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_band
 my_nspinor=max(1,hamk%nspinor/mpi_enreg%nproc_spinor)
 nspinortot=hamk%nspinor

!Check sizes
 npw=hamk%npw_k;if (.not.do_transpose) npw=hamk%npw_fft_k
 if (size(cwavef)/=2*npw*my_nspinor*blocksize) then
   msg = 'Incorrect size for cwavef!'
   MSG_BUG(msg)
 end if
 if(choice/=0.and.signs==2) then
   if (paw_opt/=3) then
     if (size(gvnlc)/=2*npw*my_nspinor*blocksize) then
       msg = 'Incorrect size for gvnlc!'
       MSG_BUG(msg)
     end if
   end if
   if(paw_opt>=3) then
     if (size(gsc)/=2*npw*my_nspinor*blocksize) then
       msg = 'Incorrect size for gsc!'
       MSG_BUG(msg)
     end if
   end if
 end if
 if(cpopt>=0) then
   if (size(cwaveprj)/=hamk%natom*my_nspinor*mpi_enreg%bandpp) then
     msg = 'Incorrect size for cwaveprj!'
     MSG_BUG(msg)
   end if
 end if

 ABI_ALLOCATE(sendcountsloc,(nproc_band))
 ABI_ALLOCATE(sdisplsloc   ,(nproc_band))
 ABI_ALLOCATE(recvcountsloc,(nproc_band))
 ABI_ALLOCATE(rdisplsloc   ,(nproc_band))

 ikpt_this_proc=bandfft_kpt_get_ikpt()

 recvcounts   => bandfft_kpt(ikpt_this_proc)%recvcounts(:)
 sendcounts   => bandfft_kpt(ikpt_this_proc)%sendcounts(:)
 rdispls      => bandfft_kpt(ikpt_this_proc)%rdispls   (:)
 sdispls      => bandfft_kpt(ikpt_this_proc)%sdispls   (:)
 ndatarecv    =  bandfft_kpt(ikpt_this_proc)%ndatarecv

 ABI_ALLOCATE(cwavef_alltoall2,(2,ndatarecv*my_nspinor*bandpp))
 ABI_ALLOCATE(gsc_alltoall2,(2,ndatarecv*my_nspinor*(paw_opt/3)*bandpp))
 ABI_ALLOCATE(gvnlc_alltoall2,(2,ndatarecv*my_nspinor*bandpp))

 if(do_transpose .and. (bandpp/=1 .or. (bandpp==1 .and. mpi_enreg%paral_spinor==0.and.nspinortot==2)))then
   ABI_ALLOCATE(cwavef_alltoall1,(2,ndatarecv*my_nspinor*bandpp))
   if(signs==2)then
     if (paw_opt/=3) then
       ABI_ALLOCATE(gvnlc_alltoall1,(2,ndatarecv*my_nspinor*bandpp))
     end if
     if (paw_opt==3.or.paw_opt==4) then
       ABI_ALLOCATE(gsc_alltoall1,(2,ndatarecv*my_nspinor*bandpp))
     end if
   end if
 end if

 ABI_ALLOCATE(enlout,(nnlout*bandpp))
 enlout = zero

 recvcountsloc(:)=recvcounts(:)*2*my_nspinor*bandpp
 rdisplsloc(:)=rdispls(:)*2*my_nspinor*bandpp
 sendcountsloc(:)=sendcounts(:)*2*my_nspinor
 sdisplsloc(:)=sdispls(:)*2*my_nspinor

 if(do_transpose) then
   call timab(581,1,tsec)
   if(bandpp/=1 .or. (bandpp==1 .and. mpi_enreg%paral_spinor==0.and.nspinortot==2))then
     call xmpi_alltoallv(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall1,&
&     recvcountsloc,rdisplsloc,spaceComm,ier)
   else
     call xmpi_alltoallv(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall2,&
&     recvcountsloc,rdisplsloc,spaceComm,ier)
   end if
   call timab(581,2,tsec)
 else
   ! Here, we cheat, and use DCOPY to bypass some compiler's overzealous bound-checking
   ! (ndatarecv*my_nspinor*bandpp might be greater than the declared size of cwavef)
   call DCOPY(2*ndatarecv*my_nspinor*bandpp,cwavef,1,cwavef_alltoall2,1)
 end if

 if(hamk%istwf_k==2) then
   old_me_g0=mpi_enreg%me_g0
   if (mpi_enreg%me_fft==0) then
     mpi_enreg%me_g0=1
   else
     mpi_enreg%me_g0=0
   end if
 end if

!=====================================================================
 if (bandpp==1) then


   if (do_transpose .and. mpi_enreg%paral_spinor==0.and.nspinortot==2) then !Sort WF by spin
     call prep_sort_wavef_spin(nproc_band,my_nspinor,ndatarecv,recvcounts,rdispls,index_wavef_band)
     cwavef_alltoall2(:, :)=cwavef_alltoall1(:,index_wavef_band)
   end if

   if (paw_opt==2) lambda_nonlop(1)=lambdablock(mpi_enreg%me_band+1)
   call nonlop(choice,cpopt,cwaveprj,enlout,hamk,idir,lambda_nonlop,mpi_enreg,1,nnlout,paw_opt,&
&   signs,gsc_alltoall2,tim_nonlop,cwavef_alltoall2,gvnlc_alltoall2)

   if (do_transpose .and. mpi_enreg%paral_spinor == 0 .and. nspinortot==2.and.signs==2) then
     if (paw_opt/=3) gvnlc_alltoall1(:,index_wavef_band)=gvnlc_alltoall2(:,:)
     if (paw_opt==3.or.paw_opt==4) gsc_alltoall1(:,index_wavef_band)=gsc_alltoall2(:,:)
   end if

 else   ! bandpp/=1

!  -------------------------------------------------------------
!  Computation of the index used to sort the waves functions below bandpp
!  -------------------------------------------------------------
   if(do_transpose) then
     call prep_index_wavef_bandpp(nproc_band,bandpp,&
&     my_nspinor,ndatarecv,recvcounts,rdispls,index_wavef_band)

!  -------------------------------------------------------
!  Sorting of the waves functions below bandpp
!  -------------------------------------------------------
     cwavef_alltoall2(:,:) = cwavef_alltoall1(:,index_wavef_band)
   end if

!  -------------------------------------------------------
!  Call nonlop
!  -------------------------------------------------------
   if(paw_opt == 2) lambda_nonlop(1:bandpp) = lambdablock((mpi_enreg%me_band*bandpp)+1:((mpi_enreg%me_band+1)*bandpp))
   call nonlop(choice,cpopt,cwaveprj,enlout,hamk,idir,lambda_nonlop,mpi_enreg,bandpp,nnlout,paw_opt,&
&   signs,gsc_alltoall2,tim_nonlop,cwavef_alltoall2,gvnlc_alltoall2)

!  -----------------------------------------------------
!  Sorting of waves functions below the processors
!  -----------------------------------------------------
   if(do_transpose.and.signs==2) then
     if (paw_opt/=3) gvnlc_alltoall1(:,index_wavef_band)=gvnlc_alltoall2(:,:)
     if (paw_opt==3.or.paw_opt==4) gsc_alltoall1(:,index_wavef_band)=gsc_alltoall2(:,:)
   end if

 end if

!=====================================================================
!  -------------------------------------------------------
!  Deallocation
!  -------------------------------------------------------
 if (allocated(index_wavef_band)) then
   ABI_DEALLOCATE(index_wavef_band)
 end if

!Transpose the gsc_alltoall or gvlnc_alltoall tabs
!according to the paw_opt and signs values
 if(do_transpose) then
   if (signs==2) then
     call timab(581,1,tsec)
     if(bandpp/=1 .or. (bandpp==1 .and. mpi_enreg%paral_spinor==0.and.nspinortot==2))then
       if (paw_opt/=3) then
         call xmpi_alltoallv(gvnlc_alltoall1,recvcountsloc,rdisplsloc,gvnlc,&
&         sendcountsloc,sdisplsloc,spaceComm,ier)
       end if
       if (paw_opt==3.or.paw_opt==4) then
         call xmpi_alltoallv(gsc_alltoall1,recvcountsloc,rdisplsloc,gsc,&
&         sendcountsloc,sdisplsloc,spaceComm,ier)
       end if
     else
       if (paw_opt/=3) then
         call xmpi_alltoallv(gvnlc_alltoall2,recvcountsloc,rdisplsloc,gvnlc,&
&         sendcountsloc,sdisplsloc,spaceComm,ier)
       end if
       if (paw_opt==3.or.paw_opt==4) then
         call xmpi_alltoallv(gsc_alltoall2,recvcountsloc,rdisplsloc,gsc,&
&         sendcountsloc,sdisplsloc,spaceComm,ier)
       end if
     end if
     call timab(581,2,tsec)
   end if
 else
   ! TODO check other usages, maybe
   call DCOPY(2*ndatarecv*my_nspinor*bandpp, gsc_alltoall2, 1, gsc, 1)
 end if
 if (hamk%istwf_k==2) mpi_enreg%me_g0=old_me_g0

 if (nnlout>0) then
   call xmpi_allgather(enlout,nnlout*bandpp,enlout_block,spaceComm,ier)
 end if
 ABI_DEALLOCATE(enlout)
 ABI_DEALLOCATE(sendcountsloc)
 ABI_DEALLOCATE(sdisplsloc)
 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)
 ABI_DEALLOCATE(cwavef_alltoall2)
 ABI_DEALLOCATE(gvnlc_alltoall2)
 ABI_DEALLOCATE(gsc_alltoall2)
 if(do_transpose .and. (bandpp/=1 .or. (bandpp==1 .and. mpi_enreg%paral_spinor==0.and.nspinortot==2)))then
   ABI_DEALLOCATE(cwavef_alltoall1)
   if(signs==2)then
     if (paw_opt/=3) then
       ABI_DEALLOCATE(gvnlc_alltoall1)
     end if
     if (paw_opt==3.or.paw_opt==4) then
       ABI_DEALLOCATE(gsc_alltoall1)
     end if
   end if
 end if

 call timab(570,2,tsec)

 DBG_EXIT('COLL')

end subroutine prep_nonlop
!!***
