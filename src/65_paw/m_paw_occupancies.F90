!{\src2tex{textfont=tt}}
!!****m* m_paw_occupancies/m_paw_occupancies
!! NAME
!!  m_paw_occupancies
!!
!! FUNCTION
!!  This module contains routines related to the computation of PAW on-site occupancies (rhoij).
!!
!! COPYRIGHT
!! Copyright (C) 2018-2019 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_occupancies

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi

 use defs_abitypes, only : MPI_type
 use m_pawtab,     only : pawtab_type
 use m_pawrhoij,   only : pawrhoij_type,pawrhoij_init_unpacked,pawrhoij_mpisum_unpacked, &
&                         pawrhoij_alloc,pawrhoij_free,pawrhoij_inquire_dim
 use m_pawcprj,    only : pawcprj_type,pawcprj_alloc,pawcprj_get, &
&                         pawcprj_gather_spin, pawcprj_free, pawcprj_mpi_send, &
&                         pawcprj_mpi_recv, pawcprj_copy, pawcprj_unpack, pawcprj_pack
 use m_paw_io,     only : pawio_print_ij
 use m_paral_atom, only : get_my_atmtab,free_my_atmtab
 use m_paw_dmft,   only : paw_dmft_type
 use m_mpinfo,     only : proc_distrb_cycle

 implicit none

 private

!public procedures.
 public :: pawmkrhoij  ! Compute the PAW occupancies rhoij
 public :: pawaccrhoij ! Accumulate the contribution of one band to the PAW occupancies rhoij
 public :: initrhoij   ! Initialize the PAW occupancies rhoij from atomic data

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_occupancies/pawmkrhoij
!!
!! NAME
!! pawmkrhoij
!!
!! FUNCTION
!! Calculate the PAW quantities rhoij (augmentation occupancies)
!! Remember:for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cprj(natom,mcprj)= wave functions projected with non-local projectors:
!!                     cprj_nk(i)=<p_i|Cnk> where p_i is a non-local projector.
!!  dimcprj(natom)=array of dimensions of array cprj (ordered by atom-type)
!!  istwfk(nkpt)=parameter that describes the storage of wfs
!!  kptopt=option for the generation of k points
!!  mband=maximum number of bands
!!  mband_cprj=maximum number of bands used in the dimensioning of cprj array (usually mband/nproc_band)
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mkmem =number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell
!!  nband=number of bands for all k points
!!  nkpt=number of k points
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol)=occupation number for each band for each k
!!  paral_kgb=Flag related to the kpoint-band-fft parallelism
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  unpaw=unit number for cprj PAW data (if used)
!!  wtk(nkpt)=weight assigned to each k point
!!
!! SIDE EFFECTS
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  On input: arrays dimensions
!!  On output:
!!    pawrhoij(:)%rhoij_(lmn2_size,nspden)=
!!          Sum_{n,k} {occ(n,k)*conjugate[cprj_nk(ii)].cprj_nk(jj)} (non symetrized)
!!
!! PARENTS
!!      afterscfloop,scfcv,vtorho
!!
!! CHILDREN
!!      pawaccrhoij,pawcprj_alloc,pawcprj_free,pawcprj_gather_spin,pawcprj_get
!!      pawrhoij_free,pawrhoij_init_unpacked
!!      pawrhoij_mpisum_unpacked,wrtout
!!
!! NOTES
!!  The cprj are distributed over band processors.
!!  Only the mod((iband-1)/mpi_enreg%bandpp,mpi_enreg%nproc_band) projectors
!!  are stored on each proc.
!!
!! SOURCE

 subroutine pawmkrhoij(atindx,atindx1,cprj,dimcprj,istwfk,kptopt,mband,mband_cprj,mcprj,mkmem,mpi_enreg,&
&                      natom,nband,nkpt,nspinor,nsppol,occ,paral_kgb,paw_dmft,pawrhoij,unpaw,usewvl,wtk)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: kptopt,mband,mband_cprj,mcprj,mkmem,natom,nkpt,nspinor,nsppol
 integer,intent(in) :: paral_kgb,unpaw,usewvl
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),dimcprj(natom),istwfk(nkpt)
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),wtk(nkpt)
 type(pawcprj_type),target,intent(in) :: cprj(natom,mcprj)
 type(paw_dmft_type),intent(in) :: paw_dmft
 type(pawrhoij_type),intent(inout),target:: pawrhoij(:)

!Local variables ---------------------------------------
!scalars
 integer,parameter :: max_nband_cprj=100
 integer :: bdtot_index,cplex
 integer :: iatom,iatom_tot,ib,ib1,iband,ibc1,ibg,ib_this_proc,ierr
 integer :: ikpt,iorder_cprj,isppol,jb_this_proc,jbg,me,my_nspinor,nband_k,nband_k_cprj
 integer :: nbandc1,nband_k_cprj_read,nband_k_cprj_used,nprocband,nrhoij
 integer :: option,spaceComm,use_nondiag_occup_dmft
 logical :: locc_test,paral_atom,usetimerev
 integer :: ib1_this_proc, ib_loop, proc_sender, proc_recver
 real(dp) :: wtk_k
 character(len=500) :: msg

!arrays
 integer :: n2buff
 integer, allocatable :: req_correl(:,:,:)
 real(dp) :: occup(2)
 real(dp) ABI_ASYNC, allocatable :: buffer_cprj_correl(:,:,:)
 character(len=8),parameter :: dspin(6)=(/"up      ","down    ","dens (n)","magn (x)","magn (y)","magn (z)"/)
 type(pawcprj_type),allocatable :: cprj_tmp(:,:),cwaveprj(:,:),cwaveprjb(:,:)
 type(pawcprj_type),pointer :: cprj_ptr(:,:)
 type(pawrhoij_type),pointer :: pawrhoij_all(:)

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(mkmem/=0,"mkmem==0 not supported anymore!")

!Init MPI data
! spaceComm=mpi_enreg%comm_cell
! if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_kpt
 spaceComm=mpi_enreg%comm_kpt
 me=mpi_enreg%me_kpt

!Check size of cprj
 my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)
 if (mcprj/=my_nspinor*mband_cprj*mkmem*nsppol) then
   msg=' wrong size for cprj !'
   MSG_BUG(msg)
 end if

!Check if cprj is distributed over bands
 nprocband=(mband/mband_cprj)
 if (paral_kgb==1.and.nprocband/=mpi_enreg%nproc_band) then
   msg='mband/mband_cprj must be equal to nproc_band!'
   MSG_BUG(msg)
 end if

 if( usewvl==1 .and. (nprocband/=1)) then
   write(msg,'(2a)') ch10,&
&   'Parallelization over bands is not compatible with WAVELETS!'
   MSG_ERROR(msg)
 end if

!Initialise and check dmft variables
 if(paw_dmft%use_sc_dmft/=0) then
   nbandc1=paw_dmft%mbandc
 else
   nbandc1=1
 end if

!Size of pawrhoij datastructure
 nrhoij=size(pawrhoij)

!Check if pawrhoij is distributed over atomic sites
 paral_atom=(nrhoij/=natom.and.mpi_enreg%nproc_atom>1)
 if (paral_atom.and.nrhoij/=mpi_enreg%my_natom) then
   msg='Size of pawrhoij should be natom or my_natom!'
   MSG_BUG(msg)
 end if

!Allocate temporary cwaveprj storage
 ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,nspinor))
 call pawcprj_alloc(cwaveprj,0,dimcprj)
 if(paw_dmft%use_sc_dmft/=0) then
   ABI_DATATYPE_ALLOCATE(cwaveprjb,(natom,nspinor))
   call pawcprj_alloc(cwaveprjb,0,dimcprj)
 end if

 if (paw_dmft%use_sc_dmft /= 0 .and. mpi_enreg%paral_kgb /= 0) then
   if(paw_dmft%use_bandc(mpi_enreg%me_band+1)) then
     n2buff = nspinor*sum(dimcprj)
     ABI_ALLOCATE(buffer_cprj_correl,(2,n2buff,nbandc1))
     ABI_ALLOCATE(req_correl,(nbandc1, nkpt, nsppol))
     req_correl(:,:,:) = 0
   end if
 end if

!Initialize temporary file (if used)
 iorder_cprj=0

!Build and initialize unpacked rhoij (to be computed here)
 call pawrhoij_init_unpacked(pawrhoij)

!If pawrhoij is MPI-distributed over atomic sites, gather it
 if (paral_atom) then
   ABI_DATATYPE_ALLOCATE(pawrhoij_all,(natom))
 else
   pawrhoij_all => pawrhoij
 end if

!LOOP OVER SPINS
 option=1
 usetimerev=(kptopt>0.and.kptopt<3)
 bdtot_index=0;ibg=0;jbg=0
 do isppol=1,nsppol

!  LOOP OVER k POINTS
   do ikpt=1,nkpt

     nband_k=nband(ikpt+(isppol-1)*nkpt)
     nband_k_cprj=nband_k/nprocband
     wtk_k=wtk(ikpt)

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k
       cycle
     end if

     cplex=2;if (istwfk(ikpt)>1) cplex=1

!    In case of spinors parallelism, need some extra storage
     if (mpi_enreg%paral_spinor==1) then
       nband_k_cprj_used=min(max_nband_cprj,nband_k_cprj)
       ABI_DATATYPE_ALLOCATE(cprj_tmp,(natom,my_nspinor*nband_k_cprj_used))
       ABI_DATATYPE_ALLOCATE(cprj_ptr,(natom,   nspinor*nband_k_cprj_used))
       call pawcprj_alloc(cprj_tmp,0,dimcprj)
       call pawcprj_alloc(cprj_ptr,0,dimcprj)
     else
       cprj_ptr => cprj
     end if

!    In case of band parallelism combined with self consistent DMFT, need to
!    exchange bands cprj
     if (paw_dmft%use_sc_dmft /= 0 .and. mpi_enreg%paral_kgb /= 0) then
       if (paw_dmft%use_bandc(mpi_enreg%me_band+1)) then
! only proc using correlated band have to do this
         do ibc1=1,nbandc1
           proc_sender = paw_dmft%bandc_proc(ibc1)
           if(proc_sender == mpi_enreg%me_band) then

!            get the index of band local to this proc
             ib1 = paw_dmft%include_bands(ibc1)
             ib1_this_proc = 0
             do ib_loop=1,ib1-1
               if (mod((ib_loop-1)/mpi_enreg%bandpp,mpi_enreg%nproc_band) == mpi_enreg%me_band) then
                 ib1_this_proc = ib1_this_proc+1
               end if
             end do
             ib1_this_proc = ib1_this_proc+1

!            extract the band
             call pawcprj_get(atindx1,cwaveprjb,cprj_ptr,natom,ib1_this_proc,ibg,ikpt,&
&                             iorder_cprj,isppol,mband_cprj,mkmem,natom,1,nband_k_cprj,nspinor,nsppol,&
&                             unpaw,mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)

             call pawcprj_pack(dimcprj,cwaveprjb,buffer_cprj_correl(:,:,ibc1))
             do proc_recver=0,mpi_enreg%nproc_band-1
               if (proc_sender /= proc_recver .and. paw_dmft%use_bandc(proc_recver+1)) then
!                locc_test = At least one of the bands used by proc_recver have a non neglectable occnd
                 locc_test = .false.
                 do ib_loop=1,nbandc1
                   if(proc_recver == paw_dmft%bandc_proc(ib_loop)) then
                     ib = paw_dmft%include_bands(ib_loop)
                     locc_test = locc_test .or. (abs(paw_dmft%occnd(1,ib,ib1,ikpt,isppol))+&
&                                                abs(paw_dmft%occnd(2,ib,ib1,ikpt,isppol))>tol8)
                   end if
                 end do
                 if(locc_test) then
                   ! send to proc_recver
                   ierr = 0
                   call xmpi_isend(buffer_cprj_correl(:,:,ibc1),proc_recver,&
&                                  10000+ibc1+nbandc1*(ikpt+nsppol*isppol),mpi_enreg%comm_band,&
&                                  req_correl(ibc1,ikpt,isppol),ierr)
!                  force sending or buffering
                   call xmpi_wait(req_correl(ibc1,ikpt,isppol), ierr)
                 end if
               end if
             end do
           else
!            locc_test = At least one of the bands used by this proc have a non neglectable occnd
             locc_test = .false.
             do ib_loop=1,nbandc1
               if(mpi_enreg%me_band == paw_dmft%bandc_proc(ib_loop)) then
                 ib = paw_dmft%include_bands(ib_loop)
                 ib1 = paw_dmft%include_bands(ibc1)
                 locc_test = locc_test .or. (abs(paw_dmft%occnd(1,ib,ib1,ikpt,isppol))+&
&                                            abs(paw_dmft%occnd(2,ib,ib1,ikpt,isppol))>tol8)
               end if
             end do
             if(locc_test) then
               ! recv from proc_sender
               ierr = 0
               call xmpi_irecv(buffer_cprj_correl(:,:,ibc1),proc_sender,&
&                              10000+ibc1+nbandc1*(ikpt+nsppol*isppol),mpi_enreg%comm_band,&
&                              req_correl(ibc1,ikpt,isppol),ierr)
             end if
           end if
         end do
       end if
     end if

     ierr = 0
!    LOOP OVER BANDS
     ib_this_proc=0;jb_this_proc=0
     do ib=1,nband_k
       iband=bdtot_index+ib

!      Parallelization: treat only some bands
       if(xmpi_paral==1)then
         if (paral_kgb==1) then
           if (mod((ib-1)/mpi_enreg%bandpp,mpi_enreg%nproc_band)/=mpi_enreg%me_band) cycle
         else
           if (mpi_enreg%proc_distrb(ikpt,ib,isppol)/=me) cycle
         end if
       end if
       ib_this_proc=ib_this_proc+1

!      In case of spinors parallelism, gather cprj because we need both components together
!      We do that nband_k_cprj_used by nband_k_cprj_used bands
       if (mpi_enreg%paral_spinor==1) then
         jb_this_proc=jb_this_proc+1
         if (mod(jb_this_proc,nband_k_cprj_used)==1) then
           ib_this_proc=1
           nband_k_cprj_read=nband_k_cprj_used
           if (nband_k_cprj<jb_this_proc+nband_k_cprj_used-1) nband_k_cprj_read=nband_k_cprj-jb_this_proc+1
           call pawcprj_get(atindx1,cprj_tmp,cprj,natom,jb_this_proc,jbg,ikpt,iorder_cprj,isppol,&
&           mband_cprj,mkmem,natom,nband_k_cprj_read,nband_k_cprj,my_nspinor,nsppol,unpaw,&
&           mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
           call pawcprj_gather_spin(cprj_tmp,cprj_ptr,natom,nband_k_cprj_read,my_nspinor,nspinor,&
&           mpi_enreg%comm_spinor,ierr)
         end if
       end if

!      DMFT: LOOP ON ADDITIONAL BANDS
       do ibc1=1,nbandc1
!        check if dmft and occupations
!        write(std_out,*) 'ib,ibc1          ',ib,ibc1

!        if ib is not part a band correlated in dmft do not repeat the following
         if(ibc1 /= 1) then
           if (paw_dmft%use_sc_dmft == 0) then
             cycle
           else
             if (.not.(paw_dmft%band_in(ib))) cycle
           end if
         end if

!        DMFT stuff: extract cprj and occupations for additional band
         if(paw_dmft%use_sc_dmft /= 0) then
           if(paw_dmft%band_in(ib)) then
!            write(std_out,*) 'use_sc_dmft=1 ib,ib1',ib,ib1
!            write(std_out,*) 'ib, ib1          ',paw_dmft%band_in(ib),paw_dmft%band_in(ib1)

             ib1 = paw_dmft%include_bands(ibc1) ! indice reel de la bande

             use_nondiag_occup_dmft = 1
             locc_test = abs(paw_dmft%occnd(1,ib,ib1,ikpt,isppol))+abs(paw_dmft%occnd(2,ib,ib1,ikpt,isppol))>tol8

             occup(1) = paw_dmft%occnd(1,ib,ib1,ikpt,isppol)
             occup(2) = paw_dmft%occnd(2,ib,ib1,ikpt,isppol)

!            write(std_out,*) 'use_sc_dmft=1,band_in(ib)=1, ib,ibc1',ib,ib1,locc_test
!
             if (locc_test .or. mkmem == 0) then

               if (paral_kgb==1) then
!                cprj have already been extracted
                 if (paw_dmft%bandc_proc(ibc1) /= mpi_enreg%me_band) then
!                  if the band is not on this proc, wait for the recv to complete
                   ierr = 0
                   call xmpi_wait(req_correl(ibc1,ikpt,isppol), ierr)
                 end if
                 call pawcprj_unpack(dimcprj,cwaveprjb,buffer_cprj_correl(:,:,ibc1))
               else ! paral_kgb /= 0
                 call pawcprj_get(atindx1,cwaveprjb,cprj_ptr,natom,ib1,ibg,ikpt,iorder_cprj,isppol,&
&                                 mband_cprj,mkmem,natom,1,nband_k_cprj,nspinor,nsppol,unpaw,&
&                                 mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
               end if
             end if
           else  ! nbandc1=1
             use_nondiag_occup_dmft=0
             locc_test = (abs(occ(iband))>tol8)
             occup(1) = occ(iband)
           end if
         else  ! nbandc1=1
           use_nondiag_occup_dmft=0
           locc_test = (abs(occ(iband))>tol8)
           occup(1) = occ(iband)
         end if

!        Extract cprj for current band
!        Must read cprj when mkmem=0 (even if unused) to have right pointer inside _PAW file
         if (locc_test.or.mkmem==0) then
           call pawcprj_get(atindx1,cwaveprj,cprj_ptr,natom,ib_this_proc,ibg,ikpt,iorder_cprj,isppol,&
&                           mband_cprj,mkmem,natom,1,nband_k_cprj,nspinor,nsppol,unpaw,&
&                           mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         end if

!        Accumulate contribution from (occupied) current band
         if (locc_test) then
           if(use_nondiag_occup_dmft == 1) then
             call pawaccrhoij(atindx,cplex,cwaveprj,cwaveprjb,0,isppol,nrhoij,natom,&
&                             nspinor,occup(1),option,pawrhoij_all,usetimerev,wtk_k,&
&                             occ_k_2=occup(2))
           else
             call pawaccrhoij(atindx,cplex,cwaveprj,cwaveprj ,0,isppol,nrhoij,natom,&
&                             nspinor,occup(1),option,pawrhoij_all,usetimerev,wtk_k)
           end if
         end if
       end do ! ib1c
     end do ! ib

     if (mpi_enreg%paral_spinor==1) then
       call pawcprj_free(cprj_tmp)
       call pawcprj_free(cprj_ptr)
       ABI_DATATYPE_DEALLOCATE(cprj_tmp)
       ABI_DATATYPE_DEALLOCATE(cprj_ptr)
     else
       nullify(cprj_ptr)
     end if

     bdtot_index=bdtot_index+nband_k
     if (mkmem/=0) then
       if (mpi_enreg%paral_spinor==0) then
         ibg=ibg+   nspinor*nband_k_cprj
       else
         jbg=jbg+my_nspinor*nband_k_cprj
       end if
     end if

   end do ! ikpt
 end do ! isppol

!call xmpi_barrier(mpi_enreg%comm_band)

!deallocate temporary cwaveprj/cprj storage
 call pawcprj_free(cwaveprj)
 ABI_DATATYPE_DEALLOCATE(cwaveprj)

 if(paw_dmft%use_sc_dmft/=0) then
   call pawcprj_free(cwaveprjb)
   ABI_DATATYPE_DEALLOCATE(cwaveprjb)
 end if

 if (allocated(buffer_cprj_correl)) then
   ABI_DEALLOCATE(buffer_cprj_correl)
   ABI_DEALLOCATE(req_correl)
 end if

!MPI: need to exchange rhoij_ between procs
 if (paral_kgb==1.and.nprocband>1) then
   call pawrhoij_mpisum_unpacked(pawrhoij_all,spaceComm,comm2=mpi_enreg%comm_band)
 else
   call pawrhoij_mpisum_unpacked(pawrhoij_all,spaceComm)
 end if

!In case of distribution over atomic sites, dispatch rhoij
 if (paral_atom) then
   do iatom=1,nrhoij
     iatom_tot=mpi_enreg%my_atmtab(iatom)
     pawrhoij(iatom)%rhoij_(:,:)=pawrhoij_all(iatom_tot)%rhoij_(:,:)
   end do
   call pawrhoij_free(pawrhoij_all)
   ABI_DATATYPE_DEALLOCATE(pawrhoij_all)
 end if

 DBG_EXIT("COLL")

end subroutine pawmkrhoij
!!***

!----------------------------------------------------------------------

!!****f* m_paw_occupancies/pawaccrhoij
!!
!! NAME
!! pawaccrhoij
!!
!! FUNCTION
!! Accumulate the PAW quantities rhoij (augmentation occupancies)
!! or their 1-st order change or their gradient vs r
!! Add the contribution of a given k-point and band
!! Remember: for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (sorted-->random), inverse of atindx.
!!  cplex: if 1, WFs (or 1st-order WFs) are REAL, if 2, COMPLEX
!!  cwaveprj(natom,nspinor)= LEFT wave function at given n,k
!!                         projected with non-local projectors: cwaveprj=<p_i|Cnk>
!!  cwaveprj1(natom,nspinor)= RIGHT wave function at n,k,q
!!                          projected with non-local projectors: cwaveprj1=<p_i|C1nk,q>
!!                          * USED for RF  : C1nk is the first-order wave function
!!                          * USED for DMFT: C1nk is the RIGHT wave function
!!                          * NOT USED in usual GS case; can be set to cwaveprj in that case
!!  ipert=index of perturbation (RF only, i.e. option=2)
!!  isppol=index of current spin component
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell
!!  nspinor=number of spinorial components (on current proc)
!!  occ_k=occupation number for current band n,k
!!  option: choice of calculation:
!!          1: update rhoij (Ground-State)
!!          2: update 1-st order rhoij (Response Function) according to ipert
!!          3: update gradients of rhoij with respect to r,strain of both
!!  usetimerev=.TRUE. if time-reversal symmetry is used (WF(-k)=Conjg[WF(k)])
!!  wtk_k=weight assigned to current k-point
!!
!! SIDE EFFECTS
!!  pawrhoij(natom) <type(pawrhoij_type)>= GS: paw rhoij occupancies and related data
!!                                         RF: 1-st order paw rhoij occupancies and related data
!!  On output, has been updated with the contribution of current n,k
!!    === option=1:
!!        pawrhoij(:)%rhoij_(lmn2_size,nspden)=      (non symetrized)
!!            Sum_{n,k} {occ(n,k)*conjugate[cprj_nk(ii)].cprj_nk(jj)}
!!    === option=2:
!!        pawrhoij(:)%rhoij_(lmn2_size,nspden)=      (non symetrized)
!!            Sum_{n,k} {occ(n,k)*(conjugate[cprj_nk(ii)].cprj1_nk,q(jj)
!!                                 conjugate[cprj_nk(jj)].cprj1_nk,q(ii)}
!!          + Sum_{n,k} {occ(n,k)*(conjugate[dcprj_nk(ii)/dlambda].cprj_nk(jj)
!!                                +conjugate[cprj_nk(ii)].dcprj_nk(jj)/dlambda)}
!!    === option=3:
!!        pawrhoij(:)%grhoij(lmn2_size,mu,nspden)=   (non symetrized)
!!            Sum_{n,k} {occ(n,k)*(conjugate[dcprj_nk(ii)/dr_mu].cprj_nk(jj)
!!                                +conjugate[cprj_nk(ii)].dcprj_nk(jj)/dr_mu)}
!!
!! PARENTS
!!      d2frnl,dfpt_accrho,energy,pawmkrhoij,posdoppler,wfd_pawrhoij
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab
!!
!! SOURCE

 subroutine pawaccrhoij(atindx,cplex,cwaveprj,cwaveprj1,ipert,isppol,my_natom,natom,&
&                       nspinor,occ_k,option,pawrhoij,usetimerev,wtk_k,occ_k_2, &
&                       comm_atom,mpi_atmtab ) ! optional (parallelism)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,ipert,isppol,my_natom,natom,nspinor,option
 integer,optional,intent(in) :: comm_atom
 logical,intent(in) :: usetimerev
 real(dp),intent(in) :: occ_k,wtk_k
 real(dp),optional,intent(in) :: occ_k_2
!arrays
 integer,intent(in) :: atindx(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(pawcprj_type),intent(in) :: cwaveprj(natom,nspinor),cwaveprj1(natom,nspinor)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)

!Local variables ---------------------------------------
!scalars
 integer :: cplex_rhoij,iatm,iatom,iatom1,ilmn,iplex,iq0,j0lmn,jlmn,klmn,klmn_im,klmn_re
 integer :: mu,my_comm_atom,ncpgr,nspden_rhoij
 logical :: compute_impart,compute_impart_cplex,substract_diagonal
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: ro11_im,ro11_re,ro12_im,ro12_re,ro21_im,ro21_re,ro22_im,ro22_re,weight,weight_2
 character(len=500) :: message
!arrays
 integer,pointer :: my_atmtab(:)
 real(dp) :: cpi0(2,nspinor),cpi1(2,nspinor),cpj0(2,nspinor),cpj1(2,nspinor)
 real(dp) :: dcpi0(2,nspinor,9),dcpj0(2,nspinor,9)

! ***********************************************************************

 DBG_ENTER("COLL")

 if (my_natom==0) return

 ncpgr=0
 if (option==2.and.(ipert<=natom.or.ipert==natom+3.or.ipert==natom+4)) ncpgr=1
 if (option==3) ncpgr=cwaveprj(1,1)%ncpgr

!Tests
 if(option==2.and.(ipert==natom+1.or.ipert==natom+10.or.ipert==natom+11)) then
   message = 'Not relevant for ipert=natom+1 or ipert=natom+10 or ipert=natom+11!'
   MSG_BUG(message)
 end if
 if(option==2.and.cwaveprj(1,1)%ncpgr<ncpgr) then
   message = 'Error on cwaveprj1 factors derivatives!'
   MSG_BUG(message)
 end if
 if(option==3.and.cwaveprj(1,1)%ncpgr/=ncpgr) then
   message = 'Error on cwaveprj factors derivatives!'
   MSG_BUG(message)
 end if
 if (pawrhoij(1)%qphase==2.and.option/=2) then
   message = 'pawaccrhoij: qphase=2 only allowed with option=2 (1st-order rhoij)!'
   MSG_BUG(message)
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
& my_natom_ref=my_natom)

 weight=wtk_k*occ_k
 weight_2=zero
 if(present(occ_k_2)) then
   if (cplex==1) then
     !DMFT need complex cprj
     message='DMFT computation must be done with complex cprj. Check that istwfk = 1!'
     MSG_ERROR(message)
   end if
   weight_2=wtk_k*occ_k_2
 end if
 if (pawrhoij(1)%nspden==2.and.pawrhoij(1)%nsppol==1.and.nspinor==1) weight=half*weight
 if (pawrhoij(1)%nspden==2.and.pawrhoij(1)%nsppol==1.and.nspinor==1.and.present(occ_k_2)) weight_2=half*weight_2

 if (option==1) then

!  ==================================================================
!  === OPTION 1: Accumulate (n,k) contribution to rhoij =============
!  ==================================================================
   compute_impart=((.not.usetimerev).and.(pawrhoij(1)%cplex_rhoij==2))
   compute_impart_cplex=((compute_impart).and.(cplex==2))
   if (nspinor==1) then
     cplex_rhoij=pawrhoij(1)%cplex_rhoij
     if (cplex_rhoij==1) then
       do iatom=1,my_natom
         iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
         iatm=atindx(iatom1)
         do jlmn=1,pawrhoij(iatom)%lmn_size
           j0lmn=jlmn*(jlmn-1)/2
           cpj0(1:cplex,1)=cwaveprj(iatm,1)%cp(1:cplex,jlmn)
           do ilmn=1,jlmn
             klmn=j0lmn+ilmn
             cpi0(1:cplex,1)=cwaveprj1(iatm,1)%cp(1:cplex,ilmn)
             ro11_re=zero
             do iplex=1,cplex
               ro11_re=ro11_re+cpi0(iplex,1)*cpj0(iplex,1)
             end do
             pawrhoij(iatom)%rhoij_(klmn,isppol)=pawrhoij(iatom)%rhoij_(klmn,isppol)+weight*ro11_re
             if (present(occ_k_2)) then
               ro11_im=cpi0(1,1)*cpj0(2,1)-cpi0(2,1)*cpj0(1,1)
               pawrhoij(iatom)%rhoij_(klmn,isppol)=pawrhoij(iatom)%rhoij_(klmn,isppol)-weight_2*ro11_im
             end if
           end do
         end do
       end do
     else
       do iatom=1,my_natom
         iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
         iatm=atindx(iatom1)
         do jlmn=1,pawrhoij(iatom)%lmn_size
           j0lmn=jlmn*(jlmn-1)/2
           cpj0(1:cplex,1)=cwaveprj(iatm,1)%cp(1:cplex,jlmn)
           do ilmn=1,jlmn
             klmn=j0lmn+ilmn
             klmn_re=cplex_rhoij*(klmn-1)+1
             cpi0(1:cplex,1)=cwaveprj1(iatm,1)%cp(1:cplex,ilmn)
             ro11_re=zero
             do iplex=1,cplex
               ro11_re=ro11_re+cpi0(iplex,1)*cpj0(iplex,1)
             end do
             pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)+weight*ro11_re
             if (present(occ_k_2)) then
               ro11_im=cpi0(1,1)*cpj0(2,1)-cpi0(2,1)*cpj0(1,1)
               pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)-weight_2*ro11_im
             end if
             if (compute_impart_cplex) then
               klmn_im=klmn_re+1
               ro11_im=cpi0(1,1)*cpj0(2,1)-cpi0(2,1)*cpj0(1,1)
               pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight*ro11_im
               if (present(occ_k_2)) then
                 pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight_2*ro11_re
               end if
             end if
           end do
         end do
       end do
     end if
   else ! nspinor=2
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       cplex_rhoij=pawrhoij(iatom)%cplex_rhoij
       nspden_rhoij=pawrhoij(iatom)%nspden
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         cpj0(1:cplex,1)=cwaveprj(iatm,1)%cp(1:cplex,jlmn)
         cpj0(1:cplex,2)=cwaveprj(iatm,2)%cp(1:cplex,jlmn)
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           cpi0(1:cplex,1)=cwaveprj1(iatm,1)%cp(1:cplex,ilmn)
           cpi0(1:cplex,2)=cwaveprj1(iatm,2)%cp(1:cplex,ilmn)
           ro11_re=zero;ro22_re=zero
           ro12_re=zero;ro21_re=zero
           ro12_im=zero;ro21_im=zero
           do iplex=1,cplex
             ro11_re=ro11_re+cpi0(iplex,1)*cpj0(iplex,1)
             ro22_re=ro22_re+cpi0(iplex,2)*cpj0(iplex,2)
           end do
           pawrhoij(iatom)%rhoij_(klmn_re,1)=pawrhoij(iatom)%rhoij_(klmn_re,1)+weight*(ro11_re+ro22_re)
           if (nspden_rhoij>1) then
             do iplex=1,cplex
               ro12_re=ro12_re+cpi0(iplex,2)*cpj0(iplex,1)
               ro21_re=ro21_re+cpi0(iplex,1)*cpj0(iplex,2)
             end do
             pawrhoij(iatom)%rhoij_(klmn_re,4)=pawrhoij(iatom)%rhoij_(klmn_re,4)+weight*(ro11_re-ro22_re)
             pawrhoij(iatom)%rhoij_(klmn_re,2)=pawrhoij(iatom)%rhoij_(klmn_re,2)+weight*(ro12_re+ro21_re)
             if (cplex==2) then
               !Important note: the present implementation follows eq(15) in Hobbs et al, PRB 62, 11556(2000)
               ! rho^alpha,beta_ij = Sum[<Psi^beta|pi><pj|Psi^alpha]  (alpha and beta exponents inverted)
               ro12_im=cpi0(1,2)*cpj0(2,1)-cpi0(2,2)*cpj0(1,1)
               ro21_im=cpi0(1,1)*cpj0(2,2)-cpi0(2,1)*cpj0(1,2)
               pawrhoij(iatom)%rhoij_(klmn_re,3)=pawrhoij(iatom)%rhoij_(klmn_re,3)+weight*(ro21_im-ro12_im)
             end if
           end if
           if (present(occ_k_2)) then
             ro11_im=cpi0(1,1)*cpj0(2,1)-cpi0(2,1)*cpj0(1,1)
             ro22_im=cpi0(1,2)*cpj0(2,2)-cpi0(2,2)*cpj0(1,2)
             pawrhoij(iatom)%rhoij_(klmn_re,1)=pawrhoij(iatom)%rhoij_(klmn_re,1)+weight_2*(-ro11_im-ro22_im)
             pawrhoij(iatom)%rhoij_(klmn_re,2)=pawrhoij(iatom)%rhoij_(klmn_re,2)+weight_2*(-ro21_im-ro12_im)
             pawrhoij(iatom)%rhoij_(klmn_re,3)=pawrhoij(iatom)%rhoij_(klmn_re,3)+weight_2*(-ro12_re+ro21_re)
             pawrhoij(iatom)%rhoij_(klmn_re,4)=pawrhoij(iatom)%rhoij_(klmn_re,4)+weight_2*(-ro11_im+ro22_im)
           end if
           if (compute_impart) then
             klmn_im=klmn_re+1
             if (nspden_rhoij>1) pawrhoij(iatom)%rhoij_(klmn_im,3)=pawrhoij(iatom)%rhoij_(klmn_im,3)+weight*(ro12_re-ro21_re)
             if (cplex==2) then
               ro11_im=cpi0(1,1)*cpj0(2,1)-cpi0(2,1)*cpj0(1,1)
               ro22_im=cpi0(1,2)*cpj0(2,2)-cpi0(2,2)*cpj0(1,2)
               pawrhoij(iatom)%rhoij_(klmn_im,1)=pawrhoij(iatom)%rhoij_(klmn_im,1)+weight*(ro11_im+ro22_im)
               if (nspden_rhoij>1) then
                 pawrhoij(iatom)%rhoij_(klmn_im,4)=pawrhoij(iatom)%rhoij_(klmn_im,4)+weight*(ro11_im-ro22_im)
                 pawrhoij(iatom)%rhoij_(klmn_im,2)=pawrhoij(iatom)%rhoij_(klmn_im,2)+weight*(ro12_im+ro21_im)
               end if
               if (present(occ_k_2)) then
                 pawrhoij(iatom)%rhoij_(klmn_im,1)=pawrhoij(iatom)%rhoij_(klmn_im,1)+weight_2*( ro11_re+ro22_re)
                 pawrhoij(iatom)%rhoij_(klmn_im,2)=pawrhoij(iatom)%rhoij_(klmn_im,2)+weight_2*( ro21_re+ro12_re)
                 pawrhoij(iatom)%rhoij_(klmn_im,3)=pawrhoij(iatom)%rhoij_(klmn_im,3)+weight_2*(-ro12_im+ro21_im)
                 pawrhoij(iatom)%rhoij_(klmn_im,4)=pawrhoij(iatom)%rhoij_(klmn_im,4)+weight_2*( ro11_re-ro22_re)
               end if
             end if
           end if
         end do
       end do
     end do
   end if

 else if (option==2) then

!  ==================================================================
!  === OPTION 2: Accumulate (n,k) contribution to 1st-order rhoij ===
!  ==================================================================

!  Accumulate (n,k) contribution to rhoij1
!  due to derivative of wave-function
   compute_impart=(pawrhoij(1)%qphase==2)
   compute_impart_cplex=(compute_impart.and.(cplex==2))
   if (nspinor==1) then
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       cplex_rhoij=pawrhoij(iatom)%cplex_rhoij
       iq0=cplex_rhoij*pawrhoij(iatom)%lmn2_size
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         cpj0(1:2,1)=cwaveprj (iatm,1)%cp(1:2,jlmn)
         cpj1(1:2,1)=cwaveprj1(iatm,1)%cp(1:2,jlmn)
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           cpi0(1:2,1)=cwaveprj (iatm,1)%cp(1:2,ilmn)
           cpi1(1:2,1)=cwaveprj1(iatm,1)%cp(1:2,ilmn)
           ro11_re=zero
           do iplex=1,cplex
             ro11_re=ro11_re+cpi0(iplex,1)*cpj1(iplex,1)+cpj0(iplex,1)*cpi1(iplex,1)
           end do
           pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)+weight*ro11_re
           if (compute_impart_cplex) then
             klmn_im=klmn_re+iq0
             ro11_im=cpi0(1,1)*cpj1(2,1)-cpi0(2,1)*cpj1(1,1)+cpj0(1,1)*cpi1(2,1)-cpj0(2,1)*cpi1(1,1)
             pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight*ro11_im
           end if
         end do
       end do
     end do
   else ! nspinor=2
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       cplex_rhoij=pawrhoij(iatom)%cplex_rhoij
       nspden_rhoij=pawrhoij(iatom)%nspden
       iq0=cplex_rhoij*pawrhoij(iatom)%lmn2_size
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         cpj0(1:2,1)=cwaveprj (iatm,1)%cp(1:2,jlmn)
         cpj0(1:2,2)=cwaveprj (iatm,2)%cp(1:2,jlmn)
         cpj1(1:2,1)=cwaveprj1(iatm,1)%cp(1:2,jlmn)
         cpj1(1:2,2)=cwaveprj1(iatm,2)%cp(1:2,jlmn)
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           cpi0(1:2,1)=cwaveprj (iatm,1)%cp(1:2,ilmn)
           cpi0(1:2,2)=cwaveprj (iatm,2)%cp(1:2,ilmn)
           cpi1(1:2,1)=cwaveprj1(iatm,1)%cp(1:2,ilmn)
           cpi1(1:2,2)=cwaveprj1(iatm,2)%cp(1:2,ilmn)
           ro11_re=zero;ro22_re=zero
           ro12_re=zero;ro21_re=zero
           ro12_im=zero;ro21_im=zero
           do iplex=1,cplex
             ro11_re=ro11_re+cpj0(iplex,1)*cpi1(iplex,1)+cpi0(iplex,1)*cpj1(iplex,1)
             ro22_re=ro22_re+cpj0(iplex,2)*cpi1(iplex,2)+cpi0(iplex,2)*cpj1(iplex,2)
           end do
           pawrhoij(iatom)%rhoij_(klmn_re,1)=pawrhoij(iatom)%rhoij_(klmn_re,1)+weight*(ro11_re+ro22_re)
           if (nspden_rhoij>1) then
             do iplex=1,cplex
               ro12_re=ro12_re+cpj0(iplex,1)*cpi1(iplex,2)+cpi0(iplex,2)*cpj1(iplex,1)
               ro21_re=ro21_re+cpj0(iplex,2)*cpi1(iplex,1)+cpi0(iplex,1)*cpj1(iplex,2)
             end do
             pawrhoij(iatom)%rhoij_(klmn_re,4)=pawrhoij(iatom)%rhoij_(klmn_re,4)+weight*(ro11_re-ro22_re)
             pawrhoij(iatom)%rhoij_(klmn_re,2)=pawrhoij(iatom)%rhoij_(klmn_re,2)+weight*(ro12_re+ro21_re)
             if (cplex==2) then
               !Important note: the present implementation follows eq(15) in Hobbs et al, PRB 62, 11556(2000)
               ! rho^alpha,beta_ij = Sum[<Psi^beta|pi><pj|Psi^alpha]  (alpha and beta exponents inverted)
               ro12_im=cpj0(2,1)*cpi1(1,2)-cpi1(2,2)*cpj0(1,1)+cpi0(1,2)*cpj1(2,1)-cpj1(1,1)*cpi0(2,2)
               ro21_im=cpj0(2,2)*cpi1(1,1)-cpi1(2,1)*cpj0(1,2)+cpi0(1,1)*cpj1(2,2)-cpj1(1,2)*cpi0(2,1)
               pawrhoij(iatom)%rhoij_(klmn_re,3)=pawrhoij(iatom)%rhoij_(klmn_re,3)+weight*(ro21_im-ro12_im)
             end if
           end if
           if (compute_impart) then
             klmn_im=klmn_re+iq0
             if (nspden_rhoij>1) pawrhoij(iatom)%rhoij_(klmn_re,3)=pawrhoij(iatom)%rhoij_(klmn_re,3)+weight*(ro12_re-ro21_re)
             if (cplex==2) then
               ro11_im=cpj0(2,1)*cpi1(1,1)-cpi1(2,1)*cpj0(1,1)+cpi0(1,1)*cpj1(2,1)-cpj1(1,1)*cpi0(2,1)
               ro22_im=cpj0(2,2)*cpi1(1,2)-cpi1(2,2)*cpj0(1,2)+cpi0(1,2)*cpj1(2,2)-cpj1(1,2)*cpi0(2,2)
               pawrhoij(iatom)%rhoij_(klmn_im,1)=pawrhoij(iatom)%rhoij_(klmn_im,1)+weight*(ro11_im+ro22_im)
               if (nspden_rhoij>1) then
                 pawrhoij(iatom)%rhoij_(klmn_im,4)=pawrhoij(iatom)%rhoij_(klmn_im,4)+weight*(ro11_im-ro22_im)
                 pawrhoij(iatom)%rhoij_(klmn_im,2)=pawrhoij(iatom)%rhoij_(klmn_re,2)+weight*(ro12_im+ro21_im)
               end if
             end if
           end if
         end do
       end do
     end do
   end if

!  Accumulate (n,k) contribution to rhoij1
!  due to derivative of projectors
   if (ipert/=natom+2) then
     compute_impart=(pawrhoij(1)%cplex_rhoij==2)
     compute_impart_cplex=(compute_impart.and.(cplex==2))
     substract_diagonal=(ipert==natom+3)
     if (nspinor==1) then
       do iatom=1,my_natom
         iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
         iatm=atindx(iatom1)
         if (ipert<=natom.and.iatom/=ipert) cycle
         cplex_rhoij=pawrhoij(iatom)%cplex_rhoij
         do jlmn=1,pawrhoij(iatom)%lmn_size
           j0lmn=jlmn*(jlmn-1)/2
           cpj0 (1:2,1)  =cwaveprj(iatm,1)%cp (1:2  ,jlmn)
           dcpj0(1:2,1,1)=cwaveprj(iatm,1)%dcp(1:2,1,jlmn)
           do ilmn=1,jlmn
             klmn=j0lmn+ilmn
             klmn_re=cplex_rhoij*(klmn-1)+1
             cpi0 (1:2,1)  =cwaveprj(iatm,1)%cp (1:2  ,ilmn)
             dcpi0(1:2,1,1)=cwaveprj(iatm,1)%dcp(1:2,1,ilmn)
             ro11_re=zero
             do iplex=1,cplex
               ro11_re=ro11_re+dcpi0(iplex,1,1)*cpj0(iplex,1)+cpi0(iplex,1)*dcpj0(iplex,1,1)
             end do
             if (substract_diagonal) then
               do iplex=1,cplex
                 ro11_re=ro11_re-cpi0(iplex,1)*cpj0(iplex,1)
               end do
             end if
             pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol)+weight*ro11_re
!            This imaginary part does not have to be computed
!            It is cancelled because rho_ij+rho_ji is stored in rho_ij
             if (compute_impart_cplex) then
               klmn_im=klmn_re+1
               ro11_im=dcpi0(1,1,1)*cpj0(2,1)-dcpi0(2,1,1)*cpj0(1,1)+cpi0(1,1)*dcpj0(2,1,1)-cpi0(2,1)*dcpj0(1,1,1)
               if (substract_diagonal) then
                 ro11_im=ro11_im-cpi0(1,1)*cpj0(2,1)+cpi0(2,1)*cpj0(1,1)
               end if
               pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol)+weight*ro11_im
             end if
           end do
         end do
       end do
     else ! nspinor=2
       do iatom=1,my_natom
         iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
         iatm=atindx(iatom1)
         if (ipert<=natom.and.iatom/=ipert) cycle
         cplex_rhoij=pawrhoij(iatom)%cplex_rhoij
         nspden_rhoij=pawrhoij(iatom)%nspden
         do jlmn=1,pawrhoij(iatom)%lmn_size
           j0lmn=jlmn*(jlmn-1)/2
           cpj0 (1:2,1)  =cwaveprj(iatm,1)%cp (1:2  ,jlmn)
           dcpj0(1:2,1,1)=cwaveprj(iatm,1)%dcp(1:2,1,jlmn)
           cpj0 (1:2,2)  =cwaveprj(iatm,2)%cp (1:2  ,jlmn)
           dcpj0(1:2,2,1)=cwaveprj(iatm,2)%dcp(1:2,1,jlmn)
           do ilmn=1,jlmn
             klmn=j0lmn+ilmn
             klmn_re=cplex_rhoij*(klmn-1)+1
             cpi0 (1:2,1)  =cwaveprj(iatm,1)%cp (1:2  ,ilmn)
             dcpi0(1:2,1,1)=cwaveprj(iatm,1)%dcp(1:2,1,ilmn)
             cpi0 (1:2,2)  =cwaveprj(iatm,2)%cp (1:2  ,ilmn)
             dcpi0(1:2,2,1)=cwaveprj(iatm,2)%dcp(1:2,1,ilmn)
             ro11_re=zero;ro22_re=zero
             ro12_re=zero;ro21_re=zero
             ro12_im=zero;ro21_im=zero
             do iplex=1,cplex
               ro11_re=ro11_re+dcpi0(iplex,1,1)*cpj0(iplex,1)+cpi0(iplex,1)*dcpj0(iplex,1,1)
               ro22_re=ro22_re+dcpi0(iplex,2,1)*cpj0(iplex,2)+cpi0(iplex,2)*dcpj0(iplex,2,1)
             end do
             if (substract_diagonal) then
               do iplex=1,cplex
                 ro11_re=ro11_re-cpi0(iplex,1)*cpj0(iplex,1)
                 ro22_re=ro22_re-cpi0(iplex,2)*cpj0(iplex,2)
               end do
             end if
             pawrhoij(iatom)%rhoij_(klmn_re,1)=pawrhoij(iatom)%rhoij_(klmn_re,1)+weight*(ro11_re+ro22_re)
             if (nspden_rhoij>1) then
               do iplex=1,cplex
                 ro12_re=ro12_re+dcpi0(iplex,2,1)*cpj0(iplex,1)+cpi0(iplex,2)*dcpj0(iplex,1,1)
                 ro21_re=ro21_re+dcpi0(iplex,1,1)*cpj0(iplex,2)+cpi0(iplex,1)*dcpj0(iplex,2,1)
               end do
               if (substract_diagonal) then
                 do iplex=1,cplex
                   ro12_re=ro12_re-cpi0(iplex,2)*cpj0(iplex,1)
                   ro21_re=ro21_re-cpi0(iplex,1)*cpj0(iplex,2)
                 end do
               end if
               pawrhoij(iatom)%rhoij_(klmn_re,4)=pawrhoij(iatom)%rhoij_(klmn_re,4)+weight*(ro11_re-ro22_re)
               pawrhoij(iatom)%rhoij_(klmn_re,2)=pawrhoij(iatom)%rhoij_(klmn_re,2)+weight*(ro12_re+ro21_re)
               if (cplex==2) then
                 !Important note: the present implementation follows eq(15) in Hobbs et al, PRB 62, 11556(2000)
                 ! rho^alpha,beta_ij = Sum[<Psi^beta|pi><pj|Psi^alpha]  (alpha and beta exponents inverted)
                 ro12_im=dcpi0(1,2,1)*cpj0(2,1)-dcpi0(2,2,1)*cpj0(1,1)+cpi0(1,2)*dcpj0(2,1,1)-cpi0(2,2)*dcpj0(1,1,1)
                 ro21_im=dcpi0(1,1,1)*cpj0(2,2)-dcpi0(2,1,1)*cpj0(1,2)+cpi0(1,1)*dcpj0(2,2,1)-cpi0(2,1)*dcpj0(1,2,1)
                 if (substract_diagonal) then
                   ro12_im=ro12_im-cpi0(1,2)*cpj0(2,1)+cpi0(2,2)*cpj0(1,1)
                   ro21_im=ro21_im-cpi0(1,1)*cpj0(2,2)+cpi0(2,1)*cpj0(1,2)
                 end if
                 pawrhoij(iatom)%rhoij_(klmn_re,3)=pawrhoij(iatom)%rhoij_(klmn_re,3)+weight*(ro21_im-ro12_im)
               end if
             end if
             if (compute_impart) then
               klmn_im=klmn_re+1
               if (nspden_rhoij>1) pawrhoij(iatom)%rhoij_(klmn_im,3)=pawrhoij(iatom)%rhoij_(klmn_im,3)+weight*(ro12_re-ro21_re)
               if (cplex==2) then
                 ro11_im=dcpi0(1,1,1)*cpj0(2,1)-dcpi0(2,1,1)*cpj0(1,1)+cpi0(1,1)*dcpj0(2,1,1)-cpi0(2,1)*dcpj0(1,1,1)
                 ro22_im=dcpi0(1,2,1)*cpj0(2,2)-dcpi0(2,2,1)*cpj0(1,2)+cpi0(1,2)*dcpj0(2,2,1)-cpi0(2,2)*dcpj0(1,2,1)
                 if (substract_diagonal) then
                   ro11_im=ro11_im-cpi0(1,1)*cpj0(2,1)+cpi0(2,1)*cpj0(1,1)
                   ro22_im=ro22_im-cpi0(1,2)*cpj0(2,2)+cpi0(2,2)*cpj0(1,2)
                 end if
                 pawrhoij(iatom)%rhoij_(klmn_im,1)=pawrhoij(iatom)%rhoij_(klmn_im,1)+weight*(ro11_im+ro22_im)
                 if (nspden_rhoij>1) then
                   pawrhoij(iatom)%rhoij_(klmn_im,4)=pawrhoij(iatom)%rhoij_(klmn_im,4)+weight*(ro11_im-ro22_im)
                   pawrhoij(iatom)%rhoij_(klmn_im,2)=pawrhoij(iatom)%rhoij_(klmn_im,2)+weight*(ro12_im+ro21_im)
                 end if
               end if
             end if
           end do
         end do
       end do
     end if
   end if

 else if (option==3) then

!  ==================================================================
!  === OPTION 3: Accumulate (n,k) contribution to drhoij/dr =========
!  ==================================================================

   compute_impart=((.not.usetimerev).and.(pawrhoij(1)%cplex_rhoij==2))
   compute_impart_cplex=((compute_impart).and.(cplex==2))
   if (nspinor==1) then
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       cplex_rhoij=pawrhoij(iatom)%cplex_rhoij
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         cpj0(1:cplex,1)         =cwaveprj(iatm,1)%cp (1:cplex,jlmn)
         dcpj0(1:cplex,1,1:ncpgr)=cwaveprj(iatm,1)%dcp(1:cplex,1:ncpgr,jlmn)
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           cpi0(1:cplex,1)         =cwaveprj(iatm,1)%cp (1:cplex,ilmn)
           dcpi0(1:cplex,1,1:ncpgr)=cwaveprj(iatm,1)%dcp(1:cplex,1:ncpgr,ilmn)
           do mu=1,ncpgr
             ro11_re=zero
             do iplex=1,cplex
               ro11_re=ro11_re+dcpi0(iplex,1,mu)*cpj0(iplex,1)+cpi0(iplex,1)*dcpj0(iplex,1,mu)
             end do
             pawrhoij(iatom)%grhoij(mu,klmn_re,isppol)=pawrhoij(iatom)%grhoij(mu,klmn_re,isppol)+weight*ro11_re
           end do
           if (compute_impart_cplex) then
             klmn_im=klmn_re+1
             do mu=1,ncpgr
               ro11_im=dcpi0(1,1,mu)*cpj0(2,1)+cpi0(1,1)*dcpj0(2,1,mu)-dcpi0(2,1,mu)*cpj0(1,1)-cpi0(2,1)*dcpj0(1,1,mu)
               pawrhoij(iatom)%grhoij(mu,klmn_im,isppol)=pawrhoij(iatom)%grhoij(mu,klmn_im,isppol)+weight*ro11_im
             end do
           end if
         end do
       end do
     end do
   else ! nspinor=2
     do iatom=1,my_natom
       iatom1=iatom;if (paral_atom) iatom1=my_atmtab(iatom)
       iatm=atindx(iatom1)
       cplex_rhoij=pawrhoij(iatom)%cplex_rhoij
       nspden_rhoij=pawrhoij(iatom)%nspden
       do jlmn=1,pawrhoij(iatom)%lmn_size
         j0lmn=jlmn*(jlmn-1)/2
         cpj0(1:cplex,1)     =cwaveprj(iatm,1)%cp (1:cplex,jlmn)
         cpj0(1:cplex,2)     =cwaveprj(iatm,2)%cp (1:cplex,jlmn)
         dcpj0(1:cplex,1,1:ncpgr)=cwaveprj(iatm,1)%dcp(1:cplex,1:ncpgr,jlmn)
         dcpj0(1:cplex,2,1:ncpgr)=cwaveprj(iatm,2)%dcp(1:cplex,1:ncpgr,jlmn)
         do ilmn=1,jlmn
           klmn=j0lmn+ilmn
           klmn_re=cplex_rhoij*(klmn-1)+1
           klmn_im=klmn_re+1
           cpi0(1:cplex,1)     =cwaveprj(iatm,1)%cp (1:cplex,ilmn)
           cpi0(1:cplex,2)     =cwaveprj(iatm,2)%cp (1:cplex,ilmn)
           dcpi0(1:cplex,1,1:ncpgr)=cwaveprj(iatm,1)%dcp(1:cplex,1:ncpgr,ilmn)
           dcpi0(1:cplex,2,1:ncpgr)=cwaveprj(iatm,2)%dcp(1:cplex,1:ncpgr,ilmn)
           do mu=1,ncpgr
             ro11_re=zero;ro22_re=zero
             ro12_re=zero;ro21_re=zero
             ro12_im=zero;ro21_im=zero
             do iplex=1,cplex
               ro11_re=ro11_re+dcpi0(iplex,1,mu)*cpj0(iplex,1)+cpi0(iplex,1)*dcpj0(iplex,1,mu)
               ro22_re=ro22_re+dcpi0(iplex,2,mu)*cpj0(iplex,2)+cpi0(iplex,2)*dcpj0(iplex,2,mu)
             end do
             pawrhoij(iatom)%grhoij(mu,klmn_re,1)=pawrhoij(iatom)%grhoij(mu,klmn_re,1)+weight*(ro11_re+ro22_re)
             if (nspden_rhoij>1) then
               do iplex=1,cplex
                 ro12_re=ro12_re+dcpi0(iplex,2,mu)*cpj0(iplex,1)+cpi0(iplex,2)*dcpj0(iplex,1,mu)
                 ro21_re=ro21_re+dcpi0(iplex,1,mu)*cpj0(iplex,2)+cpi0(iplex,1)*dcpj0(iplex,2,mu)
               end do
               pawrhoij(iatom)%grhoij(mu,klmn_re,4)=pawrhoij(iatom)%grhoij(mu,klmn_re,4)+weight*(ro11_re-ro22_re)
               pawrhoij(iatom)%grhoij(mu,klmn_re,2)=pawrhoij(iatom)%grhoij(mu,klmn_re,2)+weight*(ro12_re+ro21_re)
               if (cplex==2) then
                 !Important note: the present implementation follows eq(15) in Hobbs et al, PRB 62, 11556(2000)
                 ! rho^alpha,beta_ij = Sum[<Psi^beta|pi><pj|Psi^alpha]  (alpha and beta exponents inverted)
                 ro12_im=dcpi0(1,2,mu)*cpj0(2,1)+cpi0(1,2)*dcpj0(2,1,mu)-dcpi0(2,2,mu)*cpj0(1,1)-cpi0(2,2)*dcpj0(1,1,mu)
                 ro21_im=dcpi0(1,1,mu)*cpj0(2,2)+cpi0(1,1)*dcpj0(2,2,mu)-dcpi0(2,1,mu)*cpj0(1,2)-cpi0(2,1)*dcpj0(1,2,mu)
                 pawrhoij(iatom)%grhoij(mu,klmn_re,3)=pawrhoij(iatom)%grhoij(mu,klmn_re,3)+weight*(ro21_im-ro12_im)
               end if
             end if
             if (compute_impart) then
               if (nspden_rhoij>1) then
                 pawrhoij(iatom)%grhoij(mu,klmn_im,3)=pawrhoij(iatom)%grhoij(mu,klmn_im,3)+weight*(ro12_re-ro21_re)
               end if
               if (cplex==2) then
                 ro11_im=dcpi0(1,1,mu)*cpj0(2,1)+cpi0(1,1)*dcpj0(2,1,mu)-dcpi0(2,1,mu)*cpj0(1,1)-cpi0(2,1)*dcpj0(1,1,mu)
                 ro22_im=dcpi0(1,2,mu)*cpj0(2,2)+cpi0(1,2)*dcpj0(2,2,mu)-dcpi0(2,2,mu)*cpj0(1,2)-cpi0(2,2)*dcpj0(1,2,mu)
                 pawrhoij(iatom)%grhoij(mu,klmn_im,1)=pawrhoij(iatom)%grhoij(mu,klmn_im,1)+weight*(ro11_im+ro22_im)
                 if (nspden_rhoij>1) then
                   pawrhoij(iatom)%grhoij(mu,klmn_im,4)=pawrhoij(iatom)%grhoij(mu,klmn_im,4)+weight*(ro11_im-ro22_im)
                   pawrhoij(iatom)%grhoij(mu,klmn_im,2)=pawrhoij(iatom)%grhoij(mu,klmn_im,2)+weight*(ro12_im+ro21_im)
                 end if
               end if
             end if
           end do
         end do
       end do
     end do
   end if

!  End
 end if ! option

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine pawaccrhoij
!!***

!----------------------------------------------------------------------

!!****f* m_paw_occupancies/initrhoij
!! NAME
!! initrhoij
!!
!! FUNCTION
!! Initialize PAW rhoij occupancies (in packed storage)
!! from atomic ones
!!
!! INPUTS
!!  cpxocc=1 if rhoij are real, 2 if they are complex
!!  lexexch(ntypat)=l on which local exact-exchange is applied for a given type of atom
!!  lpawu(ntypat)=l on which U is applied for a given type of atom (PAW+U)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of atom types
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated in PAW augmentation regions
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!                                     (containing initial rhoij)
!!  qphase=2 if rhoij have a exp(iqR) phase, 1 if not (typical use: response function at q<>0)
!!  spinat(3,natom)=initial spin of each atom, in unit of hbar/2.
!!  typat(natom)=type of each atom
!!  === Optional arguments
!!    ngrhoij=number of gradients to be allocated (OPTIONAL, default=0)
!!    nlmnmix=number of rhoij elements to be mixed during SCF cycle (OPTIONAL, default=0)
!!    use_rhoij_=1 if pawrhoij(:)%rhoij_ has to be allocated (OPTIONAL, default=0)
!!    use_rhoijres=1 if pawrhoij(:)%rhoijres has to be allocated (OPTIONAL, default=0)

!!
!! OUTPUT
!!  pawrhoij(natom) <type(pawrhoij_type)>=rhoij quantities for each atom
!!                                        in packed storage
!!
!! PARENTS
!!      gstate,respfn,setup_positron
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawrhoij_alloc
!!
!! SOURCE

subroutine initrhoij(cpxocc,lexexch,lpawu,my_natom,natom,nspden,nspinor,nsppol,&
&                    ntypat,pawrhoij,pawspnorb,pawtab,qphase,spinat,typat,&
&                    ngrhoij,nlmnmix,use_rhoij_,use_rhoijres,& ! optional arguments
&                    mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cpxocc,my_natom,natom,nspden,nspinor,nsppol,ntypat,pawspnorb,qphase
 integer,intent(in),optional :: comm_atom,ngrhoij,nlmnmix,use_rhoij_,use_rhoijres
 character(len=500) :: message
!arrays
 integer,intent(in) :: lexexch(ntypat),lpawu(ntypat)
 integer,intent(in) :: typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: spinat(3,natom)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!Arrays
!scalars
 integer :: cplex_rhoij,iatom,iatom_rhoij,ilmn,ispden,itypat,j0lmn,jl,jlmn,jspden
 integer :: klmn,klmn1,ln,lnspinat0,my_comm_atom
 integer :: ngrhoij0,nlmnmix0,nselect,nselect1,nspden_rhoij,qphase_rhoij
 integer :: use_rhoij_0,use_rhoijres0
 real(dp) :: ratio,ro,roshift,zratio,zz
 logical :: my_atmtab_allocated,paral_atom,spinat_zero,test_exexch,test_pawu,test_lnspinat
!arrays
 integer,pointer :: my_atmtab(:),lnspinat(:)
 real(dp),allocatable :: occ(:)
!************************************************************************

 DBG_ENTER("COLL")

!PAW+U and local exact-exchange restriction
 do itypat=1,ntypat
   if (lpawu(itypat)/=lexexch(itypat).and. lpawu(itypat)/=-1.and.lexexch(itypat)/=-1) then
     message = ' lpawu must be equal to lexexch !'
     MSG_ERROR(message)
   end if
 end do

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,qphase_rhoij=qphase_rhoij,nspden_rhoij=nspden_rhoij,&
&                          nspden=nspden,spnorb=pawspnorb,cpxocc=cpxocc,cplex=qphase)

 ratio=one;if (nspden_rhoij==2) ratio=half
 spinat_zero=all(abs(spinat(:,:))<tol10)

 if (my_natom>0) then
   ngrhoij0=0;if (present(ngrhoij)) ngrhoij0=ngrhoij
   nlmnmix0=0;if (present(nlmnmix)) nlmnmix0=nlmnmix
   use_rhoij_0=0;if (present(use_rhoij_)) use_rhoij_0=use_rhoij_
   use_rhoijres0=0;if (present(use_rhoijres)) use_rhoijres0=use_rhoijres
   if (paral_atom) then
     call pawrhoij_alloc(pawrhoij,cplex_rhoij,nspden_rhoij,nspinor,nsppol,typat,&
&     ngrhoij=ngrhoij0,nlmnmix=nlmnmix0,use_rhoij_=use_rhoij_0,use_rhoijres=use_rhoijres0,&
&     qphase=qphase_rhoij,pawtab=pawtab,comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
   else
     call pawrhoij_alloc(pawrhoij,cplex_rhoij,nspden_rhoij,nspinor,nsppol,typat,qphase=qphase_rhoij,&
&     pawtab=pawtab,ngrhoij=ngrhoij0,nlmnmix=nlmnmix0,use_rhoij_=use_rhoij_0,use_rhoijres=use_rhoijres0)
   end if
 end if

 do iatom_rhoij=1,my_natom
   iatom=iatom_rhoij;if (paral_atom) iatom=my_atmtab(iatom_rhoij)
   itypat=typat(iatom)
   nselect=0
   ABI_ALLOCATE(lnspinat,(pawtab(itypat)%basis_size))
   lnspinat=-1
! Determine occupancies of each orbital
   if (nspden_rhoij==2) then
     ABI_ALLOCATE(occ,(pawtab(itypat)%basis_size))
     occ=zero
     do jlmn=1,pawtab(itypat)%lmn_size
       ln=pawtab(itypat)%indlmn(5,jlmn)
       klmn=jlmn*(jlmn+1)/2
       occ(ln)=occ(ln)+pawtab(itypat)%rhoij0(klmn)
     end do
     do ln=1,pawtab(itypat)%basis_size
       if(pawtab(itypat)%orbitals(ln)==0.and.occ(ln)==1) lnspinat(ln)=ln
       if(pawtab(itypat)%orbitals(ln)==1.and.(occ(ln)>=1.and.occ(ln)<=5)) lnspinat(ln)=ln
       if(pawtab(itypat)%orbitals(ln)==2.and.(occ(ln)>=1.and.occ(ln)<=9)) lnspinat(ln)=ln
       if(pawtab(itypat)%orbitals(ln)==3.and.(occ(ln)>=1.and.occ(ln)<=13)) lnspinat(ln)=ln
     end do
     ABI_DEALLOCATE(occ)
   end if
   lnspinat0=maxval(lnspinat)
   lnspinat0=-1

!  Determine Z (trace of rhoij0 or part of it)
   zz=zero
   do jlmn=1,pawtab(itypat)%lmn_size
     jl=pawtab(itypat)%indlmn(1,jlmn)
     ln=pawtab(itypat)%indlmn(5,jlmn)
     j0lmn=jlmn*(jlmn-1)/2
     test_lnspinat=(lnspinat0==-1.or.lnspinat(ln)==ln)
     test_pawu=(lpawu(itypat)==-1.or.lpawu(itypat)==jl)
     test_exexch=(lexexch(itypat)==-1.or.lexexch(itypat)==jl)
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       if ((ilmn==jlmn).and.test_pawu.and.test_exexch.and.test_lnspinat) &
&       zz=zz+pawtab(itypat)%rhoij0(klmn)
     end do
   end do

!  Compute rhoij from tabulated value and magnetization
   do ispden=1,nspden_rhoij

     zratio=zero
     roshift=one
     ratio=one
     if (nspden_rhoij==2) then
       ratio=half
       if ((spinat(3,iatom)>zero.and.ispden==1).or.&
&       (spinat(3,iatom)<zero.and.ispden==2)) then
         if(abs(zz)>tol12)then
           zratio=two*abs(spinat(3,iatom))/zz
         else
           zratio=zero
         end if
       end if
     else if (nspden_rhoij==4.and.ispden>=2) then
       roshift=zero
       if(abs(zz)>tol12)then
         zratio=spinat(ispden-1,iatom)/zz
       else
         zratio=zero
       end if
     end if

     nselect=0;nselect1=1-cpxocc
     do jlmn=1,pawtab(itypat)%lmn_size
       jl=pawtab(itypat)%indlmn(1,jlmn)
       ln=pawtab(itypat)%indlmn(5,jlmn)
       j0lmn=jlmn*(jlmn-1)/2
       test_lnspinat=(lnspinat0==-1.or.lnspinat(ln)==ln)
       test_pawu=(lpawu(itypat)==-1.or.lpawu(itypat)==jl)
       test_exexch=(lexexch(itypat)==-1.or.lexexch(itypat)==jl)
       do ilmn=1,jlmn
         klmn=j0lmn+ilmn
         ro=pawtab(itypat)%rhoij0(klmn)
         if ((ilmn==jlmn).and.test_pawu.and.test_exexch.and.test_lnspinat) then
           ro=ro*ratio*(roshift+zratio)
         else
           ro=ro*ratio*roshift
         end if

         klmn1=cpxocc*(klmn-1)+1
         if (abs(ro)>tol10) then
           pawrhoij(iatom_rhoij)%rhoijp(klmn1,ispden)=ro
         else
           pawrhoij(iatom_rhoij)%rhoijp(klmn1,ispden)=zero
         end if

         if (ispden==nspden_rhoij) then
           if (any(abs(pawrhoij(iatom_rhoij)%rhoijp(klmn1,:))>tol10)) then
             nselect=nselect+1;nselect1=nselect1+cpxocc
             pawrhoij(iatom_rhoij)%rhoijselect(nselect)=klmn
             do jspden=1,nspden_rhoij
               pawrhoij(iatom_rhoij)%rhoijp(nselect1,jspden)=pawrhoij(iatom_rhoij)%rhoijp(klmn1,jspden)
             end do
           end if
         end if

       end do
     end do

   end do
   pawrhoij(iatom_rhoij)%nrhoijsel=nselect
   if (nselect<pawrhoij(iatom_rhoij)%lmn2_size) &
&    pawrhoij(iatom_rhoij)%rhoijselect(nselect+1:pawrhoij(iatom_rhoij)%lmn2_size)=0

!  Non-collinear magnetism: avoid zero magnetization, because it produces numerical instabilities
!    Add a small real to the magnetization ; not yet activated => must be tested.
!   if (pawrhoij(iatom_rhoij)%nspden==4.and.spinat_zero) then
!     pawrhoij(iatom_rhoij)%rhoijp(:,4)=pawrhoij(iatom_rhoij)%rhoijp(:,4)+tol10
!   end if
   ABI_DEALLOCATE(lnspinat)
 end do ! iatom_rhoij

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine initrhoij
!!***

!----------------------------------------------------------------------

END MODULE m_paw_occupancies
!!***
