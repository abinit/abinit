!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawmkrhoij
!!
!! NAME
!! pawmkrhoij
!!
!! FUNCTION
!! Calculate the PAW quantities rhoij (augmentation occupancies)
!! Remember:for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!!  pawprtvol=control print volume and debugging output for PAW
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
!!      vtorho
!!
!! CHILDREN
!!      pawaccrhoij,pawcprj_alloc,pawcprj_free,pawcprj_gather_spin,pawcprj_get
!!      pawio_print_ij,pawrhoij_free,pawrhoij_init_unpacked
!!      pawrhoij_mpisum_unpacked,wrtout
!!
!! NOTES
!!  The cprj are distributed over band processors.
!!  Only the mod((iband-1)/mpi_enreg%bandpp,mpi_enreg%nproc_band) projectors
!!  are stored on each proc.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawmkrhoij(atindx,atindx1,cprj,dimcprj,istwfk,kptopt,mband,mband_cprj,mcprj,mkmem,mpi_enreg,&
&                      natom,nband,nkpt,nspinor,nsppol,occ,paral_kgb,paw_dmft,&
&                      pawprtvol,pawrhoij,unpaw,usewvl,wtk)


 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_pawrhoij, only : pawrhoij_type, pawrhoij_init_unpacked, pawrhoij_mpisum_unpacked, pawrhoij_free
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_get, &
&                       pawcprj_gather_spin, pawcprj_free
 use m_paw_io,   only : pawio_print_ij
 use m_paw_dmft, only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawmkrhoij'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_65_paw, except_this_one => pawmkrhoij
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: kptopt,mband,mband_cprj,mcprj,mkmem,natom,nkpt,nspinor,nsppol
 integer,intent(in) :: paral_kgb,pawprtvol,unpaw,usewvl
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
 integer :: iatom,iatom_tot,ib,ib1,iband,iband1,ibc1,ibg,ib_this_proc,ierr
 integer :: ikpt,iorder_cprj,isppol,jb_this_proc,jbg,me,my_nspinor,nband_k,nband_k_cprj
 integer :: nbandc1,nband_k_cprj_read,nband_k_cprj_used,nprocband,nrhoij,nsp2
 integer :: option,spaceComm,use_nondiag_occup_dmft
 logical :: locc_test,paral_atom,usetimerev
 real(dp) :: wtk_k
 character(len=4) :: wrt_mode
 character(len=500) :: msg

!arrays
 integer :: idum(0)
 real(dp) :: occup(2)
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
   msg=' mband/mband_cprj must be equal to mband !'
   MSG_BUG(msg)
 end if
 if (paw_dmft%use_sc_dmft/=0.and.nprocband/=1) then
   write(msg,'(4a,e14.3,a)') ch10,&
&   ' Parallelization over bands is not yet compatible with self-consistency in DMFT !',ch10,&
&   ' Calculation is thus restricted to nstep =1.'
   MSG_WARNING(msg)
 end if

 if( usewvl==1 .and. (nprocband/=1)) then
   write(msg,'(2a)') ch10,&
&   '  ERROR: parallelization over bands is not compatible with WAVELETS'
   MSG_ERROR(msg)
 end if

!Initialise and check dmft variables
 if(paw_dmft%use_sc_dmft/=0) then
   nbandc1=(paw_dmft%mbandc-1)*paw_dmft%use_sc_dmft+1
 else
   nbandc1=1
 end if

!Size of pawrhoij datastructure
 nrhoij=size(pawrhoij)

!Check if pawrhoij is distributed over atomic sites
 paral_atom=(nrhoij/=natom.and.mpi_enreg%nproc_atom>1)
 if (paral_atom.and.nrhoij/=mpi_enreg%my_natom) then
   msg=' Size of pawrhoij should be natom or my_natom !'
   MSG_BUG(msg)
 end if

!Allocate temporary cwaveprj storage
 ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,nspinor))
 call pawcprj_alloc(cwaveprj,0,dimcprj)
 if(paw_dmft%use_sc_dmft/=0) then
   ABI_DATATYPE_ALLOCATE(cwaveprjb,(natom,nspinor))
   call pawcprj_alloc(cwaveprjb,0,dimcprj)
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

!        DMFT stuff: extract cprj and occupations for additional band
         if(paw_dmft%use_sc_dmft /= 0) then
           ib1 = paw_dmft%include_bands(ibc1)
!          write(std_out,*) 'use_sc_dmft=1 ib,ib1',ib,ib1
           iband1 = bdtot_index+ib1
!          write(std_out,*) 'ib, ib1          ',paw_dmft%band_in(ib),paw_dmft%band_in(ib1)
           if(paw_dmft%band_in(ib)) then
             if(.not.paw_dmft%band_in(ib1))  stop
             use_nondiag_occup_dmft = 1
             occup(1) = paw_dmft%occnd(1,ib,ib1,ikpt,isppol)
             if(nspinor==2) occup(2) = paw_dmft%occnd(2,ib,ib1,ikpt,isppol)
             if(nspinor==1) occup(2) = zero
             locc_test = abs(paw_dmft%occnd(1,ib,ib1,ikpt,isppol))+abs(paw_dmft%occnd(2,ib,ib1,ikpt,isppol))>tol8
!            write(std_out,*) 'use_sc_dmft=1,band_in(ib)=1, ib,ibc1',ib,ib1,locc_test
             if (locc_test .or. mkmem == 0) then
               call pawcprj_get(atindx1,cwaveprjb,cprj_ptr,natom,ib1,ibg,ikpt,iorder_cprj,isppol,&
&               mband_cprj,mkmem,natom,1,nband_k_cprj,nspinor,nsppol,unpaw,&
&               mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
             end if
           else
             use_nondiag_occup_dmft = 0
             locc_test = (abs(occ(iband))>tol8)
             occup(1) = occ(iband)
             if(ibc1 /= 1 .and. .not.(paw_dmft%band_in(ib))) cycle
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
&           mband_cprj,mkmem,natom,1,nband_k_cprj,nspinor,nsppol,unpaw,&
&           mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         end if

!        Accumulate contribution from (occupied) current band
         if (locc_test) then
           if(use_nondiag_occup_dmft == 1) then
             call pawaccrhoij(atindx,cplex,cwaveprj,cwaveprjb,0,isppol,nrhoij,natom,&
&             nspinor,occup(1),option,pawrhoij_all,usetimerev,wtk_k,&
&             occ_k_2=occup(2))
           else
             call pawaccrhoij(atindx,cplex,cwaveprj,cwaveprj ,0,isppol,nrhoij,natom,&
&             nspinor,occup(1),option,pawrhoij_all,usetimerev,wtk_k)
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

!deallocate temporary cwaveprj/cprj storage
 call pawcprj_free(cwaveprj)
 ABI_DATATYPE_DEALLOCATE(cwaveprj)
 if(paw_dmft%use_sc_dmft/=0) then
   call pawcprj_free(cwaveprjb)
   ABI_DATATYPE_DEALLOCATE(cwaveprjb)
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

!Print info
 if (abs(pawprtvol)>=1) then
   wrt_mode='COLL';if (paral_atom) wrt_mode='PERS'
   do iatom=1,nrhoij
     iatom_tot=iatom;if (paral_atom) iatom_tot=mpi_enreg%my_atmtab(iatom)
     if (pawprtvol>=0.and.iatom_tot/=1.and.iatom_tot/=natom) cycle
     nsp2=pawrhoij(iatom)%nsppol;if (pawrhoij(iatom)%nspden==4) nsp2=4
     write(msg, '(4a,i3,a)') ch10," PAW TEST:",ch10,&
&     ' ====== Values of RHOIJ in pawmkrhoij (iatom=',iatom_tot,') ======'
     if (pawrhoij(iatom)%nspden==2.and.pawrhoij(iatom)%nsppol==1) write(msg,'(3a)') trim(msg),ch10,&
&     '      (antiferromagnetism case: only one spin component)'
     call wrtout(std_out,msg,wrt_mode)
     do isppol=1,nsp2
       if (pawrhoij(iatom)%nspden/=1) then
         write(msg,'(3a)') '   Component ',trim(dspin(isppol+2*(pawrhoij(iatom)%nspden/4))),':'
         call wrtout(std_out,msg,wrt_mode)
       end if
       option=2;if (pawrhoij(iatom)%cplex==2.and.pawrhoij(iatom)%nspinor==1) option=1
       call pawio_print_ij(std_out,pawrhoij(iatom)%rhoij_(:,isppol),pawrhoij(iatom)%lmn2_size,&
&       pawrhoij(iatom)%cplex,pawrhoij(iatom)%lmn_size,-1,idum,0,pawprtvol,idum,&
&       -1._dp,1,opt_sym=option,mode_paral=wrt_mode)
     end do
   end do
 end if

 DBG_EXIT("COLL")

end subroutine pawmkrhoij
!!***
