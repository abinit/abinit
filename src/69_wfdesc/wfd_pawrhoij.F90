!!****f* ABINIT/wfd_pawrhoij
!!
!! NAME
!! wfd_pawrhoij
!!
!! FUNCTION
!! Calculate the PAW quantities rhoij (augmentation occupancies)
!! Remember:for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= wave functions projected with non-local projectors:
!!                                   cprj_nk(i)=<p_i|Cnk> where p_i is a non-local projector.
!!  istwfk(nkpt)=parameter that describes the storage of wfs
!!  kptopt=option for the generation of k points
!!  mband=maximum number of bands
!!  natom=number of atoms in cell
!!  nkpt=number of k points
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol)=occupation number for each band for each k
!!  pawprtvol=control print volume and debugging output for PAW
!!
!! SIDE EFFECTS
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  On input: arrays dimensions
!!  On output:
!!    pawrhoij(:)%rhoij_(lmn2_size,nspden)=
!!          Sum_{n,k} {occ(n,k)*conjugate[cprj_nk(ii)].cprj_nk(jj)} (non symetrized)
!!
!! PARENTS
!!      paw_qpscgw
!!
!! CHILDREN
!!      pawaccrhoij,pawcprj_alloc,pawcprj_free,pawio_print_ij
!!      pawrhoij_mpisum_unpacked,wfd_bks_distrb,wfd_get_cprj,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wfd_pawrhoij(Wfd,Cryst,Bst,kptopt,pawrhoij,pawprtvol)

 use defs_basis
 use defs_datatypes
 use m_profiling_abi
 use m_errors

 use m_crystal,  only : crystal_t
 use m_wfd,      only : wfd_t, wfd_get_cprj, wfd_bks_distrb
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_mpisum_unpacked
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_paw_io,   only : pawio_print_ij

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_pawrhoij'
 use interfaces_14_hidewrite
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: kptopt,pawprtvol
 type(crystal_t),intent(in) :: Cryst
 type(wfd_t),intent(inout) :: Wfd
 type(ebands_t),intent(in) :: Bst
!arrays
 type(pawrhoij_type),intent(inout) :: pawrhoij(Wfd%natom)

!Local variables ---------------------------------------
!scalars
 integer :: cplex,iatom,band,ik_ibz
 integer :: spin,natinc,nband_k,nsp2,option,rhoij_cplex,lmn2_size,nspden
 logical :: usetimerev
 real(dp) :: occup,wtk_k
 character(len=500) :: msg
!arrays
 integer,allocatable :: idum(:)
 !real(dp) :: tsec(2)
 character(len=8),parameter :: dspin(6)=(/"up      ","down    ","dens (n)","magn (x)","magn (y)","magn (z)"/)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)
 integer :: bks_distrb(Wfd%mband,Wfd%nkibz,Wfd%nsppol)
 integer :: got(Wfd%nproc)
 logical :: bks_mask(Wfd%mband,Wfd%nkibz,Wfd%nsppol)

!************************************************************************

 DBG_ENTER("COLL")

 ! Allocate temporary cwaveprj storage (sorted by atom type)
 ABI_DATATYPE_ALLOCATE(cwaveprj,(Wfd%natom,Wfd%nspinor))
 call pawcprj_alloc(cwaveprj,0,Wfd%nlmn_sort)

 ! Initialize output quantities if not already done.
 do iatom=1,Wfd%natom
   if (pawrhoij(iatom)%use_rhoij_==0) then
     rhoij_cplex     = pawrhoij(iatom)%cplex
     lmn2_size = pawrhoij(iatom)%lmn2_size
     nspden    = pawrhoij(iatom)%nspden
     ABI_ALLOCATE(pawrhoij(iatom)%rhoij_,(rhoij_cplex*lmn2_size,nspden))
     pawrhoij(iatom)%use_rhoij_=1
   end if
   pawrhoij(iatom)%rhoij_=zero
 end do

 option=1; usetimerev=(kptopt>0.and.kptopt<3)

 ! Distribute (b,k,s).
 where (ABS(Bst%occ)>tol8)
   bks_mask=.TRUE.
 else where
   bks_mask=.FALSE.
 end where
 got=zero

 call wfd_bks_distrb(Wfd,bks_distrb,got,bks_mask)

 do spin=1,Wfd%nsppol
   do ik_ibz=1,Wfd%nkibz

     nband_k=Wfd%nband(ik_ibz,spin)
     wtk_k=Bst%wtk(ik_ibz)

     cplex=2; if (Wfd%istwfk(ik_ibz)>1) cplex=1

     do band=1,nband_k

       if (bks_distrb(band,ik_ibz,spin) == Wfd%my_rank) then
         !locc_test = (abs(Bst%occ(band,ik_ibz,spin))>tol8)
         occup = Bst%occ(band,ik_ibz,spin)

          ! Extract cprj for current band cwaveprj are sorted by atom type.
          call wfd_get_cprj(Wfd,band,ik_ibz,spin,Cryst,cwaveprj,sorted=.TRUE.)

          ! Accumulate contribution from (occupied) current band
          !if (locc_test) then
           call pawaccrhoij(Cryst%atindx,cplex,cwaveprj,cwaveprj ,0,spin,Wfd%natom,Wfd%natom,&
&            Wfd%nspinor,occup,option,pawrhoij,usetimerev,wtk_k)
          !end if
       end if
     end do !band

   end do !ik_ibz
 end do !spin
 !
 ! Free temporary cwaveprj storage.
 call pawcprj_free(cwaveprj)
 ABI_DATATYPE_DEALLOCATE(cwaveprj)
 !
 !==========================================
 ! MPI: need to exchange arrays between procs
 ! TODO it should be tested.
 call pawrhoij_mpisum_unpacked(pawrhoij,Wfd%comm)
 !
 ! Print info.
 if (abs(pawprtvol)>=1) then
   natinc=1; if(Wfd%natom>1.and.pawprtvol>=0) natinc=Wfd%natom-1
   do iatom=1,Cryst%natom,natinc
     nsp2=pawrhoij(iatom)%nsppol;if (pawrhoij(iatom)%nspden==4) nsp2=4
     write(msg, '(4a,i3,a)') ch10," PAW TEST:",ch10,&
&     ' ====== Values of RHOIJ in wfd_pawrhoij (iatom=',iatom,') ======'
     if (pawrhoij(iatom)%nspden==2.and.pawrhoij(iatom)%nsppol==1) write(msg,'(3a)') trim(msg),ch10,&
&     '      (antiferromagnetism case: only one spin component)'
     call wrtout(std_out,msg,'COLL')
     do spin=1,nsp2
       if (pawrhoij(iatom)%nspden/=1) then
         write(msg, '(3a)') '   Component ',trim(dspin(spin+2*(pawrhoij(iatom)%nspden/4))),':'
         call wrtout(std_out,msg,'COLL')
       end if
       call pawio_print_ij(std_out,pawrhoij(iatom)%rhoij_(:,spin),pawrhoij(iatom)%lmn2_size,&
&       pawrhoij(iatom)%cplex,pawrhoij(iatom)%lmn_size,-1,idum,0,pawprtvol,idum,-1.d0,1)
     end do
   end do
 end if

 DBG_EXIT("COLL")

end subroutine wfd_pawrhoij
!!***
