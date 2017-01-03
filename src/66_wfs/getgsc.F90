!{\src2tex{textfont=tt}}
!!****f* ABINIT/getgsc
!!
!! NAME
!! getgsc
!!
!! FUNCTION
!! Compute <G|S|C> for all input vectors |Cnk> at a given k-point,
!!              OR for one input vector |Cnk>.
!! |Cnk> are expressed in reciprocal space.
!! S is the overlap operator between |Cnk> (used for PAW).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mcg)=planewave coefficients of wavefunctions
!!  cprj(natom,mcprj)= wave functions projected with non-local projectors: cprj=<p_i|Cnk>
!!  gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  ibg=shift to be applied on the location of data in the array cprj (beginning of current k-point)
!!  icg=shift to be applied on the location of data in the array cg (beginning of current k-point)
!!  igsc=shift to be applied on the location of data in the array gsc (beginning of current k-point)
!!  ikpt,isppol=indexes of current (spin.kpoint)
!!  mcg=second dimension of the cg array
!!  mcprj=second dimension of the cprj array
!!  mgsc=second dimension of the gsc array
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in unit cell.
!!  nband= if positive: number of bands at this k point for that spin polarization
!!         if negative: abs(nband) is the index of the only band to be computed
!!  npw_k=number of planewaves in basis for given k point.
!!  nspinor=number of spinorial components of the wavefunctions
!! [select_k]=optional, option governing the choice of k points to be used.
!!             gs_ham datastructure contains quantities needed to apply overlap operator
!!             in reciprocal space between 2 kpoints, k and k^prime (equal in most cases);
!!             if select_k=1, <k^prime|S|k>       is applied [default]
!!             if select_k=2, <k|S|k^prime>       is applied
!!             if select_k=3, <k|S|k>             is applied
!!             if select_k=4, <k^prime|S|k^prime> is applied
!!
!! OUTPUT
!!  gsc(2,mgsc)= <g|S|Cnk> or <g|S^(1)|Cnk> (S=overlap)
!!
!! PARENTS
!!      dfpt_vtowfk
!!
!! CHILDREN
!!      nonlop,pawcprj_alloc,pawcprj_copy,pawcprj_free,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getgsc(cg,cprj,gs_ham,gsc,ibg,icg,igsc,ikpt,isppol,&
&                 mcg,mcprj,mgsc,mpi_enreg,natom,nband,npw_k,nspinor,select_k)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors

 use m_hamiltonian, only : gs_hamiltonian_type
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_copy
 use m_nonlop

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getgsc'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: ibg,icg,igsc,ikpt,isppol,mcg,mcprj
 integer,intent(in) :: mgsc,natom,nband,npw_k,nspinor
 integer,intent(in),optional :: select_k
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout),target :: gs_ham
!arrays
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(out) :: gsc(2,mgsc)
 type(pawcprj_type),intent(in) :: cprj(natom,mcprj)

!Local variables-------------------------------
!scalars
 integer :: choice,cpopt,dimenl1,dimenl2,iband,iband1,iband2,ierr,index_cg,index_cprj
 integer :: index_gsc,me,my_nspinor,paw_opt,select_k_,signs,tim_nonlop,useylm
 character(len=500) :: msg
!arrays
 real(dp) :: enlout_dum(1),tsec(2)
 real(dp),allocatable :: cwavef(:,:),scwavef(:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)
 if(gs_ham%usepaw==0) then
   msg='Only compatible with PAW (usepaw=1) !'
   MSG_BUG(msg)
 end if
 if(nband<0.and.(mcg<npw_k*my_nspinor.or.mgsc<npw_k*my_nspinor.or.mcprj<my_nspinor)) then
   msg='Invalid value for mcg, mgsc or mcprj !'
   MSG_BUG(msg)
 end if

!Keep track of total time spent in getgsc:
 call timab(565,1,tsec)

 gsc = zero

!Prepare some data
 ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
 ABI_ALLOCATE(scwavef,(2,npw_k*my_nspinor))
 if (mcprj>0) then
   ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,my_nspinor))
   call pawcprj_alloc(cwaveprj,0,gs_ham%dimcprj)
 end if
 dimenl1=gs_ham%dimekb1;dimenl2=natom;tim_nonlop=0
 choice=1;signs=2;cpopt=-1+3*gs_ham%usecprj;paw_opt=3;useylm=1
 select_k_=1;if (present(select_k)) select_k_=select_k
 me=mpi_enreg%me_kpt

!Loop over bands
 index_cprj=ibg;index_cg=icg;index_gsc=igsc
 if (nband>0) then
   iband1=1;iband2=nband
 else if (nband<0) then
   iband1=abs(nband);iband2=iband1
   index_cprj=index_cprj+(iband1-1)*my_nspinor
   index_cg  =index_cg  +(iband1-1)*npw_k*my_nspinor
   index_gsc =index_gsc +(iband1-1)*npw_k*my_nspinor
 end if

 do iband=iband1,iband2

   if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me.and.nband>0) then
     gsc(:,1+index_gsc:npw_k*my_nspinor+index_gsc)=zero
     index_cprj=index_cprj+my_nspinor
     index_cg=index_cg+npw_k*my_nspinor
     index_gsc=index_gsc+npw_k*my_nspinor
     cycle
   end if

!  Retrieve WF at (n,k)
   cwavef(:,1:npw_k*my_nspinor)=cg(:,1+index_cg:npw_k*my_nspinor+index_cg)
   if (gs_ham%usecprj==1) then
     call pawcprj_copy(cprj(:,1+index_cprj:my_nspinor+index_cprj),cwaveprj)
   end if

!  Compute <g|S|Cnk>
   call nonlop(choice,cpopt,cwaveprj,enlout_dum,gs_ham,0,(/zero/),mpi_enreg,1,1,paw_opt,&
&   signs,scwavef,tim_nonlop,cwavef,cwavef,select_k=select_k_)

   gsc(:,1+index_gsc:npw_k*my_nspinor+index_gsc)=scwavef(:,1:npw_k*my_nspinor)

!  End of loop over bands
   index_cprj=index_cprj+my_nspinor
   index_cg=index_cg+npw_k*my_nspinor
   index_gsc=index_gsc+npw_k*my_nspinor
 end do

!Reduction in case of parallelization
 if ((xmpi_paral==1)) then
   call timab(48,1,tsec)
   call xmpi_sum(gsc,mpi_enreg%comm_band,ierr)
   call timab(48,2,tsec)
 end if

!Memory deallocation
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(scwavef)
 if (gs_ham%usecprj==1) then
   call pawcprj_free(cwaveprj)
   ABI_DATATYPE_DEALLOCATE(cwaveprj)
 end if

 call timab(565,2,tsec)

 DBG_EXIT("COLL")

end subroutine getgsc
!!***
