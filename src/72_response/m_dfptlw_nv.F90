!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dfptlw_nv
!! NAME
!!  m_dfptlw_nv
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2022 ABINIT group (MR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

module m_dfptlw_nv
    
 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_mpinfo
 use m_dtset

 use m_dfpt_elt,    only : dfpt_ewalddq, dfpt_ewalddqdq

 implicit none

 private
!!***

 public :: dfptlw_nv
!!***

! *************************************************************************

contains 
!!***

!!****f* ABINIT/m_dfptlw_nv/dfptlw_nv
!! NAME
!!  dfptlw_nv
!!
!! FUNCTION
!!  This routine calculates the nonvariational contributions to the 
!!  spatial-dispersion third-order energy derivatives.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2 
!!  mpert=maximum number of ipert
!!  my_natom=number of atoms treated by current processor
!!  rmet(3,3)=metric tensor in real space (length units squared)
!!  rfpert(3,mpert,3,mpert,3,mpert) = array defining the type of perturbations
!!       that have to be computed
!!       1   ->   element has to be computed explicitely
!!      -1   ->   use symmetry operations to obtain the corresponding element
!!  ucvol=unit cell volume in (whatever length scale units)**3
!!  xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!!  zion(ntypat)=charge on each type of atom (real number)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!
!! OUTPUT
!!  d3etot_nv(2,3,mpert,3,mpert,3,mpert)= array with the nonvariational
!!              contributions of d3etot
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfptlw_nv(d3etot_nv,dtset,gmet,mpert,my_natom,rfpert,rmet,ucvol,xred,zion, &
&                 mpi_atmtab,comm_atom ) ! optional arguments (parallelism))
    
 implicit none

!Arguments ------------------------------------
!scalars
 integer , intent(in)  :: mpert,my_natom
 real(dp) :: ucvol
 type(dataset_type),intent(in) :: dtset
 integer,optional,intent(in) :: comm_atom

!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 integer,intent(in) :: rfpert(3,mpert,3,mpert,3,mpert)
 real(dp), intent(out) :: d3etot_nv(2,3,mpert,3,mpert,3,mpert)
 real(dp), intent(in) :: gmet(3,3),rmet(3,3),xred(3,dtset%natom),zion(*)

!Local variables-------------------------------
!scalars
 integer :: i1dir,i2dir,i3dir,ii,i1pert,i2pert,i3pert,natom,sumg0 
 real(dp) :: tmpim,tmpre
 character(len=500) :: msg

!arrays
 real(dp),allocatable :: dyewdq(:,:,:,:,:,:)
 real(dp) :: qphon(3)
 
! *************************************************************************

 DBG_ENTER("COLL")

!Anounce start of calculation
 write(msg, '(a,a)' ) ' -- Compute nonvariational contributions -- ',ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

!Initialiations
 natom=dtset%natom
 d3etot_nv(:,:,:,:,:,:,:)=zero
 
 if (dtset%lw_flexo==1.or.dtset%lw_flexo==3) then

   !1st q-gradient of Ewald contribution to the IFCs
   ABI_MALLOC(dyewdq,(2,3,natom,3,natom,3))
   sumg0=0;qphon(:)=zero
   call dfpt_ewalddq(dyewdq,gmet,my_natom,natom,qphon,rmet,sumg0,dtset%typat,ucvol,xred,zion,&
& mpi_atmtab=mpi_atmtab,comm_atom=comm_atom)

   i3pert=natom+8
   do i1pert=1,natom
     do i1dir=1,3
       do i2pert=1,natom
         do i2dir=1,3
           do i3dir=1,3
             tmpre=dyewdq(1,i1dir,i1pert,i2dir,i2pert,i3dir)
             tmpim=dyewdq(2,i1dir,i1pert,i2dir,i2pert,i3dir)
             if (abs(tmpre)>=tol8) d3etot_nv(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)= tmpre
             if (abs(tmpim)>=tol8) d3etot_nv(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)= tmpim
           end do
         end do
       end do
     end do
   end do 

 end if

 !Print results
 if (dtset%prtvol>=10) then
   write(msg,'(3a)') ch10,'LONGWAVE NONVARIATIONAL D3ETOT: ',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   do i1pert=1,mpert
     do i1dir=1,3
       do i2pert=1,mpert
         do i2dir=1,3
           do i3pert=1,mpert
             do i3dir=1,3
               if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then 
                 write(msg,'(3(a,i2,a,i1),2f18.8)') &
         ' perts : ',i1pert,'.',i1dir,' / ',i2pert,'.',i2dir,' / ',i3pert,'.',i3dir,&
                 d3etot_nv(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert), &
                 d3etot_nv(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
                 call wrtout(std_out,msg,'COLL')
                 call wrtout(ab_out,msg,'COLL')
               end if
             end do
           end do
         end do
       end do
     end do
   end do 
   write(msg,'(a)') ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   
 end if

 DBG_EXIT("COLL")

end subroutine dfptlw_nv
!!***

end module m_dfptlw_nv
!!***
