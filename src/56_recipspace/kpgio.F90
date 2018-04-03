!{\src2tex{textfont=tt}}
!!****f* ABINIT/kpgio
!! NAME
!! kpgio
!!
!! FUNCTION
!! Do initialization of kg information.
!! Includes opening disk files for kpgsph i/o.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, AR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ecut=kinetic energy planewave cutoff (hartree)
!!  exchn2n3d=if 1, n2 and n3 are exchanged
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kptns(3,nkpt)=reduced coords of k points
!!  mkmem =number of k points treated by this node.
!!  character(len=4) : mode_paral=either 'COLL' or 'PERS', tells whether
!!   the loop over k points must be done by all processors or not,
!!   in case of parallel execution.
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum number of planewaves as dimensioned in calling routine
!!  nband(nkpt*nsppol)=number of bands at each k point
!!  nkpt=number of k points
!!  nsppol=1 for unpolarized, 2 for polarized
!!
!! OUTPUT
!!  npwarr(nkpt)=array holding npw for each k point, taking into account
!!   the effect of istwfk, and the spreading over processors
!!  npwtot(nkpt)=array holding the total number of plane waves for each k point,
!!  kg(3,mpw*mkmem)=dimensionless coords of G vecs in basis sphere at k point
!!
!! NOTES
!! Note that in case of band parallelism, the number of spin-up
!! and spin-down bands must be equal at each k points
!!
!! PARENTS
!!      dfpt_looppert,dfptff_initberry,gstate,initmv,inwffil,m_cut3d,nonlinear
!!      respfn,scfcv
!!
!! CHILDREN
!!      kpgsph,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine kpgio(ecut,exchn2n3d,gmet,istwfk,kg,kptns,mkmem,nband,nkpt,&
& mode_paral,mpi_enreg,mpw,npwarr,npwtot,nsppol)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_fftcore, only : kpgsph

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kpgio'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: exchn2n3d,mkmem,mpw,nkpt,nsppol
 real(dp),intent(in) :: ecut
 character(len=4),intent(in) :: mode_paral
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol)
 integer,intent(out) :: kg(3,mpw*mkmem),npwarr(nkpt),npwtot(nkpt)
 real(dp),intent(in) :: gmet(3,3),kptns(3,nkpt)

!Local variables-------------------------------
!scalars
 integer :: ierr,ikg,ikpt,istwf_k,me,nband_down,nband_k,npw1
 logical :: test_npw
 character(len=500) :: message
!arrays
 real(dp) :: kpoint(3)

! *************************************************************************

!DEBUG
!write(std_out,*)' kpgio : enter '
!ENDDEBUG

!Define me
 me=mpi_enreg%me_kpt

 if((mpi_enreg%paralbd==1) .and. (mode_paral=='PERS')) then
   if(nsppol==2)then
     do ikpt=1,nkpt
       nband_k=nband(ikpt)
       nband_down=nband(ikpt+nkpt)
       if(nband_k/=nband_down)then
         write(message,'(a,a,a,a,a,a,a,a,i4,a,i4,a,a,a,i4,a,a,a)')ch10,&
&         ' kpgio : ERROR -',ch10,&
&         '  Band parallel case, one must have same number',ch10,&
&         '  of spin up and spin down bands, but input is :',ch10,&
&         '  nband(up)=',nband_k,', nband(down)=',nband_down,',',ch10,&
&         '  for ikpt=',ikpt,'.',ch10,&
&         '  Action : correct nband in your input file.'
!        MG: Tests v3(10,11,17) and v6(67) fail if this test is enabled
!        call wrtout(std_out,message,mode_paral)
       end if
     end do
   end if
 end if

 npwarr(:)=0
 npwtot(:)=0

 kg=0
 ikg=0
!Find (k+G) sphere for each k.

 do ikpt=1,nkpt

   nband_k = nband(ikpt)

   if(mode_paral=='PERS')then
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,-1,me)) cycle
   end if

   kpoint(:)=kptns(:,ikpt)
   istwf_k=istwfk(ikpt)
   call kpgsph(ecut,exchn2n3d,gmet,ikg,ikpt,istwf_k,kg,kpoint,mkmem,mpi_enreg,mpw,npw1)

   test_npw=.true.
   if (xmpi_paral==1)then
     if (mode_paral=='PERS')then
       test_npw=(minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,1:nsppol))==me)
     end if
   end if
   if (test_npw) npwarr(ikpt)=npw1

!  Make sure npw < nband never happens:
!  if (npw1<nband(ikpt)) then
!  write(message, '(a,a,a,a,i5,a,3f8.4,a,a,i10,a,i10,a,a,a,a)' )ch10,&
!  &   ' kpgio : ERROR -',ch10,&
!  &   '  At k point number',ikpt,' k=',(kptns(mu,ikpt),mu=1,3),ch10,&
!  &   '  npw=',npw1,' < nband=',nband(ikpt),ch10,&
!  &   '  Indicates not enough planewaves for desired number of bands.',ch10,&
!  &   '  Action : change either ecut or nband in input file.'
!  MSG_ERROR(message)
!  end if

!  Find boundary of G sphere for efficient zero padding,
!    Shift to next section of each array kg
   ikg=ikg+npw1
 end do !  End of the loop over k points

 if(mode_paral == 'PERS') then
   call xmpi_sum(npwarr,mpi_enreg%comm_kpt,ierr)
 end if

 if (mpi_enreg%nproc>1) then
   call wrtout(std_out,' kpgio: loop on k-points done in parallel','COLL')
 end if

!XG030513 MPIWF : now, one should sum npwarr over all processors
!of the WF group, to get npwtot (to be spread on all procs of the WF group
 npwtot(:)=npwarr(:)

!Taking into account istwfk
 do ikpt=1,nkpt
   if(istwfk(ikpt)>1)then
     if(istwfk(ikpt)==2)then
       npwtot(ikpt)=2*npwtot(ikpt)-1
     else
       npwtot(ikpt)=2*npwtot(ikpt)
     end if
   end if
 end do

!DEBUG
!write(std_out,*)' kpgio : exit '
!ENDDEBUG

end subroutine kpgio

!!***
