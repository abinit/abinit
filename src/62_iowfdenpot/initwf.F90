!{\src2tex{textfont=tt}}
!!****f* ABINIT/initwf
!!
!! NAME
!! initwf
!!
!! FUNCTION
!! Initialization of wavefunctions.
!! If formeig==1, and partially filled case, I am not sure that the eig_k are initialized properly ...
!! formeig option (format of the eigenvalues and eigenvector) :
!!   0 => ground-state format (initialisation of
!!        eigenvectors with random numbers, vector of eigenvalues)
!!   1 => respfn format (initialisation of
!!        eigenvectors with 0 s, hermitian matrix of eigenvalues)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! formeig=see above
!! headform=header format (might be needed to read the block of wfs)
!! icg=shift to be given to the location of the data in the array cg
!! ikpt= number of the k point of which the wf is initialised
!! spin=spin index
!! mcg=dimension of the cg array
!! mpi_enreg=information about MPI parallelization
!! nband_k=number of bands at this particular k point
!! nkpt=number of k points
!! npw=number of plane waves
!! nspinor=number of spinorial components of the wavefunctions (on current proc)
!! wff1=structure info for file containing wavefunctions (when needed)
!!
!! OUTPUT
!! cg(2,mcg)=complex wf array
!!  if ground state format (formeig=0):
!!    eig_k(nband_k)=list of eigenvalues (input or init to large number), hartree
!!  if respfn format (formeig=1):
!!    eig_k(2*nband_k*nband_k)= matrix of eigenvalues (input or init to large number), hartree
!!
!! SIDE EFFECTS
!! Input/output:
!! occ_k(nband_k)=list of occupations (input or left to their initial value)
!! ikptsp_old=number of the previous spin-k point, or 0 if first call of present file
!!
!! PARENTS
!!      wfsinp
!!
!! CHILDREN
!!      rwwf,timab,wffreadskipk,wfk_close,wfk_open_read,wfk_read_band_block
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initwf(cg,eig_k,formeig,headform,icg,ikpt,ikptsp_old,&
&  spin,mcg,mpi_enreg,nband_k,nkpt,npw,nspinor,occ_k,wff1)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_wfk
 use m_wffile

 use m_io_tools,     only : get_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initwf'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_56_io_mpi
 use interfaces_62_iowfdenpot, except_this_one => initwf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: formeig,headform,icg,ikpt,spin,mcg,nband_k,nkpt,npw,nspinor
 integer,intent(inout) :: ikptsp_old
 type(MPI_type),intent(inout) :: mpi_enreg
 type(wffile_type),intent(inout) :: wff1
!arrays
 real(dp),intent(inout) :: occ_k(nband_k)
 real(dp),intent(inout) :: cg(2,mcg),eig_k((2*nband_k)**formeig*nband_k) !vz_i

!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt_max=50
 integer :: ikpt0,nband_disk,tim_rwwf 
 character(len=500) :: msg
!arrays
 integer,allocatable :: kg_dum(:,:)
 real(dp) :: tsec(2)
#if 0
 integer :: iomode,comm,funt,ierr
 type(wfk_t) :: Wfk
#endif

! *************************************************************************

!DEBUG
!write(std_out,*)' initwf : enter, ikptsp_old,ikpt,spin,nkpt= ',ikptsp_old,ikpt,spin,nkpt
!stop
!ENDDEBUG

#if 0
 MSG_WARNING("Entering new IO section")
!call WffClose(wff1,ierr)
 comm   = MPI_enreg%comm_cell
 iomode = iomode_from_fname(wff1%fname)
 call wfk_open_read(Wfk,wff1%fname,formeig,iomode,get_unit(),comm)
 call wfk_read_band_block(Wfk,(/1,nband_k/),ikpt,spin,xmpio_at,cg_k=cg(1:,icg:),eig_k=eig_k,occ_k=occ_k)
 call wfk_close(Wfk)
!call clsopn(wff1)
 RETURN
#endif

 call timab(770,1,tsec)
 call timab(771,1,tsec)

 ABI_ALLOCATE(kg_dum,(3,0))

!Skip wavefunctions for k-points not treated by this proc.
!(from ikptsp_old+1 to ikpt+(spin-1)*nkpt-1)
 if (ikptsp_old<ikpt+(spin-1)*nkpt-1) then
   do ikpt0=ikptsp_old+1,ikpt+(spin-1)*nkpt-1
     call WffReadSkipK(formeig,headform,ikpt0,spin,mpi_enreg,wff1)
   end do
 end if

!DEBUG
!write(std_out,*)' initwf : before rwwf'
!write(std_out,*)' formeig,icg,ikpt,spin=',formeig,icg,ikpt,spin
!write(std_out,*)' nband_k,nband_disk,npw,nspinor=',nband_k,nband_disk,npw,nspinor
!write(std_out,*)' unwff1=',unwff1
!stop
!ENDDEBUG

 if(mpi_enreg%paralbd==0)tim_rwwf=2
 if(mpi_enreg%paralbd==1)tim_rwwf=20

 call timab(771,2,tsec)

 call rwwf(cg,eig_k,formeig,headform,icg,ikpt,spin,kg_dum,nband_k,mcg,mpi_enreg,nband_k,nband_disk,&
& npw,nspinor,occ_k,1,0,tim_rwwf,wff1)

 call timab(772,1,tsec)

 if (ikpt<=nkpt_max) then
   write(msg,'(3(a,i0))')' initwf: disk file gives npw= ',npw,' nband= ',nband_disk,' for kpt number= ',ikpt
   call wrtout(std_out,msg,'PERS')
 else if (ikpt==nkpt_max+1) then
   call wrtout(std_out,' initwf: the number of similar message is sufficient... stop printing them','PERS')
 end if

!Check the number of bands on disk file against desired number. These are not required to agree)
 if (nband_disk/=nband_k) then
   write(msg,'(2(a,i0),3a,i0,3a)')&
&   'For kpt number ',ikpt,' disk file has ',nband_disk,' bands',ch10,&
&   'but input file gave nband= ',nband_k,'.',ch10,&
&   'This is not fatal. Bands are skipped or filled with random numbers.'
   MSG_COMMENT(msg)
 end if

 if (ikpt<=nkpt_max) then
   write(msg,'(a,i0,a)')' initwf: ',nband_disk,' bands have been initialized from disk'
   call wrtout(std_out,msg,'PERS')
 end if

 ikptsp_old=ikpt+(spin-1)*nkpt

 ABI_DEALLOCATE(kg_dum)

 call timab(772,2,tsec)
 call timab(770,2,tsec)

end subroutine initwf
!!***
