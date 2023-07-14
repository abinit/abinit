!!****m* ABINIT/m_out_spg_anal
!! NAME
!!  m_out_spg_anal
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2022 ABINIT group (DCA, XG, GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_out_spg_anal

 use defs_basis
 use m_results_out
 use m_dtset
 use m_abicore
 use m_errors
 use m_xomp
 use m_xmpi
#if defined HAVE_NETCDF
 use netcdf
#endif
 use m_outvar_a_h
 use m_outvar_i_n
 use m_outvar_o_z

 use m_parser,    only : ab_dimensions
 use m_nctk,      only : create_nc_file
 use m_symfind,   only : symfind_expert, symanal, symlatt
 use m_geometry,  only : metric, mkrdim
 use m_spgdata,   only : prtspgroup

 implicit none

 private
!!***

 public :: out_spg_anal
!!***

contains
!!***

!!****f* ABINIT/out_spg_anal
!! NAME
!! out_spg_anal
!!
!! FUNCTION
!! Perform final spacegroup analysis of the results of ABINIT, for each dataset.
!!
!! INPUTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number for echoed output
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!   one data set. Use for most dimensioned arrays.
!!  results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including evolving variables
!!
!! OUTPUT
!!  Only writing
!!
!! NOTES
!! Note that this routine is called only by the processor me==0 .
!! In consequence, no use of message and wrtout routine.
!!
!! SOURCE

subroutine out_spg_anal(dtsets,iout,ndtset,ndtset_alloc,results_out)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout
 integer,intent(in) :: ndtset,ndtset_alloc
!arrays
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)

!Local variables-------------------------------
!scalars
 integer :: idtset,iimage,invar_z,jdtset,msym,mu,natom,nimage,nptsym,nsym,ptgroupma,spgroup
 real(dp) :: tolsym,ucvol
 character(len=500) :: msg
!arrays
 integer :: bravais(11)
 integer, allocatable :: ptsymrel(:,:,:),symafm(:),symrel(:,:,:)
 real(dp) :: acell(3),genafm(3),gmet(3,3),gprimd(3,3),rmet(3,3),rprim(3,3),rprimd(3,3)
 real(dp), allocatable :: tnons(:,:),xred(:,:)

! *************************************************************************

 msym=dtsets(1)%maxnsym
 if(ndtset_alloc>1)then
   do idtset=2,ndtset_alloc
     msym=max(dtsets(idtset)%maxnsym,msym)
   end do
 end if
 ABI_MALLOC(ptsymrel,(3,3,msym))
 ABI_MALLOC(symafm,(msym))
 ABI_MALLOC(symrel,(3,3,msym))
 ABI_MALLOC(tnons,(3,msym))

 do idtset=1,ndtset_alloc

   tolsym=dtsets(idtset)%tolsym
   natom=dtsets(idtset)%natom
   nimage=results_out(idtset)%nimage
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
   ABI_MALLOC(xred,(3,natom))

   do iimage=1,nimage 

     acell=results_out(idtset)%acell(:,iimage)
     rprim=results_out(idtset)%rprim(:,:,iimage)
     xred=results_out(idtset)%xred(:,:,iimage)
     call mkrdim(acell,rprim,rprimd)
     call metric(gmet,gprimd,dev_null,rmet,rprimd,ucvol)

     !From rprimd and tolsym, compute bravais, nptsym and ptsymrel (with maximum size msym).
     call symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)

     invar_z=0 ; if(dtsets(idtset)%jellslab/=0 .or. dtsets(idtset)%nzchempot/=0)invar_z=2
     call symfind_expert(gprimd,msym,natom,nptsym,dtsets(idtset)%nspden,nsym,&
       dtsets(idtset)%pawspnorb,dtsets(idtset)%prtvol,ptsymrel,dtsets(idtset)%spinat,symafm,symrel,&
       tnons,tolsym,dtsets(idtset)%typat,dtsets(idtset)%usepaw,xred,&
       chrgat=dtsets(idtset)%chrgat,nucdipmom=dtsets(idtset)%nucdipmom,&
       invardir_red=dtsets(idtset)%field_xred,invar_z=invar_z)

     call symanal(bravais,dtsets(idtset)%chkprim,genafm,msym,nsym,ptgroupma,rprimd,spgroup,symafm,symrel,tnons,tolsym)

!    Should modify prtspgroup to allow echo of iimage, as optional argument
     call prtspgroup(bravais,genafm,iout,jdtset,ptgroupma,spgroup)

   enddo ! iimage

   ABI_FREE(xred)

 enddo ! idtset

!###########################################################
!## Deallocations and cleaning

 ABI_FREE(ptsymrel)
 ABI_FREE(symafm)
 ABI_FREE(symrel)
 ABI_FREE(tnons)

 write(msg,'(a,80a)')ch10,('=',mu=1,80)
 call wrtout(iout,msg,'COLL')

!**************************************************************************

end subroutine out_spg_anal
!!***

end module m_out_spg_anal
!!***
