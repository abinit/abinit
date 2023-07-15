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
!! Compare with the initial one, and perform analysis, with adequate warning if there was a change.
!! Possibly echo spacegroup for all dtsets and possibly all images
!!
!! INPUTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number for echoed output
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!   one data set. Use for most dimensioned arrays.
!!  echo_spgroup = option for analysis 
!!      1 => write ;  2 => echo of spacegroup for all dtsets and possibly all images
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

subroutine out_spg_anal(dtsets,echo_spgroup,iout,ndtset,ndtset_alloc,results_out)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: echo_spgroup,iout
 integer,intent(in) :: ndtset,ndtset_alloc
!arrays
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)

!Local variables-------------------------------
!scalars
 integer, save :: counter0=1, counter1=1
 integer :: idtset,iimage,invar_z,jdtset,msym,mu,natom,nimage,nptsym,nsym
 integer :: ptgroupma,spgroup,symmetry_changed
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
     call symlatt(bravais,std_out,msym,nptsym,ptsymrel,rprimd,tolsym)

     invar_z=0 ; if(dtsets(idtset)%jellslab/=0 .or. dtsets(idtset)%nzchempot/=0)invar_z=2
     call symfind_expert(gprimd,msym,natom,nptsym,dtsets(idtset)%nspden,nsym,&
       dtsets(idtset)%pawspnorb,dtsets(idtset)%prtvol,ptsymrel,dtsets(idtset)%spinat,symafm,symrel,&
       tnons,tolsym,dtsets(idtset)%typat,dtsets(idtset)%usepaw,xred,&
       chrgat=dtsets(idtset)%chrgat,nucdipmom=dtsets(idtset)%nucdipmom,&
       invardir_red=dtsets(idtset)%field_xred,invar_z=invar_z)

     !Set chkprim to 0, to allow detecting increase of multiplicity
     call symanal(bravais,0,genafm,msym,nsym,ptgroupma,rprimd,spgroup,symafm,symrel,tnons,tolsym)

     symmetry_changed=0
     if( nsym/=dtsets(idtset)%nsym .or. spgroup/=dtsets(idtset)%spgroup .or. &
&        ptgroupma/=dtsets(idtset)%ptgroupma)then
       symmetry_changed=1
     endif

     if(symmetry_changed==1)then
       if(echo_spgroup==0 .and. counter0==1)then
         write(msg,'(8a)')ch10,' The spacegroup number, the magnetic point group, and/or the number of symmetries',ch10,&
&         ' have changed between the initial recognition based on the input file',ch10,&
&         ' and a postprocessing based on the final acell, rprim, and xred.',ch10,&
&         ' More details in the log file.'
         call wrtout(iout,msg,'COLL')
         counter0=0
       endif
       if(echo_spgroup==1 .and. counter1==1)then
         write(msg,'(10a)')ch10,' The spacegroup number, the magnetic point group, and/or the number of symmetries',ch10,&
&         ' have changed between the initial recognition based on the input file',ch10,&
&         ' and a postprocessing based on the final acell, rprim, and xred.',ch10,&
&         ' These modifications are detailed below.',ch10,&
&         ' The updated tnons, symrel or symrel have NOT been reported in the final echo of variables after computation.'
         call wrtout(iout,msg,'COLL')
         write(msg,'(5a)')' Such change of spacegroup, or magnetic point group might happen in several cases.',ch10,&
&         ' (1) If spgroup (+ptgroupma) defined in the input file, but the actual groups are supergroups of these; ',ch10,&
&         ' (2) If symrel, tnons (+symafm) defined in the input file, while the system is more symmetric; '
         call wrtout(iout,msg,'COLL')
         write(msg,'(5a)')&
&         ' (3) If the geometry has been optimized and the final structure is more symmetric than the initial one;',ch10,
&         ' (4) In case of GW of BSE calculation with inversion symmetry, as nsym has been reduced in such',ch10,&
&         '       dataset, excluding the improper symmetry operations (with determinant=-1), but not in the postprocessing.'
         call wrtout(iout,msg,'COLL')
         write(msg,'(5a)')' In some case, the recognition of symmetries strongly depends on the value of tolsym.',ch10,&
          ' You might investigate its effect by restarting abinit based on the final acell, rprim and xred,',ch10,&
&         ' and different values for tolsym.'
         counter1=0
       endif
     endif

!    Echo the spacegroup (and ptgroupma) if requested, and give more information if the symmetry changed.
     if(echo_spgroup==1)then

       if(symmetry_changed==0)then
         if(nimage==1)then
           call prtspgroup(bravais,genafm,iout,jdtset,ptgroupma,spgroup)
         else
           call prtspgroup(bravais,genafm,iout,jdtset,ptgroupma,spgroup,iimage=iimage)
         endif
       else 
         write(msg,'(2a,3i5)')ch10,' Initial data. jdtset, iimage, nsym=',jdtset,iimage,dtsets(idtset)%nsym
         call wrtout(iout,msg,'COLL')
         if(nimage==1)then
           call prtspgroup(dtsets(idtset)%bravais,dtsets(idtset)%genafm,iout,jdtset,&
&           dtsets(idtset)%ptgroupma,dtsets(idtset)%spgroup)
         else
           call prtspgroup(dtsets(idtset)%bravais,dtsets(idtset)%genafm,iout,jdtset,&
&           dtsets(idtset)%ptgroupma,dtsets(idtset)%spgroup,iimage=iimage)
         endif
         write(msg,'(a,3i5)')' Final data.   jdtset, iimage, nsym=',jdtset,iimage,nsym
         call wrtout(iout,msg,'COLL')
         if(nimage==1)then
           call prtspgroup(bravais,genafm,iout,jdtset,ptgroupma,spgroup)
         else
           call prtspgroup(bravais,genafm,iout,jdtset,ptgroupma,spgroup,iimage=iimage)
         endif
       endif

     end if ! echo_spgroup==1

   enddo ! iimage

   ABI_FREE(xred)

 enddo ! idtset

!###########################################################
!## Deallocations and cleaning

 ABI_FREE(ptsymrel)
 ABI_FREE(symafm)
 ABI_FREE(symrel)
 ABI_FREE(tnons)

 if(echo_spgroup==1)then
   write(msg,'(a,80a)')ch10,('=',mu=1,80)
   call wrtout(iout,msg,'COLL')
 endif

!**************************************************************************

end subroutine out_spg_anal
!!***

end module m_out_spg_anal
!!***
