!{\src2tex{textfont=tt}}
!!****f* ABINIT/invars2m
!! NAME
!! invars2m
!!
!! FUNCTION
!! Initialisation phase - main input routine.
!! Big loop on the datasets :
!! - for each of the datasets, write one line about the crystallographic data
!! - call invars2, that read the eventual single dataset input values ;
!! - compute mgfft,mpw,nfft,... for this data set ;
!! - compute quantities for the susceptibility computation
!! - compute the memory needs for this data set.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  bravais_(11,0:ndtset_alloc)=characteristics of Bravais lattice
!!  iout=unit number of output file
!!  lenstr=actual length of string
!!  mband_upper_(0:ndtset_alloc)=list of mband_upper values
!!  msym=default maximal number of symmetries
!!  ndtset= number of datasets to be read; if 0, no multi-dataset mode
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!      one data set.
!!  npsp=number of pseudopotentials
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!  string*(*)=character string containing all the input data.
!!   Initialized previously in instrng.
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are initialized here, while other were already
!!   initialized previously.
!!
!! NOTES
!! The outputs of this routine are the values of input variables,
!! their default value is stored at the index 0 of the last dimension
!! of their multi-dataset representation.
!!
!! PARENTS
!!      m_ab7_invars_f90
!!
!! CHILDREN
!!      invars2,metric,mkrdim,setshells
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine invars2m(dtsets,iout,lenstr,mband_upper_,msym,ndtset,ndtset_alloc,npsp,pspheads,string)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi

 use m_geometry,   only : mkrdim, metric
 use m_gsphere,    only : setshells

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'invars2m'
 use interfaces_57_iovars, except_this_one => invars2m
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,lenstr,msym,ndtset,ndtset_alloc,npsp
 character(len=*),intent(in) :: string
!arrays
 integer,intent(in) :: mband_upper_(0:ndtset_alloc)
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables -------------------------------
!scalars
 integer :: idtset,jdtset,mband_upper,nsheps,nshsigx,nshwfn,usepaw
 real(dp) :: ucvol
!arrays
 integer :: bravais(11)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: zionpsp(:)

!*************************************************************************

 do idtset=1,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
   bravais(:)=dtsets(idtset)%bravais(:)
   mband_upper  =mband_upper_(idtset)
   usepaw=dtsets(idtset)%usepaw
!  Allocate arrays
   ABI_ALLOCATE(zionpsp,(npsp))
   zionpsp(:)=pspheads(1:npsp)%zionpsp

!  Here, nearly all the remaining input variables are initialized
   call invars2(bravais,dtsets(idtset),iout,jdtset,lenstr,&
&   mband_upper,msym,npsp,string,usepaw,zionpsp)

   ABI_DEALLOCATE(zionpsp)

   call mkrdim(dtsets(idtset)%acell_orig(1:3,1),dtsets(idtset)%rprim_orig(1:3,1:3,1),rprimd)
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!  For GW or BSE calculations, we only use (npwwfn|ecutwfn) G-vectors read from the KSS file,
!  therefore the FFT box for the density should be defined according to ecut=ecutwfn.
   if ( ANY(dtsets(idtset)%optdriver == (/RUNL_SCREENING,RUNL_SIGMA,RUNL_BSE/)) ) then

     nshwfn=0
     call setshells(dtsets(idtset)%ecutwfn,dtsets(idtset)%npwwfn,nshwfn,&
&     dtsets(idtset)%nsym,gmet,gprimd,dtsets(idtset)%symrel,'wfn',ucvol)

!    MG: Hack to avoid portability problems under gfortran and g95:
!    getng and getmpw are indeed quite sensitive if ecut is small
!    and, in the GW tests, mpw and ngfft might depend on the compiler used.
!    the problem shows up if we use npwwfn instead of ecutwfn, a good
!    reason for removing npwwfn!
     dtsets(idtset)%ecutwfn=dtsets(idtset)%ecutwfn-tol14
!    MG: This is a kind of a hack, but the problem is ecutwfn that is too much redundant!
     dtsets(idtset)%ecut=dtsets(idtset)%ecutwfn

!    Close the shell for (W|chi0)
     nsheps=0
     call setshells(dtsets(idtset)%ecuteps,dtsets(idtset)%npweps,nsheps,&
&     dtsets(idtset)%nsym,gmet,gprimd,dtsets(idtset)%symrel,'eps',ucvol)

!    Close the shell for the exchange term.
     nshsigx=0
     call setshells(dtsets(idtset)%ecutsigx,dtsets(idtset)%npwsigx,nshsigx,&
&     dtsets(idtset)%nsym,gmet,gprimd,dtsets(idtset)%symrel,'sigx',ucvol)

   end if ! (SIGMA|SCREENING|SCGW|BSE)

 end do

end subroutine invars2m
!!***
