!{\src2tex{textfont=tt}}
!!****f* ABINIT/macroin
!! NAME
!! macroin
!!
!! FUNCTION
!! Treat "macro" input variables, that can :
!! - initialize several other input variables for one given dataset
!! - initialize several other input variables for a set of datasets.
!! Note that the treatment of these different types of macro input variables is different.
!! Documentation of such input variables is very important, including the
!! proper echo, in the output file, of what such input variables have done. 
!!
!! Important information : all the "macro" input variables should be properly
!! identifiable to be so, and it is proposed to make them start with the string "macro".
!!
!! COPYRIGHT
!! Copyright (C) 2009-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ndtset_alloc=number of datasets, corrected for allocation of at
!!               least one data set.
!!  ecut_tmp(3,2,10)= possible ecut values as read in psp files
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are given a value here.
!!   The dataset with number 0 should NOT be modified in the present routine.
!!
!! PARENTS
!!      m_ab7_invars_f90
!!
!! CHILDREN
!!      intagm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine macroin(dtsets,ecut_tmp,lenstr,ndtset_alloc,string)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'macroin'
 use interfaces_42_parser
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndtset_alloc,lenstr
 character(len=*),intent(inout) :: string
!arrays
 real(dp),intent(in) :: ecut_tmp(3,2,10)
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc) !vz_i ziontypat

!Local variables -------------------------------
!scalars
 integer :: idtset,iatom,jdtset,marr,tread
!!arrays
 integer,allocatable :: intarr(:)
 real(dp) :: ecutmax(3),ecutdgmax(3)
 real(dp),allocatable :: dprarr(:)
 character(len=500) :: message
!******************************************************************

 do idtset=1,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset 
   if (dtsets(idtset)%macro_uj>0) then
     dtsets(idtset)%irdwfk   = 1        ! preconverged wave function compulsory 
!    dtsets(idtset)%nline    = maxval((/ int(dtsets(idtset)%natom/2) , 6 /))   ! using default value: \DeltaU< 1%
!    dtsets(idtset)%nnsclo   = 4        ! using default value: \DeltaU< 1% 
     dtsets(idtset)%tolvrs   = 10d-8    ! convergence on the potential; 10d-8^= 10d-5 on occupation 
     dtsets(idtset)%diemix   = 0.45_dp  ! fastest convergence: dn= E^(-istep * 0.229 )
     dtsets(idtset)%dmatpuopt= 3        ! normalization of the occupation operator 
!    dtsets(idtset)%nstep    = 255      ! expected convergence after 10 \pm 3, 30 as in default normally suficient
!    dtsets(idtset)%iscf     = 17       ! mixing on potential, 17: default for PAW
   end if ! macro_uj

  !Read parameters
   marr=dtsets(idtset)%npsp;if (dtsets(idtset)%npsp<3) marr=3
   marr=max(marr,dtsets(idtset)%nimage)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
   
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),"accuracy",tread,'INT')

   ecutmax=-one
   ecutdgmax=-one
   do iatom=1,dtsets(idtset)%natom
     ecutmax(:)=max(ecutmax(:),ecut_tmp(:,1,dtsets(idtset)%typat(iatom)))
     ecutdgmax(:)=max(ecutdgmax(:),ecut_tmp(:,2,dtsets(idtset)%typat(iatom)))
   end do

   if(tread==1) then
     dtsets(idtset)%accuracy=intarr(1)
     if (dtsets(idtset)%accuracy==1) then
       if (ecutmax(1)>zero) dtsets(idtset)%ecut=ecutmax(1)
       if (ecutdgmax(1)>zero.and.dtsets(idtset)%usepaw==1) dtsets(idtset)%pawecutdg=ecutdgmax(1)
       dtsets(idtset)%boxcutmin=1.5_dp
       if (dtsets(idtset)%usepaw==1) then
         dtsets(idtset)%bxctmindg=1.5_dp
         dtsets(idtset)%pawxcdev=1
         dtsets(idtset)%pawmixdg=0
         dtsets(idtset)%pawovlp=10
         dtsets(idtset)%pawnhatxc=0
         dtsets(idtset)%mqgriddg=0
       end if
       dtsets(idtset)%mqgrid=0
       dtsets(idtset)%tolimg=5.0d-5
       dtsets(idtset)%tolvrs=tol3
       dtsets(idtset)%tolmxf=1.0d-3
       dtsets(idtset)%toldff=zero
       dtsets(idtset)%optforces=1
       dtsets(idtset)%timopt=0
       dtsets(idtset)%npulayit=4
       dtsets(idtset)%nstep=30
       dtsets(idtset)%prteig=0
       dtsets(idtset)%prtden=0
     else if (dtsets(idtset)%accuracy==2) then
       if (ecutmax(2)>zero) dtsets(idtset)%ecut=ecutmax(2)
       if (ecutdgmax(2)>zero.and.dtsets(idtset)%usepaw==1) dtsets(idtset)%pawecutdg=ecutdgmax(2)
       dtsets(idtset)%boxcutmin=1.8_dp
       if (dtsets(idtset)%usepaw==1) then
         dtsets(idtset)%bxctmindg=1.8_dp
         dtsets(idtset)%pawxcdev=1
         dtsets(idtset)%pawmixdg=0
         dtsets(idtset)%pawovlp=7
         dtsets(idtset)%pawnhatxc=1
         dtsets(idtset)%mqgriddg=0
       end if
       dtsets(idtset)%mqgrid=0
       dtsets(idtset)%tolimg=5.0d-5
       dtsets(idtset)%tolvrs=tol5
       dtsets(idtset)%tolmxf=5.0d-4
       dtsets(idtset)%toldff=zero
       dtsets(idtset)%optforces=1
       dtsets(idtset)%timopt=0
       dtsets(idtset)%npulayit=7
       dtsets(idtset)%nstep=30
       dtsets(idtset)%prteig=0
       dtsets(idtset)%prtden=0 
     else if (dtsets(idtset)%accuracy==3) then
       if (ecutmax(2)>zero) dtsets(idtset)%ecut=ecutmax(2)
       if (ecutdgmax(2)>zero.and.dtsets(idtset)%usepaw==1) dtsets(idtset)%pawecutdg=ecutdgmax(2)
       dtsets(idtset)%boxcutmin=1.8_dp
       if (dtsets(idtset)%usepaw==1) then
         dtsets(idtset)%bxctmindg=1.8_dp
         dtsets(idtset)%pawxcdev=1
         dtsets(idtset)%pawmixdg=0
         dtsets(idtset)%pawovlp=7
         dtsets(idtset)%pawnhatxc=1
         dtsets(idtset)%mqgriddg=0
       end if
       dtsets(idtset)%mqgrid=0
       dtsets(idtset)%tolimg=5.0d-5
       dtsets(idtset)%tolvrs=tol7
       dtsets(idtset)%tolmxf=1.0d-4
       dtsets(idtset)%toldff=zero
       dtsets(idtset)%optforces=2
       dtsets(idtset)%timopt=1
       if(xmpi_paral==1) dtsets(idtset)%timopt = 0
       dtsets(idtset)%npulayit=7
       dtsets(idtset)%nstep=30
       dtsets(idtset)%prteig=1
       dtsets(idtset)%prtden=1
     else if (dtsets(idtset)%accuracy==4) then
       if (ecutmax(3)>zero) dtsets(idtset)%ecut=ecutmax(3)
       if (ecutdgmax(3)>zero.and.dtsets(idtset)%usepaw==1) dtsets(idtset)%pawecutdg=ecutdgmax(3)
       dtsets(idtset)%boxcutmin=two
       if (dtsets(idtset)%usepaw==1) then
         dtsets(idtset)%bxctmindg=two
         dtsets(idtset)%pawxcdev=1
         dtsets(idtset)%pawmixdg=0
         dtsets(idtset)%pawovlp=5
         dtsets(idtset)%pawnhatxc=1
         dtsets(idtset)%mqgriddg=0
       end if
       dtsets(idtset)%mqgrid=0
       dtsets(idtset)%tolimg=5.0d-5
       dtsets(idtset)%tolvrs=tol9
       dtsets(idtset)%tolmxf=5.0d-5
       dtsets(idtset)%toldff=zero
       dtsets(idtset)%optforces=2
       dtsets(idtset)%timopt=1
       if(xmpi_paral==1) dtsets(idtset)%timopt = 0
       dtsets(idtset)%npulayit=7
       dtsets(idtset)%nstep=30
       dtsets(idtset)%prteig=1
       dtsets(idtset)%prtden=1
     else if (dtsets(idtset)%accuracy==5) then
       if (ecutmax(2)>zero) dtsets(idtset)%ecut=ecutmax(2)
       if (ecutdgmax(2)>zero.and.dtsets(idtset)%usepaw==1) dtsets(idtset)%pawecutdg=ecutdgmax(2)
       dtsets(idtset)%boxcutmin=two
       if (dtsets(idtset)%usepaw==1) then
         dtsets(idtset)%bxctmindg=two
         dtsets(idtset)%pawxcdev=2
         dtsets(idtset)%pawmixdg=1
         dtsets(idtset)%pawovlp=5
         dtsets(idtset)%pawnhatxc=1
         dtsets(idtset)%mqgriddg=0
       end if
       dtsets(idtset)%mqgrid=0
       dtsets(idtset)%tolimg=5.0d-5
       dtsets(idtset)%tolvrs=tol10
       dtsets(idtset)%tolmxf=1.0d-6
       dtsets(idtset)%toldff=zero
       dtsets(idtset)%optforces=2
       dtsets(idtset)%timopt=1
       if(xmpi_paral==1) dtsets(idtset)%timopt = 0
       dtsets(idtset)%npulayit=15
       dtsets(idtset)%nstep=50
       dtsets(idtset)%prteig=1
       dtsets(idtset)%prtden=1
     else if (dtsets(idtset)%accuracy==6) then
       if (ecutmax(3)>zero) dtsets(idtset)%ecut=ecutmax(3)
       if (ecutdgmax(3)>zero.and.dtsets(idtset)%usepaw==1) dtsets(idtset)%pawecutdg=ecutdgmax(3)
       dtsets(idtset)%boxcutmin=two
       if (dtsets(idtset)%usepaw==1) then
         dtsets(idtset)%bxctmindg=two
         dtsets(idtset)%pawxcdev=2
         dtsets(idtset)%pawmixdg=1
         dtsets(idtset)%pawovlp=5
         dtsets(idtset)%pawnhatxc=1
         dtsets(idtset)%mqgriddg=0
       end if
       dtsets(idtset)%mqgrid=0
       dtsets(idtset)%tolimg=5.0d-5
       dtsets(idtset)%tolvrs=tol12
       dtsets(idtset)%tolmxf=1.0d-6
       dtsets(idtset)%toldff=zero
       dtsets(idtset)%optforces=2
       dtsets(idtset)%timopt=1
       if(xmpi_paral==1) dtsets(idtset)%timopt = 0
       dtsets(idtset)%npulayit=15
       dtsets(idtset)%nstep=50
       dtsets(idtset)%prteig=1
       dtsets(idtset)%prtden=1
     elseif(dtsets(idtset)%accuracy>6)then
       write(message, '(a,a,a)' )&
&       'accuracy >6 is forbiden !',ch10,&
&       'Action : check your input data file.'
       MSG_ERROR(message)
     end if
   else
     if (ecutmax(3)>zero) dtsets(idtset)%ecut=ecutmax(3)
   end if
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
 end do

end subroutine macroin
!!***
