!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_uj
!! NAME
!!  m_paw_uj
!!
!! FUNCTION
!!  This module contains several routines relevant only for automatic determination of U
!!    in PAW+U context (linear response method according to Phys. Rev. B 71, 035105)
!!
!! COPYRIGHT
!! Copyright (C) 2018-2019 ABINIT group (DJA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_uj

 use defs_basis
 use defs_abitypes
 use m_abicore
 use m_errors
 use m_linalg_interfaces
 use m_xmpi

 use m_pptools,       only : prmat
 use m_special_funcs, only : iradfnh
 use m_geometry,      only : shellstruct,ioniondist
 use m_parser,        only : prttagm
 use m_supercell,     only : mksupercell
 use m_pawrad,        only : pawrad_type
 use m_pawtab,        only : pawtab_type
 use m_paw_ij,        only : paw_ij_type
 use m_paral_atom,    only : get_my_atmtab, free_my_atmtab

 implicit none

 private

!public procedures.
 public :: pawuj_ini  ! Initialize dtpawuj datastructure
 public :: pawuj_free ! Deallocate dtpawuj datastructure
 public :: pawuj_det  ! Determine U (or J) parameter
 public :: pawuj_red  ! Store atomic occupancies, potential shift, positions

!private procedures.
!  chiscwrt:   distributes values of chi_org on chi_sc according to ion-ion distances
!  linvmat:    inverts real matrix inmat
!  lprtmat:    prints out the real matrix mmat
!  lcalcu:     prints out real the real matrice mmat
!  blow_pawuj: reads a real nxn matrice and appends lines n+1 and clumn n+1

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_uj/pawuj_ini
!! NAME
!! pawuj_ini
!!
!! FUNCTION
!!  Initialize dtpawuj datastructure for one SCF run
!!  Relevant only for automatic determination of U in PAW+U context
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dtpawuj(0:ndtpawuj) (initialization of fields vsh, occ, iuj,nnat)
!!
!! PARENTS
!!      pawuj_drive,ujdet
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawuj_ini(dtpawuj,ndtset)

 implicit none

!Arguments ------------------------------------
!scalars
 integer                           :: ndtset
 type(macro_uj_type),intent(inout) :: dtpawuj(0:ndtset)

!Local variables -------------------------
!Variables for partial dos calculation
!scalars
 integer, parameter :: nwfchr=6
 integer            :: iuj,im1
! *********************************************************************

 DBG_ENTER("COLL")

 do iuj=0,ndtset
   !write(std_out,*)'pawuj_ini iuj ',iuj
   dtpawuj(iuj)%diemix=0.45_dp
   dtpawuj(iuj)%iuj=0
   dtpawuj(iuj)%nat=0
   dtpawuj(iuj)%ndtset=1
   dtpawuj(iuj)%nspden=1
   dtpawuj(iuj)%macro_uj=0
   dtpawuj(iuj)%option=1
   dtpawuj(iuj)%pawujat=1
   dtpawuj(iuj)%pawujga=one
   dtpawuj(iuj)%pawprtvol=1
   dtpawuj(iuj)%ph0phiint=one
   dtpawuj(iuj)%dmatpuopt=3
   dtpawuj(iuj)%pawujrad=3.0_dp
   dtpawuj(iuj)%pawrad=20.0_dp
   !Allocate arrays
   !write(std_out,*)'pawuj_ini before arrays'
   ABI_ALLOCATE(dtpawuj(iuj)%rprimd,(3,3))
   ABI_ALLOCATE(dtpawuj(iuj)%scdim,(3))
   ABI_ALLOCATE(dtpawuj(iuj)%wfchr,(nwfchr))
   dtpawuj(iuj)%rprimd=reshape((/ 1,0,0,0,1,0,0,0,1/),(/ 3,3 /))
   dtpawuj(iuj)%scdim=reshape((/ 250,0,0 /),(/3 /))
   dtpawuj(iuj)%wfchr=(/ (0,im1=1,nwfchr) /)
   if (iuj>0) then
     dtpawuj(iuj)%iuj=-1
   end if
 end do

 DBG_EXIT("COLL")

end subroutine pawuj_ini
!!***

!----------------------------------------------------------------------

!!****f* m_paw_uj/pawuj_free
!! NAME
!!  pawuj_free
!!
!! FUNCTION
!!   deallocate pawuj stuff
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      pawuj_drive,ujdet
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawuj_free(dtpawuj)

 implicit none

!Arguments -------------------------------
 type(macro_uj_type),intent(inout) :: dtpawuj

! *********************************************************************

 if (allocated(dtpawuj%scdim))    then
   ABI_DEALLOCATE(dtpawuj%scdim)
 end if
 if (allocated(dtpawuj%occ))      then
   ABI_DEALLOCATE(dtpawuj%occ)
 end if
 if (allocated(dtpawuj%rprimd))   then
   ABI_DEALLOCATE(dtpawuj%rprimd)
 end if
 if (allocated(dtpawuj%vsh))      then
   ABI_DEALLOCATE(dtpawuj%vsh)
 end if
 if (allocated(dtpawuj%xred))     then
   ABI_DEALLOCATE(dtpawuj%xred)
 end if
 if (allocated(dtpawuj%wfchr))    then
   ABI_DEALLOCATE(dtpawuj%wfchr)
 end if
 if (allocated(dtpawuj%zioneff))  then
   ABI_DEALLOCATE(dtpawuj%zioneff)
 end if

end subroutine pawuj_free
!!***

!----------------------------------------------------------------------
!!****f* m_paw_uj/pawuj_det
!! NAME
!!  pawuj_det
!!
!! FUNCTION
!!  From the complete dtpawuj-dataset determines U (or J) parameter for
!!  PAW+U calculations
!!  Relevant only for automatic determination of U in PAW+U context
!!
!! INPUTS
!!  dtpawuj=potential shifts (vsh) and atomic occupations (occ)
!!  ujdet_filename= Filename for write (Using NetCDF format)
!!
!! OUTPUT
!!  only printing
!!  (among other things a section in the ab.out that can be used for input in ujdet)
!!
!! PARENTS
!!      pawuj_drive,ujdet
!!
!! CHILDREN
!!      chiscwrt,lcalcu,mksupercell,prmat,prttagm,shellstruct,wrtout
!!
!! SOURCE

subroutine pawuj_det(dtpawuj,ndtpawuj,ujdet_filename,ures)

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 integer                        :: ndtpawuj
 type(macro_uj_type),intent(in) :: dtpawuj(0:ndtpawuj)
 real(dp),intent(out)           :: ures
 character(len=*),intent(in)    :: ujdet_filename

!Local variables-------------------------------
!scalars
 integer,parameter           :: natmax=2,nwfchr=6
 integer                     :: ii,jj,nat_org,jdtset,nspden,macro_uj,kdtset,marr
 integer                     :: im1,ndtuj,idtset, nsh_org, nsh_sc,nat_sc,maxnat
 integer                     :: pawujat,pawprtvol,pawujoption
 integer                     :: dmatpuopt,invopt
 real(dp)                    :: pawujga,ph0phiint,intg,fcorr,eyp

 character(len=500)          :: message
 character(len=2)            :: hstr
!arrays
 integer                     :: ext(3)
 real(dp)                    :: rprimd_sc(3,3),vsh(ndtpawuj),a(5),b(5)
 integer,allocatable         :: narrm(:)
 integer,allocatable         :: idum2(:,:),jdtset_(:),smult_org(:),smult_sc(:)
 real(dp),allocatable        :: chih(:,:),dqarr(:,:),dqarrr(:,:),dparr(:,:),dparrr(:,:),xred_org(:,:),drarr(:,:)
 real(dp),allocatable        :: magv_org(:),magv_sc(:),chi_org(:),chi0_org(:),chi0_sc(:), chi_sc(:), xred_sc(:,:)
 real(dp),allocatable        :: sdistv_org(:),sdistv_sc(:),distv_org(:),distv_sc(:)
 integer                     :: ncid=0
! *********************************************************************

 DBG_ENTER("COLL")

!write(std_out,*) 'pawuj 01'
!###########################################################
!### 01. Allocations

!Initializations
 ndtuj=count(dtpawuj(:)%iuj/=-1)-1 ! number of datasets initialized by pawuj_red
 ABI_ALLOCATE(jdtset_,(0:ndtuj))
 jdtset_(0:ndtuj)=pack(dtpawuj(:)%iuj,dtpawuj(:)%iuj/=-1)
 jdtset=maxval(dtpawuj(:)%iuj)

!DEBUG
!write(message,'(10(a,i3))')'pawuj_det jdtset ',jdtset,&
!& ' ndtuj ', ndtuj,' ndtpawuj ',ndtpawuj
!call wrtout(std_out,message,'COLL')
!call flush_unit(6)
!END DEBUG

 nspden=dtpawuj(jdtset)%nspden
 nat_org=dtpawuj(jdtset)%nat
 macro_uj=dtpawuj(jdtset)%macro_uj
 pawujat=dtpawuj(jdtset)%pawujat
 pawprtvol=dtpawuj(jdtset)%pawprtvol
 pawujga=dtpawuj(jdtset)%pawujga
 pawujoption=dtpawuj(jdtset)%option
 ph0phiint=dtpawuj(jdtset)%ph0phiint
 dmatpuopt=dtpawuj(jdtset)%dmatpuopt
 marr=maxval((/ 9, nspden*nat_org ,nat_org*3 /))
 eyp=2.5_dp ! for dmatpuopt==2 and 3
 if (dmatpuopt==1) eyp=eyp+3.0_dp
 if (dmatpuopt>=3) eyp=(eyp+3.0_dp-dmatpuopt)

 ABI_ALLOCATE(chih,(ndtpawuj,nat_org))
 ABI_ALLOCATE(idum2,(marr,0:ndtuj))
 ABI_ALLOCATE(drarr,(marr,0:ndtuj))
 ABI_ALLOCATE(magv_org,(nat_org))
 ABI_ALLOCATE(xred_org,(3,nat_org))
 ABI_ALLOCATE(chi0_org,(nat_org))
 ABI_ALLOCATE(chi_org,(nat_org))
 ABI_ALLOCATE(dparr,(marr,0:ndtuj))
 ABI_ALLOCATE(dparrr,(marr,0:ndtuj))
 ABI_ALLOCATE(dqarr,(marr,0:ndtuj))
 ABI_ALLOCATE(dqarrr,(marr,0:ndtuj))
 ABI_ALLOCATE(distv_org,(nat_org))
 ABI_ALLOCATE(narrm,(0:ndtuj))
 dparr=-one ;  dparrr=-one ;  dqarr=-one ;  dqarrr=-one
!DEBUG
!write(message,fmt='((a,i3,a))')'pawuj_det init sg'
!call wrtout(std_out,message,'COLL')
!END DEBUG
 idum2=1 ; drarr=one
 chih=zero

!write(std_out,*) 'pawuj 02'
!###########################################################
!### 02. Create the file UJDET.nc

 if(.false.)write(std_out,*)ujdet_filename ! This is for the abirules

!write(std_out,*) 'pawuj 03'
!###########################################################
!### 03. Write out the Input for UJDET

 write(message, '(3a)' ) ch10,&
& ' # input for ujdet, cut it using ''sed -n "/MARK/,/MARK/p" abi.out > ujdet.in ''------- ',ch10
 call wrtout(ab_out,message,'COLL')

 if (ndtuj/=4) then
   write (hstr,'(I0)') jdtset              ! convert integer to string
 else
   hstr=''
 end if
 if (ndtuj>=2.or.jdtset==1) then
   idum2(1,1:ndtuj)=4 !dtpawuj(:)%ndtuj
   call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid,ndtuj,'ndtset','INT',0)
 end if
 idum2(1,0:ndtuj)=pack(dtpawuj(:)%nat,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid,ndtuj,'nat'//trim(hstr),'INT',0)

 idum2(1,0:ndtuj)=pack(dtpawuj(:)%nspden,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid,ndtuj,'nspden'//trim(hstr),'INT',0)

 idum2(1,0:ndtuj)=pack(dtpawuj(:)%macro_uj,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid,ndtuj,'macro_uj'//trim(hstr),'INT',0)

 idum2(1,0:ndtuj)=pack(dtpawuj(:)%pawujat,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid,ndtuj,'pawujat'//trim(hstr),'INT',0)

 dparr(1,0:ndtuj)=pack(dtpawuj(:)%pawujga,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid,ndtuj,'pawujga'//trim(hstr),'DPR',0)

 dparr(1,0:ndtuj)=pack(dtpawuj(:)%pawujrad,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid,ndtuj,'pawujrad'//trim(hstr),'DPR',0)

 dparr(1,0:ndtuj)=pack(dtpawuj(:)%pawrad,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid,ndtuj,'pawrad'//trim(hstr),'DPR',0)

 dparr(1,0:ndtuj)=pack(dtpawuj(:)%ph0phiint,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid,ndtuj,'ph0phiint'//trim(hstr),'DPR',0)

 idum2(1,0:ndtuj)=pack(dtpawuj(:)%pawprtvol,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid,ndtuj,'pawprtvol'//trim(hstr),'INT',0)

 idum2(1,0:ndtuj)=pack(dtpawuj(:)%option,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid,ndtuj,'pawujopt'//trim(hstr),'INT',0)

 idum2(1,0:ndtuj)=pack(dtpawuj(:)%dmatpuopt,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid,ndtuj,'dmatpuopt'//trim(hstr),'INT',0)

 kdtset=0

 do idtset=0,ndtpawuj
!  DEBUG
!  write(message,fmt='((a,i3,a))')'pawuj_det m2, idtset ',idtset,ch10
!  call wrtout(std_out,message,'COLL')
!  call flush_unit(6)
!  END DEBUG
   if (dtpawuj(idtset)%iuj/=-1) then
     dparr(1:nspden*nat_org,kdtset)=reshape(dtpawuj(idtset)%vsh,(/nspden*nat_org/))
     dparrr(1:nspden*nat_org,kdtset)=reshape(dtpawuj(idtset)%occ,(/nspden*nat_org/))
!    DEBUG
!    write(message,fmt='((a,i3,a))')'pawuj_det m3, idtset ',idtset,ch10
!    call wrtout(std_out,message,'COLL')
!    write(std_out,*)' marr,narrm,ncid,ndtuj,nat_org,kdtset,idtset=',marr,narrm,ncid,ndtuj,nat_org,kdtset,idtset
!    write(std_out,*)' dtpawuj(idtset)%xred=',dtpawuj(idtset)%xred
!    call flush_unit(6)
!    END DEBUG
     dqarr(1:nat_org*3,kdtset)=reshape(dtpawuj(idtset)%xred,(/nat_org*3/))
     dqarrr(1:3*3,kdtset)=reshape(dtpawuj(idtset)%rprimd,(/3*3/))
     idum2(1:3,kdtset)=reshape(dtpawuj(idtset)%scdim,(/3/))
     drarr(1:nwfchr,kdtset)=reshape(dtpawuj(idtset)%wfchr,(/nwfchr/))
!    DEBUG
!    write(message,fmt='((a,i3,a))')'pawuj_det m4, idtset ',idtset,ch10
!    call wrtout(std_out,message,'COLL')
!    call flush_unit(6)
!    END DEBUG
     kdtset=kdtset+1
   end if
 end do

!DEBUG
!write(message,fmt='((a,i3,a))')'pawuj_det m5'
!call wrtout(std_out,message,'COLL')
!END DEBUG
 call prttagm(dparr,idum2,ab_out,jdtset_,2,marr,nspden*nat_org,narrm,ncid,ndtuj,'vsh'//trim(hstr),'DPR',0)
 call prttagm(dparrr,idum2,ab_out,jdtset_,2,marr,nspden*nat_org,narrm,ncid,ndtuj,'occ'//trim(hstr),'DPR',0)
!DEBUG
!write(message,fmt='((a,i3,a))')'pawuj_det m6'
!call wrtout(std_out,message,'COLL')
!END DEBUG
 call prttagm(dqarr,idum2,ab_out,jdtset_,2,marr,nat_org*3,narrm,ncid,ndtuj,'xred'//trim(hstr),'DPR',0)
 call prttagm(dqarrr,idum2,ab_out,jdtset_,2,marr,3*3,narrm,ncid,ndtuj,'rprimd'//trim(hstr),'DPR',0)
 call prttagm(dqarrr,idum2,ab_out,jdtset_,2,marr,3,narrm,ncid,ndtuj,'scdim'//trim(hstr),'INT',0)
 call prttagm(drarr,idum2,ab_out,jdtset_,2,marr,nwfchr,narrm,ncid,ndtuj,'wfchr'//trim(hstr),'DPR',0)
 ABI_DEALLOCATE(narrm)

 write(message, '( 15a )'  ) ch10,' # further possible options: ',ch10,&
& ' #    scdim    2 2 2 ',ch10,&
& ' #    mdist    10.0  ',ch10,&
& ' #  pawujga    2 '    ,ch10,&
& ' # pawujopt    2 '    ,ch10,&
& ' # pawujrad    3.0'   ,ch10,&
& ' # ------- end input for ujdet: end-MARK  -------- ',ch10
 call wrtout(ab_out,message,'COLL')

!write(std_out,*) 'pawuj 04'
!###########################################################
!### 04. Test

 if (ndtuj/=4)  return

 write(message, '(3a)' ) ch10,' ---------- calculate U, (J) start ---------- ',ch10
 call wrtout(ab_out,message,'COLL')

 if (all(dtpawuj(1:ndtpawuj)%pawujat==pawujat)) then
   write (message,fmt='(a,i3)') ' All pawujat  ok and equal to ',pawujat
   call wrtout(ab_out,message,'COLL')
 else
   write (message,fmt='(a,4i3,2a)') ' Differing values of pawujat were found: ',dtpawuj(1:ndtuj)%pawujat,ch10,&
&   'No determination of U.'
   call wrtout(ab_out,message,'COLL')
   return
 end if

 if (all(dtpawuj(1:ndtpawuj)%macro_uj==macro_uj)) then
   if (nspden==1) then
     write(message,fmt='(2a)') ' pawuj_det found nspden==1, determination',&
&     ' of U-parameter for unpol. struct. (non standard)'
   else if (macro_uj==1.and.nspden==2) then
     write(message,fmt='(2a)') ' pawuj_det: found macro_uj=1 and nspden=2:',&
&     ' standard determination of U-parameter'
   else if (macro_uj==2.and.nspden==2) then
     write(message,fmt='(2a)') ' pawuj_det: found macro_uj=2 and nspden=2:',&
&     ' determination of U on single spin channel (experimental)'

   else if (macro_uj==3.and.nspden==2) then
     write(message,fmt='(2a)') ' pawuj_det: found macro_uj=3 and nspden=2,',&
&     ' determination of J-parameter on single spin channel (experimental)'
   end if
   write (message,fmt='(a,i3,a,a)') ' All macro_uj ok and equal to ',macro_uj,ch10,trim(message)
   write (message,'(a,i3)') ' All macro_uj ok and equal to ',macro_uj
   call wrtout(ab_out,message,'COLL')
 else
   write (message,fmt='(a,10i3)') ' Differing values of macro_uj were found: ',dtpawuj(:)%macro_uj
   write (message,fmt='(3a)')trim(message),ch10,' No determination of U.'
   call wrtout(ab_out,message,'COLL')
   return
 end if

 if (macro_uj>1.and.nspden==1) then
   write (message,'(4a,2a)') ' U on a single spin channel (or J) can only be determined for nspden=2 ,',ch10,&
&   'No determination of U.'
   call wrtout(ab_out,message,'COLL')
   return
 end if

!Calculation of response matrix

 do jdtset=1,4
   if (nspden==1) then
     chih(jdtset,1:nat_org)=dtpawuj(jdtset)%occ(1,:)
   else if (macro_uj==1.and.nspden==2) then
     chih(jdtset,1:nat_org)=dtpawuj(jdtset)%occ(1,:)+dtpawuj(jdtset)%occ(2,:)
   else if (macro_uj==2.and.nspden==2) then
     chih(jdtset,1:nat_org)=dtpawuj(jdtset)%occ(1,:)
   else if (macro_uj==3.and.nspden==2) then
     chih(jdtset,1:nat_org)=dtpawuj(jdtset)%occ(2,:)
   end if
   vsh(jdtset)=dtpawuj(jdtset)%vsh(1,pawujat)
   if (pawprtvol==3) then
     write(message,fmt='(2a,i3,a,f15.12)') ch10,' Potential shift vsh(',jdtset,') =',vsh(jdtset)
     call wrtout(std_out,message,'COLL')
     write(message,fmt='( a,i3,a,120f15.9)') ' Occupations occ(',jdtset,') ',chih(jdtset,1:nat_org)
     call wrtout(std_out,message,'COLL')
   end if
 end do

 if (any(abs((/(vsh(ii)-vsh(ii+2), ii=1,2) /))<0.00000001)) then
   write(message, '(2a,18f10.7,a)' )  ch10,' vshift is too small: ',abs((/(vsh(ii)-vsh(ii+2), ii=1,2) /))
   call wrtout(ab_out,message,'COLL')
   return
 end if

!DEBUG
!write(message,fmt='(a)')'pawuj_det: after test vsh'
!call wrtout(std_out,message,'COLL')
!END DEBUG

 chi0_org=(chih(1,1:nat_org)-chih(3,1:nat_org))/(vsh(1)-vsh(3))/dtpawuj(1)%diemix
 chi_org=(chih(2,1:nat_org)-chih(4,1:nat_org))/(vsh(2)-vsh(4))

 if (pawprtvol==3) then
   write(message,fmt='(2a, 150f15.10)') ch10,' Chi_0n ',chi0_org
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a, 150f15.10)') ' Chi_n ',chi_org
   call wrtout(std_out,message,'COLL')
 end if

 write(message,fmt='(a)')': '
 if (nspden==2) then
   magv_org=dtpawuj(1)%occ(1,:)-dtpawuj(1)%occ(2,:)
   if (all(abs(magv_org)<0.001)) then
     magv_org=(/(1,im1=1,nat_org)/)
   else
     magv_org=abs(magv_org)/magv_org
     if (magv_org(1).lt.0) magv_org=magv_org*(-1_dp)
     if (all(magv_org(2:nat_org).lt.0)) then
       magv_org=abs(magv_org)
       write(message,'(a)')', (reset to fm): '
     end if
   end if
 else
   magv_org=(/(1,im1=1,nat_org)/)
 end if

 if (pawprtvol==3) then
   write(message,fmt='(3a, 150f4.1)') ch10,' Magnetization',trim(message),magv_org
   call wrtout(std_out,message,'COLL')
 end if

!Case of extrapolation to larger r_paw: calculate intg

 if (all(dtpawuj(1)%wfchr(:)/=0).and.ph0phiint/=1) then
   if (dtpawuj(1)%pawujrad<20.and.dtpawuj(1)%pawujrad>dtpawuj(1)%pawrad) then
     fcorr=(1-ph0phiint)/(IRadFnH(dtpawuj(1)%pawrad,20.0_dp,nint(dtpawuj(1)%wfchr(2)),&
&     nint(dtpawuj(1)%wfchr(3)),dtpawuj(1)%wfchr(1)))
     intg=ph0phiint/(1-fcorr*IRadFnH(dtpawuj(1)%pawujrad,20.0_dp,nint(dtpawuj(1)%wfchr(2)),&
&     nint(dtpawuj(1)%wfchr(3)),dtpawuj(1)%wfchr(1)))
     write(message, fmt='(a,f12.5,a,f12.5)') ' pawuj_det: met2 extrapolation to ', dtpawuj(1)%pawujrad,' using factor ',intg
     call wrtout(std_out,message,'COLL')
   else if (dtpawuj(1)%pawujrad<dtpawuj(1)%pawrad) then
     a=0 ; a(1:3)=dtpawuj(1)%wfchr(4:6)
     a=(/a(1)**2,a(1)*a(2),a(2)**2/3.0_dp+2.0_dp/3.0_dp*a(1)*a(3),a(2)*a(3)/2.0_dp,a(3)**2/5.0_dp/)
     b=(/(dtpawuj(1)%pawujrad**(im1)-dtpawuj(1)%pawrad**(im1),im1=1,5)/)
     intg=dot_product(a,b)
     intg=ph0phiint/(ph0phiint+intg)
     write(message, fmt='(a,f12.5,a,f12.5)') ' pawuj_det: met1 extrapolation to ', dtpawuj(1)%pawujrad,' using factor ',intg
     call wrtout(std_out,message,'COLL')
   else if (dtpawuj(1)%pawujrad==dtpawuj(1)%pawrad) then
     intg=one
     write(message, fmt='(a,2i7,3f12.5)') ' pawuj_det: no extrapolation (pawujrad=pawrad)'
     call wrtout(std_out,message,'COLL')
   else
     intg=ph0phiint
     write(message, fmt='(a,3f12.5)') ' pawuj_det: extrapolation to r->\inf using factor ',intg
     call wrtout(std_out,message,'COLL')
   end if
 else
   write(message, fmt='(a,2i7,3f12.5)') ' pawuj_det: no extrapolation (ph0phiint=1 or wfchr=0)'
   call wrtout(std_out,message,'COLL')
   intg=one
 end if


!Determine U in primitive cell

 write(message,fmt='(a)')' pawuj_det: determine U in primitive cell'
 call wrtout(std_out,message,'COLL')

 call lcalcu(int(magv_org),nat_org,dtpawuj(1)%rprimd,dtpawuj(1)%xred,chi_org,chi0_org,pawujat,ures,pawprtvol,pawujga,pawujoption)

!Begin calculate U in supercell

!Analize shell structure of primitive cell
!and atomic distances in primitive cell
 ABI_ALLOCATE(smult_org,(nat_org))
 ABI_ALLOCATE(sdistv_org,(nat_org))
 call shellstruct(dtpawuj(1)%xred,dtpawuj(1)%rprimd,nat_org,&
& int(magv_org),distv_org,smult_org,sdistv_org,nsh_org,pawujat,pawprtvol)

 ii=1
 write(message, fmt='(8a)') ' URES ','     ii','    nat','       r_max','    U(J)[eV]','   U_ASA[eV]','   U_inf[eV]',ch10
 write(message, fmt='(a,2i7,4f12.5)') trim(message)//' URES ',ii,nat_org,maxval(abs(distv_org)),ures,ures*exp(log(intg)*eyp),&
& ures*exp(log(ph0phiint)*eyp)
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 if (pawprtvol>1) then
   write(message,fmt='(a, 150f10.5)')' pawuj_det: ionic distances in original cell (distv_org) ', distv_org
   call wrtout(std_out,message,'COLL')
 end if

!Construct supercell, calculate limit dimensions of supercell
 ii=1
 maxnat=product(dtpawuj(1)%scdim(:))*nat_org
 if (maxnat==0) then
   maxnat=dtpawuj(1)%scdim(1)
   ext=(/ii, ii, ii/)
 else
   jj=1
   do while (jj<=3)
     ext(jj)=minval( (/ii, dtpawuj(1)%scdim(jj) /) )
     jj=jj+1
   end do
 end if
 ext=ext+(/ 1,1,1 /)
 ii=ii+1


 nat_sc=product(ext)*nat_org

!DEBUG
!write(message,fmt='(a,3i4)')'pawuj_det: ext ',ext
!call wrtout(std_out,message,'COLL')
!END DEBUG

 do while (nat_sc<=maxnat)
   ABI_ALLOCATE(chi0_sc,(nat_sc))
   ABI_ALLOCATE(chi_sc,(nat_sc))
   ABI_ALLOCATE(distv_sc,(nat_sc))
   ABI_ALLOCATE(magv_sc,(nat_sc))
   ABI_ALLOCATE(sdistv_sc,(nat_sc))
   ABI_ALLOCATE(smult_sc,(nat_sc))
   ABI_ALLOCATE(xred_sc,(3,nat_sc))

!  Determine positions=xred_sc and supercell dimensions=rpimd_sc
   call mksupercell(dtpawuj(1)%xred,int(magv_org),dtpawuj(1)%rprimd,nat_org,&
&   nat_sc,xred_sc,magv_sc,rprimd_sc,ext,pawprtvol)

!  Determine shell structure of supercell: multiplicities (smult_sc), radii of shells (sdistv_sc)
!  number of shells (nsh_sc) and atomic distances in supercell (distv_sc)

   write(message,fmt='(a,3i3,a)')' pawuj_det: determine shell structure of ',ext(1:3),' supercell'
   call wrtout(std_out,message,'COLL')

   call shellstruct(xred_sc,rprimd_sc,nat_sc,int(magv_sc),distv_sc,smult_sc,sdistv_sc,nsh_sc,pawujat,pawprtvol)

   if (pawprtvol>=2) then
     write(message,fmt='(a)')' pawuj_det: ionic distances in supercell (distv_sc) '
     call wrtout(std_out,message,'COLL')
     call prmat(distv_sc,1,nat_sc,1,std_out)
   end if

!  Determine chi and chi0 in supercell (chi0_sc, chi_sc)
!  DEBUG
!  write(message,fmt='(a)')'pawuj_det:  chi and chi0 in supercell'
!  call wrtout(std_out,message,'COLL')
!  END DEBUG

   if (pawujoption>2) then
     invopt=2
   else
     invopt=1
   end if

   call chiscwrt(chi0_org,distv_org,nat_org,sdistv_org,smult_org,nsh_org,&
&   chi0_sc,distv_sc,nat_sc,smult_sc,nsh_sc,invopt,pawprtvol)
   call chiscwrt(chi_org,distv_org,nat_org,sdistv_org,smult_org,nsh_org,&
&   chi_sc,distv_sc,nat_sc,smult_sc,nsh_sc,invopt,pawprtvol)

!  Calculate U in supercell
!  DEBUG
!  write(message,fmt='(a)')'pawuj_det:   U in supercell'
!  call wrtout(std_out,message,'COLL')
!  END DEBUG
   call lcalcu(int(magv_sc),nat_sc,rprimd_sc,xred_sc,chi_sc,chi0_sc,pawujat,ures,pawprtvol,pawujga,pawujoption)

   write(message, fmt='(a,2i7,4f12.5)') ' URES ',ii,nat_sc,maxval(abs(distv_sc)),ures,ures*exp(log(intg)*eyp),&
&   ures*exp(log(ph0phiint)*eyp)
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   ABI_DEALLOCATE(chi0_sc)
   ABI_DEALLOCATE(chi_sc)
   ABI_DEALLOCATE(distv_sc)
   ABI_DEALLOCATE(magv_sc)
   ABI_DEALLOCATE(sdistv_sc)
   ABI_DEALLOCATE(smult_sc)
   ABI_DEALLOCATE(xred_sc)

   if (product(dtpawuj(1)%scdim(:))==0) then
     ext=(/ii, ii, ii/)
   else
     jj=1
     do while (jj<=3)
       ext(jj)=minval( (/ii, dtpawuj(1)%scdim(jj) /) )
       jj=jj+1
     end do
   end if
   ext=ext+(/ 1,1,1 /)
   ii=ii+1

   nat_sc=product(ext)*nat_org

 end do

 write(message, '(3a)' )ch10,' ---------- calculate U, (J) end -------------- ',ch10
 call wrtout(ab_out,message,'COLL')

 ABI_DEALLOCATE(dparrr)
 ABI_DEALLOCATE(dqarr)
 ABI_DEALLOCATE(dqarrr)
 ABI_DEALLOCATE(jdtset_)
 ABI_DEALLOCATE(chi_org)
 ABI_DEALLOCATE(chi0_org)
 ABI_DEALLOCATE(smult_org)
 ABI_DEALLOCATE(sdistv_org)
 ABI_DEALLOCATE(chih)
 ABI_DEALLOCATE(idum2)
 ABI_DEALLOCATE(drarr)
 ABI_DEALLOCATE(dparr)
 ABI_DEALLOCATE(magv_org)
 ABI_DEALLOCATE(xred_org)
 ABI_DEALLOCATE(distv_org)

 DBG_EXIT("COLL")

end subroutine pawuj_det
!!***

!----------------------------------------------------------------------
!!****f* m_paw_uj/pawuj_red
!! NAME
!!  pawuj_red
!!
!! FUNCTION
!!  Store atomic occupancies, potential shift, positions in dtpawuj datastructure.
!!
!! INPUTS
!!  fatvshift=factor that multiplies atvshift
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell
!!  ntypat = number of atom types
!!  paw_ij(my_natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawprtvol= printing volume
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  dtpawuj(0:ndtpawuj) (initialization of fields vsh, occ, iuj,nnat)
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,linvmat,prmat,wrtout,xmpi_lor,xmpi_sum
!!
!! SOURCE

subroutine pawuj_red(dtset,dtpawuj,fatvshift,my_natom,natom,ntypat,paw_ij,pawrad,pawtab,ndtpawuj,&
&                    mpi_atmtab,comm_atom) ! optional arguments (parallelism)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)                 :: my_natom,natom,ntypat,ndtpawuj
 integer,optional,intent(in)        :: comm_atom
 real(dp),intent(in)                :: fatvshift
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(paw_ij_type),intent(in)       :: paw_ij(my_natom)
 type(pawtab_type),intent(in)       :: pawtab(ntypat)
 type(pawrad_type),intent(in)       :: pawrad(ntypat)
 type(dataset_type),intent(in)      :: dtset
 type(macro_uj_type),intent(inout)  :: dtpawuj(0:ndtpawuj)

!Local variables-------------------------------
!scalars
 integer,parameter           :: natmax=2,ncoeff=3
 integer                     :: iatom,iatom_tot,ierr,im1,im2,ispden,itypat,ll,nspden,nsppol,iuj
 integer                     :: my_comm_atom,nnat,natpawu,natvshift,pawujat,ndtset,typawujat
 logical                     :: usepawu !antiferro,
 logical                     :: my_atmtab_allocated,paral_atom
 character(len=1000)         :: message,hstr
 character(len=500)          :: messg
!arrays
 logical                     :: musk(3,natom)
 integer,pointer             :: my_atmtab(:)
 real(dp)                    :: rrtab(ncoeff),wftab(ncoeff),a(ncoeff,ncoeff),b(ncoeff,ncoeff)! ,a(ncoeff,ncoeff)
 real(dp),allocatable        :: nnocctot(:,:) !,magv(:)
 real(dp),allocatable        :: atvshift(:,:,:) ! atvshift(natvshift,2,natom)
 logical,allocatable         :: dmusk(:,:),atvshmusk(:,:,:) !atvshmusk(natvshift,2,natom)

! *********************************************************************

!Initializations
 nspden=1;nsppol=1
 if (my_natom>0) then
   nspden=paw_ij(1)%nspden ; nsppol=paw_ij(1)%nsppol
 end if

 natvshift=dtset%natvshift
 pawujat=dtset%pawujat
 natpawu=dtset%natpawu   ; ndtset=dtset%ndtset
 ABI_ALLOCATE(atvshift,(natvshift,nspden,natom))
 ABI_ALLOCATE(atvshmusk,(natvshift,nspden,natom))
 ABI_ALLOCATE(dmusk,(nspden,natom))
 ABI_ALLOCATE(nnocctot,(nspden,natom))
 musk=.false.; dmusk=.false.
 atvshift=fatvshift*dtset%atvshift
 atvshmusk=.false.
 nnocctot=zero
 typawujat=dtset%typat(pawujat)
 usepawu=(count(pawtab(:)%usepawu/=0)>0)

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 nnat=0
 if (usepawu) then
   write(message,'(3a)') ch10, '---------- pawuj_red ------ ',ch10
   call wrtout(std_out,  message,'COLL');
   do iatom=1,my_natom
     iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
     itypat=paw_ij(iatom)%itypat;ll=pawtab(itypat)%lpawu
     if ((ll>=0).and.(pawtab(itypat)%usepawu/=0).and.itypat==typawujat) then
       musk(:,iatom_tot)=(/.true., .true., .true. /)
       atvshmusk(:,:,iatom_tot)=reshape((/ (( (im1==1), im1=1,natvshift)  ,im2=1,nspden ) /),(/natvshift,nspden/))
       do ispden=1,nspden
         nnocctot(ispden,iatom_tot)=paw_ij(iatom)%nocctot(ispden)
         dmusk(ispden,iatom_tot)=.true.
       end do
       nnat=nnat+1
     end if
   end do

!  Reduction in case of parallelism
   if (paral_atom) then
     call xmpi_sum(nnocctot ,my_comm_atom,ierr)
     call xmpi_lor(atvshmusk,my_comm_atom) ! dim=natpawu ???
     call xmpi_lor(dmusk    ,my_comm_atom)
     call xmpi_lor(musk     ,my_comm_atom)
     call xmpi_sum(nnat     ,my_comm_atom,ierr)
   end if

   iuj=maxval(dtpawuj(:)%iuj)
   !write(std_out,*)' pawuj_red: iuj',iuj
   !write(std_out,*)' pawuj_red: dtpawuj(:)%iuj ',dtpawuj(:)%iuj

   if (iuj==1.or.iuj==3) then  ! 1 and 3: non-scf steps
     dtpawuj(iuj+1)%iuj=iuj+1
   end if

   ! TODO: check that this is correct: this point is passed several times
   ! for a given value of iuj - should the stuff be accumulated instead of replaced?
   if(allocated(dtpawuj(iuj)%vsh))  then
     ABI_DEALLOCATE(dtpawuj(iuj)%vsh)
   end if
   ABI_ALLOCATE(dtpawuj(iuj)%vsh,(nspden,nnat))
   if(allocated(dtpawuj(iuj)%occ))  then
     ABI_DEALLOCATE(dtpawuj(iuj)%occ)
   end if
   ABI_ALLOCATE(dtpawuj(iuj)%occ,(nspden,nnat))
   if(allocated(dtpawuj(iuj)%xred))  then
     ABI_DEALLOCATE(dtpawuj(iuj)%xred)
   end if
   ABI_ALLOCATE(dtpawuj(iuj)%xred,(3,nnat))

   if (iuj==1) then
     ABI_ALLOCATE(dtpawuj(0)%vsh,(nspden,nnat))
     ABI_ALLOCATE(dtpawuj(0)%occ,(nspden,nnat))
     ABI_ALLOCATE(dtpawuj(0)%xred,(3,nnat))
     dtpawuj(0)%vsh=0
     dtpawuj(0)%occ=0
     dtpawuj(0)%xred=0
   end if

   rrtab=(/0.75_dp,0.815_dp,1.0_dp/)*pawtab(typawujat)%rpaw
   wftab=pawtab(typawujat)%phi(pawtab(typawujat)%mesh_size,pawtab(typawujat)%lnproju(1))

!  DEBUG
!  write(std_out,*)' pawuj_red: rrtab ',rrtab
!  write(std_out,*)' pawuj_red: wftab ',wftab
!  END DEBUG

   do im1=1,ncoeff
!    DEBUG
!    write(std_out,*)' pawuj_red: ncoeff ',ncoeff,' im1 ',im1
!    END DEBUG
     if (pawrad(typawujat)%mesh_type==1) then
       im2=nint(rrtab(im1)/pawrad(typawujat)%rstep+1)
     else if (pawrad(typawujat)%mesh_type==2) then
       im2=nint(log(rrtab(im1)/pawrad(typawujat)%rstep+1)/pawrad(typawujat)%lstep+1)
     else if (pawrad(typawujat)%mesh_type==3) then
       im2=nint(log(rrtab(im1)/pawrad(typawujat)%rstep)/pawrad(typawujat)%lstep+1)
     else if (pawrad(typawujat)%mesh_type==4) then
       im2=nint(pawtab(typawujat)%mesh_size*(1-exp((-one)*rrtab(im1)/pawrad(typawujat)%rstep))+1)
     end if

!    DEBUG
!    write(std_out,*)' pawuj_red: im2 ',im2
!    END DEBUG

     rrtab(im1)=pawrad(typawujat)%rad(im2)
     wftab(im1)=pawtab(typawujat)%phi(im2,pawtab(typawujat)%lnproju(1))
   end do
   write(message,fmt='(a,i3,a,10f10.5)')' pawuj_red: mesh_type',pawrad(typawujat)%mesh_type,' rrtab:', rrtab
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a,10f10.5)')' pawuj_red: wftab', wftab
   call wrtout(std_out,message,'COLL')
   a=reshape((/ (( rrtab(im2)**(im1-1), im1=1,3)  ,im2=1,3 )/),(/ncoeff,ncoeff/))
   write(messg,fmt='(a)')'A'
   call linvmat(a,b,ncoeff,messg,2,0.0_dp,3) ! linvmat(inmat,oumat,nat,nam,option,gam,prtvol)
   write(std_out,*) 'pawuj_red: a,b ', a,b
   wftab=matmul(wftab,b)
   write(std_out,*) 'pawuj_red: wftab ', wftab
   dtpawuj(iuj)%wfchr(4:6)=wftab

   dtpawuj(iuj)%nat=nnat
   write(std_out,*) 'pawuj_red: m1'
   dtpawuj(iuj)%vsh=reshape(pack(atvshift,atvshmusk),(/ nspden,nnat /))
!  factor in next line to compensate nocctot contains just occ of 1 spin channel for nspden=1
   write(std_out,*) 'pawuj_red: m2'
   dtpawuj(iuj)%occ=reshape(pack(nnocctot,dmusk),(/nspden,nnat/))*(3-nspden)
   write(std_out,*) 'pawuj_red: m3'
!  dtpawuj(iuj)%occ=dtpawuj(iuj)%occ/pawtab(typawujat)%ph0phiint(1)

   write(std_out,*) 'pawuj_red: occ ', dtpawuj(iuj)%occ

   dtpawuj(iuj)%xred=reshape(pack(dtset%xred_orig(:,:,1),musk),(/3,nnat/))
   dtpawuj(iuj)%ph0phiint=pawtab(typawujat)%ph0phiint(1)
   dtpawuj(iuj)%wfchr(1:3)=(/ pawtab(typawujat)%zioneff(1)*(dtset%lpawu(typawujat)+2),&
&   one*(dtset%lpawu(typawujat)+1),one*(dtset%lpawu(typawujat))/)
   dtpawuj(iuj)%pawrad=pawtab(typawujat)%rpaw

   write(std_out,*) 'pawuj_red: wfchr ',dtpawuj(iuj)%wfchr


   write (hstr,'(I0)') iuj
   write(message,'(a,a,I3,I3,a)') ch10, '---------- MARK ------ ',iuj,maxval(dtpawuj(:)%iuj) ,ch10
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a)') 'vsh'//trim(hstr)
   call wrtout(std_out,message,'COLL')
   call prmat(dtpawuj(iuj)%vsh(:,:),1,nnat*nspden,1)
   write(message,fmt='(a)') 'occ'//trim(hstr)
   call wrtout(std_out,message,'COLL')
   call prmat(dtpawuj(iuj)%occ(:,:),1,nnat*nspden,1)
   write(message, '(3a)' )'---------- MARK ---------- ',ch10
   call wrtout(std_out,message,'COLL')
 end if !usepawu

 ABI_DEALLOCATE(nnocctot)
 ABI_DEALLOCATE(dmusk)
 ABI_DEALLOCATE(atvshift)
 ABI_DEALLOCATE(atvshmusk)

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine pawuj_red
!!***

!----------------------------------------------------------------------

!!****f* m_paw_uj/chiscwrt
!! NAME
!!  chiscwrt
!!
!! FUNCTION
!!  PRIVATE ROUTINE
!!  Distributes values of chi_org on chi_sc according to ion-ion distances and multiplicities in shells
!!
!! INPUTS
!!  chi_org  = response matrix in original unit cell
!!  disv_org = distances (multiplied by magntization) in original unit cell
!!  nat_org = number of atoms in original unit cell
!!  sdisv_org = radii of shells in original unit cell
!!  smult_org = multiplicities of shells in original unit cell
!!  nsh_org   = number of shells in original unit cell
!!  disv_sc   = distances (multiplied by magntization) in super-cell
!!  nat_sc  = number of atoms in super-cell
!!  sdisv_sc = radii of shells in super-cell (was unused, so suppressed)
!!  smult_sc = multiplicities of shells in super-cell
!!  nsh_sc = number of shells in super-cell
!!  opt =
!!
!! OUTPUT
!!  chi_sc = first line of reponse matrix in super-cell
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pawuj_det
!!
!! CHILDREN
!!      prmat,wrtout
!!
!! SOURCE

subroutine chiscwrt(chi_org,disv_org,nat_org,sdisv_org,smult_org,nsh_org,chi_sc,&
& disv_sc,nat_sc,smult_sc,nsh_sc,opt,prtvol)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)              :: nat_org,nsh_org,nsh_sc
 integer,intent(in)              :: nat_sc
 integer,intent(in),optional     :: prtvol
 integer,intent(in),optional     :: opt
!arrays
 real(dp),intent(in)             :: chi_org(nat_org),disv_org(nat_org),sdisv_org(nsh_org)
 integer,intent(in)              :: smult_org(nsh_org), smult_sc(nsh_sc)
 real(dp),intent(in)             :: disv_sc(nat_sc)
 real(dp),intent(out)            :: chi_sc(nat_sc)

!Local variables-------------------------------
!scalars
 integer                      :: iatom,jatom,jsh,optt
 character(len=500)           :: message
!arrays
 real(dp)                     :: chi_orgl(nat_org)

! *************************************************************************

 if (present(opt)) then
   optt=opt
 else
   optt=1
 end if


 do iatom=1,nat_org
   do jsh=1,nsh_org
     if (disv_org(iatom)==sdisv_org(jsh)) then
       if (opt==2) then
         chi_orgl(iatom)=chi_org(iatom)
       else if (opt==1) then
         chi_orgl(iatom)=chi_org(iatom)*smult_org(jsh)/smult_sc(jsh)
       end if
       exit
     end if
   end do !iatom
 end do  !jsh

 if (prtvol>1) then
   write(message,fmt='(a,150f10.5)')' chiscwrt: chi at input ',chi_org
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a,150f10.5)')' chiscwrt: chi after division ',chi_orgl
   call wrtout(std_out,message,'COLL')
 end if

 do iatom=1,nat_sc
   do jatom=1,nat_org
     if (disv_org(jatom)==disv_sc(iatom)) then
       chi_sc(iatom)=chi_orgl(jatom)
       exit
     else if (jatom==nat_org) then
       chi_sc(iatom)=0_dp
     end if
   end do
 end do

 if (prtvol>1) then
   write(message,'(a)')' chiscwrt, chi_sc '
   call wrtout(std_out,message,'COLL')
   call prmat(chi_sc,1,nat_sc,1,std_out)
 end if

end subroutine chiscwrt
!!***

!----------------------------------------------------------------------

!!****f* m_paw_uj/linvmat
!! NAME
!!  linvmat
!!
!! FUNCTION
!!  inverts real matrix inmat
!!
!! INPUTS
!!  inmat(1:nat,1:nat)=matrix to be inverted
!!  nat=dimension of inmat
!!  nam=comment specifiying the input matrix (to be printed in output)
!!  option=how to invert inmat
!!      option=1 or 3 add charge bath to matrix and add gam for inversion
!!      option=2 simply invert matrix
!!  gam=gamma added to inmat before inversion in case charge bath is used (allows inversion of otherwise
!!               singular matrix)
!!  prtvol=controls output to files (see subroutine lprtmat)
!!
!! OUTPUT
!!  oumat(nnat,nnat)=inverse of inmat, nnat=nat+1 for option=1 or option=3; nnat=nat for option=2
!!
!! PARENTS
!!      pawuj_red,pawuj_utils
!!
!! CHILDREN
!!
!! SOURCE

subroutine linvmat(inmat,oumat,nat,nam,option,gam,prtvol)

 implicit none

!Arguments -------------------------------

 integer,intent(in)              :: nat
 real(dp),intent(in)             :: gam,inmat(nat,nat)
 real(dp),intent(inout)          :: oumat(:,:)         ! nat+1,nat+1 for option=1 or 2
                                                       ! nat,nat for option=2
 character(len=500),intent(in)   :: nam
 integer,intent(in),optional     :: prtvol,option

!Local variables -------------------------
 character(len=500)             :: message
 character(len=500)             :: bastrin,gastrin
 integer                        :: info,nnat,optionn
 integer,allocatable            :: ipvt(:)
 real(dp),allocatable           :: hma(:,:),work(:)

! *********************************************************************

 if (present(option)) then
   optionn=option
 else
   optionn=1
 end if

 if (option==1.or.option==3) then
   write(bastrin,'(a)')'+ charge bath '
   write(gastrin,'(a,d10.2,a)')'+ gamma  (=',gam,') '
 else
   write(bastrin,'(a)')''
   write(gastrin,'(a)')''
 end if

 write(message,fmt='(a)')' matrix '//trim(nam)
 call lprtmat(message,1,prtvol,inmat,nat)

 if (option==1.or.option==3) then
   call blow_pawuj(inmat,nat,oumat)
   write(message,fmt='(a,a)')' ',trim(nam)//trim(bastrin)
   call lprtmat(message,1,prtvol,oumat,nat+1)
   oumat=oumat+gam
   write(message,fmt='(a,a)')' ',trim(nam)//trim(bastrin) ! //trim(gastrin)
   call lprtmat(message,1,prtvol,oumat-gam,nat+1)
   nnat=nat+1
 else
   nnat=nat
   oumat=inmat
   oumat(1,1)=inmat(1,1)
 end if


 ABI_ALLOCATE(hma,(nnat,nnat))
 ABI_ALLOCATE(work,(nnat))
 ABI_ALLOCATE(ipvt,(nnat))
 work=0_dp
 hma(:,:)=oumat

 call dgetrf(nnat,nnat,hma,nnat,ipvt,info)
 if (.not.info==0) then
   write(message, '(3a)' ) 'Matrix '//trim(nam)//' is singular',ch10,'Probably too many symmetries kept'
   call wrtout(ab_out,message,'COLL')
   return
 end if

 call dgetri(nnat,hma,nnat,ipvt,work,nnat,info)
 oumat=hma(:,:)

 write(message,fmt='(2a,a)')' ('//trim(nam)//trim(bastrin)//trim(gastrin)//')^(-1)'
 call lprtmat(message,1,prtvol,oumat,nnat)

 ABI_DEALLOCATE(hma)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(ipvt)

end subroutine linvmat
!!***

!----------------------------------------------------------------------

!!****f* m_paw_uj/lprtmat
!! NAME
!!  lprtmat
!!
!! FUNCTION
!!  prints out the real matrix mmat
!!
!! INPUTS
!!  mmat(nat,nat)=matrix to be printed
!!  nat=dimension of mmat
!!  prtvol specifies the volume of printing
!!   3: print the whole matrix
!!   2: print the first line
!!   1: do not print anything
!!  chan specifies the output files
!!   1: output only to std_out
!!   2: output also to ab_out
!!  commnt=comment specifying matirix
!!
!! OUTPUT
!!  oumat(nat+1,nat+1)=inverse of inmat
!!
!! PARENTS
!!      pawuj_utils
!!
!! CHILDREN
!!
!! SOURCE

subroutine lprtmat(commnt,chan,prtvol,mmat,nat)

 implicit none

!Arguments -------------------------------
 integer,intent(in)              :: nat,chan,prtvol
 real(dp),intent(in)             :: mmat(nat,nat)
 character(len=500),intent(in)  :: commnt

!Local variables -------------------------
 character(len=500)             :: message
! *********************************************************************

 if (prtvol==3) then
   write(message,fmt='(a)') trim(commnt)
   call wrtout(std_out,message,'COLL')
   call prmat(mmat,nat,nat,nat,std_out)
   if (chan==2) then
     call wrtout(ab_out,message,'COLL')
     call prmat(mmat,nat,nat,nat,ab_out)
   end if
   write(message,*)ch10
   call wrtout(std_out,message,'COLL')
   if (chan==2) then
     call wrtout(ab_out,message,'COLL')
   end if
 end if

 if (prtvol==2) then
   write(message,fmt='(a)') trim(commnt)
   call wrtout(std_out,message,'COLL')
   call prmat(mmat,1,nat,nat,std_out)
   if (chan==2) then
     call wrtout(ab_out,message,'COLL')
     call prmat(mmat,1,nat,nat,ab_out)
   end if
   write(message,*)ch10
   call wrtout(std_out,message,'COLL')
   if (chan==2) then
     call wrtout(ab_out,message,'COLL')
   end if
 end if

end subroutine lprtmat
!!***

!----------------------------------------------------------------------

!!****f* m_paw_uj/lcalcu
!! NAME
!!  lcalcu
!!
!! FUNCTION
!!  prints out real the real matrice mmat
!!
!! INPUTS
!!  magv=magnetic ordering of the ions (-1 of down/1 for up)
!!  natom=number of atoms
!!  rprimd(3,3)=lattic vectors of unit cell
!!  xred(3,natom)=positions of atoms
!!  chi(natom)=full response of atoms due to shift on atom pawujat
!!  chi0(natom)= response of atoms due to shift on atom pawujat
!!  pawujat=specifies on which atom the potential shift was done
!!  prtvol=controls output to files (see subroutine lprtmat)
!!  gam=gamma to be used for inversion of matrices (see subroutine livmat)
!!  opt=wether to use charge bath (1 or 3) or not (else)
!!
!! OUTPUT
!!  ures=resulting U (in eV) on atom pawujat
!!
!! PARENTS
!!      pawuj_det
!!
!! CHILDREN
!!
!! SOURCE

subroutine lcalcu(magv,natom,rprimd,xred,chi,chi0,pawujat,ures,prtvol,gam,opt)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)          :: natom
 integer,intent(in),optional :: opt,pawujat,prtvol
 real(dp),intent(in),optional:: gam
 real(dp),intent(out)        :: ures
!arrays
 integer,intent(in)          :: magv(natom)
 real(dp),intent(in)         :: rprimd(3,3),xred(3,natom),chi(natom),chi0(natom)

!Local variables-------------------------------
!scalars
 character(len=500)          :: message
 integer                     :: optt,prtvoll,nnatom
 real(dp)                    :: gamm
!arrays
 real(dp),allocatable        :: tab(:,:,:)

! *********************************************************************

 if (present(opt)) then
   optt=opt
 else
   optt=1
 end if

 if (present(prtvol)) then
   prtvoll=prtvol
 else
   prtvoll=1
 end if

 if (present(gam)) then
   gamm=gam
 else
   gamm=1_dp
 end if

 if (optt==1.or.optt==3) then
   nnatom=natom+1
 else
   nnatom=natom
 end if

 ABI_ALLOCATE(tab,(4,nnatom,nnatom))

 call ioniondist(natom,rprimd,xred,tab(1,1:natom,1:natom),3,chi0,magv,pawujat,prtvoll)
 call ioniondist(natom,rprimd,xred,tab(2,1:natom,1:natom),3,chi,magv,pawujat,prtvoll)


 write(message,fmt='(a)')'response chi_0'
 call linvmat(tab(1,1:natom,1:natom),tab(3,1:nnatom,1:nnatom),natom,message,optt,gamm,prtvoll)

 write(message,fmt='(a)')'response chi'
 call linvmat(tab(2,1:natom,1:natom),tab(4,1:nnatom,1:nnatom),natom,message,optt,gamm,prtvoll)

 tab(1,1:nnatom,1:nnatom)=(tab(3,1:nnatom,1:nnatom)-tab(4,1:nnatom,1:nnatom))*Ha_eV

 write(message,fmt='(a,i3,a)')' (chi_0)^(-1)-(chi)^(-1) (eV)'
 call lprtmat(message,2,prtvoll,tab(1,1:nnatom,1:nnatom),nnatom)

 ures=tab(1,1,pawujat)

 ABI_DEALLOCATE(tab)

end subroutine lcalcu
!!***

!----------------------------------------------------------------------

!!****f* m_paw_uj/blow_pawuj
!!
!! NAME
!! blow_pawuj
!!
!! FUNCTION
!! This subroutine reads a real nxn matrice and appends lines n+1 and clumn n+1 containing
!! the sum of the lines
!!
!! INPUTS
!!  mat(nj,nj) matrix to be completed
!!
!! OUTPUT
!!  matt(nj+1,nj+1) completed matrix
!!
!! PARENTS
!!      pawuj_utils
!!
!! CHILDREN
!!
!! SOURCE

subroutine blow_pawuj(mat,nj,matt)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)        :: nj
!arrays
 real(dp),intent(in)       :: mat(nj,nj)
 real(dp),intent(out)      :: matt(nj+1,nj+1)

!Local variables-------------------------------
!scalars
 integer                   :: ii
!arrays

! *************************************************************************

 matt(1:nj,1:nj)=mat
 do  ii = 1,nj
   matt(ii,nj+1)=-sum(matt(ii,1:nj))
 end do

 do  ii = 1,nj+1
   matt(nj+1,ii)=-sum(matt(1:nj,ii))
 end do

end subroutine blow_pawuj
!!***

!----------------------------------------------------------------------

END MODULE m_paw_uj
!!***
