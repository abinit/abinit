!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawuj_det
!! NAME
!!  pawuj_det
!!
!! FUNCTION
!!  From the complete dtpawuj-dataset determines U (or J) parameter for 
!!  PAW+U calculations
!!  Relevant only for automatic determination of U in PAW+U context
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DJA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawuj_det(dtpawuj,ndtpawuj,ujdet_filename,ures)

 use defs_basis
 use defs_abitypes
 use defs_parameters
 use m_profiling_abi
 use m_errors

 use m_special_funcs,  only : iradfnh
 use m_pptools,        only : prmat

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawuj_det'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_57_iovars
 use interfaces_65_paw, except_this_one => pawuj_det
!End of the abilint section

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
