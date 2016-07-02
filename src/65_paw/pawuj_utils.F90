!{\src2tex{textfont=tt}}
!!****f* ABINIT/linvmat 
!! NAME
!!  linvmat
!!
!! FUNCTION
!!  inverts real matrix inmat 
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DJA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine linvmat(inmat,oumat,nat,nam,option,gam,prtvol)

 use defs_basis
 use m_profiling_abi
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'linvmat'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_65_paw, except_this_one => linvmat
!End of the abilint section

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

!!****f* ABINIT/lprtmat
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

 use m_profiling_abi
 use defs_basis

 use m_pptools,   only : prmat

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lprtmat'
 use interfaces_14_hidewrite
!End of the abilint section

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

!!****f* ABINIT/lcalcu 
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

 use defs_basis
 use defs_parameters
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lcalcu'
 use interfaces_41_geometry
 use interfaces_65_paw, except_this_one => lcalcu
!End of the abilint section

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

!!****f* ABINIT/pawuj_free 
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

 use m_profiling_abi
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawuj_free'
!End of the abilint section

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

