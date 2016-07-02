!{\src2tex{textfont=tt}}
!!****f* ABINIT/shellstruct 
!! NAME
!!  shellstruct 
!!
!! FUNCTION calculates shell structure (multiplicities, radii) of an atomic configuration 
!!  
!! INPUTS
!!  natom=number of atoms in unit cell
!!  xred=reduced coordinates of atoms
!!  rprimd=unit cell vectors
!!  magv = magnetic ordering of atoms given as 1 and -1, if not given fm is assumed
!!  atp = atom on which the perturbation was done
!!
!! OUTPUT
!!  sdisv(nat)= distance of each shell to central atom (only the first nsh entries are relevant)
!!  nsh    = number of shells
!!  mult(nat) = number of atoms on shell      (only the first nsh entries are relevant)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pawuj_det
!!
!! CHILDREN
!!      ioniondist,prmat,sort_dp,sort_int,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine shellstruct(xred,rprimd,natom,magv,distv,smult,sdisv,nsh,atp,prtvol) 

 use defs_basis
 use m_profiling_abi

 use m_pptools,       only : prmat

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'shellstruct'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_41_geometry, except_this_one => shellstruct
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)              :: natom
 integer,intent(in),optional     :: atp  
 integer,intent(in),optional     :: prtvol
 integer,intent(out)             :: nsh          
!arrays
 real(dp),intent(in)             :: rprimd(3,3)
 real(dp),intent(in)             :: xred(3,natom)
 integer,intent(out)             :: smult(natom) 
 integer,intent(in),optional     :: magv(natom)
 real(dp),intent(out)            :: sdisv(natom)
 real(dp),intent(out)            :: distv(natom)

!Local variables-------------------------------
!scalars
 integer                      :: iatom,atpp,ish,prtvoll
 character(len=500)           :: message
 real(dp),parameter           :: rndfact=10000_dp
!arrays
 integer                      :: iperm(natom),jperm(natom)
 real(dp)                     :: distvh(natom,natom)
 real(dp)                     :: magvv(natom)

! *************************************************************************

 if (present(magv)) then
   magvv=magv
 else
   magvv=(/ (1, iatom=1,natom)  /) 
 end if

 if (present(atp)) then
   atpp=atp 
 else
   atpp=1 
 end if

 if (present(prtvol)) then
   prtvoll=prtvol 
 else
   prtvoll=1
 end if

!DEBUB
 write(std_out,*)'shellstruct start'
!END DEBUG 

!Calculate ionic distances  
 call ioniondist(natom,rprimd,xred,distvh,1,magv=int(magvv),atp=atpp)
 distv=distvh(1,:)

 if (prtvol>2) then
   write(std_out,'(a)')' shellstruct ionic distances in cell (distv) : '
   call prmat(distv(1:natom),1,natom,1,std_out)
 end if

 iperm=(/ (iatom, iatom=1,natom ) /)
 jperm=iperm
 distv=anint(distv*rndfact)/rndfact
!Sort distances
 call sort_dp(natom,distv,iperm,10d-5)
 call sort_int(natom,iperm,jperm)
 
 smult=0
 sdisv=dot_product(rprimd(1,:),rprimd(1,:))+dot_product(rprimd(2,:),rprimd(2,:))+dot_product(rprimd(3,:),rprimd(3,:))

 nsh=1
 smult(1)=1
 sdisv(1)=distv(1)

 do iatom=2,natom
   do ish=1,natom
     if (distv(iatom)>sdisv(ish)) then
       cycle  
     else if (distv(iatom)==sdisv(ish)) then
       smult(ish)=smult(ish)+1
       exit    
     else if (distv(iatom)<sdisv(ish)) then
       smult(ish+1:natom)=smult(ish:natom-1)
       sdisv(ish+1:natom)=sdisv(ish:natom-1)
       smult(ish)=1
       sdisv(ish)=distv(iatom)
       nsh=nsh+1
       exit
     end if
   end do  
 end do

 distv=(/ ( distv(jperm(iatom)),iatom=1,natom ) /)
 
 if (prtvoll>2) then
   write(message,'(a,i4,a)')' shellstruct found ',nsh,' shells at distances (sdisv) '
   call wrtout(std_out,message,'COLL')
   call prmat(sdisv(1:nsh),1,nsh,1,std_out)  
   write(message,fmt='(a,150i4)')' and multiplicities (smult) ', smult(1:nsh)
   call wrtout(std_out,message,'COLL')
 end if
 
!DEBUB
 write(std_out,*)'shellstruct leave'
!END DEBUG

end subroutine shellstruct 
!!***
