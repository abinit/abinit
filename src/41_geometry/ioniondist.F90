!{\src2tex{textfont=tt}}
!!****f* ABINIT/ioniondist
!! NAME
!! ioniondist
!!
!! FUNCTION
!!  Compute ion-ion distances
!!
!! INPUTS
!!  natom= number of atoms in unit cell
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  xred(3,natom)=dimensionless reduced coordinates of atoms
!!  inm(natom,natom)=index (m,n) of the atom
!!  option= 1 output ion-ion distances / 2 output ordering of ion-ion
!!          distances / 3 output variables in varlist
!!          according to ion-ion distances * magnetic ordering
!!          magv magnetic ordering of atoms given als 1 and -1, if not
!!          given fm is assumed
!!  varlist=List of variables
!!  magv(natom)= magnetic ordering of atoms
!!  atp=atom on which the perturbation was done
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pawuj_utils,shellstruct
!!
!! CHILDREN
!!      prmat,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ioniondist(natom,rprimd,xred,inm,option,varlist,magv,atp,prtvol)

 use defs_basis
 use m_profiling_abi

 use m_pptools,    only : prmat

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ioniondist'
 use interfaces_14_hidewrite
 use interfaces_41_geometry, except_this_one => ioniondist
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)              :: natom,option
 integer,intent(in),optional     :: atp                   !atom on which the perturbation was done
!arrays
 real(dp),intent(in)             :: rprimd(3,3)
 real(dp),intent(in)             :: xred(3,natom)
 real(dp),intent(out)            :: inm(natom,natom)      
 integer,intent(in),optional     :: magv(natom)
 real(dp),intent(in),optional    :: varlist(natom)
 integer,intent(in),optional     :: prtvol

!Local variables-------------------------------
!scalars
 integer                      :: iatom,jatom,katom,kdum,atpp,prtvoll
 !character(len=500)           :: message
!arrays
 integer                      :: interq(natom)
 real(dp)                     :: hxcart(3,natom),distm(natom,natom)
 real(dp)                     :: magvv(natom)

! *************************************************************************

 hxcart=matmul(rprimd,xred)
 interq=(/(iatom,iatom=1,natom)/)
 inm=0

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

 if (option==3.and.(.not.present(varlist))) then
   call  wrtout(std_out,'ioniondist error: option=3 but no variable list provided for symmetrization','COLL')
   return
 end if


!DEBUG
!write(message, '(a,a)' ) ch10,' ioniondist start '
!call wrtout(std_out,message,'COLL')
!END DEBUG

 distm=0
 katom=atpp-1
 do iatom=1,natom
   katom=katom+1
   if (katom > natom) katom=1
   distm(iatom,iatom)=0
   do jatom=iatom,natom
     distm(iatom,jatom)=dist2(xred(:,katom),xred(:,jatom),rprimd,1)*magvv(katom)*magvv(jatom)
     distm(jatom,iatom)=distm(iatom,jatom)
   end do
 end do

 if (prtvoll>=3) then
   call  wrtout(std_out,'ioniondist: ionic distances:','COLL')
   call prmat(distm,natom,natom,natom,std_out)
 end if

 distm=anint(distm*10000_dp)/10000_dp           ! rounding needed else distm(iatom,jatom)/= distm(1,kdum) sometimes fails

 do iatom=1,natom
   if (option==1) then
     inm(iatom,:)=distm(iatom,:)
   else
     do jatom=iatom,natom
       kdum=1
       do while ( (kdum <= natom) .and. (distm(iatom,jatom)/= distm(1,kdum)) )
         kdum=kdum+1
       end do
       if (option==2) then
         inm(iatom,jatom)=interq(kdum)
       else if (option==3) then
         inm(iatom,jatom)=varlist(kdum)
       end if
       inm(jatom,iatom)=inm(iatom,jatom)
     end do
   end if
 end do

 if (prtvoll==2) then
   call  wrtout(std_out,'ioniondist: symmetrized matrix:','COLL')
   call prmat(distm,1,natom,natom,std_out)
 else if (prtvoll>=3) then
   call  wrtout(std_out,'ioniondist: symmetrized matrix:','COLL')
   call prmat(distm,natom,natom,natom,std_out)
 end if

end subroutine ioniondist
!!***
