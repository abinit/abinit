!{\src2tex{textfont=tt}}
!!****f* ABINIT/dos_hdr_write
!!
!! NAME
!! dos_hdr_write
!!
!! FUNCTION
!! Write the header of the DOS files, for both
!! smearing and tetrahedron methods. Also compute the
!! minimum, maximum energies, the energy increment
!! and the number of points for the DOS.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG, AF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! buffer=approximative buffer energy for the output of the DOS
!!  (beyond the max and min energy values).
!! dosdeltae=DOS delta of Energy (if zero, take default values)
!! eigen(mband*nkpt*nsppol)=eigenvalues (input or init to large number), hartree
!! fermie=fermi energy useful for band alignment...
!! mband=maximum number of bands
!! nband(nkpt*nsppol)=number of bands at each k point
!! nkpt=number of k points
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! occopt=option for occupancies, or re-smearing scheme if dblsmr /= 0
!! prtdos=1 for smearing technique, 2 or 3 for tetrahedron technique
!! tphysel="physical" electronic temperature with FD occupations
!! tsmear=smearing width (or temperature)
!! unitdos=unit number of output of the DOS.
!!
!! OUTPUT
!! deltaene=increment of DOS energy arguments
!! enemax=maximal value of the DOS energy argument
!! enemin=minimal value of the DOS energy argument
!! nene=number of DOS energy argument
!!
!! NOTES
!!
!! PARENTS
!!      gaus_dos,getnel,tetrahedron
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dos_hdr_write(buffer,deltaene,dosdeltae,&
&  eigen,enemax,enemin,fermie,mband,nband,nene,&
&  nkpt,nsppol,occopt,prtdos,tphysel,tsmear,unitdos)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dos_hdr_write'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkpt,nsppol,occopt,prtdos,unitdos
 integer,intent(out) :: nene
 real(dp),intent(in) :: buffer,dosdeltae,fermie,tphysel,tsmear
 real(dp),intent(out) :: deltaene,enemax,enemin
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer :: bantot
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,*) ' dos_hdr_write : enter '
!stop
!ENDDEBUG

 bantot=sum(nband(:))

!Choose the lower and upper energies
 enemax=maxval(eigen(1:bantot))+buffer
 enemin=minval(eigen(1:bantot))-buffer

!Extend the range to a nicer value
 enemax=0.1_dp*ceiling(enemax*10._dp)
 enemin=0.1_dp*floor(enemin*10._dp)

!Choose the energy increment
 if(abs(dosdeltae)<tol10)then
   deltaene=0.001_dp
   if(prtdos>=2)deltaene=0.0005_dp ! Higher resolution possible (and wanted) for tetrahedron
 else
   deltaene=dosdeltae
 end if
 nene=nint((enemax-enemin)/deltaene)+1

!Write the DOS file
 write(message, '(7a,i2,a,i5,a,i4)' ) "#",ch10, &
& '# ABINIT package : DOS file  ',ch10,"#",ch10,&
& '# nsppol =',nsppol,', nkpt =',nkpt,', nband(1)=',nband(1)
 call wrtout(unitdos,message,'COLL')

 if(prtdos==1.or.prtdos==4)then
   write(message, '(a,i2,a,f6.3,a,f6.3,a)' )  &
&   '# Smearing technique, occopt =',occopt,', tsmear=',tsmear,' Hartree, tphysel=',tphysel,' Hartree'
 else
   write(message, '(a)' ) &
&   '# Tetrahedron method '
 end if
 call wrtout(unitdos,message,'COLL')

 if(mband*nkpt*nsppol>=3) then
   write(message, '(a,3f8.3,2a)' ) &
&   '# For identification : eigen(1:3)=',eigen(1:3),ch10,"#"
 else
   write(message, '(a,3f8.3)' ) &
&   '# For identification : eigen=',eigen
   write(message, '(3a)')trim(message),ch10,"#"
 end if
 call wrtout(unitdos,message,'COLL')

 write(message, '(a,f16.8)' ) &
& '# Fermi energy : ', fermie
 call wrtout(unitdos,message,'COLL')

 if(prtdos==1)then
   write(message, '(5a)' ) "#",ch10,&
&   '# The DOS (in electrons/Hartree/cell) and integrated DOS (in electrons/cell),',&
&   ch10,'# as well as the DOS with tsmear halved and doubled, are computed,'
 else if(prtdos==2)then
   write(message, '(3a)' ) "#",ch10,&
&   '# The DOS (in electrons/Hartree/cell) and integrated DOS (in electrons/cell) are computed,'
 else if(prtdos==3.or.prtdos==4)then
   write(message, '(5a)' ) "#",ch10,&
&   '# The local DOS (in electrons/Hartree for one atomic sphere)',ch10,&
&   '# and integrated local DOS (in electrons for one atomic sphere) are computed.'
 else if(prtdos==5)then
   write(message, '(9a)' ) "#",ch10,&
&   '# The spin component DOS (in electrons/Hartree/cell)',ch10,&
&   '# and integrated spin component DOS (in electrons/cell) are computed.',ch10,&
&   '# Remember that the wf are eigenstates of S_z and S^2, not S_x and S_y',ch10,&
&   '#   so the latter will not always sum to 0 for paired electronic states.'
 end if
 call wrtout(unitdos,message,'COLL')

 write(message, '(a,i5,a,a,a,f9.4,a,f9.4,a,f8.5,a,a,a)' )&
& '# at ',nene,' energies (in Hartree) covering the interval ',ch10,&
& '# between ',enemin,' and ',enemax,' Hartree by steps of ',&
& deltaene,' Hartree.',ch10,"#"
 call wrtout(unitdos,message,'COLL')

 if(prtdos==1)then
   write(message, '(a,a)' )&
&   '#       energy        DOS       Integr. DOS   ',&
&   '     DOS           DOS    '
   call wrtout(unitdos,message,'COLL')

   write(message, '(a)' )&
&   '#                                              (tsmear/2)    (tsmear*2) '
   call wrtout(unitdos,message,'COLL')
 else
   write(message, '(a)' ) '#       energy        DOS '
 end if

!DEBUG
!write(std_out,*) ' dos_hdr_write : exit '
!ENDDEBUG

end subroutine dos_hdr_write
!!***
