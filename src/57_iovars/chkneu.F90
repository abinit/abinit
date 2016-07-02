!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkneu
!! NAME
!! chkneu
!!
!! FUNCTION
!! Check neutrality of system based on band occupancies and
!! valence charges of pseudo-atoms.
!! Eventually initialize occ if occopt==1 or 3...8
!! Also return nelect, the number of valence electron per unit cell
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  charge=number of electrons missing (+) or added (-) to system (usually 0)
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | iscf= if>0, SCF calculation ; if<=0, non SCF calculation (wtk might
!!   |  not be defined)
!!   | natom=number of atoms in unit cell
!!   | nband(nkpt*nsppol)=number of bands at each k point
!!   | nkpt=number of k points
!!   | nspinor=number of spinorial components of the wavefunctions
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | ntypat=number of pseudopotentials
!!   | positron=0 if electron GS calculation
!!   |          1 if positron GS calculation
!!   |          2 if electron GS calcultaion in presence of the positron
!!   | typat(natom)=atom type (integer) for each atom
!!   | wtk(nkpt)=k point weights (defined if iscf>0 or iscf==-3)
!!   | ziontypat(ntypat)=ionic charge of each pseudoatom
!!  occopt=option for occupancies
!!
!! OUTPUT
!!  Writes warning and/or aborts if error condition exists
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | nelect=number of valence electrons per unit cell
!!   |  (from counting valence electrons in psps, and taking into
!!   |   account the input variable "charge")
!!
!! SIDE EFFECTS
!! Input/Output :
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | occ_orig(dtset%nband(1)*nkpt*nsppol)=occupation numbers for each band and k point
!!   |   must be input for occopt==0 or 2,
!!   |   will be an output for occopt==1 or 3 ... 8
!!
!! NOTES
!!
!! PARENTS
!!      invars2
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine chkneu(charge,dtset,occopt)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkneu'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: occopt
 real(dp),intent(in) :: charge
 type(dataset_type),intent(inout) :: dtset

!Local variables-------------------------------
!scalars
 integer :: bantot,iatom,iband,ii,ikpt,isppol,nocc
 real(dp) :: maxocc,nelect_occ,nelect_spin,occlast,sign_spin,zval
 character(len=500) :: message
!arrays
 real(dp),allocatable :: tmpocc(:)

! *************************************************************************

!(1) count nominal valence electrons according to ziontypat
 zval=zero
 do iatom=1,dtset%natom
   zval=zval+dtset%ziontypat(dtset%typat(iatom))
 end do
 if (dtset%positron/=1) then
   dtset%nelect=zval-charge
 else
   dtset%nelect=one
 end if

!DEBUG
 write(std_out,*)ch10
 write(std_out,*)' chkneu : enter, dtset%nelect=',dtset%nelect
 write(std_out,*)' occopt,dtset%nsppol,dtset%nspden=',occopt,dtset%nsppol,dtset%nspden
!ENDDEBUG

!(2) Optionally initialize occ with semiconductor occupancies
!(even for a metal : at this stage, the eigenenergies are unknown)
!Note that nband(1)=nband(2) in this section, as occopt=2 is avoided.
 if(occopt==1 .or. (occopt>=3 .and. occopt<=8) )then
!  Here, initialize a real(dp) variable giving the
!  maximum occupation number per band
   maxocc=2.0_dp/real(dtset%nsppol*dtset%nspinor,dp)

!  Determine the number of bands fully or partially occupied
   nocc=(dtset%nelect-1.0d-8)/maxocc + 1
!  Occupation number of the highest level
   occlast=dtset%nelect-maxocc*(nocc-1)

!DEBUG
   write(std_out,*)' maxocc,nocc,occlast=',maxocc,nocc,occlast
!ENDDEBUG


!  The number of allowed bands must be sufficiently large
   if( nocc<=dtset%nband(1)*dtset%nsppol .or. dtset%iscf==-2) then

     if(dtset%iscf==-2 .and. nocc>dtset%nband(1)*dtset%nsppol)nocc=dtset%nband(1)*dtset%nsppol

!    First treat the case where the spin magnetization is not imposed, is zero with nspden==1, or has sufficient flexibility
!    for a target not to be matched at the initialisation, but later
     if(abs(dtset%spinmagntarget+99.99_dp)<tol8 .or. (dtset%nspden==4) .or. &
&     (abs(dtset%spinmagntarget)<tol8.and.dtset%nspden==1))then

!      Use a temporary array for defining occupation numbers
       ABI_ALLOCATE(tmpocc,(dtset%nband(1)*dtset%nsppol))
!      First do it for fully occupied bands
       if (1<nocc) tmpocc(1:nocc-1)=maxocc
!      Then, do it for highest occupied band
       if (1<=nocc) tmpocc(nocc)=occlast
!      Finally do it for eventual unoccupied bands
       if ( nocc<dtset%nband(1)*dtset%nsppol ) tmpocc(nocc+1:dtset%nband(1)*dtset%nsppol)=0.0_dp

!      Now copy the tmpocc array in the occ array, taking into account the spin
       if(dtset%nsppol==1)then
         do ikpt=1,dtset%nkpt
           dtset%occ_orig(1+(ikpt-1)*dtset%nband(1):ikpt*dtset%nband(1))=tmpocc(:)
         end do
       else
         do ikpt=1,dtset%nkpt
           do iband=1,dtset%nband(1)
             do isppol=1,dtset%nsppol
               dtset%occ_orig(iband+dtset%nband(1)*(ikpt-1+dtset%nkpt*(isppol-1))) =  &
&               tmpocc(isppol+dtset%nsppol*(iband-1))
             end do
           end do
         end do
       end if
       ABI_DEALLOCATE(tmpocc)

!      Second, treat the case in which one imposes the spin magnetization (only possible for nspden==2)
!      Also treat antiferromagnetic case (nsppol==1, nspden==2), although spinmagntarget must be zero
     else if(abs(dtset%spinmagntarget+99.99_dp)>1.0d-10)then
       do isppol=1,dtset%nsppol
         sign_spin=real(3-2*isppol,dp)
         nelect_spin=half*(dtset%nelect*maxocc+sign_spin*dtset%spinmagntarget)

!        DEBUG 
         write(std_out,*)' isppol,sign_spin,nelect_spin=',isppol,sign_spin,nelect_spin
!        ENDDEBUG
!        Determines the last state, and its occupation
         if(abs(nint(nelect_spin)-nelect_spin)<tol10)then
           nocc=nint(nelect_spin/maxocc)
           occlast=maxocc
         else
           nocc=ceiling(nelect_spin/maxocc)
           occlast=nelect_spin-(real(nocc,dp)-one)*maxocc
         end if
!        DEBUG 
         write(std_out,*)' dtset%nband(1),maxocc,occlast=',dtset%nband(1),maxocc,occlast
!        ENDDEBUG
         if(dtset%nband(1)*nint(maxocc)<nocc)then
           write(message, '(a,i4,a, a,2i6,a, a,es16.6,a, a,es16.6,6a)' )&
&           'Initialization of occ, with nspden=',dtset%nspden,ch10,&
&           'number of bands=',dtset%nband(1:2),ch10,&
&           'number of electrons=',dtset%nelect,ch10,&
&           'and spinmagntarget=',dtset%spinmagntarget,ch10,&
&           'This combination is not possible, because of a lack of bands.',ch10,&
&           'Action : modify input file ... ',ch10,&
&           '(you should likely increase nband, but also check nspden, nspinor, nsppol, and spinmagntarget)'
           MSG_ERROR(message)
         end if
         do ikpt=1,dtset%nkpt
!          Fill all bands, except the upper one
           if(dtset%nband(1)>1)then
             do iband=1,nocc-1
               dtset%occ_orig(iband+dtset%nband(1)*(ikpt-1+dtset%nkpt*(isppol-1)))=maxocc
             end do
           end if
!          Fill the upper occupied band
           dtset%occ_orig(nocc+dtset%nband(1)*(ikpt-1+dtset%nkpt*(isppol-1)))=occlast
         end do
!        DEBUG 
         write(std_out,*)' dtset%occ_orig(1)=',dtset%occ_orig(1)
!        ENDDEBUG
       end do

     else
       write(message, '(a,i4,a,a,es16.6,6a)' )&
&       'Initialization of occ, with nspden=',dtset%nspden,ch10,&
&       'and spinmagntarget=',dtset%spinmagntarget,ch10,&
&       'This combination is not possible.',ch10,&
&       'Action : modify input file ... ',ch10,&
&       '(check nspden, nspinor, nsppol and spinmagntarget)'
       MSG_ERROR(message)
     end if

!    Now print the values
     if(dtset%nsppol==1)then

       write(message, '(a,i4,a,a)' ) &
&       ' chkneu : initialized the occupation numbers for occopt= ',occopt,', spin-unpolarized or antiferromagnetic case : '
!      Add iout as argument to uncomment this line.
!      call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
       do ii=0,(dtset%nband(1)-1)/12
         write(message,'(12f6.2)') dtset%occ_orig( 1+ii*12 : min(12+ii*12,dtset%nband(1)) )
!        call wrtout(iout,message,'COLL')
         call wrtout(std_out,message,'COLL')
       end do

     else

       write(message, '(a,i4,a,a)' ) &
&       ' chkneu : initialized the occupation numbers for occopt= ',occopt,&
&       ch10,'    spin up   values : '
!      call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
       do ii=0,(dtset%nband(1)-1)/12
         write(message,'(12f6.2)') dtset%occ_orig( 1+ii*12 : min(12+ii*12,dtset%nband(1)) )
!        call wrtout(iout,message,'COLL')
         call wrtout(std_out,message,'COLL')
       end do
       write(message, '(a)' ) '    spin down values : '
!      call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
       do ii=0,(dtset%nband(1)-1)/12
         write(message,'(12f6.2)') &
&         dtset%occ_orig( 1+ii*12+dtset%nkpt*dtset%nband(1) : min(12+ii*12,dtset%nband(1))+dtset%nkpt*dtset%nband(1) )
!        call wrtout(iout,message,'COLL')
         call wrtout(std_out,message,'COLL')
       end do

     end if

!    Here, treat the case when the number of allowed bands is not large enough
   else
     write(message, '(a,i4,a,a,a,a,a,a,a,a)' )&
&     'Initialization of occ, with occopt=',occopt,ch10,&
&     'There are not enough bands to get charge balance right',ch10,&
&     'Action: modify input file ... ',ch10,&
&     '(check the pseudopotential charges, the variable charge,',ch10,&
&     'and the declared number of bands, nband)'
     MSG_ERROR(message)
   end if
 end if

!The remaining of the routine is for SCF runs and special options
 if(dtset%iscf>0 .or. dtset%iscf==-1 .or. dtset%iscf==-3)then

!  (3) count electrons in bands (note : in case occ has just been
!  initialized, point (3) and (4) is a trivial test
   nelect_occ=0.0_dp
   bantot=0
   do isppol=1,dtset%nsppol
     do ikpt=1,dtset%nkpt
       do iband=1,dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
         bantot=bantot+1
         nelect_occ=nelect_occ+dtset%wtk(ikpt)*dtset%occ_orig(bantot)
       end do
     end do
   end do

!  (4) if dtset%iscf/=-3, dtset%nelect must equal nelect_occ
!  if discrepancy exceeds tol11, give warning;  tol8, stop with error

   if (abs(nelect_occ-dtset%nelect)>tol11 .and. dtset%iscf/=-3) then

!    There is a discrepancy
     write(message, &
&     '(a,a,e16.8,a,e16.8,a,a,a,e22.14,a,a,a,a,a,a,a)' ) ch10,&
&     ' chkneu: nelect_occ=',nelect_occ,', zval=',zval,',',ch10,&
&     '         and input value of charge=',charge,',',ch10,&
&     '   nelec_occ is computed from occ and wtk',ch10,&
&     '   zval is nominal charge of all nuclei, computed from zion (read in psp),',ch10,&
&     '   charge is an input variable (usually 0).'
     call wrtout(std_out,message,'COLL')

     if (abs(nelect_occ-dtset%nelect)>tol8) then
!      The discrepancy is severe
       write(message,'(a,a,e9.2,a,a)')ch10,&
&       'These must obey zval-nelect_occ=charge to better than ',tol8,ch10,&
&       ' This is not the case. '
     else
!      The discrepancy is not so severe
       write(message, '(a,a,e9.2)' )ch10,'These should obey zval-nelect_occ=charge to better than ',tol11
     end if
     MSG_WARNING(message)

     write(message, '(a,a,a,a,a,a)' ) &
&     'Action : check input file for occ,wtk, and charge.',ch10,&
&     'Note that wtk is NOT automatically normalized when occopt=2,',ch10,&
&     'but IS automatically normalized otherwise.',ch10
     call wrtout(std_out,message,'COLL')

!    If the discrepancy is severe, stop
     if (abs(nelect_occ-dtset%nelect)>tol8)then
       MSG_ERROR(message)
     end if

   end if

!  End the condition dtset%iscf>0 or -1 or -3 .
 end if

end subroutine chkneu
!!***
