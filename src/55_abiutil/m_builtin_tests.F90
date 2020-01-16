!!****m* ABINIT/m_builtin_tests
!! NAME
!!  m_builtin_tests
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA,XG,GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_builtin_tests

 use defs_basis
 use m_errors
 use m_abicore

 use m_io_tools,     only : open_file

 implicit none

 private
!!***

 public :: testfi
!!***

contains
!!***

!!****f* ABINIT/testfi
!!
!! NAME
!! testfi
!!
!! FUNCTION
!! Routine "Final test" for generation of the test report in the status file:
!! if it appears that the run was a "Built-in Test", then
!! compare the final values of different quantities to the reference
!! values, here hard-coded.
!!
!! INPUTS
!!  builtintest=number of the builtintest, from the input file.
!!  etotal=total energy (sum of 7 contributions) (hartree)
!!  filstat=name of the status file
!!  fred(3,natom)=symmetrized gradient of etotal with respect to tn
!!  natom=number of atoms in cell.
!!  strten(6)=components of the stress tensor (hartree/bohr^3)
!!  xred(3,natom)=reduced atomic coordinates
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine testfi(builtintest,etotal,filstat,fred,natom,strten,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: builtintest,natom
 real(dp),intent(in) :: etotal
 character(len=fnlen),intent(in) :: filstat
!arrays
 real(dp),intent(in) :: fred(3,natom),strten(6),xred(3,natom)

!Local variables-------------------------------
 character(len=fnlen) :: testname(7)='         '
 character(len=*), parameter :: format01000="(a,d22.14,a,d12.4)"
 character(len=500) :: msg
!scalars
 integer,parameter :: mtest=7
 integer :: iatom,ii,problem,tok,temp_unit
 real(dp) :: etot_mxdev,etot_ref,fred_mxdev,strten_mxdev,xred_mxdev
!arrays
 integer,parameter :: natom_test(mtest)=(/2,1,2,1,1,1,2/)
 real(dp) :: fred_ref(3,2),strten_ref(6),xred_ref(3,2)

! ***********************************************************************

 testname(1)='fast'
 testname(2)='v1'
 testname(3)='v5'
 testname(4)='bigdft'
 testname(5)='etsf_io'
 testname(6)='libxc'
 testname(7)='wannier90'

!---------------------------------------------------------

!Now, open the status file, and either delete it, or produce a report
 if (open_file(filstat,msg,newunit=temp_unit,form='formatted',status='unknown') /= 0) then
   MSG_ERROR(msg)
 end if

 if(builtintest==0)then

   close (temp_unit,status='delete')

 else

!  Note: all processors have their own file, so no special attention must be paid to the parallel case.
   write(temp_unit,*)
   write(temp_unit,*)'Status file, reporting on built-in test ',trim(testname(builtintest))
   write(temp_unit,*)

!  Define reference values, as well as maximum tolerable deviation
   if(builtintest==1)then

     etot_ref=-1.05814441948188d+00
     etot_mxdev=1.0d-9
     xred_ref(1:3,1)=(/ -0.65048430042634D-01 , 0.0_dp , 0.0_dp /)
     xred_ref(1:3,2)=(/  0.65048430042634D-01 , 0.0_dp , 0.0_dp /)
!    xred(*,*) are reduced coordinates
     xred_mxdev=1.0d-6
     fred_ref(1:3,1:2)= 0.0_dp
!    fred(*,*) are gradients with respect to reduced coordinates
     fred_mxdev=5.0d-4
     strten_ref(1:6)=(/ 0.149D-04  , 0.560D-04 , 0.560D-04 ,&
&     0.0_dp , 0.0_dp , 0.0_dp /)
     strten_mxdev=1.0d-5

   else if(builtintest==2)then

!    This value of etot is accurate to about 1.0d-12
     etot_ref=-.70811958266295D+02
!    Initialisation conditions might give fluctuations
     etot_mxdev=3.0d-7
     xred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     xred_mxdev=1.0d-12
     fred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     fred_mxdev=1.0d-12
!    This value of strten is accurate to at least 1.0d-8
     strten_ref(1:3)= 5.09324870E-03
     strten_ref(4:6)= 0.0_dp
     strten_mxdev=1.0d-6

   else if(builtintest==3)then

     etot_ref=-.892746696311772D+01
     etot_mxdev=1.0d-8
     xred_ref(1:3,1)=(/ -0.125_dp , 0.0_dp , 0.0_dp /)
     xred_ref(1:3,2)=(/  0.125_dp , 0.0_dp , 0.0_dp /)
     xred_mxdev=1.0d-12
     fred_ref(1:3,1)=(/ -.140263620278D+00 , 0.0_dp , 0.0_dp /)
     fred_ref(1:3,2)=(/  .140013483725D+00 , 0.0_dp , 0.0_dp /)
     fred_mxdev=1.0d-3
     strten_ref(1:6)=(/  1.3949d-3 ,  1.3643d-3 , 1.3643d-3 ,.0_dp ,  .0_dp  ,  .0_dp    /)
     strten_mxdev=1.0d-5

!    Bigdft
   else if(builtintest==4)then

!    This value of etot is accurate to about 1.0d-12
     etot_ref=-0.56231990141D+00
!    Initialisation conditions might give fluctuations
     etot_mxdev=3.0d-7
     xred_ref(1:3,1)=(/ 0.5_dp , 0.5_dp , 0.5_dp /)
     xred_mxdev=1.0d-12
     fred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     fred_mxdev=1.0d-12
     strten_ref(1)=-0.22537238382D-01
     strten_ref(2)=-0.22536232141D-01
     strten_ref(3)=-0.22529038043D-01
     strten_ref(4)= 0.39104928264D-06
     strten_ref(5)= 0.62823846408D-06
     strten_ref(6)= 0.38414125548E-09
     strten_mxdev=1.0d-7

!    ETSF_IO
   else if(builtintest==5)then

!    This value of etot is accurate to about 1.0d-12
     etot_ref=-0.33307825915795D+02
!    Initialisation conditions might give fluctuations
     etot_mxdev=3.0d-7
     xred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     xred_mxdev=1.0d-12
     fred_ref(1:3,1)=(/ 0.0_dp , -0.000000120877_dp , -0.000000164287_dp /)
     fred_mxdev=1.0d-7
!    This value of strten is accurate to at least 1.0d-8
     strten_ref(1)= 0.13569940015175D-01
     strten_ref(2)= 0.33822108610352D-01
     strten_ref(3)= 0.37991262607028D-01
     strten_ref(4:6)= 0.0_dp
     strten_mxdev=1.0d-6

!    LibXC
   else if(builtintest==6)then

     etot_ref=-0.55196688523897D+01
     etot_mxdev=5.0d-7
     xred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     xred_mxdev=1.0d-12
     fred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     fred_mxdev=1.0d-12
     strten_ref(1:3)= 0.13246699375127D-04
     strten_ref(4:6)= 0.0_dp
     strten_mxdev=1.0d-8

!    Wannier90
   else if(builtintest==7)then

!    This value of etot is accurate to about 1.0d-12
     etot_ref=-0.10620085133544D+02
!    Initialisation conditions might give fluctuations
     etot_mxdev=3.0d-7
     xred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     xred_ref(1:3,2)=(/ 0.25_dp , 0.25_dp , 0.25_dp /)
     xred_mxdev=1.0d-12
     fred_ref(1:3,1)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     fred_ref(1:3,2)=(/ 0.0_dp , 0.0_dp , 0.0_dp /)
     fred_mxdev=1.0d-12
!    This value of strten is accurate to at least 1.0d-8
     strten_ref(1:3)=0.31413922197317D-03
     strten_ref(4:6)= 0.0_dp
     strten_mxdev=1.0d-6

!    End of the reference value set up, for different tests
   end if

!  Compare reference and actual values

   problem=0

!  Take care of total energy
   if(abs(etot_ref-etotal)>etot_mxdev)then
     problem=1
     write(temp_unit,'(a)')' Error for total energy : '
     write(temp_unit,format01000)'        expected ',etot_ref,'  with maximum   deviation',etot_mxdev
     write(temp_unit,format01000)'        computed ',etotal,'  with effective deviation',abs(etotal-etot_ref)
   end if

!  Take care of nuclei positions
   tok=1
   do iatom=1,natom
     do ii=1,3
       if(abs(xred_ref(ii,iatom)-xred(ii,iatom))>xred_mxdev)then
         tok=0
         write(temp_unit, '(a,i1,a,i1,a)' )' Error for nuclei position xred(',ii,',',iatom,')'
         write(temp_unit,format01000)'        expected ',xred_ref(ii,iatom),'  with maximum   deviation',xred_mxdev
         write(temp_unit,format01000)'        computed ',xred(ii,iatom),'  with effective deviation',&
&         abs( xred(ii,iatom)-xred_ref(ii,iatom) )
       end if
     end do
   end do
   if(tok==0)problem=1

!  Take care of forces
   tok=1
   do iatom=1,natom
     do ii=1,3
       if(abs(fred_ref(ii,iatom)-fred(ii,iatom))&
&       >fred_mxdev)then
         tok=0
         write(temp_unit, '(a,i1,a,i1,a)' )' Error for force fred(',ii,',',iatom,')'
         write(temp_unit,format01000)'        expected ',fred_ref(ii,iatom),'  with maximum   deviation',fred_mxdev
         write(temp_unit,format01000)'        computed ',fred(ii,iatom),'  with effective deviation',&
&         abs( fred(ii,iatom)-fred_ref(ii,iatom) )
       end if
     end do
   end do
   if(tok==0)problem=1

!  Take care of stress
   tok=1
   do ii=1,6
     if(abs(strten_ref(ii)-strten(ii))>strten_mxdev)then
       tok=0
       write(temp_unit,'(a,i1,a)')' Error for stress strten(',ii,')'
       write(temp_unit,format01000)'        expected ',strten_ref(ii),'  with maximum   deviation',strten_mxdev
       write(temp_unit,format01000)'        computed ',strten(ii),'  with effective deviation',&
&       abs( strten(ii)-strten_ref(ii) )
     end if
   end do
   if(tok==0)problem=1

   if(problem==0)then
     write(temp_unit,'(a)')' ==> The run finished cleanly.'
     write(temp_unit,'(a,a)')'     Moreover, comparison of the total energy, and other (few) ',&
&     'relevant quantities with reference values has been successful.'
     write(temp_unit,'(a)')'     This does not mean that no problem is present, however. '
     write(temp_unit,'(a)')'     Please run the complete set of ABINIT tests to gain a better confidence in your installation.'
   end if

   write(temp_unit,*)

   close(temp_unit)

 end if !  End of the choice between produce a report, and produce no report

end subroutine testfi
!!***

end module m_builtin_tests
!!***
