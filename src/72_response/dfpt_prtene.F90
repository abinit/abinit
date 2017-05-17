!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_prtene
!!
!! NAME
!! dfpt_prtene
!!
!! FUNCTION
!! Print components of second derivative of total energy in nice format
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (XG, DRH, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! eberry=energy associated with Berry phase
!! edocc=correction to 2nd-order total energy coming from changes of occupation
!! eeig0=0th-order eigenenergies part of 2nd-order total energy
!! eew=Ewald part of 2nd-order total energy
!! efrhar=hartree frozen-wavefunction part of 2nd-order tot. en.
!! efrkin=kinetic frozen-wavefunction part of 2nd-order tot. en.
!! efrloc=local psp. frozen-wavefunction part of 2nd-order tot. en.
!! efrnl=nonlocal psp. frozen-wavefunction part of 2nd-order tot. en
!! efrx1=xc core corr.(1) frozen-wavefunction part of 2nd-order tot. en
!! efrx2=xc core corr.(2) frozen-wavefunction part of 2nd-order tot. en
!! ehart01=inhomogeneous 1st-order Hartree part of 2nd-order total energy
!!   for strain perturbation only (zero otherwise, and not used)
!! ehart1=1st-order Hartree part of 2nd-order total energy
!! eii=pseudopotential core part of 2nd-order total energy
!! ek0=0th-order kinetic energy part of 2nd-order total energy.
!! ek1=1st-order kinetic energy part of 2nd-order total energy.
!! eloc0=0th-order local (psp+vxc+Hart) part of 2nd-order total energy
!! elpsp1=1st-order local pseudopot. part of 2nd-order total energy.
!! enl0=0th-order nonlocal pseudopot. part of 2nd-order total energy.
!! enl1=1st-order nonlocal pseudopot. part of 2nd-order total energy.
!! eovl1=1st-order change of wave-functions overlap, part of 2nd-order energy
!!       PAW only - Eq(79) and Eq(80) of PRB 78, 035105 (2008)
!! epaw1=1st-order PAW on-site part of 2nd-order total energy.
!! evdw=DFT-D semi-empirical part of 2nd-order total energy
!! exc1=1st-order exchange-correlation part of 2nd-order total energy
!! iout=unit number to which output is written
!! ipert=type of the perturbation
!! natom=number of atoms in unit cell
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! usevdw= flag set to 1 if vdw DFT-D semi-empirical potential is in use
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! all energies in Hartree
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_prtene(berryopt,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&
&  ehart01,ehart1,eii,ek0,ek1,eloc0,elpsp1,enl0,enl1,eovl1,epaw1,evdw,exc1,iout,ipert,natom,&
&  usepaw,usevdw)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_prtene'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: berryopt,iout,ipert,natom,usepaw,usevdw
 real(dp),intent(in) :: eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1
 real(dp),intent(in) :: efrx2,ehart01,ehart1,eii,ek0,ek1,eloc0,elpsp1,enl0,enl1
 real(dp),intent(in) :: eovl1,epaw1,evdw,exc1

!Local variables -------------------------
!scalars
 integer :: nn
 logical :: berry_activated
 real(dp) :: enl1_effective,erelax,etotal
 character(len=10) :: numb
 character(len=10),parameter :: numbstr(20) = &
&  (/'One       ','Two       ','Three     ','Four      ','Five      ', &
&    'Six       ','Seven     ','Eight     ','Nine      ','Ten       ', &
&    'Eleven    ','Twelve    ','Thirteen  ','Fourteen  ','Fifteen   ', &
&    'Sixteen   ','Seventeen ','Eighteen  ','Nineteen  ','Twenty    '/)
 character(len=500) :: message

! *********************************************************************

!Count and print the number of components of 2nd-order energy
!MT feb 2015: this number is wrong! Should change it but
!             need to change a lot of ref. files
 berry_activated=(berryopt== 4.or.berryopt== 6.or.berryopt== 7.or. &
& berryopt==14.or.berryopt==16.or.berryopt==17)
 if (ipert==natom+1) nn=8
 if (ipert==natom+2) nn=7
 if (ipert>=1.and.ipert<=natom) nn=13
 if (ipert==natom+3.or.ipert==natom+4) nn=17
 if (ipert==natom+2.and.berry_activated) nn=nn+1
 if (usepaw==1) nn=nn+1
 if (usevdw==1) nn=nn+1
 if (ipert==natom+10.or.ipert==natom+11) nn=1 ! means nothing,
! because we do not compute derivatives of the energy in this case
 write(message, '(4a)' ) ch10,&
& ' ',trim(numbstr(nn)),' components of 2nd-order total energy (hartree) are '
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 numb='1,2,3'
 write(message, '(3a)' )&
& ' ',trim(numb),': 0th-order hamiltonian combined with 1st-order wavefunctions'
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')
 write(message, '(a,es17.8,a,es17.8,a,es17.8)' )&
& '     kin0=',ek0,   ' eigvalue=',eeig0,'  local=',eloc0
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 numb='4,5,6';if( ipert==natom+3.or.ipert==natom+4) numb='4,5,6,7'
 write(message, '(3a)' )&
& ' ',trim(numb),': 1st-order hamiltonian combined with 1st and 0th-order wfs'
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')
 if(ipert/=natom+1.and.ipert/=natom+2)then
   write(message, '(a,es17.8,a,es17.8,a,es17.8,a,a)' ) &
&   ' loc psp =',elpsp1,'  Hartree=',ehart1,'     xc=',exc1,ch10,&
&   ' note that "loc psp" includes a xc core correction that could be resolved'
 else if(ipert==natom+1) then
   write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&   '     kin1=',ek1,   '  Hartree=',ehart1,'     xc=',exc1
 else if(ipert==natom+2) then
   write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&   '    dotwf=',enl1,  '  Hartree=',ehart1,'     xc=',exc1
 end if
 if(ipert==natom+3 .or. ipert==natom+4) then
   write(message, '(a,es17.8,a,es17.8,a,es17.8,a,a,es17.8)' ) &
&   ' loc psp =',elpsp1,'  Hartree=',ehart1,'     xc=',exc1,ch10,&
&   '     kin1=',ek1
 end if
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 enl1_effective=enl1;if (ipert==natom+2) enl1_effective=zero
 numb='7,8,9';if( ipert==natom+3.or.ipert==natom+4) numb='8,9,10'
 write(message, '(5a,es17.8,a,es17.8,a,es17.8)' )&
& ' ',trim(numb),': eventually, occupation + non-local contributions',ch10,&
& '    edocc=',edocc,'     enl0=',enl0,'   enl1=',enl1_effective
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 if (usepaw==1) then
   numb='10';if( ipert==natom+3.or.ipert==natom+4) numb='11'
   write(message, '(3a,es17.8)' )&
&   ' ',trim(numb),': eventually, PAW "on-site" Hxc contribution: epaw1=',epaw1
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 if(ipert/=natom+10 .and.ipert/=natom+11) then
   erelax=0.0_dp
   if(ipert>=1.and.ipert<=natom)then
     erelax=ek0+edocc+eeig0+eloc0+elpsp1+ehart1+exc1+enl0+enl1+epaw1
   else if(ipert==natom+1.or.ipert==natom+2)then
     erelax=ek0+edocc+eeig0+eloc0+ek1+ehart1+exc1+enl0+enl1+epaw1
   else if(ipert==natom+3.or.ipert==natom+4)then
     erelax=ek0+edocc+eeig0+eloc0+ek1+elpsp1+ehart1+exc1+enl0+enl1+epaw1
   end if
   enl1_effective=enl1
   if (ipert==natom+1.or.ipert==natom+2) then
     if (1.0_dp+enl1/10.0_dp==1.0_dp) enl1_effective=zero
   end if

   numb='1-9';if (usepaw==1) numb='1-10'
   if( ipert==natom+3.or.ipert==natom+4) then
     numb='1-10';if (usepaw==1) numb='1-11'
   end if
   write(message, '(5a,es17.8)' )&
&   ' ',trim(numb),' gives the relaxation energy (to be shifted if some occ is /=2.0)',&
&   ch10,'   erelax=',erelax
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 if(ipert>=1.and.ipert<=natom)then

   numb='10,11,12';if (usepaw==1) numb='11,12,13'
   write(message, '(4a)' )&
&   ' ',trim(numb),' Non-relaxation  contributions : ',&
&   'frozen-wavefunctions and Ewald'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&   ' fr.local=',efrloc,' fr.nonlo=',efrnl,'  Ewald=',eew
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')

   write(message, '(a,es16.6)' )' dfpt_prtene : non-relax=',efrloc+efrnl+eew
   call wrtout(std_out,message,'COLL')

   numb='13,14';if (usepaw==1) numb='14,15'
   write(message, '(3a)' )&
&   ' ',trim(numb),' Frozen wf xc core corrections (1) and (2)'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es17.8,a,es17.8)' ) &
&   ' frxc 1  =',efrx1,'  frxc 2 =',efrx2
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   if (usepaw==1) then
     numb='16'
     write(message, '(5a,es17.8)' )&
&     ' ',trim(numb),' Contribution from 1st-order change of wavefunctions overlap',&
&     ch10,' eovl1 =',eovl1
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if
   if (usevdw==1) then
     numb='15';if (usepaw==1) numb='17'
     write(message, '(3a,es17.8)' )&
&     ' ',trim(numb),' Contribution from van der Waals DFT-D: evdw =',evdw
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if

   write(message, '(a)' )' Resulting in : '
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   etotal=erelax+eew+efrloc+efrnl+efrx1+efrx2+evdw
   write(message, '(a,e20.10,a,e22.12,a)' ) &
&   ' 2DEtotal=',etotal,' Ha. Also 2DEtotal=',etotal*Ha_eV,' eV'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es20.10,a,es20.10,a)' ) &
&   '    (2DErelax=',erelax,' Ha. 2DEnonrelax=',etotal-erelax,' Ha)'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es20.10,a,a)' ) &
&   '    (  non-var. 2DEtotal :',&
&   0.5_dp*(elpsp1+enl1)+eovl1+eew+efrloc+efrnl+efrx1+efrx2+evdw,' Ha)',ch10
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')

 else if(ipert==natom+1.or.ipert==natom+2)then
   if (ipert==natom+1.and.usepaw==1) then
     numb='11'
     write(message, '(5a,es17.8)' )&
&     ' ',trim(numb),' Contribution from 1st-order change of wavefunctions overlap',&
&     ch10,' eovl1 =',eovl1
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if
   write(message,*)' No Ewald or frozen-wf contrib.:',&
&   ' the relaxation energy is the total one'
   if(berry_activated) then
     numb='10';
     write(message,'(3a,es20.10)')' ',trim(numb),' Berry phase energy :',eberry
   end if
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   etotal=erelax
   write(message, '(a,e20.10,a,e22.12,a)' ) &
&   ' 2DEtotal=',etotal,' Ha. Also 2DEtotal=',etotal*Ha_eV,' eV'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es20.10,a)' ) &
&   '    (  non-var. 2DEtotal :',0.5_dp*(ek1+enl1_effective)+eovl1,' Ha)'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')

 else if(ipert==natom+3 .or. ipert==natom+4) then
   numb='11,12,13';if (usepaw==1) numb='12,13,14'
   write(message, '(4a)' )&
&   ' ',trim(numb),' Non-relaxation  contributions : ',&
&   'frozen-wavefunctions and Ewald'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&   '  fr.hart=',efrhar,'   fr.kin=',efrkin,' fr.loc=',efrloc
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')

   numb='14,15,16';if (usepaw==1) numb='15,16,17'
   write(message, '(4a)' )&
&   ' ',trim(numb),' Non-relaxation  contributions : ',&
&   'frozen-wavefunctions and Ewald'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&   '  fr.nonl=',efrnl,'    fr.xc=',efrx1,'  Ewald=',eew
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')

   numb='17';if (usepaw==1) numb='18'
   write(message, '(4a)' )&
&   ' ',trim(numb),' Non-relaxation  contributions : ',&
&   'pseudopotential core energy'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es17.8)' ) &
&   '  pspcore=',eii
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   if (usepaw==1) then
     numb='19'
     write(message, '(5a,es17.8)' )&
&     ' ',trim(numb),' Contribution from 1st-order change of wavefunctions overlap',&
&     ch10,' eovl1 =',eovl1
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if
   if (usevdw==1) then
     numb='18';if (usepaw==1) numb='20'
     write(message, '(3a,es17.8)' )&
&     ' ',trim(numb),' Contribution from van der Waals DFT-D: evdw =',evdw
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if

   write(message, '(a,es16.6)' )' dfpt_prtene : non-relax=',&
&   efrhar+efrkin+efrloc+efrnl+efrx1+eew+evdw
   call wrtout(std_out,message,'COLL')
   write(message, '(a)' )' Resulting in : '
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   etotal=erelax+efrhar+efrkin+efrloc+efrnl+efrx1+eew+eii+evdw
   write(message, '(a,e20.10,a,e22.12,a)' ) &
&   ' 2DEtotal=',etotal,' Ha. Also 2DEtotal=',etotal*Ha_eV,' eV'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es20.10,a,es20.10,a)' ) &
&   '    (2DErelax=',erelax,' Ha. 2DEnonrelax=',etotal-erelax,' Ha)'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es20.10,a,a)' ) &
&   '    (  non-var. 2DEtotal :',&
&   0.5_dp*(elpsp1+enl1+ek1+ehart01)+eovl1+&
&   efrhar+efrkin+efrloc+efrnl+efrx1+eew+eii+evdw,' Ha)',ch10
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

end subroutine dfpt_prtene
!!***
