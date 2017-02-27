!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_etot
!! NAME
!! dfpt_etot
!!
!! FUNCTION
!! Assemble different contributions to the variational part of the
!! 2nd derivative of total energy
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (XG, DRH, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  berryopt= 4/14: electric field is on; berryopt = 6/7/16/17: electric displacement field is on;
!!  eberry=energy associated with Berry phase
!!  edocc=correction to 2nd-order total energy coming from changes of occupation
!!  ehart1=1st-order Hartree part of 2nd-order total energy
!!  eeig0=0th-order eigenenergies part of 2nd-order total energy
!!  eew=2nd derivative of Ewald energy (hartree)
!!  efrhar=contrib. from frozen-wavefunction, hartree energy, to the 2nd-derivative of total energy
!!  efrkin=contrib. from frozen-wavefunction, kinetic energy, to the 2nd-derivative of total energy
!!  efrloc=contrib. from frozen-wavefunction, local potential, to the 2nd-derivative of total energy
!!  efrnl=contribution from frozen-wavefunction, non-local potential, to the 2nd-derivative of total energy
!!  efrx1=contrib. from frozen-wavefunction, xc core correction(1), to the 2nd-derivative of total energy
!!  efrx2=contribution from frozen-wavefunction, xc core correction(2),
!!           to the second-derivative of total energy.
!!  ek0=0th-order kinetic energy part of 2nd-order total energy.
!!  ek1=1st-order kinetic energy part of 2nd-order total energy.
!!  eii=2nd derivative of pseudopotential core energy (hartree)
!!  eloc0=0th-order local (psp+vxc+Hart) part of 2nd-order total energy
!!  elpsp1=1st-order local pseudopot. part of 2nd-order total energy.
!!  enl0=0th-order nonlocal pseudopot. part of 2nd-order total energy.
!!  enl1=1st-order nonlocal pseudopot. part of 2nd-order total energy.
!!  epaw1=1st-order PAW on-sitew part of 2nd-order total energy.
!!  evdw=DFT-D semi-empirical part of 2nd-order total energy
!!  exc1=1st-order exchange-correlation part of 2nd-order total energy
!!  ipert=type of the perturbation
!!  natom=number of atoms
!!  optene=option for the computation of 2nd-order total energy
!!         (-1=no computation; 0=direct scheme; 1=double-counting scheme)
!!
!! OUTPUT
!!  deltae=change in energy between the previous and present SCF cycle
!!         and previous SCF cycle.
!!  etotal=2nd-order total energy
!!  evar=variational part of the 2nd-order total energy
!!
!! SIDE EFFECTS
!! input/output
!! elast=previous value of the 2nd-order total energy, needed to compute deltae,
!!      then updated (cannot simply be saved, because set to zero
!!      at each new call of dfpt_scfcv).
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_etot(berryopt,deltae,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,&
&                efrnl,efrx1,efrx2,ehart1,ek0,ek1,eii,elast,eloc0,elpsp1,&
&                enl0,enl1,epaw1,etotal,evar,evdw,exc1,ipert,natom,optene)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_etot'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: berryopt,ipert,natom,optene
 real(dp),intent(in) :: eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1
 real(dp),intent(in) :: efrx2,ehart1,eii,ek0,ek1,eloc0,elpsp1,enl0,enl1,epaw1
 real(dp),intent(in) :: evdw,exc1
 real(dp),intent(inout) :: elast
 real(dp),intent(out) :: deltae,etotal,evar

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *********************************************************************

 if (optene==1) then
   message = 'Double-counting scheme not yet allowed !'
   MSG_BUG(message)
 end if

 if (optene>-1) then

!  Compute 2nd-order variational energy by direct scheme
   if (optene==0) then

!    Atomic displ. perturbation
     if ( ipert>=1 .and. ipert<=natom  ) then
       evar=ek0+edocc+eeig0+eloc0+elpsp1+ehart1+exc1+enl0+enl1+epaw1

     else if (ipert==natom+1) then
       evar=ek0+edocc+eeig0+eloc0+ek1+ehart1+exc1+enl0+enl1

     else if (ipert==natom+10 .or. ipert==natom+11) then
       evar=ek0+edocc+eeig0+eloc0+enl0+ek1 ! here ek1 contains a lot of contributions

!      For ipert==natom+2, some contributions vanishes, noticeably ek1
     else if (ipert==natom+2) then
       evar=ek0+edocc+eeig0+eloc0+ek1+ehart1+exc1+enl0+enl1+epaw1

!      All terms enter for strain perturbation
     else if ( ipert==natom+3 .or. ipert==natom+4 ) then
       evar=ek0+ek1+edocc+eeig0+eloc0+elpsp1+ehart1+exc1+enl0+enl1+epaw1
     end if
   end if

!  Compute energy residual
   deltae=evar-elast
   elast=evar

!  Compute 2nd-order total energy by direct scheme
   if (optene==0) then
     if (berryopt==4 .or. berryopt==6 .or. berryopt==7 .or. berryopt==14 .or. berryopt==16 .or. berryopt==17) then
       if (ipert<=natom) then
         etotal=evar+eew+evdw+eii+efrhar+efrkin+efrloc+efrnl+efrx1+efrx2+two*eberry
       else if (ipert==natom+2) then
         etotal=half*evar+eew+evdw+eii+efrhar+efrkin+efrloc+efrnl+efrx1+efrx2+two*eberry
       end if
     else
       if (ipert/=natom+10 .and. ipert/=natom+11) then
         etotal=evar+eew+evdw+eii+efrhar+efrkin+efrloc+efrnl+efrx1+efrx2
       else
         etotal=evar ! For 2nd order sternheimer equations, the total (4th order) energy is not used (yet)
       end if
     end if
   end if

 end if

end subroutine dfpt_etot
!!***
