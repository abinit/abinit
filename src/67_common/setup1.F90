!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup1
!!
!! NAME
!! setup1
!!
!! FUNCTION
!! Call near top of main routine to handle setup of various arrays,
!! filenames, checking of input data, etc.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  acell(3)=length scales (bohr)
!!  amu(ntypat)=mass of each atom type
!!  ecut_eff=effective energy cutoff (hartree) for planewave basis sphere
!!  ecutc_eff=- PAW only - effective energy cutoff (hartree) for the coarse grid
!!  natom=number of atoms
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ngfftc(18)=contain all needed information about 3D FFT for the coarse grid
!!  nkpt=number of k points
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms
!!  response=0 if called by gstate, =1 if called by respfn
!!  rprim(3,3)=dimensionless real space primitive translations
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!  amass(natom)=atomic masses for each atom (in atomic units, where the electron mass is one)
!!  bantot=total number of bands at all k points
!!  gmet(3,3)=metric for reciprocal space inner products (bohr^-2)
!!  gprimd(3,3)=dimens. primitive translations for reciprocal space (bohr**-1)
!!  gsqcut_eff=Fourier cutoff on G^2 for "large sphere" of radius double
!!  gsqcutc_eff=(PAW) Fourier cutoff on G^2 for "large sphere" of radius double for the coarse FFT grid
!!   that of the basis sphere--appropriate for charge density rho(G),
!!   Hartree potential, and pseudopotentials, corresponding to ecut_eff
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol=unit cell volume (bohr^3)
!!
!! NOTES
!! SHOULD BE CLEANED !
!!
!! PARENTS
!!      gstate,nonlinear,respfn
!!
!! CHILDREN
!!      getcut,metric,mkrdim,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine setup1(acell,amass,amu,bantot,dtset,ecut_eff,ecutc_eff,gmet,&
&  gprimd,gsqcut_eff,gsqcutc_eff,natom,ngfft,ngfftc,nkpt,nsppol,&
&  response,rmet,rprim,rprimd,ucvol,usepaw)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_geometry,   only : mkrdim, metric
 use m_kg,         only : getcut

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup1'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: natom,nkpt,nsppol
 integer,intent(in) :: response,usepaw
 integer,intent(out) :: bantot
 real(dp),intent(in) :: ecut_eff,ecutc_eff
 real(dp),intent(out) :: gsqcut_eff,gsqcutc_eff,ucvol
!arrays
 integer,intent(in) :: ngfft(18),ngfftc(18)
 real(dp),intent(in) :: acell(3),amu(dtset%ntypat),rprim(3,3)
 real(dp),intent(out) :: amass(natom),gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),intent(out) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: iatom,ikpt,isppol
 real(dp) :: boxcut,boxcutc
 character(len=500) :: message
!arrays
 real(dp) :: k0(3)

! ************************************************************************

!Compute bantot
 bantot=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     bantot=bantot+dtset%nband(ikpt+(isppol-1)*nkpt)
   end do
 end do

 if(dtset%nqpt>1.or.dtset%nqpt<0) then
   write(message,'(a,i0,5a)')&
&   '  nqpt =',dtset%nqpt,' is not allowed',ch10,&
&   '  (only 0 or 1 are allowed).',ch10,&
&   '  Action : correct your input file.'
   call wrtout(ab_out,message,'COLL')
   MSG_ERROR(message)
 end if

!Compute dimensional primitive translations rprimd
 call mkrdim(acell,rprim,rprimd)

!Obtain dimensional translations in reciprocal space gprimd,
!metrics and unit cell volume, from rprimd.
!Also output rprimd, gprimd and ucvol
 call metric(gmet,gprimd,ab_out,rmet,rprimd,ucvol)

!Assign masses to each atom (for MD)
 do iatom=1,natom
   amass(iatom)=amu_emass*amu(dtset%typat(iatom))
 end do

!Get boxcut for given acell, gmet, ngfft, and ecut_eff
!(center at 000 for groundstate, center at q for respfn):
!boxcut=ratio of basis sphere diameter to fft box side
 k0(:)=0.0_dp
 if(response==1 .and. dtset%nqpt==1)then
   k0(:)=dtset%qptn(:)
   write(message, '(a)' )' setup1 : take into account q-point for computing boxcut.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if
 if (usepaw==1) then
   write(message,'(2a)') ch10,' Coarse grid specifications (used for wave-functions):'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call getcut(boxcutc,ecutc_eff,gmet,gsqcutc_eff,dtset%iboxcut,ab_out,k0,ngfftc)
   write(message,'(2a)') ch10,' Fine grid specifications (used for densities):'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call getcut(boxcut,ecut_eff,gmet,gsqcut_eff,dtset%iboxcut,ab_out,k0,ngfft)
 else
   call getcut(boxcut,ecut_eff,gmet,gsqcut_eff,dtset%iboxcut,ab_out,k0,ngfft)
   gsqcutc_eff=gsqcut_eff
 end if

!Check that boxcut>=2 if dtset%intxc=1; otherwise dtset%intxc must be set=0
 if (boxcut<2.0_dp.and.dtset%intxc==1) then
   write(message, '(a,es12.4,a,a,a,a,a)' )&
&   '  boxcut=',boxcut,' is < 2.0  => intxc must be 0;',ch10,&
&   '  Need larger ngfft to use intxc=1.',ch10,&
&   '  Action : you could increase ngfft, or decrease ecut, or put intxcn=0.'
   MSG_ERROR(message)
 end if

end subroutine setup1
!!***
