!{\src2tex{textfont=tt}}
!!****f* ABINIT/echo_xc_name
!! NAME
!! echo_xc_name
!!
!! FUNCTION
!!  Write to log and output the xc functional which will be used for this dataset
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (MJV,CE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ixc = internal code for xc functional
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine echo_xc_name (ixc)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use libxc_functionals

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'echo_xc_name'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
 integer, intent(in) :: ixc

!Local variables -------------------------
 integer :: l_citation
 character(len=500) :: message, citation

! *********************************************************************

 message =''
 citation =''

!normal case (not libxc)
 if (ixc >= 0) then

   select case (ixc)
   case (0)
     message = 'No xc applied (usually for testing) - ixc=0'
     citation = ''
!      LDA,LSD
   case (1) 
     message = 'LDA: new Teter (4/93) with spin-polarized option - ixc=1'
     citation = 'S. Goedecker, M. Teter, J. Huetter, PRB 54, 1703 (1996)'
   case (2) 
     message = 'LDA: Perdew-Zunger-Ceperley-Alder - ixc=2'
     citation = 'J.P.Perdew and A.Zunger, PRB 23, 5048 (1981) '
   case (3) 
     message = 'LDA: old Teter (4/91) fit to Ceperley-Alder data - ixc=3'
     citation = ''
   case (4) 
     message = 'LDA: Wigner - ixc=4'
     citation = 'E.P.Wigner, Trans. Faraday Soc. 34, 67 (1938)'
   case (5) 
     message = 'LDA: Hedin-Lundqvist - ixc=5'
     citation = 'L.Hedin and B.I.Lundqvist, J. Phys. C4, 2064 (1971)'
   case (6) 
     message = 'LDA: "X-alpha" xc - ixc=6'
     citation = 'Slater J. C., Phys. Rev. 81, 385 (1951)'
   case (7) 
     message = 'LDA: Perdew-Wang 92 LSD fit to Ceperley-Alder data - ixc=7'
     citation = 'J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)'
   case (8) 
     message = 'LDA: Perdew-Wang 92 LSD , exchange-only - ixc=8'
     citation = 'J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)'
   case (9) 
     message = 'LDA: Perdew-Wang 92 Ex+Ec_RPA  energy - ixc=9'
     citation = 'J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)'
   case (10) 
     message = 'LDA: RPA LSD energy (only the energy !!) - ixc=10'
     citation = ''
!      GGA
   case (11) 
     message = 'GGA: Perdew-Burke-Ernzerhof functional - ixc=11'
     citation = 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)'
   case (12) 
     message = 'GGA: x-only Perdew-Burke-Ernzerhof functional - ixc=12'
     citation = 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)'
   case (13) 
     message = 'GGA: LDA (ixc==7) energy, and the xc _potential_ is given by van Leeuwen-Baerends GGA - ixc=13'
     citation = 'R. van Leeuwen and E. J. Baerends PRA 49, 2421 (1994)'
   case (14) 
     message = 'GGA: revPBE functional - ixc=14'
     citation = 'Zhang and Yang, PRL 80, 890 (1998)'
   case (15) 
     message = 'GGA: RPBE functional - ixc=15'
     citation = 'Hammer, L. B. Hansen, and J. K. Norskov, PRB 59, 7413 (1999)'
   case (16) 
     message = 'GGA: HCTH93 functional - ixc=16'
     citation = 'F.A. Hamprecht, A.J. Cohen, D.J. Tozer, N.C. Handy, JCP 109, 6264 (1998)'
   case (17)
     message = 'GGA: HCTH120 functional - ixc=17'
     citation = 'A.D. Boese, N.L. Doltsinis, N.C. Handy, and M. Sprik, JCP 112, 1670 (1998)'
   case (23) 
     message = 'GGA: Wu Cohen functional - ixc=23'
     citation = 'Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)'
   case (24) 
     message = 'GGA: C09x exchange functional - ixc=24'
     citation = 'Valentino R. Cooper, PRB 81, 161104(R) (2010)'
   case (26)
     message = 'GGA: HCTH147 functional - ixc=26'
     citation = 'A.D. Boese, N.L. Doltsinis, N.C. Handy, and M. Sprik, JCP 112, 1670 (1998)'
   case (27)
     message = 'GGA: HCTH407 functional - ixc=27'
     citation = 'A.D. Boese, and N.C. Handy, JCP 114, 5497 (2001)'
!      Fermi-Amaldi
   case (20) 
     message = 'Fermi-Amaldi correction - ixc=20'
     citation = ''
   case (21) 
     message = 'Fermi-Amaldi correction with LDA(ixc=1) kernel - ixc=21'
     citation = ''
   case (22) 
     message = 'Fermi-Amaldi correction with hybrid BPG kernel - ixc=22'
     citation = ''
   case (31)
     message = 'Meta-GGA fake1 - ixc=31'
     citation = ''
   case (32)
     message = 'Meta-GGA fake2 - ixc=32'
     citation = ''
   case (33)
     message = 'Meta-GGA fake3 - ixc=33'
     citation = ''
   case (34)
     message = 'Meta-GGA fake4 - ixc=34'
     citation = ''
   case (40)
     message = 'Hartree-Fock with mixing coefficient alpha=1'
     citation = ''
   case (41)
     message = 'PBE0 with alpha=0.25'
     citation = ''
   case (42)
     message = 'modified PBE0 with alpha=0.33'
     citation = ''
   case (50)
     message = 'LDA at finite T Ichimaru-Iyetomy-Tanaka - ixc=50'
     citation = 'Ichimaru S., Iyetomi H., Tanaka S., Phys. Rep. 149, 91-205 (1987) '
   case default 
     write(message,'(a,i0)')" echo_xc_name does not know how to handle ixc = ",ixc
     MSG_WARNING(message)
   end select

   message = " Exchange-correlation functional for the present dataset will be:" // ch10 &
&   // "  " // trim(message)

   l_citation=len_trim(citation)
   citation = " Citation for XC functional:" // ch10 // "  " // trim(citation)

   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   if(l_citation/=0)then
     call wrtout(ab_out,citation,'COLL')
     call wrtout(std_out,citation,'COLL')
   end if

   message =' ' 
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

 end if ! end libxc if

end subroutine echo_xc_name
!!***
