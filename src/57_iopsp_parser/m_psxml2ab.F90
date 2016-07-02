!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_psxml2ab
!! NAME
!! m_psxml2ab
!!
!! FUNCTION
!!  From a SIESTA XML format pseudopotential file 
!!  convert to abinit internal datastructures for pspheader.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2016 ABINIT group (MJV).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! psxml  = pseudopotential data structure
!!
!! OUTPUT
!! psphead = psp information structure 
!!
!! PARENTS
!!      inpspheads,pspatm_abinit
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_psxml2ab

implicit none

private

#ifdef HAVE_TRIO_PSML
public :: psxml2abheader
!public :: psxml2abfull

CONTAINS


subroutine psxml2abheader(psxmlfile, psphead, iwrite)

 use defs_basis
 use defs_datatypes
 use m_profiling_abi
 use m_errors
 use m_psml

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psxml2abheader'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: iwrite
 character(len=fnlen), intent(in) :: psxmlfile
 type(pspheader_type),intent(inout) :: psphead


!Local variables-------------------------------
!scalars
 integer :: il,lmm,ii
 integer :: iproj, ishell, ll, ll_previous
 integer,parameter :: n1xccc_default=2501
 character(len=500) :: message, frmt_str
 type(ps_t) :: psxml
!arrays
 integer, allocatable :: idx_sr(:), idx_so(:)
 real(dp),allocatable :: zeld(:),zelu(:)

! *********************************************************************

 call ps_destroy(psxml)
 call psml_reader(psxmlfile, psxml, debug=.true.)


 psphead%znuclpsp  = ps_AtomicNumber(psxml)
 psphead%zionpsp   = ps_Zpseudo(psxml)

 psphead%pspcod = 9

! impose libxc coding for pspxc
 psphead%pspxc = 0
 do ii=1, ps_NLibxcFunctionals(psxml)
   psphead%pspxc = ps_LibxcId(psxml,ii) + psphead%pspxc*1000
 end do
 psphead%pspxc = -psphead%pspxc

 if (iwrite == 1) then
   write (message,'(a,a)') 'ps_PseudoFlavor ', ps_PseudoFlavor(psxml)
!  call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write (message,'(a,I5)') 'ps_NLibxcFunctionals ', ps_NLibxcFunctionals(psxml)
!  call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write (message,'(a,I5)') 'ps_NValenceShells ', ps_NValenceShells(psxml)
!  call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

 ABI_ALLOCATE(zeld, (ps_NValenceShells(psxml)))
 ABI_ALLOCATE(zelu, (ps_NValenceShells(psxml)))
 zeld = zero
 zelu = zero
 frmt_str="(a"
 do ishell = 1, ps_NValenceShells(psxml)
   if (iwrite == 1) then
     write (message,'(a,I5)') 'ps_ValenceShellN ', ps_ValenceShellN(psxml,ishell)
!    call wrtout(ab_out,  message,'COLL')
     call wrtout(std_out,  message,'COLL')
     write (message,'(a,I5)') 'ps_ValenceShellL ', ps_ValenceShellL(psxml,ishell)
!    call wrtout(ab_out,  message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if
   if (ps_IsSpinPolarized(psxml)) then
     zelu(ishell) = ps_ValenceShellOccupation(psxml,ishell, channel='u')
     zeld(ishell) = ps_ValenceShellOccupation(psxml,ishell, channel='d')
   else
     zeld(ishell) = ps_ValenceShellOccupation(psxml,ishell)
   end if
   frmt_str = trim(frmt_str) // ", F10.3"
 end do
 frmt_str = trim(frmt_str) // ")"
 if (iwrite == 1) then
   write (message,frmt_str) 'zeld ', zeld
!  call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write (message,frmt_str) 'zelu ', zelu
!  call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

!    Find the number of projectors per angular momentum shell
 psphead%nproj(:)=0
 if (ps_Number_Of_Projectors(psxml, SET_NONREL) > 0) then
   idx_sr = ps_Projector_Indexes(psxml, SET_NONREL)
   do iproj = 1, ps_Number_Of_Projectors(psxml, SET_NONREL)
     if (iwrite == 1) then
       write (message,'(a,2I5)') 'iproj, idx for nonrel ', iproj, idx_sr(iproj)
!      call wrtout(ab_out,  message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end if
     il = ps_Projector_L(psxml, idx_sr(iproj))
     psphead%nproj(il) = psphead%nproj(il) + 1
   end do
 else if (ps_Number_Of_Projectors(psxml, SET_SREL) > 0) then
   idx_sr = ps_Projector_Indexes(psxml, SET_SREL)
   do iproj = 1, ps_Number_Of_Projectors(psxml, SET_SREL)
     if (iwrite == 1) then
       write (message,'(a,2I5)') 'iproj, idx for srel ', iproj, idx_sr(iproj)
!      call wrtout(ab_out,  message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end if
     il = ps_Projector_L(psxml, idx_sr(iproj))
     psphead%nproj(il) = psphead%nproj(il) + 1
   end do
 else
   MSG_BUG('Your psml potential should have either scalar- or non- relativistic projectors') 
 end if

 psphead%nprojso(:)=0
 idx_so = ps_Projector_Indexes(psxml, SET_SO)
 do iproj = 1, ps_Number_of_Projectors(psxml, SET_SO)
   if (iwrite == 1) then
     write (message,'(a,2I5)') 'iproj, idx for soc ', iproj, idx_so(iproj)
!    call wrtout(ab_out,  message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if
   il = ps_Projector_L(psxml, idx_so(iproj))
   psphead%nprojso(il) = psphead%nprojso(il) + 1
 end do
 if (iwrite == 1) then
   write (message,'(a,5I5)') 'nproj ', psphead%nproj
!  call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write (message,'(a,5I5)') 'nprojso ', psphead%nprojso
!  call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

 if( ps_HasCoreCorrections(psxml)) then
   psphead%xccc  = n1xccc_default
 else
   psphead%xccc  = 0
 end if

 psphead%pspso = 0
 if (sum(abs(psphead%nprojso(:))) > 0) psphead%pspso = 2


 psphead%lmax = 0
 if (iwrite == 1) then
   write (message,'(a,I5)') ' ps_Number_of_Projectors not relativistic ',&
&        ps_Number_of_Projectors(psxml, SET_NONREL)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write (message,'(a,I5)') ' ps_Number_of_Projectors scalar relativistic ', &
&        ps_Number_of_Projectors(psxml, SET_SREL)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if
 ll_previous=-1
 do iproj = 1, max(ps_Number_of_Projectors(psxml, SET_SREL), ps_Number_of_Projectors(psxml, SET_NONREL))
   ll = ps_Projector_L(psxml, idx_sr(iproj))
   if (iwrite == 1) then
     if(ll/=ll_previous)then
       write (message,'(a,I5)') ' ps_Projector_L ', ll
       call wrtout(ab_out,  message,'COLL')
       call wrtout(std_out,  message,'COLL')
       ll_previous=ll
     endif
     write (message,'(a,E20.10)') ' ps_Projector_Ekb ', ps_Projector_Ekb(psxml, idx_sr(iproj))
     call wrtout(ab_out,  message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if
   psphead%lmax = max( psphead%lmax, ll)
 end do

 if (iwrite == 1) then
   write (message,'(a,I5)') ' ps_Number_of_Projectors SOC ', ps_Number_of_Projectors(psxml, SET_SO)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if
 ll_previous=-1
 do iproj = 1, ps_Number_of_Projectors(psxml, SET_SO)
   ll = ps_Projector_L(psxml, idx_so(iproj))
   if (iwrite == 1) then
     if(ll/=ll_previous)then
       write (message,'(a,I5)') ' ps_Projector_L ', ll
       call wrtout(ab_out,  message,'COLL')
       call wrtout(std_out,  message,'COLL')
       ll_previous=ll
     endif
     write (message,'(a,E20.10)') ' ps_Projector_Ekb ', ps_Projector_Ekb(psxml, idx_so(iproj))
     call wrtout(ab_out,  message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if
   psphead%lmax = max( psphead%lmax, ll)
 end do

 lmm     = max( ps_Number_of_Projectors(psxml, SET_SREL), ps_Number_of_Projectors(psxml, SET_NONREL))
 lmm     = max( lmm, ps_Number_of_Projectors(psxml, SET_SO))

 write(std_out,'(a,f5.1,a,i4,a,i4)' ) '  read the values zionpsp=',&
& psphead%zionpsp,' , pspcod=',psphead%pspcod,' , lmax=',psphead%lmax

 if ( iwrite .eq. 1 ) then

   write(message,'(a,a)') &
&   ' psxml2ab: Atomic Label:                      ', &
&  ps_AtomicLabel(psxml)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,f12.5)') &
&   ' psxml2ab: Atomic Number:                     ', &
&   ps_AtomicNumber(psxml)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,f12.5)') &
&   ' psxml2ab: Valence charge:                    ', &
&   ps_ZPseudo(psxml)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,a)') &
&   ' psxml2ab: Pseudopotential generator code :    ', &
&   ps_Creator(psxml)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,a)') &
&   ' psxml2ab: Date of pseudopotential generation: ', &
&   ps_Date(psxml)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,a)') &
&   ' psxml2ab: Pseudopotential flavor:             ', &
&   ps_PseudoFlavor(psxml)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   do ii = 1, ps_NLibxcFunctionals(psxml)
     write(message,'(a,I4,2a)') &
&     ' psxml2ab: Exchange-correlation functional ',ii,' :    ', &
&     ps_LibxcName(psxml,ii)
     call wrtout(ab_out,  message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end do

   write(message,'(a,l3)') &
&   ' psxml2ab: Relativistically generated pseudopotential (not necessarily SOC!):   ', &
&   ps_Relativity(psxml)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,l3)') &
&   ' psxml2ab: Spin-polarized pseudopotential:     ', &
&   ps_IsSpinPolarized(psxml)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   select case(ps_HasCoreCorrections(psxml))
     case(.true.)
       write(message, '(a)' ) &
&       '- psxml2ab: XC core correction read in from XML file.'
     case(.false.)
       write(message, '(a)' ) &
&       '- psxml2ab: No core corrections.'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
   end select
 end if

 ABI_DEALLOCATE(zeld)
 ABI_DEALLOCATE(zelu)

 if (allocated(idx_sr)) then
   ABI_FREE_NOCOUNT(idx_sr)
 end if
 if (allocated(idx_so)) then
   ABI_FREE_NOCOUNT(idx_so)
 end if

 call ps_destroy(psxml)

end subroutine psxml2abheader
!!***
! end test on compiling with LIBPSML enabled
#endif   

!{\src2tex{textfont=tt}}
!!****f* ABINIT/psxml2abfull
!! NAME
!! psxml2abfull
!!
!! FUNCTION
!!  read in all data from psml file. Call header reader first, then local potential and NL projectors
!!  
!!
!! COPYRIGHT
!! Copyright (C) 2005-2016 ABINIT group (MJV).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   str = string with error message
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE




!!***
end module m_psxml2ab

!{\src2tex{textfont=tt}}
!!****f* ABINIT/psml_die
!! NAME
!! psml_die
!!
!! FUNCTION
!!  auxiliary die function for calling libPSML. Needed at link time
!!  allows calling software to decide how fatal the PSML die call actually is.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2016 ABINIT group (MJV).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   str = string with error message
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine psml_die(str)
  use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psml_die'
!End of the abilint section

  implicit none

  character(len=*), intent(in) :: str

  MSG_BUG(str)

end subroutine psml_die
!!***

