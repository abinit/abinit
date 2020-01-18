!!****m* ABINIT/m_psxml2ab
!! NAME
!! m_psxml2ab
!!
!! FUNCTION
!!  From a SIESTA XML format pseudopotential file
!!  convert to abinit internal datastructures for pspheader.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2020 ABINIT group (MJV).
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
!! atmsymb = atomic symbol
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

 use defs_basis
 use defs_datatypes
 use m_abicore
 use m_errors
#ifdef HAVE_PSML
 use m_psml
 use m_psml_api
#endif

 use defs_datatypes, only : pspheader_type
 use m_fstrings,     only : yesno

implicit none

private

#ifdef HAVE_PSML
public :: psxml2abheader
!public :: psxml2abfull

CONTAINS

subroutine psxml2abheader(psxmlfile, psphead, atmsymb, creator, iwrite)

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: iwrite
 character(len=fnlen), intent(in) :: psxmlfile
 type(pspheader_type),intent(inout) :: psphead
 character(len=3), intent(out) :: atmsymb
 character(len=30), intent(out) :: creator

!Local variables-------------------------------
!scalars
 character(len=500) :: flavor, frmt_str, label, message, relat, xc_name
 character(len=1) :: g1,g2
 integer :: dd,dm,dy
 integer :: il,lmm,ii
 integer :: nxc
 integer :: ishell, nvshells, shell_l, shell_n
 integer :: iproj, ll, ll_previous, nprojs, nprojsr, nprojso
 integer,parameter :: n1xccc_default=2501
 logical :: has_nlcc, has_spin
 real(dp) :: ekb
 type(ps_t) :: psxml
!arrays
 integer, allocatable :: idx_sr(:), idx_so(:)
 real(dp),allocatable :: zeld(:),zelu(:)

! *********************************************************************

 call ps_destroy(psxml)
 call psml_reader(psxmlfile, psxml, debug=.true.)

 psphead%pspdat = 0
 call ps_Provenance_Get(psxml, 1, creator=creator, date=message)
 read (message(1:10), '(I4,A1,I2,A1,I2)', err=10) &
&  dy, g1, dm, g2, dd
 psphead%pspdat = MODULO(dy,100) * 10000 + dm * 100 + dd

10 continue

 call ps_PseudoAtomSpec_Get(psxml, &
&  atomic_symbol=atmsymb, atomic_label=label, &
&  atomic_number=psphead%znuclpsp, z_pseudo=psphead%zionpsp, &
&  pseudo_flavor=flavor, relativity=relat, &
&  spin_dft=has_spin, core_corrections=has_nlcc)

 psphead%pspcod = 9

! impose libxc coding for pspxc
 psphead%pspxc = 0
 call ps_ExchangeCorrelation_Get(psxml, n_libxc_functionals=nxc)
 do ii=1, min(2,nxc)
   call ps_LibxcFunctional_Get(psxml, ii, code=dd)
   psphead%pspxc = dd + psphead%pspxc*1000
 end do
 psphead%pspxc = -psphead%pspxc

 call ps_ValenceConfiguration_Get(psxml, nshells=nvshells)

 if (iwrite == 1) then
   write (message,'(a,a)') '- psxml2ab: ps_PseudoFlavor ', trim(message)
!  call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write (message,'(a,I5)') '- psxml2ab: ps_NLibxcFunctionals ', nxc
!  call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write (message,'(a,I5)') '- psxml2ab: ps_NValenceShells ', nvshells
!  call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

 ABI_ALLOCATE(zeld, (nvshells))
 ABI_ALLOCATE(zelu, (nvshells))
 zeld = zero
 zelu = zero
 frmt_str="(a"
 do ishell = 1, nvshells
   if (has_spin) then
     call ps_ValenceShell_Get(psxml, ishell, n=shell_n, l=shell_l, &
&    occ_up=zelu(ishell), occ_down=zeld(ishell))
   else
     call ps_ValenceShell_Get(psxml, ishell, n=shell_n, l=shell_l, &
&    occupation=zeld(ishell))
   end if
   if (iwrite == 1) then
     write (message,'(a,I5)') '- psxml2ab: ps_ValenceShellN ', shell_n
!    call wrtout(ab_out,  message,'COLL')
     call wrtout(std_out,  message,'COLL')
     write (message,'(a,I5)') '- psxml2ab: ps_ValenceShellL ', shell_l
!    call wrtout(ab_out,  message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if
   frmt_str = trim(frmt_str) // ", F10.3"
 end do
 frmt_str = trim(frmt_str) // ")"
 if (iwrite == 1) then
   write (message,frmt_str) '- psxml2ab: zeld ', zeld
!  call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write (message,frmt_str) '- psxml2ab: zelu ', zelu
!  call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

!    Find the number of projectors per angular momentum shell
 nprojs = 0
 nprojsr = 0
 nprojso = 0
 psphead%nproj(:)=0
 call ps_NonlocalProjectors_Filter(psxml, set=SET_NONREL, number=nprojs)
 if (nprojs > 0) then
   call ps_NonlocalProjectors_Filter(psxml, set=SET_NONREL, indexes=idx_sr)
   do iproj = 1, nprojs
     if (iwrite == 1) then
       write (message,'(a,2I5)') '- psxml2ab: iproj, idx for nonrel ', iproj, idx_sr(iproj)
!      call wrtout(ab_out,  message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end if
     call ps_Projector_Get(psxml, idx_sr(iproj), l=il)
     psphead%nproj(il) = psphead%nproj(il) + 1
   end do
 else
   call ps_NonlocalProjectors_Filter(psxml, set=SET_SREL, number=nprojsr)
   if (nprojsr > 0) then
     call ps_NonlocalProjectors_Filter(psxml, set=SET_SREL, indexes=idx_sr)
     do iproj = 1, nprojsr
       if (iwrite == 1) then
         write (message,'(a,2I5)') '- psxml2ab: iproj, idx for srel ', iproj, idx_sr(iproj)
!        call wrtout(ab_out,  message,'COLL')
         call wrtout(std_out,  message,'COLL')
       end if
       call ps_Projector_Get(psxml, idx_sr(iproj), l=il)
       psphead%nproj(il) = psphead%nproj(il) + 1
     end do
   else
     MSG_BUG('Your psml potential should have either scalar- or non- relativistic projectors')
   end if
 end if

 psphead%nprojso(:)=0
 call ps_NonlocalProjectors_Filter(psxml, set=SET_SO, number=nprojso, &
&  indexes=idx_so)
 do iproj = 1, nprojso
   if (iwrite == 1) then
     write (message,'(a,2I5)') '- psxml2ab: iproj, idx for soc ', iproj, idx_so(iproj)
!    call wrtout(ab_out,  message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if
   call ps_Projector_Get(psxml, idx_so(iproj), l=il)
   psphead%nprojso(il) = psphead%nprojso(il) + 1
 end do
 if (iwrite == 1) then
   write (message,'(a,5I5)') '- psxml2ab: nproj ', psphead%nproj
!  call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write (message,'(a,5I5)') '- psxml2ab: nprojso ', psphead%nprojso
!  call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

 if( has_nlcc) then
   psphead%xccc  = n1xccc_default
 else
   psphead%xccc  = 0
 end if

 psphead%pspso = 0
 if (sum(abs(psphead%nprojso(:))) > 0) psphead%pspso = 2


 psphead%lmax = 0
 if (iwrite == 1) then
   write (message,'(a,I5)') '- psxml2ab: ps_Number_of_Projectors not relativistic ',&
&        nprojs
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write (message,'(a,I5)') '- psxml2ab: ps_Number_of_Projectors scalar relativistic ', &
&        nprojsr
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if
 ll_previous=-1
 do iproj = 1, max(nprojs, nprojsr)
   call ps_Projector_Get(psxml, idx_sr(iproj), l=ll)
   if (iwrite == 1) then
     if(ll/=ll_previous)then
       write (message,'(a,I5)') '- psxml2ab: ps_Projector_L ', ll
       call wrtout(ab_out,  message,'COLL')
       call wrtout(std_out,  message,'COLL')
       ll_previous=ll
     endif
     call ps_Projector_Get(psxml, idx_sr(iproj), ekb=ekb)
     write (message,'(a,E20.10)') '- psxml2ab: ps_Projector_Ekb ', ekb
     call wrtout(ab_out,  message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if
   psphead%lmax = max( psphead%lmax, ll)
 end do

 if (iwrite == 1) then
   write (message,'(a,I5)') '- psxml2ab: ps_Number_of_Projectors SOC ', nprojso
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if
 ll_previous=-1
 do iproj = 1, nprojso
   call ps_Projector_Get(psxml, idx_so(iproj), l=ll)
   if (iwrite == 1) then
     if(ll/=ll_previous)then
       write (message,'(a,I5)') '- psxml2ab: ps_Projector_L ', ll
       call wrtout(ab_out,  message,'COLL')
       call wrtout(std_out,  message,'COLL')
       ll_previous=ll
     endif
     call ps_Projector_Get(psxml, idx_so(iproj), ekb=ekb)
     write (message,'(a,E20.10)') '- psxml2ab: ps_Projector_Ekb ', ekb
     call wrtout(ab_out,  message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if
   psphead%lmax = max( psphead%lmax, ll)
 end do

 lmm     = max( nprojs, nprojsr)
 lmm     = max( lmm, nprojso)

 write(std_out,'(a,f5.1,a,i4,a,i4)' ) '- psxml2ab:  read the values zionpsp=',&
& psphead%zionpsp,' , pspcod=',psphead%pspcod,' , lmax=',psphead%lmax

 if ( iwrite .eq. 1 ) then

   write(message,'(a,a)') &
&   '- psxml2ab: Atomic Label:                      ', &
&  trim(label)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,f12.5)') &
&   '- psxml2ab: Atomic Number:                     ', &
&   psphead%znuclpsp
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,f12.5)') &
&   '- psxml2ab: Valence charge:                    ', &
&   psphead%zionpsp
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,a)') &
&   '- psxml2ab: Pseudopotential generator code :    ', &
&   creator
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,i8)') &
&   '- psxml2ab: Date of pseudopotential generation: ', &
&   psphead%pspdat
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(a,a)') &
&   '- psxml2ab: Pseudopotential flavor:             ', &
&   trim(flavor)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   do ii = 1, nxc
     call ps_LibxcFunctional_Get(psxml, ii, name=xc_name, code=dd)
     write(message,'(a,I4,2a," (",I4,")")') &
&     '- psxml2ab: Exchange-correlation functional ',ii,' :    ', &
&     trim(xc_name), dd
     call wrtout(ab_out,  message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end do

   write(message,'(2a)') &
&   '- psxml2ab: Relativistically generated pseudopotential (not necessarily SOC!):   ', &
&   trim(relat)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   write(message,'(2a)') &
&   '- psxml2ab: Spin-polarized pseudopotential:     ', &
&   yesno(has_spin)
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')

   select case(has_nlcc)
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

end module m_psxml2ab
!!***

!!****f* ABINIT/psml_die
!! NAME
!! psml_die
!!
!! FUNCTION
!!  auxiliary die function for calling libPSML. Needed at link time
!!  allows calling software to decide how fatal the PSML die call actually is.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2020 ABINIT group (MJV).
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
  implicit none

  character(len=*), intent(in) :: str

  MSG_BUG(str)

end subroutine psml_die
!!***

