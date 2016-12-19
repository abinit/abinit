!{\src2tex{textfont=tt}}
!!****f* ABINIT/printbxsf
!! NAME
!! printbxsf
!!
!! FUNCTION
!!  Print band structure energies in XCrysDen format.
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2016 ABINIT group (MVerstraete,MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  eigen(mband,nkpt,nsppol) = eigenvalues in hartree
!!  ewind = energy window around the fermi level.
!!          if ewind /= 0 ==> a band is considered in the plot of FSurf
!!                            only if it is inside [ ef-ewind, ef+ewind ] for some k point
!!          if ewind == 0 ==> all bands will be keept in the _BXSF file
!!  fermie = Fermi energy (Hartree)
!!  gprimd(3,3) = dimensional primitive translations for reciprocal space (bohr^-1)
!!  kptrlatt(3,3) = reciprocal of lattice vectors for full kpoint grid
!!  mband = maximum number of bands
!!  nsppol = 1 for unpolarized, 2 for spin-polarized
!!  shiftk(3,nshiftk) =shift vector for k point grid
!!  fname = filename for the fortran file
!!  symafm(nsym)=(Anti)ferromagnetic symmetries.
!!  use_afm=.TRUE. if (anti)ferromagnetic symmetries are used.
!!
!! OUTPUT
!!  ierr=Status error.
!!  BXSF file.
!!
!! PARENTS
!!      m_ebands,m_ifc
!!
!! CHILDREN
!!      destroy_kptrank,get_rank_1kpt,listkk,mkkptrank
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine printbxsf(eigen,ewind,fermie,gprimd,kptrlatt,mband,&
& nkptirred,kptirred,nsym,use_afm,symrec,symafm,use_tr,nsppol,shiftk,nshiftk,fname,ierr)

 use defs_basis
 use m_errors
 use m_kptrank
 use m_profiling_abi

 use m_io_tools, only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'printbxsf'
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkptirred,nshiftk,nsppol,nsym
 integer,intent(out) :: ierr
 real(dp),intent(in) :: ewind,fermie
 logical,intent(in) :: use_afm,use_tr
 character(len=*),intent(in) :: fname
!arrays
 integer,intent(in) :: kptrlatt(3,3),symafm(nsym),symrec(3,3,nsym)
 real(dp),intent(in) :: eigen(mband,nkptirred,nsppol),gprimd(3,3)
 real(dp),intent(in) :: kptirred(3,nkptirred),shiftk(3,nshiftk)

!Local variables-------------------------------
!scalars
 integer :: iband,ik1,ik2,ik3,ikgrid,ikpt,indx
 integer :: isppol,isym,maxband,minband,nk1,nk2,nk3,nkptfull,ubxsf,timrev
 integer :: symkptrank, nsymfm, isymfm
 real(dp) :: ene,dksqmax
 character(len=500) :: msg
 type(kptrank_type) :: kptrank_t
!arrays
 integer :: indkk_kq(1,6)
 integer,allocatable :: fulltoirred(:),symrecfm(:,:,:)
 real(dp) :: kptgrid(3),gmet(3,3)

! *************************************************************************

 ierr = 0

! Error if klatt is no simple orthogonal lattice (in red space)
! for generalization to MP grids, need new version of XCrysDen

 if ( kptrlatt(1,2)/=0 .or. kptrlatt(1,3)/=0 .or. kptrlatt(2,1)/=0 .or. &
&     kptrlatt(2,3)/=0 .or. kptrlatt(3,1)/=0 .or. kptrlatt(3,2)/=0 ) then
   write(msg,'(3a)')&
&   'kptrlatt should be diagonal, for the FS calculation ',ch10,&
&   'Action: use an orthogonal k-grid for the GS calculation '
   MSG_COMMENT(msg)
   ierr = ierr + 1
 end if

! Error if there are not at least 2 kpts in each direction:
! kptrank will fail for the intermediate points below
 if ( abs(kptrlatt(1,1))<2 .or. abs(kptrlatt(2,2))<2 .or. abs(kptrlatt(3,3))<2) then
   write(msg,'(3a)')&
&   'You need at least 2 points in each direction in k space to output BXSF files ',ch10,&
&   'Action: use an augmented k-grid for the GS calculation (at least 2x2x2) '
   MSG_COMMENT(msg)
   ierr = ierr + 1
 end if

 if (ANY(ABS(shiftk) > tol10)) then
   write(msg,'(3a)')&
&   'Origin of the k-grid should be (0,0,0) for the FS calculation ',ch10,&
&   'Action: use a non-shifted k-grid for the GS calculation. Returning '
   MSG_COMMENT(msg)
   ierr = ierr + 1
 end if

 if (ierr /= 0) return

 ! Compute reciprocal space metric.
 gmet = MATMUL(TRANSPOSE(gprimd),gprimd)

 if (use_afm) then
   nsymfm = 0
   do isym = 1, nsym
     if (symafm(isym) == 1) nsymfm = nsymfm+1
   end do
   ABI_MALLOC(symrecfm,(3,3,nsymfm))
   isymfm = 0
   do isym = 1, nsym
     if (symafm(isym) == 1) then
       isymfm = isymfm + 1
       symrecfm(:,:,isymfm) = symrec(:,:,isym)
     end if
   end do
 else
   nsymfm = nsym
   ABI_MALLOC(symrecfm,(3,3,nsymfm))
   symrecfm = symrec
 end if

!Xcrysden uses aperiodical data-grid
 nk1 = kptrlatt(1,1)
 nk2 = kptrlatt(2,2)
 nk3 = kptrlatt(3,3)
 nkptfull=(nk1+1)*(nk2+1)*(nk3+1)

 ABI_MALLOC(fulltoirred,(nkptfull))
 timrev=0; if (use_tr) timrev=1

 call mkkptrank (kptirred,nkptirred,kptrank_t, nsym=nsymfm, symrec=symrecfm, time_reversal=use_tr)

!Xcrysden employs the C-ordering for the Fermi Surface.
 ikgrid=0
 do ik1=0,nk1
   do ik2=0,nk2
     do ik3=0,nk3

       ikgrid=ikgrid+1
       kptgrid(1)=DBLE(ik1)/kptrlatt(1,1)
       kptgrid(2)=DBLE(ik2)/kptrlatt(2,2)
       kptgrid(3)=DBLE(ik3)/kptrlatt(3,3)

       ! Find correspondence between the Xcrysden grid and the IBZ ===
#if 1
       call get_rank_1kpt (kptgrid, symkptrank, kptrank_t)
       fulltoirred(ikgrid) = kptrank_t%invrank(symkptrank)

       if (fulltoirred(ikgrid) < 1) then
         write(msg,'(a,3es16.8,2a,I8,2a)')&
&         'kpt = ',kptgrid,ch10,' with rank ', symkptrank, ch10,&
&         'has no symmetric among the k-points used in the GS calculation '
         ierr=ierr + 1
         MSG_WARNING(msg)
       end if

#else
       ! TODO:: symafm are not passed in a consisten way.
       call listkk(dksqmax,gmet,indkk_kq,kptirred,kptgrid,nkptirred,1,nsymfm,&
         1,symafm,symrecfm,timrev,use_symrec=.True.)

       if (dksqmax > tol12) then
         write(msg, '(7a,es16.6,4a)' )&
          'The WFK file cannot be used to start thee present calculation ',ch10,&
          'It was asked that the wavefunctions be accurate, but',ch10,&
          'at least one of the k points could not be generated from a symmetrical one.',ch10,&
          'dksqmax=',dksqmax,ch10,&
          'Action: check your WFK file and k point input variables',ch10,&
          '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
         MSG_WARNING(msg)
         ierr = ierr + 1
       end if
       fulltoirred(ikgrid) = indkk_kq(1,1)

       !ikq_ibz = indkk_kq(1,1); isym_kq = indkk_kq(1,2)
       !trev_kq = indkk_kq(1, 6); g0_kq = indkk_kq(1, 3:5)
       !isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
       !kq_ibz = ebands%kptns(:,ikq_ibz)
#endif

     end do !ik1
   end do !ik2
 end do !ik3

 call destroy_kptrank(kptrank_t)

 if (ierr/=0) then
   ABI_FREE(fulltoirred)
   MSG_ERROR("Bug")
   RETURN
 end if

 if (abs(ewind) < tol12 ) then ! Keep all bands.
   minband=1
   maxband=mband
 else ! Select a subset of bands.
   minband = mband
   maxband = 0
   ene=abs(ewind)
   do isppol=1,nsppol
     do iband=1,mband
       if(minval(eigen(iband,:,isppol))-fermie < -ene) then
         minband = iband
       end if
     end do
     do iband=mband,1,-1
       if (maxval(eigen(iband,:,isppol))-fermie > ene) then
         maxband = iband
       end if
     end do
   end do ! isppol

 end if ! abs(energy_window)

 ! Dump the results on file
 if (open_file(fname,msg,newunit=ubxsf,status='unknown',form='formatted') /=0) then
   MSG_WARNING(msg)
   ierr=ierr +1; RETURN
 end if

! Write header
 write(ubxsf,*)' BEGIN_INFO'
 write(ubxsf,*)'   #'
 write(ubxsf,*)'   # this is a Band-XCRYSDEN-Structure-File for Visualization of Fermi Surface'
 write(ubxsf,*)'   # generated by the ABINIT package'
 write(ubxsf,*)'   #'
 write(ubxsf,*)'   #  bands between ',minband,' and ',maxband
 write(ubxsf,*)'   #'
 if (nsppol == 2 ) then
   write(ubxsf,*)'   # NOTE: the first band is relative to spin-up electrons,'
   write(ubxsf,*)'   # the second band to spin-down and so on .. '
   write(ubxsf,*)'   #'
 end if
 write(ubxsf,*)'   # Launch as: xcrysden --bxsf '
 write(ubxsf,*)'   #'
 write(ubxsf,'(a,es16.8)')'   Fermi Energy: ',fermie
 write(ubxsf,*)' END_INFO'
 write(ubxsf,*)' '
 write(ubxsf,*)' BEGIN_BLOCK_BANDGRID_3D'
 write(ubxsf,*)' band_energies'
 write(ubxsf,*)' BEGIN_BANDGRID_3D'

 write(ubxsf,*)' ',(maxband-minband+1)*nsppol
 write(ubxsf,*)' ',nk1+1,nk2+1,nk3+1
 write(ubxsf,*)' ',shiftk(:,1)
!NOTE : Angstrom units are used in the BXSF format
 write(ubxsf,*)' ',gprimd(:,1)/Bohr_Ang
 write(ubxsf,*)' ',gprimd(:,2)/Bohr_Ang
 write(ubxsf,*)' ',gprimd(:,3)/Bohr_Ang

!print out data for all relevant bands and full kpt grid (redundant, yes)
!for each kpt in full zone, find equivalent irred kpt and print eigenval
 indx=0
 do iband=minband,maxband
   do isppol=1,nsppol
     write(ubxsf,*)' BAND: ',indx+minband
     write(ubxsf,'(7(es16.8))')(eigen(iband,fulltoirred(ikpt),isppol),ikpt=1,nkptfull)
     indx=indx+1
   end do
 end do

 write(ubxsf,*)'  END_BANDGRID_3D'
 write(ubxsf,*)' END_BLOCK_BANDGRID_3D'
 close (ubxsf)

 ABI_FREE(fulltoirred)
 ABI_FREE(symrecfm)

end subroutine printbxsf
!!***
