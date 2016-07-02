!{\src2tex{textfont=tt}}
!!****f* ABINIT/printvtk
!! NAME
!! printvtk
!!
!! FUNCTION
!!  Print band structure energies and velocities in VTK format.
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2016 ABINIT group (BXu)
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
!!      elphon
!!
!! CHILDREN
!!      wrap2_pmhalf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine printvtk(eigen,v_surf,ewind,fermie,gprimd,kptrlatt,mband,&
& nkptirred,kptirred,nsym,use_afm,symrec,symafm,use_tr,nsppol,shiftk,nshiftk,fname,ierr)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_io_tools,        only : open_file
 use m_numeric_tools,   only : wrap2_pmhalf

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'printvtk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkptirred,nshiftk,nsppol,nsym
 integer,intent(out) :: ierr
 real(dp),intent(in) :: ewind,fermie
 logical,intent(in) :: use_afm,use_tr
 character(len=fnlen),intent(in) :: fname
!arrays
 integer,intent(in) :: kptrlatt(3,3),symafm(nsym),symrec(3,3,nsym)
 real(dp),intent(in) :: eigen(mband,nkptirred,nsppol),gprimd(3,3)
 real(dp),intent(in) :: v_surf(mband,kptrlatt(1,1)+1,kptrlatt(2,2)+1,kptrlatt(3,3)+1,3,nsppol)
 real(dp),intent(in) :: kptirred(3,nkptirred),shiftk(3,nshiftk)

!Local variables-------------------------------
!scalars
 integer :: iband,ikgrid,ikpt1,indx
 integer :: ikpt,jkpt,kkpt,ikpt_fine, ik1, ik2, ik3
 integer :: isppol,isym,itim,maxband,minband,nk1,nk2,nk3,nkptfull,uvtk,timrev
 real(dp) :: ene,res,ss,timsign
 logical :: found
 character(len=500) :: msg, format_str
!arrays
 integer,allocatable :: fulltoirred(:)
 real(dp) :: kconv(3),kpt(3),kptgrid(3),kptsym(3)

! *************************************************************************

 ierr=0

!Error if klatt is no simple orthogonal lattice (in red space)
!for generalization to MP grids, need new version of XCrysDen

 if ( kptrlatt(1,2)/=0 .or. kptrlatt(1,3)/=0 .or. kptrlatt(2,1)/=0 .or. &
& kptrlatt(2,3)/=0 .or. kptrlatt(3,1)/=0 .or. kptrlatt(3,2)/=0 ) then
   write(msg,'(3a)')&
&   'kptrlatt should be diagonal, for the FS calculation ',ch10,&
&   'action: use an orthogonal k-grid for the GS calculation '
   MSG_COMMENT(msg)
   ierr=ierr+1
 end if

 if (ANY(ABS(shiftk(:,:))>tol10)) then
   write(msg,'(3a)')&
&   'Origin of the k-grid should be (0,0,0) for the FS calculation ',ch10,&
&   'Action: use a non-shifted k-grid for the GS calculation. Returning '
   MSG_COMMENT(msg)
   ierr=ierr+1
 end if

 if (ierr/=0) RETURN

!Xcrysden uses aperiodical data-grid
 nk1 = kptrlatt(1,1)
 nk2 = kptrlatt(2,2)
 nk3 = kptrlatt(3,3)
 nkptfull=(nk1+1)*(nk2+1)*(nk3+1)

 ABI_ALLOCATE(fulltoirred,(nkptfull))
 timrev=0; if (use_tr) timrev=1

!Xcrysden employs the C-ordering for the Fermi Surface.
 ierr = 0
 ikgrid=0
 do ik1=0,nk1
   do ik2=0,nk2
     do ik3=0,nk3

       ikgrid=ikgrid+1
       kptgrid(1)=DBLE(ik1)/kptrlatt(1,1)
       kptgrid(2)=DBLE(ik2)/kptrlatt(2,2)
       kptgrid(3)=DBLE(ik3)/kptrlatt(3,3)
       call wrap2_pmhalf(kptgrid(1),kpt(1),res)
       call wrap2_pmhalf(kptgrid(2),kpt(2),res)
       call wrap2_pmhalf(kptgrid(3),kpt(3),res)

!      === Find correspondence between the Xcrysden grid and the IBZ ===
!      * If AFM case, use only Ferromagetic symmetries.
       found=.FALSE.
       irred: do ikpt1=1,nkptirred
         do itim=0,timrev
           do isym=1,nsym
             if (use_afm.and.symafm(isym)==-1) CYCLE
             timsign = one-two*itim
             kptsym(:) = timsign*(symrec(:,1,isym)*kptirred(1,ikpt1) + &
&             symrec(:,2,isym)*kptirred(2,ikpt1) + &
&             symrec(:,3,isym)*kptirred(3,ikpt1))
             call wrap2_pmhalf(kptsym(1),kconv(1),res)
             call wrap2_pmhalf(kptsym(2),kconv(2),res)
             call wrap2_pmhalf(kptsym(3),kconv(3),res)
!            * is kconv equivalent to kpt?
             ss= (kpt(1)-kconv(1))**2 + (kpt(2)-kconv(2))**2 + (kpt(3)-kconv(3))**2
             if (ss < tol6) then
               found=.TRUE.
               fulltoirred(ikgrid)=ikpt1
               exit irred
             end if

           end do !itim
         end do !isym
       end do irred

       if (.not.found) then
         write(msg,'(a,3es16.8,2a)')&
&         ' kpt = ',kpt,ch10,' has no symmetric among the irred k-points used in the GS calculation '
         ierr=ierr+1
         MSG_ERROR(msg)
       end if

     end do !ik1
   end do !ik2
 end do !ik3


 if (ierr/=0) then
   ABI_DEALLOCATE(fulltoirred)
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

!=== Dump the results on file ===
 if (open_file(fname,msg,newunit=uvtk,status='unknown',form='formatted') /= 0) then
   ABI_DEALLOCATE(fulltoirred)
   MSG_WARNING(msg)
   ierr=ierr +1; RETURN
 end if

!write header
 write(uvtk,"(a)") '# vtk DataFile Version 2.0'
 write(uvtk,"(a)") 'Eigen values for the Fermi surface'
 write(uvtk,"(a)") 'ASCII'
 write(uvtk,*) ''
 write(uvtk,"(a)") 'DATASET STRUCTURED_GRID'
 write(uvtk,"(a,3i6)") 'DIMENSIONS', nk1+1,nk2+1,nk3+1
 write(uvtk,"(a,i6,a)") 'POINTS',nkptfull,' float'

 do ik3 = 0, nk3
   do ik2 = 0, nk2
     do ik1 = 0, nk1
       write(uvtk,'(3es16.8)') dble(ik1)/nk1*gprimd(1,1)+ &
&       dble(ik2)/nk2*gprimd(1,2)+ &
&       dble(ik3)/nk3*gprimd(1,3), &
&       dble(ik1)/nk1*gprimd(2,1)+ &
&       dble(ik2)/nk2*gprimd(2,2)+ &
&       dble(ik3)/nk3*gprimd(2,3), &
&       dble(ik1)/nk1*gprimd(3,1)+ &
&       dble(ik2)/nk2*gprimd(3,2)+ &
&       dble(ik3)/nk3*gprimd(3,3)
     end do
   end do
 end do

!print out data for all relevant bands and full kpt grid (redundant, yes)
!for each kpt in full zone, find equivalent irred kpt and print eigenval
 write(uvtk,*) ''
 write(uvtk,"(a,i6)") 'POINT_DATA',nkptfull
 indx=0
 do iband=minband,maxband
   do isppol=1,nsppol
     if (minband+indx < 10) then
       format_str="(a14,i1,1X,a)"
     else
       format_str="(a14,i2,1X,a)"
     end if
     write(uvtk,format_str) 'SCALARS eigval', minband+indx, 'float 1'
     write(uvtk,"(a)") 'LOOKUP_TABLE default'
     write(uvtk,*) ' '
     do kkpt = nk3/2+1, nk3+nk3/2+1
       do jkpt = nk2/2+1, nk2+nk2/2+1
         do ikpt = nk1/2+1, nk1+nk1/2+1
           ik1 = ikpt
           ik2 = jkpt
           ik3 = kkpt
           if (ikpt > nk1+1) ik1 = ikpt - nk1
           if (jkpt > nk2+1) ik2 = jkpt - nk2
           if (kkpt > nk3+1) ik3 = kkpt - nk3
!          get the index with zyx order
           ikpt_fine = (ik1-1)*(nk2+1)*(nk3+1) + (ik2-1)*(nk3+1) + ik3
           write(uvtk,'(es16.8)') eigen(iband,fulltoirred(ikpt_fine),isppol)
         end do
       end do
     end do
     indx=indx+1
   end do
 end do

 write(uvtk,*) ''
 indx=0
 do iband=minband,maxband
   do isppol=1,nsppol
     if (minband+indx < 10) then
       format_str="(a10,i1,1X,a)"
     else
       format_str="(a10,i2,1X,a)"
     end if
     write(uvtk,format_str) 'SCALARS ve', minband+indx, 'float'
     write(uvtk,"(a)") 'LOOKUP_TABLE default'
     write(uvtk,*) ' '
     do kkpt = nk3/2+1, nk3+nk3/2+1
       do jkpt = nk2/2+1, nk2+nk2/2+1
         do ikpt = nk1/2+1, nk1+nk1/2+1
           ik1 = ikpt
           ik2 = jkpt
           ik3 = kkpt
           if (ikpt > nk1+1) ik1 = ikpt - nk1
           if (jkpt > nk2+1) ik2 = jkpt - nk2
           if (kkpt > nk3+1) ik3 = kkpt - nk3
!          write(uvtk,'(3i6,3es16.8)') ik1,ik2,ik3,v_surf(iband,ik1,ik2,ik3,1,isppol), &
!          &                                                 v_surf(iband,ik1,ik2,ik3,2,isppol), &
!          &                                                 v_surf(iband,ik1,ik2,ik3,3,isppol)
           write(uvtk,'(es16.8)') sqrt(v_surf(iband,ik1,ik2,ik3,1,isppol)*v_surf(iband,ik1,ik2,ik3,1,isppol)+ &
&           v_surf(iband,ik1,ik2,ik3,2,isppol)*v_surf(iband,ik1,ik2,ik3,2,isppol)+ &
&           v_surf(iband,ik1,ik2,ik3,3,isppol)*v_surf(iband,ik1,ik2,ik3,3,isppol))
         end do
       end do
     end do
     write(uvtk,format_str) 'SCALARS vx', minband+indx, 'float'
     write(uvtk,"(a)") 'LOOKUP_TABLE default'
     write(uvtk,*) ' '
     do kkpt = nk3/2+1, nk3+nk3/2+1
       do jkpt = nk2/2+1, nk2+nk2/2+1
         do ikpt = nk1/2+1, nk1+nk1/2+1
           ik1 = ikpt
           ik2 = jkpt
           ik3 = kkpt
           if (ikpt > nk1+1) ik1 = ikpt - nk1
           if (jkpt > nk2+1) ik2 = jkpt - nk2
           if (kkpt > nk3+1) ik3 = kkpt - nk3
           write(uvtk,'(es16.8)') v_surf(iband,ik1,ik2,ik3,1,isppol)
         end do
       end do
     end do
     write(uvtk,format_str) 'SCALARS vy', minband+indx, 'float'
     write(uvtk,"(a)") 'LOOKUP_TABLE default'
     write(uvtk,*) ' '
     do kkpt = nk3/2+1, nk3+nk3/2+1
       do jkpt = nk2/2+1, nk2+nk2/2+1
         do ikpt = nk1/2+1, nk1+nk1/2+1
           ik1 = ikpt
           ik2 = jkpt
           ik3 = kkpt
           if (ikpt > nk1+1) ik1 = ikpt - nk1
           if (jkpt > nk2+1) ik2 = jkpt - nk2
           if (kkpt > nk3+1) ik3 = kkpt - nk3
           write(uvtk,'(es16.8)') v_surf(iband,ik1,ik2,ik3,2,isppol)
         end do
       end do
     end do
     write(uvtk,format_str) 'SCALARS vz', minband+indx, 'float'
     write(uvtk,"(a)") 'LOOKUP_TABLE default'
     write(uvtk,*) ' '
     do kkpt = nk3/2+1, nk3+nk3/2+1
       do jkpt = nk2/2+1, nk2+nk2/2+1
         do ikpt = nk1/2+1, nk1+nk1/2+1
           ik1 = ikpt
           ik2 = jkpt
           ik3 = kkpt
           if (ikpt > nk1+1) ik1 = ikpt - nk1
           if (jkpt > nk2+1) ik2 = jkpt - nk2
           if (kkpt > nk3+1) ik3 = kkpt - nk3
           write(uvtk,'(es16.8)') v_surf(iband,ik1,ik2,ik3,3,isppol)
         end do
       end do
     end do
     indx=indx+1
   end do
 end do

 close (uvtk)
 ABI_DEALLOCATE(fulltoirred)

end subroutine printvtk
!!***
