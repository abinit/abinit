!{\src2tex{textfont=tt}}
!!****f* ABINIT/prteigrs
!!
!! NAME
!! prteigrs
!!
!! FUNCTION
!! Print out eigenvalues band by band and k point by k point.
!! If option=1, do it in a standard way, for self-consistent calculations.
!! If option=2, print out residuals and eigenvalues, in a format
!! adapted for nonself-consistent calculations, within the loops.
!! If option=3, print out eigenvalues, in a format
!! adapted for nonself-consistent calculations, at the end of the job.
!! If option=4, print out derivatives of eigenvalues (same format as option==3, except header that is printed)
!! If option=5, print out Fan contribution to zero-point motion correction to eigenvalues (averaged)
!!                  (same format as option==3, except header that is printed)
!! If option=6, print out DDW contribution to zero-point motion correction to eigenvalues (averaged)
!!                  (same format as option==3, except header that is printed)
!! If option=7, print out Fan+DDW contribution to zero-point motion correction to eigenvalues (averaged)
!!                  (same format as option==3, except header that is printed)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  eigen(mband*nkpt*nsppol)=eigenvalues (hartree)
!!   or, if option==4, diagonal of derivative of eigenvalues
!!   or, if option==5...7, zero-point motion correction to eigenvalues (averaged)
!!  enunit=choice parameter: 0=>output in hartree; 1=>output in eV;
!!   2=> output in both hartree and eV
!!  fermie=fermi energy (Hartree)
!!  fname_eig=filename for printing of the eigenenergies
!!  iout=unit number for formatted output file
!!  iscf=option for self-consistency
!!  kptns(3,nkpt)=k points in reduced coordinates
!!  kptopt=option for the generation of k points
!!  mband=maximum number of bands
!!  nband(nkpt)=number of bands at each k point
!!  nkpt=number of k points
!!  nnsclo_now=number of non-self-consistent loops for the current vtrial
!!    (often 1 for SCF calculation, =nstep for non-SCF calculations)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(maxval(nband(:))*nkpt*nsppol)=occupancies for each band and k point
!!  occopt=option for occupancies
!!  option= (see above)
!!  prteig=control print eigenenergies
!!  prtvol=control print volume and debugging
!!  resid(mband*nkpt*nsppol)=residuals (hartree**2)
!!  tolwfr=tolerance on band residual of wf, hartrees**2 (needed when option=2)
!!  vxcavg=average of vxc potential
!!  wtk(nkpt)=k-point weights
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      clnup1,dfpt_looppert,respfn,scprqt,vtorho
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prteigrs(eigen,enunit,fermie,fname_eig,iout,iscf,kptns,kptopt,mband,nband,&
&  nkpt,nnsclo_now,nsppol,occ,occopt,option,prteig,prtvol,resid,tolwfr,vxcavg,wtk)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_io_tools,  only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prteigrs'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: enunit,iout,iscf,kptopt,mband,nkpt,nnsclo_now,nsppol
 integer,intent(in) :: occopt,option,prteig,prtvol
 real(dp),intent(in) :: fermie,tolwfr,vxcavg
 character(len=*),intent(in) :: fname_eig
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),kptns(3,nkpt)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),resid(mband*nkpt*nsppol)
 real(dp),intent(in) :: wtk(nkpt)

!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt_max=50
 integer :: band_index,iband,ienunit,ii,ikpt,isppol,nband_index,nband_k,nkpt_eff,tmagnet,tmetal,temp_unit
 real(dp) :: convrt,magnet,residk,rhodn,rhoup
 character(len=2) :: ibnd_fmt,ikpt_fmt
 character(len=7) :: strunit1,strunit2
 character(len=39) :: kind_of_output
 character(len=500) :: message

! *************************************************************************

 if (nsppol<1.or.nsppol>2) then
   write(message, '(a,i0)' )' nsppol must be 1 or 2. Argument was ',nsppol
   MSG_BUG(message)
 end if

 if (enunit<0.or.enunit>2) then
   write(message, '(a,i0)' )' enunit must be 0, 1 or 2. Argument was ',enunit
   MSG_BUG(message)
 end if

 if (prteig > 0) then
   write(message, '(a,a)' ) ' prteigrs : about to open file ',TRIM(fname_eig)
   call wrtout(iout,message,'COLL')
   if (open_file(fname_eig, message, newunit=temp_unit, status='unknown', form='formatted') /= 0) then
     MSG_ERROR(message)
   end if
   rewind(temp_unit) ! always rewind disk file and print latest eigenvalues
 end if

 kind_of_output=              ' Eigenvalues                          '
 if(option==4) kind_of_output=' Expectation of eigenvalue derivatives'
 if(option==5) kind_of_output=' Fan corrections to eigenvalues at T=0'
 if(option==6) kind_of_output=' DDW corrections to eigenvalues at T=0'
 if(option==7) kind_of_output=' Fan+DDW corrs   to eigenvalues at T=0'

 nkpt_eff=nkpt

!DEBUG
!write(message,'(a,5i5)')' prtvol,iscf,kptopt,nkpt_eff,nkpt_max ',prtvol,iscf,kptopt,nkpt_eff,nkpt_max
!call wrtout(iout,message,'COLL')
!ENDDEBUG

 if( (prtvol==0.or.prtvol==1) .and. (iscf/=-2 .or. kptopt>0) .and. nkpt_eff>nkpt_max)nkpt_eff=nkpt_max
 if( (prtvol==0.or.prtvol==1) .and. (iscf/=-2 .or. kptopt>0) .and. nkpt_eff>1 .and. iout==ab_out)nkpt_eff=1

 if(option==1 .or. (option>=3 .and. option<=7))then

   do ienunit=0,1

     if (enunit==1 .and. ienunit==0)cycle
     if (enunit==0 .and. ienunit==1)cycle
!  Print eigenvalues in hartree for enunit=0 or 2
!  The definition of two different strings is quite ridiculous. Historical reasons ...
     
     if (ienunit==0)then
       convrt=one
       strunit1='hartree'
       strunit2='hartree'
     end if
     if (ienunit==1)then
       convrt=Ha_eV
       strunit1='   eV  '
       strunit2='eV     '
     end if

     band_index=0

     if(ienunit==0)then  ! XG20140730 I do not know why this is only done when ienunit==0
       tmetal=0
       if(option==1 .and. occopt>=3 .and. occopt<=8)tmetal=1
       tmagnet=0
       if(tmetal==1 .and. nsppol==2)then
         tmagnet=1
         rhoup = 0._dp
         rhodn = 0._dp
         nband_index = 1
         do isppol=1,nsppol
           do ikpt=1,nkpt
             nband_k=nband(ikpt+(isppol-1)*nkpt)
             do iband=1,nband_k
               if(isppol==1) rhoup = rhoup + wtk(ikpt)*occ(nband_index)
               if(isppol==2) rhodn = rhodn + wtk(ikpt)*occ(nband_index)
               nband_index = nband_index + 1
             end do
           end do
         end do
         magnet = abs(rhoup - rhodn)
       end if
     end if
     
     if(iscf>=0 .and. (ienunit==0 .or. option==1))then
       write(message, '(3a,f10.5,3a,f10.5)' ) &
&       ' Fermi (or HOMO) energy (',trim(strunit2),') =',convrt*fermie,'   Average Vxc (',trim(strunit2),')=',convrt*vxcavg
       call wrtout(iout,message,'COLL')
       if (prteig > 0) then
         call wrtout(temp_unit,message,'COLL')
       end if
     end if

     
!    if( (iscf>=0 .or. iscf==-3) .and. ienunit==0)then     ! This is the most correct
     if(iscf>=0 .and. ienunit==0)then ! For historical reasons
       if(tmagnet==1)then
         write(message, '(a,es16.8,a,a,es16.8,a,es16.8)' )&
&         ' Magnetization (Bohr magneton)=',magnet,ch10,&
&         ' Total spin up =',rhoup,'   Total spin down =',rhodn
         call wrtout(iout,message,'COLL')
         if (prteig > 0) then
           call wrtout(temp_unit,message,'COLL')
         end if
       end if
     end if

!    Loop over spins (suppress spin data if nsppol not 2)
     do isppol=1,nsppol

       ikpt_fmt="i4" ; if(nkpt>=10000)ikpt_fmt="i6" ; if(nkpt>=1000000)ikpt_fmt="i9"
       if (nsppol==2.and.isppol==1) then
         write(message, '(4a,'//ikpt_fmt//',2x,a)' ) &
&         trim(kind_of_output),' (',strunit1,') for nkpt=',nkpt,'k points, SPIN UP:'
       else if (nsppol==2.and.isppol==2) then
         write(message, '(4a,'//ikpt_fmt//',2x,a)' ) &
&         trim(kind_of_output),' (',strunit1,') for nkpt=',nkpt,'k points, SPIN DOWN:'
       else
         write(message, '(4a,'//ikpt_fmt//',2x,a)' ) &
&         trim(kind_of_output),' (',strunit1,') for nkpt=',nkpt,'k points:'
       end if
       call wrtout(iout,message,'COLL')
       if (prteig > 0) then
         call wrtout(temp_unit,message,'COLL')
       end if

       if(ienunit==0)then
         if(option>=4 .and. option<=7)then
           message = '  (in case of degenerate eigenvalues, averaged derivative)'
           call wrtout(iout,message,'COLL')
           if (prteig > 0) then
             call wrtout(temp_unit,message,'COLL')
           end if
         end if
       end if

       do ikpt=1,nkpt
         nband_k=nband(ikpt+(isppol-1)*nkpt)
         ikpt_fmt="i4" ; if(nkpt>=10000)ikpt_fmt="i6" ; if(nkpt>=1000000)ikpt_fmt="i9"
         ibnd_fmt="i3" ; if(nband_k>=1000)ibnd_fmt="i6" ; if(nband_k>=1000000)ibnd_fmt="i9"
         if(ikpt<=nkpt_eff)then
           write(message, '(a,'//ikpt_fmt//',a,'//ibnd_fmt//',a,f9.5,a,3f8.4,a)' ) &
&           ' kpt#',ikpt,', nband=',nband_k,', wtk=',wtk(ikpt)+tol10,', kpt=',&
&           kptns(1:3,ikpt)+tol10,' (reduced coord)'
           call wrtout(iout,message,'COLL')
           if (prteig > 0) then
             call wrtout(temp_unit,message,'COLL')
           end if
           do ii=0,(nband_k-1)/8
!            write(message, '(8f15.10)' ) (convrt*eigen(iband+band_index),&
             write(message, '(8(f10.5,1x))' ) (convrt*eigen(iband+band_index),&
&             iband=1+ii*8,min(nband_k,8+ii*8))
             call wrtout(iout,message,'COLL')
             if (prteig > 0) then
               call wrtout(temp_unit,message,'COLL')
             end if
           end do
           if(ienunit==0 .and. option==1 .and. occopt>=3 .and. occopt<=8)then
             write(message, '(5x,a,'//ikpt_fmt//')' )  ' occupation numbers for kpt#',ikpt
             call wrtout(iout,message,'COLL')
             do ii=0,(nband_k-1)/8
               write(message, '(8(f10.5,1x))' ) (occ(iband+band_index),&
&               iband=1+ii*8,min(nband_k,8+ii*8))
               call wrtout(iout,message,'COLL')
             end do
           end if

         else
           if(ikpt==nkpt_eff+1)then
             write(message, '(a,a)' ) &
&             ' prteigrs : prtvol=0 or 1, do not print more k-points.',ch10
             call wrtout(iout,message,'COLL')
           end if
           if (prteig > 0) then
             write(message, '(a,'//ikpt_fmt//',a,'//ibnd_fmt//',a,f9.5,a,3f8.4,a)' ) &
&             ' kpt#',ikpt,', nband=',nband_k,', wtk=',wtk(ikpt)+tol10,', kpt=',&
&             kptns(1:3,ikpt)+tol10,' (reduced coord)'
             call wrtout(temp_unit,message,'COLL')
             do ii=0,(nband_k-1)/8
               write(message, '(8(f10.5,1x))' ) (convrt*eigen(iband+band_index),&
&               iband=1+ii*8,min(nband_k,8+ii*8))
               call wrtout(temp_unit,message,'COLL')
             end do
           end if
         end if
         band_index=band_index+nband_k
       end do ! do ikpt=1,nkpt
     end do ! do isppol=1,nsppol

   end do ! End loop over Hartree or eV

 else if(option==2)then

   band_index=0
   do isppol=1,nsppol

     if(nsppol==2)then
       if(isppol==1)write(message, '(2a)' ) ch10,' SPIN UP channel '
       if(isppol==2)write(message, '(2a)' ) ch10,' SPIN DOWN channel '
       call wrtout(iout,message,'COLL')
       if(prteig>0) then
         call wrtout(temp_unit,message,'COLL')
       end if
     end if

     do ikpt=1,nkpt
       nband_k=nband(ikpt+(isppol-1)*nkpt)
       ikpt_fmt="i5" ; if(nkpt>=10000)ikpt_fmt="i7" ; if(nkpt>=1000000)ikpt_fmt="i9"

       if(ikpt<=nkpt_eff)then
         write(message, '(1x,a,'//ikpt_fmt//',a,f9.5,2f9.5,a)' ) &
&         'Non-SCF case, kpt',ikpt,' (',(kptns(ii,ikpt),ii=1,3),'), residuals and eigenvalues='
         call wrtout(iout,message,'COLL')
         if (prteig > 0) then
           write(message, '(1x,a,'//ikpt_fmt//',a,f9.5,2f9.5,a)' ) &
&           'Non-SCF case, kpt',ikpt,' eig(',(kptns(ii,ikpt),ii=1,3),') '
           call wrtout(temp_unit,message,'COLL')
         end if
         do ii=0,(nband_k-1)/8
           write(message, '(1p,8e10.2)' )(resid(iband+band_index),iband=1+8*ii,min(8+8*ii,nband_k))
           call wrtout(iout,message,'COLL')
         end do
         do ii=0,(nband_k-1)/6
           write(message, '(1p,6e12.4)' )(eigen(iband+band_index),iband=1+6*ii,min(6+6*ii,nband_k))
           call wrtout(iout,message,'COLL')
           if (prteig > 0) then
             call wrtout(temp_unit,message,'COLL')
           end if
         end do
       else
         if(ikpt==nkpt_eff+1)then
           write(message, '(a,a)' )' prteigrs : prtvol=0 or 1, do not print more k-points.',ch10
           call wrtout(iout,message,'COLL')
         end if
         if (prteig > 0) then
           write(message, '(1x,a,i5,a,f9.5,2f9.5,a)' ) &
&           'Non-SCF kpt',ikpt,' eig(',(kptns(ii,ikpt),ii=1,3),') '
           call wrtout(temp_unit,message,'COLL')
           do ii=0,(nband_k-1)/6
             write(message, '(1p,6e12.4)' )(eigen(iband+band_index),iband=1+6*ii,min(6+6*ii,nband_k))
             call wrtout(temp_unit,message,'COLL')
           end do
         end if
       end if

       residk=maxval(resid(band_index+1:band_index+nband_k))
       if (residk>tolwfr) then
         write(message, '(1x,a,2i5,a,1p,e13.5)' ) &
&         ' prteigrs : nnsclo,ikpt=',nnsclo_now,ikpt,' max resid (incl. the buffer)=',residk
         call wrtout(iout,message,'COLL')
       end if

       band_index=band_index+nband_k
     end do
   end do
   call wrtout(iout," ",'COLL')

 else
   write(message, '(a,i0,a)' )' option = ',option,', is not an allowed value.'
   MSG_BUG(message)
 end if

 if (prteig > 0) close (temp_unit)

end subroutine prteigrs
!!***
