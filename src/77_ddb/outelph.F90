!{\src2tex{textfont=tt}}
!!****f* ABINIT/outelph
!! NAME
!! outelph
!!
!! FUNCTION
!!  Output to stdout and file the data for electron phonon coupling, 
!!  on the q-points which were really calculated by abinit (no interpolation yet)
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2017 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  elph_ds  the elph_type structured variable
!!  enunit   from the anaddb dataset 0 ==> Hartree and cm-1;
!!                                   1 ==> meV and Thz;
!!
!! OUTPUT
!!  only write
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      bfactor,destroy_kptrank,mkkptrank,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine outelph(elph_ds,enunit,fname)

 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_errors
 use m_kptrank

 use m_io_tools,   only : open_file
 use m_nesting,    only : bfactor

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outelph'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: enunit
 character(len=fnlen),intent(in) :: fname
 type(elph_type),intent(in) :: elph_ds

!Local variables-------------------------------
!scalars
 integer :: ibranch,ii,iqfull,iqirr,isppol,jj,nfile,qmax,qnest_max,qnest_min
 integer :: nbranch,nsppol,nqptirred
 real(dp) :: lambda_q_max,lambda_qbranch_max,lambda_tot,nest_max,nest_min
 real(dp) :: omegalog_q,omegalog_qgrid,tc_macmill
 character(len=500) :: msg
 type(kptrank_type) :: kptrank_t
!arrays
 integer :: qbranch_max(2)
 real(dp),allocatable :: lambda_q(:,:),nestfactor(:),qirred(:,:)

! *************************************************************************

 if ( ALL (enunit /= (/0,1,2/)) )  then
   write(msg,'(a,i0)')' enunit should be 0 or 1 or 2 while it is ',enunit
   MSG_BUG(msg)
 end if

 nbranch   = elph_ds%nbranch
 nsppol    = elph_ds%nsppol
 nqptirred = elph_ds%nqptirred 

!==========================================================
!write header
!==========================================================
 if (open_file(fname,msg,newunit=nfile,form="formatted",status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if 

 write(msg,'(2a,80a,4a,80a)')ch10,' ',('=',ii=1,80),ch10,&
& ' Values of the parameters that define the electron-phonon calculation',ch10,&
& ' ',('=',ii=1,80)
 call wrtout(nfile,msg,'COLL')

 write(msg,'(a,i10,a,i10,a,i10)')&
& ' nkpt_phon    = ',elph_ds%k_phon%nkpt,   ' nkpt_phonirred = ',elph_ds%k_phon%nkptirr,&
& ' nqpt      = ',elph_ds%nqpt_full
 call wrtout(nfile,msg,'COLL')

 if (nsppol==1) then
   write(msg,'(2a,f10.7,a,f10.6,a,f10.7)')ch10,&
&   ' Fermi DOS = ',elph_ds%n0(1),       ' Fermi level = ',elph_ds%fermie,&
&   ' mustar    = ',elph_ds%mustar
   call wrtout(nfile,msg,'COLL')
 else if (nsppol==2) then
   write(msg,'(2a,f10.7,f10.7,a,f10.6,a,f10.7)')ch10,&
&   ' Fermi DOS (up/dn) = ',elph_ds%n0(1),elph_ds%n0(2),       ' Fermi level = ',elph_ds%fermie,&
&   ' mustar    = ',elph_ds%mustar
   call wrtout(nfile,msg,'COLL')
 else 
   MSG_BUG("bad value for nsppol")
 end if

 write(msg,'(2a,i10,a,i10,a,i10)')ch10,&
& ' minFSband = ',elph_ds%minFSband,' maxFSband   = ',elph_ds%maxFSband,&
& ' ngkkband  = ',elph_ds%ngkkband
 call wrtout(nfile,msg,'COLL')

 write(msg,'(80a,a)')('=',ii=1,80),ch10
 call wrtout(nfile,msg,'COLL')

!==========================================================
!evaluate lambda and omega_log as a weighted sum over the q grid
!NOTE: in this part of the code atomic units are used
!==========================================================

 ABI_ALLOCATE(lambda_q,(nqptirred,nsppol))
 lambda_q=zero
 lambda_tot=zero ; lambda_q_max=zero
 qmax=0          ; lambda_qbranch_max=zero
 qbranch_max(:)=1; omegalog_qgrid=zero

 do iqirr=1,nqptirred
   omegalog_q=zero

   do isppol=1,nsppol
     do ibranch=1,nbranch
!      find Max lambda(q,n)
       if (elph_ds%qgrid_data(iqirr,ibranch,isppol,3) > lambda_qbranch_max) then
         lambda_qbranch_max=elph_ds%qgrid_data(iqirr,ibranch,isppol,3)
         qbranch_max(1)=iqirr
         qbranch_max(2)=ibranch
       end if
       lambda_q(iqirr,isppol)=lambda_q(iqirr,isppol)+elph_ds%qgrid_data(iqirr,ibranch,isppol,3)
       if (abs(elph_ds%qgrid_data(iqirr,ibranch,isppol,1)) <= tol10) cycle
       omegalog_q=omegalog_q + elph_ds%qgrid_data(iqirr,ibranch,isppol,3)*log(abs(elph_ds%qgrid_data(iqirr,ibranch,isppol,1)))
     end do

     lambda_tot=lambda_tot+elph_ds%wtq(elph_ds%qirredtofull(iqirr))*lambda_q(iqirr,isppol)
     omegalog_qgrid=omegalog_qgrid+elph_ds%wtq(elph_ds%qirredtofull(iqirr))*omegalog_q


!    find Max lambda(q)
     if (lambda_q(iqirr,isppol) > lambda_q_max) then
       lambda_q_max=lambda_q(iqirr,isppol)
       qmax=iqirr
     end if
   end do

 end do !iqirr

 omegalog_qgrid=exp(omegalog_qgrid/lambda_tot)

 write (msg,'(3a,2(a,es16.8))')                                                                              &
& ' Values of Lambda, Omega_log and Tc obtained using the weighted sum over the input Q-grid',ch10,ch10,&
& ' Isotropic Lambda = ',lambda_tot,'  Input mustar     = ',elph_ds%mustar
 call wrtout(nfile,msg,'COLL')

 if (enunit==0) then !use hartree and cm-1
   write (msg,'(2a,es16.8,a,es16.8,a)')ch10,&
&   ' Omega_log        = ',omegalog_qgrid,' (Ha) ',omegalog_qgrid*Ha_cmm1,' (cm-1)'
   call wrtout(nfile,msg,'COLL')
 else if (enunit==1) then !mev Thz
   write (msg,'(2a,es16.8,a,es16.8,a)')ch10,&
&   ' Omega_log        = ',omegalog_qgrid*Ha_eV/1000._dp,' (meV) ',omegalog_qgrid*Ha_THz,' (THz)'
   call wrtout(nfile,msg,'COLL')
 else !hartree,cm-1,mev,Thz,kelvin
   write (msg,'(2a,es16.8,a,es16.8,3a,es16.8,a,es16.8,3a,es16.8,a)')ch10,                              &
&   ' Omega_log        = ',omegalog_qgrid,' (Ha)  ',omegalog_qgrid*Ha_cmm1,' (cm-1)',ch10,             &
&   '                  = ',omegalog_qgrid*Ha_eV/1000._dp,' (meV) ',omegalog_qgrid*Ha_THz,' (THz)',ch10,&
&   '                  = ',omegalog_qgrid*Ha_K,' (K) '
   call wrtout(nfile,msg,'COLL')
 end if

 tc_macmill = omegalog_qgrid/1.2_dp&
& *exp((-1.04_dp*(one+lambda_tot)) / (lambda_tot-elph_ds%mustar*(one+0.62_dp*lambda_tot)))

 if (enunit==0) then !use hartree and cm-1
   write (msg,'(2a,es16.8,a,es16.8,2a)')ch10,&
&   ' MacMillan Tc     = ',tc_macmill,' (Ha) ',tc_macmill*Ha_cmm1,' (cm-1) ',ch10
   call wrtout(nfile,msg,'COLL')
 else if (enunit==1) then !use mev and Thz
   write (msg,'(2a,es16.8,a,es16.8,2a)')ch10,&
&   ' MacMillan Tc     = ',tc_macmill*Ha_eV/1000._dp,' (meV) ',tc_macmill*Ha_THz,' (THz) ',ch10
   call wrtout(nfile,msg,'COLL')
 else !use hartree,cm-1,mev,Thz,kelvin
   write (msg,'(2a,es16.8,a,es16.8,3a,es16.8,a,es16.8,3a,es16.8,2a)')ch10,                 &
&   ' MacMillan Tc     = ',tc_macmill,' (Ha)  ',tc_macmill*Ha_cmm1,' (cm-1) ',ch10,            &
&   '                  = ',tc_macmill*Ha_eV/1000._dp,' (meV) ',tc_macmill*Ha_THz,' (THz) ',ch10,&
&   '                  = ',tc_macmill*Ha_K,' (K) ',ch10
   call wrtout(nfile,msg,'COLL')
 end if

!==========================================================
!output lambda(q) values for each q point in the irred grid
!==========================================================

 write(msg,'(2a)')' Irreducible q-points and corresponding Lambda(q)',ch10
 call wrtout(nfile,msg,'COLL')

 do isppol=1,nsppol
   write(msg,'(a,i3,2a)')'  === isppol ', isppol,' === ',ch10
   call wrtout(nfile,msg,'COLL')
!  
   do iqirr=1,nqptirred
     iqfull=elph_ds%qirredtofull(iqirr)
     write(msg,'(i5,a,3(es16.8,1x),a,es16.8,a)')&
&     iqfull,') ',elph_ds%qpt_full(:,iqfull),'(',lambda_q(iqirr,isppol),'  )'
     call wrtout(nfile,msg,'COLL')
   end do
!  
 end do

!use same indexing as that used for the full q-grid
 qmax=elph_ds%qirredtofull(qmax)
 qbranch_max(1)=elph_ds%qirredtofull(qbranch_max(1))

 write (msg,'(2a,es16.8,a,i6,3a,es16.8,a,i6,a,i4)')ch10,            &
& ' Max lambda(q)      = ',lambda_q_max,      ' at qpt ',qmax,')',ch10, &
& ' Max lambda(q,n)    = ',lambda_qbranch_max,' at qpt ',qbranch_max(1),&
& ') and Mode number ',qbranch_max(2)
 call wrtout(nfile,msg,'COLL')

!==========================================================
!evaluation of the nesting-factor over the irreducible q grid.
!==========================================================

!fill irreducile q-grid
 ABI_ALLOCATE(qirred,(3,nqptirred))
 qirred(:,:)=zero

 do iqirr=1,nqptirred
   qirred(:,iqirr)=elph_ds%qpt_full(:,elph_ds%qirredtofull(iqirr))
 end do

 call mkkptrank (elph_ds%k_phon%kpt,elph_ds%k_phon%nkpt,kptrank_t)

 ABI_ALLOCATE(nestfactor,(nqptirred))

!NOTE: weights are not normalised, the normalisation factor in reintroduced in bfactor
 call bfactor(elph_ds%k_phon%nkpt,elph_ds%k_phon%kpt,nqptirred,qirred,kptrank_t,&
& elph_ds%k_phon%nkpt,elph_ds%k_phon%wtk,elph_ds%nFSband,nestfactor)

 ABI_DEALLOCATE(qirred)
 call destroy_kptrank (kptrank_t)


!find Max and min of the nesting factor
!NOTE maxloc and minloc are arrays so they cannot be used in the formatted output
!anyway the size of nestfactor is not so huge!!!
 nest_max=maxval(nestfactor); nest_min=minval(nestfactor)

 qnest_max=0
 do iqirr=1,nqptirred
   if (nestfactor(iqirr)==nest_max) then
     qnest_max=iqirr
     exit
   end if
 end do

 qnest_min=0
 do iqirr=1,nqptirred
   if (nestfactor(iqirr)==nest_min) then
     qnest_min=iqirr
     exit
   end if
 end do

 write (std_out,*) maxloc(nestfactor),minloc(nestfactor)
 write(msg,'(a,(a,es16.8,a,i6,a),a,(a,es16.8,a,i6,a))')ch10,  &
& ' Max nesting factor = ',nest_max,' at qpt ',qnest_max,') ',ch10,&
& ' min nesting factor = ',nest_min,' at qpt ',qnest_min,') '
 call wrtout(nfile,msg,'COLL')

!==========================================================
!Write ph-linewidths and lambda(q,n) obtained before the
!Fourier interpolation
!==========================================================

 write (msg,'(2a)')ch10,&
& ' Phonon frequencies, linewidths and e-ph coefficients for each irreducible q point '
 call wrtout(nfile,msg,'COLL')

 do isppol=1,nsppol
   write (msg,'(a,i3,a)') '========= quantities for isppol = ', isppol, ' ================='
   call wrtout(nfile,msg,'COLL')
   do iqirr=1,nqptirred
!    same numbering as that used for irred q points
     iqfull=elph_ds%qirredtofull(iqirr)
!    write(std_out,*) 'iqfull = ', iqfull
     write(msg,'(64a,i6,a,3(es16.8),3a,es16.8,a,es16.8,2a,es16.8,a,f8.3,65a)')ch10,&
&     ' ',('=',jj=1,60),ch10,&
&     ' qpt ',iqfull,') ',elph_ds%qpt_full(:,iqfull),ch10,ch10,&
&     ' Weight    = ',elph_ds%wtq(iqfull),'    Lambda(q,isppol) = ',lambda_q(iqirr,isppol),ch10,&
&     ' Nest fact = ',nestfactor(iqirr),'    (',100*nestfactor(iqirr)/nest_max,' % of max_value )',ch10,&
&     ' ',('=',jj=1,60),ch10,' Mode number    Frequency       Linewidth        Lambda(q,n)'
     call wrtout(nfile,msg,'COLL')

!    use units according to enunit
     if (enunit==0 .or. enunit==2) then !hartree and cm-1
       write(msg,'(63a)')' ',('-',jj=1,60),ch10,&
       '                  (Ha)             (Ha)'
       call wrtout(nfile,msg,'COLL')
       do ibranch=1,nbranch
!        branch index, frequency, linewidth, lamda(q,n) (hartree units)
         write(msg,'(i6,5x,3(es16.8,1x))' )ibranch,(elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,3)
         call wrtout(nfile,msg,'COLL')
       end do
       write(msg,'(63a)')' ',('-',jj=1,60),ch10,&
&       '                 (cm-1)           (cm-1)'
       call wrtout(nfile,msg,'COLL')
       do ibranch=1,nbranch
!        branch index, frequency, linewidth (in cm-1)
         write(msg,'(i6,5x,2(es16.8,1x))' )ibranch,(Ha_cmm1*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2)
         call wrtout(nfile,msg,'COLL')
       end do
     end if !hartree and cm-1

     if (enunit==2 .or. enunit==1) then !write also meV Thz and Kelvin
       write(msg,'(63a)')' ',('-',jj=1,60),ch10,&
&       '                 (meV)             (meV)'
       call wrtout(nfile,msg,'COLL')
       if (enunit == 1 ) then !write also lambda values
         do ibranch=1,nbranch
!          branch index, frequency, linewidth, lamda(q,n) (mev units)
           write(msg,'(i6,5x,3(es16.8,1x))' )ibranch,((Ha_eV/1000._dp)*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2),&
&           elph_ds%qgrid_data(iqirr,ibranch,isppol,3)
           call wrtout(nfile,msg,'COLL')
         end do
       else !do not write lambda values
         do ibranch=1,nbranch
!          branch index, frequency, linewidth (in meV)
           write(msg,'(i6,5x,2(es16.8,1x))' )ibranch,((Ha_eV/1000._dp)*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2)
           call wrtout(nfile,msg,'COLL')
         end do
       end if

       write(msg,'(63a)')' ',('-',jj=1,60),ch10,&
&       '                 (Thz)             (Thz)'
       call wrtout(nfile,msg,'COLL')
       do ibranch=1,nbranch
!        branch index, frequency, linewidth (in Thz)
         write(msg,'(i6,5x,2(es16.8,1x))' )ibranch,(Ha_THz*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2)
         call wrtout(nfile,msg,'COLL')
       end do

       if (enunit == 2 ) then !kelvin
         write(msg,'(63a)')' ',('-',jj=1,60),ch10,&
&         '                  (K)               (K)'
         call wrtout(nfile,msg,'COLL')
         do ibranch=1,nbranch
!          branch index, frequency, linewidth (in Kelvin)
           write(msg,'(i6,5x,2(es16.8,1x))' )ibranch,(Ha_K*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2)
           call wrtout(nfile,msg,'COLL')
         end do
       end if !kelvin

     end if  !end write also meV Thz and Kelvin

     write(msg,'(62a)')' ',('=',jj=1,60),ch10
     call wrtout(nfile,msg,'COLL')

   end do !nqptirred
 end do !nsppol

 ABI_DEALLOCATE(nestfactor)
 ABI_DEALLOCATE(lambda_q)

 close (nfile)

end subroutine outelph
!!***
