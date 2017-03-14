!{\src2tex{textfont=tt}}
!!****f* ABINIT/ddb_diel
!!
!! NAME
!! ddb_diel
!!
!! FUNCTION
!! Get the frequency-dependent dielectric matrix, as well as the
!! oscillator strengths and mode effective charges,
!! and reflectivities (without damping)
!! See the definitions Eq.(53-54) in PRB55, 10355 (1997).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (XG,XW)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! anaddb_dtset= (derived datatype) contains all the input variables
!! matrix (diagonal in the atoms)
!! displ(2,3*natom,3*natom)=
!!  the displacements of atoms in cartesian coordinates.
!!  The first index means either the real or the imaginary part,
!!  The second index runs on the direction and the atoms displaced
!!  The third index runs on the modes.
!! d2cart(2,3,mpert,3,mpert)=dynamical matrix, effective charges, dielectric tensor,... all in cartesian coordinates
!! iout=unit number for outputs
!! lst(3*nph2l)=log. of product of frequencies**2, needed to calculate
!!  the generalized Lyddane-Sachs-Teller relation at zero frequency
!! mpert =maximum number of ipert
!! natom=number of atoms in unit cell
!! nph2l=input variable from anaddb_dtset, needed to dimension lst
!! phfrq(3*natom)=phonon frequencies (square root of the dynamical
!!  matrix eigenvalues, except if these are negative, and in this
!!  case, give minus the square root of the absolute value
!!  of the matrix eigenvalues). Hartree units.
!! comm=MPI communicator.
!!
!! OUTPUT
!! fact_oscstr(2,3,3*natom)=oscillator strengths for the different eigenmodes,
!!  for different direction of the electric field;
!! dielt_rlx(3,3) relaxed ion (zero frequency) dielectric tensor.
!!
!! NOTES
!! 1. The phonon frequencies phfrq should correspond to the
!! wavevector at Gamma, without any non-analyticities.
!! 2. Should clean for no imaginary part ...
!! This routine should be used only by one processor.
!! 3. frdiel(3,3,nfreq)= frequency-dependent dielectric tensor
!! mode effective charges for the different eigenmodes,
!! for different direction of the electric field
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      alignph,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ddb_diel(Crystal,amu,anaddb_dtset,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&
& iout,lst,mpert,natom,nph2l,phfrq,comm)

 use defs_basis
 use m_errors
 use m_xmpi
 use m_profiling_abi
 use m_ddb

 use m_anaddb_dataset, only : anaddb_dataset_type
 use m_crystal,        only : crystal_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_diel'
 use interfaces_14_hidewrite
 use interfaces_77_ddb, except_this_one => ddb_diel
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iout,mpert,natom,nph2l,comm
 type(crystal_t),intent(in) :: Crystal
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset
!arrays
 real(dp),intent(in) :: amu(Crystal%ntypat),d2cart(2,3,mpert,3,mpert),lst(nph2l)
 real(dp),intent(in) :: phfrq(3*natom)
 real(dp),intent(inout) :: displ(2,3*natom,3*natom)
 real(dp),intent(out) :: dielt_rlx(3,3),epsinf(3,3),fact_oscstr(2,3,3*natom)

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 integer :: dieflag,i1,idir1,idir2,ifreq,ii,imode,ipert1,iphl2,nfreq
 integer :: nprocs,my_rank
 real(dp) :: afreq,difffr,eps,lst0,q2,usquare,ucvol
 character(len=500) :: message
 logical :: t_degenerate
!arrays
 real(dp) :: qphon(3),refl(3)
 real(dp),allocatable :: frdiel(:,:,:),modez(:,:,:),oscstr(:,:,:,:)
 character(len=1),allocatable :: metacharacter(:)

! *********************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 dieflag=anaddb_dtset%dieflag
 nfreq=anaddb_dtset%nfreq

 ucvol = Crystal%ucvol

!Check the possibility of asking the frequency-dependent
!dielectric tensor (there should be more than one atom in the unit cell)
 if(natom==1)then
   write(message, '(6a)' )&
&   ' ddb_diel : WARNING -',ch10,&
&   '  When there is only one atom in the unit cell',ch10,&
&   '  cell, the dielectric tensor is frequency-independent.',&
&   '  Consequently, dieflag has been reset to 2 . '
   call wrtout(std_out,message,'COLL')
   call wrtout(iout,message,'COLL')
   dieflag=2
 end if

!frdiel(3,3,nfreq)= frequency-dependent dielectric tensor
!modez(2,3,3*natom)=mode effective charges for the different eigenmodes,
!for different directions of the electric field, following
!the definition Eq.(53) in PRB55, 10355 (1997)
!fact_oscstr(2,3,3*natom)=factors of the oscillator strengths
!for the different eigenmodes,
!for different direction of the electric field
 ABI_ALLOCATE(frdiel,(3,3,nfreq))
 ABI_ALLOCATE(modez,(2,3,3*natom))
!oscstr(2,3,3,3*natom)=oscillator strengths, following
!the definition Eq.(54) in PRB55, 10355 (1997)
 ABI_ALLOCATE(oscstr,(2,3,3,3*natom))

! write(std_out,'(a)')' Enter ddb_diel : displ ='
! do imode=1,3*natom
!   write(std_out,'(a,i4,a,12es16.6)')'imode=',imode,' displ(:,:,imode)',displ(:,:,imode)
! end do

!In case the frequency-dependent dielectric tensor is asked
 if(dieflag==1 .or. dieflag==3 .or. dieflag==4)then

   if (anaddb_dtset%alphon > 0) then
     write(message, '(3a)' )&
&     ' The alphon input variable is non-zero, will mix the degenerate phonon modes',ch10,&
&     ' in order to align the mode effective charges with the cartesian axes.'
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')
     call alignph(amu,displ,d2cart,mpert,natom,Crystal%ntypat,phfrq,Crystal%typat)
   end if

!  Get the factors of the oscillator strength, and the mode effective charge for each mode
   do imode=1,3*natom
     usquare=zero
     do idir1=1,3
       do ipert1=1,natom
         i1=idir1+(ipert1-1)*3
         usquare=usquare+displ(1,i1,imode)*displ(1,i1,imode)+displ(2,i1,imode)*displ(2,i1,imode)
       end do
     end do
     do idir2=1,3
       fact_oscstr(:,idir2,imode)=zero
       modez(:,idir2,imode)=zero
       do idir1=1,3
         do ipert1=1,natom
           i1=idir1+(ipert1-1)*3
           fact_oscstr(:,idir2,imode)=fact_oscstr(:,idir2,imode)+&
&           displ(:,i1,imode)*d2cart(1,idir1,ipert1,idir2,natom+2)
           modez(:,idir2,imode)=modez(:,idir2,imode)+&
&           displ(:,i1,imode)*&
&           d2cart(1,idir1,ipert1,idir2,natom+2)/sqrt(usquare)
         end do
       end do
     end do

!DEBUG
!           write(std_out,'(a,i4,a,12es16.6)')'imode=',imode,&
!&           ' displ(:,:,imode)',displ(:,:,imode)
!           write(std_out,'(a,i4,a,6es16.6)')'imode=',imode,&
!&           ' fact_oscstr(:,:,imode)=',fact_oscstr(:,:,imode)
!ENDDEBUG

   end do

!  Examine the degeneracy of each mode. The portability of the echo of the mode effective
!  charges and oscillator strengths is very hard to guarantee. On the contrary,
!  the scalar reductions of these quantities are OK.
   ABI_ALLOCATE(metacharacter,(3*natom))
   do imode=1,3*natom
!    The degenerate modes are not portable
     t_degenerate=.false.
     if(imode>1)then
       if(phfrq(imode)-phfrq(imode-1)<tol6)t_degenerate=.true.
     end if
     if(imode<3*natom)then
       if(phfrq(imode+1)-phfrq(imode)<tol6)t_degenerate=.true.
     end if
     metacharacter(imode)=';'
     if(t_degenerate)metacharacter(imode)='-'
   end do

   if (my_rank == master) then
     !  Write the mode effective charge for each mode
     write(iout, '(a)' )'  '
     write(iout, '(a)' )' Mode effective charges '
     write(iout, '(a)' )' Mode number.     x               y               z '
     do imode=1,3*natom
       write(iout, '(a,i4,3f16.3)' )metacharacter(imode),imode,(modez(1,idir1,imode),idir1=1,3)
     end do

     !  Write the mode effective charge length for each mode
     write(iout, '(a)' )'  '
     write(iout, '(a)' )' Length of mode effective charge for each phonon mode :'
     do imode=1,3*natom,5
       if (3*natom-imode<5) then
         write(iout, '(1x,5es14.6)') (sqrt(modez(1,1,ii)**2+modez(1,2,ii)**2+modez(1,3,ii)**2),ii=imode,3*natom)
       else
         write(iout, '(1x,5es14.6)') (sqrt(modez(1,1,ii)**2+modez(1,2,ii)**2+modez(1,3,ii)**2),ii=imode,imode+4)
       end if
     end do
   end if ! master

!  Get the oscillator strengths
   do imode=1,3*natom
     do idir1=1,3
       do idir2=1,3
         oscstr(1,idir1,idir2,imode)= &
&         fact_oscstr(1,idir1,imode)*fact_oscstr(1,idir2,imode) +&
&         fact_oscstr(2,idir1,imode)*fact_oscstr(2,idir2,imode)
         if(abs(oscstr(1,idir1,idir2,imode))<tol14)oscstr(1,idir1,idir2,imode)=zero

!DEBUG
!         if(idir1==1 .and. idir2==2)then
!           write(std_out,'(a,i4,a,5es16.6)')'imode=',imode,&
!&           ' oscstr(1,idir1,idir2,imode), fact_oscstr(:,idir1,imode),fact_oscstr(:,idir2,imode)=',&
!&            oscstr(1,idir1,idir2,imode), fact_oscstr(:,idir1,imode),fact_oscstr(:,idir2,imode)
!         endif
!ENDDEBUG

         oscstr(2,idir1,idir2,imode)= &
&         fact_oscstr(1,idir1,imode)*fact_oscstr(2,idir2,imode) -&
&         fact_oscstr(2,idir1,imode)*fact_oscstr(1,idir2,imode)
         if(abs(oscstr(2,idir1,idir2,imode))<tol14)oscstr(2,idir1,idir2,imode)=zero
       end do
     end do
   end do

   if (my_rank == master) then
     !  Write the oscillator strength for each mode
     write(iout, '(a)' )'  '
     write(iout, '(a)' )' Oscillator strengths (in a.u. ; 1 a.u.=253.2638413 m3/s2). Set to zero if abs()<tol14.'
     write(iout, '(a)' )' Mode number.       xx          yy          zz          xy          xz          yz '
     do imode=1,3*natom
       write(iout, '(a,i4,a,6es12.4)' )&
&       metacharacter(imode),imode,'     Real  ',(oscstr(1,idir1,idir1,imode),idir1=1,3),&
&       oscstr(1,1,2,imode), oscstr(1,1,3,imode),oscstr(1,2,3,imode)
       write(iout, '(a,a,6es12.4)' )&
&       metacharacter(imode),'         Imag  ',(oscstr(2,idir1,idir1,imode),idir1=1,3),&
&       oscstr(2,1,2,imode),oscstr(2,1,3,imode),oscstr(2,2,3,imode)
     end do

     !  Write the trace of oscillator strength (real part only) for each mode
     write(iout, '(a)' )'  '
     write(iout, '(a)' )' Trace of oscillator strength, for each phonon mode :'
     do imode=1,3*natom,5
       if (3*natom-imode<5) then
         write(iout, '(1x,5es14.6)') ((oscstr(1,1,1,ii)+oscstr(1,2,2,ii)+oscstr(1,3,3,ii)),ii=imode,3*natom)
       else
         write(iout, '(1x,5es14.6)') ((oscstr(1,1,1,ii)+oscstr(1,2,2,ii)+oscstr(1,3,3,ii)),ii=imode,imode+4)
       end if
     end do
   end if

   ABI_DEALLOCATE(metacharacter)

 end if !  end the condition on frequency-dependent dielectric tensor

!In case the electronic dielectric tensor is needed
 if(dieflag==1.or.dieflag==2.or.dieflag==3 .or. dieflag==4)then

   write(message, '(a,a)' ) ch10,' Electronic dielectric tensor'
   call wrtout(std_out,message,'COLL')
   call wrtout(iout,message,'COLL')

   do idir1=1,3
     do idir2=1,3
       epsinf(idir1,idir2)=d2cart(1,idir1,natom+2,idir2,natom+2)
     end do
     write(message, '(3f16.8)' )(epsinf(idir1,idir2),idir2=1,3)
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')
   end do
   call wrtout(iout, " ",'COLL')
   call wrtout(std_out, " ",'COLL')
 end if

!Only in case the frequency-dependent dielectric tensor is needed
 if(dieflag==1 .or. dieflag==3 .or. dieflag==4) then
!  Check the acousticity of the three lowest modes, assuming
!  that they are ordered correctly
   if (abs(phfrq(1))>abs(phfrq(4)))then
!    This means that there is at least one mode with truly negative frequency
     write(message, '(12a,4es16.8)' )&
&     'The lowest mode appears to be a "true" negative mode,',ch10,&
&     'and not an acoustic mode. This precludes the computation',ch10,&
&     'of the frequency-dependent dielectric tensor.',ch10,&
&     'Action : likely there is no action to be taken, although you,',ch10,&
&     'could try to raise your convergence parameters (ecut and k-points).',ch10,&
&     'For your information, here are the four lowest frequencies :',ch10,&
&     (phfrq(ii),ii=1,4)
     MSG_ERROR(message)
   end if

!  Extract the relaxed ion dielectric tensor
   do idir1=1,3
     do idir2=1,3
       dielt_rlx(idir1,idir2)=epsinf(idir1,idir2)
       do imode=4,3*natom
!        Note that the acoustic modes are not included : their
!        oscillator strength should be exactly zero
!        Also, only the real part of oscstr is taken into account:
!        the possible imaginary parts of degenerate modes
!        will cancel.
         dielt_rlx(idir1,idir2)=dielt_rlx(idir1,idir2)+&
&         oscstr(1,idir1,idir2,imode) / (phfrq(imode)**2)*four_pi/ucvol
!DEBUG
!         if(idir1==1 .and. idir2==2)then
!           write(std_out,'(a,i4,a,3es16.6)')'imode=',imode,' dielt_rlx(idir1,idir2),oscstr(1,idir1,idir2,imode),phfrq(imode)=',&
!&            dielt_rlx(idir1,idir2),oscstr(1,idir1,idir2,imode),phfrq(imode)
!         endif
!ENDDEBUG
       end do
     end do
   end do

   write(message,'(a,a)') ch10,' Relaxed ion dielectric tensor'
   call wrtout(std_out,message,'COLL')
   call wrtout(iout,message,'COLL')

   do idir1=1,3
     write(message,'(3f16.8)')(dielt_rlx(idir1,idir2),idir2=1,3)
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')
   end do
   call wrtout(iout, " ",'COLL')
   call wrtout(std_out, " ",'COLL')
 end if

!Only in case the frequency-dependent dielectric tensor is needed
 if(dieflag==1) then

   difffr=zero
   if(nfreq>1)difffr=(anaddb_dtset%frmax-anaddb_dtset%frmin)/(nfreq-1)

   if (nfreq>10 .and. my_rank == master) then
     write(iout, '(a,a,a,a,a,a,a,a)' )&
&     ' ddb_diel : the number of frequencies is larger',&
&     ' than 10 => I will consider only',ch10,&
&     ' the three principal directions, assume that the tensor',ch10,&
&     ' is diagonalized, and give dielectric constant and ',ch10,&
&     ' reflectivities.'
     write(iout, '(a,a)' )&
&     ' Frequency(Hartree)    Dielectric constant   ',&
&     '             Reflectivity    '
     write(iout, '(a,a)' )&
&     '                     x           y          z',&
&     '          x        y        z'
   end if

!  Loop on frequencies
   do ifreq=1,nfreq
     afreq=anaddb_dtset%frmin+difffr*(ifreq-1)
     do idir1=1,3
       do idir2=1,3
         frdiel(idir1,idir2,ifreq)=epsinf(idir1,idir2)
         do imode=4,3*natom
!          Note that the acoustic modes are not included : their
!          oscillator strength should be exactly zero
!          Also, only the real part of oscstr is taken into account:
!          the possible imaginary parts of degenerate modes
!          will cancel.
           frdiel(idir1,idir2,ifreq)=frdiel(idir1,idir2,ifreq)+&
&           oscstr(1,idir1,idir2,imode) / (phfrq(imode)**2-afreq**2)*four_pi/ucvol
         end do
       end do
     end do

     ! Write all this information (actually, there should be a choice of units for the frequencies ...
     if (nfreq>10) then
       do idir1=1,3
         if(frdiel(idir1,idir1,ifreq)<=zero)then
           refl(idir1)=one
         else
!          See Gervais and Piriou PRB11,3944(1975).
           refl(idir1)=( (sqrt(frdiel(idir1,idir1,ifreq)) -one) /(sqrt(frdiel(idir1,idir1,ifreq)) +one) )**2
         end if
       end do
       if (my_rank == master) then
         write(iout, '(7es12.4)' )afreq,(frdiel(idir1,idir1,ifreq),idir1=1,3),(refl(idir1),idir1=1,3)
       end if

     else
       if (my_rank == master) then
         write(iout, '(a,es12.4,a)' )' Full dielectric tensor at frequency',afreq,' Hartree'
         do idir1=1,3
           write(iout, '(3es16.8)' ) (frdiel(idir1,idir2,ifreq),idir2=1,3)
         end do
         write(iout, '(a)' )' '
       end if
     end if

   end do ! End of the loop on frequencies
 end if ! End the condition on frequency-dependent dielectric tensor

!Calculation of the Lyddane-Sachs-Teller value of the dielectric constant at zero frequency
 if(anaddb_dtset%nph2l/=0 .and.dieflag==1)then

!  Get the log of product of the square of the frequencies without non-analyticities.
   lst0=zero
   do imode=4,3*natom
     lst0=lst0+2*log(phfrq(imode))
   end do

!  Prepare the output
   write(message, '(a,a,a,a)' ) ch10,&
&   ' Generalized Lyddane-Sachs-Teller relation at zero frequency :',ch10,&
&   ' Direction                     Dielectric constant'
   call wrtout(std_out,message,'COLL')
   call wrtout(iout,message,'COLL')

!  Examine every wavevector in the phonon list
   do iphl2=1,anaddb_dtset%nph2l
     qphon(1:3)=anaddb_dtset%qph2l(1:3,iphl2)
     if( abs(qphon(1))>DDB_QTOL .or. abs(qphon(2))>DDB_QTOL .or. abs(qphon(3))>DDB_QTOL)then
       q2=qphon(1)**2+qphon(2)**2+qphon(3)**2
       eps=qphon(1)**2*epsinf(1,1)+qphon(2)**2*epsinf(2,2)+&
&       qphon(3)**2*epsinf(3,3)+ 2* ( qphon(1)*qphon(2)*epsinf(1,2)+&
&       qphon(1)*qphon(3)*epsinf(1,3)+qphon(2)*qphon(3)*epsinf(2,3))
       eps=eps/q2*exp(lst(iphl2)-lst0)
       if (my_rank == master) then
         write(iout, '(3f10.5,es18.8)' )qphon,eps
         write(std_out,'(3f10.5,es18.8)' )qphon,eps
       end if
     end if
   end do
 end if ! End of the condition of nph2l does not vanish

 ABI_DEALLOCATE(frdiel)
 ABI_DEALLOCATE(modez)
 ABI_DEALLOCATE(oscstr)

end subroutine ddb_diel
!!***
