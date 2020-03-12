!!****m* ABINIT/m_ddb_diel
!! NAME
!! m_ddb_diel
!!
!! FUNCTION
!! This module provides routines for the calculation of the dielectric constant (anaddb)
!! 
!! COPYRIGHT
!!  Copyright (C) 1999-2020 ABINIT group (XG,XW, MVeithen, EB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_ddb_diel

 use defs_basis
 use m_errors
 use m_xmpi
 use m_abicore
 use m_ddb
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_anaddb_dataset, only : anaddb_dataset_type
 use m_crystal,        only : crystal_t

 implicit none

 private
!!***

 public :: ddb_diel
! public :: ddb_diel_elec  ! electronic dielectric constant calculation
!!***

contains
!!***

!!****f* ABINIT/ddb_diel
!!
!! NAME
!! ddb_diel
!!
!! FUNCTION
!! Get the frequency-dependent dielectric matrix, as well as the
!! oscillator strengths and mode effective charges,
!! and reflectivities (without damping)
!! See the definitions Eq.(53-54) in PRB55, 10355 (1997) [[cite:Gonze1997a]].
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
!! ncid=the id of the open NetCDF file. Set to nctk_noid if netcdf output is not wanted.
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
!!      alignph,ddb_diel_elec,ddb_oscstr,wrtout
!!
!! SOURCE

subroutine ddb_diel(Crystal,amu,anaddb_dtset,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&
& iout,lst,mpert,natom,nph2l,phfrq,comm,ncid)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iout,mpert,natom,nph2l,comm,ncid
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
 integer :: nprocs,my_rank,ncerr
 real(dp) :: afreq,difffr,eps,lst0,q2,usquare,ucvol
 character(len=500) :: message
 logical :: t_degenerate
!arrays
 real(dp) :: qphon(3),refl(3)
 real(dp),allocatable :: frdiel(:,:,:),modez(:,:,:),oscstr(:,:,:,:),dielt_modedecompo(:,:,:)
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
     ' ddb_diel : WARNING -',ch10,&
     '  When there is only one atom in the unit cell',ch10,&
     '  cell, the dielectric tensor is frequency-independent.',&
     '  Consequently, dieflag has been reset to 2 . '
   call wrtout(std_out,message,'COLL')
   call wrtout(iout,message,'COLL')
   dieflag=2
 end if

!frdiel(3,3,nfreq)= frequency-dependent dielectric tensor
!fact_oscstr(2,3,3*natom)=factors of the oscillator strengths
!for the different eigenmodes,
!for different direction of the electric field
 ABI_ALLOCATE(frdiel,(3,3,nfreq))
 ABI_ALLOCATE(oscstr,(2,3,3,3*natom))
 ABI_ALLOCATE(modez,(2,3,3*natom))
 ABI_ALLOCATE(dielt_modedecompo,(3,3,3*natom))

! write(std_out,'(a)')' Enter ddb_diel : displ ='
! do imode=1,3*natom
!   write(std_out,'(a,i4,a,12es16.6)')'imode=',imode,' displ(:,:,imode)',displ(:,:,imode)
! end do

!In case the ionic contribution to the dielectric tensor is asked
 if(dieflag==1 .or. dieflag==3 .or. dieflag==4 .or. dieflag==5)then
! if(dieflag==1 .or. dieflag==3 .or. dieflag==4)then

! Check if the alignement of phonon modes eigenvector is requested from the input flag alphon;
! useful in case of degenerate modes 
   if (anaddb_dtset%alphon > 0) then
     write(message, '(3a)' )&
      ' The alphon input variable is non-zero, will mix the degenerate phonon modes',ch10,&
      ' in order to align the mode effective charges with the cartesian axes.'
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')
     call alignph(amu,displ,d2cart,mpert,natom,Crystal%ntypat,phfrq,Crystal%typat)
   end if ! alignment of the phonon eigenvectors

!  Compute the mode effective charge and oscillator strength
   call ddb_oscstr(displ,d2cart,fact_oscstr,oscstr,modez,iout,mpert,natom,phfrq,ncid)
 end if  ! ionic contrib to dielectric tensor

!In case the electronic contribution to the dielectric tensor is needed
 if(dieflag==1.or.dieflag==2.or.dieflag==3 .or. dieflag==4 .or. dieflag==5)then
! if(dieflag==1.or.dieflag==2.or.dieflag==3 .or. dieflag==4)then
   call ddb_diel_elec(d2cart,natom,mpert,iout,epsinf)
 end if

!Only in case the frequency-dependent dielectric tensor is needed
! EB: This not only for freq-dep epsilon but for ionic contrib to epsilon!
 if(dieflag==1 .or. dieflag==3 .or. dieflag==4 .or. dieflag==5) then
! if(dieflag==1 .or. dieflag==3 .or. dieflag==4) then
!  Check the acousticity of the three lowest modes, assuming
!  that they are ordered correctly
   if (abs(phfrq(1))>abs(phfrq(4)))then
!    This means that there is at least one mode with truly negative frequency
     write(message, '(12a,4es16.8)' )&
      'The lowest mode appears to be a "true" negative mode,',ch10,&
      'and not an acoustic mode. This precludes the computation',ch10,&
      'of the frequency-dependent dielectric tensor.',ch10,&
      'Action : likely there is no action to be taken, although you,',ch10,&
      'could try to raise your convergence parameters (ecut and k-points).',ch10,&
      'For your information, here are the four lowest frequencies :',ch10,&
       (phfrq(ii),ii=1,4)
     MSG_ERROR(message)
   end if

!  Compute the relaxed ion dielectric tensor
   do idir1=1,3
     do idir2=1,3
!      The electronic contribution to epsilon is added 
       dielt_rlx(idir1,idir2)=epsinf(idir1,idir2)
!      calculation of the phonon contribution (ionic) to epsilon
       do imode=4,3*natom
!        Note that the acoustic modes are not included : their
!        oscillator strength should be exactly zero
!        Also, only the real part of oscstr is taken into account:
!        the possible imaginary parts of degenerate modes
!        will cancel.
         dielt_rlx(idir1,idir2)=dielt_rlx(idir1,idir2)+&
&         oscstr(1,idir1,idir2,imode) / (phfrq(imode)**2)*four_pi/ucvol
!        Mode decomposition of epsilon
         if (dieflag==5)then
           dielt_modedecompo(idir1,idir2,imode)=oscstr(1,idir1,idir2,imode) / (phfrq(imode)**2)*four_pi/ucvol
         endif ! mode decompo of epsilon 
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
   
!  Mode decompo of epsilon
   if (dieflag==5) then
     write(message,'(a,a,a,a)') ch10,' Mode by mode decomposition of the ionic dielectric tensor',&
                                ch10,' (the electronic contribution is not included)'
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')
     do imode=4,3*natom
       write(message,'(a,a,i4)') ch10,' Mode number ',imode
       call wrtout(std_out,message,'COLL')
       call wrtout(iout,message,'COLL')
       do idir1=1,3
        ! write(message,'(a,a,i4)') ch10,' Mode number 2',imode
         write(message,'(3f16.8)')(dielt_modedecompo(idir1,idir2,imode),idir2=1,3)
         call wrtout(std_out,message,'COLL')
         call wrtout(iout,message,'COLL')
       end do 
     end do
   endif ! dieflag=5 mode decompo of epsilon

   ! write the relaxed ion dielectric tensor to the netcdf
#ifdef HAVE_NETCDF
   if (ncid /= nctk_noid) then
     ncerr = nctk_def_arrays(ncid, [nctkarr_t("emacro_cart_rlx", "dp", &
     "number_of_cartesian_directions, number_of_cartesian_directions")],defmode=.True.)
     NCF_CHECK(ncerr)
     NCF_CHECK(nctk_set_datamode(ncid))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "emacro_cart_rlx"), dielt_rlx))
   end if
#endif

 end if

!Only in case the frequency-dependent dielectric tensor is needed
 if(dieflag==1) then

   difffr=zero
   if(nfreq>1)difffr=(anaddb_dtset%frmax-anaddb_dtset%frmin)/(nfreq-1)

   if (nfreq>10 .and. my_rank == master) then
     write(iout, '(a,a,a,a,a,a,a,a)' )&
      ' ddb_diel : the number of frequencies is larger',&
      ' than 10 => I will consider only',ch10,&
      ' the three principal directions, assume that the tensor',ch10,&
      ' is diagonalized, and give dielectric constant and ',ch10,&
      ' reflectivities.'
     write(iout, '(a,a)' )&
      ' Frequency(Hartree)    Dielectric constant   ',&
      '             Reflectivity    '
     write(iout, '(a,a)' )&
      '                     x           y          z',&
      '          x        y        z'
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
!          See Gervais and Piriou PRB11,3944(1975) [[cite:Gervais1975]].
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
    ' Generalized Lyddane-Sachs-Teller relation at zero frequency :',ch10,&
    ' Direction                     Dielectric constant'
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
 end if ! End of the condition of nph2l does not vanish for Lyddane-Sachs-Teller

 ABI_DEALLOCATE(frdiel)
 ABI_DEALLOCATE(modez)
 ABI_DEALLOCATE(oscstr)
 ABI_DEALLOCATE(dielt_modedecompo)

end subroutine ddb_diel
!!***



!!****f* ABINIT/ddb_diel_elec
!!
!! NAME
!! ddb_diel_elec
!!
!! FUNCTION
!! Compute the electronic response of the dielectric constant
!!
!! INPUTS
!! d2cart(2,3,mpert,3,mpert)=
!!  dynamical matrix, effective charges, dielectric tensor,....
!!  all in cartesian coordinates
!! iout=unit number for outputs
!! mpert =maximum number of ipert
!! natom=number of atoms in unit cell
!!
!! OUTPUT
!! epsinf(3,3)= epsilon^infty = electronic contribution to the
!!  dielectric tensor
!! 
!! PARENTS
!!      ddb_diel
!!
!! CHILDREN
!!

subroutine ddb_diel_elec(d2cart,natom,mpert,iout,epsinf)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iout,mpert,natom
!arrays
 real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
 real(dp),intent(out) :: epsinf(3,3)

!Local variables -------------------------
!scalars
 integer :: idir1,idir2
 character(len=500) :: message
!arrays

 write(message, '(a,a)' ) ch10,' Electronic dielectric tensor'
 call wrtout(std_out,message,'COLL')
 call wrtout(iout,message,'COLL')

 !Compute the electronic contribution to the dielectric tensor
 !Needs only the perturbations with E-field from the DDB 
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

end subroutine ddb_diel_elec
!!***



!!****f* ABINIT/ddb_oscstr
!!
!! NAME
!! ddb_oscstr
!!
!! FUNCTION
!! Compute the oscillator strength and the mode effective charge
!!
!! INPUTS
!! displ(2,3*natom,3*natom)=
!!  the displacements of atoms in cartesian coordinates.
!!  The first index means either the real or the imaginary part,
!!  The second index runs on the direction and the atoms displaced
!!  The third index runs on the modes.
!! d2cart(2,3,mpert,3,mpert)=
!!  dynamical matrix, effective charges, dielectric tensor,....
!!  all in cartesian coordinates
!! iout=unit number for outputs
!! mpert =maximum number of ipert
!! natom=number of atoms in unit cell
!! phfrq(3*natom)=phonon frequencies (square root of the dynamical
!!  matrix eigenvalues, except if these are negative, and in this
!!  case, give minus the square root of the absolute value
!!  of the matrix eigenvalues). Hartree units.
!!
!! OUTPUT
!! fact_oscstr(2,3,3*natom)=oscillator strengths for the different eigenmodes,
!!  for different direction of the electric field.
!! modez(2,3,3*natom)=mode effective charges for the different eigenmodes,
!!  for different directions of the electric field, following
!!  the definition Eq.(53) in PRB55, 10355 (1997) [[cite:Gonze1997a]]
!! oscstr(2,3,3,3*natom)=oscillator strengths, following
!!  the definition Eq.(54) in PRB55, 10355 (1997) [[cite:Gonze1997a]]
!! 
!! PARENTS
!!      ddb_diel
!!
!! CHILDREN
!!

subroutine ddb_oscstr(displ,d2cart,fact_oscstr,oscstr,modez,iout,mpert,natom,phfrq,ncid)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iout,mpert,natom,ncid
!arrays
 real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
 real(dp),intent(in) :: phfrq(3*natom)
 real(dp),intent(inout) :: displ(2,3*natom,3*natom)
 real(dp),intent(out) :: fact_oscstr(2,3,3*natom),oscstr(2,3,3,3*natom),modez(2,3,3*natom)

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 integer :: i1,ii,idir1,idir2,imode,ipert1,iphl2,nfreq
 integer :: my_rank,ncerr
 real(dp) :: usquare
 character(len=500) :: message
 logical :: t_degenerate
!arrays
! real(dp),allocatable :: modez(:,:,:),oscstr(:,:,:,:)
 character(len=1),allocatable :: metacharacter(:)

! *********************************************************************

! ABI_ALLOCATE(modez,(2,3,3*natom))
!oscstr(2,3,3,3*natom)=oscillator strengths, following
!the definition Eq.(54) in PRB55, 10355 (1997) [[cite:Gonze1997a]]
! ABI_ALLOCATE(oscstr,(2,3,3,3*natom))


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

    !write(std_out,'(a,i4,a,12es16.6)')'imode=',imode,' displ(:,:,imode)',displ(:,:,imode)
    !write(std_out,'(a,i4,a,6es16.6)')'imode=',imode,' fact_oscstr(:,:,imode)=',fact_oscstr(:,:,imode)

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

     ! write the oscillator strength to the netcdf
#ifdef HAVE_NETCDF
     if (ncid /= nctk_noid) then
       ncerr = nctk_def_arrays(ncid, [nctkarr_t("oscillator_strength", "dp", &
       "complex, number_of_cartesian_directions, number_of_cartesian_directions, number_of_phonon_modes")], &
       defmode=.True.)
       NCF_CHECK(ncerr)
       NCF_CHECK(nctk_set_datamode(ncid))
       NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "oscillator_strength"), oscstr))
     end if
#endif

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
   
end subroutine ddb_oscstr
!!***



!!****f* ABINIT/alignph
!!
!! NAME
!! alignph
!!
!! FUNCTION
!! Construct linear combinations of the phonon eigendisplacements
!! of degenerate modes in order to align the mode effective charges
!! along the axes of the cartesian frame.
!!
!! INPUTS
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! displ(2,3*natom,3*natom)=
!! the displacements of atoms in cartesian coordinates.
!! The first index means either the real or the imaginary part,
!! The second index runs on the direction and the atoms displaced
!! The third index runs on the modes.
!! d2cart(2,3,mpert,3,mpert)=
!!  dynamical matrix, effective charges, dielectric tensor,....
!!  all in cartesian coordinates
!! mpert =maximum number of ipert
!! natom=number of atoms in unit cell
!! ntypat=number of types of atoms
!! phfrq(3*natom)=phonon frequencies (square root of the dynamical
!!  matrix eigenvalues, except if these are negative, and in this
!!  case, give minus the square root of the absolute value
!!  of the matrix eigenvalues). Hartree units.
!! typat(natom)=integer label of each type of atom (1,2,...)
!!
!! OUTPUT
!! displ(2,3*natom,3*natom)=
!!  the displacements of atoms in cartesian coordinates.
!!  The eigendisplacements of degenerate modes have been aligned along
!!  the cartesian axes.
!!
!! PARENTS
!!      ddb_diel
!!
!! CHILDREN
!!
!! SOURCE

subroutine alignph(amu,displ,d2cart,mpert,natom,ntypat,phfrq,typat)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,ntypat
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: amu(ntypat),d2cart(2,3,mpert,3,mpert),phfrq(3*natom)
 real(dp),intent(inout) :: displ(2,3*natom,3*natom)

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 integer :: i1,idir1,idir2,ii,imode,imodex,imodey,imodez,ipert1
 real(dp) :: theta
!arrays
 integer,allocatable :: deg(:)
 real(dp) :: zvec(3,3),zvect(3,3)
 real(dp),allocatable :: modez(:,:,:),modezabs(:),oscstr(:,:,:),vec(:,:),vect(:,:)

! *********************************************************************

!Get the oscillator strength and mode effective charge for each mode
 ABI_ALLOCATE(oscstr,(2,3,3*natom))
 ABI_ALLOCATE(modez,(2,3,3*natom))
 ABI_ALLOCATE(modezabs,(3*natom))
 ABI_ALLOCATE(vec,(3*natom,3))
 ABI_ALLOCATE(vect,(3*natom,3))
 ABI_ALLOCATE(deg,(3*natom))

 write(std_out,'(a,a)')ch10,' alignph : before modifying the eigenvectors, mode number and mode effective charges :'
 do imode=1,3*natom
   modezabs(imode)=zero
   do ii=1,2
     do idir2=1,3
       oscstr(ii,idir2,imode)=zero
       modez(ii,idir2,imode)=zero
       do idir1=1,3
         do ipert1=1,natom
           i1=idir1+(ipert1-1)*3
           oscstr(ii,idir2,imode)=oscstr(ii,idir2,imode)+&
&           displ(ii,i1,imode)*&
&           d2cart(1,idir1,ipert1,idir2,natom+2)
           modez(ii,idir2,imode)=modez(ii,idir2,imode)+&
&           displ(ii,i1,imode)*&
&           d2cart(1,idir1,ipert1,idir2,natom+2)*&
&           sqrt(amu(typat(ipert1))*amu_emass)
         end do
       end do
       if(abs(modez(ii,idir2,imode))>modezabs(imode))modezabs(imode)=abs(modez(ii,idir2,imode))
     end do
   end do
   write(std_out,'(i4,3f16.6)')imode,modez(1,:,imode)
 end do

!Find degenerate modes with non-zero mode effective charge
 imode = 0
 do while (imode < 3*natom)
   imode = imode + 1
   if (imode == 3*natom) then
     deg(imode) = 1
   else if (abs(phfrq(imode) - phfrq(imode+1)) > tol6 .or. modezabs(imode)<tol8 .or. modezabs(imode+1)<tol8) then
!    Differ by phonon frequency or zero mode effective charge
     deg(imode) = 1
   else
     deg(imode) = 2
     if (imode < 3*natom - 1) then
       if (abs(phfrq(imode) - phfrq(imode+2)) < tol6 .and. modezabs(imode+2)>tol8) then
         deg(imode) = 3
         imode = imode + 1
       end if
     end if
     imode = imode + 1
   end if
 end do


!In case of a degenerate mode, with non-zero mode effective charge, align the mode effective charge vector along
!the axes of the cartesian frame
 imode = 1
 do while (imode <= 3*natom)

   write(std_out,'(a,a,i4,a,i2)')ch10,' Mode number ',imode,' has degeneracy ',deg(imode)
   write(std_out,'(a,3es16.6)') ' Mode effective charge of this mode =',modez(1,:,imode)

   if (deg(imode) == 2) then

!    Optimize on the x direction
     write(std_out,'(a,3es16.6)') ' Mode effective charge of next mode =',modez(1,:,imode+1)
     if (abs(modez(1,1,imode)) > tol8) then
       theta = atan(-modez(1,1,imode+1)/modez(1,1,imode))
       vec(:,1) = displ(1,:,imode)
       vec(:,2) = displ(1,:,imode+1)
       displ(1,:,imode) = cos(theta)*vec(:,1) - sin(theta)*vec(:,2)
       displ(1,:,imode+1) = sin(theta)*vec(:,1) + cos(theta)*vec(:,2)
     end if

   else if (deg(imode) == 3) then

     write(std_out,'(a,3es16.6)') ' Mode effective charge of next mode =',modez(1,:,imode+1)
     write(std_out,'(a,3es16.6)') ' Mode effective charge of next-next mode =',modez(1,:,imode+2)

!    Before mixing them, select the mode-effective charge vectors as being predominently "x", "y" or "z" type.
     if(abs(modez(1,1,imode))>abs(modez(1,2,imode))-tol12 .and. &
&     abs(modez(1,1,imode))>abs(modez(1,3,imode))-tol12) then
       imodex=imode
       if(abs(modez(1,2,imode+1))>abs(modez(1,3,imode+1))-tol12)then
         imodey=imode+1 ; imodez=imode+2
       else
         imodez=imode+1 ; imodey=imode+2
       end if
     else if(abs(modez(1,2,imode))>abs(modez(1,1,imode))-tol12 .and. &
&       abs(modez(1,2,imode))>abs(modez(1,3,imode))-tol12) then
       imodey=imode
       if(abs(modez(1,1,imode+1))>abs(modez(1,3,imode+1))-tol12)then
         imodex=imode+1 ; imodez=imode+2
       else
         imodez=imode+1 ; imodex=imode+2
       end if
     else
       imodez=imode
       if(abs(modez(1,1,imode+1))>abs(modez(1,2,imode+1))-tol12)then
         imodex=imode+1 ; imodey=imode+2
       else
         imodey=imode+1 ; imodex=imode+2
       end if
     end if
     vec(:,1)=displ(1,:,imodex)
     vec(:,2)=displ(1,:,imodey)
     vec(:,3)=displ(1,:,imodez)
     zvec(:,1)=modez(1,:,imodex)
     zvec(:,2)=modez(1,:,imodey)
     zvec(:,3)=modez(1,:,imodez)


!    Optimize along x : does the first vector has a component along x ?
     if (abs(zvec(1,1)) > tol8) then
!      Optimize on the (1,2) pair of modes along x
       theta = atan(-zvec(1,2)/zvec(1,1))
       zvect(:,:)=zvec(:,:)
       zvec(:,1) = cos(theta)*zvect(:,1) - sin(theta)*zvect(:,2)
       zvec(:,2) = sin(theta)*zvect(:,1) + cos(theta)*zvect(:,2)
       vect(:,:)=vec(:,:)
       vec(:,1) = cos(theta)*vect(:,1) - sin(theta)*vect(:,2)
       vec(:,2) = sin(theta)*vect(:,1) + cos(theta)*vect(:,2)
!      Optimize on the (1,3) pair of modes along x
       theta = atan(-zvec(1,3)/zvec(1,1))
       zvect(:,:)=zvec(:,:)
       zvec(:,1) = cos(theta)*zvect(:,1) - sin(theta)*zvect(:,3)
       zvec(:,3) = sin(theta)*zvect(:,1) + cos(theta)*zvect(:,3)
       vect(:,:)=vec(:,:)
       vec(:,1) = cos(theta)*vect(:,1) - sin(theta)*vect(:,3)
       vec(:,3) = sin(theta)*vect(:,1) + cos(theta)*vect(:,3)
       if (abs(zvec(2,2)) > tol8) then
!        Optimize on the (2,3) pair of modes along y
         theta = atan(-zvec(2,3)/zvec(2,2))
         zvect(:,:)=zvec(:,:)
         zvec(:,2) = cos(theta)*zvect(:,2) - sin(theta)*zvect(:,3)
         zvec(:,3) = sin(theta)*zvect(:,2) + cos(theta)*zvect(:,3)
         vect(:,:)=vec(:,:)
         vec(:,2) = cos(theta)*vect(:,2) - sin(theta)*vect(:,3)
         vec(:,3) = sin(theta)*vect(:,2) + cos(theta)*vect(:,3)
       end if
!    Likely, the remaining is not needed ... because the vectors have been ordered in x, y, and z major component ...
!    Optimize along x : does the second vector has a component along x ?
     else if(abs(zvec(1,2)) > tol8) then
!      Optimize on the (2,3) pair of modes along x
       theta = atan(-zvec(1,3)/zvec(1,2))
       zvect(:,:)=zvec(:,:)
       zvec(:,2) = cos(theta)*zvect(:,2) - sin(theta)*zvect(:,3)
       zvec(:,3) = sin(theta)*zvect(:,2) + cos(theta)*zvect(:,3)
       vect(:,:)=vec(:,:)
       vec(:,2) = cos(theta)*vect(:,2) - sin(theta)*vect(:,3)
       vec(:,3) = sin(theta)*vect(:,2) + cos(theta)*vect(:,3)
!      Optimize on the (1,3) pair of modes along y
       if (abs(zvec(2,1)) > tol8) then
         theta = atan(-zvec(2,3)/zvec(2,1))
         zvect(:,:)=zvec(:,:)
         zvec(:,1) = cos(theta)*zvect(:,1) - sin(theta)*zvect(:,3)
         zvec(:,3) = sin(theta)*zvect(:,1) + cos(theta)*zvect(:,3)
         vect(:,:)=vec(:,:)
         vec(:,1) = cos(theta)*vect(:,1) - sin(theta)*vect(:,3)
         vec(:,3) = sin(theta)*vect(:,1) + cos(theta)*vect(:,3)
       end if
!    We are left with the pair of vectors (2,3)
     else if (abs(zvec(2,2)) > tol8) then
!      Optimize on the (2,3) pair of modes along y
       theta = atan(-zvec(2,3)/zvec(2,2))
       zvect(:,:)=zvec(:,:)
       zvec(:,2) = cos(theta)*zvect(:,2) - sin(theta)*zvect(:,3)
       zvec(:,3) = sin(theta)*zvect(:,2) + cos(theta)*zvect(:,3)
       vect(:,:)=vec(:,:)
       vec(:,2) = cos(theta)*vect(:,2) - sin(theta)*vect(:,3)
       vec(:,3) = sin(theta)*vect(:,2) + cos(theta)*vect(:,3)
     end if

     displ(1,:,imodex)=vec(:,1)
     displ(1,:,imodey)=vec(:,2)
     displ(1,:,imodez)=vec(:,3)

!    Previous coding, from Marek. Apparently, break the orthogonalization of vectors ...
!    do ii = 1,3
!      coeff(:) = 0._dp
!      if (ii == 1) then
!        jj = 2 ; kk = 3
!      else if (ii == 2) then
!        jj = 1 ; kk = 3
!      else
!        jj = 1 ; kk = 2
!      end if
!      coeff(ii) = 1._dp
!      c1 = modez(1,jj,imode+ii-1)
!      c2 = modez(1,jj,imode+jj-1)
!      c3 = modez(1,jj,imode+kk-1)
!      c4 = modez(1,kk,imode+ii-1)
!      c5 = modez(1,kk,imode+jj-1)
!      c6 = modez(1,kk,imode+kk-1)
!      dtm = c2*c6 - c3*c5
!      if (abs(dtm) > tol8) then
!        coeff(jj) = (c3*c4 - c1*c6)/dtm
!        coeff(kk) = (c1*c5 - c2*c4)/dtm
!      end if
!      mod_ = sqrt(1._dp + coeff(jj)*coeff(jj) + coeff(kk)*coeff(kk))
!      coeff(:) = coeff(:)/mod_
!      displ(1,:,imode+ii-1) = coeff(1)*vec(1,:) + coeff(2)*vec(2,:) + &
!&       coeff(3)*vec(3,:)
!    end do

   end if ! if deg mode

   imode = imode + deg(imode)

 end do

 write(std_out,'(a,a)')ch10,' alignph : after modifying the eigenvectors, mode number and mode effective charges :'
 do imode=1,3*natom
   do ii=1,2
     do idir2=1,3
       modez(ii,idir2,imode)=zero
       do idir1=1,3
         do ipert1=1,natom
           i1=idir1+(ipert1-1)*3
           modez(ii,idir2,imode)=modez(ii,idir2,imode)+&
&           displ(ii,i1,imode)*&
&           d2cart(1,idir1,ipert1,idir2,natom+2)*&
&           sqrt(amu(typat(ipert1))*amu_emass)
         end do
       end do
     end do
   end do
   write(std_out,'(i4,3f16.6)')imode,modez(1,:,imode)
 end do

 ABI_DEALLOCATE(deg)
 ABI_DEALLOCATE(oscstr)
 ABI_DEALLOCATE(modez)
 ABI_DEALLOCATE(modezabs)
 ABI_DEALLOCATE(vec)
 ABI_DEALLOCATE(vect)

end subroutine alignph
!!***

end module m_ddb_diel
!!***
