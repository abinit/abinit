!{\src2tex{textfont=tt}}
!!****f* ABINIT/ddb_compare
!!
!! NAME
!! ddb_compare
!!
!! FUNCTION
!! Compare the temporary DDB and input DDB preliminary information,
!! as well as psp information.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! 1. All the variables have their usual meaning.
!! 2. Note that fullinit==0  means that the input DDB has been
!! initialized by a ground state input file. Some comparison are
!! then not required.
!! 3. All variables with 8 appended are from the new input DDB
!!
!! INPUTS
!!  acell, acell8 = lattice parameters
!!  amu, amu8 = atomic masses
!!  dimekb = dimension of KB projector set (only used for NCPP)
!!  ecut, ecut8 = cutoff energy
!!  ekb, ekb8 = KB energies for pseudopotentials
!!  fullinit, fullmrgddb_init = flags (see notes)
!!  iscf, iscf8 = SCF algorithm
!!  ixc, ixc8 = XC functional
!!  kpt, kpt8 = kpoint array
!!  kptnrm, kptnr8 = normalization factor for kpt
!!  natom, natom8 = number of atoms
!!  nband, nband8 = number of bands at each kpt
!!  ngfft, ngfft8 = FFT grid sizes
!!  nkpt, nkpt8 = number of kpoints
!!  nsppol, nsppo8 = number of spin polarization (1 or 2)
!!  nsym, nsym8 = number of symmetry operations
!!  ntypat, ntypat8 = number of types of atoms
!!  occ, occ8 = occupation numbers
!!  occopt, occop8 = occupation style (metal, insulator, smearing...)
!!  pawecutdg,pawecutdg8= cutoff energy used for the fine "double grid" (PAW only)
!!  pawtab,pawtab8= PAW tabulated data (PAW dataset)
!!  rprim, rprim8 = primitive vectors of unit cell (cartesian coordinates)
!!  dfpt_sciss, dfpt_sciss8 = scissor correction (Ha)
!!  symrel, symre8 = symmetry operations in reciprocal space
!!  tnons, tnons8 = translations associated to symrel
!!  tolwfr, tolwf8 = tolerance on convergence of wavefunctions
!!  typat, typat8 = array of atom types
!!  usepaw = flag for utilization of PAW
!!  wtk, wtk8 = weights of kpoints
!!  xred, xred8 = reduced coordinates of atoms
!!  zion, zion8 = ionic charges of nuclei
!!
!! OUTPUT (corresponding values, checked and/or set)
!!  acell, amu, dimekb, ecut, ekb, fullinit, iscf, ixc, kpt, kptnrm,
!!  natom, nband, ngfft, nkpt, nsppol, nsym, ntypat, occ, occopt,
!!  rprim, dfpt_sciss, symrel, tnons, tolwfr, typat, usepaw, wtk, xred, zion
!!
!! PARENTS
!!      mblktyp1,mblktyp5
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ddb_compare (acell,acell8,amu,amu8,dimekb,ecut,ecut8,ekb,ekb8,&
& fullinit,fullmrgddb_init,iscf,iscf8,ixc,ixc8,kpt,kpt8,kptnrm,kptnr8,&
& natom,natom8,nband,nband8,ngfft,ngfft8,nkpt,nkpt8,&
& nsppol,nsppo8,nsym,nsym8,ntypat,ntypat8,occ,occ8,&
& occopt,occop8,pawecutdg,pawecutdg8,pawtab,pawtab8,&
& rprim,rprim8,dfpt_sciss,dfpt_sciss8,symrel,symre8,&
& tnons,tnons8,tolwfr,tolwf8,typat,typat8,usepaw,wtk,wtk8,xred,xred8,zion,zion8)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_pawtab, only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_compare'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: dimekb,fullmrgddb_init,iscf8,ixc8,natom8,nkpt8,nsppo8,nsym8
 integer,intent(in) :: ntypat8,occop8,usepaw
 integer,intent(inout) :: fullinit,iscf,ixc,natom,nkpt,nsppol,nsym,ntypat
 integer,intent(inout) :: occopt
 real(dp),intent(in) :: ecut8,kptnr8,pawecutdg8,dfpt_sciss8,tolwf8
 real(dp),intent(inout) :: ecut,kptnrm,pawecutdg,dfpt_sciss,tolwfr
!arrays
 integer,intent(in) :: nband8(*),ngfft8(18),symre8(3,3,*),typat8(*)
 integer,intent(inout) :: nband(*),ngfft(18),symrel(3,3,*),typat(*)
 real(dp),intent(in) :: acell8(3),amu8(*),ekb8(dimekb,*),kpt8(3,*),occ8(*)
 real(dp),intent(in) :: rprim8(3,3),tnons8(3,*),wtk8(*),xred8(3,*),zion8(*)
 real(dp),intent(inout) :: acell(3),amu(*),ekb(dimekb,*),kpt(3,*),occ(*)
 real(dp),intent(inout) :: rprim(3,3),tnons(3,*),wtk(*),xred(3,*),zion(*)
 type(pawtab_type),intent(in) :: pawtab8(*)
 type(pawtab_type),intent(inout) :: pawtab(*)

!Local variables -------------------------
!scalars
 integer :: bantot,ii,ij,isym,itypat
 real(dp) :: ekbcm8,ekbcmp 
 real(dp),parameter :: tol=2.0d-14
 character(len=500) :: msg

! *********************************************************************

!DEBUG
!write(std_out,*)' ddb_compare : rprim=',rprim
!write(std_out,*)' ddb_compare : rprim8=',rprim8
!ENDDEBUG

!Compare all the preliminary information
!1. natom
 call chki8(natom,natom8,' natom')
!2. nkpt
!Compares the input and transfer values only if the input has not
!been initialized by a ground state input file
!There can also be the case of perturbation at Gamma, that
!only need half of the number of k points.
 if(fullinit/=0)then
   if(nkpt/=2*nkpt8 .and. 2*nkpt/=nkpt8)then
     call chki8(nkpt,nkpt8,'  nkpt')
   else
     write(std_out,*)' compar8 : assume that one of the DDB to be',&
&     ' merged use Time-Reversal to'
     write(std_out,*)' decrease the number of k-points'
   end if
 else
!  Otherwise, takes the meaningful value
   nkpt=nkpt8
 end if
!3a. occopt
!Because the program will stop if the bloks
!do not compare well, take here the most favorable case.
 if(occop8==0)occopt=0
!3b. nband
!Compares the input and transfer values only if the input has not
!been initialized by a ground state input file
!There can also be the case of perturbation at Gamma, that
!only need half of the number of k points.
 if(fullinit==0 .or. nkpt8==2*nkpt)then
   bantot=0
   do ii=1,nkpt8
     nband(ii)=nband8(ii)
     bantot=bantot+nband(ii)
   end do
 else
   bantot=0
   do ii=1,nkpt
     if(nkpt==nkpt8)then
       call chki8(nband(ii),nband8(ii),' nband')
     end if
     bantot=bantot+nband(ii)
   end do
 end if
!9. nsppol
 call chki8(nsppol,nsppo8,'nsppol')
!4. nsym
 if(nsym/=1 .and. nsym8/=1)then
   call chki8(nsym,nsym8,'  nsym')
 end if
!5. ntypat
 call chki8(ntypat,ntypat8,'ntypat')
!6. acell
 do ii=1,3
   call chkr8(acell(ii),acell8(ii),' acell',tol)
 end do
!7. amu
 do ii=1,ntypat
   call chkr8(amu(ii),amu8(ii),'   amu',tol)
 end do
!9. date
!Compare the two dates, and put in the variable date the most
!recent one.
!(the selection on the year is valid up to 2090 only ....)
!yynow=1991
!yy=date-100*((date-yynow)/100)
!yy8=date8-100*((date8-yynow)/100)
!if(yy<yy8)then
!date=date8
!else if(yy==yy8)then
!mm=date/10000
!mm8=date8/10000
!if(mm<mm8)then
!date=date8
!else if(mm==mm8)then
!dd=date/100-100*(date/10000)
!dd8=date8/100-100*(date8/10000)
!if(dd<dd8)then
!date=date8
!end if
!end if
!end if
!10. ecut
 call chkr8(ecut,ecut8,'  ecut',tol)
!10b. pawecutdg (PAW only)
 if (usepaw==1) then
   call chkr8(pawecutdg,pawecutdg8,'  ecut',tol)
 end if
!11. iscf
!Compares the input and transfer values only if the input has not
!been initialized by a ground state input file
 if(fullinit/=0)then
   call chki8(iscf,iscf8,'  iscf')
 else
!  Otherwise, takes the meaningful value
   iscf=iscf8
 end if
!12. ixc
 call chki8(ixc,ixc8,'   ixc')
!13. kpt and 14. kptnrm
!Compares the input and transfer values only if the input
!has not been initialized by a ground state input file
!and if the number of k points is identical
 if(nkpt8 == 2*nkpt .or. fullinit==0)then
!  Copy the largest number of k points in the right place
   do ij=1,nkpt8
     do ii=1,3
       kpt(ii,ij)=kpt8(ii,ij)
     end do
   end do
   kptnrm=kptnr8
 else if (nkpt==nkpt8)then
   do ij=1,nkpt
     do ii=1,3
!      Compares the input and transfer values only if the input
!      has not been initialized by a ground state input file
       call chkr8(kpt(ii,ij)/kptnrm,&
&       kpt8(ii,ij)/kptnr8,'   kpt',tol)
     end do
   end do
 end if
!16. ngfft
!MT dec 2013: deactivate the stop on ngfft to allow for
! (nfft-converged) DFPT calculations with GS WFK obtained with a different ngfft
 do ii=1,3
   write(msg,'(3a,i10,3a,i10,a)') &
&   'Comparing integers for variable ngfft.',ch10,&
&   'Value from input DDB is',ngfft(ii),' and',ch10,&
&   'from transfer DDB is',ngfft8(ii),'.'
   MSG_WARNING(msg)
!  call chki8(ngfft(ii),ngfft8(ii),' ngfft')
 end do
!17. occ
!Compares the input and transfer values only if the input has not
!been inititialized by a ground state input file
 do ii=1,bantot
   if (fullinit==0 .or. nkpt8==2*nkpt) then
     occ(ii)=occ8(ii)
   else if(nkpt==nkpt8)then
     call chkr8(occ(ii),occ8(ii),'   occ',tol)
   end if
 end do
!18. rprim
 do ii=1,3
   do ij=1,3
     call chkr8(rprim(ii,ij),rprim8(ii,ij),' rprim',tol)
   end do
 end do
!19. dfpt_sciss
!Compares the input and transfer values only if the input has not
!been inititialized by a ground state input file
 if(fullinit/=0)then
   call chkr8(dfpt_sciss,dfpt_sciss8,' dfpt_sciss',tol)
 else
!  Otherwise, takes the meaningful value
   dfpt_sciss=dfpt_sciss8
 end if
!20. symrel
!If nsym == nsym8, compares the symmetry operations,
!otherwise, one of nsym or nsym8 is 1, and thus take the
!symrel corresponding to the largest set.
!nsym will be changed later
 if(nsym==nsym8)then
   do isym=1,nsym
     do ii=1,3
       do ij=1,3
         call chki8(symrel(ii,ij,isym),symre8(ii,ij,isym),'symrel')
       end do
     end do
   end do
 else if(nsym8/=1)then
   symrel(:,:,1:nsym8)=symre8(:,:,1:nsym8)
 end if
!21. tnons (see symrel)
 if(nsym==nsym8)then
   do isym=1,nsym
     do ii=1,3
       call chkr8(tnons(ii,isym),tnons8(ii,isym),' tnons',tol)
     end do
   end do
 else if(nsym8/=1)then
   tnons(:,1:nsym8)=tnons8(:,1:nsym8)
   nsym=nsym8
 end if
!22. tolwfr
!Take the less converged value...
 tolwfr=max(tolwfr,tolwf8)
!23. typat
 do ii=1,ntypat
   call chki8(typat(ii),typat8(ii),' typat')
 end do
!24. wtk
!Compares the input and transfer values only if the input has not
!been initialized by a ground state input file and the
!number of k-points is identical.
 if(nkpt8==2*nkpt .or. fullinit==0)then
   do ii=1,nkpt8
     wtk(ii)=wtk8(ii)
   end do
 else if(nkpt==nkpt8)then
   do ii=1,nkpt
     call chkr8(wtk(ii),wtk8(ii),'   wtk',tol)
   end do
 end if
!25.xred
 do ij=1,natom
   do ii=1,3
     call chkr8(xred(ii,ij),xred8(ii,ij),'  xred',tol)
   end do
 end do
!26. zion
 do ii=1,ntypat
   call chkr8(zion(ii),zion8(ii),'  zion',tol)
 end do

!Finally, put the correct value of nkpt in the case
!of the use of the time-reversal symmetry
 if(2*nkpt==nkpt8)then
   nkpt=nkpt8
 end if

!Now compare the NC pseudopotential information
 if (usepaw==0) then
   if(dimekb/=0 .and. fullinit/=0 .and. fullmrgddb_init/=0 )then
     do ii=1,dimekb
       do itypat=1,ntypat
         ekbcmp=ekb(ii,itypat)
         ekbcm8=ekb8(ii,itypat)
         call chkr8(ekbcmp,ekbcm8,'   ekb',tol)
       end do
     end do
   else if(dimekb/=0 .and. fullmrgddb_init/=0)then
     do ii=1,dimekb
       do itypat=1,ntypat
         ekb(ii,itypat)=ekb8(ii,itypat)
       end do
     end do
   end if
 end if

!Now compare several PAW dataset informations
 if (usepaw==1) then
   if (fullinit/=0 .and. fullmrgddb_init/=0) then
     do itypat=1,ntypat
       call chki8(pawtab(itypat)%basis_size,pawtab8(itypat)%basis_size,'bas_sz')
       call chki8(pawtab(itypat)%lmn_size,pawtab8(itypat)%lmn_size,'lmn_sz')
       call chki8(pawtab(itypat)%lmn2_size,pawtab8(itypat)%lmn2_size,'lmn2sz')
       call chkr8(pawtab(itypat)%rpaw,pawtab8(itypat)%rpaw,'  rpaw',tol3)
       call chkr8(pawtab(itypat)%rshp,pawtab8(itypat)%rshp,'rshape',tol3)
       call chki8(pawtab(itypat)%shape_type,pawtab8(itypat)%shape_type,'shp_tp')
       if (pawtab(itypat)%lmn2_size>0) then
         do ii=1,pawtab(itypat)%lmn2_size
           call chkr8(pawtab(itypat)%dij0(ii),pawtab8(itypat)%dij0(ii),'  dij0',tol)
         end do
       end if
     end do
   else if (fullmrgddb_init/=0) then
     do itypat=1,ntypat
       pawtab(itypat)%basis_size =pawtab8(itypat)%basis_size
       pawtab(itypat)%lmn_size   =pawtab8(itypat)%lmn_size
       pawtab(itypat)%rpaw       =pawtab8(itypat)%rpaw
       pawtab(itypat)%rshp       =pawtab8(itypat)%rshp
       pawtab(itypat)%shape_type =pawtab8(itypat)%shape_type
       if (pawtab8(itypat)%lmn2_size>0) then
         if (pawtab(itypat)%lmn2_size==0)  then
           ABI_ALLOCATE(pawtab(itypat)%dij0,(pawtab8(itypat)%lmn2_size))
         end if
         do ii=1,pawtab8(itypat)%lmn2_size
           pawtab(itypat)%dij0(ii)=pawtab8(itypat)%dij0(ii)
         end do
       end if
       pawtab(itypat)%lmn2_size  =pawtab8(itypat)%lmn2_size
     end do
   end if
 end if

 contains 
!!***

!!****f* ABINIT/chkr8
!!
!! NAME
!! chkr8
!!
!! FUNCTION
!! This small subroutine check the identity of reali and realt,
!! who are integers, and eventually send a message and stop
!! if they are found unequal by more than tol
!!
!! INPUTS
!! reali=first real number
!! intt=second  real number
!! character(len=6) name=name of the variable in the calling routine, to be echoed
!! tol=tolerance
!!
!! OUTPUT
!!  (only checking)
!!
!! PARENTS
!!      ddb_compare
!!
!! CHILDREN
!!
!! SOURCE

subroutine chkr8(reali,realt,name,tol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkr8'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 real(dp),intent(in) :: reali,realt,tol
 character(len=6),intent(in) :: name

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *********************************************************************

   if(abs(reali-realt)>tol) then
     write(message, '(a,a,a,a,a,es16.6,a,a,a,es16.6,a,a,a)' )&
&     'Comparing reals for variable',name,'.',ch10,&
&     'Value from input DDB is',reali,' and',ch10,&
&     'from transfer DDB is',realt,'.',ch10,&
&     'Action: check your DDBs.'
     MSG_ERROR(message)
   end if

 end subroutine chkr8
!!***

!!****f* ABINIT/chki8
!!
!! NAME
!! chki8
!!
!! FUNCTION
!! This small subroutine check the identity of inti and intt,
!! who are integers, and eventually send a message and stop
!! if they are found unequal
!!
!! INPUTS
!! inti=first integer
!! intt=second integer
!! character(len=6) name=name of the variable in the calling routine, to be echoed
!!
!! OUTPUT
!!  (only checking)
!!
!! PARENTS
!!      ddb_compare
!!
!! CHILDREN
!!
!! SOURCE

subroutine chki8(inti,intt,name)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chki8'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: inti,intt
 character(len=6),intent(in) :: name

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *********************************************************************

   if(inti/=intt) then
     write(message, '(a,a,a,a,a,i10,a,a,a,i10,a,a,a)' )&
&     'Comparing integers for variable',name,'.',ch10,&
&     'Value from input DDB is',inti,' and',ch10,&
&     'from transfer DDB is',intt,'.',ch10,&
&     'Action: check your DDBs.'
     MSG_ERROR(message)
   end if

 end subroutine chki8
!!***

end subroutine ddb_compare
!!***
