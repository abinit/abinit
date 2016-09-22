!{\src2tex{textfont=tt}}
!!****f* ABINIT/mlwfovlp_setup
!! NAME
!! mlwfovlp_setup
!!
!! FUNCTION
!! Routine which creates table g1 and ovikp  necessary to compute
!! overlap for Wannier code (www.wannier.org f90 version).
!!
!! COPYRIGHT
!! Copyright (C) 2005-2016 ABINIT group (BAmadon,FJollet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atom_symbols(natom)= table of symbol for each atom
!!                                          and each |p_lmn> non-local projector
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  filew90_win(nsppol) secondary input files for w90
!!  lwanniersetup= flag: only 1 is fully working.
!!  natom              =number of atoms in cell.
!!  mband=maximum number of bands
!!  natom=number of atoms in cell.
!!  nkpt=number of k points.
!!  num_bands(isppol)=number of bands actually used to construct the wannier function
!!  nwan(isppol)= number of wannier fonctions (read in wannier90.win).
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  real_lattice(3,3)=dimensional primitive translations for real space
!!                 in format required by wannier90
!!  recip_lattice(3,3)=dimensional primitive translations for reciprocal space
!!                 in format required by wannier90
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  seed_name=character string for generating wannier90 filenames
!!  spin = just used for nsppol>1 ; 0 both, 1 just spin up, 2 just spin down
!!  xcart(3,natom)=atomic coordinates in bohr
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  band_in(mband,nsppol)   = band to take into account for wannier calculation
!!  g1(3,nkpt,nntot) = G vector shift which is necessary to obtain k1+b
!!                     from k2 in the case where k1+b does not belong to the 1st BZ.
!!  nband_inc(nsppol) = # of included bands
!!  nntot            = number of k-point neighbour
!!  ovikp(nkpt,nntot)= gives  nntot value of k2 (in the BZ) for each k1  (k2=k1+b mod(G))
!!  
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      mlwfovlp
!!
!! CHILDREN
!!      atomdata_from_znucl,wannier_setup,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


 subroutine mlwfovlp_setup(atom_symbols,band_in,dtset,filew90_win,gamma_only,&
& g1,lwanniersetup,mband,natom,nband_inc,nkpt,&
& nntot,num_bands,num_nnmax,nsppol,nwan,ovikp,&
& proj_l,proj_m,proj_radial,proj_site,proj_x,proj_z,proj_zona,&
& real_lattice,recip_lattice,rprimd,seed_name,spin,spinors,xcart,xred)

 use defs_basis
 use defs_abitypes
 use defs_wannier90
 use m_errors
 use m_profiling_abi
 use m_atomdata

 use m_io_tools,   only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mlwfovlp_setup'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments---------------------------
! scalars
!scalars
 integer,intent(in) :: lwanniersetup,mband,natom,nkpt,nsppol
 integer,intent(in) :: num_nnmax,spin
 integer,intent(out) :: nband_inc(nsppol),nntot,num_bands(nsppol),nwan(nsppol)
 logical,intent(in) :: gamma_only,spinors
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(out) :: g1(3,nkpt,num_nnmax),ovikp(nkpt,num_nnmax)
 integer,intent(out) :: proj_l(mband,nsppol),proj_m(mband,nsppol),proj_radial(mband,nsppol)
 real(dp),intent(in) :: real_lattice(3,3)
 real(dp),intent(in) :: recip_lattice(3,3),rprimd(3,3),xred(3,natom)
 real(dp),intent(out) :: proj_site(3,mband,nsppol),proj_x(3,mband,nsppol),proj_z(3,mband,nsppol)
 real(dp),intent(out) :: proj_zona(mband,nsppol),xcart(3,natom)
 logical,intent(out) :: band_in(mband,nsppol)
 character(len=3),intent(out) :: atom_symbols(natom)
 character(len=fnlen),intent(in) :: seed_name(nsppol),filew90_win(nsppol)

!Local variables---------------------------
!scalars
 integer :: iatom,icb,ikpt,ikpt1,intot,isppol,itypat,jj,mband_,unt
 real(dp) :: znucl1
 character(len=2) :: symbol
 character(len=500) :: message
 character(len=fnlen) :: filew90_nnkp
 type(atomdata_t) :: atom
!arrays
 integer :: exclude_bands(mband,nsppol),ngkpt(3)

! *************************************************************************

!^^^^^^^^^^^^^^^^read wannier90.nnkp^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 if(lwanniersetup==0) then  !this part is not coded for nsppol>1
   isppol=1
   filew90_nnkp=trim(seed_name(isppol))//'.nnkp'
   if (open_file(filew90_nnkp,message,newunit=unt,form='formatted',status='old') /= 0) then
     MSG_ERROR(message)
   end if
   read(unt,*)
   read(unt,*) nntot , mband_, nwan(1)
   write(message, '(a,a,i6,i6,i6)' )ch10,&
&   ' mlwfovlp_setup nntot,mband,nwan ', nntot,mband_,nwan(1)
   call wrtout(std_out,message,'COLL')
   if(mband_.ne.mband) then
     write(message, '(4a)' )&
&     'mband_ is not equal to mband ',ch10,&
&     'Action: check ',trim(filew90_nnkp)
     MSG_ERROR(message)
   end if
   if(nwan(1)>mband) then
     write(message, '(4a)' )&
&     'nwan > mband ',ch10,&
&     'Action: check ',trim(filew90_nnkp)
     MSG_ERROR(message)
   end if
   if(nwan(1)==0) then
     write(message, '(4a)' )&
&     'nwan = 0 ',ch10,&
&     'Action: check ',trim(filew90_nnkp)
     MSG_ERROR(message)
   end if
   do ikpt=1,nkpt
     do intot=1,nntot
!      ikpt1: k point  (ikpt=ikpt1)
!      ovikp(intot,ikpt): neighbour number intot for ikpt
!      g1(1:3,intot,ikpt): non reciprocal space vector between the 2 k-points
       read(unt,*)  &
&       ikpt1,ovikp(ikpt,intot),(g1(jj,ikpt,intot),jj=1,3)
       if(ikpt1.ne.ikpt) write(std_out,*) "warning: ikpt1 .ne ikpt : ?"
     end do
   end do
   close(unt)
   write(message, '(3a)' )ch10,&
&   trim(filew90_nnkp),'wannier90.nnkp has been read !'
   call wrtout(std_out,message,'COLL')

   message = ' exclude bands is not given in this case (not implemented) '
   MSG_ERROR(message)

!  ^^^^^^^^^^^^^^^^^^^^^^^ call wannier_setup begin^^^^^^^^^^^^^^^^^^^^^^^^
 else if (lwanniersetup==1) then
   num_bands(:)=mband
!  num_nnmax=12 !limit fixed for compact structure in wannier_setup.
   ovikp=0.d0
!  "When nshiftk=1, kptrlatt is initialized as a diagonal (3x3) matrix, whose diagonal 
!  elements are the three values ngkpt(1:3)"
   ngkpt(1)=dtset%kptrlatt(1,1)
   ngkpt(2)=dtset%kptrlatt(2,2) !  have to verif kptrlatt is diagonal
   ngkpt(3)=dtset%kptrlatt(3,3)
   do iatom=1,natom
     itypat=dtset%typat(iatom)
     znucl1=dtset%znucl(itypat)
     call atomdata_from_znucl(atom, znucl1) 
     symbol=trim(adjustl(atom%symbol))
!    write(309,*) symbol
     atom_symbols(iatom)=symbol
     xcart(:,iatom)=rprimd(:,1)*xred(1,iatom)+&
&     rprimd(:,2)*xred(2,iatom)+&
&     rprimd(:,3)*xred(3,iatom)
   end do ! iatom
!  write(std_out,*) xcart
!  write(std_out,*) Bohr_Ang
!  write(std_out,*) rprimd*Bohr_Ang
!  write(std_out,*) seed_name
!  write(std_out,*) ngkpt
!  write(std_out,*) nkpt
!  write(std_out,*) mband
!  write(std_out,*) natom
!  write(std_out,*) atom_symbols
   write(message, '(a,a)' )ch10,&
&   '** mlwfovlp_setup:  call wannier90 library subroutine wannier_setup'
   call wrtout(std_out,message,'COLL')
#if defined HAVE_WANNIER90
   nwan(:)=0
   num_bands(:)=0
   do isppol=1,nsppol
     if(spin.ne.0 .and. spin.ne.isppol) cycle
     call wannier_setup(seed_name(isppol),ngkpt,nkpt&                    !input
&    ,real_lattice,recip_lattice,dtset%kpt&                      !input
&    ,mband,natom,atom_symbols,xcart*Bohr_Ang&                   !input
&    ,gamma_only,spinors&                                        !input
&    ,nntot,ovikp,g1,num_bands(isppol),nwan(isppol)&                     !output
&    ,proj_site(:,:,isppol),proj_l(:,isppol)&                    !output
&    ,proj_m(:,isppol),proj_radial(:,isppol)&                    !output
&    ,proj_z(:,:,isppol),proj_x(:,:,isppol)&                     !output
&    ,proj_zona(:,isppol),exclude_bands(:,isppol) )               !output
   end do !isppol


#endif
!  do isppol=1,nsppol
!  if(spin.ne.0 .and. spin.ne.isppol) cycle
!  write(std_out,*)  "1", nntot,nwan(isppol)
!  write(std_out,*)  "2",num_bands(isppol)  ! states on which wannier functions are computed
!  write(std_out,*)  "3", proj_site(:,1:nwan(isppol),isppol)
!  write(std_out,*)  "4",proj_l(1:nwan(isppol),isppol)
!  write(std_out,*)  "5",proj_m(1:nwan(isppol),isppol)
!  write(std_out,*)  "6", proj_radial(1:nwan(isppol),isppol)
!  write(std_out,*)  "7", proj_z(:,1:nwan(isppol),isppol)
!  write(std_out,*)  "8", proj_x(:,1:nwan(isppol),isppol)
!  write(std_out,*)  "9",proj_zona(1:nwan(isppol),isppol)
!  write(std_out,*)  "10",exclude_bands(:,isppol)
!  end do!isppol
!  testdebug:  ovikp(1,1)=1
!  ^^^^^^^^^^^^^^^^^^^^^^^ end ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 end if  ! lwanniersetup
 do isppol=1,nsppol
   if(spin.ne.0 .and. spin.ne.isppol) cycle
   band_in(:,isppol)=.true.
   do icb=1,mband
     if(exclude_bands(icb,isppol).ne.0)  band_in(exclude_bands(icb,isppol),isppol)=.false. 
   end do
   nband_inc(isppol)=0
   do icb=1, mband
     if (band_in(icb,isppol)) then
       nband_inc(isppol)=nband_inc(isppol)+1
     end if
   end do
 end do !isppol
 if(any(mband.gt.num_bands(:))) then
   write(message, '(a,a)' )ch10,&
&   '   Following bands are excluded from the calculation of wannier functions:'
   call wrtout(std_out,message,'COLL')
   
   do isppol=1,nsppol
     if(spin.ne.0 .and. spin.ne.isppol) cycle
     if(nsppol==2) then
       write(message,'("For spin",i4)')isppol
!      write(message,'(a,i)')'For spin=',isppol
       call wrtout(std_out,message,'COLL')
     end if !nsppol
     do jj=1,mband-num_bands(isppol),10
       write(message,'(10i7)') exclude_bands(jj:min(jj+9,mband-num_bands(isppol)),isppol)
       call wrtout(std_out,message,'COLL')
     end do
   end do !isppol
 end if

 do isppol=1,nsppol
   if(spin.ne.0 .and. spin.ne.isppol) cycle
   if(nsppol==2) then
     write(message,'("For spin",i4)')isppol
     call wrtout(std_out,message,'COLL')
   end if !nsppol
   write(message, '(a,i6,3a)' )ch10,&
&   nwan(isppol),' wannier functions will be computed (see ',trim(filew90_win(isppol)),')'
   call wrtout(std_out,message,'COLL') 
!  write(std_out,*) exclude_bands(icb),band_in(icb)
!  ^^^^^^^^^^^^^^^END OF READING
   write(message, '(a,i6,a)' )ch10,&
&   num_bands(isppol),' bands will be used to extract wannier functions'
   call wrtout(std_out,message,'COLL')
   if(num_bands(isppol).lt.nwan(isppol)) then
     write(message, '(4a)' )&
&     ' number of bands is lower than the number of wannier functions',ch10,&
&     ' Action : check input file and ',trim(filew90_win(isppol))
     MSG_ERROR(message)
   else if (num_bands(isppol)==nwan(isppol)) then
     write(message, '(a,a,a,a)' )ch10,&
&     '   Number of bands is equal to the number of wannier functions',ch10,&
&     '   Disentanglement will not be necessary'
     call wrtout(std_out,message,'COLL')
   else if  (num_bands(isppol).gt.nwan(isppol)) then
     write(message, '(a,a,a,a)' )ch10,&
&     '   Number of bands is larger than the number of wannier functions',ch10,&
&     '   Disentanglement will be necessary'
     call wrtout(std_out,message,'COLL')
   end if
   write(message, '(2x,a,a,i3,1x,a)' )ch10,&
&   '   Each k-point has', nntot,'neighbours'
   call wrtout(std_out,message,'COLL')

 end do !isppol

end subroutine mlwfovlp_setup
!!***
