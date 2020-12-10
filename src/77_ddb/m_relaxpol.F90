!!****m* ABINIT/m_relaxpol
!! NAME
!!  m_relaxpol
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2020 ABINIT group (MVeithen)
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

module m_relaxpol

 use defs_basis
 use m_abicore
 use m_errors

 use m_fstrings,  only : sjoin, itoa
 use m_symtk,     only : matr3inv
 use m_berrytk,   only : polcart
 use m_hide_lapack,   only : dzgedi, dzgefa
 use m_geometry,  only : xcart2xred
 use m_dynmat,    only : symdyma
 use m_crystal,   only : crystal_t

 implicit none

 private
!!***

 public :: relaxpol
!!***

contains
!!***

!!****f* ABINIT/relaxpol
!! NAME
!! relaxpol
!!
!! FUNCTION
!! 1) Compute polarization in cartesian coordinates
!! 2) Structural relaxation at fixed polarization: this routine
!!    solves the linear system of equations Eq.(13)
!!    of Na Sai et al., PRB 66, 104108 (2002) [[cite:Sai2002]].
!!
!! INPUTS
!! blkflg(msize) = flag for every matrix element (0=> the element
!!   is not in the data block), (1=> the element is in the data blok)
!! blkval(2,msize) = matrix that contains the second-order energy derivatives
!! etotal = Kohn-Sham energy at zero electric field
!! fred(3,natom) = -1 times the forces in reduced coordinates
!! iatfix(natom) = indices of the atoms that are held fixed in the relaxation
!! iout = unit number for output
!! istrfix(6) = indices of the elements of the strain tensor that
!!   are held fixed in the relaxation
!!      1 = xx
!!      2 = yy
!!      3 = zz
!!      4 = yz & zy
!!      5 = xz & zx
!!      6 = xy & yx
!! mpert = maximum number of ipert
!! msize = dimension of blkflg and blkval
!! natfix = number of atoms that are held fixed in the relaxation
!! natom = number of atoms in the unit cell
!! nstrfix = number of elements of the strain tensor that are held fixed in the relaxation
!! pel(3) = electronic polarization not taking into account the factor 1/ucvol
!!  red_ptot(3) = total polarization reduced units   !!REC
!! relaxat = 1: relax atomic positions
!!         = 0: do not relax atomic positions
!! relaxstr = 1: relax cell parameters
!!          = 0: do not relax cell parameters
!! strten(6) = stress tensor in cartesian coordinates
!! targetpol(3) = target value of the polarization
!! usepaw = 1 if PAW in use, 0 else (needed for polcart)
!!
!! OUTPUT
!!
!! NOTES
!! - The elements of the dynamical matrix stored in blkval
!!   are symmetrized before computing the new atomic positions and cell parameters.
!! - In case relaxat = 0 and relaxstr = 0, the routine only
!!   computes the polarization in cartesian coordinates.
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      dzgedi,dzgefa,matr3inv,polcart,symdyma,xcart2xred
!!
!! SOURCE

subroutine relaxpol(Crystal,blkflg,blkval,etotal,fred,iatfix,iout,istrfix,&
& mpert,msize,natfix,natom,nstrfix,pel,red_ptot,relaxat,relaxstr,&
& strten,targetpol,usepaw)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iout,mpert,msize,natfix,natom,nstrfix
 integer,intent(in) :: relaxat,relaxstr,usepaw
 real(dp),intent(in) :: etotal
 type(crystal_t),intent(in) :: Crystal
!arrays
 integer,intent(in) :: blkflg(msize),iatfix(natom)
 integer,intent(in) :: istrfix(6)
 real(dp),intent(in) :: fred(3,natom),pel(3),strten(6)
 real(dp),intent(in) :: red_ptot(3)
 real(dp),intent(inout) :: blkval(2,msize),targetpol(3)

!Local variables -------------------------
!scalars
 integer :: flag,iatom,idir,ii,index,index1,index_tild,info,ipert,istrain
 integer :: itypat,jdir,job,jpert,polunit,posi,posj,sizef
 logical :: iwrite
 real(dp) :: e1,fmax,poltmp,sigmax,tol,value,ucvol
 character(len=500) :: message
!arrays
 integer :: irelaxstrain(6)
 integer,allocatable :: ipvt(:),irelaxat(:),rfpert(:,:)
 real(dp) :: acell_new(3),delta_eta(6),delta_xcart(3,natom),det(2,2),diffpol(3),rprimd(3,3)
 real(dp) :: diffsig(6),favg(3),gprimd(3,3),lambda(3),pel_cart(3),pelev(3)
 real(dp) :: pion(3),pion_cart(3),ptot_cart(3),qphon(3),rprim(3,3)
 real(dp) :: rprimd_new(3,3),sigelfd(6),strainmat(3,3)
 real(dp) :: xcart_new(3,natom),xred_new(3,natom)
 real(dp),allocatable :: cfac(:,:),delta(:),dymati(:),fcart(:,:),fcmat(:,:,:)
 real(dp),allocatable :: fdiff(:,:),felfd(:,:),ifcmat(:,:,:),vec(:),zgwork(:,:)

! *********************************************************************

 rprimd = Crystal%rprimd
 ucvol = Crystal%ucvol
 iwrite = iout > 0

!Check if some degrees of freedom remain fixed during the optimization

 ABI_MALLOC(irelaxat,(natom))
 irelaxat(:) = 1   ; irelaxstrain(:) = 1
 if (natfix > 0) then
   do ii = 1, natfix
     iatom = iatfix(ii)
     if ((iatom > natom).or.(iatom < 0)) then
       write(message, '(a,i0,a,i0,a,a,a,a,a)')&
&       'The value of iatfix(',ii,') is ',iatom,', which is not allowed.',ch10,&
&       'iatfix must be larger than 0 and smaller than natom.',ch10,&
&       'Action: correct iatfix in your input file.'
       ABI_ERROR(message)
     end if
     irelaxat(iatom) = 0
   end do
 end if

 if (nstrfix > 0) then
   do ii = 1, nstrfix
     istrain = istrfix(ii)
     if ((istrain > 6).or.(istrain < 0)) then
       write(message, '(a,i0,a,i0,a,a,a,a,a)')&
&       'istrfix(',ii,') is',istrain,', which is not allowed.',ch10,&
&       'istrfix must be larger than 0 and smaller than 6.',ch10,&
&       'Action : correct istrfix in your input file.'
       ABI_ERROR(message)
     end if
     irelaxstrain(istrain) = 0
   end do
 end if


 ABI_MALLOC(rfpert,(mpert,3))
 ABI_MALLOC(cfac,(mpert,mpert))
 call matr3inv(rprimd,gprimd)

!Compute the size of the matrix that contains the second-order derivatives

 sizef = 3
 rfpert(:,:) = 0
 rfpert(natom+2,1:3) = 1
 if (relaxat == 1) then
   do iatom = 1, natom
     if (irelaxat(iatom) == 1) then
       sizef = sizef + 3
       rfpert(iatom,1:3) = 1
     end if
   end do
 end if
 ii = natom + 2
 if (relaxstr == 1) then
   istrain = 0
   do ipert = (natom+3), (natom+4)
     do idir = 1, 3
       istrain = istrain + 1
       if (irelaxstrain(istrain) == 1) then
         sizef = sizef + 1
         rfpert(ipert,idir) = 1
       end if
     end do
   end do
 end if

 ABI_MALLOC(fcmat,(2,sizef,sizef))
 ABI_MALLOC(ifcmat,(2,sizef,sizef))
 ABI_MALLOC(vec,(sizef))
 ABI_MALLOC(delta,(sizef))
 ABI_MALLOC(ipvt,(sizef))
 ABI_MALLOC(zgwork,(2,sizef))
 ABI_MALLOC(fcart,(3,natom))
 ABI_MALLOC(felfd,(3,natom))
 ABI_MALLOC(fdiff,(3,natom))

!Build the vector that stores the forces, sigma and the polarization

 vec(:) = zero
 posi = 0

 if (relaxat == 1) then

!  Note conversion to cartesian coordinates (bohr) AND
!  negation to make a force out of a gradient
!  Also subtract off average force from each force
!  component to avoid spurious drifting of atoms across cell.
   favg(:) = zero
   do iatom = 1, natom
     do idir = 1, 3
       fcart(idir,iatom) = -(gprimd(idir,1)*fred(1,iatom) + &
&       gprimd(idir,2)*fred(2,iatom) + &
&       gprimd(idir,3)*fred(3,iatom))
       favg(idir) = favg(idir) + fcart(idir,iatom)
     end do
   end do
   favg(:) = favg(:)/dble(natom)
   do iatom = 1, natom
     fcart(:,iatom) = fcart(:,iatom) - favg(:)
   end do

   do iatom = 1, natom
     if (irelaxat(iatom) == 1) then
       do idir = 1, 3
         posi = posi + 1
         vec(posi) = fcart(idir,iatom)
       end do
     end if
   end do

 end if    ! relaxat == 1

!DEBUG
!write(std_out,*)'Forces in cartesian coords'
!do iatom = 1, natom
!write(std_out,'(3(2x,e16.9))')(fcart(idir,iatom),idir = 1, 3)
!end do
!stop
!ENDDEBUG

!Transform target polarization to atomic units
 targetpol(:) = targetpol(:)*((Bohr_Ang*1.0d-10)**2)/e_Cb

!Compute ionic polarization
 pion(:) = zero
 do iatom = 1, natom
   itypat = Crystal%typat(iatom)
   do idir = 1, 3
     poltmp = Crystal%zion(itypat) * Crystal%xred(idir,iatom)
     poltmp = poltmp - two*nint(poltmp/two)   ! fold into [-1,1]
     pion(idir) = pion(idir) + poltmp
   end do
 end do
 do idir = 1, 3
   pion(idir) = pion(idir) - two*nint(pion(idir)/two) ! fold into [-1,1]
 end do

!Transform the polarization to cartesian coordinates
 polunit = 3
 pelev=zero
 call polcart(red_ptot,pel,pel_cart,pelev,pion,pion_cart,polunit,&
& ptot_cart,rprimd,ucvol,iout,usepaw)

 do idir = 1, 3
   posi = posi + 1
   vec(posi) = ptot_cart(idir) - targetpol(idir)
 end do


 if (relaxstr == 1) then
   do istrain = 1, 6
     if (irelaxstrain(istrain) == 1) then
       posi = posi + 1
       vec(posi) = -1._dp*strten(istrain)*ucvol
     end if
   end do
 end if


!Symmetrize the dynamical matrix

 ABI_MALLOC(dymati,(2*3*natom*3*natom))
!by the symdyma routine
 do ipert = 1, natom
   do idir = 1, 3
     do jpert = 1, natom
       do jdir = 1, 3
         index  = jdir +3*((jpert - 1) + mpert*((idir - 1) + 3*(ipert - 1)))
         index1 = jdir +3*((jpert - 1) + natom*((idir - 1) + 3*(ipert - 1)))
         dymati(2*index1 - 1) = blkval(1,index)
         dymati(2*index1    ) = blkval(2,index)
       end do
     end do
   end do
 end do

 qphon(:) = zero
 call symdyma(dymati,Crystal%indsym,natom,Crystal%nsym,qphon,rprimd,Crystal%symrel,Crystal%symafm)

 do ipert = 1, natom
   do idir = 1, 3
     do jpert = 1, natom
       do jdir = 1, 3
         index  = jdir +3*((jpert - 1) + mpert*((idir - 1) + 3*(ipert - 1)))
         index1 = jdir +3*((jpert - 1) + natom*((idir - 1) + 3*(ipert - 1)))
         blkval(1,index) = dymati(2*index1 - 1)
         blkval(2,index) = dymati(2*index1    )
       end do
     end do
   end do
 end do

 ABI_FREE(dymati)

!Define conversion factors for blkval
 cfac(:,:) = 1._dp
 cfac(1:natom,natom+2) = -1._dp/ucvol
 cfac(natom+2,1:natom) = -1._dp/ucvol
 cfac(natom+3:natom+4,natom+2) = -1._dp
 cfac(natom+2,natom+3:natom+4) = -1._dp


!Build the matrix that contains the second-order derivatives
!ipert = natom + 1 corresponds to the ddk perturbation, that
!is not needed; so skip it

 fcmat(:,:,:) = zero

 posi = 0
 flag = 0
! When fcmat has been build, flag = 0 if all elements were available.
! Otherwise, it will be 1. In case one element is missing, check if
! it can be obtained by changing the order of the perturbations

 do ipert = 1, mpert
   do idir = 1, 3
     if (rfpert(ipert,idir) == 1) then
       posi = posi + 1
       posj = 0

       do jpert = 1, mpert
         do jdir = 1, 3
           if (rfpert(jpert,jdir) == 1) then
             index = jdir +3*((jpert - 1) + mpert*((idir - 1) + 3*(ipert - 1)))
             index_tild = idir +3*((ipert - 1) + mpert*((jdir - 1) + 3*(jpert - 1)))
             posj = posj + 1
             if ((ipert /= natom + 2).or.(jpert /= natom + 2)) then
               if (blkflg(index) == 1) then
                 fcmat(:,posi,posj) = blkval(:,index)*cfac(ipert,jpert)
               else if (blkflg(index_tild) == 1) then
                 fcmat(:,posi,posj) = blkval(:,index_tild)*cfac(ipert,jpert)
                 blkval(:,index) = blkval(:,index_tild)
               else
                 flag = 1
                 write(std_out,'(a,4(2x,i3))')'relaxpol: could not find element:',idir,ipert,jdir,jpert
               end if
             end if
!            DEBUG
!            write(100,'(4(2x,i3),5x,f16.9)')idir,ipert,jdir,jpert,fcmat(1,posi,posj)
!            ENDDEBUG
           end if
         end do
       end do

     end if
   end do
 end do

 if (flag == 1) then
   write(message, '(a,a,a,i0,a,i0,a,a,a,a)' )&
&   'Some of the second order derivatives required to deal with the case',ch10,&
&   'relaxat = ',relaxat,', relaxstr = ', relaxstr, ch10,&
&   'are missing in the DDB.',ch10,&
&   'Action: correct your DDB or change your input file.'
   ABI_ERROR(message)
 end if


!Compute the inverse of the force constant matrix

 if ((relaxat /= 0).or.(relaxstr /= 0)) then

   job = 1          ! compute inverse only
   ifcmat(:,:,:) = fcmat(:,:,:)

   call dzgefa(ifcmat,sizef,sizef,ipvt,info)
   ABI_CHECK(info == 0, sjoin("dzgefa returned:", itoa(info)))
   call dzgedi(ifcmat,sizef,sizef,ipvt,det,zgwork,job)

!  DEBUG
!  write(100,*)'relaxat = ',relaxat
!  write(100,*)'relaxstr = ',relaxstr
!  write(100,*)'irelaxat = '
!  write(100,*)irelaxat(:)
!  write(100,*)'irelaxstrain = '
!  write(100,*)irelaxstrain(:)
!  write(100,*)'sizef = ',sizef
!  write(100,*)'targetpol ='
!  write(100,*)targetpol(:)
!  do ipert = 1, sizef
!  do jpert = 1, sizef
!  write(100,'(2(2x,i3),2x,e16.9)')ipert,jpert,fcmat(1,ipert,jpert)
!  end do
!  end do
!  stop
!  ENDDEBUG

!  Compute \delta R, \delta \eta and \lambda
   delta(:) = zero
   do ipert = 1, sizef
     do jpert = 1, sizef
       delta(ipert) = delta(ipert) + ifcmat(1,ipert,jpert)*vec(jpert)
     end do
   end do


!  Update atomic positions
   posi = 0
   if (relaxat == 1) then

     delta_xcart(:,:) = zero
     xcart_new(:,:) = zero
     do iatom = 1, natom
       if (irelaxat(iatom) == 1) then
         do idir = 1, 3
           posi = posi + 1
           delta_xcart(idir,iatom) = delta(posi)
         end do
       end if
     end do

!    Drop unsignificant digits in order to eleminate numerical noise
     tol = 10000000._dp
     do iatom = 1, natom
       do idir = 1, 3
         value = delta_xcart(idir,iatom)
         ii = log10(abs(value))
         if (ii <= 0) then
           ii = abs(ii) + 1
           value = one*int(tol*value*10.0_dp**ii)/(tol*10.0_dp**ii) !vz_d
         else
           value = one*int(tol*value/(10.0_dp**ii))*(10.0_dp**ii)/tol !vz_d
         end if
         delta_xcart(idir,iatom) = value
       end do
     end do

     xcart_new(:,:) = Crystal%xcart(:,:) + delta_xcart(:,:)
     call xcart2xred(natom,rprimd,xcart_new,xred_new)
   end if  ! relaxat == 1

! Compute lambda and the value of the energy functional F - \lambda \cdot P$

   e1 = etotal
   do idir = 1, 3
     posi = posi + 1
     lambda(idir) = delta(posi)
     e1 = e1 - lambda(idir)*ptot_cart(idir)
   end do

!  Update cell parameters
   if (relaxstr == 1) then
     delta_eta(:) = zero
     do istrain = 1, 6
       if (irelaxstrain(istrain) == 1) then
         posi = posi + 1
         delta_eta(istrain) = delta(posi)
       end if
     end do

     do istrain = 1, 3
       strainmat(istrain,istrain) = delta_eta(istrain)
     end do
     strainmat(2,3) = delta_eta(4)/2._dp ; strainmat(3,2) = delta_eta(4)/2._dp
     strainmat(1,3) = delta_eta(5)/2._dp ; strainmat(3,1) = delta_eta(5)/2._dp
     strainmat(2,1) = delta_eta(6)/2._dp ; strainmat(1,2) = delta_eta(6)/2._dp

     rprimd_new(:,:) = 0._dp
     do idir = 1, 3
       do jdir = 1, 3
         do ii = 1, 3
           rprimd_new(jdir,idir) = rprimd_new(jdir,idir) + &
&           rprimd(ii,idir)*strainmat(ii,jdir)
         end do
       end do
     end do
     rprimd_new(:,:) = rprimd_new(:,:) + rprimd(:,:)

     acell_new(:) = zero
     do idir = 1, 3
       do jdir = 1, 3
         acell_new(idir) = acell_new(idir) + &
&         rprimd_new(jdir,idir)*rprimd_new(jdir,idir)
       end do
       acell_new(idir) = sqrt(acell_new(idir))
       rprim(:,idir) = rprimd_new(:,idir)/acell_new(idir)
     end do

   end if          ! relaxstr == 1

!  Write out the results

   if (iwrite) then
     write(iout,*)
     write(iout,'(a,80a,a)') ch10,('=',ii=1,80),ch10
     write(iout,*)
     write(iout,*)'Relaxation of the geometry at fixed polarization:'
     write(iout,*)
     write(iout,'(a,3(2x,f16.9))')' Lambda = ',(lambda(idir),idir = 1, 3)
     write(iout,'(a,e16.9)')' Value of the energy functional E_1 = ',e1
     write(iout,*)
     write(iout,*)'Difference between actual value of the Polarization (C/m^2)'
     write(iout,*)'and the target value:'
   end if
   diffpol(:) = (ptot_cart(:) - targetpol(:))*e_Cb/((Bohr_Ang*1.0d-10)**2)
   if (iwrite) write(iout,'(3(3x,f16.9))')(diffpol(idir),idir = 1, 3)

   if (relaxat == 1) then
!    Compute the forces induced on the atoms by the electric field
!    The strength of the field is determined by lambda
     felfd(:,:) = zero
     do iatom = 1, natom
       do idir = 1, 3
         do jdir = 1, 3
           index = idir +3*((iatom - 1) + mpert*((jdir - 1) + 3*(natom + 1)))
           felfd(idir,iatom) = felfd(idir,iatom) - lambda(jdir)*blkval(1,index)/ucvol
         end do
       end do
     end do

!    Compute remaining forces and write them out

     fdiff(:,:) = fcart(:,:) - felfd(:,:)
     if (iwrite) then
       write(iout,*)
       write(iout,*)'Difference between the Hellmann-Feynman forces'
       write(iout,*)'and the forces induced by the electric field'
       write(iout,*)'(cartesian coordinates, hartree/bohr)'
     end if
     fmax = zero
     do iatom = 1, natom
       if (iwrite) write(iout,'(3(3x,es16.9))')(fdiff(idir,iatom),idir = 1, 3)
       do idir = 1, 3
         if (abs(fdiff(idir,iatom)) > fmax) fmax = abs(fdiff(idir,iatom))
       end do
     end do

     if (iwrite) then
       write(iout,'(a,3x,es16.9)')' fmax = ',fmax
       write(iout,*)
       write(iout,*)'Change of cartesian coordinates (delta_xcart):'
       do iatom = 1, natom
         write(iout,'(5x,i3,3(2x,f16.9))')iatom,(delta_xcart(idir,iatom),idir = 1, 3)
       end do
       write(iout,*)
       write(iout,*)'New cartesian coordinates (xcart_new):'
       write(iout,*)'  xcart'
       do iatom = 1, natom
         write(iout,'(3(3x,d22.14))')(xcart_new(idir,iatom),idir = 1, 3)
       end do
       write(iout,*)
       write(iout,*)'New reduced coordinates (xred_new):'
       write(iout,*)'  xred'
       do iatom = 1, natom
         write(iout,'(3(3x,d22.14))')(xred_new(idir,iatom),idir = 1, 3)
       end do
     end if

   end if         ! relaxat == 1

   if (relaxstr == 1) then

!    Compute the stresses induced by the electric field
     sigelfd(:) = zero
     istrain = 0
     do ipert = 1, 2
       jpert = natom + 2 + ipert
       do idir = 1, 3
         istrain = istrain + 1
         do jdir = 1, 3
           index = idir +3*((jpert - 1) + mpert*((jdir - 1) + 3*(natom + 1)))
           sigelfd(istrain) = sigelfd(istrain) + lambda(jdir)*blkval(1,index)
         end do
         sigelfd(istrain) = sigelfd(istrain)/ucvol
       end do
     end do

!    Compute the remaining stresses and write them out
     diffsig(:) = strten(:) - sigelfd(:)
     sigmax = zero
     do istrain = 1, 6
       if (abs(diffsig(istrain)) > sigmax) sigmax = abs(diffsig(istrain))
     end do
     if (iwrite) then
       write(iout,*)
       write(iout,*)'Difference between the Hellmann-Feynman stresses'
       write(iout,*)'and the stresses induced by the electric field'
       write(iout,*)'(cartesian coordinates, hartree/bohr^3)'
       write(iout,'(2x,a,f16.9,5x,a,f16.9)')'diffsig(1) = ',diffsig(1),'diffsig(4) = ',diffsig(4)
       write(iout,'(2x,a,f16.9,5x,a,f16.9)')'diffsig(2) = ',diffsig(2),'diffsig(5) = ',diffsig(5)
       write(iout,'(2x,a,f16.9,5x,a,f16.9)')'diffsig(3) = ',diffsig(3),'diffsig(6) = ',diffsig(6)
       write(iout,'(a,3x,es16.9)')' sigmax = ',sigmax
       write(iout,*)
       write(iout,*)'Induced strain (delta_eta):'
       write(iout,'(2x,a,f16.9,5x,a,f16.9)')'delta_eta(1) = ',delta_eta(1),'delta_eta(4) = ',delta_eta(4)
       write(iout,'(2x,a,f16.9,5x,a,f16.9)')'delta_eta(2) = ',delta_eta(2),'delta_eta(5) = ',delta_eta(5)
       write(iout,'(2x,a,f16.9,5x,a,f16.9)')'delta_eta(3) = ',delta_eta(3),'delta_eta(6) = ',delta_eta(6)
       write(iout,*)
       write(iout,*)'New lattice constants (acell_new):'
       write(iout,*)'  acell'
       write(iout,'(3(2x,d22.14))')(acell_new(idir),idir = 1, 3)
       write(iout,*)
       write(iout,*)'New primitive vectors (rprim_new):'
       write(iout,*)'  rprim'
       write(iout,'(3(2x,d22.14))')(rprim(idir,1),idir = 1, 3)
       write(iout,'(3(2x,d22.14))')(rprim(idir,2),idir = 1, 3)
       write(iout,'(3(2x,d22.14))')(rprim(idir,3),idir = 1, 3)
     end if
   end if         ! relaxstr /= 0

 end if    !  (relaxat /= 0).or.(relaxstr /= 0)

 ABI_FREE(cfac)
 ABI_FREE(fdiff)
 ABI_FREE(felfd)
 ABI_FREE(delta)
 ABI_FREE(fcart)
 ABI_FREE(fcmat)
 ABI_FREE(ifcmat)
 ABI_FREE(ipvt)
 ABI_FREE(rfpert)
 ABI_FREE(vec)
 ABI_FREE(zgwork)
 ABI_FREE(irelaxat)

end subroutine relaxpol
!!***

end module m_relaxpol
!!***
