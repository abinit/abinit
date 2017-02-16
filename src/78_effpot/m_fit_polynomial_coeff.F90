!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_fit_polynomial_coeff
!!
!! NAME
!! m_fit_polynomial_coeff
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2010-2015 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_fit_polynomial_coeff

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_polynomial_coeff
 use m_atomdata
 use m_xmpi
 use m_sort
 use m_phonon_supercell
 use m_effective_potential, only :  effective_potential_type
 use m_io_tools,   only : open_file
 use m_abihist, only : abihist

 implicit none

 public :: fit_polynomial_coeff_free
 public :: fit_polynomial_coeff_init
 public :: fit_polynomial_coeff_getList
 public :: fit_polynomial_coeff_get
 public :: fit_polynomial_getOrder1
 public :: fit_polynomial_getOrder2
 public :: fit_polynomial_printSystemFiles
!!***

!!****t* m_fit_polynomial_coeff/fit_polynomial_coeff_type
!! NAME
!! fit_polynomial_coeff_type
!!
!! FUNCTION
!! structure for a polynomial coefficient
!! contains the value of the coefficient and a 
!! list of terms (displacement) relating to the coefficient
!!
!! SOURCE

 type, public :: fit_polynomial_coeff_type

   character(200) :: name
!     Name of the polynomial_coeff (Sr_y-O1_y)^3) for example

   type(polynomial_coeff_type),dimension(:),allocatable :: coefficients
!     terms(nterm)
!     contains all the displacements for this coefficient      

 end type fit_polynomial_coeff_type
!!***

CONTAINS  !===========================================================================================


!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_init
!!
!! NAME
!! fit_polynomial_coeff_init
!!
!! FUNCTION
!! Initialize scell structure, from unit cell vectors, qpoint chosen, and atoms
!!
!! INPUTS
!!  name     = Name of the polynomial_coeff (Sr_y-O1_y)^3) for example
!!  nterm   = Number of terms (short range interaction) for this polynomial_coeff
!!  coefficient  = Value of the coefficient of this term
!!  termstype(nterm) = Polynomial_term_type contains all the displacements for this coefficient 
!!
!! OUTPUT
!! polynomial_coeff = polynomial_coeff structure to be initialized
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_polynomial_coeff_init(coefficient,name,nterm,polynomial_coeff,terms)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nterm
 real(dp),intent(in) :: coefficient
!arrays
 character(200),intent(in) :: name
 type(polynomial_term_type),intent(in) :: terms(nterm)
 type(fit_polynomial_coeff_type), intent(out) :: polynomial_coeff
!Local variables-------------------------------
!scalar
!arrays

! *************************************************************************
 
!First free before initilisation
 call fit_polynomial_coeff_free(polynomial_coeff)

!Initilisation

end subroutine fit_polynomial_coeff_init
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_free
!!
!! NAME
!! fit_polynomial_coeff_free
!!
!! FUNCTION
!! Free polynomial_coeff
!!
!! INPUTS
!!
!! OUTPUT
!! polynomial_coeff = polynomial_coeff structure to be free
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_polynomial_coeff_free(polynomial_coeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(fit_polynomial_coeff_type), intent(inout) :: polynomial_coeff
!Local variables-------------------------------
!scalar
!arrays

! *************************************************************************

end subroutine fit_polynomial_coeff_free
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_get
!!
!! NAME
!! fit_polynomial_coeff_get
!!
!! FUNCTION
!! Get the list of the coefficient of the polynome
!!
!! INPUTS
!!
!! OUTPUT
!! polynomial_coeff = polynomial_coeff structure to be free
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_polynomial_coeff_get(cut_off,eff_pot,option)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_get'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option
 real(dp),intent(in):: cut_off
!arrays
 type(effective_potential_type), intent(in) :: eff_pot
!Local variables-------------------------------
!scalar
 integer :: ia,ib,icoeff,ii,irpt,irpt_ref,itypat
 integer :: jj,lim1,lim2,lim3
 integer :: natom,ncoeff1,ncoeff2,ncoeff_sym,nrpt,nsym
 integer :: r1,r2,r3
!arrays
 integer :: ncell(3)
 integer,allocatable :: cell(:,:)
 integer,allocatable :: list_symcoeff(:,:,:)
 ! integer,allocatable :: powers(:)
 ! integer,allocatable :: indsym(:,:,:) ,symrec(:,:,:)
 real(dp) :: rprimd(3,3)
 real(dp),allocatable :: dist(:,:,:),rpt(:,:)
 real(dp),allocatable :: xcart(:,:),xred(:,:)
 character(len=1) :: powerchar
 character(len=5),allocatable :: symbols(:)
 ! character(len=500) :: message
  type(polynomial_coeff_type),dimension(:),allocatable :: coeffs1,coeffs2
 ! type(polynomial_coeff_type),dimension(:),allocatable :: coeffs2,coeffs3,coeffs4
 ! type(polynomial_term_type),dimension(:),allocatable :: terms

! *************************************************************************

!Initialisation of variables
 natom  = eff_pot%crystal%natom
 nsym   = eff_pot%crystal%nsym
 rprimd = eff_pot%crystal%rprimd

 ABI_ALLOCATE(xcart,(3,natom))
 ABI_ALLOCATE(xred,(3,natom))
 xcart(:,:) = eff_pot%crystal%xcart(:,:)
 xred(:,:)  = eff_pot%crystal%xred(:,:)

!Set the size of the interaction
 ncell = (/anint(cut_off/rprimd(1,1))+1,&
&          anint(cut_off/rprimd(2,2))+1,&
&          anint(cut_off/rprimd(3,3))+1/)

 lim1=((ncell(1)/2))
 lim2=((ncell(2)/2))
 lim3=((ncell(3)/2))
 if(mod(ncell(1),2)/=0) lim1=lim1+1
 if(mod(ncell(2),2)/=0) lim2=lim2+1
 if(mod(ncell(3),2)/=0) lim3=lim3+1
 nrpt=(2*lim1+1)*(2*lim2+1)*(2*lim3+1)

!compute new ncell
 ncell(1) = 2*lim1+1
 ncell(2) = 2*lim2+1
 ncell(3) = 2*lim3+1

!Build the rpt point 
 ABI_ALLOCATE(rpt,(3,nrpt))
 ABI_ALLOCATE(cell,(3,nrpt))
 
!WARNING:
!Put the reference cell into the first element
!the code will first deal with the atoms of the first cell
 irpt = one
 irpt_ref = one 
 rpt(:,1) = zero
 cell(:,irpt)=zero
!Fill other rpt:
 do r1=lim1,-lim1,-1
   do r2=lim2,-lim2,-1
     do r3=lim3,-lim3,-1
       if(r1==0.and.r2==0.and.r3==0) then
         cycle
       end if
       irpt=irpt+1
       rpt(1,irpt)=r1*rprimd(1,1)+r2*rprimd(1,2)+r3*rprimd(1,3)
       rpt(2,irpt)=r1*rprimd(2,1)+r2*rprimd(2,2)+r3*rprimd(2,3)
       rpt(3,irpt)=r1*rprimd(3,1)+r2*rprimd(3,2)+r3*rprimd(3,3)
       cell(1,irpt)=r1;cell(2,irpt)=r2;cell(3,irpt)=r3
     end do
   end do
 end do

 ABI_ALLOCATE(symbols,(natom))
 do ia=1,natom
   symbols(ia) = adjustl(znucl2symbol(eff_pot%crystal%znucl(eff_pot%crystal%typat(ia))))
 end do

!Check the atoms of the same type and add numerotation
 itypat = zero
 do itypat =1,eff_pot%crystal%ntypat
   ii = zero
   do ia=1,natom
     if(eff_pot%crystal%typat(ia)==itypat) then
       ii = ii + 1
     end if
   end do
   if(ii>1)then
     jj=1
     do ia=1,natom
       if(eff_pot%crystal%typat(ia)==itypat) then
         write(powerchar,'(I0)') jj
         symbols(ia) = trim(symbols(ia))//trim(powerchar)
         jj=jj+1
       end if
     end do
   end if
 end do

!Compute the distances between atoms
!Now dist(ia,ib,irpt) contains the distance from atom ia to atom ib in unit cell irpt.
 ABI_ALLOCATE(dist,(natom,natom,nrpt))
 do ia=1,natom
   do ib=1,natom
     do irpt=1,nrpt
       dist(ia,ib,irpt) = ((xcart(1,ib)-xcart(1,ia)+rpt(1,irpt))**2+&
&                          (xcart(2,ib)-xcart(2,ia)+rpt(2,irpt))**2+&
&                          (xcart(3,ib)-xcart(3,ia)+rpt(3,irpt))**2)**0.5
     end do
   end do
 end do


 call fit_polynomial_coeff_getList(cell,cut_off,dist,eff_pot,list_symcoeff,&
&                                  natom,ncoeff_sym,nrpt)

 call fit_polynomial_getOrder1(cell,coeffs1,cut_off,list_symcoeff,natom,ncoeff1,ncoeff_sym,&
&                              nrpt,nsym,rprimd,symbols,xcart)

 call fit_polynomial_getOrder2(cell,coeffs2,cut_off,list_symcoeff,natom,ncoeff2,ncoeff_sym,&
&                              nrpt,nsym,rprimd,symbols,xcart)

 do icoeff=1,ncoeff1
   call polynomial_coeff_free(coeffs1(icoeff))
 end do
 if(allocated(coeffs1)) then
   ABI_DEALLOCATE(coeffs1)
 end if

 do icoeff=1,ncoeff2
   call polynomial_coeff_free(coeffs2(icoeff))
 end do
 if(allocated(coeffs2)) then
   ABI_DEALLOCATE(coeffs2)
 end if

 ABI_DEALLOCATE(cell)
 ABI_DEALLOCATE(dist)
 ABI_DEALLOCATE(list_symcoeff)
 ABI_DEALLOCATE(rpt)
 ABI_DEALLOCATE(symbols)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xred)

end subroutine fit_polynomial_coeff_get
!!***


!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_getList
!!
!! NAME
!! fit_polynomial_coeff_getList
!!
!! FUNCTION
!! Get the list of the coefficient of the polynome
!!
!! INPUTS
!!
!! OUTPUT
!! polynomial_coeff = polynomial_coeff structure to be free
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_polynomial_coeff_getList(cell,cut_off,dist,eff_pot,list_symcoeff,&
&                                       natom,ncoeff_sym,nrpt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_getList'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrpt
 integer,intent(out) :: ncoeff_sym
 real(dp),intent(in):: cut_off
!arrays
 integer,intent(in) :: cell(3,nrpt)
 real(dp),intent(in):: dist(natom,natom,nrpt)
 type(effective_potential_type), intent(in) :: eff_pot
 integer,allocatable,intent(out) :: list_symcoeff(:,:,:)
!Local variables-------------------------------
!scalar
 integer :: ia,ib,icoeff,icoeff2,icoeff_tot,icoeff_tmp,idisy1,idisy2,ii
 integer :: ipesy1,ipesy2,isym,irpt,irpt3,irpt_ref,irpt_sym
 integer :: jj,mu
 integer :: ncoeff,ncoeff_max,nu
 integer :: nsym,shift_atm1(3)
 integer :: shift_atm2(3)
 real(dp):: tolsym8
 logical :: found
!arrays
 integer :: sym(3,3)
 integer :: transl(3)
 integer,allocatable :: list(:),list_symcoeff_tmp(:,:,:)
 integer,allocatable :: indsym(:,:,:) ,symrec(:,:,:)
 real(dp),allocatable :: blkval(:,:,:,:,:),tnons(:,:)
 real(dp),allocatable :: wkdist(:),xcart(:,:),xred(:,:)
 real(dp) :: difmin(3)
 real(dp) :: rprimd(3,3)
 real(dp) :: tratom(3)
 character(len=500) :: message

! *************************************************************************


!Initialisation of variables
 nsym   = eff_pot%crystal%nsym
 rprimd = eff_pot%crystal%rprimd
 ABI_ALLOCATE(xcart,(3,natom))
 ABI_ALLOCATE(xred,(3,natom))
 xcart(:,:) = eff_pot%crystal%xcart(:,:)
 xred(:,:)  = eff_pot%crystal%xred(:,:)
 ncoeff_max = nrpt*natom*natom*3*3

!Found the ref cell
 irpt_ref = one 
 do irpt=1,nrpt
!   write(100,*)"tutu",irpt,'=>',cell(:,irpt)
   if(all(cell(:,irpt)==0))then
     irpt_ref = irpt
     exit
   end if
 end do

!Obtain a list of rotated atom labels:
 ABI_ALLOCATE(indsym,(4,nsym,natom))
 ABI_ALLOCATE(symrec,(3,3,nsym))
 ABI_ALLOCATE(tnons,(3,nsym))
 symrec = eff_pot%crystal%symrec
 tnons  = eff_pot%crystal%tnons

 tolsym8=tol8 
 call symatm(indsym,natom,nsym,symrec,tnons,&
&            tolsym8,eff_pot%crystal%typat,eff_pot%crystal%xred)


 ABI_ALLOCATE(blkval,(3,natom,3,natom,nrpt))
 ABI_ALLOCATE(list,(natom*nrpt))
 ABI_ALLOCATE(list_symcoeff_tmp,(5,ncoeff_max,nsym))
 ABI_ALLOCATE(wkdist,(natom*nrpt))

 write(message,'(a,F6.3,a)') " Cut-off of ",cut_off," Angstrom is imposed"
 call wrtout(std_out,message,'COLL') 

 
!Set to one blkval, all the coeff have to be compute
 blkval = one 
 icoeff = one
 icoeff_tot = one
 list_symcoeff_tmp = zero

!Big loop over generic atom 
 do ia=1,natom
   wkdist(:)=reshape(dist(ia,:,:),(/natom*nrpt/))
   do ii=1,natom*nrpt
     list(ii)=ii
   end do
   call sort_dp(natom*nrpt,wkdist,list,tol14)
   
   do ii=1,natom*nrpt
!    Get the irpt and ib
     irpt=(list(ii)-1)/natom+1     
     ib=list(ii)-natom*(irpt-1)

     if(dist(ia,ib,irpt) > cut_off ) then
!      If this distance is superior to the cutoff, we don't compute
       blkval(:,ia,:,ib,irpt)= zero
       if(irpt==irpt_ref)blkval(:,ib,:,ia,irpt)= zero 
!      Stop the loop
       exit
     end if
     do mu=1,3
       do nu=1,3
!      Check if : - The coefficient is not yet compute
!                 - The directions are the same
!                 - The atoms are not equivalent
         if (mu/=nu) then 
           blkval(mu,ia,nu,ib,irpt)=zero
           blkval(nu,ia,mu,ib,irpt)=zero
           cycle
         end if
!        Pass if the atoms are identical and in the ref cell
         if(irpt==irpt_ref.and.ia==ib) then
           blkval(mu,ia,nu,ib,irpt)=zero
           blkval(nu,ib,mu,ia,irpt)=zero
           cycle
         end if
         
         if(blkval(mu,ia,nu,ib,irpt)==1)then
!          Loop over symmetries 
           do isym=1,nsym
!            Get the symmetry matrix 
             sym(:,:) = eff_pot%crystal%symrel(:,:,isym)

!            Get the corresponding atom and shift with the symetries 
!            For atom 1
             ipesy1 = indsym(4,isym,ia)
             shift_atm1 = indsym(1:3,isym,ia)
!            And atom 2
             do jj=1,3 ! Apply transformation to original coordinates.
               tratom(jj) = dble(sym(1,jj))*(xred(1,ib)+cell(1,irpt)-tnons(1,isym))&
&                          +dble(sym(2,jj))*(xred(2,ib)+cell(2,irpt)-tnons(2,isym))&
&                          +dble(sym(3,jj))*(xred(3,ib)+cell(3,irpt)-tnons(3,isym))
             end do

!            Find symmetrically equivalent atom
             call symchk(difmin,ipesy2,natom,tratom,transl,eff_pot%crystal%typat(ib),&
&                        eff_pot%crystal%typat,xred(:,:))

!            Put information into array indsym: translations and label
             shift_atm2(:)= transl(:) - shift_atm1(:)

             found = .false.
             do irpt3=1,nrpt
               if(cell(1,irpt3)==shift_atm2(1).and.&
&                 cell(2,irpt3)==shift_atm2(2).and.&
&                 cell(3,irpt3)==shift_atm2(3))then
                 found = .true.
                 irpt_sym = irpt3
               end if
             end do
             
!            Now that a symmetric perturbation has been obtained,
!            including the expression of the symmetry matrix, see
!            if the symmetric perturbations are available
             do idisy1=1,3
               do idisy2=1,3              

                 if (idisy1/=idisy2) then 
!                  Remove this term (is not computed)
!                  Also remove opposite term... (Srx-Tix) = (Ti-Srx)
                   blkval(idisy1,ipesy1,idisy2,ipesy2,irpt_sym) = 0
                   blkval(idisy2,ipesy1,idisy1,ipesy2,irpt_sym) = 0
                   cycle
                 else
                   if(sym(mu,idisy1)/=0.and.sym(nu,idisy2)/=0)then
                     if(.not.found.or.(irpt_sym==irpt_ref.and.ipesy1==ipesy2)) then
!                      Remove this term (is not computed)
!                      Also remove oposite term... (Srx-Tix) = (Ti-Srx)
                       blkval(idisy1,ipesy1,idisy2,ipesy2,irpt_sym) = 0
                       blkval(idisy2,ipesy2,idisy1,ipesy1,irpt_sym) = 0
                       cycle
                     else
!                      Fill the list with the coeff and symetric (need all symetrics)
                       list_symcoeff_tmp(1:4,icoeff,isym)=(/idisy1,ipesy1,ipesy2,irpt_sym/)
                       list_symcoeff_tmp(5,icoeff,isym)= sym(mu,idisy1)
                     end if
                   end if
                 end if
               end do
             end do
           end do ! end loop sym
           icoeff = icoeff + 1 
         end if
!        This coeff is now computed 
         blkval(mu,ia,nu,ib,irpt)= zero 
       end do ! end loop nu
     end do ! end loop mu
   end do ! end loop ii
 end do ! end loop ia

!Reset the output
 if(allocated(list_symcoeff))then
   ABI_DEALLOCATE(list_symcoeff)
 end if

!Transfert the final array with all the coefficients
!With this array, we can access to all the terms presents
!ncoeff1 + symetrics
!first dimension is 4 (mu,ia,ib,irpt)
!irpt is the index of the cell of the atom ib in the cell array
!example cell(:,irpt=12) can be (-1 0 -2). The cell of ia is 
!always 0 0 0
!Transfert the final array for the list of irreductible coeff and symetries
!With this array, we can access to the irretuctible coefficients (ncoeff1) and  
!all the symetrics of these coefficients (nsym)
!first dimension is 5 (mu,ia,ib,irpt,icoeff)
!icoeff is the position of this coefficients in the list_fullcoeff array

!Count the number of coefficient total and symetrics
 ncoeff = zero
 do icoeff_tmp = 1,ncoeff_max
   if(.not.(all(list_symcoeff_tmp(:,icoeff_tmp,1)==zero)))then
!    Don't take into acount if it's opposite
     ncoeff = ncoeff + 1
   end if
 end do

 ABI_ALLOCATE(list_symcoeff,(6,ncoeff,nsym))

 list_symcoeff = zero
 icoeff = zero
 do icoeff_tmp = 1,ncoeff_max
   if(.not.(all(list_symcoeff_tmp(:,icoeff_tmp,1)==zero)))then
!    Don't take into acount if it's opposite
     icoeff = icoeff + 1
     list_symcoeff(1:5,icoeff,:) = list_symcoeff_tmp(1:5,icoeff_tmp,:)
   end if
 end do

!set the dimension six of list_symcoeff(6,icoeffs,1)
 ncoeff_max = zero
 do icoeff = 1,ncoeff
!  found the index of each coeff in list_fullcoeff
   do isym = 1,nsym
     icoeff2 = getCoeffFromList(list_symcoeff(:,:,1),&
&                               list_symcoeff(2,icoeff,isym),&
&                               list_symcoeff(3,icoeff,isym),&
&                               list_symcoeff(4,icoeff,isym),&
&                               list_symcoeff(1,icoeff,isym),&
&                               real(list_symcoeff(5,icoeff,isym),dp),ncoeff)
     list_symcoeff(6,icoeff,isym) = icoeff2

!    Check if the opposite is not already compute
!    In this case, we set the 6 dimention to the opposite
     icoeff_tmp = getCoeffFromList(list_symcoeff(:,:,1),&
&                                  list_symcoeff(3,icoeff,isym),&
&                                  list_symcoeff(2,icoeff,isym),&
&                                  list_symcoeff(4,icoeff,isym),&
&                                  list_symcoeff(1,icoeff,isym),&
&                                  real(list_symcoeff(5,icoeff,isym),dp),ncoeff)
     if(icoeff_tmp/=zero.and.icoeff_tmp<icoeff2) list_symcoeff(6,icoeff,isym) = icoeff_tmp
     if(icoeff2 > ncoeff_max) ncoeff_max = icoeff2
!     write(100,*),"lala",icoeff,isym,"=>",list_symcoeff(:,icoeff,isym)
   end do
 end do

! close(100)


!Set the max number of coeff inside list_symcoeff
 ncoeff_sym = ncoeff_max

!do somes checks
 do icoeff = 1,ncoeff_sym 
   do isym = 1,nsym
     if(list_symcoeff(6,icoeff,isym)==0)then
       write(message, '(a,i0,a,I0,4a)' )&
&           'The coefficient number ',icoeff,' with the symetrie ',isym,ch10,&
&           'have no equivalent',ch10,&
&           'Action: Contact abinit group'
       MSG_BUG(message)
     end if
     if(icoeff /= list_symcoeff(6,icoeff,isym))then
       if(list_symcoeff(1,icoeff,isym)/=list_symcoeff(1,list_symcoeff(6,icoeff,isym),1))then
         write(message, '(a,i0,a,I0,2a,I0,4a)' )&
&          'The coefficient number ',icoeff,' with the symetrie ',isym,ch10,&
&          'does not refer to the same coefficient ',list_symcoeff(6,icoeff,1),ch10,&
&          'because the direction is different:',ch10,&
&          'Action: Contact abinit group'
         MSG_BUG(message)
       end if
       if(list_symcoeff(4,icoeff,isym)/=list_symcoeff(4,list_symcoeff(6,icoeff,isym),1))then
         write(message, '(a,i0,a,I0,2a,I0,4a)' )&
&          'The coefficient number ',icoeff,' with the symetrie ',isym,ch10,&
&          'does not refer to the same coefficient ',list_symcoeff(6,icoeff,1),ch10,&
&          'because the cell is different',ch10,&
&          'Action: Contact abinit group'
         MSG_BUG(message)
       end if
       if((list_symcoeff(2,icoeff,isym)/=list_symcoeff(2,list_symcoeff(6,icoeff,isym),1).and.&
&           list_symcoeff(3,icoeff,isym)/=list_symcoeff(3,list_symcoeff(6,icoeff,isym),1)).and.&
&          (list_symcoeff(2,icoeff,isym)/=list_symcoeff(3,list_symcoeff(6,icoeff,isym),1).and.&
&           list_symcoeff(3,icoeff,isym)/=list_symcoeff(2,list_symcoeff(6,icoeff,isym),1)))then
         write(message, '(a,i0,a,I0,2a,I0,4a)' )&
&          'The coefficient number ',icoeff,' with the symetrie ',isym,ch10,&
&          'does not refer to the same coefficient ',list_symcoeff(6,icoeff,1),ch10,&
&          'because the atoms different',ch10,&
&          'Action: Contact abinit group'
         MSG_BUG(message)
       end if
     end if
   end do
 end do

!Deallocation
 ABI_DEALLOCATE(blkval)
 ABI_DEALLOCATE(list)
 ABI_DEALLOCATE(list_symcoeff_tmp)
 ABI_DEALLOCATE(indsym) 
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(tnons)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xred )
 ABI_DEALLOCATE(wkdist)

end subroutine fit_polynomial_coeff_getList
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_getOrder1
!!
!! NAME
!! fit_polynomial_getOrder1
!!
!! FUNCTION
!! Free polynomial_coeff
!!
!! INPUTS
!!
!! OUTPUT
!! polynomial_coeff = polynomial_coeff structure to be free
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_polynomial_getOrder1(cell,coeffs_out,cut_off,list_symcoeff,&
&                                   natom,ncoeff_out,ncoeff,nrpt,nsym,&
&                                   rprimd,symbols,xcart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_getOrder1'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: natom,ncoeff,nsym,nrpt
 integer,intent(out) :: ncoeff_out
 real(dp),intent(in) :: cut_off
!arrays
 integer,intent(in) :: cell(3,nrpt)
 integer,intent(in) :: list_symcoeff(6,ncoeff,nsym)
 real(dp),intent(in) :: xcart(3,natom),rprimd(3,3)
 character(len=5),intent(in) :: symbols(natom)
 type(polynomial_coeff_type),allocatable,intent(inout) :: coeffs_out(:)
!Local variables-------------------------------
!scalar
 integer :: ia,ib,icoeff,icoeff1_opp,icoeff_tmp,irpt,irpt_ref
 integer :: isym,iterm,mu,ncoeff_max,ndisp,nterm_max
 real(dp):: coefficient,weight
!arrays
 integer,allocatable :: atindx(:,:),cells(:,:,:),dir_int(:)
 integer,allocatable :: blkval(:),powers(:)
 character(len=1) :: dir_char(3)
 character(len=1) :: powerchar
 character(len=1) :: mutodir(3) = (/"x","y","z"/)
 character(len=100):: name,text
 character(len=500) :: message
 character(len=fnlen) :: filename
 type(polynomial_term_type),dimension(:),allocatable :: terms
 type(polynomial_coeff_type),allocatable :: coeffs_tmp(:)

! *************************************************************************

!Initialisation of variables
 nterm_max  = ncoeff
 ncoeff_max = ncoeff
 ndisp = 1
 ABI_ALLOCATE(blkval,(ncoeff))
 ABI_ALLOCATE(coeffs_tmp,(ncoeff_max))
 ABI_ALLOCATE(terms,(nterm_max))
 blkval(:) = one
 icoeff_tmp = zero 
 icoeff1_opp = zero
 ABI_ALLOCATE(atindx,(2,ndisp))
 ABI_ALLOCATE(cells,(3,2,ndisp))
 ABI_ALLOCATE(dir_int,(ndisp))
 ABI_ALLOCATE(powers,(ndisp))

!Found the ref cell
 irpt_ref = one 
 do irpt=1,nrpt
   if(all(cell(:,irpt)==0))then
     irpt_ref = irpt
     exit
   end if
 end do

 write(message,'(a)') " Irreductible coefficient and associated atom 1, atom 2 and direction:"
 call wrtout(std_out,message,'COLL') 

  do icoeff=1,ncoeff
    if(blkval(icoeff)==1)then
!     Reset counter
      iterm = zero
      coefficient = one
      do isym=1,nsym
!       Get index of this displacement term
        mu   = list_symcoeff(1,icoeff,isym)
        ia   = list_symcoeff(2,icoeff,isym)
        ib   = list_symcoeff(3,icoeff,isym)
        irpt = list_symcoeff(4,icoeff,isym)
        weight = list_symcoeff(5,icoeff,isym)
!       And fill arrays for the initialisation
        atindx(1,:) = ia; atindx(2,:) = ib; dir_char(:) = mutodir(mu);
        dir_int(:)  = mu
        ndisp  = 1 
        powers(:)   = one
        cells(:,1,1) = (/0,0,0/)
        cells(:,2,1) = cell(:,irpt)
        if(blkval(list_symcoeff(6,icoeff,isym))==1)then
          iterm = iterm + 1
          call polynomial_term_init(atindx,cells,dir_int,ndisp,terms(iterm),powers,weight)
        end if
      end do!end do sym

      if(iterm > 0)then
!      increase coefficients and set it
        icoeff_tmp = icoeff_tmp + 1
        call polynomial_coeff_init(coefficient,iterm,coeffs_tmp(icoeff_tmp),terms(1:iterm))
      end if

!     Deallocate the terms
      do iterm=1,nterm_max
        call polynomial_term_free(terms(iterm))
      end do

    end if!end if blkval==1

!   Set this coeff and all symetric to 0 
!   because already computed
    do isym=1,nsym
      blkval(list_symcoeff(6,icoeff,isym))=zero
      mu   = list_symcoeff(1,icoeff,isym)
      ia   = list_symcoeff(2,icoeff,isym)
      ib   = list_symcoeff(3,icoeff,isym)
      irpt = list_symcoeff(4,icoeff,isym)
!     find the opposite if the to atoms are in the same cell (ref cell)
!     Can not have the opposite of (Sr-O1[100]) because 1st atom is always 
!     in the unit cell
      if(irpt==irpt_ref) then
        icoeff1_opp = getCoeffFromList(list_symcoeff(:,:,1),ib,ia,irpt,mu,weight,ncoeff)
        blkval(icoeff1_opp) = zero
      end if
    end do
  end do!end do coeff_sym

 ABI_DEALLOCATE(terms)
 ABI_DEALLOCATE(blkval)
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(cells)
 ABI_DEALLOCATE(dir_int)
 ABI_DEALLOCATE(powers)

!Count the number of terms
 ncoeff_out = zero
 do icoeff_tmp=1,ncoeff_max
   if (coeffs_tmp(icoeff_tmp)%coefficient/=zero)then
     ncoeff_out = ncoeff_out + 1
   end if
 end do

!Transfer in the final array
 ABI_ALLOCATE(coeffs_out,(ncoeff_out))
 icoeff = zero
 do icoeff_tmp=1,ncoeff_max
   if (coeffs_tmp(icoeff_tmp)%coefficient/=zero)then
     icoeff = icoeff + 1
!    Get the name of this coefficient if the term is the first
     name = ""
     write(powerchar,'(I0)') 1
     call polynomial_coeff_getName(text,atm1=symbols(coeffs_tmp(icoeff_tmp)%terms(1)%atindx(1,1)),&
&                                       atm2=symbols(coeffs_tmp(icoeff_tmp)%terms(1)%atindx(2,1)),&
&                                       dir=mutodir(coeffs_tmp(icoeff_tmp)%terms(1)%direction(1)),&
&                                       power=trim(powerchar),&
&                                       cell_atm1=coeffs_tmp(icoeff_tmp)%terms(1)%cell(:,1,1),&
&                                       cell_atm2=coeffs_tmp(icoeff_tmp)%terms(1)%cell(:,2,1))
     name = trim(name)//trim(text)
!    Set the coefficient
     call polynomial_coeff_init(one,coeffs_tmp(icoeff_tmp)%nterm,&
&                               coeffs_out(icoeff),coeffs_tmp(icoeff_tmp)%terms,&
&                               name=name)

     write(message,'(2a)')' ',trim(name)
     call wrtout(std_out,message,'COLL') 
     
     do iterm = 1,coeffs_tmp(icoeff_tmp)%nterm
       write(message,'(a,I0,a,I0,2a)') '    Atom ',coeffs_tmp(icoeff_tmp)%terms(iterm)%atindx(1,1),&
&                       ' and atom ',coeffs_tmp(icoeff_tmp)%terms(iterm)%atindx(2,1),&
&                       ' in the direction ',mutodir(coeffs_tmp(icoeff_tmp)%terms(iterm)%direction(1))
       if(any(coeffs_tmp(icoeff_tmp)%terms(iterm)%cell(:,2,1)/=zero))then
         write(message,'(2a,I0,a,I0,a,I0,a)') trim(message),' in the cell ',&
&                                       coeffs_tmp(icoeff_tmp)%terms(iterm)%cell(1,2,1),' ',&
&                                       coeffs_tmp(icoeff_tmp)%terms(iterm)%cell(2,2,1),' ',&
&                                       coeffs_tmp(icoeff_tmp)%terms(iterm)%cell(3,2,1),'.'
       end if
       call wrtout(std_out,message,'COLL') 
     end do
   end if
 end do

 filename = "terms_1st_order.xml"
 call polynomial_coeff_writeXML(coeffs_out,ncoeff_out,filename=filename)
 write(message,'(a,I0,a)')&
&       ' with ',ncoeff_out,' for the 1st order '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL') 

!Deallocation
 do icoeff=1,ncoeff_max
   call polynomial_coeff_free(coeffs_tmp(icoeff))
 end do
 ABI_DEALLOCATE(coeffs_tmp)

end subroutine fit_polynomial_getOrder1
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_getOrder2
!!
!! NAME
!! fit_polynomial_getOrder2
!!
!! FUNCTION
!! Free polynomial_coeff
!!
!! INPUTS
!!
!! OUTPUT
!! polynomial_coeff = polynomial_coeff structure to be free
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_polynomial_getOrder2(cell,coeffs_out,cut_off,list_coeff,&
&                                   natom,ncoeff_out,ncoeff,nrpt,nsym,&
&                                   rprimd,symbols,xcart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_getOrder2'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: natom,ncoeff,nsym,nrpt
 integer,intent(out) :: ncoeff_out
 real(dp),intent(in) :: cut_off
!arrays
 integer,intent(in) :: cell(3,nrpt)
 integer,intent(in) :: list_coeff(6,ncoeff,nsym)
 real(dp),intent(in) :: xcart(3,natom),rprimd(3,3)
 character(len=5),intent(in) :: symbols(natom)
 type(polynomial_coeff_type),allocatable,intent(inout) :: coeffs_out(:)
!Local variables-------------------------------
!scalar
 integer :: ia,ib,icoeff1,icoeff2,icoeff_tmp
 integer :: ii,jj,kk,idisp,irpt,irpt_ref,isym,iterm
 integer :: mu,ncoeff_max,ndisp,nterm_max,power
 real(dp):: coefficient,weight
 logical :: compatible
!arrays
 integer,allocatable :: atindx(:,:),coeffs(:),coeffs_sym(:),coeffs_opp(:)
 integer,allocatable :: cells(:,:,:),dir_int(:)
 integer,allocatable :: blkval(:,:),powers(:)
 character(len=1) :: dir_char(3)
 character(len=1) :: powerchar
 character(len=1) :: mutodir(3) = (/"x","y","z"/)
 character(len=100):: name,text
 character(len=500) :: message
 character(len=fnlen) :: filename
 type(polynomial_term_type),dimension(:),allocatable :: terms
 type(polynomial_coeff_type),allocatable :: coeffs_tmp(:)

! *************************************************************************

!Initialisation of variables
 power = 2
 nterm_max  = ncoeff
 ncoeff_max = ncoeff**power
 ndisp = power
 ABI_ALLOCATE(blkval,(ncoeff,ncoeff))
 ABI_ALLOCATE(coeffs_tmp,(ncoeff_max))
 ABI_ALLOCATE(terms,(nterm_max))
 ABI_ALLOCATE(atindx,(2,ndisp))
 ABI_ALLOCATE(cells,(3,2,ndisp))
 ABI_ALLOCATE(coeffs,(ndisp))
 ABI_ALLOCATE(coeffs_sym,(ndisp))
 ABI_ALLOCATE(coeffs_opp,(ndisp))
 ABI_ALLOCATE(dir_int,(ndisp))
 ABI_ALLOCATE(powers,(ndisp))

 blkval(:,:) = one
 icoeff_tmp = zero 

!Found the ref cell
 irpt_ref = one 
 do irpt=1,nrpt
   if(all(cell(:,irpt)==0))then
     irpt_ref = irpt
     exit
   end if
 end do
 
 do icoeff1=1,ncoeff
!TEST_AM
!   do isym=1,nsym
!     if(list_coeff(6,icoeff1,isym)==0)then
!       print*,"problem....",icoeff1,isym,":",list_coeff(:,icoeff1,isym)
!       stop
!     end if
!   end do
!TEST_AM
   do icoeff2=icoeff1,ncoeff

     if(blkval(icoeff1,icoeff2)==1)then 
!      Reset counter
       iterm = zero
       coefficient = one
!      Set the coeffs
       coeffs(:) = (/icoeff1,icoeff2/)

       do isym=1,nsym

         if(blkval(list_coeff(6,coeffs(1),isym),list_coeff(6,coeffs(2),isym))==1)then 
!          Treat this coeff
           weight = 1
           do idisp=1,ndisp
!            Get index of this displacement term
             mu   = list_coeff(1,coeffs(idisp),isym)
             ia   = list_coeff(2,coeffs(idisp),isym)
             ib   = list_coeff(3,coeffs(idisp),isym)
             irpt = list_coeff(4,coeffs(idisp),isym)
             weight = weight*list_coeff(5,coeffs(idisp),isym)
!            Fill First term arrays 
             atindx(1,idisp) = ia; atindx(2,idisp) = ib;
             dir_char(idisp) = mutodir(mu); dir_int(idisp) = mu

             powers(idisp)   = 1
             cells(:,1,idisp) = (/0,0,0/)
             cells(:,2,idisp) = cell(:,irpt)
!            Set the coefficient number of this symetric and opposite
             coeffs_sym(idisp) = list_coeff(6,coeffs(idisp),isym)
             if(irpt==irpt_ref)then
               coeffs_opp(idisp) = getCoeffFromList(list_coeff(:,:,1),ib,ia,irpt,mu,weight,ncoeff)
             else
               coeffs_opp(idisp) = list_coeff(6,coeffs(idisp),isym)
             end if
           end do
           compatible = .true.
!          Check the cut off and if the coeff is valid
           do ii=1,power
             do jj=ii+1,power
                do kk=2,3
                  if(dist(xcart(:,list_coeff(kk,coeffs_sym(ii),1)),&
&                         xcart(:,list_coeff(2,coeffs_sym(jj),1)),rprimd,&
&                          cell(:,list_coeff(4,coeffs_sym(ii),1)),&
&                          cell(:,list_coeff(4,coeffs_sym(jj),1)))>cut_off.or.&
&                    dist(xcart(:,list_coeff(kk,coeffs_sym(ii),1)),&
&                         xcart(:,list_coeff(3,coeffs_sym(jj),1)),rprimd,&
&                          cell(:,list_coeff(4,coeffs_sym(ii),1)),&
&                          cell(:,list_coeff(4,coeffs_sym(jj),1)))>cut_off)then
                    compatible =.false.
                  end if
                end do
              end do
            end do

           if(compatible)then
             iterm = iterm + 1
             call polynomial_term_init(atindx,cells,dir_int,ndisp,terms(iterm),powers,weight)
           end if

         end if!end if blkval==1
       end do!end do sym

       if(iterm > 0)then
!        Get the name of this coefficient if the term is the first         
!         write(name,'(I0,a,I0)') icoeff1," and ",icoeff2
!        increase coefficients and set it
         icoeff_tmp = icoeff_tmp + 1
         call polynomial_coeff_init(coefficient,iterm,coeffs_tmp(icoeff_tmp),terms(1:iterm))
       end if

!      Deallocate the terms
       do iterm=1,nterm_max
         call polynomial_term_free(terms(iterm))
       end do
     end if!end if blkval==1
!    Set this coeff and all symetric to 0 
!    because already computed
     do isym=1,nsym
       do idisp=1,ndisp
!        Get index of this displacement term
         mu   = list_coeff(1,coeffs(idisp),isym)
         ia   = list_coeff(2,coeffs(idisp),isym)
         ib   = list_coeff(3,coeffs(idisp),isym)
         irpt = list_coeff(4,coeffs(idisp),isym)
!        Set the coefficient number of this symetric and opposite
         coeffs_sym(idisp) = list_coeff(6,coeffs(idisp),isym)
         if(irpt==irpt_ref)then
           coeffs_opp(idisp) = getCoeffFromList(list_coeff(:,:,1),ib,ia,irpt,mu,weight,ncoeff)
         else
           coeffs_opp(idisp) = list_coeff(6,coeffs(idisp),isym)
         end if
       end do
       blkval(coeffs_sym(1),coeffs_sym(2)) = 0
       blkval(coeffs_sym(2),coeffs_sym(1)) = 0
       blkval(coeffs_opp(1),coeffs_opp(2)) = 0
       blkval(coeffs_opp(2),coeffs_opp(1)) = 0
       blkval(coeffs_sym(1),coeffs_opp(2)) = 0
       blkval(coeffs_opp(1),coeffs_sym(2)) = 0
       blkval(coeffs_sym(2),coeffs_opp(1)) = 0
       blkval(coeffs_opp(2),coeffs_sym(1)) = 0
     end do
   end do!end do coeff_sym1
 end do!end do coeff_sym2

 ABI_DEALLOCATE(terms)
 ABI_DEALLOCATE(blkval)
 ABI_DEALLOCATE(coeffs)
 ABI_DEALLOCATE(coeffs_sym)
 ABI_DEALLOCATE(coeffs_opp)
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(cells)
 ABI_DEALLOCATE(dir_int)
 ABI_DEALLOCATE(powers)
  
!Count the number of terms
 ncoeff_out = zero
 do icoeff_tmp=1,ncoeff_max
   if (coeffs_tmp(icoeff_tmp)%coefficient/=zero)then
     ncoeff_out = ncoeff_out + 1
   end if
 end do

!Transfer in the final array
 ABI_ALLOCATE(coeffs_out,(ncoeff_out))
 icoeff1 = zero
 do icoeff_tmp=1,ncoeff_max
   if (coeffs_tmp(icoeff_tmp)%coefficient/=zero)then
     name = ''
     do idisp=1,coeffs_tmp(icoeff_tmp)%terms(1)%ndisp
       write(powerchar,'(I0)') coeffs_tmp(icoeff_tmp)%terms(1)%power(idisp)
       call polynomial_coeff_getName(text,atm1=symbols(coeffs_tmp(icoeff_tmp)%terms(1)%atindx(1,idisp)),&
&                                         atm2=symbols(coeffs_tmp(icoeff_tmp)%terms(1)%atindx(2,idisp)),&
&                                         dir=mutodir(coeffs_tmp(icoeff_tmp)%terms(1)%direction(idisp)),&
&                                         power=trim(powerchar),&
&                                         cell_atm1=coeffs_tmp(icoeff_tmp)%terms(1)%cell(:,1,idisp),&
&                                         cell_atm2=coeffs_tmp(icoeff_tmp)%terms(1)%cell(:,2,idisp))
       name = trim(name)//trim(text)          
     end do
     icoeff1 = icoeff1 + 1
     call polynomial_coeff_init(one,coeffs_tmp(icoeff_tmp)%nterm,&
&                               coeffs_out(icoeff1),coeffs_tmp(icoeff_tmp)%terms,&
&                               name=name)
!TEST_AM
!     write(300,*)icoeff1,coeffs_tmp(icoeff_tmp)%nterm
!TEST_AM
   end if
 end do

 filename = "terms_2nd_order.xml"
 call polynomial_coeff_writeXML(coeffs_out,ncoeff_out,filename=filename)
 write(message,'(a,I0,a)')&
&       ' with ',ncoeff_out,' for the 2st order '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL') 

!Deallocation
 do icoeff1=1,ncoeff_max
   call polynomial_coeff_free(coeffs_tmp(icoeff1))
 end do
 ABI_DEALLOCATE(coeffs_tmp)

!TEST_AM
! close(300)
!TEST_AM
end subroutine fit_polynomial_getOrder2
!!***


!!****f* m_fit_polynomial/fit_polynomial_dist
!! NAME
!!
!! FUNCTION
!! compute the distance betwen 2 atoms in different cell
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

function dist(xcart1,xcart2,rprimd,cell1,cell2) result(distance)

!Arguments ------------------------------------
!scalar
!array

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dist'
!End of the abilint section

  real(dp),intent(in):: rprimd(3,3)
  real(dp),intent(in):: xcart1(3),xcart2(3)
  integer,intent(in) :: cell1(3),cell2(3)
  real(dp) :: distance
!Local variables -------------------------------
  real(dp) :: rpt1(3),rpt2(3)
  integer  :: mu
!! *************************************************************************
  do mu=1,3
    rpt1(mu) = cell1(1)*rprimd(mu,1)+cell1(2)*rprimd(mu,2)+cell1(3)*rprimd(mu,3)
    rpt2(mu) = cell2(1)*rprimd(mu,1)+cell2(2)*rprimd(mu,2)+cell2(3)*rprimd(mu,3)
  end do
  
  distance = ((xcart2(1)+rpt2(1)-xcart1(1)-rpt1(1))**2+&
&             (xcart2(2)+rpt2(2)-xcart1(2)-rpt1(2))**2+&
&             (xcart2(3)+rpt2(3)-xcart1(3)-rpt1(3))**2)**0.5

end function dist
!!***

!!****f* m_fit_polynomial_coeff/getCoeffFromList
!!
!! NAME
!! getCoeffFromList
!!
!! FUNCTION
!! Free polynomial_coeff
!!
!! INPUTS
!!
!! OUTPUT
!! polynomial_coeff = polynomial_coeff structure to be free
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!
!! SOURCE

function getCoeffFromList(list_coeff,ia,ib,irpt,mu,weight,ncoeff) result(coeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getCoeffFromList'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalar
 integer,intent(in) :: ia,ib,irpt,mu,ncoeff
 real(dp),intent(in):: weight
 integer :: coeff
!arrays
 integer,intent(in) :: list_coeff(6,ncoeff)
!Local variables-------------------------------
!scalar
 integer :: icoeff
!arrays

! *************************************************************************
 coeff = zero
 do icoeff = 1,ncoeff
   if(mu==list_coeff(1,icoeff).and.&
&     ia==list_coeff(2,icoeff).and.&
&     ib==list_coeff(3,icoeff).and.&
&     irpt==list_coeff(4,icoeff)) then
     coeff = icoeff
   end if
 end do

end function getCoeffFromList
!!***


!!****f* m_fit_polynomial_coeff/fit_polynomial_printSystemFiles
!!
!! NAME
!! fit_polynomial_printSystemFiles
!!
!! FUNCTION
!! Free polynomial_coeff
!!
!! INPUTS
!!
!! OUTPUT
!! polynomial_coeff = polynomial_coeff structure to be free
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_polynomial_printSystemFiles(eff_pot,hist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_printSystemFiles'
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(effective_potential_type), intent(in) :: eff_pot
 type(abihist),intent(in) :: hist
!Local variables-------------------------------
!scalar
 integer :: ia,ib,ii,jj,irpt,kk,ll,mu,nu,nstep,nshift
 integer :: unit_born=22,unit_epsiloninf=23,unit_md=24
 integer :: unit_harmonic=25,unit_ref=26,unit_strain=27,unit_sym=28
!arrays
 integer, dimension(3)  :: A,ncell
 real(dp), allocatable :: xred(:,:)
 character(len=500) :: msg
 type(supercell_type) :: supercell
! *************************************************************************

!Create new supercell corresponding to the MD
 ncell = (/2,2,2/)
 call init_supercell(eff_pot%crystal%natom, 0, real(ncell,dp), eff_pot%crystal%rprimd,&
&                    eff_pot%crystal%typat,eff_pot%crystal%xcart, supercell)

!Convert in reduced coordinates 
 ABI_ALLOCATE(xred,(3,supercell%natom_supercell))
 call xcart2xred(supercell%natom_supercell,supercell%rprimd_supercell,&
&                supercell%xcart_supercell,xred)

 A = (/ 3, 2, 1/)
 nshift = product(ncell) 

! BORN CHARGES FILE
 if (open_file('Born_Charges',msg,unit=unit_born,form="formatted",&
&    status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 do ii=1,eff_pot%crystal%ntypat
   jj = A(ii)
   do ia=1,eff_pot%crystal%natom
     if(eff_pot%crystal%typat(ia)==jj)then
       write(unit_born,'(i2,a,1F10.5)') ia,"    ",eff_pot%crystal%amu(eff_pot%crystal%typat(ia))
       do mu=1,3
         WRITE(unit_born,'(a,3(F23.14))') "     ",eff_pot%harmonics_terms%zeff(:,mu,ia)
       end do
     end if
   end do
 end do

!DIELECTRIC TENSOR FILE
 if (open_file('Dielectric_Tensor',msg,unit=unit_epsiloninf,form="formatted",&
&    status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 do mu=1,3
   WRITE(unit_epsiloninf,'(3(F23.14))') eff_pot%harmonics_terms%epsilon_inf(:,mu)
 end do


!REFERENCE STRUCTURE FILE
 if (open_file('Reference_structure',msg,unit=unit_ref,form="formatted",&
&    status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write(unit_ref,'("Energy (Hartree)")')
 write(unit_ref,'("================")')
 write(unit_ref,'(F23.14)') (hist%histE(1))
 write(unit_ref,'("")')
 write(unit_ref,'("Cell vectors")')
 write(unit_ref,'("============")')
 write(unit_ref,'(3(F23.14))') (hist%histR(:,:,1))
 write(unit_ref,'("")')
 write(unit_ref,'("Atomic positions (Bohr radius)")')
 write(unit_ref,'("==============================")')

! do ii=1,eff_pot%crystal%ntypat
!   jj = A(ii)
   do kk=1,supercell%natom_supercell
!   do kk=1,eff_pot%crystal%natom
!     if(supercell%typat_supercell(kk)==jj)then
!       do ishift=1,nshift
!         ia = eff_pot%crystal%natom*(ishift-1)+kk
         ia = kk
!        In the carlos script the atom must be in the order 11111 222222 33333 ..
!        and the order of the atom can not be change in the fit script,
!        we transform into the format of the script
         write(unit_ref,'(3(F23.14))')  hist%histXF(:,ia,1,1)
!       end do
!     end if
   end do
! end do

!Harmonic XML file
 if (open_file('harmonic.xml',msg,unit=unit_harmonic,form="formatted",&
&     status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

!Write header
 write(unit_harmonic,'("<?xml version=""1.0"" ?>")')
 write(unit_harmonic,'("<name>")')
 do irpt=1,eff_pot%harmonics_terms%ifcs%nrpt
   if(any(abs(eff_pot%harmonics_terms%ifcs%short_atmfrc(1,:,:,:,:,irpt))>tol9)) then
     write(unit_harmonic,'("  <local_force_constant units=""hartree/bohrradius**2"">")')
     write(unit_harmonic,'("    <data>")')
     do ia=1,eff_pot%crystal%natom
       do mu=1,3
         do ib=1,eff_pot%crystal%natom
           do  nu=1,3
             write(unit_harmonic,'(F22.14)', advance="no")&
&                 (eff_pot%harmonics_terms%ifcs%short_atmfrc(1,mu,ia,nu,ib,irpt))
           end do
         end do
         write(unit_harmonic,'(a)')''
       end do
     end do
     write(unit_harmonic,'("    </data>")')
     write(unit_harmonic,'("    <cell>")')
     write(unit_harmonic,'(3(I4))') (eff_pot%harmonics_terms%ifcs%cell(:,irpt))
     write(unit_harmonic,'("    </cell>")')
     write(unit_harmonic,'("  </local_force_constant>")')
   end if
 end do
 write(unit_harmonic,'("</name>")')


!STRAIN FILE
 if (open_file('Strain_Tensor',msg,unit=unit_strain,form="formatted",&
&     status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(unit_strain,'(6(F23.14))') (eff_pot%harmonics_terms%elastic_constants)

! SYM FILE
 if (open_file('symmetry_operations',msg,unit=unit_sym,form="formatted",&
&     status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(unit_sym,'("(x,y,z)  (y,-x,z) (z,x,y) (y,z,x) (x,z,y) (y,x,z) (z,y,x) (x,-y,-z) (z,-x,-y)",&
&                " (y,-z,-x) (x,-z,-y) (y,-x,-z) (z,-y,-x) (-x,y,-z) (-z,x,-y) (-y,z,-x) (-x,z,-y)",&
&                " (-y,x,-z) (-z,y,-x) (-x,-y,z) (-z,-x,y) (-y,-z,x) (-x,-z,y) (-y,-x,z) (-z,-y,x)",&
&                " (-x,-y,-z) (-z,-x,-y) (-y,-z,-x) (-x,-z,-y) (-y,-x,-z) (-z,-y,-x) (-x,y,z)",&
&                " (-z,x,y) (-y,z,x) (-x,z,y) (-y,x,z) (-z,y,x) (x,-y,z) (z,-x,y) (y,-z,x) (x,-z,y)",&
&                " (z,-y,x) (x,y,-z) (z,x,-y) (y,z,-x) (x,z,-y) (y,x,-z) (z,y,-x)")')


!MD file
 nstep = hist%mxhist
 if (open_file('Molecular_dynamic',msg,unit=unit_md,form="formatted",&
&     status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 do ii=1,nstep
   write(unit_md,'(I5)') ii-1
   write(unit_md,'(F22.14)') hist%histE(ii)
   write(unit_md,'(3(F22.14))') (hist%histR(:,:,ii))

!   do jj=1,eff_pot%crystal%ntypat
!     kk = A(jj)
     do ll=1,supercell%natom_supercell
!     do ll=1,eff_pot%crystal%natom
!       if(supercell%typat_supercell(ll)==kk)then
!       if(eff_pot%crystal%typat(ll)==kk)then
!         do ishift=1,nshift
!           ia = ll+eff_pot%crystal%natom*(ishift-1)
!           ia = kk
       ia=ll
!          In the carlos script the atom must be in the order 11111 222222 33333 ..
           write(unit_md,'(3(E22.14),3(E22.14))') hist%histXF(:,ia,1,ii),hist%histXF(:,ia,3,ii)
!         end do
!       end if
!     end do
   end do
   write(unit_md,'(6(E22.14))') hist%histS(:,ii)
 end do

!Close files
 close(unit_ref)
 close(unit_born)
 close(unit_harmonic)
 close(unit_epsiloninf)
 close(unit_md)
 close(unit_strain)
 close(unit_sym)

!Deallocation array 
 ABI_DEALLOCATE(xred)
 call destroy_supercell(supercell)

end subroutine fit_polynomial_printSystemFiles
!!***


end module m_fit_polynomial_coeff
!!***
