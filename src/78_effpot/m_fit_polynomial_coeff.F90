!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_fit_polynomial_coeff
!!
!! NAME
!! m_fit_polynomial_coeff
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2010-2017 ABINIT group (AM)
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
 use m_crystal,only : symbols_crystal
 use m_strain,only : strain_type,strain_get
 use m_effective_potential, only : effective_potential_type,effective_potential_evaluate
 use m_effective_potential, only : effective_potential_setSupercell,effective_potential_getDisp
 use m_effective_potential, only : effective_potential_GetIndexPeriodic,effective_potential_setCoeffs
 use m_io_tools,   only : open_file
 use m_abihist, only : abihist,abihist_init,abihist_free,abihist_copy

 implicit none

 public :: fit_polynomial_coeff_computeMSE
 public :: fit_polynomial_coeff_fit
 public :: fit_polynomial_coeff_free
 public :: fit_polynomial_coeff_init
 public :: fit_polynomial_coeff_getList
 public :: fit_polynomial_coeff_get
 public :: fit_polynomial_coeff_getOrder1
 public :: fit_polynomial_coeff_getOrder2
 public :: fit_polynomial_coeff_getOrder3
 public :: fit_polynomial_coeff_getOrder4
 public :: fit_polynomial_coeff_getOrder5
 public :: fit_polynomial_coeff_getFS
 public :: fit_polynomial_coeff_mapHistToRef
 public :: fit_polynomial_coeff_solve
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
!!      destroy_supercell,init_supercell,xred2xcart
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
!!      destroy_supercell,init_supercell,xred2xcart
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
!!      multibinit
!!
!! CHILDREN
!!      destroy_supercell,init_supercell,xred2xcart
!!
!! SOURCE

subroutine fit_polynomial_coeff_get(cut_off,coefficients,eff_pot,ncoeff,option)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_get'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option
 integer,intent(out):: ncoeff
 real(dp),intent(in):: cut_off
!arrays
 type(effective_potential_type), intent(inout) :: eff_pot
 type(polynomial_coeff_type),allocatable :: coefficients(:)
!Local variables-------------------------------
!scalar
 integer :: ia,ib,icoeff,icoeff2,ii,irpt,irpt_ref
 integer :: lim1,lim2,lim3
 integer :: natom,ncoeff1,ncoeff2,ncoeff3,ncoeff4,ncoeff5,ncoeff_sym,nrpt,nsym
 integer :: r1,r2,r3
!arrays
 integer :: ncell(3)
 integer,allocatable :: cell(:,:)
 integer,allocatable :: list_symcoeff(:,:,:)
 real(dp) :: rprimd(3,3)
 real(dp),allocatable :: dist(:,:,:),rpt(:,:)
 real(dp),allocatable :: xcart(:,:),xred(:,:)
 character(len=5),allocatable :: symbols(:)
 type(polynomial_coeff_type),dimension(:),allocatable :: coeffs1,coeffs2,coeffs3
 type(polynomial_coeff_type),dimension(:),allocatable :: coeffs4,coeffs5

! *************************************************************************
!Free the output
  if(allocated(coefficients))then
    do ii =1,size(coefficients)
      call polynomial_coeff_free(coefficients(ii))
    end do
    ABI_DEALLOCATE(coefficients)
  end if

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
 call symbols_crystal(eff_pot%crystal%natom,eff_pot%crystal%ntypat,eff_pot%crystal%npsp,&
&                     symbols,eff_pot%crystal%typat,eff_pot%crystal%znucl)

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

 call fit_polynomial_coeff_getOrder1(cell,coeffs1,cut_off,list_symcoeff,natom,ncoeff1,ncoeff_sym,&
&                              nrpt,nsym,rprimd,symbols,xcart)

 call fit_polynomial_coeff_getOrder2(cell,coeffs2,cut_off,list_symcoeff,natom,ncoeff2,ncoeff_sym,&
&                              nrpt,nsym,rprimd,symbols,xcart)

 call fit_polynomial_coeff_getOrder3(cell,coeffs3,cut_off,list_symcoeff,natom,ncoeff3,ncoeff_sym,&
&                              nrpt,nsym,rprimd,symbols,xcart)

 call fit_polynomial_coeff_getOrder4(cell,coeffs4,cut_off,list_symcoeff,natom,ncoeff4,ncoeff_sym,&
&                              nrpt,nsym,rprimd,symbols,xcart)

 call fit_polynomial_coeff_getOrder5(cell,coeffs5,cut_off,list_symcoeff,natom,ncoeff5,ncoeff_sym,&
&                              nrpt,nsym,rprimd,symbols,xcart)

!Final tranfert
!1- count the total number of coefficient
 ncoeff = zero
 do icoeff=1,ncoeff3
   if (coeffs3(icoeff)%coefficient /= zero) then
     ncoeff = ncoeff + 1
   end if
 end do
 do icoeff=1,ncoeff4
   if (coeffs4(icoeff)%coefficient /= zero) then
     ncoeff = ncoeff + 1
   end if
 end do
 do icoeff=1,ncoeff5
   if (coeffs5(icoeff)%coefficient /= zero) then
     ncoeff = ncoeff + 1
   end if
 end do

!2- Transfer
 ABI_ALLOCATE(coefficients,(ncoeff))
 icoeff2 = zero
 do icoeff=1,ncoeff3
   if (coeffs3(icoeff)%coefficient /= zero) then
     icoeff2 = icoeff2 + 1
     call polynomial_coeff_init(one,coeffs3(icoeff)%nterm,coefficients(icoeff2),coeffs3(icoeff)%terms,&
&                               name=coeffs3(icoeff)%name)
   end if
 end do
 do icoeff=1,ncoeff4
   if (coeffs4(icoeff)%coefficient /= zero) then
     icoeff2 = icoeff2 + 1
     call polynomial_coeff_init(one,coeffs4(icoeff)%nterm,coefficients(icoeff2),coeffs4(icoeff)%terms,&
&                               name=coeffs4(icoeff)%name)
   end if
 end do
 do icoeff=1,ncoeff5
   if (coeffs5(icoeff)%coefficient /= zero) then
     icoeff2 = icoeff2 + 1
     call polynomial_coeff_init(one,coeffs5(icoeff)%nterm,coefficients(icoeff2),coeffs5(icoeff)%terms,&
&                               name=coeffs5(icoeff)%name)
   end if
 end do

!Free them all
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

 do icoeff=1,ncoeff3
   call polynomial_coeff_free(coeffs3(icoeff))
 end do
 if(allocated(coeffs3)) then
   ABI_DEALLOCATE(coeffs3)
 end if

 do icoeff=1,ncoeff4
   call polynomial_coeff_free(coeffs4(icoeff))
 end do
 if(allocated(coeffs4)) then
   ABI_DEALLOCATE(coeffs4)
 end if

 do icoeff=1,ncoeff5
   call polynomial_coeff_free(coeffs5(icoeff))
 end do
 if(allocated(coeffs5)) then
   ABI_DEALLOCATE(coeffs5)
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
!!      destroy_supercell,init_supercell,xred2xcart
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
 integer :: ncoeff,ncoeff2,ncoeff_max,nu
 integer :: nsym,shift_atm1(3)
 integer :: shift_atm2(3)
 real(dp):: tolsym8
 logical :: found
!arrays
 integer :: sym(3,3)
 integer :: transl(3)
 integer,allocatable :: list(:),list_symcoeff_tmp(:,:,:),list_symcoeff_tmp2(:,:,:)
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
   if(all(cell(:,irpt)==0))then
     irpt_ref = irpt
!     exit
   end if
   write(100,*) irpt,":",cell(:,irpt)
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
!                      Remove this term (is not computed) Sr-Sr or not include in the cell
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

!1/ step remove the zero coeff in this array
 ncoeff = zero
 do icoeff = 1,ncoeff_max
   if(.not.(all(list_symcoeff_tmp(:,icoeff,1)==zero)))then
     ncoeff = ncoeff + 1
   end if
 end do

 ABI_ALLOCATE(list_symcoeff_tmp2,(6,ncoeff,nsym))
 list_symcoeff_tmp2 = zero
 icoeff = zero
 do icoeff_tmp = 1,ncoeff_max
   if(.not.(all(list_symcoeff_tmp(:,icoeff_tmp,1)==zero)))then
     icoeff = icoeff + 1
     list_symcoeff_tmp2(1:5,icoeff,:) = list_symcoeff_tmp(1:5,icoeff_tmp,:)
   end if
 end do


!2/ set the dimension six of list_symcoeff_tmp2(6,icoeffs,1)
!   and check is a symetric coeff is not coresspondig to an other
!   one, in this case we set this coeff to 0
! ncoeff2 = zero
 do icoeff = 1,ncoeff
!  found the index of each coeff in list_fullcoeff
   do isym = 1,nsym
     icoeff2 = getCoeffFromList(list_symcoeff_tmp2(:,:,1),&
&                               list_symcoeff_tmp2(2,icoeff,isym),&
&                               list_symcoeff_tmp2(3,icoeff,isym),&
&                               list_symcoeff_tmp2(4,icoeff,isym),&
&                               list_symcoeff_tmp2(1,icoeff,isym),&
&                               real(list_symcoeff_tmp2(5,icoeff,isym),dp),ncoeff)
     list_symcoeff_tmp2(6,icoeff,isym) = icoeff2
   end do
 end do

!2.5/do checks
 do icoeff = 1,ncoeff
   do isym = 1,nsym
     if(list_symcoeff_tmp2(6,icoeff,isym)==0)then
       write(message, '(a,i0,a,I0,4a)' )&
&           'The coefficient number ',icoeff,' with the symetrie ',isym,ch10,&
&           'have no equivalent',ch10,&
&           'Action: Contact abinit group'
       MSG_BUG(message)
     else
       if(icoeff /= list_symcoeff_tmp2(6,icoeff,isym))then
         if(list_symcoeff_tmp2(1,icoeff,isym)/=&
&           list_symcoeff_tmp2(1,list_symcoeff_tmp2(6,icoeff,isym),1))then
           write(message, '(a,i0,a,I0,2a,I0,4a)' )&
&          'The coefficient number ',icoeff,' with the symetrie ',isym,ch10,&
&          'does not refer to the same coefficient ',list_symcoeff_tmp2(6,icoeff,1),ch10,&
&          'because the direction is different:',ch10,&
&          'Action: Contact abinit group'
           MSG_BUG(message)
         end if
         if(list_symcoeff_tmp2(4,icoeff,isym)/=&
&           list_symcoeff_tmp2(4,list_symcoeff_tmp2(6,icoeff,isym),1))then
           write(message, '(a,i0,a,I0,2a,I0,4a)' )&
&          'The coefficient number ',icoeff,' with the symetrie ',isym,ch10,&
&          'does not refer to the same coefficient ',list_symcoeff_tmp2(6,icoeff,1),ch10,&
&          'because the cell is different',ch10,&
&          'Action: Contact abinit group'
           MSG_BUG(message)
         end if
         if((list_symcoeff_tmp2(2,icoeff,isym)/=&
&            list_symcoeff_tmp2(2,list_symcoeff_tmp2(6,icoeff,isym),1).and.&
&            list_symcoeff_tmp2(3,icoeff,isym)/=&
&            list_symcoeff_tmp2(3,list_symcoeff_tmp2(6,icoeff,isym),1)).and.&
&           (list_symcoeff_tmp2(2,icoeff,isym)/=&
&            list_symcoeff_tmp2(3,list_symcoeff_tmp2(6,icoeff,isym),1).and.&
&            list_symcoeff_tmp2(3,icoeff,isym)/=&
&            list_symcoeff_tmp2(2,list_symcoeff_tmp2(6,icoeff,isym),1)))then
           write(message, '(a,i0,a,I0,2a,I0,4a)' )&
&          'The coefficient number ',icoeff,' with the symetrie ',isym,ch10,&
&          'does not refer to the same coefficient ',list_symcoeff_tmp2(6,icoeff,1),ch10,&
&          'because the atoms different',ch10,&
&          'Action: Contact abinit group'
           MSG_BUG(message)
         end if
       end if
     end if
   end do
 end do

!3/ Remove useless terms like opposites
 do icoeff = 1,ncoeff
   do isym = 1,nsym
     icoeff2 = list_symcoeff_tmp2(6,icoeff,isym)
     if (icoeff2> icoeff)then
       list_symcoeff_tmp2(:,icoeff2,1) = zero
     end if
     icoeff2 = getCoeffFromList(list_symcoeff_tmp2(:,:,1),&
&                               list_symcoeff_tmp2(3,icoeff,isym),&
&                               list_symcoeff_tmp2(2,icoeff,isym),&
&                               list_symcoeff_tmp2(4,icoeff,isym),&
&                               list_symcoeff_tmp2(1,icoeff,isym),&
&                               real(list_symcoeff_tmp2(5,icoeff,isym),dp),ncoeff)
     if (icoeff2> icoeff)then
      list_symcoeff_tmp2(:,icoeff2,1) = zero
     end if
   end do
 end do

!4/ Recount the number of coeff after step 3
 ncoeff2 = zero
 do icoeff = 1,ncoeff
   if(.not.(all(list_symcoeff_tmp2(:,icoeff,1)==zero)))then
     ncoeff2 = ncoeff2 + 1
   end if
 end do

 ABI_DEALLOCATE(list_symcoeff_tmp)
 ABI_ALLOCATE(list_symcoeff_tmp,(6,ncoeff2,nsym))

 list_symcoeff_tmp = zero
 icoeff = zero

 do icoeff_tmp = 1,ncoeff
   if(.not.(all(list_symcoeff_tmp2(:,icoeff_tmp,1)==zero)))then
     icoeff = icoeff + 1
     list_symcoeff_tmp(:,icoeff,:) = list_symcoeff_tmp2(:,icoeff_tmp,:)
   end if
 end do

!5/ Final transfert
 ABI_ALLOCATE(list_symcoeff,(6,ncoeff2,nsym))
 list_symcoeff = zero
 icoeff = zero
 do icoeff = 1,ncoeff2
   list_symcoeff(1:6,icoeff,:) = list_symcoeff_tmp(1:6,icoeff,:)
   do isym=1,nsym
     write(100,*) icoeff,isym,list_symcoeff(1:6,icoeff,isym)
   end do
 end do

!6/ reset the dimension six of list_symcoeff_tmp2(6,icoeffs,1)
!   and check is a symetric coeff is not coresspondig to an other
!   one, in this case we set this coeff to 0

 do icoeff = 1,ncoeff2
!  found the index of each coeff in list_fullcoeff
   do isym = 1,nsym
     icoeff2 = getCoeffFromList(list_symcoeff(:,:,1),&
&                               list_symcoeff(2,icoeff,isym),&
&                               list_symcoeff(3,icoeff,isym),&
&                               list_symcoeff(4,icoeff,isym),&
&                               list_symcoeff(1,icoeff,isym),&
&                               real(list_symcoeff(5,icoeff,isym),dp),ncoeff)
     list_symcoeff(6,icoeff,isym) = icoeff2
   end do
 end do

!Set the max number of coeff inside list_symcoeff
 ncoeff_sym = ncoeff2

!Deallocation
 ABI_DEALLOCATE(blkval)
 ABI_DEALLOCATE(list)
 ABI_DEALLOCATE(list_symcoeff_tmp)
 ABI_DEALLOCATE(list_symcoeff_tmp2)
 ABI_DEALLOCATE(indsym) 
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(tnons)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xred )
 ABI_DEALLOCATE(wkdist)

 close(100)
end subroutine fit_polynomial_coeff_getList
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_getOrder1
!!
!! NAME
!! fit_polynomial_coeff_getOrder1
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
!!      destroy_supercell,init_supercell,xred2xcart
!!
!! SOURCE

subroutine fit_polynomial_coeff_getOrder1(cell,coeffs_out,cut_off,list_symcoeff,&
&                                   natom,ncoeff_out,ncoeff,nrpt,nsym,&
&                                   rprimd,symbols,xcart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_getOrder1'
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
 integer,allocatable :: atindx(:,:),blkval(:),cells(:,:,:),dir_int(:)
 integer,allocatable :: powers(:)
 character(len=1) :: dir_char(3)
 character(len=1) :: mutodir(9) = (/"x","y","z","1","2","3","4","5","6"/)
 character(len=100):: name
 character(len=500) :: message
 character(len=fnlen) :: filename
 type(polynomial_term_type),dimension(:),allocatable :: terms
 type(polynomial_coeff_type),allocatable :: coeffs_tmp(:)

! *************************************************************************

!Initialisation of variables
 nterm_max  = nsym
 ncoeff_max = ncoeff
 ndisp = 1

 ABI_ALLOCATE(blkval,(ncoeff_max))
 blkval = one
!  do icoeff=1,ncoeff
!    do isym=1,nsym
!      if(list_symcoeff(6,icoeff,isym) > 0.and.&
! &       list_symcoeff(6,icoeff,isym) > icoeff.and.&
! &       list_symcoeff(6,icoeff,isym) < ncoeff) blkval(list_symcoeff(6,icoeff,isym)) = zero
!    end do
!  end do
 
 ABI_ALLOCATE(coeffs_tmp,(ncoeff_max))
 ABI_ALLOCATE(terms,(nterm_max))


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
    if(blkval(icoeff)==zero)cycle
!   Reset counter
    iterm = zero
    coefficient = one
    do isym=1,nsym
!     Get index of this displacement term
      mu   = list_symcoeff(1,icoeff,isym)
      ia   = list_symcoeff(2,icoeff,isym)
      ib   = list_symcoeff(3,icoeff,isym)
      irpt = list_symcoeff(4,icoeff,isym)
      weight = list_symcoeff(5,icoeff,isym)
!     And fill arrays for the initialisation
      atindx(1,:) = ia; atindx(2,:) = ib; dir_char(:) = mutodir(mu);
      dir_int(:)  = mu
      ndisp  = 1 
      powers(:)   = one
      cells(:,1,1) = (/0,0,0/)
      cells(:,2,1) = cell(:,irpt)
      iterm = iterm + 1
      call polynomial_term_init(atindx,cells,dir_int,ndisp,terms(iterm),powers,weight,check=.true.)
    end do!end do sym

    if(iterm > 0)then
!   increase coefficients and set it
      icoeff_tmp = icoeff_tmp + 1
      call polynomial_coeff_init(coefficient,iterm,coeffs_tmp(icoeff_tmp),terms(1:iterm),check=.true.)
    end if

!   Deallocate the terms
    do iterm=1,nterm_max
      call polynomial_term_free(terms(iterm))
    end do
  end do!end do coeff_sym

 ABI_DEALLOCATE(terms)
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(cells)
 ABI_DEALLOCATE(dir_int)
 ABI_DEALLOCATE(powers)
 ABI_DEALLOCATE(blkval)

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
!    Get the name of this coefficient
     call polynomial_coeff_getName(name,natom,coeffs_tmp(icoeff_tmp),symbols)
!    Increase icoeff and fill the coeffs_out array
     icoeff = icoeff + 1
      call polynomial_coeff_init(one,coeffs_tmp(icoeff_tmp)%nterm,&
 &                               coeffs_out(icoeff),coeffs_tmp(icoeff_tmp)%terms,name=name)

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
! call polynomial_coeff_writeXML(coeffs_out,ncoeff_out,filename=filename)
 write(message,'(a,1x,I0,a)') ch10,&
&       ncoeff_out,' coefficients for the 1st order '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL') 

!Deallocation
 do icoeff=1,ncoeff_max
   call polynomial_coeff_free(coeffs_tmp(icoeff))
 end do
 ABI_DEALLOCATE(coeffs_tmp)

end subroutine fit_polynomial_coeff_getOrder1
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_getOrder2
!!
!! NAME
!! fit_polynomial_coeff_getOrder2
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
!!      destroy_supercell,init_supercell,xred2xcart
!!
!! SOURCE

subroutine fit_polynomial_coeff_getOrder2(cell,coeffs_out,cut_off,list_coeff,&
&                                   natom,ncoeff_out,ncoeff,nrpt,nsym,&
&                                   rprimd,symbols,xcart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_getOrder2'
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
 integer,allocatable :: atindx(:,:),coeffs(:)
 integer,allocatable :: cells(:,:,:),dir_int(:)
 integer,allocatable :: powers(:)
 character(len=1) :: dir_char(3)
 character(len=1) :: mutodir(9) = (/"x","y","z","1","2","3","4","5","6"/)
 character(len=100):: name
 character(len=500) :: message
 character(len=fnlen) :: filename
 type(polynomial_term_type),dimension(:),allocatable :: terms
 type(polynomial_coeff_type),allocatable :: coeffs_tmp(:)

! *************************************************************************

!Initialisation of variables
 power = 2
 nterm_max  = nsym
 ncoeff_max = ncoeff**power
 ndisp = power

 ABI_ALLOCATE(coeffs_tmp,(ncoeff_max))
 ABI_ALLOCATE(terms,(nterm_max))
 ABI_ALLOCATE(atindx,(2,ndisp))
 ABI_ALLOCATE(cells,(3,2,ndisp))
 ABI_ALLOCATE(coeffs,(ndisp))
 ABI_ALLOCATE(dir_int,(ndisp))
 ABI_ALLOCATE(powers,(ndisp))


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
   do icoeff2=icoeff1,ncoeff
     iterm = zero
     coefficient = one
!    Set the coeffs
     coeffs(:) = (/icoeff1,icoeff2/)     
     do isym=1,nsym
!      Treat this coeff
       weight = 1
       do idisp=1,ndisp
!      Get index of this displacement term
         mu   = list_coeff(1,coeffs(idisp),isym)
         ia   = list_coeff(2,coeffs(idisp),isym)
         ib   = list_coeff(3,coeffs(idisp),isym)
         irpt = list_coeff(4,coeffs(idisp),isym)
         weight = weight*list_coeff(5,coeffs(idisp),isym)
!        Fill First term arrays 
         atindx(1,idisp) = ia; atindx(2,idisp) = ib;
         dir_char(idisp) = mutodir(mu); dir_int(idisp) = mu
         powers(idisp)   = 1
         cells(:,1,idisp) = (/0,0,0/)
         cells(:,2,idisp) = cell(:,irpt)
       end do
       compatible = .true.
!      Check the cut off and if the coeff is valid
       do ii=1,power
         do jj=ii+1,power
           do kk=2,3
             if(dist(xcart(:,list_coeff(kk,icoeff1,1)),&
&                    xcart(:,list_coeff(2,icoeff2,1)),rprimd,&
&                    cell(:,list_coeff(4,icoeff1,1)),&
&                    cell(:,list_coeff(4,icoeff2,1)))>cut_off.or.&
&               dist(xcart(:,list_coeff(kk,icoeff1,1)),&
&                    xcart(:,list_coeff(3,icoeff2,1)),rprimd,&
&                    cell(:,list_coeff(4,icoeff1,1)),&
&                    cell(:,list_coeff(4,icoeff2,1)))>cut_off)then
               compatible =.false.
             end if
           end do
         end do
       end do

       if(compatible)then
         iterm = iterm + 1
         call polynomial_term_init(atindx,cells,dir_int,ndisp,terms(iterm),powers,&
&                                  weight,check=.true.)
       end if
     end do!end do sym

     if(iterm > 0)then
!    increase coefficients and set it
       icoeff_tmp = icoeff_tmp + 1
       call polynomial_coeff_init(coefficient,iterm,coeffs_tmp(icoeff_tmp),terms(1:iterm),check=.true.)
     end if

!    Deallocate the terms
     do iterm=1,nterm_max
       call polynomial_term_free(terms(iterm))
     end do
   end do!end do coeff_sym1
 end do!end do coeff_sym2

 ABI_DEALLOCATE(terms)

 ABI_DEALLOCATE(coeffs)
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
!    Get the name of this coefficient
     call polynomial_coeff_getName(name,natom,coeffs_tmp(icoeff_tmp),symbols)
!    Increase icoeff and fill the coeffs_out array
     icoeff1 = icoeff1 + 1
     call polynomial_coeff_init(one,coeffs_tmp(icoeff_tmp)%nterm,&
&                               coeffs_out(icoeff1),coeffs_tmp(icoeff_tmp)%terms,&
&                               name=name)
    end if
 end do

 filename = "terms_2nd_order.xml"
! call polynomial_coeff_writeXML(coeffs_out,ncoeff_out,filename=filename)
 write(message,'(1x,I0,a)')&
&       ncoeff_out,' coefficients for the 2st order '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL') 

!Deallocation
 do icoeff1=1,ncoeff_max
   call polynomial_coeff_free(coeffs_tmp(icoeff1))
 end do
 ABI_DEALLOCATE(coeffs_tmp)

end subroutine fit_polynomial_coeff_getOrder2
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_getOrder3
!!
!! NAME
!! fit_polynomial_coeff_getOrder3
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

subroutine fit_polynomial_coeff_getOrder3(cell,coeffs_out,cut_off,list_coeff,&
&                                   natom,ncoeff_out,ncoeff,nrpt,nsym,&
&                                   rprimd,symbols,xcart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_getOrder3'
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
 integer :: ia,ib,icoeff1,icoeff2,icoeff3,icoeff_tmp
 integer :: ii,jj,kk,idisp,irpt,irpt_ref,isym,iterm
 integer :: mu,ncoeff_max,ndisp,nterm_max,power
 real(dp):: coefficient,weight
 logical :: compatible
!arrays
 integer,allocatable :: atindx(:,:),coeffs(:)
 integer,allocatable :: cells(:,:,:),dir_int(:)
 integer,allocatable :: powers(:)
 character(len=1) :: dir_char(3)
 character(len=1) :: mutodir(9) = (/"x","y","z","1","2","3","4","5","6"/)
 character(len=100):: name
 character(len=500) :: message
 character(len=fnlen) :: filename
 type(polynomial_term_type),dimension(:),allocatable :: terms
 type(polynomial_coeff_type),allocatable :: coeffs_tmp(:)

! *************************************************************************

!Initialisation of variables
 power = 3
 nterm_max  = nsym
 ncoeff_max = ncoeff**power
 ndisp = power
 ABI_ALLOCATE(coeffs_tmp,(ncoeff_max))
 ABI_ALLOCATE(terms,(nterm_max))
 ABI_ALLOCATE(atindx,(2,ndisp))
 ABI_ALLOCATE(cells,(3,2,ndisp))
 ABI_ALLOCATE(coeffs,(ndisp))
 ABI_ALLOCATE(dir_int,(ndisp))
 ABI_ALLOCATE(powers,(ndisp))

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
   do icoeff2=icoeff1,ncoeff
     do icoeff3=icoeff2,ncoeff
       iterm = zero
       coefficient = one
!      Set the coeffs
       coeffs(:) = (/icoeff1,icoeff2,icoeff3/) 
       do isym=1,nsym
!        Treat this coeff
         weight = 1
         do idisp=1,ndisp
!          Get index of this displacement term
           mu   = list_coeff(1,coeffs(idisp),isym)
           ia   = list_coeff(2,coeffs(idisp),isym)
           ib   = list_coeff(3,coeffs(idisp),isym)
           irpt = list_coeff(4,coeffs(idisp),isym)
           weight = weight*list_coeff(5,coeffs(idisp),isym)
!          Fill First term arrays 
           atindx(1,idisp) = ia; atindx(2,idisp) = ib;
           dir_char(idisp) = mutodir(mu); dir_int(idisp) = mu
           powers(idisp)   = 1
           cells(:,1,idisp) = (/0,0,0/)
           cells(:,2,idisp) = cell(:,irpt)
         end do
         compatible = .true.
!      Check the cut off and if the coeff is valid
         do ii=1,power
           do jj=ii+1,power
             do kk=2,3
               if(dist(xcart(:,list_coeff(kk,coeffs(ii),1)),&
&                      xcart(:,list_coeff(2,coeffs(jj),1)),rprimd,&
&                      cell(:,list_coeff(4,coeffs(ii),1)),&
&                      cell(:,list_coeff(4,coeffs(jj),1)))>cut_off.or.&
&                 dist(xcart(:,list_coeff(kk,coeffs(ii),1)),&
&                      xcart(:,list_coeff(3,coeffs(jj),1)),rprimd,&
&                      cell(:,list_coeff(4,coeffs(ii),1)),&
&                      cell(:,list_coeff(4,coeffs(jj),1)))>cut_off)then
                 compatible =.false.
               end if
             end do
           end do
         end do
           
         if(compatible)then
           iterm = iterm + 1
           call polynomial_term_init(atindx,cells,dir_int,ndisp,terms(iterm),powers,&
&                                    weight,check=.true.)
         end if
       end do!end do sym

       if(iterm > 0)then
!        increase coefficients and set it
         icoeff_tmp = icoeff_tmp + 1         
         call polynomial_coeff_init(coefficient,iterm,coeffs_tmp(icoeff_tmp),terms(1:iterm),check=.true.)
       end if

!      Deallocate the terms
       do iterm=1,nterm_max
         call polynomial_term_free(terms(iterm))
       end do
     end do!end do coeff_sym1
   end do!end do coeff_sym2
 end do!end do coeff_sym3

 ABI_DEALLOCATE(terms)
 ABI_DEALLOCATE(coeffs)
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
!    Get the name of this coefficient
     call polynomial_coeff_getName(name,natom,coeffs_tmp(icoeff_tmp),symbols)
!    Increase icoeff and fill the coeffs_out array
     icoeff1 = icoeff1 + 1
     call polynomial_coeff_init(one,coeffs_tmp(icoeff_tmp)%nterm,&
&                               coeffs_out(icoeff1),coeffs_tmp(icoeff_tmp)%terms,&
&                               name=name)
    end if
 end do

 filename = "terms_3rd_order.xml"
! call polynomial_coeff_writeXML(coeffs_out,ncoeff_out,filename=filename)
 write(message,'(1x,I0,a)')&
&       ncoeff_out,' coefficient for the 3rd order '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL') 

!Deallocation
 do icoeff1=1,ncoeff_max
   call polynomial_coeff_free(coeffs_tmp(icoeff1))
 end do
 ABI_DEALLOCATE(coeffs_tmp)

end subroutine fit_polynomial_coeff_getOrder3
!!***


!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_getOrder4
!!
!! NAME
!! fit_polynomial_coeff_getOrder4
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

subroutine fit_polynomial_coeff_getOrder4(cell,coeffs_out,cut_off,list_coeff,&
&                                   natom,ncoeff_out,ncoeff,nrpt,nsym,&
&                                   rprimd,symbols,xcart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_getOrder4'
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
 integer :: ia,ib,icoeff1,icoeff2,icoeff3,icoeff4,icoeff_tmp
 integer :: ii,jj,kk,idisp,irpt,irpt_ref,isym,iterm
 integer :: mu,ncoeff_max,ndisp,nterm_max,power
 real(dp):: coefficient,weight
 logical :: compatible
!arrays
 integer,allocatable :: atindx(:,:),coeffs(:)
 integer,allocatable :: cells(:,:,:),dir_int(:)
 integer,allocatable :: powers(:)
 character(len=1),allocatable :: dir_char(:)
 character(len=1) :: mutodir(9) = (/"x","y","z","1","2","3","4","5","6"/)
 character(len=100):: name
 character(len=500) :: message
 character(len=fnlen) :: filename
 type(polynomial_term_type),dimension(:),allocatable :: terms
 type(polynomial_coeff_type),allocatable :: coeffs_tmp(:)

! *************************************************************************

!Initialisation of variables
 power = 4
 nterm_max  = nsym
 ncoeff_max = ncoeff**power
 ndisp = power
 ABI_ALLOCATE(coeffs_tmp,(ncoeff_max))
 ABI_ALLOCATE(terms,(nterm_max))
 ABI_ALLOCATE(atindx,(2,ndisp))
 ABI_ALLOCATE(cells,(3,2,ndisp))
 ABI_ALLOCATE(coeffs,(ndisp))
 ABI_ALLOCATE(dir_int,(ndisp))
 ABI_ALLOCATE(powers,(ndisp))
 ABI_ALLOCATE(dir_char,(ndisp))

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
   do icoeff2=icoeff1,ncoeff
     do icoeff3=icoeff2,ncoeff
       do icoeff4=icoeff3,ncoeff
         iterm = zero
         coefficient = one
!        Set the coeffs
         coeffs(:) = (/icoeff1,icoeff2,icoeff3,icoeff4/) 
         do isym=1,nsym
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
           end do
           compatible = .true.
!        Check the cut off and if the coeff is valid
           do ii=1,power
             do jj=ii+1,power
               do kk=2,3
                 if(dist(xcart(:,list_coeff(kk,coeffs(ii),1)),&
&                        xcart(:,list_coeff(2,coeffs(jj),1)),rprimd,&
&                        cell(:,list_coeff(4,coeffs(ii),1)),&
&                        cell(:,list_coeff(4,coeffs(jj),1)))>cut_off.or.&
&                   dist(xcart(:,list_coeff(kk,coeffs(ii),1)),&
&                        xcart(:,list_coeff(3,coeffs(jj),1)),rprimd,&
&                        cell(:,list_coeff(4,coeffs(ii),1)),&
&                        cell(:,list_coeff(4,coeffs(jj),1)))>cut_off)then
                   compatible =.false.
                 end if
               end do
             end do
           end do
           
           if(compatible)then
             iterm = iterm + 1
             call polynomial_term_init(atindx,cells,dir_int,ndisp,terms(iterm),powers,&
&                                    weight,check=.true.)
           end if
         end do!end do sym

         if(iterm > 0)then
!        increase coefficients and set it
           icoeff_tmp = icoeff_tmp + 1         
           call polynomial_coeff_init(coefficient,iterm,coeffs_tmp(icoeff_tmp),terms(1:iterm),&
&                                     check=.true.)
         end if

!      Deallocate the terms
         do iterm=1,nterm_max
           call polynomial_term_free(terms(iterm))
         end do
       end do!end do coeff_sym1
     end do!end do coeff_sym2
   end do!end do coeff_sym3
 end do!end do coeff4
 ABI_DEALLOCATE(terms)
 ABI_DEALLOCATE(coeffs)
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(cells)
 ABI_DEALLOCATE(dir_int)
 ABI_DEALLOCATE(dir_char)
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
!    Get the name of this coefficient
     call polynomial_coeff_getName(name,natom,coeffs_tmp(icoeff_tmp),symbols)
!    Increase icoeff and fill the coeffs_out array
     icoeff1 = icoeff1 + 1
     call polynomial_coeff_init(one,coeffs_tmp(icoeff_tmp)%nterm,&
&                               coeffs_out(icoeff1),coeffs_tmp(icoeff_tmp)%terms,&
&                               name=name)
   end if
 end do

 filename = "terms_4th_order.xml"
! call polynomial_coeff_writeXML(coeffs_out,ncoeff_out,filename=filename)
 write(message,'(1x,I0,a)')&
&       ncoeff_out,' coefficients for the 4th order '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL') 

!Deallocation
 do icoeff1=1,ncoeff_max
   call polynomial_coeff_free(coeffs_tmp(icoeff1))
 end do
 ABI_DEALLOCATE(coeffs_tmp)

end subroutine fit_polynomial_coeff_getOrder4
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_getOrder5
!!
!! NAME
!! fit_polynomial_coeff_getOrder5
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

subroutine fit_polynomial_coeff_getOrder5(cell,coeffs_out,cut_off,list_coeff,&
&                                   natom,ncoeff_out,ncoeff,nrpt,nsym,&
&                                   rprimd,symbols,xcart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_getOrder5'
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
 integer :: ia,ib,icoeff1,icoeff2,icoeff3,icoeff4,icoeff5,icoeff_tmp
 integer :: ii,jj,kk,idisp,irpt,irpt_ref,isym,iterm
 integer :: mu,ncoeff_max,ndisp,nterm_max,power
 real(dp):: coefficient,weight
 logical :: compatible
!arrays
 integer,allocatable :: atindx(:,:),coeffs(:)
 integer,allocatable :: cells(:,:,:),dir_int(:)
 integer,allocatable :: powers(:)
 character(len=1),allocatable :: dir_char(:)
 character(len=1) :: mutodir(9) = (/"x","y","z","1","2","3","4","5","6"/)
 character(len=100):: name
 character(len=500) :: message
 character(len=fnlen) :: filename
 type(polynomial_term_type),dimension(:),allocatable :: terms
 type(polynomial_coeff_type),allocatable :: coeffs_tmp(:)

! *************************************************************************

!Initialisation of variables
 power = 5
 nterm_max  = nsym
 ncoeff_max = ncoeff**power
 ndisp = power
 ABI_ALLOCATE(coeffs_tmp,(ncoeff_max))
 ABI_ALLOCATE(terms,(nterm_max))
 ABI_ALLOCATE(atindx,(2,ndisp))
 ABI_ALLOCATE(cells,(3,2,ndisp))
 ABI_ALLOCATE(coeffs,(ndisp))
 ABI_ALLOCATE(dir_int,(ndisp))
 ABI_ALLOCATE(powers,(ndisp))
 ABI_ALLOCATE(dir_char,(ndisp))

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
   do icoeff2=icoeff1,ncoeff
     do icoeff3=icoeff2,ncoeff
       do icoeff4=icoeff3,ncoeff
         do icoeff5=icoeff4,ncoeff
           iterm = zero
           coefficient = one
!          Set the coeffs
           coeffs(:) = (/icoeff1,icoeff2,icoeff3,icoeff4,icoeff5/) 
           do isym=1,nsym
!            Treat this coeff
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
             end do
             compatible = .true.
!            Check the cut off and if the coeff is valid
             do ii=1,power
               do jj=ii+1,power
                 do kk=2,3
                   if(dist(xcart(:,list_coeff(kk,coeffs(ii),1)),&
&                          xcart(:,list_coeff(2,coeffs(jj),1)),rprimd,&
&                          cell(:,list_coeff(4,coeffs(ii),1)),&
&                          cell(:,list_coeff(4,coeffs(jj),1)))>cut_off.or.&
&                     dist(xcart(:,list_coeff(kk,coeffs(ii),1)),&
&                          xcart(:,list_coeff(3,coeffs(jj),1)),rprimd,&
&                          cell(:,list_coeff(4,coeffs(ii),1)),&
&                          cell(:,list_coeff(4,coeffs(jj),1)))>cut_off)then
                     compatible =.false.
                   end if
                 end do
               end do
             end do
           
             if(compatible)then
               iterm = iterm + 1
               call polynomial_term_init(atindx,cells,dir_int,ndisp,terms(iterm),powers,&
&                                        weight,check=.true.)
             end if
           end do!end do sym

           if(iterm > 0)then
!          increase coefficients and set it
             icoeff_tmp = icoeff_tmp + 1         
             call polynomial_coeff_init(coefficient,iterm,coeffs_tmp(icoeff_tmp),terms(1:iterm),&
&                                       check=.true.)
           end if

!          Deallocate the terms
           do iterm=1,nterm_max
             call polynomial_term_free(terms(iterm))
           end do
         end do!end do coeff_1
       end do!end do coeff_2
     end do!end do coeff_3
   end do!end do coeff4
 end do!end do coeff5
 ABI_DEALLOCATE(terms)
 ABI_DEALLOCATE(coeffs)
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(cells)
 ABI_DEALLOCATE(dir_int)
 ABI_DEALLOCATE(dir_char)
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
!    Get the name of this coefficient
     call polynomial_coeff_getName(name,natom,coeffs_tmp(icoeff_tmp),symbols)
!     Increase icoeff and fill the coeffs_out array
      icoeff1 = icoeff1 + 1
      call polynomial_coeff_init(one,coeffs_tmp(icoeff_tmp)%nterm,&
&                               coeffs_out(icoeff1),coeffs_tmp(icoeff_tmp)%terms,&
&                               name=name)
    end if
 end do

 filename = "terms_5th_order.xml"
! call polynomial_coeff_writeXML(coeffs_out,ncoeff_out,filename=filename)
 write(message,'(1x,I0,a)')&
&       ncoeff_out,' coefficients for the 5th order '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL') 

!Deallocation
 do icoeff1=1,ncoeff_max
   call polynomial_coeff_free(coeffs_tmp(icoeff1))
 end do
 ABI_DEALLOCATE(coeffs_tmp)

end subroutine fit_polynomial_coeff_getOrder5
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_fit
!!
!! NAME
!! fit_polynomial_coeff_fit
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

subroutine fit_polynomial_coeff_fit(cut_off,eff_pot,hist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_fit'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: cut_off
!arrays
 type(effective_potential_type),intent(inout) :: eff_pot
 type(abihist),intent(in) :: hist
!Local variables-------------------------------
!scalar
 integer :: ii,icycle,index_min,itime,jj
 integer :: ncoeff_max,natom_sc,ncycle,ntime
 real(dp) :: energy,ffact,sfact,mingf,mse,msef,mses,ucvol
!arrays
 character(len=500) :: message
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 integer,allocatable :: list_coeffs(:)
 real(dp),allocatable :: coeff_values(:)
 real(dp),allocatable :: displacement(:,:,:),fcart_fixed(:,:,:)
 real(dp),allocatable :: fred_fixed(:,:,:),fcart_HIST(:,:,:)
 real(dp),allocatable :: fcart_coeffs(:,:,:,:),gf_values(:),strain(:,:),strten_coeffs(:,:,:)
 real(dp),allocatable :: strten_fixed(:,:),strten_HIST(:,:),sqomega(:)
 type(polynomial_coeff_type),allocatable :: coeffs_tmp(:),coeffs_in(:)
 type(strain_type) :: strain_t
! *************************************************************************

 write(message,'(a,(80a))') ch10,('=',ii=1,80)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 write(message,'(2a)') ch10,' Starting Fit process'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 write(message,'(a,(80a))') ch10,('-',ii=1,80)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 if(eff_pot%anharmonics_terms%ncoeff > zero)then
   write(message, '(4a)' )ch10,' The coefficients present in the effective',&
&      ' potential will be used for the fit'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 else

!  Fit the coefficients, need MD file...
   write(message, '(4a)' )ch10,' The coefficients for the fit must',&
&                              ' be generate... ',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   call fit_polynomial_coeff_get(cut_off,eff_pot%anharmonics_terms%coefficients,&
&                                eff_pot,eff_pot%anharmonics_terms%ncoeff,1)
 end if


!Initialisation of constants
 ncoeff_max = eff_pot%anharmonics_terms%ncoeff
 ntime      = hist%mxhist
 natom_sc   = eff_pot%supercell%natom_supercell
 ffact      = 1.0/(ntime*natom_sc)
 sfact      = 1.0/(ntime*6)

!Copy the initial coefficients
 ABI_DATATYPE_ALLOCATE(coeffs_in,(ncoeff_max))
 do ii=1,ncoeff_max
   call polynomial_coeff_init(eff_pot%anharmonics_terms%coefficients(ii)%coefficient,&
&                             eff_pot%anharmonics_terms%coefficients(ii)%nterm,&
&                             coeffs_in(ii),&
&                             eff_pot%anharmonics_terms%coefficients(ii)%terms,&
&                             eff_pot%anharmonics_terms%coefficients(ii)%name,&
&                             check=.false.)
 
 end do

!MOVE IT TO INPUT IN THE FUTURE
 ncycle = 5
 ncycle = min(ncycle,ncoeff_max)

!Initialisation of arrays
 ABI_DATATYPE_ALLOCATE(coeffs_tmp,(ncycle))
 ABI_ALLOCATE(list_coeffs,(ncoeff_max))
 ABI_ALLOCATE(coeff_values,(ncoeff_max))
 ABI_ALLOCATE(displacement,(3,natom_sc,ntime))
 ABI_ALLOCATE(gf_values,(ncoeff_max))
 ABI_ALLOCATE(strain,(6,ntime))
 ABI_ALLOCATE(fcart_fixed,(3,natom_sc,ntime))
 ABI_ALLOCATE(fred_fixed,(3,natom_sc,ntime))
 ABI_ALLOCATE(strten_fixed,(6,ntime))
 ABI_ALLOCATE(fcart_HIST,(3,natom_sc,ntime))
 ABI_ALLOCATE(strten_HIST,(6,ntime))
 ABI_ALLOCATE(fcart_coeffs,(3,natom_sc,ntime,ncoeff_max))
 ABI_ALLOCATE(strten_coeffs,(6,ntime,ncoeff_max))
 ABI_ALLOCATE(sqomega,(ntime))

 mse  = zero
 msef = zero
 mses = zero
 list_coeffs  = zero
 displacement = zero
 strain       = zero
 fcart_fixed  = zero
 fred_fixed   = zero
 strten_fixed = zero
 fcart_HIST   = zero
 strten_HIST  = zero
 sqomega      = zero


!Before the fit, compute constants
!Conpute the strain of each configuration
!Compute the displacmeent of each configurations.
!Compute fixed forces and stresse and get the standard deviation
!Compute Shepard and al Factors  \Omega^{2} see J.Chem Phys 136, 074103 (2012)
 fcart_HIST(:,:,:) = hist%fcart(:,:,:)
 strten_HIST(:,:)= hist%strten(:,:)

 do itime=1,ntime
!  Get strain
   call strain_get(strain_t,rprim=eff_pot%supercell%rprimd_supercell,rprim_def=hist%rprimd(:,:,itime))
    if (strain_t%name /= "reference")  then
      do ii=1,3
        strain(ii,itime) = strain_t%strain(ii,ii)
      end do
      strain(4,itime) = strain_t%strain(2,3) + strain_t%strain(3,2)
      strain(5,itime) = strain_t%strain(3,1) + strain_t%strain(1,3)
      strain(6,itime) = strain_t%strain(2,1) + strain_t%strain(1,2)
    else
      strain(:,itime) = zero
    end if

!  Get displacement
   call effective_potential_getDisp(displacement(:,:,itime),natom_sc,hist%rprimd(:,:,itime),&
&                                         eff_pot%supercell%rprimd_supercell,&
&                                         xred_hist=hist%xred(:,:,itime),&
&                                         xcart_ref=eff_pot%supercell%xcart_supercell)

!  Get forces and stresses from harmonic part (fixed part)     
   call effective_potential_evaluate(eff_pot,energy,fcart_fixed(:,:,itime),fred_fixed(:,:,itime),&
&                                    strten_fixed(:,itime),natom_sc,hist%rprimd(:,:,itime),&
&                                    hist%xred(:,:,itime),displacement=displacement(:,:,itime),&
&                                    compute_anharmonic=.FALSE.,verbose=.FALSE.)

   mse = mse + abs(hist%etot(itime) - energy)

!  Compute \Omega^{2}
   call metric(gmet,gprimd,-1,rmet,hist%rprimd(:,:,itime),ucvol)
   sqomega(itime) = ((ucvol * (natom_sc)**0.5)**(-1./3))**2
 end do

!Print the standard deviation before the fit
 write(message, '(2a,E24.16,a)' )ch10,' Standard deviation of the energy at the begining (meV/f.u.): ',&
&               mse* 1000*27.21138386 / ntime / product(eff_pot%supercell%qphon(:)),ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')


!Get the decomposition for each coefficients of the forces and stresses for 
!each atoms and each step  equations 11 & 12 of  PRB95,094115(2017) 
 call fit_polynomial_coeff_getFS(eff_pot%anharmonics_terms%coefficients,displacement,ffact,&
&                                fcart_coeffs,natom_sc,eff_pot%crystal%natom,ncoeff_max,ntime,&
&                                int(eff_pot%supercell%qphon(:)),strain,strten_coeffs,sfact,&
&                                eff_pot%mpi_coeff%my_cells,eff_pot%mpi_coeff%my_ncell,&
&                                eff_pot%mpi_coeff%my_index_cells,eff_pot%mpi_coeff%comm)

!Start fit process
 do icycle=1,ncycle

!  Fill the coeffs list with the previous coeffs
   if(icycle > 1) then
     write(message, '(5a)')'--',ch10,' Coefficient numbers from the previous cycle:',ch10,' ['
     do ii=1,icycle-1
       if(ii<icycle-1)then
         write(message, '(a,I0,a)') trim(message),list_coeffs(ii),','
       else
         write(message, '(a,I0)') trim(message),list_coeffs(ii)
       end if
     end do
     write(message, '(2a)') trim(message),']'
     call wrtout(std_out,message,'COLL')
   else
     write(message, '(4a)')'--',ch10,' No coefficient numbers from the previous cycle.',&
&                                    ' Start from scratch.'
     call wrtout(std_out,message,'COLL')
   end if
   
!  Reset gf_values
   gf_values(:) = zero

   do jj=1,ncoeff_max
     if(any(list_coeffs==jj)) cycle
     list_coeffs(icycle) = jj

!    call the fit process routine
!    This routine solves the linear system proposed by C.Escorihuela-Sayalero see PRB95,094115(2017)
     call fit_polynomial_coeff_solve(coeff_values(1:icycle),fcart_coeffs,fcart_fixed,fcart_HIST,ffact,&
&                                    list_coeffs(1:icycle),natom_sc,icycle,ncoeff_max,ntime,&
&                                    strten_coeffs,strten_fixed,strten_HIST,sfact,sqomega)


     call fit_polynomial_coeff_ComputeGF(coeff_values(1:icycle),fcart_coeffs,fcart_fixed,fcart_HIST,&
&                                       ffact,gf_values(icycle),list_coeffs(1:icycle),natom_sc,icycle,&
&                                       ncoeff_max,ntime,strten_coeffs,strten_fixed,strten_HIST,&
&                                       sfact,sqomega)

     write(message, '(a,I0,1a,E24.16)')' Testing coefficient ',jj,' GF value(meV/A^2): ',&
&                                       gf_values(icycle)*27.21138386*1000/(0.529177**2)
     call wrtout(std_out,message,'COLL')

   end do

!  Find a way to keep the bestStore the best coeff for the next step
   mingf    = 9D99
   index_min= one
   do jj=1,ncoeff_max
     if(gf_values(jj) < zero) cycle
     if(gf_values(jj) == zero) cycle
     if(gf_values(jj) < mingf ) then
       mingf = gf_values(jj)
       index_min = jj
     end if
   end do
   list_coeffs(icycle) = index_min

!  Store the best coeff 
   call polynomial_coeff_init(coeff_values(icycle),&
&                             coeffs_in(list_coeffs(icycle))%nterm,&
&                             coeffs_tmp(icycle),&
&                             coeffs_in(list_coeffs(icycle))%terms,&
&                             coeffs_in(list_coeffs(icycle))%name,&
&                             check=.false.)

!  Set the new set of coefficients into the eff_pot type
   call effective_potential_setCoeffs(coeffs_tmp(1:icycle),eff_pot,icycle)

!  Compute the MSE,MSEF,MSES for this new anharmonic part
   call fit_polynomial_coeff_computeMSE(eff_pot,hist,mse,msef,mses,natom_sc,compute_anharmonic=.TRUE.)
   write(message, '(2a,I0,a,E24.16)' )' Standard deviation of the energy for',&
&                                      ' the iteration ',icycle,' (meV/f.u.): ',&
&                        mse* 1000*27.21138386 / ntime / product(eff_pot%supercell%qphon(:))
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 end do

!Fit the last model
!This routine solves the linear system proposed by C.Escorihuela-Sayalero see PRB95,094115(2017)
 call fit_polynomial_coeff_solve(coeff_values(1:ncycle),fcart_coeffs,fcart_fixed,fcart_HIST,ffact,&
&                                list_coeffs(1:ncycle),natom_sc,ncycle,ncoeff_max,ntime,&
&                                strten_coeffs,strten_fixed,strten_HIST,sfact,sqomega)

 write(message, '(2a)') ch10,' Fitted coefficients at the end of the fit process: '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 do ii = 1,ncoeff_max   
   if(list_coeffs(ii) ==0) cycle
   write(message, '(a,I0,a,E19.10,2a)') " ",list_coeffs(ii)," =>",coeff_values(ii),&
&                              " ",trim(coeffs_in(ii)%name)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end do


!Deallocation of arrays
!Deallocate the temporary coefficient
 do ii=1,ncycle
   call polynomial_coeff_free(coeffs_tmp(ii))
 end do
 do ii=1,ncoeff_max
   call polynomial_coeff_free(coeffs_in(ii))
 end do
 ABI_DATATYPE_DEALLOCATE(coeffs_tmp)
 ABI_DATATYPE_DEALLOCATE(coeffs_in)
 ABI_DEALLOCATE(gf_values)
 ABI_DEALLOCATE(coeff_values)
 ABI_DEALLOCATE(displacement)
 ABI_DEALLOCATE(fcart_fixed)
 ABI_DEALLOCATE(fred_fixed)
 ABI_DEALLOCATE(fcart_HIST)
 ABI_DEALLOCATE(fcart_coeffs)
 ABI_DEALLOCATE(strain)
 ABI_DEALLOCATE(strten_HIST)
 ABI_DEALLOCATE(strten_fixed)
 ABI_DEALLOCATE(strten_coeffs)
 ABI_DEALLOCATE(list_coeffs)
 ABI_DEALLOCATE(sqomega)

end subroutine fit_polynomial_coeff_fit
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_solve
!!
!! NAME
!! fit_polynomial_coeff_solve
!!
!! FUNCTION
!!
!!
!! INPUTS
!!
!! OUTPUT
!!
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_polynomial_coeff_solve(coefficients,fcart_coeffs,fcart_fixed,fcart_HIST,ffact,&
&                                     list_coeffs,natom,ncoeff_fit,ncoeff_max,ntime,strten_coeffs,&
&                                     strten_fixed,strten_HIST,sfact,sqomega)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_solve'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: natom,ncoeff_fit,ncoeff_max,ntime
 real(dp),intent(in) :: ffact,sfact
!arrays
 integer,intent(in)  :: list_coeffs(ncoeff_fit)
 real(dp),intent(in) :: fcart_coeffs(3,natom,ntime,ncoeff_max)
 real(dp),intent(in) :: fcart_fixed(3,natom,ntime),fcart_HIST(3,natom,ntime)
 real(dp),intent(in) :: strten_coeffs(6,ntime,ncoeff_max),strten_fixed(6,ntime)
 real(dp),intent(in) :: strten_HIST(6,ntime),sqomega(ntime)
 real(dp),intent(out):: coefficients(ncoeff_fit)
!Local variables-------------------------------
!scalar
 integer :: ia,itime,icoeff,jcoeff,icoeff_tmp,jcoeff_tmp,mu,LDA,LDB,LDX,LDAF,N,NRHS
 real(dp):: ftmpA,stmpA,ftmpB,stmpB,fmu,fnu,smu,snu
 integer :: INFO
 real(dp):: RCOND
 real(dp),allocatable:: AF(:,:),BERR(:),FERR(:),WORK(:),C(:),R(:)
 integer,allocatable :: IPIV(:),IWORK(:)
!arrays
 real(dp),allocatable :: A(:,:),B(:,:)
! character(len=500) :: message
! *************************************************************************

!0-Set variables for the 
 N    = ncoeff_fit; LDA  = ncoeff_fit; LDB  = ncoeff_fit; LDX  = ncoeff_fit
 LDAF = ncoeff_fit; NRHS = 1; RCOND = zero; INFO  = zero

!0-Allocation
 ABI_ALLOCATE(A,(LDA,N))
 ABI_ALLOCATE(B,(LDB,NRHS))
 ABI_ALLOCATE(AF,(LDAF,N))
 ABI_ALLOCATE(IPIV,(N))
 ABI_ALLOCATE(R,(N))
 ABI_ALLOCATE(C,(N))
 ABI_ALLOCATE(FERR,(NRHS))
 ABI_ALLOCATE(BERR,(NRHS))
 ABI_ALLOCATE(WORK,(4*N))
 ABI_ALLOCATE(IWORK,(N))
 A=zero; B=zero;
 AF = zero; IPIV = one; 
 R = one; C = one; 
 FERR = zero; BERR = zero
 IWORK = zero; WORK = zero

!1-Get forces and stresses from the model and fill A
!  Fill alsor B with the forces and stresses from 
!  the DFT snapshot and the model
 do icoeff=1,ncoeff_fit
   do jcoeff=1,ncoeff_fit
     icoeff_tmp = list_coeffs(icoeff)
     jcoeff_tmp = list_coeffs(jcoeff)
     ftmpA= zero; ftmpB = zero
     stmpA= zero; stmpB = zero
!    loop over the configuration
     do itime=1,ntime
!      Fill forces
       do ia=1,natom
         do mu=1,3          
           fmu = fcart_coeffs(mu,ia,itime,icoeff_tmp)
           fnu = fcart_coeffs(mu,ia,itime,jcoeff_tmp)
           ftmpA = ftmpA + fmu*fnu
           ftmpB = ftmpB + (fcart_HIST(mu,ia,itime)-fcart_fixed(mu,ia,itime))*fmu 
         end do !End loop dir
       end do !End loop natom

!      Fill stresses
       do mu=1,6
         smu = strten_coeffs(mu,itime,icoeff_tmp)
         snu = strten_coeffs(mu,itime,jcoeff_tmp)
         stmpA = stmpA + sqomega(itime)*smu*snu
         stmpB = stmpB + sqomega(itime)*(strten_HIST(mu,itime)-strten_fixed(mu,itime))*smu 
       end do !End loop stress dir
     end do ! End loop time
     A(icoeff,jcoeff) = A(icoeff,jcoeff) + ffact*ftmpA + sfact*stmpA
     B(icoeff,1) = B(icoeff,1) + ffact*ftmpB + sfact*stmpB
   end do ! End loop jcoeff
 end do ! End loop icoeff

!2-Solve Ax=B
 call dgesvx('F','N',N,NRHS,A,LDA,A,LDA,IPIV,'N',R,C,B,LDB,coefficients,LDX,&
             RCOND,FERR,BERR,WORK,IWORK,INFO)

 ABI_DEALLOCATE(AF)
 ABI_DEALLOCATE(IPIV)
 ABI_DEALLOCATE(R)
 ABI_DEALLOCATE(C)
 ABI_DEALLOCATE(FERR)
 ABI_DEALLOCATE(BERR)
 ABI_DEALLOCATE(WORK)
 ABI_DEALLOCATE(IWORK)
 ABI_DEALLOCATE(A)
 ABI_DEALLOCATE(B)

end subroutine fit_polynomial_coeff_solve
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_ComputeGF
!!
!! NAME
!! fit_polynomial_coeff_ComputeGF
!!
!! FUNCTION
!!
!!
!! INPUTS
!!
!! OUTPUT
!!
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_polynomial_coeff_ComputeGF(coefficients,fcart_coeffs,fcart_fixed,fcart_HIST,ffact,&
&                                         gf_value,list_coeffs,natom,ncoeff_fit,ncoeff_max,ntime,&
&                                         strten_coeffs,strten_fixed,strten_HIST,sfact,sqomega)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_ComputeGF'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: natom,ncoeff_fit,ncoeff_max,ntime
 real(dp),intent(in) :: ffact,sfact
!arrays
 integer,intent(in)  :: list_coeffs(ncoeff_fit)
 real(dp),intent(in) :: fcart_coeffs(3,natom,ntime,ncoeff_max)
 real(dp),intent(in) :: fcart_fixed(3,natom,ntime),fcart_HIST(3,natom,ntime)
 real(dp),intent(in) :: strten_coeffs(6,ntime,ncoeff_max),strten_fixed(6,ntime)
 real(dp),intent(in) :: strten_HIST(6,ntime),sqomega(ntime)
 real(dp),intent(in) :: coefficients(ncoeff_fit)
 real(dp),intent(out) :: gf_value
!Local variables-------------------------------
!scalar
 integer :: ia,itime,icoeff,icoeff_tmp,mu
 real(dp):: fmu,ftmp,smu,stmp
!arrays
! *************************************************************************

!1-Compute the value of the goal function
! see equation 9 of PRB 95 094115(2017

 gf_value = zero
 ftmp     = zero
 stmp     = zero

! loop over the configuration
 do itime=1,ntime
! Fill forces
   do ia=1,natom
     do mu=1,3          
       fmu = zero
       do icoeff=1,ncoeff_fit
         icoeff_tmp = list_coeffs(icoeff)
         fmu = coefficients(icoeff)*fcart_coeffs(mu,ia,itime,icoeff_tmp)
       end do
       ftmp = ftmp + (fcart_HIST(mu,ia,itime)-fcart_fixed(mu,ia,itime)-fmu)
     end do !End loop dir
   end do !End loop natom
   do mu=1,6
     do icoeff=1,ncoeff_fit
       icoeff_tmp = list_coeffs(icoeff)
       smu = coefficients(icoeff)*strten_coeffs(mu,itime,icoeff_tmp)
     end do
     stmp = stmp + sqomega(itime)*(strten_HIST(mu,itime)-strten_fixed(mu,itime)-smu) 
   end do !End loop stress dir
 end do ! End loop time
 
 gf_value   =  ffact*ftmp + sfact*stmp
 


end subroutine fit_polynomial_coeff_ComputeGF
!!***


!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_getFS
!!
!! NAME
!! fit_polynomial_coeff_getFS
!!
!! FUNCTION
!!
!!
!! INPUTS
!!
!! OUTPUT
!!
!!
!! PARENTS
!!      m_fit_polynomial_coeff
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_polynomial_coeff_getFS(coefficients,displacement,ffact,fcart_out,&
&                                     natom_sc,natom_uc,ncoeff,ntime,sc_size,strain,strten_out,sfact,&
&                                     cells,ncell,index_cells,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_getFS'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: ffact,sfact
 integer,intent(in) :: natom_sc,natom_uc,ncoeff,ntime
 integer,intent(in) :: comm,ncell
!arrays
 integer,intent(in) :: sc_size(3)
 integer,intent(in) :: cells(ncell),index_cells(ncell,3)
 real(dp),intent(in) :: displacement(3,natom_sc,ntime)
 real(dp),intent(in) :: strain(6,ntime)
 real(dp),intent(out) :: fcart_out(3,natom_sc,ntime,ncoeff)
 real(dp),intent(out) :: strten_out(6,ntime,ncoeff)
 type(polynomial_coeff_type), intent(in) :: coefficients(ncoeff)
!Local variables-------------------------------
!scalar
 integer :: i1,i2,i3,ia1,ia2,ib1,ib2,ii,icell,icoeff
 integer :: idir1,idir2,idisp1,idisp2,ierr,iterm,itime,power
 real(dp):: disp1,disp2,tmp2,tmp3,weight
!arrays
 integer :: cell_atoma1(3),cell_atoma2(3)
 integer :: cell_atomb1(3),cell_atomb2(3)

! *************************************************************************


!1-Get forces and stresses from the model and fill A
!  Fill alsor B with the forces and stresses from 
!  the DFT snapshot and the model
! Initialisation of variables
 fcart_out(:,:,:,:) = zero
 strten_out(:,:,:)  = zero

 do icell = 1,ncell
   ii = (cells(icell)-1)*natom_uc
   i1=index_cells(icell,1); i2=index_cells(icell,2); i3=index_cells(icell,3)
!  Loop over configurations
   do itime=1,ntime
!    Loop over coefficients
     do icoeff=1,ncoeff
!      Loop over terms of this coefficient
       do iterm=1,coefficients(icoeff)%nterm
!        Set the weight of this term
         weight =coefficients(icoeff)%terms(iterm)%weight
!        Loop over displacement and strain
         do idisp1=1,coefficients(icoeff)%terms(iterm)%ndisp

!          Set to one the acculation of forces and strain
           tmp2 = one
           tmp3 = one

!          Set the power of the displacement:
           power = coefficients(icoeff)%terms(iterm)%power(idisp1)

!          Get the direction of the displacement or strain
           idir1 = coefficients(icoeff)%terms(iterm)%direction(idisp1)

!          Strain case idir => -6, -5, -4, -3, -2 or -1
           if (idir1 < zero)then

             if(abs(strain(abs(idir1),itime)) > tol10)then
               if(power > 1) then
!                Accumulate stress for each strain (\sum (Y(eta_2)^Y-1(eta_2)^Z+...))
                 tmp3 = tmp3 *  power*(strain(abs(idir1),itime))**(power-1)
               end if
             else
               if(power > 1) then
                 tmp3 = zero
               end if
             end if
           else
!            Displacement case idir = 1, 2  or 3
!            indexes of the cell of the atom a
             cell_atoma1 = coefficients(icoeff)%terms(iterm)%cell(:,1,idisp1)
             if(cell_atoma1(1)/=0.or.cell_atoma1(2)/=0.or.cell_atoma1(3)/=0) then
!                if the cell is not 0 0 0 we apply PBC:
               cell_atoma1(1) =  (i1-1) + cell_atoma1(1)
               cell_atoma1(2) =  (i2-1) + cell_atoma1(2)
               cell_atoma1(3) =  (i3-1) + cell_atoma1(3)
               call effective_potential_GetIndexPeriodic(cell_atoma1(1:3),sc_size(1:3))
!              index of the first atom (position in the supercell if the cell is not 0 0 0)
               ia1 = cell_atoma1(1)*sc_size(2)*sc_size(3)*natom_uc+&
&                    cell_atoma1(2)*sc_size(3)*natom_uc+&
&                    cell_atoma1(3)*natom_uc+&
&                    coefficients(icoeff)%terms(iterm)%atindx(1,idisp1)
             else
!              index of the first atom (position in the supercell if the cell is 0 0 0)
               ia1 = ii + coefficients(icoeff)%terms(iterm)%atindx(1,idisp1)
             end if

!            indexes of the cell of the atom b  (with PBC) same as ia1
             cell_atomb1 = coefficients(icoeff)%terms(iterm)%cell(:,2,idisp1)
             if(cell_atomb1(1)/=0.or.cell_atomb1(2)/=0.or.cell_atomb1(3)/=0) then
               cell_atomb1(1) =  (i1-1) + cell_atomb1(1)
               cell_atomb1(2) =  (i2-1) + cell_atomb1(2)
               cell_atomb1(3) =  (i3-1) + cell_atomb1(3)
               call effective_potential_GetIndexPeriodic(cell_atomb1(1:3),sc_size(1:3))

!              index of the second atom in the (position in the supercell  if the cell is not 0 0 0) 
               ib1 = cell_atomb1(1)*sc_size(2)*sc_size(3)*natom_uc+&
&                    cell_atomb1(2)*sc_size(3)*natom_uc+&
&                    cell_atomb1(3)*natom_uc+&
&                    coefficients(icoeff)%terms(iterm)%atindx(2,idisp1)
             else
!              index of the first atom (position in the supercell if the cell is 0 0 0)
               ib1 = ii + coefficients(icoeff)%terms(iterm)%atindx(2,idisp1)
             end if

!            Get the displacement for the both atoms
             disp1 = displacement(idir1,ia1,itime)
             disp2 = displacement(idir1,ib1,itime)

             if(abs(disp1) > tol10 .or. abs(disp2)> tol10)then
               if(power > 1) then
!                Accumulate forces for each displacement (\sum (Y(A_x-O_x)^Y-1(A_y-O_c)^Z+...))
                 tmp2 = tmp2 * power*(disp1-disp2)**(power-1)
               end if
             else
               if(power > 1) then
                 tmp2 = zero
               end if
             end if
           end if
             
           do idisp2=1,coefficients(icoeff)%terms(iterm)%ndisp
               
             if(idisp2 /= idisp1) then
               idir2 = coefficients(icoeff)%terms(iterm)%direction(idisp2)
               if (idir2 < zero)then
!               Strain case
!               Set the power of the strain:
                 power = coefficients(icoeff)%terms(iterm)%power(idisp2)
!                Accumulate energy forces
                 tmp2 = tmp2 * (strain(abs(idir2),itime))**power
!                Accumulate stress for each strain (\sum (Y(eta_2)^Y-1(eta_2)^Z+...))
                 tmp3 = tmp3 * (strain(abs(idir2),itime))**power
               else
                 cell_atoma2=coefficients(icoeff)%terms(iterm)%cell(:,1,idisp2)
                 if(cell_atoma2(1)/=0.or.cell_atoma2(2)/=0.or.cell_atoma2(3)/=0) then
                   cell_atoma2(1) =  (i1-1) + cell_atoma2(1)
                   cell_atoma2(2) =  (i2-1) + cell_atoma2(2)
                   cell_atoma2(3) =  (i3-1) + cell_atoma2(3)
                   call effective_potential_GetIndexPeriodic(cell_atoma2(1:3),sc_size(1:3))
!                  index of the first atom (position in the supercell and direction)
!                  if the cell of the atom a is not 0 0 0 (may happen)
                   ia2 = cell_atoma2(1)*sc_size(2)*sc_size(3)*natom_uc+&
&                        cell_atoma2(2)*sc_size(3)*natom_uc+&
&                        cell_atoma2(3)*natom_uc+&
&                        coefficients(icoeff)%terms(iterm)%atindx(1,idisp2)
                 else
!                  index of the first atom (position in the supercell and direction)
                   ia2 = ii + coefficients(icoeff)%terms(iterm)%atindx(1,idisp2)
                 end if

                 cell_atomb2= coefficients(icoeff)%terms(iterm)%cell(:,2,idisp2)
                 
                 if(cell_atomb2(1)/=0.or.cell_atomb2(2)/=0.or.cell_atomb2(3)/=0) then
!                  indexes of the cell2 (with PBC)
                   cell_atomb2(1) =  (i1-1) + cell_atomb2(1)
                   cell_atomb2(2) =  (i2-1) + cell_atomb2(2)
                   cell_atomb2(3) =  (i3-1) + cell_atomb2(3)
                   call effective_potential_GetIndexPeriodic(cell_atomb2(1:3),sc_size(1:3))

!                  index of the second atom in the (position in the supercell) 
                   ib2 = cell_atomb2(1)*sc_size(2)*sc_size(3)*natom_uc+&
&                        cell_atomb2(2)*sc_size(3)*natom_uc+&
&                        cell_atomb2(3)*natom_uc+&
&                        coefficients(icoeff)%terms(iterm)%atindx(2,idisp2)
                 else
                   ib2 = ii + coefficients(icoeff)%terms(iterm)%atindx(2,idisp2)
                 end if

                 disp1 = displacement(idir2,ia2,itime)
                 disp2 = displacement(idir2,ib2,itime)

!                Set the power of the displacement:
                 power = coefficients(icoeff)%terms(iterm)%power(idisp2)
                   
                 tmp2 = tmp2 * (disp1-disp2)**power
                 tmp3 = tmp3 * (disp1-disp2)**power

               end if
             end if
           end do

           if(idir1<zero)then
!            Accumule stress tensor
             strten_out(abs(idir1),itime,icoeff) = strten_out(abs(idir1),itime,icoeff) +  weight * tmp3
           else
!            Accumule  forces
             fcart_out(idir1,ia1,itime,icoeff) = fcart_out(idir1,ia1,itime,icoeff) + weight * tmp2
             fcart_out(idir1,ib1,itime,icoeff) = fcart_out(idir1,ib1,itime,icoeff) - weight * tmp2
           end if
         end do
        
       end do
     end do!End do coeff
   end do!End time
 end do!End do cell

!MPI_SUM
 call xmpi_sum(fcart_out, comm, ierr)
 call xmpi_sum(strten_out, comm, ierr)

!multiply derivative by -1 
 fcart_out  = -1 * fcart_out
 strten_out = -1 * strten_out

end subroutine fit_polynomial_coeff_getFS
!!***


!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_computeMSE
!!
!! NAME
!! fit_polynomial_coeff_computeMSE
!!
!! FUNCTION

!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_polynomial_coeff_computeMSE(eff_pot,hist,mse,msef,mses,natom,compute_anharmonic)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_computeMSE'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 real(dp),intent(out):: mse,msef,mses
 logical,optional,intent(in) :: compute_anharmonic
!arrays
 type(effective_potential_type),intent(in) :: eff_pot
 type(abihist),intent(in) :: hist
!Local variables-------------------------------
!scalar
 integer :: ii,nstep,ntime
 real(dp):: energy
 logical :: need_anharmonic = .TRUE.
!arrays
 real(dp):: fcart(3,natom),fred(3,natom),strten(6),rprimd(3,3),xred(3,natom)

! *************************************************************************

 mse  = zero
 msef = zero
 mses = zero

 nstep= 1
 ntime = hist%mxhist

 if(present(compute_anharmonic))then
   need_anharmonic = compute_anharmonic
 end if

 do ii=1,ntime
   xred(:,:)   = hist%xred(:,:,ii)
   rprimd(:,:) = hist%rprimd(:,:,ii)
   call effective_potential_evaluate(eff_pot,energy,fcart,fred,strten,natom,rprimd,&
&                                    xred,compute_anharmonic=need_anharmonic,verbose=.false.)

   mse = mse + abs(hist%etot(ii) - energy)
  
 end do

 mse = mse  / ntime 

end subroutine fit_polynomial_coeff_computeMSE
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_mapHistToRef
!!
!! NAME
!! fit_polynomial_coeff_mapHistToRef
!!
!! FUNCTION
!! Generate the supercell in the effective potential according to the size of the 
!! supercell in the hist file
!! Check if the hist file match to reference supercell in the effective potential
!! If not, the hist file is reordering 
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine fit_polynomial_coeff_mapHistToRef(eff_pot,hist,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_mapHistToRef'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
!arrays
 type(effective_potential_type),intent(inout) :: eff_pot
 type(abihist),intent(inout) :: hist
!Local variables-------------------------------
!scalar
 integer :: ia,ib,ii,jj,natom_hist,nstep_hist
 real(dp):: factor
 logical :: revelant_factor,need_map
!arrays
 real(dp) :: rprimd_hist(3,3),rprimd_ref(3,3),scale_cell(3)
 integer :: n_cell(3)
 integer,allocatable  :: blkval(:),list(:)
 real(dp),allocatable :: xred_hist(:,:),xred_ref(:,:)
 character(len=500) :: msg
 type(abihist) :: hist_tmp
! *************************************************************************

 natom_hist = size(hist%xred,2)
 nstep_hist = size(hist%xred,3)

!Try to set the supercell according to the hist file
 rprimd_ref(:,:)  = eff_pot%crystal%rprimd
 rprimd_hist(:,:) = hist%rprimd(:,:,1)
 do ia=1,3
   scale_cell(:) = zero
   do ii=1,3
     if(abs(rprimd_ref(ii,ia)) > tol10)then
       scale_cell(ii) = rprimd_hist(ii,ia) / rprimd_ref(ii,ia)
     end if
   end do
!  Check if the factor for the supercell is revelant
   revelant_factor = .TRUE.
   do ii=1,3
     if(abs(scale_cell(ii)) < tol10) cycle
     factor = abs(scale_cell(ii))
     do jj=ii,3
       if(abs(scale_cell(jj)) < tol10) cycle
       if(abs(abs(scale_cell(ii))-abs(scale_cell(jj))) > tol10) revelant_factor = .FALSE.
     end do
   end do
   if(.not.revelant_factor)then
     write(msg, '(3a)' )&
&         'unable to map the hist file ',ch10,&
&         'Action: check/change your MD file'
     MSG_ERROR(msg)
   else
     n_cell(ia) = nint(factor)
   end if
 end do
 
!Set the new supercell structure into the effective potential reference
 call effective_potential_setSupercell(eff_pot,comm,n_cell)

!allocation
 ABI_ALLOCATE(blkval,(natom_hist))
 ABI_ALLOCATE(list,(natom_hist))
 ABI_ALLOCATE(xred_hist,(3,natom_hist))
 ABI_ALLOCATE(xred_ref,(3,natom_hist))
 blkval = one

 call xcart2xred(eff_pot%supercell%natom_supercell,eff_pot%supercell%rprimd_supercell,&
&                eff_pot%supercell%xcart_supercell,xred_ref)

 xred_hist = hist%xred(:,:,1)

 write(msg,'(2a,I2,a,I2,a,I2)') ch10,&
&       ' The size of the supercell for the fit is ',n_cell(1),' ',n_cell(2),' ',n_cell(3)
 call wrtout(std_out,msg,'COLL') 
 call wrtout(ab_out,msg,'COLL') 

!try to map
 do ia=1,natom_hist
   do ib=1,natom_hist
     if(blkval(ib)==1)then
       if(abs(xred_ref(1,ia)-xred_hist(1,ib)) < tol10 .and.abs(xred_ref(2,ia)-xred_hist(2,ib)) < tol10 &
&      .and.abs(xred_ref(3,ia)-xred_hist(3,ib)) < tol10) then
         blkval(ib) = zero
         list(ib) = ia
       end if
     end if
   end do
 end do

!Check before transfert
 if(.not.all(blkval==zero))then
   write(msg, '(5a)' )&
&         'The hist file does correspond ',ch10,&
&         'to the reference supercell structure',ch10,&
&         'Action: change MD file'
     MSG_ERROR(msg)
 end if

 do ia=1,natom_hist
   if(.not.any(list(:)==ia))then
     write(msg, '(5a)' )&
&         'The hist file does correspond ',ch10,&
&         'to the reference supercell structure',ch10,&
&         'Action: change MD file'
     MSG_ERROR(msg)
   end if
 end do

 need_map = .FALSE.
 do ia=1,natom_hist
   if(list(ia) /= ia) need_map = .TRUE.
 end do
 if(need_map)then
   write(msg, '(11a)' )ch10,&
&      ' --- !WARNING',ch10,&
&      '     The ordering of the atoms in the hist file is different,',ch10,&
&      '     of the one built by multibinit. The hist file will be map,',ch10,&
&      '     on the ordering of multibinit.',ch10,&
&      ' ---',ch10
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')


! Allocate hist structure 
   call abihist_init(hist_tmp,natom_hist,nstep_hist,.false.,.false.)
! copy all the information
   do ia=1,nstep_hist
     hist%ihist = ia
     hist_tmp%ihist = ia
     call abihist_copy(hist,hist_tmp)
   end do
   hist_tmp%mxhist = nstep_hist

! reoder array
   do ia=1,natom_hist
     hist_tmp%xred(:,list(ia),:)=hist%xred(:,ia,:)
     hist_tmp%fcart(:,list(ia),:)=hist%fcart(:,ia,:)
     hist_tmp%vel(:,list(ia),:)=hist%vel(:,ia,:)
   end do

! free the old hist and reinit
   call abihist_free(hist)
   call abihist_init(hist,natom_hist,nstep_hist,.false.,.false.)
! copy the temporary hist into output
   do ia=1,nstep_hist
     hist%ihist = ia
     hist_tmp%ihist = ia
     call abihist_copy(hist_tmp,hist)
   end do
   hist_tmp%mxhist = nstep_hist

   call abihist_free(hist_tmp)
 end if

!deallocation
 ABI_DEALLOCATE(blkval)
 ABI_DEALLOCATE(list)
 ABI_DEALLOCATE(xred_hist)
 ABI_DEALLOCATE(xred_ref)

end subroutine fit_polynomial_coeff_mapHistToRef
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
!!      multibinit
!!
!! CHILDREN
!!      destroy_supercell,init_supercell,xred2xcart
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
 integer :: ia,ib,ib1,ii,jj,irpt,kk,ll,mu,nu,nstep,nshift
 integer :: natom_uc
 integer :: unit_born=22,unit_epsiloninf=23,unit_md=24
 integer :: unit_harmonic=25,unit_ref=26,unit_strain=27,unit_sym=28
!arrays
 integer,allocatable :: typat_order(:),typat_order_uc(:)
 integer, dimension(3)  :: A,ncell
 real(dp), allocatable :: xcart(:,:),fcart(:,:)
 character(len=500) :: msg
 type(supercell_type) :: supercell
! *************************************************************************

!Create new supercell corresponding to the MD
 ncell = (/2,2,2/)
 call init_supercell(eff_pot%crystal%natom, 0, real(ncell,dp), eff_pot%crystal%rprimd,&
&                    eff_pot%crystal%typat,eff_pot%crystal%xcart, supercell)

!allocation of array
 ABI_ALLOCATE(xcart,(3,supercell%natom_supercell))
 ABI_ALLOCATE(fcart,(3,supercell%natom_supercell))
 ABI_ALLOCATE(typat_order,(supercell%natom_supercell))
 ABI_ALLOCATE(typat_order_uc,(eff_pot%crystal%natom))

 A = (/ 2, 3, 1/)
 
 nshift = product(ncell) 
 natom_uc = eff_pot%crystal%natom
!Fill the typat_order array:
!In the fit script the atom must be in the order 11111 222222 33333 ..
!and the order of the atom can not be change in the fit script,
!we transform into the format of the script
 ib = 1
 ib1= 1
 do ii=1,eff_pot%crystal%ntypat
   jj = A(ii)
   do kk=1,natom_uc
     if(supercell%typat_supercell(kk)==jj)then
       typat_order_uc(ib1) = kk
       ib1 = ib1 + 1 
       do ll=1,nshift
         ia = (ll-1)*natom_uc + kk
         typat_order(ib) = ia
         ib = ib + 1
       end do
     end if
   end do
 end do

! BORN CHARGES FILE
 if (open_file('system/Born_Charges',msg,unit=unit_born,form="formatted",&
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
 if (open_file('system/Dielectric_Tensor',msg,unit=unit_epsiloninf,form="formatted",&
&    status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 do mu=1,3
   WRITE(unit_epsiloninf,'(3(F23.14))') eff_pot%harmonics_terms%epsilon_inf(:,mu)
 end do


!REFERENCE STRUCTURE FILE
 if (open_file('system/Reference_structure',msg,unit=unit_ref,form="formatted",&
&    status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write(unit_ref,'("Energy (Hartree)")')
 write(unit_ref,'("================")')
 write(unit_ref,'(F23.14)') (hist%etot(1)/nshift)
 write(unit_ref,'("")')
 write(unit_ref,'("Cell vectors")')
 write(unit_ref,'("============")')
 write(unit_ref,'(3(F23.14))') (hist%rprimd(:,:,1))
 write(unit_ref,'("")')
 write(unit_ref,'("Atomic positions (Bohr radius)")')
 write(unit_ref,'("==============================")')

 call xred2xcart(supercell%natom_supercell,supercell%rprimd_supercell,&
&                  xcart,hist%xred(:,:,1))

 do ia=1,supercell%natom_supercell
   write(unit_ref,'(3(F23.14))') xcart(:,typat_order(ia))
 end do

!Harmonic XML file
 if (open_file('system/harmonic.xml',msg,unit=unit_harmonic,form="formatted",&
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
&                 (eff_pot%harmonics_terms%ifcs%short_atmfrc(1,mu,typat_order_uc(ia),&
&                                                              nu,typat_order_uc(ib),irpt))
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
 if (open_file('system/Strain_Tensor',msg,unit=unit_strain,form="formatted",&
&     status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write(unit_strain,'(6(F23.14))') (eff_pot%harmonics_terms%elastic_constants)

! SYM FILE
 if (open_file('system/symmetry_operations',msg,unit=unit_sym,form="formatted",&
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
 if (open_file('system/Molecular_dynamic',msg,unit=unit_md,form="formatted",&
&     status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 do ii=2,nstep
   write(unit_md,'(I5)') ii-2
   write(unit_md,'(F22.14)') hist%etot(ii)/nshift
   write(unit_md,'(3(F22.14))') (hist%rprimd(:,:,ii))

!  Set xcart and fcart for this step
   call xred2xcart(supercell%natom_supercell,supercell%rprimd_supercell,&
&                  xcart,hist%xred(:,:,ii))

   fcart(:,:) = hist%fcart(:,:,ii)

   do ia=1,supercell%natom_supercell
     write(unit_md,'(3(E22.14),3(E22.14))') xcart(:,typat_order(ia)),fcart(:,typat_order(ia))
   end do
   write(unit_md,'(6(E22.14))') hist%strten(:,ii)
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
 ABI_DEALLOCATE(typat_order)
 ABI_DEALLOCATE(typat_order_uc)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(fcart)
 call destroy_supercell(supercell)

end subroutine fit_polynomial_printSystemFiles
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
     exit
   end if
 end do

end function getCoeffFromList
!!***

end module m_fit_polynomial_coeff
!!***
