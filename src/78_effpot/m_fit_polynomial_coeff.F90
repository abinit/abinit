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
 use m_phonon_supercell
 use m_effective_potential, only :  effective_potential_type
 use m_io_tools,   only : open_file
 use m_abihist, only : abihist

 implicit none

 public :: fit_polynomial_coeff_free
 public :: fit_polynomial_coeff_init
 public :: fit_polynomial_coeff_get
 public :: fit_polynomial_getNcoeffNextOrder
 public :: fit_polynomial_getNextOrder
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

subroutine fit_polynomial_coeff_get(eff_pot,option)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_get'
 use interfaces_14_hidewrite
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option
!arrays
 type(effective_potential_type), intent(in) :: eff_pot
! type(fit_polynomial_coeff_type), intent(inout) :: polynomial_coeff
!Local variables-------------------------------
!scalar
 integer :: found,ia,ib,icoeff,idisy1,idisy2,ii
 integer :: ipesy1,ipesy2,isym,iterm,itypat
 integer :: jj,lim1,lim2,lim3,mu,natom
 integer :: ncoeff1,ncoeff2,ncoeff3,ncoeff4,nu,ndisp
 integer :: npower,power,nsym,nterm
 real(dp):: tolsym8,weight
!arrays
 integer :: ncell(3),sym1(3,3),sym2(3,3)
 integer,allocatable :: all_powers(:),atindx(:,:),cell(:,:,:)
 integer,allocatable :: dir_int(:),nterms(:),powers(:)
 integer,allocatable :: pertsy(:,:),indsym(:,:,:) ,sym_coeff(:,:,:,:),sym_coeff_tmp(:,:,:,:)
 character(len=1) :: dir_char(3),powerchar
 character(len=1) :: mutodir(3) = (/"x","y","z"/)
 character(len=5),allocatable :: symbols(:),srting_sym(:,:)
 character(len=40):: name,text
 character(len=fnlen) :: filename
 character(len=500) :: message
 character(len=80) :: tmpstring
 type(polynomial_coeff_type),dimension(:),allocatable :: coeffs1
 type(polynomial_coeff_type),dimension(:),allocatable :: coeffs2,coeffs3,coeffs4
 type(polynomial_term_type),dimension(:),allocatable :: terms
 type(supercell_type) :: supercell

! *************************************************************************

!set variables
 natom = eff_pot%crystal%natom
 nsym = eff_pot%crystal%nsym

!Allocation of the symbol array of the atoms
 ABI_ALLOCATE(symbols,(natom))
 do ia=1,natom
   symbols(ia) = adjustl(znucl2symbol(eff_pot%crystal%znucl(eff_pot%crystal%typat(ia))))
 end do

!check the atoms of the same type and add numerotation
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

!Get the symetries
 ABI_ALLOCATE(srting_sym,(3,nsym))
 srting_sym = ""
 do isym = 1, eff_pot%crystal%nsym
   call  symrel2string(eff_pot%crystal%symrel(:,:,isym), eff_pot%crystal%tnons(:,isym), tmpstring)
   srting_sym(1,isym) = tmpstring(1:2)
   srting_sym(2,isym) = tmpstring(4:5)
   srting_sym(3,isym) = tmpstring(7:8)
!  Check if the values are correct
   if (.not.(any((/'+x','+y','+z'/) ==  srting_sym(1,isym)).or.&
&           (any((/'-x','-y','-z'/) ==  srting_sym(1,isym)))).and.&
&           (any((/'+x','+y','+z'/) ==  srting_sym(2,isym)).or.&
&           (any((/'-x','-y','-z'/) ==  srting_sym(2,isym)))).and.&
&           (any((/'+x','+y','+z'/) ==  srting_sym(3,isym)).or.&
&           (any((/'-x','-y','-z'/) ==  srting_sym(3,isym)))))then
     write(message, '(a,I0,9a)' )&
&          ' There is a problem with the symmetry ',isym,ch10,&
&          ' The values of the sym are:',srting_sym(:,isym),ch10,&
&          ' Action: Contact abinit group'
     MSG_ERROR(message)
   end if
 end do
 
!Init supercell

 ncell = (/2,2,2/)
 lim1=((ncell(1)/2)+1)
 lim2=((ncell(2)/2)+1)
 lim3=((ncell(3)/2)+1)

 call init_supercell(natom, 0,real(ncell,dp),&
&                    eff_pot%crystal%rprimd,eff_pot%crystal%typat,&
&                    eff_pot%crystal%xcart,supercell)

!Obtain a list of rotated atom labels:
 ABI_ALLOCATE(indsym,(4,nsym,natom))
 tolsym8=tol8
 call symatm(indsym,natom,nsym,eff_pot%crystal%symrec,eff_pot%crystal%tnons,&
&            tolsym8,eff_pot%crystal%typat,eff_pot%crystal%xred)


!Get The number of coefficient with the symmetries
 ABI_ALLOCATE(sym_coeff,(3,natom,3,natom))
 ABI_ALLOCATE(pertsy,(3,natom))
!Zero pertsy 
 sym_coeff = zero
 pertsy = zero
 do ia=1,natom
   do mu=1,3
     do ib=ia,natom
       do nu=1,3     
         do isym=1,nsym
!          reset flag
           found = 1
!          Fill sym1
           ipesy1=indsym(4,isym,ia)
           ipesy2=indsym(4,isym,ib)
           sym1(:,:)=eff_pot%crystal%symrec(:,:,isym)
           sym2(:,:)=eff_pot%crystal%symrec(:,:,isym)
!          Now that a symmetric perturbation has been obtained,
!          including the expression of the symmetry matrix, see
!          if the symmetric perturbations are available
           do idisy1=1,3
             do idisy2=1,3
               if(sym1(mu,idisy1)/=0.and.sym2(nu,idisy2)/=0)then
                 if(sym_coeff(idisy1,ipesy1,idisy2,ipesy2) == 0)then
                   found=0
                   exit
                 end if
               end if
             end do
           end do
!          Now, if still found, then it is a symmetric
!          of some linear combination of existing perturbations
           if(found==1)then
             sym_coeff(mu,ia,nu,ib)=-1
             exit ! Exit loop on symmetry operations
           end if
         end do
!        Now that all symmetries have been examined,
!        if still not symmetric of a linear combination
!        of basis perturbations, then it is a basis perturbation
!        if(pertsy(mu,ia)/=-1) pertsy(mu,ia)=1
         if(sym_coeff(mu,ia,nu,ib)/=-1) then
           sym_coeff(mu,ia,nu,ib)=1
         end if
       end do
     end do
   end do
 end do
 
!Count number of coefficients and max term (without symetries)
 nterm   = zero
 ncoeff1 = zero
 do ia=1,natom
   do mu=1,3
     do ib=ia,natom
       do nu=1,3
         if(ia/=ib.and.mu==nu)then
           nterm = nterm + 1
           if(sym_coeff(mu,ia,nu,ib)==1) then
             ncoeff1 = ncoeff1 + 1
           end if
         end if
       end do
     end do
   end do
 end do


!Count the number of term by coefficent (with the symetries)
 ABI_ALLOCATE(nterms,(ncoeff1))
 ABI_ALLOCATE(sym_coeff_tmp,(3,natom,3,natom))
 sym_coeff_tmp(:,:,:,:) = sym_coeff(:,:,:,:) ! Temporary
 nterms = zero
 icoeff = zero
 do ia=1,natom
   do mu=1,3
     do ib=1,natom
       do nu=1,3
         if(ia/=ib .and. mu==nu)then
           if(sym_coeff_tmp(mu,ia,nu,ib)==1)then
             icoeff = icoeff +1 
             nterms(icoeff) = one
!           else
             do isym=1,nsym
!              Reset flag
               found = 1
!              Fill sym1 and sym2
               ipesy1=indsym(4,isym,ia)
               ipesy2=indsym(4,isym,ib)             
               do ii=1,3
                 do jj=1,3
                   sym1(ii,jj)=eff_pot%crystal%symrec(ii,jj,isym)
                   sym2(ii,jj)=eff_pot%crystal%symrec(ii,jj,isym)
                 end do
               end do
!              Now that a symmetric perturbation has been obtained,
!              including the expression of the symmetry matrix, see
!              if the symmetric perturbations are available
               do idisy1=1,3
                 do idisy2=1,3
                   if(sym1(mu,idisy1)/=0.and.sym2(nu,idisy2)/=0)then
                     if(sym_coeff_tmp(idisy1,ipesy1,idisy2,ipesy2) == -1.and.&
&                       idisy1==idisy2.and.(ipesy1/=ipesy2)) then
                       sym_coeff_tmp(idisy1,ipesy1,idisy2,ipesy2) = 0
                       nterms(icoeff) = nterms(icoeff) +1
                       found=1
                       exit
                     end if
                   end if
                 end do
               end do
!               if(found==1)then
!                 exit ! Exit loop on symmetry operations
!               end if
             end do
           end if
         end if
       end do
     end do
   end do
 end do
 ABI_DEALLOCATE(sym_coeff_tmp)

 write(message,'(a)')," Irreductible coefficient and associated atom 1, atom 2 and direction:"
 call wrtout(std_out,message,'COLL') 
 
!Construct the coeffcients
!Allcation of the arrays
 ABI_DATATYPE_ALLOCATE(coeffs1,(ncoeff1))
 ABI_ALLOCATE(all_powers,(power))
 all_powers(:)=(/3/)
 
 ndisp = 1
 ABI_ALLOCATE(atindx,(2,ndisp))
 ABI_ALLOCATE(cell,(3,2,ndisp))
 ABI_ALLOCATE(dir_int,(ndisp))
 ABI_ALLOCATE(powers,(ndisp))

 npower  = one
 icoeff  = zero
 iterm   = zero
 do ia=1,natom
   do mu=1,3
     do ib=1,natom
       do nu=1,3
         if(ia/=ib.and.mu==nu)then
           if(sym_coeff(mu,ia,nu,ib)==1)then
             icoeff = icoeff + 1 
             iterm = 1
             ABI_DATATYPE_ALLOCATE(terms,(nterms(icoeff)))

             name = ""
             atindx(1,:) = ia; atindx(2,:) = ib;
             cell(:,:,:) = zero
             dir_char(:) = mutodir(mu)
             dir_int(:)  = mu
             ndisp  = 1 
             weight = one           
             powers(:) = one
             write(powerchar,'(I0)') 1
             call polynomial_coeff_getName(text,atm1=symbols(ia),atm2=symbols(ib),&
&                                        dir=dir_char(1),power=trim(powerchar))
             name = trim(name)//trim(text)
             write(message,'(2a)'),' ',trim(name)
             call wrtout(std_out,message,'COLL') 
             call polynomial_term_init(atindx,cell,dir_int,ndisp,terms(1),powers,weight)
             write(message,'(a,I0,a,I0,2a)'),'    Atom ',ia,' and atom ',ib&
&                                             ," in the direction ",dir_char(1)
             call wrtout(std_out,message,'COLL')
!           else
             do isym=1,nsym
!              reset flag
               found = 1

!              Fill sym1
               ipesy1=indsym(4,isym,ia)
               ipesy2=indsym(4,isym,ib)
             
               sym1(:,:)=eff_pot%crystal%symrec(:,:,isym)
               sym2(:,:)=eff_pot%crystal%symrec(:,:,isym)

!              Now that a symmetric perturbation has been obtained,
!              including the expression of the symmetry matrix, see
!              if the symmetric perturbations are available
               do idisy1=1,3
                 do idisy2=1,3
                   if(sym1(mu,idisy1)/=0.and.sym2(nu,idisy2)/=0)then
                     if(sym_coeff(idisy1,ipesy1,idisy2,ipesy2) == -1.and.&
&                       idisy1==idisy2.and.(ipesy1/=ipesy2)) then
                       sym_coeff(idisy1,ipesy1,idisy2,ipesy2) = 0
                       iterm = iterm + 1
                       atindx(1,:) = ipesy1; atindx(2,:) = ipesy2; cell(:,:,:) = zero
                       dir_char(:) = mutodir(idisy1); dir_int(:)  = idisy1
                       ndisp  = one; weight = one           
                       powers(:) = one
                       write(powerchar,'(I0)') 1
                       call polynomial_term_init(atindx,cell,dir_int,ndisp,terms(iterm),powers,weight)
                       write(message,'(a,I0,a,I0,3a)'),'    Atom ',ipesy1,' and atom ',ipesy2&
&                                                     ," in the direction ",dir_char(1)
                       call wrtout(std_out,message,'COLL')
                       found=0
                       exit
                     end if
                   end if
                 end do
               end do
!              Now, if still found, then it is a symmetric
!              of some linear combination of existing perturbations
               if(found==1)then
 !                exit ! Exit loop on symmetry operations
               end if
             end do
             if(iterm == nterms(icoeff))then
!            Initialisatoin of the coefficients
               call polynomial_coeff_init(zero,name,nterms(icoeff),coeffs1(icoeff),terms)
               do iterm=1,nterms(icoeff)
                 call polynomial_term_free(terms(iterm))
               end do
               iterm = zero
               ABI_DEALLOCATE(terms)
             else
               write (message, '(a,I0,a,I0,a,I0,3a)')&
&                     'Probleme with the total number of term for the coefficient number ',icoeff,&
&                     ch10,iterm,' are computed but ',nterms(icoeff),' are needed',ch10,&
&                     'action: contact abinit groupe'
               MSG_BUG(message)
             end if
           end if
         end if
       end do
     end do
   end do
 end do




  filename = "test1.xml"
  call polynomial_coeff_writeXML(coeffs1,ncoeff1,filename=filename)
  write(message,'(a,I0,a)'),&
 &       ' with ',ncoeff1,' for the 1st order '
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL') 

!Compute higher orders
 call fit_polynomial_getNcoeffNextOrder(coeffs1,coeffs1,ncoeff1,ncoeff1,ncoeff2)
 ABI_DATATYPE_ALLOCATE(coeffs2,(ncoeff2))
 call fit_polynomial_getNextOrder(coeffs1,coeffs1,coeffs2,natom,ncoeff1,ncoeff1,ncoeff2,symbols)
 filename = "test2.xml"
 call polynomial_coeff_writeXML(coeffs2,ncoeff2,filename=filename)
 write(message,'(a,I0,a)'),&
&       ' with ',ncoeff2,' for the 2nd order '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL') 

 call fit_polynomial_getNcoeffNextOrder(coeffs2,coeffs1,ncoeff2,ncoeff1,ncoeff3)
 ABI_DATATYPE_ALLOCATE(coeffs3,(ncoeff3))
 call fit_polynomial_getNextOrder(coeffs2,coeffs1,coeffs3,natom,ncoeff2,ncoeff1,ncoeff3,symbols)
 filename = "test3.xml"
 call polynomial_coeff_writeXML(coeffs3,ncoeff3,filename=filename)
 write(message,'(a,I0,a)'),&
&       ' with ',ncoeff3,' for the 3rd order '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL') 

 call fit_polynomial_getNcoeffNextOrder(coeffs3,coeffs1,ncoeff3,ncoeff1,ncoeff4)
 ABI_DATATYPE_ALLOCATE(coeffs4,(ncoeff4))
 call fit_polynomial_getNextOrder(coeffs3,coeffs1,coeffs4,natom,ncoeff3,ncoeff1,ncoeff4,symbols)
 filename = "test4.xml"
 call polynomial_coeff_writeXML(coeffs4,ncoeff4,filename=filename)
 write(message,'(a,I0,a)'),&
&       ' with ',ncoeff4,' for the 4th order '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL') 


!Deallocation 
 ABI_DEALLOCATE(all_powers)
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(cell)
 ABI_DEALLOCATE(dir_int)
 ABI_DEALLOCATE(nterms)
 ABI_DEALLOCATE(pertsy)
 ABI_DEALLOCATE(powers)
 ABI_DEALLOCATE(sym_coeff)
 ABI_DEALLOCATE(symbols)
 ABI_DEALLOCATE(srting_sym)
 ABI_DEALLOCATE(indsym)
 do ii=1,ncoeff1
   call polynomial_coeff_free(coeffs1(ii))
 end do
 ABI_DEALLOCATE(coeffs1)
 do ii=1,ncoeff2
   call polynomial_coeff_free(coeffs2(ii))
 end do
 ABI_DEALLOCATE(coeffs2)
 do ii=1,ncoeff3
   call polynomial_coeff_free(coeffs3(ii))
 end do
 ABI_DEALLOCATE(coeffs3)
 do ii=1,ncoeff4
   call polynomial_coeff_free(coeffs4(ii))
 end do
 ABI_DEALLOCATE(coeffs4)
 
 
 
 call destroy_supercell(supercell)


end subroutine fit_polynomial_coeff_get
!!***


!!****f* m_fit_polynomial_coeff/fit_polynomial_getNextOrder
!!
!! NAME
!! fit_polynomial_getNextOrder
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

subroutine fit_polynomial_getNextOrder(coeffs1,coeffs2,coeffs3,natom,ncoeff1,ncoeff2,ncoeff3,symbols)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_getNextOrder'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ncoeff1,ncoeff2,ncoeff3
!arrays
 character(len=5) :: symbols(natom)
 type(polynomial_coeff_type),intent(in) :: coeffs1(ncoeff1),coeffs2(ncoeff2)
 type(polynomial_coeff_type),intent(inout) :: coeffs3(:)
!Local variables-------------------------------
!scalar
 integer :: ii,icoeff1,icoeff2,icoeff3
 integer :: idisp1,idisp2,idisp3,iterm1,iterm2,iterm3
 integer :: ndisp
 integer :: nterm
 real(dp):: weight
 logical :: compatible,found
!arrays
 integer,allocatable :: atindx(:,:),cell(:,:,:),dir_int(:),powers(:)
 character(len=1) :: powerchar
 character(len=1) :: mutodir(3) = (/"x","y","z"/)
 character(len=40):: name,text
 character(len=500) :: message
 type(polynomial_term_type),dimension(:),allocatable :: terms

! *************************************************************************

 icoeff3 = one
 do icoeff1=1,ncoeff1
   do icoeff2=1,ncoeff2

     compatible = .true.

     if (compatible)then
!      Get the number of term for the new coefficient
!      and allocate the new array terms
       nterm = coeffs1(icoeff1)%nterm*coeffs2(icoeff2)%nterm
       ABI_DATATYPE_ALLOCATE(terms,(nterm))
!      Loops over the terms
       iterm3 = 1
       do iterm1=1,coeffs1(icoeff1)%nterm
         do iterm2=1,coeffs2(icoeff2)%nterm

!          Count the number of displacement
           ndisp = 1
           do idisp1=1,coeffs1(icoeff1)%terms(iterm1)%ndisp
             do idisp2=1,coeffs2(icoeff2)%terms(iterm2)%ndisp
               if (coeffs1(icoeff1)%terms(iterm1)%atindx(1,idisp1)/=&
&                  coeffs2(icoeff2)%terms(iterm2)%atindx(1,idisp2).or.&
&                  coeffs1(icoeff1)%terms(iterm1)%atindx(2,idisp1)/=&
&                  coeffs2(icoeff2)%terms(iterm2)%atindx(2,idisp2).or.&
&                  coeffs1(icoeff1)%terms(iterm1)%direction(idisp1)/=&
&                  coeffs2(icoeff2)%terms(iterm2)%direction(idisp2))  then
                 ndisp = ndisp + 1
               end if
             end do
           end do

!          Allocation of the new array for this term
           ABI_ALLOCATE(atindx,(2,ndisp))
           ABI_ALLOCATE(cell,(3,2,ndisp))
           ABI_ALLOCATE(dir_int,(ndisp))
           ABI_ALLOCATE(powers,(ndisp))

           idisp3 = 1
!          1-copy the displacement from the first coefficient
           do idisp1=1,coeffs1(icoeff1)%terms(iterm1)%ndisp
             atindx(:,idisp3) = coeffs1(icoeff1)%terms(iterm1)%atindx(:,idisp1)
             cell(:,:,idisp3) = coeffs1(icoeff1)%terms(iterm1)%cell(:,:,idisp1)
             dir_int( idisp3) = coeffs1(icoeff1)%terms(iterm1)%direction(idisp1)
             powers( idisp3)  = coeffs1(icoeff1)%terms(iterm1)%power(idisp1)
             weight           = coeffs1(icoeff1)%terms(iterm1)%weight
             idisp3 = idisp3 + 1 
           end do

!         In this case, copy  displacement from the second coefficient             
           do idisp2=1,coeffs2(icoeff2)%terms(iterm2)%ndisp
!            Test if all the cofficients are different
             if (ndisp == (coeffs1(icoeff1)%terms(iterm1)%ndisp+&
&                          coeffs2(icoeff2)%terms(iterm2)%ndisp))then
               atindx(:,idisp3) = coeffs2(icoeff2)%terms(iterm2)%atindx(:,idisp2)
               cell(:,:,idisp3) = coeffs2(icoeff2)%terms(iterm2)%cell(:,:,idisp2)
               dir_int( idisp3) = coeffs2(icoeff2)%terms(iterm2)%direction(idisp2)
               powers( idisp3)  = coeffs2(icoeff2)%terms(iterm2)%power(idisp2)
               weight           = coeffs2(icoeff2)%terms(iterm2)%weight
               idisp3 = idisp3 + 1 
             else
!            There is identical coefficients
               found = .false.
               do idisp1=1,coeffs1(icoeff1)%terms(iterm1)%ndisp
                 if (coeffs1(icoeff1)%terms(iterm1)%atindx(1,idisp1)==&
&                    coeffs2(icoeff2)%terms(iterm2)%atindx(1,idisp2).and.&
&                    coeffs1(icoeff1)%terms(iterm1)%atindx(2,idisp1)==&
&                    coeffs2(icoeff2)%terms(iterm2)%atindx(2,idisp2).and.&
&                    coeffs1(icoeff1)%terms(iterm1)%direction(idisp1)==&
&                    coeffs2(icoeff2)%terms(iterm2)%direction(idisp2))  then                  
                   powers(idisp1)  = powers(idisp1)+1                   
                   found = .true.
                 end if
               end do
               if(.not.found)then
                 atindx(:,idisp3) = coeffs2(icoeff2)%terms(iterm2)%atindx(:,idisp2)
                 cell(:,:,idisp3) = coeffs2(icoeff2)%terms(iterm2)%cell(:,:,idisp2)
                 dir_int( idisp3) = coeffs2(icoeff2)%terms(iterm2)%direction(idisp2)
                 powers( idisp3)  = coeffs2(icoeff2)%terms(iterm2)%power(idisp2)
                 weight           = coeffs2(icoeff2)%terms(iterm2)%weight
                 idisp3 = idisp3 + 1 
               end if
             end if
           end do

           if (ndisp /= idisp3-1) then
             write(message, '(3a)' )&
&           ' The number of displacement does not correspond. ',ch10,&
&           ' Action: Contact abinit group'
             MSG_ERROR(message)
           end if

           call polynomial_term_init(atindx,cell,dir_int,ndisp,terms(iterm3),powers,weight)
!          Deallocation of the  array
           ABI_DEALLOCATE(atindx)
           ABI_DEALLOCATE(cell)
           ABI_DEALLOCATE(dir_int)
           ABI_DEALLOCATE(powers)

           iterm3 = iterm3 + 1
         end do
       end do
       name=''

       do ii=1,terms(1)%ndisp
         write(powerchar,'(I0)') terms(1)%power(ii)
         call polynomial_coeff_getName(text,&
&                                      atm1=symbols(terms(1)%atindx(1,ii)),&
&                                      atm2=symbols(terms(1)%atindx(2,ii)),&
&                                      dir=mutodir(terms(1)%direction(ii)),&
&                                      power=trim(powerchar))
         name = trim(name)//trim(text)
       end do
       call polynomial_coeff_init(zero,name,nterm,coeffs3(icoeff3),terms)
       icoeff3 = icoeff3 + 1

!      Deallocate the array terms
       do ii=1,nterm
         call polynomial_term_free(terms(ii))
       end do
       ABI_DEALLOCATE(terms)
     else
       cycle
     end if
   end do
 end do

end subroutine fit_polynomial_getNextOrder

!!****f* m_fit_polynomial_coeff/fit_polynomial_getNcoeffNextOrder
!!
!! NAME
!! fit_polynomial_getNcoeffNextOrder
!!
!! FUNCTION
!! Return the number of coefficients for the next order 
!! of the polynome
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

subroutine fit_polynomial_getNcoeffNextOrder(coeffs1,coeffs2,ncoeff1,ncoeff2,ncoeff3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_getNcoeffNextOrder'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncoeff1,ncoeff2
 integer,intent(out):: ncoeff3
!arrays
 type(polynomial_coeff_type),intent(in) :: coeffs1(ncoeff1),coeffs2(ncoeff2)
!Local variables-------------------------------
!scalar
 integer :: icoeff1,icoeff2
!arrays
! *************************************************************************

 ncoeff3 = zero

 do icoeff1=1,ncoeff1
   do icoeff2=1,ncoeff2
     ncoeff3 = ncoeff3 + 1
   end do
 end do

end subroutine fit_polynomial_getNcoeffNextOrder

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
 use interfaces_14_hidewrite
 use interfaces_32_util
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

















