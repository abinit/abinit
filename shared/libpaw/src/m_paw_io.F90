!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_io
!! NAME
!!  m_paw_io
!!
!! FUNCTION
!!  PAW I/O related operations
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2019 ABINIT group (MT, TR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

module m_paw_io

 USE_DEFS
 USE_MSG_HANDLING
 USE_MEMORY_PROFILING

 implicit none

 private

 public :: pawio_print_ij
!!***

CONTAINS !===========================================================

!!****f* m_paw_io/pawio_print_ij
!! NAME
!! pawio_print_ij
!!
!! FUNCTION
!! Print ij_ square matrixes in a "suitable" format.
!! Data are "energy-like" in Hartree units.
!! Devoted to the printing of rhoij, Dij -like PAW matrixes.
!!
!! INPUTS
!!  a_ij(cplex*adim)= input square matrix
!!  asym_ij(cplex*adim)= -OPTIONAL ARGUMENT-
!!                       When present, A(j,i) is deduced from asym_ij
!!                                     instead of a_ij
!!  adim= dimension of array a_ij:
!!        adim=ndim*(ndim+1)/2                   if opt_pack= 0
!!        adim=number of non-zero values of a_ij if opt_pack=+1
!!  cplex=1 if a_ij is real, 2 if it is complex
!!  [mode_paral]= --optional argument, default='COLL'--
!!   'COLL' if all procs are calling the routine with the same message to be written once only.
!!   'PERS' if the procs are calling the routine with different messages each to be written,
!!          or if one proc is calling the routine
!!  ndim= dimension of input square matrix
!!  opt_l= if <0  all parts of a_ij are printed
!!         if >=0 only parts of a_ij corresponding to li=lj=opt_l are printed
!!  opt_l_index(ndim)= array giving l quantum number for each 1<=ilmn<=ndim
!!                     not used if opt_l<0
!!  opt_pack= 0 if a_ij is given as A(j(j-1)/2+i), i<=j
!!           +1 if a_ij is given as A(j(j-1)/2+i) and is in "packed storage"
!!                                  (i.e. only non-zero values are stored)
!!  opt_prtvol= >=0 if up to 12 components of _ij matrix have to be printed
!!               <0 if all components of ij_ matrix have to be printed
!!  opt_sym= -OPTIONAL ARGUMENT- (default if not present: opt_sym=2)
!!          Define the symmetry of a_ij matrix:
!!            opt_sym=1 : A(j,i)= A(i,j)
!!            opt_sym=2 : A(j,i)= Conjg[A(i,j)]
!!            opt_sym=3 : A(j,i)=-A(i,j)
!!            opt_sym=4 : A(j,i)=-Conjg[A(i,j)]
!!            When asym_ij argument is present, A[i,j] is taken from it.
!!  pack2ij(adim)= gives the (i,j) index of of packed value of rhoij
!!                 used only if opt_packed=+1
!!  test_value= (real number) if positive, print a warning when the
!!              magnitude of a_ij is greater than opt_test
!!              No test when test_value<0
!!  unit=the unit number for output
!!  Ha_or_eV= 1: output in hartrees, 2: output in eV
!!
!! PARENTS
!!      m_paw_ij,m_paw_slater,m_pawdij,m_pawrhoij,pawmkrhoij,pawprt
!!      setrhoijpbe0,wfd_pawrhoij
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine pawio_print_ij(unit,a_ij,adim,cplex,ndim,opt_l,opt_l_index, &
&                         opt_pack,opt_prtvol,pack2ij,test_value,Ha_or_eV, &
&                         mode_paral,opt_sym,asym_ij) ! Optional arguments

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: adim,cplex,ndim,opt_l,opt_pack,opt_prtvol,unit,Ha_or_eV
 integer,intent(in),optional :: opt_sym
 character(len=*),optional,intent(in) :: mode_paral
 real(dp),intent(in) :: test_value
!arrays
 integer,intent(in) :: opt_l_index(ndim*min(1+opt_l,1)),pack2ij(adim*opt_pack)
 real(dp),intent(in) :: a_ij(cplex*adim)
 real(dp),intent(in),optional :: asym_ij(cplex*adim)

!Local variables ---------------------------------------
! Adjust format bellow according to maxprt
!scalars
 integer,parameter :: maxprt_default=12
 integer :: dplex,ilmn,ilmn1,j0lmn,jlmn,jlmn1,klmn,klmn1,klmn2,maxprt,nhigh
 integer :: nmin,optsym
 real(dp) :: testval
 logical :: use_asym
 character(len=4) :: mode_paral_
 character(len=500) :: msg=''
!arrays
 real(dp),parameter :: fact_re(4)=(/one,one,-one,-one/),fact_im(4)=(/one,-one,-one,one/)
 real(dp) :: tabmax(cplex),tabmin(cplex)
 real(dp),allocatable :: b_ij(:),bsym_ij(:),prtab(:,:,:)

! *************************************************************************

 10 format(100(1x,f9.5))
 11 format(12(1x,f9.5),a) !Change this format according to variable "maxprt"


!DEBUG
!write(std_out,*)' pawio_print_ij : enter '
!ENDDEBUG

!Optional arguments
 mode_paral_='COLL';if (present(mode_paral)) mode_paral_=mode_paral
 use_asym=present(asym_ij)
 if (present(opt_sym)) then
   optsym=opt_sym
 else
   optsym=2
 end if

!Define size of square matrix
 if (opt_prtvol>=0) then
   maxprt=maxprt_default
 else
   maxprt=ndim
 end if
 nmin=min(ndim,maxprt)

 if (opt_l>=0) nmin=count(opt_l_index(:)==opt_l)
 LIBPAW_ALLOCATE(prtab,(cplex,nmin,nmin))
 dplex=cplex-1

!Eventually unpack input matrix(es)
 LIBPAW_ALLOCATE(b_ij,(cplex*ndim*(ndim+1)/2))
 if (opt_pack==0) then
   b_ij=a_ij
 else if (opt_pack==1) then
   b_ij=zero
   do klmn=1,adim
     klmn1=cplex*klmn-dplex
     klmn2=cplex*pack2ij(klmn)-dplex
     b_ij(klmn2:klmn2+dplex)=a_ij(klmn1:klmn1+dplex)
   end do
 end if
 if (opt_prtvol<0.and.opt_l<0) then
   if (cplex==1) then
     tabmax(1)=maxval(abs(b_ij))
     tabmin(1)=minval(abs(b_ij))
   else
     tabmax(1:2)=zero;tabmin(1:2)=1.d20
     do klmn=1,size(b_ij)/cplex
       klmn2=2*klmn
       tabmax(1)=max(tabmax(1),b_ij(klmn2-1))
       tabmin(1)=min(tabmin(1),b_ij(klmn2-1))
       tabmax(2)=max(tabmax(2),b_ij(klmn2  ))
       tabmin(2)=min(tabmin(2),b_ij(klmn2  ))
     end do
   end if
 end if
 if (use_asym) then
   LIBPAW_ALLOCATE(bsym_ij,(cplex*ndim*(ndim+1)/2))
   if (opt_pack==0) then
     bsym_ij=asym_ij
   else if (opt_pack==1) then
     bsym_ij=zero
     do klmn=1,adim
       klmn1=cplex*klmn-dplex
       klmn2=cplex*pack2ij(klmn)-dplex
       bsym_ij(klmn2:klmn2+dplex)=asym_ij(klmn1:klmn1+dplex)
     end do
   end if
   if (opt_prtvol<0.and.opt_l<0) then
     if (cplex==1) then
       tabmax(1)=max(tabmax(1),maxval(abs(bsym_ij)))
       tabmin(1)=min(tabmin(1),minval(abs(bsym_ij)))
     else
       do klmn=1,ndim
         klmn2=2*klmn
         tabmax(1)=max(tabmax(1),bsym_ij(klmn2-1))
         tabmin(1)=min(tabmin(1),bsym_ij(klmn2-1))
         tabmax(2)=max(tabmax(2),bsym_ij(klmn2  ))
         tabmin(2)=min(tabmin(2),bsym_ij(klmn2  ))
       end do
     end if
   end if
 end if

!Transfer triangular matrix to rectangular one
 jlmn1=0
 do jlmn=1,ndim
   if (opt_l<0) then
     jlmn1=jlmn;if (jlmn1>nmin) cycle
   else if (opt_l_index(jlmn)==opt_l) then
     jlmn1=jlmn1+1
   else
     cycle
   end if
   ilmn1=0;j0lmn=jlmn*(jlmn-1)/2
   do ilmn=1,jlmn
     if (opt_l<0) then
       ilmn1=ilmn
     else if (opt_l_index(ilmn)==opt_l) then
       ilmn1=ilmn1+1
     else
       cycle
     end if
     klmn=j0lmn+ilmn
     if (cplex==1) then
       prtab(1,ilmn1,jlmn1)=b_ij(klmn)
       if (use_asym) then
         prtab(1,jlmn1,ilmn1)=fact_re(optsym)*bsym_ij(klmn)
       else
         prtab(1,jlmn1,ilmn1)=fact_re(optsym)*b_ij(klmn)
       end if
     else
       klmn=2*klmn
       prtab(1:2,ilmn1,jlmn1)=b_ij(klmn-1:klmn)
       if (use_asym) then
         prtab(1,jlmn1,ilmn1)=fact_re(optsym)*bsym_ij(klmn-1)
         prtab(2,jlmn1,ilmn1)=fact_im(optsym)*bsym_ij(klmn  )
       else
         prtab(1,jlmn1,ilmn1)=fact_re(optsym)*b_ij(klmn-1)
         prtab(2,jlmn1,ilmn1)=fact_im(optsym)*b_ij(klmn  )
       end if
     end if
   end do
 end do
 LIBPAW_DEALLOCATE(b_ij)

 if (use_asym)  then
   LIBPAW_DEALLOCATE(bsym_ij)
 end if

 if (Ha_or_eV==2) then
   prtab=prtab*Ha_eV
   if (opt_prtvol<0.and.opt_l<0) then
     tabmax=tabmax*Ha_eV
     tabmin=tabmin*Ha_eV
   end if
 end if

 if (cplex==2) then
   write(msg,'(3x,a)') '=== REAL PART:'
   call wrtout(unit,msg,mode_paral_)
 end if

 if (ndim<=maxprt.or.opt_l>=0) then
   do ilmn=1,nmin
     write(msg,fmt=10) prtab(1,1:nmin,ilmn)
     call wrtout(unit,msg,mode_paral_)
   end do
 else
   do ilmn=1,nmin
     write(msg,fmt=11) prtab(1,1:nmin,ilmn),' ...'
     call wrtout(unit,msg,mode_paral_)
   end do
   write(msg,'(3x,a,i2,a)') '...  only ',maxprt,'  components have been written...'
   call wrtout(unit,msg,mode_paral_)
 end if
 if (opt_prtvol<0.and.opt_l<0) then
   write(msg,'(3x,2(a,es9.2))') 'max. value= ',tabmax(1),', min. value= ',tabmin(1)
   call wrtout(unit,msg,mode_paral_)
 end if

 if (cplex==2) then
   write(msg,'(3x,a)') '=== IMAGINARY PART:'
   call wrtout(unit,msg,mode_paral_)
   if (ndim<=maxprt.or.opt_l>=0) then
     do ilmn=1,nmin
       write(msg,fmt=10) prtab(2,1:nmin,ilmn)
       call wrtout(unit,msg,mode_paral_)
     end do
   else
     do ilmn=1,nmin
       write(msg,fmt=11) prtab(2,1:nmin,ilmn),' ...'
       call wrtout(unit,msg,mode_paral_)
     end do
     write(msg,'(3x,a,i2,a)') '...  only ',maxprt,'  components have been written...'
     call wrtout(unit,msg,mode_paral_)
   end if
   if (opt_prtvol<0.and.opt_l<0) then
     write(msg,'(3x,2(a,es9.2))') 'max. value= ',tabmax(2),', min. value= ',tabmin(2)
     call wrtout(unit,msg,mode_paral_)
   end if
 end if

 if (test_value>zero) then
   testval=test_value;if (Ha_or_eV==2) testval=testval*Ha_eV
   nhigh=0;nhigh=count(abs(prtab(:,:,:))>=testval)
   if (nhigh>0) then
     write(msg,'(5a,i3,a,f6.1,7a)')&
&     ' pawio_print_ij: WARNING -',ch10,&
&     '  The matrix seems to have high value(s) !',ch10,&
&     '  (',nhigh,' components have a value greater than ',testval,').',ch10,&
&     '  It can cause instabilities during SCF convergence.',ch10,&
&     '  Action: you should check your atomic dataset (psp file)',ch10,&
&     '          and look for "high" projector functions...'
     call wrtout(unit,msg,mode_paral_)
   end if
 end if

 LIBPAW_DEALLOCATE(prtab)

!DEBUG
!write(std_out,*)' pawio_print_ij : exit '
!ENDDEBUG

end subroutine pawio_print_ij
!!***

end module m_paw_io
!!***
