!!****m* ABINIT/m_oper
!! NAME
!!  m_oper
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2020 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

MODULE m_oper

 use defs_basis
 use m_abicore
 use m_errors

 use m_matlu,    only : matlu_type
 use m_hide_lapack,  only : xginv

 implicit none

 private

 public :: init_oper
 public :: diff_oper
 public :: destroy_oper
 public :: print_oper
 public :: inverse_oper
 public :: loc_oper
 public :: identity_oper
 public :: copy_oper
 public :: trace_oper
 public :: upfold_oper
 public :: prod_oper
!!***

!!****t* m_oper/oper_type
!! NAME
!!  oper_type
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! SOURCE

 type, public :: oper_type ! for each atom

!  integer :: maxlpawu         ! Number of correlated atoms
!
!  integer :: mband
!  ! Number of bands
!
!
  integer :: nkpt
  ! Number of k-point in the IBZ.
!
  integer :: natom
!
  integer :: mbandc
!  ! Total number of bands in the Kohn-Sham Basis for PAW+DMFT
!
  integer :: nspinor
!
  integer :: nsppol

  integer :: has_operks

  integer :: has_opermatlu

  character(len=12) :: whichoper
  ! describe the type of operator computed (DFT, DMFT, KS..)

!  ! Polarisation
  type(matlu_type), allocatable :: matlu(:)
!   Local projection on correlated orbitals

  complex(dpc), allocatable :: ks(:,:,:,:)
!   In the KS basis  (nsppol,nkpt,nband,nband)

  real(dp), pointer :: wtk(:) => null()

 end type oper_type
!!***

!----------------------------------------------------------------------


CONTAINS  !========================================================================================
!!***

!!****f* m_oper/init_oper
!! NAME
!! init_oper
!!
!! FUNCTION
!!  Allocate variables used in type oper_type.
!!
!! INPUTS
!!
!! OUTPUTS
!! oper  = operator of type oper_type
!!
!! PARENTS
!!      m_datafordmft,m_forctqmc,m_green,m_hubbard_one,m_outscfcv,m_self
!!      m_vtorho
!!
!! CHILDREN
!!      prod_matlu
!!
!! SOURCE

subroutine init_oper(paw_dmft,oper,nkpt,wtk,opt_ksloc)

 use defs_basis
 use m_matlu, only : init_matlu
 use m_paw_dmft, only : paw_dmft_type
 use m_errors

!Arguments ------------------------------------
!scalars
 integer, optional, intent(in) :: nkpt
 integer, optional, intent(in) :: opt_ksloc
!type
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(oper_type), intent(inout) :: oper
!arrays
 real(dp), pointer, optional :: wtk(:)
!oper variables ------------------------------------
 integer :: iatom,optksloc

!************************************************************************
 DBG_ENTER("COLL")
 if(present(opt_ksloc)) then
   optksloc=opt_ksloc
 else
   optksloc=3
 endif

 if(optksloc/=3) then
    ! FIXME: empty line!
 endif

 oper%has_operks=0
 oper%has_opermatlu=0
! ===================
!  Integers
! ===================
 oper%nsppol=paw_dmft%nsppol
 oper%nspinor=paw_dmft%nspinor
 oper%mbandc=paw_dmft%mbandc
 oper%natom=paw_dmft%natom

! ===================
!  KS variables
! ===================
 if(optksloc==1.or.optksloc==3) then
   if(.not.present(nkpt)) then
     oper%nkpt=paw_dmft%nkpt
   else
     oper%nkpt=nkpt
   endif
! allocate(oper%wtk(oper%nkpt))
   if(.not.present(wtk)) then
     oper%wtk=>paw_dmft%wtk
   else
     oper%wtk=>wtk
   endif
   ABI_ALLOCATE(oper%ks,(paw_dmft%nsppol,oper%nkpt,paw_dmft%mbandc,paw_dmft%mbandc))
   oper%has_operks=1
   oper%ks=czero
 endif

! ===================
!  matlu variables
! ===================
 if(optksloc==2.or.optksloc==3) then
   oper%has_opermatlu=0
   ABI_DATATYPE_ALLOCATE(oper%matlu,(oper%natom))
   oper%has_opermatlu=1
   call init_matlu(oper%natom,paw_dmft%nspinor,paw_dmft%nsppol,paw_dmft%lpawu,oper%matlu)
   do iatom=1,oper%natom
    oper%matlu(iatom)%mat=czero
   enddo
 endif

 DBG_EXIT("COLL")

end subroutine init_oper
!!***

!!****f* m_oper/destroy_oper
!! NAME
!! destroy_oper
!!
!! FUNCTION
!!  deallocate oper
!!
!! INPUTS
!!  oper
!!
!! OUTPUT
!!
!! PARENTS
!!      m_datafordmft,m_forctqmc,m_green,m_hubbard_one,m_outscfcv,m_self
!!      m_vtorho
!!
!! CHILDREN
!!      prod_matlu
!!
!! SOURCE

subroutine destroy_oper(oper)

 use defs_basis
 use m_crystal, only : crystal_t
 use m_matlu, only : destroy_matlu
 use m_errors

!Arguments ------------------------------------
!scalars
 type(oper_type),intent(inout) :: oper
!local variables-------------------------------
 character(len=500) :: message

!! *********************************************************************
 DBG_ENTER("COLL")
 if(oper%has_opermatlu==1) then
   call destroy_matlu(oper%matlu,oper%natom)
 else
   message = " Operator is not defined to be used in destroy_oper"
   MSG_ERROR(message)
 endif
 if ( allocated(oper%matlu))  then
   ABI_DATATYPE_DEALLOCATE(oper%matlu)
   oper%has_opermatlu=0
 endif
 if ( allocated(oper%ks)) then
   ABI_DEALLOCATE(oper%ks)
   oper%has_operks=0
 endif
 oper%wtk => null()
!  no deallocation for wtk: wtk is an explicit pointer

 DBG_EXIT("COLL")
end subroutine destroy_oper
!!***

!!****f* m_oper/copy_oper
!! NAME
!! copy_oper
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_datafordmft,m_green
!!
!! CHILDREN
!!      prod_matlu
!!
!! SOURCE

subroutine copy_oper(oper1,oper2)

 use defs_basis
 use m_matlu, only : copy_matlu
 use m_errors

!Arguments ------------------------------------
!type
 type(oper_type),intent(in) :: oper1
 type(oper_type),intent(inout) :: oper2 !vz_i

!oper variables-------------------------------
 integer ::  ib, ib1, ikpt, isppol
! *********************************************************************
 DBG_ENTER("COLL")
 if(oper1%has_opermatlu==1.and.oper2%has_opermatlu==1)  then
   call copy_matlu(oper1%matlu,oper2%matlu,oper1%natom)
 endif

 if(oper1%has_operks==1.and.oper2%has_operks==1) then
   do isppol=1,oper1%nsppol
     do ikpt=1,oper1%nkpt
       do ib=1,oper1%mbandc
         do ib1=1,oper1%mbandc
           oper2%ks(isppol,ikpt,ib,ib1)=oper1%ks(isppol,ikpt,ib,ib1)
         enddo
       enddo
     enddo
   enddo
 endif

 DBG_EXIT("COLL")
end subroutine copy_oper
!!***

!!****f* m_oper/print_oper
!! NAME
!! print_oper
!!
!! FUNCTION
!!
!! INPUTS
!!
!! option= 1
!!         2
!!         3
!!         4
!!  below  5: write diagonal part of KS occupation matrix
!!         6
!!         7
!!  above  8: write all elements of KS occup. matrix.
!!         9
!!
!! OUTPUT
!!
!! PARENTS
!!      m_green,m_self
!!
!! CHILDREN
!!      prod_matlu
!!
!! SOURCE

subroutine print_oper(oper,option,paw_dmft,prtopt)

 use defs_basis
 use m_matlu, only : print_matlu
 use m_paw_dmft, only : paw_dmft_type
 use m_errors

!Arguments ------------------------------------
!type
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(oper_type),intent(in) :: oper
 integer, intent(in) :: option,prtopt

!oper variables-------------------------------
 integer :: ib,ib1,ikpt,isppol,nkptr,iband1,iband2
 character(len=2000) :: message
 logical :: ximag
 real(dp) :: maximag(3)
! *********************************************************************
 DBG_ENTER("COLL")
 maximag(:)=zero

 if(oper%has_opermatlu==1) then
   write(message,'(2a)') ch10,'   = In the atomic basis'
   call wrtout(std_out,message,'COLL')
   call print_matlu(oper%matlu,oper%natom,prtopt)
 endif
 if(oper%has_operks==1) then
   write(message,'(2a)') ch10,'   = In the KS basis'
   call wrtout(std_out,message,'COLL')

!todo_ba complete print_out
   iband1=1
   iband2=oper%mbandc
!   do ib=1,oper%mbandc
!     if(-(paw_dmft%eigen_dft(1,1,ib)+paw_dmft%fermie).ge.0.3) iband1=ib
!     if( (paw_dmft%eigen_dft(1,1,ib)-paw_dmft%fermie).le.0.3) iband2=ib
!   enddo

   ximag=.false.
   if((abs(prtopt)>=3.and.((option<5).or.(option>8))).and.oper%has_operks==1) then
!     write(message,'(x,a,a,i4,2x,a)') ch10,'  -KS states'
!     call wrtout(std_out,message,'COLL')
     do isppol=1,paw_dmft%nsppol
       write(message, '(a,3x,a,i4)') ch10,"--isppol--",isppol
       call wrtout(std_out,message,'COLL')
       nkptr=oper%nkpt
       nkptr=min(oper%nkpt,4)
       write(message, '(2a)') ch10,&
&       "   - ( in the following only the value for the first k-points are printed)"
       call wrtout(std_out,message,'COLL')
       do ikpt=1,nkptr
         if(option<5) then
           write(message, '(2a,i4,2x,f14.5,a)') ch10,&
&           "   -k-pt--",ikpt,oper%wtk(ikpt),"(<-weight(k-pt))"
           call wrtout(std_out,message,'COLL')
         else if(abs(prtopt)>=4.or.option>8) then
           write(message, '(2a,i5,a,i5,a,i5)') ch10,"  Writes occupations for k-pt",&
&           ikpt, "and between bands",&
&           iband1," and",iband2
           call wrtout(std_out,message,'COLL')
         endif
         do ib=1,oper%mbandc
           if(option<5) then
             if(abs(aimag(oper%ks(isppol,ikpt,ib,ib))).ge.tol10) then
               write(message, '(a,i4,e14.5,3x,e14.5,3x,e21.14)') "   -iband--",ib,&
&               paw_dmft%eigen_dft(isppol,ikpt,ib),oper%ks(isppol,ikpt,ib,ib)
               call wrtout(std_out,message,'COLL')
             else
               write(message, '(a,i4,e14.5,3x,e14.5)') "   -iband--",ib,&
&               paw_dmft%eigen_dft(isppol,ikpt,ib),real(oper%ks(isppol,ikpt,ib,ib))
               call wrtout(std_out,message,'COLL')
             endif
           endif
           if(abs(prtopt)>=4.or.option>8.and.ib>=iband1.and.ib<=iband2) then
             write(message, '(2000(f8.3))') &
&               (real(oper%ks(isppol,ikpt,ib,ib1)),ib1=iband1,iband2)
             call wrtout(std_out,message,'COLL')
             write(message, '(2000(f8.3))') &
&               (aimag(oper%ks(isppol,ikpt,ib,ib1)),ib1=iband1,iband2)
             call wrtout(std_out,message,'COLL')
             write(message, '(2000(f8.3))') &
&               (sqrt(real(oper%ks(isppol,ikpt,ib,ib1))**2+aimag(oper%ks(isppol,ikpt,ib,ib1))**2),ib1=iband1,iband2)
             call wrtout(std_out,message,'COLL')
!   to write imaginary part
!             write(message, '(1000(2f9.3,2x))') &
!&               (real(oper%ks(isppol,ikpt,ib,ib1)),imag(oper%ks(isppol,ikpt,ib,ib1)),ib1=iband1,iband2)
!             call wrtout(std_out,message,'COLL')
           endif ! prtopt>=20
           do ib1=1,oper%mbandc
             if(abs(aimag(oper%ks(isppol,ikpt,ib,ib1)))>max(tol10,maximag(1))) then
               ximag=.true.
               maximag(1)=aimag(oper%ks(isppol,ikpt,ib,ib1))
             endif
           enddo
         enddo ! ib
       enddo ! ikpt
     enddo ! isppol
   else
    write(message,'(5x,a,i10,a)') '(not written)'
    call wrtout(std_out,message,'COLL')
   endif
   if(ximag) then
     write(message, '(3a,e12.4,a)')"Occupations are imaginary !",ch10, &
&    "  Maximal value is ", maximag(1), ch10
     MSG_WARNING(message)
   endif
 else if(abs(prtopt)>=3.and.((option<5).or.(option>8))) then
   write(message, '(2a)') ch10," Prb with options and has_operks in print_oper"
   call wrtout(std_out,message,'COLL')
 endif ! if oper%has_operks
! write(message, '(2a)') ch10," end print_oper"
!     call wrtout(std_out,message,'COLL')


 DBG_EXIT("COLL")
end subroutine print_oper
!!***

!!****f* m_oper/inverse_oper
!! NAME
!! inverse_oper
!!
!! FUNCTION
!!  Compute the inverse of the operator either in the KS space or in the
!!  correlated subspace.
!!
!! INPUTS
!!  oper <type(oper_type)>= operator
!!  paw_dmft <type(paw_dmft_type)>
!!  option= integer
!!
!! OUTPUT
!!  oper <type(oper_type)>= operator inverted
!!
!! PARENTS
!!      m_dmft,m_forctqmc,m_green
!!
!! CHILDREN
!!      prod_matlu
!!
!! SOURCE

subroutine inverse_oper(oper,option,prtopt,procb,iproc)

 use defs_basis
 use m_crystal, only : crystal_t
 use m_paw_dmft, only : paw_dmft_type
 use m_matlu, only : inverse_matlu
 use m_errors

!Arguments ------------------------------------
!type
 integer, intent(in):: option
 integer, intent(in) :: prtopt
 type(oper_type),intent(inout) :: oper
 integer, optional, intent(in) ::  iproc
 integer, optional, intent(in) :: procb(oper%nkpt)
!oper variables-------------------------------
 integer :: ikpt,isppol
 complex(dpc), allocatable :: matrix(:,:)
 character(len=500) :: message
 integer :: paral
 integer, allocatable :: procb2(:)
!todo_ba: prb with gwpc here: necessary for matcginv but should be dpc
! *********************************************************************
 DBG_ENTER("COLL")
 if(option==2.or.option==3) then
   ABI_ALLOCATE(procb2,(oper%nkpt))
 endif
!  if option=1 do inversion in local space
!  if option=2 do inversion in KS band space
!  if option=3 do both
 if(present(procb).and.present(iproc)) then
   paral=1
   procb2=procb
 else
   paral=0
 endif

 if(((option==1.or.option==3).and.(oper%has_opermatlu==0)).or.&
&   ((option==2.or.option==3).and.(oper%has_operks==0))) then
   message = " Options are not coherent with definitions of this operator"
   MSG_ERROR(message)
 endif

 if(option==1.or.option==3) then
   call inverse_matlu(oper%matlu,oper%natom,prtopt)
 else if(option==2.or.option==3) then
   ABI_ALLOCATE(matrix,(oper%mbandc,oper%mbandc))
     do isppol=1,oper%nsppol
       do ikpt=1,oper%nkpt
        if ((paral==1.and.(procb2(ikpt)==iproc)).or.(paral==0)) then
          matrix(:,:)=oper%ks(isppol,ikpt,:,:)
!          write(std_out,*) "isppol,ikpt",isppol,ikpt,m
!          write(std_out,*) "isppol,ikpt",matrix
         !call matcginv_dpc(matrix,oper%mbandc,oper%mbandc)
         call xginv(matrix,oper%mbandc)
         oper%ks(isppol,ikpt,:,:)=matrix(:,:)
        endif
       enddo ! ikpt
     enddo ! isppol
   ABI_DEALLOCATE(matrix)
 endif

 if(option==2.or.option==3) then
   ABI_DEALLOCATE(procb2)
 endif
 DBG_EXIT("COLL")
end subroutine inverse_oper
!!***

!!****f* m_oper/loc_oper
!! NAME
!! loc_oper
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_datafordmft,m_green
!!
!! CHILDREN
!!      prod_matlu
!!
!! SOURCE

subroutine loc_oper(oper,paw_dmft,option,jkpt,procb,iproc)

 use defs_basis
 use m_paw_dmft, only : paw_dmft_type
 use m_errors

!Arguments ------------------------------------
!type
 integer, intent(in):: option
 integer, optional, intent(in):: jkpt
 type(oper_type),intent(inout) :: oper
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, optional, intent(in) ::  iproc
 integer, optional, intent(in) :: procb(oper%nkpt)
!oper variables-------------------------------
 integer :: iatom,ib,ib1,ikpt,ikpt1,ispinor,ispinor1,isppol,im1,im
 integer :: natom,mbandc,ndim,nkpt,nspinor,nsppol,paral
 character(len=500) :: message
 integer, allocatable :: procb2(:)
 logical lvz  !vz_d
! *********************************************************************
 DBG_ENTER("COLL")
 ABI_ALLOCATE(procb2,(oper%nkpt))
 if((oper%has_opermatlu==0).or.(oper%has_operks==0)) then
   message = " Operator is not defined to be used in loc_oper"
   MSG_ERROR(message)
 endif
 if(present(procb).and.present(iproc)) then
   paral=1
   procb2=procb
 else
   paral=0
 endif

 if(option<0) then
 endif
 nkpt=oper%nkpt
 natom=oper%natom
 nsppol=oper%nsppol
 mbandc=oper%mbandc
 nspinor=oper%nspinor


 do iatom=1,natom
   oper%matlu(iatom)%mat=czero
 enddo
 do isppol=1,nsppol
  do ikpt=1,nkpt
   ikpt1=ikpt
   if(present(jkpt)) ikpt1=jkpt
   lvz=paral==0  !vz_d
   if(present(iproc)) lvz=lvz.or.(paral==1.and.(procb2(ikpt1)==iproc))  !vz_d
!  if ((paral==1.and.(procb2(ikpt1)==iproc)).or.(paral==0)) then    !vz_d
   if(lvz) then !vz_d
   do ib=1,mbandc
    do ib1=1,mbandc
     do iatom=1,natom
      if(oper%matlu(iatom)%lpawu.ne.-1) then
      ndim=2*oper%matlu(iatom)%lpawu+1
       do im=1,ndim
        do im1=1,ndim
         do ispinor=1,nspinor
          do ispinor1=1,nspinor
            if (im1 == im .and. im1 == 1 .and. ib1 == ib .and. iatom == 1) then

            end if
            oper%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=     &
&            oper%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)+    &
&            paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)*        &
&            conjg(paw_dmft%psichi(isppol,ikpt1,ib1,ispinor1,iatom,im1))* &
!false&            paw_dmft%psichi(isppol,ikpt1,ib1,ispinor1,iatom,im1)*        &
!false&            conjg(paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im))* &
&            oper%ks(isppol,ikpt,ib,ib1)*oper%wtk(ikpt)
! one  could suppress wtk here if present(jkpt)
! ks(ib,ib1)=ks(ib1,ib) -> ib and ib1 can be underchanged !
          enddo ! ispinor1
         enddo ! ispinor
        enddo ! im1
       enddo ! im
      endif
     enddo ! iatom
    enddo ! ib
   enddo ! ib
   endif
  enddo ! ikpt
 enddo ! isppol
 ABI_DEALLOCATE(procb2)




 DBG_EXIT("COLL")
end subroutine loc_oper
!!***

!!****f* m_oper/upfold_oper
!! NAME
!! upfold_oper
!!
!! FUNCTION
!!
!! INPUTS
!!  psichi=<chi|psi> !
!!
!! OUTPUT
!!
!! PARENTS
!!      m_green
!!
!! CHILDREN
!!      prod_matlu
!!
!! SOURCE

subroutine upfold_oper(oper,paw_dmft,option,procb,iproc,prt)

 use defs_basis
 use m_paw_dmft, only : paw_dmft_type
 use m_errors

!Arguments ------------------------------------
!type
 integer, intent(in):: option
 type(oper_type),intent(inout) :: oper
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, optional, intent(in) :: iproc
 integer, optional, intent(in) :: procb(oper%nkpt)
 integer, optional, intent(in) :: prt
!oper variables-------------------------------
 integer :: iatom,ib,ib1,ikpt,ispinor,ispinor1,isppol,im1,im
 integer :: natom,mbandc,ndim,nkpt,nspinor,nsppol,paral
 integer, allocatable :: procb2(:)
 character(len=500) :: message
! *********************************************************************

 ABI_UNUSED(prt)
 ABI_ALLOCATE(procb2,(oper%nkpt))
 if(present(procb).and.present(iproc)) then
   paral=1
   procb2=procb
!   write(6,*) "upfold_oper procb",procb
!   write(6,*) "iproc",iproc
!   write(6,*) size(procb)
!   write(6,*) size(procb2)
!   write(6,*) procb2(1),procb2(16)
 else
   paral=0
 endif



 DBG_ENTER("COLL")
 if((oper%has_opermatlu==0).or.(oper%has_operks==0)) then
   message = " Operator is not defined to be used in upfold_oper"
   MSG_ERROR(message)
 endif
 if(option<0) then
 endif
 nkpt=oper%nkpt
 natom=paw_dmft%natom
 nsppol=paw_dmft%nsppol
 mbandc=paw_dmft%mbandc
 nspinor=paw_dmft%nspinor

 do isppol=1,nsppol
   do ikpt=1,nkpt
    if ((paral==1.and.(procb2(ikpt)==iproc)).or.(paral==0)) then
     do ib=1,mbandc
       do ib1=1,mbandc
!               if(ib==1.and.ib1==3) write(std_out,*) "IKPT=",ikpt
        oper%ks(isppol,ikpt,ib,ib1)=czero

         do iatom=1,natom
           if(oper%matlu(iatom)%lpawu.ne.-1) then
             ndim=2*oper%matlu(iatom)%lpawu+1
             do im=1,ndim
               do im1=1,ndim
                 do ispinor=1,nspinor
                   do ispinor1=1,nspinor

! psichi(isppol,ikpt,ib,ispinor,iatom,im)=<\chi_{m,R,ispinor)|\Psi(s,k,nu)>
                     oper%ks(isppol,ikpt,ib,ib1)= oper%ks(isppol,ikpt,ib,ib1) &
&                     + ( paw_dmft%psichi(isppol,ikpt,ib1,ispinor1,iatom,im1)        &
&                     * oper%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)    &
&                     * conjg(paw_dmft%psichi(isppol,ikpt,ib,ispinor,iatom,im)))
              ! if(present(prt).and.(ib==1.and.ib1==1)) then
              !   write(6,*) "im,im1",im,im1
              !   write(6,*) "ispinor,ispinor1",ispinor,ispinor1
              !   write(6,*) "psichi",paw_dmft%psichi(isppol,ikpt,ib1,ispinor1,iatom,im1)
              !   write(6,*) "psichi 2",paw_dmft%psichi(isppol,ikpt,ib,ispinor,iatom,im1)
              !   write(6,*) "oper%matlu", oper%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)
              ! endif

                   enddo ! ispinor1
                 enddo ! ispinor
               enddo ! im1
             enddo ! im
           endif
         enddo ! iatom

       enddo ! ib
     enddo ! ib
    endif
   enddo ! ikpt
 enddo ! isppol
 ABI_DEALLOCATE(procb2)

 DBG_EXIT("COLL")
end subroutine upfold_oper
!!***

!!****f* m_oper/identity_oper
!! NAME
!! identity_oper
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_datafordmft
!!
!! CHILDREN
!!      prod_matlu
!!
!! SOURCE

subroutine identity_oper(oper,option)

 use defs_basis
 use m_crystal, only : crystal_t
 use m_paw_dmft, only : paw_dmft_type
 use m_errors

!Arguments ------------------------------------
!type
 integer, intent(in):: option
 type(oper_type),intent(inout) :: oper
!oper variables-------------------------------
 integer :: iatom,ib,ikpt,ispinor,isppol,im
 integer :: natom,mbandc,ndim,nkpt,nspinor,nsppol
 character(len=500) :: message
! *********************************************************************

 DBG_ENTER("COLL")

 if(((option==1.or.option==3).and.(oper%has_opermatlu==0)).or.&
&   ((option==2.or.option==3).and.(oper%has_operks==0))) then
   message = " Options in identity_oper are not coherent with definitions of this operator"
   MSG_ERROR(message)
 endif
 nkpt=oper%nkpt
 nsppol=oper%nsppol
 mbandc=oper%mbandc
 natom=oper%natom
 nspinor=oper%nspinor
 oper%ks=czero
 do iatom=1,natom
  oper%matlu(iatom)%mat= czero
 enddo

 if(option==1.or.option==3) then
   do isppol=1,nsppol
     do ikpt=1,nkpt
       do ib=1,mbandc
         oper%ks(isppol,ikpt,ib,ib)=cone
       enddo ! ib
     enddo ! ikpt
   enddo ! isppol

 else if (option==2.or.option==3) then
   do iatom=1,natom
     if(oper%matlu(iatom)%lpawu.ne.-1) then
       ndim=2*oper%matlu(iatom)%lpawu+1
       do isppol=1,nsppol
         do im=1,ndim
           do ispinor=1,nspinor
             oper%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)= cone
           enddo ! ispinor
         enddo ! im
       enddo
     endif ! lpawu
   enddo ! iatom
 endif

 DBG_EXIT("COLL")
end subroutine identity_oper
!!***

!!****f* m_oper/diff_oper
!! NAME
!! diff_oper
!!
!! FUNCTION
!! Compute a norm of the differences between two occupations matrices.
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  occup1 <type(oper_type)>= occupations
!!  occup2 <type(oper_type)>= occupations
!!  option : option for printing (if 1 assume data are related to lda only)
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!
!! PARENTS
!!      m_datafordmft,m_dmft
!!
!! CHILDREN
!!      prod_matlu
!!
!! SOURCE

subroutine diff_oper(char1,char2,occup1,occup2,option,toldiff)

 use defs_basis
 use m_paw_dmft, only : paw_dmft_type
 use m_crystal, only : crystal_t
 use m_matlu, only : diff_matlu
 use m_errors

!Arguments ------------------------------------
!type
 type(oper_type), intent(in) :: occup1,occup2
 integer, intent(in) :: option
 real(dp), intent(in) :: toldiff
 character(len=*), intent(in) :: char1,char2
!local variables-------------------------------
 integer :: mbandc,nkpt
 character(len=500) :: message
! *********************************************************************

 DBG_ENTER("COLL")
 if(occup1%has_opermatlu==0.or.occup2%has_opermatlu==0) then
   message = " Operators are not defined to be used in diff_oper"
   MSG_ERROR(message)
 endif
 mbandc   = occup1%mbandc
 nkpt    = occup1%nkpt
 if(occup1%nkpt/=occup2%nkpt) then
  write(message,'(a,2x,2i9)')' Operators are not equals',occup1%nkpt,occup2%nkpt
  MSG_ERROR(message)
 endif

 call diff_matlu(char1,char2,occup1%matlu,&
& occup2%matlu,occup1%natom,option,toldiff)
! if(option==1) then
!  toldiff=tol4
!  if( matludiff < toldiff ) then
!   write(message,'(6a,e12.4,a,e12.4)') ch10,&
!&   '   Differences between ',trim(char1),' and ',trim(char2),' is small enough:',&
!&   matludiff,'is lower than',toldiff
!   call wrtout(std_out,message,'COLL')
!  else
!   write(message,'(6a,e12.4,a,e12.4)') ch10,&
!&   '   Error: Differences between ',trim(char1),' and ',trim(char2),' is too large:',&
!&   matludiff,'is large than',toldiff
!   call wrtout(std_out,message,'COLL')
!   call abi_abort('COLL')
!  endif
! endif
! call abi_abort('COLL')

 DBG_EXIT("COLL")
end subroutine diff_oper
!!***

!!****f* m_oper/trace_oper
!! NAME
!! trace_oper
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dmft,m_green
!!
!! CHILDREN
!!      prod_matlu
!!
!! SOURCE

subroutine trace_oper(oper,trace_ks,trace_loc,opt_ksloc)

 use defs_basis
 use m_matlu, only : trace_matlu
 use m_errors

!Arguments ------------------------------------
!type
 type(oper_type),intent(in) :: oper
 real(dp), intent(out) :: trace_ks  !vz_i
 real(dp), intent(inout) :: trace_loc(oper%natom,oper%nsppol+1) !vz_i
 integer, intent(in) :: opt_ksloc

!oper variables-------------------------------
 integer ::  ib, ikpt, isppol
 character(len=500) :: message
 real(dp) :: temp1
! *********************************************************************

 DBG_ENTER("COLL")
 if(((opt_ksloc==1.or.opt_ksloc==3).and.(oper%has_opermatlu==0)).or.&
&   ((opt_ksloc==2.or.opt_ksloc==3).and.(oper%has_operks==0))) then
   message = " Options in trace_oper are not coherent with definitions of this operator"
   MSG_ERROR(message)
 endif

 if(opt_ksloc==1.or.opt_ksloc==3) then
   trace_ks=zero
   temp1=zero
   do isppol=1,oper%nsppol
     do ikpt=1,oper%nkpt
       do ib=1,oper%mbandc
         trace_ks=trace_ks+real(oper%ks(isppol,ikpt,ib,ib)*oper%wtk(ikpt))
         temp1=temp1+oper%wtk(ikpt)
       enddo
     enddo
   enddo
   if(oper%nsppol==1.and.oper%nspinor==1) trace_ks=two*trace_ks
!   write(std_out,*) "temp1",temp1
 endif
 if(opt_ksloc==2.or.opt_ksloc==3) then
   trace_loc=zero
   call trace_matlu(oper%matlu,oper%natom,trace_loc)
 endif

 DBG_EXIT("COLL")
end subroutine trace_oper
!!***

!!****f* m_oper/prod_oper
!! NAME
!! prod_oper
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_datafordmft
!!
!! CHILDREN
!!      prod_matlu
!!
!! SOURCE

subroutine prod_oper(oper1,oper2,oper3,opt_ksloc)

 use defs_basis
 use m_errors
 use m_matlu, only : prod_matlu

!Arguments ------------------------------------
!type
 type(oper_type),intent(in) :: oper1
 type(oper_type),intent(in) :: oper2
 type(oper_type),intent(inout) :: oper3
 integer :: opt_ksloc

!oper variables-------------------------------
 integer ::  ib, ib1, ib2, ikpt, isppol
! *********************************************************************
 DBG_ENTER("COLL")
 if(opt_ksloc==2) then
   if(oper1%has_opermatlu==1.and.oper2%has_opermatlu==1.and.oper3%has_opermatlu==1)  then
     call prod_matlu(oper1%matlu,oper2%matlu,oper3%matlu,oper1%natom) !
   endif
 endif

 if(opt_ksloc==1) then
   if(oper1%has_operks==1.and.oper2%has_operks==1.and.oper3%has_operks==1) then
     do isppol=1,oper1%nsppol
       do ikpt=1,oper1%nkpt
         do ib=1,oper1%mbandc
           do ib1=1,oper1%mbandc
             oper3%ks(isppol,ikpt,ib,ib1)=czero
             do ib2=1,oper1%mbandc
               oper3%ks(isppol,ikpt,ib,ib1)=oper3%ks(isppol,ikpt,ib,ib1)+ oper1%ks(isppol,ikpt,ib,ib2)*oper2%ks(isppol,ikpt,ib2,ib1)
             enddo
           enddo
         enddo
       enddo
     enddo
   endif
 endif

 DBG_EXIT("COLL")
end subroutine prod_oper
!!***

END MODULE m_oper
!!***
