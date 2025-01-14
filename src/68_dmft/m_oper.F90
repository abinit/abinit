!!****m* ABINIT/m_oper
!! NAME
!!  m_oper
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2024 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

 use m_matlu, only : matlu_type

 implicit none

 private

 public :: init_oper
 public :: diff_oper
 public :: destroy_oper
 public :: print_oper
 public :: inverse_oper
 public :: downfold_oper
 public :: identity_oper
 public :: copy_oper
 public :: trace_oper
 public :: upfold_oper
 public :: prod_oper
 public :: trace_prod_oper
 public :: gather_oper
 public :: gather_oper_ks
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

  integer :: has_operks
  ! Is the operator allocated in the KS basis ?

  integer :: has_opermatlu
  ! Is the operator allocated in the local basis ?
!
  integer :: mbandc
  ! Total number of correlated bands

  integer :: natom
  ! Number of atoms

  integer :: nkpt
  ! Number of k-point in the IBZ.
!
  integer :: nspinor
  ! Number of spinors
!
  integer :: nsppol
  ! Number of spin polarizations

  integer :: paral
  ! =1 if the operator has been memory-parallelized over kpt, 0 otherwise

  integer :: shiftk
  ! Shift to get the physical kpt index (when the operator is memory-parallelized over kpt)

  !character(len=12) :: whichoper
  ! describe the type of operator computed (DFT, DMFT, KS..)

!  ! Polarisation
  type(matlu_type), allocatable :: matlu(:)
  ! Local projection on correlated orbitals

  complex(dpc), allocatable :: ks(:,:,:,:)
  ! In the KS basis  (mbandc,mbandc,nkpt,nsppol)

  real(dp), ABI_CONTIGUOUS pointer :: wtk(:) => null()
  ! Weights for each kpt

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
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  nkpt = number of k-pts
!!  wtk = weights for each k-pt
!!  shiftk = shift for the kpt index
!!  opt_ksloc = 1: initialize in KS space only
!!            = 2: initialize in local space only
!!            = 3: initialize in both KS and local space
!!
!! OUTPUTS
!!  oper <type(oper_type)>= operator
!!
!! SOURCE

subroutine init_oper(paw_dmft,oper,nkpt,wtk,shiftk,opt_ksloc)

 use m_matlu, only : init_matlu
 use m_paw_dmft, only : paw_dmft_type

!Arguments ------------------------------------
 integer, optional, intent(in) :: nkpt,opt_ksloc,shiftk
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(oper_type), intent(inout) :: oper
 real(dp), target, optional :: wtk(paw_dmft%nkpt)
!Local variables ------------------------------------
 integer :: optksloc
!************************************************************************

 DBG_ENTER("COLL")

 optksloc = 3
 if (present(opt_ksloc)) optksloc = opt_ksloc

 !if(optksloc/=3) then
    ! FIXME: empty line!
 !endif

 oper%has_operks    = 0
 oper%has_opermatlu = 0

! ===================
!  Integers
! ===================
 oper%mbandc  = paw_dmft%mbandc
 oper%natom   = paw_dmft%natom
 oper%nspinor = paw_dmft%nspinor
 oper%nsppol  = paw_dmft%nsppol
 oper%paral   = 0
 oper%shiftk  = 0

 if (present(shiftk)) oper%shiftk = shiftk

 oper%nkpt = paw_dmft%nkpt
 if (present(nkpt)) oper%nkpt = nkpt

 if (present(shiftk) .or. oper%nkpt /= paw_dmft%nkpt) oper%paral = 1

! allocate(oper%wtk(oper%nkpt))
 if (present(wtk)) then
   oper%wtk => wtk(:)
 else
   oper%wtk => paw_dmft%wtk(:)
 end if ! present(wtk)

! ===================
!  KS variables
! ===================
 if (optksloc == 1 .or. optksloc == 3) then

   ABI_MALLOC(oper%ks,(oper%mbandc,oper%mbandc,oper%nkpt,oper%nsppol))
   oper%has_operks  = 1
   oper%ks(:,:,:,:) = czero

 end if ! optksloc=1 or optksloc=3

! ===================
!  matlu variables
! ===================
 if (optksloc == 2 .or. optksloc == 3) then
   ABI_MALLOC(oper%matlu,(oper%natom))
   oper%has_opermatlu = 1
   call init_matlu(oper%natom,oper%nspinor,oper%nsppol,paw_dmft%lpawu(:),oper%matlu(:))
 end if ! optksloc=2 or optksloc=3

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
!!  oper <type(oper_type)>= operator
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_oper(oper)

 use m_matlu, only : destroy_matlu

!Arguments ------------------------------------
 type(oper_type), intent(inout) :: oper
!Local variables-------------------------------
!! *********************************************************************

 DBG_ENTER("COLL")

 if (oper%has_opermatlu == 1) then
   call destroy_matlu(oper%matlu(:),oper%natom)
 !else
 !  message = " Operator is not defined to be used in destroy_oper"
 !  ABI_ERROR(message)
 end if ! has_opermatlu=1

 if (allocated(oper%matlu)) then
   ABI_FREE(oper%matlu)
   oper%has_opermatlu = 0
 end if

 if (allocated(oper%ks)) then
   ABI_FREE(oper%ks)
   oper%has_operks = 0
 end if

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
!!  Copy oper1 into oper2
!!
!! INPUTS
!!  oper1 <type(oper_type)>= operator
!!
!! OUTPUT
!!  oper2 <type(oper_type)>= operator
!!
!! SOURCE

subroutine copy_oper(oper1,oper2)

 use m_matlu, only : copy_matlu

!Arguments ------------------------------------
 type(oper_type), intent(in) :: oper1
 type(oper_type), intent(inout) :: oper2 !vz_i
!Local variables-------------------------------
! *********************************************************************

 DBG_ENTER("COLL")

 if (oper1%has_opermatlu == 1 .and. oper2%has_opermatlu == 1) then
   call copy_matlu(oper1%matlu(:),oper2%matlu(:),oper1%natom)
 end if

 if (oper1%has_operks == 1 .and. oper2%has_operks == 1) &
    & oper2%ks(:,:,:,:) = oper1%ks(:,:,:,:)

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
!! oper <type(oper_type)>= operator
!! option= < 5: write diagonal part of KS occupation matrix
!!         > 8: write all elements of KS occup. matrix.
!! paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!! prtopt= in local space: print option for print_matlu
!!         in KS space: only prints if abs(prtopt)>=3
!!                      print off-diagonal elements if abs(prtopt)>=4
!!
!! OUTPUT
!!
!! SOURCE

subroutine print_oper(oper,option,paw_dmft,prtopt)

 use m_matlu, only : print_matlu
 use m_paw_dmft, only : paw_dmft_type

!Arguments ------------------------------------
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(oper_type), intent(in) :: oper
 integer, intent(in) :: option,prtopt
!Local variables-------------------------------
 integer :: ib,ib1,iband1,iband2,ikpt,isppol,mbandc,nkpt,nkptr
 character(len=100000) :: message
 logical  :: ximag
 real(dp) :: maximag
! *********************************************************************

 DBG_ENTER("COLL")

 if (oper%has_opermatlu == 1) then
   write(message,'(2a)') ch10,'   = In the atomic basis'
   call wrtout(std_out,message,'COLL')
   call print_matlu(oper%matlu(:),oper%natom,prtopt)
 end if ! has_opermatlu=1

 if (oper%has_operks == 1) then
   write(message,'(2a)') ch10,'   = In the Kohn-Sham basis'
   call wrtout(std_out,message,'COLL')

!todo_ba complete print_out
   mbandc  = oper%mbandc
   iband1  = 1
   iband2  = mbandc
   maximag = zero
   nkpt    = oper%nkpt
   ximag   = .false.
!   do ib=1,oper%mbandc
!     if(-(paw_dmft%eigen_dft(1,1,ib)+paw_dmft%fermie).ge.0.3) iband1=ib
!     if( (paw_dmft%eigen_dft(1,1,ib)-paw_dmft%fermie).le.0.3) iband2=ib
!   enddo

   if (abs(prtopt) >= 3 .and. ((option < 5) .or. (option > 8))) then
!     write(message,'(x,a,a,i4,2x,a)') ch10,'  -KS states'
!     call wrtout(std_out,message,'COLL')
     nkptr = min(nkpt,4)
     do isppol=1,paw_dmft%nsppol
       write(message,'(a,3x,a,1x,i1)') ch10,"--isppol--",isppol
       call wrtout(std_out,message,'COLL')
       write(message,'(2a)') ch10,&
         & "   - (in the following only the values for the correlated bands and the first k-points are printed)"
       call wrtout(std_out,message,'COLL')
       do ikpt=1,nkptr
         write(message,'(2a,i4,2x,f14.5,a)') ch10,&
             & "   -k-pt--",ikpt,oper%wtk(ikpt),"(<-weight(k-pt))"
         call wrtout(std_out,message,'COLL')
         if (option < 5) then
           write(message,'(19x,a,6x,a)') "Eigenvalues","Occupations"
           call wrtout(std_out,message,'COLL')
         else if (abs(prtopt) >= 4 .or. option > 8) then
           write(message,'(a,10x,2000(i5,12x))') ch10,(paw_dmft%include_bands(ib),ib=iband1,iband2)
           call wrtout(std_out,message,'COLL')
         end if ! option
         do ib=1,mbandc
           if (option < 5) then
             if (abs(aimag(oper%ks(ib,ib,ikpt,isppol))) >= tol10) then
               write(message,'(a,i5,e14.5,3x,e14.5,3x,e21.14)') "   -iband--",paw_dmft%include_bands(ib),&
                 & paw_dmft%eigen_dft(ib,ikpt,isppol),oper%ks(ib,ib,ikpt,isppol)
             else
               write(message,'(a,i5,e14.5,3x,e14.5)') "   -iband--",paw_dmft%include_bands(ib),&
                 & paw_dmft%eigen_dft(ib,ikpt,isppol),dble(oper%ks(ib,ib,ikpt,isppol))
             end if ! imaginary part
             call wrtout(std_out,message,'COLL')
           end if ! option<5
           if (abs(prtopt) >= 4 .or. option > 8 .and. ib >= iband1 .and. ib <= iband2) then

             write(message,'(i5,1x,2000(2f7.3,3x))') paw_dmft%include_bands(ib),(dble(oper%ks(ib,ib1,ikpt,isppol)), &
                 & aimag(oper%ks(ib,ib1,ikpt,isppol)),ib1=iband1,iband2)
             call wrtout(std_out,message,'COLL')

!   to write imaginary part
!             write(message, '(1000(2f9.3,2x))') &
!&               (real(oper%ks(isppol,ikpt,ib,ib1)),imag(oper%ks(isppol,ikpt,ib,ib1)),ib1=iband1,iband2)
!             call wrtout(std_out,message,'COLL')
           end if ! prtopt>=20
           if (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7) then ! only need to check diagonal elements
             if (abs(aimag(oper%ks(ib,ib,ikpt,isppol))) > max(tol10,maximag)) then
               ximag   = .true.
               maximag = aimag(oper%ks(ib,ib,ikpt,isppol))
             end if
           else
             do ib1=1,mbandc
               if (abs(aimag(oper%ks(ib1,ib,ikpt,isppol))) > max(tol10,maximag)) then
                 ximag   = .true.
                 maximag = aimag(oper%ks(ib1,ib,ikpt,isppol))
               end if
             end do ! ib1
           end if
         end do ! ib
       end do ! ikpt
     end do ! isppol
   else
    write(message,'(5x,a,i10,a)') '(not written)'
    call wrtout(std_out,message,'COLL')
   end if ! abs(prtopt)>=3 and (option<5 or option>8)
   if (ximag) then
     write(message,'(3a,e12.4,a)') "Occupations are imaginary !",ch10,&
        & "  Maximal value is ",maximag,ch10
     ABI_WARNING(message)
   end if ! ximag
 else if (abs(prtopt) >= 3 .and. ((option < 5) .or. (option > 8))) then
   write(message, '(2a)') ch10," Prb with options and has_operks in print_oper"
   call wrtout(std_out,message,'COLL')
 end if ! if oper%has_operks
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
!!  option=1 do inversion in KS band space
!!        =2 do inversion in local space
!!        =3 do both
!!  procb(ikpt)=for kpt parallelization; gives the rank (in the kpt communicator) of the CPU handling each ikpt
!!  iproc=rank of the current process in the kpt communicator
!!
!! OUTPUT
!!  oper <type(oper_type)>= operator inverted
!!
!! SOURCE

subroutine inverse_oper(oper,option,procb,iproc)

 use m_matlu, only : inverse_matlu
 use m_hide_lapack, only : xginv

!Arguments ------------------------------------
 integer, intent(in) :: option
 type(oper_type), intent(inout) :: oper
 integer, optional, intent(in) :: iproc
 integer, optional, intent(in) :: procb(oper%nkpt)
!Local variables-------------------------------
 integer :: ikpt,isppol,paral
!todo_ba: prb with gwpc here: necessary for matcginv but should be dpc
! *********************************************************************

 DBG_ENTER("COLL")

 paral = 0
 if (present(procb) .and. present(iproc) .and. oper%paral == 0) paral = 1

 !if (((option == 1 .or. option == 3) .and. (oper%has_operks == 0)) .or. &
 !  & ((option == 2 .or. option == 3) .and. (oper%has_opermatlu == 0))) then
 !  message = " Options are not coherent with definitions of this operator"
 !  ABI_ERROR(message)
 !end if

 if (option == 2 .or. option == 3) then
   call inverse_matlu(oper%matlu(:),oper%natom)
 end if

 if (option == 1 .or. option == 3) then
   do isppol=1,oper%nsppol
     do ikpt=1,oper%nkpt
       if (paral == 1) then
         if (procb(ikpt) /= iproc) cycle
       end if
!          write(std_out,*) "isppol,ikpt",isppol,ikpt,m
!          write(std_out,*) "isppol,ikpt",matrix
         !call matcginv_dpc(matrix,oper%mbandc,oper%mbandc)
       call xginv(oper%ks(:,:,ikpt,isppol),oper%mbandc)
     end do ! ikpt
   end do ! isppol
 end if ! option

 DBG_EXIT("COLL")

end subroutine inverse_oper
!!***

!!****f* m_oper/downfold_oper
!! NAME
!! downfold_oper
!!
!! FUNCTION
!!  Downfold an operator from KS space to local space.
!!
!! INPUTS
!!  oper <type(oper_type)>= operator
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  procb(ikpt)=for kpt parallelization; gives the rank (in the kpt communicator) of the CPU handling each ikpt
!!  iproc=rank of the current process in the kpt communicator
!!  option = 1 (default) : downfold an operator represented by a matrix in KS space
!!         = 2 : downfold the identity
!!         = 3 : downfold a diagonal KS operator
!!         = 4 : computes downfold(upfold)
!!  op_ks_diag = when option=3, you can provide the diagonal KS operator in this variable
!!               with the format mbandc*nkpt*nsppol, instead of storing it in the mband*mband matrix
!!               of oper%ks
!!
!! OUTPUT
!!
!! SOURCE

subroutine downfold_oper(oper,paw_dmft,procb,iproc,option,op_ks_diag)

 use m_paw_dmft, only : paw_dmft_type
 use m_abi_linalg, only : abi_xgemm

!Arguments ------------------------------------
 type(oper_type), intent(inout) :: oper
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, optional, intent(in) :: iproc,option
 integer, optional, intent(in) :: procb(oper%nkpt)
 real(dp), optional, intent(in) :: op_ks_diag(oper%mbandc,oper%nkpt,oper%nsppol)
!oper variables-------------------------------
 integer :: iatom,ib,ik,ikpt,isppol,lpawu,mbandc,ndim
 integer :: ndim_max,nspinor,opt,paral,shift
 character(len=500) :: message
 complex(dpc), allocatable :: mat_temp(:,:),mat_temp2(:,:),mat_temp3(:,:)
! *********************************************************************

 DBG_ENTER("COLL")

 if (oper%has_opermatlu == 0) then
   message = " Operator is not defined to be used in downfold_oper"
   ABI_ERROR(message)
 end if

 mbandc   = oper%mbandc
 nspinor  = oper%nspinor
 ndim_max = nspinor * (2*paw_dmft%maxlpawu+1)
 paral    = 0
 shift    = oper%shiftk

 if (present(procb) .and. present(iproc) .and. oper%paral == 0) paral = 1

 opt = 1
 if (present(option)) opt = option

 do iatom=1,oper%natom
   lpawu = oper%matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   oper%matlu(iatom)%mat(:,:,:) = czero
   ndim = nspinor * (2*lpawu+1)
   ABI_MALLOC(mat_temp,(ndim,mbandc))
   ABI_MALLOC(mat_temp2,(ndim,ndim))
   ABI_MALLOC(mat_temp3,(ndim,ndim))
   do isppol=1,oper%nsppol
     do ikpt=1,oper%nkpt ! index of kpt on the current CPU

       if (paral == 1) then
         if (procb(ikpt) /= iproc) cycle
       end if

       ik = ikpt + shift ! true kpt index (needed for chipsi)

       if (opt == 1 .or. opt == 3) then

         if (opt == 1) then

           call abi_xgemm("n","n",ndim,mbandc,mbandc,cone,paw_dmft%chipsi(:,:,ik,isppol,iatom),&
                        & ndim_max,oper%ks(:,:,ikpt,isppol),mbandc,czero,mat_temp(:,:),ndim)

         else if (opt == 3) then

           do ib=1,mbandc
             if (present(op_ks_diag)) then
               mat_temp(:,ib) = paw_dmft%chipsi(1:ndim,ib,ik,isppol,iatom) * op_ks_diag(ib,ikpt,isppol)
             else
               mat_temp(:,ib) = paw_dmft%chipsi(1:ndim,ib,ik,isppol,iatom) * oper%ks(ib,ib,ikpt,isppol)
             end if ! present(op_ks_diag)
           end do ! ib

         end if ! opt=1 or 3

         call abi_xgemm("n","c",ndim,ndim,mbandc,cone,mat_temp(:,:),ndim,&
                      & paw_dmft%chipsi(:,:,ik,isppol,iatom),ndim_max,czero,mat_temp2(:,:),ndim)

       else if (opt == 2) then

         call abi_xgemm("n","c",ndim,ndim,mbandc,cone,paw_dmft%chipsi(:,:,ik,isppol,iatom),&
                      & ndim_max,paw_dmft%chipsi(:,:,ik,isppol,iatom),ndim_max,czero,mat_temp2(:,:),ndim)

       else if (opt == 4) then

         call abi_xgemm("n","c",ndim,ndim,mbandc,cone,paw_dmft%chipsi(:,:,ik,isppol,iatom),&
                      & ndim_max,paw_dmft%chipsi(:,:,ik,isppol,iatom),ndim_max,czero,mat_temp3(:,:),ndim)

         call abi_xgemm("n","n",ndim,ndim,ndim,cone,mat_temp3(:,:),ndim,&
                      & mat_temp3(:,:),ndim,czero,mat_temp2(:,:),ndim)

       end if ! opt

       oper%matlu(iatom)%mat(:,:,isppol) = oper%matlu(iatom)%mat(:,:,isppol) + mat_temp2(:,:)*oper%wtk(ik)

     end do ! ikpt
   end do ! isppol
   ABI_FREE(mat_temp)
   ABI_FREE(mat_temp2)
   ABI_FREE(mat_temp3)
 end do ! iatom

!do isppol=1,nsppol
 ! do ikpt=1,nkpt
 !  ikpt1=ikpt
 !  if(present(jkpt)) ikpt1=jkpt
 !  lvz=paral==0  !vz_d
 !  if(present(iproc)) lvz=lvz.or.(paral==1.and.(procb2(ikpt1)==iproc))  !vz_d
!!  if ((paral==1.and.(procb2(ikpt1)==iproc)).or.(paral==0)) then    !vz_d
 !  if(lvz) then !vz_d
 !  do ib=1,mbandc
 !   do ib1=1,mbandc
 !    do iatom=1,natom
 !     if(oper%matlu(iatom)%lpawu.ne.-1) then
 !     ndim=2*oper%matlu(iatom)%lpawu+1
 !      do im=1,ndim
 !       do im1=1,ndim
 !        do ispinor=1,nspinor
 !         do ispinor1=1,nspinor
 !           if (im1 == im .and. im1 == 1 .and. ib1 == ib .and. iatom == 1) then

  !          end if
  !          oper%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=     &
!&            oper%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)+    &
!&            paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)*        &
!&            conjg(paw_dmft%psichi(isppol,ikpt1,ib1,ispinor1,iatom,im1))* &
!false&            paw_dmft%psichi(isppol,ikpt1,ib1,ispinor1,iatom,im1)*
!&
!false&            conjg(paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im))* &
!&            oper%ks(isppol,ikpt,ib,ib1)*oper%wtk(ikpt)
! one  could suppress wtk here if present(jkpt)
! ks(ib,ib1)=ks(ib1,ib) -> ib and ib1 can be underchanged !
!          enddo ! ispinor1
!         enddo ! ispinor
!        enddo ! im1
!       enddo ! im
!      endif
!     enddo ! iatom
 !   enddo ! ib
 !  enddo ! ib
 !  endif
 ! enddo ! ikpt
 !enddo ! isppol

 DBG_EXIT("COLL")

end subroutine downfold_oper
!!***

!!****f* m_oper/upfold_oper
!! NAME
!! upfold_oper
!!
!! FUNCTION
!!  Upfold an operator from local space to KS space
!!
!! INPUTS
!!  oper <type(oper_type)>= operator
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  procb(ikpt)=for kpt parallelization; gives the rank (in the kpt communicator) of the CPU handling each ikpt
!!  iproc=rank of the current process in the kpt communicator
!!
!! OUTPUT
!!
!! SOURCE

subroutine upfold_oper(oper,paw_dmft,procb,iproc)

 use m_paw_dmft, only : paw_dmft_type
 use m_abi_linalg, only : abi_xgemm

!Arguments ------------------------------------
 type(oper_type), intent(inout)  :: oper
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, optional, intent(in)   :: iproc
 integer, optional, intent(in)   :: procb(oper%nkpt)
!Local variables-------------------------------
 integer :: iatom,ik,ikpt,isppol,lpawu,mbandc
 integer :: ndim,ndim_max,nspinor,paral,shift
 complex(dpc), allocatable :: mat_temp(:,:),mat_temp2(:,:)
! *********************************************************************

!   write(6,*) "upfold_oper procb",procb
!   write(6,*) "iproc",iproc
!   write(6,*) size(procb)
!   write(6,*) size(procb2)
!   write(6,*) procb2(1),procb2(16)

 DBG_ENTER("COLL")

 !if ((oper%has_opermatlu == 0) .or. (oper%has_operks == 0)) then
 !  message = " Operator is not defined to be used in upfold_oper"
 !  ABI_ERROR(message)
 !end if

 mbandc   = paw_dmft%mbandc
 nspinor  = paw_dmft%nspinor
 ndim_max = nspinor * (2*paw_dmft%maxlpawu+1)
 paral    = 0
 shift    = oper%shiftk

 if (present(procb) .and. present(iproc) .and. oper%paral == 0) paral = 1

 oper%ks(:,:,:,:) = czero

 ABI_MALLOC(mat_temp,(mbandc,ndim_max))
 ABI_MALLOC(mat_temp2,(mbandc,mbandc))

 do iatom=1,oper%natom
   lpawu = oper%matlu(iatom)%lpawu
   if (lpawu == -1) cycle
   ndim = (2*lpawu+1) * nspinor
   do isppol=1,oper%nsppol
     do ikpt=1,oper%nkpt ! index of kpt on the current CPU

       if (paral == 1) then
         if (procb(ikpt) /= iproc) cycle
       end if

       ik = ikpt + shift ! true kpt index (needed for chipsi)

       call abi_xgemm("c","n",mbandc,ndim,ndim,cone,paw_dmft%chipsi(:,:,ik,isppol,iatom),&
                    & ndim_max,oper%matlu(iatom)%mat(:,:,isppol),ndim,czero,mat_temp(:,1:ndim),mbandc)

       call abi_xgemm("n","n",mbandc,mbandc,ndim,cone,mat_temp(:,1:ndim),mbandc,&
                    & paw_dmft%chipsi(:,:,ik,isppol,iatom),ndim_max,czero,mat_temp2(:,:),mbandc)

       oper%ks(:,:,ikpt,isppol) = oper%ks(:,:,ikpt,isppol) + mat_temp2(:,:)

     end do ! ikpt
   end do ! isppol
 end do ! iatom

 ABI_FREE(mat_temp)
 ABI_FREE(mat_temp2)

!do isppol=1,nsppol
 !  do ikpt=1,nkpt
 !   if ((paral==1.and.(procb2(ikpt)==iproc)).or.(paral==0)) then
 !    do ib=1,mbandc
 !      do ib1=1,mbandc
!               if(ib==1.and.ib1==3) write(std_out,*) "IKPT=",ikpt
 !       oper%ks(isppol,ikpt,ib,ib1)=czero

  !       do iatom=1,natom
  !         if(oper%matlu(iatom)%lpawu.ne.-1) then
  !           ndim=2*oper%matlu(iatom)%lpawu+1
  !           do im=1,ndim
  !             do im1=1,ndim
  !               do ispinor=1,nspinor
  !                 do ispinor1=1,nspinor

! psichi(isppol,ikpt,ib,ispinor,iatom,im)=<\chi_{m,R,ispinor)|\Psi(s,k,nu)>
   !                  oper%ks(isppol,ikpt,ib,ib1)= oper%ks(isppol,ikpt,ib,ib1) &
!&                     + ( paw_dmft%psichi(isppol,ikpt,ib1,ispinor1,iatom,im1)
!&
!&                     * oper%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)
!&
!&                     *
!conjg(paw_dmft%psichi(isppol,ikpt,ib,ispinor,iatom,im)))
              ! if(present(prt).and.(ib==1.and.ib1==1)) then
              !   write(6,*) "im,im1",im,im1
              !   write(6,*) "ispinor,ispinor1",ispinor,ispinor1
              !   write(6,*)
              !   "psichi",paw_dmft%psichi(isppol,ikpt,ib1,ispinor1,iatom,im1)
              !   write(6,*) "psichi
              !   2",paw_dmft%psichi(isppol,ikpt,ib,ispinor,iatom,im1)
              !   write(6,*) "oper%matlu",
              !   oper%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)
              ! endif

   !                enddo ! ispinor1
   !              enddo ! ispinor
   !            enddo ! im1
   !          enddo ! im
   !        endif
   !      enddo ! iatom

   !    enddo ! ib
   !  enddo ! ib
   ! endif
   !enddo ! ikpt
 !enddo ! isppol

 DBG_EXIT("COLL")

end subroutine upfold_oper
!!***

!!****f* m_oper/identity_oper
!! NAME
!! identity_oper
!!
!! FUNCTION
!!  Construct the identity operator
!!
!! INPUTS
!!  oper <type(oper_type)>= operator
!!  option = 1: in KS space
!!         = 2: in local space
!!         = 3: both
!!
!! OUTPUT
!!
!! SOURCE

subroutine identity_oper(oper,option)

 use m_matlu, only : identity_matlu,zero_matlu

!Arguments ------------------------------------
 integer, intent(in) :: option
 type(oper_type), intent(inout) :: oper
!Local variables-------------------------------
 integer :: ib,ikpt,isppol,natom
 character(len=500) :: message
! *********************************************************************

 DBG_ENTER("COLL")

 if (((option == 1 .or. option == 3) .and. (oper%has_operks == 0)) .or. &
   & ((option == 2 .or. option == 3) .and. (oper%has_opermatlu == 0))) then
   message = " Options in identity_oper are not coherent with definitions of this operator"
   ABI_ERROR(message)
 end if

 if (option == 1 .or. option == 3) then

   oper%ks(:,:,:,:) = czero
   do isppol=1,oper%nsppol
     do ikpt=1,oper%nkpt
       do ib=1,oper%mbandc
         oper%ks(ib,ib,ikpt,isppol) = cone
       end do ! ib
     end do ! ikpt
   end do ! isppol

 end if ! option=1 or 3

 if (option == 2 .or. option == 3) then
   natom = oper%natom
   call zero_matlu(oper%matlu(:),natom)
   call identity_matlu(oper%matlu(:),natom)
 end if ! option=2 or 3

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
!!  char1 = character describing occup1
!!  char2 = character describing occup2
!!  occup1 <type(oper_type)>= occupations
!!  occup2 <type(oper_type)>= occupations
!!  option : option for printing (if 1 assume data are related to lda only)
!!  toldiff : tolerance for the difference
!!
!! OUTPUT
!!
!! SOURCE

subroutine diff_oper(char1,char2,occup1,occup2,option,toldiff)

 use m_matlu, only : diff_matlu

!Arguments ------------------------------------
 type(oper_type), intent(in) :: occup1,occup2
 integer, intent(in) :: option
 real(dp), intent(in) :: toldiff
 character(len=*), intent(in) :: char1,char2
!Local variables-------------------------------
 character(len=500) :: message
! *********************************************************************

 DBG_ENTER("COLL")

 if (occup1%has_opermatlu == 0 .or. occup2%has_opermatlu == 0) then
   message = " Operators are not defined to be used in diff_oper"
   ABI_ERROR(message)
 end if

 if (occup1%nkpt /= occup2%nkpt) then
   write(message,'(a,2x,2i9)')' Operators are not equal',occup1%nkpt,occup2%nkpt
   ABI_ERROR(message)
 end if

 call diff_matlu(char1,char2,occup1%matlu(:),occup2%matlu(:),occup1%natom,option,toldiff)
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
!!  Computes the trace of an operator
!!
!! INPUTS
!!  oper <type(oper_type)>= operator
!!  opt_ksloc = 1: trace in KS space
!!            = 2: trace in local space
!!            = 3: both
!!
!! OUTPUT
!!  trace_ks  :: trace in KS space
!!  trace_loc :: trace in local space
!!  trace_ks_cmplx :: complex trace in KS space
!!
!! SOURCE

subroutine trace_oper(oper,trace_ks,trace_loc,opt_ksloc,trace_ks_cmplx)

 use m_matlu, only : trace_matlu

!Arguments ------------------------------------
 type(oper_type), intent(in) :: oper
 real(dp), intent(out) :: trace_ks  !vz_i
 real(dp), intent(inout) :: trace_loc(oper%nsppol+1,oper%natom) !vz_i
 integer, intent(in) :: opt_ksloc
 complex(dpc), optional, intent(out) :: trace_ks_cmplx
!Local variables-------------------------------
 integer :: ib,ikpt,isppol
 complex(dpc) :: trace
 character(len=500) :: message
! *********************************************************************

 DBG_ENTER("COLL")

 if (((opt_ksloc == 1 .or. opt_ksloc == 3) .and. (oper%has_operks == 0)) .or. &
   & ((opt_ksloc == 2 .or. opt_ksloc == 3) .and. (oper%has_opermatlu == 0))) then
   message = " Options in trace_oper are not coherent with definitions of this operator"
   ABI_ERROR(message)
 end if

 if (opt_ksloc == 1 .or. opt_ksloc == 3) then
   trace = czero
   !temp1=zero
   do isppol=1,oper%nsppol
     do ikpt=1,oper%nkpt
       do ib=oper%mbandc,1,-1 ! Never change this order
         trace = trace + oper%ks(ib,ib,ikpt,isppol)*oper%wtk(ikpt+oper%shiftk)
         !temp1=temp1+oper%wtk(ikpt)
       end do ! ib
     end do ! ikpt
   end do ! isppol
   if (oper%nsppol == 1 .and. oper%nspinor == 1) trace = two * trace
   if (present(trace_ks_cmplx)) trace_ks_cmplx = trace
   trace_ks = dble(trace)
!   write(std_out,*) "temp1",temp1
 end if ! opt_ksloc

 if (opt_ksloc == 2 .or. opt_ksloc == 3) then
   call trace_matlu(oper%matlu(:),oper%natom,trace_loc=trace_loc(:,:))
 end if

 DBG_EXIT("COLL")

end subroutine trace_oper
!!***

!!****f* m_oper/prod_oper
!! NAME
!! prod_oper
!!
!! FUNCTION
!!  Computes the matrix product of oper1 and oper2
!!
!! INPUTS
!!  oper1,oper2 <type(oper_type)>= operator
!!  opt_ksloc = 1 : in KS space
!!            = 2 : in local space
!!  opt_diag = 1 if oper1 and oper2 are diagonal in KS space, 0 otherwise (default)
!!
!! OUTPUT
!!  oper3 <type(oper_type)>= matrix product of oper1 and oper2
!!
!! SOURCE

subroutine prod_oper(oper1,oper2,oper3,opt_ksloc,opt_diag)

 use m_matlu, only : prod_matlu
 use m_abi_linalg, only : abi_xgemm

!Arguments ------------------------------------
 type(oper_type), intent(in) :: oper1,oper2
 type(oper_type), intent(inout) :: oper3
 integer, intent(in) :: opt_ksloc
 integer, optional, intent(in) :: opt_diag
!Local variables-------------------------------
 integer :: ib,ikpt,isppol,mbandc
 logical :: diag
! *********************************************************************

 DBG_ENTER("COLL")

 if (opt_ksloc == 2 .and. oper1%has_opermatlu == 1 .and. &
   & oper2%has_opermatlu == 1 .and. oper3%has_opermatlu == 1) then
   call prod_matlu(oper1%matlu(:),oper2%matlu(:),oper3%matlu(:),oper1%natom)
 end if

 if (opt_ksloc == 1 .and. oper1%has_operks == 1 .and. &
    & oper2%has_operks == 1 .and. oper3%has_operks == 1) then
   mbandc = oper1%mbandc
   diag   = .false.
   if (present(opt_diag)) then
     if (opt_diag == 1) diag = .true.
   end if
   do isppol=1,oper1%nsppol
     do ikpt=1,oper1%nkpt
       if (diag) then
         do ib=1,mbandc
           oper3%ks(ib,ib,ikpt,isppol) = oper1%ks(ib,ib,ikpt,isppol) * oper2%ks(ib,ib,ikpt,isppol)
         end do ! ib
       else
         call abi_xgemm("n","n",mbandc,mbandc,mbandc,cone,oper1%ks(:,:,ikpt,isppol),mbandc,&
                      & oper2%ks(:,:,ikpt,isppol),mbandc,czero,oper3%ks(:,:,ikpt,isppol),mbandc)
         end if ! diag
       end do ! ikpt
     end do ! isppol
 end if ! opt_ksloc=1

 DBG_EXIT("COLL")

end subroutine prod_oper
!!***

!!****f* m_oper/trace_prod_oper
!! NAME
!! trace_prod_oper
!!
!! FUNCTION
!!  Computes Tr(oper1*oper2) in KS space
!!
!! INPUTS
!!  oper1,oper2 <type(oper_type)>= operator
!!
!! OUTPUT
!!  trace = Tr(op1*op2)
!!
!! SOURCE

subroutine trace_prod_oper(oper1,oper2,trace)

!Arguments ------------------------------------
 type(oper_type), intent(in) :: oper1,oper2
 complex(dpc), intent(out) :: trace
!Local variables-------------------------------
 integer :: ikpt,isppol
 character(len=500) :: message
! *********************************************************************

 if (oper1%shiftk /= oper2%shiftk) then
   message = "Inconsistency in trace_prod_oper: oper1%shiftk should be equal to oper2%shiftk"
   ABI_ERROR(message)
 end if

 trace = czero

 do isppol=1,oper1%nsppol
   do ikpt=1,oper1%nkpt
     trace = trace + sum(oper1%ks(:,:,ikpt,isppol)*transpose(oper2%ks(:,:,ikpt,isppol)))*&
       & oper1%wtk(ikpt+oper1%shiftk)
   end do ! ikpt
 end do ! isppol

 if (oper1%nsppol == 1 .and. oper1%nspinor == 1) trace = trace * two

end subroutine trace_prod_oper
!!***

!!****f* m_oper/gather_oper
!! NAME
!! gather_oper
!!
!! FUNCTION
!!  Gather the contributions from all CPUs, for a frequency-dependent
!!  operator, and for both levels of parallelization (kpt and then frequency,
!!  and frequency only).
!!
!! INPUTS
!!  oper <type(oper_type)>= operator for each frequency
!!  distrib <type(mpi_distrib_dmft_type)> = mpi related data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  opt_ksloc = 1 : gather the KS operator on the kpt and frequency communicator
!!              2 : gather the local operator (the exact behavior can be defined via opt_commkpt)
!!  master = if present, only gather on the master node
!!  opt_diag = 1 if the operator is diagonal in KS space, 0 (default) otherwise
!!  opt_commkpt (only meaningful for the local quantity)
!!              = 0 (default) : frequency-only parallelization
!!                              -> xmpi_allgatherv on the whole communicator
!!              = 1 : kpt and then frequency parallelization (CAREFUL: here, frequencies
!!                     are not distributed in the same way as the frequency-only parallelization scheme)
!!                    -> xmpi_sum on the kpt-communicator, and then xmpi_allgatherv
!!                      on the frequency communicator (useful after a downfold for instance)
!!
!! OUTPUT
!!
!! SOURCE

subroutine gather_oper(oper,distrib,paw_dmft,opt_ksloc,master,opt_diag,opt_commkpt)

 use m_paw_dmft, only : mpi_distrib_dmft_type,paw_dmft_type
 use m_xmpi, only : xmpi_allgatherv,xmpi_gatherv,xmpi_sum,xmpi_sum_master

!Arguments ------------------------------------
 type(mpi_distrib_dmft_type), target, intent(in) :: distrib
 type(oper_type), intent(inout) :: oper(distrib%nw)
 type(paw_dmft_type) :: paw_dmft
 integer, intent(in) :: opt_ksloc
 integer, optional, intent(in) :: master,opt_commkpt,opt_diag
!Local variables-------------------------------
 integer :: comm,iatom,ib,ib1,ibuf,ierr,ifreq,ikpt,im,im1,irank,irank1,irank2,isppol,lpawu,mbandc
 integer :: myproc,myproc2,natom,ndim,nkpt,nproc,nproc_freq,nproc_kpt,nproc2,nspinor,nsppol,nw,optcommkpt,siz_buf
 logical :: diag
 integer, allocatable :: displs(:),recvcounts(:)
 complex(dpc), allocatable :: buffer(:),buffer_tot(:)
! *********************************************************************

 comm    = paw_dmft%spacecomm
 mbandc  = paw_dmft%mbandc
 myproc  = paw_dmft%myproc
 natom   = paw_dmft%natom
 nkpt    = paw_dmft%nkpt
 nproc   = paw_dmft%nproc
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol
 nw      = distrib%nw

 nproc_kpt  = min(nkpt,nproc)

 optcommkpt = 0
 if (present(opt_commkpt)) optcommkpt = opt_commkpt

 if (opt_ksloc == 1) then

   nproc_freq = max(1,nproc/nkpt)

   ABI_MALLOC(recvcounts,(nproc))
   ABI_MALLOC(displs,(nproc))

   diag = .false.
   if (present(opt_diag)) then
     if (opt_diag == 1) diag = .true.
   end if

   irank2 = 1
   do irank=0,nproc_kpt-1
     do irank1=0,nproc_freq-1
       recvcounts(irank2) = distrib%nkpt_mem(irank+1) * distrib%nw_mem_kptparal(irank1+1)
       irank2 = irank2 + 1
     end do ! irank1
   end do ! irank
   if (nproc > nproc_freq*nproc_kpt) recvcounts(nproc_freq*nproc_kpt+1:nproc) = 0

   recvcounts(:) = mbandc * recvcounts(:)
   if (.not. diag) recvcounts(:) = recvcounts(:) * mbandc
   displs(1) = 0
   do irank=2,nproc
     displs(irank) = displs(irank-1) + recvcounts(irank-1)
   end do ! irank

   ABI_MALLOC(buffer,(recvcounts(myproc+1)))
   ABI_MALLOC(buffer_tot,(displs(nproc)+recvcounts(nproc)))

   do isppol=1,nsppol

     ibuf = 0
     do ikpt=1,nkpt
       if (distrib%procb(ikpt) /= distrib%me_kpt) cycle
       do ifreq=1,nw
         if (distrib%proct(ifreq) /= distrib%me_freq) cycle
         do ib1=1,mbandc
           if (diag) then
             ibuf = ibuf + 1
             buffer(ibuf) = oper(ifreq)%ks(ib1,ib1,ikpt,isppol)
           else
             do ib=1,mbandc
               ibuf = ibuf + 1
               buffer(ibuf) = oper(ifreq)%ks(ib,ib1,ikpt,isppol)
             end do ! ib
           end if ! diag
         end do ! ib1
       end do ! ifreq
     end do ! ikpt

     if (present(master)) then
       call xmpi_gatherv(buffer(:),recvcounts(myproc+1),buffer_tot(:),recvcounts(:),displs(:),master,comm,ierr)
     else
       call xmpi_allgatherv(buffer(:),recvcounts(myproc+1),buffer_tot(:),recvcounts(:),displs(:),comm,ierr)
     end if  ! present(master)

     ibuf = 0
     do ikpt=1,nkpt
       do ifreq=1,nw
         do ib1=1,mbandc
           if (diag) then
             ibuf = ibuf + 1
             oper(ifreq)%ks(ib1,ib1,ikpt,isppol) = buffer_tot(ibuf)
           else
             do ib=1,mbandc
               ibuf = ibuf + 1
               oper(ifreq)%ks(ib,ib1,ikpt,isppol) = buffer_tot(ibuf)
             end do ! ib
           end if ! diag
         end do ! ib1
       end do ! ifreq
     end do ! ikpt

   end do ! isppol

 else if (opt_ksloc == 2) then

   nproc_freq = nproc / nkpt

   if (optcommkpt == 1) then
     myproc2 = distrib%me_freq
     nproc2  = nproc_freq + 1
   else
     myproc2 = myproc
     nproc2  = nproc
   end if

   ABI_MALLOC(recvcounts,(nproc2))
   ABI_MALLOC(displs,(nproc2))

   siz_buf = 0

   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     siz_buf = siz_buf + (2*lpawu+1)**2
   end do ! iatom

   siz_buf = siz_buf * (nspinor**2) * nsppol
   if (optcommkpt == 1) then
     recvcounts(:) = siz_buf * distrib%nw_mem_kptparal(:)
   else if (optcommkpt == 0) then
     recvcounts(:) = siz_buf * distrib%nw_mem(:)
   end if ! optcommkpt
   displs(1) = 0
   do irank=2,nproc2
     displs(irank) = displs(irank-1) + recvcounts(irank-1)
   end do ! irank

   if (optcommkpt == 1 .and. recvcounts(myproc2+1) == 0) then
     if (nproc_freq > 1) then
       siz_buf = siz_buf * distrib%nw_mem_kptparal(mod(paw_dmft%myproc,nproc_freq)+1)
     else
       siz_buf = siz_buf * nw
     end if
   else
     siz_buf = recvcounts(myproc2+1)
   end if

   ABI_MALLOC(buffer,(siz_buf))
   ABI_MALLOC(buffer_tot,(recvcounts(nproc2)+displs(nproc2)))

   buffer(:) = czero

   ibuf = 0
   do ifreq=1,nw
     if (optcommkpt == 1) then
       if (distrib%proct(ifreq) /= distrib%me_freq) cycle
     else if (optcommkpt == 0) then
       if (distrib%procf(ifreq) /= myproc) cycle
     end if ! optcommkpt
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       ndim = (2*lpawu+1) * nspinor
       do isppol=1,nsppol
         do im1=1,ndim
           do im=1,ndim
             ibuf = ibuf + 1
             buffer(ibuf) = oper(ifreq)%matlu(iatom)%mat(im,im1,isppol)
           end do ! im
         end do ! im1
       end do ! isppol
     end do ! iatom
   end do ! ifreq

   if (optcommkpt == 1) comm = distrib%comm_freq
   if (present(master)) then
     if (optcommkpt == 1) then
       call xmpi_sum_master(buffer(:),master,distrib%comm_kpt,ierr)
     end if
     call xmpi_gatherv(buffer(:),recvcounts(myproc2+1),buffer_tot(:),recvcounts(:),displs(:),master,comm,ierr)
   else
     if (optcommkpt == 1) then
       call xmpi_sum(buffer(:),distrib%comm_kpt,ierr)
     end if
     call xmpi_allgatherv(buffer(:),recvcounts(myproc2+1),buffer_tot(:),recvcounts(:),displs(:),comm,ierr)
   end if  ! present(master)

   ibuf = 0
   do ifreq=1,nw
     do iatom=1,natom
       lpawu = paw_dmft%lpawu(iatom)
       if (lpawu == -1) cycle
       ndim = (2*lpawu+1) * nspinor
       do isppol=1,nsppol
         do im1=1,ndim
           do im=1,ndim
             ibuf = ibuf + 1
             oper(ifreq)%matlu(iatom)%mat(im,im1,isppol) = buffer_tot(ibuf)
           end do ! im
         end do ! im1
       end do ! isppol
     end do ! iatom
   end do ! ifreq

 end if ! opt_ksloc

 ABI_FREE(recvcounts)
 ABI_FREE(displs)
 ABI_FREE(buffer)
 ABI_FREE(buffer_tot)

end subroutine gather_oper
!!***

!!****f* m_oper/gather_oper_ks
!! NAME
!! gather_oper_ks
!!
!! FUNCTION
!!  For a single KS operator, performs a xmpi_sum on the frequency communicator,
!!  and a xmpi_allgatherv on the kpt-communicator
!!
!! INPUTS
!!  oper <type(oper_type)>= operator
!!  distrib <type(mpi_distrib_dmft_type)> = mpi related data
!!  opt_diag = 1 if the operator is diagonal in KS space, 0 (default) otherwise
!!
!! OUTPUT
!!
!! SOURCE

subroutine gather_oper_ks(oper,distrib,paw_dmft,opt_diag)

 use m_paw_dmft, only : mpi_distrib_dmft_type,paw_dmft_type
 use m_xmpi, only : xmpi_allgatherv,xmpi_sum

!Arguments ------------------------------------
 type(oper_type), intent(inout) :: oper
 type(mpi_distrib_dmft_type), intent(in) :: distrib
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, optional, intent(in) :: opt_diag
!Local variables-------------------------------
 integer :: ib,ib1,ibuf,ierr,ikpt,irank,isppol,mbandc
 integer :: me_kpt,nkpt,nproc,nproc_freq,nsppol,nw,siz_buf
 logical :: diag
 integer, allocatable :: displs(:),recvcounts(:)
 complex(dpc), allocatable :: buffer(:),buffer_tot(:)
! *********************************************************************

 mbandc = paw_dmft%mbandc
 me_kpt = distrib%me_kpt
 nkpt   = paw_dmft%nkpt
 nproc  = paw_dmft%nproc
 nsppol = paw_dmft%nsppol
 nw     = distrib%nw

 nproc_freq = nproc / nkpt

 diag = .false.
 if (present(opt_diag)) then
   if (opt_diag == 1) diag = .true.
 end if

 ABI_MALLOC(recvcounts,(nproc))
 ABI_MALLOC(displs,(nproc))

 recvcounts(:) = mbandc * distrib%nkpt_mem(:)
 if (.not. diag) recvcounts(:) = recvcounts(:) * mbandc

 displs(1) = 0
 do irank=2,nproc
   displs(irank) = displs(irank-1) + recvcounts(irank-1)
 end do ! irank

 siz_buf = recvcounts(me_kpt+1)
 if (siz_buf == 0) then
   siz_buf = mbandc
   if (nproc_freq <= 1) siz_buf = siz_buf * distrib%nkpt_mem(mod(paw_dmft%myproc,min(nkpt,nproc))+1)
   if (.not. diag) siz_buf = siz_buf * mbandc
 end if
 ABI_MALLOC(buffer,(siz_buf))
 ABI_MALLOC(buffer_tot,(recvcounts(nproc)+displs(nproc)))

 do isppol=1,nsppol

   buffer(:) = czero

   ibuf = 0
   do ikpt=1,nkpt
     if (distrib%procb(ikpt) /= me_kpt) cycle
     do ib1=1,mbandc
       if (diag) then
         ibuf = ibuf + 1
         buffer(ibuf) = oper%ks(ib1,ib1,ikpt,isppol)
       else
         do ib=1,mbandc
           ibuf = ibuf + 1
           buffer(ibuf) = oper%ks(ib,ib1,ikpt,isppol)
         end do ! ib
       end if ! diag
     end do ! ib1
   end do ! ikpt

   call xmpi_sum(buffer(:),distrib%comm_freq,ierr)

   call xmpi_allgatherv(buffer(:),recvcounts(me_kpt+1),&
      & buffer_tot(:),recvcounts(:),displs(:),distrib%comm_kpt,ierr)

   ibuf = 0
   do ikpt=1,nkpt
     do ib1=1,mbandc
       if (diag) then
         ibuf = ibuf + 1
         oper%ks(ib1,ib1,ikpt,isppol) = buffer_tot(ibuf)
       else
         do ib=1,mbandc
           ibuf = ibuf + 1
           oper%ks(ib,ib1,ikpt,isppol) = buffer_tot(ibuf)
         end do ! ib
       end if ! diag
     end do ! ib
   end do ! ikpt

 end do ! isppol

 ABI_FREE(buffer)
 ABI_FREE(buffer_tot)
 ABI_FREE(recvcounts)
 ABI_FREE(displs)

end subroutine gather_oper_ks
!!***

END MODULE m_oper
!!***
