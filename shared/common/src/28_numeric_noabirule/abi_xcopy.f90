!{\src2tex{textfont=tt}}
!!****f* m_abi_linalg/abi_xcopy
!! NAME
!!  abi_xcopy
!!
!! FUNCTION
!! abi_xcopy is the generic function for vector copy
!! It performs the data copy: dst(1:n:incdst) = src(1:n:incsrc)
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2019 ABINIT group (LNguyen,FDahm (CS))
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/Infos/copyright
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

!!***

!!****f* m_abi_linalg/abi_zcopy
!! NAME
!! abi_zcopy
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE
!!
  subroutine abi_zcopy(size,tsrc,incsrc,tdest,incdest)

 implicit none

!Arguments-------------------------------------
 integer, intent(in) :: size
 integer, intent(in) :: incsrc
 integer, intent(in) :: incdest
 complex(dpc),intent(in) :: tsrc
 complex(dpc),intent(inout) :: tdest

 !Local variables-------------------------------
#ifdef DEV_LINALG_TIMING
 real(dp) :: tsec(2)
 call timab(TIMAB_XCOPY,1,tsec)
#endif

 call zcopy(size,tsrc,incsrc,tdest,incdest)

#ifdef DEV_LINALG_TIMING
 call timab(TIMAB_XCOPY,2,tsec)
#endif

end subroutine abi_zcopy
!!***

!!****f* m_abi_linalg/abi_zcopy_1d
!! NAME
!! abi_zcopy_1d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

subroutine abi_zcopy_1d(size,tsrc,incsrc,tdest,incdest)

 implicit none

!Arguments-------------------------------------
 integer, intent(in) :: size
 integer, intent(in) :: incsrc
 integer, intent(in) :: incdest
 complex(dpc), intent(in) :: tsrc(*)
 complex(dpc), intent(inout) :: tdest(*)

 call abi_zcopy(size,tsrc(1),incsrc,tdest(1),incdest)

end subroutine abi_zcopy_1d
!!***

!!****f* m_abi_linalg/abi_dcopy
!! NAME
!! abi_dcopy
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

  subroutine abi_dcopy(size,tsrc,incsrc,tdest,incdest,x_cplx)

 implicit none

!Arguments-------------------------------------
 integer, intent(in) :: size
 integer, intent(in) :: incsrc
 integer, intent(in) :: incdest
 real(dp),intent(in) :: tsrc
 real(dp),intent(inout) :: tdest
 integer, intent(in), optional :: x_cplx

!Local variables-------------------------------
 integer  :: cplx_

#ifdef DEV_LINALG_TIMING
 real(dp) :: tsec(2)
 call timab(TIMAB_XCOPY,1,tsec)
#endif

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx
 if(cplx_ == 2) then
    call zcopy(size,tsrc,incsrc,tdest,incdest)
 else
    call dcopy(size,tsrc,incsrc,tdest,incdest)
 end if

#ifdef DEV_LINALG_TIMING
 call timab(TIMAB_XCOPY,2,tsec)
#endif

end subroutine abi_dcopy
!!***

!!****f* m_abi_linalg/abi_dcopy_1d
!! NAME
!! abi_dcopy_1d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

subroutine abi_dcopy_1d(size,tsrc,incsrc,tdest,incdest,x_cplx)

 implicit none

!Arguments-------------------------------------
 integer, intent(in) :: size
 integer, intent(in) :: incsrc
 integer, intent(in) :: incdest
 real(dp), intent(in)  :: tsrc(*)
 real(dp), intent(inout) :: tdest(*)
!Only for lobpcgwf
 integer, intent(in), optional :: x_cplx

 call abi_dcopy(size,tsrc(1),incsrc,tdest(1),incdest,x_cplx)

end subroutine abi_dcopy_1d
!!***

!!****f* m_abi_linalg/abi_dcopy_2d
!! NAME
!! abi_dcopy_2d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

subroutine abi_dcopy_2d(size,tsrc,incsrc,tdest,incdest,x_cplx)

 implicit none

!Arguments-------------------------------------
 integer, intent(in) :: size
 integer, intent(in) :: incsrc
 integer, intent(in) :: incdest
 real(dp), DEV_CONTARRD intent(in) :: tsrc(:,:)
 real(dp), DEV_CONTARRD intent(inout) :: tdest(:,:)
!Only for lobpcgwf
 integer, intent(in), optional :: x_cplx

 ! write(*,*) "dcopy2D size=",size
 call abi_dcopy(size,tsrc(1,1),incsrc,tdest(1,1),incdest,x_cplx)

end subroutine abi_dcopy_2d
!!***

!!****f* m_abi_linalg/abi_dcopy_0d_1d
!! NAME
!! abi_dcopy_0d_1d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

subroutine abi_dcopy_0d_1d(size,tsrc,incsrc,tdest,incdest,x_cplx)

 implicit none

!Arguments-------------------------------------
 integer, intent(in) :: size
 integer, intent(in) :: incsrc
 integer, intent(in) :: incdest
 real(dp),intent(in) :: tsrc
 real(dp), intent(inout) :: tdest(*)
 integer,intent(in), optional :: x_cplx !only lobpcgwf

 call abi_dcopy(size,tsrc,incsrc,tdest(1),incdest,x_cplx)

end subroutine abi_dcopy_0d_1d
!!***

!!****f* m_abi_linalg/abi_dcopy_1d_0d
!! NAME
!! abi_dcopy_1d_0d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

subroutine abi_dcopy_1d_0d(size,tsrc,incsrc,tdest,incdest,x_cplx)


 implicit none
 !Arguments-------------------------------------
 integer, intent(in) :: size
 integer, intent(in) :: incsrc
 integer, intent(in) :: incdest
 real(dp),intent(in) :: tsrc(*)
 real(dp),intent(inout) :: tdest
 integer,intent(in), optional :: x_cplx !only lobpcgwf

 call abi_dcopy(size,tsrc(1),incsrc,tdest,incdest,x_cplx)

end subroutine abi_dcopy_1d_0d
!!***

!!****f* m_abi_linalg/abi_d2zcopy_2d
!! NAME
!! abi_d2zcopy_2d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

  subroutine abi_d2zcopy_2d(size,tsrc,incsrc,tdest,incdest,x_cplx)


 implicit none

 !Arguments-------------------------------------
 integer, intent(in) :: size
 integer, intent(in) :: incsrc
 integer, intent(in) :: incdest
 real(dp),     DEV_CONTARRD intent(in) :: tsrc(:,:)
 complex(dpc), DEV_CONTARRD intent(inout) :: tdest(:,:)
 !only in lobpcgwf
 integer, intent(in),optional :: x_cplx

!Local variables-------------------------------
 integer  :: cplx_

#ifdef DEV_LINALG_TIMING
 real(dp) :: tsec(2)
 call timab(TIMAB_XCOPY,1,tsec)
#endif

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx
 if(cplx_ == 2) then
    call zcopy(size,tsrc,incsrc,tdest,incdest)
 else
    call dcopy(size,tsrc,incsrc,tdest,incdest)
 end if

#ifdef DEV_LINALG_TIMING
 call timab(TIMAB_XCOPY,2,tsec)
#endif

end subroutine abi_d2zcopy_2d
!!***

!!****f* m_abi_linalg/abi_z2dcopy_2d
!! NAME
!! abi_z2dcopy_2d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! SOURCE

subroutine abi_z2dcopy_2d(size,tsrc,incsrc,tdest,incdest,x_cplx)

 implicit none

!Arguments-------------------------------------
 integer, intent(in) :: size
 integer, intent(in) :: incsrc
 integer, intent(in) :: incdest
 complex(dpc), DEV_CONTARRD intent(in) :: tsrc(:,:)
 real(dp),     DEV_CONTARRD intent(inout) :: tdest(:,:)
 !only in lobpcgwf
 integer,intent(in), optional :: x_cplx

!Local variables-------------------------------
 integer  :: cplx_
#ifdef DEV_LINALG_TIMING
 real(dp) :: tsec(2)
 call timab(TIMAB_XCOPY,1,tsec)
#endif

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx
 if(cplx_ == 2) then
    call zcopy(size,tsrc,incsrc,tdest,incdest)
 else
    call dcopy(size,tsrc,incsrc,tdest,incdest)
 end if

#ifdef DEV_LINALG_TIMING
 call timab(TIMAB_XCOPY,2,tsec)
#endif

end subroutine abi_z2dcopy_2d
!!***
