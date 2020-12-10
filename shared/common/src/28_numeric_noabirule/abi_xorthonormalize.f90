!{\src2tex{textfont=tt}}
!!****f* m_abi_linalg/abi_xorthonormalize
!! NAME
!!  abi_xorthonormalize
!!
!! FUNCTION
!!  abi_xorthonormalize is the generic function for computing the
!!  overlap of two complex wavefunctions (for a given number of bands)
!!  and orthonormalizes it:
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2020 ABINIT group (LNguyen,FDahm (CS), FBottin, GZ, AR, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/Infos/copyright
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

!!***

!!****f* m_abi_linalg/xorthonormalize
!! NAME
!!  xorthonormalize
!!
!! FUNCTION
!!  This routine computes the overlap of two complex wavefunctions (for a given number of bands)
!!  and orthonormalizes it:
!!      - Computes the products of two rectangular matrices
!!         containing the wavefunctions psi and S.psi (where S is the
!!         overlap (with the PAW terms if necessary)).
!!      - Does a Cholesky decomposition of this overlap
!!      - rotates the initial matrix blockvectorx by the triangular matrix to
!!         have an orthonormal set of wavefunctions
!!
!! INPUTS
!!  blockvectorbx = matrix of dimension (blocksize,vectsize)
!!                  (e.g. block of overlap*wavefunction)
!!  blocksize     = dimension of matrices (e.g number of bands)
!!  spaceComm     = communicator used for MPI parallelization
!!  vectsize      = dimension of matrices (e.g number of G vector)
!!
!! OUTPUT
!!  sqgram        = Choleski decomposition of transpose(blockvector)*blockvectorx
!!
!! SIDE EFFECTS
!!  blockvectorx  = on input, matrix of dimension (vectsize,blocksize)
!!                  (e.g block of wavefunction)
!!  blockvectorx  = on output, orthonormalized wavefunction.
!!
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine xorthonormalize(blockvectorx,blockvectorbx,blocksize,spaceComm,sqgram,vectsize,&
&                          x_cplx,timopt,tim_xortho) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,vectsize,spaceComm,x_cplx
 integer, intent(in), optional :: timopt,tim_xortho
 !arrays
 real(dp),intent(in) :: blockvectorbx(vectsize,blocksize)
 real(dp),intent(inout) :: blockvectorx(vectsize,blocksize)
 real(dp),intent(out) :: sqgram(x_cplx*blocksize,blocksize)

!Local variables-------------------------------
 real(dp) :: tsec(2)
 integer  :: ierr,info
 character(len=500) :: message
 character, dimension(2) :: cparam

 ! *********************************************************************

 if (present(tim_xortho).and.present(timopt)) then
   if(abs(timopt)==3) then
     call timab(tim_xortho,1,tsec)
   end if
 end if

 cparam(1)='t'
 cparam(2)='c'

 call abi_xgemm(cparam(x_cplx),'n',blocksize,blocksize,vectsize,cone,blockvectorx,&
&   vectsize,blockvectorbx,vectsize,czero,sqgram,blocksize,x_cplx=x_cplx)

 call xmpi_sum(sqgram,spaceComm,ierr)

 !Cholesky factorization of sqgram (ouside upper Triangular of sqgram)
 call abi_xpotrf('u',blocksize,sqgram,blocksize,info,x_cplx=x_cplx)

 if (info /= 0 )  then
   write(message,'(a,i0)')'abi_xpotrf, info=',info
   ABI_ERROR(message)
 end if

 !Find X  X*sqgram=blockvectorx
 call abi_xtrsm('r','u','n','n',vectsize,blocksize,cone,sqgram,blocksize,&
&   blockvectorx,vectsize,x_cplx=x_cplx)

 if (present(tim_xortho).and.present(timopt)) then
   if(abs(timopt)==3) then
     call timab(tim_xortho,2,tsec)
   end if
 end if

end subroutine xorthonormalize
!!***

!!****f* ABINIT/ortho_reim
!! NAME
!! ortho_reim
!!
!! FUNCTION
!! This routine computes the overlap of two wavefunctions (for a given number of bands)
!! and orthonormalizes it:
!!      - Computes the products of two rectangular matrices
!!         containing the wavefunctions psi and S.psi (where S is the
!!         overlap (with the PAW terms if necessary)).
!!      - Does a Cholesky decomposition of this overlap
!!      - rotates the initial matrix blockvectorx by the triangular matrix to
!!         have an orthonormal set of wavefunctions
!!
!! This version operates on arrays in which the real and the imaginary part
!! are packed together (real parts first, them imaginary parts), used when istwfk=2
!!
!! INPUTS
!!  blockvectorbx = matrix of dimension (blocksize,vectsize)
!!                  (e.g. block of overlap*wavefunction)
!!  blocksize     = dimension of matrices (e.g number of bands)
!!  spaceComm     = communicator used for MPI parallelization
!!  vectsize      = dimension of matrices (e.g number of G vector)
!!
!! OUTPUT
!!  sqgram        = Choleski decomposition of transpose(blockvector)*blockvectorx
!!
!! SIDE EFFECTS
!!  blockvectorx  = on input, matrix of dimension (vectsize,blocksize)
!!                  (e.g block of wavefunction)
!!  blockvectorx  = on output, orthonormalized wavefunction.
!!
!! PARENTS
!!      lobpcgIIwf,m_lobpcg,m_lobpcgIIIwf,pw_orthon
!!
!! CHILDREN
!!      abi_xgemm,abi_xpotrf,abi_xtrsm,wrtout,xmpi_sum
!!
!! SOURCE

subroutine ortho_reim(blockvectorx,blockvectorbx,blocksize,spaceComm,sqgram,vectsize)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,vectsize,spaceComm
!arrays
 real(dp),intent(in) :: blockvectorbx(vectsize,blocksize)
 real(dp),intent(inout) :: blockvectorx(vectsize,blocksize)
 real(dp),intent(out) :: sqgram(blocksize,blocksize)

!Local variables-------------------------------
!scalars
 integer :: ierr,info
 character(len=500) :: message

! *********************************************************************

 call abi_xgemm('t','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
&   vectsize,blockvectorbx,vectsize,czero,sqgram,blocksize)

 call xmpi_sum(sqgram,spaceComm,ierr)

 !Cholesky factorization of sqgram (ouside upper Triangular of sqgram)
 call abi_d2zpotrf('u',blocksize,sqgram,blocksize,info) !vz_d

 if (info /= 0 )  then
   write(message,'(a,i0)')'dpotrf, info=',info
   ABI_ERROR(message)
 end if

!Find X  X*sqgram=blockvectorx
 call abi_xtrsm('r','u','n','n',vectsize,blocksize,one,sqgram,blocksize,blockvectorx,vectsize)

end subroutine ortho_reim
!!***


!!****f* ABINIT/zorthonormalize
!! NAME
!! zorthonormalize
!!
!! FUNCTION
!! This routine computes the overlap of two complex wavefunctions (for a given number of bands)
!! and orthonormalizes it:
!!      - Computes the products of two rectangular matrices
!!         containing the wavefunctions psi and S.psi (where S is the
!!         overlap (with the PAW terms if necessary)).
!!      - Does a Cholesky decomposition of this overlap
!!      - rotates the initial matrix blockvectorx by the triangular matrix to
!!         have an orthonormal set of wavefunctions
!!
!! INPUTS
!!  blockvectorbx = matrix of dimension (blocksize,vectsize)
!!                  (e.g. block of overlap*wavefunction)
!!  blocksize     = dimension of matrices (e.g number of bands)
!!  spaceComm     = communicator used for MPI parallelization
!!  vectsize      = dimension of matrices (e.g number of G vector)
!!
!! OUTPUT
!!  sqgram        = Choleski decomposition of transpose(blockvector)*blockvectorx
!!
!! SIDE EFFECTS
!!  blockvectorx  = on input, matrix of dimension (vectsize,blocksize)
!!                  (e.g block of wavefunction)
!!  blockvectorx  = on output, orthonormalized wavefunction.
!!
!!
!! PARENTS
!!      lobpcgccIIIwf,lobpcgccIIwf,m_lobpcg,pw_orthon
!!
!! CHILDREN
!!      wrtout,xmpi_sum,zgemm,zpotrf,ztrsm
!!
!! SOURCE

subroutine zorthonormalize(blockvectorx,blockvectorbx,blocksize,spaceComm,sqgram,vectsize)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,spaceComm,vectsize
!arrays
 complex(dpc),intent(in) :: blockvectorbx(vectsize,blocksize)
 complex(dpc),intent(inout) :: blockvectorx(vectsize,blocksize)
 complex(dpc),intent(out) :: sqgram(blocksize,blocksize)

!Local variables-------------------------------
!scalars
 integer :: ierr,info
 character(len=500) :: message

! *********************************************************************

 call abi_xgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
& vectsize,blockvectorbx,vectsize,czero,sqgram,blocksize)

 call xmpi_sum(sqgram,spaceComm,ierr)

 call abi_xpotrf('u',blocksize,sqgram,blocksize,info)

 if (info /= 0 )  then
   write(message,'(a,i0)')'zpotrf, info=',info
   ABI_ERROR(message)
 end if

 call abi_xtrsm('r','u','n','n',vectsize,blocksize,cone,sqgram,blocksize,blockvectorx,vectsize)

end subroutine zorthonormalize
!!***

