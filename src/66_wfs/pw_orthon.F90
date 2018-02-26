!{\src2tex{textfont=tt}}
!!****f* ABINIT/pw_orthon
!! NAME
!! pw_orthon
!!
!! FUNCTION
!! Normalize nvec complex vectors each of length nelem and then orthogonalize by modified Gram-Schmidt.
!! Two orthogonality conditions are available:
!!  Simple orthogonality: ${<Vec_{i}|Vec_{j}>=Delta_ij}$
!!  Orthogonality with overlap S: ${<Vec_{i}|S|Vec_{j}>=Delta_ij}$
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA,XG,GMR,FF,MT)
!! This file is distributed under the terms of the
!! GNU General Public License,see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt.
!!
!! INPUTS
!!  icg=shift to be given to the location of the data in cg(=vecnm)
!!  igsc=shift to be given to the location of the data in gsc(=ovl_vecnm)
!!  istwf_k=option parameter that describes the storage of wfs
!!  mcg=maximum size of second dimension of cg(=vecnm)
!!  mgsc=maximum size of second dimension of gsc(=ovl_vecnm)
!!  nelem=number of complex elements in each vector
!!  nvec=number of vectors to be orthonormalized
!!  ortalgo= option for the choice of the algorithm
!!         -1: no orthogonalization (direct return)
!!          0 or 2: old algorithm (use of buffers)
!!          1: new algorithm (use of blas)
!!          3: new new algorithm (use of lapack without copy)
!!  useoverlap=select the orthogonality condition
!!               0: no overlap between vectors
!!               1: vectors are overlapping
!!  me_g0=1 if this processor has G=0, 0 otherwise
!!  comm=MPI communicator
!!
!! SIDE EFFECTS
!!  vecnm= input: vectors to be orthonormalized; array of nvec column
!!                vectors,each of length nelem,shifted by icg
!!                This array is complex or else real(dp) of twice length
!!         output: orthonormalized set of vectors
!!  if (useoverlap==1) only:
!!    ovl_vecnm= input: product of overlap and input vectors:
!!                      S|vecnm>,where S is the overlap operator
!!               output: updated S|vecnm> according to vecnm
!!
!! NOTES
!! Note that each vector has an arbitrary phase which is not fixed in this routine.
!!
!! WARNING : not yet suited for nspinor=2 with istwfk/=1
!!
!! PARENTS
!!      lapackprof,vtowfk,wfconv
!!
!! CHILDREN
!!      abi_xcopy,abi_xorthonormalize,abi_xtrsm,cgnc_cholesky,cgnc_gramschmidt
!!      cgpaw_cholesky,cgpaw_gramschmidt,ortho_reim,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pw_orthon(icg,igsc,istwf_k,mcg,mgsc,nelem,nvec,ortalgo,ovl_vecnm,useoverlap,vecnm,me_g0,comm)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_cgtools
 use m_abi_linalg

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pw_orthon'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,igsc,istwf_k,mcg,mgsc,nelem,nvec,ortalgo,useoverlap,me_g0,comm
!arrays
 real(dp),intent(inout) :: ovl_vecnm(2,mgsc*useoverlap),vecnm(2,mcg)

!Local variables-------------------------------
!scalars
 integer :: ierr,ii,ii0,ii1,ii2,ivec,ivec2
 integer :: rvectsiz,vectsize,cg_idx,gsc_idx
 real(dp) :: doti,dotr,sum,xnorm
 character(len=500) :: msg
!arrays
 integer :: cgindex(nvec),gscindex(nvec)
 real(dp) :: buffer2(2),tsec(2)
 real(dp),allocatable :: rblockvectorbx(:,:),rblockvectorx(:,:),rgramxbx(:,:)
 complex(dpc),allocatable :: cblockvectorbx(:,:),cblockvectorx(:,:)
 complex(dpc),allocatable :: cgramxbx(:,:)

! *************************************************************************

#ifdef DEBUG_MODE
!Make sure imaginary part at G=0 vanishes
 if (istwf_k==2) then
   do ivec=1,nvec
     if(abs(vecnm(2,1+nelem*(ivec-1)+icg))>zero)then
!      if(abs(vecnm(2,1+nelem*(ivec-1)+icg))>tol16)then
       write(msg,'(2a,3i0,2es16.6,a,a)')&
&       ' For istwf_k=2,observed the following element of vecnm :',ch10,&
&       nelem,ivec,icg,vecnm(1:2,1+nelem*(ivec-1)+icg),ch10,&
&       '  with a non-negligible imaginary part.'
       MSG_BUG(msg)
     end if
   end do
 end if
#endif

!Nothing to do if ortalgo=-1
 if(ortalgo==-1) return

 do ivec=1,nvec
   cgindex(ivec)=nelem*(ivec-1)+icg+1
   gscindex(ivec)=nelem*(ivec-1)+igsc+1
 end do

 if (ortalgo==3) then
!  =========================
!  First (new new) algorithm
!  =========================
!  NEW VERSION: avoid copies, use ZHERK for NC
   cg_idx = cgindex(1)
   if (useoverlap==1) then
     gsc_idx = gscindex(1)
     call cgpaw_cholesky(nelem,nvec,vecnm(1,cg_idx),ovl_vecnm(1,gsc_idx),istwf_k,me_g0,comm)
   else
     call cgnc_cholesky(nelem,nvec,vecnm(1,cg_idx),istwf_k,me_g0,comm,use_gemm=.FALSE.)
   end if

 else if(ortalgo==1) then
!  =======================
!  Second (new) algorithm
!  =======================
!  This first algorithm seems to be more efficient especially in the parallel band-FFT mode.

   if(istwf_k==1) then
     vectsize=nelem
     ABI_ALLOCATE(cgramxbx,(nvec,nvec))
     ABI_ALLOCATE(cblockvectorx,(vectsize,nvec))
     ABI_ALLOCATE(cblockvectorbx,(vectsize,nvec))
     call abi_xcopy(nvec*vectsize,vecnm(:,cgindex(1):cgindex(nvec)-1),1,cblockvectorx,1,x_cplx=2)
     if (useoverlap == 1) then
       call abi_xcopy(nvec*vectsize,ovl_vecnm(:,gscindex(1):gscindex(nvec)-1),1,cblockvectorbx,1,x_cplx=2)
     else
       call abi_xcopy(nvec*vectsize,vecnm(:,cgindex(1):cgindex(nvec)-1),1,cblockvectorbx,1,x_cplx=2)
     end if
     call abi_xorthonormalize(cblockvectorx,cblockvectorbx,nvec,comm,cgramxbx,vectsize)
     call abi_xcopy(nvec*vectsize,cblockvectorx,1,vecnm(:,cgindex(1):cgindex(nvec)-1),1,x_cplx=2)
     if (useoverlap == 1) then
       call abi_xtrsm('r','u','n','n',vectsize,nvec,cone,cgramxbx,nvec,cblockvectorbx,vectsize)
       call abi_xcopy(nvec*vectsize,cblockvectorbx,1,ovl_vecnm(:,gscindex(1):gscindex(nvec)-1),1,x_cplx=2)
     end if
     ABI_DEALLOCATE(cgramxbx)
     ABI_DEALLOCATE(cblockvectorx)
     ABI_DEALLOCATE(cblockvectorbx)

   else if(istwf_k==2) then
     ! Pack real and imaginary part of the wavefunctions.
     rvectsiz=nelem
     vectsize=2*nelem; if(me_g0==1) vectsize=vectsize-1
     ABI_ALLOCATE(rgramxbx,(nvec,nvec))
     ABI_ALLOCATE(rblockvectorx,(vectsize,nvec))
     ABI_ALLOCATE(rblockvectorbx,(vectsize,nvec))
     do ivec=1,nvec
       if (me_g0 == 1) then
         call abi_xcopy(1,vecnm(1,cgindex(ivec)),1,rblockvectorx (1,ivec),1)
         call abi_xcopy(rvectsiz-1,vecnm(1,cgindex(ivec)+1),2,rblockvectorx(2,ivec),1)
         call abi_xcopy(rvectsiz-1,vecnm(2,cgindex(ivec)+1),2,rblockvectorx(rvectsiz+1,ivec),1)
         if (useoverlap == 1) then
           call abi_xcopy(1,ovl_vecnm(1,gscindex(ivec)),1,rblockvectorbx(1,ivec),1)
           call abi_xcopy(rvectsiz-1,ovl_vecnm(1,gscindex(ivec)+1),2,rblockvectorbx(2,ivec),1)
           call abi_xcopy(rvectsiz-1,ovl_vecnm(2,gscindex(ivec)+1),2,rblockvectorbx(rvectsiz+1,ivec),1)
         else
           call abi_xcopy(1,vecnm(1,cgindex(ivec)),1,rblockvectorbx(1,ivec),1)
           call abi_xcopy(rvectsiz-1,vecnm(1,cgindex(ivec)+1),2,rblockvectorbx(2,ivec),1)
           call abi_xcopy(rvectsiz-1,vecnm(2,cgindex(ivec)+1),2,rblockvectorbx(rvectsiz+1,ivec),1)
         end if
         rblockvectorx (2:vectsize,ivec)=rblockvectorx (2:vectsize,ivec)*sqrt2
         rblockvectorbx(2:vectsize,ivec)=rblockvectorbx(2:vectsize,ivec)*sqrt2
       else
         call abi_xcopy(rvectsiz,vecnm(1,cgindex(ivec)),2,rblockvectorx(1,ivec),1)
         call abi_xcopy(rvectsiz,vecnm(2,cgindex(ivec)),2,rblockvectorx(rvectsiz+1,ivec),1)
         if (useoverlap == 1) then
           call abi_xcopy(rvectsiz,ovl_vecnm(1,gscindex(ivec)),2,rblockvectorbx(1,ivec),1)
           call abi_xcopy(rvectsiz,ovl_vecnm(2,gscindex(ivec)),2,rblockvectorbx(rvectsiz+1,ivec),1)
         else
           call abi_xcopy(rvectsiz,vecnm(1,cgindex(ivec)),2,rblockvectorbx(1,ivec),1)
           call abi_xcopy(rvectsiz,vecnm(2,cgindex(ivec)),2,rblockvectorbx(rvectsiz+1,ivec),1)
         end if
         rblockvectorx (1:vectsize,ivec)=rblockvectorx (1:vectsize,ivec)*sqrt2
         rblockvectorbx(1:vectsize,ivec)=rblockvectorbx(1:vectsize,ivec)*sqrt2
       end if
     end do

     call ortho_reim(rblockvectorx,rblockvectorbx,nvec,comm,rgramxbx,vectsize)

     do ivec=1,nvec
       ! Unpack results
       if (me_g0 == 1) then
         call abi_xcopy(1,rblockvectorx(1,ivec),1,vecnm(1,cgindex(ivec)),1)
         vecnm(2,cgindex(ivec))=zero
         rblockvectorx(2:vectsize,ivec)=rblockvectorx(2:vectsize,ivec)/sqrt2
         call abi_xcopy(rvectsiz-1,rblockvectorx(2,ivec),1,vecnm(1,cgindex(ivec)+1),2)
         call abi_xcopy(rvectsiz-1,rblockvectorx(rvectsiz+1,ivec),1,vecnm(2,cgindex(ivec)+1),2)
       else
         rblockvectorx(1:vectsize,ivec)=rblockvectorx(1:vectsize,ivec)/sqrt2
         call abi_xcopy(rvectsiz,rblockvectorx(1,ivec),1,vecnm(1,cgindex(ivec)),2)
         call abi_xcopy(rvectsiz,rblockvectorx(rvectsiz+1,ivec),1,vecnm(2,cgindex(ivec)),2)
       end if

       if(useoverlap == 1) then
         call abi_xtrsm('r','u','n','n',vectsize,nvec,one,rgramxbx,nvec,rblockvectorbx,vectsize)
         if (me_g0 == 1) then
           call abi_xcopy(1,rblockvectorbx(1,ivec),1,ovl_vecnm(1,gscindex(ivec)),1)
           ovl_vecnm(2,gscindex(ivec))=zero
           rblockvectorbx(2:vectsize,ivec)=rblockvectorbx(2:vectsize,ivec)/sqrt2
           call abi_xcopy(rvectsiz-1,rblockvectorbx(2,ivec),1,ovl_vecnm(1,gscindex(ivec)+1),2)
           call abi_xcopy(rvectsiz-1,rblockvectorbx(rvectsiz+1,ivec),1,ovl_vecnm(2,gscindex(ivec)+1),2)
         else
           rblockvectorbx(1:vectsize,ivec)=rblockvectorbx(1:vectsize,ivec)/sqrt2
           call abi_xcopy(rvectsiz,rblockvectorbx(1,ivec),1,ovl_vecnm(1,gscindex(ivec)),2)
           call abi_xcopy(rvectsiz,rblockvectorbx(rvectsiz+1,ivec),1,ovl_vecnm(2,gscindex(ivec)),2)
         end if
       end if
     end do
     ABI_DEALLOCATE(rgramxbx)
     ABI_DEALLOCATE(rblockvectorx)
     ABI_DEALLOCATE(rblockvectorbx)
   end if

 else if (ortalgo==4) then 
!  else if (ANY(ortalgo==(/0,2/))) then 

   cg_idx = cgindex(1)
   if (useoverlap==0) then 
     call cgnc_gramschmidt(nelem,nvec,vecnm(1,cg_idx),istwf_k,me_g0,comm)
   else 
     gsc_idx = gscindex(1)
     call cgpaw_gramschmidt(nelem,nvec,vecnm(1,cg_idx),ovl_vecnm(1,gsc_idx),istwf_k,me_g0,comm)
   end if

 else if (ANY(ortalgo==(/0,2/))) then 
!  =======================
!  Third (old) algorithm
!  =======================

   do ivec=1,nvec
!    Normalize each vecnm(n,m) in turn:

     if (useoverlap==1) then ! Using overlap S...
       if(istwf_k/=2)then
         sum=zero;ii0=1
       else
         if (me_g0 ==1) then
           sum=half*ovl_vecnm(1,1+nelem*(ivec-1)+igsc)*vecnm(1,1+nelem*(ivec-1)+icg)
           ii0=2
         else
           sum=zero;ii0=1
         end if
       end if
!$OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:sum) SHARED(icg,ivec,nelem,vecnm)
       do ii=ii0+nelem*(ivec-1),nelem*ivec
         sum=sum+vecnm(1,ii+icg)*ovl_vecnm(1,ii+igsc)+vecnm(2,ii+icg)*ovl_vecnm(2,ii+igsc)
       end do

     else ! Without overlap...
       if(istwf_k/=2)then
         sum=zero;ii0=1
       else
         if (me_g0 ==1) then
           sum=half*vecnm(1,1+nelem*(ivec-1)+icg)**2
           ii0=2
         else
           sum=zero;ii0=1
         end if
       end if
!$OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:sum) SHARED(icg,ivec,nelem,vecnm)
       do ii=ii0+nelem*(ivec-1)+icg,nelem*ivec+icg
         sum=sum+vecnm(1,ii)**2+vecnm(2,ii)**2
       end do
     end if

     call timab(48,1,tsec)
     call xmpi_sum(sum,comm,ierr)
     call timab(48,2,tsec)

     if(istwf_k>=2)sum=two*sum
     xnorm = sqrt(abs(sum)) ;  sum=1.0_dp/xnorm
!$OMP PARALLEL DO PRIVATE(ii) SHARED(icg,ivec,nelem,sum,vecnm)
     do ii=1+nelem*(ivec-1)+icg,nelem*ivec+icg
       vecnm(1,ii)=vecnm(1,ii)*sum
       vecnm(2,ii)=vecnm(2,ii)*sum
     end do
     if (useoverlap==1) then
!$OMP PARALLEL DO PRIVATE(ii) SHARED(icg,ivec,nelem,sum,ovl_vecnm)
       do ii=1+nelem*(ivec-1)+igsc,nelem*ivec+igsc
         ovl_vecnm(1,ii)=ovl_vecnm(1,ii)*sum
         ovl_vecnm(2,ii)=ovl_vecnm(2,ii)*sum
       end do
     end if

!    Remove projection in all higher states.
     if (ivec<nvec) then

       if(istwf_k==1)then
!        Cannot use time-reversal symmetry

         if (useoverlap==1) then ! Using overlap.
           do ivec2=ivec+1,nvec
!            First compute scalar product
             dotr=zero ; doti=zero
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+igsc
!$OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:doti,dotr) SHARED(ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               dotr=dotr+vecnm(1,ii1+ii)*ovl_vecnm(1,ii2+ii)+vecnm(2,ii1+ii)*ovl_vecnm(2,ii2+ii)
               doti=doti+vecnm(1,ii1+ii)*ovl_vecnm(2,ii2+ii)-vecnm(2,ii1+ii)*ovl_vecnm(1,ii2+ii)
             end do

             call timab(48,1,tsec)
             buffer2(1)=doti;buffer2(2)=dotr
             call xmpi_sum(buffer2,comm,ierr)
             call timab(48,2,tsec)
             doti=buffer2(1)
             dotr=buffer2(2)

!            Then subtract the appropriate amount of the lower state
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+icg
#ifdef FC_INTEL
!            DIR$ ivdep
#endif
!$OMP PARALLEL DO PRIVATE(ii) SHARED(doti,dotr,ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               vecnm(1,ii2+ii)=vecnm(1,ii2+ii)-dotr*vecnm(1,ii1+ii)+doti*vecnm(2,ii1+ii)
               vecnm(2,ii2+ii)=vecnm(2,ii2+ii)-doti*vecnm(1,ii1+ii)-dotr*vecnm(2,ii1+ii)
             end do

             ii1=nelem*(ivec-1)+igsc;ii2=nelem*(ivec2-1)+igsc
             do ii=1,nelem
               ovl_vecnm(1,ii2+ii)=ovl_vecnm(1,ii2+ii)&
&               -dotr*ovl_vecnm(1,ii1+ii)&
&               +doti*ovl_vecnm(2,ii1+ii)
               ovl_vecnm(2,ii2+ii)=ovl_vecnm(2,ii2+ii)&
               -doti*ovl_vecnm(1,ii1+ii)&
&               -dotr*ovl_vecnm(2,ii1+ii)
             end do
           end do
         else
!          ----- No overlap -----
           do ivec2=ivec+1,nvec
!            First compute scalar product
             dotr=zero ; doti=zero
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+icg
!$OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:doti,dotr) SHARED(ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               dotr=dotr+vecnm(1,ii1+ii)*vecnm(1,ii2+ii)+&
&               vecnm(2,ii1+ii)*vecnm(2,ii2+ii)
               doti=doti+vecnm(1,ii1+ii)*vecnm(2,ii2+ii)-&
&               vecnm(2,ii1+ii)*vecnm(1,ii2+ii)
             end do
!            Init mpi_comm
             buffer2(1)=doti
             buffer2(2)=dotr
             call timab(48,1,tsec)
             call xmpi_sum(buffer2,comm,ierr)
!            call xmpi_sum(doti,spaceComm,ierr)
!            call xmpi_sum(dotr,spaceComm,ierr)
             call timab(48,2,tsec)
             doti=buffer2(1)
             dotr=buffer2(2)

!            Then subtract the appropriate amount of the lower state
#ifdef FC_INTEL
!            DIR$ ivdep
#endif
!$OMP PARALLEL DO PRIVATE(ii) SHARED(doti,dotr,ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               vecnm(1,ii2+ii)=vecnm(1,ii2+ii)-dotr*vecnm(1,ii1+ii)+&
&               doti*vecnm(2,ii1+ii)
               vecnm(2,ii2+ii)=vecnm(2,ii2+ii)-doti*vecnm(1,ii1+ii)-&
&               dotr*vecnm(2,ii1+ii)
             end do
           end do

         end if  ! Test on useoverlap

       else if(istwf_k==2)then
!        At gamma point use of time-reversal symmetry saves cpu time.

         if (useoverlap==1) then
!          ----- Using overlap -----
           do ivec2=ivec+1,nvec
!            First compute scalar product
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+igsc
             if (me_g0 ==1) then
               dotr=half*vecnm(1,ii1+1)*ovl_vecnm(1,ii2+1)
!              Avoid double counting G=0 contribution
!              Imaginary part of vecnm at G=0 should be zero,so only take real part
!$OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:dotr) SHARED(ii1,ii2,nelem,vecnm)
               do ii=2,nelem
                 dotr=dotr+vecnm(1,ii1+ii)*ovl_vecnm(1,ii2+ii)+&
&                 vecnm(2,ii1+ii)*ovl_vecnm(2,ii2+ii)
               end do
             else
               dotr=0._dp
!$OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:dotr) SHARED(ii1,ii2,nelem,vecnm)
               do ii=1,nelem
                 dotr=dotr+vecnm(1,ii1+ii)*ovl_vecnm(1,ii2+ii)+&
&                 vecnm(2,ii1+ii)*ovl_vecnm(2,ii2+ii)
               end do
             end if

             dotr=two*dotr

             call timab(48,1,tsec)
             call xmpi_sum(dotr,comm,ierr)
             call timab(48,2,tsec)

!            Then subtract the appropriate amount of the lower state
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+icg
#ifdef FC_INTEL
!            DIR$ ivdep
#endif
!$OMP PARALLEL DO PRIVATE(ii) SHARED(dotr,ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               vecnm(1,ii2+ii)=vecnm(1,ii2+ii)-dotr*vecnm(1,ii1+ii)
               vecnm(2,ii2+ii)=vecnm(2,ii2+ii)-dotr*vecnm(2,ii1+ii)
             end do
             ii1=nelem*(ivec-1)+igsc;ii2=nelem*(ivec2-1)+igsc
             do ii=1,nelem
               ovl_vecnm(1,ii2+ii)=ovl_vecnm(1,ii2+ii)-dotr*ovl_vecnm(1,ii1+ii)
               ovl_vecnm(2,ii2+ii)=ovl_vecnm(2,ii2+ii)-dotr*ovl_vecnm(2,ii1+ii)
             end do
           end do
         else
!          ----- No overlap -----
           do ivec2=ivec+1,nvec
!            First compute scalar product
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+icg
             if (me_g0 ==1) then
!              Avoid double counting G=0 contribution
!              Imaginary part of vecnm at G=0 should be zero,so only take real part
               dotr=half*vecnm(1,ii1+1)*vecnm(1,ii2+1)
!$OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:dotr) SHARED(ii1,ii2,nelem,vecnm)
               do ii=2,nelem
                 dotr=dotr+vecnm(1,ii1+ii)*vecnm(1,ii2+ii)+vecnm(2,ii1+ii)*vecnm(2,ii2+ii)
               end do
             else
               dotr=0._dp
!$OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:dotr) SHARED(ii1,ii2,nelem,vecnm)
               do ii=1,nelem
                 dotr=dotr+vecnm(1,ii1+ii)*vecnm(1,ii2+ii)+vecnm(2,ii1+ii)*vecnm(2,ii2+ii)
               end do
             end if
             dotr=two*dotr

             call timab(48,1,tsec)
             call xmpi_sum(dotr,comm,ierr)
             call timab(48,2,tsec)

!            Then subtract the appropriate amount of the lower state
#ifdef FC_INTEL
!            DIR$ ivdep
#endif
!$OMP PARALLEL DO PRIVATE(ii) SHARED(dotr,ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               vecnm(1,ii2+ii)=vecnm(1,ii2+ii)-dotr*vecnm(1,ii1+ii)
               vecnm(2,ii2+ii)=vecnm(2,ii2+ii)-dotr*vecnm(2,ii1+ii)
             end do
           end do
         end if  ! Test on useoverlap

       else
!        At other special points,use of time-reversal symmetry saves cpu time.

         if (useoverlap==1) then
!          ----- Using overlap -----
           do ivec2=ivec+1,nvec
!            First compute scalar product
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+igsc
!            Avoid double counting G=0 contribution
!            Imaginary part of vecnm at G=0 should be zero,so only take real part
             dotr=zero
!$OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:dotr) SHARED(ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               dotr=dotr+vecnm(1,ii1+ii)*ovl_vecnm(1,ii2+ii)+vecnm(2,ii1+ii)*ovl_vecnm(2,ii2+ii)
             end do
             dotr=two*dotr

             call timab(48,1,tsec)
             call xmpi_sum(dotr,comm,ierr)
             call timab(48,2,tsec)
             
!            Then subtract the appropriate amount of the lower state
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+icg
#ifdef FC_INTEL
!            DIR$ ivdep
#endif
!$OMP PARALLEL DO PRIVATE(ii) SHARED(dotr,ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               vecnm(1,ii2+ii)=vecnm(1,ii2+ii)-dotr*vecnm(1,ii1+ii)
               vecnm(2,ii2+ii)=vecnm(2,ii2+ii)-dotr*vecnm(2,ii1+ii)
             end do
             ii1=nelem*(ivec-1)+igsc;ii2=nelem*(ivec2-1)+igsc
             do ii=1,nelem
               ovl_vecnm(1,ii2+ii)=ovl_vecnm(1,ii2+ii)-dotr*ovl_vecnm(1,ii1+ii)
               ovl_vecnm(2,ii2+ii)=ovl_vecnm(2,ii2+ii)-dotr*ovl_vecnm(2,ii1+ii)
             end do
           end do
         else
!          ----- No overlap -----
           do ivec2=ivec+1,nvec
!            First compute scalar product
             ii1=nelem*(ivec-1)+icg;ii2=nelem*(ivec2-1)+icg
!            Avoid double counting G=0 contribution
!            Imaginary part of vecnm at G=0 should be zero,so only take real part
             dotr=zero
!$OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:dotr) SHARED(ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               dotr=dotr+vecnm(1,ii1+ii)*vecnm(1,ii2+ii)+vecnm(2,ii1+ii)*vecnm(2,ii2+ii)
             end do
             dotr=two*dotr

             call timab(48,1,tsec)
             call xmpi_sum(dotr,comm,ierr)
             call timab(48,2,tsec)

!            Then subtract the appropriate amount of the lower state
!$OMP PARALLEL DO PRIVATE(ii) SHARED(dotr,ii1,ii2,nelem,vecnm)
             do ii=1,nelem
               vecnm(1,ii2+ii)=vecnm(1,ii2+ii)-dotr*vecnm(1,ii1+ii)
               vecnm(2,ii2+ii)=vecnm(2,ii2+ii)-dotr*vecnm(2,ii1+ii)
             end do
           end do
         end if

!        End use of time-reversal symmetry
       end if

     end if  ! Test on "ivec"

!    end loop over vectors (or bands) with index ivec :
   end do

 else 
   write(msg,'(a,i0)')"Wrong value for ortalgo: ",ortalgo
   MSG_ERROR(msg)
 end if ! End of the second algorithm

end subroutine pw_orthon
!!***
