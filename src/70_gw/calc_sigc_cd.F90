!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_sigc_cd
!! NAME
!! calc_sigc_cd
!!
!! FUNCTION
!! Calculate contributions to the self-energy operator with the contour deformation method.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (FB, GMR, VO, LR, RWG, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nomega=Total number of frequencies where $\Sigma_c$ matrix elements are evaluated.
!!  nomegae=Number of frequencies where $\epsilon^{-1}$ has been evaluated.
!!  nomegaei=Number of imaginary frequencies for $\epsilon^{-1}$ (non zero).
!!  nomegaer=Number of real frequencies for $\epsilon^{-1}$
!!  npwc=Number of G vectors for the correlation part.
!!  npwx=Number of G vectors in rhotwgp for each spinorial component.
!!  nspinor=Number of spinorial components.
!!  theta_mu_minus_e0i=1 if e0i is occupied, 0 otherwise. Fractional occupancy in case of metals. 
!!  omegame0i(nomega)=Contains $\omega-\epsilon_{k-q,b1,\sigma}$
!!  epsm1q(npwc,npwc,nomegae)=Symmetrized inverse dielectric matrix (exchange part is subtracted).
!!  omega(nomegae)=Set of frequencies for $\epsilon^{-1}$.
!!  rhotwgp(npwx*nspinor)=Matrix elements: $<k-q,b1,\sigma|e^{-i(q+G)r} |k,b2,\sigma>*vc_sqrt$
!!
!! OUTPUT
!! ket(npwc,nomega)=Contains \Sigma_c(\omega)|\phi> in reciprocal space. 
!!
!! SIDE EFFECTS
!! npoles_missing=Incremented with the number of poles whose contribution has not been taken into account due to
!!  limited frequency mesh used for W.
!!
!! PARENTS
!!      calc_sigc_me,m_screen
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine calc_sigc_cd(npwc,npwx,nspinor,nomega,nomegae,nomegaer,nomegaei,rhotwgp,&
& omega,epsm1q,omegame0i,theta_mu_minus_e0i,ket,plasmafreq,npoles_missing,calc_poles,method)

 use defs_basis
 use m_profiling_abi
 use m_splines
 use m_xomp

 use m_gwdefs, only : czero_gw, cone_gw, Kron15N, Kron15W, Gau7W, &
&                     Kron23N, Kron23W, Gau11W, Kron31N, Kron31W, Gau15W
 use m_blas,   only : xgemv, xgemm

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_sigc_cd'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,nomegae,nomegaei,nomegaer,npwc,npwx,nspinor
 integer,intent(inout) :: npoles_missing
 real(dp),intent(in) :: theta_mu_minus_e0i,plasmafreq
!arrays
 real(dp),intent(in) :: omegame0i(nomega)
 complex(dpc),intent(in) :: omega(nomegae)
 complex(gwpc),intent(in) :: epsm1q(npwc,npwc,nomegae) 
 complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
 complex(gwpc),intent(inout) :: ket(nspinor*npwc,nomega)
 logical, intent(in), optional :: calc_poles(nomega)
 integer, intent(in), optional :: method

!Local variables-------------------------------
!scalars
 integer, parameter :: FABIEN=1,TRAPEZOID=2,NSPLINE=3 
 integer :: ii,ig,io,ios,ispinor,spadc,spadx,my_err,ierr,GK_LEVEL,INTMETHOD
 integer :: i,j
 real(dp) :: rt_imag,rt_real,local_one,local_zero
 real(dp) :: intsign,temp1,temp2,temp3,temp4
 real(dp) :: alph,inv_alph,beta,alphsq,betasq,inv_beta
 real(dp) :: re_intG,re_intK,im_intG,im_intK,GKttab,tau,ttil
 real(dp) :: ref,imf,r,s,r2,s2
 complex(dpc) :: ct,domegaleft,domegaright
 complex(gwpc) :: fact
!arrays
 real(dp) :: omegame0i_tmp(nomega),tmp_x(2),tmp_y(2)
 real(dp) :: left(nomega),right(nomega)
 real(dp) :: tbeta(nomega),tinv_beta(nomega),tbetasq(nomega)
 real(dp) :: atermr(nomega),aterml(nomega),logup(nomega),logdown(nomega)
 real(dp) :: rtmp_r(nomegaer),rtmp_i(nomegaer)
 real(dp) :: ftab(nomegaei+2),ftab2(nomegaei+2),xtab(nomegaei+2),y(3,nomegaei+2)
 real(dp) :: work(nomegaei+2),work2(nomegaei+2),y2(3,nomegaei+2)
 complex(dpc) :: omega_imag(nomegaei+1)
 complex(gwpc) :: epsrho(npwc,nomegae),epsrho_imag(npwc,nomegaei+1)
 complex(gwpc) :: tfone(npwc,nomegaei+1),tftwo(npwc,nomegaei+1)
 complex(gwpc) :: weight(nomegaei+1,nomega)
 complex(gwpc) :: weight2(nomegaei,nomega)
 logical :: my_calc_poles(nomega)
 real(dp), allocatable :: KronN(:),KronW(:),GaussW(:),fint(:),fint2(:)
!*************************************************************************

 my_calc_poles=.TRUE.
 my_err=0

 ! Set integration method for imaginary axis
 INTMETHOD = FABIEN
 if (present(method)) then
   if (method==1) INTMETHOD = FABIEN
   if (method==2) INTMETHOD = TRAPEZOID
   if (method>2) then
     INTMETHOD = NSPLINE
     if (method==3) then
       GK_LEVEL = 15
       ABI_ALLOCATE(KronN,(GK_LEVEL))
       ABI_ALLOCATE(KronW,(GK_LEVEL))
       ABI_ALLOCATE(GaussW,(GK_LEVEL-8))
       ABI_ALLOCATE(fint,(GK_LEVEL))
       ABI_ALLOCATE(fint2,(GK_LEVEL))
       KronN(:) = Kron15N(:); KronW(:) = Kron15W(:); GaussW(:) = Gau7W(:) 
     else if (method==4) then
       GK_LEVEL = 23
       ABI_ALLOCATE(KronN,(GK_LEVEL))
       ABI_ALLOCATE(KronW,(GK_LEVEL))
       ABI_ALLOCATE(GaussW,(GK_LEVEL-12))
       ABI_ALLOCATE(fint,(GK_LEVEL))
       ABI_ALLOCATE(fint2,(GK_LEVEL))
       KronN(:) = Kron23N(:); KronW(:) = Kron23W(:); GaussW(:) = Gau11W(:)
     else if (method>4) then
       GK_LEVEL = 31
       ABI_ALLOCATE(KronN,(GK_LEVEL))
       ABI_ALLOCATE(KronW,(GK_LEVEL))
       ABI_ALLOCATE(GaussW,(GK_LEVEL-16))
       ABI_ALLOCATE(fint,(GK_LEVEL))
       ABI_ALLOCATE(fint2,(GK_LEVEL))
       KronN(:) = Kron31N(:); KronW(:) = Kron31W(:); GaussW(:) = Gau15W(:) 
     end if
   end if
 end if

 ! Avoid divergences in $\omega - \omega_s$.
 omegame0i_tmp(:)=omegame0i(:)
 do ios=1,nomega
   if (ABS(omegame0i_tmp(ios))<tol6) omegame0i_tmp(ios)=sign(tol6,omegame0i_tmp(ios))
 end do

 do ispinor=1,nspinor
   spadx=(ispinor-1)*npwx
   spadc=(ispinor-1)*npwc
   !
   ! Calculate $ \sum_{Gp} (\epsilon^{-1}_{G Gp}(\omega)-\delta_{G Gp}) \rhotwgp(Gp) $
!$omp parallel do
   do io=1,nomegae
     call XGEMV('N',npwc,npwc,cone_gw,epsm1q(:,:,io),npwc,rhotwgp(1+spadx:),1,czero_gw,epsrho(:,io),1)
   end do

   ! Integrand along the imaginary axis.
   epsrho_imag(:,1)=epsrho(:,1)
   epsrho_imag(:,2:nomegaei+1)=epsrho(:,nomegaer+1:nomegae)

   ! Frequency mesh for integral along the imaginary axis.    
   omega_imag(1)=omega(1)
   omega_imag(2:nomegaei+1)=omega(nomegaer+1:nomegae)

   ! Original implementation -- saved here for reference during development
   ! === Perform integration along the imaginary axis ===
   !do io=1,nomegaei+1
   ! 
   !  if (io==1) then
   !    domegaleft  = omega_imag(io)
   !    domegaright =(omega_imag(io+1)-omega_imag(io  ))*half
   !  else if (io==nomegaei+1) then
   !    domegaleft  =(omega_imag(io  )-omega_imag(io-1))*half
   !    domegaright =(omega_imag(io  )-omega_imag(io-1))*half
   !  else
   !    domegaleft  =(omega_imag(io  )-omega_imag(io-1))*half
   !    domegaright =(omega_imag(io+1)-omega_imag(io  ))*half
   !  end if
   !  do ios=1,nomega
   !    omg2 = -AIMAG(omega_imag(io)+domegaright)/REAL(omegame0i_tmp(ios))
   !    omg1 = -AIMAG(omega_imag(io)-domegaleft )/REAL(omegame0i_tmp(ios))
   !    fact = ATAN(omg2)-ATAN(omg1)
   !    ket(spadc+1:spadc+npwc,ios)=ket(spadc+1:spadc+npwc,ios)+epsrho_imag(:,io)*fact
   !  end do
   !end do !io
   
   !ket(spadc+1:spadc+npwc,:)=ket(spadc+1:spadc+npwc,:)/pi
   !
   ! ---------------- end of original implementation -----------------------

 select case(INTMETHOD)
 case(FABIEN)
   ! Hopefully more effective implementation MS 04.08.2011
   ! === Perform integration along imaginary axis using BLAS ===
   ! First calculate first and last weights
   weight(1,:) = ATAN(-half*AIMAG(omega_imag(2))/REAL(omegame0i_tmp(:)))
   domegaleft  = (three*omega_imag(nomegaei+1)-omega_imag(nomegaei))
   domegaright = (omega_imag(nomegaei+1)+omega_imag(nomegaei))
   right(:)    = -AIMAG(omega_imag(nomegaei+1)-omega_imag(nomegaei))*REAL(omegame0i_tmp(:))
   left(:)     = quarter*AIMAG(domegaleft)*AIMAG(domegaright) &
&                 +REAL(omegame0i_tmp(:))*REAL(omegame0i_tmp(:))
   weight(nomegaei+1,:) = ATAN(right(:)/left(:))
   ! Calculate the rest of the weights
   do io=2,nomegaei
     domegaleft  = (omega_imag(io  )+omega_imag(io-1))
     domegaright = (omega_imag(io+1)+omega_imag(io  ))
     right(:)    = -half*AIMAG(omega_imag(io+1)-omega_imag(io-1))*REAL(omegame0i_tmp(:))
     left(:)     = REAL(omegame0i_tmp(:))*REAL(omegame0i_tmp(:)) &
&     +quarter*AIMAG(domegaleft)*AIMAG(domegaright)
     weight(io,:) = ATAN(right(:)/left(:))
   end do

   ! Use BLAS call to perform matrix-matrix multiplication and accumulation
   fact = CMPLX(piinv,zero) 

   call xgemm('N','N',npwc,nomega,nomegaei+1,fact,epsrho_imag,npwc,&
&   weight,nomegaei+1,cone_gw,ket(spadc+1:spadc+npwc,:),npwc)

!-----------------------------------------------------------
! Trapezoidal rule
 case(TRAPEZOID)

!  Transform omega coordinates
   alph     = plasmafreq
   alphsq   = alph*alph
   inv_alph = one/alph

   xtab(1:nomegaei+1) = AIMAG(omega_imag(:))/(AIMAG(omega_imag(:)) + alph)
   xtab(nomegaei+2)   = one

! Efficient trapezoidal rule with BLAS calls
   tbeta(:)     = REAL(omegame0i_tmp(:))
   tbetasq(:)   = tbeta(:)*tbeta(:)
   tinv_beta(:) = one/tbeta(:)
  
   do io=1,nomegaei
     atermr(:)    = inv_alph*tinv_beta(:)*((alphsq+tbetasq(:))*xtab(io+1)-tbetasq(:)) 
     aterml(:)    = inv_alph*tinv_beta(:)*((alphsq+tbetasq(:))*xtab(io  )-tbetasq(:))
     right(:)     = ATAN((atermr(:)-aterml(:))/(one+atermr(:)*aterml(:)))
     logup(:)     = ABS(((alphsq+tbetasq(:))*xtab(io+1)-two*tbetasq(:)) &
&                   *xtab(io+1)+tbetasq(:))
     logdown(:)   = ABS(((alphsq+tbetasq(:))*xtab(io  )-two*tbetasq(:)) &
&                   *xtab(io  )+tbetasq(:))
     ! Trapezoid integration weights
     weight(io,:)  = CMPLX(-(half*alph*tbeta(:)*LOG(logup(:)/logdown(:)) + tbetasq(:) &
&                       *right(:))/(alphsq+tbetasq(:)),zero)
     weight2(io,:) = CMPLX(-right(:),zero)
     ! Linear interpolation coefficients for each section (sum over ig)
     tfone(:,io)   = (epsrho_imag(:,io+1)-epsrho_imag(:,io)) &
&                    /(xtab(io+1)-xtab(io))
     tftwo(:,io)   = epsrho_imag(:,io) - tfone(:,io)*xtab(io)
   end do
   ! Calculate weights for asymptotic behaviour
   atermr(:)   = alph*tinv_beta(:) 
   aterml(:)   = inv_alph*tinv_beta(:)*((alphsq+tbetasq(:))*xtab(nomegaei+1)-tbetasq(:))
   logup(:)    = alphsq*xtab(nomegaei+1)*xtab(nomegaei+1)
   logdown(:)  = ABS(((alphsq+tbetasq(:))*xtab(nomegaei+1)-two*tbetasq(:)) &
&                 *xtab(nomegaei+1)+tbetasq(:))
   right(:)     = ATAN((atermr(:)-aterml(:))/(one+atermr(:)*aterml(:)))
   weight (nomegaei+1,:) = CMPLX(-(half*(alphsq*tinv_beta(:)*LOG(logdown(:)/logup(:)) &
&   - tbeta(:)*LOG(xtab(nomegaei+1)*xtab(nomegaei+1))) - alph*right(:)),zero)
   tfone(:,nomegaei+1) = -(zero-epsrho_imag(:,nomegaei+1)*AIMAG(omega_imag(nomegaei+1))) &
                         /(one-xtab(nomegaei+1))

   ! Use BLAS call to perform matrix-matrix multiplication and accumulation
   fact = CMPLX(piinv,zero) 

   call xgemm('N','N',npwc,nomega,nomegaei+1,fact,tfone,npwc,&
&   weight ,nomegaei+1,cone_gw,ket(spadc+1:spadc+npwc,:),npwc)
   call xgemm('N','N',npwc,nomega,nomegaei  ,fact,tftwo,npwc,&
&   weight2,nomegaei  ,cone_gw,ket(spadc+1:spadc+npwc,:),npwc)


!-----------------------------------------------------------
! Natural spline followed by Gauss-Kronrod
 case(NSPLINE)

!  Transform omega coordinates
   alph     = plasmafreq
   alphsq   = alph*alph
   inv_alph = one/alph

   xtab(1:nomegaei+1) = AIMAG(omega_imag(:))/(AIMAG(omega_imag(:)) + alph)
   xtab(nomegaei+2)   = one

! Gauss-Kronrod integration of spline fit of f(t)/(1-t) in transformed space
! *** OPENMP SECTION *** Added by MS
!!$OMP PARALLEL PRIVATE(ig,ftab,ftab2,s,s2,r,r2,y,y2,work,work2,beta,betasq,inv_beta, &
!!$OMP  intsign,io,ii,i,j,re_intG,re_intK,im_intG,im_intK,temp1,temp2,temp3,temp4, &
!!$OMP  ttil,tau,ref,fint,imf,fint2,GKttab) 
!!$OMP DO 
   do ig=1,npwc
     ! Spline fit
     ftab (1:nomegaei+1) =  REAL(epsrho_imag(ig,1:nomegaei+1))/(one-xtab(1:nomegaei+1))
     ftab2(1:nomegaei+1) = AIMAG(epsrho_imag(ig,1:nomegaei+1))/(one-xtab(1:nomegaei+1))
     ftab (nomegaei+2)   = zero; ftab2(nomegaei+2) = zero
     ! Explicit calculation of spline coefficients
     s  = zero; s2 = zero
     do i = 1, nomegaei+2-1
       r  = ( ftab (i+1) - ftab (i) ) / ( xtab(i+1) - xtab(i) )
       r2 = ( ftab2(i+1) - ftab2(i) ) / ( xtab(i+1) - xtab(i) )
       y (2,i) = r  - s; y2(2,i) = r2 - s2
       s  = r; s2 = r2
     end do 
     s = zero; s2 = zero
     r = zero; r2 = zero
     y(2,1) = zero; y2(2,1) = zero
     y(2,nomegaei+2) = zero; y2(2,nomegaei+2) = zero
     do i = 2, nomegaei+2-1
       y (2,i) = y (2,i) + r  * y (2,i-1)
       y2(2,i) = y2(2,i) + r2 * y2(2,i-1)
       work (i) = two * ( xtab(i-1) - xtab(i+1) ) - r  * s
       work2(i) = two * ( xtab(i-1) - xtab(i+1) ) - r2 * s2
       s = xtab(i+1) - xtab(i)
       s2 = s
       r  = s  / work (i)
       r2 = s2 / work2(i)
     end do 
     do j = 2, nomegaei+2-1
       i = nomegaei+2+1-j
       y (2,i) = ( ( xtab(i+1) - xtab(i) ) * y (2,i+1) - y (2,i) ) / work (i)
       y2(2,i) = ( ( xtab(i+1) - xtab(i) ) * y2(2,i+1) - y2(2,i) ) / work2(i)
     end do 
     do i = 1, nomegaei+2-1
       s = xtab(i+1) - xtab(i); s2 = s;
       r = y(2,i+1) - y(2,i); r2 = y2(2,i+1) - y2(2,i);
       y(3,i) = r / s; y2(3,i) = r2 / s2;
       y(2,i) = three * y(2,i); y2(2,i) = three * y2(2,i);
       y (1,i) = ( ftab (i+1) - ftab (i) ) / s  - ( y (2,i) + r  ) * s
       y2(1,i) = ( ftab2(i+1) - ftab2(i) ) / s2 - ( y2(2,i) + r2 ) * s2
     end do
     ! End of spline interpolation
     do ios=1,nomega
       beta     = REAL(omegame0i_tmp(ios))
       betasq   = beta*beta
       inv_beta = one/beta
       intsign = sign(half*piinv,beta)
       beta = ABS(beta)
       io = 1; re_intG = zero; re_intK = zero; im_intG = zero; im_intK = zero
       do ii=1,GK_LEVEL
         do
           GKttab = two*alph*xtab(io+1)/(beta-(beta-alph)*xtab(io+1))-one
           if (GKttab > KronN(ii)) EXIT
           io = io + 1
         end do
         temp1     = half*(KronN(ii)+one)
         temp2     = temp1 - half
         temp3     = temp2*temp2
         temp4     = half/(temp3 + quarter)
         ttil      = beta*temp1/(alph-(alph-beta)*temp1)
         tau       = ttil - xtab(io)
         ref       = ftab (io) + tau*(y (1,io)+tau*(y (2,io)+tau*y (3,io)))
         fint (ii) = -ref*(one-ttil)*temp4
         imf       = ftab2(io) + tau*(y2(1,io)+tau*(y2(2,io)+tau*y2(3,io))) 
         fint2(ii) = -imf*(one-ttil)*temp4
         re_intK   = KronW(ii)*fint (ii)  
         im_intK   = KronW(ii)*fint2(ii)
         ket(spadc+ig,ios) = ket(spadc+ig,ios)+intsign*CMPLX(re_intK,im_intK)
         end do ! ii
     end do !ios
   end do !ig
!!$OMP END DO
!!$OMP END PARALLEL

 end select !intmethod

   local_one = one
   local_zero = zero

   ! ============================================
   ! ==== Add contribution coming from poles ====
   ! ============================================
   ! First see if the contribution has been checked before the routine is entered
   if (present(calc_poles)) then
     my_calc_poles = calc_poles
   else ! Otherwise check locally if there is a contribution
     do ios=1,nomega
       if (omegame0i_tmp(ios)>tol12) then
         if ((local_one-theta_mu_minus_e0i)<tol12) my_calc_poles(ios) = .FALSE.
       end if
       if (omegame0i_tmp(ios)<-tol12) then
         if (theta_mu_minus_e0i<tol12) my_calc_poles(ios) = .FALSE.
       end if
     end do !ios
   end if

   if (ANY(my_calc_poles(:))) then ! Make sure we only enter if necessary

! *** OPENMP SECTION *** Added by MS
!!OMP !write(std_out,'(a,i0)') ' Entering openmp loop. Number of threads: ',xomp_get_num_threads()
!$OMP PARALLEL SHARED(npwc,nomega,nomegaer,theta_mu_minus_e0i,spadc,local_one,local_zero, &
!$OMP                    omega,epsrho,omegame0i_tmp,ket,my_calc_poles) &
!$OMP PRIVATE(ig,ios,rtmp_r,rtmp_i,tmp_x,tmp_y,rt_real,rt_imag,ct,ierr) REDUCTION(+:my_err) 
!!OMP $ write(std_out,'(a,i0)') ' Entering openmp loop. Number of threads: ',xomp_get_num_threads()
!$OMP DO 
     do ig=1,npwc
       !
       ! * Prepare the spline interpolation by filling at once the arrays rtmp_r, rtmp_i
       call spline(DBLE(omega(1:nomegaer)),DBLE(epsrho(ig,1:nomegaer)),nomegaer,local_zero,local_zero,rtmp_r)
       call spline(DBLE(omega(1:nomegaer)),DBLE(AIMAG(epsrho(ig,1:nomegaer))),nomegaer,local_zero,local_zero,rtmp_i)

       !! call spline_complex( DBLE(omega(1:nomegaer)), epsrho(ig,1:nomegaer), nomegaer, zero, zero, rtmp )

       do ios=1,nomega

         if (.NOT.my_calc_poles(ios)) CYCLE

         ! * Interpolate real and imaginary part of epsrho at |omegame0i_tmp|.
         tmp_x(1) = ABS(omegame0i_tmp(ios))
         call splint(nomegaer,DBLE(omega(1:nomegaer)),DBLE(epsrho(ig,1:nomegaer)),rtmp_r,1,tmp_x,tmp_y,ierr=ierr)
         if (ig==1.and.ispinor==1) my_err = my_err + ierr
         rt_real = tmp_y(1)

         tmp_x(1) = ABS(omegame0i_tmp(ios))
         call splint(nomegaer,DBLE(omega(1:nomegaer)),DBLE(AIMAG(epsrho(ig,1:nomegaer))),rtmp_i,1,tmp_x,tmp_y)
         rt_imag = tmp_y(1)

         !!call splint_complex(nomegaer,DBLE(omega(1:nomegaer)),epsrho(ig,1:nomegaer),rtmp,1,tmp_x,yfit)

         ct=DCMPLX(rt_real,rt_imag)

         if (omegame0i_tmp(ios)>tol12) then
           ket(spadc+ig,ios)=ket(spadc+ig,ios)+ct*(local_one-theta_mu_minus_e0i)
         end if
         if (omegame0i_tmp(ios)<-tol12) then
           ket(spadc+ig,ios)=ket(spadc+ig,ios)-ct*theta_mu_minus_e0i
         end if

       end do !ios
     end do !ig
!$OMP END DO
!$OMP END PARALLEL
   end if ! ANY(my_calc_poles)
 end do !ispinor

 npoles_missing = npoles_missing + my_err

 if (INTMETHOD>2) then
   ABI_DEALLOCATE(KronN)
   ABI_DEALLOCATE(KronW)
   ABI_DEALLOCATE(GaussW)
   ABI_DEALLOCATE(fint)
   ABI_DEALLOCATE(fint2)
 end if

end subroutine calc_sigc_cd
!!***
