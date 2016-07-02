!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_sigc_pole_cd
!! NAME
!! calc_sigc_pole_cd
!!
!! FUNCTION
!! Calculate contributions to the self-energy operator with the contour deformation method,
!! using a pole-fit screening.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (FB, GMR, VO, LR, RWG, MG)
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
!! npoles_missing=Incremeented with the number of poles whose contribution has not been taken into account due to
!!  limited frequency mesh used for W.
!!
!! PARENTS
!!
!! CHILDREN
!!      cgemm,re_and_im_screening_with_phase,spline,splint,xgemv,zgemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine calc_sigc_pole_cd(npwc,npwx,nspinor,ncoeff,nomega,nomegae,nomegaer,nomegaei,rhotwgp,&
& omega,epsm1q,omegame0i,theta_mu_minus_e0i,ket,npoles_missing,calc_poles)

 use m_profiling_abi

 use defs_basis
 use m_splines

 use m_gwdefs, only : czero_gw, cone_gw
 use m_blas,   only : xgemv
 use m_model_screening, only : re_and_im_screening_with_phase
!$ use omp_lib, only: omp_get_thread_num, omp_get_num_threads

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_sigc_pole_cd'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,ncoeff,nomegae,nomegaei,nomegaer,npwc,npwx,nspinor
 integer,intent(inout) :: npoles_missing
 real(dp),intent(in) :: theta_mu_minus_e0i
!arrays
 real(dp),intent(in) :: omegame0i(nomega)
 complex(dpc),intent(in) :: omega(nomegae)
 real(gwp) :: epsm1q(npwc,npwc,ncoeff) 
 complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
 complex(gwpc),intent(inout) :: ket(nspinor*npwc,nomega)
 logical, intent(in), optional :: calc_poles(nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,io,ios,ispinor,spadc,spadx,my_err,ig1,ig2 
 real(dp) :: rt_imag,rt_real,local_one,local_zero
! real(dp) :: intresult,intsign
 complex(dpc) :: ct,domegaleft,domegaright
#if defined HAVE_GW_DPC
 complex(dpc) :: fact
#else
 complex(spc) :: fact
#endif
!arrays
 real(dp) :: omegame0i_tmp(nomega),tmp_x(2),tmp_y(2)
 real(dp) :: rtmp_r(nomegaer),rtmp_i(nomegaer)
! real(dp) :: ftab(nomegaei+2),xtab(nomegaei+2),y(3,nomegaei+2)
! real(dp) :: work(nomegaei+2),e(nomegaei+2)
 complex(dpc) :: omega_imag(nomegaei+1)
 complex(gwpc) :: epsrho(npwc,nomegae),epsrho_imag(npwc,nomegaei+1)
 complex(gwpc) :: weight(nomegaei+1,nomega)
 complex(gwpc) :: work(npwc,npwc,nomegae)
 logical :: my_calc_poles(nomega)
!*************************************************************************

 my_err=0
 my_calc_poles=.TRUE.

 ! Avoid divergences in $\omega - \omega_s$.
 omegame0i_tmp(:)=omegame0i(:)
 do ios=1,nomega
   if (ABS(omegame0i_tmp(ios))<tol6) omegame0i_tmp(ios)=tol6
 end do

 do ispinor=1,nspinor
   spadx=(ispinor-1)*npwx
   spadc=(ispinor-1)*npwc
   !
   ! Calculate $ \sum_{Gp} (\epsilon^{-1}_{G Gp}(\omega)-\delta_{G Gp}) \rhotwgp(Gp) $
   ! First calculate the values from the poles
   do ig2=1,npwc
     do ig1=1,npwc
       call re_and_im_screening_with_phase(omega,work(ig1,ig2,:),nomegae, &
&       epsm1q(ig1,ig2,:),ncoeff)
     end do
   end do
    
   do io=1,nomegae
     call XGEMV('N',npwc,npwc,cone_gw,work(:,:,io),npwc,rhotwgp(1+spadx:),1,czero_gw,epsrho(:,io),1)
   end do

   !write(std_out,*) ch10,' Epsrho: ',epsrho

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

   ! Hopefully more effective implementation MS 04.08.2011
   ! === Perform integration along imaginary axis using BLAS ===
   ! First calculate first and last weights
   domegaleft  = omega_imag(1)
   domegaright = (omega_imag(2)-omega_imag(1))*half
   weight(1,:) = ATAN(-AIMAG(omega_imag(1)+domegaright)/REAL(omegame0i_tmp(:)))&
&               -ATAN(-AIMAG(omega_imag(1)-domegaleft)/REAL(omegame0i_tmp(:)))
   domegaleft  = (omega_imag(nomegaei+1)-omega_imag(nomegaei))*half
   domegaright = (omega_imag(nomegaei+1)-omega_imag(nomegaei))*half
   weight(nomegaei+1,:) = ATAN(-AIMAG(omega_imag(nomegaei+1)+domegaright)/REAL(omegame0i_tmp(:)))&
&               -ATAN(-AIMAG(omega_imag(nomegaei+1)-domegaleft)/REAL(omegame0i_tmp(:)))
   ! Calculate the rest of the weights
   do io=2,nomegaei
     domegaleft  = (omega_imag(io  )-omega_imag(io-1))*half
     domegaright = (omega_imag(io+1)-omega_imag(io  ))*half
     weight(io,:) = ATAN(-AIMAG(omega_imag(io)+domegaright)/REAL(omegame0i_tmp(:)))&
&                  -ATAN(-AIMAG(omega_imag(io)-domegaleft)/REAL(omegame0i_tmp(:))) 
   end do

   ! Use BLAS call to perform matrix-matrix multiplication and accumulation
   fact = CMPLX(piinv,zero) 
#if defined HAVE_GW_DPC
   call zgemm('N','N',npwc,nomega,nomegaei+1,fact,epsrho_imag,npwc,&
&   weight,nomegaei+1,cone,ket(spadc+1:spadc+npwc,:),npwc)
#else
   call cgemm('N','N',npwc,nomega,nomegaei+1,fact,epsrho_imag,npwc,&
&   weight,nomegaei+1,CMPLX(1.0_sp,0.0_sp),ket(spadc+1:spadc+npwc,:),npwc)
#endif
   
!#else
! Use spline integration
!do ios=1,nomega
!  if (ABS(omegame0i(ios))<tol12) then
!    intsign = sign(one,-REAL(omegame0i(ios)))
!    ket(spadc+1:spadc+npwc,ios) = ket(spadc+1:spadc+npwc,ios) + half*intsign*epsrho_imag(:,1)
!  else
!    xtab(1:nomegaei+1) = ABS(two*piinv*ATAN(-AIMAG(omega_imag(:))/REAL(omegame0i(ios))))
!    intsign = sign(one,-REAL(omegame0i(ios)))
!    xtab(nomegaei+2)   = one
!    do ig=1,npwc
!      ftab(1:nomegaei+1) = epsrho_imag(ig,:)
!      ftab(nomegaei+2)   = zero
!      call cspint(ftab,xtab,nomegaei+2,zero,one,y,e,work,intresult)
!      ket(spadc+ig,ios) = ket(spadc+ig,ios) + half*intsign*intresult
!    end do
!  end if
!end do
!#endif

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

!!OMP *** OPENMP SECTION *** Added by MS
!!OMP!$ write(std_out,'(a,i0)') ' Entering openmp loop. Number of threads: ',omp_get_num_threads()
!$OMP PARALLEL SHARED(npwc,nomega,nomegaer,theta_mu_minus_e0i,spadc,local_one,local_zero, &
!$OMP                    omega,epsrho,omegame0i_tmp,ket,my_calc_poles) &
!$OMP          PRIVATE(ig,ios,rtmp_r,rtmp_i,tmp_x,tmp_y,rt_real,rt_imag,ct)
!!OMP!$ write(std_out,'(a,i0)') ' Entering openmp loop. Number of threads: ',omp_get_num_threads()
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
         call splint(nomegaer,DBLE(omega(1:nomegaer)),DBLE(epsrho(ig,1:nomegaer)),rtmp_r,1,tmp_x,tmp_y,ierr=my_err)
         if (ig==1.and.ispinor==1) npoles_missing = npoles_missing + my_err
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

end subroutine calc_sigc_pole_cd
!!***
