!{\src2tex{textfont=tt}}
!!****f* ABINIT/pspnl_hgh_rec
!! NAME
!! pspnl_hgh_rec
!!
!! FUNCTION
!! This routine computes the matrices S_kk'^{l,A} of the projectors
!! (it is the exponential of the overlap matrix).
!! it coorresponds to the matrix:
!!   $$\left[(U_l)^{-1}*Exp(-temperature D_l )*U_l* (g_l)^{-1} -Identity\right]_kk'
!!   where (U_l)^-1* D_l* U_l = h^lg_l.
!! $g_l = <f^l_k|f^l_{k'}>$ is the overlap matrix between projectors
!! and $h^l_{kk'}$ is the strength matrix of the projectors.
!! It calulates also the strength eigenvalues and eigenvectors of $h$,
!! used in the calculus of non-local energy
!!
!! COPYRIGHT
!! Copyright (C) 2008-2016 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  temperature=4*rtrotter/beta=4*rtrotter*tsmear: the effective temp. in  recursion
!!  psps <type(pseudopotential_type)>=variables related to pseudo-potentials 
!!  debug=debug variable
!! 
!! OUTPUT
!!  
!!  nlrec%mat_exp_psp_nl=the matrix of the exponential of the projectors:
!!   for any psp, for any angular moment:
!!   h_ij=strength; g_ij=ovelap => exp(-h.g/temp/4p).g^(-1)
!!  nlrec%radii=Local radii of nl psp
!!  nlrec%pspinfo(:,:) for any typat:  (momang,typat)=number of projectors
!!  nlrec%eival(:,:,:) for any psp, any momang, the eigenvalues of the
!!    strength matrix H: eival(:,mang,psp)
!!  nlrec%eivec(:,:,:,:)for any psp, any momang, the eigenvectors of the
!!    strength matrix H: eivec(:,:,mang,psp)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_rec
!!
!! CHILDREN
!!      dgetrf,dgetri,dsyev,exp_mat,gamma_function,set2unit,wrtout
!!
!! NOTES
!!  at this time :
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pspnl_hgh_rec(psps,temperature,nlrec,debug)

 use m_profiling_abi
  
 use defs_basis
 use defs_datatypes
 use defs_rectypes
 use m_exp_mat,       only : exp_mat
 use m_numeric_tools, only : set2unit
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pspnl_hgh_rec'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none
   
!Arguments -----------------------------------
!scalars
 real(dp),intent(in) :: temperature
 logical,intent(in) :: debug
 type(pseudopotential_type),intent(in) :: psps
 type(nlpsprec_type),intent(inout) :: nlrec
!arrays
!Local variables-------------------------------
!scalars
 integer,parameter :: maxsize=3
 integer,parameter :: lwork=(1+32)*maxsize
 integer :: iangol,ipseudo,info,nproj 
 integer :: g_mat_size,ii,nproj2
 real(dp) :: denom_1,denom_2,tot_proj
 character(len=500) :: msg
!arrays
 integer :: ipvt(1)
 !real(dp) :: rwork(2*maxsize)
 real(dp) :: h_mat_init(3,3), rework(lwork)
 real(dp), allocatable :: g_mat(:,:),h_mat(:,:),eig_val_h(:)
 real(dp), allocatable :: identity(:,:),inv_g_mat(:,:),u_mat(:,:)

 complex(dpc),allocatable :: hg_mat(:,:)
! *************************************************************************

 if(debug)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' pspnl_hgh_rec : enter '
   call wrtout(std_out,msg,'COLL')
 end if

!--For any pseudo potential:
 do ipseudo = 1, psps%npsp !--Loop on the pseudos
   write(msg,'(a,80a)')' pseudo file',('-',ii=1,10)
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a)') psps%filpsp(ipseudo)
   call wrtout(std_out,msg,'COLL')
   
!  --For any angular moment:
   do iangol = 0,psps%mpsang-1 !--Loop on the angular moment
     
!    --Local radius
     nlrec%radii(iangol+1,ipseudo) = psps%gth_params%psppar(iangol+1,0,ipseudo)

!    --Strenghts of non-local projectors  (matrix h)
!    --Diagonal part:
     h_mat_init = zero
     h_mat_init(1,1) = psps%gth_params%psppar(iangol+1,1,ipseudo) 
     h_mat_init(2,2) = psps%gth_params%psppar(iangol+1,2,ipseudo) 
     h_mat_init(3,3) = psps%gth_params%psppar(iangol+1,3,ipseudo) 
!    --Out-diagonal part 
!    --Depending on angular moment the projectors 
!    strength is calculated differently
     select case(iangol)
     case(0) 
       h_mat_init(1,2) = -half*sqrt(3.d0/5.d0)*h_mat_init(2,2)
       h_mat_init(1,3) =  half*sqrt(5.d0/21.d0)*h_mat_init(3,3)
       h_mat_init(2,3) = -half*sqrt(100.d0/63.d0)*h_mat_init(3,3)
     case(1) 
       h_mat_init(1,2) = -half*sqrt(5.d0/7.d0)*h_mat_init(2,2)
       h_mat_init(1,3) =  sixth*sqrt(35.d0/11.d0)*h_mat_init(3,3)
       h_mat_init(2,3) = -14.d0/six/sqrt(11.d0) *h_mat_init(3,3)
     case(2)
       h_mat_init(1,2) = -half*sqrt(7.d0/9.d0)*h_mat_init(2,2)
       h_mat_init(1,3) =  half*sqrt(63.d0/143.d0)*h_mat_init(3,3)
       h_mat_init(2,3) = -nine/sqrt(143.d0)*h_mat_init(3,3)
     case(3)
       h_mat_init(1,2) = zero;  h_mat_init(1,3) = zero;  h_mat_init(2,3) = zero; 
     case default 
       write(msg,'(a)')' error angular: moment component'
       call wrtout(std_out,msg,'COLL')
     end select
     
     

!    --Real dimensions of projectors.
     g_mat_size = count(abs((/ (h_mat_init(ii,ii),ii=1,3) /))>1.d-8)
     nlrec%pspinfo(iangol+1,ipseudo) = g_mat_size   
     write(msg,'(a,i2,a,i2)')' ang. moment=',iangol,', N projectors=',g_mat_size
     call wrtout(std_out,msg,'COLL')
     if (g_mat_size>0) then
!      --Identity matrix
       ABI_ALLOCATE(identity,(g_mat_size,g_mat_size))
       call set2unit(identity)
!      identity = zero
!      identity(:,1) = one
!      identity(:,:) = cshift(array=identity,shift=(/ (-ii,ii=0,g_mat_size) /), dim=2 )

       
!      ############## CALCULOUS OF THE EIGEN_SPACE OF THE PROJECTORS STRENGTHS ##################
!      --Inverse of the matrix h
       ABI_ALLOCATE(eig_val_h,(g_mat_size))
       ABI_ALLOCATE(u_mat,(g_mat_size,g_mat_size))
!      --u-mat will contain the eigenvectors of h_mat_init
       u_mat = h_mat_init(:g_mat_size,:g_mat_size)

!      write(std_out,*)'hmat_init'
!      do ii=1,g_mat_size
!      write(std_out,*)h_mat_init(ii,:)
!      end do
       call DSYEV('v','u',g_mat_size,u_mat,g_mat_size,eig_val_h,rework,lwork,info)
       
!      --THE DIAGONAL MATRIX IS GIVEN BY D=U^t.H.U 
!      (eival=transpose(eivec).h_mat_init.eivec)
       write(msg,'(a,3d10.3)')'  eigenvalues=',eig_val_h
       call wrtout(std_out,msg,'COLL')
!      write(std_out,*)'autovec';write(std_out,*)u_mat
       
       nlrec%eival(:g_mat_size,1+iangol,ipseudo) = eig_val_h
       nlrec%eivec(:g_mat_size,:g_mat_size,1+iangol,ipseudo) = u_mat
       ABI_DEALLOCATE(eig_val_h)
       ABI_DEALLOCATE(u_mat)

!      ##########END  CALCULOUS OF THE EIGEN_SPACE OF THE PROJECTORS STRENGTHS ##################

       ABI_ALLOCATE(g_mat,(g_mat_size,g_mat_size))
       ABI_ALLOCATE(inv_g_mat,(g_mat_size,g_mat_size))
       ABI_ALLOCATE(h_mat,(g_mat_size,g_mat_size))
       ABI_ALLOCATE(hg_mat,(g_mat_size,g_mat_size))
       
       g_mat(:,:) = one
       h_mat(:,:) = zero
       h_mat(:,:) = h_mat_init(:g_mat_size,:g_mat_size)

!      -------------------------------------------------------
!      --Matrix  of the overlap between projetors (matrix g) 
!      and the h matrix of strength
       do nproj = 1,g_mat_size-1
         do nproj2 = 1+nproj,g_mat_size
           tot_proj = zero            
!          --Analytic value of overlap
!          g_ij=Gamma[-1/2+i+j+l]/Sqrt(Gamma[-1/2+i+iangol]*Gamma[-1/2+j+iangol])
           call gamma_function(-half+real(nproj+nproj2+iangol,dp),tot_proj)
           call gamma_function(-half+real(iangol+2*nproj,dp),denom_1)
           call gamma_function(-half+real(iangol+2*nproj2,dp),denom_2)

           g_mat(nproj,nproj2) = tot_proj/sqrt(denom_1*denom_2)
           g_mat(nproj2,nproj) = g_mat(nproj,nproj2)
           
           h_mat(nproj,nproj2) = h_mat_init(nproj,nproj2)
           h_mat(nproj2,nproj) = h_mat_init(nproj,nproj2)
         end do
       end do

!      --Inverse of the overlap matrix g
       inv_g_mat = g_mat   
       call DGETRF(g_mat_size,g_mat_size,inv_g_mat,g_mat_size,ipvt,info)
       call DGETRI(g_mat_size,inv_g_mat,g_mat_size,ipvt,rework,lwork,info)  

       
!      -----------------------------------------------------------
!      --Now it calculates the exponential of the matrix h.g
       hg_mat = matmul(h_mat,g_mat)
       
       call exp_mat(hg_mat,g_mat_size,-one/temperature)
       
!      --(exp(h.g)-Identity).(g^-1)
       hg_mat = hg_mat-identity(:,:)


!      --results on output
       nlrec%mat_exp_psp_nl(:g_mat_size,:g_mat_size,1+iangol,ipseudo) = matmul(real(hg_mat,dp),inv_g_mat)

!      write(std_out,*) nlrec%mat_exp_psp_nl(:g_mat_size,:g_mat_size,1+iangol,ipseudo)

       ABI_DEALLOCATE(g_mat)
       ABI_DEALLOCATE(hg_mat)
       ABI_DEALLOCATE(h_mat)
       ABI_DEALLOCATE(inv_g_mat)
       ABI_DEALLOCATE(identity)
     end if

   end do !enddo on angular moment
 end do !enddo on pseudos


 if(debug)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' pspnl_hgh_rec : exit '
   call wrtout(std_out,msg,'COLL')
 end if


end subroutine pspnl_hgh_rec
!!***
