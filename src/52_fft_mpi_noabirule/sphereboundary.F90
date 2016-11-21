!{\src2tex{textfont=tt}}
!!****f* ABINIT/sphereboundary
!! NAME
!! sphereboundary
!!
!! FUNCTION
!! Finds the boundary of the basis sphere of G vectors (at a given
!! k point) for use in improved zero padding of ffts in 3 dimensions.
!! Provides data to be used by subroutine fourwf, in the form of
!! an array gbound(2*mgfft+8,2).
!!
!! The first component (for use when mod(fftalg,10)==2))
!! provides integer values g1min,g1max,g2min,g2max
!! and then for g2 in the
!! sequence g2=0,1,2,...,g2max,g2min,g2min+1,...,-1, provides g1min, g1max.
!!
!! The second component (for use when mod(fftalg,10)==1))
!! provides integer values g1min,g1max,g3min,g3max,
!! where g3min and g3max have been corrected in case of time-reversal
!! and then for g3 in the sequence
!! g3=0,1,2,...,g3max,g3min,g3min+1,...,-1, provides g2min, g2max.
!! (also corrected in case of time-reversal)
!!
!! These values are stored in the above order in array gbound.
!! Debug mode, if fftalg is between 000 and 099
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  istwf_k=option parameter that describes the storage of wfs
!!  kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!!  mgfft=maximum size of 1D FFTs (only for dimensioning purposes)
!!  npw=number of G vectors in basis at this k point
!!
!! OUTPUT
!!  gbound(2*mgfft+8,2)=defined above
!!
!! PARENTS
!!      dens_in_sph,dfpt_eltfrkin,dfpt_mkrho,fock_getghc,m_bandfft_kpt,m_cut3d
!!      m_fft,m_fft_prof,m_fock,m_gsphere,m_hamiltonian,m_wfd,mkrho,mlwfovlp
!!      pawmkaewf,posdoppler,scfcv,spin_current,suscep_stat,susk,tddft,wfconv
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine sphereboundary(gbound,istwf_k,kg_k,mgfft,npw)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sphereboundary'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,mgfft,npw
!arrays
 integer,intent(in) :: kg_k(3,npw)
 integer,intent(out) :: gbound(2*mgfft+8,2)

!Local variables-------------------------------
!scalars
 integer :: dim_a,dim_b,fftalgc,g_a,gmax_a,gmax_b,gmax_b1,gmax_b2,gmin_a,gmin_b
 integer :: gmin_b1,gmin_b2,igb,ii,iloop,ipw,testm,testp,kgk
 character(len=500) :: message
!arrays
 integer :: gmax(2),gmin(2)

! *************************************************************************
!
!DEBUG
!write(std_out,*)' sphereboundary : enter'
!write(std_out, '(a)' )' sphereboundary : list of plane waves coordinates for k point '
!write(std_out, '(a)' )'       ipw       kg_k(1:3,ipw) '
!do ipw=1,npw
!write(std_out, '(i10,a,3i6)' )ipw,'  ',kg_k(1:3,ipw)
!end do
!gbound=-999
!ENDDEBUG

!Determine cube boundaries
 gbound(1,1)=minval(kg_k(1,:))
 gbound(2,1)=maxval(kg_k(1,:))
 gbound(1:2,2)=gbound(1:2,1)

!Treat differently the fftalgc cases
 do ii=1,2

   fftalgc=3-ii

   if(fftalgc/=2)then
     dim_a=3
     dim_b=2
   else
     dim_a=2
     dim_b=1
   end if

!  Relevant boundaries
   gbound(3,ii)=minval(kg_k(dim_a,:))
   gbound(4,ii)=maxval(kg_k(dim_a,:))
   gmin_a=gbound(3,ii)
   gmax_a=gbound(4,ii)

!  Must complete the sphere for fftalgc==1 and special storage modes.
!  Explanation : sg_fftpad is not able to take into account
!  the time-reversal symmetry, so that the boundaries will not be delimited
!  by the kg_k set, but by their symmetric also.
   if(istwf_k>=2 .and. fftalgc==1)then
     if( istwf_k==2 .or. istwf_k==3 .or. istwf_k==6 .or. istwf_k==7 )then
       gbound(4,2)=max(gmax_a,-gmin_a)
       gbound(3,2)=-gbound(4,2)
     else if( istwf_k==4 .or. istwf_k==5 .or. istwf_k==8 .or. istwf_k==9 )then
       gbound(4,2)=max(gmax_a,-gmin_a-1)
       gbound(3,2)=-gbound(4,2)-1
     end if
     gmax_a=gbound(4,2) ; gmin_a=gbound(3,2)
   end if

   igb=5

!  Consider first every g_a in range 0 ... gmax_a, then gmin_a ... -1
   gmin(1)=0         ; gmax(1)=gmax_a
   gmin(2)=gmin_a    ; gmax(2)=-1

   do iloop=1,2

     if( gmin(iloop) <= gmax(iloop) )then

       do g_a=gmin(iloop),gmax(iloop)

         if(istwf_k==1 .or. fftalgc/=1)then
           ! Select the minimal and maximal values, in the selected plane
           gmin_b=mgfft+1 ! Initialized with a value larger than all possible ones
           gmax_b=-mgfft-1 ! Initialized with a value smaller than all possible ones
           do ipw=1,npw
             if(kg_k(dim_a,ipw)==g_a)then
               kgk=kg_k(dim_b,ipw)
               if(kgk<=gmin_b)gmin_b=kgk
               if(kgk>=gmax_b)gmax_b=kgk
             end if
           end do

         else if(istwf_k>=2 .and. fftalgc==1)then

           ! Here, must take into account time-reversal symmetry explicitely

           ! Determine the boundaries for the plane g_a
           testp=0
           if(g_a<=gmax_a)then
             ! Select the minimal and maximal values, in the selected plane
             gmin_b1=mgfft+1 ! Initialized with a value larger than all possible ones
             gmax_b1=-mgfft-1 ! Initialized with a value smaller than all possible ones
             do ipw=1,npw
               if(kg_k(dim_a,ipw)==g_a)then
                 kgk=kg_k(dim_b,ipw)
                 if(kgk<=gmin_b1)gmin_b1=kgk
                 if(kgk>=gmax_b1)gmax_b1=kgk
               end if
             end do


             testp=1
           end if

           ! Determine the boundaries for the plane -g_a or -g_a-1
           testm=0
           if( istwf_k==2 .or. istwf_k==3 .or. istwf_k==6 .or. istwf_k==7 )then

             if(-g_a>=gmin_a)then
               ! Select the minimal and maximal values, in the selected plane
               ! Warning : there is an inversion of search (might be confusing)
               gmax_b2=mgfft+1 ! Initialized with a value larger than all possible ones
               gmin_b2=-mgfft-1 ! Initialized with a value smaller than all possible ones
               do ipw=1,npw
                 if(kg_k(dim_a,ipw)==-g_a)then
                   kgk=kg_k(dim_b,ipw)
                   if(kgk<=gmax_b2)gmax_b2=kgk
                   if(kgk>=gmin_b2)gmin_b2=kgk
                 end if
               end do
               testm=1
             end if

           else if( istwf_k==4 .or. istwf_k==5 .or. istwf_k==8 .or. istwf_k==9 )then

             if(-g_a-1>=gmin_a)then
               ! Select the minimal and maximal values, in the selected plane
               ! Warning : there is an inversion of search (might be confusing)
               gmax_b2=mgfft+1 ! Initialized with a value larger than all possible ones
               gmin_b2=-mgfft-1 ! Initialized with a value smaller than all possible ones
               do ipw=1,npw
                 if(kg_k(dim_a,ipw)==-g_a-1)then
                   kgk=kg_k(dim_b,ipw)
                   if(kgk<=gmax_b2)gmax_b2=kgk
                   if(kgk>=gmin_b2)gmin_b2=kgk
                 end if
               end do
               testm=1
             end if

           end if

           !  Must invert the boundaries, to use them for plane g_a
           if(testm==1)then
             ! This is needed to avoid border effect
             ! if the search did not lead to any element
             gmin_b2=max(gmin_b2,-mgfft) ; gmax_b2=min(gmax_b2,mgfft)
             if(istwf_k<=5)then
               gmax_b2=-gmax_b2 ; gmin_b2=-gmin_b2
             else
               gmax_b2=-gmax_b2-1 ; gmin_b2=-gmin_b2-1
             end if
           end if

           if( testp==1 .and. testm==1)then
             gmin_b=min(gmin_b1,gmin_b2) ; gmax_b=max(gmax_b1,gmax_b2)
           else if( testp==1 )then
             gmin_b=gmin_b1 ; gmax_b=gmax_b1
           else if( testm==1 )then
             gmin_b=gmin_b2 ; gmax_b=gmax_b2
           end if

         end if ! Endif take into account time-reversal symmetry

         if (igb+1>2*mgfft+4) then
           write(message, '(2a, 4(a,3(i0,1x),a))' )&
             "About to overwrite gbound array (FFT mesh too small) ",ch10, &
             "   iloop, igb, mgb = ",iloop,igb,2*mgfft+4, ch10, &
             "   istwfk, mgfft, npw = ",istwf_k, mgfft, npw, ch10, &
             "   minval(kg_k) = ",minval(kg_k, dim=2), ch10, &
             "   maxval(kg_k) = ",maxval(kg_k, dim=2), ch10
           MSG_BUG(message)
         end if

         gbound(igb,ii)=gmin_b
         gbound(igb+1,ii)=gmax_b

         if( iloop==1 .and. istwf_k>=2 .and. istwf_k<=5 .and. fftalgc==2 .and. g_a==0)then
!          If k_y=0 , for fftalgc==2, the g_a==0 plane must be completed
           if(istwf_k==2 .or. istwf_k==4)then
             gbound(igb+1,ii)=max(gmax_b,-gmin_b)
             gbound(igb,ii)=-gbound(igb+1,ii)
           else if(istwf_k==3 .or. istwf_k==5)then
             gbound(igb+1,ii)=max(gmax_b,-gmin_b-1)
             gbound(igb,ii)=-gbound(igb+1,ii)-1
           end if

         end if

         igb=igb+2

       end do ! g_a
     end if
   end do  ! iloop
 end do ! ii (fftalgc)

!DEBUG
!write(std_out,'(a)')' sphereoundary : list of plane waves coordinates for 1st k point '
!write(std_out,'(a)')'       ipw       kg_k(1:3,ipw) '
!do ipw=1,npw
!write(std_out, '(i10,a,3i6)' )ipw,'  ',kg_k(1:3,ipw)
!end do
!write(std_out, '(a)' )' sphereboundary : list of boundaries '
!do igb=1,2*mgfft+8
!write(std_out, '(i10,a,2i6)' )igb,'  ',gbound(igb,1),gbound(igb,2)
!end do
!write(std_out,*)' sphereboundary : exit '
!ENDDEBUG

end subroutine sphereboundary
!!***
