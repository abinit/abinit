!!****m* ABINIT/predict_ga
!! NAME
!! predict_ga
!!
!! FUNCTION
!! Given a given set of images, which represent a population, it predicts a new set of images.
!! The implementation is based on a Genetic Algorithm idea, where the best fit candidates are passed
!! to the next generation. Those are chosen from ga_opt_percent% best fit and (1-ga_opt_percent)% from Genetic rules
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT group (XG, AHR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! itimimage=time index for image propagation (itimimage+1 is to be predicted here)
!! itimimage_eff=time index in the history
!! list_dynimage(nimage)=list of dynamical images.
!! This is quite useful when ground states of the A and B states is known
!! natom=dimension of vel_timimage and xred_timimage
!! ndynimage=number of dynamical images
!! nimage= population size
!! ntimimage_stored=number of time steps stored in the history
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! results_img(ntimimage_stored,nimage)=datastructure that holds the history of previous computations.
!!   results_img(:,:)%acell(3)
!!    at input, history of the values of acell for all images
!!    at output, the predicted values of acell for all images
!!   results_img(:,:)%results_gs
!!    at input, history of the values of energies and forces for all images
!!   results_img(:,:)%rprim(3,3)
!!    at input, history of the values of rprim for all images
!!    at output, the predicted values of rprim for all images
!!   results_img(:,:)%vel(3,natom)
!!    at input, history of the values of vel for all images
!!    at output, the predicted values of vel for all images
!!   results_img(:,:)%vel_cell(3,3)
!!    at input, history of the values of vel_cell for all images
!!    at output, the predicted values of vel_cell for all images
!!   results_img(:,:)%xred(3,natom)
!!    at input, history of the values of xred for all images
!!    at output, the predicted values of xred for all images
!!
!! PARENTS
!!      predictimg
!!
!! CHILDREN
!!      convert_coortogen,convert_gentocoor,initialize_perm,metric,mkradim
!!      mkrdim,randomize_parent,sort_dp,swap,symanal,symfind,symlatt
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_use_ga

 use defs_basis
 use m_abicore
 use m_ga
 use m_sort

 use m_symfind,        only : symfind, symanal, symlatt
 use m_geometry,       only : mkradim, mkrdim, metric, dist2
 use m_results_img,    only : results_img_type,gather_array_img
 use m_numeric_tools,  only : uniformrandom
 use m_symtk,          only : matr3inv
 implicit none

 private

 public :: predict_ga
 public :: checksymmetrygroup


CONTAINS

subroutine predict_ga(itimimage_eff,idum,ga_param,natom,nimage,&
&                     ntimimage_stored,results_img)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)     :: itimimage_eff,natom,nimage,ntimimage_stored
 integer,intent(inout)  :: idum
!arrays
 type(results_img_type) :: results_img(nimage,ntimimage_stored)
 type(ga_type),intent(inout) :: ga_param

!Local variables------------------------------
!scalars
 integer :: ii,jj,kk,itmp,iimage,indiv,oper,iparent1,iparent2,ndimen
 integer :: next_itimimage,nsurvivor,nspinor
!character(len=500)   :: message
!arrays
 integer,allocatable  :: ieperm(:),ihperm(:),ibgoptperm(:,:),ibgperm(:,:)
 integer,allocatable  :: zperm1(:),zperm2(:)
 real(dp) :: gprimd(3,3),rmet(3,3),gmet(3,3)
 real(dp) :: rprimdparent1(3,3),rprimdparent2(3,3)
 real(dp) :: rprimson1(3,3),rprimson2(3,3)
 real(dp) :: rprimdson1(3,3),rprimdson2(3,3)
 real(dp) :: acellson1(3),acellson2(3)
 real(dp) :: son1(3,natom),son2(3,natom)
 real(dp) :: parent1(3,natom),parent2(3,natom)
! store the energy and the enthalpy of all elements of the population
 real(dp),allocatable :: etotal_img(:),enthalpy_img(:),bg_img(:,:),bg_opt_img(:,:)
 real(dp),allocatable :: acell(:,:),acell_old(:,:),rprim(:,:,:),rprim_old(:,:,:),rprimd(:,:,:)
 real(dp),allocatable :: fitness(:),zcoor1(:),zcoor2(:)
 real(dp),allocatable :: coor(:,:),coor_old(:,:)
 real(dp),allocatable :: vson1(:),vson2(:)

!real quantities
 real(dp) :: denom,sumH,Hmin,Hmax,ucvol,rtmp(3,3)

! *************************************************************************

!DEBUG
!write(std_out,*)' MODULE predict_ga : enter '
!ENDDEBUG

! gen dimension

 ndimen=3*natom
 nspinor=results_img(nimage,itimimage_eff)%nsppol

 ABI_ALLOCATE(coor,(ndimen,nimage))
 ABI_ALLOCATE(coor_old,(ndimen,nimage))
 ABI_ALLOCATE(acell,(3,nimage))
 ABI_ALLOCATE(acell_old,(3,nimage))
 ABI_ALLOCATE(rprim,(3,3,nimage))
 ABI_ALLOCATE(rprimd,(3,3,nimage))
 ABI_ALLOCATE(rprim_old,(3,3,nimage))

 ABI_ALLOCATE(zperm1,(natom))
 ABI_ALLOCATE(zperm2,(natom))
 ABI_ALLOCATE(ieperm,(nimage))
 ABI_ALLOCATE(ihperm,(nimage))
 ABI_ALLOCATE(ibgperm,(nspinor,nimage))
 ABI_ALLOCATE(ibgoptperm,(nspinor,nimage))

 ABI_ALLOCATE(etotal_img,(nimage))
 ABI_ALLOCATE(enthalpy_img,(nimage))

!to store locally the calculated band gaps
 ABI_ALLOCATE(bg_img,(nspinor,nimage))
 ABI_ALLOCATE(bg_opt_img,(nspinor,nimage))

 ABI_ALLOCATE(fitness,(nimage))

 ABI_ALLOCATE(zcoor1,(natom))
 ABI_ALLOCATE(zcoor2,(natom))
 ABI_ALLOCATE(vson1,(ndimen))
 ABI_ALLOCATE(vson2,(ndimen))

 call initialize_perm(ieperm,nimage)
 call initialize_perm(ihperm,nimage)
 do ii=1,nspinor
   call initialize_perm(ibgperm(ii,:),nimage)
   call initialize_perm(ibgoptperm(ii,:),nimage)
 enddo

!from gs_results_image take energies, rprim and acell and build energy and enthalpy vectors reordered from larger to smaller

 do iimage=1,nimage
   etotal_img(iimage)=results_img(iimage,itimimage_eff)%results_gs%etotal
   call convert_coortogen(results_img(iimage,itimimage_eff)%xred(:,:),coor_old(:,iimage),natom)
   acell_old(:,iimage)=results_img(iimage,itimimage_eff)%acell(:)
   rprim_old(:,:,iimage)=results_img(iimage,itimimage_eff)%rprim(:,:)
   if (results_img(iimage,itimimage_eff)%results_gs%gaps(3,1) > 0.0_dp) then
       bg_img(:,iimage)=results_img(iimage,itimimage_eff)%results_gs%gaps(1,:)
       bg_opt_img(:,iimage)=results_img(iimage,itimimage_eff)%results_gs%gaps(2,:)
   else
       bg_img(:,iimage)=-100.0_dp
       bg_opt_img(:,iimage)=-100.0_dp
   endif
   call mkrdim(acell_old(:,iimage),rprim_old(:,:,iimage),rprimd(:,:,iimage))
   call metric(gmet,gprimd,-1,rmet,rprimd(:,:,iimage),ucvol)
   enthalpy_img(iimage)=etotal_img(iimage) &
&    +sum(results_img(iimage,itimimage_eff)%results_gs%strten(1:3))*ucvol
 enddo

!sort energies

 call sort_dp(nimage,etotal_img,ieperm,tol9)

!sort enthalpies

 call sort_dp(nimage,enthalpy_img,ihperm,tol9)

!sort band gaps

 do ii=1,nspinor
   call sort_dp(nimage,bg_img(ii,:),ibgperm(ii,:),tol9)
   call sort_dp(nimage,bg_opt_img(ii,:),ibgoptperm(ii,:),tol9)
 enddo

! Fitness is calculated

 sumH=0.d0
 Hmin=minval(enthalpy_img)
 Hmax=maxval(enthalpy_img)
 fitness = zero

 select case (ga_param%ga_fitness)
 case ( 1 )

! Function weighted over the difference with respect to the minimum

    do iimage=1,nimage
        sumh=sumh+enthalpy_img(iimage); fitness(iimage)=sumh
    enddo
    denom=nimage*enthalpy_img(nimage)-sumh
    do iimage=1,nimage
        fitness(iimage)=(iimage*enthalpy_img(nimage)-fitness(iimage))/denom
    enddo

 case  ( 2 )

! Botzmann like, using the enthalpy difference.

   do iimage=1,nimage
      fitness(iimage)=exp(-one*(enthalpy_img(iimage)-enthalpy_img(1)))
      sumH = sumH + fitness(iimage)
   enddo
   fitness = fitness/sumH
   do iimage=2,nimage
     fitness(iimage)=fitness(iimage)+fitness(iimage-1)
   enddo

  case ( 3 )

! weighted ove the position in the ordered list, with probability 1/i (non uniform)

   do iimage=1,nimage
      fitness(iimage)=one/float(iimage)
      sumH = sumH + fitness(iimage)
   enddo
   fitness = fitness/sumH
   do iimage=2,nimage
     fitness(iimage)=fitness(iimage)+fitness(iimage-1)
   enddo

  end select

!  do a single GA boocle

 indiv=0

! Selection over the best ga_opt_percent of the population

 nsurvivor=int(ga_param%ga_opt_percent*nimage)

 if (nsurvivor < one) nsurvivor=1

! pass coordinates,rprim, acell of survivors to next generation

 do iimage=1,nsurvivor
   indiv=indiv+1
   coor(:,iimage)=coor_old(:,ihperm(iimage))
   acell(:,iimage)=acell_old(:,ihperm(iimage))
   rprim(:,:,iimage)=rprim_old(:,:,ihperm(iimage))
 enddo

! complete the number of individuals of the generation by choosing them through GA

 do while(indiv<nimage)

! ga_n_rules corresponds to the number of chosen Genetic rules

   oper=ga_param%ga_rules(int(ga_param%ga_n_rules*uniformrandom(idum)+1))

   select case(oper)
     case(1) ! cut and splice elements after a randomization
       iparent1=choosefather(fitness,nimage,idum)
       iparent2=choosefather(fitness,nimage,idum)
       call convert_gentocoor(parent1,coor_old(:,ihperm(iparent1)),natom)
       call convert_gentocoor(parent2,coor_old(:,ihperm(iparent2)),natom)
! randomize the cell: random rotation and translation
       call randomize_parent(parent1,natom,idum)
       call randomize_parent(parent2,natom,idum)
! choose direction of cutting plane
       itmp = int(3*uniformrandom(idum)+1)
! order coordinates from small to large along that random direction axis of both parents
       zcoor1(:)=parent1(itmp,:)
       zcoor2(:)=parent2(itmp,:)
       call initialize_perm(zperm1,natom)
       call initialize_perm(zperm2,natom)
       call sort_dp(natom,zcoor1,zperm1,tol9)
       call sort_dp(natom,zcoor2,zperm2,tol9)
! choose the atom position to take the cut and take atoms below from one parent and above from the other parent
       itmp=int(natom*uniformrandom(idum)+1)
       do ii=1,itmp
          son1(:,ii)=parent1(:,zperm1(ii))
          son2(:,ii)=parent2(:,zperm2(ii))
       enddo
       do ii=itmp+1,natom
          son1(:,ii)=parent2(:,zperm2(ii))
          son2(:,ii)=parent1(:,zperm1(ii))
       enddo
! random combination of rprimd from parents
       call mkrdim(acell_old(:,ihperm(iparent1)),rprim_old(:,:,ihperm(iparent1)),rprimdparent1)
       call mkrdim(acell_old(:,ihperm(iparent2)),rprim_old(:,:,ihperm(iparent2)),rprimdparent2)
       do ii=1,3
         do jj=1,3
           rtmp(ii,jj)=uniformrandom(idum)
         enddo
       enddo
       rprimdson1=(one-rtmp)*rprimdparent1+rtmp*rprimdparent2
       do ii=1,3
         do jj=1,3
           rtmp(ii,jj)=uniformrandom(idum)
         enddo
       enddo
       rprimdson2=(one-rtmp)*rprimdparent1+rtmp*rprimdparent2
!create acell and rprim from rprimd of springoffs
       call mkradim(acellson1,rprimson1,rprimdson1)
       call mkradim(acellson2,rprimson2,rprimdson2)
! check distances of atoms of every springoff
       if (checkatomicdist(natom,son1,rprimdson1)==0 .and. indiv<nimage) then
         indiv=indiv+1
!if any fix coordinate restore the parent coordinate without modifications
         do ii=1,natom
           do jj=1,3
            if (ga_param%ga_iatfix(jj,ii) == 1) son1(jj,ii)=parent1(jj,ii)
           enddo
         enddo
         call convert_coortogen(son1,coor(:,indiv),natom)
         acell(:,indiv)=acellson1
         rprim(:,:,indiv)=rprimson1
       endif
       if (checkatomicdist(natom,son2,rprimdson2)==0 .and. indiv<nimage) then
         indiv=indiv+1
!if any fix coordinate restore the parent coordinate without modifications
         do ii=1,natom
           do jj=1,3
            if (ga_param%ga_iatfix(jj,ii) == 1) son2(jj,ii)=parent1(jj,ii)
           enddo
         enddo
         call convert_coortogen(son2,coor(:,indiv),natom)
         acell(:,indiv)=acellson2
         rprim(:,:,indiv)=rprimson2
       endif
     case(2)! vector flip mutation
       iparent1=choosefather(fitness,nimage,idum)
       ii=int(ndimen*uniformrandom(idum)+1)
       jj=int(ndimen*uniformrandom(idum)+1)
       if (ii>jj) then
         call swap(ii,jj)
       end if
       vson1(1:ii)=coor_old(1:ii,ihperm(iparent1))
       if (jj<ndimen) vson1(jj+1:ndimen)=coor_old(jj+1:ndimen,ihperm(iparent1))
       do kk=ii,jj
         vson1(kk)=coor_old(jj+1-kk,ihperm(iparent1))
       enddo
       rprimson1=rprim_old(:,:,ihperm(iparent1))
       acellson1=acell_old(:,ihperm(iparent1))
       call convert_gentocoor(son1,vson1,natom)
       call mkrdim(acellson1,rprimson1,rprimdson1)
!if any fix coordinate restore the parent coordinate without modifications
       do ii=1,natom
         do jj=1,3
          if (ga_param%ga_iatfix(jj,ii) == 1) son1(jj,ii)=parent1(jj,ii)
         enddo
       enddo
       call convert_coortogen(son1,vson1,natom)
       if (checkatomicdist(natom,son1,rprimdson1)==0 .and. indiv<nimage) then
         indiv=indiv+1
         coor(:,indiv)=vson1
         acell(:,indiv)=acellson1
         rprim(:,:,indiv)=rprimson1
       endif
     case(3) ! random strain - nondiagonal
       iparent1=choosefather(fitness,nimage,idum)
       rtmp(1,1)=one+gaussian_random(idum,0.1_dp)
       rtmp(2,2)=one+gaussian_random(idum,0.1_dp)
       rtmp(3,3)=one+gaussian_random(idum,0.1_dp)
       rtmp(1,2)=gaussian_random(idum,0.1_dp)*half
       rtmp(1,3)=gaussian_random(idum,0.1_dp)*half
       rtmp(2,3)=gaussian_random(idum,0.1_dp)*half
       rtmp(2,1)=rtmp(1,2)
       rtmp(3,1)=rtmp(1,3)
       rtmp(3,2)=rtmp(2,3)
       rprimson1=matmul(rtmp,rprim_old(:,:,ihperm(iparent1)))
       vson1=coor_old(:,ihperm(iparent1))
       acellson1=acell_old(:,ihperm(iparent1))
       call convert_gentocoor(son1,vson1,natom)
       call mkrdim(acellson1,rprimson1,rprimdson1)
!if any fix coordinate restore the parent coordinate without modifications
       do ii=1,natom
         do jj=1,3
          if (ga_param%ga_iatfix(jj,ii) == 1) son1(jj,ii)=parent1(jj,ii)
         enddo
       enddo
       call convert_coortogen(son1,vson1,natom)
       if (checkatomicdist(natom,son1,rprimdson1)==0 .and. indiv<nimage) then
         indiv=indiv+1
         coor(:,indiv)=vson1
         acell(:,indiv)=acellson1
         rprim(:,:,indiv)=rprimson1
       endif
     case(4) ! coordinates mutation
       iparent1=choosefather(fitness,nimage,idum)
       vson1=coor_old(:,ihperm(iparent1))
       itmp=ndimen/4
       if (itmp<1) itmp=1
       do jj=1,itmp
         ii=int(ndimen*uniformrandom(idum)+1)
         vson1(ii)=vson1(ii)+0.15*uniformrandom(idum)
         if (vson1(ii)>one) vson1(ii)=vson1(ii)-one
         if (vson1(ii)<zero) vson1(ii)=vson1(ii)+one
       enddo
       rprimson1=rprim_old(:,:,ihperm(iparent1))
       acellson1=acell_old(:,ihperm(iparent1))
       call convert_gentocoor(son1,vson1,natom)
       call mkrdim(acellson1,rprimson1,rprimdson1)
!if any fix coordinate restore the parent coordinate without modifications
       do ii=1,natom
         do jj=1,3
          if (ga_param%ga_iatfix(jj,ii) == 1) son1(jj,ii)=parent1(jj,ii)
         enddo
       enddo
       call convert_coortogen(son1,vson1,natom)
       if (checkatomicdist(natom,son1,rprimdson1)==0 .and. indiv<nimage) then
         indiv=indiv+1
         coor(:,indiv)=vson1
         acell(:,indiv)=acellson1
         rprim(:,:,indiv)=rprimson1
       endif
    end select
  enddo

 next_itimimage=itimimage_eff+1
 if (next_itimimage>ntimimage_stored) next_itimimage=1

 do iimage=1,nimage
   results_img(iimage,next_itimimage)%acell = acell(:,iimage)
   results_img(iimage,next_itimimage)%rprim = rprim(:,:,iimage)
   results_img(iimage,next_itimimage)%vel = results_img(iimage,itimimage_eff)%vel
   results_img(iimage,next_itimimage)%vel_cell = results_img(iimage,itimimage_eff)%vel_cell
   call convert_gentocoor(results_img(iimage,next_itimimage)%xred,coor(:,iimage),natom)
 enddo

 ABI_DEALLOCATE(coor)
 ABI_DEALLOCATE(acell)
 ABI_DEALLOCATE(acell_old)
 ABI_DEALLOCATE(rprim)
 ABI_DEALLOCATE(rprimd)
 ABI_DEALLOCATE(rprim_old)
 ABI_DEALLOCATE(zcoor1)
 ABI_DEALLOCATE(zcoor2)
 ABI_DEALLOCATE(coor_old)
 ABI_DEALLOCATE(zperm1)
 ABI_DEALLOCATE(zperm2)
 ABI_DEALLOCATE(ieperm)
 ABI_DEALLOCATE(ihperm)
 ABI_DEALLOCATE(ibgperm)
 ABI_DEALLOCATE(ibgoptperm)
 ABI_DEALLOCATE(etotal_img)
 ABI_DEALLOCATE(bg_img)
 ABI_DEALLOCATE(bg_opt_img)
 ABI_DEALLOCATE(enthalpy_img)
 ABI_DEALLOCATE(fitness)
 ABI_DEALLOCATE(vson1)
 ABI_DEALLOCATE(vson2)

end subroutine predict_ga

!!*************** local subroutines

INTEGER FUNCTION choosefather(fitf,n,idum)

 implicit none

 integer,intent(in) :: n
 integer,intent(inout) :: idum
 real(dp), dimension(:), intent(in) :: fitf

 real(dp) :: x1
 integer :: ii

 x1=uniformrandom(idum);
 choosefather=1
 do ii=2,n
    if(fitf(ii-1)<x1.and.x1<=fitf(ii))then
       choosefather=ii
       exit
    endif
 enddo
end FUNCTION choosefather

!!

SUBROUTINE swap(a,b)

 integer, intent(inout) :: a,b
 integer :: dum
   dum=a; a=b;  b=dum
END SUBROUTINE swap

!!

INTEGER FUNCTION comp_indiv(distances,indiv,natom,nimage)

!! comparing individuals from the same generation and check they
!! are not too close.  We compare all individuals with individual: indiv.
!! We assume, all distances for a given individual are ordered and
!! we define a metric from the total difference distances between individuals.
!! if comp_indiv is 0, means that two individuals are two close.
 implicit none

 integer, intent(in) :: indiv,natom,nimage
 real(dp),intent(in) ::distances(natom,nimage)
 real(dp) :: diff

 integer :: ii,jj

 comp_indiv=1

 do ii=1,indiv-1
   diff=0.0_dp
   do jj=1,natom
     diff=(distances(jj,indiv)-distances(jj,ii))**2
   enddo
   if (diff.le.0.001_dp) comp_indiv=0
 enddo

end FUNCTION comp_indiv

!!

SUBROUTINE randomize_parent(parent,natom,idum)

! take a parent and randomize the positions of the atoms
implicit none

integer,intent(in)     :: natom
integer,intent(inout)     :: idum
real(dp),intent(inout) :: parent(3,natom)

real(dp)               :: tmp(3),rot(3,3)
integer                :: ii,jj

!random rotation

do ii=1,natom
  if (uniformrandom(idum)>half) then
    do jj=1,3
      tmp(jj)=uniformrandom(idum)*two_pi
    enddo
    rot(1,1)=cos(tmp(3))*cos(tmp(2))
    rot(1,2)=cos(tmp(3))*sin(tmp(1))*sin(tmp(2))-cos(tmp(1))*sin(tmp(3))
    rot(1,3)=cos(tmp(1))*cos(tmp(3))*sin(tmp(2))+sin(tmp(1))*sin(tmp(3))
    rot(2,1)=cos(tmp(2))*sin(tmp(3))
    rot(2,2)=cos(tmp(1))*cos(tmp(3))+sin(tmp(1))*sin(tmp(2))*sin(tmp(3))
    rot(2,3)=cos(tmp(1))*sin(tmp(2))*sin(tmp(3))-cos(tmp(3))*sin(tmp(1))
    rot(3,1)=-sin(tmp(2))
    rot(3,2)=cos(tmp(2))*sin(tmp(1))
    rot(3,3)=cos(tmp(1))*cos(tmp(2))
    parent(:,ii)=matmul(rot(:,:),parent(:,ii))
  endif

! random reflection
  if (uniformrandom(idum)>half) then
      parent(:,ii)=-parent(:,ii)
  endif

! fold back all coordinates to reduced coordinates in cube [0,1]
  parent(:,ii)=parent(:,ii)-anint((parent(:,ii)-one)/two)

enddo

! shift all atoms along a random direction

do ii=1,3
   tmp(ii)=uniformrandom(idum)
enddo

do ii=1,natom
   parent(:,ii)=parent(:,ii)+tmp(:)
   parent(:,ii)=parent(:,ii)-anint(parent(:,ii)/two)
enddo

end SUBROUTINE randomize_parent

!!

SUBROUTINE convert_gentocoor(parent,coor,natom)

! convert gene (single vector) to coordinates (3,natom)
implicit none

integer,intent(in) :: natom
double precision,intent(out)    :: parent(3,natom)
double precision,intent(in) :: coor(3*natom)

integer            :: ii,jj


do ii=1,natom
  jj=(ii-1)*3
  parent(:,ii)=coor(jj+1:jj+3)
enddo

end SUBROUTINE convert_gentocoor

!!

SUBROUTINE convert_coortogen(parent,coor,natom)

! convert coordinates in gene notation (single vector with all coordinates)
implicit none

integer,intent(in) :: natom
double precision,intent(in)    :: parent(3,natom)
double precision,intent(out) :: coor(3*natom)

integer            :: ii,jj

do ii=1,natom
  jj=(ii-1)*3
  coor(jj+1:jj+3)=parent(:,ii)
enddo

end SUBROUTINE convert_coortogen

!!

SUBROUTINE initialize_perm(iperm,nimage)

!! initialize the vector iperm with corresponding indices.
 implicit none

 integer, intent(in) :: nimage
 integer, intent(inout) :: iperm(nimage)
 integer :: ii


 do ii=1,nimage
   iperm(ii)=ii
 enddo

end SUBROUTINE initialize_perm

! if after a genetic rule, check if two atoms in the same gen
! are too close


INTEGER FUNCTION checkatomicdist(natom,coord,rprimd)

!! check if two atoms are two close.
!! if they are, checkatomicdist=0
 implicit none

 integer, intent(in) :: natom
 real(dp), intent(in) :: coord(3,natom),rprimd(3,3)
 real(dp) :: d,v1(3),v2(3)
 integer :: ii,jj

 checkatomicdist=0
 do ii=1,natom
   v1=coord(:,ii)
   do jj=ii+1,natom
     v2=coord(:,jj)
     d=dist2(v1,v2,rprimd,1)
     if (d<1.0d0) then
       checkatomicdist=1
       EXIT
     endif
   enddo
 enddo

END FUNCTION checkatomicdist

SUBROUTINE checksymmetrygroup(rprimd,xred,typat,msym,natom,ptgroupma,spgroup,symrel_out,tnons_out)

  implicit none

!Arguments ------------------------------------
!scalars
  integer,intent(in) :: msym,natom
  integer,intent(in)  :: typat(natom)
  integer,intent(out) :: ptgroupma,spgroup
! Arrays
  real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
  integer,optional,intent(out) :: symrel_out(3,3,msym)
  real(dp),optional,intent(out) :: tnons_out(3,msym)

!Local variables ---------------------------------------
!scalars
  integer :: berryopt,jellslab=0,noncoll=0,nptsym,nzchempot=0,use_inversion
  integer :: chkprim,nsym
! Arrays
  integer :: bravais(11),ptsymrel(3,3,msym)
  integer :: symafm(msym),symrel(3,3,msym)
  real(dp) :: efield(3),gprimd(3,3),spinat(3,natom)
  real(dp) :: tnons(3,msym)
  real(dp) :: genafm(3)

! given the acel, rprim and coor
! this suroutine find the symmetry group

write(std_out,*) 'symlatt'
  call symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tol3)

call matr3inv(rprimd,gprimd)

write(std_out,*) 'symfind'
  call symfind(berryopt,efield,gprimd,jellslab,msym,natom,noncoll,nptsym,nsym,&
&           nzchempot,0,ptsymrel,spinat,symafm,symrel,tnons,tol3,typat,use_inversion,xred)

write(std_out,*) 'symanal'
  call symanal(bravais,chkprim,genafm,msym,nsym,ptgroupma,rprimd,spgroup,symafm,symrel,tnons,tol3)

  if(present(symrel_out))symrel_out = symrel 
  if(present(tnons_out)) tnons_out  = tnons


END SUBROUTINE checksymmetrygroup

DOUBLE PRECISION FUNCTION gaussian_random(idum,sigma)

  implicit none

  integer,intent(inout) :: idum
  real(dp), intent(in) :: sigma

  real(dp) :: r1,r2,w
!! The polar form of the Box-Muller transformation
  w=two
  do while (w >= one .or. w == zero)
       r1=two*uniformrandom(idum)-one
       r2=two*uniformrandom(idum)-one
       w=r1*r1+r2*r2
  enddo
  w=sqrt(-two*log(w)/w)
  gaussian_random=r2*w*sigma

END FUNCTION gaussian_random

end MODULE m_use_ga
