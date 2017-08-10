!{\src2tex{textfont=tt}}
!!****f* ABINIT/predictimg
!! NAME
!! predictimg
!!
!! FUNCTION
!! Given the past history of images, predict the new set of images
!!
!! COPYRIGHT
!! Copyright (C) 2009-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! deltae=averaged energy difference used to control convergence over images
!! imagealgo_str=name of the algorithm (with images) used
!! imgmov=gives the algorithm to be used for prediction of new set of images
!! itimimage=time index for image propagation (itimimage+1 is to be predicted here)
!! itimimage_eff=time index in the history
!! list_dynimage(nimage)=list of dynamical images. The non-dynamical ones will not change.
!!       Example : in the NEB method, or in the string method, one expect the two end images to be fixed.
!! mep_param=several parameters for Minimal Energy Path (MEP) search
!! mpi_enreg=MPI-parallelisation information
!! natom= number of atoms
!! ndynimage=number of dynamical images
!! nimage=number of images (treated by current proc)
!! nimage_tot=total number of images
!! ntimimage_stored=number of time steps stored in the history
!! pimd_param=several parameters for Path-Integral MD
!! prtvolimg=printing volume
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
!!      gstateimg
!!
!! CHILDREN
!!      predict_copy,predict_ga,predict_neb,predict_pimd,predict_steepest
!!      predict_string,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine predictimg(deltae,imagealgo_str,imgmov,itimimage,itimimage_eff,list_dynimage,&
&                     ga_param,mep_param,mpi_enreg,natom,ndynimage,nimage,nimage_tot,&
&                     ntimimage_stored,pimd_param,prtvolimg,results_img)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_mep
 use m_ga
 use m_pimd
 use m_results_img
 use m_use_ga

 use m_results_gs , only : results_gs_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'predictimg'
 use interfaces_14_hidewrite
 use interfaces_45_geomoptim, except_this_one => predictimg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: imgmov,itimimage,itimimage_eff,natom,ndynimage
 integer,intent(in) :: nimage,nimage_tot,ntimimage_stored,prtvolimg
 character(len=60),intent(in) :: imagealgo_str
 real(dp),intent(in) :: deltae
 type(mep_type),intent(inout) :: mep_param
 type(ga_type),intent(inout) :: ga_param
 type(pimd_type),intent(inout) :: pimd_param
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 type(results_img_type) :: results_img(nimage,ntimimage_stored)

!Local variables-------------------------------
!scalars
 integer,save :: idum=5
 logical :: is_pimd
 character(len=500) :: msg
!arrays

! *************************************************************************

 is_pimd=(imgmov==9.or.imgmov==10.or.imgmov==13)

!Write convergence info
 write(msg,'(3a)') ch10,&
& '------------------------------------------------------------',ch10
 if (prtvolimg<2) write(msg,'(5a)') trim(msg),' ',trim(imagealgo_str),':',ch10

!Specific case of 4th-order RK algorithm
 if (mep_param%mep_solver==4) then
   if (mod(itimimage,4)==0) then
     write(msg,'(2a,i1,2a)') trim(msg),&
&     ' Fourth-order Runge-Kutta algorithm - final step',ch10
     if (itimimage>4) write(msg,'(2a,es11.3,2a)') trim(msg),&
&     ' Average[Abs(Etotal(t)-Etotal(t-dt))]=',deltae,' Hartree',ch10
     write(msg,'(2a)') trim(msg),' Moving images of the cell...'
   else
     write(msg,'(2a,i1,2a)') trim(msg),&
&     ' Fourth-order Runge-Kutta algorithm - intermediate step ',mod(itimimage,4),ch10
     write(msg,'(2a)') trim(msg),' Computing new intermediate positions...'
   end if
 else if (is_pimd) then

!  PIMD
   write(msg,'(2a)') trim(msg),' Moving images of the cell...'
 else

!  Other cases
   if (itimimage>1) write(msg,'(2a,es11.3,2a)') trim(msg),&
&   ' Average[Abs(Etotal(t)-Etotal(t-dt))]=',deltae,' Hartree',ch10
   write(msg,'(2a)') trim(msg),' Moving images of the cell...'

 end if

 call wrtout(ab_out ,msg,'COLL')
 call wrtout(std_out,msg,'COLL')


 select case(imgmov)

 case(0)

   call predict_copy(itimimage,itimimage_eff,list_dynimage,ndynimage,nimage,&
&   ntimimage_stored,results_img)

 case(1)

   call predict_steepest(itimimage,itimimage_eff,list_dynimage,mep_param,natom,ndynimage,nimage,&
&   ntimimage_stored,results_img)

 case(2)

   call predict_string(itimimage,itimimage_eff,list_dynimage,mep_param,mpi_enreg,natom,&
&   ndynimage,nimage,nimage_tot,ntimimage_stored,results_img)

 case(4)

   call predict_ga(itimimage,itimimage_eff,idum,ga_param,natom,nimage,&
&   ntimimage_stored,results_img)

 case(5)

   call predict_neb(itimimage,itimimage_eff,list_dynimage,mep_param,mpi_enreg,natom,&
&   ndynimage,nimage,nimage_tot,ntimimage_stored,results_img)

 case(9, 10, 13)
!    Path Integral Molecular Dynamics
   call predict_pimd(imgmov,itimimage,itimimage_eff,mpi_enreg,natom,nimage,nimage_tot,&
&   ntimimage_stored,pimd_param,prtvolimg,results_img)

 case default

 end select

end subroutine predictimg
!!***
