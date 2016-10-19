!{\src2tex{textfont=tt}}
!!****f* ABINIT/predictimg
!! NAME
!! predictimg
!!
!! FUNCTION
!! Given the past history of images, predict the new set of images
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! deltae=averaged energy difference used to control convergence over images
!! imagealgo_str=name of the algorithm (with images) used
!! imgmov=gives the algorithm to be used for prediction of new set of images
!! itimimage=number of the current time for image propagation (itimimage+1 is to be predicted here)
!! list_dynimage(nimage)=list of dynamical images. The non-dynamical ones will not change.
!!       Example : in the NEB method, or in the string method, one expect the two end images to be fixed.
!! mep_param=several parameters for Minimal Energy Path (MEP) search
!! mpi_enreg=MPI-parallelisation information
!! natom= number of atoms
!! ndynimage=number of dynamical images
!! nimage=number of images (treated by current proc)
!! nimage_tot=total number of images
!! ntimimage=dimension of several arrays
!! pimd_param=several parameters for Path-Integral MD
!! prtvolimg=printing volume
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! results_img(ntimimage,nimage)=datastructure that hold all the history of previous computations.
!!   results_img(:,:)%acell(3)
!!    at input, history of the values of acell for all images, up to itimimage
!!    at output, the predicted values of acell for all images
!!   results_img(:,:)%results_gs
!!    at input, history of the values of energies and forces for all images, up to itimimage
!!   results_img(:,:)%rprim(3,3)
!!    at input, history of the values of rprim for all images, up to itimimage
!!    at output, the predicted values of rprim for all images
!!   results_img(:,:)%vel(3,natom)
!!    at input, history of the values of vel for all images, up to itimimage
!!    at output, the predicted values of vel for all images
!!   results_img(:,:)%vel_cell(3,3)
!!    at input, history of the values of vel_cell for all images, up to itimimage
!!    at output, the predicted values of vel_cell for all images
!!   results_img(:,:)%xred(3,natom)
!!    at input, history of the values of xred for all images, up to itimimage
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


subroutine predictimg(deltae,imagealgo_str,imgmov,itimimage,list_dynimage,&
&                     ga_param,mep_param,mpi_enreg,natom,ndynimage,nimage,nimage_tot,ntimimage,&
&                     pimd_param,prtvolimg,results_img)

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
 integer,intent(in) :: imgmov,itimimage,natom,ndynimage,nimage,nimage_tot,ntimimage,prtvolimg
 character(len=60),intent(in) :: imagealgo_str
 real(dp),intent(in) :: deltae
 type(mep_type),intent(inout) :: mep_param
 type(ga_type),intent(inout) :: ga_param
 type(pimd_type),intent(inout) :: pimd_param
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: list_dynimage(ndynimage)
 type(results_img_type) :: results_img(nimage,ntimimage)

!Local variables-------------------------------
!scalars
 integer,save :: idum=5
 character(len=500) :: msg
!arrays

! *************************************************************************

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

   call predict_copy(itimimage,list_dynimage,ndynimage,nimage,ntimimage,results_img)

 case(1)

   call predict_steepest(itimimage,list_dynimage,mep_param,natom,ndynimage,nimage,&
&   ntimimage,results_img)

 case(2)

   call predict_string(itimimage,list_dynimage,mep_param,mpi_enreg,natom,&
&   ndynimage,nimage,nimage_tot,ntimimage,results_img)

 case(4)

   call predict_ga(itimimage,idum,ga_param,natom,nimage,ntimimage,results_img)

 case(5)

   call predict_neb(itimimage,list_dynimage,mep_param,mpi_enreg,natom,&
&   ndynimage,nimage,nimage_tot,ntimimage,results_img)

 case(9, 10, 13)
!    Path Integral Molecular Dynamics
   call predict_pimd(imgmov,itimimage,mpi_enreg,natom,nimage,nimage_tot,&
&   ntimimage,pimd_param,prtvolimg,results_img)

 case default

 end select

end subroutine predictimg
!!***
