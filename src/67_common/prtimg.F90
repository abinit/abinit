!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtimg
!! NAME
!! prtimg
!!
!! FUNCTION
!! Print out results obtained by as ground-state calculation of
!! an image of the cell. The printing format is condensed in order
!! to facilitate the reading.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2017 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dynimage(nimagetot)=flags defining static/dynamic state of images
!!  imagealgo_str=name of the algorithm (with images) used
!!  imgmov=index of algorithm (with images) used
!!  iout=unit number for output
!!  mpi_enreg=MPI-parallelisation information
!!  nimage=number of images stored on current proc
!!  nimage_tot=total number of images (should be dtset%nimage)
!!  prt_all_images=true if all images have to be printed out (ignoring dynimage)
!!  prtvolimg=printing volume for each image
!!           <0 : nothing
!!            0 : only a title
!!            1 : energy, residuals, forces, stresses, velocities, atomic positions
!!            2 : energy, residuals
!!  resimg(nimage) <type(results_img_type)>=results of the ground-state computations
!!                                          for all images treated by current proc
!!
!!
!! OUTPUT
!!  (data written to unit iout)
!!
!! PARENTS
!!      gstateimg
!!
!! CHILDREN
!!      destroy_results_img,gather_results_img,metric,prtxvf,wrtout,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prtimg(dynimage,imagealgo_str,imgmov,iout,mpi_enreg,nimage,nimage_tot,&
&                 prt_all_images,prtvolimg,resimg)

 use defs_basis
 use defs_abitypes
 use m_results_img
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtimg'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_67_common, except_this_one => prtimg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nimage_tot,dynimage(nimage_tot),imgmov,iout,nimage,prtvolimg !vz_d
 logical,intent(in) :: prt_all_images
 character(len=60),intent(in) :: imagealgo_str
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 type(results_img_type),target,intent(inout) :: resimg(nimage)

!Local variables-------------------------------
!scalars
 integer :: ii,prtvel
 logical :: test_img
 real(dp) :: ucvol_img
 character(len=500) :: msg
!arrays
 integer,allocatable :: iatfix_img(:,:)
 real(dp),allocatable :: gmet_img(:,:),gprimd_img(:,:),rmet_img(:,:),xcart_img(:,:)
 type(results_img_type),pointer :: resimg_all(:)

! ****************************************************************

 DBG_ENTER('COLL')

 if (prtvolimg<=0) return
 if (mpi_enreg%me_cell/=0) return

!Gather data
 if (prtvolimg==1.or.prtvolimg==2) then
   test_img=(nimage_tot/=1.and.mpi_enreg%paral_img==1)
   if (test_img) then
     if (mpi_enreg%me==0)  then
       ABI_DATATYPE_ALLOCATE(resimg_all,(nimage_tot))
     end if
     call gather_results_img(mpi_enreg,resimg,resimg_all,master=0,&
&     allgather=.false.,only_one_per_img=.true.)
   else
     resimg_all => resimg
   end if
 end if

!===== First option for the printing volume ===
 if (prtvolimg==1.and.mpi_enreg%me==0) then

   prtvel=0;if (imgmov==0.or.imgmov==9.or.imgmov==10.or.imgmov==13) prtvel=1

   do ii=1,nimage_tot
     if (dynimage(ii)==1.or.prt_all_images) then

!      Title
       write(msg,'(6a,i4,a,i4,2a)') ch10,&
&       '----------------------------------------------------------------------',ch10,&
&       ' ',trim(imagealgo_str),' - CELL # ',ii,'/',nimage_tot,ch10,&
&       '----------------------------------------------------------------------'
       call wrtout(iout,msg,'COLL')

!      Total energy
       write(msg,'(2a,es20.12)') ch10,' Total energy for the cell [Ha]: ',resimg_all(ii)%results_gs%etotal
       call wrtout(iout,msg,'COLL')

!      Residuals of the SCF cycle
       write(msg,'(3a,4(a,es16.8,a))') ch10,&
&       ' Residuals from SCF cycle: ',ch10,&
&       '    Total energy difference        =',resimg_all(ii)%results_gs%deltae,ch10,&
&       '    Maximal forces difference      =',resimg_all(ii)%results_gs%diffor,ch10,&
&       '    Max. residual of wave-functions=',resimg_all(ii)%results_gs%residm,ch10,&
&       '    Density/potential residual (^2)=',resimg_all(ii)%results_gs%res2,ch10
       call wrtout(iout,msg,'COLL')

!      Cell parameters
       ABI_ALLOCATE(rmet_img,(3,3))
       ABI_ALLOCATE(gmet_img,(3,3))
       ABI_ALLOCATE(gprimd_img,(3,3))
       call metric(gmet_img,gprimd_img,iout,rmet_img,resimg_all(ii)%rprim,ucvol_img)
       ABI_DEALLOCATE(rmet_img)
       ABI_DEALLOCATE(gmet_img)
       ABI_DEALLOCATE(gprimd_img)

!      Positions, forces and velocities
       ABI_ALLOCATE(iatfix_img,(3,resimg_all(ii)%natom))
       ABI_ALLOCATE(xcart_img,(3,resimg_all(ii)%natom))
       iatfix_img=0
       call xred2xcart(resimg_all(ii)%natom,resimg_all(ii)%rprim,xcart_img,resimg_all(ii)%xred)
       call prtxvf(resimg_all(ii)%results_gs%fcart,resimg_all(ii)%results_gs%fred,&
&       iatfix_img,iout,resimg_all(ii)%natom,prtvel,&
&       resimg_all(ii)%vel,xcart_img,resimg_all(ii)%xred)
       ABI_DEALLOCATE(iatfix_img)
       ABI_DEALLOCATE(xcart_img)

!      Stress tensor
       write(msg, '(a,es12.4,a)' ) &
&       '-Cartesian components of stress tensor (GPa)         [Pressure=',&
&       -(resimg_all(ii)%results_gs%strten(1)+resimg_all(ii)%results_gs%strten(2) &
&       +resimg_all(ii)%results_gs%strten(3))*HaBohr3_GPa/three,' GPa]'
       call wrtout(iout,msg,'COLL')
       write(msg, '(2(a,1p,e16.8))' ) '- sigma(1 1)=',resimg_all(ii)%results_gs%strten(1)*HaBohr3_GPa,&
&       '  sigma(3 2)=',resimg_all(ii)%results_gs%strten(4)*HaBohr3_GPa
       call wrtout(iout,msg,'COLL')
       write(msg, '(2(a,1p,e16.8))' ) '- sigma(2 2)=',resimg_all(ii)%results_gs%strten(2)*HaBohr3_GPa,&
&       '  sigma(3 1)=',resimg_all(ii)%results_gs%strten(5)*HaBohr3_GPa
       call wrtout(iout,msg,'COLL')
       write(msg, '(2(a,1p,e16.8))' ) '- sigma(3 3)=',resimg_all(ii)%results_gs%strten(3)*HaBohr3_GPa,&
&       '  sigma(2 1)=',resimg_all(ii)%results_gs%strten(6)*HaBohr3_GPa
       call wrtout(iout,msg,'COLL')
     end if
   end do
 end if


!===== 2nd option for the printing volume ===
 if (prtvolimg==2.and.mpi_enreg%me==0) then
   write(msg,'(a,1x,a)') ch10,&
&   'Cell   Total_energy[Ha]     deltae       diffor       residm         res2'
   call wrtout(iout,msg,'COLL')
   do ii=1,nimage_tot
     if (dynimage(ii)==1.or.prt_all_images) then
       write(msg,'(1x,i4,2x,es16.8,4(1x,es13.5))') &
&       ii,resimg_all(ii)%results_gs%etotal,resimg_all(ii)%results_gs%deltae,&
&       resimg_all(ii)%results_gs%diffor,resimg_all(ii)%results_gs%residm,&
&       resimg_all(ii)%results_gs%res2
       call wrtout(iout,msg,'COLL')
     end if
   end do
 end if

!=====
 if (prtvolimg==1.or.prtvolimg==2) then
   if (test_img.and.mpi_enreg%me==0) then
     call destroy_results_img(resimg_all)
     ABI_DATATYPE_DEALLOCATE(resimg_all)
   end if
   nullify(resimg_all)
 end if

 DBG_EXIT('COLL')

end subroutine prtimg
!!***
