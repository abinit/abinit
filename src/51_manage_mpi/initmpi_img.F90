!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_img
!! NAME
!!  initmpi_img
!!
!! FUNCTION
!!  Initializes the mpi informations for parallelism over images of the cell (npimage>1).
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2016 ABINIT group (MT,GG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mpi_enreg= informations about MPI parallelization
!!  option= see below
!!
!! OUTPUT
!!  mpi_enreg%my_nimage= number of images of the cell treated by current proc
!!  ===== if option==1 or option==-1
!!    mpi_enreg%my_imgtab= indexes of images of the cell treated by current proc
!!  ===== if option==2 or option==3 or option==-1
!!    mpi_enreg%comm_cell=Communicator over all processors treating the same image
!!    mpi_enreg%nproc_cell=size of comm_cell
!!    mpi_enreg%me_cell=my rank in comm_cell
!!  ===== if option==3 or option==-1
!!    mpi_enreg%comm_img=Communicator over all images
!!    mpi_enreg%nproc_img=size of comm_img
!!    mpi_enreg%me_img=my rank in comm_img
!!    mpi_enreg%distrb_img(:)=index of processor treating each image (in comm_img communicator)
!!
!! PARENTS
!!      mpi_setup
!!
!! CHILDREN
!!      sort_int
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initmpi_img(dtset,mpi_enreg,option)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_xmpi
 use m_profiling_abi
 use m_sort

 use m_io_tools,  only: flush_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmpi_img'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: option
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
 integer :: imod,irank,iprocmax,iprocmin,jrank
 integer :: ndynimage_eff,nimage_eff,nproc_per_image,nrank
 logical,parameter :: debug=.false.
 integer,allocatable :: ranks(:)
 character(len=500) :: msg

!integer :: group_cell,ierr

! ***********************************************************************

 DBG_ENTER("COLL")
 if (option/=0) then
   mpi_enreg%comm_img=xmpi_comm_self
   mpi_enreg%comm_cell=mpi_enreg%comm_world
 end if

 if (xmpi_paral==1.and.dtset%npimage>1.and.dtset%optdriver==RUNL_GSTATE) then

!  Activate flag for parallelization over images
   mpi_enreg%paral_img=1

   ndynimage_eff=dtset%ndynimage;if (dtset%ntimimage<=1) ndynimage_eff=0

!  Print several warnings
   if (option==0) then
     nimage_eff=max(ndynimage_eff,dtset%nimage-ndynimage_eff)
     if (dtset%npimage>nimage_eff) then
       write(unit=msg,fmt='(3a,i4,a,i4,4a)') &
&       'The number of processors used for the parallelization',ch10,&
&       ' over images (npimage=',dtset%npimage,&
&       ') is greater than the number of dynamic (or static) images (',nimage_eff,') !',ch10,&
&       ' This is unefficient.',ch10
       MSG_WARNING(msg)
     end if
     if (dtset%npimage>mpi_enreg%nproc) then
       write(unit=msg,fmt='(3a,i6,a,i4,4a)') &
&       'The number of processors used for the parallelization',ch10,&
&       ' over images (nproc=',mpi_enreg%nproc,&
&       ') is smaller than npimage in input file (',dtset%npimage,&
&       ')!',ch10,' This is unconsistent.',ch10
       MSG_ERROR(msg)
     end if
     if (mod(nimage_eff,dtset%npimage)/=0) then
       write(unit=msg,fmt='(3a,i4,a,i4,4a)') &
&       'The number of processors used for the parallelization',ch10,&
&       ' over images (npimage=',dtset%npimage,&
&       ') does not divide the number of dynamic images (',nimage_eff,&
&       ') !',ch10,' This is unefficient (charge unbalancing).',ch10
       MSG_WARNING(msg)
     end if
   end if

!  # of images treated by current proc
   nproc_per_image=mpi_enreg%nproc/dtset%npimage
   iprocmax=nproc_per_image*dtset%npimage-1
   if (mpi_enreg%me<=iprocmax) then
     mpi_enreg%my_nimage=(ndynimage_eff/dtset%npimage)+((dtset%nimage-ndynimage_eff)/dtset%npimage)
     imod=mod(ndynimage_eff,dtset%npimage)-1
     if (mpi_enreg%me/nproc_per_image<=imod) mpi_enreg%my_nimage=mpi_enreg%my_nimage+1
     imod=mod((dtset%nimage-ndynimage_eff),dtset%npimage)-1
     if (mpi_enreg%me/nproc_per_image<=imod) mpi_enreg%my_nimage=mpi_enreg%my_nimage+1
   else
     mpi_enreg%my_nimage=0
   end if
   if (option==1.or.option==-1) then
!    Indexes of images treated by current proc
     if (mpi_enreg%me<=iprocmax) then
       ABI_ALLOCATE(mpi_enreg%my_imgtab,(mpi_enreg%my_nimage))
       nrank=0
       imod=mpi_enreg%me/nproc_per_image+1;imod=mod(imod,dtset%npimage)
!      Dynamic images
       irank=0
       do jrank=1,dtset%nimage
         if (dtset%dynimage(jrank)/=0.and.dtset%ntimimage>1) then
           irank=irank+1
           if (mod(irank,dtset%npimage)==imod) then
             nrank=nrank+1
             mpi_enreg%my_imgtab(nrank)=jrank
           end if
         end if
       end do
!      Static images
       irank=0
       do jrank=1,dtset%nimage
         if (dtset%dynimage(jrank)==0.or.dtset%ntimimage<=1) then
           irank=irank+1
           if (mod(irank,dtset%npimage)==imod) then
             nrank=nrank+1
             mpi_enreg%my_imgtab(nrank)=jrank
           end if
         end if
       end do
       if (nrank/=mpi_enreg%my_nimage) then
         MSG_BUG('Error on nrank !')
       end if
!      Sort images by increasing index (this step is MANDATORY !!)
       ABI_ALLOCATE(ranks,(nrank))
       call sort_int(nrank,mpi_enreg%my_imgtab,ranks)
       ABI_DEALLOCATE(ranks)
     else
       ABI_ALLOCATE(mpi_enreg%my_imgtab,(0))
     end if
   end if
   if (option==2.or.option==3.or.option==-1) then
!    Communicator over one image
     if (mpi_enreg%me<=iprocmax) then
       ABI_ALLOCATE(ranks,(nproc_per_image))
       iprocmin=(mpi_enreg%me/nproc_per_image)*nproc_per_image
       ranks=(/((iprocmin+irank-1),irank=1,nproc_per_image)/)
       mpi_enreg%comm_cell=xmpi_subcomm(mpi_enreg%comm_world,nproc_per_image,ranks)
       ABI_DEALLOCATE(ranks)
       mpi_enreg%me_cell=xmpi_comm_rank(mpi_enreg%comm_cell)
       mpi_enreg%nproc_cell=nproc_per_image
       if (mpi_enreg%me_cell==0.and.mod(mpi_enreg%me,nproc_per_image)/=0) then
         MSG_BUG('Error on me_cell !')
       end if
     else
       mpi_enreg%comm_img=xmpi_comm_null
       mpi_enreg%nproc_cell=0
       mpi_enreg%me_cell=-1
     end if
   end if
   if (option==3.or.option==-1) then
!    Communicator over all images
     if (mpi_enreg%me<=iprocmax) then
       ABI_ALLOCATE(ranks,(dtset%npimage))
       iprocmin=mod(mpi_enreg%me,nproc_per_image)
       ranks=(/((iprocmin+(irank-1)*nproc_per_image),irank=1,dtset%npimage)/)
       mpi_enreg%comm_img=xmpi_subcomm(mpi_enreg%comm_world,dtset%npimage,ranks)
       ABI_DEALLOCATE(ranks)
       mpi_enreg%me_img=xmpi_comm_rank(mpi_enreg%comm_img)
       mpi_enreg%nproc_img=dtset%npimage
       if (iprocmin==0.and.mpi_enreg%me_img==0.and.mpi_enreg%me/=0) then
         MSG_BUG('Error on me_img!')
       end if
       ABI_ALLOCATE(mpi_enreg%distrb_img,(dtset%nimage))
!      Dynamic images
       nrank=0
       do irank=1,dtset%nimage
         if (dtset%dynimage(irank)/=0.and.dtset%ntimimage>1) then
           nrank=nrank+1
           mpi_enreg%distrb_img(irank)=mod(nrank,dtset%npimage)-1
           if (mpi_enreg%distrb_img(irank)==-1) mpi_enreg%distrb_img(irank)=dtset%npimage-1
         end if
       end do
!      Static images
       nrank=0
       do irank=1,dtset%nimage
         if (dtset%dynimage(irank)==0.or.dtset%ntimimage<=1) then
           nrank=nrank+1
           mpi_enreg%distrb_img(irank)=mod(nrank,dtset%npimage)-1
           if (mpi_enreg%distrb_img(irank)==-1) mpi_enreg%distrb_img(irank)=dtset%npimage-1
         end if
       end do
     else
       mpi_enreg%comm_img=xmpi_comm_null
       mpi_enreg%nproc_img=0
       mpi_enreg%me_img=-1
       ABI_ALLOCATE(mpi_enreg%distrb_img,(0))
     end if
   end if

!  if (debug) then
!  write(200+mpi_enreg%me,*) "=================================="
!  write(200+mpi_enreg%me,*) "DEBUGGING STATMENTS IN INITMPI_IMG"
!  write(200+mpi_enreg%me,*) "=================================="
!  write(200+mpi_enreg%me,*) "option         =",option
!  write(200+mpi_enreg%me,*) "MPI_UNDEFINED  =",MPI_UNDEFINED
!  write(200+mpi_enreg%me,*) "MPI_IDENT      =",MPI_IDENT
!  write(200+mpi_enreg%me,*) "MPI_CONGRUENT  =",MPI_CONGRUENT
!  write(200+mpi_enreg%me,*) "MPI_SIMILAR    =",MPI_SIMILAR
!  write(200+mpi_enreg%me,*) "MPI_UNEQUAL    =",MPI_UNEQUAL
!  write(200+mpi_enreg%me,*) "null_comm      =",MPI_COMM_NULL
!  write(200+mpi_enreg%me,*) "self_comm      =",xmpi_comm_self
!  write(200+mpi_enreg%me,*) "world_comm     =",mpi_enreg%comm_world
!  write(200+mpi_enreg%me,*) "empty_group    =",MPI_GROUP_EMPTY
!  write(200+mpi_enreg%me,*) "nimage         =",mpi_enreg%my_nimage
!  write(200+mpi_enreg%me,*) "nproc_per_image=",nproc_per_image
!  call MPI_COMM_SIZE(mpi_enreg%comm_world,irank,ierr)
!  write(200+mpi_enreg%me,*) "Size of world_comm    =",irank
!  call MPI_COMM_RANK(mpi_enreg%comm_world,irank,ierr)
!  write(200+mpi_enreg%me,*) "My rank in world_comm =",irank
!  if (option==1.or.option==-1) then
!  write(200+mpi_enreg%me,*) "index_img=",mpi_enreg%my_imgtab(:)
!  end if
!  if (option==2.or.option==3.or.option==-1) then
!  write(200+mpi_enreg%me,*) "nproc_cell  =",mpi_enreg%nproc_cell
!  write(200+mpi_enreg%me,*) "me_cell     =",mpi_enreg%me_cell
!  call xmpi_comm_group(mpi_enreg%comm_cell,group_cell,ierr)
!  write(200+mpi_enreg%me,*) "group_cell  =",group_cell
!  write(200+mpi_enreg%me,*) "comm_cell   =",mpi_enreg%comm_cell
!  if (group_cell/=MPI_GROUP_EMPTY) then
!  call MPI_GROUP_SIZE(group_cell,irank,ierr)
!  write(200+mpi_enreg%me,*) "Size of group_cell   =",irank
!  call MPI_GROUP_RANK(group_cell,irank,ierr)
!  write(200+mpi_enreg%me,*) "My rank in group_cell=",irank
!  else
!  write(200+mpi_enreg%me,*) "Size of group_cell   =",0
!  write(200+mpi_enreg%me,*) "My rank in group_cell=",-1
!  end if
!  if (mpi_enreg%comm_cell/=MPI_COMM_NULL) then
!  call MPI_COMM_SIZE(mpi_enreg%comm_cell,irank,ierr)
!  write(200+mpi_enreg%me,*) "Size of comm_cell   =",irank
!  call MPI_COMM_RANK(mpi_enreg%comm_cell,irank,ierr)
!  write(200+mpi_enreg%me,*) "My rank in comm_cell=",irank
!  call MPI_COMM_COMPARE(mpi_enreg%comm_world,mpi_enreg%comm_cell,irank,ierr)
!  write(200+mpi_enreg%me,*) "Comparison world_comm/comm_cell=",irank
!  call MPI_COMM_COMPARE(xmpi_comm_self,mpi_enreg%comm_cell,irank,ierr)
!  write(200+mpi_enreg%me,*) "Comparison self_comm/comm_cell =",irank
!  else
!  write(200+mpi_enreg%me,*) "Size of comm_cell   =",0
!  write(200+mpi_enreg%me,*) "My rank in comm_cell=",-1
!  write(200+mpi_enreg%me,*) "Comparison world_comm/comm_cell=",MPI_UNEQUAL
!  write(200+mpi_enreg%me,*) "Comparison self_comm/comm_cell =",MPI_UNEQUAL
!  end if
!  end if
!  if (option==3.or.option==-1) then
!  write(200+mpi_enreg%me,*) "nproc_img  =",mpi_enreg%nproc_img
!  write(200+mpi_enreg%me,*) "me_img     =",mpi_enreg%me_img
!  write(200+mpi_enreg%me,*) "img_comm   =",mpi_enreg%comm_img
!  if (mpi_enreg%comm_img/=MPI_COMM_NULL) then
!  call MPI_COMM_SIZE(mpi_enreg%comm_img,irank,ierr)
!  write(200+mpi_enreg%me,*) "Size of img_comm   =",irank
!  call MPI_COMM_RANK(mpi_enreg%comm_img,irank,ierr)
!  write(200+mpi_enreg%me,*) "My rank in img_comm=",irank
!  call MPI_COMM_COMPARE(mpi_enreg%comm_world,mpi_enreg%comm_img,irank,ierr)
!  write(200+mpi_enreg%me,*) "Comparison world_comm/img_comm=",irank
!  call MPI_COMM_COMPARE(xmpi_comm_self,mpi_enreg%comm_img,irank,ierr)
!  write(200+mpi_enreg%me,*) "Comparison self_comm/img_comm =",irank
!  else
!  write(200+mpi_enreg%me,*) "Size of img_comm   =",0
!  write(200+mpi_enreg%me,*) "My rank in img_comm=",-1
!  write(200+mpi_enreg%me,*) "Comparison world_comm/img_comm=",MPI_UNEQUAL
!  write(200+mpi_enreg%me,*) "Comparison self_comm/img_comm =",MPI_UNEQUAL
!  end if
!  write(200+mpi_enreg%me,*) "distrb_img=",mpi_enreg%distrb_img(:)
!  end if
!  write(200+mpi_enreg%me,*)
!  call flush_unit(200+mpi_enreg%me)
!  if (option==-1) stop
!  end if

 else

!  Do not activate flag for parallelization over images
   mpi_enreg%paral_img=0
!  # of images treated by current proc
   if (dtset%optdriver==RUNL_GSTATE) then
     mpi_enreg%my_nimage=dtset%nimage
   else
     mpi_enreg%my_nimage=1
   end if
!  Indexes of images treated by current proc
   if (option==1.or.option==-1) then
     ABI_ALLOCATE(mpi_enreg%my_imgtab,(mpi_enreg%my_nimage))
     mpi_enreg%my_imgtab=(/(irank,irank=1,mpi_enreg%my_nimage)/)
   end if
!  Communicator over all images
   if (option==2.or.option==3.or.option==-1) then
!    Communicator for one image
     mpi_enreg%nproc_cell=mpi_enreg%nproc
     mpi_enreg%me_cell=mpi_enreg%me
   end if
   if (option==3.or.option==-1) then
!    Communicator over all images
     mpi_enreg%nproc_img=1
     mpi_enreg%comm_img=xmpi_comm_self
     mpi_enreg%me_img=0
     ABI_ALLOCATE(mpi_enreg%distrb_img,(dtset%nimage))
     mpi_enreg%distrb_img(:)=0
   end if
 end if

 DBG_EXIT("COLL")

end subroutine initmpi_img
!!***
