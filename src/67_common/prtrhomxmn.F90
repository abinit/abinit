!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtrhomxmn
!! NAME
!! prtrhomxmn
!!
!! FUNCTION
!! If option==1, compute the maximum and minimum of the density (and spin-polarization
!! if nspden==2), and print it.
!! If option==2, also compute and print the second maximum or minimum
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  iout=unit for output file
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  option, see above
!!  optrhor=option for rhor (If optrhor==0, rhor is expected to be the electron density)
!!                          (If optrhor==1, rhor is expected to be the kinetic energy density (taur))
!!                          (If optrhor==2, rhor is expected to be the gradient of the electron density (grhor))
!!                          (If optrhor==3, rhor is expected to be the laplacian of the electron density (lrhor))
!!                          (If optrhor==4, rhor is expected to be the ELF (elfr))
!!  rhor(nfft,nspden)=electron density (electrons/bohr^3)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  The tolerance tol12 aims at giving a machine-independent ordering.
!!  (this trick is used in bonds.f, listkk.f, prtrhomxmn.f and rsiaf9.f)
!!
!! PARENTS
!!      afterscfloop,bethe_salpeter,clnup1,mkrho,screening,sigma,vtorho
!!
!! CHILDREN
!!      wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prtrhomxmn(iout,mpi_enreg,nfft,ngfft,nspden,option,rhor,optrhor,ucvol)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtrhomxmn'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,nfft,nspden,option
 integer,intent(in),optional :: optrhor
 real(dp),intent(in),optional :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,ierr,ifft,ii,iisign,iitems,index1,ioptrhor
 integer :: index2,indsign,iproc,istart,me,n1,n2,n3,nitems
 integer :: nfft_,nfftot,nproc,spaceComm
 real(dp) :: temp,value1,value2
 character(len=500) :: message,txt1_in_mssg,txt2_in_mssg,txt3_in_mssg
 logical :: reduce=.false.
!arrays
 integer,allocatable :: iindex(:,:,:),index_fft(:,:,:,:)
 real(dp) :: rhomn1(4),rhomn2(4),rhomx1(4),rhomx2(4),ri_rhomn1(3,4)
 real(dp) :: ri_rhomn2(3,4),ri_rhomx1(3,4),ri_rhomx2(3,4),ri_zetmn1(3,2)
 real(dp) :: ri_zetmn2(3,2),ri_zetmx1(3,2),ri_zetmx2(3,2),zetmn1(2)
 real(dp) :: zetmn2(2),zetmx1(2),zetmx2(2)
 real(dp),allocatable :: array(:),coord(:,:,:,:),value(:,:,:),integrated(:)
 real(dp),allocatable :: value_fft(:,:,:)

! *************************************************************************

 if(.not.(present(optrhor))) then
   ioptrhor=0
 else
   ioptrhor=optrhor
 end if

 if(option/=1 .and. option/=2)then
   write(message, '(a,i0)' )' Option must be 1 or 2, while it is ',option
   MSG_BUG(message)
 end if

 if (mpi_enreg%nproc_wvl>1) then
!  nfft is always the potential size (in GGA, the density has buffers).
   nfft_ = ngfft(1) * ngfft(2) * mpi_enreg%nscatterarr(mpi_enreg%me_wvl, 2)
   n1 = ngfft(1)
   n2 = ngfft(2)
   n3 = sum(mpi_enreg%nscatterarr(:, 2))
   istart = mpi_enreg%nscatterarr(mpi_enreg%me_wvl, 4)
 else
   nfft_ = nfft
   n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
   istart = 0
 end if

!--------------------------------------------------------------------------
!One has to determine the maximum and minimum (etc...) values
!over all space, and then output it, as well as to identify
!the point at which it occurs ...
!This will require a bit of data exchange, and correct indirect indexing ...

!For the local processor, find different items :
!maximum and minimum total electron density and locations
!and also spin-polarisation and magnetization
!also keep the second maximal or minimal value
 if(nspden==1)nitems=1   ! Simply the total density
 if(nspden==2)nitems=5   ! Total density, spin up, spin down, magnetization, zeta
 if(nspden==4)nitems=6   ! Total density, x, y, z, magnetization, zeta

 ABI_ALLOCATE(value,(2,2,nitems))
 ABI_ALLOCATE(iindex,(2,2,nitems))
 ABI_ALLOCATE(array,(nfft))
 ABI_ALLOCATE(integrated,(nitems))

 do iitems=1,nitems

!  Copy the correct values into the array
!  First set of items : the density, for each spin component
   if(iitems<=nspden)then
     array(:)=rhor(:,iitems)
   end if
!  Case nspden==2, some computation to be done
   if(nspden==2)then
     if(iitems==3)then ! Spin down
       array(:)=rhor(:,1)-rhor(:,2)
     else if(iitems==4)then  ! Magnetization
       array(:)=2*rhor(:,2)-rhor(:,1)
     else if(iitems==5)then  ! zeta = relative magnetization
       ! Avoid 0/0: the limit of (x - y) / (x+ y) depends on the direction.
       array(:)=zero
       where (abs(rhor(:,1)) > tol12) array(:)=(2*rhor(:,2)-rhor(:,1))/rhor(:,1)
     end if
!    Case nspden==4, some other computation to be done
   else if(nspden==4)then
     if(iitems==5)then ! Magnetization
       array(:)=sqrt(rhor(:,2)**2+rhor(:,3)**2+rhor(:,4)**2)
     else if(iitems==6)then ! zeta = relative magnetization
       array(:)=(sqrt(rhor(:,2)**2+rhor(:,3)**2+rhor(:,4)**2))/rhor(:,1)
     end if
   end if

!  Zero all the absolute values that are lower than tol8, for portability reasons.
   do ifft = 1, nfft_
     if(abs(array(ifft))<tol8)array(ifft)=zero
   end do

!  DEBUG
!  write(std_out,*) ' iitems,array(1:2)=',iitems,array(1:2)
!  ENDDEBUG

   do indsign=1,2 ! Find alternatively the maximum and the minimum
     iisign=3-2*indsign

     if (nfft_ > 1) then
!      Initialize the two first values
       value1=array(istart + 1) ; value2=array(istart + 2)
       index1=1 ; index2=2

!      Ordering, if needed
       if( iisign*(value2+tol12) > iisign*(value1)) then
         temp=value2 ; value2=value1 ; value1=temp
         index1=2 ; index2=1
       end if

!      Integration, if relevant
       if(present(ucvol).and. indsign==1)then
         integrated(iitems) = array(istart + 1)+array(istart + 2)
       end if
     else
       value1 = zero; value2 = zero
       index1 = 0;    index2 = 0
     end if

!    DEBUG
!    write(std_out,*) ' value1,value2,index1,index2=',value1,value2,index1,index2
!    ENDDEBUG

!    Loop over all points
     do ifft = 3, nfft_

       temp=array(istart + ifft)
       if(present(ucvol).and. indsign==1)integrated(iitems) = integrated(iitems)+temp
!      Compares it to the second value
       if( iisign*(temp+tol12) > iisign*value2 ) then
!        Compare it to the first value
         if( iisign*(temp+tol12) > iisign*value1 ) then
           value2=value1 ; index2=index1
           value1=temp   ; index1=ifft
         else
           value2=temp   ; index2=ifft
         end if
       end if

     end do ! ifft

     value(1,indsign,iitems)=value1
     value(2,indsign,iitems)=value2
     iindex(1,indsign,iitems)=index1
     iindex(2,indsign,iitems)=index2

!    DEBUG
!    write(std_out,*) ' it,v1,i1=',iitems, value1,index1
!    write(std_out,*) ' it,v2,i2=',iitems, value2,index2
!    ENDDEBUG

   end do ! indsign

   if(present(ucvol))then
     nfftot=ngfft(1) * ngfft(2) * ngfft(3)
     integrated(iitems)=integrated(iitems)*ucvol/nfftot
   end if

!  Integrate the array
!  integrated(iitems)=zero
!  do ifft=1,nfft_
!  integrated(iitems) = integrated(iitems) + array(istart + ifft)
!  enddo
!  if(present(ucvol))integrated(iitems) = integrated(iitems)*ucvol/nfft_
!  write(std_err,*)present(ucvol)
!  if(present(ucvol))then
!  write(std_err,*)ucvol
!  endif

 end do ! iitems

 ABI_DEALLOCATE(array)

!-------------------------------------------------------------------
!Enter section for FFT parallel case
!if(mpi_enreg%paral_kgb>1) spaceComm=mpi_enreg%comm_fft; reduce=.true.
 spaceComm=mpi_enreg%comm_fft; reduce=.false.
 if(mpi_enreg%nproc_fft>1) then
   spaceComm=mpi_enreg%comm_fft; reduce=.true.
 else if(mpi_enreg%nproc_wvl>1) then
   spaceComm=mpi_enreg%comm_wvl; reduce=.true.
 end if
 nproc=xmpi_comm_size(spaceComm)
 me=xmpi_comm_rank(spaceComm)

 if (reduce) then

!  Communicate all data to all processors with only two global communications
   ABI_ALLOCATE(value_fft,(5,nitems,nproc))
   ABI_ALLOCATE(index_fft,(2,2,nitems,nproc))
   value_fft(:,:,:)=zero
   index_fft(:,:,:,:)=zero
   value_fft(1,:,me + 1)=value(1,1,:)
   value_fft(2,:,me + 1)=value(2,1,:)
   value_fft(3,:,me + 1)=value(1,2,:)
   value_fft(4,:,me + 1)=value(2,2,:)
   if(present(ucvol))value_fft(5,:,me + 1)=integrated(:)
   index_fft(:,:,:,me + 1)=iindex(:,:,:)
   call xmpi_sum(value_fft,spaceComm,ierr)
   call xmpi_sum(index_fft,spaceComm,ierr)

!  Determine the global optimum and second optimum for each item
!  Also, the integrated quantities, if relevant.
   do iitems=1,nitems

     if(present(ucvol))integrated(iitems)=sum(value_fft(5,iitems,1:nproc))

     do indsign=1,2 ! Find alternatively the maximum and the minimum
       iisign=3-2*indsign

!      Initialisation
       value1=value_fft(2*indsign-1,iitems,1)
       value2=value_fft(2*indsign  ,iitems,1)
       index1=index_fft(1,indsign,iitems,1)
       index2=index_fft(2,indsign,iitems,1)

!      Loop
       do iproc=1, nproc, 1
         do ii=1,2
           if(iproc>1 .or. ii==2)then

             temp=value_fft(ii+2*(indsign-1),iitems,iproc)
!            Compares it to the second value
             if( iisign*(temp+tol12) > iisign*value2 ) then
!              Compare it to the first value
               if( iisign*(temp+tol12) > iisign*value1 ) then
                 value2=value1 ; index2=index1
                 value1=temp   ; index1=index_fft(ii,indsign,iitems,iproc)
               else
                 value2=temp   ; index2=index_fft(ii,indsign,iitems,iproc)
               end if
             end if

           end if ! if(iproc>1 .or. ii==2)
         end do ! ii
       end do ! iproc

       value(1,indsign,iitems)=value1
       value(2,indsign,iitems)=value2
       iindex(1,indsign,iitems)=index1
       iindex(2,indsign,iitems)=index2

     end do ! iisign
   end do ! iitems

   ABI_DEALLOCATE(value_fft)
   ABI_DEALLOCATE(index_fft)

 end if !if(reduce)

!-------------------------------------------------------------------

!Determines the reduced coordinates of the min and max for each item
 ABI_ALLOCATE(coord,(3,2,2,nitems))
 do iitems=1,nitems
   do indsign=1,2
     do ii=1,2
       index1=iindex(ii,indsign,iitems)
       i3=(index1-1)/n1/n2
       i2=(index1-1-i3*n1*n2)/n1
       i1=index1-1-i3*n1*n2-i2*n1
       coord(1,ii,indsign,iitems)=dble(i1)/dble(n1)+tol12
       coord(2,ii,indsign,iitems)=dble(i2)/dble(n2)+tol12
       coord(3,ii,indsign,iitems)=dble(i3)/dble(n3)+tol12
!      DEBUG
!      write(std_out,*)' ii,indsign,iitems,coord(1:3)=',ii,indsign,iitems,coord(:,ii,indsign,iitems)
!      write(std_out,*)' value ', value(ii, indsign, iitems)
!      ENDDEBUG
     end do
   end do
 end do

!-------------------------------------------------------------------------
!Output
 if (mpi_enreg%paral_kgb==0.or.mpi_enreg%me_fft==0) then
   if(.true.)then
     do iitems=1,nitems

       if(ioptrhor==4 .and. iitems>2)exit

       select case (ioptrhor)
       case(0)

         if(iitems==1) write(message,'(a)')' Total charge density [el/Bohr^3]'
         if(nspden==2)then
           if(iitems==2) write(message,'(a)')' Spin up density      [el/Bohr^3]'
           if(iitems==3) write(message,'(a)')' Spin down density    [el/Bohr^3]'
           if(iitems==4) write(message,'(a)')' Magnetization (spin up - spin down) [el/Bohr^3]'
           if(iitems==5) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         else if(nspden==4)then
           if(iitems==2) write(message,'(a)')' x component of magnetization [el/Bohr^3]'
           if(iitems==3) write(message,'(a)')' y component of magnetization [el/Bohr^3]'
           if(iitems==4) write(message,'(a)')' z component of magnetization [el/Bohr^3]'
           if(iitems==5) write(message,'(a)')' Magnetization (absolute value) [el/Bohr^3]'
           if(iitems==6) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         end if

       case(1)

         if(iitems==1) write(message,'(a)')' Total kinetic energy density [Ha/Bohr^3]'
         if(nspden==2)then
           if(iitems==2) write(message,'(a)')' Spin up density      [Ha/Bohr^3]'
           if(iitems==3) write(message,'(a)')' Spin down density    [Ha/Bohr^3]'
           if(iitems==4) write(message,'(a)')' Magnetization (spin up - spin down) [Ha/Bohr^3]'
           if(iitems==5) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         else if(nspden==4)then
           if(iitems==2) write(message,'(a)')' x component of magnetization [Ha/Bohr^3]'
           if(iitems==3) write(message,'(a)')' y component of magnetization [Ha/Bohr^3]'
           if(iitems==4) write(message,'(a)')' z component of magnetization [Ha/Bohr^3]'
           if(iitems==5) write(message,'(a)')' Magnetization (absolute value) [Ha/Bohr^3]'
           if(iitems==6) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         end if

       case(2)

         if(iitems==1) write(message,'(a)')' Gradient of the electronic density [el/Bohr^4]'
         if(nspden==2)then
           if(iitems==2) write(message,'(a)')' Spin up density      [el/Bohr^4]'
           if(iitems==3) write(message,'(a)')' Spin down density    [el/Bohr^4]'
           if(iitems==4) write(message,'(a)')' Magnetization (spin up - spin down) [el/Bohr^4]'
           if(iitems==5) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         else if(nspden==4)then
           if(iitems==2) write(message,'(a)')' x component of magnetization [el/Bohr^4]'
           if(iitems==3) write(message,'(a)')' y component of magnetization [el/Bohr^4]'
           if(iitems==4) write(message,'(a)')' z component of magnetization [el/Bohr^4]'
           if(iitems==5) write(message,'(a)')' Magnetization (absolute value) [el/Bohr^4]'
           if(iitems==6) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         end if

       case(3)

         if(iitems==1) write(message,'(a)')' Laplacian of the electronic density [el/Bohr^5]'
         if(nspden==2)then
           if(iitems==2) write(message,'(a)')' Spin up density      [el/Bohr^5]'
           if(iitems==3) write(message,'(a)')' Spin down density    [el/Bohr^5]'
           if(iitems==4) write(message,'(a)')' Magnetization (spin up - spin down) [el/Bohr^5]'
           if(iitems==5) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         else if(nspden==4)then
           if(iitems==2) write(message,'(a)')' x component of magnetization [el/Bohr^5]'
           if(iitems==3) write(message,'(a)')' y component of magnetization [el/Bohr^5]'
           if(iitems==4) write(message,'(a)')' z component of magnetization [el/Bohr^5]'
           if(iitems==5) write(message,'(a)')' Magnetization (absolute value) [el/Bohr^5]'
           if(iitems==6) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         end if

       case(4)

         if(iitems==1) write(message,'(a)')' Electron Localization Function (ELF) [min:0;max:1]'
         if(nspden==2)then
           if(iitems==2) write(message,'(a)')' Spin up ELF      [min:0;max:1]'
!            if(iitems==3) write(message,'(a)')' Spin down ELF    [min:0;max:1]'
!            if(iitems==4) write(message,'(a)')' Magnetization (spin up - spin down) [el/Bohr^4]'
!            if(iitems==5) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         else if(nspden==4)then
!            if(iitems==2) write(message,'(a)')' x component of magnetization [el/Bohr^4]'
!            if(iitems==3) write(message,'(a)')' y component of magnetization [el/Bohr^4]'
!            if(iitems==4) write(message,'(a)')' z component of magnetization [el/Bohr^4]'
!            if(iitems==5) write(message,'(a)')' Magnetization (spin up - spin down) [el/Bohr^4]'
!            if(iitems==6) write(message,'(a)')' Relative magnetization (=zeta, between -1 and 1)   '
         end if


       end select

       call wrtout(iout,message,'COLL')

       write(message,'(a,es13.4,a,3f10.4)') '      Maximum= ',&
&       value(1,1,iitems),'  at reduced coord.',coord(:,1,1,iitems)
       call wrtout(iout,message,'COLL')
       if(option==2)then
         write(message,'(a,es13.4,a,3f10.4)')' Next maximum= ',&
&         value(2,1,iitems),'  at reduced coord.',coord(:,2,1,iitems)
         call wrtout(iout,message,'COLL')
       end if
       write(message,'(a,es13.4,a,3f10.4)') '      Minimum= ',&
&       value(1,2,iitems),'  at reduced coord.',coord(:,1,2,iitems)
       call wrtout(iout,message,'COLL')
       if(option==2)then
         write(message,'(a,es13.4,a,3f10.4)')' Next minimum= ',&
&         value(2,2,iitems),'  at reduced coord.',coord(:,2,2,iitems)
         call wrtout(iout,message,'COLL')
       end if
       if(present(ucvol))then
         if(.not.(nspden==2.and.iitems==5) .and. .not.(nspden==4.and.iitems==6))then
           if(abs(integrated(iitems))<tol10)integrated(iitems)=zero
           write(message,'(a,es13.4)')'   Integrated= ',integrated(iitems)
           call wrtout(iout,message,'COLL')
         end if
       end if

     end do ! iitems
   end if

   if(.false.)then

     select case(optrhor)
     case(0)
       write(txt1_in_mssg, '(a)')" Min el dens="
       write(txt2_in_mssg, '(a)')" el/bohr^3 at reduced coord."
       write(txt3_in_mssg, '(a)')" Max el dens="
     case(1)
       write(txt1_in_mssg, '(a)')" Min kin energy dens="
       write(txt2_in_mssg, '(a)')" bohr^(-5) at reduced coord."
       write(txt3_in_mssg, '(a)')" Max kin energy dens="
     end select

     write(message, '(a,a,1p,e12.4,a,0p,3f8.4)' ) ch10,&
&     trim(txt1_in_mssg),value(1,2,1),&
&     trim(txt2_in_mssg),coord(:,1,2,1)
     call wrtout(iout,message,'COLL')
     if(option==2)then
       write(message, '(a,1p,e12.4,a,0p,3f8.4)' ) &
&       ',   next min=',value(2,2,1),&
&       trim(txt2_in_mssg),coord(:,2,2,1)
       call wrtout(iout,message,'COLL')
     end if
     write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&     trim(txt3_in_mssg),value(1,1,1),&
&     trim(txt2_in_mssg),coord(:,1,1,1)
     call wrtout(iout,message,'COLL')
     if(option==2)then
       write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&       ',   next max=',value(2,1,1),&
&       trim(txt2_in_mssg),coord(:,2,1,1)
       call wrtout(iout,message,'COLL')
     end if

     if(nspden>=2)then
       write(message, '(a,a,1p,e12.4,a,0p,3f8.4)' ) ch10,&
&       ',Min spin pol zeta=',value(1,2,4+nspden/2),&
&       ' at reduced coord.',coord(:,1,2,4+nspden/2)
       call wrtout(iout,message,'COLL')
       if(option==2)then
         write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&         ',         next min=',value(2,2,4+nspden/2),&
&         ' at reduced coord.',coord(:,2,2,4+nspden/2)
         call wrtout(iout,message,'COLL')
       end if
       write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&       ',Max spin pol zeta=',value(1,1,4+nspden/2),&
&       ' at reduced coord.',coord(:,1,1,4+nspden/2)
       call wrtout(iout,message,'COLL')
       if(option==2)then
         write(message, '(a,1p,e12.4,a,0p,3f8.4)' )&
&         ',         next max=',value(2,1,4+nspden/2),&
&         ' at reduced coord.',coord(:,2,1,4+nspden/2)
         call wrtout(iout,message,'COLL')
       end if
     end if ! nspden

   end if ! second section always true

   if(nspden==2 .and. .false.)then
     write(message,'(a)')&
&     '                               Position in reduced coord.       (  x         y         z )'
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Total  el-den) : [el/Bohr^3]',&
&     rhomn1(1),'  at',ri_rhomn1(1,1),ri_rhomn1(2,1),ri_rhomn1(3,1)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Spin-up   den) : [el/Bohr^3]',&
&     rhomn1(2),'  at',ri_rhomn1(1,2),ri_rhomn1(2,2),ri_rhomn1(3,2)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Spin-down den) : [el/Bohr^3]',&
&     zetmn1(1),'  at',ri_zetmn1(1,1),ri_zetmn1(2,1),ri_zetmn1(3,1)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Spin pol zeta) :   [m/|m|]  ',&
&     zetmn1(2),'  at',ri_zetmn1(1,2),ri_zetmn1(2,2),ri_zetmn1(3,2)
     call wrtout(iout,message,'COLL')
     if(option==2)then
       write(message,'(a,es13.4,a,3f10.4)')' Next minimum (Total  el-den) : [el/Bohr^3]',&
&       rhomn2(1),'  at',ri_rhomn2(1,1),ri_rhomn2(2,1),ri_rhomn2(3,1)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next minimum (Spin-up   den) : [el/Bohr^3]',&
&       rhomn2(2),'  at',ri_rhomn2(1,2),ri_rhomn2(2,2),ri_rhomn2(3,2)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next minimum (Spin-down den) : [el/Bohr^3]',&
&       zetmn2(1),'  at',ri_zetmn2(1,1),ri_zetmn2(2,1),ri_zetmn2(3,1)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next minimum (Spin pol zeta) :   [m/|m|]  ',&
&       zetmn2(2),'  at',ri_zetmn2(1,2),ri_zetmn2(2,2),ri_zetmn2(3,2)
       call wrtout(iout,message,'COLL')
     end if
     write(message,*)' '
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Total  el-den) : [el/Bohr^3]',&
&     rhomx1(1),'  at',ri_rhomx1(1,1),ri_rhomx1(2,1),ri_rhomx1(3,1)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Spin-up   den) : [el/Bohr^3]',&
&     rhomx1(2),'  at',ri_rhomx1(1,2),ri_rhomx1(2,2),ri_rhomx1(3,2)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Spin-down den) : [el/Bohr^3]',&
&     zetmx1(1),'  at',ri_zetmx1(1,1),ri_zetmx1(2,1),ri_zetmx1(3,1)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Spin pol zeta) :   [m/|m|]  ',&
&     zetmx1(2),'  at',ri_zetmx1(1,2),ri_zetmx1(2,2),ri_zetmx1(3,2)
     call wrtout(iout,message,'COLL')
     if(option==2)then
       write(message,'(a,es13.4,a,3f10.4)')' Next maximum (Total  el-den) : [el/Bohr^3]',&
&       rhomx2(1),'  at',ri_rhomx2(1,1),ri_rhomx2(2,1),ri_rhomx2(3,1)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next maximum (Spin-up   den) : [el/Bohr^3]',&
&       rhomx2(2),'  at',ri_rhomx2(1,2),ri_rhomx2(2,2),ri_rhomx2(3,2)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next maximum (Spin-down den) : [el/Bohr^3]',&
&       zetmx2(1),'  at',ri_zetmx2(1,1),ri_zetmx2(2,1),ri_zetmx2(3,1)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next maximum (Spin pol zeta) :   [m/|m|]  ',&
&       zetmx2(2),'  at',ri_zetmx2(1,2),ri_zetmx2(2,2),ri_zetmx2(3,2)
       call wrtout(iout,message,'COLL')
     end if
   end if

   if(nspden==4 .and. .false.)then
     write(message,'(a)')&
&     '                               Position in reduced coord.       (  x         y         z )'
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Total  el-den) : [el/Bohr^3]',&
&     rhomn1(1),'  at',ri_rhomn1(1,1),ri_rhomn1(2,1),ri_rhomn1(3,1)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Magnetizat.-x) :   [m/|m|]  ',&
&     rhomn1(2),'  at',ri_rhomn1(1,2),ri_rhomn1(2,2),ri_rhomn1(3,2)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Magnetizat.-y) :   [m/|m|]  ',&
&     rhomn1(3),'  at',ri_rhomn1(1,3),ri_rhomn1(2,3),ri_rhomn1(3,3)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Magnetizat.-z) :   [m/|m|]  ',&
&     rhomn1(4),'  at',ri_rhomn1(1,4),ri_rhomn1(2,4),ri_rhomn1(3,4)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Minimum (Spin pol zeta) :   [m/|m|]  ',&
&     zetmn1(1),'  at',ri_zetmn1(1,1),ri_zetmn1(2,1),ri_zetmn1(3,1)
     call wrtout(iout,message,'COLL')
     if(option==2)then
       write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Total  el-den) : [el/Bohr^3]',&
&       rhomn2(1),'  at',ri_rhomn2(1,1),ri_rhomn2(2,1),ri_rhomn2(3,1)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Magnetizat.-x) :   [m/|m|]  ',&
&       rhomn2(2),'  at',ri_rhomn2(1,2),ri_rhomn2(2,2),ri_rhomn2(3,2)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Magnetizat.-y) :   [m/|m|]  ',&
&       rhomn2(3),'  at',ri_rhomn2(1,3),ri_rhomn2(2,3),ri_rhomn2(3,3)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Magnetizat.-z) :   [m/|m|]  ',&
&       rhomn2(4),'  at',ri_rhomn2(1,4),ri_rhomn2(2,4),ri_rhomn2(3,4)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Minimum (Spin pol zeta) :   [m/|m|]  ',&
&       zetmn2(1),'  at',ri_zetmn2(1,1),ri_zetmn2(2,1),ri_zetmn2(3,1)
       call wrtout(iout,message,'COLL')
     end if
     write(message,*)' '
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Total  el-den) : [el/Bohr^3]',&
&     rhomx1(1),'  at',ri_rhomx1(1,1),ri_rhomx1(2,1),ri_rhomx1(3,1)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Magnetizat.-x) :   [m/|m|]  ',&
&     rhomx1(2),'  at',ri_rhomx1(1,2),ri_rhomx1(2,2),ri_rhomx1(3,2)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Magnetizat.-y) :   [m/|m|]  ',&
&     rhomx1(3),'  at',ri_rhomx1(1,3),ri_rhomx1(2,3),ri_rhomx1(3,3)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Magnetizat.-z) :   [m/|m|]  ',&
&     rhomx1(4),'  at',ri_rhomx1(1,4),ri_rhomx1(2,4),ri_rhomx1(3,4)
     call wrtout(iout,message,'COLL')
     write(message,'(a,es13.4,a,3f10.4)')'      Maximum (Spin pol zeta) :   [m/|m|]  ',&
&     zetmx1(1),'  at',ri_zetmx1(1,1),ri_zetmx1(2,1),ri_zetmx1(3,1)
     call wrtout(iout,message,'COLL')
     if(option==2)then
       write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Total  el-den) : [el/Bohr^3]',&
&       rhomx2(1),'  at',ri_rhomx2(1,1),ri_rhomx2(2,1),ri_rhomx2(3,1)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Magnetizat.-x) :   [m/|m|]  ',&
&       rhomx2(2),'  at',ri_rhomx2(1,2),ri_rhomx2(2,2),ri_rhomx2(3,2)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Magnetizat.-y) :   [m/|m|]  ',&
&       rhomx2(3),'  at',ri_rhomx2(1,3),ri_rhomx2(2,3),ri_rhomx2(3,3)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Magnetizat.-z) :   [m/|m|]  ',&
&       rhomx2(4),'  at',ri_rhomx2(1,4),ri_rhomx2(2,4),ri_rhomx2(3,4)
       call wrtout(iout,message,'COLL')
       write(message,'(a,es13.4,a,3f10.4)')' Next-Maximum (Spin pol zeta) :   [m/|m|]  ',&
&       zetmx2(1),'  at',ri_zetmx2(1,1),ri_zetmx2(2,1),ri_zetmx2(3,1)
       call wrtout(iout,message,'COLL')
     end if
   end if
 end if

 ABI_DEALLOCATE(coord)
 ABI_DEALLOCATE(value)
 ABI_DEALLOCATE(iindex)
 ABI_DEALLOCATE(integrated)

end subroutine prtrhomxmn
!!***
