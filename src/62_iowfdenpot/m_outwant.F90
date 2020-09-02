!!****m* ABINIT/m_outwant
!! NAME
!! m_outwant
!!
!! FUNCTION
!! Interface with want code.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2020 ABINIT group (CMorari)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_outwant

 use defs_basis
 use m_errors
 use m_abicore
 use m_hdr
 use m_dtset

 use m_io_tools,   only : open_file
 use m_symtk,      only : matr3inv

 implicit none

 private
 public ::  outwant

contains
!!***

!!****f* m_outwant/outwant
!! NAME
!! outwant
!!
!! FUNCTION
!! This routine creates an output file containing all the
!! information needed to run WanT as a post-processing program
!! The resulting file is 'launch.dat'.
!!
!! The routine writes to the disk (unformatted file unitwnt) the following information:
!!
!!     alat - lattice parameter
!!     rprim - primitive translation vectors
!!     ntypat - nr of atom types in elementary cell
!!     tpat - nr of types of atoms in the elementary cell
!!     xcart - cartesian coordinates of the atoms in the elem. cell
!!     ecut - energy cut-off
!!     mband - # of bands taken in calculation (same for each K-pt)
!!     nk(3) - # of k-pts for each direction (uniform grid in the WHOLE BZ)
!!     s0(3) - the origin of the K-space
!!     kg_tmp(3,mpw*mkmem ) - reduced planewave coordinates
!!     imax - Maximum index  of a G vector among all k points (see explanation bellow)
!!     nkpt - total no of K-pts
!!     nsppol - nr of spin polarisations (1 or 2)
!!     eig(mband, nkp_tot) - eigenvalues/band/K_point
!!     ngfft(3) - nr of points used for FFT in each direction
!!     wfc(i)- cmplx(cg(1,i),cg(2,i)) - wavefunction
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eig(mband*nkpt*nsppol) = array for holding eigenvalues (Hartree)
!!  cg(2,mcg) = planewave coefficients of wavefunction
!!  kg(3, mpw*mkmem) = reduced planewave coordinates
!!  npwarr(nkpt) = number of planewaves in basis at this k-point
!!  mband = maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  nkpt = number of k - points
!!  nsppol = 1 for unpolarized, 2 for spin polarized
!!  nspinor = number of spinorial components of the wavefunction (on current proc)
!!  mkmem = number of k points treated by this node.
!!  mpw = maximum dimensioned size of npw
!!  prtwant = if set to 1, print 0 in S0 output
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      m_outscfcv
!!
!! CHILDREN
!!      matr3inv,wrtout
!!
!! SOURCE

subroutine outwant(dtset,eig,cg,kg,npwarr,mband,mcg,nkpt,nsppol,mkmem,mpw,prtwant)

!Arguments ------------------------------------
!scalars
 integer :: mband,mcg,mkmem,mpw,nkpt,nsppol,prtwant
 type(dataset_type),intent(in) :: dtset
!arrays
 integer :: kg(3,mpw*mkmem),npwarr(nkpt)
 real(dp) :: cg(2,mcg),eig(mband*nkpt*nsppol)

!Local variables-------------------------------
! the following variables are not used; they are written to 'launch.dat'
! in order to be compatible with the WANT format
!scalars
 integer :: bandtot,i,icount,ifind,ig,ii,iij,ik,ik_,imax
 integer :: index,index1,ispin_,iunit,iwf,iwf_k,j,k
 integer :: maxat,ngm,ngw_,nk_,nkp,nspin_
 integer :: unitwnt
 real(dp) :: alat,scal,scal_,tt
 logical :: twrite=.true.
 character(len=20) :: section_name
 character(len=3) :: nameat
 character(len=500) :: message
 character(len=fnlen) :: filewnt
!arrays
 integer :: ikg(3),nk(3)
 integer,allocatable :: iwfi(:,:),kg_tmp(:,:),tpat(:)
 real(dp) :: drprim(3,3),gmat(3,3),gmod(3),s0(3),t1(3),t2(3)
 real(dp),allocatable :: xcoord(:,:,:)
 complex,allocatable :: wfc(:)

! ***************************************************************************

!WARNING: not tested for nsppol,nspinor >1
!
!Initialisations
 nameat = ' '
 bandtot=mband*nkpt*nsppol*dtset%nspinor
 filewnt='launch.dat'

 write(message,'(3a)')ch10,' Opening file for WanT input: ',trim(filewnt)
 call wrtout(std_out,message,'COLL')

!Open the file
 if (open_file(filewnt,message,newunit=unitwnt,form='unformatted', status='unknown') /=0) then
   MSG_ERROR(message)
 end if

!Comments
 if(prtwant>1) then
   write(std_out,*) 'Wrong value for prtwant. Reseting to 1'
   prtwant=1
 elseif(prtwant==1) then
   do i=1,3
     s0(i)=0._dp
   end do
 end if

!Discussion of 'alat' ABINIT/ WanT
 if(dtset%acell_orig(1,1)==dtset%acell_orig(2,1).and.&
 dtset%acell_orig(1,1)==dtset%acell_orig(3,1)) then
   alat=dtset%acell_orig(1,1)
   do i=1,3
     do j=1,3
       drprim( i, j) = dtset%rprim_orig( i, j, 1 )
     end do
   end do
 else
!  Redefining the drprim( i, j)
   alat=dtset%acell_orig(1,1)
   do i=1,3
     do j=1,3
       drprim( i, j) = dtset%rprim_orig( i, j, 1 )*dtset%acell_orig(j, 1)/alat
     end do
   end do
 end if

!Now finding the no of k-pt for each direction PARALEL with the
!generators of the first B.Z.
!First decide if we have the Gamma point in the list; its index in the list is ... index
 nk(:)=1
 ifind=0
 icount=2
 do i=1,nkpt
   index1=0
   do j=1,3
     if(dtset%kptns(j,i)<tol8) index1=index1+1
   end do
   if(index1==3) then
     index=i
     ifind=1
     cycle
   end if
 end do
 if(ifind==0) then
   write(std_out,*) 'GAMMA POINT NOT IN THE LIST OF KPTS?'
   do ii=1,nkpt
     write(std_out,*) (dtset%kptns(j,ii),j=1,3)
   end do
   MSG_ERROR("fatal error")
 end if

 call matr3inv(drprim,gmat)

!Modules for each vector in recip. space; nb: g(index coord, index point)
 do j=1,3
   gmod(j)=0.D0
   do i=1,3
     gmod(j)=gmod(j)+gmat(i,j)**2
   end do
   gmod(j)=sqrt(gmod(j))
 end do
 if(nkpt==2) then
   do j=1,3
     do ii=1,3
       t1(ii)=dtset%kptns(ii,1)-dtset%kptns(ii,2)
     end do
     tt=0._dp
     do iij=1,3
       t2(iij)=0._dp
       do ii=1,3
         t2(iij)=t2(iij)+t1(ii)*gmat(ii,iij)
       end do
       tt=tt + t2(iij)**2
     end do
     tt=sqrt(tt)
     scal=0._dp
     do ii=1,3
       scal=scal+t2(ii)*gmat(j,ii)
     end do
     scal=abs(scal)
!    Compare scal(tt,gmat) with simple product of modules -> paralel or not
     if(abs(scal-tt*gmod(j))<tol8) nk(j)=2
   end do

 elseif(nkpt>2) then

   do i=1,nkpt
     if(i.ne.index) then
       do ii=1,3
         t1(ii)=dtset%kptns(ii,index)-dtset%kptns(ii,i)
       end do
       tt=0._dp
       do iij=1,3
         t2(iij)=0._dp
         do ii=1,3
           t2(iij)=t2(iij)+t1(ii)*gmat(ii,iij)
         end do
         tt=tt + t2(iij)**2
       end do
       tt=sqrt(tt)
!      check for each direction in the BZ
       do j=1,3
         scal=0._dp
         do ii=1,3
           scal=scal+t2(ii)*gmat(j,ii)
         end do
         scal=abs(scal)
!        Compare scal(t1,gmat) with simple product of modules -> paralel or not
         if(abs(scal-tt*gmod(j))<tol8) nk(j)=nk(j)+1
       end do
     end if
   end do
 end if
 index=1
 do i=1,3
   index=index*nk(i)
 end do

 if(index.ne.nkpt) then
   write(message,'(a,2i0)')' OutwanT: Wrong assignemt of kpts', index,nkpt
   MSG_ERROR(message)
 end if

!End counting/assigning no of kpts/direction
!Reordering the coordinates of all atoms - xcoord array
 ABI_ALLOCATE(tpat,(dtset%ntypat))
 tpat(:)=zero
 do i=1,dtset%natom
   do j=1,dtset%ntypat
     if(dtset%typat(i)==j) tpat(j)=tpat(j)+1
   end do
 end do
 maxat=maxval(tpat(:))
 ABI_ALLOCATE(xcoord,(3,maxat,dtset%ntypat))
 index=1
 do i=1, dtset%ntypat
   do k=1,tpat(i)
     do j=1,3
       xcoord(j,k,i)=dtset%xred_orig(j,index,1)
     end do
     index=index+1
   end do
 end do
!
!Defining the kg_tmp list
!Preparing the output of reduced coords., in a single list (kg_tmp(3,imax))
!We start with kg_tmp(:,i)=kg(:,i=1,npwarr(1)) then the new coordinates are added
!ONLY if they are not allready in the list. An index is associated
!for each kg_tmp which allow us to recover kg(3,mpw*nkpt) from
!the smaller list kg_tmp(3, imax)
 ABI_ALLOCATE(kg_tmp,(3,mpw*nkpt))
 ABI_ALLOCATE(iwfi,(nkpt,mpw))
 kg_tmp(:,:)=zero
 iwfi(:,:)=zero
 imax=npwarr(1)
 index=0

 do i=1,  nkpt
   if(i>1) then
     index=index+npwarr(i-1)
   end if
   do j=1, npwarr(i)
     if(i.eq.1) then
       iwfi(i,j)=j
       if(mkmem>0) kg_tmp(:,j)=kg(:,j)
     else
       ifind=0
       if(mkmem>0) ikg(:)=kg(:,index+j)

       do k=1,imax
         if(ikg(1)==kg_tmp(1,k)) then
           if(ikg(2)==kg_tmp(2,k)) then
             if(ikg(3)==kg_tmp(3,k)) then
               ifind=1
               iwfi(i,j)=k
             end if
           end if
         end if
       end do

       if(ifind==0) then
         imax=imax+1
         kg_tmp(:,imax)=ikg(:)
         iwfi(i,j)=imax
       end if
     end if
   end do
 end do
 ngm=imax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!PART ONE: writing the header
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(unitwnt) alat
 write(unitwnt) ( drprim( i, 1 ), i = 1, 3 )  ! save A1
 write(unitwnt) ( drprim( i, 2 ), i = 1, 3 )  ! save A2
 write(unitwnt) ( drprim( i, 3 ), i = 1, 3 )  ! save A3
!write(std_out,* ) ( drprim( i, 1 ), i = 1, 3 )  ! save A1
!write(std_out,* ) ( drprim( i, 2 ), i = 1, 3 )  ! save A2
!write(std_out,* ) ( drprim( i, 3 ), i = 1, 3 )  ! save A3

 write(unitwnt) dtset%ntypat
!write(std_out,*) dtset%ntypat, 'NTYPAT', tpat

 do i = 1, dtset%ntypat
   write(unitwnt) tpat(i), nameat
   write(unitwnt) ((xcoord(j,k,i),j=1,3), k=1, tpat(i))
!  write(std_out,*) tpat(i), nameat
!  write(std_out,*) ((xcoord(j,k,i),j=1,3),k=1,tpat(i)), 'XCART'
 end do
 ABI_DEALLOCATE(tpat)
 ABI_DEALLOCATE(xcoord)

!energy cut-off in Rydberg (WANT option)
 write (unitwnt) 2._dp*dtset%ecut, mband
!write(std_out,*)   2._dp*dtset%ecut, mband
 write (unitwnt) ( nk(i), i = 1, 3 ), ( s0(j), j = 1, 3 ),ngm
!write(std_out,*) ( nk(i), i = 1, 3 ), ( s0(j), j = 1, 3 ),imax
 write (unitwnt) ( kg_tmp( 1, i ), kg_tmp( 2, i ), kg_tmp( 3, i ), i = 1, ngm )
 write (unitwnt) mpw, mband, dtset%nkpt/dtset%nsppol
!write(std_out,*) mpw, mband,  dtset%nkpt/dtset%nsppol

 do i=1, nkpt
   write(unitwnt) (iwfi(i,j), j=1,mpw)
 end do
 ABI_DEALLOCATE(kg_tmp)

!Eigenvalues in HARTREE
 write (unitwnt)  ( eig( i ), i = 1, bandtot)
 write (unitwnt) ( npwarr( ik ), ik = 1, nkpt )
 write (unitwnt) ( mband, ik = 1, nkpt )
 write (unitwnt) (dtset%ngfft(i),i=1,3), imax, imax
!write(std_out,*)  ( eig( i ), i = 1, bandtot )
!write(std_out,*) ( npwarr( ik ), ik = 1, nkpt )
!write(std_out,*) ( mband, ik = 1, nkpt )
!write(std_out,*) (dtset%ngfft(i),i=1,3), imax ,imax
!a list with the band structure; usefull for 'windows' and 'disentangle' programs
!from WanT distribution

 if (open_file('band.gpl',message,newunit=iunit,status='unknown') /=0) then
   MSG_ERROR(message)
 end if

 index=1
 do i=1,mband
   index=1
   do j=1,nkpt
     write(iunit,*) index, Ha_eV*eig(i+(j-1)*mband), eig(i+(j-1)*mband)
     index=index+1
   end do
   write(iunit,*)
 end do

 close(iunit)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!PART TWO: Writing the wavefunction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!not used
 ngw_=0
 ik_=0
 nk_=0
 ispin_=0
 nspin_=0
 scal_=1._dp
!!!!!!!!!!!!!!!!!!!!!!!!!!
 iwf = 1
 iwf_k=1
 ABI_ALLOCATE(wfc,(imax))

!Loop over k-pt
 do nkp=1,nkpt
!  Not relevant
   write(unitwnt) twrite, ik_, section_name
!  Only 'mband' is relevant here
   write(unitwnt) ngw_, mband, ik_, nk_, nk_,ispin_, nspin_, scal_
   write(unitwnt) imax
!  Not relevant
   write(unitwnt) twrite
!  Loop over bands

!  Preparing WF
   do k=1,mband
     if(mkmem >0) then
       wfc(:)=zero
!      From cg to wf:
       do i=iwf, iwf+npwarr(nkp)-1
         index=i-iwf+1
         wfc(iwfi(nkp,index))=cmplx(cg(1,i), cg(2,i), kind(0._dp))
       end do
       iwf=iwf+npwarr(nkp)
     else
       message = 'Wrong mkmem in outwant'
       MSG_ERROR(message)
     end if
     write(unitwnt) (wfc(ig), ig=1,imax)
!    End loop over bands
   end do

!  Not relevant
   write(unitwnt) twrite
!  Not relevant
   do i=1,mband
     write(unitwnt) i
   end do

!  End loop over k-pts
 end do

 ABI_DEALLOCATE(iwfi)
 ABI_DEALLOCATE(wfc)

 call wrtout(std_out,'Closing file','COLL')
 close(unit=unitwnt)

end subroutine outwant
!!***

end module m_outwant
!!***
