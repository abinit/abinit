!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pawang
!! NAME
!!  m_pawang
!!
!! FUNCTION
!!  This module contains the definition of the pawang_type structured datatype,
!!  as well as related functions and methods.
!!  pawang_type variables define ANGular mesh discretization of PAW augmentation
!!  regions and related data.
!!
!! COPYRIGHT
!! Copyright (C) 2013-2019 ABINIT group (MT, FJ, BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

MODULE m_pawang

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 use m_paw_sphharm, only : initylmr, mat_mlms2jmj, mat_slm2ylm

 implicit none

 private

!public procedures.
 public :: pawang_init        ! Constructor
 public :: pawang_free        ! Free memory
 public :: pawang_lsylm       ! Compute the LS operator in the real spherical harmonics basis
 public :: initang            ! Initialize angular mesh for PAW calculations
 public :: make_angular_mesh  ! Build angular mesh from ntheta, nphi.

 ! MGPAW: Private?
 public :: realgaunt          ! compute "real Gaunt coefficients" with "real spherical harmonics"
 public :: gaunt              ! Gaunt coeffients for complex Yml
 public :: gauleg
 public :: mat_mlms2jmj
 public :: mat_slm2ylm

!!***

!----------------------------------------------------------------------


!!****t* m_pawang/pawang_type
!! NAME
!! pawang_type
!!
!! FUNCTION
!! For PAW: ANGular mesh discretization of PAW augmentation regions and related data.
!!
!! SOURCE

 type,public :: pawang_type

!Integer scalars

  integer :: angl_size=0
   ! Dimension of paw angular mesh
   ! angl_size=ntheta*nphi

  integer :: l_max=-1
   ! Maximum value of angular momentum l+1

  integer :: l_size_max=-1
   ! Maximum value of angular momentum +1
   ! leading to non-zero Gaunt coefficients
   ! l_size_max = 2*l_max-1

  integer :: ngnt=0
   ! Number of non zero Gaunt coefficients

  integer :: ntheta, nphi
   ! Dimensions of paw angular mesh

  integer :: nsym
   ! Number of symmetry elements in space group

  integer :: gnt_option=-1
   ! Option for Gaunt coefficients:
   !  gnt_option==0, Gaunt coeffs are not computed (and not allocated)
   !  gnt_option==1, Gaunt coeffs are computed up to 2*l_size_max-1
   !  gnt_option==2, Gaunt coeffs are computed up to l_size_max

  integer :: use_ls_ylm=0
   ! Flag: use_ls_ylm=1 if pawang%ls_ylm is allocated

  integer :: ylm_size=0
   ! Size of ylmr/ylmrgr arrays

!Integer arrays

  integer, allocatable :: gntselect(:,:)
   ! gntselect(l_size_max**2,l_max**2*(l_max**2+1)/2)
   ! Selection rules for Gaunt coefficients stored as (LM,ij) where ij is in packed form.
   ! (if gntselect>0, Gaunt coeff. is non-zero)

!Real (real(dp)) arrays

  real(dp), allocatable :: anginit(:,:)
   ! anginit(3,angl_size)
   ! For each point of the angular mesh, gives the coordinates
   ! of the corresponding point on an unitary sphere
   ! Not used in present version (5.3)

  real(dp), allocatable :: angwgth(:)
   ! angwgth(angl_size)
   ! For each point of the angular mesh, gives the weight
   ! of the corresponding point on an unitary sphere

  real(dp), allocatable :: ls_ylm(:,:,:)
   ! ls_ylm(2,l_max**2*(l_max**2+1)/2,2)
   ! LS operator in the real spherical harmonics basis
   ! ls_ylm(ilm1m2,ispin)= <sigma, y_lm1| LS |y_lm2, sigma_prime>

  real(dp), allocatable :: realgnt(:)
   ! realgnt(ngnt)
   ! Non zero real Gaunt coefficients

  real(dp), allocatable :: ylmr(:,:)
   ! ylmr(ylm_size,angl_size)
   ! Real Ylm calculated in real space

  real(dp), allocatable :: ylmrgr(:,:,:)
   ! ylmrgr(1:3,ylm_size,angl_size)
   ! First gradients of real Ylm calculated in real space (cart. coordinates)

  real(dp), allocatable :: zarot(:,:,:,:)
   !  zarot(l_size_max,l_size_max,l_max,nsym)
   !  Coeffs of the transformation of real spherical
   !  harmonics under the symmetry operations symrec.

 end type pawang_type
!!***

CONTAINS

!===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_pawang/pawang_init
!! NAME
!! pawang_init
!!
!! FUNCTION
!!  Initialize a pawang datastructure
!!
!! INPUTS
!!  gnt_option=flag activated if pawang%gntselect and pawang%realgnt have to be allocated
!!             also determine the size of these pointers
!!  lmax=maximum value of angular momentum l
!!  nphi,ntheta=dimensions of paw angular mesh
!!  nsym=number of symetries
!!  pawxcdev=choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  use_ls_ylm=flag activated if pawang%ls_ylm has to be allocated
!!  use_ylm=flag activated if pawang%ylmr and pawang%ylmrgr have to be allocated
!!  xclevel=XC functional level (1=LDA, 2=GGA)
!!
!! OUTPUT
!!  Pawang <type(pawang_type)>=ANGular mesh discretization and related data
!!
!! PARENTS
!!      dfpt_looppert,pawinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawang_init(Pawang,gnt_option,lmax,nphi,nsym,ntheta,pawxcdev,use_ls_ylm,use_ylm,xclevel)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: gnt_option,lmax,nphi,nsym,ntheta
 integer,intent(in) :: pawxcdev,use_ls_ylm,use_ylm,xclevel
 type(Pawang_type),intent(inout) :: Pawang

!Local variables-------------------------------
!scalars
 integer :: ll,sz1,sz2,sz3
!arrays
 real(dp),allocatable :: rgnt_tmp(:)

! *************************************************************************

 !@Pawang_type

 Pawang%l_max=lmax+1
 Pawang%l_size_max=2*Pawang%l_max-1
 Pawang%nsym=nsym

 if (pawxcdev==0) then
   Pawang%nphi=nphi
   Pawang%ntheta=ntheta
   Pawang%angl_size=ntheta*nphi
 else
   Pawang%nphi=0
   Pawang%ntheta=0
   Pawang%angl_size=0
 end if

 if (Pawang%angl_size>0) then
   LIBPAW_ALLOCATE(Pawang%anginit,(3,Pawang%angl_size))
   LIBPAW_ALLOCATE(Pawang%angwgth,(Pawang%angl_size))
   call initang(pawang)
 end if

 if (use_ylm>0.and.Pawang%angl_size>0) then
   if (xclevel>=0) ll=Pawang%l_size_max
   if (xclevel>=2) ll=ll+1
   Pawang%ylm_size=ll**2
   LIBPAW_ALLOCATE(Pawang%ylmr,(Pawang%ylm_size,Pawang%angl_size))
   if (xclevel==2) then
     LIBPAW_ALLOCATE(Pawang%ylmrgr,(3,Pawang%ylm_size,Pawang%angl_size))
     call initylmr(ll,0,pawang%angl_size,pawang%angwgth,2,pawang%anginit,pawang%ylmr,&
&                  ylmr_gr=pawang%ylmrgr)
   else
     call initylmr(ll,0,pawang%angl_size,pawang%angwgth,1,pawang%anginit,pawang%ylmr)
   end if
 else
   Pawang%ylm_size=0
 end if

 Pawang%gnt_option=gnt_option
 if (Pawang%gnt_option==1.or.Pawang%gnt_option==2) then
   if (Pawang%gnt_option==1) then
     sz1=(2*Pawang%l_max-1)**2*(Pawang%l_max)**4
     sz2=(Pawang%l_size_max)**2
     sz3=(Pawang%l_max**2)*(Pawang%l_max**2+1)/2
     LIBPAW_ALLOCATE(rgnt_tmp,(sz1))
     LIBPAW_ALLOCATE(pawang%gntselect,(sz2,sz3))
     call realgaunt(Pawang%l_max,Pawang%ngnt,Pawang%gntselect,rgnt_tmp)
   else if (Pawang%gnt_option==2) then
     sz1=(4*Pawang%l_max-3)**2*(2*Pawang%l_max-1)**4
     sz2=(2*Pawang%l_size_max-1)**2
     sz3=((2*Pawang%l_max-1)**2)*((2*Pawang%l_max-1)**2+1)/2
     LIBPAW_ALLOCATE(rgnt_tmp,(sz1))
     LIBPAW_ALLOCATE(pawang%gntselect,(sz2,sz3))
     call realgaunt(2*Pawang%l_max-1,Pawang%ngnt,Pawang%gntselect,rgnt_tmp)
   end if
   if (allocated(pawang%realgnt))  then
     LIBPAW_DEALLOCATE(pawang%realgnt)
   end if
   LIBPAW_ALLOCATE(Pawang%realgnt,(Pawang%ngnt))
   Pawang%realgnt(1:Pawang%ngnt)=rgnt_tmp(1:Pawang%ngnt)
   LIBPAW_DEALLOCATE(rgnt_tmp)
 end if

 Pawang%use_ls_ylm=use_ls_ylm
 if (use_ls_ylm>0) then
   LIBPAW_ALLOCATE(pawang%ls_ylm,(2,Pawang%l_max**2*(Pawang%l_max**2+1)/2,2))
   call pawang_lsylm(pawang)
 end if

 if (nsym>0) then
   LIBPAW_ALLOCATE(Pawang%zarot,(Pawang%l_size_max,Pawang%l_size_max,Pawang%l_max,nsym))
 end if

end subroutine pawang_init
!!***

!----------------------------------------------------------------------

!!****f* m_pawang/pawang_free
!! NAME
!! pawang_free
!!
!! FUNCTION
!!  Free all dynamic memory and reset all flags stored in a pawang datastructure
!!
!! SIDE EFFECTS
!!  Pawang <type(pawang_type)>=ANGular mesh discretization and related data
!!
!! PARENTS
!!      dfpt_looppert,driver,pawinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawang_free(Pawang)

!Arguments ------------------------------------
!scalars
 type(Pawang_type),intent(inout) :: Pawang

! *************************************************************************

 !@Pawang_type
 if (allocated(pawang%angwgth))    then
   LIBPAW_DEALLOCATE(pawang%angwgth)
 end if
 if (allocated(pawang%anginit))    then
   LIBPAW_DEALLOCATE(pawang%anginit)
 end if
 if (allocated(pawang%zarot))      then
   LIBPAW_DEALLOCATE(pawang%zarot)
 end if
 if (allocated(pawang%gntselect))  then
   LIBPAW_DEALLOCATE(pawang%gntselect)
 end if
 if (allocated(pawang%realgnt))    then
   LIBPAW_DEALLOCATE(pawang%realgnt)
 end if
 if (allocated(pawang%ylmr))       then
   LIBPAW_DEALLOCATE(pawang%ylmr)
 end if
 if (allocated(pawang%ylmrgr))     then
   LIBPAW_DEALLOCATE(pawang%ylmrgr)
 end if
 if (allocated(pawang%ls_ylm))     then
   LIBPAW_DEALLOCATE(pawang%ls_ylm)
 end if

 pawang%angl_size =0
 pawang%ylm_size  =0
 pawang%use_ls_ylm=0
 pawang%l_max=-1
 pawang%l_size_max=-1
 pawang%gnt_option=-1
 pawang%ngnt=0

end subroutine pawang_free
!!***

!----------------------------------------------------------------------

!!****f* m_pawang/pawang_lsylm
!! NAME
!! pawang_lsylm
!!
!! FUNCTION
!! Compute the LS operator in the real spherical harmonics basis
!! ls_ylm(ilm1,ilm2,ispin)= <sigma, S_lm1| L.S |S_lm2, sigma_prime>
!!   ilm,1m2=(l,m1,m2) with -l<=m1<=l, -l<=m2<=l and 0<l<=lmax
!!   ispin=(sigma,sigma_prime) 1=(up,up), 2=(up,dn), 3=(dn,up), 4=(dn,dn)
!!
!! INPUTS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!
!! OUTPUT
!!  pawang%ls_ylm(2,l_max**2*(l_max**2+1)/2,2)=LS operator in the real spherical harmonics basis
!!        ls_ylm(:,:,1)=<up, S_lm1| L.S |S_lm2, up>
!!        ls_ylm(:,:,2)=<up, S_lm1| L.S |S_lm2, down>
!!        One can deduce:
!!        <down, S_lm1| L.S |S_lm2, down>=-<up, S_lm1| L.S |S_lm2, up>
!!        <down, S_lm1| L.S |S_lm2, up>  =-Conjg[<up, S_lm1| L.S |S_lm2, down>]
!!        Also, only ilm1<=ilm2 terms are stored, because:
!!         <sigma, S_lm1| L.S |S_lm2, sigma_prime>=-<sigma_prime, S_lm1| L.S |S_lm2, sigma>
!!
!! PARENTS
!!      m_pawang
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawang_lsylm(pawang)

!Arguments ---------------------------------------------
!scalars
 type(pawang_type),intent(inout) :: pawang

!Local variables ---------------------------------------
!scalars
 integer :: ii,ilm,im,j0lm,jj,jlm,jm,klm,l_max,ll,lm0,mm,ispden
 real(dp),parameter :: invsqrt2=one/sqrt2
 real(dp) :: onem
 character(len=500) :: msg
 logical,parameter :: tso=.false. ! use true to Test Spin Orbit and
!                                   write the matrix of L.S in different basis
!arrays
 complex(dpc) :: tmp(2)
 complex(dpc),allocatable :: ls_cplx(:,:,:),slm2ylm(:,:)
 complex(dpc),allocatable :: mat_inp_c(:,:,:),mat_out_c(:,:,:)
 complex(dpc),allocatable :: mat_ls_ylm(:,:,:),mat_jmj(:,:)
 character(len=9),parameter :: dspin2(2)=(/"up-up    ","up-dn    "/)
 character(len=9),parameter :: dspin6(6)=(/"dn       ","up       ","dn-dn    ","up-up    ","dn-up    ","up-dn    "/)
 character(len=9),parameter :: dspinm(6)=(/"dn       ","up       ","n        ","mx       ","my       ","mz       "/)

! *************************************************************************

 if (pawang%use_ls_ylm==0) then
   msg='  ls_ylm pointer is not allocated !'
   MSG_BUG(msg)
 end if

!Initialization
 pawang%ls_ylm=zero
 l_max=pawang%l_max-1

!Nothing to do if lmax=0
 if (l_max<=0) return

!Loop on l quantum number
 do ll=1,l_max

!  Transformation matrixes: real->complex spherical harmonics
   LIBPAW_ALLOCATE(slm2ylm,(2*ll+1,2*ll+1))
   slm2ylm=czero
   do im=1,2*ll+1
     mm=im-ll-1;jm=-mm+ll+1
     onem=dble((-1)**mm)
     if (mm> 0) then
       slm2ylm(im,im)= cmplx(onem*invsqrt2,zero,kind=dp)
       slm2ylm(jm,im)= cmplx(invsqrt2,     zero,kind=dp)
     end if
     if (mm==0) then
       slm2ylm(im,im)=cone
     end if
     if (mm< 0) then
       slm2ylm(im,im)= cmplx(zero,     invsqrt2,kind=dp)
       slm2ylm(jm,im)=-cmplx(zero,onem*invsqrt2,kind=dp)
     end if
   end do

!  Compute <sigma, Y_lm1|L.S|Y_lm2, sigma_prime> (Y_lm=complex spherical harmonics)
!  1= <up|L.S|up>  ;  2= <up|L.S|dn>
   LIBPAW_ALLOCATE(ls_cplx,(2*ll+1,2*ll+1,2))
   ls_cplx=czero
   if(tso)  then
     LIBPAW_ALLOCATE(mat_ls_ylm,(2*ll+1,2*ll+1,4))
     if(tso) mat_ls_ylm=czero
   end if
   if(tso)  then
     LIBPAW_ALLOCATE(mat_jmj,(2*(2*ll+1),2*(2*ll+1)))
     if(tso) mat_jmj=czero
   end if
   do im=1,2*ll+1
     mm=im-ll-1
     ls_cplx(im,im,1)=half*mm
     if(tso) mat_ls_ylm(im,im,1)=-half*mm ! dn dn
     if(tso) mat_ls_ylm(im,im,2)=half*mm  ! up up
     if ((mm+1)<= ll) then
       ls_cplx(im,im+1,2)=half*sqrt(real((ll-mm)*(ll+mm+1),kind=dp))
       if(tso) mat_ls_ylm(im,im+1,4)=half*sqrt(real((ll-mm)*(ll+mm+1),kind=dp))  ! up dn
       if(tso) mat_ls_ylm(im+1,im,3)=half*sqrt(real((ll-mm)*(ll+mm+1),kind=dp))  ! dn up
     end if
     if ((mm-1)>=-ll) then
       ls_cplx(im-1,im,2)=half*sqrt(real((ll+mm)*(ll-mm+1),kind=dp))
       if(tso) mat_ls_ylm(im-1,im,4)=half*sqrt(real((ll+mm)*(ll-mm+1),kind=dp))  ! up dn
       if(tso) mat_ls_ylm(im,im-1,3)=half*sqrt(real((ll+mm)*(ll-mm+1),kind=dp))  ! dn up
     end if
   end do

!  test : print LS in J,M_J basis
   if(tso) then
     do ispden=1,4
       write(msg,'(3a)') ch10,"value of LS in the Ylm basis for " ,trim(dspin6(ispden+2*(4/4)))
       call wrtout(std_out,msg,'COLL')
       do im=1,ll*2+1
         write(msg,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))') (mat_ls_ylm(im,jm,ispden),jm=1,ll*2+1)
         call wrtout(std_out,msg,'COLL')
       end do
     end do
     call mat_mlms2jmj(ll,mat_ls_ylm,mat_jmj,4,1,2,3,std_out,'COLL')  ! optspin=2 : dn spin are first
   end if

!  Compute <sigma, S_lm1|L.S|S_lm2, sigma_prime> (S_lm=real spherical harmonics)
!  1= <up|L.S|up>  ;  2= <up|L.S|dn>
   if(tso) then
     LIBPAW_ALLOCATE(mat_inp_c,(2*ll+1,2*ll+1,4))
     LIBPAW_ALLOCATE(mat_out_c,(2*ll+1,2*ll+1,4))
   end if
   lm0=ll**2
   do jm=1,2*ll+1
     jlm=lm0+jm;j0lm=jlm*(jlm-1)/2
     do im=1,jm
       ilm=lm0+im;klm=j0lm+ilm
       tmp(:)=czero
       do ii=1,2*ll+1
         do jj=1,2*ll+1
           tmp(:)=tmp(:)+ls_cplx(ii,jj,:)*CONJG(slm2ylm(ii,im))*slm2ylm(jj,jm)
         end do
       end do
       pawang%ls_ylm(1,klm,:)=REAL(tmp(:),kind=dp)
       pawang%ls_ylm(2,klm,:)=AIMAG(tmp(:))
     end do
   end do

!  Test: print LS in Slm basis
   if(tso) then
     call mat_slm2ylm(ll,mat_ls_ylm,mat_inp_c,4,2,2,3,std_out,'COLL') ! from Ylm to Slm, and dn spin are first
     do ispden=1,4
       write(msg,'(3a)') ch10,"value of LS in the Slm basis for " ,trim(dspin6(ispden+2*(4/4)))
       call wrtout(std_out,msg,'COLL')
       do im=1,ll*2+1
         write(msg,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))') (mat_inp_c(im,jm,ispden),jm=1,ll*2+1)
         call wrtout(std_out,msg,'COLL')
       end do
     end do
!    change into n,m basis
     mat_ls_ylm(:,:,1)=(mat_inp_c(:,:,1)+mat_inp_c(:,:,2))
     mat_ls_ylm(:,:,2)=(mat_inp_c(:,:,3)+mat_inp_c(:,:,4))
     mat_ls_ylm(:,:,3)=-cmplx(0.d0,1.d0)*(mat_inp_c(:,:,4)-mat_inp_c(:,:,3))
     mat_ls_ylm(:,:,4)=(mat_inp_c(:,:,1)-mat_inp_c(:,:,2))
     do ispden=1,4
       write(msg,'(3a)') ch10,"value of LS in the Slm basis for " ,trim(dspinm(ispden+2*(4/4)))
       call wrtout(std_out,msg,'COLL')
       do im=1,ll*2+1
         write(msg,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))') (mat_ls_ylm(im,jm,ispden),jm=1,ll*2+1)
         call wrtout(std_out,msg,'COLL')
       end do
     end do
     LIBPAW_DEALLOCATE(mat_inp_c)
     LIBPAW_DEALLOCATE(mat_ls_ylm)
     LIBPAW_DEALLOCATE(mat_jmj)
     LIBPAW_DEALLOCATE(mat_out_c)
   end if ! tso

   LIBPAW_DEALLOCATE(ls_cplx)
   LIBPAW_DEALLOCATE(slm2ylm)

!  End loop on l
 end do

 end subroutine pawang_lsylm
!!***

!----------------------------------------------------------------------

!!****f* m_pawang/initang
!! NAME
!! initang
!!
!! FUNCTION
!! Initialize angular mesh for PAW calculations
!!
!! INPUTS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!       pawang%angl_size  - Total number of sample points in the angular mesh
!!       pawang%ntheta     - Number of sample points in the theta dir
!!       pawang%nphi       - Number of sample points in the phi dir
!!
!! OUTPUT
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!       pawang%anginit    - (3 x angl_size) array, the ntheta*nphi
!!                           dimensional arrays ax, ay, and az
!!       pawang%angwgth    - (angl_size) array, the weight factor of the
!!                           point (ax, ay, az)
!!
!! PARENTS
!!      m_pawang
!!
!! CHILDREN
!!
!! SOURCE

 subroutine initang(pawang)

!Arguments ------------------------------------
!scalars
 type(pawang_type),intent(inout) :: pawang

!Local variables-------------------------------
!scalars
 integer :: ip,it,npoints
 real(dp) :: ang,con,cos_phi,cos_theta,sin_phi,sin_theta
 character(len=500) :: msg
!arrays
 real(dp) :: th(pawang%ntheta),wth(pawang%ntheta)

! ***********************************************************************

 if (pawang%angl_size==0) return

!Initializations
 npoints=0
 con=two_pi / pawang%nphi
 call gauleg(-one,one,th,wth,pawang%ntheta)

!We now open two nested do-loops. The first loops through the number
!of theta angles, the second through the number of phi angles (?).
!The two together initialize anginit.

 do it = 1, pawang%ntheta

   cos_theta = th(it)
   sin_theta = sqrt(one - cos_theta*cos_theta)

   do ip = 1, pawang%nphi

     ang = con * (ip-1)
     cos_phi = cos(ang)
     sin_phi = sin(ang)

     npoints = npoints + 1

     pawang%anginit(1, npoints) = sin_theta * cos_phi
     pawang%anginit(2, npoints) = sin_theta * sin_phi
     pawang%anginit(3, npoints) = cos_theta

!    Normalization required
     pawang%angwgth(npoints) = wth(it) / (2 * pawang%nphi)

   end do
 end do

!The following is an error statement that will be generated
!if npoints exceeds nang...
 if (npoints > pawang%angl_size) then
   write(msg, '(a,i4,a,a,i4)' ) &
&   '  anginit%npoints =',npoints,ch10,&
&   '        angl_size =',pawang%angl_size
   MSG_BUG(msg)
 end if

end subroutine initang
!!***

!----------------------------------------------------------------------

!!****f* m_pawang/make_angular_mesh
!! NAME
!! make_angular_mesh
!!
!! FUNCTION
!!  Build angular mesh from ntheta, nphi.
!!
!! INPUTS
!!   ntheta: Number of sample points in the theta dir
!!   nphi: Number of sample points in the phi dir
!!
!! OUTPUT
!!   angl_size: Total number of sample points in the angular mesh. (ntheta * nphi)
!!   vers_cart(3, angl_size): For each point of the angular mesh, gives the Cartesian coordinates
!!     of the corresponding point on an unitary sphere.
!!  angwgth(angl_size):
!!     For each point of the angular mesh, gives the weight of the corresponding point on an unitary sphere.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_angular_mesh(ntheta, nphi, angl_size, vers_cart, angwgth)

!Arguments ------------------------------------
 integer,intent(in) :: ntheta, nphi
 integer,intent(out) :: angl_size
 real(dp),allocatable,intent(out) :: vers_cart(:,:)
 real(dp),allocatable,intent(out) :: angwgth(:)

!Local variables ------------------------------
!scalars
 integer :: it, ip, npoints
 real(dp) :: ang, con, cos_phi, cos_theta, sin_phi, sin_theta
!arrays
 real(dp),allocatable :: th(:),wth(:)

! *************************************************************************

 LIBPAW_ALLOCATE(th, (ntheta))
 LIBPAW_ALLOCATE(wth, (ntheta))

 con = two_pi / nphi
 call gauleg(-one, one, th, wth, ntheta)

 ! Initialize vers_cart and angular weights angwgth
 ! NB: summing over f * angwgth gives the spherical average 1/(4pi) \int domega f(omega)
 angl_size = ntheta * nphi
 LIBPAW_ALLOCATE(vers_cart, (3, angl_size))
 LIBPAW_ALLOCATE(angwgth, (angl_size))
 npoints = 0
 do it = 1, ntheta
   cos_theta = th(it)
   sin_theta = sqrt(one - cos_theta*cos_theta)
   do ip = 1, nphi
     ang = con * (ip-1)
     cos_phi = cos(ang); sin_phi = sin(ang)
     npoints = npoints + 1
     vers_cart(1, npoints) = sin_theta * cos_phi
     vers_cart(2, npoints) = sin_theta * sin_phi
     vers_cart(3, npoints) = cos_theta
     ! Normalization required
     angwgth(npoints) = wth(it) / (two * nphi)
   end do
 end do
 !write(std_out, *)"Sum angwgth: ", sum(angwgth)
 !write(std_out, *)"int sig(theta): ", sum(angwgth * sqrt(one - vers_cart(3, :) **2))
 !write(std_out, *)pi/4

 LIBPAW_DEALLOCATE(th)
 LIBPAW_DEALLOCATE(wth)

end subroutine make_angular_mesh
!!***

!----------------------------------------------------------------------

!!****f* m_pawang/realgaunt
!! NAME
!! realgaunt
!!
!! FUNCTION
!! This routine compute "real Gaunt coefficients", i.e. gaunt
!! coefficients according to "real spherical harmonics"
!!
!! INPUTS
!!  l_max= max. value of ang. momentum l+1;  Gaunt coeffs up to
!!          [(2*l_max-1,m),(l_max,m),(l_max,m)] are computed
!!
!! OUTPUT
!!  gntselect((2*l_max-1)**2,l_max**2*(l_max**2+1)/2)=
!!          selection rules for Gaunt coefficients
!!          if Gaunt coeff. is zero, gntselect=0
!!          if Gaunt coeff. is non-zero, gntselect is the index of
!!                           the coeff. in realgnt(:) array
!!  ngnt= number of non-zero Gaunt coefficients
!!  realgnt((2*l_max-1)**2*l_max**4)= non-zero real Gaunt coefficients
!!
!! PARENTS
!!      m_paw_slater,m_pawang,m_pawpwij
!!
!! CHILDREN
!!
!! SOURCE

subroutine realgaunt(l_max,ngnt,gntselect,realgnt)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: l_max
 integer,intent(out) :: ngnt
!arrays
 integer,intent(out) :: gntselect((2*l_max-1)**2,l_max**2*(l_max**2+1)/2)
 real(dp),intent(out) :: realgnt((2*l_max-1)**2*(l_max)**4)

!Local variables ------------------------------
!scalars
 integer :: ilm1,ilm2,ilmp1,k0lm1,klm1,l1,l2,ll,lp1,m1,m2,mm,mm1,mm2,mm3,mp1
 real(dp) :: c11,c12,c21,c22,c31,c32,fact,realgnt_tmp
!arrays
 integer,allocatable :: ssgn(:)
 type(coeff3_type), allocatable :: coeff(:)

!************************************************************************

! Initialize output arrays with zeros.
gntselect = 0; realgnt = zero

!Compute matrix cc where Sl=cc*Yl (Sl=real sph. harm.)
!------------------------------------------------
 LIBPAW_DATATYPE_ALLOCATE(coeff,(4*l_max-3))
 do ll=1,4*l_max-3
   LIBPAW_ALLOCATE(coeff(ll)%value,(2,2*ll-1,2*ll-1))
   coeff(ll)%value(:,:,:)=zero
   coeff(ll)%value(1,ll,ll)=one
   do mm=1,ll-1
     coeff(ll)%value(1,ll+mm,ll+mm)= (-1._dp)**mm/sqrt(2._dp)
     coeff(ll)%value(1,ll-mm,ll+mm)= ( 1._dp)    /sqrt(2._dp)
     coeff(ll)%value(2,ll+mm,ll-mm)=-(-1._dp)**mm/sqrt(2._dp)
     coeff(ll)%value(2,ll-mm,ll-mm)= ( 1._dp)    /sqrt(2._dp)
   end do
 end do

 LIBPAW_ALLOCATE(ssgn,(l_max**2))
 ssgn(:)=1
 if (l_max>0) then
   do l1=1,l_max-1
     ilm1=1+l1**2+l1
     do m1=-l1,-1
       ssgn(ilm1+m1)=-1
     end do
   end do
 end if

 ngnt=0

!Loop on (lp1,mp1)
!------------------------------------------------
 do lp1=0,l_max-1
   do mp1=-lp1,lp1
     ilmp1=1+lp1**2+lp1+mp1
     k0lm1=ilmp1*(ilmp1-1)/2

!    Loop on (l1,m1)<=(lp1,mp1)
!    ------------------------------------------------
     do l1=0,l_max-1
       do m1=-l1,l1
         ilm1=1+l1**2+l1+m1

         if (ilm1<=ilmp1) then

           klm1=k0lm1+ilm1
           gntselect(:,klm1)=0

!          Loop on (l2,m2)
!          ------------------------------------------------
           do l2=abs(l1-lp1),l1+lp1,2
             do m2=-l2,l2
               ilm2=1+l2**2+l2+m2

!              Real Gaunt coeffs selection rules
!              ------------------------------------------------
               if ((l2<=l1+lp1).and.&
&               (((m1== mp1).and.((m2==0).or.(m2==2*abs(mp1)))).or.&
&               ((m1==-mp1).and.(m2==-abs(m1)-abs(mp1))).or.&
&               ((abs(m1)/=(abs(mp1)).and.&
&               ((m2==ssgn(ilm1)*ssgn(ilmp1)*   (abs(m1)+abs(mp1))).or.&
&               (m2==ssgn(ilm1)*ssgn(ilmp1)*abs(abs(m1)-abs(mp1)))&
               ))))) then

!                Compute selected real Gaunt coefficient
!                ------------------------------------------------
                 realgnt_tmp=zero
                 do mm1=-l1,l1
                   c11=coeff(l1+1)%value(1,l1+mm1+1,l1+m1+1)
                   c12=coeff(l1+1)%value(2,l1+mm1+1,l1+m1+1)
                   do mm2= -lp1,lp1
                     c21=coeff(lp1+1)%value(1,lp1+mm2+1,lp1+mp1+1)
                     c22=coeff(lp1+1)%value(2,lp1+mm2+1,lp1+mp1+1)
                     do mm3= -l2,l2
                       c31=coeff(l2+1)%value(1,l2+mm3+1,l2+m2+1)
                       c32=coeff(l2+1)%value(2,l2+mm3+1,l2+m2+1)
                       fact=c11*c21*c31  -  c12*c22*c31&
&                       -c11*c22*c32  -  c12*c21*c32
                       if((abs(fact)>=tol12).and.(mm3==-mm2-mm1)) &
&                       realgnt_tmp=realgnt_tmp+fact*(-1)**mm2 &
&                       *gaunt(l2,mm3,l1,mm1,lp1,-mm2)
                     end do
                   end do
                 end do

!                Count additional non-zero real Gaunt coeffs
!                ------------------------------------------------
                 if (abs(realgnt_tmp)>=tol12) then
                   ngnt=ngnt+1
                   gntselect(ilm2,klm1)=ngnt
                   realgnt(ngnt)=realgnt_tmp/sqrt(four_pi)
                 end if

!                End loops
!                ------------------------------------------------
               end if
             end do
           end do
         end if
       end do
     end do
   end do
 end do

!Deallocate memory
!------------------------------------------------
 do ll=1,4*l_max-3
   LIBPAW_DEALLOCATE(coeff(ll)%value)
 end do
 LIBPAW_DATATYPE_DEALLOCATE(coeff)
 LIBPAW_DEALLOCATE(ssgn)

end subroutine realgaunt
!!***

!----------------------------------------------------------------------

!!****f* m_pawang/gaunt
!! NAME
!! gaunt
!!
!! FUNCTION
!! Returns gaunt coefficient, i.e.
!!   the integral of Sqrt[4 \pi] Y*(l_i,m_i) Y*(ll,mm) Y(l_j,m_j)
!!   See the 3-j and 6-j symbols by Rotenberg, etc., (Technology Press, 1959), pg.5.
!!
!! INPUTS
!!   ll,mm,l1,l2,m1,m2= six quantum numbers defining the Gaunt coef.
!!
!! OUTPUT
!!   gaunt(ll,mm,l1,l2,m1,m2)=the value of the integral
!!
!! CHILDREN
!!
!! SOURCE

function gaunt(ll,mm,l1,m1,l2,m2)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: l1,l2,ll,m1,m2,mm
 real(dp) :: gaunt

!Local variables ------------------------------
!scalars
 integer :: i1,i2,j1,j1half,j2,j2half,j3,j3half,j_half,jj,k1,k2,n1,n2
 real(dp) :: argument,sign,sum,xx,yy
 logical :: ok

!************************************************************************

 gaunt=zero;sum=zero;ok =.true.

 if((-m1-mm+m2) /= 0) ok = .false.
 if(abs(m1) > l1) ok = .false.
 if(abs(mm) > ll) ok = .false.
 if(abs(m2) > l2) ok = .false.

 jj = l1 + ll + l2
 if (mod(jj,2)/=0) ok = .false.
 j1 = jj-2*l2
 j2 = jj-2*ll
 j3 = jj-2*l1

 if (j1<0 .or. j2<0 .or. j3<0) ok = .false.

 if (ok) then

   xx = (2 * l1 + 1) * (2 * ll + 1) * (2 * l2 + 1)

   j1half = j1/2
   j2half = j2/2
   j3half = j3/2
   j_half = jj/2

   gaunt = (-1)**j1half * sqrt(xx)
   gaunt = gaunt * rfactorial(j2)*rfactorial(j3)/rfactorial(jj+1)
   gaunt = gaunt * rfactorial(j_half)/(rfactorial(j1half)&
&                * rfactorial(j2half)*rfactorial(j3half))

   yy = rfactorial(l2 + m2) * rfactorial(l2 - m2)

   if (mm>=0) then
     yy = yy * perms(ll+mm,2*mm)
   else
     yy = yy / perms(ll-mm,-2*mm)
   end if

   if (m1>=0) then
     yy = yy / perms(l1+m1,2*m1)
   else
     yy = yy * perms(l1-m1,-2*m1)
   end if

   gaunt = gaunt * sqrt(yy)

   i1 = l2 - ll - m1
   i2 = l2 - l1 + mm
   k1 = -min(0, i1, i2)
   n1 = l1 + m1
   n2 = ll - mm
   k2 = min(j1, n1, n2)

   sign = 1._dp
   if(k1>0) sign = (-1._dp)**k1

   argument = sign     * perms(n1,k1)/rfactorial(k1)
   argument = argument * perms(n2,k1)/rfactorial(i1 + k1)
   argument = argument * perms(j1,k1)/rfactorial(i2 + k1)
   sum = sum + argument

   sign = -sign
   k1 = k1 + 1
   do while(k1 <= k2)
     argument = sign     * perms(n1, k1)/rfactorial(k1)
     argument = argument * perms(n2, k1)/rfactorial(i1 + k1)
     argument = argument * perms(j1, k1)/rfactorial(i2 + k1)
     sum = sum + argument
     sign = -sign
     k1 = k1 + 1
   end do

 end if

 gaunt = gaunt * sum

 end function gaunt
!!***

!----------------------------------------------------------------------

!!****f* m_pawang/gauleg
!! NAME
!! gauleg
!!
!! FUNCTION
!! Compute the coefficients (supports and weights) for Gauss-Legendre integration
!!
!! INPUTS
!!  xmin=lower bound of integration
!!  xmax=upper bound of integration
!!  n=order of integration
!!
!! OUTPUT
!!  x(n)=array of support points
!!  weights(n)=array of integration weights
!!
!! PARENTS
!!      m_pawang
!!
!! CHILDREN
!!
!! SOURCE

 subroutine gauleg(xmin,xmax,x,weights,n)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: n
 real(dp),intent(in) :: xmax,xmin
!arrays
 real(dp),intent(out) :: weights(n),x(n)

!Local variables ------------------------------
!scalars
 integer :: ii,jj
 real(dp),parameter :: tol=1.d-13
 real(dp) :: p1,p2,p3,pi,xl,pp,xmean,z,z1

!************************************************************************

 pi=4._dp*atan(1._dp)
 xl=(xmax-xmin)*0.5_dp
 xmean=(xmax+xmin)*0.5_dp

 do ii=1,(n+1)/2
   z=cos(pi*(ii-0.25_dp)/(n+0.5_dp))
   do
     p1=1._dp
     p2=0._dp
     do jj=1,n
       p3=p2
       p2=p1
       p1=((2._dp*jj-1._dp)*z*p2-(jj-1._dp)*p3)/jj
     end do
     pp=n*(p2-z*p1)/(1._dp-z**2)
     z1=z
     z=z1-p1/pp
     if(abs(z-z1) < tol) exit
   end do
   x(ii)=xmean-xl*z
   x(n+1-ii)=xmean+xl*z
   weights(ii)=2._dp*xl/((1._dp-z**2)*pp**2)
   weights(n+1-ii)=weights(ii)
 end do

 end subroutine gauleg
!!***

!----------------------------------------------------------------------

!!****f* m_pawang/rfactorial
!! NAME
!! rfactorial
!!
!! FUNCTION
!! Private function
!! Calculates N! as a double precision real.
!!
!! INPUTS
!!   nn=number to use
!!
!! OUTPUT
!!   factorial= n! (real)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental function rfactorial(nn)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp) :: rfactorial

!Local variables ---------------------------------------
!scalars
 integer :: ii

! *********************************************************************

 rfactorial=one
 do ii=2,nn
   rfactorial=rfactorial*ii
 end do

end function rfactorial
!!***

!----------------------------------------------------------------------

!!****f* m_pawang/perms
!! NAME
!! perms
!!
!! FUNCTION
!! Private function
!! Returns N!/(N-k)!  if N>=0 and N>k ; otherwise 0 is returned
!!
!! INPUTS
!!   kk=number k to use
!!   nn=number N to use
!!
!! OUTPUT
!!   perms= n!/(n-k)!
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function perms(nn,kk)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: kk,nn
 real(dp) :: perms

!Local variables ---------------------------------------
!scalars
 integer :: ii
 real(dp) :: pp

! *********************************************************************

 if (nn>=0.and.nn>=kk) then
   pp=1._dp
   do ii=nn-kk+1,nn
     pp=pp*ii
   end do
 else
   pp=0._dp
 end if

 perms=pp

end function perms
!!***

!----------------------------------------------------------------------

END MODULE m_pawang
!!***
