!!****m* ABINIT/m_bader
!! NAME
!! m_bader
!!
!! FUNCTION
!! Procedures used by AIM code.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (PCasek,FF,XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_bader

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_sort
 use m_hdr
 use m_splines
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_time,          only : timein
 use m_geometry,      only : metric
 use m_parser,        only : inread
 use m_numeric_tools, only : coeffs_gausslegint
 use m_hide_lapack,   only : jacobi, lubksb, ludcmp

 implicit none

 !private
 public
!!***

!!****t* m_bader/aim_dataset_type
!! NAME
!! aim_dataset_type
!!
!! FUNCTION
!! The aim_dataset_type structured datatype
!! gathers all the input variables for the aim code
!!
!! SOURCE

 type aim_dataset_type

! Since all these input variables are described in the aim_help.html
! file, they are not described in length here ...

! Integer
  integer :: crit
  integer :: denout
  integer :: dltyp
  integer :: gpsurf
  integer :: irho
  integer :: ivol
  integer :: lapout
  integer :: nsa
  integer :: nsb
  integer :: nsc

  integer :: batom  ! Warning : corresponds to the input variable atom
  integer :: foll   ! Warning : corresponds to the input variable follow
  integer :: isurf  ! Warning : corresponds to the input variable surf
  integer :: irsur  ! Warning : corresponds to the input variable rsurf
  integer :: nph    ! Warning : corresponds to the input variable nphi
  integer :: npt    ! Warning : corresponds to the input variable inpt
  integer :: nth    ! Warning : corresponds to the input variable ntheta
  integer :: plden  ! Warning : not documented in help file ?!

  integer :: ngrid(3)

! Real
  real(dp) :: atrad
  real(dp) :: coff1
  real(dp) :: coff2
  real(dp) :: dpclim
  real(dp) :: folstp
  real(dp) :: lgrad
  real(dp) :: lgrad2
  real(dp) :: lstep
  real(dp) :: lstep2
  real(dp) :: maxatd
  real(dp) :: maxcpd
  real(dp) :: phimax
  real(dp) :: phimin

  real(dp) :: dr0    ! Warning : correspond to the input variable radstp
  real(dp) :: phi0   ! Warning : correspond to the input variable rsurdir(2)
  real(dp) :: rmin   ! Warning : correspond to the input variable ratmin
  real(dp) :: th0    ! Warning : correspond to the input variable rsurdir(1)
  real(dp) :: themax ! Warning : correspond to the input variable thetamax
  real(dp) :: themin ! Warning : correspond to the input variable thetamin

  real(dp) :: foldep(3)
  real(dp) :: scal(3)
  real(dp) :: vpts(3,4)

 end type aim_dataset_type
!!***

 public :: adini
 public :: drvaim
 public :: inpar
 public :: defad
 public :: aim_shutdown

 ! Global from defs_aimfields
 integer, save :: ngfft(3),nmax
 integer, allocatable, save :: ndat(:)
 real(dp), allocatable, target, save :: dig1(:),llg1(:),dig2(:),llg2(:),dig3(:),llg3(:)
 real(dp), allocatable, target, save :: cdig1(:),cdig2(:),cdig3(:)
 real(dp), save :: dix(3)
 real(dp), allocatable, target, save :: dvl(:,:,:),ddx(:,:,:),ddy(:,:,:),ddz(:,:,:),rval(:,:,:)
 real(dp), allocatable, save :: rrad(:,:),crho(:,:),sp2(:,:),sp3(:,:),sp4(:,:),pdd(:),pd(:)

 ! Global from defs_aimprom

 ! UNITS
 integer, save :: unt0,unto,unt,untc,unts,untd,untl,untg,unta,untad,untp,untout
 integer,save :: aim_iomode
 ! DRIVER VARIABLES
 real(dp), save :: maxatdst,maxcpdst
 integer, parameter :: ndif=45,ngaus=200,npos=1000
 integer, allocatable, save :: typat(:), corlim(:)
 integer, save :: ntypat,nnpos,natom
 integer, save :: nsimax,batcell,npc,nbcp,nrcp,nccp
 integer, save :: icpc(npos*ndif),npcm3,slc
 real(dp), save :: rprimd(3,3),ivrprim(3,3),trivrp(3,3)
 real(dp), allocatable, save :: xred(:,:),xatm(:,:),rminl(:)
 real(dp), save :: tpi,sqfp,fpi,sqpi,sqtpi,atp(3,npos)
 real(dp), save :: h0,hmin,r0,ttsrf,ttcp,tttot
 real(dp), save :: cth(ngaus),th(ngaus),ph(ngaus),wcth(ngaus),wph(ngaus),rs(ngaus,ngaus)
 real(dp), save :: pc(3,npos*ndif), evpc(3,npos*ndif),zpc(3,3,npos*ndif), pcrb(3,npos*ndif)
 logical, save :: deb,ldeb
!!! interface chgbas
!!!    subroutine bschg1(vv,dir)
!!!      implicit none
!!!      integer, intent(in) :: dir
!!!      real(dp),intent(inout) :: vv(3)
!!!    end subroutine bschg1
!!!    subroutine bschg2(aa,dir)
!!!      implicit none
!!!      integer, intent(in) :: dir
!!!      real(dp),intent(inout) :: aa(3,3)
!!!    end subroutine bschg2
!!! end interface chgbas

!- Set of parameters for the aim utility -----------------------------------
 real(dp), parameter :: aim_rhocormin=1.d-10  ! the minimal core density
 real(dp), parameter :: aim_epstep=0.5
 real(dp), parameter :: aim_rhomin=1.d-5,aim_dgmin=1.d-9,aim_dmaxcrit=5.d-2
 real(dp), parameter :: aim_dmin=1.d-3,aim_hmax=2.d7,aim_fac0=2.1_dp,aim_facmin=1.d-3
 real(dp), parameter :: aim_hmult=15._dp,aim_tiny=1.d-4,aim_snull=1.d-6
 real(dp), parameter :: aim_deltarmin=1.d-7
!the minimal length of one step following the gradient line
 real(dp), parameter :: aim_fac=1.2_dp,aim_drmin=1.d-5
 real(dp), parameter :: aim_dlimit=1.d-4,aim_dmaxcs=3.d-1
 real(dp), parameter :: aim_dpc0=1.d-2
 integer, parameter :: aim_maxstep=100
 real(dp), parameter :: aim_xymin=1.d-10
 integer, parameter :: aim_npmaxin=17
 real(dp), parameter :: aim_stmax=0.05
 real(dp), parameter :: aim_dmaxc1=1.d-1, aim_dmaxcl=5.d-2

!----------------------------------------------------------------------

!!****t* m_bader/bcp_type
!! NAME
!! bcp_type
!!
!! FUNCTION
!! a "bonding critical point" for aim
!!
!! SOURCE

 type, private :: bcp_type

! Integer
  integer :: iat       ! number of the bonding atom inside a primitive cell
  integer :: ipos      ! number of the primitive cell of the bonding atom

! Real
  real(dp) :: chg      ! charge at the critical point
  real(dp) :: diff(3)  ! three distances : AT-CP,BAT-CP,AT-BAT
  real(dp) :: ev(3)    ! eigenvalues of the Hessian
  real(dp) :: pom(3)   ! position of the bonding atom
  real(dp) :: rr(3)    ! position of the bcp
  real(dp) :: vec(3,3) ! eigenvectors of the Hessian
  real(dp) :: vv(3)    ! position of the bcp relative to the central atom

 end type bcp_type
!!***

 contains
!!***

!----------------------------------------------------------------------

!!****f* defs_aimprom/aim_shutdown
!! NAME
!!  aim_shutdown
!!
!! FUNCTION
!!  Free memory allocated in the module. Close units. Mainly used to pass the abirules
!!
!! PARENTS
!!      aim
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

 subroutine aim_shutdown()

!Local variables-------------------------------
 integer :: ii
 logical :: is_open
 integer :: all_units(12)

 ! *********************************************************************

 !if (allocated(typat)) then
 !  ABI_FREE(typat)
 !end if
 !if (allocated(corlim)) then
 !  ABI_FREE(corlim)
 !end if
 !if (allocated(xred)) then
 !  ABI_FREE(xred)
 !end if
 !if (allocated(xatm)) then
 !  ABI_FREE(rminl)
 !end if

 all_units(:) = [unt0,unto,unt,untc,unts,untd,untl,untg,unta,untad,untp,untout]
 do ii=1,size(all_units)
   inquire(unit=all_units(ii), opened=is_open)
   if (is_open) close(all_units(ii))
 end do

end subroutine aim_shutdown
!!***

!!****f* m_bader/adini
!! NAME
!! adini
!!
!! FUNCTION
!! Analysis of the input string "inpstr" (the content of input file)
!! and setting of the corresponding input variables
!!
!! INPUTS
!!  inpstr=character string containing the input data, to be treated
!!  lenstr=actual length of the string contained in inpstr
!!
!! OUTPUT
!!  aim_dtset=the structured entity containing all input variables
!!
!! PARENTS
!!      aim
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine adini(aim_dtset,inpstr,lenstr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lenstr
 character(len=*),intent(in) :: inpstr
!no_abirules
 type(aim_dataset_type), intent(inout) :: aim_dtset !vz_i

!Local variables ------------------------------
!scalars
 integer :: errcod,ii,inxh,ipos,jj,lenc,ll,outi,tstngr=0,tstvpt=0 !vz_z
 real(dp) :: outr
 logical :: nbtst,try
 character(len=20) :: cmot

! *********************************************************************

 if (iachar(inpstr(1:1)) < 32) then
   ipos=2
 else
   ipos=1
 end if

 write(std_out,*) 'ECHO of the INPUT'
 write(std_out,*) '************************'
 write(untout,*) 'ECHO of the INPUT'
 write(untout,*) '************************'

 mread:  do ii=1,lenstr
   try=.false.
   nbtst=.true.
   inxh=index(inpstr(ipos:lenstr),' ')
   if ((ipos >= lenstr)) exit
   if ((inxh==2).or.(inxh==1)) then
     ipos=ipos+inxh
     cycle
   end if
   lenc=inxh-1
   cmot(1:lenc)=inpstr(ipos:ipos+inxh-2)
   ipos=ipos+inxh
!  write(std_out,*) cmot(1:lenc), lenc

   select case (cmot(1:lenc))

!    DRIVER SPECIFICATIONS

   case ('SURF')
     inxh=index(inpstr(ipos:lenstr),' ')
     if ((inxh /= 2).and.(inpstr(ipos:ipos)/='-')) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%isurf=outi
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%isurf
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%isurf
     ipos=ipos+inxh

   case ('CRIT')
     inxh=index(inpstr(ipos:lenstr),' ')
     if ((inxh /= 2).and.(inpstr(ipos:ipos)/='-')) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%crit=outi
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%crit
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%crit
     ipos=ipos+inxh

   case ('RSURF')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%irsur=outi
     write(std_out,*) cmot(1:lenc),'     ', aim_dtset%irsur
     write(untout,*) cmot(1:lenc),'     ', aim_dtset%irsur
     ipos=ipos+inxh

   case ('FOLLOW')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%foll=outi
     write(std_out,*) cmot(1:lenc),'    ', aim_dtset%foll
     write(untout,*) cmot(1:lenc),'    ', aim_dtset%foll
     ipos=ipos+inxh

   case ('IRHO')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%irho=outi
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%irho
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%irho
     ipos=ipos+inxh

   case ('PLDEN')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%plden=outi
     write(std_out,*) cmot(1:lenc),'     ', aim_dtset%plden
     write(untout,*) cmot(1:lenc),'     ', aim_dtset%plden
     ipos=ipos+inxh


   case ('IVOL')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%ivol=outi
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%ivol
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%ivol
     ipos=ipos+inxh

   case ('DENOUT')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%denout=outi
     if ((aim_dtset%denout < -1).or.(aim_dtset%denout>3)) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     write(std_out,*) cmot(1:lenc),'    ', aim_dtset%denout
     write(untout,*) cmot(1:lenc),'    ', aim_dtset%denout
     ipos=ipos+inxh

   case ('LAPOUT')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%lapout=outi
     if ((aim_dtset%lapout < -1).or.(aim_dtset%lapout>3)) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     write(std_out,*) cmot(1:lenc),'    ', aim_dtset%lapout
     write(untout,*) cmot(1:lenc),'    ', aim_dtset%lapout
     ipos=ipos+inxh

   case ('DLTYP')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%dltyp=outi
     write(std_out,*) cmot(1:lenc),'     ', aim_dtset%dltyp
     write(untout,*) cmot(1:lenc),'     ', aim_dtset%dltyp
     ipos=ipos+inxh

   case ('GPSURF')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%gpsurf=outi
     write(std_out,*) cmot(1:lenc),'    ', aim_dtset%gpsurf
     write(untout,*) cmot(1:lenc),'    ', aim_dtset%gpsurf
     ipos=ipos+inxh


!      END OF THE DRIVER SPECIFICATIONS

   case ('ATOM')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%batom=outi
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%batom
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%batom
     ipos=ipos+inxh

   case ('NSA')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%nsa=outi
     write(std_out,*) cmot(1:lenc),'       ', aim_dtset%nsa
     write(untout,*) cmot(1:lenc),'       ', aim_dtset%nsa
     ipos=ipos+inxh

   case ('NSB')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%nsb=outi
     write(std_out,*) cmot(1:lenc),'       ', aim_dtset%nsb
     write(untout,*) cmot(1:lenc),'       ', aim_dtset%nsb
     ipos=ipos+inxh

   case ('NSC')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%nsc=outi
     write(std_out,*) cmot(1:lenc),'       ', aim_dtset%nsc
     write(untout,*) cmot(1:lenc),'       ', aim_dtset%nsc
     ipos=ipos+inxh

   case ('INPT')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%npt=outi
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%npt
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%npt
     ipos=ipos+inxh

   case ('NTHETA')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%nth=outi
     write(std_out,*) cmot(1:lenc),'    ', aim_dtset%nth
     write(untout,*) cmot(1:lenc),'    ', aim_dtset%nth
     ipos=ipos+inxh

   case ('NPHI')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%nph=outi
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%nph
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%nph
     ipos=ipos+inxh

   case ('THETAMIN')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%themin=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'  ', aim_dtset%themin
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'  ', aim_dtset%themin
     ipos=ipos+inxh

   case ('THETAMAX')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%themax=outr
     write(std_out, '(1x,a,a,es17.10)' ) cmot(1:lenc),'  ', aim_dtset%themax
     write(untout,'(1x,a,a,es17.10)') cmot(1:lenc),'  ', aim_dtset%themax
     ipos=ipos+inxh

   case ('PHIMIN')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%phimin=outr
     write(std_out, '(1x,a,a,es17.10)' ) cmot(1:lenc),'    ', aim_dtset%phimin
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%phimin
     ipos=ipos+inxh

   case ('PHIMAX')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%phimax=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%phimax
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%phimax
     ipos=ipos+inxh

   case ('ATRAD')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%atrad=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'     ', aim_dtset%atrad
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'     ', aim_dtset%atrad
     ipos=ipos+inxh

   case ('RADSTP')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%dr0=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%dr0
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%dr0
     ipos=ipos+inxh

   case ('FOLSTP')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%folstp=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%folstp
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%folstp
     ipos=ipos+inxh

   case ('RATMIN')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%rmin=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%rmin
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%rmin
     ipos=ipos+inxh

   case ('COFF1')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%coff1=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%coff1
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%coff1
     ipos=ipos+inxh

   case ('COFF2')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%coff2=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%coff2
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%coff2
     ipos=ipos+inxh

   case ('DPCLIM')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%dpclim=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%dpclim
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%dpclim
     ipos=ipos+inxh

   case ('LGRAD')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%lgrad=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lgrad
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lgrad
     ipos=ipos+inxh

   case ('LGRAD2')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%lgrad2=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lgrad2
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lgrad2
     ipos=ipos+inxh

   case ('LSTEP')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%lstep=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lstep
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lstep
     ipos=ipos+inxh

   case ('LSTEP2')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%lstep2=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lstep2
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lstep2
     ipos=ipos+inxh

   case ('RSURDIR')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%th0=outr
     ipos=ipos+inxh
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%phi0=outr
     ipos=ipos+inxh
     write(std_out, '(1x,a,a,2es17.10)') cmot(1:lenc),'   ', aim_dtset%th0, aim_dtset%phi0
     write(untout, '(1x,a,a,2es17.10)') cmot(1:lenc),'   ', aim_dtset%th0, aim_dtset%phi0

   case ('FOLDEP')
     do jj=1,3
       inxh=index(inpstr(ipos:lenstr),' ')
       call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
       aim_dtset%foldep(jj)=outr
       ipos=ipos+inxh
     end do
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%foldep
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%foldep

   case ('SCAL')
     do jj=1,3
       inxh=index(inpstr(ipos:lenstr),' ')
       call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
       aim_dtset%scal(jj)=outr
       ipos=ipos+inxh
     end do
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%scal
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%scal

   case ('NGRID')
     try=.true.
     do jj=1,3
       inxh=index(inpstr(ipos:lenstr),' ')
       call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
       aim_dtset%ngrid(jj)=outi
       if (.not.nbtst) then
         tstngr=jj-1
         cycle mread
       end if
       if (inxh==0) then
         tstvpt=jj
         exit mread
       end if
       ipos=ipos+inxh
       if (ipos==lenstr-1) then
         tstngr=jj
         exit mread
       end if
     end do
!      Why no echo ?? XG 030218
     tstngr=3

   case ('VPTS')
     do jj=1,4
       do ll=1,3
         try=.true.
         inxh=index(inpstr(ipos:lenstr),' ')
         call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
         aim_dtset%vpts(ll,jj)=outr
         if (.not.nbtst) then
           tstvpt=jj-1
           cycle mread
         end if
         ipos=ipos+inxh
         if (ipos>=lenstr) then
           tstvpt=jj
           exit mread
         end if
       end do
     end do
!      Why no echo ?? XG 030218
     tstvpt=4

   case ('MAXATD')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%maxatd=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%maxatd
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%maxatd
     ipos=ipos+inxh

   case ('MAXCPD')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%maxcpd=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%maxcpd
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%maxcpd
     ipos=ipos+inxh

   case default
     write(std_out,*) 'ERROR Bad key word ! ',cmot(1:lenc)
   end select
 end do mread

 write(std_out,*) '************************'

 call consist(aim_dtset,tstngr,tstvpt)

end subroutine adini
!!***

!!****f* m_bader/addout
!! NAME
!! addout
!!
!! FUNCTION
!! Output density and laplacian (see input variables denout and lapout)
!!
!! INPUTS
!!  aim_dtset=the structured entity containing all input variables
!!  also, uses the variables saved in the module "defs_aimprom"
!!
!! OUTPUT
!!  (print)
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet) : the use
!! of a module to transfer data should be avoided
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine addout(aim_dtset)

!Arguments ------------------------------------
!scalars
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: cod,dims,iat,ii,ipos,jj,nn,tgrd
 real(dp) :: alfa,rho,rr,xx,yy
!arrays
 real(dp) :: grho(3),hrho(3,3),orig(3),vv(3)
 real(dp),allocatable :: dfld(:),lfld(:),nr(:),stp(:),uu(:,:)

!************************************************************************
 orig(:)=aim_dtset%vpts(:,1)
 if (aim_dtset%denout > 0) then
   dims=aim_dtset%denout
 elseif (aim_dtset%lapout > 0) then
   dims=aim_dtset%lapout
 end if

 select case (aim_dtset%dltyp)
 case (1)
   cod=1
 case (2)
   cod=2
 case default
   cod=0
 end select

 ABI_ALLOCATE(uu,(3,dims))
 ABI_ALLOCATE(nr,(dims))
 ABI_ALLOCATE(stp,(dims))

 write(std_out,*) 'grid:', aim_dtset%ngrid(1:dims)
 write(std_out,*) 'kod :', cod
 tgrd=1
 do ii=1,dims
   tgrd=tgrd*aim_dtset%ngrid(ii)
   uu(:,ii)=aim_dtset%vpts(:,ii+1)-aim_dtset%vpts(:,1)
   nr(ii)=vnorm(uu(:,ii),0)
   stp(ii)=nr(ii)/(aim_dtset%ngrid(ii)-1)
   uu(:,ii)=uu(:,ii)/nr(ii)
 end do
 write(std_out,*) 'tgrd :', tgrd
 do ii=1,dims
   write(std_out,*) 'uu :', uu(1:3,ii)
 end do

 if (aim_dtset%denout > 0) then
   ABI_ALLOCATE(dfld,(tgrd+1))
   dfld(:)=0._dp
 end if
 if (aim_dtset%lapout > 0)  then
   ABI_ALLOCATE(lfld,(tgrd+1))
 end if

 select case (dims)
 case (1)
   nn=0
   do ii=0,aim_dtset%ngrid(1)-1
     nn=nn+1
     vv(:)=orig(:)+ii*stp(1)*uu(:,1)
     call vgh_rho(vv,rho,grho,hrho,rr,iat,ipos,cod)
     if (aim_dtset%denout > 0) dfld(nn)=rho
     if (aim_dtset%lapout > 0) lfld(nn)=hrho(1,1)+hrho(2,2)+hrho(3,3)
   end do
   if (aim_dtset%denout==1) then
     do ii=0,aim_dtset%ngrid(1)-1
       xx=ii*stp(1)
       write(untd,'(2E16.8)') xx, dfld(ii+1)
     end do
   end if
   if (aim_dtset%lapout==1) then
     do ii=0,aim_dtset%ngrid(1)-1
       xx=ii*stp(1)
       write(untl,'(2E16.8)') xx, lfld(ii+1)
     end do
   end if
 case (2)
   nn=0
   alfa=dot_product(uu(:,1),uu(:,2))
   alfa=acos(alfa)
   do ii=0,aim_dtset%ngrid(2)-1
     do jj=0,aim_dtset%ngrid(1)-1
       nn=nn+1
       vv(:)=orig(:)+jj*uu(:,2)*stp(2)+ii*stp(1)*uu(:,1)
       call vgh_rho(vv,rho,grho,hrho,rr,iat,ipos,cod)
       if (aim_dtset%denout > 0) dfld(nn)=rho
       if (aim_dtset%lapout > 0) lfld(nn)=hrho(1,1)+hrho(2,2)+hrho(3,3)
     end do
   end do
   write(std_out,*) 'generace hotova', nn
   nn=0
   if (aim_dtset%denout==2) then
     do ii=0,aim_dtset%ngrid(2)-1
       do jj=0,aim_dtset%ngrid(1)-1
         nn=nn+1
         xx=jj*stp(1)+cos(alfa)*ii*stp(2)
         yy=sin(alfa)*ii*stp(2)
         write(untd,'(3E16.8)') xx, yy, dfld(nn)
       end do
       write(untd,*) ' '
     end do
   end if
   nn=0
   if (aim_dtset%lapout==2) then
     write(std_out,*) 'lezes sem?'
     do ii=0,aim_dtset%ngrid(2)-1
       do jj=0,aim_dtset%ngrid(1)-1
         nn=nn+1
         xx=jj*stp(1)+cos(alfa)*ii*stp(2)
         yy=sin(alfa)*ii*stp(2)
         write(untl,'(3E16.8)') xx, yy, lfld(nn)
       end do
       write(untl,*) ' '
     end do
   end if
 end select
 ABI_DEALLOCATE(uu)
 ABI_DEALLOCATE(stp)
 ABI_DEALLOCATE(nr)
 if(aim_dtset%denout>0) then
   ABI_DEALLOCATE(dfld)
 end if
 if(aim_dtset%lapout>0) then
   ABI_DEALLOCATE(lfld)
 end if

end subroutine addout
!!***

!!****f* m_bader/aim_follow
!! NAME
!! aim_follow
!!
!! FUNCTION
!! This routine follows the gradient line starting from the point
!! vv. It stop when it arrives to the atom (nearer than rminl(iat))
!! or - if srch=true - also if it arrives under the already known
!! part of Bader surface
!!
!! INPUTS
!! aim_dtset= the structured entity containing all input variables
!! iatinit,iposinit= indexes of initial atom
!! npmax= maximum number of division in each step
!!
!! OUTPUT
!! iat,ipos= index of final atom
!! nstep= returns the number of step needed
!!
!! SIDE EFFECTS
!! srch=  (true/false) check if the line is outside or
!!             inside the atomic surface.
!! vv(3)= initial point in orthogonal coordinates
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine aim_follow(aim_dtset,vv,npmax,srch,iatinit,iposinit,iat,ipos,nstep)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iatinit,iposinit,npmax
 integer,intent(out) :: iat,ipos,nstep
 logical,intent(inout) :: srch
 type(aim_dataset_type),intent(in) :: aim_dtset
!arrays
 real(dp),intent(inout) :: vv(3)

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3,ii,iph,ires,ith,jj,kk,nit,np,nph,nsi,nth
 real(dp) :: deltar,dg,dist,dph,dth,fac2,facf,h0old,hh,hold,rho,rr,rsmed
 real(dp) :: t1,t2,t3,vcth,vph,vth,wall,xy,xyz
 logical :: fin,ldebold,srchold,stemp,stemp2
 character(len=50) :: formpc
 character(len=500) :: msg
!arrays
 real(dp) :: ev(3),grho(3),hrho(3,3),pom(3),vold(3),vt(3),vt1(3)
 real(dp) :: zz(3,3)

!************************************************************************
 formpc='(":CP",2I5,3F12.8,3E12.4,I4,2E12.4)'


 fin=.false.

 srchold=srch
 ldebold=ldeb
 h0old=h0

 nth=aim_dtset%nth
 nph=aim_dtset%nph

 if (slc==0) then
   rminl(:)=aim_dtset%rmin
 end if

 if (deb) then
   ldeb=.true.
 end if

 call vgh_rho(vv,rho,grho,hrho,rr,iat,ipos,slc)

!Initial tests

 if (iat/=0) then
   if (rr<rminl(iat)) then
     fin=.true.
     write(std_out,*) 'rr < rmin iat=',iat,' ipos=',ipos
   elseif (rho<aim_rhomin) then
     fin=.true.
     write(std_out,*) 'CHARGE LT rhomin ',rho,' < ',aim_rhomin
     if (rho<zero) then
       MSG_ERROR('RHO < 0 !!!')
     end if
   end if
 end if

 facf=aim_fac0
 hh=aim_hmax

 call timein(t1,wall)
 nstep=0
 nsi=0

!the principal cycle

 madw : do while(.not.fin)
   hold=hh

   dg=vnorm(grho,0)
   if (ldeb.or.deb) write(std_out,*) 'dg= ',dg

!  the time test

   call timein(t3,wall)
   t2=t3-t1
   if (t2>300.0) then
     write(std_out,*) 'TIME EXCEEDED 5 min IN FOLLOW'
     write(std_out,*) 'h0 =',h0,'  h =',hh,'  h0old =',h0old,'  dg =',dg
     write(std_out,*) 'facf =',facf
     msg =  'TIME EXCEEDED 5 min IN FOLLOW'
     MSG_ERROR(msg)
   end if

   if (dg<aim_dgmin) then
     write(std_out,*) 'gradient < dgmin ',dg,' < ',aim_dgmin
     fin=.true.
     iat=0
     ipos=0
!    testing for the CP
     if (npc>0) then
       call critic(aim_dtset,vv,ev,zz,aim_dmaxcrit,ires,0)
       if (ires==0) then
         do jj=1,npc
           pom(:)=pc(:,jj)-vv(:)+xatm(:,aim_dtset%batom)
           dist=vnorm(pom,0)
           if (dist<aim_tiny) cycle madw
         end do
         write(std_out,*) 'C.P. found !!'
         npc=npc+1
         do jj=1,3
           pc(jj,npc)=vv(jj)
           evpc(jj,npc)=ev(jj)
           do kk=1,3
             zpc(kk,jj,npc)=zz(kk,jj)
           end do
         end do
         i1=ev(1)/abs(ev(1))
         i2=ev(2)/abs(ev(2))
         i3=ev(3)/abs(ev(3))
         icpc(npc)=i1+i2+i3
         if (icpc(npc)==-3) then           ! pseudoatom handling
           npcm3=npcm3+1
           write(std_out,*) 'Pseudo-atom found !!'
         end if

         call vgh_rho(vv,rho,grho,hrho,rr,iat,ipos,slc)
         write(22,formpc) 0,0,(pcrb(jj,npc),jj=1,3),(ev(jj),jj=1,3),icpc(npc),&
&         ev(1)+ev(2)+ev(3),rho
         write(std_out,formpc) 0,0,(pcrb(jj,npc),jj=1,3),(ev(jj),jj=1,3),icpc(npc),&
&         ev(1)+ev(2)+ev(3),rho
       else
         write(std_out,*) 'C.P. not found !!'
       end if
     end if

     cycle madw
   end if

   hh=h0/dg
   if (ldeb.or.deb) write(std_out,*) 'h= ',hh,' h0= ',h0,' dg= ',dg
   if (hh>aim_hmax) hh=aim_hmax
!  step modifications

   hh=hh*facf
   if (hh>(hold*aim_hmult)) then
     hh=hold*aim_hmult
   end if

   do ii=1,3
     vold(ii)=vv(ii)
   end do

   nit=0
   hold=hh

!  one step following the gradient line
!
   call onestep(vv,rho,grho,hh,np,npmax,deltar)
   do while (((np>npmax).or.(deltar>aim_stmax)).and.(deltar>aim_dmin))
     nit=nit+1
     if (nit>5) then
       if (deltar>aim_stmax) then
         write(std_out,*) 'nit > 5 and deltar > stmax   nit=',nit
       else
         write(std_out,*) 'nit > 5 and np > npmax   nit=',nit
       end if
     end if
     do ii=1,3
       vv(ii)=vold(ii)
     end do
     hh=hh*0.3
     call onestep(vv,rho,grho,hh,np,npmax,deltar)
   end do


   nstep=nstep+1
   if (ldeb.or.deb) write(std_out,*) 'h= ',hh

   fac2=hh/hold
   if (fac2>=1._dp) then
     facf=facf*1.2
   else
     if (fac2>=aim_facmin) then
       facf=fac2
     else
       facf=aim_facmin
     end if
   end if

   if (deb.or.ldeb) then
     write(std_out,*) ':POS ',vv
     write(std_out,*) ':RBPOS ',vt1
     write(std_out,*) ':GRAD ',grho
   end if

   call vgh_rho(vv,rho,grho,hrho,rr,iat,ipos,slc)
   dg=vnorm(grho,0)
   pom(:)=vv(:)-xatm(:,iatinit)-atp(:,iposinit)

   if (iat /= 0) then
     fin=.true.
     write(std_out,*) 'r < rmin iat=',iat,' ipos=',ipos
     cycle madw
   end if

   if (rho<aim_rhomin) then
     fin=.true.
     write(std_out,*) 'charge < rhomin ',rho,' < ',aim_rhomin
     if (rho<zero) then
       MSG_ERROR('RHO < 0 !!!')
     end if
     iat=0
     ipos=0
     cycle madw
   end if

   if (npcm3>0) then
     do jj=1,npc
       if (icpc(jj)==(-3)) then
         pom(:)=pc(:,jj)-vv(:)+xatm(:,aim_dtset%batom)
         dist=vnorm(pom,0)
         if (dist<(aim_dtset%rmin**2*0.1)) then
           iat=0
           ipos=0
           fin=.true.
           write(std_out,*) 'We are inside a pseudo-atom'
           cycle madw
         end if
       end if
     end do
   end if

   nsi=nsi+1

!  surface checking

   if (srch.and.(nsi>=nsimax)) then
     nsi=0
     ith=0
     iph=0
     do ii=1,3
       vt(ii)=vv(ii)-xatm(ii,iatinit)
     end do
     xy=vt(1)*vt(1)+vt(2)*vt(2)
     xyz=xy+vt(3)*vt(3)
     xyz=sqrt(xyz)
     if (xy<aim_snull) then
       vcth=1._dp
       if (vt(3)<0._dp) vcth=-vcth
       vph=0._dp
     else
       vcth=vt(3)/xyz
       vph=atan2(vt(2),vt(1))
     end if
     vth=acos(vcth)
     if (vth<th(1)) then
       ith=0
     else
       if (vth>th(nth)) then
         ith=nth
       else
         do ii=2,nth
           if (vth<th(ii)) then
             ith=ii-1
             exit
           end if
         end do
       end if
     end if

     if (vph<ph(1)) then
       iph=0
     else
       if (vph>ph(nph)) then
         iph=nph
       else
         do ii=2,nph
           if (vph<ph(ii)) then
             iph=ii-1
             exit
           end if
         end do
       end if
     end if

     stemp=(iph>0).and.(iph<nph)
     stemp=stemp.and.(ith>0).and.(ith<nth)

     if (stemp) then
       stemp2=rs(ith,iph)>0._dp
       stemp2=stemp2.and.(rs(ith+1,iph)>0._dp)
       stemp2=stemp2.and.(rs(ith+1,iph+1)>0._dp)
       stemp2=stemp2.and.(rs(ith,iph+1)>0._dp)
       if (stemp2) then
         dth=th(ith+1)-th(ith)
         dph=ph(iph+1)-ph(iph)
         rsmed=rs(ith,iph)*(th(ith+1)-vth)/dth*(ph(iph+1)-vph)/dph
         rsmed=rsmed+rs(ith+1,iph)*(vth-th(ith))/dth*(ph(iph+1)-vph)/dph
         rsmed=rsmed+rs(ith+1,iph+1)*(vth-th(ith))/dth*(vph-ph(iph))/dph
         rsmed=rsmed+rs(ith,iph+1)*(th(ith+1)-vth)/dth*(vph-ph(iph))/dph
         if (rsmed>xyz) then
           write(std_out,*) 'We are inside the surface'
           iat=iatinit
           ipos=iposinit
         else
           write(std_out,*) 'We are outside the surface'
           iat=0
           ipos=0
         end if
         fin=.true.
         cycle madw
       end if
     end if
   end if

 end do madw


 srch=srchold
 ldeb=ldebold
 h0=h0old


end subroutine aim_follow
!!***

!!****f* m_bader/consist
!! NAME
!! consist
!!
!! FUNCTION
!! Checking of the consistency between the values of input variables
!!
!! INPUTS
!!  aim_dtset= the structured entity containing all input variables
!!  tstngr= information about the test on the ngrid input variable
!!  tstvpt= information about the test on the vpts input variable
!!
!! OUTPUT
!!  (only checking : print error message and stop if there is a problem)
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine consist(aim_dtset,tstngr,tstvpt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: tstngr,tstvpt
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------

! *********************************************************************

!write(std_out,*) tstngr, tstvpt

 if (((aim_dtset%denout/=0).or.(aim_dtset%lapout/=0)).and.((tstngr < 1).or.(tstvpt < 2))) then
   MSG_ERROR('in input1 - I cannot do the output !')
 end if
 if ((aim_dtset%denout > 0).and.(aim_dtset%lapout>0)) then
   if (aim_dtset%denout/=aim_dtset%lapout) then
     write(std_out,*) 'ERROR in input - when both denout and lapout are positive non-zero,'
     write(std_out,*) 'they must be equal.'
     MSG_ERROR("Aborting now")
   end if
   if ((tstvpt < aim_dtset%denout+1).or.(tstngr < aim_dtset%denout)) then
     write(std_out,*) 'ERROR in input2 - I cannot do the output !'
     MSG_ERROR("Aborting now")
   end if
 elseif (aim_dtset%denout > 0) then
   if ((tstvpt < aim_dtset%denout+1).or.(tstngr < aim_dtset%denout)) then
     write(std_out,*) 'ERROR in input - I cannot do the output !'
     MSG_ERROR("Aborting now")
   end if
 elseif (aim_dtset%lapout > 0) then
   if ((tstvpt < aim_dtset%lapout+1).or.(tstngr < aim_dtset%lapout)) then
     write(std_out,*) 'ERROR in input - I cannot do the output !'
     MSG_ERROR("Aborting now")
   end if
 end if

 if ((aim_dtset%isurf==1).and.(aim_dtset%crit==0)) then
   write(std_out,*) 'ERROR in input - must have crit/=0 for isurf==1'
   MSG_ERROR("Aborting now")
 end if

 if (((aim_dtset%ivol/=0).or.(aim_dtset%irho/=0)).and.(aim_dtset%isurf==0)) then
   MSG_ERROR('in input - I cannot integrate without surface !')
 end if

end subroutine consist
!!***

!!****f* m_bader/cpdrv
!! NAME
!! cpdrv
!!
!! FUNCTION
!! Critical points (CPs) searching driver
!! First Bond CPs are searched for each pair atom-its neighbor
!! (distance cutoff=maxatdst)
!! then Ring CPs for each pair of BCPs
!! and finally Cage CPs for each pair of RCPs.
!!
!! INPUTS
!! aim_dtset= the structured entity containing all input variables
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  this routine treat information contained in the aim_prom module
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! TODO
!! Should combine parts of code that are similar ...
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine cpdrv(aim_dtset)

!Arguments ------------------------------------
!scalars
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: iat,iatinit,ii,inxat,inxcell,ipair,ipos,iposinit,ires,jj,kk,nb,nb_now
 integer :: nn,nstep,nvs,me,nproc,ierr
 real(dp) :: candidate,chg,diff1,diff2,diff3,dist,prj,rtdiff,ss,tt0,wall
 logical :: srch=.false.
!arrays
 integer :: ibat(nnpos*natom),inatm(nnpos*natom),incell(nnpos*natom)
 integer :: ipibat(nnpos*natom)
 integer,allocatable :: indexcp(:),nr(:)
 real(dp) :: bmin(natom),dif(3),dists(nnpos*natom),ev(3),evec(3,3),grho(3)
 real(dp) :: hrho(3,3),pom(3),rr(3),uu(3),vv(3),xorig(3)
 real(dp),allocatable :: buffer(:,:),sortguide(:)
!no_abirules
!Warning : bcp_type should be transformed to cp_type
 type(bcp_type),allocatable :: bcp(:),ccp(:),cp_tmp(:),rcp(:)

!************************************************************************

 me=xmpi_comm_rank(xmpi_world)
 nproc=xmpi_comm_size(xmpi_world)

!Consider the critical points starting from atom #batom
 inxat=aim_dtset%batom
 slc=-1
 rminl(:)=aim_dtset%rmin
 bmin(:)=0._dp
 ttcp=0._dp

 write(std_out,*)
 write(std_out,*) "CRITICAL POINTS ANALYSIS"
 write(std_out,*) "========================"
 write(std_out,*)

 write(untout,*)
 write(untout,*) "CRITICAL POINTS ANALYSIS"
 write(untout,*) "========================"
 write(untout,*)


 xorig(:)=xatm(:,inxat)

 call timein(tt0,wall)

!Searching the neighbouring atoms

 if (aim_dtset%crit > 0) then
   nvs=0
   do ii=1,nnpos
     do jj=1,natom
       dist=0._dp
       dif(:)=xatm(:,inxat)-xatm(:,jj)-atp(:,ii)
       dif(:)=dif(:)/aim_dtset%scal(:)
       dist=vnorm(dif,0)
       if (dist < tol6 ) then
         inxcell=ii
       elseif (dist < maxatdst) then
         nvs=nvs+1
         dists(nvs)=dist
         inatm(nvs)=jj
         incell(nvs)=ii
       end if
     end do
   end do

   write(std_out,*) "ATOM:"
   write(std_out,*) 'inxat :', inxat, 'inxcell :', inxcell
   write(std_out, '(3es16.6)' ) (xorig(ii),ii=1,3)
   write(std_out,*)

   write(untout,*) "ATOM:"
   write(untout,*) 'inxat :', inxat, 'inxcell :', inxcell
   write(untout, '(3es16.6)') (xorig(ii),ii=1,3)
   write(untout,*)

   ABI_ALLOCATE(nr,(nvs))
   do ii=1,nvs
     nr(ii)=ii
   end do

!  Ordering of the nearest neighbouring atoms
   call sort_dp(nvs,dists,nr,tol14)

   nb=0
   write(std_out,*) "NEIGHBORING ATOMS (atindex,cellindex,distance(in bohr)):"
   write(untout,*) "NEIGHBORING ATOMS (atindex,cellindex,distance(in bohr)):"
   do ii=1,nvs
     nn=nr(ii)
     if (dists(ii) < maxatdst) then
       nb=nb+1
       ibat(nb)=inatm(nn)
       ipibat(nb)=incell(nn)
       write(std_out,*) ':NEIG ',inatm(nn),incell(nn),dists(ii)
       write(untout,'("      ",2I6,F16.8)')inatm(nn),incell(nn),dists(ii)
     else
       exit
     end if
   end do

!  SEARCHING BCP
   ABI_DATATYPE_ALLOCATE(bcp,(nb))
   nbcp=0
   iatinit=inxat
   iposinit=inxcell
   bcp(:)%iat=0
   bcp(:)%ipos=0

   write(std_out,*)
   write(std_out,*) "BONDING CRITICAL POINTS (BCP)"
   write(std_out,*) "============================="
   write(std_out,*)

   write(untout,*)
   write(untout,*) "BONDING CRITICAL POINTS (BCP)"
   write(untout,*) "============================="
   write(untout,*)

   srbcp: do ii=1,nb

!    Start the search for BCP from the midistance between the atom
!    and his neighbor.
     vv(:)=(xatm(:,inxat)+xatm(:,ibat(ii))+atp(:,ipibat(ii)))/2._dp

     call critic(aim_dtset,vv,ev,evec,aim_dmaxcs,ires,-1)

     if (ires==0) then
!      Testing if CP is already known
       if (nbcp > 0) then
         do jj=1,nbcp
           pom(:)=vv(:)-bcp(jj)%rr(:)-xorig(:)
           dist=vnorm(pom,0)
           if (dist < aim_dtset%dpclim) then
             write(std_out,*) 'BCP already known  !'
             cycle srbcp
           end if
         end do
       end if
       rr(:)=vv(:)-xorig(:)
       ss=vnorm(rr,0)
       if (ss > maxcpdst) then
         write(std_out, '(a,es16.6,a,es16.6)' ) 'BCP distance from atom,',ss,', exceed maxcpdst =',maxcpdst
         cycle srbcp
       end if
       nn=0
       do jj=1,3
         nn=nn+ev(jj)/abs(ev(jj))
       end do
       write(std_out, '(a,3es16.6,i4)') ' vv(1:3), nn',(vv(jj), jj=1,3), nn
       write(std_out, '(a,3es16.6)') 'ev: ', (ev(jj), jj=1,3)
       if (nn /= -1) then
         write(std_out,*) ' The trial critical point is not a BCP !'
         cycle srbcp
       end if
       write(std_out, '(a,3es16.6)' ) 'evec(:,1): ',(evec(jj,1), jj=1,3)
       pom(:)=evec(:,1)
       dist=vnorm(pom,0)
       prj=dot_product(evec(:,1),rr)
       write(std_out,*) 'prj:', prj, vnorm(evec(:,1),0)
       dist=vnorm(evec(:,1),0)
       uu(:)=vv(:)-sign(aim_epstep,prj)*evec(:,1)/dist

!      Testing whether this BCP "is bonded" to the considered atom
       call aim_follow(aim_dtset,uu,aim_npmaxin,srch,iatinit,iposinit,iat,ipos,nstep)
!      write(std_out,*) 'do', iat, ipos
!      if ((iat==0).or.(ipos==0)) cycle
!      write(std_out,*) 'APOS: ',(xatm(jj,iat)+atp(jj,ipos), jj=1,3)
       if ((iat/=inxat).or.(inxcell/=ipos)) then
         write(std_out,*) ' The trial BCP is not bonded to the Bader atom'
         cycle srbcp
       end if

!      A new BCP has been found !
       nbcp=nbcp+1

!      Searching for the second bonded atom
       ss=vnorm(rr,0)
       diff1=ss
       diff3=dists(ii)
       uu(:)=vv(:)+sign(aim_epstep,prj)*evec(:,1)/dist
       if ((abs(bmin(iat))<1.0d-12).or.( ss<bmin(iat))) then
         bmin(iat)=ss
       end if
       call aim_follow(aim_dtset,uu,aim_npmaxin,srch,iatinit,iposinit,iat,ipos,nstep)
       if ((iat==0).or.(ipos==0)) then
         write(std_out,*) ' The trial BCP is not bonded to a bonding atom !'
!        cycle srbcp
       end if
       pom(:)=vv(:)-xatm(:,iat)-atp(:,ipos)
       ss=vnorm(pom,0)
       diff2=ss
       pom(:)=xorig(:)-xatm(:,iat)-atp(:,ipos)
       diff3=vnorm(pom,0)
       rtdiff=diff1/diff3
       if ((abs(bmin(iat))<1.0d-12).or.(ss<bmin(iat))) then
         bmin(iat)=ss
       end if
       pom(:)=xatm(:,iat)+atp(:,ipos)

!      Store more results, for coherent, and portable output
       bcp(nbcp)%iat=iat
       bcp(nbcp)%ipos=ipos
       bcp(nbcp)%chg=chg
       bcp(nbcp)%diff(1)=diff1
       bcp(nbcp)%diff(2)=diff2
       bcp(nbcp)%diff(3)=diff3
       bcp(nbcp)%ev(:)=ev(:)
       bcp(nbcp)%pom(:)=pom(:)
       bcp(nbcp)%rr(:)=rr(:)
       bcp(nbcp)%vec(:,:)=evec(:,:)
       bcp(nbcp)%vv(:)=vv(:)
!      Warning : iat, ipos might be modified by this call
       call vgh_rho(vv,chg,grho,hrho,dist,iat,ipos,0)
       bcp(nbcp)%chg=chg

     end if ! ires==0
   end do srbcp

   if(nbcp>0)then

!    Order the BCP. CPs should appear by increasing values of x,y,z , the latter
!    varying the fastest
     ABI_ALLOCATE(sortguide,(nbcp))
     ABI_ALLOCATE(indexcp,(nbcp))
     ABI_DATATYPE_ALLOCATE(cp_tmp,(nbcp))
     do ii=3,1,-1
!      DEBUG
!      write(std_out,*)' cpdrv : sort on index ii=',ii
!      ENDDEBUG

       do jj=1,nbcp
!        DEBUG
!        write(std_out,*)bcp(jj)%vv(:)
!        ENDDEBUG
         sortguide(jj)=bcp(jj)%vv(ii)
         indexcp(jj)=jj
       end do

!      Try to be platform-independent. Might need a larger tolerance.
       call sort_dp(nbcp,sortguide,indexcp,tol3)
       do jj=1,nbcp
         cp_tmp(jj)=bcp(indexcp(jj))
       end do
       do jj=1,nbcp
         bcp(jj)=cp_tmp(jj)
       end do
     end do
!    DEBUG
!    write(std_out,*)' cpdrv : after the sort '
!    do jj=1,nbcp
!    write(std_out,*)bcp(jj)%vv(:)
!    end do
!    ENDDEBUG


!    Output the info about the BCP
     do jj=1,nbcp
       write(untout,'(" Bonded atom (BAT) (indxatm,indxcell,position): ",/,2I6,3F16.8)')&
&       bcp(jj)%iat,bcp(jj)%ipos,bcp(jj)%pom(:)
       write(untout,'("%Bonding CP: ",3F16.8)') bcp(jj)%vv(:)
       write(untout,'("%Eigenval. of Hessian: ",3F16.8)') bcp(jj)%ev(:)
       write(untout,'(a,a,a,3f16.8,a,a,3f16.8,a,a,3f16.8,a)') &
&       ' Eigenvec. of Hessian:',char(10),&
&       '-',bcp(jj)%vec(1,:),char(10),&
&       '-',bcp(jj)%vec(2,:),char(10),&
&       '-',bcp(jj)%vec(3,:),char(10)
       write(untout,'("%Density and laplacian in CP: ",2F16.8)') &
&       bcp(jj)%chg, bcp(jj)%ev(1)+bcp(jj)%ev(2)+bcp(jj)%ev(3)
       write(untout,'("%Relative position of BCP (AT-CP,BAT-CP,AT-BAT,relative(AT): ",/,4F16.8)') &
&       bcp(jj)%diff(:),bcp(jj)%diff(1)/bcp(jj)%diff(3)
       write(untout,*) "********************************************************************"
       write(std_out,'(/," BCP: ",3F10.6,3E12.4,E12.4,/)') &
&       bcp(jj)%rr(:),bcp(jj)%ev(:),bcp(jj)%ev(1)+bcp(jj)%ev(2)+bcp(jj)%ev(3)
       write(std_out,'(":DISPC ",4F12.6)') bcp(jj)%diff(:),bcp(jj)%diff(1)/bcp(jj)%diff(3)
     end do

     ABI_DATATYPE_DEALLOCATE(cp_tmp)
     ABI_DEALLOCATE(indexcp)
     ABI_DEALLOCATE(sortguide)

   end if ! nbcp>0

   if (abs(bmin(inxat))>1.0d-12) then
     rminl(inxat)=aim_dtset%coff1*bmin(inxat)
     r0=bmin(inxat)
   else
     r0=0._dp
   end if

!  !AD-HOC PARAMETER

   do ii=1,natom
     if ((abs(bmin(ii))>1.0d-12).and.(ii /= inxat)) rminl(ii)=aim_dtset%coff2*bmin(ii)
   end do

!  END WARNING

   write(std_out,*) ' number of BCP:', nbcp
   write(untout,'(" Number of BCP found: ",I4)') nbcp
   nn=nbcp*(nbcp-1)*(nbcp-2)/6
   if (bit_size(ii) <= nbcp+1) then
     MSG_ERROR("b-test!")
   end if

!  SEARCHING RCP

   write(std_out,*)
   write(std_out,*) "RING CRITICAL POINTS (RCP)"
   write(std_out,*) "============================="
   write(std_out,*)

   write(untout,*)
   write(untout,*) "RING CRITICAL POINTS (RCP)"
   write(untout,*) "============================="
   write(untout,*)

   nrcp=0
   if(aim_dtset%crit==1)nb_now=nbcp
   if(aim_dtset%crit==2)nb_now=nb
!  DEBUG
!  nb_now=nbcp
!  ENDDEBUG
   nn=nb_now*(nb_now-1)/2
   ABI_DATATYPE_ALLOCATE(rcp,(nn))

!  Loop on pairs of BCP or atoms
   ipair=0
   ABI_ALLOCATE(buffer,(16,nn))
   buffer=zero

!  DEBUG
!  write(std_out,*)ch10,ch10,' drvcpr : enter loop to search for RCPs,nb_now,nn=',nb_now,nn
!  ENDDEBUG

   do ii=1,nb_now-1
     srcp1: do jj=ii+1,nb_now
       ipair=ipair+1
       if(mod(ipair,nproc)==me)then
         if (aim_dtset%crit==1) then
           vv(:)=xorig(:)+(bcp(ii)%rr(:)+bcp(jj)%rr(:))/2._dp
         else if (aim_dtset%crit==2) then
           vv(:)=xorig(:)*half+(xatm(:,ibat(ii))+atp(:,ipibat(ii))+xatm(:,ibat(jj))+atp(:,ipibat(jj)))*quarter
         end if

         call critic(aim_dtset,vv,ev,evec,aim_dmaxcs,ires,1)

         if(ires==1)then
           cycle srcp1
         end if

!        Check that it is within the maximum allowed distance for a CP
         rr(:)=vv(:)-xorig(:)
         ss=vnorm(rr,0)
         if (ss > maxcpdst) then
           write(std_out,*) 'RCP distance from atom exceed maxcpdst !'
           cycle srcp1
         end if
!        Check that it is a RCP
         nn=0
         do kk=1,3
           nn=nn+ev(kk)/abs(ev(kk))
         end do
         if (nn /= 1) then
           write(std_out,*) ' the critical point that is found is not a RCP '
           cycle srcp1
         end if
!        Might be the same RCP than one already found on the same processor
         if (nrcp > 0) then
           do kk=1,nrcp
             pom(:)=vv(:)-rcp(kk)%rr(:)-xorig(:)
             dist=vnorm(pom,0)
             if (dist < aim_dtset%dpclim) then
               write(std_out,*) ':RCP already known'
               cycle srcp1
             end if
           end do
         end if
!        If crit==2, check that it is on the Bader surface
         if (aim_dtset%crit==2) then
           uu(:)=vv(:)-aim_epstep*rr(:)/ss
           call aim_follow(aim_dtset,uu,aim_npmaxin,srch,iatinit,iposinit,iat,ipos,nstep)
           if ((iat/=inxat).or.(inxcell/=ipos))then
             write(std_out,*) ' RCP is not on the Bader surface (outside of it)'
             cycle srcp1
           end if
         end if
         nrcp=nrcp+1
         rcp(nrcp)%rr(:)=vv(:)-xorig(:)

         buffer(1:3,ipair)=vv
         buffer(4:6,ipair)=ev
         buffer(7:9,ipair)=evec(:,1)
         buffer(10:12,ipair)=evec(:,2)
         buffer(13:15,ipair)=evec(:,3)
         buffer(16,ipair)=one

!        DEBUG
!        write(std_out,*)ch10,ch10,' drvcpr : ipair,candidate=',ipair,candidate
!        ENDDEBUG
       end if
     end do srcp1
   end do
   call xmpi_sum(buffer,xmpi_world,ierr)

   nrcp=0
   ipair=0
   do ii=1,nb_now-1
     srcp: do jj=ii+1,nb_now
       ipair=ipair+1
       candidate=buffer(16,ipair)

!      One CP has been found, must make tests to see whether it is a new RCP
       if (nint(candidate)==1) then

         vv=buffer(1:3,ipair)
         ev=buffer(4:6,ipair)
         evec(:,1)=buffer(7:9,ipair)
         evec(:,2)=buffer(10:12,ipair)
         evec(:,3)=buffer(13:15,ipair)

!        Check that it is not the same as a previous one
         if (nrcp > 0) then
           do kk=1,nrcp
             pom(:)=vv(:)-rcp(kk)%rr(:)-xorig(:)
             dist=vnorm(pom,0)
             if (dist < aim_dtset%dpclim) then
               write(std_out,*) ':RCP already known'
               cycle srcp
             end if
           end do
         end if

!        A new RCP has been found !
         nrcp=nrcp+1

!        DEBUG
!        write(std_out,*)' drvcpr : A new RCP has been found, for kk=',kk
!        ENDDEBUG


         rcp(nrcp)%iat=iat
         rcp(nrcp)%ipos=ipos
         rcp(nrcp)%rr(:)=vv(:)-xorig(:)
         rcp(nrcp)%vec(:,:)=evec(:,:)
         rcp(nrcp)%ev(:)=ev(:)
         rcp(nrcp)%vv(:)=vv(:)
         call vgh_rho(vv,chg,grho,hrho,dist,iat,ipos,0)
         rcp(nrcp)%chg=chg

       end if ! ires==0
     end do srcp ! jj=ii+2,nb_now
   end do ! ii=1,nb_now-1

   ABI_DEALLOCATE(buffer)

   if(nrcp>0)then

!    Order the RCP. CPs should appear by increasing values of x,y,z , the latter
!    varying the fastest
     ABI_ALLOCATE(sortguide,(nrcp))
     ABI_ALLOCATE(indexcp,(nrcp))
     ABI_DATATYPE_ALLOCATE(cp_tmp,(nrcp))
     do ii=3,1,-1
!      DEBUG
!      write(std_out,*)' cpdrv : sort on index ii=',ii
!      ENDDEBUG
       do jj=1,nrcp

!        DEBUG
!        write(std_out,*)rcp(jj)%vv(:)
!        ENDDEBUG

!        Try to be platform-independent. Might need a larger tolerance.
         sortguide(jj)=rcp(jj)%vv(ii)
         indexcp(jj)=jj
       end do
       call sort_dp(nrcp,sortguide,indexcp,tol3)
       do jj=1,nrcp
         cp_tmp(jj)=rcp(indexcp(jj))
       end do
       do jj=1,nrcp
         rcp(jj)=cp_tmp(jj)
       end do
     end do

!    DEBUG
!    write(std_out,*)' cpdrv : after the sort '
!    do jj=1,nrcp
!    write(std_out,*)rcp(jj)%vv(:)
!    end do
!    ENDDEBUG


!    Write the Ring Critical Point information
     do jj=1,nrcp
       write(untout,'(";Ring CP: ",3F16.8)') rcp(jj)%vv(:)
       write(untout,'("%Eigenval. of Hessian: ",3F16.8)') rcp(jj)%ev(:)
       write(untout,'(a,a,a,3f16.8,a,a,3f16.8,a,a,3f16.8,a)') &
&       ' Eigenvec. of Hessian:',char(10),&
&       '-',rcp(jj)%vec(1,:),char(10),&
&       '-',rcp(jj)%vec(2,:),char(10),&
&       '-',rcp(jj)%vec(3,:),char(10)
       write(untout,'("%Density and laplacian in CP: ",2F16.8)') &
&       rcp(jj)%chg, rcp(jj)%ev(1)+rcp(jj)%ev(2)+rcp(jj)%ev(3)
       write(untout,*) "********************************************************************"
       write(std_out,'(/," RCP: ",3F10.6,3E12.4,E12.4,/)') &
&       rcp(jj)%rr(:),rcp(jj)%ev(:),rcp(jj)%ev(1)+rcp(jj)%ev(2)+rcp(jj)%ev(3)
     end do

     ABI_DATATYPE_DEALLOCATE(cp_tmp)
     ABI_DEALLOCATE(indexcp)
     ABI_DEALLOCATE(sortguide)

   end if ! nrcp>0

   write(untout,'(" Number of RCP found: ",I4)') nrcp
   write(std_out,*) ' Number of RCP:', nrcp

!  SEARCHING CCP

   write(std_out,*)
   write(std_out,*) "CAGE CRITICAL POINTS (CCP)"
   write(std_out,*) "============================="
   write(std_out,*)

   write(untout,*)
   write(untout,*) "CAGE CRITICAL POINTS (CCP)"
   write(untout,*) "============================="
   write(untout,*)


   nn=nrcp*(nrcp-1)/2
   ABI_DATATYPE_ALLOCATE(ccp,(nn))

   nccp=0
   do ii=1,nrcp-1
     srccp: do jj=ii+1,nrcp
       vv(:)=xorig(:)+(rcp(ii)%rr(:)+rcp(jj)%rr(:))/2._dp
       call critic(aim_dtset,vv,ev,evec,aim_dmaxcs,ires,3)
       if (ires==0) then
         rr(:)=vv(:)-xorig(:)
         ss=vnorm(rr,0)
         if (ss > maxcpdst) then
           write(std_out,*) 'CCP distance from atom exceed maxcpdst !'
           cycle srccp
         end if
         nn=0
         do kk=1,3
           nn=nn+ev(kk)/abs(ev(kk))
         end do
         if (nn /= 3) then
           write(std_out,*) ' the critical point that is found is not a CCP '
           cycle srccp
         end if

         if (nccp > 0) then
           do kk=1,nccp
             pom(:)=vv(:)-ccp(kk)%rr(:)-xorig(:)
             dist=vnorm(pom,0)
             if (dist < aim_dtset%dpclim) then
               write(std_out,*) ':CCP already known'
               cycle srccp
             end if
           end do
         end if
         if (aim_dtset%crit==2) then
           uu(:)=vv(:)-aim_epstep*rr(:)/ss
           call aim_follow(aim_dtset,uu,aim_npmaxin,srch,iatinit,iposinit,iat,ipos,nstep)
           if ((iat/=inxat).or.(inxcell/=ipos)) then
             write(std_out,*) ' This CCP is not on the Bader surface (outside of it)'
             cycle srccp
           end if
         end if

         nccp=nccp+1

         ccp(nccp)%iat=iat
         ccp(nccp)%ipos=ipos
         ccp(nccp)%rr(:)=vv(:)-xorig(:)
         ccp(nccp)%vec(:,:)=evec(:,:)
         ccp(nccp)%ev(:)=ev(:)
         ccp(nccp)%vv(:)=vv(:)
         call vgh_rho(vv,chg,grho,hrho,dist,iat,ipos,0)
         ccp(nccp)%chg=chg

       end if
     end do srccp
   end do

   if(nccp>0)then

!    Order the CCP. CPs should appear by increasing values of x,y,z , the latter
!    varying the fastest
     ABI_ALLOCATE(sortguide,(nccp))
     ABI_ALLOCATE(indexcp,(nccp))
     ABI_DATATYPE_ALLOCATE(cp_tmp,(nccp))
     do ii=3,1,-1
       do jj=1,nccp
!        Try to be platform-independent. Might need a larger tolerance.
         sortguide(jj)=ccp(jj)%vv(ii)
         indexcp(jj)=jj
       end do
       call sort_dp(nccp,sortguide,indexcp,tol3)
       do jj=1,nccp
         cp_tmp(jj)=ccp(indexcp(jj))
       end do
       do jj=1,nccp
         ccp(jj)=cp_tmp(jj)
       end do
     end do

!    Write the Cage Critical Point information
     do jj=1,nccp
       write(untout,'("%Cage CP: ",3F16.8)') ccp(jj)%vv(:)
       write(untout,'("%Eigenval. of Hessian: ",3F16.8)') ccp(jj)%ev(:)
       write(untout,'(a,a,a,3f16.8,a,a,3f16.8,a,a,3f16.8,a)') &
&       ' Eigenvec. of Hessian:',char(10),&
&       '-',ccp(jj)%vec(1,:),char(10),&
&       '-',ccp(jj)%vec(2,:),char(10),&
&       '-',ccp(jj)%vec(3,:),char(10)
       write(untout,'("%Density and laplacian in CP: ",2F16.8)') &
&       ccp(jj)%chg, ccp(jj)%ev(1)+ccp(jj)%ev(2)+ccp(jj)%ev(3)
       write(untout,*) "********************************************************************"
       write(std_out,'(/," CCP: ",3F10.6,3E12.4,E12.4,/)') &
&       ccp(jj)%rr(:),ccp(jj)%ev(:),ccp(jj)%ev(1)+ccp(jj)%ev(2)+ccp(jj)%ev(3)
     end do

     ABI_DEALLOCATE(sortguide)
     ABI_DEALLOCATE(indexcp)
     ABI_DATATYPE_DEALLOCATE(cp_tmp)

   end if ! nccp>0

   write(untout,'(" Number of CCP found: ",I4)') nccp
   write(std_out,*) 'Number of CCP:', nccp
   write(std_out,*)
   write(untout,*)
   write(std_out, '(a,3i8)' ) 'BCP-RCP-CCP', nbcp,nrcp,nccp
   write(untout, '(a,3i8)' ) 'BCP-RCP-CCP', nbcp,nrcp,nccp

   write(std_out,*)
   write(std_out,*) "==============================="
   write(std_out,*) "END OF CRITICAL POINTS ANALYSIS"
   write(std_out,*)

   write(untout,*)
   write(untout,*) "==============================="
   write(untout,*) "END OF CRITICAL POINTS ANALYSIS"
   write(untout,*)


!  Output of the CPs

   write(untc,'(I4, " :BCP''s, coordinates, laplacian eigs, type of bonding at., sum of lap.eigs., density")') nbcp
   do ii=1,nbcp
     write(untc,'(3F10.6,3E12.4,I4,2E12.4)') &
&     bcp(ii)%rr(:),bcp(ii)%ev(:),bcp(ii)%iat,bcp(ii)%ev(1)+bcp(ii)%ev(2)+bcp(ii)%ev(3),bcp(ii)%chg
   end do

   write(untc,'(I4, " :RCP''s, coordinates, laplacian eigenvalues, sum of these, density")') nrcp
   do ii=1,nrcp
     vv(:)=rcp(ii)%rr(:)+xorig(:)
     call vgh_rho(vv,chg,grho,hrho,dist,iat,ipos,0)
     write(untc,'(3F10.6,3E12.4,2E12.4)') &
&     rcp(ii)%rr(:),rcp(ii)%ev(:),rcp(ii)%ev(1)+rcp(ii)%ev(2)+rcp(ii)%ev(3),rcp(ii)%chg
   end do

   write(untc,'(I4, " :CCP''s coordinates, laplacian eigenvalues, sum of these, density")') nccp
   do ii=1,nccp
     vv(:)=ccp(ii)%rr(:)+xorig(:)
     call vgh_rho(vv,chg,grho,hrho,dist,iat,ipos,0)
     write(untc,'(3F10.6,3E12.4,2E12.4)') &
&     ccp(ii)%rr(:),ccp(ii)%ev(:),ccp(ii)%ev(1)+ccp(ii)%ev(2)+ccp(ii)%ev(3),ccp(ii)%chg
   end do

 end if ! End the condition on aim_dtset%crit > 0

!Reading of the CPs from the file

 if (aim_dtset%crit==-1) then
   read(untc,*) nbcp
   ABI_DATATYPE_ALLOCATE(bcp,(nbcp))
   do ii=1,nbcp
     read(untc,*) bcp(ii)%rr(:)
   end do
   read(untc,*) nrcp
   ABI_DATATYPE_ALLOCATE(rcp,(nrcp))
   do ii=1,nrcp
     read(untc,*) rcp(ii)%rr(:)
   end do
   read(untc,*) nccp
   ABI_DATATYPE_ALLOCATE(ccp,(nccp))
   do ii=1,nccp
     read(untc,*) ccp(ii)%rr(:)
   end do
 end if

 do ii=1,nbcp
   pc(:,ii)=bcp(ii)%rr(:)
   icpc(ii)=-1
 end do
 do ii=1,nrcp
   pc(:,nbcp+ii)=rcp(ii)%rr(:)
   icpc(nbcp+ii)=1
 end do
 do ii=1,nccp
   pc(:,nbcp+nrcp+ii)=ccp(ii)%rr(:)
   icpc(nbcp+nrcp+ii)=3
 end do
 npc=nbcp+nrcp+nccp

!Checking

 if (allocated(bcp)) then
   do ii=1,nbcp
     do jj=1,npc
       iat=bcp(ii)%iat
       ipos=bcp(ii)%ipos
       if ((iat/=0).and.(ipos/=0)) then
         pom(:)=pc(:,jj)+xorig(:)-xatm(:,iat)-atp(:,ipos)
         ss=aim_dtset%coff2*vnorm(pom,0)
         if (rminl(iat) >= ss) rminl(iat)=ss
       end if
     end do
   end do
   ABI_DATATYPE_DEALLOCATE(bcp)
 end if
 do ii=1,natom
   write(std_out,*) 'atom: ', ii, rminl(ii)
 end do

 if(allocated(rcp)) then
   ABI_DATATYPE_DEALLOCATE(rcp)
 end if
 if(allocated(ccp)) then
   ABI_DATATYPE_DEALLOCATE(ccp)
 end if

!END CP ANALYSIS

 call timein(ttcp,wall)
 ttcp=ttcp-tt0

end subroutine cpdrv
!!***

!!****f* m_bader/critic
!! NAME
!! critic
!!
!! FUNCTION
!!     Search for a critical point starting from point vv
!!
!! INPUTS
!! aim_dtset= the structured entity containing all input variables
!! dmax= maximal step
!! sort= 0(default) general CP searching (Newton-Raphson)
!!                  -1,1,3 searching of specific type CP (Popelier)
!!
!! OUTPUT
!! ev= eigenvalues (ordered) of the Hessian in the final point
!! zz=  eigenvectors of the Hessian in the final point
!! ires= if ires==0 => CP found
!!       if ires==1 => CP not found within the maximum steps
!!
!! SIDE EFFECTS
!! vv(3)= starting point and final point
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine critic(aim_dtset,vv,ev,zz,dmax,ires,sort)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: sort
 integer,intent(out) :: ires
 real(dp),intent(in) :: dmax
!arrays
 real(dp),intent(inout) :: vv(3)
 real(dp),intent(out) :: ev(3),zz(3,3)
!no_abirules
 type(aim_dataset_type), intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: iat,id,ii,info,ipos,istep,jii,jj,nrot
 real(dp),parameter :: evol=1.d-3
 real(dp) :: chg,dg,dltcmax,dv,dvold,rr,ss
 logical :: oscl,outof
!arrays
 integer :: ipiv(3)
 real(dp) :: dc(3),ff(3),grho(3),hrho(3,3),lp(3),vold(3),vt(3),yy(3,3)
 real(dp),allocatable :: lamb(:),pom(:,:),pom2(:,:)

!************************************************************************

!DEBUG
!write(std_out,*)' critic : enter '
!ENDDEBUG
 oscl=.false.
 if (sort==3) then
   ABI_ALLOCATE(pom,(4,4))
   ABI_ALLOCATE(pom2,(4,4))
   ABI_ALLOCATE(lamb,(4))
 elseif (sort/=0) then
   ABI_ALLOCATE(pom,(3,3))
   ABI_ALLOCATE(pom2,(3,3))
   ABI_ALLOCATE(lamb,(3))
 end if


 deb=.false.
 istep=0
 ires=0

!DEBUG
!write(std_out,'(":POSIN ",3F16.8)') vv
!do jj=1,3
!vt(jj)=rprimd(1,jj)*vv(1)+rprimd(2,jj)*vv(2)+rprimd(3,jj)*vv(3)
!end do
!write(std_out,'(":RBPOSIN ",3F16.8)') vt
!ENDDEBUG

 call vgh_rho(vv,chg,grho,hrho,rr,iat,ipos,0)

!write(std_out,'(":GRAD ",3F16.8)') grho
!write(std_out,'(":HESSIAN ",/,3F16.8,/,3F16.8,/,3F16.8)') ((hrho(ii,jj),jj=1,3),ii=1,3)

 dg=1.0_dp
 dv=1.0_dp
 dg = vnorm(grho,0)

 if (chg < aim_rhomin) then
   ires=1
!  DEBUG
!  write(std_out,*)' critic : exit, ires=1'
!  ENDDEBUG
   return
 end if

!main cycle => limits (adhoc):
!aim_dtset%lstep - minimal step
!aim_dtset%lgrad - minimal norm of gradient
!aim_maxstep - max number of steps

 do while ((dv>aim_dtset%lstep).and.(dg>aim_dtset%lgrad).and.(istep<aim_maxstep))
   istep=istep+1
   vold(:)=vv(:)
   dvold=dv
   ev(:)=0._dp
   yy(:,:)=0._dp
   call jacobi(hrho,3,3,ev,yy,nrot)   ! eigenval of Hessian
   call ordr(ev,yy,3,-1)  ! ordering

!  modification of the Newton-Raphson step to searching
!  specific type of CP (Popelier algorithm)

   ff(:)=0._dp
   lp(:)=0._dp
   dc(:)=0._dp
   outof=.false.
   do ii=1,3
     do jj=1,3
       ff(ii)=ff(ii)+yy(jj,ii)*grho(jj)
     end do
   end do
   id=sign(1._dp,ev(1))+sign(1._dp,ev(2))+sign(1._dp,ev(3))
   if (id /= sort) then
     outof=.true.
     select case (sort)
     case (-1)
       lp(3)=0.5_dp*(ev(3)-sqrt(ev(3)*ev(3)+4._dp*ff(3)*ff(3)))
       pom(:,:)=0._dp
       pom2(:,:)=0._dp
       lamb(:)=0._dp
       do ii=1,2
         pom(ii,ii)=ev(ii)
         pom(ii,3)=ff(ii)
         pom(3,ii)=ff(ii)
       end do
       call jacobi(pom,3,3,lamb,pom2,nrot)
       call ordr(lamb,pom2,3,1)
       do ii=1,3
         lp(1)=lamb(ii)
         if (abs(pom2(3,ii))>1.0d-24) exit
       end do
       lp(2)=lp(1)

!        write(std_out,*) (ev(ii),ii=1,3)
!        write(std_out,*) (lamb(ii),ii=1,3)
!        write(std_out,*) ':ID  ',id,lp(1),lp(3)

     case (1)
       lp(1)=0.5_dp*(ev(1)+sqrt(ev(1)*ev(1)+4._dp*ff(1)*ff(1)))
       pom(:,:)=0._dp
       pom2(:,:)=0._dp
       lamb(:)=0._dp
       do ii=2,3
         pom(ii-1,ii-1)=ev(ii)
         pom(ii-1,3)=ff(ii)
         pom(3,ii-1)=ff(ii)
       end do
       call jacobi(pom,3,3,lamb,pom2,nrot)
       call ordr(lamb,pom2,3,1)
       do ii=3,1,-1
         lp(2)=lamb(ii)
         if (abs(pom2(3,ii))>1.0d-24) exit
       end do
       lp(3)=lp(2)

     case (3)
       pom(:,:)=0._dp
       pom2(:,:)=0._dp
       lamb(:)=0._dp
       do ii=1,3
         pom(ii,ii)=ev(ii)
         pom(ii,4)=ff(ii)
         pom(4,ii)=ff(ii)
       end do
       call jacobi(pom,4,4,lamb,pom2,nrot)
       call ordr(lamb,pom2,4,1)
       do ii=4,1,-1
         lp(1)=lamb(ii)
         if (abs(pom2(4,ii))>1.0d-24) exit
       end do
       lp(2)=lp(1); lp(3)=lp(1)
     case default
       lp(:)=0._dp
     end select
   end if

   do ii=1,3
     if (abs(ev(ii)-lp(ii))<1.0d-24) then
       outof=.false.
       exit
     end if
   end do
   do ii=1,3                      ! SEARCHING STEP
     do jj=1,3
       if (outof) then
         dc(ii)=dc(ii)+ff(jj)*yy(ii,jj)/(ev(jj)-lp(jj))
       elseif (abs(ev(jj))>1.0d-24) then
         dc(ii)=dc(ii)+ff(jj)*yy(ii,jj)/ev(jj)
       else
         MSG_ERROR("zero eigval of Hessian")
       end if
     end do
   end do

   dltcmax = vnorm(dc,0)
   if (dltcmax>dmax) then                 ! STEP RESTRICTION
     do ii=1,3
       dc(ii)=dc(ii)*dmax/dltcmax
     end do
   end if                                  ! primitive handling of oscillations
   ss=vnorm(dc,0)                          ! usually not needed
   ss=abs(ss-dv)/ss
   if ((ss < evol).and.(oscl)) then
     dc(:)=dc(:)/2._dp
   end if


   do ii=1,3
     vv(ii) = vv(ii) - dc(ii)
   end do

!  DEBUG
!  write(std_out,'(":POSIN ",3F16.8)') vv
!  ENDDEBUG

   call vgh_rho(vv,chg,grho,hrho,rr,iat,ipos,0)
   dg = vnorm(grho,0)

   if (deb) then                 !  DEBUGG OUTPUT
     write(std_out,'("AFTER STEP ===================================")')
     write(std_out,'(":HESSIAN^(-1) ",/,3F16.8,/,3F16.8,/,3F16.8)') ((yy(ii,jii),jii=1,3),ii=1,3)
     write(std_out,'(":DC ",3F16.8)') dc
     write(std_out,*) 'STEP ',istep
     write(std_out,'(":POS ",3F16.8)') vv
     write(std_out,'(":GRAD ",3F16.8)') grho
     write(std_out,*) ':DGRAD,CHG ',dg,chg
     write(std_out,'(":HESSIAN ",/,3F16.8,/,3F16.8,/,3F16.8)') ((hrho(ii,jii),jii=1,3),ii=1,3)
   end if
   vt(:)=vv(:)-vold(:)
   dv=vnorm(vt,0)
   ss=abs(dvold-dv)/dv
   if (ss < evol) oscl=.true.
 end do

!end of main cycle

!the final output

 write(std_out,*) 'iste:',istep, dv, dg
 if (istep>=aim_maxstep)then
   write(std_out,*) ' istep=MAXSTEP ! Examine lstep2 and lgrad2 .'
   if ( (dv>aim_dtset%lstep2) .and. (dg>aim_dtset%lgrad2 )) then
     write(std_out,'(":POSOUT ",3F16.8)') vv
     ires=1
   end if
 end if

 vt(:)=vv(:)

!write(std_out,'(":POSOUT ",3F16.8)') vv
!write(std_out,'(":RBPOSOUT ",3F16.8)') vt

 call vgh_rho(vv,chg,grho,hrho,rr,iat,ipos,0)

!write(std_out,'(":GRAD ",3F16.8)') grho
!write(std_out,'(":HESSIAN ",/,3F16.8,/,3F16.8,/,3F16.8)')&
!& ((hrho(ii,jii),jii=1,3),ii=1,3)


!FINAL INVERSION OF HESSIAN

 call ludcmp(hrho,3,3,ipiv,id,info)
 if (info /= 0) then
   write(std_out,*) 'Error inverting hrho:'
   do ii=1,3
     write(std_out,*) (hrho(ii,jii),jii=1,3)
   end do
   ires=1
!  DEBUG
!  write(std_out,*)' critic : exit, ires=1'
!  ENDDEBUG
   return
!  stop 'ERROR INVERTING HESSIAN'
 end if
 do ii=1,3
   yy(ii,1:3)=0.
   yy(ii,ii)=1.
 end do
 do jii=1,3
   call lubksb(hrho,3,3,ipiv,yy(1,jii))
 end do


!write(std_out,'(":HESSIAN^(-1) ",/,3F16.8,/,3F16.8,/,3F16.8)') ((y(ii,jii),jii=1,3),ii=1,3)

 call vgh_rho(vv,chg,grho,hrho,rr,iat,ipos,0)

!write(std_out,'("LAPLAC:",F16.8)') hrho(1,1)+hrho(2,2)+hrho(3,3)

 call jacobi(hrho,3,3,ev,yy,nrot)
 call ordr(ev,yy,3,1)
 zz(:,:)=yy(:,:)

!do ii=1,3
!do jii=1,3
!zz(ii,jii)=yy(jii,ii)
!end do
!end do

!write(std_out,'(":AUTOVAL ",3F16.8)') (ev(ii),ii=1,3)
!write(std_out,'(":AUTOVEC ",/,3F16.8,/,3F16.8,/,3F16.8)') ((zz(ii,jii),ii=1,3),jii=1,3)

 if (sort/=0)  then
   ABI_DEALLOCATE(pom)
   ABI_DEALLOCATE(pom2)
   ABI_DEALLOCATE(lamb)
 end if

!DEBUG
!write(std_out,*)' critic : exit, ires= ',ires
!ENDDEBUG
end subroutine critic
!!***

!!****f* m_bader/ordr
!! NAME
!! ordr
!!
!! FUNCTION
!!
!! INPUTS
!!  (to be filled)
!!
!! OUTPUT
!!  (to be filled)
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE
!!

subroutine ordr(aa,dd,nn,cff)

!Arguments ----------------------------
!scalars
 integer,intent(in) :: cff,nn
!arrays
 real(dp),intent(inout) :: aa(nn),dd(nn,nn)

!Local variables ----------------------
!scalars
 integer :: ii,jj,kk
 real(dp) :: uu

! *********************************************************************

 do ii=1,nn-1
   kk=ii
   uu=aa(ii)
   do jj=ii+1,nn
     if (cff==1) then
       if (aa(jj) >= uu+tol12) then
         kk=jj
         uu=aa(jj)
       end if
     else
       if (aa(jj) <= uu-tol12) then
         kk=jj
         uu=aa(jj)
       end if
     end if
   end do
   if (kk /= ii) then
     aa(kk)=aa(ii)
     aa(ii)=uu
     do jj=1,nn
       uu=dd(jj,ii)
       dd(jj,ii)=dd(jj,kk)
       dd(jj,kk)=uu
     end do
   end if
 end do
end subroutine ordr
!!***

!!****f* m_bader/critics
!! NAME
!! critics
!!
!! FUNCTION
!! Search for critical points starting between
!!    atom inxat and its neighbors.
!!
!! INPUTS
!!  aim_dtset= the structured entity containing all input variables
!!  dstmax=maximum distance to search for neighbors
!!  stwo, sthree, sfour: logical switches (TRUE/FALSE) indicating
!!                          to search CP starting in the middle point
!!                          of two, three or four atoms. One of these
!!                          atoms is inxat.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  This routines acts primarily on the data contained in the aim_prom module
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine  critics(aim_dtset,inxat,stwo,sthree,sfour,dstmax)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: inxat
 real(dp),intent(in) :: dstmax
 logical,intent(in) :: sfour,sthree,stwo
!no_abirules
 type(aim_dataset_type), intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3,iat,ii,ipos,ires,jii,jj,kjj,kk,ll,n1,n2,n3,nb,nc
 integer :: nshell
 real(dp) :: chg,dif1,dif2,diff,dist,olddist,rr
! real(dp) :: ss,uu
 logical :: found,inter
!arrays
 integer :: ibat(nnpos*natom),inat(nnpos*natom),ipibat(nnpos*natom)
 integer :: nnat(nnpos*natom),nr(nnpos*natom)
 real(dp) :: dif(3),dists(nnpos*natom),ev(3),grho(3),hrho(3,3)
 real(dp) :: pom(3),v1(3),v2(3),v3(3),v4(3),vi(3),vt(3),zz(3,3)

!************************************************************************
 vi(:)=xatm(:,inxat)

 nc=0
 do jii=1,nnpos
   do kjj=1,natom
     dist=0._dp
     dif(:)=xatm(:,inxat)-xatm(:,kjj)-atp(:,jii)

!    do ii=1,3
!    dif(ii)=xatm(ii,inxat)-xatm(ii,kjj)-atp(ii,jii)
!    end do
     dist=vnorm(dif,0)
     if (.not.((dist>dstmax).or.(dist<0.001))) then
       nc=nc+1
       dists(nc)=dist
       nnat(nc)=kjj
       inat(nc)=jii
     end if
   end do
 end do
 do n1=1,nc
   nr(n1)=n1
 end do
 call sort_dp(nc,dists,nr,tol14)
 nb=0
 olddist=0._dp
 nshell=0
!write(std_out,*) ':ORIAT ', (xatm(ii,inxat),ii=1,3)
 do n1=1,nc
   n2=nr(n1)
   n3=nnat(n2)
   if (dists(n1)<(2*dists(1))) then
     if ((dists(n1)-olddist)>aim_dlimit) then
       nshell=nshell+1
       olddist=dists(n1)
       if (nshell==5) exit
     end if
     nb=nb+1
     ibat(nb)=n3
     ipibat(nb)=inat(n2)
     write(std_out,*) ':NEIG ',inxat,n3,inat(n2),dists(n1)
!    write(std_out,*) ':POSAT',(xatm(ii,ibat(nb))+atp(ii,ipibat(nb)),ii=1,3)
   else
     exit
   end if
 end do

 npc=0
 npcm3=0

!
!.....SEARCH BETWEEN EACH PAIR OF ATOMS
!

 if (stwo) then
   do jii=1,nb
     do ii=1,3
       v1(ii)=xatm(ii,inxat)
       v2(ii)=xatm(ii,ibat(jii))+atp(ii,ipibat(jii))
       vt(ii)=(v1(ii)+v2(ii))/2._dp
     end do
     inter=.true.
     diff=0._dp
     pom(:)=vt(:)
     pom(:)=pom(:)-vi(:)
     diff=vnorm(pom,0)
     if (diff > maxcpdst) inter=.false.
     if (inter) then
       call critic(aim_dtset,vt,ev,zz,aim_dmaxcs,ires,0)
       if (ires==0) then
         found=.false.
         if (npc > 0) then
           do jj=1,npc
             pom(:)=vt(:)-pc(:,jj)
             dist=vnorm(pom,0)
             if (dist < aim_dtset%dpclim) found=.true.
           end do
         end if
         if (.not.found) then
           pom(:)=vt(:)
           call bschg1(pom,-1)
           pcrb(:,npc+1)=pom(:)
           pom(:)=pom(:)-vi(:)
           diff=vnorm(pom,0)
           if (abs(diff) > maxcpdst) found=.true.
         end if
         if (.not.found) then
           npc=npc+1
           do jj=1,3
             pc(jj,npc)=vt(jj)
             evpc(jj,npc)=ev(jj)
             do kk=1,3
               zpc(kk,jj,npc)=zz(kk,jj)
             end do
           end do
           i1=ev(1)/abs(ev(1))
           i2=ev(2)/abs(ev(2))
           i3=ev(3)/abs(ev(3))
           icpc(npc)=i1+i2+i3
           if (icpc(npc)==-3) then
             npcm3=npcm3+1
           end if
           write(std_out,*) 'New critical point found'
           write(std_out,'("POS: ",3F16.8)') (pc(ii,npc),ii=1,3)
           write(std_out,'("POS in base: ",3F16.8)') (pcrb(ii,npc),ii=1,3)
           write(std_out,'("AUTOVAL: ",3F16.8)') ev
           write(std_out,'("AUTOVEC: ",/,3F16.8,/,3F16.8,/,3F16.8)') &
&           ((zpc(ii,jj,npc),ii=1,3),jj=1,3)
           call vgh_rho(vt,chg,grho,hrho,rr,iat,ipos,0)
           write(22,'(":PC2",3F10.6,3E12.4,I4,2E12.4)') &
&           (pc(jj,npc),jj=1,3),(ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg
           write(std_out,'(":PC2",3F10.6,3E12.4,I4,2E12.4)')  &
&           (pc(jj,npc),jj=1,3),(ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg
           pom(:)=vt(:)-v1(:)
           dif1=vnorm(pom,0)
           pom(:)=vt(:)-v2(:)
           dif2=vnorm(pom,0)
           write(std_out,'(":DISPC ",2F12.8)') dif1,dif2
         end if
       end if
     end if
   end do
 end if
!
!.....SEARCH BETWEEN EACH THREE ATOMS
!
 if(sthree) then
   do jii=1,nb
     do kjj=jii+1,nb
       do ii=1,3
         v1(ii)=xatm(ii,inxat)
         v2(ii)=xatm(ii,ibat(jii))+atp(ii,ipibat(jii))
         v3(ii)=xatm(ii,ibat(kjj))+atp(ii,ipibat(kjj))
         vt(ii)=(v1(ii)+v2(ii)+v3(ii))/3._dp
       end do
       inter=.true.
       pom(:)=vt(:)
       pom(:)=pom(:)-vi(:)
       dist=vnorm(pom,0)
       if (abs(diff)>maxcpdst) then
         inter=.false.
         exit
       end if
       if (inter) then
         do jj=1,npc
           pom(:)=pc(:,jj)-vt(:)
           diff=vnorm(pom,0)
           if (diff<aim_dpc0) then
             inter=.false.
             exit
           end if
         end do
       end if
       if (inter) then
         call critic(aim_dtset,vt,ev,zz,aim_dmaxcs,ires,0)
         if (ires==0) then
           found=.false.
           if (npc>0) then
             do jj=1,npc
               pom(:)=vt(:)-pc(:,jj)
               dist=vnorm(pom,0)
               if (dist<aim_dtset%dpclim) then
                 found=.true.
                 exit
               end if
             end do
           end if
           if (.not.found) then
             pom(:)=vt(:)
             call bschg1(pom,-1)
             pcrb(:,npc+1)=pom(:)
             pom(:)=pom(:)-vi(:)
             diff=vnorm(pom,0)
             if (abs(diff)>maxcpdst) found=.true.
           end if
           if (.not.found) then
             npc=npc+1
             do jj=1,3
               pc(jj,npc)=vt(jj)
               evpc(jj,npc)=ev(jj)
               do kk=1,3
                 zpc(kk,jj,npc)=zz(kk,jj)
               end do
             end do
             i1=ev(1)/abs(ev(1))
             i2=ev(2)/abs(ev(2))
             i3=ev(3)/abs(ev(3))
             icpc(npc)=i1+i2+i3
             if (icpc(npc)==-3) then
               npcm3=npcm3+1
             end if
             write(std_out,*) 'New critical point found'
             write(std_out,'("POS: ",3F16.8)') (pc(ii,npc),ii=1,3)
             write(std_out,'("POS in base: ",3F16.8)') (pcrb(ii,npc),ii=1,3)
             write(std_out,'("AUTOVAL: ",3F16.8)') ev
             write(std_out,'("AUTOVEC: ",/,3F16.8,/,3F16.8,/,3F16.8)') &
&             ((zpc(ii,jj,npc),ii=1,3),jj=1,3)
             call vgh_rho(vt,chg,grho,hrho,rr,iat,ipos,0)
             write(22,'(":PC3",3F10.6,3E12.4,I4,2E12.4)') &
&             (pc(jj,npc),jj=1,3),(ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg
             write(std_out,'(":PC3",3F10.6,3E12.4,I4,2E12.4)') &
&             (pc(jj,npc),jj=1,3),(ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg
           end if
         end if
       end if
     end do
   end do
 end if

!
!.....SEARCH BETWEEN EACH FOUR ATOMS
!
 if (sfour) then
   do jii=1,nb
     do kjj=jii+1,nb
       do ll=jii+1,nb
         do ii=1,3
           v1(ii)=xatm(ii,inxat)
           v2(ii)=xatm(ii,ibat(jii))+atp(ii,ipibat(jii))
           v3(ii)=xatm(ii,ibat(kjj))+atp(ii,ipibat(kjj))
           v4(ii)=xatm(ii,ibat(ll))+atp(ii,ipibat(ll))
           vt(ii)=(v1(ii)+v2(ii)+v3(ii)+v4(ii))/4._dp
         end do
         inter=.true.
         pom(:)=vt(:)
         pom(:)=pom(:)-vi(:)
         diff=vnorm(pom,0)
         if (abs(diff)>maxcpdst) then
           inter=.false.
           exit
         end if
         if (inter) then
           do jj=1,npc
             pom(:)=pc(:,jj)-vt(:)
             diff=vnorm(pom,0)
             if (diff < aim_dpc0) then
               inter=.false.
               exit
             end if
           end do
         end if
         if (inter) then
           call critic(aim_dtset,vt,ev,zz,aim_dmaxcs,ires,0)
           if (ires==0) then
             found=.false.
             if (npc>0) then
               do jj=1,npc
                 pom(:)=vt(:)-pc(:,jj)
                 dist=vnorm(pom,0)
                 if (dist < aim_dtset%dpclim) found=.true.
               end do
             end if
             if (.not.found) then
               pom(:)=vt(:)
               pcrb(:,npc+1)=pom(:)
               pom(:)=pom(:)-vi(:)
               diff=vnorm(pom,0)
               if (abs(diff)>maxcpdst) found=.true.
             end if
             if (.not.found) then
               npc=npc+1
               do jj=1,3
                 pc(jj,npc)=vt(jj)
                 evpc(jj,npc)=ev(jj)
                 do kk=1,3
                   zpc(kk,jj,npc)=zz(kk,jj)
                 end do
               end do
               i1=ev(1)/abs(ev(1))
               i2=ev(2)/abs(ev(2))
               i3=ev(3)/abs(ev(3))
               icpc(npc)=i1+i2+i3
               if (icpc(npc)==-3) then
                 npcm3=npcm3+1
               end if
               write(std_out,*) 'New critical point found'
               write(std_out,'("POS: ",3F16.8)') (pc(ii,npc),ii=1,3)
               write(std_out,'("POS in base: ",3F16.8)') (pcrb(ii,npc),ii=1,3)
               write(std_out,'("AUTOVAL: ",3F16.8)') ev
               write(std_out,'("AUTOVEC: ",/,3F16.8,/,3F16.8,/,3F16.8)') &
&               ((zpc(ii,jj,npc),ii=1,3),jj=1,3)
               call vgh_rho(vt,chg,grho,hrho,rr,iat,ipos,0)
               write(22,'(":PC4",3F10.6,3E12.4,I4,2E12.4)') &
&               (pc(jj,npc),jj=1,3), (ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg
               write(std_out,'(":PC4",3F10.6,3E12.4,I4,2E12.4)') &
&               (pc(jj,npc),jj=1,3), (ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg
             end if
           end if
         end if
       end do
     end do
   end do
 end if

 write(std_out,*) npc
end subroutine critics
!!***

!!****f* m_bader/defad
!! NAME
!! defad
!!
!! FUNCTION
!! Initialisation of aim input variables to their default values.
!!
!! INPUTS
!!  (no input : initialisation by default values)
!!
!! OUTPUT
!! aim_dtset = the structured entity containing all input variables
!!
!! PARENTS
!!      aim
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine defad(aim_dtset)

!Arguments ------------------------------------
!scalars
 type(aim_dataset_type),intent(out) :: aim_dtset

!Local variables ------------------------------

! *********************************************************************

 aim_dtset%isurf=0
 aim_dtset%crit=0
 aim_dtset%irsur=0
 aim_dtset%foll=0
 aim_dtset%irho=0
 aim_dtset%ivol=0
 aim_dtset%denout=0
 aim_dtset%lapout=0
 aim_dtset%gpsurf=0
 aim_dtset%plden=0
 aim_dtset%dltyp=0

 aim_dtset%batom=1
 aim_dtset%nsa=3
 aim_dtset%nsb=3
 aim_dtset%nsc=3
 aim_dtset%npt=100
 aim_dtset%nth=32
 aim_dtset%nph=48

 aim_dtset%themax=pi
 aim_dtset%themin=zero
 aim_dtset%phimin=zero
 aim_dtset%phimax=two_pi
 aim_dtset%phi0=zero
 aim_dtset%th0=zero
 aim_dtset%folstp=5.d-2
 aim_dtset%dr0=5.d-2
 aim_dtset%atrad=one
 aim_dtset%rmin=one

 aim_dtset%foldep(:)=zero
 aim_dtset%vpts(:,:)=zero
 aim_dtset%ngrid(:)=30
 aim_dtset%scal(:)=one
 aim_dtset%maxatd=1.d1
 aim_dtset%maxcpd=5.d1

 aim_dtset%dpclim=1.d-2
 aim_dtset%lstep=1.d-10
 aim_dtset%lstep2=1.d-5
 aim_dtset%lgrad=1.d-12
 aim_dtset%lgrad2=1.d-5
 aim_dtset%coff1=0.98_dp
 aim_dtset%coff2=0.95_dp

end subroutine defad
!!***

!!****f* m_bader/drvaim
!! NAME
!! drvaim
!!
!! FUNCTION
!! Main driver for the Bader analysis
!! it looks the values of the input variables
!! and calls corresponding procedures
!!
!! INPUTS
!! aim_dtset = the structured entity containing all input variables
!! tcpui=initial CPU time
!! twalli=initial wall clock time
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  this routine acts primarily on the data contained in the aimprom module
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      aim
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine drvaim(aim_dtset,tcpui,twalli)

!Arguments ------------------------------------
!scalars
 real(dp) :: tcpui,twalli
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: iat,iatinit,inxat,ipos,iposinit
 integer :: me,npmax,nproc,nstep
 real(dp) :: dstlim,rr,ss,t1,t2,tf,wall
 real(dp) :: tcpu,twall,znucl_batom
 logical :: debold,sfour,srch,sthree,stwo
!arrays
 real(dp) :: tsec(2)
 real(dp) :: grho(3),xstart(3)

! *********************************************************************

 me=xmpi_comm_rank(xmpi_world)
 nproc=xmpi_comm_size(xmpi_world)

!These input variables might be modified during what follows,
!so, they are copied outside of aim_dtset.
 inxat=aim_dtset%batom
 r0=aim_dtset%atrad
 h0=aim_dtset%folstp
 maxatdst=aim_dtset%maxatd
 maxcpdst=aim_dtset%maxcpd

 dstlim=maxcpdst

!Flags from the old version
!to be remove later
 deb=.false.
 stwo=.true.
 sthree=.true.
 sfour=.false.
 srch=.false.

 npmax=aim_npmaxin

!Main initialisation procedure -
!- it reads ABINIT density file and files
!with core densities and initialises the fields for
!spline interpolation

 call initaim(aim_dtset,znucl_batom)


!CP SEARCHING

 if (aim_dtset%crit /= 0) then

   call timein(tcpu,twall)
   tsec(1)=tcpu-tcpui ; tsec(2)=twall-twalli
   write(std_out, '(5a,f13.1,a,f13.1)' ) &
&   '-',ch10,'- Before searching the CP ',ch10,&
&   '- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)

   if (aim_dtset%crit==3) then
!    old version of the driver for searching CPs (original code)
     call critics(aim_dtset,inxat,stwo,sthree,sfour,dstlim)
   else
!    driver for searching CPs with Popellier algorithm
     call cpdrv(aim_dtset)
   end if

   call timein(tcpu,twall)
   tsec(1)=tcpu-tcpui ; tsec(2)=twall-twalli
   write(std_out, '(5a,f13.1,a,f13.1)' ) &
&   '-',ch10,'- After searching the CP ',ch10,&
&   '- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)

 end if

!
!BADER SURFACE CALCULATION
!

 if (aim_dtset%isurf==1) then
!  driver for determination of the Bader surface

   call timein(tcpu,twall)
   tsec(1)=tcpu-tcpui ; tsec(2)=twall-twalli
   write(std_out, '(5a,f13.1,a,f13.1)' ) &
&   '-',ch10,'- Before determinating the Bader surface ',ch10,&
&   '- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)

   call surf(aim_dtset)

   call timein(tcpu,twall)
   tsec(1)=tcpu-tcpui ; tsec(2)=twall-twalli
   write(std_out, '(5a,f13.1,a,f13.1)' ) &
&   '-',ch10,'- After determinating the Bader surface ',ch10,&
&   '- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)

 end if

!
!CHARGE INTEGRATIOM
!

 if (aim_dtset%irho==1) then
   call integrho(aim_dtset,znucl_batom)
 end if

!
!VOLUME INTEGRATION OF THE BADER ATOM
!

 if (aim_dtset%ivol==1) then
   call integvol()
 end if

!
!ONE RADIUS OF THE BADER SURFACE
!

 if (aim_dtset%irsur==1) then
   if (aim_dtset%isurf/=0) srch=.true.
   iat=aim_dtset%batom
   ss=r0
   call timein(t1,wall)
   call rsurf(aim_dtset,rr,grho,aim_dtset%th0,aim_dtset%phi0,ss,iat,npmax,srch)
   call timein(t2,wall)
   t2=t2-t1
   write(unts,'(2F12.8,F15.10)') aim_dtset%th0,aim_dtset%phi0,rr
   write(std_out,'(":RSUR ",2F12.8,2F15.10)') aim_dtset%th0,aim_dtset%phi0,rr,t2
 end if

!
!FOLLOW THE GRADIENT PATH FROM ONE POINT
!

 if (aim_dtset%foll==1) then
   iatinit=aim_dtset%batom
   iposinit=batcell
   if (aim_dtset%isurf/=0) srch=.true.
   debold=deb
   xstart(:)=aim_dtset%foldep(:)
   call timein(t1,wall)
   call aim_follow(aim_dtset,xstart,npmax,srch,iatinit,iposinit,iat,ipos,nstep)
   call timein(t2,wall)
   tf=t2-t1
   write(std_out,'(":TIME in aim_follow:", F12.4)') tf
 end if

 if (aim_dtset%plden == 1) then
!  profile of the density integrated in plane xy
!  belong the z-axes - not finished - cut3d better !
   call plint()
 end if

 if ((aim_dtset%denout > 0).or.(aim_dtset%lapout > 0)) then
!  additional outputs of density and laplacian fields
!  in the plane or line
   call addout(aim_dtset)
 end if

 if (aim_dtset%gpsurf == 1) then
!  script for gnuplot - simple demonstration of the
!  computed surface
   call graph(unts,untg)
 end if

!Deallocation of global variables allocated in initaim
!and declared in defs_aimfields.
 ABI_DEALLOCATE(dig1)
 ABI_DEALLOCATE(dig2)
 ABI_DEALLOCATE(dig3)
 ABI_DEALLOCATE(llg1)
 ABI_DEALLOCATE(llg2)
 ABI_DEALLOCATE(llg3)
 ABI_DEALLOCATE(cdig1)
 ABI_DEALLOCATE(cdig2)
 ABI_DEALLOCATE(cdig3)
 ABI_DEALLOCATE(ddx)
 ABI_DEALLOCATE(ddy)
 ABI_DEALLOCATE(ddz)
 ABI_DEALLOCATE(rrad)
 ABI_DEALLOCATE(crho)
 ABI_DEALLOCATE(sp2)
 ABI_DEALLOCATE(sp3)
 ABI_DEALLOCATE(sp4)
 ABI_DEALLOCATE(corlim)
 ABI_DEALLOCATE(dvl)
 ABI_DEALLOCATE(ndat)
 ABI_DEALLOCATE(rminl)
!Deallocation of global variables allocated in initaim
!and declared in defs_aimprom.
 ABI_DEALLOCATE(typat)
 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(xatm)

end subroutine drvaim
!!***

!!****f* m_bader/graph
!! NAME
!! graph
!!
!! FUNCTION
!! Writing  of the gnuplot script to show the computed part
!! of Bader surface with lines
!!
!! INPUTS
!!  untg = unit number of the file on which the info is written
!!  unts = unit number of the file from which the Bader surface is read
!!
!! OUTPUT
!!  (written in the untg file)
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine graph(unts,untg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: untg,unts

!Local variables ------------------------------
!scalars
 integer :: ii,indx,jj,nphi,nth
 real(dp),parameter :: snull=1.d-6
 real(dp) :: phimax,phimin,ss,thmax,thmin
!arrays
 real(dp) :: xorig(3)
 real(dp),allocatable :: phi(:),rr(:,:),th(:)

! *********************************************************************

 rewind(unts)
 read(unts,*) indx, xorig(1:3)
 read(unts,*) nth, thmin, thmax
 read(unts,*) nphi, phimin, phimax
 ABI_ALLOCATE(th,(nth))
 ABI_ALLOCATE(phi,(nphi))
 ABI_ALLOCATE(rr,(nth,nphi))
 do ii=1,nth
   do jj=1,nphi
     read(unts,*) th(ii),phi(jj),rr(ii,jj),ss
   end do
 end do

!end of reading

 write(untg,*) 'reset'
 write(untg,*) 'set st d l'
 write(untg,*) 'set ticslevel 0'
 write(untg,*) 'set title ''Bader surface'' '
 write(untg,*) 'splot ''-'' using ($3*sin($1)*cos($2)):($3*sin($1)*sin($2)):($3*cos($1)) notitle'
 do ii=1,nth
   do jj=1,nphi
     write(untg,'(2F12.8,E16.8)') th(ii),phi(jj),rr(ii,jj)
   end do
   if ((ii==nth).and.(jj==nphi)) then
     cycle
   else
     write(untg,*)
   end if
 end do

end subroutine graph
!!***

!!****f* m_bader/initaim
!! NAME
!! initaim
!!
!! FUNCTION
!! Initialization for the 3D interpolation for the AIM code:
!!  - this procedure reads the charge density of the electrons of valence on
!!    the equidistant 3D grid (*_DEN output file of ABINIT) and the core charge
!!    density of electrons from *.fc files (fhi package)
!!  - the Cholesky decomposition  of the general matrix for
!!    the computation of the 1D spline coeficients in each direction is done.
!!    Warning - the procedure is modified to use periodic boundary conditions
!!    already during the decomposition
!!  - the second derivations of valence density in three directions are computed
!!    and stored in the real space grid of the density for interpolation.
!!  - the core density is stored separately in the radial grid together with th
!!    second radial derivation
!!
!! INPUTS
!! aim_dtset= the structured entity containing all input variables
!!
!! OUTPUT
!! znucl_batom= the nuclear charge of the Bader atom
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  thie routine works on the data contained in the aim_fields and aim_prom modules
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine initaim(aim_dtset,znucl_batom)

!Arguments ------------------------------------
!scalars
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: fform0,id,ierr,ii,info,jj,kk,kod,mm,ndtmax,nn,nsa,nsb,nsc,nsym,me,nproc,npsp
 integer :: unth,comm
#ifdef HAVE_NETCDF
 integer :: den_id
#endif
 real(dp) :: ss,ucvol,znucl_batom
 real(dp) :: zz
 type(hdr_type) :: hdr
!arrays
 integer :: ipiv(3)
 integer,allocatable :: symrel(:,:,:)
 real(dp) :: aa(3),bb(3),gmet(3,3),gprimd(3,3),rmet(3,3),yy(3,3)
 real(dp),allocatable :: tnons(:,:),znucl(:),zionpsp(:)
 real(dp),pointer :: ptc(:),ptd(:),ptf(:),ptp(:),ptsd(:)

! *********************************************************************

!DEBUG
!write(std_out,*) ' initaim : enter '
!ENDDEBUG

 comm = xmpi_world
 me=xmpi_comm_rank(comm)
 nproc=xmpi_comm_size(comm)

 slc=0    ! code for follow

!The use of the "hdr" routines is much better for the future
!maintenance of the code. Indeed, the content of the header
!will continue to change from time to time, and the associated
!changes should be done in one place only.

!Read ABINIT header ----------------------------------------------------------
 if(me==master)then
   if (aim_iomode == IO_MODE_ETSF) then
     call hdr_ncread(hdr, untad, fform0)
   else
     call hdr_fort_read(hdr, untad, fform0)
   end if
 end if
 ABI_CHECK(fform0 /= 0, "hdr_read returned fform == 0")
 call hdr%bcast(master, me, comm)

!Echo part of the header
 call hdr%echo(fform0, 4, unit=std_out)
 call hdr%echo(fform0, 4, unit=untout)

 natom=hdr%natom
 ngfft(1:3)=hdr%ngfft(:)
 nsym=hdr%nsym
 npsp=hdr%npsp
 ntypat=hdr%ntypat
 rprimd(:,:)=hdr%rprimd(:,:)

 ABI_ALLOCATE(zionpsp,(npsp))
 ABI_ALLOCATE(znucl,(ntypat))
 ABI_ALLOCATE(typat,(natom))
 ABI_ALLOCATE(xred,(3,natom))
 ABI_ALLOCATE(symrel,(3,3,nsym))
 ABI_ALLOCATE(tnons,(3,nsym))
 ABI_ALLOCATE(xatm,(3,natom))

 symrel(:,:,:)=hdr%symrel(:,:,:)
 typat(:)=hdr%typat(:)
 tnons(:,:)=hdr%tnons(:,:)
 znucl(:)=hdr%znucltypat(:)
 zionpsp(:)=hdr%zionpsp(:)
 xred(:,:)=hdr%xred(:,:)

 call hdr%free()

!-------------------------------------------------------------------------------

 ABI_ALLOCATE(dvl,(ngfft(1),ngfft(2),ngfft(3)))

 if(me==master)then
   if (aim_iomode == IO_MODE_ETSF) then
#ifdef HAVE_NETCDF
     ! netcdf array has shape [cplex, n1, n2, n3, nspden]), here we read only the total density.
     NCF_CHECK(nf90_inq_varid(untad, "density", den_id))
     NCF_CHECK(nf90_get_var(untad, den_id, dvl, start=[1,1,1,1], count=[1, ngfft(1), ngfft(2), ngfft(3), 1]))
#endif
   else
     read(untad,iostat=nn) dvl(1:ngfft(1),1:ngfft(2),1:ngfft(3))
     ABI_CHECK(nn==0,"error of reading !")
   end if
 end if
 call xmpi_bcast(dvl, master, comm, ierr)

 write(std_out,*)ch10,' initaim : the valence density has been read' ,ch10

!INITIALISATION OF SOME IMPORTANT FIELDS

!Only interpolation is computed (inside vgh_rho) in reduced
!coordinates. In all other routines the cart. coordinates (CC) are used.

!transformation of the atom positions to CC
 do ii=1,natom
   xatm(:,ii)=xred(:,ii)
   call bschg1(xatm(:,ii),1)
 end do

!Generation of the neighbouring cells + transf to CC
 nn=0
 nsa=aim_dtset%nsa ; nsb=aim_dtset%nsb ; nsc=aim_dtset%nsc
 do ii=-nsa,nsa
   do jj=-nsb,nsb
     do kk=-nsc,nsc
       nn=nn+1
       atp(1,nn)=ii*1._dp
       atp(2,nn)=jj*1._dp
       atp(3,nn)=kk*1._dp
       call bschg1(atp(:,nn),1)
     end do
   end do
 end do
 nnpos=nn

!DEBUG
!write(std_out,*)' initaim : nnpos=',nnpos
!ENDDEBUG

 batcell=nsa*(2*nsb+1)*(2*nsc+1)+(2*nsc+1)*nsb+nsc+1
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 maxatdst=min(maxatdst, nsa*sqrt(rmet(1,1)), nsb*sqrt(rmet(2,2)), nsc*sqrt(rmet(3,3)) )
 if (maxcpdst > maxatdst) maxcpdst=0.75*maxatdst


!RPRIM ITS INVERSE AND TRANSPOSE

 do ii=1,3
   do jj=1,3
     yy(ii,jj)=rprimd(ii,jj)
   end do
 end do
 call ludcmp(yy,3,3,ipiv,id,info)
 ABI_CHECK(info==0,'Error inverting rprimd')

 do  ii=1,3
   do jj=1,3
     ivrprim(ii,jj)=0._dp
   end do
   ivrprim(ii,ii)=1._dp
 end do
 do ii=1,3
   call lubksb(yy,3,3,ipiv,ivrprim(:,ii))
 end do
 do ii=1,3
   do jj=1,3
     trivrp(ii,jj)=ivrprim(jj,ii)
   end do
 end do

 write(std_out,'(" INVERSE OF RPRIMD: ",/,3F16.8,/,3F16.8,/,3F16.8,/)') &
& ((ivrprim(ii,jj), jj=1,3), ii=1,3)
 write(untout,'(" INVERSE OF RPRIMD: ",/,3F16.8,/,3F16.8,/,3F16.8,/)') &
& ((ivrprim(ii,jj), jj=1,3), ii=1,3)

 write(std_out,*) "ATOMS (index,at.number,Zionic,position(xcart.))"
 write(std_out,*) "======================================="
 do ii=1,natom
   jj=typat(ii)
   write(std_out,'(I4,2F10.6,3F16.8)') ii, znucl(jj), zionpsp(jj), (xatm(kk,ii),kk=1,3)
 end do
 write(untout,*) "ATOMS (index,at.number,Zionic,position(xcart.))"
 write(untout,*) "======================================="
 do ii=1,natom
   jj=typat(ii)
   write(untout,'(I4,2F10.6,3F16.8)') ii, znucl(jj), zionpsp(jj), (xatm(kk,ii),kk=1,3)
 end do

!STEPS IN REAL SPACE GRID (REDUCED)
 do ii=1,3
   dix(ii)=1._dp/ngfft(ii)
 end do

!READING OF THE CORE DENSITY
 write(std_out,*)ch10,' initaim : will read the core densities' ,ch10

 ABI_ALLOCATE(ndat,(ntypat))
 ABI_ALLOCATE(rminl,(natom))
 ndtmax=0
 if(me==master)then
   do ii=1,ntypat
     unth=unt+ii
!    DEBUG
!    write(std_out,*)' read from unit ',unth
!    call flush(std_out)
!    stop
!    ENDDEBUG
     read(unth,*) ndat(ii),ss
     if (ndat(ii)>ndtmax) ndtmax=ndat(ii)
   end do
 end if
 call xmpi_bcast(ndat,master,comm,ierr)
 call xmpi_bcast(ndtmax,master,comm,ierr)
 call xmpi_bcast(ss,master,comm,ierr)

!FIELDS FOR STORING CORE DENSITY

 ABI_ALLOCATE(rrad,(ndtmax,ntypat))
 ABI_ALLOCATE(crho,(ndtmax,ntypat))
 ABI_ALLOCATE(sp2,(ndtmax,ntypat))
 ABI_ALLOCATE(sp3,(ndtmax,ntypat))
 ABI_ALLOCATE(sp4,(ndtmax,ntypat))
 ABI_ALLOCATE(corlim,(ntypat))

 sp2(:,:)=zero
 sp3(:,:)=zero
 sp4(:,:)=zero

!Reading of the core densities
 corlim(:)=0
 kod=0
 if(me==master)then
   do ii=1,ntypat
     unth=unt+ii
     do jj=1,ndat(ii)
       read(unth,*) rrad(jj,ii),crho(jj,ii),sp2(jj,ii),sp3(jj,ii)
       ! this is the integral of the core charge read in
       crho(jj,ii) = crho(jj,ii)/4._dp/pi
       if ((crho(jj,ii) < aim_rhocormin) .and. (corlim(ii)==0)) corlim(ii)=jj
       sp2(jj,ii)=sp2(jj,ii)/4._dp/pi
       sp3(jj,ii)=sp3(jj,ii)/4._dp/pi   ! ATENTION!!! in sp3 is just second derivation
     end do
     do jj=1,ndat(ii)-1
       sp4(jj,ii)=(sp3(jj+1,ii)-sp3(jj,ii))/(6._dp*(rrad(jj+1,ii)-rrad(jj,ii)))
     end do
     !
     zz = crho(1,ii) * rrad(1,ii)**2 * (rrad(2,ii)-rrad(1,ii))
     do jj=2,ndat(ii)-1
       zz = zz + crho(jj,ii) * rrad(jj,ii)**2 * (rrad(jj+1,ii)-rrad(jj-1,ii))
     end do
     zz = zz * half * 4._dp * pi
     if (corlim(ii)==0) corlim(ii)=ndat(ii)

     ! add check on zion wrt FHI .fc file
     ! compare zion to zionpsp(typat(aim_dtset%batom))
     if (abs(znucl(ii) - zz - zionpsp(ii)) > 1.e-1_dp) then
       write (std_out,*) 'error: your core charge ', zz, ' does not correspond to the correct number'
       write (std_out,*) ' of valence electrons', zionpsp(ii), ' and the nuclear charge ', znucl(ii)
       write (std_out,*) ' You have probably used a pseudopotential which has more valence electrons than the'
       write (std_out,*) ' original FHI ones. ACTION: make a .fc file with the correct core charge'
       stop
     end if

   end do
 end if
 call xmpi_bcast(rrad,master,comm,ierr)
 call xmpi_bcast(crho,master,comm,ierr)
 call xmpi_bcast(sp2,master,comm,ierr)
 call xmpi_bcast(sp3,master,comm,ierr)
 call xmpi_bcast(sp4,master,comm,ierr)
 call xmpi_bcast(corlim,master,comm,ierr)

 write(std_out,*)ch10,' initaim : the core densities have been read' ,ch10


!CORRECTION OF THE CORE DENSITY NORMALISATION
 crho(:,:)=1.0003*crho(:,:)
 sp2(:,:)=1.0003*sp2(:,:)
 sp3(:,:)=1.0003*sp3(:,:)
 sp4(:,:)=1.0003*sp4(:,:)

!FIELDS FOR INTERPOLATIONS OF THE VALENCE DENSITY

 ABI_ALLOCATE(dig1,(ngfft(1)))
 ABI_ALLOCATE(dig2,(ngfft(2)))
 ABI_ALLOCATE(dig3,(ngfft(3)))
 ABI_ALLOCATE(llg1,(ngfft(1)))
 ABI_ALLOCATE(llg2,(ngfft(2)))
 ABI_ALLOCATE(llg3,(ngfft(3)))
 ABI_ALLOCATE(cdig1,(ngfft(1)-1))
 ABI_ALLOCATE(cdig2,(ngfft(2)-1))
 ABI_ALLOCATE(cdig3,(ngfft(3)-1))
 ABI_ALLOCATE(ddx,(ngfft(1),ngfft(2),ngfft(3)))
 ABI_ALLOCATE(ddy,(ngfft(1),ngfft(2),ngfft(3)))
 ABI_ALLOCATE(ddz,(ngfft(1),ngfft(2),ngfft(3)))

!DECOMPOSITION OF THE MATRIX FOR THE DETERMINATION OF COEFFICIENTS
!FOR CUBIC SPLINE INTERPOLATION (using the periodic boundary conditions)

!MAIN DIAGONAL (aa) AND SECONDARY DIAGONAL (bb) MATRIX ELEMENTS

 nmax=ngfft(1)
 do ii=2,3
   if (ngfft(ii) > nmax) nmax=ngfft(ii)
 end do
 nullify(ptf,ptsd)
 nullify(ptd,ptc,ptp)
 aa(:)=2.0*dix(:)**2/3.0
 bb(:)=dix(:)**2/6.0

 do ii=1,3
   if(ii==1) then
     ptd=>dig1;ptc=>cdig1;ptp=>llg1
   elseif (ii==2) then
     ptd=>dig2;ptc=>cdig2;ptp=>llg2
   else
     ptd=>dig3;ptc=>cdig3;ptp=>llg3
   end if
   ptd(1)=sqrt(aa(ii))
   ptc(1)=bb(ii)/ptd(1)
   ptp(1)=ptc(1)
   do jj=2,ngfft(ii)-1
     ptd(jj)=aa(ii)-ptc(jj-1)**2
     if(ptd(jj)<zero) then
       MSG_ERROR('Matrix is not positive definite !')
     end if
     ptd(jj)=sqrt(ptd(jj))
     if (jj==ngfft(ii)-1) then
       ptc(jj)=(bb(ii)-ptp(jj-1)*ptc(jj-1))/ptd(jj)
       ptp(jj)=ptc(jj)
       exit
     end if
     ptc(jj)=bb(ii)/ptd(jj)
     ptp(jj)=-ptp(jj-1)*ptc(jj-1)/ptd(jj)
   end do
   ss=0._dp
   do jj=1,ngfft(ii)-1
     ss=ss+ptp(jj)**2
   end do
   ss=aa(ii)-ss
   if(ss<zero) then
     MSG_ERROR('Matrix is not positive definite !')
   end if
   ptd(ngfft(ii))=sqrt(ss)
   ptp(ngfft(ii))=ptd(ngfft(ii))


!  INITIALISATION OF THE SECOND DERIVATIVE FIELDS

   nn=ii+1
   if (nn>3) nn=nn-3
   mm=ii+2
   if (mm>3) mm=mm-3
   do jj=1,ngfft(nn)
     do kk=1,ngfft(mm)
!      The calcul of the second derivations on the grid
       call inspln(ii,jj,kk)
     end do
   end do
   nullify(ptd,ptc,ptp)
 end do
 nullify(ptd,ptc,ptp)

 znucl_batom=znucl(typat(aim_dtset%batom))

 ABI_DEALLOCATE(znucl)
 ABI_DEALLOCATE(zionpsp)
 ABI_DEALLOCATE(symrel)
 ABI_DEALLOCATE(tnons)

!the pointers are obsolete - to remove later

end subroutine initaim
!!***

!!****f* m_bader/inpar
!! NAME
!! inpar
!!
!! FUNCTION
!! Parser for the aim utility (shorter than the one of ABINIT)
!!
!! INPUTS
!!  This routine uses data from the defs_aimprom module
!!
!! OUTPUT
!!  instr=string of character containing the input data
!!  lenstr=actual length of the character string
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      aim
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine inpar(instr,lenstr)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: lenstr
 character(len=*),intent(out) :: instr

!Local variables ------------------------------
 character(len=1),parameter :: space=' '
 character(len=26),parameter :: uplett='ABCDEFGHIJKLMNOPQRSTUVWXYZ', lolett='abcdefghijklmnopqrstuvwxyz'
!scalars
 integer,parameter :: nline=100
 integer :: ii,inxh,inxl,ios,jj,kk,ll
 character(len=fnlen) :: line

! *********************************************************************

 lenstr=0

 do ii=1,26
   inxh=index(lolett,uplett(ii:ii))
   if (inxh > 0) then
     write(std_out,*) 'ERROR The ', uplett(ii:ii) ,' is considered come lowcase !'
     MSG_ERROR("Aborting now")
   end if
 end do
 rewind(unt0)
 do ii=1,nline
   read(unt0,'(A)',iostat=ios) line(1:fnlen)
   if (ios/=0) exit
   inxh=index(line,'#')
   if (inxh == 1) then
     cycle
   elseif (inxh > 0) then
     inxl=inxh-1
     line(inxh:inxh)=space
   else
     inxl=len_trim(line)
     if (inxl==0) cycle
   end if
   inxh=index(line(1:inxl),char(9))
   if (inxh/=0) line(inxh:inxh)=space
   do ll=1,inxl
     if (iachar(line(ll:ll)) < 32) line(ll:ll)=space
   end do
   inxh=index(line(1:inxl),'- ')
   if (inxh/=0) then
     write(std_out,*) 'ERROR sign minus with white space in input file'
     MSG_ERROR("Aborting now")
   end if
   line(1:inxl)=adjustl(line(1:inxl))
   inxl=len_trim(line(1:inxl))+1
   jj=2;kk=0
   line(1:inxl)=adjustl(line(1:inxl))
   kk=len_trim(line(1:inxl))+1
   do ll=1,inxl
     inxh=index(line(jj:kk),space)
     if ((inxh==0).or.((jj+inxh-1)==kk)) exit
     line(inxh+jj:kk)=adjustl(line(inxh+jj:kk))
     kk=len_trim(line(1:inxl))
     if (kk == inxl) then
       exit
     end if
     jj=jj+inxh
   end do
   inxl=len_trim(line(1:inxl))+1
   do ll=1,inxl-1
     inxh=index(lolett,line(ll:ll))
     if (inxh/=0) line(ll:ll)=uplett(inxh:inxh)
   end do
   if ((lenstr+inxl) > strlen ) then
     write(std_out,*) 'ERROR Too large input !'
     MSG_ERROR("Aborting now")
   else
     instr(lenstr+1:lenstr+inxl)=line(1:inxl)
     lenstr=lenstr+inxl
   end if
 end do
end subroutine inpar
!!***

!!****f* m_bader/inspln
!! NAME
!! inspln
!!
!! FUNCTION
!! This procedure gives the values of the spline coefficients
!! (second derivatives) in the 1D grid with periodic boundary
!! conditions at rsid - the values of the unknown functions specified
!! in the vector valf of direction idir
!!
!! INPUTS
!!  idir= direction following which the derivatives are evaluated
!!  snn, tnn=remaining bi-dimensional coordinates of the line along which
!!        the derivative is to be computed
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  This routine works on the data contained in the aimfields module
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine inspln(idir,snn,tnn)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,snn,tnn

!Local variables-------------------------------
!scalars
 integer :: dim,ii
 real(dp) :: ss
!arrays
 real(dp) :: rsid(ngfft(idir)),valf(ngfft(idir))
 real(dp),pointer :: ptc(:),ptd(:),ptp(:)

! *************************************************************************

!POINTER INITIALIZATION

 if (idir==1) then
   valf(:)=dvl(:,snn,tnn)
 elseif (idir==2) then
   valf(:)=dvl(tnn,:,snn)
 else
   valf(:)=dvl(snn,tnn,:)
 end if

 nullify(ptd,ptc,ptp)
 if(idir==1) then
   ptd=>dig1;ptc=>cdig1;ptp=>llg1
 elseif (idir==2) then
   ptd=>dig2;ptc=>cdig2;ptp=>llg2
 else
   ptd=>dig3;ptc=>cdig3;ptp=>llg3
 end if

 dim=ngfft(idir)

!FIRST CYCLE OF RECURRENCE

 rsid(1)=valf(2)+valf(dim)-2.*valf(1)
 rsid(1)=rsid(1)/ptd(1)
 do ii=2,dim-1
   rsid(ii)=valf(ii+1)+valf(ii-1)-2.*valf(ii)
   rsid(ii)=(rsid(ii)-ptc(ii-1)*rsid(ii-1))/ptd(ii)
 end do
 ss=0._dp
 do ii=1,dim-1
   ss=ss+rsid(ii)*ptp(ii)
 end do
 rsid(dim)=valf(1)+valf(dim-1)-2.*valf(dim)
 rsid(dim)=(rsid(dim)-ss)/ptd(dim)

!SECOND CYCLE WITH TRANSPOSED MATRIX

 rsid(dim)=rsid(dim)/ptd(dim)
 rsid(dim-1)=(rsid(dim-1)-ptc(dim-1)*rsid(dim))/ptd(dim-1)
 do ii=dim-2,1,-1
   rsid(ii)=(rsid(ii)-ptc(ii)*rsid(ii+1)-ptp(ii)*rsid(dim))/ptd(ii)
 end do

 if (idir==1) then
   ddx(:,snn,tnn)=rsid(:)
 elseif (idir==2) then
   ddy(tnn,:,snn)=rsid(:)
 else
   ddz(snn,tnn,:)=rsid(:)
 end if

end subroutine inspln
!!***

!!****f* m_bader/integrho
!! NAME
!! integrho
!!
!! FUNCTION
!! This routine integrates the electron density inside the
!! atomic surface already calculated - it reads the file *.surf
!! The radial integration is always performed with splines and
!! the two angular integrations with Gauss quadrature
!!
!! INPUTS
!! aim_dtset = the structured entity containing all input variables
!! znucl_batom=the nuclear charge of the Bader atom
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  This routine works primarily on the data contained in the aimfields and aimprom modules
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine integrho(aim_dtset,znucl_batom)

!Arguments ------------------------------------
!scalars
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: batom,chs,iat,ii,inx,inxf,ipos,jj,kk,ll,nn,nph,nth
 real(dp) :: chg,chgint,cintr,ct1,ct2,lder,nsphe,phimax,phimin,rder
 real(dp) :: rsmax,rsmin,ss,stp,themax,themin,uu
 real(dp) :: znucl_batom,zz
 logical :: gaus,weit
!arrays
 real(dp) :: grho(3),hrho(3,3),shift(3),unvec(3),vv(3)
 real(dp),allocatable :: ncrho(:),nsp2(:),nsp3(:),nsp4(:),rdint(:,:),rr(:)
 real(dp),allocatable :: vdd(:),vrho(:),wgrs(:,:),work(:)

! *********************************************************************

 gaus=.true.
 weit=.true.

 write(std_out,*) 'npt = ',aim_dtset%npt

 rewind(unts)
 read(unts,*) batom,shift  ! Warning : batom is read, instead of coming from aim_dtset
 read(unts,*) nth,themin,themax ! Warning : these numbers are read, instead of coming from aim_dtset
 read(unts,*) nph,phimin,phimax ! Warning : these numbers are read, instead of coming from aim_dtset

 write(std_out,*) 'NTH NPH ',nth,nph

 ABI_ALLOCATE(wgrs,(nth,nph))
 ABI_ALLOCATE(rdint,(nth,nph))

 do ii=1,nth
   do jj=1,nph
     if (weit) then
       read(unts,*) th(ii),ph(jj),rs(ii,jj),wgrs(ii,jj)
     else
       read(unts,*) th(ii),ph(jj),rs(ii,jj)
     end if
   end do
 end do
 read(unts,*) rsmin,rsmax


 if (gaus) then
   ct1=cos(themin)
   ct2=cos(themax)
   call coeffs_gausslegint(ct1,ct2,cth,wcth,nth)
   call coeffs_gausslegint(phimin,phimax,ph,wph,nph)
 end if

 do ii=1,nth
   do jj=1,nph
     if (.not.weit) then
       if (gaus) then
         wgrs(ii,jj)=wcth(ii)*wph(jj)
       else
         wgrs(ii,jj)=1._dp
       end if
     end if
   end do
 end do


 do ii=1,nth
   do jj=1,nph
     if (rs(ii,jj) < rsmin) rsmin=rs(ii,jj)
   end do
 end do


!INTEGRATION OF THE CORE DENSITY

 nn=typat(batom)
 kk=ndat(nn)


!spherical integration of the core density in the sphere
!of the minimal Bader radius

!COEF. FOR SPHERICAL INTEGRATION

 ABI_ALLOCATE(nsp2,(kk))
 ABI_ALLOCATE(nsp3,(kk))
 ABI_ALLOCATE(nsp4,(kk))
 ABI_ALLOCATE(ncrho,(kk))

 do ii=1,kk
   ncrho(ii)=crho(ii,nn)*4._dp*pi*rrad(ii,nn)*rrad(ii,nn)
   nsp3(ii)=4._dp*pi*(2._dp*crho(ii,nn)+2._dp*rrad(ii,nn)*sp2(ii,nn)+&
&   rrad(ii,nn)*rrad(ii,nn)*sp3(ii,nn))
 end do

 if (rsmin < rrad(ndat(nn),nn)) then        ! search index
   inx=0
   if (rsmin < rrad(1,nn)) then
     MSG_ERROR('absurd')
   elseif (rsmin > rrad(ndat(nn),nn)) then
     inx=ndat(nn)
   else
     do while (rsmin >= rrad(inx+1,nn))
       inx=inx+1
     end do
   end if
 else
   inx=ndat(nn)
 end if

 cintr=4._dp/3._dp*pi*rrad(1,nn)**3*crho(1,nn)

!spline integration

 do ii=1,inx-1
   uu=rrad(ii+1,nn)-rrad(ii,nn)
   cintr=cintr+(ncrho(ii)+ncrho(ii+1))*uu/2._dp-uu*uu*uu/2.4d1*(nsp3(ii)+nsp3(ii+1))
 end do
 if (inx/=ndat(nn)) then
   uu=rsmin-rrad(inx,nn)
   zz=rrad(inx+1,nn)-rsmin
   ss=rrad(inx+1,nn)-rrad(inx,nn)
   cintr=cintr+ncrho(inx)/2._dp*(ss-zz*zz/ss)+ncrho(inx+1)/2._dp*uu*uu/ss+&
   nsp3(inx)/1.2d1*(zz*zz*ss-zz*zz*zz*zz/2._dp/ss-ss*ss*ss/2._dp)+&
   nsp3(inx+1)/1.2d1*(uu*uu*uu*uu/2._dp/ss-uu*uu*ss)
 end if


!INTEGRATION OF THE REST OF THE CORE DENSITY
!(for gauss quadrature)
!For the Gauss quadrature it is added
!to the radial integrated valence density

 rdint(:,:)=0._dp
 nsphe=0._dp
 do ii=1,nth
   do jj=1,nph
     if (inx==ndat(nn)) cycle
     inxf=inx
     if (rs(ii,jj) < rsmin) then
       write(std_out,*) rs(ii,jj),rsmin
       MSG_ERROR('in surface')
     elseif (rs(ii,jj) > rrad(ndat(nn),nn)) then
       inxf=ndat(nn)
     else
       do while (rs(ii,jj) >= rrad(inxf+1,nn))
         inxf=inxf+1
       end do
     end if

     if (inxf==inx) then
       uu=rrad(inx+1,nn)-rs(ii,jj)
       zz=rrad(inx+1,nn)-rsmin
       ss=rrad(inx+1,nn)-rrad(inx,nn)

       rdint(ii,jj)=(ncrho(inx)/2._dp/ss-nsp3(inx)/1.2d1*ss)*(zz*zz-uu*uu)+&
       nsp3(inx)/2.4d1/ss*(zz**4-uu**4)
       uu=rs(ii,jj)-rrad(inx,nn)
       zz=rsmin-rrad(inx,nn)
       rdint(ii,jj)=rdint(ii,jj)+(uu*uu-zz*zz)*(ncrho(inx+1)/2._dp/ss-nsp3(inx+1)/1.2d1*ss)+&
       nsp3(inx+1)/2.4d1/ss*(uu**4-zz**4)
     else
       uu=rrad(inx+1,nn)-rsmin
       zz=rsmin-rrad(inx,nn)

       rdint(ii,jj)=ncrho(inx)/2._dp/ss*uu*uu+ncrho(inx+1)/2._dp*(ss-zz*zz/ss)+&
       nsp3(inx)/1.2d1*(uu**4/2._dp/ss-uu*uu*ss)+nsp3(inx+1)/1.2d1*(zz*zz*ss-ss**3/2._dp-zz**4/2._dp/ss)
       if (inxf > inx+1) then
         do kk=inx+1,inxf-1
           uu=rrad(kk+1,nn)-rrad(kk,nn)
           rdint(ii,jj)=rdint(ii,jj)+(ncrho(kk)+ncrho(kk+1))*uu/2._dp-uu*uu*uu/2.4d1*(nsp3(kk)+nsp3(kk+1))
         end do
       end if

       if (inxf/=ndat(nn)) then
         uu=rs(ii,jj)-rrad(inxf,nn)
         zz=rrad(inxf+1,nn)-rs(ii,jj)
         ss=rrad(inxf+1,nn)-rrad(inxf,nn)
         rdint(ii,jj)=rdint(ii,jj)+ncrho(inxf)/2._dp*(ss-zz*zz/ss)+ncrho(inxf+1)/2._dp*uu*uu/ss+&
         nsp3(inxf)/1.2d1*(zz*zz*ss-zz*zz*zz*zz/2._dp/ss-ss*ss*ss/2._dp)+&
         nsp3(inxf+1)/1.2d1*(uu*uu*uu*uu/2._dp/ss-uu*uu*ss)
       end if
     end if
     rdint(ii,jj)=rdint(ii,jj)/4._dp/pi
     nsphe=nsphe+rdint(ii,jj)*wgrs(ii,jj)
   end do
 end do
 nsphe=nsphe*(pi/(themin-themax))*(two_pi/(phimax-phimin))

 write(untout,*)
 write(untout,*) "CHARGE INTEGRATION"
 write(untout,*) "=================="
 write(untout,'(" Core density contribution: ",/,/,"    ",F16.8)') cintr+nsphe

 write(std_out,*) ':INTECOR ', cintr+nsphe

 ABI_DEALLOCATE(ncrho)
 ABI_DEALLOCATE(nsp2)
 ABI_DEALLOCATE(nsp3)
 ABI_DEALLOCATE(nsp4)

!INTEGRATION OF THE VALENCE DENSITY

 ABI_ALLOCATE(rr,(aim_dtset%npt+1))
 ABI_ALLOCATE(vrho,(aim_dtset%npt+1))
 ABI_ALLOCATE(vdd,(aim_dtset%npt+1))

!in the case of the only irho appelation

 nn=0
 do ii=-3,3
   do jj=-3,3
     do kk=-3,3
       nn=nn+1
       atp(1,nn)=ii*1._dp
       atp(2,nn)=jj*1._dp
       atp(3,nn)=kk*1._dp
       call bschg1(atp(:,nn),1)
       if ((ii==0).and.(jj==0).and.(kk==0)) ipos=nn
     end do
   end do
 end do
 nnpos=nn
 iat=batom

!XG020629 There is a problem with this routine
!(or vgh_rho), when one uses the PGI compiler :
!The following line is needed, otherwise, iat and ipos
!are set to 0 inside vgh_now. Why ????
 write(std_out,*)' integrho : iat,ipos=',iat,ipos
!

 nsphe=0._dp
 ABI_ALLOCATE(work,(aim_dtset%npt+1))
 do ii=1,nth
   do jj=1,nph

     stp=rs(ii,jj)/aim_dtset%npt
     unvec(1)=sin(th(ii))*cos(ph(jj))
     unvec(2)=sin(th(ii))*sin(ph(jj))
     unvec(3)=cos(th(ii))
     do kk=0,aim_dtset%npt
       rr(kk+1)=kk*stp
       vv(:)=xatm(:,batom)+kk*stp*unvec(:)
       chs=-2
       call vgh_rho(vv,chg,grho,hrho,uu,iat,ipos,chs)
       vrho(kk+1)=chg*rr(kk+1)*rr(kk+1)
       if (kk==aim_dtset%npt) then
         rder=0._dp
         do ll=1,3
           rder=rder+grho(ll)*unvec(ll)
         end do
         rder=rder*rr(kk+1)*rr(kk+1)+2._dp*rr(kk+1)*chg
       end if
     end do
     lder=0._dp
     kk=aim_dtset%npt+1
     call spline(rr,vrho,kk,lder,rder,vdd)

!    INTEGRATION

     do kk=1,aim_dtset%npt
       rdint(ii,jj)=rdint(ii,jj)+stp/2._dp*(vrho(kk)+vrho(kk+1))&
&       -stp*stp*stp/24._dp*(vdd(kk)+vdd(kk+1))
     end do
     nsphe=nsphe+rdint(ii,jj)*wgrs(ii,jj)
   end do
 end do
 ABI_DEALLOCATE(work)

 if (gaus.or.weit) then
   nsphe=nsphe*(pi/(themin-themax))*(two_pi/(phimax-phimin))
 else
   nsphe=nsphe/(nth*nph)*2.0*two_pi
 end if
 chgint=cintr+nsphe

 write(untout,'(/," Different density contributions: Core (only spherical part) and the rest ",/,/,"      ",2F16.8)') &
& cintr, nsphe
 write(untout,'(/,a,i4,a,f14.8)') ' For atom number ',batom,', the number of electrons in the Bader volume is ',chgint
 write(untout,'(a,f15.7,a,f17.8)') ' The nuclear charge is',znucl_batom,', so that the Bader charge is ',znucl_batom-chgint
 write(untout,*)
 write(std_out,*) ':INTEPAR ', cintr, nsphe
 write(std_out,*) ':RHOTOT ',batom,chgint

end subroutine integrho
!!***

!!****f* m_bader/integvol
!! NAME
!! integvol
!!
!! FUNCTION
!! This routine integrates the volume of the Bader atom
!!
!! INPUTS
!!  (see side effects)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  This routine works on the data contained in the aimfields and aimprom modules
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine integvol()

!Arguments ------------------------------------

!Local variables ------------------------------
!scalars
 integer :: batom,ii,jj,nph,nth
 real(dp) :: chgint,ct1,ct2,nsphe,phimax,phimin
 real(dp) :: rsmax,rsmin,themax,themin
 logical :: gaus,weit
!arrays
 real(dp) :: shift(3)
 real(dp),allocatable :: rdint(:,:)
 real(dp),allocatable :: wgrs(:,:)

! *********************************************************************

 tpi=two_pi
 gaus=.true.
 weit=.true.


 rewind(unts)
 read(unts,*) batom,shift
 read(unts,*) nth,themin,themax
 read(unts,*) nph,phimin,phimax

 write(std_out,*) 'NTH NPH ',nth,nph

 ABI_ALLOCATE(wgrs,(nth,nph))
 ABI_ALLOCATE(rdint,(nth,nph))

 do ii=1,nth
   do jj=1,nph
     if (weit) then
       read(unts,*) th(ii),ph(jj),rs(ii,jj),wgrs(ii,jj)
     else
       read(unts,*) th(ii),ph(jj),rs(ii,jj)
     end if
   end do
 end do
 read(unts,*) rsmin,rsmax


 if (gaus) then
   ct1=cos(themin)
   ct2=cos(themax)
   call coeffs_gausslegint(ct1,ct2,cth,wcth,nth)
   call coeffs_gausslegint(phimin,phimax,ph,wph,nph)
 end if

 do ii=1,nth
   do jj=1,nph
     if (.not.weit) then
       if (gaus) then
         wgrs(ii,jj)=wcth(ii)*wph(jj)
       else
         wgrs(ii,jj)=1._dp
       end if
     end if
   end do
 end do

 nsphe=0._dp
 do ii=1,nth
   do jj=1,nph
     nsphe=nsphe+rs(ii,jj)**3/3._dp*wgrs(ii,jj)
   end do
 end do
 if (gaus.or.weit) then
   nsphe=nsphe*(pi/(themin-themax))*(tpi/(phimax-phimin))
 else
   nsphe=nsphe/(nth*nph)*2.0*tpi
 end if
 chgint=nsphe

 write(std_out,*) ':VOLTOT ',batom,chgint
 write(untout,'("Volume of the Bader atom: ", I6, F16.8)') batom,chgint

end subroutine integvol
!!***

!!****f* m_bader/onestep
!! NAME
!! onestep
!!
!! FUNCTION
!! Advance one step following the gradient from vv(3).
!! It returns a new point in vv(3) and the value and gradient of the
!! electron density at this point in chg and grho(3)
!!
!! INPUTS
!!  npmax= maximum number of divisions
!!  hh= determines the initial value of the step (to be multiplied by grho)
!!
!! OUTPUT
!!  chg= value of electron density
!!  deltar= the length of the step thaty was needed
!!  grho(3)= gradient of electron density
!!  np= returns the number of divisions that were needed
!!
!! SIDE EFFECTS
!!  vv(3)=starting and updated point
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine onestep(vv,chg,grho,hh,np,npmax,deltar)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npmax
 integer,intent(out) :: np
 real(dp),intent(in) :: hh
 real(dp),intent(out) :: chg,deltar
!arrays
 real(dp),intent(inout) :: vv(3)
 real(dp),intent(out) :: grho(3)

!Local variables ------------------------------
!scalars
 integer :: iat,ii,ipos,jj
 real(dp) :: dt,rr
!arrays
 real(dp) :: hrho(3,3),pom(3),vinter(3,200),vk(3),vkold(3)

!************************************************************************
 dt=hh
 np=1
 deltar=1._dp
 vk(1:3)=vv(1:3)


 do while((np<3).or.((np<=npmax).and.(deltar>aim_deltarmin)))
   np=np*2
   dt=dt*0.5_dp
   vkold(1:3)=vk(1:3)
   call vgh_rho(vk,chg,grho,hrho,rr,iat,ipos,0)
   vinter(1:3,1)=vv(1:3)+dt*grho(1:3)
   do jj=2,np
     call vgh_rho(vinter(1,jj-1),chg,grho,hrho,rr,iat,ipos,0)
     if(jj.eq.2) then
       vinter(1:3,2)=vv(1:3)+2.0*dt*grho(1:3)
     else
       vinter(1:3,jj)=vinter(1:3,jj-2)+2.0*dt*grho(1:3)
     end if
   end do

   call vgh_rho(vinter(1,np),chg,grho,hrho,rr,iat,ipos,0)
   vinter(1:3,np+1)=vinter(1:3,np-1)+dt*grho(1:3)

   deltar=0._dp
   do ii=1,3
     vk(ii)=(vinter(ii,np)+vinter(ii,np+1))*0.5_dp
     deltar=deltar+(vkold(ii)-vk(ii))*(vkold(ii)-vk(ii))
   end do
 end do

 pom(:)=vk(:)-vv(:)
 deltar=vnorm(pom,0)
 vv(1:3)=vk(1:3)

 call vgh_rho(vv,chg,grho,hrho,rr,iat,ipos,0)
 if(deb) write(std_out,*) ':VKf ',np,vk

end subroutine onestep
!!***

!!****f* m_bader/plint
!! NAME
!! plint
!!
!! FUNCTION
!! This simple routine gives the profile of the density
!! integrated in xy plane belong the z-axes (it works only
!! for orthogonal coordinates at present - it is better to use cut3d)
!! integration in plane - with equilateral triangles (not really
!! finished and not tested!)
!!
!! INPUTS
!!  (this routine works on the data in the aimprom module)
!!
!! OUTPUT
!!  (this routine works on the data in the aimprom module)
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine plint()

!Arguments ------------------------------------

!Local variables ------------------------------
!scalars
 integer,parameter :: nd=150,ng=300
 integer :: cod,iat,ii,ipos,jj,kk,nn
 real(dp) :: dd,ee,ff,gg,hh,igr,rho,ss
 logical :: prep
!arrays
 real(dp) :: grho(3),hrho(3,3),vv(3),xl(nd+1),xs(nd)
 real(dp),allocatable :: uu(:)

! *********************************************************************

 ff=rprimd(1,1)/nd
 ss=2._dp/sqrt(3._dp)*rprimd(2,2)/rprimd(1,1)*nd
 nn=int(ss)
 gg=sqrt(3._dp)/2.*ff
 hh=rprimd(2,2)-nn/nd*sqrt(3._dp)/2.*rprimd(1,1)
 ee=hh/sqrt(3._dp)
 hh=hh/2.
 ss=sqrt(3._dp)*ff*ff/4.
 dd=ee*ff/2.

 do ii=1,nd
   xl(ii)=ii*ff
   xs(ii)=ff/2.+ii*ff
 end do
 xl(nd+1)=rprimd(1,1)

 ABI_ALLOCATE(uu,(nn+3))

 uu(1)=0._dp
 uu(nn+3)=rprimd(2,2)
 do ii=2,nn+2
   uu(ii)=hh+(ii-1)*gg
 end do
 igr=0._dp
 prep=.true.
 do kk=1,ng
   igr=0._dp
   vv(3)=(kk-1)*rprimd(3,3)/ng
   do ii=1,nn+3
     vv(2)=uu(ii)
     do jj=1,nd
       if (prep) then
         vv(1)=xl(jj)
         prep=.false.
       else
         vv(1)=xs(jj)
         prep=.true.
       end if
       call vgh_rho(vv,rho,grho,hrho,dd,iat,ipos,cod)
       if ((ii==1).or.(ii==nn+3)) then
         igr=igr+dd*rho
       elseif ((ii==2).or.(ii==nn+2)) then
         igr=igr+(dd+ss)*rho
       else
         igr=igr+ss*2*rho
       end if
     end do
   end do
   write(untp,'(2E16.8)') vv(3), igr
 end do
 ABI_DEALLOCATE(uu)

end subroutine plint
!!***

!!****f* m_bader/rsurf
!! NAME
!! rsurf
!!
!! FUNCTION
!! Basic routine for determination of the radius of Bader surface
!! for spherical rayon theta,phi
!! the bassin is tested by following the gradient line
!! If srch==true (in general for calls from surf) the routine aim_follow
!! is called to stop when it arrives under already known part of surface
!! Simple bissection method is used to obtain the radius
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! INPUTS
!! aim_dtset= the structured entity containing all input variables
!! rr0= starting radius
!! theta,phi = the spherical direction
!! iatinit= the atom index
!! srch= see above
!! npmax= maximum number of divisions in one step for follow
!!
!! OUTPUT
!! rr= radius
!! grho(3)= gradient on the surface
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine rsurf(aim_dtset,rr,grho,theta,phi,rr0,iatinit,npmax,srch)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iatinit,npmax
 real(dp),intent(in) :: phi,rr0,theta
 real(dp),intent(out) :: rr
 logical,intent(in) :: srch
!arrays
 real(dp),intent(out) :: grho(3)
!no_abirules
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: iat,ii,ipos,iposinit,jj,nstep
 real(dp),parameter :: mfkt=1.d1
 real(dp) :: aa,dmax,dr,drr,rho,rr1,rr2,t1,t2,wall
 logical :: cross,deb_tmp,in,in1,in2,low,srch_tmp
!arrays
 real(dp) :: hrho(3,3),unvec(3),vv(3)

! *********************************************************************

 srch_tmp=srch
 deb_tmp=deb

!unity vecteur in the direction (theta,phi)

 unvec(1)=sin(theta)*cos(phi)
 unvec(2)=sin(theta)*sin(phi)
 unvec(3)=cos(theta)


 rr=rr0
 rr1=rr
 rr2=rr
 drr=1._dp
 if (abs(rr0-r0)<1.0d-12) then
   dr=aim_dtset%dr0*mfkt
 else
   dr=aim_dtset%dr0
 end if

 vv(1)=xatm(1,aim_dtset%batom)
 vv(2)=xatm(2,aim_dtset%batom)
 vv(3)=xatm(3,aim_dtset%batom)


 iposinit=batcell
 write(std_out,'("ATOM iat=",i4," ipos=",i4)') aim_dtset%batom,batcell
 jj=0

 cross=.false.

 in=.true.
 low=.false.

 dmax=h0

 in1=.true.
 in2=in1

 do while((drr>aim_drmin).or.(jj<2))
   call timein(t1,wall)
   jj=jj+1
   do ii=1,3
     vv(ii)=xatm(ii,aim_dtset%batom)+rr*unvec(ii)
   end do

!  VACUUM CONDITION

   call vgh_rho(vv,rho,grho,hrho,aa,iat,ipos,0)
   if (rho < aim_rhomin) exit

   ldeb=.false.

   call aim_follow(aim_dtset,vv,npmax,srch_tmp,iatinit,iposinit,iat,ipos,nstep)

   call timein(t2,wall)
   t2=t2-t1

   write(std_out,'(a,i4,a,f12.8,a,i4,a,i4,a,f10.5,a,i4)') &
&   ' :STEP ',jj,' r=',rr,' iat=',iat,' ipos=',ipos,' time(sec)=',t2,' nstep=',nstep

   if ((iat.eq.iatinit).and.(ipos.eq.iposinit)) then
     in=.true.
   else
     in=.false.
   end if

!
!  NEW RADIUS
!

   if ((jj.eq.1).or.((in1.eqv.in).and.(.not.cross))) then
     if (in) then
       rr2=rr1
       rr1=rr
       rr=rr+dr
     else
       rr2=rr1
       rr1=rr
       rr=rr-dr
     end if
     if ((jj>2).and.(dr<(0.6))) then
!      modification of the step
       dr=dr*aim_fac
       if (deb_tmp) write(std_out,*) ':DR ',dr
     end if
   else
     if (.not.cross) then
       cross=.true.
       rr2=rr1
     else
       if (in2) then
         if (in) then
           rr2=rr1
         else
           in1=in2
         end if
       else
         if (in) then
           in1=in2
         else
           rr2=rr1
         end if
       end if
     end if
     rr1=rr
     rr=(rr2+rr1)/2.0
   end if

   in2=in1
   in1=in
   drr=abs(rr2-rr1)/rr
   if (deb_tmp) write(std_out,*) ':DRR ',jj,rr2,rr1,drr
 end do

end subroutine rsurf
!!***

!!****f* m_bader/surf
!! NAME
!! surf
!!
!! FUNCTION
!! Determination of the Bader surface.
!! Use rsurf to determine radius for one direction
!! simple bisection method is used
!! the bassin is tested following the gradient (follow) =
!! = the most time consuming
!! follow stops if the gradient line is near the atom
!! or if it is under already known part of surface - this is why
!! the surface is not computed row by row.
!!
!! INPUTS
!! aim_dtset= the structured entity containing all input variables
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  This routine works primarily on the data contained in the defs_aimprom module
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine surf(aim_dtset)

!Arguments ------------------------------------
!scalars
 type(aim_dataset_type) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: ierr,ii,ijj,ijj_exist,incr,init,iph,iph2,ith,ith2,jj,jj_exist,kk,level,me,mm,nn,nph,npmax,nproc,nth,comm
 real(dp) :: ct1,ct2,phi,rr,rsmax,rsmin,rthe,rthe0,t1,t2,theta,tt0,vcth,vph,vth
 real(dp) :: wall,xy,xyz
 logical :: srch,stemp
!arrays
 real(dp) :: grho(3),vr(3),vv(3)
 real(dp),allocatable :: rs_computed(:,:)

!************************************************************************

 comm = xmpi_world
 me=xmpi_comm_rank(comm)
 nproc=xmpi_comm_size(comm)

 ttsrf=zero

 rewind(unts)

 nth=aim_dtset%nth
 nph=aim_dtset%nph

!Coefficients for spherical Gauss quadrature

 ct1=cos(aim_dtset%themin)
 ct2=cos(aim_dtset%themax)
 call coeffs_gausslegint(ct1,ct2,cth,wcth,nth)
 call coeffs_gausslegint(aim_dtset%phimin,aim_dtset%phimax,ph,wph,nph)

!DEBUG
!write(std_out,*)' surf : wcth=',wcth(1:nth)
!write(std_out,*)' surf : wph=',wph(1:nth)
!ENDDEBUG

 do ijj=1,nth
   th(ijj)=acos(cth(ijj))
   if (aim_dtset%isurf/=-1) then
     do jj=1,nph
       rs(ijj,jj)=zero
     end do
   end if
 end do

 npmax=aim_npmaxin
 rsmax=0.0
 rsmin=100.0
 rthe0=r0
 srch=.false.

 do ijj=1,3
   vv(ijj)=xatm(ijj,aim_dtset%batom)
 end do


 write(std_out,*)
 write(std_out,*) "BADER SURFACE DETERMINATION"
 write(std_out,*) "==========================="
 write(std_out,*)

 write(untout,*)
 write(untout,*) "BADER SURFACE DETERMINATION"
 write(untout,*) "==========================="
 write(untout,*)

 write(std_out,'(" Atom:  ",i3,3F15.10)') aim_dtset%batom,vv
 write(std_out,'(" Theta: ",i3,2F15.10)') nth,aim_dtset%themin,aim_dtset%themax
 write(std_out,'(" Phi:   ",i3,2F15.10)') nph,aim_dtset%phimin,aim_dtset%phimax

 write(untout,'(" Atom:  ",i3,3F15.10)') aim_dtset%batom,vv
 write(untout,'(" Theta: ",i3,2F15.10)') nth,aim_dtset%themin,aim_dtset%themax
 write(untout,'(" Phi:   ",i3,2F15.10)') nph,aim_dtset%phimin,aim_dtset%phimax

 write(unts,'(i3,3F15.10)') aim_dtset%batom,vv
 write(unts,'(i3,2F15.10)') nth,aim_dtset%themin,aim_dtset%themax
 write(unts,'(i3,2F15.10)') nph,aim_dtset%phimin,aim_dtset%phimax

!write(std_out,*) 'npmax in surf= ',npmax

 ith=0
 iph=0
 tt0=0._dp
 call timein(tt0,wall)

 write(untout,*)
 write(untout,*) "DEVELOPMENT OF THE RADII DETERMINATIONS"
 write(untout,*) "========================================"
 write(untout,*)
 write(untout,*) "Determination near the CPs:"

!Determination of the CP neighbouring radii

 if (aim_dtset%isurf/=-1) then

!  Precomputation of the value of the radii (for parallelisation)
!  To make the output independent of the number of processors, but still
!  cut down the CPU time, use a multigrid technique
   srch=.true.
   ABI_ALLOCATE(rs_computed,(nth,nph))
   rs(:,:)=zero
   rs_computed(:,:)=zero
   kk=0 ; init=0
   do level=3,0,-1
     incr=2**level
     if(incr<nth .and. incr<nph)then
       rs_computed(:,:)=rs(1:nth,1:nph)
       rs(1:nth,1:nph)=zero
       do ijj=1,nth,incr
         do jj=1,nph,incr
           if(rs_computed(ijj,jj)<1.0d-12) then
             kk=kk+1
             if(mod(kk,nproc)==me)then
!              Find an approximate starting radius, from the already computed ones
               if(init==0)then
                 rthe=r0
               else
                 ijj_exist=ijj ; if(mod(ijj-1,2*incr)>=incr)ijj_exist=ijj-incr
                 jj_exist=jj ; if(mod(jj-1,2*incr)>=incr)jj_exist=jj-incr
                 rthe=rs_computed(ijj_exist,jj_exist)
                 if(rthe<1.0d-12)then
                   write(std_out,*)' surf : there is a bug ! rthe=',rthe
                   MSG_ERROR("Aborting now")
                 end if
               end if
               call timein(t1,wall) ; t2=zero
               call rsurf(aim_dtset,rr,grho,th(ijj),ph(jj),rthe,aim_dtset%batom,npmax,srch)
               rs(ijj,jj)=rr
               if (deb) then
                 call timein(t2,wall) ; t2=t2-t1
                 write(std_out,*) ':CALCULATED NP',ijj,jj,th(ijj),ph(jj),rthe,npmax,rs(ijj,jj),t2
               end if
             end if
           end if
         end do ! jj
       end do ! ijj
       call xmpi_sum(rs,comm,ierr)
!      Combine the set of already computed radii and the set of the newly computed, to obtain all computed.
       rs(1:nth,1:nph)=rs(1:nth,1:nph)+rs_computed(:,:)
       init=1
     end if
   end do
   ABI_DEALLOCATE(rs_computed)

   srch=.true.

   do ijj=1,nbcp
!    if ((icpc(ijj) == -1)) then
     rthe0=vnorm(pc(:,ijj),0)
     do jj=1,3
       vr(jj)=pc(jj,ijj)-vv(jj)+xatm(jj,aim_dtset%batom)
     end do
     xy=vr(1)*vr(1)+vr(2)*vr(2)
     xyz=xy+vr(3)*vr(3)
     xyz=sqrt(xyz)

     if (xy < aim_xymin) then
       vcth=1._dp
       if (vr(3) < 0._dp) vcth=-vcth
       vph=0._dp
     else
       vcth=vr(3)/xyz
       vph=atan2(vr(2),vr(1))
     end if

     vth=acos(vcth)
     write(untout,'(/," BCP: (index,theta,phi)",I4,2E16.8)') ijj,vth,vph

     if (vth < th(1)) then
       ith=0
     else
       if (vth > th(nth)) then
         ith=nth
       else
         do ii=2,nth
           if (vth < th(ii)) then
             ith=ii-1
             exit
           end if
         end do
       end if
     end if

     if (vph < ph(1)) then
       iph=0
     else
       if (vph > ph(nph)) then
         iph=nph
       else
         do ii=2,nph
           if (vph < ph(ii)) then
             iph=ii-1
             exit
           end if
         end do
       end if
     end if

     write(untout,*) "ATOMIC RADII (ith,iph,theta,phi,radius)"
     do jj=-1,2
       do kk=-1,2
         ith2=ith+jj
         iph2=iph+kk
         stemp=(iph2 > 0).and.(iph2 < nph+1)
         stemp=(stemp.and.((ith2 > 0).and.(ith2 < nth+1)))
         if (stemp) then
           theta=th(ith2)
           phi=ph(iph2)
           if (abs(rs(ith2,iph2))<1.0d-12) then
             rthe=rthe0
             if (deb) write(std_out,*) ':CALCULATING NP',theta,phi,rthe,npmax
             call timein(t1,wall)
             call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
             call timein(t2,wall)
             t2=t2-t1
             rs(ith2,iph2)=rr
           end if
           rr=rs(ith2,iph2)
!          write(unts,'(2F12.8,2E16.8)') theta,phi,rr,wcth(ijj)*wph(jj)
           write(std_out,'(":RSUR PC ",3i3,4E16.8,F10.4)') ijj,jj,kk,theta,phi,rr,wcth(ith2)*wph(iph2),t2
           write(untout,'(a,2i3,3E16.8)') '-  ',jj,kk,theta,phi,rr
           rthe0=rr
         end if

       end do ! kk
     end do ! jj

!    end if

   end do ! ijj (loop on BCP)

!  DEBUG
!  write(std_out,*)' surf : near BCP '
!  do ijj=1,nth
!  do jj=1,nph
!  write(std_out,*)ijj,jj,rs(ijj,jj)
!  end do
!  end do
!  ENDDEBUG


   srch=.true.
   do ijj=nbcp+1,nbcp+nrcp     ! Loop on RCP
!    if ((icpc(ijj) == 1)) then
     rthe0=max(rminl(aim_dtset%batom),r0)
     do jj=1,3
       vr(jj)=pc(jj,ijj)-vv(jj)+xatm(jj,aim_dtset%batom)
     end do
     xy=vr(1)*vr(1)+vr(2)*vr(2)
     xyz=xy+vr(3)*vr(3)
     xyz=sqrt(xyz)

     if (xy < aim_xymin) then
       vcth=1._dp
       if (vr(3) < 0._dp) vcth=-vcth
       vph=0._dp
     else
       vcth=vr(3)/xyz
       vph=atan2(vr(2),vr(1))
     end if
     vth=acos(vcth)
     write(untout,'(/,";RCP: (index,theta,phi)",I4,2E16.8)') ijj-nbcp,vth,vph

     if (vth < th(1)) then
       ith=0
     else
       if (vth > th(nth)) then
         ith=nth
       else
         do ii=2,nth
           if (vth < th(ii)) then
             ith=ii-1
             exit
           end if
         end do
       end if
     end if

     if (vph < ph(1)) then
       iph=0
     else
       if (vph > ph(nph)) then
         iph=nph
       else
         do ii=2,nph
           if (vph < ph(ii)) then
             iph=ii-1
             exit
           end if
         end do
       end if
     end if

     write(untout,*) "ATOMIC RADIUS (ith,iph,theta,phi,radius)"
     do jj=-1,2
       do kk=-1,2
         ith2=ith+jj
         iph2=iph+kk
         stemp=(iph2 > 0).and.(iph2 < nph+1)
         stemp=stemp.and.(ith2 > 0).and.(ith2 < nth+1)

         if (stemp) then
           theta=th(ith2)
           phi=ph(iph2)
           if ((abs(rs(ith2,iph2))<1.0d-12)) then
             rthe=rthe0
             if (deb) write(std_out,*) ':CALCULATING NP',theta,phi,rthe,npmax
             call timein(t1,wall)
             call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
             call timein(t2,wall)
             t2=t2-t1
             rs(ith2,iph2)=rr
           end if
           rr=rs(ith2,iph2)
!          write(unts,'(2F12.8,2E16.8)') theta,phi,rr,wcth(ijj)*wph(jj)
           write(std_out,'(":RSUR PC ",3i3,4E16.8,F10.4)') ijj,jj,kk,theta,phi,rr,wcth(ith2)*wph(iph2),t2
           write(untout,'(a,2i3,3E16.8)') '-  ',jj,kk,theta,phi,rr
           rthe0=rr
         end if

       end do ! kk
     end do ! jj
!    end if

   end do ! ijj (Loop on RCP)

!  DEBUG
!  write(std_out,*)' surf : near RCP '
!  do ijj=1,nth
!  do jj=1,nph
!  write(std_out,*)ijj,jj,rs(ijj,jj)
!  end do
!  end do
!  ENDDEBUG

!  Boundary angles
   rthe0=r0
   srch=.true.
   write(untout,*)
   write(untout,*) "The boundary angles:"
   write(untout,*) "===================="
   write(untout,*) "ATOMIC RADIUS (ith,iph,theta,phi,radius)"

!  Must have sufficient angular sampling
   if ((nth > 8).and.(nph > 8)) then
     rthe=r0
     do ijj=1,2
       theta=th(ijj)
       if (ijj==2) rthe=rs(1,1)
       do jj=1,nph
         phi=ph(jj)
         call timein(t1,wall)
         if (abs(rs(ijj,jj))<1.0d-12) then
           if (deb) write(std_out,*) ':CALC NP',theta,phi,rthe,npmax
           call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
           rs(ijj,jj)=rr
         end if
         rr=rs(ijj,jj)
         call timein(t2,wall)
         t2=t2-t1
         write(std_out,'(":RSUR ",2F12.8,2E16.8,F10.4)') theta,phi,rr,wcth(ijj)*wph(jj),t2
         write(untout,'(a,2i3,3E16.8)') '-  ',ijj,jj,theta,phi,rr
         rthe=rs(ijj,jj)
       end do ! jj
     end do ! ijj

     write(untout,*)

     rthe=rs(2,1)
     do jj=1,2
       phi=ph(jj)
       if (jj==2) rthe=rs(2,2)
       do ijj=3,nth
         theta=th(ijj)
         t2=0.0
         call timein(t1,wall)
         if (abs(rs(ijj,jj))<1.0d-12) then
           if (deb) write(std_out,*) ':CALC NP',theta,phi,rthe,npmax
           call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
           rs(ijj,jj)=rr
         end if
         rr=rs(ijj,jj)
         call timein(t2,wall)
         t2=t2-t1
         write(std_out,'(":RSUR ",2F12.8,2E16.8,F10.4)') theta,phi,rr,wcth(ijj)*wph(jj),t2
         write(untout,'(2i3,3E16.8)') ijj,jj,theta,phi,rr
         rthe=rs(ijj,jj)
       end do ! ijj
     end do ! jj

     write(untout,*)

     rthe=rs(nth-1,2)
     do ijj=nth-1,nth
       theta=th(ijj)
       if (ijj==nth) rthe=rs(nth,2)
       do jj=3,nph
         phi=ph(jj)
         call timein(t1,wall)
         if (abs(rs(ijj,jj))<1.0d-12) then
           if (deb) write(std_out,*) ':CALC NP',theta,phi,rthe,npmax
           call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
           rs(ijj,jj)=rr
         end if
         rr=rs(ijj,jj)
         call timein(t2,wall)
         t2=t2-t1
         write(std_out,'(":RSUR ",2F12.8,2E16.8,F10.4)') theta,phi,rr,wcth(ijj)*wph(jj),t2
         write(untout,'(2i3,3E16.8)') ijj,jj,theta,phi,rr
         rthe=rs(ijj,jj)
       end do ! jj
     end do ! ijj

     rthe=rs(2,nph-1)
     do jj=nph-1,nph
       phi=ph(jj)
       if (jj==nph) rthe=rs(2,nph)
       do ijj=3,nth-2
         theta=th(ijj)
         t2=0.0
         call timein(t1,wall)
         if (abs(rs(ijj,jj))<1.0d-12) then
           if (deb) write(std_out,*) ':CALC NP',theta,phi,rthe,npmax
           call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
           rs(ijj,jj)=rr
         end if
         rr=rs(ijj,jj)
         call timein(t2,wall)
         t2=t2-t1
         write(std_out,'(":RSUR ",2F12.8,2E16.8,F10.4)') theta,phi,rr,wcth(ijj)*wph(jj),t2
         write(untout,'(2i3,3E16.8)') ijj,jj,theta,phi,rr
         rthe=rs(ijj,jj)
       end do ! ijj
     end do ! jj
     write(untout,*)

!    Complementary bands for boundary angles
     nn=int(real(nth)/1.4d1)
     if (nn > 1) then
       do ii=1,nn-1
         mm=int(nth/nn)*ii
         do kk=0,1
           mm=mm+kk
           theta=th(mm)
           rthe=rs(mm,2)
           do jj=3,nph-2
             phi=ph(jj)
             call timein(t1,wall)
             if (abs(rs(mm,jj))<1.0d-12) then
               if (deb) write(std_out,*) ':CALC NP',theta,phi,rthe,npmax
               call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
               rs(mm,jj)=rr
             end if
             rr=rs(mm,jj)
             call timein(t2,wall)
             t2=t2-t1
             write(std_out,'(":RSUR ",2F12.8,2E16.8,F10.4)') theta,phi,rr,wcth(mm)*wph(jj),t2
             write(untout,'(2i3,3E16.8)') mm,jj,theta,phi,rr
             rthe=rs(mm,jj)
           end do ! jj
         end do ! kk
       end do ! ii
     end if ! nn>1

     write(untout,*)

     nn=nint(real(nph)/1.2d1)
     if (nn > 1) then
       do ii=1,nn-1
         mm=int(nph/nn)*ii
         do kk=0,1
           mm=mm+kk
           phi=ph(mm)
           rthe=rs(2,mm)

           do jj=3,nth-2
             theta=th(jj)
             call timein(t1,wall)
             if (abs(rs(jj,mm))<1.0d-12) then
               if (deb) write(std_out,*) ':CALC NP',theta,phi,rthe,npmax
               call rsurf(aim_dtset,rr,grho,theta,phi,rthe,aim_dtset%batom,npmax,srch)
               rs(jj,mm)=rr
             end if
             rr=rs(mm,jj)
             call timein(t2,wall)
             t2=t2-t1
             write(std_out,'(":RSUR ",2F12.8,2E16.8,F10.4)') theta,phi,rr,wcth(jj)*wph(mm),t2
             write(untout,'(2i3,3E16.8)') jj,mm,theta,phi,rr
             rthe=rs(jj,mm)
           end do ! jj

         end do ! kk
       end do ! ii
     end if  ! nn>1

   end if ! sufficient sampling to determine boundary angles

   write(untout,*)

!  DEBUG
!  write(std_out,*)' surf : after boundary angles '
!  do ijj=1,nth
!  do jj=1,nph
!  write(std_out,*)ijj,jj,rs(ijj,jj)
!  end do
!  end do
!  ENDDEBUG

!  Output the complete Bader surface

   write(untout,*) "The complete Bader surface:"
   write(untout,*) "==========================="
   write(untout,*) "ATOMIC RADIUS (ith,iph,theta,phi,radius)"
   rthe0=r0
   srch=.true.

!  Write all the values

   do ijj=1,nth
     theta=th(ijj)
     do jj=1,nph
       phi=ph(jj)
       rr=rs(ijj,jj)
       write(unts,'(2F12.8,2E16.8)') theta,phi,rr,wcth(ijj)*wph(jj)
       write(std_out,'(":RSUR ",2F12.8,2E16.8,F10.4)') theta,phi,rr,wcth(ijj)*wph(jj),t2
       write(untout,'(a,2i3,3E16.8)') '   ',ijj,jj,theta,phi,rr
       if (rr < rsmin) rsmin=rr
       if (rr> rsmax) rsmax=rr
     end do ! jj
   end do ! ijj
   write(unts,'(2F15.10)') rsmin,rsmax
   write(untout,'(/," The minimal and maximal radii:",/,/,"     ",2F15.10)') rsmin,rsmax

!  DEBUG
!  write(std_out,*)' surf : final output '
!  do ijj=1,nth
!  do jj=1,nph
!  write(std_out,*)ijj,jj,rs(ijj,jj)
!  end do
!  end do
!  ENDDEBUG

 end if ! determination of the critical surface

 call timein(ttsrf,wall)
 ttsrf=ttsrf-tt0

end subroutine surf
!!***

!!****f* m_bader/vgh_rho
!! NAME
!! vgh_rho
!!
!! FUNCTION
!! The general procedure to obtain the value, the gradient and the hessian
!! of the density of electrons in the point vv (in cart.coord).
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! INPUTS
!! vv(3)=position
!! chs  1 only valence density
!!      2 only core density
!!      0 total density
!!     -2 iat, ipos are nulify and ignored
!!     -1 iat,ipos = index of atom if vv is inside
!!         the "core sphere (rminl)", 0 otherwise
!!
!! OUTPUT
!! rho,grho(3),hrho(3,3) - density, gradient of density, hessian of density
!!                                 (cart. coord)
!! iat, ipos - index of the nearest atom (except chs < 0 see above)
!! rdmin  - the distance to the nearest atom
!!
!! SIDE EFFECTS
!!  This routine also works on the data contained in the defs_aimprom and defs_aimfields modules
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine vgh_rho(vv,rho,grho,hrho,rdmin,iat,ipos,chs)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chs
 integer,intent(inout) :: iat,ipos
 real(dp),intent(out) :: rdmin,rho
!arrays
 real(dp),intent(in) :: vv(3)
 real(dp),intent(out) :: grho(3),hrho(3,3)

!Local variables ------------------------------
!scalars
 integer :: ii,inmax,inmin,inx,jj,kk,ll,nn,oii,omm,onn
 integer :: selct
! real(dp),save :: cumul_cpu=0.0_dp,cumul_cpu_old=0.0_dp
 real(dp),save :: tcpui,tcpuo,twalli
 real(dp),save :: twallo
 real(dp) :: aa,bb,cc,cgrad1_rr_inv,coeff,dd,rr,rr2,rr_inv
 real(dp) :: rrad2_nn,rrad_nn,ss,uu,uu_inv,val,vt1,vt2,vt3,vw1,vw2
! real(dp) :: ss_inv
 real(dp) :: vw3
!arrays
 integer :: indx(3),inii(4,3)
 real(dp) :: cgrad(3),ches(3,3),cof(2,3),ddstar(6),ddu(2),grd(4)
 real(dp) :: hh(4,2),hrh(2),lder(4),pom2sq(2,3),pomsq(2,3)
 real(dp) :: rhstar(6),sqder(6,4),sqvlr(6,4),trsf(3,3),xx(3)
 real(dp),pointer :: ptddx(:,:,:),ptddy(:,:,:),ptddz(:,:,:),ptrho(:,:,:)

!************************************************************************
 tcpui=0.0_dp
 tcpuo=0.0_dp
 twalli=0.0_dp
 twallo=0.0_dp

 nullify(ptddx,ptddy,ptddz,ptrho)

 selct=chs

 if (selct/=2) then

!  call timein(tcpui,twalli)

!  TRANSFORMATION TO THE REDUCED COORD.

   xx(:)=vv(:)
   call bschg1(xx,-1)

!  call timein(tcpuo,twallo)
!  cumul_cpu=cumul_cpu+(tcpuo-tcpui)

!  REDUCTION TO THE PRIMITIVE CELL

   do ii=1,3
     if (xx(ii) >= one-tol12 ) then
       xx(ii)=xx(ii)-aint(xx(ii))
     elseif (xx(ii) < -tol12 ) then
       xx(ii)=xx(ii)-floor(xx(ii))
     end if
   end do


!  DETERMINATION OF THE INDEX IN THE GRID

   do ii=1,3
     indx(ii)=aint(xx(ii)*ngfft(ii))
     bb=(xx(ii)-indx(ii)*dix(ii))*ngfft(ii)
     if (indx(ii)==ngfft(ii)) then
       indx(ii)=1
       xx(ii)=0._dp
     else
       indx(ii)=indx(ii)+1
     end if

!    Explicit handling to avoid numeric problems

     if (bb > 1._dp+tol12 ) then
       cof(1,ii)=0._dp
       cof(2,ii)=1._dp
     elseif (bb < -tol12 ) then
       cof(1,ii)=1._dp
       cof(2,ii)=0._dp
     else
       cof(1,ii)=1._dp-bb
       cof(2,ii)=bb
     end if
   end do

!  3D INTERPOLATION OF THE VALENCE DENSITY

!  determination of the values of density and of its second derivative
!  at the "star" = constructed at vv with primitive directions
!  To interpolation the values at the faces of the grid cell are needed

   rhstar(:)=0._dp
   sqder(:,:)=0._dp
   sqvlr(:,:)=0._dp
   ddstar(:)=0._dp
   pomsq(:,:)=0._dp
   pom2sq(:,:)=0._dp

   oii=1; onn=1; omm=1
   if (indx(1)==ngfft(1)) oii=1-ngfft(1)
   if (indx(2)==ngfft(2)) onn=1-ngfft(2)
   if (indx(3)==ngfft(3)) omm=1-ngfft(3)

!  the values in the corners of the grid cell

   ptddx=>ddx(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm)
   ptddy=>ddy(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm)
   ptddz=>ddz(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm)
   ptrho=>dvl(indx(1):indx(1)+oii:oii,indx(2):indx(2)+onn:onn,indx(3):indx(3)+omm:omm)

!  the coefficients for spline interpolation of density and its derivation
   do ii=1,3
     do jj=1,2
       pomsq(jj,ii)=(cof(jj,ii)*cof(jj,ii)*cof(jj,ii)-cof(jj,ii))/6._dp*dix(ii)*dix(ii)
       pom2sq(jj,ii)=(3._dp*cof(jj,ii)*cof(jj,ii)-1._dp)/6._dp*dix(ii)
       if (jj==1) pom2sq(jj,ii)=-pom2sq(jj,ii)
     end do
   end do


   do ii=1,2
     do jj=1,2
       do kk=1,2
         ddstar(ii)=ddstar(ii)+cof(jj,2)*cof(kk,3)*ptddx(ii,jj,kk)
         ddstar(ii+2)=ddstar(ii+2)+cof(jj,3)*cof(kk,1)*ptddy(kk,ii,jj)
         ddstar(ii+4)=ddstar(ii+4)+cof(jj,1)*cof(kk,2)*ptddz(jj,kk,ii)
         sqder(ii,jj)=sqder(ii,jj)+cof(kk,2)*ptddz(ii,kk,jj)
         sqder(ii,jj+2)=sqder(ii,jj+2)+cof(kk,3)*ptddy(ii,jj,kk)
         sqder(ii+2,jj)=sqder(ii+2,jj)+cof(kk,3)*ptddx(jj,ii,kk)
         sqder(ii+2,jj+2)=sqder(ii+2,jj+2)+cof(kk,1)*ptddz(kk,ii,jj)
         sqder(ii+4,jj)=sqder(ii+4,jj)+cof(kk,1)*ptddy(kk,jj,ii)
         sqder(ii+4,jj+2)=sqder(ii+4,jj+2)+cof(kk,2)*ptddx(jj,kk,ii)
         sqvlr(ii,jj)=sqvlr(ii,jj)+cof(kk,2)*ptrho(ii,kk,jj)+pomsq(kk,2)*ptddy(ii,kk,jj)
         sqvlr(ii,jj+2)=sqvlr(ii,jj+2)+cof(kk,3)*ptrho(ii,jj,kk)+pomsq(kk,3)*ptddz(ii,jj,kk)
         sqvlr(ii+2,jj+2)=sqvlr(ii+2,jj+2)+cof(kk,1)*ptrho(kk,ii,jj)+pomsq(kk,1)*ptddx(kk,ii,jj)
       end do
     end do
   end do

   do ii=1,2
     do jj=1,2
       sqvlr(ii+2,jj)=sqvlr(jj,ii+2)
       sqvlr(ii+4,jj)=sqvlr(jj+2,ii+2)
       sqvlr(ii+4,jj+2)=sqvlr(jj,ii)
     end do
   end do

   do ii=1,2
     do jj=1,2
       rhstar(ii)=rhstar(ii)+cof(jj,3)*sqvlr(ii,jj)+pomsq(jj,3)*sqder(ii,jj)+&
&       cof(jj,2)*sqvlr(ii,jj+2)+pomsq(jj,2)*sqder(ii,jj+2)
       rhstar(ii+2)=rhstar(ii+2)+cof(jj,1)*sqvlr(ii+2,jj)+pomsq(jj,1)*sqder(ii+2,jj)+&
&       cof(jj,3)*sqvlr(ii+2,jj+2)+pomsq(jj,3)*sqder(ii+2,jj+2)
       rhstar(ii+4)=rhstar(ii+4)+cof(jj,2)*sqvlr(ii+4,jj)+pomsq(jj,2)*sqder(ii+4,jj)+&
&       cof(jj,1)*sqvlr(ii+4,jj+2)+pomsq(jj,1)*sqder(ii+4,jj+2)
     end do
   end do
   rhstar(:)=rhstar(:)/2._dp

   rho=0._dp
   grho(:)=0._dp
   hrho(:,:)=0._dp
   kk=1; nn=1
   do ii=1,5,2
     do jj=1,2
       nn=-nn
       rho=rho+cof(jj,kk)*rhstar(ii+jj-1)+pomsq(jj,kk)*ddstar(ii+jj-1)
       grho(kk)=grho(kk)+pom2sq(jj,kk)*ddstar(ii+jj-1)
       hrho(kk,kk)=hrho(kk,kk)+cof(jj,kk)*ddstar(ii+jj-1)
       grho(kk)=grho(kk)+nn*rhstar(ii+jj-1)/dix(kk)
     end do
     kk=kk+1
   end do
   rho=rho/3._dp

!  Off-diagonal elements of the hessian

!  for the speed reasons the polynomial interpolation
!  for second derivation fields is used in this case
!  but the last step is always done by spline interpolation.


   do ii=1,3
     do jj=-1,2
       inii(jj+2,ii)=indx(ii)+jj
       if (inii(jj+2,ii) < 1) inii(jj+2,ii)=inii(jj+2,ii)+ngfft(ii)
       if (inii(jj+2,ii) > ngfft(ii)) inii(jj+2,ii)=inii(jj+2,ii)-ngfft(ii)
     end do
   end do

!  Not very nice

   do ii=1,3
     select case (ii)
     case (1)
       do jj=1,4
         ddu(1)=cof(1,2)*ddz(inii(jj,1),inii(2,2),inii(2,3))+cof(2,2)*ddz(inii(jj,1),inii(3,2),inii(2,3))
         ddu(2)=cof(1,2)*ddz(inii(jj,1),inii(2,2),inii(3,3))+cof(2,2)*ddz(inii(jj,1),inii(3,2),inii(3,3))
         hrh(1)=cof(1,2)*dvl(inii(jj,1),inii(2,2),inii(2,3))+cof(2,2)*dvl(inii(jj,1),inii(3,2),inii(2,3))+&
&         pomsq(1,2)*ddy(inii(jj,1),inii(2,2),inii(2,3))+pomsq(2,2)*ddy(inii(jj,1),inii(3,2),inii(2,3))
         hrh(2)=cof(1,2)*dvl(inii(jj,1),inii(2,2),inii(3,3))+cof(2,2)*dvl(inii(jj,1),inii(3,2),inii(3,3))+&
&         pomsq(1,2)*ddy(inii(jj,1),inii(2,2),inii(3,3))+pomsq(2,2)*ddy(inii(jj,1),inii(3,2),inii(3,3))
         hh(jj,2)=(hrh(2)-hrh(1))/dix(3)+pom2sq(1,3)*ddu(1)+pom2sq(2,3)*ddu(2)

         ddu(1)=cof(1,3)*ddy(inii(jj,1),inii(2,2),inii(2,3))+cof(2,3)*ddy(inii(jj,1),inii(2,2),inii(3,3))
         ddu(2)=cof(1,3)*ddy(inii(jj,1),inii(3,2),inii(2,3))+cof(2,3)*ddy(inii(jj,1),inii(3,2),inii(3,3))
         hrh(1)=cof(1,3)*dvl(inii(jj,1),inii(2,2),inii(2,3))+cof(2,3)*dvl(inii(jj,1),inii(2,2),inii(3,3))+&
&         pomsq(1,3)*ddz(inii(jj,1),inii(2,2),inii(2,3))+pomsq(2,3)*ddz(inii(jj,1),inii(2,2),inii(3,3))
         hrh(2)=cof(1,3)*dvl(inii(jj,1),inii(3,2),inii(2,3))+cof(2,3)*dvl(inii(jj,1),inii(3,2),inii(3,3))+&
&         pomsq(1,3)*ddz(inii(jj,1),inii(3,2),inii(2,3))+pomsq(2,3)*ddz(inii(jj,1),inii(3,2),inii(3,3))
         hh(jj,1)=(hrh(2)-hrh(1))/dix(2)+pom2sq(1,2)*ddu(1)+pom2sq(2,2)*ddu(2)
       end do
     case (2)
       do jj=1,4
         ddu(1)=cof(1,3)*ddx(inii(2,1),inii(jj,2),inii(2,3))+cof(2,3)*ddx(inii(2,1),inii(jj,2),inii(3,3))
         ddu(2)=cof(1,3)*ddx(inii(3,1),inii(jj,2),inii(2,3))+cof(2,3)*ddx(inii(3,1),inii(jj,2),inii(3,3))
         hrh(1)=cof(1,3)*dvl(inii(2,1),inii(jj,2),inii(2,3))+cof(2,3)*dvl(inii(2,1),inii(jj,2),inii(3,3))+&
&         pomsq(1,3)*ddz(inii(2,1),inii(jj,2),inii(2,3))+pomsq(2,3)*ddz(inii(2,1),inii(jj,2),inii(3,3))
         hrh(2)=cof(1,3)*dvl(inii(3,1),inii(jj,2),inii(2,3))+cof(2,3)*dvl(inii(3,1),inii(jj,2),inii(3,3))+&
&         pomsq(1,3)*ddz(inii(3,1),inii(jj,2),inii(2,3))+pomsq(2,3)*ddz(inii(3,1),inii(jj,2),inii(3,3))
         hh(jj,2)=(hrh(2)-hrh(1))/dix(1)+pom2sq(1,1)*ddu(1)+pom2sq(2,1)*ddu(2)

         ddu(1)=cof(1,1)*ddz(inii(2,1),inii(jj,2),inii(2,3))+cof(2,1)*ddz(inii(3,1),inii(jj,2),inii(2,3))
         ddu(2)=cof(1,1)*ddz(inii(2,1),inii(jj,2),inii(3,3))+cof(2,1)*ddz(inii(3,1),inii(jj,2),inii(3,3))
         hrh(1)=cof(1,1)*dvl(inii(2,1),inii(jj,2),inii(2,3))+cof(2,1)*dvl(inii(3,1),inii(jj,2),inii(2,3))+&
&         pomsq(1,1)*ddx(inii(2,1),inii(jj,2),inii(2,3))+pomsq(2,1)*ddx(inii(3,1),inii(jj,2),inii(2,3))
         hrh(2)=cof(1,1)*dvl(inii(2,1),inii(jj,2),inii(3,3))+cof(2,1)*dvl(inii(3,1),inii(jj,2),inii(3,3))+&
&         pomsq(1,1)*ddx(inii(2,1),inii(jj,2),inii(3,3))+pomsq(2,1)*ddx(inii(3,1),inii(jj,2),inii(3,3))
         hh(jj,1)=(hrh(2)-hrh(1))/dix(3)+pom2sq(1,3)*ddu(1)+pom2sq(2,3)*ddu(2)
       end do
     case (3)
       do jj=1,4
         ddu(1)=cof(1,1)*ddy(inii(2,1),inii(2,2),inii(jj,3))+cof(2,1)*ddy(inii(3,1),inii(2,2),inii(jj,3))
         ddu(2)=cof(1,1)*ddy(inii(2,1),inii(3,2),inii(jj,3))+cof(2,1)*ddy(inii(3,1),inii(3,2),inii(jj,3))
         hrh(1)=cof(1,1)*dvl(inii(2,1),inii(2,2),inii(jj,3))+cof(2,1)*dvl(inii(3,1),inii(2,2),inii(jj,3))+&
&         pomsq(1,1)*ddx(inii(2,1),inii(2,2),inii(jj,3))+pomsq(2,1)*ddx(inii(3,1),inii(2,2),inii(jj,3))
         hrh(2)=cof(1,1)*dvl(inii(2,1),inii(3,2),inii(jj,3))+cof(2,1)*dvl(inii(3,1),inii(3,2),inii(jj,3))+&
&         pomsq(1,1)*ddx(inii(2,1),inii(3,2),inii(jj,3))+pomsq(2,1)*ddx(inii(3,1),inii(3,2),inii(jj,3))
         hh(jj,2)=(hrh(2)-hrh(1))/dix(2)+pom2sq(1,2)*ddu(1)+pom2sq(2,2)*ddu(2)

         ddu(1)=cof(1,2)*ddx(inii(2,1),inii(2,2),inii(jj,3))+cof(2,2)*ddx(inii(2,1),inii(3,2),inii(jj,3))
         ddu(2)=cof(1,2)*ddx(inii(3,1),inii(2,2),inii(jj,3))+cof(2,2)*ddx(inii(3,1),inii(3,2),inii(jj,3))
         hrh(1)=cof(1,2)*dvl(inii(2,1),inii(2,2),inii(jj,3))+cof(2,2)*dvl(inii(2,1),inii(3,2),inii(jj,3))+&
&         pomsq(1,2)*ddy(inii(2,1),inii(2,2),inii(jj,3))+pomsq(2,2)*ddy(inii(2,1),inii(3,2),inii(jj,3))
         hrh(2)=cof(1,2)*dvl(inii(3,1),inii(2,2),inii(jj,3))+cof(2,2)*dvl(inii(3,1),inii(3,2),inii(jj,3))+&
&         pomsq(1,2)*ddy(inii(3,1),inii(2,2),inii(jj,3))+pomsq(2,2)*ddy(inii(3,1),inii(3,2),inii(jj,3))
         hh(jj,1)=(hrh(2)-hrh(1))/dix(1)+pom2sq(1,1)*ddu(1)+pom2sq(2,1)*ddu(2)
       end do
     end select
     do jj=-2,1
       grd(jj+3)=(indx(ii)+jj)*dix(ii)
     end do

!    write(std_out,'("hh: ",/,4F16.8,/,4F16.8)') ((hh(kk,jj),kk=1,4),jj=1,2)
!    write(std_out,'("grad: ",3F16.8)') (grho(kk),kk=1,3)
!    write(std_out,'("dix: ",3F16.8)') (dix(kk),kk=1,3)
!    write(std_out,'("grd: ",4F16.8)') (grd(kk),kk=1,4)
!    write(std_out,'("inii: ",4I4)') (inii(kk,ii),kk=1,4)

     do jj=1,2

!      polynomial interpolation

       do kk=1,3
         do ll=4,kk+1,-1
           hh(ll,jj)=(hh(ll,jj)-hh(ll-1,jj))/(grd(ll)-grd(ll-1))
         end do
       end do
       lder(4)=hh(4,jj)
       do kk=3,1,-1
         lder(kk)=hh(kk,jj)+(xx(ii)-grd(kk))*lder(kk+1)
       end do
       do kk=1,2
         do ll=3,kk+1,-1
           lder(ll)=lder(ll)+(xx(ii)-grd(ll-kk))*lder(ll+1)
         end do
       end do
       nn=ii+jj
       if (nn > 3) nn=nn-3
       hrho(ii,nn)=hrho(ii,nn)+lder(2)
       hrho(nn,ii)=hrho(nn,ii)+lder(2)
     end do
   end do

!  averaging of the mixed derivations obtained in different order

   do ii=1,3
     do jj=1,3
       if (ii /= jj) hrho(ii,jj)=hrho(ii,jj)/2._dp
     end do
   end do


!  write(std_out,'("xx:",3F16.8)') (xx(ii),ii=1,3)
!  write(std_out,'("hrho: ",/,3F16.8,/,3F16.8,/,3F16.8)') &
!  & ((hrho(ii,jj),ii=1,3),jj=1,3)
!  stop
!  write(std_out,'("xx:",3F16.8)') (xx(ii),ii=1,3)
!  write(std_out,'(":GRAD pred tr ",3F16.8)') grho
!  write(std_out,'(":HESSIAN pred tr",/,3F16.8,/,3F16.8,/,3F16.8)') ((hrho(ii,jj),jj=1,3),ii=1,3)


!  Transformation back to Cart. coordonnes

   call bschg1(grho,2)
   call bschg2(hrho,2)

!  write(std_out,'("hrho: ",/,3F16.8,/,3F16.8,/,3F16.8)') &
!  & ((hrho(ii,jj),ii=1,3),jj=1,3)
!  stop

   nullify(ptddx,ptddy,ptddz,ptrho)

   if (selct==1) return

 end if

!write(51,'(":GRADv ",3F16.8)') grho
!write(52,'(":LAPv ",F16.8)') hrho(1,1)+hrho(2,2)+hrho(3,3)
!write(52,'(":HESNv ",/,3F16.8,/,3F16.8,/,3F16.8)') ((hrho(ii,jj),jj=1,3),ii=1,3)

!INTERPOLATION OF THE CORE DENSITY

 if (selct/=1) then

   if (selct==2) then
     grho(:)=0._dp
     hrho(:,:)=0._dp
     rho=0._dp
   end if

!  SEARCH OF THE NEIGHBOUR ATOMS

   if (selct /= -2) then
     iat=0
     ipos=0
   end if
   rdmin=20._dp

   do jj=1,natom
     nn=typat(jj)
     rrad_nn=rrad(corlim(nn),nn)
     rrad2_nn=rrad_nn*rrad_nn
     vw1=vv(1)-xatm(1,jj)
     vw2=vv(2)-xatm(2,jj)
     vw3=vv(3)-xatm(3,jj)

     do kk=1,nnpos

       vt1=vw1-atp(1,kk)
       vt2=vw2-atp(2,kk)
       vt3=vw3-atp(3,kk)
       rr2=vt1*vt1+vt2*vt2+vt3*vt3

!      rr=vnorm(vt,0)

!      Only contribution > rhocormin (adhoc.f90) are considered

       if (rr2 < rrad2_nn .and.(.not.((selct==-2).and.(iat==jj).and.(ipos==kk)))) then
!        if (rr /= 0.0_dp) then    ! XG020629 : never test a real number against zero (not portable)
         if (rr2 > 1.0d-28) then         ! SEARCH INDEX

           rr=sqrt(rr2)
           rr_inv=1.0_dp/rr

           if (rr < rrad(1,nn)) then
             inx=-1
           elseif (rr > rrad(ndat(nn),nn)) then
             inx=ndat(nn)
           else
!            Find the index of the radius by bissection
             inmin=1
             inmax=ndat(nn)
             inx=1
             do
               if(inmax-inmin==1)exit
               inx=(inmin+inmax)/2
               if(rr>=rrad(inx,nn))then
                 inmin=inx
               else
                 inmax=inx
               end if
             end do
             inx=inmin

!            XG020629 : old coding, slower
!            inx=0
!            do while (rr >= rrad(inx+1,nn))
!            inx=inx+1
!            end do

           end if

!          Transformation matrix radial -> cart. coord
           ss=sqrt(vt1*vt1+vt2*vt2)
!          if (ss /=0._dp) then    ! XG020629 : never test a real number against zero (not portable)
           if (ss*ss > 1.0d-28) then  ! ss non-zero
!            XG020629 : very strange : only trsf(:,1) is needed in what follows ? !
!            ss_inv=1.0_dp/ss
             trsf(1,1)=vt1*rr_inv
!            trsf(1,2)=-vt2*ss_inv
!            trsf(1,3)=vt3*vt1*rr_inv*ss_inv
             trsf(2,1)=vt2*rr_inv
!            trsf(2,2)=vt1*ss_inv
!            trsf(2,3)=vt3*vt2*rr_inv*ss_inv
             trsf(3,1)=vt3*rr_inv
!            trsf(3,2)=0._dp
!            trsf(3,3)=-ss*rr_inv
!            XG020629 Not needed
!            do  ii=1,3
!            do ll=1,3
!            ches(ii,ll)=0._dp
!            end do
!            cgrad(ii)=0._dp
!            end do
           else                      ! ss zero
             do ii=1,3
               do ll=1,3
                 trsf(ii,ll)=0._dp
               end do
               trsf(ii,4-ii)=1._dp
             end do
           end if ! ss zero or non-zero

           if (inx == -1) then   ! LEFT EXTRAPOLATION y=a*x^2+b (a<0)
             val=sp2(1,nn)*0.5_dp*rr*rr/rrad(1,nn)+crho(1,nn)-sp2(1,nn)*rrad(1,nn)*0.5_dp
             cgrad(1)=sp2(1,nn)*rr/rrad(1,nn)
             ches(1,1)=sp2(1,nn)/rrad(1,nn)
           elseif (inx == ndat(nn) ) then  ! RIGHT EXTRAPOLATION y=a*exp(b*x)
             val=rrad(ndat(nn),nn)*exp(sp2(ndat(nn),nn)*(rr-rrad(ndat(nn),nn))/crho(ndat(nn),nn))
             cgrad(1)=val*sp2(ndat(nn),nn)/crho(ndat(nn),nn)
             ches(1,1)=cgrad(1)*sp2(ndat(nn),nn)/crho(ndat(nn),nn)
           else                    ! INTERPOLATION
             uu=rrad(inx+1,nn)-rrad(inx,nn)
             uu_inv=1.0_dp/uu
             aa=(rrad(inx+1,nn)-rr)*uu_inv
             bb=(rr-rrad(inx,nn))*uu_inv
             cc=(aa*aa*aa-aa)*uu*uu*0.16666666666666666_dp
             dd=(bb*bb*bb-bb)*uu*uu*0.16666666666666666_dp
             val=aa*crho(inx,nn)+bb*crho(inx+1,nn)+cc*sp3(inx,nn)+dd*sp3(inx+1,nn)
             cgrad(1)=(crho(inx+1,nn)-crho(inx,nn))*uu_inv&
&             -(3._dp*aa*aa-1._dp)*uu*0.16666666666666666_dp*sp3(inx,nn)+&
&             (3._dp*bb*bb-1._dp)*uu*0.16666666666666666_dp*sp3(inx+1,nn)
             ches(1,1)=aa*sp3(inx,nn)+bb*sp3(inx+1,nn)

           end if     ! TRANSFORMATION TO CARTEZ. COORD.

           cgrad1_rr_inv=cgrad(1)*rr_inv
           coeff=(ches(1,1)-cgrad1_rr_inv)*rr_inv*rr_inv
           cgrad(3)=trsf(3,1)*cgrad(1)
           cgrad(2)=trsf(2,1)*cgrad(1)
           cgrad(1)=trsf(1,1)*cgrad(1)
           ches(1,1)=coeff*vt1*vt1+cgrad1_rr_inv
           ches(2,2)=coeff*vt2*vt2+cgrad1_rr_inv
           ches(3,3)=coeff*vt3*vt3+cgrad1_rr_inv
           ches(1,2)=coeff*vt1*vt2 ; ches(2,1)=coeff*vt1*vt2
           ches(1,3)=coeff*vt1*vt3 ; ches(3,1)=coeff*vt1*vt3
           ches(2,3)=coeff*vt2*vt3 ; ches(3,2)=coeff*vt2*vt3

         else                                            ! case rr==0

           val=crho(1,nn)-sp2(1,nn)*rrad(1,nn)/2._dp
           do ii=1,3
             do ll=1,3
               ches(ii,ll)=0._dp
             end do
             cgrad(ii)=0._dp
             ches(ii,ii)=sp2(1,nn)/rrad(1,nn)
           end do

         end if ! rr>0 or rr==0

         do ii=1,3
           do ll=1,3
             hrho(ii,ll)=hrho(ii,ll)+ches(ii,ll)
           end do
           grho(ii)=grho(ii)+cgrad(ii)
         end do
         rho=rho+val

       end if ! rr2< rrad_nn*rrad_nn

       if (selct==-1) then
         if (rr2 < rminl(jj)*rminl(jj) ) then
           iat=jj
           ipos=kk
           rdmin=sqrt(rr2)
         end if
       elseif (selct==-2) then
         cycle
       else
         if (rr2 < rdmin*rdmin) then
           iat=jj
           ipos=kk
           rdmin=sqrt(rr2)
         end if
       end if

     end do
   end do

 end if

!write(51,'(":GRADt ",3F16.8)') grho
!write(52,'(":LAPt ",F16.8)') hrho(1,1)+hrho(2,2)+hrho(3,3)
!write(52,'(":HESNt ",/,3F16.8,/,3F16.8,/,3F16.8)') ((hrho(ii,jj),jj=1,3),ii=1,3)

!if(abs(cumul_cpu-cumul_cpu_old)>0.499)then
!write(std_out,'(a,f7.1)' )' vgh_rho : cumul_cpu=',cumul_cpu
!cumul_cpu_old=cumul_cpu
!end if

end subroutine vgh_rho
!!***

!!****f* m_bader/vnorm
!! NAME
!! vnorm
!!
!! FUNCTION
!! Default declarations, and interfaces for the aim.f utility.
!!
!! INPUTS
!! vector norm ->dir==1: vector in reduced coordinates
!!               dir==0: vector in cartes. coordinates
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! vv = on entry, vector to normalized
!!    = on exit, normalized vector
!!
!! SOURCE

function vnorm(vv,dir)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dir
 real(dp) :: vnorm
!arrays
 real(dp),intent(in) :: vv(3)

!Local variables-------------------------------
!scalars
 integer :: ii
!arrays
 real(dp) :: vt(3)

! *************************************************************************

 vnorm=zero
 if (dir==1) then
   do ii=1,3
     vt(ii)=rprimd(ii,1)*vv(1)+rprimd(ii,2)*vv(2)+rprimd(ii,3)*vv(3)
     vnorm=vnorm+vt(ii)*vt(ii)
   end do
 elseif (dir==0) then
   do ii=1,3
     vnorm=vnorm+vv(ii)*vv(ii)
   end do
 else
   MSG_ERROR('vnorm calcul')
 end if
 vnorm=sqrt(vnorm)
end function vnorm
!!***

!!****f* m_bader/vec_prod
!! NAME
!! vec_prod
!!
!! FUNCTION
!! Vector product
!!
!! INPUTS
!!  vv,uu = vectors to compute vector product
!!
!! OUTPUT
!!  (return the value of the vector product)
!!
!! SOURCE

function vec_prod(uu,vv)

!Arguments ------------------------------------
!arrays
 real(dp) :: vec_prod(3)
 real(dp),intent(in) :: uu(3),vv(3)

!Local variables-------------------------------

! *************************************************************************

 vec_prod(1)=uu(2)*vv(3)-vv(2)*uu(3)
 vec_prod(2)=uu(3)*vv(1)-vv(3)*uu(1)
 vec_prod(3)=uu(1)*vv(2)-vv(1)*uu(2)

end function vec_prod
!!***

!!****f* m_bader/mprod
!! NAME
!! mprod
!!
!! FUNCTION
!! Matrix multiplication cc=aa*bb
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine mprod(aa,bb,cc)

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: aa(3,3),bb(3,3)
 real(dp),intent(out) :: cc(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,kk

! *************************************************************************

 do ii=1,3
   do jj=1,3
     cc(ii,jj)=0._dp
     do kk=1,3
       cc(ii,jj)=cc(ii,jj)+aa(ii,kk)*bb(kk,jj)
     end do
   end do
 end do

end subroutine mprod
!!***

!!****f* m_bader/bschg1
!! NAME
!! bschg1
!!
!! FUNCTION
!! bschg1: Vector transformation of coordinates
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine bschg1(vv,dir)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dir
!arrays
 real(dp),intent(inout) :: vv(3)

!Local variables ------------------------------
!scalars
 integer :: ii
!arrays
 real(dp) :: vt(3)

! *********************************************************************

 if (dir==1) then
   do ii=1,3
     vt(ii)=rprimd(ii,1)*vv(1)+rprimd(ii,2)*vv(2)+rprimd(ii,3)*vv(3)
   end do
 elseif (dir==-1) then
   do ii=1,3
     vt(ii)=ivrprim(ii,1)*vv(1)+ivrprim(ii,2)*vv(2)+ivrprim(ii,3)*vv(3)
   end do
 elseif (dir==2) then
   do ii=1,3
     vt(ii)=trivrp(ii,1)*vv(1)+trivrp(ii,2)*vv(2)+trivrp(ii,3)*vv(3)
   end do
 else
   MSG_ERROR('Transformation of coordinates')
 end if
 vv(:)=vt(:)

end subroutine bschg1
!!***

!!****f* m_bader/bschg2
!! NAME
!! bschg2
!!
!! FUNCTION
!! bschg2: Matrix transformation of coordinates
!!
!! PARENTS
!!      m_bader
!!
!! CHILDREN
!!      mprod
!!
!! SOURCE

subroutine bschg2(aa,dir)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dir
!arrays
 real(dp),intent(inout) :: aa(3,3)

!Local variables ------------------------------
!arrays
 real(dp) :: bb(3,3)

! *********************************************************************

 if (dir==1) then
   call mprod(aa,ivrprim,bb)
   call mprod(rprimd,bb,aa)
 elseif (dir==2) then
   call mprod(aa,ivrprim,bb)
   call mprod(trivrp,bb,aa)
 elseif (dir==-1) then
   call mprod(aa,rprimd,bb)
   call mprod(ivrprim,bb,aa)
 else
   MSG_ERROR("transformation of coordinates")
 end if
end subroutine bschg2
!!***

end module m_bader
!!***
