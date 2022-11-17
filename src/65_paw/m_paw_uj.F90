!!****m* ABINIT/m_paw_uj
!! NAME
!!  m_paw_uj
!!
!! FUNCTION
!!  This module contains several routines relevant only for automatic determination of U
!!    in PAW+U context (linear response method according to Phys. Rev. B 71, 035105)
!!
!! COPYRIGHT
!! Copyright (C) 2018-2022 ABINIT group (DJA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_uj

 use defs_basis
 use m_abicore
 use m_errors
 use m_linalg_interfaces
 use m_xmpi
 use m_dtset
 use netcdf
 use m_nctk

 use m_fstrings,      only : strcat
 use m_pptools,       only : prmat
 use m_special_funcs, only : iradfnh
 use m_geometry,      only : shellstruct,ioniondist
 use m_parser,        only : prttagm
 use m_supercell,     only : mksupercell
 use m_pawrad,        only : pawrad_type
 use m_pawtab,        only : pawtab_type
 use m_paw_ij,        only : paw_ij_type
 use m_paral_atom,    only : get_my_atmtab, free_my_atmtab
 use m_dtfil,         only : datafiles_type
 use m_crystal,       only : crystal_t

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_paw_uj/macro_uj_type
!! NAME
!! dtmacro_uj
!!
!! FUNCTION
!! This data type contains the potential shifts and the occupations
!! for the determination of U and J for the DFT+U calculations.
!! iuj=1,3: non-selfconsistent calculations. iuj=2,4 selfconsistent calculations.
!!
!! SOURCE

 type, public :: macro_uj_type

! Integer
  integer :: iuj        ! dataset treated
  integer :: nat        ! number of atoms U (J) is determined on
  integer :: ndtset     ! total number of datasets
  integer :: nspden     ! number of densities treated
  integer :: macro_uj   ! which mode the determination runs in
  integer :: pawujat    ! which atom U (J) is determined on
  integer :: pawprtvol  ! controlling amount of output
  integer :: option     ! controls the determination of U (1 with compensating charge bath, 2 without)
  integer :: dmatpuopt  ! controls the renormalisation of the PAW projectors

! Real
  real(dp) :: diemix    ! mixing parameter
  real(dp) :: diemixmag ! magnetic mixing parameter
  real(dp) :: mdist     ! maximal distance of ions
  real(dp) :: pawujga   ! gamma for inversion of singular matrices
  real(dp) :: ph0phiint ! integral of phi(:,1)*phi(:,1)
  real(dp) :: pawujrad  ! radius to which U should be extrapolated.
  real(dp) :: pawrad    ! radius of the paw atomic data

! Integer arrays
  integer , allocatable  :: scdim(:)
  ! size of supercell

! Real arrays
  real(dp) , allocatable :: occ(:,:)
  ! occupancies after a potential shift: occ(ispden,nat)

  real(dp) , allocatable :: rprimd(:,:)
  ! unit cell for symmetrization

  real(dp) , allocatable :: vsh(:,:)
  ! potential shifts on atoms, dimensions: nspden,nat

  real(dp) , allocatable :: xred(:,:)
  ! atomic position for symmetrization

  real(dp) , allocatable :: wfchr(:)
  ! wfchr(1:3): zion, n and l of atom on which projection is done
  ! wfchr(4:6): coefficients ai of a0+a1*r+a2*r^2, fit to the wfc for r< r_paw

  real(dp), allocatable :: zioneff(:)
  ! zioneff(ij_proj), "Effective charge"*n "seen" at r_paw, deduced from Phi at r_paw, n:
  ! pricipal quantum number; good approximation to model wave function outside PAW-sphere through

 end type macro_uj_type
!!***

!public procedures.
 public :: pawuj_ini  ! Initialize dtpawuj datastructure
 public :: pawuj_free ! Deallocate dtpawuj datastructure
 public :: pawuj_det  ! Determine U (or J) parameter
 public :: pawuj_red  ! Store atomic occupancies, potential shift, positions

!private procedures.
!  chiscwrt:   distributes values of chi_org on chi_sc according to ion-ion distances
!  linvmat:    inverts real matrix inmat
!  lprtmat:    prints out the real matrix mmat
!  lcalcu:     prints out real the real matrice mmat
!  blow_pawuj: reads a real nxn matrice and appends lines n+1 and clumn n+1

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_uj/pawuj_ini
!! NAME
!! pawuj_ini
!!
!! FUNCTION
!!  Initialize dtpawuj datastructure for one SCF run
!!  Relevant only for automatic determination of U in PAW+U context
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dtpawuj(0:ndtpawuj) (initialization of fields vsh, occ, iuj,nnat)
!!
!! SOURCE

subroutine pawuj_ini(dtpawuj,ndtset)

!Arguments ------------------------------------
!scalars
 integer                           :: ndtset
 type(macro_uj_type),intent(inout) :: dtpawuj(0:ndtset)

!Local variables -------------------------
!Variables for partial dos calculation
!scalars
 integer, parameter :: nwfchr=6
 integer            :: iuj,im1
! *********************************************************************

 DBG_ENTER("COLL")
 do iuj=0,ndtset
   !write(std_out,*)'pawuj_ini iuj ',iuj
   dtpawuj(iuj)%diemix=0.45_dp
   dtpawuj(iuj)%diemixmag=0.45_dp
   dtpawuj(iuj)%iuj=0
   dtpawuj(iuj)%nat=0
   dtpawuj(iuj)%ndtset=1
   dtpawuj(iuj)%nspden=1
   dtpawuj(iuj)%macro_uj=0
   dtpawuj(iuj)%option=1
   dtpawuj(iuj)%pawujat=1
   dtpawuj(iuj)%pawujga=one
   dtpawuj(iuj)%pawprtvol=1
   dtpawuj(iuj)%ph0phiint=one
   dtpawuj(iuj)%dmatpuopt=2
   dtpawuj(iuj)%pawujrad=3.0_dp
   dtpawuj(iuj)%pawrad=20.0_dp
   !Allocate arrays
   !write(std_out,*)'pawuj_ini before arrays'
   ABI_MALLOC(dtpawuj(iuj)%rprimd,(3,3))
   ABI_MALLOC(dtpawuj(iuj)%scdim,(3))
   ABI_MALLOC(dtpawuj(iuj)%wfchr,(nwfchr))
   dtpawuj(iuj)%rprimd=reshape((/ 1,0,0,0,1,0,0,0,1/),(/ 3,3 /))
   dtpawuj(iuj)%scdim=reshape((/ 250,0,0 /),(/3 /))
   dtpawuj(iuj)%wfchr=(/ (0,im1=1,nwfchr) /)
   if (iuj>0) then
     dtpawuj(iuj)%iuj=-1
   end if
 end do

 DBG_EXIT("COLL")

end subroutine pawuj_ini
!!***

!----------------------------------------------------------------------

!!****f* m_paw_uj/pawuj_free
!! NAME
!!  pawuj_free
!!
!! FUNCTION
!!   deallocate pawuj stuff
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine pawuj_free(dtpawuj)

!Arguments -------------------------------
 type(macro_uj_type),intent(inout) :: dtpawuj

! *********************************************************************

 if (allocated(dtpawuj%scdim))    then
   ABI_FREE(dtpawuj%scdim)
 end if
 if (allocated(dtpawuj%occ))      then
   ABI_FREE(dtpawuj%occ)
 end if
 if (allocated(dtpawuj%rprimd))   then
   ABI_FREE(dtpawuj%rprimd)
 end if
 if (allocated(dtpawuj%vsh))      then
   ABI_FREE(dtpawuj%vsh)
 end if
 if (allocated(dtpawuj%xred))     then
   ABI_FREE(dtpawuj%xred)
 end if
 if (allocated(dtpawuj%wfchr))    then
   ABI_FREE(dtpawuj%wfchr)
 end if
 if (allocated(dtpawuj%zioneff))  then
   ABI_FREE(dtpawuj%zioneff)
 end if

end subroutine pawuj_free
!!***

!----------------------------------------------------------------------
!!****f* m_paw_uj/pawuj_det
!! NAME
!!  pawuj_det
!!
!! FUNCTION
!!  From the complete dtpawuj-dataset determines U (or J) parameter for
!!  PAW+U calculations
!!  Relevant only for automatic determination of Hubbard Parameters in
!!  PAW+U context.
!!  Hubbard U = diemix/chi0-1/chi
!!  Hund's J = 1/chi-diemixmag/chi0
!!
!! INPUTS
!!  dtpawuj=potential shifts (vsh) and atomic occupations (occ)
!!
!! OUTPUT
!!  only printing
!!  (among other things a section in the ab.out that can be used for input in ujdet)
!!
!! SOURCE

subroutine pawuj_det(dtpawuj,ndtpawuj,dtset,dtfil,ures,comm)

!Arguments ------------------------------------
!scalars
!arrays
 integer                        :: ndtpawuj, comm
 type(macro_uj_type),intent(in) :: dtpawuj(0:ndtpawuj)
 type(dataset_type),intent(in)      :: dtset
 type(datafiles_type),intent(in) :: dtfil
 real(dp),intent(out)           :: ures

!Local variables-------------------------------
!scalars
 integer,parameter           :: natmax=2,nwfchr=6
 integer                     :: ii,jj,nat_org,jdtset,nspden,macro_uj,kdtset,marr
 integer                     :: im1,ndtuj,idtset, nsh_org, nsh_sc,nat_sc,maxnat
 integer                     :: pawujat,pawprtvol,pawujoption
 integer                     :: dmatpuopt,invopt,ipert
 integer                     :: my_rank, ncid, ncerr
 real(dp)                    :: pawujga,ph0phiint,intg,fcorr,eyp
 real(dp)                    :: diem,signum,scalarHP !LMac quantities

 character(len=500)          :: message
 character(len=2)            :: hstr
 character(len=5)            :: pertname
 character(len=1)            :: parname
 character(len=14)           :: occmag
 type(crystal_t) :: cryst
!arrays
 integer                     :: ext(3)
 real(dp)                    :: rprimd_sc(3,3),vsh(ndtpawuj),a(5),b(5)
 integer,allocatable         :: narrm(:)
 integer,allocatable         :: idum2(:,:),jdtset_(:),smult_org(:),smult_sc(:)
 real(dp),allocatable        :: luocc(:,:),dqarr(:,:),dqarrr(:,:),dparr(:,:),dparrr(:,:),xred_org(:,:),drarr(:,:)
 real(dp),allocatable        :: magv_org(:),magv_sc(:),chi(:),chi0(:),chi0_sc(:), chi_sc(:), xred_sc(:,:)
 real(dp),allocatable        :: sdistv_org(:),sdistv_sc(:),distv_org(:),distv_sc(:)
 integer,parameter           :: ncid0 = 0, master = 0
! *********************************************************************

 DBG_ENTER("COLL")

 my_rank = xmpi_comm_rank(comm)

!write(std_out,*) 'pawuj 01'
!###########################################################
!### 01. Allocations

!Initializations
 ndtuj=count(dtpawuj(:)%iuj/=-1)-1 ! number of datasets initialized by pawuj_red
 ABI_MALLOC(jdtset_,(0:ndtuj))
 jdtset_(0:ndtuj)=pack(dtpawuj(:)%iuj,dtpawuj(:)%iuj/=-1)
 jdtset=maxval(dtpawuj(:)%iuj)

!DEBUG
write(message,'(10(a,i3))')'pawuj_det jdtset ',jdtset,&
& ' ndtuj ', ndtuj,' ndtpawuj ',ndtpawuj
call wrtout(std_out,message,'COLL')
!call flush_unit(6)
!END DEBUG

 nspden=dtpawuj(jdtset)%nspden
 nat_org=dtpawuj(jdtset)%nat
 macro_uj=dtpawuj(jdtset)%macro_uj
 pawujat=dtpawuj(jdtset)%pawujat
 pawprtvol=dtpawuj(jdtset)%pawprtvol
 pawujga=dtpawuj(jdtset)%pawujga
 pawujoption=dtpawuj(jdtset)%option
 ph0phiint=dtpawuj(jdtset)%ph0phiint
 dmatpuopt=dtpawuj(jdtset)%dmatpuopt
 marr=maxval((/ 9, nspden*nat_org ,nat_org*3 /))
 eyp=2.5_dp ! for dmatpuopt==2 and 3
 if (dmatpuopt==1) eyp=eyp+3.0_dp
 if (dmatpuopt>=3) eyp=(eyp+3.0_dp-dmatpuopt)

 ABI_MALLOC(luocc,(ndtpawuj,nat_org))
 ABI_MALLOC(idum2,(marr,0:ndtuj))
 ABI_MALLOC(drarr,(marr,0:ndtuj))
 ABI_MALLOC(magv_org,(nat_org))
 ABI_MALLOC(xred_org,(3,nat_org))
 ABI_MALLOC(chi0,(nat_org))
 ABI_MALLOC(chi,(nat_org))
 ABI_MALLOC(dparr,(marr,0:ndtuj))
 ABI_MALLOC(dparrr,(marr,0:ndtuj))
 ABI_MALLOC(dqarr,(marr,0:ndtuj))
 ABI_MALLOC(dqarrr,(marr,0:ndtuj))
 ABI_MALLOC(distv_org,(nat_org))
 ABI_MALLOC(narrm,(0:ndtuj))
 dparr=-one ;  dparrr=-one ;  dqarr=-one ;  dqarrr=-one
!DEBUG
!write(message,fmt='((a,i3,a))')'pawuj_det init sg'
!call wrtout(std_out,message,'COLL')
!END DEBUG
 idum2=1 ; drarr=one
 luocc=zero

!write(std_out,*) 'pawuj 03'
!###########################################################
!### 03. Write out the input for the ujdet utility

 write(message, '(3a)' ) ch10,&
& ' # input for ujdet, cut it using ''sed -n "/MARK/,/MARK/p" abi.out > ujdet.in ''------- ',ch10
 call wrtout(ab_out,message,'COLL')

 if (ndtuj/=4) then
   write (hstr,'(I0)') jdtset              ! convert integer to string
 else
   hstr=''
 end if
 if (ndtuj>=2.or.jdtset==1) then
   idum2(1,1:ndtuj)=4 !dtpawuj(:)%ndtuj
   call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid0,ndtuj,'ndtset','INT',0)
 end if
 idum2(1,0:ndtuj)=pack(dtpawuj(:)%nat,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid0,ndtuj,'nat'//trim(hstr),'INT',0)

 idum2(1,0:ndtuj)=pack(dtpawuj(:)%nspden,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid0,ndtuj,'nspden'//trim(hstr),'INT',0)

 idum2(1,0:ndtuj)=pack(dtpawuj(:)%macro_uj,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid0,ndtuj,'macro_uj'//trim(hstr),'INT',0)

 idum2(1,0:ndtuj)=pack(dtpawuj(:)%pawujat,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid0,ndtuj,'pawujat'//trim(hstr),'INT',0)

 dparr(1,0:ndtuj)=pack(dtpawuj(:)%pawujga,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid0,ndtuj,'pawujga'//trim(hstr),'DPR',0)

 dparr(1,0:ndtuj)=pack(dtpawuj(:)%pawujrad,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid0,ndtuj,'pawujrad'//trim(hstr),'DPR',0)

 dparr(1,0:ndtuj)=pack(dtpawuj(:)%pawrad,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid0,ndtuj,'pawrad'//trim(hstr),'DPR',0)

 dparr(1,0:ndtuj)=pack(dtpawuj(:)%ph0phiint,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid0,ndtuj,'ph0phiint'//trim(hstr),'DPR',0)

 idum2(1,0:ndtuj)=pack(dtpawuj(:)%pawprtvol,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid0,ndtuj,'pawprtvol'//trim(hstr),'INT',0)

 idum2(1,0:ndtuj)=pack(dtpawuj(:)%option,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid0,ndtuj,'pawujopt'//trim(hstr),'INT',0)

 idum2(1,0:ndtuj)=pack(dtpawuj(:)%dmatpuopt,dtpawuj(:)%iuj/=-1)
 call prttagm(dparr,idum2,ab_out,jdtset_,1,marr,1,narrm,ncid0,ndtuj,'dmatpuopt'//trim(hstr),'INT',0)

 kdtset=0

 do idtset=0,ndtpawuj
   if (dtpawuj(idtset)%iuj/=-1) then
     dparr(1:nspden*nat_org,kdtset)=reshape(dtpawuj(idtset)%vsh,(/nspden*nat_org/))
     dparrr(1:nspden*nat_org,kdtset)=reshape(dtpawuj(idtset)%occ,(/nspden*nat_org/))
     dqarr(1:nat_org*3,kdtset)=reshape(dtpawuj(idtset)%xred,(/nat_org*3/))
     dqarrr(1:3*3,kdtset)=reshape(dtpawuj(idtset)%rprimd,(/3*3/))
     idum2(1:3,kdtset)=reshape(dtpawuj(idtset)%scdim,(/3/))
     drarr(1:nwfchr,kdtset)=reshape(dtpawuj(idtset)%wfchr,(/nwfchr/))
     kdtset=kdtset+1
   end if
 end do

 call prttagm(dparr,idum2,ab_out,jdtset_,2,marr,nspden*nat_org,narrm,ncid0,ndtuj,'vsh'//trim(hstr),'DPR',0)
 call prttagm(dparrr,idum2,ab_out,jdtset_,2,marr,nspden*nat_org,narrm,ncid0,ndtuj,'occ'//trim(hstr),'DPR',0)
 call prttagm(dqarr,idum2,ab_out,jdtset_,2,marr,nat_org*3,narrm,ncid0,ndtuj,'xred'//trim(hstr),'DPR',0)
 call prttagm(dqarrr,idum2,ab_out,jdtset_,2,marr,3*3,narrm,ncid0,ndtuj,'rprimd'//trim(hstr),'DPR',0)
 call prttagm(dqarrr,idum2,ab_out,jdtset_,2,marr,3,narrm,ncid0,ndtuj,'scdim'//trim(hstr),'INT',0)
 call prttagm(drarr,idum2,ab_out,jdtset_,2,marr,nwfchr,narrm,ncid0,ndtuj,'wfchr'//trim(hstr),'DPR',0)
 ABI_FREE(narrm)

 write(message, '( 15a )'  ) ch10,' # further possible options: ',ch10,&
& ' #    scdim    2 2 2 ',ch10,&
& ' #    mdist    10.0  ',ch10,&
& ' #  pawujga    2 '    ,ch10,&
& ' # pawujopt    2 '    ,ch10,&
& ' # pawujrad    3.0'   ,ch10,&
& ' # ------- end input for ujdet: end-MARK  -------- ',ch10
 call wrtout(ab_out,message,'COLL')

!###########################################################
!### 04. Testing consistency of parameters and outputting info

!LMac
! if (ndtuj/=4)  return

 write(message, '(3a)' ) ch10,' ---------- calculate U, (J) start ---------- ',ch10
 call wrtout(ab_out,message,'COLL')

!Tests if perturbed atom is consistent with ujdet atom
 if (all(dtpawuj(1:ndtpawuj)%pawujat==pawujat)) then
   write (message,fmt='(a,i3)') ' All pawujat  ok and equal to ',pawujat
   call wrtout(ab_out,message,'COLL')
 else
   write (message,fmt='(a,4i3,2a)') ' Differing values of pawujat were found: ',dtpawuj(1:ndtuj)%pawujat,ch10,&
&   'No determination of U.'
   call wrtout(ab_out,message,'COLL')
   return
 end if

!Tests consistency of macro_uj, then writes message about macro_uj procedure selected. LMac
 if (nspden==1) then
   write(message,fmt='(2a)') ' nspden==1, determination',&
&   ' of U-parameter for unpolarized structure (non standard)'
 else if (macro_uj==1.and.nspden==2) then
   write(message,fmt='(2a)') ' macro_uj=1 and nspden=2:',&
&   ' standard determination of Hubbard U-parameter'
 else if (macro_uj==2.and.nspden==2) then
   write(message,fmt='(2a)') ' macro_uj=2 and nspden=2:',&
&   ' determination of parameter on single spin channel (experimental)'
 else if (macro_uj==3.and.nspden==2) then
   write(message,fmt='(2a)') ' macro_uj=3 and nspden=2,',&
&   ' determination of (not Hunds) J-parameter on single spin channel (experimental)'
 else if (macro_uj==4.and.nspden==2) then
   write(message,fmt='(2a)') ' macro_uj=4 and nspden=2,',&
&   ' Hunds J determination – L. MacEnulty August 2021'
 end if
 call wrtout(ab_out,message,'COLL')

!Tests compatibility of nspden and macro_uj
 if (macro_uj>1.and.nspden==1) then
   write (message,'(4a,2a)') ' U on a single spin channel (or J) can only be determined for nspden=2 ,',ch10,&
&   'No determination of U.'
   call wrtout(ab_out,message,'COLL')
   return
 end if

!Calculation of response matrix
 diem=dtpawuj(1)%diemix !Unscreened response in Hubbard U impacted by diemix
 pertname='alpha' !Hubbard U perturbation; applied equally to spins up and down
 parname='U'
 occmag='  Occupations'
 signum=1.0d0 !Hubbard U is signum*(1/chi0-1/chi)
 do jdtset=1,4
   if (nspden==1) then
     luocc(jdtset,1:nat_org)=dtpawuj(jdtset)%occ(1,:)
   !Hubbard U: uses sum of spin up and spin down occupancies
   else if (macro_uj==1.and.nspden==2) then
     luocc(jdtset,1:nat_org)=dtpawuj(jdtset)%occ(1,:)+dtpawuj(jdtset)%occ(2,:) !Total occupation
   else if (macro_uj==2.and.nspden==2) then
     luocc(jdtset,1:nat_org)=dtpawuj(jdtset)%occ(1,:)
   else if (macro_uj==3.and.nspden==2) then
     luocc(jdtset,1:nat_org)=dtpawuj(jdtset)%occ(2,:)
     parname='J'
   !Hund's J: uses difference of spin up and spin down occupancies
   else if (macro_uj==4.and.nspden==2) then
     luocc(jdtset,1:nat_org)=dtpawuj(jdtset)%occ(1,:)-dtpawuj(jdtset)%occ(2,:) !Magnetization
     diem=dtpawuj(1)%diemixmag !Unscreened response in Hund's J impacted by diemixmag
     pertname='beta ' !Hund's J perturbation: +beta to spin up, -beta to down
     parname='J'
     occmag='Magnetizations'
     signum=-1.0d0 !Hund's J is -1*(1/chi0-1/chi)
   end if
   vsh(jdtset)=dtpawuj(jdtset)%vsh(1,pawujat)
!   if (pawprtvol==-3) then
     write(message,fmt='(2a,i3,a,f15.12)') ch10,' Potential shift vsh(',jdtset,') =',vsh(jdtset)
     call wrtout(std_out,message,'COLL')
     write(message,fmt='( a,i3,a,120f15.9)') ' Occupations occ(',jdtset,') ',luocc(jdtset,1:nat_org)
     call wrtout(std_out,message,'COLL')
!   end if
 end do

 write(message,fmt='( a)') 'Occupations assigned.'
 call wrtout(std_out,message,'COLL')

 !Two-point linear regression of response matrices.
 chi0=(luocc(1,1:nat_org)-luocc(3,1:nat_org))/(vsh(1)-vsh(3))/diem
 chi=(luocc(2,1:nat_org)-luocc(4,1:nat_org))/(vsh(2)-vsh(4))

 if ((abs(chi0(pawujat))<0.0000001).or.(abs(chi(pawujat))<0.0000001)) then
   write(message, '(2a,2f12.5,a)' ) ch10,'Chi0 or Chi is too small for inversion.',&
     &chi0(pawujat),chi(pawujat),ch10
   call wrtout(ab_out,message,'COLL')
   return
 end if

 !LMac: Scalar Hubbard Parameter
 scalarHP=signum*(1.0d0/chi0(pawujat)-1.0d0/chi(pawujat))*Ha_eV

 write(message,fmt='(a)')': '
 if (nspden==2) then
   magv_org=dtpawuj(1)%occ(1,:)-dtpawuj(1)%occ(2,:)
   if (all(abs(magv_org)<0.001)) then
     magv_org=(/(1,im1=1,nat_org)/)
   else
     magv_org=abs(magv_org)/magv_org
     if (magv_org(1).lt.0) magv_org=magv_org*(-1_dp)
     if (all(magv_org(2:nat_org).lt.0)) then
       magv_org=abs(magv_org)
       write(message,'(a)')', (reset to fm): '
     end if
   end if
 else
   magv_org=(/(1,im1=1,nat_org)/)
 end if

 if (pawprtvol==-3) then
   write(message,fmt='(3a, 150f4.1)') ch10,' Magnetization',trim(message),magv_org
   call wrtout(std_out,message,'COLL')
 end if

!Case of extrapolation to larger r_paw: calculate intg
 if (all(dtpawuj(1)%wfchr(:)/=0).and.ph0phiint/=1) then
   if (dtpawuj(1)%pawujrad<20.and.dtpawuj(1)%pawujrad>dtpawuj(1)%pawrad) then
     fcorr=(1-ph0phiint)/(IRadFnH(dtpawuj(1)%pawrad,20.0_dp,nint(dtpawuj(1)%wfchr(2)),&
&     nint(dtpawuj(1)%wfchr(3)),dtpawuj(1)%wfchr(1)))
     intg=ph0phiint/(1-fcorr*IRadFnH(dtpawuj(1)%pawujrad,20.0_dp,nint(dtpawuj(1)%wfchr(2)),&
&     nint(dtpawuj(1)%wfchr(3)),dtpawuj(1)%wfchr(1)))
     write(message, fmt='(a,f12.5,a,f12.5)') ' pawuj_det: met2 extrapolation to ', dtpawuj(1)%pawujrad,' using factor ',intg
     call wrtout(std_out,message,'COLL')
   else if (dtpawuj(1)%pawujrad<dtpawuj(1)%pawrad) then
     a=0 ; a(1:3)=dtpawuj(1)%wfchr(4:6)
     a=(/a(1)**2,a(1)*a(2),a(2)**2/3.0_dp+2.0_dp/3.0_dp*a(1)*a(3),a(2)*a(3)/2.0_dp,a(3)**2/5.0_dp/)
     b=(/(dtpawuj(1)%pawujrad**(im1)-dtpawuj(1)%pawrad**(im1),im1=1,5)/)
     intg=dot_product(a,b)
     intg=ph0phiint/(ph0phiint+intg)
     write(message, fmt='(a,f12.5,a,f12.5)') ' pawuj_det: met1 extrapolation to ', dtpawuj(1)%pawujrad,' using factor ',intg
     call wrtout(std_out,message,'COLL')
   else if (dtpawuj(1)%pawujrad==dtpawuj(1)%pawrad) then
     intg=one
     write(message, fmt='(a,2i7,3f12.5)') ' pawuj_det: no extrapolation (pawujrad=pawrad)'
     call wrtout(std_out,message,'COLL')
   else
     intg=ph0phiint
     write(message, fmt='(a,3f12.5)') ' pawuj_det: extrapolation to r->\inf using factor ',intg
     call wrtout(std_out,message,'COLL')
   end if
 else
   write(message, fmt='(a,2i7,3f12.5)') ' pawuj_det: no extrapolation (ph0phiint=1 or wfchr=0)'
   call wrtout(std_out,message,'COLL')
   intg=one
 end if

!Determine U in primitive cell
 write(message,fmt='(a)')' pawuj_det: determine U in primitive cell'
 call wrtout(std_out,message,'COLL')

 call lcalcu(int(magv_org),nat_org,dtpawuj(1)%rprimd,dtpawuj(1)%xred,chi,chi0,pawujat,ures,pawprtvol,pawujga,pawujoption)

!Begin calculate U in supercell

!Analyze shell structure of primitive cell
!and atomic distances in primitive cell
 ABI_MALLOC(smult_org,(nat_org))
 ABI_MALLOC(sdistv_org,(nat_org))
 call shellstruct(dtpawuj(1)%xred,dtpawuj(1)%rprimd,nat_org,&
& int(magv_org),distv_org,smult_org,sdistv_org,nsh_org,pawujat,pawprtvol)


!LMac: Printing relevant information about the Hubbard parameter just calculated.
 write(message,'(3a)') ch10,ch10,'*********************************************************************'
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 write(message,'(4a)') '************************  Linear Response ',parname,'  ************************',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 write(message, '(a,i4,a)' ) ' Info printed for perturbed atom: ',pawujat,ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 write(message, fmt='(10a)')'  Perturbations         ',occmag,ch10,&
' --------------- -----------------------------',ch10,&
'    ',pertname,' [eV]     Unscreened      Screened',ch10,&
' --------------- -----------------------------'
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 do ipert=1,2
   write(message, fmt='(3f15.10)') vsh(ipert*2-1)*Ha_eV,luocc(ipert*2-1,pawujat),luocc(ipert*2,pawujat)
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 end do
 write(message,'(2a)') ch10,'                    Scalar response functions:'
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 write(message,fmt='(a,f12.5)') '                    Chi0 [eV^-1]: ',chi0(pawujat)/Ha_eV
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 write(message,fmt='(a,f12.5)') '                    Chi [eV^-1]:  ',chi(pawujat)/Ha_eV
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 write(message,'(4a,f9.5,a)') ch10,' The scalar ',parname,' from the two-point regression scheme is ',scalarHP,' eV.'
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 write(message,'(3a)') '*********************************************************************',ch10,&
'*********************************************************************'
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 write(message,'(7a)') 'Note: For more reliable linear regressions of the response',ch10,&
'matrices, it is advised that you have more than two points.',ch10,&
'See the LRUJ protocol for more information.',ch10,ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

!###########################################################
!### Create the file LRUJ.nc for the LRUJ ujdet utility

 if (my_rank == master) then
   ! First call:
   !  - Create NC file, define dimensions, scalars and arrays.
   !  - Add crystal structure and metadata required to post-process the data.
   NCF_CHECK(nctk_open_create(ncid, strcat(dtfil%filnam_ds(4), "_LRUJ.nc"), xmpi_comm_self))

   cryst = dtset%get_crystal(1)
   NCF_CHECK(cryst%ncwrite(ncid))
   call cryst%free()

   ! Define dimensions.
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("nnat", nat_org), &
     nctkdim_t("ndtpawuj", ndtpawuj), &
     nctkdim_t("nspden", dtset%nspden), &
     nctkdim_t("nsppol", dtset%nsppol) ], &
     defmode=.True.)
   NCF_CHECK(ncerr)

   ! Define integer scalars
   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
     "usepaw", "macro_uj", "pawujat", "nspden", "dmatpuopt"  &
   ])

   ! Define double precision scalars
   ! @lmacenul: Here I write pawujv so that one can order the results by alpha in the post-processing tool.
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
     "diemix", "diemixmag", "ph0phiint", "uj_pert" &
   ])
   NCF_CHECK(ncerr)

   ! Define arrays with results.
   ! TODO: Here I need an extra dimension to store results with different iuj and/or different names.
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("luocc", "dp", "ndtpawuj, nnat") &
     !nctkarr_t("vsh", "dp", "nspden, nnat") &
   ])
   NCF_CHECK(ncerr)

   ! ===========================================
   ! Write metadata that does not depend on icyc
   ! ===========================================
   NCF_CHECK(nctk_set_datamode(ncid))

   ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
     "usepaw", "macro_uj", "pawujat", "nspden", "dmatpuopt"],  &
     [dtset%usepaw, macro_uj, pawujat, nspden, dmatpuopt])
   NCF_CHECK(ncerr)

   associate (dt => dtpawuj(1))
   ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
     "diemix", "diemixmag", "ph0phiint", "uj_pert" ], &
     [dt%diemix, dt%diemixmag, ph0phiint, vsh(3)])
   NCF_CHECK(ncerr)
   end associate

   ! Write arrays
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "luocc"), luocc))

   NCF_CHECK(nf90_close(ncid))
 end if

 ii=1
 write(message, fmt='(8a)') ' URES ','     ii','    nat','       r_max','    U(J)[eV]','   U_ASA[eV]','   U_inf[eV]',ch10
 write(message, fmt='(a,2i7,4f12.5)') trim(message)//' URES ',ii,nat_org,maxval(abs(distv_org)),signum*ures,signum*ures*exp(log(intg)*eyp),&
& signum*ures*exp(log(ph0phiint)*eyp)
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 if (pawprtvol>1) then
   write(message,fmt='(a, 150f10.5)')' pawuj_det: ionic distances in original cell (distv_org) ', distv_org
   call wrtout(std_out,message,'COLL')
 end if

!Construct supercell, calculate limit dimensions of supercell
 ii=1
 maxnat=product(dtpawuj(1)%scdim(:))*nat_org
 if (maxnat==0) then
   maxnat=dtpawuj(1)%scdim(1)
   ext=(/ii, ii, ii/)
 else
   jj=1
   do while (jj<=3)
     ext(jj)=minval( (/ii, dtpawuj(1)%scdim(jj) /) )
     jj=jj+1
   end do
 end if
 ext=ext+(/ 1,1,1 /)
 ii=ii+1


 nat_sc=product(ext)*nat_org

!DEBUG
!write(message,fmt='(a,3i4)')'pawuj_det: ext ',ext
!call wrtout(std_out,message,'COLL')
!END DEBUG

 do while (nat_sc<=maxnat)
   ABI_MALLOC(chi0_sc,(nat_sc))
   ABI_MALLOC(chi_sc,(nat_sc))
   ABI_MALLOC(distv_sc,(nat_sc))
   ABI_MALLOC(magv_sc,(nat_sc))
   ABI_MALLOC(sdistv_sc,(nat_sc))
   ABI_MALLOC(smult_sc,(nat_sc))
   ABI_MALLOC(xred_sc,(3,nat_sc))

!  Determine positions=xred_sc and supercell dimensions=rpimd_sc
   call mksupercell(dtpawuj(1)%xred,int(magv_org),dtpawuj(1)%rprimd,nat_org,&
&   nat_sc,xred_sc,magv_sc,rprimd_sc,ext,pawprtvol)

!  Determine shell structure of supercell: multiplicities (smult_sc), radii of shells (sdistv_sc)
!  number of shells (nsh_sc) and atomic distances in supercell (distv_sc)

   write(message,fmt='(a,3i3,a)')' pawuj_det: determine shell structure of ',ext(1:3),' supercell'
   call wrtout(std_out,message,'COLL')

   call shellstruct(xred_sc,rprimd_sc,nat_sc,int(magv_sc),distv_sc,smult_sc,sdistv_sc,nsh_sc,pawujat,pawprtvol)

   if (pawprtvol>=2) then
     write(message,fmt='(a)')' pawuj_det: ionic distances in supercell (distv_sc) '
     call wrtout(std_out,message,'COLL')
     call prmat(distv_sc,1,nat_sc,1,std_out)
   end if

!  Determine chi and chi0 in supercell (chi0_sc, chi_sc)
!  DEBUG
!  write(message,fmt='(a)')'pawuj_det:  chi and chi0 in supercell'
!  call wrtout(std_out,message,'COLL')
!  END DEBUG

   if (pawujoption>2) then
     invopt=2
   else
     invopt=1
   end if

   call chiscwrt(chi0,distv_org,nat_org,sdistv_org,smult_org,nsh_org,&
&   chi0_sc,distv_sc,nat_sc,smult_sc,nsh_sc,invopt,pawprtvol)
   call chiscwrt(chi,distv_org,nat_org,sdistv_org,smult_org,nsh_org,&
&   chi_sc,distv_sc,nat_sc,smult_sc,nsh_sc,invopt,pawprtvol)

!  Calculate U in supercell
!  DEBUG
!  write(message,fmt='(a)')'pawuj_det:   U in supercell'
!  call wrtout(std_out,message,'COLL')
!  END DEBUG
   call lcalcu(int(magv_sc),nat_sc,rprimd_sc,xred_sc,chi_sc,chi0_sc,pawujat,ures,pawprtvol,pawujga,pawujoption)

   write(message, fmt='(a,2i7,4f12.5)') ' URES ',ii,nat_sc,maxval(abs(distv_sc)),signum*ures,signum*ures*exp(log(intg)*eyp),&
&   signum*ures*exp(log(ph0phiint)*eyp)
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   ABI_FREE(chi0_sc)
   ABI_FREE(chi_sc)
   ABI_FREE(distv_sc)
   ABI_FREE(magv_sc)
   ABI_FREE(sdistv_sc)
   ABI_FREE(smult_sc)
   ABI_FREE(xred_sc)

   if (product(dtpawuj(1)%scdim(:))==0) then
     ext=(/ii, ii, ii/)
   else
     jj=1
     do while (jj<=3)
       ext(jj)=minval( (/ii, dtpawuj(1)%scdim(jj) /) )
       jj=jj+1
     end do
   end if
   ext=ext+(/ 1,1,1 /)
   ii=ii+1

   nat_sc=product(ext)*nat_org

 end do

 write(message, '(3a)' )ch10,' ---------- calculate U, (J) end -------------- ',ch10
 call wrtout(ab_out,message,'COLL')

 ABI_FREE(dparrr)
 ABI_FREE(dqarr)
 ABI_FREE(dqarrr)
 ABI_FREE(jdtset_)
 ABI_FREE(chi)
 ABI_FREE(chi0)
 ABI_FREE(smult_org)
 ABI_FREE(sdistv_org)
 ABI_FREE(luocc)
 ABI_FREE(idum2)
 ABI_FREE(drarr)
 ABI_FREE(dparr)
 ABI_FREE(magv_org)
 ABI_FREE(xred_org)
 ABI_FREE(distv_org)

 DBG_EXIT("COLL")

end subroutine pawuj_det
!!***

!----------------------------------------------------------------------
!!****f* m_paw_uj/pawuj_red
!! NAME
!!  pawuj_red
!!
!! FUNCTION
!!  Store atomic occupancies, potential shift, positions in dtpawuj datastructure.
!!
!! INPUTS
!!  istep: SCF iteration step
!!  pert_state: 0 if the routine is called with unperturbed occupancies
!!              1 if the routine is called after having applied the perturbation.
!!  fatvshift=factor that multiplies atvshift
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell
!!  ntypat = number of atom types
!!  paw_ij(my_natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawprtvol= printing volume
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  comm=MPI commuicator
!!
!! OUTPUT
!!  dtpawuj(0:ndtpawuj) (initialization of fields vsh, occ, occ0, iuj,nnat)
!!
!! SOURCE

subroutine pawuj_red(istep, pert_state, dtfil, &
                     dtset,dtpawuj,fatvshift,my_natom,natom,ntypat,paw_ij,pawrad,pawtab,ndtpawuj, comm, &
&                    mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ------------------------------------
!scalars
 integer,intent(in)                 :: istep, pert_state, my_natom,natom,ntypat,ndtpawuj, comm
 integer,optional,intent(in)        :: comm_atom
 real(dp),intent(in)                :: fatvshift
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(paw_ij_type),intent(in)       :: paw_ij(my_natom)
 type(pawtab_type),intent(in)       :: pawtab(ntypat)
 type(pawrad_type),intent(in)       :: pawrad(ntypat)
 type(dataset_type),intent(in)      :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(macro_uj_type),intent(inout)  :: dtpawuj(0:ndtpawuj)

!Local variables-------------------------------
!scalars
 integer,parameter           :: natmax=2,ncoeff=3, master = 0
 integer                     :: iatom,iatom_tot,ierr,im1,im2,ispden,itypat,ll,nspden,nsppol,iuj,ncyc,icyc
 integer                     :: my_comm_atom,nnat,natpawu,natvshift,pawujat,ndtset,typawujat
 !integer                     :: my_rank, ncid, ncerr
 logical                     :: usepawu !antiferro,
 logical                     :: my_atmtab_allocated,paral_atom
 character(len=1000)         :: message,hstr
 character(len=500)          :: messg
!arrays
 logical                     :: musk(3,natom)
 integer,pointer             :: my_atmtab(:)
 real(dp)                    :: rrtab(ncoeff),wftab(ncoeff),a(ncoeff,ncoeff),b(ncoeff,ncoeff)! ,a(ncoeff,ncoeff)
 real(dp),allocatable        :: nnocctot(:,:) !,magv(:)
 real(dp),allocatable        :: atvshift(:,:,:) ! atvshift(natvshift,2,natom)
 logical,allocatable         :: dmusk(:,:),atvshmusk(:,:,:) !atvshmusk(natvshift,2,natom)

! *********************************************************************

!Initializations
 nspden=1;nsppol=1
 if (my_natom>0) then
   nspden=paw_ij(1)%nspden ; nsppol=paw_ij(1)%nsppol
 end if

 !my_rank = xmpi_comm_rank(comm)

 natvshift=dtset%natvshift
 pawujat=dtset%pawujat
 natpawu=dtset%natpawu   ; ndtset=dtset%ndtset
 ABI_MALLOC(atvshift,(natvshift,nspden,natom))
 ABI_MALLOC(atvshmusk,(natvshift,nspden,natom))
 ABI_MALLOC(dmusk,(nspden,natom))
 ABI_MALLOC(nnocctot,(nspden,natom))
 musk=.false.; dmusk=.false.
 atvshift=fatvshift*dtset%atvshift
 atvshmusk=.false.
 nnocctot=zero
 typawujat=dtset%typat(pawujat)
 usepawu=(count(pawtab(:)%usepawu/=0)>0)

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 nnat=0
 if (usepawu) then
   write(message,'(3a)') ch10, '---------- pawuj_red ------ ',ch10
   call wrtout(std_out,  message,'COLL');
   do iatom=1,my_natom
     iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
     itypat=paw_ij(iatom)%itypat;ll=pawtab(itypat)%lpawu
     if ((ll>=0).and.(pawtab(itypat)%usepawu/=0).and.itypat==typawujat) then
       musk(:,iatom_tot)=(/.true., .true., .true. /)
       atvshmusk(:,:,iatom_tot)=reshape((/ (( (im1==1), im1=1,natvshift)  ,im2=1,nspden ) /),(/natvshift,nspden/))
       do ispden=1,nspden
         nnocctot(ispden,iatom_tot)=paw_ij(iatom)%nocctot(ispden)
         dmusk(ispden,iatom_tot)=.true.
       end do
       nnat=nnat+1
     end if
   end do

!  Reduction in case of parallelism
   if (paral_atom) then
     call xmpi_sum(nnocctot ,my_comm_atom,ierr)
     call xmpi_lor(atvshmusk,my_comm_atom) ! dim=natpawu ???
     call xmpi_lor(dmusk    ,my_comm_atom)
     call xmpi_lor(musk     ,my_comm_atom)
     call xmpi_sum(nnat     ,my_comm_atom,ierr)
   end if

   iuj=maxval(dtpawuj(:)%iuj)
   write(std_out,*)' pawuj_red: iuj',iuj
   write(std_out,*)' pawuj_red: dtpawuj(:)%iuj ',dtpawuj(:)%iuj

   !If this is the unperturbed state, then unscreened and screened occupancies
   !are the same. Also set perturbation vsh to zero.
   !If this is the unperturbed case, then do everything twice; once for iuj=1
   !and the same for iuj=2. If this is the perturbed case, do everything only
   !once. LMac
   if (iuj==1) then
     ABI_MALLOC(dtpawuj(0)%vsh,(nspden,nnat))
     ABI_MALLOC(dtpawuj(0)%occ,(nspden,nnat))
     ABI_MALLOC(dtpawuj(0)%xred,(3,nnat))
     dtpawuj(0)%vsh=0
     dtpawuj(0)%occ=0
     dtpawuj(0)%xred=0
     ncyc=2
   else
     ncyc=iuj
   end if

   do icyc=iuj,ncyc

     if (icyc==1) then  ! 1 and 3: non-scf steps
       dtpawuj(icyc+1)%iuj=icyc+1
       dtpawuj(icyc+2)%iuj=icyc+2
     else if (icyc==3) then
       dtpawuj(icyc+1)%iuj=icyc+1
     end if

     ! TODO: check that this is correct: this point is passed several times
     ! for a given value of iuj - should the stuff be accumulated instead of replaced?
     if(allocated(dtpawuj(icyc)%vsh))  then
       ABI_FREE(dtpawuj(icyc)%vsh)
     end if
     ABI_MALLOC(dtpawuj(icyc)%vsh,(nspden,nnat))
     !Allocate screened occupancy array: LMac
     if(allocated(dtpawuj(icyc)%occ))  then
       ABI_FREE(dtpawuj(icyc)%occ)
     end if
     ABI_MALLOC(dtpawuj(icyc)%occ,(nspden,nnat))
     if(allocated(dtpawuj(icyc)%xred))  then
       ABI_FREE(dtpawuj(icyc)%xred)
     end if
     ABI_MALLOC(dtpawuj(icyc)%xred,(3,nnat))


     rrtab=(/0.75_dp,0.815_dp,1.0_dp/)*pawtab(typawujat)%rpaw
     wftab=pawtab(typawujat)%phi(pawtab(typawujat)%mesh_size,pawtab(typawujat)%lnproju(1))

     do im1=1,ncoeff
       if (pawrad(typawujat)%mesh_type==1) then
         im2=nint(rrtab(im1)/pawrad(typawujat)%rstep+1)
       else if (pawrad(typawujat)%mesh_type==2) then
         im2=nint(log(rrtab(im1)/pawrad(typawujat)%rstep+1)/pawrad(typawujat)%lstep+1)
       else if (pawrad(typawujat)%mesh_type==3) then
         im2=nint(log(rrtab(im1)/pawrad(typawujat)%rstep)/pawrad(typawujat)%lstep+1)
       else if (pawrad(typawujat)%mesh_type==4) then
         im2=nint(pawtab(typawujat)%mesh_size*(1-exp((-one)*rrtab(im1)/pawrad(typawujat)%rstep))+1)
       end if

       rrtab(im1)=pawrad(typawujat)%rad(im2)
       wftab(im1)=pawtab(typawujat)%phi(im2,pawtab(typawujat)%lnproju(1))
     end do
     write(message,fmt='(a,i3,a,10f10.5)')' pawuj_red: mesh_type',pawrad(typawujat)%mesh_type,' rrtab:', rrtab
     call wrtout(std_out,message,'COLL')
     write(message,fmt='(a,10f10.5)')' pawuj_red: wftab', wftab
     call wrtout(std_out,message,'COLL')
     a=reshape((/ (( rrtab(im2)**(im1-1), im1=1,3)  ,im2=1,3 )/),(/ncoeff,ncoeff/))
     write(messg,fmt='(a)')'A'
     call linvmat(a,b,ncoeff,messg,2,0.0_dp,3) ! linvmat(inmat,oumat,nat,nam,option,gam,prtvol)
     write(std_out,*) 'pawuj_red: a,b ', a,b
     wftab=matmul(wftab,b)
     write(std_out,*) 'pawuj_red: wftab ', wftab
     dtpawuj(icyc)%wfchr(4:6)=wftab

     dtpawuj(icyc)%nat=nnat
     write(std_out,*) 'pawuj_red: m1'
     dtpawuj(icyc)%vsh=reshape(pack(atvshift,atvshmusk),(/ nspden,nnat /))
     !factor in next line to compensate nocctot contains just occ of 1 spin channel for nspden=1
     write(std_out,*) 'pawuj_red: m2'
     dtpawuj(icyc)%occ=reshape(pack(nnocctot,dmusk),(/nspden,nnat/))*(3-nspden)
     write(std_out,*) 'pawuj_red: m3'

     write(std_out,*) 'pawuj_red: occ ', dtpawuj(icyc)%occ

     dtpawuj(icyc)%xred=reshape(pack(dtset%xred_orig(:,:,1),musk),(/3,nnat/))
     dtpawuj(icyc)%ph0phiint=pawtab(typawujat)%ph0phiint(1)
     dtpawuj(icyc)%wfchr(1:3)=(/ pawtab(typawujat)%zioneff(1)*(dtset%lpawu(typawujat)+2),&
&     one*(dtset%lpawu(typawujat)+1),one*(dtset%lpawu(typawujat))/)
     dtpawuj(icyc)%pawrad=pawtab(typawujat)%rpaw

     write(std_out,*) 'pawuj_red: wfchr ',dtpawuj(icyc)%wfchr

     if (icyc.le.2) dtpawuj(icyc)%vsh=0.0d0

     write (hstr,'(I0)') icyc
     write(message,'(a,a,I3,a)') ch10, '---------- MARK ------ ',icyc,ch10
     call wrtout(std_out,message,'COLL')
     write(message,fmt='(a)') 'vsh'//trim(hstr)
     call wrtout(std_out,message,'COLL')
     call prmat(dtpawuj(icyc)%vsh(:,:),1,nnat*nspden,1)
     write(message,fmt='(a)') 'occ'//trim(hstr)
     call wrtout(std_out,message,'COLL')
     call prmat(dtpawuj(icyc)%occ(:,:),1,nnat*nspden,1)
     write(message, '(3a)' )'---------- MARK ---------- ',ch10
     call wrtout(std_out,message,'COLL')
   end do
 end if !usepawu

 ABI_FREE(nnocctot)
 ABI_FREE(dmusk)
 ABI_FREE(atvshift)
 ABI_FREE(atvshmusk)

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine pawuj_red
!!***

!----------------------------------------------------------------------

!!****f* m_paw_uj/chiscwrt
!! NAME
!!  chiscwrt
!!
!! FUNCTION
!!  PRIVATE ROUTINE
!!  Distributes values of chi_org on chi_sc according to ion-ion distances and multiplicities in shells
!!
!! INPUTS
!!  chi_org  = response matrix in original unit cell
!!  disv_org = distances (multiplied by magntization) in original unit cell
!!  nat_org = number of atoms in original unit cell
!!  sdisv_org = radii of shells in original unit cell
!!  smult_org = multiplicities of shells in original unit cell
!!  nsh_org   = number of shells in original unit cell
!!  disv_sc   = distances (multiplied by magntization) in super-cell
!!  nat_sc  = number of atoms in super-cell
!!  sdisv_sc = radii of shells in super-cell (was unused, so suppressed)
!!  smult_sc = multiplicities of shells in super-cell
!!  nsh_sc = number of shells in super-cell
!!  opt =
!!
!! OUTPUT
!!  chi_sc = first line of reponse matrix in super-cell
!!
!! SIDE EFFECTS
!!
!! SOURCE

subroutine chiscwrt(chi_org,disv_org,nat_org,sdisv_org,smult_org,nsh_org,chi_sc,&
& disv_sc,nat_sc,smult_sc,nsh_sc,opt,prtvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in)              :: nat_org,nsh_org,nsh_sc
 integer,intent(in)              :: nat_sc
 integer,intent(in),optional     :: prtvol
 integer,intent(in),optional     :: opt
!arrays
 real(dp),intent(in)             :: chi_org(nat_org),disv_org(nat_org),sdisv_org(nsh_org)
 integer,intent(in)              :: smult_org(nsh_org), smult_sc(nsh_sc)
 real(dp),intent(in)             :: disv_sc(nat_sc)
 real(dp),intent(out)            :: chi_sc(nat_sc)

!Local variables-------------------------------
!scalars
 integer                      :: iatom,jatom,jsh,optt
 character(len=500)           :: message
!arrays
 real(dp)                     :: chi_orgl(nat_org)

! *************************************************************************

 if (present(opt)) then
   optt=opt
 else
   optt=1
 end if


 do iatom=1,nat_org
   do jsh=1,nsh_org
     if (disv_org(iatom)==sdisv_org(jsh)) then
       if (opt==2) then
         chi_orgl(iatom)=chi_org(iatom)
       else if (opt==1) then
         chi_orgl(iatom)=chi_org(iatom)*smult_org(jsh)/smult_sc(jsh)
       end if
       exit
     end if
   end do !iatom
 end do  !jsh

 if (prtvol>1) then
   write(message,fmt='(a,150f10.5)')' chiscwrt: chi at input ',chi_org
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a,150f10.5)')' chiscwrt: chi after division ',chi_orgl
   call wrtout(std_out,message,'COLL')
 end if

 do iatom=1,nat_sc
   do jatom=1,nat_org
     if (disv_org(jatom)==disv_sc(iatom)) then
       chi_sc(iatom)=chi_orgl(jatom)
       exit
     else if (jatom==nat_org) then
       chi_sc(iatom)=0_dp
     end if
   end do
 end do

 if (prtvol>1) then
   write(message,'(a)')' chiscwrt, chi_sc '
   call wrtout(std_out,message,'COLL')
   call prmat(chi_sc,1,nat_sc,1,std_out)
 end if

end subroutine chiscwrt
!!***

!----------------------------------------------------------------------

!!****f* m_paw_uj/linvmat
!! NAME
!!  linvmat
!!
!! FUNCTION
!!  inverts real matrix inmat
!!
!! INPUTS
!!  inmat(1:nat,1:nat)=matrix to be inverted
!!  nat=dimension of inmat
!!  nam=comment specifiying the input matrix (to be printed in output)
!!  option=how to invert inmat
!!      option=1 or 3 add charge bath to matrix and add gam for inversion
!!      option=2 simply invert matrix
!!  gam=gamma added to inmat before inversion in case charge bath is used (allows inversion of otherwise
!!               singular matrix)
!!  prtvol=controls output to files (see subroutine lprtmat)
!!
!! OUTPUT
!!  oumat(nnat,nnat)=inverse of inmat, nnat=nat+1 for option=1 or option=3; nnat=nat for option=2
!!
!! SOURCE

subroutine linvmat(inmat,oumat,nat,nam,option,gam,prtvol)

!Arguments -------------------------------

 integer,intent(in)              :: nat
 real(dp),intent(in)             :: gam,inmat(nat,nat)
 real(dp),intent(inout)          :: oumat(:,:)         ! nat+1,nat+1 for option=1 or 2
                                                       ! nat,nat for option=2
 character(len=500),intent(in)   :: nam
 integer,intent(in),optional     :: prtvol,option

!Local variables -------------------------
 character(len=500)             :: message
 character(len=500)             :: bastrin,gastrin
 integer                        :: info,nnat,optionn
 integer,allocatable            :: ipvt(:)
 real(dp),allocatable           :: hma(:,:),work(:)

! *********************************************************************

 if (present(option)) then
   optionn=option
 else
   optionn=1
 end if

 if (option==1.or.option==3) then
   write(bastrin,'(a)')'+ charge bath '
   write(gastrin,'(a,d10.2,a)')'+ gamma  (=',gam,') '
 else
   write(bastrin,'(a)')''
   write(gastrin,'(a)')''
 end if

 write(message,fmt='(a)')' matrix '//trim(nam)
 call lprtmat(message,1,prtvol,inmat,nat)

 if (option==1.or.option==3) then
   call blow_pawuj(inmat,nat,oumat)
   write(message,fmt='(a,a)')' ',trim(nam)//trim(bastrin)
   call lprtmat(message,1,prtvol,oumat,nat+1)
   oumat=oumat+gam
   write(message,fmt='(a,a)')' ',trim(nam)//trim(bastrin) ! //trim(gastrin)
   call lprtmat(message,1,prtvol,oumat-gam,nat+1)
   nnat=nat+1
 else
   nnat=nat
   oumat=inmat
   oumat(1,1)=inmat(1,1)
 end if


 ABI_MALLOC(hma,(nnat,nnat))
 ABI_MALLOC(work,(nnat))
 ABI_MALLOC(ipvt,(nnat))
 work=0_dp
 hma(:,:)=oumat

 call dgetrf(nnat,nnat,hma,nnat,ipvt,info)
 if (.not.info==0) then
   write(message, '(3a)' ) 'Matrix '//trim(nam)//' is singular',ch10,'Probably too many symmetries kept'
   call wrtout(ab_out,message,'COLL')
   return
 end if

 call dgetri(nnat,hma,nnat,ipvt,work,nnat,info)
 oumat=hma(:,:)

 write(message,fmt='(2a,a)')' ('//trim(nam)//trim(bastrin)//trim(gastrin)//')^(-1)'
 call lprtmat(message,1,prtvol,oumat,nnat)

 ABI_FREE(hma)
 ABI_FREE(work)
 ABI_FREE(ipvt)

end subroutine linvmat
!!***

!----------------------------------------------------------------------

!!****f* m_paw_uj/lprtmat
!! NAME
!!  lprtmat
!!
!! FUNCTION
!!  prints out the real matrix mmat
!!
!! INPUTS
!!  mmat(nat,nat)=matrix to be printed
!!  nat=dimension of mmat
!!  prtvol specifies the volume of printing
!!   3: print the whole matrix
!!   2: print the first line
!!   1: do not print anything
!!  chan specifies the output files
!!   1: output only to std_out
!!   2: output also to ab_out
!!  commnt=comment specifying matirix
!!
!! OUTPUT
!!  oumat(nat+1,nat+1)=inverse of inmat
!!
!! SOURCE

subroutine lprtmat(commnt,chan,prtvol,mmat,nat)

!Arguments -------------------------------
 integer,intent(in)              :: nat,chan,prtvol
 real(dp),intent(in)             :: mmat(nat,nat)
 character(len=500),intent(in)  :: commnt

!Local variables -------------------------
 character(len=500)             :: message
! *********************************************************************

 if (prtvol==3) then
   write(message,fmt='(a)') trim(commnt)
   call wrtout(std_out,message,'COLL')
   call prmat(mmat,nat,nat,nat,std_out)
   if (chan==2) then
     call wrtout(ab_out,message,'COLL')
     call prmat(mmat,nat,nat,nat,ab_out)
   end if
   write(message,*)ch10
   call wrtout(std_out,message,'COLL')
   if (chan==2) then
     call wrtout(ab_out,message,'COLL')
   end if
 end if

 if (prtvol==2) then
   write(message,fmt='(a)') trim(commnt)
   call wrtout(std_out,message,'COLL')
   call prmat(mmat,1,nat,nat,std_out)
   if (chan==2) then
     call wrtout(ab_out,message,'COLL')
     call prmat(mmat,1,nat,nat,ab_out)
   end if
   write(message,*)ch10
   call wrtout(std_out,message,'COLL')
   if (chan==2) then
     call wrtout(ab_out,message,'COLL')
   end if
 end if

end subroutine lprtmat
!!***

!----------------------------------------------------------------------

!!****f* m_paw_uj/lcalcu
!! NAME
!!  lcalcu
!!
!! FUNCTION
!!  prints out real the real matrice mmat
!!
!! INPUTS
!!  magv=magnetic ordering of the ions (-1 of down/1 for up)
!!  natom=number of atoms
!!  rprimd(3,3)=lattic vectors of unit cell
!!  xred(3,natom)=positions of atoms
!!  chi(natom)=full response of atoms due to shift on atom pawujat
!!  chi0(natom)= response of atoms due to shift on atom pawujat
!!  pawujat=specifies on which atom the potential shift was done
!!  prtvol=controls output to files (see subroutine lprtmat)
!!  gam=gamma to be used for inversion of matrices (see subroutine livmat)
!!  opt=wether to use charge bath (1 or 3) or not (else)
!!
!! OUTPUT
!!  ures=resulting U (in eV) on atom pawujat
!!
!! SOURCE

subroutine lcalcu(magv,natom,rprimd,xred,chi,chi0,pawujat,ures,prtvol,gam,opt)

!Arguments ------------------------------------
!scalars
 integer,intent(in)          :: natom
 integer,intent(in),optional :: opt,pawujat,prtvol
 real(dp),intent(in),optional:: gam
 real(dp),intent(out)        :: ures
!arrays
 integer,intent(in)          :: magv(natom)
 real(dp),intent(in)         :: rprimd(3,3),xred(3,natom),chi(natom),chi0(natom)

!Local variables-------------------------------
!scalars
 character(len=500)          :: message
 integer                     :: optt,prtvoll,nnatom
 real(dp)                    :: gamm
!arrays
 real(dp),allocatable        :: tab(:,:,:)

! *********************************************************************

 if (present(opt)) then
   optt=opt
 else
   optt=1
 end if

 if (present(prtvol)) then
   prtvoll=prtvol
 else
   prtvoll=1
 end if

 if (present(gam)) then
   gamm=gam
 else
   gamm=1_dp
 end if

 if (optt==1.or.optt==3) then
   nnatom=natom+1
 else
   nnatom=natom
 end if

 ABI_MALLOC(tab,(4,nnatom,nnatom))

 call ioniondist(natom,rprimd,xred,tab(1,1:natom,1:natom),3,chi0,magv,pawujat,prtvoll)
 call ioniondist(natom,rprimd,xred,tab(2,1:natom,1:natom),3,chi,magv,pawujat,prtvoll)


 write(message,fmt='(a)')'response chi_0'
 call linvmat(tab(1,1:natom,1:natom),tab(3,1:nnatom,1:nnatom),natom,message,optt,gamm,prtvoll)
 call wrtout(std_out,message,'COLL')

 write(message,fmt='(a)')'response chi'
 call linvmat(tab(2,1:natom,1:natom),tab(4,1:nnatom,1:nnatom),natom,message,optt,gamm,prtvoll)
 call wrtout(std_out,message,'COLL')

 tab(1,1:nnatom,1:nnatom)=(tab(3,1:nnatom,1:nnatom)-tab(4,1:nnatom,1:nnatom))*Ha_eV

 write(message,fmt='(a,i3,a)')' (chi_0)^(-1)-(chi)^(-1) (eV)'
 call lprtmat(message,2,prtvoll,tab(1,1:nnatom,1:nnatom),nnatom)
 call wrtout(std_out,message,'COLL')

 ures=tab(1,1,pawujat)

 ABI_FREE(tab)

end subroutine lcalcu
!!***

!----------------------------------------------------------------------

!!****f* m_paw_uj/blow_pawuj
!!
!! NAME
!! blow_pawuj
!!
!! FUNCTION
!! This subroutine reads a real nxn matrice and appends lines n+1 and clumn n+1 containing
!! the sum of the lines
!!
!! INPUTS
!!  mat(nj,nj) matrix to be completed
!!
!! OUTPUT
!!  matt(nj+1,nj+1) completed matrix
!!
!! SOURCE

subroutine blow_pawuj(mat,nj,matt)

!Arguments ------------------------------------
!scalars
 integer,intent(in)        :: nj
!arrays
 real(dp),intent(in)       :: mat(nj,nj)
 real(dp),intent(out)      :: matt(nj+1,nj+1)

!Local variables-------------------------------
!scalars
 integer                   :: ii
!arrays

! *************************************************************************

 matt(1:nj,1:nj)=mat
 do  ii = 1,nj
   matt(ii,nj+1)=-sum(matt(ii,1:nj))
 end do

 do  ii = 1,nj+1
   matt(nj+1,ii)=-sum(matt(1:nj,ii))
 end do

end subroutine blow_pawuj
!!***

!----------------------------------------------------------------------

END MODULE m_paw_uj
!!***
