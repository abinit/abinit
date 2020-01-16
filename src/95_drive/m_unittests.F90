!!****m* ABINIT/m_unittests
!! NAME
!! m_unittests
!!
!! FUNCTION
!! Module to implement unit tests
!!
!! COPYRIGHT
!! Copyright (C) 1999-2019 ABINIT group (HM)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_unittests

 use defs_basis
 use m_errors
 use m_abicore
 use m_crystal
 use m_defs_ptgroups
 use m_ptgroups
 use m_tetrahedron
 use m_htetra
 use m_krank
 use m_hashtable
 use m_symkpt
 use m_sort

 use m_time,            only : cwtime, cwtime_report
 use m_numeric_tools,   only : linspace, ctrap
 use m_symtk,           only : matr3inv
 use m_io_tools,        only : open_file
 use m_kpts,            only : kpts_ibz_from_kptrlatt, listkk
 use m_geometry,        only : normv

 implicit none

 public :: tetra_unittests
 public :: kptrank_unittests
!!***

contains

!!****f* m_unittests/rprim_from_ptgroup
!! NAME
!!  rprim_from_ptgroup
!!
!! FUNCTION
!!  Get an rprim compatible with the point group
!!  This is only for testing porposes
!!
function rprim_from_ptgroup(ptgroup) result(rprim)
 character(len=*),intent(in) :: ptgroup
 real(dp) :: rprim(3,3)
 real(dp) :: tx,ty,tz
 real(dp) :: a,b,c
 real(dp) :: tmp1, tmp2
 real(dp) :: alpha,beta,gamma

 a = one
 b = two
 c = three
 alpha =     pi/two
 beta  =     pi/three
 gamma = two_pi/three

 select case(ptgroup)
 case ('1','-1')
   ! Triclinic
   tmp1 = c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
   tmp2 = c*sqrt( 1 + 2*cos(alpha)*cos(beta)*cos(gamma)-&
            cos(alpha)**2-cos(beta)**2-cos(gamma)**2 )/sin(gamma)
   rprim(:,1) = [a,                  0.0_dp, 0.0_dp]
   rprim(:,2) = [b*cos(gamma), b*sin(gamma), 0.0_dp]
   rprim(:,3) = [c*cos(beta),          tmp1,   tmp2]
 case ('2','m','-2','2/m')
   ! Monoclinic
   rprim(:,1) = [           a,       0.0_dp, 0.0_dp]
   rprim(:,2) = [b*cos(gamma), b*sin(gamma), 0.0_dp]
   rprim(:,3) = [      0.0_dp,       0.0_dp,      c]
 case ('222','mm2','mmm')
   ! Orthorhombic
   rprim(:,1) = [     a,0.0_dp,0.0_dp]
   rprim(:,2) = [0.0_dp,     b,0.0_dp]
   rprim(:,3) = [0.0_dp,0.0_dp,     c]
 case ('4','-4','4/m','422','4mm','-42m','4/mmm')
   ! Tetragonal
   rprim(:,1) = a*[1.0_dp, 0.0_dp, 0.0_dp]
   rprim(:,2) = a*[0.0_dp, 1.0_dp, 0.0_dp]
   rprim(:,3) = a*[0.0_dp, 0.0_dp,    c/a]
 case ('3','-3','32','3m','-3m')
   ! Trigonal
   tx = sqrt((1-c)/2)
   ty = sqrt((1-c)/6)
   tz = sqrt((1+2*c)/3)
   rprim(:,1) = a*[     tx,       -ty, tz]
   rprim(:,2) = a*[ 0.0_dp, 2.0_dp*ty, tz]
   rprim(:,3) = a*[    -tx,       -ty, tz]
 case ('6','-6','6/m','622','6mm','-62m','6/mmm')
   ! Hexagonal
   rprim(:,1) = a*[1.0_dp,        0.0_dp,  0.0_dp]
   rprim(:,2) = a*[ -half,sqrt(3.0_dp)/2,  0.0_dp]
   rprim(:,3) = a*[0.0_dp,        0.0_dp,     c/a]
 case ('23','m-3','432','-43m','m-3m')
   ! cubic
   rprim(:,1) = a*[1.0_dp, 0.0_dp, 0.0_dp]
   rprim(:,2) = a*[0.0_dp, 1.0_dp, 0.0_dp]
   rprim(:,3) = a*[0.0_dp, 0.0_dp, 1.0_dp]
 end select

end function rprim_from_ptgroup
!!***

!----------------------------------------------------------------------

!!****f* m_unittests/crystal_from_ptgroup
!! NAME
!!  crystal_from_ptgroup
!!
!! FUNCTION
!!  Create a crystal structure from a user defined point group
!!
type(crystal_t) function crystal_from_ptgroup(ptgroup) result(crystal)

!Arguments -------------------------------
 character(len=*),intent(in) :: ptgroup

!Variables -------------------------------
 type(irrep_t),allocatable :: irr(:)
 logical,parameter :: use_antiferro_true=.true.,remove_inv_false=.false.
 integer,parameter :: npsp1=1,space_group0=0,timrev1=1
 integer :: natom,nclass,ntypat,nsym
 integer :: typat(1)
 real(dp) :: rprimd(3,3)
 real(dp) :: amu(1),xred(3,1),znucl(1),zion(1)
 real(dp),allocatable :: tnons(:,:)
 character(len=5),allocatable :: class_names(:)
 integer,allocatable :: symrel(:,:,:),symafm(:),class_ids(:,:)

 rprimd = rprim_from_ptgroup(ptgroup)

 call get_point_group(ptgroup,nsym,nclass,symrel,class_ids,class_names,irr)
 ABI_FREE(class_ids)
 ABI_FREE(class_names)
 call irrep_free(irr)
 ABI_DT_FREE(irr)

 ABI_MALLOC(symafm,(nsym))
 symafm = 1
 ABI_MALLOC(tnons,(3,nsym))
 tnons = 0

 amu = 1
 natom = 1
 ntypat = 1
 typat = 1
 znucl=1
 zion=1
 xred(:,1) = [0,0,0]
 call crystal_init(amu,crystal,space_group0,natom,npsp1,ntypat,nsym,rprimd,typat,xred,&
                   zion,znucl,timrev1,use_antiferro_true,remove_inv_false,"test",&
                   symrel=symrel,symafm=symafm,tnons=tnons)
 ABI_FREE(symrel)
 ABI_FREE(tnons)
 ABI_FREE(symafm)

end function crystal_from_ptgroup
!!***

!----------------------------------------------------------------------

!!****f* m_unittests/tetra_unittests
!! NAME
!!  tetra_unittests
!!
!! FUNCTION
!!  Unit tests for the tetrahedron routines
!!
subroutine tetra_unittests(comm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: comm

!Local variables -------------------------
!scalars
 type(crystal_t) :: crystal
 integer,parameter :: timrev1=1,npsp1=1
 integer,parameter :: brav1=1,bcorr0=0,bcorr1=1,qptopt1=1,nqshft1=1,space_group0=0
 real(dp),parameter :: max_occ1=1.d0
 integer :: nqibz,iqibz,nqbz,ierr,nw
 real(dp) :: cpu, wall, gflops
 real(dp) :: dosdeltae, emin, emax, qnorm, dos_int
 character(len=80) :: errstr
 type(t_tetrahedron) :: tetraq
 type(htetra_t) :: htetraq
 integer :: in_qptrlatt(3,3),new_qptrlatt(3,3)
 integer,allocatable :: bz2ibz(:,:)
 real(dp) :: dos_qshift(3,nqshft1)
 real(dp) :: rlatt(3,3),qlatt(3,3)
 real(dp),allocatable :: tweight(:,:),dweight(:,:)
 real(dp),allocatable :: qbz(:,:),qibz(:,:)
 real(dp),allocatable :: wtq_ibz(:), wdt(:,:)
 real(dp),allocatable :: energies(:),eig(:),mat(:)
 real(dp),allocatable :: dos(:),idos(:)
 complex(dp),allocatable :: cenergies(:)
 complex(dp),allocatable :: cweight(:,:)

! *********************************************************************

 call cwtime(cpu, wall, gflops, "start")
 call wrtout(std_out,'DOS unit tests')
 !
 ! 0. Intialize
 !
 call wrtout(std_out,'0. Initialize')

 ! Create fake crystal
 crystal = crystal_from_ptgroup('m-3m')

 ! Create a regular grid
 in_qptrlatt(:,1) = [ 20, 0, 0]
 in_qptrlatt(:,2) = [ 0, 20, 0]
 in_qptrlatt(:,3) = [ 0, 0, 20]
 dos_qshift(:,1) = [0.0,0.0,0.0]
 call kpts_ibz_from_kptrlatt(crystal, in_qptrlatt, qptopt1, nqshft1, dos_qshift, &
                             nqibz, qibz, wtq_ibz, nqbz, qbz, new_kptrlatt=new_qptrlatt, bz2ibz=bz2ibz)
 call cwtime_report(" kpts_ibz_from_kptrlatt", cpu, wall, gflops)

 ! Initialize old tetrahedra
 rlatt = new_qptrlatt; call matr3inv(rlatt, qlatt)
 call init_tetra(bz2ibz(1,:), crystal%gprimd, qlatt, qbz, nqbz, tetraq, ierr, errstr, comm)
 call cwtime_report(" init_tetra ", cpu, wall, gflops)

 ! Initialize new tetrahedra
 call htetra_init(htetraq, bz2ibz(1,:), crystal%gprimd, qlatt, qbz, nqbz, qibz, &
                  nqibz, ierr, errstr, comm)
 call cwtime_report(" init_htetra", cpu, wall, gflops)

 !
 ! 1. Compute energies of a parabolic band
 !
 call wrtout(std_out,'1. Parabolic Band')
 ABI_MALLOC(eig,(nqibz))
 ABI_MALLOC(mat,(nqibz))
 do iqibz=1,nqibz
   qnorm = normv(qibz(:,iqibz),crystal%gmet,'R')
   ! The DOS for this function goes as sqrt(w-0.5)*two_pi
   eig(iqibz) = qnorm**2+half
   mat(iqibz) = abs(one/eig(iqibz))
   mat(iqibz) = 1!abs(one/eig(iqibz))
 end do

 ! Prepare DOS calculation
 emax = maxval(eig)+1.d0
 emin = 0.d0
 nw = 500
 dosdeltae = (emax-emin)/(nw-1)
 ABI_MALLOC(energies,(nw))
 ABI_MALLOC(dos,(nw))
 ABI_MALLOC(idos,(nw))
 ABI_MALLOC(wdt,(nw,2))
 energies = linspace(emin,emax,nw)
 call cwtime_report(" init", cpu, wall, gflops)

 ! Compute DOS using tetrahedron implementation
 ABI_CALLOC(tweight,(nw,nqibz))
 ABI_CALLOC(dweight,(nw,nqibz))
 call tetra_blochl_weights(tetraq,eig,emin,emax,max_occ1,nw,&
                           nqibz,bcorr0,tweight,dweight,comm)
 do iqibz=1,nqibz
   dweight(:,iqibz) = dweight(:,iqibz)*mat(iqibz)
   tweight(:,iqibz) = tweight(:,iqibz)*mat(iqibz)
 end do
 dos(:)  = sum(dweight,2)
 idos(:) = sum(tweight,2)
 call ctrap(nw,dos,dosdeltae,dos_int)
 write(std_out,*) "dos_int", dos_int, idos(nw)
 call cwtime_report(" tetra_blochl   ", cpu, wall, gflops)
 call write_file('parabola_tetra.dat', nw, energies, idos, dos)

 ! Compute blochl weights
 call htetraq%blochl_weights(eig,emin,emax,max_occ1,nw,&
                            nqibz,bcorr0,tweight,dweight,comm)
 do iqibz=1,nqibz
   dweight(:,iqibz) = dweight(:,iqibz)*mat(iqibz)
   tweight(:,iqibz) = tweight(:,iqibz)*mat(iqibz)
 end do
 dos(:)  = sum(dweight,2)
 idos(:) = sum(tweight,2)
 call ctrap(nw,dos,dosdeltae,dos_int)
 write(std_out,*) "dos_int", dos_int, idos(nw)
 call cwtime_report(" htetra_blochl   ", cpu, wall, gflops)
 call write_file('parabola_htetra.dat', nw, energies, idos, dos)

 ! Compute blochl weights
 call htetraq%blochl_weights(eig,emin,emax,max_occ1,nw,&
                            nqibz,1,tweight,dweight,comm)
 do iqibz=1,nqibz
   dweight(:,iqibz) = dweight(:,iqibz)*mat(iqibz)
   tweight(:,iqibz) = tweight(:,iqibz)*mat(iqibz)
 end do
 dos(:)  = sum(dweight,2)
 idos(:) = sum(tweight,2)
 call ctrap(nw,dos,dosdeltae,dos_int)
 write(std_out,*) "dos_int", dos_int, idos(nw)
 call cwtime_report(" htetra_blochl_corr   ", cpu, wall, gflops)
 call write_file('parabola_htetra_corr.dat', nw, energies, idos, dos)

 ! Compute weights using LV integration from TDEP
 call htetraq%blochl_weights(eig,emin,emax,max_occ1,nw,&
                              nqibz,2,tweight,dweight,comm)
 do iqibz=1,nqibz
   dweight(:,iqibz) = dweight(:,iqibz)*mat(iqibz)
   tweight(:,iqibz) = tweight(:,iqibz)*mat(iqibz)
 end do
 dos(:)  = sum(dweight,2)
 idos(:) = sum(tweight,2)
 call ctrap(nw,dos,dosdeltae,dos_int)
 write(std_out,*) "dos_int", dos_int, idos(nw)
 call cwtime_report(" htetra_blochl_lv   ", cpu, wall, gflops)
 call write_file('parabola_htetra_lv.dat', nw, energies, idos, dos)

 ! Compute DOS using new tetrahedron implementation
 ABI_MALLOC(cenergies,(nw))
 ABI_MALLOC(cweight,(nw,nqibz))
 cenergies = energies

 ! Use SIMTET routines
 call htetraq%weights_wvals_zinv(eig,nw,cenergies,max_occ1,&
                                 nqibz,1,cweight,comm)
 dos(:)  = -sum(aimag(cweight(:,:)),2)/pi
 idos(:) =  sum(real(cweight(:,:)),2)
 call ctrap(nw,dos,dosdeltae,dos_int)
 write(std_out,*) "dos_int", dos_int, idos(nw)
 call cwtime_report(" htetra_weights_wvals_zinv_simtet   ", cpu, wall, gflops)
 call write_file('parabola_zinv_simtet.dat', nw, energies, idos, dos)

 ! Use LV integration from TDEP
 call htetraq%weights_wvals_zinv(eig,nw,cenergies,max_occ1,&
                                 nqibz,2,cweight,comm)
 dos(:)  = -sum(aimag(cweight(:,:)),2)/pi
 idos(:) =  sum(real(cweight(:,:)),2)
 call ctrap(nw,dos,dosdeltae,dos_int)
 write(std_out,*) "dos_int", dos_int, idos(nw)
 call cwtime_report(" htetra_weights_wvals_zinv_lv   ", cpu, wall, gflops)
 call write_file('parabola_zinv_lv.dat', nw, energies, idos, dos)

 dos = zero; idos = zero
 do iqibz=1,nqibz
   call htetraq%get_onewk_wvals(iqibz,bcorr0,nw,energies,max_occ1,nqibz,eig,wdt)
   wdt(:,:) = wdt(:,:)*mat(iqibz)
   dos(:)  = dos(:)  + wdt(:,1)*wtq_ibz(iqibz)
   idos(:) = idos(:) + wdt(:,2)*wtq_ibz(iqibz)
 end do
 call ctrap(nw,dos,dosdeltae,dos_int)
 write(std_out,*) "dos_int", dos_int, idos(nw)
 call cwtime_report(" htetra_get_onewk_wvals", cpu, wall, gflops)
 call write_file('parabola_htetra_onewk_wvals.dat', nw, energies, idos, dos)

 dos = zero; idos = zero
 do iqibz=1,nqibz
   call htetraq%get_onewk(iqibz,bcorr0,nw,nqibz,eig,emin,emax,max_occ1,wdt)
   wdt(:,:) = wdt(:,:)*mat(iqibz)
   dos(:)  = dos(:)  + wdt(:,1)*wtq_ibz(iqibz)
   idos(:) = idos(:) + wdt(:,2)*wtq_ibz(iqibz)
 end do
 call ctrap(nw,dos,dosdeltae,dos_int)
 write(std_out,*) "dos_int", dos_int, idos(nw)
 call cwtime_report(" htetra_get_onewk", cpu, wall, gflops)
 call write_file('parabola_htetra_onewk.dat', nw, energies, idos, dos)
 call htetraq%print()

 !
 ! 2. Compute energies for a flat band
 !
 call wrtout(std_out,'2. Flat Band')
 eig = 0.5

 ! Compute DOS using old tetrahedron implementation
 call tetra_blochl_weights(tetraq,eig,emin,emax,max_occ1,nw,&
                           nqibz,bcorr0,tweight,dweight,comm)
 dos(:)  = sum(dweight,2)
 idos(:) = sum(tweight,2)
 call cwtime_report(" tetra_blochl", cpu, wall, gflops)
 call write_file('flat_tetra.dat', nw, energies, idos, dos)

 ! Compute DOS using new tetrahedron implementation
 call htetraq%blochl_weights(eig,emin,emax,max_occ1,nw,&
                             nqibz,bcorr0,tweight,dweight,comm)
 dos(:)  = sum(dweight,2)
 idos(:) = sum(tweight,2)
 call cwtime_report(" htetra_blochl", cpu, wall, gflops)
 call write_file('flat_htetra.dat', nw, energies, idos, dos)

 dos = zero; idos = zero
 do iqibz=1,nqibz
   call htetraq%get_onewk_wvals(iqibz,bcorr0,nw,energies,max_occ1,nqibz,eig,wdt)
   dos(:)  = dos(:)  + wdt(:,1)*wtq_ibz(iqibz)
   idos(:) = idos(:) + wdt(:,2)*wtq_ibz(iqibz)
 end do
 call cwtime_report(" htetra_get_onewk_wvals", cpu, wall, gflops)
 call write_file('flat_htetra_onewk_wvals.dat', nw, energies, idos, dos)

 dos = zero; idos = zero
 do iqibz=1,nqibz
   call htetraq%get_onewk(iqibz,bcorr0,nw,nqibz,eig,emin,emax,max_occ1,wdt)
   dos(:)  = dos(:)  + wdt(:,1)*wtq_ibz(iqibz)
   idos(:) = idos(:) + wdt(:,2)*wtq_ibz(iqibz)
 end do
 call cwtime_report(" htetra_get_onewk", cpu, wall, gflops)
 call write_file('flat_htetra_onewk.dat', nw, energies, idos, dos)

 !
 ! 3. Compute tetrahedron for simple TB FCC
 !
 !TODO

 ! Free memory
 ABI_SFREE(energies)
 ABI_SFREE(eig)
 ABI_SFREE(mat)
 ABI_SFREE(wdt)
 ABI_SFREE(tweight)
 ABI_SFREE(dweight)
 ABI_SFREE(wtq_ibz)
 ABI_SFREE(dos)
 ABI_SFREE(idos)
 ABI_SFREE(qbz)
 ABI_SFREE(qibz)
 ABI_SFREE(bz2ibz)
 ABI_SFREE(cenergies)
 ABI_SFREE(cweight)

 call crystal%free()
 call htetraq%free()
 call destroy_tetra(tetraq)

 contains

 subroutine write_file(fname,nw,energies,idos,dos)
  integer,intent(in) :: nw
  character(len=*), intent(in) :: fname
  real(dp),intent(in) :: energies(nw),dos(nw),idos(nw)
  character(len=500) :: msg
  integer :: iw,funit

  if (open_file(fname,msg,newunit=funit,action="write") /= 0) then
    MSG_ERROR(msg)
  end if

  write(funit, *)"# Energies, IDOS, DOS"
  do iw=1,nw
    write(funit,*) energies(iw), idos(iw), dos(iw)
  end do
  close(unit=funit)

 end subroutine write_file

end subroutine tetra_unittests
!!***

!!****f* m_unittests/kptrank_unittests
!!
!! NAME
!!  kptrank_unittests
!!
!! FUNCTION
!!  Test the krank routines
!!
subroutine kptrank_unittests(comm)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: comm

!Variables -------------------------------
 type(crystal_t) :: crystal
 type(krank_t) :: krank
 integer,parameter :: qptopt1=1,nqshft1=1,iout0=0,chksymbreak0=0,sppoldbl1=1
 integer :: nqibz,iqbz,iqibz,iqbz_rank,nqbz,nqibz_symkpt,nqibz_symkpt_new
 integer :: in_qptrlatt(3,3),new_qptrlatt(3,3)
 real(dp) :: cpu, gflops, wall
 real(dp) :: dksqmax
 real(dp) :: dos_qshift(3,nqshft1)
 character(len=500) :: msg
 integer,allocatable :: bz2ibz(:,:)
 integer,allocatable :: bz2ibz_symkpt(:,:), bz2ibz_symkpt_new(:,:)
 integer,allocatable :: bz2ibz_listkk(:,:)
 integer,allocatable :: ibz2bz(:), ibz2bz_new(:)
 real(dp),allocatable :: wtq_fullbz(:), wtq_folded(:)
 real(dp),allocatable :: wtq_ibz(:)
 real(dp),allocatable :: qbz(:,:),qibz(:,:)

 ABI_UNUSED(comm)

 ! Create fake crystal
 !crystal = crystal_from_ptgroup('6')
 crystal = crystal_from_ptgroup('m-3m')

 ! Create a regular grid
 in_qptrlatt(:,1)=[ 100, 0, 0]
 in_qptrlatt(:,2)=[ 0, 100, 0]
 in_qptrlatt(:,3)=[ 0, 0, 100]
 dos_qshift(:,1) =[0.0,0.0,0.0]
 call cwtime(cpu, wall, gflops, "start")
 call kpts_ibz_from_kptrlatt(crystal, in_qptrlatt, qptopt1, nqshft1, dos_qshift, &
                             nqibz, qibz, wtq_ibz, nqbz, qbz, new_kptrlatt=new_qptrlatt, bz2ibz=bz2ibz)
 call cwtime_report(" kpts_ibz_from_kptrlatt", cpu, wall, gflops)

 ! Test krank
 krank = krank_new(nqbz, qbz)
 do iqbz=1,nqbz
   iqbz_rank = krank%get_index(qbz(:,iqbz))
   ABI_CHECK(iqbz==iqbz_rank,'wrong q-point')
 end do
 call cwtime_report(" krank", cpu, wall, gflops)
 call krank%free()

 ABI_MALLOC(wtq_fullbz,(nqbz))
 ABI_MALLOC(wtq_folded,(nqbz))
 ABI_MALLOC(ibz2bz,(nqbz))
 ABI_MALLOC(ibz2bz_new,(nqbz))
 ABI_MALLOC(bz2ibz_symkpt,(6,nqbz))
 ABI_MALLOC(bz2ibz_symkpt_new,(6,nqbz))
 wtq_fullbz = one/nqbz

 ! Test symkpt (note that the above call to kpts_ibz_from_kptrlatt already involves calling this routine)
 call symkpt(chksymbreak0,crystal%gmet,ibz2bz,iout0,qbz,nqbz,&
             nqibz_symkpt,crystal%nsym,crystal%symrec,crystal%timrev,&
             wtq_fullbz,wtq_folded, bz2ibz_symkpt, comm)
 call cwtime_report(" symkpt", cpu, wall, gflops)

 wtq_fullbz = one/nqbz
 call symkpt_new(chksymbreak0,crystal%gmet,ibz2bz_new,iout0,qbz,nqbz,&
                 nqibz_symkpt_new,crystal%nsym,crystal%symrec,crystal%timrev,&
                 bz2ibz_symkpt_new, comm)
 call cwtime_report(" symkpt_new", cpu, wall, gflops)
 ABI_CHECK(nqibz_symkpt==nqibz_symkpt_new,'Wrong number of qpoints in the IBZ')

 ! check if ibz is the same
 do iqibz=1,nqibz
   if (ibz2bz(iqibz) == ibz2bz_new(iqibz)) cycle
   MSG_ERROR("The IBZ is different")
 end do

 ! check if mapping is the same
 do iqbz=1,nqbz
   if (bz2ibz_symkpt(1,iqbz) == bz2ibz_symkpt_new(1,iqbz)) cycle
   write(msg,*) "Inconsistent mapping", iqbz, bz2ibz_symkpt(1,iqbz), bz2ibz_symkpt_new(1,iqbz)
   MSG_ERROR(msg)
 end do

 ! call listkk
 ABI_MALLOC(bz2ibz_listkk,(nqbz,6))
 call listkk(dksqmax,crystal%gmet,bz2ibz_listkk,qibz,qbz,nqibz,nqbz,crystal%nsym,&
             sppoldbl1,crystal%symafm,crystal%symrec,crystal%timrev,comm,&
             exit_loop=.True., use_symrec=.True.)
 call cwtime_report(" listkk", cpu, wall, gflops)

 ! check if indkk is the same
 do iqbz=1,nqbz
   if (bz2ibz_listkk(iqbz,1) /= bz2ibz_symkpt_new(1,iqbz)) then
     write(std_out,*) "Inconsistent ikpt", iqbz, bz2ibz_listkk(iqbz,1), bz2ibz_symkpt_new(1,iqbz)
   end if
   if (bz2ibz_listkk(iqbz,2) /= bz2ibz_symkpt_new(2,iqbz)) then
     write(std_out,*) "Inconsistent isym", iqbz, bz2ibz_listkk(iqbz,2), bz2ibz_symkpt_new(2,iqbz)
   end if
   if (bz2ibz_listkk(iqbz,6) /= bz2ibz_symkpt_new(3,iqbz)) then
     write(std_out,*) "Inconsistent itim", iqbz, bz2ibz_listkk(iqbz,6), bz2ibz_symkpt_new(3,iqbz)
   end if
   if (.not.all(bz2ibz_listkk(iqbz,3:5) == bz2ibz_symkpt_new(4:,iqbz))) then
     write(std_out,*) "Inconsistent shift", iqbz
     write(std_out,*) bz2ibz_listkk(iqbz,3:5)
     write(std_out,*) bz2ibz_symkpt_new(4:,iqbz)
   end if
 end do

 ABI_SFREE(bz2ibz_symkpt)
 ABI_SFREE(bz2ibz_symkpt_new)
 ABI_SFREE(bz2ibz_listkk)
 ABI_SFREE(wtq_fullbz)
 ABI_SFREE(wtq_folded)
 ABI_SFREE(ibz2bz)
 ABI_SFREE(ibz2bz_new)
 ABI_SFREE(qbz)
 ABI_SFREE(qibz)
 ABI_SFREE(bz2ibz)
 ABI_SFREE(wtq_ibz)
 call crystal%free()

end subroutine kptrank_unittests
!!***

end module m_unittests
