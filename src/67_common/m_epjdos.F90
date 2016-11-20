!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_epjdos
!! NAME
!!  m_epjdos
!!
!! FUNCTION
!!  Tools for the computiation of electronic PJDOSes
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2016 ABINIT group (MVer, XG, SM, MT, BAmadon, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_epjdos

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_tetrahedron
 use m_splines
 use m_cgtools
 use m_atomdata
 use m_crystal
 use m_crystal_io
 use m_ebands
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr

 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_time,           only : cwtime
 use m_io_tools,       only : open_file
 use m_numeric_tools,  only : simpson, simpson_int
 use m_fstrings,       only : int2char4
 use m_pawtab,         only : pawtab_type

 implicit none

 private
!!***

 public :: tetrahedron         ! Calculate DOS by tetrahedron method.
 public :: gaus_dos            ! Calculate DOS by gauss method.
 public :: recip_ylm           ! Project input wavefunctions (real space) on to Ylm
 public :: dens_in_sph         ! Calculate integrated density in sphere around each atom
!!***

!----------------------------------------------------------------------

!!****t* m_epjdos/epjdos_t
!! NAME
!! epjdos_t
!!
!! FUNCTION
!!  Stores different contributions to the electronic DOS.
!!
!! NOTES
!!  Please contact gmatteo if you plan to change the internal implementation
!!  or add new DOSes. These results are saved in a netcdf file (see fatbands_ncwrite)
!!  so that one can read it with python and plot fatbands and PJDOSEs.
!!  The python version is able to handle the different cases (L, LM, Spin ...) but
!!  any change in the internal Abinit implementation is likely to break the python interface.
!!
!! SOURCE

 type,public :: epjdos_t

   integer :: mbesslang
   ! Max L+1 used in LM-DOS  (Bessel function expansion)

   integer :: ndosfraction
   ! Defines the last dimension of the dos arrays.
   ! Actual value depends on the other variables.

   integer :: prtdos
   ! 2 --> Standard DOS with tetra.
   ! 3 --> L-DOS with tetra (prtdosm>0 if LM is wanted in Ylm/Slm basis).
   ! 4 --> L-DOS with gaussian (prtdosm if LM is wanted in Ylm/Slm basis).
   ! 5 --> Spin-DOS

   integer :: prtdosm
   ! Option for the m-contributions to the partial DOS
   ! 1 if LM-projection is done onto complex Ylm
   ! 2 if LM-projection is done onto real Slm

   integer :: partial_dos_flag

   integer :: paw_dos_flag
   ! 1 if both PAW contributions are evaluated AND stored

   !integer :: pawfatbnd
   integer :: fatbands_flag

   integer :: nkpt, mband, nsppol
   ! Used to dimension arrays

   integer,allocatable :: mlang_type(:)
   ! mlang_type(ntypat + natsph_extra)
   ! Max L+1 used in LM-DOS for each atom type

   real(dp),allocatable :: fractions(:,:,:,:)
   ! fractions(nkpt,mband,nsppol,ndosfraction))

   real(dp),allocatable :: fractions_m(:,:,:,:)
   ! fractions_m(nkpt,mband,nsppol,ndosfraction*mbesslang)

   !real(dp),allocatable :: fractions_average_m(:,:,:,:)
   ! fractions_average_m(nkpt,mband,nsppol,ndosfraction*mbesslang)

   real(dp),allocatable :: fractions_paw1(:,:,:,:)
   ! fractions_paw1(nkpt,mband,nsppol,ndosfraction)

   real(dp),allocatable :: fractions_pawt1(:,:,:,:)
   ! fractions_pawt1(nkpt,mband,nsppol,ndosfraction))

 end type epjdos_t

 public :: epjdos_new
 public :: epjdos_free

 public :: prtfatbands
 public :: fatbands_ncwrite

!----------------------------------------------------------------------

contains  !============================================================
!!***

!!****f* m_epjdos/epjdos_new
!! NAME
!!  epjdos_new
!!
!! FUNCTION
!!  Create new object from dataset input variables.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! PARENTS
!!
!! SOURCE

type(epjdos_t) function epjdos_new(dtset, psps, pawtab) result(new)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'epjdos_new'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: ierr,itypat,iat

! *********************************************************************

 new%nkpt = dtset%nkpt; new%mband = dtset%mband; new%nsppol = dtset%nsppol

 new%prtdos = dtset%prtdos
 new%partial_dos_flag = 0
 if (new%prtdos==2) new%partial_dos_flag = 0 ! Standard DOS with tetra.
 if (new%prtdos==3) new%partial_dos_flag = 1 ! L-DOS with tetra (prtdosm>0 if LM is wanted in Ylm/Slm basis).
 if (new%prtdos==4) new%partial_dos_flag = 1 ! L-DOS with gaussian (prtdosm if LM is wanted in Ylm/Slm basis).
 if (new%prtdos==5) new%partial_dos_flag = 2 ! Spin DOS

 new%prtdosm=0
 if (new%partial_dos_flag==1) new%prtdosm=dtset%prtdosm
 ! paw_dos_flag= 1 if both PAW contributions are evaluated AND stored
 new%paw_dos_flag=0
 if (dtset%usepaw==1 .and. new%partial_dos_flag==1 .and. dtset%pawprtdos==1) new%paw_dos_flag=1

 new%fatbands_flag=0
 if (dtset%pawfatbnd>0 .and. new%prtdosm==0) new%fatbands_flag=1
 if (new%prtdosm==1.and.dtset%pawfatbnd>0)then
!  because they compute quantities in real and complex harmonics respectively
   MSG_ERROR('pawfatbnd>0  and prtdosm=1 are not compatible')
 end if

 ! mjv : initialization is needed as mbesslang is used for allocation below
 ! NOTE: 10/5/2010 the whole of this could be looped over ndosfraction,
 ! to store much less in memory. The DOS is accumulated in an array
 ! and then printed to file at the end.
 new%mbesslang = 1
 if (new%partial_dos_flag==1 .or. new%fatbands_flag==1) then

   ABI_MALLOC(new%mlang_type, (dtset%ntypat + dtset%natsph_extra))
   new%mlang_type = 0

   ! TODO: Could use mbesslang = 4 or compute it from psps/pawtab
   ! Increment by one (could underestimate if vloc = vlmax)
   if (dtset%usepaw == 0) then
     do iat=1,dtset%natsph
       itypat = dtset%typat(dtset%iatsph(iat))
       new%mlang_type(itypat) = 1 + maxval(psps%indlmn(1, :, itypat))
     end do
   else
     do iat=1,dtset%natsph
       itypat= dtset%typat(dtset%iatsph(iat))
       new%mlang_type(itypat) = 1 + (pawtab(itypat)%l_size - 1) / 2
     end do
   end if

   ! Up to l=g if we have natsph_extra.
   if (dtset%natsph_extra > 0) new%mlang_type(dtset%ntypat+1:) = 5

   new%mlang_type = 5  ! This is to preserve the old implementation
   new%mbesslang = maxval(new%mlang_type)
   new%ndosfraction = (dtset%natsph + dtset%natsph_extra) * new%mbesslang

 else if (new%partial_dos_flag == 2) then
   new%ndosfraction = 7

 else
   new%ndosfraction = 1
   new%mbesslang = 0
 end if

 ! Check allocations status as these arrays are not distributed and the wavefunctions are still in memory.
 ABI_STAT_MALLOC(new%fractions, (dtset%nkpt,dtset%mband,dtset%nsppol,new%ndosfraction), ierr)
 ABI_CHECK(ierr==0, "out of memory in new%fractions")
 new%fractions = zero

 if (new%prtdosm>=1 .or. new%fatbands_flag==1) then
   ABI_STAT_MALLOC(new%fractions_m,(dtset%nkpt,dtset%mband,dtset%nsppol,new%ndosfraction*new%mbesslang), ierr)
   ABI_CHECK(ierr==0, "out of memory in new%fractions_m")
   new%fractions_m = zero
 end if

 if (dtset%usepaw==1 .and. new%partial_dos_flag==1) then
   ABI_STAT_MALLOC(new%fractions_paw1,(dtset%nkpt,dtset%mband,dtset%nsppol,new%ndosfraction), ierr)
   ABI_CHECK(ierr==0, "out of memory in new%fraction_paw1")
   ABI_STAT_MALLOC(new%fractions_pawt1,(dtset%nkpt,dtset%mband,dtset%nsppol,new%ndosfraction), ierr)
   ABI_CHECK(ierr==0, "out of memory in new%fraction_pawt1")
   new%fractions_paw1 = zero; new%fractions_pawt1 = zero
 end if

end function epjdos_new
!!***

!!****f* m_epjdos/epjdos_free
!! NAME
!!  epjdos_free
!!
!! FUNCTION
!!  deallocate memory
!!
!! PARENTS
!!
!! SOURCE

subroutine epjdos_free(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'epjdos_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(epjdos_t),intent(inout) :: self

! *********************************************************************

 ! integer
 if (allocated(self%mlang_type)) then
   ABI_FREE(self%mlang_type)
 end if

 ! real
 if (allocated(self%fractions)) then
   ABI_FREE(self%fractions)
 end if
 if (allocated(self%fractions_m)) then
   ABI_FREE(self%fractions_m)
 end if
 if (allocated(self%fractions_paw1)) then
   ABI_FREE(self%fractions_paw1)
 end if
 if (allocated(self%fractions_pawt1)) then
   ABI_FREE(self%fractions_pawt1)
 end if

end subroutine epjdos_free
!!***

!!****f* m_epjdos/tetrahedron
!! NAME
!! tetrahedron
!!
!! FUNCTION
!! calculate DOS by tetrahedron method
!!
!! INPUTS
!!  dos_fractions= projections of wavefunctions on each angular momentum Ylm
!!     which is the weight going into the DOS for an l-decomposed dos
!!  dos_fractions_m= same as dos_fractions, but m-decomposed not just l-
!!  dos_fractions_paw1= contribution to dos fractions from the PAW partial waves (phi)
!!  dos_fractions_pawt1= contribution to dos fractions from the PAW pseudo partial waves (phi_tild)
!!  dtset     structured datatype, from which one uses :
!!   kpt(3,nkpt)  =irreducible kpoints
!!   kptrlatt(3,3)=lattice vectors for full kpoint grid
!!   nband        =number of electronic bands for each kpoint
!!   nkpt         =number of irreducible kpoints
!!   nshiftk      =number of kpoint grid shifts
!!   nsym         =number of symmetries
!!   pawprtdos    =option to output the individual contributions to the partial DOS (0, 1 or 2)
!!   shiftk(3,nshiftk)=kpoint shifts
!!   symrel(3,3,nsym)=symmetry matrices in real space
!!   typat(ntypat)=types of atoms
!!   usepaw       =option for PAW
!!  crystal<crystal_t>=Object defining the unit cell and its symmetries.
!!  ebands<ebands_t>=Band structure data.
!!  fermie=Fermi energy
!!  fildata=name of the DOS output file
!!  mbesslang=maximum angular momentum for Bessel function expansion
!!  prtdosm=option for the m-contributions to the partial DOS
!!  ndosfraction= number of types of DOS we are calculating, e.g. the number
!!    of l channels. Could be much more general, for other types of partial DOS
!!  paw_dos_flag= option for partial dos in PAW
!!
!! OUTPUT
!!  (no explicit output)
!!  note: could have routine return DOS for other purposes, and separate printing
!!  in another routine (MJV 8/2010)
!!
!! NOTES
!!  This routine should be called by master only
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      destroy_tetra,dos_hdr_write,get_dos_1band,get_dos_1band_m
!!      get_full_kgrid,get_tetra_weight,init_tetra,int2char4,matr3inv,metric
!!      wrtout
!!
!! SOURCE

subroutine tetrahedron(dos,dtset,crystal,ebands,fildata)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tetrahedron'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_56_recipspace
 use interfaces_61_occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: fildata
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: crystal
 type(ebands_t),intent(in) :: ebands
 type(epjdos_t),intent(in) :: dos

!Local variables-------------------------------
!scalars
 integer,parameter :: bcorr0=0
 integer :: iat,iband,iene,ifract,ikpt,isppol,natsph,natsph_extra,nkpt,nsppol,i1,i2
 integer :: nene,nkpt_fullbz,prtdos,unitdos,ierr,prtdosm,paw_dos_flag,mbesslang,ndosfraction
 real(dp),parameter :: dos_max=9999.9999_dp
 real(dp) :: buffer,deltaene,enemax,enemin,enex,integral_DOS,max_occ,rcvol
 real(dp) :: cpu,wall,gflops
 logical :: bigDOS
 character(len=10) :: tag
 character(len=500) :: frmt,frmt_extra,message
 character(len=fnlen) :: tmpfil
 character(len=80) :: errstr
 type(t_tetrahedron) :: tetrahedra
!arrays
 integer,allocatable :: indkpt(:),unitdos_arr(:)
 real(dp) :: gprimd(3,3),klatt(3,3),rlatt(3,3)
 real(dp),allocatable :: dtweightde(:,:),integ_dos(:,:),integ_dos_m(:,:)
 real(dp),allocatable :: kpt_fullbz(:,:),partial_dos(:,:)
 real(dp),allocatable :: partial_dos_m(:,:),tmp_eigen(:),total_dos(:,:)
 real(dp),allocatable :: total_dos_m(:,:),total_dos_paw1(:,:)
 real(dp),allocatable :: total_dos_pawt1(:,:),total_integ_dos(:,:)
 real(dp),allocatable :: total_integ_dos_m(:,:),tweight(:,:)
 real(dp),allocatable :: work_ndos(:,:),work_ndosmbessl(:,:)

! *********************************************************************
 prtdosm = dos%prtdosm; paw_dos_flag = dos%paw_dos_flag
 mbesslang = dos%mbesslang; ndosfraction = dos%ndosfraction

 nkpt = dtset%nkpt; nsppol = dtset%nsppol

!m-decomposed DOS not compatible with PAW-decomposed DOS
 if(prtdosm>=1.and.paw_dos_flag==1) then
   message = 'm-decomposed DOS (prtdosm>=1) not compatible with PAW-decomposed DOS (pawprtdos=1) !'
   MSG_ERROR(message)
 end if

!Refuse only 1 kpoint: the algorithms are no longer valid. DOH !
 if (nkpt == 1) then
   MSG_WARNING('You need at least 2 kpoints to use the tetrahedron method. DOS cannot be computed')
   return
 end if

!Do not support nshiftk > 1: lattice must be decomposed into boxes
!and this is not always possible (I think) with bizzare shiftks
!normally at this point we have incorporated everything into
!kptrlatt, and only 1 shift is needed (in particular for MP grids).
 if (dtset%nshiftk > 1) then
   write(std_out,*) 'tetrahedron: skip subroutine.'
   write(std_out,*) 'Problem with a composite k-point grid.'
   write(std_out,*) 'Only simple lattices are supported. Action: use nshiftk=1.'
   write(std_out,*) 'dtset%nshiftk, dtset%shiftk = ', dtset%nshiftk, dtset%shiftk
   write(std_out,*) 'dtset%kptrlatt= ', dtset%kptrlatt
   MSG_WARNING('tetrahedron: skip subroutine. See message above')
   return
 end if

!Refuse nband different for different kpoints
!Note: This means we can pass ebands%eig(:,:,:) instead of eigen(mband*nkpt*nsppol) in packed form
 do isppol=1,nsppol
   do ikpt=1,nkpt
     if ( dtset%nband(nkpt*(isppol-1) + ikpt) /= dtset%nband(1) ) then
       write(std_out,*) 'tetrahedron: skip subroutine.'
       write(std_out,*) 'nband must be the same for all kpoints'
       write(std_out,*) 'nband=', dtset%nband
       MSG_WARNING('tetrahedron: skip subroutine. See message above')
       return
     end if
   end do
 end do

 call cwtime(cpu, wall, gflops, "start")

 ! TODO
 !call tetra_from_kptrlatt(tetra, crystal, dtset%kptopt, dtset%kptrlatt, &
 !                         dtset%nshiftk, dtset%shiftk, dtset%nkpt, dtset%kpt)

! Calculate nkpt_fullbz
 nkpt_fullbz= dtset%kptrlatt(1,1)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,3) &
& +dtset%kptrlatt(1,2)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,1) &
& +dtset%kptrlatt(1,3)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,2) &
& -dtset%kptrlatt(1,2)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,3) &
& -dtset%kptrlatt(1,3)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,1) &
& -dtset%kptrlatt(1,1)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,2)
 nkpt_fullbz = nkpt_fullbz*dtset%nshiftk

 if (nkpt_fullbz==0) then
   write(std_out,*)'tetrahedron: skip subroutine.'
   write(std_out,*)'no homogeneous grid  of k-points is defined ...'
   write(std_out,*)'in order to obtain the DOS using the tetrahedron method,'
   write(std_out,*)'you need to re-define ngkpt or kptrlatt.'
   MSG_WARNING('tetrahedron: skip subroutine. See message above')
   return
 end if

!Make klatt
 rlatt = dtset%kptrlatt(:,:)
 call matr3inv(rlatt,klatt)

!Get metric tensors
 gprimd = crystal%gprimd
 rcvol = abs (gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
& -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
& +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))

 ABI_ALLOCATE(indkpt,(nkpt_fullbz))
 ABI_ALLOCATE(kpt_fullbz,(3,nkpt_fullbz))

!Make full kpoint grid and get equivalence to irred kpoints
 call get_full_kgrid(indkpt,dtset%kpt,kpt_fullbz,dtset%kptrlatt,&
& nkpt,nkpt_fullbz,dtset%nshiftk,dtset%nsym,dtset%shiftk,dtset%symrel)

!Get tetrahedra, ie indexes of the full kpoints at their summits
 call init_tetra (indkpt,gprimd,klatt,kpt_fullbz,nkpt_fullbz,tetrahedra, ierr, errstr)
 ABI_CHECK(ierr==0,errstr)

 natsph=dtset%natsph; natsph_extra=dtset%natsph_extra

 ! Open the DOS file
 if (dtset%prtdos == 2 .or. dtset%prtdos == 5) then
   if (open_file(fildata,message,newunit=unitdos,status='unknown',form='formatted') /= 0) then
     MSG_ERROR(message)
   end if

 else if (dtset%prtdos == 3) then
   ABI_ALLOCATE(unitdos_arr,(natsph+natsph_extra))
   do iat=1,natsph
     call int2char4(dtset%iatsph(iat),tag)
     ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
     tmpfil = trim(fildata)//'_AT'//trim(tag)
     if (open_file(tmpfil, message, newunit=unitdos_arr(iat), status='unknown',form='formatted') /= 0) then
       MSG_ERROR(message)
     end if
   end do
!  do extra spheres in vacuum too. Use _ATEXTRA[NUM] suffix
   do iat=1,natsph_extra
     call int2char4(iat,tag)
     ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
     tmpfil = trim(fildata)//'_ATEXTRA'//trim(tag)
     if (open_file(tmpfil, message, newunit=unitdos_arr(natsph+iat), status='unknown',form='formatted') /= 0) then
       MSG_ERROR(message)
     end if
   end do
 end if

!Write the header of the DOS file, and determine the energy range and spacing
 prtdos=dtset%prtdos
 buffer=0.01_dp ! Size of the buffer around the min and max ranges
 if (dtset%prtdos == 2 .or. dtset%prtdos == 5) then
   call dos_hdr_write(buffer,deltaene,dtset%dosdeltae,ebands%eig,enemax,enemin,ebands%fermie,dtset%mband,&
&   dtset%nband,nene,nkpt,nsppol,dtset%occopt,prtdos,&
&   dtset%tphysel,dtset%tsmear,unitdos)
 else if (dtset%prtdos == 3) then
   do iat=1,natsph+natsph_extra
     call dos_hdr_write(buffer,deltaene,dtset%dosdeltae,ebands%eig,enemax,enemin,ebands%fermie,dtset%mband,&
&     dtset%nband,nene,nkpt,nsppol,dtset%occopt,prtdos,&
&     dtset%tphysel,dtset%tsmear,unitdos_arr(iat))
   end do
 end if

 ABI_ALLOCATE(partial_dos,(nene,ndosfraction))
 ABI_ALLOCATE(integ_dos,(nene,ndosfraction))
 ABI_ALLOCATE(total_dos,(nene,ndosfraction))
 ABI_ALLOCATE(total_integ_dos,(nene,ndosfraction))
 ABI_ALLOCATE(tweight,(nene, nkpt))
 ABI_ALLOCATE(dtweightde,(nene, nkpt))

 if (paw_dos_flag==1) then
   ABI_ALLOCATE(total_dos_paw1,(nene,ndosfraction))
   ABI_ALLOCATE(total_dos_pawt1,(nene,ndosfraction))
 end if
 if (prtdosm>=1) then
   ABI_ALLOCATE(partial_dos_m,(nene,ndosfraction*mbesslang))
   ABI_ALLOCATE(integ_dos_m,(nene,ndosfraction*mbesslang))
   ABI_ALLOCATE(total_dos_m,(nene,ndosfraction*mbesslang))
   ABI_ALLOCATE(total_integ_dos_m,(nene,ndosfraction*mbesslang))
 end if

!Get maximum occupation value (2 or 1)
 max_occ = one; if (dtset%nspinor == 1 .and. nsppol == 1) max_occ = two

!-------------------------------------------------------------------
!For each spin polarisation and band, interpolate band over kpoints
!calculate integration weights and DOS contib from
!-------------------------------------------------------------------

 ! Workspace arrays.
 ABI_MALLOC(work_ndos, (nkpt, ndosfraction))
 ABI_MALLOC(work_ndosmbessl, (nkpt, ndosfraction*mbesslang))
 ABI_ALLOCATE(tmp_eigen,(nkpt))

 do isppol=1,nsppol
   total_dos = zero; total_integ_dos = zero
   if (prtdosm>=1) total_dos_m = zero
   if (paw_dos_flag==1) then
     total_dos_paw1(:,:)=zero;total_dos_pawt1(:,:)=zero
   end if

   if (nsppol==2) then
     if(isppol==1) write(message,'(a,16x,a)')  '#','Spin-up DOS'
     if(isppol==2) write(message,'(2a,16x,a)')  ch10,'#','Spin-dn DOS'
     ! NB: dtset%prtdos == 5 should not happen for nsppol==2
     if (dtset%prtdos == 2 .or. dtset%prtdos == 5) then
       write(unitdos, "(a)")trim(message)
     else if (dtset%prtdos == 3) then
       do iat=1,natsph+natsph_extra
         write(unitdos_arr(iat), "(a)")trim(message)
       end do
     end if
   end if

   do iband=1,dtset%nband(1)
     ! For each band get its contribution
     tmp_eigen(:) = ebands%eig(iband, :, isppol)

     ! calculate general integration weights at each irred kpoint as in Blochl et al PRB 49 16223
     call tetra_blochl_weights(tetrahedra,tmp_eigen,enemin,enemax,max_occ,nene,nkpt,&
       bcorr0,tweight,dtweightde,xmpi_comm_self)

     ! calculate DOS and integrated DOS projected with the input dos_fractions
     if (paw_dos_flag==1) then
       work_ndos = dos%fractions_paw1(:,iband,isppol,:)
       call get_dos_1band(work_ndos,integ_dos,nene,nkpt,ndosfraction,partial_dos,tweight,dtweightde)

       total_dos_paw1 = total_dos_paw1 + partial_dos

       work_ndos = dos%fractions_pawt1(:,iband,isppol,:)
       call get_dos_1band(work_ndos,integ_dos,nene,nkpt,ndosfraction,partial_dos,tweight,dtweightde)

       total_dos_pawt1 = total_dos_pawt1 + partial_dos
     end if

     work_ndos = dos%fractions(:,iband,isppol,:)
     call get_dos_1band(work_ndos,integ_dos,nene,nkpt,ndosfraction,partial_dos,tweight,dtweightde)

     ! Add to total dos
     total_dos = total_dos + partial_dos
     total_integ_dos = total_integ_dos + integ_dos

     if (prtdosm>=1) then
       work_ndosmbessl = dos%fractions_m(:,iband,isppol,:)
       call get_dos_1band_m(work_ndosmbessl,integ_dos_m,nene,nkpt,ndosfraction*mbesslang,&
&       partial_dos_m,tweight,dtweightde)

        total_dos_m = total_dos_m + partial_dos_m
        total_integ_dos_m = total_integ_dos_m + integ_dos_m
     end if

   end do ! iband
   bigDOS=(maxval(total_dos)>999._dp)

   !call wrtout(std_out,'about to write to the DOS file ',"COLL")
   ! header lines depend on the type of DOS (projected etc...) which is output
   enex=enemin
   integral_DOS=zero

   if(prtdos==2)then
     write(unitdos, '(a)' )'# energy(Ha)     DOS  integrated DOS'

   else if(prtdos==3)then
     if (paw_dos_flag/=1.or.dtset%pawprtdos==2) then
       do iat=1,natsph
         write(unitdos_arr(iat), '(3a,i5,a,i5,a,a,es16.6,3a)' ) &
&         '# Local DOS (columns 2-6) and integrated local DOS (columns 7-11),',ch10,&
&         '# for atom number iat=',iat,'  iatom=',dtset%iatsph(iat),ch10,&
&         '# inside sphere of radius ratsph=',dtset%ratsph(dtset%typat(dtset%iatsph(iat))),' Bohr.',ch10,"#"

         if (dtset%usepaw==1.and.dtset%pawprtdos==2) then
           write(unitdos_arr(iat), '(3a)' ) &
&           '# PAW: note that only all-electron on-site part has been used to compute DOS !',ch10,"#"
         end if
         if (bigDOS) then
           write(message, '(a,a)' ) &
&           '# energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&           '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
         else
           write(message, '(a,a)' ) &
&           '# energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&           '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
         end if
         if (prtdosm>=1) then
           write(message, '(7a)' ) trim(message),'          ',&
&           '  lm=0 0',&
&           '  lm=1-1  lm=1 0  lm=1 1',&
&           '  lm=2-2  lm=2-1  lm=2 0  lm=2 1  lm=2 2',&
&           '  lm=3-3  lm=3-2  lm=3-1  lm=3 0  lm=3 1  lm=3 2  lm=3 3',&
&           '  lm=4-4  lm=4-3  lm=4-2  lm=4-1  lm=4 0  lm=4 1  lm=4 2  lm=4 3  lm=4 4'
         end if
         write(unitdos_arr(iat), "(a)")trim(message)
       end do
     else
       do iat=1,natsph
         write(unitdos_arr(iat), '(9a,i5,a,i5,a,a,es16.6,3a)' ) &
&         '# Local DOS (columns 2-6),',ch10,&
&         '#  plane-waves contrib. to DOS (columns 7-11),',ch10,&
&         '#  AE on-site  contrib. to DOS (columns 12-16),',ch10,&
&         '# -PS on-site  contrib. to DOS (columns 17-21),',ch10,&
&         '# for atom number iat=',iat,'  iatom=',dtset%iatsph(iat),ch10,&
&         '# inside sphere of radius ratsph=',dtset%ratsph(dtset%typat(dtset%iatsph(iat))),' Bohr.',ch10,"#"
         if (bigDOS) then
           write(message, '(4a)' ) &
&           '#energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&           '       (PW)  l=0       l=1       l=2       l=3       l=4',&
&           '      (Phi)  l=0       l=1       l=2       l=3       l=4',&
&           '     (tPhi)  l=0       l=1       l=2       l=3       l=4'
         else
           write(message, '(4a)' ) &
&           '#energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&           '       (PW) l=0      l=1      l=2      l=3      l=4',&
&           '      (Phi) l=0      l=1      l=2      l=3      l=4',&
&           '     (tPhi) l=0      l=1      l=2      l=3      l=4'
         end if
         write(unitdos_arr(iat), "(a)")trim(message)
       end do
     end if
     do iat=1,natsph_extra
       write(unitdos_arr(natsph+iat), '(3a,i5,2a,es16.6,3a)' ) &
&       '# Local DOS (columns 2-6) and integrated local DOS (columns 7-11),',ch10,&
&       '# for non-atomic sphere number iat=',iat,ch10,&
&       '# of radius ratsph=',dtset%ratsph_extra,' Bohr.',ch10,"#"
       if (bigDOS) then
         write(message, '(a,a)' ) &
&         '# energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&         '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
       else
         write(message, '(a,a)' ) &
&         '# energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&         '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
       end if
       if (prtdosm>=1) then
         write(message, '(7a)' ) trim(message),'          ',&
&         '  lm=0 0',&
&         '  lm=1-1  lm=1 0  lm=1 1',&
&         '  lm=2-2  lm=2-1  lm=2 0  lm=2 1  lm=2 2',&
&         '  lm=3-3  lm=3-2  lm=3-1  lm=3 0  lm=3 1  lm=3 2  lm=3 3',&
&         '  lm=4-4  lm=4-3  lm=4-2  lm=4-1  lm=4 0  lm=4 1  lm=4 2  lm=4 3  lm=4 4'
       end if
       write(unitdos_arr(natsph+iat), "(a)")trim(message)
     end do
   else if(prtdos==5)then
     write(unitdos, '(a)' )&
&      '# energy(Ha)     DOS up,up  up,dn  dn,up  dn,dn  sigma_x sigma_y sigma_z  and integrated DOS components'
   end if ! prtdos value
   !  finished with header printing

!  Write the DOS value in the DOS file
!  Print the data for this energy. Note the upper limit (dos_max), to be consistent with the format.
!  The use of "E" format is not adequate, for portability of the self-testing procedure.
   if (prtdos==2) then
     ! E, DOS, IDOS
     do iene=1,nene
       write(unitdos, '(f11.5,1x,2(f10.4,1x))') &
&        enex, min(total_dos(iene,:), dos_max), total_integ_dos(iene,:)
       enex=enex+deltaene
     end do
   else if (prtdos==3) then
     ! E, DOS(L=1,LMAX), IDOS(L=1,LMAX)
     ! Here we assume mpsang = 5 in the format.
     if (paw_dos_flag/=1.or.dtset%pawprtdos==2) then
       frmt = '(f11.5,1x,5(f9.4,1x),10x,5(f8.2,1x),10x,25(f8.2,1x))'
       if (bigDOS) frmt = '(f11.5,1x,5(f10.4,1x),10x,5(f8.2,1x),10x,25(f8.2,1x))'
       ! for extra atoms in vacuum need more precision
       frmt_extra = '(f11.5,1x,5(f20.16,1x),10x,5(f20.16,1x),10x,25(f20.16,1x))'

       do iene=1,nene
         do iat=1,natsph
           if (prtdosm==0) then
             i1 = (iat-1)*mbesslang+1; i2 = iat*mbesslang
             write(unitdos_arr(iat), fmt=frmt) enex, &
&             min(total_dos(iene, i1:i2), dos_max), total_integ_dos(iene, i1:i2)
           else
             i1 = (iat-1)*mbesslang+1; i2 = iat*mbesslang
             write(unitdos_arr(iat), fmt=frmt) enex, &
&             min(total_dos(iene, i1:i2), dos_max),&
&             total_integ_dos(iene, i1:i2),&
&             min(total_dos_m(iene,(iat-1)*mbesslang**2+1:iat*mbesslang**2), dos_max)
           end if
         end do

         do iat=natsph+1,natsph+natsph_extra
           if (prtdosm==0) then
             i1 = (iat-1)*mbesslang+1; i2 = iat*mbesslang
             write(unitdos_arr(iat), fmt=frmt_extra) enex, &
&             total_dos(iene, i1:i2), &
&             total_integ_dos(iene, i1:i2)
           else
             i1 = (iat-1)*mbesslang+1; i2 = iat*mbesslang
             write(unitdos_arr(iat), fmt=frmt_extra) enex, &
&             total_dos(iene, i1:i2),&
&             total_integ_dos(iene, i1:i2),&
&             total_dos_m(iene,(iat-1)*mbesslang**2+1:iat*mbesslang**2)
           end if
         end do

         enex=enex+deltaene
       end do

     else
       frmt = '(f11.5,1x,5(f9.4,1x),3(6x,5f9.4))'
       if (bigDOS) frmt = '(f11.5,1x,5(f10.4,1x),3(6x,5f10.4))'
       ! for extra atom spheres in vacuum need more precision
       frmt_extra = '(f11.5,1x,5(f20.16,1x),3(6x,5f20.16))'

       do iene=1,nene
         do iat=1,natsph
           i1 = iat*5-4; i2 = iat*5
           write(unitdos_arr(iat), fmt=frmt) enex, &
&           min(total_dos(iene, i1:i2), dos_max),&
&           min(total_dos(iene,i1:i2) - total_dos_paw1(iene,i1:i2) + total_dos_pawt1(iene,i1:i2), dos_max),&
&           min(total_dos_paw1(iene,i1:i2), dos_max),&
&           min(total_dos_pawt1(iene,i1:i2), dos_max)
         end do
         do iat=natsph+1,natsph+natsph_extra
           i1 = iat*5-4; i2 = iat*5
           write(unitdos_arr(iat), fmt=frmt_extra) enex, &
&           min(total_dos(iene,i1:i2), dos_max),&
&           min(total_dos(iene,i1:i2) - total_dos_paw1(iene,i1:i2) + total_dos_pawt1(iene,i1:i2), dos_max),&
&           min(total_dos_paw1(iene,i1:i2), dos_max),&
&           min(total_dos_pawt1(iene,i1:i2), dos_max)
         end do
         enex=enex+deltaene
       end do
     end if

   else if(prtdos==5)then
     ! E, SPIN-DOS
     frmt = '(f11.5,1x,7(f9.4,1x),10x,7(f8.2,1x))'
     if (bigDOS) frmt = '(f11.5,1x,7(f10.4,1x),10x,7(f8.2,1x))'
     do iene=1,nene
       write(unitdos, fmt=frmt) enex, min(total_dos(iene,1:7), dos_max), total_integ_dos(iene,1:7)
       enex=enex+deltaene
     end do
   end if

   integral_DOS=sum(total_integ_dos(nene,:))
   write(message, '(a,es16.8)' ) ' tetrahedron : integrate to',integral_DOS
   call wrtout(std_out,message,'COLL')
 end do ! isppol

 ABI_DEALLOCATE(tmp_eigen)
 ABI_FREE(work_ndos)
 ABI_FREE(work_ndosmbessl)

 if(prtdos==2 .or. prtdos==5) then
   close(unitdos)
 else if(prtdos==3) then
   do iat=1,natsph+natsph_extra
     close(unitdos_arr(iat))
   end do
   ABI_DEALLOCATE(unitdos_arr)
 end if

 ABI_DEALLOCATE(indkpt)
 ABI_DEALLOCATE(kpt_fullbz)
 ABI_DEALLOCATE(partial_dos)
 ABI_DEALLOCATE(integ_dos)
 ABI_DEALLOCATE(total_dos)
 ABI_DEALLOCATE(total_integ_dos)
 ABI_DEALLOCATE(tweight)
 ABI_DEALLOCATE(dtweightde)

 if (prtdosm>=1)  then
   ABI_DEALLOCATE(partial_dos_m)
   ABI_DEALLOCATE(integ_dos_m)
   ABI_DEALLOCATE(total_dos_m)
   ABI_DEALLOCATE(total_integ_dos_m)
 end if

 if (paw_dos_flag==1)  then
   ABI_DEALLOCATE(total_dos_paw1)
   ABI_DEALLOCATE(total_dos_pawt1)
 end if

 call destroy_tetra(tetrahedra)

 call cwtime(cpu,wall,gflops,"stop")
 write(message,'(2(a,f8.2),a)')" tetrahedron: cpu_time: ",cpu,"[s], walltime: ",wall," [s]"
 call wrtout(std_out,message,"PERS")

end subroutine tetrahedron
!!***

!----------------------------------------------------------------------

!!****f* m_epjdos/gaus_dos
!! NAME
!! gaus_dos
!!
!! FUNCTION
!! Calculate DOS by gauss method
!!
!! INPUTS
!!  dos_fractions= projections of wavefunctions on each angular momentum Ylm
!!     which is the weight going into the DOS for an l-decomposed dos
!!  dos_fractions_paw1= contribution to dos fractions from the PAW partial waves (phi)
!!  dos_fractions_pawt1= contribution to dos fractions from the PAW pseudo partial waves (phi_tild)
!!  dtset     structured datatype, from which one uses :
!!   kpt(3,nkpt)  =irreducible kpoints
!!   kptrlatt(3,3)=lattice vectors for full kpoint grid
!!   nband        =number of electronic bands for each kpoint
!!   nkpt         =number of irreducible kpoints
!!   nshiftk      =number of kpoint grid shifts
!!   nsym         =number of symmetries
!!   pawprtdos    =option to output the individual contributions to the partial DOS (0, 1 or 2)
!!   shiftk(3,nshiftk)=kpoint shifts
!!   symrel(3,3,nsym)=symmetry matrices in real space
!!   typat(ntypat)=types of atoms
!!   usepaw       =option for PAW
!!  eigen(mband*nkpt*nsppol)=eigenvalues at irred kpoints
!!  fermie=Fermi energy
!!  fildata=name of the DOS output file
!!  mbesslang=maximum angular momentum for Bessel function expansion
!!  prtdosm=option for the m-contributions to the partial DOS
!!  ndosfraction= number of types of DOS we are calculating, e.g. the number
!!    of l channels. Could be much more general, for other types of partial DOS
!!  paw_dos_flag= option for partial dos in PAW
!!
!! OUTPUT
!!  (no explicit output)
!!
!! NOTES
!!  This routine should be called by master only
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      dos_hdr_write,int2char4,splfit,wrtout
!!
!! SOURCE

subroutine gaus_dos(dos, dtset, ebands, eigen,fildata)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gaus_dos'
 use interfaces_14_hidewrite
 use interfaces_61_occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: fildata
 type(dataset_type),intent(in) :: dtset
 type(epjdos_t),intent(in) :: dos
 type(ebands_t),intent(in) :: ebands
!arrays
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: nptsdiv2=6000
 integer :: bantot,ii,iat,iband,iene,ifract,ikpt,index,isppol,natsph,nkpt,nsppol,nene,prtdos
 integer :: mbesslang,ndosfraction
 real(dp) :: buffer,deltaene,enemax,enemin,integral_DOS,max_occ,enex,tsmearinv,tsmear,limit_occ,tratio
 real(dp) :: dblsmr,dsqrpi,increm,xx
 real(dp) :: cpu,wall,gflops
 logical :: bigDOS
 character(len=10) :: tag
 character(len=500) :: frmt,message
 character(len=fnlen) :: tmpfil
!scalars
 integer,allocatable :: unitdos_arr(:)
 real(dp) :: xgrid(-nptsdiv2:nptsdiv2),smdfun(-nptsdiv2:nptsdiv2,2)
 real(dp),allocatable :: vals(:),arg(:),derfun(:),total_dos(:,:)
 real(dp),allocatable :: total_integ_dos(:,:),total_dos_paw1(:,:),total_dos_pawt1(:,:)
 !real(dp),allocatable :: integ_dos_m(:,:,:),partial_dos_m(:,:,:),partial_dos(:,:,:),
 real(dp),allocatable :: total_dos_m(:,:) !,total_integ_dos_m(:,:)

! *********************************************************************

 natsph=dtset%natsph
 nsppol=dtset%nsppol
 nkpt=dtset%nkpt
 prtdos=dtset%prtdos
 deltaene=dtset%dosdeltae
 tsmear=dtset%tsmear
 tsmearinv=one/tsmear
 mbesslang = dos%mbesslang
 ndosfraction = dos%ndosfraction

 call cwtime(cpu, wall, gflops, "start")

 ! Open the DOS file
 ABI_ALLOCATE(unitdos_arr,(natsph))
 do iat=1,natsph
   call int2char4(dtset%iatsph(iat),tag)
   ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
   tmpfil = trim(fildata)//'_AT'//trim(tag)
   if (open_file(tmpfil, message, newunit=unitdos_arr(iat), status='unknown',form='formatted') /= 0) then
     MSG_ERROR(message)
   end if
 end do

 ! Write the header of the DOS file, and determine the energy range and spacing
 buffer=0.01_dp ! Size of the buffer around the min and max ranges

 do iat=1,natsph
   call dos_hdr_write(buffer,deltaene,dtset%dosdeltae,eigen,enemax,enemin,ebands%fermie,dtset%mband,&
&   dtset%nband,nene,dtset%nkpt,dtset%nsppol,dtset%occopt,prtdos,&
&   dtset%tphysel,dtset%tsmear,unitdos_arr(iat))
 end do

 bantot=sum(dtset%nband(:))
 ABI_ALLOCATE(arg,(bantot))
 ABI_ALLOCATE(derfun,(bantot))
 ABI_ALLOCATE(vals,(bantot))

 !write(std_out,*) ' ndosfraction,dtset%mband,dtset%nkpt,nene'ndosfraction,dtset%mband,dtset%nkpt,nene
 !ABI_ALLOCATE(partial_dos,(nene,ndosfraction,dtset%mband))
 ABI_ALLOCATE(total_dos,(nene,ndosfraction))
 ABI_ALLOCATE(total_integ_dos,(nene,ndosfraction))

 if (dos%paw_dos_flag==1) then
   ABI_ALLOCATE(total_dos_paw1,(nene,ndosfraction))
   ABI_ALLOCATE(total_dos_pawt1,(nene,ndosfraction))
   total_dos_paw1=zero
   total_dos_pawt1=zero
 end if
 if (dos%prtdosm>=1) then
 !  ABI_ALLOCATE(partial_dos_m,(nene,ndosfraction*mbesslang,dtset%mband))
 !  ABI_ALLOCATE(integ_dos_m,(nene,ndosfraction*mbesslang,dtset%mband))
    ! TODO: This is not used
    ABI_ALLOCATE(total_dos_m,(nene,ndosfraction*mbesslang))
 !  ABI_ALLOCATE(total_integ_dos_m,(nene,ndosfraction*mbesslang))
 end if

!Get maximum occupation value (2 or 1)
 max_occ = one
 if (dtset%nspinor == 1 .and. dtset%nsppol == 1) max_occ = two

!Define xgrid and gaussian on smdfun following getnel
 limit_occ=6.0_dp; dblsmr=zero
 if (abs(dtset%tphysel)>tol12 .and. abs(dtset%tsmear)>tol12) dblsmr = one

 if(dtset%occopt==3)limit_occ=30.0_dp
 if(dblsmr /= zero) then
   tratio = dtset%tsmear / dtset%tphysel
   limit_occ=30.0_dp + 6.0_dp*tratio
 end if

 increm=limit_occ/nptsdiv2
 do ii=-nptsdiv2,nptsdiv2
   xgrid(ii)=ii*increm
 end do
 dsqrpi=one/sqrt(pi)
 do ii=0,nptsdiv2
   xx=xgrid(ii)
   smdfun( ii,1)=dsqrpi*exp(-xx**2)
   smdfun(-ii,1)=smdfun(ii,1)
 end do

!calculate DOS and integrated DOS projected with the input dos_fractions
 total_dos=zero
 enex=enemin
 do iene=1,nene
   arg(:)=(enex-eigen(1:bantot))*tsmearinv
   call splfit(xgrid,derfun,smdfun,0,arg,vals,(2*nptsdiv2+1),bantot)
   index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       do iband=1,dtset%nband(ikpt)
         index=index+1
         do ifract=1,ndosfraction
           if (dos%paw_dos_flag==1) then
             total_dos_paw1(iene,ifract)=total_dos_paw1(iene,ifract)+&
&             dos%fractions_paw1(ikpt,iband,isppol,ifract)*dtset%wtk(ikpt) * max_occ*vals(index)*tsmearinv
             total_dos_pawt1(iene,ifract)=total_dos_pawt1(iene,ifract)+&
&             dos%fractions_pawt1(ikpt,iband,isppol,ifract)*dtset%wtk(ikpt) * max_occ*vals(index)*tsmearinv
           end if
           total_dos(iene,ifract) = total_dos(iene,ifract) + &
&           dos%fractions(ikpt,iband,isppol,ifract)*dtset%wtk(ikpt)*max_occ * vals(index)*tsmearinv
         end do
       end do ! iband
     end do ! ikpt
   end do ! isppol
   enex=enex+deltaene
 end do   ! iene

 do ifract=1,ndosfraction
   call simpson_int(nene,deltaene,total_dos(:,ifract), total_integ_dos(:, ifract))
 end do

!Write the DOS value in the DOS file
 do isppol=1,nsppol
   enex=enemin
   integral_DOS=zero

   bigDOS=(maxval(total_dos)>999._dp)
   if (dos%paw_dos_flag/=1.or.dtset%pawprtdos==2) then
     do iat=1,natsph
       write(unitdos_arr(iat), '(3a,i5,a,i5,a,a,es16.6,3a)' ) &
&       '# Local DOS (columns 3-7) and integrated local DOS (columns 8-12),',ch10,&
&       '# for atom number iat=',iat,'  iatom=',dtset%iatsph(iat),ch10,&
&       '# inside sphere of radius ratsph=',dtset%ratsph(dtset%typat(dtset%iatsph(iat))),' Bohr.',ch10,"#"
       if (dtset%usepaw==1.and.dtset%pawprtdos==2) then
         write(unitdos_arr(iat), '(3a)' ) &
&         '# PAW: note that only all-electron on-site part has been used to compute DOS !',ch10,"#"
       end if
       if (bigDOS) then
         write(message, '(a,a)' ) &
&         '# energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&         '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
       else
         write(message, '(a,a)' ) &
&         '# energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&         '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
       end if
       if (dos%prtdosm>=1) then
         write(unitdos_arr(iat), '(7a)' ) trim(message),'          ',&
&         '  lm=0 0',&
&         '  lm=1-1  lm=1 0  lm=1 1',&
&         '  lm=2-2  lm=2-1  lm=2 0  lm=2 1  lm=2 2',&
&         '  lm=3-3  lm=3-2  lm=3-1  lm=3 0  lm=3 1  lm=3 2  lm=3 3',&
&         '  lm=4-4  lm=4-3  lm=4-2  lm=4-1  lm=4 0  lm=4 1  lm=4 2  lm=4 3  lm=4 4'
       end if
     end do
   else
     do iat=1,natsph
       write(unitdos_arr(iat), '(9a,i5,a,i5,a,a,es16.6,3a)' ) &
&       '# Local DOS (columns 3-7),',ch10,&
&       '#  plane-waves contrib. to DOS (columns 8-12),',ch10,&
&       '#  AE on-site  contrib. to DOS (columns 13-17),',ch10,&
&       '# -PS on-site  contrib. to DOS (columns 18-22),',ch10,&
&       '# for atom number iat=',iat,'  iatom=',dtset%iatsph(iat),ch10,&
&       '# inside sphere of radius ratsph=',dtset%ratsph(dtset%typat(dtset%iatsph(iat))),' Bohr.',ch10,"#"
       if (bigDOS) then
         write(message, '(4a)' ) &
&         '#energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&         '       (PW)  l=0       l=1       l=2       l=3       l=4',&
&         '      (Phi)  l=0       l=1       l=2       l=3       l=4',&
&         '     (tPhi)  l=0       l=1       l=2       l=3       l=4'
       else
         write(message, '(4a)' ) &
&         '#energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&         '       (PW) l=0      l=1      l=2      l=3      l=4',&
&         '      (Phi) l=0      l=1      l=2      l=3      l=4',&
&         '     (tPhi) l=0      l=1      l=2      l=3      l=4'
       end if
       write(unitdos_arr(iat), "(a)")trim(message)
     end do
   end if

   if (dos%paw_dos_flag/=1.or.dtset%pawprtdos==2) then
     frmt='(f11.5,5f9.4 ,10x,5f8.2,10x,25f8.2)'
     if (bigDOS) frmt='(f11.5,5f10.4,10x,5f8.2,10x,25f8.2)'

     do iene=1,nene
       do iat=1,natsph
!        write(message, '(i6,5f9.4,10x,5f7.2))') iene-1, total_dos(iene,iat*5-2)
!        Note the upper limit, to be consistent with the format. The use of "E" format is not adequate,
!        for portability of the self-testing procedure.
         if (dos%prtdosm==0) then
           write(unitdos_arr(iat), fmt=frmt) enex, &
&           min(total_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang),9999.9999_dp), &
&           total_integ_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang)
         else
           write(unitdos_arr(iat), fmt=frmt) enex, &
&           min(total_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang),9999.9999_dp),&
&           total_integ_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang),&
&           min(total_dos_m(iene,(iat-1)*mbesslang**2+1:iat*mbesslang**2),9999.9999_dp)
         end if
       end do
       enex=enex+deltaene
     end do
   else
     frmt='(f11.5,5f9.4 ,3(6x,5f9.4 ))'
     if (bigDOS) frmt='(f11.5,5f10.4,3(6x,5f10.4))'

     do iene=1,nene
       do iat=1,natsph
         write(unitdos_arr(iat), fmt=frmt) enex, &
&         min(total_dos(iene,iat*5-4:iat*5),9999.9999_dp),&
&         min(total_dos(iene,iat*5-4:iat*5)-total_dos_paw1(iene,iat*5-4:iat*5)&
&         +total_dos_pawt1(iene,iat*5-4:iat*5),9999.9999_dp),&
&         min(total_dos_paw1(iene,iat*5-4:iat*5),9999.9999_dp),&
&         min(total_dos_pawt1(iene,iat*5-4:iat*5),9999.9999_dp)
       end do
       enex=enex+deltaene
     end do
   end if

   ! integral_DOS=integral_DOS+deltaene*sum(total_dos(iene,:))
   integral_DOS=sum(total_integ_dos(nene,:))
   !write(message, '(a,es16.8)' ) ' gaus_dos: integrate to',integral_DOS
   !call wrtout(std_out,message,'COLL')
 end do ! isppol

 do iat=1,natsph
   close(unitdos_arr(iat))
 end do

 call cwtime(cpu,wall,gflops,"stop")
 write(message,'(2(a,f8.2),a)')" gaus_dos: cpu_time: ",cpu,"[s], walltime: ",wall," [s]"
 call wrtout(std_out,message,"PERS")

 ABI_DEALLOCATE(unitdos_arr)
 ABI_FREE(arg)
 ABI_FREE(derfun)
 ABI_FREE(vals)
 !ABI_DEALLOCATE(partial_dos)
 ABI_DEALLOCATE(total_dos)
 ABI_DEALLOCATE(total_integ_dos)

 if (dos%paw_dos_flag==1)  then
   ABI_DEALLOCATE(total_dos_paw1)
   ABI_DEALLOCATE(total_dos_pawt1)
 end if
 if (dos%prtdosm>=1)  then
 !  ABI_DEALLOCATE(partial_dos_m)
 !  ABI_DEALLOCATE(integ_dos_m)
   ABI_DEALLOCATE(total_dos_m)
 !  ABI_DEALLOCATE(total_integ_dos_m)
 end if

end subroutine gaus_dos
!!***

!!****f* m_epjdos/get_dos_1band
!! NAME
!! get_dos_1band
!!
!! FUNCTION
!! calculate DOS from tetrahedron method for 1 band and 1 sppol
!!
!! INPUTS
!! dos_fractions=fractional DOS at each irred kpoint
!! nene=number of energies for DOS
!! nkpt=number of irreducible kpoints
!! ndosfraction=number of different fractional DOSs
!! tweight=sum of tetrahedron weights for each irred kpoint
!! dtweightde=energy derivative of tweight
!!
!! OUTPUT
!!  partial_dos(nene,ndosfraction)=partial DOS, for each different channel
!!  integ_dos(nene,ndosfraction)=integrated DOS, for each different channel
!!
!! PARENTS
!!      tetrahedron
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_dos_1band(dos_fractions,integ_dos,nene,nkpt,ndosfraction,&
&            partial_dos,tweight,dtweightde)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_dos_1band'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndosfraction,nene,nkpt
!arrays
 real(dp),intent(in) :: dos_fractions(nkpt,ndosfraction),dtweightde(nene,nkpt)
 real(dp),intent(in) :: tweight(nene,nkpt)
 real(dp),intent(out) :: integ_dos(nene,ndosfraction)
 real(dp),intent(out) :: partial_dos(nene,ndosfraction)

!Local variables-------------------------------
!scalars
 integer :: ieps,ifract,ikpt

! *********************************************************************

 partial_dos = zero; integ_dos = zero

 ! Calculate parameters of DOS at each point eps in [epsmin,epsmax]
 do ifract=1,ndosfraction
   do ikpt=1,nkpt
     do ieps=1,nene
       partial_dos(ieps,ifract) = partial_dos(ieps,ifract) + dtweightde(ieps,ikpt)*dos_fractions(ikpt,ifract)
       integ_dos(ieps,ifract) = integ_dos(ieps,ifract) + tweight(ieps,ikpt)*dos_fractions(ikpt,ifract)
     end do
   end do
 end do

end subroutine get_dos_1band
!!***

!!****f* m_epjdos/get_dos_1band_m
!! NAME
!! get_dos_1band_m
!!
!! FUNCTION
!! calculate DOS from tetrahedron method for 1 band and 1 sppol
!!
!! INPUTS
!! dos_fractions_m=fractional DOS at each irred kpoint
!! nene=number of energies for DOS
!! nkpt=number of irreducible kpoints
!! ndosfraction_m=number of different fractional DOSs
!! tweight=sum of tetrahedron weights for each irred kpoint
!! dtweightde=energy derivative of tweight
!!
!! OUTPUT
!!  partial_dos_m(nene,ndosfraction_m)=partial DOS, for each different channel
!!  integ_dos_m_m(nene,ndosfraction_m)=integrated DOS, for each different channel
!!
!! PARENTS
!!      tetrahedron
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_dos_1band_m(dos_fractions_m,integ_dos_m,nene,nkpt,ndosfraction_m,&
&                          partial_dos_m,tweight,dtweightde)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_dos_1band_m'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndosfraction_m,nene,nkpt
!arrays
 real(dp),intent(in) :: dos_fractions_m(nkpt,ndosfraction_m)
 real(dp),intent(in) :: dtweightde(nene,nkpt),tweight(nene,nkpt)
 real(dp),intent(out) :: integ_dos_m(nene,ndosfraction_m)
 real(dp),intent(out) :: partial_dos_m(nene,ndosfraction_m)

!Local variables-------------------------------
!scalars
 integer :: ieps,ifract,ikpt

! *********************************************************************

 partial_dos_m = zero; integ_dos_m = zero

 ! Calculate parameters of DOS at each point eps in [epsmin,epsmax]
 do ifract=1,ndosfraction_m
   do ikpt=1,nkpt
     do ieps=1,nene
       partial_dos_m(ieps,ifract) = partial_dos_m(ieps,ifract) + dtweightde(ieps,ikpt)*dos_fractions_m(ikpt,ifract)
       integ_dos_m(ieps,ifract) = integ_dos_m(ieps,ifract) + tweight(ieps,ikpt)*dos_fractions_m(ikpt,ifract)
     end do
   end do
 end do

end subroutine get_dos_1band_m
!!***

!!****f* m_epjdos/recip_ylm
!! NAME
!! recip_ylm
!!
!! FUNCTION
!! Project input wavefunctions in reciprocal space on to Ylm
!! (real or complex harmonics depending on rc_ylm).
!!
!! INPUTS
!!  bess_fit(mpw,nradintmax,ll) = Bessel functions for L, splined
!!   with arguments $2 \pi |k+G| \Delta r$, for all G vectors in sphere
!!   and all points on radial grid.
!!  cg_1band(2,npw_k)=wavefunction in recip space (note that nspinor is missing, see Notes).
!!  comm_pw=MPI communicator over plane waves (all npw-dependent data are distributed)
!!  istwfk= storage mode of cg_1band
!!  nradint(natsph)=number of points on radial real-space grid for a given atom.
!!  nradintmax=dimension of rint array.
!!  me_g0=1 if this processor has G=0, 0 otherwise
!!  mlang=maximum angular momentum in Bessel functions.
!!  mpw=Maximum number of planewaves. Used to dimension bess_fit
!!  natsph=number of atoms around which ang mom projection has to be done
!!  typat_extra(natsph)=Type of each atom. ntypat + 1 if empty sphere
!!  mlang_type(ntypat + natsph_extra)=Max L+1 for each atom type
!!  npw_k=number of plane waves for this kpt
!!  ph3d(2,npw_k,natsph)=3-dim structure factors, for each atom and plane wave.
!!  prtsphere= if 1, print a complete analysis of the angular momenta in atomic spheres
!!  rint(nradintmax) = points on radial real-space grid for integration
!!  rmax(natsph)=maximum radius for real space integration sphere
!!  rc_ylm= 1 for real spherical harmonics. 2 for complex spherical harmonics,
!!  ucvol=unit cell volume in bohr**3.
!!  ylm_k(npw_k,mlang**2)=real spherical harmonics for each G and LM.
!!  znucl_sph(natsph)=gives the nuclear number for each type of atom
!!
!! OUTPUT
!!  sum_1ll_1atom(mlang,natsph)= projected scalars for each atom and ang. mom.
!!  sum_1lm_1atom(mlang*mlang,natsph)= projected scalars for each atom and LM component.
!!
!! NOTES
!!  * ph3d atoms are ordered with natsph and must be provided by the caller in the correct order!
!!
!!  * spinor components are not treated here. This facilitates the implementation of spinor parallelism
!!    because the caller can easily call the routine inside a loop over spinors and then sum the
!!    different contributions outside the loop thus reducing the number of MPI calls.
!!
!! PARENTS
!!      m_cut3d,partial_dos_fractions
!!
!! CHILDREN
!!      atomdata_from_znucl,dotprod_g
!!
!! SOURCE

subroutine recip_ylm (bess_fit,cg_1band,comm_pw,istwfk,nradint,nradintmax,me_g0,mlang,&
&  mpw,natsph, typat_extra, mlang_type, npw_k,ph3d,prtsphere,rint,rmax,&
&  rc_ylm,sum_1ll_1atom,sum_1lm_1atom,ucvol,ylm_k,znucl_sph)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'recip_ylm'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm_pw,istwfk,me_g0,mlang,mpw,natsph,npw_k,nradintmax
 integer,intent(in) :: prtsphere,rc_ylm
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: nradint(natsph),typat_extra(natsph),mlang_type(:)
 real(dp),intent(in) :: bess_fit(mpw,nradintmax,mlang),cg_1band(2,npw_k)
 real(dp),intent(in) :: ph3d(2,npw_k,natsph),rint(nradintmax)
 real(dp),intent(in) :: rmax(natsph),ylm_k(npw_k,mlang*mlang)
 real(dp),intent(in) :: znucl_sph(natsph)
 real(dp),intent(out) :: sum_1ll_1atom(mlang,natsph)
 real(dp),intent(out) :: sum_1lm_1atom(mlang*mlang,natsph)

!Local variables-------------------------------
!scalars
 integer :: ilm,iat,ipw,ixint,ll,mm,il,jlm,ierr,lm_size,itypat
 real(dp),parameter :: invsqrt2=one/sqrt2
 real(dp) :: doti, dotr, sum_all, dr, fact
 type(atomdata_t) :: atom
 character(len=500) :: msg
!arrays
 integer :: ilang(mlang**2)
 real(dp) :: c1(2),c2(2)
 real(dp) :: sum_1atom(natsph),sum_1ll(mlang),sum_1lm(mlang**2)
 real(dp) :: func(nradintmax)
 real(dp) :: tmppsia(2,npw_k),tmppsim(2,npw_k),vect(2,npw_k)
 real(dp),allocatable :: values(:,:,:,:)

! *************************************************************************

 ! Workspace array (used to reduce the number of MPI communications)
 ! One could reduce a bit the memory requirement by using non-blocking operations ...
 ABI_STAT_MALLOC(values, (2, nradintmax, mlang**2, natsph), ierr)
 ABI_CHECK(ierr==0, "oom in values")
 values = zero

 sum_1lm_1atom = zero

 do ll=0,mlang-1
   do mm=-ll,ll
     ilm = (ll+1)**2-ll+mm
     ilang(ilm) = ll+1
   end do
 end do

 ! Big loop on all atoms
 do iat=1,natsph
   itypat = typat_extra(iat)
   lm_size = mlang_type(itypat) ** 2
   dr = rmax(iat) / (nradint(iat)-1)

   ! u(G) e^{i(k+G).Ra}
   ! Temporary array for part which depends only on iat
   do ipw=1,npw_k
     tmppsia(1,ipw) = cg_1band(1,ipw) * ph3d(1,ipw,iat) - cg_1band(2,ipw) * ph3d(2,ipw,iat)
     tmppsia(2,ipw) = cg_1band(1,ipw) * ph3d(2,ipw,iat) + cg_1band(2,ipw) * ph3d(1,ipw,iat)
   end do

   ! tmppsim = temporary arrays for part of psi which doesnt depend on ixint
   ! u(G) Y_LM^*(k+G) e^{i(k+G).Ra}
   ! Take into account the fact that ylm_k are REAL spherical harmonics, see initylmg.f
   ! For time-reversal states, detailed treatment show that only the real or imaginary
   ! part of tmppsia is needed here, depending on l being even or odd: only one of the coef is 1, the other 0
   do ilm=1,lm_size
     il = ilang(ilm)
     ll = ilang(ilm) - 1
     mm = ilm - (ll+1)**2 + ll

     select case (rc_ylm)
     case (1)
       ! to get PDOS for real spherical harmonics, multiply here by ylm_k instead of linear combination
       do ipw=1,npw_k
         tmppsim(1,ipw) = tmppsia(1,ipw) * ylm_k(ipw,ilm)
         tmppsim(2,ipw) = tmppsia(2,ipw) * ylm_k(ipw,ilm)
       end do

       ! Handle time-reversal
       if (istwfk /= 1) then
         if (mod(ll, 2) == 0) then
            tmppsim(2,:) = zero
         else
            tmppsim(1,:) = tmppsim(2,:)
            tmppsim(2,:) = zero
         end if
       end if

     case (2)
       ! to get PDOS for complex spherical harmonics, build linear combination of real ylm_k
       jlm = (ll+1)**2-ll-mm ! index of (l, -m)
       if (mm == 0) then
         vect(1,:) = ylm_k(1:npw_k,ilm)
         vect(2,:) = zero
       else if (mm > 0) then
          !vect(1,:) =  invsqrt2 * ylm_k(1:npw_k,ilm) * (-1)**mm
          !vect(2,:) = +invsqrt2 * ylm_k(1:npw_k,jlm) * (-1)**mm
          c1 = sy(ll, mm, mm)
          c2 = sy(ll,-mm, mm)
          vect(1,:) = c1(1) * ylm_k(1:npw_k,ilm) + c2(1) * ylm_k(1:npw_k,jlm)
          vect(2,:) = c1(2) * ylm_k(1:npw_k,ilm) + c2(2) * ylm_k(1:npw_k,jlm)

       else if (mm < 0) then
          !vect(1,:) =  invsqrt2 * ylm_k(1:npw_k,jlm) !* (-1)**mm
          !vect(2,:) = -invsqrt2 * ylm_k(1:npw_k,ilm) !* (-1)**mm
          c1 = sy(ll, mm,  mm)
          c2 = sy(ll,-mm,  mm)
          vect(1,:) = c1(1) * ylm_k(1:npw_k,ilm) + c2(1) * ylm_k(1:npw_k,jlm)
          vect(2,:) = c1(2) * ylm_k(1:npw_k,ilm) + c2(2) * ylm_k(1:npw_k,jlm)
       end if
       vect(2,:) = -vect(2,:)

       if (istwfk == 1) then
         do ipw=1,npw_k
           tmppsim(1, ipw) = tmppsia(1, ipw) * vect(1, ipw) - tmppsia(2, ipw) * vect(2, ipw)
           tmppsim(2, ipw) = tmppsia(1, ipw) * vect(2, ipw) + tmppsia(2, ipw) * vect(1, ipw)
         end do
       else
         ! Handle time-reversal
         if (mod(ll, 2) == 0) then
           do ipw=1,npw_k
             tmppsim(1, ipw) = tmppsia(1, ipw) * vect(1, ipw)
             tmppsim(2, ipw) = tmppsia(1, ipw) * vect(2, ipw)
           end do
         else
           do ipw=1,npw_k
             tmppsim(1, ipw) = tmppsia(2, ipw) * vect(1, ipw)
             tmppsim(2, ipw) = tmppsia(2, ipw) * vect(2, ipw)
           end do
         end if
       end if

     case default
       MSG_ERROR("Wrong value for rc_ylm")
     end select

     ! Compute integral $ \int_0^{rc} dr r**2 ||\sum_G u(G) Y_LM^*(k+G) e^{i(k+G).Ra} j_L(|k+G| r)||**2 $
     do ixint=1,nradint(iat)
       dotr = zero; doti = zero
       do ipw=1,npw_k
         dotr = dotr + bess_fit(ipw, ixint, il) * tmppsim(1, ipw)
         doti = doti + bess_fit(ipw, ixint, il) * tmppsim(2, ipw)
       end do
       if (istwfk /= 1) then
         dotr = two * dotr; doti = two * doti
         if (istwfk == 2 .and. me_g0 == 1) then
           dotr = dotr - bess_fit(1, ixint, il) * tmppsim(1, 1)
           doti = doti - bess_fit(1, ixint, il) * tmppsim(2, 1)
         end if
       end if

       ! Store results to reduce number of xmpi_sum calls if MPI
       values(1, ixint, ilm, iat) = dotr
       values(2, ixint, ilm, iat) = doti
     end do ! ixint

   end do ! ilm
 end do ! iat

 ! Collect results in comm_pw (data are distributed over plane waves)
 call xmpi_sum(values, comm_pw, ierr)

 ! Multiply by r**2 and take norm, integrate
 do iat=1,natsph
   itypat = typat_extra(iat)
   lm_size = mlang_type(itypat) ** 2
   do ilm=1,lm_size
     do ixint=1,nradint(iat)
       func(ixint) = rint(ixint)**2 * (values(1, ixint, ilm, iat)**2 + values(2, ixint, ilm, iat)**2)
     end do
     ! Here I should treat the case in which the last point /= rcut
     sum_1lm_1atom(ilm, iat) = simpson(dr, func(1:nradint(iat)))
   end do
 end do

 ! Normalize with unit cell volume and include 4pi term coming from Rayleigh expansion.
 fact = four_pi**2 / ucvol
 sum_1lm_1atom = fact * sum_1lm_1atom
 sum_1ll_1atom = zero
 do iat=1,natsph
   itypat = typat_extra(iat)
   lm_size = mlang_type(itypat) ** 2
   do ilm=1,lm_size
     il = ilang(ilm)
     sum_1ll_1atom(il, iat) = sum_1ll_1atom(il, iat) + sum_1lm_1atom(ilm, iat)
   end do
 end do

 ABI_FREE(values)

 ! Output
 if (prtsphere == 1) then
   sum_1ll = zero
   sum_1lm = zero
   sum_1atom = zero
   do iat=1,natsph
     sum_1atom(iat) = sum(sum_1lm_1atom(:,iat))
     sum_1ll(:)=sum_1ll(:)+sum_1ll_1atom(:,iat)
     sum_1lm(:)=sum_1lm(:)+sum_1lm_1atom(:,iat)
   end do
   sum_all = sum(sum_1atom)

   if (rc_ylm == 1) msg = " Angular analysis (real spherical harmonics)"
   if (rc_ylm == 2) msg =" Angular analysis (complex spherical harmonics)"
   call wrtout(std_out, msg)
   do iat=1,natsph
     call atomdata_from_znucl(atom, znucl_sph(iat))
     call wrtout(std_out, " ")
     write(msg,'(a,i3,a,a,a,f10.6)' )' Atom # ',iat, ' is  ',  atom%symbol,', in-sphere charge =',sum_1atom(iat)
     call wrtout(std_out, msg)
     do ll=0,mlang-1
       write(msg,'(a,i1,a,f9.6,a,9f6.3)' )&
&       ' l=',ll,', charge=',sum_1ll_1atom(ll+1,iat),&
&       ', m=-l,l splitting:',sum_1lm_1atom(1+ll**2:(ll+1)**2,iat)
       call wrtout(std_out, msg)
     end do ! ll
   end do ! iat
   write(msg,'(a,a)') ch10,' Sum of angular contributions for all atomic spheres '
   call wrtout(std_out, msg)
   do ll=0,mlang-1
     write(msg,'(a,i1,a,f9.6,a,f9.6)' )&
&     ' l=',ll,', charge =',sum_1ll(ll+1),' proportion =',sum_1ll(ll+1)/sum_all
     call wrtout(std_out, msg)
   end do
   write(msg,'(a,a,f10.6)' ) ch10,' Total over all atoms and l=0 to 4 :',sum_all
   call wrtout(std_out, msg)
   call wrtout(std_out, " ")
 end if

contains

 function sy(ll, mm, mp)
   use  m_paw_sphharm, only : ys
   ! Computes the matrix element <Slm|Ylm'>

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sy'
!End of the abilint section

   integer,intent(in) :: ll,mm, mp

   real(dp) :: sy(2)
   complex(dpc) :: ys_val

   ! Computes the matrix element <Yl'm'|Slm>
   call ys(ll,mp,ll,mm,ys_val)
   !call ys(ll,mm,ll,mp,ys_val)
   sy(1) = real(ys_val)
   sy(2) = -aimag(ys_val)

 end function sy

end subroutine recip_ylm
!!***

!!****f* m_epjdos/dens_in_sph
!! NAME
!! dens_in_sph
!!
!! FUNCTION
!!  Calculate integrated density in sphere around each atom
!!
!! INPUTS
!!  cg      = wavefunction coefficitents in recip space
!!  gmet    = metric in recip space
!!  istwfk  = storage mode for cg coefficients
!!  kg_k    = G vector indices
!!  natom   = number of atoms
!!  mpi_enreg=information about MPI parallelization
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  npw_k   = number of plane waves for this kpoint
!!  ph1d    = phase factors for different atoms for all G vectors
!!  rmax(natom) = max radius to integrate to (in bohr)
!!
!! OUTPUT
!!  cmax = integrated density for each atom for a rmax-radius sphere
!!
!! WARNING
!!  cg should not be modified by fourwf.
!!
!! PARENTS
!!      m_cut3d
!!
!! CHILDREN
!!      dotprod_v,fftpac,fourdp,fourwf,ph1d3d,sphereboundary,sphericaldens
!!      sqnorm_g
!!
!! SOURCE

subroutine dens_in_sph(cmax,cg,gmet,istwfk,kg_k,natom,ngfft,mpi_enreg,npw_k,&
&                       paral_kgb,ph1d,rmax,ucvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dens_in_sph'
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwfk,natom,npw_k,paral_kgb
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: kg_k(3,npw_k),ngfft(18)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: ph1d(2,(2*ngfft(1)+1+2*ngfft(2)+1+2*ngfft(3)+1)*natom)
 real(dp),intent(in) :: rmax(natom)
 real(dp),intent(inout) :: cg(2,npw_k)
 real(dp),intent(out) :: cmax(natom)

!Local variables -------------------------
!scalars
 integer,parameter :: tim_fourwf=0
 integer :: cplex,i1,i2,i3,iatom,id1,id2,id3,ifft,mgfft,n1,n2,n3,n4,n5,n6,nfft,nfftot
 real(dp) :: cmaxr,g1,g2,g3,norm,weight
!arrays
 integer :: ngfft_here(18)
 integer,allocatable :: garr(:,:),gbound(:,:)
 real(dp),allocatable :: denpot(:,:,:),fofgout(:,:),fofr(:,:,:,:),gnorm(:)
 real(dp),allocatable :: ph3d(:,:,:),phkxred(:,:),rhog(:,:),rhor(:)
 real(dp),allocatable :: sphrhog(:,:)

! *********************************************************************

 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
 n4=ngfft(4)
 n5=ngfft(5)
 n6=ngfft(6)
 nfftot = n1*n2*n3
 nfft=n1*n2*n3
 ngfft_here(:) = ngfft(:)
!fourwf doesnt work with other options for mode 0 (fft G -> r)
 ngfft_here(7)=111
 ngfft_here(8)=256
 mgfft=maxval(ngfft_here(1:3))

 call sqnorm_g(norm,istwfk,npw_k,cg,mpi_enreg%me_g0,mpi_enreg%comm_fft)

 if (abs(one-norm) > tol6) then
   write(std_out,'(a,f8.5)' ) ' dens_in_sph : this state is not normalized : norm=',norm
 end if

!-----------------------------------------------------------------
!inverse FFT of wavefunction to real space => density in real space
!-----------------------------------------------------------------
 ABI_ALLOCATE(gbound,(2*mgfft+8,2))
 call sphereboundary(gbound,istwfk,kg_k,mgfft,npw_k)

 weight = one
 cplex=1
 ABI_ALLOCATE(denpot,(cplex*n4,n5,n6))
 denpot(:,:,:)=zero
 ABI_ALLOCATE(fofgout,(2,npw_k))
 ABI_ALLOCATE(fofr,(2,n4,n5,n6))
 call fourwf(cplex,denpot,cg,fofgout,fofr,gbound,gbound, &
& istwfk,kg_k,kg_k,mgfft,mpi_enreg,1,ngfft_here,npw_k,&
& npw_k,n4,n5,n6,1,paral_kgb,tim_fourwf,weight,weight)
 ABI_DEALLOCATE(fofgout)
 ABI_DEALLOCATE(fofr)
 ABI_DEALLOCATE(gbound)

 norm = sum(denpot(:,:,:))/nfftot
 if (abs(one-norm) > tol6) then
   write(std_out,'(a,f8.5)') ' dens_in_sph : this state is not normalized in real space : norm=',norm
 end if

!-----------------------------------------------------------------
!FFT of new density: we obtain n(G) in rhog(1,:)
!-----------------------------------------------------------------

!Change the packing of the reciprocal space density
 ABI_ALLOCATE(rhor,(nfft))
 call fftpac(1,mpi_enreg,1,n1,n2,n3,n4,n5,n6,ngfft,rhor,denpot,1)

 ABI_ALLOCATE(rhog,(2,nfft))
 call fourdp(1,rhog,rhor,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)

 ABI_DEALLOCATE(rhor)
 ABI_DEALLOCATE(denpot)

 do ifft=1,nfft
   rhog(:,ifft) = rhog(:,ifft) / ucvol
 end do

!-----------------------------------------------------------------
!calculate norms of G vectors
!-----------------------------------------------------------------

 ABI_ALLOCATE(garr,(3,nfft))
 ABI_ALLOCATE(gnorm,(nfft))
 id3=ngfft(3)/2+2 ; id2=ngfft(2)/2+2 ; id1=ngfft(1)/2+2
 do i3=1,n3
   g3=i3-(i3/id3)*ngfft(3)-1
   do i2=1,n2
     g2=i2-(i2/id2)*ngfft(2)-1
     do i1=1,n1
       g1=i1-(i1/id1)*ngfft(1)-1
       ifft=i1+(i2-1)*n1+(i3-1)*n1*n2
       garr(1,ifft)=g1
       garr(2,ifft)=g2
       garr(3,ifft)=g3
       gnorm(ifft)=sqrt(gmet(1,1)*g1*g1 + &
&       two*gmet(2,1)*g2*g1 + &
&       two*gmet(3,1)*g3*g1 + &
&       gmet(2,2)*g2*g2 + &
&       gmet(3,2)*g3*g2 + &
&       gmet(3,3)*g3*g3)
     end do
   end do
 end do

!-----------------------------------------------------------------
!For each atom call sphericaldens to calculate
!n(G) * 1/|G|^3  *  int_0^2*\pi*r_{max}*|G| 4 \pi y^2 j_0 (y) dy
!for all G vectors put into array sphrhog
!scalar product of phase factors with spherically convoluted density
!-----------------------------------------------------------------

!largest mem occupation = nfft * (2(sphrog) +2*1(ph3d) +3(garr) +2(rhog) +1(gnorm)) = nfft * 10
 ABI_ALLOCATE(sphrhog,(2,nfft))
 ABI_ALLOCATE(phkxred,(2,natom))
 phkxred(1,:)=one
 phkxred(2,:)=zero
 ABI_ALLOCATE(ph3d,(2,nfft,1))

 do iatom=1,natom

   call sphericaldens(rhog,gnorm,nfft,rmax(iatom),sphrhog)
!  -----------------------------------------------------------------
!  Compute the phases for the whole set of fft vectors
!  -----------------------------------------------------------------

   call ph1d3d(iatom,iatom,garr,natom,natom,nfft,ngfft(1),ngfft(2),ngfft(3),&
&   phkxred,ph1d,ph3d)

!  For the phase factors, take the compex conjugate, before evaluating the scalar product
   do ifft=1,nfft
     ph3d(2,ifft,1)=-ph3d(2,ifft,1)
   end do
   cplex=2
   call dotprod_v(cplex,cmaxr,nfft,1,0,ph3d,sphrhog,mpi_enreg%comm_fft)
   cmax(iatom) = cmaxr
!  write(std_out,'(a,i4,a,es14.6,a,es12.6)' )' dens_in_sph : At ', iatom, ' has ',cmaxr, ' el.s in a sphere of rad ', rmax
 end do

 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(gnorm)
 ABI_DEALLOCATE(garr)
 ABI_DEALLOCATE(sphrhog)
 ABI_DEALLOCATE(ph3d)
 ABI_DEALLOCATE(phkxred)

end subroutine dens_in_sph
!!***

!!****f* m_epjdos/sphericaldens
!! NAME
!! sphericaldens
!!
!! FUNCTION
!! Compute the convolution of a function with
!! the unity constant function over a sphere of radius rmax .
!! The function is to be given in reciprocal space,
!! the resulting function is also given in reciprocal space.
!! The routine needs the norm of the reciprocal space vectors.
!!
!! The resulting function in reciprocal space can give the
!! integral of the density in any sphere of that radius, centered
!! on any point, by a simple scalar product.
!!
!! INPUTS
!!  fofg(2,nfft)=initial function, in reciprocal space
!!  gnorm(nfft)=norm of the reciprocal space vectors
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  rmax=radius of the sphere
!!
!! OUTPUT
!!  sphfofg(2,nfft)=convoluted function, in reciprocal space
!!
!! PARENTS
!!      dens_in_sph
!!
!! CHILDREN
!!
!! SOURCE

subroutine sphericaldens(fofg,gnorm,nfft,rmax,sphfofg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sphericaldens'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
 real(dp),intent(in) :: rmax
!arrays
 real(dp),intent(in) :: fofg(2,nfft),gnorm(nfft)
 real(dp),intent(out) :: sphfofg(2,nfft)

!Local variables-------------------------------
!scalars
 integer :: ifft
 real(dp) :: factor,int0yy,rmax_2pi,yy

! *************************************************************************

 rmax_2pi=two_pi*rmax
 factor=four_pi/(two_pi)**3

 do ifft=1,nfft
   if(abs(gnorm(ifft)) < tol12)then
     sphfofg(1,ifft)=fofg(1,ifft)*four_pi*third*rmax**3
     sphfofg(2,ifft)=fofg(2,ifft)*four_pi*third*rmax**3
   else
     yy=gnorm(ifft)*rmax_2pi
     int0yy=factor*(sin(yy)-yy*cos(yy))/(gnorm(ifft)**3)
     sphfofg(1,ifft)=fofg(1,ifft)*int0yy
     sphfofg(2,ifft)=fofg(2,ifft)*int0yy
   end if
 end do

end subroutine sphericaldens
!!***

!!****f* m_epjdos/prtfatbands
!! NAME
!! prtfatbands
!!
!! FUNCTION
!! Print dos_fractions_m in order to plot easily fatbands
!! if pawfatbnd=1  1 : fatbands are resolved in L.
!! if pawfatbnd=1  2 : fatbands are resolved in L and M.
!!
!! INPUTS
!!  dos_fractions_m(nkpt,mband,nsppol,ndosfraction*mbesslang*m_dos_flag)
!!               = m-resolved projected dos inside PAW sphere.
!!  dtset        = Input variables
!!  ebands<ebands_t>=Band structure data.
!!  pawfatbnd    = keyword for fatbands
!!  mbesslang    =maximum angular momentum for Bessel function expansion
!!  m_dos_flag   =option for the m-contributions to the partial DOS
!!  ndosfraction =natsph*mbesslang
!!
!! OUTPUT
!! (only writing)
!!
!! NOTES
!!  This routine should be called by master only
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      atomdata_from_znucl,int2char4,wrtout
!!
!! SOURCE

subroutine prtfatbands(dos,dtset,ebands,fildata,pawfatbnd,pawtab)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtfatbands'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: pawfatbnd
 type(epjdos_t),intent(in) :: dos
 type(ebands_t),intent(in) :: ebands
 type(dataset_type),intent(in) :: dtset
 character(len=fnlen),intent(in) :: fildata
!arrays
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer :: iall,il,iat,natsph,inbfatbands,iband,mband,ixfat,isppol,nkpt,lmax,ll,mm
 integer :: ikpt,nband_k,ndosfraction,mbesslang
 real(dp) :: xfatband,cpu,wall,gflops
 character(len=1) :: tag_l,tag_1m,tag_is
 character(len=2) :: tag_2m
 character(len=10) :: tag_il,tag_at,tag_grace
 character(len=500) :: message
 character(len=fnlen) :: tmpfil
 type(atomdata_t) :: atom
!arrays
 integer,allocatable :: unitfatbands_arr(:,:)
 real(dp),allocatable :: eigenvalues(:,:,:)
 character(len=2) :: symbol

!*************************************************************************

 DBG_ENTER("COLL")

 ndosfraction = dos%ndosfraction; mbesslang = dos%mbesslang

 if(dos%prtdosm.ne.0) then
   write(message,'(3a)')&
&   'm decomposed dos is activated',ch10, &
&   'Action: deactivate it with prtdosm=0 !'
   MSG_ERROR(message)
 end if

 if(dtset%nspinor==2) then
   MSG_WARNING("Fatbands are not yet available in the case nspinor==2!")
 end if

 ABI_CHECK(allocated(dos%fractions_m), "dos%fractions_m is not allocated!")

 natsph=dtset%natsph
 nkpt=dtset%nkpt
 mband=dtset%mband

 if(natsph>1000) then
   write(message,'(3a)')&
&   'Too big number of fat bands!',ch10, &
&   'Action: decrease natsph in input file !'
   MSG_ERROR(message)
 end if

!--------------  PRINTING IN LOG
 call cwtime(cpu, wall, gflops, "start")
 write(message,'(a,a,a,a,i5,a,a,1000i5)') ch10," ***** Print of fatbands activated ****** ",ch10,&
& "  Number of atom: natsph = ",natsph,ch10, &
& "  atoms  are             = ",(dtset%iatsph(iat),iat=1,natsph)
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 iall=0;inbfatbands=0

 if(pawfatbnd==1) then
   inbfatbands=mbesslang-1
   write(message,'(3a)')"  (fatbands are in eV and are given for each value of L)",ch10
 else if(pawfatbnd==2) then
   write(message,'(3a)')"  (fatbands are in eV and are given for each value of L and M)",ch10
   inbfatbands=(mbesslang-1)**2
 end if
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 write(message,'(a,e12.5,a,e12.5,a)') "  Fermi energy is ",ebands%fermie*Ha_eV," eV = ",ebands%fermie," Ha"
 call wrtout(std_out,message,'COLL')

!--------------  OPEN AND NAME FILES FOR FATBANDS
 ABI_ALLOCATE(unitfatbands_arr,(natsph*inbfatbands,dtset%nsppol))
 unitfatbands_arr = -3

 do iat=1,natsph
   lmax=(pawtab(dtset%typat(dtset%iatsph(iat)))%l_size-1)/2
   call int2char4(dtset%iatsph(iat),tag_at)
   ABI_CHECK((tag_at(1:1)/='#'),'Bug: string length too short!')
   call atomdata_from_znucl(atom,dtset%znucl(dtset%typat(dtset%iatsph(iat))))
   symbol = atom%symbol
   do il=1,inbfatbands
     iall=iall+1
     ll=int(sqrt(float(il-1)))  ! compute l
     if(ll.le.lmax) then  ! print only angular momentum included in the PAW data
       do isppol=1,dtset%nsppol
         write(tag_is,'(i1)')isppol
         if(pawfatbnd==1) then
           call int2char4(il-1,tag_il)
           ABI_CHECK((tag_il(1:1)/='#'),'Bug: string length too short!')
           tmpfil = trim(fildata)// &
&           '_at'//trim(tag_at)//'_'//trim(adjustl(symbol))//'_is'//tag_is//'_l'//trim(tag_il)
         else if (pawfatbnd==2) then
           write(tag_l,'(i1)') ll
           mm=il-(ll**2+ll+1)      ! compute m
           if(mm<0) write(tag_2m,'(i2)') mm
           if(mm>=0) write(tag_1m,'(i1)') mm
           if(mm<0) tmpfil = trim(fildata)// &
&           '_at'//trim(tag_at)//'_'//trim(adjustl(symbol))//'_is'//tag_is//'_l'//tag_l//'_m'//tag_2m
           if(mm>=0) tmpfil = trim(fildata)// &
&           '_at'//trim(tag_at)//'_'//trim(adjustl(symbol))//'_is'//tag_is//'_l'//tag_l//'_m+'//tag_1m
         end if
         !unitfatbands_arr(iall,isppol)=tmp_unit+100+iall-1+(natsph*inbfatbands)*(isppol-1)
         !open (unit=unitfatbands_arr(iall,isppol),file=trim(tmpfil),status='unknown',form='formatted')
         if (open_file(tmpfil, message, newunit=unitfatbands_arr(iall,isppol), status='unknown',form='formatted') /= 0) then
           MSG_ERROR(message)
         end if

         write(message,'(a,a,a,i4)') 'opened file : ', trim(tmpfil), ' unit', unitfatbands_arr(iall,isppol)
         call wrtout(std_out,message,'COLL')
         write(message,'(9a)') "# ",ch10,"# ABINIT package : FATBAND file ", ch10,&
&         "# It contains, for each band: the eigenvalues in eV (and the character of the band) as a function of the k-point",&
&         ch10,"# This file can be read with xmgrace (http://plasma-gate.weizmann.ac.il/Grace/)  ",ch10,"#  "
         write(unitfatbands_arr(iall,isppol), "(a)")trim(message)
         do iband=1,mband
           call int2char4(iband-1,tag_grace)
           ABI_CHECK((tag_grace(1:1)/='#'),'Bug: string length too short!')
           write(message,'(16a)') ch10,"@    s",trim(tag_grace)," line color 1",&
&           ch10,"@    s",trim(tag_grace)," errorbar color 2",&
&           ch10,"@    s",trim(tag_grace)," errorbar riser linewidth 5.0", &
&           ch10,"@    s",trim(tag_grace)," errorbar linestyle 0"
           write(unitfatbands_arr(iall,isppol), "(a)")trim(message)
         end do  !iband
         write(unitfatbands_arr(iall,isppol), '(a,a)') ch10,'@type xydy'
       end do   ! isppol
     end if ! ll=<lmax
   end do   ! il
 end do  ! iat

 if(iall.ne.(natsph*inbfatbands)) then
   MSG_ERROR("error1 ")
 end if



!--------------  WRITE FATBANDS IN FILES
 if (pawfatbnd>0) then
   ! Store eigenvalues with nkpt as first dimension for efficiency reasons
   ABI_ALLOCATE(eigenvalues,(nkpt,mband,dtset%nsppol))
   do isppol=1,dtset%nsppol
     do ikpt=1,nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
       do iband=1,mband
         eigenvalues(ikpt,iband,isppol)= ebands%eig(iband, ikpt, isppol) - ebands%fermie
       end do
     end do
   end do
   iall=0
   do iat=1,natsph
     lmax=(pawtab(dtset%typat(dtset%iatsph(iat)))%l_size-1)/2
     do il=1,inbfatbands
       iall=iall+1
       ll=int(sqrt(float(il-1)))
       if(ll.le.lmax) then
         do isppol=1,dtset%nsppol
           do iband=1,mband
             write(unitfatbands_arr(iall,isppol),'(a,a,i8)') ch10,"# BAND number :",iband
             do ikpt=1,nkpt
               if(pawfatbnd==1) then
                 xfatband=0.d0
                 do ixfat=(il-1)**2+1,il**2
                   xfatband=xfatband+dos%fractions_m(ikpt,iband,isppol,(iat-1)*mbesslang**2+ixfat)
                 end do ! ixfat
               else if (pawfatbnd==2) then
                 xfatband=dos%fractions_m(ikpt,iband,isppol,(iat-1)*mbesslang**2+il)
               end if
               write(unitfatbands_arr(iall,isppol),'(i5,e20.5,e20.5)')&
                 ikpt-1,eigenvalues(ikpt,iband,isppol)*Ha_eV,xfatband
             end do ! ikpt
           end do  !iband
           write(unitfatbands_arr(iall,isppol),'(a)') '&'
           !close(unitfatbands_arr(iall,isppol))
         end do  !isppol
       end if
     end do ! il
   end do ! iat
   ABI_DEALLOCATE(eigenvalues)
 end if

 do isppol=1,size(unitfatbands_arr, dim=2)
   do iat=1,size(unitfatbands_arr, dim=1)
     if (unitfatbands_arr(iat, isppol) /= -3) close (unitfatbands_arr(iat, isppol))
   end do
 end do

 ABI_DEALLOCATE(unitfatbands_arr)

 call cwtime(cpu,wall,gflops,"stop")
 write(message,'(2(a,f8.2),a)')" prtfatbands: cpu_time: ",cpu,"[s], walltime: ",wall," [s]"
 call wrtout(std_out,message,"PERS")

 DBG_EXIT("COLL")

end subroutine prtfatbands
!!***

!----------------------------------------------------------------------

!!****f* m_epjdos/fatbands_ncwrite
!! NAME
!! fatbands_ncwrite
!!
!! FUNCTION
!!
!! INPUTS
!!  crystal<crystal_t>=Object defining the unit cell and its symmetries.
!!  ncid=NC file handle.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE



subroutine fatbands_ncwrite(dos, crystal, ebands, hdr, dtset, psps, pawtab, ncid)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fatbands_ncwrite'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 type(epjdos_t),intent(in) :: dos
 type(crystal_t),intent(in) :: crystal
 type(ebands_t),intent(in) :: ebands
 type(hdr_type),intent(in) :: hdr
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)

#ifdef HAVE_NETCDF
!Local variables-------------------------------
!scalars
 integer :: itype,ncerr,fform
 real(dp) :: cpu,wall,gflops
 character(len=500) :: msg
!arrays
 integer :: lmax_type(crystal%ntypat)

!*************************************************************************

 ABI_CHECK(dtset%natsph > 0, "natsph <=0")
 call cwtime(cpu, wall, gflops, "start")

 fform = fform_from_ext("FATBANDS.nc")
 ABI_CHECK(fform /= 0, "Cannot find fform associated to FATBANDS.nc")

 ! Write header, crystal structure and band energies.
 NCF_CHECK(hdr_ncwrite(hdr, ncid, fform, nc_define=.True.))
 NCF_CHECK(crystal_ncwrite(crystal, ncid))
 NCF_CHECK(ebands_ncwrite(ebands, ncid))

 ! Add fatband-specific quantities
 ncerr = nctk_def_dims(ncid, [ &
   nctkdim_t("natsph", dtset%natsph), &
   nctkdim_t("ndosfraction", dos%ndosfraction)], defmode=.True.)
 NCF_CHECK(ncerr)

 if (dos%ndosfraction*dos%mbesslang > 0) then
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("mbesslang", dos%mbesslang), &
     nctkdim_t("dos_fractions_m_lastsize", dos%ndosfraction*dos%mbesslang)])
   NCF_CHECK(ncerr)
 end if
 if (dtset%natsph_extra /= 0) then
   NCF_CHECK(nctk_def_dims(ncid, [nctkdim_t("natsph_extra", dtset%natsph_extra)]))
 end if

 ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "prtdos", "pawprtdos", "prtdosm"])
 NCF_CHECK(ncerr)
 ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "ratsph_extra"])
 NCF_CHECK(ncerr)

 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t("lmax_type", "int", "number_of_atom_species"), &
   nctkarr_t("iatsph", "int", "natsph"), &
   nctkarr_t("ratsph", "dp", "number_of_atom_species"), &
   nctkarr_t("dos_fractions", "dp", "number_of_kpoints, max_number_of_states, number_of_spins, ndosfraction") &
 ])
 NCF_CHECK(ncerr)

 if (allocated(dos%fractions_m)) then
   ncerr = nctk_def_arrays(ncid, &
     nctkarr_t("dos_fractions_m", "dp", &
               "number_of_kpoints, max_number_of_states, number_of_spins, dos_fractions_m_lastsize"))
   NCF_CHECK(ncerr)
 end if

 if (allocated(dos%fractions_paw1)) then
   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t("dos_fractions_paw1", "dp", "number_of_kpoints, max_number_of_states, number_of_spins, ndosfraction"), &
     nctkarr_t("dos_fractions_pawt1", "dp", "number_of_kpoints, max_number_of_states, number_of_spins, ndosfraction") &
   ])
   NCF_CHECK(ncerr)
 end if

 if (dtset%natsph_extra /= 0) then
   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t("xredsph_extra", "dp", "number_of_reduced_dimensions, natsph_extra") &
   ])
   NCF_CHECK(ncerr)
 end if

 ! Write variables
 NCF_CHECK(nctk_set_datamode(ncid))

 ! scalars
 NCF_CHECK(nf90_put_var(ncid, vid("pawprtdos"), dtset%pawprtdos))
 NCF_CHECK(nf90_put_var(ncid, vid("prtdos"), dos%prtdos))
 NCF_CHECK(nf90_put_var(ncid, vid("prtdosm"), dos%prtdosm))

 ! arrays
 if (dtset%usepaw == 1) then
   lmax_type = (pawtab(:)%l_size - 1) / 2
 else
   do itype=1,crystal%ntypat
     lmax_type(itype) = maxval(psps%indlmn(1, :, itype))
   end do
 end if
 NCF_CHECK(nf90_put_var(ncid, vid("lmax_type"), lmax_type))
 NCF_CHECK(nf90_put_var(ncid, vid("dos_fractions"), dos%fractions))

 if (dos%prtdos == 3) then
   NCF_CHECK(nf90_put_var(ncid, vid("iatsph"), dtset%iatsph(1:dtset%natsph)))
   NCF_CHECK(nf90_put_var(ncid, vid("ratsph"), dtset%ratsph(1:dtset%ntypat)))
   NCF_CHECK(nf90_put_var(ncid, vid("ratsph_extra"), dtset%ratsph_extra))
   if (dtset%natsph_extra /= 0) then
     NCF_CHECK(nf90_put_var(ncid, vid("xredsph_extra"), dtset%xredsph_extra(:, 1:dtset%natsph_extra)))
   end if
 end if

 if (allocated(dos%fractions_m)) then
   NCF_CHECK(nf90_put_var(ncid, vid("dos_fractions_m"), dos%fractions_m))
 end if
 if (allocated(dos%fractions_paw1)) then
   NCF_CHECK(nf90_put_var(ncid, vid("dos_fractions_paw1"), dos%fractions_paw1))
   NCF_CHECK(nf90_put_var(ncid, vid("dos_fractions_pawt1"), dos%fractions_pawt1))
 end if

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2),a)')" fatbands_ncwrite: cpu_time: ",cpu,"[s], walltime: ",wall," [s]"
 call wrtout(std_out,msg,"PERS")
#endif

contains
 integer function vid(vname)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vid'
!End of the abilint section

   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine fatbands_ncwrite
!!***

end module m_epjdos
!!***
