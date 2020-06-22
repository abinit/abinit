!!****m* ABINIT/m_paw_atomorb
!! NAME
!!  m_paw_atomorb
!!
!! FUNCTION
!!  This module provides the definition of the atomorb_type used
!!  to store atomic orbitals on a radial mesh as well
!!  as methods to operate on it.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_atomorb

 use defs_basis
 use m_abicore
 use m_errors
 use m_splines

 use m_io_tools,      only : open_file
 use m_fstrings,      only : tolower
 use m_paw_lmn,       only : make_indlmn, make_indklmn, make_kln2ln
 use m_pawrad,        only : pawrad_type, pawrad_init, &
&                            pawrad_free, pawrad_print, pawrad_isame, pawrad_ifromr, simp_gen

 implicit none

 private
!!***

!!****t* m_paw_atomorb/atomorb_type
!! NAME
!!
!! FUNCTION
!!  Defines the atomorb_type datastructure type.
!!  It contains the atomic orbitals on a radial mesh for a given type of atom.
!!
!! NOTES
!!  * Should the radial mesh included in the object or not?
!!  * We might have different meshes, useful for deep states in heavy atoms!
!!  * Methods to be added: corekin
!!
!! SOURCE

 type, public :: atomorb_type

!scalars
  integer :: ixc
   ! Exchange and correlation functional used to generate the orbitals

  integer :: method
  ! 1 for restricted, compatible only with nsppol=1.
  ! 2 for spin unrestricted, compatible only with nsppol=2.

  integer :: nspden
   ! Number of spin-density components.

  integer :: nsppol
   ! Number of independent spin-components.
   ! FIXME: here a lot of quantities might depend on nsppol in the
   ! case of magnetic atoms!

  integer :: nspinor
   ! Number of spinorial components
   ! TODO this is a quite delicate issue, for S.O. one should use J = L+S instead of L!
   ! If we use scalar relativistic then...

  integer :: l_max
   ! Maximum value of angular momentum l+1

  integer :: l_size
  ! Maximum value of l+1 leading to non zero Gaunt coeffs
  ! l_size=2*l_max-1

  integer :: ln_size
  ! Number of (l,n) components.

  integer :: ln2_size
  ! ln2_size=ln_size*(ln_size+1)/2
  ! where ln_size is the number of (l,n) elements for core orbitals.

  integer :: lmn_size
  ! Number of (l,m,n) elements.

  integer :: lmn2_size
   ! lmn2_size=lmn_size*(lmn_size+1)/2
   ! where lmn_size is the number of (l,m,n) elements for core orbitals.

  integer :: mesh_size
  ! Dimension of the radial mesh.

  real(dp) :: rcore
  ! Radius of the sphere used to describe core electrons.
  ! It should be <= rpaw

  real(dp) :: zion
   ! zionpsp
   ! The ionic pseudo-charge, (giving raise to a long-range coulomb potential)

  ! TODO alchemy?
  !real(dp) :: ziontypat
   ! ziontypat
   !  For each type of atom (might be alchemy wrt psps), the ionic pseudo-charge
   ! (giving raise to a long-range coulomb potential)

  real(dp) :: znucl
   ! The atomic number of the atom.

  ! TODO alchemy?
  !real(dp) :: znucltypat
   ! znucltypat
   ! The atomic number of each type of atom (might be alchemy wrt psps)

  character(len=fnlen) :: fname
   ! The filename for temporary storage.

!arrays
  integer, allocatable :: indlmn(:,:)
  ! indlmn(6,lmn_size)
  ! Array giving l,m,n,lm,ln,spin for i=lmn.

  integer, allocatable :: indln(:,:)
  ! indln(2,ln_size)
  ! Array giving l and n for i=ln

  integer, allocatable :: indklmn(:,:)
   ! indklmn(8,lmn2_size)
   ! Array giving klm, kln, abs(il-jl), (il+jl), ilm and jlm, ilmn and jlmn for each klmn=(ilmn,jlmn)
   ! Note: ilmn=(il,im,in) and ilmn<=jlmn

  !integer, allocatable :: klm2lm TODO add
   !  klm2lm(6,lm2_size)=Table giving il, jl ,im, jm, ilm, jlm for each klm=(ilm,jlm)
   !  where ilm=(il,im) and ilm<=jlm. NB: klm2lm is an application and not a bijection.

  integer, allocatable :: klm_diag(:)
   ! klm_diag(lmn2_size)
   ! 1 il==jl and im==jm, 0 otherwise.

  integer, allocatable :: klmntomn(:,:)
   ! klmntomn(4,lmn2_size)
   ! Array giving im, jm ,in, and jn for each klmn=(ilmn,jlmn)
   ! Note: ilmn=(il,im,in) and ilmn<=jlmn
   ! NB: klmntomn is an application and not a bijection

  integer, allocatable :: kln2ln(:,:)
   ! kln2ln(6,ln2_size)
   ! Table giving il, jl ,in, jn, iln, jln for each kln=(iln,jln)
   ! where iln=(il,in) and iln<=jln. NB: kln2ln is an application and not a bijection

  integer, allocatable :: mode(:,:)
  ! mode(ln_size,nsppol)
  ! Flag defining how the orbital is treated.
  ! During the pseudopotential generation we can have: ORB_FROZEN or ORB_VALENCE
  ! For calculations in extended systems we can have:  ORB_FROZEN or ORB_RELAXED_CORE
  ! Namely different treatment depending of the degree of localization.
  ! For example the 1s in plutonium might be treated as ORB_FROZEN during
  ! a relaxed core calculation.
  ! TODO define function to test the type, much safer!

  real(dp), allocatable :: eig(:,:)
  ! eig(ln_size,nsppol)
  ! Eigenvalues for each ln channel and spin.

  real(dp), allocatable :: occ(:,:)
  ! occ(ln_size,nsppol)
  ! Occupation for each ln channel and spin.

  real(dp), allocatable :: phi(:,:,:)
  ! phi(mesh_size,ln_size,nsppol)
  ! Here we might have different meshes, useful for deep states in heavy atoms!

  ! this might be retrieved with  a method get_atomden
  !real(dp), allocatable :: density(:)
   ! density(mesh_size,nspden)
   ! Gives the core density of the atom for each spin channel
   ! Total charge in first dimension,up component in second one (if present)

 end type atomorb_type

! public procedures.
 public :: init_atomorb
 public :: destroy_atomorb
 public :: print_atomorb
 public :: get_overlap
!!***

!----------------------------------------------------------------------

 integer,public,parameter :: ORB_FROZEN       =1
 integer,public,parameter :: ORB_RELAXED_CORE =2
 integer,public,parameter :: ORB_VALENCE      =4


CONTAINS  !=========================================================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_atomorb/destroy_atomorb
!! NAME
!!  destroy_atomorb
!!
!! FUNCTION
!!  Free the dynamic memory allocated in a structure of type atomorb_type.
!!
!! SIDE EFFECTS
!!  Atm <type(atomorb_type)>=datastructure containing atomic orbitals for a given type of atom.
!!
!! PARENTS
!!      m_paw_slater
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine destroy_atomorb(Atm)

 implicit none

!Arguments ------------------------------------
!scalars
 type(atomorb_type),intent(inout) :: Atm

!************************************************************************

 !@atomorb_type

 ! integers
 if (allocated(Atm%indlmn)) then
   ABI_FREE(Atm%indlmn)
 end if
 if (allocated(Atm%indln)) then
   ABI_FREE(Atm%indln)
 end if
 if (allocated(Atm%indklmn)) then
   ABI_FREE(Atm%indklmn)
 end if
 if (allocated(Atm%klm_diag)) then
   ABI_FREE(Atm%klm_diag)
 end if
 if (allocated(Atm%klmntomn)) then
   ABI_FREE(Atm%klmntomn)
 end if
 if (allocated(Atm%kln2ln)) then
   ABI_FREE(Atm%kln2ln)
 end if
 if (allocated(Atm%mode)) then
   ABI_FREE(Atm%mode)
 end if

 !real
 if (allocated(Atm%eig)) then
   ABI_FREE(Atm%eig)
 end if
 if (allocated(Atm%occ)) then
   ABI_FREE(Atm%occ)
 end if
 if (allocated(Atm%phi)) then
   ABI_FREE(Atm%phi)
 end if

end subroutine destroy_atomorb
!!***

!----------------------------------------------------------------------

!!****f* m_paw_atomorb/init_atomorb
!! NAME
!!  init_atomorb
!!
!! FUNCTION
!!  Initialize a structure of type atomorb_type from file.
!!
!! INPUTS
!!  filename=Name of the file containing core electrons
!!
!! OUTPUT
!!  Atmrad<pawrad_type>=Info oh the Radial mesh used for core electrons.
!!  Atm<atomorb_type>=Structure defining the set of core orbitals.
!!  ierr=Status error.
!!   * 1 if error during the opening of the file.
!!   * 2 for generic error during the reading.
!!
!! PARENTS
!!      m_paw_slater
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine init_atomorb(Atm,Atmrad,rcut,filename,prtvol,ierr)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: prtvol
 integer,intent(out) :: ierr
 real(dp),intent(in) :: rcut
 character(len=*),intent(in) :: filename
 type(atomorb_type),intent(out) :: Atm
 type(pawrad_type),intent(out) :: Atmrad

!Local variables-------------------------------
!scalars
 integer :: creatorid
 integer :: imainmesh,jlmn,jl,k0lmn,ilmn,il
 integer :: lcutdens,version,klmn,nn,ll
 integer :: lmax,pspdat
 integer :: ii,unt,iln,ir,nmesh,imsh
 integer :: iread1,pspcod,msz,isppol
 integer :: mesh_type,mesh_size
 real(dp) :: charge
 real(dp) :: my_rcore
 real(dp) :: rstep,lstep
 real(dp):: occ,ene
 character(len=80) :: line,title
 character(len=500) :: msg
!arrays
 integer,allocatable :: orbitals(:)
 real(dp),allocatable :: phitmp(:,:,:),radens(:,:)
 type(pawrad_type),allocatable :: Radmesh(:)

! ************************************************************************

 ! File format of formatted core orbitals.

 !(1) title (character) line
 !(2) method, nspinor,nsppol
 !(3) znucl, zion, pspdat
 !(4) ixc, lmax
 !(5) version, creatorID
 !(6) ln_size, lmn_size
 !(7) orbitals (for i=1 to ln_size)
 !(8) number_of_meshes
 !
 !For imsh=1 to number_of_meshes
 !(9)  mesh_index, mesh_type ,mesh_size, rad_step[, log_step]
 !(10) rcore (SPH)
 !
 !For iln=1 to basis_size
 !  (11) comment(character)
 !  (12) radial mesh index for phi
 !  (13) nn, ll,
 !  (14) phi(r) (for ir=1 to phi_meshsz)
 !
 !Comments:
 ! * allowed values for method are:
 !    1 for restricted, compatible only with nsppol=1.
 !    2 for spin unrestricted, compatible only with nsppol=2.
 !* psp_version= ID of PAW_psp version
 !    4 characters string of the form 'pawn' (with n varying)
 !* creatorID= ID of psp generator
 !  creatorid=1xyz : psp generated from Holzwarth AtomPAW generator version x.yz
 !  creatorid=2xyz : psp generated from Vanderbilt ultra-soft generator version x.yz
 !  creatorid=-1: psp for tests (for developpers only)
 !* mesh_type= type of radial mesh
 !    mesh_type=1 (regular grid): rad(i)=(i-1)*AA
 !    mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
 !    mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
 !    mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
 ! --------------------------------------------------------------------------------

 !@atomorb_type
 ierr=0
 if (open_file(filename,msg,newunit=unt,form="formatted",status="old",action="read") /=0) then
   MSG_WARNING(msg)
   ierr=1; RETURN
 end if

 Atm%fname = filename
 imainmesh = 1

 !1) title
 read(unt,'(a80)',ERR=10)title
 write(msg,'(a,1x,a)')ch10,TRIM(title)
 call wrtout(std_out,msg,'COLL')

 read(unt,*,ERR=10)Atm%method, Atm%nspinor, Atm%nsppol
 write(msg,'(3(i2,2x),22x,a)' )Atm%method, Atm%nspinor, Atm%nsppol,' method, nspinor, nsppol.'
 call wrtout(std_out,msg,'COLL')

 Atm%nspden   = 1 !FIXME pass it as input

 !2)
 read(unt,*,ERR=10) Atm%znucl, Atm%zion, pspdat
 write(msg,'(2f10.5,2x,i8,2x,a)' )Atm%znucl, Atm%zion, pspdat,' znucl, zion, pspdat'
 call wrtout(std_out,msg,'COLL')

 !3)
 read(unt,*,ERR=10)pspcod,Atm%ixc,lmax
 write(msg,'(2i5,2x,2x,a)')Atm%ixc,lmax,'ixc,lmax'
 call wrtout(std_out,msg,'COLL')

 Atm%l_max  =  lmax+1
 Atm%l_size =2*lmax+1

 !4)
!Read psp version in line 4 of the header
 version=1
 read(unt,'(a80)',ERR=10) line
 line=ADJUSTL(line)

 if (tolower(line(1:4))=="core") read(unit=line(5:80),fmt=*) version
 if (version/=1) then
   write(msg,'(a,i2,3a)')&
&   'This version of core psp file (',version,') is not compatible with',ch10,&
&   'the current version of Abinit.'
   MSG_ERROR(msg)
 end if

 read(unit=line(6:80),fmt=*,ERR=10) creatorid

 !==========================================================

 ! 5)
 read(unt,*,ERR=10)Atm%ln_size, Atm%lmn_size

 Atm%ln2_size  = Atm%ln_size *(Atm%ln_size +1)/2
 Atm%lmn2_size = Atm%lmn_size*(Atm%lmn_size+1)/2

 ! 6)
 ABI_MALLOC(orbitals,(Atm%ln_size))
 read(unt,*,ERR=10) (orbitals(iln), iln=1,Atm%ln_size)

 lmax = MAXVAL(orbitals)
 if (lmax+1/=Atm%l_max) then
   write(msg,'(a)')" lmax read from file does not agree with orbitals. "
   MSG_ERROR(msg)
 end if

 ! 7)
 read(unt,*)nmesh
 ABI_MALLOC(Radmesh,(nmesh))
 do imsh=1,nmesh
   lstep=zero
   read(unt,'(a80)') line
   read(unit=line,fmt=*,err=20,end=20) ii,mesh_type,mesh_size,rstep,lstep
   20 continue
   ABI_CHECK(ii<=nmesh,"Index of mesh out of range!")
   call pawrad_init(Radmesh(ii),mesh_size,mesh_type,rstep,lstep,-one)
 end do

 ! 8)
 read(unt,*,ERR=10) my_rcore

 Atm%rcore = my_rcore
 if (rcut > tol6) then
   Atm%rcore = rcut
   write(msg,'(a,f8.5,a,f8.5,a)')&
&    " Truncating radial mesh for core orbitals using new rcore = ",Atm%rcore," (",my_rcore,")"
   call wrtout(std_out,msg,'COLL')
   if (rcut > my_rcore) then
     MSG_ERROR("rcut should not exceed my_rcore")
   end if
 end if

 !==========================================================
 ! Mirror pseudopotential parameters to the output and log files

 write(msg,'(a,i1)')' Pseudopotential format is: core ',version
 call wrtout(std_out,msg,'COLL')
 write(msg,'(2(a,i3),a,64i4)')' (ln_size)= ',Atm%ln_size,' (lmn_size= ',Atm%lmn_size,'), orbitals= ',orbitals(1:Atm%ln_size)
 call wrtout(std_out,msg,'COLL')
 write(msg,'(a,f11.8)')' Radius used for core orbitals: rc_sph= ',Atm%rcore
 call wrtout(std_out,msg,'COLL')
 write(msg,'(a,i1,a)')' ',nmesh,' radial meshes are used:'
 call wrtout(std_out,msg,'COLL')

 do imsh=1,nmesh
   call pawrad_print(Radmesh(imsh),prtvol=prtvol)
 end do

 !==========================================================

 !---------------------------------
 ! (11-12-13) Read wave-functions (phi)
 ! * Initialize also the mesh.
 ABI_MALLOC(Atm%indln,(2,Atm%ln_size))
 ABI_MALLOC(Atm%eig,(Atm%ln_size,Atm%nsppol))
 ABI_MALLOC(Atm%occ,(Atm%ln_size,Atm%nsppol))
 Atm%occ = zero

 do isppol=1,Atm%nsppol
   do iln=1,Atm%ln_size
     read(unt,*,ERR=10) line
     read(unt,*,ERR=10) iread1
     if (iln==1.and.isppol==1) then
       imainmesh=iread1
       Atm%mesh_size = Radmesh(iread1)%mesh_size
       ABI_MALLOC(Atm%phi,(Atm%mesh_size,Atm%ln_size,Atm%nsppol))
     else if (iread1/=imainmesh) then
       write(msg,'(3a)')&
&       ' All Phi core must be given on the same radial mesh !',ch10,&
&       ' Action: check your pseudopotential file.'
       MSG_ERROR(msg)
     end if
     read(unt,*,ERR=10)nn,ll,ii
     ABI_CHECK(ii==isppol,"Wrong spin index")
     read(unt,*,ERR=10)ene,occ
     Atm%indln(1,iln) = ll
     Atm%indln(2,iln) = nn
     Atm%occ(iln,isppol) = occ
     Atm%eig(iln,isppol) = ene
     read(unt,*,ERR=10) (Atm%phi(ir,iln,isppol),ir=1,Radmesh(imainmesh)%mesh_size)
     !do ir=1,Radmesh(imainmesh)%mesh_size
     !  write(100+iln,*)Radmesh(imainmesh)%rad(ir),Atm%phi(ir,iln,isppol)
     !end do
   end do
 end do !isppol

 close(unt)
! -------------------- END OF READING -------------------------------

 if ( rcut < tol16) then
   ! * Use full mesh reported on file.
   call pawrad_init(Atmrad,mesh_size=Radmesh(imainmesh)%mesh_size,mesh_type=Radmesh(imainmesh)%mesh_type,&
&                   rstep=Radmesh(imainmesh)%rstep,lstep=Radmesh(imainmesh)%lstep, &
!&                  r_for_intg=Atm%rcore) ! check this
&                   r_for_intg=-one)
 else
   ! * Truncate orbitals and radial mesh within rcut.
   msz = MIN(pawrad_ifromr(Radmesh(imainmesh),rcut)+6, Radmesh(imainmesh)%mesh_size) ! add six more points

   write(msg,'(a,f8.5,a,f8.5,a)')&
&   " Truncating radial mesh for core orbitals ",Radmesh(imainmesh)%rad(msz),"(",my_rcore,")"
   call wrtout(std_out,msg,'COLL')
   Atm%rcore = Radmesh(imainmesh)%rad(msz)

   ABI_MALLOC(phitmp,(msz,Atm%ln_size,Atm%nsppol))
   phitmp = Atm%phi(1:msz,:,:)

   ABI_FREE(Atm%phi)
   ABI_MALLOC(Atm%phi,(msz,Atm%ln_size,Atm%nsppol))
   Atm%phi = phitmp
   ABI_FREE(phitmp)

   ! * Compute new mesh for core orbitals.
   mesh_type = Radmesh(imainmesh)%mesh_type
   mesh_size = msz !radmesh(imainmesh)%mesh_size
   rstep     = Radmesh(imainmesh)%rstep
   lstep     = Radmesh(imainmesh)%lstep

   call pawrad_init(Atmrad,mesh_size=mesh_size,mesh_type=mesh_type,rstep=rstep,lstep=lstep,&
!&                  r_for_intg=Atm%rcore)
&                   r_for_intg=-one)

   !do isppol=1,Atm%nsppol
   !  do iln=1,Atm%ln_size
   !    do ir=1,Atmrad%mesh_size
   !      write(200+iln,*)Atmrad%rad(ir),Atm%phi(ir,iln,isppol)
   !    end do
   !  end do
   !end do

 end if

 Atm%mesh_size = Atmrad%mesh_size
 call pawrad_print(Atmrad,header="Final mesh",prtvol=prtvol)

 ABI_MALLOC(radens,(Atm%mesh_size,Atm%nspden))
 call get_atomorb_charge(Atm,Atmrad,charge,radens=radens)

 write(std_out,*)"core charge  = ",charge
 !do ii=1,Atmrad%mesh_size
 ! write(77,*)Atmrad%rad(ii),(radens(ii,isppol),isppol=1,Atm%nspden)
 !end do
 ABI_FREE(radens)

 !==========================================================
 ! Free temporary allocated space

 call pawrad_free(Radmesh)
 ABI_FREE(Radmesh)

 ! * Setup of indlmn
 ABI_MALLOC(Atm%indlmn,(6,Atm%lmn_size))
 call make_indlmn(Atm%ln_size, Atm%lmn_size, orbitals, Atm%indlmn)

 ABI_FREE(orbitals)

 ! * Setup of indklmn and klm_diag.
 lcutdens=HUGE(1)
 ABI_MALLOC(Atm%indklmn,(8,Atm%lmn2_size))
 ABI_MALLOC(Atm%klm_diag,(Atm%lmn2_size))

 call make_indklmn(lcutdens, Atm%lmn_size, Atm%lmn2_size, Atm%indlmn, Atm%indklmn, Atm%klm_diag)

 ! * Setup of klmntomn.
 ABI_MALLOC(Atm%klmntomn,(4,Atm%lmn2_size))

 do jlmn=1,Atm%lmn_size
   jl= Atm%indlmn(1,jlmn)
   k0lmn=jlmn*(jlmn-1)/2
   do ilmn=1,jlmn
     il= Atm%indlmn(1,ilmn)
     klmn=k0lmn+ilmn
     Atm%klmntomn(1,klmn) = Atm%indlmn(2,ilmn)+il+1 ! im
     Atm%klmntomn(2,klmn) = Atm%indlmn(2,jlmn)+jl+1 ! jm
     Atm%klmntomn(3,klmn) = Atm%indlmn(3,ilmn)      ! in
     Atm%klmntomn(4,klmn) = Atm%indlmn(3,jlmn)      ! jn
     !write(std_out,*) jlmn,ilmn,Atm%klmntomn(:,klmn)
   end do
 end do

 ! * Setup of kln2ln.
 !TODO this has to be tested
 ABI_MALLOC(Atm%kln2ln,(6,Atm%ln2_size))
 call make_kln2ln(Atm%lmn_size,Atm%lmn2_size,Atm%ln2_size,Atm%indlmn,Atm%indklmn,Atm%kln2ln)

 ! * Treat all states as core.
 ABI_MALLOC(Atm%mode,(Atm%ln_size,Atm%nsppol))
 Atm%mode = ORB_FROZEN

 !call print_atomorb(Atm,prtvol=prtvol)

 return
 !
 ! === Propagate the error ===
10 continue
 ierr=2
 MSG_WARNING("Wrong file format")
 return

end subroutine init_atomorb
!!***

!----------------------------------------------------------------------

!!****f* m_paw_atomorb/get_atomorb_charge
!! NAME
!!  get_atomorb_charge
!!
!! FUNCTION
!!  Get core charge from a structure of type atomorb_type
!!  and optionally core density.
!!
!! INPUTS
!!  Atm<atomorb_type>=Structure defining the set of core orbitals.
!!  Radmesh<pawrad_type>=Info oh the Radial mesh used for core electrons.
!!
!! OUTPUT
!!  nele=core charge
!!  raddens(mesh_size)=core density (optional)
!!
!! PARENTS
!!      m_paw_atomorb
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine get_atomorb_charge(Atm,Radmesh,nele,radens)

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: nele
 type(atomorb_type),intent(in) :: Atm
 type(pawrad_type),intent(in) :: Radmesh
!arrays
 real(dp),optional,intent(out) :: radens(Atm%mesh_size,Atm%nspden)

!Local variables-------------------------------
!scalars
 integer :: iln,isppol
 real(dp) :: intg,focc
 real(dp),allocatable :: phi2nl(:)

!************************************************************************

 if (Atm%nsppol==2) then
   MSG_ERROR("nsppol==2 is Working in progress")
 end if

 ABI_MALLOC(phi2nl,(Atm%mesh_size))
 if (PRESENT(radens)) radens = zero

 nele = zero
 do isppol=1,Atm%nsppol
   do iln=1,Atm%ln_size
     !Atm%mode(iln,isppol) TODO add option to select particular states
     focc   = Atm%occ(iln,isppol)
     if (ABS(focc) > tol16) then
       phi2nl = Atm%phi(:,iln,isppol)**2
       call simp_gen(intg,phi2nl,Radmesh)
       nele = nele + focc*intg
!       if (PRESENT(radens)) then  !FIXME maybe it is better to rr**2 radens
!        radens(2:Atm%mesh_size) = radens(2:Atm%mesh_size) &
!&         + focc * phi2nl(2:Atm%mesh_size)/(four_pi*Radmesh%rad(2:Atm%mesh_size)**2)
!       end if
     end if
   end do
 end do

 ABI_FREE(phi2nl)

end subroutine get_atomorb_charge
!!***

!----------------------------------------------------------------------

!!****f* m_paw_atomorb/get_overlap
!! NAME
!!  get_overlap
!!
!! FUNCTION
!!  Get overlap between core and valence states
!!
!! INPUTS
!!  Atm<atomorb_type>=Structure defining the set of core states
!!  Atmesh<pawrad_type>=Info oh the Radial mesh used for core states
!!  isppol=index for spin component
!!  nphi=number of core states
!!  phi(Radmesh2%mesh_size,nphi)=valence states
!!  phi_indln(nphi)=Array giving l and and n for i=1,nphi
!!  Radmesh2<pawrad_type>=Info oh the Radial mesh used for valence states
!!
!! OUTPUT
!!  overlap(ln_size,nphi)=core-valence overlap matrix
!!
!! PARENTS
!!      m_paw_slater
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine get_overlap(Atm,Atmesh,Radmesh2,isppol,nphi,phi,phi_indln,overlap)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nphi,isppol
 type(atomorb_type),intent(in) :: Atm
 type(pawrad_type),target,intent(in) :: Atmesh,Radmesh2
!arrays
 integer,intent(in) :: phi_indln(2,nphi)
 real(dp),target,intent(in) :: phi(Radmesh2%mesh_size,nphi)
 real(dp),intent(out) :: overlap(Atm%ln_size,nphi)

!Local variables-------------------------------
!scalars
 integer :: iln_atm,iphi,ll_phi,ll_atm,do_spline,iln
 integer :: whichdenser,size4spl,my_mesh_size
 real(dp) :: ybcbeg,ybcend,intg
 logical :: hasameq
!arrays
 real(dp),pointer :: ff_spl(:,:)
 real(dp),allocatable :: der(:),ypp(:),func(:)
 real(dp),pointer :: rad4spl(:),my_pts(:)

!************************************************************************

 ABI_CHECK(isppol>0.and.isppol<=Atm%nsppol,"Wrong isppol")

 call pawrad_isame(Atmesh,Radmesh2,hasameq,whichdenser)

 do_spline= 0; if (.not.hasameq) do_spline=1

 my_mesh_size = MIN(Atmesh%mesh_size,Radmesh2%mesh_size)
 ff_spl => phi

 ! === Spline valence onto Atom mesh (natural spline) ===
 if (do_spline==1) then
   MSG_COMMENT("Splining in overlap")

   my_mesh_size  =  Atmesh%mesh_size
   my_pts        => Atmesh%rad(1:my_mesh_size)
   ABI_MALLOC(ff_spl,(my_mesh_size,nphi))

   size4spl =  Radmesh2%mesh_size
   rad4spl  => Radmesh2%rad
   ABI_MALLOC(der,(size4spl))
   ABI_MALLOC(ypp,(size4spl))

   do iln=1,nphi
     ypp(:) = zero; ybcbeg = zero; ybcend = zero
     call spline(rad4spl,phi(:,iln),size4spl,ybcbeg,ybcend,ypp)

     call splint(size4spl,rad4spl,phi(:,iln),ypp,my_mesh_size,my_pts,ff_spl(:,iln))
   end do

   ABI_FREE(der)
   ABI_FREE(ypp)
 end if

 ABI_MALLOC(func,(my_mesh_size))
 overlap = zero

 do iphi=1,nphi
   ll_phi = phi_indln(1,iphi)
   do iln_atm=1,Atm%ln_size
     ll_atm = Atm%indln(1,iln_atm)

     if (ll_atm == ll_phi) then ! selection rule on l
       func(:) = Atm%phi(1:my_mesh_size,iln_atm,isppol) * ff_spl(1:my_mesh_size,iphi)
       call simp_gen(intg,func,Atmesh)
       overlap(iln_atm,iphi)=intg
       write(std_out,*)"overlap <phic_i|phi_j> for ll_phi",ll_phi,"ll_phic",ll_atm,"=",intg
     end if

   end do
 end do
 ABI_FREE(func)

 if (do_spline==1)  then
   ABI_FREE(ff_spl)
 end if

end subroutine get_overlap
!!***

!----------------------------------------------------------------------

!!****f* m_paw_atomorb/print_atomorb
!! NAME
!!  print_atomorb
!!
!! FUNCTION
!!  Reports info on a structure of type atomorb_type.
!!
!! INPUTS
!!  Atm <type(atomorb_type)>=datastructure containing atomic orbitals for a given type of atom.
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_slater
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine print_atomorb(Atm,header,unit,prtvol,mode_paral)

 implicit none

!Arguments ------------------------------------
!scalars
 type(atomorb_type),intent(in) :: Atm
 integer,optional,intent(in) :: prtvol,unit
 character(len=*),optional,intent(in) :: header
 character(len=4),optional,intent(in) :: mode_paral

!Local variables-------------------------------
 integer :: my_unt,my_prtvol,iln,ll,nn,isppol
 character(len=4) :: my_mode
 character(len=500) :: msg
! ************************************************************************

 !@atomorb_type
 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Info on the atomorb_type ==== '
 if (PRESENT(header)) msg=header
 call wrtout(my_unt,msg,my_mode)

 select case (Atm%method)
 case (1)
   msg = "  Spin restricted"
 case(2)
   msg = "  Spin unrestricted"
 case default
   write(msg,'(a,i3)')" Wrong method= ",Atm%method
   MSG_BUG(msg)
 end select
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(7(a,i5,a),(a,f8.5,a))')&
& '  Number of spinorial components ...... ',Atm%nspinor,ch10,&
& '  Number of ind. spin polarizations ... ',Atm%nsppol,ch10,&
& '  Number of spin-density components ... ',Atm%nspden,ch10,&
& '  Maximum angular momentum + 1 ........ ',Atm%l_max,ch10,&
& '  Number of (l,n) orbitals  ........... ',Atm%ln_size,ch10,&
& '  Number of (l,m,n) orbitals  ......... ',Atm%lmn_size,ch10,&
& '  Dimensions of radial mesh ........... ',Atm%mesh_size,ch10,&
& '  Core Radius  ........................ ',Atm%rcore,ch10
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(2(a,f8.5,a))')&
& '  Ionic charge ........................ ',Atm%zion,ch10,&
& '  Atomic number ....................... ',Atm%znucl,ch10
 call wrtout(my_unt,msg,my_mode)

 do isppol=1,Atm%nsppol
   do iln=1,Atm%ln_size
     ll = Atm%indln(1,iln)
     nn = Atm%indln(2,iln)
     write(msg,'(" n=",i2,", l=",i2,", spin=",i2,", nocc=",f15.7,", energy=",f15.7,2x,"(",a,")")')&
&      nn,ll,isppol,Atm%occ(iln,isppol),Atm%eig(iln,isppol),TRIM(my_mode2str(Atm%mode(iln,isppol)))
     call wrtout(my_unt,msg,my_mode)
   end do
 end do

end subroutine print_atomorb
!!***

!----------------------------------------------------------------------

!!****f* m_paw_atomorb/my_mode2str
!! NAME
!!  my_mode2str
!!
!! FUNCTION
!!  Converts an integer flags defining the way an orbital is treated to a string.
!!
!! INPUTS
!!  mode=Integer
!!
!! OUTPUT
!!  str=mode. Either "Frozen", "Relazed Core", "Valence"
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function my_mode2str(mode) result(str)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mode
 character(len=50) :: str

!Local variables
 character(len=500) :: msg

!************************************************************************

 select case (mode)
 case (ORB_FROZEN)
   str="Frozen Orbital"
 case (ORB_RELAXED_CORE)
   str="Relaxed Core Orbital"
 case (ORB_VALENCE)
   str="Valence Orbital"
 case default
   write(msg,'(a,i3)')" Wrong mode= ",mode
   MSG_BUG(msg)
 end select

end function my_mode2str
!!***

END MODULE m_paw_atomorb
!!***
