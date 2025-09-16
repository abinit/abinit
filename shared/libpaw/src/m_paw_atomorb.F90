!!****m* ABINIT/m_paw_atomorb
!! NAME
!!  m_paw_atomorb
!!
!! FUNCTION
!!  This module provides the definition of the atomorb_type used
!!  to store atomic orbitals on a radial mesh as well
!!  as methods to operate on it.
!!
!! Copyright (C) 2008-2025 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#include "libpaw.h"

MODULE m_paw_atomorb

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 use m_paw_numeric
 use m_libpaw_tools,  only : libpaw_basename, libpaw_get_free_unit
 use m_pawrad,        only : pawrad_type, pawrad_init, bound_deriv, &
&                            pawrad_free, pawrad_print, pawrad_isame, pawrad_ifromr, simp_gen,pawrad_deducer0

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

  integer :: mult
  ! Number of atoms of the same typat

  logical :: dirac
  ! Dirac relativism or not

  logical :: nc_conv
  ! nc has converged ?

  logical :: zcore_conv
  ! zcore has converged ?

  real(dp) :: edcc
  ! DC core energy

  real(dp) :: eeigc
  ! Core eigenvalue energy contribution

  real(dp) :: ehnzc
  ! Core Hartree nc+Z energy

  real(dp) :: ekinc
  ! Core kinetic energy

  real(dp) :: min_eigv
  ! Minimal eigenvalue of the  valence orbitals

  real(dp) :: nresid_c
  ! Residual error on core density

  real(dp) :: rcore
  ! Radius of the sphere used to describe core electrons.
  ! It should be <= rpaw

  real(dp) :: zion
   ! zionpsp
   ! The ionic pseudo-charge, (giving raise to a long-range coulomb potential)

  real(dp) :: zcore
   ! Number of core electrons

  real(dp) :: zcore_orig
   ! original zcore
   ! This is used in RCPAW

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

  type(pawrad_type) :: radmesh
  ! Radial mesh

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

  integer, allocatable :: kappa(:)
  ! Kappa for dirac relativism

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

  real(dp), allocatable :: max_occ(:,:)
  ! max_occ(ln_size,nsppol)
  ! Maximal occupancy for each, used in RCPAW

  real(dp), allocatable :: occ(:,:)
  ! occ(ln_size,nsppol)
  ! Occupation for each ln channel and spin.

  real(dp), allocatable :: occ_res(:,:)
  ! occ_res(ln_size,nsppol)
  ! Occupation residue for each ln channel and spin, used in RCPAW

  real(dp), allocatable :: occ_respc(:,:)
  ! occ_respc(ln_size,nsppol)
  ! Occupation preconditionned residue for each ln channel and spin, used in RCPAW.

  real(dp), allocatable :: phi(:,:,:)
  ! phi(mesh_size,ln_size,nsppol)
  ! Here we might have different meshes, useful for deep states in heavy atoms!

  ! this might be retrieved with  a method get_atomden
  !real(dp), allocatable :: density(:)
   ! density(mesh_size,nspden)
   ! Gives the core density of the atom for each spin channel
   ! Total charge in first dimension,up component in second one (if present)

  real(dp), allocatable :: vhtnzc_orig(:)
   ! vhtnzc_orig(size(pawtab(itypat)%vhtnzc))
   ! Original vhtnzc, used in RCPAW

 end type atomorb_type

! public procedures.
 public :: destroy_atomorb
 public :: print_atomorb
 public :: get_overlap
!!***

!----------------------------------------------------------------------

 integer,public,parameter :: ORB_FROZEN       =1
 integer,public,parameter :: ORB_RELAXED_CORE =0
 integer,public,parameter :: ORB_VALENCE      =2


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
!! SOURCE

subroutine destroy_atomorb(Atm)

!Arguments ------------------------------------
!scalars
 type(atomorb_type),intent(inout) :: Atm

!************************************************************************

 !@atomorb_type

 ! integers
 if (allocated(Atm%indlmn)) then
   LIBPAW_DEALLOCATE(Atm%indlmn)
 end if
 if (allocated(Atm%indln)) then
   LIBPAW_DEALLOCATE(Atm%indln)
 end if
 if (allocated(Atm%indklmn)) then
   LIBPAW_DEALLOCATE(Atm%indklmn)
 end if
 if (allocated(Atm%klm_diag)) then
   LIBPAW_DEALLOCATE(Atm%klm_diag)
 end if
 if (allocated(Atm%klmntomn)) then
   LIBPAW_DEALLOCATE(Atm%klmntomn)
 end if
 if (allocated(Atm%kln2ln)) then
   LIBPAW_DEALLOCATE(Atm%kln2ln)
 end if
 if (allocated(Atm%kappa)) then
   LIBPAW_DEALLOCATE(Atm%kappa)
 end if
 if (allocated(Atm%mode)) then
   LIBPAW_DEALLOCATE(Atm%mode)
 end if

 !real
 if (allocated(Atm%eig)) then
   LIBPAW_DEALLOCATE(Atm%eig)
 end if
 if (allocated(Atm%max_occ)) then
   LIBPAW_DEALLOCATE(Atm%max_occ)
 end if
 if (allocated(Atm%occ)) then
   LIBPAW_DEALLOCATE(Atm%occ)
 end if
 if (allocated(Atm%occ_res)) then
   LIBPAW_DEALLOCATE(Atm%occ_res)
 end if
if (allocated(Atm%occ_respc)) then
   LIBPAW_DEALLOCATE(Atm%occ_respc)
 end if
 if (allocated(Atm%phi)) then
   LIBPAW_DEALLOCATE(Atm%phi)
 end if
 if(allocated(atm%vhtnzc_orig)) then
   LIBPAW_DEALLOCATE(atm%vhtnzc_orig)
 endif
 call pawrad_free(atm%radmesh)

end subroutine destroy_atomorb
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
!! SOURCE

subroutine get_atomorb_charge(Atm,Radmesh,nele,radens)

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
   LIBPAW_ERROR("nsppol==2 is Working in progress")
 end if

 LIBPAW_ALLOCATE(phi2nl,(Atm%mesh_size))
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

 LIBPAW_DEALLOCATE(phi2nl)

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
!! SOURCE

subroutine get_overlap(Atm,Atmesh,Radmesh2,isppol,nphi,phi,phi_indln,overlap)

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

 if(isppol<=0.or.isppol>Atm%nsppol) LIBPAW_ERROR("Wrong isppol")

 call pawrad_isame(Atmesh,Radmesh2,hasameq,whichdenser)

 do_spline= 0; if (.not.hasameq) do_spline=1

 my_mesh_size = MIN(Atmesh%mesh_size,Radmesh2%mesh_size)
 ff_spl => phi

 ! === Spline valence onto Atom mesh (natural spline) ===
 if (do_spline==1) then
   LIBPAW_COMMENT("Splining in overlap")

   my_mesh_size  =  Atmesh%mesh_size
   my_pts        => Atmesh%rad(1:my_mesh_size)
   LIBPAW_ALLOCATE(ff_spl,(my_mesh_size,nphi))

   size4spl =  Radmesh2%mesh_size
   rad4spl  => Radmesh2%rad
   LIBPAW_ALLOCATE(der,(size4spl))
   LIBPAW_ALLOCATE(ypp,(size4spl))

   do iln=1,nphi
     ypp(:) = zero; ybcbeg = zero; ybcend = zero
     call paw_spline(rad4spl,phi(:,iln),size4spl,ybcbeg,ybcend,ypp)

     call paw_splint(size4spl,rad4spl,phi(:,iln),ypp,my_mesh_size,my_pts,ff_spl(:,iln))
   end do

   LIBPAW_DEALLOCATE(der)
   LIBPAW_DEALLOCATE(ypp)
 end if

 LIBPAW_ALLOCATE(func,(my_mesh_size))
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
 LIBPAW_DEALLOCATE(func)

 if (do_spline==1)  then
   LIBPAW_DEALLOCATE(ff_spl)
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
!! SOURCE

subroutine print_atomorb(Atm,header,unit,prtvol,mode_paral)

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
   LIBPAW_BUG(msg)
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
!! SOURCE

function my_mode2str(mode) result(str)

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
   ABI_BUG(msg)
 end select

end function my_mode2str
!!***

END MODULE m_paw_atomorb
!!***
