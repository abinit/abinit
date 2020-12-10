!!****m* ABINIT/m_pawrad
!! NAME
!! m_pawrad
!!
!! FUNCTION
!! Module containing all the functions related to the PAW radial meshes
!!
!! COPYRIGHT
!! Copyright (C) 2013-2020 ABINIT group (MT,FJ,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  * Routines tagged with "@type_name" are strongly connected to the definition of the data type.
!!    Strongly connected means that the proper functioning of the implementation relies on the
!!    assumption that the tagged procedure is consistent with the type declaration.
!!    Every time a developer changes the structure "type_name" adding new entries, he/she has to make sure
!!    that all the strongly connected routines are changed accordingly to accommodate the modification of the data type
!!    Typical examples of strongly connected routines are creation, destruction or reset methods.
!!
!!  * FOR DEVELOPERS: in order to preserve the portability of libPAW library,
!!    please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

MODULE m_pawrad

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 use m_paw_numeric, only : paw_derfc

 implicit none

 private

!public procedures.
 public :: pawrad_init          ! Main creation method
 public :: pawrad_free          ! Free the allocated memory
 public :: pawrad_print         ! Printout of the basic info
 public :: pawrad_isame         ! Checks whether two meshes are equivalent or have the same equation.
 public :: pawrad_copy          ! Returns a copy of the mesh.
 public :: pawrad_ifromr        ! Retrieve the Index FROM a given R value in a radial grid.
 public :: pawrad_deducer0      ! Extrapolate r=0 value of a function from values near r=0.
 public :: pawrad_bcast         ! Broadcast pawrad datastructure over a given MPI communicator
 public :: simp_gen             ! Performs integral on a given (generalized) grid using Simpson rule.
 public :: nderiv_gen           ! Do corrected first (and 2nd) derivation on a given (generalized) grid.
 public :: nderiv_lin           ! Do corrected first (and 2nd) derivation on a given linear grid.
 public :: bound_deriv          ! Computes derivatives at boundaries of the mesh
 public :: poisson              ! Solves Poisson eq. for angularly dependent charge distribution of angular momentum l
 public :: screened_coul_kernel ! Kernel used to compute short-range screened Coulomb integrals
 public :: calc_slatradl        ! Calculates the radial part of Slater integrals.

 interface pawrad_free
   module procedure pawrad_free_0D
   module procedure pawrad_free_1D
 end interface pawrad_free

 ! TODO: Might use bit flags, but all radmesh stuff should be encapsulated here!
 integer,private,parameter :: RMESH_LINEAR = 1
 integer,private,parameter :: RMESH_LOG1   = 2
 integer,private,parameter :: RMESH_LOG2   = 3
 integer,private,parameter :: RMESH_LOG3   = 4
 integer,private,parameter :: RMESH_NL     = 5
!!***

!-------------------------------------------------------------------------

!!****t* m_pawrad/pawrad_type
!! NAME
!! pawrad_type
!!
!! FUNCTION
!! For PAW, RADial mesh discretization and related data
!!
!! SOURCE

 type, public :: pawrad_type

!Integer scalars

  integer :: int_meshsz=0
   ! Mesh size used in integrals computation
   ! Integrals will be computed up to r(int_meshsz)

  integer :: mesh_size=0
   ! Dimension of radial mesh

  integer :: mesh_type=-1
   ! Type of mesh
   !     1=regular grid: r(i)=(i-1)*AA
   !     2=logarithmic grid: r(i)=AA*(exp[BB*(i-1)]-1)
   !     3=logarithmic grid: r(i>1)=AA*exp[BB*(i-1)] and r(1)=0
   !     4=logarithmic grid: r(i)=-AA*ln[1-BB*(i-1)] with BB=1/n

!Real (real(dp)) scalars

  real(dp) :: lstep=zero
   ! Exponential step of the mesh (BB parameter above)
   ! Defined only if mesh type is logarithmic

  real(dp) :: rmax=zero
   ! Max. value of r = rad(mesh_size)

  real(dp) :: rstep=zero
   ! Radial step of the mesh (AA parameter above)

  real(dp) :: stepint=zero
   ! Radial step used to convert any function from the
   ! present grid onto a regular grid in order to
   ! integrate it using trapeze method

!Real (real(dp)) arrays

  real(dp), allocatable :: rad(:)
   ! rad(mesh_size)
   ! Coordinates of all the points of the mesh

  real(dp), allocatable :: radfact(:)
   ! radfact(mesh_size)
   ! Factor used to compute radial integrals
   ! Before being integrated on the present mesh,
   ! any function is multiplied by this factor

  real(dp), allocatable :: simfact(:)
   ! simfact(mesh_size)
   ! Factor used to compute radial integrals by the a Simpson scheme
   ! Integral[f] = Sum_i [simfact(i)*f(i)]

 end type pawrad_type
!!***

CONTAINS
!===========================================================
!!***

!!****f* m_pawrad/pawrad_init
!! NAME
!!  pawrad_init
!!
!! FUNCTION
!!  Creation method for radial meshes.
!!  Compute all points (and related weights) of a radial mesh.
!!  Grid can be regular or logarithimc.
!!
!! INPUTS
!!  [mesh_size]=Dimension of the radial mesh
!!              If not present, take mesh%mesh_size
!!  [mesh_type]=Type of mesh
!!              If not present, take mesh%mesh_type
!!  [rstep]=Radial step of the mesh (AA parameter above)
!!          If not present, take mesh%rstep
!!  [lstep]=Exponential step of the mesh (BB parameter above)
!!          If not present, take mesh%lstep
!!          Needed only if mesh type is logarithmic.
!!  [r_for_intg]=Mesh size used in integrals computation
!!               If not present, take mesh%r_for_intg
!!    Integrals will be computed up to rr(r_for_intg)
!!    (can be negative for an integration over the whole grid)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  mesh<pawrad_type>=The object completely initialized (containing radial grid information).
!!  The following quantities are calculated inside the routine:
!!    %stepint = Radial step used to convert any function from the
!!               present grid onto a regular grid in order to integrate it using trapeze method
!!    %rad(mesh_size) = Coordinates of all the points of the mesh.
!!    %radfact(mesh_size) = Factors used to compute radial integrals.
!!    %int_meshsz = Integrals will be computed up to r(int_meshsz)
!!    %simfact(mesh_size) = Factor used to compute radial integrals by the a Simpson scheme
!!                          Integral[f] = Sum_i [simfact(i)*f(i)]
!!    %rmax = Max. value of r = rad(mesh_size)
!!
!! PARENTS
!!      m_dfpt_elt,m_mkcore,m_paw_atomorb,m_paw_gaussfit,m_pawpsp,m_pawpwij
!!      m_pawxmlps,m_psp8,m_psp9,m_psps,m_wvl_rho,mkcore_wvl
!!
!! CHILDREN
!!      poisson,simp_gen
!!
!! NOTES
!!  Possible mesh types (mesh%mesh_type)
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!   mesh_type=5 ( grid): rad(i)=AA*i/(n-i)
!!
!! SOURCE

subroutine pawrad_init(mesh,mesh_size,mesh_type,rstep,lstep,r_for_intg)

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: mesh_size,mesh_type
 real(dp),intent(in),optional :: rstep,lstep
 real(dp),intent(in),optional :: r_for_intg
 type(pawrad_type),intent(inout) :: mesh

!Local variables-------------------------------
!scalars
 integer :: ir,ir_last,isim,mesh_size_,mesh_type_
 real(dp) :: hh,lstep_,rstep_,r_for_intg_
 character(len=500) :: msg

! *************************************************************************

 !@pawrad_type

!Retrieve mesh data
 mesh_size_ = mesh%mesh_size ; if (present(mesh_size)) mesh_size_ = mesh_size
 mesh_type_ = mesh%mesh_type ; if (present(mesh_type)) mesh_type_ = mesh_type
 rstep_ = mesh%rstep ; if (present(rstep)) rstep_ = rstep
 lstep_ = mesh%lstep ; if (present(lstep)) lstep_ = lstep

 r_for_intg_=-1._dp;if (present(r_for_intg)) r_for_intg_=r_for_intg

 mesh%mesh_size = mesh_size_
 mesh%mesh_type = mesh_type_
 mesh%rstep     = rstep_
 mesh%lstep     = lstep_
 LIBPAW_ALLOCATE(mesh%rad    ,(mesh%mesh_size))
 LIBPAW_ALLOCATE(mesh%radfact,(mesh%mesh_size))
 LIBPAW_ALLOCATE(mesh%simfact,(mesh%mesh_size))
 mesh%simfact=zero
 if (mesh%mesh_type==1) then
   isim=3
   mesh%stepint=mesh%rstep
   mesh%rad(1)=zero;mesh%radfact(1)=one
   do ir=2,mesh%mesh_size
     mesh%rad(ir) =mesh%rstep*dble(ir-1)
     mesh%radfact(ir)=one
   end do
 else if (mesh%mesh_type==2) then
   isim=3
   mesh%stepint=mesh%lstep
   mesh%rad(1)=zero;mesh%radfact(1)=mesh%rstep
   do ir=2,mesh%mesh_size
     mesh%rad(ir) =mesh%rstep*(exp(mesh%lstep*dble(ir-1))-one)
     mesh%radfact(ir)=mesh%rad(ir)+mesh%rstep
   end do
 else if (mesh%mesh_type==3) then
   isim=4
   mesh%stepint=mesh%lstep
   mesh%rad(1)=zero;mesh%radfact(1)=zero
   do ir=2,mesh%mesh_size
     mesh%rad(ir) =mesh%rstep*exp(mesh%lstep*dble(ir-2))
     mesh%radfact(ir)=mesh%rad(ir)
   end do
 else if (mesh%mesh_type==4) then
   isim=3
   mesh%lstep=one/dble(mesh%mesh_size)
   mesh%stepint=mesh%lstep
   mesh%rad(1)=zero;mesh%radfact(1)=mesh%rstep
   do ir=2,mesh%mesh_size
     mesh%rad(ir) =-mesh%rstep*log(one-mesh%lstep*dble(ir-1))
     mesh%radfact(ir)=mesh%rstep/(one-mesh%lstep*dble(ir-1))
   end do
 else if (mesh%mesh_type==5) then
   isim=3
   mesh%stepint=mesh%rstep
   mesh%rad(1)=zero;mesh%radfact(1)=1/mesh%lstep
   do ir=2,mesh%mesh_size
     mesh%rad(ir) =mesh%rstep*dble(ir-1)/(mesh%lstep-dble(ir-1))
     mesh%radfact(ir)=(mesh%rstep+mesh%rad(ir))/(mesh%lstep-dble(ir-1))/mesh%rstep
   end do

 else !  Other values of mesh_type are not allowed (see psp7in.F90)
  write(msg,'(a,i0)')" Unknown value of mesh_type: ",mesh%mesh_type
  LIBPAW_ERROR(msg)
 end if

 mesh%int_meshsz=mesh%mesh_size
 if (r_for_intg_>0.d0) then
   ir=min(pawrad_ifromr(mesh,r_for_intg_),mesh%mesh_size)
   if (ir<mesh%mesh_size) then
     if (abs(mesh%rad(ir+1)-r_for_intg_)<abs(mesh%rad(ir)-r_for_intg_)) ir=ir+1
   end if
   if (ir>1) then
     if (abs(mesh%rad(ir-1)-r_for_intg_)<abs(mesh%rad(ir)-r_for_intg_)) ir=ir-1
   end if
   mesh%int_meshsz=ir
 end if

 hh=mesh%stepint/3.d0
 mesh%simfact(mesh%int_meshsz)=hh*mesh%radfact(mesh%int_meshsz)
 mesh%simfact(1:isim-2)=zero
 ir_last=1
 do ir=mesh%int_meshsz,isim,-2
   mesh%simfact(ir-1)=4.d0*hh*mesh%radfact(ir-1)
   mesh%simfact(ir-2)=2.d0*hh*mesh%radfact(ir-2)
   ir_last=ir-2
 end do
 mesh%simfact(ir_last)=half*mesh%simfact(ir_last)
 if (mesh%int_meshsz<mesh%mesh_size) mesh%simfact(mesh%int_meshsz+1:mesh%mesh_size)=zero

 mesh%rmax=mesh%rad(mesh%mesh_size)

end subroutine pawrad_init
!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/pawrad_free_0D
!! NAME
!!  pawrad_free_0D
!!
!! FUNCTION
!!  Frees all memory allocated in the object
!!
!! PARENTS
!!      m_pawrad
!!
!! CHILDREN
!!      poisson,simp_gen
!!
!! SOURCE

subroutine pawrad_free_0D(Rmesh)

!Arguments ------------------------------------
!arrays
 type(pawrad_type),intent(inout) :: Rmesh

!Local variables-------------------------------

! *************************************************************************

 !@Pawrad_type

 if (allocated(Rmesh%rad    ))  then
   LIBPAW_DEALLOCATE(Rmesh%rad)
 end if
 if (allocated(Rmesh%radfact))  then
   LIBPAW_DEALLOCATE(Rmesh%radfact)
 end if
 if (allocated(Rmesh%simfact))  then
   LIBPAW_DEALLOCATE(Rmesh%simfact)
 end if
 Rmesh%int_meshsz=0
 Rmesh%mesh_size=0
 Rmesh%mesh_type=-1

end subroutine pawrad_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/pawrad_free_1D
!! NAME
!!  pawrad_free_1D
!!
!! FUNCTION
!!  Destroy all objects in an array of pawrad data structures
!!
!! PARENTS
!!
!! CHILDREN
!!      poisson,simp_gen
!!
!! SOURCE

subroutine pawrad_free_1D(Rmesh)

!Arguments ------------------------------------
 type(pawrad_type),intent(inout) :: Rmesh(:)

!Local variables-------------------------------
 integer :: ii,nn

! *************************************************************************

 !@pawrad_type

 nn=size(Rmesh)
 if (nn==0) return

 do ii=1,nn
   call pawrad_free_0D(Rmesh(ii))
 end do

end subroutine pawrad_free_1D
!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/pawrad_print
!! NAME
!!  pawrad_print
!!
!! FUNCTION
!!  Reports basic info on the object.
!!
!! INPUTS
!!  Rmesh<pawrad_type>=Object defining the radial mesh
!!  header=String for the header provided by the user.
!!  [unit]=Unit number for output, defaults to std_out
!!  [prtvol]=Verbosity level, minimal if not specified.
!!  [mode_paral]=Either "COLL" or "PERS". Passed to wrtout. Defaults to "COLL"
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_paw_atomorb
!!
!! CHILDREN
!!      poisson,simp_gen
!!
!! SOURCE

subroutine pawrad_print(Rmesh,header,unit,prtvol,mode_paral)

!Arguments ------------------------------------
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 character(len=*),intent(in),optional :: header
 type(pawrad_type),intent(in) :: Rmesh

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 !@pawrad_type
 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=ch10//' ==== Info on the Radial Mesh ==== '
 if (PRESENT(header)) msg=ch10//' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 SELECT CASE (Rmesh%mesh_type)

 CASE (RMESH_LINEAR)
   write(msg,'(a,i4,a,g12.5)')&
&    ' - Linear mesh: r(i)=step*(i-1), size=',Rmesh%mesh_size,', step=',Rmesh%rstep

 CASE (RMESH_LOG1)
   write(msg,'(a,i4,2(a,g12.5))')&
&    ' - Logarithimc mesh: r(i)=AA*[exp(BB*(i-1))-1], size=',Rmesh%mesh_size,', AA=',Rmesh%rstep,' BB=',Rmesh%lstep

 CASE (RMESH_LOG2)
   write(msg,'(a,i4,2(a,g12.5))')&
&    ' - Logarithimc mesh: r(i)=AA*exp(BB*(i-2)), size=',Rmesh%mesh_size,', AA=',Rmesh%rstep,' BB=',Rmesh%lstep

 CASE (RMESH_LOG3)
   write(msg,'(a,i1,a,i4,a,g12.5)')&
&    ' - Logarithimc mesh: r(i)=-AA*ln(1-(i-1)/n), n=size=',Rmesh%mesh_size,', AA=',Rmesh%rstep

 CASE (RMESH_NL)
   write(msg,'(a,i1,a,i4,a,g12.5)')&
&    ' - Non-linear mesh: r(i)=-AA*i/(n-i), n=size=',Rmesh%mesh_size,', AA=',Rmesh%rstep

 CASE DEFAULT
   msg = ' Unknown mesh type! Action : check your pseudopotential or input file.'
   LIBPAW_ERROR(msg)
 END SELECT

 call wrtout(my_unt,msg,my_mode)

 if (my_prtvol>1) then
   write(msg,'(a,i4)')' Mesh size for integrals = ',Rmesh%int_meshsz
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,g12.5)')' rmax=rad(mesh_size)     = ',Rmesh%rmax
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(a,g12.5)')' Value of stepint        = ',Rmesh%stepint
   call wrtout(my_unt,msg,my_mode)
 end if

end subroutine pawrad_print
!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/pawrad_isame
!! NAME
!!  pawrad_isame
!!
!! FUNCTION
!!  Check two radial meshes, returns a logical flag defining
!!  whether the meshes have the same equation and an integer
!!  flag
!!
!! INPUTS
!!  Rmesh1,Rmesh2<pawrad_type>=The two radial meshes.
!!
!! OUTPUT
!!  hasameq=.true. if the two meshes are defined by the same equation.
!!  whichdenser=
!!    * 0 if meshes are not compatible
!!    * 1 if Rmesh1 is denser than Rmesh2
!!    * 2 if Rmesh2 is denser than Rmesh1
!!
!! PARENTS
!!      m_paw_atomorb,m_paw_slater
!!
!! CHILDREN
!!      poisson,simp_gen
!!
!! SOURCE

subroutine pawrad_isame(Rmesh1,Rmesh2,hasameq,whichdenser)

!Arguments ------------------------------------
 integer,intent(out) :: whichdenser
 logical,intent(out) :: hasameq
 type(pawrad_type),intent(in) :: Rmesh1
 type(pawrad_type),intent(in) :: Rmesh2

!Local variables-------------------------------
 character(len=50) :: msg

! *************************************************************************

 !@pawrad_type

 whichdenser =0 ; hasameq=.FALSE.

 if (Rmesh1%mesh_type /= Rmesh2%mesh_type) RETURN

 SELECT CASE (Rmesh1%mesh_type)

 CASE (RMESH_LINEAR) !check for linear meshes
   hasameq = (Rmesh1%rstep == Rmesh2%rstep)

 CASE (RMESH_LOG1,&  !check for logarithmic meshes
&      RMESH_LOG2,&
&      RMESH_LOG3)

   hasameq = (    Rmesh1%rstep == Rmesh2%rstep &
&            .and.Rmesh1%lstep == Rmesh2%lstep )

 CASE (RMESH_NL) !check for linear meshes
   hasameq = (Rmesh1%rstep == Rmesh2%rstep)

 CASE DEFAULT
   msg='Unknown mesh type'
   LIBPAW_ERROR(msg)

 END SELECT

 ! === If meshes have same equation, check whether they are equal ===
 ! * Note that also int_meshsz must be equal
 if (hasameq) then
   whichdenser= 1
   if (Rmesh2%mesh_size  > Rmesh1%mesh_size ) whichdenser = 2
   !if (Rmesh1%int_meshsz == Rmesh2%int_meshsz) whichdenser = 2
 end if

end subroutine pawrad_isame
!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/pawrad_copy
!! NAME
!! pawrad_copy
!!
!! FUNCTION
!! Copy one radial mesh (in a generalized format) to another
!!
!! INPUTS
!!  mesh1 <type(pawrad_type)>=data containing radial grid information of input mesh
!!
!! OUTPUT
!!  mesh2 <type(pawrad_type)>=data containing radial grid information of output mesh
!!
!! PARENTS
!!      m_pawpsp,m_pawpwij
!!
!! CHILDREN
!!      poisson,simp_gen
!!
!! NOTES
!!  Possible mesh types (mesh%mesh_type)
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!   mesh_type=5 ( grid): rad(i)=AA*i/(n-i)
!!
!! SOURCE

subroutine pawrad_copy(mesh1,mesh2)

!Arguments ------------------------------------
!scalars
 type(pawrad_type),intent(in) :: mesh1
 type(pawrad_type),intent(out) :: mesh2

!Local variables-------------------------------
!scalars
 integer :: ir

! *************************************************************************

 mesh2%mesh_type =mesh1%mesh_type
 mesh2%mesh_size =mesh1%mesh_size
 mesh2%int_meshsz=mesh1%int_meshsz
 mesh2%lstep     =mesh1%lstep
 mesh2%rstep     =mesh1%rstep
 mesh2%stepint   =mesh1%stepint
 mesh2%rmax      =mesh1%rmax

 LIBPAW_ALLOCATE(mesh2%rad,(mesh1%mesh_size))
 LIBPAW_ALLOCATE(mesh2%radfact,(mesh1%mesh_size))
 LIBPAW_ALLOCATE(mesh2%simfact,(mesh1%mesh_size))
 do ir=1,mesh1%mesh_size
   mesh2%rad(ir)    =mesh1%rad(ir)
   mesh2%radfact(ir)=mesh1%radfact(ir)
   mesh2%simfact(ir)=mesh1%simfact(ir)
 end do

end subroutine pawrad_copy
!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/pawrad_deducer0
!! NAME
!! pawrad_deducer0
!!
!! FUNCTION
!! Extrapolate r=0 value of a function from values near r=0
!! using a 3 points formula
!!
!! INPUTS
!!  funcsz=size of array func
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!
!! SIDE EFFECTS
!!  func(funcsz)=array containing values of function to extrapolate
!!
!! PARENTS
!!      m_orbmag,m_paw_atom,m_paw_denpot,m_paw_gaussfit,m_paw_init,m_paw_mkrho
!!      m_paw_nmr,m_paw_onsite,m_paw_pwaves_lmn,m_paw_slater,m_pawdij,m_pawpsp
!!      m_pawrad,m_pawxc
!!
!! CHILDREN
!!      poisson,simp_gen
!!
!! SOURCE

subroutine pawrad_deducer0(func,funcsz,radmesh)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: funcsz
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(inout) :: func(funcsz)

! *************************************************************************

 if (radmesh%mesh_type==1.or.radmesh%mesh_type==2.or.radmesh%mesh_type==4.or.radmesh%mesh_type==5) then
   func(1)=func(4)+3*(func(2)-func(3))
 else if (radmesh%mesh_type==3) then
   func(1)=func(4)+exp(two*radmesh%lstep)/(exp(radmesh%lstep)-one)*(func(2)-func(3))
 end if

end subroutine pawrad_deducer0
!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/pawrad_bcast
!! NAME
!! pawrad_bcast
!!
!! FUNCTION
!! Communicate pawrad data over a given MPI communicator
!!
!! INPUTS
!! comm_mpi= communicator used to broadcast data
!!
!! SIDE EFFECTS
!!  pawrad=<type pawrad_type>= a radial mesh datastructure for PAW
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!      poisson,simp_gen
!!
!! SOURCE

subroutine pawrad_bcast(pawrad,comm_mpi)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm_mpi
 type(pawrad_type),intent(inout) :: pawrad

!Local variables-------------------------------
!scalars
 integer :: ierr,indx,me,nn,isz1
 integer :: if_rad,if_radfact,if_simfact !flags used to communicate
 character(len=500) :: msg
!arrays
 integer,allocatable :: list_int(:)
 real(dp),allocatable :: list_dpr(:)

!*************************************************************************

 me=xmpi_comm_rank(comm_mpi)

!Initializations
 if_rad=0; if_radfact=0; if_simfact=0

!calculate the size of the reals
 if(me==0) then
   if (allocated(pawrad%rad)) then
     if_rad=1 !communicate rad
     isz1=size(pawrad%rad)
     if(isz1/=pawrad%mesh_size) then
       msg='rad: sz1 /= pawrad%mesh_size (1)'
       LIBPAW_BUG(msg)
     end if
   end if
   if (allocated(pawrad%radfact)) then
     if_radfact=1 !communicate radfact
     isz1=size(pawrad%radfact)
     if(isz1/=pawrad%mesh_size) then
       msg='radfact: sz1 /= pawrad%mesh_size (2)'
       LIBPAW_BUG(msg)
     end if
   end if
   if (allocated(pawrad%simfact)) then
     if_simfact=1 !communicate simfact
     isz1=size(pawrad%simfact)
     if(isz1/=pawrad%mesh_size) then
       msg='simfact: sz1 /= pawrad%mesh_size (3)'
       LIBPAW_BUG(msg)
     end if
   end if
 end if

!Brodcast the integers
 LIBPAW_ALLOCATE(list_int,(6))
 if(me==0) then
   list_int(1)=pawrad%int_meshsz
   list_int(2)=pawrad%mesh_size
   list_int(3)=pawrad%mesh_type
   list_int(4)=if_rad
   list_int(5)=if_radfact
   list_int(6)=if_simfact
 end if
 call xmpi_bcast(list_int,0,comm_mpi,ierr)
 if(me/=0) then
   pawrad%int_meshsz =list_int(1)
   pawrad%mesh_size =list_int(2)
   pawrad%mesh_type =list_int(3)
   if_rad=list_int(4)
   if_radfact=list_int(5)
   if_simfact=list_int(6)
 end if
 LIBPAW_DEALLOCATE(list_int)

!Broadcast the reals
 nn=4+pawrad%mesh_size*(if_rad+if_radfact+if_simfact)
 LIBPAW_ALLOCATE(list_dpr,(nn))
 if(me==0) then
   list_dpr(1)=pawrad%lstep
   list_dpr(2)=pawrad%rmax
   list_dpr(3)=pawrad%rstep
   list_dpr(4)=pawrad%stepint
   indx=5
   if (if_rad==1) then
     isz1=pawrad%mesh_size
     list_dpr(indx:indx+isz1-1)=pawrad%rad(1:isz1)
     indx=indx+isz1
   end if
   if (if_radfact==1) then
     isz1=pawrad%mesh_size
     list_dpr(indx:indx+isz1-1)=pawrad%radfact(1:isz1)
     indx=indx+isz1
   end if
   if (if_simfact==1) then
     isz1=pawrad%mesh_size
     list_dpr(indx:indx+isz1-1)=pawrad%simfact(1:isz1)
     indx=indx+isz1
   end if
 end if
 call xmpi_bcast(list_dpr,0,comm_mpi,ierr)
 if(me/=0) then
   pawrad%lstep=list_dpr(1)
   pawrad%rmax=list_dpr(2)
   pawrad%rstep=list_dpr(3)
   pawrad%stepint=list_dpr(4)
   indx=5
!  Deallocate all arrays:
   if (allocated(pawrad%rad)) then
     LIBPAW_DEALLOCATE(pawrad%rad)
   end if
   if (allocated(pawrad%radfact)) then
     LIBPAW_DEALLOCATE(pawrad%radfact)
   end if
   if (allocated(pawrad%simfact)) then
     LIBPAW_DEALLOCATE(pawrad%simfact)
   end if
!  Communicate if flag is set to 1:
   if(if_rad==1) then
     isz1=pawrad%mesh_size
     LIBPAW_ALLOCATE(pawrad%rad,(isz1))
     pawrad%rad(1:isz1)=list_dpr(indx:indx+isz1-1)
     indx=indx+isz1
   end if
   if(if_radfact==1) then
     isz1=pawrad%mesh_size
     LIBPAW_ALLOCATE(pawrad%radfact,(isz1))
     pawrad%radfact(1:isz1)=list_dpr(indx:indx+isz1-1)
     indx=indx+isz1
   end if
   if(if_simfact==1) then
     isz1=pawrad%mesh_size
     LIBPAW_ALLOCATE(pawrad%simfact,(isz1))
     pawrad%simfact(1:isz1)=list_dpr(indx:indx+isz1-1)
     indx=indx+isz1
   end if
 end if
 LIBPAW_DEALLOCATE(list_dpr)

end subroutine pawrad_bcast
!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/simp_gen
!! NAME
!! simp_gen
!!
!! FUNCTION
!! Do integral on a given (generalized) grid using Simpson rule
!!
!! INPUTS
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!  func(:)=integrand values
!!  r_for_intg=upper boundary for (future) integration over the radial grid
!!            (can be negative for an integration over the whole grid)
!!
!! OUTPUT
!!  intg=resulting integral by Simpson rule
!!
!! PARENTS
!!      m_epjdos,m_mlwfovlp,m_orbmag,m_paw_atom,m_paw_atomorb
!!      m_paw_correlations,m_paw_denpot,m_paw_dfptnl,m_paw_hr,m_paw_init
!!      m_paw_nmr,m_paw_onsite,m_paw_overlap,m_paw_slater,m_pawdij,m_pawpsp
!!      m_pawpwij,m_pawrad,m_pawxc,m_plowannier,m_positron,m_psps
!!
!! CHILDREN
!!      poisson,simp_gen
!!
!! NOTES
!!  Possible mesh types (radmesh%mesh_type)
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!   mesh_type=5 ( grid): rad(i)=AA*i/(n-i)
!!
!! SOURCE

subroutine simp_gen(intg,func,radmesh,r_for_intg)

#if defined HAVE_AVX_SAFE_MODE
!DEC$ NOOPTIMIZE
#endif

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: intg
 real(dp),intent(in),optional :: r_for_intg
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(in) :: func(:)

!Local variables-------------------------------
!scalars
 integer :: ii,int_meshsz,ir,ir_last,isim,nn
 real(dp) :: hh,resid,simp
 real(dp),allocatable :: simfact(:)
 character(len=500) :: msg

! *************************************************************************

 if (present(r_for_intg)) then
   if (r_for_intg>0.d0) then
     ir=min(pawrad_ifromr(radmesh,r_for_intg),radmesh%mesh_size)
     if (ir<radmesh%mesh_size) then
       if (abs(radmesh%rad(ir+1)-r_for_intg)<abs(radmesh%rad(ir)-r_for_intg)) ir=ir+1
     end if
     if (ir>1) then
       if (abs(radmesh%rad(ir-1)-r_for_intg)<abs(radmesh%rad(ir)-r_for_intg)) ir=ir-1
     end if
     int_meshsz=ir
   else
     int_meshsz=radmesh%mesh_size
   end if
   if (int_meshsz>radmesh%mesh_size.or.int_meshsz>size(func)) then
     write(msg,'(3(a,i4))')"int_meshsz= ",int_meshsz," > mesh_size=",radmesh%mesh_size,&
&                          ", size(func)=",size(func)
     LIBPAW_BUG(msg)
   end if
   isim=3; if (radmesh%mesh_type==3)isim=4
   LIBPAW_ALLOCATE(simfact,(radmesh%mesh_size))
   hh=radmesh%stepint/3.d0
   simfact(int_meshsz)=hh*radmesh%radfact(int_meshsz)
   simfact(1:isim-2)=zero
   ir_last=1
   do ir=int_meshsz,isim,-2
     simfact(ir-1)=4.d0*hh*radmesh%radfact(ir-1)
     simfact(ir-2)=2.d0*hh*radmesh%radfact(ir-2)
     ir_last=ir-2
   end do
   simfact(ir_last)=half*simfact(ir_last)
   if (int_meshsz<radmesh%mesh_size) simfact(int_meshsz+1:radmesh%mesh_size)=zero

   nn=int_meshsz
   simp=zero
   do ii=1,nn
     simp=simp+func(ii)*simfact(ii)
   end do
   LIBPAW_DEALLOCATE(simfact)

 else
   if (radmesh%int_meshsz>size(func)) then
     write(msg,'(2(a,i4))')"int_meshsz= ",int_meshsz," > size(func)=",size(func)
     LIBPAW_BUG(msg)
   end if
   nn=radmesh%int_meshsz
   simp=zero
   do ii=1,nn
     simp=simp+func(ii)*radmesh%simfact(ii)
   end do
 end if

 resid=zero
 if (radmesh%mesh_type==3) then
   resid=half*(func(2)+func(1))*(radmesh%rad(2)-radmesh%rad(1))
   if (mod(nn,2)==1) resid=resid+radmesh%stepint/3.d0*(1.25d0*func(2)*radmesh%radfact(2) &
&   +2.d0*func(3)*radmesh%radfact(3)-0.25d0*func(4)*radmesh%radfact(4))
 else if (mod(nn,2)==0) then
   resid=radmesh%stepint/3.d0*(1.25d0*func(1)*radmesh%radfact(1)+2.d0*func(2)*radmesh%radfact(2) &
&   -0.25d0*func(3)*radmesh%radfact(3))
 end if

 intg=simp+resid

end subroutine simp_gen
!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/nderiv_gen
!! NAME
!! nderiv_gen
!!
!! FUNCTION
!! Do corrected first (and -if requested- second) derivation on a given (generalized) grid.
!! This routine interfaces nderiv_lin (derivation on a linear grid).
!!
!! INPUTS
!!  func(:)=input function
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!
!! OUTPUT
!!  der(:)= 1st derivative of input function
!!  [der2(:)]= -- optional -- 2nd derivative of input function
!!
!! PARENTS
!!      m_paw_init,m_paw_onsite,m_paw_pwaves_lmn,m_pawdij,m_pawpsp,m_pawxc
!!      m_positron
!!
!! CHILDREN
!!      poisson,simp_gen
!!
!! NOTES
!!  Possible mesh types (radmesh%mesh_type)
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!
!! SOURCE

subroutine nderiv_gen(der,func,radmesh,der2)

!Arguments ------------------------------------
!scalars
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(in) :: func(:)
 real(dp),intent(out) :: der(:)
 real(dp),optional,intent(out) :: der2(:)

!Local variables-------------------------------
!scalars
 integer :: msz
 logical :: compute_2der
 character(len=500) :: msg

! *************************************************************************

 msz=size(func)
 if (size(der)/=msz.or.msz>radmesh%mesh_size) then
   msg='wrong sizes for in/out arrays!'
   LIBPAW_BUG(msg)
 end if

 compute_2der=(present(der2))

 if (radmesh%mesh_type==1) then

   call nderiv_lin(radmesh%rstep,func,der,msz,1)
   if (compute_2der) then
     call nderiv_lin(radmesh%rstep,func,der2,msz,2)
   end if

 else if (radmesh%mesh_type==2) then

   call nderiv_lin(radmesh%lstep,func,der,msz,1)
   der(1:msz)=der(1:msz)/radmesh%radfact(1:msz)
   if (compute_2der)then
     call nderiv_lin(radmesh%lstep,func,der2,msz,2)
     der2(1:msz)=(der2(1:msz)/radmesh%radfact(1:msz)-der(1:msz))/radmesh%radfact(1:msz)
   end if

 else if (radmesh%mesh_type==3) then

   call nderiv_lin(radmesh%lstep,func(2:msz),der(2:msz),msz-1,1)
   der(2:msz)=der(2:msz)/radmesh%radfact(2:msz)
   call pawrad_deducer0(der,msz,radmesh)
   if (compute_2der)then
     call nderiv_lin(radmesh%lstep,func(2:msz),der2(2:msz),msz-1,2)
     der2(2:msz)=(der2(2:msz)/radmesh%radfact(2:msz)-der(2:msz))/radmesh%radfact(2:msz)
     call pawrad_deducer0(der2,msz,radmesh)
   end if

 else if (radmesh%mesh_type==4) then

   call nderiv_lin(radmesh%lstep,func,der,msz,1)
   der(1:msz)=der(1:msz)/radmesh%radfact(1:msz)
   if (compute_2der)then
     call nderiv_lin(radmesh%lstep,func,der2,msz,2)
     der2(1:msz)=der2(1:msz)/radmesh%radfact(1:msz)**2-der(1:msz)/radmesh%rstep
   end if

 else if (radmesh%mesh_type==5) then

   call nderiv_lin(one,func,der,msz,1)
   der(1:msz)=der(1:msz)/(radmesh%radfact(1:msz)*radmesh%rstep)
   if (compute_2der)then
     call nderiv_lin(one,func,der2,msz,2)
     der2(1:msz)=der2(1:msz)/(radmesh%radfact(1:msz)*radmesh%rstep)**2-two*der(1:msz)/&
&                            (radmesh%rstep+radmesh%rad(1:msz))
   end if
 end if

end subroutine nderiv_gen
!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/nderiv_lin
!! NAME
!! nderiv_lin
!!
!! FUNCTION
!! Do corrected first (and -if requested- second) derivation on a given LINEAR grid.
!!
!! INPUTS
!!  hh= radial step
!!  ndim= radial mesh size
!!  yy(ndim)= input function
!!  norder= order of derivation (1 or 2)
!!
!! OUTPUT
!!  zz(ndim)= first or second derivative of y
!!
!! PARENTS
!!      m_pawrad
!!
!! CHILDREN
!!      poisson,simp_gen
!!
!! SOURCE

subroutine nderiv_lin(hh,yy,zz,ndim,norder)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ndim,norder
 real(dp),intent(in) :: hh
!arrays
 real(dp),intent(in) :: yy(ndim)
 real(dp),intent(out) :: zz(ndim)

!Local variables ---------------------------------------
!scalars
 integer :: ier,ii
 real(dp) :: aa,bb,cc,h1,y1

! *************************************************************************

!Initialization (common to 1st and 2nd derivative)
 h1=one/(12.d0*hh)
 y1=yy(ndim-4)

!FIRST DERIVATIVE
!================
 if (norder==1) then

!  Prepare differentiation loop
   bb=h1*(-25.d0*yy(1)+48.d0*yy(2)-36.d0*yy(3)+16.d0*yy(4)-3.d0*yy(5))
   cc=h1*(-3.d0*yy(1)-10.d0*yy(2)+18.d0*yy(3)-6.d0*yy(4)+yy(5))
!  Start differentiation loop
   do ii=5,ndim
     aa=bb;bb=cc
     cc=h1*(yy(ii-4)-yy(ii)+8.d0*(yy(ii-1)-yy(ii-3)))
     zz(ii-4)=aa
   end do
!  Normal exit
   ier=0
   aa=h1*(-y1+6.d0*yy(ndim-3)-18.d0*yy(ndim-2)+10.d0*yy(ndim-1)+3.d0*yy(ndim))
   zz(ndim)=h1*(3.d0*y1-16.d0*yy(ndim-3)+36.d0*yy(ndim-2) -48.d0*yy(ndim-1)+25.d0*yy(ndim))
   zz(ndim-1)=aa
   zz(ndim-2)=cc
   zz(ndim-3)=bb

!  SECOND DERIVATIVE
!  =================
 else
   h1=h1/hh
!  Prepare differentiation loop
   bb=h1*(35.d0*yy(1)-104.d0*yy(2)+114.d0*yy(3)-56.d0*yy(4)+11.d0*yy(5))
   cc=h1*(11.d0*yy(1)-20.d0*yy(2)+6.d0*yy(3)+4.d0*yy(4)-yy(5))
!  Start differentiation loop
   do ii=5,ndim
     aa=bb;bb=cc
     cc=h1*(-yy(ii-4)-yy(ii)+16.d0*(yy(ii-1)+yy(ii-3))-30.d0*yy(ii-2))
     zz(ii-4)=aa
   end do
!  Normal exit
   ier=0
   aa=h1*(-y1+4.d0*yy(ndim-3)+6.d0*yy(ndim-2)-20.d0*yy(ndim-1)+11.d0*yy(ndim))
   zz(ndim)=h1*(11.d0*y1-56.d0*yy(ndim-3)+114.d0*yy(ndim-2) -104.d0*yy(ndim-1)+35.d0*yy(ndim))
   zz(ndim-1)=aa
   zz(ndim-2)=cc
   zz(ndim-3)=bb

 end if !norder

end subroutine nderiv_lin
!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/bound_deriv
!! NAME
!! bound_deriv
!!
!! FUNCTION
!! Computes derivatives of a function a boundaries of interval (first and last derivative)
!!
!! INPUTS
!!  func(n)= array containing function
!!  mesh <type(pawrad_type)>= radial mesh and related data
!!  nn= size of intervall
!!
!! OUTPUT
!!  yp1,ypn= derivatives of func at r(1) and r(n)
!!
!! PARENTS
!!      m_outscfcv,m_paw_atom,m_pawpsp,m_pawxmlps
!!
!! CHILDREN
!!      poisson,simp_gen
!!
!! SOURCE

 subroutine bound_deriv(func,mesh,nn,yp1,ypn)

!Arguments----------------------
 integer, intent(in) :: nn
 real(dp), intent(in) :: func(nn)
 real(dp), intent(out) :: yp1,ypn
 type(pawrad_type),intent(in) :: mesh

!*************************************************************************

 if (mesh%radfact(1)>zero) then
   yp1=1._dp/12._dp/mesh%stepint/mesh%radfact(1) &
&   *(-25._dp*func(1)+48._dp*func(2)-36._dp*func(3)+16._dp*func(4)-3._dp*func(5))
 else
   yp1=(func(2)-func(1))/(mesh%rad(2)-mesh%rad(1))
 end if
 ypn=1._dp/12._dp/mesh%stepint &
& *( 3._dp*func(nn-4)-16._dp*func(nn-3)+36._dp*func(nn-2)-48._dp*func(nn-1) &
& +25._dp*func(nn))/mesh%radfact(nn)

end subroutine bound_deriv

!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/poisson
!! NAME
!! poisson
!!
!! FUNCTION
!!  Solve poisson equation for angularly dependent charge
!!  distribution of angular momentum l
!!  Densities and potentials are given on a (generalized) radial grid
!!
!! INPUTS
!!  den(:)= electron density * (4*pi*r**2) appropriate for l
!!  ll= l quantum number
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!  [screened_sr_separation]= --optional-- separation for screened short-range kernel
!!                            no screening by default
!!
!! OUTPUT
!!  [qq]= --optional-- lth moment of the charge; not compatible with screened Coulomb interaction
!!  rv(:)= electrostatic potential * r in (Hartree*Bohr) units
!!          where v(r)=\frac{1}{2l+1}(\frac{int[(r''^(l+2))g(r'')dr'']} {r^(l+1)}
!!                                   +(r^l) int[r''^(1-l)g(r'')dr''])
!!
!! PARENTS
!!      m_paw_atom,m_paw_correlations,m_paw_denpot,m_paw_init,m_pawpsp,m_pawrad
!!
!! CHILDREN
!!      poisson,simp_gen
!!
!! SOURCE

subroutine poisson(den,ll,radmesh,rv,screened_sr_separation,qq)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll
 real(dp),intent(out),optional :: qq
 real(dp),intent(in),optional :: screened_sr_separation
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(in) :: den(:)
 real(dp),intent(out) :: rv(:)

!Local variables ---------------------------------------
!scalars
 integer :: ir,jr,mesh_size,mm,nn
 logical :: use_numerov,use_screened
 real(dp) :: angm,hh,intg,qq_,ri,rj,rr,sr_omega
!arrays
 real(dp) :: ee(4)
 real(dp),allocatable :: aa(:),bb(:),cc(:),dd(:),radl(:),radl1(:)

! ***************************************************************************

 mesh_size=size(den)
 if (size(rv)/=mesh_size.or.mesh_size>radmesh%mesh_size) then
   LIBPAW_BUG('wrong sizes!')
 end if

 use_numerov=(radmesh%mesh_type==1)

 use_screened=.false.;sr_omega=zero
 if (present(screened_sr_separation)) then
   sr_omega=screened_sr_separation
   use_screened=(sr_omega>tol8)
!  Numerov method not coded for screened Coulomb interaction
   if (use_numerov.and.use_screened) use_numerov=.false.
 end if

!==============================================
!UNIFORM GRID - NUMEROV ALGORITHM
!==============================================
 if (use_numerov) then
   hh=radmesh%rstep
   nn=radmesh%int_meshsz-1
   LIBPAW_ALLOCATE(aa,(nn))
   LIBPAW_ALLOCATE(bb,(nn))
   LIBPAW_ALLOCATE(cc,(nn+1))
   do ir=1,nn
     aa(ir)=two*hh*den(ir+1)/(ir)
     bb(ir)=den(ir+1)*((ir*hh)**ll)
   end do
   cc(1)=zero
   cc(2:nn+1)=bb(1:nn)
   call simp_gen(qq_,cc,radmesh)
   qq_=qq_/dble(2*ll+1)
   rv(1)=aa(1)+0.1_dp*aa(2)
   do ir=2,nn-1
     rv(ir)=aa(ir)+0.1_dp*(aa(ir+1)+aa(ir-1))
   end do
   rv(nn)=aa(nn)+0.1_dp*aa(nn-1)
   angm=dble(ll*(ll+1))
   rr=(nn+1)*hh
   rv(nn)=rv(nn)+(2.4_dp-0.2_dp*angm/((nn+1)**2))*qq_/(rr**ll)
   do ir=1,nn
     aa(ir)=angm/(ir*ir)
     bb(ir)=2.4_dp+aa(ir)
   end do
   do ir=1,nn-1
     cc(ir)=-1.2_dp+0.1_dp*aa(ir+1)
   end do
   do ir=nn,2,-1
     aa(ir)=-1.2_dp+0.1_dp*aa(ir-1)
   end do
   if (nn.eq.1) then
     rv(2)=rv(1)/bb(1)
     rv(1)=zero
   else
     do ir=2,nn
       bb(ir)=bb(ir)-aa(ir)*cc(ir-1)/bb(ir-1)
     end do
     rv(1)=rv(1)/bb(1)
     do ir=2,nn
       rv(ir)=(rv(ir)-aa(ir)*rv(ir-1))/bb(ir)
     end do
     do ir=nn-1,1,-1
       rv(ir)=rv(ir)-cc(ir)*rv(ir+1)/bb(ir)
     end do
     do ir=nn+1,2,-1
       rv(ir)=rv(ir-1)
     end do
     rv(1)=zero
     if (nn+1<mesh_size) rv(nn+2:mesh_size)=zero
   end if
   rv(:)=half*rv(:)
   if (present(qq)) qq=qq_
   LIBPAW_DEALLOCATE(aa)
   LIBPAW_DEALLOCATE(bb)
   LIBPAW_DEALLOCATE(cc)

!  ==============================================
!  ANY OTHER GRID - SIMPSON ALGORITHM
!  ==============================================
 else
   nn=mesh_size;hh=third*radmesh%stepint
   do while (abs(den(nn))<tol16.and.nn>radmesh%int_meshsz)
     nn=nn-1
   end do
   mm=nn;if (radmesh%mesh_type==3) mm=mm-1
   LIBPAW_ALLOCATE(aa,(nn))
   LIBPAW_ALLOCATE(bb,(nn))
   LIBPAW_ALLOCATE(cc,(nn))
   LIBPAW_ALLOCATE(dd,(nn))

   if (.not.use_screened) then

!    Standard Coulomb integral
     LIBPAW_ALLOCATE(radl,(nn))
     LIBPAW_ALLOCATE(radl1,(nn))
     do jr=nn,2,-1
       ir=nn-jr+1
       radl(jr) =radmesh%rad(jr)**ll
       radl1(jr)=radmesh%rad(jr)*radl(jr)
       aa(ir)=den(jr)*radmesh%radfact(jr)*radl(jr)
       bb(ir)=den(jr)*radmesh%radfact(jr)/radl1(jr)
     end do
     radl(1)=zero;radl1(1)=zero
     ee(2)=aa(nn-1);ee(3)=aa(nn-2);ee(4)=aa(nn-3)
     call pawrad_deducer0(ee,4,radmesh)
     aa(nn)=ee(1)
     ee(2)=bb(nn-1);ee(3)=bb(nn-2);ee(4)=bb(nn-3)
     call pawrad_deducer0(ee,4,radmesh)
     bb(nn)=ee(1)
     cc(1)=zero;dd(1)=zero
     do ir=3,mm,2
       cc(ir)  =cc(ir-2)+hh*(aa(ir-2)+four*aa(ir-1)+aa(ir))
       cc(ir-1)=cc(ir-2)+hh*(1.25_dp*aa(ir-2)+two*aa(ir-1)-quarter*aa(ir))
       dd(ir)  =dd(ir-2)+hh*(bb(ir-2)+four*bb(ir-1)+bb(ir))
       dd(ir-1)=dd(ir-2)+hh*(1.25_dp*bb(ir-2)+two*bb(ir-1)-quarter*bb(ir))
     end do
     if (mod(mm,2)==0) then
!      cc(mm)=cc(mm-2)+hh*(aa(mm-2)+four*aa(mm-1)+aa(mm))
!      dd(mm)=dd(mm-2)+hh*(bb(mm-2)+four*bb(mm-1)+bb(mm))
       cc(mm)=cc(mm-1)+hh*(1.25_dp*aa(mm-2)+two*aa(mm-1)-quarter*aa(mm))
       dd(mm)=dd(mm-1)+hh*(1.25_dp*bb(mm-2)+two*bb(mm-1)-quarter*bb(mm))
     end if
     if (mm<nn) then
       cc(nn)=cc(mm)+half*(aa(mm)+aa(nn))*(radmesh%rad(1+nn-mm)-radmesh%rad(1))
       dd(nn)=dd(mm)+half*(bb(mm)+bb(nn))*(radmesh%rad(1+nn-mm)-radmesh%rad(1))
     end if
     rv(1)=zero
     do ir=2,nn
       jr=nn-ir+1
       rv(ir)=(dd(jr)*radl1(ir)+(cc(nn)-cc(jr))/radl(ir))/(two*ll+one)
     end do
     if (nn<mesh_size) rv(nn+1:mesh_size)=rv(nn)
     if (present(qq)) qq=cc(nn)/(two*ll+one)
     LIBPAW_DEALLOCATE(radl)
     LIBPAW_DEALLOCATE(radl1)
   else

!    Short-range screened Coulomb integral
     rv(1)=zero
     do ir=2,nn
       ri=sr_omega*radmesh%rad(ir)
       do jr=2,nn
         rj=sr_omega*radmesh%rad(jr)
         aa(jr)=den(jr)*radmesh%radfact(jr)*sr_omega*screened_coul_kernel(ll,ri,rj)
       end do
       call pawrad_deducer0(aa(1:4),4,radmesh)
       intg=zero
!      Compute the integral in 2 parts because the function is not differentiable at r(ir)
!      Integral from zero to r(ir)
       if (ir>2) then
         do jr=ir-2,1,-2
           intg=intg+hh*(aa(jr)+four*aa(jr+1)+aa(jr+2))
         end do
         if (mod(ir,2)==0) intg=intg+hh*(-quarter*aa(1)+two*aa(2)+1.25_dp*aa(3))
       else if (ir==2) then
        if (nn >2) intg=intg+hh*(+1.25_dp*aa(1)+two*aa(2)-quarter*aa(3))
        if (nn==2) intg=intg+(aa(1)+aa(2))*hh*half
       end if
       if (mm<nn) intg=intg+half*(aa(1+nn-mm)+aa(nn))*(radmesh%rad(1+nn-mm)-radmesh%rad(1))
!      Integral from r(ir) to rmax
       if (ir<nn-2) then
         do jr=ir+2,nn,2
           intg=intg+hh*(aa(jr-2)+four*aa(jr-1)+aa(jr))
         end do
         if (mod(nn-ir,2)==1) then
           if (ir<=nn-4) then
             !Use a high-order formula because points are spaced
             intg=intg+hh*(251._dp*aa(nn)+646._dp*aa(nn-1)-264._dp*aa(nn-2) &
&                         +106._dp*aa(nn-3)-19._dp*aa(nn-4))/240._dp
           else
             intg=intg+hh*(1.25_dp*aa(nn)+two*aa(nn-1)-quarter*aa(nn-2))
           endif
         end if
       else if (ir==nn-1) then
         intg=intg+(aa(nn-1)+aa(nn))*hh*half
       end if
       rv(ir)=intg/(two*ll+one)*radmesh%rad(ir)
     end do
     if (nn<mesh_size) rv(nn+1:mesh_size)=rv(nn)
     if (present(qq)) qq=zero ! Not relevant

   end if
   LIBPAW_DEALLOCATE(aa)
   LIBPAW_DEALLOCATE(bb)
   LIBPAW_DEALLOCATE(cc)
   LIBPAW_DEALLOCATE(dd)

 end if

end subroutine poisson
!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/pawrad_ifromr
!! NAME
!! pawrad_ifromr
!!
!! FUNCTION
!! Retreive Index FROM a given R value in a radial grid
!! Grid can be regular or logarithimc
!!
!! INPUTS
!!  rr=given input r value
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!
!! OUTPUT
!!  pawrad_ifromr=index of rr in radial grid
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!  Possible mesh types (radmesh%mesh_type)
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!   mesh_type=5 ( grid): rad(i)=AA*i/(n-i)
!!
!! SOURCE

function pawrad_ifromr(radmesh,rr)

!Arguments ------------------------------------
!scalars
 integer :: pawrad_ifromr
 real(dp),intent(in) :: rr
 type(pawrad_type),intent(in) :: radmesh

!Local variables-------------------------------
 character(len=500) :: msg

! *************************************************************************

 if (radmesh%mesh_type==1) then
   pawrad_ifromr=int(tol8+rr/radmesh%rstep)+1
 else if (radmesh%mesh_type==2) then
   pawrad_ifromr=int(tol8+log(1.d0+rr/radmesh%rstep)/radmesh%lstep)+1
 else if (radmesh%mesh_type==3) then
   if (rr<radmesh%rstep) then
     pawrad_ifromr=1
   else
     pawrad_ifromr=int(tol8+log(rr/radmesh%rstep)/radmesh%lstep)+2
   end if
 else if (radmesh%mesh_type==4) then
   pawrad_ifromr=int(tol8+(1.d0-exp(-rr/radmesh%rstep))/radmesh%lstep)+1
 else if (radmesh%mesh_type==5) then
   pawrad_ifromr=int(tol8+(radmesh%lstep*rr)/(radmesh%rstep+rr))+1
 else
!  Other values of mesh_type are not allowed (see psp7in.F90)
   write(msg,'(a,i0)')" Unknown value of %mesh_type ",radmesh%mesh_type
   LIBPAW_ERROR(msg)
 end if

end function pawrad_ifromr
!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/screened_coul_kernel
!! NAME
!! screened_coul_kernel
!!
!! FUNCTION
!!  Evaluates the kernel function used to compute the integral
!!  of a screened Coulomb potential (with erfc function)
!!  See Angyan, Gerber, Marsman, J. Phys. A: Math. Gen. 39, 8613 (2006) [[cite:Angyan2006]]
!!
!! INPUTS
!!  [formula]=optional; used to force formula (1=full; 2=dev. near r2=0; 3=dev near r1*r2=0)
!!  order= order of the function (typically l quantum number)
!!  r1,r2=input arguments
!!
!! OUTPUT
!!  screened_coul_kernel=output radial function
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function screened_coul_kernel(order,r1,r2,formula)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: order
 integer,optional :: formula
 real(dp),intent(in) :: r1,r2
 real(dp) :: screened_coul_kernel
!Local variables-------------------------------
 integer :: formula_
 real(dp) :: difexp,erfcp,erfcm,f0,f1,f2,f3,f4,f5,f6,hh,sqrtpi,sumexp,xx,xy,yy
 character(len=100) :: msg

!******************************************************************

 if (order>6) then
   msg='PAW screened exchange not coded for l>2!'
   LIBPAW_ERROR(msg)
 end if

!Use max and min of arguments
 xx=max(r1,r2) ; yy=min(r1,r2)

!Unscreened Coulomb interaction for very small (or negative) arguments
 if (xx<tol8) then
  screened_coul_kernel=yy**order/xx**(order+1) ; return
 end if

!Choice of formula
!Empirical criterion, inspired by J. Phys. A 39, pp8624 [[cite:Angyan2006]] and adjusted
 formula_=1;if (xx<0.25_dp.and.yy<0.25_dp) formula_=2
 if (present(formula)) formula_=formula
 select case (formula_)

!Full formula
!J. Phys. A 39, 8613 (2006) - Eq. (21), (22), (24) [[cite:Angyan2006]]
!   Note: typo in the paper: in eq (22), the sum begins at p=0
!------------------------------------------------------------------
 case(1)
   xy=xx*yy ; sqrtpi=sqrt(pi)
   sumexp=exp(-(xx+yy)**2)+exp(-(xx-yy)**2)
   difexp=exp(-(xx+yy)**2)-exp(-(xx-yy)**2)
   erfcp=paw_derfc(xx+yy) ; erfcm=paw_derfc(xx-yy)
   if (order==0) then
     f0=-difexp/(two*sqrtpi*xy)
     hh=((xx+yy)*erfcp-(xx-yy)*erfcm)/(two*xy)
     screened_coul_kernel = hh + f0
   else if(order==1) then
     f0=-difexp/(two*sqrtpi*xy)
     f1=(difexp+two*xy*sumexp)/(sqrtpi*(two*xy)**2)
     hh=((xx**3+yy**3)*erfcp-(xx**3-yy**3)*erfcm)/(two*xy**2)
     screened_coul_kernel = hh + f1 + f0*(xx**2+yy**2)/xy
   else if (order==2) then
     f0=-difexp/(two*sqrtpi*xy)
     f1=(difexp+two*xy*sumexp)/(sqrtpi*(two*xy)**2)
     f2= -((four*xy**2+three)*difexp+6._dp*xy*sumexp)/(sqrtpi*(two*xy)**3)
     hh=((xx**5+yy**5)*erfcp-(xx**5-yy**5)*erfcm)/(two*xy**3)
     screened_coul_kernel = hh + f2 + f1*(xx**2+yy**2)/xy + f0*(xx**4+yy**4)/xy**2
   else if (order==3) then
     f0=-difexp/(two*sqrtpi*xy)
     f1=(difexp+two*xy*sumexp)/(sqrtpi*(two*xy)**2)
     f2= -((four*xy**2+three)*difexp+6._dp*xy*sumexp)/(sqrtpi*(two*xy)**3)
     f3=((24._dp*xy**2+15._dp)*difexp+(8._dp*xy**3+30._dp*xy)*sumexp)/(sqrtpi*(two*xy)**4)
     hh=((xx**7+yy**7)*erfcp-(xx**7-yy**7)*erfcm)/(two*xy**4)
     screened_coul_kernel = hh + f3 + f2*(xx**2+yy**2)/xy    + f1*(xx**4+yy**4)/xy**2 &
&                                   + f0*(xx**6+yy**6)/xy**3
   else if (order==4) then
     f0=-difexp/(two*sqrtpi*xy)
     f1=(difexp+two*xy*sumexp)/(sqrtpi*(two*xy)**2)
     f2= -((four*xy**2+three)*difexp+6._dp*xy*sumexp)/(sqrtpi*(two*xy)**3)
     f3=((24._dp*xy**2+15._dp)*difexp+(8._dp*xy**3+30._dp*xy)*sumexp)/(sqrtpi*(two*xy)**4)
     f4=-((16._dp*xy**4+180._dp*xy**2+105._dp)*difexp &
&        +(80._dp*xy**3+210._dp*xy)*sumexp)/(sqrtpi*(two*xy)**5)
     hh=((xx**9+yy**9)*erfcp-(xx**9-yy**9)*erfcm)/(two*xy**5)
     screened_coul_kernel = hh + f4 + f3*(xx**2+yy**2)/xy    + f2*(xx**4+yy**4)/xy**2 &
&                                   + f1*(xx**6+yy**6)/xy**3 + f0*(xx**8+yy**8)/xy**4
   else if (order==5) then
     f0=-difexp/(two*sqrtpi*xy)
     f1=(difexp+two*xy*sumexp)/(sqrtpi*(two*xy)**2)
     f2= -((four*xy**2+three)*difexp+6._dp*xy*sumexp)/(sqrtpi*(two*xy)**3)
     f3=((24._dp*xy**2+15._dp)*difexp+(8._dp*xy**3+30._dp*xy)*sumexp)/(sqrtpi*(two*xy)**4)
     f4=-((16._dp*xy**4+180._dp*xy**2+105._dp)*difexp &
&        +(80._dp*xy**3+210._dp*xy)*sumexp)/(sqrtpi*(two*xy)**5)
     f5=((240._dp*xy**4+1680._dp*xy**2+945._dp)*difexp &
&       +(32._dp*xy**5+840._dp*xy**3+1890._dp*xy)*sumexp)/(sqrtpi*(two*xy)**6)
     hh=((xx**11+yy**11)*erfcp-(xx**11-yy**11)*erfcm)/(two*xy**6)
     screened_coul_kernel = hh + f5 + f4*(xx**2+yy**2)  /xy    + f3*(xx**4+yy**4)/xy**2 &
&                                   + f2*(xx**6+yy**6)  /xy**3 + f1*(xx**8+yy**8)/xy**4 &
&                                   + f0*(xx**10+yy**10)/xy**5
   else if (order==6) then
     f0=-difexp/(two*sqrtpi*xy)
     f1=(difexp+two*xy*sumexp)/(sqrtpi*(two*xy)**2)
     f2= -((four*xy**2+three)*difexp+6._dp*xy*sumexp)/(sqrtpi*(two*xy)**3)
     f3=((24._dp*xy**2+15._dp)*difexp+(8._dp*xy**3+30._dp*xy)*sumexp)/(sqrtpi*(two*xy)**4)
     f4=-((16._dp*xy**4+180._dp*xy**2+105._dp)*difexp &
&        +(80._dp*xy**3+210._dp*xy)*sumexp)/(sqrtpi*(two*xy)**5)
     f5=((240._dp*xy**4+1680._dp*xy**2+945._dp)*difexp &
&       +(32._dp*xy**5+840._dp*xy**3+1890._dp*xy)*sumexp)/(sqrtpi*(two*xy)**6)
     f6=-((64._dp*xy**6+3360._dp*xy**4+18900._dp*xy**2+10395._dp)*difexp+ &
&         (672._dp*xy**5+10080._dp*xy**3+20790._dp*xy)*sumexp)/(sqrtpi*(two*xy)**7)
     hh=((xx**13+yy**13)*erfcp-(xx**13-yy**13)*erfcm)/(two*xy**7)
     screened_coul_kernel = hh + f6 + f5*(xx**2+yy**2)  /xy    + f4*(xx**4+yy**4)  /xy**2 &
&                                   + f3*(xx**6+yy**6)  /xy**3 + f1*(xx**8+yy**8)  /xy**4 &
&                                   + f1*(xx**10+yy**10)/xy**5 + f0*(xx**12+yy**12)/xy**6
   end if

!Development for yy->0
!J. Phys. A 39, 8613 (2006) - Eq. (28), (29), (30) [[cite:Angyan2006]]
!------------------------------------------------------------------
 case(2)
   sqrtpi=sqrt(pi)
   if (order==0) then
     screened_coul_kernel = paw_derfc(xx)/xx + exp(-xx**2)/sqrtpi * &
&            (two/three              *yy**2 &
&           +(two*xx**2-three)/15._dp*yy**4)
   else if(order==1) then
     screened_coul_kernel = paw_derfc(xx)*yy/xx**2 + exp(-xx**2)/sqrtpi * &
&            (two/xx                         *yy &
&            +four*xx/5._dp                  *yy**3 &
&            +two*xx*(two*xx**2-5._dp)/35._dp*yy**5)
   else if (order==2) then
     screened_coul_kernel = paw_derfc(xx)*yy**2/xx**3 + exp(-xx**2)/sqrtpi * &
&            (two*(two*xx**2+three)/(three*xx**2) *yy**2 &
&            +8._dp*xx**2/21._dp                  *yy**4 &
&            +four*xx**2*(two*xx**2-7._dp)/189._dp*yy**6)
   else if (order==3) then
     screened_coul_kernel = paw_derfc(xx)*yy**3/xx**4 + exp(-xx**2)/sqrtpi * &
&            (two*(four*xx**4+10._dp*xx**2+15._dp)/(15._dp*xx**3)*yy**3 &
&            +16._dp*xx**3/135._dp                               *yy**5 &
&            +8._dp*xx**3*(two*xx**2-9._dp)/1485._dp             *yy**7)
   else if (order==4) then
     screened_coul_kernel = paw_derfc(xx)*yy**4/xx**5 + exp(-xx**2)/sqrtpi * &
&            (two*(8._dp*xx**6+28._dp*xx**4+70._dp*xx**2+105._dp)/(105._dp*xx**4)*yy**4 &
&            +32._dp*xx**4/1155._dp                                              *yy**6 &
&            +16._dp*xx**4*(two*xx**2-11._dp)/15015._dp                          *yy**8)
   else if (order==5) then
     screened_coul_kernel = paw_derfc(xx)*yy**5/xx**6 + exp(-xx**2)/sqrtpi * &
&            (two*(16._dp*xx**8+72._dp*xx**6+252._dp*xx**4+630._dp*xx**2+945._dp)/(945._dp*xx**5)*yy**5 &
&            +64._dp*xx**5/12285._dp                                                             *yy**7 &
&            +32._dp*xx**5*(two*xx**2-13._dp)/184275._dp                                         *yy**9)
   else if (order==6) then
     screened_coul_kernel = paw_derfc(xx)*yy**6/xx**7 + exp(-xx**2)/sqrtpi * &
&            (two*(32._dp*xx**10+176._dp*xx**8+792._dp*xx**6+2772._dp*xx**4+6930._dp*xx**2+10395._dp)/(10395._dp*xx**6)*yy**6 &
&            +128._dp*xx**6/155925._dp                                                                                 *yy**8 &
&            +64._dp*xx**6*(two*xx**2-15._dp)/2650725._dp                                                              *yy**10)
   end if

!Development for xx*yy->0
!J. Phys. A 39, 8613 (2006) - Eq. (24), (26) [[cite:Angyan2006]]
!   Note: typo in the paper: (2n+3)! should be (2n+3)!!
!------------------------------------------------------------------
 case(3)
   xy=xx*yy ; sqrtpi=sqrt(pi)
   sumexp=exp(-(xx**2+yy**2))
   erfcp=paw_derfc(xx+yy) ; erfcm=paw_derfc(xx-yy)
   if (order==0) then
     f0 = two          *(three +two*xy**2)/(three      *sqrtpi)*sumexp
     hh=((xx+yy)*erfcp-(xx-yy)*erfcm)/(two*xy)
     screened_coul_kernel = hh + f0
   else if(order==1) then
     f0 = two          *(three +two*xy**2)/(three      *sqrtpi)*sumexp
     f1 = four   *xy   *(5._dp +two*xy**2)/(15._dp     *sqrtpi)*sumexp
     hh=((xx**3+yy**3)*erfcp-(xx**3-yy**3)*erfcm)/(two*xy**2)
     screened_coul_kernel = hh + f1 + f0*(xx**2+yy**2)/xy
   else if (order==2) then
     f0 = two          *(three +two*xy**2)/(three      *sqrtpi)*sumexp
     f1 = four   *xy   *(5._dp +two*xy**2)/(15._dp     *sqrtpi)*sumexp
     f2 = 8._dp  *xy**2*(7._dp +two*xy**2)/(105._dp    *sqrtpi)*sumexp
     hh=((xx**5+yy**5)*erfcp-(xx**5-yy**5)*erfcm)/(two*xy**3)
     screened_coul_kernel = hh + f2 + f1*(xx**2+yy**2)/xy + f0*(xx**4+yy**4)/xy**2
   else if (order==3) then
     f0 = two          *(three +two*xy**2)/(three      *sqrtpi)*sumexp
     f1 = four   *xy   *(5._dp +two*xy**2)/(15._dp     *sqrtpi)*sumexp
     f2 = 8._dp  *xy**2*(7._dp +two*xy**2)/(105._dp    *sqrtpi)*sumexp
     f3 = 16._dp *xy**2*(9._dp +two*xy**2)/(945._dp    *sqrtpi)*sumexp
     hh=((xx**7+yy**7)*erfcp-(xx**7-yy**7)*erfcm)/(two*xy**4)
     screened_coul_kernel = hh + f3 + f2*(xx**2+yy**2)/xy    + f1*(xx**4+yy**4)/xy**2 &
&                                   + f0*(xx**6+yy**6)/xy**3
   else if (order==4) then
     f0 = two          *(three +two*xy**2)/(three      *sqrtpi)*sumexp
     f1 = four   *xy   *(5._dp +two*xy**2)/(15._dp     *sqrtpi)*sumexp
     f2 = 8._dp  *xy**2*(7._dp +two*xy**2)/(105._dp    *sqrtpi)*sumexp
     f3 = 16._dp *xy**2*(9._dp +two*xy**2)/(945._dp    *sqrtpi)*sumexp
     f4 = 32._dp *xy**2*(11._dp+two*xy**2)/(10395._dp  *sqrtpi)*sumexp
     hh=((xx**9+yy**9)*erfcp-(xx**9-yy**9)*erfcm)/(two*xy**5)
     screened_coul_kernel = hh + f4 + f3*(xx**2+yy**2)/xy    + f2*(xx**4+yy**4)/xy**2 &
&                                   + f1*(xx**6+yy**6)/xy**3 + f0*(xx**8+yy**8)/xy**4
   else if (order==5) then
     f0 = two          *(three +two*xy**2)/(three      *sqrtpi)*sumexp
     f1 = four   *xy   *(5._dp +two*xy**2)/(15._dp     *sqrtpi)*sumexp
     f2 = 8._dp  *xy**2*(7._dp +two*xy**2)/(105._dp    *sqrtpi)*sumexp
     f3 = 16._dp *xy**2*(9._dp +two*xy**2)/(945._dp    *sqrtpi)*sumexp
     f4 = 32._dp *xy**2*(11._dp+two*xy**2)/(10395._dp  *sqrtpi)*sumexp
     f5 = 64._dp *xy**2*(13._dp+two*xy**2)/(135135._dp *sqrtpi)*sumexp
     hh=((xx**11+yy**11)*erfcp-(xx**11-yy**11)*erfcm)/(two*xy**6)
     screened_coul_kernel = hh + f5 + f4*(xx**2+yy**2)  /xy    + f3*(xx**4+yy**4)/xy**2 &
&                                   + f2*(xx**6+yy**6)  /xy**3 + f1*(xx**8+yy**8)/xy**4 &
&                                   + f0*(xx**10+yy**10)/xy**5
   else if (order==6) then
     f0 = two          *(three +two*xy**2)/(three      *sqrtpi)*sumexp
     f1 = four   *xy   *(5._dp +two*xy**2)/(15._dp     *sqrtpi)*sumexp
     f2 = 8._dp  *xy**2*(7._dp +two*xy**2)/(105._dp    *sqrtpi)*sumexp
     f3 = 16._dp *xy**2*(9._dp +two*xy**2)/(945._dp    *sqrtpi)*sumexp
     f4 = 32._dp *xy**2*(11._dp+two*xy**2)/(10395._dp  *sqrtpi)*sumexp
     f5 = 64._dp *xy**2*(13._dp+two*xy**2)/(135135._dp *sqrtpi)*sumexp
     f6 = 128._dp*xy**2*(15._dp+two*xy**2)/(1027025._dp*sqrtpi)*sumexp
     hh=((xx**13+yy**13)*erfcp-(xx**13-yy**13)*erfcm)/(two*xy**7)
     screened_coul_kernel = hh + f6 + f5*(xx**2+yy**2)  /xy    + f4*(xx**4+yy**4)  /xy**2 &
&                                   + f3*(xx**6+yy**6)  /xy**3 + f1*(xx**8+yy**8)  /xy**4 &
&                                   + f1*(xx**10+yy**10)/xy**5 + f0*(xx**12+yy**12)/xy**6
   end if

 end select

end function screened_coul_kernel
!!***

!----------------------------------------------------------------------

!!****f* m_pawrad/calc_slatradl
!! NAME
!!  calc_slatradl
!!
!! FUNCTION
!!  Calculate the radial part of Slater integrals. See below.
!!
!! INPUTS
!!  ll= l quantum number in the expansion of the Coulomb term.
!!  mesh_size=Number of points on the radial mesh.
!!  ff1(radmesh), ff2(radmesh)= The two functions to be integrated.
!!  Pawrad <type(pawrad_type)>=Structure containing radial grid information.
!!
!! OUTPUT
!!  integral =
!!   $ \dfrac{4\pi}{2L+1} \int ff1(r1) \dfrac{r_<^L}{r_>^{L+1}} ff2(r2) dr1 dr2 $
!!  where $r_< = min(r1,r2)$ and $r_> = Max(r1,r2)$.
!!
!! PARENTS
!!      m_paw_slater
!!
!! CHILDREN
!!      poisson,simp_gen
!!
!! SOURCE

subroutine calc_slatradl(ll,mesh_size,ff1,ff2,Pawrad,integral)

!scalars
 integer,intent(in) :: mesh_size,ll
 real(dp),intent(out) :: integral
!arrays
 real(dp),intent(in) :: ff1(mesh_size),ff2(mesh_size)
 type(pawrad_type),intent(in) :: Pawrad

!Local variables ---------------------------------------
!scalars
 integer :: int_meshsz
 character(len=100) :: msg
!arrays
 real(dp),allocatable :: hh(:),gg(:)

! *************************************************************************

 if (mesh_size > Pawrad%mesh_size) then
   msg='mesh_size > pawrad%mesh_size!'
   LIBPAW_BUG(msg)
 end if

 LIBPAW_ALLOCATE(hh,(mesh_size))
 LIBPAW_ALLOCATE(gg,(mesh_size))
 !hh = zero
 !gg = zero

 int_meshsz=Pawrad%int_meshsz
 !$int_meshsz=Pawrad%mesh_size

 ! the line below requires hh as work array.
 hh = ff2

 ! TODO find where int_meshsz is calculated and if it can affects the results.
 if (int_meshsz<mesh_size) hh(int_meshsz+1:mesh_size)=zero

 call poisson(hh,ll,Pawrad,gg)

 gg(2:mesh_size) = gg(2:mesh_size)/Pawrad%rad(2:mesh_size)

 hh(1)          = zero
 hh(2:mesh_size)= ff1(2:mesh_size) * gg(2:mesh_size)
 LIBPAW_DEALLOCATE(gg)

 call simp_gen(integral,hh,Pawrad)
 integral = four_pi * integral

 LIBPAW_DEALLOCATE(hh)

end subroutine calc_slatradl
!!***

!----------------------------------------------------------------------

END MODULE m_pawrad
!!***
