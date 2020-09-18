!!****m* ABINIT/m_gsphere
!! NAME
!!  m_gsphere
!!
!! FUNCTION
!!   The Gsphere data type defines the set of G-vectors
!!   centered on Gamma used to describe (chi0|epsilon|W) in the GW code.
!!   Note that, unlike the kg_k arrays used for wavefunctions, here the
!!   G-vectors are ordered in shells (increasing length). Moreover
!!   the sphere can be enlarged to take into account umklapps for which
!!   one need the knowledge of several quantities at G-G0.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group (MG, GMR, VO, LR, RWG, MT, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_gsphere

 use defs_basis
 use m_abicore
 use m_errors
 use m_distribfft
 use m_sort

 use defs_abitypes,   only : MPI_type
 use m_fstrings,      only : sjoin, itoa
 use m_numeric_tools, only : bisect
 use m_geometry,      only : normv
 use m_crystal,       only : crystal_t
 use m_fftcore,       only : kpgsph, kgindex, sphereboundary
 use m_mpinfo,        only : destroy_mpi_enreg, initmpi_seq

 implicit none

 private

! Low-level procedures.
 public :: merge_and_sort_kg   ! Merges a set of k-centered G-spheres of cutoff ecut. Return a Gamma-centered G-spheres.
 public :: table_gbig2kg       ! Associate the kg_k set of G-vectors with Gamma-centered G-sphere.
 public :: get_irredg          ! Given a set of G vectors, find the set of G"s generating the others by symmetry.
 public :: merge_kgirr         ! Merge a list of irreducible G vectors (see routine for more info)
 public :: setshells           ! Set consistently the number of shells, the number of plane-waves,  and the energy cut-off
 public :: kg_map              ! Compute the mapping between two lists of g-vectors.
 public :: make_istwfk_table
 public :: getkpgnorm          ! Compute the norms of the k+G vectors
 public :: symg
!!***

!----------------------------------------------------------------------

!!****t* m_gsphere/gsphere_t
!! NAME
!! gsphere_t
!!
!! FUNCTION
!! The gsphere_t data type contains information related to the set of G vectors
!! used during a screening or a GW calculation, as well as symmetry tables relating
!! these vectors. Presently the following quantities are stored
!!
!! 1) The reduced coordinates of the G vectors (arrays gvec)
!! 2) Tables giving the correspondence between a G-vector and its rotated image
!!    through a symmetry operation in reciprocal space.
!! 3) List of the irreducible G pairs
!! 4) Tables giving, for each pair in the full reciprocal space, the corresponding
!!    irreducible pair as well as the symmetry operation in reciprocal space
!!
!! Note that, unlike the GS part, the basis set does not depend on the k-point.
!!
!! NOTES
!! To indicate the indices in the arrays grottb, grottbm1 we use the following notation :
!!
!!  g defines the index of the reciprocal lattice vector in the array gvec
!!  s  indicates the index of the symmetry operation in reciprocal space
!!  i  can  be one or two. 1 is used to indicate the identity operator
!!
!! SOURCE

 type,public :: gsphere_t

  integer :: ng
  ! Total number of G vectors in the sphere taking into account umklapps
  ! it defines the size of the array gvec and it accounts for possible umklapps for which
  ! one has to shift the sphere.

  ! TODO: The sphere should be enlarged in gpairs_init using mg0 and (ecut|ng) as input.
  ! For the time being we keep the old implementation.
  !%integer :: ng_eff
  ! Effective number of G vectors, i.e. the number of G in the smaller sphere without umklapps.
  ! ng_eff<=ng and should be used to loop over the elements of (chi0|epsilon|W).
  !
  ! TODO: Add info on FFT including zero-padding algorithm.
  ! table must be recalculated for each G0 and rho_tw_g should accept gpshere in input.

  integer :: nsh                   ! Number of shells
  integer :: nsym                  ! The number of symmetry operations
  integer :: timrev                ! 2 if time-reversal is used, 1 otherwise
  integer :: istwfk=1              ! Storage mode. At present time-reversal is not used.

  !integer :: mg0(3)=0
  ! For each reduced direction gives the max G0 component to account for umklapp processes.

  real(dp) :: ecut
  ! Cutoff energy of the sphere.

  real(dp) :: gmet(3,3)
  ! Reciprocal space metric ($\textrm{bohr}^{-2}$).

  real(dp) :: gprimd(3,3)
  ! Dimensional reciprocal space primitive translations (bohr^{-1})

  integer,allocatable :: g2sh(:)
  ! g2sh(ng)
  ! For each G, it gives the index of the shell to which it belongs.

  integer,allocatable :: gvec(:,:)
  ! gvec(3,ng)
  ! Reduced coordinates of G vectors.

  integer,allocatable :: g2mg(:)
  ! g2mg(ng)
  ! Correspondence G --> -G

  integer,allocatable :: rottb(:,:,:)
  ! rottb(ng,timrev,nsym)
  ! rottb(G,I,S) is the index of (SI) G in the array gvec
  ! where I is either the identity or the inversion.

  integer,allocatable :: rottbm1(:,:,:)
  ! rottb(ng,timrev,nsym)
  ! rottbm1(G,I,S) is the index of IS{^-1} G in the array gvec

  integer,allocatable :: shlim(:)
  ! shlim(nsh+1)
  ! Index of the first G vector in each shell, =ng+1 for nsh+1

  real(dp),allocatable :: shlen(:)
  ! shlen(nsh)
  ! Radius of each shell.

  !TODO switch to dpc
  complex(gwpc),allocatable :: phmGt(:,:)
  ! phmGt(ng,nsym)
  ! Phase factor e^{-i2\pi(G.\tau)} where $\tau$ is the fractional translation associated to isym.

  complex(gwpc),allocatable :: phmSGt(:,:)
  ! phmSGt(ng,nsym)
  ! Phase factor e^{-i2\pi(SG.\tau)} where S is one of the symmetry properties in reciprocal space.

 end type gsphere_t

 public :: gsph_init          ! Initialize the G-sphere.
 public :: gsph_fft_tabs      ! Returns useful tables for FFT (with or without padding).
 public :: gsph_in_fftbox     ! Initialize the largest Gsphere contained in the FFT box.
 public :: print_gsphere      ! Printout of basic dimensions.
 public :: gsph_free          ! Free memory allocated in the object.
 public :: gsph_g_idx         ! Returns the index of G from its reduced coordinates.
 public :: gsph_gmg_idx       ! Returns the index of G1-G2 from their indeces
 public :: gsph_gmg_fftidx    ! Returns the index of G1-G2 in the FFT mesh defined by ngfft.
 public :: gsph_extend        ! Construct a new gsphere_t with a larger cutoff energy
!!***

CONTAINS  !=================================================================================
!!***

!!****f* m_gsphere/setup_G_rotation
!! NAME
!! setup_G_rotation
!!
!! FUNCTION
!! Set up tables indicating rotation of G-vectors.
!!
!! INPUTS
!! nsym=Number of symmetry operations.
!! symrec(3,3,nsym)=Symmetry operations in reciprocal space.
!! timrev=2 if time reversal can be used, 1 otherwise.
!! npw=Number of planewaves in the sphere.
!! gvec(3,npw)=Coordinates of plane waves, supposed to be ordered in increasing modulus
!! g2sh(npw)=For each G, it gives the index of the shell to which it belongs.
!! nsh=Number of shells
!! shlim(nsh+1)=Index of the first G vector in each shell, =npw+1 for nsh+1
!!
!! OUTPUT
!!  grottb  (npw,2,nsym)= grottb(G,I,S) is the index of (SI) G in the array gvec.
!!  grottbm1(npw,2,nsym)= index of IS^{-1} G.
!!
!! NOTES:
!!  I is either the identity or the inversion (time reversal in reciprocal space).
!!  S is one of the symmetry operation in reciprocal space belonging to the Space group.
!!
!! PARENTS
!!      m_gsphere
!!
!! CHILDREN
!!
!! SOURCE

subroutine setup_G_rotation(nsym,symrec,timrev,npw,gvec,g2sh,nsh,shlim,grottb,grottbm1)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,nsh,nsym,timrev
!arrays
 integer,intent(in) :: g2sh(npw),gvec(3,npw),shlim(nsh+1),symrec(3,3,nsym)
 integer,intent(inout) :: grottb  (npw,timrev,nsym)
 integer,intent(inout) :: grottbm1(npw,timrev,nsym)

!Local variables ------------------------------
!scalars
 integer :: ee,ig1,ig2,ish1,isym,itim,ss
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: gbase(3),grot(3)

!************************************************************************
 !
 ! === Set up G-rotation table ===
 do ig1=1,npw
   ish1=g2sh(ig1) ; ss=shlim(ish1) ; ee=shlim(ish1+1)-1
   gbase(:)=gvec(:,ig1)

   do itim=1,timrev
     do isym=1,nsym
       grot=(3-2*itim)*MATMUL(symrec(:,:,isym),gbase)
       found=.FALSE.
       ! * Loop on the shell of ig1 to speed up the search.
       do ig2=ss,ee
         if (ALL(ABS(grot(:)-gvec(:,ig2))==0)) then
           found=.TRUE.
           grottb  (ig1,itim,isym)=ig2
           grottbm1(ig2,itim,isym)=ig1
         end if
       end do
       if (.not.found) then
         write(msg,'(3a,i5,a,i5,1x,2(3i5,a),a,i3,a,i3)')&
&         'G-shell not closed',ch10,&
&         '  Initial G vector ',ig1,'/',npw,gbase(:),' Rotated G vector ',grot(:),ch10,&
&         '  Through sym ',isym,' and itim ',itim
         MSG_ERROR(msg)
       end if
     end do
   end do

 end do !ig1

end subroutine setup_G_rotation
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/gsph_init
!! NAME
!! gsph_init
!!
!! FUNCTION
!!  Main creation method for the Gvectors data type
!!
!! INPUTS
!!  Cryst<crystal_t> = Info on unit cell and its symmetries
!!     %nsym=number of symmetry operations
!!     %symrec(3,3,nsym)=symmetry operations in reciprocal space
!!     %tnons(3,nsym)=fractional translations
!!     %gmet(3,3)=reciprocal space metric (bohr**-2).
!!     %gprimd(3,3)=dimensional reciprocal space primitive translations
!!  ng=number of G vectors, needed only if gvec is passed.
!!  [gvec(3,ng)]=coordinates of G vectors
!!  [ecut]=Cutoff energy for G-sphere. gvec and ecut are mutually exclusive.
!!
!! OUTPUT
!!  Gsph<gsphere_t>=Data type containing information related to the set of G vectors
!!   completetly initialized in output.
!!
!! NOTES
!!  gvec are supposed to be ordered with increasing norm.
!!
!! PARENTS
!!      m_bethe_salpeter,m_gsphere,m_gwls_hamiltonian,m_screening_driver
!!      m_sigma_driver,mrgscr
!!
!! CHILDREN
!!
!! SOURCE

subroutine gsph_init(Gsph,Cryst,ng,gvec,ecut)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng
 real(dp),optional,intent(in) :: ecut
 type(crystal_t),target,intent(in) :: Cryst
 type(gsphere_t),intent(out) :: Gsph
!arrays
 integer,optional,intent(in) :: gvec(3,ng)
!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt1=1
 integer :: ig,isearch,img,ish,isym,nsh,nsym,timrev,pinv,g1,g2,g3,ss,ee
 real(dp) :: eps,norm,norm_old,max_ecut,gsq
!arrays
 real(dp),parameter :: k_gamma(3)=(/zero,zero,zero/)
 integer :: sg(3),gsearch(3)
 integer,allocatable :: shlim(:)
 integer,pointer :: symrec(:,:,:),gvec_ptr(:,:)
 real(dp) :: kptns1(3,nkpt1)
 real(dp),allocatable :: shlen(:)
 real(dp),pointer :: tnons(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 !@gsphere_t
 ! === Copy info on symmetries ===
 nsym   =  Cryst%nsym
 timrev =  Cryst%timrev
 symrec => Cryst%symrec
 tnons  => Cryst%tnons
 !
 ! === Initialize the object ===
 Gsph%istwfk = 1           ! Time reversal is not used here.
 Gsph%nsym   = nsym
 Gsph%timrev = timrev

 Gsph%gmet   = Cryst%gmet
 Gsph%gprimd = Cryst%gprimd

 if (PRESENT(gvec)) then
   if (PRESENT(ecut)) then
     MSG_BUG("ecut cannot be present when gvec is used")
   end if
   Gsph%ng= ng
   ABI_MALLOC(Gsph%gvec,(3,ng))
   Gsph%gvec=gvec
   !
   ! Calculate cutoff energy.of the sphere.
   max_ecut=-one
   do ig=1,ng
     g1=gvec(1,ig)
     g2=gvec(2,ig)
     g3=gvec(3,ig)
     gsq=       Cryst%gmet(1,1)*g1**2+Cryst%gmet(2,2)*g2**2+Cryst%gmet(3,3)*g3**2+ &
&          two*(Cryst%gmet(1,2)*g1*g2+Cryst%gmet(1,3)*g1*g3+Cryst%gmet(2,3)*g2*g3)
     max_ecut=MAX(max_ecut,gsq)
   end do
   max_ecut=two*max_ecut*pi**2
   Gsph%ecut= max_ecut

 else
   ! To be consistent with the previous implementation.
   MSG_WARNING("Init from ecut has to be tested")
   !call setshells(ecut,npw,nsh,nsym,Cryst%gmet,Cryst%gprimd,Cryst%symrel,tag,Cryst%ucvol)
   Gsph%ecut = ecut
   pinv=+1; kptns1(:,1)=k_gamma
   call merge_and_sort_kg(nkpt1,kptns1,ecut,Cryst%nsym,pinv,Cryst%symrel,Cryst%gprimd,gvec_ptr,0)
   Gsph%ng = SIZE(gvec_ptr,DIM=2)
   ABI_MALLOC(Gsph%gvec, (3,Gsph%ng))
   Gsph%gvec = gvec_ptr
   ABI_FREE(gvec_ptr)
 end if
 !
 ! === Calculate phase exp{-i2\pi G.\tau} ===
 ABI_MALLOC(Gsph%phmGt,(Gsph%ng,nsym))
 do ig=1,Gsph%ng
   do isym=1,nsym
    Gsph%phmGt(ig,isym)=EXP(-j_dpc*two_pi*DOT_PRODUCT(Gsph%gvec(:,ig),tnons(:,isym)))
   end do
 end do
 !
 ! === Calculate phase phsgt= exp{-i2\pi SG\cdot t} ===
 ! TODO Here we can store only one of this arrays but I have to rewrite screeening!
 ABI_MALLOC(Gsph%phmSGt,(Gsph%ng,nsym))
 do ig=1,Gsph%ng
   do isym=1,nsym
     sg=MATMUL(symrec(:,:,isym),Gsph%gvec(:,ig))
     Gsph%phmSGt(ig,isym)=EXP(-j_dpc*two_pi*DOT_PRODUCT(sg,tnons(:,isym)))
   end do
 end do
 !
 ! === Calculate number of shells and corresponding starting index ===
 ! * Shells are useful to speed up search algorithms see e.g setup_G_rotation.
 ! * The last shell ends at ng+1, thus gvec is supposed to be closed.

 ABI_CHECK(ALL(Gsph%gvec(1:3,1)==0),'First G must be 0')

 ABI_MALLOC(Gsph%g2sh,(Gsph%ng))
 Gsph%g2sh(1)=1 ! This table is useful if we dont loop over shell

 ! For each shell, gives the index of the initial G-vector.
 ABI_MALLOC(shlim,(Gsph%ng+1))
 shlim(1)=1

 ! For each shell, gives the radius of the shell.
 ABI_MALLOC(shlen,(Gsph%ng))
 shlen(1)=zero

 nsh=1; norm_old=zero
 do ig=2,Gsph%ng
   norm=two_pi*SQRT(DOT_PRODUCT(Gsph%gvec(:,ig),MATMUL(Cryst%gmet,Gsph%gvec(:,ig))))
   eps=norm*tol8
   if (ABS(norm-norm_old)>eps) then
     norm_old = norm; nsh = nsh + 1
     shlim(nsh)=ig
     shlen(nsh)=norm
   end if
   Gsph%g2sh(ig)=nsh
 end do
 shlim(nsh+1)=Gsph%ng+1

 ! === Save info on the shells ===
 Gsph%nsh=nsh
 ABI_MALLOC(Gsph%shlim,(nsh+1))
 Gsph%shlim=shlim(1:nsh+1)
 ABI_MALLOC(Gsph%shlen,(nsh  ))
 Gsph%shlen=shlen(1:nsh)
 ABI_FREE(shlim)
 ABI_FREE(shlen)
 !
 ! === Calculate tables for rotated G"s ===
 ABI_MALLOC(Gsph%rottb  ,(Gsph%ng,timrev,nsym))
 ABI_MALLOC(Gsph%rottbm1,(Gsph%ng,timrev,nsym))

 call setup_G_rotation(nsym,symrec,timrev,Gsph%ng,Gsph%gvec,&
&  Gsph%g2sh,Gsph%nsh,Gsph%shlim,Gsph%rottb,Gsph%rottbm1)

 ! Store Mapping G --> -G
 ! (we use a specialized table instead of rootb since rottb assumes time-reversal symmetry.
 ABI_MALLOC(gsph%g2mg, (gsph%ng))

 do ig=1,gsph%ng
   ish=gsph%g2sh(ig)
   ss=gsph%shlim(ish); ee=gsph%shlim(ish+1)-1
   gsearch = -gsph%gvec(:,ig)
   img = 0
   ! Loop on shells to speed up the search.
   do isearch=ss,ee
     if (all(gsph%gvec(:,isearch) == gsearch)) then
       img = isearch; exit
     end if
   end do
   if (img==0) MSG_ERROR("Cannot find -G in G-sphere!")
   gsph%g2mg(ig) = img
 end do

 !call print_gsphere(Gsph,unit=std_out,prtvol=1)

 DBG_EXIT("COLL")

end subroutine gsph_init
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/gsph_fft_tabs
!! NAME
!! gsph_fft_tabs
!!
!! FUNCTION
!!
!! INPUTS
!!  Gsph<gsphere_t>=Info on the G-sphere
!!  g0(3)
!!  mgfft=MAXVAL(ngfft(1:3))
!!  ngfftf(18)=Info on the FFT mesh.
!!
!! OUTPUT
!!  use_padfft=1 if padded FFT can be used, 0 otherwise.
!!  gmg0_gbound(2*mgfft+8,2)=Tables for improved zero-padded FFTS. Calculated only if use_padfft==1
!!  gmg0_ifft(Gsph%ng)=Index of G-G0 in the FFT mesh defined by ngfft.
!!
!! NOTES
!!  The routine will stop if any G-G0 happens to be outside the FFT box.
!!
!! PARENTS
!!      m_chi0,m_cohsex,m_exc_build,m_prep_calc_ucrpa,m_sigc,m_sigx
!!
!! CHILDREN
!!
!! SOURCE

subroutine gsph_fft_tabs(Gsph,g0,mgfft,ngfft,use_padfft,gmg0_gbound,gmg0_ifft)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft
 integer,intent(out) :: use_padfft
 type(gsphere_t),intent(in) :: Gsph
!arrays
 integer,intent(in) :: g0(3),ngfft(18)
 integer,intent(out) :: gmg0_gbound(2*mgfft+8,2),gmg0_ifft(Gsph%ng)

!Local variables-------------------------------
!scalars
 integer :: ig,ng,ierr
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer,allocatable :: gmg0(:,:)
 logical,allocatable :: kg_mask(:)

! *************************************************************************

 if (mgfft/=MAXVAL(ngfft(1:3))) then
   MSG_ERROR("mgfft/-MAXVAL(ngfft(1:3)")
 end if

 ng = Gsph%ng

 ierr=0; use_padfft=0
 ABI_MALLOC(gmg0,(3,ng))
 do ig=1,ng
   gmg0(:,ig) = Gsph%gvec(:,ig)-g0
   ! Consider possible wrap around errors.
   if ( ANY(gmg0(:,ig)>ngfft(1:3)/2) .or. ANY(gmg0(:,ig)<-(ngfft(1:3)-1)/2) ) then
     !gmg0_ifft(ig,ig01+mg0(1)+1,ig02+mg0(2)+1,ig03+mg0(3)+1) = 0
     write(std_out,*)" outside FFT box ",gmg0(:,ig)
     ierr=ierr+1
   end if
   if (ALL(gmg0(:,ig) == 0)) use_padfft=1
 end do

 if (ierr/=0) then
   write(msg,'(a,i0,a)')'Found ',ierr,' G-G0 vectors falling outside the FFT box. This is not allowed '
   MSG_ERROR(msg)
 end if
 !
 ! Evaluate the tables needed for the padded FFT performed in rhotwg. Note that we have
 ! to pass G-G0 to sphereboundary instead of G as we need FFT results on the shifted G-sphere,
 ! If Gamma is not inside G-G0 one has to disable FFT padding as sphereboundary will give wrong tables.
 if (use_padfft == 1) call sphereboundary(gmg0_gbound,1,gmg0,mgfft,ng)

 call initmpi_seq(MPI_enreg_seq) ! No FFT parallelism.
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfft(2),ngfft(3),'all')

 ABI_MALLOC(kg_mask, (ng))
 call kgindex(gmg0_ifft, gmg0, kg_mask, MPI_enreg_seq, ngfft, ng)

 ABI_CHECK(ALL(kg_mask),"FFT para not yet implemented")
 ABI_FREE(kg_mask)

 ABI_FREE(gmg0)
 call destroy_mpi_enreg(MPI_enreg_seq)

end subroutine gsph_fft_tabs
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/gsph_in_fftbox
!! NAME
!! gsph_in_fftbox
!!
!! FUNCTION
!!  Initialize the largest Gsphere contained in the FFT box.
!!
!! INPUTS
!!  Cryst<crystal_t> = Info on unit cell and its symmetries.
!!  ngfft(18)=Info on the FFT box.
!!
!! OUTPUT
!!  Gsph<gsphere_t>=Data type containing information related to the set of G vectors
!!   completetly initialized in output.
!!
!! PARENTS
!!      m_chi0
!!
!! CHILDREN
!!
!! SOURCE

subroutine gsph_in_fftbox(Gsph,Cryst,ngfft)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: Cryst
 type(gsphere_t),intent(out) :: Gsph
!arrays
 integer,intent(in) :: ngfft(18)

!Local variables-------------------------------
!scalars
 integer :: dir1,dir2,dir3,npw,ig,ist
 real(dp) :: ecut,trial_ene
!arrays
 integer :: n1_max(3),n2_max(3),n3_max(3),vec(3)
 integer,allocatable :: gvec(:,:)

!************************************************************************

 ! Find ecut for the largest G-sphere contained in the FFT box.
 n1_max(1) = -(ngfft(1)-1)/2
 n2_max(1) = -(ngfft(2)-1)/2
 n3_max(1) = -(ngfft(3)-1)/2

 n1_max(2) = 0
 n2_max(2) = 0
 n3_max(2) = 0

 n1_max(3) = ngfft(1)/2
 n2_max(3) = ngfft(2)/2
 n3_max(3) = ngfft(3)/2

 ecut = HUGE(one)
 do dir3=1,3
   vec(3) = n1_max(dir3)
   do dir2=1,3
     vec(2) = n2_max(dir2)
     do dir1=1,3
       vec(1) = n1_max(dir1)
       if (ANY(vec/=0)) then
         trial_ene = half * normv(vec,Cryst%gmet,"G")**2
         ecut = MIN(ecut,trial_ene)
         !write(std_out,*)vec(:),trial_ene
       end if
     end do
   end do
 end do
 !
 ! Init sphere from ecut.
 call gsph_init(Gsph,Cryst,0,ecut=ecut)
 !
 ! Make sure that Gsph does not contain G vectors outside the FFT box.
 ! kpgsph might return G whose energy is larger than the input ecut.
 npw = Gsph%ng
 star_loop: do ist=1,Gsph%nsh-1
   do ig=Gsph%shlim(ist),Gsph%shlim(ist+1)
     if ( ANY(Gsph%gvec(:,ig)>ngfft(1:3)/2) .or. ANY(Gsph%gvec(:,ig)<-(ngfft(1:3)-1)/2) ) then
       npw = Gsph%shlim(ist)-1  ! Gsph exceeds the FFT box. Only the shells up to npw will be used.
       EXIT star_loop
     end if
   end do
 end do star_loop

 if (npw<Gsph%ng) then
   MSG_COMMENT("Have to reinit Gpshere")
   ABI_MALLOC(gvec,(3,npw))
   gvec =Gsph%gvec(:,1:npw)
   call gsph_free(Gsph)
   call gsph_init(Gsph,Cryst,npw,gvec=gvec)
   ABI_FREE(gvec)
 end if

end subroutine gsph_in_fftbox
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/print_gsphere
!! NAME
!! print_gsphere
!!
!! FUNCTION
!!  Print the content of a gvectors data type
!!
!! INPUTS
!!  Gsph<gsphere_t>=Info on the G-sphere
!!  unit=the unit number for output
!!  prtvol = verbosity level
!!  mode_paral =either "COLL" or "PERS"
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_bethe_salpeter,m_chi0,m_gwls_hamiltonian
!!
!! CHILDREN
!!
!! SOURCE

subroutine print_gsphere(Gsph,unit,prtvol,mode_paral)

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 type(gsphere_t),intent(in) :: Gsph

!Local variables-------------------------------
!scalars
 integer :: ish,nsc,my_unt,my_prtvol
 real(dp) :: fact,kin
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt    =std_out; if (PRESENT(unit      )) my_unt    =unit
 my_prtvol=0       ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode   ='COLL' ; if (PRESENT(mode_paral)) my_mode   =mode_paral

 write(msg,'(3a,2(a,i8,a))')ch10,&
& ' ==== Info on the G-sphere ==== ',ch10,&
& '  Number of G vectors ... ',Gsph%ng,ch10,&
& '  Number of shells ...... ',Gsph%nsh,ch10
 call wrtout(my_unt,msg,my_mode)

 SELECT CASE (Gsph%timrev)
 CASE (1)
   call wrtout(my_unt,' Time reversal symmetry cannot be used',my_mode)
 CASE (2)
   call wrtout(my_unt,' Time reversal symmetry is used',my_mode)
 CASE DEFAULT
   MSG_BUG("Wrong timrev")
 END SELECT

 if (my_prtvol/=0) then
   fact=half*two_pi**2
   write(msg,'(a)')
   call wrtout(my_unt,' Shell   Tot no. of Gs   Cutoff [Ha]',my_mode)
   do ish=1,Gsph%nsh
     nsc=Gsph%shlim(ish+1)-1
     kin=half*Gsph%shlen(ish)**2
     write(msg,'(2x,i4,10x,i6,5x,f8.3)')ish,nsc,kin
     call wrtout(my_unt,msg,'COLL')
   end do
   call wrtout(my_unt,ch10,my_mode)
 end if

end subroutine print_gsphere
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/gsph_free
!! NAME
!! gsph_free
!!
!! FUNCTION
!!  Deallocate the memory in a gsphere_t data type.
!!
!! INPUTS
!!   Gsph = datatype to be freed
!!
!! PARENTS
!!      m_bethe_salpeter,m_chi0,m_gsphere,m_gwls_hamiltonian,m_screening_driver
!!      m_sigma_driver,mrgscr
!!
!! CHILDREN
!!
!! SOURCE

subroutine gsph_free(Gsph)

!Arguments ------------------------------------
!scalars
 type(gsphere_t),intent(inout) :: Gsph

! *************************************************************************

 DBG_ENTER("COLL")

 !@gsphere_t

! integer arrays.
 ABI_SFREE(Gsph%g2sh)
 ABI_SFREE(Gsph%gvec)
 ABI_SFREE(Gsph%g2mg)
 ABI_SFREE(Gsph%rottb)
 ABI_SFREE(Gsph%rottbm1)
 ABI_SFREE(Gsph%shlim)

 ABI_SFREE(Gsph%shlen)

! complex arrays
 ABI_SFREE(Gsph%phmGt)
 ABI_SFREE(Gsph%phmSGt)

 DBG_EXIT("COLL")

end subroutine gsph_free
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/gsph_g_idx
!! NAME
!! gsph_g_idx
!!
!! FUNCTION
!! Return the index of G in the sphere. zero if not in the sphere
!!
!! INPUTS
!!  Gsph<gsphere_t>=Info on the G-sphere
!!  gg(3)=Reduced coordinates of the G-vector.
!!
!! NOTES
!!  The function assumes that the G-vectors are ordered with increasing length.
!!
!! PARENTS
!!
!! SOURCE

pure function gsph_g_idx(Gsph,gg) result(g_idx)

!Arguments ------------------------------------
!scalars
 type(gsphere_t),intent(in) :: Gsph
 integer :: g_idx
!arrays
 integer,intent(in) :: gg(3)

!Local variables-------------------------------
!scalars
 integer :: ishbsc,igs,ige
 real(dp) :: glen
 logical :: found

! *************************************************************************

 ! * Use shells and bisect to find the star and stop index thus avoiding the storage of a table (ig1,ig2)
 glen = two_pi*SQRT(DOT_PRODUCT(gg,MATMUL(Gsph%gmet,gg)))

 ishbsc = bisect(Gsph%shlen,glen)
 if ( ANY(ishbsc==(/0,Gsph%nsh/)) ) then ! glen out of range.
   g_idx=0; RETURN
 end if

 igs = Gsph%shlim(ishbsc)
 ige = Gsph%shlim(MIN(ishbsc+2,Gsph%nsh+1))-1

 g_idx=igs-1; found=.FALSE.
 do while (.not.found .and. g_idx<ige)
   g_idx=g_idx+1
   found=(ALL(Gsph%gvec(:,g_idx)==gg(:)))
 end do
 if (.not.found) g_idx=0

end function gsph_g_idx
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/gsph_gmg_idx
!! NAME
!! gsph_gmg_idx
!!
!! FUNCTION
!! Return the index of G1-G2 in the sphere. zero if not in the sphere
!!
!! INPUTS
!!  Gsph<gsphere_t>=Info on the G-sphere
!!  ig1,ig2 index of g1 and g2 in the G-sphere.
!!
!! NOTES
!!  The function assumes that the G-vectors are ordered with increasing length.
!!
!! PARENTS
!!
!! SOURCE

pure function gsph_gmg_idx(Gsph,ig1,ig2) result(ig1mg2)

!Arguments ------------------------------------
!scalars
 type(gsphere_t),intent(in) :: Gsph
 integer,intent(in) :: ig1,ig2
 integer :: ig1mg2

!Local variables-------------------------------
!scalars
 integer :: ishbsc,igs,ige
 real(dp) :: difflen
 logical :: found
!arrays
 integer :: g1mg2(3)

! *************************************************************************

 g1mg2 = Gsph%gvec(:,ig1)-Gsph%gvec(:,ig2)

 ! * Use shells and bisect to find the star and stop index thus avoiding the storage of a table (ig1,ig2)
 difflen = two_pi*SQRT(DOT_PRODUCT(g1mg2,MATMUL(Gsph%gmet,g1mg2)))

 ! FIXME It seems bisect is not portable, on my laptop test v5/t72 the number of skipped G-vectors is > 0
 ishbsc = bisect(Gsph%shlen,difflen)
 if ( ANY(ishbsc==(/0,Gsph%nsh/)) ) then ! difflen out of range.
   ig1mg2=0; RETURN
 end if

 igs = Gsph%shlim(ishbsc)
 ige = Gsph%shlim(MIN(ishbsc+2,Gsph%nsh+1))-1

 ig1mg2=igs-1; found=.FALSE.
 do while (.not.found .and. ig1mg2<ige)
   ig1mg2=ig1mg2+1
   found=(ALL(Gsph%gvec(:,ig1mg2)==g1mg2(:)))
 end do
 if (.not.found) ig1mg2=0

end function gsph_gmg_idx
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/gsph_gmg_fftidx
!! NAME
!! gsph_gmg_fftidx
!!
!! FUNCTION
!! Return the index of G1-G2 in the FFT mesh defined by ngfft. zero if not found.
!!
!! INPUTS
!!  Gsph<gsphere_t>=Info on the G-sphere
!!  ig1,ig2 index of g1 and g2 in the G-sphere.
!!  ngfft(18)=Info on the FFT mesh.
!!
!! PARENTS
!!
!! SOURCE

pure function gsph_gmg_fftidx(Gsph,ig1,ig2,ngfft) result(fft_idx)

!Arguments ------------------------------------
!scalars
 type(gsphere_t),intent(in) :: Gsph
 integer,intent(in) :: ig1,ig2
 integer :: fft_idx
!arrays
 integer,intent(in) :: ngfft(18)

!Local variables-------------------------------
!scalars
 integer :: id1,id2,id3
!arrays
 integer :: g1mg2(3)

! *************************************************************************

 g1mg2(:)=Gsph%gvec(:,ig1)-Gsph%gvec(:,ig2)

 ! Make sure G1-G2 is still in the FFT mesh.
 ! MODULO wraps G1-G2 in the FFT box but the Fourier components are not periodic!
 if (ANY(g1mg2(:)>ngfft(1:3)/2) .or. ANY(g1mg2(:)<-(ngfft(1:3)-1)/2)) then
   fft_idx=0; RETURN
 end if

 id1=MODULO(g1mg2(1),ngfft(1))
 id2=MODULO(g1mg2(2),ngfft(2))
 id3=MODULO(g1mg2(3),ngfft(3))
 fft_idx= 1 + id1 + id2*ngfft(1) + id3*ngfft(1)*ngfft(2)

end function gsph_gmg_fftidx
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/merge_and_sort_kg
!! NAME
!!  merge_and_sort_kg
!!
!! FUNCTION
!!  This routine merges a set of k-centered G-spheres of cutoff energy ecut and
!!  returns a Gamma-centered G-spheres. The elements in the final G-spheres are packed with increasing module.
!!
!! INPUTS
!!  nkpt=Number of k-points
!!  kptns(3,nkpt)=The k-points in reduced coordinates defining the k-centered G-spheres.
!!  ecut=Cutoff energy for the k-centered G-spheres.
!!  nsym2=Number of symmetry operations.
!!  pinv=-1 if time-reversal can be used, 1 otherwise
!!  symrel2(3,3,nsym2)=symmetry operations in real space.
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  prtvol=Flag defining the verbosity level.
!!
!! SIDE EFFECTS
!!  gbig(:,:)
!!    in input : pointer to NULL
!!    in output: gbig(3,1:npw) contains the set of G-vectors ordered by shell obtained by
!!               merging the k-centered sphere.
!!  shlim_p(:)
!!    in input : pointer to NULL
!!    in output: shlim_p(nbase)=Cumulative number of G-vectors for each shell.
!!               where nbase is the number of irreducible G"s found.
!!
!! PARENTS
!!      m_gsphere,m_io_kss,m_sigma_driver
!!
!! CHILDREN
!!
!! SOURCE

subroutine merge_and_sort_kg(nkpt,kptns,ecut,nsym2,pinv,symrel2,gprimd,gbig,prtvol,shlim_p)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nsym2,pinv,prtvol
 real(dp),intent(in) :: ecut
!arrays
 integer,intent(in) :: symrel2(3,3,nsym2)
 real(dp),intent(in) :: kptns(3,nkpt),gprimd(3,3)
 integer,pointer :: gbig(:,:)
 integer,optional,pointer :: shlim_p(:)

!Local variables-------------------------------
!scalars
 integer,parameter :: mkmem_=1
 integer :: ikg,ig,ikpt,nbase,sizepw,in,maxpw,is,iinv,ish,ilim,mpw
 integer :: exchn2n3d,istwf_k,onpw_k,ierr,npw_k,ii,isym,sizeold
 logical :: found
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: gcur(3),geq(3)
 integer :: symrec2t(3,3,nsym2)
 integer :: dum_kg(3,0)
 integer,allocatable :: gbase(:,:),gbasek(:,:,:)
 integer,allocatable :: gcurr(:,:),gshell(:,:),insort(:),gtmp(:,:)
 integer,allocatable :: nbasek(:),nshell(:),shlim(:)
 integer,allocatable :: npwarr(:)
 real(dp) :: kpoint(3),gmet(3,3)
 real(dp),allocatable :: cnorm(:),cnormk(:,:),ctmp(:)

! *********************************************************************

 ! * Fake MPI_type for the sequential part.
 ! This routine should not be parallelized as communicating gbig and other
 ! tables takes more time than recalculating them in sequential.
 call initmpi_seq(MPI_enreg_seq)

!Compute reciprocal space metrics
 do ii=1,3
   gmet(ii,:)=gprimd(1,ii)*gprimd(1,:)+&
&   gprimd(2,ii)*gprimd(2,:)+&
&   gprimd(3,ii)*gprimd(3,:)
 end do

 ! * Here we use TRANSPOSE(symrel2) instead of the more intuitive symrel2^{-1t} for historical reasons
 ! It does not affect the results since in the code below we only check the module of G
 do isym=1,nsym2
   symrec2t(:,:,isym)=TRANSPOSE(symrel2(:,:,isym))
 end do
 !
 ! ==============================================
 ! ==== Find irreducible G-vectors at each k ====
 ! ==============================================

 ABI_MALLOC(npwarr,(nkpt))
 exchn2n3d=0; ikg=0
 do ikpt=1,nkpt
   kpoint=kptns(:,ikpt); istwf_k=1
   call kpgsph(ecut,exchn2n3d,gmet,ikg,0,istwf_k,dum_kg,kpoint,0,MPI_enreg_seq,0,npwarr(ikpt))
 end do
 mpw = MAXVAL(npwarr)

 ABI_MALLOC(nbasek,(nkpt))
 ABI_MALLOC(gbasek,(3,mpw,nkpt))
 ABI_MALLOC(cnormk,(mpw,nkpt))
 nbasek=0     ! # of irreducible G at each k.
 cnormk=zero  ! Norm of each irreducible G.
 gbasek=0     ! The set of irreducible G"s at each k.

 do ikpt=1,nkpt

   kpoint = kptns(:,ikpt)
   npw_k  = npwarr(ikpt)

   exchn2n3d=0; ikg=0; istwf_k=1
   ABI_MALLOC(gcurr,(3,npw_k))
   call kpgsph(ecut,exchn2n3d,gmet,ikg,0,istwf_k,gcurr,kpoint,mkmem_,MPI_enreg_seq,npw_k,onpw_k)

   if (ANY(gcurr(:,1)/=0)) then
     MSG_BUG("gcurr(:,1)/=0")
   end if
   !
   ! * Search for the G"s generating the others by symmetry.
   !  NB: Here we use symrec2t=TRANSPOSE(symrel2) for historical reasons, see note above
   call get_irredg(npw_k,nsym2,pinv,gprimd,symrec2t,gcurr,nbasek(ikpt),gbasek(:,:,ikpt),cnormk(:,ikpt))

   ABI_FREE(gcurr)
 end do
 !
 ! === Reduce info over k-points ===
 ! * Here symrec2t=TRANSPOSE(symrel2) for historical reasons, see note above
 sizepw=2*mpw
 ABI_MALLOC(gbase,(3,sizepw))
 ABI_MALLOC(cnorm,(sizepw))
 nbase=0                    ! # of irred G found.

 call merge_kgirr(nsym2,pinv,nkpt,mpw,sizepw,symrec2t,nbasek,cnormk,gbasek,nbase,gbase,cnorm,ierr)
 if (ierr/=0) then
   MSG_ERROR('merge_kgirr returned a non-zero status error')
 end if

 ABI_FREE(nbasek)
 ABI_FREE(cnormk)
 ABI_FREE(gbasek)
 !
 !=== Reorder base G-vectors in order of increasing module ===
 !
 !Generate all shells of G-vectors: star of a g==set of all symetrics of this g
 ABI_MALLOC(gshell,(3,2*nsym2))
 ABI_MALLOC(shlim,(nbase))

 ABI_MALLOC(gbig,(3,sizepw))
!
!TODO
#if 0
!* Here symrec2t=TRANSPOSE(symrel2) for historical reasons, see note above
 call getfullg(nbase,nsym2,pinv,sizepw,gbase,symrec2t,cnorm,maxpw,gbig,shlim,ierr)
 if (ierr/0) RETURN

#else
 ABI_MALLOC(insort,(nbase))
 ABI_MALLOC(nshell,(nbase))
 do in=1,nbase
   insort(in)=in
 end do
 call sort_dp(nbase,cnorm,insort,tol14)
!
!Loop over all different modules of g''s (=shells):
 maxpw=0
 do in=1,nbase
   nshell(in)=0
   gcur(:)=gbase(:,insort(in))

   do is=1,nsym2 !  Loop over all symetries:
     do iinv=pinv,1,2
       geq(:)=iinv*(symrel2(1,:,is)*gcur(1)+symrel2(2,:,is)*gcur(2)+symrel2(3,:,is)*gcur(3))

       found=.FALSE.; ish=1
       do while ((.not.found) .and. (ish<=nshell(in))) ! Search for symetric of g and eventually add it:
         found=ALL(geq(:)==gshell(:,ish))
         ish=ish+1
       end do
       if (.not.found) then
         nshell(in)=nshell(in)+1
         gshell(:,nshell(in))=geq(:)
       end if
     end do
   end do

   if ((maxpw+nshell(in)) > sizepw) then
     ! We need to increase the size of the gbase, gbig and cnorm arrays while still keeping their content.
     ! This is done using two temporary arrays gtmp and ctmp
     MSG_WARNING("Had to reallocate gbase, gbig, cnorm")
     ABI_MALLOC(ctmp,(sizepw))
     ABI_MALLOC(gtmp,(3,sizepw))
     sizeold=sizepw
     sizepw=maxpw+nshell(in)

     ctmp(:)=cnorm(:)
     gtmp(:,:)=gbase(:,:)

     ABI_FREE(cnorm)
     ABI_MALLOC(cnorm,(sizepw))
     cnorm(1:sizeold)=ctmp(1:sizeold)
     cnorm(sizeold+1:sizepw)=zero
     ABI_FREE(ctmp)

!    MG why this? gbase should not be changed!
     ABI_FREE(gbase)
     ABI_MALLOC(gbase,(3,sizepw))
     gbase(:,:sizeold)=gtmp(:,:sizeold)
     gbase(:,sizeold+1:sizepw)=0
     gtmp(:,:)=gbig(:,:)

     ABI_FREE(gbig)
     ABI_MALLOC(gbig,(3,sizepw))
     gbig(:,:sizeold)=gtmp(:,:sizeold)
     gbig(:,sizeold+1:sizepw)=0
     ABI_FREE(gtmp)
   end if
   !
   ! Store this shell of g''s in a big array of g (gbig):
   do ig=1,nshell(in)
     gbig(:,ig+maxpw)=gshell(:,ig)
   end do
   maxpw=maxpw+nshell(in)

 end do ! End loop over shells
 !
 ! * Compute shell limits
 ilim=0
 do in=1,nbase
   ilim=ilim+nshell(in)
   shlim(in)=ilim
 end do

 if (PRESENT(shlim_p)) then ! Return shlim_p
  ABI_MALLOC(shlim_p,(nbase))
  shlim_p = shlim
 end if

 ! Re-allocate gbig with correct sizes so that caller can inquire the size
 ABI_MALLOC(gtmp,(3,ilim))
 gtmp = gbig(:,1:ilim)
 ABI_FREE(gbig)
 ABI_MALLOC(gbig,(3,ilim))
 gbig=gtmp
 ABI_FREE(gtmp)

 if (prtvol>10) then ! Print out shell limits
   write(msg,'(3a)')&
&    ' Shells found:',ch10,&
&    ' number of shell    number of G vectors      cut-off energy [Ha} '
   call wrtout(std_out,msg,'COLL')
   do in=1,nbase
     write(msg,'(12x,i4,17x,i6,12x,f8.3)')in,shlim(in),2*pi**2*cnorm(in)
     call wrtout(std_out,msg,'COLL')
   end do
   call wrtout(std_out,ch10,'COLL')
 end if

 ABI_FREE(gshell)
 ABI_FREE(insort)
 ABI_FREE(nshell)
#endif

 call destroy_mpi_enreg(MPI_enreg_seq)
 ABI_FREE(gbase)
 ABI_FREE(shlim)
 ABI_FREE(cnorm)
 ABI_FREE(npwarr)

end subroutine merge_and_sort_kg
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/getfullg
!! NAME
!! getfullg
!!
!! FUNCTION
!!  Reconstruct a G-sphere starting from a set of irreducible lattice vectors
!!
!! INPUTS
!!  pinv=-1 if time-reversal can be used, 1 otherwise
!!  nsym=number of symmetry operations
!!  sizepw=Max expected number of G vectors in the shere
!!  symrec(3,3,nsym)=symmetry operation in reciprocal space
!!  nbase=number of irreducible G vectors
!!  gbase(3,nbase)=irreducible G-vectors
!!  cnorm(nbase)=norm of the irreducible G vectors (supposed not yet sorted)
!!
!! OUTPUT
!!  maxpw=Number of G vectors found
!!  gbig(3,sizepw)=G vectors in the sphere packed in the first maxpw columns
!!  shlim(nbase)=number of G vectors within each shell
!!  ierr= Exit status, if /=0 the number of G vectors found exceeds sizepw
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  cnorm is a bit redundant since it can be calculated from gbase. However this procedure
!!  is called by outkss in which cnorm is already calculated and we dont want to do it twice
!!
!! PARENTS
!!      m_gsphere
!!
!! CHILDREN
!!
!! SOURCE

subroutine getfullg(nbase,nsym,pinv,sizepw,gbase,symrec,cnorm,maxpw,gbig,shlim,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbase,nsym,pinv,sizepw
 integer,intent(out) :: ierr,maxpw
!arrays
 integer,intent(in) :: gbase(3,nbase),symrec(3,3,nsym)
 integer,intent(out) :: gbig(3,sizepw),shlim(nbase)
 real(dp),intent(inout) :: cnorm(nbase) !sort_dp can change cnorm

!Local variables-------------------------------
!scalars
 integer :: ibase,ig,ilim,ish,isym,itim
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: gcur(3),geq(3)
 integer,allocatable :: gshell(:,:),insort(:),nshell(:)

! *************************************************************************

 if (pinv/=1.and.pinv/=-1) then
   write(msg,'(a,i6)')&
&   ' The argument pinv should be -1 or 1, however, pinv =',pinv
   MSG_BUG(msg)
 end if
 !
 ! === Reorder base g-vectors in order of increasing module ===
 ABI_MALLOC(insort,(nbase))
 do ibase=1,nbase
   insort(ibase)=ibase
 end do
 call sort_dp(nbase,cnorm,insort,tol14)
 !
 ! === Generate all stars of G-vectors ===
 ! Star of G is the set of all symetrical images of the vector
 ! gshell contains the symmetrical G at fixed gbase. No need to add an additional dimension
 ! or initialize to zero the array inside the loop over nbase as we loop over (ish<=nshell(ibase))
 ABI_MALLOC(nshell,(nbase))
 ABI_MALLOC(gshell,(3,2*nsym))
 !
 ! === Start with zero number of G vectors found ===
 maxpw=0 ; ierr=0
 do ibase=1,nbase
   !
   ! === Loop over all different modules of G ===
   ! * Start with zero G vectors found in this star
   nshell(ibase)=0
   gcur(:)=gbase(:,insort(ibase))
   !
   !  === Loop over symmetries ===
   do isym=1,nsym
     do itim=pinv,1,2
       geq(:)=itim*MATMUL(symrec(:,:,isym),gcur)
       !
       ! * Search for symetric of g and eventually add it:
       found=.FALSE. ; ish=1
       do while ((.not.found).and. (ish<=nshell(ibase)))
         found=ALL(geq(:)==gshell(:,ish))
         ish=ish+1
       end do
       if (.not.found) then
         nshell(ibase)=nshell(ibase)+1
         gshell(:,nshell(ibase))=geq(:)
       end if
     end do
   end do
   !
   ! * Was sizepw large enough?
   if ((maxpw+nshell(ibase))>sizepw) then
     write(msg,'(a,i6,2a)')&
&     ' Number of G in sphere exceeds maximum allowed value =',sizepw,ch10,&
&     ' check the value of sizepw in calling routine '
     MSG_WARNING(msg)
     ierr=1; RETURN
   end if
   !
   ! === Store this shell of Gs in a big array (gbig) ===
   do ig=1,nshell(ibase)
     gbig(:,ig+maxpw)=gshell(:,ig)
   end do
   maxpw=maxpw+nshell(ibase)
 end do ! ibase
 !
 ! === Compute number of G"s within each shell ===
 ilim=0
 do ibase=1,nbase
   ilim=ilim+nshell(ibase)
   shlim(ibase)=ilim
 end do
 !
 ! === Print out shell limits ===
 write(msg,'(3a)')&
& ' Shells found:',ch10,&
& ' number of shell    number of G vectors      cut-off energy [Ha] '
 call wrtout(std_out,msg,'COLL')

 do ibase=1,nbase
   write(msg,'(12x,i4,17x,i6,12x,f8.3)')ibase,shlim(ibase),two*pi**2*cnorm(ibase)
   call wrtout(std_out,msg,'COLL')
 end do
 write(msg,'(a)')ch10
 call wrtout(std_out,msg,'COLL')
 ABI_FREE(gshell)
 ABI_FREE(insort)
 ABI_FREE(nshell)

end subroutine getfullg
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/get_irredg
!! NAME
!! get_irredg
!!
!! FUNCTION
!!  Given a set of reciprocal lattice vectors, find the set of G"s generating the others by symmetry.
!!
!! INPUTS
!!  nsym=number of symmetry operations
!!  pinv=-1 if time-reversal can be used, 1 otherwise
!!  npw_k=number of G vectors (for this k-point, as the set of G is k-centered)
!!  gcurr(3,npw_k)=the list of G vectors
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  symrec(3,3,nsym)=symmetry operations in terms of reciprocal space primitive translations.
!!
!! OUTPUT
!!  nbasek=number of irreducible G vectors found
!!  cnormk(npw_k)=first nbasek elements are the norm of each irreducible G-vector
!!  gbasek(3,npw_k)=first nbasek elements are the irreducible G vectors
!!
!! NOTES
!!  The search can be optimized by looping over shells. See m_skw for a faster algo
!!
!! PARENTS
!!      m_gsphere,m_skw
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_irredg(npw_k,nsym,pinv,gprimd,symrec,gcurr,nbasek,gbasek,cnormk)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,nsym,pinv
 integer,intent(out) :: nbasek
!arrays
 integer,intent(in) :: gcurr(3,npw_k),symrec(3,3,nsym)
 integer,intent(out) :: gbasek(3,npw_k)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(out) :: cnormk(npw_k)

!Local variables-------------------------------
!scalars
 integer :: ig,irr,isym,jj
 real(dp) :: eps,norm
 logical :: found
!arrays
 integer :: gbas(3),gcur(3),geq(3)
 real(dp) :: gcar(3)

! *************************************************************************

 DBG_ENTER("COLL")

 if (pinv/=1.and.pinv/=-1) then
   MSG_BUG(sjoin('pinv should be -1 or 1, however, pinv =', itoa(pinv)))
 end if

 ! Zero irred G vectors found, zeroing output arrays.
 nbasek = 0; cnormk(:) = zero; gbasek(:,:) = 0

 do ig=1,npw_k
   gcur(:) = gcurr(:,ig); norm = zero
   do jj=1,3
     gcar(jj)=gcur(1)*gprimd(jj,1)+gcur(2)*gprimd(jj,2)+gcur(3)*gprimd(jj,3)
     norm=norm+gcar(jj)**2
   end do
   eps = tol8 * norm; found = .False.; irr = 1
   do while (.not.found .and. irr <= nbasek)  ! This loop can be optimized by looping inside the shell.
     if (abs(norm - cnormk(irr)) <= eps) then
       gbas(:) = gbasek(:,irr); isym = 1
       do while (.not.found .and. isym <= nsym)
         geq(:) = matmul(symrec(:,:,isym),gcur)
         found = all(geq(:) == gbas(:))
         if (pinv == -1) found = (found .or. all(geq == -gbas)) ! For time-reversal
         isym = isym + 1
       end do
     end if
     irr = irr + 1
   end do
   if (.not. found) then
     nbasek = nbasek + 1; cnormk(nbasek) = norm; gbasek(:,nbasek) = gcur(:)
   end if
 end do

 DBG_EXIT("COLL")

end subroutine get_irredg
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/merge_kgirr
!! NAME
!! merge_kgirr
!!
!! FUNCTION
!!  Given a list of irreducible reciprocal vectors associated to different k-centered spheres,
!!  this subroutine finds the minimal set of G vectors needed to reconstruct the union of the spheres
!!  through symmetry operations.
!!
!! INPUTS
!! nsym=number of symmetry operations
!! pinv=-1 if time-reversal can be used, 0 otherwise
!! nkpt=number of k-points for k-centered spheres
!! mpw=Max number of G vectors for each k-point
!! sizepw=Max expected number of G vectors found
!! symrec(3,3,nsym)=symmetries in reciprocal space given in reduced coordinates
!! nbasek(nkpt)=number of irred G for each k-point
!! cnormk(mpw,nkpt)=the norm of each k-centered G (only 1:nbase(ik)) is used
!! gbasek(3,mpw,nkpt)
!!
!! OUTPUT
!! nbase=number of irreducible G needed to reconstruct the initial set of spheres
!! gbase(3,sizepw)=irreducible G found in reciprocal coordinates
!! cnorm(sizepw)=Norm of each irred G vector
!! ierr= Exit status, if /=0 the number of G vectors found exceeds sizepw
!!
!! PARENTS
!!      m_gsphere
!!
!! CHILDREN
!!
!! SOURCE

subroutine merge_kgirr(nsym,pinv,nkpt,mpw,sizepw,symrec,nbasek,cnormk,gbasek,nbase,gbase,cnorm,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpw,nkpt,nsym,pinv,sizepw
 integer,intent(out) :: ierr,nbase
!arrays
 integer,intent(in) :: gbasek(3,mpw,nkpt),nbasek(nkpt),symrec(3,3,nsym)
 integer,intent(inout) :: gbase(3,sizepw) !vz_i
 real(dp),intent(in) :: cnormk(mpw,nkpt)
 real(dp),intent(inout) :: cnorm(sizepw) !vz_i

!Local variables-------------------------------
!scalars
 integer :: ikpt,inb,irgk,isym
 real(dp) :: eps,norm
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: gbas(3),gcur(3),geq(3)

! *************************************************************************

 DBG_ENTER("COLL")

 if (pinv/=1.and.pinv/=-1) then
   write(msg,'(a,i6)')' The argument pinv should be -1 or 1, however, pinv =',pinv
   MSG_BUG(msg)
 end if
 !
 ! === Start with zero number of G found ===
 nbase=0 ; ierr=0
 do ikpt=1,nkpt
   do irgk=1,nbasek(ikpt)
     gcur(:)=gbasek(:,irgk,ikpt)
     norm=cnormk(irgk,ikpt) ; eps=tol8*norm
     found=.FALSE. ; inb=1
     do while ((.not.found).and.(inb<=nbase))
       if (ABS(norm-cnorm(inb))<=eps) then
         gbas(:)=gbase(:,inb)
         isym=1
         do while ((.not.found).and.(isym<=nsym))
           geq(:)=MATMUL(symrec(:,:,isym),gcur)
           found=ALL(geq(:)==gbas(:))
           if (pinv==-1) found= (found.or.ALL(geq(:)==-gbas(:)) ) ! For time-reversal
           isym=isym+1
         end do
       end if
       inb=inb+1
     end do
     if (.not.found) then
       ! === Add to the list ===
       nbase=nbase+1
       if (nbase>sizepw) then
         write(msg,'(2(a,i5),a)')&
&         ' nbase (',nbase,') became greater than sizepw = ',sizepw,' returning ierr=1 '
         MSG_WARNING(msg)
         ierr=1; RETURN
       end if
       cnorm(nbase)=cnormk(irgk,ikpt)
       gbase(:,nbase)=gcur(:)
     end if
   end do
 end do

 DBG_EXIT("COLL")

end subroutine merge_kgirr
!!***

!----------------------------------------------------------------------

!!****f* m_gpshere/setshells
!! NAME
!! setshells
!!
!! FUNCTION
!! Set consistently the number of shells, the number of plane-waves, and the energy cut-off
!!
!! INPUTS
!!  nsym=number of symmetry operations
!!  gmet(3,3)=metric tensor in reciprocal space
!!  gprimd(3,3)=dimensional primitive vectors in reciprocal space
!!  symrel(3,3,nsym)=symmetry operations in real space
!!  tag=suffix to account for the different possibilities for these variables (npw, ecut or nsh ..)
!!  ucvol=unit cell volume
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  ecut,npw,nsh=one of them is an input, the two others are output
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  npw=number of plane waves
!!  nsh=number of shells
!!
!! PARENTS
!!      m_invars2,m_screening_driver,m_sigma_driver
!!
!! CHILDREN
!!
!! SOURCE

subroutine setshells(ecut,npw,nsh,nsym,gmet,gprimd,symrel,tag,ucvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(inout) :: npw,nsh
 real(dp),intent(in) :: ucvol
 real(dp),intent(inout) :: ecut
 character(len=*),intent(in) :: tag
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: exchn2n3d,ifound,ig,ii,ish,isym,npw_found,npwave
 integer :: npwwrk,nsh_found,pad=50
 real(dp) :: ecut_found,ecut_trial,eps,scale=1.3_dp
 logical :: found
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: geq(3)
 integer,allocatable :: gvec(:,:),gvec_sh(:,:),insort(:),npw_sh(:)
 real(dp) :: gctr(3)
 real(dp),allocatable :: gnorm(:),gnorm_sh(:)

!******************************************************************

 ! Check coherence of input variables ecut, npw, and nsh.
 ! 1-> one at least should be non-null
 if (npw==0.and.nsh==0.and.ecut<=tol6) then
   write(msg,'(8a)')&
    'One of the three variables ecut',TRIM(tag),', npw',TRIM(tag),', or nsh',TRIM(tag),ch10,&
    'must be non-null. Returning.'
   MSG_COMMENT(msg)
   RETURN
 end if
 ! 2-> one and only one should be non-null
 if (npw/=0.and.nsh/=0) then
   write(msg,'(6a)')&
    'Only one of the two variables npw',TRIM(tag),' and nsh',TRIM(tag),ch10,&
    'can be non-null. Modify the value of one of these in input file.'
   MSG_ERROR(msg)
 end if
 if (ecut>tol6.and.npw/=0) then
   write(msg,'(6a)')&
    'Only one of the two variables ecut',TRIM(tag),' and npw',TRIM(tag),ch10,&
    'can be non-null. Modify the value of one of these in input file.'
   MSG_ERROR(msg)
 end if
 if (ecut>tol6.and.nsh/=0) then
   write(msg,'(6a)')&
    'Only one of the two variables ecut',TRIM(tag),' and nsh',TRIM(tag),ch10,&
    'can be non-null Action : modify the value of one of these in input file.'
   MSG_ERROR(msg)
 end if

 ! Calculate an upper bound for npw.
 ! gctr is center of the g-vector sphere
 gctr(:)= [zero,zero,zero]
 if (ecut>tol6) then
   ! The average number of plane-waves in the cutoff sphere is given by:
   ! npwave = (2*ecut)**(3/2)*ucvol/(6*pi**2)
   ! The upper bound is calculated as npwwrk=int(scale * npwave) + pad
   npwave=NINT(ucvol*(two*ecut)**1.5_dp/(six*pi**2))
   npwwrk=NINT(DBLE(npwave)*scale)+pad
   ecut_trial=ecut
 else if (npw/=0) then
   ! npw is given in the input
   npwwrk=NINT(DBLE(npw)*scale)+pad
   ecut_trial=(six*pi**2*npw/ucvol)**two_thirds/two
 else
   ! If nsh is given in the input
   npwwrk=nsh*18+2*pad
   ecut_trial=(six*pi**2*nsh*18/ucvol)**two_thirds/two
 end if

 call initmpi_seq(MPI_enreg_seq)

 ABI_MALLOC(gvec,(3,npwwrk))
 ifound=0
 do while (ifound==0)
   !write(msg,'(a,f8.2)')' setshells : ecut_trial = ',ecut_trial
   !call wrtout(std_out,msg,'COLL')
   exchn2n3d=0 ! For the time being, no exchange of n2 and n3

   call kpgsph(ecut_trial,exchn2n3d,gmet,0,1,1,gvec,gctr,1,MPI_enreg_seq,npwwrk,npw_found)

   ABI_MALLOC(gnorm,(npw_found))
   ABI_MALLOC(insort,(npw_found))

   do ig=1,npw_found
     insort(ig)=ig
     gnorm(ig)=zero
     do ii=1,3
       gnorm(ig)=gnorm(ig)+(gvec(1,ig)*gprimd(ii,1)+&
                            gvec(2,ig)*gprimd(ii,2)+&
                            gvec(3,ig)*gprimd(ii,3))**2
     end do
   end do
   call sort_dp(npw_found,gnorm,insort,tol14)

   ABI_MALLOC(npw_sh,(npw_found))
   ABI_MALLOC(gnorm_sh,(npw_found))
   ABI_MALLOC(gvec_sh,(3,npw_found))
   npw_sh(:)=0
   gnorm_sh(:)=zero
   gvec_sh(:,:)=0
   ! Count the number of shells:
   ! (search for the G-vectors generating the others by symmetry)
   nsh_found=0

   do ig=1,npw_found
     eps=1.d-8*gnorm(ig)
     found=.FALSE.
     ish=1
     do while ((.not.found).and.(ish<=nsh_found))
       if (ABS(gnorm(ig)-gnorm_sh(ish))<=eps) then
         isym=1
         do while ((.not.found).and.(isym<=nsym))
           geq(:)=(symrel(1,:,isym)*gvec(1,insort(ig))+&
                  symrel(2,:,isym)*gvec(2,insort(ig))+&
                  symrel(3,:,isym)*gvec(3,insort(ig)))

           found=((geq(1)==gvec_sh(1,ish)).and.&
                  (geq(2)==gvec_sh(2,ish)).and.&
                  (geq(3)==gvec_sh(3,ish)))
           isym=isym+1
         end do
       end if
       ish=ish+1
     end do
     if (.not.found) then
       nsh_found=nsh_found+1
       gnorm_sh(nsh_found)=gnorm(ig)
       gvec_sh(:,nsh_found)=gvec(:,insort(ig))
       npw_sh(nsh_found)=1
     else
       ish=ish-1
       npw_sh(ish)=npw_sh(ish)+1
     end if
   end do

   ecut_found=two*pi**2*gnorm(npw_found)

   if(ecut>tol6) then
     ! ecut is given in the input
     if (ecut_found<ecut-0.1) then
       write(msg,'(3a,e14.6,9a,e14.6,3a)')&
        'The value ecut',TRIM(tag),'=',ecut,' given in the input file leads to',ch10,&
        'the same values for nsh',TRIM(tag),' and npw',TRIM(tag),' as ecut',TRIM(tag),'=',ecut_found,ch10!,&
!        'This value will be adopted for the calculation.',ch10
       MSG_COMMENT(msg)
     end if
     ifound=1
   else if (npw/=0) then
     ! If npw is given in the input
     if (npw_found==npw) then
       ecut_found=two*pi**2*gnorm(npw_found)
       ifound=1
     else if (npw_found>npw) then
       npw_found=0
       nsh_found=0
       do while (npw_found<npw)
         nsh_found=nsh_found+1
         npw_found=npw_found+npw_sh(nsh_found)
       end do
       ! check that the shell is closed
       if(npw_found>npw) then
         ! shell not closed
         npw_found=npw_found-npw_sh(nsh_found)
         nsh_found=nsh_found-1
         do while (ABS(gnorm_sh(nsh_found)-gnorm_sh(nsh_found+1))<0.000001)
           npw_found=npw_found-npw_sh(nsh_found)
           nsh_found=nsh_found-1
         end do
         write(msg,'(3a,i6,5a,i6,3a)')&
          'The value npw',TRIM(tag),'=',npw,' given in the input file does not close the shell',ch10,&
          'The lower closed-shell is obtained for a value npw',TRIM(tag),'=',npw_found,ch10,&
          'This value will be adopted for the calculation.',ch10
         MSG_WARNING(msg)
       end if
       ecut_found=two*pi**2*gnorm(npw_found)
       ifound=1
     end if
   else if (nsh/=0) then
     ! If nsh is given in the input
     if (nsh_found==nsh) then
       ecut_found=two*pi**2*gnorm(npw_found)
       ifound=1
     else if (nsh_found>nsh) then
       npw_found=0
       nsh_found=0
       do ish=1,nsh
         npw_found=npw_found+npw_sh(ish)
         nsh_found=nsh_found+1
       end do
       if (ABS(gnorm_sh(nsh_found)-gnorm_sh(nsh_found+1))<0.000001) then
         do while (ABS(gnorm_sh(nsh_found)-gnorm_sh(nsh_found+1))<0.000001)
           nsh_found=nsh_found+1
           npw_found=npw_found+npw_sh(nsh_found)
         end do
         write(msg,'(3a,i6,5a,i6,3a)')&
          'The value nsh',TRIM(tag),'=',nsh,' given in the input file corresponds to the same',ch10,&
          'cut-off energy as for closed-shell upto nsh',TRIM(tag),'=',nsh_found,ch10,&
          'This value will be adopted for the calculation.',ch10
         MSG_WARNING(msg)
       end if
       ecut_found=two*pi**2*gnorm(npw_found)
       ifound=1
     end if
   end if

   if (ifound==0) then
     ecut_trial=1.1*ecut_trial
     ABI_FREE(gnorm)
     ABI_FREE(gnorm_sh)
     ABI_FREE(gvec_sh)
     ABI_FREE(insort)
     ABI_FREE(npw_sh)
   else
     if (ecut<tol6) then
       ecut=ecut_found
     else
       MSG_COMMENT("MRM: ecut and/or ecutwfn was given in the input file and will not be resized.")
     end if
     npw=npw_found
     nsh=nsh_found
   end if

 end do ! while(ifound==0)

 call destroy_mpi_enreg(MPI_enreg_seq)

 ABI_FREE(gnorm)
 ABI_FREE(gnorm_sh)
 ABI_FREE(gvec)
 ABI_FREE(gvec_sh)
 ABI_FREE(insort)
 ABI_FREE(npw_sh)

end subroutine setshells
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/kg_map
!! NAME
!!  kg_map
!!
!! FUNCTION
!!  Compute the mapping between two lists of g-vectors.
!!
!! INPUTS
!!   npw1, kg1(3,npw1)=First list of G-vectors
!!   npw2, kg2(3,npw2)=Second list of G-vectors
!!
!! OUTPUT
!!   g2g1(npw2) = Mapping kg2 index --> kg1 index.
!!                Set to 0 if kg2(:,ig) not in kg1
!!   nmiss = Number of G-vectors in kg2 not found in kg1
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine kg_map(npw1,kg1,npw2,kg2,g2g1,nmiss)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw1,npw2
 integer,intent(out) :: nmiss
!arrays
 integer,intent(in) :: kg1(3,npw1),kg2(3,npw2)
 integer,intent(out) :: g2g1(npw2)

!Local variables ------------------------------
!scalars
 integer :: ii,ipw,i1,i2,i3,n1,n2,n3
!arrays
 integer :: gmax(3),g1_max(3),g2_max(3)
 integer,allocatable :: iwork(:,:,:)

!************************************************************************

 g1_max = maxval(abs(kg1))
 g2_max = maxval(abs(kg2))
 do ii=1,3
   gmax(ii) = max(g1_max(ii), g2_max(ii))
 end do
 gmax = 2*gmax + 1
 n1 = gmax(1); n2 = gmax(2); n3 = gmax(3)

 ABI_MALLOC(iwork, (n1, n2, n3))

 ! Insert kg1 into work with extra 0 s around outside:
 iwork = 0
 do ipw=1,npw1
   i1 = kg1(1,ipw); if (i1<0) i1=i1+n1; i1=i1+1
   i2 = kg1(2,ipw); if (i2<0) i2=i2+n2; i2=i2+1
   i3 = kg1(3,ipw); if (i3<0) i3=i3+n3; i3=i3+1
   iwork(i1,i2,i3) = ipw
 end do

 g2g1 = 0; nmiss = 0
 do ipw=1,npw2
   i1 = kg2(1,ipw); if (i1<0) i1=i1+n1; i1=i1+1
   i2 = kg2(2,ipw); if (i2<0) i2=i2+n2; i2=i2+1
   i3 = kg2(3,ipw); if (i3<0) i3=i3+n3; i3=i3+1
   g2g1(ipw) = iwork(i1,i2,i3)
   if (g2g1(ipw) == 0) nmiss = nmiss + 1
 end do

 ABI_FREE(iwork)

end subroutine kg_map
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/make_istwk_table
!! NAME
!! make_istwfk_table
!!
!! FUNCTION
!!
!! INPUTS
!!  ng1,ng2,ng3
!!
!! OUTPUT
!!
!! NOTES
!!   Useful relations:
!!     u_k(G) = u_{k+G0}(G-G0); u_{-k}(G) = u_k(G)^*
!!   and therefore:
!!     u_{G0/2}(G) = u_{G0/2}(-G-G0)^*.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_istwfk_table(istwf_k,ng1,ng2,ng3,ig1_inver,ig2_inver,ig3_inver)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng1,ng2,ng3,istwf_k
!arrays
 integer,intent(out) :: ig1_inver(ng1),ig2_inver(ng2),ig3_inver(ng3)

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3
 character(len=500) :: msg

!************************************************************************

! Initialize the inverse coordinates
 select case (istwf_k)

 case (1)
   ig1_inver(1)=1
   do i1=2,ng1
     ig1_inver(i1)=ng1+2-i1
   end do
   ig2_inver(1)=1
   do i2=2,ng2
     ig2_inver(i2)=ng2+2-i2
   end do
   ig3_inver(1)=1
   do i3=2,ng3
     ig3_inver(i3)=ng3+2-i3
   end do

 case (2:8)
   if (istwf_k==2 .or. istwf_k==4 .or. istwf_k==6 .or. istwf_k==8) then
     ig1_inver(1)=1
     do i1=2,ng1
       ig1_inver(i1)=ng1+2-i1
     end do
   else
     do i1=1,ng1
       ig1_inver(i1)=ng1+1-i1
     end do
   end if
   if (istwf_k>=2 .and. istwf_k<=5) then
     ig2_inver(1)=1
     do i2=2,ng2
       ig2_inver(i2)=ng2+2-i2
     end do
   else
     do i2=1,ng2
       ig2_inver(i2)=ng2+1-i2
     end do
   end if
   if (istwf_k==2 .or. istwf_k==3 .or. istwf_k==6 .or. istwf_k==7) then
     ig3_inver(1)=1
     do i3=2,ng3
       ig3_inver(i3)=ng3+2-i3
     end do
   else
     do i3=1,ng3
       ig3_inver(i3)=ng3+1-i3
     end do
   end if

 case default
   write(msg,'(a,i0)')" Wrong value for istwf_k: ",istwf_k
   MSG_ERROR(msg)
 end select

end subroutine make_istwfk_table
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/table_gbig2kg
!! NAME
!!  table_gbig2kg
!!
!! FUNCTION
!!  Associate the kg_k set of g-vectors with the big array of gbig
!!  The array gbig(3,maxpw) contains all g-vectors used for all k-points, in order of
!!  increasing shells. For a each k-point, the wave-functions are defined only on a particular set
!!  of g-vectors kg_k (included in gbig). This set is defined by array gamma2k:
!!  The array gamma2k(ig=1,maxpw) translates the index of the gbig (from 1 to maxpw) into the corresponding
!!  index in array kg_k. If gbig(ig) does not exist in kg_k, gamma2k(ig) contains npw_k+1.
!!
!! INPUTS
!!  npw_k=Number of planewaves in the k-centered basis set
!!  kg_k(3,npw_k)=The k-centered basis set
!!  maxpw=Number of G in gbig
!!  gbig(3,maxpw)=The union of the G-spheres at different k-points.
!!
!! OUTPUT
!!  ierr=Status error. It gives the number of G of kg_k not contained in gbig.
!!  gamma2k(maxpw)=Mapping gbig -> kg_k
!!
!! PARENTS
!!      m_io_kss
!!
!! CHILDREN
!!      gsph_free,gsph_init
!!
!! SOURCE

pure subroutine table_gbig2kg(npw_k,kg_k,maxpw,gbig,gamma2k,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,maxpw
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 integer,intent(in) :: gbig(3,maxpw)
 integer,intent(out) :: gamma2k(maxpw)

!Local variables-------------------------------
!scalars
 integer :: ig,igp
 logical :: found
!arrays
 integer :: gcur(3)

! *********************************************************************

 ierr=0
 gamma2k(:)=npw_k+1  ! Initialize array gamma2k

 do ig=1,npw_k       ! Loop over g-vectors, for this k point.
   gcur(:)=kg_k(:,ig)
   igp=0; found=.FALSE.
   do while ((.not.found) .and. igp<maxpw) ! Search selected vector in array gbig: TODO this part can be optimized
     igp=igp+1
     found=ALL(gcur(:)==gbig(:,igp))
   end do
   if (found) then ! Store it if found:
     gamma2k(igp)=ig
   else
     ierr=ierr+1
   end if
 end do

end subroutine table_gbig2kg
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/gsph_extend
!! NAME
!!  gsph_extend
!!
!! FUNCTION
!!  Construct a new gsphere_t with a larger cutoff energy
!!  while preserving the ordering of the first G-vectors stored in in_Gsph
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_bethe_salpeter
!!
!! CHILDREN
!!
!! SOURCE

subroutine gsph_extend(in_Gsph,Cryst,new_ecut,new_Gsph)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: new_ecut
 type(crystal_t),intent(in) :: Cryst
 type(gsphere_t),intent(in) :: in_Gsph
 type(gsphere_t),intent(out) :: new_Gsph
!arrays

!Local variables-------------------------------
!scalars
 integer :: new_ng,in_ng,ig,ierr,sh
!arrays
 integer,allocatable :: new_gvec(:,:)

! *********************************************************************

 call gsph_init(new_Gsph,Cryst,0,ecut=new_ecut)

 if (new_Gsph%ng > in_Gsph%ng) then

   new_ng = new_Gsph%ng
   in_ng  = in_Gsph%ng

   ierr = 0
   do ig=1,in_ng
     if (ANY(new_Gsph%gvec(:,ig) /= in_Gsph%gvec(:,ig)) ) then
       ierr = ierr + 1
       write(std_out,*)" new_gvec, in_gvec",ig,new_Gsph%gvec(:,ig),in_Gsph%gvec(:,ig)
     end if
   end do

   if (ierr==0) RETURN

   ierr = 0
   do sh=1,in_Gsph%nsh
     if ( new_Gsph%shlim(sh) /= in_Gsph%shlim(sh) .or. &
&         ABS(new_Gsph%shlen(sh)-in_Gsph%shlen(sh)) > tol12 ) then
       ierr = ierr + 1
       write(std_out,*)"new_shlim, in_shlim",sh,new_Gsph%shlim(sh),in_Gsph%shlim(sh)
       write(std_out,*)"new_shlen, in_shlen",sh,new_Gsph%shlen(sh),in_Gsph%shlen(sh)
     end if
   end do
   ABI_CHECK(ierr==0,"Wrong shells")

   ABI_MALLOC(new_gvec,(3,new_ng))
   new_gvec = new_Gsph%gvec
   new_gvec(:,1:in_ng) = in_Gsph%gvec

   call gsph_free(new_Gsph)
   call gsph_init(new_Gsph,Cryst,new_ng,gvec=new_gvec)
   ABI_FREE(new_gvec)

 else
   ierr = 0
   do ig=1,MIN(new_Gsph%ng,in_Gsph%ng)
     if (ANY(new_Gsph%gvec(:,ig) /= in_Gsph%gvec(:,ig)) ) then
       ierr = ierr + 1
       write(std_out,*)" new_gvec, in_gvec",ig,new_Gsph%gvec(:,ig),in_Gsph%gvec(:,ig)
     end if
   end do
   ABI_CHECK(ierr==0,"Fatal error")
 end if

end subroutine gsph_extend
!!***

!----------------------------------------------------------------------

!!****f* m_gsphere/getkpgnorm
!! NAME
!! getkpgnorm
!!
!! FUNCTION
!!  compute the norms of the k+G vectors
!!
!! INPUTS
!!  gprimd(3,3)=metric tensor
!!  kg_k(3,npw_k)= G vectors, in reduced coordinates
!!  kpt(3)=k vector, in reduced coordinates
!!  npw_k=size of the G-vector set
!!
!! OUTPUT
!!  kpgnorm(npw_k)=norms of the k+G vectors
!!
!! PARENTS
!!      m_cut3d,m_epjdos
!!
!! CHILDREN
!!
!! SOURCE

subroutine getkpgnorm(gprimd,kpt,kg_k,kpgnorm,npw_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(in) :: gprimd(3,3),kpt(3)
 real(dp),intent(out) :: kpgnorm(npw_k)

!Local variables-------------------------------
!scalars
 integer :: ipw
 real(dp) :: g11,g12,g13,g21,g22,g23,g31,g32,g33,k1,k2,k3,kpg1,kpg2,kpg3,rr,xx
 real(dp) :: yy,zz

! *************************************************************************

 k1=kpt(1) ; k2=kpt(2) ; k3=kpt(3)
 g11=gprimd(1,1)
 g12=gprimd(1,2)
 g13=gprimd(1,3)
 g21=gprimd(2,1)
 g22=gprimd(2,2)
 g23=gprimd(2,3)
 g31=gprimd(3,1)
 g32=gprimd(3,2)
 g33=gprimd(3,3)

!Loop over all k+G
 do ipw=1,npw_k

!  Load k+G
   kpg1=k1+dble(kg_k(1,ipw))
   kpg2=k2+dble(kg_k(2,ipw))
   kpg3=k3+dble(kg_k(3,ipw))

!  Calculate module of k+G
   xx=g11*kpg1+g12*kpg2+g13*kpg3
   yy=g21*kpg1+g22*kpg2+g23*kpg3
   zz=g31*kpg1+g32*kpg2+g33*kpg3
   rr=sqrt(xx**2+yy**2+zz**2)
   kpgnorm(ipw) = rr

 end do ! ipw

end subroutine getkpgnorm
!!***

!!****f* m_gsphere/symg
!! NAME
!! symg
!!
!! FUNCTION
!! Treat symmetries applied to the G vectors, in view of the application
!! to symmetrization of the dielectric matrix.
!! Generate a list of time-reversed G vectors, as well as a list
!! of spatially-symmetric G vectors.
!!
!! INPUTS
!! kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!! npwdiel=number of planewaves for the dielectric matrix
!! nsym=number of symmetry
!! symrel(3,3,nsym)=symmetry matrices in real space (integers)
!! tnons(3,nsym)=reduced nonsymmorphic translations
!! (symrel and tnons are in terms of real space primitive translations)
!!
!! OUTPUT
!! phdiel(2,npwdiel,nsym)=phase associated with point symmetries applied to G
!! sym_g(npwdiel,nsym)=index list of symmetric G vectors
!! (could save a bit of space by suppressing isym=1, since the
!! corresponding symmetry is the identity)
!! tmrev_g(npwdiel)=index list of inverted G vectors (time-reversed)
!!
!! PARENTS
!!      m_suscep_stat
!!
!! CHILDREN
!!
!! SOURCE

subroutine symg(kg_diel,npwdiel,nsym,phdiel,sym_g,symrel,tmrev_g,tnons)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwdiel,nsym
!arrays
 integer,intent(in) :: kg_diel(3,npwdiel),symrel(3,3,nsym)
 integer,intent(out) :: sym_g(npwdiel,nsym),tmrev_g(npwdiel)
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),intent(out) :: phdiel(2,npwdiel,nsym)

!Local variables-------------------------------
!scalars
 integer :: g1,g2,g3,ipw,isym,j1,j2,j3,m1m,m1p,m2m,m2p,m3m,m3p,symmg,trevg
 real(dp) :: arg,tau1,tau2,tau3
 !character(len=500) :: msg
!arrays
 integer,allocatable :: grid(:,:,:)

! *************************************************************************

!Determines maximal bounds of the zone spanned by the planewaves
 m1m=0 ; m2m=0 ; m3m=0 ; m1p=0 ; m2p=0 ; m3p=0
 do ipw=1,npwdiel
   g1=kg_diel(1,ipw)
   g2=kg_diel(2,ipw)
   g3=kg_diel(3,ipw)
   if(g1<m1m)m1m=g1 ; if(g1>m1p)m1p=g1
   if(g2<m2m)m2m=g2 ; if(g2>m2p)m2p=g2
   if(g3<m3m)m3m=g3 ; if(g3>m3p)m3p=g3
 end do

!Set up grid, that associate to each point the index of the
!corresponding planewave, if there is one
 ABI_MALLOC(grid, (m1m:m1p,m2m:m2p,m3m:m3p))
 grid(:,:,:)=0
 do ipw=1,npwdiel
   g1=kg_diel(1,ipw)
   g2=kg_diel(2,ipw)
   g3=kg_diel(3,ipw)
   grid(g1,g2,g3)=ipw
 end do

!Set up tmrev_g and sym_g arrays
 do ipw=1,npwdiel
   g1=kg_diel(1,ipw)
   g2=kg_diel(2,ipw)
   g3=kg_diel(3,ipw)

!  Treat first time-reversal symmetry
   trevg=grid(-g1,-g2,-g3)
   if(trevg==0)then
     MSG_BUG('Do not find the time-reversed symmetric of a G-vector.')
   end if
   tmrev_g(ipw)=trevg

!  Treat now spatial symmetries
   do isym=1,nsym

!    Get rotated G vector Gj for each symmetry element
!    -- here we use the TRANSPOSE of symrel; assuming symrel expresses
!    the rotation in real space, the transpose is then appropriate
!    for G space symmetrization (according to Doug : see routine irrzg.f)
     j1=symrel(1,1,isym)*g1+&
&     symrel(2,1,isym)*g2+symrel(3,1,isym)*g3
     j2=symrel(1,2,isym)*g1+&
&     symrel(2,2,isym)*g2+symrel(3,2,isym)*g3
     j3=symrel(1,3,isym)*g1+&
&     symrel(2,3,isym)*g2+symrel(3,3,isym)*g3
     symmg=grid(j1,j2,j3)
     if(symmg==0)then
       MSG_BUG('Do not find the spatially symmetric of a G-vector.')
     end if
     sym_g(ipw,isym)=symmg

!    Get associated phase
     tau1=tnons(1,isym)
     tau2=tnons(2,isym)
     tau3=tnons(3,isym)
     if (abs(tau1)>tol12.or.abs(tau2)>tol12.or.abs(tau3)>tol12) then
!      compute exp(-2*Pi*I*G dot tau) using original G
       arg=two_pi*(dble(g1)*tau1+dble(g2)*tau2+dble(g3)*tau3)
       phdiel(1,ipw,isym)=cos(arg)
       phdiel(2,ipw,isym)=-sin(arg)
     else
       phdiel(1,ipw,isym)=1._dp
       phdiel(2,ipw,isym)=0._dp
     end if

   end do
 end do

 ABI_FREE(grid)

end subroutine symg
!!***

END MODULE m_gsphere
!!***
