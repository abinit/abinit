!!****m* ABINIT/m_gwdefs
!! NAME
!! m_gwdefs
!!
!! FUNCTION
!! This module contains definitions for a number of named constants used in the GW part of abinit
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG, FB, GMR, VO, LR, RWG)
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

MODULE m_gwdefs

 use defs_basis
 use m_abicore
 use m_errors

 implicit none

 private

! Unit number for formatted files produced by GW calculations.
! These files are not supposed to be read by abinit therefore
! their names and unit numbers are not defined in the Dtfil% structure.
 integer,public,parameter :: unt_gw  = 21  ! GW corrections
 integer,public,parameter :: unt_sig = 22  ! Self-energy as a function of frequency
 integer,public,parameter :: unt_sgr = 23  ! Derivative wrt omega of the Self-energy
 integer,public,parameter :: unt_sgm = 20  ! Sigma on the Matsubara axis
 integer,public,parameter :: unt_gwdiag  = 40 ! GW diagonal

 real(dp),public,parameter :: GW_TOLQ =0.0001_dp
 ! Tolerance below which two BZ points are considered equal within a RL vector:
 ! for each red. direct. the abs value of the difference btw the two coord must be smaller that tolq.

 real(dp),public,parameter :: GW_TOLQ0=0.001_dp
 ! Tolerance below which a q-point is treated as zero (long wavelength limit)

 real(dp),public,parameter :: GW_TOL_DOCC=0.01_dp
 ! Tolerance on the difference between two occupation numbers.
 ! below this value, the contribution of the transition is neglected in the evaluation of chi0

 real(dp),public,parameter :: GW_TOL_W0=0.001_dp
 ! Tolerance on the real part of the frequency appearing in the denominator of the
 ! non-interacting Green function G0. Above this value, a small purely imaginary
 ! complex shift is added to the denominator during the evaluation of chi0.

 real(gwp),public,parameter :: one_gw =1._gwp
 real(gwp),public,parameter :: zero_gw=0._gwp

 complex(gwpc),public,parameter :: czero_gw=(0._gwp,0._gwp)
 complex(gwpc),public,parameter :: cone_gw =(1._gwp,0._gwp)
 complex(gwpc),public,parameter :: j_gw    =(0._gwp,1._gwp)

!arrays
 real(dp),public,parameter :: GW_Q0_DEFAULT(3) = (/0.00001_dp, 0.00002_dp, 0.00003_dp/)

! Weights and nodes for Gauss-Kronrod integration rules
! Gauss 7 Kronrod 15
 real(dp),public,parameter :: Kron15N(15) = (/-0.991455371120812639206854697526328516642040_dp, &
& -0.949107912342758524526189684047851262400770_dp,-0.864864423359769072789712788640926201210972_dp, &
& -0.741531185599394439863864773280788407074150_dp,-0.586087235467691130294144838258729598436780_dp, &
& -0.405845151377397166906606412076961463347382_dp,-0.207784955007898467600689403773244913479780_dp, &
&  0.000000000000000000000000000000000000000000_dp, 0.207784955007898467600689403773244913479780_dp, &
&  0.405845151377397166906606412076961463347382_dp, 0.586087235467691130294144838258729598436780_dp, &
&  0.741531185599394439863864773280788407074150_dp, 0.864864423359769072789712788640926201210972_dp, &
&  0.949107912342758524526189684047851262400770_dp, 0.991455371120812639206854697526328516642040_dp/)

 real(dp),public,parameter :: Kron15W(15) = (/0.022935322010529224963732008058969591993561_dp, &
& 0.063092092629978553290700663189204286665070_dp,0.104790010322250183839876322541518017443757_dp, &
& 0.140653259715525918745189590510237920399890_dp,0.169004726639267902826583426598550284106245_dp, &
& 0.190350578064785409913256402421013682826078_dp,0.204432940075298892414161999234649084716518_dp, &
& 0.209482141084727828012999174891714263697760_dp,0.204432940075298892414161999234649084716518_dp, &
& 0.190350578064785409913256402421013682826078_dp,0.169004726639267902826583426598550284106245_dp, &
& 0.140653259715525918745189590510237920399890_dp,0.104790010322250183839876322541518017443757_dp, &
& 0.063092092629978553290700663189204286665070_dp,0.022935322010529224963732008058969591993561_dp/)

 real(dp),public,parameter :: Gau7W(7) = (/0.12948496616886969327061143267908201832859_dp, &
& 0.279705391489276667901467771423779582486925_dp,0.38183005050511894495036977548897513387837_dp, &
& 0.417959183673469387755102040816326530612245_dp,0.38183005050511894495036977548897513387837_dp, &
& 0.279705391489276667901467771423779582486925_dp,0.12948496616886969327061143267908201832859_dp/)

! Gauss 11 Kronrod 23
 real(dp),public,parameter :: Kron23N(23) = (/-0.996369613889542634360164573335160773367030_dp, &
& -0.978228658146056992803938001122857390771420_dp,-0.941677108578067946455396730352134962188750_dp, &
& -0.887062599768095299075157769303927266631680_dp,-0.816057456656220942392261355192625879277860_dp, &
& -0.730152005574049324093416252031153458049643_dp,-0.630599520161965092168263312405390318822425_dp, &
& -0.519096129206811815925725669458609554480227_dp,-0.397944140952377573675073943298232277259112_dp, &
& -0.269543155952344972331531985400861524679620_dp,-0.136113000799361815798364355994952934344660_dp, &
&  0.000000000000000000000000000000000000000000_dp, 0.136113000799361815798364355994952934344660_dp, &
&  0.269543155952344972331531985400861524679620_dp, 0.397944140952377573675073943298232277259112_dp, &
&  0.519096129206811815925725669458609554480227_dp, 0.630599520161965092168263312405390318822425_dp, &
&  0.730152005574049324093416252031153458049643_dp, 0.816057456656220942392261355192625879277860_dp, &
&  0.887062599768095299075157769303927266631680_dp, 0.941677108578067946455396730352134962188750_dp, &
&  0.978228658146056992803938001122857390771420_dp, 0.996369613889542634360164573335160773367030_dp/)

 real(dp),public,parameter :: Kron23W(23) = (/0.00976544104596075802247917260964352216369_dp, &
& 0.027156554682104262051721401617851679412810_dp,0.04582937856442641598526161637154832107050_dp, &
& 0.063097424750374906584540530495371467781318_dp,0.07866457193222732928421712412387330212537_dp, &
& 0.092953098596900827769293667912429161939839_dp,0.10587207448138939648189189823198112055496_dp, &
& 0.116739502461047270810811060893282832324909_dp,0.125158799100319505060067189609770147044437_dp, &
& 0.131280684229805644255977537339957161670610_dp,0.135193572799884533184261853141533217156135_dp, &
& 0.136577794711118301018953895305516133510830_dp,0.135193572799884533184261853141533217156135_dp, &
& 0.131280684229805644255977537339957161670610_dp,0.125158799100319505060067189609770147044437_dp, &
& 0.116739502461047270810811060893282832324909_dp,0.10587207448138939648189189823198112055496_dp, &
& 0.092953098596900827769293667912429161939839_dp,0.07866457193222732928421712412387330212537_dp, &
& 0.063097424750374906584540530495371467781318_dp,0.04582937856442641598526161637154832107050_dp, &
& 0.027156554682104262051721401617851679412810_dp,0.00976544104596075802247917260964352216369_dp/)

 real(dp),public,parameter :: Gau11W(11) = (/0.0556685671161736664827537204425485787285_dp, &
& 0.125580369464904624634694299223940100197616_dp,0.186290210927734251426097641431655891691285_dp, &
& 0.233193764591990479918523704843175139431800_dp,0.262804544510246662180688869890509195372765_dp, &
& 0.272925086777900630714483528336342189156042_dp,0.262804544510246662180688869890509195372765_dp, &
& 0.233193764591990479918523704843175139431800_dp,0.186290210927734251426097641431655891691285_dp, &
& 0.125580369464904624634694299223940100197616_dp,0.055668567116173666482753720442548578728500_dp/)

! Gauss 15 Kronrod 31
 real(dp),public,parameter :: Kron31N(31) = (/-0.998002298693397060285172840152271209073410_dp, &
& -0.987992518020485428489565718586612581146970_dp,-0.967739075679139134257347978784337225283360_dp, &
& -0.937273392400705904307758947710209471244000_dp,-0.897264532344081900882509656454495882831780_dp, &
& -0.848206583410427216200648320774216851366260_dp,-0.790418501442465932967649294817947346862140_dp, &
& -0.724417731360170047416186054613938009630900_dp,-0.650996741297416970533735895313274692546948_dp, &
& -0.570972172608538847537226737253910641238390_dp,-0.485081863640239680693655740232350612866339_dp, &
& -0.394151347077563369897207370981045468362750_dp,-0.299180007153168812166780024266388962661603_dp, &
& -0.201194093997434522300628303394596207812836_dp,-0.101142066918717499027074231447392338787451_dp, &
&  0.000000000000000000000000000000000000000000_dp, 0.101142066918717499027074231447392338787451_dp, &
&  0.201194093997434522300628303394596207812836_dp, 0.299180007153168812166780024266388962661603_dp, &
&  0.394151347077563369897207370981045468362750_dp, 0.485081863640239680693655740232350612866339_dp, &
&  0.570972172608538847537226737253910641238390_dp, 0.650996741297416970533735895313274692546948_dp, &
&  0.724417731360170047416186054613938009630900_dp, 0.790418501442465932967649294817947346862140_dp, &
&  0.848206583410427216200648320774216851366260_dp, 0.897264532344081900882509656454495882831780_dp, &
&  0.937273392400705904307758947710209471244000_dp, 0.967739075679139134257347978784337225283360_dp, &
&  0.987992518020485428489565718586612581146970_dp, 0.998002298693397060285172840152271209073410_dp/)

 real(dp),public,parameter :: Kron31W(31) = (/0.0053774798729233489877920514301276498183100_dp, &
& 0.0150079473293161225383747630758072680946390_dp, 0.0254608473267153201868740010196533593972700_dp, &
& 0.0353463607913758462220379484783600481226300_dp, 0.0445897513247648766082272993732796902232570_dp, &
& 0.0534815246909280872653431472394302967715500_dp, 0.0620095678006706402851392309608029321904000_dp, &
& 0.0698541213187282587095200770991474757860450_dp, 0.0768496807577203788944327774826590067221100_dp, &
& 0.0830805028231330210382892472861037896015540_dp, 0.0885644430562117706472754436937743032122700_dp, &
& 0.0931265981708253212254868727473457185619300_dp, 0.0966427269836236785051799076275893351366570_dp, &
& 0.0991735987217919593323931734846031310595673_dp, 0.1007698455238755950449466626175697219163500_dp, &
& 0.1013300070147915490173747927674925467709270_dp, 0.1007698455238755950449466626175697219163500_dp, &
& 0.0991735987217919593323931734846031310595673_dp, 0.0966427269836236785051799076275893351366570_dp, &
& 0.0931265981708253212254868727473457185619300_dp, 0.0885644430562117706472754436937743032122700_dp, &
& 0.0830805028231330210382892472861037896015540_dp, 0.0768496807577203788944327774826590067221100_dp, &
& 0.0698541213187282587095200770991474757860450_dp, 0.0620095678006706402851392309608029321904000_dp, &
& 0.0534815246909280872653431472394302967715500_dp, 0.0445897513247648766082272993732796902232570_dp, &
& 0.0353463607913758462220379484783600481226300_dp, 0.0254608473267153201868740010196533593972700_dp, &
& 0.0150079473293161225383747630758072680946390_dp, 0.0053774798729233489877920514301276498183100_dp/)

 real(dp),public,parameter :: Gau15W(15) = (/0.030753241996117268354628393577204417721700_dp, &
0.070366047488108124709267416450667338466710_dp, 0.107159220467171935011869546685869303415544_dp, &
0.139570677926154314447804794511028322520850_dp, 0.166269205816993933553200860481208811130900_dp, &
0.186161000015562211026800561866422824506226_dp, 0.198431485327111576456118326443839324818693_dp, &
0.202578241925561272880620199967519314838662_dp, 0.198431485327111576456118326443839324818693_dp, &
0.186161000015562211026800561866422824506226_dp, 0.166269205816993933553200860481208811130900_dp, &
0.139570677926154314447804794511028322520850_dp, 0.107159220467171935011869546685869303415544_dp, &
0.070366047488108124709267416450667338466710_dp, 0.030753241996117268354628393577204417721700_dp/)


! Flags for self-consistent GW calculations used in gw_driver and for parsing the input file.
 integer,public,parameter :: GWSC_one_shot      =1
 integer,public,parameter :: GWSC_only_W        =2
 integer,public,parameter :: GWSC_only_G        =3
 integer,public,parameter :: GWSC_both_G_and_W  =4

! Flags defining the approximation used for the self-energy (used in csigme).
 integer,public,parameter :: SIG_GW_PPM      =0  ! standard GW with PPM
 integer,public,parameter :: SIG_GW_AC       =1  ! standard GW without PPM (analytical continuation)
 integer,public,parameter :: SIG_GW_CD       =2  ! standard GW without PPM (contour deformation)
 integer,public,parameter :: SIG_HF          =5  ! Hartree-Fock calculation
 integer,public,parameter :: SIG_SEX         =6  ! Screened Exchange calculation
 integer,public,parameter :: SIG_COHSEX      =7  ! COHSEX calculation
 integer,public,parameter :: SIG_QPGW_PPM    =8  ! model GW with PPM
 integer,public,parameter :: SIG_QPGW_CD     =9  ! model GW without PPM

 public :: sigma_type_from_key
 public :: sigma_is_herm
 public :: sigma_needs_w
 public :: sigma_needs_ppm
 !public :: sigma_sc_on_wfs
 !public :: sigma_sc_on_ene
 public :: g0g0w

! Private variables
 integer,private,parameter :: STR_LEN=500
!!***

!----------------------------------------------------------------------

!!****t* m_gwdefs/em1params_t
!! NAME
!! em1params_t
!!
!! FUNCTION
!! For the GW part of ABINIT, the  em1params_t structured datatype
!! gather different parameters used to calculate the inverse dielectric matrices.
!!
!! SOURCE

 type,public :: em1params_t

!scalars
  integer :: awtr                   ! If 1 the Adler-Wiser expression for Chi_0 is evaluated
                                    !  taking advantage of time-reversal symmetry
  integer :: gwcalctyp              ! Calculation type (see input variable)
  integer :: gwcomp                 ! 1 if extrapolar technique is used. 0 otherwise.
  integer :: inclvkb                ! Integer flag related to the evaluation of the commutator for q-->0
  integer :: spmeth                 ! Method used to approximate the delta function in the expression for Im Chi_0
  integer :: nI                     ! Number of components (rows) in the chi0 matrix.
  integer :: nJ                     ! Number of components (columns) in the chi0 matrix.
  integer :: npwvec                 ! Max between npwe and npwwfn, used to pass the dimension of arrays e.g gvec
  integer :: npwwfn                 ! Number of planewaves for wavefunctions
  integer :: npwe                   ! Number of planewaves for $\tilde \epsilon$
  integer :: npwepG0                ! Number of planewaves in the enlarged sphere G-G0, to account for umklapp G0 vectors
  integer :: nbnds                  ! Number of bands used to evaluate $\tilde \epsilon$
  integer :: nkibz                  ! Number of k-points in the IBZ
  integer :: nsppol                 ! 1 for spin unpolarized, 2 for collinear spin polarized
  integer :: nqcalc                 ! Number of q-points that are calculated (subset of qibz)
  integer :: nqibz                  ! Number of q-points in the IBZ
  integer :: nqlwl                  ! Number of directions to analyze the non analytical behavior for q-->0
  integer :: nomega                 ! Number of frequencies where evaluate $\tilde \epsilon (\omega)$
  integer :: nomegaer,nomegaei      ! Number of real and imaginary frequencies, respectively
  integer :: nomegaec               ! Number of frequencies on a grid in the complex plane nomegaec = nomegaei*(nomegaer-1)
  integer :: nomegasf               ! Number of frequencies used for the spectral function
  integer :: symchi                 ! 0 ==> do not use symmetries to reduce the k-points summed over in chi0
                                    ! 1 ==> take advantage of point group symmetries as well as time-reversal

  real(dp) :: gwencomp              ! Extrapolar energy used if gwcomp==1.
  real(dp) :: omegaermin            ! Minimum real frequency used in the contour deformation method
  real(dp) :: omegaermax            ! Maximum real frequency used in the contour deformation method
  real(dp) :: mbpt_sciss              ! Scissor energy used in chi0
  real(dp) :: spsmear               ! Smearing of the delta in case of spmeth==2
  real(dp) :: zcut                  ! Small imaginary shift to avoid poles in chi0

  logical :: analytic_continuation  ! if true calculate chi0 only along the imaginary axis
  logical :: contour_deformation    ! if true calculated chi0 both along the real and the imaginary axis
  logical :: plasmon_pole_model     ! if true a plasmonpole model is used (only 1 or 2 frequencies are calculated)

!arrays
  integer :: mG0(3)
  ! For each reduced direction gives the max G0 component to account for umklapp processes

  real(dp),allocatable :: qcalc(:,:)
  ! qcalc(3,nqcalc)
  ! q-points that are explicitely calculated (subset of qibz).

  real(dp),allocatable :: qibz(:,:)
  ! qibz(3,nqibz)
  ! q-points in the IBZ.

  real(dp),allocatable :: qlwl(:,:)
  ! qlwl(3,nqlwl)
  ! q-points used for the long-wavelength limit.

  real(dp),allocatable :: omegasf(:)
  ! omegasf(nomegasf)
  ! real frequencies used to calculate the imaginary part of chi0.

  complex(dpc),allocatable :: omega(:)
  ! omega(nomegasf)
  ! real and imaginary frequencies in chi0,epsilon and epsilonm1.

 end type em1params_t

 public :: em1params_free
!!***

 type,public :: sigij_col_t
   integer :: size1
   integer,allocatable :: bidx(:)
 end type sigij_col_t

 type,public :: sigijtab_t
   type(sigij_col_t),allocatable :: col(:)
 end type sigijtab_t

 public :: sigijtab_free

!----------------------------------------------------------------------

!!****t* m_gwdefs/sigparams_t
!! NAME
!! sigparams_t
!!
!! FUNCTION
!! For the GW part of ABINIT, the sigparams_t structured datatype
!! gather different parameters that characterize the calculation of the matrix
!! elements of the self-energy operator.
!!
!! SOURCE

 type,public :: sigparams_t

  integer :: gwcalctyp                   ! Calculation type

  integer :: gwgamma                     ! If 1 include vertex correction (GWGamma)
  integer :: gwcomp                      ! 1 if the extrapolar technique is used.

  integer :: minbdgw,maxbdgw             ! Minimum and maximum band index (considering the spin) defining
                                         ! The set of bands where GW corrections are evaluated

  integer :: mG0(3)                      ! For each reduced direction gives the max G0 component
                                         ! to account for umklapp processes

  integer :: npwvec                      ! Max betwenn npwe and npwwfn, used to pass the dimension of arrays e.g gvec
  integer :: npwwfn                      ! No. of planewaves for wavefunctions
  integer :: npwx                        ! No. of planewaves for $\Sigma_x$
  integer :: npwc                        ! No. of planewaves for $\Sigma_c$ and W
  integer :: nbnds                       ! No. of bands summed over.
  integer :: nomegasr                    ! No. of frequencies on the real axis to evaluate the spectral function
  integer :: nomegasrd                   ! No. of frequencies on the real axis to evaluate $\Sigma(E)$
  integer :: nomegasi                    ! No. of frequencies along the imaginary axis for Sigma in case of AC
  integer :: nsig_ab                     ! No. of components in the self-energy operator (1 if nspinor==1, 4 if nspinor==2)
  integer :: nspinor                     ! No. of spinorial components.
  integer :: nsppol                      ! 1 for unpolarized, 2 for spin-polarized calculation
  integer :: nkptgw                      ! No. of k-points where GW corrections have been calculated
  integer :: ppmodel                     ! Integer defining the plasmon pole model used, 0 for None.
  integer :: symsigma                    ! 0 ==> do not use symmetries to reduce the k-points summed over in sigma
                                         ! 1 ==> take advantage of space group symmetries as well as time-reversal
  integer :: use_sigxcore                ! 1 if core contribution to sigma is estimated by using Hartree-Fock

  real(dp) :: deltae                     ! Energy step used to evaluate numerically the derivative of the self energy
                                         ! $\frac{\partial \Re \Sigma(E)}{\partial E_o}$
  real(dp) :: ecutwfn                    ! cutoff energy for the wavefunctions.
  real(dp) :: ecutsigx                   ! cutoff energy for the the exchange parth of Sigma.
  real(dp) :: ecuteps                    ! cutoff energy for W

  real(dp) :: gwencomp                   ! Extrapolar energy used if gwcomp==1.

  real(dp) :: mbpt_sciss                 ! Scissor energy used in G0

  real(dp) :: minomega_r                 ! Minimum real frequency for the evaluation of the spectral function
  real(dp) :: maxomega_r                 ! Maximum real frequency for the evaluation of the spectral function
  real(dp) :: maxomega4sd                ! Maximum displacement around the KS energy where evaluate the diagonal
                                         ! Elements of $ \Sigma(E)$
  real(dp) :: omegasimax                 ! Max omega for Sigma along the imag axis in case of analytic continuation
  real(dp) :: omegasimin                 ! min omega for Sigma along the imag axis in case of analytic continuation

  real(dp) :: sigma_mixing               ! Global factor that multiplies Sigma to give the final matrix element.
                                         ! Usually one, except for the hybrid functionals.

  real(dp) :: zcut                       ! Value of $\delta$ used to avoid the divergences (see related input variable)

  integer,allocatable :: kptgw2bz(:)
  ! kptgw2bz(nkptgw)
  ! For each k-point where GW corrections are calculated, the corresponding index in the BZ.

  integer,allocatable :: minbnd(:,:), maxbnd(:,:)
  ! minbnd(nkptgw,nsppol), maxbnd(nkptgw,nsppol)
  ! For each k-point at which GW corrections are calculated, the min and Max band index considered
  ! (see also input variable dtset%bdgw).

  real(dp),allocatable :: kptgw(:,:)
  ! kptgw(3,nkptgw)
  ! k-points for the GW corrections in reduced coordinates.

  !TODO should be removed, everything should be in Sr%

  complex(dpc),allocatable :: omegasi(:)
  ! omegasi(nomegasi)
  ! Frequencies along the imaginary axis used for the analytical continuation.

  complex(dpc),allocatable :: omega_r(:)
  ! omega_r(nomegasr)
  ! Frequencies used to evaluate the spectral function.

  type(sigijtab_t),allocatable :: Sigcij_tab(:,:)
  ! Sigcij_tab(nkptgw,nsppol)%col(kb)%bidx(ii)  gived the index of the left wavefunction.
  ! in the <i,kgw,s|\Sigma_c|j,kgw,s> matrix elements that has to be calculated in cisgme.
  ! in the case of self-consistent GW on wavefunctions.

  type(sigijtab_t),allocatable :: Sigxij_tab(:,:)
  ! Save as Sigcij_tab but for the Hermitian \Sigma_x where only the upper triangle is needed.

 end type sigparams_t

 public :: sigparams_free
!!***

CONTAINS  !==============================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/em1params_free
!! NAME
!! em1params_free
!!
!! FUNCTION
!!  Free dynamic memory allocated in the structure.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_screening_driver
!!
!! CHILDREN
!!      sigijtab_free
!!
!! SOURCE

subroutine em1params_free(Ep)

!Arguments ------------------------------------
!scalars
 type(em1params_t),intent(inout) :: Ep

! *************************************************************************

 !@em1params_t

!real
 ABI_SFREE(Ep%qcalc)
 ABI_SFREE(Ep%qibz)
 ABI_SFREE(Ep%qlwl)
 ABI_SFREE(Ep%omegasf)
!complex
 ABI_SFREE(Ep%omega)

end subroutine em1params_free
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/sigijtab_free
!! NAME
!! sigijtab_free
!!
!! FUNCTION
!!   deallocate all memory in a sigijtab_t datatype.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwdefs,m_sigma_driver
!!
!! CHILDREN
!!      sigijtab_free
!!
!! SOURCE

subroutine sigijtab_free(Sigijtab)

!Arguments ------------------------------------
!scalars
 type(sigijtab_t),intent(inout) :: Sigijtab(:,:)

!Local variables
 integer :: ii,jj,kk,ilow,iup
! *************************************************************************

 !@sigijtab_t
  do jj=1,SIZE(Sigijtab,DIM=2)
    do ii=1,SIZE(Sigijtab,DIM=1)

      ilow=LBOUND(Sigijtab(ii,jj)%col,DIM=1)
      iup =UBOUND(Sigijtab(ii,jj)%col,DIM=1)
      do kk=ilow,iup
        ABI_FREE(Sigijtab(ii,jj)%col(kk)%bidx)
      end do
      ABI_FREE(Sigijtab(ii,jj)%col)

    end do
  end do

end subroutine sigijtab_free
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/sigparams_free
!! NAME
!! sigparams_free
!!
!! FUNCTION
!!  Free dynamic memory allocated in the structure.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_sigma_driver
!!
!! CHILDREN
!!      sigijtab_free
!!
!! SOURCE

subroutine sigparams_free(Sigp)

!Arguments ------------------------------------
!scalars
 type(sigparams_t),intent(inout) :: Sigp

! *************************************************************************

 !@sigparams_t

!integer
 ABI_SFREE(Sigp%kptgw2bz)
 ABI_SFREE(Sigp%minbnd)
 ABI_SFREE(Sigp%maxbnd)
!real
 ABI_FREE(Sigp%kptgw)
!complex
 ABI_SFREE(Sigp%omegasi)
 ABI_SFREE(Sigp%omega_r)

!types
 if (allocated(Sigp%Sigcij_tab)) then
   call sigijtab_free(Sigp%Sigcij_tab)
   ABI_FREE(Sigp%Sigcij_tab)
 end if

 if (allocated(Sigp%Sigxij_tab)) then
   call sigijtab_free(Sigp%Sigxij_tab)
   ABI_FREE(Sigp%Sigxij_tab)
 end if

end subroutine sigparams_free
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/sigma_type_from_key
!! NAME
!! sigma_type_from_key
!!
!! FUNCTION
!!  Return a string definining the particular approximation used for the self-energy.
!!  Stops if the key is not in the list of allowed possibilities.
!!
!! INPUTS
!!  key=Integer
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function sigma_type_from_key(key) result(sigma_type)

 integer,intent(in) :: key
 character(len=STR_LEN) :: sigma_type

!Local variables ------------------------------
!scalars
 character(len=500) :: msg

!************************************************************************

 sigma_type = "None"
 if (key==SIG_GW_PPM  )  sigma_type = ' standard GW with PPM'
 if (key==SIG_GW_AC   )  sigma_type = ' standard GW without PPM (analytical continuation)'
 if (key==SIG_GW_CD   )  sigma_type = ' standard GW without PPM (contour deformation)'
 if (key==SIG_HF      )  sigma_type = ' Hartree-Fock calculation'
 if (key==SIG_SEX     )  sigma_type = ' Screened Exchange calculation'
 if (key==SIG_COHSEX  )  sigma_type = ' COHSEX calculation'
 if (key==SIG_QPGW_PPM)  sigma_type = ' model GW with PPM'
 if (key==SIG_QPGW_CD )  sigma_type = ' model GW without PPM'

 if (sigma_type == "None") then
   write(msg,'(a,i0)')" Unknown value for key= ",key
   MSG_ERROR(msg)
 end if

end function sigma_type_from_key
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/sigma_is_herm
!! NAME
!! sigma_is_herm
!!
!! FUNCTION
!!  Return .TRUE. if the approximated self-energy is hermitian.
!!
!! INPUTS
!!  Sigp<sigparams_t>=datatype gathering data and info on the self-energy run.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure logical function sigma_is_herm(Sigp)

!Arguments ------------------------------
!scalars
 type(sigparams_t),intent(in) :: Sigp

!Local variables ------------------------------
 integer :: mod10

!************************************************************************

 mod10=MOD(Sigp%gwcalctyp,10)
 sigma_is_herm = ANY(mod10 == [SIG_HF, SIG_SEX, SIG_COHSEX])

end function sigma_is_herm
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/sigma_needs_w
!! NAME
!! sigma_needs_w
!!
!! FUNCTION
!!  Return .TRUE. if self-energy requires the screened interaction W.
!!  For example HF does not need the SCR file.
!!
!! INPUTS
!!  Sigp<sigparams_t>=datatype gathering data and info on the self-energy run.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure logical function sigma_needs_w(Sigp)

!Arguments ------------------------------
!scalars
 type(sigparams_t),intent(in) :: Sigp

!Local variables ------------------------------
 integer :: mod10

!************************************************************************

 mod10=MOD(Sigp%gwcalctyp,10)
 sigma_needs_w = (mod10/=SIG_HF)

end function sigma_needs_w
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/sigma_needs_ppm
!! NAME
!! sigma_needs_ppm
!!
!! FUNCTION
!!  Return .TRUE. if the self-energy run requires a plasmon-pole model.
!!
!! INPUTS
!!  Sigp<sigparams_t>=datatype gathering data and info on the self-energy run.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure logical function sigma_needs_ppm(Sigp)

!Arguments ------------------------------
!scalars
 type(sigparams_t),intent(in) :: Sigp

!Local variables ------------------------------
 integer :: mod10

!************************************************************************

 mod10=MOD(Sigp%gwcalctyp,10)
 sigma_needs_ppm = ( ANY(mod10 == (/SIG_GW_PPM, SIG_QPGW_PPM/)) .or. &
&                   Sigp%gwcomp==1                                   &
&                  )

end function sigma_needs_ppm
!!***

!----------------------------------------------------------------------

!!****f* m_gwdefs/g0g0w
!! NAME
!! g0g0w
!!
!! FUNCTION
!!  Calculates the frequency-dependent part of the RPA polarizability G0G0.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

function g0g0w(omega,numerator,delta_ene,zcut,TOL_W0,opt_poles)

!Arguments ------------------------------------
!scalars
 integer,intent(in):: opt_poles
 real(dp),intent(in) :: TOL_W0,delta_ene,numerator,zcut
 complex(dpc) :: g0g0w
 complex(dpc),intent(in) :: omega

!Local variables ------------------------------
!scalars
 real(dp) :: sgn
 character(len=500) :: msg

!************************************************************************

 if (delta_ene**2>tol14) then
   sgn=SIGN(1.0_dp,delta_ene)
   !
   if (opt_poles == 2) then ! Resonant and anti-resonant contributions.
     if (DABS(REAL(omega))>TOL_W0) then
       g0g0w =  numerator / (omega + delta_ene - j_dpc*sgn*zcut)&
&              -numerator / (omega - delta_ene + j_dpc*sgn*zcut)
     else
       g0g0w =  numerator / (omega + delta_ene)&
&              -numerator / (omega - delta_ene)
     end if

   else if (opt_poles == 1) then ! Only resonant contribution is included.
     if (DABS(REAL(omega))>TOL_W0) then
       g0g0w =  numerator / (omega + delta_ene - j_dpc*sgn*zcut)
     else
       g0g0w =  numerator / (omega + delta_ene)
     end if

   else
     write(msg,'(a,i0)')" Wrong value for opt_poles: ",opt_poles
     MSG_ERROR(msg)
   end if ! opt_poles

 else ! delta_ene**2<tol14
   g0g0w = czero
 end if

end function g0g0w
!!***

END MODULE m_gwdefs
!!***
