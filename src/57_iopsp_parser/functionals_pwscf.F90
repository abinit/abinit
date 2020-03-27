!!****m* ABINIT/functionals_pwscf.F90
!! NAME
!!  functionals_pwscf
!!
!! FUNCTION
!! This module contains data defining the DFT functional in use
!! and a number of functions and subroutines to manage them.
!! Data are PRIVATE and are accessed and set only by function calls.
!! Basic drivers to compute XC quantities are also included.
!!
!! imported into abinit by MJV 20/6/2009
!    changed module name
!    removed DP references
!    removed functions for explicit calculation of xc (E,V,fxc)
!!
!! COPYRIGHT
!! Copyright (C) 2004 PWSCF group
!! This file is distributed under the terms of the
!! GNU General Public License. See the file `License'
!! in the root directory of the present distribution,
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!
!-------------------------------------------------------------------

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module funct_pwscf
!-------------------------------------------------------------------
!  
!  setting routines:   set_dft_from_name (previously which_dft)
!                      set_dft_from_indices
!                      enforce_input_dft
!                      start_exx
!                      stop_exx
!  retrive functions:  get_dft_name
!                      get_iexch
!                      get_icorr
!                      get_igcx
!                      get_igcc
!                      get_exx_fraction
!                      dft_name
!                      write_dft_name
!  logical functions:  dft_is_gradient
!                      dft_is_meta
!                      dft_is_hybrid
!                      exx_is_active
!
!  XC computation drivers: xc, xc_spin, gcxc, gcx_spin, gcc_spin, gcc_spin_more
!  derivatives of XC computation drivers: dmxc, dmxc_spin, dmxc_nc
!
  use flib_pwscf

  IMPLICIT NONE
  PRIVATE
  SAVE
  ! subroutines/functions managing dft name and indices
  PUBLIC  :: set_dft_from_indices, set_dft_from_name
  PUBLIC  :: enforce_input_dft, write_dft_name, dft_name
  PUBLIC  :: get_dft_name, get_iexch, get_icorr, get_igcx, get_igcc
  PUBLIC  :: dft_is_gradient, dft_is_meta, dft_is_hybrid
  ! additional subroutines/functions for hybrid functionale
  PUBLIC  :: start_exx, stop_exx, get_exx_fraction, exx_is_active
  !
  ! PRIVATE variables defining the DFT functional
  !
  PRIVATE :: dft, iexch, icorr, igcx, igcc
  PRIVATE :: discard_input_dft
  PRIVATE :: isgradient, ismeta, ishybrid
  PRIVATE :: exx_fraction, exx_started
  !PRIVATE :: dft_shortname
  !
  character (len=50) :: dft = 'not set'
  !character (len=4)  :: dft_shortname = ' '
  !
  ! dft is the exchange-correlation functional, described by
  ! any nonconflicting combination of the following keywords
  ! (case-insensitive):
  !
  ! Exchange:    "nox"    none                           iexch=0
  !              "sla"    Slater (alpha=2/3)             iexch=1 (default)
  !              "sl1"    Slater (alpha=1.0)             iexch=2
  !              "rxc"    Relativistic Slater            iexch=3
  !              "oep"    Optimized Effective Potential  iexch=4
  !              "hf"     Hartree-Fock                   iexch=5
  !              "pb0x"   PBE0                           iexch=6
  !
  ! Correlation: "noc"    none                           icorr=0
  !              "pz"     Perdew-Zunger                  icorr=1 (default)
  !              "vwn"    Vosko-Wilk-Nusair              icorr=2
  !              "lyp"    Lee-Yang-Parr                  icorr=3
  !              "pw"     Perdew-Wang                    icorr=4
  !              "wig"    Wigner                         icorr=5
  !              "hl"     Hedin-Lunqvist                 icorr=6
  !              "obz"    Ortiz-Ballone form for PZ      icorr=7
  !              "obw"    Ortiz-Ballone form for PW      icorr=8
  !              "gl"     Gunnarson-Lunqvist             icorr=9
  !
  ! Gradient Correction on Exchange:
  !              "nogx"   none                           igcx =0 (default)
  !              "b88"    Becke88 (beta=0.0042)          igcx =1
  !              "ggx"    Perdew-Wang 91                 igcx =2
  !              "pbx"    Perdew-Burke-Ernzenhof exch    igcx =3
  !              "rpb"    revised PBE by Zhang-Yang      igcx =4
  !              "hcth"   Cambridge exch, Handy et al    igcx =5
  !              "optx"   Handy's exchange functional    igcx =6
  !              "meta"   meta-gga                       igcx =7
  !              "pb0x"   PBE0                           igcx =8
  !
  ! Gradient Correction on Correlation:
  !              "nogc"   none                           igcc =0 (default)
  !              "p86"    Perdew86                       igcc =1
  !              "ggc"    Perdew-Wang 91 corr.           igcc =2
  !              "blyp"   Lee-Yang-Parr                  igcc =3
  !              "pbc"    Perdew-Burke-Ernzenhof corr    igcc =4
  !              "hcth"   Cambridge corr, Handy et al    igcc =5
  !              "meta"   meta-gga                       igcc =6
  !
  ! Special cases (dft_shortnames):
  !              "bp"   = "b88+p86"         = Becke-Perdew grad.corr.
  !              "pw91" = "pw +ggx+ggc"     = PW91 (aka GGA)
  !              "blyp" = "sla+b88+lyp+blyp"= BLYP
  !              "pbe"  = "sla+pw+pbx+pbc"  = PBE
  !              "revpbe"="sla+pw+rpb+pbc"  = revPBE (Zhang-Yang)
  !              "hcth" = "nox+noc+hcth+hcth"=HCTH/120
  !              "olyp" = "nox+lyp+optx+blyp" !!! UNTESTED !!!
  !
  ! References:
  !              pz      J.P.Perdew and A.Zunger, PRB 23, 5048 (1981) [[cite:Perdew1981]]
  !              vwn     S.H.Vosko, L.Wilk, M.Nusair, Can.J.Phys. 58,1200(1980) [[cite:Vosko1980]]
  !              wig     E.P.Wigner, Trans. Faraday Soc. 34, 67 (1938) [[cite:Wigner1938]]
  !              hl      L.Hedin and B.I.Lundqvist, J. Phys. C4, 2064 (1971) [[cite:Hedin1971]]
  !              gl      O.Gunnarsson and B.I.Lundqvist, PRB 13, 4274 (1976) [[cite:Gunnarsson1976]]
  !              pw      J.P.Perdew and Y.Wang, PRB 45, 13244 (1992) [[cite:Perdew1992a]]
  !              obpz    G.Ortiz and P.Ballone, PRB 50, 1391 (1994) [[cite:Ortiz1994]]
  !              obpw    as above
  !              b88     A.D.Becke, PRA 38, 3098 (1988) [[cite:Becke1988]]
  !              p86     J.P.Perdew, PRB 33, 8822 (1986) [[cite:Perdew1986]] 
  !              pbe     J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996) [[cite:Perdew1996]]
  !              pw91    J.P.Perdew and Y. Wang, PRB 46, 6671 (1992) [[cite:Perdew1992]]
  !              blyp    C.Lee, W.Yang, R.G.Parr, PRB 37, 785 (1988) [[cite:Lee1988]]
  !              hcth    Hamprecht et al, JCP 109, 6264 (1998) [[cite:Hamprecht1998]]
  !              olyp    Handy and Cohen, JCP 116, 5411 (2002) [[cite:Handy2002]]
  !              revPBE  Zhang and Yang, PRL 80, 890 (1998) [[cite:Zhang1998]]
  !              oep

  integer, parameter:: notset = -1
  !
  integer :: iexch = notset
  integer :: icorr = notset
  integer :: igcx  = notset
  integer :: igcc  = notset
  real(8):: exx_fraction = 0.0d0
  logical :: isgradient  = .false.
  logical :: ismeta      = .false.
  logical :: ishybrid    = .false.
  logical :: exx_started = .false.

  logical :: discard_input_dft = .false.
  !
  ! internal indices for exchange-correlation
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlation
  !
  !    ismeta: .TRUE. if gradient correction is of meta-gga type
  !    ishybrid: .TRUE. if the xc finctional is an HF+DFT hybrid like
  !              PBE0 or B3LYP or HF itself
  !
  ! see comments above and routine "set_dft_from_name" below 
  !
  ! data
  integer :: nxc, ncc, ngcx, ngcc
  parameter (nxc = 7, ncc =10, ngcx = 9, ngcc = 7)
  character (len=4) :: exc, corr
  character (len=4) :: gradx, gradc
  dimension exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0: ngcc)

  data exc / 'NOX', 'SLA', 'SL1', 'RXC', 'OEP', 'HF', 'PB0X', 'B3LP' /
  data corr / 'NOC', 'PZ', 'VWN', 'LYP', 'PW', 'WIG', 'HL', 'OBZ', &
              'OBW', 'GL' , 'B3LP' /
  data gradx / 'NOGX', 'B88', 'GGX', 'PBX',  'RPB', 'HCTH', 'OPTX', 'META', 'PB0X', 'B3LP'  /
  data gradc / 'NOGC', 'P86', 'GGC', 'BLYP', 'PBC', 'HCTH', 'META', 'B3LP' /

CONTAINS
!!***

!!****f* functionals_pwscf/set_dft_from_name
!!
!! NAME
!! set_dft_from_name
!!
!! FUNCTION
!! translates a string containing the exchange-correlation name
!! into internal indices iexch, icorr, igcx, igcc
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  !-----------------------------------------------------------------------
  subroutine set_dft_from_name( dft_ )
  !-----------------------------------------------------------------------

    use flib_pwscf
    implicit none
    ! input
    character(len=*)               :: dft_
    ! local
    integer :: len, l, i
    character (len=50):: dftout
    !
    !
    ! if 
    !
    if ( discard_input_dft ) return
    !
    ! convert to uppercase
    len = len_trim(dft_)
    dftout = ' '
    do l = 1, len
       dftout (l:l) = capital (dft_(l:l) )
    enddo

    !  exchange
    iexch = notset
    do i = 0, nxc
       if (matches (exc (i), dftout) ) then
         call set_dft_value (iexch, i)
       end if
    enddo

    !  correlation
    icorr = notset
    do i = 0, ncc
       if (matches (corr (i), dftout) ) then
         call set_dft_value (icorr, i)
       end if
    enddo

    !  gradient correction, exchange
    igcx = notset
    do i = 0, ngcx
       if (matches (gradx (i), dftout) ) then
         call set_dft_value (igcx, i)
       end if
    enddo

    !  gradient correction, correlation
    igcc = notset
    do i = 0, ngcc
       if (matches (gradc (i), dftout) ) then
         call set_dft_value (igcc, i)
       end if
    enddo

    ! special case : BLYP => B88 for gradient correction on exchange
    if (matches ('BLYP', dftout) ) then
      call set_dft_value (igcx, 1)
    end if

    ! special case : revPBE
    if (matches ('REVPBE', dftout) ) then
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 4)
       call set_dft_value (igcc, 4)
    else if (matches('RPBE',dftout)) then
         call errore('set_dft_from_name', &
     &   'RPBE (Hammer-Hansen-Norskov) not implemented (revPBE is)',1)
   else if (matches ('PBE0', dftout) ) then
    ! special case : PBE0
       call set_dft_value (iexch,6)
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 8)
       call set_dft_value (igcc, 4)
   else if (matches ('PBE', dftout) ) then
    ! special case : PBE
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 3)
       call set_dft_value (igcc, 4)
   endif

    if (matches ('PBC', dftout) ) then
    ! special case : PBC  = PW + PBC 
       call set_dft_value (icorr,4)
       call set_dft_value (igcc, 4)
    endif

    ! special case : BP = B88 + P86
    if (matches ('BP', dftout) ) then
       call set_dft_value (igcx, 1)
       call set_dft_value (igcc, 1)
    endif

    ! special case : PW91 = GGX + GGC
    if (matches ('PW91', dftout) ) then
       call set_dft_value (igcx, 2)
       call set_dft_value (igcc, 2)
    endif

    ! special case : HCTH already contains LDA exchange and correlation

    if (matches('HCTH',dftout)) then
       call set_dft_value(iexch,0)
       call set_dft_value(icorr,0)
    end if

    ! special case : OPTX already contains LDA exchange
     
    if (matches('OPTX',dftout)) then
       call set_dft_value(iexch,0)
    end if

    ! special case : OLYP = OPTX + LYP

    if (matches('OLYP',dftout)) then
       call set_dft_value(iexch,0)
       call set_dft_value(icorr,3)
       call set_dft_value(igcx,6)
       call set_dft_value(igcc,3)
    end if
    !
    ! ... special case : TPSS meta-GGA Exc
    !
    IF ( matches( 'TPSS', dftout ) ) THEN
       !
       CALL set_dft_value( iexch, 1 )
       CALL set_dft_value( icorr, 4 )
       CALL set_dft_value( igcx,  7 )
       CALL set_dft_value( igcc,  6 )
       !
    END IF
    !
    ! ... special cases : OEP and HF need not GC part (nor LDA...)
    !                     and include no correlation by default
    !
    IF ( matches( 'OEP', dftout ) .OR. matches( 'HF', dftout )) THEN
       !
       CALL set_dft_value( igcx,  0 )
       if (icorr == notset) then
         call set_dft_value (icorr, 0)
       end if
       !
    END IF


    if (igcx == 6) &
         call errore('set_dft_from_name','OPTX untested! please test',-igcx)
    ! Default value: Slater exchange
    if (iexch == notset) then
      call set_dft_value (iexch, 1)
    end if

    ! Default value: Perdew-Zunger correlation
    if (icorr == notset) then
      call set_dft_value (icorr, 1)
    end if

    ! Default value: no gradient correction on exchange
    if (igcx == notset) then
      call set_dft_value (igcx, 0)
    end if

    ! Default value: no gradient correction on correlation
    if (igcc == notset) then
      call set_dft_value (igcc, 0)
    end if

    dft = dftout

    dftout = exc (iexch) //'-'//corr (icorr) //'-'//gradx (igcx) //'-' &
         &//gradc (igcc)

    call set_auxiliary_flags

    return
  end subroutine set_dft_from_name
!!***

!!****f* functionals_pwscf/set_auxiliary_flags
!!
!! NAME
!! set_auxiliary_flags
!!
!! FUNCTION
!! set logical flags describing the complexity of the xc functional
!! define the fraction of exact exchange used by hybrid fuctionals
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

  !-----------------------------------------------------------------------
  subroutine set_auxiliary_flags
  !-----------------------------------------------------------------------

    use flib_pwscf
    isgradient =  (igcx > 0) .or. (igcc > 0) 
    ismeta     =  (igcx == 7) .or. (igcx == 6 )

    ! PBE0
    IF ( iexch==6 .or. igcx==8 ) exx_fraction = 0.25d0
    ! HF or OEP
    IF ( iexch==4 .or. iexch==5 ) exx_fraction = 1.0d0
    !B3LYP
    IF ( matches( 'B3LP',dft ) ) exx_fraction = 0.2d0
    ishybrid = ( exx_fraction /= 0.0d0 )

    return
  end subroutine set_auxiliary_flags
!!***

!!****f* functionals_pwscf/set_dft_value
!!
!! NAME
!! set_dft_value
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  !-----------------------------------------------------------------------
  subroutine set_dft_value (m, i)
  !-----------------------------------------------------------------------
    use flib_pwscf
    implicit none
    integer :: m, i
    ! local

    if ( m /= notset .and. m /= i) &
         call errore ('set_dft_value', 'two conflicting matching values', 1)
    m = i
    return

  end subroutine set_dft_value
!!***

!!****f* functionals_pwscf/enforce_input_dft
!!
!! NAME
!! enforce_input_dft
!!
!! FUNCTION
!! translates a string containing the exchange-correlation name
!! into internal indices and force any subsequent call to set_dft_from_name
!! to return without changing them
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

  !-----------------------------------------------------------------------
  subroutine enforce_input_dft (dft_)
  !
    use defs_basis, only : std_out,std_out_default
    use flib_pwscf
    implicit none
    ! input
    character(len=*) :: dft_
    ! data

     call set_dft_from_name (dft_)
     if (dft == 'not set') then
       call errore('enforce_input_dft','cannot fix unset dft',1)
     end if
     discard_input_dft = .true.

     write(std_out,'(/,5x,a)') "!!! XC functional enforced from input :"
     call write_dft_name
     write(std_out,'(5x,a)') "!!! Any further DFT definition will be discarded"
     write(std_out,'(5x,a)') "!!! Please, verify this is what you really want !"

     return
  end subroutine enforce_input_dft
!!***

!!****f* functionals_pwscf/start_exx
!!
!! NAME
!! start_exx
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

  subroutine start_exx 

    use flib_pwscf
     if (.not. ishybrid) &
        call errore('start_exx','dft is not hybrid, wrong call',1)
     exx_started = .true.
  end subroutine start_exx
!!***

!!****f* functionals_pwscf/stop_exx
!!
!! NAME
!! stop_exx
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

  !-----------------------------------------------------------------------
  subroutine stop_exx 

    use flib_pwscf
     if (.not. ishybrid) &
        call errore('stop_exx','dft is not hybrid, wrong call',1)
     exx_started = .false.
  end subroutine stop_exx
!!***

!!****f* functionals_pwscf/exx_is_active
!!
!! NAME
!! exx_is_active
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

  function exx_is_active ()

     logical exx_is_active
     exx_is_active = exx_started
  end function exx_is_active
!!***

!!****f* functionals_pwscf/get_iexch
!!
!! NAME
!! get_iexch
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  !-----------------------------------------------------------------------
  function get_iexch ()

     integer get_iexch
     get_iexch = iexch
     return
  end function get_iexch
!!***

!!****f* functionals_pwscf/get_icorr
!!
!! NAME
!! get_icorr
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  !-----------------------------------------------------------------------
  function get_icorr ()

     integer get_icorr
     get_icorr = icorr
     return
  end function get_icorr
!!***

!!****f* functionals_pwscf/get_igcx
!!
!! NAME
!! get_igcx
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  !-----------------------------------------------------------------------
  function get_igcx ()

     integer get_igcx
     get_igcx = igcx
     return
  end function get_igcx
!!***

!!****f* functionals_pwscf/get_igcc
!!
!! NAME
!! get_igcc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  !-----------------------------------------------------------------------
  function get_igcc ()

     integer get_igcc
     get_igcc = igcc
     return
  end function get_igcc
!!***

!!****f* functionals_pwscf/get_exx_fraction
!!
!! NAME
!! get_exx_fraction
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  !-----------------------------------------------------------------------
  function get_exx_fraction ()

     real(8):: get_exx_fraction
     get_exx_fraction = exx_fraction
     return
  end function get_exx_fraction
!!***

!!****f* functionals_pwscf/get_dft_name
!!
!! NAME
!! get_dft_name
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  !-----------------------------------------------------------------------
  function get_dft_name ()

     character (len=50) :: get_dft_name
     get_dft_name = dft
     return
  end function get_dft_name
!!***

!!****f* functionals_pwscf/dft_is_gradient
!!
!! NAME
!! dft_is_gradient
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  !-----------------------------------------------------------------------
  function dft_is_gradient ()

     logical :: dft_is_gradient
     dft_is_gradient = isgradient
     return
  end function dft_is_gradient
!!***

!!****f* functionals_pwscf/dft_is_meta
!!
!! NAME
!! dft_is_meta
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  !-----------------------------------------------------------------------
  function dft_is_meta ()

     logical :: dft_is_meta
     dft_is_meta = ismeta
     return
  end function dft_is_meta
!!***

!!****f* functionals_pwscf/dft_is_hybrid
!!
!! NAME
!! dft_is_hybrid
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  !-----------------------------------------------------------------------
  function dft_is_hybrid ()

     logical :: dft_is_hybrid
     dft_is_hybrid = ishybrid
     return
  end function dft_is_hybrid
!!***

!!****f* functionals_pwscf/set_dft_from_indices
!!
!! NAME
!! set_dft_from_indices
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  !-----------------------------------------------------------------------
  subroutine set_dft_from_indices(iexch_,icorr_,igcx_,igcc_)

    use defs_basis, only : std_out,std_out_default
    use flib_pwscf
    implicit none
     integer :: iexch_, icorr_, igcx_, igcc_
     if ( discard_input_dft ) return
     if (iexch == notset) iexch = iexch_
     if (iexch /= iexch_) then
        write(std_out,*) iexch, iexch_
        call errore('set_dft',' conflicting values for iexch',1)
     end if
     if (icorr == notset) icorr = icorr_
     if (icorr /= icorr_) then
        write(std_out,*) icorr, icorr_
        call errore('set_dft',' conflicting values for icorr',1)
     end if
     if (igcx  == notset) igcx = igcx_
     if (igcx /= igcx_) then
        write(std_out,*) igcx, igcx_
        call errore('set_dft',' conflicting values for igcx',1)
     end if
     if (igcc  == notset) igcc = igcc_
     if (igcc /= igcc_) then
        write(std_out,*) igcc, igcc_
        call errore('set_dft',' conflicting values for igcc',1)
     end if
     dft = exc (iexch) //'-'//corr (icorr) //'-'//gradx (igcx) //'-' &
           &//gradc (igcc)
     ! write(std_out,'(a)') dft
     call set_auxiliary_flags
     return
  end subroutine set_dft_from_indices
!!***

!!****f* functionals_pwscf/dft_name
!!
!! NAME
!! dft_name
!!
!! FUNCTION
!! convert the four indices iexch, icorr, igcx, igcc
!! into user-readable strings
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
  !---------------------------------------------------------------------
  subroutine dft_name(iexch_, icorr_, igcx_, igcc_, longname_, shortname_)
  !---------------------------------------------------------------------
  implicit none
  integer iexch_, icorr_, igcx_, igcc_
  character (len=4) :: shortname_
  character (len=20):: longname_
  !
  if (iexch_==1.and.igcx_==0.and.igcc_==0) then
     shortname_ = corr(icorr_)
  else if (iexch_==1.and.icorr_==3.and.igcx_==1.and.igcc_==3) then
     shortname_ = 'BLYP'
  else if (iexch_==1.and.icorr_==1.and.igcx_==1.and.igcc_==0) then
     shortname_ = 'B88'
  else if (iexch_==1.and.icorr_==1.and.igcx_==1.and.igcc_==1) then
     shortname_ = 'BP'
  else if (iexch_==1.and.icorr_==4.and.igcx_==2.and.igcc_==2) then
     shortname_ = 'PW91'
  else if (iexch_==1.and.icorr_==4.and.igcx_==3.and.igcc_==4) then
     shortname_ = 'PBE'
  else if (iexch_==6.and.icorr_==4.and.igcx_==8.and.igcc_==4) then
     shortname_ = 'PBE0'
  else
     shortname_ = ' '
  end if
  write(longname_,'(4a5)') exc(iexch_),corr(icorr_),gradx(igcx_),gradc(igcc_)
  
  return
end subroutine dft_name
!!***

!!****f* functionals_pwscf/write_dft_name
!!
!! NAME
!! write_dft_name
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_dft_name
!-----------------------------------------------------------------------
   use defs_basis, only : std_out,std_out_default
   implicit none

   !write(std_out,'(5X,"Exchange-correlation      = ",A, &
   !     &  " (",4I1,")")') TRIM( dft ), iexch, icorr, igcx, igcc

   write(std_out,'(5X,a,A,a,4I1,a)') "Exchange-correlation      = ", TRIM( dft ), " (", iexch, icorr, igcx, igcc, ")"
   return
end subroutine write_dft_name


end module funct_pwscf
!!***
