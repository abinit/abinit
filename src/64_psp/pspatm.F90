!{\src2tex{textfont=tt}}
!!****f* ABINIT/pspatm
!! NAME
!! pspatm
!!
!! FUNCTION
!! Open atomic pseudopotential data file for a given atom,
!! read the three first lines, make some checks, then
!! call appropriate subroutine for the reading of
!! V(r) and wf R(r) data for each angular momentum, and subsequent
!! Fourier and Bessel function transforms for local and nonlocal potentials.
!! Close psp file at end.
!!
!! Handles pseudopotential files produced by (pspcod=1 or 4) Teter code,
!! or from the Goedecker-Teter-Hutter paper (pspcod=2),
!! or from the Hartwigsen-Goedecker-Hutter paper (pspcod=3 or 10)
!! or "Phoney pseudopotentials" (Hamman grid in real space) (pspcod=5)
!! or "Troullier-Martins pseudopotentials" from the FHI (pspcod=6)
!! or "XML format" (pspcod=9)
!! or "UPF PWSCF format" (pspcod=11)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, FrD, AF, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dq= spacing of the q-grid
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | ixc=exchange-correlation choice from main routine data file
!!   | pawxcdev=choice of XC development in PAW formalism
!!   | usexcnhat_orig=choice for use of nhat in Vxc in PAW formalism
!!   | xclevel= XC functional level
!!  ipsp=id in the array of the currently read pseudo.
!!
!! OUTPUT
!!  ekb(dimekb)=
!!    ->NORM-CONSERVING PSPS ONLY (pspcod/=7):
!!      (Real) Kleinman-Bylander energies (hartree)
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for number of basis functions (l,n) (dimekb=lnmax)
!!             If any, spin-orbit components begin at l=mpsang+1
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$(hartree)
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  vlspl(mqgrid_vl,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  ffspl(mqgrid_ff,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector; if any, spin-orbit components begin at l=mpsang+1
!!  xcccrc=XC core correction cutoff radius (bohr) from psp file
!!  xccc1d(n1xccc*(1-usepaw),6)=1D core charge function and five derivatives, from psp file (used in NC only)
!!  nctab=<nctab_t>
!!    has_tvale=True if the pseudo provides the valence density (used in NC only)
!!    tvalespl(mqgrid_vl(1-usepaw),2)=the pseudo valence density and 2nd derivative in reciprocal space on a regular grid 
!!                                     (used in NC only)
!!
!! SIDE EFFECTS
!! Input/Output :
!!  psps <type(pseudopotential_type)>=at output, values depending on the read
!!                                    pseudo are set.
!!   | dimekb(IN)=dimension of ekb (see module defs_psp.f)
!!   | filpsp(IN)=name of formatted external file containing atomic psp data.
!!   | lmnmax(IN)=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!   |           =if useylm=0, max number of (l,n)   comp. over all type of psps
!!   | lnmax(IN)=max. number of (l,n) components over all type of psps
!!   |           angular momentum of nonlocal pseudopotential
!!   | mpsang(IN)= 1+maximum angular momentum for nonlocal pseudopotentials
!!   | mpssoang(IN)= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!   | mqgrid_ff(IN)=dimension of q (or G) grid for nl form factors (array ffspl)
!!   | mqgrid_vl(IN)=dimension of q (or G) grid for Vloc (array vlspl)
!!   | n1xccc(IN)=dimension of xccc1d ; 0 if no XC core correction is used
!!   | optnlxccc(IN)=option for nl XC core correction
!!   | positron(IN)=0 if electron GS calculation
!!   |              1 if positron GS calculation
!!   |              2 if electron GS calculation in presence of the positron
!!   | pspso(INOUT)=spin-orbit characteristics, govern the content of ffspl and ekb
!!   |          if =0 : this input requires NO spin-orbit characteristics of the psp
!!   |          if =2 : this input requires HGH characteristics of the psp
!!   |          if =3 : this input requires HFN characteristics of the psp
!!   |          if =1 : this input will be changed at output to 1, 2, 3, according
!!   |                  to the intrinsic characteristics of the psp file
!!   | qgrid_ff(mqgrid_ff)(IN)=values of q on grid from 0 to qmax (bohr^-1) for nl form factors
!!   | qgrid_vl(mqgrid_vl)(IN)=values of q on grid from 0 to qmax (bohr^-1) for Vloc
!!   | usepaw(IN)= 0 for non paw calculation; =1 for paw calculation
!!   | useylm(IN)=governs the way the nonlocal operator is to be applied:
!!   |            1=using Ylm, 0=using Legendre polynomials
!!   | vlspl_recipSpace(IN)=.true. if pseudo are expressed in reciprocal space.
!!   | znuclpsp(IN)=atomic number of atom as specified in input file to main routine
!!
!! NOTES
!!  Format expected for the three first lines of pseudopotentials
!!  (1) title (character) line
!!  (2) znucl,zion,pspdat
!!  (3) pspcod,pspxc,lmax,lloc,mmax,r2well
!!
!!  Dimensions of form factors and Vloc q grids must be the same in Norm-Conserving case
!!
!! PARENTS
!!      pspini
!!
!! CHILDREN
!!      nctab_eval_tcorespl,pawpsp_17in,pawpsp_7in,pawpsp_bcast
!!      pawpsp_read_header_xml,pawpsp_read_pawheader,pawpsp_wvl,psp10in,psp1in
!!      psp2in,psp3in,psp5in,psp6in,psp8in,psp9in,psxml2abheader
!!      test_xml_xmlpaw_upf,timab,upf2abinit,wrtout,wvl_descr_psp_fill
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pspatm(dq,dtset,dtfil,ekb,epsatm,ffspl,indlmn,ipsp,pawrad,pawtab,&
&  psps,vlspl,dvlspl,xcccrc,xccc1d,nctab,comm_mpi)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_psps

 use m_io_tools, only : open_file
 use m_pawrad,   only : pawrad_type
 use m_pawtab,   only : pawtab_type
 use m_pawpsp,   only : pawpsp_bcast, pawpsp_read_pawheader, pawpsp_read_header_xml,&
&                       pawpsp_header_type, pawpsp_wvl, pawpsp_7in, pawpsp_17in
 use m_pawxmlps,  only : paw_setup, ipsp2xml
 use m_psxml2ab
#if defined HAVE_DFT_BIGDFT
 use BigDFT_API, only : dictionary, atomic_info, dict_init, dict_free, UNINITIALIZED
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pspatm'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_43_wvl_wrappers
 use interfaces_64_psp, except_this_one => pspatm
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ipsp
 integer, optional,intent(in) :: comm_mpi
 real(dp),intent(in) :: dq
 real(dp),intent(out) :: epsatm,xcccrc
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(pawrad_type),intent(inout) :: pawrad
 type(pawtab_type),intent(inout) :: pawtab
 type(nctab_t),intent(inout) :: nctab
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 integer,intent(out) :: indlmn(6,psps%lmnmax)
 real(dp),intent(out) :: dvlspl(psps%mqgrid_vl,2)
 real(dp),intent(inout) :: ekb(psps%dimekb*(1-psps%usepaw))
 real(dp),intent(inout) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
 real(dp),intent(out) :: vlspl(psps%mqgrid_vl,2)
 real(dp),intent(inout) :: xccc1d(psps%n1xccc*(1-psps%usepaw),6)

!Local variables ---------------------------------------
!scalars
 integer :: ii,il,ilmn,iln,iln0,lloc,lmax,me,mmax
 integer :: paral_mode,pspcod,pspdat,pspxc,useupf,usexml,xmlpaw,unt
 real(dp) :: maxrad,qchrg,r2well,zion,znucl
 character(len=500) :: message,errmsg
 character(len=fnlen) :: title
 character(len=fnlen) :: filnam
 type(pawpsp_header_type):: pawpsp_header
!arrays
 integer,allocatable :: nproj(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: e990(:),e999(:),ekb1(:),ekb2(:),epspsp(:),rcpsp(:)
 real(dp),allocatable :: rms(:)
#if defined HAVE_TRIO_PSML
!!  usexml= 0 for non xml ps format ; =1 for xml ps format
 type(pspheader_type) :: psphead
#endif

! ******************************************************************************

!paral_mode defines how we access to the psp file
!  paral_mode=0: all processes access to the file (sequentially)
!  paral_mode=1: only proc. 0 access to the file and then broadcast
 paral_mode=0
 if (present(comm_mpi)) then
   if (psps%usepaw==1.and.xmpi_comm_size(comm_mpi)>1) paral_mode=1
 end if
 me=0;if (paral_mode==1) me=xmpi_comm_rank(comm_mpi)

 if (paral_mode == 1) then
   ABI_CHECK(psps%usepaw==1, "paral_mode==1 is only compatible with PAW, see call to pawpsp_bcast below")
 end if

 nctab%has_tvale = .False.; nctab%has_tcore = .False. 

 if (me==0) then
!  Dimensions of form factors and Vloc q grids must be the same in Norm-Conserving case
   if (psps%usepaw==0 .and. psps%mqgrid_ff/=psps%mqgrid_vl) then
     write(message, '(a,a,a,a,a)' )&
&     'Dimension of q-grid for nl form factors (mqgrid_ff)',ch10,&
&     'is different from dimension of q-grid for Vloc (mqgrid_vl) !',ch10,&
&     'This is not allowed for norm-conserving psp.'
     MSG_ERROR(message)
   end if

   write(message, '(a,t38,a)' )'- pspatm: opening atomic psp file',trim(psps%filpsp(ipsp))
   call wrtout(ab_out,  message,'COLL')
   call wrtout(std_out,  message,'COLL')
   !write(message, "(2a)")"- md5: ",trim(psps%md5_pseudos(ipsps))

   !  Check if the file pseudopotential file is written in (XML| XML-PAW | UPF)
   call test_xml_xmlpaw_upf(psps%filpsp(ipsp), usexml, xmlpaw, useupf)

!  ----------------------------------------------------------------------------
!  allocate nproj here: can be read in now for UPF
   ABI_ALLOCATE(nproj,(psps%mpssoang))
   nproj(:)=0

   if (usexml /= 1 .and. useupf /= 1) then

!    Open the atomic data file, and read the three first lines
!    These three first lines have a similar format in all allowed psp files

!    Open atomic data file (note: formatted input file)
     if (open_file(psps%filpsp(ipsp), message, unit=tmp_unit, form='formatted', status='old') /= 0) then
       MSG_ERROR(message)
     end if 
     rewind (unit=tmp_unit,err=10,iomsg=errmsg)

!    Read and write some description of file from first line (character data)
     read (tmp_unit,'(a)',err=10,iomsg=errmsg) title
     write(message, '(a,a)' ) '- ',trim(title)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

!    Read and write more data describing psp parameters
     read (tmp_unit,*,err=10,iomsg=errmsg) znucl,zion,pspdat
     write(message,'(a,f9.5,f10.5,2x,i8,t47,a)')'-',znucl,zion,pspdat,'znucl, zion, pspdat'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

     read (tmp_unit,*,err=10,iomsg=errmsg) pspcod,pspxc,lmax,lloc,mmax,r2well
     if(pspxc<0) then
       write(message, '(i5,i8,2i5,i10,f10.5,t47,a)' ) &
&       pspcod,pspxc,lmax,lloc,mmax,r2well,'pspcod,pspxc,lmax,lloc,mmax,r2well'
     else
       write(message, '(4i5,i10,f10.5,t47,a)' ) &
&       pspcod,pspxc,lmax,lloc,mmax,r2well,'pspcod,pspxc,lmax,lloc,mmax,r2well'
     end if
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

   else if (usexml == 1 .and. xmlpaw == 0) then

! the following is probably useless - already read in everything in inpspheads
#if defined HAVE_TRIO_PSML
     write(message,'(a,a)') &
&     '- pspatm: Reading pseudopotential header in XML form from ', trim(psps%filpsp(ipsp))
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

     lloc   = 0 ! does this mean s? in psml case the local potential can be different from any l channel
     r2well = 0

     call psxml2abheader( psps%filpsp(ipsp), psphead, 1 )
     znucl = psphead%znuclpsp
     zion = psphead%zionpsp
     pspcod = psphead%pspcod 
     pspxc =  psphead%pspxc
     lmax = psphead%lmax
     
#else
     write(message,'(a,a)')  &
&     'ABINIT is not compiled with XML support for reading this type of pseudopotential ', &
&     trim(psps%filpsp(ipsp))
     MSG_BUG(message)
     
#endif
! END useless
   else if (usexml == 1 .and. xmlpaw == 1) then
     write(message,'(a,a)')  &
&     '- pspatm : Reading pseudopotential header in XML form from ', trim(psps%filpsp(ipsp))
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

!PENDING: These routines will be added to libpaw:
!    Return header informations
     call pawpsp_read_header_xml(lloc,lmax,pspcod,&
&     pspxc,paw_setup(ipsp2xml(ipsp)),r2well,zion,znucl)
!    Fill in pawpsp_header object:
     call pawpsp_read_pawheader(pawpsp_header%basis_size,&
&     lmax,pawpsp_header%lmn_size,&
&     pawpsp_header%l_size,pawpsp_header%mesh_size,&
&     pawpsp_header%pawver,paw_setup(ipsp2xml(ipsp)),&
&     pawpsp_header%rpaw,pawpsp_header%rshp,pawpsp_header%shape_type)

   else if (useupf == 1) then
     if (psps%usepaw /= 0) then
       MSG_ERROR("UPF format not allowed with PAW (USPP part not read yet)")
     end if

     pspcod = 11
     r2well = 0

!    should initialize znucl,zion,pspxc,lmax,lloc,mmax
     call upf2abinit (psps%filpsp(ipsp), znucl, zion, pspxc, lmax, lloc, mmax, &
&     psps, epsatm, xcccrc, indlmn, ekb, ffspl, nproj, vlspl, xccc1d)

   else
     MSG_ERROR("You should not be here! erroneous type or pseudopotential input")
   end if

!  ------------------------------------------------------------------------------
!  Check data for consistency against main routine input

!  Does required spin-orbit characteristics agree with format
!  TODO: in case of pspcod 5 (phoney) and 8 (oncvpsp) this is not specific enough.
!  they can be non-SOC as well.
!  HGH is ok - can always turn SOC on or off.
!  write(std_out,*) pspso
   if((pspcod/=3).and.(pspcod/=5).and.(pspcod/=8).and.(pspcod/=10))then
!    If pspso requires internal characteristics, set it to 1 for non-HGH psps
     if(psps%pspso(ipsp)==1) psps%pspso(ipsp)=0
     if(psps%pspso(ipsp)/=0)then
       write(message, '(3a,i0,3a)' )&
&       'Pseudopotential file cannot give spin-orbit characteristics,',ch10,&
&       'while pspso(itypat)= ',psps%pspso(ipsp),'.',ch10,&
&       'Action: check your pseudopotential and input files for consistency.'
       MSG_ERROR(message)
     end if
   end if

!  Does nuclear charge znuclpsp agree with psp input znucl
   !write(std_out,*)znucl,ipsp,psps%znuclpsp(ipsp)
   !MGNAG: v5[66] gives NAG in %znuclpsp if -nan
   if (abs(psps%znuclpsp(ipsp)-znucl)>tol8) then
     write(message, '(a,f10.5,2a,f10.5,5a)' )&
&     'Pseudopotential file znucl=',znucl,ch10,&
&     'does not equal input znuclpsp=',psps%znuclpsp(ipsp),' better than 1e-08 .',ch10,&
&     'znucl is read from the psp file in pspatm, while',ch10,&
&     'znuclpsp is read in iofn2.'
     MSG_BUG(message)
   end if

!  Is the highest angular momentum within limits?
!  Recall mpsang is 1+highest l for nonlocal correction.
!  Nonlocal corrections for s, p, d, and f are supported.
   if (lmax+1>psps%mpsang) then
     write(message, '(a,i0,a,i0,a,a)' )&
&     'input lmax+1= ',lmax+1,' exceeds mpsang= ',psps%mpsang,ch10,&
&     'indicates input lmax too large for dimensions.'
     MSG_BUG(message)
   end if

!  Check several choices for ixc against pspxc
!  ixc is from ABINIT code; pspxc is from atomic psp file
   if (dtset%ixc==0) then
     MSG_WARNING('Note that input ixc=0 => no xc is being used.')
   else if(dtset%ixc/=pspxc) then
     write(message, '(a,i0,3a,i0,9a)' )&
&     'Pseudopotential file pspxc= ',pspxc,',',ch10,&
&     'not equal to input ixc= ',dtset%ixc,'.',ch10,&
&     'These parameters must agree to get the same xc ',ch10,&
&     'in ABINIT code as in psp construction.',ch10,&
&     'Action: check psp design or input file.',ch10,&
&     'Assume experienced user. Execution will continue.'
     MSG_WARNING(message)
   end if

   if (lloc>lmax .and. pspcod/=4 .and. pspcod/=8 .and. pspcod/=10) then
     write(message, '(a,2i12,a,a,a,a)' )&
&     'lloc,lmax=',lloc,lmax,ch10,&
&     'chosen l of local psp exceeds range from input data.',ch10,&
&     'Action: check pseudopotential input file.'
     MSG_ERROR(message)
   end if

!  Does the pspcod agree with type of calculation (paw or not)?
   if (((pspcod/=7.and.pspcod/=17).and.psps%usepaw==1).or.((pspcod==7.or.pspcod==17).and.psps%usepaw==0)) then
     write(message, '(a,i2,a,a,i0,a)' )&
&     'In reading atomic psp file, finds pspcod= ',pspcod,ch10,&
&     'This is not an allowed value with usepaw= ',psps%usepaw,'.'
     MSG_BUG(message)
   end if

   if (.not.psps%vlspl_recipSpace .and. &
&   (pspcod /= 2 .and. pspcod /= 3 .and. pspcod /= 10 .and. pspcod /= 7)) then
!    The following "if" statement can substitute the one just before once libBigDFT
!    has been upgraded to include pspcod 10
!    if (.not.psps%vlspl_recipSpace .and. (pspcod /= 2 .and. pspcod /= 3 .and. pspcod /= 10)) then
     write(message, '(a,i2,a,a)' )&
&     'In reading atomic psp file, finds pspcod=',pspcod,ch10,&
&     'This is not an allowed value with real space computation.'
     MSG_BUG(message)
   end if

!  MJV 16/6/2009 added pspcod 11 for upf format
   if( pspcod<1 .or. (pspcod>11.and.pspcod/=17) ) then
     write(message, '(a,i0,4a)' )&
&     'In reading atomic psp file, finds pspcod= ',pspcod,ch10,&
&     'This is not an allowed value. Allowed values are 1 to 11 .',ch10,&
&     'Action: check pseudopotential input file.'
     MSG_ERROR(message)
   end if

!  -----------------------------------------------------------------------
!  Set various terms to 0 in case not defined below
   ABI_ALLOCATE(e990,(psps%mpssoang))
   ABI_ALLOCATE(e999,(psps%mpssoang))
   ABI_ALLOCATE(rcpsp,(psps%mpssoang))
   ABI_ALLOCATE(rms,(psps%mpssoang))
   ABI_ALLOCATE(epspsp,(psps%mpssoang))
   ABI_ALLOCATE(ekb1,(psps%mpssoang))
   ABI_ALLOCATE(ekb2,(psps%mpssoang))
   e990(:)=zero ;e999(:)=zero
   rcpsp(:)=zero;rms(:)=zero
   ekb1(:)=zero ;ekb2(:)=zero
   epspsp(:)=zero
   qchrg=0

!  ----------------------------------------------------------------------
   if(pspcod==1 .or. pspcod==4)then

!    Teter pseudopotential (pspcod=1 or 4)
     call psp1in(dq,ekb,ekb1,ekb2,epsatm,epspsp,&
&     e990,e999,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,&
&     mmax,psps%mpsang,psps%mqgrid_ff,nproj,psps%n1xccc,pspcod,qchrg,psps%qgrid_ff,&
&     rcpsp,rms,psps%useylm,vlspl,xcccrc,xccc1d,zion,psps%znuclpsp(ipsp))

   else if (pspcod==2)then

!    GTH pseudopotential
     call psp2in(dtset,ekb,epsatm,ffspl,indlmn,ipsp,lmax,nproj,psps,vlspl,dvlspl,zion)
     xccc1d(:,:)=0.0d0 ; qchrg=0.0d0 ; xcccrc=0.0d0

   else if (pspcod==3)then

!    HGH pseudopotential
     call psp3in(dtset,ekb,epsatm,ffspl,indlmn,ipsp,lmax,nproj,psps, psps%pspso(ipsp), &
&     vlspl,zion)
     xccc1d(:,:)=0.0d0 ; qchrg=0.0d0 ; xcccrc=0.0d0

   else if (pspcod==5)then

!    Old phoney pseudopotentials
     call psp5in(ekb,ekb1,ekb2,epsatm,epspsp,&
&     e990,e999,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,&
&     mmax,psps%mpsang,psps%mpssoang,psps%mqgrid_ff,nproj,psps%n1xccc,psps%pspso(ipsp),qchrg,psps%qgrid_ff,&
&     rcpsp,rms,psps%useylm,vlspl,xcccrc,xccc1d,zion,psps%znuclpsp(ipsp))

   else if (pspcod==6)then

!    FHI pseudopotentials
     call psp6in(ekb,epsatm,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,mmax,&
&     psps%mpsang,psps%mqgrid_ff,nproj,psps%n1xccc,psps%optnlxccc,psps%positron,qchrg,psps%qgrid_ff,psps%useylm,vlspl,&
&     xcccrc,xccc1d,zion,psps%znuclpsp(ipsp))

   else if (pspcod==7)then
!    PAW "pseudopotentials"
     call pawpsp_7in(epsatm,ffspl,dtset%icoulomb,dtset%ixc,&
&     lmax,psps%lnmax,mmax,psps%mqgrid_ff,psps%mqgrid_vl,&
&     pawrad,pawtab,dtset%pawxcdev,psps%qgrid_ff,psps%qgrid_vl,&
&     dtset%usewvl,dtset%usexcnhat_orig,vlspl,xcccrc,dtset%xclevel,&
&     dtset%xc_denpos,zion,psps%znuclpsp(ipsp))

   else if (pspcod==8)then

!    DRH pseudopotentials
     call psp8in(ekb,epsatm,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,mmax,&
&     psps%mpsang,psps%mpssoang,psps%mqgrid_ff,psps%mqgrid_vl,nproj,psps%n1xccc,psps%pspso(ipsp),&
&     qchrg,psps%qgrid_ff,psps%qgrid_vl,psps%useylm,vlspl,xcccrc,xccc1d,zion,psps%znuclpsp(ipsp),nctab)

   else if (pspcod==9)then

#if defined HAVE_TRIO_PSML
     call psp9in(psps%filpsp(ipsp),ekb,epsatm,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,mmax,&
&     psps%mpsang,psps%mpssoang,psps%mqgrid_ff,nproj,psps%n1xccc, &
&     psps%pspso(ipsp),qchrg,psps%qgrid_ff,psps%useylm,vlspl,&
&     xcccrc,xccc1d,zion,psps%znuclpsp(ipsp))
#else
     write(message,'(2a)')  &
&     'ABINIT is not compiled with XML support for reading this type of pseudopotential ', &
&     trim(psps%filpsp(ipsp))
     MSG_BUG(message)
#endif

   else if (pspcod==10)then

!    HGH pseudopotential, full h/k matrix read
     call psp10in(dtset,ekb,epsatm,ffspl,indlmn,ipsp,lmax,nproj,psps, psps%pspso(ipsp), &
&     vlspl,zion)
     xccc1d(:,:)=0.0d0 ; qchrg=0.0d0 ; xcccrc=0.0d0

!    NB for pspcod 11 the reading has already been done above.
   else if (pspcod==17)then
!    PAW XML "pseudopotentials"
     call pawpsp_17in(epsatm,ffspl,dtset%icoulomb,ipsp,dtset%ixc,lmax,&
&     psps%lnmax,mmax,psps%mqgrid_ff,psps%mqgrid_vl,pawpsp_header,pawrad,pawtab,&
&     dtset%pawxcdev,psps%qgrid_ff,psps%qgrid_vl,dtset%usewvl,&
&     dtset%usexcnhat_orig,vlspl,xcccrc,&
&     dtset%xclevel,dtset%xc_denpos,zion,psps%znuclpsp(ipsp))
   end if

   close (unit=tmp_unit)

!  ----------------------------------------------------------------------
   if (pspcod==2 .or. pspcod==3 .or. pspcod==10)then
     write(message, '(a,a,a,a,a,a,a,a,a,a)' )ch10,&
&     ' pspatm : COMMENT -',ch10,&
&     '  the projectors are not normalized,',ch10,&
&     '  so that the KB energies are not consistent with ',ch10,&
&     '  definition in PRB44, 8503 (1991). ',ch10,&
&     '  However, this does not influence the results obtained hereafter.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
!    The following lines are added to keep backward compatibilty
     maxrad=zero
#if defined HAVE_DFT_BIGDFT
     do ii=1,size(psps%gth_params%psppar,1)-1 ! psppar first dim begins at 0
       if (psps%gth_params%psppar(ii,0,ipsp)/=zero) maxrad=max(maxrad,psps%gth_params%psppar(ii,0,ipsp))
     end do
     if (abs(maxrad)<=tol12) then
       psps%gth_params%radii_cf(ipsp,3)=zero
     else
       psps%gth_params%radii_cf(ipsp,3)=max( &
&       min(dtset%wvl_crmult*psps%gth_params%radii_cf(ipsp,1),15._dp*maxrad)/dtset%wvl_frmult, &
&       psps%gth_params%radii_cf(ipsp,2))
     end if
#endif
   end if

   if (pspcod/=7.and.pspcod/=17) then
     write(message, '(a,f14.8,a,a)' ) ' pspatm: epsatm=',epsatm,ch10,&
&     '         --- l  ekb(1:nproj) -->'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     iln0=0
     do ilmn=1,psps%lmnmax
       iln=indlmn(5,ilmn)
       if (iln>iln0) then
         il=indlmn(1,ilmn)
         if (indlmn(6,ilmn)==1) then
           iln0=iln0+nproj(il+1)
           write(message, '(13x,i1,4f12.6)' ) il,(ekb(iln+ii),ii=0,nproj(il+1)-1)
         else
           iln0=iln0+nproj(il+psps%mpsang)
           write(message, '(2x,a,i1,4f12.6)' ) 'spin-orbit ',il,(ekb(iln+ii),ii=0,nproj(il+psps%mpsang)-1)
         end if
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
       end if
     end do
   end if

   ! NC: Evalute spline-fit of the model core charge in reciprocal space.
   ! TODO: Be careful, because we will be using the PAW part in which tcore is always avaiable!
   ! Should add a test with 2 NC pseudos: one with NLCC and the other without!
   if (psps%usepaw == 0) then
     call nctab_eval_tcorespl(nctab, psps%n1xccc, xcccrc, xccc1d, psps%mqgrid_vl, psps%qgrid_vl)
   end if

   write(message,'(3a)') ' pspatm: atomic psp has been read ',' and splines computed',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   ABI_DEALLOCATE(e990)
   ABI_DEALLOCATE(e999)
   ABI_DEALLOCATE(rcpsp)
   ABI_DEALLOCATE(rms)
   ABI_DEALLOCATE(ekb1)
   ABI_DEALLOCATE(ekb2)
   ABI_DEALLOCATE(epspsp)
   ABI_DEALLOCATE(nproj)

   if (dtset%prtvol > 9 .and. psps%usepaw==0 .and. psps%lmnmax>3) then

     write (filnam, '(a,i0,a)') trim(dtfil%fnameabo_pspdata), ipsp, ".dat"
     if (open_file(filnam, message, newunit=unt) /= 0) then
       MSG_ERROR(message)
     end if
     write (unt,*) '# Pseudopotential data in reciprocal space as used by ABINIT'
     write (unt,'(a)', ADVANCE='NO') '# index       vlocal   '
     if (psps%lnmax > 0) &
&     write (unt,'(a,I3)', ADVANCE='NO')   '           1st proj(l=', indlmn(1,1)
     if (psps%lnmax > 1) &
&     write (unt,'(a,I3)', ADVANCE='NO')   ')            2nd(l=', indlmn(1,2)
     if (psps%lnmax > 2) &
&     write (unt,'(a,I3,a)', ADVANCE='NO') ')            3rd(l=', indlmn(1,3), ')'
     write (unt,*)

     do ii = 1, psps%mqgrid_vl
       write(unt, '(I5,E24.16)', ADVANCE='NO') ii, vlspl(ii,1)
       if (psps%lnmax > 0) write(unt, '(E24.16)', ADVANCE='NO') ffspl(ii,1,1)
       if (psps%lnmax > 1) write(unt, '(E24.16)', ADVANCE='NO') ffspl(ii,1,2)
       if (psps%lnmax > 2) write(unt, '(E24.16)', ADVANCE='NO') ffspl(ii,1,3)
       write(unt, *)
     end do
     close(unt)

     write (filnam, '(a,i0,a)') trim(dtfil%fnameabo_nlcc_derivs), ipsp, ".dat"
     if (open_file(filnam, message, newunit=unt) /= 0) then
       MSG_ERROR(message)
     end if 
     write (unt,*) '# Non-linear core corrections'
     write (unt,*) '#  r, pseudocharge, 1st, 2nd, 3rd, 4th, 5th derivatives'
     do ii = 1, psps%n1xccc
       write (unt,*) xcccrc*(ii-1)/(psps%n1xccc-1), xccc1d(ii,1), xccc1d(ii,2), &
&       xccc1d(ii,3), xccc1d(ii,4), &
&       xccc1d(ii,5), xccc1d(ii,6)
     end do
     close(unt)
   end if

 end if !me=0

 if (paral_mode==1) then
   call timab(48,1,tsec)
   call pawpsp_bcast(comm_mpi,epsatm,ffspl,pawrad,pawtab,vlspl,xcccrc)
   call timab(48,2,tsec)
 end if

 if (psps%usepaw==1) then
   indlmn(:,:)=0
   indlmn(1:6,1:pawtab%lmn_size)=pawtab%indlmn(1:6,1:pawtab%lmn_size)
 end if

!--------------------------------------------------------------------
!WVL+PAW: 
 if (dtset%usepaw==1 .and. (dtset%icoulomb /= 0 .or. dtset%usewvl==1)) then
#if defined HAVE_DFT_BIGDFT
   psps%gth_params%psppar(:,:,ipsp) = UNINITIALIZED(1._dp)
   psps%gth_params%radii_cf(ipsp,:) = UNINITIALIZED(1._dp)
   call wvl_descr_psp_fill(psps%gth_params, ipsp, psps%pspxc(1), int(psps%zionpsp(ipsp)), int(psps%znuclpsp(ipsp)), 0)
#endif

!  The following lines are added to keep backward compatibilty
   maxrad=zero
#if defined HAVE_DFT_BIGDFT
   do ii=1,size(psps%gth_params%psppar,1)-1 ! psppar first dim begins at 0
     if (psps%gth_params%psppar(ii,0,ipsp)/=zero) maxrad=max(maxrad,psps%gth_params%psppar(ii,0,ipsp))
   end do
   if (abs(maxrad)<=tol12) then
!== MT COMMENT
!    Damien wants to activate this (in order to directly compare to bigDFT):
     psps%gth_params%radii_cf(ipsp,3)= psps%gth_params%radii_cf(ipsp,2)
!    But, this changes strongly file references.
!    So, I keep this, waiting for Tonatiuh s validation
     psps%gth_params%radii_cf(ipsp,3)= &
&     (psps%gth_params%radii_cf(ipsp,1)+psps%gth_params%radii_cf(ipsp,2))*half
!== MT COMMENT
   else
     psps%gth_params%radii_cf(ipsp,3)=max( &
&     min(dtset%wvl_crmult*psps%gth_params%radii_cf(ipsp,1),15._dp*maxrad)/dtset%wvl_frmult, &
&     psps%gth_params%radii_cf(ipsp,2))
   end if
   if(present(comm_mpi)) then
     call pawpsp_wvl(psps%filpsp(ipsp),pawrad,pawtab,dtset%usewvl,dtset%wvl_ngauss,comm_mpi)
   else
     call pawpsp_wvl(psps%filpsp(ipsp),pawrad,pawtab,dtset%usewvl,dtset%wvl_ngauss)
   end if
#endif
 end if
 
!end of WVL+PAW section
!----------------------------------------------------

 return 

 ! Handle IO error
 10 continue
 MSG_ERROR(errmsg)

end subroutine pspatm
!!***
