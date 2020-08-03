!!****m* ABINIT/m_timana
!! NAME
!!  m_timana
!!
!! FUNCTION
!! Analyse the timing, and print in unit ab_out. Some discussion of the
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (XG, GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_timana

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_xomp

 use m_time,        only : time_accu, timab, TIMER_SIZE
 use defs_abitypes, only : MPI_type

 implicit none

 private
!!***

 public :: timana
!!***

contains
!!***

!!****f* ABINIT/timana
!! NAME
!! timana
!!
!! FUNCTION
!! Analyse the timing, and print in unit ab_out. Some discussion of the
!! number of calls to different routines is also provided, as comments,
!! at the end of the routine, as well as, in the single dataset mode (ndtset<2),
!! a detailed analysis of the time-consuming routines.
!!
!! INPUTS
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell.
!!  nband(nkpt*nsppol)=number of bands at each k point, for each polarization
!!  ndtset=number of datasets
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nkpt=number of k points
!!  npwtot(nkpt)=number of planewaves in basis at this k point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  timopt= if >0, write short analysis, if <0, write full analysis
!!          if timopt>=2, or timopt==-2 do not time the timer
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! *) One can suppress the cpu timer call in timein.f, if line 315 of the present routine is uncommented.
!!
!! *) The number of fourwf and nonlop calls can be computed as follows, in the
!!    groud-state case, with no reading of wavefunctions (irdwfk==0 and the like),
!!    and iscf>0 :
!!
!!    1) For fourwf.f
!!
!!    In each cgwf call, there will be
!!    1 call (isign=+1 and -1) for the first gradient calculation,
!!    and iline calls for the line minimizations,
!!    minus the number of ffts skipped because some wfs are sufficiently converged
!!    (there is a counter for that, see the log file)
!!
!!    There are nband*nkpt*(nstep+2) calls to cgwf presently, where the
!!    (nstep+2) comes from the number of the presence of 2 nonscf loops
!!    in the first 2 steps.
!!    Thus, the number of fourwf calls in cgwf is
!!    nband*nkpt*(nstep+2)*(1+iline) - nskip_fourwf_in_cgwf
!!
!!    To compute the density (either in vtowfk or in vtorho - by a mkrho call - )
!!    at each step, there will be nband*nkpt one-way calls,
!!    minus the number of bands skipped because the occupation number
!!    is too small (smaller than 1.0d-14). There is another counter for that.
!!    Thus, the number of fourwf calls for the density is
!!    nband*nkpt*nstep - nskip_fourwf_for_density
!!
!!    For example, for Si with nline=3, nkpt=2, nband=4, nstep=10, and supposing
!!    no fourwf calls are skipped, there will be
!!    at most 4*2*12=96 calls to cgwf, with 4 two-way fft,
!!    that is 384 two-way ffts,
!!    and 4*2*10=80 one-way ffts to make the density.
!!    Altogether 464-nskip one-way ffts at most.
!!
!!    2) For nonlop.f
!!
!!    Presently, there are three different types of call to nonlop :
!!    for energy and gradient wrt wavefunctions (choice=1), for forces (choice=2),
!!    and for stresses (choice=3).
!!
!!    In each cgwf call, there will be one nonlop call for two fourwf calls
!!    (independently of the number of skipped fourwf calls, since
!!    nonlop is also skipped then). These are the only calls with choice=1.
!!    Thus the number will be
!!    nband*nkpt*(nstep+2)*(1+iline) - nskip_fourwf_in_cgwf
!!
!!    The number of choice=2 nonlop calls is equal to the number of fourwf calls
!!    to make the density, that is
!!    nband*nkpt*nstep - nskip_fourwf_for_density
!!
!!    The number of choice=8 calls is equal to the number of occupied bands
!!    at the end of the calculation :
!!    nband(occupied)*nkpt
!!    The number of bands skipped then is not counted.
!!
!!    NOTE : the number of fourwf calls is equal to
!!    the # of nonlop (choice=1) calls + the # of nonlop (choice=2) calls
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      timab,time_accu,wrtout,xmpi_sum
!!
!! SOURCE

subroutine timana(mpi_enreg,natom,nband,ndtset,nfft,nkpt,npwtot,nsppol,timopt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ndtset,nfft,nkpt,nsppol,timopt
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: nband(nkpt*nsppol),npwtot(nkpt)

!Local variables-------------------------------
!scalars
 integer :: aslot,bslot,cslot,flag_count,flag_write,ierr,ii,ikpt,ipart
 integer :: ilist,isort,islot,isppol,itim,itimab,ltimab,maxii,me
 integer :: npart,nlist,nothers,nproc,nthreads,return_ncount
 integer(i8b) :: npwmean,npwnbdmean
 integer :: spaceworld,temp_list,totcount,tslot,utimab,ount
 real(dp) :: cpunm,lflops,other_cpu,other_wal,percent_limit,subcpu,subwal,timab_cpu,timab_wall,wallnm
 character(len=500) :: msg
!arrays
 integer(i8b) :: basic(TIMER_SIZE),ndata(TIMER_SIZE),tslots(TIMER_SIZE)
 integer :: ncount(TIMER_SIZE)
 integer,allocatable :: list(:)
 real(dp) :: ftimes(2,TIMER_SIZE),ftsec(2),mflops(TIMER_SIZE),nflops(TIMER_SIZE),times(2,TIMER_SIZE),tsec(2),my_tsec(2)
 character(len=32) :: names(-1:TIMER_SIZE),entry_name
 character(len=*),parameter :: format01040 ="('- ',a24,f12.3,f6.1,f11.3,f6.1,i15,16x,f7.2,1x,f10.2)"
 character(len=*),parameter :: format01041 ="('- ',a24,f12.3,f6.1,f11.3,f6.1,i15,3x,g12.3,1x,f7.2,1x,f10.2)"
 character(len=*),parameter :: format01042 ="('- ',a24,f12.3,f6.1,f11.3,g12.3,i15)"
 character(len=*),parameter :: format01045 ="('-',i3,a19,f15.3,f6.1,f11.3,f6.1)"
 !character(len=*),parameter ::  format01200 ="('- subtotal     ',f15.3,f6.1,f11.3,f6.1)"

! *************************************************************************

 01200 format('- subtotal             ',f15.3,f6.1,f11.3,f6.1,31x,f7.2,1x,f10.2)
 01201 format(/,'- subtotal             ',f15.3,f6.1,f11.3,f6.1,31x,f7.2,1x,f10.2)

 ount = ab_out

 call timab(49,1,tsec)

!The means are computed as integers, for later compatibility
 npwmean=0; npwnbdmean=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     npwmean=npwmean+npwtot(ikpt)
     npwnbdmean=npwnbdmean+npwtot(ikpt)*nband(ikpt+(isppol-1)*nkpt)
   end do
 end do

 ! initialize ftime, valgrind complains on line 832 = sum up of all Gflops
 ftimes=zero

 npwmean=dble(npwmean)/dble(nkpt*nsppol)
 npwnbdmean=dble(npwnbdmean)/dble(nkpt*nsppol)

!List of timed subroutines, eventual initialisation of the number of data, and declaration of a slot as being "basic"
!Channels 1 to 299 are for optdriver=0 (GS), 1 (RF) and 2 (Suscep), at random
!Channels 300 to 399 are for optdriver=3 (Screening)
!Channels 400 to 499 are for optdriver=4 (Sigma)
!Channels 500 to 529 are for optdriver=5 (Nonlinear)
!Channels 530 to 549 are for various counters
!Channels 550 to 599 are for PAW
!Channels 600 to 619 are for Recursion Method
!Channels 620 to 639 are for DMFT
!Channels 650 to 699 are for bethe_salpeter code.
!Channels 700 to 799 are for optdriver=0 (again ...)
!Channels 800 to 899 are for the detailed analysis of fourwf
!Channels 900 to 1499 are for optdriver=0 (again ...)
!Channels 1500 to 1519 are for Hartree-Fock.
!Channels 1700 to 1747 are for GWLS.
!Channels 1520 and beyond are not yet attributed.

 names(1:TIMER_SIZE)='***                             '
!Basic slots are not overlapping. Their sum should cover most of the code.
!WARNING : the slots from 1 to 99 should be avoided in the future ... They are hard to track.
 basic(1:TIMER_SIZE)=0
 names(1)='abinit                          '
 names(5)='ewald                           ' ; basic(5)=1
 names(6)='setsym                          ' ; basic(6)=1
 names(9)='fourdp                          ' ; basic(9)=1 ;    ndata(9)=nfft
 names(10)='hartre                          '
 names(11)='xc:pot/=fourdp                  '; basic(11)=1;    ndata(11)=nfft*nsppol
 names(12)='mkcore                          '; basic(12)=1
 names(13)='mkresi                          '
 names(14)='rwwf                            '; basic(13)=1
 names(15)='pspini                          '; basic(15)=1
 names(16)='mkffnl                          '; basic(16)=1
 names(17)='symrhg(no FFT)                  '; basic(17)=1
 names(19)='inwffil                         '

 names(22)='cgwf                            '
 names(23)='kpgsph                          '; basic(23)=1   ! Actually, should not be basic ... too complicated, too much overlap ...
 names(28)='vtowfk                          '
 names(30)='vtowfk  (afterloop)             '
 names(31)='vtowfk  (1)                     '; basic(31)=1
 names(32)='gstate                          '
 names(33)='gstate->kpgsph                  '
 names(34)='gstate  (2)                     '
 names(35)='gstate(...scfcv)                '
 names(36)='gstate  (3)                     '
 names(37)='stress                          '; basic(37)=1   ! Actually, should not be basic !
 names(38)='ewald2 (+vdw_dftd)              '; basic(38)=1
 names(39)='vtowfk (loop)                   '
 names(40)='cgwf-O(npw)                     '
 names(41)='abinit(1)                       '
 names(42)='abinit(2)                       '; basic(42)=1
 names(43)='indefo+macroin+invars2m         '
 names(44)='abinit(4)                       '
 names(45)='abinit(5)                       '
 names(46)='abinit(6)                       '
 names(47)='ingeo/symgroup                  '
 names(48)='communic.MPI                    '
 names(49)='timana(1)                       '
 names(50)='timing timab                    '; basic(50)=1
 names(51)='total timab                     '
 names(52)='scfcv-scprqt                    '; basic(52)=1
 names(53)='forces-mkcore                   '
 names(54)='scfcv_core(1)                   '
 names(55)='stress-mkcore                   '
 names(56)='scfcv_core-read                 '
 names(57)='rhotov                          '
 names(59)='energy                          '
 names(60)='scfcv_core(etotfor)             '
 names(61)='scfcv_core :synchro             '
 names(62)='kpgio :synchro                  '
 names(63)='mkrho :synchro                  '
 names(64)='strkin:synchro                  '
 names(65)='forstrnps:synchr                '
 names(66)='vtorho:synchro                  '; basic(66)=1
 names(67)='wfsinp:synchro                  '
 names(68)='scfcv_core(mix den - newrho)    '
 names(69)='forces                          '; basic(69)=1 ! Actually, should not be basic !
 names(70)='vtorho(symrhg)                  '
 names(71)='mkrho :MPIrhor                  '
 names(72)='mklocl(2)                       '
 names(73)='status                          '; basic(73)=1
 names(74)='newocc                          '
 names(75)='nonlop(apply)                   '; basic(75)=1; ndata(75)=npwmean*natom
 names(76)='nonlop(forces)                  '; basic(76)=1; ndata(76)=npwmean*natom
 names(77)='nonlop(forstr)                  '; basic(77)=1; ndata(77)=npwmean*natom
 names(78)='nonlop(dyfrnl)                  '
 names(79)='nonlop(ddk)                     '
 names(80)='etotfor/=forces                 '
 names(81)='xc:pot                          ' ! rhotoxc_coll, except the call to hartre.f
 names(82)='xc:fourdp                       '
 names(83)='newvtr/rho(3):io                '; basic(83)=1
 names(84)='suscep                          '
 names(85)='suscep:MPI                      '; basic(85)=1
 names(86)='suscep:synchro                  '; basic(86)=1
 names(87)='suskXX:loop(1)                  '
 names(88)='suskXX:loop(2)                  '
 names(89)='suscep:other                    '
 names(90)='dielmt                          '; basic(90)=1
 names(91)='setvtr                          '
 names(92)='setvtr:mkcore                   '
 names(93)='newvtr                          '
 names(94)='newrho                          '
 names(95)='tddft                           '
 names(96)='dieltcel                        '; basic(96)=1
 names(97)='nonlop(total)                   '
 names(98)='getghc-other                    '; basic(98)=1

 names(101)='dfpt_nstdy                      '
 names(102)='dfpt_nstwf                      '
 names(108)='dfpt_vtowfk(contrib)            '; basic(108)=1
 names(118)='dfpt_vtorho (1)                 '; basic(118)=1
 names(120)='dfpt_scfcv                      '
 names(121)='dfpt_vtorho                     '
 names(122)='dfpt_cgwf                       '
 names(124)='dfpt_vtorho (1)(2)              '
 names(125)='dfpt_vtorho (2)                 '
 names(126)='dfpt_vtorho-kpt loop            '; basic(126)=1
 names(127)='dfpt_vtorho (4)                 '
 names(128)='dfpt_vtowfk                     '
 names(129)='dfpt_vtorho:MPI                 '; basic(129)=1
 names(130)='dfpt_vtowfk (3)                 '; basic(130)=1
 names(131)='dfpt_vtowfk (1)                 '; basic(131)=1
 names(132)='respfn                          '
 names(133)='respfn(kpgio)                   '
 names(134)='respfn(pspini)                  '
 names(135)='respfn(inwffil)                 '
 names(136)='respfn(frozen)                  '
 names(137)='respfn(dfpt_dyxc1+bef.dfpt_lop) '
 names(138)='respfn(after dfpt_loper)        '
 names(139)='dfpt_vtowfk (loop)              '
 names(140)='dfpt_cgwf-O(npw)                '; basic(140)=1
 names(141)='dfpt_loper                      '
 names(142)='dfpt_loper(kpgio)               '
 names(143)='dfpt_loper(getmpw)              '
 names(144)='dfpt_loper(inwffil)             '
 names(146)='dfpt_loper(outwf)               '
 names(147)='dfpt_loper(eig2tot)             '
 names(148)='eig2tot                         '; basic(148)=1
 names(150)='dfpt_nselt/nstdy/nstpaw         '
 names(152)='dfpt_scfcv-scprqt               '
 names(154)='dfpt_scfcv  (1)                 '; basic(154)=1
 names(157)='dfpt_rhotov                     '
 names(158)='dfpt_newvtr                     '
 names(159)='d2frnl                          '
 names(160)='dfpt_scfcv (6)                  '
 names(161)='dfpt_nstdy:synchro              '; basic(161)=1
 names(166)='dfpt_vtorho:synchro             '; basic(166)=1
 names(181)='dfpt_mkvxc                      '
 names(182)='dfpt_dyxc1                      '; basic(182)=1
!names(184)='dfpt_dyxc1(analysis)            '

 names(191)='invars2                         '; basic(191)=1
 names(192)='inkpts                          '
 names(193)='fresid                          '

 names(195)='getgh1c_setup'; basic(195) = 1
 names(196)='getgh1c'; basic(196) = 1
 names(197)='getgh1c%dfpt_cgwf               '
 names(198)='getgh1c%dfpt_nstwf              '
 names(199)='getgh1c%dfpt_nstpaw             '

 names(200)='getghc                          '
 names(201)='getghc%cgwf                     '
 names(202)='getghc%dfpt_cgwf                '
 names(203)='getghc%mkresi                   '
 names(204)='getghc%kss_ddiago               '
 names(205)='getghc%lobpcgwf                 '
 names(206)='getghc%prep_getghc              '
 names(207)='getghc%other lobpcg             '
 names(208)='getghc%update_mmat              '

 names(210)='projbd                          '; basic(210)=1;    ndata(210)=npwnbdmean
 names(211)='projbd%cgwf                     '
 names(212)='projbd%dfpt_cgwf                '
 names(213)='projbd%dfpt_nstpaw              '
 names(214)='corrmetalwf1%dfpt_vtowfk        '

 names(220)='nonlop%(other)                  '
 names(221)='nonlop%getghc                   '
 names(222)='nonlop%vtowfk                   '
 names(223)='nonlop%energy                   '
 names(224)='nonlop%forstrnps                '
 names(225)='nonlop%dfpt_nstwf               '
 names(226)='nonlop%d2frnl                   '
 names(227)='nonlop%dfpt_cgwf !2             '
 names(228)='nonlop%dfpt_cgwf !5             '
 names(229)='nonlop%outkss                   '
 names(230)='nonlop%vtowfk(rhoij)            '
 names(231)='nonlop%prep_nonl%vtowfk         '
 names(232)='nonlop%prep_nonl%forstrn        '
 names(233)='nonlop%appinvovl                '
 names(234)='nonlop%prep_nonl%energy         '

 names(238)='scfcv_core                      '
 names(239)='scfcv_core(Berry)               '
 names(240)='scfcv_core(iniloop, setvtr  )   '
 names(241)='scfcv_core(loop, PAW)           '
 names(242)='scfcv_core(vtorho(f))           '
 names(243)='scfcv_core(rhotov)              '
 names(244)='scfcv_core(qui loop)            '
 names(245)='scfcv_core(mix pot)             '
 names(246)='scfcv_core(just after scf)      '
 names(247)='scfcv_core(afterscfloop)        '
 names(248)='scfcv_core(outscfcv)            '
 names(249)='scfcv_core(free)                '

 names(250)='afterscfloop                    '
 names(251)='afterscfloop(wvl)               '
 names(252)='afterscfloop(pol/magn)          '
 names(253)='afterscfloop(grad/lapl)         '
 names(254)='afterscfloop(kin.en.den)        '
 names(255)='afterscfloop(elf)               '
 names(256)='afterscfloop(forstr)            '
 names(257)='afterscfloop(final)             '

 names(260)='fourdp%(other)                  '
 names(261)='fourdp%rhotwg%ch                '
 names(262)='fourdp%rhotwg%si                '
 names(263)='fourdp%ckxcldag                 '
 names(264)='fourdp%fftwfn%ch                '
 names(265)='fourdp%fftwfn%si                '
 names(266)='fourdp%rec%rho                  '
 names(267)='fourdp%rec%ek                   '
 names(268)='fourdp%newvtr                   '
 names(269)='fourdp%newrho                   '

 names(270)='rwwf%(other)                    '
 names(271)='rwwf%vtorho                     '
 names(272)='rwwf%initwf(GS)                 '
 names(273)='rwwf%energy                     '
 names(274)='rwwf%wfsinp(GS)                 '
 names(275)='rwwf%mkrho                      '
 names(276)='rwwf%outwf                      '
 names(277)='rwwf%strnps                     '
 names(278)='rwwf%tddft                      '
 names(279)='rwwf%suscep                     '
 names(281)='rwwf%wfsinp(RF)                 '
 names(282)='rwwf%mkrho2                     '
 names(283)='rwwf%outwf2                     '
 names(284)='rwwf%dfpt_dyfnl                 '
 names(285)='rwwf%dfpt_mkrho                 '
 names(286)='rwwf%dfpt_nstwf                 '
 names(287)='rwwf%dfpt_vtorho                '
 names(288)='rwwf%dfpt_vtowfk                '
 names(289)='rwwf%dfpt_nstdy                 '
 names(290)='rwwf%initwf(RF)                 '
 names(291)='rwwf%newkpt(GS)                 '
 names(292)='rwwf%newkpt(RF)                 '

 names(301)='screening                       '
 names(302)='screening(init1)                '
 names(304)='screening(KS=>QP[wfrg])         '
 names(305)='screening(density)              '
 names(306)='screening(q-loop,init )         '
 names(307)='screening(cchi0q0)              '
 names(308)='screening(cchi0)                '
 names(309)='screening(q-loop,end)           '
 names(310)='screening(wrt scr files)        '
 names(315)='screening(pawin)                '
 names(316)='screening(wfs)                  '
 names(319)='screening(1)                    '
 names(320)='screening(paw)                  '; basic(320)=1
 names(321)='screening(2)                    '

 names(331)='cchi0                           '
 names(332)='cchi0(rho_tw_g)                 '
 names(333)='cchi0(assembly)                 '

 names(401)='sigma                           '
 names(402)='sigma(Init1)                    '
 names(403)='setup_sigma                     '
 names(404)='sigma(rdkss)                    '
 names(405)='sigma(Init2)                    '
 names(406)='sigma(make_vhxc)                '
 names(407)='sigma(vHxc_me)                  '
 names(408)='sigma(hqp_init)                 '
 names(409)='sigma(getW)                     '
 names(410)='sigma/=fourdp                   '; basic(410)=1

 names(421)='sigma(calc_sigx_me)             '
 names(423)='sigma(cohsex_me)                '
 names(424)='sigma(calc_sigc_me)             '
 names(425)='sigma(solve_dyson)              '
 names(426)='sigma(finalize)                 '

 names(430)='calc_sigx_me                    '

 names(431)='calc_sigc_me                    '
 names(432)='calc_sigc_me(Init)              '
 names(433)='calc_sigc_me(Init spin)         '
 names(434)='calc_sigc_me(Init q)            '
 names(435)='calc_sigc_me(eet_sigma)         '
 names(436)='calc_sigc_me(1)                 '
 names(437)='calc_sigc_me(rho_tw_g)          '
 names(438)='calc_sigc_me(2)                 '
 names(439)='calc_sigc_me(sigma_me)          '
 names(440)='calc_sigc_me(wfd_barrier        '
 names(441)='calc_sigc_me(xmpi_sum)          '
 names(442)='calc_sigc_me(final ops)         '

 names(445)='calc_sigc_me(loop)              '

 names(490)='solve_dyson                     '
 names(491)='cohsex_me                       '

 names(501)='nonlinear                       '
 names(502)='pead_nl_loop                    '
 names(503)='dfptnl_loop                     '
 names(511)='dfptnl_mv                       '; basic(511)=1
 names(512)='pead_nl_resp                    '; basic(512)=1
 names(513)='dfptnl_pert                     '
 names(514)='rf2_init                        '

 names(520)='lobpcgwf(init)                  '; if(abs(timopt)==4)basic(520)=1
 names(521)='lobpcgwf(bef.getghc 1           '; if(abs(timopt)==4)basic(521)=1
 names(522)='lobpcgwf(aft.getghc 1           '; if(abs(timopt)==4)basic(522)=1
 names(523)='lobpcgwf(bef.getghc 2           '; if(abs(timopt)==4)basic(523)=1
 names(524)='lobpcgwf(aft.getghc 2           '; if(abs(timopt)==4)basic(524)=1
 names(525)='lobpcgwf(aft.loop)              '; if(abs(timopt)==4)basic(525)=1
 names(526)='lobpcgwf(prep-getghc)           '

 names(530)='lobpcgwf                        '
 names(532)='xgemm%lobpcg                    '
 names(533)='xmpi_sum%lobpcg                 '
 names(535)='xorthon-xtrsm                   '
 names(536)='xprecon%lobpcg                  '
 names(537)='prep_fourwf%vtow                '
 names(538)='prep_fourwf%mkrh                '
 names(539)='prep_fourwf                     '

 names(540)='sg_fourwf%fourwf                '
 names(541)='back_wf%sg_fourw                '
 names(542)='forw_wf%sg_fourw                '
 names(543)='alltoall%back_wf                '
 names(544)='alltoall%forw_wf                '
 names(545)='prep_getghc(alltoall)           '
 names(547)='alltoall%prep_fo                '
 names(548)='allgather%prep_f                '
 names(549)='symrhg%mkrho                    '

 names(550)='forces:pawatm2ff                '
 names(551)='stress:pawatm2ff                '
 names(552)='setvtr:pawatm2ff                '
 names(553)='pawinit                         '; basic(553)=1
 names(554)='vtowfk:rhoij                    '
 names(555)='vtorho:pawmkrhoij               '; basic(555)=1
 names(556)='pawmkrho                        '; basic(556)=1
 names(557)='pawmkrho:symrhoij               '; basic(557)=1
 names(558)='scfcv_core:mknhat               '
 names(559)='nhatgrid                        '; basic(559)=1
 names(560)='pawdenpot                       '; basic(560)=1
 names(561)='pawdij/symdij                   '; basic(561)=1
 names(562)='respfn:pawatm2ff                '; basic(562)=1
 names(563)='dfpt_dyfro:pawatm2ff            '; basic(563)=1
 names(564)='dfpt_scfcv:dfpt_mknhat          '; basic(564)=1
 names(565)='getgsc                          '
 names(566)='dfpt_nstpaw                     '; basic(566)=1
 names(567)='pawnstd2e                       '
 names(568)='stress%strhar                   '

 names(570)='prep_nonlop                     '
 names(572)='prep_nonlop%vtowfk              '
 names(573)='prep_nonlop%forstrnps           '

 names(575)='prep_bandfft_tabs               '; basic(575)=1

 names(581)='prep_nonlop(alltoall)           '
 names(583)='vtowfk(pw_orthon)               '
 names(584)='xcopy%lobpcg                    '
 names(585)='vtowfk(subdiago)                '
 names(586)='vtowfk(nonlocalpart)            '
 names(587)='zheegv-dsyegv                   '

 names(588)='vtowfk(ssdiag)                  '; basic(588)=1
 names(589)='vtowfk(contrib)                 '; basic(589)=1
 names(590)='vtowfk(2)                       '
 names(591)='vtowfk(3)                       '

 names(593)='set_paw_pert                    '
 names(594)='get_exchange_atom               '
 names(595)='pawrhoij_redistribute           '
 names(596)='paw_ij_redistribute             '
 names(597)='paw_an_redistribute             '
 names(598)='pawfgrtab_redistribute          '

 names(600)='vtorhorec                       '
 names(601)='Definitions                     '
 names(602)='getngrec                        '
 names(603)='green_kernel                    '
 names(604)='transgrid (c->f)                '
 names(605)='recursion (other)               '
 names(606)='recursion (den)                 '
 names(607)='recursion (cuda)                '
 names(608)='recursion_nl                    '
 names(609)='fermisolverec                   '
 names(610)='entropyrec                      '
 names(611)='gran_potrec                     '
 names(612)='nonlocal-energy                 '
 names(613)='sync. cpu (wait)                '
 names(614)='sync. gpu (wait)                '
 names(615)='vn_nl_rec                       '
 names(616)='null recursion                  '
 names(617)='recursion (other_cuda)          '

 names(620)='datafordmft                     '
 names(621)='initialize dmft loop            '
 names(622)='impurity_solve                  '
 names(623)='Dyson                           '
 names(624)='compute_green                   '
 names(625)='integrate_green                 '
 names(626)='dmft-other                      '
 names(627)='Print/Read self                 '

 names(630)='prep_getghc                     '
 names(631)='prep_getghc(before if)          '
 names(632)='prep_getghc(bef. getghc)        '
 names(633)='prep_getghc(betw getghc)        '
 names(634)='prep_getghc(aft. getghc)        '
 names(635)='prep_getghc(getghc - 1 )        '
 names(636)='prep_getghc(getghc - 2 )        '
 names(637)='prep_getghc(getghc - 3 )        '
 names(638)='prep_getghc(getghc - 4 )        '

 names(640)='driver                          '
 names(641)='driver(bef. loop dtset)         '
 names(642)='driver(bef. select case)        '
 names(643)='driver(aft. select case)        '
 names(644)='driver(aft. loop dtset)         '

 names(650)='bse                             '
 names(651)='bse(Init1)                      '; basic(651)=1
 names(652)='setup_bse                       '; basic(652)=1
 names(653)='bse(rdkss)                      '; basic(653)=1
 names(654)='bse(rdmkeps^-1)                 '; basic(654)=1
 names(655)='bse(mkrho)                      '; basic(655)=1
 names(656)='bse(mkexcham)                   '; basic(656)=1
 names(657)='bse(mkexceps)                   '; basic(657)=1
 names(658)='bse(wfd_wave_free)              '; basic(658)=1
 names(659)='bse(mk_pawhur_t)                '; basic(659)=1
 names(660)='bse(exc_diago_driver)           '; basic(660)=1
 names(661)='bse(exc_haydock_driver)         '; basic(661)=1


 names(670)='exc_build_ham                   '
 names(671)='exc_build_ham(q=0)              '
 names(672)='exc_build_ham(block-res)        '
 names(673)='exc_build_ham(block-coupling)   '

 names(680)='exc_build_block                 '
 names(681)='exc_build_block(init,read)      '
 names(682)='exc_build_block(Coulomb)        '
 names(683)='exc_build_block(exchange)       '
 names(684)='exc_build_block(synchro)        '
 names(685)='exc_build_block(write_ha        '
 names(686)='exc_build_block(exch.spi        '

 names(690)='exc_haydock_driver              '
 names(691)='exc_haydock_driver(read)        '
 names(692)='exc_haydock_driver(prep)        '
 names(693)='exc_haydock_driver(wo lf        '
 names(694)='exc_haydock_driver(apply)       '
 names(695)='exc_haydock_driver(end)         '
 names(696)='exc_haydock_driver(inter        '
 names(697)='exc_haydock_driver(matmul)      '
!Slots up to 699 are reserved for bethe_salpeter code.

 names(700)='gstateimg                       '
 names(701)='gstate(pspini)                  '
 names(702)='gstateimg(leave_test)           '
 names(703)='gstateimg(init)                 '
 names(704)='gstateimg(bef. loop img)        '
 names(705)='gstateimg(bef. gstate)          '
 names(706)='gstateimg(aft. gstate)          '
 names(707)='gstateimg(aft. loop img)        '
 names(708)='gstateimg(finalize)             '


 names(710)='inwffil                         '
 names(711)='inwffil(read header)            '
 names(712)='inwffil(init params)            '
 names(713)='inwffil(prepa wfsinp)           '
 names(714)='inwffil(call wfsinp)            '
 names(715)='inwffil(after wfsinp)           '
 names(716)='inwffil(spin convert)           '
 names(717)='inwffil(call newkpt)            '
 names(718)='inwffil(excl. calls)            '; basic(718)=1

 names(720)='wfsinp                          '
 names(721)='wfsinp(before loop)             '
 names(722)='wfsinp(find kpt)                '
 names(723)='wfsinp(prepa initwf)            '
 names(724)='wfsinp(call  initwf)            '
 names(725)='wfsinp(transfer of wfs)         '
 names(726)='wfsinp(call rwwf)               '
 names(727)='wfsinp(wfconv section)          '
 names(728)='wfsinp(excl. calls)             '; basic(728)=1

 names(740)='suscep_stat                     '
 names(741)='suscep_stat(init)               '
 names(742)='suscep_stat(bef. susk-mm        '
 names(743)='suscep_stat(susk-mm)            '
 names(744)='suscep_stat(extrapol)           '
 names(745)='suscep_stat:synchro             '
 names(746)='suscep_stat:MPI                 '
 names(747)='suscep_stat(symmetries)         '

 names(750)='susk                            '
 names(751)='susk (init)                     '; basic(751)=1
 names(752)='susk (loop)                     '
 names(753)='susk:MPI (1)                    '; basic(753)=1
 names(754)='susk (accumul.)                 '
 names(755)='susk:MPI (2)                    '; basic(755)=1
 names(756)='susk (loop except FFT)          '; basic(756)=1
 names(757)='susk (accumul.except FFT        '; basic(757)=1

 names(760)='suskmm                          '
 names(761)='suskmm (init)                   '; basic(761)=1
 names(762)='suskmm (loop : part1)           '
 names(763)='suskmm (loop : part2)           '
 names(764)='suskmm(loop1 except FFT)        '; basic(764)=1
 names(765)='suskmm(loop2 except FFT)        '; basic(765)=1

 names(770)='initwf                          '
 names(771)='initwf(before rwwf)             '; basic(771)=1
 names(772)='initwf(after rwwf)              '; basic(772)=1

 names(780)='newkpt                          '
 names(781)='newkpt(before loop)             '
 names(782)='newkpt(before rwwf)             '
 names(783)='newkpt(after rwwf)              '
 names(784)='newkpt(call wfconv)             '
 names(785)='newkpt(finalize loop)           '
 names(786)='newkpt(after loop   )           '
 names(787)='newkpt:synchro                  '
 names(788)='newkpt(excl. rwwf   )           '; basic(788)=1

 names(790)='mkrho                           '
 names(791)='mkrho%gstate                    '
 names(792)='mkrho%vtorho                    '
 names(793)='mkrho%energy                    '
 names(794)='mkrho%respfn                    '
 names(795)='mkrho%afterscfloop              '
 names(796)='mkrho%scfcv_core                '
 names(798)='mkrho/=                         '; basic(798)=1
 names(799)='mkrho/=+fourwf                  '

 names(801)='fourwf                          '
 names(802)='fourwf%(pot)                    '; basic(802)=1;    ndata(802)=2*nfft
 names(803)='fourwf%(den)                    '; basic(803)=1;    ndata(803)=nfft
 names(804)='fourwf%(G->r)                   '; basic(804)=1
 names(805)='fourwf%(r->G)                   '; basic(805)=1


 names(840)='fourwf%(other)                  '
 names(841)='fourwf%getghc                   '
 names(842)='fourwf%vtowfk                   '
 names(843)='fourwf%mkrho                    '
 names(844)='fourwf%dfpt_cgwf                '
 names(845)='fourwf%dfpt_accrho%dfpt_vtowfk  '
 names(846)='fourwf%mkrho2                   '
 names(847)='fourwf%dfpt_mkrho               '
 names(854)='fourwf%tddft                    '
 names(855)='fourwf%outkss                   '
 names(856)='fourwf%prep_four                '
 names(858)='fourwf%dfpt_accrho%idfpt_nstpaw '
 names(861)='fourwf%suskmm !0 part 1         '
 names(862)='fourwf%suskmm !0 part 2         '
 names(871)='fourwf%suskmm !3 part 1         '
 names(872)='fourwf%suskmm !3 part 2         '

 names(901)='newvtr(before selection)        '
 names(902)='newvtr(bef. prcref_PMA)         '
 names(903)='newvtr(call prcref_PMA)         '
 names(904)='newvtr(aft. prcref_PMA)         '
 names(905)='newvtr(mean potential)          '

 names(910)='forstr                          '
 names(911)='forstr(forstrnps)               '
 names(912)='forstr(pawgrnl)                 '
 names(913)='forstr(forces)                  '
 names(914)='forstr(stress)                  '

 names(920)='forstrnps                       '
 names(921)='forstrnps(bef.loop spin)        '
 names(922)='forstrnps(bef.loop band)        '
 names(923)='forstrnps(copy)                 '
 names(924)='forstrnps(kinetic contr)        '
 names(925)='forstrnps(aft.loop kptsp        '
 names(926)='forstrnps(nonlop+prep_ba        '
 names(927)='forstrnps(bef.loop kpt)         '

 names(933)='outkss                          '
 names(934)='outkss(Gsort+hd)                '
 names(935)='outkss(k-loop)                  '
 names(936)='outkss(diago)                   '; basic(936)=1
 names(937)='outkss(MPI_exch)                '; basic(937)=1
 names(938)='outkss(write)                   '

 names(940)='rhotov                          '
 names(941)='rhotov(rhotoxc)                 '
 names(942)='rhotov(dotprod_vn)              '
 names(943)='rhotov(PSolver_rhohxc)          '
 names(944)='rhotov(rhohxcpositron)          '
 names(945)='rhotov(other)                   '

 names(950)='outscfcv                        '
 names(951)='outscfcv(mlwfovlp)              '
 names(952)='outscfcv([PAW]prtden)           '
 names(953)='outscfcv(prtelf)                '
 names(954)='outscfcv(prt[g,k,l]den)         '
 names(955)='outscfcv(prtwf)                 '
 names(956)='outscfcv(prtpot)                '
 names(957)='outscfcv(prt geo misc.)         '
 names(958)='outscfcv(prt stm,vha,..)        '
 names(959)='outscfcv(prtdos)                '
 names(960)='outscfcv(calcdenmagsph)         '
 names(961)='outscfcv(pawprt)                '
 names(962)='outscfcv(optics)                '
 names(963)='outscfcv(pawmkaewf)             '
 names(964)='outscfcv(outkss)                '
 names(965)='outscfcv(poslifetime)           '
 names(966)='outscfcv(outwant)               '
 names(967)='outscfcv(cal[cs,efg,fc])        '
 names(968)='outscfcv(prt[surf,nest])        '
 names(969)='outscfcv(misc.)                 '

 names(980)='vtorho                          '
 names(981)='vtorho(bef. spin loop)          '
 names(982)='vtorho(bef. kpt  loop)          '
 names(983)='vtorho(Berry)                   '
 names(984)='vtorho(bef. vtowfk)             '
 names(985)='vtorho(aft. vtowfk)             '
 names(986)='vtorho(aft. kpt loop)           '
 names(987)='vtorho(leave_test)              '; basic(987)=1
 names(988)='vtorho(aft. spin loop)          '
 names(989)='vtorho(MPI)                     '; basic(989)=1
 names(990)='vtorho(newocc)                  '
 names(991)='vtorho(DMFT)                    '
 names(992)='vtorho(mkrho 1)                 '
 names(993)='vtorho(highest occ. eig)        '
 names(994)='vtorho(mkrho 2)                 '
 names(995)='vtorho(tddft)                   '
 names(996)='vtorho(suscep_stat)             '
 names(997)='vtorho(init kpt loop)           '

 names(1001)='initberry                       '; basic(1001)=1
 names(1002)='initberry(before listkk)        '
 names(1003)='initberry(call listkk)          '
 names(1004)='initberry(after listkk)         '
 names(1005)='initberry(find neighb.)         '
 names(1006)='initberry(build strings)        '
 names(1007)='initberry(PAW on-site)          '
 names(1008)='initberry(pwind)                '
 names(1009)='initberry(MPI stuff)            '

 names(1021)='listkk                          '; basic(1021) = 1

! CMartins: TEST for HF
 names(1501)='HF_init                         '; basic(1501)=1
 names(1502)='HF_updatecgocc                  '; basic(1502)=1
 !===START===containers used for testing purposes
 names(1503)='HF_updatecgocc-MPI              ';! basic(1503)=1  ! 100 % nested inside 1502
 names(1504)='HF_getghc                       ';! basic(1504)=1  !1504 = 1505 + 1506 + 1511 + extra_loops
 names(1505)='HF_getghc-init                  ';! basic(1505)=1  ! 100 % nested inside 1504
 names(1506)='HF_getghc-kmu_loop              ';! basic(1506)=1  ! 1506 = 1508 + 1509 + 1510 + 1507 + extra_loops
 names(1507)='HF_getghc-calc_vlocpsi          ';! basic(1507)=1  ! 100 % nested inside 1506
 names(1508)='HF_getghc-mult-cwf*cwocc        ';! basic(1508)=1  ! 100 % nested inside 1506
 names(1509)='HF_getghc-calc_rhog_munu        ';! basic(1509)=1  ! 100 % nested inside 1506
 names(1510)='HF_getghc-calc_vloc             ';! basic(1510)=1  ! 100 % nested inside 1506
 names(1511)='HF_getghc-calc_ghc              ';! basic(1511)=1  ! 100 % nested inside 1504
 !===STOP===containers used for testing purposes
 names(1515)='HF_getghc_main                  ';  basic(1515)=1  ! ulterior slot for test

 ! Chebfi
 names(1600) = 'chebfi                        '; basic(1601) = 1
 names(1601) = 'chebfi(alltoall)              '
 names(1602) = 'chebfi(appinvovl)             '
 names(1603) = 'chebfi(rotation)              '
 names(1604) = 'chebfi(subdiago)              '
 names(1605) = 'chebfi(subham)                '
 names(1606) = 'chebfi(ortho)                 '
 names(1607) = 'chebfi(getghc)                '
 names(1608) = 'chebfi(residuals)             '
 names(1609) = 'chebfi(update_eigens)         '
 names(1610) = 'chebfi(sync)'

 names(1630) = 'chebfi(opernla)               '
 names(1631) = 'chebfi(opernlb)               '
 names(1632) = 'chebfi(inv_s)                 '

 names(1620) = 'mkinvovl                      '
 names(1621) = 'mkinvovl(build_d)             '
 names(1622) = 'mkinvovl(build_ptp)           '

 ! lobpcg2
 names(1650) = 'lobpcgwf2                     '; basic(1650) = 1
 names(1651) = 'lobpcg_init                    '
 names(1652) = 'lobpcg_free                    '
 names(1653) = 'lobpcg_run                     '
 names(1654) = 'lobpcg_getAX_BX                '
 names(1655) = 'lobpcg_orthoWrtPrev            '
 names(1656) = 'lobpcg_Bortho                  '
 names(1657) = 'lobpcg_RayleighRitz            '
 names(1658) = 'lobpcg_maxResidu               '
 names(1659) = 'lobpcg_run@getAX_BX            '
 names(1660) = 'lobpcg_pcond                   '
 names(1661) = 'lobpcg_RayleighRitz@hegv       '

 ! xg_t
 names(1662) = 'xgTransposer_transpose@ColsRows'
 names(1663) = 'xgTransposer_transpose@Linalg  '
 names(1664) = 'xgTransposer_*@all2all         '
 names(1665) = 'xgTransposer_*@gatherv         '
 names(1666) = 'xgTransposer_@reorganize       '
 names(1667) = 'xgTransposer_init              '
 names(1668) = 'xgTransposer_free              '
 names(1669) = 'xgTransposer_transpose         '
 names(1670) = 'xgBlock_potrf                  '
 names(1671) = 'xgBlock_trsm                   '
 names(1672) = 'xgBlock_gemm                   '
 names(1673) = 'xgBlock_set                    '
 names(1674) = 'xgBlock_get                    '
 names(1675) = 'xgBlock_heev                   '
 names(1676) = 'xgBlock_heevd                  '
 names(1677) = 'xgBlock_hpev                   '
 names(1678) = 'xgBlock_hpevd                  '
 names(1679) = 'xgBlock_hegv                   '
 names(1680) = 'xgBlock_hegvx                  '
 names(1681) = 'xgBlock_hegvd                  '
 names(1682) = 'xgBlock_hpgv                   '
 names(1683) = 'xgBlock_hpgvx                  '
 names(1684) = 'xgBlock_hpgvd                  '
 names(1685) = 'xgBlock_copy                   '
 names(1686) = 'xgBlock_cshift                 '
 names(1687) = 'xgBlock_pack                   '
 names(1690) = 'xgScalapack_init               '
 names(1691) = 'xgScalapack_free               '
 names(1692) = 'xgScalapack_heev               '
 names(1693) = 'xgScalapack_hegv               '
 names(1694) = 'xgScalapack_scatter            '

 ! GWLS GW code
 names(1701)='gwls_sternheimer                ';basic(1701)=1
 names(1702)='exchange and correlation        '
 names(1703)='correl. shift lanczos           '
 names(1704)='Dielectric matrix               '
 names(1705)='Model Dielectric matrix         '
 names(1706)='setup proj. sternheimer         '
 names(1707)='compute proj.sternheimer        '
 names(1708)='eps^{-1} - eps_m^{-1}           '
 names(1709)='eps_m^{-1} - 1                  '
 names(1710)='Modify Lbasis Coulomb           '
 names(1711)='Diag eps^{-1}-eps_m^{-1}        '
 names(1712)='exact  AT shift lanczos         '
 names(1713)='model  AT shift lanczos         '
 names(1714)='exact  BT shift lanczos         '
 names(1715)='model  BT shift lanczos         '
 names(1716)='compute poles                   '
 names(1717)='Sigma_A Lanczos                 '
 names(1718)='Sigma_B num. integrands         '


 names(1719)='gwls: extract_QR                ';basic(1719)=1
 names(1720)='gwls: extract_SVD               ';basic(1720)=1

 ! these entry are not in a logical order.
 names(1721)='gwls: gstateimg                 '
 names(1722)='prepareValenceWfk               '

 names(1723)='gwls: sqmr                      ';basic(1723)=1


 names(1724)='gwls: Pk                        ';basic(1724)=1
 names(1725)='Pk- allocating                  '
 names(1726)='Pk- wfk to denpot               '
 names(1727)='Pk- wfk product with val        '
 names(1728)='Pk- pc_k                        '
 names(1729)='Pk- sqmr case 1                 '
 names(1730)='Pk- sqmr case 2                 '
 names(1731)='Pk- sqmr case 3                 '
 names(1732)='Pk-  qmr case 4                 '
 names(1733)='Pk- apply H (case 2)            '


 names(1734)='gwls: Pk_model                  ';basic(1734)=1
 names(1735)='Pk_model- allocating            '
 names(1736)='Pk_model- wfk to denpot         '
 names(1737)='Pk_model- wfk x val             '
 names(1738)='Pk_model- pc_k                  '
 names(1739)='Pk_model- act with Y            '
 names(1740)='Pk_model- add contrib.          '


 names(1741)='gwls: calc eps_m^-1(w)-1        ';basic(1741)=1
 names(1742)='Allocating                      '
 names(1743)='modifying Lanczos basis         '
 names(1744)='calc <mod_L_1|Y|mod_L_2>        '
 names(1745)='    make array hermitian        '
 names(1746)='               xsum_mpi         '
 names(1747)='inv eps_m and subtract 1        '

 ! IFC object
 names(1748)='ifc_fourq'; basic(1748) = 1
 !names(1749)='ewald9'; basic(1749) = 1
 !names(1750)='gtdyn9'; basic(1750) = 1
 !names(1751)='dfpt_phfrq'; basic(1751) = 1

 names(1780)='cgtk_rotate'; basic(1780) = 1

 ! DVDB object
 names(1800)='dvdb_new'; basic(1800) = 1
 names(1801)='dvdb_qcache_read'; basic(1801) = 1
 names(1802)='dvdb_readsym_qbz'; basic(1802) = 1
 names(1803)='dvdb_rotate_fqg'; basic(1803) = 1
 names(1804)='v1phq_rotate'; basic(1804) = 1
 names(1805)='dvdb_readsym_allv1'; basic(1805) = 1
 names(1806)='dvdb_xmpi_sum'; basic(1806) = 1
 names(1807)='dvdb_qcache_update'; basic(1807) = 1
 names(1808)='dvdb_ftqcache_build'; basic(1808) = 1
 names(1809)='dvdb_get_ftqbz'; basic(1809) = 1

 ! SIGEPH
 !names(1900)='sigph_pre_qloop'; basic(1900) = 1
 !names(1901)='sigph_qloop_preamble'; basic(1901) = 1
 !names(1902)='sigph_qloop_cg_and_h1'; basic(1902) = 1
 names(1903)='sigph_bsum'; basic(1903) = 1

 !names(1904)='rf_transgrid_and_pack'; basic(1904) = 1
 !names(1905)='splfit'; basic(1905) = 1
 !names(1906)='mkkin'; basic(1906) = 1

 names(TIMER_SIZE)='(other)                         ' ! This is a generic slot, to compute a complement

!==================================================================================

 spaceworld= mpi_enreg%comm_world
 nproc     = mpi_enreg%nproc
 me        = mpi_enreg%me
 nthreads  = xomp_get_num_threads(open_parallel=.true.)

 call timab(49,2,tsec)

 if(abs(timopt)==1 .or. timopt==-3 .or. timopt==-4)then ! Time the timing routine (precision should be better than 3%)
   ltimab=1
   utimab=1000
   maxii=20
!  maxii=1    ! Uncomment this line if no timer is provided in timein.f
   do ii=1,20

     call timab(50,1,tsec)
     do itimab=ltimab,utimab
!      The channel 51 is here used as a dummy channel
       call timab(51,1,tsec)
       call timab(51,2,tsec)
     end do
     call timab(50,2,tsec)
     call time_accu(50,return_ncount,tsec, lflops, ftsec)
!    Exit the timing loop if the CPU time is bigger than 0.10 second
!    of if the number of calls is too large.
!    Since the accuracy of the timing is expected to be better than 0.01 sec,
!    gives about 10% accuracy
     if(tsec(1)>0.10_dp)then
       exit
     else
       ltimab=utimab+1
!      Increase the number of timab calls in a block.
!      This small factor of increase allows to have less than
!      0.15 second for this testing
       utimab=(3*utimab)/2
     end if
   end do
!  Get the time per combined call timab(*,1,tsec) + timab(*,2,tsec)
   timab_cpu=tsec(1)/utimab
   timab_wall=tsec(2)/utimab
   if(timopt<0 .and. me==0 .and. timopt/=-2)then
     write(ount,*)
     write(ount,*)'Test the timer : '
     write(ount,*)' a combined call timab(*,1,tsec) + timab(*,2,tsec) is '
     write(ount, '(a,es14.4,a,es14.4,a)' )'- CPU time =',timab_cpu,' sec,    Wall time =',timab_wall,' sec'
   end if
 else
   timab_cpu=zero; timab_wall=zero
 end if

!Eventually reenable the timab routine
 call timab(1,5,tsec)

!Get overall elapsed cpu and wall clock time
 call timab(1,2,tsec)
 call time_accu(1,return_ncount,tsec,lflops,ftsec)
 ncount(1)=return_ncount

!Sum over all procs
 my_tsec(:)=tsec(:)
 call xmpi_sum(my_tsec,tsec,2,spaceworld,ierr)

!Only the world master writes
 if (me==0) then
   write(ount,'(/,a,f13.1,f12.2,f11.3)')'- Total cpu        time (s,m,h):',tsec(1),tsec(1)/60._dp,tsec(1)/3600._dp
   write(ount,'(a,f13.1,f12.2,f11.3)')  '- Total wall clock time (s,m,h):',tsec(2),tsec(2)/60._dp,tsec(2)/3600._dp
 end if

!Get separate time reports from all timed sections
 totcount=0
 do itim=1,TIMER_SIZE
   call time_accu(itim,return_ncount,times(:,itim),nflops(itim),ftimes(:,itim))
   ncount(itim)=return_ncount
   totcount=totcount+return_ncount
 end do

!Estimate additional timings.

!Estimate the values associated with timab, put it in channel 51
 ncount(51)=totcount
 times(1,51)=timab_cpu*totcount
 times(2,51)=timab_wall*totcount

!Gather the different parts of selected time slots
!Or, alternatively, deduce the value of the complement of some time slots.
!This loop is finished when the default case is hit (see below)
 do ii=1,TIMER_SIZE

   tslots(:)=0

!  List first the time slot in which the result will be accumulated.
!  If this number is negative, the positive value will be used for the time slot, but the ncount will be set to -1 .
!  Then, list the time slots whose value will be either accumulate or subtracted. The latter is obtained by
!  entering a minus sign in front of the time slot number ...
!  If a negative number is present in the list, while the accumulated time slot is positive,
!  then the number of counts will be set to the value of the first routine to be accumulated.
   select case(ii)
!    Gather the different parts of nonlop  (SHOULD BE REEXAMINED !)
   case(1)
     tslots(:5)=(/75, 221,223,229,233/)
   case(2)
     tslots(:4)=(/76, 222,225,227/)
   case(3)
     tslots(:2)=(/77, 224/)
   case(4)
     tslots(:2)=(/78, 226/)
   case(5)
     tslots(:2)=(/79, 228/)
   case(6)
!      Gather the different parts of selected time channels
     tslots(:10)=(/97, 75,76,77,78,79,220,230,231,232/)
   case(7)
!      Gather the different parts of fourwf (NOTE : should attribute the channel 840 to one of the 4 modes !!!)
     tslots(:3)=(/802, 841,844/)
   case(8)
     tslots(:4)=(/803, 842,843,846/)
   case(9)
     tslots(:10)=(/804, 845,847,848,850,854,858,859,861,862/)
   case(10)
     tslots(:6)=(/805, 849,851,857,871,872/)
   case(11)
!      In the following, the part coming from the prep_fourwf interface is added to the total.
     tslots(:7)=(/801, 802,803,804,805,840,856/)
   case(13)
!      Gather the different parts of prep_fourwf
     tslots(:3)=(/539, 537,538/)
   case(14)
!      Gather the different parts of fourdp
     tslots(:11)=(/9, 260,261,262,263,264,265,266,267,268,269/)
   case(15)
!      Gather the different parts of getghc
     tslots(:9)=(/200, 201,202,203,204,205,206,207,208/)
   case(16)
!      Gather the different parts of projbd
     tslots(:3)=(/210, 211,212/)
   case(17)
!      Gather the different parts of rwwf (wavefunctions read/write)
     tslots(:24)=&
&     (/14, 270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292/)
   case(18)
!      Estimate the complement of getghc (non fourwf, non nonlop)
     tslots(:4)=(/-98, 200,-841,-221/)
   case(19)
!      Estimate the complement of cgwf (non getghc,projbd)
     tslots(:5)=(/-40, 22,530,-201,-211/)
   case(20)
!      Estimate the complement of dfpt_cgwf (non getghc,projbd,nonlop,fourwf)
     tslots(:8)=(/-140, 122,-202,-197,-212,-227,-228,-844/)
   case(21)
!      Estimate different complements in vtowfk
!      vtowfk(ssdiag) (= vtowfk(loop)    - cgwf )
     tslots(:5)=(/-588, 39,-22,-530, -1600/)
   case(22)
!      vtowfk(contrib) (= vtowfk (afterloop) - nonlop%vtowfk - fourwf%vtowfk )
     tslots(:4)=(/589, 30,-222,-842/)
   case(23)
!      vtowfk (1) = vtowfk - vtowfk(loop) - vtowfk(afterloop)
     tslots(:4)=(/31, 28,-39,-30/)
   case(24)
!      Estimate different complements in dfpt_vtowfk
!      dfpt_vtowfk(contrib) (= vtowfk3(loop) - cgwf - fourwf%vtowfk3 - rwwf%vtowfk3 - corrmetalwf1)
     tslots(:6)=(/-108, 139,-122,-845,-288,-214/)
   case(25)
!      vtowfk (1) = dfpt_vtowfk - vtowfk3(loop) - vtowfk3 (3)
     tslots(:4)=(/ 131, 128,-139,-130/)
   case(28)
!      dfpt_vtorho-kpt loop (= dfpt_vtowfk (2) - vtowfk3 - rwwf)
     tslots(:4)=(/126,125,-128,-287/)
   case(29)
!      Estimate complement in mkrho
     tslots(:3)=(/798,799,-843/)
   case(30)
!      Estimate complement in dfpt_looppert
!      dfpt_looppert(other) (= loper3 - loper3(kpgio) - loper3(getmpw) - loper3(inwffil)
!      dfpt_scfcv - dfpt_looppert(outwf) -loper3(eigt2tot)
     tslots(:8)=(/145,141,-142,-143,-144,-120,-146,-147/)
   case(31)
!      Estimate complement in sigma
!      sigma/=fourdp = sigma - fourdp%rhotwg%si - fourdp%fftwfn%si
     tslots(:4)=(/410,401,-262,-265/)
   case(32)
!      Estimate complement in bethe_salpeter
     tslots(:2)=(/699,650/)
   case(33)
!      Estimate complement in susk
!      NOTE : fourwf%susk _PAW should actually be split between susk (loop except FFT)
!      and susk (accumul.except FFT . But a renumbering of the fourwf splitting should be done ...
!      susk (loop except FFT) = susk (loop) - fourwf%susk !0 - fourwf%susk !3
     tslots(:4)=(/756,752,-848,-849/)
   case(34)
!      susk (accumul.except FFT = susk (accumul) - fourwf%susk !3bis - fourwf%susk _PAW
     tslots(:4)=(/757,754,-859,-857/)
   case(35)
!      Estimate complement in suskmm
!      NOTE : fourwf%susk _PAW should actually be split between susk (loop except FFT)
!      and suskmm (accum.except FFT . But a renumbering of the fourwf splitting should be done ...
!      suskmm (loop except FFT) = suskmm (loop) - fourwf%suskmm !0 part 1 - fourwf%suskmm !3 part 1
     tslots(:4)=(/764,762,-861,-871/)
   case(36)
!      suskmm (accum.except FFT = suskmm (accumul) - fourwf%suskmm !0 part 2 - fourwf%suskmm !3 part 2 - fourwf%susk _PAW
     tslots(:5)=(/765,763,-862,-872,-857/)
   case(37)
!      inwffil(excl. calls) = inwffil - inwffil(call wfsinp) - inwffil(call newkpt);
     tslots(:4)=(/718,710,-714,-717/)
   case(38)
!      wfsinp(excl. calls) = wfsinp - wfsinp(call  initwf) - wfsinp(call rwwf)
     tslots(:4)=(/728,720,-724,-727/)
   case(39)
!      newkpt(excl. rwwf   )=newkpt(before loop) + newkpt(before rwwf) + newkpt(after rwwf)
!      newkpt(call wfconv) + newkpt(finalize loop) + newkpt(after loop   )
     tslots(:7)=(/-788,781,782,783,784,785,786/)
   case(40)
!      More complements in vtowfk
!      vtowfk (2) = vtowfk (loop) - cgwf - lobpcg - subdiago - pw_orthon
     tslots(:7)=(/-590,39,-22,-530,-585,-583, -1600/)
   case(41)
!      vtowfk (3) = vtowfk (afterloop) - nonlop%vtowfk - prep_nonlop%vtowfk - fourwf%vtowfk - prep_fourwf%vtowfk - vtowfk(nonlocalpart)
     tslots(:7)=(/-591,30,-222,-572,-842,-537,-586/)
   case(43)
!      mkrho = mkrho%gstate + mkrho%vtorho + mkrho%energy + mkrho%respfn + mkrho%afterscfloop + mkrho%scfcv_core
     tslots(:7)=(/790,791,792,793,794,795,796/)
   case(44)
!      Estimate the complement of dmft (in vtorho, only)
     tslots(:9)=(/-626, 991,-620,-621,-622,-623,-624,-625,-627/)

   case default
     cycle
   end select

   tslot=tslots(1)
   aslot=abs(tslot)
   ncount(    aslot)=0 ; if (tslot<0)ncount(aslot)=-1
   times(1:2, aslot)=zero
   nflops(    aslot)=zero
   ftimes(1:2,aslot)=zero
   flag_count=1
   do islot=2,TIMER_SIZE
     bslot=tslots(islot)
     cslot=abs(bslot)
     if(bslot>0)then
       if(tslot>0)ncount(aslot)=ncount(aslot)+ncount(cslot)
       times(1:2, aslot)=times(1:2, aslot)+times(1:2,cslot)
       nflops(    aslot)=nflops(    aslot)+nflops(   cslot)
       ftimes(1:2,aslot)=ftimes(1:2,aslot)+ftimes(1:2,cslot)
     else if(bslot<0)then
       if(tslot>0)flag_count=-1
       times(1:2, aslot)=times(1:2, aslot)-times(1:2,cslot)
       nflops(    aslot)=nflops(    aslot)-nflops(   cslot)
       ftimes(1:2,aslot)=ftimes(1:2,aslot)-ftimes(1:2,cslot)
     else if(bslot==0)then
       exit
     end if
   end do
   if(flag_count==-1)ncount(aslot)=ncount(abs(tslots(2)))
 end do

!For the following sections, the number of counts is non standard, and thus these sections have not been placed
!in the previous doloop.

!Compute xc part of rhotoxc and dfpt_mkvxc, minus the calls to fourdp inside that part
 ncount(11)=ncount(81)+ncount(181)
 times(1:2,11)=times(1:2,81)+times(1:2,181)-times(1:2,82)
 ftimes(1:2,11)=ftimes(1:2,81)+ftimes(1:2,181)-ftimes(1:2,82)
 nflops(11)=nflops(81)+nflops(181)-nflops(82)

!Estimate different complements in dfpt_vtorho
!dfpt_vtorho (1) (= vtorho3 (1,2) - vtorho3(2) - vtorho3:synchro )
 ncount(118)=ncount(121)
 times(1:2,118)=times(1:2,124)-times(1:2,125)-times(1:2,166)
 ftimes(1:2,118)=ftimes(1:2,124)-ftimes(1:2,125)-ftimes(1:2,166)
 nflops(118)=nflops(124)-nflops(125)-nflops(166)

!Calculating Gigaflops for all cases
 do itim=1,TIMER_SIZE
   mflops(itim)=-2
   if(abs(ftimes(1,itim)) > tol10) then ! VALGRIND complains that here there is a jump on uninitialized values
     mflops(itim)=nflops(itim)*1.e-9/ftimes(1,itim)
   else
     mflops(itim)=-1
   end if
 end do

!Warning if the time is negative
 do itim=1,TIMER_SIZE
   if(times(1,itim)<-tol6 .or. times(2,itim)<-tol6 .or. ncount(itim)<-1 )then
     write(msg, '(6a,i4,4a,es16.6,a,es16.6,a,i6,a,es16.6)' ) ch10,&
      ' timana: WARNING -',ch10,&
      '  One among cpu, wall and ncount is negative.',ch10,&
      '  Timing section #',itim,', name :  ',names(itim),ch10,&
      '  CPU =',times(1,itim),', Wall=',times(2,itim),' ncount=',ncount(itim),' flops=',nflops(itim)
     call wrtout(std_out,msg,'PERS')
   end if
 end do

!List of major independent code sections
 ABI_MALLOC(list, (TIMER_SIZE))
 list(:)=0
 nlist=0
 do itim=1,TIMER_SIZE
   if(basic(itim)/=0)then
     nlist=nlist+1
     list(nlist)=itim
   end if
 end do

 percent_limit=0.5_dp; if (timopt<0) percent_limit=0.0001_dp
 !percent_limit=tol12

!In case there is parallelism, report times for node 0
!if (me==0 .and. nproc>1) then
 if (me==0) then

!  Find normalization to report timing as % total time
   cpunm=100._dp/tsec(1)
   wallnm=100._dp/tsec(2)

!  (0) Take care of major independent code sections for this account of node 0 timing

   write(ount,  '(a,a,a,a,/,a,a,a)' ) '-',ch10,&
    '- For major independent code sections,',' cpu and wall times (sec),',&
    '-  as well as % of the time and number of calls for node 0',&
    '-'

   write(ount,"(3(a,i0),a)")&
    "-<BEGIN_TIMER mpi_nprocs = ",nproc,", omp_nthreads = ",nthreads,", mpi_rank = ",me,">"

!  write(ount,"(2(a,f13.1))")"- tot_cpu_time = ",tsec(1),   ", tot_wall_time = ",tsec(2)
   write(ount,"(2(a,f13.1))")"- cpu_time =  ",my_tsec(1),", wall_time =  ",my_tsec(2)
   write(ount,"(a)")"-"

   write(ount, '(a,t34,a,t42,a,t50,a,t59,a,t65,a,t82,a,3x,a7,1x,a10)' )&
     '- routine','cpu','%','wall','%',' number of calls ',' Gflops ', 'Speedup', 'Efficacity'
   write(ount,'(a,t35,a,t43,a,t51,a,t60,a,t66,a,t78,a)')&
     '-                ','   ',' ','    ',' ','  (-1=no count)'

!  Sort the list by decreasing CPU time
   do ii=1,nlist
     do ilist=1,nlist-1
       if (times(1,list(ilist))<times(1,list(ilist+1))) then
         temp_list=list(ilist)
         list(ilist)=list(ilist+1)
         list(ilist+1)=temp_list
       end if
     end do
   end do

   subcpu=zero; subwal=zero; other_cpu=zero; other_wal=zero; nothers=0

   do ilist=1,nlist
     isort = list(ilist)

     if ( (times(1,isort)*cpunm  > percent_limit .and. &
           times(2,isort)*wallnm > percent_limit) .and. ncount(isort)/=0 ) then ! Timing analysis

       write(ount,format01041)names(isort),&
         times(1,isort),times(1,isort)*cpunm,times(2,isort),times(2,isort)*wallnm,ncount(isort),mflops(isort), &
         times(1,isort)/times(2,isort),times(1,isort)/times(2,isort)/nthreads

     else
       nothers=nothers+1
       other_cpu=other_cpu+times(1,isort)
       other_wal=other_wal+times(2,isort)
     end if

     subcpu=subcpu+times(1,isort)
     subwal=subwal+times(2,isort)
   end do

   other_wal = other_wal + tol14
   write(entry_name,"(a,i0,a)")"others (",nothers,")"
   write(ount,format01041)entry_name,other_cpu,other_cpu*cpunm,other_wal,other_wal*wallnm,-1,-1.0, &
     other_cpu/other_wal,other_cpu/other_wal/nthreads
   write(ount,"(a)")"-<END_TIMER>"

   write(ount,'(a)' ) '-'
   subwal = subwal + tol14
   write(ount,01200) subcpu,subcpu*cpunm,subwal,subwal*wallnm,subcpu/subwal,subcpu/subwal/nthreads
 end if

!Now, gather all information
 call xmpi_sum(times,spaceworld,ierr)
 call xmpi_sum(ncount,spaceworld,ierr)
 call xmpi_sum(ftimes,spaceworld,ierr)
 call xmpi_sum(nflops,spaceworld,ierr)

 if (me==0) then ! Only the world master writes

!  Find normalization to report timing as % total time
   cpunm=100._dp/tsec(1)
   wallnm=100._dp/tsec(2)

!  Calculating Gigaflops for all process
   do itim=1,TIMER_SIZE
     mflops(itim)=-2
     if(abs(ftimes(1,itim)) > tol10) then ! VALGRIND complains that here there is a jump on uninitialized values
       mflops(itim)=nflops(itim)*1.e-9/ftimes(1,itim)
     else
       mflops(itim)=-1
     end if
   end do

!  _______________________________________

!  Write timing output for cpu times

!  (1) Take care of major independent code sections
   write(ount,'(/,a,/,a,/)' )&
     '- For major independent code sections, cpu and wall times (sec),',&
     '- as well as % of the total time and number of calls '

   write(ount,"(2(a,i0),a)")&
     "-<BEGIN_TIMER mpi_nprocs = ",nproc,", omp_nthreads = ",nthreads,", mpi_rank = world>"

   write(ount,"(2(a,f13.1))")"- cpu_time = ",tsec(1),   ", wall_time = ",tsec(2)
!  write(ount,"(2(a,f13.1))")"- my_cpu_time =  ",my_tsec(1),", my_wall_time =  ",my_tsec(2)
   write(ount,"(a)")"-"

   write(ount,'(a,t35,a,t43,a,t51,a,t60,a,t66,a,t82,a,3x,a7,1x,a10)')&
    '- routine        ','cpu','%','wall','%', ' number of calls ',' Gflops ', &
    'Speedup', 'Efficacity'
   write(ount,'(a,t35,a,t43,a,t51,a,t60,a,t66,a,t78,a)')&
    '-                ','   ',' ','    ',' ','  (-1=no count)'

!  Sort the list by decreasing CPU time
   do ii=1,nlist
     do ilist=1,nlist-1
       if(times(1,list(ilist))<times(1,list(ilist+1)))then
         temp_list=list(ilist)
         list(ilist)=list(ilist+1)
         list(ilist+1)=temp_list
       end if
     end do
   end do

   subcpu=zero; subwal=zero; other_cpu=zero; other_wal=zero; nothers=0

   do ilist=1,nlist
     isort = list(ilist)
     if( (times(1,isort)*cpunm > percent_limit .and. times(2,isort)*wallnm> percent_limit) .and. ncount(isort)/=0 )then

       write(ount,format01041)names(isort),&
         times(1,isort),times(1,isort)*cpunm,times(2,isort),times(2,isort)*wallnm,ncount(isort),mflops(isort), &
         times(1,isort)/times(2,isort),times(1,isort)/times(2,isort)/nthreads
     else
       nothers=nothers+1
       other_cpu=other_cpu+times(1,isort)
       other_wal=other_wal+times(2,isort)
     end if
     subcpu=subcpu+times(1,isort)
     subwal=subwal+times(2,isort)
   end do

   other_wal = other_wal + tol14
   write(entry_name,"(a,i0,a)")"others (",nothers,")"
   write(ount,format01041)entry_name,other_cpu,other_cpu*cpunm,other_wal,other_wal*wallnm,-1,-1.0, &
     other_cpu/other_wal,other_cpu/other_wal/nthreads

   write(ount,"(a)")"-<END_TIMER>"

   subwal = subwal + tol14
   write(ount,01201) subcpu,subcpu*cpunm,subwal,subwal*wallnm,subcpu/subwal,subcpu/subwal/nthreads

!  (2) Partitionings
   if (timopt<0) then

     npart=1000
     do ipart=1,npart
       list(:)=0
       select case(ipart)

       case(1)
         list(:11)=(/1,41,42,43,44,45,640,46,49,50,TIMER_SIZE/)      ; msg='abinit '
       case(2)
         list(:13)=(/640,641,642,700,132,84,301,401,501,650,643,644,TIMER_SIZE/)  ; msg='driver '
       case(3)
         list(:13)=(/700,703,704,705,33,701,34,35,36,706,702,707,708/)       ; msg='gstateimg+gstate '
       case(4)
         list(:19)=(/238,54,240,241,56,242,60,52,68,239,243,244,245,246,247,248,61,249,TIMER_SIZE/); msg='scfcv_core '
       case(5)
         list(:7)=(/940,941,942,943,944,945,TIMER_SIZE/)             ; msg= 'rhotov '
       case(6)
         list(:22)=(/980,981,982,983,984,28,985,271,986,987,988,989,990,991,992,993,994,995,996,997,1620,TIMER_SIZE/)
         msg= 'vtorho '
       case(7)
         list(:15)=(/28,31,22,530,585,583,590,222,572,842,537,586,591,1600,TIMER_SIZE/) ; msg='vtowfk '
       case(8)
         if(abs(timopt)==3)then
           list(:11)=(/530,204,205,571,532,533,630,535,536,584,587/)  ; msg='lobpcgwf (abs(timopt)==3)'
         else if(abs(timopt)==4)then
           list(:8)=(/530,520,521,522,523,524,525,526/)               ; msg='lobpcgwf (abs(timopt)==4)'
!            else
!            list(:3)=(/530,204,205/)
!            msg='lobpcgwf (light analysis: for a deeper one, use abs(timopt)=3 or 4)'
         end if
       case(9)
         list(:4)=(/22,201,40,211/)                            ; msg='cgwf '
       case(10)
         list(:8)=(/132,133,134,135,136,137,138,141/)          ; msg='respfn '
       case(11)
         list(:8)=(/141,142,143,144,120,146,147,TIMER_SIZE/)         ; msg='dfpt_looppert '
       case(12)
         list(:9)=(/120,154,121,157,152,158,160,150,564/) ; msg='dfpt_scfcv '
       case(13)
         list(:9)=(/121,118,128,126,287,166,129,127,556/)      ; msg='dfpt_vtorho '
       case(14)
         list(:9)=(/128,131,122,845,288,214,108,130,565/)      ; msg='dfpt_vtowfk '
       case(15)
         list(:8)=(/122,140,202,197,212,227,228,844/)          ; msg='dfpt_cgwf '
       case(16)
         list(:4)=(/200,841,221,98/)                           ; msg='getghc '
       case(17)
         list(:20)=(/801,840,841,842,843,844,845,846,847,848,849,850,851,852,853,854,855,856,857,858/)
         msg='fourwf (upwards partitioning)'
       case(18)
         list(:5)=(/933,934,936,937,938/)                      ; msg='outkss '
       case(19)
         list(:14)=(/301,302,315,316,319,304,305,320,321,306,307,308,309,310/)
         msg='screening '
       case(20)
         list(:13)=(/401,402,403,404,405,406,407,408,409,421,423,424,425/); msg='sigma  '
       case(21)
         list(:9)=(/431,432,433,434,435,445,440,441,442/)     ; msg='calc_sigc_me '
       case(23)
         list(:11)=(/630,631,632,633,634,545,635,636,637,638,TIMER_SIZE/)         ; msg='prep_getghc '
       case(24)
         list(:4)=(/539,856,547,548/)                          ; msg='prep_fourwf '
       case(25)
         list(:5)=(/570,231,232,581,TIMER_SIZE/)                     ; msg='prep_nonlop '
       case(26)
         list(:6)=(/790,791,792,793,794,795/)                  ; msg='mkrho (upwards partitioning)'
!          Disabled (temporarily ?) because the partitioning was not correct
!          case(27);list(:17)=(/600,601,602,603,604,605,617,606,607,608,609,610,611,612,613,614,615/)
!          msg='vtorhorec '
       case(28)
         list(:10)=(/650,651,653,654,655,656,658,659,660,661/)
         msg='bethe_salpeter '
       case(29)
         list(:8)=(/740,741,742,743,744,745,746,747/)          ; msg='suscep_stat '
       case(30)
         list(:9)=(/750,751,848,849,753,756,859,757,755/)      ; msg='susk '
       case(31)
         list(:8)=(/760,761,764,861,871,765,862,872/)          ; msg='suskmm '
       case(32)
         list(:8)=(/710,711,712,713,714,715,716,717/)          ; msg='inwffil '
       case(33)
         list(:10)=(/720,721,722,723,724,725,726,727,67,TIMER_SIZE/)  ; msg='wfsinp '
       case(34)
         list(:5)=(/770,771,772,272,290/)                      ; msg='initwf '
       case(35)
         list(:9)=(/780,781,782,783,784,785,786,291,292/)      ; msg='newkpt '
       case(36)
         list(:8)=(/93,901,902,903,904,905,268,TIMER_SIZE/)          ; msg='newvtr '
       case(37)
         list(:2)=(/94,269/)                                   ; msg='newrho '
       case(38)
         list(:11)=(/9,260,261,262,263,264,265,266,267,268,269/) ; msg=' fourdp (upwards partitioning)'
       case(39)
         list(:8)=(/250,251,252,253,254,255,256,257/)          ; msg='afterscfloop '
       case(40)
         list(:5)=(/910,911,912,913,914/)                      ; msg='forstr '
       case(41)
         list(:10)=(/920,921,927,922,923,926,924,65,925,TIMER_SIZE/) ; msg='forstrnps '
       case(42)
         list(:4)=(/670,671,672,673/)                          ; msg='exc_build_ham '
       case(43)
         list(:7)=(/680,681,682,683,684,685,686/)              ; msg='exc_build_block'
       case(44)
         list(:8)=(/690,691,692,693,694,695,696,697/)                  ; msg='exc_haydock_driver '
       case(45)
         list(:20)=(/950,951,952,953,954,955,956,957,958,959,960,961,962,963,964,965,966,967,968,969/)
         msg='outscfcv '
       case(46)
         list(:8)=(/620,621,622,623,624,625,626,627/)          ; msg='dmft '
       case(47)
         list(:9)=(/1001,1002,1003,1004,1005,1006,1007,1008,1009/)
         msg='initberry '
       case(50)
         list(:3)=(/1501,1502,1515/)                           ; msg='hartreefock '
       case(60)
         list(:13) = (/1600,1607,1630,1631,1632,1601,1603,1604,1605,1606,1608,1609,1610/)
         msg = 'chebfi'
       case(61)
         list(:3) = (/1620,1621,1622/)
         msg = 'mkinvovl'
       case(70)
         list(:5)=(/1701,1702,1703,1721,1722/)
         msg='gwls GW code'
       case(71)
         list(:16)=(/1703,1704,1705,1706,1707,1708,1709,1710,1711,1712,1713,1714,1715,1716,1717,1718/)
         msg='gwls: compute_correlations_shift_lanczos'
       case(72)
         list(:10)=(/1724,1725,1726,1727,1728,1729,1730,1731,1732,1733/)
         msg='gwls: Applying the susceptibility Pk'
       case(73)
         list(:7)=(/1734,1735,1736,1737,1738,1739,1740/)
         msg='gwls: Applying the model susceptibility Pk_model'
       case(74)
         list(:7)=(/1741,1742,1743,1744,1745,1746,1747/)
         msg='gwls: computing the matrix elements of eps_model^{-1}(w) -1 '
       case(75)
         list(:12)=(/1650,1651,1652,1653,1654,1655,1656,1657,1658,1659,1660,1661/)
         msg='lobpcgwf2 core engine '
       case(76)
         list(:18)=(/1670,1671,1672,1673,1674,1675,1676,1677,1678,1679,1680,1681,1682,1683,1684,1685,1686,1687/)
         msg='low-level xgBlock type '
       case(77)
         list(:5)=(/1690,1691,1692,1693,1694/)
         msg='low-level xgScalapack type '
       case(78)
         list(:8)=(/1662,1663,1664,1665,1666,1667,1668,1669/)
         msg='low-level xgTransposer type '
       case default
         cycle ! This allows to disable temporarily some partitionings

       end select

       nlist=0
       do itim=1,TIMER_SIZE
         if(list(itim)/=0)then
           nlist=nlist+1
         else
           exit
         end if
       end do

       if(nlist==0)then
         cycle
       end if

       if(ncount(list(1))/=0)then
         write(ount,'(/,a,a)')' Partitioning of ',trim(msg)
         subcpu=zero
         subwal=zero
         do ilist=1,nlist
           isort = list(ilist)
!          When the LAST item is TIMER_SIZE, a complement is evaluated (count number set to -1)
           if(ilist==nlist .and. list(nlist)==TIMER_SIZE)then
             times(1,TIMER_SIZE)=times(1,list(1))-subcpu
             times(2,TIMER_SIZE)=times(2,list(1))-subwal
             ncount(TIMER_SIZE)=-1
             ftimes(1,TIMER_SIZE)=zero
             mflops(TIMER_SIZE)=0
#if defined HAVE_TEST_TIME_PARTITIONING
             if(times(2,TIMER_SIZE)>1.2d0 .and. wallnm*times(2,TIMER_SIZE)>3.d0)then
               write(ount, '(3a,es16.6,4a,es16.6,2a)')&
                ' Note : the partitioning does not work well for this routine.',ch10,&
                '   The (other) Wall time            ',times(2,TIMER_SIZE),ch10,&
                '   is bigger than 1.2 secs. ',ch10,&
                '   The (other) Wall time percentage ',wallnm*times(2,TIMER_SIZE),ch10,&
                '   is bigger than 3% '
             else if (times(2,TIMER_SIZE)<0.2d0 .and. wallnm*times(2,TIMER_SIZE)<-0.2d0)then
               write(ount, '(3a,es16.6,2a)')&
                ' Note : the partitioning does not work well for this routine.',ch10,&
                '   The (other) Wall time percentage ',wallnm*times(2,TIMER_SIZE),ch10,&
                '   is negative '
             end if
#endif
           end if
           if(ncount(isort)/=0)then
             if(times(2,isort)*wallnm>0.02d0 .or. ilist==1)then   ! Does not write a slot if the wall time ratio is below a threshold
               if ( times(2,isort) < 0.0001 ) times(2,isort) = -1.d0
               write(ount,format01040)names(isort),&
                 times(1,isort),times(1,isort)*cpunm,&
                 times(2,isort),times(2,isort)*wallnm,ncount(isort), &
                 times(1,isort)/times(2,isort),times(1,isort)/times(2,isort)/nthreads
             end if
             if(ilist/=1)then
               subcpu=subcpu+times(1,isort)
               subwal=subwal+times(2,isort)
             else
               write(ount, '(a)' ) ' '
             end if
           end if
         end do

         subwal = subwal + tol14
         write(ount, 01201 ) subcpu,subcpu*cpunm,subwal,subwal*wallnm, subcpu/subwal,subcpu/subwal/nthreads
#ifdef HAVE_TEST_TIME_PARTITIONING
         if( wallnm*abs(subwal-times(2,list(1)))>1.d0 .and. abs(subwal-times(2,list(1)))>0.2d0 )then
           write(ount, '(3a,es16.6,2a,es16.6,4a,es16.6,2a,es16.6,6a,i4)')&
            ' Note : the partitioning does not work well for this routine.',ch10,&
            '   The subtotal Wall time            ',subwal,ch10,&
            '   differs from the total Wall time  ',times(2,list(1)),ch10,&
            '   by more than 0.2 secs.',ch10,&
            '   The subtotal Wall time percentage ',wallnm*subwal,ch10,&
            '   differs from the total Wall time %',wallnm*times(2,list(1)),ch10,&
            '   by more than 1%. ',ch10,&
            '   The partitioning might not have been coded properly.',ch10,&
            '   nlist=',nlist
           do ilist=1,nlist
             write(ount, '(a,i4,i4,es16.6,i8)' )&
              ' ilist,list(ilist),wallnm*times(2,list(ilist)),ncount(list(ilist))=',&
              ilist,isort,wallnm*times(2,isort),ncount(isort)
           end do
         end if
#endif
       end if

     end do ! End of loop on partitionings

!    For parallel case
     if(xmpi_paral==1)then
       write(ount, '(a,/,a)' )'-','-Synchronisation (=leave_test) and MPI calls '
       nlist=14
       list(:14)=(/48,61,62,63,64,65,66,67,71,85,86,543,544,787/)
       subcpu=zero; subwal=zero
       if(ncount(list(1))/=0)then
         do ilist=1,nlist
           isort = list(ilist)
!
           if (ncount(isort)/=0) then
             write(ount,format01040)names(isort),&
              times(1,isort),times(1,isort)*cpunm,&
              times(2,isort),times(2,isort)*wallnm,ncount(isort), &
              times(1,isort)/(tol14+times(2,isort)),times(1,isort)/(times(2,isort)+tol14)/nthreads

             if(ilist/=1)then
               subcpu=subcpu+times(1,isort)
               subwal=subwal+times(2,isort)
             else
               write(ount, '(a)' ) '-'
             end if
           end if !ncount
         end do !ilist

         subwal = subwal + tol14
         write(ount, 01200 ) subcpu,subcpu*cpunm,subwal,subwal*wallnm, subcpu/subwal,subcpu/subwal/nthreads
       end if !ncount
     end if !xmpi_paral

     nlist=23
     list(:23)=(/47,49,51,801,72,73,74,77,78,79,97,82,87,88,436,437,438,439,804,805,331,332,333/)
     flag_write=1
     do ilist=1,nlist
       isort = list(ilist)
       if(ncount(isort)/=0)then
         if(flag_write==1)then
           write(ount, '(/,a)' ) ' Additional information'
           flag_write=0
         end if
         write(ount,format01040)names(isort),&
           times(1,isort),times(1,isort)*cpunm,times(2,isort),times(2,isort)*wallnm,ncount(isort), &
           times(1,isort)/(tol14+times(2,isort)),times(1,isort)/(tol14+times(2,isort))/nthreads
       end if
     end do

     nlist=23
     list(:23)=(/550,551,552,553,554,555,556,558,559,560,561,562,563,564,565,566,567,593,594,595,596,597,598/)
     flag_write=1
     do ilist=1,nlist
       isort = list(ilist)
       if(ncount(isort)/=0)then
         if(flag_write==1)then
           write(ount, '(/,a)' ) ' Additional information about PAW segments'
           flag_write=0
         end if
         write(ount,format01040)names(isort),&
           times(1,isort),times(1,isort)*cpunm,times(2,isort),times(2,isort)*wallnm,ncount(isort), &
           times(1,isort)/(tol14+times(2,isort)),times(1,isort)/(tol14+times(2,isort))/nthreads
       end if
     end do

!    The detailed analysis cannot be done in the multidataset mode
     if(ndtset<2)then
       write(ount, '(/,/,a,/,a,/,a)' ) &
        ' Detailed analysis of some time consuming routines ',&
        '                                  tcpu    ncalls  tcpu/ncalls    ndata tcpu/ncalls/ndata',&
        '                                 (sec)                (msec)              (microsec)'
       nlist=8
       list(:8)=(/802,803,9,75,76,77,210,11/)
       do ilist=1,nlist
         isort = list(ilist)
         if(ncount(isort)/=0)then
           write(ount, '(a,a24,f12.3,i10,f12.3,i10,f12.3)' )'- ',names(isort),&
             times(1,isort),ncount(isort),&
             1000.0_dp*times(1,isort)/dble(ncount(isort)),ndata(isort),&
             1000000.0_dp*times(1,isort)/dble(ncount(isort)*dble(ndata(isort)))
         else
           write(ount, '(a,a24,f12.3,i10)' )'- ',names(isort),times(1,isort),ncount(isort)
         end if
       end do !ilist
     else
       write(ount,'(/,a)') ' timana : in multi dataset mode, the more detailed analysis is not done.'
     end if !ndtset

   end if ! End the condition of timopt<0

 end if ! me==0

 ABI_FREE(list)

end subroutine timana
!!***

end module m_timana
!!***
