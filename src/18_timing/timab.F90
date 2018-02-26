!{\src2tex{textfont=tt}}
!!****f* ABINIT/timab
!! NAME
!!  timab
!!
!! FUNCTION
!!  Timing subroutine.  Calls machine-dependent "timein" which returns elapsed cpu and wall clock times in sec.
!!
!!  Depending on value of "option" routine will:
!!  (0) zero all accumulators
!!  (1) start with new incremental time slice for accumulator n using explicit call to timein (or PAPI)
!!  (2) stop time slice; add time to accumulator n also increase by one the counter for this accumulator
!!  (3) start with new incremental time slice for accumulator n
!!        using stored values for cpu, wall, and PAPI infos ( ! do not use for stop )
!!  (4) report time slice for accumlator n (not full time accumlated)
!!  (5) option to suppress timing (nn should be 0) or reenable it (nn /=0)
!!
!!  If, on first entry, subroutine is not being initialized, it
!!  will automatically initialize as well as rezero accumulator n.
!!  However, initialization SHOULD be done explicitly by the user
!!  so that it can be done near the top of his/her main routine.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nn=index of accumulator (distinguish what is being timed); NOT used if option=0
!!  option=see comment above
!!
!! OUTPUT
!!  on option=4 :
!!    tottim(2,nn)=accumulated time for accumulator nn; otherwise
!!     tottim is a dummy variable.
!!    option gives the number of times that the
!!     accumulator has been incremented
!!
!! PARENTS
!!      abinit,afterscfloop,atm2fft,bethe_salpeter,calc_sigc_me,calc_sigx_me
!!      calcdensph,cchi0,cgq_builder,cgwf,chebfi,cohsex_me,corrmetalwf1,d2frnl
!!      density_rec,dfpt_cgwf,dfpt_dyfro,dfpt_dyxc1,dfpt_eltfrhar,dfpt_eltfrkin
!!      dfpt_eltfrloc,dfpt_eltfrxc,dfpt_ewald,dfpt_looppert,dfpt_mkrho
!!      dfpt_mkvxc,dfpt_mkvxc_noncoll,dfpt_mkvxcstr,dfpt_newvtr,dfpt_nstdy
!!      dfpt_nstpaw,dfpt_nstwf,dfpt_rhofermi,dfpt_rhotov,dfpt_scfcv,dfpt_vtorho
!!      dfpt_vtowfk,dfpt_wfkfermi,dfptnl_loop,dielmt,dieltcel,dmft_solve
!!      dotprodm_v,dotprodm_vn,driver,dyson,eig2stern,eig2tot,elt_ewald
!!      eltxccore,energy,entropyrec,etotfor,exc_build_block,exc_build_ham
!!      fermisolverec,first_rec,fock2ACE,fock_getghc,forces,forstr,forstrnps
!!      fourdp,fourwf,fxphas,getgh1c,getghc,getgsc,getngrec,gran_potrec
!!      green_kernel,gstate,gstateimg,gwls_ComputeCorrelationEnergy
!!      gwls_DielectricArray,gwls_QR_factorization,gwls_lineqsolver
!!      gwls_model_polarisability,gwls_polarisability,gwls_sternheimer,hartre
!!      impurity_solve,initberry,initorbmag,initwf,inkpts,invars2,inwffil
!!      listkk,lobpcgwf,m_ab7_invars_f90,m_ab7_mixing,m_cgtools,m_dyson_solver
!!      m_fftcore,m_fftw3,m_fock,m_green,m_haydock,m_hexc,m_invovl,m_iowf
!!      m_lobpcg,m_lobpcg2,m_lobpcgwf,m_paral_pert,m_sg2002,m_wfutils,m_xg
!!      mag_constr,mkcore,mkcore_paw,mkcore_wvl,mkffnl,mklocl_realspace
!!      mklocl_recipspace,mkresi,mkrho,newkpt,newocc,newrho,newvtr,nhatgrid
!!      nlenergyrec,nonlinear,nonlop,odamix,opernla_ylm,optics_paw
!!      optics_paw_core,optics_vloc,outkss,outscfcv,pareigocc
!!      partial_dos_fractions_paw,pawdenpot,pawdfptenergy,pawinit,pawmknhat
!!      pawmknhat_psipsi,pawmkrho,pawpolev,prep_bandfft_tabs,prep_calc_ucrpa
!!      prep_fourwf,prep_getghc,prep_nonlop,pspatm,pspheads_comm,pspini
!!      pw_orthon,rayleigh_ritz,recursion,recursion_nl,respfn,rhotov,rhotoxc
!!      rwwf,scfcv,screening,setsym,setvtr,sigma,sqnormm_v,status,stress,strhar
!!      suscep_stat,susk,suskmm,symrhg,symsgcube,tddft,timana,vn_nl_rec,vtorho
!!      vtorhorec,vtorhotf,vtowfk,wf_mixing,wfconv,wfk_analyze,wfsinp
!!      wvl_nhatgrid,xcden,xcpot
!!
!! CHILDREN
!!      papif_flops,papif_perror,timein
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine timab(nn,option,tottim)

 use defs_basis
 use defs_time
 use m_profiling_abi
 use m_errors
 use iso_c_binding

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'timab'
 use interfaces_18_timing, except_this_one => timab
!End of the abilint section

 implicit none

#ifdef HAVE_PAPI
#include "f90papi.h"
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn,option
!arrays
 real(dp),intent(out) :: tottim(2)

!Local variables-------------------------------
!scalars
 real(dp),save :: cpu,wall
 character(len=500) :: message
#ifdef HAVE_PAPI
 integer(C_INT) :: check 
 integer(C_LONG_LONG),save :: flops1
 real(C_FLOAT),save :: real_time,proc_time
 real(C_FLOAT) :: mflops1
 character(len=PAPI_MAX_STR_LEN) :: papi_errstr
#endif
! *************************************************************************

 if (option==5) timopt=nn

!If timopt was set to zero by a call with option=5, suppress
!all action of this routine (might as well return at this point !)
 if(timopt/=0 .and. option/=5)then
!  
!  Check that nn lies in sensible bounds
   if (nn<1.or.nn>mtim) then
     write(message,'(a,i0,a,i0)')'  mtim = ',mtim,' but input nn = ',nn
     MSG_BUG(message)
   end if

#ifdef HAVE_PAPI
!  for all active options for time if papi analysis has been selected.
   if (option/=3.and.time_get_papiopt()==1) then 
     call PAPIf_flops(real_time, proc_time, flops1, mflops1, check)
     if (check /= PAPI_OK) then
       call papif_perror(check,papi_errstr,check)
       write(std_out,*) 'Problem to initialize papi high level inteface'
       write(std_out,*) 'Error code', papi_errstr
     end if
     if (flops1 < 0) then  
       MSG_WARNING("Number of floating point instruction Overflow")
       papi_flops(:)=-1            
     end if
   end if
#endif
   
   select case (option)
   case (0)  
       ! Zero out all accumulators of time and init timers
     acctim(:,:)      = 0.0d0
     tzero(:,:)       = 0.0d0
     ncount(:)        = 0
     papi_flops(:)    = 0
     papi_acctim(:,:) = 0. 
     papi_accflops(:) = 0. 
     papi_tzero(:,:)  = 0. 

   case (1)  
       ! Initialize timab for nn
     call timein(cpu,wall)
     tzero(1,nn)=cpu
     tzero(2,nn)=wall
#ifdef HAVE_PAPI
     papi_flops(nn)   = flops1       ! Initialize megaflops for nn
     papi_tzero(1,nn) = proc_time
     papi_tzero(2,nn) = real_time
#endif

   case (2)  
       ! Accumulate time for nn (also keep the values of cpu, wall, proc_time, real_time, flops1)
     call timein(cpu,wall)
     acctim(1,nn)=acctim(1,nn)+cpu -tzero(1,nn)
     acctim(2,nn)=acctim(2,nn)+wall-tzero(2,nn)
     ncount(nn)=ncount(nn)+1
#ifdef HAVE_PAPI
!      accumulate time and flops for nn Difference between 2 calls to Papif_flops 
     papi_acctim(1,nn)=papi_acctim(1,nn)+ proc_time - papi_tzero(1,nn)
     papi_acctim(2,nn)=papi_acctim(2,nn)+ real_time - papi_tzero(2,nn)
     papi_accflops(nn)=papi_accflops(nn)+ flops1- papi_flops(nn) 
#endif

   case (3) 
       ! Use previously obtained values to initialize timab for nn
     tzero(1,nn)=cpu
     tzero(2,nn)=wall
#ifdef HAVE_PAPI
     papi_flops(nn)=flops1
     papi_tzero(1,nn) = proc_time
     papi_tzero(2,nn) = real_time
#endif

   case (4) 
       ! Return elapsed time for nn (do not accumulate)
     call timein(cpu,wall)
     tottim(1)=cpu-tzero(1,nn)
     tottim(2)=wall-tzero(2,nn)
#ifdef HAVE_PAPI
!      return elapsed floating point operationfor nn (do not accumulate)
     papi_tottim(1,nn)= proc_time - papi_tzero(1,nn)
     papi_tottim(2,nn)= real_time - papi_tzero(2,nn)
     papi_totflops(nn)= flops1 - papi_flops(nn) 
#endif

   case default
     write(message,'(a,i10,a)')'  Input option not valid, =',option,'.'
     MSG_BUG(message)
   end select 
 end if

end subroutine timab
!!***
