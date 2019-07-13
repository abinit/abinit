!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_common
!! NAME
!!  m_common
!!
!! FUNCTION
!!  This module gathers routines used by higher-level procedures.
!!  Mainly printing routines.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, AF, GMR, LBoeri, MT)
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

module m_common

 use defs_basis
 use m_errors
 use m_abicore
 use m_exit
 use m_fftcore
 use m_fock
 use m_io_tools
#if defined DEV_YP_VDWXC
 use m_xc_vdw
#endif
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_nctk
 use m_crystal
 use m_wfk
 use m_ebands
 use m_hdr
 use m_xmpi
 use m_dtset
 use m_xpapi
 use m_stream_string
 use m_yaml_out
 use m_invars2

 use m_fstrings,          only : indent, endswith, sjoin
 use m_electronpositron,  only : electronpositron_type
 use m_energies,          only : energies_type, energies_eval_eint
 use m_pair_list,         only : pair_list
 use m_neat,              only : neat_open_etot, neat_finish_etot, neat_etot_add_line
 use m_geometry,          only : mkrdim, metric
 use m_kg,                only : getcut
 use m_parser,            only : parsefile
 use m_invars1,           only : invars0, invars1m, indefo
 use m_time,              only : timab, time_set_papiopt
 use defs_abitypes,       only : dataset_type, ab_dimensions, hdr_type, MPI_type
 use defs_datatypes,      only : pspheader_type, ebands_t
 use m_pspheads,          only : inpspheads, pspheads_comm

 implicit none

 private
!!***

 public :: scprqt
 public :: setup1
 public :: prteigrs
 public :: prtene
 public :: get_dtsets_pspheads     ! Parse input file, get list of pseudos for files file and build list of datasets
                                   ! pseudopotential headers, maxval of dimensions needed in outvars
 public :: ebands_from_file        ! Build an ebands_t object from file. Supports Fortran and netcdf files
 public :: crystal_from_file       ! Build a crystal_t object from netcdf file or Abinit input file
                                   ! with file extension in [".abi", ".in"]

 type(stream_string),private,save :: etot_yaml_doc
!!***

contains
!!***

!!****f* ABINIT/scprqt
!! NAME
!! scprqt
!!
!! FUNCTION
!! Conducts printing inside the scfcv.F90 routine, according to the value of choice.
!! Also checks the convergence with respect to the different criteria.
!! Eventually send a signal to quit the SCF cycle.
!!
!! INPUTS
!!  choice= if 1 => called at the initialisation of scfcv.f
!!          if 2 => called during the loop in scfcv.f
!!          if 3 => called at the end of scfcv.f
!!  cpus=cpu time limit in seconds
!!  deltae=change in energy between the previous and present SCF cycle
!!  diffor=maximum absolute change in component of fcart between present
!!          and previous SCF cycle.
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | chkexit= if non-zero, check whether the user wishes to exit
!!   | enunit=parameter determining units of output energies
!!   | ionmov=governs the movement of atoms (see help file)
!!   | kptopt=option for the generation of k points
!!   | mband=maximum number of bands
!!   | natom=number of atoms in cell.
!!   | nnsclo_now=number of non-self-consistent loops for the current vtrial
!!   |  (often 1 for SCF calculation, =nstep for non-SCF calculations)
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | occopt=option for occupancies
!!   | prtxml=1 if values have to be stored in an XML file.
!!   | prteig=
!!   | prtstm=print STM input variable
!!   | prtvol= control print volume
!!   | usedmatpu=LDA+U: number of SCF steps keeping occ. matrix fixed
!!   | usefock=1 if Fock operator is present (hence possibility of a double loop)
!!   | usepawu=0 if no LDA+U; /=0 if LDA+U
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  etotal=total energy (hartree)
!!  favg(3)=average of forces (ha/bohr)
!!  fcart(3,natom)=cartesian forces (hartree/bohr)
!!  fermie=fermi energy (Hartree)
!!  fname_eig=filename for printing of the eigenenergies
!!  fock <type(fock_type)>=quantities for the fock operator (optional argument)
!!  character(len=fnlen) :: filnam1=character strings giving input file name
!!  initGS= 1 if one GS SCF cycle has already be done
!!  iscf=( <= 0 =>non-SCF), >0 => SCF)
!!   iscf =1 => determination of the largest eigenvalue of the SCF cycle
!!   iscf =2 => SCF cycle, simple mixing
!!   iscf =3 => SCF cycle, anderson mixing
!!   iscf =5 => SCF cycle, CG based on estimations of gradients of the energy
!!   iscf =6 => SCF cycle, CG based on true minimization of the energy
!!   iscf =-3, although non-SCF, the energy is computed, so print it here.
!!  istep=number of the SCF iteration (needed if choice=2)
!!  istep_fock_outer=number of outer SCF iteration in the double loop approach
!!  istep_mix=number of inner SCF iteration in the double loop approach
!!  kpt(3,nkpt)=reduced coordinates of k points.
!!  maxfor=maximum absolute value of fcart
!!  moved_atm_inside: if==1, the atoms are allowed to move.
!!  mpi_enreg=information about MPI parallelization
!!  nband(nkpt*nsppol)=number of bands at each k point, for each polarization
!!  nkpt=number of k points
!!  nstep=number of steps expected in iterations.
!!  occ(mband*nkpt*nsppol)=occupation number for each band at each k point.
!!  optres=0 if the residual (res2) is a POTENTIAL residual
!!         1 if the residual (res2) is a DENSITY residual
!!  prtfor=1 only if forces have to be printed (0 otherwise)
!!  prtxml=1 if XML file has to be output
!!  res2=square of the density/potential residual
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points and spins
!!  residm=maximum value from resid array (except for nbdbuf highest bands)
!!         in Wavelets mode, it is used as the maximum value for the gradient norm.
!!  response= if 0, GS case, if 1, RF case.
!!  tollist(12)=tolerance list. Presently, the following are defined :
!!    tollist(1)=tolmxf ; tollist(2)=tolwfr ; tollist(3)=toldff
!!    tollist(4)=toldfe ; tollist(5)=toleig ; tollist(6)=tolvrs
!!    tollist(7)=tolrff
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vxcavg=mean of the vxc potential
!!  wtk(nkpt)=weight assigned to each k point.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  quit= 0 if the SCF cycle is not finished; 1 otherwise.
!!  conv_retcode=Only if choice==3, != 0 if convergence is not achieved.
!!
!! PARENTS
!!      afterscfloop,dfpt_scfcv,scfcv
!!
!! CHILDREN
!!      exit_check,flush_unit,prteigrs,wrtout,xc_vdw_trigger
!!
!! SOURCE

subroutine scprqt(choice,cpus,deltae,diffor,dtset,&
&  eigen,etotal,favg,fcart,fermie,fname_eig,filnam1,initGS,&
&  iscf,istep,istep_fock_outer,istep_mix,kpt,maxfor,moved_atm_inside,mpi_enreg,&
&  nband,nkpt,nstep,occ,optres,&
&  prtfor,prtxml,quit,res2,resid,residm,response,tollist,usepaw,&
&  vxcavg,wtk,xred,conv_retcode,&
&  electronpositron, fock) ! optional arguments)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,initGS,iscf,istep,istep_fock_outer,istep_mix
 integer,intent(in) :: moved_atm_inside,nkpt,nstep
 integer,intent(in) :: optres,prtfor,prtxml,response,usepaw
 integer,intent(out) :: quit,conv_retcode
 real(dp),intent(in) :: cpus,deltae,diffor,etotal,fermie,maxfor,res2,residm
 real(dp),intent(in) :: vxcavg
 character(len=fnlen),intent(in) :: fname_eig,filnam1
 type(electronpositron_type),pointer,optional :: electronpositron
 type(fock_type),pointer,optional :: fock
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: nband(nkpt*dtset%nsppol)
 real(dp),intent(in) :: eigen(dtset%mband*nkpt*dtset%nsppol),favg(3)
 real(dp),intent(in) :: fcart(3,dtset%natom),kpt(3,nkpt)
 real(dp),intent(in) :: occ(dtset%mband*nkpt*dtset%nsppol)
 real(dp),intent(in) :: resid(dtset%mband*nkpt*dtset%nsppol),tollist(12)
 real(dp),intent(in) :: wtk(nkpt),xred(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer,save :: toldfe_ok,toldff_ok,tolrff_ok,ttoldfe,ttoldff,ttolrff,ttolvrs,ttolwfr
 integer :: iatom,iband,iexit,ikpt,ii,ishift,isppol,my_rank
 integer :: nband_index,nband_k,nnsclohf
 integer :: openexit,option,tmagnet,usefock
#if defined DEV_YP_VDWXC
 integer :: ivdw
#endif
 real(dp),save :: toldfe,toldff,tolrff,tolvrs,tolwfr,vdw_df_threshold
 real(dp) :: diff_e,diff_f,magnet,rhodn,rhoup
 logical :: noquit,use_dpfft
 character(len=500) :: message, message2, message3
 character(len=2) :: format_istep
 character(len=5) :: format_magnet
 character(len=8) :: colname
 character(len=1) :: firstchar
!arrays
 real(dp) :: residm_band(dtset%mband,dtset%nsppol)
 real(dp) :: f_tmp(3)

! *********************************************************************

 DBG_ENTER("COLL")

 my_rank = mpi_enreg%me_cell

 quit=0; conv_retcode=0
 usefock=dtset%usefock
 nnsclohf=dtset%nnsclohf
 use_dpfft = .False.

 tmagnet=0
 if(response==0.and.(iscf>0.or.iscf==-3).and.dtset%nsppol==2.and.dtset%occopt>2)tmagnet=1

 ishift=0
 residm_band = zero
 do isppol=1, dtset%nsppol
   do ikpt=1, nkpt
     do iband=1, nband(ikpt+(isppol-1)*nkpt)
       ishift = ishift+1
       residm_band(iband, isppol) = max (resid(ishift), residm_band(iband, isppol))
     end do
   end do
 end do

 select case (choice)

 case (1)
   ! Examine tolerance criteria
   ! NB: The tests on tolwfr and the presence of tolerances in the SCF case are
   ! also done at the level of the parser in chkinp.
   tolwfr=tollist(2)
   toldff=tollist(3)
   toldfe=tollist(4)
   tolvrs=tollist(6)
   tolrff=tollist(7)
   vdw_df_threshold=tollist(8)
   ttolwfr=0 ; ttoldff=0 ; ttoldfe=0 ; ttolvrs=0; ttolrff=0;
   if(abs(tolwfr)>tiny(zero))ttolwfr=1
   if(abs(toldff)>tiny(zero))ttoldff=1
   if(abs(tolrff)>tiny(zero))ttolrff=1
   if(abs(toldfe)>tiny(zero))ttoldfe=1
   if(abs(tolvrs)>tiny(zero))ttolvrs=1
   !  If non-scf calculations, tolwfr must be defined
   if(ttolwfr /= 1 .and. (iscf<0 .and. iscf/=-3) )then
     write(message,'(a,a,a,es14.6,a,a)')&
&     'when iscf <0 and /= -3, tolwfr must be strictly',ch10,&
&     'positive, while it is ',tolwfr,ch10,&
&     'Action: change tolwfr in your input file and resubmit the job.'
     MSG_ERROR(message)
   end if
   ! toldff only allowed when prtfor==1
   ! FIXME: this test should be done on input, not during calculation
   if((ttoldff == 1 .or. ttolrff == 1) .and. prtfor==0 )then
     MSG_ERROR('toldff only allowed when prtfor=1!')
   end if
   ! If SCF calculations, one and only one of these can differ from zero
   if(ttolwfr+ttoldff+ttoldfe+ttolvrs+ttolrff /= 1 .and. (iscf>0 .or. iscf==-3))then
     write(message,'(6a,es14.6,a,es14.6,a,es14.6,a,es14.6,a,a,es14.6,a,a,a)' )&
&     'For the SCF case, one and only one of the input tolerance criteria ',ch10,&
&     'tolwfr, toldff, tolrff, toldfe or tolvrs ','must differ from zero, while they are',ch10,&
&     'tolwfr=',tolwfr,', toldff=',toldff,', tolrff=',tolrff,', toldfe=',toldfe,ch10,&
&     'and tolvrs=',tolvrs,' .',ch10,&
&     'Action: change your input file and resubmit the job.'
     MSG_ERROR(message)
   end if

   if (dtset%usewvl == 1) then
     write(colname, "(A)") "grdnorm "
   else
     write(colname, "(A)") "residm  "
   end if
   if (nstep>0 .and. (iscf>=0 .or.iscf==-3) .and. dtset%prtstm==0) then
     if(tmagnet==1)then
       if (prtfor==0) then
         if (optres==0) then
           write(message, '(4a)' ) ch10,&
&           '     iter   Etot(hartree)     deltaE(h) ',colname,'  vres2    magn'
         else
           write(message, '(4a)' ) ch10,&
&           '     iter   Etot(hartree)     deltaE(h) ',colname,'  nres2    magn'
         end if
       else
         if (optres==0) then
           write(message, '(4a)' ) ch10,&
&           '     iter   Etot(hartree)     deltaE(h) ',colname,'  vres2   diffor   maxfor   magn'
         else
           write(message, '(4a)' ) ch10,&
&           '     iter   Etot(hartree)     deltaE(h) ',colname,'  nres2   diffor   maxfor   magn'
         end if
       end if
     else
       if(response==0)then
         if (prtfor==0) then
           if (optres==0) then
             write(message, '(4a)' ) ch10,&
&             '     iter   Etot(hartree)      deltaE(h)  ', colname, '   vres2'
           else
             write(message, '(4a)' ) ch10,&
&             '     iter   Etot(hartree)      deltaE(h)  ', colname, '   nres2'
           end if
         else
           if (optres==0) then
             write(message, '(4a)' ) ch10,&
&             '     iter   Etot(hartree)      deltaE(h)  ',colname,'   vres2    diffor    maxfor '
           else
             write(message, '(4a)' ) ch10,&
&             '     iter   Etot(hartree)      deltaE(h)  ',colname,'   nres2    diffor    maxfor '
           end if
         end if
       else
         if (optres==0) then
           write(message, '(4a)' ) ch10,&
&           '     iter   2DEtotal(Ha)        deltaE(Ha) ', colname, '  vres2'
         else
           write(message, '(4a)' ) ch10,&
&           '     iter   2DEtotal(Ha)        deltaE(Ha) ', colname, '  nres2'
         end if
       end if
     end if
     ! Will save iterations in this global variables.
     call neat_open_etot(etot_yaml_doc, '', message)
     call wrtout(ab_out,message,'COLL')
   end if

 case (2)

   ! Examine tolerance criteria
   tolwfr=tollist(2)
   toldff=tollist(3)
   toldfe=tollist(4)
   tolvrs=tollist(6)
   tolrff=tollist(7)
   vdw_df_threshold=tollist(8)
   ttolwfr=0 ; ttoldff=0 ; ttoldfe=0 ; ttolvrs=0; ttolrff=0;
   if(abs(tolwfr)>tiny(0.0_dp))ttolwfr=1
   if(abs(toldff)>tiny(0.0_dp))ttoldff=1
   if(abs(tolrff)>tiny(0.0_dp))ttolrff=1
   if(abs(toldfe)>tiny(0.0_dp))ttoldfe=1
   if(abs(tolvrs)>tiny(0.0_dp))ttolvrs=1
   ! Conduct printing. If extra output follows, then put a blank line into the output here
   if (dtset%prtvol>=10) then
     call wrtout([std_out, ab_out], ' ')
   end if

   ! Calculate up and down charge and magnetization
   if(tmagnet==1) then
     rhoup = zero
     rhodn = zero
     nband_index = 1
     do isppol=1,dtset%nsppol
       do ikpt=1,nkpt
         nband_k=nband(ikpt+(isppol-1)*nkpt)
         do iband=1,nband_k
           if(isppol==1) rhoup = rhoup + wtk(ikpt)*occ(nband_index)
           if(isppol==2) rhodn = rhodn + wtk(ikpt)*occ(nband_index)
           nband_index = nband_index + 1
         end do
       end do
     end do
     magnet = abs(rhoup - rhodn)
   end if

   if (prtxml == 1) then
     write(ab_xml_out, "(A)", advance = "NO") '      <scfcvStep'
     write(message, "(es22.10)") etotal
     write(ab_xml_out, "(A,A,A)", advance = "NO") ' eTotal="', trim(message) ,'"'
     write(message, "(es20.8)") deltae
     write(ab_xml_out, "(A,A,A)", advance = "NO") ' deltaETotal="', trim(message) ,'"'
     write(message, "(es20.8)") residm
     write(ab_xml_out, "(A,A,A)", advance = "NO") ' maxResid="', trim(message) ,'"'
     write(message, "(es20.8)") res2
     if (optres == 0) then
       write(ab_xml_out, "(A,A,A)", advance = "NO") ' potResid="', trim(message) ,'"'
     else
       write(ab_xml_out, "(A,A,A)", advance = "NO") ' denResid="', trim(message) ,'"'
     end if
     if (tmagnet== 1) then
       write(message, "(es20.8)") magnet
       write(ab_xml_out, "(A,A,A)", advance = "NO") ' magn="', trim(message) ,'"'
     end if
     if (prtfor == 1) then
       write(message, "(es20.8)") diffor
       write(ab_xml_out, "(A,A,A)", advance = "NO") ' deltaForces="', trim(message) ,'"'
       write(message, "(es20.8)") maxfor
       write(ab_xml_out, "(A,A,A)", advance = "NO") ' maxForces="', trim(message) ,'"'
     end if
     write(ab_xml_out, "(A)") " />"
   end if

   ! Print total (free) energy (hartree) and other convergence measures
   if(dtset%prtstm==0)then
     format_istep='i3'
     if(istep>99)format_istep='i5'
     if(istep>9999)format_istep='i7'
     if(tmagnet==1)then
       if(magnet<10)then
         format_magnet='f6.3)'
       else if(magnet<100)then
         format_magnet='f6.2)'
       else
         format_magnet='f6.1)'
       end if
       if (prtfor==0) then
         write(message, '(a,'//format_istep//',1p,g22.14,3es9.2,0p,'//format_magnet ) &
&         ' ETOT',istep,etotal,deltae,residm,res2,magnet
       else
         write(message, '(a,'//format_istep//',1p,g22.14,3es9.2,es8.1,es9.2,0p,'//format_magnet ) &
&         ' ETOT',istep,etotal,deltae,residm,res2,diffor,maxfor,magnet
       end if
     else
       firstchar=' '
       if (response/=0.and.istep==1) firstchar="-"
       if (response==0) then
         if (prtfor==0) then
           write(message, '(2a,'//format_istep//',1p,g22.14,3es10.3)' ) &
&           firstchar,'ETOT',istep,etotal,deltae,residm,res2
         else
           write(message, '(2a,'//format_istep//',1p,g22.14,5es10.3)' ) &
&           firstchar,'ETOT',istep,etotal,deltae,residm,res2,diffor,maxfor
         end if
       else
         write(message, '(2a,'//format_istep//',1p,g22.14,1x,3es10.3)' ) &
&         firstchar,'ETOT',istep,etotal,deltae,residm,res2
       end if
     end if
     if (etot_yaml_doc%length /= 0) then
       call neat_etot_add_line(etot_yaml_doc, message)
     end if
     call wrtout(ab_out,message,'COLL')

     if(mpi_enreg%paral_pert==1) then
       call wrtout(std_out,  message,'PERS')
     elseif(mpi_enreg%paral_pert==0) then
       call wrtout(std_out,  message,'COLL')
     end if

   end if ! dtset%prtstm==0

   ! Print positions/forces every step if dtset%prtvol>=10 and iscf>0 or -3 and GS case
   if (dtset%prtvol>=10.and.(iscf>=0.or.iscf==-3).and.response==0.and.dtset%prtstm==0) then
     call wrtout(ab_out," ",'COLL')

     ! Print up and down charge and magnetization
     if(tmagnet==1) then
       write(message,'(a,f11.6,a,f11.6,a,f10.6)')&
&       ' #electrons spin up=',rhoup,', spin down=',rhodn,', magnetization=',magnet
       call wrtout([std_out, ab_out], message)
     end if

     ! Moreover, print atomic positions if dtset%ionmov==4, and moved_atm_inside==1
     if (dtset%ionmov==4 .and. moved_atm_inside==1)then
       call wrtout([std_out, ab_out], ' reduced coordinates :')
       do iatom=1,dtset%natom
         write(message, '(i5,1x,3es21.11)' ) iatom,xred(:,iatom)
         call wrtout([std_out, ab_out], message)
       end do
     end if

     ! Slightly change favg for printing reasons
     if (prtfor>0) then
       f_tmp(:)=favg(:)
       if(abs(favg(1))<1.0d-13)f_tmp(1)=zero
       if(abs(favg(2))<1.0d-13)f_tmp(2)=zero
       if(abs(favg(3))<1.0d-13)f_tmp(3)=zero
       write(message, '(a,3es10.2)' )' cartesian forces (ha/bohr); non-corrected avg=',f_tmp(:)
       call wrtout([std_out, ab_out], message)
       do iatom=1,dtset%natom
         f_tmp(:)=fcart(:,iatom)
         if(abs(fcart(1,iatom))<1.0d-13)f_tmp(1)=zero
         if(abs(fcart(2,iatom))<1.0d-13)f_tmp(2)=zero
         if(abs(fcart(3,iatom))<1.0d-13)f_tmp(3)=zero
         write(message, '(i5,1x,3es21.11)' ) iatom,f_tmp(:)
         call wrtout([std_out, ab_out], message)
       end do
     end if

   end if

   ! Print eigenvalues every step if dtset%prtvol>=10 and GS case
   if (my_rank == master .and. (dtset%prtvol>=10 .and. response==0 .and. dtset%tfkinfunc==0 .and. dtset%usewvl==0)) then
     option=1
     call prteigrs(eigen,dtset%enunit,fermie,fname_eig,ab_out,iscf,kpt,dtset%kptopt,dtset%mband,&
&     nband,nkpt,dtset%nnsclo,dtset%nsppol,occ,dtset%occopt,option,dtset%prteig,dtset%prtvol,resid,tolwfr,vxcavg,wtk)

     call prteigrs(eigen,dtset%enunit,fermie,fname_eig,std_out,iscf,kpt,dtset%kptopt,dtset%mband,&
&     nband,nkpt,dtset%nnsclo,dtset%nsppol,occ,dtset%occopt,option,dtset%prteig,dtset%prtvol,resid,tolwfr,vxcavg,wtk)
   end if

   if(response==0)then
     write(message, '(a,1p,e15.7,a)'  ) ' scprqt: <Vxc>=',vxcavg,' hartree'
     call wrtout(std_out,message,'COLL')
   end if

   ! Check whether exiting was required by the user.
   openexit=1 ; if(dtset%chkexit==0) openexit=0
   call exit_check(cpus,filnam1,iexit,ab_out,mpi_enreg%comm_cell,openexit)
   if (iexit/=0) quit=1

   ! In special cases, do not quit even if convergence is reached
   noquit=((istep<nstep).and.(usepaw==1).and.(dtset%usepawu/=0).and.&
&   (dtset%usedmatpu/=0).and.(istep<=abs(dtset%usedmatpu)).and.&
&   (dtset%usedmatpu<0.or.initGS==0))

   ! Additional stuff for electron/positron
   if (present(electronpositron)) then
     if (associated(electronpositron)) then
       if (electronpositron%istep_scf==1) then
         toldff_ok=0;tolrff_ok=0;toldfe_ok=0
       end if
     end if
   end if

   ! Stopping criteria in the SCF case
   if(iscf>1 .or. iscf==-3 .or. iscf == 0) then
     ! Here treat the vdw_df_threshold criterion : if the change of energy is less than
     ! input vdw_df_threshold, trigger the calculation of vdW interactions
     ! write(message,'(1x,a,e10.3,1x,a,e10.3,1x,l1,a)') &
     ! &      '[vdW-DF][DEBUG] deltae=',deltae,'vdw_df_threshold=',vdw_df_threshold, &
     ! &      (abs(deltae)<vdw_df_threshold),ch10
     ! call wrtout(std_out,message,'COLL')
#if defined DEV_YP_VDWXC
     call xc_vdw_trigger( (abs(deltae)<vdw_df_threshold) )
#endif
     ! Here treat the tolwfr criterion: if maximum residual is less than
     ! input tolwfr, stop steps (exit loop here)
     if (ttolwfr == 1 .and. .not. noquit) then
       if (residm < tolwfr) then
         if (dtset%usewvl == 0) then
           write(message, '(a,a,i5,a,1p,e10.2,a,e10.2,a,a)' )ch10, &
           ' At SCF step',istep,'   max residual=',residm,' < tolwfr=',tolwfr,' =>converged.'
         else
           write(message, '(a,a,i5,a,1p,e10.2,a,e10.2,a,a)' )ch10, &
           ' At SCF step',istep,'   max grdnorm=',residm,' < tolwfr=',tolwfr,' =>converged.'
         end if
         call wrtout([std_out, ab_out], message)
         quit=1
       else
         use_dpfft = residm < tol7
       end if
     end if

     ! Here treat the toldff criterion: if maximum change of fcart is less than
     ! input toldff twice consecutively, stop steps (exit loop here)
     if (ttoldff==1) then
       if (istep==1) then
         toldff_ok=0
       else if (diffor < toldff) then
         toldff_ok=toldff_ok+1
         ! add warning for forces which are 0 by symmetry. Also added Matteo check below that the wave
         ! functions are relatively converged as well
         if (diffor < tol12) then
           write (message,'(3a)') ' toldff criterion is satisfied, but your forces are suspiciously low.', ch10,&
&           ' Check if the forces are 0 by symmetry: in that case you can not use the toldff convergence criterion!'
           MSG_WARNING(message)
           if (maxfor < tol16 .and. res2 > tol9) tolrff_ok=0
         end if
       else
         toldff_ok=0
         use_dpfft = diffor < tol6
       end if

       if(toldff_ok==2 .and. .not.noquit)then
         write(message, '(a,a,i5,a,a,a,es11.3,a,es11.3)' ) ch10, &
&         ' At SCF step',istep,', forces are converged : ',ch10,&
&         '  for the second time, max diff in force=',diffor,' < toldff=',toldff
         call wrtout([std_out, ab_out], message)
         quit=1
       end if
     end if

     ! Here treat the tolrff criterion: if maximum change of fcart is less than
     ! input tolrff times fcart itself twice consecutively, stop steps (exit loop here)
     if (ttolrff==1) then
       if (istep==1) then
         tolrff_ok=0
         ! 27/7/2009: added test for absolute value of maxfor, otherwise if it is 0 this never exits the scf loop.
       else if (diffor < tolrff*maxfor .or. (maxfor < tol16 .and. diffor < tol16)) then
         tolrff_ok=tolrff_ok+1
           ! Thu Mar 12 19:01:40 MG: added additional check on res2 to make sure the SCF cycle is close to convergence.
           ! Needed for structural relaxations otherwise the stress tensor is wrong and the relax algo makes wrong moves.
         if (maxfor < tol16 .and. res2 > tol9) tolrff_ok=0
       else
         tolrff_ok=0
         use_dpfft = diffor < tolrff * maxfor * five
       end if
       if(tolrff_ok==2 .and. (.not.noquit))then
         write(message, '(a,a,i5,a,a,a,es11.3,a,es11.3,a)' ) ch10, &
         ' At SCF step',istep,', forces are sufficiently converged : ',ch10,&
         '  for the second time, max diff in force=',diffor,&
         ' is less than < tolrff=',tolrff, ' times max force'
         call wrtout([std_out, ab_out], message)
         quit=1
       end if
     end if

     ! Here treat the toldfe criterion: if the change of energy is less than
     ! input toldfe twice consecutively, stop steps (exit loop here)
     if (ttoldfe==1) then
       if (istep==1) then
         toldfe_ok=0
       else if (abs(deltae)<toldfe) then
         toldfe_ok=toldfe_ok+1
       else
         toldfe_ok=0
         use_dpfft = abs(deltae) < tol8
       end if
       if(toldfe_ok==2 .and. (.not.noquit))then
         if(usefock==0 .or. nnsclohf<2)then
           write(message, '(a,a,i5,a,a,a,es11.3,a,es11.3)' ) ch10, &
            ' At SCF step',istep,', etot is converged : ',ch10,&
            '  for the second time, diff in etot=',abs(deltae),' < toldfe=',toldfe
         else
           write(message, '(a,i3,a,i3,a,a,a,es11.3,a,es11.3)' ) &
            ' Outer loop step',istep_fock_outer,' - inner step',istep_mix,' - frozen Fock etot converged : ',ch10,&
            '  for the second time, diff in etot=',abs(deltae),' < toldfe=',toldfe
         endif
         call wrtout([std_out, ab_out], message)
         quit=1
       end if
       if(usefock==1 .and. nnsclohf>1)then
         if(istep_mix==1 .and. (.not.noquit))then
!          The change due to the update of the Fock operator is sufficiently small. No need to meet it a second times.
           if (abs(deltae)<toldfe) then
             write(message, '(a,i3,a,i3,a,a,a,es11.3,a,es11.3)' ) &
             ' Outer loop step',istep_fock_outer,' - inner step',istep_mix,' - etot converged : ',ch10,&
             '  update of Fock operator yields diff in etot=',abs(deltae),' < toldfe=',toldfe
             call wrtout([std_out, ab_out], message)
             fock%fock_common%fock_converged=.true.
             quit=1
           endif
         endif
         if(istep_mix==nnsclohf .and. quit==0)then
           write(message, '(a,i3,a,i3,a,a,a,es11.3,a,es11.3)' ) &
           ' Outer loop step',istep_fock_outer,' - inner step',istep_mix,' - frozen Fock etot NOT converged : ',ch10,&
           '  diff in etot=',abs(deltae),' > toldfe=',toldfe
           call wrtout([std_out, ab_out], message)
         endif
       endif

!      Here treat the vdw_df_threshold criterion for non-SCF vdW-DF
!      calculations: If input vdw_df_threshold is lesss than toldfe
!      then the vdW-DF is triggered once selfconsistency criteria is
!      reached for the first time.
!      write(message,'(1x,a,e10.3,1x,a,e10.3,1x,l1,a)') &
!      &      '[vdW-DF][DEBUG] deltae=',deltae,'vdw_df_threshold=',vdw_df_threshold, &
!      &      (abs(deltae)<toldfe),ch10
!      call wrtout(std_out,message,'COLL')
#if defined DEV_YP_VDWXC
       ivdw = 0
       if ( toldfe > vdw_df_threshold ) then
         ivdw = ivdw + 1
       end if
       call xc_vdw_trigger((toldfe_ok==1 .and. toldfe>vdw_df_threshold))
       if ( ivdw == 2) then
         quit=1
       end if
#endif
     end if

     ! Here treat the tolvrs criterion: if density/potential residual (squared)
     ! is less than input tolvrs, stop steps (exit loop here)
     if (ttolvrs==1 .and. .not. noquit) then
       if (res2 < tolvrs) then
         if (optres==0) then
           write(message, '(a,a,i5,a,1p,e10.2,a,e10.2,a)' ) ch10,&
            ' At SCF step',istep,'       vres2   =',res2,' < tolvrs=',tolvrs,' =>converged.'
         else
           write(message, '(a,a,i5,a,1p,e10.2,a,e10.2,a)' ) ch10,&
            ' At SCF step',istep,'       nres2   =',res2,' < tolvrs=',tolvrs,' =>converged.'
         end if
         call wrtout([std_out, ab_out], message)
         quit=1
       else
         use_dpfft = res2 < tol5
       end if
     end if

     if (quit==1.and.noquit) then
       write(message, '(a,a,a)' ) ch10, &
        ' SCF cycle will continue as it is in an initialization stage',' (occ. matrix was kept constant)...'
       call wrtout([std_out, ab_out], message)
     end if

   end if

   ! Activate FFT in double-precision.
   if (use_dpfft) then
     if (fftcore_mixprec == 1) call wrtout(std_out, " Approaching convergence. Activating FFT in double-precision")
     ii = fftcore_set_mixprec(0)
   end if

 case (3)
   ! If wavefunction convergence was not reached (for nstep>0) print a warning and return conv_retcode
   conv_retcode = 0
   if(nstep>0) then
     if (.not. converged()) then
       conv_retcode = 1

       if(iscf>=1 .or. iscf==-3 .or. iscf == 0)then
         write(message, '(a,a,a,a,i5,a)' ) ch10,&
         ' scprqt:  WARNING -',ch10,&
         '  nstep=',nstep,' was not enough SCF cycles to converge;'

         write(std_out,'(6a,i0,3a)')ch10,&
         "--- !ScfConvergenceWarning",ch10,&
         "message: |",ch10,&
         '    nstep ',nstep,' was not enough SCF cycles to converge.',ch10,&
         "..."
           !MSG_WARNING_CLASS(message, "ScfConvergenceWarning")
       else
         write(message, '(a,a,a,a,i5,a)' ) ch10,&
         ' scprqt:  WARNING -',ch10,&
         '  nstep=',nstep,' was not enough non-SCF iterations to converge;'

         write(std_out,'(8a)')ch10,&
         "--- !NscfConvergenceWarning",ch10,&
         "message: |",ch10,TRIM(indent(message)),ch10,&
         "..."
           !MSG_WARNING_CLASS(message, "NScfConvergenceWarning")
       end if
       call wrtout([std_out, ab_out], message)

       if (ttolwfr==1) then
         if (dtset%usewvl == 0) then
           write(message, '(a,es11.3,a,es11.3,a)' ) &
           '  maximum residual=',residm,' exceeds tolwfr=',tolwfr,ch10

           write(message2, '(a,es11.3,2a)' ) &
           '  maximum residual each band. tolwfr= ',tolwfr,ch10,&
           '  iband, isppol, individual band residuals (max over all k-points):'
           call wrtout(std_out, message2,'COLL')
           do isppol = 1, dtset%nsppol
             do iband = 1, dtset%mband
               write(message3, '(2i6, es11.3)') iband, isppol, residm_band(iband,isppol)
               call wrtout(std_out,message3,'COLL')
             end do
           end do

         else
           write(message, '(a,es11.3,a,es11.3,a)' ) &
           '  maximum grdnorm=',residm,' exceeds tolwfr=',tolwfr,ch10
         end if

       else if (ttoldff==1) then
         write(message, '(a,es11.3,a,es11.3,a)' ) &
         '  maximum force difference=',diffor,' exceeds toldff=',toldff,ch10

       else if (ttolrff==1) then
         write(message, '(a,es11.3,a,es11.3,a)' ) &
         '  maximum force difference=',diffor,' exceeds tolrff*maxfor=',tolrff*maxfor,ch10

       else if (ttoldfe==1) then
         write(message, '(a,es11.3,a,es11.3,a)' ) &
         '  maximum energy difference=',abs(deltae),' exceeds toldfe=',toldfe,ch10

       else if(ttolvrs==1)then
         if (optres==0) then
           write(message, '(a,es11.3,a,es11.3,a)' ) &
           '  potential residual=',res2,' exceeds tolvrs=',tolvrs,ch10
         else
           write(message, '(a,es11.3,a,es11.3,a)' ) &
           '  density residual=',res2,' exceeds tolvrs=',tolvrs,ch10
         end if
       end if
       call wrtout([std_out, ab_out], message)

       if (prtxml == 1) then
         write(ab_xml_out, "(A)", advance = "NO") '      <status cvState="Failed"'
       end if

     else
       ! Convergence is OK
       if (prtxml == 1) then
         write(ab_xml_out, "(A)", advance = "NO") '      <status cvState="Ok"'
       end if
     end if ! test for convergence reached or not

     if (prtxml == 1) then
       if (ttoldfe == 1) then
         write(ab_xml_out, "(A)") ' stop-criterion="toldfe" />'
       else if (ttoldff == 1) then
         write(ab_xml_out, "(A)") ' stop-criterion="toldff" />'
       else if (ttolrff == 1) then
         write(ab_xml_out, "(A)") ' stop-criterion="tolrff" />'
       else if (ttolvrs == 1) then
         write(ab_xml_out, "(A)") ' stop-criterion="tolvrs" />'
       else if (ttolwfr == 1) then
         write(ab_xml_out, "(A)") ' stop-criterion="tolwfr" />'
       else
         write(ab_xml_out, "(A)") ' />'
       end if
     end if

     ! If enabled, output a YAML document with the ETOT iterations
     call neat_finish_etot(etot_yaml_doc, ab_out)
   end if ! nstep == 0 : no output

 case default
   write(message, '(a,i0,a)' )' choice = ',choice,' is not an allowed value.'
   MSG_BUG(message)
 end select

 ! Additional stuff for the Fock+SCF cycle
 if (present(fock)) then
   if (associated(fock)) then
     fock%fock_common%scf_converged=(quit==1)
     ! At present, the decision that the Fock loop is converged is not taken here
     if (.not.fock%fock_common%fock_converged)quit=0
   end if
 end if

 ! Additional stuff for the two-component DFT SCF cycle (electrons+positron)
 if (present(electronpositron)) then
   if (associated(electronpositron)) then
     electronpositron%scf_converged=(quit==1)
     if (dtset%positron<0) then
       diff_e=abs(etotal-electronpositron%etotal_prev)
       diff_f=abs(maxfor-electronpositron%maxfor_prev)
     end if
     if (choice==1) then
       ttoldff=0;ttoldfe=0
       if(abs(dtset%postoldff)>tiny(0.0_dp))ttoldff=1
       if(abs(dtset%postoldfe)>tiny(0.0_dp))ttoldfe=1
       if (dtset%positron<0.and.ttoldff+ttoldfe/=1.and.iscf>0) then
         MSG_ERROR('one and only one of toldff or toldfe must differ from zero !')
       end if
     end if
     if (choice==2) then
       if (dtset%positron<0.and.istep<=nstep) then
         if (electronpositron%scf_converged) then
           if (electronpositron%istep/=electronpositron%nstep) then
             if ((.not.noquit).and.&
&             (diff_e<electronpositron%postoldfe.or.diff_f<electronpositron%postoldff).and.&
&             (mod(electronpositron%calctype,2)==0.or.(dtset%positron>-20.and.dtset%positron/=-2))) then
               if (diff_e<electronpositron%postoldfe) then
                 write(message, '(2a,i5,5a,es11.3,a,es11.3)' ) ch10, &
&                 ' At SCF step',istep,', the difference between',ch10,&
&                 ' etotal from electronic calculation and etotal from positronic calculation',ch10,&
&                 ' is converged :  diff(etot_el-etot_pos)=',diff_e,' < postoldfe=',electronpositron%postoldfe
               else
                 write(message, '(2a,i5,5a,es11.3,a,es11.3)' ) ch10, &
&                 ' At SCF step',istep,', the difference between',ch10,&
&                 ' max. force from electronic calculation and max. force from positronic calculation',ch10,&
&                 ' is converged :  diff(maxfor_el-maxfor_pos)=',diff_f,' < postoldff=',electronpositron%postoldff
               end if
               call wrtout([std_out, ab_out], message)
             else
               quit=0
             end if
           end if
         end if
       end if
     end if
     if (choice==3) then
       if (dtset%positron<0.and.nstep>0)then
         if (diff_e>=electronpositron%postoldfe.and.abs(dtset%postoldfe)>tiny(0.0_dp)) then
           write(message, '(4a,i5,5a,es11.3,a,es11.3)' ) ch10,&
&           ' scprqt:  WARNING -',ch10,&
&           '  posnstep=',dtset%posnstep,' was not enough SCF cycles to converge difference between',ch10,&
&           '  etotal from electronic calculation and etotal from positronic calculation;',ch10,&
&           '  diff=',diff_e,' exceeds postoldfe=',electronpositron%postoldfe
           call wrtout([std_out, ab_out], message)
         end if
         if (diff_f>=electronpositron%postoldff.and.abs(dtset%postoldff)>tiny(0.0_dp)) then
           write(message, '(4a,i5,5a,es11.3,a,es11.3)' ) ch10,&
&           ' scprqt:  WARNING -',ch10,&
&           '  posnstep=',dtset%posnstep,' was not enough SCF cycles to converge difference between',ch10,&
&           '  max. force from electronic calculation and max. force from positronic calculation;',ch10,&
&           '  diff=',diff_e,' exceeds postoldff=',electronpositron%postoldff
           call wrtout([std_out, ab_out], message)
         end if
       end if
     end if
   end if
 end if

 call flush_unit(ab_out)

 DBG_EXIT("COLL")

 contains

   logical function converged()

   ! LB-02/01/2017:
   ! This code avoids evaluation of undefined variables (which could happen in respfn, apparently)
   logical :: loc_conv
   loc_conv = .true.
   if (ttolwfr==1) then
     if (residm > tolwfr) loc_conv=.false.
   end if
   if (ttoldff==1) then
     if (diffor > toldff) loc_conv=.false.
   end if
   if (ttolrff==1) then
     if (diffor > tolrff*maxfor .and. maxfor > tol16) loc_conv=.false.
   end if
   if (ttoldfe==1) then
     if (abs(deltae) > toldfe) loc_conv=.false.
   end if
   if (ttolvrs==1) then
     if (res2  > tolvrs) loc_conv=.false.
   end if
   converged = loc_conv

 end function converged

end subroutine scprqt
!!***

!!****f* ABINIT/setup1
!! NAME
!! setup1
!!
!! FUNCTION
!! Call near top of main routine to handle setup of various arrays,
!! filenames, checking of input data, etc.
!!
!! INPUTS
!!  acell(3)=length scales (bohr)
!!  ecut_eff=effective energy cutoff (hartree) for planewave basis sphere
!!  ecutc_eff=- PAW only - effective energy cutoff (hartree) for the coarse grid
!!  natom=number of atoms
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ngfftc(18)=contain all needed information about 3D FFT for the coarse grid
!!  nkpt=number of k points
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms
!!  response=0 if called by gstate, =1 if called by respfn
!!  rprim(3,3)=dimensionless real space primitive translations
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!  bantot=total number of bands at all k points
!!  gmet(3,3)=metric for reciprocal space inner products (bohr^-2)
!!  gprimd(3,3)=dimens. primitive translations for reciprocal space (bohr**-1)
!!  gsqcut_eff=Fourier cutoff on G^2 for "large sphere" of radius double
!!  gsqcutc_eff=(PAW) Fourier cutoff on G^2 for "large sphere" of radius double for the coarse FFT grid
!!   that of the basis sphere--appropriate for charge density rho(G),
!!   Hartree potential, and pseudopotentials, corresponding to ecut_eff
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol=unit cell volume (bohr^3)
!!
!! NOTES
!! SHOULD BE CLEANED !
!!
!! PARENTS
!!      gstate,nonlinear,respfn
!!
!! CHILDREN
!!      getcut,metric,mkrdim,wrtout
!!
!! SOURCE

subroutine setup1(acell,bantot,dtset,ecut_eff,ecutc_eff,gmet,&
&  gprimd,gsqcut_eff,gsqcutc_eff,ngfft,ngfftc,nkpt,nsppol,&
&  response,rmet,rprim,rprimd,ucvol,usepaw)

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: nkpt,nsppol
 integer,intent(in) :: response,usepaw
 integer,intent(out) :: bantot
 real(dp),intent(in) :: ecut_eff,ecutc_eff
 real(dp),intent(out) :: gsqcut_eff,gsqcutc_eff,ucvol
!arrays
 integer,intent(in) :: ngfft(18),ngfftc(18)
 real(dp),intent(in) :: acell(3),rprim(3,3)
 real(dp),intent(out) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),intent(out) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ikpt,isppol
 real(dp) :: boxcut,boxcutc
 character(len=500) :: message
!arrays
 real(dp) :: k0(3)

! ************************************************************************

!Compute bantot
 bantot=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     bantot=bantot+dtset%nband(ikpt+(isppol-1)*nkpt)
   end do
 end do

 if(dtset%nqpt>1.or.dtset%nqpt<0) then
   write(message,'(a,i0,5a)')&
   'nqpt =',dtset%nqpt,' is not allowed',ch10,&
   '(only 0 or 1 are allowed).',ch10,&
   'Action: correct your input file.'
   MSG_ERROR(message)
 end if

 ! Compute dimensional primitive translations rprimd
 call mkrdim(acell,rprim,rprimd)

 ! Obtain dimensional translations in reciprocal space gprimd,
 ! metrics and unit cell volume, from rprimd.
 ! Also output rprimd, gprimd and ucvol
 call metric(gmet,gprimd,ab_out,rmet,rprimd,ucvol)

 ! Get boxcut for given acell, gmet, ngfft, and ecut_eff
 ! (center at 000 for groundstate, center at q for respfn):
 ! boxcut=ratio of basis sphere diameter to fft box side
 k0(:)=0.0_dp
 if(response==1 .and. dtset%nqpt==1)then
   k0(:)=dtset%qptn(:)
   write(message, '(a)' )' setup1 : take into account q-point for computing boxcut.'
   call wrtout([std_out, ab_out], message)
 end if
 if (usepaw==1) then
   write(message,'(2a)') ch10,' Coarse grid specifications (used for wave-functions):'
   call wrtout([std_out, ab_out], message)
   call getcut(boxcutc,ecutc_eff,gmet,gsqcutc_eff,dtset%iboxcut,ab_out,k0,ngfftc)
   write(message,'(2a)') ch10,' Fine grid specifications (used for densities):'
   call wrtout([std_out, ab_out], message)
   call getcut(boxcut,ecut_eff,gmet,gsqcut_eff,dtset%iboxcut,ab_out,k0,ngfft)
 else
   call getcut(boxcut,ecut_eff,gmet,gsqcut_eff,dtset%iboxcut,ab_out,k0,ngfft)
   gsqcutc_eff=gsqcut_eff
 end if

 ! Check that boxcut>=2 if dtset%intxc=1; otherwise dtset%intxc must be set=0
 if (boxcut<2.0_dp.and.dtset%intxc==1) then
   write(message, '(a,es12.4,a,a,a,a,a)' )&
   'boxcut= ',boxcut,' is < 2.0  => intxc must be 0;',ch10,&
   'Need larger ngfft to use intxc=1.',ch10,&
   'Action: you could increase ngfft, or decrease ecut, or put intxcn=0.'
   MSG_ERROR(message)
 end if

end subroutine setup1
!!***

!!****f* ABINIT/prteigrs
!!
!! NAME
!! prteigrs
!!
!! FUNCTION
!! Print out eigenvalues band by band and k point by k point.
!! If option=1, do it in a standard way, for self-consistent calculations.
!! If option=2, print out residuals and eigenvalues, in a format
!! adapted for nonself-consistent calculations, within the loops.
!! If option=3, print out eigenvalues, in a format
!! adapted for nonself-consistent calculations, at the end of the job.
!! If option=4, print out derivatives of eigenvalues (same format as option==3, except header that is printed)
!! If option=5, print out Fan contribution to zero-point motion correction to eigenvalues (averaged)
!!                  (same format as option==3, except header that is printed)
!! If option=6, print out DDW contribution to zero-point motion correction to eigenvalues (averaged)
!!                  (same format as option==3, except header that is printed)
!! If option=7, print out Fan+DDW contribution to zero-point motion correction to eigenvalues (averaged)
!!                  (same format as option==3, except header that is printed)
!!
!! INPUTS
!!  eigen(mband*nkpt*nsppol)=eigenvalues (hartree)
!!   or, if option==4, diagonal of derivative of eigenvalues
!!   or, if option==5...7, zero-point motion correction to eigenvalues (averaged)
!!  enunit=choice parameter: 0=>output in hartree; 1=>output in eV;
!!   2=> output in both hartree and eV
!!  fermie=fermi energy (Hartree)
!!  fname_eig=filename for printing of the eigenenergies
!!  iout=unit number for formatted output file
!!  iscf=option for self-consistency
!!  kptns(3,nkpt)=k points in reduced coordinates
!!  kptopt=option for the generation of k points
!!  mband=maximum number of bands
!!  nband(nkpt)=number of bands at each k point
!!  nkpt=number of k points
!!  nnsclo_now=number of non-self-consistent loops for the current vtrial
!!    (often 1 for SCF calculation, =nstep for non-SCF calculations)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(maxval(nband(:))*nkpt*nsppol)=occupancies for each band and k point
!!  occopt=option for occupancies
!!  option= (see above)
!!  prteig=control print eigenenergies
!!  prtvol=control print volume and debugging
!!  resid(mband*nkpt*nsppol)=residuals (hartree**2)
!!  tolwfr=tolerance on band residual of wf, hartrees**2 (needed when option=2)
!!  vxcavg=average of vxc potential
!!  wtk(nkpt)=k-point weights
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      clnup1,dfpt_looppert,respfn,scprqt,vtorho
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine prteigrs(eigen,enunit,fermie,fname_eig,iout,iscf,kptns,kptopt,mband,nband,&
&  nkpt,nnsclo_now,nsppol,occ,occopt,option,prteig,prtvol,resid,tolwfr,vxcavg,wtk)

 use m_io_tools,  only : open_file

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: enunit,iout,iscf,kptopt,mband,nkpt,nnsclo_now,nsppol
 integer,intent(in) :: occopt,option,prteig,prtvol
 real(dp),intent(in) :: fermie,tolwfr,vxcavg
 character(len=*),intent(in) :: fname_eig
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),kptns(3,nkpt)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),resid(mband*nkpt*nsppol)
 real(dp),intent(in) :: wtk(nkpt)

!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt_max=50
 integer :: band_index,iband,ienunit,ii,ikpt,isppol,nband_index,nband_k,nkpt_eff,tmagnet,tmetal,temp_unit
 real(dp) :: convrt,magnet,residk,rhodn,rhoup
 character(len=2) :: ibnd_fmt,ikpt_fmt
 character(len=7) :: strunit1,strunit2
 character(len=39) :: kind_of_output
 character(len=500) :: msg

! *************************************************************************

 if (enunit<0.or.enunit>2) then
   write(msg, '(a,i0)' )' enunit must be 0, 1 or 2. Argument was ',enunit
   MSG_BUG(msg)
 end if

 if (prteig > 0) then
   write(msg, '(2a)' ) ' prteigrs : about to open file ',TRIM(fname_eig)
   call wrtout(iout,msg,'COLL')
   if (open_file(fname_eig, msg, newunit=temp_unit, status='unknown', form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if
   rewind(temp_unit) ! always rewind disk file and print latest eigenvalues
 end if

 kind_of_output=              ' Eigenvalues                          '
 if(option==4) kind_of_output=' Expectation of eigenvalue derivatives'
 if(option==5) kind_of_output=' Fan corrections to eigenvalues at T=0'
 if(option==6) kind_of_output=' DDW corrections to eigenvalues at T=0'
 if(option==7) kind_of_output=' Fan+DDW corrs   to eigenvalues at T=0'

 nkpt_eff=nkpt

!write(msg,'(a,5i5)')' prtvol,iscf,kptopt,nkpt_eff,nkpt_max ',prtvol,iscf,kptopt,nkpt_eff,nkpt_max
!call wrtout(iout,msg,'COLL')

 if( (prtvol==0.or.prtvol==1) .and. (iscf/=-2 .or. kptopt>0) .and. nkpt_eff>nkpt_max)nkpt_eff=nkpt_max
 if( (prtvol==0.or.prtvol==1) .and. (iscf/=-2 .or. kptopt>0) .and. nkpt_eff>1 .and. iout==ab_out)nkpt_eff=1

 if(option==1 .or. (option>=3 .and. option<=7))then

   do ienunit=0,1

     if (enunit==1 .and. ienunit==0)cycle
     if (enunit==0 .and. ienunit==1)cycle
     ! Print eigenvalues in hartree for enunit=0 or 2
     ! The definition of two different strings is quite ridiculous. Historical reasons ...

     if (ienunit==0)then
       convrt=one
       strunit1='hartree'
       strunit2='hartree'
     end if
     if (ienunit==1)then
       convrt=Ha_eV
       strunit1='   eV  '
       strunit2='eV     '
     end if

     band_index=0

     if(ienunit==0)then  ! XG20140730 I do not know why this is only done when ienunit==0
       tmetal=0
       if(option==1 .and. occopt>=3 .and. occopt<=8)tmetal=1
       tmagnet=0
       if(tmetal==1 .and. nsppol==2)then
         tmagnet=1
         rhoup = 0._dp
         rhodn = 0._dp
         nband_index = 1
         do isppol=1,nsppol
           do ikpt=1,nkpt
             nband_k=nband(ikpt+(isppol-1)*nkpt)
             do iband=1,nband_k
               if(isppol==1) rhoup = rhoup + wtk(ikpt)*occ(nband_index)
               if(isppol==2) rhodn = rhodn + wtk(ikpt)*occ(nband_index)
               nband_index = nband_index + 1
             end do
           end do
         end do
         magnet = abs(rhoup - rhodn)
       end if
     end if

     if(iscf>=0 .and. (ienunit==0 .or. option==1))then
       write(msg, '(3a,f10.5,3a,f10.5)' ) &
        ' Fermi (or HOMO) energy (',trim(strunit2),') =',convrt*fermie,'   Average Vxc (',trim(strunit2),')=',convrt*vxcavg
       call wrtout(iout,msg,'COLL')
       if (prteig > 0) call wrtout(temp_unit,msg,'COLL')
     end if

     ! if( (iscf>=0 .or. iscf==-3) .and. ienunit==0)then     ! This is the most correct
     if(iscf>=0 .and. ienunit==0)then ! For historical reasons
       if(tmagnet==1)then
         write(msg, '(a,es16.8,a,a,es16.8,a,es16.8)' )&
&         ' Magnetization (Bohr magneton)=',magnet,ch10,&
&         ' Total spin up =',rhoup,'   Total spin down =',rhodn
         call wrtout(iout,msg,'COLL')
         if (prteig > 0) call wrtout(temp_unit,msg,'COLL')
       end if
     end if

     ! Loop over spins (suppress spin data if nsppol not 2)
     do isppol=1,nsppol

       ikpt_fmt="i4" ; if(nkpt>=10000)ikpt_fmt="i6" ; if(nkpt>=1000000)ikpt_fmt="i9"
       if (nsppol==2.and.isppol==1) then
         write(msg, '(4a,'//ikpt_fmt//',2x,a)' ) &
&         trim(kind_of_output),' (',strunit1,') for nkpt=',nkpt,'k points, SPIN UP:'
       else if (nsppol==2.and.isppol==2) then
         write(msg, '(4a,'//ikpt_fmt//',2x,a)' ) &
&         trim(kind_of_output),' (',strunit1,') for nkpt=',nkpt,'k points, SPIN DOWN:'
       else
         write(msg, '(4a,'//ikpt_fmt//',2x,a)' ) &
&         trim(kind_of_output),' (',strunit1,') for nkpt=',nkpt,'k points:'
       end if
       call wrtout(iout,msg,'COLL')
       if (prteig > 0) call wrtout(temp_unit,msg,'COLL')

       if(ienunit==0)then
         if(option>=4 .and. option<=7)then
           msg = '  (in case of degenerate eigenvalues, averaged derivative)'
           call wrtout(iout,msg,'COLL')
           if (prteig > 0) call wrtout(temp_unit,msg,'COLL')
         end if
       end if

       do ikpt=1,nkpt
         nband_k=nband(ikpt+(isppol-1)*nkpt)
         ikpt_fmt="i4" ; if(nkpt>=10000)ikpt_fmt="i6" ; if(nkpt>=1000000)ikpt_fmt="i9"
         ibnd_fmt="i3" ; if(nband_k>=1000)ibnd_fmt="i6" ; if(nband_k>=1000000)ibnd_fmt="i9"
         if(ikpt<=nkpt_eff)then
           write(msg, '(a,'//ikpt_fmt//',a,'//ibnd_fmt//',a,f9.5,a,3f8.4,a)' ) &
&           ' kpt#',ikpt,', nband=',nband_k,', wtk=',wtk(ikpt)+tol10,', kpt=',&
&           kptns(1:3,ikpt)+tol10,' (reduced coord)'
           call wrtout(iout,msg,'COLL')
           if (prteig > 0) call wrtout(temp_unit,msg,'COLL')
           do ii=0,(nband_k-1)/8
!            write(msg, '(8f15.10)' ) (convrt*eigen(iband+band_index),&
             write(msg, '(8(f10.5,1x))' ) (convrt*eigen(iband+band_index),&
&             iband=1+ii*8,min(nband_k,8+ii*8))
             call wrtout(iout,msg,'COLL')
             if (prteig > 0) call wrtout(temp_unit,msg,'COLL')
           end do
           if(ienunit==0 .and. option==1 .and. occopt>=3 .and. occopt<=8)then
             write(msg, '(5x,a,'//ikpt_fmt//')' )  ' occupation numbers for kpt#',ikpt
             call wrtout(iout,msg,'COLL')
             do ii=0,(nband_k-1)/8
               write(msg, '(8(f10.5,1x))' ) (occ(iband+band_index),iband=1+ii*8,min(nband_k,8+ii*8))
               call wrtout(iout,msg,'COLL')
             end do
           end if

         else
           if(ikpt==nkpt_eff+1)then
             write(msg, '(a,a)' )' prteigrs : prtvol=0 or 1, do not print more k-points.',ch10
             call wrtout(iout,msg,'COLL')
           end if
           if (prteig > 0) then
             write(msg, '(a,'//ikpt_fmt//',a,'//ibnd_fmt//',a,f9.5,a,3f8.4,a)' ) &
&             ' kpt#',ikpt,', nband=',nband_k,', wtk=',wtk(ikpt)+tol10,', kpt=',&
&             kptns(1:3,ikpt)+tol10,' (reduced coord)'
             call wrtout(temp_unit,msg,'COLL')
             do ii=0,(nband_k-1)/8
               write(msg, '(8(f10.5,1x))' ) (convrt*eigen(iband+band_index),iband=1+ii*8,min(nband_k,8+ii*8))
               call wrtout(temp_unit,msg,'COLL')
             end do
           end if
         end if
         band_index=band_index+nband_k
       end do ! do ikpt=1,nkpt
     end do ! do isppol=1,nsppol

   end do ! End loop over Hartree or eV

 else if(option==2)then

   band_index=0
   do isppol=1,nsppol

     if(nsppol==2)then
       if(isppol==1)write(msg, '(2a)' ) ch10,' SPIN UP channel '
       if(isppol==2)write(msg, '(2a)' ) ch10,' SPIN DOWN channel '
       call wrtout(iout,msg,'COLL')
       if(prteig>0) call wrtout(temp_unit,msg,'COLL')
     end if

     do ikpt=1,nkpt
       nband_k=nband(ikpt+(isppol-1)*nkpt)
       ikpt_fmt="i5" ; if(nkpt>=10000)ikpt_fmt="i7" ; if(nkpt>=1000000)ikpt_fmt="i9"

       if(ikpt<=nkpt_eff)then
         write(msg, '(1x,a,'//ikpt_fmt//',a,f9.5,2f9.5,a)' ) &
&         'Non-SCF case, kpt',ikpt,' (',(kptns(ii,ikpt),ii=1,3),'), residuals and eigenvalues='
         call wrtout(iout,msg,'COLL')
         if (prteig > 0) then
           write(msg, '(1x,a,'//ikpt_fmt//',a,f9.5,2f9.5,a)' ) &
&           'Non-SCF case, kpt',ikpt,' eig(',(kptns(ii,ikpt),ii=1,3),') '
           call wrtout(temp_unit,msg,'COLL')
         end if
         do ii=0,(nband_k-1)/8
           write(msg, '(1p,8e10.2)' )(resid(iband+band_index),iband=1+8*ii,min(8+8*ii,nband_k))
           call wrtout(iout,msg,'COLL')
         end do
         do ii=0,(nband_k-1)/6
           write(msg, '(1p,6e12.4)' )(eigen(iband+band_index),iband=1+6*ii,min(6+6*ii,nband_k))
           call wrtout(iout,msg,'COLL')
           if (prteig > 0) call wrtout(temp_unit,msg,'COLL')
         end do
       else
         if(ikpt==nkpt_eff+1)then
           write(msg, '(a,a)' )' prteigrs : prtvol=0 or 1, do not print more k-points.',ch10
           call wrtout(iout,msg,'COLL')
         end if
         if (prteig > 0) then
           write(msg, '(1x,a,i5,a,f9.5,2f9.5,a)' )'Non-SCF kpt',ikpt,' eig(',(kptns(ii,ikpt),ii=1,3),') '
           call wrtout(temp_unit,msg,'COLL')
           do ii=0,(nband_k-1)/6
             write(msg, '(1p,6e12.4)' )(eigen(iband+band_index),iband=1+6*ii,min(6+6*ii,nband_k))
             call wrtout(temp_unit,msg,'COLL')
           end do
         end if
       end if

       ! MG: I don't understand why we should include the buffer in the output.
       ! It's already difficult to make the tests pass for the residuals without the buffer if nband >> nbocc
       residk=maxval(resid(band_index+1:band_index+nband_k))
       if (residk>tolwfr) then
         write(msg, '(1x,a,2i5,a,1p,e13.5)' ) &
&         ' prteigrs : nnsclo,ikpt=',nnsclo_now,ikpt,' max resid (incl. the buffer)=',residk
         call wrtout(iout,msg,'COLL')
       end if

       band_index=band_index+nband_k
     end do
   end do
   call wrtout(iout," ",'COLL')

 else
   write(msg, '(a,i0,a)' )' option = ',option,', is not an allowed value.'
   MSG_BUG(msg)
 end if

 if (prteig > 0) close (temp_unit)

end subroutine prteigrs
!!***

!!****f* ABINIT/prtene
!!
!! NAME
!! prtene
!!
!! FUNCTION
!! Print components of total energy in nice format
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryphase
!!   | kptopt
!!   | occopt
!!   | positron=option for electron-positron calculation
!!   | tphysel="physical" electronic temperature with FD occupations
!!   | tsmear=smearing energy or temperature (if metal)
!!  energies <type(energies_type)>=values of parts of total energy
!!  iout=unit number to which output is written
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      gstate,scfcv
!!
!! CHILDREN
!!      energies_eval_eint,wrtout
!!
!! SOURCE

subroutine prtene(dtset,energies,iout,usepaw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,usepaw
 type(dataset_type),intent(in) :: dtset
 type(energies_type),intent(in) :: energies

!Local variables-------------------------------
!scalars
 integer :: ipositron,mu,optdc
 logical :: directE_avail,testdmft
 real(dp) :: eent,enevalue,etotal,etotaldc,exc_semilocal
 ! Do not modify the length of these strings
 character(len=14) :: eneName
 character(len=500) :: msg
 type(pair_list) :: e_components, e_components_dc
 type(stream_string) :: stream
!arrays
 character(len=10) :: EPName(1:2)=(/"Positronic","Electronic"/)

! *************************************************************************

 directE_avail=(usepaw==0.or.dtset%pawspnorb==0.or.dtset%pawcpxocc==2.or.dtset%kptopt==1.or.dtset%kptopt==2)

!============= Evaluate some parts of the energy ===========

 optdc=-1;ipositron=merge(0,2,dtset%positron==0)
 if (abs(energies%e_ewald)<1.e-15_dp.and.abs(energies%e_hartree)<1.e-15_dp) ipositron=1
 call energies_eval_eint(energies,dtset,usepaw,optdc,etotal,etotaldc)

!Here, treat the case of metals
!In re-smeared case the free energy is defined with tphysel
 if(dtset%occopt>=3 .and. dtset%occopt<=8)then
   if (abs(dtset%tphysel) < tol10) then
     eent=-dtset%tsmear * energies%entropy
   else
     eent=-dtset%tphysel * energies%entropy
   end if
 else
   eent=zero
 end if
! If DMFT is used and DMFT Entropy is not computed, then do not print
! non interacting entropy
 testdmft=(dtset%dmftcheck>=0.and.dtset%usedmft>=1.and.(sum(dtset%upawu(:,1))>=tol8.or.  &
& sum(dtset%jpawu(:,1))>tol8).and.dtset%dmft_entropy==0)
 if(testdmft) eent=zero

 etotal   = etotal   + eent
 etotaldc = etotaldc + eent

 write(msg,'(a,80a)') ch10,('-',mu=1,80)
 call wrtout(iout,msg,'COLL')

!============= Printing of Etotal by direct scheme ===========

 if (dtset%icoulomb == 1) then
   write(eneName, "(A)") "Ion-ion energy"
 else
   write(eneName, "(A)") "Ewald energy"
 end if
 enevalue = energies%e_ewald


 if (optdc==0.or.optdc==2) then

   if (directE_avail) then
     write(msg, '(2a)' ) ' Components of total free energy (in Hartree) :',ch10
     call wrtout(iout,msg,'COLL')
     call e_components%set('comment', s='Components of total free energy in Hartree')
     write(msg, '(a,es21.14)' ) '    Kinetic energy  = ',energies%e_kinetic
     call wrtout(iout,msg,'COLL')
     call e_components%set('kinetic', r=energies%e_kinetic)
     if (ipositron/=1) then
       exc_semilocal=energies%e_xc+energies%e_hybcomp_E0-energies%e_hybcomp_v0+energies%e_hybcomp_v
!XG20181025 This should NOT be a part of the semilocal XC energy, but treated separately.
!      At present, there is still a problem with the variational formulation for the Fock term with PAW.
!      So, for the time being, keep it inside.
       if(usepaw==1)exc_semilocal=exc_semilocal+energies%e_fock
       write(msg, '(3(a,es21.14,a),a,es21.14)' ) &
&       '    Hartree energy  = ',energies%e_hartree,ch10,&
&       '    XC energy       = ',exc_semilocal,ch10,&
&       '    '//eneName//'  = ',enevalue,ch10,&
&       '    PspCore energy  = ',energies%e_corepsp
       call wrtout(iout,msg,'COLL')
       call e_components%set('hartree', r=energies%e_hartree)
       call e_components%set('xc', r=exc_semilocal)
       call e_components%set(eneName, r=enevalue)
       call e_components%set('psp_core', r=energies%e_corepsp)
#if defined DEV_YP_VDWXC
       if ( (dtset%vdw_xc > 0) .and. (dtset%vdw_xc < 10) .and. (xc_vdw_status()) ) then
         write(msg, '(a,es21.14)' )'    vdW-DF energy   = ',energies%e_xc_vdw
         call wrtout(iout,msg,'COLL')
         call e_components%set('VdWaals_df', r=energies%e_xc_vdw)
       end if
#endif
     end if
     write(msg, '(a,es21.14)' ) '    Loc. psp. energy= ',energies%e_localpsp
     call wrtout(iout,msg,'COLL')
     call e_components%set('local_psp', r=energies%e_localpsp)
     if (usepaw==0) then
       if(abs(energies%e_fock0)<tol8)then
         write(msg, '(a,es21.14)' ) &
&         '    NL   psp  energy= ',energies%e_nlpsp_vfock
         call e_components%set('non_local_psp', r=energies%e_nlpsp_vfock)
       else
         write(msg, '(a,es21.14)' ) &
&         '    NL(psp+X) energy= ',energies%e_nlpsp_vfock-energies%e_fock0
         call e_components%set('non_local_psp+x', r=energies%e_nlpsp_vfock-energies%e_fock0)
       endif
       call wrtout(iout,msg,'COLL')
     else
       write(msg, '(a,es21.14)' ) &
&       '    Spherical terms = ',energies%e_paw
       call wrtout(iout,msg,'COLL')
       call e_components%set('spherical_terms', r=energies%e_paw)
       !XG20181025 Does not work (yet)...
       !if(abs(energies%e_nlpsp_vfock)>tol8)then
       !  write(msg, '(a,es21.14)' )'    Fock-type term  = ',energies%e_nlpsp_vfock
       !  call wrtout(iout,msg,'COLL')
       !  write(msg, '(a,es21.14)' ) '    -frozen Fock en.= ',-energies%e_fock0
       !  call wrtout(iout,msg,'COLL')
       !endif
     end if
     if ((dtset%vdw_xc>=5.and.dtset%vdw_xc<=7).and.ipositron/=1) then
       write(msg, '(a,es21.14)' ) '    Vd Waals DFT-D = ',energies%e_vdw_dftd
       call wrtout(iout,msg,'COLL')
       call e_components%set('VdWaals_dft_d', r=energies%e_vdw_dftd)
     end if
     if (dtset%nzchempot>=1) then
       write(msg, '(a,es21.14)' ) '    Chem. potential = ',energies%e_chempot
       call wrtout(iout,msg,'COLL')
       call e_components%set('chem_potential', r=energies%e_chempot)
     end if
     if(dtset%occopt>=3.and.dtset%occopt<=8.and.ipositron==0) then
       call e_components%set('internal', r=etotal-eent)
       if(.not.testdmft) then
         write(msg, '(a,es21.14,a,a,a,es21.14)' ) &
&         '    >>>>> Internal E= ',etotal-eent,ch10,ch10,&
&         '    -kT*entropy     = ',eent
         call wrtout(iout,msg,'COLL')
         call e_components%set('-kT*entropy', r=eent)
       else
         write(msg, '(a,es21.14,a)' ) &
&         '    >>>>> Internal E= ',etotal-eent,ch10
         call wrtout(iout,msg,'COLL')
       end if
     else if (ipositron/=0) then
       if (dtset%occopt>=3.and.dtset%occopt<=8) then
         write(msg, '(a,es21.14)' ) '    -kT*entropy     = ',eent
         call wrtout(iout,msg,'COLL')
         call e_components%set('-kT*entropy', r=eent)
       end if
       write(msg, '(3a,es21.14,a)' ) &
&       '    >>> ',EPName(ipositron),' E= ',etotal-energies%e0_electronpositron &
&       -energies%e_electronpositron,ch10
       call wrtout(iout,msg,'COLL')
       write(msg, '(3a,es21.14,2a,es21.14)' ) &
&       '    ',EPName(3-ipositron),' ener.= ',energies%e0_electronpositron,ch10,&
&       '    EP interaction E= '             ,energies%e_electronpositron
       call wrtout(iout,msg,'COLL')
       if(ipositron == 1) then
        call e_components%set('positronic', r=etotal- &
&                                         energies%e0_electronpositron-energies%e_electronpositron)
        call e_components%set('electronic', r=energies%e0_electronpositron)
       else
        call e_components%set('electronic', r=etotal- &
&                                         energies%e0_electronpositron-energies%e_electronpositron)
        call e_components%set('positronic', r=energies%e0_electronpositron)
       end if
       call e_components%set('electron_positiron_interaction', r=energies%e_electronpositron)
     end if
     if ((dtset%berryopt==4 .or.  dtset%berryopt==6 .or. dtset%berryopt==7 .or.  &
&     dtset%berryopt==14 .or. dtset%berryopt==16 .or. dtset%berryopt==17) .and.ipositron/=1) then
       write(msg, '(a,es21.14)' ) '    Electric energy = ',energies%e_elecfield
       call wrtout(iout,msg,'COLL')
       write(msg, '(a,es21.14)' ) '    Kohn-Sham energy= ',etotal-energies%e_elecfield
       call wrtout(iout,msg,'COLL')
       call e_components%set('electric', r=energies%e_elecfield)
       call e_components%set('kohn_sham', r=etotal-energies%e_elecfield)
     end if
     write(msg, '(a,es21.14)' ) '    >>>>>>>>> Etotal= ',etotal
     call wrtout(iout,msg,'COLL')
     call e_components%set('total_energy', r=etotal)

   else
     write(msg, '(9a)' ) &
&     ' COMMENT: ',ch10,&
&     '  "Direct" decomposition of total free energy cannot be printed out !!!',ch10,&
&     '  PAW contribution due to spin-orbit coupling cannot be evaluated',ch10,&
&     '  without the knowledge of imaginary part of Rhoij atomic occupancies',ch10,&
&     '  (computed only when pawcpxocc=2).'
     call wrtout(iout,msg,'COLL')
     call e_components%set('comment', s='"Direct" decomposition of total free energy cannot be printed out !!!'//ch10// &
&                       'PAW contribution due to spin-orbit coupling cannot be evaluated'//ch10// &
&                       'without the knowledge of imaginary part of Rhoij atomic occupancies'//ch10// &
&                       '(computed only when pawcpxocc=2).')
   end if
 end if
!============= Printing of Etotal by double-counting scheme ===========

 if (optdc>=1) then

   write(msg, '(4a,es21.14)' ) ch10,&
&   ' "Double-counting" decomposition of free energy:',ch10,&
&   '    Band energy     = ',energies%e_eigenvalues
   call wrtout(iout,msg,'COLL')
   call e_components_dc%set('comment', s='"Double-counting" decomposition of free energy')
   call e_components_dc%set('band_energy', r=energies%e_eigenvalues)
   if (ipositron/=1) then
     write(msg, '(2(a,es21.14,a),a,es21.14)' ) &
&     '    '//eneName//'  =',enevalue,ch10,&
&     '    PspCore energy  = ',energies%e_corepsp-energies%e_corepspdc,ch10,&
&     '    Dble-C XC-energy= ',-energies%e_hartree+energies%e_xc-energies%e_xcdc&
&     -energies%e_fock0+&
&     energies%e_hybcomp_E0-energies%e_hybcomp_v0
     call wrtout(iout,msg,'COLL')
     call e_components_dc%set(eneName, r=enevalue)
     call e_components_dc%set('psp_core', r=energies%e_corepsp-energies%e_corepspdc)
     call e_components_dc%set('xc_dc', r=-energies%e_hartree+energies%e_xc-energies%e_xcdc)
   end if
   if ((dtset%berryopt==4 .or.  dtset%berryopt==6 .or. dtset%berryopt==7 .or.  &
&   dtset%berryopt==14 .or. dtset%berryopt==16 .or. dtset%berryopt==17).and.ipositron/=1) then
     write(msg, '(a,es21.14)' ) '    Electric field  = ',energies%e_elecfield
     call wrtout(iout,msg,'COLL')
     call e_components_dc%set('electric_field', r=energies%e_elecfield)
   end if
   if (usepaw==1) then
     write(msg, '(a,es21.14)' ) '    Spherical terms = ',energies%e_pawdc
     call wrtout(iout,msg,'COLL')
     call e_components_dc%set('spherical_terms', r=energies%e_pawdc)
   end if
   if ((dtset%vdw_xc>=5.and.dtset%vdw_xc<=7).and.ipositron/=1) then
     write(msg, '(a,es21.14)' ) '    Vd Waals DFT-D = ',energies%e_vdw_dftd
     call wrtout(iout,msg,'COLL')
     call e_components_dc%set('VdWaals_dft_d', r=energies%e_vdw_dftd)
   end if
   if (dtset%nzchempot>=1) then
     write(msg, '(a,es21.14)' ) '    Chem. potential = ',energies%e_chempot
     call wrtout(iout,msg,'COLL')
     call e_components_dc%set('chem_potential', r=energies%e_chempot)
   end if
   if(dtset%occopt>=3.and.dtset%occopt<=8.and.ipositron==0) then
     if(.not.testdmft) then
       write(msg, '(a,es21.14,a,a,a,es21.14)' ) &
&       '    >>>>> Internal E= ',etotaldc-eent,ch10,ch10,&
&       '    -kT*entropy     = ',eent
       call wrtout(iout,msg,'COLL')
       call e_components_dc%set('internal', r=etotaldc-eent)
       call e_components_dc%set('-kT*entropy', r=eent)
     else
       write(msg, '(a,es21.14,a)' ) '    >>>>> Internal E= ',etotaldc-eent,ch10
       call wrtout(iout,msg,'COLL')
       call e_components_dc%set('internal', r=etotaldc-eent)
     end if
   else if (ipositron/=0) then
     if (dtset%occopt>=3 .and. dtset%occopt<=8) then
       write(msg, '(a,es21.14)' ) '    -kT*entropy     = ',eent
       call wrtout(iout,msg,'COLL')
       call e_components_dc%set('-kT*entropy', r=eent)
     end if
     write(msg, '(a,es21.14,4a,es21.14,a)' ) &
&     '    - EP dble-ct En.= ',-energies%edc_electronpositron,ch10,&
&     '    >>> ',EPName(ipositron),' E= ',etotaldc-energies%e0_electronpositron &
&     -energies%e_electronpositron,ch10
     call wrtout(iout,msg,'COLL')
     write(msg, '(3a,es21.14,2a,es21.14)' ) &
&     '    ',EPName(3-ipositron),' ener.= ',energies%e0_electronpositron,ch10,&
&     '    EP interaction E= '            ,energies%e_electronpositron
     call wrtout(iout,msg,'COLL')
     call e_components_dc%set('electron_positron_dc', r=-energies%edc_electronpositron)
     if(ipositron == 1) then
       call e_components_dc%set('positronic', r=etotaldc-energies%e0_electronpositron-energies%e_electronpositron)
       call e_components_dc%set('electronic', r=energies%e0_electronpositron)
     else
       call e_components_dc%set('electronic', r=etotaldc-energies%e0_electronpositron-energies%e_electronpositron)
       call e_components_dc%set('positronic', r=energies%e0_electronpositron)
     end if
     call e_components_dc%set('electron_positron_interaction', r=energies%e_electronpositron)
   end if
   write(msg, '(a,es21.14)' ) '    >>>> Etotal (DC)= ',etotaldc
   call wrtout(iout,msg,'COLL')
   call e_components_dc%set('total_energy_dc', r=etotaldc)

 end if

!======= Additional printing for compatibility  ==========

 if (usepaw==0.and.optdc==0) then
   write(msg, '(a,a,a,a,es21.14,a,es18.10)' ) ch10,&
&   ' Other information on the energy :',ch10,&
&   '    Total energy(eV)= ',etotal*Ha_eV,' ; Band energy (Ha)= ',energies%e_eigenvalues
   call wrtout(iout,msg,'COLL')
   call e_components%set('band_energy', r=energies%e_eigenvalues)
   call e_components%set('total_energy_eV', r=etotal*Ha_eV)
 end if

 if ((optdc==0.or.optdc==2).and.(.not.directE_avail)) then
   write(msg, '(a,a,es18.10)' ) ch10,' Band energy (Ha)= ',energies%e_eigenvalues
   call wrtout(iout,msg,'COLL')
   call e_components%set('band_energy', r=energies%e_eigenvalues)
 end if

 if (usepaw==1) then
   if ((optdc==0.or.optdc==2).and.(directE_avail)) then
     write(msg, '(a,a,es21.14)' ) ch10,'  >Total energy in eV           = ',etotal*Ha_eV
     call wrtout(iout,msg,'COLL')
     call e_components%set('total_energy_eV', r=etotal*Ha_eV)
   end if
   if (optdc>=1) then
     if (optdc==1) write(msg, '(a,a,es21.14)' ) ch10,&
&     '  >Total DC energy in eV        = ',etotaldc*Ha_eV
     if (optdc==2) write(msg, '(a,es21.14)' ) &
&     '  >Total DC energy in eV        = ',etotaldc*Ha_eV
     call wrtout(iout,msg,'COLL')
     call e_components_dc%set('total_energy_dc_eV', r=etotal*Ha_eV)
   end if
 end if

 if( dtset%icoulomb/=1.and.abs(dtset%charge)>tol8) then
   write(msg, '(6a)' ) &
&   ch10,' Calculation was performed for a charged system with PBC',&
&   ch10,' You may consider including the monopole correction to the total energy',&
&   ch10,' The correction is to be divided by the dielectric constant'
   call wrtout(iout,msg,'COLL')
   write(msg, '(a,es21.14)' ) '    Monopole correction (Ha)=',energies%e_monopole
   call wrtout(iout,msg,'COLL')
   write(msg, '(a,es21.14)' ) '    Monopole correction (eV)=',energies%e_monopole*Ha_eV
   call wrtout(iout,msg,'COLL')
   call e_components%set('monopole_correction', r=energies%e_monopole)
   call e_components%set('monopole_correction_eV', r=energies%e_monopole*Ha_eV)
 end if

 write(msg,'(a,80a)')('-',mu=1,80)
 call wrtout(iout,msg,'COLL')

 call wrtout(iout, ch10, 'COLL')

 ! Write components of total energies in a structured way
 call yaml_single_dict('EnergyTerms', '', e_components, 35, 500, width=20, stream=stream, real_fmt='(es21.14)')
 call stream%dump(iout)
 call e_components%free()

 if(e_components_dc%length() > 1) then
   call wrtout(iout, ch10, 'COLL')
   call yaml_single_dict('EnergyTermsDC', '', e_components_dc, 35, 500, width=20, stream=stream, real_fmt='(es21.14)')
   call stream%dump(iout)
   call e_components_dc%free()
 end if

end subroutine prtene
!!***

!!****f* ABINIT/get_dtsets_pspheads
!! NAME
!! get_dtsets_pspheads
!!
!! FUNCTION
!!  Parse input file, get list of pseudos for files file and build list of datasets
!!  pseudopotential headers, maxval of dimensions needed in outvars
!!
!! INPUTS
!!  path: Input Filename
!!  comm: MPI communicator
!!
!! OUTPUT
!!  lenstr= the length of the resulting string.
!!  ndtset= the number of declared datasets.
!!  string= contains on output the content of the file, ready for parsing.
!!  dtsets(0:ndtset): List of datasets
!!  dmatpuflag=flag controlling the use of an initial density matrix in PAW+U (max. value over datasets)
!!  mx<ab_dimensions>=datatype storing the maximal dimensions.
!!  pspheads(npsp)=<type pspheader_type>=all the important information from the
!!   pseudopotential file headers, as well as the psp file names
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_dtsets_pspheads(path, ndtset, lenstr, string, timopt, dtsets, pspheads, mx, dmatpuflag, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: lenstr, ndtset
 type(ab_dimensions),intent(out) :: mx
 character(len=strlen), intent(out) :: string
 character(len=*),intent(in) :: path
 integer,intent(in) :: comm
 integer,intent(out) :: timopt, dmatpuflag
!arrays
 type(dataset_type),allocatable,intent(out)  :: dtsets(:)
 type(pspheader_type),allocatable,intent(out):: pspheads(:)

!Local variables-------------------------------
!scalars
 integer :: ipsp,ios, me, ndtset_alloc, nprocs
 integer :: istatr,istatshft, papiopt, npsp, ii, idtset, msym, usepaw
 character(len=fnlen) :: filpsp
 character(len=500) :: msg
!arrays
 integer,allocatable :: mband_upper_(:)
 real(dp) :: ecut_tmp(3,2,10),tsec(2)
 real(dp),allocatable :: zionpsp(:)
 character(len=fnlen), allocatable :: pspfilnam_(:)

!************************************************************************

 ! Call the parser from the parser module.
 me = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 ! Read the file, stringify it and return the number of datasets.
 call parsefile(path, lenstr, ndtset, string, comm)

 ndtset_alloc = ndtset; if (ndtset == 0) ndtset_alloc=1
 ABI_DATATYPE_ALLOCATE(dtsets, (0:ndtset_alloc))

 timopt = 1; if (xmpi_paral==1) timopt = 0

 ! Continue to analyze the input string, get upper dimensions, and allocate the remaining arrays.
 call invars0(dtsets, istatr, istatshft, lenstr, msym, mx%natom, mx%nimage, mx%ntypat, &
              ndtset, ndtset_alloc, npsp, papiopt, timopt, string, comm)

 ! Enable PAPI timers
 call time_set_papiopt(papiopt)

 dtsets(:)%timopt = timopt
 dtsets(0)%timopt = 1
 if (xmpi_paral == 1) dtsets(0)%timopt = 0

 call timab(41,2,tsec)
 call timab(timopt,5,tsec)

 ! Finish to read the "file" file completely, as npsp is known,
 ! and also initialize pspheads, that contains the important information
 ! from the pseudopotential headers, as well as the psp filename

 call timab(42,1,tsec)

 usepaw = 0
 ABI_DATATYPE_ALLOCATE(pspheads,(npsp))
 if (npsp > 10) then
   MSG_BUG('ecut_tmp is not well defined.')
 end if
 ecut_tmp = -one

 pspheads(:)%usewvl = dtsets(1)%usewvl
 if (me == 0) then
    !if (.not. present(pspfilnam)) then
    ! Read the name of the psp file
    ABI_MALLOC(pspfilnam_,(npsp))
    do ipsp=1,npsp
      write(std_out,'(/,a)' )' Please give name of formatted atomic psp file'
      read (std_in, '(a)' , iostat=ios ) filpsp
      ! It might be that a file name is missing
      if (ios/=0) then
        write(msg, '(7a)' )&
        'There are not enough names of pseudopotentials',ch10,&
        'provided in the files file.',ch10,&
        'Action: check first the variable ntypat (and/or npsp) in the input file;',ch10,&
        'if they are correct, complete your files file.'
        MSG_ERROR(msg)
      end if
      pspfilnam_(ipsp) = trim(filpsp)
      write(std_out,'(a,i0,2a)' )' For atom type ',ipsp,', psp file is ',trim(filpsp)
    end do ! ipsp=1,npsp

    call inpspheads(pspfilnam_, npsp, pspheads, ecut_tmp)
    ABI_FREE(pspfilnam_)
    !else
    !   call inpspheads(pspfilnam, npsp, pspheads, ecut_tmp)
    !end if
    if (minval(abs(pspheads(1:npsp)%pspcod - 7)) == 0) usepaw=1
    if (minval(abs(pspheads(1:npsp)%pspcod - 17)) == 0) usepaw=1
 end if

 ! Communicate pspheads to all processors
 call pspheads_comm(npsp, pspheads, usepaw)

 ! If (all) pspcod are 7 then this is a PAW calculation. Initialize (default) the value of ratsph
 do idtset=0,ndtset_alloc
    dtsets(idtset)%usepaw = usepaw
    if (usepaw == 0) then
      dtsets(idtset)%ratsph(:)=two
    else
      ! Note that the following coding assumes that npsp=ntypat for PAW, which is true as of now (XG20101024).
      ! dtsets(idtset)%ratsph(1:npsp)=token%pspheads(1:npsp)%pawheader%rpaw
      do ipsp=1,npsp
        dtsets(idtset)%ratsph(ipsp) = pspheads(ipsp)%pawheader%rpaw
      end do
    endif
 end do

 !Take care of other dimensions, and part of the content of dtsets that is or might be needed early.
 !zion_max=maxval(pspheads(1:npsp)%zionpsp) ! This might not work properly with HP compiler

! zion_max=token%pspheads(1)%zionpsp
! do ii=1,npsp
!    zion_max=max(token%pspheads(ii)%zionpsp,zion_max)
! end do
 ABI_MALLOC(zionpsp,(npsp))
 do ii=1,npsp
  zionpsp(ii) = pspheads(ii)%zionpsp
 end do

 ABI_MALLOC(mband_upper_, (0:ndtset_alloc))

 ! Get MAX dimension over datasets
 call invars1m(dmatpuflag, dtsets, ab_out, lenstr, mband_upper_, mx,&
               msym, ndtset, ndtset_alloc, string, npsp, zionpsp, comm)

 ABI_FREE(zionpsp)
 call timab(42,2,tsec)
 call timab(43,3,tsec)

 ! Provide defaults for the variables that have not yet been initialized.
 call indefo(dtsets, ndtset_alloc, nprocs)
 call macroin(dtsets, ecut_tmp, lenstr, ndtset_alloc, string)

 ! Perform some global initialization, depending on the value of
 ! pseudopotentials, parallelism variables, or macro input variables

 ! If all the pseudopotentials have the same pspxc, override the default value for dtsets 1 to ndtset
 if (minval(abs((pspheads(1:npsp)%pspxc - pspheads(1)%pspxc)))==0) then
   dtsets(1:ndtset_alloc)%ixc = pspheads(1)%pspxc
 end if

 ! Call the main input routine.
 call invars2m(dtsets,ab_out,lenstr,mband_upper_,msym,ndtset,ndtset_alloc,npsp,pspheads,string, comm)

 call macroin2(dtsets, ndtset_alloc)

 mx%mband = dtsets(1)%mband
 do ii=1,ndtset_alloc
    mx%mband = max(dtsets(ii)%mband, mx%mband)
 end do

 call timab(43,2,tsec)

 ABI_FREE(mband_upper_)

end subroutine get_dtsets_pspheads
!!***


!!****f* ABINIT/ebands_from_file
!! NAME
!! ebands_from_file
!!
!! FUNCTION
!!  Build and ebands_t object from file. Supports Fortran and netcdf files
!!  provided they have a Abinit header and obviously GS eigenvalues
!!
!! INPUTS
!!  path: File name.
!!  comm: MPI communicator.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE


type(ebands_t) function ebands_from_file(path, comm) result(new)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path
 integer,intent(in) :: comm

!Local variables-------------------------------
!scalars
 integer :: fform, ncid
 type(hdr_type) :: hdr
!arrays
 real(dp),pointer :: gs_eigen(:,:,:)

! *************************************************************************

 ! NOTE: Assume file with header. Must use wfk_read_eigenvalues to handle Fortran WFK
 if (endswith(path, "_WFK") .or. endswith(path, "_WFK.nc")) then
   call hdr_read_from_fname(hdr, path, fform, comm)
   ABI_CHECK(fform /= 0, "fform == 0")
   call wfk_read_eigenvalues(path, gs_eigen, hdr, comm)
   new = ebands_from_hdr(hdr, maxval(hdr%nband), gs_eigen)

 else if (endswith(path, ".nc")) then
#ifdef HAVE_NETCDF
   NCF_CHECK(nctk_open_read(ncid, path, comm))
   call hdr_ncread(hdr, ncid, fform)
   ABI_MALLOC(gs_eigen, (hdr%mband, hdr%nkpt, hdr%nsppol))
   NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "eigenvalues"), gs_eigen))
   new = ebands_from_hdr(hdr, maxval(hdr%nband), gs_eigen)
   NCF_CHECK(nf90_close(ncid))
#endif
 else
   MSG_ERROR(sjoin("Don't know how to construct crystal structure from: ", path, ch10, "Supported extensions: _WFK or .nc"))
 end if

 ABI_FREE(gs_eigen)
 call hdr_free(hdr)

end function ebands_from_file
!!***


!!****f* ABINIT/crystal_from_file
!! NAME
!! crystal_from_file
!!
!! FUNCTION
!!  Build crystal_t object from netcdf file or Abinit input file with file extension in [".abi", ".in"]
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


type(crystal_t) function crystal_from_file(path, comm) result(new)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path
 integer,intent(in) :: comm

!Local variables-------------------------------
!scalars
 integer :: fform, timrev !, lenstr, ndtset, timopt, dmatpuflag
 !character(len=strlen) :: string
 type(hdr_type) :: hdr
 !type(ab_dimensions) :: mx
!arrays
 !type(dataset_type),allocatable :: dtsets(:)
 !type(pspheader_type),allocatable :: pspheads(:)

! *************************************************************************

 if (endswith(path, ".abi") .or. endswith(path, ".in")) then
   NOT_IMPLEMENTED_ERROR()

   ! TODO
   ! This routine prompts for the list of pseudos! One should get rid of the files file
   ! before activating this part.
   !call get_dtsets_pspheads(path, ndtset, lenstr, string, timopt, dtsets, pspheads, mx, dmatpuflag, comm)
   !call crystal_init(dtset%amu_orig(:,1), new, dtset%spgroup, dtset%natom, dtset%npsp, &
   !  psps%ntypat,dtset%nsym, rprimd, dtset%typat, xred, dtset%ziontypat, dtset%znucl,1, &
   !  dtset%nspden==2.and.dtset%nsppol==1, remove_inv, psps%title, &
   !  symrel=dtset%symrel, tnons=dtset%tnons, symafm=dtset%symafm)
   !ABI_FREE(dtsets)
   !ABI_FREE(pspheads)

 else
    ! Assume file header
    call hdr_read_from_fname(hdr, path, fform, comm)
    ABI_CHECK(fform /= 0, "fform == 0")
    timrev = 2 !; (if kpts_timrev_from_kptopt(hdr%kptopt) == 0) timrev = 1
    new = hdr_get_crystal(hdr, timrev)
    call hdr_free(hdr)
 end if

end function crystal_from_file
!!***

end module m_common
!!***
