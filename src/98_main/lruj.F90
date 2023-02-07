!!****p* ABINIT/lruj
!! NAME
!! lruj
!!
!! FUNCTION
!!  Linear Response U and J:
!!  Determines Hubbard U or Hund's J from series of *DS*_LRUJ.nc
!!  files containing information regarding the perturbation applied to a particular
!!  atom and the resulting occupations/magnetizations. The procedure implemented
!!  is that of the SCF linear response for the Hubbard U (Phys. Rev. B 71,035105)
!!  and the Hund's J (Phys. Rev. B 98, 235157) parameters.
!!  This protocol was coded up in November 2022 by Lorien MacEnulty (macenulty.com),
!!  doctoral researcher in the Quantum Theory of Materials group (theoryofmaterials.com)
!!  at Trinity College Dublin, headed by Dr. David O'Regan.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2022 ABINIT group (LMac)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! INPUTS
!!  -d <n> = Command line argument: highest degree of intended polynomial fits
!!  *DS*_LRUJ.nc files = gives data from perturbative Abinit calculations
!! OUTPUT
!!  std_out = log file
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


program lruj

 use defs_basis
 use m_xmpi
 use m_abicore
 use m_build_info
 use m_errors
 use m_argparse
 use m_crystal
 use netcdf
 use m_nctk

 use m_fstrings,    only : itoa, sjoin, ltoa
 use m_specialmsg,  only : specialmsg_getcount, herald
 use m_sort,        only : sort_dp
 use m_common,      only : crystal_from_file
 use m_mpinfo,      only : destroy_mpi_enreg, initmpi_seq
 use m_paw_uj,      only : pawuj_ini,pawuj_free,pawuj_det, macro_uj_type

 implicit none

!Local variables-------------------------------

!scalars
 integer,parameter                  :: master=0
 integer                            :: nproc,my_rank,comm
 logical                            :: iam_master

 integer                            :: ncid,nnat,prtvol,nargs,nfiles,ndtpawuj,degarg
 integer                            :: ndata,nspden,macro_uj,pawujat,dmatpuopt
 integer                            :: degree,mdegree,ii,ipert
 real(dp)                           :: diem,diemix,diemixmag,ph0phiint,signum,Ha2eV
 type(crystal_t)                    :: cryst

!arrays
 integer, allocatable               :: iperm(:),pawujat_file(:),macrouj_file(:),dmatpuopt_file(:)
 real(dp), allocatable              :: diem_file(:),ph0phiint_file(:),nspden_file(:)

 real(dp), allocatable              :: perts(:),occs0(:),occs(:)
 real(dp), allocatable              :: uj_perts(:),luocc(:,:),luocc_nnat(:,:)
 real(dp), allocatable              :: chi0coeffs(:),chicoeffs(:),chi0(:),chi(:),hubpar(:)
 real(dp), allocatable              :: chi0err(:),chierr(:),hubparerr(:)

!characters
 character(len=1)                   :: parname
 character(len=5)                   :: pertname,degreename
 character(len=12)                  :: regname
 character(len=14)                  :: occmag
 character(len=24)                  :: codename
 character(len=30)                  :: diem_token
 character(len=500)                 :: message,command,arg,msg
 character(len=fnlen),allocatable   :: file_paths(:)


!##########################################################################################################
!##################################  Set up MPI architecture (unused)  ####################################

 !Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 !Initialize MPI (not used but necessary)
 call xmpi_init()
 comm = xmpi_world
 nproc = xmpi_comm_size(comm)
 my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

 !Initialize memory profiling if it is activated
 !if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
 !note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 !No MPI functionality needed for main procedure.
 if (my_rank /= master) goto 100

!##########################################################################################################
!######################################  Read command line options  #######################################

 !Syntax: ./lruj  FILE1 FILE2 FILE3 ... [-d 5]
 !i.e. list of netcdf files come first, followed by options

 !Count arguments and number of files (= #perturbations)
 nargs = command_argument_count()
 ABI_MALLOC(file_paths, (nargs))
 nfiles = 0
 do ii=1,nargs
   call get_command_argument(ii, arg)
   if (arg(1:1) == "-") exit
   nfiles = nfiles + 1
   file_paths(nfiles) = trim(arg)
 end do

 !Assess options
 do ii=1,command_argument_count()
   call get_command_argument(ii, arg)
   if (arg == "--version") then
     write(std_out,"(a)") trim(abinit_version); goto 100
   else if (arg == "-h" .or. arg == "--help") then
     !Document the options.
     call lruj_show_help()
     goto 100
   end if
 end do

 !If no files found, exit program.
 if (nfiles == 0) then
   write(std_out, *) "Empty file list!"
   goto 100
 end if

 !Get other options from the CLI. e.g. --prtvol 0 -d 3.0
 !Should be documented in lruj_show_help
 ABI_CHECK(get_arg("prtvol", prtvol, msg, default=0) == 0, msg)
 ABI_CHECK(get_arg("d", degarg, msg, default=1) == 0, msg)

 !Print header
 codename='LRUJ'//repeat(' ',18)
 call herald(codename, abinit_version, std_out)

!##########################################################################################################
!######################################  Read *DSi*_LRUJ.nc files  ########################################

 ABI_MALLOC(uj_perts,(nfiles))
 ABI_MALLOC(macrouj_file, (nfiles))

 !Read perturbation strengths, parameter type, #spins and #atoms
 !from each file.
 do ii=1,nfiles
   NCF_CHECK(nctk_open_read(ncid, file_paths(ii), xmpi_comm_self))
   NCF_CHECK(nf90_get_var(ncid, vid("uj_pert"), uj_perts(ii)))
   NCF_CHECK(nf90_get_var(ncid, vid("macro_uj"), macrouj_file(ii)))
   macro_uj=macrouj_file(1)
   !Make sure ndtpawuj is always 4.
   NCF_CHECK(nctk_get_dim(ncid, "ndtpawuj", ndtpawuj))
   ABI_CHECK_IEQ(ndtpawuj, 4, "Wrong ndtpawuj")
   NCF_CHECK(nctk_get_dim(ncid, "nspden", nspden))
   NCF_CHECK(nctk_get_dim(ncid, "nnat", nnat))
   NCF_CHECK(nf90_close(ncid))
 end do

 !Sort files by perturbation magnitude.
 ABI_MALLOC(iperm, (nfiles))
 iperm = [(ii, ii=1,nfiles)]
 call sort_dp(nfiles, uj_perts, iperm, tol12)
 file_paths(1:nfiles) = file_paths(iperm(:))
 ABI_FREE(iperm)

 !Allocate main data-holding arrays.
 ABI_MALLOC(luocc, (ndtpawuj, nfiles))
 ABI_MALLOC(luocc_nnat, (ndtpawuj, nnat))
 ABI_MALLOC(pawujat_file, (nfiles))
 ABI_MALLOC(diem_file, (nfiles))
 ABI_MALLOC(dmatpuopt_file, (nfiles))
 ABI_MALLOC(ph0phiint_file, (nfiles))
 ABI_MALLOC(nspden_file, (nfiles))

 !Set macro_uj-specific variables, strings and constants.
 Ha2eV=27.2113961318d0
 if (macro_uj==4) then         !Calculation of the Hunds J parameter
   diem_token="diemixmag"     !Unscreened response in Hund's J impacted by diemixmag
   pertname='beta '           !Hund's J perturbation: +beta to spin up, -beta to down
   parname='J'
   occmag='Magnetizations'    !Magnetic moments are monitored.
   signum=-1.0d0              !Hund's J is -1*(1/chi0-1/chi)
 else
   diem_token="diemix"        !Unscreened response in Hubbard U impacted by diemix
   pertname='alpha'           !Hubbard U perturbation; applied equally to spins up and down
   parname='U'
   occmag='  Occupations'     !Total occupation is monitored.
   signum=1.0d0               !Hubbard U is 1*(1/chi0-1/chi)
 end if

 !Allocate perturbation and occupation arrays. Set the first
 !to the unperturbed case (i.e., when perturbation=0.0d0).
 ABI_MALLOC(perts,(0:nfiles))
 ABI_MALLOC(occs0,(0:nfiles))
 ABI_MALLOC(occs,(0:nfiles))
 perts(0)=0.0d0               !Unperturbed case.

 write(std_out,'(a,i2)') ' Number of perturbations detected: ',nfiles

 !Open _LRUJ.nc files and read in the relevant data.
 do ii=1,nfiles
   NCF_CHECK(nctk_open_read(ncid, file_paths(ii), xmpi_comm_self))
   NCF_CHECK(nf90_get_var(ncid, vid("pawujat"), pawujat_file(ii)))
   pawujat=pawujat_file(1)
   NCF_CHECK(nf90_get_var(ncid, vid("nspden"), nspden_file(ii)))
   nspden=nspden_file(1)
   NCF_CHECK(nf90_get_var(ncid, vid("luocc"), luocc_nnat))
   luocc(:,ii) = luocc_nnat(:,pawujat)
   NCF_CHECK(nf90_get_var(ncid, vid(diem_token), diem_file(ii)))
   diem=diem_file(1)
   NCF_CHECK(nf90_get_var(ncid, vid("dmatpuopt"), dmatpuopt_file(ii)))
   dmatpuopt=dmatpuopt_file(1)
   NCF_CHECK(nf90_get_var(ncid, vid("ph0phiint"), ph0phiint_file(ii)))
   ph0phiint=ph0phiint_file(1)
   NCF_CHECK(nf90_close(ncid))
   !Testing if the unperturbed occupancies are equal across each run. If they
   !aren't, then exit. If they are equal, then save them in appropriate arrays.
   if ((ii>1).and.((occs0(0)/=luocc(1,ii)).or.(occs(0)/=luocc(2,ii)))) then
     write(std_out,'(a)') "ERROR: Unperturbed ground state occupations across LRUJ datasets are not equal."
     write(std_out,'(a)') "Check the consistency of input variables in your perturbative calculations:"
     write(std_out,'(a)') "    1. Are they each reading in the same WFK file?"
     write(std_out,'(a)') "    2. Are macro-uj, pawujat, dmatpuopt, diemix(mag) consistent"
     write(std_out,'(a)') "       across all perturbations?"
     write(std_out,'(a)') "If not, relaunch perturbative Abinit calculations, then"
     write(std_out,'(2a)') "reexecute lruj utility. Exiting.",ch10
     goto 100
   else
     perts(ii)=uj_perts(ii)*Ha2eV
     occs0(0)=luocc(1,ii)
     occs(0)=luocc(2,ii)
     occs0(ii)=luocc(3,ii)
     occs(ii)=luocc(4,ii)
   end if
 end do

!##########################################################################################################
!####################################  Tests on input information  ########################################

 !Tests if we have enough data points (at least 3) to conduct distinct regression.
 ndata=nfiles+1
 write(std_out,'(a,i2,a)') ' Including unperturbed state, we have ',ndata,' data points.'
 if (ndata==0) then
   ABI_ERROR('No linear response data points found.')
 else if (ndata==1) then
   msg = sjoin('Only one data point found. This utility needs',ch10,&
    'at least three (3) data points (two non-zero perturbations and one unperturbed) to compute',ch10,&
    'the Hubbard parameter.')
  ABI_ERROR(msg)
 else if (ndata==2) then
   msg = sjoin('Only two data points found. The scalar Hubbard Parameter from',ch10,&
    'the two-point linear regression scheme has already been printed in your .abo file. Try bashing',ch10,&
    '==>   grep "two-point regression" <run_name.abo> ',ch10,'to find the result of this calculation.')
  ABI_ERROR(msg)
 end if

 !pawujat consistency check.
 if (any(pawujat_file /= pawujat_file(1))) then
   msg = sjoin("Found different values of pawujat in files: ", ltoa(pawujat_file),ch10,&
         "Perturbed atom has to be consistent across perturbations to compute U or J.")
  ABI_ERROR(msg)
 end if

 !dmatpuopt consistency check
 if (any(dmatpuopt_file /= dmatpuopt_file(1))) then
   msg = sjoin("Found different values of dmatpuopt in files: ", ltoa(dmatpuopt_file),ch10,&
           "PAW projector must be consistent across perturbations to compute U or J.")
  ABI_ERROR(msg)
 end if

 !macro_uj consistency check
 if (any(macrouj_file /= macrouj_file(1))) then
   msg = sjoin("Found different values of macro_uj in files: ",ltoa(macrouj_file),ch10,&
           "Perturbation protocol and occupancy monitoring must be consistent to compute U or J.")
  ABI_ERROR(msg)
 end if

 !diemix/diemixmag consistency check
 if (any(diem_file /= diem_file(1))) then
   msg = sjoin("Found different values of mixing constant in files: ",ltoa(diem_file),ch10,&
          "Unscreened response functions will factor into U (J) incorrectly.")
  ABI_ERROR(msg)
 end if

 !Tests consistency of macro_uj, then writes message about macro_uj procedure selected.
 !Also assigns Hubbard parameter-specific variables.
 if (nspden==1) then
   write(message,'(a)') ' Determination of U-parameter for unpolarized structure (non standard)'
 else if (macro_uj==1.and.nspden==2) then
   write(message,'(a)') ' Standard determination of the Hubbard U parameter.'
 else if (macro_uj==2.and.nspden==2) then
   write(message,'(a)') ' Determination of parameter on single spin channel (experimental)'
   pertname='Pert. '
 else if (macro_uj==3.and.nspden==2) then
   parname='J'
   pertname='Pert. '
   write(message,'(a)') ' Determination of (not Hunds) J-parameter on single spin channel (experimental)'
 else if (macro_uj==4.and.nspden==2) then
   write(message,'(a)') ' Hunds J determination, implemented by L. MacEnulty August 2021'
 end if
 call wrtout(std_out,message,'COLL')

 !Tests compatibility of nspden and macro_uj
 if (macro_uj>1.and.nspden==1) then
   msg = sjoin('U on a single spin channel (or J) can only be determined for nspden=2 ,',ch10,&
    'Cannot calculate the chosen Hubbard parameter.')
   ABI_ERROR(msg)
 end if

 !Tests if perturbations are too small.
 if (maxval(abs(uj_perts))<0.00000001) then
   msg = sjoin('Perturbation magnitudes are too small.',ch10,&
     'Rerun perturbative Abinit calculations with pawujv >> 1d-8.')
   ABI_ERROR(msg)
 end if

!##########################################################################################################
!###############################  Calculation of the Response Functions  ##################################

 !Test compatibility of polynomial degree (if present as an argument) with
 !number of data points. Otherwise, default to the following:
 !If we have 3 data points, conduct at maximum a linear regression.
 !If 4 data points, conduct linear and quadratic regressions.
 !If more than 5 data points, conduct linear, quadratic and cubic regressions.
 if (degarg/=1) then
   if (degarg>ndata-2) then
     write(std_out,'(4a,i2,3a,i2,2a)') ch10,' ERROR: Your chosen polynomial degree is too large. The resulting',ch10,&
&       ' parameters will certainly be overfitted. Either conduct ',degarg+1,' perturbations,',ch10,&
&       ' or execute this utility again with --d',ndata-2,' or smaller. Exiting program.',ch10
     goto 100
   else
     mdegree=degarg
   end if
 else !Default max polynomial degree
   if (ndata<=4) then
     mdegree=ndata-2
   else
     mdegree=3
   end if
 end if
 write(std_out,'(a,i2)') ' Maximum degree of polynomials analyzed: ',mdegree

 !Write warning about response matrices
 write(std_out,'(5a)') ' NOTE: Unlike the ujdet utility, lruj treats the ',ch10,&
' response functions as scalars, not matrices!',ch10,&
' See lruj tutorial for more information.'

 !Allocate the response and error arrays
 ABI_MALLOC(chi0err,(mdegree))
 ABI_MALLOC(chierr,(mdegree))
 ABI_MALLOC(chi0,(mdegree))
 ABI_MALLOC(chi,(mdegree))
 ABI_MALLOC(hubpar,(mdegree))
 ABI_MALLOC(hubparerr,(mdegree))

 !For all regressions, call subroutine to calculate polynomial fit for chi0 and chi.
 do degree=1,mdegree
   ABI_MALLOC(chi0coeffs,(degree+1))
   ABI_MALLOC(chicoeffs,(degree+1))
   call polynomial_regression(ndata,perts,occs0,degree,chi0coeffs,chi0err(degree))
   call polynomial_regression(ndata,perts,occs,degree,chicoeffs,chierr(degree))
   chi0(degree)=chi0coeffs(2)/diem         !The derivative of all polynomial regressions
   chi(degree)=chicoeffs(2)                !at pert=0.0 is just the second coefficient.
   chi0err(degree)=chi0err(degree)/diem    !Chi0 error divided by diem also.
   hubpar(degree)=signum*(1.0d0/chi0(degree)-1.0d0/chi(degree))
   hubparerr(degree)=sqrt(chi0err(degree)/chi0(degree)**2+chierr(degree)/chi(degree)**2)
   ABI_FREE(chi0coeffs)
   ABI_FREE(chicoeffs)
 end do

!##########################################################################################################
!#############################  Printing information on Hubbard Parameters ################################

 !Printing relevant information about the Hubbard parameter just calculated.
 write(message,'(3a)') ch10,ch10,'*************************************************************************************************'
 call wrtout(std_out,message,'COLL')
 write(message,'(4a)') '**************************************  Linear Response ',parname,'  **************************************',ch10
 call wrtout(std_out,message,'COLL')
 write(message, '(a,i4)' ) ' Index of perturbed atom: ',pawujat
 call wrtout(std_out,message,'COLL')
 write(message, '(a,i4)' ) ' Value of macro_uj:  ',macro_uj
 call wrtout(std_out,message,'COLL')
 write(message, '(a,i4)' ) ' Value of dmatpuopt:  ',dmatpuopt
 call wrtout(std_out,message,'COLL')
 write(message, '(a,f6.3)' ) ' Mixing constant factored out of Chi0:  ',diem
 call wrtout(std_out,message,'COLL')
 write(message, '(a,f12.5,2a)' ) ' Percentage of AE orbital within the PAW sphere of perturbed subspace: ',ph0phiint*100.00,'%',ch10
 call wrtout(std_out,message,'COLL')

 write(message, fmt='(10a)')'  Perturbations         ',occmag,ch10,&
' --------------- -----------------------------',ch10,&
'    ',pertname,' [eV]     Unscreened      Screened',ch10,&
' --------------- -----------------------------'
 call wrtout(std_out,message,'COLL')
 do ipert=0,nfiles
   write(message, fmt='(3f15.10)') perts(ipert),occs0(ipert),occs(ipert)
   call wrtout(std_out,message,'COLL')
 end do

 write(message, fmt='(11a)') '                                                                       RMS Errors',&
&   ch10,'                                                         ---------------------------------------',ch10,&
&   ' Regression   Chi0 [eV^-1]   Chi [eV^-1]      ',parname,' [eV]    | Chi0 [eV^-1]  Chi [eV^-1]     ',parname,&
&   ' [eV]',ch10,'--------------------------------------------------------|---------------------------------------'
 call wrtout(std_out,message,'COLL')
 do degree=1,mdegree
   if (degree==1) then
     regname=' Linear:    '
   else if (degree==2) then
     regname=' Quadratic: '
   else if (degree==3) then
     regname=' Cubic:     '
   else
     write(degreename,'(i2)') degree
     regname=' Degree'//degreename//': '
   end if
   write(message,fmt='(a,3f14.7,a,3f13.7)') regname,chi0(degree),chi(degree),hubpar(degree),&
&     '  |',chi0err(degree),chierr(degree),hubparerr(degree)
   call wrtout(std_out,message,'COLL')
 end do

 write(message,'(3a)') '*************************************************************************************************',ch10,&
'*************************************************************************************************'
 call wrtout(std_out,message,'COLL')


!##########################################################################################################
!############################################  Deallocations ##############################################


 ABI_FREE(perts)
 ABI_FREE(occs0)
 ABI_FREE(occs)
 ABI_FREE(luocc_nnat)
 ABI_FREE(pawujat_file)
 ABI_FREE(diem_file)
 ABI_FREE(dmatpuopt_file)
 ABI_FREE(ph0phiint_file)
 ABI_FREE(macrouj_file)
 ABI_FREE(nspden_file)
 ABI_FREE(chi0err)
 ABI_FREE(chierr)
 ABI_FREE(chi0)
 ABI_FREE(chi)
 ABI_FREE(hubpar)
 ABI_FREE(hubparerr)
 ABI_FREE(luocc)
 ABI_FREE(uj_perts)
 ABI_FREE(file_paths)

 !Ending herald.
 write(std_out,*) ch10,'Linear Response UJ (LRUJ) program complete. Live long and prosper. ~LMac',ch10
 goto 100

 100 call xmpi_end()



!##########################################################################################################
!#####################################  Subroutines and Functions  ########################################

contains

!!****f* abitk/lruj_show_help
!! NAME
!! lruj_show_help
!!
!! FUNCTION
!!  Show command line help
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine lruj_show_help()

  write(std_out,"(a)")" "
  write(std_out,"(a)")"           Linear Response Hubbard U and Hund's J (LRUJ) Utility"
  write(std_out,"(a)")"-----------------------------------------------------------------------------"
  write(std_out,"(a)")"To execute the LRUJ utility, execute: "
  write(std_out,"(a)")"   ./lruj  FILE1 FILE2 FILE3 ... [options]"
  write(std_out,"(2a)")"             ^ input files must be _LRUJ.nc from Abinit run",ch10
  write(std_out,"(a)")" --version              Show version number and exit."
  write(std_out,"(a)")" -h, --help             Show this help and exit."
  write(std_out,"(a)")" --d <n>                Set the maximum degree n polynomial calculated for"
  write(std_out,"(a)")"                           the response functions chi and chi0."
  write(std_out,"(a)")"                           (i.e., 1=linear, 2=quadratic, 3=cubic, etc.)"
  write(std_out,"(a)")"                           NOTE: a degree n polynomial will require at minimum"
  write(std_out,"(a)")"                           n+2 points (n+1 perturbations and the unperturbed"
  write(std_out,"(a)")"                           case) or more so as to avoid overfitting."

end subroutine lruj_show_help
!!***


 !Function to simplify reading in of variables from netcdf files.
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid



!----------------------------------------------------------------------
!!****f* polynomial_regression
!! NAME
!!  polynomial_regression
!!
!! FUNCTION
!!  Perform a polynomial regression on incoming data points, the
!!  x-values of which are stored in array xvals and the y-values
!!  stored in array yvals. Returns a one dimensional array with
!!  fit coefficients (coeffs) and the unbiased RMS error of the
!!  fit as a scalar (RMSerr).
!!
!! INPUTS
!!  npoints = number of data points
!!  xvals(npoints) = x-values of those data points
!!  yvals(npoints) = y-values of those data points
!!  degree = order of the polynomial
!!
!! OUTPUT
!!  coeffs(degree+1) = coefficients of the polynomial regression
!!  RMSerr = unbiased RMS error on the fit
!!            RMSerr=\sqrt{\frac{1}{npoints-1}*
!!                      \sum_i^npoints{(fitval-yvals(i))**2}}
!!
!! SOURCE
!!  Polynomial regression algorithm from Rosetta Code under Creative Commons
!!  and GNU Free Documentation License.
!!  Link: https://rosettacode.org/wiki/Polynomial_regression#Fortranf
!!  Some variables changed to simplify.

subroutine polynomial_regression(npoints,xvals,yvals,degree,coeffs,RMSerr)

!Arguments ------------------------------------

!scalars
 integer                     :: npoints,degree
 real(dp),intent(out)        :: RMSerr
!arrays
 real(dp),intent(in)         :: xvals(1:npoints),yvals(1:npoints)
 real(dp),intent(out)        :: coeffs(degree+1)

!Local variables-------------------------------

!scalars
 integer                     :: ncoeffs,icoeff,ipoint,info
 real(dp)                    :: residual,fitval
!arrays
 integer,allocatable         :: ipiv(:)
 real(dp),allocatable        :: work(:)
 real(dp),allocatable        :: A(:,:),AT(:,:),ATA(:,:)
!characters
 character(len=500)          :: message

!####################################################################
!#####################  Get Polynomial Fit  #########################

  ncoeffs=degree+1

  ABI_MALLOC(ipiv,(ncoeffs))
  ABI_MALLOC(work,(ncoeffs))
  ABI_MALLOC(A,(size(xvals),ncoeffs))
  ABI_MALLOC(AT,(ncoeffs,size(xvals)))
  ABI_MALLOC(ATA,(ncoeffs,ncoeffs))

  !Prepare the matrix A
  do icoeff=0,ncoeffs-1
    do ipoint=1,size(xvals)
       A(ipoint,icoeff+1) = xvals(ipoint)**icoeff
    end do
  end do

  AT  = transpose(A)
  ATA = matmul(AT,A)

  !Call LAPACK subroutines DGETRF and DGETRI
  call DGETRF(ncoeffs,ncoeffs,ATA,ncoeffs,ipiv,info)
  if (info/=0) then
    ABI_ERROR('Problem returned from LAPACK DGETRF in polynomial regression routine.')
    return
  end if
  call DGETRI(ncoeffs,ATA,ncoeffs,ipiv,work,ncoeffs,info)
  if (info/=0) then
    ABI_ERROR('Problem returned from LAPACK DGETRI in polynomial regression routine.')
    return
  end if

  coeffs = matmul(matmul(ATA,AT),yvals)

!####################################################################
!##############  RMS error on the polynomial fit  ###################

  residual=0.0d0
  do ipoint=1,npoints
    fitval=0.0d0
    do icoeff=1,ncoeffs
      fitval=fitval+coeffs(icoeff)*xvals(ipoint)**(icoeff-1)
    end do
    residual=residual+(fitval-yvals(ipoint))**2
  end do
  RMSerr=sqrt(residual/(real(npoints-1,8)))


!####################################################################
!########################  Deallocations  ###########################

  ABI_FREE(ipiv)
  ABI_FREE(work)
  ABI_FREE(A)
  ABI_FREE(AT)
  ABI_FREE(ATA)

end subroutine polynomial_regression
!!***


end program lruj
!!***
