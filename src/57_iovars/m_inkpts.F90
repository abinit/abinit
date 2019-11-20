!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_inkpts
!! NAME
!!  m_inkpts
!!
!! FUNCTION
!!  Routines to initialize k-point and q-point sampling from input file.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR)
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

module m_inkpts

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr

 use m_time,      only : timab
 use m_fstrings,  only : sjoin, itoa
 use m_numeric_tools, only : isdiagmat
 use m_geometry,  only : metric
 use m_symfind,   only : symfind, symlatt
 use m_cgtools,   only : set_istwfk
 use m_parser,    only : intagm
 use m_kpts,      only : getkgrid, testkgrid, mknormpath

 implicit none

 private
!!***

 public :: inkpts      ! Initialize k points (list of k points, weights, storage) for one particular dataset.
 public :: inqpt       ! Initialize the q point for one particular dataset, characterized by jdtset.
!!***

contains
!!***

!!****f* m_inkpts/inkpts
!! NAME
!! inkpts
!!
!! FUNCTION
!! Initialize k points (list of k points, weights, storage)
!! for one particular dataset, characterized by jdtset.
!! Note that nkpt (and nkpthf) can be computed by calling this routine with input value of nkpt=0, provided kptopt /= 0.
!!
!! INPUTS
!! bravais(11): bravais(1)=iholohedry
!!              bravais(2)=center
!!              bravais(3:11)=coordinates of rprim in the axes of the conventional bravais lattice (*2 if center/=0)
!! chksymbreak= if 1, will check whether the k point grid is symmetric, and stop if not.
!! [impose_istwf_1]= (optional argument):
!!                 0: no restriction on istwfk
!!                 1: impose istwfk=1 for all k points
!!                 2: impose istwfk=1 for all k points non equal to zero
!! iout=unit number for echoed output
!! iscf= <= 0 => non-SCF, >0 => SCF.
!! jdtset=number of the dataset looked for
!! lenstr=actual length of the string
!! kptopt=option for the generation of k points
!! msym=default maximal number of symmetries
!! getkerange_path= Path of KERANGE.nc file used to initialize k-point sampling if kptopt == 0 and string != ABI_NOFILE
!! nqpt=number of q points (0 or 1)
!! nsym=number of symetries
!! occopt=option for occupation numbers
!! qptn(3)=reduced coordinates of eventual q point shift (already normalized).
!! response=0 if GS case, =1 if RF case.
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! string*(*)=character string containing all the input data. Initialized previously in instrng.
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)=symmetry operations in real space in terms of primitive translations
!! vacuum(3)=for each direction, 0 if no vacuum, 1 if vacuum
!! comm= MPI communicator
!!
!! OUTPUT
!! fockdownsampling(3)=echo of input variable fockdownsampling(3)
!! kptnrm=normalisation of k points
!! kptrlatt_orig(3,3)=Original value of kptrlatt as specified in the input file (if kptopt/=0)
!! kptrlatt(3,3)=k-point lattice specification (if kptopt/=0)
!! kptrlen=length of the smallest real space supercell vector
!! nshiftk_orig=Original number of k-point shifts (0 if not read)
!! nshiftk=actual number of k-point shifts in shiftk (if kptopt/=0)
!! shiftk(3,MAX_NSHIFTK)=shift vectors for k point generation (if kptopt/=0)
!! If nkpt/=0  the following arrays are also output:
!!   istwfk(nkpt)=option parameters that describes the storage of wfs
!!   kpt(3,nkpt)=reduced coordinates of k points.
!!   kpthf(3,nkpthf)=reduced coordinates of k points for Fock operator.
!!   wtk(nkpt)=weight assigned to each k point.
!! ngkpt(3)=Number of divisions along the three reduced directions
!!   (0 signals that this variable has not been used.
!! shiftk_orig(3,MAX_NSHIFTK)=Original shifts read from the input file
!!   (0 signals that this variable has not been read).
!!
!! SIDE EFFECTS
!! Input/output:
!! nkpt=number of k points
!!  if non-zero at input, is only an input variable
!!  if zero at input, its actual value will be computed
!! nkpthf=number of k points for Fock operator, computed if nkpt=0 at input
!!
!! NOTES
!! Warning: this routine can be called with nkpt=0 (in which
!! case it returns the true value of nkpt), which can lead
!! to strange bugs in the debugging procedure, if one tries to print wtk or istwfk, in this case!
!!
!! PARENTS
!!      invars1,invars2
!!
!! CHILDREN
!!      getkgrid,intagm,metric,mknormpath,testkgrid,timab,wrtout
!!
!! SOURCE

subroutine inkpts(bravais,chksymbreak,fockdownsampling,iout,iscf,istwfk,jdtset,&
& kpt,kpthf,kptopt,kptnrm,kptrlatt_orig,kptrlatt,kptrlen,lenstr,msym, getkerange_path, &
& nkpt,nkpthf,nqpt,ngkpt,nshiftk,nshiftk_orig,shiftk_orig,nsym,&
& occopt,qptn,response,rprimd,shiftk,string,symafm,symrel,vacuum,wtk,comm,&
& impose_istwf_1) ! Optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chksymbreak,iout,iscf,jdtset,kptopt,lenstr,msym,nqpt,nsym,occopt
 integer,intent(in) :: response, comm
 integer,intent(in),optional :: impose_istwf_1
 integer,intent(inout) :: nkpt,nkpthf
 integer,intent(out) :: nshiftk,nshiftk_orig
 integer,intent(out) :: fockdownsampling(3)
 real(dp),intent(out) :: kptnrm,kptrlen
 character(len=*),intent(in) :: string
 character(len=*),intent(in) :: getkerange_path
!arrays
 integer,intent(in) :: bravais(11),symafm(msym),symrel(3,3,msym),vacuum(3)
 integer,intent(out) :: istwfk(nkpt),kptrlatt(3,3),kptrlatt_orig(3,3),ngkpt(3)
 real(dp),intent(in) :: rprimd(3,3),qptn(3)
 real(dp),intent(out) :: kpt(3,nkpt),kpthf(3,nkpthf),shiftk(3,MAX_NSHIFTK),wtk(nkpt),shiftk_orig(3,MAX_NSHIFTK)

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: dkpt,ii,ikpt,jkpt,marr,ndiv_small,nkpt_computed,my_rank,nprocs
 integer :: nsegment,prtkpt,tread,tread_kptrlatt,tread_ngkpt, ncid, fform, ierr
 real(dp) :: fraction,norm,ucvol,wtksum
 character(len=500) :: msg
 type(hdr_type) :: hdr
!arrays
 integer,allocatable :: ndivk(:),intarr(:), krange2ibz(:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),rmet(3,3),tsec(2)
 real(dp),allocatable :: kptbounds(:,:),dprarr(:)

! *************************************************************************

 call timab(192,1,tsec)

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 ! Compute the maximum size of arrays intarr and dprarr
 marr = max(3*nkpt,3*MAX_NSHIFTK)
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

 ! Use zero to signal that these values have not been read.
 ngkpt = 0
 shiftk_orig = zero
 kptrlatt_orig = 0; kptrlatt = 0
 nshiftk_orig = 1; nshiftk = 1

 !fockdownsampling(:)=1
 !kptnrm = one
 !kpthf = zero

 ! MG: FIXME These values should be initialized because they are intent(out)
 ! but several tests fails. So we keep this bug to avoid problems somewhere else
 ! The initialization of the kpoints should be rewritten in a cleaner way
 ! without all these side effects!
 !shiftk = zero
 ! Initializing these three variables is OK but we keep the bug to preserve the old behavior
 !wtk = one
 !kpt = zero
 !istwfk = 1

 ! Initialize kptrlen
 kptrlen=30.0_dp
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'kptrlen',tread,'DPR')
 if(tread==1)kptrlen=dprarr(1)

 ! Initialize kpt, kptnrm and wtk according to kptopt.
 if (kptopt == 0 .and. getkerange_path == ABI_NOFILE) then
   ! For kptopt==0, one must have nkpt defined.
   kpt(:,:)=zero
   call intagm(dprarr,intarr,jdtset,marr,3*nkpt,string(1:lenstr),'kpt',tread,'DPR')
   if(tread==1) kpt(:,:)=reshape( dprarr(1:3*nkpt), [3,nkpt])

   kptnrm=one
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'kptnrm',tread,'DPR')
   if(tread==1) kptnrm=dprarr(1)

   ! Only read wtk when iscf >0 or iscf=-1 or iscf=-3 or (iscf=-2 and response=1)
   ! (this last option is for Zach Levine)
   ! Normalize the k-point weights when occopt/=2
   ! Check that k point weights add to 1 when occopt==2
   if  (iscf>0.or.iscf==-1.or.iscf==-3.or.(iscf==-2.and.response==1))  then
     call intagm(dprarr,intarr,jdtset,marr,nkpt,string(1:lenstr),'wtk',tread,'DPR')
     if(tread==1) wtk(1:nkpt)=dprarr(1:nkpt)

     wtksum=sum(wtk(:))
     write(msg,'(a,i0,a,f12.6)')' inkpts: Sum of ',nkpt,' k point weights is',wtksum
     call wrtout(std_out,msg,'COLL')

     if (wtksum < tol6) then
       write(msg, '(3a)' )&
       'This sum is too close to zero. ',ch10,&
       'Action: correct the array wtk in the input file.'
       MSG_ERROR(msg)
     end if
     if (abs(wtksum - one) > tol6) then
       if (occopt==2) then
         write(msg, '(a,1p,e18.8,a,a,a)' )&
          'wtksum= ',wtksum,' /= 1.0 means wts do not add to 1 , while occopt=2.',ch10,&
          'Action: correct the array wtk in input file.'
         MSG_ERROR(msg)
       else
         write(msg,'(a,i0,a)')' With present occopt= ',occopt,', renormalize it to one'
         call wrtout(std_out,msg,'COLL')
         norm=one/wtksum
         wtk(1:nkpt)=wtk(1:nkpt)*norm
       end if
     end if
   end if

 else if (kptopt == 0 .and. getkerange_path /= ABI_NOFILE) then
   ! Initialize kpts from kerange_path file.
   ABI_MALLOC(krange2ibz, (nkpt))
   if (my_rank == master) then
#ifdef HAVE_NETCDF
     NCF_CHECK(nctk_open_read(ncid, getkerange_path, xmpi_comm_self))
     call hdr_ncread(hdr, ncid, fform)
     ABI_CHECK(fform == fform_from_ext("KERANGE.nc"), sjoin("Error while reading:", getkerange_path, ", fform:", itoa(fform)))
     ! TODO Add code for consistency check
     !kptopt, nsym, occopt
     !ABI_CHECK(nkpt == hdr%nkpt, "nkpt from kerange != nkpt")
     NCF_CHECK(nf90_get_var(ncid, nctk_idname(ncid, "krange2ibz"), krange2ibz))
     NCF_CHECK(nf90_close(ncid))
#else
     MSG_ERROR("getkerange_path requires NETCDF support")
#endif
   end if
   call xmpi_bcast(krange2ibz, master, comm, ierr)
   call hdr%bcast(master, my_rank, comm)
   ! Hdr contains kpts in the IBZ. Extract data corresponding to pockets via krange2ibz.
   nshiftk = hdr%nshiftk; nshiftk_orig = hdr%nshiftk_orig
   istwfk = hdr%istwfk(krange2ibz(:))
   kptrlatt = hdr%kptrlatt; kptrlatt_orig = hdr%kptrlatt_orig
   ABI_CHECK(isdiagmat(hdr%kptrlatt), "kptrlatt is not diagonal!")
   ngkpt(1) = hdr%kptrlatt(1, 1); ngkpt(2) = hdr%kptrlatt(2, 2); ngkpt(3) = hdr%kptrlatt(3, 3)
   kpt = hdr%kptns(:, krange2ibz(:)) !; kpthf(3,nkpthf)
   shiftk(:,1:nshiftk) = hdr%shiftk; shiftk_orig(:, 1:nshiftk_orig) = hdr%shiftk_orig
   wtk = hdr%wtk(krange2ibz(:))
   call hdr%free()
   kptnrm = one
   ABI_FREE(krange2ibz)

 else if (kptopt < 0) then
   ! Band structure calculation
   nsegment=abs(kptopt)

   if (iscf /= -2)then
     write(msg,  '(3a,i0,3a)' ) &
      'For a negative kptopt, iscf must be -2,',ch10,&
      'while it is found to be ',iscf,'.',ch10,&
      'Action: change the value of iscf in your input file, or change kptopt.'
     MSG_ERROR(msg)
   end if

   if(marr<3*nsegment+3)then
     marr=3*nsegment+3
     ABI_DEALLOCATE(dprarr)
     ABI_DEALLOCATE(intarr)
     ABI_ALLOCATE(dprarr,(marr))
     ABI_ALLOCATE(intarr,(marr))
   end if

   ABI_ALLOCATE(kptbounds,(3,nsegment+1))
   ABI_ALLOCATE(ndivk,(nsegment))

   call intagm(dprarr,intarr,jdtset,marr,3*nsegment+3,string(1:lenstr),'kptbounds',tread,'DPR')

   if(tread==1)then
     kptbounds(:,:)=reshape( dprarr(1:3*nsegment+3), [3,nsegment+1])
   else
     write(msg,'(5a)') &
     'When kptopt is negative, kptbounds must be initialized ',ch10,&
     'in the input file, which is not the case.',ch10,&
     'Action: initialize kptbounds in your input file, or change kptopt.'
     MSG_ERROR(msg)
   end if

   call intagm(dprarr,intarr,jdtset,marr,nsegment,string(1:lenstr),'ndivk',tread,'INT')
   if(tread==1)then
     ndivk(1:nsegment)=intarr(1:nsegment)
     ! The 1 stand for the first point
     nkpt_computed=1+sum(ndivk(1:nsegment))

     ! ndivk and ndivsm are mutually exclusive
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ndivsm',tread,'INT')
     if (tread == 1) then
       MSG_ERROR("ndivk and ndivsm are mutually exclusive. Choose only one variable")
     end if

   else
     ! Calculate ndivk such as the path is normalized
     ! Note that if both ndivk and ndivsm are defined in in input file, only ndivk is used !
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ndivsm',tread,'INT')
     if(tread==1)then
       ndiv_small=intarr(1)
       call metric(gmet,gprimd,std_out,rmet,rprimd,ucvol)
       call mknormpath(nsegment+1,kptbounds,gmet,ndiv_small,ndivk,nkpt_computed)
     else
       write(msg,'(5a)') &
        'When kptopt is negative, ndivsm or ndivk must be initialized ',ch10,&
        'in the input file, which is not the case.',ch10,&
        'Action: initialize ndivsm or ndivk in your input file, or change kptopt.'
       MSG_ERROR(msg)
     end if
   end if

   ! Check that the argument nkpt is coherent with nkpt_computed
   if (nkpt/=0 .and. nkpt /= nkpt_computed) then
     write(msg,  '(a,i0,5a,i0,7a)' ) &
     'The argument nkpt = ',nkpt,', does not match',ch10,&
     'the number of k points generated by kptopt, ndivk, kptbounds,',ch10,&
     'and the eventual symmetries, that is, nkpt= ',nkpt_computed,'.',ch10,&
     'However, note that it might due to the user,',ch10,&
     'if nkpt is explicitely defined in the input file.',ch10,&
     'In this case, please check your input file.'
     MSG_ERROR(msg)
   end if

   if (nkpt/=0) then
     ! The array kpt has the right dimension and we can generate the k-path
     call intagm(dprarr,intarr,jdtset,marr,3*nsegment+3,string(1:lenstr),'kptbounds',tread,'DPR')
     if(tread==1)then
       kptbounds(:,:)=reshape( dprarr(1:3*nsegment+3), [3,nsegment+1])
     else
       write(msg, '(5a)') &
       'When kptopt is negative, kptbounds must be initialized ',ch10,&
       'in the input file, which is not the case.',ch10,&
       'Action: initialize kptbounds in your input file, or change kptopt.'
       MSG_ERROR(msg)
     end if

     ! First k point
     jkpt=1
     kpt(:,1)=kptbounds(:,1)
     do ii=1,nsegment
       dkpt=ndivk(ii)
       do ikpt=1,dkpt
         fraction=dble(ikpt)/dble(dkpt)
         kpt(:,ikpt+jkpt)=fraction *kptbounds(:,ii+1)+(one-fraction)*kptbounds(:,ii)
       end do
       jkpt=jkpt+dkpt
     end do

   end if

   kptnrm=one
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'kptnrm',tread,'DPR')
   if(tread==1) kptnrm=dprarr(1)

   ABI_DEALLOCATE(kptbounds)
   ABI_DEALLOCATE(ndivk)

 else if (kptopt>=1 .and. kptopt<=4) then
   ! Read ngkpt
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ngkpt',tread_ngkpt,'INT')
   if(tread_ngkpt==1)then
     ngkpt(1:3)=intarr(1:3)
     do ii=1,3
       if(ngkpt(ii)<1)then
         write(msg,'(a,i0,3a,i0,3a)') &
         'The input variable ngkpt(',ii,') must be strictly positive,',ch10,&
         'while it is found to be',ngkpt(ii),'.',ch10,&
         'Action: change it in your input file, or change kptopt.'
         MSG_ERROR(msg)
       end if
     end do
   end if

   call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),'kptrlatt',tread_kptrlatt,'INT')
   if(tread_kptrlatt==1) kptrlatt(:,:)=reshape(intarr(1:9), [3,3])

   if(tread_ngkpt==1 .and. tread_kptrlatt==1)then
     write(msg, '(5a)') &
     'The input variables ngkpt and kptrlatt cannot both ',ch10,&
     'be defined in the input file.',ch10,&
     'Action: change one of ngkpt or kptrlatt in your input file.'
     MSG_ERROR(msg)
   else if(tread_ngkpt==1)then
     kptrlatt(:,:)=0
     kptrlatt(1,1)=ngkpt(1)
     kptrlatt(2,2)=ngkpt(2)
     kptrlatt(3,3)=ngkpt(3)
     ! Save kptrlatt for reference.
     kptrlatt_orig = kptrlatt
   end if

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nshiftk',tread,'INT')
   if(tread==1)nshiftk=intarr(1)

   if (nshiftk < 1 .or. nshiftk > MAX_NSHIFTK )then
     write(msg,  '(a,i0,2a,i0,3a)' )&
     'The only allowed values of nshiftk are between 1 and ',MAX_NSHIFTK,ch10,&
     'while it is found to be',nshiftk,'.',ch10,&
     'Action: change the value of nshiftk in your input file, or change kptopt.'
     MSG_ERROR(msg)
   end if

   call intagm(dprarr,intarr,jdtset,marr,3*nshiftk,string(1:lenstr),'shiftk',tread,'DPR')
   if(tread==1)then
     shiftk(:,1:nshiftk)=reshape( dprarr(1:3*nshiftk), [3,nshiftk])
     ! Save input shifts as they will be changes in getkgrid.
     nshiftk_orig = nshiftk
     shiftk_orig(:,1:nshiftk) = shiftk(:,1:nshiftk)
   else
     if(nshiftk/=1)then
       write(msg,  '(3a,i0,2a)' )&
       'When nshiftk is not equal to 1, shiftk must be defined in the input file.',ch10,&
       'However, shiftk is not defined, while nshiftk=',nshiftk,ch10,&
       'Action: change the value of nshiftk in your input file, or define shiftk.'
       MSG_ERROR(msg)
     end if
     ! Default values used in indefo
     nshiftk_orig = 1
     shiftk_orig(:,1:nshiftk) = half
   end if

   prtkpt=0
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtkpt',tread,'INT')
   if(tread==1)prtkpt=intarr(1)

   if(sum(abs(kptrlatt(:,:)))==0)then
     kptrlatt_orig = 0
     do ii=1,3
       kptrlatt_orig(ii,ii) = ngkpt(ii)
     end do
     ! The parameters of the k lattice are not known, compute kptrlatt, nshiftk, shiftk.
     call testkgrid(bravais,iout,kptrlatt,kptrlen, msym,nshiftk,nsym,prtkpt,rprimd,shiftk,symafm,symrel,vacuum)
   end if

   fockdownsampling(:)=1
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'fockdownsampling',tread,'INT')
   if(tread==1)fockdownsampling=intarr(1:3)

   call getkgrid(chksymbreak,0,iscf,kpt,kptopt,kptrlatt,kptrlen,&
    msym,nkpt,nkpt_computed,nshiftk,nsym,rprimd,&
    shiftk,symafm,symrel,vacuum,wtk,nkpthf=nkpthf,kpthf=kpthf,downsampling=fockdownsampling)

   kptnrm=one

 else
   write(msg,'(3a,i0,3a)' ) &
   'The only values of kptopt allowed are smaller than 4.',ch10,&
   'The input value of kptopt is: ',kptopt,'.',ch10,&
   'Action: change kptopt in your input file.'
   MSG_ERROR(msg)
 end if

 if (kptnrm < tol10) then
   write(msg, '(5a)' )&
   'The input variable kptnrm is lower than 1.0d-10,',ch10,&
   'while it must be a positive, non-zero number.   ',ch10,&
   'Action: correct the kptnrm in the input file.'
   MSG_ERROR(msg)
 end if

 ! The k point number has been computed, and, if nkpt/=0, also the list of k points.
 ! Also nkpthf has been computed, and, if nkpt/=0, also the list kpthf.
 ! Now, determine istwfk, and eventually shift the k points by the value of qptn.
 if (nkpt /= 0) then
   istwfk(1:nkpt)=0
   call intagm(dprarr,intarr,jdtset,marr,nkpt,string(1:lenstr),'istwfk',tread,'INT')
   if(tread==1) istwfk(1:nkpt)=intarr(1:nkpt)

   ! Impose istwfk=1 for RF calculations
   if(response==1)istwfk(1:nkpt)=1

   ! Also impose istwfk=1 for spinor calculations
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nspinor',tread,'INT')
   if(tread/=0 .and. intarr(1)/=1)istwfk(1:nkpt)=1

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawspnorb',tread,'INT')
   if(tread/=0 .and. intarr(1)/=0)istwfk(1:nkpt)=1

   do ikpt=1,nkpt
     if(istwfk(ikpt)==0)then
       kpoint=kpt(:,ikpt)/kptnrm; if (nqpt/=0.and.response==0) kpoint=kpoint+qptn
       istwfk(ikpt) = set_istwfk(kpoint)
     end if
     if (present(impose_istwf_1)) then
       if (impose_istwf_1==1) then
         istwfk(ikpt)=1
       else if (impose_istwf_1==2.and.any(kpt(:,ikpt)>tol10)) then
         istwfk(ikpt)=1
       end if
     end if
   end do
 end if

 ! If nkpt was to be computed, transfer it from nkpt_computed
 if (nkpt == 0) nkpt = nkpt_computed

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

 call timab(192,2,tsec)

end subroutine inkpts
!!***

!!****f* m_inkpts/inqpt
!! NAME
!! inqpt
!!
!! FUNCTION
!! Initialize the q point for one particular dataset, characterized by jdtset.
!!
!! INPUTS
!!  chksymbreak= if 1, will check whether the k point grid is symmetric, and stop if not.
!!  iout=unit number for echoed output
!!  jdtset=number of the dataset looked for
!!  lenstr=actual length of the string
!!  msym=default maximal number of symmetries
!!  natom=number of atoms
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  spinat(3,1:natom)=spin-magnetization of the atoms
!!  string*(*)=character string containing all the input data. Initialized previously in instrng.
!!  typat(natom)=type for each atom
!!  vacuum(3)=for each direction, 0 if no vacuum, 1 if vacuum
!!  xred(3,natom,nimage) =reduced coordinates of atoms
!!
!! OUTPUT
!!  qptn(3)=reduced coordinates of eventual q point (normalisation is already included)
!!  kptrlatt(3,3)=q-point lattice specification (if kptopt/=0)
!!  wtqc=weigth of the eventual current q point
!!
!! PARENTS
!!      invars1
!!
!! CHILDREN
!!      getkgrid,intagm,metric,symfind,symlatt
!!
!! SOURCE

subroutine inqpt(chksymbreak,iout,jdtset,lenstr,msym,natom,qptn,wtqc,rprimd,spinat,string,typat,vacuum,xred,qptrlatt)

!Arguments ------------------------------------
!scalars
 integer,intent(in)   :: chksymbreak,iout,jdtset,lenstr,msym,natom
 real(dp),intent(inout) :: wtqc
 character(len=*),intent(in) :: string
!arrays
 integer,intent(in) :: typat(natom),vacuum(3)
 real(dp),intent(out) :: qptn(3)
 integer,intent(inout) :: qptrlatt(3,3) !vz_i
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: spinat(3,natom)
 real(dp),intent(in) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: ii,iqpt,iscf_fake,marr,nptsym,nqpt,nqpt_computed,nshiftq,nsym_new,qptopt
 integer :: tread,tread_qptrlatt,tread_ngqpt,use_inversion
 real(dp) :: qptnrm,qptrlen,tolsym,ucvol
 character(len=500) :: msg
!arrays
 integer :: bravais(11), ngqpt(3)
 integer, allocatable :: symafm_new(:), ptsymrel(:,:,:),symrel_new(:,:,:), intarr(:)
 real(dp) :: gmet(3,3),gprimd(3,3),qpt(3),rmet(3,3),shiftq(3,MAX_NSHIFTK)
 real(dp),allocatable :: qpts(:,:),tnons_new(:,:),wtq(:), dprarr(:)

! *************************************************************************

 ! Compute the maximum size of arrays intarr and dprarr (nshiftq is MAX_NSHIFTK at maximum)
 marr=630
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

 ! Find the method to generate the q-points
 qptopt=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'qptopt',tread,'INT')
 if(tread==1)qptopt=intarr(1)

 if(qptopt==0)then
   ! Read qpt and qptnrm
   qpt=zero
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'qpt',tread,'DPR')
   if(tread==1) qpt(1:3)=dprarr(1:3)

   qptnrm=1.0_dp
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'qptnrm',tread,'DPR')

   if(tread==1) qptnrm=dprarr(1)
   if(qptnrm<tol10)then
     write(msg, '(5a)' )&
     'The input variable qptnrm is lower than 1.0d-10,',ch10,&
     'while it must be a positive, non-zero number.   ',ch10,&
     'Action: correct the qptnrm in the input file.'
     MSG_ERROR(msg)
   end if

   qptn(:)=qpt(:)/qptnrm

   ! DBSP: one could want ot define wtq in order to reproduce what is obtained
   ! with ngqpt but without having to do initialize the qgrid (extremly slow in case of large grid > 50x50x50
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'wtq',tread,'DPR')
   if(tread==1) wtqc=dprarr(1)

 else if (qptopt>=1 .and. qptopt<=4) then
   ngqpt(:)=0
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ngqpt',tread_ngqpt,'INT')

   if(tread_ngqpt==1)then
     ngqpt(1:3)=intarr(1:3)
     do ii=1,3
       if(ngqpt(ii)<1)then
         write(msg,  '(a,i0,a,a,a,i0,a,a,a)' ) &
         'The input variable ngqpt(',ii,') must be strictly positive,',ch10,&
         'while it is found to be',ngqpt(ii),'.',ch10,&
         'Action: change it in your input file, or change qptopt.'
         MSG_ERROR(msg)
       end if
     end do
   end if

   call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),'qptrlatt',tread_qptrlatt,'INT')
   if(tread_qptrlatt==1) qptrlatt(:,:)=reshape(intarr(1:9), (/3,3/) )

   if(tread_ngqpt==1 .and. tread_qptrlatt==1)then
     write(msg, '(5a)' ) &
     'The input variables ngqpt and qptrlatt cannot both ',ch10,&
     'be defined in the input file.',ch10,&
     'Action: change one of ngqpt or qptrlatt in your input file.'
     MSG_ERROR(msg)
   else if(tread_ngqpt==1)then
     qptrlatt(:,:)=0
     qptrlatt(1,1)=ngqpt(1)
     qptrlatt(2,2)=ngqpt(2)
     qptrlatt(3,3)=ngqpt(3)
   end if

   nshiftq=1
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nshiftq',tread,'INT')
   if(tread==1)nshiftq=intarr(1)

   if (nshiftq<1 .or. nshiftq>MAX_NSHIFTK) then
     write(msg,  '(a,i0,2a,i0,3a)' )&
     'The only allowed values of nshiftq are between 1 and,',MAX_NSHIFTK,ch10,&
     'while it is found to be',nshiftq,'.',ch10,&
     'Action: change the value of nshiftq in your input file, or change qptopt.'
     MSG_ERROR(msg)
   end if

   shiftq=zero
   call intagm(dprarr,intarr,jdtset,marr,3*nshiftq,string(1:lenstr),'shiftq',tread,'DPR')

   if(tread==1)then
     shiftq(:,1:nshiftq)=reshape( dprarr(1:3*nshiftq), (/3,nshiftq/) )
   else
     if(nshiftq/=1)then
       write(msg,  '(3a,i0,2a)' )&
       'When nshiftq is not equal to 1, shiftq must be defined in the input file.',ch10,&
       'However, shiftq is not defined, while nshiftq=',nshiftq,ch10,&
       'Action: change the value of nshiftq in your input file, or define shiftq.'
       MSG_ERROR(msg)
     end if
   end if

   ! Re-generate symmetry operations from the lattice and atomic coordinates
   ! This is a fundamental difference with respect to the k point generation.
   tolsym=tol8
   ABI_ALLOCATE(ptsymrel,(3,3,msym))
   ABI_ALLOCATE(symafm_new,(msym))
   ABI_ALLOCATE(symrel_new,(3,3,msym))
   ABI_ALLOCATE(tnons_new,(3,msym))
   call symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
   use_inversion=1
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
   call symfind(0,(/zero,zero,zero/),gprimd,0,msym,natom,0,nptsym,nsym_new,0,0,&
    ptsymrel,spinat,symafm_new,symrel_new,tnons_new,tolsym,typat,use_inversion,xred)

   ! Prepare to compute the q-point grid in the ZB or IZB
   iscf_fake=0 ! Do not need the weights

   ! Compute the maximum number of q points
   nqpt=0
   ABI_ALLOCATE(qpts,(3,nqpt))
   ABI_ALLOCATE(wtq,(nqpt))
   call getkgrid(chksymbreak,0,iscf_fake,qpts,qptopt,qptrlatt,qptrlen,&
     msym,nqpt,nqpt_computed,nshiftq,nsym_new,rprimd,&
     shiftq,symafm_new,symrel_new,vacuum,wtq)

   nqpt=nqpt_computed
   ABI_DEALLOCATE(qpts)
   ABI_DEALLOCATE(wtq)

   ! Find the index of the q point within the set of q points that will be generated
   iqpt=0
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'iqpt',tread,'INT')
   if(tread==1)iqpt=intarr(1)

   ! Checks that iqpt is among the computed q points
   if(iqpt<0)then
     write(msg, '(a,i0,3a)' )&
      'The input variable iqpt,',iqpt,' is negative, while it should be 0 or positive.',ch10,&
      'Action: correct iqpt in the input file.'
     MSG_ERROR(msg)
   end if

   if (iqpt > nqpt_computed) then
     write(msg, '(a,i0,3a,i0,7a)' )&
      'The input variable iqpt,',iqpt,' is bigger than the computed number of q-points in the grid,',ch10,&
      'which is ',nqpt,'.',ch10,&
      'The latter has been computed from the input variables qptrlatt, ngqpt, nshiftq,',ch10,&
      'shiftq, as well as qptopt, the symmetries of the lattice, and spinat.',ch10,&
      'Action: correct iqpt in the input file, or correct the computed q-point grid.'
     MSG_ERROR(msg)
   end if

   ! Compute the q-point grid in the BZ or the IBZ
   ABI_ALLOCATE(qpts,(3,nqpt))
   ABI_ALLOCATE(wtq,(nqpt))

   call getkgrid(chksymbreak,iout,iscf_fake,qpts,qptopt,qptrlatt,qptrlen,&
    msym,nqpt,nqpt_computed,nshiftq,nsym_new,rprimd,&
    shiftq,symafm_new,symrel_new,vacuum,wtq)

   ! Transfer to qptn, and deallocate
   qptn(:)=zero
   if(iqpt/=0)then
     qptn(:)=qpts(:,iqpt)
     wtqc = wtq(iqpt)
   end if

   ABI_DEALLOCATE(ptsymrel)
   ABI_DEALLOCATE(symafm_new)
   ABI_DEALLOCATE(symrel_new)
   ABI_DEALLOCATE(tnons_new)
   ABI_DEALLOCATE(qpts)
   ABI_DEALLOCATE(wtq)

 else
   write(msg, '(3a,i0,3a)' ) &
    'The only values of qptopt allowed are smaller than 4.',ch10,&
    'The input value of qptopt is',qptopt,'.',ch10,&
    'Action: change qptopt in your input file.'
   MSG_ERROR(msg)
 end if

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

end subroutine inqpt
!!***

end module m_inkpts
!!***
