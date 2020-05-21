!!****m* ABINIT/m_invars1
!! NAME
!!  m_invars1
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, AR, MKV, FF, MM)
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

module m_invars1

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_atomdata
 use m_dtset
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_fstrings, only : inupper, itoa, endswith, strcat, sjoin, startswith
 use m_geometry, only : mkrdim
 use m_parser,   only : intagm, chkint_ge, ab_dimensions, geo_t, geo_from_abivar_string
 use m_inkpts,   only : inkpts, inqpt
 use m_ingeo,    only : ingeo, invacuum
 use m_symtk,    only : mati3det

#if defined HAVE_GPU_CUDA
 use m_gpu_toolbox
#endif

 implicit none

 private
!!***

 public :: invars0
 public :: invars1
 public :: invars1m
 public :: indefo
!!***

contains
!!***

!!****f* ABINIT/invars0
!! NAME
!! invars0
!!
!! FUNCTION
!! Initialisation phase: prepare the main input subroutine call by
!! reading most of the NO MULTI variables, as well as natom, nimage, and ntypat,
!! needed for allocating some input arrays in abinit, and also useri
!! and userr. The variable usewvl is also read here for later reading
!! of input path for the atomic orbital file (if required).
!!
!! INPUTS
!!  lenstr=actual length of string
!!  ndtset= number of datasets to be read; if 0, no multi-dataset mode
!!  ndtset_alloc=number of datasets, corrected for allocation of at least one data set.
!!  string*(*)=string of characters containing all input variables and data
!!  comm= MPI communicator
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are initialized here:
!!   cpus,jdtset,natom,nimage,npsp,ntypat,useri*,userr*
!!  istatr=repetition rate for status file
!!  istatshft=shift of the repetition rate for status file
!!  msym=maximal value of input msym for all the datasets
!!  mxnatom=maximal value of input natom for all the datasets
!!  mxnimage=maximal value of input nimage for all the datasets
!!  mxntypat=maximal value of input ntypat for all the datasets
!!  npsp=number of pseudopotentials
!!  pseudo_paths(npsp): List of paths to pseudopotential files as read from input file.
!!   List of empty strings if we are legacy "files file" mode. Allocated here, caller should free memory.
!!
!! PARENTS
!!      m_ab7_invars_f90
!!
!! CHILDREN
!!      get_ndevice,intagm
!!
!! SOURCE

subroutine invars0(dtsets, istatr, istatshft, lenstr, msym, mxnatom, mxnimage, mxntypat, ndtset, ndtset_alloc, &
    npsp, pseudo_paths, papiopt, timopt, string, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lenstr,ndtset,ndtset_alloc, comm
 integer,intent(out) :: istatr,istatshft,msym,mxnatom,mxnimage,mxntypat,npsp,papiopt
 integer,intent(inout) :: timopt
 character(len=*),intent(in) :: string
!arrays
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc) !vz_i
 character(len=fnlen),allocatable,intent(out) :: pseudo_paths(:)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,idtset,ii,jdtset,marr,multiplicity,tjdtset,tread,treadh,treadm,tread_pseudos,cnt, tread_geo
 integer :: treads, use_gpu_cuda, ierr
 real(dp) :: cpus
 character(len=500) :: msg
 character(len=fnlen) :: pp_dirpath, shell_var
 character(len=20*fnlen) :: pseudos_string ! DO NOT decrease len
 character(len=len(string)) :: geo_string
 type(geo_t) :: geo
!arrays
 integer,allocatable :: intarr(:), sidx(:)
 real(dp),allocatable :: dprarr(:)

!******************************************************************

 !write(std_out,"(3a)")"invars1 with string:", ch10, trim(string)

 marr=max(9,ndtset_alloc,2)
 ABI_ALLOCATE(dprarr,(marr))
 ABI_ALLOCATE(intarr,(marr))

 ! Set up jdtset
 if (ndtset/=0) then
   ! Default values
   dtsets(0)%jdtset = -1 ! unused value
   dtsets(1:ndtset_alloc)%jdtset=(/ (ii,ii=1,ndtset_alloc) /)

   ! Read explicitly the jdtset array
   call intagm(dprarr,intarr,0,marr,ndtset,string(1:lenstr),'jdtset',tjdtset,'INT')
   if(tjdtset==1) dtsets(1:ndtset)%jdtset=intarr(1:ndtset)

   ! Read the udtset array
   call intagm(dprarr,intarr,0,marr,2,string(1:lenstr),'udtset',tread,'INT')

   ! jdtset and udtset cannot be defined together
   if(tjdtset==1 .and. tread==1)then
     write(msg, '(3a)' )&
     'jdtset and udtset cannot be defined both in the input file.',ch10,&
     'Action: remove one of them from your input file.'
     MSG_ERROR(msg)
   end if

   ! Check values of udtset
   if(tread==1)then
     if(intarr(1)<1 .or. intarr(1)>999)then
       write(msg, '(a,i0,3a)' )&
       'udtset(1) must be between 1 and 999, but it is ',intarr(1),'.',ch10,&
       'Action: change the value of udtset(1) in your input file.'
       MSG_ERROR(msg)
     end if
     if(intarr(2)<1 .or. intarr(2)>9)then
       write(msg, '(a,i0,3a)' )&
       'udtset(2) must be between 1 and 9, but it is ',intarr(2),'.',ch10,&
       'Action: change the value of udtset(2) in your input file.'
       MSG_ERROR(msg)
     end if
     if(intarr(1)*intarr(2) /= ndtset)then
       write(msg, '(3a,i0,3a,i0,a,i0,3a,i0,3a)' )&
       'udtset(1)*udtset(2) must be equal to ndtset,',ch10,&
       'but it is observed that udtset(1) = ',intarr(1),',',ch10,&
       'and udtset(2) = ',intarr(2),' so that their product is ',intarr(1)*intarr(2),',',ch10,&
       'while ndtset is ',ndtset,'.',ch10,&
       'Action: change udtset or ndtset in your input file.'
       MSG_ERROR(msg)
     end if
     idtset=0
     do i1=1,intarr(1)
       do i2=1,intarr(2)
         idtset=idtset+1
         dtsets(idtset)%jdtset=i1*10+i2
       end do
     end do
   end if

   ! Final check on the jdtset values
   do idtset=1,ndtset
     if(dtsets(idtset)%jdtset<1 .or. dtsets(idtset)%jdtset>9999)then
       write(msg, '(3a,i0,a,i0,a,a)' )&
       'The components of jdtset must be between 1 and 9999.',ch10,&
       'However, the input value of the component ',idtset,' of jdtset is ',dtsets(idtset)%jdtset,ch10,&
       'Action: correct jdtset in your input file.'
       MSG_ERROR(msg)
     end if
   end do

 else
   dtsets(1)%jdtset=0
 end if

 papiopt = 0
 call intagm(dprarr,intarr,0,1,1,string(1:lenstr),'papiopt',tread,'INT')
 if(tread==1) papiopt=intarr(1)

 ! Read timopt and pass it to timab
 call intagm(dprarr,intarr,0,1,1,string(1:lenstr),'timopt',tread,'INT')
 if(tread==1) timopt=intarr(1)

 istatr=0
 dtsets(0)%istatr=istatr
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),'istatr',tread,'INT')
 if(tread==1) istatr=intarr(1)
 dtsets(1:)%istatr=istatr

 istatshft=1
 dtsets(0)%istatshft=istatshft
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),'istatshft',tread,'INT')
 if(tread==1) istatshft=intarr(1)
 dtsets(1:)%istatshft=istatshft

 cpus=zero
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),'cpus ',treads,'DPR')
 if(treads==1) cpus=dprarr(1)
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),'cpum ',treadm,'DPR')
 if(treadm==1) cpus=dprarr(1)*60.0_dp
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),'cpuh ',treadh,'DPR')

 if(treadh==1) cpus=dprarr(1)*3600.0_dp
 if(treads+treadm+treadh>1)then
   write(msg, '(5a)' )&
   'More than one input variable is used to defined the CPU time limit.',ch10,&
   'This is not allowed.',ch10,&
   'Action: in the input file, suppress either cpus, cpum or cpuh.'
   MSG_ERROR(msg)
 end if
 dtsets(:)%cpus=cpus

 ! Default for natom, nimage, ntypat, useri and userr
 dtsets(:)%natom=1
 dtsets(:)%nimage=1
 dtsets(:)%ntypat=1 ; dtsets(0)%ntypat=0    ! Will always echo ntypat
 dtsets(:)%macro_uj=0
 dtsets(:)%maxnsym=384
 dtsets(:)%useria=0
 dtsets(:)%userib=0
 dtsets(:)%useric=0
 dtsets(:)%userid=0
 dtsets(:)%userie=0
 dtsets(:)%userra=zero
 dtsets(:)%userrb=zero
 dtsets(:)%userrc=zero
 dtsets(:)%userrd=zero
 dtsets(:)%userre=zero
 dtsets(:)%usewvl = 0
 dtsets(:)%plowan_compute=0

 ! Loop on datasets, to find natom and mxnatom, as well as useri and userr
 do idtset=1,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0

   ! proposal: supercell generation in input string before it is read in
   ! call expand_supercell_input(jdtset, lenstr, string)
   !  find supercell, else exit
   !  determinant = ncells
   !  copy rprim,    acell,    xred,    xcart,    vel,    typat,   to
   !       rprim_uc, acell_uc, xred_uc, xcart_uc, vel_uc, typat_uc
   !     NB: also rprim and angdeg need to be updated in non diagonal case!!!
   !  generate supercell info for each of these copying out with translation vectors etc...
   !  set chkprim to 0
   !  done!

   !  Generate the supercell if supercell_latt is specified and update string
   dtsets(idtset)%supercell_latt(:) = 0
   do ii=1,3
     dtsets(idtset)%supercell_latt(ii) = 1
   end do
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),"supercell_latt",tread,'INT')
   if (tread==1) dtsets(idtset)%supercell_latt(:)=intarr(1:3)
   !This test should be update if in the future we allow non-diagonal supercell
   if (any(dtsets(idtset)%supercell_latt(:) < tol10 )) then
     write(msg, '(5a)' )&
      'supercell_latt must have positive parameters and diagonal part',ch10,&
      'This is not allowed.  ',ch10,&
      'Action: modify supercell_latt in the input file.'
     MSG_ERROR(msg)
   end if
   ! Compute the multiplicity of the supercell
   multiplicity=dtsets(idtset)%supercell_latt(1)  &
&   *dtsets(idtset)%supercell_latt(2)  & 
&   *dtsets(idtset)%supercell_latt(3)  
!  call mati3det(dtsets(idtset)%supercell_latt,multiplicity)

   ! Read natom from string
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natom',tread,'INT')

   ! or get it from the structure variable
   call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'structure', tread_geo, &
               'KEY', key_value=geo_string)

   if (tread_geo /= 0) then
     geo = geo_from_abivar_string(geo_string, comm)
     if (tread /= 0) then
       ABI_CHECK(intarr(1) == geo%natom, "natom from variable and from structure do not agree with each other")
     end if
     intarr(1) = geo%natom
     tread = 1
   end if

   !  Might also initialize natom from XYZ file
   if (tread==0) call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'_natom',tread,'INT')

   if (tread==1) then
     dtsets(idtset)%natom=intarr(1)
   else
     write(msg, '(a,i0,2a)' )&
      'Input natom must be defined, but was absent for dataset ',jdtset,ch10,&
      'Action: check the input file.'
     MSG_ERROR(msg)
   end if

   ! Check that natom is greater than 0
   if (dtsets(idtset)%natom<=0) then
     write(msg, '(a,i0,2a,i0,3a)' )&
      'Input natom must be > 0, but was ',dtsets(idtset)%natom,ch10,&
      'for dataset ',jdtset,'. This is not allowed.',ch10,&
      'Action: check the input file.'
     MSG_ERROR(msg)
   end if

   if(multiplicity > 1)then
     dtsets(idtset)%natom = dtsets(idtset)%natom * multiplicity
   end if

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nimage',tread,'INT')
   if(tread==1) dtsets(idtset)%nimage=intarr(1)

   ! Check that nimage is greater than 0
   if (dtsets(idtset)%nimage<=0) then
     write(msg, '(a,i0,4a)' )&
      'nimage must be > 0, but was ',dtsets(idtset)%nimage,ch10,&
      'This is not allowed.',ch10,&
      'Action: check the input file.'
     MSG_ERROR(msg)
   end if

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ntypat',tread,'INT')
   if (tread==1) dtsets(idtset)%ntypat=intarr(1)

   if (tread_geo /= 0) then
     if (tread == 1) then
       ABI_CHECK(geo%ntypat == dtsets(idtset)%ntypat, "ntypat and geo%ntypat do not agree with each other")
     end if
     dtsets(idtset)%ntypat = geo%ntypat
   end if

   ! Check that ntypat is greater than 0
   if (dtsets(idtset)%ntypat<=0) then
     write(msg, '(a,i0,2a,i0,3a)' )&
      'Input ntypat must be > 0, but was ',dtsets(idtset)%ntypat,ch10,&
      'for dataset ',jdtset,'. This is not allowed.',ch10,&
      'Action: check the input file.'
     MSG_ERROR(msg)
   end if

   ! Read msym from string
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'maxnsym',tread,'INT')
   if(tread==1)dtsets(idtset)%maxnsym=intarr(1)
   !  Check that maxnsym is greater than 1
   if (dtsets(idtset)%maxnsym<1) then
     write(msg, '(a,i0,2a,i0,3a)' )&
      'Input maxnsym must be > 1, but was ',dtsets(idtset)%maxnsym,ch10,&
      'for dataset ',jdtset,'. This is not allowed.',ch10,&
      'Action: check the input file.'
     MSG_ERROR(msg)
   end if

   ! Read plowan_compute
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'plowan_compute',tread,'INT')
   if(tread==1) dtsets(idtset)%plowan_compute=intarr(1)
   

   ! Read user* variables
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'useria',tread,'INT')
   if(tread==1) dtsets(idtset)%useria=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userib',tread,'INT')
   if(tread==1) dtsets(idtset)%userib=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'useric',tread,'INT')
   if(tread==1) dtsets(idtset)%useric=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userid',tread,'INT')
   if(tread==1) dtsets(idtset)%userid=intarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userie',tread,'INT')
   if(tread==1) dtsets(idtset)%userie=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userra',tread,'DPR')
   if(tread==1) dtsets(idtset)%userra=dprarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userrb',tread,'DPR')
   if(tread==1) dtsets(idtset)%userrb=dprarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userrc',tread,'DPR')
   if(tread==1) dtsets(idtset)%userrc=dprarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userrd',tread,'DPR')
   if(tread==1) dtsets(idtset)%userrd=dprarr(1)
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'userre',tread,'DPR')
   if(tread==1) dtsets(idtset)%userre=dprarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'usewvl',tread,'INT')
   if(tread==1) dtsets(idtset)%usewvl=intarr(1)

   call geo%free()
 end do ! idtset

!mxnatom =maxval(dtsets(1:ndtset_alloc)%natom)
!mxntypat =maxval(dtsets(1:ndtset_alloc)%ntypat)
!msym =maxval(dtsets(1:ndtset_alloc)%maxnsym)
!There is a bug in the HP compiler, the following should execute properly
 mxnatom=dtsets(1)%natom ; mxnimage=dtsets(1)%nimage
 mxntypat=dtsets(1)%ntypat ; msym=dtsets(1)%maxnsym
 if(ndtset_alloc>1)then
   do idtset=2,ndtset_alloc
     mxnatom =max(dtsets(idtset)%natom,mxnatom)
     mxnimage=max(dtsets(idtset)%nimage,mxnimage)
     mxntypat=max(dtsets(idtset)%ntypat,mxntypat)
     msym    =max(dtsets(idtset)%maxnsym,msym)
   end do
 end if

 if(mxnimage>1)then
   do idtset=2,ndtset_alloc
     if(mxnatom/=dtsets(idtset)%natom)then
       write(msg,'(5a,i0,a,i0,3a,i0,a)')&
       'When there exist one dataset with more than one image,',ch10,&
       'the number of atoms in each dataset must be the same.',ch10,&
       'However, it has been found that for dataset= ',idtset,ch10,&
       'natom= ',dtsets(idtset)%natom,' differs from the maximum number',ch10,&
       'of atoms, mxnatom= ',mxnatom,&
       'Action: check the input variables natom for different datasets.'
       MSG_ERROR(msg)
     end if
   end do
 end if

 ! Set up npsp
 npsp=mxntypat   ! Default value
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),'npsp',tread,'INT')

 if(tread==1)then
   npsp=intarr(1)
 else
   if(ndtset_alloc>1)then
     do idtset=1,ndtset_alloc
       if(dtsets(idtset)%ntypat/=mxntypat)then
         write(msg, '(5a,i0,a,i0,2a,i0,2a)' )&
          ' When npsp is not defined, the input variable ntypat must be',ch10,&
          ' the same for all datasets. However, it has been found that for',ch10,&
          ' jdtset: ',dtsets(idtset)%jdtset,', ntypat= ',dtsets(idtset)%ntypat,ch10,&
          ' differs from the maximum value of ntypat= ',mxntypat,ch10,&
          ' Action: check the input variables npsp and ntypat.'
         MSG_ERROR(msg)
       end if
       if(dtsets(idtset)%ntypat>npsp)then
         write(msg, '(5a,i0,a,i0,a,i0,2a)' )&
          ' The number of pseudopotentials, npsp, must never be smaller than ntypat.',ch10,&
          ' However, it has been found that for',ch10,&
          ' jdtset: ',dtsets(idtset)%jdtset,', ntypat= ',dtsets(idtset)%ntypat,' and npsp=',npsp,ch10,&
          ' Action: check the input variables npsp and ntypat.'
         MSG_ERROR(msg)
       endif
     end do
   end if
 end if
 dtsets(0)%npsp = mxntypat   ! Default value
 dtsets(1:ndtset_alloc)%npsp = npsp

 ! Read pseudopotential directory and pseudo paths from input.
 ! Remember that in "files file mode", this info is passed through the files file so these variables are optional
 pp_dirpath = ""
 call intagm(dprarr, intarr, 0, marr, 1, string(1:lenstr), 'pp_dirpath', tread, 'KEY', key_value=pp_dirpath)
 if (tread == 1) then
   if (pp_dirpath(1:1) == "$") then
     shell_var = pp_dirpath(2:)
     call get_environment_variable(shell_var, pp_dirpath, status=ierr)
     if (ierr == -1) MSG_ERROR(sjoin(shell_var, "is present but string too short for the environment variable"))
     if (ierr == +1) MSG_ERROR(sjoin(shell_var, "does not exist"))
     if (ierr == +2) MSG_ERROR(sjoin(shell_var, "used in input file but processor does not support environment variables"))
     call wrtout(std_out, sjoin(shell_var, "found in env. Assuming pseudos located in:",  pp_dirpath))
   end if
   if (.not. endswith(pp_dirpath, "/")) pp_dirpath = strcat(pp_dirpath, "/")
 end if

 ! String must be large enough to contain ntypat filepaths.
 pseudos_string = ""
 call intagm(dprarr, intarr, 0, marr, 1, string(1:lenstr), "pseudos", tread_pseudos, 'KEY', key_value=pseudos_string)

 ABI_MALLOC(pseudo_paths, (npsp))
 pseudo_paths = ""

 if (tread_pseudos == 1) then
   ! Split pseudos_string using comma and transfer results to pseudos_paths
   ! Make sure string lenght is large enough and input string is consistent with npsp
   ! Lot of checks must be done here!
   !print *, "pseudos_string: ", trim(pseudos_string)
   ABI_ICALLOC(sidx, (npsp + 1))
   sidx(1) = 1; sidx(npsp + 1) = len(pseudos_string)
   cnt = 1
   do ii=1,len(pseudos_string)
     if (pseudos_string(ii:ii) == ",") then
       pseudos_string(ii:ii) = " "
       cnt = cnt + 1
       sidx(cnt) = ii
       ABI_CHECK(cnt <= npsp, "Too many commas in pseudos string!")
     end if
   end do
   if (cnt /= npsp) then
     MSG_ERROR(sjoin("Not enough pseudos in input `pseudos` string, expecting npsp:", itoa(npsp)))
   end if

   do ii=1,npsp
     i1 = sidx(ii)
     i2 = sidx(ii + 1)
     cnt = len(adjustl(trim(pseudos_string(i1:i2))))
     ABI_CHECK(cnt <= fnlen, "pseudo path too small, increase fnlen")
     pseudo_paths(ii) = adjustl(trim(pseudos_string(i1:i2)))
     if (len_trim(pp_dirpath) > 0) then
       if (len_trim(pp_dirpath) + len_trim(pseudo_paths(ii)) > fnlen) then
         MSG_ERROR(sjoin("String of len fnlen:", itoa(fnlen), " too small to contain full pseudo path"))
       end if
       pseudo_paths(ii) = strcat(pp_dirpath, pseudo_paths(ii))
     end if
   end do
   ABI_FREE(sidx)
   !print *, "pp_dirpath: ", trim(pp_dirpath), "pseudos: ", trim(pseudos_string)
 end if

 ! KGB parallelism information (needed at this stage)
 dtsets(:)%paral_kgb=0
 do idtset=1,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'paral_kgb',tread,'INT')
   if(tread==1)dtsets(idtset)%paral_kgb=intarr(1)

   if (dtsets(idtset)%paral_kgb<0 .or. dtsets(idtset)%paral_kgb>1) then
     write(msg,'(a,i0,2a,i0,3a)')&
      'Input paral_kgb must be 0 or 1, but was ',dtsets(idtset)%paral_kgb,ch10,&
      'for dataset ',jdtset,'. This is not allowed.',ch10,&
      'Action: check the input file.'
     MSG_ERROR(msg)
   end if


 end do

!GPU information
 use_gpu_cuda=0
 dtsets(:)%use_gpu_cuda=0
#if defined HAVE_GPU_CUDA && defined HAVE_GPU_CUDA_DP
 call Get_ndevice(ii)
 if (ii>0) then
   do i1=1,ndtset_alloc
     dtsets(i1)%use_gpu_cuda=-1
   end do
 end if
#endif
 do idtset=1,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'use_gpu_cuda',tread,'INT')
   if(tread==1)dtsets(idtset)%use_gpu_cuda=intarr(1)
   if (dtsets(idtset)%use_gpu_cuda==1) use_gpu_cuda=1
 end do
 if (use_gpu_cuda==1) then
#if defined HAVE_GPU_CUDA && defined HAVE_GPU_CUDA_DP
   if (ii<=0) then
     write(msg,'(5a)')&
&     'Input variables use_gpu_cuda is on',ch10,&
&     'but no available GPU device has been detected !',ch10,&
&     'Action: change the input variable use_gpu_cuda.'
     MSG_ERROR(msg)
   end if
#else
   write(msg,'(7a)')&
&   'Input variables use_gpu_cuda is on but abinit hasn''t been built',ch10,&
&   'with (double precision) gpu mode enabled !',ch10,&
&   'Action: change the input variable use_gpu_cuda',ch10,&
&   '        or re-compile ABINIT with double-precision Cuda enabled.'
   MSG_ERROR(msg)
#endif
 end if

 ABI_DEALLOCATE(dprarr)
 ABI_DEALLOCATE(intarr)

 ! We allocate the internal array, depending on the computed values.
 ! WARNING: do not forget to deallocate these arrays in the routine dtset_free
 ! (should make a separate subroutine for allocating/deallocating these records)
 do idtset=0,ndtset_alloc
   ABI_ALLOCATE(dtsets(idtset)%acell_orig,(3,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%algalch,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%amu_orig,(mxntypat,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%chrgat,(mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%constraint_kind,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%corecs,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%densty,(mxntypat,4))
   ABI_ALLOCATE(dtsets(idtset)%dynimage,(mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%iatfix,(3,mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%f4of2_sla,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%f6of2_sla,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%jpawu,(mxntypat,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%kberry,(3,20))
   ABI_ALLOCATE(dtsets(idtset)%lexexch,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%ldaminushalf,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%lpawu,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%mixalch_orig,(npsp,mxntypat,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%mixesimgf,(mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%nucdipmom,(3,mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%pimass,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%ptcharge,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%prtatlist,(mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%quadmom,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%ratsph,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%rprim_orig,(3,3,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%rprimd_orig,(3,3,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%so_psp,(npsp))
   ABI_ALLOCATE(dtsets(idtset)%spinat,(3,mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%shiftk,(3,MAX_NSHIFTK))
   ABI_ALLOCATE(dtsets(idtset)%typat,(mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%upawu,(mxntypat,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%plowan_iatom,(mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%plowan_it,(100*3))
   ABI_ALLOCATE(dtsets(idtset)%plowan_nbl,(mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%plowan_lcalc,(12*mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%plowan_projcalc,(12*mxnatom))
   ABI_ALLOCATE(dtsets(idtset)%vel_orig,(3,mxnatom,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%vel_cell_orig,(3,3,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%xred_orig,(3,mxnatom,mxnimage))
   ABI_ALLOCATE(dtsets(idtset)%ziontypat,(mxntypat))
   ABI_ALLOCATE(dtsets(idtset)%znucl,(npsp))
 end do

!DEBUG
!write(std_out,*)' invars0 : nimage, mxnimage = ',dtsets(:)%nimage, mxnimage
!write(std_out,*)' invars0 : natom = ',dtsets(:)%natom
!write(std_out,*)' invars0 : mxnatom = ',mxnatom
!ENDDEBUG

end subroutine invars0
!!***

!!****f* ABINIT/invars1m
!! NAME
!! invars1m
!!
!! FUNCTION
!! Initialisation phase: prepare the main input subroutine call by
!! reading all the NO MULTI variables, as well as the dimensions
!! needed for allocating the input arrays in abinit.
!!
!! INPUTS
!!  iout=unit number of output file
!!  lenstr=actual length of string
!!  msym=default maximal number of symmetries
!!  mxnatom=maximal value of input natom for all the datasets
!!  mxnimage=maximal value of input nimage for all the datasets
!!  ndtset= number of datasets to be read; if 0, no multi-dataset mode
!!  ndtset_alloc=number of datasets, corrected for allocation of at least one data set.
!!  npsp= number of pseudopotential files
!!  string*(*)=string of characters containing all input variables and data
!!  zionpsp(npsp)= valence charge over all psps
!!  comm=MPI communicator
!!
!! OUTPUT
!!  dmatpuflag=flag controlling the use of an initial density matrix in PAW+U (max. value over datasets)
!!  mband_upper_(0:ndtset_alloc)=list of mband_upper values
!!
!! SIDE EFFECTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are initialized here (see invars1.f for more details on the initialized records)
!!  mx<ab_dimensions>=datatype storing the maximal dimensions. Partly initialized in input.
!!
!! PARENTS
!!
!! CHILDREN
!!      indefo1,invars1
!!
!! SOURCE

subroutine invars1m(dmatpuflag, dtsets, iout, lenstr, mband_upper_, mx,&
& msym, ndtset, ndtset_alloc, string, npsp, zionpsp, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,lenstr,msym,ndtset,ndtset_alloc,npsp, comm
 integer,intent(out) :: dmatpuflag
 character(len=*),intent(inout) :: string
 type(ab_dimensions),intent(inout) :: mx
!arrays
 integer,intent(out) :: mband_upper_(0:ndtset_alloc)
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
 real(dp),intent(in) :: zionpsp(npsp)

!Local variables-------------------------------
!scalars
 integer :: idtset,ii,jdtset,lpawu,mband_upper,iatom,nat,nsp
!arrays
 integer,allocatable :: symafm_(:,:),symrel_(:,:,:,:)
 integer,allocatable :: symafm(:),symrel(:,:,:)
 real(dp),allocatable :: tnons_(:,:,:),tnons(:,:)

!******************************************************************

 ! Here, allocation of the arrays that depend on msym.
 ABI_ALLOCATE(symrel_,(3,3,msym,0:ndtset_alloc))
 ABI_ALLOCATE(symafm_,(msym,0:ndtset_alloc))
 ABI_ALLOCATE(tnons_,(3,msym,0:ndtset_alloc))
 ABI_ALLOCATE(symafm,(msym))
 ABI_ALLOCATE(symrel,(3,3,msym))
 ABI_ALLOCATE(tnons,(3,msym))

 ! Set up default values (note that the default acell, amu mkmem, mkmem1,mkqmem, and nkpt must be overcome
 do idtset=0,ndtset_alloc
   call indefo1(dtsets(idtset))
 end do

 ! natom and nimage are already initialized in invars0
 dtsets(0)%natom=-1
 dtsets(0)%nimage=1

!Initialization for parallelization data has changed
!these lines aim to keep old original default values
 dtsets(0)%npimage=1
 dtsets(0)%npkpt=1
 dtsets(0)%npspinor=1
 dtsets(0)%npfft=1
 dtsets(0)%npband=1
 dtsets(0)%bandpp=1

 symafm_(:,0)=1
 symrel_(:,:,:,0)=0
 symrel_(1,1,:,0)=1 ; symrel_(2,2,:,0)=1 ; symrel_(3,3,:,0)=1
 tnons_(:,:,0)=0.0_dp

 ! Loop on datasets
 do idtset=1,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
   write(std_out,'(2a)') ch10,'======================================================= '
   write(std_out,'(a,i0)') 'invars1m : enter jdtset= ',jdtset

   ! Input default values
   dtsets(idtset)%bravais(:)=0
   symafm(:)=symafm_(:,0)
   symrel(:,:,:)=symrel_(:,:,:,0)
   tnons(:,:)=tnons_(:,:,0)

   call invars1(dtsets(idtset)%bravais,dtsets(idtset),iout,jdtset,lenstr,&
&   mband_upper,msym,npsp,string,symafm,symrel,tnons,zionpsp, comm)

   mband_upper_ (idtset)=mband_upper
   symafm_(:,idtset)=symafm(:)
   symrel_(:,:,:,idtset)=symrel(:,:,:)
   tnons_(:,:,idtset)=tnons(:,:)
 end do

 mx%mband_upper = maxval(mband_upper_ (1:ndtset_alloc))

 dmatpuflag = 0; mx%natpawu = 0; mx%lpawu = 0
 mx%natsph = dtsets(1)%natsph
 mx%natsph_extra = dtsets(1)%natsph_extra
 mx%natvshift = dtsets(1)%natvshift
 mx%nconeq = dtsets(1)%nconeq
 mx%n_efmas_dirs=0
 mx%ga_n_rules = dtsets(1)%ga_n_rules
 mx%gw_nqlwl = dtsets(1)%gw_nqlwl
 mx%nimfrqs = 0
 mx%nfreqsp = 0
 mx%n_projection_frequencies = 0
 mx%nkpt  = dtsets(1)%nkpt
 mx%nkptgw = dtsets(1)%nkptgw
 mx%nkpthf = dtsets(1)%nkpthf
 mx%nnos  = dtsets(1)%nnos
 mx%nqptdm = dtsets(1)%nqptdm
 mx%nspinor = dtsets(1)%nspinor
 mx%nsppol = dtsets(1)%nsppol
 mx%ntypat = dtsets(1)%ntypat
 mx%nzchempot = dtsets(1)%nzchempot
 mx%nberry = 20   ! This is presently a fixed value. Should be changed.

 ! Get MAX dimension over datasets
 do ii=1,ndtset_alloc
   mx%natsph = max(dtsets(ii)%natsph, mx%natsph)
   mx%natsph_extra=max(dtsets(ii)%natsph_extra, mx%natsph_extra)
   mx%nconeq=max(dtsets(ii)%nconeq, mx%nconeq)
   mx%n_efmas_dirs = max(dtsets(ii)%efmas_n_dirs, mx%n_efmas_dirs)
   mx%ga_n_rules = max(dtsets(ii)%ga_n_rules,mx%ga_n_rules)
   mx%gw_nqlwl = max(dtsets(ii)%gw_nqlwl,mx%gw_nqlwl)
   mx%nimfrqs = max(dtsets(ii)%cd_customnimfrqs, mx%nimfrqs)
   mx%nfreqsp = max(dtsets(ii)%gw_customnfreqsp, mx%nfreqsp)
   mx%n_projection_frequencies = max(dtsets(ii)%gwls_n_proj_freq, mx%n_projection_frequencies)
   mx%nkpt  = max(dtsets(ii)%nkpt, mx%nkpt)
   mx%nkptgw = max(dtsets(ii)%nkptgw, mx%nkptgw)
   mx%nkpthf = max(dtsets(ii)%nkpthf, mx%nkpthf)
   mx%nnos  = max(dtsets(ii)%nnos, mx%nnos)
   mx%nqptdm = max(dtsets(ii)%nqptdm, mx%nqptdm)
   mx%nspinor = max(dtsets(ii)%nspinor, mx%nspinor)
   mx%nsppol = max(dtsets(ii)%nsppol, mx%nsppol)
   mx%ntypat = max(dtsets(ii)%ntypat, mx%ntypat)
   mx%nzchempot = max(dtsets(ii)%nzchempot, mx%nzchempot)
   if (dtsets(ii)%usepawu/=0) then
     if (dtsets(ii)%usepawu>0.and.dtsets(ii)%usedmatpu/=0) dmatpuflag=1
     lpawu=maxval(dtsets(ii)%lpawu(:))
     mx%lpawu=max(lpawu,mx%lpawu)
     !dtsets(ii)%natpawu=count(dtsets(ii)%lpawu(dtsets(ii)%typat((/(i1,i1=1,dtsets(ii)%natom)/)))/=-1)
     ! Old fashon way that should do fine
     dtsets(ii)%natpawu = 0
     do iatom=1, dtsets(ii)%natom
       if (dtsets(ii)%lpawu(dtsets(ii)%typat(iatom)) /= -1 ) dtsets(ii)%natpawu = dtsets(ii)%natpawu + 1
     end do
     mx%natpawu = max(dtsets(ii)%natpawu, mx%natpawu)
     if (dtsets(ii)%macro_uj/=0) dtsets(ii)%natvshift=lpawu*2+1
   end if
   mx%natvshift = max(dtsets(ii)%natvshift, mx%natvshift)
 end do

!mx%nsym=maxval(dtsets(1:ndtset_alloc)%nsym) ! This might not work properly with HP compiler
 mx%nsym=dtsets(1)%nsym
 do idtset=1,ndtset_alloc
   mx%nsym = max(dtsets(idtset)%nsym, mx%nsym)
 end do

 do idtset=0,ndtset_alloc
   ABI_ALLOCATE(dtsets(idtset)%atvshift, (mx%natvshift, mx%nsppol, mx%natom))
   ABI_ALLOCATE(dtsets(idtset)%bs_loband,(mx%nsppol))
   ABI_ALLOCATE(dtsets(idtset)%bdgw,(2, mx%nkptgw, mx%nsppol))
   ABI_ALLOCATE(dtsets(idtset)%cd_imfrqs,(mx%nimfrqs))
   ABI_ALLOCATE(dtsets(idtset)%chempot,(3, mx%nzchempot, mx%ntypat))
   nsp = max(mx%nsppol, mx%nspinor); nat = mx%natpawu*dmatpuflag
   ABI_ALLOCATE(dtsets(idtset)%dmatpawu,(2*mx%lpawu+1,2*mx%lpawu+1,nsp,nat, mx%nimage))
   ABI_ALLOCATE(dtsets(idtset)%efmas_bands,(2, mx%nkpt))
   ABI_ALLOCATE(dtsets(idtset)%efmas_dirs,(3, mx%n_efmas_dirs))
   ABI_ALLOCATE(dtsets(idtset)%gw_freqsp, (mx%nfreqsp))
   ABI_ALLOCATE(dtsets(idtset)%gwls_list_proj_freq, (mx%n_projection_frequencies))
   ABI_ALLOCATE(dtsets(idtset)%gw_qlwl,(3,mx%gw_nqlwl))
   ABI_ALLOCATE(dtsets(idtset)%kpt,(3, mx%nkpt))
   ABI_ALLOCATE(dtsets(idtset)%kptgw,(3, mx%nkptgw))
   ABI_ALLOCATE(dtsets(idtset)%kptns,(3, mx%nkpt))
   ABI_ALLOCATE(dtsets(idtset)%kptns_hf,(3, mx%nkpthf))
   ABI_ALLOCATE(dtsets(idtset)%iatsph,(mx%natsph))
   ABI_ALLOCATE(dtsets(idtset)%istwfk, (mx%nkpt))
   ABI_ALLOCATE(dtsets(idtset)%nband, (mx%nkpt*mx%nsppol))
   ABI_ALLOCATE(dtsets(idtset)%occ_orig,(mx%mband_upper*mx%nkpt*mx%nsppol, mx%nimage))
   ABI_ALLOCATE(dtsets(idtset)%qmass, (mx%nnos))
   ABI_ALLOCATE(dtsets(idtset)%qptdm,(3, mx%nqptdm))
   ABI_ALLOCATE(dtsets(idtset)%symafm, (mx%nsym))
   ABI_ALLOCATE(dtsets(idtset)%symrel,(3,3,mx%nsym))
   ABI_ALLOCATE(dtsets(idtset)%tnons,(3,mx%nsym))
   ABI_ALLOCATE(dtsets(idtset)%wtatcon,(3,mx%natom, mx%nconeq))
   ABI_ALLOCATE(dtsets(idtset)%wtk, (mx%nkpt))
   ABI_ALLOCATE(dtsets(idtset)%xredsph_extra,(3, mx%natsph_extra))
   dtsets(idtset)%symrel(:,:,:)=symrel_(:,:,1:mx%nsym,idtset)
   dtsets(idtset)%symafm(:)    =symafm_(1:mx%nsym,idtset)
   dtsets(idtset)%tnons (:,:)  =tnons_ (:,1:mx%nsym,idtset)
 end do

 ABI_DEALLOCATE(symafm_)
 ABI_DEALLOCATE(symrel_)
 ABI_DEALLOCATE(tnons_)
 ABI_DEALLOCATE(symafm)
 ABI_DEALLOCATE(symrel)
 ABI_DEALLOCATE(tnons)

end subroutine invars1m
!!***

!!****f* ABINIT/indefo1
!! NAME
!! indefo1
!!
!! FUNCTION
!! Initialisation phase : defaults values for a first batch of input variables
!! (especially dimensions, needed to allocate other parts of dtsets, as well
!!  as other input variables whose existence is needed for other initialisations to proceed).
!!
!! INPUTS
!!
!! OUTPUT
!!  dtset=<type datafiles_type>contains all input variables for one dataset,
!!   some of which are given a default value here.
!!
!! PARENTS
!!      invars1m
!!
!! CHILDREN
!!
!! SOURCE

subroutine indefo1(dtset)

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(inout) :: dtset

!Local variables -------------------------------
!scalars
 !integer :: ii

!******************************************************************

!Set up default values. All variables to be output in outvars.f
!should have a default, even if a nonsensible one can be chosen to garantee print in that routine.

 DBG_ENTER("COLL")

!Use alphabetic order

!A
 dtset%acell_orig(:,:)=zero
 dtset%algalch(:)=1
 dtset%amu_orig(:,:)=-one
 dtset%autoparal=0
!B
 dtset%bandpp=1
 dtset%berryopt=0
 dtset%berrysav=0
 dtset%bfield(:)=zero
!C
 dtset%cd_customnimfrqs=0
 dtset%chkprim=1
 dtset%chrgat(:)=zero
 dtset%constraint_kind(:)=0
!D
 dtset%densty(:,:)=zero
 dtset%dfield(:)=zero    !!HONG
 dtset%dynimage(:)=1
!E
 dtset%efield(:)=zero
 dtset%efmas_calc_dirs=0
 dtset%efmas_n_dirs=0
!F
!G
 dtset%ga_n_rules=1
 dtset%gw_customnfreqsp=0
 dtset%gw_nqlwl=0
 dtset%gwls_n_proj_freq=0
!I
 dtset%iatfix(:,:)=0
 dtset%icoulomb=0
 dtset%imgmov=0
!J
 dtset%jellslab=0
 dtset%jfielddir(:)=0
!K
 dtset%kptopt=0
!L
 dtset%lexexch(:)=-1
 dtset%ldaminushalf(:)=0
 dtset%lpawu(:)=-1
!M
 dtset%maxestep=0.005d0
 dtset%mixalch_orig(:,:,:)=zero
 dtset%mkmem=-1
 dtset%mkqmem=-1
 dtset%mk1mem=-1
!N
 dtset%natpawu=0
 dtset%natsph=0
 dtset%natsph_extra=0
 dtset%natvshift=0
 dtset%nconeq=0
 dtset%ndynimage=1
 dtset%nkpt=-1
 dtset%nkptgw=0
 dtset%nkpthf=0
 dtset%nnos=0
 dtset%npband=1
 dtset%npfft=1
 dtset%nphf=1
 dtset%npimage=1
 dtset%npkpt=1
 dtset%nppert=1
 dtset%npspalch=0
 dtset%npspinor=1
 dtset%np_slk=1000000
 dtset%nqptdm=0
 dtset%nspden=1
 dtset%nspinor=1
 dtset%nsppol=1
 dtset%nsym=0     ! Actually, this default value is not used : it is to be reimposed before each call to ingeo in invars1
 dtset%ntimimage=1
 dtset%ntypalch=0
 dtset%ntyppure=-1
 dtset%nucdipmom(:,:)=zero
 dtset%nzchempot=0
!O
 dtset%optdriver=0
!P
 dtset%paral_rf=0
!dtset%paral_kgb ! Is even initialized earlier.
 dtset%pawspnorb=0  ! will be changed to 1 as soon as usepaw==1 and nspinor==2
 dtset%pimass(:)=-one
!Q
 dtset%qptn=zero
!R
 dtset%red_efield(:)=zero
 dtset%red_dfield(:)=zero
 dtset%red_efieldbar(:)=zero
 dtset%rprim_orig(:,:,:)=zero
 dtset%rprim_orig(1,1,:)=one
 dtset%rprim_orig(2,2,:)=one
 dtset%rprim_orig(3,3,:)=one
!S
 dtset%slabzbeg=zero
 dtset%slabzend=zero
 dtset%so_psp(:)=1
 dtset%spinat(:,:)=zero
 dtset%symmorphi=1
!T
 dtset%tfkinfunc=0
 dtset%typat(:)=0  ! This init is important because dimension of typat is mx%natom (and not natom).
!U
 dtset%usedmatpu=0
 dtset%usedmft=0
 dtset%useexexch=0
 dtset%usepawu=0
 dtset%usepotzero=0
 dtset%use_slk=0
!V
 dtset%vel_orig(:,:,:)=zero
 dtset%vel_cell_orig(:,:,:)=zero
!W
 dtset%wtq=zero
 if (dtset%usepaw==0) dtset%wfoptalg=0
 if (dtset%usepaw/=0) dtset%wfoptalg=10
 if (dtset%optdriver==RUNL_GSTATE.and.dtset%paral_kgb>0) dtset%wfoptalg=14
 dtset%wvl_bigdft_comp=1

!X
 dtset%xred_orig(:,:,:)=zero
!Y
!Z
 dtset%zeemanfield(:)=zero

 DBG_EXIT("COLL")

end subroutine indefo1
!!***

!!****f* ABINIT/invars1
!! NAME
!! invars1
!!
!! FUNCTION
!! Initialize the dimensions needed to allocate the input arrays
!! for one dataset characterized by jdtset, by taking from string the necessary data.
!! Perform some preliminary checks and echo these dimensions.
!!
!! INPUTS
!!  iout=unit number of output file
!!  jdtset=number of the dataset looked for
!!  lenstr=actual length of string
!!  msym=default maximal number of symmetries
!!  npsp1= number of pseudopotential files
!!  zionpsp(npsp1)= valence charge over all psps
!!  comm= MPI communicator
!!
!! OUTPUT
!!  mband_upper=estimation of the maximum number of bands for any k-point
!!
!! SIDE EFFECTS
!! Input/Output (the default value is given in the calling routine)
!!  dtset=<type datafiles_type>contains all input variables,
!!   some of which are initialized here, while other were already
!!   initialized, while some others will still be initialized later.
!!   The list of records of dtset initialized in the present routine is:
!!
!!       acell_orig,chrgat,densty,iatfix,kptopt,kptrlatt,
!!       mkmem,mkqmem,mk1mem,natsph,natvshift,nconeq,nkpt,nkptgw,nkpthf,
!!       nqptdm,nshiftk,nucdipmom,nzchempot,optdriver,
!!       rprim_orig,rprimd_orig,shiftk,
!!       spgroup,spinat,typat,vel_orig,vel_cell_orig,xred_orig
!!
!!  bravais(11)=characteristics of Bravais lattice (see symlatt.F90)
!!  symafm(1:msym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,1:msym)=symmetry operations in real space in terms of primitive translations
!!  tnons(3,1:msym)=nonsymmorphic translations for symmetry operations
!!  string*(*)=string of characters containing all input variables and data
!!
!! NOTES
!! Must set up the geometry of the system, needed to compute k point grids in an automatic fashion.
!! Treat separately mband_upper, since fband, charge and zionpsp must be known for being able to initialize it.
!!
!! Defaults are provided in the calling routine.
!! Defaults are also provided here for the following variables:
!!
!!      mband_upper, occopt, fband, charge
!!
!! They should be kept consistent with defaults of the same variables provided to the invars routines.
!!
!! PARENTS
!!      invars1m
!!
!! CHILDREN
!!      atomdata_from_znucl,chkint_ge,ingeo,inkpts,inqpt,intagm,inupper
!!      invacuum,mkrdim,wrtout
!!
!! SOURCE

subroutine invars1(bravais,dtset,iout,jdtset,lenstr,mband_upper,msym,npsp1,&
& string,symafm,symrel,tnons,zionpsp, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,jdtset,lenstr,msym,npsp1, comm
 integer,intent(out) :: mband_upper
 character(len=*),intent(inout) :: string
 type(dataset_type),intent(inout) :: dtset
!arrays
 integer,intent(inout) :: bravais(11),symafm(msym),symrel(3,3,msym)
 real(dp),intent(inout) :: tnons(3,msym)
 real(dp),intent(in) :: zionpsp(npsp1)

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: chksymbreak,found,ierr,iatom,ii,ikpt,iimage,index_blank,index_lower, tread_geo
 integer :: index_typsymb,index_upper,ipsp,iscf,intimage,itypat,leave,marr
 integer :: natom,nkpt,nkpthf,npsp,npspalch, ncid
 integer :: nqpt,nspinor,nsppol,ntypat,ntypalch,ntyppure,occopt,response
 integer :: rfddk,rfelfd,rfphon,rfstrs,rfuser,rf2_dkdk,rf2_dkde,rfmagn
 integer :: tfband,tnband,tread,tread_alt, my_rank, nprocs
 real(dp) :: charge,fband,kptnrm,kptrlen,zelect,zval
 character(len=1) :: blank=' ',string1
 character(len=2) :: string2,symbol
 character(len=500) :: msg
 type(atomdata_t) :: atom
!arrays
 integer :: cond_values(4),vacuum(3)
 integer,allocatable :: iatfix(:,:),intarr(:),istwfk(:),nband(:),typat(:)
 real(dp) :: acell(3),rprim(3,3)
!real(dp) :: field(3)
 real(dp),allocatable :: amu(:),chrgat(:),dprarr(:),kpt(:,:),kpthf(:,:),mixalch(:,:),nucdipmom(:,:)
 real(dp),allocatable :: ratsph(:),reaalloc(:),spinat(:,:)
 real(dp),allocatable :: vel(:,:),vel_cell(:,:),wtk(:),xred(:,:),znucl(:)
 character(len=32) :: cond_string(4)
 character(len=fnlen) :: key_value
 character(len=len(string)) :: geo_string
 type(geo_t) :: geo

!************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 ! This counter is incremented when we find a non-critical error.
 ! The code outputs a warning and stops at end.
 leave = 0

 ! Some initialisations
 ierr=0
 cond_string(1:4)=' '
 cond_values(1:4)=(/0,0,0,0/)

 ! Read parameters
 marr=dtset%npsp;if (dtset%npsp<3) marr=3
 marr=max(marr,dtset%nimage)
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

!---------------------------------------------------------------------------

 rfddk=0; rfelfd=0; rfphon=0; rfmagn=0; rfstrs=0; rfuser=0; rf2_dkdk=0; rf2_dkde=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfddk',tread,'INT')
 if(tread==1) rfddk=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfelfd',tread,'INT')
 if(tread==1) rfelfd=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfmagn',tread,'INT')
 if(tread==1) rfmagn=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfphon',tread,'INT')
 if(tread==1) rfphon=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfstrs',tread,'INT')
 if(tread==1) rfstrs=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfuser',tread,'INT')
 if(tread==1) rfuser=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf2_dkdk',tread,'INT')
 if(tread==1) rf2_dkdk=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf2_dkde',tread,'INT')
 if(tread==1) rf2_dkde=intarr(1)

 response=0
 if(rfddk/=0.or.rf2_dkdk/=0.or.rf2_dkde/=0.or.rfelfd/=0.or.rfphon/=0.or.rfstrs/=0.or.rfuser/=0.or.rfmagn/=0)response=1

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'optdriver',tread,'INT')
 if (tread==1) then
   dtset%optdriver=intarr(1)
 else
   ! If optdriver was not read, while response=1, set optdriver to 1
   if(response==1)dtset%optdriver=1
 end if

!---------------------------------------------------------------------------
!For now, waiting express parallelisation for recursion
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tfkinfunc',tread,'INT')
 if(tread==1) dtset%tfkinfunc=intarr(1)

!---------------------------------------------------------------------------
! wvl_bigdft_comp, done here since default values of nline, nwfshist and iscf depend on its value (see indefo)
 if(dtset%usewvl==1) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'wvl_bigdft_comp',tread,'INT')
   if(tread==1) dtset%wvl_bigdft_comp=intarr(1)
 end if

!---------------------------------------------------------------------------

 natom=dtset%natom
 npsp=dtset%npsp
 ntypat=dtset%ntypat

 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'structure', tread_geo, &
             'KEY', key_value=geo_string)

 if (tread_geo == 0) then
   ! No default value for znucl
   call intagm(dprarr,intarr,jdtset,marr,dtset%npsp,string(1:lenstr),'znucl',tread,'DPR')
   if(tread==1) dtset%znucl(1:dtset%npsp)=dprarr(1:dtset%npsp)

   if(tread/=1)then
     write(msg, '(3a)' )&
     'The array znucl MUST be initialized in the input file while this is not done.',ch10,&
     'Action: initialize znucl in your input file.'
     MSG_ERROR(msg)
   end if

 else
   call wrtout(std_out, sjoin(" Initializing lattice and positions from:", geo_string))
   geo = geo_from_abivar_string(geo_string, comm)
   dtset%znucl(1:dtset%ntypat) = geo%znucl
   call geo%free()
 end if

 ! The default for ratsph has already been initialized
 call intagm(dprarr,intarr,jdtset,marr,dtset%ntypat,string(1:lenstr),'ratsph',tread,'LEN')
 if(tread==1)then
   do ii=1,dtset%ntypat
     dtset%ratsph(ii)=dprarr(ii)
   end do
 end if
 ABI_ALLOCATE(ratsph,(dtset%ntypat))
 do ii=1,dtset%ntypat
   ratsph(ii)=dtset%ratsph(ii)
 end do

!Special treatment of _TYPAX (from a XYZ file), taking into account
!the fact that znucl does NOT depend on the dataset
!Examine all occurences of '_TYPAX'

 do
   index_typsymb=index(string(1:lenstr),'_TYPAX')
   if(index_typsymb==0)exit
!  Replace '_TYPAX' by '_TYPAT'
   string(index_typsymb:index_typsymb+5)='_TYPAT'
   index_upper=index_typsymb+5
!  Must start from the first blank after the tag (including possible dtset_char)
   index_upper=index(string(index_upper:lenstr),blank)+index_upper-1
   index_lower=index_upper

!  Examine all atoms (the end of the symbol string is delimited by a XX )
   do
     index_blank=index(string(index_upper:lenstr),blank)+index_upper-1
     string2=string(index_blank+1:index_blank+2)
     if(string2=="XX")exit
     found=0
!    Find the matching symbol
     do ipsp=1,dtset%npsp
       call atomdata_from_znucl(atom,dtset%znucl(ipsp))
       symbol = atom%symbol
       call inupper(symbol)
       call inupper(string2)
!      write(std_out,'(a)')' invars1 : before test, trim(adjustl(symbol)),trim(adjustl(string2))'
!      write(std_out,'(5a)' )'"',trim(adjustl(symbol)),'","',trim(adjustl(string2)),'"'
       if(trim(adjustl(symbol))==trim(adjustl(string2)))then
         found=1
         index_upper=index_blank+1
         ! Cannot deal properly with more that 9 psps
         if(ipsp>=10)then
           MSG_ERROR('Need to use a pseudopotential with number larger than 9. Not allowed yet.')
         end if

         ! write(std_out,*)' invars1 : found ipsp=',ipsp
         write(string1,'(i1)')ipsp
         string(index_lower:index_lower+1)=blank//string1
         index_lower=index_lower+2
       end if
     end do ! ipsp
!    if not found ...
     if(found==0)then
       write(msg,'(6a)' )&
&       'Did not find matching pseudopotential for XYZ atomic symbol,',ch10,&
&       'with value ',string2,ch10,&
&       'Action: check that the atoms required by the XYZ file correspond to one psp file.'
       MSG_ERROR(msg)
     end if
   end do ! Loop on atoms
!  One should find blanks after the last significant type value
   string(index_lower:index_blank+2)=blank
 end do ! loop to identify _TYPAX

!---------------------------------------------------------------------------

! Here, set up quantities that are related to geometrical description of the system (acell,rprim,xred), as well as
! initial velocity(vel), charge and spin of atoms (chrgat,spinat), nuclear dipole moments of atoms (nucdipmom),
! the symmetries (symrel,symafm, and tnons) and the list of fixed atoms (iatfix,iatfixx,iatfixy,iatfixz).
! Arrays have already been dimensioned thanks to the knowledge of msym and mx%natom

!ji: We need to read the electric field before calling ingeo
!****** Temporary ******

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'berryopt',tread,'INT')
 if(tread==1) dtset%berryopt=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'berrysav',tread,'INT')
 if(tread==1) dtset%berrysav=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'bfield',tread,'DPR')
 if (tread==1) dtset%bfield(1:3) = dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'dfield',tread,'DPR')
 if (tread==1) dtset%dfield(1:3) = dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'efield',tread,'DPR')
 if (tread==1) dtset%efield(1:3) = dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'red_dfield',tread,'DPR')
 if (tread==1) dtset%red_dfield(1:3) = dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'red_efield',tread,'DPR')
 if (tread==1) dtset%red_efield(1:3) = dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'red_efieldbar',tread,'DPR')
 if (tread==1) dtset%red_efieldbar(1:3) = dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'jfielddir',tread,'INT')
 if(tread==1) dtset%jfielddir(1:3)=intarr(1:3)

 ! We need to know nsppol/nspinor/nspden before calling ingeo
 nsppol=dtset%nsppol
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nsppol',tread,'INT')
 if(tread==1) nsppol=intarr(1)

!Alternate SIESTA definition of nsppol
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'SpinPolarized',tread_alt,'LOG')
 if(tread_alt==1)then
   if(tread==1)then
     msg = 'nsppol and SpinPolarized cannot be specified simultaneously for the same dataset.'
     MSG_ERROR_NOSTOP(msg, leave)
   else
!    Note that SpinPolarized is a logical input variable
     nsppol=1
     if(intarr(1)==1)nsppol=2
     tread=1
   end if
 end if
 dtset%nsppol=nsppol

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nspinor',tread,'INT')
 if(tread==1) dtset%nspinor=intarr(1)

!Has to read pawspnorb now, in order to adjust nspinor
!Also, if nspinor=1, turn on spin-orbit coupling by default, here for the PAW case. NC case is treated elsewhere.
 if (dtset%usepaw>0)then
!  Change the default value
   if(dtset%nspinor==2)dtset%pawspnorb=1
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawspnorb',tread,'INT')
   if(tread==1)then
     dtset%pawspnorb=intarr(1)
     if(dtset%pawspnorb>0) dtset%nspinor=2
   else
     if(dtset%nspinor==2)then
       write(msg, '(4a)' ) ch10,&
&       ' invars1: COMMENT -',ch10,&
&       '  With nspinor=2 and usepaw=1, pawspnorb=1 has been switched on by default.'
       call wrtout(iout, msg,'COLL')
     end if
   end if
 end if
 nspinor=dtset%nspinor

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nspden',tread,'INT')
 if(tread==1) then
   dtset%nspden=intarr(1)
 else
   dtset%nspden=dtset%nsppol
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ntypalch',tread,'INT')
 if(tread==1) dtset%ntypalch=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nzchempot',tread,'INT')
 if(tread==1) dtset%nzchempot=intarr(1)

 ntypalch=dtset%ntypalch
 if(ntypalch>ntypat)then
   write(msg, '(3a,i0,a,i0,a,a)' )&
    'The input variable ntypalch must be smaller than ntypat, while it is',ch10,&
    'ntypalch=',dtset%ntypalch,', and ntypat=',ntypat,ch10,&
    'Action: check ntypalch vs ntypat in your input file.'
   MSG_ERROR(msg)
 end if

 ntyppure=ntypat-ntypalch
 dtset%ntyppure=ntyppure
 npspalch=npsp-ntyppure
 dtset%npspalch=npspalch
 if(npspalch<0)then
   write(msg, '(a,i0,2a,i0,a,a)' )&
    'The number of available pseudopotentials, npsp=',npsp,ch10,&
    'is smaller than the requested number of types of pure atoms, ntyppure=',ntyppure,ch10,&
    'Action: check ntypalch versus ntypat and npsp in your input file.'
   MSG_ERROR(msg)
 end if

 if(ntypalch>0)then
   call intagm(dprarr,intarr,jdtset,marr,ntypalch,string(1:lenstr),'algalch',tread,'INT')
   if(tread==1) dtset%algalch(1:ntypalch)=intarr(1:ntypalch)
   if (tread_geo /= 0) then
     MSG_ERROR("Alchemical mixing cannot be used with geo variable, use typat, znucl etc.")
   end if
 end if

!Read the Zeeman field
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'zeemanfield',tread,'BFI')
 if(tread==1) then
   if(dtset%nspden == 2)then
     write(msg,'(7a)')&
      'A Zeeman field has been specified without noncollinear spins.',ch10,&
      'Only the z-component of the magnetic field will be used.'
     MSG_WARNING(msg)
   else if (dtset%nspden == 1)then
     write(msg, '(a,a,a)' )&
      'A Zeeman field has been specified for a non-spin-polarized calculation.',ch10,&
      'Action: check the input file.'
     MSG_ERROR(msg)
   end if

   dtset%zeemanfield(1:3) = dprarr(1:3)
 end if

 ABI_ALLOCATE(amu,(ntypat))
 ABI_ALLOCATE(mixalch,(npspalch,ntypalch))
 ABI_ALLOCATE(vel,(3,natom))
 ABI_ALLOCATE(vel_cell,(3,3))
 ABI_ALLOCATE(xred,(3,natom))
 intimage=2 ; if(dtset%nimage==1)intimage=1
 do ii=1,dtset%nimage+1
   iimage=ii
   if(dtset%nimage==1 .and. ii==2)exit
   if(dtset%nimage==2 .and. ii==3)exit
   if(dtset%nimage> 2 .and. ii==intimage)cycle ! Will do the intermediate reference image at the last reading
   if(dtset%nimage>=2 .and. ii==dtset%nimage+1)iimage=intimage

   if (dtset%nimage /= 1) call wrtout(std_out, sjoin(' invars1: treat image number: ',itoa(iimage)))

!  Need to reset nsym to default value for each image
   dtset%nsym=0

!  Call ingeo for each image in turn, with the possible default values
   acell=dtset%acell_orig(1:3,iimage)
   amu=dtset%amu_orig(1:ntypat,iimage)
   mixalch=dtset%mixalch_orig(1:npspalch,1:ntypalch,iimage)
   rprim=dtset%rprim_orig(1:3,1:3,iimage)
   vel=dtset%vel_orig(1:3,1:natom,iimage)
   vel_cell=dtset%vel_cell_orig(1:3,1:3,iimage)
   xred=dtset%xred_orig(1:3,1:natom,iimage)
   ABI_ALLOCATE(chrgat,(natom))
   ABI_ALLOCATE(iatfix,(3,natom))
   ABI_ALLOCATE(nucdipmom,(3,natom))
   ABI_ALLOCATE(spinat,(3,natom))
   ABI_ALLOCATE(typat,(natom))
   ABI_ALLOCATE(znucl,(dtset%npsp))
   chrgat(1:natom)=dtset%chrgat(1:natom)
   nucdipmom(1:3,1:natom)=dtset%nucdipmom(1:3,1:natom)
   spinat(1:3,1:natom)=dtset%spinat(1:3,1:natom)
   znucl(1:dtset%npsp)=dtset%znucl(1:dtset%npsp)

   call ingeo(acell,amu,bravais,chrgat,dtset,dtset%genafm(1:3),iatfix,&
    dtset%icoulomb,iimage,iout,jdtset,dtset%jellslab,lenstr,mixalch,&
    msym,natom,dtset%nimage,dtset%npsp,npspalch,dtset%nspden,dtset%nsppol,&
    dtset%nsym,ntypalch,dtset%ntypat,nucdipmom,dtset%nzchempot,&
    dtset%pawspnorb,dtset%ptgroupma,ratsph,&
    rprim,dtset%slabzbeg,dtset%slabzend,dtset%spgroup,spinat,&
    string,dtset%supercell_latt,symafm,dtset%symmorphi,symrel,tnons,dtset%tolsym,&
    typat,vel,vel_cell,xred,znucl, comm)

   dtset%chrgat(1:natom)=chrgat(1:natom)
   dtset%iatfix(1:3,1:natom)=iatfix(1:3,1:natom)
   dtset%nucdipmom(1:3,1:natom)=nucdipmom(1:3,1:natom)
   dtset%spinat(1:3,1:natom)=spinat(1:3,1:natom)
   dtset%typat(1:natom)=typat(1:natom)
   ABI_DEALLOCATE(chrgat)
   ABI_DEALLOCATE(iatfix)
   ABI_DEALLOCATE(nucdipmom)
   ABI_DEALLOCATE(spinat)
   ABI_DEALLOCATE(typat)
   ABI_DEALLOCATE(znucl)
   dtset%acell_orig(1:3,iimage)=acell
   dtset%amu_orig(1:ntypat,iimage)=amu
   dtset%mixalch_orig(1:npspalch,1:ntypalch,iimage)=mixalch
   dtset%rprim_orig(1:3,1:3,iimage)=rprim
   dtset%vel_orig(1:3,1:natom,iimage)=vel
   dtset%vel_cell_orig(1:3,1:3,iimage)=vel_cell
   dtset%xred_orig(1:3,1:natom,iimage)=xred
   call mkrdim(dtset%acell_orig(1:3,iimage),dtset%rprim_orig(1:3,1:3,iimage),dtset%rprimd_orig(1:3,1:3,iimage))
 end do

 ABI_DEALLOCATE(amu)
 ABI_DEALLOCATE(mixalch)
 ABI_DEALLOCATE(vel)
 ABI_DEALLOCATE(vel_cell)
 ABI_DEALLOCATE(xred)

 ! Examine whether there is some vacuum space in the unit cell
 call invacuum(jdtset,lenstr,natom,dtset%rprimd_orig(1:3,1:3,intimage),string,vacuum,&
& dtset%xred_orig(1:3,1:natom,intimage))

!write(std_out,*)' invars1: before inkpts, dtset%mixalch_orig(1:npspalch,1:ntypalch,:)=',&
!dtset%mixalch_orig(1:npspalch,1:ntypalch,1:dtset%nimage)

!---------------------------------------------------------------------------

!Set up k point grid number
!First, get additional information
 dtset%kptopt=1
 if(dtset%nspden==4)dtset%kptopt=4
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'kptopt',tread,'INT')
 if(tread==1) dtset%kptopt=intarr(1)

 dtset%qptopt=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'qptopt',tread,'INT')
 if(tread==1) dtset%qptopt=intarr(1)

 iscf=5
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'iscf',tread,'INT')
 if(tread==1) iscf=intarr(1)

 dtset%natsph=dtset%natom
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natsph',tread,'INT')
 if(tread==1) dtset%natsph=intarr(1)

 dtset%natsph_extra=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natsph_extra',tread,'INT')
 if(tread==1) dtset%natsph_extra=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natvshift',tread,'INT')
 if(tread==1) dtset%natvshift=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nconeq',tread,'INT')
 if(tread==1) dtset%nconeq=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nkptgw',tread,'INT')
 if(tread==1) dtset%nkptgw=intarr(1)
 if (dtset%nkptgw<0) then
   write(msg, '(a,i0,4a)' )&
   'Input nkptgw must be >= 0, but was ',dtset%nkptgw,ch10,&
   'This is not allowed.',ch10,'Action: check the input file.'
   MSG_ERROR(msg)
 end if

 ! Number of points for long wavelength limit. Default is dtset%gw_nqlwl=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw_nqlwl',tread,'INT')
 if(tread==1) dtset%gw_nqlwl=intarr(1)
 if (dtset%gw_nqlwl<0) then
   write(msg, '(a,i0,4a)' )&
   'Input gw_nqlwl must be > 0, but was ',dtset%gw_nqlwl,ch10,&
   'This is not allowed.',ch10,'Action: check the input file.'
   MSG_ERROR(msg)
 end if

 ! Read number of k-points from input file (if specified)
 nkpt=0
 if(dtset%kptopt==0)nkpt=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nkpt',tread,'INT')
 if(tread==1) nkpt=intarr(1)

 ! or from KERANGE file.
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr), "getkerange_filepath", tread, 'KEY', key_value=key_value)
 if (tread==1) dtset%getkerange_filepath = key_value

#ifdef HAVE_NETCDF
 if (dtset%getkerange_filepath /= ABI_NOFILE) then
   ! Get number of k-points in sigma_erange energy windows.
   !dtset%kptopt = 0
   if (my_rank == master) then
     NCF_CHECK(nctk_open_read(ncid, dtset%getkerange_filepath, xmpi_comm_self))
     NCF_CHECK(nctk_get_dim(ncid, "nkpt_inerange", nkpt, datamode=.True.))
     NCF_CHECK(nf90_close(ncid))
   end if
   call xmpi_bcast(nkpt, master, comm, ierr)
 end if
#endif

 dtset%nkpt = nkpt

 call chkint_ge(0,0,cond_string,cond_values,ierr,'nkpt',nkpt,0,iout)
 if (dtset%kptopt==0) then
   cond_string(1)='kptopt'; cond_values(1)=0
   call chkint_ge(1,1,cond_string,cond_values,ierr,'nkpt',nkpt,1,iout)
 end if

 nkpthf=nkpt
 dtset%nkpthf=nkpt

 ! Will compute the actual value of nkpt, if needed. Otherwise,
 ! test that the value of nkpt is OK, if kptopt/=0
 ! Set up dummy arrays istwfk, kpt, wtk

 if(nkpt/=0 .or. dtset%kptopt/=0)then
   ABI_ALLOCATE(istwfk,(nkpt))
   ABI_ALLOCATE(kpt,(3,nkpt))
   ABI_ALLOCATE(kpthf,(3,nkpthf))
   ABI_ALLOCATE(wtk,(nkpt))
   ! Here, occopt is also a dummy argument
   occopt=1; dtset%nshiftk=1; dtset%kptrlatt(:,:)=0

   kptrlen=20.0_dp ; wtk(:)=1.0_dp
   dtset%shiftk(:,:)=half

   nqpt=0
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nqpt',tread,'INT')
   if(tread==1) nqpt=intarr(1)

   chksymbreak=1
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chksymbreak',tread,'INT')
   if(tread==1) chksymbreak=intarr(1)

   ! Use the first image to predict k and/or q points, except if an intermediate image is available
   intimage=1; if(dtset%nimage>2)intimage=(1+dtset%nimage)/2

   ! Find the q-point, if any.
   if(nqpt/=0)then
     call inqpt(chksymbreak,std_out,jdtset,lenstr,msym,natom,dtset%qptn,dtset%wtq,&
       dtset%rprimd_orig(1:3,1:3,intimage),dtset%spinat,string,dtset%typat,&
       vacuum,dtset%xred_orig(1:3,1:natom,intimage),dtset%qptrlatt)
   endif

   ! Find the k point grid
   call inkpts(bravais,chksymbreak,dtset%fockdownsampling,iout,iscf,istwfk,jdtset,&
     kpt,kpthf,dtset%kptopt,kptnrm,dtset%kptrlatt_orig,dtset%kptrlatt,kptrlen,lenstr,msym, dtset%getkerange_filepath, &
     nkpt,nkpthf,nqpt,dtset%ngkpt,dtset%nshiftk,dtset%nshiftk_orig,dtset%shiftk_orig,dtset%nsym,&
     occopt,dtset%qptn,response,dtset%rprimd_orig(1:3,1:3,intimage),dtset%shiftk,&
     string,symafm,symrel,vacuum,wtk,comm)

   ABI_DEALLOCATE(istwfk)
   ABI_DEALLOCATE(kpt)
   ABI_DEALLOCATE(kpthf)
   ABI_DEALLOCATE(wtk)

   ! nkpt and nkpthf have been computed, as well as the k point grid, if needed
   dtset%nkpt=nkpt
   dtset%nkpthf=nkpthf
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nqptdm',tread,'INT')
 if(tread==1) dtset%nqptdm=intarr(1)

 if (dtset%nqptdm<-1) then
   write(msg, '(a,i0,4a)' )&
    'Input nqptdm must be >= 0, but was ',dtset%nqptdm,ch10,&
    'This is not allowed.',ch10,'Action: check the input file.'
   MSG_ERROR(msg)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nzchempot',tread,'INT')
 if(tread==1) dtset%nzchempot=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'cd_customnimfrqs',tread,'INT')
 if(tread==1) dtset%cd_customnimfrqs=intarr(1)

 if (dtset%cd_customnimfrqs<0) then
   write(msg, '(a,i0,4a)' )&
    'Input cd_customnimfrqs must be >= 0, but was ',dtset%cd_customnimfrqs,ch10,&
    'This is not allowed.',ch10,'Action: check the input file.'
   MSG_ERROR(msg)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw_customnfreqsp',tread,'INT')
 if(tread==1) dtset%gw_customnfreqsp=intarr(1)

 if (dtset%gw_customnfreqsp<0) then
   write(msg, '(a,i0,4a)' )&
    'Input gw_customnfreqsp must be >= 0, but was ',dtset%gw_customnfreqsp,ch10,&
    'This is not allowed.',ch10,'Action: check the input file.'
   MSG_ERROR(msg)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_n_proj_freq',tread,'INT')
 if(tread==1) dtset%gwls_n_proj_freq=intarr(1)

 if (dtset%gwls_n_proj_freq<0) then
   write(msg, '(a,i0,4a)' )&
   'Input gwls_n_proj_freq must be >= 0, but was ',dtset%gwls_n_proj_freq,ch10,&
   'This is not allowed.',ch10,'Action: check the input file.'
   MSG_ERROR(msg)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'efmas_calc_dirs',tread,'INT')
 if(tread==1) dtset%efmas_calc_dirs=intarr(1)

 if (ABS(dtset%efmas_calc_dirs)>3) then
   write(msg, '(a,i0,4a)' )&
   'Input efmas_calc_dirs must be between -3 and 3, but was ',dtset%efmas_calc_dirs,ch10,&
   'This is not allowed.',ch10,'Action: check the input file.'
   MSG_ERROR(msg)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'efmas_n_dirs',tread,'INT')
 if(tread==1) dtset%efmas_n_dirs=intarr(1)

 if (dtset%efmas_n_dirs<0) then
   write(msg, '(a,i0,4a)' )&
   'Input efmas_n_dirs must be >= 0, but was ',dtset%efmas_n_dirs,ch10,&
   'This is not allowed.',ch10,'Action: check the input file.'
   MSG_ERROR(msg)
 end if

!---------------------------------------------------------------------------

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nnos',tread,'INT')
 if(tread==1) dtset%nnos=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ga_n_rules',tread,'INT')
 if(tread==1) dtset%ga_n_rules=intarr(1)

 ! Perform the first checks
 ! Check that nkpt is greater than 0
 if (nkpt<=0) then
   write(msg, '(a,i0)' )'After inkpts, nkpt must be > 0, but was ',nkpt
   MSG_ERROR_NOSTOP(msg, leave)
 end if

 ! Check that nsppol is 1 or 2
 if (nsppol/=1 .and. nsppol/=2) then
   write(msg, '(a,i0)' )'Input nsppol must be 1 or 2, but was ',nsppol
   MSG_ERROR_NOSTOP(msg, leave)
 end if

 ! Check that nspinor is 1 or 2
 if (nspinor/=1 .and. nspinor/=2) then
   write(msg, '(a,i0)' )'Input nspinor must be 1 or 2, but was ',nspinor
   MSG_ERROR_NOSTOP(msg, leave)
 end if

 ! Check that nspinor and nsppol are not 2 together
 if (nsppol==2 .and. nspinor==2) then
   MSG_ERROR_NOSTOP('nspinor and nsppol cannot be 2 together!', leave)
 end if

 ! Here, leave if an error has been detected earlier
 if (leave /= 0) then
   MSG_ERROR('Errors are present in the input file. See ABOVE messages')
 end if

 ! Now, take care of mband_upper
 mband_upper=1
 occopt=1
 fband=0.5_dp

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'occopt',tread,'INT')
 if(tread==1) occopt=intarr(1)

 ! Also read fband, that is an alternative to nband. The default
 ! is different for occopt==1 and for metallic occupations.
 if(occopt==1)fband=0.125_dp
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fband',tfband,'DPR')
 if(tfband==1)fband=dprarr(1)

 ! fband cannot be used when occopt==0 or occopt==2
 if(tfband==1 .and. (occopt==0 .or. occopt==2) )then
   write(msg, '(3a)' )&
   'fband cannot be used if occopt==0 or occopt==2 ',ch10,&
   'Action: correct your input file, suppress fband, or change occopt.'
   MSG_ERROR(msg)
 end if

 ABI_ALLOCATE(nband,(nkpt*nsppol))
 tnband=0

 ! Compute ziontypat
 ! When the pseudo-atom is pure, simple copy
 if(ntyppure>0)then
   do itypat=1,ntyppure
     dtset%ziontypat(itypat)=zionpsp(itypat)
   end do
 end if

 ! When the pseudo-atom is alchemical, must make mixing
 if(ntypalch>0)then
   do itypat=ntyppure+1,ntypat
     dtset%ziontypat(itypat)=zero
     do ipsp=ntyppure+1,npsp
       dtset%ziontypat(itypat)=dtset%ziontypat(itypat) &
&       +dtset%mixalch_orig(ipsp-ntyppure,itypat-ntyppure,1)*zionpsp(ipsp)
     end do
   end do
 end if

 if (occopt==0 .or. occopt==1 .or. (occopt>=3 .and. occopt<=8) ) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nband',tnband,'INT')
   ! Note: mband_upper is initialized, not nband
   if(tnband==1) mband_upper=intarr(1)

   if(tfband==1 .and. tnband==1)then
     write(msg, '(3a)' )&
     'fband and nband cannot be used together. ',ch10,&
     'Action: correct your input file, suppress either fband or nband.'
     MSG_ERROR(msg)
   end if

   ! In case nband was not read, use fband, either read, or the default,
   ! to provide an upper limit for mband_upper
   if(tnband==0)then

     charge=0.0_dp
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'charge',tread,'DPR')
     if(tread==1) charge=dprarr(1)

     ! Only take into account negative charge, to compute maximum number of bands
     if(charge > 0.0_dp)charge=0.0_dp

!     mband_upper=nspinor*((nint(zion_max)*natom+1)/2 - floor(charge/2.0_dp)&
!&     + ceiling(fband*natom-1.0d-10))
     zval=0.0_dp
     do iatom=1,natom
       zval=zval+dtset%ziontypat(dtset%typat(iatom))
     end do
     zelect=zval-charge
     mband_upper=nspinor * ((ceiling(zelect-1.0d-10)+1)/2 + ceiling( fband*natom - 1.0d-10 ))

     nband(:)=mband_upper

!    write(std_out,*)' invars1 : zion_max,natom,fband,mband_upper '
!    write(std_out,*)zion_max,natom,fband,mband_upper
   end if

   nband(:)=mband_upper

 else if (occopt==2) then
   ABI_ALLOCATE(reaalloc,(nkpt*nsppol))
   call intagm(reaalloc,nband,jdtset,nkpt*nsppol,nkpt*nsppol,string(1:lenstr),'nband',tnband,'INT')
   if(tnband==1)then
     do ikpt=1,nkpt*nsppol
       if (nband(ikpt)>mband_upper) mband_upper=nband(ikpt)
     end do
   end if
   ABI_DEALLOCATE(reaalloc)
 else
   write(msg, '(a,i0,3a)' )'occopt=',occopt,' is not an allowed value.',ch10,'Action: correct your input file.'
   MSG_ERROR(msg)
 end if

 ! Check that mband_upper is greater than 0
 if (mband_upper<=0) then
   write(msg, '(a,i0,4a)' )&
   'Maximal nband must be > 0, but was ',mband_upper,ch10,&
   'This is not allowed.',ch10,'Action: check the input file.'
   MSG_ERROR(msg)
 end if

 ! The following 3 values are needed to dimension the parallelism over images
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'imgmov',tread,'INT')
 if(tread==1) dtset%imgmov=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ntimimage',tread,'INT')
 if(tread==1) dtset%ntimimage=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,dtset%nimage,string(1:lenstr),'dynimage',tread,'INT')
 if(tread==1)then
   dtset%dynimage(1:dtset%nimage)=intarr(1:dtset%nimage)
 else if (dtset%imgmov==2.or.dtset%imgmov==5) then
   dtset%dynimage(1)=0;dtset%dynimage(dtset%nimage)=0
 end if
 dtset%ndynimage=count(dtset%dynimage(1:dtset%nimage)/=0)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'wfoptalg',tread,'INT')
 if(tread==1) then
   dtset%wfoptalg=intarr(1)
 else
   if (dtset%usepaw==0)    dtset%wfoptalg=0
   if (dtset%usepaw/=0)    dtset%wfoptalg=10
   if (dtset%optdriver==RUNL_GSTATE) then
     if (dtset%paral_kgb/=0) dtset%wfoptalg=14
   end if
 end if

!---------------------------------------------------------------------------
!Some PAW+DMFT keywords
 dtset%usedmft=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'usedmft',tread,'INT')
 if(tread==1) dtset%usedmft=intarr(1)

!Some ucrpa keywords
 dtset%ucrpa=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ucrpa',tread,'INT')
 if(tread==1) dtset%ucrpa=intarr(1)

 if (dtset%ucrpa > 0 .and. dtset%usedmft > 0) then
   write(msg, '(9a)' )&
   'usedmft and ucrpa are both activated in the input file ',ch10,&
   'In the following, abinit assume you are doing a ucrpa calculation and ',ch10,&
   'you define Wannier functions as in DFT+DMFT calculation',ch10,&
   'If instead, you want to do a full dft+dmft calculation and not only the Wannier construction, use ucrpa=0',ch10,&
   'This keywords are depreciated, please use the new keywords to perform cRPA calculation'
   MSG_WARNING(msg)
 end if

!Some PAW+U keywords
 dtset%usepawu=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'usepawu',tread,'INT')
 if(tread==1) dtset%usepawu=intarr(1)
 if ( dtset%usedmft > 0 .and. dtset%usepawu >= 0 ) dtset%usepawu = 1

 dtset%usedmatpu=0
 dtset%lpawu(1:dtset%ntypat)=-1
 if (dtset%usepawu/=0.or.dtset%usedmft>0) then
   call intagm(dprarr,intarr,jdtset,marr,dtset%ntypat,string(1:lenstr),'lpawu',tread,'INT')
   if(tread==1) dtset%lpawu(1:dtset%ntypat)=intarr(1:dtset%ntypat)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'usedmatpu',tread,'INT')
   if(tread==1) dtset%usedmatpu=intarr(1)
 end if

!Some PAW+Exact exchange keywords
 dtset%useexexch=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'useexexch',tread,'INT')
 if(tread==1) dtset%useexexch=intarr(1)

 dtset%lexexch(1:dtset%ntypat)=-1

 if (dtset%useexexch/=0) then
   call intagm(dprarr,intarr,jdtset,marr,dtset%ntypat,string(1:lenstr),'lexexch',tread,'INT')
   if(tread==1) dtset%lexexch(1:dtset%ntypat)=intarr(1:dtset%ntypat)
 end if

!LDA minus half keyword
 call intagm(dprarr,intarr,jdtset,marr,dtset%ntypat,string(1:lenstr),'ldaminushalf',tread,'INT')
 if(tread==1) dtset%ldaminushalf(1:dtset%ntypat)=intarr(1:dtset%ntypat)

!Some plowan data
 dtset%plowan_natom=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'plowan_natom',tread,'INT')
 if(tread==1) dtset%plowan_natom=intarr(1)

 dtset%plowan_nt=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'plowan_nt',tread,'INT')
 if(tread==1) dtset%plowan_natom=intarr(1)

 !if (dtset%ucrpa > 0 .and. dtset%plowan_compute==0) then
   !dtset%plowan_natom=1
   !dtset%plowan_nt=1
 !endif

!PAW potential zero keyword
 dtset%usepotzero=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'usepotzero',tread,'INT')
 if(tread==1) dtset%usepotzero=intarr(1)

!Macro_uj (determination of U in PAW+U), governs also allocation of atvshift
 dtset%macro_uj = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'macro_uj',tread,'INT')
 if(tread==1) dtset%macro_uj=intarr(1)

!Constraint DFT keyword
 call intagm(dprarr,intarr,jdtset,marr,dtset%ntypat,string(1:lenstr),'constraint_kind',tread,'INT')
 if(tread==1) dtset%constraint_kind(1:dtset%ntypat)=intarr(1:dtset%ntypat)

 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(ratsph)
 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

end subroutine invars1
!!***

!!****f* ABINIT/indefo
!! NAME
!! indefo
!!
!! FUNCTION
!! Initialisation phase: default values for most input variables
!! (some are initialized earlier, see indefo1 routine)
!!
!! INPUTS
!!  ndtset_alloc=number of datasets, corrected for allocation of at least one data set.
!!  nprocs=Number of MPI processors available.
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are given a default value here.
!!   The dataset with number 0 should be the reference default value in the remaining of the code.
!!
!! NOTES
!! The outputs of this routine are the defaults values of input
!! variables, stored at the index 0 of the last dimension of their multi-dataset representation.
!!
!!  Scalars and static arrays can be initialized directly at the level of the datatype declaration
!!  provided the value does not depend on runtime conditions.
!!
!! PARENTS
!!      m_ab7_invars_f90
!!
!! CHILDREN
!!
!! SOURCE

subroutine indefo(dtsets,ndtset_alloc,nprocs)

 use m_gwdefs
#if defined DEV_YP_VDWXC
 use m_xc_vdw
#endif

 use m_fftcore,      only : get_cache_kb, fftalg_for_npfft

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndtset_alloc,nprocs
!arrays
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc) !vz_i

!Local variables -------------------------------
!scalars
 integer :: idtset,ii,jdtset,paral_atom_default
 logical :: wvl_bigdft
#if defined DEV_YP_VDWXC
 type(xc_vdw_type) :: vdw_defaults
#endif

!******************************************************************

 DBG_ENTER("COLL")

!Set up default values. All variables to be output in outvars.f
!should have a default, even if a nonsensible one can be chosen to garantee print in that routine.

!These variables have already been initialized, for idtset/=0
 dtsets(0)%istatr=0
 dtsets(0)%istatshft=1
 dtsets(0)%kptrlatt(1:3,1:3)=0
 !dtsets(0)%kptrlatt_orig=0
 dtsets(0)%qptrlatt(1:3,1:3)=0
 dtsets(0)%ptgroupma=0
 dtsets(0)%spgroup=0
 dtsets(0)%shiftk(:,:)=half
 dtsets(0)%tolsym=tol8
 dtsets(0)%znucl(:)=zero
 dtsets(0)%ucrpa=0
 dtsets(0)%usedmft=0

 paral_atom_default=0
 if (nprocs>1.and.maxval(dtsets(:)%usepaw)>0) paral_atom_default=1

!WARNING: set default in all datasets, including idtset=0 !!!
!Use alphabetic order

 do idtset=0,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset

   wvl_bigdft=.false.
   if(dtsets(idtset)%usewvl==1 .and. dtsets(idtset)%wvl_bigdft_comp==1) wvl_bigdft=.true.
!  Special case of use_gpu_cuda (can be undertermined at this point)
!  use_gpu_cuda=-1 means undetermined ; here impose its value due to some restrictions
   if (dtsets(idtset)%use_gpu_cuda==-1) then
     if (dtsets(idtset)%optdriver/=0.or. dtsets(idtset)%tfkinfunc/=0.or. dtsets(idtset)%nspinor/=1) then
       dtsets(idtset)%use_gpu_cuda=0
     else
       dtsets(idtset)%use_gpu_cuda=1
     end if
   end if

!  A
!  Here we change the default value of iomode according to the configuration options.
!  Ideally, all the sequential tests should pass independently of the default value.
!  The parallel tests may require IO_MODE_MPI or, alternatively, IO_MODE_ETSF with HDF5 support.
!  MG FIXME Sun Sep 6 2015: Many tests fail if IO_MODE_MPI is used as default. IO errors in v1, v2 ...
!  with np=1 and wonderful deadlocks if np>1.

!  Note that this default value might be overriden for specific datasets later, in case of parallelism
   dtsets(idtset)%iomode=IO_MODE_FORTRAN
#ifdef HAVE_NETCDF_DEFAULT
   dtsets(idtset)%iomode=IO_MODE_ETSF
#endif
#ifdef HAVE_MPI_IO_DEFAULT
   dtsets(idtset)%iomode=IO_MODE_MPI
#endif

   dtsets(idtset)%adpimd=0
   dtsets(idtset)%adpimd_gamma=one
   dtsets(idtset)%accuracy=0
   dtsets(idtset)%atvshift(:,:,:)=zero
   dtsets(idtset)%auxc_ixc=11
   dtsets(idtset)%auxc_scal=one
   dtsets(idtset)%awtr=1
!  B
   dtsets(idtset)%bdberry(1:4)=0
   dtsets(idtset)%bdeigrf=-1
   dtsets(idtset)%bdgw=0
   dtsets(idtset)%berrystep=1
   dtsets(idtset)%bmass=ten
   dtsets(idtset)%boxcenter(1:3)=half
   dtsets(idtset)%boxcutmin=two
   dtsets(idtset)%brvltt=0
   dtsets(idtset)%bs_nstates=0
!  dtsets(idtset)%bs_hayd_term=0
   dtsets(idtset)%bs_hayd_term=1
   dtsets(idtset)%builtintest=0
   dtsets(idtset)%bxctmindg=two
!  C
   dtsets(idtset)%cd_halfway_freq=3.674930883_dp !(100 eV)
   dtsets(idtset)%cd_imfrqs(:) = zero
   dtsets(idtset)%cd_max_freq=36.74930883_dp     !(1000 eV)
   dtsets(idtset)%cd_subset_freq(1:2)=0
   dtsets(idtset)%cd_frqim_method=1
   dtsets(idtset)%cd_full_grid=0
   dtsets(idtset)%charge=zero
   dtsets(idtset)%chempot(:,:,:)=zero
   dtsets(idtset)%chkdilatmx=1
   dtsets(idtset)%chkexit=0
   dtsets(idtset)%chksymbreak=1
   dtsets(idtset)%cineb_start=7
   dtsets(idtset)%corecs(:) = zero
!  D
   dtsets(idtset)%ddamp=0.1_dp
   dtsets(idtset)%delayperm=0
   dtsets(idtset)%densfor_pred=2
   if (dtsets(idtset)%paral_kgb>0.and.idtset>0) dtsets(idtset)%densfor_pred=6 ! Recommended for band-FFT parallelism
!XG170502 : This section is completely useless, as ionmov is NOT know at present !
!#ifdef HAVE_LOTF
!   if (dtsets(idtset)%ionmov==23) dtsets(idtset)%densfor_pred=2 ! Recommended for LOTF
!#endif
   dtsets(idtset)%dfpt_sciss=zero
   dtsets(idtset)%diecut=2.2_dp
   dtsets(idtset)%dielng=1.0774841_dp
   dtsets(idtset)%diemac=1.0d6
   if (dtsets(idtset)%usepaw==0) then
     dtsets(idtset)%diemix=one
   else
     dtsets(idtset)%diemix=0.7_dp
   end if
   dtsets(idtset)%diemixmag=dtsets(idtset)%diemix
   dtsets(idtset)%diegap=0.1_dp
   dtsets(idtset)%dielam=half
   dtsets(idtset)%diismemory=8
   dtsets(idtset)%dilatmx=one
   dtsets(idtset)%dmatpuopt=2
   if (size(dtsets(idtset)%dmatpawu,4)>0) dtsets(idtset)%dmatpawu=-10._dp
   dtsets(idtset)%dmatudiag=0
   dtsets(idtset)%dmft_entropy=0
   dtsets(idtset)%dmft_dc  =1
   dtsets(idtset)%dmft_iter=0
   dtsets(idtset)%dmft_kspectralfunc=0
   dtsets(idtset)%dmft_nlambda=6
   dtsets(idtset)%dmft_nwli=0
   dtsets(idtset)%dmft_nwlo=0
   dtsets(idtset)%dmft_mxsf=0.3_dp
   dtsets(idtset)%dmft_occnd_imag=1
   dtsets(idtset)%dmft_read_occnd=0
   dtsets(idtset)%dmft_rslf=0
   dtsets(idtset)%dmft_solv=5
   if(dtsets(idtset)%ucrpa>0.and.dtsets(idtset)%usedmft==1) dtsets(idtset)%dmft_solv=0
   dtsets(idtset)%dmft_t2g=0
!  dtsets(idtset)%dmft_x2my2d=0
   dtsets(idtset)%dmft_tolfreq=tol4
   dtsets(idtset)%dmft_tollc=tol5
   dtsets(idtset)%dmft_charge_prec=tol6
   dtsets(idtset)%dmftbandi=0
   dtsets(idtset)%dmftbandf=0
   dtsets(idtset)%dmftcheck=0
   dtsets(idtset)%dmftctqmc_basis =1
   dtsets(idtset)%dmftctqmc_check =0
   dtsets(idtset)%dmftctqmc_correl=0
   dtsets(idtset)%dmftctqmc_grnns =0
   dtsets(idtset)%dmftctqmc_meas  =1
   dtsets(idtset)%dmftctqmc_mrka  =0
   dtsets(idtset)%dmftctqmc_mov   =0
   dtsets(idtset)%dmftctqmc_order =0
   dtsets(idtset)%dmftctqmc_triqs_nleg=30
   dtsets(idtset)%dmftqmc_l=0
   dtsets(idtset)%dmftqmc_n=0.0_dp
   dtsets(idtset)%dmftqmc_seed=jdtset
   dtsets(idtset)%dmftqmc_therm=1000
   dtsets(idtset)%dmftctqmc_gmove = dtsets(idtset)%dmftqmc_therm / 10
   dtsets(idtset)%dosdeltae=zero
   dtsets(idtset)%dtion=100.0_dp
   dtsets(idtset)%d3e_pert1_atpol(1:2)=1
   dtsets(idtset)%d3e_pert1_dir(1:3)=0
   dtsets(idtset)%d3e_pert1_elfd=0
   dtsets(idtset)%d3e_pert1_phon=0
   dtsets(idtset)%d3e_pert2_atpol(1:2)=1
   dtsets(idtset)%d3e_pert2_dir(1:3)=0
   dtsets(idtset)%d3e_pert2_elfd=0
   dtsets(idtset)%d3e_pert2_phon=0
   dtsets(idtset)%d3e_pert2_strs=0
   dtsets(idtset)%d3e_pert3_atpol(1:2)=1
   dtsets(idtset)%d3e_pert3_dir(1:3)=0
   dtsets(idtset)%d3e_pert3_elfd=0
   dtsets(idtset)%d3e_pert3_phon=0
!  E
   dtsets(idtset)%ecut=-one
   dtsets(idtset)%ecuteps=zero
   dtsets(idtset)%ecutsigx=zero ! The true default value is ecut . This is defined in invars2.F90
   dtsets(idtset)%ecutsm=zero
   dtsets(idtset)%ecutwfn=zero ! The true default value is ecut . This is defined in invars2.F90
   dtsets(idtset)%effmass_free=one
   dtsets(idtset)%efmas=0
   dtsets(idtset)%efmas_bands=0 ! The true default is nband. This is defined in invars2.F90
   dtsets(idtset)%efmas_deg=1
   dtsets(idtset)%efmas_deg_tol=tol5
   dtsets(idtset)%efmas_dim=3
   dtsets(idtset)%efmas_dirs=zero
   dtsets(idtset)%efmas_ntheta=1000
   dtsets(idtset)%elph2_imagden=zero
   dtsets(idtset)%enunit=0
   dtsets(idtset)%eshift=zero
   dtsets(idtset)%esmear=0.01_dp
   dtsets(idtset)%exchn2n3d=0
   dtsets(idtset)%extrapwf=0
   dtsets(idtset)%exchmix=quarter
!  F
   dtsets(idtset)%fermie_nest=zero
   dtsets(idtset)%fftgw=21
   dtsets(idtset)%focktoldfe=zero
   dtsets(idtset)%fockoptmix=0
   dtsets(idtset)%fockdownsampling(:)=1
   dtsets(idtset)%freqim_alpha=five
   dtsets(idtset)%freqremin=zero
   dtsets(idtset)%freqremax=zero
   dtsets(idtset)%freqspmin=zero
   dtsets(idtset)%freqspmax=zero
   dtsets(idtset)%friction=0.001_dp
   dtsets(idtset)%frzfermi=0
   dtsets(idtset)%fxcartfactor=one ! Should be adjusted to the H2 conversion factor
!  G
   dtsets(idtset)%ga_algor =1
   dtsets(idtset)%ga_fitness =1
   dtsets(idtset)%ga_opt_percent =0.2_dp
   dtsets(idtset)%ga_rules(:) =1
   dtsets(idtset)%goprecon =0
   dtsets(idtset)%goprecprm(:)=0
   dtsets(idtset)%gpu_devices=(/-1,-1,-1,-1,-1/)
   dtsets(idtset)%gpu_linalg_limit=2000000
   if (dtsets(idtset)%gw_customnfreqsp/=0) dtsets(idtset)%gw_freqsp(:) = zero
   dtsets(idtset)%gw_nstep =30
   dtsets(idtset)%gwgamma =0
   if ( dtsets(idtset)%gw_nqlwl > 0 ) then
     dtsets(idtset)%gw_qlwl(:,:)=zero
     dtsets(idtset)%gw_qlwl(1,1)=0.00001_dp
     dtsets(idtset)%gw_qlwl(2,1)=0.00002_dp
     dtsets(idtset)%gw_qlwl(3,1)=0.00003_dp
   end if
   dtsets(idtset)%gw_frqim_inzgrid=0
   dtsets(idtset)%gw_frqre_inzgrid=0
   dtsets(idtset)%gw_frqre_tangrid=0
   dtsets(idtset)%gw_invalid_freq=0
   dtsets(idtset)%gw_qprange=0
   dtsets(idtset)%gw_sigxcore=0
   dtsets(idtset)%gw_sctype = GWSC_one_shot
   dtsets(idtset)%gw_toldfeig=0.1/Ha_eV
   dtsets(idtset)%getbseig=0
   dtsets(idtset)%getbsreso=0
   dtsets(idtset)%getbscoup=0
   dtsets(idtset)%getcell =0
   dtsets(idtset)%getddb  =0
   dtsets(idtset)%getddk  =0
   dtsets(idtset)%getdelfd=0
   dtsets(idtset)%getdkdk =0
   dtsets(idtset)%getdkde =0
   dtsets(idtset)%getden  =0
   dtsets(idtset)%getefmas=0
   dtsets(idtset)%getgam_eig2nkq  =0
   dtsets(idtset)%gethaydock=0
   dtsets(idtset)%getocc  =0
   dtsets(idtset)%getpawden=0
   dtsets(idtset)%getqps  =0
   dtsets(idtset)%getscr  =0
   dtsets(idtset)%getsuscep=0
   dtsets(idtset)%getvel  =0
   dtsets(idtset)%getwfk  =0
   dtsets(idtset)%getwfkfine = 0
   dtsets(idtset)%getwfq  =0
   dtsets(idtset)%getxcart=0
   dtsets(idtset)%getxred =0
   dtsets(idtset)%get1den =0
   dtsets(idtset)%get1wf  =0
   dtsets(idtset)%gwcalctyp=0
   dtsets(idtset)%gwcomp=0
   dtsets(idtset)%gwencomp=2.0_dp
   dtsets(idtset)%gwmem=11
   dtsets(idtset)%gwpara=2
   dtsets(idtset)%gwrpacorr=0
   dtsets(idtset)%gwls_stern_kmax=1
   dtsets(idtset)%gwls_model_parameter=1.0_dp
   dtsets(idtset)%gwls_npt_gauss_quad=10
   dtsets(idtset)%gwls_diel_model=2
   dtsets(idtset)%gwls_print_debug=0
   if (dtsets(idtset)%gwls_n_proj_freq/=0) dtsets(idtset)%gwls_list_proj_freq(:) = zero
   dtsets(idtset)%gwls_nseeds=1
   dtsets(idtset)%gwls_recycle=2
   dtsets(idtset)%gwls_kmax_complement=1
   dtsets(idtset)%gwls_kmax_poles=4
   dtsets(idtset)%gwls_kmax_analytic=8
   dtsets(idtset)%gwls_kmax_numeric=16
   dtsets(idtset)%gwls_band_index=1
   dtsets(idtset)%gwls_exchange=1
   dtsets(idtset)%gwls_correlation=3
   dtsets(idtset)%gwls_first_seed=0
!  H
   dtsets(idtset)%hmcsst=3
   dtsets(idtset)%hmctt=4
   dtsets(idtset)%hyb_mixing=-999.0_dp
   dtsets(idtset)%hyb_mixing_sr=-999.0_dp
   dtsets(idtset)%hyb_range_dft=-999.0_dp
   dtsets(idtset)%hyb_range_fock=-999.0_dp
!  I
   if(dtsets(idtset)%natsph/=0) then
!    do not use iatsph(:) but explicit boundaries
!    to avoid to read to far away in the built array (/ ... /)
     dtsets(idtset)%iatsph(1:dtsets(idtset)%natsph)=(/ (ii,ii=1,dtsets(idtset)%natsph) /)
   else
     dtsets(idtset)%iatsph(:)=0
   end if
   dtsets(idtset)%iboxcut=0
   dtsets(idtset)%icsing=6
   dtsets(idtset)%icutcoul=6
   dtsets(idtset)%ieig2rf=0
   dtsets(idtset)%imgwfstor=0
   dtsets(idtset)%inclvkb=2
   dtsets(idtset)%intxc=0
!  if (dtsets(idtset)%paral_kgb>0.and.idtset>0) dtsets(idtset)%intxc=0
   dtsets(idtset)%ionmov=0
   dtsets(idtset)%densfor_pred=2
   if (dtsets(idtset)%paral_kgb>0.and.idtset>0) dtsets(idtset)%densfor_pred=6 ! Recommended for band-FFT parallelism
!This section is completely useless, as ionmov is NOT know at present !
!#ifdef HAVE_LOTF
!   if (dtsets(idtset)%ionmov==23) dtsets(idtset)%densfor_pred=2 ! Recommended for LOTF
!#endif
   dtsets(idtset)%iprcel=0
   dtsets(idtset)%iprcfc=0
   dtsets(idtset)%irandom=3
   dtsets(idtset)%irdbseig=0
   dtsets(idtset)%irdbsreso=0
   dtsets(idtset)%irdbscoup=0
   dtsets(idtset)%irdddb=0
   dtsets(idtset)%irdddk=0
   dtsets(idtset)%irdden=0
   dtsets(idtset)%irdefmas=0
   dtsets(idtset)%irdhaydock=0
   dtsets(idtset)%irdpawden=0
   dtsets(idtset)%irdqps=0
   dtsets(idtset)%irdscr=0
   dtsets(idtset)%irdsuscep=0
   dtsets(idtset)%irdvdw=0
   dtsets(idtset)%irdwfk=0
   dtsets(idtset)%irdwfkfine=0
   dtsets(idtset)%irdwfq=0
   dtsets(idtset)%ird1den=0
   dtsets(idtset)%ird1wf=0
!iscf
   if(wvl_bigdft) then
     dtsets(idtset)%iscf=0
   else
     if(dtsets(idtset)%usepaw==0) then
       dtsets(idtset)%iscf=7
     else
       dtsets(idtset)%iscf=17
     end if
   end if
   dtsets(idtset)%isecur=0
   dtsets(idtset)%istatimg = 1
   dtsets(idtset)%istwfk(:)=0
   dtsets(idtset)%ixc=1
   dtsets(idtset)%ixc_sigma=1
   dtsets(idtset)%ixcpositron=1
   dtsets(idtset)%ixcrot=1
!  J
   dtsets(idtset)%f4of2_sla(:)=-one
   dtsets(idtset)%f6of2_sla(:)=-one
   dtsets(idtset)%jpawu(:,:)=zero
!  K
   dtsets(idtset)%kberry(1:3,:)=0
   dtsets(idtset)%kpt(:,:)=zero
   dtsets(idtset)%kptgw(:,:)=zero
   dtsets(idtset)%kptnrm=one
   dtsets(idtset)%kptns_hf(:,:)=zero
   dtsets(idtset)%kptopt=1
   if(dtsets(idtset)%nspden==4)dtsets(idtset)%kptopt=4
   dtsets(idtset)%kptrlen=30.0_dp
   dtsets(idtset)%kssform=1
!  L
   dtsets(idtset)%localrdwf=1

#if defined HAVE_LOTF
   dtsets(idtset)%lotf_classic=5
   dtsets(idtset)%lotf_nitex=10
   dtsets(idtset)%lotf_nneigx=40
   dtsets(idtset)%lotf_version=2
#endif
   dtsets(idtset)%lw_qdrpl=0
   dtsets(idtset)%lw_flexo=0
!  M
   dtsets(idtset)%magconon = 0
   dtsets(idtset)%magcon_lambda = 0.01_dp
   dtsets(idtset)%max_ncpus = 0
   dtsets(idtset)%mbpt_sciss=zero
   dtsets(idtset)%mband = -1
   dtsets(idtset)%mdf_epsinf = zero
   dtsets(idtset)%mdtemp(:)=300.0_dp
   dtsets(idtset)%mdwall=10000_dp
   dtsets(idtset)%mep_mxstep=100._dp
   dtsets(idtset)%mep_solver=0
   dtsets(idtset)%mffmem=1
   dtsets(idtset)%mgfft = -1
   dtsets(idtset)%mgfftdg = -1
   dtsets(idtset)%mixesimgf(:)=zero
   dtsets(idtset)%mpw = -1
   dtsets(idtset)%mqgrid=0
   dtsets(idtset)%mqgriddg=0
!  N
   dtsets(idtset)%natrd = -1
   dtsets(idtset)%nband(:)=0
   dtsets(idtset)%nbandhf=0
   dtsets(idtset)%nbdblock=1
   dtsets(idtset)%nbdbuf=0
   dtsets(idtset)%nberry=1
   if (dtsets(idtset)%usepaw==0) then
     dtsets(idtset)%nc_xccc_gspace=0
   else
     dtsets(idtset)%nc_xccc_gspace=1
   end if
   dtsets(idtset)%nbandkss=0
   dtsets(idtset)%nctime=0
   dtsets(idtset)%ndtset = -1
   dtsets(idtset)%neb_algo=1
   dtsets(idtset)%neb_spring(1:2)=(/0.05_dp,0.05_dp/)
   dtsets(idtset)%npwkss=0
   dtsets(idtset)%nfft = -1
   dtsets(idtset)%nfftdg = -1

   dtsets(idtset)%nfreqim=-1
   dtsets(idtset)%nfreqre=-1
   dtsets(idtset)%nfreqsp=0

   dtsets(idtset)%npulayit=7

!  ngfft is a special case
   dtsets(idtset)%ngfft(1:8)=0
   dtsets(idtset)%ngfft(7) = fftalg_for_npfft(1)
!  fftcache=ngfft(8) is machine-dependent.
   dtsets(idtset)%ngfft(8) = get_cache_kb()

   dtsets(idtset)%ngfftdg(:)=dtsets(idtset)%ngfft(:)
!
   !nline
   dtsets(idtset)%nline=4
   if(dtsets(idtset)%usewvl==1 .and. .not. wvl_bigdft) then
     if(dtsets(idtset)%usepaw==1) then
       dtsets(idtset)%nline=4
     else
       dtsets(idtset)%nline=2
     end if
   end if

!  nloalg is also a special case
   dtsets(idtset)%nloalg(1)=4
   dtsets(idtset)%nloalg(2)=1
   dtsets(idtset)%nloalg(3)=dtsets(idtset)%usepaw
   !if (dtsets(idtset)%optdriver == RUNL_EPH) dtsets(idtset)%nloalg(3) = 1
   dtsets(idtset)%ngkpt=0
   dtsets(idtset)%nnsclo=0
   dtsets(idtset)%nnsclohf=0
   dtsets(idtset)%nomegasf=100
   dtsets(idtset)%nomegasrd=9
   dtsets(idtset)%nomegasi=12
   dtsets(idtset)%nonlinear_info=0
   dtsets(idtset)%noseinert=1.0d5
   dtsets(idtset)%npvel=0
   dtsets(idtset)%npweps=0
   dtsets(idtset)%npwsigx=0
   dtsets(idtset)%npwwfn=0
   dtsets(idtset)%nqpt=0
   dtsets(idtset)%nscforder=16
   dtsets(idtset)%nshiftk=1
   dtsets(idtset)%nshiftk_orig=1
   dtsets(idtset)%nstep=30
   dtsets(idtset)%ntime=1
   dtsets(idtset)%nwfshist=0
   if(dtsets(idtset)%usewvl==1 .and. .not. wvl_bigdft) then
     if(dtsets(idtset)%usepaw==1) then
       dtsets(idtset)%nwfshist=4
     else
       dtsets(idtset)%nwfshist=2
     end if
   end if
!  O
   dtsets(idtset)%occopt=1
   dtsets(idtset)%occ_orig(:,:)=zero
   dtsets(idtset)%omegasrdmax=1.0_dp/Ha_eV  ! = 1eV
   dtsets(idtset)%omegasimax=50/Ha_eV
   dtsets(idtset)%optcell=0
   dtsets(idtset)%optforces=2
   if(dtsets(idtset)%usedmft>0) dtsets(idtset)%optforces=0
   dtsets(idtset)%optstress=1
   dtsets(idtset)%optnlxccc=1
   dtsets(idtset)%orbmag=0
   if (dtsets(idtset)%usepaw==0) then
     dtsets(idtset)%ortalg=2
!    dtsets(idtset)%ortalg=999
   else
     dtsets(idtset)%ortalg=-2
!    dtsets(idtset)%ortalg=999
   end if
!  P
   dtsets(idtset)%paral_atom=paral_atom_default
   dtsets(idtset)%pawcpxocc=1
   dtsets(idtset)%pawcross=0
   dtsets(idtset)%pawecutdg=-one
   dtsets(idtset)%pawfatbnd=0
   dtsets(idtset)%pawlcutd=10
   dtsets(idtset)%pawlmix=10
   dtsets(idtset)%pawmixdg=0 ! Will be set to 1 when npfft>1
   dtsets(idtset)%pawnhatxc=1
   dtsets(idtset)%pawntheta=12
   dtsets(idtset)%pawnphi=13
   dtsets(idtset)%pawnzlm=1
   dtsets(idtset)%pawoptmix=0
   dtsets(idtset)%pawoptosc=0
   dtsets(idtset)%pawovlp=5._dp
   dtsets(idtset)%pawprtdos=0
   dtsets(idtset)%pawprtvol=0
   dtsets(idtset)%pawprtwf=0
   dtsets(idtset)%pawprt_k=0
   dtsets(idtset)%pawprt_b=0
   dtsets(idtset)%pawstgylm=1
   dtsets(idtset)%pawsushat=0
   dtsets(idtset)%pawujat=1
   dtsets(idtset)%pawujrad=20.0_dp
   dtsets(idtset)%pawujv=0.1_dp/Ha_eV
   dtsets(idtset)%pawusecp=1
   dtsets(idtset)%pawxcdev=1
   dtsets(idtset)%ph_nqshift = 0
   if(dtsets(idtset)%ph_nqshift > 0)then
     dtsets(idtset)%ph_qshift = zero
   end if
   dtsets(idtset)%pimd_constraint=0
   dtsets(idtset)%pitransform=0
   dtsets(idtset)%ptcharge(:) = zero
   !dtsets(idtset)%plowan_compute=0
   dtsets(idtset)%plowan_bandi=0
   dtsets(idtset)%plowan_bandf=0
   !if(dtsets(idtset)%plowan_compute>0) then
   dtsets(idtset)%plowan_it(:)=0
   dtsets(idtset)%plowan_iatom(:)=0
   dtsets(idtset)%plowan_lcalc(:)=-1
   dtsets(idtset)%plowan_projcalc(:)=0
   dtsets(idtset)%plowan_nbl(:)=0
   !end if
   dtsets(idtset)%plowan_natom=0
   dtsets(idtset)%plowan_nt=0
   dtsets(idtset)%plowan_realspace=0
   dtsets(idtset)%pol(:)=zero
   dtsets(idtset)%polcen(:)=zero
   dtsets(idtset)%posdoppler=0
   dtsets(idtset)%positron=0
   dtsets(idtset)%posnstep=50
   dtsets(idtset)%posocc=one
   dtsets(idtset)%postoldfe=0.000001_dp
   dtsets(idtset)%postoldff=zero
   dtsets(idtset)%ppmodel=1
   dtsets(idtset)%ppmfrq=zero
   dtsets(idtset)%prepalw=0
   dtsets(idtset)%prepanl=0
   dtsets(idtset)%prepgkk=0
   dtsets(idtset)%prtbbb=0
   dtsets(idtset)%prtbltztrp=0
   dtsets(idtset)%prtcif=0
   dtsets(idtset)%prtden=1;if (dtsets(idtset)%nimage>1) dtsets(idtset)%prtden=0
   dtsets(idtset)%prtdensph=1
   dtsets(idtset)%prtdipole=0
   dtsets(idtset)%prtdos=0
   dtsets(idtset)%prtdosm=0
   dtsets(idtset)%prtebands=1;if (dtsets(idtset)%nimage>1) dtsets(idtset)%prtebands=0
   dtsets(idtset)%prtefg=0
   dtsets(idtset)%prtefmas=1
   dtsets(idtset)%prteig=1;if (dtsets(idtset)%nimage>1) dtsets(idtset)%prteig=0
   dtsets(idtset)%prtelf=0
   dtsets(idtset)%prtfc=0
   dtsets(idtset)%prtfull1wf=0
   dtsets(idtset)%prtfsurf=0
   dtsets(idtset)%prtgden=0
   dtsets(idtset)%prtgeo=0
   dtsets(idtset)%prtgkk=0
   dtsets(idtset)%prtkden=0
   dtsets(idtset)%prtkpt = -1
   dtsets(idtset)%prtlden=0
   dtsets(idtset)%prtnabla=0
   dtsets(idtset)%prtnest=0
   dtsets(idtset)%prtphdos=1
   dtsets(idtset)%prtphsurf=0
   dtsets(idtset)%prtposcar=0
   dtsets(idtset)%prtprocar=0
   dtsets(idtset)%prtpot=0
   dtsets(idtset)%prtpsps=0
   dtsets(idtset)%prtspcur=0
   dtsets(idtset)%prtsuscep=0
   dtsets(idtset)%prtstm=0
   dtsets(idtset)%prtvclmb=0
   dtsets(idtset)%prtvdw=0
   dtsets(idtset)%prtvha=0
   dtsets(idtset)%prtvhxc=0
   dtsets(idtset)%prtvxc=0
   dtsets(idtset)%prtvol=0
   dtsets(idtset)%prtvolimg=0
   dtsets(idtset)%prtvpsp=0
   dtsets(idtset)%prtwant=0
   dtsets(idtset)%prtwf=1; if (dtsets(idtset)%nimage>1) dtsets(idtset)%prtwf=0
   !if (dtset%(idtset)%optdriver == RUNL_RESPFN and all(dtsets(:)%optdriver /= RUNL_NONLINEAR) dtsets(idtset)%prtwf = -1
   dtsets(idtset)%prtwf_full=0
   dtsets(idtset)%prtxml = 0
   do ii=1,dtsets(idtset)%natom,1
     dtsets(idtset)%prtatlist(ii)=ii
   end do
   dtsets(idtset)%prt1dm=0
   dtsets(idtset)%pvelmax(:)=one
   dtsets(idtset)%pw_unbal_thresh=40._dp
!  Q
   dtsets(idtset)%qmass(:)=ten
   dtsets(idtset)%qprtrb(1:3)=0
   dtsets(idtset)%qptdm(:,:)=zero
   dtsets(idtset)%quadmom(:) = zero
!  R
   dtsets(idtset)%random_atpos=0
   dtsets(idtset)%ratsm=zero
   if (any(dtsets(idtset)%constraint_kind(1:dtsets(idtset)%ntypat)>0)) dtsets(idtset)%ratsm=0.05_dp
   dtsets(idtset)%ratsph_extra=two
   dtsets(idtset)%recefermi=zero
   dtsets(idtset)%recgratio=1
   dtsets(idtset)%recnpath=500
   dtsets(idtset)%recnrec=10
   dtsets(idtset)%recrcut=zero
   dtsets(idtset)%recptrott=0
   dtsets(idtset)%rectesteg=0
   dtsets(idtset)%rectolden=zero
   dtsets(idtset)%rcut=zero
   dtsets(idtset)%restartxf=0
   dtsets(idtset)%rfasr=0
   dtsets(idtset)%rfatpol(1:2)=1
   dtsets(idtset)%rfddk=0
   dtsets(idtset)%rfdir(1:3)=0
   dtsets(idtset)%rfelfd=0
   dtsets(idtset)%rfmagn=0
   dtsets(idtset)%rfmeth=1
   dtsets(idtset)%rfphon=0
   dtsets(idtset)%rfstrs=0
   dtsets(idtset)%rfuser=0
   dtsets(idtset)%rf2_dkdk=0
   dtsets(idtset)%rf2_dkde=0
   dtsets(idtset)%rf2_pert1_dir(1:3)=0
   dtsets(idtset)%rf2_pert2_dir(1:3)=0
   dtsets(idtset)%rhoqpmix=one
!  S
   dtsets(idtset)%shiftk_orig(:,:)=one
   dtsets(idtset)%signperm=1
   dtsets(idtset)%slabwsrad=zero
   dtsets(idtset)%slk_rankpp=1000
   dtsets(idtset)%smdelta=0
   dtsets(idtset)%spbroad=0.1_dp
   dtsets(idtset)%spgaxor = -1
   dtsets(idtset)%spgorig = -1
   dtsets(idtset)%spinmagntarget=-99.99_dp
   dtsets(idtset)%spmeth=0
   dtsets(idtset)%spnorbscl=one
   dtsets(idtset)%stmbias=zero
   dtsets(idtset)%strfact=100.0_dp
   dtsets(idtset)%string_algo=1
   dtsets(idtset)%strprecon=one
   dtsets(idtset)%strtarget(1:6)=zero
   dtsets(idtset)%symchi=1
   dtsets(idtset)%symsigma=1
!  T
   dtsets(idtset)%td_maxene=zero
   dtsets(idtset)%td_mexcit=0
   dtsets(idtset)%tfw_toldfe=0.000001_dp
   dtsets(idtset)%tim1rev = 1
   dtsets(idtset)%tl_nprccg = 30
   dtsets(idtset)%tl_radius = zero
   dtsets(idtset)%tphysel=zero
   dtsets(idtset)%toldfe=zero
   dtsets(idtset)%tolmxde=zero
   dtsets(idtset)%toldff=zero
   dtsets(idtset)%tolimg=5.0d-5
   dtsets(idtset)%tolrde=0.005_dp
   dtsets(idtset)%tolrff=zero
   dtsets(idtset)%tolmxf=5.0d-5
   dtsets(idtset)%tolvrs=zero
   dtsets(idtset)%tolwfr=zero
   dtsets(idtset)%tmesh=[5._dp, 59._dp, 6._dp]
   dtsets(idtset)%tsmear=0.01_dp
!  U
   dtsets(idtset)%ucrpa_bands(:)=-1
   dtsets(idtset)%ucrpa_window(:)=-1.0_dp
   dtsets(idtset)%upawu(:,:)=zero
   dtsets(idtset)%usepead=1
   dtsets(idtset)%usefock=0
   dtsets(idtset)%usekden=0
   dtsets(idtset)%use_gemm_nonlop=0
   dtsets(idtset)%use_nonscf_gkk=0 !1 ! deactivate by default, for now 6 Oct 2013
   dtsets(idtset)%userec=0
   dtsets(idtset)%usexcnhat_orig=-1
   dtsets(idtset)%useylm=0
!  V
   dtsets(idtset)%vacnum = -1
   dtsets(idtset)%vcutgeo(:)=zero
   dtsets(idtset)%vdw_nfrag = 1
#if defined DEV_YP_VDWXC
   dtsets(idtset)%vdw_df_acutmin = vdw_defaults%acutmin
   dtsets(idtset)%vdw_df_aratio = vdw_defaults%aratio
   dtsets(idtset)%vdw_df_damax = vdw_defaults%damax
   dtsets(idtset)%vdw_df_damin = vdw_defaults%damin
   dtsets(idtset)%vdw_df_dcut = vdw_defaults%dcut
   dtsets(idtset)%vdw_df_dratio = vdw_defaults%dratio
   dtsets(idtset)%vdw_df_dsoft = vdw_defaults%dsoft
   dtsets(idtset)%vdw_df_gcut = vdw_defaults%gcut
   dtsets(idtset)%vdw_df_ndpts = vdw_defaults%ndpts
   dtsets(idtset)%vdw_df_ngpts = vdw_defaults%ngpts
   dtsets(idtset)%vdw_df_nqpts = vdw_defaults%nqpts
   dtsets(idtset)%vdw_df_nrpts = vdw_defaults%nrpts
   dtsets(idtset)%vdw_df_nsmooth = vdw_defaults%nsmooth
   dtsets(idtset)%vdw_df_phisoft = vdw_defaults%phisoft
   dtsets(idtset)%vdw_df_qcut = vdw_defaults%qcut
   dtsets(idtset)%vdw_df_qratio = vdw_defaults%qratio
   dtsets(idtset)%vdw_df_rcut = vdw_defaults%rcut
   dtsets(idtset)%vdw_df_rsoft = vdw_defaults%rsoft
   dtsets(idtset)%vdw_df_tolerance = vdw_defaults%tolerance
   dtsets(idtset)%vdw_df_tweaks = vdw_defaults%tweaks
   dtsets(idtset)%vdw_df_zab = vdw_defaults%zab
   dtsets(idtset)%vdw_df_threshold = 1.0d-2
#endif
   dtsets(idtset)%vdw_supercell(:) = 0
   dtsets(idtset)%vdw_tol = tol10
   dtsets(idtset)%vdw_tol_3bt = -1
   dtsets(idtset)%vdw_typfrag(:) = 1
   dtsets(idtset)%vdw_xc = 0
   dtsets(idtset)%vis=100.0_dp
   dtsets(idtset)%vprtrb(1:2)=zero
!  W
   dtsets(idtset)%wtatcon(:,:,:)=zero
   dtsets(idtset)%wfmix=one
   dtsets(idtset)%wfk_task=0
   dtsets(idtset)%wtk=one
   dtsets(idtset)%wvl_crmult  = 6._dp
   dtsets(idtset)%wvl_frmult  = 10._dp
   dtsets(idtset)%wvl_hgrid   = 0.5_dp
   dtsets(idtset)%wvl_ngauss  =(/1,100/)
   dtsets(idtset)%wvl_nprccg  = 10
   dtsets(idtset)%w90iniprj   = 1
   dtsets(idtset)%w90prtunk   = 0

!  X
   dtsets(idtset)%xclevel  = 0
   dtsets(idtset)%xc_denpos = tol14
   dtsets(idtset)%xc_tb09_c = 99.99_dp
   dtsets(idtset)%xredsph_extra(:,:)=zero
!  Y
!  Z
   dtsets(idtset)%zcut=3.67493260d-03  ! = 0.1eV
   if(dtsets(idtset)%optdriver == RUNL_GWLS) dtsets(idtset)%zcut=zero
   !if(dtsets(idtset)%optdriver == RUNL_EPH) dtsets(idtset)%zcut = 0.01 * eV_Ha
   dtsets(idtset)%ziontypat(:)=zero

!  BEGIN VARIABLES FOR @Bethe-Salpeter
   dtsets(idtset)%bs_algorithm    =2
   dtsets(idtset)%bs_haydock_niter=100
   dtsets(idtset)%bs_exchange_term=1
   dtsets(idtset)%bs_coulomb_term=11
   dtsets(idtset)%bs_calctype=1
   dtsets(idtset)%bs_coupling=0

   dtsets(idtset)%bs_haydock_tol=(/0.02_dp,zero/)

   dtsets(idtset)%bs_loband=0
!  Take big absolute value numbers, but the the biggest ones, otherwise overflow can happen
   dtsets(idtset)%bs_eh_cutoff = [smallest_real*tol6,greatest_real*tol6]
   dtsets(idtset)%bs_freq_mesh = [zero,zero,0.01_dp/Ha_eV]

!  Interpolation
   dtsets(idtset)%bs_interp_method = 1 ! YG interpolation
   dtsets(idtset)%bs_interp_mode = 0 ! No interpolation
   dtsets(idtset)%bs_interp_prep = 0 ! Do not prepare interp
   dtsets(idtset)%bs_interp_kmult = 0
   dtsets(idtset)%bs_interp_m3_width = one
   dtsets(idtset)%bs_interp_rl_nb = 1

!  END VARIABLES FOR @Bethe-Salpeter.

! EPH variables
   dtsets(idtset)%asr = 1
   dtsets(idtset)%dipdip = 1
   dtsets(idtset)%chneut = 1
   dtsets(idtset)%symdynmat = 1

   dtsets(idtset)%ph_ndivsm = 20
   dtsets(idtset)%ph_nqpath = 0
   dtsets(idtset)%ph_ngqpt = [20, 20, 20]

   dtsets(idtset)%eph_mustar = 0.1_dp
   dtsets(idtset)%eph_intmeth = 2
   dtsets(idtset)%eph_extrael = zero
   dtsets(idtset)%eph_fermie = zero
   dtsets(idtset)%eph_frohlichm = 0
   dtsets(idtset)%eph_fsmear = 0.01
   dtsets(idtset)%eph_fsewin = 0.04
   dtsets(idtset)%eph_ngqpt_fine = [0, 0, 0]
   dtsets(idtset)%eph_task = 1
   dtsets(idtset)%eph_transport  = 0

   dtsets(idtset)%ph_wstep = 0.1/Ha_meV
   dtsets(idtset)%ph_intmeth = 2
   dtsets(idtset)%ph_nqshift = 1
   dtsets(idtset)%ph_smear = 0.00002_dp
   dtsets(idtset)%ddb_ngqpt = [0, 0, 0]
   dtsets(idtset)%dvdb_ngqpt = [0, 0, 0]
   dtsets(idtset)%ddb_shiftq(:) = zero

! JB:UNINITIALIZED VALUES (not found in this file neither indefo1)
! They might be initialized somewhereelse, I don't know.
! That might cause unitialized error with valgrind depending on the compilo
! chkprim
! maxnsym
! nsym
! macro_uj
! prtpmp
! timopt
! useria
! userib
! useric
! userid
! userie
! bravais
! symafm
! symrel
! fband
! nelect
! userra
! userrb
! userrc
! userrd
! userre
! vacwidth
! genafm
! kptns
! rprimd_orig
! tnons

 end do

 DBG_EXIT("COLL")

end subroutine indefo
!!***

end module m_invars1
!!***
