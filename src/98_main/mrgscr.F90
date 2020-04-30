!!****p* ABINIT/mrgscr
!! NAME
!! mrgscr
!!
!! FUNCTION
!! This code reads partial (SCR|SUSC) files for different q points creating a single file that
!! can be used to perform a sigma calculation.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2020 ABINIT group (RS, MG, MS)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! If the number of SCR files to be merged is equal to 1, the program checks
!! the integrity of the file reporting the list of missing q-points.
!! Note that the list of required q-points depends on the k-mesh
!! used during the calculation of the KSS file. We assume indeed that the same k-mesh
!! is used during the calculation of the matrix elements of sigma.
!!
!! INPUTS
!!  (Main program)
!!
!! OUTPUT
!!  Only checking and writing
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,abinit_doctor,cqratio,crystal_free
!!      crystal_from_hdr,cspint,decompose_epsm1,destroy_mpi_enreg
!!      em1results_free,em1results_print,find_qmesh,flush_unit,fourdp
!!      get_ppm_eigenvalues,getem1_from_ppm_one_ggp,getng,gsph_free,gsph_init
!!      hdr_free,herald,hscr_free,hscr_from_file,hscr_print,init_er_from_file
!!      initmpi_seq,int2char4,ioscr_qmerge,ioscr_qrecover,ioscr_wmerge
!!      ioscr_wremove,kmesh_free,kmesh_init,kmesh_print,metric,mkdump_er
!!      pawrhoij_free,ppm_free,ppm_init,prompt,read_rhor,read_screening
!!      remove_phase,setup_ppmodel,test_charge,timein,vcoul_free,vcoul_init
!!      wrtout,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program mrgscr

 use defs_basis
 use m_xmpi
 use m_abicore
 use m_build_info
 use m_errors
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr
 use m_crystal
 use m_pawrhoij
 use m_dtset

 use defs_abitypes,         only : MPI_type
 use m_specialmsg,          only : herald
 use m_time,                only : timein
 use m_gwdefs,              only : GW_TOLQ, GW_TOLQ0, GW_Q0_DEFAULT
 use m_io_tools,            only : prompt, flush_unit, open_file
 use m_fstrings,            only : int2char4, endswith, itoa, sjoin
 use m_fft_mesh,            only : g2ifft
 use m_fftcore,             only : get_cache_kb, getng
 use m_fft,                 only : fourdp
 use m_numeric_tools,       only : iseven, cspint
 use m_mpinfo,              only : destroy_mpi_enreg, initmpi_seq
 use m_geometry,            only : normv, metric
 use m_gsphere,             only : gsph_init, gsph_free, gsphere_t
 use m_bz_mesh,             only : kmesh_t, find_qmesh, kmesh_init, kmesh_print, kmesh_free
 use m_vcoul,               only : vcoul_t, vcoul_init, vcoul_free
 use m_ioarr,               only : read_rhor
 use m_io_screening,        only : hscr_print, read_screening, hscr_free, hscr_t, hscr_from_file,&
                                   ioscr_qmerge, ioscr_qrecover, ioscr_wmerge, ioscr_wremove
 use m_ppmodel,             only : ppm_init, ppm_free, setup_ppmodel, getem1_from_PPm_one_ggp, &
                                   get_PPm_eigenvalues, ppmodel_t, cqratio
 use m_model_screening,     only : remove_phase
 use m_screening,           only : mkdump_er, em1results_free, em1results_print, decompose_epsm1, &
                                   init_er_from_file, Epsilonm1_results
 use m_wfd,                 only : test_charge

 implicit none

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,paral_kgb0=0,rdwr2=2,prtvol=0,cplex1=1
 integer :: iomode,fform1,ifile,ierr,ii,ios,iqibz,iqf,nfiles,timrev
 integer :: unt_dump,idx,ig1,ig2,iomega,ppmodel,npwe_asked,mqmem,io,unt_dump2
 integer :: id_required,ikxc,approx_type,option_test,dim_kxcg,usexcnhat,usefinegrid
 integer :: mgfft,nqlwl,nfft,igmax,comm,nq_selected,kptopt
 integer :: choice,nfreq_tot,nfreqre,nfreqim,nfreqc,ifrq,imax
 integer :: ig1_start,ig1_end,ig2_start,ig2_end,gmgp_idx,orig_npwe
 real(dp) :: ucvol,boxcutmin,ecut,drude_plsmf,compch_fft,compch_sph
 real(dp) :: nelectron_exp,freqremax,eps_diff,eps_norm,eps_ppm_norm
 real(dp) :: value1,value2,factor,GN_drude_plsmf
 real(dp) :: tcpu,tcpui,twall,twalli
 real(gwp) :: phase
 logical :: is_sus,is_scr,same_freqs,calc_epsilon,only_diag
 character(len=1) :: ans
 character(len=10) :: tagq
 character(len=24) :: codename
 character(len=500) :: msg
 character(len=nctk_slen) :: varname
 character(len=fnlen) :: fname_out,fname,fname_dump,fname_rho,prefix,fname_eigen,fname_dump2
 type(hdr_type) :: hdr_rhor
 type(abifile_t) :: abifile
 type(hscr_t),pointer :: Hscr0
 type(hscr_t),target :: Hscr_merge
 type(MPI_type) :: MPI_enreg
 type(kmesh_t) :: Kmesh,Qmesh
 type(crystal_t) :: Cryst
 type(gsphere_t)  :: Gsphere
 type(ppmodel_t) :: PPm
 type(Epsilonm1_results) :: Er
 type(vcoul_t),target :: Vcp
 type(Dataset_type) :: Dtset
!arrays
 integer :: ngfft(18)
 integer,allocatable :: foundq(:),freq_indx(:,:)
 real(dp),parameter :: k0(3) = [zero,zero,zero]
 real(dp) :: gmet(3,3),gprimd(3,3),qdiff(3),rmet(3,3),qtmp(3),tsec(2)
 real(dp),allocatable :: qlwl(:,:),real_omega(:),rhor(:,:),rhog(:,:),nhat(:,:)
 real(dp),allocatable :: work(:),ftab(:),ysp(:,:),eint(:),qratio(:,:)
 complex(gwpc),pointer :: vc_sqrt(:)
 complex(gwpc),allocatable :: epsm1(:,:,:,:),kxcg(:,:)
 complex(dpc),allocatable :: omega(:),em1_ppm(:),epsm1_eigen(:,:),ppm_eigen(:,:),rhoggp(:,:)
 character(len=fnlen),allocatable :: filenames(:)
 type(pawrhoij_type),allocatable :: pawrhoij(:)
 type(hscr_t),target,allocatable :: Hscr_file(:) ! Cannot use allocatable as pathscale5 miscompiles the code

! *************************************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 ! Initialize MPI
 call xmpi_init()

 ! Initialize memory profiling if it is activated
 ! if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
 ! note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 call timein(tcpui,twalli)

 ! Default for sequential use
 call initmpi_seq(MPI_enreg); comm = MPI_enreg%comm_world

 is_sus=.FALSE.; is_scr=.FALSE.

 ! Write greetings, and read the number of files ===
 codename='MRGSCR'//REPEAT(' ',18)
 call herald(codename,abinit_version,std_out)

 call prompt(' Enter the number of files to merge: ',nfiles)
 ABI_CHECK(nfiles > 0, 'nfiles must be >0')

 ABI_MALLOC(filenames,(nfiles))
 ABI_MALLOC(Hscr_file,(nfiles))

 if (nfiles == 1) then
   call prompt(' Enter the name of the file to be analyzed: ',filenames(1))
   write(msg,'(7a)')ch10,&
    ' Running single-file mode:',ch10,&
    ' Checking the integrity of file: ',TRIM(filenames(1)),ch10,&
    ' reporting the list of q-points that are missing. '
   call wrtout(std_out,msg,'COLL')

   if (nctk_try_fort_or_ncfile(filenames(1), msg) /= 0) then
     MSG_ERROR(msg)
   end if

 else if (nfiles > 1) then
   ! Read name of files to be merged and check for existence.
   call prompt(' Enter the prefix for the final output file: ',fname_out)

   do ifile=1,nfiles
     write(msg,'(a,i4)')' Enter the name for the partial screening file no.',ifile
     call prompt(msg,filenames(ifile))

     if (nctk_try_fort_or_ncfile(filenames(ifile), msg) /= 0) then
       MSG_ERROR(msg)
     end if

   end do
 end if

 ! Read the header of each file.
 do ifile=1,nfiles
   iomode = IO_MODE_FORTRAN; if (endswith(filenames(ifile), ".nc")) iomode = IO_MODE_ETSF

   call hscr_from_file(Hscr_file(ifile), filenames(ifile), fform1, comm)
   ABI_CHECK(fform1 /= 0, sjoin("fform == 0 in", filenames(ifile)))

   abifile = abifile_from_fform(fform1)
   if (abifile%fform == 0) then
     MSG_ERROR(sjoin("Cannot find any abifile object associated to fform1:", itoa(fform1)))
   end if
   if (abifile%class /= "polariz" .and. abifile%class /= "epsm1") then
     MSG_ERROR(sjoin('Error while reading header, fform= ',itoa(fform1)))
   end if
   is_scr = abifile%class == "epsm1"
   is_sus = abifile%class == "polariz"

   call hscr_print(Hscr_file(ifile),unit=std_out,prtvol=1)

   if (ifile == 1) call metric(gmet,gprimd,-1,rmet,Hscr_file(ifile)%Hdr%rprimd,ucvol)
 end do !ifile

 if (nfiles > 1) then
   ! Put the correct ending on the output file
   if (is_scr) fname_out=TRIM(fname_out)//'_SCR'
   if (is_sus) fname_out=TRIM(fname_out)//'_SUS'
 end if

 ! Produce output file in netcdf format we are merging netcdf files.
 if (iomode == IO_MODE_ETSF .and. .not. endswith(fname_out, ".nc")) fname_out = nctk_ncify(fname_out)

 !============================
 !=== Merge multiple files ===
 !============================
 if (nfiles > 1) then

   ! Check what kind of merging is to be performed
   write(std_out,'(2(a))') ch10,' Do you want to merge q-points        (= 1) ?'
   write(std_out,'(a)')         '  or do you want to merge frequencies (= 2) ?'
   read(std_in,*)choice

   select case(choice)
   case(1)
     write(std_out,'(3a)') ch10,' 1 => merging q-points',ch10
     call ioscr_qmerge(nfiles, filenames, hscr_file, fname_out, hscr_merge)

   case (2)
     ! Merge frequencies
     write(std_out,'(3a)') ch10,' 2 => merging frequency grids',ch10
     !MSG_WARNING("Advanced user option, consistency in fform etc. will not be checked.")

     write(std_out,'(2a)') ch10,' Enter freqremax [eV] for the merged file (Enter 0 to use all freq. found):'
     read(std_in,*)freqremax

     freqremax = freqremax/Ha_eV; if (freqremax<tol16) freqremax = HUGE(freqremax)
     call ioscr_wmerge(nfiles, filenames, hscr_file, freqremax, fname_out, hscr_merge)

   case default
     MSG_ERROR("Invalid choice!")
   end select

 end if ! nfiles>1

 ! Now check if the list of q-points is complete
 ! Here we assume that the k-mesh reported in the header is the same as that used during the sigma calculation.
 write(msg,'(3a)') ch10,' Checking if the list of q-points is complete. ',ch10
 call wrtout(std_out,msg,'COLL')

 !call hscr_check_qpoints(hscr0)

 Hscr0 => Hscr_file(1)
 fname =filenames(1)
 if (nfiles>1) then
   Hscr0 => Hscr_merge
   fname = fname_out
 end if

 timrev=2 ! This should be read from kptopt
 cryst = HScr0%Hdr%get_crystal(timrev,remove_inv=.FALSE.)

 kptopt=1
 call kmesh_init(Kmesh,Cryst,HScr0%Hdr%nkpt,Hscr0%Hdr%kptns,kptopt)
 call kmesh_print(Kmesh,"K-mesh for the wavefunctions",prtvol=prtvol)

 call find_qmesh(Qmesh,Cryst,Kmesh)
 call kmesh_print(Qmesh,"Q-mesh for the screening function",prtvol=prtvol)

 ABI_MALLOC(foundq,(Qmesh%nibz))
 foundq(:)=0
 do iqibz=1,Qmesh%nibz
   do iqf=1,Hscr0%nqibz
     qdiff(:)=Qmesh%ibz(:,iqibz)-Hscr0%qibz(:,iqf)
     if (normv(qdiff,gmet,'G')<GW_TOLQ) foundq(iqibz)=foundq(iqibz)+1
   end do
 end do

 if (ANY(foundq==0)) then
   write(msg,'(6a)')ch10,&
    ' File ',TRIM(fname),' is not complete ',ch10,&
    ' The following q-points are missing:'
   call wrtout(std_out,msg,'COLL')
   ii=0
   do iqibz=1,Qmesh%nibz
     if (foundq(iqibz)==0) then
       ii=ii+1
       write(msg,'(i3,a,3f12.6)')ii,') ',Qmesh%ibz(:,iqibz)
       call wrtout(std_out,msg,'COLL')
     end if
   end do
 end if

 if (ANY(foundq>1)) then
   write(msg,'(6a)')ch10,&
    ' File ',TRIM(fname),' is overcomplete ',ch10,&
    ' The following q-points are present more than once:'
   call wrtout(std_out,msg,'COLL')
   ii=0
   do iqibz=1,Qmesh%nibz
     if (foundq(iqibz)>1) then
       ii=ii+1
       write(msg,'(i3,a,3f12.6)')ii,') ',Qmesh%ibz(:,iqibz)
       call wrtout(std_out,msg,'COLL')
     end if
   end do
 end if

 if (ALL(foundq==1)) then
   write(msg,'(5a)')ch10,&
    '.File ',TRIM(fname),' contains a complete list of q-points ',ch10
   call wrtout(std_out,msg,'COLL')
 end if

 !=====================
 !=== Recovery mode ===
 !=====================
 if (nfiles==1) then

   write(std_out,'(2(a))') ch10,' Do you want to recover a subset of q-points    (= 1) ?'
   write(std_out,'(a)')         '  or extract the contents of the file           (= 2) ?'
   write(std_out,'(a)')         '  or create dielectric function (SCR file)'
   write(std_out,'(a)')         '    and/or extract plasmon-pole parameters      (= 3) ?'
   write(std_out,'(a)')         '  or remove real frequencies                    (= 4) ?'
   write(std_out,'(a)')         '  or remove imaginary frequencies               (= 5) ?'
   write(std_out,'(a)')         '  or calculate a model screening                (= 6) ?'
   !write(std_out,'(a)')         '  or interpolate the screening in k-space       (= 8) ?'
   !write(std_out,'(a)')         '  or convert a netcd file to Fortran            (= 9) ?'
   read(std_in,*)choice

   select case(choice)
   case(1)
       ! Recover subset of q-points --------------------------------------------------
     write(std_out,'(a)') ' 1 => Recovering subset of q-points'
     call prompt(' Enter the number of q-points to be extracted: ',nq_selected)
     call prompt(' Enter the name of the final output file: ',fname_out)

     if (endswith(filenames(1), ".nc") .and. .not. endswith(fname_out, ".nc")) then
       fname_out = nctk_ncify(fname_out)
       call wrtout(std_out,"- Added .nc extension to output file as input data is in netcdf format.")
     end if

     call ioscr_qrecover(filenames(1), nq_selected, fname_out)

   case(2)
     ! Analyse file ----------------------------------------------------------------
     ABI_CHECK(iomode==IO_MODE_FORTRAN, "netcdf output not coded")

     write(std_out,'(a)') ' 2 => Extraction of file contents'

     ! Initialize the G-sphere.
     call gsph_init(Gsphere,Cryst,Hscr0%npwe,gvec=Hscr0%gvec)

     ABI_MALLOC_OR_DIE(epsm1,(Hscr0%npwe,Hscr0%npwe,Hscr0%nomega,1), ierr)

     ! Give option to output epsilon instead of chi0
     calc_epsilon = .FALSE.
     if (is_sus) then
       write(std_out,'(2a)') ch10,' You have provided a chi_0 file for analysis. Would you like to output'
       write(std_out,'(2a)',advance='no') ' the dielectric function epsilon_GG'' ', '= delta_GG'' - v_G*chi0_GG''[Y/N] ? '
       read(std_in,*)ans

       if (ans=='Y'.or.ans=='y') then
         ! Initialise Coulomb terms
         if (Er%Hscr%nqlwl==0) then
           nqlwl=1
           ABI_MALLOC(qlwl,(3,nqlwl))
           qlwl(:,1)= GW_Q0_DEFAULT
         else
           nqlwl=Er%Hscr%nqlwl
           ABI_MALLOC(qlwl,(3,nqlwl))
           qlwl(:,:)=Er%Hscr%qlwl(:,1:nqlwl)
         end if

         Dtset%icsing=3; Dtset%rcut=zero
         Dtset%vcutgeo=(/zero,zero,zero/);
         Dtset%boxcenter=(/zero,zero,zero/)

         write(std_out,'(2a)',advance='no') ch10,' Was a Coulomb cutoff technique used [Y/N] ? '
         read(std_in,*)ans
         if (ans=='Y'.or.ans=='y') then
           write(std_out,'(2a)',advance='no') ' Enter icsing: '
           read(std_in,*)Dtset%icsing
           write(std_out,'(2a)',advance='no') ' Enter vcutgeo: '
           read(std_in,*)Dtset%vcutgeo
           write(std_out,'(2a)',advance='no') ' Enter boxcenter: '
           read(std_in,*)Dtset%boxcenter
         end if
         dtset%ecutsigx = -one

         call vcoul_init(Vcp,Gsphere,Cryst,Qmesh,Kmesh,Dtset%rcut,Dtset%icsing,&
           Dtset%vcutgeo,Dtset%ecutsigx,Hscr0%npwe,nqlwl,qlwl,ngfft,comm)
         ABI_FREE(qlwl)

         calc_epsilon = .TRUE.
       end if
     end if

     ig1 = 0; ig2 = 0
     write(std_out,'(2(a),I0,a)',advance='NO') ch10,' Enter the starting index for G (1 - ',Hscr0%npwe,' ): '
     read(std_in,*)ig1_start
     if (ig1_start<1.OR.ig1_start>Hscr0%npwe) then
       MSG_ERROR(' Starting index out of bounds')
     end if
     write(std_out,'(a,I0,a,I0,a)',advance='NO')    ' Enter the ending index for G ( ',ig1_start,' - ',Hscr0%npwe,' ): '
     read(std_in,*)ig1_end
     if (ig1_end<ig1_start.OR.ig1_end>Hscr0%npwe) then
       MSG_ERROR(' Ending index out of bounds')
     end if
     write(std_out,'(a,I0,a)',advance='NO')         ' Enter the starting index for G'' (1 - ',Hscr0%npwe,' ): '
     read(std_in,*)ig2_start
     if (ig2_start<1.OR.ig2_start>Hscr0%npwe) then
       MSG_ERROR(' Starting index out of bounds')
     end if
     write(std_out,'(a,I0,a,I0,a)',advance='NO')    ' Enter the ending index for G'' ( ',ig2_start,' - ',Hscr0%npwe,' ): '
     read(std_in,*)ig2_end
     if (ig2_end<ig2_start.OR.ig2_end>Hscr0%npwe) then
       MSG_ERROR(' Ending index out of bounds')
     end if

     only_diag = .FALSE.
     write(std_out,'(a)',advance='no') ' Would you like to output only the diagonal [Y/N] ? '
     read(std_in,*)ans
     if (ans=='Y'.or.ans=='y') only_diag = .TRUE.

     do iqibz=1,Hscr0%nqibz
       ! In the long wavelength limit we set q==0, because we still can use symmetries for the Body.
       qtmp(:)=Hscr0%qibz(:,iqibz); if (normv(qtmp,Cryst%gmet,'G')<GW_TOLQ0) qtmp(:)=zero

       ! FIXME
       varname = "none"
       call read_screening(varname,fname,Hscr0%npwe,1,Hscr0%nomega,epsm1,iomode,comm,iqiA=iqibz)

       if (calc_epsilon) then ! Calculate epsilon
         do iomega=1,Hscr0%nomega
           if (iqibz==1) then
             if (nqlwl>1) then
               MSG_ERROR('nqlwl>1 not coded yet!')
             end if
             vc_sqrt => Vcp%vcqlwl_sqrt(:,iqibz)  ! Use Coulomb term for q-->0
           else
             vc_sqrt => Vcp%vc_sqrt(:,iqibz)
           end if
           do ig2=ig2_start,ig2_end
             do ig1=ig1_start,ig1_end
               epsm1(ig1,ig2,iomega,1) = -(vc_sqrt(ig1)**2)*epsm1(ig1,ig2,iomega,1)
             end do ! ig1
             epsm1(ig2,ig2,iomega,1) = one + epsm1(ig2,ig2,iomega,1)
           end do ! ig2
         end do ! iomega
       end if ! Do we calculate epsilon

       ! Find out the total number of frequencies along real/imaginary axes
       ! and possibly in the z-plane
       nfreqre=0; nfreqim=0; nfreqc=0;
       do iomega=1,Hscr0%nomega
         if (ABS(REAL(Hscr0%omega(iomega)))<tol8.AND. ABS(AIMAG(Hscr0%omega(iomega)))<tol8) nfreqre = nfreqre + 1
         if (ABS(REAL(Hscr0%omega(iomega)))>tol8.AND. ABS(AIMAG(Hscr0%omega(iomega)))<tol8) nfreqre = nfreqre + 1
         if (ABS(REAL(Hscr0%omega(iomega)))<tol8.AND. ABS(AIMAG(Hscr0%omega(iomega)))>tol8) nfreqim = nfreqim + 1
       end do
       if (Hscr0%nomega-nfreqre-nfreqim/=0) then
         write(std_out,'(/,a)') ' WARNING: There are frequencies in the full complex plane.'
         write(std_out,'(a)')   '          The _SCR or _SUS file might not be suitable'
         write(std_out,'(a,/)') '          for self-energy calculations.'
         nfreqc = Hscr0%nomega-nfreqre-nfreqim
       end if
       write(std_out,'(2a,I0,a)') ch10,' Found ',Hscr0%nomega,' frequencies.'
       write(std_out,'(2(a,I0),2a)') ' ',nfreqre,' real, and ',nfreqim,' imaginary.',ch10
       if (nfreqc>0) then
         write(std_out,'(a,I0)') ' There is a grid in the complex plane with ',nfreqc
         write(std_out,'(2a)')   '  extra frequencies in the list.',ch10
       end if

       ! Get Q index for name
       call int2char4(iqibz,tagq)
       ABI_CHECK((tagq(1:1)/='#'),'Bug: string length too short!')

       if (nfreqre>0) then ! Output real frequency axis
         if (calc_epsilon) then
           fname_dump=TRIM(fname)//'_EPS_Q'//TRIM(tagq)
         else
           fname_dump=TRIM(fname)//'_Q'//TRIM(tagq)
         end if

         if (open_file(fname_dump, msg, newunit=unt_dump, status='replace', form='formatted') /= 0) then
           MSG_ERROR(msg)
         end if

         do ig1=ig1_start,ig1_end
           do ig2=ig2_start,ig2_end
             if (only_diag.AND.ig1/=ig2) CYCLE
             write(unt_dump,'(2(a,i8),/,a,3f12.6,/,a,3i6,a,3i6,/,a,/)')&
             '# ig1= ',ig1,'    ig2= ',ig2,&
             '# q = ',Hscr0%qibz(:,iqibz),&
             '# G = ',Hscr0%gvec(:,ig1),'  G''= ',Hscr0%gvec(:,ig2),&
             '#   omega [eV]           Re             Im '
             do iomega=1,nfreqre
               write(unt_dump,'(f8.2,4x,2es16.8)') REAL(Hscr0%omega(iomega))*Ha_eV,&
                 REAL(epsm1(ig1,ig2,iomega,1)),AIMAG(epsm1(ig1,ig2,iomega,1))
             end do
             write(unt_dump,*)
             write(unt_dump,*)
           end do !ig2
         end do !ig1
         close(unt_dump)
       end if ! Output real frequency axis

       if (nfreqim>0) then ! output imaginary frequency axis
         if (calc_epsilon) then
           fname_dump=TRIM(fname)//'_EPS_Imfrq_Q'//TRIM(tagq)
         else
           fname_dump=TRIM(fname)//'_Imfrq_Q'//TRIM(tagq)
         end if
         if (open_file(fname_dump,msg,newunit=unt_dump,status='replace',form='formatted') /= 0) then
           MSG_ERROR(msg)
         end if
         do ig1=ig1_start,ig1_end
           do ig2=ig2_start,ig2_end
             if (only_diag.AND.ig1/=ig2) CYCLE
             write(unt_dump,'(a,i4,2(a,i8),/,a,3f12.6,/,a,3i6,a,3i6,/,a,/)')&
               '# index= ',idx,'    ig1= ',ig1,'    ig2= ',ig2,&
               '# q = ',Hscr0%qibz(:,iqibz),&
               '# G = ',Hscr0%gvec(:,ig1),'  G''= ',Hscr0%gvec(:,ig2),&
               '#   omega [eV]           Re             Im '
             do iomega=nfreqre+1,nfreqre+nfreqim
               write(unt_dump,'(f8.2,4x,2es16.8)') AIMAG(Hscr0%omega(iomega))*Ha_eV,epsm1(ig1,ig2,iomega,1)
             end do
             write(unt_dump,*)
             write(unt_dump,*)
           end do !ig2
         end do !ig1
         close(unt_dump)
       end if ! Check for imaginary frequencies

       ! Check for complex plane values
       if (nfreqc>0) then
         if (calc_epsilon) then
           fname_dump=TRIM(fname)//'_EPS_ZPLANE_Q'//TRIM(tagq)
         else
           fname_dump=TRIM(fname)//'_ZPLANE_Q'//TRIM(tagq)
         end if

         if (open_file(fname_dump,msg,newunit=unt_dump,status='replace',form='formatted') /= 0) then
           MSG_ERROR(msg)
         end if

         do ig1=ig1_start,ig1_end
           do ig2=ig2_start,ig2_end
             if (only_diag.AND.ig1/=ig2) CYCLE
             write(unt_dump,'(a,i4,2(a,i8),/,a,3f12.6,/,a,3i6,a,3i6,/,a,/)')&
              '# index= ',idx,'    ig1= ',ig1,'    ig2= ',ig2,&
              '# q = ',Hscr0%qibz(:,iqibz),&
              '# G = ',Hscr0%gvec(:,ig1),'  G''= ',Hscr0%gvec(:,ig2),&
              '#   omega [eV]           Re             Im '
             do iomega=1,nfreqre
               write(unt_dump,'(2(f8.2),4x,2es16.8)') REAL(Hscr0%omega(iomega))*Ha_eV,&
                 AIMAG(Hscr0%omega(iomega))*Ha_eV,epsm1(ig1,ig2,iomega,1)
             end do
             write(unt_dump,*)
             do ios=1,nfreqim
               do iomega=1,nfreqre
                 if (iomega==1) then
                   io = nfreqre + ios
                 else
                   io = nfreqre + nfreqim + (ios-1)*(nfreqre-1) + (iomega-1)
                 end if
                 write(unt_dump,'(2(f8.2),4x,2es16.8)') REAL(Hscr0%omega(io))*Ha_eV,&
                   AIMAG(Hscr0%omega(io))*Ha_eV,epsm1(ig1,ig2,io,1)
               end do
               write(unt_dump,*)
             end do
             write(unt_dump,*)
             write(unt_dump,*)
           end do !ig2
         end do !ig1
         close(unt_dump)
       end if ! Check for complex plane freqs

     end do !iqibz

     ABI_FREE(epsm1)
     call gsph_free(Gsphere)

   case(3)
       ! Extract dielectric function and plasmon-pole stuff --------------------------
     ABI_CHECK(iomode==IO_MODE_FORTRAN, "netcdf output not coded")
     write(std_out,'(a)') ' 3 => Calculation of dielectric function and plasmon-pole model'

     npwe_asked=Hscr0%npwe; mqmem=Hscr0%nqibz
     call init_Er_from_file(Er,fname,mqmem,npwe_asked,comm)

     ! Initialize the G-sphere ===
     call gsph_init(Gsphere,Cryst,Hscr0%npwe,gvec=Hscr0%gvec)

     boxcutmin=two; igmax=Gsphere%shlim(Gsphere%nsh)
     ecut=Er%Hscr%Hdr%ecutdg

     call getng(boxcutmin,ecut,Gsphere%gmet,k0,MPI_enreg%me_fft,&
       mgfft,nfft,ngfft,MPI_enreg%nproc_fft,Cryst%nsym,paral_kgb0,Cryst%symrel)

     ! I am using standard valued, it would be better to call indefo
     ! ngfft(1:3)=Er%Hscr%Hdr%ngfft(1:3)
     ngfft(7)=112
     ngfft(8)=get_cache_kb()
     nfft = PRODUCT(ngfft(1:3))

     Dtset%icsing=3; Dtset%rcut=zero
     Dtset%vcutgeo=(/zero,zero,zero/); Dtset%boxcenter=(/zero,zero,zero/)
     Dtset%ecutsigx = -1

     if (Er%Hscr%nqlwl==0) then
       nqlwl=1
       ABI_MALLOC(qlwl,(3,nqlwl))
       qlwl(:,1)= GW_Q0_DEFAULT
     else
       nqlwl=Er%Hscr%nqlwl
       ABI_MALLOC(qlwl,(3,nqlwl))
       qlwl(:,:)=Er%Hscr%qlwl(:,1:nqlwl)
     end if

     call vcoul_init(Vcp,Gsphere,Cryst,Qmesh,Kmesh,Dtset%rcut,Dtset%icsing,Dtset%vcutgeo,Dtset%ecutsigx,Hscr0%npwe,nqlwl,&
       qlwl,ngfft,comm)
     ABI_FREE(qlwl)

     ! Get the density from an external file ===
     ! If meshes are not the same, do an FFT interpolation to have rhor on ngfft.
     call prompt(' Enter name for external DEN (or PAWDEN) file: ', fname_rho)

     ABI_MALLOC(rhor,(nfft,Hscr0%Hdr%nspden))
     ABI_MALLOC(pawrhoij,(Hscr0%Hdr%natom*Hscr0%Hdr%usepaw))

     call read_rhor(fname_rho, cplex1, nfft, Hscr0%Hdr%nspden, ngfft, 1, MPI_enreg, rhor, hdr_rhor, pawrhoij, comm)

     call hdr_rhor%free()
     call pawrhoij_free(pawrhoij)
     ABI_FREE(pawrhoij)

     ABI_MALLOC(rhog,(2,nfft))
     call fourdp(1,rhog,rhor(:,1),-1,MPI_enreg,nfft,1,ngfft,0)

     ABI_MALLOC(nhat,(nfft,Hscr0%Hdr%nspden*Hscr0%Hdr%usepaw))
     compch_sph=greatest_real; compch_fft=greatest_real
     usexcnhat=0; usefinegrid=0

     nelectron_exp = Hscr0%Hdr%nelect

     call test_charge(nfft,nelectron_exp,Hscr0%Hdr%nspden,rhor,Cryst%ucvol,&
       Hscr0%Hdr%usepaw,usexcnhat,usefinegrid,compch_sph,compch_fft,drude_plsmf)
     GN_drude_plsmf = drude_plsmf

     ! Read and in case make Epsilon^{-1} according the the options specified
     id_required=4; ikxc=0; approx_type=0; option_test=0; dim_kxcg=0
     ABI_MALLOC(kxcg,(nfft,dim_kxcg))

     call prompt(' Enter prefix for output files: ',prefix)
     fname_dump=TRIM(prefix)//'_SCR'

     orig_npwe = Er%npwe
     write(std_out,'(2a,I0)') ch10,' Number of plane waves is: ',Er%npwe
     write(std_out,'(a)',advance='no') ' Would you like to change it [Y/N] ?'
     read(std_in,*) ans
     if (ans=='Y'.or.ans=='y') then
       write(std_out,'(a)',advance='no') ' Enter new no. of plane waves (0 means use old value): '
       read(std_in,*) ii
       if (ii>0.or.ii<=Er%npwe) Er%npwe = ii
       if (ii<0.or.ii>Er%npwe) then
         MSG_ERROR(' Wrong value for no. of plane waves!')
       end if
     end if

     if (is_scr) Er%mqmem=1
     if (is_sus) Er%mqmem=0
     call mkdump_Er(Er,Vcp,Er%npwe,Gsphere%gvec,dim_kxcg,kxcg,id_required,approx_type,ikxc,option_test,&
       fname_dump,iomode,nfft,ngfft,comm)
     Er%mqmem=1

     call em1results_print(Er)

     write(std_out,'(2a)',advance='no') ch10,&
     ' Would you like to calculate the eigenvalues of eps^{-1}_GG''(omega) [Y/N] ? '
     read(std_in,*) ans

     if (ans=='Y'.or.ans=='y') then
       ABI_MALLOC(epsm1_eigen,(Er%npwe,Er%nomega))
       imax = 10
       if (Er%npwe < imax) imax = Er%npwe
       do iqibz=1,Er%nqibz
         call int2char4(iqibz,tagq)
         ABI_CHECK((tagq(1:1)/='#'),'Bug: string length too short!')
         fname_eigen=TRIM(prefix)//'_EM1_EIG_Q'//TRIM(tagq)
         if (open_file(fname_eigen,msg,newunit=unt_dump,status='replace',form='formatted') /= 0) then
           MSG_ERROR(msg)
         end if
         call decompose_epsm1(Er,iqibz,epsm1_eigen)
         write(unt_dump,'(a)')       '# First (max 10) eigenvalues of eps^{-1}(omega)'
         write(unt_dump,'(a,3f12.6)')'# q = ',Hscr0%qibz(:,iqibz)
         write(unt_dump,'(a)')       '# REAL omega [eV]  REAL(eigen(esp^-1(1,w)))  AIMAG(eigen(esp^-1(1,w))  ...'
         do iomega=1,Er%nomega_r
           write(unt_dump,'(21(es16.8))')REAL(Er%omega(iomega))*Ha_eV,&
             (REAL(epsm1_eigen(ii,iomega)),ii=1,imax),(AIMAG(epsm1_eigen(ii,iomega)),ii=1,imax)
         end do
         close(unt_dump)
       end do
       ABI_FREE(epsm1_eigen)
       ABI_FREE(kxcg)
     end if ! Calculate eigenvalues

     ! Analyze the PPmodel.
     write(std_out,'(2a)') ch10,' Would you like to analyse plasmon-pole models [Y/N] ? '
     read(std_in,*)ans

     if (ans=='Y'.or.ans=='y') then
       write(std_out,'(2a,f6.2,a)') ch10,' Plasma frequency for GN PPM is: ',GN_drude_plsmf*Ha_eV, ' eV'
       write(std_out,'(a)',advance='no') ' Would you like to change it [Y/N] ?'
       read(std_in,*) ans
       if (ans=='Y'.or.ans=='y') then
         write(std_out,'(2a)',advance='no') ch10,' Enter plasma frequency [eV]: '
         read(std_in,*) GN_drude_plsmf
         GN_drude_plsmf = GN_drude_plsmf/Ha_eV
       end if

       write(std_out,'(2a)') ch10,' Would you like to calculate the plasmon-pole model'
       write(std_out,'(a)',advance='no')       '        eigenvalues of eps^{-1}_GG''(omega) [Y/N] ? '
       read(std_in,*) ans

       if (ans=='Y'.or.ans=='y') then
         ABI_MALLOC(ppm_eigen,(PPm%npwc,Er%nomega))
         imax = 10; if (Er%npwe < imax) imax = Er%npwe
         do iqibz=1,Er%nqibz
           do ppmodel=1,2

             call int2char4(iqibz,tagq)
             ABI_CHECK((tagq(1:1)/='#'),'Bug: string length too short!')
             if (ppmodel==1) fname_dump=TRIM(prefix)//'_PPM_GN_EM1_EIG_Q'//TRIM(tagq)
             if (ppmodel==2) fname_dump=TRIM(prefix)//'_PPM_HL_EM1_EIG_Q'//TRIM(tagq)
             if (ppmodel==3) fname_dump=TRIM(prefix)//'_PPM_vdLH_EM1_EIG_Q'//TRIM(tagq)
             if (ppmodel==4) fname_dump=TRIM(prefix)//'_PPM_EF_EM1_EIG_Q'//TRIM(tagq)

             if (open_file(fname_eigen,msg,newunit=unt_dump,status='new',form='formatted') /= 0) then
               MSG_ERROR(msg)
             end if

             call ppm_free(PPm)
             if (ppmodel==1) then
               call ppm_init(PPm,Er%mqmem,Er%nqibz,Er%npwe,ppmodel,GN_drude_plsmf,Dtset%gw_invalid_freq)
             else
               call ppm_init(PPm,Er%mqmem,Er%nqibz,Er%npwe,ppmodel,drude_plsmf,Dtset%gw_invalid_freq)
             end if
             call setup_ppmodel(PPm,Cryst,Qmesh,Er%npwe,Er%nomega,Er%omega,Er%epsm1,nfft,Gsphere%gvec,ngfft,rhor(:,1),iqibz)

             call get_PPm_eigenvalues(PPm,iqibz,Er%Hscr%zcut,Er%nomega,Er%omega,Vcp,ppm_eigen)

             write(unt_dump,'(a)')       '# First (max 10) eigenvalues of eps^{-1}(omega) from Plasmon-pole model'
             write(unt_dump,'(a,3f12.6)')'# q = ',Hscr0%qibz(:,iqibz)
             select case(ppmodel)
             case(1)
               write(unt_dump,'(a)')     '# ppmodel = 1 : Godby - Needs'
             case(2)
               write(unt_dump,'(a)')     '# ppmodel = 2 : Hybertsen - Louie'
             case(3)
               write(unt_dump,'(a)')     '# ppmodel = 3 : von der Linden - Horsch'
             case(4)
               write(unt_dump,'(a)')     '# ppmodel = 4 : Engel - Farid'
             end select
             write(unt_dump,'(a)')       '# REAL omega [eV]  REAL(eigen(ppm_eps^-1(1,w)))  AIMAG(eigen(ppm_eps^-1(1,w))  ...'
             do iomega=1,Er%nomega_r
               write(unt_dump,'(21(es16.8))')REAL(Er%omega(iomega))*Ha_eV,&
                 (REAL(ppm_eigen(ii,iomega)),ii=1,imax),(AIMAG(ppm_eigen(ii,iomega)),ii=1,imax)
             end do
             close(unt_dump)

           end do !ppmodel
         end do ! iqibz
         ABI_FREE(ppm_eigen)
       end if ! Calculate PPM eigenvalues

       ! Optionally output eps^{-1}_GG''(w) for a given set of GG' and gridpoints
       write(std_out,'(2a)',advance='no') ch10,' Would you like to extract eps^{-1}_GG''(omega) for the PPM [Y/N] ?'
       read(std_in,*) ans

       if (ans=='Y'.or.ans=='y') then
         ! Reconstruct e^{-1}_GG'(w) according to PPmodel for statistical analysis.
         write(std_out,'(a)') ' Enter the number of frequency points in the'
         write(std_out,'(a)') '  interval 0 - freqremax (0 means same as input file ): '
         read(std_in,*) nfreqre
         if (nfreqre==0) then
           nfreqre   = Er%nomega_r
           nfreqim   = Er%nomega_i
           nfreq_tot = Er%nomega
           freqremax = REAL(Er%omega(Er%nomega_r))
           ABI_MALLOC(omega,(nfreq_tot))
           omega(:) = Er%omega(:)
           same_freqs = .TRUE.
         else
           write(std_out,'(a)') ' Enter the value of freqremax (in eV): '
           read(std_in,*) freqremax
           nfreqim   = Er%nomega_i
           nfreq_tot = nfreqre+nfreqim
           ABI_MALLOC(omega,(nfreqre+Er%nomega_i))
           do iomega=1,nfreqre
             omega(iomega) =  CMPLX((freqremax/REAL((nfreqre-1)))*(iomega-1),zero)
           end do
           omega(nfreqre+1:nfreq_tot) = Er%omega(Er%nomega_r+1:Er%nomega)
           same_freqs = .FALSE.
         end if ! frequencies

         do iqibz=1,Er%nqibz
           qtmp(:)=Er%qibz(:,iqibz); if (normv(qtmp,Cryst%gmet,'G')<GW_TOLQ0) qtmp(:)=zero
           call int2char4(iqibz,tagq)
           ABI_CHECK((tagq(1:1)/='#'),'Bug: string length too short!')

           ! At this time only the Godby-Needs and Hybertsen-Louie models
           ! TODO: Check the results from the others
           do ppmodel=1,2

             call ppm_free(PPm)
             if (ppmodel==1) then
               call ppm_init(PPm,Er%mqmem,Er%nqibz,Er%npwe,ppmodel,GN_drude_plsmf,Dtset%gw_invalid_freq)
             else
               call ppm_init(PPm,Er%mqmem,Er%nqibz,Er%npwe,ppmodel,drude_plsmf,Dtset%gw_invalid_freq)
             end if
             call setup_ppmodel(PPm,Cryst,Qmesh,Er%npwe,Er%nomega,Er%omega,Er%epsm1,&
               nfft,Gsphere%gvec,ngfft,rhor(:,1),iqibz)

             ! Prepare file for data on real omega axis
             if (ppmodel==1) fname_dump=TRIM(prefix)//'_PPM_w_GN_Q'//TRIM(tagq)
             if (ppmodel==2) fname_dump=TRIM(prefix)//'_PPM_w_HL_Q'//TRIM(tagq)
             if (ppmodel==3) fname_dump=TRIM(prefix)//'_PPM_w_vdLH_Q'//TRIM(tagq)
             if (ppmodel==4) fname_dump=TRIM(prefix)//'_PPM_w_EF_Q'//TRIM(tagq)

             if (open_file(fname_dump,msg,newunit=unt_dump,status='replace',form='formatted') /= 0) then
               MSG_ERROR(msg)
             end if

             ! Prepare file for data on imaginary omega axis
             if (ppmodel==1) fname_dump2=TRIM(prefix)//'_PPM_iw_GN_Q'//TRIM(tagq)
             if (ppmodel==2) fname_dump2=TRIM(prefix)//'_PPM_iw_HL_Q'//TRIM(tagq)
             if (ppmodel==3) fname_dump2=TRIM(prefix)//'_PPM_iw_vdLH_Q'//TRIM(tagq)
             if (ppmodel==4) fname_dump2=TRIM(prefix)//'_PPM_iw_EF_Q'//TRIM(tagq)

             if (open_file(fname_dump2,msg,newunit=unt_dump2,status='replace',form='formatted') /= 0) then
               MSG_ERROR(msg)
             end if

             ABI_MALLOC(em1_ppm,(nfreq_tot))

             ig1 = 0; ig2 = 0
             write(std_out,'(3a,I0,a,I0)') ch10,' Enter indices for G and G''.',&
               'Entering 0 exits the loop. iqibz = ',iqibz,' ppmodel = ',ppmodel

             do
               write(std_out,'(2(a),I0,a)',advance='NO') ch10,' Enter index for G (1 - ',Er%npwe,' ): '
               read(std_in,*)ig1
               if (ig1==0) EXIT
               if (ig1<0.OR.ig1>Er%npwe) MSG_ERROR(' index out of bounds')
               write(std_out,'(2(a),I0,a)',advance='NO') ch10,' Enter index for G'' (1 - ',Er%npwe,' ): '
               read(std_in,*)ig2
               if (ig2==0) EXIT
               if (ig2<0.OR.ig2>Er%npwe) MSG_ERROR(' index out of bounds')

               ! Generate the PPM representation of epsilon^-1
               call getem1_from_PPm_one_ggp(PPm,iqibz,Er%Hscr%zcut,nfreq_tot,omega,Vcp,em1_ppm,ig1,ig2)

               write(unt_dump,'(a,I1)') '# epsilon^-1_GG''(omega) from ppmodel = ',ppmodel
               write(unt_dump,'(2(a,i8),/,a,3f12.6,/,a,3i6,a,3i6,/,a,3F9.4,a,3F9.4,a,/a,f9.4,a,f9.4,a,/,a,/)')&
                 '# ig1= ',ig1,'    ig2= ',ig2,&
                 '# q = ',Er%qibz(:,iqibz),&
                 '# G = ',Er%gvec(:,ig1),'  G''= ',Er%gvec(:,ig2),&
                 '# G = (',MATMUL(two_pi*Cryst%gmet,Er%gvec(:,ig1)),&
                 ')  G''= (',MATMUL(two_pi*Cryst%gmet,Er%gvec(:,ig2)),')',&
                 '# 1/2|G|^2 =',half*normv(Er%gvec(:,ig1),Cryst%gmet,'G')**2,&
                 ' Ha 1/2|G''|^2 =',half*normv(Er%gvec(:,ig2),Cryst%gmet,'G')**2,' Ha',&
                 '#   omega [eV]           Re             Im '
               write(unt_dump2,'(a,I1)') '# epsilon^-1_GG''(iomega) from ppmodel = ',ppmodel
               write(unt_dump2,'(2(a,i8),/,a,3f12.6,/,a,3i6,a,3i6,/,a,3F9.4,a,3F9.4,a,/a,f9.4,a,f9.4,a,/,a,/)')&
                 '# ig1= ',ig1,'    ig2= ',ig2,&
                 '# q = ',Er%qibz(:,iqibz),&
                 '# G = ',Er%gvec(:,ig1),'  G''= ',Er%gvec(:,ig2),&
                 '# G = (',MATMUL(two_pi*Cryst%gmet,Er%gvec(:,ig1)),&
                 ')  G''= (',MATMUL(two_pi*Cryst%gmet,Er%gvec(:,ig2)),')',&
                 '# 1/2|G|^2 =',half*normv(Er%gvec(:,ig1),Cryst%gmet,'G')**2,&
                 ' Ha 1/2|G''|^2 =',half*normv(Er%gvec(:,ig2),Cryst%gmet,'G')**2,' Ha',&
                 '#   iomega [eV]           Re             Im '

               do iomega=1,nfreqre
                 if (same_freqs) then
                   write(unt_dump,'(f8.2,4x,4es16.8)') REAL(omega(iomega))*Ha_eV,em1_ppm(iomega),&
                     Er%epsm1(ig1,ig2,iomega,iqibz)
                 else
                   write(unt_dump,'(f8.2,4x,2es16.8)') REAL(omega(iomega))*Ha_eV,em1_ppm(iomega)
                 end if
               end do
               ! First output the iomega = 0 point
               write(unt_dump2,'(f8.2,4x,4es16.8)') AIMAG(omega(1))*Ha_eV,em1_ppm(1),&
                 Er%epsm1(ig1,ig2,1,iqibz)
               ! Then the rest
               do iomega=nfreqre+1,nfreq_tot
                 write(unt_dump2,'(f8.2,4x,4es16.8)') AIMAG(omega(iomega))*Ha_eV,em1_ppm(iomega),&
                   Er%epsm1(ig1,ig2,iomega,iqibz)
               end do
               write(unt_dump,*)
               write(unt_dump,*)
               write(unt_dump2,*)
               write(unt_dump2,*)
             end do ! Empty
             ABI_FREE(em1_ppm)
             close(unt_dump)
             close(unt_dump2)

           end do ! ppmodel
         end do ! iqibz
         ABI_FREE(omega)
       end if ! Output epsilon for PPM

       ! Optionally statistics for all PPMs
       write(std_out,'(2a)',advance='no') ch10,' Would you like to output statistics for all PPMs [Y/N] ?'
       read(std_in,*) ans

       if (ans=='Y'.or.ans=='y') then
         nfreqre   = Er%nomega_r
         nfreq_tot = Er%nomega
         freqremax = REAL(Er%omega(Er%nomega_r))
         ABI_MALLOC(real_omega,(nfreqre))
         real_omega(:) = REAL(Er%omega(1:nfreqre))

         do iqibz=1,Er%nqibz
           do ppmodel=1,2

             qtmp(:)=Er%qibz(:,iqibz)
             if (normv(qtmp,Cryst%gmet,'G')<GW_TOLQ0) qtmp(:)=zero

             call ppm_free(PPm)
             if (ppmodel==1) then
               call ppm_init(PPm,Er%mqmem,Er%nqibz,Er%npwe,ppmodel,GN_drude_plsmf,Dtset%gw_invalid_freq)
             else
               call ppm_init(PPm,Er%mqmem,Er%nqibz,Er%npwe,ppmodel,drude_plsmf,Dtset%gw_invalid_freq)
             end if
             call setup_ppmodel(PPm,Cryst,Qmesh,Er%npwe,Er%nomega,Er%omega,Er%epsm1,&
             nfft,Gsphere%gvec,ngfft,rhor(:,1),iqibz)

             ! Prepare ratios and density for the f-sum rule
             ABI_MALLOC_OR_DIE(qratio,(orig_npwe,orig_npwe), ierr)
             ABI_MALLOC_OR_DIE(rhoggp,(Er%npwe,Er%npwe), ierr)

             call cqratio(orig_npwe,Gsphere%gvec,qtmp,Cryst%gmet,Cryst%gprimd,qratio)
             ! Arrange n(G-G')->n(G,G')
             ierr=0
             do ig1=1,Er%npwe
               do ig2=1,Er%npwe
                 gmgp_idx = g2ifft(Gsphere%gvec(:,ig1)-Gsphere%gvec(:,ig2),ngfft)
                 if (gmgp_idx/=0) then
                   rhoggp(ig1,ig2)=CMPLX(rhog(1,gmgp_idx),rhog(2,gmgp_idx))
                 else
                   ierr=ierr+1
                   rhoggp(ig1,ig2)=czero
                 end if
               end do
             end do
             if (ierr/=0) then
               write(std_out,'(a,i0,a)')' Found ',ierr,' G1-G2 vectors falling outside the FFT box. '
             end if

             ! Prepare files
             call int2char4(iqibz,tagq)
             ABI_CHECK((tagq(1:1)/='#'),'Bug: string length too short!')
             if (ppmodel==1) fname_dump=TRIM(prefix)//'_norms_GN_Q'//TRIM(tagq)
             if (ppmodel==2) fname_dump=TRIM(prefix)//'_norms_HL_Q'//TRIM(tagq)
             if (open_file(fname_dump,msg, newunit=unt_dump, status='replace',form='formatted') /= 0) then
               MSG_ERROR(msg)
             end if
             write(unt_dump,'(a)') '# Various norms integrated through spline interpolation'
             write(unt_dump,'(a)') '# over all frequencies in the input file,'
             write(unt_dump,'(a)') '# for all G and G'' vectors.'
             write(unt_dump,'(a,I0)') '#               ppmodel: ',ppmodel
             write(unt_dump,'(a,I0)') '# Number of frequencies: ',nfreqre
             write(unt_dump,'(a,f12.6)') '# Maximum frequency    : ',freqremax
             write(unt_dump,'(a)') '# Columns:'
             write(unt_dump,'(2a)') '#  ig1      ig2   |eps-eps_PPM|/|eps|',&
               '   |eps-eps_PPM|    |eps|     |eps_PPM|            G                  G'''
             if (ppmodel==1) fname_dump2=TRIM(prefix)//'_f_sumrule_GN_Q'//TRIM(tagq)
             if (ppmodel==2) fname_dump2=TRIM(prefix)//'_f_sumrule_HL_Q'//TRIM(tagq)

             if (open_file(fname_dump2,msg,newunit=unt_dump2,status='replace',form='formatted') /= 0) then
               MSG_ERROR(msg)
             end if

             write(unt_dump2,'(a)') '# The fulfillment of the f-sum rule: I(epsilon) ='
             write(unt_dump2,'(a)') '#   int_0^{inf}{omega*Im[epsilon_G,G''(omega)]}/C_qGG'''
             write(unt_dump2,'(a)') '# C_qGG'' = '
             write(unt_dump2,'(a)') '#   -Pi/2*omega_p^2*(q+G)*(q+G'')/|q+G|^2*n(G-G'')/n(0)'
             write(unt_dump2,'(a)') '# for all G and G'' vectors.'
             write(unt_dump2,'(a,I0)') '#               ppmodel: ',ppmodel
             write(unt_dump2,'(a,I0)') '# Number of frequencies: ',nfreqre
             write(unt_dump2,'(a,f12.6)') '# Maximum frequency    : ',freqremax
             write(unt_dump2,'(a)') '# Columns:'
             write(unt_dump2,'(3a)') '#  ig1      ig2   I(epsilon)',&
               '   I(eps_PPM)   Re[n(G-G'')]    Im[n(G-G'')]    qratio      I1*C_qGG''',&
               ' Re[Omegatwsq] Im[Omegatwsq]   Re[omegatw]   Im[omegatw]    |G|    1/2|G|^2'

             ABI_MALLOC(em1_ppm,(nfreq_tot))
             ABI_MALLOC(ftab,(nfreqre))
             ABI_MALLOC(ysp,(3,nfreqre))
             ABI_MALLOC(work,(nfreqre))
             ABI_MALLOC(eint,(nfreqre))

             do ig1=1,Er%npwe
               write(std_out,'(2(a,I0))') ' ig1= ',ig1, ' of ',Er%npwe
               do ig2=1,Er%npwe
                 !ig2 = ig1
                 call getem1_from_PPm_one_ggp(PPm,iqibz,Er%Hscr%zcut,nfreq_tot,Er%omega,Vcp,em1_ppm,ig1,ig2)

                 ! Calculate norms in real
                 eps_diff=0; eps_norm=0; eps_ppm_norm=0
                 ftab(1:nfreqre) = ABS(Er%epsm1(ig1,ig2,1:nfreqre,iqibz)-em1_ppm(1:nfreqre))
                 call cspint(ftab,real_omega,nfreqre,real_omega(1),real_omega(nfreqre),ysp,eint,work,eps_diff)
                 ftab(1:nfreqre) = ABS(Er%epsm1(ig1,ig2,1:nfreqre,iqibz))
                 call cspint(ftab,real_omega,nfreqre,real_omega(1),real_omega(nfreqre),ysp,eint,work,eps_norm)
                 ftab(1:nfreqre) = ABS(em1_ppm(1:nfreqre))
                 call cspint(ftab,real_omega,nfreqre,real_omega(1),real_omega(nfreqre),ysp,eint,work,eps_ppm_norm)
                 write(unt_dump,'(2i6,f12.4,3es14.4,6i4)') ig1,ig2,eps_diff/eps_norm,eps_diff,&
                   eps_norm,eps_ppm_norm,Er%gvec(:,ig1),Er%gvec(:,ig2)

                 ! Evaluate the f-sum rule
                 if (ig1==ig2) then
                   ftab(1:nfreqre) = real_omega(1:nfreqre)*AIMAG(Er%epsm1(ig1,ig2,1:nfreqre,iqibz))
                 else
                   ! Dephase first - HERE epsm1 is changed!
                   call remove_phase(epsm1(ig1,ig2,:,1),Hscr_file(1)%nomega,phase)
                   ftab(1:nfreqre) = real_omega(1:nfreqre)*AIMAG(Er%epsm1(ig1,ig2,1:nfreqre,iqibz))
                 end if

                 call cspint(ftab,real_omega,nfreqre,real_omega(1),&
                   real_omega(nfreqre),ysp,eint,work,eps_diff)

                 if (ig1==ig2) then
                   factor = -two*pi*pi*REAL(rhoggp(ig1,ig2))*qratio(ig1,ig2)
                 else
                   rhoggp(ig1,ig2) = CMPLX(COS(phase),-SIN(phase))*rhoggp(ig1,ig2)
                   factor = -two*pi*pi*REAL(rhoggp(ig1,ig2))*qratio(ig1,ig2)
                 end if

                 if (ABS(qratio(ig1,ig2))>zero) then
                   value1 = eps_diff/factor
                   if (ppmodel==1) then
                     value2 = -pi*half*(REAL(PPm%bigomegatwsq(iqibz)%vals(ig1,ig2))&
                       /(REAL(PPm%omegatw(iqibz)%vals(ig1,ig2))))&
                       /factor*(2*sqrt(pi*rhoggp(1,1)))
                   else
                     value2 = -pi*half*(SQRT(REAL(PPm%bigomegatwsq(iqibz)%vals(ig1,ig2))))&
                      /factor*(2*sqrt(pi*rhoggp(1,1)))
                   end if
                 else
                   value1 = zero
                   value2 = zero
                 end if

                 write(unt_dump2,'(2i6,12es14.4)') ig1,ig2,value1,value2,&
                   REAL(rhoggp(ig1,ig2)),AIMAG(rhoggp(ig1,ig2)),qratio(ig1,ig2),&
                   eps_diff,REAL(PPm%bigomegatwsq(iqibz)%vals(ig1,ig2)),&
                   AIMAG(PPm%bigomegatwsq(iqibz)%vals(ig1,ig2)),&
                   REAL(PPm%omegatw(iqibz)%vals(ig1,ig2)),&
                   AIMAG(PPm%omegatw(iqibz)%vals(ig1,ig2)),&
                   normv(Er%gvec(:,ig1),Cryst%gmet,'G'),&
                   half*normv(Er%gvec(:,ig1),Cryst%gmet,'G')**2

               end do !ig2
             end do !ig1

             ABI_FREE(em1_ppm)
             ABI_FREE(ftab)
             ABI_FREE(ysp)
             ABI_FREE(work)
             ABI_FREE(eint)
             ABI_FREE(qratio)
             ABI_FREE(rhoggp)
             close(unt_dump); close(unt_dump2)

           end do ! ppmodel
         end do ! iqibz

         ABI_FREE(real_omega)
       end if ! Output statistics

       call ppm_free(PPm)
     end if ! If ppmodel>0

     ABI_FREE(rhor)
     ABI_FREE(rhog)
     ABI_FREE(nhat)

     call vcoul_free(Vcp)
     call em1results_free(Er)
     call gsph_free(Gsphere)

   case(4)
     ! Remove real frequencies ----------------------------------------------------------
     write(std_out,'(2(a))') ch10,' Do you want to remove every other real frequency  (= 1) ?'
     write(std_out,'(a)')         '  or specify for each real frequency individually  (= 2) ?'
     write(std_out,'(a)')         '  or remove ALL real frequencies                   (= 3) ?'
     read(std_in,*)choice

     ! Calculate the total number of real freq
     nfreqre = 0; nfreqim = 0
     do ifrq=1,Hscr_file(1)%nomega
       ! If frequency is not imaginary, count.
       if (AIMAG(Hscr_file(1)%omega(ifrq)) < tol8) nfreqre = nfreqre + 1
       if (REAL(Hscr_file(1)%omega(ifrq)) < tol8 .and. AIMAG(Hscr_file(1)%omega(ifrq))>tol8)  nfreqim = nfreqim + 1
     end do

     nfreq_tot = nfreqre + nfreqim ! Here nfreq_tot becomes the *true* number of freq
     write(std_out,'(2a,I0,a)') ch10,' Found ',nfreq_tot,' frequencies.'
     write(std_out,'(2(a,I0),2a)') ' ',nfreqre,' real, and ',nfreqim,' imaginary.',ch10

     ! Array with the index of frequencies to be kept.
     ABI_MALLOC(freq_indx,(nfreq_tot,1))
     freq_indx = 0

     select case(choice)
     case(1)
       ! Remove every other frequency
       write(std_out,'(2(a))') ch10,' Removing every other real frequency, i.e. every even one.'
       write(std_out,'(a)')         ' If the total number of frequencies is odd, the first and last one will be kept.'
       write(std_out,'(a)')         ' If the total number is even, the first one will still be in the final set.'

       ! Test for no real frequencies
       ABI_CHECK(nfreqre /= 0, "No real frequencies in file!")

       ii=nfreqre; nfreqre = 0
       do ifrq=1,ii
         if (.not. iseven(ifrq)) then
           nfreqre = nfreqre + 1; freq_indx(nfreqre,1) = ifrq
         end if
       end do
       write(std_out,'(2a,I0,a)') ch10,' ',nfreqre,' real frequencies will be kept.'

     case(2)
       ! Specify freq. individually
       ii = nfreqre; nfreqre = 0
       do ifrq=1,ii
         write(std_out,'(a,f12.6,a)') ' Would you like to keep freq. at: ',REAL(Hscr_file(1)%omega(ifrq))*Ha_eV,' eV? [y/n]'
         read(std_in,*) ans
         if (ans=='Y'.or.ans=='y') then
           nfreqre = nfreqre + 1; freq_indx(nfreqre,1) = ifrq
         end if
       end do
       write(std_out,'(2a,I0,a)') ch10,' ',nfreqre,' real frequencies will be kept.'

     case(3)
       ! Remove all real freq.
       nfreqre = 0

     case default
       MSG_ERROR("Invalid choice!")
     end select

     ! Add imaginary frequencies if any
     if (nfreqim > 0) then
       nfreqim = 0
       do ifrq=1,Hscr_file(1)%nomega
         if (AIMAG(Hscr_file(1)%omega(ifrq)) > tol8) then
           nfreqim = nfreqim + 1; freq_indx(nfreqre+nfreqim,1) = ifrq
         end if
       end do
     end if

     nfreq_tot = nfreqre + nfreqim ! Here nfreq_tot becomes the *true* number of freq
     write(std_out,'(3(a,i0),a)')' Finally, we have ',nfreq_tot,' frequencies. ',nfreqre,' real, and ',nfreqim,' imaginary.'

     call prompt(' Enter the full name of the final output file: ', fname_out)

     if (endswith(filenames(1), ".nc") .and. .not. endswith(fname_out, ".nc")) then
       fname_out = nctk_ncify(fname_out)
       call wrtout(std_out,"- Added .nc extension to output file as input data is in netcdf format.")
     end if

     call ioscr_wremove(filenames(1), hscr_file(1), fname_out, nfreq_tot, freq_indx, hscr_merge)
     ABI_FREE(freq_indx)

   case(5)
     ! Remove imaginary frequencies ------------------------------------------------

     ! Calculate the total number of freq
     nfreqre = 0; nfreqim = 0
     do ifrq=1,Hscr_file(1)%nomega
         ! If frequency is not imaginary, count.
       if (AIMAG(Hscr_file(1)%omega(ifrq))<tol8) nfreqre = nfreqre + 1
       if (REAL(Hscr_file(1)%omega(ifrq))<tol8.and.AIMAG(Hscr_file(1)%omega(ifrq))>tol8)  nfreqim = nfreqim + 1
     end do ! ifrq

     ! Test for no real frequencies
     if (nfreqim == 0) then
       MSG_ERROR("No imaginary frequencies in file!")
     end if

     nfreq_tot = nfreqre + nfreqim ! Here nfreq_tot becomes the *true* number of freq
     write(std_out,'(2a,I0,a)') ch10,' Found ',nfreq_tot,' frequencies.'
     write(std_out,'(2(a,I0),2a)') ' ',nfreqre,' real, and ',nfreqim,' imaginary.',ch10

     ABI_MALLOC(freq_indx,(nfreq_tot,1))
     freq_indx = 0

     ! Specify freq. individually
     ii=nfreq_tot; nfreqim = 0
     do ifrq=nfreqre+1,ii
       write(std_out,'(a,f12.6,a)')&
        ' Would you like to keep imaginary freq. at: ',AIMAG(Hscr_file(1)%omega(ifrq))*Ha_eV,' eV? [y/n]'
       read(std_in,*) ans
       if (ans=='Y'.or.ans=='y') then
         nfreqim = nfreqim + 1; freq_indx(nfreqre+nfreqim,1) = ifrq
       end if
     end do ! ifrq
     write(std_out,'(2a,I0,a)') ch10,' ',nfreqim,' imaginary frequencies will be kept.'

     ! Add real frequencies if any
     if (nfreqre > 0) then
       nfreqre = 0
       do ifrq=1,Hscr_file(1)%nomega
         if (AIMAG(Hscr_file(1)%omega(ifrq)) < tol8) then
           nfreqre = nfreqre + 1; freq_indx(nfreqre,1) = ifrq
         end if
       end do
     end if

     nfreq_tot = nfreqre + nfreqim ! Here nfreq_tot becomes the *true* number of freq
     write(std_out,'(2a,I0,a)') ch10,' Finally, we have ',nfreq_tot,' frequencies.'
     write(std_out,'(2(a,I0),2a)') ' ',nfreqre,' real, and ',nfreqim,' imaginary.',ch10

     call prompt(' Enter the full name of the final output file: ',fname_out)

     if (endswith(filenames(1), ".nc") .and. .not. endswith(fname_out, ".nc")) then
       fname_out = nctk_ncify(fname_out)
       call wrtout(std_out,"- Added .nc extension to output file as input data is in netcdf format.")
     end if

     call ioscr_wremove(filenames(1), hscr_file(1), fname_out, nfreq_tot, freq_indx, hscr_merge)

     ABI_FREE(freq_indx)

   case(6)
     ! Model screening -------------------------------------------------------------
     MSG_ERROR("Model screening has been removed")

   !case(9)
   !  TODO
   !  ! netcdf --> Fortran converter -------------------------------------------------------------
   !  call prompt(' Enter the name of the final output Fortran file: ',fname_out)
   !  ! fname_out extension should be consistent with filenames(1)
   !  call ioscr_nc2fort(filenames(1), fname_out)

   case default
     ! Bail if choice is wrong
     write(std_out,*) ' Invalid choice! Exiting...'
     goto 100
   end select

 end if ! Single file mode

 call timein(tcpu,twall)

 tsec(1)=tcpu-tcpui
 tsec(2)=twall-twalli

 write(std_out, '(a,a,a,f13.1,a,f13.1)' )  '-',ch10,'- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)

 ! =====================
 ! ==== Free memory ====
 ! =====================
 ABI_FREE(filenames)
 ABI_SFREE(kxcg)
 ABI_SFREE(foundq)

 call cryst%free()
 call kmesh_free(Kmesh)
 call kmesh_free(Qmesh)
 call destroy_mpi_enreg(MPI_enreg)

 nullify(Hscr0)
 call hscr_free(Hscr_merge)

 do ifile=1,nfiles
   call hscr_free(Hscr_file(ifile))
 end do
 ABI_FREE(Hscr_file)

 call flush_unit(std_out)

 call abinit_doctor("__mrgscr")

 100 call xmpi_end()

 end program mrgscr
!!***
