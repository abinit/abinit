!{\src2tex{textfont=tt}}
!!****f* ABINIT/mlwfovlp
!! NAME
!! mlwfovlp
!!
!! FUNCTION
!! Routine which computes overlap M_{mn}(k,b) and projection A_{mn}(k)
!! for Wannier code (www.wannier.org f90 version).
!! Various file are written (wannier90.*) which can be used to run a
!! separate wannier calculation with the wannier90 code.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2017 ABINIT group (BAmadon,CEspejo,FJollet,TRangel)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cg(2,mcg)=planewave coefficients of wavefunctions.
!!  cprj(natom,mcprj)= <p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  dtfil <type(datafiles_type)>=variables related to files
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mgfft=maximum size of 1D FFTs
!!  mgfftc=maximum size of 1D FFTs (coarse grid)
!!  mkmem =number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw.
!   natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nfft=(effective) number of FFT grid points (for this processor) (see NOTES at beginning of scfcv)
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  occ(mband*nkpt*nsppol) Occupation number for each band (often 2) for each k point.
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  ucvol=unit cell volume (bohr**3)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine mlwfovlp(atindx1,cg,cprj,dtset,dtfil,eigen,gprimd,hdr,kg,&
& mband,mcg,mcprj,mgfftc,mkmem,mpi_enreg,mpw,natom,&
& nattyp,nfft,ngfft,nkpt,npwarr,nsppol,ntypat,occ,&
& pawang,pawrad,pawtab,prtvol,psps,rprimd,ucvol,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wannier90
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_io_tools, only : delete_file, get_unit
 use m_pawang,   only : pawang_type
 use m_pawrad,   only : pawrad_type
 use m_pawtab,   only : pawtab_type
 use m_pawcprj,  only : pawcprj_type
#ifdef FC_NAG
 use f90_unix_dir
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mlwfovlp'
 use interfaces_14_hidewrite
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts
 use interfaces_65_paw
 use interfaces_67_common, except_this_one => mlwfovlp
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mcprj,mgfftc,mkmem,mpw,natom,nfft,nkpt
 integer,intent(in) :: nsppol,ntypat,prtvol
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(hdr_type),intent(in) :: hdr
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx1(natom)
 integer :: kg(3,mpw*mkmem),nattyp(ntypat),ngfft(18),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),gprimd(3,3),rprimd(3,3)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: xred(3,natom)
 type(pawcprj_type) :: cprj(natom,mcprj)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: band_index,cplex,i,iatom,iband,iband1,iband2,icgtemp
 integer :: ig,ii,ikg,ierr
 integer :: ikpt,ikpt1,ikpt2,ilmn,intot,isppol,itypat
 integer :: iun(nsppol),iun_plot,iwan,jband,jband1,jband2,jj,jj1,jj2,jj3
 integer :: kk,lmn_size,lproj,lwanniersetup,mwan,mgfft,n1
 integer :: n1tmp,n2,n2tmp,n3,n3tmp,n4,n5,n6,nband_k
 integer :: nntot,npw_k,num_nnmax,spacing
 integer :: tim_fourwf
 integer :: master,max_num_bands,nprocs,spaceComm,spin,rank
!integer :: j,k,l
 integer  :: nwan(nsppol),nband_inc(nsppol),num_bands(nsppol)
 real(dp) :: corrvdw
 real(dp) :: weight
 complex(dpc) :: caux,caux2,caux3
 logical :: gamma_only,leig,lmmn,lwannierrun,spinors !,have_disentangled
 character(len=20) :: wfnname
 character(len=500) :: message
 character(len=fnlen) :: seed_name(nsppol)
 character(len=fnlen) :: fname,filew90_win(nsppol),filew90_wout(nsppol),filew90_amn(nsppol),filew90_ramn(nsppol)
 character(len=fnlen) :: filew90_mmn(nsppol),filew90_eig(nsppol)
!arrays
 integer :: g1temp(3),ngkpt(3)
 integer,allocatable :: g1(:,:,:),gbound(:,:),icg(:,:)
 integer,allocatable:: iwav(:,:,:),kg_k(:,:),ovikp(:,:)
 integer,allocatable :: proj_l(:,:),proj_m(:,:),proj_radial(:,:)
 real(dp) :: real_lattice(3,3)
 real(dp) :: recip_lattice(3,3),spreadw(3)
 real(dp),allocatable :: cm1(:,:,:,:,:,:),cm2_paw(:,:,:),csix(:,:,:,:),cwavef(:,:)
 real(dp),allocatable :: denpot(:,:,:)
 real(dp),allocatable :: eigenvalues_w(:,:,:),fofgout(:,:),fofr(:,:,:,:)
 real(dpc),allocatable :: occ_arr(:,:,:),occ_wan(:,:,:)
 real(dp),allocatable :: proj_site(:,:,:),proj_x(:,:,:),proj_z(:,:,:),proj_zona(:,:)
 real(dp),allocatable :: tdocc_wan(:,:)
 real(dp),allocatable :: wann_centres(:,:,:),wann_spreads(:,:),xcart(:,:)
 complex(dpc),allocatable :: A_paw(:,:,:,:)
 complex(dpc),allocatable :: M_matrix(:,:,:,:,:),U_matrix(:,:,:,:)
 complex(dpc),allocatable :: U_matrix_opt(:,:,:,:)
 complex(dpc),pointer :: A_matrix(:,:,:,:)
 logical,allocatable :: band_in(:,:),lwindow(:,:,:)
 character(len=3),allocatable :: atom_symbols(:)
 logical,allocatable::just_augmentation(:,:)

!************************************************************************

 ABI_UNUSED(hdr%natom)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!1) Initialize variables and allocations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Some initialization and checks
!
 lwanniersetup=1 ! 1 is mandatory ( 0 is for debug)
!to use lwanniersetup=0, one would need
!to define which bands to exclude.
 lwannierrun=.true.   ! .false. and .true. are possible
 lmmn=.true.          ! .false. and .true. are possible
 leig=.true.          ! .false. and .true. are possible
!
 gamma_only=.false. !not yet implemented
 spinors=.false. !not yet implemented
 spin=dtset%useria !spin polarization, just used in case of nsppol==2
!0=> both (default)
!1=> spin up
!2=> spin down
!
!mpi initialization
!
 spaceComm=MPI_enreg%comm_cell
 nprocs=xmpi_comm_size(spaceComm)
 rank=MPI_enreg%me_kpt
 master=0

!write(std_out,'("master ",i3," rank ",i3," nprocs ",i3)') master,rank,nprocs
!
!Generate seed names for wannier90 files, and file names
!
 call mlwfovlp_seedname(dtfil%fnameabo_w90,filew90_win,filew90_wout,filew90_amn,&
& filew90_ramn,filew90_mmn,filew90_eig,nsppol,seed_name)
!
!Check the validity of input variables
!FIXME: this is not a check, and prints a warning even if the input is fine!
!must be changed to not print anything if kptopt 3 and istwfk 1 (the latter is easier to check)
!
 if(rank==master) then
   write(message, '(a,a,a,a)' ) ch10,&
&   '   mlwfovlp:  you should give k-point in the full brillouin zone ',ch10,&
&   '   with explicit k-points (or kptopt=3) and istwfk 1'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if
!
 if(MPI_enreg%paral_spinor==1) then
   message = ' Parallelization over spinorial components not yet available !'
   MSG_ERROR(message)
 end if

 if(nsppol==2) then
   if(spin==1) then
     write(message, '(3a)' ) ch10,&
&     '   mlwfovlp:  Calculating matrices for spin-up polarization  ',ch10
   elseif(spin==2) then
     write(message, '(3a)' ) ch10,&
&     '   mlwfovlp:  Calculating matrices for spin-down polarization  ',ch10
   end if
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   if(spin <0 .or. spin>2) then
     message = '  mlwfovlp:  spin variable should be equal to 0, 1 or 2 '
     MSG_ERROR(message)
   end if
 end if
!
!get lattice parameters in wannier90 format
!
 do i=1, 3
   real_lattice(:,i)=Bohr_Ang*rprimd(i,:)
   recip_lattice(:,i)=two_pi*gprimd(i,:)/Bohr_Ang
 end do
!
 if(psps%npsp/=psps%ntypat) then
   MSG_ERROR("prb npsp")
 end if
!
!Allocations.
!
 num_nnmax=12 !limit fixed for compact structure in wannier_setup.
 ABI_ALLOCATE(g1,(3,nkpt,num_nnmax))
 ABI_ALLOCATE(ovikp,(nkpt,num_nnmax))
 ABI_ALLOCATE(atom_symbols,(natom))
 ABI_ALLOCATE(xcart,(3,natom))
 ABI_ALLOCATE(band_in,(mband,nsppol))
 ABI_ALLOCATE(proj_site,(3,mband,nsppol))
 ABI_ALLOCATE(proj_l,(mband,nsppol))
 ABI_ALLOCATE(proj_m,(mband,nsppol))
 ABI_ALLOCATE(proj_radial,(mband,nsppol))
 ABI_ALLOCATE(proj_x,(3,mband,nsppol))
 ABI_ALLOCATE(proj_z,(3,mband,nsppol))
 ABI_ALLOCATE(proj_zona,(mband,nsppol))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!2) Call to  Wannier setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 nullify(A_matrix)

 call mlwfovlp_setup(atom_symbols,band_in,dtset,filew90_win,gamma_only,&
& g1,lwanniersetup,mband,natom,nband_inc,nkpt,&
& nntot,num_bands,num_nnmax,nsppol,nwan,ovikp,&
& proj_l,proj_m,proj_radial,proj_site,proj_x,proj_z,proj_zona,&
& real_lattice,recip_lattice,rprimd,seed_name,spin,spinors,xcart,xred)



!
 do isppol=1, nsppol
   if(spin.ne.0 .and. spin.ne.isppol) cycle
   write(message, '(6a)' ) ch10,&
&   '   mlwfovlp :  mlwfovlp_setup done -',ch10,&
&   '-  see ',trim(filew90_wout(isppol)),' for details.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end do
!
!some allocations after wannier90 setup
!
 max_num_bands=maxval(num_bands(:))
 mwan=maxval(nwan(:))
 ABI_ALLOCATE(eigenvalues_w,(max_num_bands,nkpt,nsppol))
 ABI_ALLOCATE(M_matrix,(max_num_bands,max_num_bands,nntot,nkpt,nsppol))
 ABI_ALLOCATE(A_matrix,(max_num_bands,mwan,nkpt,nsppol))
 ABI_ALLOCATE(iwav,(mband,nkpt,nsppol))

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!3) Write Eigenvalues (file seed_name.eig)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
 if(leig) then
!  
!  Assign file unit numbers
   if(rank==master) then
     iun(1)=444
     if(nsppol==2) iun(2)=455
     do isppol=1,nsppol
       if(spin.ne.0 .and. spin.ne.isppol) cycle
       open(unit=iun(isppol),file=filew90_eig(isppol),form='formatted',status='unknown')
     end do
   end if !rank==master
!  Loop to write eigenvalues
   band_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
       jband=0
       do iband=1,mband
         if(band_in(iband,isppol)) then
           jband=jband+1
           if(spin.eq.0 .or. spin.eq.isppol)then
!            Writing data
             if(rank==master) write(iun(isppol), '(2i6,4x,f10.5)' ) jband,ikpt,Ha_eV*eigen(iband+band_index)
!            Finish writing, now save eigenvalues
             eigenvalues_w(jband,ikpt,isppol)=Ha_eV*eigen(iband+band_index)
           end if !spin
         end if
       end do !iband
       band_index=band_index+nband_k
     end do !ikpt
   end do  !nsppol
   if(rank==master) then
     do isppol=1,nsppol
       if(spin.ne.0 .and. spin.ne.isppol) cycle
       close(iun(isppol))
     end do
     write(message, '(a,a)' ) ch10,&
&     '   mlwfovlp :  eigenvalues written'
     call wrtout(std_out,  message,'COLL')
   end if !master
 end if !leig
!else if( leig . and. lwannierun) then
!read .eig file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!4) Calculate overlaps (file seed_name.mmn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!First calculate indices and shift
!
!write(std_out,*) "Computes shift for cg"
 write(message, '(a,a)' ) ch10,&
& '   mlwfovlp : compute shifts for g-points '
 call wrtout(std_out,  message,'COLL')
!----------------------------------------------------------------------
!Compute shifts for g points (icg,iwav)
!(here mband is not used, because shifts are internal
!variables of abinit)
!----------------------------------------------------------------------
!write(std_out,*) mpw*dtset%nspinor*mband*mkmem*nsppol
 ABI_ALLOCATE(icg,(nsppol,nkpt))
 icg=0
 icgtemp=0
 iwav(:,:,:)=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
!    
!    MPI:cycle over k-points not treated by this node
!    

     if (nprocs>1 ) then !sometimes we can have just one processor
       if ( ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-rank)  /=0) CYCLE
     end if
     
!    write(std_out,*)'rank',rank,'ikpt',ikpt,'isppol',isppol
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
!    write(std_out,*) ikpt+(isppol-1)*nkpt,nkpt
     npw_k=npwarr(ikpt)
     do iband=1,nband_k
       if(iband.gt.mband) then
         write(message,'(a,3i0)')" mband",iband,mband,nband_k
         MSG_ERROR(message)
       end if
       iwav(iband,ikpt,isppol)= &
&       (iband-1)*npw_k*dtset%nspinor+icgtemp
     end do ! iband
     icgtemp=icgtemp+ npw_k*dtset%nspinor*nband_k
!    icg(isppol,ikpt)=icgtemp
!    write(std_out,*) "icg", isppol,ikpt,icg(isppol,ikpt)
   end do  ! ikpt
 end do   ! isppol
!write(std_out,*) "shift for cg computed"
 ABI_DEALLOCATE(icg)
!
!Shifts computed.
!
 if( lmmn) then
!  
!  In case of parallelization write out cg for all k-points
!  
   if (nprocs > 1) then
!    
     if(prtvol>0) then
       write(message, '(3a)' ) ch10,&
&       '   mlwfovlp :  Creating temporary files with cg and cprj (PAW)',ch10
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end if
!    
     do isppol=1,nsppol
       if(spin.ne.0 .and. spin.ne.isppol) cycle
       do ikpt=1,nkpt
!        
!        MPI:cycle over k-points not treated by this node
!        
         if (nprocs>1 ) then !sometimes we can have just one processor
           if ( ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-rank)  /=0) CYCLE
         end if
         
!        write(std_out,*)'writing kpt ',ikpt,'isppol',isppol,' by node ', rank
         write(wfnname,'(a,I5.5,".",I1)') trim(dtfil%fnametmp_cg),ikpt,isppol
         iun_plot=1000+ikpt+ikpt*(isppol-1)

         open (unit=iun_plot, file=wfnname,form='unformatted')
         npw_k=npwarr(ikpt)
         do iband=1,mband
           do ig=1,npw_k*dtset%nspinor
             write(iun_plot) (cg(i,ig+iwav(iband,ikpt,isppol)),i=1,2)
           end do
         end do
         close(iun_plot)
       end do !ikpt
     end do !isppol
!    
!    In the PAW case we also need to write out cprj into files
!    
     if(psps%usepaw==1) then
!      
!      big loop on atoms, kpts, bands and lmn
!      
       ikpt2=0
       do isppol=1,nsppol
         if(spin.ne.0 .and. spin.ne.isppol) cycle
         do ikpt=1,nkpt
!          
!          MPI:cycle over k-points not treated by this node
!          
           if (nprocs>1 ) then !sometimes we can have just one processor
             if ( ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-MPI_enreg%me)  /=0) CYCLE
           end if
           
           ikpt2=ikpt2+1 !sums just on the k-points treated by this node
!          
           write(wfnname,'(a,I5.5,".",I1)') trim(dtfil%fnametmp_cprj),ikpt,isppol
           iun_plot=1000+ikpt
           open (unit=iun_plot, file=wfnname,form='unformatted')
!          
           do iband=1,mband*dtset%nspinor
             ig=iband+(ikpt2-1)*mband*dtset%nspinor +(isppol-1)*nkpt*mband*dtset%nspinor !index for cprj(:,ig)
!            
             do iatom=1,natom
               itypat=dtset%typat(iatom)
               lmn_size=pawtab(itypat)%lmn_size
!              
               do ilmn=1,lmn_size
                 write(iun_plot) (( cprj(iatom,ig)%cp(i,ilmn)),i=1,2)
               end do !ilmn
             end do !iatom
           end do !iband

           close(iun_plot)
         end do !ikpt
       end do !isppol
     end if !usepaw==1

!    
!    
   end if !MPI nprocs>1
!  
!  End of MPI preliminarities
!  Calculate PW contribution of overlaps
!  
   ABI_ALLOCATE(cm1,(2,mband,mband,nntot,nkpt,nsppol))
   call mlwfovlp_pw(cg,cm1,g1,iwav,kg,mband,&
&   mkmem,mpi_enreg,mpw,nfft,ngfft,nkpt,nntot,&
&   npwarr,dtset%nspinor,nsppol,ovikp,dtfil%fnametmp_cg,spin)
   write(message, '(a,a)' ) ch10,&
&   '   mlwfovlp : PW part of overlap computed   '
   call wrtout(std_out,  message,'COLL')
!  
!  compute PAW Contribution and add it to PW contribution
!  
   if(psps%usepaw==1) then
     write(message, '(a,a)' ) ch10,&
&     '** smatrix_pawinit : PAW part of overlap  '
     call wrtout(std_out,  message,'COLL')
     ABI_ALLOCATE(cm2_paw,(2,mband,mband))
     do isppol=1,nsppol
       if(spin.ne.0 .and. spin.ne.isppol) cycle
       do ikpt1=1,nkpt
!        
!        MPI:cycle over k-points not treated by this node
!        
         if (nprocs>1 ) then !sometimes we can have just one processor
           if ( ABS(MPI_enreg%proc_distrb(ikpt1,1,isppol)-rank)  /=0) CYCLE
         end if
         
         write(message, '(a,i6,a,2i6)' ) &
&         '   processor',rank,' computes PAW part for kpt and spin',ikpt1,isppol
         call wrtout(std_out,  message,'COLL')

         do intot=1,nntot
           ikpt2= ovikp(ikpt1,intot)
           g1temp(:)=g1(:,ikpt1,intot)
           call smatrix_pawinit(atindx1,cm2_paw,cprj,ikpt1,ikpt2,isppol,&
&           g1temp,gprimd,dtset%kpt,mband,mband,mkmem,mpi_enreg,&
&           natom,dtset%nband,nkpt,dtset%nspinor,nsppol,dtset%ntypat,pawang,pawrad,pawtab,rprimd,&
&           dtfil%fnametmp_cprj,dtset%typat,xred)
!          cm1(:,:,:,intot,ikpt1,isppol)=four_pi*cm2_paw(:,:,:)
!           write(6,*) "ikpt1=",ikpt1
!           do iband=1,mband
!             write(6,*) "iband=",iband
!             write(6,*) "Wannier PW       overlap",cm1(:,iband,iband,intot,ikpt1,isppol)
!             write(6,*) "Wannier PAW      overlap",four_pi*cm2_paw(:,iband,iband)
!             write(6,*) "Wannier PW+PAW   overlap",cm1(:,iband,iband,intot,ikpt1,isppol)+four_pi*cm2_paw(:,iband,iband)
!           enddo
           cm1(:,:,:,intot,ikpt1,isppol)=cm1(:,:,:,intot,ikpt1,isppol)+four_pi*cm2_paw(:,:,:)
         end do ! intot
       end do ! ikpt1
     end do ! isppol
     ABI_DEALLOCATE(cm2_paw)
     write(message, '(a,a)' ) ch10,&
&     '   mlwfovlp : PAW part of overlap computed '
     call wrtout(std_out,  message,'COLL')
   end if ! usepaw
!  
   call xmpi_barrier(spaceComm)
   call xmpi_sum(cm1,spaceComm,ierr)
!  
!  write overlap for separate calculation of wannier functions
!  
   if(rank==master) then
     do isppol=1,nsppol !we write separate output files for each isppol
       if(spin.ne.0 .and. spin.ne.isppol) cycle
       iun(isppol)=220+isppol
       open(unit=iun(isppol),file=filew90_mmn(isppol),form='formatted',status='unknown')
       write(iun(isppol),*) "nnkp version 90"
       write(iun(isppol),*) num_bands(isppol),nkpt,nntot
     end do
   end if ! rank==master

   do isppol=1,nsppol
     if(spin.ne.0 .and. spin.ne.isppol) cycle
     do ikpt1=1,nkpt
       do intot=1,nntot
         if( rank==master) write(iun(isppol),'(2i6,3x,3x,3i5)') ikpt1,ovikp(ikpt1,intot),(g1(jj,ikpt1,intot),jj=1,3)
         jband2=0
         do iband2=1,mband ! the first index is faster
           if(band_in(iband2,isppol)) then
             jband2=jband2+1
             jband1=0
             do iband1=1,mband
               if(band_in(iband1,isppol)) then
                 jband1=jband1+1
                 if(rank==master) write(iun(isppol),*) &
&                 cm1(1,iband1,iband2,intot,ikpt1,isppol),cm1(2,iband1,iband2,intot,ikpt1,isppol)
                 M_matrix(jband1,jband2,intot,ikpt1,isppol)=&
&                 cmplx(cm1(1,iband1,iband2,intot,ikpt1,isppol),cm1(2,iband1,iband2,intot,ikpt1,isppol))
!                write(2211,*) ikpt1,intot,iband1,iband2
!                write(2211,*) cm1(1,iband1,iband2,intot,ikpt1,isppol),cm1(2,iband1,iband2,intot,ikpt1,isppol)
               end if ! band_in(iband1)
             end do ! iband1
           end if ! band_in(iband2)
         end do ! iband2
       end do !intot
     end do !ikpt
     if( rank==master ) then
       close(iun(isppol))
       write(message, '(3a)' )  '   ',trim(filew90_mmn(isppol)),' written'
       call wrtout(std_out,  message,'COLL')
     end if !rank==master
   end do !isppol
!  
   ABI_DEALLOCATE(cm1)
!  
!  Write down part of the matrix to the output file
!  This is for the automatic tests
!  
   if(rank==master) then
     write(message, '(4a)' ) ch10,&
&     '   Writing top of the overlap matrix: M_mn(ikb,ik)',ch10,&
&     '   m=n=1:3, ikb=1, ik=1'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
!    
!    just write down the first 3 elements
!    
     do isppol=1,nsppol
       if(spin.ne.0 .and. spin.ne.isppol) cycle
       write(message, '( " " )')
       if (nsppol>1 ) then
         if (isppol==1) write(message,'(2a)')trim(message),'   spin up:'
         if (isppol==2) write(message,'(2a)')trim(message),'   spin down:'
       end if
       do ii=1,3
         if(ii>num_bands(isppol)) cycle
         write(message,'(3a)') trim(message),ch10,';   ( '
         do jj=1,3
           if(jj>num_bands(isppol))cycle
           write(message, '(a,2f11.6,a)') trim(message),&
&           M_matrix(ii,jj,1,1,isppol),' , '
         end do
         write(message,'(2a)') trim(message),'    ) '
       end do
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end do
!    
!    Now write down bottom of the matrix
!    
     write(message, '(4a)' ) ch10,&
&     '   Writing bottom of the overlap matrix: M_mn(ikb,ik)',ch10,&
&     '   m=n=num_bands-2:num_bands, ikb=nntot, ik=nkpt'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
!    
     do isppol=1,nsppol
       if(spin.ne.0 .and. spin.ne.isppol) cycle
       write(message, '( " " )')
       if (nsppol>1 ) then
         if (isppol==1) write(message,'(2a)')trim(message),'   spin up:'
         if (isppol==2) write(message,'(2a)')trim(message),'   spin down:'
       end if
       do ii=num_bands(isppol)-2,num_bands(isppol)
         if(ii<1) cycle
         write(message,'(3a)') trim(message),ch10,';   ( '
         do jj=num_bands(isppol)-2,num_bands(isppol)
           if(jj<1)cycle
           write(message, '(a,2f11.6,a)') trim(message),&
&           M_matrix(ii,jj,nntot,nkpt,isppol),' , '
         end do !j
         write(message,'(2a)') trim(message),'    ) '
       end do !ii
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end do !isppol
   end if !rank==master
!  
!  erase temporary files created for parallel runs
!  
   if (nprocs > 1) then
!    
     if(prtvol>0) then
       write(message, '(3a)' ) ch10,&
&       '   mlwfovlp :  Removing temporary files with cg and cprj (PAW)',ch10
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end if
!    
!    Just master  node will remove the files
!    
     if(rank==master) then
       do isppol=1,nsppol
         if(spin.ne.0 .and. spin.ne.isppol) cycle
         do ikpt=1,nkpt
           write(wfnname,'(a,I5.5,".",I1)') trim(dtfil%fnametmp_cg),ikpt,isppol
           call delete_file(wfnname,ierr)
           if(psps%usepaw==1) then
             write(wfnname,'(a,I5.5,".",I1)') trim(dtfil%fnametmp_cprj),ikpt,isppol
             call delete_file(wfnname,ierr)
           end if
         end do !ikpt
       end do !isppol
     end if
   end if !MPI nprocs>1
!  
 end if !lmmn
!if ( lmmn== .false. .and. lwannierun ) the
!read .mmn file
!
!Deallocate arrays  no longer used
!
 ABI_DEALLOCATE(ovikp)
 ABI_DEALLOCATE(g1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!5) Calculate initial projections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if(dtset%w90iniprj/=0 )  then


!  
!  Set value for lproj (type of projections to be computed)
!  In PAW, options 5 and 6 are not in use.
!  5 means that there will be a contribution from inside the spheres and another from the PW part
!  6 means that we take into account just the inside-spheres contribution
!  2 means that PW part will be calculated

!  
   lproj=dtset%w90iniprj
   if(dtset%w90iniprj == 5 ) lproj=2 ! Necessary to calculate PW contribution
!  
   ABI_ALLOCATE(just_augmentation,(mwan,nsppol))
   just_augmentation(:,:)=.false.

   if( psps%usepaw==1 .and. (dtset%w90iniprj==2 .or. dtset%w90iniprj>4)) then
     if (dtset%w90iniprj==6) just_augmentation(:,:)=.true.
     if (dtset%w90iniprj==5) then
       do isppol=1,nsppol
         if(spin.ne.0 .and. spin.ne.isppol) cycle
         do iwan=1,nwan(isppol)
!          
!          Trick to skip the planewave contribution for some Wannier functions
!          (Not in production).
!          
           if(proj_radial(iwan,isppol) > 4) then
             just_augmentation(iwan,isppol)=.true.
             proj_radial(iwan,isppol)=proj_radial(iwan,isppol)-3
             write(message, '(2a,2i4)' ) &
&             '   ','Skiping planewave contribution for iwan, ispin=',iwan,isppol
             call wrtout(std_out,  message,'COLL')
           end if !proj_radial>4
         end do !iwan
       end do !isppol
     end if !w90iniprj == 5
   end if !paw
!  
!  Call mlwfovlp_proj (plane waves part of projections)
!  
   if (dtset%w90iniprj/=6) then ! option 6 not yet in use
     call mlwfovlp_proj(A_matrix,band_in,cg,cprj,dtset,gprimd,just_augmentation,kg,&
&     lproj,max_num_bands,mband,mkmem,mpi_enreg,mpw,mwan,natom,&
&     nattyp,nkpt,npwarr,&
&     dtset%nspinor,nsppol,ntypat,num_bands,nwan,pawtab,proj_l,proj_m,&
&     proj_radial,proj_site,proj_x,proj_z,proj_zona,psps,spin,ucvol)
     write(message, '(a,a,a,a)' ) ch10,&
&     '   mlwfovlp:  mlwfovlp_proj done -',ch10,&
&     '   Projectors computed.'
     call wrtout(std_out,  message,'COLL')
   end if !w90proj/=6
!  
!  Calculate inside-sphere part of projections (PAW)
!  
   if (psps%usepaw ==1 .and. ( dtset%w90iniprj>4)) then
     ABI_ALLOCATE(A_paw,(max_num_bands,mwan,nkpt,nsppol))
     call mlwfovlp_projpaw(A_paw,band_in,cprj,just_augmentation,max_num_bands,mband,mkmem,&
&     mwan,natom,dtset%nband,nkpt,&
&     dtset%nspinor,nsppol,dtset%ntypat,nwan,pawrad,pawtab,&
&     proj_l,proj_m,proj_radial,proj_site,proj_x,proj_z,proj_zona,psps,&
&     rprimd,spin,dtset%typat,xred)
!    
     write(message, '(a,a,a,a)' ) ch10,&
&     '   mlwfovlp:  mlwfovlp_proj_paw done -',ch10,&
&     '   Inside-spheres part of projectors computed.'
     call wrtout(std_out,  message,'COLL')
!    
!    Add in-sphere contribution to A_matrix
!    
!    
!    w90iniprj==5. Plane waves + augmentation contributions
!    
     if(dtset%w90iniprj==5) A_matrix(:,:,:,:)=A_matrix(:,:,:,:)+A_paw(:,:,:,:)
!    
!    w90iniprj==6. Just augmentation contribution
!    
     if(dtset%w90iniprj==6) A_matrix(:,:,:,:)=A_paw(:,:,:,:)
!    
!    deallocations
!    
     ABI_DEALLOCATE(A_paw)
   end if !usepaw==1

   ABI_DEALLOCATE(just_augmentation)
!  
   call xmpi_barrier(spaceComm)
   call xmpi_sum(A_matrix,spaceComm,ierr)

!  
!  write      projections  to a file
!  
   if(rank==master) then
     if(dtset%w90iniprj==1) then
       do isppol=1,nsppol
         if(spin.ne.0 .and. spin.ne.isppol) cycle
         iun(isppol)=219+isppol
         open(unit=iun(isppol),file=trim(filew90_ramn(isppol)),form='formatted',status='unknown')
         write(iun(isppol),*) 'Projections from Abinit : mband,nkpt,nwan. indices: iband1,iwan,ikpt'
         write(iun(isppol),*) num_bands(isppol),nkpt,nwan(isppol)
       end do
     else
       do isppol=1,nsppol
         if(spin.ne.0 .and. spin.ne.isppol) cycle
         iun(isppol)=220+isppol
         open(unit=iun(isppol),file=trim(filew90_amn(isppol)),form='formatted',status='unknown')
         write(iun(isppol),*) 'Projections from Abinit : mband,nkpt,nwan. indices: iband1,iwan,ikpt'
         write(iun(isppol),*) num_bands(isppol),nkpt,nwan(isppol)
       end do
     end if
!    
     do isppol=1,nsppol
       if(spin.ne.0 .and. spin.ne.isppol) cycle
       do ikpt=1,nkpt
         do iwan=1,nwan(isppol)
           jband=0
           do iband=1,mband
             if(band_in(iband,isppol)) then
               jband=jband+1
               write(iun(isppol),'(3i6,13x,3x,2f18.14)')jband,iwan,ikpt,A_matrix(jband,iwan,ikpt,isppol)
             end if !band_in
           end do !iband
         end do !iwan
       end do !ikpt
     end do !isppol
!    
     if(dtset%w90iniprj==1) then
       do isppol=1,nsppol
         if(spin.ne.0 .and. spin.ne.isppol) cycle
         close(iun(isppol))
         write(message, '(3a)' ) &
&         '   ',trim(filew90_ramn(isppol)),' written'
         call wrtout(std_out,  message,'COLL')
       end do
     else
       do isppol=1,nsppol
         if(spin.ne.0 .and. spin.ne.isppol) cycle
         close(iun(isppol))
         write(message, '(3a)' ) &
&         '   ',trim(filew90_amn(isppol)),' written'
         call wrtout(std_out,  message,'COLL')
       end do
     end if
   end if !rank==master

!  
!  
!  Write down part of the matrix to the output file
!  This is for the automatic tests
!  
   if(rank==master) then
     write(message, '(4a)' ) ch10,&
&     '   Writing top of the initial projections matrix: A_mn(ik)',ch10,&
&     '   m=1:3, n=1:3, ik=1'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
!    
!    just write down the first 3 elements
!    
     do isppol=1,nsppol
       if(spin.ne.0 .and. spin.ne.isppol) cycle
       write(message, '( " " )')
       if (nsppol>1 ) then
         if (isppol==1) write(message,'(2a)')trim(message),'   spin up:'
         if (isppol==2) write(message,'(2a)')trim(message),'   spin down:'
       end if
       do ii=1,3
         if(ii>num_bands(isppol)) cycle
         write(message,'(3a)') trim(message),ch10,';   ( '
         do jj=1,3
           if(jj>nwan(isppol))cycle
           write(message, '(a,2f11.6,a)') trim(message),&
&           A_matrix(ii,jj,1,isppol),' , '
         end do
         write(message,'(2a)') trim(message),'    ) '
       end do
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end do
!    
!    Now write down bottom of the matrix
!    
     write(message, '(4a)' ) ch10,&
&     '   Writing bottom of the initial projections matrix: A_mn(ik)',ch10,&
&     '   m=num_bands-2:num_bands, n=nwan-2:nwan, ik=nkpt'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
!    
     do isppol=1,nsppol
       if(spin.ne.0 .and. spin.ne.isppol) cycle
       write(message, '( " " )')
       if (nsppol>1 ) then
         if (isppol==1) write(message,'(2a)')trim(message),'   spin up:'
         if (isppol==2) write(message,'(2a)')trim(message),'   spin down:'
       end if
       do ii=num_bands(isppol)-2,num_bands(isppol)
         if(ii<1) cycle
         write(message,'(3a)') trim(message),ch10,';   ( '
         do jj=nwan(isppol)-2,nwan(isppol)
           if(jj<1)cycle
           write(message, '(a,2f11.6,a)') trim(message),&
&           A_matrix(ii,jj,nkpt,isppol),' , '
         end do
         write(message,'(2a)') trim(message),'    ) '
       end do
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end do !isppol
   end if !rank==master
 end if !dtset%w90iniprj/=0
!
!Deallocations
!
 ABI_DEALLOCATE(proj_site)
 ABI_DEALLOCATE(proj_l)
 ABI_DEALLOCATE(proj_m)
 ABI_DEALLOCATE(proj_radial)
 ABI_DEALLOCATE(proj_x)
 ABI_DEALLOCATE(proj_z)
 ABI_DEALLOCATE(proj_zona)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!6) write files for wannier function plot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if( dtset%w90prtunk>0) then
   if(psps%usepaw==1) then
     write(message, '( a,a,a,a,a,a,a,a,a)')ch10,&
&     "   WARNING: The UNK matrices will not contain the correct wavefunctions ",ch10,&
&     "   since we are just writing the plane wave contribution.",ch10,&
&     "   The contribution from inside the spheres is missing. ",ch10,&
&     "   However, these files can be used for plotting purposes",ch10
     call wrtout(std_out,  message,'COLL')
   end if
!  
   spacing = dtset%w90prtunk
   write(message, '( 8a,i3,2a)')ch10,&
&   "   UNK files will be written.",ch10,&
&   "   According to the chosen value of w90prtunk",ch10,&
&   "   the wavefunctions are to be written ",ch10, &
&   "   at every ", spacing," records.",ch10
   call wrtout(std_out,  message,'COLL')
!  
   ABI_ALLOCATE(kg_k,(3,mpw))
   n1=ngfft(1)
   n2=ngfft(2)
   n3=ngfft(3)
   n4=ngfft(4)
   n5=ngfft(5)
   n6=ngfft(6)
   cplex=1
   mgfft=mgfftc ! error
   do isppol=1,nsppol
     ikg=0
     do ikpt=1,nkpt
!      
!      MPI:cycle over k-points not treated by this node
!      
       if (nprocs>1 ) then !sometimes we can have just one processor
         if ( ABS(MPI_enreg%proc_distrb(ikpt,1,isppol)-rank)  /=0) CYCLE
       end if
!      
       if(spin.eq.0 .or. spin.eq.isppol) then
         npw_k=npwarr(ikpt)
         kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
         ABI_ALLOCATE(denpot,(cplex*n4,n5,n6))
         ABI_ALLOCATE(cwavef,(2,npw_k))
         ABI_ALLOCATE(fofr,(2,n4,n5,n6))
         ABI_ALLOCATE(gbound,(2*mgfft+8,2))
         ABI_ALLOCATE(fofgout,(2,npw_k))
         iun_plot=1000+ikpt+ikpt*(isppol-1)
         write(wfnname,'("UNK",I5.5,".",I1)') ikpt, isppol
!        open (unit=iun_plot, file=wfnname,form='formatted')
         open(unit=iun_plot, file=wfnname,form='unformatted')
!        optimizing grid for UNK files
         n1tmp = n1/spacing
         n2tmp = n2/spacing
         n3tmp = n3/spacing
         if( mod(n1,spacing) /= 0) then
           n1tmp = n1tmp + 1
         end if
         if( mod(n2,spacing) /= 0) then
           n2tmp = n2tmp + 1
         end if
         if( mod(n3,spacing) /= 0) then
           n3tmp = n3tmp + 1
         end if
!        write(iun_plot,*) n1tmp,n2tmp,n3tmp,ikpt,nband_inc
         write(iun_plot) n1tmp,n2tmp,n3tmp,ikpt,nband_inc(isppol)
!        gbound=zero
         call sphereboundary(gbound,dtset%istwfk(ikpt),kg_k,mgfft,npw_k)
         write(std_out,*) "  writes UNK file for ikpt, spin=",ikpt,isppol
         denpot(:,:,:)=zero
         weight = one
         do iband=1,mband
           if(band_in(iband,isppol)) then
             do ig=1,npw_k*dtset%nspinor
               cwavef(1,ig)=cg(1,ig+iwav(iband,ikpt,isppol))
               cwavef(2,ig)=cg(2,ig+iwav(iband,ikpt,isppol))
             end do
             tim_fourwf=0
             call fourwf(cplex,denpot,cwavef,fofgout,fofr,&
&             gbound,gbound,dtset%istwfk(ikpt),kg_k,kg_k,mgfft,&
&             mpi_enreg,1,ngfft,npw_k,npw_k,n4,n5,n6,0,dtset%paral_kgb,&
&             tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)
!            do jj3=1,n3,spacing
!            do jj2=1,n2,spacing
!            do jj1=1,n1,spacing
!            write(iun_plot,*) fofr(1,jj1,jj2,jj3),&
!            & fofr(2,jj1,jj2,jj3)
!            end do !jj1
!            end do !jj2
!            end do !jj3
!            unformatted (must be one record)
             write(iun_plot) (((fofr(1,jj1,jj2,jj3),fofr(2,jj1,jj2,jj3),&
&             jj1=1,n1,spacing),jj2=1,n2,spacing),jj3=1,n3,spacing)
           end if !iband
         end do ! iband
         ABI_DEALLOCATE(cwavef)
         ABI_DEALLOCATE(fofr)
         ABI_DEALLOCATE(gbound)
         ABI_DEALLOCATE(denpot)
         ABI_DEALLOCATE(fofgout)
       end  if !spin
       ikg=ikg+npw_k
       if(spin.eq.0 .or. spin.eq.isppol) close(iun_plot)
     end do  ! ikpt
   end do  ! nsppol
   ABI_DEALLOCATE(kg_k)
!  
   write(message, '(4a)' )ch10, &
&   '   ','UNK files written',ch10
   call wrtout(std_out,  message,'COLL')
 end if !dtset%w90prtunk
!
 ABI_DEALLOCATE(iwav)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!7) Call to  Wannier90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if(lwannierrun) then
   if(lwanniersetup.ne.1) MSG_ERROR("lwanniersetup.ne.1")
   ABI_ALLOCATE(U_matrix,(mwan,mwan,nkpt,nsppol))
   ABI_ALLOCATE(U_matrix_opt,(max_num_bands,mwan,nkpt,nsppol))
   ABI_ALLOCATE(lwindow,(max_num_bands,nkpt,nsppol))
   ABI_ALLOCATE(wann_centres,(3,mwan,nsppol))
   ABI_ALLOCATE(wann_spreads,(mwan,nsppol))
!  Initialize
   U_matrix(:,:,:,:)=czero
   U_matrix_opt(:,:,:,:)=czero
   lwindow(:,:,:)=.false.
   wann_centres(:,:,:)=zero
   wann_spreads(:,:)=zero
!  
!  write(std_out,*) seed_name
!  write(std_out,*) ngkpt
   ngkpt(1)=dtset%kptrlatt(1,1)
   ngkpt(2)=dtset%kptrlatt(2,2) ! ajouter test de verif que kptrlatt est bien diagonal
   ngkpt(3)=dtset%kptrlatt(3,3)
!  write(std_out,*) nkpt
!  write(std_out,*) rprimd*Bohr_Ang
!  write(std_out,*) two_pi*gprimd/Bohr_Ang
!  write(std_out,*) mband
!  write(std_out,*) "nwan",nwan
!  write(std_out,*) nntot
!  write(std_out,*) natom
!  write(std_out,*) atom_symbols
!  write(std_out,*) xcart
!  write(std_out,*) num_bands,num_bands,nntot,nkpt
!  write(std_out,*) wann_spreads
!  wann_spreads=2
!  do i=1, nkpt
!  do j=1, nntot
!  write(std_out,*) i,j
!  do k=1, num_bands
!  do l=1, num_bands
!  write(std_out,*) "m",M_matrix(l,k,j,i,1)
!  enddo
!  enddo
!  enddo
!  enddo



#if defined HAVE_WANNIER90
   do isppol=1,nsppol
     if(spin.ne.0 .and. spin.ne.isppol) cycle
!    when nsppol>1, master runs isppol 1 and rank==1 runs isppol 2
     if(nprocs>1 .and. isppol==1.and.rank.ne.master) cycle
     if(nprocs>1 .and. isppol==2.and.rank.ne.1) cycle

     write(message, '(8a)' ) ch10,&
&     '** mlwfovlp :   call wannier90 library subroutine wannier_run ',ch10,&
&     '   Calculation is running         ',ch10,&
&     '-  see ',trim(filew90_wout(isppol)),' for details.'
     call wrtout(std_out,  message,'COLL')
!    
     call wannier_run(trim(seed_name(isppol)),ngkpt,nkpt,&            !input
&    real_lattice,recip_lattice,dtset%kpt,num_bands(isppol),& !input
&    nwan(isppol),nntot,natom,atom_symbols,&                  !input
&    xcart*Bohr_Ang,gamma_only,M_matrix(:,:,:,:,isppol),A_matrix(:,:,:,isppol),eigenvalues_w(:,:,isppol),& !input
&    U_matrix(1:nwan(isppol),1:nwan(isppol),:,isppol),& !output
&    U_matrix_opt(1:num_bands(isppol),1:nwan(isppol),:,isppol),& !output
&    lwindow_loc=lwindow(1:num_bands(isppol),:,isppol),& !output
&    wann_centres_loc=wann_centres(:,1:nwan(isppol),isppol),&     !output
&    wann_spreads_loc=wann_spreads(1:nwan(isppol),isppol),spread_loc=spreadw)                            !output

!    ----------------------------------------------------------------------------------------------

     write(message, '(7a)' ) ch10,&
&     '   mlwfovlp :  mlwfovlp_run completed -',ch10,&
&     '-  see ',trim(filew90_wout(isppol)),' for details.',ch10
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

   end do !isppol

!  collect output of  wannier90 from different processors
   call xmpi_barrier(spaceComm)

   call xmpi_sum(U_matrix,spaceComm,ierr)
   call xmpi_sum(U_matrix_opt,spaceComm,ierr)
   call xmpi_lor(lwindow,spaceComm)
   call xmpi_sum(wann_centres,spaceComm,ierr)
   call xmpi_sum(wann_spreads,spaceComm,ierr)


!  CALL SILVESTRELLI'S APPROACH TO EVALUATE vdW INTERACTION ENERGY USING MLWF!!
!  ----------------------------------------------------------------------------------------------
   if (dtset%vdw_xc==10.or.dtset%vdw_xc==11.or.dtset%vdw_xc==12.or.dtset%vdw_xc==14.and.rank==master) then 
!    vdw_xc==10,11,12,14 starts the 
!    vdW interaction using MLWFs
     write(std_out,*) 'nwan(nsppol)=',ch10
     do ii=1,nsppol
       write(std_out,*) 'nsppol=',ii, 'nwan(nsppol)=',nwan(ii),ch10
     end do
     write(std_out,*) 'mwan=', mwan, ch10

     ABI_ALLOCATE(occ_arr,(mband,nkpt,isppol))
     ABI_ALLOCATE(occ_wan,(mwan,nkpt,nsppol))
     ABI_ALLOCATE(tdocc_wan,(mwan,nsppol))

     occ_arr(:,:,:)=zero
     occ_wan(:,:,:)=zero      
     tdocc_wan(:,:)=zero 
     jj = 0
     do isppol=1,nsppol
       do ikpt=1,nkpt
         do iband=1,num_bands(isppol)   
           jj = jj + 1 
           occ_arr(iband,ikpt,isppol) = occ(jj)  
         end do
       end do
     end do 
     
     do isppol=1,nsppol
       do ikpt=1,nkpt
         do iwan=1,nwan(isppol)
           caux=czero
           caux2=czero
           caux3=czero
           do iband=1,num_bands(isppol) !nband_inc(isppol) !nwan(isppol)
             do ii=1,nwan(isppol) 
               caux=U_matrix(ii,iwan,ikpt,isppol)*U_matrix_opt(iband,ii,ikpt,isppol)
!              DEBUG
!              if(ISNAN(dble(caux))) then
!              write(std_out,*) 'NaN: caux(ikpt,iwan,iband,ii):',ikpt,iwan,iband,ii,ch10
!              end if
!              END DEBUG 
               do kk=1,nwan(isppol)
                 caux2=conjg(U_matrix(kk,iwan,ikpt,isppol))*conjg(U_matrix_opt(iband,kk,ikpt,isppol))
                 caux3= caux3+caux*caux2*occ_arr(iband,ikpt,isppol) !take care here as exclude_bands case is not well
!                DEBUG
!                if(ISNAN(dble(caux2))) then
!                write(std_out,*) 'NaN: caux2(ikpt,iwan,iband,kk):',ikpt,iwan,iband,kk,ch10
!                end if
!                if(ISNAN(dble(caux3))) then
!                write(std_out,*) 'NaN: caux3(ikpt,iwan,iband,kk,jj):',ikpt,iwan,iband,kk,jj
!                write(std_out,*) 'caux,caux2 and occ=',caux,caux2,occ(jj),ch10 
!                end if
!                END DEBUG 
               end do
             end do
           end do
           occ_wan(iwan,ikpt,isppol) = dble(caux3) 
!          DEBUG
!          write(std_out,*) occ_wan(iwan,ikpt,isppol)
!          END DEBUG      
!          end do
         end do
       end do
     end do

     write(std_out,*) ch10,'MLWFs Occupation Matrix diagonal terms:',ch10
     
     do jj=1,nsppol     
       forall(iwan=1:nwan(jj)) tdocc_wan(iwan,jj) = sum(occ_wan(iwan,1:nkpt,jj)) / real(nkpt,dp) 
       write(std_out,*) 'tdocc_wan(iwan),isppol:',ch10
       write(std_out,*) (tdocc_wan(iwan,jj),iwan=1,nwan(jj)),jj
     end do

     ABI_ALLOCATE(csix,(mwan,mwan,nsppol,nsppol))

     call evdw_wannier(csix,corrvdw,mwan,natom,nsppol,nwan,tdocc_wan,dtset%vdw_nfrag,&
&     dtset%vdw_supercell,dtset%vdw_typfrag,dtset%vdw_xc,rprimd,wann_centres,wann_spreads,xcart)

     ABI_DEALLOCATE(csix)
     ABI_DEALLOCATE(occ_arr)
     ABI_DEALLOCATE(occ_wan)
     ABI_DEALLOCATE(tdocc_wan)
   end if

#endif
!  FIXME: looks like there is no automatic test which goes through here: g95 bot did not catch
!  the missing deallocations
   ABI_DEALLOCATE(wann_centres)
   ABI_DEALLOCATE(wann_spreads)
   ABI_DEALLOCATE(U_matrix)
   ABI_DEALLOCATE(U_matrix_opt)
   ABI_DEALLOCATE(lwindow)

 end if !lwannierrun
!
!deallocation
!
 ABI_DEALLOCATE(band_in)
 ABI_DEALLOCATE(atom_symbols)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(eigenvalues_w)
 ABI_DEALLOCATE(M_matrix)
 ABI_DEALLOCATE(A_matrix)

 contains

!!***
!!****f* mlwfovlp/read_chkunit
!! NAME
!! read_chkunit
!!
!! FUNCTION
!! Function which reads the .chk file produced by Wannier90
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
 subroutine read_chkunit(seed_name,nkpt,ndimwin,ierr)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'read_chkunit'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt
 character(len=*),intent(in) :: seed_name
 integer,intent(out) :: ierr
!arrays
 integer,intent(out) :: ndimwin(nkpt)

!Local variables-------------------------------
!scalars
 integer :: chk_unit,ios,ikpt
 logical :: have_disentangled

!************************************************************************

   chk_unit=get_unit()
   fname=TRIM(seed_name)//'.chk'
   open(unit=chk_unit,file=fname,form='unformatted',status='old',iostat=ios)

   ierr=0
   read(chk_unit) ! header                                   ! Date and time
   read(chk_unit) ! ((real_lattice(i,j),i=1,3),j=1,3)        ! Real lattice
   read(chk_unit) ! ((recip_lattice(i,j),i=1,3),j=1,3)       ! Reciprocal lattice
   read(chk_unit) ! num_kpts
   read(chk_unit) ! ((kpt_latt(i,nkp),i=1,3),nkp=1,num_kpts) ! K-points
   read(chk_unit) ! nntot                  ! Number of nearest k-point neighbours
   read(chk_unit) ! num_wann               ! Number of wannier functions
   read(chk_unit) ! chkpt1                 ! Position of checkpoint
   read(chk_unit) have_disentangled        ! Whether a disentanglement has been performed
   if (have_disentangled) then
!    read(chk_unit) ! omega_invariant     ! Omega invariant
!    read(chk_unit) ((lwindow(i,nkp),i=1,num_bands),nkp=1,num_kpts)
     read(chk_unit) (ndimwin(ikpt),ikpt=1,nkpt)
!    read(chk_unit) (((u_matrix_opt(i,j,nkp),i=1,num_bands),j=1,num_wann),nkp=1,num_kpts)
   else
!    this is not expected. we should have disentanglement. Report the error.
     ierr=-1
   end if
!  read(chk_unit)  (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann),k=1,num_kpts)               ! U_matrix
!  read(chk_unit)  ((((m_matrix(i,j,k,l,1),i=1,num_wann),j=1,num_wann),k=1,nntot),l=1,num_kpts) ! M_matrix
!  read(chk_unit)  ((wannier_centres(i,j),i=1,3),j=1,num_wann)
   close(chk_unit)

 end subroutine read_chkunit

end subroutine mlwfovlp
!!***
