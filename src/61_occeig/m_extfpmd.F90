!!****m* ABINIT/m_extfpmd
!! NAME
!! m_extfpmd
!!
!! FUNCTION
!! This module provides routines to run computations at very high temperature
!! with reduced number of bands. High energy orbitals are represented as
!! pure plane waves.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2022 ABINIT group (A. Blanchet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!! 1) Add contribution to conductivity.
!! 2) Smooth the contributions.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_extfpmd
  use defs_basis
  use defs_abitypes
  use m_io_tools
  use m_errors
  use m_geometry
  use m_special_funcs
  use m_specialmsg
  use m_xmpi
  use m_dtset,          only : dataset_type
  use m_hdr,            only : hdr_type, hdr_init
  use m_cgtools,        only : pw_orthon, cgnc_cholesky
  use m_energies,       only : energies_type
  use m_gsphere,        only : getkpgnorm
  use m_kg,             only : kpgio, getmpw, mkkin
  use m_mpinfo,         only : ptabs_fourdp, proc_distrb_cycle, destroy_mpi_enreg, initmpi_seq, copy_mpi_enreg, distrb2
  use m_numeric_tools,  only : simpson,simpson_int
  use m_spacepar,       only : meanvalue_g

  implicit none
  public :: extfpmd_dos,extfpmd_e_fg,extfpmd_i_fg,extfpmd_chkinp
  !!***

  !----------------------------------------------------------------------

  !!****t* m_extfpmd/extfpmd_type
  !! NAME
  !! extfpmd_type
  !!
  !! FUNCTION
  !! Store extfpmd functions and parameters.
  !!
  !! SOURCE
  type,public :: extfpmd_type
    logical :: truecg
    integer :: bcut,nbcut,nbdbuf,nfft,nspden,version
    real(dp) :: ebcut,edc_kinetic,e_kinetic,entropy
    real(dp) :: nelect,eshift,ucvol,bandshift
    real(dp),allocatable :: vtrial(:,:)
    real(dp),allocatable :: nelectarr(:,:)
    real(dp),allocatable :: bandshiftk(:)
    !! Scalars and arrays for numerical extended PW method
    integer :: mpw,mcg,mband,mkmem
    real(dp) :: ecut,ecut_eff
    type(hdr_type) :: hdr
    type(MPI_type) :: mpi_enreg
    integer,allocatable :: kg(:,:),npwarr(:),npwtot(:),nband(:)
    real(dp),allocatable :: cg(:,:),eigen(:),occ(:),doccde(:),resid(:)
  contains
    procedure :: compute_e_kinetic
    procedure :: compute_entropy
    procedure :: compute_nelect
    procedure :: generate_extpw
    procedure :: compute_eshift
    procedure :: extpw_orthon
    procedure :: init
    procedure :: destroy
  end type extfpmd_type
  !!***

contains

  !!****f* ABINIT/m_extfpmd/init
  !! NAME
  !!  init
  !!
  !! FUNCTION
  !!  Initialize extfpmd_type object, memory allocation of arrays...
  !!
  !! INPUTS
  !!  this=extfpmd_type object concerned
  !!  mband=maximum number of bands
  !!  nbcut=number of states used to average the constant potential value
  !!  nbdbuf=Number of bands in the buffer to converge scf cycle with extfpmd models
  !!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
  !!  version=extfpmd implementation version
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! SOURCE
  subroutine init(this,mband,extpw_eshift,nbcut,nbdbuf,nfft,nspden,nsppol,nkpt,rprimd,&
  & version,ecut,exchn2n3d,istwfk,kptns,mpi_enreg,mkmem,dilatmx,extfpmd_ecut,extfpmd_mband,&
  & nspinor,truecg)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: mband,nbcut,nbdbuf,nfft,nspden,nsppol,nkpt,version
    integer,intent(in) :: exchn2n3d,mkmem,extfpmd_mband,nspinor,truecg
    real(dp),intent(in) :: ecut,extfpmd_ecut,extpw_eshift,dilatmx
    type(MPI_type),intent(inout) :: mpi_enreg
    ! Arrays
    integer,intent(in) :: istwfk(nkpt)
    real(dp),intent(in) :: rprimd(3,3)
    real(dp),intent(in) :: kptns(3,nkpt)

    ! Local variables -------------------------
    ! Scalars
    integer :: nb_per_proc,my_nspinor
    ! Arraysz
    real(dp) :: gprimd(3,3),rmet(3,3),gmet(3,3)

    ! *********************************************************************
    ! write(0,*) "DEBUG: Init extfpmd object"
    
    this%truecg=.false.
    this%bcut=mband-nbdbuf
    this%nbcut=nbcut
    this%mband=extfpmd_mband
    this%nbdbuf=nbdbuf
    this%version=version
    ABI_MALLOC(this%vtrial,(nfft,nspden))
    this%vtrial(:,:)=zero
    this%nfft=nfft
    this%nspden=nspden
    this%ebcut=zero
    this%edc_kinetic=zero
    this%e_kinetic=zero
    this%entropy=zero
    this%nelect=zero
    ABI_MALLOC(this%nelectarr,(nfft,nspden))
    this%nelectarr(:,:)=zero
    this%bandshift=zero
    ABI_MALLOC(this%bandshiftk,(nkpt*nsppol))
    this%bandshiftk(:)=zero
    this%eshift=extpw_eshift
    call metric(gmet,gprimd,-1,rmet,rprimd,this%ucvol)
    
    if(this%version==11) then
      if(truecg==1) this%truecg=.true.
      if(extfpmd_ecut==zero) then
        ! Adding dsqrt(...) to extended cutoff energy to make
        ! sure extended pw basis set is large enough.
        this%ecut=extfpmd_e_fg(one*extfpmd_mband,this%ucvol)+&
        & dsqrt((2*PI*gprimd(1,1))**2+(2*PI*gprimd(2,1))**2+(2*PI*gprimd(3,1))**2)
        ! Force extended plane wave kinetic energy to be geq Kohn-Sham ecut.
        this%ecut=max(this%ecut,ecut)
      else
        ! Automatically determine extended plane waves energy cutoff
        this%ecut=extfpmd_ecut
      end if
      this%ecut_eff=this%ecut*dilatmx**2
      ABI_MALLOC(this%nband,(nkpt*nsppol))
      this%nband(:)=extfpmd_mband
      ABI_MALLOC(this%eigen,(extfpmd_mband*nkpt*nsppol))
      this%eigen(:)=zero
      ABI_MALLOC(this%resid,(extfpmd_mband*nkpt*nsppol))
      this%resid(:)=zero
      ABI_MALLOC(this%occ,(extfpmd_mband*nkpt*nsppol))
      this%occ(:)=zero
      ABI_MALLOC(this%doccde,(extfpmd_mband*nkpt*nsppol))
      this%doccde(:)=zero
      
      ! Initialize mpi_enreg variable for seq use.
      this%mkmem=nkpt
      call initmpi_seq(this%mpi_enreg)
      ABI_MALLOC(this%mpi_enreg%proc_distrb,(nkpt,this%mband,nsppol))
      this%mpi_enreg%proc_distrb=0
      this%mpi_enreg%me_g0=1

      ! WARNING: using sequential full array for this%kg. Otherwise: memory leak when using 
      ! kpgsph subroutine in m_fftcore "integer, save :: alloc_size=0" ?
      call getmpw(this%ecut_eff,exchn2n3d,gmet,istwfk,kptns,mpi_enreg,this%mpw,nkpt)
      ABI_MALLOC(this%kg,(3,this%mpw*nkpt))
      this%kg(:,:)=0
      ABI_MALLOC(this%npwarr,(nkpt))
      this%npwarr(:)=0
      ABI_MALLOC(this%npwtot,(nkpt))
      this%npwtot(:)=0
      call kpgio(this%ecut_eff,exchn2n3d,gmet,istwfk,this%kg, &
        & kptns,this%mkmem,this%nband,nkpt,'PERS',this%mpi_enreg,&
        & this%mpw,this%npwarr,this%npwtot,nsppol)
      call destroy_mpi_enreg(this%mpi_enreg)
      
      ! Setup MPI parallelization
      this%mkmem=mkmem
      call copy_mpi_enreg(mpi_enreg,this%mpi_enreg)
      ! Free proc_distrib and reallocate it with extended values
      ABI_FREE(this%mpi_enreg%proc_distrb)
      ABI_MALLOC(this%mpi_enreg%proc_distrb,(nkpt,this%mband,nsppol))
      call distrb2(this%mband,nb_per_proc,this%nband,nkpt,this%mpi_enreg%nproc,nsppol,this%mpi_enreg)
      
      ! Get number of ext pw coefficients
      my_nspinor=max(1,nspinor/this%mpi_enreg%nproc_spinor)
      this%mcg=this%mpw*my_nspinor*this%mband*this%mkmem*nsppol
      ABI_MALLOC(this%cg,(2,this%mcg))
      this%cg(:,:)=zero
    else if(this%version==5) then
      ! Make a copy of mpi_enreg in order to cycle.
      call copy_mpi_enreg(mpi_enreg,this%mpi_enreg)
    end if

  end subroutine init
  !!***

  !!****f* ABINIT/m_extfpmd/destroy
  !! NAME
  !!  destroy
  !!
  !! FUNCTION
  !!  Destroy extfpmd_type object, memory deallocation of arrays...
  !!
  !! INPUTS
  !!  this=extfpmd_type object concerned
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! SOURCE
  subroutine destroy(this)

    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this

    ! *********************************************************************
    if(this%version==11) then
      this%cg(:,:)=zero
      ABI_FREE(this%cg)
      this%npwtot(:)=0
      ABI_FREE(this%npwtot)
      this%npwarr(:)=0
      ABI_FREE(this%npwarr)
      this%kg(:,:)=0
      ABI_FREE(this%kg)
      this%doccde(:)=zero
      ABI_FREE(this%doccde)
      this%occ(:)=zero
      ABI_FREE(this%occ)
      this%eigen(:)=zero
      ABI_FREE(this%eigen)
      this%nband(:)=0
      ABI_FREE(this%nband)
      this%ecut=zero
      this%ecut_eff=zero
      this%mkmem=0
      this%mpw=0
      this%truecg=.false.
    end if
    
    if(this%version==5.or.this%version==11) then
      call destroy_mpi_enreg(this%mpi_enreg)
    end if
    
    this%vtrial(:,:)=zero
    ABI_FREE(this%vtrial)
    this%nelectarr(:,:)=zero
    ABI_FREE(this%nelectarr)
    this%bandshiftk(:)=zero
    ABI_FREE(this%bandshiftk)
    this%nfft=0
    this%nspden=0
    this%bcut=0
    this%mband=0
    this%nbcut=0
    this%nbdbuf=0
    this%version=1
    this%ebcut=zero
    this%edc_kinetic=zero
    this%e_kinetic=zero
    this%entropy=zero
    this%nelect=zero
    this%bandshift=zero
    this%eshift=zero
    this%ucvol=zero
  end subroutine destroy
  !!***

  !!****f* ABINIT/m_extfpmd/compute_eshift
  !! NAME
  !!  compute_eshift
  !!
  !! FUNCTION
  !!  Computes the energy shift factor $U_0$ corresponding to constant
  !!  potential contribution.
  !!
  !! INPUTS
  !!  this=extfpmd_type object concerned
  !!  eigen(mband*nkpt*nsppol)=eigenvalues (hartree)
  !!  eknk(mband*nkpt*nsppol)=kinetic energies (hartree)
  !!  mband=maximum number of bands
  !!  nband(nkpt*nsppol)=desired number of bands at each k point
  !!  nfft=(effective) number of FFT grid points (for this processor)
  !!  nkpt=number of k points
  !!  nsppol=1 for unpolarized, 2 for spin-polarized
  !!  nspden=number of spin-density components
  !!  wtk(nkpt)=k point weights
  !!  vtrial(nfftf,nspden)=GS potential (Hartree)
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! SOURCE
  subroutine compute_eshift(this,eigen,eknk,mband,me,nband,nfft,nkpt,nsppol,nspden,wtk,vtrial)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: mband,me,nfft,nkpt,nsppol,nspden
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
    real(dp),intent(in) :: eknk(mband*nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)
    real(dp),intent(in) :: vtrial(nfft,nspden)

    ! Local variables -------------------------
    ! Scalars
    integer :: band_index,ii,ikpt,isppol,nband_k

    ! *********************************************************************
    this%vtrial=vtrial

    if(this%version==2) then
      ! Computes U_0^{HEG} from the difference between
      ! eigenvalues and Fermi gas energies, averaged
      ! over lasts nbcut bands.
      this%eshift=zero
      band_index=0
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          do ii=nband_k-this%nbdbuf-this%nbcut+1,nband_k-this%nbdbuf
            this%eshift=this%eshift+&
            & wtk(ikpt)*(eigen(band_index+ii)-extfpmd_e_fg(dble(ii),this%ucvol))
          end do
          band_index=band_index+nband_k
        end do
      end do
      this%eshift=this%eshift/this%nbcut
    else if(this%version==3) then
      ! Computes U_0^K from the difference between
      ! eigenvalues and kinetic energies, averaged
      ! over lasts nbcut bands.
      this%eshift=zero
      band_index=0
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          do ii=nband_k-this%nbdbuf-this%nbcut+1,nband_k-this%nbdbuf
            this%eshift=this%eshift+&
            & wtk(ikpt)*(eigen(band_index+ii)-eknk(band_index+ii))
          end do
          band_index=band_index+nband_k
        end do
      end do
      this%eshift=this%eshift/this%nbcut
    else
      ! Computes U_0 from the sum of local
      ! potentials (vtrial), averaging over all space.
      ! Simplest and most precise way to evaluate U_0.
      this%eshift=sum(this%vtrial)/(this%nfft*this%nspden)
    end if

    ! Get extended FPMD band energy cutoff
    this%ebcut=zero
    this%bandshift=zero
    this%bandshiftk(:)=zero
    band_index=0
    do isppol=1,nsppol
      do ikpt=1,nkpt
        nband_k=nband(ikpt+(isppol-1)*nkpt)
        this%ebcut=this%ebcut+wtk(ikpt)*eigen(band_index+nband_k-this%nbdbuf)/nsppol
        this%bandshift=this%bandshift+wtk(ikpt)*&
        & (extfpmd_i_fg(eigen(band_index+nband_k-this%nbdbuf)-this%eshift,this%ucvol)-(nband_k-this%nbdbuf))/nsppol
        this%bandshiftk(ikpt+(isppol-1)*nkpt)=extfpmd_i_fg(eigen(band_index+nband_k-this%nbdbuf)-this%eshift,this%ucvol)-(nband_k-this%nbdbuf)
        band_index=band_index+nband_k
      end do
    end do
  end subroutine compute_eshift
  !!***eigen(band_index+nband_k-this%nbdbuf)

  !!****f* ABINIT/m_extfpmd/generate_extpw
  !! NAME
  !!  generate_extpw
  !!
  !! FUNCTION
  !!  Computes the extended plane waves basis set vectors
  !!
  !! INPUTS
  !!  this=extfpmd_type object concerned
  !!  exchn2n3d=if 1, n2 and n3 are exchanged
  !!  gmet(3,3)=reciprocal space metric (bohr^-2)
  !!  istwfk(nkpt)=input option parameter that describes the storage of wfs
  !!  kptns(3,nkpt)=reduced coords of k points
  !!  mkmem =number of k points treated by this node.
  !!  character(len=4) : mode_paral=either 'COLL' or 'PERS', tells whether
  !!   the loop over k points must be done by all processors or not,
  !!   in case of parallel execution.
  !!  mpi_enreg=information about MPI parallelization
  !!  mpw=maximum number of planewaves as dimensioned in calling routine
  !!  nband(nkpt*nsppol)=number of bands at each k point
  !!  nkpt=number of k points
  !!  nsppol=1 for unpolarized, 2 for polarized
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! SOURCE
  subroutine generate_extpw(this,exchn2n3d,effmass_free,gmet,istwfk,kptns,mkmem,nband,nkpt,&
    & mode_paral,mpi_enreg,nsppol,dilatmx,nspinor,cg,mcg,npwarr,kg,mpw,eigen,mband,ecut,ecutsm,&
    & resid)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: exchn2n3d,mkmem,nkpt,nsppol,nspinor,mcg,mpw,mband
    real(dp),intent(in) :: dilatmx,effmass_free,ecut,ecutsm
    character(len=4),intent(in) :: mode_paral
    type(MPI_type),intent(inout) :: mpi_enreg
    ! Arrays
    integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt),kg(3,mpw*mkmem)
    real(dp),intent(in) :: gmet(3,3),kptns(3,nkpt)
    real(dp),intent(in) :: cg(2,mcg),eigen(mband*nkpt*nsppol),resid(mband*nkpt*nsppol)
    
    ! Local variables -------------------------
    ! Scalars
    integer :: isppol,ikpt,iband,icg,ext_icg,ipw,ext_ipw,ikg,ext_ikg,ierr
    integer :: nband_k,ext_nband_k,istwf_k,npw_k,ext_npw_k
    integer :: my_nspinor,bdtot_index,ext_bdtot_index
    integer :: index_below,index_above,g0_count_below,g0_count_above
    real(dp) :: ecut_eff,ekin_max,fg_kin,closest_below,closest_above,prop_below,prop_above
    real(dp) :: norm,nband_k_kin,count_below,count_above
    real(dp) :: tmp,bandshift,dotr,phase,random_value
    ! Arrays
    real(dp) :: kpoint(3)
    integer,allocatable :: kg_k(:,:),ext_kg_k(:,:),indices_below(:),indices_above(:)
    real(dp),allocatable :: ext_kinpw(:)

    ! *********************************************************************
    
    ! write(0,*) mpi_enreg%me_kpt,'DEBUG: Generating extended plane wave basis set...'
    this%eigen(:)=zero
    this%resid(:)=zero
    my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)
    
    ! Get ext pw coefficients and kinetic and eigenvalues
    ! Loop over spins
    bdtot_index=0
    ext_bdtot_index=0
    icg=0
    ext_icg=0
    ! write(0,*) "DEBUG: Starting loops"
    do isppol=1,nsppol
      ikg=0
      ext_ikg=0
      ! Loop over k points
      do ikpt=1,nkpt
        nband_k=nband(ikpt+(isppol-1)*nkpt)
        ext_nband_k=this%mband
        istwf_k=istwfk(ikpt)
        npw_k=npwarr(ikpt)
        ext_npw_k=this%npwarr(ikpt)
        ! Skip this k-point if not the proper processor
        if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,mpi_enreg%me_kpt)) then
          bdtot_index=bdtot_index+nband_k
          ext_bdtot_index=ext_bdtot_index+this%mband
          cycle
        end if
        ! Temporary allocate k-point dependant arrays
        ABI_MALLOC(kg_k,(3,npw_k))
        ABI_MALLOC(ext_kg_k,(3,ext_npw_k))
        ABI_MALLOC(ext_kinpw,(ext_npw_k))

        kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
        ext_kg_k(:,1:ext_npw_k)=this%kg(:,1+ext_ikg:ext_npw_k+ext_ikg)
        kpoint(:)=kptns(:,ikpt)
        
        ! Copy non extended plane waves coefficients and eigenvalues here
        do iband=1,nband_k-this%nbdbuf
          do ipw=1,npw_k*my_nspinor
            do ext_ipw=1,ext_npw_k*my_nspinor
              if(all(kg_k(:,ipw)==ext_kg_k(:,ext_ipw))) then
                this%cg(:,ext_ipw+(iband-1)*ext_npw_k*my_nspinor+ext_icg)=cg(:,ipw+(iband-1)*npw_k*my_nspinor+icg)
              end if
            end do
          end do
          this%eigen(iband+ext_bdtot_index)=eigen(iband+bdtot_index)
          this%resid(iband+ext_bdtot_index)=resid(iband+bdtot_index)
        end do
        
        ! bandshift=this%bandshift ! Converges faster, but transition less smooth
        bandshift=this%bandshiftk(ikpt+(isppol-1)*nkpt)
        
        ! Compute kinetic energy of extended pw.
        ! /!\ We don't keep this in memory because it can be heavy.
        ext_kinpw(:)=zero
        call mkkin(this%ecut,zero,effmass_free,gmet,ext_kg_k,ext_kinpw,kpoint,ext_npw_k,0,0)
        ! Set extended plane waves coefficients here
        do iband=nband_k-this%nbdbuf+1,this%mband
          ! Get fermi gas energy and set eigenvalue
          fg_kin=extfpmd_e_fg(one*iband+bandshift,this%ucvol)
          this%eigen(iband+ext_bdtot_index)=fg_kin+this%eshift

          ! Find closest values in ext_kinpw
          closest_below=zero
          closest_above=this%ecut
          count_below=zero
          count_above=zero
          g0_count_below=0
          g0_count_above=0

          do ext_ipw=1,ext_npw_k*my_nspinor
            if(ext_kinpw(ext_ipw)<=fg_kin.and.ext_kinpw(ext_ipw)>=closest_below) then
              if(ext_ipw==1) then
                g0_count_below=1
              else
                g0_count_below=0
              end if
              if(ext_kinpw(ext_ipw)==closest_below) then
                count_below=count_below+one
              else
                count_below=one
                closest_below=ext_kinpw(ext_ipw)
              end if
            else if(ext_kinpw(ext_ipw)>=fg_kin.and.ext_kinpw(ext_ipw)<=closest_above) then
              if(ext_ipw==1) then
                g0_count_above=1
              else
                g0_count_above=0
              end if
              if(ext_kinpw(ext_ipw)==closest_above) then
                count_above=count_above+one
              else
                count_above=one
                closest_above=ext_kinpw(ext_ipw)
              end if
            end if
          end do
          
          ! Taking time-reversal symmetry into account excluding g0.
          if(istwf_k >= 2) then
            count_above=count_above*(two-g0_count_above)
            count_below=count_below*(two-g0_count_below)
          end if
          prop_below=(fg_kin-closest_above)/(count_below*(closest_below-closest_above))
          prop_above=(one-prop_below*count_below)/count_above
          
          ! Do a second loop to set cg coefficients (and reset old ones).
          ! Generating random phase between 0 and 2pi.
          ! random_value=zero
          do ext_ipw=1,ext_npw_k*my_nspinor
            ! call random_number(random_value)
            phase=random_value*two*PI
            this%cg(:,ext_ipw+(iband-1)*ext_npw_k*my_nspinor+ext_icg)=zero
            if(ext_kinpw(ext_ipw)==closest_below) then
              this%cg(1,ext_ipw+(iband-1)*ext_npw_k*my_nspinor+ext_icg)=dsqrt(prop_below)*cos(phase)
              this%cg(2,ext_ipw+(iband-1)*ext_npw_k*my_nspinor+ext_icg)=dsqrt(prop_below)*sin(phase)
            end if
            if(ext_kinpw(ext_ipw)==closest_above) then
              this%cg(1,ext_ipw+(iband-1)*ext_npw_k*my_nspinor+ext_icg)=dsqrt(prop_above)*cos(phase)
              this%cg(2,ext_ipw+(iband-1)*ext_npw_k*my_nspinor+ext_icg)=dsqrt(prop_above)*sin(phase)
            end if
          end do
        end do

        ! write(0,*) "DEBUG: Checking extended plane waves coefficients normalization"
        ! do iband=1,this%mband
        !   fg_kin=extfpmd_e_fg(one*iband+bandshift,this%ucvol)
        !   norm=zero
        !   do ext_ipw=1,ext_npw_k*my_nspinor
        !     norm=norm+(this%cg(1,ext_ipw+(iband-1)*ext_npw_k*my_nspinor+ext_icg)**2+this%cg(2,ext_ipw+(iband-1)*ext_npw_k*my_nspinor+ext_icg)**2)
        !   end do
        !   write(99,*) ikpt, iband, norm
        ! end do

        ! Increment indexes
        bdtot_index=bdtot_index+nband_k
        ext_bdtot_index=ext_bdtot_index+this%mband
        if(mkmem/=0) then
          icg=icg+npw_k*my_nspinor*nband_k
          ext_icg=ext_icg+ext_npw_k*my_nspinor*this%mband
          ikg=ikg+npw_k
          ext_ikg=ext_ikg+ext_npw_k
        end if

        ABI_FREE(ext_kinpw)
        ABI_FREE(ext_kg_k)
        ABI_FREE(kg_k)
      end do
    end do
    ! write(0,*) mpi_enreg%me_kpt,'DEBUG: End...'
    call xmpi_sum(this%eigen,xmpi_world,ierr)
  end subroutine generate_extpw

  !!****f* ABINIT/m_extfpmd/extpw_orthon
  !! NAME
  !!  extpw_orthon
  !!
  !! FUNCTION
  !!  Orthonormalize extended plane wave wave functions
  !!
  !! INPUTS
  !!  this=extfpmd_type object concerned
  !!  exchn2n3d=if 1, n2 and n3 are exchanged
  !!  gmet(3,3)=reciprocal space metric (bohr^-2)
  !!  istwfk(nkpt)=input option parameter that describes the storage of wfs
  !!  kptns(3,nkpt)=reduced coords of k points
  !!  mkmem =number of k points treated by this node.
  !!  character(len=4) : mode_paral=either 'COLL' or 'PERS', tells whether
  !!   the loop over k points must be done by all processors or not,
  !!   in case of parallel execution.
  !!  mpi_enreg=information about MPI parallelization
  !!  mpw=maximum number of planewaves as dimensioned in calling routine
  !!  nband(nkpt*nsppol)=number of bands at each k point
  !!  nkpt=number of k points
  !!  nsppol=1 for unpolarized, 2 for polarized
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! SOURCE
  subroutine extpw_orthon(this,effmass_free,gmet,istwfk,kptns,mkmem,nkpt,&
    & nsppol,nspinor,usepaw)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: mkmem,nkpt,nsppol,nspinor,usepaw
    real(dp),intent(in) :: effmass_free
    ! Arrays
    integer,intent(in) :: istwfk(nkpt)
    real(dp),intent(in) :: gmet(3,3),kptns(3,nkpt)
    
    ! Local variables -------------------------
    ! Scalars
    integer :: isppol,ikpt,iband,ext_icg,ext_ikg,ierr
    integer :: ext_nband_k,istwf_k,ext_npw_k,ortalgo
    integer :: my_nspinor,ext_bdtot_index
    real(dp) :: dotr
    ! Arrays
    real(dp) :: kpoint(3)
    integer,allocatable :: ext_kg_k(:,:)
    real(dp),allocatable :: ext_kinpw(:)
    real(dp),allocatable :: ext_cwavef(:,:)

    ! *********************************************************************
    
    write(0,*) "DEBUG: Orthogonalization...."
    ortalgo=this%mpi_enreg%paral_kgb
    my_nspinor=max(1,nspinor/this%mpi_enreg%nproc_spinor)
    
    ! Get ext pw coefficients and kinetic and eigenvalues
    ext_bdtot_index=0
    ext_icg=0
    ! Loop over spins
    do isppol=1,nsppol
      ext_ikg=0
      ! Loop over k points
      do ikpt=1,nkpt
        ext_nband_k=this%mband
        istwf_k=istwfk(ikpt)
        ext_npw_k=this%npwarr(ikpt)
        ! Skip this k-point if not the proper processor
        if(proc_distrb_cycle(this%mpi_enreg%proc_distrb,ikpt,1,this%mband,isppol,this%mpi_enreg%me_kpt)) then
          ext_bdtot_index=ext_bdtot_index+this%mband
          cycle
        end if
        ! Temporary allocate k-point dependant arrays
        ABI_MALLOC(ext_kg_k,(3,ext_npw_k))
        ABI_MALLOC(ext_kinpw,(ext_npw_k))
        ABI_MALLOC(ext_cwavef,(2,ext_npw_k*my_nspinor*this%mband))

        ext_kg_k(:,1:ext_npw_k)=this%kg(:,1+ext_ikg:ext_npw_k+ext_ikg)
        kpoint(:)=kptns(:,ikpt)
        ! Compute kinetic energy of extended pw.
        ! /!\ We don't keep this in memory because it can be heavy.
        ext_kinpw(:)=zero
        call mkkin(this%ecut,zero,effmass_free,gmet,ext_kg_k,ext_kinpw,kpoint,ext_npw_k,0,0)
        
        call pw_orthon(0,0,istwf_k,this%mcg,this%mcg,ext_npw_k*my_nspinor,this%mband,4,this%cg,usepaw,this%cg,&
        & this%mpi_enreg%me_g0,this%mpi_enreg%comm_bandspinorfft)
        ! call cgnc_cholesky(ext_npw_k*my_nspinor,this%mband,this%cg,istwf_k,this%mpi_enreg%me_g0,this%mpi_enreg%comm_bandspinorfft,use_gemm=.False.)
        
        ! Update extpw eigen values with orthogonalized ones (should not be different)
        do iband=1,this%mband
          ext_cwavef(1:2,1:ext_npw_k*my_nspinor)= &
          & this%cg(:,1+(iband-1)*ext_npw_k*my_nspinor+ext_icg:iband*ext_npw_k*my_nspinor+ext_icg)
          call meanvalue_g(dotr,ext_kinpw,0,istwf_k,this%mpi_enreg,ext_npw_k,my_nspinor,ext_cwavef,ext_cwavef,0)
          this%eigen(iband+ext_bdtot_index)=dotr+this%eshift
        end do

        ! Increment indexes
        ext_bdtot_index=ext_bdtot_index+this%mband
        if(mkmem/=0) then
          ext_icg=ext_icg+ext_npw_k*my_nspinor*this%mband
          ext_ikg=ext_ikg+ext_npw_k
        end if

        ABI_FREE(ext_cwavef)
        ABI_FREE(ext_kinpw)
        ABI_FREE(ext_kg_k)
      end do
    end do
    call xmpi_sum(this%eigen,xmpi_world,ierr)
    write(0,*) this%mpi_enreg%me_kpt,'DEBUG: End...'
  end subroutine extpw_orthon

  !!****f* ABINIT/m_extfpmd/compute_nelect
  !! NAME
  !!  compute_nelect
  !!
  !! FUNCTION
  !!  Computes the value of the integral corresponding to the missing
  !!  free electrons contribution after band cut, with an order 1/2
  !!  incomplete Fermi-Dirac integral.
  !!
  !! INPUTS
  !!  this=extfpmd_type object concerned
  !!  fermie=chemical potential (Hartree)
  !!  nelect=number of electrons per unit cell
  !!  tsmear=smearing width (or temperature)
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! SOURCE
  subroutine compute_nelect(this,fermie,nband,nelect,nkpt,nspinor,nsppol,tsmear,wtk)
    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: nkpt,nsppol,nspinor
    real(dp),intent(in) :: fermie,tsmear
    real(dp),intent(inout) :: nelect
    class(extfpmd_type),intent(inout) :: this
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)

    ! Local variables -------------------------
    ! Scalars
    integer :: ifft,ispden,isppol,ikpt,iband,ext_bdtot_index,nband_k,ierr
    real(dp) :: factor,gamma,xcut,fn,maxocc
    ! Arrays
    real(dp),allocatable :: gamma_hybrid_tf(:,:)
    real(dp),allocatable :: xcut_hybrid_tf(:,:)

    ! *********************************************************************
    
    maxocc=two/(nsppol*nspinor)
    factor=dsqrt(two)/(PI*PI)*this%ucvol*tsmear**(1.5)
    gamma=(fermie-this%eshift)/tsmear

    ! Computes extfpmd contribution to nelect integrating
    ! over accessible states from bcut to infinity with
    ! order 1/2 incomplete Fermi-Dirac integral.
    if(this%version==2.or.this%version==4) then
      xcut=extfpmd_e_fg(one*this%bcut+this%bandshift,this%ucvol)/tsmear
      if(one*this%bcut+this%bandshift.lt.zero) xcut=zero
      nelect=nelect+factor*djp12(xcut,gamma)
    end if

    ! Computes extfpmd contribution to nelect integrating
    ! over energy from ebcut to infinity with order 1/2
    ! incomplete Fermi-Dirac integral.
    if(this%version==1.or.this%version==3) then
      xcut=(this%ebcut-this%eshift)/tsmear
      if(this%ebcut.lt.this%eshift) xcut=zero
      nelect=nelect+factor*djp12(xcut,gamma)
    end if
    
    ! Computes extfpmd contribution to nelect summing
    ! over accessible states from bcut to infinity, with
    ! integer band numbers. Total number of bands
    ! is controlled with the input variable extfpmd_nband.
    if(this%version==5) then
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          if(proc_distrb_cycle(this%mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,this%mpi_enreg%me_kpt)) cycle
          do iband=nband_k-this%nbdbuf+1,this%mband
            fn=fermi_dirac(extfpmd_e_fg(one*iband+this%bandshiftk(ikpt+(isppol-1)*nkpt),this%ucvol)+this%eshift,fermie,tsmear)
            nelect=nelect+wtk(ikpt)*maxocc*fn
          end do
        end do
      end do
      call xmpi_sum(nelect,xmpi_world,ierr)
    end if

    ! Computes extended pw contribution to nelect summing
    ! over ext pw states from nband_k to this%mband
    if(this%version==11) then
      ext_bdtot_index=0
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          ! Skip this k-point if not the proper processor
          if(proc_distrb_cycle(this%mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,this%mpi_enreg%me_kpt)) then
            ext_bdtot_index=ext_bdtot_index+this%mband
            cycle
          end if
          do iband=nband_k-this%nbdbuf+1,this%mband
            nelect=nelect+wtk(ikpt)*this%occ(iband+ext_bdtot_index)
          end do
          ext_bdtot_index=ext_bdtot_index+this%mband
        end do
      end do
      call xmpi_sum(nelect,xmpi_world,ierr)
    end if

    ! Computes extfpmd contribution to nelect using a sum
    ! of Fermi gas contributions for each point of the fftf grid.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==10) then
      ABI_MALLOC(gamma_hybrid_tf,(this%nfft,this%nspden))
      ABI_MALLOC(xcut_hybrid_tf,(this%nfft,this%nspden))
      gamma_hybrid_tf(:,:)=(fermie-this%vtrial(:,:))/tsmear
      xcut_hybrid_tf(:,:)=(this%ebcut-this%vtrial(:,:))/tsmear
      if(ANY(this%ebcut.lt.this%vtrial(:,:))) xcut_hybrid_tf(:,:)=zero

      !$OMP PARALLEL DO
      do ifft=1,this%nfft
        do ispden=1,this%nspden
          this%nelectarr(ifft,ispden)=factor*djp12(xcut_hybrid_tf(ifft,ispden),gamma_hybrid_tf(ifft,ispden))
        end do
      end do
      !$OMP END PARALLEL DO

      nelect=nelect+sum(this%nelectarr)/(this%nfft*this%nspden)
      gamma_hybrid_tf(:,:)=zero
      xcut_hybrid_tf(:,:)=zero
      ABI_FREE(gamma_hybrid_tf)
      ABI_FREE(xcut_hybrid_tf)
    end if
  end subroutine compute_nelect
  !!***

  !!****f* ABINIT/m_extfpmd/compute_e_kinetic
  !! NAME
  !!  compute_e_kinetic
  !!
  !! FUNCTION
  !!  Computes the value of the integral corresponding to the missing
  !!  kinetic energy contribution of free electrons after band cut,
  !!  with an order 3/2 incomplete Fermi-Dirac integral.
  !!
  !! INPUTS
  !!  this=extfpmd_type object concerned
  !!  fermie=chemical potential (Hartree)
  !!  tsmear=smearing width (or temperature)
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! SOURCE
  subroutine compute_e_kinetic(this,fermie,tsmear,effmass_free,gmet,kptns,nkpt,mkmem,istwfk,&
  & nspinor,nsppol,nband,wtk,total_e_kinetic,total_e_eigenvalues)
    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: nkpt,mkmem,nspinor,nsppol
    class(extfpmd_type),intent(inout) :: this
    real(dp),intent(in) :: fermie,tsmear,effmass_free
    ! Arrays
    integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol)
    real(dp),intent(in) :: gmet(3,3),kptns(3,nkpt),wtk(nkpt)
    real(dp),intent(inout) :: total_e_kinetic,total_e_eigenvalues

    ! Local variables -------------------------
    ! Scalars
    logical :: cut_warn=.false.
    integer :: ikpt,isppol,istwf_k,my_nspinor,nband_k,iband,ierr,bandstart
    integer :: ifft,ispden,ext_bdtot_index,ext_icg,ext_ikg,ext_npw_k
    real(dp) :: factor,gamma,xcut,dotr,e_eigenvalues
    real(dp) :: e_kinetic_hybrid_tf,maxocc,fn
    character(len=500) :: msg
    ! Arrays
    real(dp) :: kpoint(3)
    real(dp),allocatable :: gamma_hybrid_tf(:,:),xcut_hybrid_tf(:,:)
    real(dp),allocatable :: ext_cwavef(:,:),ext_kinpw(:)
    integer,allocatable :: ext_kg_k(:,:)

    ! *********************************************************************

    dotr=zero
    maxocc=two/(nsppol*nspinor)
    this%e_kinetic=zero
    e_eigenvalues=zero
    factor=dsqrt(two)/(PI*PI)*this%ucvol*tsmear**(2.5)
    gamma=(fermie-this%eshift)/tsmear

    ! Computes extfpmd contribution to kinetic energy integrating
    ! over accessible states from bcut to infinity with
    ! order 3/2 incomplete Fermi-Dirac integral.
    if(this%version==2.or.this%version==4) then
      xcut=extfpmd_e_fg(one*this%bcut+this%bandshift,this%ucvol)/tsmear
      if(one*this%bcut+this%bandshift.lt.zero) then
        cut_warn=.true.
        xcut=zero
      end if
      this%e_kinetic=this%e_kinetic+factor*djp32(xcut,gamma)
    end if

    ! Computes extfpmd contribution to kinetic energy integrating
    ! over energy from ebcut to infinity with order 3/2
    ! incomplete Fermi-Dirac integral.
    if(this%version==1.or.this%version==3) then
      xcut=(this%ebcut-this%eshift)/tsmear
      if(this%ebcut.lt.this%eshift) then
        cut_warn=.true.
        xcut=zero
      end if
      this%e_kinetic=this%e_kinetic+factor*djp32(xcut,gamma)
    end if
    
    ! Computes extfpmd contribution to kinetic energy summing
    ! over accessible states from bcut to infinity, with
    ! integer band numbers. Total number of bands
    ! is controlled with the input variable extfpmd_nband.
    if(this%version==5) then
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          if(proc_distrb_cycle(this%mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,this%mpi_enreg%me_kpt)) cycle
          if(one*this%bcut+this%bandshiftk(ikpt+(isppol-1)*nkpt).lt.zero) then
            cut_warn=.true.
          end if
          do iband=nband_k-this%nbdbuf+1,this%mband
            dotr=extfpmd_e_fg(one*iband+this%bandshiftk(ikpt+(isppol-1)*nkpt),this%ucvol)
            fn=fermi_dirac(dotr+this%eshift,fermie,tsmear)
            this%e_kinetic=this%e_kinetic+wtk(ikpt)*maxocc*fn*dotr
          end do
        end do
      end do
      call xmpi_sum(this%e_kinetic,xmpi_world,ierr)
    end if

    ! Computes extfpmd contribution to kinetic energy using a sum
    ! of Fermi gas contributions for each point of the fftf grid.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==10) then
      ABI_MALLOC(gamma_hybrid_tf,(this%nfft,this%nspden))
      ABI_MALLOC(xcut_hybrid_tf,(this%nfft,this%nspden))
      gamma_hybrid_tf(:,:)=(fermie-this%vtrial(:,:))/tsmear
      xcut_hybrid_tf(:,:)=(this%ebcut-this%vtrial(:,:))/tsmear
      e_kinetic_hybrid_tf=zero

      !$OMP PARALLEL DO REDUCTION (+:e_kinetic_hybrid_tf)
      do ifft=1,this%nfft
        do ispden=1,this%nspden
          e_kinetic_hybrid_tf=e_kinetic_hybrid_tf+factor*djp32(xcut_hybrid_tf(ifft,ispden),gamma_hybrid_tf(ifft,ispden))/&
          & (this%nfft*this%nspden)
        end do
      end do
      !$OMP END PARALLEL DO

      this%e_kinetic=e_kinetic_hybrid_tf
      gamma_hybrid_tf(:,:)=zero
      xcut_hybrid_tf(:,:)=zero
      ABI_FREE(gamma_hybrid_tf)
      ABI_FREE(xcut_hybrid_tf)
    end if

    ! Computes extended pw contribution to kinetic energy summing
    ! over ext pw states from nband_k to this%mband
    if(this%version==11) then
      my_nspinor=max(1,nspinor/this%mpi_enreg%nproc_spinor)
      ! Loop over spins
      ext_bdtot_index=0
      ext_icg=0
      do isppol=1,nsppol
        ext_ikg=0
        ! Loop over k points
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          istwf_k=istwfk(ikpt)
          ext_npw_k=this%npwarr(ikpt)
          ! Skip this k-point if not the proper processor
          if(proc_distrb_cycle(this%mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,this%mpi_enreg%me_kpt)) then
            ext_bdtot_index=ext_bdtot_index+this%mband
            cycle
          end if
          
          bandstart=nband_k-this%nbdbuf+1
          if(this%truecg) then
            kpoint(:)=kptns(:,ikpt)
            ABI_MALLOC(ext_kg_k,(3,ext_npw_k))
            ABI_MALLOC(ext_kinpw,(ext_npw_k))
            ABI_MALLOC(ext_cwavef,(2,ext_npw_k*my_nspinor*this%mband))
            ext_kg_k(:,1:ext_npw_k)=this%kg(:,1+ext_ikg:ext_npw_k+ext_ikg)
            ext_kinpw(:)=zero
            call mkkin(this%ecut,zero,effmass_free,gmet,ext_kg_k,ext_kinpw,kpoint,ext_npw_k,0,0)
            if(this%truecg) bandstart=1
          end if
          
          ! Compute kinetic energy of each band
          do iband=bandstart,this%mband
            if(this%truecg) then
              ext_cwavef(1:2,1:ext_npw_k*my_nspinor)= &
              & this%cg(:,1+(iband-1)*ext_npw_k*my_nspinor+ext_icg:iband*ext_npw_k*my_nspinor+ext_icg)
              call meanvalue_g(dotr,ext_kinpw,0,istwf_k,this%mpi_enreg,ext_npw_k,my_nspinor,ext_cwavef,ext_cwavef,0)
            else
              dotr=extfpmd_e_fg(one*iband+this%bandshift,this%ucvol)
            end if
            this%e_kinetic=this%e_kinetic+wtk(ikpt)*this%occ(iband+ext_bdtot_index)*dotr
            e_eigenvalues=e_eigenvalues+wtk(ikpt)*this%occ(iband+ext_bdtot_index)*this%eigen(iband+ext_bdtot_index)
          end do

          ! Increment indexes
          ext_bdtot_index=ext_bdtot_index+this%mband
          if(mkmem/=0) then
            ext_icg=ext_icg+ext_npw_k*my_nspinor*this%mband
            ext_ikg=ext_ikg+ext_npw_k
          end if

          if(this%truecg) then
            ABI_FREE(ext_cwavef)
            ABI_FREE(ext_kinpw)
            ABI_FREE(ext_kg_k)
          end if
        end do
      end do
      call xmpi_sum(this%e_kinetic,xmpi_world,ierr)
      call xmpi_sum(e_eigenvalues,xmpi_world,ierr)
      if(this%truecg) then
        total_e_eigenvalues=e_eigenvalues
        total_e_kinetic=this%e_kinetic
      end if
    end if

    ! Computes the double counting term from the eshift, and
    ! from the contributions to the kinetic energy and
    ! the number of electrons
    if(this%version==10) then
      this%edc_kinetic=this%e_kinetic+sum(this%nelectarr(:,:)*this%vtrial(:,:)/(this%nfft*this%nspden))
    else
      this%edc_kinetic=this%e_kinetic+this%nelect*this%eshift
    end if

    if(cut_warn) then
      write(msg,'(11a)')&
      & 'Extended FPMD could not properly compute the contribution to the energy.',ch10,&
      & 'This can be due to a too low number of bands in the calculation.',ch10,&
      & 'This can also happen when restarting from a previous calculation.',ch10,&
      & 'Poor prediction of the electron density based on forces may results in this error.',ch10,&
      & 'Action: slightly increase nband if the electron density is supposed to be converged.',ch10,&
      & 'Otherwise: wait for the density to be converged.'
      ABI_WARNING(msg)
    end if
  end subroutine compute_e_kinetic
  !!***

  !!****f* ABINIT/m_extfpmd/compute_entropy
  !! NAME
  !!  compute_entropy
  !!
  !! FUNCTION
  !!  Computes the value of the integral corresponding to the missing
  !!  entropy contribution of free electrons after band cut using
  !!  incomplete Fermi-Dirac integrals.
  !!
  !! INPUTS
  !!  this=extfpmd_type object concerned
  !!  fermie=chemical potential (Hartree)
  !!  tsmear=smearing width (or temperature)
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! SOURCE
  subroutine compute_entropy(this,fermie,tsmear,nkpt,nsppol,nspinor,wtk,nband,mband,occ)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: nkpt,nsppol,nspinor,mband
    real(dp),intent(in) :: fermie,tsmear
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)
    real(dp),intent(in) :: occ(mband*nkpt*nsppol)

    ! Local variables -------------------------
    ! Scalars
    integer :: ii,ifft,ispden,ext_bdtot_index,isppol,ikpt,iband,nband_k,ierr
    real(dp) :: ix,step,factor,fn,gamma,maxocc
    ! Arrays
    real(dp),dimension(:),allocatable :: valuesent
    real(dp),dimension(:,:),allocatable :: gamma_hybrid_tf
    real(dp),dimension(:,:),allocatable :: step_hybrid_tf

    ! *********************************************************************
    maxocc=two/(nsppol*nspinor)
    this%entropy=zero
    factor=dsqrt(two)/(PI*PI)*this%ucvol*tsmear**(2.5)
    gamma=(fermie-this%eshift)/tsmear
    ABI_MALLOC(valuesent,(this%bcut+1))

    ! Computes extfpmd contribution to the entropy integrating
    ! over accessible states with Fermi-Dirac complete integrals and
    ! substracting 0 to bcut contribution with numeric integration.
    if(this%version==2.or.this%version==4) then
      step=(dble(this%bcut)+this%bandshift)/(this%bcut)
      !$OMP PARALLEL DO PRIVATE(fn,ix) SHARED(valuesent)
      do ii=1,this%bcut+1
        ix=(dble(ii)-one)*step
        fn=fermi_dirac(extfpmd_e_fg(ix,this%ucvol)+this%eshift,fermie,tsmear)
        if(fn>tol16.and.(one-fn)>tol16) then
          valuesent(ii)=-maxocc*(fn*log(fn)+(one-fn)*log(one-fn))
        else
          valuesent(ii)=zero
        end if
      end do
      !$OMP END PARALLEL DO

      ! We need at least 6 elements in valuesent to call simpson function.
      if(size(valuesent)>=6) then
        this%entropy=5./3.*factor*dip32(gamma)/tsmear-&
        & gamma*factor*dip12(gamma)/tsmear-&
        simpson(step,valuesent)
      end if
    end if

    ! Computes extfpmd contribution to the entropy integrating
    ! over energy with Fermi-Dirac complete integrals and
    ! substracting 0 to bcut contribution with numeric integration.
    if(this%version==1.or.this%version==3) then
      step=(this%ebcut-this%eshift)/(this%bcut)
      !$OMP PARALLEL DO PRIVATE(fn,ix) SHARED(valuesent)
      do ii=1,this%bcut+1
        ix=this%eshift+(dble(ii)-one)*step
        fn=fermi_dirac(ix,fermie,tsmear)
        if(fn>tol16.and.(one-fn)>tol16) then
          valuesent(ii)=-(fn*log(fn)+(one-fn)*log(one-fn))*&
          & extfpmd_dos(ix,this%eshift,this%ucvol)
        else
          valuesent(ii)=zero
        end if
      end do
      !$OMP END PARALLEL DO

      ! We need at least 6 elements in valuesent to call simpson function.
      if(size(valuesent)>=6) then
        this%entropy=5./3.*factor*dip32(gamma)/tsmear-&
        & gamma*factor*dip12(gamma)/tsmear-simpson(step,valuesent)
      end if
    end if

    ! Computes extfpmd contribution to the entropy summing
    ! over accessible states from bcut to infinity, with
    ! integer band numbers. Total number of bands
    ! is controlled with the input variable extfpmd_nband.
    if(this%version==5) then
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          if(proc_distrb_cycle(this%mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,this%mpi_enreg%me_kpt)) cycle
          do iband=nband_k-this%nbdbuf+1,this%mband
            fn=fermi_dirac(extfpmd_e_fg(one*iband+this%bandshiftk(ikpt+(isppol-1)*nkpt),this%ucvol)+this%eshift,fermie,tsmear)
            this%entropy=this%entropy-wtk(ikpt)*maxocc*(fn*log(fn)+(one-fn)*log(one-fn))/nsppol
          end do
        end do
      end do
      call xmpi_sum(this%entropy,xmpi_world,ierr)
    end if

    ! Computes extended pw contribution to nelect summing
    ! over ext pw states from nband_k to this%mband.
    ! In this case, entropy is computed just for information as
    ! ext pw entropy is already included while computing occupations.
    if(this%version==11) then
      ext_bdtot_index=0
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          ! Skip this k-point if not the proper processor
          if(proc_distrb_cycle(this%mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,this%mpi_enreg%me_kpt)) then
            ext_bdtot_index=ext_bdtot_index+this%mband
            cycle
          end if

          do iband=nband_k-this%nbdbuf+1,this%mband
            fn=this%occ(iband+ext_bdtot_index)/maxocc
            this%entropy=this%entropy-wtk(ikpt)*maxocc*(fn*log(fn)+(one-fn)*log(one-fn))/nsppol
          end do

          ext_bdtot_index=ext_bdtot_index+this%mband
        end do
      end do
      call xmpi_sum(this%entropy,xmpi_world,ierr)
    end if

    ! Computes extfpmd contribution to the entropy using a sum
    ! of Fermi gas contributions for each point of the fftf grid,
    ! as we do for version=1 and version=2.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==10) then
      ABI_MALLOC(gamma_hybrid_tf,(this%nfft,this%nspden))
      ABI_MALLOC(step_hybrid_tf,(this%nfft,this%nspden))
      gamma_hybrid_tf(:,:)=(fermie-this%vtrial(:,:))/tsmear
      step_hybrid_tf(:,:)=(this%ebcut-this%vtrial(:,:))/(this%bcut)
      this%entropy=zero

      do ifft=1,this%nfft
        do ispden=1,this%nspden
          !$OMP PARALLEL DO PRIVATE(ix,fn) SHARED(valuesent)
          do ii=1,this%bcut+1
            ix=this%vtrial(ifft,ispden)+(dble(ii)-one)*step_hybrid_tf(ifft,ispden)
            fn=fermi_dirac(ix,fermie,tsmear)
            if(fn>tol16.and.(one-fn)>tol16) then
              valuesent(ii)=-(fn*log(fn)+(one-fn)*log(one-fn))*&
              & extfpmd_dos(ix,this%vtrial(ifft,ispden),this%ucvol)
            else
              valuesent(ii)=zero
            end if
          end do
          !$OMP END PARALLEL DO

          ! We need at least 6 elements in valuesent to call simpson function.
          if(size(valuesent)>=6) then
            this%entropy=this%entropy+(5./3.*factor*dip32(gamma_hybrid_tf(ifft,ispden))/tsmear-&
            & gamma_hybrid_tf(ifft,ispden)*factor*dip12(gamma_hybrid_tf(ifft,ispden))/tsmear-&
            & simpson(step_hybrid_tf(ifft,ispden),valuesent))/(this%nfft*this%nspden)
          end if
        end do
      end do

      gamma_hybrid_tf(:,:)=zero
      step_hybrid_tf(:,:)=zero
      ABI_FREE(step_hybrid_tf)
      ABI_FREE(gamma_hybrid_tf)
    end if
    ABI_FREE(valuesent)
  end subroutine compute_entropy
  !!***

  !----------------------------------------------------------------------

  !!****f* ABINIT/m_extfpmd/extfpmd_dos
  !! NAME
  !!  extfpmd_dos
  !!
  !! FUNCTION
  !!  Returns the free particle density of states for a given energy.
  !!
  !! INPUTS
  !!  energy=get the value of the free particle density of states at this energy
  !!  eshift=energy shift factor
  !!  ucvol=unit cell volume (bohr^3)
  !!
  !! OUTPUT
  !!  extfpmd_dos=value of free particle density of states at given energy
  !!
  !! SOURCE
  function extfpmd_dos(energy,eshift,ucvol)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: energy,eshift,ucvol
    real(dp) :: extfpmd_dos

    ! *********************************************************************

    extfpmd_dos=dsqrt(two)*ucvol*dsqrt(energy-eshift)/(PI*PI)
  end function extfpmd_dos
  !!***

  !!****f* ABINIT/m_extfpmd/extfpmd_e_fg
  !! NAME
  !!  extfpmd_e_fg
  !!
  !! FUNCTION
  !!  Returns the energy of the Fermi gas for a given number of
  !!  accessible states.
  !!
  !! INPUTS
  !!  iband=number of accessible states
  !!  ucvol=unit cell volume (bohr^3)
  !!
  !! OUTPUT
  !!  extfpmd_e_fg=energy of homogeneous electron gas for a given number of accessible states
  !!
  !! SOURCE
  function extfpmd_e_fg(iband,ucvol)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: iband,ucvol
    real(dp) :: extfpmd_e_fg

    ! *********************************************************************

    extfpmd_e_fg=half*(iband*six*PI*PI/ucvol)**(two/three)
  end function extfpmd_e_fg

  !!***
  !!****f* ABINIT/m_extfpmd/extfpmd_i_fg
  !! NAME
  !!  extfpmd_i_fg
  !!
  !! FUNCTION
  !!  Returns the number of doubly occupied orbitals of the
  !!  for a Fermi Gas for a given kinetic energy.
  !!
  !! INPUTS
  !!  ekin=kinetic energy
  !!  ucvol=unit cell volume (bohr^3)
  !!
  !! OUTPUT
  !!  extfpmd_i_fg=number of doubly occupied states of the Fermi gas
  !!
  !! SOURCE
  function extfpmd_i_fg(ekin,ucvol)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: ekin,ucvol
    real(dp) :: extfpmd_i_fg

    ! *********************************************************************

    extfpmd_i_fg=(two*ekin)**(three/two)*ucvol/(six*PI*PI)
  end function extfpmd_i_fg
  !!***
  
  !!***
  !!****f* ABINIT/m_extfpmd/extfpmd_chkinp
  !! NAME
  !!  extfpmd_chkinp
  !!
  !! FUNCTION
  !!  Returns true if extfpmd input variables are valid and false otherwise.
  !!
  !! INPUTS
  !!  dtset=<type datafiles_type>contains all input variables.
  !!
  !! OUTPUT
  !!  extfpmd_chkinp=allocate extfpmd object or not
  !!
  !! SOURCE
  function extfpmd_chkinp(dtset)
    ! Arguments -------------------------------
    ! Scalars
    class(dataset_type),intent(in) :: dtset
    logical :: extfpmd_chkinp
    ! Local variables -------------------------
    ! Scalars
    character(len=500) :: msg
    
    ! *********************************************************************
    
    extfpmd_chkinp=.false.
    if(dtset%extfpmd_nbcut>dtset%mband) then
      write(msg,'(3a,i0,a,i0,3a)') "Not enough bands to activate ExtFPMD routines.",ch10,&
      & "extfpmd_nbcut = ",dtset%extfpmd_nbcut," should be less than or equal to nband = ",dtset%mband,".",ch10,&
      & "Action: Increase nband or decrease extfpmd_nbcut."
      ABI_ERROR(msg)
    else if(dtset%extfpmd_nbdbuf+dtset%extfpmd_nbcut>dtset%mband) then
      write(msg,'(a,i0,a,i0,a,i0,2a,i0,a)') "(extfpmd_nbdbuf = ",dtset%extfpmd_nbdbuf," + extfpmd_nbcut = ",&
      & dtset%extfpmd_nbcut,") = ",dtset%extfpmd_nbdbuf+dtset%extfpmd_nbcut,ch10,&
      & "should be less than or equal to nband = ",dtset%mband,"."
      ABI_ERROR(msg)
    else if(dtset%extfpmd_nbdbuf>dtset%mband) then
      write(msg,'(a,i0,a,i0,a)') "extfpmd_nbdbuf = ",dtset%extfpmd_nbdbuf,&
      & " should be less than or equal to nband = ",dtset%mband,"."
      ABI_ERROR(msg)
    else if((dtset%useextfpmd==5.or.dtset%useextfpmd==11).and.(dtset%extfpmd_nband<=dtset%mband)) then
      write(msg,'(3a,i0,a,i0,3a)') "Not enough bands to activate ExtFPMD routines.",ch10,&
      & "extfpmd_nband = ",dtset%extfpmd_nband," should be strictly greater than nband = ",dtset%mband,".",ch10,&
      & "Action: Increase extfpmd_nband or decrease nband."
      ABI_ERROR(msg)
    else
      extfpmd_chkinp=.true.
    end if
  end function extfpmd_chkinp
  !!***

end module m_extfpmd
!!***
