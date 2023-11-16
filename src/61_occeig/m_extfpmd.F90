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
  use m_energies,       only : energies_type
  use m_gsphere,        only : getkpgnorm
  use m_kg,             only : kpgio, getmpw, mkkin
  use m_mpinfo,         only : ptabs_fourdp,proc_distrb_cycle
  use m_numeric_tools,  only : simpson,simpson_int
  use m_spacepar,       only : meanvalue_g

  implicit none
  public :: extfpmd_dos,extfpmd_e_fg
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
    integer :: bcut,nbcut,nbdbuf,nfft,nspden,version
    real(dp) :: e_bcut,edc_kinetic,e_kinetic,entropy
    real(dp) :: nelect,shiftfactor,ucvol
    real(dp),allocatable :: vtrial(:,:)
    real(dp),allocatable :: nelectarr(:,:)
    !! Scalars and arrays for numerical extended PW method
    integer :: mpw,mcg,mband
    real(dp) :: ecut,ecut_eff
    integer,allocatable :: kg(:,:),npwarr(:),npwtot(:),nband(:)
    real(dp),allocatable :: cg(:,:),eigen(:),occ(:),doccde(:)
  contains
    procedure :: compute_e_kinetic
    procedure :: compute_entropy
    procedure :: compute_nelect
    procedure :: generate_extpw
    procedure :: update_ks_occ
    procedure :: compute_shiftfactor
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
  !!  ecut=extended plane waves basis set energy cut off
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! SOURCE
  subroutine init(this,mband,nbcut,nbdbuf,nfft,nspden,nsppol,nkpt,rprimd,version,&
  & exchn2n3d,istwfk,kptns,mpi_enreg,mkmem,dilatmx,extfpmd_mband)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: mband,nbcut,nbdbuf,nfft,nspden,nsppol,nkpt,version
    integer,intent(in) :: exchn2n3d,mkmem,extfpmd_mband
    real(dp),intent(in) :: dilatmx
    type(MPI_type),intent(inout) :: mpi_enreg
    ! Arrays
    integer,intent(in) :: istwfk(nkpt)
    real(dp),intent(in) :: rprimd(3,3)
    real(dp),intent(in) :: kptns(3,nkpt)

    ! Local variables -------------------------
    ! Arrays
    real(dp) :: gprimd(3,3),rmet(3,3),gmet(3,3)

    ! *********************************************************************

    this%bcut=mband-nbdbuf
    this%nbcut=nbcut
    this%nbdbuf=nbdbuf
    this%version=version
    ABI_MALLOC(this%vtrial,(nfft,nspden))
    this%vtrial(:,:)=zero
    this%nfft=nfft
    this%nspden=nspden
    this%e_bcut=zero
    this%edc_kinetic=zero
    this%e_kinetic=zero
    this%entropy=zero
    this%nelect=zero
    ABI_MALLOC(this%nelectarr,(nfft,nspden))
    this%nelectarr(:,:)=zero
    this%shiftfactor=zero
    this%mband=extfpmd_mband
    call metric(gmet,gprimd,-1,rmet,rprimd,this%ucvol)
    ! Add sqrt(...) to extended cutoff energy to make
    ! sure extended pw basis set is large enough.
    this%ecut=extfpmd_e_fg(one*extfpmd_mband,this%ucvol)+&
    & sqrt((2*PI*gprimd(1,1))**2+(2*PI*gprimd(2,1))**2+(2*PI*gprimd(3,1))**2)
    this%ecut_eff=this%ecut*dilatmx**2
    ABI_MALLOC(this%nband,(nkpt*nsppol))
    this%nband(:)=extfpmd_mband
    ABI_MALLOC(this%npwarr,(nkpt))
    this%npwarr(:)=0
    ABI_MALLOC(this%npwtot,(nkpt))
    this%npwtot(:)=0
    call getmpw(this%ecut_eff,exchn2n3d,gmet,istwfk,kptns,mpi_enreg,this%mpw,nkpt)
    ABI_MALLOC(this%kg,(3,this%mpw*mkmem))
    this%kg(:,:)=0
    this%mcg=0
    ABI_MALLOC(this%eigen,(extfpmd_mband*nkpt*nsppol))
    this%eigen(:)=zero
    ABI_MALLOC(this%occ,(extfpmd_mband*nkpt*nsppol))
    this%occ(:)=zero
    ABI_MALLOC(this%doccde,(extfpmd_mband*nkpt*nsppol))
    this%doccde(:)=zero
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
    if (allocated(this%cg)) then
      this%cg(:,:)=zero
      ABI_FREE(this%cg)
    end if
    this%doccde(:)=zero
    ABI_FREE(this%doccde)
    this%eigen(:)=zero
    ABI_FREE(this%eigen)
    this%occ(:)=zero
    ABI_FREE(this%occ)
    this%kg(:,:)=0
    ABI_FREE(this%kg)
    this%npwtot(:)=0
    ABI_FREE(this%npwtot)
    this%npwarr(:)=0
    ABI_FREE(this%npwarr)
    this%nband(:)=0
    ABI_FREE(this%nband)
    this%vtrial(:,:)=zero
    ABI_FREE(this%vtrial)
    this%nelectarr(:,:)=zero
    ABI_FREE(this%nelectarr)
    this%mpw=0
    this%nfft=0
    this%nspden=0
    this%bcut=0
    this%mband=0
    this%nbcut=0
    this%nbdbuf=0
    this%version=1
    this%e_bcut=zero
    this%edc_kinetic=zero
    this%e_kinetic=zero
    this%entropy=zero
    this%nelect=zero
    this%shiftfactor=zero
    this%ecut=zero
    this%ecut_eff=zero
    this%ucvol=zero
  end subroutine destroy
  !!***

  !!****f* ABINIT/m_extfpmd/compute_shiftfactor
  !! NAME
  !!  compute_shiftfactor
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
  !!  nkpt=number of k points
  !!  nsppol=1 for unpolarized, 2 for spin-polarized
  !!  wtk(nkpt)=k point weights
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! SOURCE
  subroutine compute_shiftfactor(this,eigen,eknk,mband,me,nband,nkpt,nsppol,wtk)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: mband,me,nkpt,nsppol
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
    real(dp),intent(in) :: eknk(mband*nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)

    ! Local variables -------------------------
    ! Scalars
    integer :: band_index,ii,ikpt,isppol,nband_k
    real(dp) :: abs_err,rel_err
    character(len=500) :: msg

    ! *********************************************************************

    ! Computes U_0 from the sum of local
    ! potentials (vtrial), averaging over all space.
    ! Simplest and most precise way to evaluate U_0.
    if(this%version==1.or.this%version==4.or.this%version==5.or.this%version==10) then
      this%shiftfactor=sum(this%vtrial)/(this%nfft*this%nspden)

      ! Computes the relative error of the model vs last eigenvalues
      if(me==0.and.this%version==4) then
        band_index=0
        rel_err=zero
        abs_err=zero
        do isppol=1,nsppol
          do ikpt=1,nkpt
            nband_k=nband(ikpt+(isppol-1)*nkpt)
            rel_err=rel_err+wtk(ikpt)*abs((eigen(band_index+nband_k-this%nbdbuf)-&
            & extfpmd_e_fg(dble(nband_k-this%nbdbuf),this%ucvol)-this%shiftfactor)/&
            & eigen(band_index+nband_k-this%nbdbuf))
            abs_err=abs_err+wtk(ikpt)*abs(eigen(band_index+nband_k-this%nbdbuf)-&
            & extfpmd_e_fg(dble(nband_k-this%nbdbuf),this%ucvol)-this%shiftfactor)
            band_index=band_index+nband_k
          end do
        end do
        if(rel_err.gt.tol1) then
          write(msg,'(a,es8.2,3a,es8.2,a,es8.2,7a)')&
          & 'Relative difference between eigenvalues and Fermi gas energy (',rel_err,')',ch10,&
          & 'is over ',tol1,' at band cut. Absolute difference is ',abs_err,' Ha.',ch10,&
          & 'Execution will continue as the code will still add contributions in the right',ch10,&
          & 'direction, but you should likely increase nband to make sure electrons of last',ch10,&
          & 'band can be considered as free fermions in a constant potential.'
          ABI_WARNING(msg)
        end if
      end if
    end if

    ! Computes U_0^{HEG} from the difference between
    ! eigenvalues and Fermi gas energies, averaged
    ! over lasts nbcut bands.
    if(this%version==2) then
      this%shiftfactor=zero
      band_index=0
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          do ii=nband_k-this%nbdbuf-this%nbcut+1,nband_k-this%nbdbuf
            this%shiftfactor=this%shiftfactor+&
            & wtk(ikpt)*(eigen(band_index+ii)-extfpmd_e_fg(dble(ii),this%ucvol))
          end do
          band_index=band_index+nband_k
        end do
      end do
      this%shiftfactor=this%shiftfactor/this%nbcut
    end if

    ! Computes U_0^K from the difference between
    ! eigenvalues and kinetic energies, averaged
    ! over lasts nbcut bands.
    if(this%version==3) then
      this%shiftfactor=zero
      band_index=0
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          do ii=nband_k-this%nbdbuf-this%nbcut+1,nband_k-this%nbdbuf
            this%shiftfactor=this%shiftfactor+&
            & wtk(ikpt)*(eigen(band_index+ii)-eknk(band_index+ii))
          end do
          band_index=band_index+nband_k
        end do
      end do
      this%shiftfactor=this%shiftfactor/this%nbcut
    end if

    ! Get extended FPMD band energy cutoff for version 1, 3 and 10.
    if(this%version==1.or.this%version==3.or.this%version==10) then
      this%e_bcut=0
      band_index=0
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          this%e_bcut=this%e_bcut+wtk(ikpt)*eigen(band_index+nband_k-this%nbdbuf)
          band_index=band_index+nband_k
        end do
      end do
    end if
  end subroutine compute_shiftfactor
  !!***

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
    & mode_paral,mpi_enreg,nsppol,dilatmx,nspinor,cg,mcg,npwarr,kg,mpw,eigen,mband)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: exchn2n3d,mkmem,nkpt,nsppol,nspinor,mcg,mpw,mband
    real(dp),intent(in) :: dilatmx,effmass_free
    character(len=4),intent(in) :: mode_paral
    type(MPI_type),intent(inout) :: mpi_enreg
    ! Arrays
    integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt),kg(3,mpw*mkmem)
    real(dp),intent(in) :: gmet(3,3),kptns(3,nkpt)
    real(dp),intent(in) :: cg(2,mcg),eigen(mband*nkpt*nsppol)
    
    ! Local variables -------------------------
    ! Scalars
    integer :: isppol,ikpt,iband,icg,ext_icg,ipw,ext_ipw,ikg,ext_ikg
    integer :: nband_k,ext_nband_k,istwf_k,npw_k,ext_npw_k,nblockbd
    integer :: my_nspinor,my_ikpt,blocksize,bdtot_index,ext_bdtot_index
    integer :: count_below,count_above,index_below,index_above
    real(dp) :: ecut_eff,ekin_max,fg_kin,closest_below,closest_above,prop_below,prop_above
    real(dp) :: norm
    ! Arrays
    real(dp) :: kpoint(3)
    integer,allocatable :: kg_k(:,:),ext_kg_k(:,:),indices_below(:),indices_above(:)
    real(dp),allocatable :: ext_kinpw(:),eig_k(:)

    ! *********************************************************************

    ! Generating ext pw basis set (this%kg,this%npwarr,this%npwtot)
    if(.not.allocated(this%cg)) then
      write(0,*) 'DEBUG: Generating extended plane wave basis set...'
      call kpgio(this%ecut_eff,exchn2n3d,gmet,istwfk,this%kg, &
      & kptns,mkmem,nband,nkpt,'PERS',mpi_enreg,&
      & this%mpw,this%npwarr,this%npwtot,nsppol)
      
      ! Get number of ext pw coefficients
      my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)
      this%mcg=this%mpw*my_nspinor*this%mband*mkmem*nsppol
      ABI_MALLOC(this%cg,(2,this%mcg))
      this%cg(:,:)=zero
      
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
            cycle
          end if

          kpoint(:)=kptns(:,ikpt)

          ABI_MALLOC(eig_k,(nband_k))
          ABI_MALLOC(kg_k,(3,npw_k))
          ABI_MALLOC(ext_kg_k,(3,ext_npw_k))
          ABI_MALLOC(ext_kinpw,(ext_npw_k))
          eig_k(:)=eigen(1+bdtot_index:nband_k+bdtot_index)
          if (minval(eig_k)>1.d100) eig_k=zero
          kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
          ext_kg_k(:,1:ext_npw_k)=this%kg(:,1+ext_ikg:ext_npw_k+ext_ikg)

          ! Compute kinetic energy of extended pw.
          ! /!\ We don't keep this in memory because it can be heavy.
          ext_kinpw(:)=zero
          call mkkin(this%ecut,zero,effmass_free,gmet,ext_kg_k,ext_kinpw,kpoint,ext_npw_k,0,0)

          ! Copy non extended plane waves coefficients and eigenvalues here
          do iband=1,nband_k
            do ipw=1,npw_k*my_nspinor
              do ext_ipw=1,ext_npw_k*my_nspinor
                if (all(kg_k(:,ipw)==ext_kg_k(:,ext_ipw))) then
                  this%cg(:,ext_ipw+(iband-1)*ext_npw_k*my_nspinor+ext_icg)=cg(:,ipw+(iband-1)*npw_k*my_nspinor+icg)
                  this%eigen(iband+ext_bdtot_index)=eig_k(iband)
                end if
              end do
            end do
          end do

          ! Set extended plane waves coefficients here
          do iband=nband_k+1,this%mband
            ! Get fermi gas energy and set eigenvalue
            fg_kin=extfpmd_e_fg(one*iband,this%ucvol)
            this%eigen(iband+ext_bdtot_index)=fg_kin+this%shiftfactor
            ! Find closest values in ext_kinpw

            closest_below=zero
            closest_above=this%ecut
            count_below=0
            count_above=0

            do ext_ipw=1,ext_npw_k*my_nspinor
              if (ext_kinpw(ext_ipw)<=fg_kin.and.ext_kinpw(ext_ipw)>=closest_below) then
                if (ext_kinpw(ext_ipw)==closest_below) then
                  count_below=count_below+1
                else
                  count_below=1
                  closest_below=ext_kinpw(ext_ipw)
                end if
              else if (ext_kinpw(ext_ipw)>=fg_kin.and.ext_kinpw(ext_ipw)<=closest_above) then
                if (ext_kinpw(ext_ipw)==closest_above) then
                  count_above=count_above+1
                else
                  count_above=1
                  closest_above=ext_kinpw(ext_ipw)
                end if
              end if
            end do

            prop_below=(fg_kin-closest_above)/(count_below*(closest_below-closest_above))
            prop_above=(1-prop_below*count_below)/count_above

            ! Do a second loop to set cg coefficients
            do ext_ipw=1,ext_npw_k*my_nspinor
              if (ext_kinpw(ext_ipw)==closest_below) then
                this%cg(1,ext_ipw+(iband-1)*ext_npw_k*my_nspinor+ext_icg)=sqrt(prop_below)
              end if
              if (ext_kinpw(ext_ipw)==closest_above) then
                this%cg(1,ext_ipw+(iband-1)*ext_npw_k*my_nspinor+ext_icg)=sqrt(prop_above)
              end if
            end do
          end do

          ! write(0,*) "DEBUG: Checking extended plane waves coefficients normalization"
          ! do iband=1,this%mband
          !   fg_kin=extfpmd_e_fg(one*iband,this%ucvol)
          !   norm=zero
          !   do ext_ipw=1,ext_npw_k*my_nspinor
          !     norm=norm+(this%cg(1,ext_ipw+(iband-1)*ext_npw_k*my_nspinor+ext_icg)**2+this%cg(2,ext_ipw+(iband-1)*ext_npw_k*my_nspinor+ext_icg)**2)
          !   end do
          !   write(0,*) ikpt, iband, norm, this%eigen(iband+ext_bdtot_index), fg_kin
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
          ABI_FREE(eig_k)
        end do
      end do
      write(0,*) "DEBUG: Exiting routine..."
    else

      my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)

      ! Loop over spins
      bdtot_index=0
      ext_bdtot_index=0
      icg=0
      ext_icg=0
      do isppol=1,nsppol
        ikg=0
        ext_ikg=0
        ! Loop over k points
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          npw_k=npwarr(ikpt)
          ext_npw_k=this%npwarr(ikpt)
          ! Skip this k-point if not the proper processor
          if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,mpi_enreg%me_kpt)) then
            bdtot_index=bdtot_index+nband_k
            cycle
          end if

          ABI_MALLOC(eig_k,(nband_k))
          ABI_MALLOC(kg_k,(3,npw_k))
          ABI_MALLOC(ext_kg_k,(3,ext_npw_k))
          eig_k(:)=eigen(1+bdtot_index:nband_k+bdtot_index)
          if (minval(eig_k)>1.d100) eig_k=zero
          kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
          ext_kg_k(:,1:ext_npw_k)=this%kg(:,1+ext_ikg:ext_npw_k+ext_ikg)

          ! Copy updated non extended plane waves coefficients and eigenvalues here
          do iband=1,nband_k
            do ipw=1,npw_k*my_nspinor
              do ext_ipw=1,ext_npw_k*my_nspinor
                if (all(kg_k(:,ipw)==ext_kg_k(:,ext_ipw))) then
                  this%cg(:,ext_ipw+(iband-1)*ext_npw_k*my_nspinor+ext_icg)=cg(:,ipw+(iband-1)*npw_k*my_nspinor+icg)
                  this%eigen(iband+ext_bdtot_index)=eig_k(iband)
                end if
              end do
            end do
          end do

          ! Increment indexes
          bdtot_index=bdtot_index+nband_k
          ext_bdtot_index=ext_bdtot_index+this%mband
          if(mkmem/=0) then
            icg=icg+npw_k*my_nspinor*nband_k
            ext_icg=ext_icg+ext_npw_k*my_nspinor*this%mband
            ikg=ikg+npw_k
            ext_ikg=ext_ikg+ext_npw_k
          end if

          ABI_FREE(ext_kg_k)
          ABI_FREE(kg_k)
          ABI_FREE(eig_k)
        end do
      end do
    end if
  end subroutine generate_extpw

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
  subroutine compute_nelect(this,fermie,nband,nelect,nkpt,nsppol,tsmear,wtk)
    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: nkpt,nsppol
    real(dp),intent(in) :: fermie,tsmear
    real(dp),intent(inout) :: nelect
    class(extfpmd_type),intent(inout) :: this
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)

    ! Local variables -------------------------
    ! Scalars
    integer :: ifft,ispden,isppol,ikpt,iband,index_tot,nband_k
    real(dp) :: factor,gamma,xcut
    ! Arrays
    real(dp),allocatable :: gamma_hybrid_tf(:,:)
    real(dp),allocatable :: xcut_hybrid_tf(:,:)

    ! *********************************************************************

    factor=sqrt(2.)/(PI*PI)*this%ucvol*tsmear**(1.5)

    ! Computes extfpmd contribution to nelect integrating
    ! over accessible states from bcut to infinity with
    ! order 1/2 incomplete Fermi-Dirac integral.
    if(this%version==4.or.this%version==2) then
      gamma=(fermie-this%shiftfactor)/tsmear
      xcut=extfpmd_e_fg(dble(this%bcut),this%ucvol)/tsmear
      nelect=nelect+factor*djp12(xcut,gamma)
    end if

    ! Computes extfpmd contribution to nelect integrating
    ! over energy from e_bcut to infinity with order 1/2
    ! incomplete Fermi-Dirac integral.
    if(this%version==1.or.this%version==3) then
      gamma=(fermie-this%shiftfactor)/tsmear
      xcut=(this%e_bcut-this%shiftfactor)/tsmear
      if(this%e_bcut.lt.this%shiftfactor) xcut=zero
      nelect=nelect+factor*djp12(xcut,gamma)
    end if

    ! Computes extended pw contribution to nelect summing
    ! over ext pw states from nband_k to this%mband
    if(this%version==5) then
      index_tot=0
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          do iband=nband_k+1,this%mband
            nelect=nelect+wtk(ikpt)*this%occ(iband+index_tot)
          end do
          index_tot=index_tot+this%mband
        end do
      end do
    end if

    ! Computes extfpmd contribution to nelect using a sum
    ! of Fermi gas contributions for each point of the fftf grid.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==10) then
      ABI_MALLOC(gamma_hybrid_tf,(this%nfft,this%nspden))
      ABI_MALLOC(xcut_hybrid_tf,(this%nfft,this%nspden))
      gamma_hybrid_tf(:,:)=(fermie-this%vtrial(:,:))/tsmear
      xcut_hybrid_tf(:,:)=(this%e_bcut-this%vtrial(:,:))/tsmear
      if(ANY(this%e_bcut.lt.this%vtrial(:,:))) xcut_hybrid_tf(:,:)=zero

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

  !!****f* ABINIT/m_extfpmd/update_ks_occ
  !! NAME
  !!  update_ks_occ
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
  subroutine update_ks_occ(this,occ,doccde,mband,nkpt,nsppol,mpi_enreg,nband)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: mband,nkpt,nsppol
    type(MPI_type),intent(inout) :: mpi_enreg
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(out) :: doccde(mband*nkpt*nsppol)
    real(dp),intent(inout) :: occ(mband*nkpt*nsppol)

    ! Local variables -------------------------
    ! Scalars
    integer :: bdtot_index,ext_bdtot_index,isppol,ikpt,iband,nband_k
    ! Arrays

    ! *********************************************************************
    
    bdtot_index=0
    ext_bdtot_index=0

    do isppol=1,nsppol
      do ikpt=1,nkpt
        nband_k=nband(ikpt+(isppol-1)*nkpt)
        ! Skip this k-point if not the proper processor
        if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,mpi_enreg%me_kpt)) then
          bdtot_index=bdtot_index+nband_k
          cycle
        end if

        do iband=1,nband_k
          occ(iband+bdtot_index)=this%occ(iband+ext_bdtot_index)
          doccde(iband+bdtot_index)=this%doccde(iband+ext_bdtot_index)
        end do

        ! Increment indexes
        bdtot_index=bdtot_index+nband_k
        ext_bdtot_index=ext_bdtot_index+this%mband
      end do
    end do

  end subroutine update_ks_occ

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
  subroutine compute_e_kinetic(this,fermie,tsmear,mpi_enreg,effmass_free,gmet,kptns,nkpt,mkmem,istwfk,&
  & nspinor,nsppol,nband,wtk)
    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: nkpt,mkmem,nspinor,nsppol
    class(extfpmd_type),intent(inout) :: this
    real(dp),intent(in) :: fermie,tsmear,effmass_free
    type(MPI_type),intent(inout) :: mpi_enreg
    ! Arrays
    integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol)
    real(dp),intent(in) :: gmet(3,3),kptns(3,nkpt),wtk(nkpt)

    ! Local variables -------------------------
    ! Scalars
    integer :: ikpt,isppol,istwf_k,my_nspinor,nband_k,iband
    integer :: ifft,ispden,ext_bdtot_index,ext_icg,ext_ikg,ext_npw_k
    real(dp) :: factor,gamma,xcut,dotr
    real(dp) :: e_kinetic_hybrid_tf
    character(len=500) :: msg
    ! Arrays
    real(dp) :: kpoint(3)
    real(dp),allocatable :: gamma_hybrid_tf(:,:),xcut_hybrid_tf(:,:)
    real(dp),allocatable :: cwavef(:,:),ext_kinpw(:)
    integer,allocatable :: ext_kg_k(:,:)

    ! *********************************************************************

    this%e_kinetic=zero
    factor=sqrt(2.)/(PI*PI)*this%ucvol*tsmear**(2.5)

    ! Computes extfpmd contribution to kinetic energy integrating
    ! over accessible states from bcut to infinity with
    ! order 3/2 incomplete Fermi-Dirac integral.
    if(this%version==4.or.this%version==2) then
      gamma=(fermie-this%shiftfactor)/tsmear
      xcut=extfpmd_e_fg(dble(this%bcut),this%ucvol)/tsmear
      this%e_kinetic=this%e_kinetic+factor*djp32(xcut,gamma)
    end if

    ! Computes extfpmd contribution to kinetic energy integrating
    ! over energy from e_bcut to infinity with order 3/2
    ! incomplete Fermi-Dirac integral.
    if(this%version==1.or.this%version==3) then
      gamma=(fermie-this%shiftfactor)/tsmear
      xcut=(this%e_bcut-this%shiftfactor)/tsmear
      this%e_kinetic=this%e_kinetic+factor*djp32(xcut,gamma)
    end if

    ! Computes extended pw contribution to kinetic energy summing
    ! over ext pw states from nband_k to this%mband
    if(this%version==5) then
      my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)
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

          kpoint(:)=kptns(:,ikpt)

          ABI_MALLOC(ext_kg_k,(3,ext_npw_k))
          ABI_MALLOC(ext_kinpw,(ext_npw_k))
          ABI_MALLOC(cwavef,(2,ext_npw_k*my_nspinor*this%mband))
          ext_kg_k(:,1:ext_npw_k)=this%kg(:,1+ext_ikg:ext_npw_k+ext_ikg)
          ext_kinpw(:)=zero
          call mkkin(this%ecut,zero,effmass_free,gmet,ext_kg_k,ext_kinpw,kpoint,ext_npw_k,0,0)

          ! Compute kinetic energy of each band
          do iband=nband_k+1,this%mband
            cwavef(1:2,1:ext_npw_k*my_nspinor)= &
            & this%cg(:,1+(iband-1)*ext_npw_k*my_nspinor+ext_icg:iband*ext_npw_k*my_nspinor+ext_icg)
            call meanvalue_g(dotr,ext_kinpw,0,istwf_k,mpi_enreg,ext_npw_k,my_nspinor,cwavef,cwavef,0)
            this%e_kinetic=this%e_kinetic+wtk(ikpt)*this%occ(iband+ext_bdtot_index)*dotr
          end do

          ! Increment indexes
          ext_bdtot_index=ext_bdtot_index+this%mband
          if(mkmem/=0) then
            ext_icg=ext_icg+ext_npw_k*my_nspinor*this%mband
            ext_ikg=ext_ikg+ext_npw_k
          end if

          ABI_FREE(cwavef)
          ABI_FREE(ext_kinpw)
          ABI_FREE(ext_kg_k)
        end do
      end do
    end if

    ! Computes extfpmd contribution to kinetic energy using a sum
    ! of Fermi gas contributions for each point of the fftf grid.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==10) then
      ABI_MALLOC(gamma_hybrid_tf,(this%nfft,this%nspden))
      ABI_MALLOC(xcut_hybrid_tf,(this%nfft,this%nspden))
      gamma_hybrid_tf(:,:)=(fermie-this%vtrial(:,:))/tsmear
      xcut_hybrid_tf(:,:)=(this%e_bcut-this%vtrial(:,:))/tsmear
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


    ! Computes the double counting term from the shiftfactor, and
    ! from the contributions to the kinetic energy and
    ! the number of electrons
    if(this%version==1.or.this%version==2.or.this%version==3.or.this%version==4.or.this%version==5) then
      this%edc_kinetic=this%e_kinetic+this%nelect*this%shiftfactor
    else if(this%version==10) then
      this%edc_kinetic=this%e_kinetic+sum(this%nelectarr(:,:)*this%vtrial(:,:)/(this%nfft*this%nspden))
    end if

    if((this%e_bcut.lt.this%shiftfactor).and.&
    & (this%version==1.or.this%version==3.or.this%version==10)) then
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
  subroutine compute_entropy(this,fermie,tsmear,nkpt,nsppol,nspinor,wtk,nband)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: nkpt,nsppol,nspinor
    real(dp),intent(in) :: fermie,tsmear
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)

    ! Local variables -------------------------
    ! Scalars
    integer :: ii,ifft,ispden,index_tot,isppol,ikpt,iband,nband_k
    real(dp) :: ix,step,factor,fn,gamma,maxocc
    ! Arrays
    real(dp),dimension(:),allocatable :: valuesent
    real(dp),dimension(:,:),allocatable :: gamma_hybrid_tf
    real(dp),dimension(:,:),allocatable :: step_hybrid_tf

    ! *********************************************************************
    maxocc=two/(nsppol*nspinor)
    this%entropy=zero

    ! Computes extfpmd contribution to the entropy integrating
    ! over accessible states with Fermi-Dirac complete integrals and
    ! substracting 0 to bcut contribution with numeric integration.
    if(this%version==2.or.this%version==4) then
      factor=sqrt(2.)/(PI*PI)*this%ucvol*tsmear**(2.5)
      gamma=(fermie-this%shiftfactor)/tsmear
      ABI_MALLOC(valuesent,(this%bcut+1))

      !$OMP PARALLEL DO PRIVATE(fn,ix) SHARED(valuesent)
      do ii=1,this%bcut+1
        ix=dble(ii)-one
        fn=fermi_dirac(extfpmd_e_fg(ix,this%ucvol)+this%shiftfactor,fermie,tsmear)
        if(fn>tol16.and.(one-fn)>tol16) then
          valuesent(ii)=-two*(fn*log(fn)+(one-fn)*log(one-fn))
        else
          valuesent(ii)=zero
        end if
      end do
      !$OMP END PARALLEL DO

      ! We need at least 6 elements in valuesent to call simpson function.
      if(size(valuesent)>=6) then
        this%entropy=5./3.*factor*dip32(gamma)/tsmear-&
        & gamma*factor*dip12(gamma)/tsmear-&
        simpson(one,valuesent)
      end if
      ABI_FREE(valuesent)
    end if

    ! Computes extfpmd contribution to the entropy integrating
    ! over energy with Fermi-Dirac complete integrals and
    ! substracting 0 to bcut contribution with numeric integration.
    if(this%version==1.or.this%version==3) then
      factor=sqrt(2.)/(PI*PI)*this%ucvol*tsmear**(2.5)
      gamma=(fermie-this%shiftfactor)/tsmear
      ABI_MALLOC(valuesent,(this%bcut+1))

      step=(this%e_bcut-this%shiftfactor)/(this%bcut)
      !$OMP PARALLEL DO PRIVATE(fn,ix) SHARED(valuesent)
      do ii=1,this%bcut+1
        ix=this%shiftfactor+(dble(ii)-one)*step
        fn=fermi_dirac(ix,fermie,tsmear)
        if(fn>tol16.and.(one-fn)>tol16) then
          valuesent(ii)=-(fn*log(fn)+(one-fn)*log(one-fn))*&
          & extfpmd_dos(ix,this%shiftfactor,this%ucvol)
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
      ABI_FREE(valuesent)
    end if

    ! Computes extended pw contribution to nelect summing
    ! over ext pw states from nband_k to this%mband
    if(this%version==5) then
      index_tot=0
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          do iband=nband_k+1,this%mband
            fn=this%occ(iband+index_tot)/maxocc
            this%entropy=this%entropy-two*(fn*log(fn)+(one-fn)*log(one-fn))
          end do
          index_tot=index_tot+this%mband
        end do
      end do
    end if

    ! Computes extfpmd contribution to the entropy using a sum
    ! of Fermi gas contributions for each point of the fftf grid,
    ! as we do for version=1 and version=2.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==10) then
      ABI_MALLOC(valuesent,(this%bcut+1))
      ABI_MALLOC(gamma_hybrid_tf,(this%nfft,this%nspden))
      ABI_MALLOC(step_hybrid_tf,(this%nfft,this%nspden))
      gamma_hybrid_tf(:,:)=(fermie-this%vtrial(:,:))/tsmear
      step_hybrid_tf(:,:)=(this%e_bcut-this%vtrial(:,:))/(this%bcut)
      this%entropy=zero
      factor=sqrt(2.)/(PI*PI)*this%ucvol*tsmear**(2.5)

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
      ABI_FREE(valuesent)
    end if
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
  !!  shiftfactor=energy shift factor
  !!  ucvol=unit cell volume (bohr^3)
  !!
  !! OUTPUT
  !!  extfpmd_dos=value of free particle density of states at given energy
  !!
  !! SOURCE
  function extfpmd_dos(energy,shiftfactor,ucvol)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: energy,shiftfactor,ucvol
    real(dp) :: extfpmd_dos

    ! *********************************************************************

    extfpmd_dos=sqrt(2.)*ucvol*sqrt(energy-shiftfactor)/(PI*PI)
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

    extfpmd_e_fg=.5*(iband*6*PI*PI/ucvol)**(2./3.)
  end function extfpmd_e_fg
  !!***

end module m_extfpmd
!!***
