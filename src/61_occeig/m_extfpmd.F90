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
!! Copyright (C) 2018-2025 ABINIT group (A. Blanchet)
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
  use m_energies,       only : energies_type
  use m_gsphere,        only : getkpgnorm
  use m_kg,             only : mkkin,kpgio
  use m_mpinfo,         only : ptabs_fourdp,proc_distrb_cycle,copy_mpi_enreg,destroy_mpi_enreg
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
    integer :: bcut,mband,nbcut,nbdbuf,nfftf,nspden,version
    real(dp) :: ebcut,edc_kinetic,e_kinetic,entropy
    real(dp) :: nelect,eshift,ucvol,el_temp,bandshift
    real(dp) :: nelect_res, nelect_respc
    real(dp),allocatable :: vtrial(:,:)
    real(dp),allocatable :: nelectarr(:,:)
    real(dp),allocatable :: bandshiftk(:)
    type(MPI_type) :: mpi_enreg
  contains
    procedure :: compute_e_kinetic
    procedure :: compute_entropy
    procedure :: compute_nelect
    procedure :: compute_eshift
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
  !!  extfpmd_eshift=pre-defined extfpmd energy shift
  !!  nbcut=number of states used to average the constant potential value
  !!  nbdbuf=Number of bands in the buffer to converge scf cycle with extfpmd models
  !!  nfftf=number of FFT fine grid points
  !!  nspden=number of spin-density components
  !!  nsppol=number of independent spin WF components
  !!  nkpt=number of k-points
  !!  occopt=option for occupancies
  !!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
  !!  tphysel="physical" electronic temperature with FD occupations
  !!  tsmear=smearing energy or temperature (if metal)
  !!  version=extfpmd implementation version
  !!  mpi_enreg=information about MPI parallelization
  !!  extfpmd_mband=number of extfpmd bands
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! SOURCE
  subroutine init(this,mband,extfpmd_eshift,nbcut,nbdbuf,nfftf,nspden,&
  & nsppol,nkpt,occopt,rprimd,tphysel,tsmear,version,mpi_enreg,extfpmd_mband)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: mband,nbcut,nbdbuf,nfftf,nspden
    integer,intent(in) :: nsppol,nkpt,version,extfpmd_mband,occopt
    real(dp),intent(in) :: extfpmd_eshift,tphysel,tsmear
    type(MPI_type),intent(in) :: mpi_enreg
    ! Arrays
    real(dp),intent(in) :: rprimd(3,3)

    ! Local variables -------------------------
    ! Arrays
    real(dp) :: gprimd(3,3),rmet(3,3),gmet(3,3)

    ! *********************************************************************

    this%bcut=mband-nbdbuf
    this%nbcut=nbcut
    this%mband=extfpmd_mband
    this%nbdbuf=nbdbuf
    this%version=version
    ABI_MALLOC(this%vtrial,(nfftf,nspden))
    this%vtrial(:,:)=zero
    this%nfftf=nfftf
    this%nspden=nspden
    this%ebcut=zero
    this%edc_kinetic=zero
    this%e_kinetic=zero
    this%entropy=zero
    this%nelect=zero
    this%nelect_res=zero
    this%nelect_respc=zero
    this%bandshift=zero
    ABI_MALLOC(this%bandshiftk,(nkpt*nsppol))
    this%bandshiftk(:)=zero
    this%eshift=extfpmd_eshift
    call metric(gmet,gprimd,-1,rmet,rprimd,this%ucvol)
    this%el_temp=merge(tphysel,tsmear,tphysel>tol8.and.occopt/=3.and.occopt/=9)

    if(this%version==5) then
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

    if(this%version==5) then
      call destroy_mpi_enreg(this%mpi_enreg)
    end if

    this%vtrial(:,:)=zero
    ABI_FREE(this%vtrial)
    if(allocated(this%nelectarr)) then
      this%nelectarr(:,:)=zero
      ABI_FREE(this%nelectarr)
    end if
    this%bandshiftk(:)=zero
    ABI_FREE(this%bandshiftk)
    this%nfftf=0
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
    this%nelect_res=zero
    this%nelect_respc=zero
    this%bandshift=zero
    this%eshift=zero
    this%ucvol=zero
    this%el_temp=zero
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
  !!  nfftf=number of FFT fine grid points
  !!  nkpt=number of k points
  !!  nsppol=1 for unpolarized, 2 for spin-polarized
  !!  nspden=number of spin-density components
  !!  wtk(nkpt)=k point weights
  !!  vtrial(nfftf,nspden)=GS potential on the fine grid (Hartree)
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! SOURCE
  subroutine compute_eshift(this,eigen,eknk,mband,nband,nfftf,nkpt,nsppol,nspden,wtk,vtrial)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: mband,nfftf,nkpt,nsppol,nspden
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
    real(dp),intent(in) :: eknk(mband*nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)
    real(dp),intent(in) :: vtrial(nfftf,nspden)

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
      this%eshift=sum(this%vtrial)/(nfftf*nspden)
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
  !!  nband(nkpt*nsppol)=desired number of bands at each k point
  !!  nelect=number of electrons per unit cell
  !!  nkpt=number of k points
  !!  nspinor=number of spinor components
  !!  nsppol=1 for unpolarized, 2 for spin-polarized
  !!  wtk(nkpt)=k point weights
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!  nelect=number of electrons per unit cell
  !!
  !! SOURCE
  subroutine compute_nelect(this,fermie,nband,nelect,nkpt,nspinor,nsppol,wtk)
    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: nkpt,nsppol,nspinor
    real(dp),intent(in) :: fermie
    real(dp),intent(inout) :: nelect
    class(extfpmd_type),intent(inout) :: this
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)

    ! Local variables -------------------------
    ! Scalars
    integer :: ifft,ispden,isppol,ikpt,iband,nband_k,ierr
    real(dp) :: factor,gamma,xcut,fn,maxocc,nelect_tmp
    ! Arrays
    real(dp),allocatable :: gamma_hybrid_tf(:,:)
    real(dp),allocatable :: xcut_hybrid_tf(:,:)

    ! *********************************************************************

    maxocc=two/(nsppol*nspinor)
    factor=dsqrt(two)/(PI*PI)*this%ucvol*this%el_temp**(1.5)
    gamma=(fermie-this%eshift)/this%el_temp
    nelect_tmp=zero

    ! Computes extfpmd contribution to nelect integrating
    ! over accessible states from bcut to infinity with
    ! order 1/2 incomplete Fermi-Dirac integral.
    if(this%version==2.or.this%version==4) then
      xcut=extfpmd_e_fg(one*this%bcut+this%bandshift,this%ucvol)/this%el_temp
      if(one*this%bcut+this%bandshift.lt.zero) xcut=zero
      nelect=nelect+factor*djp12(xcut,gamma)
    end if

    ! Computes extfpmd contribution to nelect integrating
    ! over energy from ebcut to infinity with order 1/2
    ! incomplete Fermi-Dirac integral.
    if(this%version==1.or.this%version==3) then
      xcut=(this%ebcut-this%eshift)/this%el_temp
      if(this%ebcut.lt.this%eshift) xcut=zero
      nelect=nelect+factor*djp12(xcut,gamma)
    end if

    ! Computes extfpmd contribution to nelect summing
    ! over accessible states from bcut to mband, with
    ! integer band numbers. Total number of bands
    ! is controlled with the input variable extfpmd_nband.
    if(this%version==5) then
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          if(proc_distrb_cycle(this%mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,this%mpi_enreg%me_kpt)) cycle
          do iband=nband_k-this%nbdbuf+1,this%mband
            fn=fermi_dirac(extfpmd_e_fg(one*iband+this%bandshiftk(ikpt+(isppol-1)*nkpt),this%ucvol)+this%eshift,fermie,this%el_temp)
            nelect_tmp=nelect_tmp+wtk(ikpt)*maxocc*fn
          end do
        end do
      end do
      call xmpi_sum(nelect_tmp,this%mpi_enreg%comm_kpt,ierr)
      nelect=nelect+nelect_tmp
    end if

    ! Computes extfpmd contribution to nelect using a sum
    ! of Fermi gas contributions for each point of the fftf grid.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==10) then
      ABI_MALLOC(gamma_hybrid_tf,(this%nfftf,this%nspden))
      ABI_MALLOC(xcut_hybrid_tf,(this%nfftf,this%nspden))
      if(.not.allocated(this%nelectarr)) ABI_MALLOC(this%nelectarr,(this%nfftf,this%nspden))
      this%nelectarr(:,:)=zero
      gamma_hybrid_tf(:,:)=(fermie-this%vtrial(:,:))/this%el_temp
      xcut_hybrid_tf(:,:)=(this%ebcut-this%vtrial(:,:))/this%el_temp
      if(ANY(this%ebcut.lt.this%vtrial(:,:))) xcut_hybrid_tf(:,:)=zero

      !$OMP PARALLEL DO
      do ifft=1,this%nfftf
        do ispden=1,this%nspden
          this%nelectarr(ifft,ispden)=factor*djp12(xcut_hybrid_tf(ifft,ispden),gamma_hybrid_tf(ifft,ispden))
        end do
      end do
      !$OMP END PARALLEL DO

      nelect=nelect+sum(this%nelectarr)/(this%nfftf*this%nspden)
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
  !!  nkpt=number of k points
  !!  nspinor=number of spinor components
  !!  nsppol=1 for unpolarized, 2 for spin-polarized
  !!  nband(nkpt*nsppol)=desired number of bands at each k point
  !!  wtk(nkpt)=k point weights
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! SOURCE
  subroutine compute_e_kinetic(this,fermie,nkpt,nspinor,nsppol,nband,wtk)
    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: nkpt,nspinor,nsppol
    class(extfpmd_type),intent(inout) :: this
    real(dp),intent(in) :: fermie
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)

    ! Local variables -------------------------
    ! Scalars
    logical :: cut_warn=.false.
    integer :: ikpt,isppol,nband_k,iband,ierr,ifft,ispden
    real(dp) :: factor,gamma,xcut,dotr
    real(dp) :: e_kinetic_hybrid_tf,maxocc,fn
    character(len=500) :: msg
    ! Arrays
    real(dp),allocatable :: gamma_hybrid_tf(:,:),xcut_hybrid_tf(:,:)

    ! *********************************************************************

    dotr=zero
    maxocc=two/(nsppol*nspinor)
    this%e_kinetic=zero
    factor=dsqrt(two)/(PI*PI)*this%ucvol*this%el_temp**(2.5)
    gamma=(fermie-this%eshift)/this%el_temp

    ! Computes extfpmd contribution to kinetic energy integrating
    ! over accessible states from bcut to infinity with
    ! order 3/2 incomplete Fermi-Dirac integral.
    if(this%version==2.or.this%version==4) then
      xcut=extfpmd_e_fg(one*this%bcut+this%bandshift,this%ucvol)/this%el_temp
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
      xcut=(this%ebcut-this%eshift)/this%el_temp
      if(this%ebcut.lt.this%eshift) then
        cut_warn=.true.
        xcut=zero
      end if
      this%e_kinetic=this%e_kinetic+factor*djp32(xcut,gamma)
    end if

    ! Computes extfpmd contribution to kinetic energy summing
    ! over accessible states from bcut to mband, with
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
            fn=fermi_dirac(dotr+this%eshift,fermie,this%el_temp)
            this%e_kinetic=this%e_kinetic+wtk(ikpt)*maxocc*fn*dotr
          end do
        end do
      end do
      call xmpi_sum(this%e_kinetic,this%mpi_enreg%comm_kpt,ierr)
    end if

    ! Computes extfpmd contribution to kinetic energy using a sum
    ! of Fermi gas contributions for each point of the fftf grid.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==10) then
      ABI_MALLOC(gamma_hybrid_tf,(this%nfftf,this%nspden))
      ABI_MALLOC(xcut_hybrid_tf,(this%nfftf,this%nspden))
      gamma_hybrid_tf(:,:)=(fermie-this%vtrial(:,:))/this%el_temp
      xcut_hybrid_tf(:,:)=(this%ebcut-this%vtrial(:,:))/this%el_temp
      e_kinetic_hybrid_tf=zero

      !$OMP PARALLEL DO REDUCTION (+:e_kinetic_hybrid_tf)
      do ifft=1,this%nfftf
        do ispden=1,this%nspden
          e_kinetic_hybrid_tf=e_kinetic_hybrid_tf+factor*djp32(xcut_hybrid_tf(ifft,ispden),gamma_hybrid_tf(ifft,ispden))/&
          & (this%nfftf*this%nspden)
        end do
      end do
      !$OMP END PARALLEL DO

      this%e_kinetic=e_kinetic_hybrid_tf
      gamma_hybrid_tf(:,:)=zero
      xcut_hybrid_tf(:,:)=zero
      ABI_FREE(gamma_hybrid_tf)
      ABI_FREE(xcut_hybrid_tf)
    end if

    ! Computes the double counting term from the eshift, and
    ! from the contributions to the kinetic energy and
    ! the number of electrons
    if(this%version==10) then
      this%edc_kinetic=this%e_kinetic+sum(this%nelectarr(:,:)*this%vtrial(:,:)/(this%nfftf*this%nspden))
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
  !!  nkpt=number of k points
  !!  nsppol=1 for unpolarized, 2 for spin-polarized
  !!  nspinor=number of spinor components
  !!  wtk(nkpt)=k point weights
  !!  nband(nkpt*nsppol)=desired number of bands at each k point
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!  entropy_extfpmd=extfpmd contribution to the entropy
  !!
  !! SOURCE
  subroutine compute_entropy(this,entropy_extfpmd,fermie,nkpt,nsppol,nspinor,wtk,nband)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: nkpt,nsppol,nspinor
    real(dp),intent(in) :: fermie
    real(dp),intent(out) :: entropy_extfpmd
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)

    ! Local variables -------------------------
    ! Scalars
    integer :: ii,ifft,ispden,isppol,ikpt,iband,nband_k,ierr
    real(dp) :: ix,step,factor,fn,gamma,maxocc
    ! Arrays
    real(dp),dimension(:),allocatable :: valuesent
    real(dp),dimension(:,:),allocatable :: gamma_hybrid_tf
    real(dp),dimension(:,:),allocatable :: step_hybrid_tf

    ! *********************************************************************
    maxocc=two/(nsppol*nspinor)
    this%entropy=zero
    factor=dsqrt(two)/(PI*PI)*this%ucvol*this%el_temp**(2.5)
    gamma=(fermie-this%eshift)/this%el_temp
    ABI_MALLOC(valuesent,(this%bcut+1))

    ! Computes extfpmd contribution to the entropy integrating
    ! over accessible states with Fermi-Dirac complete integrals and
    ! substracting 0 to bcut contribution with numeric integration.
    if(this%version==2.or.this%version==4) then
      step=(dble(this%bcut)+this%bandshift)/(this%bcut)
      !$OMP PARALLEL DO PRIVATE(fn,ix) SHARED(valuesent)
      do ii=1,this%bcut+1
        ix=(dble(ii)-one)*step
        fn=fermi_dirac(extfpmd_e_fg(ix,this%ucvol)+this%eshift,fermie,this%el_temp)
        if(fn>tol16.and.(one-fn)>tol16) then
          valuesent(ii)=-maxocc*(fn*log(fn)+(one-fn)*log(one-fn))
        else
          valuesent(ii)=zero
        end if
      end do
      !$OMP END PARALLEL DO

      ! We need at least 6 elements in valuesent to call simpson function.
      if(size(valuesent)>=6) then
        this%entropy=5./3.*factor*dip32(gamma)/this%el_temp-&
        & gamma*factor*dip12(gamma)/this%el_temp-&
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
        fn=fermi_dirac(ix,fermie,this%el_temp)
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
        this%entropy=5./3.*factor*dip32(gamma)/this%el_temp-&
        & gamma*factor*dip12(gamma)/this%el_temp-simpson(step,valuesent)
      end if
    end if

    ! Computes extfpmd contribution to the entropy summing
    ! over accessible states from bcut to mband, with
    ! integer band numbers. Total number of bands
    ! is controlled with the input variable extfpmd_nband.
    if(this%version==5) then
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          if(proc_distrb_cycle(this%mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,this%mpi_enreg%me_kpt)) cycle
          do iband=nband_k-this%nbdbuf+1,this%mband
            fn=fermi_dirac(extfpmd_e_fg(one*iband+this%bandshiftk(ikpt+(isppol-1)*nkpt),this%ucvol)+this%eshift,fermie,this%el_temp)
            this%entropy=this%entropy-wtk(ikpt)*maxocc*(fn*log(fn)+(one-fn)*log(one-fn))/nsppol
          end do
        end do
      end do
      call xmpi_sum(this%entropy,this%mpi_enreg%comm_kpt,ierr)
    end if

    ! Computes extfpmd contribution to the entropy using a sum
    ! of Fermi gas contributions for each point of the fftf grid,
    ! as we do for version=1 and version=2.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==10) then
      ABI_MALLOC(gamma_hybrid_tf,(this%nfftf,this%nspden))
      ABI_MALLOC(step_hybrid_tf,(this%nfftf,this%nspden))
      gamma_hybrid_tf(:,:)=(fermie-this%vtrial(:,:))/this%el_temp
      step_hybrid_tf(:,:)=(this%ebcut-this%vtrial(:,:))/(this%bcut)
      this%entropy=zero

      do ifft=1,this%nfftf
        do ispden=1,this%nspden
          !$OMP PARALLEL DO PRIVATE(ix,fn) SHARED(valuesent)
          do ii=1,this%bcut+1
            ix=this%vtrial(ifft,ispden)+(dble(ii)-one)*step_hybrid_tf(ifft,ispden)
            fn=fermi_dirac(ix,fermie,this%el_temp)
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
            this%entropy=this%entropy+(5./3.*factor*dip32(gamma_hybrid_tf(ifft,ispden))/this%el_temp-&
            & gamma_hybrid_tf(ifft,ispden)*factor*dip12(gamma_hybrid_tf(ifft,ispden))/this%el_temp-&
            & simpson(step_hybrid_tf(ifft,ispden),valuesent))/(this%nfftf*this%nspden)
          end if
        end do
      end do

      gamma_hybrid_tf(:,:)=zero
      step_hybrid_tf(:,:)=zero
      ABI_FREE(step_hybrid_tf)
      ABI_FREE(gamma_hybrid_tf)
    end if
    ABI_FREE(valuesent)
    entropy_extfpmd=this%entropy
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
    if(.not.(dtset%occopt>=3.and.dtset%occopt<=9)) then
      write(msg,'(3a)') "ExtFPMD routines need metallic occupation option.",ch10,&
      & "Action: Set occopt input variable to a value >= 3 and <= 9."
      ABI_ERROR(msg)
    else if((dtset%useextfpmd==2.or.dtset%useextfpmd==3).and.(dtset%extfpmd_nbcut>dtset%mband)) then
      write(msg,'(3a,i0,a,i0,3a)') "Not enough bands to activate ExtFPMD routines.",ch10,&
      & "extfpmd_nbcut = ",dtset%extfpmd_nbcut," should be less than or equal to nband = ",dtset%mband,".",ch10,&
      & "Action: Increase nband or decrease extfpmd_nbcut."
      ABI_ERROR(msg)
    else if((dtset%useextfpmd==2.or.dtset%useextfpmd==3).and.(dtset%extfpmd_nbdbuf+dtset%extfpmd_nbcut>dtset%mband)) then
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
