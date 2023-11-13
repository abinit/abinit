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
  use m_kg,             only : mkkin,kpgio
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
    real(dp) :: ecut
    real(dp),allocatable :: kg(:,:)
  contains
    procedure :: compute_e_kinetic
    procedure :: compute_entropy
    procedure :: compute_nelect
    procedure :: compute_kg
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
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! SOURCE
  subroutine init(this,mband,nbcut,nbdbuf,nfft,nspden,rprimd,version)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: mband,nbcut,nbdbuf,nfft,nspden,version
    ! Arrays
    real(dp),intent(in) :: rprimd(3,3)

    ! Local variables -------------------------
    ! Arrays
    real(dp) :: gprimd(3,3),rmet(3,3), gmet(3,3)

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
    this%ecut=zero
    call metric(gmet,gprimd,-1,rmet,rprimd,this%ucvol)
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

    this%vtrial(:,:)=zero
    ABI_FREE(this%vtrial)
    this%nelectarr(:,:)=zero
    ABI_FREE(this%nelectarr)
    this%nfft=0
    this%nspden=0
    this%bcut=0
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

  !!****f* ABINIT/m_extfpmd/compute_kg
  !! NAME
  !!  compute_kg
  !!
  !! FUNCTION
  !!  Computes the extended plane wave basis set, computing extended
  !!  plane wave cutoff energy
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
  subroutine compute_kg(this)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    ! integer,intent(in) :: mband,me,nkpt,nsppol
    ! Arrays
    ! integer,intent(in) :: nband(nkpt*nsppol)
    ! real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
    ! real(dp),intent(in) :: eknk(mband*nkpt*nsppol)
    ! real(dp),intent(in) :: wtk(nkpt)

    ! Local variables -------------------------
    ! Scalars
    ! integer :: band_index,ii,ikpt,isppol,nband_k
    ! real(dp) :: abs_err,rel_err
    ! character(len=500) :: msg

    ! *********************************************************************
    if(this%version==5) then
      
    end if
  end subroutine compute_kg

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
  subroutine compute_nelect(this,fermie,nelect,tsmear)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: fermie,tsmear
    real(dp),intent(inout) :: nelect
    class(extfpmd_type),intent(inout) :: this

    ! Local variables -------------------------
    ! Scalars
    integer :: ifft,ispden
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
  subroutine compute_e_kinetic(this,fermie,tsmear)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    real(dp),intent(in) :: fermie,tsmear

    ! Local variables -------------------------
    ! Scalars
    integer :: ifft,ispden
    real(dp) :: factor,gamma,xcut
    real(dp) :: e_kinetic_hybrid_tf
    character(len=500) :: msg
    ! Arrays
    real(dp),allocatable :: gamma_hybrid_tf(:,:)
    real(dp),allocatable :: xcut_hybrid_tf(:,:)

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
    if(this%version==1.or.this%version==2.or.this%version==3.or.this%version==4) then
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
  subroutine compute_entropy(this,fermie,tsmear)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    real(dp),intent(in) :: fermie,tsmear

    ! Local variables -------------------------
    ! Scalars
    integer :: ii,ifft,ispden
    real(dp) :: ix,step,factor,fn,gamma
    ! Arrays
    real(dp),dimension(:),allocatable :: valuesent
    real(dp),dimension(:,:),allocatable :: gamma_hybrid_tf
    real(dp),dimension(:,:),allocatable :: step_hybrid_tf

    ! *********************************************************************

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
