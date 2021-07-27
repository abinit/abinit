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
!! Copyright (C) 2018-2021 ABINIT group (A. Blanchet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! TODO
!! 1) Add contribution to conductivity.
!! 2) Smooth the contributions.
!!
!! CHILDREN
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
  public :: dip12,djp12,dip32,djp32,extfpmd_dos,extfpmd_e_fg
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
    integer :: bcut,nbcut,nfftf,nspden,version
    real(dp) :: e_bcut,edc_kinetic,e_kinetic,entropy
    real(dp) :: nelect,shiftfactor,ucvol
    real(dp),allocatable :: vtrial(:,:)
  contains
    procedure :: compute_e_kinetic
    procedure :: compute_entropy
    procedure :: compute_nelect
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
  !!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
  !!  version=extfpmd implementation version
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! PARENTS
!!
  !! CHILDREN
!!
  !! SOURCE
  subroutine init(this,mband,nbcut,nfftf,nspden,rprimd,version)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    integer,intent(in) :: mband,nbcut,nfftf,nspden,version
    ! Arrays
    real(dp),intent(in) :: rprimd(3,3)

    ! Local variables -------------------------
    ! Arrays
    real(dp) :: gprimd(3,3),rmet(3,3), gmet(3,3)

    ! *********************************************************************

    this%bcut=mband
    this%nbcut=nbcut
    this%version=version
    ABI_MALLOC(this%vtrial,(nfftf,nspden))
    this%vtrial(:,:)=zero
    this%nfftf=nfftf
    this%nspden=nspden
    this%e_bcut=zero
    this%edc_kinetic=zero
    this%e_kinetic=zero
    this%entropy=zero
    this%nelect=zero
    this%shiftfactor=zero
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
  !! PARENTS
!!
  !! CHILDREN
!!
  !! SOURCE
  subroutine destroy(this)

    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this

    ! *********************************************************************

    this%vtrial(:,:)=zero
    ABI_FREE(this%vtrial)
    this%nfftf=0
    this%nspden=0
    this%bcut=0
    this%nbcut=0
    this%version=1
    this%e_bcut=zero
    this%edc_kinetic=zero
    this%e_kinetic=zero
    this%entropy=zero
    this%nelect=zero
    this%shiftfactor=zero
    this%ucvol=zero
  end subroutine destroy
  !!***

  !!****f* ABINIT/m_extfpmd/compute_shiftfactor
  !! NAME
  !!  compute_shiftfactor
  !!
  !! FUNCTION
  !!  Compute the energy shift factor $U_0$ corresponding to constant
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
  !! PARENTS
!!
  !! CHILDREN
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

    ! Compute U_0 from the sum of local
    ! potentials (vtrial), averaging over all space.
    ! Simplest and most precise way to evaluate U_0.
    if(this%version==1) then
      this%shiftfactor=sum(this%vtrial)/(this%nfftf*this%nspden)

      ! Compute the relative error of the model vs last eigenvalues
      if(me==0) then
        band_index=0
        rel_err=zero
        abs_err=zero
        do isppol=1,nsppol
          do ikpt=1,nkpt
            nband_k=nband(ikpt+(isppol-1)*nkpt)
            rel_err=rel_err+wtk(ikpt)*abs((eigen(band_index+nband_k)-&
            & extfpmd_e_fg(dble(nband_k),this%ucvol)-this%shiftfactor)/&
            & eigen(band_index+nband_k))
            abs_err=abs_err+wtk(ikpt)*abs(eigen(band_index+nband_k)-&
            & extfpmd_e_fg(dble(nband_k),this%ucvol)-this%shiftfactor)
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

    ! Compute U_0^{HEG} from the difference between
    ! eigenvalues and Fermi gas energies, averaged
    ! over lasts nbcut bands.
    if(this%version==2) then
      this%shiftfactor=zero
      band_index=0
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          do ii=nband_k-this%nbcut+1,nband_k
            this%shiftfactor=this%shiftfactor+&
            & wtk(ikpt)*(eigen(band_index+ii)-extfpmd_e_fg(dble(ii),this%ucvol))
          end do
          band_index=band_index+nband_k
        end do
      end do
      this%shiftfactor=this%shiftfactor/this%nbcut
    end if

    ! Compute U_0^K from the difference between
    ! eigenvalues and kinetic energies, averaged
    ! over lasts nbcut bands.
    if(this%version==3) then
      this%shiftfactor=zero
      this%e_bcut=0
      band_index=0
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          do ii=nband_k-this%nbcut+1,nband_k
            this%shiftfactor=this%shiftfactor+&
            & wtk(ikpt)*(eigen(band_index+ii)-eknk(band_index+ii))
          end do
          this%e_bcut=this%e_bcut+wtk(ikpt)*eigen(band_index+nband_k)
          band_index=band_index+nband_k
        end do
      end do
      this%shiftfactor=this%shiftfactor/this%nbcut
    end if
  end subroutine compute_shiftfactor
  !!***

  !!****f* ABINIT/m_extfpmd/compute_nelect
  !! NAME
  !!  compute_nelect
  !!
  !! FUNCTION
  !!  Compute the value of the integral corresponding to the missing
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
  !! PARENTS
!!
  !! CHILDREN
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
    integer :: ifftf,ispden
    real(dp) :: factor,gamma,xcut

    ! *********************************************************************

    factor=sqrt(2.)/(PI*PI)*this%ucvol*tsmear**(1.5)

    ! Compute extfpmd contribution to nelect integrating
    ! over accessible states from bcut to infinity with
    ! order 1/2 incomplete Fermi-Dirac integral.
    if(this%version==1.or.this%version==2) then
      gamma=(fermie-this%shiftfactor)/tsmear
      xcut=extfpmd_e_fg(dble(this%bcut),this%ucvol)/tsmear
      nelect=nelect+factor*djp12(xcut,gamma)
    end if

    ! Compute extfpmd contribution to nelect integrating
    ! over energy from e_bcut to infinity with order 1/2
    ! incomplete Fermi-Dirac integral.
    if(this%version==3) then
      gamma=(fermie-this%shiftfactor)/tsmear
      xcut=(this%e_bcut-this%shiftfactor)/tsmear
      nelect=nelect+factor*djp12(xcut,gamma)
    end if

    ! Compute extfpmd contribution to nelect using a sum
    ! of Fermi gas contributions for each point of the fftf grid.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==4) then
      do ifftf=1,this%nfftf
        do ispden=1,this%nspden
          gamma=(fermie-this%vtrial(ifftf,ispden))/tsmear
          xcut=extfpmd_e_fg(dble(this%bcut),this%ucvol)/tsmear
          nelect=nelect+factor*djp12(xcut,gamma)/(this%nfftf*this%nspden)
        end do
      end do
    end if
  end subroutine compute_nelect
  !!***

  !!****f* ABINIT/m_extfpmd/compute_e_kinetic
  !! NAME
  !!  compute_e_kinetic
  !!
  !! FUNCTION
  !!  Compute the value of the integral corresponding to the missing
  !!  kinetic energy contribution of free electrons after band cut,
  !!  with an order 3/2 incomplete Fermi-Dirac integral.
  !!
  !! INPUTS
  !!  this=extfpmd_type object concerned
  !!  fermie=chemical potential (Hartree)
  !!  nfftf=number of points in the fine FFT mesh (for this processor)
  !!  nspden=number of spin-density components
  !!  tsmear=smearing width (or temperature)
  !!  vtrial(nfft,nspden)= trial potential (Hartree)
  !!
  !! OUTPUT
  !!  this=extfpmd_type object concerned
  !!
  !! PARENTS
!!
  !! CHILDREN
!!
  !! SOURCE
  subroutine compute_e_kinetic(this,fermie,nfftf,nspden,tsmear,vtrial)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    real(dp),intent(in) :: fermie,tsmear
    integer,intent(in) :: nfftf,nspden
    ! Arrays
    real(dp),intent(in) :: vtrial(nfftf,nspden)

    ! Local variables -------------------------
    ! Scalars
    integer :: ifftf,ispden
    real(dp) :: factor,gamma,xcut

    ! *********************************************************************

    this%e_kinetic=zero
    factor=sqrt(2.)/(PI*PI)*this%ucvol*tsmear**(2.5)

    ! Compute extfpmd contribution to kinetic energy integrating
    ! over accessible states from bcut to infinity with
    ! order 3/2 incomplete Fermi-Dirac integral.
    if(this%version==1.or.this%version==2) then
      gamma=(fermie-this%shiftfactor)/tsmear
      xcut=extfpmd_e_fg(dble(this%bcut),this%ucvol)/tsmear
      this%e_kinetic=this%e_kinetic+factor*djp32(xcut,gamma)
    end if

    ! Compute extfpmd contribution to kinetic energy integrating
    ! over energy from e_bcut to infinity with order 3/2
    ! incomplete Fermi-Dirac integral.
    if(this%version==3) then
      gamma=(fermie-this%shiftfactor)/tsmear
      xcut=(this%e_bcut-this%shiftfactor)/tsmear
      this%e_kinetic=this%e_kinetic+factor*djp32(xcut,gamma)
    end if

    ! Compute extfpmd contribution to kinetic energy using a sum
    ! of Fermi gas contributions for each point of the fftf grid.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==4) then
      do ifftf=1,this%nfftf
        do ispden=1,this%nspden
          gamma=(fermie-this%vtrial(ifftf,ispden))/tsmear
          xcut=extfpmd_e_fg(dble(this%bcut),this%ucvol)/tsmear
          this%e_kinetic=this%e_kinetic+factor*djp32(xcut,gamma)/&
          & (this%nfftf*this%nspden)
        end do
      end do
    end if

    ! Compute the double counting term of the contribution to
    ! kinetic energy from the vtrial potential.
    this%edc_kinetic=zero
    do ispden=1,nspden
      do ifftf=1,nfftf
        this%edc_kinetic=this%edc_kinetic+vtrial(ifftf,ispden)
      end do
    end do
    this%edc_kinetic=this%edc_kinetic*this%nelect/nspden/nfftf/nspden
  end subroutine compute_e_kinetic
  !!***

  !!****f* ABINIT/m_extfpmd/compute_entropy
  !! NAME
  !!  compute_entropy
  !!
  !! FUNCTION
  !!  Compute the value of the integral corresponding to the missing
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
  !! PARENTS
!!
  !! CHILDREN
!!
  !! SOURCE
  subroutine compute_entropy(this,fermie,tsmear)
    ! Arguments -------------------------------
    ! Scalars
    class(extfpmd_type),intent(inout) :: this
    real(dp),intent(in) :: fermie,tsmear

    ! Local variables -------------------------
    ! Scalars
    integer :: ii,ifftf,ispden
    real(dp) :: ix,step,factor,fn,gamma,minocc
    ! Arrays
    real(dp),dimension(:),allocatable :: valuesent

    ! *********************************************************************

    this%entropy=zero

    ! Compute extfpmd contribution to the entropy integrating
    ! over accessible states with Fermi-Dirac complete integrals and
    ! substracting 0 to bcut contribution with numeric integration.
    if(this%version==1.or.this%version==2) then
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

    ! Compute extfpmd contribution to the entropy integrating
    ! over energy with Fermi-Dirac complete integrals and
    ! substracting 0 to bcut contribution with numeric integration.
    if(this%version==3) then
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
        & gamma*factor*dip12(gamma)/tsmear-&
        simpson(step,valuesent)
      end if
      ABI_FREE(valuesent)
    end if

    ! Compute extfpmd contribution to the entropy using a sum
    ! of Fermi gas contributions for each point of the fftf grid,
    ! as we do for version=1 and version=2.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==4) then
      step=one
      do ifftf=1,this%nfftf
        do ispden=1,this%nspden
          ! Dynamic array find size
          ix=dble(this%bcut)
          ii=0
          fn=fermi_dirac(extfpmd_e_fg(ix,this%ucvol)+&
          & this%vtrial(ifftf,ispden),fermie,tsmear)
          minocc=tol16
          do while(fn>minocc)
            fn=fermi_dirac(extfpmd_e_fg(ix,this%ucvol)+&
            & this%vtrial(ifftf,ispden),fermie,tsmear)
            ii=ii+1
            ix=ix+step
          end do
          ABI_MALLOC(valuesent,(ii))
          ix=dble(this%bcut)
          ii=0
          fn=fermi_dirac(extfpmd_e_fg(ix,this%ucvol)+&
          & this%vtrial(ifftf,ispden),fermie,tsmear)
          do while(fn>minocc)
            fn=fermi_dirac(extfpmd_e_fg(ix,this%ucvol)+&
            & this%vtrial(ifftf,ispden),fermie,tsmear)
            ii=ii+1
            valuesent(ii)=-2*(fn*log(fn)+(1.-fn)*log(1.-fn))
            ix=ix+step
          end do
          if (ii>1) then
            this%entropy=this%entropy+simpson(step,valuesent)/&
            & (this%nfftf*this%nspden)
          end if
          ABI_FREE(valuesent)
        end do
      end do
    end if
  end subroutine compute_entropy
  !!***

  !----------------------------------------------------------------------

  !!****f* ABINIT/m_extfpmd/dip12
  !! NAME
  !!  dip12
  !!
  !! FUNCTION
  !!  Returns the complete Fermi integral of order 1/2.
  !!
  !! INPUTS
  !!  gamma=complete Fermi integral argument
  !!
  !! OUTPUT
  !!  dip12=resulting function
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function dip12(gamma)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: gamma

    ! Local variables -------------------------
    ! Scalars
    real(dp) :: d,dip12,dy

    ! *********************************************************************

    if (gamma.lt.3.) then
      dy=exp(gamma)
      if (gamma+1.9375.LE.0) then
        dip12=dy*&
        & (1.-dy*(0.35355283-dy*(0.19242767-dy*(0.12456909-dy*&
        & (8.5114507E-02-dy*4.551794E-02)))))
      else
        d=gamma-0.5
        dip12=dy*(0.677695804-d*(0.187773135+d*(2.16197521E-02-d*&
        & (9.23703807E-03+d*&
        & (1.71735167E-03-d*(6.07913775E-04+d*&
        & (1.1448629E-04-d*&
        & (4.544432E-05+d*(6.4719368E-06-d*(3.794983E-06+d*&
        & (1.7338029E-07-d*&
        & (3.5546516E-07-d*(3.7329191E-08+d*&
        & (3.3097822E-08-d*&
        & (8.3190193E-09+d*(2.2752769E-09-d*(7.836005E-10+d*&
        & (7.519551E-11-d*2.960006E-11))))))))))))))))))
      end if
    else if (gamma.lt.20.) then
      if (gamma.lt.10.) then
        d=gamma-6.5
        dip12=12.839811+d*&
        & (2.844774+d*(0.114920926-d*(3.43733039E-03-d*&
        & (2.3980356E-04-d*&
        & (2.0201888E-05-d*(1.5219883E-06-d*&
        & (6.2770524E-08+d*&
        & (4.8830336E-09-d*(2.1031164E-09-d*(5.785753E-10-d*&
        & (7.233066E-11-d*1.230727E-12)))))))))))
      else
        d=gamma-14.5
        dip12=41.7799227+d*&
        & (4.2881461+d*(7.45407825E-02-d*(8.79243296E-04-d*&
        & (2.38288861E-05-d*&
        & (8.82474867E-07-d*(3.82865217E-08-d*&
        & (1.9274292E-09-d*(1.42248669E-10-d*8.17019813E-12))))))))
      end if
    else
      d=1./gamma
      dy=gamma*sqrt(gamma)/1.329340388
      dip12=dy*(1.-d*(9.354E-07-d*(1.2338391-d*(6.77931E-03-d*&
      & 1.17871643))))
    end if
    dip12=dip12*0.88622692
  end function dip12
  !!***

  !!****f* ABINIT/m_extfpmd/djp12
  !! NAME
  !!  djp12
  !!
  !! FUNCTION
  !!  Returns the incomplete Fermi integral of order 1/2.
  !!
  !! INPUTS
  !!  xcut=lower bound of the incomplete Fermi integral
  !!  gamma=incomplete Fermi integral argument
  !!
  !! OUTPUT
  !!  djp12=resulting function
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function djp12(xcut,gamma)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: xcut,gamma
    real(dp) :: djp12

    ! Local variables -------------------------
    ! Scalars
    real(dp) :: d2h,db,dc,dd,de,df1,df2,df3,dh
    real(dp) :: ds,dt,dv,dw,dxm,dxp
    integer :: i,ind,iq,k,nm,np,nq
    ! Arrays
    real(dp) :: dq(5),df(101),dy(101)

    ! *********************************************************************

    dh=0.2D+0
    d2h=0.4D+0
    nm=101
    ind=0
    dq=(/1.D+0,2.828427124D+0,5.196152423D+0,8.D+0,1.118033989D+1/)

    djp12=0.D+0
    dxm=gamma-1.5D+1
    if (xcut.gt.dxm) then
      if (ind.eq.0) then
        do i=1,nm
          dy(i)=-1.5D+1+(i-1)*dh
          df(i)=1.D+0+dexp(dy(i))
        end do
        ind=1
      end if
      dxp=gamma+5.D+0
      if (xcut.lt.dxp) then
        dc=dxp
      else
        dc=xcut
      end if
      db=dexp(gamma-dc)
      dt=db
      do iq=1,5
        dd=iq*dc
        ds=dsqrt(dd)
        dw=1.+.3275911*ds
        dw=1.D+0/dw
        dv=dw*(.2258368458D+0+&
        & dw*(-.2521286676D+0+dw*(1.2596951294D+0+&
        & dw*(-1.2878224530D+0+dw*(.9406460699D+0)))))
        dv=dv+ds
        de=dt*dv/dq(iq)
        djp12=djp12+de
        if (dabs(de).lt.(1.D-07*djp12)) exit
        dt=-dt*db
      end do
      if (xcut.ge.dxp) return
      np=(dxp-xcut)/dh
      np=2*(np/2)
      np=nm-np
      nq=(15.-gamma)/dh
      nq=1+2*(nq/2)
      if (np.lt.nq) np=nq
      if (np.le.nm) then
        df3=0.D+0
        dt=dy(np)+gamma
        dv=(dt-xcut)/2.D+0
        df3=0.D0
        if (dt.ge.1.D-13) df3=dsqrt(dt)
        if (dabs(dv).ge.1.D-13) then
          df1=dsqrt(xcut)
          dt=df1+df3
          dw=(dv+dv)/(dt*dt)
          df2=dw*dw
          df2=df2+df2
          db=df2*(df2+7.D+0)+7.D+1
          dc=7.D+0*(1.D+1-df2)
          dc=dc*dw
          dd=-df2*(df2-2.8D+1)+1.4D+2
          dd=dd+dd
          ds=dt*((db-dc)/(1.D+0+dexp(xcut-gamma))+&
          & dd/(1.D+0+dexp(xcut+dv-gamma))+(db+dc)/df(np))
          ds=ds*dv/4.2D+2
          djp12=djp12+ds
        end if
        if (np.ne.nm) then
          ds=0.D+0
          np=np+2
          do k=np,nm,2
            df1=df3
            df3=dsqrt(dy(k)+gamma)
            dt=df1+df3
            dw=d2h/(dt*dt)
            df2=dw*dw
            df2=df2+df2
            db=df2*(df2+7.D+0)+7.D+1
            dc=7.D+0*(1.D+1-df2)
            dc=dc*dw
            dd=-df2*(df2-2.8D+1)+1.4D+2
            dd=dd+dd
            ds=ds+dt*((db-dc)/df(k-2)+dd/df(k-1)+(db+dc)/df(k))
          end do
          ds=ds*dh/4.2D+2
          djp12=djp12+ds
        end if
        if (xcut.ge.dxm) return
      end if
    end if
    djp12=dip12(gamma)-xcut*dsqrt(xcut)/1.5D+0
  end function djp12
  !!***

  !!****f* ABINIT/m_extfpmd/dip32
  !! NAME
  !!  dip32
  !!
  !! FUNCTION
  !!  Returns the complete Fermi integral of order 3/2.
  !!
  !! INPUTS
  !!  gamma=complete Fermi integral argument
  !!
  !! OUTPUT
  !!  dip32=resulting function
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function dip32(gamma)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: gamma

    ! Local variables -------------------------
    ! Scalars
    real(dp) :: d,dip32,dval

    ! *********************************************************************

    if (gamma.GT.1.75) then
      dval=gamma*gamma*SQRT(gamma)
      if (gamma.LT.4.5) then
        d=gamma-3.125
        dip32=(1.27623+0.596065*gamma+0.3*dval)*&
        & (1.0055558385-d*(5.23889494E-03+d*&
        & (3.13523144E-03-d*(3.06124286E-03-d*&
        & (1.3644667E-03-d*&
        & (4.1528384E-04-d*(8.901188E-05-d*(1.079979E-05+d*&
        & (2.29058E-06-d*(2.58985E-06-d*7.30909E-07))))))))))
      else if (gamma.LT.12.) then
        if (gamma.LT.8.) then
          d=gamma-6.25
          dip32=(2.01508+0.425775*gamma+0.3*dval)*&
          & (1.000387131-d*(3.93626295E-04+d*&
          & (2.55710115E-04-d*&
          & (1.57383494E-04-d*(5.0286036E-05-d*&
          & (1.2073559865E-05-d*&
          & (2.4909523213E-06-d*(5.244328548E-07-d*&
          & 8.0884033896E-08))))))))
        else
          d=gamma-10.
          dip32=0.3*dval*&
          & (1.064687247-d*(1.22972303E-02-d*(1.8362121E-03-d*&
          & (2.433558E-04-d*(3.018186E-05-d*(3.5694E-06-d*&
          & (4.11212E-07-d*(5.2151E-08-d*5.8424E-09))))))))
        end if
      else
        d=1./gamma
        dip32=0.30090111127*dval*&
        & (1.-d*(2.863E-06-d*(6.168876549-d*&
        & (1.740553E-02+d*(1.425257+d*2.95887)))))
      end if
    else if (gamma+0.75.LE.0) then
      d=EXP(gamma)
      dip32=d*&
      & (1.-d*(1.76775246E-01-d*(6.4124584E-02-d*&
      & (3.1027055E-02-d*(1.6797637E-02-d*&
      & (8.212636E-03-d*(2.384106E-03)))))))
    else
      d=gamma-0.5
      dip32=EXP(gamma)*(0.846691-0.128948*gamma)*&
      & (1.034064158+d*(2.778947E-03-d*&
      & (3.572502805E-02+d*(3.0411645E-03-d*&
      & (1.7380548E-03+d*(2.7756776E-04-d*&
      & (8.08302E-05+d*(1.59606E-05-d*&
      & (3.8144E-06+d*7.4446E-07)))))))))
    end if
    dip32=dip32*1.32934038
  end function dip32
  !!***

  !!****f* ABINIT/m_extfpmd/djp32
  !! NAME
  !!  djp32
  !!
  !! FUNCTION
  !!  Returns the incomplete Fermi integral of order 3/2.
  !!
  !! INPUTS
  !!  xcut=lower bound of the incomplete Fermi integral
  !!  gamma=incomplete Fermi integral argument
  !!
  !! OUTPUT
  !!  djp32=resulting function
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function djp32(xcut,gamma)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: xcut,gamma
    real(dp) :: djp32

    ! Local variables -------------------------
    ! Scalars
    real(dp) :: d2h,db,dc,dd,de,df1,df2,df3,dh
    real(dp) :: ds,dt,dv,dw,dx1,dx2
    real(dp) :: dx3,dxm,dxp
    integer :: i,ind,iq,k,nm,np,nq
    ! Arrays
    real(dp) :: dq(5),df(101),dy(101)

    ! *********************************************************************

    dh=0.2D+0
    d2h=0.4D+0
    nm=101
    ind=0
    dq=(/1.D+0,5.656854228D+0,1.558845727D+1,3.2D+1,5.590169945D+0/)

    djp32=0.D+0
    dxm=gamma-1.5D+1
    if (xcut.GT.dxm) then
      if (ind.EQ.0) then
        do i=1,nm
          dy(i)=-1.5D+1+(i-1)*dh
          df(i)=1.D+0+DEXP(dy(i))
        end do
        ind=1
      end if
      dxp=gamma+5.D+0
      if (xcut.LT.dxp) then
        dc=dxp
      else
        dc=xcut
      end if
      db=DEXP(gamma-dc)
      dt=db
      do iq=1,5
        dd=iq*dc
        ds=DSQRT(dd)
        dw=1.+.3275911*ds
        dw=1.D+0/dw
        dv=dw*(.2258368458D+0+&
        & dw*(-.2521286676D+0+dw*(1.2596951294D+0+&
        & dw*(-1.2878224530D+0+dw*(.9406460699D+0)))))
        dv=dv+ds
        dv=1.5D+0*dv+ds*dd
        de=dt*dv/dq(iq)
        djp32=djp32+de
        if (DABS(de).LT.(1.D-07*djp32)) exit
        dt=-dt*db
      end do
      if (xcut.GE.dxp) return
      np=(dxp-xcut)/dh
      np=2*(np/2)
      np=nm-np
      nq=(15.-gamma)/dh
      nq=1+2*(nq/2)
      if (np.LT.nq) np=nq
      if (np.LE.nm) then
        df3=0.D+0
        dt=dy(np)+gamma
        dv=(dt-xcut)/2.D+0
        df3=DSQRT(dt)
        dx3=dt
        if (DABS(dv).GE.1.D-13) then
          df1=DSQRT(xcut)
          dt=df1+df3
          dw=(dv+dv)/(dt*dt)
          df2=dw*dw
          df2=df2+df2
          db=df2*(df2+7.D+0)+7.D+1
          dc=7.D+0*(1.D+1-df2)
          dc=dc*dw
          dd=-df2*(df2-2.8D+1)+1.4D+2
          dd=dd+dd
          ds=dt*((db-dc)*xcut/(1.D+0+DEXP(xcut-gamma))+dd*(xcut+dv)&
          & /(1.D+0+DEXP(xcut+dv-gamma))+(db+dc)*(dy(np)+gamma)/df(np))
          ds=ds*dv/4.2D+2
          djp32=djp32+ds
        end if
        if (np.NE.nm) then
          ds=0.D+0
          np=np+2
          do k=np,nm,2
            dx1=dx3
            df1=df3
            dx2=dy(k-1)+gamma
            dx3=dy(k)+gamma
            df3=DSQRT(dx3)
            dt=df1+df3
            dw=d2h/(dt*dt)
            df2=dw*dw
            df2=df2+df2
            db=df2*(df2+7.D+0)+7.D+1
            dc=7.D+0*(1.D+1-df2)
            dc=dc*dw
            dd=-df2*(df2-2.8D+1)+1.4D+2
            dd=dd+dd
            ds=ds+dt*((db-dc)*dx1/df(k-2)+dd*dx2/df(k-1)&
            & +(db+dc)*dx3/df(k))
          end do
          ds=ds*dh/4.2D+2
          djp32=djp32+ds
        end if
        if (xcut.GE.dxm) return
      end if
    end if
    djp32=dip32(gamma)-xcut*xcut*DSQRT(xcut)/2.5D+0
  end function djp32
  !!***

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
  !! PARENTS
  !!
  !! CHILDREN
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
  !! PARENTS
  !!
  !! CHILDREN
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
