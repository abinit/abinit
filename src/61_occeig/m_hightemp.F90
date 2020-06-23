!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_hightemp
!! NAME
!!  m_hightemp
!!
!! FUNCTION
!!  This module provides routines to run computations at very high temperatures
!!  with orbitals.
!!
!! COPYRIGHT
!!  Copyright (C) 2018-2019 ABINIT group (A. Blanchet)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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

module m_hightemp
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
  use m_kg,             only : mkkin
  use m_mpinfo,         only : ptabs_fourdp, proc_distrb_cycle
  use m_numeric_tools,  only : simpson, simpson_int
  use m_spacepar,       only : meanvalue_g

  implicit none

  type,public :: hightemp_type
    integer :: bcut,nbcut,version
    real(dp) :: ebcut,edc_kin_freeel,e_kin_freeel,e_ent_freeel
    real(dp) :: std_init,nfreeel,e_shiftfactor,ucvol
    logical :: prt_cg
  contains
    procedure :: compute_e_ent_freeel,compute_e_kin_freeel_approx
    procedure :: compute_nfreeel,compute_pw_avg_std
    procedure :: compute_e_shiftfactor,init,destroy
  end type hightemp_type

  ! type(hightemp_type),save,pointer :: hightemp=>null()
  public :: dip12,djp12,dip32,djp32,hightemp_e_heg
  public :: hightemp_gaussian_jintegral,hightemp_gaussian_kintegral
  public :: hightemp_dosfreeel,hightemp_get_e_shiftfactor
  public :: hightemp_get_nfreeel_approx,hightemp_prt_eigocc
contains

  !!****f* ABINIT/m_hightemp/init
  !! NAME
  !! init
  !!
  !! FUNCTION
  !! Initialize hightemp_type object, memory allocation of arrays...
  !!
  !! INPUTS
  !! this=hightemp_type object concerned
  !! activate=condition to activate or not the hightemp methods
  !! mband=maximum number of bands
  !! nbcut=input value of on how many states would you like to average the band energy
  !! rprimd(3,3)=dimensional primitive translations in real space (bohr)
  !!
  !! OUTPUT
  !! this=hightemp_type object concerned
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine init(this,mband,nbcut,prt_cg,rprimd,version)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    integer,intent(in) :: mband,nbcut,version,prt_cg
    ! Arrays
    real(dp),intent(in) :: rprimd(3,3)

    ! Local variables -------------------------
    ! Arrays
    real(dp) :: gprimd(3,3),rmet(3,3), gmet(3,3)

    ! *********************************************************************

    this%bcut=mband
    this%nbcut=nbcut
    this%version=version
    this%ebcut=zero
    this%edc_kin_freeel=zero
    this%e_kin_freeel=zero
    this%e_ent_freeel=zero
    this%std_init=zero
    this%nfreeel=zero
    this%e_shiftfactor=zero
    if(prt_cg==1) then
      this%prt_cg=.true.
    else
      this%prt_cg=.false.
    end if
    call metric(gmet,gprimd,-1,rmet,rprimd,this%ucvol)
  end subroutine init

  subroutine destroy(this)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this

    ! *********************************************************************

    this%bcut=0
    this%nbcut=0
    this%version=1
    this%ebcut=zero
    this%edc_kin_freeel=zero
    this%e_kin_freeel=zero
    this%e_ent_freeel=zero
    this%std_init=zero
    this%nfreeel=zero
    this%e_shiftfactor=zero
    this%prt_cg=.false.
    this%ucvol=zero
  end subroutine destroy

  !!****f* ABINIT/m_hightemp/compute_e_shiftfactor
  !! NAME
  !! compute_e_shiftfactor
  !!
  !! FUNCTION
  !! Compute the shift factor $U_0$ that appears in the density of states of free electrons.
  !!
  !! INPUTS
  !! this=hightemp_type object concerned
  !! eigen(mband*nkpt*nsppol)=eigenvalues (hartree)
  !! eknk(mband*nkpt*nsppol)=kinetic energies (hartree)
  !! mband=maximum number of bands
  !! nkpt=number of k points
  !! nsppol=1 for unpolarized, 2 for spin-polarized
  !!
  !! OUTPUT
  !! this=hightemp_type object concerned
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine compute_e_shiftfactor(this,eigen,eknk,mband,mpi_enreg,nband,nkpt,nsppol,wtk)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    integer,intent(in) :: mband,nkpt,nsppol
    type(MPI_type), intent(inout) :: mpi_enreg
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
    real(dp),intent(in) :: eknk(mband*nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)

    ! Local variables -------------------------
    ! Scalars
    integer :: band_index,ii,ikpt,isppol,krow,nband_k,niter
    ! Arrays
    real(dp) :: eigentemp(mband*nkpt*nsppol)
    real(dp) :: eknktemp(mband*nkpt*nsppol)
    logical :: mk(mband*nkpt*nsppol)

    ! *********************************************************************


    ! U_{K_0}
    ! this%e_shiftfactor=zero
    ! band_index=0
    ! do isppol=1,nsppol
    !   do ikpt=1,nkpt
    !     nband_k=nband(ikpt+(isppol-1)*nkpt)
    !     do ii=nband_k-this%nbcut+1,nband_k
    !       this%e_shiftfactor=this%e_shiftfactor+wtk(ikpt)*(eigen(band_index+ii)-eknk(band_index+ii))
    !     end do
    !     band_index=band_index+nband_k
    !   end do
    ! end do
    ! this%e_shiftfactor=this%e_shiftfactor/this%nbcut
    ! write(0,*) 'eknk_new', this%e_shiftfactor

    ! U_{HEG_0}
    if(this%version==1) then
      this%e_shiftfactor=zero
      band_index=0
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          do ii=nband_k-this%nbcut+1,nband_k
            this%e_shiftfactor=this%e_shiftfactor+wtk(ikpt)*(eigen(band_index+ii)-hightemp_e_heg(dble(ii),this%ucvol))
          end do
          band_index=band_index+nband_k
        end do
      end do
      this%e_shiftfactor=this%e_shiftfactor/this%nbcut
      this%ebcut=hightemp_e_heg(dble(this%bcut),this%ucvol)+this%e_shiftfactor
      ! write(0,*) this%ebcut, eigen(this%bcut*nkpt*nsppol)
      ! write(0,*) 'eheg_new', this%e_shiftfactor
    end if

    ! U_{LEGACY_0}
    if(this%version==2) then
      eigentemp(:)=zero
      eknktemp(:)=zero
      mk(:)=.true.
      ! Sorting eigen and eknk in ascending energy order.
      do ii=1,this%bcut*nkpt*nsppol
        krow=minloc(eigen,dim=1,mask=mk)
        mk(minloc(eigen,dim=1,mask=mk))=.false.
        eigentemp(ii)=eigen(krow)
        eknktemp(ii)=eknk(krow)
      end do
      ! Doing the average over the this%nbcut lasts states...
      niter=0
      this%e_shiftfactor=zero
      do ii=this%bcut*nkpt*nsppol-this%nbcut+1,this%bcut*nkpt*nsppol
        this%e_shiftfactor=this%e_shiftfactor+(eigentemp(ii)-eknktemp(ii))
        niter=niter+1
      end do
      this%e_shiftfactor=this%e_shiftfactor/niter
      this%ebcut=eigentemp(this%bcut*nkpt*nsppol)
      ! write(0,*) 'eknk_legacy', this%e_shiftfactor
    end if
  end subroutine compute_e_shiftfactor

  !!****f* ABINIT/m_hightemp/compute_nfreeel
  !! NAME
  !! compute_nfreeel
  !!
  !! FUNCTION
  !! Compute the value of the integral corresponding to the missing free electrons after energy cut.
  !! $I = \int_{Ec}^{\Infty}f(\epsilon)\frac{\sqrt{2}\Omega}{\pi^2}\sqrt{\epsilon - U_0}d \epsilon$
  !!
  !! INPUTS
  !! this=hightemp_type object concerned
  !! fermie=fermi energy (Hartree)
  !! mrgrid=number of grid points to compute the integral
  !! tsmear=smearing width (or temperature)
  !!
  !! OUTPUT
  !! this=hightemp_type object concerned
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine compute_nfreeel(this,fermie,tsmear)

    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: fermie,tsmear
    class(hightemp_type),intent(inout) :: this

    ! *********************************************************************

    call hightemp_get_nfreeel_approx(this%e_shiftfactor,this%ebcut,&
    & fermie,this%bcut,this%nfreeel,tsmear,this%ucvol,this%version)
    ! write(0,*) "nfreeel=",this%nfreeel
    ! write(0,*) this%nfreeel, factor*djp12(xcut,gamma), abs(factor*djp12(xcut,gamma)-this%nfreeel)
  end subroutine compute_nfreeel

  !!****f* ABINIT/m_hightemp/compute_e_kin_freeel
  !! NAME
  !! compute_e_kin_freeel
  !!
  !! FUNCTION
  !! Compute the value of the integral corresponding to the residual of energy after the band cut
  !! $I = \Omega\int_{Ec}^{\Infty}f(\epsilon)\frac{\sqrt{2}}{\pi^2}(\epsilon - U_0)^{3/2}d \epsilon$
  !!
  !! INPUTS
  !! this=hightemp_type object concerned
  !! fermie=fermi energy (Hartree)
  !! mrgrid=number of grid points to compute the integral
  !! tsmear=smearing width (or temperature)
  !!
  !! OUTPUT
  !! this=hightemp_type object concerned
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  ! subroutine compute_e_kin_freeel(this,fermie,mrgrid,nfftf,nspden,tsmear,vtrial)
  !
  !   ! Arguments -------------------------------
  !   ! Scalars
  !   class(hightemp_type),intent(inout) :: this
  !   integer,intent(in) :: mrgrid,nfftf,nspden
  !   real(dp),intent(in) :: fermie,tsmear
  !   ! Arrays
  !   real(dp),intent(in) :: vtrial(nfftf,nspden)
  !
  !   ! Local variables -------------------------
  !   ! Scalars
  !   integer :: ii,ifft,ispden
  !   real(dp) :: ix,step
  !   ! Arrays
  !   real(dp),dimension(:),allocatable :: values
  !
  !   ! *********************************************************************
  !   step=1e-5
  !
  !   ! Dynamic array find size
  !   ix=this%ebcut
  !   ii=0
  !   do while(fermi_dirac(ix,fermie,tsmear)>tol16)
  !     ii=ii+1
  !     ix=ix+step
  !   end do
  !   ABI_ALLOCATE(values,(ii))
  !
  !   ! Computation of e_kin_freeel
  !   ix=this%ebcut
  !   ii=0
  !   do while(fermi_dirac(ix,fermie,tsmear)>tol16)
  !     ii=ii+1
  !     values(ii)=fermi_dirac(ix,fermie,tsmear)*&
  !     & hightemp_dosfreeel(ix,this%e_shiftfactor,this%ucvol)*(ix-this%e_shiftfactor)
  !     ix=ix+step
  !   end do
  !   if (ii>1) then
  !     this%e_kin_freeel=simpson(step,values)
  !   end if
  !
  !   ! Computation of edc_kin_freeel
  !   this%edc_kin_freeel=zero
  !   do ispden=1,nspden
  !     do ifft=1,nfftf
  !       this%edc_kin_freeel=this%edc_kin_freeel+vtrial(ifft,ispden)
  !     end do
  !   end do
  !   ! Verifier la constante (/nspden**2)
  !   this%edc_kin_freeel=this%edc_kin_freeel*this%nfreeel/nspden/nfftf/nspden
  ! end subroutine compute_e_kin_freeel

  !!****f* ABINIT/m_hightemp/compute_e_kin_freeel_approx
  !! NAME
  !! compute_e_kin_freeel_approx
  !!
  !! FUNCTION
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
  subroutine compute_e_kin_freeel_approx(this,fermie,nfftf,nspden,tsmear,vtrial)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    real(dp),intent(in) :: fermie,tsmear
    integer,intent(in) :: nfftf,nspden
    ! Arrays
    real(dp),intent(in) :: vtrial(nfftf,nspden)

    ! Local variables -------------------------
    ! Scalars
    real(dp) :: factor,xcut,gamma
    integer :: ispden,ifft

    ! *********************************************************************

    factor=sqrt(2.)/(PI*PI)*this%ucvol*tsmear**(2.5)
    if(this%version==2) then
      xcut=(this%ebcut-this%e_shiftfactor)/tsmear
    else if(this%version==1) then
      xcut=hightemp_e_heg(dble(this%bcut),this%ucvol)/tsmear
    end if
    gamma=(fermie-this%e_shiftfactor)/tsmear

    this%e_kin_freeel=factor*djp32(xcut,gamma)
    ! write(0,*) "e_kin_freeel=",this%e_kin_freeel

    ! Computation of edc_kin_freeel
    this%edc_kin_freeel=zero
    do ispden=1,nspden
      do ifft=1,nfftf
        this%edc_kin_freeel=this%edc_kin_freeel+vtrial(ifft,ispden)
      end do
    end do
    ! Verifier la constante (/nspden**2)
    this%edc_kin_freeel=this%edc_kin_freeel*this%nfreeel/nspden/nfftf/nspden
    if(this%e_kin_freeel<tol16) this%e_kin_freeel=zero
  end subroutine compute_e_kin_freeel_approx

  !!****f* ABINIT/m_hightemp/compute_e_ent_freeel
  !! NAME
  !! compute_e_ent_freeel
  !!
  !! FUNCTION
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
  subroutine compute_e_ent_freeel(this,fermie,tsmear)
    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    real(dp),intent(in) :: fermie,tsmear

    ! Local variables -------------------------
    ! Scalars
    integer :: ii
    real(dp) :: ix,step
    ! Arrays
    real(dp),dimension(:),allocatable :: valuesent

    ! *********************************************************************
    step=1e-4
    this%e_ent_freeel=zero

    ! Dynamic array find size
    ix=this%ebcut
    ii=0
    do while(fermi_dirac(ix,fermie,tsmear)>tol14)
      ii=ii+1
      ix=ix+step
    end do

    ABI_ALLOCATE(valuesent,(ii))

    ix=this%ebcut
    ii=0
    do while(fermi_dirac(ix,fermie,tsmear)>tol14)
      ii=ii+1
      valuesent(ii)=(fermi_dirac(ix,fermie,tsmear)*log(fermi_dirac(ix,fermie,tsmear))+&
      & (1.-fermi_dirac(ix,fermie,tsmear))*log(1.-fermi_dirac(ix,fermie,tsmear)))*&
      & hightemp_dosfreeel(ix,this%e_shiftfactor,this%ucvol)
      ix=ix+step
    end do
    if (ii>1) then
      this%e_ent_freeel=simpson(step,valuesent)
    end if
  end subroutine compute_e_ent_freeel

  !****f* ABINIT/m_hightemp/compute_pw_avg_std
  ! NAME
  ! compute_pw_avg_std
  !
  ! FUNCTION
  ! Compute the planewaves average standard deviation over lasts bands.
  !
  ! INPUTS
  ! this=hightemp_type object concerned
  ! eigen(mband*nkpt*nsppol)=eigenvalues (hartree)
  ! eknk(mband*nkpt*nsppol)=kinetic energies (hartree)
  ! mband=maximum number of bands
  ! nkpt=number of k points
  ! nsppol=1 for unpolarized, 2 for spin-polarized
  !
  ! OUTPUT
  ! this=hightemp_type object concerned
  !
  ! PARENTS
  !
  ! CHILDREN
  !
  ! SOURCE
  subroutine compute_pw_avg_std(this,cg,ckpt,ecut,eig_k,ek_k,exchn2n3d,fnameabo,istwfk,kg_k,kpt,&
  & mcg,mpi_enreg,mpw,nband,nkpt,npw_k,nsppol,rprimd,wtk)

    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: ckpt,mcg,mpw,nkpt,npw_k,nsppol,exchn2n3d
    real(dp),intent(in) :: ecut,wtk
    class(hightemp_type),intent(inout) :: this
    type(MPI_type),intent(inout) :: mpi_enreg
    ! Arrays
    integer,intent(in) :: istwfk(nkpt),kg_k(3,npw_k),nband(nkpt)
    real(dp),intent(in) :: kpt(3,nkpt),rprimd(3,3)
    real(dp),intent(in) :: cg(2,mcg)
    real(dp),intent(in) :: eig_k(nband(ckpt)),ek_k(nband(ckpt))
    character(len=*),intent(in) :: fnameabo

    ! Local variables -------------------------
    ! Scalars
    integer :: cgshift,cgshift_tot,iband,iout,ipw,ipw_tot
    integer :: mcg_tot,mpierr,npw_tot
    real(dp) :: average,ucvol,kpg1_tot,kpg2_tot,kpg3_tot,tmp_std1,tmp_std2
    character(len=50) :: filenameoutpw
    character(len=4) :: mode_paral
    ! Arrays
    integer :: nppw_tot(mpi_enreg%nproc_band)
    real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
    integer,allocatable :: kg_tot(:,:)
    real(dp),allocatable :: cg_tot(:,:),kpgnorm(:),kpgnorm_tot(:),tempwfk(:,:)

    ! *********************************************************************

    write(std_out,'(a,I5.5,a)') 'Writing plane waves coefs of kpt=',ckpt,'...'
    mode_paral='PERS'
    iout=-1
    call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

    ABI_ALLOCATE(kpgnorm,(npw_k))
    call getkpgnorm(gprimd,kpt(:,ckpt),kg_k,kpgnorm,npw_k)

    mcg_tot=mcg
    nppw_tot(:)=0
    nppw_tot(mpi_enreg%me_band + 1)=npw_k
    call xmpi_sum(nppw_tot,mpi_enreg%comm_world,mpierr)
    call xmpi_sum(mcg_tot,mpi_enreg%comm_world,mpierr)

    npw_tot=sum(nppw_tot)
    ABI_ALLOCATE(cg_tot,(2,mcg_tot))
    ABI_ALLOCATE(kpgnorm_tot,(npw_tot))
    ABI_ALLOCATE(kg_tot,(3,npw_tot))

    cg_tot(:,:)=zero
    kpgnorm_tot(:)=zero
    kg_tot(:,:)=0
    do ipw=1,npw_k
      ipw_tot=sum(nppw_tot(1:mpi_enreg%me_band))+ipw
      kpgnorm_tot(ipw_tot)=kpgnorm(ipw)
      kg_tot(:,ipw_tot)=kg_k(:,ipw)

      do iband=1,nband(ckpt)
        cgshift=(iband-1)*npw_k
        cgshift_tot=(iband-1)*npw_tot
        cg_tot(:,cgshift_tot+ipw_tot)=cg(:,cgshift+ipw)
      end do
    end do

    call xmpi_sum(cg_tot,mpi_enreg%comm_world,mpierr)
    call xmpi_sum(kpgnorm_tot,mpi_enreg%comm_world,mpierr)
    call xmpi_sum(kg_tot,mpi_enreg%comm_world,mpierr)

    if(mpi_enreg%me==mpi_enreg%me_kpt) then
      if(this%prt_cg) then
        ! Writting plane waves vector coordinates
        write(filenameoutpw,'(A,I5.5)') '_PW_MESH_k',ckpt

        open(file=trim(fnameabo)//trim(filenameoutpw), unit=23)
        do ipw=1, npw_tot
          kpg1_tot=kpt(1,ckpt)+dble(kg_tot(1,ipw))
          kpg2_tot=kpt(2,ckpt)+dble(kg_tot(2,ipw))
          kpg3_tot=kpt(3,ckpt)+dble(kg_tot(3,ipw))

          write(23,'(i14,ES13.5,ES13.5,ES13.5,ES13.5)')&
          & ipw,&
          & gprimd(1,1)*kpg1_tot+gprimd(1,2)*kpg2_tot+gprimd(1,3)*kpg3_tot,&
          & gprimd(2,1)*kpg1_tot+gprimd(2,2)*kpg2_tot+gprimd(2,3)*kpg3_tot,&
          & gprimd(3,1)*kpg1_tot+gprimd(3,2)*kpg2_tot+gprimd(3,3)*kpg3_tot,&
          & kpgnorm_tot(ipw)
        end do
        close(23)

        ! Writting Eigen energies
        write(filenameoutpw,'(A,I5.5)') '_PW_EIG_k',ckpt
        open(file=trim(fnameabo)//trim(filenameoutpw), unit=23)
        write(23,'(ES12.5)') eig_k
        close(23)

        ! Writting Kinetic energies
        write(filenameoutpw,'(A,I5.5)') '_PW_KIN_k',ckpt
        open(file=trim(fnameabo)//trim(filenameoutpw), unit=23)
        write(23,'(ES12.5)') ek_k
        close(23)
      end if

      ABI_ALLOCATE(tempwfk,(npw_tot, 3))
      tmp_std1=zero
      do iband=1, nband(ckpt)
        tempwfk(:,:)=zero
        cgshift_tot=(iband-1)*npw_tot
        do ipw=1, npw_tot
          tempwfk(ipw,2)=ABS(dcmplx(cg_tot(1,cgshift_tot+ipw), cg_tot(2,cgshift_tot+ipw)))
          ! Hartree to eV * Kinetic energy of the pw
          tempwfk(ipw,1)=(2*pi*kpgnorm_tot(ipw))**2/2.
          tempwfk(ipw,3)=ipw
        end do

        ! Printing Plane wave in _PW output file on request
        if(this%prt_cg) then
          write(filenameoutpw, '(A,I5.5,A,I5.5)') '_PW_k',ckpt,'_b',iband
          open(file=trim(fnameabo)//trim(filenameoutpw), unit=23)
          do ipw=1, npw_tot
            write(23,'(i14,ES14.6,ES14.6,i14)') int(tempwfk(ipw, 3)),27.2114*tempwfk(ipw, 1),tempwfk(ipw,2)
          end do
          close(23)
        end if


        ! Computing standard deviation for lasts bands
        write(0,*) ckpt,sum(tempwfk(:,2)**2)
        if(iband>=nband(ckpt)-this%nbcut+1) then
          tmp_std2=zero
          average=zero
          cgshift_tot=(iband-1)*npw_tot
          do ipw=1, npw_tot
            average=average+tempwfk(ipw,1)*tempwfk(ipw,2)**2
          end do
          do ipw=1, npw_tot
            tmp_std2=tmp_std2+(tempwfk(ipw,1)-average)**2*tempwfk(ipw,2)**2
          end do
          tmp_std1=tmp_std1+wtk*0.5*sqrt(tmp_std2)/this%nbcut
        end if
      end do
      ! Add kpt weighted contribution to standard deviation
      this%std_init=this%std_init+tmp_std1
      write(0,*) 'tmp = ',ckpt,tmp_std1
      ABI_DEALLOCATE(tempwfk)
    end if
    ABI_DEALLOCATE(kpgnorm)
    ABI_DEALLOCATE(kg_tot)
    ABI_DEALLOCATE(cg_tot)
    ABI_DEALLOCATE(kpgnorm_tot)
  end subroutine compute_pw_avg_std

  ! *********************************************************************

  !!****f* ABINIT/m_hightemp/dip12
  !! NAME
  !! dip12
  !!
  !! FUNCTION
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

  !!****f* ABINIT/m_hightemp/djp12
  !! NAME
  !! djp12
  !!
  !! FUNCTION
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

  !!****f* ABINIT/m_hightemp/dip32
  !! NAME
  !! dip32
  !!
  !! FUNCTION
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

  !!****f* ABINIT/m_hightemp/djp32
  !! NAME
  !! djp32
  !!
  !! FUNCTION
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

  !!****f* ABINIT/m_hightemp/hightemp_dosfreeel
  !! NAME
  !! hightemp_dosfreeel
  !!
  !! FUNCTION
  !! Returns the free particle density of states for a given energy (in Hartree)
  !!
  !! INPUTS
  !! energy=get the value of the free particle density of states at this energy
  !! e_shiftfactor=energy shift factor
  !! ucvol=unit cell volume (bohr^3)
  !!
  !! OUTPUT
  !! freedos=value of free particle density of states at given energy
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function hightemp_dosfreeel(energy,e_shiftfactor,ucvol)

    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: energy,e_shiftfactor,ucvol
    real(dp) :: hightemp_dosfreeel

    ! *********************************************************************

    hightemp_dosfreeel=sqrt(2.)*ucvol*sqrt(energy-e_shiftfactor)/(PI*PI)
  end function hightemp_dosfreeel

  !!****f* ABINIT/m_hightemp/hightemp_e_heg
  !! NAME
  !! hightemp_e_heg
  !!
  !! FUNCTION
  !! Returns the free particle energy for a given orbital (in Hartree)
  !!
  !! INPUTS
  !! iband=considered orbital
  !! ucvol=unit cell volume (bohr^3)
  !!
  !! OUTPUT
  !! hightemp_e_heg=value of free particle density of states at given energy
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function hightemp_e_heg(iband,ucvol)

    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: iband,ucvol
    real(dp) :: hightemp_e_heg

    ! *********************************************************************
    hightemp_e_heg=.5*(iband*6*PI*PI/ucvol)**(2./3.)
  end function hightemp_e_heg

  !!****f* ABINIT/m_hightemp/hightemp_gaussian_jintegral
  !! NAME
  !! hightemp_dosfreeel
  !!
  !! FUNCTION
  !! Returns the value of the first integral J(\sigma,\rho_{0_n}) =
  !! \int_0^\infty \rho^2 e^{-(\rho_{0_n}- \rho)^2/(2\sigma^2)}d\rho
  !!
  !! INPUTS
  !! sigma=standard deviation
  !! alpha=center of the gaussian
  !!
  !! OUTPUT
  !! hightemp_gaussian_jintegral=value of the integral
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function hightemp_gaussian_jintegral(sigma,alpha)

    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: sigma,alpha
    real(dp) :: hightemp_gaussian_jintegral

    ! *********************************************************************

    hightemp_gaussian_jintegral=alpha*sigma**2*exp(-alpha**2/(2*sigma**2))+.5*sigma*sqrt(2*PI)&
    & *(alpha**2+sigma**2)*(1+erf(alpha/(sigma*sqrt(2.))))
  end function hightemp_gaussian_jintegral

  !!****f* ABINIT/m_hightemp/hightemp_gaussian_kintegral
  !! NAME
  !! hightemp_dosfreeel
  !!
  !! FUNCTION
  !! Returns the value of the first integral J(\sigma,\rho_{0_n}) =
  !! \int_0^\infty \rho^4 e^{-(\rho_{0_n}- \rho)^2/(2\sigma^2)}d\rho
  !!
  !! INPUTS
  !! sigma=standard deviation
  !! alpha=center of the gaussian
  !!
  !! OUTPUT
  !! hightemp_gaussian_kintegral=value of the integral
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function hightemp_gaussian_kintegral(sigma,alpha)

    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: sigma,alpha
    real(dp) :: hightemp_gaussian_kintegral

    ! *********************************************************************

    hightemp_gaussian_kintegral=sigma**2*(alpha**3+5*alpha*sigma**2)*exp(-alpha**2/(2*sigma**2))&
    & +.5*sigma*sqrt(2*PI)*(alpha**4+6*alpha**2*sigma**2+3*sigma**4)*(1+erf(alpha/(sigma*sqrt(2.))))
  end function hightemp_gaussian_kintegral

  !!****f* ABINIT/m_hightemp/hightemp_get_e_shiftfactor
  !! NAME
  !! hightemp_get_e_shiftfactor
  !!
  !! FUNCTION
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
  subroutine hightemp_get_e_shiftfactor(cg,ecut,ecutsm,effmass_free,eigen,gmet,hightemp,&
  & istwfk,kg,kptns,mband,mcg,mkmem,mpi_enreg,mpw,my_nspinor,nband,nkpt,nsppol,npwarr,wtk)
    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: mband,mcg,mpw,mkmem,my_nspinor,nkpt,nsppol
    real(dp),intent(in) :: ecut,ecutsm,effmass_free
    type(MPI_type),intent(inout) :: mpi_enreg
    type(hightemp_type),pointer,intent(inout) :: hightemp
    ! Arrays
    integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt)
    integer,intent(in) :: kg(3,mpw*mkmem)
    real(dp),intent(in) :: cg(2,mcg),eigen(mband*nkpt*nsppol)
    real(dp),intent(in) :: kptns(3,nkpt),gmet(3,3)
    real(dp),intent(in) :: wtk(nkpt)

    ! Local variables -------------------------
    ! Scalars
    integer :: bdtot_index,blocksize,iband,iblock,iblocksize
    integer :: ikpt,isppol,mpierr,nband_k,nblockbd,npw_k
    integer :: ii
    real(dp) :: ar
    ! Arrays
    integer,allocatable :: kg_k(:,:)
    real(dp),allocatable :: eknk(:),ek_k(:),kinpw(:)

    ! *********************************************************************

    bdtot_index=0

    ABI_ALLOCATE(eknk,(mband*nkpt*nsppol))
    do isppol=1,nsppol
      do ikpt=1,nkpt
        npw_k=npwarr(ikpt)
        nband_k=nband(ikpt+(isppol-1)*nkpt)

        ABI_ALLOCATE(kg_k,(3,npw_k))
        ABI_ALLOCATE(kinpw,(npw_k))
        ABI_ALLOCATE(ek_k,(nband_k))

        kg_k(:,1:npw_k)=kg(:,1:npw_k)
        nblockbd=nband_k/(mpi_enreg%nproc_band*mpi_enreg%bandpp)
        blocksize=nband_k/nblockbd

        call mkkin(ecut,ecutsm,effmass_free,gmet,kg_k,kinpw,kptns(:,ikpt),npw_k,0,0)

        do iblock=1,nblockbd
          do iblocksize=1,blocksize
            iband=(iblock-1)*blocksize+iblocksize
            call meanvalue_g(ar,kinpw,0,istwfk(ikpt),mpi_enreg,npw_k,my_nspinor,&
            & cg(:,1+(iband-1)*npw_k*my_nspinor:iband*npw_k*my_nspinor),&
            & cg(:,1+(iband-1)*npw_k*my_nspinor:iband*npw_k*my_nspinor),0)

            ek_k(iband)=ar
          end do
        end do

        eknk (1+bdtot_index : nband_k+bdtot_index)=ek_k (:)
        ABI_DEALLOCATE(ek_k)
        ABI_DEALLOCATE(kinpw)
        ABI_DEALLOCATE(kg_k)

        bdtot_index=bdtot_index+nband_k
      end do
    end do

    call hightemp%compute_e_shiftfactor(eigen,eknk,mband,mpi_enreg,&
    & nband,nkpt,nsppol,wtk)

    ABI_DEALLOCATE(eknk)
  end subroutine hightemp_get_e_shiftfactor

  !!****f* ABINIT/m_hightemp/hightemp_getnfreeel
  !! NAME
  !! hightemp_getnfreeel
  !!
  !! FUNCTION
  !! Compute the value of the integral corresponding to the missing free electrons after energy cut.
  !! $I = \int_{Ec}^{\Infty}f(\epsilon)\frac{\sqrt{2}\Omega}{\pi^2}\sqrt{\epsilon - U_0}d \epsilon$
  !!
  !! INPUTS
  !! this=hightemp_type object concerned
  !! fermie=fermi energy (Hartree)
  !! mrgrid=number of grid points to compute the integral
  !! tsmear=smearing width (or temperature)
  !! e_shiftfactor=energy shift factor
  !! ucvol=unit cell volume (bohr^3)
  !!
  !! OUTPUT
  !! nfreeel=number of free electrons after the energy cut with given fermi level
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  ! subroutine hightemp_getnfreeel(ebcut,entropy,fermie,mrgrid,nfreeel,tsmear,e_shiftfactor,ucvol)
  !
  !   ! Arguments -------------------------------
  !   ! Scalars
  !   integer,intent(in) :: mrgrid
  !   real(dp),intent(in) :: ebcut,fermie,tsmear,e_shiftfactor,ucvol
  !   real(dp),intent(out) :: entropy,nfreeel
  !
  !   ! Local variables -------------------------
  !   ! Scalars
  !   integer :: ii
  !   real(dp) :: ix,step
  !   ! Arrays
  !   real(dp),dimension(:),allocatable :: valuesnel,valuesent
  !
  !   ! *********************************************************************
  !   step=1e-5
  !   nfreeel=zero
  !   entropy=zero
  !
  !   ! Dynamic array find size
  !   ix=ebcut
  !   ii=0
  !   do while(fermi_dirac(ix,fermie,tsmear)>1e-14)
  !     ii=ii+1
  !     ix=ix+step
  !   end do
  !
  !   ABI_ALLOCATE(valuesnel,(ii))
  !   ABI_ALLOCATE(valuesent,(ii))
  !
  !   ix=ebcut
  !   ii=0
  !   do while(fermi_dirac(ix,fermie,tsmear)>1e-14)
  !     ii=ii+1
  !     valuesnel(ii)=fermi_dirac(ix,fermie,tsmear)*hightemp_dosfreeel(ix,e_shiftfactor,ucvol)
  !     valuesent(ii)=(fermi_dirac(ix,fermie,tsmear)*log(fermi_dirac(ix,fermie,tsmear))+&
  !     & (1.-fermi_dirac(ix,fermie,tsmear))*log(1.-fermi_dirac(ix,fermie,tsmear)))*&
  !     & hightemp_dosfreeel(ix,e_shiftfactor,ucvol)
  !     ix=ix+step
  !   end do
  !   if (ii>1) then
  !     nfreeel=simpson(step,valuesnel)
  !     entropy=simpson(step,valuesent)
  !   end if
  ! end subroutine hightemp_getnfreeel

  !!****f* ABINIT/m_hightemp/hightemp_get_nfreeel_approx
  !! NAME
  !! hightemp_get_nfreeel_approx
  !!
  !! FUNCTION
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
  subroutine hightemp_get_nfreeel_approx(e_shiftfactor,ebcut,fermie,bcut,&
  & nfreeel,tsmear,ucvol,version)

    ! Arguments -------------------------------
    ! Scalars
    integer :: bcut,version
    real(dp),intent(inout) :: nfreeel
    real(dp),intent(in) :: e_shiftfactor,ebcut,fermie,tsmear,ucvol

    ! Local variables -------------------------
    ! Scalars
    real(dp) :: factor,xcut,gamma

    ! *********************************************************************

    nfreeel=0.
    factor=sqrt(2.)/(PI*PI)*ucvol*tsmear**(1.5)
    if(version==0) then
      xcut=(ebcut-e_shiftfactor)/tsmear
    else if(version==1) then
      xcut=hightemp_e_heg(dble(bcut),ucvol)/tsmear
    end if
    gamma=(fermie-e_shiftfactor)/tsmear

    nfreeel=factor*djp12(xcut,gamma)
    if(nfreeel<tol16) nfreeel=zero
  end subroutine hightemp_get_nfreeel_approx

  !!****f* ABINIT/m_hightemp/hightemp_prt_eigocc
  !! NAME
  !! hightemp_prt_eigocc
  !!
  !! FUNCTION
  !! Printing in a _EIGOCC file many informations like Fermi energy, Unit cell volume,
  !! electronic temperature (tsmear when occopt=3), eigenvalues and occupation for each k-point,
  !! wheight of the k-point, number of bands, number of k-points and more...
  !! This file is intended to be used for custom DOS computation with external tools for example.
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
  subroutine hightemp_prt_eigocc(e_kin_freeel,e_shiftfactor,eigen,etotal,energies,fnameabo,&
  & iout,iter,kptns,mband,nband,nfreeel,nkpt,nsppol,occ,rprimd,tsmear,usepaw,wtk,&
  & strten,istep) ! Optional arguments

    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: iout,iter,mband,nkpt,nsppol,usepaw
    real(dp),intent(in) :: e_kin_freeel,e_shiftfactor,etotal,nfreeel,tsmear
    character(len=*),intent(in) :: fnameabo
    type(energies_type),intent(in) :: energies
    integer,intent(in),optional :: istep
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: eigen(mband*nkpt*nsppol),kptns(3,nkpt)
    real(dp),intent(in) :: rprimd(3,3)
    real(dp),intent(in) :: occ(mband*nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)
    real(dp),intent(in),optional :: strten(6)

    ! Local variables -------------------------
    ! Scalars
    integer :: band_index,iband,ii,ikpt,isppol,nband_k,temp_unit
    real(dp) :: ucvol,pressure
    character(len=200) :: fnameabo_eigocc
    character(len=39) :: kind_of_output
    character(len=500) :: msg
    ! Arrays
    real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)

    ! *********************************************************************

    band_index=0
    call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

    fnameabo_eigocc=trim(fnameabo) // trim('_EIGOCC')
    write(msg,'(a,a)') ' prt_eigocc : about to open file ',trim(fnameabo_eigocc)
    call wrtout(iout,msg,'COLL')
    if (open_file(fnameabo_eigocc,msg,newunit=temp_unit,status='unknown',form='formatted') /= 0) then
      MSG_ERROR(msg)
    end if
    rewind(temp_unit)

    ! write(msg, '(a,f12.6,a,f12.6,a,ES16.8,a)') &
    !   & ' Fermi (or HOMO) energy (hartree)= ',fermie,'   Average Vxc (hartree)=',vxcavg,&
    !   & ', Total energy= ',etotal, ' Ha'
    ! call wrtout(temp_unit,msg,'COLL')
    ! write(msg, '(a,e16.8,a,f16.8,a,a,ES12.4,a)') &
    !   & ' Unit cell volume ucvol= ',ucvol,' bohr^3, Electronic temperature= ',tsmear,' Ha',&
    !   & ', Pressure= ',-(strten(1)+strten(2)+strten(3))*HaBohr3_GPa/3.0_dp,' GPa'
    ! call wrtout(temp_unit,msg,'COLL')

    write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
      & ' Chemical potential = ',energies%e_fermie,' Ha         XC energy          = ',energies%e_xc,&
      & ' Ha         Total free energy  = ',etotal, ' Ha'
    call wrtout(temp_unit,msg,'COLL')

    write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
    & ' Kinetic energy     = ',energies%e_kinetic,' Ha         Kin. free el. E    = ',e_kin_freeel,&
    & ' Ha         Total kin. energy  = ',energies%e_kinetic+e_kin_freeel, ' Ha'
    call wrtout(temp_unit,msg,'COLL')

    write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
    & ' Hartree energy     = ',energies%e_hartree,' Ha         Ewald energy       = ',energies%e_ewald,&
    & ' Ha         E. Shift factor U0 = ',e_shiftfactor, ' Ha'
    call wrtout(temp_unit,msg,'COLL')

    if (usepaw==0) then
      write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
        & ' PsPCore energy     = ',energies%e_corepsp,' Ha         Loc. Psp. energy   = ',energies%e_localpsp,&
        & ' Ha         NL PsP energy      = ',energies%e_nlpsp_vfock, ' Ha'
    else if (usepaw==1) then
      write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
        & ' PsPCore energy     = ',energies%e_corepsp,' Ha         Loc. Psp. energy   = ',energies%e_localpsp,&
        & ' Ha         PAW PsP energy     = ',energies%e_paw, ' Ha'
    end if
    call wrtout(temp_unit,msg,'COLL')

    write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
      & ' Entropy energy     = ',energies%e_entropy,' Ha         Band energy        = ',energies%e_eigenvalues,&
      & ' Ha         Internal energy    = ',etotal-energies%e_entropy, ' Ha'
    call wrtout(temp_unit,msg,'COLL')
    if(present(strten)) then
      pressure=-(strten(1)+strten(2)+strten(3))*HaBohr3_GPa/3.0_dp
    else
      pressure=zero
    end if
    write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
    & ' Unit cell vol      = ',ucvol,' Bohr^3     Elec. temperature  = ',tsmear,&
    & ' Ha         Pressure           = ',pressure,' GPa'
    call wrtout(temp_unit,msg,'COLL')

    write(msg, '(a,i12,a,i12,a,ES15.8)') &
    & ' N. elec. iter      = ',istep,'            MD iteration       = ',iter,&
    & '            N. free electrons  = ',nfreeel
    call wrtout(temp_unit,msg,'COLL')

    ! Loop over spins
    do isppol=1,nsppol
      ! write(msg, '(a,i6,a)') ' Eigenvalues (hartree) for nkpt=',nkpt,'k points:'
      write(msg, '(a,i12,a,i12)') ' Number of kpts     = ',nkpt,'            Number of bands    = ',&
      & mband
      call wrtout(temp_unit,msg,'COLL')

      ! Loop over k-points
      do ikpt=1,nkpt
        nband_k=nband(ikpt+(isppol-1)*nkpt)
        write(msg, '(a,i6,a,i6,a,ES12.5,a,3f8.4,a)') &
          & ' ------- ikpt = ',ikpt,'      nband = ',nband_k,'     ikpt weight = ',wtk(ikpt)+tol10,&
          & '   Reduced coordinates : ',kptns(1:3,ikpt)+tol10,' -------'
        call wrtout(temp_unit,msg,'COLL')

        ! Loop over bands
        do ii=0,(nband_k-1)/4
          write(msg, '(4(i6,ES12.4,ES13.5,a,1x))') &
            & (iband,eigen(iband+band_index),occ(iband+band_index),',',iband=1+ii*4,min(nband_k,4+ii*4))
          call wrtout(temp_unit,msg,'COLL')
        end do

        band_index=band_index+nband_k
      end do ! do ikpt=1,nkpt
    end do ! do isppol=1,nsppol

    close(temp_unit)
  end subroutine hightemp_prt_eigocc

end module m_hightemp
