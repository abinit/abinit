
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_phdos

  use defs_basis
  use m_errors
  use m_abicore
  use m_phonons
  use m_ifc,              only : ifc_type
  use m_crystal,          only : crystal_t
  use m_ddb,              only : ddb_type
  use m_tdep_qpt,         only : Qpoints_type
  use m_tdep_readwrite,   only : Input_Variables_type
  use m_tdep_latt,        only : Lattice_Variables_type
  use m_tdep_sym,         only : Symetries_Variables_type
  use m_tdep_shell,       only : Shell_Variables_type
  use m_tdep_abitypes,    only : tdep_ifc2phij, tdep_read_ifc, tdep_write_ifc
  use m_xmpi

  implicit none

  public :: tdep_calc_phdos
  public :: tdep_calc_thermo
  public :: tdep_calc_elastic

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tdep_calc_phdos(Crystal,Ifc,InVar,Lattice,natom,natom_unitcell,Phij_NN,PHdos,Qpt,Rlatt4abi,Shell2at,Sym)

  integer :: prtdos,ii,iqpt,iatom
  integer :: natom,natom_unitcell,iomega
  integer :: dos_ngqpt(3)
  integer :: count_wminmax(2)
  character (len=25):: phdos_fname
  character(len=500) :: msg
  double precision :: dossmear,integ,domega
  double precision :: Phij_NN(3*natom,3*natom)
  double precision :: Rlatt4abi(3,natom_unitcell,natom)
  double precision :: dos_qshift(3)
  real(dp) :: wminmax(2)
  double precision, allocatable :: displ(:,:),omega(:,:)
  type(Input_Variables_type),intent(in) :: InVar
  type(phonon_dos_type),intent(out) :: PHdos
  type(ifc_type),intent(inout) :: Ifc
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Symetries_Variables_type),intent(in) :: Sym
  type(crystal_t),intent(in) :: Crystal
  type(Qpoints_type),intent(in) :: Qpt
  type(Shell_Variables_type),intent(in) :: Shell2at

  write(InVar%stdout,*)' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '################### vibrational Density OF States (vDOS) ####################'
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,'(a)') ' See the vdos.dat and TDEP_PHDOS* files'

! Copy Phij_NN to Ifc%atmfrc
! ==========================
  call tdep_ifc2phij(Ifc%dipdip,Ifc,InVar,Lattice,natom_unitcell,0,Phij_NN,Rlatt4abi,Shell2at,Sym)

! Write Ifc%atmfrc in the ifc.tdep file
! =====================================
  call tdep_write_ifc(Crystal,Ifc,InVar,natom_unitcell,0)

! Read the previous IFC from ifc.tdep input file and copy to Phij
! ===============================================================
  if (InVar%ReadIFC.eq.2) then
!   Read IFC from ifc.tdep (ReadIFC=2)
    call tdep_read_ifc(Ifc,InVar,natom_unitcell)
!   Copy Ifc%atmfrc to Phij_NN
    call tdep_ifc2phij(Ifc%dipdip,Ifc,InVar,Lattice,natom_unitcell,1,Phij_NN,Rlatt4abi,Shell2at,Sym)
!   Copy Phij_NN to Ifc%atmfrc
    call tdep_ifc2phij(Ifc%dipdip,Ifc,InVar,Lattice,natom_unitcell,0,Phij_NN,Rlatt4abi,Shell2at,Sym)
!   Write IFC in ifc.out (for check)
    call tdep_write_ifc(Crystal,Ifc,InVar,natom_unitcell,1)

!   Write the Phij_NN-new.dat file
    if (InVar%debug) then
      write(InVar%stdout,'(a)') ' See the Phij_NN-new.dat file corresponding to the ifc.tdep/Phij_NN file'
      open(unit=55,file=trim(InVar%output_prefix)//'Phij_NN-new.dat')
      do iatom=1,3*natom
        write(55,'(10000(f10.6,1x))') Phij_NN(iatom,:)
      end do
      close(55)
    end if
  end if

! Compute the DOS
! ===============
  prtdos=1 !Gaussian
!  prtdos=2 !Tetra
  dossmear=4.5d-6
!  dossmear=4.5d-5
  dos_qshift(:)=     zero
  dos_ngqpt(:)=InVar%ngqpt2(:)
  write(InVar%stdout,'(a)') ' Compute the vDOS'
  ! Only 1 shift in q-mesh
  wminmax = zero
  do
    call mkphdos(PHdos,Crystal,Ifc,prtdos,InVar%dosdeltae,dossmear,dos_ngqpt,1,dos_qshift, &
      "freq_displ", wminmax, count_wminmax, XMPI_WORLD)
     if (all(count_wminmax == 0)) exit
     wminmax(1) = wminmax(1) - abs(wminmax(1)) * 0.05
     wminmax(2) = wminmax(2) + abs(wminmax(2)) * 0.05
     call phdos%free()
     write(msg, "(a, 2f8.5)")"Initial frequency mesh not large enough. Recomputing PHDOS with wmin, wmax: ",wminmax
     call wrtout(std_out, msg)
  end do

  write(InVar%stdout,'(a)') ' ------- achieved'
  write(InVar%stdout,'(a)') ' (Please, pay attention to convergency wrt the BZ mesh : the ngqpt2 input variable)'

! Compute the frequencies
! =======================
  ABI_MALLOC(displ,(2*3*natom_unitcell*3*natom_unitcell,Qpt%nqpt)); displ(:,:)=zero
  ABI_MALLOC(omega,(3*natom_unitcell,Qpt%nqpt)); omega(:,:)=zero
  open(unit=53,file=trim(InVar%output_prefix)//'omega-abinit.dat')
  do iqpt=1,Qpt%nqpt
    call ifc%fourq(Crystal,Qpt%qpt_red(:,iqpt),omega(:,iqpt),displ(:,iqpt))
    if (iqpt.le.Qpt%nqpt) then
      if (InVar%Enunit.eq.0) write(53,'(i5,1x,100(f15.6,1x))') iqpt,(omega(ii,iqpt)*Ha_eV*1000,ii=1,3*natom_unitcell)
      if (InVar%Enunit.eq.1) write(53,'(i5,1x,100(f15.6,1x))') iqpt,(omega(ii,iqpt)*Ha_cmm1   ,ii=1,3*natom_unitcell)
      if (InVar%Enunit.eq.2) write(53,'(i5,1x,100(f15.6,1x))') iqpt,(omega(ii,iqpt)           ,ii=1,3*natom_unitcell)
    end if
  end do
  close(53)

! Print the DOS
! =============
  phdos_fname = trim(InVar%output_prefix)//"_PHDOS"
  call phdos%print(phdos_fname)
  domega=(InVar%dosdeltae*Ha_meV)
  integ=0.d0
  do iomega=1,PHdos%nomega
    integ=integ + domega*PHdos%phdos(iomega)
  end do
  PHdos%phdos(:)=PHdos%phdos(:)/integ
  open(unit=56,file=trim(InVar%output_prefix)//'vdos.dat')
  do iomega=1,PHdos%nomega
    if (InVar%Enunit.eq.0) write(56,'(2(f18.6,1x))') PHdos%omega(iomega)*Ha_eV*1000,PHdos%phdos(iomega)
    if (InVar%Enunit.eq.1) write(56,'(2(f18.6,1x))') PHdos%omega(iomega)*Ha_cmm1   ,PHdos%phdos(iomega)
    if (InVar%Enunit.eq.2) write(56,'(2(f18.6,1x))') PHdos%omega(iomega)           ,PHdos%phdos(iomega)
  end do
  close(56)

end subroutine tdep_calc_phdos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tdep_calc_thermo(DeltaFree_AH2,InVar,Lattice,PHdos,U0)

  integer :: iomega,itemp,iatom,itypat
  double precision :: k_B,wovert,heatcapa,entropy,internalE,freeE,expm2x,ln2shx,cothx,xx
  double precision :: Ftot,domega,MSD,Omega_m2,mass_amu,vdos
  double precision, intent(in) :: U0,DeltaFree_AH2
  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type), intent(in) :: Lattice
  type(phonon_dos_type),intent(in) :: PHdos

  write(InVar%stdout,*)' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '################# Thermodynamic quantities: Free energy,...##################'
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,'(a)') ' See the thermo.dat file'

! The heat capacity, entropy, internal and free energies (direct calculation)
! ===========================================================================
  mass_amu=zero
  do iatom=1,InVar%natom_unitcell
    itypat=InVar%typat_unitcell(iatom)
    mass_amu=mass_amu+InVar%amu(itypat)
  end do
  mass_amu=mass_amu*amu_emass/real(InVar%natom_unitcell)

  open(unit=20,file=trim(InVar%output_prefix)//'thermo.dat')
!FB  k_B=8.617343d-5 !in eV/K
  k_B=kb_HaK*Ha_eV
  domega=(InVar%dosdeltae*Ha_meV)
  wovert=1.d0/(2*InVar%temperature*k_B)
  heatcapa=0.d0 ; entropy=0.d0 ; internalE=0.d0 ; freeE=0.d0 ; MSD=0.d0 ; Omega_m2=0.d0 ; vdos=0.d0
  do iomega=1,PHdos%nomega
    xx=PHdos%omega(iomega)*Ha_eV
    if (xx.lt.tol8) cycle
    expm2x=exp(-2.d0*wovert*xx)
    ln2shx=wovert*xx+log(1.d0-expm2x)
    cothx=(1.d0+expm2x)/(1.d0-expm2x)
    heatcapa =heatcapa  + (wovert*xx/sinh(wovert*xx))**2*PHdos%phdos(iomega)*domega
    internalE=internalE + (wovert*xx)*cothx*PHdos%phdos(iomega)*domega
    entropy  =entropy   + ((wovert*xx)*cothx-ln2shx)*PHdos%phdos(iomega)*domega
    freeE    =freeE     + log(2*sinh(wovert*xx))*PHdos%phdos(iomega)*domega
    MSD      =MSD       + (cothx/PHdos%omega(iomega))*PHdos%phdos(iomega)*domega
    Omega_m2 =Omega_m2  + (1.d0/PHdos%omega(iomega))**2*PHdos%phdos(iomega)*domega
    vdos     =vdos      + PHdos%phdos(iomega)*domega
  end do
  heatcapa=heatcapa*3
  entropy=entropy*3
  internalE=internalE*3*k_B*InVar%temperature
  freeE=freeE*3*k_B*InVar%temperature
  MSD=MSD*3.d0/mass_amu/2.d0
  Omega_m2=Omega_m2*3.d0
  write(20,'(a)')'============= Direct results (without any inter/extrapolation) =================='
  write(20,'(1x,a,f10.5)')'For present temperature (in Kelvin): T= ',InVar%temperature
  write(20,'(1x,a,f12.5)')'  The cold contribution (in eV/atom): U_0 =',U0*Ha_eV
  write(20,'(1x,a,f10.5)')'  The specific heat (in k_b/atom): C_v=',heatcapa
  write(20,'(1x,a,f10.5)')'  The vibrational entropy (in k_b/atom): S_vib =',entropy
  write(20,'(1x,a,f10.5)')'  The internal energy (in eV/atom): U_vib =',internalE
  write(20,'(1x,a,f10.5)')'  The vibrational contribution (in eV/atom): F_vib = U_vib -T.S_vib =',freeE
  write(20,'(1x,a,f10.5)')'  The anharmonic contribution (in eV/atom): DeltaF_AH =',DeltaFree_AH2*Ha_eV
  Ftot=U0*Ha_eV+freeE+DeltaFree_AH2*Ha_eV
  write(20,'(1x,a)')'  So the free energy (in eV/atom) is equal to:'
  write(20,'(1x,a,f12.5)')'     Harmonic only -->  F_tot^HA = U_0 + F_vib =',U0*Ha_eV+freeE
  write(20,'(1x,a,f12.5)')'     With anharmonic contribution -->  F_tot^AH = U_0 + F_vib + DeltaF_AH =',Ftot
  write(20,'(1x,a)')'  Useful quantities for melting :'
  write(20,'(1x,a,f10.5)')'     The mean square displacement (in a.u.): sqrt(<u^2>) =',(MSD)**0.5
  write(20,'(1x,a,f10.5)')'     The <Omega^(-2)> factor (in THz^(-2)) =',Omega_m2/(Ha_THz)**2
  write(20,'(1x,a,f10.5)')'     The Wigner-Seitz radius (in a.u.) : d_at =',(6*Lattice%ucvol/pi)**(1./3.)
  write(20,'(1x,a,f10.5)')'     The average mass / proton-electron mass ratio (in a.u.) =', mass_amu/amu_emass
  write(20,'(1x,a,f10.5)')'     The Lindemann constant : sqrt(<u^2>)/d_at =',(MSD)**0.5/(6*Lattice%ucvol/pi)**(1./3.)
  write(20,'(1x,a,f10.5)')'     The integral of vDOS =',vdos
  write(20,'(a)')' '

! The free energy (extrapolation)
! ===============================
  write(20,'(a)')'============= Quasi-Harmonic Approximation (QHA) =================='
  write(20,'(1x,a)')'  Note that the following results come from an EXTRAPOLATION:'
  write(20,'(1x,a,i5,a)')'    1/ F_vib^QHA(T) is computed for each T using vDOS(T=',int(InVar%temperature),')'
  write(20,'(1x,a)')'    2/ F_tot^QHA(T) = F_vib^QHA(T) + U_0'
  write(20,'(1x,a)')'    3/ We assume that DeltaF_AH^QHA(T)=a(V)*T**2'
  write(20,'(1x,a)')'    4/ F_tot^QHA+AH(T) = U_0 + F_vib^QHA(T) + DeltaF_AH^QHA(T)'
  write(20,'(a)')'   T      F_vib^QHA(T)   F_tot^QHA(T)           C_v(T)  '&
&   //'       S_vib(T)        U_vib(T)   DeltaF_AH^QHA(T)   F_tot^QHA+AH(T)   MSD(T)'
  do itemp=1,100
    wovert=1.d0/(2*real(itemp)*100*k_B)
    freeE=0.d0
    heatcapa=0.d0
    entropy=0.d0
    internalE=0.d0
    MSD=0.d0
    do iomega=1,PHdos%nomega
      xx=PHdos%omega(iomega)*Ha_eV
      if (xx.lt.tol8) cycle
      expm2x=exp(-2.d0*wovert*xx)
      ln2shx=wovert*xx+log(1.d0-expm2x)
      cothx=(1.d0+expm2x)/(1.d0-expm2x)
      heatcapa =heatcapa  + (wovert*xx/sinh(wovert*xx))**2*PHdos%phdos(iomega)*domega
      internalE=internalE + (wovert*xx)*cothx*PHdos%phdos(iomega)*domega
      entropy  =entropy   + ((wovert*xx)*cothx-ln2shx)*PHdos%phdos(iomega)*domega
      freeE    =freeE     + log(2*sinh(wovert*xx))*PHdos%phdos(iomega)*domega
      MSD      =MSD       + (cothx/PHdos%omega(iomega))*PHdos%phdos(iomega)*domega
    end do
    heatcapa=heatcapa*3
    entropy=entropy*3
    internalE=internalE*3*k_B*itemp*100
    freeE=freeE*3*k_B*itemp*100
    MSD=MSD*3.d0/mass_amu/2.d0
    Ftot=U0*Ha_eV+freeE+DeltaFree_AH2*Ha_eV*(itemp*100)**2/(InVar%temperature*100)**2
    write(20,'(1x,i5,8(1x,f15.5))') itemp*100,freeE,U0*Ha_eV+freeE,heatcapa,entropy,internalE,&
&                          DeltaFree_AH2*Ha_eV*(itemp*100)**2/(InVar%temperature)**2,Ftot,(MSD)**0.5
  end do
  close(20)

end subroutine tdep_calc_thermo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tdep_calc_elastic(Phij_NN,distance,InVar,Lattice)

  integer :: iatom,ii,jj,kk,ll,iatcell,itypat
! integer :: istep
  double precision :: BH,BR,BV,GR,GV,GH,Eaverage,Nuaverage,Laverage,Vp,Vs,Vphi
  double precision :: rho,E1,E2,E3,Nu12,Nu13,Nu23,Nu21,Nu31,Nu32,G23,G13,G12
  double precision :: mass_amu,bohr
! real(dp) :: sigma_11,sigma_21,sigma_22,sigma_31,sigma_32,sigma_33
  double precision, allocatable :: Sij(:,:),Cij(:,:),aijkl(:,:,:,:),cijkl(:,:,:,:)
  type(Input_Variables_type), intent(in) :: InVar
  type(Lattice_Variables_type), intent(inout) :: Lattice
  double precision, intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision, intent(in) :: Phij_NN(3*InVar%natom,3*InVar%natom)

  write(InVar%stdout,*)' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '######################### Elastic constants #################################'
  write(InVar%stdout,*) '################ Bulk and Shear modulus--Sound velocities ###################'
  write(InVar%stdout,*) '#############################################################################'

  bohr=0.5291772108e-10
! Define atomic mass average
  mass_amu=zero
  do iatom=1,InVar%natom_unitcell
    itypat=InVar%typat_unitcell(iatom)
    mass_amu=mass_amu+InVar%amu(itypat)
  end do
  mass_amu=mass_amu/real(InVar%natom_unitcell)

  rho=(mass_amu/1e3)*InVar%natom_unitcell/Lattice%ucvol/bohr**3/6.022e23

!==========================================================================================
!===================== Elastic constants ==================================================
!==========================================================================================
!FB! Stress and Pressure
!FB  sigma_11=zero; sigma_21=zero
!FB  do istep=1,InVar%nstep
!FB    sigma_11=sigma_11+InVar%sigma(1,istep)
!FB    sigma_21=sigma_21+InVar%sigma(2,istep)
!FB  end do
!FB  sigma_11=sigma_11/real(InVar%nstep)
!FB  sigma_21=sigma_21/real(InVar%nstep)

! New calculation of elastic constants using the formula (12.28 and 12.29 of
! Wallace, Statistical physics of crystals and liquids, Worl Scientific)
  ABI_MALLOC(aijkl,(3,3,3,3)); aijkl(:,:,:,:)=0.d0
  ABI_MALLOC(cijkl,(3,3,3,3)); cijkl(:,:,:,:)=0.d0
  do ii=1,3
    do jj=1,3
      do kk=1,3
        do ll=1,3
          do iatom=1,InVar%natom
            do iatcell=1,InVar%natom_unitcell
              aijkl(ii,jj,kk,ll)=aijkl(ii,jj,kk,ll)-Phij_NN(ii+(iatcell-1)*3,3*(iatom-1)+jj)&
&               *distance(iatcell,iatom,kk+1)*distance(iatcell,iatom,ll+1)/2.d0/Lattice%ucvol
            end do
          end do
        enddo
      enddo
    enddo
  enddo

  do ii=1,3
    do jj=1,3
      do kk=1,3
        do ll=1,3
          cijkl(ii,jj,kk,ll)=aijkl(ii,kk,jj,ll)+aijkl(jj,kk,ii,ll)-aijkl(ii,jj,kk,ll)
!FB          cijkl(ii,jj,kk,ll)=aijkl(ii,kk,jj,ll)+aijkl(jj,kk,ii,ll)-aijkl(jj,ii,kk,ll)
        enddo
      enddo
    enddo
  enddo
  ABI_FREE(aijkl)

  cijkl(:,:,:,:)=cijkl(:,:,:,:)*29421.033d0

  ABI_MALLOC(Cij,(6,6)) ; Cij(:,:)=0.d0
  Cij(1,1)=cijkl(1,1,1,1) ; Cij(1,2)=cijkl(1,1,2,2) ; Cij(1,3)=cijkl(1,1,3,3)
  Cij(1,4)=cijkl(1,1,2,3) ; Cij(1,5)=cijkl(1,1,1,3) ; Cij(1,6)=cijkl(1,1,1,2)
  Cij(2,1)=cijkl(2,2,1,1) ; Cij(2,2)=cijkl(2,2,2,2) ; Cij(2,3)=cijkl(2,2,3,3)
  Cij(2,4)=cijkl(2,2,2,3) ; Cij(2,5)=cijkl(2,2,1,3) ; Cij(2,6)=cijkl(2,2,1,2)
  Cij(3,1)=cijkl(3,3,1,1) ; Cij(3,2)=cijkl(3,3,2,2) ; Cij(3,3)=cijkl(3,3,3,3)
  Cij(3,4)=cijkl(3,3,2,3) ; Cij(3,5)=cijkl(3,3,1,3) ; Cij(3,6)=cijkl(3,3,1,2)
  Cij(4,1)=cijkl(2,3,1,1) ; Cij(4,2)=cijkl(2,3,2,2) ; Cij(4,3)=cijkl(2,3,3,3)
  Cij(4,4)=cijkl(2,3,2,3) ; Cij(4,5)=cijkl(2,3,1,3) ; Cij(4,6)=cijkl(2,3,1,2)
  Cij(5,1)=cijkl(1,3,1,1) ; Cij(5,2)=cijkl(1,3,2,2) ; Cij(5,3)=cijkl(1,3,3,3)
  Cij(5,4)=cijkl(1,3,2,3) ; Cij(5,5)=cijkl(1,3,1,3) ; Cij(5,6)=cijkl(1,3,1,2)
  Cij(6,1)=cijkl(1,2,1,1) ; Cij(6,2)=cijkl(1,2,2,2) ; Cij(6,3)=cijkl(1,2,3,3)
  Cij(6,4)=cijkl(1,2,2,3) ; Cij(6,5)=cijkl(1,2,1,3) ; Cij(6,6)=cijkl(1,2,1,2)
! Remove the rounding errors before writing (for non regression testing purposes)
  do ii=1,6
    do jj=1,6
      if (abs(Cij(ii,jj)).lt.tol8) Cij(ii,jj)=zero
    end do
  end do
  write(InVar%stdout,'(a)') ' '
  write(InVar%stdout,'(a)') ' ========== Using the formulation proposed by Wallace (using the IFC) ========='
  write(InVar%stdout,'(a)') ' Cijkl='
  write(InVar%stdout,'(a,6(f8.3,1x))') ' | C11 C12 C13 C14 C15 C16 |   ',Cij(1,1),Cij(1,2),Cij(1,3),Cij(1,4),Cij(1,5),Cij(1,6)
  write(InVar%stdout,'(a,6(f8.3,1x))') ' | C21 C22 C23 C24 C25 C26 |   ',Cij(2,1),Cij(2,2),Cij(2,3),Cij(2,4),Cij(2,5),Cij(2,6)
  write(InVar%stdout,'(a,6(f8.3,1x))') ' | C31 C32 C33 C34 C35 C36 |   ',Cij(3,1),Cij(3,2),Cij(3,3),Cij(3,4),Cij(3,5),Cij(3,6)
  write(InVar%stdout,'(a,6(f8.3,1x))') ' | C41 C42 C43 C44 C45 C46 | = ',Cij(4,1),Cij(4,2),Cij(4,3),Cij(4,4),Cij(4,5),Cij(4,6)
  write(InVar%stdout,'(a,6(f8.3,1x))') ' | C51 C52 C53 C54 C55 C56 |   ',Cij(5,1),Cij(5,2),Cij(5,3),Cij(5,4),Cij(5,5),Cij(5,6)
  write(InVar%stdout,'(a,6(f8.3,1x))') ' | C61 C62 C63 C64 C65 C66 |   ',Cij(6,1),Cij(6,2),Cij(6,3),Cij(6,4),Cij(6,5),Cij(6,6)

! Mean value of the off-diagonal elements
  Cij(1,2)=(Cij(1,2)+Cij(2,1))/2.d0 ; Cij(2,1)=Cij(1,2)
  Cij(1,3)=(Cij(1,3)+Cij(3,1))/2.d0 ; Cij(3,1)=Cij(1,3)
  Cij(2,3)=(Cij(2,3)+Cij(3,2))/2.d0 ; Cij(3,2)=Cij(2,3)

! Young's modulus
  E1=(Cij(1,1)*Cij(2,2)*Cij(3,3)+2.d0*Cij(2,3)*Cij(1,2)*Cij(1,3)-Cij(1,1)*Cij(2,3)**2-Cij(2,2)*Cij(1,3)**2-Cij(3,3)*Cij(1,2)**2)&
&   /(Cij(2,2)*Cij(3,3)-Cij(2,3)**2)
  E2=(Cij(1,1)*Cij(2,2)*Cij(3,3)+2.d0*Cij(2,3)*Cij(1,2)*Cij(1,3)-Cij(1,1)*Cij(2,3)**2-Cij(2,2)*Cij(1,3)**2-Cij(3,3)*Cij(1,2)**2)&
&   /(Cij(1,1)*Cij(3,3)-Cij(1,3)**2)
  E3=(Cij(1,1)*Cij(2,2)*Cij(3,3)+2.d0*Cij(2,3)*Cij(1,2)*Cij(1,3)-Cij(1,1)*Cij(2,3)**2-Cij(2,2)*Cij(1,3)**2-Cij(3,3)*Cij(1,2)**2)&
&   /(Cij(1,1)*Cij(2,2)-Cij(1,2)**2)
  write(InVar%stdout,'(a,3(f8.3,1x))') ' Young modulus E1, E2 and E3=',E1,E2,E3

! Poisson Ratio
  Nu21=(Cij(1,2)*Cij(3,3)-Cij(1,3)*Cij(2,3))/(Cij(1,1)*Cij(3,3)-Cij(1,3)**2)
  Nu31=(Cij(1,3)*Cij(2,2)-Cij(1,2)*Cij(2,3))/(Cij(1,1)*Cij(2,2)-Cij(1,2)**2)
  Nu23=(Cij(1,1)*Cij(2,3)-Cij(1,2)*Cij(1,3))/(Cij(1,1)*Cij(3,3)-Cij(1,3)**2)
  Nu12=(Cij(1,2)*Cij(3,3)-Cij(1,3)*Cij(2,3))/(Cij(2,2)*Cij(3,3)-Cij(2,3)**2)
  Nu13=(Cij(2,2)*Cij(1,3)-Cij(1,2)*Cij(2,3))/(Cij(2,2)*Cij(3,3)-Cij(2,3)**2)
  Nu32=(Cij(1,1)*Cij(2,3)-Cij(1,2)*Cij(1,3))/(Cij(1,1)*Cij(2,2)-Cij(1,2)**2)
  write(InVar%stdout,'(a,6(f8.3,1x))') ' Poisson ratio Nu21, Nu31, Nu23, Nu12, Nu13 and Nu32=',Nu21,Nu31,Nu23,Nu12,Nu13,Nu32

! Shear modulus
  G23=Cij(4,4) ; G13=Cij(5,5) ; G12=Cij(6,6)
  write(InVar%stdout,'(a,3(f8.3,1x))') ' Shear modulus G23, G13 and G12=',G23,G13,G12

! Compliance matrix
  ABI_MALLOC(Sij,(6,6)) ; Sij(:,:)=0.d0
  Sij(1,1)= 1.d0/E1 ; Sij(1,2)=-Nu21/E2 ; Sij(1,3)=-Nu31/E3 ; Sij(1,4)=0.d0     ; Sij(1,5)=0.d0     ; Sij(1,6)=0.d0
  Sij(2,1)=-Nu12/E1 ; Sij(2,2)= 1.d0/E2 ; Sij(2,3)=-Nu32/E3 ; Sij(2,4)=0.d0     ; Sij(2,5)=0.d0     ; Sij(2,6)=0.d0
  Sij(3,1)=-Nu13/E1 ; Sij(3,2)=-Nu23/E2 ; Sij(3,3)= 1.d0/E3 ; Sij(3,4)=0.d0     ; Sij(3,5)=0.d0     ; Sij(3,6)=0.d0
  Sij(4,1)= 0.d0    ; Sij(4,2)= 0.d0    ; Sij(4,3)= 0.d0    ; Sij(4,4)=1.d0/G23 ; Sij(4,5)=0.d0     ; Sij(4,6)=0.d0
  Sij(5,1)= 0.d0    ; Sij(5,2)= 0.d0    ; Sij(5,3)= 0.d0    ; Sij(5,4)=0.d0     ; Sij(5,5)=1.d0/G13 ; Sij(5,6)=0.d0
  Sij(6,1)= 0.d0    ; Sij(6,2)= 0.d0    ; Sij(6,3)= 0.d0    ; Sij(6,4)=0.d0     ; Sij(6,5)=0.d0     ; Sij(6,6)=1.d0/G12
! Remove the rounding errors before writing (for non regression testing purposes)
  do ii=1,6
    do jj=1,6
      if (abs(Sij(ii,jj)).lt.tol8) Sij(ii,jj)=zero
    end do
  end do
  write(InVar%stdout,'(a)') ' '
  write(InVar%stdout,'(a)') ' Sijkl='
  write(InVar%stdout,'(a,6(f8.3,1x))') ' | S11 S12 S13 S14 S15 S16 |   ',Sij(1,1),Sij(1,2),Sij(1,3),Sij(1,4),Sij(1,5),Sij(1,6)
  write(InVar%stdout,'(a,6(f8.3,1x))') ' | S21 S22 S23 S24 S25 S26 |   ',Sij(2,2),Sij(2,2),Sij(2,3),Sij(2,4),Sij(2,5),Sij(2,6)
  write(InVar%stdout,'(a,6(f8.3,1x))') ' | S31 S32 S33 S34 S35 S36 |   ',Sij(3,1),Sij(3,2),Sij(3,3),Sij(3,4),Sij(3,5),Sij(3,6)
  write(InVar%stdout,'(a,6(f8.3,1x))') ' | S41 S42 S43 S44 S45 S46 | = ',Sij(4,1),Sij(4,2),Sij(4,3),Sij(4,4),Sij(4,5),Sij(4,6)
  write(InVar%stdout,'(a,6(f8.3,1x))') ' | S51 S52 S53 S54 S55 S56 |   ',Sij(5,1),Sij(5,2),Sij(5,3),Sij(5,4),Sij(5,5),Sij(5,6)
  write(InVar%stdout,'(a,6(f8.3,1x))') ' | S61 S62 S63 S64 S65 S66 |   ',Sij(6,1),Sij(6,2),Sij(6,3),Sij(6,4),Sij(6,5),Sij(6,6)

!==========================================================================================
!===================== Bulk and Shear modulus--Sound velocities ===========================
!==========================================================================================
! Voigt notation
  write(InVar%stdout,'(a,f9.3)')' For density rho=',rho
  write(InVar%stdout,*)' '
  write(InVar%stdout,*)' ========================= Voigt average (constant strain) ==================='
  BV=((Cij(1,1)+Cij(2,2)+Cij(3,3))+2.d0*(Cij(1,2)+Cij(1,3)+Cij(2,3)))/9.d0
  GV=((Cij(1,1)+Cij(2,2)+Cij(3,3))-     (Cij(1,2)+Cij(1,3)+Cij(2,3))+3.d0*(Cij(4,4)+Cij(5,5)+Cij(6,6)))/15.d0
  write(InVar%stdout,'(2(a,f9.3))')' ISOTHERMAL modulus: Bulk Kt=',BV,' and Shear G=',GV
  Eaverage=9.d0*BV*GV/(3*BV+GV)
  Nuaverage=0.5*(1.d0-(3.d0*GV)/(3.d0*BV+GV) )
  Laverage=(3.d0*BV-2.d0*GV)/3.d0
  write(InVar%stdout,'(3(a,f9.3))')' Average of Young modulus E=',Eaverage,' Lame modulus Lambda=',Laverage,&
&   ' and Poisson ratio Nu=',Nuaverage
  Vp=dsqrt(1.d9*(BV+4.d0*GV/3.d0)/rho)
  Vs=dsqrt(1.d9*GV/rho)
  Vphi=dsqrt(1.d9*BV/rho)
  write(InVar%stdout,'(3(a,f9.3,1x))')' Velocities: compressional Vp=',Vp,' shear Vs=',Vs,' and bulk Vphi=',Vphi

! Reuss notation
  write(InVar%stdout,*)' '
  write(InVar%stdout,*)' ========================= Reuss average (constant stress) ==================='
  BR=1.d0/(Sij(1,1)+Sij(2,2)+Sij(3,3)+2.d0*(Sij(1,2)+Sij(1,3)+Sij(2,3)))
  GR=15.d0/(4.d0*(Sij(1,1)+Sij(2,2)+Sij(3,3))-4.d0*(Sij(1,2)+Sij(1,3)+Sij(2,3))+3.d0*(Sij(4,4)+Sij(5,5)+Sij(6,6)))
  write(InVar%stdout,'(2(a,f9.3))')' ISOTHERMAL modulus: Bulk Kt=',BR,' and Shear G=',GR
  Eaverage=9.d0*BR*GR/(3*BR+GR)
  Nuaverage=0.5*(1.d0-(3.d0*GR)/(3.d0*BR+GR) )
  Laverage=(3.d0*BR-2.d0*GR)/3.d0
  write(InVar%stdout,'(3(a,f9.3))')' Average of Young modulus E=',Eaverage,' Lame modulus Lambda=',Laverage,&
&   ' and Poisson ratio Nu=',Nuaverage
  Vp=dsqrt(1.d9*(BR+4.d0*GR/3.d0)/rho)
  Vs=dsqrt(1.d9*GR/rho)
  Vphi=dsqrt(1.d9*BR/rho)
  write(InVar%stdout,'(3(a,f9.3,1x))')' Velocities: compressional Vp=',Vp,' shear Vs=',Vs,' and bulk Vphi=',Vphi

! Voigt-Reuss-Hill notation
  write(InVar%stdout,*)' '
  write(InVar%stdout,*)' ============================== Hill average ================================='
  BH=(BR+BV)/2.d0
  GH=(GR+GV)/2.d0
  write(InVar%stdout,'(2(a,f9.3))')' ISOTHERMAL modulus: Bulk Kt=',BH,' and Shear G=',GH
  Eaverage=9.d0*BH*GH/(3*BH+GH)
  Nuaverage=0.5*(1.d0-(3.d0*GH)/(3.d0*BH+GH) )
  Laverage=(3.d0*BH-2.d0*GH)/3.d0
  write(InVar%stdout,'(3(a,f9.3))')' Average of Young modulus E=',Eaverage,' Lame modulus Lambda=',Laverage,&
&   ' and Poisson ratio Nu=',Nuaverage
  Vp=dsqrt(1.d9*(BH+4.d0*GH/3.d0)/rho)
  Vs=dsqrt(1.d9*GH/rho)
  Vphi=dsqrt(1.d9*BH/rho)
  write(InVar%stdout,'(3(a,f9.3,1x))')' Velocities: compressional Vp=',Vp,' shear Vs=',Vs,' and bulk Vphi=',Vphi
! Store the RVH value of the bulk modulus (will be useful for the Gruneisen)
  Lattice%BulkModulus=BH

  ABI_FREE(cijkl)
  ABI_FREE(Cij)
  ABI_FREE(Sij)
end subroutine tdep_calc_elastic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module m_tdep_phdos
