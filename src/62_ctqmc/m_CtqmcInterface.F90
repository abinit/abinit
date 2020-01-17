
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!****m* ABINIT/m_CtqmcInterface
!! NAME
!!  m_CtqmcInterface
!! 
!! FUNCTION 
!!  Manage a ctqmc simulation. 
!!  friendly interface for the user
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

#include "defs.h"
MODULE m_CtqmcInterface
USE m_Ctqmc

IMPLICIT NONE

!!***

PRIVATE

!!****t* m_CtqmcInterface/CtqmcInterface
!! NAME
!!  CtqmcInterface
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

TYPE, PUBLIC :: CtqmcInterface
  TYPE(Ctqmc)         :: Hybrid
  INTEGER _PRIVATE :: opt_fk       = 0
  INTEGER _PRIVATE :: opt_order    = 0
  INTEGER _PRIVATE :: opt_movie    = 0
  INTEGER _PRIVATE :: opt_analysis = 0
  INTEGER _PRIVATE :: opt_check    = 0
  INTEGER _PRIVATE :: opt_spectra  = 0
  INTEGER _PRIVATE :: opt_noise    = 0
  INTEGER _PRIVATE :: opt_gMove    = 0
END TYPE CtqmcInterface
!!***

PUBLIC :: CtqmcInterface_init
PUBLIC :: CtqmcInterface_setOpts
PUBLIC :: CtqmcInterface_run
PUBLIC :: CtqmcInterface_setSweeps
PUBLIC :: CtqmcInterface_finalize

CONTAINS
!!***

!!****f* ABINIT/m_CtqmcInterface/CtqmcInterface_init
!! NAME
!!  CtqmcInterface_init
!!
!! FUNCTION
!!  Initialize with permanent parameters
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmcinterface
!!  iseed=seed for rng
!!  sweeps=number of sweeps for the run
!!  thermalization=number of sweeps to thermalize
!!  measurements=how often we measure (modulo)
!!  flavors=number of orbitals (spin degenerated)
!!  samples=imaginary time slices
!!  beta=inverse temperature
!!  U=interaction parameter
!!  ostream=where to write output
!!  MPI_COMM=mpi communicator for the run
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE CtqmcInterface_init(this,iseed,sweeps,thermalization,measurements,flavors,samples,beta,U,ostream,MPI_COMM)

!Arguments ------------------------------------
  TYPE(CtqmcInterface), INTENT(INOUT) :: this
  INTEGER, OPTIONAL, INTENT(IN) :: MPI_COMM
  INTEGER, INTENT(IN) :: iseed
  DOUBLE PRECISION, INTENT(IN) :: sweeps
  INTEGER, INTENT(IN) :: thermalization
  INTEGER, INTENT(IN) :: measurements
  INTEGER, INTENT(IN) :: flavors
  INTEGER, INTENT(IN) :: samples
  !INTEGER, INTENT(IN) :: Wmax
  INTEGER, INTENT(IN) :: ostream
  DOUBLE PRECISION, INTENT(IN) :: beta
  DOUBLE PRECISION, INTENT(IN) :: u
  !DOUBLE PRECISION, INTENT(IN) :: mu
!Local arguements -----------------------------
  INTEGER          :: ifstream!,opt_nondiag
  DOUBLE PRECISION, DIMENSION(1:9) :: buffer
 ! opt_nondiag=0

  ifstream = 42

  buffer(1)=DBLE(iseed)
  buffer(2)=sweeps
  buffer(3)=DBLE(thermalization)
  buffer(4)=DBLE(measurements)
  buffer(5)=DBLE(flavors)
  buffer(6)=DBLE(samples)
  buffer(7)=beta
  buffer(8)=U
  buffer(9)=GREENHYB_TAU
 ! buffer(10)=DBLE(opt_nondiag)
  !buffer(9)=0.d0!mu
  !buffer(9)=DBLE(Wmax)

  IF ( PRESENT( MPI_COMM ) ) THEN
    CALL Ctqmc_init(this%Hybrid, ostream, ifstream, .FALSE., MY_COMM=MPI_COMM,iBuffer=buffer)
  ELSE
    CALL Ctqmc_init(this%Hybrid, ostream, ifstream, .FALSE.,iBuffer=buffer)
  END IF
  this%opt_fk       = 0
  this%opt_order    = 0
  this%opt_movie    = 0
  this%opt_analysis = 0
  this%opt_check    = 0
  this%opt_noise    = 0
  this%opt_spectra  = 0
END SUBROUTINE CtqmcInterface_init
!!***

!!****f* ABINIT/m_CtqmcInterface/CtqmcInterface_setOpts
!! NAME
!!  CtqmcInterface_setOpts
!!
!! FUNCTION
!!  Set and save options for many runs
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmcinterface
!!  opt_Fk=0 if we give us Gw0 and 1 if 1/Gw0+iwn
!!  opt_order=maximal perturbation order to scope(>0)
!!  opt_movie=print a latex file (0 or 1)
!!  opt_analysis=measure correlations (0 or 1)
!!  opt_check=check fast calculation :0 nothing
!!                                    1 Impurity
!!                                    2 Bath
!!                                    3 Both
!!  opt_noise=calculate noise ofr green functions(0 ro 1)
!!  opt_spectra=fourier transform of time evolution of number of electrons
!!               (0 or 1)
!!  opt_gMove=number of global moves (>0)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE CtqmcInterface_setOpts(this,opt_Fk,opt_order,opt_movie,opt_analysis,opt_check, opt_noise, opt_spectra, opt_gMove) 

!Arguments ------------------------------------
  TYPE(CtqmcInterface), INTENT(INOUT) :: this
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_Fk
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_order
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_movie
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_analysis
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_check
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_noise
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_spectra
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_gMove

  IF ( PRESENT(opt_Fk) ) &
    this%opt_Fk = opt_fk
  IF ( PRESENT(opt_order) ) &
    this%opt_order = opt_order
  IF ( PRESENT(opt_analysis) ) &
    this%opt_analysis = opt_analysis
  IF ( PRESENT(opt_check) ) &
    this%opt_check = opt_check
  IF ( PRESENT(opt_movie) ) &
    this%opt_movie = opt_movie
  IF ( PRESENT(opt_noise) ) &
    this%opt_noise = opt_noise
  IF ( PRESENT(opt_spectra) ) &
    this%opt_spectra = opt_spectra
  IF ( PRESENT(opt_gMove) ) &
    this%opt_gMove = opt_gMove

END SUBROUTINE CtqmcInterface_setOpts
!!***

!!****f* ABINIT/m_CtqmcInterface/CtqmcInterface_run
!! NAME
!!  CtqmcInterface_run
!!
!! FUNCTION
!!  run a ctqmc simu and get results
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmcinterface
!!  G0omega=Gw0 (according to opt_Fk)
!!  matU=interaction matrice
!!  opt_sym=weight factors to symmetrise G
!!  opt_levels=energy for each level (with respect to fermi level)
!!
!! OUTPUT
!!  Gtau=G(tau)
!!  Gw=fourier transform of Gtau
!!  D=full double occupancy
!!  E=Interaction energy
!!  Noise=Noise on E
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE CtqmcInterface_run(this,G0omega, Gtau, Gw, D,E,Noise,matU,opt_sym,opt_levels) 

!Arguments ------------------------------------
  TYPE(CtqmcInterface), INTENT(INOUT) :: this
  COMPLEX(KIND=8) , DIMENSION(:,:)          , INTENT(IN   ) :: G0omega
  DOUBLE PRECISION, DIMENSION(:,:), OPTIONAL, INTENT(  OUT) :: Gtau
  COMPLEX(KIND=8) , DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: Gw
  DOUBLE PRECISION, DIMENSION(:,:), OPTIONAL, INTENT(  OUT) :: D
  DOUBLE PRECISION                , OPTIONAL, INTENT(  OUT) :: E
  DOUBLE PRECISION                , OPTIONAL, INTENT(  OUT) :: Noise
  DOUBLE PRECISION, DIMENSION(:,:), OPTIONAL, INTENT(IN   ) :: matU
  DOUBLE PRECISION, DIMENSION(:,:), OPTIONAL, INTENT(IN   ) :: opt_sym
  DOUBLE PRECISION, DIMENSION(:  ), OPTIONAL, INTENT(IN   ) :: opt_levels

  CALL Ctqmc_reset(this%Hybrid)

!  ifstream = 42
!
!  OPEN(UNIT=ifstream, FILE="Gw.dat")
!  CALL Ctqmc_setG0w(Hybrid, ifstream)
!  CLOSE(ifstream)
!  

  IF ( PRESENT(opt_levels) ) &
    CALL Ctqmc_setMu(this%Hybrid, opt_levels)

  CALL Ctqmc_setG0wTab(this%Hybrid, G0omega,this%opt_fk)

  IF ( PRESENT(matU) ) &
    CALL Ctqmc_setU(this%Hybrid, matU)

  CALL Ctqmc_run(this%Hybrid,opt_order=this%opt_order, &
                           opt_movie=this%opt_movie, &
                           opt_analysis=this%opt_analysis, &
                           opt_check=this%opt_check, &
                           opt_noise=this%opt_noise, &
                           opt_spectra=this%opt_spectra, &
                           opt_gMove=this%opt_gMove)

  CALL Ctqmc_getResult(this%Hybrid)

  IF ( PRESENT(opt_sym) ) THEN
    CALL Ctqmc_symmetrizeGreen(this%Hybrid,opt_sym)
  END IF

  IF ( PRESENT(Gtau) .AND. PRESENT(Gw) ) THEN
    CALL Ctqmc_getGreen(this%Hybrid, Gtau=Gtau, Gw=Gw)
  !      write(6,*) "size",size(Gw,dim=1),size(gw,dim=2)
  !      IF ( hybrid%rank .EQ. 0 ) write(389,*) Gw(:,hybrid%flavors+1)
  !      call flush(389)
  ELSE IF ( PRESENT(Gtau) .AND. .NOT. PRESENT(Gw) ) THEN
    CALL Ctqmc_getGreen(this%Hybrid, Gtau=Gtau)
  ELSE IF ( .NOT. PRESENT(Gtau) .AND. PRESENT(Gw) ) THEN
    CALL Ctqmc_getGreen(this%Hybrid, Gw=Gw)
  END IF

  IF ( PRESENT(D) ) &
  CALL Ctqmc_getD(this%Hybrid, D)

  IF ( PRESENT(E) .AND. PRESENT(Noise) ) THEN
    CALL Ctqmc_getE(this%Hybrid, E=E, Noise=Noise)
  ELSE IF ( PRESENT(E) .AND. .NOT. PRESENT(Noise) ) THEN
    CALL Ctqmc_getE(this%Hybrid, E=E)
  ELSE IF ( .NOT. PRESENT(E) .AND. PRESENT(Noise) ) THEN
    CALL Ctqmc_getE(this%Hybrid, Noise=noise)
  END IF


  CALL Ctqmc_printAll(this%Hybrid)
  !CALL Ctqmc_printQMC(this%Hybrid)

END SUBROUTINE CtqmcInterface_run
!!***

!!****f* ABINIT/m_CtqmcInterface/CtqmcInterface_setSweeps
!! NAME
!!  CtqmcInterface_setSweeps
!!
!! FUNCTION
!!  change sweeps on the fly
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmcinterface
!!  sweeps=new number of sweeps
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE CtqmcInterface_setSweeps(this, sweeps)

!Arguments ------------------------------------
  TYPE(CtqmcInterface), INTENT(INOUT) :: this
  DOUBLE PRECISION, INTENT(IN) :: sweeps

  CALL Ctqmc_setSweeps(this%Hybrid,sweeps)
END SUBROUTINE CtqmcInterface_setSweeps
!!***

!!****f* ABINIT/m_CtqmcInterface/CtqmcInterface_finalize
!! NAME
!!  CtqmcInterface_finalize
!!
!! FUNCTION
!!  Destroy simulation
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  this=ctqmcinterface
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE CtqmcInterface_finalize(this)

!Arguments ------------------------------------
  TYPE(CtqmcInterface), INTENT(INOUT) :: this

  !IF ( this%Hybrid%init .EQV. .TRUE. ) THEN
!    CALL Ctqmc_printAll(this%Hybrid)
    CALL Ctqmc_destroy(this%Hybrid)
  !END IF

END SUBROUTINE CtqmcInterface_finalize
!!***

END MODULE m_CtqmcInterface
!!***
