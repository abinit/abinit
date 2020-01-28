
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!****m* ABINIT/m_CtqmcoffdiagInterface
!! NAME
!!  m_CtqmcoffdiagInterface
!! 
!! FUNCTION 
!!  Manage a ctqmc simulation. 
!!  friendly interface for the user
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder, B. Amadon, J. Denier)
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
MODULE m_CtqmcoffdiagInterface
USE m_Ctqmcoffdiag

IMPLICIT NONE

!!***

!!****t* m_CtqmcoffdiagInterface/CtqmcoffdiagInterface
!! NAME
!!  CtqmcoffdiagInterface
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

TYPE CtqmcoffdiagInterface
  TYPE(Ctqmcoffdiag) :: Hybrid
  INTEGER :: opt_fk       = 0
  INTEGER :: opt_order    = 0
  INTEGER :: opt_movie    = 0
  INTEGER :: opt_analysis = 0
  INTEGER :: opt_check    = 0
  INTEGER :: opt_spectra  = 0
  INTEGER :: opt_noise    = 0
  INTEGER :: opt_gMove    = 0
END TYPE CtqmcoffdiagInterface
!!***

CONTAINS
!!***

!!****f* ABINIT/m_CtqmcoffdiagInterface/CtqmcoffdiagInterface_init
!! NAME
!!  CtqmcoffdiagInterface_init
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
!!  op=ctqmcinterface
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

SUBROUTINE CtqmcoffdiagInterface_init(op,iseed,sweeps,thermalization,&
&measurements,flavors,samples,beta,U,ostream,MPI_COMM,opt_nondiag)

!Arguments ------------------------------------
  TYPE(CtqmcoffdiagInterface), INTENT(INOUT) :: op
  INTEGER, OPTIONAL, INTENT(IN) :: MPI_COMM
  INTEGER, INTENT(IN) :: iseed
  DOUBLE PRECISION, INTENT(IN) :: sweeps
  INTEGER, INTENT(IN) :: thermalization
  INTEGER, INTENT(IN) :: measurements
  INTEGER, INTENT(IN) :: flavors
  INTEGER, INTENT(IN) :: samples
  !INTEGER, INTENT(IN) :: Wmax
  INTEGER, INTENT(IN) :: ostream
  INTEGER, INTENT(IN) :: opt_nondiag
  DOUBLE PRECISION, INTENT(IN) :: beta
  DOUBLE PRECISION, INTENT(IN) :: u
  !DOUBLE PRECISION, INTENT(IN) :: mu
!Local arguements -----------------------------
  INTEGER          :: ifstream
  DOUBLE PRECISION, DIMENSION(1:10) :: buffer

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
  buffer(10)=DBLE(opt_nondiag)
  !buffer(9)=0.d0!mu
  !buffer(9)=DBLE(Wmax)

  IF ( PRESENT( MPI_COMM ) ) THEN
    CALL Ctqmcoffdiag_init(op%Hybrid, ostream, ifstream, .FALSE., MPI_COMM,buffer)
  ELSE
    CALL Ctqmcoffdiag_init(op%Hybrid, ostream, ifstream, .FALSE.,iBuffer=buffer)
  END IF
  op%opt_fk       = 0
  op%opt_order    = 0
  op%opt_movie    = 0
  op%opt_analysis = 0
  op%opt_check    = 0
  op%opt_noise    = 0
  op%opt_spectra  = 0
END SUBROUTINE CtqmcoffdiagInterface_init
!!***

!!****f* ABINIT/m_CtqmcoffdiagInterface/CtqmcoffdiagInterface_setOpts
!! NAME
!!  CtqmcoffdiagInterface_setOpts
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
!!  op=ctqmcinterface
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

SUBROUTINE CtqmcoffdiagInterface_setOpts(op,opt_Fk,opt_order,opt_movie,opt_analysis,opt_check, opt_noise, opt_spectra, opt_gMove) 

!Arguments ------------------------------------
  TYPE(CtqmcoffdiagInterface), INTENT(INOUT) :: op
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_Fk
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_order
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_movie
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_analysis
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_check
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_noise
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_spectra
  INTEGER , OPTIONAL  , INTENT(IN   ) :: opt_gMove

  IF ( PRESENT(opt_Fk) ) &
    op%opt_Fk = opt_fk
  IF ( PRESENT(opt_order) ) &
    op%opt_order = opt_order
  IF ( PRESENT(opt_analysis) ) &
    op%opt_analysis = opt_analysis
  IF ( PRESENT(opt_check) ) &
    op%opt_check = opt_check
  IF ( PRESENT(opt_movie) ) &
    op%opt_movie = opt_movie
  IF ( PRESENT(opt_noise) ) &
    op%opt_noise = opt_noise
  IF ( PRESENT(opt_spectra) ) &
    op%opt_spectra = opt_spectra
  IF ( PRESENT(opt_gMove) ) &
    op%opt_gMove = opt_gMove

END SUBROUTINE CtqmcoffdiagInterface_setOpts
!!***

!!****f* ABINIT/m_CtqmcoffdiagInterface/CtqmcoffdiagInterface_run
!! NAME
!!  CtqmcoffdiagInterface_run
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
!!  op=ctqmcinterface
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

SUBROUTINE CtqmcoffdiagInterface_run(op,G0omega, Gtau, Gw, D,E,Noise,matU,Docc,opt_sym,opt_levels,hybri_limit) 

!Arguments ------------------------------------
  TYPE(CtqmcoffdiagInterface), INTENT(INOUT) :: op
  COMPLEX(KIND=8)      , DIMENSION(:,:,:), INTENT(IN ) :: G0omega
  DOUBLE PRECISION, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT) :: Gtau
  COMPLEX(KIND=8)      , DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: Gw
  DOUBLE PRECISION, OPTIONAL      , INTENT(OUT) :: D
  DOUBLE PRECISION, OPTIONAL      , INTENT(OUT) :: E
  DOUBLE PRECISION, OPTIONAL      , INTENT(OUT) :: Noise
  DOUBLE PRECISION, DIMENSION(:,:),OPTIONAL,  INTENT(IN ) :: matU
  DOUBLE PRECISION, DIMENSION(:,:),OPTIONAL,  INTENT(OUT ) :: Docc
  DOUBLE PRECISION, DIMENSION(:,:),OPTIONAL,  INTENT(IN ) :: opt_sym
  DOUBLE PRECISION, DIMENSION(:),OPTIONAL,  INTENT(IN ) :: opt_levels
  COMPLEX(KIND=8) , DIMENSION(:,:),OPTIONAL,  INTENT(IN ) :: hybri_limit

  CALL Ctqmcoffdiag_reset(op%Hybrid)

!  ifstream = 42
!
!  OPEN(UNIT=ifstream, FILE="Gw.dat")
!  CALL Ctqmcoffdiag_setG0w(Hybrid, ifstream)
!  CLOSE(ifstream)
!  

  IF ( PRESENT(opt_levels)) &
    CALL Ctqmcoffdiag_setMu(op%Hybrid, opt_levels)

  IF ( PRESENT(hybri_limit)) &
    CALL Ctqmcoffdiag_sethybri_limit(op%Hybrid, hybri_limit)

  CALL Ctqmcoffdiag_setG0wTab(op%Hybrid, G0omega,op%opt_fk)

  IF ( PRESENT(matU) ) &
    CALL Ctqmcoffdiag_setU(op%Hybrid, matU)

  CALL Ctqmcoffdiag_run(op%Hybrid,opt_order=op%opt_order, &
                           opt_movie=op%opt_movie, &
                           opt_analysis=op%opt_analysis, &
                           opt_check=op%opt_check, &
                           opt_noise=op%opt_noise, &
                           opt_spectra=op%opt_spectra, &
                           opt_gMove=op%opt_gMove)

 ! write(6,*) "op%Hybrid%stats",op%Hybrid%stats
 ! write(6,*) "opt_gMove",op%opt_gMove

  CALL Ctqmcoffdiag_getResult(op%Hybrid)

  IF ( PRESENT(opt_sym) ) THEN
    CALL Ctqmcoffdiag_symmetrizeGreen(op%Hybrid,opt_sym)
  END IF

 ! write(6,*) "op%Hybrid%stats",op%Hybrid%stats

  IF ( PRESENT(Gtau) .AND. PRESENT(Gw) ) THEN
    CALL Ctqmcoffdiag_getGreen(op%Hybrid, Gtau=Gtau, Gw=Gw)
  !      write(6,*) "size",size(Gw,dim=1),size(gw,dim=2)
  !      IF ( hybrid%rank .EQ. 0 ) write(389,*) Gw(:,hybrid%flavors+1)
  !      call flush(389)
  ELSE IF ( PRESENT(Gtau) .AND. .NOT. PRESENT(Gw) ) THEN
    CALL Ctqmcoffdiag_getGreen(op%Hybrid, Gtau=Gtau)
  ELSE IF ( .NOT. PRESENT(Gtau) .AND. PRESENT(Gw) ) THEN
    CALL Ctqmcoffdiag_getGreen(op%Hybrid, Gw=Gw)
  END IF

  !write(6,*) "op%Hybrid%stats",op%Hybrid%stats

  IF ( PRESENT(D) ) &
  CALL Ctqmcoffdiag_getD(op%Hybrid, D)
  Docc=op%Hybrid%measDE

  !write(6,*) "op%Hybrid%stats",op%Hybrid%stats

  IF ( PRESENT(E) .AND. PRESENT(Noise) ) &
  CALL Ctqmcoffdiag_getE(op%Hybrid, E, Noise)

  CALL Ctqmcoffdiag_printAll(op%Hybrid)
  !CALL Ctqmcoffdiag_printQMC(op%Hybrid)

END SUBROUTINE CtqmcoffdiagInterface_run
!!***

!!****f* ABINIT/m_CtqmcoffdiagInterface/CtqmcoffdiagInterface_setSweeps
!! NAME
!!  CtqmcoffdiagInterface_setSweeps
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
!!  op=ctqmcinterface
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

SUBROUTINE CtqmcoffdiagInterface_setSweeps(op, sweeps)

!Arguments ------------------------------------
  TYPE(CtqmcoffdiagInterface), INTENT(INOUT) :: op
  DOUBLE PRECISION, INTENT(IN) :: sweeps

  CALL Ctqmcoffdiag_setSweeps(op%Hybrid,sweeps)
END SUBROUTINE CtqmcoffdiagInterface_setSweeps
!!***

!!****f* ABINIT/m_CtqmcoffdiagInterface/CtqmcoffdiagInterface_finalize
!! NAME
!!  CtqmcoffdiagInterface_finalize
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
!!  op=ctqmcinterface
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

SUBROUTINE CtqmcoffdiagInterface_finalize(op)

!Arguments ------------------------------------
  TYPE(CtqmcoffdiagInterface), INTENT(INOUT) :: op

  !IF ( op%Hybrid%init .EQV. .TRUE. ) THEN
!    CALL Ctqmcoffdiag_printAll(op%Hybrid)
         !write(6,*) "before ctqmc_destroy in Ctqmcoffdiaginterface_finalize"
    CALL Ctqmcoffdiag_destroy(op%Hybrid)
  !END IF

END SUBROUTINE CtqmcoffdiagInterface_finalize
!!***

END MODULE m_CtqmcoffdiagInterface
!!***
