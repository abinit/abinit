#if defined HAVE_CONFIG_H
#include "config.h"

#include "abi_common.h"

#define MALLOC(ARR,SIZE)     ABI_ALLOCATE(ARR,SIZE)
#define FREE(ARR)            ABI_DEALLOCATE(ARR) 
#define FREEIF(ARR)          IF(ALLOCATED(ARR)) THEN NEWLINE ABI_DEALLOCATE(ARR) NEWLINE END IF

#define DT_MALLOC(ARR,SIZE)  ABI_DATATYPE_ALLOCATE(ARR,SIZE)
#define DT_FREE(ARR)         ABI_DATATYPE_DEALLOCATE(ARR) 
#define DT_FREEIF(ARR)       IF(ALLOCATED(ARR)) THEN NEWLINE ABI_DATATYPE_DEALLOCATE(ARR) NEWLINE END IF

#define myWARNALL(msg)       ABI_WARNING(msg)
#define myWARN(msg)          call msg_hndl(msg,"WARNING","PERS")
#define myERROR(msg)         ABI_ERROR(msg) 
#define MY_WORLD             xmpi_world

#define _PRIVATE              ABI_PRIVATE


#else

#define MALLOC(ARR,SIZE)     ALLOCATE(ARR SIZE)
#define FREE(ARR)            DEALLOCATE(ARR)
#define FREEIF(ARR)          IF(ALLOCATED(ARR)) DEALLOCATE(ARR)
#define DT_MALLOC(ARR,SIZE)  ALLOCATE(ARR SIZE) 
#define DT_FREE(ARR)         DEALLOCATE(ARR) 
#define DT_FREEIF(ARR)       IF(ALLOCATED(ARR)) DEALLOCATE(ARR)

#define std_err                6
#define myWARNALL(msg)         WRITE(std_err,'(A)') msg
#define myWARN(msg)            WRITE(std_err,'(A)') msg
#define myERROR(msg)           WRITE(std_err,'(A)') msg

#define MY_WORLD               MPI_COMM_WORLD

#define _PRIVATE                ,PRIVATE

#ifdef HAVE_MPI
#define HAVE_MPI2
#endif
  
#endif


#define Global_SIZE 100
#define MODCYCLE(a,b,c) c=a; IF(c .GT. b) c = c-b;
#define Vector_QuickResize(a,b) IF( b .GT. a%size ) CALL Vector_enlarge(a,MAX(b-a%size,Global_SIZE)); a%tail = b
#define VectorInt_QuickResize(a,b) IF( b .GT. a%size ) CALL VectorInt_enlarge(a,MAX(b-a%size,Global_SIZE)); a%tail = b
#define ImpurityOperator_QuickActivation(a,b) a%activeFlavor = b
#define BathOperator_QuickActivation(a,b) a%activeFlavor = b; a%MAddFlag = .FALSE.;a%MRemoveFlag = .FALSE.
#define Cdag_ 1
#define C_    2
#define GREENHYB_TAU    0
#define GREENHYB_OMEGA  1


#ifdef HAVE_MPI2
#define MYMPI2 use mpi
#else
#define MYMPI2 !I know it is useless
#endif
