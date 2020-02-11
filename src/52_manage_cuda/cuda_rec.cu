/* cuda_rec.cu */

/*
 * Copyright (C) 2008-2020 ABINIT Group (MMancini)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACE WITH FORTRAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <stdlib.h>
#include <unistd.h>
#include "cuda_header.h"

void 
recursion_no_bth(
		 const int ,const int ,const int ,const int ,      
		 const int ,int* ,const cureal ,const cureal ,
		 const cureal ,const cureal ,const int3* ,
		 const int3* , int* , const cureal* ,const cureal* ,       
		 cureal* ,cureal* );  
  
void 
recursion_bth(
	      const int, const int, const int, const int,           
	      const int, int* ,const cureal ,const cureal ,
	      const cureal ,const cureal ,const int3*,          
	      const int3*,int*, const cureal*, 
	      const cureal* , cureal* ,cureal* );
				      
extern "C"
void cuda_rec_cal_(const int* trotter, 
		   const int* gratio,
		   const int* npt,
		   const int* nrec,              //- Max number of recursion
		   const int* height_max,        //- Max number of points at the same time (depends on GPU)
		   int* max_rec,                 //- Out-max recursion to have required precision
		   const cureal* beta,
		   const cureal* fermie,
		   const cureal* tolrec,         //- Tollerance of recursion
		   const cureal* inf_vol,        //- Infinitesimal volume
		   const int3* pt0,              //- Coord of initial point
		   const int3* pt1,              //- Coord of final point
		   int* ngfftrec,                //- Linear sizes of the grid
		   const cureal* T_p,            //- The green kernel
		   const cureal* pot,            //- Potential
		   cureal* an,cureal* bn2)       //- Rec coefficients an and bn2
{
  printf("\n-------cudarec---------- \n"); 
  /*    printf("pt0 %d %d %d \n",pt0->x,pt0->y,pt0->z); */
  /*    printf("pt1 %d %d %d \n",pt1->x,pt1->y,pt1->z); */
  /*    printf("ngfftrec %d %d \n",ngfftrec[0],ngfftrec[2]); */
  /*    printf("nrec %d \n",*nrec); */
  /*    printf("size cureal %d \n",int(sizeof(cureal))); */
  /*    printf("tolrec %e \n",*tolrec); */
  /*    printf("inf_vol %f \n",*inf_vol); */
  /*    printf("maxpt %d \n",*height_max); */

  if((ngfftrec[1]*ngfftrec[2]*ngfftrec[3])%16==0){

    recursion_bth(*trotter,*gratio,*npt,
		  *nrec,*height_max,max_rec,
		  *beta,*fermie,*tolrec,*inf_vol,
		  pt0,pt1,ngfftrec,T_p,pot,an,bn2);
  }
  else{
     
    recursion_no_bth(*trotter,*gratio,*npt,
		     *nrec,*height_max,max_rec,
		     *beta,*fermie,*tolrec,*inf_vol,
		     pt0,pt1,ngfftrec,T_p,pot,an,bn2);
  }

  return;
}
