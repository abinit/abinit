/* recursion_cut_no_bth.cu */

/*
 * Copyright (C) 2008-2020 ABINIT Group (MMancini)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#include "cuda_common.h"
#include "cuda_header.h"
#include "cuda_rec_head.h"
//#include "rec_kernels.cu"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* This CUDA module contains the functions (kernels) to perform Recursion
   Method DFT on GPU devices when recrcut!=0.

   Kernels function:
  
   Host function:
   recursion_cut_no_bth
*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void 
recursion_cut_no_bth(
		     const int trotter, 
		     const int gratio, 
		     const int npt, 
		     const int nrec,           //- Max number of recursion
		     const int nptrec,         //- Max number of points at the same time (depends on GPU)
		     int* max_rec,             //- Out-max recursion to have required precision
		     const cureal beta,
		     const cureal fermie,
		     const cureal tolrec,      //- Tollerance of recursion
		     const cureal inf_vol,     //- Infinitesimal volume
		     const int3* pt0,          //- Coord of initial point
		     const int3* pt1,	       //- Coord of final point
		     const int* ngfft,         //- Total spatial grid
		     int* ngfftrec,            //- Linear sizes of the grid
		     const cureal* T_p,        //- The green kernel
		     const cureal* pot,        //- Potential
		     cureal* an,cureal* bn2)   //- Rec coefficients an and bn2

/* NOTES:
   -The vectors for which to calculate the recursion are set in a
   matrix of size (pth_size X height_max).
   -pth_size= is calculated starting of nfftrec (pth_size>=nfftrec) to make
   GPU memory copy efficient.
   -height_max= is the max number of vector(of size nfftrec) which can be
   putted in the matrix  to make the recursion and it depends on the device.
   -At any call, the recursion is computed for npt points, by a single allocation on the device. 
   -The calculation is made only on min(height_max,pos[1]-pos[0])
   points where pos is the current number of point to compute.
   -oldtonew, for any step in recursion, for any point, the un,unold
   inthe next step are obtaind.
*/
	        
{
 /*------------------------------------------------------------------------------------*/
 /*---------------------------------INITIALIZATION-------------------------------------*/
 /*------------------------------------------------------------------------------------*/
 
 /*-------------- Time setting -----------*/
 float* timing = (float*)calloc(DEBUGLEN,sizeof(float));
 cudaEvent_t start; cudaEventCreate(&start);
 cudaEvent_t stop;  cudaEventCreate(&stop);

 starttime(&start);//-Start timing: memory allocation+setting

 /*--------------- Sizes of Vectors ---------------*/
 int nfftrec = ngfftrec[0]*ngfftrec[1]*ngfftrec[2];         //- Size of the 3d-grid
 int nfft = ngfft[0]*ngfft[1]*ngfft[2];         //- Size of the total 3d-grid   
 int target = (ngfftrec[0]>>1)*(1+ngfftrec[1]*(1+ngfftrec[2]));

 size_t height_max = (size_t) nptrec;
 size_t totallarg  = (size_t)nfft*sizeof(cureal);   //- Size of real vectors
 size_t largeur  = (size_t)nfftrec*sizeof(cureal);  //- Size of real vectors
 size_t clargeur = (size_t)nfftrec*sizeof(cucmplx); //- Size of Complex vectors
 size_t un_pitch = largeur;  //- Pitch to put multi-vectors in a matrix: intial guess
  
 /*------------------ Grids and Blocks ---------------------------*/
 //-Multi vector grid
 dim3 block(320, 1);
 dim3 grid(((size_t)nfftrec+block.x-1)/block.x,(height_max+block.y-1)/block.y);

 //-Grid to select potential
 dim3 block_pot(16,16);
 dim3 grid_pot(((size_t)ngfftrec[0]+block_pot.x-1)/block_pot.x,((size_t)ngfftrec[1]+block_pot.y-1)/block_pot.y);

 /*--------------------- FFT Planes ------------------------------*/
 cufftHandle plan_dir;
 cufftPlan3d(&plan_dir,ngfftrec[0],ngfftrec[1],ngfftrec[2],FFT_C2C);	  
  
 /*----------------- FFT of the Green Kernel ---------------------*/ 
 //-Get the green kernel  from host 
 cureal *T_p_gpu = NULL;
 cudaMalloc((void**)&T_p_gpu,largeur);
 cudaMemcpy(T_p_gpu,T_p,largeur,cudaMemcpyHostToDevice);
 //-Compute the FFT 
 cucmplx *ZT_p_gpu = NULL;
 cudaMalloc((void**)&ZT_p_gpu,clargeur);
 /*Obtain the FFT of the Green Kernel on device */ 
 realtocmplx <<< ((size_t)nfftrec+320-1)/320,320 >>>(T_p_gpu, ZT_p_gpu,nfftrec);
 cudaFree(T_p_gpu);
 FFTEXECC2C(plan_dir,ZT_p_gpu ,ZT_p_gpu,CUFFT_FORWARD);

 /*------------- Allocation of Matrices on Device ----------------*/
 cureal* vn_gpu    = NULL;
 cudaMallocPitch((void**) &vn_gpu,&un_pitch,largeur,height_max);  

 cureal* un_gpu    = NULL;
 cudaMallocPitch((void**) &un_gpu,&un_pitch,largeur,height_max);

 cureal* unold_gpu = NULL;
 cudaMallocPitch((void**) &unold_gpu,&un_pitch,largeur,height_max);

 cureal* an_gpu    = NULL;
 cudaMalloc((void**)&an_gpu,height_max*sizeof(cureal));

 cureal* bn2_gpu   = NULL;
 cudaMalloc((void**)&bn2_gpu,height_max*sizeof(cureal));

 /*----------------- Get the Potential from Host -----------------*/
 cureal* pot_gpu   = NULL;
 cudaMalloc((void**)&pot_gpu,totallarg);
 cudaMemcpy(pot_gpu,pot,totallarg,cudaMemcpyHostToDevice);

 /*----------------- Local Potential -----------------*/
 cureal* locpot_gpu   = NULL;
 cudaMallocPitch((void**) &locpot_gpu,&un_pitch,largeur,height_max);

  
 /*------------ Local coordinates of points to calculate ---------------------*/
 int delta = pt0->x+ngfft[0]*(pt0->y+pt0->z*ngfft[1]); //-(virtual) linear initial point
 int final = pt1->x+ngfft[0]*(pt1->y+pt1->z*ngfft[1]); //-(virtual) linear final point
 int pth_size = (int)(un_pitch/sizeof(cureal));  //-Pitched size of vectors 
 int ntranche = (final-delta)+1;//-How many (virtual) pts to compute

 /*------------ Auxiliary Complex Vector -----------*/
 cucmplx* cvn_gpu  = NULL;
 cudaMallocPitch((void**) &cvn_gpu,&un_pitch,clargeur,height_max);
 int cvpthsz = int(un_pitch/sizeof(cucmplx)); //-Pitched size of complex auxiliar vec.
 printf("CUFFT PITCHED SIZES %d %d\n",pth_size,cvpthsz);

 /*------------ Copy Variable in cnst Memory -------------------*/
 copytoconstmem(nfftrec, nptrec, pth_size, cvpthsz);

 /*-------- Initialization points to compute in the First Loop ---*/
 int pos0 = 0;
 int loctranc = min(nptrec,npt);
 int pos1 = pos0+loctranc; 
 *max_rec = 0;


 //-Trotter
 int loctrott = 4*trotter;
 if(trotter == 0) loctrott = 2;

 /*-------- Arrays for Exit Criteria (Density Calculation) ---*/
 cucmplx* ND      = (cucmplx*) malloc(npt*loctrott*sizeof(cucmplx));
 cucmplx* NDold   = (cucmplx*) malloc(npt*loctrott*sizeof(cucmplx));
 cucmplx* NDnew   = (cucmplx*) malloc(npt*loctrott*sizeof(cucmplx));
 cucmplx* acc_rho = (cucmplx*) malloc(npt*sizeof(cucmplx));
 cureal*  erreur  = (cureal*)  malloc(2*npt*sizeof(cureal));
 cureal*  prod_b2 = (cureal*)  malloc(npt*sizeof(cureal));

 calctime(&stop,start,timing,0); //-End Timing: memory allocation+setting

 /*---------------------------------------------------------------------------------*/
 /*---------------- MAIN LOOP on pos1>pos0 -----------------------------------------*/
 /*---------------------------------------------------------------------------------*/
 printf(" Start  %10d\n End    %10d\n Npt    %10d\n gratio %10d\n Tot pt %10d\n",delta,final,npt,gratio,ntranche);

 int ipt;
 int npoint = delta;

 do{starttime(&start);
   int contrec = 0;
   printf("now: from %d to %d, so %d pts of %d \n",pos0+delta,pos1+delta,loctranc,npt);    
   int3 trasl;
   ipt = 0;
   for(int kk=pt0->z;kk<=pt1->z;kk+=gratio){
     trasl.z=kk-(ngfftrec[2]>>1);
     for(int jj=0;jj<ngfft[1];jj+=gratio){
       trasl.y=jj-(ngfftrec[1]>>1);
       for(int ii=0;ii<ngfft[0];ii+=gratio){
	 trasl.x=ii-(ngfftrec[0]>>1);
	 int ipoint = ii+ngfft[0]*(jj+kk*ngfft[1]);
	 if(ipoint<npoint) continue;
	 get_loc_potent <<< grid_pot,block_pot >>> (pot_gpu, locpot_gpu, trasl,ipt, ngfftrec[0],ngfft[0]);
	 //prt_dbg_arr(&(locpot_gpu[ipt*pth_size]),largeur,6,target-3,"potloc");
	 ipt++;
	 if(ipt==loctranc || ipoint == final) {npoint = ipoint+1; goto end_3_loop;}
       }}}
   end_3_loop:  calctime(&stop,start,timing,1);
   check_err(0);
    
   /*--------- Setting arrays Un,Unold on the Device --------------*/  
   starttime(&start);
   setting_un_cut <<< grid, block >>>(un_gpu,unold_gpu,vn_gpu,an_gpu, bn2_gpu,rsqrt(inf_vol),target);
   calctime(&stop,start,timing,1);
   //prt_dbg_arr(un_gpu,largeur,10,0,"un0");
   check_err(0);
  

   /*------------------ Loop on nrec ------------------------------*/
   int irec;
   for(irec=0; irec<nrec+1; irec++){
#ifdef  HAVE_GPU_CUDA_DEBUG
     printf("IREC------------%d\n",irec);
#endif
     starttime(&start);
     //prt_dbg_arr(un_gpu,largeur,6,target-3+pth_size,"un_gpu");
     un_x_pot_cut <<< grid, block >>> (cvn_gpu,un_gpu, locpot_gpu, loctranc);
     calctime(&stop,start,timing,2);
     //prt_dbg_arr(vn_gpu,largeur,6,target-3+pth_size,"vn=un*pot");

     /*-------------- Loop on loctranc: CONVOLUTION by FFT -----*/
      
     /*----- FFT -----*/
     starttime(&start);
     for( int ii=0; ii<loctranc; ii++)
     {
       cucmplx* a=&(cvn_gpu[ii*cvpthsz]);      
       FFTEXECC2C(plan_dir,a,a,CUFFT_FORWARD);
     };     
     calctime(&stop,start,timing,3);

     /*---- Moltiplication of the FFT with the Green Kernel ------*/
     starttime(&start); 
     complex_prod_tot <<< grid, block >>> (cvn_gpu, ZT_p_gpu,loctranc,nfftrec);
     calctime(&stop,start,timing,4);

     /*---- Inverse FFT -----*/
     starttime(&start); 
     for( int ii=0; ii<loctranc; ii++){
       cucmplx* a=&(cvn_gpu[ii*cvpthsz]);      
       FFTEXECC2C(plan_dir,a,a,CUFFT_INVERSE);
     }//-end loop:loctranc
     calctime(&stop,start,timing,3);

     /*-------------- Compute Vn = Pot*Vn -------------*/
     starttime(&start);
     vn_x_pot_dv_cut <<< grid, block >>>(cvn_gpu,vn_gpu, locpot_gpu,inf_vol/cureal(nfftrec), loctranc);
     calctime(&stop,start,timing,5);   
     //prt_dbg_arr(vn_gpu,largeur,6,target-3,"vn=vn*pot*infvol2");

     /*-------------- Compute An = Un*Vn -------------*/
     starttime(&start);
     scalarProdGPU<<< 128,256 >>>(an_gpu, vn_gpu, un_gpu, inf_vol);
     //-Copying An on the host
     cudaMemcpy(&(an[(irec)*npt+pos0]),an_gpu,(size_t)loctranc*sizeof(cureal),cudaMemcpyDeviceToHost);
     calctime(&stop,start,timing,6);
     //prt_dbg_arr(an_gpu,(height_max)*sizeof(cureal),10,0,"an");

   
     /*-------------- PREPARING NEXT ITERATION IN IREC -----------------*/
     if(irec<nrec){
       /*---------- Compute Un,Vn,Unold: Old to New -------*/
       starttime(&start);
       oldtonew <<< grid,block >>> (un_gpu,vn_gpu,unold_gpu,an_gpu,bn2_gpu,loctranc);
       calctime(&stop,start,timing,7);

       /*---------- Compute Bn = Un*Un ---------------*/
       starttime(&start);
       scalarProdGPU <<< 128,256 >>> (bn2_gpu, un_gpu, un_gpu, inf_vol);
       //-Copying An on the host
       cudaMemcpy(&(bn2[(irec+1)*npt+pos0]),bn2_gpu,loctranc*sizeof(cureal),cudaMemcpyDeviceToHost);
       calctime(&stop,start,timing,8);
       //prt_dbg_arr(bn2_gpu,(height_max)*sizeof(cureal),10,0,"bn2");
	
       /*---------- Compute Un = Un/Sqrt(Bn) ---------------*/
       starttime(&start);
       un_invsqrt_scale <<< grid,block >>> (un_gpu, bn2_gpu, loctranc);
       calctime(&stop,start,timing,9);
       //prt_dbg_arr(un_gpu,largeur,6,target-3+pth_size,"unnew rescaled");

 
       /*--------- Exit Criterium: Density and Error Calculations ----------*/
       starttime(&start);
       density_calc( beta*fermie, 2./inf_vol,tolrec,irec,trotter,npt,loctranc,
		     pos0,&(contrec), bn2, an,
		     erreur,prod_b2, 
		     acc_rho, ND, NDold, NDnew);
       calctime(&stop,start,timing,10);
       if(contrec==loctranc) break;
     }     
   }//End loop on nrec

   *max_rec = max(*max_rec,irec+1);     
    
    /*-------- Points to compute in the Next Loop ---*/
   loctranc = min(nptrec,npt-pos1);
   pos0 = pos1; 
   pos1 = min(pos0+loctranc,npt-1);
   
 }while(pos0 <pos1);

#if defined HAVE_GPU_CUDA3
 cudaThreadSynchronize();
#else
 cudaDeviceSynchronize();
#endif   

 /*--------------Free Memory on Device and Host----------*/
 free(ND);
 free(NDold);
 free(NDnew);
 free(acc_rho);
 free(erreur);
 free(prod_b2);

 cufftDestroy(plan_dir);
 cudaFree(cvn_gpu);
 cudaFree(ZT_p_gpu);
 cudaFree(pot_gpu);
 cudaFree(locpot_gpu);
 cudaFree(vn_gpu);
 cudaFree(unold_gpu);
 cudaFree(un_gpu);
 cudaFree(an_gpu);
 cudaFree(bn2_gpu);

 /*____________________Printing Time_____________________*/
 prt_device_timing(timing,DEBUGLEN);
 free(timing);
 cudaEventDestroy(start);
 cudaEventDestroy(stop);
 
 printf("\n--end--cudarec cut------ \n"); 
 return;
}


