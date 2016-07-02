// FUNCTION
// For a given wave-function |c>, get all projected scalars
// <p_lmn|c> where |p_lmn> are non-local projectors
//   With:
//   <p_lmn|c>=4pi/sqr(vol) (i)^l Sum_g[c(g).f_nl(g).Y_lm(g).exp(2pi.i.g.R)]

#include <stdio.h>
#include <stdlib.h>

//Constant conversion indices
__constant__ unsigned char indalpha_c3[]={0,1,2,3,2,1,1};
__constant__ unsigned char  indbeta_c3[]={0,1,2,3,3,3,2};

/*******************************************************************/
/**********                                      *******************/
/**********          kernel definition           *******************/
/**********                                      *******************/
/*******************************************************************/

//Note that we assume in a first time we have only one block in X direction,
// SO we hava gridDim.x = 1 and blockIdx.x=0
// Blocks in Y direction represente a couple (iaton,ilmn) for which we have
// to compute C_ilmn^iaton
//The primary compute and the reduction are made over the threads of this block
__global__ void kernel_compute_nl_projections(double2 *proj,double2 *vectin,double2 *ph3din,double *ffnlin,
					      const int *indlmn,const unsigned short int* atoms,
					      const unsigned char *lmn,const unsigned char * typat,
					      const int dimffnlin,const int npw,
					      const int lmnmax,const char cplex,const double four_pi_by_squcvol
					      ){
  
  //Definition of locals
  int ilmn,ipw,itypat;
  int dec_iatom_ph3d,dec_ffnl;
  double2 vect,ph3d;
  double ffnl;
  
  //Shared memory areas to compute and reduce
  extern __shared__ double sh_mem[];
  double *sh_proj_x = sh_mem ;
  double *sh_proj_y = &(sh_mem[blockDim.x]);
  
  //Get the couple (iatom,ilmn) of the block
  //iatom = atoms[blockIdx.y];// we don't need it because we use it only in dec_iatom_ph3d
  ilmn =lmn[blockIdx.y];
  itypat=typat[blockIdx.y];
  dec_iatom_ph3d = npw * atoms[blockIdx.y];
  dec_ffnl= npw*(0 + dimffnlin*(ilmn + lmnmax*itypat));
  
  //Initialisation of Shared mem
  sh_proj_x[threadIdx.x]=0.;
  sh_proj_y[threadIdx.x]=0.;
  
  //Step 1: Compute value for each plane wave and reduce by thread in sh mem
  for(ipw=threadIdx.x ;ipw<npw; ipw+=blockDim.x){
    //ffnl= ffnlin[ipw + npw*(0 + dimffnlin*(ilmn + lmnmax*itypat))];
    ffnl= ffnlin[ipw + dec_ffnl];
    vect= vectin[ipw];
    //ph3d= ph3din[ipw + npw*iatom];
    ph3d= ph3din[ipw + dec_iatom_ph3d];
    
    //We add the complex product in shared mem
    sh_proj_x[threadIdx.x] += ffnl*(vect.x*ph3d.x - vect.y*ph3d.y) ;
    sh_proj_y[threadIdx.x] += ffnl*(vect.y*ph3d.x + vect.x*ph3d.y) ;
  }
  
  //Step 2: Reduce between threads
  for(int decalage=blockDim.x>>1;decalage>0;decalage=decalage>>1){
    //We ensure every load on shared mem is accomplished in the block
    __syncthreads();
    if( threadIdx.x <  decalage) 
      sh_proj_x[threadIdx.x] += sh_proj_x[threadIdx.x + decalage];
    else if(threadIdx.x >= (blockDim.x - decalage))
      sh_proj_y[threadIdx.x] += sh_proj_y[threadIdx.x - decalage];
  }
  
  __syncthreads();
  
  //Write results for the block
  if(threadIdx.x == 0){
    //We need to know l;
    int l=indlmn[ 0 + 6*(ilmn + lmnmax*itypat)];
    //We write the result as a complex in register before storing it in glob mem
    switch (l%4){
    case 0: //i^l = 1
      vect.x = four_pi_by_squcvol * sh_proj_x[0] ;
      vect.y = four_pi_by_squcvol * sh_proj_y[blockDim.x -1] ;
      break;
    case 1: //i^l = i
      vect.x = -four_pi_by_squcvol * sh_proj_y[blockDim.x -1] ;
      vect.y = four_pi_by_squcvol * sh_proj_x[0] ;
      break;
    case 2: //i^l = -1
      vect.x = -four_pi_by_squcvol * sh_proj_x[0] ;
      vect.y = -four_pi_by_squcvol * sh_proj_y[blockDim.x -1] ;
      break;
    case 3: //i^l = -i
      vect.x = four_pi_by_squcvol * sh_proj_y[blockDim.x -1] ;
      vect.y = -four_pi_by_squcvol * sh_proj_x[0] ;
      break;
    default:
      break;
    }
    if(cplex==1){
      vect.x *= 2.;
      vect.y = 0.;
    }
    proj[blockIdx.y] = vect;
  }
  
}//end of kernel_compute_nl_projections


//Simultaneous computing of proj and Dproj
__global__ void kernel_compute_nl_projections_and_derivates_choice2(double2 *proj,double2 *dproj,
								    double2 *vectin,double2 *ph3din,
								    double *ffnlin,double *kpgin,
								    const int *indlmn,const unsigned short int* atoms,
								    const unsigned char *lmn,const unsigned char * typat,
								    const int dimffnlin,const int npw, const int nb_projections,
								    const int lmnmax,const char cplex,
								    const double four_pi_by_squcvol,const double height_pi2_by_squcvol
								    //,double *debug
								    ){
  
  //Definition of locals
  int iatom,ilmn,ipw,itypat;
  double2 vect,ph3d;//,proj_loc;
  double2 dproj0,dproj1,dproj2,dproj3;
  double ffnl,kpg1,kpg2,kpg3;
  //int idebug=0;
  //Shared memory areas to compute and reduce
  extern __shared__ double sh_mem[];
  double *sh_proj_x = sh_mem ;
  double *sh_proj_y = &(sh_mem[blockDim.x]);
  
  //Get the couple (iatom,ilmn) of the block
  iatom=atoms[blockIdx.y];
  ilmn =lmn[blockIdx.y];
  itypat=typat[blockIdx.y];
  
  //Initialisation of Shared mem
  sh_proj_x[threadIdx.x]=0.;
  sh_proj_y[threadIdx.x]=0.;
  dproj1.x=0.;dproj2.x=0.;dproj3.x=0.;
  dproj1.y=0.;dproj2.y=0.;dproj3.y=0.;
  
  //Step 1: Compute value for each plane wave and reduce by thread in sh mem
  for(ipw=threadIdx.x + blockIdx.x*gridDim.x;ipw<npw; ipw+=blockDim.x*gridDim.x){
    
    ffnl= ffnlin[ipw + npw*(0 + dimffnlin*(ilmn + lmnmax*itypat))];
    vect= vectin[ipw];
    ph3d= ph3din[ipw + npw*iatom];
    kpg1= kpgin[ipw + 0*npw];
    kpg2= kpgin[ipw + 1*npw];
    kpg3= kpgin[ipw + 2*npw];
    
    //We add the complex product in shared mem and derivates part in register
    dproj0.x = ffnl*(vect.x*ph3d.x - vect.y*ph3d.y);
    dproj0.y = ffnl*(vect.y*ph3d.x + vect.x*ph3d.y);
    sh_proj_x[threadIdx.x] += dproj0.x ;
    sh_proj_y[threadIdx.x] += dproj0.y ;
    dproj1.x  += kpg1*dproj0.x ;
    dproj1.y  += kpg1*dproj0.y ;
    dproj2.x  += kpg2*dproj0.x ;
    dproj2.y  += kpg2*dproj0.y ;
    dproj3.x  += kpg3*dproj0.x ;
    dproj3.y  += kpg3*dproj0.y ;
    //     if(blockIdx.y==0)
    //       debug[ipw]=-kpg1*dproj0.y;
  }
    
  //Step 2: Reduce between threads
  for(int decalage=blockDim.x>>1;decalage>0;decalage=decalage>>1){
    //We ensure every load on shared mem is accomplished in the block
    __syncthreads();
    if( threadIdx.x <  decalage) 
      sh_proj_x[threadIdx.x] += sh_proj_x[threadIdx.x + decalage];
    else if(threadIdx.x >= (blockDim.x - decalage))
      sh_proj_y[threadIdx.x] += sh_proj_y[threadIdx.x - decalage];
  }
  __syncthreads();
  //Store locally the sum pro projection
  if(threadIdx.x==0){
    dproj0.x=sh_proj_x[0] ;
    dproj0.y=sh_proj_y[blockDim.x -1] ;
  }
  
  //Step 3: Reduce for first derivate
  __syncthreads();
  sh_proj_x[threadIdx.x]= dproj1.x;
  sh_proj_y[threadIdx.x]= dproj1.y;;
  for(int decalage=blockDim.x>>1;decalage>0;decalage=decalage>>1){
    //We ensure every load on shared mem is accomplished in the block
    __syncthreads();
    if( threadIdx.x <  decalage) 
      sh_proj_x[threadIdx.x] += sh_proj_x[threadIdx.x + decalage];
    else if(threadIdx.x >= (blockDim.x - decalage))
      sh_proj_y[threadIdx.x] += sh_proj_y[threadIdx.x - decalage];
  }
  __syncthreads();
  
  //Store locally the sum for first derivate
  if(threadIdx.x==0){
    dproj1.x=sh_proj_x[0] ;
    dproj1.y=sh_proj_y[blockDim.x -1];
  }
  
  //Step 4: Reduce for second derivate
  __syncthreads();
  sh_proj_x[threadIdx.x]= dproj2.x;
  sh_proj_y[threadIdx.x]= dproj2.y;
  for(int decalage=blockDim.x>>1;decalage>0;decalage=decalage>>1){
    //We ensure every load on shared mem is accomplished in the block
    __syncthreads();
    if( threadIdx.x <  decalage) 
      sh_proj_x[threadIdx.x] += sh_proj_x[threadIdx.x + decalage];
    else if(threadIdx.x >= (blockDim.x - decalage))
      sh_proj_y[threadIdx.x] += sh_proj_y[threadIdx.x - decalage];
  }
  __syncthreads();
  
  //Store locally the sum for second derivate
  if(threadIdx.x==0){
    dproj2.x=sh_proj_x[0] ;
    dproj2.y=sh_proj_y[blockDim.x -1];
  }
  
  //Step 5: Reduce for third derivate
  __syncthreads();
  sh_proj_x[threadIdx.x]= dproj3.x;
  sh_proj_y[threadIdx.x]= dproj3.y;;
  for(int decalage=blockDim.x>>1;decalage>0;decalage=decalage>>1){
    //We ensure every load on shared mem is accomplished in the block
    __syncthreads();
    if( threadIdx.x <  decalage) 
      sh_proj_x[threadIdx.x] += sh_proj_x[threadIdx.x + decalage];
    else if(threadIdx.x >= (blockDim.x - decalage))
      sh_proj_y[threadIdx.x] += sh_proj_y[threadIdx.x - decalage];
  }
  __syncthreads();
  
   //Write results for the block
  //At this moment, proj,dproj1 and dproj2 are in register for tx0 and dproj3 in ShMem
  //We use the dproj3 varaible name to store tmp permutation
  if(threadIdx.x == 0){
    //We need to know l;
    int l=indlmn[ 0 + 6*(ilmn + lmnmax*itypat)];
    //We write the result as a complex in registers before storing it in glob mem
    double tmp;
    switch (l%4){
    case 0: //i^l = 1 ==> i^(l+1) = i 
      vect.x = four_pi_by_squcvol * dproj0.x;    
      vect.y = four_pi_by_squcvol * dproj0.y;
      
      tmp = dproj1.x;
      dproj1.x = -height_pi2_by_squcvol * dproj1.y;
      dproj1.y = height_pi2_by_squcvol * tmp;
      
      tmp = dproj2.x;
      dproj2.x = -height_pi2_by_squcvol * dproj2.y;
      dproj2.y = height_pi2_by_squcvol * tmp;
      
      dproj3.x = -height_pi2_by_squcvol *sh_proj_y[blockDim.x -1] ;
      dproj3.y = height_pi2_by_squcvol * sh_proj_x[0] ;
      break;
    case 1: //i^l = i ==> i^(l+1)=-1
      vect.x = -four_pi_by_squcvol * dproj0.y;    
      vect.y = four_pi_by_squcvol * dproj0.x;
      
      dproj1.x = -height_pi2_by_squcvol * dproj1.x;
      dproj1.y = -height_pi2_by_squcvol * dproj1.y;
      
      dproj2.x = -height_pi2_by_squcvol * dproj2.x;
      dproj2.y = -height_pi2_by_squcvol * dproj2.y;
      
      dproj3.x = -height_pi2_by_squcvol * sh_proj_x[0] ;
      dproj3.y = -height_pi2_by_squcvol * sh_proj_y[blockDim.x -1] ;
      break;
    case 2: //i^l = -1 ==> i^(l+1)=-i
      vect.x = -four_pi_by_squcvol * dproj0.x;    
      vect.y = -four_pi_by_squcvol * dproj0.y;
      
      tmp = dproj1.x;
      dproj1.x = height_pi2_by_squcvol * dproj1.y;
      dproj1.y = -height_pi2_by_squcvol * tmp;

      tmp = dproj2.x;
      dproj2.x = height_pi2_by_squcvol * dproj2.y;
      dproj2.y = -height_pi2_by_squcvol * tmp;
      
      dproj3.x = height_pi2_by_squcvol *sh_proj_y[blockDim.x -1] ;
      dproj3.y = -height_pi2_by_squcvol * sh_proj_x[0] ;
      break;
    case 3: //i^l = -i == i^(l+1)=1
      vect.x = four_pi_by_squcvol * dproj0.y;    
      vect.y = -four_pi_by_squcvol * dproj0.x;
      
      dproj1.x = height_pi2_by_squcvol * dproj1.x;
      dproj1.y = height_pi2_by_squcvol * dproj1.y;
      
      dproj2.x = height_pi2_by_squcvol * dproj2.x;
      dproj2.y = height_pi2_by_squcvol * dproj2.y;
      
      dproj3.x = height_pi2_by_squcvol * sh_proj_x[0] ;
      dproj3.y = height_pi2_by_squcvol * sh_proj_y[blockDim.x -1] ;
      break;
    default:
      break;
    }
    if(cplex==1){
      vect.x *= 2.;
      vect.y = 0.;
      dproj1.x *= 2;
      dproj1.y = 0.;
      dproj2.x *= 2;
      dproj2.y = 0.;
      dproj3.x *= 2;
      dproj3.y = 0.;
    }
    proj[blockIdx.y] = vect;
    dproj[blockIdx.y + 0*nb_projections] = dproj1;;
    dproj[blockIdx.y + 1*nb_projections] = dproj2;
    dproj[blockIdx.y + 2*nb_projections] = dproj3;
  }
  
}//end of kernel_compute_nl_projections_and_derivates_choice2


//Simultaneous computing of proj and Dproj second version 
//One thread's block compute one derivative of one projection
//A order 0 on derivation is the projection itself  
//Each block is indiced by (iproj,ialphabeta) where iproj is the projection's indice
//ialphabeta correspond to a unique couple (ialpha,ibeta) 
//We derivate ffnlin relatively to ialpha and kpg to ibeta
//The reduction over planewaves for each derivative is split among threads of the block
__global__ void kernel_compute_nl_projections_and_derivates_all(double2 *dproj,
								double2 *vectin,double2 *ph3din,
								double *ffnlin,double *kpgin,
								const int *indlmn,const unsigned short int* atoms,
								const unsigned char *lmn,const unsigned char * typat,
								const int dimffnlin,const int npw, const int nb_projections,
								const int lmnmax,const char cplex,const int choice,
								const double four_pi_by_squcvol,const double two_pi
								){
  
  //Definition of locals
  int iatom,ilmn,ipw,itypat;
  int iproj,ialphabeta;
  unsigned char ialpha,ibeta;
  double2 vect,ph3d;
  double ffnl_alpha,kpg_beta,fact_alphabeta;
  
  //Shared memory areas to compute and reduce
  extern __shared__ double sh_mem[];
  double *sh_proj_x = sh_mem ;
  double *sh_proj_y = &(sh_mem[blockDim.x]);
  
  //Get indices
  iproj=blockIdx.y;
  ialphabeta = blockIdx.x;

  //Compute alpha and beta from ialphabeta
  if(choice==2){
    ialpha=0;
    ibeta=ialphabeta;
  }
  else if(choice==3){
    ialpha = indalpha_c3[ialphabeta];
    ibeta  =  indbeta_c3[ialphabeta];
  }
  else if(choice==23){
    if(ialphabeta<4){
      ialpha=0;
      ibeta=ialphabeta;
    }
    else{
      ialpha = indalpha_c3[ialphabeta - 3];
      ibeta  =  indbeta_c3[ialphabeta - 3];
    }
  }
  
  //Get the couple (iatom,ilmn) of the block
  iatom=atoms[iproj];
  ilmn =lmn[iproj];
  itypat=typat[iproj];
  
  //Initialisation of Shared mem
  sh_proj_x[threadIdx.x]=0.;
  sh_proj_y[threadIdx.x]=0.;
  
  //Step 1: Compute value for each plane wave and reduce by thread in sh mem
  for(ipw=threadIdx.x ;ipw<npw; ipw+=blockDim.x){
    
    ffnl_alpha = ffnlin[ipw + npw*(ialpha + dimffnlin*(ilmn + lmnmax*itypat))];
    vect= vectin[ipw];
    ph3d= ph3din[ipw + npw*iatom];
    if(ibeta>0)
      kpg_beta = kpgin[ipw + (ibeta-1)*npw];
    else
      kpg_beta = 1.;
    
    fact_alphabeta = kpg_beta*ffnl_alpha; 
    
    if( ((choice==3)&&(ialphabeta>3)) || ((choice==23)&&(ialphabeta>6)) ){
      //We switch the alpha beta symetry
      ffnl_alpha = ffnlin[ipw + npw*(ibeta + dimffnlin*(ilmn + lmnmax*itypat))];
      kpg_beta = kpgin[ipw + (ialpha-1)*npw];
      fact_alphabeta = (fact_alphabeta + kpg_beta*ffnl_alpha)/2 ;
    }
    
    //We add the complex product in shared mem and derivates part in register
    sh_proj_x[threadIdx.x] += fact_alphabeta*(vect.x*ph3d.x - vect.y*ph3d.y);
    sh_proj_y[threadIdx.x] += fact_alphabeta*(vect.y*ph3d.x + vect.x*ph3d.y);
  }
  
  //Step 2: Reduce between threads
  for(int decalage=blockDim.x>>1;decalage>0;decalage=decalage>>1){
    //We ensure every load on shared mem is accomplished in the block
    __syncthreads();
    if( threadIdx.x <  decalage) 
      sh_proj_x[threadIdx.x] += sh_proj_x[threadIdx.x + decalage];
    else if(threadIdx.x >= (blockDim.x - decalage))
      sh_proj_y[threadIdx.x] += sh_proj_y[threadIdx.x - decalage];
  }
  __syncthreads();
  
  //Step 3: Write results for the block
  if(threadIdx.x == 0){
    //We need to know l;
    int l=indlmn[ 0 + 6*(ilmn + lmnmax*itypat)];
    if( ((ialphabeta>0)&&(choice==2)) || ((choice==23)&&(ialphabeta>0)&&(ialphabeta<4)) )
      l++;
    //We write the result as a complex in registers before storing it in glob mem
    switch (l%4){
    case 0: //i^l = 1 
      vect.x = four_pi_by_squcvol * sh_proj_x[0] ;
      vect.y = four_pi_by_squcvol * sh_proj_y[blockDim.x -1] ;
      break;
    case 1: //i^l = i
      vect.x = -four_pi_by_squcvol * sh_proj_y[blockDim.x -1] ;
      vect.y = four_pi_by_squcvol * sh_proj_x[0] ;
      break;
    case 2: //i^l = -1
      vect.x = -four_pi_by_squcvol * sh_proj_x[0] ;
      vect.y = -four_pi_by_squcvol * sh_proj_y[blockDim.x -1] ;
      break;
    case 3: //i^l = -i
      vect.x = four_pi_by_squcvol * sh_proj_y[blockDim.x -1] ;
      vect.y = -four_pi_by_squcvol * sh_proj_x[0] ;
      break;
    default:
      break;
    }
    
    if(cplex==1){
      vect.x *= 2.;
      vect.y = 0.;
    }
    if( ((ialphabeta>0)&&(choice==2)) || ((choice==23)&&(ialphabeta>0)&&(ialphabeta<4)) ){
      vect.x *= two_pi;
      vect.y *= two_pi;
    }
    if( ((ialphabeta>0)&&(choice==3)) || ((choice==23)&&(ialphabeta>3)) ){
      vect.x = -vect.x;
      vect.y = -vect.y;
    }
    dproj[iproj + ialphabeta*nb_projections] = vect;
  }
  
}//end of kernel_compute_nl_projections_and_derivates_all



/*******************************************************************/
/**********                                      *******************/
/**********          calling fonction            *******************/
/**********                                      *******************/
/*******************************************************************/

extern "C" void gpu_compute_nl_projections_(double2 *proj_gpu,double2 *dproj_gpu,
					    double2 *vectin_gpu,double2 *ph3din_gpu,
					    double *ffnlin_gpu,double *kpgin_gpu,
					    int *indlmn_gpu,unsigned short int *atoms_gpu,
					    unsigned char *lmn_gpu, unsigned char *typat_gpu,
					    int *nb_proj_to_compute,int *npw,int *choice,
					    int *dimffnlin,int *lmnmax,
					    const char *cplex,const double *pi,const double *ucvol
					    ){
  
  //Configuration of the cuda grid
  dim3 grid,block;
  double four_pi_by_squcvol = 4*(*pi)/sqrt(*ucvol);
  if((*choice)<2){
    //One block by projection to compute and plane waves split among a batch of threads by block   
    block.x=64;
    grid.x=1;
    grid.y=*nb_proj_to_compute;
    
    kernel_compute_nl_projections<<<grid,block,block.x*2*sizeof(double),0>>>(proj_gpu,vectin_gpu,ph3din_gpu,ffnlin_gpu,
									     indlmn_gpu,atoms_gpu,lmn_gpu,typat_gpu,
									     *dimffnlin,*npw,*lmnmax,
									     *cplex,four_pi_by_squcvol);
  }
  else if((*choice)==2){
    
    block.x=64;
    grid.x=1;
    grid.y=*nb_proj_to_compute;
    
    double height_pi2_by_squcvol  = 2*(*pi)*four_pi_by_squcvol; 
    
    kernel_compute_nl_projections_and_derivates_choice2<<<grid,block,block.x*2*sizeof(double),0>>>(proj_gpu,dproj_gpu,
												   vectin_gpu,ph3din_gpu,
												   ffnlin_gpu,kpgin_gpu,
												   indlmn_gpu,atoms_gpu,lmn_gpu,typat_gpu,
												   *dimffnlin,*npw,*nb_proj_to_compute,
												   *lmnmax,*cplex,
												   four_pi_by_squcvol,height_pi2_by_squcvol
												   );
  }
  else if( ((*choice)==3) || ((*choice)==23) ){
    
    double two_pi = 2*(*pi);
    block.x=64;
    
    //Grid.x = Nb derivates to compute (order 0 included)
    if((*choice)==3)
      grid.x = 7 ;
    else
      grid.x = 10;
    
    grid.y = *nb_proj_to_compute;
    
    kernel_compute_nl_projections_and_derivates_all<<<grid,block,block.x*2*sizeof(double),0>>>(dproj_gpu,
											       vectin_gpu,ph3din_gpu,
											       ffnlin_gpu,kpgin_gpu,
											       indlmn_gpu,atoms_gpu,lmn_gpu,typat_gpu,
											       *dimffnlin,*npw,*nb_proj_to_compute,
											       *lmnmax,*cplex,*choice,
											       four_pi_by_squcvol,two_pi
											       );
  }
}//end of routine gpu_compute_nl_projections_
