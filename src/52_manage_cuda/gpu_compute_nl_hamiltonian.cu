//{\src2tex{textfont=tt}}
//****f* ABINIT/gpu_compute_nl_hamiltonian
// NAME
// gpu_compute_nl_hamiltonian
//
// FUNCTION
//  Kernels definitions and calling functions for
//  the application of non-local hamiltonian
//
// COPYRIGHT
// Copyright (C) 1998-2019 ABINIT group (FDahm)
// This file is distributed under the terms of the
// GNU General Public License, see ~abinit/COPYING
// or http://www.gnu.org/copyleft/gpl.txt .
// For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
//
// INPUTS
//
// OUTPUT
//
// SIDE EFFECTS
//
// SOURCE

#include <stdio.h>
#include <stdlib.h>

/*******************************************************************/
/**********                                      *******************/
/**********         kernels definitions          *******************/
/**********                                      *******************/
/*******************************************************************/

__device__ static inline void reduce_complex_64(volatile double *sh_vect_x,volatile double *sh_vect_y){
  //We have only 2 warps, the first reduce the real part, the second the imaginary one
  //No more synchronization is needed
  if(threadIdx.x < 32)
    sh_vect_x[threadIdx.x] += sh_vect_x[threadIdx.x + 32];
  else//if(threadIdx.x > 31)
    sh_vect_y[threadIdx.x] += sh_vect_y[threadIdx.x - 32];

  if(threadIdx.x < 16)
    sh_vect_x[threadIdx.x] += sh_vect_x[threadIdx.x + 16];
  if(threadIdx.x > 47)
    sh_vect_y[threadIdx.x] += sh_vect_y[threadIdx.x - 16];

  if(threadIdx.x < 8)
    sh_vect_x[threadIdx.x] += sh_vect_x[threadIdx.x + 8];
  if(threadIdx.x > 55)
    sh_vect_y[threadIdx.x] += sh_vect_y[threadIdx.x - 8];

  if(threadIdx.x < 4)
    sh_vect_x[threadIdx.x] += sh_vect_x[threadIdx.x + 4];
  if(threadIdx.x > 59)
    sh_vect_y[threadIdx.x] += sh_vect_y[threadIdx.x - 4];

  if(threadIdx.x < 2)
    sh_vect_x[threadIdx.x] += sh_vect_x[threadIdx.x + 2];
  if(threadIdx.x > 61)
    sh_vect_y[threadIdx.x] += sh_vect_y[threadIdx.x - 2];

  if(threadIdx.x < 1)
    sh_vect_x[threadIdx.x] += sh_vect_x[threadIdx.x + 1];
  if(threadIdx.x > 62)
    sh_vect_y[threadIdx.x] += sh_vect_y[threadIdx.x - 1];
}

__device__ static inline  void reduce_double_32(volatile double *sh_vect){
  //We have only 1 active warp,
  //No more synchronization is needed
  if(threadIdx.x < 16)
    sh_vect[threadIdx.x] += sh_vect[threadIdx.x + 16];
  if(threadIdx.x < 8)
    sh_vect[threadIdx.x] += sh_vect[threadIdx.x + 8];
  if(threadIdx.x < 4)
    sh_vect[threadIdx.x] += sh_vect[threadIdx.x + 4];
  if(threadIdx.x < 2)
    sh_vect[threadIdx.x] += sh_vect[threadIdx.x + 2];
  if(threadIdx.x < 1)
    sh_vect[threadIdx.x] += sh_vect[threadIdx.x + 1];
}

__global__ void kernel_compute_proj_factor(double2 *proj,double2 *dproj,
					   double2 *val_ajlmn,double2 *val_sajlmn,
					   double *enl,double *sij,
					   double *rdlmn,double *rdproj,
					   const int *atindx1,const unsigned char *nlmn,
					   const int *indlmn,const unsigned short int* atoms,
					   const unsigned char *lmn,const unsigned char *typat,
					   const int paw_opt,const int dimenl1,const int lmnmax,
					   const int nb_projections,const int natoms,
					   const int choice, const int signs, const double lambda
					   ){

  int iproj=threadIdx.x + blockIdx.x*blockDim.x;
  double2 a_jlmn,sa_jlmn;
  double tmp_loc;
  unsigned short int iatom,jatom;
  unsigned char jlmn,itypat;

  if(iproj<nb_projections){

    //Get the couple (iatom,ilmn) of the thread's projection
    iatom=atoms[iproj];//atoms's indice sorted by type
    jatom=atindx1[iatom] - 1; //atoms's indice not sorted (-1 because of fortran cpu storage)
    jlmn =lmn[iproj];
    itypat=typat[iproj];
    int l=indlmn[ 0 + 6*(jlmn + lmnmax*itypat)];
    //Norm Conserving Pseudo Potential
    if(paw_opt==0){
      double2 proj_ilmn_R;
      double val_enl; //enl is real when norm conserving
      int iln=indlmn[ 4 + 6*(jlmn + lmnmax*itypat)]; // iln=indlmn(5,ilmn) (pour un type fixe)
      val_enl = enl[(iln-1) + dimenl1*itypat]; //enl_(1)=enl(iln,itypat,ispinor)
      proj_ilmn_R = proj[iproj];
      a_jlmn.x = val_enl * proj_ilmn_R.x ;//gxfac(1:cplex,ilmn,ia,ispinor)=enl_(1)*gx(1:cplex,ilmn,ia,ispinor)
      a_jlmn.y = val_enl * proj_ilmn_R.y ;
    }
    else{ //PAW
      a_jlmn.x=0.;
      a_jlmn.y=0.;
      sa_jlmn.x=0.;
      sa_jlmn.y=0.;
      //Accumulate projection;
      for(int ilmn=0;ilmn<nlmn[itypat];ilmn++){
	int ijlmn= (ilmn<jlmn? ilmn +  jlmn*(jlmn+1)/2 : jlmn +  ilmn*(ilmn+1)/2 );
	double val_enl = enl[ijlmn + dimenl1*jatom];
	double2 proj_ilmn_R = proj[ iproj + ilmn - jlmn];
	a_jlmn.x +=  val_enl*proj_ilmn_R.x;
	a_jlmn.y +=  val_enl*proj_ilmn_R.y;
	double val_sij= sij[ijlmn + dimenl1*itypat];
	sa_jlmn.x +=  val_sij*proj_ilmn_R.x;
	sa_jlmn.y +=  val_sij*proj_ilmn_R.y;
      }
    }
    if(paw_opt==2){
      a_jlmn.x -= lambda*sa_jlmn.x;
      a_jlmn.y -= lambda*sa_jlmn.y;
    }

    //Cas choice==1 && signs==2 :  a_jlmn = a_jlmn*i^(-l)
    if(choice==1){
      tmp_loc = a_jlmn.x;
      switch (l%4){
	//case 0: //i^(-l) = 1 : Nothing to do
      case 1: //i^(-l) = -i
	a_jlmn.x = a_jlmn.y;
	a_jlmn.y = -tmp_loc;
	tmp_loc = sa_jlmn.x;
	sa_jlmn.x = sa_jlmn.y;
	sa_jlmn.y = -tmp_loc;
	break;
      case 2: //i^(-l) = -1
	a_jlmn.x =  -a_jlmn.x ;
	a_jlmn.y =  -a_jlmn.y ;
	sa_jlmn.x = -sa_jlmn.x;
	sa_jlmn.y = -sa_jlmn.y;
	break;
      case 3: //i^(-l) = i
	a_jlmn.x = -a_jlmn.y;
	a_jlmn.y = tmp_loc;
	tmp_loc = sa_jlmn.x;
	sa_jlmn.x = -sa_jlmn.y;
	sa_jlmn.y = tmp_loc;
	break;
      default:
	break;
      }
      //Store computed values
      val_ajlmn[iproj]=a_jlmn;
      val_sajlmn[iproj]=sa_jlmn;
    }//End choice==1 && signs==2

    //Case choice==2 && signs==1: a_jlmn_alpha = a_jlmn*conjuguate(dproj(jlmn,alpha))
    else if(choice==2){
      if(paw_opt==3){
	a_jlmn.x = sa_jlmn.x;
	a_jlmn.y = sa_jlmn.y;
      }
      double2 val_dproj;
      //Alpha1 : ialpha=0
      val_dproj=dproj[iproj + 0*nb_projections];
      rdlmn[jlmn + lmnmax*(iatom + 0*natoms)] = val_dproj.x*a_jlmn.x + val_dproj.y*a_jlmn.y;
      //Alpha2 : ialpha=1
      val_dproj=dproj[iproj + 1*nb_projections];
      rdlmn[jlmn + lmnmax*(iatom + 1*natoms)] = val_dproj.x*a_jlmn.x + val_dproj.y*a_jlmn.y;
      //Alpha3 : ialpha=2
      val_dproj=dproj[iproj + 2*nb_projections];
      rdlmn[jlmn + lmnmax*(iatom + 2*natoms)] = val_dproj.x*a_jlmn.x + val_dproj.y*a_jlmn.y;
    }

    //Cas choice==3 && signs==1: a_jlmn_R_alphabeta = a_jlmn_R*conjuguate(dproj(iproj,alphabeta))
    else if(choice==3){
      if(paw_opt==3){
	a_jlmn.x = sa_jlmn.x;
	a_jlmn.y = sa_jlmn.y;
      }
      double2 val_dproj;
      for(int ialphabeta=0;ialphabeta<7;ialphabeta++){
	//ialphabeta=0
	val_dproj=dproj[iproj + ialphabeta*nb_projections];
	rdproj[iproj + ialphabeta*nb_projections] = val_dproj.x*a_jlmn.x + val_dproj.y*a_jlmn.y;
      }
    }
    //Case choice==23: both choice 2 and choice 3
    else if(choice==23){
      if(paw_opt==3){
	a_jlmn.x = sa_jlmn.x;
	a_jlmn.y = sa_jlmn.y;
      }
      //Choice 2 part
      double2 val_dproj;
      //Alpha1 : ialpha=0
      val_dproj=dproj[iproj + 1*nb_projections];
      rdlmn[jlmn + lmnmax*(iatom + 0*natoms)] = val_dproj.x*a_jlmn.x + val_dproj.y*a_jlmn.y;
      //Alpha2 : ialpha=1
      val_dproj=dproj[iproj + 2*nb_projections];
      rdlmn[jlmn + lmnmax*(iatom + 1*natoms)] = val_dproj.x*a_jlmn.x + val_dproj.y*a_jlmn.y;
      //Alpha3 : ialpha=2
      val_dproj=dproj[iproj + 3*nb_projections];
      rdlmn[jlmn + lmnmax*(iatom + 2*natoms)] = val_dproj.x*a_jlmn.x + val_dproj.y*a_jlmn.y;


      //choice 3 part
      //ialphabeta =0
      val_dproj=dproj[iproj];
      rdproj[iproj] = val_dproj.x*a_jlmn.x + val_dproj.y*a_jlmn.y;
      //others
      for(int ialphabeta=1;ialphabeta<7;ialphabeta++){
	val_dproj=dproj[iproj + (ialphabeta+3)*nb_projections];
	rdproj[iproj + ialphabeta*nb_projections] = val_dproj.x*a_jlmn.x + val_dproj.y*a_jlmn.y;
      }
    }//end choice 23
  }//end if thread's global id < nb_projection
}//End of kernel_compute_proj_factor


//Note that we assume in a first time we have only one block in X direction,
// So we have gridDim.x = 1 and blockIdx.x=0
// Blocks in Y direction represent a plane wave couple indexed by jpw
// threads of the blocks take care of a couple (iaton,ilmn)
//The primary compute and the reduction are made over the threads of this block
__global__ void kernel_compute_nl_hamiltonian(double2 *vectin,double2 *vectout,double2 *svectout,
					      double2 *val_ajlmn,double2 *val_sajlmn,
					      double2 *ph3dout,double *ffnlout,
					      const unsigned short int* atoms,
					      const unsigned char *lmn,const unsigned char *typat,
					      const int paw_opt,const int dimffnlout,const int npwout,
					      const int nb_projections,const int lmnmax,
					      const double four_pi_by_ucvol,const double lambda
					      ){
  //Definition of locals
  unsigned short int jpw,iatom;
  unsigned char jlmn,itypat;
  double2 vect_loc,svect_loc;

  //Shared memory areas to compute and reduce
  extern __shared__ double sh_mem[];
  double *sh_vect_x = sh_mem ;
  double *sh_vect_y = &(sh_mem[blockDim.x]);

  jpw=blockIdx.y;

  vect_loc.x = 0.;
  vect_loc.y = 0.;
  svect_loc.x = 0.;
  svect_loc.y = 0.;

  //Step 1: Compute value for each plane wave and reduce by thread in sh mem
  for(int iproj=threadIdx.x ; iproj<nb_projections ; iproj+=blockDim.x){
    double2 a_jlmn,sa_jlmn,ph3d;
    double ffnl;

    //Get the couple (iatom,ilmn) of the thread's projection
    iatom=atoms[iproj];//atoms's indice sorted by type
    jlmn =lmn[iproj];
    itypat=typat[iproj];

    //Read the projections factor
    a_jlmn =val_ajlmn[iproj];
    sa_jlmn=val_sajlmn[iproj];

    //Accumulate the contribution of the projection for this thread in register
    ffnl = ffnlout[jpw + npwout*(0 + dimffnlout*(jlmn + lmnmax*itypat))]; //ffnlout(npwout,dimffnlout,lmnmax,ntypat)
    ph3d = ph3dout[jpw + npwout*iatom];

    ph3d.x *= ffnl;
    ph3d.y *= ffnl;

    //Warning: the product is between a_jlmn and conjuguate(ph3d)
    vect_loc.x += a_jlmn.x*ph3d.x + a_jlmn.y*ph3d.y;
    vect_loc.y += a_jlmn.y*ph3d.x - a_jlmn.x*ph3d.y;
    svect_loc.x += sa_jlmn.x*ph3d.x + sa_jlmn.y*ph3d.y;
    svect_loc.y += sa_jlmn.y*ph3d.x - sa_jlmn.x*ph3d.y;

  }//End loop to performs all projections

  if(paw_opt!=3){
    //Step2: We reduce in Shared Memory
    sh_vect_x[threadIdx.x]=vect_loc.x;
    sh_vect_y[threadIdx.x]=vect_loc.y;
    for(int decalage=blockDim.x>>1;decalage>0;decalage=decalage>>1){
      //We ensure every access in shared mem is accomplished in the block
      __syncthreads();
      if( threadIdx.x <  decalage)
	sh_vect_x[threadIdx.x] += sh_vect_x[threadIdx.x + decalage];
      else if(threadIdx.x >= (blockDim.x - decalage))
	sh_vect_y[threadIdx.x] += sh_vect_y[threadIdx.x - decalage];
    }
    __syncthreads();
    //Step3: Thread 0 writes the results for the block (ie for the plane wave)
    if(threadIdx.x==0){
      if(paw_opt==2){
	vect_loc.x = -lambda*vectin[jpw].x + four_pi_by_ucvol * sh_vect_x[0];
	vect_loc.x = -lambda*vectin[jpw].y + four_pi_by_ucvol * sh_vect_y[blockDim.x-1];
      }
      else{
	vect_loc.x = four_pi_by_ucvol * sh_vect_x[0];
	vect_loc.y = four_pi_by_ucvol * sh_vect_y[blockDim.x-1];
      }
      vectout[jpw] = vect_loc;
    }
  }

  if(paw_opt>2){
    //Step4: We reduce for overlap
    __syncthreads();
    sh_vect_x[threadIdx.x]=svect_loc.x;
    sh_vect_y[threadIdx.x]=svect_loc.y;
    for(int decalage=blockDim.x>>1;decalage>0;decalage=decalage>>1){
      //We ensure every access in shared mem is accomplished in the block
      __syncthreads();
      if( threadIdx.x <  decalage)
	sh_vect_x[threadIdx.x] += sh_vect_x[threadIdx.x + decalage];
      else if(threadIdx.x >= (blockDim.x - decalage))
	sh_vect_y[threadIdx.x] += sh_vect_y[threadIdx.x - decalage];
    }
    __syncthreads();
    //Step5: Thread 0 writes the results for the block (ie for the plane wave)
    if(threadIdx.x==0){
      svect_loc = vectin[jpw];
      svect_loc.x += four_pi_by_ucvol * sh_vect_x[0];
      svect_loc.y += four_pi_by_ucvol * sh_vect_y[blockDim.x-1];
      svectout[jpw] = svect_loc;
    }
  }//End of overlap calculation

}//end of kernel_compute_nl_hamiltonian

__global__ void kernel_compute_nl_hamiltonian_64(double2 *vectin,double2 *vectout,double2 *svectout,
						 double2 *val_ajlmn,double2 *val_sajlmn,
						 double2 *ph3dout,double *ffnlout,
						 const unsigned short int* atoms,
						 const unsigned char *lmn,const unsigned char *typat,
						 const int paw_opt,const int dimffnlout,const int npwout,
						 const int nb_projections,const int lmnmax,
						 const double four_pi_by_ucvol,const double lambda
						 ){

  //Definition of locals
  unsigned short int jpw,iatom;
  unsigned char jlmn,itypat;

  double2 vect_loc,svect_loc;

  //Shared memory areas to compute and reduce
  extern __shared__ double sh_mem[];
  double *sh_vect_x = sh_mem ;
  double *sh_vect_y = &(sh_mem[64]);

  jpw=blockIdx.y;

  vect_loc.x = 0.;
  vect_loc.y = 0.;
  svect_loc.x = 0.;
  svect_loc.y = 0.;

  //Step 1: Compute value for each plane wave and reduce by thread in sh mem
  for(int iproj=threadIdx.x ; iproj<nb_projections ; iproj+=64){
    double2 a_jlmn,sa_jlmn,ph3d;
    double ffnl;

    //Get the couple (iatom,ilmn) of the thread's projection
    iatom=atoms[iproj];//atoms's indice sorted by type
    jlmn =lmn[iproj];
    itypat=typat[iproj];

    //Read the projections factor
    a_jlmn =val_ajlmn[iproj];
    sa_jlmn=val_sajlmn[iproj];

    //Accumulate the contribution of the projection for this thread in register
    //These 2 loads are'nt coalesced and may be costly
    ffnl = ffnlout[jpw + npwout*(0 + dimffnlout*(jlmn + lmnmax*itypat))]; //ffnlout(npwout,dimffnlout,lmnmax,ntypat)
    ph3d = ph3dout[jpw + npwout*iatom];

    ph3d.x *= ffnl;
    ph3d.y *= ffnl;

    //Warning: the product is between a_jlmn and conjuguate(ph3d)
    vect_loc.x += a_jlmn.x*ph3d.x + a_jlmn.y*ph3d.y;
    vect_loc.y += a_jlmn.y*ph3d.x - a_jlmn.x*ph3d.y;
    svect_loc.x += sa_jlmn.x*ph3d.x + sa_jlmn.y*ph3d.y;
    svect_loc.y += sa_jlmn.y*ph3d.x - sa_jlmn.x*ph3d.y;

  }//End loop to performs all projections

  //if(paw_opt!=3){
    //Step2: We reduce in Shared Memory
    sh_vect_x[threadIdx.x]=vect_loc.x;
    sh_vect_y[threadIdx.x]=vect_loc.y;
    __syncthreads();
    reduce_complex_64(sh_vect_x,sh_vect_y);
    __syncthreads();
    //Step3: Thread 0 writes the results for the block (ie for the plane wave)
    if(threadIdx.x==0){
      if(paw_opt==2){
	vect_loc.x = -lambda*vectin[jpw].x + four_pi_by_ucvol * sh_vect_x[0];
	vect_loc.x = -lambda*vectin[jpw].y + four_pi_by_ucvol * sh_vect_y[63];
      }
      else{
	vect_loc.x = four_pi_by_ucvol * sh_vect_x[0];
	vect_loc.y = four_pi_by_ucvol * sh_vect_y[63];
      }
      vectout[jpw] = vect_loc;
    }
    //}

    //  if(paw_opt>2){
    //Step4: We reduce for overlap
    __syncthreads();
    sh_vect_x[threadIdx.x]=svect_loc.x;
    sh_vect_y[threadIdx.x]=svect_loc.y;
     __syncthreads();
     reduce_complex_64(sh_vect_x,sh_vect_y);
     __syncthreads();

    //Step5: Thread 0 writes the results for the block (ie for the plane wave)
    if(threadIdx.x==0){
      svect_loc = vectin[jpw];
      svect_loc.x += four_pi_by_ucvol * sh_vect_x[0];
      svect_loc.y += four_pi_by_ucvol * sh_vect_y[63];
      svectout[jpw] = svect_loc;
    }
    //  }//End of overlap calculation

}//end of kernel_compute_nl_hamiltonian_64


//Accumulate previously calculate factors to get forces
__global__ void kernel_compute_forces(double *rdlmn,double *enlout,
				      const int decalage_enlout,
				      const int natoms,const int lmnmax
				      ){

  extern __shared__ double enl_sh_sum[];
  enl_sh_sum[threadIdx.x] = 0.;
  int ialpha=blockIdx.x;
  int iatom=blockIdx.y;

  //Load in shared memory
  for(int jlmn=threadIdx.x; jlmn<lmnmax; jlmn+=blockDim.x)
    //rdlmn is set to 0 if jlmn >= nlmn(type(iatoms))
    enl_sh_sum[threadIdx.x] += rdlmn[jlmn + lmnmax*(iatom + ialpha*natoms)];

  //Reduce in shared memory
  __syncthreads();
  reduce_double_32(enl_sh_sum);
  // for(int decalage=blockDim.x/2;decalage>0;decalage=decalage/2){
  //     __syncthreads();
  //     if(threadIdx.x < decalage)
  //       enl_sh_sum[threadIdx.x] += enl_sh_sum[threadIdx.x + decalage];
  //   }

  //Write results
  if(threadIdx.x==0)
    enlout[decalage_enlout + ialpha + 3*iatom] = 2*enl_sh_sum[0];

}//end of kernel_compute_forces


//Compute the stress tensor in reduced coordinates
__global__ void kernel_reduce_stress(double *rdproj,double *d_enlk,
				     const int nb_projections){
  extern __shared__ double sh_mem[];
  int ialphabeta=blockIdx.x;

  sh_mem[threadIdx.x]=0.;

  //Acumulate in Shared Mem
  for(int iproj= threadIdx.x; iproj < nb_projections ; iproj+=blockDim.x){
    sh_mem[threadIdx.x] += rdproj[iproj + ialphabeta*nb_projections];
  }

  //Reduce in shared memory
  for(int decalage=blockDim.x/2;decalage>0;decalage=decalage/2){
    __syncthreads();
    if(threadIdx.x < decalage)
      sh_mem[threadIdx.x] += sh_mem[threadIdx.x + decalage];
  }

  //Write results
  if(threadIdx.x==0)
    d_enlk[ialphabeta] = sh_mem[0];
}

__constant__ unsigned char imunu[] = {0,1,2,3,4,5,3,4,5};
__constant__ unsigned char imu[] = {0,1,2,1,0,0,2,2,1};
__constant__ unsigned char inu[] = {0,1,2,2,2,1,1,0,0};
__constant__ unsigned char ialpha[] = {0,1,2,1,0,0};
__constant__ unsigned char ibeta[]= {0,1,2,2,2,1};

//Convert stress tensor to cartesian coordinates
__global__ void kernel_stress_convert(double *gprimd,double *d_enlk,double *enlout){

  __shared__ double sh_gprim[9];
  __shared__ double sh_stress[6];
  __shared__ double sh_enl[9];

  sh_gprim[threadIdx.x] = gprimd[threadIdx.x];

  if(threadIdx.x < 6)
    sh_stress[threadIdx.x] = d_enlk[threadIdx.x+1];

  sh_enl[threadIdx.x] = sh_stress[imunu[threadIdx.x]] * sh_gprim[ialpha[blockIdx.x] + 3*imu[threadIdx.x]] * sh_gprim[ibeta[blockIdx.x] + 3*inu[threadIdx.x]];

  __syncthreads();

  if(threadIdx.x < 3 )
    sh_enl[threadIdx.x] += sh_enl[threadIdx.x + 3] + sh_enl[threadIdx.x + 6];

  if(threadIdx.x == 0 )
    sh_enl[0] += sh_enl[1] + sh_enl[2];

  if(threadIdx.x==0){
    if(blockIdx.x < 3 )
      enlout[blockIdx.x] = 2*sh_enl[0] - d_enlk[0];
    else
      enlout[blockIdx.x] = 2*sh_enl[0];
  }
}


/*******************************************************************/
/**********                                      *******************/
/**********          calling fonction            *******************/
/**********                                      *******************/
/*******************************************************************/

extern "C" void gpu_compute_nl_hamiltonian_(double2 *proj_gpu,double2 *dproj_gpu,
					    double2 *vectin_gpu,
					    double2 *vectout_gpu,double2 *svectout_gpu,
					    double2 *val_ajlmn_gpu,double2 *val_sajlmn_gpu,
					    double2 *ph3dout_gpu,double *ffnlout_gpu,
					    double *rdlmn_gpu, double *rdproj_gpu,
					    double *enlout_gpu, double *d_enlk_gpu,
					    double *enl_gpu,double *sij_gpu,double *gprimd_gpu,
					    int *atindx1_gpu,unsigned char *nlmn_gpu,
					    int *indlmn_gpu,unsigned short int *atoms_gpu,
					    unsigned char *lmn_gpu, unsigned char *typat_gpu,
					    int *paw_opt,int *dimffnlout,int *dimenl1,
					    int *npw,int *nb_projections,int *lmnmax,
					    int *natoms, int *choice, int *signs,
					    const double *four_pi_by_ucvol,const double *lambda
					    ){

  /************* Accumulate the projection factors ***********************/
  //Configuration of the cuda grid
  dim3 grid,block;
  block.x= 64;
  grid.x=((*nb_projections) + 64 - 1)/64 ;

  double2 *effective_proj_gpu=proj_gpu;
  //If choice = 3 or 23 proj is stored in dproj
  if( ((*choice)==3) || ((*choice)==23) )
    effective_proj_gpu=dproj_gpu;

  kernel_compute_proj_factor<<<grid,block>>>(effective_proj_gpu,dproj_gpu,
					     val_ajlmn_gpu,val_sajlmn_gpu,
					     enl_gpu,sij_gpu,
					     rdlmn_gpu,rdproj_gpu,
					     atindx1_gpu,nlmn_gpu,
					     indlmn_gpu,atoms_gpu,lmn_gpu,typat_gpu,
					     *paw_opt,*dimenl1,*lmnmax,
					     *nb_projections,*natoms,
					     *choice,*signs,*lambda
					     );

  /************* Compute output needed with signs ***********************/

  if( ((*choice)==1) && (*signs==2)) {
    //One thread in block by projection and plane waves split amongs blocks
    block.x = 64;
    grid.x = 1;
    grid.y=*npw;

    //kernel_compute_nl_hamiltonian<<<grid,block,block.x*2*sizeof(double),0>>>(vectin_gpu,vectout_gpu,svectout_gpu,
    kernel_compute_nl_hamiltonian_64<<<grid,block,block.x*2*sizeof(double),0>>>(vectin_gpu,vectout_gpu,svectout_gpu,
									       val_ajlmn_gpu,val_sajlmn_gpu,
									       ph3dout_gpu,ffnlout_gpu,
									       atoms_gpu,lmn_gpu,typat_gpu,
									       *paw_opt,*dimffnlout,*npw,
									       *nb_projections,*lmnmax,
									       *four_pi_by_ucvol,*lambda
									       );

  }//End choice==1,signs==2

  else {
    //Compute forces
    if( (((*choice)==2)&&(*signs==1)) || ((*choice)==23) ) {
      //Threads In blocks will take care of lmnmax
      block.x= 32;

      grid.x= 3; //Nb alpha directions
      grid.y= *natoms ;

      int decalage_enlout_gpu = 0;
      if( (*choice)==23)
	decalage_enlout_gpu = 6;

      kernel_compute_forces<<<grid,block,block.x*sizeof(double),0>>>(rdlmn_gpu,enlout_gpu,
								     decalage_enlout_gpu,
								     *natoms,*lmnmax);
    }//End choice2 / signs1 or choice 23

    //Compute stress tensor
    if( (((*choice)==3) && ((*signs)==1))  || ((*choice)==23) ) {

      grid.x= 7; //Nb alphabeta directions
      block.x= 512; // Enough threads to reduce on projection
      while(block.x > (unsigned int)(*nb_projections)){
	block.x /= 2;
      }
      block.x=64;

      kernel_reduce_stress<<<grid,block,block.x*sizeof(double),0>>>(rdproj_gpu,d_enlk_gpu,*nb_projections);

      kernel_stress_convert<<<6,9>>>(gprimd_gpu,d_enlk_gpu,enlout_gpu);

    }

  }//end choice!=1
}//End of gpu_compute_nl_hamiltonian

//***
