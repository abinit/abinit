//{\src2tex{textfont=tt}}
//****f* ABINIT/gpu_nonlop
//
// NAME
// gpu_nonlop
//
// FUNCTION
// This routine is a driver to compute:
// * Application of a nonlocal operator Vnl in order to get:
//    - contracted elements (energy, forces, stresses, ...), if signs=1
//    - a function in reciprocal space (|out> = Vnl|in>), if signs=2
// * Optionally, in case of PAW calculation:
//   - Application of the overlap matrix in reciprocal space
//     (<in|S|in> or (I+S)|in>).
//   - Application of (Vnl-lambda.S) in reciprocal space
//     (<in|Vnl-lambda.S|in> and derivatives or (Vnl-lambda.S)|in>).
// According to user s choice, the routine calls a subroutine, computing all quantities:
//   - using Legendre Polynomials Pl (Norm-conserving psps only)
//   - using Spherical Harmonics Ylm (N-conserving or PAW ; compulsory for PAW)
//
// COPYRIGHT
// Copyright (C) 1998-2019 ABINIT group (MT,FDahm)
// This file is distributed under the terms of the
// GNU General Public License, see ~abinit/COPYING
// or http://www.gnu.org/copyleft/gpl.txt .
// For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
//
// INPUTS
//  atindx1(natom)=index table for atoms, inverse of atindx
//  choice: chooses possible output:
//    choice=0 => do nothing (only compute WF projected with NL projectors)
//          =1 => a non-local energy contribution
//          =2 => a gradient with respect to atomic position(s)
//          =3 => a gradient with respect to strain(s)
//          =23=> a gradient with respect to atm. pos. and strain(s)
//          =4 => a gradient and 2nd derivative with respect to atomic pos.
//          =24=> a gradient and 2nd derivative with respect to atomic pos.
//          =5 => a gradient with respect to k wavevector, typically
//                $\sum_{ij}|p_i\rangle D_{ij}\langle dp_j/dk| + |dp_i/dk\rangle D_{ij}\langle p_j|$
//          =51 => the right derivative with respect to k wavevector, typically
//                $\sum_{ij}|p_i\rangle D_{ij}\langle dp_j/dk|$
//          =52 => the left derivative with respect to k wavevector, typically
//                $\sum_{ij}|dp_i/dk\rangle D_{ij}\langle p_j|$
//          =53 => the twist derivative with respect to k, typically
//                $\sum_{ij}|dp_i/dk_(idir+1)\rangle D_{ij}\langle dp_j/dk_(idir-1)| -
//                 |dp_i/dk_(idir-1)\rangle D_{ij}\langle dp_j/dk_(idir+1)|$
//          =6 => 2nd derivatives with respect to strain and atm. pos.
//  cpopt=flag defining the status of cprjin%cp(:)=<Proj_i|Cnk> scalars (see below, side effects)
//  dimenl1,dimenl2=dimensions of enl (see enl)
//  dimffnlin=second dimension of ffnlin (1+number of derivatives)
//  dimffnlout=second dimension of ffnlout (1+number of derivatives)
//  enl(dimenl1,dimenl2,nspinor**2)=
//  ->Norm conserving : ==== when paw_opt=0 ====
//                      (Real) Kleinman-Bylander energies (hartree)
//                      dimenl1=lmnmax  -  dimenl2=ntypat
//  ->PAW :             ==== when paw_opt=1 or 4 ====
//                      (Real, symmetric) Dij coefs to connect projectors
//                      dimenl1=cplex_enl*lmnmax*(lmnmax+1)/2  -  dimenl2=natom
//                      These are complex numbers if cplex_enl=2
//                        enl(:,:,1) contains Dij^up-up
//                        enl(:,:,2) contains Dij^dn-dn
//                        enl(:,:,3) contains Dij^up-dn (only if nspinor=2)
//                        enl(:,:,4) contains Dij^dn-up (only if nspinor=2)
//  ffnlin(npwin,dimffnlin,lmnmax,ntypat)=nonlocal form factors to be used
//          for the application of the nonlocal operator to the |in> vector
//  ffnlout(npwout,dimffnlout,lmnmax,ntypat)=nonlocal form factors to be used
//          for the application of the nonlocal operator to the |out> vector
//  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
//  gprimd(3,3)=dimensional reciprocal space primitive translations
//  idir=direction of the - atom to be moved in the case (choice=2,signs=2),
//                        - k point direction in the case (choice=5,51,or 52 and signs=2)
//                          for choice 53 signs=2, cross derivatives are in idir-1 and idir+1 directions
//                        - strain component (1:6) in the case (choice=3,signs=2) or (choice=6,signs=1)
//  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
//                                                  or i=lmn (if useylm=1)
//  istwf_k=option parameter that describes the storage of wfs
//  kgin(3,npwin)=integer coords of planewaves in basis sphere, for the |in> vector
//  kgout(3,npwout)=integer coords of planewaves in basis sphere, for the |out> vector
//  kpgin(npw,npkgin)= (k+G) components and related data, for the |in> vector  (only if useylm=1)
//  kpgout(npw,nkpgout)=(k+G) components and related data, for the |out> vector (only if useylm=1)
//  kptin(3)=k point in terms of recip. translations, for the |in> vector
//  kptout(3)=k point in terms of recip. translations, for the |out> vector
//  lambda=factor to be used when computing (Vln-lambda.S) - only for paw_opt=2
//         Typically lambda is the eigenvalue (or its guess)
//  lmnmax=max. number of (l,m,n) components over all types of atoms
//  matblk=dimension of the arrays ph3din and ph3dout
//  mgfft=maximum size of 1D FFTs
//  mpi_enreg=informations about MPI parallelization
//  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
//  mpssoang= 1+max(spin*angular momentum) for nonlocal pseudopotentials
//  natom=number of atoms in cell
//  nattyp(ntypat)=number of atoms of each type
//  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
//  nkpgin,nkpgout=second sizes of arrays kpgin/kpgout
//  nloalg(5)=governs the choice of the algorithm for nonlocal operator
//  nnlout=dimension of enlout (when signs=1 and choice>0):
//         ==== if paw_opt=0, 1 or 2 ====
//         choice=1=>nnlout=1           choice=2=>nnlout=3*natom
//         choice=3=>nnlout=6           choice=4=>nnlout=6*natom
//         choice=5=>nnlout=3           choice=6=>nnlout=6*(3*natom+6)
//         choice=23=>nnlout=6+3*natom  choice=24=>nnlout=9*natom
//         ==== if paw_opt=3 ====
//         choice=1 =>nnlout=1
//         ==== if paw_opt=4 ====
//         not available
//  npwin=number of planewaves for given k point, for the |in> vector
//  npwout=number of planewaves for given k point, for the |out> vector
//  nspinor=number of spinorial components of the wavefunctions
//  ntypat=number of types of atoms in cell
//  only_SO=flag to calculate only the SO part in nonlop
//  paw_opt= define the nonlocal operator concerned with:
//           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
//           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
//           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
//           paw_opt=3 : PAW overlap matrix (Sij)
//           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
//  phkxredin(2,natom)=phase factors exp(2 pi kptin.xred)
//  phkxredout(2,natom)=phase factors exp(2 pi kptout.xred)
//  ph1d(2,3*(2*mgfft+1)*natom)=1D structure factors phase information
//  ph3din(2,npwin,matblk)=3D structure factors, for each atom and plane wave (in)
//  ph3dout(2,npwout,matblk)=3-dim structure factors, for each atom and plane wave (out)
//  ----- Removed in beautification because unused MS -----------
//  pspso(ntypat)=spin-orbit characteristic for each atom type
//  -------------------------------------------------------------
//  signs= if 1, get contracted elements (energy, forces, stress, ...)
//         if 2, applies the non-local operator to a function in reciprocal space
//  sij(dimenl1,ntypat)=overlap matrix components (only if paw_opt=2, 3 or 4)
//  tim_nonlop=timing code of the calling routine (can be set to 0 if not attributed)
//  ucvol=unit cell volume (bohr^3)
//  useylm=governs the way the nonlocal operator is to be applied:
//         1=using Ylm, 0=using Legendre polynomials
//  vectin(2,npwin*nspinor)=input cmplx wavefunction coefficients <G|Cnk>
//
// OUTPUT
//  ==== if (signs==1) ====
// --If (paw_opt==0, 1 or 2)
//    enlout(nnlout)= contribution to the non-local part of the following properties:
//      if choice=1 : enlout(1)               -> the energy
//      if choice=2 : enlout(1:3*natom)       -> the forces
//      if choice=3 : enlout(1:6)             -> the stresses
//      if choice=23: enlout(1:6+3*natom)     -> the forces and the stresses
//      if choice=4 : enlout(1:6*natom)       -> the frozen wf part of dyn. mat.
//      if choice=5 : enlout(3)               -> the derivatives of energy wrt to k
//      if choice=53: enlout(3)               -> the twist derivatives of energy wrt to k, needed for magnetization
//      if choice=24: enlout(1:9*natom)       -> the forces and the frozen wf part of dyn. mat.
//      if choice=6 : enlout(1:6*(3*natom+6)) -> the frozen wf part of elastic tensor
// --If (paw_opt==3)
//    if choice=1 : enlout(nnlout)= contribution to <c|S|c>  (nnlout=1)
// --If (paw_opt==4)
//    not available
// ==== if (signs==2) ====
// --if (paw_opt=0, 1 or 4)
//    vectout(2,npwout*nspinor)=result of the aplication of the concerned operator
//                or one of its derivatives to the input vect.:
//      if (choice=1) <G|V_nonlocal|vect_start>
//      if (choice=2) <G|dV_nonlocal/d(atm coord)|vect_start>
//      if (choice=3) <G|dV_nonlocal/d(strain)|vect_start>
//      if (choice=5) <G|dV_nonlocal/dk|vect_start>
//      if (choice=51) <G|d(right)V_nonlocal/dk|vect_start>
//      if (choice=52) <G|d(left)V_nonlocal/dk|vect_start>
//      if (choice=53) <G|d(twist)V_nonlocal/dk|vect_start>
//  if (paw_opt=2)
//    vectout(2,npwout*nspinor)=final vector in reciprocal space:
//      if (choice=1) <G|V_nonlocal-lamdba.(I+S)|vect_start>
//      if (choice=2) <G|d[V_nonlocal-lamdba.(I+S)]/d(atm coord)|vect_start>
//      if (choice=3) <G|d[V_nonlocal-lamdba.(I+S)]/d(strain)|vect_start>
//      if (choice=5) <G|d[V_nonlocal-lamdba.(I+S)]/dk|vect_start>
//      if (choice=51) <G|d(right)[V_nonlocal-lamdba.(I+S)]/dk|vect_start>
//      if (choice=52) <G|d(left)[V_nonlocal-lamdba.(I+S)]/dk|vect_start>
//      if (choice=53) <G|d(twist)[V_nonlocal-lamdba.(I+S)]/dk|vect_start>
// --if (paw_opt=3 or 4)
//    svectout(2,npwout*nspinor)=result of the aplication of Sij (overlap matrix)
//                  or one of its derivatives to the input vect.:
//      if (choice=1) <G|I+S|vect_start>
//      if (choice=2) <G|dS/d(atm coord)|vect_start>
//      if (choice=3) <G|dS/d(strain)|vect_start>
//      if (choice=5) <G|dS/dk|vect_start>
//      if (choice=51) <G|d(right)S/dk|vect_start>
//      if (choice=52) <G|d(left)S/dk|vect_start>
//      if (choice=53) <G|d(twist)S/dk|vect_start>
//
// SIDE EFFECTS
//  ==== ONLY IF useylm=1
//  cprjin(natom,nspinor) <type(cprj_type)>=projected input wave function |in> on non-local projectors
//                                  =<p_lmn|in> and derivatives
//                    Treatment depends on cpopt parameter:
//                     if cpopt=-1, <p_lmn|in> (and derivatives)
//                                  are computed here (and not saved)
//                     if cpopt= 0, <p_lmn|in> are computed here and saved
//                                  derivatives are eventually computed but not saved
//                     if cpopt= 1, <p_lmn|in> and first derivatives are computed here and saved
//                                  other derivatives are eventually computed but not saved
//                     if cpopt= 2  <p_lmn|in> are already in memory;
//                                  first (and 2nd) derivatives are computed here and not saved
//                     if cpopt= 3  <p_lmn|in> are already in memory;
//                                  first derivatives are computed here and saved
//                                  other derivatives are eventually computed but not saved
//                     if cpopt= 4  <p_lmn|in> and first derivatives are already in memory;
//                                  other derivatives are not computed
//                                  This option is not compatible with choice=4,24 or 6
// (if useylm=0, should have cpopt=-1)
//
// NOTES
// * See nonlop_pl and nonlop_ylm to have more comments...
// * In the case signs=1, the array vectout is not used, nor modified
//   so that the same array as vectin can be used as a dummy argument;
//   the same is true for the pairs npwin-npwout, ffnlin-ffnlout,
//   kgin-kgout, ph3din-ph3dout, phkredin-phkxredout).
//
// PARENTS
//      bloch_interp,cgwf,dyfnl3,eig1fixed,energy,forstrnps,getgh1c
//      getghc,getgsc,ladielmt,lavnl,lobpcgIIIwf,lobpcgIIwf,lobpcgccIIIwf
//      lobpcgccIIwf,mag_loc_k,nstwf4,prctfvw1,prctfvw2,prep_nonlop,resp3dte
//      vso_realspace_nonlop,vtowfk
//
// CHILDREN
//      leave_new,nonlop_pl,nonlop_ylm,timab,wrtout
//
// SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "abi_gpu_header.h"


//GPU memory areas
static double2 *vectin_gpu,*ph3din_gpu,*proj_gpu,*dproj_gpu;
static double2 *vectout_gpu,*svectout_gpu,*ph3dout_gpu;
static double2 *val_ajlmn_gpu,*val_sajlmn_gpu;
static double *ffnlin_gpu,*kpgin_gpu,*ffnlout_gpu,*enl_gpu,*sij_gpu,*rdlmn_gpu,*rdproj_gpu,*enlout_gpu,*d_enlk_gpu,*gprimd_gpu;
static unsigned char *typat_gpu,*lmn_gpu,*nlmn_gpu;
static unsigned short int *atoms_gpu;
static int *indlmn_gpu,*atindx1_gpu,*kgin_gpu;

//CPU Transfert buffers
static unsigned char *typat,*lmn,*nlmn;
static unsigned short int *atoms;

//
static int gpu_initialization=0;
static int nb_proj_to_compute=0;
static int m_ham_used=0;
static int ffnl_ph3d_updated=0;


//Compute Routine
extern "C" void gpu_nonlop_(int *atindx1,int *choice,int *cpopt,double *proj,int *dimenl1,int *dimenl2,int *dimffnlin,int *dimffnlout,
			    double *enl,double *enlout,double *ffnlin,double *ffnlout,double *gprimd,int *idir,int *indlmn,int *istwf_k,
			    int *kgin,int *kgout,double *kpgin,double *kpgout,double *kptin,double *kptout,double *lambda,int *lmnmax,int *matblk,int *mgfft,
			    int *mpi_enreg_me_g0,int *natom,int *nattyp,int *ngfft,int *nkpgin,int *nkpgout,int *nloalg,
			    int *nnlout,int *npwin,int *npwout,int *nspinor,int *ntypat,int *paw_opt,double *phkxredin,
			    double *phkxredout,double *ph1d,double *ph3din,double *ph3dout,int *signs,double *sij,double *svectout,
			    double *pi,double *ucvol,double *vectin,double *vectout)
{

  //  use defs_basis
  //  use defs_datatypes
  //  use defs_abitypes

  // !This section has been created automatically by the script Abilint (TD).
  // !Do not modify the following lines by hand.
  //  use interfaces_14_hidewrite
  //  use interfaces_16_hideleave
  //  use interfaces_18_timing
  //  use interfaces_65_nonlocal, except_this_one => nonlop
  // !End of the abilint section

  //  implicit none

  // !Arguments ------------------------------------
  // !scalars
  //  integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimffnlin,dimffnlout,idir
  //  integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,nkpgin
  //  integer,intent(in) :: nkpgout,nnlout,npwin,npwout,nspinor,ntypat,only_SO
  //  integer,intent(in) :: paw_opt,signs,tim_nonlop,useylm
  //  real(dp),intent(in) :: lambda,ucvol
  //  type(MPI_type),intent(inout) :: mpi_enreg
  // !arrays
  //  integer,intent(in) :: atindx1(natom),indlmn(6,lmnmax,ntypat),kgin(3,npwin)
  //  integer,intent(in) :: kgout(3,npwout),nattyp(ntypat),ngfft(18),nloalg(5)
  // !integer,intent(in) :: pspso(ntypat) ! Unused
  //  real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinor**2)
  //  real(dp),intent(in) :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
  //  real(dp),intent(in) :: ffnlout(npwout,dimffnlout,lmnmax,ntypat),gmet(3,3)
  //  real(dp),intent(in) :: gprimd(3,3),kpgin(npwin,nkpgin*useylm)
  //  real(dp),intent(in) :: kpgout(npwout,nkpgout*useylm),kptin(3),kptout(3)
  //  real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),phkxredin(2,natom)
  //  real(dp),intent(in) :: phkxredout(2,natom),sij(dimenl1,ntypat*((paw_opt+1)/3))
  //  real(dp),intent(inout) :: ph3din(2,npwin,matblk),ph3dout(2,npwout,matblk)
  //  real(dp),intent(inout) :: vectin(2,npwin*nspinor)
  //  real(dp),intent(out) :: enlout(nnlout),svectout(2,npwout*nspinor*(paw_opt/3))
  //  real(dp),intent(out) :: vectout(2,npwout*nspinor)
  //  type(cprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5))

  // Local variables-------------------------------
  // scalars
  //  character(len=500) :: message
  //  integer :: i
  //arrays
  //  real(dp) :: tsec(2)

  cudaError_t cuda_state;
  double2 vectin_0={0.,0.};


  //We check for functionality
  if(((*choice)!=0)&&((*choice)!=1)&&((*choice)!=2)&&((*choice)!=3)&&((*choice)!=23)){
    printf("gpu_nonlop: Error:\n choice %d was used but only choice 0,1,2,3 and 23 are available on GPU\n",(*choice));
    printf(" Try to change your input files with correct optforces & optstress ...\n");
    fflush(stdout);
    abi_cabort();
  }

  if((*nspinor)!=1){
    printf("gpu_nonlop: Error:\n Only nspinor = 1 is allowed at the moment but you used %d\n",(*nspinor));
    fflush(stdout);
    abi_cabort();
  }

  if(((*choice)==1)&&((*signs)!=2)){
    printf("gpu_nonlop: Error:\n when choice is 1 only signs=2 is implemented. ( %d was used) \n",(*signs));
    fflush(stdout);
    abi_cabort();
  }
  if((((*choice)==2)||((*choice)==3)||((*choice)==23))&&((*signs)!=1)){
    printf("gpu_nonlop: Error:\n when choice is 2,3 or 23 only signs=1 is implemented. ( %d was used) \n",(*signs));
    fflush(stdout);
    abi_cabort();
  }

  if((*nloalg)<=0){
    printf("gpu_nonlop: Error:\n Only nloalg(1)>0 is allowed at the moment but you used %d\n",(*nloalg));
    fflush(stdout);
    abi_cabort();
  }

  if( ((*paw_opt)>0) && ((2*(*dimenl1)/((*lmnmax)*((*lmnmax)+1)))!=1) ){
    printf("gpu_nonlop: Error:\n only the real case for enl has been implemented on gpu and complex was used\n");
    fflush(stdout);
    abi_cabort();
  }

  //GPU allocation
  if(!gpu_initialization)
    alloc_nonlop_gpu_(npwin,npwout,nspinor,
		      natom,ntypat,lmnmax,
		      indlmn,nattyp,
		      atindx1,gprimd,
		      dimffnlin,
		      dimenl1,dimenl2,dimffnlout);

  //If no projections everything is zero
  if (nb_proj_to_compute==0){
    int i;
    if((*signs)==1){
      for (i=0;i<(*nnlout);i++) {enlout[i]=0.;}
    }
    if((*signs)==2){
      if (((*paw_opt)==0)||((*paw_opt)==1)||((*paw_opt)==4)){
        for (i=0;i<2*(*npwout);i++) {vectout[i]=0.;}
      }
      if (((*paw_opt)==2)&&((*choice)==1)){
        for (i=0;i<2*(*npwout);i++) {vectout[i]=-(*lambda)*vectin[i];}
      }
      if (((*paw_opt)==2)&&((*choice)>1)){
        for (i=0;i<2*(*npwout);i++) {vectout[i]=0.;}
      }
      if (((*paw_opt)==3)||((*paw_opt)==4)){
        if ((*choice)==1){
          for (i=0;i<2*(*npwout);i++) {svectout[i]=vectin[i];}
        }
        if ((*choice)>1){
          for (i=0;i<2*(*npwout);i++) {svectout[i]=0.;}
        }
      }
    }
    return;
  }

  //Trick to avoid complex branchement in cuda kernels
  if( ((*istwf_k)==2) && ((*mpi_enreg_me_g0)==1) ){
    vectin_0.x  = vectin[0];
    vectin_0.y  = vectin[1];
    vectin[0] = vectin_0.x/2;
    vectin[1] = 0.;
  }

  //Memcopies
  cuda_state=cudaMemcpy(vectin_gpu,vectin,(*npwin)*sizeof(double2),cudaMemcpyHostToDevice);

  if(ffnl_ph3d_updated!=1){
    cuda_state=cudaMemcpy(ph3din_gpu,ph3din,(*natom)*(*npwin)*sizeof(double2),cudaMemcpyHostToDevice);
    cuda_state=cudaMemcpy(ffnlin_gpu,ffnlin,(*npwin)*(*dimffnlin)*(*lmnmax)*(*ntypat)*sizeof(double),cudaMemcpyHostToDevice);
  }

  if((*choice)>=2){
    if((*nkpgin)>0)
      cuda_state=cudaMemcpy(kpgin_gpu,kpgin,(*npwin)*(*nkpgin)*sizeof(double),cudaMemcpyHostToDevice);
    else{
      cuda_state=cudaMemcpy(kgin_gpu,kgin,3*(*npwin)*sizeof(int),cudaMemcpyHostToDevice);
      gpu_mkkpg_(kgin_gpu,kpgin_gpu,kptin,npwin);
    }
  }

  if(cuda_state!=cudaSuccess){
    printf("gpu_nonlop: Error while copying data to gpu : %s \n",cudaGetErrorString(cuda_state));
    fflush(stdout);
    abi_cabort();
  }

  //Compute Projection
  char cplex;
  if((*istwf_k)>1)
    cplex=1;
  else
    cplex=2;

  gpu_compute_nl_projections_(proj_gpu,dproj_gpu,
			      vectin_gpu,ph3din_gpu,
			      ffnlin_gpu,kpgin_gpu,
			      indlmn_gpu,atoms_gpu,lmn_gpu,typat_gpu,
			      &nb_proj_to_compute,npwin,choice,
			      dimffnlin,lmnmax,&cplex,pi,ucvol);

  //Copy back projections if wanted
  if(((*cpopt)==0)||((*cpopt)==1))
    cudaMemcpy(proj , proj_gpu,nb_proj_to_compute*sizeof(double2),cudaMemcpyDeviceToHost);

  if((*choice)>0){
    const double four_pi_by_ucvol=4*(*pi)/sqrt(*ucvol);

    //Memcopies
    //No need since ffnlout=ffnlin and ph3dout=ph3din
    //     if((*choice)==1){
    //       cuda_state=cudaMemcpy(ph3dout_gpu,ph3dout,(*natom)*(*npwout)*sizeof(double2),cudaMemcpyHostToDevice);
    //       cuda_state=cudaMemcpy(ffnlout_gpu,ffnlout,(*npwout)*(*dimffnlout)*(*lmnmax)*(*ntypat)*sizeof(double),cudaMemcpyHostToDevice);
    //     }

    if(m_ham_used != 1) {
      if((*paw_opt)!=3)
	cuda_state=cudaMemcpy(enl_gpu, enl, (*dimenl1)*(*dimenl2)*(*nspinor)*(*nspinor)*sizeof(double),cudaMemcpyHostToDevice);
      if((*paw_opt)>1)
      cuda_state=cudaMemcpy(sij_gpu, sij, (*dimenl1)*(*ntypat)*sizeof(double),cudaMemcpyHostToDevice);

      if((*choice==3)||((*choice==23)))
	cuda_state=cudaMemcpy(gprimd_gpu,gprimd,9*sizeof(double),cudaMemcpyHostToDevice);

      if(cuda_state!=cudaSuccess){
	printf("gpu_nonlop: Error while copying data 2 to gpu :\n %s \n",cudaGetErrorString(cuda_state));
	fflush(stdout);
	abi_cabort();
      }
    }

    gpu_compute_nl_hamiltonian_(proj_gpu,dproj_gpu,
				vectin_gpu,
				vectout_gpu,svectout_gpu,
				val_ajlmn_gpu,val_sajlmn_gpu,
				ph3din_gpu,ffnlin_gpu,// <== ph3dout_gpu,ffnlout_gpu,
				rdlmn_gpu,rdproj_gpu,
				enlout_gpu,d_enlk_gpu,
				enl_gpu,sij_gpu,gprimd_gpu,
				atindx1_gpu,nlmn_gpu,
				indlmn_gpu,atoms_gpu,
				lmn_gpu,typat_gpu,
				paw_opt,dimffnlout,dimenl1,
				npwout,&nb_proj_to_compute,lmnmax,
				natom,choice,signs,
				&four_pi_by_ucvol,lambda);

    //We copy back results
    if((*choice)==1){
      if((*paw_opt)!=3)
	cuda_state=cudaMemcpy(vectout,vectout_gpu,(*npwout)*sizeof(double2),cudaMemcpyDeviceToHost);
      if((*paw_opt)>2)
	cuda_state=cudaMemcpy(svectout,svectout_gpu,(*npwout)*sizeof(double2),cudaMemcpyDeviceToHost);
    }
    else if(((*choice)==2) && ((*signs)==1))
      cuda_state=cudaMemcpy(enlout,enlout_gpu,3*(*natom)*sizeof(double),cudaMemcpyDeviceToHost);
    else if(((*choice)==3) && ((*signs)==1))
      cuda_state=cudaMemcpy(enlout,enlout_gpu,6*sizeof(double),cudaMemcpyDeviceToHost);
    else if(((*choice)==23) && ((*signs)==1))
      cuda_state=cudaMemcpy(enlout,enlout_gpu,(6+3*(*natom))*sizeof(double),cudaMemcpyDeviceToHost);

    if(cuda_state!=cudaSuccess){
      printf("gpu_nonlop: Error while retrieving results from gpu :\n %s \n",cudaGetErrorString(cuda_state));
      fflush(stdout);
      abi_cabort();
    }
  }

  //We put back the correct value of vectin and correct svectout if needed
  if( ((*istwf_k)==2) && ((*mpi_enreg_me_g0)==1) ){
    vectin[0] = vectin_0.x;
    vectin[1] = vectin_0.y;
    if((*paw_opt)>2){
      svectout[0] += vectin[0]/2;
      svectout[1] += vectin[1];
    }
  }

  //Impossible
  if(!gpu_initialization)
    free_nonlop_gpu_();

}// end subroutine gpu_nonlop



//Allocation routine
extern "C" void alloc_nonlop_gpu_(int *npwin,int *npwout,int *nspinor,
				  int *natom,int *ntypat,int *lmnmax,
				  int *indlmn,int *nattyp,
				  int *atindx1,double *gprimd,
				  int *dimffnlin,int *dimffnlout,
				  int *dimenl1, int *dimenl2 ){

//   printf("calling alloc_nonlop_gpu with npw=%d \n",*npwin);
//   fflush(stdout);

  //Si on avait deja alloue
  if(gpu_initialization==1)
    free_nonlop_gpu_();
  else
    gpu_initialization = 1;

  nb_proj_to_compute=0;

  cudaError_t cuda_state;

  cudaMallocHost(&nlmn,(*ntypat)*sizeof(char));
  //Calcul du nombre de projections
  for(int itypat=0;itypat<(*ntypat);itypat++){
    int count_ilmn=0;
    for(int ilmn=0;ilmn<(*lmnmax);ilmn++){
      if(indlmn[2+6*(ilmn + (*lmnmax)*(itypat))] > 0)
	count_ilmn++;
    }
    nlmn[itypat]=count_ilmn;
    nb_proj_to_compute+=count_ilmn*nattyp[itypat];
  }

  cuda_state=cudaMalloc(&vectin_gpu, (*npwin)*sizeof(double2));
  cuda_state=cudaMalloc(&vectout_gpu, (*npwout)*sizeof(double2));
  cuda_state=cudaMalloc(&svectout_gpu, (*npwout)*sizeof(double2));
  cuda_state=cudaMalloc(&ph3din_gpu, (*natom)*(*npwin)*sizeof(double2));
  cuda_state=cudaMalloc(&ffnlin_gpu, (*npwin)*(*dimffnlin)*(*lmnmax)*(*ntypat)*sizeof(double));
  cuda_state=cudaMalloc(&kpgin_gpu, (*npwin)*3*sizeof(double));
  cuda_state=cudaMalloc(&kgin_gpu, (*npwin)*3*sizeof(int));
  cuda_state=cudaMalloc(&ph3dout_gpu, (*npwout)*(*natom)*sizeof(double2));
  cuda_state=cudaMalloc(&ffnlout_gpu, (*npwout)*(*dimffnlout)*(*lmnmax)*(*ntypat)*sizeof(double));
  cuda_state=cudaMalloc(&enl_gpu,(*dimenl1)*(*dimenl2)*(*nspinor)*(*nspinor)*sizeof(double));
  cuda_state=cudaMalloc(&sij_gpu,(*dimenl1)*(*ntypat)*sizeof(double));

  cuda_state=cudaMalloc(&indlmn_gpu, 6*(*lmnmax)*(*ntypat)*sizeof(int));
  cuda_state=cudaMalloc(&atindx1_gpu, (*natom)*sizeof(int));
  cuda_state=cudaMalloc(&typat_gpu, nb_proj_to_compute*sizeof(char));
  cuda_state=cudaMalloc(&nlmn_gpu,(*ntypat)*sizeof(char));
  cuda_state=cudaMalloc(&lmn_gpu, nb_proj_to_compute*sizeof(char));
  cuda_state=cudaMalloc(&atoms_gpu, nb_proj_to_compute*sizeof(short int));

  cuda_state=cudaMalloc(&proj_gpu, nb_proj_to_compute*sizeof(double2));
  cuda_state=cudaMalloc(&dproj_gpu, 10*nb_proj_to_compute*sizeof(double2));
  cuda_state=cudaMalloc(&rdproj_gpu, 7*nb_proj_to_compute*sizeof(double));
  cuda_state=cudaMalloc(&rdlmn_gpu, 3*(*natom)*(*lmnmax)*sizeof(double));

  cuda_state=cudaMalloc(&d_enlk_gpu, 7*sizeof(double));
  cuda_state=cudaMalloc(&gprimd_gpu, 9*sizeof(double));
  cuda_state=cudaMalloc(&enlout_gpu, (6+3*(*natom))*sizeof(double));
  cuda_state=cudaMalloc(&val_ajlmn_gpu, nb_proj_to_compute*sizeof(double2));
  cuda_state=cudaMalloc(&val_sajlmn_gpu, nb_proj_to_compute*sizeof(double2));


  cuda_state=cudaMallocHost(&typat, nb_proj_to_compute*sizeof(char));
  cuda_state=cudaMallocHost(&lmn, nb_proj_to_compute*sizeof(char));
  cuda_state=cudaMallocHost(&atoms, nb_proj_to_compute*sizeof(short int));

  if(cuda_state!=cudaSuccess){
    printf("alloc_nonlop_gpu: Error during gpu memory allocation :\n %s \n",cudaGetErrorString(cuda_state));
    fflush(stdout);
    abi_cabort();
  }

  //Precompute couple (atom,lmn) for each projection
  int iproj=0;
  int iatom=0;
  for(int itypat=0;itypat<(*ntypat);itypat++){
    for(int iat=0;iat<nattyp[itypat];iat++){
      for(int ilmn=0;ilmn<nlmn[itypat];ilmn++){
	typat[iproj]=itypat;
	lmn[iproj]=ilmn;
	atoms[iproj]=iatom;
	iproj++;
      }
      iatom++;
    }
  }
  assert(iproj==nb_proj_to_compute);

  cuda_state=cudaMemcpy(indlmn_gpu,indlmn, 6*(*lmnmax)*(*ntypat)*sizeof(int),cudaMemcpyHostToDevice);
  cuda_state=cudaMemcpy(typat_gpu,typat, nb_proj_to_compute*sizeof(char),cudaMemcpyHostToDevice);
  cuda_state=cudaMemcpy(lmn_gpu,lmn, nb_proj_to_compute*sizeof(char),cudaMemcpyHostToDevice);
  cuda_state=cudaMemcpy(atoms_gpu,atoms,nb_proj_to_compute*sizeof(short int),cudaMemcpyHostToDevice);
  cuda_state=cudaMemcpy(atindx1_gpu, atindx1, (*natom)*sizeof(int),cudaMemcpyHostToDevice);
  cuda_state=cudaMemcpy(nlmn_gpu, nlmn, (*ntypat)*sizeof(char),cudaMemcpyHostToDevice);
  cuda_state=cudaMemcpy(gprimd_gpu,gprimd,9*sizeof(double),cudaMemcpyHostToDevice);
  if(cuda_state!=cudaSuccess){
    printf("alloc_nonlop_gpu: Error while copying data to gpu :\n %s \n",cudaGetErrorString(cuda_state));
    fflush(stdout);
    abi_cabort();
  }
  cuda_state=cudaMemset(rdlmn_gpu,0,3*(*natom)*(*lmnmax)*sizeof(double));
  if(cuda_state!=cudaSuccess){
    printf("alloc_nonlop_gpu: Error while set tabs to 0:\n %s \n",cudaGetErrorString(cuda_state));
    fflush(stdout);
    abi_cabort();
  }
}


//Deallocation routine
extern "C" void free_nonlop_gpu_(){

  gpu_initialization=0;

  cudaError_t cuda_state;

  //Free Memory
  cuda_state=cudaFreeHost(nlmn);

  cuda_state=cudaFree(vectin_gpu);
  cuda_state=cudaFree(vectout_gpu);
  cuda_state=cudaFree(svectout_gpu);
  cuda_state=cudaFree(ph3din_gpu);
  cuda_state=cudaFree(ffnlin_gpu);
  cuda_state=cudaFree(kpgin_gpu);
  cuda_state=cudaFree(kgin_gpu);
  cuda_state=cudaFree(ph3dout_gpu);
  cuda_state=cudaFree(ffnlout_gpu);
  cuda_state=cudaFree(enl_gpu);
  cuda_state=cudaFree(sij_gpu);

  cuda_state=cudaFree(indlmn_gpu);
  cuda_state=cudaFree(atindx1_gpu);
  cuda_state=cudaFree(typat_gpu);
  cuda_state=cudaFree(nlmn_gpu);
  cuda_state=cudaFree(lmn_gpu);
  cuda_state=cudaFree(atoms_gpu);

  cuda_state=cudaFree(proj_gpu);
  cuda_state=cudaFree(dproj_gpu);
  cuda_state=cudaFree(rdproj_gpu);
  cuda_state=cudaFree(rdlmn_gpu);

  cuda_state=cudaFree(d_enlk_gpu);
  cuda_state=cudaFree(gprimd_gpu);
  cuda_state=cudaFree(enlout_gpu);
  cuda_state=cudaFree(val_ajlmn_gpu);
  cuda_state=cudaFree(val_sajlmn_gpu);

  cuda_state=cudaFreeHost(typat);
  cuda_state=cudaFreeHost(lmn);
  cuda_state=cudaFreeHost(atoms);

  if(cuda_state!=cudaSuccess){
    printf("free_nonlop_gpu: Error while freeing gpu data:\n %s \n",cudaGetErrorString(cuda_state));
    fflush(stdout);
    abi_cabort();
  }
}

extern "C" void gpu_update_ham_data_(double *enl,int *size_enl, double *sij,int *size_sij,
                                     double *gprimd,int *size_gprimd){

  cudaError_t cuda_state;

  cuda_state=cudaMemcpy(enl_gpu, enl, (*size_enl)*sizeof(double),cudaMemcpyHostToDevice);
  if(cuda_state!=cudaSuccess){
    printf("gpu_update_ham_data: Error while copying data 1 to gpu :\n %s \n",cudaGetErrorString(cuda_state));
    fflush(stdout);
    abi_cabort();
  }
  if((*size_sij)>0){
    cuda_state=cudaMemcpy(sij_gpu, sij, (*size_sij)*sizeof(double),cudaMemcpyHostToDevice);
    if(cuda_state!=cudaSuccess){
      printf("gpu_update_ham_data: Error while copying data 2 to gpu :\n %s \n",cudaGetErrorString(cuda_state));
      fflush(stdout);
      abi_cabort();
    }
  }

  cuda_state=cudaMemcpy(gprimd_gpu,gprimd,(*size_gprimd)*sizeof(double),cudaMemcpyHostToDevice);
  if(cuda_state!=cudaSuccess){
    printf("gpu_update_ham_data: Error while copying data 3 to gpu :\n %s \n",cudaGetErrorString(cuda_state));
    fflush(stdout);
    abi_cabort();
  }
  m_ham_used = 1;
}

extern "C" void gpu_update_ffnl_ph3d_(double *ph3din,int *dimph3din,double *ffnlin,int *dimffnlin){

  cudaError_t cuda_state;

  cuda_state=cudaMemcpy(ffnlin_gpu,ffnlin,(*dimffnlin)*sizeof(double),cudaMemcpyHostToDevice);
  if(cuda_state!=cudaSuccess){
    printf("gpu_update_ffnl_ph3d: Error while copying data 1 to gpu :\n %s \n",cudaGetErrorString(cuda_state));
    fflush(stdout);
    abi_cabort();
  }
  cuda_state=cudaMemcpy(ph3din_gpu,ph3din,(*dimph3din)*sizeof(double),cudaMemcpyHostToDevice);
  if(cuda_state!=cudaSuccess){
    printf("gpu_update_ffnl_ph3d: Error while copying data 2 to gpu :\n %s \n",cudaGetErrorString(cuda_state));
    fflush(stdout);
    abi_cabort();
  }

  ffnl_ph3d_updated = 1;
}

extern "C" void gpu_finalize_ffnl_ph3d_(){
  ffnl_ph3d_updated = 0;
}


extern "C" void gpu_finalize_ham_data_(){
  m_ham_used = 0;
}
//***
