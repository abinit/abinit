!{\src2tex{textfont=tt}}
!!****f* ABINIT/nonlop_gpu
!! NAME
!! nonlop_gpu
!!
!! FUNCTION
!!  Compute application of a nonlocal operator, using GPU (NVidia Cuda)
!!  This routine is an interface to Cuda Kernel gpu_nonlop.cu
!!
!! COPYRIGHT
!! Copyright (C) 2011-2017 ABINIT group (FDahm, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  choice: chooses possible output:
!!    choice=0 => do nothing (only compute WF projected with NL projectors)
!!          =1 => a non-local energy contribution
!!          =2 => a gradient with respect to atomic position(s)
!!          =3 => a gradient with respect to strain(s)
!!          =23=> a gradient with respect to atm. pos. and strain(s)
!!  cpopt=flag defining the status of cprjin%cp(:)=<Proj_i|Cnk> scalars (see below, side effects)
!!  dimenl1,dimenl2=dimensions of enl (see enl)
!!  dimffnlin=second dimension of ffnlin (1+number of derivatives)
!!  dimffnlout=second dimension of ffnlout (1+number of derivatives)
!!  enl(dimenl1,dimenl2,nspinortot**2)=
!!  ->Norm conserving : ==== when paw_opt=0 ====
!!                      (Real) Kleinman-Bylander energies (hartree)
!!                      dimenl1=lmnmax  -  dimenl2=ntypat
!!  ->PAW :             ==== when paw_opt=1, 2 or 4 ====
!!                      (Real or complex, hermitian) Dij coefs to connect projectors
!!                      dimenl1=cplex_enl*lmnmax*(lmnmax+1)/2  -  dimenl2=natom
!!  ffnlin(npwin,dimffnlin,lmnmax,ntypat)=nonlocal form factors to be used
!!          for the application of the nonlocal operator to the |in> vector
!!  ffnlout(npwout,dimffnlout,lmnmax,ntypat)=nonlocal form factors to be used
!!          for the application of the nonlocal operator to the |out> vector
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2),
!!                        - k point direction in the case (choice=5,signs=2)
!!                          for choice 53, twisted derivative involves idir+1 and idir-1
!!                        - strain component (1:6) in the case (choice=3,signs=2) or (choice=6,signs=1)
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=lmn
!!  istwf_k=option parameter that describes the storage of wfs
!!  kgin(3,npwin)=integer coords of planewaves in basis sphere, for the |in> vector
!!  kgout(3,npwout)=integer coords of planewaves in basis sphere, for the |out> vector
!!  kpgin(npw,npkgin)= (k+G) components and related data, for the |in> vector
!!  kpgout(npw,nkpgout)=(k+G) components and related data, for the |out> vector
!!  kptin(3)=k point in terms of recip. translations, for the |in> vector
!!  kptout(3)=k point in terms of recip. translations, for the |out> vector
!!  lambda=factor to be used when computing (Vln-lambda.S) - only for paw_opt=2
!!         Typically lambda is the eigenvalue (or its guess)
!!  lmnmax=max. number of (l,m,n) components over all types of atoms
!!  matblk=dimension of the arrays ph3din and ph3dout
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell
!!  nattyp(ntypat)=number of atoms of each type
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpgin,nkpgout=second sizes of arrays kpgin/kpgout
!!  nloalg(3)=governs the choice of the algorithm for nonlocal operator
!!  nnlout=dimension of enlout (when signs=1 and choice>0):
!!         ==== if paw_opt=0, 1 or 2 ====
!!         choice=1=>nnlout=1   choice=2=>nnlout=3*natom    choice=3=>nnlout=6
!!         ==== if paw_opt=3 ====
!!         choice=1 =>nnlout=1
!!         ==== if paw_opt=4 ====
!!         not available
!!  npwin=number of planewaves for given k point, for the |in> vector
!!  npwout=number of planewaves for given k point, for the |out> vector
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nspinortot=number of spinorial components of the wavefunctions on current proc
!!  ntypat=number of types of atoms in cell
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!  phkxredin(2,natom)=phase factors exp(2 pi kptin.xred)
!!  phkxredout(2,natom)=phase factors exp(2 pi kptout.xred)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1D structure factors phase information
!!  ph3din(2,npwin,matblk)=3D structure factors, for each atom and plane wave (in)
!!  ph3dout(2,npwout,matblk)=3-dim structure factors, for each atom and plane wave (out)
!!  signs= if 1, get contracted elements (energy, forces, stress, ...)
!!         if 2, applies the non-local operator to a function in reciprocal space
!!  sij(dimenl1,ntypat*(paw_opt/3))=overlap matrix components (only if paw_opt=2, 3 or 4)
!!  ucvol=unit cell volume (bohr^3)
!!  vectin(2,npwin*nspinor)=input cmplx wavefunction coefficients <G|Cnk>
!!  [cprjin_left(natom,nspinor)]=The projected input wave function <p_nlm|in_left>
!!    for the left wavefunction. Data are assumed to be in memory, they are NOT recalculated here.
!!    Only signs==1 and choice==1 are supported.
!!
!! OUTPUT
!! ==== if (signs==1) ====
!! --If (paw_opt==0, 1 or 2)
!!    enlout(nnlout)= contribution to the non-local part of the following properties:
!!      if choice=1 : enlout(1)               -> the energy
!!      if choice=2 : enlout(1:3*natom)       -> the forces
!!      if choice=3 : enlout(1:6)             -> the stresses
!!      if choice=23: enlout(1:6+3*natom)     -> the forces and the stresses
!! --If (paw_opt==3)
!!    if choice=1 : enlout(nnlout)= contribution to <c|S|c>  (nnlout=1)
!! --If (paw_opt==4)
!!    not available
!! ==== if (signs==2) ====
!! --if (paw_opt=0, 1 or 4)
!!    vectout(2,npwout*nspinor)=result of the aplication of the concerned operator
!!                or one of its derivatives to the input vect.:
!!      if (choice=1) <G|V_nonlocal|vect_start>
!!      if (choice=2) <G|dV_nonlocal/d(atm coord)|vect_start>
!!      if (choice=3) <G|dV_nonlocal/d(strain)|vect_start>
!!  if (paw_opt=2)
!!    vectout(2,npwout*nspinor)=final vector in reciprocal space:
!!      if (choice=1) <G|V_nonlocal-lamdba.(I+S)|vect_start>
!!      if (choice=2) <G|d[V_nonlocal-lamdba.(I+S)]/d(atm coord)|vect_start>
!!      if (choice=3) <G|d[V_nonlocal-lamdba.(I+S)]/d(strain)|vect_start>
!! --if (paw_opt=3 or 4)
!!    svectout(2,npwout*nspinor)=result of the aplication of Sij (overlap matrix)
!!                  or one of its derivatives to the input vect.:
!!      if (choice=1) <G|I+S|vect_start>
!!      if (choice=2) <G|dS/d(atm coord)|vect_start>
!!      if (choice=3) <G|dS/d(strain)|vect_start>
!!
!! SIDE EFFECTS
!!  cprjin(natom,nspinor) <type(pawcprj_type)>=projected input wave function |in> on non-local projectors
!!                                  =<p_lmn|in> and derivatives
!!                    Treatment depends on cpopt parameter:
!!                     if cpopt=-1, <p_lmn|in> (and derivatives)
!!                                  are computed here (and not saved)
!!                     if cpopt= 0, <p_lmn|in> are computed here and saved
!!                                  derivatives are eventually computed but not saved
!!                     if cpopt= 1, <p_lmn|in> and first derivatives are computed here and saved
!!                                  other derivatives are eventually computed but not saved
!!
!! TODO
!! * Implementation for spinorial wave functions (nspinor=2)
!! * Implementation for response function (phonons, ddk, elastic tensor, ...)
!!
!! PARENTS
!!      nonlop
!!
!! CHILDREN
!!      dotprod_g,gpu_nonlop
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine nonlop_gpu(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimffnlin,dimffnlout,&
&                      enl,enlout,ffnlin,ffnlout,gprimd,idir,indlmn,istwf_k,&
&                      kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                      mpi_enreg,natom,nattyp,ngfft,nkpgin,nkpgout,nloalg,nnlout,&
&                      npwin,npwout,nspinor,nspinortot,ntypat,paw_opt,phkxredin,phkxredout,ph1d,&
&                      ph3din,ph3dout,signs,sij,svectout,ucvol,vectin,vectout)


 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_cgtools

 use m_pawcprj, only : pawcprj_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nonlop_gpu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimffnlin,dimffnlout,idir
 integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,natom,nkpgin,nkpgout,nnlout
 integer,intent(in) :: npwin,npwout,nspinor,nspinortot,ntypat,paw_opt,signs
 real(dp),intent(in) :: lambda,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),indlmn(6,lmnmax,ntypat),kgin(3,npwin)
 integer,intent(in) :: kgout(3,npwout),nattyp(ntypat),ngfft(18),nloalg(3)
 real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2)
 real(dp),intent(in) :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
 real(dp),intent(in) :: ffnlout(npwout,dimffnlout,lmnmax,ntypat) !,gmet(3,3)
 real(dp),intent(in) :: gprimd(3,3),kpgin(npwin,nkpgin),kpgout(npwout,nkpgout)
 real(dp),intent(in) :: kptin(3),kptout(3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: phkxredin(2,natom),phkxredout(2,natom)
 real(dp),intent(in) :: sij(dimenl1,ntypat*((paw_opt+1)/3))
 real(dp),intent(inout) :: ph3din(2,npwin,matblk),ph3dout(2,npwout,matblk)
 real(dp),intent(inout) :: vectin(2,npwin*nspinor)
 real(dp),intent(out) :: enlout(nnlout)
 real(dp),intent(out),target :: svectout(2,npwout*nspinor*(paw_opt/3))
 real(dp),intent(out),target :: vectout (2,npwout*nspinor)
 type(pawcprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5))

!Local variables-------------------------------
!scalars
 integer :: ia,iatom,ilmn,iproj,ispinor,itypat,signs_
 real(dp) :: doti
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: proj(:,:)
 real(dp),pointer :: svectout_(:,:),vectout_(:,:)

! **********************************************************************

 DBG_ENTER("COLL")

!Error on bad choice
 if ((choice<0 .or. choice>3).and. choice/=23 .and. choice/=24) then
   write(msg,'(a,i0,a)')'Does not presently support this choice=',choice,'.'
   MSG_BUG(msg)
 end if
 if (cpopt<-1.or.cpopt>1) then
   msg='  Bad value for cpopt !'
   MSG_BUG(msg)
 end if
 if (nspinor==2) then
   msg='  nspinor=2 (spinorial WF) not yet allowed !'
   MSG_ERROR(msg)
 end if

 if ((cpopt==0).or.(cpopt==1))  then
   ABI_ALLOCATE(proj,(2,lmnmax*natom))
   proj=zero;
 end if

!Workaround to get choice=1/signs=1 working
 if (choice==1.and.signs==1) then
   signs_=2
   ABI_ALLOCATE(vectout_,(2,npwin*nspinor))
   ABI_ALLOCATE(svectout_,(2,npwin*nspinor*(paw_opt/3)))
 else
   signs_=signs;vectout_=>vectout;svectout_=>svectout
 end if

#if defined HAVE_GPU_CUDA
 call gpu_nonlop(atindx1,choice,cpopt,proj,dimenl1,dimenl2,dimffnlin,dimffnlout,&
& enl,enlout,ffnlin,ffnlout,gprimd,idir,indlmn,istwf_k,&
& kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
& mpi_enreg%me_g0,natom,nattyp,ngfft,nkpgin,nkpgout,nloalg,nnlout,&
& npwin,npwout,nspinor,ntypat,paw_opt,phkxredin,phkxredout,ph1d,&
& ph3din,ph3dout,signs_,sij,svectout_,pi,ucvol,vectin,vectout_)
#endif

 if (choice==1.and.signs==1) then
   if (paw_opt/=3) then
     call dotprod_g(enlout(1),doti,istwf_k,npwin*nspinor,1,vectin,vectout_,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   else
     call dotprod_g(enlout(1),doti,istwf_k,npwin*nspinor,1,vectin,svectout_,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   end if 
   ABI_DEALLOCATE(vectout_)
   ABI_DEALLOCATE(svectout_)
 else
   nullify(vectout_,svectout_)
 end if

 if ((cpopt==0).or.(cpopt==1)) then
   iproj=0
   do ispinor=1,nspinor
     iatom=0
     do itypat=1,ntypat
       do ia=1,nattyp(itypat)
         iatom=iatom+1;cprjin(iatom,1)%nlmn=count(indlmn(3,:,itypat)>0)
         do ilmn=1,cprjin(iatom,1)%nlmn
           iproj=iproj+1;cprjin(iatom,1)%cp(:,ilmn)=proj(:,iproj)
         end do
       end do
     end do
   end do
   ABI_DEALLOCATE(proj)
 end if

 DBG_EXIT("COLL")

!Fake statements to satisfy ABI rules
#if ! defined HAVE_GPU_CUDA
 if (.false.) then
   write(std_out,*) atindx1,enl,ffnlin,ffnlout,gprimd
   write(std_out,*) idir,istwf_k,kgin,kgout,kpgin,kpgout
   write(std_out,*) kptin,kptout,lambda,mpi_enreg%me
   write(std_out,*) ngfft,nloalg,ph1d,ph3din,ph3dout
   write(std_out,*) phkxredin,phkxredout,signs,sij
   write(std_out,*) ucvol,vectin
 end if
#endif

end subroutine nonlop_gpu
!!***
