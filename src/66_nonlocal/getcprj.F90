!{\src2tex{textfont=tt}}
!!****f* ABINIT/getcprj
!! NAME
!! getcprj
!!
!! FUNCTION
!!  Compute <Proj_i|Cnk> for one wave function |Cnk> expressed in reciprocal space.
!!  Compute also derivatives of <Proj_i|Cnk>.
!!  |Proj_i> are non-local projectors (for each atom and each l,m,n)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  choice=chooses possible output:
!!    In addition to projected wave function:
!!    choice=1 => nothing else
!!          =2 => 1st gradients with respect to atomic position(s)
!!          =3 => 1st gradients with respect to strain(s)
!!          =23=> 1st gradients with respect to strain(s) and atm pos
!!          =4 => 2nd derivatives with respect to atomic pos.
!!          =24=> 1st and 2nd derivatives with respect to atomic pos.
!!          =5 => 1st gradients with respect to k wavevector
!!          =6 => 2nd derivatives with respect to strain and atm. pos.
!!  cpopt=1 if <Proj_i|Cnk> are already in memory; see below (side effects).
!!  cwavef(2,nspinor*npw_k)=input cmplx wavefunction coefficients <G|Cnk>
!!  ffnl(npw_k,dimffnl,lmnmax,ntypat)=nonlocal form factors to be used for the application of the nl operator
!!  idir=direction of the derivative, i.e. dir. of - atom to be moved  in the case choice=2
!!                                                 - strain component  in the case choice=3
!!                                                 - k point direction in the case choice=5
!!       Compatible only with choice=2,3,5; if idir=0, all derivatives are computed
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=lmn
!!  istwf_k=option parameter that describes the storage of wfs
!!  kg_k(3,npw_k)=reduced planewave coordinates
!!  kpg(npw_k,npk)=(k+G) components and related data
!!  kpoint(3)=k point in terms of recip. translations
!!  lmnmax=max. number of (l,m,n) components over all types of atoms
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nattyp(ntypat)=number of atoms of each type
!!  ngfft(18)=contain all needed information about 3D FFT, see ~ABINIT/Infos/vargs.htm#ngfft
!!  nloalg(3)=governs the choice of the algorithm for nonlocal operator
!!  npw_k=number of planewaves for given k point
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  ntypat=number of types of atoms in unit cell
!!  phkxred(2,natom)=phase factors exp(2 pi kpoint.xred)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  ph3d(2,npw_k,natom)=3D structure factors, for each atom and plane wave
!!                      only used if nloalg(2)>0
!!  ucvol= unit cell volume
!!  useylm=governs the way the nonlocal operator is to be applied
!!
!! SIDE EFFECTS
!!  cwaveprj(natom,nspinor) <type(pawcprj_type)>=projected input wave function <Proj_i|Cnk> with all NL projectors
!!                                (and derivatives)
!!                                if cpopt=1 the projected scalars have been already been computed and
!!                                           only derivatives are computed here
!!                                if cpopt=0 the projected scalars and derivatives are computed here
!!
!! TODO
!!  Spin-orbit
!!
!! PARENTS
!!      cgwf,ctocprj,debug_tools,dfpt_accrho,dfpt_nstpaw,ks_ddiago,m_wfd
!!      rf2_init
!!
!! CHILDREN
!!      mkkpg,opernla_ylm,ph1d3d
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine getcprj(choice,cpopt,cwavef,cwaveprj,ffnl,&
&                   idir,indlmn,istwf_k,kg_k,kpg,kpoint,lmnmax,mgfft,mpi_enreg,&
&                   natom,nattyp,ngfft,nloalg,npw_k,nspinor,ntypat,&
&                   phkxred,ph1d,ph3d,ucvol,useylm)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_pawcprj, only : pawcprj_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getcprj'
 use interfaces_56_recipspace
 use interfaces_66_nonlocal, except_this_one => getcprj
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: choice,cpopt,idir,istwf_k,lmnmax
 integer,intent(in) :: mgfft,natom,npw_k,nspinor,ntypat,useylm
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),kg_k(3,npw_k),nattyp(ntypat)
 integer,intent(in) :: ngfft(18),nloalg(3)
 real(dp),intent(in) :: cwavef(2,npw_k*nspinor)
 real(dp),intent(in),target :: ffnl(:,:,:,:),kpg(:,:),ph3d(:,:,:)
 real(dp),intent(in) :: kpoint(3),ph1d(2,3*(2*mgfft+1)*natom),phkxred(2,natom)
 type(pawcprj_type),intent(inout) :: cwaveprj(natom,nspinor)

!Local variables-------------------------------
!scalars
 integer :: choice_,cplex,dimffnl,ia,ia1,ia2,ia3,ia4,iatm,ic,ii,ilmn,ispinor,itypat
 integer :: jc,matblk,mincat,nd2gxdt,ndgxdt,nincat,nkpg,nkpg_,nlmn,signs
!arrays
 integer,allocatable :: cplex_dgxdt(:),cplex_d2gxdt(:),indlmn_typ(:,:)
 real(dp),allocatable :: d2gxdt(:,:,:,:,:),dgxdt(:,:,:,:,:),ffnl_typ(:,:,:)
 real(dp),allocatable :: gx(:,:,:,:)
 real(dp), ABI_CONTIGUOUS pointer :: kpg_(:,:),ph3d_(:,:,:)

! *********************************************************************

 DBG_ENTER('COLL')

!Nothing to do in that case
 if (cpopt==1.and.choice==1) return

!Not available for useylm=0
 if (useylm==0) then
   MSG_ERROR('Not available for useylm=0 !')
 end if

!Error on bad choice
 if ((choice<1.or.choice>6).and.choice/=23.and.choice/=24) then
   MSG_BUG('Does not presently support this choice !')
 end if

!Error on bad idir
 if (idir>0.and.choice/=2.and.choice/=3.and.choice/=5) then
   MSG_BUG('Does not support idir>0 for this choice')
 end if

!Error on sizes
 nkpg=size(kpg,2)
 if (nkpg>0) then
   if( (choice==2.and.nkpg<3) .or. &
&   ((choice==4.or.choice==24).and.nkpg<9) .or. &
&   ((choice==6.or.choice==3.or.choice==23).and.nkpg<3) ) then
     MSG_BUG('Incorrect size for kpg array !')
   end if
 end if
 if (size(ffnl,1)/=npw_k.or.size(ffnl,3)/=lmnmax) then
   MSG_BUG('Incorrect size for ffnl!')
 end if
 if (size(ph3d)>0) then
   if (size(ph3d,2)/=npw_k) then
     MSG_BUG('Incorrect size for ph3d!')
   end if
 end if

!Define dimensions of projected scalars
 dimffnl=size(ffnl,2)
 ndgxdt=0;nd2gxdt=0
 if (idir==0) then
   if (choice==2) ndgxdt=3
   if (choice==3) ndgxdt=6
   if (choice==23) ndgxdt=9
   if (choice==4) nd2gxdt=6
   if (choice==24) then
     ndgxdt=3;nd2gxdt=6
   end if
   if (choice==5) ndgxdt=3
   if (choice==6) then
     ndgxdt=9;nd2gxdt=54
   end if
 else
   ndgxdt=1
 end if

!Eventually re-compute (k+G) vectors (and related data)
 if (nkpg==0) then
   nkpg_=0
   if (choice==4.or.choice==24) nkpg_=9
   if (choice==2.or.choice==3.or.choice==23) nkpg_=3
   ABI_ALLOCATE(kpg_,(npw_k,nkpg_))
   if (nkpg_>0) then
     call mkkpg(kg_k,kpg_,kpoint,nkpg_,npw_k)
   end if
 else
   nkpg_=nkpg
   kpg_ => kpg
 end if

!Some other dims
 mincat=min(NLO_MINCAT,maxval(nattyp))
 cplex=2;if (istwf_k>1) cplex=1
 choice_=choice;if (cpopt==1) choice_=-choice
 signs=1;if (idir>0) signs=2

!Eventually allocate temporary array for ph3d
 if (nloalg(2)<=0) then
   matblk=mincat
   ABI_ALLOCATE(ph3d_,(2,npw_k,matblk))
 else
   matblk=size(ph3d,3)
   ph3d_ => ph3d
 end if

!Loop over atom types
 ia1=1;iatm=0
 do itypat=1,ntypat
   ia2=ia1+nattyp(itypat)-1;if (ia2<ia1) cycle
   nlmn=count(indlmn(3,:,itypat)>0)

!  Retrieve some data for this type of atom
   ABI_ALLOCATE(indlmn_typ,(6,nlmn))
   ABI_ALLOCATE(ffnl_typ,(npw_k,dimffnl,nlmn))
   indlmn_typ(:,1:nlmn)=indlmn(:,1:nlmn,itypat)
   ffnl_typ(:,:,1:nlmn)=ffnl(:,:,1:nlmn,itypat)

!  Loop on blocks of atoms inside type
   do ia3=ia1,ia2,mincat
     ia4=min(ia2,ia3+mincat-1);nincat=ia4-ia3+1
!     Prepare the phase factors if they were not already computed
     if (nloalg(2)<=0) then
       call ph1d3d(ia3,ia4,kg_k,matblk,natom,npw_k,ngfft(1),ngfft(2),ngfft(3),&
&       phkxred,ph1d,ph3d_)
     end if

!    Allocate memory for projected scalars
     ABI_ALLOCATE(gx,(cplex,nlmn,nincat,nspinor))
     ABI_ALLOCATE(dgxdt,(cplex,ndgxdt,nlmn,nincat,nspinor))
     ABI_ALLOCATE(d2gxdt,(cplex,nd2gxdt,nlmn,nincat,nspinor))
     ABI_ALLOCATE(cplex_dgxdt,(ndgxdt))
     ABI_ALLOCATE(cplex_d2gxdt,(nd2gxdt))

!    Retrieve eventually <p_i|c> coeffs
     if (cpopt==1) then
       do ispinor=1,nspinor
         do ia=1,nincat
           gx(1:cplex,1:nlmn,ia,ispinor)=cwaveprj(iatm+ia,ispinor)%cp(1:cplex,1:nlmn)
         end do
       end do
     end if

!    Compute <p_i|c> scalars (and derivatives) for this block of atoms
     call opernla_ylm(choice_,cplex,cplex_dgxdt,cplex_d2gxdt,dimffnl,d2gxdt,dgxdt,ffnl_typ,gx,&
&     ia3,idir,indlmn_typ,istwf_k,kpg_,matblk,mpi_enreg,nd2gxdt,ndgxdt,nincat,nkpg_,nlmn,&
&     nloalg,npw_k,nspinor,ph3d_,signs,ucvol,cwavef)

!    Transfer result to output variable cwaveprj
     if (cpopt==0) then
       do ispinor=1,nspinor
         do ia=1,nincat
           cwaveprj(iatm+ia,ispinor)%nlmn=nlmn
           cwaveprj(iatm+ia,ispinor)%cp(1:cplex,1:nlmn)=gx(1:cplex,1:nlmn,ia,ispinor)
           if(cplex==1) cwaveprj(iatm+ia,ispinor)%cp(2,1:nlmn)=zero
         end do
       end do
     end if
     if (cpopt>=0.and.choice>1) then
       if(cplex==2)then
         do ispinor=1,nspinor
           do ia=1,nincat
             cwaveprj(iatm+ia,ispinor)%ncpgr=ndgxdt+nd2gxdt
             if (ndgxdt>0) cwaveprj(iatm+ia,ispinor)%dcp(1:2,1:ndgxdt,1:nlmn)=&
&             dgxdt(1:2,1:ndgxdt,1:nlmn,ia,ispinor)
             if (nd2gxdt>0)cwaveprj(iatm+ia,ispinor)%dcp(1:2,ndgxdt+1:ndgxdt+nd2gxdt,1:nlmn)=&
&             d2gxdt(1:2,1:nd2gxdt,1:nlmn,ia,ispinor)
           end do
         end do
       else
!        cplex_dgxdt(i)  = 1 if dgxdt(1,i,:,:)  is real, 2 if it is pure imaginary
!        cplex_d2gxdt(i) = 1 if d2gxdt(1,i,:,:) is real, 2 if it is pure imaginary
         do ispinor=1,nspinor
           do ia=1,nincat
             cwaveprj(iatm+ia,ispinor)%ncpgr=ndgxdt+nd2gxdt
             if (ndgxdt>0) then
               do ilmn =1,nlmn
                 do ii = 1,ndgxdt
                   ic = cplex_dgxdt(ii) ; jc = 3 - ic
                   cwaveprj(iatm+ia,ispinor)%dcp(ic,ii,ilmn)=dgxdt(1,ii,ilmn,ia,ispinor)
                   cwaveprj(iatm+ia,ispinor)%dcp(jc,ii,ilmn)=zero
                 end do
               end do
             end if
             if (nd2gxdt>0) then
               do ilmn =1,nlmn
                 do ii = 1,nd2gxdt
                   ic = cplex_d2gxdt(ii) ; jc = 3 - ic
                   cwaveprj(iatm+ia,ispinor)%dcp(ic,ndgxdt+ii,ilmn)=d2gxdt(1,ii,ilmn,ia,ispinor)
                   cwaveprj(iatm+ia,ispinor)%dcp(jc,ndgxdt+ii,ilmn)=zero
                 end do
               end do
             end if
           end do
         end do
       end if
     end if

!    End loop inside block of atoms
     iatm=iatm+nincat
     ABI_DEALLOCATE(gx)
     ABI_DEALLOCATE(dgxdt)
     ABI_DEALLOCATE(d2gxdt)
     ABI_DEALLOCATE(cplex_dgxdt)
     ABI_DEALLOCATE(cplex_d2gxdt)
   end do

!  End loop over atom types
   ia1=ia2+1
   ABI_DEALLOCATE(indlmn_typ)
   ABI_DEALLOCATE(ffnl_typ)
 end do

 if (nkpg==0) then
   ABI_DEALLOCATE(kpg_)
 end if
 if (nloalg(2)<=0) then
   ABI_DEALLOCATE(ph3d_)
 end if

 DBG_EXIT('COLL')

 end subroutine getcprj
!!***
