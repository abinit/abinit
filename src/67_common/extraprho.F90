!{\src2tex{textfont=tt}}
!!****f* ABINIT/extraprho
!!
!! NAME
!! extraprho
!!
!! FUNCTION
!! Extrapolate electronic density for new ionic positions
!! from values of density of previous SCF cycle.
!! Use algorithm proposed by D. Alfe in Comp. Phys. Comm. 118 (1999), 31-33
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cg(2,mcg)= plane wave wavefunction coefficient
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | densty(ntypat,4)=parameters for initialisation of the gaussian density
!!   | jellslab,slabzbeg,slabzend,slabwsrad=parameters for jellium slab
!!   | natom=number of atoms in cell.
!!   | nspden=number of spin-density components
!!  gmet(3,3)=reciprocal space metric
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=cutoff value on G**2 for sphere inside fft box
!!  istep=number of call the routine
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mgfft=maximum size of 1D FFTs
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=informations about MPI parallelization
!!  mqgrid=number of grid pts in q array for f(q) spline.
!!  my_natom=number of atoms treated by current processor
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT
!!  npwarr=(nkpt)=number of planewaves in basis at this k point
!!  ntypat=number of types of atoms in cell
!!  pawtab(ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  psps<type(pseudopotential_type)>=variables related to pseudopotentials
!!  qgrid(mqgrid)=q grid for spline from 0 to qmax
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  ucvol=unit cell volume (bohr**3).
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  xred_new(3,natom)=new reduced coordinates for atoms in unit cell
!!  xred_old(3,natom)=old reduced coordinates for atoms in unit cell
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  zion(ntypat)=charge on each type of atom
!!  znucl(ntypat)=atomic number of each atom type
!!
!! SIDE EFFECTS
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= PAW rhoij occupancies and related data
!!                                            Value from previous SCF cycle is input
!!                                            Extrapolated value is output
!!  rhor(nfft,nspden)=the density from previous SCF cycle is input
!!                    the extrapolated density is output
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      atm2fft,extrapwf,jellium,pawrhoij_alloc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine extraprho(atindx,atindx1,cg,dtset,gmet,gprimd,gsqcut,istep,&
& kg,mcg,mgfft,mpi_enreg,mqgrid,my_natom,nattyp,nfft,ngfft,npwarr,ntypat,pawrhoij,&
& pawtab,ph1d,psps,qgrid,rhor,rprimd,scf_history,ucvol,usepaw,&
& xred_new,xred_old,ylm,zion,znucl)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_scf_history
 use m_errors

 use m_atomdata, only : atom_length
 use m_pawtab,   only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_alloc

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'extraprho'
 use interfaces_14_hidewrite
 use interfaces_64_psp
 use interfaces_67_common, except_this_one => extraprho
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,mcg,mgfft,my_natom,mqgrid,nfft,ntypat,usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(scf_history_type),intent(inout) :: scf_history
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem)
 integer,intent(in) :: nattyp(ntypat),ngfft(18),npwarr(dtset%nkpt)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(in) :: qgrid(mqgrid),rprimd(3,3)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) ::  zion(ntypat),znucl(ntypat)
 real(dp), intent(inout) :: cg(2,mcg)
 real(dp),intent(inout) :: rhor(nfft,dtset%nspden),xred_new(3,dtset%natom)
 real(dp),intent(in) :: xred_old(3,dtset%natom)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: cplex,dplex,iatom,ii,ind1,ind1new,ind2,ind2new,irhoij,ispden,itypat,jrhoij,klmn
 integer :: lmn2_size,nselect,nspden_rhoij,optatm,optdyfr,opteltfr,optgr,option,optn,optn2
 integer :: optstr,optv
 real(dp) :: a11,a12,a22,a33,alpha,b1,b2,beta,detA,fact,ratio1,ratio2
 logical :: hasmoved,usegauss
 character(len=500) :: message
!arrays
 integer :: dummy3(3)
 real(dp) :: diff_t(3),diff_tmdt(3),diff_tpdt(3),dummy2(2)
 real(dp) :: dummy_in(0)
 real(dp) :: dummy_out1(0),dummy_out2(0),dummy_out3(0),dummy_out4(0),dummy_out5(0),dummy_out6(0)
 real(dp) :: strn_dummy6(6),strv_dummy6(6)
 real(dp),allocatable :: deltarho(:),gauss(:,:),rhoijtmp(:,:),work1(:)
 real(dp),allocatable :: work2(:,:),work3(:,:),xred_tpdt(:,:)

! *************************************************************************

!---------------------------------------------------------------
!----------- Inits
!---------------------------------------------------------------

!History indexes
 ind1=scf_history%hindex(1)
 ind2=scf_history%hindex(2)

!Compatibility tests
 if (ind1==0.and.ind2>0)then
   MSG_BUG(' Incompatible history indexes !')
 end if

!Rotated values of history indexes
 if (ind1>0.and.ind2>0) then
   ind1new=ind2;ind2new=ind1
 else if (ind1>0.and.ind2==0) then
   ind1new=3-ind1;ind2new=ind1
 else if (ind1==0.and.ind2==0) then
   ind1new=1;ind2new=0
 end if

!Compute ionic positions at t+dt in red. coordinates
!Has to take the boundary conditions into account
 ABI_ALLOCATE(xred_tpdt,(3,dtset%natom))
 do iatom=1,dtset%natom
   xred_tpdt(1,iatom)=xred_old(1,iatom)+mod(xred_new(1,iatom)-xred_old(1,iatom)+1.5_dp,one)-half
   xred_tpdt(2,iatom)=xred_old(2,iatom)+mod(xred_new(2,iatom)-xred_old(2,iatom)+1.5_dp,one)-half
   xred_tpdt(3,iatom)=xred_old(3,iatom)+mod(xred_new(3,iatom)-xred_old(3,iatom)+1.5_dp,one)-half
 end do

!---------------------------------------------------------------
!----------- Compute Alpha and Beta
!----------- see (4) in Comp. Phys. Comm. 118 (1999), 31-33
!---------------------------------------------------------------

!Compute a_ij matrix
 a11=zero;a12=zero;a22=zero;a33=zero;b1=zero;b2=zero
 diff_t=zero;diff_tmdt=zero;diff_tpdt=zero
 do iatom=1,dtset%natom

   diff_tpdt(1:3)=xred_tpdt(1:3,iatom)-xred_old(1:3,iatom)
   if (ind1>0) then
     diff_t(1:3)=scf_history%xreddiff(1:3,iatom,ind1)
     if (ind2>0) diff_tmdt(1:3)=scf_history%xreddiff(1:3,iatom,ind2)
   end if
   do ii=1,3
     a11=a11+diff_t(ii)**2
     a22=a22+diff_tmdt(ii)**2
     a33=a33+diff_tpdt(ii)**2
     a12=a12+diff_t(ii)   *diff_tmdt(ii)
     b1 =b1 +diff_t(ii)   *diff_tpdt(ii)
     b2 =b2 +diff_tmdt(ii)*diff_tpdt(ii)
   end do

!  Store reduced coordinates diffs in SCF history
   scf_history%xreddiff(1:3,iatom,ind1new)=diff_tpdt(1:3)

 end do
 ABI_DEALLOCATE(xred_tpdt)
 hasmoved=(a11>=tol10.or.a22>=tol10.or.a33>=tol10)

!Compute alpha and beta
 alpha=zero;beta=zero
 if (hasmoved.and.ind1>0) then
   ratio1=one;if (abs(a33)>=tol10) ratio1=(a11+a33-two*b1)/a33
   ratio2=one;if (abs(a33)>=tol10) ratio2=(a11+a33-two*b2)/a33
   detA=a11*a22-a12**2
   if (abs(a11)>=tol10.and.(abs(a22)<tol10.or.abs(detA)<tol10)) then
     alpha=b1/a11
   else if (abs(a22)>=tol10.and.(abs(a11)<tol10.or.abs(detA)<tol10)) then
     beta=b2/a22
   else if (abs(ratio1)+abs(ratio2)<tol6) then
     if (ind2>0) then
       alpha=two;beta=-one
     else
       alpha=one
     end if
     write(message,'(3a,f4.1,a,f4.1)')&
&     'Ionic positions lead to a collinear system !',ch10,&
&     'Mixing coeffs have been set to: alpha=',alpha,' beta=',beta
     MSG_WARNING(message)
   else if (abs(a11)>=tol10.and.abs(a22)>=tol10) then
     alpha=(b1*a22-b2*a12)/detA
     beta =(b2*a11-b1*a12)/detA
   end if
 end if

!---------------------------------------------------------------
!----------- Contribution from delta_rho(t), delta_rho(t-dt)
!----------- and delta_rho(t-2dt) to predicted rho(t+dt)
!---------------------------------------------------------------

!deltarho(t+dt) <- deltarho(t) + alpha.[deltarho(t)-deltarho(t-dt)]
!+ beta .[deltarho(t-dt)-deltarho(t-2dt)]
!Note: scf_history%deltarhor is updated at the same time

 ABI_ALLOCATE(deltarho,(nfft))
 do ispden=1,dtset%nspden

   if (ispden==1) then
     deltarho(:)=rhor(:,ispden)-scf_history%atmrho_last(:)
   else if (ispden==2.and.dtset%nspden==2) then
     deltarho(:)=rhor(:,ispden)-half*scf_history%atmrho_last(:)
   end if


!  rho(t+dt) <- deltarho(t) + alpha.deltarho(t)
   if (dtset%nspden/=4.or.ispden==1) then
     rhor(:,ispden)=(one+alpha)*deltarho(:)
   else
     rhor(:,ispden)=(one+alpha)*rhor(:,ispden)
   end if

   if (hasmoved) then

!    rho(t+dt) <- -alpha.deltarho(t-dt) + beta.deltarho(t-dt)
     if (abs(beta-alpha)>tol14.and.ind1>0) then
       rhor(:,ispden)=rhor(:,ispden)+(beta-alpha)*scf_history%deltarhor(:,ispden,ind1)
     end if

!    rho(t+dt) <- -beta.deltarho(t-2dt)
     if (abs(beta)>tol14.and.ind2>0) then
       rhor(:,ispden)=rhor(:,ispden)-beta*scf_history%deltarhor(:,ispden,ind2)
     end if

   end if

!  Store deltarho(t) in history
   if (dtset%nspden/=4.or.ispden==1) then
     scf_history%deltarhor(:,ispden,ind1new)=deltarho(:)
   else
     scf_history%deltarhor(:,ispden,ind1new)=rhor(:,ispden)
   end if

 end do

 ABI_DEALLOCATE(deltarho)

!---------------------------------------------------------------
!----------- Contribution from rho_at(t+dt) to predicted rho(t+dt)
!---------------------------------------------------------------

!Determine whether a gaussian atomic density has to be used or not
!MG: Note that there's a small inconsistency between initro and extraprho because in initrho
! we use `use_gaussian(ntypat)`.
 usegauss=.true.
 if (usepaw==0) usegauss = any(.not.psps%nctab(1:ntypat)%has_tvale)
 if (usepaw==1) usegauss=(minval(pawtab(1:ntypat)%has_tvale)==0)
 if (usegauss) then
   optn2=3
   ABI_ALLOCATE(gauss,(2,ntypat))
   do itypat=1,ntypat
     gauss(1,itypat)=zion(itypat)
     gauss(2,itypat) = atom_length(dtset%densty(itypat,1),zion(itypat),znucl(itypat))
   end do
   call wrtout(std_out," Extrapolating rho(t+dt) using gaussian functions as atomic densities", "COLL")
 else
   optn2=2
   ABI_ALLOCATE(gauss,(2,0))
   call wrtout(std_out," Extrapolating rho(t+dt) using atomic densities taken from pseudos", "COLL")
 end if

!Compute rho_at(t+dt) as sum of atomic densities
!Note: scf_history%atmrho_last is updated at the same time
 optatm=1;optdyfr=0;opteltfr=0;optgr=0;optstr=0;optv=0;optn=1
 call atm2fft(atindx1,scf_history%atmrho_last,dummy_out1,dummy_out2,dummy_out3,&
& dummy_out4,gauss,gmet,gprimd,dummy_out5,dummy_out6,gsqcut,mgfft,mqgrid,dtset%natom,nattyp,&
& nfft,ngfft,ntypat,optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,psps,pawtab,ph1d,qgrid,&
& dummy3,dummy_in,strn_dummy6,strv_dummy6,ucvol,usepaw,dummy_in,dummy_in,dummy_in,dummy2,dummy_in,&
& comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
& paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)
 ABI_DEALLOCATE(gauss)

!Take eventually into account jellium slab
 if (dtset%jellslab/=0) then
   option=2
   ABI_ALLOCATE(work1,(nfft))
   ABI_ALLOCATE(work2,(nfft,1))
   ABI_ALLOCATE(work3,(2,nfft))
   work2(:,1)=scf_history%atmrho_last(:)
   call jellium(gmet,gsqcut,mpi_enreg,nfft,ngfft,1,option,dtset%paral_kgb,&
&   dtset%slabwsrad,work3,work2,rprimd,work1,dtset%slabzbeg,dtset%slabzend)
   scf_history%atmrho_last(:)=work2(:,1)
   ABI_DEALLOCATE(work1)
   ABI_DEALLOCATE(work2)
   ABI_DEALLOCATE(work3)
 end if

!Add rho_at(t+dt) to rho(t+dt)
 rhor(:,1)=rhor(:,1)+scf_history%atmrho_last(:)
 if (dtset%nspden==2) rhor(:,2)=rhor(:,2)+half*scf_history%atmrho_last(:)

!---------------------------------------------------------------
!----------- Extrapolation of PAW rhoij occupancy matrixes
!---------------------------------------------------------------

 if (usepaw==1) then

   if (ind2==0) then
     nspden_rhoij=dtset%nspden;if (dtset%pawspnorb>0.and.dtset%nspinor==2) nspden_rhoij=4
     call pawrhoij_alloc(scf_history%pawrhoij(:,ind1new),dtset%pawcpxocc,nspden_rhoij,&
&     dtset%nspinor,dtset%nsppol,dtset%typat,pawtab=pawtab,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   end if

   do iatom=1,my_natom

     nspden_rhoij=pawrhoij(iatom)%nspden
     lmn2_size=pawrhoij(iatom)%lmn2_size
     cplex=pawrhoij(iatom)%cplex;dplex=cplex-1

     if (hasmoved) then
       ABI_ALLOCATE(rhoijtmp,(cplex*lmn2_size,nspden_rhoij))
       rhoijtmp=zero

       do ispden=1,nspden_rhoij

!        rhoij(t+dt) <- rhoij(t) + alpha.rhoij(t)
         fact=one+alpha
         jrhoij=1
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=cplex*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
           rhoijtmp(klmn:klmn+dplex,ispden)=rhoijtmp(klmn:klmn+dplex,ispden) &
&           +fact*pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
           jrhoij=jrhoij+cplex
         end do

!        rhoij(t+dt) <- -alpha.rhoij(t-dt) + beta.rhoij(t-dt)
         if (abs(beta-alpha)>tol14.and.ind1>0) then
           fact=beta-alpha
           jrhoij=1
           do irhoij=1,scf_history%pawrhoij(iatom,ind1)%nrhoijsel
             klmn=cplex*scf_history%pawrhoij(iatom,ind1)%rhoijselect(irhoij)-dplex
             rhoijtmp(klmn:klmn+dplex,ispden)=rhoijtmp(klmn:klmn+dplex,ispden) &
&             +fact*scf_history%pawrhoij(iatom,ind1)%rhoijp(jrhoij:jrhoij+dplex,ispden)
             jrhoij=jrhoij+cplex
           end do
         end if

!        rho(t+dt) <- -beta.rhoij(t-2dt)
         if (abs(beta)>tol14.and.ind2>0) then
           fact=-beta
           jrhoij=1
           do irhoij=1,scf_history%pawrhoij(iatom,ind2)%nrhoijsel
             klmn=cplex*scf_history%pawrhoij(iatom,ind2)%rhoijselect(irhoij)-dplex
             rhoijtmp(klmn:klmn+dplex,ispden)=rhoijtmp(klmn:klmn+dplex,ispden) &
&             +fact*scf_history%pawrhoij(iatom,ind2)%rhoijp(jrhoij:jrhoij+dplex,ispden)
             jrhoij=jrhoij+cplex
           end do
         end if

       end do !ispden
     end if !hasmoved

!    Store rhoij(t) in history
!    (cannot use pawrhoij_copy here because update for single atom)
     nselect=pawrhoij(iatom)%nrhoijsel
     scf_history%pawrhoij(iatom,ind1new)%nrhoijsel=nselect
     scf_history%pawrhoij(iatom,ind1new)%rhoijselect(1:nselect)=pawrhoij(iatom)%rhoijselect(1:nselect)
     scf_history%pawrhoij(iatom,ind1new)%rhoijp(1:cplex*nselect,1:nspden_rhoij)= &
&     pawrhoij(iatom)%rhoijp(1:cplex*nselect,1:nspden_rhoij)

!    Select non-zero values of rhoij(t+dt)
     if (hasmoved) then
       nselect=0
       if (cplex==1) then
         do klmn=1,lmn2_size
           if (any(abs(rhoijtmp(klmn,:))>tol10)) then
             nselect=nselect+1
             pawrhoij(iatom)%rhoijselect(nselect)=klmn
             do ispden=1,nspden_rhoij
               pawrhoij(iatom)%rhoijp(nselect,ispden)=rhoijtmp(klmn,ispden)
             end do
           end if
         end do
       else
         do klmn=1,lmn2_size
           if (any(abs(rhoijtmp(2*klmn-1:2*klmn,:))>tol10)) then
             nselect=nselect+1
             pawrhoij(iatom)%rhoijselect(nselect)=klmn
             do ispden=1,nspden_rhoij
               pawrhoij(iatom)%rhoijp(2*nselect-1:2*nselect,ispden)=rhoijtmp(2*klmn-1:2*klmn,ispden)
             end do
           end if
         end do
       end if
       pawrhoij(iatom)%nrhoijsel=nselect
       ABI_DEALLOCATE(rhoijtmp)
     end if

   end do !iatom
 end if !usepaw

 scf_history%alpha=alpha
 scf_history%beta=beta

!---------------------------------------------------------------
!----------- End
!---------------------------------------------------------------

 if(scf_history%usecg==1) then
   if (hasmoved) then
     call extrapwf(atindx,atindx1,cg,dtset,istep,kg,mcg,mgfft,mpi_enreg,nattyp,&
&     ngfft,npwarr,ntypat,pawtab,psps,rprimd,scf_history,usepaw,xred_old,ylm)
   else
     scf_history%cg(:,:,2)=zero
   end if

 end if

!Rotate history indexes
 scf_history%hindex(1)=ind1new
 scf_history%hindex(2)=ind2new

end subroutine extraprho
!!***
