!!****m* ABINIT/m_extraprho
!! NAME
!!  m_extraprho
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (MT, FJ)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_extraprho

 use defs_basis
 use m_abicore
 use m_scf_history
 use m_errors
 use m_xmpi
 use m_cgtools
 use m_dtset

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes, only : MPI_type
 use m_atomdata, only : atom_length
 use m_numeric_tools,   only : hermit
 use m_geometry, only : metric
 use m_kg, only : getph
 use m_jellium,  only : jellium
 use m_atm2fft,  only : atm2fft
 use m_pawtab,   only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_alloc, pawrhoij_inquire_dim, pawrhoij_filter
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_copy, pawcprj_get, pawcprj_lincom, &
                        pawcprj_free, pawcprj_zaxpby,pawcprj_axpby, pawcprj_put, pawcprj_getdim
 use m_mpinfo,   only : proc_distrb_cycle
 use m_cgprj,    only : ctocprj

 implicit none

 private
!!***

 public :: extraprho
!!***

contains
!!***

!!****f* ABINIT/extraprho
!!
!! NAME
!! extraprho
!!
!! FUNCTION
!! Extrapolate electronic density for new ionic positions
!! from values of density of previous SCF cycle.
!! Use algorithm proposed by D. Alfe in Comp. Phys. Comm. 118 (1999), 31-33 [[cite:Alfe1999]]
!!
!! INPUTS
!!  atindx
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cg(2,mcg)= plane wave wavefunction coefficient
!!  cprj(natom,mcprj*usecprj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each NL proj |p_lmn>
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
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mpi_enreg=information about MPI parallelization
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
!!      m_scfcv_core
!!
!! CHILDREN
!!      cgcprj_cholesky,dotprod_set_cgcprj,lincom_cgcprj,pawcprj_alloc
!!      pawcprj_axpby,pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_put,zgesv
!!
!! SOURCE

subroutine extraprho(atindx,atindx1,cg,cprj,dtset,gmet,gprimd,gsqcut,istep,&
& kg,mcg,mcprj,mgfft,mpi_enreg,mqgrid,my_natom,nattyp,nfft,ngfft,npwarr,ntypat,pawrhoij,&
& pawtab,ph1d,psps,qgrid,rhor,rprimd,scf_history,ucvol,usepaw,&
& xred_new,xred_old,ylm,zion,znucl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,mcg,mcprj,mgfft,my_natom,mqgrid,nfft,ntypat,usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(scf_history_type),intent(inout) :: scf_history
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem)
 integer,intent(in) :: nattyp(ntypat),ngfft(18),npwarr(dtset%nkpt)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) ::  zion(ntypat),znucl(ntypat)
 real(dp), intent(inout) :: cg(2,mcg)
 real(dp),intent(inout) :: rhor(nfft,dtset%nspden),rprimd(3,3),xred_new(3,dtset%natom)
 real(dp),intent(in) :: xred_old(3,dtset%natom)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
 type(pawcprj_type),intent(inout) :: cprj(:,:)

!Local variables-------------------------------
!scalars
 integer :: cplex_rhoij,dplex,iatom,ii,ind1,ind1new,ind2,ind2new,iq,iq0,irhoij,ispden,itypat,jrhoij,klmn
 integer :: lmn2_size,nselect,nspden_rhoij,optatm,optdyfr,opteltfr,optgr,option,optn,optn2
 integer :: optstr,optv,qphase_rhoij
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
!----------- see (4) in Comp. Phys. Comm. 118 (1999), 31-33 [[cite:Alfe1999]]
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
& dummy3,dtset%rcut,dummy_in,rprimd,strn_dummy6,strv_dummy6,ucvol,usepaw,dummy_in,dummy_in,dummy_in,dummy2,dummy_in,&
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
   call jellium(gmet,gsqcut,mpi_enreg,nfft,ngfft,1,option,&
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
     call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,nspden_rhoij=nspden_rhoij,&
&                nspden=dtset%nspden,spnorb=dtset%pawspnorb,cpxocc=dtset%pawcpxocc)
     call pawrhoij_alloc(scf_history%pawrhoij(:,ind1new),cplex_rhoij,nspden_rhoij,&
&     dtset%nspinor,dtset%nsppol,dtset%typat,pawtab=pawtab,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   end if

   do iatom=1,my_natom

     nspden_rhoij=pawrhoij(iatom)%nspden
     lmn2_size=pawrhoij(iatom)%lmn2_size
     cplex_rhoij=pawrhoij(iatom)%cplex_rhoij;dplex=cplex_rhoij-1
     qphase_rhoij=pawrhoij(iatom)%qphase

     if (hasmoved) then
       ABI_ALLOCATE(rhoijtmp,(cplex_rhoij*qphase_rhoij*lmn2_size,nspden_rhoij))
       rhoijtmp=zero

       do ispden=1,nspden_rhoij
         do iq=1,qphase_rhoij
           iq0=merge(0,cplex_rhoij*lmn2_size,iq==1)

!          rhoij(t+dt) <- rhoij(t) + alpha.rhoij(t)
           fact=one+alpha
           jrhoij=1+iq0
           do irhoij=1,pawrhoij(iatom)%nrhoijsel
             klmn=cplex_rhoij*pawrhoij(iatom)%rhoijselect(irhoij)-dplex+iq0
             rhoijtmp(klmn:klmn+dplex,ispden)=rhoijtmp(klmn:klmn+dplex,ispden) &
&             +fact*pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
             jrhoij=jrhoij+cplex_rhoij
           end do

!          rhoij(t+dt) <- -alpha.rhoij(t-dt) + beta.rhoij(t-dt)
           if (abs(beta-alpha)>tol14.and.ind1>0) then
             fact=beta-alpha
             jrhoij=1+iq0
             do irhoij=1,scf_history%pawrhoij(iatom,ind1)%nrhoijsel
               klmn=cplex_rhoij*scf_history%pawrhoij(iatom,ind1)%rhoijselect(irhoij)-dplex+iq0
               rhoijtmp(klmn:klmn+dplex,ispden)=rhoijtmp(klmn:klmn+dplex,ispden) &
&               +fact*scf_history%pawrhoij(iatom,ind1)%rhoijp(jrhoij:jrhoij+dplex,ispden)
               jrhoij=jrhoij+cplex_rhoij
             end do
           end if

!          rho(t+dt) <- -beta.rhoij(t-2dt)
           if (abs(beta)>tol14.and.ind2>0) then
             fact=-beta
             jrhoij=1+iq0
             do irhoij=1,scf_history%pawrhoij(iatom,ind2)%nrhoijsel
               klmn=cplex_rhoij*scf_history%pawrhoij(iatom,ind2)%rhoijselect(irhoij)-dplex+iq0
               rhoijtmp(klmn:klmn+dplex,ispden)=rhoijtmp(klmn:klmn+dplex,ispden) &
&               +fact*scf_history%pawrhoij(iatom,ind2)%rhoijp(jrhoij:jrhoij+dplex,ispden)
               jrhoij=jrhoij+cplex_rhoij
             end do
           end if

         end do ! iq
       end do !ispden
     end if !hasmoved

!    Store rhoij(t) in history
!    (cannot use pawrhoij_copy here because update for single atom)
     nselect=pawrhoij(iatom)%nrhoijsel
     scf_history%pawrhoij(iatom,ind1new)%nrhoijsel=nselect
     scf_history%pawrhoij(iatom,ind1new)%rhoijselect(:)=0
     scf_history%pawrhoij(iatom,ind1new)%rhoijselect(1:nselect)=pawrhoij(iatom)%rhoijselect(1:nselect)
     scf_history%pawrhoij(iatom,ind1new)%rhoijp(1:cplex_rhoij*nselect,1:nspden_rhoij)= &
&     pawrhoij(iatom)%rhoijp(1:cplex_rhoij*nselect,1:nspden_rhoij)

!    Select non-zero values of rhoij(t+dt)
     if (hasmoved) then
       call pawrhoij_filter(pawrhoij(iatom)%rhoijp,pawrhoij(iatom)%rhoijselect,pawrhoij(iatom)%nrhoijsel,&
&                           cplex_rhoij,qphase_rhoij,lmn2_size,nspden_rhoij,&
&                           rhoij_input=rhoijtmp)
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
     if (dtset%extrapwf==1) then
       call extrapwf(atindx,atindx1,cg,dtset,istep,kg,mcg,mgfft,mpi_enreg,nattyp,&
&       ngfft,npwarr,ntypat,pawtab,psps,rprimd,scf_history,usepaw,xred_old,ylm)
     elseif(dtset%extrapwf==2) then
       scf_history%hindex(1)=ind1
       scf_history%hindex(2)=ind2
       scf_history%hindex(3)=ind1new
       call extrapwf_biortho(atindx1,cg,cprj,dtset,istep,mcg,mcprj,mpi_enreg,&
&         nattyp,npwarr,pawtab,scf_history)
     end if
   else
     scf_history%cg(:,:,2)=zero
   end if

 end if
!Rotate history indexes
 scf_history%hindex(1)=ind1new
 scf_history%hindex(2)=ind2new


end subroutine extraprho
!!***

!!****f* ABINIT/extrapwf
!!
!! NAME
!! extrapwf
!!
!! FUNCTION
!! Extrapolate wavefunctions for new ionic positions
!! from values of wavefunctions of previous SCF cycle.
!! Use algorithm proposed by T. A.  Arias et al. in PRB 45, 1538 (1992) [[cite:Arias1992]]
!!
!! INPUTS
!!  atindx(natom)=index table for atoms
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  istep=number of call the routine
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  ngfft(18)=contain all needed information about 3D FFT
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  ntypat=number of types of atoms in cell
!!  pawtab(ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps<type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  xred_old(3,natom)=old reduced coordinates for atoms in unit cell
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!
!! SIDE EFFECTS
!!  cg(2,mcg)= plane wave wavefunction coefficient
!!                          Value from previous SCF cycle is input
!!                          Extrapolated value is output
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!
!! NOTES
!!  THIS ROUTINE IS NOT USEABLE AT PRESENT.
!!  SHOULD BE CAREFULY TESTED AND DEBUGGED (ESPECIALLY WITHIN PAW).
!!
!! PARENTS
!!      m_extraprho
!!
!! CHILDREN
!!      cgcprj_cholesky,dotprod_set_cgcprj,lincom_cgcprj,pawcprj_alloc
!!      pawcprj_axpby,pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_put,zgesv
!!
!! SOURCE

subroutine extrapwf(atindx,atindx1,cg,dtset,istep,kg,mcg,mgfft,mpi_enreg,&
& nattyp,ngfft,npwarr,ntypat,pawtab,psps,rprimd,scf_history,usepaw,xred_old,ylm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,mcg,mgfft,ntypat,usepaw
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(scf_history_type),intent(inout) :: scf_history
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem),nattyp(ntypat),ngfft(18)
 integer,intent(in) :: npwarr(dtset%nkpt)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(inout) :: cg(2,mcg)
 real(dp),intent(in) :: xred_old(3,dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)


!Local variables-------------------------------
!scalars
 integer :: ia,iat,iatom,iband_max,iband_max1,iband_min,iband_min1,ibd,ibg,iblockbd,iblockbd1,icg,icgb,icgb1,icgb2
 integer :: ierr,ig,ii,ikpt,ilmn1,ilmn2,inc,ind1,ind2,iorder_cprj
 integer :: isize,isppol,istep1,istwf_k,itypat,klmn,me_distrb,my_nspinor
 integer :: nband_k,nblockbd,nprocband,npw_k,npw_nk,spaceComm_band
 real(dp) :: dotr,dotr1,doti,doti1,eigval
 !character(len=500) :: message
!arrays
 real(dp) :: alpha(2),beta(2),gmet(3,3),gprimd(3,3),rmet(3,3),ph1d(2,3*(2*mgfft+1)*dtset%natom),ucvol
 integer,allocatable :: bufsize(:),bufsize_wf(:),bufdisp(:),bufdisp_wf(:),dimcprj(:),npw_block(:),npw_disp(:)
 real(dp),allocatable :: al(:,:),anm(:),cwavef(:,:),cwavef1(:,:),cwavef_tmp(:,:),deltawf1(:,:),deltawf2(:,:)
 real(dp),allocatable :: eig(:),evec(:,:)
 real(dp),allocatable :: unm(:,:,:)
 real(dp),allocatable :: work(:,:),work1(:,:),wf1(:,:),ylmgr_k(:,:,:),zhpev1(:,:),zhpev2(:)
 complex(dpc),allocatable :: unm_tmp(:,:),anm_tmp(:,:)
 type(pawcprj_type),allocatable :: cprj(:,:),cprj_k(:,:),cprj_k1(:,:),cprj_k2(:,:),cprj_k3(:,:),cprj_k4(:,:)
!complex(dpc) :: aa

! *************************************************************************

 if (istep==0) return

!Useful array
 if (usepaw==1) then
   ABI_ALLOCATE(dimcprj,(dtset%natom))
   call pawcprj_getdim(dimcprj,dtset%natom,nattyp,ntypat,dtset%typat,pawtab,'O')
 end if

!Metric
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!History indexes
 ind1=1;ind2=2

!First step
 if (istep==1) then
   scf_history%cg(:,:,ind1)=cg(:,:)
!  scf_history%cg(:,:,ind2)=zero
   scf_history%cg(:,:,ind2)= cg(:,:)
   if(usepaw==1) then
!    WARNING: THIS SECTION IS USELESS; NOW crpj CAN BE READ FROM SCFCV
     call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred_old)
     iatom=0 ; iorder_cprj=0
     call pawcprj_alloc(scf_history%cprj(:,:,ind1),0,dimcprj)
     call pawcprj_alloc(scf_history%cprj(:,:,ind2),0,dimcprj)
     ABI_ALLOCATE(ylmgr_k,(dtset%mpw,3,0))
     call ctocprj(atindx,cg,1,scf_history%cprj(:,:,ind1),gmet,gprimd,&
&     iatom,0,iorder_cprj,dtset%istwfk,kg,dtset%kptns,mcg,scf_history%mcprj,&
&     dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,&
&     dtset%natom,nattyp,dtset%nband,dtset%natom,ngfft,dtset%nkpt,&
&     dtset%nloalg,npwarr,dtset%nspinor,dtset%nsppol,dtset%ntypat,&
&     dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,ucvol,0,&
&     xred_old,ylm,ylmgr_k)
     ABI_DEALLOCATE(ylmgr_k)
!    call pawcprj_set_zero(scf_history%cprj(:,:,ind2))
     call pawcprj_copy(scf_history%cprj(:,:,ind1),scf_history%cprj(:,:,ind2))
   end if
 else

!From 2nd step

!  Init parallelism
   me_distrb=mpi_enreg%me_kpt
   if (mpi_enreg%paral_kgb==1.or.mpi_enreg%paralbd==1) then
     spaceComm_band=mpi_enreg%comm_band
     nprocband=mpi_enreg%nproc_band
   else
     spaceComm_band=xmpi_comm_self
     nprocband=1
   end if

!  For the moment sequential part only
   nprocband=1

!  Additional statements if band-fft parallelism
   if (nprocband>1) then
     ABI_ALLOCATE(npw_block,(nprocband))
     ABI_ALLOCATE(npw_disp,(nprocband))
     ABI_ALLOCATE(bufsize,(nprocband))
     ABI_ALLOCATE(bufdisp,(nprocband))
     ABI_ALLOCATE(bufsize_wf,(nprocband))
     ABI_ALLOCATE(bufdisp_wf,(nprocband))
   end if

   icg=0
   ibg=0

   if(usepaw==1) then
!    WARNING: THIS SECTION IS USELESS; NOW cprj CAN BE READ FROM SCFCV
     call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred_old)
     ABI_DATATYPE_ALLOCATE(cprj,(dtset%natom,scf_history%mcprj))
     call pawcprj_alloc(cprj,0,dimcprj)
     iatom=0 ; iorder_cprj=0
     ABI_ALLOCATE(ylmgr_k,(dtset%mpw,3,0))
     call ctocprj(atindx,cg,1,cprj,gmet,gprimd,iatom,0,iorder_cprj,&
&     dtset%istwfk,kg,dtset%kptns,mcg,scf_history%mcprj,dtset%mgfft,&
&     dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,&
&     nattyp,dtset%nband,dtset%natom,ngfft,dtset%nkpt,dtset%nloalg,&
&     npwarr,dtset%nspinor,dtset%nsppol,dtset%ntypat,dtset%paral_kgb,&
&     ph1d,psps,rmet,dtset%typat,ucvol,0,xred_old,&
&     ylm,ylmgr_k)
     ABI_DEALLOCATE(ylmgr_k)
   end if  ! end usepaw=1

!  LOOP OVER SPINS
   my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
   do isppol=1,dtset%nsppol

!    BIG FAT k POINT LOOP
     do ikpt=1,dtset%nkpt

!      Select k point to be treated by this proc
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) cycle

       istwf_k=dtset%istwfk(ikpt)

!      Retrieve number of plane waves
       npw_k=npwarr(ikpt)
       if (nprocband>1) then
!        Special treatment for band-fft //
         call xmpi_allgather(npw_k,npw_block,spaceComm_band,ierr)
         npw_nk=sum(npw_block);npw_disp(1)=0
         do ii=2,nprocband
           npw_disp(ii)=npw_disp(ii-1)+npw_block(ii-1)
         end do
       else
         npw_nk=npw_k
       end if

!      Allocate arrays for a wave-function (or a block of WFs)
       ABI_ALLOCATE(cwavef,(2,npw_nk*my_nspinor))
       ABI_ALLOCATE(cwavef1,(2,npw_nk*my_nspinor))
       if (nprocband>1) then
         isize=2*my_nspinor;bufsize(:)=isize*npw_block(:);bufdisp(:)=isize*npw_disp(:)
         isize=2*my_nspinor*npw_k;bufsize_wf(:)=isize
         do ii=1,nprocband
           bufdisp_wf(ii)=(ii-1)*isize
         end do
       end if

!      Subspace alignment

!      Loop over bands or blocks of bands
       nblockbd=nband_k/nprocband
       icgb=icg

       if(usepaw==1) then
         ABI_DATATYPE_ALLOCATE( cprj_k,(dtset%natom,my_nspinor*nblockbd))
         call pawcprj_alloc(cprj_k,cprj(1,1)%ncpgr,dimcprj)
         call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,1,isppol,dtset%mband,&
&         dtset%mkmem,dtset%natom,nblockbd,nblockbd,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         ABI_DATATYPE_ALLOCATE( cprj_k1,(dtset%natom,my_nspinor*nblockbd))
         call pawcprj_alloc(cprj_k1,scf_history%cprj(1,1,ind1)%ncpgr,dimcprj)
         call pawcprj_get(atindx1,cprj_k1,scf_history%cprj(:,:,ind1),dtset%natom,1,ibg,ikpt,1,isppol,&
&         dtset%mband,dtset%mkmem,dtset%natom,nblockbd,nblockbd,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         ABI_DATATYPE_ALLOCATE( cprj_k2,(dtset%natom,my_nspinor*nblockbd))
         call pawcprj_alloc(cprj_k2,scf_history%cprj(1,1,ind2)%ncpgr,dimcprj)
         call pawcprj_get(atindx1,cprj_k2,scf_history%cprj(:,:,ind2),dtset%natom,1,ibg,ikpt,1,isppol,&
&         dtset%mband,dtset%mkmem,dtset%natom,nblockbd,nblockbd,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       end if  !end usepaw=1

       ABI_ALLOCATE(unm,(2,nblockbd,nblockbd))
       unm=zero
       icgb2=0

       do iblockbd=1,nblockbd
         iband_min=1+(iblockbd-1)*nprocband
         iband_max=iblockbd*nprocband

         if(xmpi_paral==1.and.mpi_enreg%paral_kgb/=1) then
           if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband_min,iband_max,isppol,me_distrb)) cycle
         end if

!        Extract wavefunction information
         if (nprocband>1) then
!          Special treatment for band-fft //
           ABI_ALLOCATE(cwavef_tmp,(2,npw_k*my_nspinor*nprocband))
           do ig=1,npw_k*my_nspinor*nprocband
             cwavef_tmp(1,ig)=cg(1,ig+icgb)
             cwavef_tmp(2,ig)=cg(2,ig+icgb)
           end do
           call xmpi_alltoallv(cwavef_tmp,bufsize_wf,bufdisp_wf,cwavef,bufsize,bufdisp,spaceComm_band,ierr)
           ABI_DEALLOCATE(cwavef_tmp)
         else
           do ig=1,npw_k*my_nspinor
             cwavef(1,ig)=cg(1,ig+icgb)
             cwavef(2,ig)=cg(2,ig+icgb)
           end do
         end if

         icgb1=icg

         do iblockbd1=1,nblockbd
           iband_min1=1+(iblockbd1-1)*nprocband
           iband_max1=iblockbd1*nprocband

           if(xmpi_paral==1.and.mpi_enreg%paral_kgb/=1) then
             if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband_min1,iband_max1,isppol,me_distrb)) cycle
           end if

!          Extract wavefunction information

           if (nprocband>1) then
!            Special treatment for band-fft //
             ABI_ALLOCATE(cwavef_tmp,(2,npw_k*my_nspinor*nprocband))
             do ig=1,npw_k*my_nspinor*nprocband
               cwavef_tmp(1,ig)=scf_history%cg(1,ig+icgb1,ind1)
               cwavef_tmp(2,ig)=scf_history%cg(2,ig+icgb1,ind1)
             end do
             call xmpi_alltoallv(cwavef_tmp,bufsize_wf,bufdisp_wf,cwavef1,bufsize,bufdisp,spaceComm_band,ierr)
             ABI_DEALLOCATE(cwavef_tmp)
           else
             do ig=1,npw_k*my_nspinor
               cwavef1(1,ig)=scf_history%cg(1,ig+icgb1,ind1)
               cwavef1(2,ig)=scf_history%cg(2,ig+icgb1,ind1)
             end do
           end if

!          Calculate Unm=<psi_nk(t)|S|psi_mk(t-dt)>
           call dotprod_g(dotr,doti,istwf_k,npw_k*my_nspinor,2,cwavef,cwavef1,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
           if(usepaw==1) then
             ia =0
             do itypat=1,ntypat
               do iat=1+ia,nattyp(itypat)+ia
                 do ilmn1=1,pawtab(itypat)%lmn_size
                   do ilmn2=1,ilmn1
                     klmn=((ilmn1-1)*ilmn1)/2+ilmn2
                     dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd1)%cp(1,ilmn2)+&
&                     cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd1)%cp(2,ilmn2))
                     doti=doti+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd1)%cp(2,ilmn2)-&
&                     cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd1)%cp(1,ilmn2))
                   end do
                   do ilmn2=ilmn1+1,pawtab(itypat)%lmn_size
                     klmn=((ilmn2-1)*ilmn2)/2+ilmn1
                     dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd1)%cp(1,ilmn2)+&
&                     cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd1)%cp(2,ilmn2))
                     doti=doti+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd1)%cp(2,ilmn2)-&
&                     cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd1)%cp(1,ilmn2))
                   end do
                 end do
               end do
               ia=ia+nattyp(itypat)
             end do
           end if
!          unm(1,iblockbd,iblockbd1)=dotr
!          unm(2,iblockbd,iblockbd1)=doti
           unm(1,iblockbd1,iblockbd)=dotr
           unm(2,iblockbd1,iblockbd)=doti
!          End loop over bands iblockbd1
           icgb1=icgb1+npw_k*my_nspinor*nprocband

         end do

!        End loop over bands iblockbd
         icgb2=icgb2+npw_k*my_nspinor*nprocband
         icgb=icgb+npw_k*my_nspinor*nprocband
       end do

!      write(std_out,*) 'UNM'
!      do iblockbd=1,nblockbd
!      write(std_out,11) (unm(1,iblockbd,iblockbd1),unm(2,iblockbd,iblockbd1),iblockbd1=1,nblockbd)
!      end do
!      11 format(12(1x,f9.5),a)
!      Compute A=tU^*U
       ABI_ALLOCATE(unm_tmp,(nblockbd,nblockbd))
       ABI_ALLOCATE(anm_tmp,(nblockbd,nblockbd))
       ABI_ALLOCATE(anm,(nblockbd*(nblockbd+1)))
       unm_tmp(:,:)=cmplx(unm(1,:,:),unm(2,:,:),kind=dp)
       call zgemm('C','N',nblockbd,nblockbd,nblockbd,dcmplx(1._dp), unm_tmp,nblockbd, &
&       unm_tmp,nblockbd,dcmplx(0._dp),anm_tmp,nblockbd)
       do iblockbd=1,nblockbd
         do iblockbd1=iblockbd,nblockbd
           ii=iblockbd1*(iblockbd1-1)+2*(iblockbd-1)+1
           anm(ii)=real(anm_tmp(iblockbd,iblockbd1))
           anm(ii+1)=aimag(anm_tmp(iblockbd,iblockbd1))
         end do
       end do
       call hermit(anm,anm,ierr,nblockbd)
!      aa=dcmplx(0._dp)
!      do iblockbd=1,nblockbd
!      aa=aa+conjg(unm_tmp(iblockbd,1))*unm_tmp(iblockbd,1)
!      end do
!      write(std_out,*) 'tU*U', aa
!      write(std_out,*) 'ANM_tmp'
!      do iblockbd=1,nblockbd
!      write(std_out,11) (anm_tmp(iblockbd,iblockbd1),iblockbd1=1,nblockbd)
!      end do
!      write(std_out,*) 'ANM'
!      do iblockbd=1,nblockbd*(nblockbd+1)
!      write(std_out,11) anm(iblockbd)
!      end do

!      Diagonalize A
       ABI_ALLOCATE(eig,(nblockbd))
       ABI_ALLOCATE(evec,(2*nblockbd,nblockbd))
       ABI_ALLOCATE(zhpev1,(2,2*nblockbd-1))
       ABI_ALLOCATE(zhpev2,(3*nblockbd-2))
       call zhpev('V','U',nblockbd,anm,eig,evec,nblockbd,zhpev1,&
&       zhpev2,ierr)
       ABI_DEALLOCATE(anm)
       ABI_DEALLOCATE(zhpev1)
       ABI_DEALLOCATE(zhpev2)
!      aa=dcmplx(0._dp)
!      do iblockbd=1,nblockbd
!      aa=aa+anm_tmp(1,iblockbd)*cmplx(evec((2*iblockbd-1),1),evec(2*iblockbd,1),kind=dp)
!      end do
!      write(std_out,*) 'EIG', aa, eig(1)*evec(1,1),eig(1)*evec(2,1)

!      Compute A'=evec*tU^/sqrt(eig)
       call zgemm('C','C',nblockbd,nblockbd,nblockbd,dcmplx(1._dp),evec,nblockbd, &
&       unm_tmp,nblockbd,dcmplx(0._dp),anm_tmp,nblockbd)
       do iblockbd=1,nblockbd
         eigval=dsqrt(eig(iblockbd))
         do iblockbd1=1,nblockbd
           anm_tmp(iblockbd,iblockbd1)=anm_tmp(iblockbd,iblockbd1)/eigval
         end do
       end do

!      Compute tA^A'to come back to the initial subspace for the cg's

       call zgemm('N','N',nblockbd,nblockbd,nblockbd,dcmplx(1._dp),evec,nblockbd, &
&       anm_tmp,nblockbd,dcmplx(0._dp),unm_tmp,nblockbd)
       anm_tmp=unm_tmp
!      write(std_out,*) 'ANM_tmp'
!      do iblockbd=1,nblockbd
!      write(std_out,11) (anm_tmp(iblockbd,iblockbd1),iblockbd1=1,nblockbd)
!      end do

!      Wavefunction alignment (istwfk=1 ?)
       ABI_ALLOCATE(work,(2,npw_nk*my_nspinor*nblockbd))
       ABI_ALLOCATE(work1,(2,my_nspinor*nblockbd*npw_nk))
       work1(:,:)=scf_history%cg(:,icg+1:icg+my_nspinor*nblockbd*npw_nk,ind1)
       call zgemm('N','N',npw_nk*my_nspinor,nblockbd,nblockbd,dcmplx(1._dp), &
&       work1,npw_nk*my_nspinor, &
&       anm_tmp,nblockbd,dcmplx(0._dp),work,npw_nk*my_nspinor)
       scf_history%cg(:,1+icg:npw_nk*my_nspinor*nblockbd+icg,ind1)=work(:,:)

       work1(:,:)=scf_history%cg(:,icg+1:icg+my_nspinor*nblockbd*npw_nk,ind2)
       call zgemm('N','N',npw_nk*my_nspinor,nblockbd,nblockbd,dcmplx(1._dp), &
&       work1,npw_nk*my_nspinor, &
&       anm_tmp,nblockbd,dcmplx(0._dp),work,npw_nk*my_nspinor)
       scf_history%cg(:,1+icg:npw_nk*my_nspinor*nblockbd+icg,ind2)=work(:,:)
       ABI_DEALLOCATE(work1)
!      If paw, must also align cprj:
       if (usepaw==1) then
!        New version (MT):
         ABI_DATATYPE_ALLOCATE(cprj_k3,(dtset%natom,my_nspinor))
         ABI_DATATYPE_ALLOCATE(cprj_k4,(dtset%natom,my_nspinor))
         call pawcprj_alloc(cprj_k3,cprj_k1(1,1)%ncpgr,dimcprj)
         call pawcprj_alloc(cprj_k4,cprj_k2(1,1)%ncpgr,dimcprj)
         ABI_ALLOCATE(al,(2,nblockbd))
         do iblockbd=1,nblockbd
           ii=(iblockbd-1)*my_nspinor
           do iblockbd1=1,nblockbd
             al(1,iblockbd1)=real (anm_tmp(iblockbd,iblockbd1))
             al(2,iblockbd1)=aimag(anm_tmp(iblockbd,iblockbd1))
           end do
           call pawcprj_lincom(al,cprj_k1,cprj_k3,nblockbd)
           call pawcprj_lincom(al,cprj_k2,cprj_k4,nblockbd)
           call pawcprj_copy(cprj_k3,cprj_k1(:,ii+1:ii+my_nspinor))
           call pawcprj_copy(cprj_k4,cprj_k2(:,ii+1:ii+my_nspinor))
         end do
         ABI_DEALLOCATE(al)
!        Old version (FJ):
!        allocate( cprj_k3(dtset%natom,my_nspinor*nblockbd))
!        call pawcprj_alloc(cprj_k3,cprj_k1(1,1)%ncpgr,dimcprj)
!        allocate( cprj_k4(dtset%natom,my_nspinor*nblockbd))
!        call pawcprj_alloc(cprj_k4,cprj_k2(1,1)%ncpgr,dimcprj)
!        beta(1)=one;beta(2)=zero
!        do iblockbd=1,nblockbd*my_nspinor
!        do iblockbd1=1,nblockbd*my_nspinor
!        alpha(1)=real(anm_tmp(iblockbd,iblockbd1));alpha(2)=aimag(anm_tmp(iblockbd,iblockbd1))
!        call pawcprj_zaxpby(alpha,beta,cprj_k1(:,iblockbd1:iblockbd1),cprj_k3(:,iblockbd:iblockbd))
!        call pawcprj_zaxpby(alpha,beta,cprj_k2(:,iblockbd1:iblockbd1),cprj_k4(:,iblockbd:iblockbd))
!        end do
!        end do
!        call pawcprj_copy(cprj_k3,cprj_k1)
!        call pawcprj_copy(cprj_k4,cprj_k2)

         call pawcprj_free(cprj_k3)
         call pawcprj_free(cprj_k4)
         ABI_DATATYPE_DEALLOCATE(cprj_k3)
         ABI_DATATYPE_DEALLOCATE(cprj_k4)
       end if
       ABI_DEALLOCATE(anm_tmp)
       ABI_DEALLOCATE(unm_tmp)
       ABI_DEALLOCATE(work)

!      Wavefunction extrapolation
       ibd=0
       inc=npw_nk*my_nspinor
       ABI_ALLOCATE(deltawf2,(2,npw_nk*my_nspinor))
       ABI_ALLOCATE(wf1,(2,npw_nk*my_nspinor))
       ABI_ALLOCATE(deltawf1,(2,npw_nk*my_nspinor))
       do iblockbd=1,nblockbd
         deltawf2(:,:)=scf_history%cg(:,1+icg+ibd:icg+ibd+inc,ind2)
         wf1(:,:)=scf_history%cg(:,1+icg+ibd:icg+ibd+inc,ind1)
!        wf1(2,1)=zero;deltawf2(2,1)=zero

         call dotprod_g(dotr,doti,istwf_k,npw_nk*my_nspinor,2,cg(:,icg+1+ibd:ibd+icg+inc),cg(:,icg+1+ibd:ibd+icg+inc),&
&         mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
         call dotprod_g(dotr1,doti1,istwf_k,npw_nk*my_nspinor,2,cg(:,icg+1+ibd:ibd+icg+inc),wf1,&
&         mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
         if(usepaw==1) then
           ia =0
           do itypat=1,ntypat
             do iat=1+ia,nattyp(itypat)+ia
               do ilmn1=1,pawtab(itypat)%lmn_size
                 do ilmn2=1,ilmn1
                   klmn=((ilmn1-1)*ilmn1)/2+ilmn2
                   dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k(iat,iblockbd)%cp(1,ilmn2)+&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k(iat,iblockbd)%cp(2,ilmn2))
                   doti=doti+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k(iat,iblockbd)%cp(2,ilmn2)-&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k(iat,iblockbd)%cp(1,ilmn2))
                   dotr1=dotr1+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd)%cp(1,ilmn2)+&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd)%cp(2,ilmn2))
                   doti1=doti1+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd)%cp(2,ilmn2)-&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd)%cp(1,ilmn2))
                 end do
                 do ilmn2=ilmn1+1,pawtab(itypat)%lmn_size
                   klmn=((ilmn2-1)*ilmn2)/2+ilmn1
                   dotr=dotr+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k(iat,iblockbd)%cp(1,ilmn2)+&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k(iat,iblockbd)%cp(2,ilmn2))
                   doti=doti+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k(iat,iblockbd)%cp(2,ilmn2)-&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k(iat,iblockbd)%cp(1,ilmn2))
                   dotr1=dotr1+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd)%cp(1,ilmn2)+&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd)%cp(2,ilmn2))
                   doti1=doti1+pawtab(itypat)%sij(klmn)*(cprj_k(iat,iblockbd)%cp(1,ilmn1)*cprj_k1(iat,iblockbd)%cp(2,ilmn2)-&
&                   cprj_k(iat,iblockbd)%cp(2,ilmn1)*cprj_k1(iat,iblockbd)%cp(1,ilmn2))
                 end do
               end do
             end do
             ia=ia+nattyp(itypat)
           end do
         end if
         dotr=sqrt(dotr**2+doti**2)
         dotr1=sqrt(dotr1**2+doti1**2)
         write(std_out,*)'DOTR, DOTR1',dotr,dotr1
         dotr=dotr1/dotr
         write(std_out,*)'DOTR',dotr
         deltawf1=zero
         if(dotr>=0.9d0) then
           deltawf1(:,:)=cg(:,icg+1+ibd:ibd+icg+inc)-wf1(:,:)
           if(usepaw==1) then
             alpha(1)=one;alpha(2)=zero
             beta(1)=-one;beta(2)=zero
             ia =0
             call pawcprj_zaxpby(alpha,beta,cprj_k(:,iblockbd:iblockbd),cprj_k1(:,iblockbd:iblockbd))
           end if
           istep1=istep
         else
           istep1=1
         end if
         scf_history%cg(:,1+icg+ibd:icg+ibd+inc,ind1)=cg(:,icg+1+ibd:ibd+icg+inc)
         scf_history%cg(:,1+icg+ibd:icg+ibd+inc,ind2)=deltawf1(:,:)
         if(usepaw==1) then
           call pawcprj_put(atindx1,cprj_k,scf_history%cprj(:,:,ind1),dtset%natom,1,ibg,ikpt,1,isppol,&
&           dtset%mband,dtset%mkmem,dtset%natom,nblockbd,nblockbd,dimcprj,my_nspinor,dtset%nsppol,0,&
&           mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
           call pawcprj_put(atindx1,cprj_k1,scf_history%cprj(:,:,ind2),dtset%natom,1,ibg,ikpt,1,isppol,&
&           dtset%mband,dtset%mkmem,dtset%natom,nblockbd,nblockbd,dimcprj,my_nspinor,dtset%nsppol,0,&
&           mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
         end if

!        if(istep1>=3) then
         cg(:,icg+1+ibd:ibd+icg+inc)=cg(:,icg+1+ibd:ibd+icg+inc)+scf_history%alpha*deltawf1(:,:) &
&         +scf_history%beta *deltawf2(:,:)

!        to be used later
!        if(usepaw==1) then
!        alpha(2)=zero
!        beta(1)=one;beta(2)=zero
!        alpha(1)=scf_history%alpha
!        call pawcprj_zaxpby(alpha,beta,cprj_k1(:,iblockbd:iblockbd),cprj_k(:,iblockbd:iblockbd))
!        alpha(1)=scf_history%beta
!        call pawcprj_zaxpby(alpha,beta,cprj_k2(:,iblockbd:iblockbd),cprj_k(:,iblockbd:iblockbd))
!        call pawcprj_put(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,1,isppol,&
!        &    dtset%mband,dtset%mkmem,dtset%natom,nblockbd,nblockbd,dimcprj,my_nspinor,dtset%nsppol,0,&
!        &    mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
!        end if
!        else if (istep1==2) then
!          cg(:,icg+1+ibd:ibd+icg+inc)=cg(:,icg+1+ibd:ibd+icg+inc)+scf_history%alpha*deltawf1(:,:)+scf_history%beta*wf1(:,:)
!       !     cg(:,icg+1+ibd:ibd+icg+inc)=cg(:,icg+1+ibd:ibd+icg+inc)+deltawf1(:,:)
!        if(usepaw==1) then
!        alpha(2)=zero
!        beta(1)=one;beta(2)=zero
!        alpha(1)=scf_history%alpha
!        call pawcprj_zaxpby(alpha,beta,cprj_k1(:,iblockbd:iblockbd),cprj_k(:,iblockbd:iblockbd))
!        alpha(1)=scf_history%beta
!        call pawcprj_zaxpby(alpha,beta,cprj_k2(:,iblockbd:iblockbd),cprj_k(:,iblockbd:iblockbd))
!        call pawcprj_put(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,1,isppol,&
!        &    dtset%mband,dtset%mkmem,dtset%natom,nblockbd,nblockbd,dimcprj,my_nspinor,dtset%nsppol,0,&
!        &    mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
!        end if
!        end if
         ibd=ibd+inc
       end do ! end loop on iblockbd

       ABI_DEALLOCATE(deltawf1)
       ABI_DEALLOCATE(deltawf2)
       ABI_DEALLOCATE(wf1)
       ABI_DEALLOCATE(cwavef)
       ABI_DEALLOCATE(cwavef1)
       ABI_DEALLOCATE(eig)
       ABI_DEALLOCATE(evec)
       ABI_DEALLOCATE(unm)
       if(usepaw==1) then
         call pawcprj_free(cprj_k)
         ABI_DATATYPE_DEALLOCATE(cprj_k)
         call pawcprj_free(cprj_k1)
         ABI_DATATYPE_DEALLOCATE(cprj_k1)
         call pawcprj_free(cprj_k2)
         ABI_DATATYPE_DEALLOCATE(cprj_k2)
       end if

       ibg=ibg+my_nspinor*nband_k
       icg=icg+my_nspinor*nband_k*npw_k

!      End big k point loop
     end do
!    End loop over spins
   end do

   if(usepaw==1) then
     call pawcprj_free(cprj)
     ABI_DATATYPE_DEALLOCATE(cprj)
   end if
   if (nprocband>1) then
     ABI_DEALLOCATE(npw_block)
     ABI_DEALLOCATE(npw_disp)
     ABI_DEALLOCATE(bufsize)
     ABI_DEALLOCATE(bufdisp)
     ABI_DEALLOCATE(bufsize_wf)
     ABI_DEALLOCATE(bufdisp_wf)
   end if

 end if ! istep>=2

 if (usepaw==1) then
   ABI_DEALLOCATE(dimcprj)
 end if

end subroutine extrapwf
!!***


!!****f* ABINIT/extrapwf_biortho
!!
!! NAME
!! extrapwf_biortho
!!
!! FUNCTION
!! Extrapolate wavefunctions for new ionic positions
!! from values of wavefunctions of previous SCF cycle.
!! Use biorthogonal algorithm proposed XG
!!
!! INPUTS
!!  atindx1(dtset%natom)=index table for atoms, inverse of atindx
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  istep=number of call the routine
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of cprj array
!!  mpi_enreg=information about MPI parallelization
!!  nattyp(dtset%ntypat)=number of atoms of each type in cell.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  pawtab(dtset%ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! SIDE EFFECTS
!!  cg(2,mcg)= plane wave wavefunction coefficient
!!                          Value from previous SCF cycle is input and stored in some form
!!                          Extrapolated value is output
!!  cprj(natom,mcprj) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors
!!                          Value from previous SCF cycle is input and stored in some form
!!                          Extrapolated value is output
!!  scf_history_wf <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!
!! PARENTS
!!      m_extraprho
!!
!! CHILDREN
!!      cgcprj_cholesky,dotprod_set_cgcprj,lincom_cgcprj,pawcprj_alloc
!!      pawcprj_axpby,pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_put,zgesv
!!
!! SOURCE

 subroutine extrapwf_biortho(atindx1,cg,cprj,dtset,istep,mcg,mcprj,mpi_enreg,&
& nattyp,npwarr,pawtab,scf_history_wf)

 !use m_scf_history
 use m_cgcprj,  only : dotprod_set_cgcprj,cgcprj_cholesky,lincom_cgcprj

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,mcg,mcprj
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(scf_history_type),intent(inout) :: scf_history_wf
!arrays
 integer,intent(in) :: atindx1(dtset%natom),nattyp(dtset%ntypat)
 integer,intent(in) :: npwarr(dtset%nkpt)
 real(dp), intent(inout) :: cg(2,mcg)
 type(pawcprj_type),intent(inout) :: cprj(dtset%natom,mcprj)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: hermitian
 integer :: ibdmix,ibg,ibg_hist,icg,icg_hist !,iband
 integer :: ierr,ikpt,indh,ind1,ind2,ind1new,inplace
 integer :: isppol,istwf_k,kk,me_distrb,mband,my_nspinor,mcprj_k
 integer :: nband_k,nbdmix,nbdmax,npw_k,ntypat
 integer :: spaceComm_band,usepaw
 real(dp) :: alpha,beta !,dotr,doti

!arrays
 integer,allocatable :: ipiv(:),dimcprj(:)
 real(dp),allocatable ::psi_ortho(:,:),mmn(:,:,:)
 real(dp),allocatable :: smn(:,:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kh(:,:)

! *************************************************************************

 if (istep==0) return

 ntypat=dtset%ntypat
 usepaw=dtset%usepaw
 mband=dtset%mband
 nbdmax=dtset%mband
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 me_distrb=mpi_enreg%me_kpt
 spaceComm_band=xmpi_comm_self

!scf_history_wf%alpha contains dtset%wfmix
 alpha=scf_history_wf%alpha
 beta=scf_history_wf%beta
 ind1=scf_history_wf%hindex(1)
 ind2=scf_history_wf%hindex(2)
 ind1new=scf_history_wf%hindex(3)
 icg=0
 icg_hist=0
 ibg=0
 ibg_hist=0

!Useful array
 ABI_ALLOCATE(dimcprj,(dtset%natom))
 if (usepaw==1) then
   call pawcprj_getdim(dimcprj,dtset%natom,nattyp,ntypat,dtset%typat,pawtab,'O')
 end if

 if(istep==1)then
   do indh=1,scf_history_wf%history_size
     call pawcprj_alloc(scf_history_wf%cprj(:,:,indh),0,dimcprj)
   end do
 end if

 mcprj_k=my_nspinor*nbdmax
 ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,mcprj_k))
 ABI_DATATYPE_ALLOCATE(cprj_kh,(dtset%natom,mcprj_k))

 if(usepaw==1) then
   call pawcprj_alloc(cprj_k,0,dimcprj)
   call pawcprj_alloc(cprj_kh,0,dimcprj)
 end if
 ABI_ALLOCATE(smn,(2,nbdmax,nbdmax))
 ABI_ALLOCATE(mmn,(2,nbdmax,nbdmax))

!Explanation for the index for the wavefunction stored in scf_history_wf
!The reference is the cg+cprj output after the wf optimization at istep 1.
!For wavefunction mixing for molecular dynamics, we use the same mixing as for the density in extraprho. To keep the same indexes,
! we choose to take indh=3 for the reference.

!First step
 if (istep==1) then

   indh=3   ! This input wavefunction is the reference

!  LOOP OVER SPINS
   do isppol=1,dtset%nsppol

!    BIG FAT k POINT LOOP
     do ikpt=1,dtset%nkpt

!      Select k point to be treated by this proc
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       nbdmix=min(nband_k,nbdmax)

       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) cycle

       npw_k=npwarr(ikpt)

       scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,indh)=cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)

       if(usepaw==1) then
!        scf_history_wf%cprj(:,ibg_hist+1:ibg_hist+my_nspinor*nbdmix,1)=cprj(:,ibg+1:ibg+my_nspinor*nbdmix)
         call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,0,isppol,mband,&
&         dtset%mkmem,dtset%natom,nbdmax,nbdmix,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         call pawcprj_put(atindx1,cprj_k,scf_history_wf%cprj(:,:,indh),dtset%natom,1,ibg_hist,ikpt,0,isppol,&
&         mband,dtset%mkmem,dtset%natom,nbdmax,nbdmix,dimcprj,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
       end if

!      Update the counters
       ibg=ibg+my_nspinor*nband_k
       ibg_hist=ibg_hist+my_nspinor*nbdmix
       icg=icg+my_nspinor*nband_k*npw_k
       icg_hist=icg_hist+my_nspinor*nbdmix*npw_k

     end do
   end do

 else
!  From istep==2
   if (istep==2) ind1=3
   if (istep==3) ind2=3
!  biorthogonalization
   indh=3   ! This input wavefunction is the reference

!  LOOP OVER SPINS
   do isppol=1,dtset%nsppol

!    BIG FAT k POINT LOOP
     do ikpt=1,dtset%nkpt

!      Select k point to be treated by this proc
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       nbdmix=min(nband_k,nbdmax)
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) cycle
       istwf_k=dtset%istwfk(ikpt)
       npw_k=npwarr(ikpt)
       ABI_ALLOCATE(psi_ortho,(2,npw_k*my_nspinor*nbdmix))
       psi_ortho=zero
!      Biorthogonalization

       if(usepaw==1) then
         call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,0,isppol,mband,&
&         dtset%mkmem,dtset%natom,nbdmax,nbdmix,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         call pawcprj_get(atindx1,cprj_kh,scf_history_wf%cprj(:,:,indh),dtset%natom,1,ibg_hist,ikpt,0,isppol,&
&         mband,dtset%mkmem,dtset%natom,nbdmax,nbdmix,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       end if  !end usepaw=1

       hermitian=0

       call dotprod_set_cgcprj(atindx1,scf_history_wf%cg(:,:,indh),cg,cprj_kh,cprj_k,dimcprj,hermitian,&
&         0,0,icg_hist,icg,ikpt,isppol,istwf_k,nbdmax,mcg,mcg,mcprj_k,mcprj_k,dtset%mkmem,&
&         mpi_enreg,dtset%natom,nattyp,nbdmix,nbdmix,npw_k,my_nspinor,dtset%nsppol,ntypat,pawtab,smn(:,1:nbdmix,1:nbdmix),usepaw)

!      Invert S matrix, that is NOT hermitian.
!      Calculate M=S^-1
       mmn=zero
       do kk=1,nbdmix
         mmn(1,kk,kk)=one
       end do

       ABI_ALLOCATE(ipiv,(nbdmix))
!      The smn is destroyed by the following inverse call
       call zgesv(nbdmix,nbdmix,smn,nbdmax,ipiv,mmn,nbdmax,ierr)
       ABI_DEALLOCATE(ipiv)
!DEBUG
       if(ierr/=0)then
         MSG_ERROR(' The call to cgesv general inversion routine failed')
       end if
!ENDDEBUG

!      The M matrix is used to compute the biorthogonalized set of wavefunctions, and to store it at the proper place
       inplace=0
       call lincom_cgcprj(mmn(:,1:nbdmix,1:nbdmix),cg,cprj_k,dimcprj,&
&         icg,inplace,mcg,mcprj_k,dtset%natom,nbdmix,nbdmix,npw_k,my_nspinor,usepaw,&
&         cgout=psi_ortho,cprjout=cprj_kh,icgout=0)

!!!TEST
!      if (usepaw==0) then
!        do iband=1,nband_k
!          call dotprod_g(dotr,doti,istwf_k,npw_k,2,scf_history_wf%cg(:,icg+1+my_nspinor*npw_k:icg+2*my_nspinor*npw_k,indh),&
!&            psi_ortho(:,1+(iband-1)*my_nspinor*npw_k:iband*my_nspinor*npw_k),mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
!          write(80+mpi_enreg%me,*) dotr,doti
!          flush(80+mpi_enreg%me)
!        end do
!      else
!        hermitian=0
!        call dotprod_set_cgcprj(atindx1,scf_history_wf%cg(:,:,indh),psi_ortho,scf_history_wf%cprj(:,:,indh)&
!&         ,cprj_kh,dimcprj,hermitian,0,0,icg_hist,icg_hist,ikpt,isppol,istwf_k,nbdmax,mcg,mcg,mcprj_k,mcprj_k,dtset%mkmem,&
!&         mpi_enreg,dtset%natom,nattyp,nbdmix,nbdmix,npw_k,my_nspinor,dtset%nsppol,ntypat,pawtab,smn(:,1:nbdmix,1:nbdmix),usepaw)
!        write(90+mpi_enreg%me,*) smn
!        flush(90+mpi_enreg%me)
!      end if
!!!TEST
!      The biorthogonalised set of wavefunctions is now stored at the proper place

!       cg(:,icg+1:icg+my_nspinor*npw_k*nband_k)=zero

!      psi(t+dt) <- psi(t) + alpha.psi(t)
       cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)=(one+alpha)*psi_ortho(:,1:my_nspinor*npw_k*nbdmix)
       if(usepaw==1) then
         do ibdmix=1,nbdmix
           call pawcprj_axpby(one+alpha,zero,cprj_kh(:,ibdmix:ibdmix),cprj_k(:,ibdmix:ibdmix))
         end do ! end loop on ibdmix
       end if
!      psi(t+dt) <- -alpha.psi(t-dt) + beta.psi(t-dt)
       if (abs(beta-alpha)>tol14.and.ind1>0) then
         cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)=cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)+&
&             (beta-alpha)*scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind1)
         if(usepaw==1) then
           do ibdmix=1,nbdmix
             call pawcprj_axpby(beta-alpha,one,scf_history_wf%cprj(:,ibg_hist+ibdmix:ibg_hist+ibdmix,ind1),&
&                                 cprj_k(:,ibdmix:ibdmix))
           end do ! end loop on ibdmix
         end if
       end if

!      psi(t+dt) <- -beta.psi(t-2dt)
       if (abs(beta)>tol14.and.ind2>0) then
         cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)=cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix)&
&                               -beta*scf_history_wf%cg(:,icg_hist+1:icg_hist+my_nspinor*npw_k*nbdmix,ind2)
         if(usepaw==1) then
           do ibdmix=1,nbdmix
             call pawcprj_axpby(-beta,one,scf_history_wf%cprj(:,ibg_hist+ibdmix:ibg_hist+ibdmix,ind2),cprj_k(:,ibdmix:ibdmix))
           end do ! end loop on ibdmix
         end if
       end if

!      Store psi(t) in history
       scf_history_wf%cg(:,icg+1:icg+my_nspinor*npw_k*nbdmix,ind1new)=psi_ortho(:,1:my_nspinor*npw_k*nbdmix)
       if(usepaw==1) then
         call pawcprj_put(atindx1,cprj_kh,scf_history_wf%cprj(:,:,ind1new),dtset%natom,1,ibg_hist,ikpt,0,isppol,&
&         nbdmix,dtset%mkmem,dtset%natom,nbdmax,nbdmix,dimcprj,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
       end if

!      Back to usual orthonormalization for the cg and cprj_k
       call cgcprj_cholesky(atindx1,cg,cprj_k,dimcprj,icg,ikpt,isppol,istwf_k,mcg,mcprj_k,dtset%mkmem,&
&       mpi_enreg,dtset%natom,nattyp,nbdmax,npw_k,my_nspinor,dtset%nsppol,ntypat,pawtab,usepaw)

!      Need to transfer cprj_k to cprj
       if(usepaw==1) then
         call pawcprj_put(atindx1,cprj_k,cprj,dtset%natom,1,ibg,ikpt,0,isppol,&
&         mband,dtset%mkmem,dtset%natom,nbdmax,nbdmix,dimcprj,my_nspinor,dtset%nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,mpi_comm_band=spaceComm_band,proc_distrb=mpi_enreg%proc_distrb)
       end if

       ibg=ibg+my_nspinor*nband_k
       ibg_hist=ibg_hist+my_nspinor*nbdmix
       icg=icg+my_nspinor*nband_k*npw_k
       icg_hist=icg_hist+my_nspinor*nbdmix*npw_k
       ABI_DEALLOCATE(psi_ortho)
!      End big k point loop
     end do
!    End loop over spins
   end do


 end if !istep>1


 if(usepaw==1) then
   call pawcprj_free(cprj_k)
   call pawcprj_free(cprj_kh)
 end if
 ABI_DATATYPE_DEALLOCATE(cprj_k)
 ABI_DATATYPE_DEALLOCATE(cprj_kh)
 ABI_DEALLOCATE(dimcprj)
 ABI_DEALLOCATE(mmn)
 ABI_DEALLOCATE(smn)



end subroutine extrapwf_biortho
!!***
end module m_extraprho
!!***
