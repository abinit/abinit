!{\src2tex{textfont=tt}}
!!****f* ABINIT/forces
!!
!! NAME
!! forces
!!
!! FUNCTION
!! Assemble gradients of various total energy terms with respect
!! to reduced coordinates, including possible symmetrization,
!! in order to produce forces.
!!     fcart(i,iat) = d(Etot)/(d(r(i,iat)))
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!  dtset <type(dataset_type)>=all input variables in this dataset
!! berryopt    =  4/14: electric field is on -> add the contribution of the
!!                      -ebar_i p_i - Omega/(8*pi) (g^{-1})_ij ebar_i ebar_j  terms to the total energy
!!     = 6/16, or 7/17: electric displacement field is on  -> add the contribution of the
!!                      Omega/(8*pi) (g^{-1})_ij ebar_i ebar_j  terms to the total energy
!!   | efield = cartesian coordinates of the electric field in atomic units
!!   | dfield = cartesian coordinates of the electric displacement field in atomic units
!!   | iatfix(3,natom)=1 for frozen atom along specified direction, 0 for unfrozen
!!   | ionmov=governs the movement of atoms (see help file)
!!   | densfor_pred=governs the mixed electronic-atomic part of the preconditioner
!!   | natom=number of atoms in cell
!!   | nconeq=number of atomic constraint equations
!!   | nspden=number of spin-density components
!!   | nsym=number of symmetries in space group
!!   | prtvol=integer controlling volume of printed output
!!   | typat(natom)=type integer for each atom in cell
!!   | wtatcon(3,natom,nconeq)=weights for atomic constraints
!!  fock <type(fock_type)>= quantities to calculate Fock exact exchange
!!  grewtn(3,natom)=d(Ewald)/d(xred) (hartree)
!!  grnl(3*natom)=gradients of Etot due to nonlocal contributions
!!  grvdw(3,ngrvdw)=gradients of energy due to Van der Waals DFT-D dispersion (hartree)
!!  gsqcut=cutoff value on G**2 for (large) sphere inside FFT box.
!!                       gsqcut=(boxcut**2)*ecut/(2._dp*(Pi**2)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngrvdw=size of grvdw(:,:); can be 0 or natom according to dtset%vdw_xc
!!  ntypat=number of types of atoms
!!  pawtab(ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase (structure factor) array
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=Fourier transform of charge density (bohr^-3)
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  usefock=1 if fock operator is used; 0 otherwise.
!!  vresid(nfft,nspden)=potential residual (if non-collinear magn., only trace of it)
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree) in real space
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)=previous reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  diffor=maximal absolute value of changes in the components of
!!         force between the input and the output.
!!  favg(3)=mean of the forces before correction for translational symmetry
!!  forold(3,natom)=cartesian forces of previous SCF cycle (hartree/bohr)
!!  fred(3,natom)=symmetrized grtn = d(etotal)/d(xred)
!!  gresid(3,natom)=forces due to the residual of the density/potential
!!  grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!  grxc(9+3*natom)=d(Exc)/d(xred) if core charges are used
!!  maxfor=maximal absolute value of the output array force.
!!  synlgr(3,natom)=symmetrized d(enl)/d(xred)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  fcart(3,natom)=forces in cartesian coordinates (Ha/Bohr)
!!    Note : unlike fred, this array has been corrected by enforcing
!!    the translational symmetry, namely that the sum of force
!!    on all atoms is zero.
!!
!! NOTES
!! * Symmetrization of gradients with respect to reduced
!!   coordinates xred is conducted according to the expression
!!   [d(e)/d(t(n,a))]_symmetrized = (1/Nsym) Sum(S) symrec(n,m,S)*
!!                [d(e)/d(t(m,b))]_unsymmetrized
!!   where t(m,b)= (symrel^-1)(m,n)*(t(n,a)-tnons(n)) and tnons
!!   is a possible nonsymmorphic translation.  The label "b" here
!!   refers to the atom which gets rotated into "a" under symmetry "S".
!!   symrel is the symmetry matrix in real space, which is the inverse
!!   transpose of symrec.  symrec is the symmetry matrix in reciprocal
!!   space.  sym_cartesian = R * symrel * R^-1 = G * symrec * G^-1
!!   where the columns of R and G are the dimensional primitive translations
!!   in real and reciprocal space respectively.
!! * Note the use of "symrec" in the symmetrization expression above.
!!
!! PARENTS
!!      afterscfloop,etotfor,forstr
!!
!! CHILDREN
!!      atm2fft,constrf,dgemv,fourdp,fred2fcart,fresid,fresidrsp,metric,mkcore
!!      mklocl,sygrad,timab,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine forces(atindx1,diffor,dtefield,dtset,favg,fcart,fock,forold,fred,gresid,grewtn,&
&                  grhf,grnl,grvdw,grxc,gsqcut,indsym,&
&                  maxfor,mgfft,mpi_enreg,n1xccc,n3xccc,&
&                  nattyp,nfft,ngfft,ngrvdw,ntypat,&
&                  pawtab,ph1d,psps,rhog,rhor,rprimd,symrec,synlgr,usefock,&
&                  vresid,vxc,wvl,wvl_den,xred,&
&                  electronpositron) ! optional argument

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_efield
 use m_errors
 use m_fock,             only : fock_type
 use m_pawtab,           only : pawtab_type
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'forces'
 use interfaces_18_timing
 use interfaces_41_geometry
 use interfaces_53_ffts
 use interfaces_56_xc
 use interfaces_64_psp
 use interfaces_67_common, except_this_one => forces
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,n1xccc,n3xccc,nfft,ngrvdw,ntypat,usefock
 real(dp),intent(in) :: gsqcut
 real(dp),intent(out) :: diffor,maxfor
 type(MPI_type),intent(in) :: mpi_enreg
 type(efield_type),intent(in) :: dtefield
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer,optional :: electronpositron
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type), intent(in) :: wvl
 type(wvl_denspot_type), intent(inout) :: wvl_den
 type(fock_type),pointer, intent(inout) :: fock
!arrays
 integer,intent(in) :: atindx1(dtset%natom),indsym(4,dtset%nsym,dtset%natom)
 integer,intent(in) :: nattyp(ntypat),ngfft(18),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: grewtn(3,dtset%natom),grvdw(3,ngrvdw),grnl(3*dtset%natom)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,dtset%nspden),rprimd(3,3)
 real(dp),intent(in) :: vxc(nfft,dtset%nspden)
 real(dp),intent(inout) :: fcart(3,dtset%natom),forold(3,dtset%natom)
 real(dp),intent(inout) :: vresid(nfft,dtset%nspden),xred(3,dtset%natom)
 real(dp),intent(out) :: favg(3),fred(3,dtset%natom),gresid(3,dtset%natom)
 real(dp),intent(out) :: grhf(3,dtset%natom) !vz_i
 real(dp),intent(inout) :: grxc(3,dtset%natom) !vz_i
 real(dp),intent(out) :: synlgr(3,dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: fdir,iatom,idir,indx,ipositron,itypat,mu
 integer :: optatm,optdyfr,opteltfr,optgr,option,optn,optn2,optstr,optv
 real(dp) :: eei_dum,ucvol,vol_element
 logical :: calc_epaw3_forces, efield_flag
!arrays
 integer :: qprtrb_dum(3)
 real(dp) :: dummy6(6),ep3(3),fioncart(3),gmet(3,3),gprimd(3,3) 
 real(dp) :: rmet(3,3),strn_dummy6(6),strv_dummy6(6),tsec(2),vprtrb_dum(2)
 real(dp),allocatable :: atmrho_dum(:),atmvloc_dum(:),dyfrlo_dum(:,:,:)
 real(dp),allocatable :: dyfrn_dum(:,:,:),dyfrv_dum(:,:,:)
 real(dp),allocatable :: dyfrx2_dum(:,:,:),eltfrn_dum(:,:),gauss_dum(:,:)
 real(dp),allocatable :: epawf3red(:,:),fin(:,:),fionred(:,:),grl(:,:),grnl_tmp(:,:),grtn(:,:)
 real(dp),allocatable :: grtn_indx(:,:),v_dum(:),vxctotg(:,:)
 real(dp),allocatable :: xccc3d_dum(:)

! *************************************************************************

 call timab(69,1,tsec)

!Save input value of forces
 ABI_ALLOCATE(fin,(3,dtset%natom))
 fin(:,:)=fcart(:,:)

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!=======================================================================
!========= Local pseudopotential and core charge contributions =========
!=======================================================================

 ABI_ALLOCATE(grl,(3,dtset%natom))

!PAW or NC with nc_xccc_gspace: compute local psp and core charge contribs together
!in reciprocal space
!-----------------------------------------------------------------------
 if (psps%usepaw==1 .or. psps%nc_xccc_gspace==1) then
   call timab(550,1,tsec)
   optatm=0;optdyfr=0;optgr=1;optstr=0;optv=1;optn=n3xccc/nfft;optn2=1;opteltfr=0

   if (n3xccc>0) then
     ABI_ALLOCATE(v_dum,(nfft))
     ABI_ALLOCATE(vxctotg,(2,nfft))
     v_dum(:)=vxc(:,1);if (dtset%nspden>=2) v_dum(:)=0.5_dp*(v_dum(:)+vxc(:,2))
     call fourdp(1,vxctotg,v_dum,-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
     call zerosym(vxctotg,2,ngfft(1),ngfft(2),ngfft(3),&
&     comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     ABI_DEALLOCATE(v_dum)
   else 
     ! Allocate zero-sized array.
     ABI_ALLOCATE(vxctotg,(0,0))
   end if

! allocate (unused) dummy variables, otherwise some compilers complain
   ABI_ALLOCATE(gauss_dum,(0,0))
   ABI_ALLOCATE(atmrho_dum,(0))
   ABI_ALLOCATE(atmvloc_dum,(0))
   ABI_ALLOCATE(dyfrn_dum,(0,0,0))
   ABI_ALLOCATE(dyfrv_dum,(0,0,0))
   ABI_ALLOCATE(eltfrn_dum,(0,0))

   call atm2fft(atindx1,atmrho_dum,atmvloc_dum,dyfrn_dum,dyfrv_dum,&
&   eltfrn_dum,gauss_dum,gmet,gprimd,&
&   grxc,grl,gsqcut,mgfft,psps%mqgrid_vl,dtset%natom,nattyp,nfft,ngfft,ntypat,&
&   optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,psps,pawtab,ph1d,psps%qgrid_vl,qprtrb_dum,&
&   rhog,strn_dummy6,strv_dummy6,ucvol,psps%usepaw,vxctotg,vxctotg,vxctotg,vprtrb_dum,psps%vlspl,&
&   comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
&   paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)

   ABI_DEALLOCATE(gauss_dum)
   ABI_DEALLOCATE(atmrho_dum)
   ABI_DEALLOCATE(atmvloc_dum)
   ABI_DEALLOCATE(dyfrn_dum)
   ABI_DEALLOCATE(dyfrv_dum)
   ABI_DEALLOCATE(eltfrn_dum)
   ABI_DEALLOCATE(vxctotg)

   if (n3xccc==0) grxc=zero
   call timab(550,2,tsec)
 else

!  Norm-conserving: compute local psp contribution in reciprocal space
!  and core charge contribution in real space
!  -----------------------------------------------------------------------
   option=2
   ABI_ALLOCATE(dyfrlo_dum,(3,3,dtset%natom))
   ABI_ALLOCATE(grtn_indx,(3,dtset%natom))
   ABI_ALLOCATE(v_dum,(nfft))
   call mklocl(dtset,dyfrlo_dum,eei_dum,gmet,gprimd,grtn_indx,gsqcut,dummy6,mgfft,&
&   mpi_enreg,dtset%natom,nattyp,nfft,ngfft,dtset%nspden,ntypat,option,pawtab,ph1d,psps,&
&   qprtrb_dum,rhog,rhor,rprimd,ucvol,vprtrb_dum,v_dum,wvl,wvl_den,xred)

   do iatom=1,dtset%natom
!    Has to use the indexing array atindx1
     grl(1:3,atindx1(iatom))=grtn_indx(1:3,iatom)
   end do
   ABI_DEALLOCATE(dyfrlo_dum)
   ABI_DEALLOCATE(grtn_indx)
   ABI_DEALLOCATE(v_dum)
!  If gradients are computed in real space, we need to symetrise
!  the system before summing.
!  Rshaltaf: I changed the following line to include surfaces BC
   if (dtset%icoulomb == 1 .or. dtset%icoulomb == 2) then
     ABI_ALLOCATE(grnl_tmp,(3,dtset%natom))
     call sygrad(grnl_tmp,dtset%natom,grl,dtset%nsym,symrec,indsym)
     grl(:, :) = grnl_tmp(:, :)
     ABI_DEALLOCATE(grnl_tmp)
   end if

   if (n3xccc>0) then
     call timab(53,1,tsec)
     ABI_ALLOCATE(dyfrx2_dum,(3,3,dtset%natom))
     ABI_ALLOCATE(xccc3d_dum,(n3xccc))
     call mkcore(dummy6,dyfrx2_dum,grxc,mpi_enreg,dtset%natom,nfft,dtset%nspden,ntypat,ngfft(1),n1xccc,ngfft(2),&
&     ngfft(3),option,rprimd,dtset%typat,ucvol,vxc,psps%xcccrc,psps%xccc1d,xccc3d_dum,xred)
     ABI_DEALLOCATE(dyfrx2_dum)
     ABI_DEALLOCATE(xccc3d_dum)
     call timab(53,2,tsec)
   else
     grxc(:,:)=zero
   end if
 end if

!=======================================================================
!===================== Nonlocal contributions ==========================
!=======================================================================

!Only has to apply symmetries
 ABI_ALLOCATE(grnl_tmp,(3,dtset%natom))
 do iatom=1,dtset%natom
   indx=3*(iatom-1);grnl_tmp(1:3,atindx1(iatom))=grnl(indx+1:indx+3)
 end do
 if (dtset%usewvl == 0) then
   call sygrad(synlgr,dtset%natom,grnl_tmp,dtset%nsym,symrec,indsym)
 else
   synlgr = grnl_tmp
 end if
 ABI_DEALLOCATE(grnl_tmp)

!=======================================================================
!============ Density/potential residual contributions =================
!=======================================================================

 if (dtset%usewvl==0.and.abs(dtset%densfor_pred)>=1.and.abs(dtset%densfor_pred)<=3) then
   call fresid(dtset,gresid,mpi_enreg,nfft,ngfft,ntypat,1,&
&   pawtab,rhor,rprimd,ucvol,vresid,xred,xred,psps%znuclpsp)
 else if (dtset%usewvl==0.and.(abs(dtset%densfor_pred)==4.or.abs(dtset%densfor_pred)==6)) then
   call fresidrsp(atindx1,dtset,gmet,gprimd,gresid,gsqcut,mgfft,&
&   mpi_enreg,psps%mqgrid_vl,nattyp,nfft,ngfft,ntypat,psps,pawtab,ph1d,&
&   psps%qgrid_vl,ucvol,psps%usepaw,vresid,psps%zionpsp,psps%znuclpsp)
 else
   gresid(:,:)=zero
 end if

!=======================================================================
!======================= Other contributions ===========================
!=======================================================================

!Ewald energy contribution to forces as already been computed in "ewald"

!Potential residual contribution to forces as already been computed (forstr)

!Add Berry phase contributions (berryopt == 4/6/7/14/16/17)
!(compute the electric field force on the ion cores)
 efield_flag = (dtset%berryopt==4 .or. dtset%berryopt==6 .or. dtset%berryopt==7 .or. &
& dtset%berryopt==14 .or. dtset%berryopt==16 .or. dtset%berryopt==17)
 calc_epaw3_forces = (efield_flag .and. dtset%optforces /= 0 .and. psps%usepaw == 1)
 if ( efield_flag ) then 
   ABI_ALLOCATE(fionred,(3,dtset%natom))
   fionred(:,:)=zero
   do iatom=1,dtset%natom
     itypat=dtset%typat(iatom)
! force on ion due to electric field, cartesian representation
     fioncart(:)=psps%ziontypat(itypat)*dtset%efield(:)
! form fionred = rprimd^T * fioncart, note that forces transform 
! oppositely to coordinates, because they are derivative with respect to 
! coordinates
     call dgemv('T',3,3,one,rprimd,3,fioncart,1,zero,fionred(1:3,iatom),1)
!     do mu=1,3
!       fionred(mu,iatom)=rprimd(1,mu)*fioncart(1) &
!&       +rprimd(2,mu)*fioncart(2) &
!&       +rprimd(3,mu)*fioncart(3)
!     end do
   end do
 end if

!(compute additional F3-type force due to projectors for electric field with PAW)
 if ( efield_flag .and. calc_epaw3_forces ) then  
   ABI_ALLOCATE(epawf3red,(3,dtset%natom))
! dtefield%epawf3(iatom,idir,fdir) contains 
   epawf3red(:,:)=zero
   do iatom=1,dtset%natom
     do fdir = 1, 3
       do idir = 1, 3
! vol_element is volume/pt for integration of epawf3. volume is BZ volume
! so 1/ucvol, and number of kpts is nstr(idir)*nkstr(idir)
         vol_element=one/(ucvol*dtefield%nstr(idir)*dtefield%nkstr(idir))
         ep3(idir) = vol_element*dtefield%epawf3(iatom,idir,fdir)
       end do
       epawf3red(fdir,iatom) = -ucvol*dot_product(dtset%red_efieldbar(1:3),ep3(1:3))
     end do
   end do ! end loop over iatom
 end if

!This was incorrect coding. Bug found by Jiawang Hong
!if (dtset%berryopt==4) then
!allocate(fionred(3,dtset%natom));fionred(:,:)=zero
!iatom = 0
!do itypat=1,ntypat
!do iattyp=1,nattyp(itypat)
!iatom=iatom+1
!fioncart(:)=psps%ziontypat(itypat)*dtset%efield(:)
!do mu=1,3
!fionred(mu,iatom)=rprimd(1,mu)*fioncart(1) &
!&         +rprimd(2,mu)*fioncart(2) &
!&         +rprimd(3,mu)*fioncart(3)
!end do
!end do
!end do
!end if

!=======================================================================
!======= Assemble the various contributions to the forces ==============
!=======================================================================

!Collect grads of etot wrt reduced coordinates
!This gives non-symmetrized Hellman-Feynman reduced gradients
 ABI_ALLOCATE(grtn,(3,dtset%natom))
 grtn(:,:)=grl(:,:)+grewtn(:,:)+synlgr(:,:)+grxc(:,:)
 if (usefock==1 .and. associated(fock).and.fock%optfor) then
   grtn(:,:)=grtn(:,:)+fock%forces(:,:)
 end if
 if (ngrvdw==dtset%natom) grtn(:,:)=grtn(:,:)+grvdw(:,:)
! note that fionred is subtracted, because it really is a force and we need to
! turn it back into a gradient. The fred2fcart routine below includes the minus
! sign to convert gradients back to forces
 if ( efield_flag ) grtn(:,:)=grtn(:,:)-fionred(:,:)  
! epawf3red is added, because it actually is a gradient, not a force
 if ( efield_flag .and. calc_epaw3_forces ) grtn(:,:) = grtn(:,:) + epawf3red(:,:)

!Symmetrize explicitly for given space group and store in grhf :
 call sygrad(grhf,dtset%natom,grtn,dtset%nsym,symrec,indsym)

!If residual forces are too large, there must be a problem: cancel them !
 if (dtset%usewvl==0.and.abs(dtset%densfor_pred)>0.and.abs(dtset%densfor_pred)/=5) then
   do iatom=1,dtset%natom
     do mu=1,3
       if (abs(gresid(mu,iatom))>10000._dp*abs(grtn(mu,iatom))) gresid(mu,iatom)=zero
     end do
   end do
 end if

!Add residual potential correction
 grtn(:,:)=grtn(:,:)+gresid(:,:)

!Additional stuff for electron-positron
 ipositron=0
 if (present(electronpositron)) then
   if (associated(electronpositron)) then
     if (allocated(electronpositron%fred_ep)) ipositron=electronpositron_calctype(electronpositron)
   end if
 end if
 if (abs(ipositron)==1) then
   grtn(:,:)=grtn(:,:)-grxc(:,:)-grewtn(:,:)-gresid(:,:)-two*grl(:,:)
   grl(:,:)=-grl(:,:);grxc(:,:)=zero;gresid(:,:)=zero
   if (ngrvdw==dtset%natom) grtn(:,:)=grtn(:,:)-grvdw(:,:)
   if ( dtset%berryopt==4 .or. dtset%berryopt==6 .or. dtset%berryopt==7 .or. &
&   dtset%berryopt==14 .or. dtset%berryopt==16 .or. dtset%berryopt==17)  then     !!HONG
     grtn(:,:)=grtn(:,:)+fionred(:,:)
     fionred(:,:)=zero
   end if
 end if
 if (ipositron>0) grtn(:,:)=grtn(:,:)+electronpositron%fred_ep(:,:)

!Symmetrize all grads explicitly for given space group:
 if (dtset%usewvl == 0) then
   call sygrad(fred,dtset%natom,grtn,dtset%nsym,symrec,indsym)
 else
   fred = grtn
 end if

!Conversion to cartesian coordinates (bohr) AND
!Subtract off average force from each force component
!to avoid spurious drifting of atoms across cell.
! notice that fred2fcart multiplies fred by -1 to convert it 
! from a gradient (input) to a force (output)

 call fred2fcart(favg,(dtset%jellslab==0),fcart,fred,gprimd,dtset%natom)

!Compute maximal force and maximal difference
 maxfor=zero;diffor=zero
 do iatom=1,dtset%natom
   do mu=1,3
     if (dtset%iatfix(mu,iatom) /= 1) then
       maxfor=max(maxfor,abs(fcart(mu,iatom)))
       diffor=max(diffor,abs(fcart(mu,iatom)-fin(mu,iatom)))
     else if (dtset%ionmov==4 .or. dtset%ionmov==5) then
!      Make the force vanish on fixed atoms when ionmov=4 or 5
!      This is because fixing of atom cannot be imposed at the
!      level of a routine similar to brdmin or moldyn for these options.
       fcart(mu,iatom)=zero
     end if
   end do
 end do

!Apply any generalized constraints to the forces
 if (dtset%nconeq>0) call constrf(diffor,fcart,forold,fred,dtset%iatfix,dtset%ionmov,maxfor,&
& dtset%natom,dtset%nconeq,dtset%prtvol,rprimd,dtset%wtatcon,xred)

!=======================================================================
!Memory deallocations
 ABI_DEALLOCATE(grl)
 ABI_DEALLOCATE(grtn)
 ABI_DEALLOCATE(fin)
 if ( efield_flag )  then   
   ABI_DEALLOCATE(fionred)
   if ( calc_epaw3_forces ) then 
     ABI_DEALLOCATE(epawf3red)
   end if
 end if

 call timab(69,2,tsec)

end subroutine forces
!!***
