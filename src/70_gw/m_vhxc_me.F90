!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_vhxc_me
!! NAME
!! m_vhxc_me
!!
!! FUNCTION
!!  Evaluate the matrix elements of $v_H$ and $v_{xc}$ and $v_U$
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (MG)
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

module m_vhxc_me

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_abicore
 use m_errors
 use m_xcdata
 use libxc_functionals

 use m_pawang,      only : pawang_type
 use m_pawtab,      only : pawtab_type
 use m_paw_an,      only : paw_an_type
 use m_paw_ij,      only : paw_ij_type
 use m_pawfgrtab,   only : pawfgrtab_type
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_paw_denpot,  only : paw_mknewh0
 use m_hide_blas,   only : xdotc
 use m_wfd,         only : wfd_t
 use m_crystal,     only : crystal_t
 use m_melemts,     only : melements_init, melements_herm, melements_mpisum, melflags_t, melements_t
 use m_dtset,       only : dtset_copy, dtset_free
 use m_mpinfo,      only : destroy_mpi_enreg, initmpi_seq
 use m_kg,          only : mkkin
 use m_rhotoxc,     only : rhotoxc

 implicit none

 private
!!***

 public :: calc_vhxc_me
!!***

contains
!!***

!!****f* ABINIT/calc_vhxc_me
!! NAME
!!  calc_vhxc_me
!!
!! FUNCTION
!!  Evaluate the matrix elements of $v_H$ and $v_{xc}$ and $v_U$
!!  both in case of NC pseudopotentials and PAW (LDA+U, presently, is only available in PAW)
!!  The matrix elements of $v_{xc}$ are calculated with and without the core contribution.
!!  The later quantity is required in case of GW calculations.
!!
!! INPUTS
!!  Mflags
!!   that of the basis sphere--appropriate for charge density rho(G),Hartree potential, and pseudopotentials
!!  Dtset <type(dataset_type)>=all input variables in this dataset
!!  ngfftf(18)contain all needed information about 3D fine FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nfftf=number of points in the fine FFT mesh (for this processor)
!!  Pawtab(Dtset%ntypat*Dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  Paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  Pawang <type(pawang_type)>=paw angular mesh and related data
!!  Paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  Pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  Cryst<crystal_t>=unit cell and symmetries
!!     %natom=number of atoms in the unit cell
!!     %rprimd(3,3)=direct lattice vectors
!!     %ucvol=unit cell volume
!!     %ntypat= number of type of atoms
!!     %typat(natom)=type of each atom
!!  vhartr(nfftf)= Hartree potential in real space on the fine FFT mesh
!!  vxc(nfftf,nspden)= xc potential in real space on the fine FFT grid
!!  Wfd <type (wfd_t)>=Structure gathering information on the wavefunctions.
!!  rhor(nfftf,nspden)=density in real space (smooth part if PAW).
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!  kstab(2,Wfd%nkibz,Wfd%nsppol)=Table temporary used to be compatible with the old implementation.
!!
!! OUTPUT
!!  Mels
!!   %vxc   =matrix elements of $v_{xc}[nv+nc]$.
!!   %vxcval=matrix elements of $v_{xc}[nv]$.
!!   %vxcval_hybrid=matrix elements of $v_{xc}[nv]^{hybrid functional}$.
!!   %vhartr=matrix elements of $v_H$.
!!   %vu    =matrix elements of $v_U$.
!!
!! SIDE EFFECTS
!!  Paw_ij= In case of self-Consistency it is changed. It will contain the new H0
!!  Hamiltonian calculated using the QP densities. The valence contribution to XC
!!  is removed.
!!
!! NOTES
!!  All the quantities ($v_H$, $v_{xc}$ and $\psi$ are evaluated on the "fine" FFT mesh.
!!  In case of calculations with pseudopotials the usual mesh is defined by ecut.
!!  For PAW calculations the dense FFT grid defined bt ecutdg is used
!!  Besides, in case of PAW, the matrix elements of V_hartree do not contain the onsite
!!  contributions due to the coulombian potentials generate by ncore and tncore.
!!  These quantities, as well as the onsite kinetic terms, are stored in Paw_ij%dij0.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      destroy_mpi_enreg,get_auxc_ixc,init_distribfft_seq,initmpi_seq
!!      libxc_functionals_end,libxc_functionals_init
!!      libxc_functionals_set_hybridparams,melements_herm,melements_init
!!      melements_mpisum,mkkin,paw_mknewh0,pawcprj_alloc,pawcprj_free,rhotoxc
!!      wfd_change_ngfft,wfd_distribute_bbp,wfd_get_cprj,wfd_get_ur,wrtout
!!      xcdata_init
!!
!! SOURCE


subroutine calc_vhxc_me(Wfd,Mflags,Mels,Cryst,Dtset,nfftf,ngfftf,&
  vtrial,vhartr,vxc,Psps,Pawtab,Paw_an,Pawang,Pawfgrtab,Paw_ij,dijexc_core,&
  rhor,usexcnhat,nhat,nhatgr,nhatgrdim,kstab,&
  taur) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nhatgrdim,usexcnhat,nfftf
 type(Dataset_type),intent(in) :: Dtset
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfd_t),target,intent(inout) :: Wfd
 type(Pawang_type),intent(in) :: Pawang
 type(crystal_t),intent(in) :: Cryst
 type(melflags_t),intent(in) :: Mflags
 type(melements_t),intent(out) :: Mels
!arrays
 integer,intent(in) :: ngfftf(18)
 integer,intent(in) :: kstab(2,Wfd%nkibz,Wfd%nsppol)
 real(dp),intent(in) :: vhartr(nfftf),vxc(nfftf,Wfd%nspden),vtrial(nfftf,Wfd%nspden)
 real(dp),intent(in) :: rhor(nfftf,Wfd%nspden)
 real(dp),intent(in) :: nhat(nfftf,Wfd%nspden*Wfd%usepaw)
 real(dp),intent(in) :: nhatgr(nfftf,Wfd%nspden,3*nhatgrdim)
 real(dp),intent(in),optional :: taur(nfftf,Wfd%nspden*Dtset%usekden)
 !real(dp),intent(in) :: dijexc_core(cplex_dij*lmn2_size_max,ndij,Cryst%ntypat)
 real(dp),intent(in) :: dijexc_core(:,:,:)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(Paw_an_type),intent(in) :: Paw_an(Cryst%natom)
 type(Paw_ij_type),intent(inout) :: Paw_ij(Cryst%natom)
 type(Pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom)

!Local variables-------------------------------
!scalars
 integer :: auxc_ixc,iat,ikc,ik_ibz,ib,jb,is,b_start,b_stop,istwf_k
 integer :: itypat,lmn_size,j0lmn,jlmn,ilmn,klmn,klmn1,lmn2_size_max
 integer :: isppol,cplex_dij,npw_k
 integer :: nspinor,nsppol,nspden,nk_calc
 integer :: rank
 integer :: iab,isp1,isp2,ixc_sigma,nsploop,nkxc,option,n3xccc_,nk3xc,my_nbbp,my_nmels
 real(dp) :: nfftfm1,fact,DijH,enxc_val,enxc_hybrid_val,vxcval_avg,vxcval_hybrid_avg,h0dij,vxc1,vxc1_val,re_p,im_p,dijsigcx
 complex(dpc) :: cdot
 logical :: ltest,nmxc
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
 type(xcdata_type) :: xcdata,xcdata_hybrid
!arrays
 integer,parameter :: spinor_idxs(2,4)=RESHAPE([1,1,2,2,1,2,2,1], [2,4])
 integer :: got(Wfd%nproc)
 integer,allocatable :: kcalc2ibz(:),dimlmn(:),bbp_ks_distrb(:,:,:,:)
 integer,ABI_CONTIGUOUS pointer :: kg_k(:,:)
 real(dp) :: tmp_xc(2,Wfd%nspinor**2),tmp_xcval(2,Wfd%nspinor**2)
 real(dp) :: tmp_H(2,Wfd%nspinor**2),tmp_U(2,Wfd%nspinor**2)
 real(dp) :: tmp_h0ij(2,Wfd%nspinor**2),tmp_sigcx(2,Wfd%nspinor**2)
 real(dp) :: dijU(2),strsxc(6),kpt(3),vxc1ab(2),vxc1ab_val(2)
 real(dp),allocatable :: kxc_(:,:),xccc3d_(:),vxc_val(:,:),vxc_val_hybrid(:,:)
 real(dp),allocatable :: kinpw(:),veffh0(:,:)
 complex(dpc) :: tmp(3)
 complex(gwpc),ABI_CONTIGUOUS pointer :: ur1_up(:),ur1_dwn(:),ur2_up(:),ur2_dwn(:),cg1(:),cg2(:)
 complex(gwpc),target,allocatable :: ur1(:),ur2(:)
 complex(dpc),allocatable :: vxcab(:),vxcab_val(:),vxcab_val_hybrid(:),u1cjg_u2dpc(:),kinwf2(:),veffh0_ab(:)
 logical,allocatable :: bbp_mask(:,:)
 type(pawcprj_type),allocatable ::  Cprj_b1ks(:,:),Cprj_b2ks(:,:)
 type(libxc_functional_type) :: xc_funcs_hybrid(2)

! *************************************************************************

 DBG_ENTER("COLL")

 ABI_MALLOC(bbp_mask,(Wfd%mband,Wfd%mband))

 ! Usually FFT meshes for wavefunctions and potentials are not equal. Two approaches are possible:
 ! Either we Fourier interpolate potentials on the coarse WF mesh or we FFT the wfs on the dense mesh.
 ! The later approach is used, more CPU demanding but more accurate.
 if ( ANY(ngfftf(1:3) /= Wfd%ngfft(1:3)) ) call wfd%change_ngfft(Cryst,Psps,ngfftf)

 ! Fake MPI_type for sequential part
 rank = Wfd%my_rank
 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',ngfftf(2),ngfftf(3),'all')

 nspinor=Wfd%nspinor; nsppol =Wfd%nsppol; nspden =Wfd%nspden
 if (nspinor == 2) MSG_WARNING("Remember to ADD SO")

 ! TODO not used for the time being but it should be a standard input of the routine.
 !  bbks_mask(Wfd%mband,Wfd%mband,Wfd%nkibz,Wfd%nsppol)=Logical mask used to select
 !  the matrix elements to be calculated.
 ABI_MALLOC(kcalc2ibz,(Wfd%nkibz))
 kcalc2ibz=0

 ! Index in the IBZ of the GW k-points.
 ! Only these points will be considered.
 nk_calc=0
 do ik_ibz=1,Wfd%nkibz
   if ( ALL(kstab(1,ik_ibz,:)/=0) .and. ALL(kstab(2,ik_ibz,:)/=0) ) then
     nk_calc=nk_calc+1; kcalc2ibz(nk_calc) = ik_ibz
   end if
 end do

 call melements_init(Mels,Mflags,nsppol,nspden,Wfd%nspinor,Wfd%nkibz,Wfd%kibz,kstab)

 if (Mflags%has_lexexch==1) then
   MSG_ERROR("Local EXX not coded!")
 end if

 ! Evaluate $v_\xc$ using only the valence charge.
 call wrtout(std_out," calc_vhxc_braket : calculating v_xc[n_val] (excluding non-linear core corrections)")

 do isppol=1,nsppol
   write(msg,'(a,i2,a,e16.6)')' For spin ',isppol,' Min density rhor = ',MINVAL(rhor(:,isppol))
   call wrtout(std_out,msg,'COLL')
   if (Wfd%usepaw==1) then
     write(msg,'(a,i2,a,e16.6)')' For spin ',isppol,' Min density nhat = ',MINVAL(nhat(:,isppol))
     call wrtout(std_out,msg,'COLL')
     write(msg,'(a,i2,a,e16.6)')' For spin ',isppol,' Min density trho-nhat = ',MINVAL(rhor(:,isppol)-nhat(:,isppol))
     call wrtout(std_out,msg,'COLL')
     write(msg,'(a,i2)')' using usexcnhat = ',usexcnhat
     call wrtout(std_out,msg,'COLL')
   end if
 end do

 option = 0 ! Only exc, vxc, strsxc
 nkxc   = 0 ! No computation of XC kernel
 n3xccc_= 0 ! No core
 nk3xc  = 0 ! k3xc not needed
 nmxc=(Dtset%usepaw==1.and.mod(abs(Dtset%usepawu),10)==4)

 ABI_MALLOC(xccc3d_,(n3xccc_))
 ABI_MALLOC(kxc_,(nfftf,nkxc))
 ABI_MALLOC(vxc_val,(nfftf,nspden))

 call xcdata_init(xcdata,dtset=Dtset)

 call rhotoxc(enxc_val,kxc_,MPI_enreg_seq,nfftf,ngfftf,&
& nhat,Wfd%usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,nmxc,n3xccc_,option,rhor,Cryst%rprimd,&
& strsxc,usexcnhat,vxc_val,vxcval_avg,xccc3d_,xcdata,taur=taur)

 ! FABIEN's development
 ! Hybrid functional treatment
 if(Mflags%has_vxcval_hybrid==1) then

   call wrtout(std_out,' Hybrid functional xc potential is being set')
   ixc_sigma=Dtset%ixc_sigma
   call get_auxc_ixc(auxc_ixc,ixc_sigma)
   call xcdata_init(xcdata_hybrid,dtset=Dtset,auxc_ixc=auxc_ixc,ixc=ixc_sigma)

   if(ixc_sigma<0)then
     if(libxc_functionals_check()) then
       call libxc_functionals_init(ixc_sigma,Dtset%nspden,xc_functionals=xc_funcs_hybrid)
!      Do not forget, negative values of hyb_mixing(_sr),hyb_range_* means that they have been user-defined.
       if (dtset%ixc==-402.or.dtset%ixc==-406.or.dtset%ixc==-427.or.dtset%ixc==-428 .or. dtset%ixc==-456 .or. &
&        min(Dtset%hyb_mixing,Dtset%hyb_mixing_sr,Dtset%hyb_range_dft,Dtset%hyb_range_fock)<-tol8)then
         call libxc_functionals_set_hybridparams(hyb_range=abs(Dtset%hyb_range_dft),&
&          hyb_mixing=abs(Dtset%hyb_mixing),hyb_mixing_sr=abs(Dtset%hyb_mixing_sr),xc_functionals=xc_funcs_hybrid)
       endif
     else
       call wrtout(std_out, 'LIBXC is not present: hybrid functionals are not available')
     end if
   end if

   write(msg, '(a, f4.2)') ' Fock fraction = ', max(abs(Dtset%hyb_mixing),abs(Dtset%hyb_mixing_sr))
   call wrtout(std_out,msg,'COLL')
   write(msg, '(a, f5.2, a)') ' Fock inverse screening length = ',abs(Dtset%hyb_range_dft), ' (bohr^-1)'
   call wrtout(std_out,msg,'COLL')

   ABI_MALLOC(vxc_val_hybrid,(nfftf,nspden))

   if(ixc_sigma<0)then
     call rhotoxc(enxc_hybrid_val,kxc_,MPI_enreg_seq,nfftf,ngfftf,&
&     nhat,Wfd%usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,nmxc,n3xccc_,option,rhor,Cryst%rprimd,&
&     strsxc,usexcnhat,vxc_val_hybrid,vxcval_hybrid_avg,xccc3d_,xcdata_hybrid,xc_funcs=xc_funcs_hybrid)
     call libxc_functionals_end(xc_functionals=xc_funcs_hybrid)
   else
     call rhotoxc(enxc_hybrid_val,kxc_,MPI_enreg_seq,nfftf,ngfftf,&
&     nhat,Wfd%usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,nmxc,n3xccc_,option,rhor,Cryst%rprimd,&
&     strsxc,usexcnhat,vxc_val_hybrid,vxcval_hybrid_avg,xccc3d_,xcdata_hybrid)
   end if

 endif

 ABI_FREE(xccc3d_)
 ABI_FREE(kxc_)

 write(msg,'(a,f8.4,2a,f8.4,a)')' E_xc[n_val]  = ',enxc_val,  ' [Ha]. ','<V_xc[n_val]> = ',vxcval_avg,' [Ha]. '
 call wrtout(std_out,msg,'COLL')

 ! If PAW and qp-SCGW then update Paw_ij and calculate the matrix elements ===
 ! We cannot simply rely on gwcalctyp because I need KS vxc in sigma.
 if (Wfd%usepaw==1.and.Mflags%has_hbare==1) then
   ABI_CHECK(Mflags%only_diago==0,"Wrong only_diago")

   call paw_mknewh0(Cryst%natom,nsppol,nspden,nfftf,Dtset%pawspnorb,Dtset%pawprtvol,Cryst,&
&    Pawtab,Paw_an,Paw_ij,Pawang,Pawfgrtab,vxc,vxc_val,vtrial)

   ! Effective potential of the bare Hamiltonian: valence term is subtracted.
   ABI_MALLOC(veffh0,(nfftf,nspden))
   veffh0=vtrial-vxc_val
   !veffh0=vtrial !this is to retrieve the KS Hamiltonian
 end if

 ! Setup of the hermitian operator vxcab ===
 ! if nspden==4 vxc contains (v^11, v^22, Re[V^12], Im[V^12].
 ! Cannot use directly Re and Im since we also need off-diagonal elements.
 if (wfd%nspden == 4) then
   ABI_MALLOC(vxcab, (nfftf))
   ABI_MALLOC(vxcab_val, (nfftf))
   vxcab    (:) = DCMPLX(vxc    (:,3), vxc    (:,4))
   vxcab_val(:) = DCMPLX(vxc_val(:,3), vxc_val(:,4))
   if (Mflags%has_vxcval_hybrid==1) then
     ABI_MALLOC(vxcab_val_hybrid,(nfftf))
     vxcab_val_hybrid(:)=DCMPLX(vxc_val_hybrid(:,3),vxc_val_hybrid(:,4))
   end if
   if (Mflags%has_hbare==1) then
     ABI_MALLOC(veffh0_ab,(nfftf))
     veffh0_ab(:)=DCMPLX(veffh0(:,3),veffh0(:,4))
   end if
 end if

 ABI_MALLOC(ur1, (nfftf * nspinor))
 ABI_MALLOC(ur2, (nfftf * nspinor))
 ABI_MALLOC(u1cjg_u2dpc, (nfftf * nspinor))

 ! Create distribution table for tasks ===
 ! This section is parallelized inside wfd%comm
 ! as all processors are calling the routine with all GW wavefunctions
 ! TODO the table can be calculated at each (k,s) to save some memory.
 got=0; my_nmels=0
 ABI_MALLOC(bbp_ks_distrb,(Wfd%mband,Wfd%mband,nk_calc,nsppol))
 do is=1,nsppol
   do ikc=1,nk_calc
     ik_ibz=kcalc2ibz(ikc)
     bbp_mask=.FALSE.
     b_start=kstab(1,ik_ibz,is)
     b_stop =kstab(2,ik_ibz,is)
     if (Mflags%only_diago==1) then
       !do jb=b1,b2
       do jb=b_start,b_stop
         bbp_mask(jb,jb)=.TRUE.
       end do
     else
       bbp_mask(b_start:b_stop,b_start:b_stop)=.TRUE.
     end if

     call wfd%distribute_bbp(ik_ibz,is,"Upper",my_nbbp,bbp_ks_distrb(:,:,ikc,is),got,bbp_mask)
     my_nmels = my_nmels + my_nbbp
   end do
 end do
 ABI_FREE(bbp_mask)

 write(msg,'(a,i0,a)')" Will calculate ",my_nmels," <b,k,s|O|b',k,s> matrix elements in calc_vhxc_me."
 call wrtout(std_out,msg,'PERS')

 ! =====================================
 ! ==== Loop over required k-points ====
 ! =====================================
 nfftfm1=one/nfftf

 do is=1,nsppol
   if (ALL(bbp_ks_distrb(:,:,:,is)/=rank)) CYCLE
   do ikc=1,nk_calc
     if (ALL(bbp_ks_distrb(:,:,ikc,is)/=rank)) CYCLE

     ik_ibz = kcalc2ibz(ikc)
     b_start = kstab(1,ik_ibz,is)
     b_stop  = kstab(2,ik_ibz,is)
     npw_k = Wfd%Kdata(ik_ibz)%npw
     kpt = Wfd%kibz(:,ik_ibz)
     kg_k => Wfd%kdata(ik_ibz)%kg_k
     istwf_k = wfd%istwfk(ik_ibz)

     ! Calculate |k+G|^2 needed by hbareme
     !FIXME Here I have a problem if I use ecutwfn there is a bug somewhere in setshell or invars2m!
     ! ecutwfn is slightly smaller than the max kinetic energy in gvec. The 0.1 pad should partially solve the problem
     if (Mflags%has_hbare==1) then
       ABI_MALLOC(kinpw,(npw_k))
       ABI_MALLOC(kinwf2,(npw_k*nspinor))
       call mkkin(Dtset%ecutwfn+0.1_dp,Dtset%ecutsm,Dtset%effmass_free,Cryst%gmet,kg_k,kinpw,kpt,npw_k,0,0)
       where (kinpw>HUGE(zero)*1.d-11)
         kinpw=zero
       end where
     end if

     !do jb=b1,b2
     do jb=b_start,b_stop
       if (ALL(bbp_ks_distrb(:,jb,ikc,is)/=rank)) CYCLE

       if (Mflags%has_hbare==1) then
         cg2 => Wfd%Wave(jb,ik_ibz,is)%ug  ! Wfd contains 1:nkptgw wave functions
         kinwf2(1:npw_k)=cg2(1:npw_k)*kinpw(:)
         if (nspinor==2) kinwf2(npw_k+1:)=cg2(npw_k+1:)*kinpw(:)
       end if

       call wfd%get_ur(jb,ik_ibz,is,ur2)

       !do ib=b1,jb ! Upper triangle
       do ib=b_start,jb
         if (bbp_ks_distrb(ib,jb,ikc,is)/=rank) CYCLE
         ! Off-diagonal elements only for QPSCGW.
         if (Mflags%only_diago==1.and.ib/=jb) CYCLE

         call wfd%get_ur(ib,ik_ibz,is,ur1)
         u1cjg_u2dpc(:) = CONJG(ur1) *ur2

         if (Mflags%has_vxc == 1) then
           Mels%vxc(ib, jb, ik_ibz, is) = sum(u1cjg_u2dpc(1:nfftf) * vxc(1:nfftf, is)) * nfftfm1
           if (wfd%nspinor == 2 .and. wfd%nspden == 1) &
             Mels%vxc(ib, jb, ik_ibz, 2) = sum(u1cjg_u2dpc(nfftf+1:) * vxc(1:nfftf, is)) * nfftfm1
         end if
         if (Mflags%has_vxcval == 1) then
           Mels%vxcval(ib, jb, ik_ibz, is) = SUM(u1cjg_u2dpc(1:nfftf) * vxc_val(1:nfftf, is)) * nfftfm1
           if (wfd%nspinor == 2 .and. wfd%nspden == 1) &
             Mels%vxcval(ib, jb, ik_ibz, 2) = SUM(u1cjg_u2dpc(nfftf+1:) * vxc_val(1:nfftf, is)) * nfftfm1
         end if
         if (Mflags%has_vxcval_hybrid == 1) then
           Mels%vxcval_hybrid(ib, jb, ik_ibz, is) = SUM(u1cjg_u2dpc(1:nfftf) * vxc_val_hybrid(1:nfftf, is)) * nfftfm1
           if (wfd%nspinor == 2 .and. wfd%nspden == 1) &
             Mels%vxcval_hybrid(ib, jb, ik_ibz, 2) = SUM(u1cjg_u2dpc(nfftf+1) * vxc_val_hybrid(1:nfftf, is)) * nfftfm1
         end if
         if (Mflags%has_vhartree==1) then
           Mels%vhartree(ib, jb, ik_ibz, is) = SUM(u1cjg_u2dpc(1:nfftf) * vhartr(1:nfftf)) * nfftfm1
           if (wfd%nspinor == 2 .and. wfd%nspden == 1) &
             Mels%vhartree(ib, jb, ik_ibz, 2) = SUM(u1cjg_u2dpc(nfftf+1:) * vhartr(1:nfftf)) * nfftfm1
         end if

         if (Mflags%has_hbare==1) then
           cg1 => Wfd%Wave(ib, ik_ibz, is)%ug(1:npw_k)
           cdot = DOT_PRODUCT(cg1, kinwf2(1:npw_k))
           !if (istwf_k /= 1) then
           !  cdot = two * cdot; if (istwf_k == 2) cdot = cdot - GWPC_CONJG(cg1(1)) * kinwf2(1)
           !end if
           Mels%hbare(ib, jb, ik_ibz, is) = cdot + SUM(u1cjg_u2dpc(1:nfftf) * veffh0(1:nfftf, is)) * nfftfm1
           if (wfd%nspinor == 2 .and. wfd%nspden == 1) then
             cg1 => Wfd%Wave(ib, ik_ibz, is)%ug(npw_k+1:)
             Mels%hbare(ib, jb, ik_ibz, 2) = &
               DOT_PRODUCT(cg1, kinwf2(npw_k+1:)) + SUM(u1cjg_u2dpc(nfftf+1:) * veffh0(1:nfftf, is)) * nfftfm1
           end if
         end if

         if (nspinor == 2 .and. wfd%nspden == 4) then
           ! Here I can skip 21 if ib==jb
           ur1_up  => ur1(1:nfftf)
           ur1_dwn => ur1(nfftf+1:2*nfftf)
           ur2_up  => ur2(1:nfftf)
           ur2_dwn => ur2(nfftf+1:2*nfftf)

           if (Mflags%has_hbare==1) then
             cg1 => Wfd%Wave(ib,ik_ibz,is)%ug(npw_k+1:)
             tmp(1)=SUM(CONJG(ur1_dwn)*veffh0(:,2)*ur2_dwn)*nfftfm1 + DOT_PRODUCT(cg1,kinwf2(npw_k+1:))
             tmp(2)=SUM(CONJG(ur1_dwn)*      veffh0_ab(:) *ur2_dwn)*nfftfm1
             tmp(3)=SUM(CONJG(ur1_dwn)*CONJG(veffh0_ab(:))*ur2_dwn)*nfftfm1
             Mels%hbare(ib,jb,ik_ibz,2:4)=tmp(:)
           end if
           if (Mflags%has_vxc==1) then
             tmp(1) = SUM(CONJG(ur1_dwn)*      vxc(:,2) *ur2_dwn)*nfftfm1
             tmp(2) = SUM(CONJG(ur1_up )*      vxcab(:) *ur2_dwn)*nfftfm1
             tmp(3) = SUM(CONJG(ur1_dwn)*CONJG(vxcab(:))*ur2_up )*nfftfm1
             Mels%vxc(ib,jb,ik_ibz,2:4)=tmp(:)
           end if
           if (Mflags%has_vxcval==1) then
             tmp(1) = SUM(CONJG(ur1_dwn)*      vxc_val(:,2) *ur2_dwn)*nfftfm1
             tmp(2) = SUM(CONJG(ur1_up )*      vxcab_val(:) *ur2_dwn)*nfftfm1
             tmp(3) = SUM(CONJG(ur1_dwn)*CONJG(vxcab_val(:))*ur2_up )*nfftfm1
             Mels%vxcval(ib,jb,ik_ibz,2:4)=tmp(:)
           end if
           if (Mflags%has_vxcval_hybrid==1) then
             tmp(1) = SUM(CONJG(ur1_dwn)*      vxc_val_hybrid(:,2) *ur2_dwn)*nfftfm1
             tmp(2) = SUM(CONJG(ur1_up )*      vxcab_val_hybrid(:) *ur2_dwn)*nfftfm1
             tmp(3) = SUM(CONJG(ur1_dwn)*CONJG(vxcab_val_hybrid(:))*ur2_up )*nfftfm1
             Mels%vxcval_hybrid(ib,jb,ik_ibz,2:4)=tmp(:)
           end if
           if (Mflags%has_vhartree==1) then
             tmp(1) = SUM(CONJG(ur1_dwn)*vhartr(:)*ur2_dwn)*nfftfm1
             Mels%vhartree(ib,jb,ik_ibz,2  )=tmp(1)
             Mels%vhartree(ib,jb,ik_ibz,3:4)=czero
           end if
         end if !nspinor==2

       end do !ib
     end do !jb

     if (Mflags%has_hbare==1) then
       ABI_FREE(kinpw)
       ABI_FREE(kinwf2)
     end if

   end do !ikc
 end do !is

 ABI_FREE(ur1)
 ABI_FREE(ur2)
 ABI_FREE(vxc_val)
 ABI_FREE(u1cjg_u2dpc)
 if(Mflags%has_vxcval_hybrid==1) then
   ABI_FREE(vxc_val_hybrid)
 end if
 if (wfd%nspden == 4)  then
   ABI_FREE(vxcab)
   ABI_FREE(vxcab_val)
   if(Mflags%has_vxcval_hybrid==1) then
     ABI_FREE(vxcab_val_hybrid)
   end if
 end if

 if (Mflags%has_hbare==1) then
   ABI_FREE(veffh0)
   if (nspinor==2)  then
     ABI_FREE(veffh0_ab)
   end if
 end if

 ! ====================================
 ! ===== Additional terms for PAW =====
 ! ====================================
 if (Wfd%usepaw==1) then
   ! Tests if needed pointers in Paw_ij are allocated.
   ltest=(allocated(Paw_ij(1)%dijxc).and.allocated(Paw_ij(1)%dijxc_hat).and.allocated(Paw_ij(1)%dijxc_val))
   ABI_CHECK(ltest,"dijxc, dijxc_hat or dijxc_val not allocated")
   ABI_CHECK(nspinor == 1, "PAW with nspinor not tested")

   ! For LDA+U
   do iat=1,Cryst%natom
     itypat=Cryst%typat(iat)
     if (Pawtab(itypat)%usepawu/=0) then
       ltest=(allocated(Paw_ij(iat)%dijU))
       ABI_CHECK(ltest,"LDA+U but dijU not allocated")
     end if
   end do

   if (Dtset%pawspnorb>0) then
     ltest=(allocated(Paw_ij(1)%dijso))
     ABI_CHECK(ltest,"dijso not allocated")
   end if

   lmn2_size_max=MAXVAL(Pawtab(:)%lmn2_size)

   if (Mflags%has_sxcore==1) then
     if (     SIZE(dijexc_core,DIM=1) /= lmn2_size_max  &
&        .or. SIZE(dijexc_core,DIM=2) /= 1              &
&        .or. SIZE(dijexc_core,DIM=3) /= Cryst%ntypat ) then
       MSG_BUG("Wrong sizes in dijexc_core")
     end if
   end if

   nsploop=nspinor**2

   ! ====================================
   ! === Assemble PAW matrix elements ===
   ! ====================================
   ABI_MALLOC(dimlmn,(Cryst%natom))
   do iat=1,Cryst%natom
     dimlmn(iat)=Pawtab(Cryst%typat(iat))%lmn_size
   end do

   ABI_DATATYPE_ALLOCATE(Cprj_b1ks,(Cryst%natom,nspinor))
   ABI_DATATYPE_ALLOCATE(Cprj_b2ks,(Cryst%natom,nspinor))
   call pawcprj_alloc(Cprj_b1ks,0,dimlmn)
   call pawcprj_alloc(Cprj_b2ks,0,dimlmn)

   do is=1,nsppol
     if (ALL(bbp_ks_distrb(:,:,:,is)/=rank)) CYCLE

     ! Loop over required k-points
     do ikc=1,nk_calc
       if (ALL(bbp_ks_distrb(:,:,ikc,is)/=rank)) CYCLE
       ik_ibz=kcalc2ibz(ikc)
       b_start=kstab(1,ik_ibz,is)
       b_stop =kstab(2,ik_ibz,is)

       !do jb=b1,b2
       do jb=b_start,b_stop
         if (ALL(bbp_ks_distrb(:,jb,ikc,is)/=rank)) CYCLE

         ! Load projected wavefunctions for this k-point, spin and band ===
         ! Cprj are unsorted, full correspondence with xred. See ctocprj.F90!!
         call wfd%get_cprj(jb,ik_ibz,is,Cryst,Cprj_b2ks,sorted=.FALSE.)

         !do ib=b1,jb ! Upper triangle
         do ib=b_start,jb
           if (bbp_ks_distrb(ib,jb,ikc,is)/=rank) CYCLE

           ! * Off-diagonal elements only for QPSCGW.
           if (Mflags%only_diago==1.and.ib/=jb) CYCLE

           call wfd%get_cprj(ib,ik_ibz,is,Cryst,Cprj_b1ks,sorted=.FALSE.)
           !
           ! === Get onsite matrix elements summing over atoms and channels ===
           ! * Spin is external and fixed (1,2) if collinear.
           ! * if noncollinear loop internally over the four components ab.
           tmp_xc = zero; tmp_xcval = zero; tmp_H = zero; tmp_U = zero; tmp_h0ij = zero; tmp_sigcx = zero

           do iat=1,Cryst%natom
             itypat   =Cryst%typat(iat)
             lmn_size =Pawtab(itypat)%lmn_size
             cplex_dij=Paw_ij(iat)%cplex_dij
             klmn1=1

             do jlmn=1,lmn_size
               j0lmn=jlmn*(jlmn-1)/2
               do ilmn=1,jlmn
                 klmn=j0lmn+ilmn
                 ! TODO Be careful, here I assume that the onsite terms ij are symmetric
                 ! should check the spin-orbit case!
                 fact=one; if (ilmn==jlmn) fact=half

                 ! Loop over four components if nspinor==2 ===
                 ! If collinear nsploop==1
                 do iab=1,nsploop
                   isp1=spinor_idxs(1,iab); isp2=spinor_idxs(2,iab)

                   re_p=  Cprj_b1ks(iat,isp1)%cp(1,ilmn) * Cprj_b2ks(iat,isp2)%cp(1,jlmn) &
&                        +Cprj_b1ks(iat,isp1)%cp(2,ilmn) * Cprj_b2ks(iat,isp2)%cp(2,jlmn) &
&                        +Cprj_b1ks(iat,isp1)%cp(1,jlmn) * Cprj_b2ks(iat,isp2)%cp(1,ilmn) &
&                        +Cprj_b1ks(iat,isp1)%cp(2,jlmn) * Cprj_b2ks(iat,isp2)%cp(2,ilmn)

                   im_p=  Cprj_b1ks(iat,isp1)%cp(1,ilmn) * Cprj_b2ks(iat,isp2)%cp(2,jlmn) &
&                        -Cprj_b1ks(iat,isp1)%cp(2,ilmn) * Cprj_b2ks(iat,isp2)%cp(1,jlmn) &
&                        +Cprj_b1ks(iat,isp1)%cp(1,jlmn) * Cprj_b2ks(iat,isp2)%cp(2,ilmn) &
&                        -Cprj_b1ks(iat,isp1)%cp(2,jlmn) * Cprj_b2ks(iat,isp2)%cp(1,ilmn)

                   ! ==================================================
                   ! === Load onsite matrix elements and accumulate ===
                   ! ==================================================
                   if (nspinor==1) then

                     if (Mflags%has_hbare==1) then ! * Get new dij of h0 and accumulate.
                       h0dij=Paw_ij(iat)%dij(klmn,is)
                       tmp_h0ij(1,iab)=tmp_h0ij(1,iab) + h0dij*re_p*fact
                       tmp_h0ij(2,iab)=tmp_h0ij(2,iab) + h0dij*im_p*fact
                     end if

                     if (Mflags%has_sxcore==1) then ! * Fock operator generated by core electrons.
                       dijsigcx = dijexc_core(klmn,1,itypat)
                       tmp_sigcx(1,iab)=tmp_sigcx(1,iab) + dijsigcx*re_p*fact
                       tmp_sigcx(2,iab)=tmp_sigcx(2,iab) + dijsigcx*im_p*fact
                     end if

                     if (Mflags%has_vxc==1) then ! * Accumulate vxc[n1+nc] + vxc[n1+tn+nc].
                       vxc1 = Paw_ij(iat)%dijxc(klmn,is)+Paw_ij(iat)%dijxc_hat(klmn,is)
                       tmp_xc(1,iab)=tmp_xc(1,iab) + vxc1*re_p*fact
                       tmp_xc(2,iab)=tmp_xc(2,iab) + vxc1*im_p*fact
                     end if

                     if (Mflags%has_vxcval==1) then ! * Accumulate valence-only XC.
                       vxc1_val=Paw_ij(iat)%dijxc_val(klmn,is)
                       tmp_xcval(1,1)=tmp_xcval(1,1) + vxc1_val*re_p*fact
                       tmp_xcval(2,1)=tmp_xcval(2,1) + vxc1_val*im_p*fact
                     end if

                     if (Mflags%has_vhartree==1) then ! * Accumulate Hartree term of the PAW Hamiltonian.
                       DijH=Paw_ij(iat)%dijhartree(klmn)
                       tmp_H(1,1)=tmp_H(1,1) + DijH*re_p*fact
                       tmp_H(2,1)=tmp_H(2,1) + DijH*im_p*fact
                     end if

                     ! * Accumulate U term of the PAW Hamiltonian (only onsite AE contribution)
                     if (Mflags%has_vu==1) then
                       if (Pawtab(itypat)%usepawu/=0) then
                         dijU(1)=Paw_ij(iat)%dijU(klmn,is)
                         tmp_U(1,1)=tmp_U(1,1) + dijU(1)*re_p*fact
                         tmp_U(2,1)=tmp_U(2,1) + dijU(1)*im_p*fact
                       end if
                     end if

                   else ! Spinorial case ===

                     ! FIXME H0 + spinor not implemented
                     if (Mflags%has_hbare==1.or.Mflags%has_sxcore==1) then
                       MSG_ERROR("not implemented")
                     end if

                     if (Mflags%has_vxc==1) then ! * Accumulate vxc[n1+nc] + vxc[n1+tn+nc].
                       vxc1ab(1) = Paw_ij(iat)%dijxc(klmn1,  iab)+Paw_ij(iat)%dijxc_hat(klmn1,  iab)
                       vxc1ab(2) = Paw_ij(iat)%dijxc(klmn1+1,iab)+Paw_ij(iat)%dijxc_hat(klmn1+1,iab)
                       tmp_xc(1,iab) = tmp_xc(1,iab) + (vxc1ab(1)*re_p - vxc1ab(2)*im_p)*fact
                       tmp_xc(2,iab) = tmp_xc(2,iab) + (vxc1ab(2)*re_p + vxc1ab(1)*im_p)*fact
                     end if

                     if (Mflags%has_vxcval==1) then ! * Accumulate valence-only XC.
                       vxc1ab_val(1) = Paw_ij(iat)%dijxc_val(klmn1,  iab)
                       vxc1ab_val(2) = Paw_ij(iat)%dijxc_val(klmn1+1,iab)
                       tmp_xcval(1,iab) = tmp_xcval(1,iab) + (vxc1ab_val(1)*re_p - vxc1ab_val(2)*im_p)*fact
                       tmp_xcval(2,iab) = tmp_xcval(2,iab) + (vxc1ab_val(2)*re_p + vxc1ab_val(1)*im_p)*fact
                     end if

                     ! * In GW, dijhartree is always real.
                     if (Mflags%has_vhartree==1) then ! * Accumulate Hartree term of the PAW Hamiltonian.
                       if (iab==1.or.iab==2) then
                         DijH = Paw_ij(iat)%dijhartree(klmn)
                         tmp_H(1,iab) = tmp_H(1,iab) + DijH*re_p*fact
                         tmp_H(2,iab) = tmp_H(2,iab) + DijH*im_p*fact
                       end if
                     end if

                     ! TODO "ADD LDA+U and SO"
                     ! check this part
                     if (Mflags%has_vu==1) then
                       if (Pawtab(itypat)%usepawu/=0) then
                         ! Accumulate the U term of the PAW Hamiltonian (only onsite AE contribution)
                         dijU(1)=Paw_ij(iat)%dijU(klmn1  ,iab)
                         dijU(2)=Paw_ij(iat)%dijU(klmn1+1,iab)
                         tmp_U(1,iab) = tmp_U(1,iab) + (dijU(1)*re_p - dijU(2)*im_p)*fact
                         tmp_U(2,iab) = tmp_U(2,iab) + (dijU(2)*re_p + dijU(1)*im_p)*fact
                       end if
                     end if

                   end if
                 end do !iab

                 klmn1=klmn1+cplex_dij

               end do !ilmn
             end do !jlmn
           end do !iat

           ! ========================================
           ! ==== Add to plane wave contribution ====
           ! ========================================
           if (nspinor==1) then

             if (Mflags%has_hbare==1)    &
&              Mels%hbare(ib,jb,ik_ibz,is) = Mels%hbare(ib,jb,ik_ibz,is) + DCMPLX(tmp_h0ij(1,1),tmp_h0ij(2,1))

             if (Mflags%has_vxc==1)      &
&              Mels%vxc(ib,jb,ik_ibz,is) = Mels%vxc(ib,jb,ik_ibz,is) + DCMPLX(tmp_xc(1,1),tmp_xc(2,1))

             if (Mflags%has_vxcval==1)   &
&              Mels%vxcval(ib,jb,ik_ibz,is) = Mels%vxcval(ib,jb,ik_ibz,is) + DCMPLX(tmp_xcval(1,1),tmp_xcval(2,1))

             if (Mflags%has_vxcval_hybrid==1)   &
&              Mels%vxcval_hybrid(ib,jb,ik_ibz,is) = Mels%vxcval_hybrid(ib,jb,ik_ibz,is) + DCMPLX(tmp_xcval(1,1),tmp_xcval(2,1))

             if (Mflags%has_vhartree==1) &
&              Mels%vhartree(ib,jb,ik_ibz,is) = Mels%vhartree(ib,jb,ik_ibz,is) + DCMPLX(tmp_H (1,1),tmp_H (2,1))

             if (Mflags%has_vu==1)       &
&              Mels%vu(ib,jb,ik_ibz,is) = DCMPLX(tmp_U(1,1),tmp_U(2,1))

             if (Mflags%has_sxcore==1)   &
&              Mels%sxcore(ib,jb,ik_ibz,is) = DCMPLX(tmp_sigcx(1,1),tmp_sigcx(2,1))

           else

             if (Mflags%has_hbare==1)    &
&              Mels%hbare(ib,jb,ik_ibz,:) = Mels%hbare(ib,jb,ik_ibz,:) + DCMPLX(tmp_h0ij(1,:),tmp_h0ij(2,:))

             if (Mflags%has_vxc==1)      &
&              Mels%vxc(ib,jb,ik_ibz,:) = Mels%vxc(ib,jb,ik_ibz,:) + DCMPLX(tmp_xc(1,:),tmp_xc(2,:))

             if (Mflags%has_vxcval==1)   &
&              Mels%vxcval(ib,jb,ik_ibz,:) = Mels%vxcval(ib,jb,ik_ibz,:) + DCMPLX(tmp_xcval(1,:),tmp_xcval(2,:))

             if (Mflags%has_vxcval_hybrid==1)   &
&              Mels%vxcval_hybrid(ib,jb,ik_ibz,:) = Mels%vxcval_hybrid(ib,jb,ik_ibz,:) + DCMPLX(tmp_xcval(1,:),tmp_xcval(2,:))

             if (Mflags%has_vhartree==1) &
&              Mels%vhartree(ib,jb,ik_ibz,:) = Mels%vhartree(ib,jb,ik_ibz,:) + DCMPLX(tmp_H (1,:),tmp_H (2,:))

             if (Mflags%has_vu==1)       &
&              Mels%vu(ib,jb,ik_ibz,:) = DCMPLX(tmp_U(1,:),tmp_U(2,:))
           end if

         end do !ib
       end do !jb

     end do !is
   end do !ikc

   ABI_FREE(dimlmn)
   call pawcprj_free(Cprj_b1ks)
   ABI_DATATYPE_DEALLOCATE(Cprj_b1ks)
   call pawcprj_free(Cprj_b2ks)
   ABI_DATATYPE_DEALLOCATE(Cprj_b2ks)
 end if !PAW

 ABI_FREE(bbp_ks_distrb)

 ! Sum up contributions on each node ===
 ! Set the corresponding has_* flags to 2.
 call melements_mpisum(Mels, wfd%comm)

 ! Reconstruct lower triangle.
 call melements_herm(Mels)

 ABI_FREE(kcalc2ibz)
 call destroy_mpi_enreg(MPI_enreg_seq)

 DBG_EXIT("COLL")

end subroutine calc_vhxc_me
!!***

end module m_vhxc_me
!!***
