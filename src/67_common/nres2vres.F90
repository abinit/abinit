!{\src2tex{textfont=tt}}
!!****f* ABINIT/nres2vres
!!
!! NAME
!! nres2vres
!!
!! FUNCTION
!! Convert a density residual into a potential residual
!! using a first order formula:
!!     V^res(r)=dV/dn.n^res(r)
!!             =V_hartree(n^res)(r) + Kxc.n^res(r)
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dtset <type(dataset_type)>=all input variables in this dataset
!!  | icoulomb=0 periodic treatment of Hartree potential, 1 use of Poisson solver
!!  | natom= number of atoms in cell
!!  | nspden=number of spin-density components
!!  | ntypat=number of atom types
!!  | typat(natom)=type (integer) for each atom
!! gsqcut=cutoff value on G**2 for sphere inside fft box
!! izero=if 1, unbalanced components of Vhartree(g) are set to zero
!! kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
!! mpi_enreg=information about MPI parallelization
!! my_natom=number of atoms treated by current processor
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT
!! nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!! nkxc=second dimension of the array kxc, see rhohxc.F90 for a description
!! nresid(nfft,nspden)= the input density residual
!! n3xccc=dimension of the xccc3d array (0 or nfft).
!! optnc=option for non-collinear magnetism (nspden=4):
!!       1: the whole 2x2 Vres matrix is computed
!!       2: only Vres^{11} and Vres^{22} are computed
!! optxc=0 if LDA part of XC kernel has only to be taken into account (even for GGA)
!!       1 if XC kernel has to be fully taken into
!!      -1 if XC kernel does not have to be taken into account
!! pawang <type(pawang_type)>=paw angular mesh and related data
!! pawfgrtab(my_natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!! pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! rhor(nfft,nspden)=electron density in real space
!!                   (used only if Kxc was not computed before)
!! rprimd(3,3)=dimensional primitive translation vectors (bohr)
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!! xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!! vresid(nfft,nspden)= the output potential residual
!!
!! PARENTS
!!      etotfor,forstr
!!
!! CHILDREN
!!      dfpt_mkvxc,dfpt_mkvxc_noncoll,fourdp,hartre,metric,pawmknhat
!!      psolver_hartree,rhohxc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine nres2vres(dtset,gsqcut,izero,kxc,mpi_enreg,my_natom,nfft,ngfft,nhat,&
&                 nkxc,nresid,n3xccc,optnc,optxc,pawang,pawfgrtab,pawrhoij,pawtab,&
&                 rhor,rprimd,usepaw,vresid,xccc3d,xred)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_xmpi
 use m_xcdata

 use m_pawang,   only : pawang_type
 use m_pawtab,   only : pawtab_type
 use m_pawfgrtab,only : pawfgrtab_type
 use m_pawrhoij, only : pawrhoij_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nres2vres'
 use interfaces_41_geometry
 use interfaces_53_ffts
 use interfaces_56_xc
 use interfaces_62_poisson
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: izero,my_natom,n3xccc,nfft,nkxc,optnc,optxc,usepaw
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: kxc(nfft,nkxc),nresid(nfft,dtset%nspden)
 real(dp),intent(in) :: rhor(nfft,dtset%nspden),rprimd(3,3),xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp),intent(inout) :: nhat(nfft,dtset%nspden*usepaw)
 real(dp),intent(out) :: vresid(nfft,dtset%nspden)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom*usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: cplex,ider,idir,ipert,ispden,nhatgrdim,nkxc_cur,option,me,nproc,comm,usexcnhat
 logical :: has_nkxc_gga
 real(dp) :: dum,energy,m_norm_min,ucvol,vxcavg
 character(len=500) :: message
 type(xcdata_type) :: xcdata
!arrays
 integer :: nk3xc
 real(dp) :: dummy6(6),gmet(3,3),gprimd(3,3),qq(3),rmet(3,3)
 real(dp),allocatable :: dummy(:),kxc_cur(:,:),nhatgr(:,:,:)
 real(dp),allocatable :: nresg(:,:),rhor0(:,:),vhres(:)

! *************************************************************************

!Compatibility tests:
 has_nkxc_gga=(nkxc==7.or.nkxc==19)

 if(optxc<-1.or.optxc>1)then
   write(message,'(a,i0)')' Wrong value for optxc ',optxc
   MSG_BUG(message)
 end if

 if((optnc/=1.and.optnc/=2).or.(dtset%nspden/=4.and.optnc/=1))then
   write(message,'(a,i0)')' Wrong value for optnc ',optnc
   MSG_BUG(message)
 end if

 if(dtset%icoulomb==1.and.optxc/=-1)then
   write(message,'(a)')' This routine is not compatible with icoulomb==1 and optxc/=-1 !'
   MSG_BUG(message)
 end if

 if(dtset%nspden==4.and.dtset%xclevel==2.and.optxc==1.and.(.not.has_nkxc_gga))then
   MSG_ERROR(' Wrong values for optxc and nkxc !')
 end if

 qq=zero
 nkxc_cur=0
 m_norm_min=EPSILON(0.0_dp)**2
 usexcnhat=0;if (usepaw==1) usexcnhat=maxval(pawtab(1:dtset%ntypat)%usexcnhat)
 if (dtset%xclevel==1.or.optxc==0) nkxc_cur= 2*min(dtset%nspden,2)-1 ! LDA: nkxc=1,3
 if (dtset%xclevel==2.and.optxc==1)nkxc_cur=12*min(dtset%nspden,2)-5 ! GGA: nkxc=7,19
 ABI_ALLOCATE(vhres,(nfft))

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Compute density residual in reciprocal space
 if (dtset%icoulomb==0) then
   ABI_ALLOCATE(nresg,(2,nfft))
   ABI_ALLOCATE(dummy,(nfft))
   dummy(:)=nresid(:,1)
   call fourdp(1,nresg,dummy,-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
   ABI_DEALLOCATE(dummy)
 end if

!For GGA, has to recompute gradients of nhat
 nhatgrdim=0
 if ((nkxc==nkxc_cur.and.has_nkxc_gga).or.(optxc==-1.and.has_nkxc_gga).or.&
& (optxc/=-1.and.nkxc/=nkxc_cur)) then
   if (usepaw==1.and.dtset%xclevel==2.and.usexcnhat>0.and.dtset%pawnhatxc>0) then
     nhatgrdim=1
     ABI_ALLOCATE(nhatgr,(nfft,dtset%nspden,3))
     ider=1;cplex=1;ipert=0;idir=0
     call pawmknhat(dum,cplex,ider,idir,ipert,izero,gprimd,my_natom,dtset%natom,&
&     nfft,ngfft,nhatgrdim,dtset%nspden,dtset%ntypat,pawang,pawfgrtab,&
&     nhatgr,nhat,pawrhoij,pawrhoij,pawtab,qq,rprimd,ucvol,dtset%usewvl,xred,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&     comm_fft=mpi_enreg%comm_fft,paral_kgb=dtset%paral_kgb,me_g0=mpi_enreg%me_g0,&
&     distribfft=mpi_enreg%distribfft,mpi_comm_wvl=mpi_enreg%comm_wvl)
   else
     ABI_ALLOCATE(nhatgr,(0,0,0))
   end if
 else
   ABI_ALLOCATE(nhatgr,(0,0,0))
 end if

 ABI_ALLOCATE(dummy,(0))
!First case: Kxc has already been computed
!-----------------------------------------
 if (nkxc==nkxc_cur.or.optxc==-1) then

!  Compute VH(n^res)(r)
   if (dtset%icoulomb == 0) then
     call hartre(1,gsqcut,izero,mpi_enreg,nfft,ngfft,dtset%paral_kgb,qq,nresg,rprimd,vhres)
   else
     comm=mpi_enreg%comm_cell
     nproc=xmpi_comm_size(comm)
     me=xmpi_comm_rank(comm)
     call psolver_hartree(energy, (/ rprimd(1,1) / dtset%ngfft(1), &
&     rprimd(2,2) / dtset%ngfft(2), rprimd(3,3) / dtset%ngfft(3) /), dtset%icoulomb, &
&     me, comm, dtset%nfft, dtset%ngfft(1:3), nproc, dtset%nscforder, dtset%nspden, &
&     nresid(:,1), vhres, dtset%usewvl)
   end if

!  Compute Kxc(r).n^res(r)
   if (optxc/=-1) then

!    Collinear magnetism or non-polarized
     if (dtset%nspden/=4) then
       call dfpt_mkvxc(1,dtset%ixc,kxc,mpi_enreg,nfft,ngfft,nhat,usepaw,nhatgr,nhatgrdim,&
&       nkxc,dtset%nspden,0,2,dtset%paral_kgb,qq,nresid,rprimd,1,vresid,dummy)
     else
!FR      call routine for Non-collinear magnetism
       ABI_ALLOCATE(rhor0,(nfft,dtset%nspden))
       rhor0(:,:)=rhor(:,:)-nresid(:,:)
       call dfpt_mkvxc_noncoll(1,dtset%ixc,kxc,mpi_enreg,nfft,ngfft,nhat,usepaw,nhatgr,nhatgrdim,&
&       nkxc,nkxc_cur,dtset%nspden,0,2,2,optxc,dtset%paral_kgb,qq,rhor0,nresid,rprimd,1,vresid,xccc3d)
       ABI_DEALLOCATE(rhor0)  
     end if

   else
     vresid=zero
   end if

 end if

!2nd case: Kxc has to be computed
!--------------------------------
 if (nkxc/=nkxc_cur.and.optxc/=-1) then

!  Has to use the "initial" density to compute Kxc
   ABI_ALLOCATE(rhor0,(nfft,dtset%nspden))
   rhor0(:,:)=rhor(:,:)-nresid(:,:)

!  Compute VH(n^res) and XC kernel (Kxc) together
   ABI_ALLOCATE(kxc_cur,(nfft,nkxc_cur))
   option=2;if (dtset%xclevel==2.and.optxc==0) option=12

   call hartre(1,gsqcut,izero,mpi_enreg,nfft,ngfft,dtset%paral_kgb,qq,nresg,rprimd,vhres)
   call xcdata_init(dtset%intxc,dtset%ixc,&
&    dtset%nelect,dtset%tphysel,dtset%usekden,dtset%vdw_xc,dtset%xc_tb09_c,dtset%xc_denpos,xcdata)

!  To be adjusted for the call to rhohxc
   nk3xc=1
   call rhohxc(energy,kxc_cur,mpi_enreg,nfft,ngfft,&
&   nhat,usepaw,nhatgr,nhatgrdim,nkxc_cur,nk3xc,dtset%nspden,n3xccc,option,dtset%paral_kgb,&
&   rhor0,rprimd,dummy6,usexcnhat,vresid,vxcavg,xccc3d,xcdata,vhartr=vhres)  !vresid=work space
   if (dtset%nspden/=4)  then
     ABI_DEALLOCATE(rhor0)
   end if

!  Compute Kxc(r).n^res(r)

!  Collinear magnetism or non-polarized
   if (dtset%nspden/=4) then
     call dfpt_mkvxc(1,dtset%ixc,kxc_cur,mpi_enreg,nfft,ngfft,nhat,usepaw,nhatgr,nhatgrdim,&
&     nkxc_cur,dtset%nspden,0,2,dtset%paral_kgb,qq,nresid,rprimd,1,vresid,dummy)
   else
!FR      call routine for Non-collinear magnetism
     ABI_ALLOCATE(rhor0,(nfft,dtset%nspden))
     rhor0(:,:)=rhor(:,:)-nresid(:,:)
     call dfpt_mkvxc_noncoll(1,dtset%ixc,kxc_cur,mpi_enreg,nfft,ngfft,nhat,usepaw,nhatgr,nhatgrdim,&
&     nkxc,nkxc_cur,dtset%nspden,0,2,2,optxc,dtset%paral_kgb,qq,rhor0,nresid,rprimd,1,vresid,xccc3d)
     ABI_DEALLOCATE(rhor0)
   end if

   ABI_DEALLOCATE(kxc_cur)
 end if

 !if (nhatgrdim>0)  then
 ABI_DEALLOCATE(nhatgr)
 !end if

!Assemble potential residual: V^res(r)=VH(n^res)(r) + Kxc(r).n^res(r)
!--------------------------------------------------------------------
 do ispden=1,dtset%nspden/optnc
   vresid(:,ispden)=vresid(:,ispden)+vhres(:)
 end do

 if (dtset%icoulomb==0)  then
   ABI_DEALLOCATE(nresg)
 end if
 ABI_DEALLOCATE(vhres)
 ABI_DEALLOCATE(dummy)

end subroutine nres2vres
!!***
