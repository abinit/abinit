!{\src2tex{textfont=tt}}
!!****f* ABINIT/nhatgrid
!! NAME
!! nhatgrid
!!
!! FUNCTION
!! Determine parts of the rectangular (fine) grid that are contained
!! inside spheres around atoms (used to compute n_hat density).
!! If corresponding option is selected, compute also g_l(r)*Y_lm(r)
!! (and derivatives) on this grid (g_l=radial shape function).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  distribfft<type(distribfft_type)>=--optional-- contains all the informations related
!!                                    to the FFT parallelism and plane sharing
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  comm_fft=--optional-- MPI communicator over FFT components
!!  my_natom=number of atoms treated by current processor
!!  natom=total number of atoms in cell
!!  nattyp(ntypat)= # atoms of each type.
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms in unit cell
!!  optcut= option for the cut-off radius of spheres:
!!          if optcut=0, cut-off radius=pawtab%rshp=cut-off radius of compensation charge
!!          if optcut=1, cut-off radius=pawtab%rpaw=radius of PAW augmentation regions
!!  optgr0= 1 if g_l(r)*Y_lm(r) are computed
!!  optgr1= 1 if first derivatives of g_l(r)*Y_lm(r) are computed
!!  optgr2= 1 if second derivatives of g_l(r)*Y_lm(r) are computed
!!  optrad= 1 if vectors (r-r_atom) on the fine grid around atoms have to be stored
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  typat(natom)=type (integer) for each atom
!!  typord=1 if the output is ordered by type of atoms, 0 otherwise
!!  ucvol=unit cell volume in bohr**3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  pawfgrtab(natom)%ifftsph(nfgd)=FFT index (fine grid) of a points in paw spheres around each atom
!!  pawfgrtab(natom)%nfgd= number of (fine grid) FFT points in paw spheres around atoms
!!  if (optgr0==1)
!!    pawfgrtab(natom)%gylm(nfgd,l_size**2)= g_l(r)*Y_lm(r) around each atom
!!  if (optgr1==1)
!!    pawfgrtab(natom)%gylmgr(3,nfgd,l_size**2)= derivatives of g_l(r)*Y_lm(r) wrt cart. coordinates
!!  if (optgr2==1)
!!    pawfgrtab(natom)%gylmgr2(6,nfgd,l_size**2)= second derivatives of g_l(r)*Y_lm(r) wrt cart. coordinates
!!  if (optrad==1)
!!    pawfgrtab(natom)%rfgd(3,nfgd)= coordinates of r-r_atom around each atom
!!
!! PARENTS
!!      afterscfloop,bethe_salpeter,classify_bands,denfgr,exc_plot,m_wfd
!!      pawmkaewf,respfn,scfcv,screening,sigma,wfk_analyze
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawgylm,pawrfgd_fft,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine nhatgrid(atindx1,gmet,my_natom,natom,nattyp,ngfft,ntypat,&
& optcut,optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,typat,ucvol,xred, &
& mpi_atmtab,comm_atom,comm_fft,distribfft,typord) ! optional arguments (parallelism)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_xmpi,         only : xmpi_comm_self,xmpi_comm_rank,xmpi_comm_size
 use m_pawtab,       only : pawtab_type
 use m_pawfgrtab,    only : pawfgrtab_type
 use m_paw_finegrid, only : pawgylm, pawrfgd_fft
 use m_paral_atom,   only : get_my_atmtab, free_my_atmtab
 use m_distribfft,   only : distribfft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nhatgrid'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: my_natom,natom,ntypat,optcut,optgr0,optgr1,optgr2,optrad
 integer,optional,intent(in) :: comm_atom,comm_fft,typord
 real(dp),intent(in) :: ucvol
 type(distribfft_type),optional,target,intent(in)  :: distribfft
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 integer,intent(in),target :: atindx1(natom),nattyp(ntypat)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gmet(3,3),rprimd(3,3),xred(3,natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ------------------------------
!scalars
 integer :: i3,iat,iatm,iatom,iatom_,iatom_tot,itypat,lm_size,me_fft,my_comm_atom,n1,n2,n3,nfgd
 logical :: grid_found,my_atmtab_allocated,paral_atom
 real(dp) :: rcut
 character(len=500) :: msg
!arrays
 integer,allocatable :: ifftsph_tmp(:)
 integer,pointer :: my_atindx1(:),my_atmtab(:),my_nattyp(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: rfgd_tmp(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(559,1,tsec)
 if (my_natom==0) return

!Set up parallelism over FFT
 me_fft=0
 if (present(comm_fft)) then
   me_fft=xmpi_comm_rank(comm_fft)
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
& my_natom_ref=my_natom)
 if (paral_atom) then
   ABI_ALLOCATE(my_atindx1,(natom))
   ABI_ALLOCATE(my_nattyp,(ntypat))
   my_atindx1(:)=0;my_nattyp(:)=0
   iat=1
   do itypat=1,ntypat
     if (my_natom>0) then
       do iatom=1,my_natom
         if(typat(my_atmtab(iatom))==itypat)then
           my_nattyp(itypat)=my_nattyp(itypat)+1
           my_atindx1(iat)=iatom
           iat=iat+1
         end if
       end do
     end if
   end do
 else
   my_atindx1 => atindx1
   my_nattyp => nattyp
 end if

!Get the distrib associated with this fft_grid
 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3)
 if (present(distribfft)) then
   grid_found=.false.
   if (n2 == distribfft%n2_coarse) then
     if (n3== size(distribfft%tab_fftdp3_distrib)) then
       fftn3_distrib => distribfft%tab_fftdp3_distrib
       ffti3_local => distribfft%tab_fftdp3_local
       grid_found=.true.
     end if
   end if
   if (n2 == distribfft%n2_fine) then
     if (n3 == size(distribfft%tab_fftdp3dg_distrib)) then
       fftn3_distrib => distribfft%tab_fftdp3dg_distrib
       ffti3_local => distribfft%tab_fftdp3dg_local
       grid_found = .true.
     end if
   end if
   if (.not.(grid_found)) then
     msg='Unable to find an allocated distrib for this fft grid!'
     MSG_BUG(msg)
   end if
 else
   ABI_ALLOCATE(fftn3_distrib,(n3))
   ABI_ALLOCATE(ffti3_local,(n3))
   fftn3_distrib=0;ffti3_local=(/(i3,i3=1,n3)/)
 end if

!Loop over types of atom
!-------------------------------------------
 iatm=0
 do itypat=1,ntypat

   if (optcut==1) then
     rcut=pawtab(itypat)%rpaw
   else
     rcut=pawtab(itypat)%rshp
   end if

!  Loop over atoms
!  -------------------------------------------
   do iat=1,my_nattyp(itypat)
     iatm=iatm+1;iatom=my_atindx1(iatm)
     iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
     iatom_=iatom;if(present(typord)) iatom_=merge(iatm,iatom,typord==1)
     lm_size=pawfgrtab(iatom_)%l_size**2

!    ------------------------------------------------------------------
!    A-Determine FFT points and r-R vectors around the atom
!    ------------------------------------------------------------------

     call pawrfgd_fft(ifftsph_tmp,gmet,n1,n2,n3,nfgd,rcut,rfgd_tmp,rprimd,ucvol,&
&     xred(:,iatom_tot),fft_distrib=fftn3_distrib,fft_index=ffti3_local,me_fft=me_fft)

!    Allocate arrays defining sphere (and related data) around current atom
     if (allocated(pawfgrtab(iatom_)%ifftsph)) then
       ABI_DEALLOCATE(pawfgrtab(iatom_)%ifftsph)
     end if
     ABI_ALLOCATE(pawfgrtab(iatom_)%ifftsph,(nfgd))
     pawfgrtab(iatom_)%nfgd=nfgd
     pawfgrtab(iatom_)%ifftsph(1:nfgd)=ifftsph_tmp(1:nfgd)

     if (optrad==1) then
       if (allocated(pawfgrtab(iatom_)%rfgd))  then
         ABI_DEALLOCATE(pawfgrtab(iatom_)%rfgd)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom_)%rfgd,(3,nfgd))
       pawfgrtab(iatom_)%rfgd_allocated=1
       pawfgrtab(iatom_)%rfgd(1:3,1:nfgd)=rfgd_tmp(1:3,1:nfgd)
     end if

     if (optgr0==1) then
       if (allocated(pawfgrtab(iatom_)%gylm))  then
         ABI_DEALLOCATE(pawfgrtab(iatom_)%gylm)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom_)%gylm,(nfgd,lm_size))
       pawfgrtab(iatom_)%gylm_allocated=1
     end if

     if (optgr1==1) then
       if (allocated(pawfgrtab(iatom_)%gylmgr))  then
         ABI_DEALLOCATE(pawfgrtab(iatom_)%gylmgr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom_)%gylmgr,(3,nfgd,lm_size))
       pawfgrtab(iatom_)%gylmgr_allocated=1
     end if

     if (optgr2==1) then
       if (allocated(pawfgrtab(iatom_)%gylmgr2))  then
         ABI_DEALLOCATE(pawfgrtab(iatom_)%gylmgr2)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom_)%gylmgr2,(6,nfgd,lm_size))
       pawfgrtab(iatom_)%gylmgr2_allocated=1
     end if

!    ------------------------------------------------------------------
!    B-Calculate g_l(r-R)*Y_lm(r-R) for each r around the atom R
!    ------------------------------------------------------------------
     if (optgr0+optgr1+optgr2>0) then
       call pawgylm(pawfgrtab(iatom_)%gylm,pawfgrtab(iatom_)%gylmgr,pawfgrtab(iatom_)%gylmgr2,&
&       lm_size,nfgd,optgr0,optgr1,optgr2,pawtab(itypat),rfgd_tmp(:,1:nfgd))
     end if

!    End loops over types/atoms
!    -------------------------------------------
     ABI_DEALLOCATE(ifftsph_tmp)
     ABI_DEALLOCATE(rfgd_tmp)
   end do
 end do

!Destroy atom tables used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)
 if (paral_atom) then
   ABI_DEALLOCATE(my_atindx1)
   ABI_DEALLOCATE(my_nattyp)
 end if

 if (.not.present(distribfft)) then
   ABI_DEALLOCATE(fftn3_distrib)
   ABI_DEALLOCATE(ffti3_local)
 end if

 call timab(559,2,tsec)

 DBG_EXIT("COLL")

end subroutine nhatgrid

!!***
