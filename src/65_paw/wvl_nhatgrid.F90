!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_nhatgrid
!! NAME
!! wvl_nhatgrid
!!
!! FUNCTION
!! Determine parts of the rectangular (fine) grid that are contained
!! inside spheres around atoms (used to compute n_hat density).
!! If corresponding option is selected, compute also g_l(r)*Y_lm(r)
!! (and derivatives) on this grid (g_l=radial shape function).
!!
!! COPYRIGHT
!! Copyright (C) 2011-2017 ABINIT group (T Rangel)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
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
!!      afterscfloop,scfcv
!!
!! CHILDREN
!!      pawgylm,pawrfgd_wvl,timab,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

!PENDING: ADD PARALELLISM OVER ATOMS:
!COPY NHATGRID

subroutine wvl_nhatgrid(atindx1,geocode,h,i3s,natom,natom_tot,&
& nattyp,ntypat,n1,n1i,n2,n2i,n3,n3pi,optcut,optgr0,optgr1,optgr2,optrad,&
& pawfgrtab,pawtab,psppar,rprimd,shift,xred)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_pawtab,       only : pawtab_type
 use m_pawfgrtab,    only : pawfgrtab_type
 use m_paw_finegrid, only : pawgylm, pawrfgd_wvl

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_nhatgrid'
 use interfaces_18_timing
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: i3s,natom,natom_tot,ntypat,optcut,optgr0,optgr1,optgr2,optrad
 integer,intent(in) :: n1,n2,n3,n1i,n2i,n3pi,shift
 real(dp),intent(in) :: h(3)
 character(1),intent(in) :: geocode
!integer,intent(in),optional :: mpi_comm_wvl
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat)
 real(dp),intent(in) :: psppar(0:4,0:6,ntypat),rprimd(3,3)
 real(dp),intent(inout) :: xred(3,natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ------------------------------
!scalars
!buffer to be added at the end of the last dimension of an array to control bounds_check
 integer :: iat,iatm,iatom,iatom_tot,itypat,lm_size,nfgd
 real(dp) :: rloc,rshp,xcart(3,natom)
!arrays
 integer,allocatable :: ifftsph_tmp(:)
 real(dp) :: hh(3) !fine grid spacing for wavelets
 real(dp) :: tsec(2)
 real(dp),allocatable :: rfgd_tmp(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

#if !defined HAVE_BIGDFT
 BIGDFT_NOTENABLED_ERROR()
#endif

 call timab(559,1,tsec)

!Set up parallelism for wvl
!for debug: use me_wvl=xmpi_comm_rank(MPI_COMM_WORLD)
!if (present(mpi_comm_wvl)) then
!me_wvl=xmpi_comm_rank(mpi_comm_wvl)
!nproc_wvl=xmpi_comm_size(mpi_comm_wvl)
!else
!me_wvl=0;nproc_wvl=1
!end if
!Pending: parallelism over atoms: see nhatgrid

 if (natom_tot<natom) then   ! This test has to be remove when natom_tot is used
   MSG_BUG(' natom_tot<natom !')
 end if

!Fine grid
 hh(:)=0.5d0*h(:)

!Compute xcart from xred
 call xred2xcart(natom,rprimd,xcart,xred)

!Loop over types of atom
 iatm=0
 do itypat=1,ntypat
   
   rloc=psppar(0,0,itypat)  
   if (optcut==1) then
     rshp=pawtab(itypat)%rpaw
   else
     rshp=pawtab(itypat)%rshp
   end if

!  Loop over atoms
   do iat=1,nattyp(itypat)
     iatm=iatm+1;iatom=atindx1(iatm)
     iatom_tot=iatom; !if (paral_atom) iatom_tot=my_atmtab(iatom)
     lm_size=pawfgrtab(iatom)%l_size**2

!    Determine FFT points and r-R vectors around the atom
     call pawrfgd_wvl(geocode,hh,ifftsph_tmp,i3s,n1,n1i,n2,n2i,n3,n3pi,nfgd,rshp,rloc,&
&     rfgd_tmp,shift,xcart(:,iatom_tot))

!    Allocate arrays defining sphere (and related data) around current atom
     if (allocated(pawfgrtab(iatom)%ifftsph)) then
       ABI_DEALLOCATE(pawfgrtab(iatom)%ifftsph)
     end if
     ABI_ALLOCATE(pawfgrtab(iatom)%ifftsph,(nfgd))
     pawfgrtab(iatom)%nfgd=nfgd
     pawfgrtab(iatom)%ifftsph(1:nfgd)=ifftsph_tmp(1:nfgd)

     if (optrad==1) then
       if (allocated(pawfgrtab(iatom)%rfgd)) then
         ABI_DEALLOCATE(pawfgrtab(iatom)%rfgd)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%rfgd,(3,nfgd))
       pawfgrtab(iatom)%rfgd_allocated=1
       pawfgrtab(iatom)%rfgd(1:3,1:nfgd)=rfgd_tmp(1:3,1:nfgd)
     end if

     if (optgr0==1) then
       if (allocated(pawfgrtab(iatom)%gylm)) then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylm)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylm,(nfgd,lm_size))
       pawfgrtab(iatom)%gylm_allocated=1
     end if

     if (optgr1==1) then
       if (allocated(pawfgrtab(iatom)%gylmgr)) then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr,(3,nfgd,lm_size))
       pawfgrtab(iatom)%gylmgr_allocated=1
     end if

     if (optgr2==1) then
       if (allocated(pawfgrtab(iatom)%gylmgr2)) then
         ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr2)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%gylmgr2,(6,nfgd,lm_size))
       pawfgrtab(iatom)%gylmgr2_allocated=1
     end if

!    Calculate g_l(r-R)*Y_lm(r-R) for each r around the atom R
     if (optgr0+optgr1+optgr2>0) then
       call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,pawfgrtab(iatom)%gylmgr2,&
&       lm_size,nfgd,optgr0,optgr1,optgr2,pawtab(itypat),rfgd_tmp(:,1:nfgd))
     end if

!    End loops over types/atoms
     ABI_DEALLOCATE(ifftsph_tmp)
     ABI_DEALLOCATE(rfgd_tmp)
   end do
 end do

 call timab(559,2,tsec)

 DBG_EXIT("COLL")

end subroutine wvl_nhatgrid
!!***
