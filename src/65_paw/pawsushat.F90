!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawsushat
!! NAME
!! pawsushat
!!
!! FUNCTION
!! PAW only, for susceptibility matrix:
!! Compute contribution to the product of two wavefunctions (exchange charge density)
!! from hat (compensation charge) density (in reciprocal space and eventually in real space):
!!    sushat_{ij,R}(g)=Sum_{L}[Q^L_ijR(g)]
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms, inverse of atindx
!!  cprj_k(natom,nspinor*nband_k)= wave functions projected with non-local projectors:
!!                                 cprj_k=<p_i|Cnk> where p_i is a non-local projector.
!!                                 WARNING: cprj(iatom,:) ARE SORTED BY ATOM TYPE !!!
!!  distribfft<type(distribfft_type)>=--optional-- contains infos related to FFT parallelism
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for small FFT sphere.
!!  gylmg_diel(npwdiel,lmax_diel**2,ntypat)= -PAW only- Fourier transform of g_l(r).Y_ml(r) shape functions
!!  iband1,iband2= indices of the bands concerned with
!!  ispinor1,ispinor2= indices of spinorial components concerned with
!!  istwf_k=input option parameter that describes the storage of wfs
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  lmax_diel=1+max. value of l angular momentum used for dielectric matrix
!!  me_g0=--optional-- 1 if the current process treat the g=0 plane-wave (only needed when comm_fft is present)
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of the dielectric matrix
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  comm_fft=--optional-- MPI communicator over FT components
!!  natom=number of atoms in cell
!!  nband=number of bands at this k point for that spin polarization
!!  ndiel4,ndiel5,ndiel6= FFT dimensions, modified to avoid cache trashing
!!  nfftdiel=number of FFT grid points for the small (diel) grid
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in unit cell.
!!  optreal=0 if WF product has to be output in reciprocal space
!!          1 if WF product has to be output in real space
!!  paral_kgb=--optional-- 1 if "band-FFT" parallelism is activated (only needed when comm_fft is present)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph3d_diel(2,npwdiel,natom*usepaw)=3-dim structure factors, for each atom and plane wave, for dielectric matrix
!!  typat(natom)=type (integer) for each atom
!!
!! SIDE EFFECTS
!!  === if optreal=0
!!  wfprod(2,npwdiel)=PAW contrib. to product of two wavefunctions (iband1,iband2):
!!                    is added (in reciprocal space)
!!  === if optreal=1
!!  wfraug(2,ndiel4,ndiel5,ndiel6)=PAW contrib. to product of two wavefunctions (iband1,iband2)
!!                                 is added (in real space)
!!
!! PARENTS
!!      susk,suskmm
!!
!! CHILDREN
!!      destroy_distribfft,fourwf,free_my_atmtab,get_my_atmtab
!!      init_distribfft_seq,initmpi_seq,set_mpi_enreg_fft,unset_mpi_enreg_fft
!!      xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawsushat(atindx,cprj_k,gbound_diel,gylmg_diel,iband1,iband2,ispinor1,ispinor2,istwf_k,kg_diel,&
&                    lmax_diel,mgfftdiel,natom,nband,ndiel4,ndiel5,ndiel6,&
&                    ngfftdiel,npwdiel,nspinor,ntypat,optreal,&
&                    pawang,pawtab,ph3d_diel,typat,wfprod,wfraug, &
&                    mpi_atmtab,comm_atom,comm_fft,me_g0,paral_kgb,distribfft) ! optional arguments (parallelism)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_distribfft
 use m_errors
 use m_xmpi, only : xmpi_comm_self,xmpi_sum

 use m_pawang,     only : pawang_type
 use m_pawtab,     only : pawtab_type
 use m_pawcprj,    only : pawcprj_type
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab
 use m_mpinfo,     only : set_mpi_enreg_fft, unset_mpi_enreg_fft

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawsushat'
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: iband1,iband2,ispinor1,ispinor2,istwf_k,lmax_diel,mgfftdiel
 integer,intent(in) :: natom,nband,ndiel4,ndiel5,ndiel6,npwdiel,nspinor
 integer,intent(in) :: ntypat,optreal
 integer,optional,intent(in) :: me_g0,comm_atom,comm_fft,paral_kgb
 type(distribfft_type),optional,intent(in),target :: distribfft
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: atindx(natom),gbound_diel(2*mgfftdiel+8,2)
 integer,intent(in) :: kg_diel(3,npwdiel),ngfftdiel(18),typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gylmg_diel(npwdiel,lmax_diel**2,ntypat)
 real(dp),intent(in) :: ph3d_diel(2,npwdiel,natom)
 real(dp),intent(inout) :: wfprod(2,npwdiel*(1-optreal))
 real(dp),intent(inout) :: wfraug(2,ndiel4,ndiel5,ndiel6*optreal)
 type(pawcprj_type),intent(in) :: cprj_k(natom,nspinor*nband)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: cplex,iatm,iatom,iatom_tot,ibsp1,ibsp2,ierr,il,ilmn,ils,ilslm,ipw
 integer :: itypat,j0lmn,jlmn,klm,klmn,lmax,lmin,mm,my_comm_atom,my_comm_fft,my_natom,paral_kgb_fft,tim_fourwf
 real(dp) :: phil1,phil2,sgn,weight_dum,wf1,wf2
 logical :: my_atmtab_allocated,parity,paral_atom
 type(distribfft_type),pointer :: my_distribfft
 type(mpi_type) :: mpi_enreg_fft
!arrays
 integer,pointer :: my_atmtab(:)
 real(dp) :: ro(2),ro_ql(2)
 real(dp),allocatable :: dummy(:,:),wfprod_paw(:,:),wfraug_paw(:,:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 if (present(comm_fft)) then
   if ((.not.present(paral_kgb)).or.(.not.present(me_g0))) then
     MSG_BUG('Need paral_kgb and me_g0 with comm_fft !')
   end if
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom)
 my_natom=natom;if (paral_atom) my_natom=size(my_atmtab)

 cplex=1;if (istwf_k>1) cplex=2
 ABI_ALLOCATE(wfprod_paw,(2,npwdiel))
 wfprod_paw(:,:)=zero
 ibsp1=(iband1-1)*nspinor+ispinor1
 ibsp2=(iband2-1)*nspinor+ispinor2

!------------------------------------------------------------------------
!----- Loop over atoms
!------------------------------------------------------------------------
 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
   iatm=atindx(iatom_tot)
   itypat=typat(iatom_tot)

!  ------------------------------------------------------------------------
!  ----- Loop over ij channels (basis components)
!  ------------------------------------------------------------------------
   do jlmn=1,pawtab(itypat)%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       klm =pawtab(itypat)%indklmn(1,klmn)
       lmin=pawtab(itypat)%indklmn(3,klmn)
       lmax=pawtab(itypat)%indklmn(4,klmn)

       ro(1)=cprj_k(iatm,ibsp1)%cp(1,ilmn)*cprj_k(iatm,ibsp2)%cp(1,jlmn)
       if (cplex==2) then
         ro(1)=ro(1)+cprj_k(iatm,ibsp1)%cp(2,ilmn)*cprj_k(iatm,ibsp2)%cp(2,jlmn)
         ro(2)=cprj_k(iatm,ibsp1)%cp(2,ilmn)*cprj_k(iatm,ibsp2)%cp(1,jlmn) &
&         -cprj_k(iatm,ibsp1)%cp(1,ilmn)*cprj_k(iatm,ibsp2)%cp(2,jlmn)
       end if
       ro(1:cplex)=ro(1:cplex)*pawtab(itypat)%dltij(klmn)

       do ils=lmin,lmax,2
         il=mod(ils,4);parity=(mod(il,2)==0)
         sgn=one;if (il>1) sgn=-one

         do mm=-ils,ils
           ilslm=ils*ils+ils+mm+1
           if (pawang%gntselect(ilslm,klm)>0) then

             ro_ql(1:cplex)=pawtab(itypat)%qijl(ilslm,klmn)*ro(1:cplex)

!            Compute: Sum_{ijR} [ cpi* cpj qij^l (-i)^l g_l(g) S_lm(g) ]

             if (cplex==1) then
               if (parity) then
                 do ipw=1,npwdiel
                   phil1= sgn*ph3d_diel(1,ipw,iatm)     ! (i)^l.exp(i.g.R)
                   phil2= sgn*ph3d_diel(2,ipw,iatm)
                   wf1= phil1*ro_ql(1)                  ! cpi* cpj qij^l (-i)^l.exp(-i.g.R)
                   wf2=-phil2*ro_ql(1)
                   wfprod_paw(1,ipw)=wfprod_paw(1,ipw)+wf1*gylmg_diel(ipw,ilslm,itypat)
                   wfprod_paw(2,ipw)=wfprod_paw(2,ipw)+wf2*gylmg_diel(ipw,ilslm,itypat)
                 end do
               else
                 do ipw=1,npwdiel
                   phil1=-sgn*ph3d_diel(2,ipw,iatm)  ! (i)^l.exp(i.g.R)
                   phil2= sgn*ph3d_diel(1,ipw,iatm)
                   wf1= phil1*ro_ql(1)               ! cpi* cpj qij^l (-i)^l.exp(-i.g.R)
                   wf2=-phil2*ro_ql(1)
                   wfprod_paw(1,ipw)=wfprod_paw(1,ipw)+wf1*gylmg_diel(ipw,ilslm,itypat)
                   wfprod_paw(2,ipw)=wfprod_paw(2,ipw)+wf2*gylmg_diel(ipw,ilslm,itypat)
                 end do
               end if

             else

               if (parity) then
                 do ipw=1,npwdiel
                   phil1= sgn*ph3d_diel(1,ipw,iatm)     ! (i)^l.exp(i.g.R)
                   phil2= sgn*ph3d_diel(2,ipw,iatm)
                   wf1=phil1*ro_ql(1)+phil2*ro_ql(2)    ! cpi* cpj qij^l (-i)^l.exp(-i.g.R)
                   wf2=phil1*ro_ql(2)-phil2*ro_ql(1)
                   wfprod_paw(1,ipw)=wfprod_paw(1,ipw)+wf1*gylmg_diel(ipw,ilslm,itypat)
                   wfprod_paw(2,ipw)=wfprod_paw(2,ipw)+wf2*gylmg_diel(ipw,ilslm,itypat)
                 end do
               else
                 do ipw=1,npwdiel
                   phil1=-sgn*ph3d_diel(2,ipw,iatm)     ! (i)^l.exp(i.g.R)
                   phil2= sgn*ph3d_diel(1,ipw,iatm)
                   wf1=phil1*ro_ql(1)+phil2*ro_ql(2)    ! cpi* cpj qij^l (-i)^l.exp(-i.g.R)
                   wf2=phil1*ro_ql(2)-phil2*ro_ql(1)
                   wfprod_paw(1,ipw)=wfprod_paw(1,ipw)+wf1*gylmg_diel(ipw,ilslm,itypat)
                   wfprod_paw(2,ipw)=wfprod_paw(2,ipw)+wf2*gylmg_diel(ipw,ilslm,itypat)
                 end do
               end if

             end if
           end if
         end do
       end do

!      ----- End loop over ij channels
     end do
   end do

!  ----- End loop over atoms
 end do

!Reduction in case of parallelism over atoms
 if (paral_atom) then
   call xmpi_sum(wfprod_paw,my_comm_atom,ierr)
 end if

 if (optreal==0) then

!  === Output in reciprocal space
   wfprod(:,:)=wfprod(:,:)+wfprod_paw(:,:)

 else
!  === Output in reciprocal space
   tim_fourwf=17;weight_dum=0
!  Create fake mpi_enreg to wrap fourdp
   if (present(distribfft)) then
     my_distribfft => distribfft
   else
     ABI_DATATYPE_ALLOCATE(my_distribfft,)
     call init_distribfft_seq(my_distribfft,'c',ngfftdiel(2),ngfftdiel(3),'fourwf')
   end if
   call initmpi_seq(mpi_enreg_fft)
   ABI_DATATYPE_DEALLOCATE(mpi_enreg_fft%distribfft)
   if (present(comm_fft)) then
     call set_mpi_enreg_fft(mpi_enreg_fft,comm_fft,my_distribfft,me_g0,paral_kgb)
     my_comm_fft=comm_fft;paral_kgb_fft=paral_kgb
   else
     my_comm_fft=xmpi_comm_self;paral_kgb_fft=0;
     mpi_enreg_fft%distribfft => my_distribfft
   end if
!  do FFT
   ABI_ALLOCATE(wfraug_paw,(2,ndiel4,ndiel5,ndiel6))
   call fourwf(1,dummy,wfprod_paw,dummy,wfraug_paw,gbound_diel,gbound_diel,&
&   istwf_k,kg_diel,kg_diel,mgfftdiel,mpi_enreg_fft,1,ngfftdiel,1,npwdiel,&
&   ndiel4,ndiel5,ndiel6,0,paral_kgb_fft,tim_fourwf,weight_dum,weight_dum)
   wfraug(:,:,:,:)=wfraug(:,:,:,:)+wfraug_paw(:,:,:,:)
   ABI_DEALLOCATE(wfraug_paw)
   call unset_mpi_enreg_fft(mpi_enreg_fft)
   if (.not.present(distribfft)) then
     call destroy_distribfft(my_distribfft)
     ABI_DATATYPE_DEALLOCATE(my_distribfft)
   end if
 end if

 ABI_DEALLOCATE(wfprod_paw)

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine pawsushat
!!***
