!{\src2tex{textfont=tt}}
!!****f* ABINIT/partial_dos_fractions_paw
!! NAME
!! partial_dos_fractions_paw
!!
!! FUNCTION
!!  Calculate PAW contributions to the partial DOS fractions (tetrahedron method)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (SM,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  cprj(natom,mcprj)= <p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!!  dimcprj(natom)=array of dimensions of array cprj (not ordered)
!!  dtset     structured datatype, from which one uses :
!!   iatsph(nasph)=number of atoms used to project dos
!!   kpt(3,nkpt)  =irreducible kpoints
!!   mband        =max number of bands per k-point
!!   mkmem        =number of kpoints in memory
!!   natom        =number of atoms in total
!!   natsph       =number of atoms ofor which the spherical decomposition must be done
!!   nband        =number of electronic bands for each kpoint
!!   nkpt         =number of irreducible kpoints
!!   nspinor      =number of spinor components
!!   nsppol       =1 or 2 spin polarization channels
!!  fatbands_flag =1 if pawfatbnd=1 or 2
!!  mbesslang=maximum angular momentum for Bessel function expansion
!!  mpi_enreg=informations about MPI parallelization
!!  prtdosm=option for the m-contributions to the partial DOS
!!  ndosfraction=natsph*mbesslang
!!  paw_dos_flag=option for the PAW contributions to the partial DOS
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data:
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  === If paw_dos_flag==1:
!!   dos%fractions_paw1(ikpt,iband,isppol,natom*mbesslang) = contribution to
!!       dos fractions from the PAW partial waves (phi)
!!   dos%fractions_pawt1(ikpt,iband,isppol,natom*mbesslang) = contribution to
!!       dos fractions from the PAW pseudo partial waves (phi_tild)
!!
!! SIDE EFFECTS
!!  dos%fractions(ikpt,iband,isppol,ndosfraction) = percentage of s, p, d..
!!    character on each atom for the wavefunction # ikpt,iband, isppol
!!    As input: contains only the pseudo contribution
!!    As output: contains pseudo contribution + PAW corrections
!!  == if prtdosm==1
!!  dos%fractions_m(ikpt,iband,isppol,ndosfraction*mbesslang*prtdosm) =
!!              m discretization of partial DOS fractions
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,simp_gen,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine partial_dos_fractions_paw(dos,cprj,dimcprj,dtset,mcprj,mkmem,mpi_enreg,pawrad,pawtab)

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_profiling_abi

 use m_pawrad,  only : pawrad_type, simp_gen
 use m_pawtab,  only : pawtab_type
 use m_pawcprj, only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_free
 use m_epjdos,  only : epjdos_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'partial_dos_fractions_paw'
 use interfaces_18_timing
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcprj,mkmem
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(epjdos_t),intent(inout) :: dos
!arrays
 integer,intent(in) :: dimcprj(dtset%natom)
 type(pawcprj_type),intent(in) :: cprj(dtset%natom,mcprj)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
 type(pawtab_type),target,intent(in) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer :: basis_size,cplex,iat,iatom,iband,ibg,ibsp,ierr,ikpt,il,ilang,ilmn,iln,im,iorder_cprj,ispinor,isppol
 integer :: itypat,j0lmn,j0ln,jl,jlmn,jln,jm,klmn,kln,lmn_size,me,my_nspinor,nband_k,spaceComm
 integer :: fatbands_flag,mbesslang,prtdosm,ndosfraction,paw_dos_flag
 real(dp) :: cpij
 character(len=500) :: message
!arrays
 integer ,allocatable :: dimcprj_atsph(:)
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: int1(:,:),int2(:,:),int1m2(:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:)

!******************************************************************************************

 DBG_ENTER("COLL")
 !return

 ABI_CHECK(mkmem/=0,"mkmem==0 not supported anymore!")

 fatbands_flag = dos%fatbands_flag 
 mbesslang = dos%mbesslang 
 prtdosm = dos%prtdosm
 ndosfraction = dos%ndosfraction
 paw_dos_flag = dos%paw_dos_flag

!m-decomposed DOS not compatible with PAW-decomposed DOS
 if(prtdosm>=1.and.paw_dos_flag==1) then
   message = 'm-decomposed DOS not compatible with PAW-decomposed DOS !'
   MSG_ERROR(message)
 end if

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

!Prepare some useful integrals
 basis_size=pawtab(1)%basis_size
 if (dtset%ntypat>1) then
   do itypat=1,dtset%ntypat
     basis_size=max(basis_size,pawtab(itypat)%basis_size)
   end do
 end if
 ABI_ALLOCATE(int1  ,(basis_size*(basis_size+1)/2,dtset%natsph))
 ABI_ALLOCATE(int2,(basis_size*(basis_size+1)/2,dtset%natsph))
 ABI_ALLOCATE(int1m2,(basis_size*(basis_size+1)/2,dtset%natsph))
 int1=zero;int2=zero;int1m2=zero
 do iat=1,dtset%natsph
   iatom=dtset%iatsph(iat)
   itypat= dtset%typat(iatom)
   do jln=1,pawtab(itypat)%basis_size
     j0ln=jln*(jln-1)/2
     do iln=1,jln
       kln=j0ln+iln
       call simp_gen(int1(kln,iat),pawtab(itypat)%phiphj(:,kln),pawrad(itypat))
       if (dtset%pawprtdos<2) then
         call simp_gen(int2(kln,iat),pawtab(itypat)%tphitphj(:,kln),pawrad(itypat))
         int1m2(kln,iat)=int1(kln,iat)-int2(kln,iat)
       else
         int2(kln,iat)=zero;int1m2(kln,iat)=int1(kln,iat)
       end if
     end do !iln
   end do !jln
 end do

!Antiferro case
 if (dtset%nspden==2.and.dtset%nsppol==1.and.dtset%nspinor==1) then
   int1m2(:,:)=half*int1m2(:,:)
   if (paw_dos_flag==1.or.fatbands_flag==1.or.prtdosm==2) then
     int1(:,:)=half*int1(:,:);int2(:,:)=half*int2(:,:)
   end if
 end if

!Init parallelism
 spaceComm=mpi_enreg%comm_cell
 if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_kpt
 me=mpi_enreg%me_kpt

 iorder_cprj=0

!LOOPS OVER SPINS,KPTS
 ibg=0
 do isppol =1,dtset%nsppol
   do ikpt=1,dtset%nkpt

     nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle
     cplex=2;if (dtset%istwfk(ikpt)>1) cplex=1
     ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natsph,my_nspinor*nband_k))
     ABI_ALLOCATE(dimcprj_atsph,(dtset%natsph))
     do iat=1,dtset%natsph
       dimcprj_atsph(iat)=dimcprj(dtset%iatsph(iat))
     end do
     call pawcprj_alloc(cprj_k,0,dimcprj_atsph)
     ABI_DEALLOCATE(dimcprj_atsph)

!    Extract cprj for this k-point.
! FIXME: Parallelism over atoms is buggy
     ibsp=0
     do iband=1,nband_k
       do ispinor=1,my_nspinor
         ibsp=ibsp+1
         do iat=1,dtset%natsph
           iatom=dtset%iatsph(iat)
           cprj_k(iat,ibsp)%cp(:,:)=cprj(iatom,ibsp+ibg)%cp(:,:)
         end do
       end do
     end do

!    LOOP OVER ATOMS
     do iat=1,dtset%natsph
       iatom=dtset%iatsph(iat)
       itypat= dtset%typat(iatom)
       lmn_size=pawtab(itypat)%lmn_size
       indlmn => pawtab(itypat)%indlmn

!      LOOP OVER BANDS
       do iband=1,nband_k
         if(abs(mpi_enreg%proc_distrb(ikpt,iband,isppol)-me)/=0) cycle
         ! FIXME: MPI-FFT parallelism is buggy
         !if (mod(iband, mpi_enreg%nproc_fft) /= mpi_enreg%me_fft) cycle
         ibsp=(iband-1)*my_nspinor
         do ispinor=1,my_nspinor
           ibsp=ibsp+1

           do ilang=1,mbesslang

             do jlmn=1,lmn_size
               jl=indlmn(1,jlmn)
               jm=indlmn(2,jlmn)
               j0lmn=jlmn*(jlmn-1)/2
               do ilmn=1,jlmn
                 il=indlmn(1,ilmn)
                 im=indlmn(2,ilmn)
                 klmn=j0lmn+ilmn
                 kln=pawtab(itypat)%indklmn(2,klmn)

                 if (il==ilang-1.and.jl==ilang-1.and.im==jm) then

                   cpij=cprj_k(iat,ibsp)%cp(1,ilmn)*cprj_k(iat,ibsp)%cp(1,jlmn)
                   if (cplex==2) cpij=cpij+cprj_k(iat,ibsp)%cp(2,ilmn)*cprj_k(iat,ibsp)%cp(2,jlmn)
                   cpij=pawtab(itypat)%dltij(klmn)*cpij

                   dos%fractions(ikpt,iband,isppol,mbesslang*(iat-1)+ilang)=  &
&                   dos%fractions(ikpt,iband,isppol,mbesslang*(iat-1)+ilang) + &
&                   cpij*int1m2(kln,iat)
                   if (prtdosm==1) then
                     dos%fractions_m(ikpt,iband,isppol,mbesslang**2*(iat-1)+il**2+il+1+im)= &
&                     dos%fractions_m(ikpt,iband,isppol,mbesslang**2*(iat-1)+il**2+il+1+im) + &
&                     cpij*int1m2(kln,iat)
                   end if
                   if (fatbands_flag==1.or.prtdosm==2) then
                     dos%fractions_m(ikpt,iband,isppol,mbesslang**2*(iat-1)+il**2+il+1+im)= &
&                     dos%fractions_m(ikpt,iband,isppol,mbesslang**2*(iat-1)+il**2+il+1+im) + &
&                     cpij*int1(kln,iat)
                   end if
                   if (paw_dos_flag==1) then
                     dos%fractions_paw1(ikpt,iband,isppol,mbesslang*(iat-1)+ilang)=  &
&                     dos%fractions_paw1(ikpt,iband,isppol,mbesslang*(iat-1)+ilang) + &
&                     cpij*int1(kln,iat)
                     dos%fractions_pawt1(ikpt,iband,isppol,mbesslang*(iat-1)+ilang)=  &
&                     dos%fractions_pawt1(ikpt,iband,isppol,mbesslang*(iat-1)+ilang) + &
&                     cpij*int2(kln,iat)
                   end if

                 end if

               end do !ilmn
             end do   !jlmn

           end do ! ilang
         end do ! ispinor
       end do ! iband

     end do !iatom

     if (mkmem/=0) ibg = ibg + my_nspinor*nband_k
     call pawcprj_free(cprj_k)
     ABI_DATATYPE_DEALLOCATE(cprj_k)
   end do ! ikpt
 end do ! isppol

 ABI_DEALLOCATE(int1)
 ABI_DEALLOCATE(int2)
 ABI_DEALLOCATE(int1m2)

!Reduce data in case of parallelism
 call timab(48,1,tsec)
 call xmpi_sum(dos%fractions,spaceComm,ierr)
 if (prtdosm>=1.or.fatbands_flag==1) then
   call xmpi_sum(dos%fractions_m,spaceComm,ierr)
 end if
 if (paw_dos_flag==1) then
   call xmpi_sum(dos%fractions_paw1,spaceComm,ierr)
   call xmpi_sum(dos%fractions_pawt1,spaceComm,ierr)
 end if
 call timab(48,2,tsec)
 if (mpi_enreg%paral_spinor==1) then
   call xmpi_sum(dos%fractions,mpi_enreg%comm_spinor,ierr)
   if (prtdosm>=1.or.fatbands_flag==1) then
     call xmpi_sum(dos%fractions_m,mpi_enreg%comm_spinor,ierr)
   end if
   if (paw_dos_flag==1) then
     call xmpi_sum(dos%fractions_paw1, mpi_enreg%comm_spinor,ierr)
     call xmpi_sum(dos%fractions_pawt1, mpi_enreg%comm_spinor,ierr)
   end if
 end if

!Averaging: A quick hack for m-decomposed LDOS:
!BA: not valid in presence of spin-orbit coupling  !
 if (prtdosm==1.and.fatbands_flag==0) then
!  if pawfatbnd is activated, one think in the cubic harmonics basis
!  whereas prtdosm=1 is in the spherical harmonics basis.
!  the following trick is done in order to have everything
!  in the complex spherical basis (not useful for pawfatbnd if we want to
!  have e.g t2g and eg d-orbitals).
   do iat=1,dtset%natsph
     do il = 0, mbesslang-1
       do im = 1, il
         dos%fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1+im) = &
&         (dos%fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1+im) + &
&         dos%fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1-im))/2
         dos%fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1-im) = &
&         dos%fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1+im)
       end do
     end do
   end do !iatom
 end if

 DBG_EXIT("COLL")

end subroutine partial_dos_fractions_paw
!!***
