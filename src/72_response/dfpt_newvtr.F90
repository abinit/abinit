!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_newvtr
!! NAME
!! dfpt_newvtr
!!
!! FUNCTION
!! Compute new first-order trial potential by mixing new and old values.
!! First, compute preconditioned residual first-order potential.
!! Then, call one of the self-consistency drivers, and  update vtrial.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | isecur=level of security of the computation
!!   | mffmem=governs the number of FFT arrays which are fit in core memory
!!   |          it is either 1, in which case the array f_fftgr is used,
!!   |          or 0, in which case the array f_fftgr_disk is used
!!   | natom=number of atoms
!!   | nspden=number of spin-density components
!!   | paral_kgb=option for (kpt,g vectors,bands) parallelism
!!   | pawoptmix= - PAW only - 1 if the computed residuals include the PAW (rhoij) part
!!  etotal=the total energy obtained from the input vtrial
!!  ffttomix(nfft*(1-nfftmix/nfft))=Index of the points of the FFT (fine) grid on the grid used for mixing (coarse)
!!  initialized= if 0 the initialization of the RF run is not yet finished
!!   iscf=( <= 0 =>non-SCF), >0 => SCF)
!!    iscf =1 => determination of the largest eigenvalue of the SCF cycle
!!    iscf =2 => SCF cycle, simple mixing
!!    iscf =3 => SCF cycle, Anderson mixing
!!    iscf =4 => SCF cycle, Anderson mixing (order 2)
!!    iscf =5 => SCF cycle, CG based on the minimization of the energy
!!    iscf =7 => SCF cycle, Pulay mixing
!!  ispmix=1 if mixing is done in real space, 2 if mixing is done in reciprocal space
!!  istep= number of the step in the SCF cycle
!!  mixtofft(nfftmix*(1-nfftmix/nfft))=Index of the points of the FFT grid used for mixing (coarse) on the FFT (fine) grid
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nfftmix=dimension of FFT grid used to mix the densities (used in PAW only)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ngfftmix(18)=contain all needed information about 3D FFT, for the grid corresponding to nfftmix
!!  npawmix=-PAW only- number of spherical part elements to be mixed
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!                                         Use here rhoij residuals (and gradients)
!!  rhor(cplex*nfft,nspden)=array for 1st-order electron density
!!    in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vresid(cplex*nfft,nspden)=array for the residual of the potential
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  dbl_nnsclo=1 if nnsclo has to be doubled to secure the convergence.
!!
!! SIDE EFFECTS
!!  vtrial(cplex*nfft,nspden)= at input, it is the "in" trial potential that gave vresid=(v_out-v_in)
!!       at output, it is an updated "mixed" trial potential
!!  ==== if usepaw==1
!!    pawrhoij(natom)%nrhoijsel,rhoijselect,rhoijp= several arrays
!!                containing new values of rhoij (augmentation occupancies)
!!
!! NOTES
!!  In case of PAW calculations:
!!    Computations are done either on the fine FFT grid or the coarse grid (depending on dtset%pawmixdg)
!!    All variables (nfft,ngfft,mgfft) refer to the fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!    Developpers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!!  Subtility in PAW and non-collinear magnetism:
!!    Potentials are stored in (up-up,dn-dn,Re[up-dn],Im[up-dn]) format
!!    On-site occupancies (rhoij) are stored in (n,mx,my,mz)
!!    This is compatible provided that the mixing factors for n and m are identical
!!    and that the residual is not a combination of V_res and rhoij_res (pawoptmix=0).
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      ab7_mixing_copy_current_step,ab7_mixing_eval,ab7_mixing_eval_allocate
!!      ab7_mixing_eval_deallocate,fourdp,metric,moddiel,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dfpt_newvtr(cplex,dbl_nnsclo,dielar,dtset,etotal,ffttomix,&
&          initialized,iscf,ispmix,istep,mix,mixtofft,&
&          mpi_enreg,my_natom,nfft,nfftmix,ngfft,ngfftmix,npawmix,pawrhoij,&
&          qphon,rhor,rprimd,usepaw,vresid,vtrial)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_ab7_mixing
 use m_errors

 use m_pawrhoij, only : pawrhoij_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_newvtr'
 use interfaces_18_timing
 use interfaces_41_geometry
 use interfaces_53_ffts
 use interfaces_67_common
!End of the abilint section

 implicit none

!Arguments-------------------------------
!scalars
 integer,intent(in) :: cplex,initialized,iscf,ispmix,istep,my_natom,nfft
 integer,intent(in) :: nfftmix,npawmix,usepaw
 integer,intent(inout) :: dbl_nnsclo !vz_i
 real(dp),intent(in) :: etotal
 type(MPI_type),intent(in) :: mpi_enreg
 type(ab7_mixing_object), intent(inout) :: mix
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ffttomix(nfft*(1-nfftmix/nfft))
 integer,intent(in) :: mixtofft(nfftmix*(1-nfftmix/nfft)),ngfft(18)
 integer,intent(in) :: ngfftmix(18)
 real(dp),intent(in) :: dielar(7),qphon(3)
 real(dp), intent(in), target :: rhor(cplex*nfft,dtset%nspden)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: vresid(cplex*nfft,dtset%nspden)
 real(dp),intent(inout) :: vtrial(cplex*nfft,dtset%nspden)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*usepaw)

!Local variables-------------------------------
!scalars
 integer :: cplex_mix,cplex_rhoij,dplex,i_vresid1,i_vrespc1,iatom,ifft,indx
 integer :: irhoij,ispden,jfft,jrhoij,klmn,kmix,moved_atm_inside,nfftot,nselect
 integer :: mpicomm,errid
 logical :: mpi_summarize,reset
 real(dp) :: fact,mixfac,mixfac_eff,mixfacmag,ucvol
 character(len=500) :: msg
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),tsec(2)
 real(dp),allocatable :: rhoijrespc(:),rhoijtmp(:,:)
 real(dp),allocatable :: vresid0(:,:),vrespc(:,:),vreswk(:,:)
 real(dp), pointer :: vtrial0(:,:),vpaw(:)
 real(dp),allocatable :: vtrialg(:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(158,1,tsec)

!Compatibility tests
 if(usepaw==1) then
   if(dtset%nspden==4.and.dtset%pawoptmix==1) then
     msg='pawoptmix=1 is not compatible with nspden=4 !'
     MSG_ERROR(msg)
   end if
   if (my_natom>0) then
     if (pawrhoij(1)%cplex<cplex) then
       msg='pawrhoij()%cplex must be >=cplex !'
       MSG_ERROR(msg)
     end if
   end if
 end if

 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 cplex_mix=max(cplex,ispmix)
 if (usepaw==1.and.my_natom>0) cplex_rhoij=pawrhoij(1)%cplex

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 moved_atm_inside=0

!Select components of potential to be mixed
 ABI_ALLOCATE(vtrial0,(cplex_mix*nfftmix,dtset%nspden))
 ABI_ALLOCATE(vresid0,(cplex_mix*nfftmix,dtset%nspden))
 if (ispmix==1.and.nfft==nfftmix) then
   vtrial0=vtrial;vresid0=vresid
 else if (nfft==nfftmix) then
   do ispden=1,dtset%nspden
     call fourdp(cplex,vtrial0(:,ispden),vtrial(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
     call fourdp(cplex,vresid0(:,ispden),vresid(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
   end do
 else
   ABI_ALLOCATE(vtrialg,(2,nfft,dtset%nspden))
   ABI_ALLOCATE(vreswk,(2,nfft))
   do ispden=1,dtset%nspden
     fact=dielar(4);if (ispden>1) fact=dielar(7)
     call fourdp(cplex,vtrialg(:,:,ispden),vtrial(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
     call fourdp(cplex,vreswk,vresid(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
     do ifft=1,nfft
       if (ffttomix(ifft)>0) then
         jfft=2*ffttomix(ifft)
         vtrial0(jfft-1,ispden)=vtrialg(1,ifft,ispden)
         vtrial0(jfft  ,ispden)=vtrialg(2,ifft,ispden)
         vresid0(jfft-1,ispden)=vreswk(1,ifft)
         vresid0(jfft  ,ispden)=vreswk(2,ifft)
       else
         vtrialg(:,ifft,ispden)=vtrialg(:,ifft,ispden)+fact*vreswk(:,ifft)
       end if
     end do
   end do
   ABI_DEALLOCATE(vreswk)
 end if

!Precondition the potential residual:
!Use a model dielectric function preconditioning, or simple mixing
 ABI_ALLOCATE(vrespc,(cplex_mix*nfftmix,dtset%nspden))
 call moddiel(cplex_mix,dielar,mpi_enreg,nfftmix,ngfftmix,dtset%nspden,ispmix,0,&
& dtset%paral_kgb,qphon,rprimd,vresid0,vrespc)

!PAW only : precondition the rhoij quantities (augmentation occupancies) residuals.
!Use a simple preconditionning with the same mixing factor
!as the model dielectric function.
 if (usepaw==1.and.my_natom>0) then
   ABI_ALLOCATE(rhoijrespc,(npawmix))
   mixfac=dielar(4);mixfacmag=abs(dielar(7))
   if (cplex_rhoij==1) then
     indx=0
     do iatom=1,my_natom
       do ispden=1,pawrhoij(iatom)%nspden
         mixfac_eff=mixfac;if (ispden>1) mixfac_eff=mixfacmag
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           indx=indx+1;klmn=pawrhoij(iatom)%kpawmix(kmix)
           rhoijrespc(indx)=mixfac_eff*pawrhoij(iatom)%rhoijres(klmn,ispden)
         end do
       end do
     end do
   else
     indx=-1
     do iatom=1,my_natom
       do ispden=1,pawrhoij(iatom)%nspden
         mixfac_eff=mixfac;if (ispden>1) mixfac_eff=mixfacmag
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           indx=indx+2;klmn=2*pawrhoij(iatom)%kpawmix(kmix)-1
           rhoijrespc(indx:indx+1)=mixfac_eff*pawrhoij(iatom)%rhoijres(klmn:klmn+1,ispden)
         end do
       end do
     end do
   end if
 end if

!------Compute new vtrial

 i_vresid1=mix%i_vresid(1)
 i_vrespc1=mix%i_vrespc(1)

!Initialise working arrays for the mixing object.
 call ab7_mixing_eval_allocate(mix, istep)

!Copy current step arrays.
 call ab7_mixing_copy_current_step(mix, vresid0, errid, msg, arr_respc = vrespc)

 if (errid /= AB7_NO_ERROR) then
   MSG_ERROR(msg)
 end if

 ABI_DEALLOCATE(vrespc)
 ABI_DEALLOCATE(vresid0)

!PAW: either use the array f_paw or the array f_paw_disk
 ABI_ALLOCATE(vpaw,(npawmix*usepaw))
 if (usepaw==1.and.my_natom>0) then
   dplex=cplex_rhoij-1
   indx=-dplex
   do iatom=1,my_natom
     do ispden=1,pawrhoij(iatom)%nspden
       ABI_ALLOCATE(rhoijtmp,(cplex_rhoij*pawrhoij(iatom)%lmn2_size,1))
       rhoijtmp=zero
       jrhoij=1
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=cplex_rhoij*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
         rhoijtmp(klmn:klmn+dplex,1)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
         jrhoij=jrhoij+cplex_rhoij
       end do
       do kmix=1,pawrhoij(iatom)%lmnmix_sz
         indx=indx+cplex_rhoij;klmn=cplex_rhoij*pawrhoij(iatom)%kpawmix(kmix)-dplex
         vpaw(indx:indx+dplex)=rhoijtmp(klmn:klmn+dplex,1)-pawrhoij(iatom)%rhoijres(klmn:klmn+dplex,ispden)
         mix%f_paw(indx:indx+dplex,i_vresid1)=pawrhoij(iatom)%rhoijres(klmn:klmn+dplex,ispden)
         mix%f_paw(indx:indx+dplex,i_vrespc1)=rhoijrespc(indx:indx+dplex)
       end do
       ABI_DEALLOCATE(rhoijtmp)
     end do
   end do
 end if

!Unlike for GS, no need to modify the mean of vtrial

 mpicomm=0;mpi_summarize=.false.
 reset=.false.;if (initialized==0) reset=.true.
 call ab7_mixing_eval(mix, vtrial0, istep, nfftot, ucvol, &
& mpicomm, mpi_summarize, errid, msg, &
& reset = reset, isecur = dtset%isecur, &
& pawopt = dtset%pawoptmix, response = 1, pawarr = vpaw, &
& etotal = etotal, potden = rhor, comm_atom=mpi_enreg%comm_atom)

 if (errid == AB7_ERROR_MIXING_INC_NNSLOOP) then
   dbl_nnsclo = 1
 else if (errid /= AB7_NO_ERROR) then
   ! MG FIXME, Why this?
   ! One should propagate the error so that we can handle it
   ! in the caller!
   MSG_ERROR(msg)
 end if

!Do here the mixing of the potential
 if(iscf==2 .or. iscf==3 .or. iscf==7)then
!  PAW: restore rhoij from compact storage
   if (usepaw==1.and.my_natom>0) then
     dplex=cplex_rhoij-1
     indx=-dplex
     do iatom=1,my_natom
       ABI_ALLOCATE(rhoijtmp,(cplex_rhoij*pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
       rhoijtmp=zero
       if (pawrhoij(iatom)%lmnmix_sz<pawrhoij(iatom)%lmn2_size) then
         do ispden=1,pawrhoij(iatom)%nspden
           jrhoij=1
           do irhoij=1,pawrhoij(iatom)%nrhoijsel
             klmn=cplex_rhoij*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
             rhoijtmp(klmn:klmn+dplex,ispden)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
             jrhoij=jrhoij+cplex_rhoij
           end do
         end do
       end if
       do ispden=1,pawrhoij(iatom)%nspden
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           indx=indx+cplex_rhoij;klmn=cplex_rhoij*pawrhoij(iatom)%kpawmix(kmix)-dplex
           rhoijtmp(klmn:klmn+dplex,ispden)=vpaw(indx:indx+dplex)
         end do
       end do
       nselect=0
       do klmn=1,pawrhoij(iatom)%lmn2_size
         if (any(abs(rhoijtmp(cplex_rhoij*klmn-dplex:cplex_rhoij*klmn,:))>tol10)) then
           nselect=nselect+1
           pawrhoij(iatom)%rhoijselect(nselect)=klmn
           do ispden=1,pawrhoij(iatom)%nspden
             pawrhoij(iatom)%rhoijp(cplex_rhoij*nselect-dplex:cplex_rhoij*nselect,ispden)=&
&             rhoijtmp(cplex_rhoij*klmn-dplex:cplex_rhoij*klmn,ispden)
           end do
         end if
       end do
       pawrhoij(iatom)%nrhoijsel=nselect
       ABI_DEALLOCATE(rhoijtmp)
     end do
   end if

 else if(iscf==5 .or. iscf==6)then
   if(ispmix/=1) then
     msg = ' Mixing on reciprocal space not allowed with iscf=5 or 6.'
     MSG_ERROR(msg)
   end if
!  PAW: apply a simple mixing to rhoij (this is temporary)
   if (usepaw==1.and.my_natom>0) then
     indx=1-cplex_rhoij
     do iatom=1,my_natom
       ABI_ALLOCATE(rhoijtmp,(cplex_rhoij*pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
       rhoijtmp=zero
       if (pawrhoij(iatom)%lmnmix_sz<pawrhoij(iatom)%lmn2_size) then
         do ispden=1,pawrhoij(iatom)%nspden
           do kmix=1,pawrhoij(iatom)%lmnmix_sz
             indx=indx+cplex_rhoij;klmn=cplex_rhoij*pawrhoij(iatom)%kpawmix(kmix)-dplex
             rhoijtmp(klmn:klmn+dplex,ispden)=rhoijrespc(indx:indx+dplex) &
&             -pawrhoij(iatom)%rhoijres(klmn:klmn+dplex,ispden)
           end do
         end do
       end if
       do ispden=1,pawrhoij(iatom)%nspden
         jrhoij=1
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=cplex_rhoij*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
           rhoijtmp(klmn:klmn+dplex,ispden)=rhoijtmp(klmn:klmn+dplex,ispden) &
&           +pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
           jrhoij=jrhoij+cplex_rhoij
         end do
       end do
       nselect=0
       do klmn=1,pawrhoij(iatom)%lmn2_size
         if (any(abs(rhoijtmp(cplex_rhoij*klmn-dplex:cplex_rhoij*klmn,:))>tol10)) then
           nselect=nselect+1
           pawrhoij(iatom)%rhoijselect(nselect)=klmn
           do ispden=1,pawrhoij(iatom)%nspden
             pawrhoij(iatom)%rhoijp(cplex_rhoij*nselect-dplex:cplex_rhoij*nselect,ispden)=&
&             rhoijtmp(cplex_rhoij*klmn-dplex:cplex_rhoij*klmn,ispden)
           end do
         end if
       end do
       pawrhoij(iatom)%nrhoijsel=nselect
       ABI_DEALLOCATE(rhoijtmp)
     end do
   end if
 end if

 ABI_DEALLOCATE(vpaw)
 if (usepaw==1.and.my_natom>0)  then
   ABI_DEALLOCATE(rhoijrespc)
 end if

!Eventually write the data on disk and deallocate f_fftgr_disk
 call ab7_mixing_eval_deallocate(mix)

!Restore potential
 if (ispmix==1.and.nfft==nfftmix) then
   vtrial=vtrial0
 else if (nfft==nfftmix) then
   do ispden=1,dtset%nspden
     call fourdp(cplex,vtrial0(:,ispden),vtrial(:,ispden),+1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
   end do
 else
   do ispden=1,dtset%nspden
     do ifft=1,nfftmix
       jfft=mixtofft(ifft)
       vtrialg(1,jfft,ispden)=vtrial0(2*ifft-1,ispden)
       vtrialg(2,jfft,ispden)=vtrial0(2*ifft  ,ispden)
     end do
     call fourdp(cplex,vtrialg(:,:,ispden),vtrial(:,ispden),+1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
   end do
   ABI_DEALLOCATE(vtrialg)
 end if
 ABI_DEALLOCATE(vtrial0)

 call timab(158,2,tsec)

 DBG_ENTER("COLL")

end subroutine dfpt_newvtr
!!***
