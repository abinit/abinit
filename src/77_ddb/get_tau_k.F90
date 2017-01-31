!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_tau_k
!! NAME
!!  get_tau_k
!!
!! FUNCTION
!!  Calculate the k-dependent relaxation time due to EPC. Impelementation based
!!  on derivation from Grmvall's book or OD Restrepo's paper (PRB 94 212103 (2009))
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2016 ABINIT group (BXu)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Cryst<crystal_t>=Info on the unit cell and on its symmetries.
!!  Ifc<ifc_type>=Object containing the interatomic force constants.
!!  elph_ds = elphon datastructure with data and dimensions
!!  eigenGS = Ground State eigenvalues
!!  max_occ = maximal occupancy for a band
!!
!! OUTPUT
!!  tau_k(nsppol,nkptirr,nband)=mode relaxation time due to electron phonono coupling
!!  rate_e(nene)= scattering rate due to electron phonono coupling vs. energy
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      dgemm,ebands_prtbltztrp_tau_out,ebands_update_occ,ep_el_weights
!!      ep_ph_weights,ftgam,ftgam_init,gam_mult_displ,ifc_fourq,matrginv
!!      mkqptequiv,phdispl_cart2red,spline,splint,wrtout,xmpi_sum,zgemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine get_tau_k(Cryst,ifc,Bst,elph_ds,elph_tr_ds,eigenGS,max_occ)
    
 use defs_basis
 use defs_elphon
 use defs_datatypes
 use m_kptrank
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_splines
 use m_ifc
 use m_ebands

 use m_io_tools,   only : open_file
 use m_geometry,   only : phdispl_cart2red
 use m_dynmat,     only : ftgam_init, ftgam
 use m_crystal,    only : crystal_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_tau_k'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_77_ddb, except_this_one => get_tau_k
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(crystal_t),intent(in) :: Cryst
 type(ifc_type),intent(in) :: ifc
 type(ebands_t),intent(inout)   :: Bst
 type(elph_type),intent(inout) :: elph_ds
 type(elph_tr_type), intent(inout) :: elph_tr_ds
 real(dp),intent(in) :: max_occ
 real(dp),intent(in) :: eigenGS(elph_ds%nband,elph_ds%k_phon%nkpt,elph_ds%nsppol)

!Local variables-------------------------------
!scalars
 character(len=500) :: message
 character(len=fnlen) :: fname
 integer :: ntemper,nsppol,nbranch,nband,natom
 integer :: nkpt,nqpt,nkptirr,nqptirr,new_nkptirr
 integer :: isppol,iFSkpt,iFSqpt,iqpt,iqpt_fullbz,imqpt_fullbz,ikpt_kpq,ikpt_kmq
 integer :: iband,jband,jpband,jbeff,ibranch,jbranch,itemp
 integer :: irec,ierr,nrpt,ik_this_proc
 integer :: unit_tau,unit_invtau
 integer :: nene,nene_all,iene,iene_fine,unit_taue,unit_mfp
 integer :: icomp,jcomp,itensor
 integer :: ikpt_irr,iomega,unit_cond,unit_therm,unit_sbk
 integer :: nskip,nspline
 real(dp) :: occ_omega,occ_e
 real(dp) :: xx,Temp,therm_factor
 real(dp) :: factor,dfermide
 real(dp) :: e_k,chu_tau,rate_e,mfp_e
 real(dp) :: ene,enemin,enemax,deltaene
 real(dp) :: omega,omega_min,omega_max,domega
 real(dp) :: diagerr
 real(dp) :: chu_mfp,chu_cond,chu_cth,chu_sbk,femto
 real(dp) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: eigval(elph_ds%nbranch),eigval2(elph_ds%nbranch)
 real(dp) :: imeigval(elph_ds%nbranch)
 real(dp) :: tmp_wtkpq, tmp_wtkmq, tol_wtk
 real(dp) :: yp1,ypn
!arrays
 integer,allocatable :: FSfullpktofull(:,:),mqtofull(:)
 integer,allocatable :: kpttokpt(:,:,:)
 real(dp) :: cond_inv(3,3)
 real(dp),allocatable :: fermie(:)
 real(dp),allocatable :: tmp_eigenGS(:,:,:)
 real(dp),allocatable :: tmp_gkk_qpt(:,:,:),tmp_gkk_rpt(:,:,:),tmp_gkk_kpt(:,:)
 real(dp),allocatable :: tmp_gkk_kpt2(:,:,:), gkk_kpt(:,:,:)
 real(dp),allocatable :: tau_k(:,:,:,:),inv_tau_k(:,:,:,:),tmp_tau_k(:,:,:,:)
 real(dp),allocatable :: phfrq(:,:),pheigvec(:,:)
 real(dp),allocatable :: displ(:,:,:,:)
 real(dp),allocatable :: a2f_2d(:),a2f_2d2(:)
 real(dp),allocatable :: tmp_wtk(:,:,:,:),tmp2_wtk(:),tmp_wtk1(:),tmp_wtk2(:)
 real(dp),allocatable :: ene_pt(:),ene_ptfine(:),ff2(:)
 real(dp),allocatable :: wtq(:,:,:),tmp_wtq(:,:,:),tmp2_wtq(:,:)
 real(dp),allocatable :: dos_e(:,:)
 real(dp),allocatable :: coskr1(:,:),sinkr1(:,:)
 real(dp),allocatable :: coskr2(:,:),sinkr2(:,:)
 real(dp),allocatable :: cond_e(:,:,:,:),cond(:,:,:,:),sbk(:,:,:,:),seebeck(:,:,:,:),cth(:,:,:,:)
 
! *************************************************************************
 
 write(std_out,*) 'get_tau_k : enter '

 nrpt = ifc%nrpt
 natom = cryst%natom

 nsppol   = elph_ds%nsppol
 nbranch  = elph_ds%nbranch
 nband    = elph_ds%ngkkband
 nkpt     = elph_ds%k_phon%nkpt
 nqpt     = elph_ds%nqpt_full
 nkptirr  = elph_ds%k_phon%nkptirr
 new_nkptirr  = elph_ds%k_phon%new_nkptirr
 nqptirr  = elph_ds%nqptirred
 ntemper  = elph_ds%ntemper
 nene = 2*elph_ds%na2f-1 ! only need e_k +- omega_max range, take deltaene=delta_oemga

 chu_tau  = 2.4188843265*1.0d-17
 chu_mfp  = 5.291772*1.0d-11
 chu_cond = 4.59988159904764*1.0d6
 chu_cth  = 1.078637439971599*1.0d4
 chu_sbk  = 8.617343101*1.0d-5
 femto    = 1.0d-15

 tol_wtk = tol7/nkptirr/nband

 ABI_ALLOCATE(fermie ,(ntemper))
 ABI_ALLOCATE(tmp_gkk_qpt ,(2,nbranch**2,nqpt))
 ABI_ALLOCATE(tmp_gkk_rpt ,(2,nbranch**2,nrpt))
 ABI_ALLOCATE(tmp_gkk_kpt ,(2,nbranch**2))
 ABI_ALLOCATE(tmp_gkk_kpt2 ,(2,nbranch,nbranch))
 ABI_ALLOCATE(gkk_kpt ,(2,nbranch,nbranch))
 ABI_ALLOCATE(a2f_2d, (nene))
 ABI_ALLOCATE(a2f_2d2, (nene))
 ABI_ALLOCATE(inv_tau_k, (ntemper,nsppol,nkpt,nband))
 ABI_ALLOCATE(tau_k, (ntemper,nsppol,nkpt,nband))
 ABI_ALLOCATE(tmp_tau_k ,(ntemper,nsppol,new_nkptirr,nband))

 if (elph_ds%gkqwrite == 0) then
   call wrtout(std_out,' get_tau_k : keeping gkq matrices in memory','COLL')
 else if (elph_ds%gkqwrite == 1) then
   fname=trim(elph_ds%elph_base_name) // '_GKKQ'
   write (message,'(2a)')' get_tau_k : reading gkq matrices from file ',trim(fname)
   call wrtout(std_out,message,'COLL')
 else
   write (message,'(a,i0)')' Wrong value for gkqwrite = ',elph_ds%gkqwrite
   MSG_BUG(message)
 end if

!=========================================================    
!Get equivalence between a kpt_phon pair and a qpt in qpt_full
!only works if the qpt grid is complete (identical to
!the kpt one, with a basic shift of (0,0,0)
!=========================================================    

!mapping of k + q onto k' for k and k' in full BZ
!for dense k grid
 ABI_ALLOCATE(FSfullpktofull,(nkpt,nkpt))
 ABI_ALLOCATE(mqtofull,(nkpt))

!kpttokpt(itim,isym,iqpt) = kpoint index which transforms to ikpt under isym and with time reversal itim.
 ABI_ALLOCATE(kpttokpt,(2,Cryst%nsym,nkpt))

 call wrtout(std_out,'get_tau_k: calling mkqptequiv to set up the FS kpoint set',"COLL")

 call mkqptequiv (FSfullpktofull,Cryst,elph_ds%k_phon%kpt,nkpt,nkpt,kpttokpt,elph_ds%k_phon%kpt,mqtofull)

!=========================================================    
!=========================================================    

 omega_max       = elph_ds%omega_max
 omega_min       = elph_ds%omega_min
 domega          = elph_ds%domega
 enemax = maxval(eigenGS(elph_ds%maxFSband,:,:))
 enemin = minval(eigenGS(elph_ds%minFSband,:,:))

 if (enemin < (elph_ds%fermie-0.2)) then
   enemin = elph_ds%fermie-0.2
 end if
 if (enemax > (elph_ds%fermie+0.2)) then
   enemax = elph_ds%fermie+0.2
 end if

 nspline = elph_ds%ep_nspline
 nene_all = INT((enemax-enemin+domega)/(nspline*domega)) + 1
 deltaene = domega
 write(std_out,*) 'E_min= ',enemin, 'E_max= ',enemax
 write(std_out,*) 'Number of energy points= ',nene_all
 write(std_out,'(a,I8)') 'scale factor for spline interpolation in RTA = ', elph_ds%ep_nspline
 write(std_out,*) 'delta_ene before spline interpolation= ',deltaene*nspline
 write(std_out,*) 'delta_ene after spline interpolation= ',deltaene
 write(std_out,*) 'Omega_min= ',omega_min, 'Omega_max= ',omega_max
 write(std_out,*) 'Number of phonon points= ',elph_ds%na2f
 write(std_out,*) 'delta_omega= ',domega
 write(std_out,*) 'number of bands= ', elph_ds%nband, nband

 ABI_ALLOCATE(tmp_wtk,(nband,nkpt,nsppol,nene_all))
 ABI_ALLOCATE(tmp2_wtk,(nene_all))
 ABI_ALLOCATE(ff2,(nene_all))
 ABI_ALLOCATE(ene_pt,(nene_all))
 ABI_ALLOCATE(ene_ptfine,(nene_all*nspline))
 ABI_ALLOCATE(tmp_wtk1,(nene_all*nspline))
 ABI_ALLOCATE(tmp_wtk2,(nene_all*nspline))
 ABI_ALLOCATE(dos_e,(nsppol,nene_all))

!Get energy points for spline interpolation
 do iene = 1, nene_all
   ene_pt(iene) = enemin + (iene-1)*nspline*deltaene
 end do

 do iene = 1, nene_all*nspline
   ene_ptfine(iene) = enemin + (iene-1)*deltaene
 end do

 ABI_ALLOCATE(tmp_wtq,(elph_ds%nbranch, elph_ds%k_phon%nkpt, elph_ds%na2f+1))
 ABI_ALLOCATE(wtq,(elph_ds%nbranch, elph_ds%k_phon%nkpt, elph_ds%na2f))
 ABI_ALLOCATE(tmp2_wtq,(elph_ds%nbranch, elph_ds%na2f))

!phonon
 ABI_ALLOCATE(phfrq,(nbranch, nkptirr))
 ABI_ALLOCATE(displ,(2, nbranch, nbranch, nkptirr))
 ABI_ALLOCATE(pheigvec,(2*nbranch*nbranch, nkptirr))

 do iFSqpt = 1, nkptirr
   call ifc_fourq(ifc,cryst,elph_ds%k_phon%kptirr(:,iFSqpt),phfrq(:,iFSqpt),displ(:,:,:,iFSqpt),out_eigvec=pheigvec(:,iFSqpt))
 end do

 omega_min = omega_min - domega

!bxu, obtain wtq for the q_fine, then condense to q_phon
 call ep_ph_weights(phfrq,elph_ds%a2fsmear,omega_min,omega_max,elph_ds%na2f+1,Cryst%gprimd,elph_ds%kptrlatt, &
& elph_ds%nbranch,elph_ds%telphint,elph_ds%k_phon,tmp_wtq)
 omega_min = omega_min + domega

 do iomega = 1, elph_ds%na2f
   wtq(:,:,iomega) = tmp_wtq(:,:,iomega+1)
   !write(1005,*) omega_min+(iomega-1)*domega, sum(tmp_wtq(:,:,iomega+1))/nkpt
 end do
 ABI_DEALLOCATE(tmp_wtq)

! electron 
 tmp_wtk =zero
 dos_e = zero
 call ep_el_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS(elph_ds%minFSband:elph_ds%minFSband+nband-1,:,:), &
& elph_ds%elphsmear, &
& enemin, enemax, nene_all, Cryst%gprimd, elph_ds%k_phon%irredtoGS, elph_ds%kptrlatt, max_occ, &
& 1, nband, elph_ds%nFSband, nsppol, elph_ds%telphint, elph_ds%k_phon, tmp_wtk)
!& elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, nsppol, elph_ds%telphint, elph_ds%k_phon, tmp_wtk)

 do isppol = 1, nsppol
   do iene = 1, nene_all
     dos_e(isppol,iene) = sum(tmp_wtk(:,:,isppol,iene))/nkpt
   end do
 end do

 ABI_ALLOCATE(coskr1, (nqpt,nrpt))
 ABI_ALLOCATE(sinkr1, (nqpt,nrpt))
 call ftgam_init(ifc%gprim, nqpt, nrpt, elph_ds%k_phon%kpt, Ifc%rpt, coskr1, sinkr1)
 ABI_ALLOCATE(coskr2, (nkptirr,nrpt))
 ABI_ALLOCATE(sinkr2, (nkptirr,nrpt))
 call ftgam_init(ifc%gprim, nkptirr, nrpt, elph_ds%k_phon%kpt, Ifc%rpt, coskr2, sinkr2)

!get fermie for itemp
 fermie = elph_ds%fermie
 do itemp=1,ntemper  ! runs over termperature in K
   Temp=elph_ds%tempermin+elph_ds%temperinc*dble(itemp)

   Bst%occopt = 3
   Bst%tsmear = Temp*kb_HaK
   call ebands_update_occ(Bst,-99.99_dp)
   write(message,'(a,f12.6,a,E20.12)')'At T=',Temp,' Fermi level is:',Bst%fermie
   call wrtout(std_out,message,'COLL')
   
   if (abs(elph_ds%fermie) < tol10) then
     fermie(itemp) = Bst%fermie
   end if
 end do

 inv_tau_k = zero
!get a2f_2d = \sum_{q,nbranch,jband'} |gkk|^2*\delta(\epsilon_{k'j'}-\epsilon')*\delta(\omega_q-\omega)
 do isppol=1,nsppol
   write (std_out,*) '##############################################'
   write (std_out,*) 'get_tau_k : Treating spin polarization ', isppol
   write (std_out,*) '##############################################'

!   do iFSkpt =1,nkpt
   do ik_this_proc =1,elph_ds%k_phon%my_nkpt
     iFSkpt = elph_ds%k_phon%my_ikpt(ik_this_proc)
     write (std_out,*) 'get_tau_k : working on kpt # ', iFSkpt, '/', nkpt
     do jband = 1, nband
!          write(*,*)'i am here 1 ', isppol,iFSkpt,jband
       a2f_2d = zero
       a2f_2d2 = zero

!sum from here
       nskip = 0
       do jpband = 1, nband
         jbeff = jpband+(jband-1)*nband

         if (elph_ds%gkqwrite == 0) then
           tmp_gkk_qpt(:,:,:) = elph_ds%gkk_qpt(:,jbeff,:,ik_this_proc,isppol,:)
         else if (elph_ds%gkqwrite == 1) then
           irec = (ik_this_proc-1)*elph_ds%k_phon%my_nkpt + iqpt
           if (iFSkpt == 1) then
             write (std_out,*) ' get_tau_k  read record ', irec
           end if
           read (elph_ds%unitgkq,REC=irec) tmp_gkk_qpt(:,:,iqpt_fullbz)
         end if

!FT to real space
         call ftgam(Ifc%wghatm,tmp_gkk_qpt,tmp_gkk_rpt,natom,nqpt,nrpt,1,coskr1,sinkr1)

!sum over irred q over k_phon, with corresponding weights
         do iFSqpt = 1, nkptirr
           iqpt_fullbz = elph_ds%k_phon%irredtoGS(iFSqpt)
           ikpt_kpq = FSfullpktofull(iFSkpt,iqpt_fullbz)

           imqpt_fullbz = mqtofull(iqpt_fullbz)
           ikpt_kmq = FSfullpktofull(iFSkpt,imqpt_fullbz)

!Do FT from real-space gamma grid to 1 kpt in k_phon%new_kptirr
           call ftgam(Ifc%wghatm,tmp_gkk_kpt,tmp_gkk_rpt,natom,1,nrpt,0,coskr2(iqpt_fullbz,:),sinkr2(iqpt_fullbz,:))
!tmp_gkk_kpt(:,:)=tmp_gkk_qpt(:,:,iFSqpt)

!if ep_scalprod==0 we have to dot in the displacement vectors here
           if (elph_ds%ep_scalprod==0) then

             call phdispl_cart2red(natom,Cryst%gprimd,displ(:,:,:,iFSqpt),displ_red)

             tmp_gkk_kpt2 = reshape (tmp_gkk_kpt(:,:), (/2,nbranch,nbranch/))
             call gam_mult_displ(nbranch, displ_red, tmp_gkk_kpt2, gkk_kpt)

             do jbranch=1,nbranch
               eigval(jbranch) = gkk_kpt(1, jbranch, jbranch)
               imeigval(jbranch) = gkk_kpt(2, jbranch, jbranch)

               if (abs(imeigval(jbranch)) > tol10) then
                 write (message,'(a,i0,a,es16.8)')" real values  branch = ",jbranch,' eigval = ',eigval(jbranch)
                 MSG_WARNING(message)
                 write (message,'(a,i0,a,es16.8)')" imaginary values  branch = ",jbranch,' imeigval = ',imeigval(jbranch)
                 MSG_WARNING(message)
               end if

             end do

!            if ep_scalprod==1 we have to diagonalize the matrix we interpolated.
           else if (elph_ds%ep_scalprod == 1) then

!            MJV NOTE : gam_now is being recast as a (3*natom)**2 matrix here
             call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, cone, tmp_gkk_kpt, 3*natom,&
&             pheigvec(:,iFSqpt), 3*natom, czero, tmp_gkk_kpt2, 3*natom)

             call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, cone, pheigvec(:,iFSqpt), 3*natom,&
&             tmp_gkk_kpt2, 3*natom, czero, gkk_kpt, 3*natom)

             diagerr = zero
             do ibranch=1,nbranch
               eigval(ibranch) = gkk_kpt(1,ibranch,ibranch)
               do jbranch=1,ibranch-1
                 diagerr = diagerr + abs(gkk_kpt(1,jbranch,ibranch))
               end do
               do jbranch=ibranch+1,nbranch
                 diagerr = diagerr + abs(gkk_kpt(1,jbranch,ibranch))
               end do
             end do

             if (diagerr > tol12) then
               write(message,'(a,es15.8)') 'get_tau_k: residual in diagonalization of gamma with phon eigenvectors: ', diagerr
               MSG_WARNING(message)
             end if

           else
             write (message,'(a,i0)')' Wrong value for ep_scalprod = ',elph_ds%ep_scalprod
             MSG_BUG(message)
           end if ! end ep_scalprod if

!For k'=k-q
!Do FT from real-space gamma grid to 1 kpt in k_phon%new_kptirr
           call ftgam(Ifc%wghatm,tmp_gkk_kpt,tmp_gkk_rpt,natom,1,nrpt,0,coskr2(imqpt_fullbz,:),sinkr2(imqpt_fullbz,:))
!tmp_gkk_kpt(:,:)=tmp_gkk_qpt(:,:,iFSqpt)

!if ep_scalprod==0 we have to dot in the displacement vectors here
           if (elph_ds%ep_scalprod==0) then

             call phdispl_cart2red(natom,Cryst%gprimd,displ(:,:,:,iFSqpt),displ_red)

             tmp_gkk_kpt2 = reshape (tmp_gkk_kpt(:,:), (/2,nbranch,nbranch/))
             call gam_mult_displ(nbranch, displ_red, tmp_gkk_kpt2, gkk_kpt)

             do jbranch=1,nbranch
               eigval2(jbranch) = gkk_kpt(1, jbranch, jbranch)
               imeigval(jbranch) = gkk_kpt(2, jbranch, jbranch)

               if (abs(imeigval(jbranch)) > tol10) then
                 write (message,'(a,i0,a,es16.8)')" real values  branch = ",jbranch,' eigval = ',eigval2(jbranch)
                 MSG_WARNING(message)
                 write (message,'(a,i0,a,es16.8)')" imaginary values  branch = ",jbranch,' imeigval = ',imeigval(jbranch)
                 MSG_WARNING(message)
               end if

             end do

!            if ep_scalprod==1 we have to diagonalize the matrix we interpolated.
           else if (elph_ds%ep_scalprod == 1) then

!            MJV NOTE : gam_now is being recast as a (3*natom)**2 matrix here
             call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, cone, tmp_gkk_kpt, 3*natom,&
&             pheigvec(:,iFSqpt), 3*natom, czero, tmp_gkk_kpt2, 3*natom)

             call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, cone, pheigvec(:,iFSqpt), 3*natom,&
&             tmp_gkk_kpt2, 3*natom, czero, gkk_kpt, 3*natom)

             diagerr = zero
             do ibranch=1,nbranch
               eigval2(ibranch) = gkk_kpt(1,ibranch,ibranch)
               do jbranch=1,ibranch-1
                 diagerr = diagerr + abs(gkk_kpt(1,jbranch,ibranch))
               end do
               do jbranch=ibranch+1,nbranch
                 diagerr = diagerr + abs(gkk_kpt(1,jbranch,ibranch))
               end do
             end do

             if (diagerr > tol12) then
               write(message,'(a,es15.8)') 'get_tau_k: residual in diagonalization of gamma with phon eigenvectors: ', diagerr
               MSG_WARNING(message)
             end if

           else
             write (message,'(a,i0)')' Wrong value for ep_scalprod = ',elph_ds%ep_scalprod
             MSG_BUG(message)
           end if ! end ep_scalprod if

           tmp2_wtk(:) = tmp_wtk(jpband,ikpt_kpq,isppol,:)
           yp1 = (tmp2_wtk(2)-tmp2_wtk(1))/nspline/deltaene
           ypn = (tmp2_wtk(nene_all)-tmp2_wtk(nene_all-1))/nspline/deltaene
           call spline(ene_pt,tmp2_wtk,nene_all,yp1,ypn,ff2)
           call splint(nene_all,ene_pt,tmp2_wtk,ff2,nene_all*nspline,ene_ptfine,tmp_wtk1)

           tmp2_wtk(:) = tmp_wtk(jpband,ikpt_kmq,isppol,:)
           yp1 = (tmp2_wtk(2)-tmp2_wtk(1))/nspline/deltaene
           ypn = (tmp2_wtk(nene_all)-tmp2_wtk(nene_all-1))/nspline/deltaene
           call spline(ene_pt,tmp2_wtk,nene_all,yp1,ypn,ff2)
           call splint(nene_all,ene_pt,tmp2_wtk,ff2,nene_all*nspline,ene_ptfine,tmp_wtk2)

           tmp2_wtq(:,:) = wtq(:,iFSqpt,:)
           do iene=1,nene 
             e_k = eigenGS(elph_ds%minFSband+jband-1,iFSkpt,isppol)
             ene = e_k - omega_max + (iene-1)*deltaene
             if (ene<enemin .or. ene>enemax) cycle
             iene_fine = NINT((ene-enemin+deltaene)/deltaene)
             tmp_wtkpq = tmp_wtk1(iene_fine) * elph_ds%k_phon%wtkirr(iFSqpt)
             tmp_wtkmq = tmp_wtk2(iene_fine) * elph_ds%k_phon%wtkirr(iFSqpt)

             if (tmp_wtkpq+tmp_wtkmq < tol_wtk ) then
               nskip = nskip +1 
               cycle
             end if

             do ibranch = 1, nbranch
               if (abs(phfrq(ibranch,iFSqpt)) < tol7) cycle

               if (ene > e_k) then
                 omega = ene - e_k
                 if (abs(omega) < tol7 .or. abs(omega) > omega_max) cycle
                 iomega = NINT((omega-omega_min+domega)/domega)

                 a2f_2d(iene) = a2f_2d(iene) +&
&                 eigval(ibranch)/phfrq(ibranch,iFSqpt)*&
&                 tmp_wtkpq * tmp2_wtq(ibranch,iomega)
               end if

               if (ene < e_k) then
                 omega = e_k - ene
                 if (abs(omega) < tol7 .or. abs(omega) > omega_max) cycle
                 iomega = NINT((omega-omega_min+domega)/domega)

                 a2f_2d2(iene) = a2f_2d2(iene) +&
&                 eigval(ibranch)/phfrq(ibranch,iFSqpt)*&
&                 tmp_wtkmq * tmp2_wtq(ibranch,iomega)
               end if

             end do ! ibranch 3
           end do ! nene  800 
         end do ! kptirr 216
       end do ! j' band 3
!      print *, ' skipped ',  nskip, ' energy points out of ', nene*nband*nkptirr

! get inv_tau_k
       do itemp=1,ntemper  ! runs over termperature in K
         Temp=elph_ds%tempermin+elph_ds%temperinc*dble(itemp)
         do iene=1,nene 
           e_k = eigenGS(elph_ds%minFSband+jband-1,iFSkpt,isppol)
           ene = e_k - omega_max + (iene-1)*deltaene
           if (ene<enemin .or. ene>enemax) cycle

           xx=(ene-fermie(itemp))/(kb_HaK*Temp)
           occ_e=1.0_dp/(exp(xx)+1.0_dp)
           if (ene > e_k .and. (ene-e_k) .le. omega_max) then
             omega = ene - e_k
             if (abs(omega) < tol7) cycle
             xx = omega/(kb_HaK*Temp)
             occ_omega=1.0_dp/(exp(xx)-1.0_dp)

             therm_factor = occ_e + occ_omega

             inv_tau_k(itemp,isppol,iFSkpt,jband) = inv_tau_k(itemp,isppol,iFSkpt,jband) +&
             a2f_2d(iene)*therm_factor*deltaene
           end if
           if (ene < e_k .and. (e_k-ene) .le. omega_max) then
             omega = e_k - ene
             if (abs(omega) < tol7) cycle
             xx = omega/(kb_HaK*Temp)
             occ_omega=1.0_dp/(exp(xx)-1.0_dp)

             therm_factor = 1 - occ_e + occ_omega

             inv_tau_k(itemp,isppol,iFSkpt,jband) = inv_tau_k(itemp,isppol,iFSkpt,jband) +&
             a2f_2d2(iene)*therm_factor*deltaene
           end if

         end do ! nene
       end do ! Temp
!          write(*,*)'i am here 2 ', isppol,iFSkpt,jband
     end do ! jband
   end do ! kpt
 end do ! nsppol

!write (300+mpi_enreg%me,*) inv_tau_k
 call xmpi_sum (inv_tau_k, xmpi_world, ierr)

 ABI_DEALLOCATE(phfrq)
 ABI_DEALLOCATE(displ)
 ABI_DEALLOCATE(pheigvec)
 ABI_DEALLOCATE(tmp2_wtk)
 ABI_DEALLOCATE(ff2)
 ABI_DEALLOCATE(ene_pt)
 ABI_DEALLOCATE(ene_ptfine)
 ABI_DEALLOCATE(tmp_wtk1)
 ABI_DEALLOCATE(tmp_wtk2)
 ABI_DEALLOCATE(tmp2_wtq)
 ABI_DEALLOCATE(wtq)
 ABI_DEALLOCATE(coskr1)
 ABI_DEALLOCATE(sinkr1)
 ABI_DEALLOCATE(coskr2)
 ABI_DEALLOCATE(sinkr2)
 ABI_DEALLOCATE(kpttokpt)
 ABI_DEALLOCATE(FSfullpktofull)
 ABI_DEALLOCATE(mqtofull)
 ABI_DEALLOCATE(tmp_gkk_qpt)
 ABI_DEALLOCATE(tmp_gkk_rpt)
 ABI_DEALLOCATE(tmp_gkk_kpt)
 ABI_DEALLOCATE(tmp_gkk_kpt2)
 ABI_DEALLOCATE(gkk_kpt)
 ABI_DEALLOCATE(a2f_2d)
 ABI_DEALLOCATE(a2f_2d2)

!output inv_tau_k and tau_k
 fname = trim(elph_ds%elph_base_name) // '_INVTAUK'
 if (open_file(fname,message,newunit=unit_invtau,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to relaxation time file
 write (unit_invtau,*) '# k-dep inverse of the relaxation time as a function of temperature.'
 write (unit_invtau,*) '# '
 write (unit_invtau,*) '# nkptirr= ', nkptirr, 'nband= ', nband
 write (unit_invtau,*) '# number of temperatures=  ', ntemper
 write (unit_invtau,*) '# tau [femtosecond^-1]     '

 fname = trim(elph_ds%elph_base_name) // '_TAUK'
 if (open_file(fname,message,newunit=unit_tau,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to relaxation time file
 write (unit_tau,*) '# k-dep relaxation time as a function of temperature.'
 write (unit_tau,*) '# '
 write (unit_tau,*) '# nkptirr= ', nkptirr, 'nband= ', nband
 write (unit_tau,*) '# number of temperatures=  ', ntemper
 write (unit_tau,*) '# tau [femtosecond]     '

 tau_k = zero
 do itemp=1,ntemper  ! runs over termperature in K
   Temp=elph_ds%tempermin+elph_ds%temperinc*dble(itemp)
   write(unit_invtau,'(a,f16.8)') '# Temperature = ', Temp
   write(unit_tau,'(a,f16.8)') '# Temperature = ', Temp
   do isppol=1,nsppol
     write(unit_invtau,'(a,i6)') '# For isppol = ', isppol
     write(unit_tau,'(a,i6)') '# For isppol = ', isppol
     do iFSkpt = 1,nkpt
!FIXME: check when tau_k is too small, whether there should be a phonon
!scattering or not, and should tau_k be zero or not.
       do jband = 1,nband
         if (abs(inv_tau_k(itemp,isppol,iFSkpt,jband)) < tol9) then
           inv_tau_k(itemp,isppol,iFSkpt,jband) = zero
           tau_k(itemp,isppol,iFSkpt,jband) = zero
         else
!no need to *nkpt due to wtkirr, as we need /nkpt for the sum
!no need to *two_pi due to the missing prefactor in gkk (see mka2f_tr_lova)
           inv_tau_k(itemp,isppol,iFSkpt,jband) = inv_tau_k(itemp,isppol,iFSkpt,jband)*elph_ds%occ_factor
           tau_k(itemp,isppol,iFSkpt,jband) = one/inv_tau_k(itemp,isppol,iFSkpt,jband)
         end if
       end do ! nband
       write(unit_invtau,'(a,i8,a,3f12.6)') '# kpt# ', iFSkpt, '   kpt=', elph_ds%k_phon%kptirr(:,iFSkpt)
       write(unit_invtau,'(100D16.8)') (inv_tau_k(itemp,isppol,iFSkpt,iband)*femto/chu_tau,iband=1,nband)
       write(unit_tau,'(a,i8,a,3f12.6)') '# kpt# ', iFSkpt, '   kpt=', elph_ds%k_phon%kptirr(:,iFSkpt)
       write(unit_tau,'(100D16.8)') (tau_k(itemp,isppol,iFSkpt,iband)*chu_tau/femto,iband=1,nband)
     end do ! nkptirr
     write(unit_invtau,*) ' '
     write(unit_tau,*) ' '
   end do ! nsppol
   write(unit_invtau,*) ' '
   write(unit_invtau,*) ' '
   write(unit_tau,*) ' '
   write(unit_tau,*) ' '
 end do ! ntemper

! Only use the irred k for eigenGS and tau_k
 ABI_ALLOCATE(tmp_eigenGS,(elph_ds%nband,elph_ds%k_phon%new_nkptirr,elph_ds%nsppol))

 do ikpt_irr = 1, new_nkptirr
   tmp_eigenGS(:,ikpt_irr,:) = eigenGS(:,elph_ds%k_phon%new_irredtoGS(ikpt_irr),:)
   tmp_tau_k(:,:,ikpt_irr,:) = tau_k(:,:,elph_ds%k_phon%new_irredtoGS(ikpt_irr),:)*chu_tau
 end do

!BoltzTraP output files in SIESTA format
 if (elph_ds%prtbltztrp == 1) then
   call ebands_prtbltztrp_tau_out (tmp_eigenGS(elph_ds%minFSband:elph_ds%maxFSband,:,:),&
&   elph_ds%tempermin,elph_ds%temperinc,ntemper,fermie, &
&   elph_ds%elph_base_name,elph_ds%k_phon%new_kptirr,natom,nband,elph_ds%nelect,new_nkptirr, &
&   elph_ds%nspinor,nsppol,Cryst%nsym,Cryst%rprimd,Cryst%symrel,tmp_tau_k)
 end if !prtbltztrp
 ABI_DEALLOCATE(tmp_eigenGS)
 ABI_DEALLOCATE(tmp_tau_k)

!Get the energy dependence of tau. 
!Eq. (6) in  Restrepo et al. Appl. Phys. Lett. 94, 212103 (2009)  

 fname = trim(elph_ds%elph_base_name) // '_TAUE'
 if (open_file(fname,message,newunit=unit_taue,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to relaxation time file
 write (unit_taue,*) '# Energy-dep relaxation time as a function of temperature.'
 write (unit_taue,*) '# '
 write (unit_taue,*) '# number of temperatures=  ', ntemper
 write (unit_taue,*) '# ene[Ha] tau [femtosecond] DOS[au]    '

 fname = trim(elph_ds%elph_base_name) // '_MFP'
 if (open_file(fname,message,newunit=unit_mfp,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

 write (unit_mfp,*) '# Energy-dep mean free path as a function of temperature.'
 write (unit_mfp,*) '# '
 write (unit_mfp,*) '# number of temperatures=  ', ntemper
 write (unit_mfp,*) '# ene[Ha] mfp [femtometer]   '

 do itemp=1,ntemper  ! runs over termperature in K
   Temp=elph_ds%tempermin+elph_ds%temperinc*dble(itemp)
   write(unit_taue,'(a,f16.8)') '# Temperature = ', Temp
   do isppol = 1, nsppol
     write(unit_taue,*) '# Tau_e for isppol = ',isppol
     do iene = 1, nene_all
       rate_e = zero
       do iFSkpt = 1, nkpt
         do jband = 1, nband
           rate_e = rate_e + inv_tau_k(itemp,isppol,iFSkpt,jband)* &
&           tmp_wtk(jband,iFSkpt,isppol,iene)
         end do ! jband
       end do ! kpt
       if (dabs(dos_e(isppol,iene)) < tol7) then
         rate_e = zero
       else
         rate_e = rate_e/nkpt/dos_e(isppol,iene)
       end if 
       write(unit_taue,"(3D16.8)") enemin+(iene-1)*deltaene*nspline, rate_e*femto/chu_tau, dos_e(isppol,iene)
     end do ! number of energies
     write(unit_taue,*) ' '
   end do ! nsppol
   write(unit_taue,*) ' '
 end do ! ntemperature

! calculate and output mean free path
 do itemp=1,ntemper  ! runs over termperature in K
   Temp=elph_ds%tempermin+elph_ds%temperinc*dble(itemp)
   write(unit_mfp,'(a,f16.8)') '# Temperature = ', Temp
   do isppol = 1, nsppol
     do icomp = 1, 3
       write(unit_mfp,*) '# Mean free path for isppol, icomp= ',isppol,icomp
       do iene = 1, nene_all
         mfp_e = zero
         do iFSkpt = 1, nkptirr
           do jband = 1, nband
             mfp_e = mfp_e + tau_k(itemp,isppol,iFSkpt,jband)* &
&             elph_tr_ds%el_veloc(iFSkpt,elph_ds%minFSband+jband-1,icomp,isppol)* &
&             tmp_wtk(jband,iFSkpt,isppol,iene)
!&                          elph_ds%k_phon%new_wtkirr(iFSqpt)
           end do ! jband
         end do ! kpt
         if (dabs(dos_e(isppol,iene)) < tol7) then
           mfp_e = zero
         else
           mfp_e = mfp_e/nkptirr/dos_e(isppol,iene)
         end if 
         write(unit_mfp,"(2D16.8)") enemin+(iene-1)*deltaene*nspline, mfp_e*chu_mfp/femto
       end do ! number of energies
       write(unit_mfp,*) ' '
     end do ! icomp
     write(unit_mfp,*) ' '
   end do ! nsppol
   write(unit_mfp,*) ' '
 end do ! ntemperature

 ABI_ALLOCATE(cond_e ,(ntemper,nsppol,nene_all,9))

!get cond_e
 cond_e = zero
 do itemp=1,ntemper  ! runs over termperature in K
   do isppol = 1, nsppol
     do iene = 1, nene_all
!       do iFSkpt =1,nkpt
       do ik_this_proc =1,elph_ds%k_phon%my_nkpt
         iFSkpt = elph_ds%k_phon%my_ikpt(ik_this_proc)
         do jband = 1, nband
           do icomp = 1, 3
             do jcomp = 1, 3
               itensor = (icomp-1)*3+jcomp
               cond_e(itemp,isppol,iene,itensor) = cond_e(itemp,isppol,iene,itensor) + &
&               tau_k(itemp,isppol,iFSkpt,jband)* &
&               elph_tr_ds%el_veloc(iFSkpt,elph_ds%minFSband+jband-1,icomp,isppol)* &
&               elph_tr_ds%el_veloc(iFSkpt,elph_ds%minFSband+jband-1,jcomp,isppol)* &
&               tmp_wtk(jband,iFSkpt,isppol,iene)
             end do
           end do
         end do ! jband
       end do ! kpt
     end do ! number of energies
   end do ! nsppol
 end do ! ntemperature

 ! MG FIXME: Why xmpi_world, besides only master should perform IO in the section below.
 call xmpi_sum (cond_e, xmpi_world, ierr)

 cond_e = cond_e/nkpt

!get transport coefficients

 fname = trim(elph_ds%elph_base_name) // '_COND'
 if (open_file(fname,message,newunit=unit_cond,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to conductivity file
 write (unit_cond,*) '#  Conductivity as a function of temperature.'
 write (unit_cond,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_cond,*) '#  '
 write (unit_cond,*) '#  Columns are: '
 write (unit_cond,*) '#  temperature[K]   cond[au]   cond [SI]    '
 write (unit_cond,*) '#  '

 fname = trim(elph_ds%elph_base_name) // '_CTH'
 if (open_file(fname,message,newunit=unit_therm,status='unknown') /= 0) then
   MSG_ERROR(message)
 end if

!print header to thermal conductivity file
 write (unit_therm,'(a)') '# Thermal conductivity as a function of temperature.'
 write (unit_therm,'(a)') '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_therm,'(a)') '#  '
 write (unit_therm,'(a)') '#  Columns are: '
 write (unit_therm,'(a)') '#  temperature[K]   thermal cond [au]   thermal cond [SI]'
 write (unit_therm,'(a)') '#  '

 fname = trim(elph_ds%elph_base_name) // '_SBK'
 if (open_file(fname,message,newunit=unit_sbk,status='unknown') /=0) then
   MSG_ERROR(message)
 end if

!print header to relaxation time file
 write (unit_sbk,*) '# Seebeck Coefficint as a function of temperature.'
 write (unit_sbk,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_sbk,*) '#  '
 write (unit_sbk,*) '#  Columns are: '
 write (unit_sbk,*) '#  temperature[K]   S [au]   S [SI]     '
 write (unit_sbk,*) '#  '

 ABI_ALLOCATE(cond ,(ntemper,nsppol,3,3))
 ABI_ALLOCATE(cth ,(ntemper,nsppol,3,3))
 ABI_ALLOCATE(sbk ,(ntemper,nsppol,3,3))
 ABI_ALLOCATE(seebeck ,(ntemper,nsppol,3,3))

 cond = zero
 cth = zero
 sbk = zero
 seebeck = zero
 do isppol=1,nsppol
   do icomp=1, 3
     do jcomp=1, 3
       itensor=(icomp-1)*3+jcomp
       do itemp=1,ntemper
         Temp=elph_ds%tempermin + elph_ds%temperinc*dble(itemp)
         do iene = 1, nene_all
           factor = (enemin+(iene-1)*deltaene*nspline - fermie(itemp))/(kb_HaK*Temp)
           if (factor < -40.0d0) then
             dfermide = zero
           else if (factor > 40.0d0) then
             dfermide = zero
           else
             dfermide = EXP(factor)/(kb_HaK*Temp*(EXP(factor)+one)**2.0d0)
           end if
           cond(itemp,isppol,icomp,jcomp) = cond(itemp,isppol,icomp,jcomp) + &
&           cond_e(itemp,isppol,iene,itensor)*dfermide*deltaene*nspline
           cth(itemp,isppol,icomp,jcomp) = cth(itemp,isppol,icomp,jcomp) + cond_e(itemp,isppol,iene,itensor)* &
&           (enemin+(iene-1)*deltaene*nspline - fermie(itemp))**2.0d0*dfermide*deltaene*nspline
           sbk(itemp,isppol,icomp,jcomp) = sbk(itemp,isppol,icomp,jcomp) + cond_e(itemp,isppol,iene,itensor)* &
&           (enemin+(iene-1)*deltaene*nspline - fermie(itemp))*dfermide*deltaene*nspline
         end do
       end do ! temperature
     end do ! jcomp
   end do ! icomp
 end do !end isppol

 do isppol=1,nsppol
   do itemp=1,ntemper
     cond_inv(:,:)=cond(itemp,isppol,:,:)
     call matrginv(cond_inv,3,3)
     call DGEMM('N','N',3,3,3,one,sbk(itemp,isppol,:,:),3,cond_inv,&
&     3,zero,seebeck(itemp,isppol,:,:),3)
   end do 
 end do

 do isppol=1,nsppol
   do icomp=1, 3
     do jcomp=1, 3
       itensor=(icomp-1)*3+jcomp
       write(unit_cond,*) '# Conductivity for isppol, itrten= ',isppol,itensor
       write(unit_therm,*) '# Thermal conductivity for isppol, itrten= ',isppol,itensor
       write(unit_sbk,*) '# Seebeck coefficient for isppol, itrten= ',isppol,itensor
       do itemp=1,ntemper
         Temp=elph_ds%tempermin + elph_ds%temperinc*dble(itemp)

         seebeck(itemp,isppol,icomp,jcomp) = -1.0d0*seebeck(itemp,isppol,icomp,jcomp)/(kb_HaK*Temp)
         cond(itemp,isppol,icomp,jcomp) = cond(itemp,isppol,icomp,jcomp)/cryst%ucvol
         cth(itemp,isppol,icomp,jcomp) = cth(itemp,isppol,icomp,jcomp)/(kb_HaK*Temp)/cryst%ucvol
         write(unit_cond,'(3D20.10)')Temp,cond(itemp,isppol,icomp,jcomp),cond(itemp,isppol,icomp,jcomp)*chu_cond
         write(unit_therm,'(3D20.10)')Temp,cth(itemp,isppol,icomp,jcomp),cth(itemp,isppol,icomp,jcomp)*chu_cth
         write(unit_sbk,'(3D20.10)')Temp,seebeck(itemp,isppol,icomp,jcomp),seebeck(itemp,isppol,icomp,jcomp)*chu_sbk
       end do ! temperature
       write(unit_cond,*)
       write(unit_therm,*)
       write(unit_sbk,*)
     end do ! jcomp
   end do ! icomp
 end do !end isppol


 ABI_DEALLOCATE(inv_tau_k)
 ABI_DEALLOCATE(tau_k)
 ABI_DEALLOCATE(tmp_wtk)
 ABI_DEALLOCATE(dos_e)
 ABI_DEALLOCATE(cond_e)
 ABI_DEALLOCATE(fermie)
 ABI_DEALLOCATE(cond)
 ABI_DEALLOCATE(sbk)
 ABI_DEALLOCATE(cth)
 ABI_DEALLOCATE(seebeck)

 close (unit=unit_tau)
 close (unit=unit_taue)
 close (unit=unit_mfp)
 close (unit=unit_invtau)
 close (unit=unit_cond)
 close (unit=unit_therm)
 close (unit=unit_sbk)

end subroutine get_tau_k
!!***
