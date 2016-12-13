!{\src2tex{textfont=tt}}
!!****f* ABINIT/read_gkk
!!
!! NAME
!! read_gkk
!!
!! FUNCTION
!! This routine reads in elphon matrix elements and completes them
!! using the appropriate symmetries
!!
!! COPYRIGHT
!! Copyright (C) 2004-2016 ABINIT group (MVer, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  elph_ds = datastructure containing elphon matrix elements
!!  Cryst<crystal_t>=Info on the crystal unit cell.
!!  Ifc<ifc_type>=Object containing the interatomic force constants.
!!  FSfullpqtofull = mapping of k+q to k
!!  n1wf = number of 1WF files to be read and analyzed
!!  nband = number of bands per kpoint
!!  unitgkk = unit of GKK file for reading
!!
!! OUTPUT
!!  elph_ds = modified gkq
!!  gkk_qpt = el-ph matrix elements for irreducible qpoints and
!!    kpoints (as a function of the reduced symmetry for the qpoint)
!!  gkk_flag = flag array:
!!       -1 -> element is missing
!!        0 -> element is from symmetric qpt (Now done in complete_gkk)
!!        1 -> element is from symmetric pert
!!        2 -> element is kptsym of gkk file
!!        3 -> element was read from gkk file
!!
!! PARENTS
!!      get_all_gkq
!!
!! CHILDREN
!!      completeperts,get_rank_1kpt,hdr_bcast,hdr_fort_read,hdr_free,ifc_fourq
!!      littlegroup_pert,littlegroup_q,mati3inv,normsq_gkq,phdispl_cart2red
!!      prt_gkk_yambo,wrap2_pmhalf,wrtout,xmpi_bcast
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine read_gkk(elph_ds,Cryst,ifc,Bst,FSfullpqtofull,gkk_flag,n1wf,nband,ep_prt_yambo,unitgkk)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_elphon
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_kptrank
 use m_hdr

 use m_numeric_tools,   only : wrap2_pmhalf
 use m_geometry,        only : phdispl_cart2red
 use m_crystal,         only : crystal_t
 use m_ifc,             only : ifc_type, ifc_fourq

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'read_gkk'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_77_ddb, except_this_one => read_gkk
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1wf,nband,unitgkk,ep_prt_yambo
 type(crystal_t),intent(in) :: Cryst
 type(ifc_type),intent(in) :: ifc
 type(ebands_t),intent(in) :: Bst
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 integer,intent(out) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol,elph_ds%nqpt_full)

!Local variables-------------------------------
!scalars
 integer :: nsppol,nbranch,nFSband,minFSband,comm,use_sym
 integer :: fform,i1wf,ikpt_phon,iatom1,iatom2
 integer :: ib,ib1,ib2,ibb,ibranch,idir,idir1,idir2,ierr,ii,ikpt1
 integer :: ipert,ipert1,ipert2,iqptirred,iqptfull,isppol,isym1
 integer :: itim1,jkpt_phon,new
 integer :: nsym1,qtimrev,syuse
 integer :: tdonecompl,test_flag,verify
 integer :: nqptirred_local
 integer :: master, me
 integer :: symrankkpt, ikpt1_phon, ik_this_proc
 real(dp) :: res,ss,timsign
 character(len=500) :: msg
 type(hdr_type) :: hdr1
!arrays
 integer :: FSirrtok(3,elph_ds%k_phon%nkpt)
 integer :: symaf1(Cryst%nsym),symq(4,2,Cryst%nsym)
 integer :: symrc1(3,3,Cryst%nsym),symrl1(3,3,Cryst%nsym)
 integer :: tmpflg(3,Cryst%natom+2,3,Cryst%natom+2)
 real(dp) :: displ_cart(2,3*Cryst%natom,3*Cryst%natom)
 real(dp) :: displ_red(2,3*Cryst%natom,3*Cryst%natom)
 real(dp) :: eigvec(2,3*Cryst%natom,3*Cryst%natom),kpt(3),phfrq_tmp(3*Cryst%natom),redkpt(3)
 real(dp) :: qptirred_local(3,n1wf)
 real(dp) :: tnons1(3,Cryst%nsym)
 real(dp),allocatable :: eigen1(:,:,:),gkk_qpt_tmp(:,:,:,:)
 real(dp),allocatable :: h1_mat_el(:,:,:,:,:),h1_mat_el_sq(:,:,:,:,:)
 real(dp),allocatable :: qdata(:,:,:),qdata_tmp(:,:,:,:)

! *************************************************************************

 ABI_UNUSED(Bst%bantot)

 use_sym   = 1
 nsppol    = elph_ds%nsppol
 nbranch   = elph_ds%nbranch
 if (ep_prt_yambo==1) then
   nFSband = nband
   minFSband = 1
 else
   nFSband   = elph_ds%nFSband
   minFSband = elph_ds%minFSband
 end if

!init values for parallelization
 comm = xmpi_world
 me = xmpi_comm_rank(comm)
 master = 0

 ABI_STAT_ALLOCATE(h1_mat_el,(2, nFSband**2, nbranch, elph_ds%k_phon%my_nkpt, nsppol), ierr)
 ABI_CHECK(ierr==0, 'trying to allocate array h1_mat_el')
 h1_mat_el= zero

 ABI_STAT_ALLOCATE(h1_mat_el_sq,(2, nFSband**2, nbranch**2,elph_ds%k_phon%my_nkpt, nsppol), ierr)
 ABI_CHECK(ierr==0, 'trying to allocate array h1_mat_el_sq')
 h1_mat_el_sq = zero

 ABI_ALLOCATE(elph_ds%qirredtofull,(elph_ds%nqptirred))

!MG array to store the e-ph quantities calculated over the input Q-grid
 ABI_ALLOCATE(qdata_tmp,(elph_ds%nqptirred,nbranch,nsppol,3))
 qdata_tmp=zero

 nqptirred_local=0 !zero number of irred q-points found
 qptirred_local(:,:)=zero

 gkk_flag = -1

 if (elph_ds%gkqwrite ==0) then
   elph_ds%gkk_qpt = zero

 else if (elph_ds%gkqwrite == 1) then
   ABI_STAT_ALLOCATE(gkk_qpt_tmp,(2,elph_ds%ngkkband**2,nbranch**2,nsppol), ierr)
   ABI_CHECK(ierr==0, 'trying to allocate array gkk_qpt_tmp')
   gkk_qpt_tmp = zero
   do iqptirred=1,elph_ds%nqptirred*elph_ds%k_phon%nkpt
     write (elph_ds%unitgkq,REC=iqptirred) gkk_qpt_tmp
   end do
   ABI_DEALLOCATE(gkk_qpt_tmp)

 else
   write (msg,'(a,i0)')' Wrong values for gkqwrite = ',elph_ds%gkqwrite
   MSG_BUG(msg)
 end if !gkqwrite

!===========================================================
!Loop over all files we have
!read in header for perturbation
!should check that all files are complete, have same header
!(taking into account the symmetries for the qpoint),
!represent the correct qpoints ...
!MG: this task should be performed in mrggkk
!===========================================================

 ABI_ALLOCATE(eigen1,(2,nband,nband))
 do i1wf=1,n1wf

   if (master == me) then
     write (msg,'(2a,i4,a,i4)')ch10,' read_gkk : reading 1WF header # ',i1wf,' /',n1wf
     call wrtout(std_out,msg,'COLL')

!    Could check for compatibility of natom, kpt grids, ecut, qpt with DDB grid...
!    MG: Also this task should be done in mrggkk

     call hdr_fort_read(hdr1, unitgkk, fform)
     if (fform == 0) then
       write (msg,'(a,i0,a)')' 1WF header number ',i1wf,' was mis-read. fform == 0'
       MSG_ERROR(msg)
     end if

     write(msg,'(a,i4)')' read_gkk : have read 1WF header #',i1wf
     call wrtout(std_out,msg,'COLL')
     write (msg,'(2a,i4,a)')ch10,' read_gkk : # of kpt for this perturbation: ',hdr1%nkpt,ch10
     call wrtout(std_out,msg,'COLL')

   end if

!  broadcast data to all nodes:
   call hdr_bcast(hdr1, master, me, comm)

!  Find qpoint in full grid
   new=1
   do iqptfull=1,elph_ds%nqpt_full
     kpt(:) = hdr1%qptn(:) - elph_ds%qpt_full(:,iqptfull)
     call wrap2_pmhalf(kpt(1),redkpt(1),res)
     call wrap2_pmhalf(kpt(2),redkpt(2),res)
     call wrap2_pmhalf(kpt(3),redkpt(3),res)
     ss=redkpt(1)**2+redkpt(2)**2+redkpt(3)**2
     if(ss < tol6) then
       new = 0
       exit !exit with iqptfull
     end if
   end do !iqptfull

   if (new == 1) then
!    Test should be at the end: dont care if there are additional
!    qpts in gkk file which are not on the main grid. Ignore them.
     write (msg,'(4a,3es16.6,2a)')ch10,&
&     ' read_gkk : WARNING-  ',ch10,&
&     ' qpoint = ',hdr1%qptn(:),ch10,&
&     ' not found in the input q-grid. Ignoring this point '
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,msg,'COLL')
     if (me == master) then
       do isppol=1,hdr1%nsppol
         do ikpt1=1,hdr1%nkpt
           read(unitgkk) ((eigen1(:,ii,ib),ii=1,nband),ib=1,nband)
         end do
       end do
     end if

     cycle !cycle the loop on i1wf
   end if !end if (new ==1)


!  Check whether other pieces of the DDB have used this qpt already
   new=1
   do iqptirred=1,nqptirred_local
     kpt(:) = qptirred_local(:,iqptirred) - hdr1%qptn(:)
     call wrap2_pmhalf(kpt(1),redkpt(1),res)
     call wrap2_pmhalf(kpt(2),redkpt(2),res)
     call wrap2_pmhalf(kpt(3),redkpt(3),res)
     ss=redkpt(1)**2+redkpt(2)**2+redkpt(3)**2
     if(ss < tol6) then
       new=0
       exit  !MG We can use this information to avoid recalculating the dynamical matrix
     end if !but we need to use a fixed format in GKK!
   end do !iqptirred

   if (new==1) then  !we have a new valid irreducible qpoint, add it!
     nqptirred_local = nqptirred_local+1
     if (nqptirred_local > elph_ds%nqptirred) then
       write (msg, '(a,a,a,i6,i6)') &
&       'found too many qpoints in GKK file wrt anaddb input ', ch10, &
&       'nqpt_anaddb nqpt_gkk = ', elph_ds%nqptirred, nqptirred_local
       MSG_ERROR(msg)
     end if
     qptirred_local(:,nqptirred_local) = hdr1%qptn(:)
     iqptirred = nqptirred_local
     tdonecompl = 0
     h1_mat_el = zero
   end if

!  now iqptirred is the index of the present qpoint in the array qptirred_local
!  and iqptfull is the index in the full qpt_full array for future reference
   elph_ds%qirredtofull(iqptirred) = iqptfull

   write (msg,'(a,i5,a,3es16.8)')&
&   ' read_gkk : full zone qpt number ',iqptfull,' is ',elph_ds%qpt_full(:,iqptfull)
   call wrtout(std_out,msg,'COLL')

!  if this perturbation has already been filled (overcomplete gkk)
!  check only 1st kpoint and spinpol, then check others
   verify = 0
   if (gkk_flag(hdr1%pertcase,hdr1%pertcase,1,1,elph_ds%qirredtofull(iqptirred)) /= -1) then
!
     do isppol=1,nsppol
       do ik_this_proc=1,elph_ds%k_phon%my_nkpt
         if (gkk_flag(hdr1%pertcase,hdr1%pertcase,ik_this_proc,isppol,elph_ds%qirredtofull(iqptirred)) == -1) then
           write (std_out,*)" hdr1%pertcase,ik_this_proc,iqptirred",hdr1%pertcase,ik_this_proc,iqptirred
           MSG_ERROR('Partially filled perturbation ')
         end if
       end do ! ikpt_phon
     end do ! isppol
!
     MSG_WARNING(' gkk perturbation is already filled')
     write(std_out,*)' hdr1%pertcase,iqptirred,iqptfull = ',hdr1%pertcase,iqptirred,iqptfull,&
&     gkk_flag(hdr1%pertcase,hdr1%pertcase,1,1,elph_ds%qirredtofull(iqptirred))
     verify = 1
     write (125,*) '# matrix elements for symmetric perturbation'
!    Instead of reading eigen1 into void, verify == 1 checks them later on wrt values in memory
   end if !gkk_flag

!  Examine the symmetries of the q wavevector
!  these will be used to complete the perturbations for other atoms and idir
   if (ep_prt_yambo==1) then
     ! If one wants to print GKKs along phonon modes, it mean mixing of
     ! perturbations with differnt jauge. Symmetries must then be disable.
     call littlegroup_q(Cryst%nsym,qptirred_local(:,iqptirred),symq,Cryst%symrec,Cryst%symafm,qtimrev,prtvol=0,use_sym=0)
   else
     call littlegroup_q(Cryst%nsym,qptirred_local(:,iqptirred),symq,Cryst%symrec,Cryst%symafm,qtimrev,prtvol=0)
   end if

   ! Determine dynamical matrix, phonon frequencies and displacement vector for qpoint
   !call wrtout(std_out,' read_gkk: calling inpphon to calculate the dynamical matrix','COLL')

   call ifc_fourq(ifc,cryst,qptirred_local(:,iqptirred),phfrq_tmp,displ_cart,out_eigvec=eigvec)

!  Get displacement vectors for all branches in reduced coordinates
!  used in scalar product with H(1)_atom,idir  matrix elements
!  Calculate $displ_red = displ_cart \cdot gprimd$ for each phonon branch

   call phdispl_cart2red(Cryst%natom,Cryst%gprimd,displ_cart,displ_red)

!  prefactors for gk+q,n\prime;k,n matrix element
!  COMMENT : in decaft there is a weird term in the mass factor, of M-zval(species)
!  dont know why. Not needed to reproduce decaft results, though...
!  weight is squared in evaluation of
!  gamma_{k,q,j} = 2 \pi omega_{q,j} sum_{nu,nu\prime} |g^{q,j}_{k+q,nu\prime; k,nu}|^2
!  normally cancels with the 2 \pi omega_{q,j} factor in front of the sum...

!  hdr1%pertcase = idir + (ipert-1)*3 where ipert=iatom in the interesting cases
   idir = mod (hdr1%pertcase-1,3)+1
   ipert = int(dble(hdr1%pertcase-idir)/three)+1

   write (msg,'(4a,i3,a,i3,a,i4,a)')ch10,&
&   ' read_gkk : calling littlegroup_pert to examine the symmetries of the full perturbation ',ch10,&
&   ' idir = ',idir,' ipert = ',ipert,' and Q point = ',iqptirred,ch10
   call wrtout(std_out,msg,'COLL')

!  Examine the symmetries of the full perturbation these will be used to complete the kpoints
!  DOESNT USE TIME REVERSAL IN littlegroup_pert except for gamma

   syuse=0

   call littlegroup_pert(Cryst%gprimd,idir,Cryst%indsym,ab_out,ipert,Cryst%natom,Cryst%nsym,nsym1,2,Cryst%symafm,symaf1,&
&   symq,Cryst%symrec,Cryst%symrel,symrl1,syuse,Cryst%tnons,tnons1)

   do isym1=1,nsym1
     call mati3inv(symrl1(:,:,isym1),symrc1(:,:,isym1))
   end do
   FSirrtok = 0

!  ========================================================
!  Loop over irred kpts in file, and fill the default gkk
!  ========================================================

!  MG NOTE : in the present implementation, if nsppol /=1 the code stops in rchkGSheader!
   do isppol=1,hdr1%nsppol !Loop over spins is trivial? Not tested.
     write (std_out,*) ' read_gkk : isppol = ', isppol

     do ikpt1=1,hdr1%nkpt   !Loop over irred kpoints, WARNING  nkpt depends on qpoint and symmetry!
!
!      this is the main read of the gkk matrix elements from the file (eigen1 arrays)
!      it has to be done exactly nsppol*nkpt times, and the kpt_phon are completed
!      where appropriate in the loop below (normally succeeding only once for each kpt)
!
       if (master == me) then
         read(unitgkk) ((eigen1(:,ii,ib),ii=1,nband),ib=1,nband)
       end if

!      MPI broadcast data to all nodes:
       call xmpi_bcast(eigen1, master, comm, ierr)

!      find place of irred k in k_phon
!      the kpoints in the file (kptns) could be ordered arbitrarily
       call get_rank_1kpt (hdr1%kptns(:,ikpt1)-qptirred_local(:,iqptirred), &
&       symrankkpt, elph_ds%k_phon%kptrank_t)
       ikpt1_phon = elph_ds%k_phon%kptrank_t%invrank(symrankkpt)
       if (ikpt1_phon < 0) then
         write (msg,'(a,3es16.6,a)')&
&         ' irred k ',hdr1%kptns(:,ikpt1),' was not found in full grid'
         MSG_ERROR(msg)
       end if
!      find correspondence between this kpt_phon and the others
!      symrc1 conserves perturbation as well as qpoint
!      add to FSirrtok list
       do isym1=1,nsym1
         do itim1=0,qtimrev
           timsign=one-two*itim1
           kpt(:) = timsign*matmul(symrc1(:,:,isym1), elph_ds%k_phon%kpt(:,ikpt1_phon))

           call get_rank_1kpt (kpt,symrankkpt,elph_ds%k_phon%kptrank_t)
           jkpt_phon = elph_ds%k_phon%kptrank_t%invrank(symrankkpt)

           if (jkpt_phon > 0) then
             FSirrtok(1,jkpt_phon) = ikpt1_phon
             FSirrtok(2,jkpt_phon) = isym1
             FSirrtok(3,jkpt_phon) = itim1
           else
             write (msg,'(a,3es16.6,a,i5,a,i4,a)')&
&             ' sym equivalent of kpt ',hdr1%kptns(:,ikpt1),' by sym ',&
&             isym1,' and itime ',itim1,' was not found'
             MSG_ERROR(msg)
           end if
         end do !itim1
       end do !isim1


       !
       !  Here check if the symmetry-copied gkk at new k point is equal to the one found in the file for non-irreducible point
       !  NB This is DEBUG code
       !
       if (verify == 1 .and. elph_ds%k_phon%my_kpt(ikpt1_phon) == me) then
         do ik_this_proc = 1, elph_ds%k_phon%my_nkpt
           if (elph_ds%k_phon%my_ikpt(ik_this_proc) == ikpt1_phon) exit
         end do
         do ib1=1,nFSband
           do ib2=1,nFSband
             ibb = (ib1-1)*nFSband+ib2
             write (125,'(2(2E16.6,2x))') h1_mat_el(:,ibb,hdr1%pertcase,ik_this_proc,isppol),&
&             eigen1(:,minFSband-1+ib2,minFSband-1+ib1)
           end do
         end do
       end if !verify end DEBUG code


       do ik_this_proc = 1, elph_ds%k_phon%my_nkpt
!        should I be dealing with this k-point?
         jkpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

!        does present ikpt1 contribute to this k-point?
         if (FSirrtok(1,jkpt_phon) /= ikpt1_phon) cycle

!        if this kpoint has already been filled (overcomplete gkk)
         if (gkk_flag(hdr1%pertcase,hdr1%pertcase,ik_this_proc,isppol,elph_ds%qirredtofull(iqptirred)) /= -1) then
           MSG_WARNING("gkk element is already filled")
           write(std_out,*)' hdr1%pertcase,ik_this_proc,isppol,iqptirred = ',&
&           hdr1%pertcase,ik_this_proc,isppol,iqptirred,&
&           gkk_flag(hdr1%pertcase,hdr1%pertcase,ik_this_proc,isppol,elph_ds%qirredtofull(iqptirred))
!           exit
         end if !gkk_flag

!        ===============================================================
!        TODO: if there is a phase factor in swapping k-points, insert it here in copy to h1_mat_el
!        as a function of symops in FSirrtok
!        complete gkk for symmetric ikpt_phon with sym1 which conserve
!        the full perturbation+qpoint
!        Not tested explicitly, but the results for Pb using reduced kpts look good
!        should do same RF calculation with nsym=1 and check
!        ===============================================================

!        save this kpoint
         do ib1=1,nFSband
           do ib2=1,nFSband
             ibb = (ib1-1)*nFSband+ib2

!            real
             res=eigen1(1,minFSband-1+ib2,minFSband-1+ib1)
             h1_mat_el(1,ibb,hdr1%pertcase,ik_this_proc,isppol) = res

!            imag
             res=eigen1(2,minFSband-1+ib2,minFSband-1+ib1)
             h1_mat_el(2,ibb,hdr1%pertcase,ik_this_proc,isppol) = res
           end do !ib2
         end do !ib1
!        if jkpt is equal to ikpt1_phon (if clause above) flag == 3
         if (FSirrtok(2,jkpt_phon) == 1) then
           gkk_flag(hdr1%pertcase,hdr1%pertcase,ik_this_proc,isppol,elph_ds%qirredtofull(iqptirred)) = 3
!          if jkpt_phon comes from ikpt1_phon flag == 2 with some symop
         else
           gkk_flag(hdr1%pertcase,hdr1%pertcase,ik_this_proc,isppol,elph_ds%qirredtofull(iqptirred)) = 2
         end if

       end do !jkpt_phon

!      ===============================================================
!      we now have contribution to g(k+q,k; \kappa,\alpha) from one
!      kpoint,and one perturbation,
!      NB: each perturbation will contribute to all the modes later!
!
!      SHOULD ONLY DO THIS FOR THE SYMS NEEDED
!      TO COMPLETE THE PERTURBATIONS!!!
!      ================================================================

     end do !ikpt1
   end do !isppol

! 14 Jan 2014 removed test on verify - in new scheme full BZ is read in and should be used to avoid phase errors
!   if (verify == 1) cycle

!  Checks on irred grid provided and on gkk_flag accumulated up to now
   if (elph_ds%tuniformgrid == 1) then  ! check if irred kpoints found reconstitute the FS kpts
     do ikpt_phon=1,elph_ds%k_phon%nkpt
       if (FSirrtok(1,ikpt_phon) == 0) then
         write(msg,'(a,3es16.6,2a)')&
&         ' kpt = ',elph_ds%k_phon%kpt(:,ikpt_phon),ch10,&
&         ' is not the symmetric of one of those found in the GKK file'
         MSG_ERROR(msg)
       end if
     end do !ikpt_phon

!    normally at this point we have used all the gkk for all kpoints on the FS
!    for the given irred perturbation: check
     do ik_this_proc = 1, elph_ds%k_phon%my_nkpt
       ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

       if (gkk_flag(hdr1%pertcase, hdr1%pertcase, ik_this_proc, 1, elph_ds%qirredtofull(iqptirred)) == -1) then
         write (msg,'(a,i3,a,3es18.6,2a,i3,a,i3,a,3es18.6,a,a,i4,a,a)')&
&         ' For irreducible qpt ', iqptirred,' = ',qptirred_local(:,iqptirred),ch10, &
&         ' the gkk element : pertcase = ',hdr1%pertcase,' ik_this_proc = ',ik_this_proc, &
&         ' kpt = ',elph_ds%k_phon%kpt(:,ikpt_phon),ch10,&
&         ' and isppol ',1,ch10,&
&         ' was not found by symmetry operations on the irreducible kpoints given'
         MSG_ERROR(msg)
       end if
     end do !ikpt_phon
   end if ! end elph_ds%tuniformgrid == 1 checks

   write(msg,'(a,i0)')' read_gkk : Done completing the kpoints for pertcase ',hdr1%pertcase
   call wrtout(std_out,msg,'COLL')

   tmpflg(:,:,:,:) = 0

   do idir1=1,3
     do iatom1=1,Cryst%natom
       ipert1 = (iatom1-1)*3+idir1
       do idir2=1,3
         do iatom2=1,Cryst%natom
           ipert2 = (iatom2-1)*3+idir2
           if (gkk_flag(ipert1,ipert1,1,1,elph_ds%qirredtofull(iqptirred)) >= 0 .and. &
&           gkk_flag(ipert2,ipert2,1,1,elph_ds%qirredtofull(iqptirred)) >= 0) then
             tmpflg(idir1,iatom1,idir2,iatom2) = 1
           end if
         end do
       end do
     end do
   end do


!  ===============================================
!  Full test: need all perturbations explicitly
!  ===============================================

   test_flag = 0
   if (sum(tmpflg(:,1:Cryst%natom,:,1:Cryst%natom)) == (3*Cryst%natom)**2 .and. tdonecompl == 0) test_flag = 1

   write(std_out,*)'read_gkk: tdonecompl = ', tdonecompl

!  de-activate completion of perts by symmetry for now.
!  Must be called when all irreducible perturbations are in memory!!!!
   if (test_flag == 1 .and. tdonecompl == 0) then

!    write(std_out,*) ' read_gkk : enter fxgkkphase before completeperts'
!    call fxgkkphase(elph_ds,gkk_flag,h1_mat_el,iqptirred)

     if (ep_prt_yambo==1) then
       if (elph_ds%k_phon%my_nkpt /= elph_ds%k_phon%nkpt) then
         write (msg, '(a)') 'prt_gkk_yambo can not handle parallel anaddb yet'
         MSG_ERROR(msg)
       end if
       call prt_gkk_yambo(displ_cart,displ_red,elph_ds%k_phon%kpt,h1_mat_el,iqptirred,&
&       Cryst%natom,nFSband,elph_ds%k_phon%my_nkpt,phfrq_tmp,hdr1%qptn)
     end if

!    ========================================================================
!    Now use more general symops to complete the other equivalent
!    perturbations: the kpoints are also shuffled by these symops
!    afterwards h1_mat_el_sq contains gamma_\tau\alpha,\tau'\alpha' in reduced coordinates
!
!    \gamma_{\tau'\alpha',\tau\alpha} =
!    <psi_{k+q,ib2}| H(1)_{\tau'\alpha'}| psi_{k,ib1}>* \cdot
!    <psi_{k+q,ib2}| H(1)_{\tau \alpha }| psi_{k,ib1}>
!
!    ========================================================================

     call completeperts(Cryst,nbranch,nFSband,elph_ds%k_phon%my_nkpt,nsppol,&
&     gkk_flag(:,:,:,:,elph_ds%qirredtofull(iqptirred)),h1_mat_el,h1_mat_el_sq,qptirred_local(:,iqptirred),symq,qtimrev)

     tdonecompl = 1
   end if

!  ==============================================================
!  if we have all the perturbations for this qpoint, proceed
!  with scalar product, norm squared, and add weight factors
!
!  SHOULD HAVE A TEST SO h1_mat_el IS NOT OVERWRITTEN
!  BEFORE PREVIOUS QPOINT IS FINISHED!!!!!
!  ==============================================================

   test_flag = 1
   do isppol=1,nsppol
     do ik_this_proc = 1, elph_ds%k_phon%my_nkpt
       do ibranch=1,nbranch
         if (gkk_flag (ibranch,ibranch,ik_this_proc,isppol,elph_ds%qirredtofull(iqptirred)) == -1) then
           test_flag = 0
           exit
         end if
       end do
     end do
   end do

   if (test_flag /= 0) then
     call wrtout(std_out,' read_gkk : enter normsq_gkq',"COLL")

!    MG temporary array to save ph-linewidths before Fourier interpolation
     ABI_ALLOCATE(qdata,(nbranch,nsppol,3))
     qdata(:,:,:)=zero

     call normsq_gkq(displ_red,eigvec,elph_ds,FSfullpqtofull,&
&     h1_mat_el_sq,iqptirred,phfrq_tmp,qptirred_local,qdata)

!    save gkk_qpt, eventually to disk, for bands up to ngkkband,
!    NB: if the sum over bands has been performed ngkkband is 1 instead of nFSband
     if (elph_ds%gkqwrite == 0) then
       elph_ds%gkk_qpt(:,:,:,:,:,iqptirred) = h1_mat_el_sq(:,1:elph_ds%ngkkband*elph_ds%ngkkband,:,:,:)
     else
!      write all kpoints to disk
       write (std_out,*) 'size of record to be written: ', 8  * 2*elph_ds%ngkkband*elph_ds%ngkkband*&
&       elph_ds%nbranch*elph_ds%nbranch*elph_ds%k_phon%my_nkpt*elph_ds%nsppol
       inquire(unit=elph_ds%unitgkq, recl=isppol)
       write (std_out,*) 'recl =', isppol
       write (std_out,*) 'iqptirred ', iqptirred
       do ik_this_proc = 1, elph_ds%k_phon%my_nkpt
         write (elph_ds%unitgkq,REC=((iqptirred-1)*elph_ds%k_phon%my_nkpt+ik_this_proc)) &
&         h1_mat_el_sq(:,1:elph_ds%ngkkband*elph_ds%ngkkband,:,ik_this_proc,:)
       end do
     end if

     qdata_tmp(iqptirred,:,:,:)=qdata(:,:,:)
     ABI_DEALLOCATE(qdata)
   end if

   call hdr_free(hdr1)

 end do !of i1wf

!got all the gkk perturbations

 ABI_DEALLOCATE(eigen1)
 ABI_DEALLOCATE(h1_mat_el)
 ABI_DEALLOCATE(h1_mat_el_sq)

 if (nqptirred_local /= elph_ds%nqptirred) then
   write (msg, '(3a,i0,i0)') &
&   ' Found wrong number of qpoints in GKK file wrt anaddb input ', ch10, &
&   ' nqpt_anaddb nqpt_gkk = ', elph_ds%nqptirred, nqptirred_local
   MSG_ERROR(msg)
 end if

!normally at this point we have the gkk for all kpoints on the FS
!for all the perturbations. Otherwise a 1WF file is missing.
!NOTE: still havent checked the qpoint grid completeness
 do iqptirred=1,elph_ds%nqptirred
   do isppol=1,nsppol
     do ik_this_proc = 1, elph_ds%k_phon%my_nkpt
       ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)
       do ipert=1,nbranch
         if (gkk_flag(ipert,ipert,ik_this_proc,isppol,elph_ds%qirredtofull(iqptirred)) == -1) then
           write (msg,'(a,i5,1x,i5,1x,i5,1x,i5,a,a)')&
&           ' gkk element',ipert,ikpt_phon,isppol,iqptirred,' was not found by symmetry operations ',&
&           ' on the irreducible perturbations and qpoints given'
           MSG_ERROR(msg)
         end if
       end do !ipert
     end do !ik_this_proc
   end do !isppol
 end do !iqptirred

 call wrtout(std_out,'read_gkk : done completing the perturbations (and checked!)','COLL')

!MG save phonon frequencies, ph-linewidths and lambda(q,n) values before Fourier interpolation
 ABI_ALLOCATE(elph_ds%qgrid_data,(elph_ds%nqptirred,nbranch,nsppol,3))

 do iqptirred=1,elph_ds%nqptirred
   elph_ds%qgrid_data(iqptirred,:,:,:)=qdata_tmp(iqptirred,:,:,:)
 end do

 ABI_DEALLOCATE(qdata_tmp)

end subroutine read_gkk
!!***
