!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_iogkk
!! NAME
!!  m_iogkk
!!
!! FUNCTION
!!  IO routines for GKK files
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (MVer)
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

module m_iogkk

 use defs_basis
 use defs_datatypes
 use defs_elphon
 use m_errors
 use m_abicore
 use m_xmpi
 use m_krank
 use m_hdr

 use defs_abitypes,     only : MPI_type
 use m_numeric_tools,   only : wrap2_pmhalf
 use m_io_tools,        only : open_file, get_unit
 use m_symtk,           only : mati3inv, littlegroup_q
 use m_geometry,        only : phdispl_cart2red, littlegroup_pert
 use m_crystal,         only : crystal_t
 use m_ifc,             only : ifc_type
 use m_dynmat,          only : d2sym3

 implicit none

 private
!!***

 public :: read_gkk
 public :: outgkk
 public :: read_el_veloc
!!***

contains
!!***

!!****f* m_iogkk/read_gkk
!!
!! NAME
!! read_gkk
!!
!! FUNCTION
!! This routine reads in elphon matrix elements and completes them
!! using the appropriate symmetries
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
!!      completeperts,get_rank,hdr_bcast,hdr_fort_read,hdr_free,ifc_fourq
!!      littlegroup_pert,littlegroup_q,mati3inv,normsq_gkq,phdispl_cart2red
!!      prt_gkk_yambo,wrap2_pmhalf,wrtout,xmpi_bcast
!!
!! SOURCE

subroutine read_gkk(elph_ds,Cryst,ifc,Bst,FSfullpqtofull,gkk_flag,n1wf,nband,ep_prt_yambo,unitgkk)

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

 ABI_MALLOC_OR_DIE(h1_mat_el,(2, nFSband**2, nbranch, elph_ds%k_phon%my_nkpt, nsppol), ierr)
 h1_mat_el= zero

 ABI_MALLOC_OR_DIE(h1_mat_el_sq,(2, nFSband**2, nbranch**2,elph_ds%k_phon%my_nkpt, nsppol), ierr)
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
   ABI_MALLOC_OR_DIE(gkk_qpt_tmp,(2,elph_ds%ngkkband**2,nbranch**2,nsppol), ierr)
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

   call ifc%fourq(cryst,qptirred_local(:,iqptirred),phfrq_tmp,displ_cart,out_eigvec=eigvec)

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
       symrankkpt = elph_ds%k_phon%krank%get_rank (hdr1%kptns(:,ikpt1)-qptirred_local(:,iqptirred))
       ikpt1_phon = elph_ds%k_phon%krank%invrank(symrankkpt)
       if (ikpt1_phon < 0) then
         write (msg,'(a,3es16.6,a)')' irred k ',hdr1%kptns(:,ikpt1),' was not found in full grid'
         MSG_ERROR(msg)
       end if
!      find correspondence between this kpt_phon and the others
!      symrc1 conserves perturbation as well as qpoint
!      add to FSirrtok list
       do isym1=1,nsym1
         do itim1=0,qtimrev
           timsign=one-two*itim1
           kpt(:) = timsign*matmul(symrc1(:,:,isym1), elph_ds%k_phon%kpt(:,ikpt1_phon))

           symrankkpt = elph_ds%k_phon%krank%get_rank (kpt)
           jkpt_phon = elph_ds%k_phon%krank%invrank(symrankkpt)

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

!!****f* m_iogkk/outgkk
!! NAME
!! outgkk
!!
!! FUNCTION
!! output gkk file for one perturbation (used for elphon calculations in anaddb)
!!
!! INPUTS
!!  bantot0 = total number of bands for all kpoints
!!  bantot1 = total number of matrix elements for 1st order eigenvalues
!!  eigen0 = GS eigenvalues
!!  eigen1 = response function 1st order eigenvalue matrix
!!  hdr0 = GS header
!!  hdr1 = RF header
!!  mpi_enreg=information about MPI parallelization
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!      hdr_fort_write,wrtout
!!
!! SOURCE

subroutine outgkk(bantot0,bantot1,outfile,eigen0,eigen1,hdr0,hdr1,mpi_enreg,phasecg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bantot0,bantot1
 character(len=fnlen),intent(in) :: outfile
 type(MPI_type),intent(in) :: mpi_enreg
 type(hdr_type),intent(inout) :: hdr0,hdr1
!arrays
 real(dp),intent(in) :: eigen0(bantot0),eigen1(2*bantot1)
 real(dp),intent(in) :: phasecg(2,bantot1)

!Local variables-------------------------------
!scalars
 integer :: fform,iband,ikpt,isppol,me,ntot,unitout
 integer :: iband_off, mband, ierr
 character(len=500) :: msg
 real(dp), allocatable :: tmpeig(:)

! *************************************************************************

!only master should be writing to disk
!Init me
 me=mpi_enreg%me_kpt
 if (me /= 0) return

 call wrtout(std_out,' writing gkk file: '//outfile,"COLL")

!initializations
 fform = 42
 ntot = 1

!open gkk file
 if (open_file(outfile, msg, newunit=unitout, form='unformatted', status='unknown', action="write") /= 0) then
   MSG_ERROR(msg)
 end if

!output GS header
 call hdr_fort_write(hdr0, unitout, fform, ierr)
 ABI_CHECK(ierr == 0 , "hdr_fort_write returned ierr != 0")

!output GS eigenvalues
 iband=0
 do isppol=1,hdr0%nsppol
   do ikpt=1,hdr0%nkpt
     write (unitout) eigen0(iband+1:iband+hdr0%nband(ikpt))
     iband=iband+hdr0%nband(ikpt)
   end do
 end do

!output number of gkk in this file (1)
 write (unitout) ntot

!output RF header
 call hdr_fort_write(hdr1, unitout, fform, ierr)
 ABI_CHECK(ierr == 0 , "hdr_fort_write returned ierr != 0")

!output RF eigenvalues
 mband = maxval(hdr1%nband(:))
 ABI_ALLOCATE(tmpeig,(2*mband**2))
 iband_off = 0
 tmpeig(1) = phasecg(1, 1)
 do isppol = 1, hdr1%nsppol
   do ikpt = 1, hdr1%nkpt
     tmpeig = zero
     do iband = 1, hdr1%nband(ikpt)**2
       tmpeig (2*(iband-1)+1) = eigen1(2*(iband_off+iband-1)+1)
       tmpeig (2*(iband-1)+2) = eigen1(2*(iband_off+iband-1)+2)
     end do
     write (unitout) tmpeig(1:2*hdr1%nband(ikpt)**2)
     iband_off = iband_off + hdr1%nband(ikpt)**2
   end do
 end do
 ABI_DEALLOCATE(tmpeig)

!close gkk file
 close (unitout)

end subroutine outgkk
!!***

!!****f* m_iogkk/prt_gkk_yambo
!!
!! NAME
!! prt_gkk_yambo
!!
!! FUNCTION
!! This routine outputs el-phon related quantities for the yambo code at 1
!!   q-point
!!
!! INPUTS
!!  displ_cart = phonon displacement vectors for this q-point in Cartesian coordinates.
!!  displ_red = phonon displacement vectors for this q-point, in reduced coordinates
!!  elph_ds = datastructure containing elphon matrix elements
!!  h1_mat_el = matrix elements of first order hamiltonian for present q-point,
!!     all perturbations
!!  iqptfull = index of present q-point in full array of q-points
!!  irredpert = index of irreducible perturbation (atom displaced)
!!  natom = number of atoms
!!  phfrq = phonon frequencies at present q-point
!!  qptn = q-point we will print for
!!
!! OUTPUT
!!  only writes to a file
!!
!! NOTES
!!
!! PARENTS
!!      read_gkk
!!
!! CHILDREN
!!
!! SOURCE

subroutine prt_gkk_yambo(displ_cart,displ_red,kpt_phon,h1_mat_el,iqpt,&
&       natom,nFSband,nkpt_phon,phfrq,qptn)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,iqpt
 integer,intent(in) :: nFSband,nkpt_phon
 !arrays
 real(dp),intent(in) :: kpt_phon(3,nkpt_phon)
 real(dp),intent(in) :: h1_mat_el(2,nFSband*nFSband,3*natom,nkpt_phon,1)
 real(dp),intent(in) :: phfrq(3*natom)
 real(dp),intent(in) :: displ_cart(2,3*natom,3*natom)
 real(dp),intent(in) :: displ_red(2,3*natom,3*natom)
 real(dp),intent(in) :: qptn(3)

!Local variables-------------------------------
 !scalars
 integer, save :: firsttime=1
 integer :: outunit,ikpt,imode,iband,ibandp,iatom,idir,ibandindex
 integer :: jmode, outunit2, outunit3
 !arrays
 real(dp) :: gkk_mode_dep(2)
! *************************************************************************

!if first time round:
 if (firsttime==1) then
   firsttime=0
!  squash file
   outunit=get_unit()
   open (unit=outunit,file="yambo_elphon_data",status="REPLACE")
   outunit2=get_unit()
   open (unit=outunit2,file="yambo_elphon_gkk_bymode",status="replace")
   outunit3=get_unit()
   open (unit=outunit3,file="yambo_elphon_gkksqtw_bymode",status="replace")

!  write dimensions
   write (outunit,'(a,I6)') 'number of el atoms ', natom
   write (outunit2,'(a,I6)') 'number of el atoms ', natom
   write (outunit3,'(a,I6)') 'number of el atoms ', natom
   write (outunit,'(a,I6)') 'number of ph modes ', 3*natom
   write (outunit2,'(a,I6)') 'number of ph modes ', 3*natom
   write (outunit3,'(a,I6)') 'number of ph modes ', 3*natom
   write (outunit,'(a,I6)') 'number of el bands ', nFSband
   write (outunit2,'(a,I6)') 'number of el bands ', nFSband
   write (outunit3,'(a,I6)') 'number of el bands ', nFSband

!  write k-points
   write (outunit,'(a,I6)') 'number of k-points ', nkpt_phon
   write (outunit2,'(a,I6)') 'number of k-points ', nkpt_phon
   write (outunit3,'(a,I6)') 'number of k-points ', nkpt_phon
   do ikpt=1,nkpt_phon
     write (outunit,'(a,I6,3E20.10)') 'reduced coord kpoint no ', ikpt, kpt_phon(:,ikpt)
     write (outunit2,'(a,I6,3E20.10)') 'reduced coord kpoint no ', ikpt, kpt_phon(:,ikpt)
     write (outunit3,'(a,I6,3E20.10)') 'reduced coord kpoint no ', ikpt, kpt_phon(:,ikpt)
   end do

!  band energies are not accessible this deep in the code: simpler to get them
!  from elsewhere

   close (outunit)
   close (outunit2)
   close (outunit3)
 end if ! first time round

!open file
 outunit=get_unit()
 open (unit=outunit,file="yambo_elphon_data",status="unknown",position="append")

!qpoint
 write (outunit,'(a,I6,3E20.10)') 'reduced coord qpoint no ', iqpt, qptn(:)

!frequencies
 do imode=1,3*natom
   write (outunit,'(a,I6,3E20.10)') 'phonon freq no ', imode, phfrq(imode)
 end do

!displacement vector
 do imode=1,3*natom
   write (outunit,'(a,I6,3E20.10)') 'phonon displ vec no ', imode
   do iatom=1,natom
     write (outunit,'(3(2E20.10,2x))') displ_cart(:,(iatom-1)*3+1:iatom*3,imode)
   end do
 end do

!the beef: matrix elements of the first order hamiltonian for displacement of
!all atoms along all reduced directions
 write (outunit,'(a)') ' matrix elements of all perturbations for this q-point'
 do ikpt=1,nkpt_phon
   write (outunit,'(a,I6)') ' kpoint ', ikpt
   imode=0
   do iatom=1,natom
     do idir=1,3
       imode=imode+1
       write (outunit,'(a,I6,I6)') ' atom, direction = ', iatom,idir
       ibandindex=0
       do iband=1,nFSband
         do ibandp=1,nFSband
           ibandindex=ibandindex+1
           write (outunit,'(a,I6,I6,2E20.10)') ' mat el for n,np ', iband,ibandp,&
&           h1_mat_el(:,ibandindex,imode,ikpt,1)
         end do !bandp
       end do !band
     end do !dir
   end do !atom
 end do

!blank line
 write (outunit,*)
 close (outunit)

 outunit2=get_unit()
 open (unit=outunit2,file="yambo_elphon_gkk_bymode",status="unknown",position="append")
 outunit3=get_unit()
 open (unit=outunit3,file="yambo_elphon_gkksqtw_bymode",status="unknown",position="append")

!qpoint
 write (outunit2,'(a,I6,3E20.10)') 'reduced coord qpoint no ', iqpt, qptn(:)
 write (outunit3,'(a,I6,3E20.10)') 'reduced coord qpoint no ', iqpt, qptn(:)

!print out mode-dependent matrix elements
 write (outunit2,'(a)') ' matrix elements of all phonon modes for this q-point'
 write (outunit3,'(a)') ' 1/w**1/2 times matrix elements of all phonon modes for this q-point'
 do ikpt=1,nkpt_phon
   write (outunit2,'(a,I6)') ' kpoint ', ikpt
   write (outunit3,'(a,I6)') ' kpoint ', ikpt
   ibandindex=0
   do iband=1,nFSband
     do ibandp=1,nFSband
       ibandindex=ibandindex+1
       write (outunit2,'(a,I6,I6)') ' el bands n,np ', iband,ibandp
       write (outunit3,'(a,I6,I6)') ' el bands n,np ', iband,ibandp
       do imode=1,3*natom
!        gkk_mode_dep = cg_zdotc(3*natom,displ_red(:,:,imode),h1_mat_el(:,ibandindex,:,ikpt,1))
         gkk_mode_dep = zero
         do jmode=1,3*natom
           gkk_mode_dep(1) = gkk_mode_dep(1) &
&           + displ_red(1,jmode,imode)*h1_mat_el(1,ibandindex,jmode,ikpt,1) &
&           + displ_red(2,jmode,imode)*h1_mat_el(2,ibandindex,jmode,ikpt,1)
           gkk_mode_dep(2) = gkk_mode_dep(2) &
&           + displ_red(1,jmode,imode)*h1_mat_el(2,ibandindex,jmode,ikpt,1) &
&           - displ_red(2,jmode,imode)*h1_mat_el(1,ibandindex,jmode,ikpt,1)
         end do
         write (outunit2,'(a,I6,2E20.10)') ' mat el for phonon mode num = ', imode, gkk_mode_dep
         write (outunit3,'(a,I6,2E20.10)') ' 1/w**1/2 * mat el for phonon mode num = ', &
&         imode, gkk_mode_dep/sqrt(two*abs(phfrq(imode))+tol10)
       end do !imode
     end do !bandp
   end do !band
 end do
!blank line
 write (outunit2,*)
 write (outunit3,*)

 close (outunit2)
 close (outunit3)

end subroutine prt_gkk_yambo
!!***

!!****f* m_iogkk/read_el_veloc
!!
!! NAME
!! read_el_veloc
!!
!! FUNCTION
!! This routine reads the velocities of the electronic GS
!! for all kpts and bands
!! then maps them into the FS kpt states
!!
!! COPYRIGHT
!! Copyright (C) 2002-2019 ABINIT group (JPCroc) based on conducti
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nkpt_in = number of kpoints according to parent routine
!! nband_in = number of bands according to parent routine
!! nsppol_in = number of spin polarizations
!!
!! OUTPUT
!! el_veloc(nkpt_in,nband_in,3)
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      destroy_kptrank,get_rank,hdr_free,inpgkk,mkkptrank
!!
!! SOURCE

subroutine read_el_veloc(nband_in,nkpt_in,kpt_in,nsppol_in,elph_tr_ds)

!Arguments -----------------------------------
!scalars
 integer, intent(in) :: nband_in,nkpt_in,nsppol_in
 type(elph_tr_type), intent(inout) :: elph_tr_ds
 real(dp), intent(in) :: kpt_in(3,nkpt_in)

!Local variables-------------------------------
!scalars
 integer :: bd2tot_index
 integer :: iband,ii,ikpt, ikpt_ddk
 integer :: isppol,l1,mband
 integer :: bantot1
 integer :: unit_ddk
 integer :: symrankkpt
 character(len=fnlen) :: filnam1,filnam2,filnam3
 character(len=500) :: msg
 type(hdr_type) :: hdr1
 type(krank_t) :: krank
!arrays
 real(dp) :: im_el_veloc(3)
 real(dp),allocatable :: eig1_k(:,:)
 real(dp),allocatable :: eigen11(:),eigen12(:),eigen13(:)

! *********************************************************************************

!Read data file name
!TODO: this should be standardized and read in anaddb always, not
!conditionally. Otherwise when new files are added to the anaddb files
!file...  Catastrophe!

 write(std_out,*)'enter read_el_veloc '

!Read data file
 if (open_file(elph_tr_ds%ddkfilename,msg,newunit=unit_ddk,form='formatted') /= 0) then
   MSG_ERROR(msg)
 end if

 rewind(unit_ddk)
 read(unit_ddk,'(a)')filnam1       ! first ddk file
 read(unit_ddk,'(a)')filnam2       ! second ddk file
 read(unit_ddk,'(a)')filnam3       ! third ddk file
 close (unit_ddk)

 bantot1 = 2*nband_in**2*nkpt_in*nsppol_in

 call inpgkk(eigen11,filnam1,hdr1)
 call hdr_free(hdr1)

 call inpgkk(eigen12,filnam2,hdr1)
 call hdr_free(hdr1)

!we use the hdr1 from the last call - should add some consistency
!testing here, we are trusting users not to mix different ddk files...
 call inpgkk(eigen13,filnam3,hdr1)

!Extract info from the header
 if(hdr1%nsppol /= nsppol_in) then
   MSG_ERROR('nsspol /= input nsppol')
 end if

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(hdr1%nband(1:hdr1%nkpt))
 if (mband /= nband_in) then
   MSG_ERROR('nband_in input to read_el_veloc is inconsistent with mband')
 end if

 write(std_out,*)
 write(std_out,*)                     'readings from read_el_veloc header'
 write(std_out,'(a,i8)')              ' natom                =',hdr1%natom
 write(std_out,'(a,3i8)')             ' nkpt,nband_in,mband  =',hdr1%nkpt,nband_in,mband
 write(std_out,'(a, f10.5,a)' )      ' ecut                 =',hdr1%ecut,' Ha'
 write(std_out,'(a,e15.5,a,e15.5,a)' )' fermie               =',hdr1%fermie,' Ha ',hdr1%fermie*Ha_eV,' eV'

 ABI_ALLOCATE(eig1_k,(2*nband_in**2,3))
 bd2tot_index = 0
 elph_tr_ds%el_veloc=zero

!need correspondence between the DDK kpoints and the kpt_phon
 krank = krank_new(hdr1%nkpt, hdr1%kptns)

 do isppol=1,nsppol_in
   im_el_veloc(:)=zero
   do ikpt=1,nkpt_in
    symrankkpt = krank%get_rank (kpt_in(:,ikpt))
     ikpt_ddk = krank%invrank(symrankkpt)
     if (ikpt_ddk == -1) then
       write(std_out,*)'read_el_veloc ******** error in correspondence between ddk and gkk kpoint sets'
       write(std_out,*)' kpt sets in gkk and ddk files must agree.'
       MSG_ERROR("Aborting now")
     end if
     bd2tot_index=2*nband_in**2*(ikpt_ddk-1)

!    first derivative eigenvalues for k-point
     eig1_k(:,1)=eigen11(1+bd2tot_index:2*nband_in**2+bd2tot_index)
     eig1_k(:,2)=eigen12(1+bd2tot_index:2*nband_in**2+bd2tot_index)
     eig1_k(:,3)=eigen13(1+bd2tot_index:2*nband_in**2+bd2tot_index)

!    turn el_veloc to cartesian coordinates
     do iband=1,nband_in
       do l1=1,3
         do ii=1,3
           elph_tr_ds%el_veloc(ikpt,iband,l1,isppol)=elph_tr_ds%el_veloc(ikpt,iband,l1,isppol)+&
&           hdr1%rprimd(l1,ii)*eig1_k(2*iband-1+(iband-1)*2*nband_in,ii)/two_pi
           im_el_veloc(l1)=im_el_veloc(l1)+&
&           hdr1%rprimd(l1,ii)*eig1_k(2*iband+(iband-1)*2*nband_in,ii)/two_pi
         end do
       end do ! l1
     end do
   end do
 end do ! end isppol

 call krank%free()
 ABI_DEALLOCATE(eig1_k)
 ABI_DEALLOCATE(eigen11)
 ABI_DEALLOCATE(eigen12)
 ABI_DEALLOCATE(eigen13)

 call hdr_free(hdr1)

 write(std_out,*)'out of read_el_veloc '

end subroutine read_el_veloc
!!***

!!****f* m_iogkk/inpgkk
!! NAME
!! inpgkk
!!
!! FUNCTION
!! read in gkk file and return eigenvalue matrix
!! Only works for a single gkk matrix (1 perturbation and qpoint) in the file
!! like the files produced by outgkk
!!
!! INPUTS
!!
!!  filegkk= filename
!!
!! OUTPUT
!!  eigen1 = response function 1st order eigenvalue matrix
!!
!! PARENTS
!!      read_el_veloc
!!
!! CHILDREN
!!      hdr_fort_read,hdr_free,wrtout
!!
!! SOURCE

subroutine inpgkk(eigen1,filegkk,hdr1)

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: filegkk
 type(hdr_type), intent(out) :: hdr1
!arrays
 real(dp),allocatable,intent(out) :: eigen1(:)

!Local variables-------------------------------
!scalars
 integer :: bantot1
 integer :: isppol, ikpt, mband, ikb
 integer :: unitgkk, fform, ierr, n1wf, i1wf
 type(hdr_type) :: hdr0
 real(dp), allocatable :: eigen(:)
 character(len=500) :: message

! *************************************************************************

 if (open_file(filegkk,message,newunit=unitgkk,form='unformatted',status='old') /= 0) then
   MSG_ERROR(message)
 end if


!read in header of GS file and eigenvalues
 call hdr_fort_read(hdr0, unitgkk, fform)
 ABI_CHECK(fform /= 0, "hdr_fort_read returned fform == 0")

 mband = maxval(hdr0%nband(:))
 ABI_ALLOCATE(eigen,(mband))
 call wrtout(std_out,'inpgkk : try to reread GS eigenvalues','COLL')

 do isppol=1,hdr0%nsppol
   do ikpt=1,hdr0%nkpt
     read (unitgkk,IOSTAT=ierr) eigen(1:hdr0%nband(ikpt))
     ABI_CHECK(ierr==0,'reading eigen from gkk file')
   end do
 end do

 read(unitgkk,IOSTAT=ierr) n1wf
 ABI_CHECK(ierr==0,"reading n1wf from gkk file")

 ABI_DEALLOCATE(eigen)
 call hdr_free(hdr0)

 if (n1wf > 1) then
   write(message,'(3a)')&
&   'several 1wf records were found in the file,',ch10, &
&   'which is not allowed for reading with this routine'
   MSG_ERROR(message)
 end if

!read in header of 1WF file
 call hdr_fort_read(hdr1, unitgkk, fform)
 if (fform == 0) then
   write(message,'(a,i0,a)')' 1WF header number ',i1wf,' was mis-read. fform == 0'
   MSG_ERROR(message)
 end if

 bantot1 = 2*hdr1%nsppol*hdr1%nkpt*mband**2
 ABI_ALLOCATE(eigen1, (bantot1))


!retrieve 1WF <psi_k+q | H | psi_k> from gkk file and echo to output
 ikb = 0
 do isppol=1,hdr1%nsppol
   do ikpt=1,hdr1%nkpt
     read (unitgkk,IOSTAT=ierr) eigen1(ikb+1:ikb+2*hdr1%nband(ikpt)**2)
     ikb = ikb + 2*hdr1%nband(ikpt)**2
     if (ierr /= 0) then
       write(message,'(a,2i0)')'reading eigen1 from gkk file, spin, kpt_idx',isppol,ikpt
       MSG_ERROR(message)
     end if
   end do
 end do

 close(unitgkk)

end subroutine inpgkk
!!***

!!****f* m_iogkk/completeperts
!!
!! NAME
!! completeperts
!!
!! FUNCTION
!!  Complete perturbations wrt atoms and reduced directions
!!  for a fixed qpoint. Normally there is a test in read_gkk which guarantees
!!  that enough irreducible perturbations are present to generate everything.
!!  h1_mat_el is first squared, making a (ipert,jpert) matrix which has the same
!!  symmetry properties as the dynamical matrix.
!!
!! INPUTS
!!  Cryst<crystal_t>=Info on the unit cell and symmetries.
!!   nbranch=Number of phonon branches.
!!   nFSband=Number of bands in H1 matrix elements.
!!   nkpt=Number of k-points in matrix elements.
!!   nsppol=Number of independent spin polarizations.
!!   gkk_flag = flags for presence of gkk matrix elements
!!   h1_mat_el = irreducible matrix elements to be completed and squared
!!   qpt = qpoint
!!   symq = flags for symmetry elements conserving the present qpoint
!!   tnons = translation vectors associated with symops
!!
!! OUTPUT
!!   h1_mat_el_sq = irreducible matrix elements squared and completed
!!   gkk_flag = changed on output
!!
!! PARENTS
!!      read_gkk
!!
!! CHILDREN
!!      d2sym3
!!
!! SOURCE

subroutine completeperts(Cryst,nbranch,nFSband,nkpt,nsppol,gkk_flag,h1_mat_el,h1_mat_el_sq,&
&   qpt,symq,qtimrev)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: qtimrev,nbranch,nFSband,nkpt,nsppol
 type(crystal_t),intent(in) :: Cryst
!arrays
 integer,intent(in) :: symq(4,2,Cryst%nsym)
 integer,intent(inout) :: gkk_flag(nbranch,nbranch,nkpt,nsppol)
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(in) :: h1_mat_el(2,nFSband**2,nbranch,nkpt,nsppol)
 real(dp),intent(out) :: h1_mat_el_sq(2,nFSband**2,nbranch**2,nkpt,nsppol)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,iatom1,iatom2,ibb,idir1,idir2,ipert1,ipert2,isppol,mpert,natom
 real(dp) :: im1,im2,re1,re2,res
 character(len=500) :: msg
!arrays
 integer,allocatable :: tmpflg(:,:,:,:)
 real(dp),allocatable :: tmpval(:,:,:,:,:)

! *************************************************************************

!WARNING! Stupid patch in d2sym3 imposes these matrices to have size natom+2
 natom = Cryst%natom
 mpert = natom+2

 ABI_ALLOCATE(tmpflg,(3,mpert,3,mpert))
 ABI_ALLOCATE(tmpval,(2,3,mpert,3,mpert))

 h1_mat_el_sq = zero

 write (std_out,*) ' completeperts: shape(h1_mat_el_sq) = ', shape(h1_mat_el_sq)

 do isppol=1,nsppol
   write(std_out,*)'completeperts: isppol = ', isppol
!
   do ikpt_phon=1,nkpt
     do ibb=1,nFSband**2
!
       tmpval = zero
       tmpflg = 0
       ! for a fixed k (q) band and sppol construct the gamma matrix for (3 natom)^2 perturbation pairs
       do iatom1=1,natom
         do idir1=1,3
           ipert1 = (iatom1-1)*3+idir1
           if (gkk_flag(ipert1,ipert1,ikpt_phon,isppol) < 0) cycle
           re1 = h1_mat_el(1,ibb,ipert1,ikpt_phon,isppol)
           im1 = h1_mat_el(2,ibb,ipert1,ikpt_phon,isppol)

           do iatom2=1,natom
             do idir2=1,3
               ipert2 = (iatom2-1)*3+idir2
               if (gkk_flag(ipert2,ipert2,ikpt_phon,isppol) < 0) cycle
               tmpflg(idir1,iatom1,idir2,iatom2) = 1
               re2 = h1_mat_el(1,ibb,ipert2,ikpt_phon,isppol)
               im2 = h1_mat_el(2,ibb,ipert2,ikpt_phon,isppol)
!
!              conjg(h1_mat_el_2) * h1_mat_el_1
               res =  re1*re2 + im1*im2
               tmpval(1,idir1,iatom1,idir2,iatom2) =  res
               res =  re1*im2 - im1*re2
               tmpval(2,idir1,iatom1,idir2,iatom2) = res

             end do !idir2
           end do !iatom2
         end do !idir1
       end do !iatom1

       ! matrix is symmetrized like a dynamical matrix. No change of band or k
       !  in here. This should be checked (if we have to restrict further the symmetry operations)
       call d2sym3(tmpflg,tmpval,Cryst%indsym,mpert,natom,Cryst%nsym,qpt,symq,Cryst%symrec,Cryst%symrel,qtimrev,1)
       if (sum(tmpflg(:,1:natom,:,1:natom)) /= 3*natom*3*natom) then
         write(msg,'(3a,4i0)')&
&         'A perturbation is missing after completion with d2sym3',ch10,&
&         'tmpflg, ikpt_phon, isppol: ',tmpflg,ikpt_phon,isppol
         MSG_ERROR(msg)
       end if
!
!      Save values for calculation of |gkk|^2
       do iatom1=1,natom
         do idir1=1,3
           ipert1 = (iatom1-1)*3+idir1
           do iatom2=1,natom
             do idir2=1,3
!
!              mjv 29/10/2007 ipert2 now contains the composite index ip1*nperts+ip2
               ipert2 = (iatom2-1)*3 + idir2 + (ipert1-1)*3*natom
               h1_mat_el_sq(1,ibb,ipert2,ikpt_phon,isppol) = pi*tmpval(1,idir2,iatom2,idir1,iatom1)
               h1_mat_el_sq(2,ibb,ipert2,ikpt_phon,isppol) = pi*tmpval(2,idir2,iatom2,idir1,iatom1)
             end do
           end do
         end do
       end do
!
     end do !end ibb band dos
!
!    Set flags.
     do ipert1=1,3*natom
       do ipert2=1,3*natom
         if (gkk_flag(ipert2,ipert1,ikpt_phon,isppol) < 0) gkk_flag(ipert2,ipert1,ikpt_phon,isppol) = 1
       end do
     end do

   end do !end kpt_phon do
 end do !end sppol do

 ABI_DEALLOCATE(tmpflg)
 ABI_DEALLOCATE(tmpval)

end subroutine completeperts
!!***

!!****f* ABINIT/normsq_gkq
!!
!! NAME
!! normsq_gkq
!!
!! FUNCTION
!! This routine takes the gkq matrix elements for a given qpoint,
!!   does the scalar product with the phonon displacement vector,
!!   squares the gkq matrix elements multiplies by the appropriate weights
!!   and puts them in a uniform (atom,icart) basis
!!
!! INPUTS
!!   displ_red = phonon mode displacement vectors in reduced coordinated.
!!   eigvec = eigenvectors of phonons (to turn to cartesian coord frame)
!!   elph_ds = datastructure with gkk matrix elements
!!   FSfullpqtofull = mapping of k + q to k
!!   h1_mat_el_sq = matrix elements $<psi_{k+q,m} | H^{1} | psi_{k,n}>$ matrix-squared
!!   iqptirred = index of present qpoint
!!   phfrq_tmp = phonon frequencies
!!   qpt_irred = array of qpoint coordinates
!!
!! OUTPUT
!!   elph_ds%gkq filled
!!   qdata(elph_ds%nbranch,elph_ds%nsppol,3) = array containing the phonon frequency, the linewidth and $\lambda_{q,\nu}$.
!!
!! PARENTS
!!      read_gkk
!!
!! CHILDREN
!!      gam_mult_displ,nmsq_gam,nmsq_gam_sumfs,nmsq_pure_gkk
!!      nmsq_pure_gkk_sumfs,wrtout,xmpi_sum,zhpev
!!
!! SOURCE

subroutine normsq_gkq(displ_red,eigvec,elph_ds,FSfullpqtofull,&
&    h1_mat_el_sq,iqptirred,phfrq_tmp,qpt_irred,qdata)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptirred
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(inout) :: &
& h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol)
 real(dp),intent(in) :: phfrq_tmp(elph_ds%nbranch),qpt_irred(3,elph_ds%nqptirred)
 real(dp),intent(out) :: qdata(elph_ds%nbranch,elph_ds%nsppol,3)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,ier,ii,isppol,jbranch,comm
 real(dp) :: lambda_tot
 character(len=500) :: message
!arrays
 real(dp) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp) :: gam_now2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: lambda(elph_ds%nsppol)
 real(dp),allocatable :: matrx(:,:),val(:),vec(:,:,:)
 real(dp),allocatable :: zhpev1(:,:),zhpev2(:)

! *************************************************************************

 DBG_ENTER("COLL")

 accum_mat  = zero
 accum_mat2 = zero
 comm = xmpi_world

 if (elph_ds%ep_scalprod == 1) then
!
   if (elph_ds%ep_keepbands == 0) then
     call wrtout(std_out,' normsq_gkq : calling nmsq_gam_sumFS',"COLL")
     call nmsq_gam_sumFS (accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&
&     h1_mat_el_sq,iqptirred)

   else if (elph_ds%ep_keepbands == 1) then
     call wrtout(std_out,' normsq_gkq : calling nmsq_gam',"COLL")
     call nmsq_gam (accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&
&     h1_mat_el_sq,iqptirred)

   else
     write (message,'(a,i0)')' Wrong value for elph_ds%ep_keepbands = ',elph_ds%ep_keepbands
     MSG_BUG(message)
   end if
!
 else if (elph_ds%ep_scalprod == 0) then  ! Interpolate on the pure "matrix of matrix elements" and do the scalar products later.
!
   if (elph_ds%ep_keepbands == 0) then
     call wrtout(std_out,' normsq_gkq : calling nmsq_pure_gkk_sumFS',"COLL")
     call nmsq_pure_gkk_sumFS (accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,&
&     h1_mat_el_sq,iqptirred)

   else if (elph_ds%ep_keepbands == 1) then
     call wrtout(std_out,' normsq_gkq : calling nmsq_pure_gkk',"COLL")

     call nmsq_pure_gkk (accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,&
&     h1_mat_el_sq,iqptirred)
   else
     write (message,'(a,i0)')' Wrong value for elph_ds%ep_keepbands = ',elph_ds%ep_keepbands
     MSG_BUG(message)
   end if
!
 else
   write (message,'(a,i0)')' Wrong value for elph_ds%ep_scalprod = ',elph_ds%ep_scalprod
   MSG_BUG(message)
 end if
!end if flag for doing scalar product now.


!MG: values without the good prefactor
 accum_mat = accum_mat * elph_ds%occ_factor/elph_ds%k_phon%nkpt

!MG: accum_mat2 contains the line-widhts before the Fourier interpolation
 accum_mat2 = accum_mat2 * elph_ds%occ_factor/elph_ds%k_phon%nkpt

!mpi sum over procs for accum_mat2
 call xmpi_sum (accum_mat, comm, ier)
 call xmpi_sum (accum_mat2, comm, ier)

!MG20060531i
!write e-ph quantities before Fourier interpolation
!save e-ph values in the temporary array qdata that will be copied into elph_ds%qgrid_data

 write (message,'(4a,3es16.6,63a)')ch10,                  &
& ' Phonon linewidths before interpolation ',ch10,        &
& ' Q point = ',qpt_irred(:,iqptirred),ch10,('=',ii=1,60),ch10,&
& ' Mode          Frequency (Ha)  Linewidth (Ha)  Lambda '
 call wrtout(std_out,message,'COLL')

 lambda_tot = zero
 do isppol=1,elph_ds%nsppol
   do ii=1,elph_ds%nbranch
     lambda(isppol)=zero
!    MG: the tolerance factor is somehow arbitrary
     if (abs(phfrq_tmp(ii)) > tol10) lambda(isppol)=accum_mat2(1,ii,ii,isppol)/&
&     (pi*elph_ds%n0(isppol)*phfrq_tmp(ii)**2)
     lambda_tot=lambda_tot+lambda(isppol)
     write(message,'(i8,es20.6,2es16.6)' )ii,phfrq_tmp(ii),accum_mat2(1,ii,ii,isppol),lambda(isppol)
     call wrtout(std_out,message,'COLL')
!    save values
     qdata(ii,isppol,1)=phfrq_tmp(ii)
     qdata(ii,isppol,2)=accum_mat2(1,ii,ii,isppol)
     qdata(ii,isppol,3)=lambda(isppol)
   end do !loop over branch
 end do !loop over sppol

!normalize for number of spins
 lambda_tot = lambda_tot / elph_ds%nsppol

 write(message,'(61a,44x,es16.6,62a)' )('=',ii=1,60),ch10,lambda_tot,ch10,('=',ii=1,60),ch10
 call wrtout(std_out,message,'COLL')
!ENDMG20060531

!immediately calculate linewidths:
 write(std_out,*) 'summed accum_mat = '
 write(std_out,'(3(2E18.6,1x))') accum_mat(:,:,:,1)
 write(std_out,*) 'summed accum_mat2 = '
 write(std_out,'(3(2E18.6,1x))')  (accum_mat2(:,ii,ii,1),ii=1,elph_ds%nbranch)
 write(std_out,*) 'displ_red  = '
 write(std_out,'(3(2E18.6,1x))') displ_red

 if (elph_ds%ep_scalprod == 1) then
   do isppol=1,elph_ds%nsppol
!    Diagonalize gamma matrix at qpoint (complex matrix). Copied from dfpt_phfrq
     ier=0
     ii=1
     ABI_ALLOCATE(matrx,(2,(elph_ds%nbranch*(elph_ds%nbranch+1))/2))
     do i2=1,elph_ds%nbranch
       do i1=1,i2
         matrx(1,ii)=accum_mat2(1,i1,i2,isppol)
         matrx(2,ii)=accum_mat2(2,i1,i2,isppol)
         ii=ii+1
       end do
     end do
     ABI_ALLOCATE(zhpev1,(2,2*elph_ds%nbranch-1))
     ABI_ALLOCATE(zhpev2,(3*elph_ds%nbranch-2))
     ABI_ALLOCATE(val,(elph_ds%nbranch))
     ABI_ALLOCATE(vec,(2,elph_ds%nbranch,elph_ds%nbranch))
     call ZHPEV ('V','U',elph_ds%nbranch,matrx,val,vec,elph_ds%nbranch,zhpev1,zhpev2,ier)

     write (std_out,*) ' normsq_gkq : accumulated eigenvalues isppol ',isppol, ' = '
     write (std_out,'(3E18.6)') val
     ABI_DEALLOCATE(matrx)
     ABI_DEALLOCATE(zhpev1)
     ABI_DEALLOCATE(zhpev2)
     ABI_DEALLOCATE(vec)
     ABI_DEALLOCATE(val)
   end do ! isppol

 else if (elph_ds%ep_scalprod == 0) then


   do isppol=1,elph_ds%nsppol
     call gam_mult_displ(elph_ds%nbranch, displ_red, accum_mat(:,:,:,isppol), gam_now2)

     write (std_out,*) ' normsq_gkq : accumulated eigenvalues isppol ', isppol, ' = '
     write (std_out,'(3(E14.6,1x))') (gam_now2(1,jbranch,jbranch), jbranch=1,elph_ds%nbranch)
     write (std_out,*) ' normsq_gkq : imag part = '
     write (std_out,'(3(E14.6,1x))') (gam_now2(2,jbranch,jbranch), jbranch=1,elph_ds%nbranch)
   end do ! isppol

 end if

 DBG_EXIT("COLL")

end subroutine normsq_gkq
!!***

!!****f* ABINIT/nmsq_gam
!!
!! NAME
!! nmsq_gam
!!
!! FUNCTION
!!  Calculate gamma matrices keeping full dependence on bands
!!  from original h1_mat_el_sq matrix elements (no averaging over
!!  bands near the Fermi surface)
!!
!! INPUTS
!!   displ_red = phonon mode displacement vectors, post-multiplied by gprim matrix
!!     (ie. turned to reduced coordinates)
!!   eigvec = phonon eigenvectors
!!   elph_ds = datastructure with gkk matrix elements
!!   FSfullpqtofull = mapping of k+q to k
!!   kpt_phon = coordinates of kpoints near to FS
!!   h1_mat_el_sq = matrix elements $<psi_{k+q,m} | H^{1} | psi_{k,n}>$ squared
!!   iqptirred = index of present qpoint
!!
!! OUTPUT
!!   accum_mat = matrix for accumulating FS average of gkk (gamma matrix -> linewidths)
!!   accum_mat2 = matrix for accumulating FS average of gamma matrix with good prefactors
!!
!! PARENTS
!!      normsq_gkq
!!
!! CHILDREN
!!      gam_mult_displ,zgemm
!!
!! SOURCE

subroutine nmsq_gam (accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&
&  h1_mat_el_sq,iqptirred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptirred
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(inout) :: &
& h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)

!Local variables-------------------------------
! tmp variables for diagonalization
!scalars
 integer :: ikpt_phon,ikpt_phonq,ib1,ib2,ibeff,ibranch,isppol,ipert1
 integer :: jbranch
 integer :: iqpt_fullbz
 integer :: ik_this_proc
 real(dp) :: sd1,sd2
 character(len=500) :: message
!arrays
 real(dp) :: gkq_1band(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_mat2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: zgemm_tmp_mat(2,elph_ds%nbranch,elph_ds%nbranch)

! *************************************************************************

 if (elph_ds%ep_keepbands == 0) then
   write (message,'(a,i0)')' elph_ds%ep_keepbands should be 1 while is ',elph_ds%ep_keepbands
   MSG_ERROR(message)
 end if

!MG20060603 NOTE:
!accum_mat and accum_mat2 are real, the imaginary part is used for debugging purpose
!accum_mat2 is used to store the phonon-linewidhts before interpolation

 iqpt_fullbz = elph_ds%qirredtofull(iqptirred)
 write(std_out,*) 'nmsq_gam : iqptirred = ', iqptirred

 do isppol=1,elph_ds%nsppol
   do ik_this_proc =1, elph_ds%k_phon%my_nkpt
     ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

     ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)

     do ib1=1,elph_ds%nFSband
       sd1 = elph_ds%k_phon%wtk(ib1,ikpt_phon,isppol) !weights for distance from the fermi surface

       do ib2=1,elph_ds%nFSband
         sd2 = elph_ds%k_phon%wtk(ib2,ikpt_phonq,isppol) !weights for distance from the fermi surface
         ibeff = ib2+elph_ds%nFSband*(ib1-1)

         gkq_1band(:,:,:) = zero

         zgemm_tmp_mat= reshape (h1_mat_el_sq(:,ibeff,:,ik_this_proc,isppol),(/2,elph_ds%nbranch,elph_ds%nbranch/))

         call gam_mult_displ(elph_ds%nbranch, displ_red, zgemm_tmp_mat, tmp_mat2)

!        sum over bands
         do ipert1=1,elph_ds%nbranch
           gkq_1band(1,ipert1,ipert1) = gkq_1band(1,ipert1,ipert1) + tmp_mat2(1,ipert1,ipert1)
         end do

!        summing over k points and bands, still diagonal in jbranch
         accum_mat(:,:,:,isppol) = accum_mat(:,:,:,isppol) + gkq_1band(:,:,:)*sd1*sd2

!        MG20060603 : summing over bands and kpoints with weights to calculate the phonon linewidth
         do jbranch=1,elph_ds%nbranch
           accum_mat2(:,jbranch,jbranch,isppol) = accum_mat2(:,jbranch,jbranch,isppol) + gkq_1band(:,jbranch,jbranch)*sd1*sd2
         end do
!        END MG


!        now turn to cartesian coordinates

!        Final Gamma matrix (hermitian) = E * D_g * E^{+}
!        Where E^{+} is the hermitian conjugate of the eigenvector matrix E
!        And D_g is the diagonal matrix of values of gamma for this qpoint

!        Here gkq_1band is indexed with real phonon modes (not atom+idir)
!        turn gkq_1band to atom+cartesian coordinates (instead of normal coordinates for qpoint)
         tmp_mat2(:,:,:) = zero
         do ibranch =1,elph_ds%nbranch
           do jbranch =1,elph_ds%nbranch
             tmp_mat2(1,ibranch,jbranch) = tmp_mat2(1,ibranch,jbranch) + &
&             eigvec(1,ibranch,jbranch) * gkq_1band(1,jbranch,jbranch)
             tmp_mat2(2,ibranch,jbranch) = tmp_mat2(2,ibranch,jbranch) + &
&             eigvec(2,ibranch,jbranch) * gkq_1band(1,jbranch,jbranch)
           end do
         end do
         gkq_1band(:,:,:) = zero

!        here eigvec is transposed and complexconjugated.
         zgemm_tmp_mat=zero
         call zgemm('n','c',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,cone,&
&         tmp_mat2,elph_ds%nbranch,eigvec,elph_ds%nbranch,czero,zgemm_tmp_mat,elph_ds%nbranch)

         gkq_1band = zgemm_tmp_mat

!        gamma matrix contribution in cartesian coordinates (ie interpolatable form)
         h1_mat_el_sq(:,ibeff,:,ik_this_proc,isppol) = reshape(gkq_1band,(/2,elph_ds%nbranch*elph_ds%nbranch/))

       end do
     end do
!    END loop over bands ib1 ib2

   end do
!  END loop over kpt_phon
 end do
!END loop over nsppol


end subroutine nmsq_gam
!!***

!!****f* ABINIT/nmsq_gam_sumfs
!!
!! NAME
!! nmsq_gam_sumfs
!!
!! FUNCTION
!!  Calculate gamma matrices from original h1_mat_el_sq matrix
!!  elements averaging over bands near the Fermi surface
!!
!! INPUTS
!!   displ_red = phonon mode displacement vectors, post-multiplied by gprim matrix
!!     (ie. turned to reduced coordinates)
!!   eigvec = eigenvectors of phonons (to turn to cartesian coord frame)
!!   elph_ds = datastructure with gkk matrix elements
!!   FSfullpqtofull = mapping of k+q to k
!!   kpt_phon = coordinates of kpoints near to FS
!!   h1_mat_el_sq = matrix elements $<psi_{k+q,m} | H^{1} | psi_{k,n}>$ squared
!!   iqptirred = index of present qpoint
!!
!! OUTPUT
!!   accum_mat = matrix for accumulating FS average of gkk (gamma matrix -> linewidths)
!!   accum_mat2 = matrix for accumulating FS average of gamma matrix with good prefactors
!!
!! PARENTS
!!      normsq_gkq
!!
!! CHILDREN
!!      gam_mult_displ,zgemm
!!
!! SOURCE

subroutine nmsq_gam_sumFS(accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&
&   h1_mat_el_sq,iqptirred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptirred
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(inout) :: &
& h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,ikpt_phonq,ib1,ib2,ibeff,ibranch,ipert1,isppol,jbranch,iqpt_fullbz
 integer :: ik_this_proc
 real(dp) :: sd1,sd2
 character(len=500) :: message
!arrays
 real(dp) :: gkq_sum_bands(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_gkq_sum_bands(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_mat2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),allocatable :: zgemm_tmp_mat(:,:,:)

! *************************************************************************

 if (elph_ds%ep_keepbands /= 0) then
   write (message,'(a,i0)')' elph_ds%ep_keepbands should be 0 in order to average over bands!',elph_ds%ep_keepbands
   MSG_ERROR(message)
 end if

 iqpt_fullbz = elph_ds%qirredtofull(iqptirred)


!MG20060603 NOTE:
!accum_mat and accum_mat2 are real, the imaginary part is used for debugging purpose
!accum_mat2 is used to store the phonon-linewidhts before interpolation

 ABI_ALLOCATE(zgemm_tmp_mat ,(2,elph_ds%nbranch,elph_ds%nbranch))

 do isppol=1,elph_ds%nsppol
   do ik_this_proc =1, elph_ds%k_phon%my_nkpt
     ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

     ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)

     gkq_sum_bands = zero
     tmp_gkq_sum_bands = zero


     do ib1=1,elph_ds%nFSband
!      weights for distance from the fermi surface
       sd1 = elph_ds%k_phon%wtk(ib1,ikpt_phon,isppol)

       do ib2=1,elph_ds%nFSband
!        weights for distance from the fermi surface
         sd2 = elph_ds%k_phon%wtk(ib2,ikpt_phonq,isppol)
         ibeff=ib2+(ib1-1)*elph_ds%nFSband

         zgemm_tmp_mat = reshape(h1_mat_el_sq(:,ibeff,:,isppol,ik_this_proc),(/2,elph_ds%nbranch,elph_ds%nbranch/))

         call gam_mult_displ(elph_ds%nbranch, displ_red, zgemm_tmp_mat, tmp_mat2)

!        sum over bands in gkq_sum_bands
         do ipert1=1,elph_ds%nbranch
           gkq_sum_bands(1,ipert1,ipert1) = gkq_sum_bands(1,ipert1,ipert1) + sd1*sd2*tmp_mat2(1,ipert1,ipert1)
         end do



       end do
     end do
!    END loop over bands

!    summing over k points, still diagonal in jbranch
     accum_mat(:,:,:,isppol) = accum_mat(:,:,:,isppol) + gkq_sum_bands(:,:,:)
     accum_mat2(:,:,:,isppol) = accum_mat2(:,:,:,isppol) + gkq_sum_bands(:,:,:)

!    summed over bands, now turn to cartesian coordinates

!    Final Gamma matrix (hermitian) = E * D_g * E^{+}
!    Where E^{+} is the hermitian conjugate of the eigenvector matrix E
!    And D_g is the diagonal matrix of values of gamma for this qpoint

!    Here gkq_sum_bands is indexed with real phonon modes (not atom+idir)
!    turn gkq_sum_bands to atom+cartesian coordinates (instead of normal coordinates for qpoint)
!    This is not a full matrix multiplication, just vector one, by
!    gkq_sum_bands(1,jbranch,jbranch)
     tmp_mat2(:,:,:) = zero
     do ibranch =1,elph_ds%nbranch
       do jbranch =1,elph_ds%nbranch
         tmp_mat2(1,ibranch,jbranch) = tmp_mat2(1,ibranch,jbranch) + &
&         eigvec(1,ibranch,jbranch) * &
&         gkq_sum_bands(1,jbranch,jbranch)
         tmp_mat2(2,ibranch,jbranch) = tmp_mat2(2,ibranch,jbranch) + &
&         eigvec(2,ibranch,jbranch) * &
&         gkq_sum_bands(1,jbranch,jbranch)
       end do
     end do

!    here eigvec is transposed and complex conjugated.
     zgemm_tmp_mat=zero
     call zgemm('n','c',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,cone,&
&     tmp_mat2,elph_ds%nbranch,eigvec,elph_ds%nbranch,czero,zgemm_tmp_mat,elph_ds%nbranch)

     gkq_sum_bands = zgemm_tmp_mat

!    ! gamma matrix contribution in cartesian coordinates (ie interpolatable form)
!    gamma matrix contribution in reduced coordinates (ie interpolatable form)
     h1_mat_el_sq(:,1,:,ik_this_proc,isppol) = reshape(gkq_sum_bands(:,:,:),(/2,elph_ds%nbranch*elph_ds%nbranch/))

!    accum_mat(:,:,:,isppol) = accum_mat(:,:,:,isppol) + gkq_sum_bands(:,:,:)
   end do
!  END loop over kpt_phon
 end do
!END loop over sppol

 ABI_DEALLOCATE(zgemm_tmp_mat)

end subroutine nmsq_gam_sumFS
!!***


!{\src2tex{textfont=tt}}
!!****f* ABINIT/nmsq_pure_gkk
!!
!! NAME
!! nmsq_pure_gkk
!!
!! FUNCTION
!!  Calculate gamma matrices for pure gkk case, ie when the
!!  scalar product with the displacement vector is done later
!!  Sum over bands is carried out later.
!!
!! INPUTS
!!   displ_red = phonon displacement in reduced coordinates (used to calculate the ph linewidth)
!!   elph_ds = datastructure with gkk matrix elements
!!   FSfullpqtofull = mapping of k+q to k
!!   kpt_phon = coordinates of kpoints near to FS
!!   h1_mat_el_sq = matrix elements $<psi_{k+q,m} | H^{1} | psi_{k,n}>$ squared
!!   iqptirred = index of present qpoint
!!
!! OUTPUT
!!   elph_ds%gkq filled
!!   accum_mat = matrix for accumulating FS average of gkk (gamma matrix -> linewidths)
!!   accum_mat2 = complex array whose real part contains the phonon linewidth
!!
!! PARENTS
!!      normsq_gkq
!!
!! CHILDREN
!!      gam_mult_displ
!!
!! SOURCE

subroutine nmsq_pure_gkk(accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,&
&   h1_mat_el_sq,iqptirred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptirred
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(inout) :: &
& h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,ikpt_phonq,ib1,ib2,ibeff,ipert1,isppol
 integer :: iqpt_fullbz
 integer :: ik_this_proc
 real(dp) :: sd1,sd2
 character(len=500) :: message
!arrays
 real(dp) :: gkq_sum_bands(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_mat2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: zgemm_tmp_mat(2,elph_ds%nbranch,elph_ds%nbranch)

! *************************************************************************

 if (elph_ds%ep_keepbands /= 1) then
   message = ' elph_ds%ep_keepbands should be 1 to keep bands!'
   MSG_ERROR(message)
 end if

 iqpt_fullbz = elph_ds%qirredtofull(iqptirred)

!h1_mat_el_sq is already fine here - nothing to do


!MG20060603 NOTE:
!accum_mat and accum_mat2 are real, the imaginary part is used for debugging purpose
!accum_mat2 is used to store the phonon-linewidhts before interpolation

!MJV 20070525 NOTE:
!in some of the nmsq routines, in particular this one, the work done to
!calculate accum_mat,accum_mat2 is completely superfluous and will be re-done
!on the interpolated values.
!MG uses them for the QPT output, however, so keep it for consistency for the
!moment.

 do isppol=1,elph_ds%nsppol
   do ik_this_proc =1, elph_ds%k_phon%my_nkpt
     ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

     ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)

     gkq_sum_bands(:,:,:) = zero

!    gkq_sum_bands = \sum_{ib1,ib2} \langle k+q \mid H^{(1)}_{q,\tau_i,\alpha_i} \mid k   \rangle
!    \cdot \langle k   \mid H^{(1)}_{q,\tau_j,\alpha_j} \mid k+q \rangle
!    where ibranch -> \tau_i,\alpha_i  and  jbranch -> \tau_j,\alpha_j

     do ib1=1,elph_ds%nFSband

       sd1 = elph_ds%k_phon%wtk(ib1,ikpt_phon,isppol)      !  weights for distance from the fermi surface

       do ib2=1,elph_ds%nFSband

         sd2 = elph_ds%k_phon%wtk(ib2,ikpt_phonq,isppol)  !  weights for distance from the fermi surface
         ibeff = ib2+(ib1-1)*elph_ds%nFSband

         gkq_sum_bands = gkq_sum_bands + &
&         sd1*sd2*reshape(h1_mat_el_sq(:,ibeff,:,ik_this_proc,isppol),(/2,elph_ds%nbranch,elph_ds%nbranch/))

       end do !ib2
     end do !ib1
!    END loops over bands


     accum_mat(:,:,:,isppol) = accum_mat(:,:,:,isppol) + gkq_sum_bands(:,:,:)
   end do
!  END loop over kpt_phon

!  MG20060603
!  do scalar product with the displ_red to calculate the ph lwdth before interpolation (stored in accum_mat2)

   zgemm_tmp_mat = accum_mat(:,:,:,isppol)

   call gam_mult_displ(elph_ds%nbranch, displ_red, zgemm_tmp_mat, tmp_mat2)

   do ipert1=1,elph_ds%nbranch
     accum_mat2(1,ipert1,ipert1,isppol) = accum_mat2(1,ipert1,ipert1,isppol) + tmp_mat2(1,ipert1,ipert1)
   end do

!  ENDMG

 end do ! isppol

end subroutine nmsq_pure_gkk
!!***

!!****f* ABINIT/nmsq_pure_gkk_sumfs
!!
!! NAME
!! nmsq_pure_gkk_sumfs
!!
!! FUNCTION
!!  Calculate gamma matrices for pure gkk case, i.e, when the
!!  scalar product with the displacement vector is done later
!!  Sum over bands is carried out now.
!!
!! INPUTS
!!   displ_red = phonon displacement in reduced coordinates (used to calculate the ph linewidth)
!!   elph_ds = datastructure with gkk matrix elements
!!   FSfullpqtofull = mapping of k+q to k
!!   kpt_phon = coordinates of kpoints near to FS
!!   h1_mat_el_sq = matrix elements $<psi_{k+q,m} | H^{1} | psi_{k,n}>$ squared
!!   iqptirred = index of present qpoint
!!
!! OUTPUT
!!   accum_mat = matrix for accumulating FS average of gkk (gamma matrix -> linewidths)
!!   accum_mat2 = complex array whose real part contains the phonon linewidth
!!
!! PARENTS
!!      normsq_gkq
!!
!! CHILDREN
!!      gam_mult_displ
!!
!! SOURCE

subroutine nmsq_pure_gkk_sumfs(accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,h1_mat_el_sq,iqptirred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptirred
 type(elph_type),intent(in) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(inout) :: &
& h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,ikpt_phonq,ib1,ib2,ibeff,ipert1,isppol,iqpt_fullbz
 integer :: nbranch,nsppol,nFSband,nkpt_phon
 integer :: ik_this_proc
 real(dp) :: sd1,sd2
 !character(len=500) :: message
!arrays
 real(dp) :: gkq_sum_bands(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_mat2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: zgemm_tmp_mat(2,elph_ds%nbranch,elph_ds%nbranch)

! *************************************************************************

 if (elph_ds%ep_keepbands /= 0) then
   MSG_BUG('ep_keepbands should be 0 to average over bands!')
 end if

 nbranch   = elph_ds%nbranch
 nsppol    = elph_ds%nsppol
 nFSband   = elph_ds%nFSband
 nkpt_phon = elph_ds%k_phon%nkpt

 iqpt_fullbz = elph_ds%qirredtofull(iqptirred)

!MG20060603 NOTE:
!accum_mat and accum_mat2 are real, the imaginary part is used for debugging purpose
!accum_mat2 is used to store the phonon-linewidhts before interpolation

 do isppol=1,nsppol
   do ik_this_proc =1, elph_ds%k_phon%my_nkpt
     ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

!
!    The index of k+q in the BZ.
     ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)
!
!    gkq_sum_bands =
!    \sum_{ib1,ib2} <k+q| H^{(1)}_{q,\tau_i,\alpha_i} |k> \cdot <k| H^{(1)}_{q,\tau_j,\alpha_j}|k+q>
!
!    where ibranch = (\tau_i,\alpha_i) and  jbranch = (\tau_j,\alpha_j).
     gkq_sum_bands(:,:,:) = zero

     do ib1=1,nFSband
       sd1 = elph_ds%k_phon%wtk(ib1,ikpt_phon,isppol)      !  weights for distance from the fermi surface

       do ib2=1,nFSband
         sd2 = elph_ds%k_phon%wtk(ib2,ikpt_phonq,isppol)  !  weights for distance from the fermi surface
         ibeff=ib2+(ib1-1)*nFSband

         gkq_sum_bands = gkq_sum_bands + &
&         sd1*sd2* reshape(h1_mat_el_sq(:,ibeff,:,ik_this_proc,isppol),(/2,nbranch,nbranch/))
       end do !ib2
     end do !ib1
!
!    gamma matrix contribution in reduced coordinates (ie interpolatable form)
!    The sum over Fermi surface bands is done here, and fed into (ib1,ib2)=(1,1)
     h1_mat_el_sq(:,1,:,ik_this_proc,isppol) = reshape(gkq_sum_bands,(/2,nbranch**2/))

     accum_mat(:,:,:,isppol) = accum_mat(:,:,:,isppol) + gkq_sum_bands(:,:,:)
   end do ! kpt_phon
 end do ! isppol
!
!MG20060603
!do scalar product wit displ_red to calculate the ph lwdth before interpolation (stored in accum_mat2)
 do isppol=1,nsppol
   zgemm_tmp_mat = accum_mat(:,:,:,isppol)
!
   call gam_mult_displ(nbranch, displ_red, zgemm_tmp_mat, tmp_mat2)

   do ipert1=1,nbranch
     accum_mat2(1,ipert1,ipert1,isppol) = accum_mat2(1,ipert1,ipert1,isppol) + tmp_mat2(1,ipert1,ipert1)
   end do
!
 end do

end subroutine nmsq_pure_gkk_sumfs
!!***

end module m_iogkk
!!***
