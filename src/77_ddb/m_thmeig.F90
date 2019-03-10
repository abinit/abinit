!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_thmeig
!! NAME
!!  m_thmeig
!!
!! FUNCTION
!! Calculate thermal corrections to the eigenvalues.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (PB, XG, GA)
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

module m_thmeig

 use defs_basis
 use m_abicore
 use m_tetrahedron
 use m_errors
 use m_ddb
 use m_ddb_hdr
 use m_xmpi
 use m_sort

 use m_geometry,       only : mkrdim, xred2xcart, metric
 use m_symfind,        only : symfind, symlatt
 use m_symtk,          only : mati3inv, matr3inv, symatm
 use m_crystal,        only : crystal_t
 use m_io_tools,       only : open_file
 use m_dynmat,         only : asria_corr, dfpt_phfrq
 use m_anaddb_dataset, only : anaddb_dataset_type
 use m_pawtab,         only : pawtab_type,pawtab_nullify,pawtab_free
 use m_kpts,           only : getkgrid

 implicit none

 private
!!***

 public :: thmeig
!!***

contains
!!***

!!****f* m_thmeig/thmeig
!! NAME
!! thmeig
!!
!! FUNCTION
!! This routine calculates the thermal corrections to the eigenvalues.
!! The output is this quantity for the input k point.
!!
!! INPUTS
!!  elph_base_name = root filename for outputs
!!  eig2_filnam = name of the eig2 database file
!!  comm=MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine thmeig(inp, ddb, crystal, &
&                 elph_base_name, eig2_filnam, ddbun, iout, &
&                 natom, mpert, msize, d2asr, &
&                 comm)

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: natom
 integer,intent(in) :: mpert,msize
 integer,intent(in) :: comm
 character(len=*),intent(in) :: elph_base_name, eig2_filnam
 integer,intent(in) :: ddbun,iout
 type(crystal_t), intent(inout) :: crystal
 type(anaddb_dataset_type),intent(inout) :: inp
 type(ddb_type),intent(inout) :: ddb
!arrays
 real(dp),intent(inout) :: d2asr(2,3,natom,3,natom)


!Local variables-------------------------------
!scalars
 integer,parameter :: msppol=2,master=0,bcorr0=0
 integer :: msym
 integer :: nkpt,mband,ntypat
 integer :: usepaw,natifc
 integer :: nsym,occopt,nblok2
 integer :: ntemper,telphint,thmflag
 integer :: brav,chksymbreak,found,gqpt,iatom1,iatom2,iband,iblok,iblok2,idir1,idir2,ii,ikpt,ilatt,imod,index
 integer :: iomega,iqpt,iqpt1,iqpt2,iqpt2_previous,iqpt3,iscf_fake,itemper
 integer :: mpert_eig2,msize2,nene,ng2f,nqshft,nsym_new,unit_g2f,nqpt,nqpt_computed,qptopt,rftyp
!integer :: mqpt,nqpt2,option
 integer :: unit_phdos,unitout
 integer :: isym
 integer :: nptsym,use_inversion
 integer :: ierr
 real(dp) :: ucvol
 real(dp) :: g2fsmear,temperinc,tempermin
 real(dp) :: bosein,deltaene,det,domega,enemax,enemin,fact2i,fact2r,factr
 real(dp) :: gaussfactor,gaussprefactor,gaussval,invdet,omega,omega_max,omega_min,qnrm,qptrlen
 real(dp) :: rcvol,tmp,tol,vec1i,vec1r,vec2i,vec2r,veci,vecr,xx
 real(dp) :: tolsym,tolsym8  !new
 character(len=500) :: message
 character(len=fnlen) :: outfile
 type(ddb_type) :: ddb_eig2
 type(ddb_hdr_type) :: ddb_hdr
!arrays
 ! FIXME now these must be allocated
 integer :: ngqpt(9),qptrlatt(3,3),rfelfd(4),rfphon(4),rfstrs(4),vacuum(3)
 integer :: bravais(11)
 integer,allocatable :: typat(:),atifc(:)
 integer,allocatable :: symrel(:,:,:),symrec(:,:,:)
 integer,allocatable :: indsym(:,:,:)
 integer,allocatable :: indqpt(:)
 integer,allocatable :: symafm(:),symafm_new(:)
 integer,allocatable :: carflg_eig2(:,:,:,:)
 integer,allocatable :: ptsymrel(:,:,:),symrel_new(:,:,:)
!integer,allocatable :: symrec_new(:,:,:)
 real(dp) :: rprim(3,3),gprim(3,3),rmet(3,3),gmet(3,3)
 real(dp) :: acell(3)
 real(dp) :: diff_qpt(3)
 real(dp) :: gprimd(3,3),mesh(3,3)
 real(dp) :: qlatt(3,3),qphnrm(3),qpt_search(3,3)
 real(dp) :: rprimd(3,3),shiftq(3,MAX_NSHIFTK),tempqlatt(3)
 real(dp) :: dummy(0),dummy2(0,0)
 real(dp),allocatable :: xcart(:,:),xred(:,:)
 real(dp),allocatable :: amu(:),zion(:)
 real(dp),allocatable :: tnons(:,:)
 real(dp),allocatable :: deigi(:,:), deigr(:,:), multi(:,:), multr(:,:)
 real(dp),allocatable :: dwtermi(:,:), dwtermr(:,:)
 real(dp),allocatable :: slope(:,:,:),thmeigen(:,:,:),zeropoint(:,:,:)
 real(dp),allocatable :: displ(:)
 real(dp),allocatable :: dos_phon(:),dtweightde(:,:),d2cart(:,:)
 real(dp),allocatable :: eigvec(:,:,:,:),eigval(:,:),g2f(:,:,:),intweight(:,:,:)
 real(dp),allocatable :: indtweightde(:,:,:),tmpg2f(:,:,:),tmpphondos(:),total_dos(:),tweight(:,:)
 real(dp),allocatable :: phfreq(:,:)
 real(dp),allocatable :: blkval2(:,:,:,:),blkval2gqpt(:,:,:,:),kpnt(:,:,:)
 real(dp),allocatable :: dedni(:,:,:,:),dednr(:,:,:,:)
 real(dp),allocatable :: eigen_in(:)
 real(dp),allocatable :: qpt_full(:,:),qptnrm(:)
 real(dp),allocatable :: spqpt(:,:),tnons_new(:,:),spinat(:,:)
 real(dp),allocatable :: wghtq(:)

 type(t_tetrahedron) :: tetrahedra
 character(len=80) :: errstr

! *********************************************************************

!DEBUG
! write(std_out,*)'-thmeig : enter '
!call flush(6)
!ENDDEBUG

 ! Only master works for the time being
 if (xmpi_comm_rank(comm) /= master) return

 write(message,'(83a)') ch10,('=',ii=1,80),ch10,&
& ' Computation of the electron-phonon changes to the electronic eigenenergies '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')


!=========================================================================
!0) Initializations
!=========================================================================


 g2fsmear = inp%a2fsmear

 telphint = inp%telphint
 temperinc = inp%temperinc
 tempermin = inp%tempermin
 thmflag = inp%thmflag

 ntemper = inp%ntemper
 natifc = inp%natifc

!Open Derivative DataBase then r/w Derivative DataBase preliminary information.

 write(std_out, '(a)' )  '- thmeig: Initialize the second-order electron-phonon file with name :'
 write(std_out, '(a,a)' )'-         ',trim(eig2_filnam)

 call ddb_hdr_open_read(ddb_hdr, eig2_filnam, ddbun, DDB_VERSION)

 mband = ddb_hdr%mband
 nkpt = ddb_hdr%nkpt
 ntypat = ddb_hdr%ntypat

 msym = ddb_hdr%msym
 nblok2 = ddb_hdr%nblok
 nsym = ddb_hdr%nsym
 occopt = ddb%occopt
 usepaw = ddb_hdr%usepaw

 ABI_ALLOCATE(typat, (natom))
 ABI_ALLOCATE(atifc, (natom))
 ABI_ALLOCATE(zion, (ntypat))
 ABI_ALLOCATE(amu, (ntypat))

 ABI_ALLOCATE(xcart,(3,natom))
 ABI_ALLOCATE(xred,(3,natom))

 ABI_ALLOCATE(symafm, (nsym))
 ABI_ALLOCATE(spinat,(3,natom))

 ABI_ALLOCATE(symrel, (3,3,nsym))
 ABI_ALLOCATE(symrec, (3,3,nsym))
 ABI_ALLOCATE(tnons, (3,nsym))
 ABI_ALLOCATE(indsym, (4,nsym,natom))

 ABI_ALLOCATE(deigi, (mband,nkpt))
 ABI_ALLOCATE(deigr, (mband,nkpt))
 ABI_ALLOCATE(dwtermi, (mband,nkpt))
 ABI_ALLOCATE(dwtermr, (mband,nkpt))
 ABI_ALLOCATE(multi, (mband,nkpt))
 ABI_ALLOCATE(multr, (mband,nkpt))
 ABI_ALLOCATE(slope, (2,mband,nkpt))
 ABI_ALLOCATE(thmeigen, (2,mband,nkpt))
 ABI_ALLOCATE(zeropoint, (2,mband,nkpt))

!At present, only atom-type perturbations are allowed for eig2 type matrix elements.
 mpert_eig2=natom
 msize2=3*mpert_eig2*3*mpert_eig2

 call ddb_malloc(ddb_eig2,msize2,nblok2,natom,ntypat)

 ABI_ALLOCATE(blkval2,(2,msize2,mband,nkpt))
 ABI_ALLOCATE(blkval2gqpt,(2,msize2,mband,nkpt))

 ABI_ALLOCATE(eigvec,(2,3,natom,3*natom))
 ABI_ALLOCATE(phfreq,(3*natom,ddb%nblok))

 atifc = inp%atifc

 !amu = ddb%amu
 amu(:) = ddb_hdr%amu(1:ntypat)
 typat(:) = ddb_hdr%typat(1:natom)
 zion(:) = ddb_hdr%zion(1:ntypat)
 symrel(:,:,1:nsym) = ddb_hdr%symrel(:,:,1:nsym)
 tnons(:,1:nsym) = ddb_hdr%tnons(:,1:nsym)

 xred(:,:) = ddb_hdr%xred(:,:)

 symafm(:) = ddb_hdr%symafm(1:nsym)
 spinat(:,:) = ddb_hdr%spinat(:,1:natom)

 !symrel = ddb_hdr%symrel  ! out
 !tnons = ddb_hdr%tnons  ! out

 !acell = ddb%acell
 !natom = ddb_hdr%natom
 acell = ddb_hdr%acell
 rprim = ddb_hdr%rprim

!Compute different matrices in real and reciprocal space, also
!checks whether ucvol is positive.
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Obtain reciprocal space primitive transl g from inverse trans of r
 call matr3inv(rprim,gprim)

!Generate atom positions in cartesian coordinates
 call xred2xcart(natom,rprimd,xcart,xred)

!Transposed inversion of the symmetry matrices, for use in
!the reciprocal space
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
 end do

!SYMATM generates for all the atoms and all the symmetries, the atom
!on which the referenced one is sent and also the translation bringing
!back this atom to the referenced unit cell
 tolsym8=tol8
 call symatm(indsym,natom,nsym,symrec(:,:,1:nsym),tnons(:,1:nsym),tolsym8,typat,xred)

!Check the correctness of some input parameters,
!and perform small treatment if needed.
 call chkin9(atifc,natifc,natom)

 blkval2gqpt(:,:,:,:)=zero

 ABI_ALLOCATE(carflg_eig2,(3,mpert_eig2,3,mpert_eig2))
 ABI_ALLOCATE(kpnt,(3,nkpt,1))

 ! Copy a bunch of stuff back into crystal (to retain old behavior)
 ! TODO comment these: doesnt make a difference
 crystal%xcart = xcart
 crystal%ucvol = ucvol
 crystal%zion = zion
 crystal%gmet = gmet
 crystal%rmet = rmet
 crystal%nsym = nsym
 crystal%symrel = symrel
 crystal%symrec = symrec
 crystal%tnons = tnons
 crystal%indsym = indsym


!=========================================================================
!1) Take care of the Gamma point for thmflag=3, 5 or 7
!=========================================================================

 if(thmflag==3 .or. thmflag==5 .or. thmflag==7) then
   found=0
   do iblok2=1,nblok2

     call read_blok8(ddb_eig2,iblok2,mband,mpert_eig2,msize2,&
&     nkpt,ddbun,blkval2(:,:,:,:),kpnt(:,:,1))


     qnrm = ddb_eig2%qpt(1,iblok2)*ddb_eig2%qpt(1,iblok2)+ &
&     ddb_eig2%qpt(2,iblok2)*ddb_eig2%qpt(2,iblok2)+ &
&     ddb_eig2%qpt(3,iblok2)*ddb_eig2%qpt(3,iblok2)
     if(qnrm < DDB_QTOL) then
       blkval2gqpt(:,:,:,:) = blkval2(:,:,:,:)
       gqpt=iblok2
       write(std_out,*)'-thmeig: found Gamma point in EIG2 DDB, blok number ',iblok2
       found=1
       exit
     end if
   end do

   if(found==0)then
     write(message,'(a,i3,2a)')&
&     'Was unable to find the blok for Gamma point in EIG2 DDB file, while thmflag= ',thmflag,ch10,&
&     'Action: compute the contribution from Gamma, and merge it in your EIG2 DDB file.'
     MSG_ERROR(message)
   end if

!  Put blkval2gqpt in cartesian coordinates
   call carttransf(ddb_eig2%flg,blkval2gqpt,carflg_eig2,gprimd,gqpt,mband,&
&   mpert_eig2,msize2,natom,nblok2,nkpt,rprimd)

 end if

 close(ddbun)

!=========================================================================
!2) Calculation of dE(n,k)/dn(Q,j) : consider all q and modes
!=========================================================================

 if(thmflag==3 .or. thmflag==4)then


!  Use the first list of q wavevectors
   nqpt=inp%nph1l
   ABI_ALLOCATE(spqpt,(3,nqpt))
   do iqpt=1,inp%nph1l
     spqpt(:,iqpt)=inp%qph1l(:,iqpt)/inp%qnrml1(iqpt)
   end do
   ABI_ALLOCATE(wghtq,(nqpt))
   wghtq(:)=one/nqpt

 else if(thmflag>=5 .and. thmflag<=8)then

!  Generates the q point grid
   ngqpt(1:3)=inp%ngqpt(1:3)
   nqshft=inp%nqshft
   qptrlatt(:,:)=0
   qptrlatt(1,1)=ngqpt(1)
   qptrlatt(2,2)=ngqpt(2)
   qptrlatt(3,3)=ngqpt(3)

   ABI_ALLOCATE(ptsymrel,(3,3,msym))
   ABI_ALLOCATE(symafm_new,(msym))
   ABI_ALLOCATE(symrel_new,(3,3,msym))
   ABI_ALLOCATE(tnons_new,(3,msym))
   if(thmflag==7 .or. thmflag==8) then
!    Re-generate symmetry operations from the lattice and atomic coordinates
     tolsym=tol8
     call symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
     use_inversion=1
     call symfind(0,(/zero,zero,zero/),gprimd,0,msym,natom,0,nptsym,nsym_new,0,0,&
&     ptsymrel,spinat,symafm_new,symrel_new,tnons_new,tolsym,typat,use_inversion,xred)
     write(std_out,*)' thmeig : found ',nsym_new,' symmetries ',ch10
     qptopt=1
   else
     nsym_new=1
     symrel_new(:,:,1)=0 ; symrel_new(1,1,1)=1 ; symrel_new(2,2,1)=1 ; symrel_new(3,3,1)=1
     tnons_new(:,1)=zero
     symafm_new(1)=1
     qptopt=3
   end if

   brav=inp%brav

   if(abs(brav)/=1)then
     message = ' The possibility to have abs(brav)/=1 for thmeig was disabled.'
     MSG_ERROR(message)
   end if

!  Prepare to compute the q-point grid in the ZB or IZB
   iscf_fake=5 ! Need the weights
   chksymbreak=0
   vacuum=0
   shiftq(:,1:nqshft)=inp%q1shft(:,1:nqshft)
!  Compute the final number of q points
   call getkgrid(chksymbreak,0,iscf_fake,dummy2,qptopt,qptrlatt,qptrlen,&
&   nsym_new,0,nqpt,nqshft,nsym_new,rprimd,&
&   shiftq,symafm_new,symrel_new,vacuum,dummy)
   ABI_ALLOCATE(spqpt,(3,nqpt))
   ABI_ALLOCATE(wghtq,(nqpt))
   call getkgrid(chksymbreak,iout,iscf_fake,spqpt,qptopt,qptrlatt,qptrlen,&
&   nsym_new,nqpt,nqpt_computed,nqshft,nsym_new,rprimd,&
&   shiftq,symafm_new,symrel_new,vacuum,wghtq)

   ABI_DEALLOCATE(ptsymrel)
   ABI_DEALLOCATE(symafm_new)
   ABI_DEALLOCATE(symrel_new)
   ABI_DEALLOCATE(tnons_new)

 end if

 call ddb_hdr_free(ddb_hdr)


 write(message,'(a,a)')ch10,' thmeig : list of q wavevectors, with integration weights '
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 do iqpt=1,nqpt
   write(message,'(i6,3es16.6,es20.6)')iqpt,spqpt(:,iqpt),wghtq(iqpt)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end do

 if(.not.allocated(indqpt))then
   ABI_ALLOCATE(indqpt,(nqpt))
 end if
 ABI_ALLOCATE(dedni,(mband,nkpt,3*natom,nqpt))
 ABI_ALLOCATE(dednr,(mband,nkpt,3*natom,nqpt))
 ABI_ALLOCATE(eigen_in,(nqpt))
 ABI_ALLOCATE(qpt_full,(3,nqpt))
 ABI_ALLOCATE(qptnrm,(nqpt))

 dednr(:,:,:,:) = zero
 dedni(:,:,:,:) = zero

!!Prepare the reading of the EIG2 files
 call ddb_hdr_open_read(ddb_hdr, eig2_filnam, ddbun, DDB_VERSION, msym=msym)

 call ddb_hdr_free(ddb_hdr)

!iqpt2 will be the index of the q point bloks inside the EIG2 file
 iqpt2=0

!Sum on all phonon wavevectors and modes
 do iqpt=1,nqpt

!  Finding the target wavevector in DDB file
   qpt_search(:,:)=0.0d0
   qpt_search(:,1)=spqpt(:,iqpt)
   qphnrm(:)=one
   rfphon(1:2)=1
!  NOTE : at present, no LO-TO splitting included !!!
   rfelfd(1:2)=0
   rfstrs(1:2)=0
   rftyp=1

   write(std_out,'(a,3es16.6)' )' Looking for spqpt=',qpt_search(:,1)

   call gtblk9(ddb,iblok,qpt_search,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

   if(iblok==0) then
     write(message,'(a,3es16.6,2a)')&
&     'Was unable to find in DDB file, the blok for point ',spqpt(:,iqpt),ch10,&
&     'Action: compute the contribution from this point, and merge it in your DDB file.'
     MSG_ERROR(message)
   end if

   ABI_ALLOCATE(d2cart,(2,msize))
!  Copy the dynamical matrix in d2cart
   d2cart(:,1:msize)=ddb%val(:,:,iblok)

!  Eventually impose the acoustic sum rule based on previously calculated d2asr
   !call asrq0_apply(asrq0, natom, mpert, msize, crystal%xcart, d2cart)
   if (inp%asr==1 .or. inp%asr==2 .or. inp%asr==5) then
     call asria_corr(inp%asr,d2asr,d2cart,mpert,natom)
   end if

!  Calculation of the eigenvectors and eigenvalues
!  of the dynamical matrix
   ABI_ALLOCATE(displ,(2*3*natom*3*natom))
   ABI_ALLOCATE(eigval,(3,natom))
   call dfpt_phfrq(amu,displ,d2cart,eigval,eigvec,indsym,&
&   mpert,msym,natom,nsym,ntypat,phfreq(:,iqpt),qphnrm(1),spqpt(:,iqpt),rprimd,inp%symdynmat,&
&   symrel,symafm,typat,ucvol)
   ABI_DEALLOCATE(displ)
   ABI_DEALLOCATE(eigval)
   ABI_DEALLOCATE(d2cart)


!  Read the next bloks to find the next q point.
   found=0 ; iqpt2_previous=iqpt2
   do while (iqpt2<nblok2)
     iqpt2=iqpt2+1
     call read_blok8(ddb_eig2,iqpt2,mband,mpert_eig2,msize2,&
&     nkpt,ddbun,blkval2(:,:,:,:),kpnt(:,:,1))
     !write (300,*) 'blkval2 _bis_ in thmeig'
     !write (300,*) blkval2
     diff_qpt(:)=ddb_eig2%qpt(1:3,iqpt2)/ddb_eig2%nrm(1,iqpt2)-spqpt(:,iqpt)
     if(diff_qpt(1)**2+diff_qpt(2)**2+diff_qpt(3)**2 < DDB_QTOL )then
       found=1
       exit
     end if
   end do

!  Usually, the q points come in the right order. However, this is not always the case...
   if(found==0)then

!    If the EIG2 database file has to be read again, close it, then search for the right q point,
!    from the beginning of the file
     close(ddbun)

     call ddb_hdr_open_read(ddb_hdr, eig2_filnam, ddbun, DDB_VERSION, msym=msym)
     call ddb_hdr_free(ddb_hdr)

!    And examine again the EIG2 file. Still, not beyond the previously examined value.
     found=0
     do iqpt2=1,iqpt2_previous
       call read_blok8(ddb_eig2,iqpt2,mband,mpert_eig2,msize2,&
&       nkpt,ddbun,blkval2(:,:,:,:),kpnt(:,:,1))
       diff_qpt(:)=ddb_eig2%qpt(1:3,iqpt2)/ddb_eig2%nrm(1,iqpt2)-spqpt(:,iqpt)
       if(diff_qpt(1)**2+diff_qpt(2)**2+diff_qpt(3)**2 < DDB_QTOL )then
         found=1
         exit
       end if
     end do

     if(found==0)then
       write(message,'(a,3es16.6,2a)')&
&       'Was unable to find in EIG2 DDB file, the blok for point ',spqpt(:,iqpt),ch10,&
&       'Action: compute the contribution from this point, and merge it in your EIG2 DDB file.'
       MSG_ERROR(message)
     end if

   end if

!  Put blkval2 in cartesian coordinates
   call carttransf(ddb_eig2%flg,blkval2,carflg_eig2,gprimd,iqpt,mband,&
&   mpert_eig2,msize2,natom,nblok2,nkpt,rprimd)

   do imod=1,3*natom

!    Calculate the derivative
     deigr(:,:) = zero
     deigi(:,:) = zero
     dwtermr(:,:)=zero
     dwtermi(:,:)=zero
     index=0
     do iatom1=1,natom
       do idir1=1,3
         do iatom2=1,natom
!          Compute factor for SE term
           if(phfreq(imod,iqpt)<tol6)then
             factr = zero
           else
             factr=one/sqrt(amu(typat(iatom1))*amu(typat(iatom2)))/phfreq(imod,iqpt)/amu_emass
           end if

           do idir2=1,3
             index = idir1 + 3*((iatom1 - 1) + natom * ((idir2-1)+3*(iatom2-1)))

!            Compute products of polarization vectors
             vecr = eigvec(1,idir1,iatom1,imod)*eigvec(1,idir2,iatom2,imod)+&
&             eigvec(2,idir1,iatom1,imod)*eigvec(2,idir2,iatom2,imod)
             veci = eigvec(2,idir1,iatom1,imod)*eigvec(1,idir2,iatom2,imod)-&
&             eigvec(1,idir1,iatom1,imod)*eigvec(2,idir2,iatom2,imod)

             vec1r = eigvec(1,idir1,iatom1,imod)*eigvec(1,idir2,iatom1,imod)+&
&             eigvec(2,idir1,iatom1,imod)*eigvec(2,idir2,iatom1,imod)
             vec1i = eigvec(2,idir1,iatom1,imod)*eigvec(1,idir2,iatom1,imod)-&
&             eigvec(1,idir1,iatom1,imod)*eigvec(2,idir2,iatom1,imod)

             vec2r = eigvec(1,idir1,iatom2,imod)*eigvec(1,idir2,iatom2,imod)+&
&             eigvec(2,idir1,iatom2,imod)*eigvec(2,idir2,iatom2,imod)
             vec2i = eigvec(2,idir1,iatom2,imod)*eigvec(1,idir2,iatom2,imod)-&
&             eigvec(1,idir1,iatom2,imod)*eigvec(2,idir2,iatom2,imod)

!            Compute factor for DW term
             if(phfreq(imod,iqpt)<tol6)then
               fact2r = zero
               fact2i = zero
             else
               fact2r = -wghtq(iqpt)*(vec1r/amu(typat(iatom1)) + vec2r/amu(typat(iatom2)))/phfreq(imod,iqpt)/&
&               amu_emass/2 !/norm(idir1)/norm(idir2)
               fact2i = -wghtq(iqpt)*(vec1i/amu(typat(iatom1)) + vec2i/amu(typat(iatom2)))/phfreq(imod,iqpt)/&
&               amu_emass/2 !/norm(idir1)/norm(idir2)
             end if

             multr(:,:) =(blkval2(1,index,:,:)*vecr - blkval2(2,index,:,:)*veci) !/(norm(idir1)*norm(idir2))
             multi(:,:) =(blkval2(1,index,:,:)*veci + blkval2(2,index,:,:)*vecr) !/(norm(idir1)*norm(idir2))


!            Debye-Waller Term
             if(thmflag==3 .or. thmflag==5 .or. thmflag==7) then
               dwtermr(1:mband,1:nkpt)=dwtermr(1:mband,1:nkpt)+fact2r*blkval2gqpt(1,index,:,:)-fact2i*blkval2gqpt(2,index,:,:)
               dwtermi(1:mband,1:nkpt)=dwtermi(1:mband,1:nkpt)+fact2r*blkval2gqpt(2,index,:,:)+fact2i*blkval2gqpt(1,index,:,:)
             end if

!            Self-energy Term (Fan)
             deigr(1:mband,1:nkpt) = deigr(1:mband,1:nkpt) + wghtq(iqpt)*factr*multr(1:mband,1:nkpt)
             deigi(1:mband,1:nkpt) = deigi(1:mband,1:nkpt) + wghtq(iqpt)*factr*multi(1:mband,1:nkpt)

           end do !idir2
         end do !iatom2
       end do !idir1
     end do !iatom1
!    Eigenvalue derivative or broadening
     if(thmflag==3 .or. thmflag==5 .or. thmflag==7) then
       dednr(1:mband,1:nkpt,imod,iqpt) = deigr(1:mband,1:nkpt) + dwtermr(1:mband,1:nkpt)
       dedni(1:mband,1:nkpt,imod,iqpt) = deigi(1:mband,1:nkpt) + dwtermi(1:mband,1:nkpt)
     else if(thmflag==4 .or. thmflag==6 .or. thmflag==8) then
       dednr(1:mband,1:nkpt,imod,iqpt) = pi*deigr(1:mband,1:nkpt)
       dedni(1:mband,1:nkpt,imod,iqpt) = pi*deigi(1:mband,1:nkpt)
     end if

   end do ! imod
 end do !iqpt

 close(ddbun)


!=============================================================================
!3) Evaluation of the Eliashberg type spectral function
!and phonon DOS via gaussian broadning
!=============================================================================

 if(telphint==1)then
   ng2f = 500  ! number of frequencies
   omega_min=zero
   omega_max=zero
   do iqpt=1,nqpt
     do imod=1,3*natom
       omega_min = min(omega_min,phfreq(imod,iqpt))
       omega_max = max(omega_max,phfreq(imod,iqpt))
     end do
   end do

   ABI_ALLOCATE(dos_phon,(ng2f))
   ABI_ALLOCATE(g2f,(mband,nkpt,ng2f))
   ABI_ALLOCATE(tmpg2f,(mband,nkpt,ng2f))
   ABI_ALLOCATE(tmpphondos,(ng2f))

   write(std_out,'(a,es13.6)') 'omega_min :', omega_min
   write(std_out,'(a,es13.6)') 'omega_max :', omega_max
   write(std_out,'(a,i8)') 'ng2f :', ng2f

   omega_max = omega_max + 0.1 * omega_max
   domega = (omega_max-omega_min)/(ng2f-one)

   gaussprefactor = sqrt(piinv) / g2fsmear
   gaussfactor = one / g2fsmear

   g2f(:,:,:) = zero
   dos_phon(:) = zero

   do iqpt=1,nqpt
     do imod=1,3*natom
       omega = omega_min
       tmpg2f(:,:,:) = zero
       tmpphondos(:) = zero
       do iomega=1,ng2f
         xx = (omega-phfreq(imod,iqpt))*gaussfactor
         gaussval = gaussprefactor*exp(-xx*xx)
         tmpg2f(:,:,iomega) = tmpg2f(:,:,iomega) + gaussval*dednr(:,:,imod,iqpt)
         tmpphondos(iomega) = tmpphondos(iomega) + gaussval
         omega = omega+domega
       end do

       g2f(:,:,:) = g2f(:,:,:) + tmpg2f(:,:,:)
       dos_phon(:) = dos_phon(:) + tmpphondos(:)

     end do !imod
   end do !iqpt

   dos_phon(:) = dos_phon(:) / nqpt

!  output the g2f
   unit_g2f = 108
   call outg2f(domega,omega_min,omega_max,elph_base_name,g2f,g2fsmear,kpnt,mband,ng2f,nkpt,nqpt,1,telphint,unit_g2f)

!  output the phonon DOS
   unit_phdos = 108
   call outphdos(domega,dos_phon,omega_min,omega_max,elph_base_name,g2fsmear,ng2f,nqpt,1,telphint,unit_g2f)


   ABI_DEALLOCATE(dos_phon)
   ABI_DEALLOCATE(g2f)
   ABI_DEALLOCATE(tmpg2f)
   ABI_DEALLOCATE(tmpphondos)

 end if !telphint

!=======================================================================
!4) Evaluation of the Eliashberg type spectral function
!and phonon DOS via improved tetrahedron method
!=======================================================================

 if(telphint==0)then

!  make dimension-ful rprimd and gprimd for transformation of derivatives to cartesian coordinates.
   call mkrdim(acell,rprim,rprimd)
   call matr3inv(rprimd,gprimd)

!  Q point Grid
   qpt_full(:,:) = ddb%qpt(1:3,:)

!  Trivial Q point index
   do iqpt=1,nqpt
     indqpt(iqpt)=iqpt
     qptnrm(iqpt)= qpt_full(1,iqpt)*qpt_full(1,iqpt)+qpt_full(2,iqpt)*qpt_full(2,iqpt)+qpt_full(3,iqpt)*qpt_full(3,iqpt)
   end do

!  Build qlatt from scratch (for 5.7)
   tol = 0.1_dp
   ilatt = 0
   call sort_dp(nqpt,qptnrm,indqpt,tol)

   do iqpt1=1,nqpt-2
     mesh(1:3,1) = qpt_full(1:3,indqpt(iqpt1))
     do iqpt2=iqpt1+1,nqpt-1
       mesh(1:3,2)= qpt_full(1:3,indqpt(iqpt2))
       do iqpt3=iqpt2+1,nqpt
         mesh(1:3,3)= qpt_full(1:3,indqpt(iqpt3))
         det = mesh(1,1)*mesh(2,2)*mesh(3,3) + mesh(1,2)*mesh(2,3)*mesh(3,1) + mesh(1,3)*mesh(2,1)*mesh(3,2) &
&         -mesh(3,1)*mesh(2,2)*mesh(1,3) - mesh(3,2)*mesh(2,3)*mesh(1,1) - mesh(3,3)*mesh(2,1)*mesh(1,2)
         invdet = one/det
         if (abs(nint(invdet))==nqpt .and. abs(invdet)-nqpt < tol) then
           ilatt = 1
           qlatt(:,:) = mesh(:,:)
           exit
         end if
       end do
       if(ilatt==1) exit
     end do
     if(ilatt==1) exit
   end do

!  error message if qlatt not found and stop
   if(ilatt==0) then
     write(message, '(a,a)' ) &
&     ' Could not find homogeneous basis vectors for Q point grid ',ch10
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
     MSG_ERROR("Aborting now")
   end if

!  test if qlatt is righthanded and possibly fixe it
   if(invdet < 0) then
     tempqlatt(:) = qlatt(:,2)
     qlatt(:,2) = qlatt(:,1)
     qlatt(:,1) = tempqlatt(:)
   end if

   write(std_out,*) 'qlatt',qlatt

!  test if qlatt generates all Q points  TO DO

!  Get tetrahedra, ie indexes of the full kpoints at their summits
   call init_tetra(indqpt,gprimd,qlatt,qpt_full,nqpt, tetrahedra, ierr, errstr, xmpi_comm_self)
   ABI_CHECK(ierr==0,errstr)

   rcvol = abs (gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
&   -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
&   +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))

!  Calculate weights for phonon DOS
!  Special precautions must be taking for Gamma point
!  because of non-analytic term.
!  Non-analyticity must be taken out and treated separatly.

   nene = 100     !nene=number of energies for DOS
   enemin = minval(phfreq)
   enemax = maxval(phfreq)
   deltaene = (enemax-enemin)/dble(nene-1)
!  redefine enemin enemax to be at rounded multiples of deltaene
!  enemin = elph_ds%fermie - dble(ifermi)*deltaene
!  enemax = elph_ds%fermie + dble(nene-ifermi-1)*deltaene

   ABI_ALLOCATE(tweight,(nqpt,nene))
   ABI_ALLOCATE(dtweightde,(nqpt,nene))
   ABI_ALLOCATE(intweight,(3*natom,nqpt,nene))
   ABI_ALLOCATE(indtweightde,(3*natom,nqpt,nene))

   do iband=1,3*natom
     eigen_in(:) = phfreq(iband,:)

!    calculate general integration weights at each irred kpoint
!    as in Blochl et al PRB 49 16223 [[cite:Bloechl1994a]]
     call get_tetra_weight(eigen_in,enemin,enemax,&
&     one,nene,nqpt,tetrahedra,bcorr0,&
&     tweight,dtweightde,xmpi_comm_self)

     intweight(iband,:,:) = tweight(:,:)
     indtweightde(iband,:,:) = dtweightde(:,:)

   end do !iband

!  intdtweightse(nband,nqpt,nene) represents the weight in each energy bin for every kpt and every band
!  So phonon DOS is calculated (neglecting the non-analyticity contribution for now !!!)

   ABI_ALLOCATE(total_dos,(nene))
   ABI_ALLOCATE(g2f,(mband,nkpt,nene))

   total_dos(:) = zero
   do iband=1,3*natom
     do iqpt=1,nqpt
       total_dos(:) = total_dos + indtweightde(iband,iqpt,:)
     end do
   end do

!  For the g2f function
!  Right now for one electronic band and one K point: dednr(1:mband,1:nkpt,imod,iqpt)
!  Once again must pay close attention to the Gamma point
   g2f(:,:,:) = zero
   do ii=1,mband
     do ikpt=1,nkpt
       do iband=1,3*natom
         do iqpt=1,nqpt
           g2f(ii,ikpt,:) = g2f(ii,ikpt,:) + dednr(ii,ikpt,iband,iqpt) * indtweightde(iband,iqpt,:)
         end do
       end do
     end do
   end do

!  output the g2f
   unit_g2f = 108
   call outg2f(deltaene,enemin,enemax,elph_base_name,g2f,g2fsmear,kpnt,mband,nene,nkpt,nqpt,tetrahedra%ntetra,telphint,unit_g2f)

!  output the phonon DOS
   unit_phdos = 108
   call outphdos(deltaene,total_dos,enemin,enemax,elph_base_name,g2fsmear,nene,nqpt,tetrahedra%ntetra,telphint,unit_g2f)

   ABI_DEALLOCATE(tweight)
   ABI_DEALLOCATE(dtweightde)
   ABI_DEALLOCATE(intweight)
   ABI_DEALLOCATE(indtweightde)
   ABI_DEALLOCATE(total_dos)
   ABI_DEALLOCATE(g2f)
 end if !telphint

!=======================================================================
!5) direct evaluation of thermal corrections
!=======================================================================

!open TBS file
 outfile = trim(elph_base_name)//"_TBS"
 if (open_file(outfile,message,newunit=unitout,form='formatted',status='unknown') /= 0) then
   MSG_ERROR(message)
 end if
 write(unitout,'(a)')'thmeig: Thermal Eigenvalue corrections (eV)'

 slope(:,:,:) = zero
 zeropoint(:,:,:) = zero
!Loop on temperatures
 do itemper= 1, ntemper
   tmp=tempermin+temperinc*float(itemper-1)
   thmeigen(:,:,:) = zero

!  Sum on all phonon wavevectors and modes
   do iqpt=1,nqpt
     do imod=1,3*natom

!      Bose-Einstein distribution
! jmb overflow with exp(). So, select bosein to be still significant wrt half
       if(phfreq(imod,iqpt)<tol6 .or. (phfreq(imod,iqpt)/(kb_HaK*tmp)) > -log(tol16))then
         bosein = zero
       else
         bosein = one/(exp(phfreq(imod,iqpt)/(kb_HaK*tmp))-one)
       end if

!      Calculate total
       thmeigen(1,1:mband,1:nkpt) = thmeigen(1,1:mband,1:nkpt) + dednr(1:mband,1:nkpt,imod,iqpt)*(bosein+half)
       thmeigen(2,1:mband,1:nkpt) = thmeigen(2,1:mband,1:nkpt) + dedni(1:mband,1:nkpt,imod,iqpt)*(bosein+half)

       if(itemper==1)then
!        Calculate slope of linear regime
         if(phfreq(imod,iqpt)<tol6)then
           slope(1,1:mband,1:nkpt) = slope(1,1:mband,1:nkpt)
           slope(2,1:mband,1:nkpt) = slope(2,1:mband,1:nkpt)
         else
           slope(1,1:mband,1:nkpt) = slope(1,1:mband,1:nkpt) + dednr(1:mband,1:nkpt,imod,iqpt)*(kb_HaK/phfreq(imod,iqpt))
           slope(2,1:mband,1:nkpt) = slope(2,1:mband,1:nkpt) + dedni(1:mband,1:nkpt,imod,iqpt)*(kb_HaK/phfreq(imod,iqpt))
         end if
!        Calculate zero-point renormalization
         zeropoint(1,1:mband,1:nkpt) = zeropoint(1,1:mband,1:nkpt) + dednr(1:mband,1:nkpt,imod,iqpt)*half
         zeropoint(2,1:mband,1:nkpt) = zeropoint(2,1:mband,1:nkpt) + dedni(1:mband,1:nkpt,imod,iqpt)*half

       end if
     end do ! imod
   end do !iqpt

!  Write temperature independent results
   if(itemper==1)then
     write(unitout,'(a)')'Temperature independent results (zero-point renormalization and slope)'
     do ikpt=1,nkpt
       write(unitout,'(a,3es16.8)')' Kpt :', kpnt(:,ikpt,1)
       do iband=1,mband
         write(unitout,'(4d22.14)') Ha_eV*zeropoint(1,iband,ikpt),Ha_eV*zeropoint(2,iband,ikpt),&
&         Ha_eV*slope(1,iband,ikpt),Ha_eV*slope(2,iband,ikpt)
       end do
     end do
     write(unitout,'(a)')'Temperature dependent corrections'
   end if
!  Write result in file for each temperature
   write(unitout,'(a,es10.3,a)')'T :', tmp,' K'
   do ikpt=1,nkpt
     write(unitout,'(a,3es16.8)')' Kpt :', kpnt(:,ikpt,1)
     do iband=1,mband
       write(unitout,'(2d22.14)') Ha_eV*thmeigen(1,iband,ikpt), Ha_eV*thmeigen(2,iband,ikpt)
     end do
   end do
 end do !itemper

 close(unitout)

!Write temperature-independent results to the main output file
 write(iout,'(a)')' '
 write(iout,'(80a)') ('-',ii=1,80)
 write(iout,'(a)')' '
 write(iout,'(a)')' Electron-phonon change of electronic structure.'
 write(iout,'(a)')' The temperature-dependent values are written in the _TBS file.'
 write(iout,'(a)')' Here follows, for each electronic wavevector and band :'
 write(iout,'(a)')'      zero-point renormalisation (Ha) and linear slope (Ha/Kelvin)'
 do ikpt=1,nkpt
   write(iout,'(2a,i6,a,3es16.6)')ch10,' Kpt number ',ikpt,', with reduced coordinates :',kpnt(:,ikpt,1)
   do iband=1,mband
     write(iout,'(i6,2es20.6)') iband,zeropoint(1,iband,ikpt),slope(1,iband,ikpt)
   end do
 end do

 ABI_DEALLOCATE(typat)
 ABI_DEALLOCATE(atifc)
 ABI_DEALLOCATE(zion)
 ABI_DEALLOCATE(amu)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(symafm)
 ABI_DEALLOCATE(spinat)
 ABI_DEALLOCATE(symrel)
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(indsym)
 ABI_DEALLOCATE(tnons)
 ABI_DEALLOCATE(deigi)
 ABI_DEALLOCATE(deigr)
 ABI_DEALLOCATE(dwtermi)
 ABI_DEALLOCATE(dwtermr)
 ABI_DEALLOCATE(multi)
 ABI_DEALLOCATE(multr)
 ABI_DEALLOCATE(slope)
 ABI_DEALLOCATE(thmeigen)
 ABI_DEALLOCATE(zeropoint)

 ABI_DEALLOCATE(dedni)
 ABI_DEALLOCATE(dednr)
 if(allocated(indqpt)) then
   ABI_DEALLOCATE(indqpt)
 end if
 ABI_DEALLOCATE(eigen_in)
 ABI_DEALLOCATE(qpt_full)
 ABI_DEALLOCATE(qptnrm)
 ABI_DEALLOCATE(wghtq)
 ABI_DEALLOCATE(spqpt)
 ABI_DEALLOCATE(eigvec)
 ABI_DEALLOCATE(phfreq)

 ABI_DEALLOCATE(blkval2)
 ABI_DEALLOCATE(blkval2gqpt)
 ABI_DEALLOCATE(kpnt)
 ABI_DEALLOCATE(carflg_eig2)



 call ddb_free(ddb_eig2)

 call destroy_tetra(tetrahedra)

end subroutine thmeig
!!***

!!****f* m_thmeig/outphdos
!! NAME
!! outphdos
!!
!! FUNCTION
!!  Print out phonon density of states
!!
!! INPUTS
!!   deltaene = step on energy/frequency grid, in Hartree
!!   dos_phon = phonon DOS calculated on a grid
!!   enemin = minimal frequency
!!   enemax = maximal frequency
!!   filnam = file name for output to disk
!!   g2fsmear = smearing width
!!   nene = number of points on energy axis
!!   nqpt = number of q-points
!!   ntetra = number of tetrahedra, if tetrahedron interpolation is used
!!   telphint = flag for el-phonon interpolation method (to indicate Gaussian or tetrahedron integration)
!!   unit_phdos = unit for phonon DOS output
!!
!!
!! OUTPUT
!!  only write
!!
!! SIDE EFFECTS
!!
!! NOTES
!!   FIXME
!!   overcomplete inputs. Eliminate unit_phdos (just filnam) and deltaene (gotten from max-min/nene)
!!
!! PARENTS
!!      thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine outphdos(deltaene,dos_phon,enemin,enemax,filnam,g2fsmear,nene,nqpt,ntetra,telphint,unit_phdos)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nene,nqpt,ntetra,telphint,unit_phdos
 character(len=fnlen),intent(in) :: filnam
 real(dp) :: deltaene,enemin,enemax,g2fsmear
!arrays
 real(dp) :: dos_phon(nene)

!Local variables-------------------------------
!scalars
 integer :: iomega,iost,step10
 real(dp) :: dos_effective,omega
 character(len=fnlen) :: outfile
 character(len=500) :: message
!arrays

! *************************************************************************

   outfile = trim(filnam) // '_PDS'
   write(message, '(3a)')ch10,&
&   ' Will write phonon DOS in file ',trim(outfile)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   write(message, '(4a)')ch10,&
&   ' For checking purposes, write ten values in the present file.',ch10,&
&   '       Index    Energy (in Ha)      DOS '
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   open (unit=unit_phdos,file=outfile,status='replace',iostat=iost)
   if (iost /= 0) then
     write (message,'(3a)')' Opening file ',trim(outfile),' as new'
     MSG_ERROR(message)
   end if

   write (unit_phdos,'(a)') '#'
   write (unit_phdos,'(a)') '# ABINIT package : phonon DOS file'
   write (unit_phdos,'(a)') '#'
   write (unit_phdos,'(a,i10)') '#   Number of Qpoints integrated over : ', nqpt
   write (unit_phdos,'(a,i10)') '#   Number of energy points : ', nene
   write (unit_phdos,'(a,es16.6,a,es16.6,a)') '#   between omega_min = ', enemin, &
&   ' Ha and omega_max = ', enemax, ' Ha'
   if(telphint==1)then
     write (unit_phdos,'(a,es16.6)') '#   The smearing width for gaussians is ', g2fsmear
   end if
   if(telphint==0)then
     write (unit_phdos,'(a,i10)') '#   Number of tetrahedrons', ntetra
   end if
   write (unit_phdos,'(a)') '#'
   write (unit_phdos,'(a)') '#      Index    Energy (in Ha)      DOS '

   omega = enemin
   do iomega=1,nene
     dos_effective=dos_phon(iomega)
     if(abs(dos_effective)<tol16)then
       dos_effective=zero
     end if
     step10=nene/10
     if(mod(iomega,step10)==1)write (std_out,'(i10,es18.6,es18.6)')iomega, omega, dos_effective
     if(mod(iomega,step10)==1)write (ab_out,'(i10,es18.6,es18.6)')iomega, omega, dos_effective
     write (unit_phdos, '(i10,es18.6,es18.6)')iomega, omega, dos_effective
     omega=omega+deltaene
   end do

   close (unit=unit_phdos)

 end subroutine outphdos
!!***

!!****f* m_thmeig/outg2f
!! NAME
!! outg2f
!!
!! FUNCTION
!!  Output g2f function to file. FIXME: Paul, please explain what g2f is.
!!  Probably a variant on the Eliashberg spectral function a2F
!!
!! INPUTS
!!
!! OUTPUT
!!  only write
!!
!! PARENTS
!!      thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine outg2f(deltaene,enemin,enemax,filnam,g2f,g2fsmear,kpnt,mband,nene,nkpt,nqpt,ntetra,telphint,unit_g2f)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nene,nkpt,nqpt,ntetra,telphint,unit_g2f
 character(len=fnlen),intent(in) :: filnam
 real(dp) :: deltaene,enemin,enemax,g2fsmear
!arrays
 real(dp) :: g2f(mband,nkpt,nene),kpnt(3,nkpt,nqpt)

!Local variables-------------------------------
!scalars
 integer :: iband,ikpt,iomega,iost
 real(dp) :: omega
 character(len=fnlen) :: outfile
 character(len=500) :: message
!arrays

! *************************************************************************

!output the g2f
   outfile = trim(filnam) // '_G2F'
   open (unit=unit_g2f,file=outfile,status='unknown',iostat=iost)
   if (iost /= 0) then
     write (message,'(3a)')' thmeig : ERROR- opening file ',trim(outfile),' as new'
     MSG_ERROR(message)
   end if

   write(std_out,*) ' g2f function'
   write (unit_g2f,'(a)') '#'
   write (unit_g2f,'(a)') '# ABINIT package : g2f file'
   write (unit_g2f,'(a)') '#'
   write (unit_g2f,'(a,I10)') '#     number of qpoints integrated over : ', nqpt
   write (unit_g2f,'(a,I10)') '#     number of energy points : ', nene
   write (unit_g2f,'(a,E16.6,a,E16.6,a)') '#       between omega_min = ', enemin, &
&   ' Ha and omega_max = ', enemax, ' Ha'
   if(telphint==1)then
     write (unit_g2f,'(a,E16.6)') '#   and the smearing width for gaussians is ', g2fsmear
     write (unit_g2f,'(a)') '#'
   end if
   if(telphint==0)then
     write (unit_g2f,'(a,I10)') '#   number of tetrahedrons', ntetra
     write (unit_g2f,'(a)') '#'
   end if

!Write only the a2f function for the first K point
!ikpt=1
   do ikpt=1,nkpt
     write(unit_g2f,'(a,3es16.8)')' Kpt :', kpnt(:,ikpt,1)
     do iband=1,mband
       write(unit_g2f,*) 'band :', iband
       omega = enemin
       do iomega=1,nene
         write (unit_g2f,*) omega*Ha_eV*1000, g2f(iband, ikpt,iomega)
         omega=omega+deltaene
       end do
     end do
   end do

   close (unit=unit_g2f)

end subroutine outg2f
!!***

end module m_thmeig
!!***
