!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_harmonic_thermo
!! NAME
!! m_harmonic_thermo
!!
!! FUNCTION
!! This routine to calculate phonon density of states,
!! thermodynamical properties, Debye-Waller factor, and atomic mean square velocity
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (CL, XG)
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

module m_harmonic_thermo

 use defs_basis
 use m_errors
 use m_abicore
 use m_sortph
 use m_xmpi

 use m_io_tools,       only : open_file
 use m_symtk,          only : matr3inv
 use m_dynmat,         only : gtdyn9
 use m_geometry,       only : mkrdim
 use m_crystal,        only : crystal_t
 use m_anaddb_dataset, only : anaddb_dataset_type
 use m_ifc,            only : ifc_type
 use m_kpts,           only : smpbz
 use m_symkpt,         only : symkpt

 implicit none

 private
!!***

 public :: harmonic_thermo
!!***

contains
!!***

!!****f* m_harmonic_thermo/harmonic_thermo
!!
!! NAME
!! harmonic_thermo
!!
!! FUNCTION
!! This routine to calculate phonon density of states,
!! thermodynamical properties, Debye-Waller factor, and atomic mean square velocity
!!
!! INPUTS
!! Crystal<crystal_t>=data type gathering info on the crystalline structure.
!! Ifc<ifc_type>=Object containing the interatomic force constants.
!!    %atmfrc(2,3,natom,3,natom,nrpt) = Interatomic Forces in real space
!!    %dyewq0(3,3,natom)=Ewald part of the dynamical matrix, at q=0.
!!    %rpt(3,nrpt)=canonical coordinates of the R points in the unit cell These coordinates are normalized (=> * acell(3)!!)
!!    %nrpt=number of R points in the Big Box
!!    %trans(3,natom)=atomic translations : xred = rcan + trans
!!    %wghatm(natom,natom,nrpt)=weights associated to a pair of atoms and to a R vector
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! anaddb_dtset= (derived datatype) contains all the input variables
!! iout =unit number for output
!! natom=number of atoms in the unit cell
!! outfilename_radix=radix of anaddb output file name: append _THERMO for thermodynamic quantities
!! comm=MPI communicator
!!
!! OUTPUT
!!
!! NOTES
!! dosinc=increment between the channels for the phonon DOS in cm-1
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      end_sortph,ifc_fourq,matr3inv,mkrdim,smpbz,sortph,symkpt,wrtout
!!
!! SOURCE

subroutine harmonic_thermo(Ifc,Crystal,amu,anaddb_dtset,iout,outfilename_radix,comm,thmflag)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iout,comm
 integer,intent(in),optional :: thmflag
 character(len=*),intent(in) :: outfilename_radix
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset
 type(crystal_t),intent(in) :: Crystal
 type(ifc_type),intent(in) :: Ifc
!arrays
 real(dp),intent(in) :: amu(Crystal%ntypat)

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 integer :: convth,facbrv,iatom,ichan,icomp,igrid,iii,iiii,ij,ind
 integer :: iqpt2,isym,itemper,iwchan,jjj,mqpt2,nchan,ngrids,natom
 integer :: nqpt2,nspqpt,ntemper,nwchan,option,timrev
 integer :: thermal_unit
 integer :: bij_unit
 integer :: vij_unit
 integer :: nomega, iomega
 real(dp) :: change,cothx,diffbb,dosinc,expm2x,factor,factorw,factorv,gerr
 real(dp) :: ggsum,ggsumsum,ggrestsum
 real(dp) :: gijerr,gijsum,gnorm,ln2shx,qphnrm,relchg,tmp,wovert,thmtol
 logical :: part1,part2
 character(len=500) :: msg
 character(len=fnlen) :: thermal_filename
 character(len=fnlen) :: bij_filename
 character(len=fnlen) :: vij_filename
!arrays
 integer :: symrec(3,3,Crystal%nsym),symrel(3,3,Crystal%nsym)
 real(dp) :: symrec_cart(3,3,Crystal%nsym)
 integer :: igqpt2(3),ii(6),jj(6),qptrlatt(3,3)
 integer,allocatable :: indqpt1(:),nchan2(:),bz2ibz_smap(:,:)
 real(dp) :: gprimd(3,3),qphon(3),rprimd(3,3),tens1(3,3),tens2(3,3)
 real(dp) :: displ(2*3*Crystal%natom*3*Crystal%natom)
 real(dp) :: eigvec(2,3,Crystal%natom,3*Crystal%natom)
 real(dp) :: phfrq(3*Crystal%natom)
 real(dp),allocatable :: bbij(:,:,:),bij(:,:,:),energy(:),energy0(:),entropy(:)
 real(dp),allocatable :: entropy0(:),free(:),free0(:),gdos(:,:),gg(:,:),gg_sum(:,:),gg_rest(:,:)
 real(dp),allocatable :: ggij(:,:,:,:),gij(:,:,:,:),gw(:,:),qpt2(:,:),spheat(:)
 real(dp),allocatable :: spheat0(:),spqpt2(:,:),wme(:),wtq(:),wtq2(:)
 real(dp),allocatable :: wtq_folded(:),vij(:,:,:)
 real(dp),allocatable :: phon_dos(:)
 logical,allocatable :: wgcnv(:),wgijcnv(:)

! *********************************************************************

 ! For the time being, this routine can use only 1 MPI process
 if (xmpi_comm_rank(comm) /= master) return

 natom = Crystal%natom
 symrel = Crystal%symrel
 symrec = Crystal%symrec

 call mkrdim(Ifc%acell,Ifc%rprim,rprimd)
 call matr3inv(rprimd,gprimd)
 do isym = 1, Crystal%nsym
   symrec_cart(:,:,isym) = matmul( gprimd, matmul(dble(symrec(:,:,isym)), transpose(rprimd)) )
! net result is tens1 = rprimd symrec^T gprimd^T   tens1   gprimd symrec rprimd^T
 end do

 thermal_filename=trim(outfilename_radix)//"_THERMO"
 if (open_file(thermal_filename, msg, newunit=thermal_unit) /= 0) then
   MSG_ERROR(msg)
 end if

 write(thermal_unit,*) '#'
 write(thermal_unit,*) '#  Thermodynamical quantities calculated by ANADDB'
 write(thermal_unit,*) '#'

 bij_filename=trim(outfilename_radix)//"_DEBYE_WALLER"
 if (open_file(bij_filename, msg, newunit=bij_unit) /= 0) then
   MSG_ERROR(msg)
 end if

 vij_filename=trim(outfilename_radix)//"_VELOC_SQUARED"
 if (open_file(vij_filename, msg, newunit=vij_unit) /= 0) then
   MSG_ERROR(msg)
 end if

 thmtol = anaddb_dtset%thmtol
 nchan=anaddb_dtset%nchan
 ntemper=anaddb_dtset%ntemper
 nwchan=anaddb_dtset%nwchan
 ngrids=anaddb_dtset%ngrids

 ABI_ALLOCATE(bbij,(6,natom,ntemper))
 ABI_ALLOCATE(bij,(6,natom,ntemper))
 ABI_ALLOCATE(vij,(6,natom,ntemper))
 ABI_ALLOCATE(energy,(ntemper))
 ABI_ALLOCATE(energy0,(ntemper))
 ABI_ALLOCATE(entropy,(ntemper))
 ABI_ALLOCATE(entropy0,(ntemper))
 ABI_ALLOCATE(free,(ntemper))
 ABI_ALLOCATE(free0,(ntemper))
!Doubling the size (nchan) of gg_sum is needed, because the maximum sum frequency is the double
!the for maximum frequency. However, for the write statement, it is better to double
!also the size of the other arrays ...
 ABI_ALLOCATE(gdos,(2*nchan,nwchan))
 ABI_ALLOCATE(gg,(2*nchan,nwchan))
 ABI_ALLOCATE(gg_sum,(2*nchan,nwchan))
 ABI_ALLOCATE(gg_rest,(2*nchan,nwchan))
 ABI_ALLOCATE(ggij,(6,natom,nchan,nwchan))
 ABI_ALLOCATE(gij,(6,natom,nchan,nwchan))
 ABI_ALLOCATE(nchan2,(nwchan))
 ABI_ALLOCATE(spheat,(ntemper))
 ABI_ALLOCATE(spheat0,(ntemper))
 ABI_ALLOCATE(wgcnv,(nwchan))
 ABI_ALLOCATE(wgijcnv,(nwchan))
 ABI_ALLOCATE(wme,(ntemper))
 ABI_ALLOCATE(gw,(nchan,nwchan))


!initialize ii and jj arrays
 ii(1)=1 ; ii(2)=2 ; ii(3)=3
 ii(4)=1 ; ii(5)=1 ; ii(6)=2
 jj(1)=1 ; jj(2)=2 ; jj(3)=3
 jj(4)=2 ; jj(5)=3 ; jj(6)=3

!Put zeros in the bins of channels of frequencies
 gdos(:,:)=0._dp
 gg_sum(:,:)=0._dp
 gg_rest(:,:)=0._dp
 gij(:,:,:,:)=0._dp

!Neither part1 nor part2 are assumed converged initially
!None of the channel widths are assumed converged initially

 part1=.false.
 part2=.false.

 wgcnv(:)=.false.
 wgijcnv(:)=.false.

!Thermodynamic quantities are put to zero.
!(If exactly zero, initial convergence tests will fail.)

 free0(:)=1.d-05
 energy0(:)=1.d-05
 entropy0(:)=1.d-05
 spheat0(:)=1.d-05
 bij(:,:,:)=1.d-05
 vij(:,:,:)=1.d-05

!Number of actual channels set
 do iwchan=1,nwchan
   nchan2(iwchan)=(nchan-1)/iwchan+1
 end do

!For different types of Bravais lattices
 facbrv=1
 if(anaddb_dtset%brav==2)facbrv=2
 if(anaddb_dtset%brav==3)facbrv=4

!Loops on the q point grids
 do igrid=1,ngrids

   igqpt2(:)=max((anaddb_dtset%ng2qpt(:)*igrid)/ngrids, 1)
   mqpt2=(igqpt2(1)*igqpt2(2)*igqpt2(3))/facbrv
   ABI_ALLOCATE(qpt2,(3,mqpt2))
   ABI_ALLOCATE(spqpt2,(3,mqpt2))

   option=1
   qptrlatt(:,:)=0
   qptrlatt(1,1)=igqpt2(1)
   qptrlatt(2,2)=igqpt2(2)
   qptrlatt(3,3)=igqpt2(3)
   call smpbz(anaddb_dtset%brav,iout,qptrlatt,mqpt2,nspqpt,1,option,anaddb_dtset%q2shft,spqpt2)

   ABI_ALLOCATE(indqpt1,(nspqpt))
   ABI_ALLOCATE(wtq,(nspqpt))
   ABI_ALLOCATE(wtq_folded,(nspqpt))
   ABI_ALLOCATE(bz2ibz_smap,(6, nspqpt))

!  Reduce the number of such points by symmetrization
   wtq(:)=1.0_dp

   timrev=1
   call symkpt(0,Crystal%gmet,indqpt1,ab_out,spqpt2,nspqpt,nqpt2,Crystal%nsym,symrec,timrev,wtq,wtq_folded, &
     bz2ibz_smap, xmpi_comm_self)

   ABI_DEALLOCATE(bz2ibz_smap)

   ABI_ALLOCATE(wtq2,(nqpt2))
   do iqpt2=1,nqpt2
     wtq2(iqpt2)=wtq_folded(indqpt1(iqpt2))
     qpt2(:,iqpt2)=spqpt2(:,indqpt1(iqpt2))
     !write(std_out,*)' harmonic_thermo : iqpt2, wtq2 :',iqpt2,wtq2(iqpt2)
   end do
   ABI_DEALLOCATE(wtq_folded)

!  Temporary counters are put zero.

   gg(:,:)=zero
   ggij(:,:,:,:)=zero
   gw(:,:)=zero

!  Sum over the sampled q points
   do iqpt2=1,nqpt2

     qphon(:)=qpt2(:,iqpt2)
     qphnrm=1.0_dp

     ! Fourier Interpolation
     call ifc%fourq(Crystal,qphon,phfrq,displ,out_eigvec=eigvec)

     if (present(thmflag)) then
       if (thmflag ==2)then
         call sortph(eigvec,displ,outfilename_radix,natom,phfrq)
       end if
     end if

!    Sum over the phonon modes
     do iii=1,3*natom

!      Slightly negative frequencies are put to zero
!      Imaginary frequencies are also put to zero
       if(phfrq(iii)<0._dp) phfrq(iii)=0._dp

!      Note: frequencies are now in cm-1

!      Sort frequencies into channels of frequencies for each channel width of frequency
       do iwchan=nwchan,1,-1
         if (.not.wgcnv(iwchan))then
           dosinc=dble(iwchan)
           ichan=int(phfrq(iii)*Ha_cmm1/dosinc)+1

           if(ichan>nchan2(iwchan)) then
             write(msg, '(a,es16.6,a,a,a,i7,a,a,a)' )&
&             'There is a phonon frequency,',phfrq(iii),' larger than the maximum one,',ch10,&
&             'as defined by the number of channels of width 1 cm-1, nchan=',nchan,'.',ch10,&
&             'Action: increase nchan (suggestion : double it).'
             MSG_ERROR(msg)
           end if

           gg(ichan,iwchan)=gg(ichan,iwchan)+wtq2(iqpt2)

           gw(ichan,iwchan)=gw(ichan,iwchan)+wtq2(iqpt2)*phfrq(iii)*Ha_cmm1

!          to calculate two phonon DOS for qshift = 0.0
           do iiii=1,3*natom
             if(phfrq(iiii)<0.d0) phfrq(iiii)=0.d0

             ichan=int(abs(phfrq(iii)+phfrq(iiii))*Ha_cmm1/dosinc)+1
             gg_sum(ichan,iwchan)=gg_sum(ichan,iwchan)+wtq2(iqpt2)

             ichan=int(abs(phfrq(iii)-phfrq(iiii))*Ha_cmm1/dosinc)+1
             gg_rest(ichan,iwchan)=gg_rest(ichan,iwchan)+wtq2(iqpt2)
           end do ! end loop for iiii

         end if
       end do

       do iwchan=nwchan,1,-1
         if (.not.wgijcnv(iwchan))then

           dosinc=dble(iwchan)
           ichan=int(phfrq(iii)*Ha_cmm1/dosinc)+1

           if(ichan>nchan2(iwchan)) then
             write(msg, '(a,a,a,a,a)' )&
&             'Phonon frequencies larger than the maximum one,',ch10,&
&             'as defined by the number of channels of width 1 cm-1.',ch10,&
&             'Action: increase nchan (suggestion : double it).'
             MSG_ERROR(msg)
           end if

           do iatom=1,natom
             do ij=1,6
               ggij(ij,iatom,ichan,iwchan)=ggij(ij,iatom,ichan,iwchan)&
&               +wtq2(iqpt2)*&
&               (eigvec(1,ii(ij),iatom,iii)*eigvec(1,jj(ij),iatom,iii)+&
&               eigvec(2,ii(ij),iatom,iii)*eigvec(2,jj(ij),iatom,iii) )
             end do
           end do

         end if
       end do

     end do ! End of the sum over the phonon modes
   end do ! End of the sum over q-points in the irreducible zone

!  deallocate sortph array
   call end_sortph()

!  Symmetrize the gij
   do ichan=1,nchan
     do iwchan=nwchan,1,-1
       do iatom=1,natom
         do ij=1,6
!          Uses bbij as work space
           bbij(ij,iatom,1)=ggij(ij,iatom,ichan,iwchan)
           ggij(ij,iatom,ichan,iwchan)=0.0_dp
         end do
       end do
       do iatom=1,natom
         do isym=1,Crystal%nsym
!          Find the atom that will be applied on atom iatom
           ind=Crystal%indsym(4,isym,iatom)
           do ij=1,6
             tens1(ii(ij),jj(ij))=bbij(ij,ind,1)
           end do
!          complete the 3x3 tensor from the upper triangle
           tens1(2,1)=tens1(1,2)
           tens1(3,1)=tens1(1,3)
           tens1(3,2)=tens1(2,3)
!          Here acomplishes the tensorial operations
!!          make this a BLAS call, or better yet batch the whole thing?
!          2) Apply the symmetry operation on both indices   USING symrec in
!          cartesian coordinates
           do iii=1,3
             do jjj=1,3
               tens2(iii,jjj)=tens1(iii,1)*symrec_cart(1,jjj,isym)&
&               +tens1(iii,2)*symrec_cart(2,jjj,isym)&
&               +tens1(iii,3)*symrec_cart(3,jjj,isym)
             end do
           end do
           do jjj=1,3
             do iii=1,3
               tens1(iii,jjj)=tens2(1,jjj)*symrec_cart(1,iii,isym)&
&               +tens2(2,jjj)*symrec_cart(2,iii,isym)&
&               +tens2(3,jjj)*symrec_cart(3,iii,isym)
             end do
           end do
! net result is tens1 = rprimd symrec^T gprimd^T   tens1   gprimd symrec rprimd^T

!          This accumulates over atoms, to account for all symmetric ones
           do ij=1,6
             ggij(ij,iatom,ichan,iwchan)=ggij(ij,iatom,ichan,iwchan) + tens1(ii(ij),jj(ij))
           end do

         end do
!        Each one will be replicated nsym times in the end:
         do ij=1,6
           ggij(ij,iatom,ichan,iwchan)=ggij(ij,iatom,ichan,iwchan)/dble(Crystal%nsym)
         end do
       end do
     end do
   end do

   call wrtout(std_out,' harmonic_thermo: g(w) and gij(k|w) calculated given a q sampling grid.','COLL')

!  Sum up the counts in the channels to check the normalization
!  and test convergence with respect to q sampling
   gnorm=(igqpt2(1)*igqpt2(2)*igqpt2(3)*3*natom)/facbrv

   if(.not.part1)then
     do iwchan=nwchan,1,-1

       !write(msg,'(a,i0)' )' harmonic_thermo : iwchan=',iwchan
       !call wrtout(std_out,msg,'COLL')

       if (wgcnv(iwchan)) cycle
       !call wrtout(std_out,' harmonic_thermo : compute g, f, e, s, c','COLL')

!      Calculate g(w) and F,E,S,C

       ggsum=0.0_dp
       ggsumsum=0.0_dp
       ggrestsum=0.0_dp
       do ichan=1,nchan2(iwchan)
         ggsum=ggsum+gg(ichan,iwchan)
         ggsumsum=ggsumsum+gg_sum(ichan,iwchan)
         ggrestsum=ggrestsum+gg_rest(ichan,iwchan)
       end do

       if(ggsum/=gnorm)then
         write(msg, '(a,es14.6,i6,a,a,es14.6,a)' )&
&         'Frequencies are missing in g(w) : ggsum,iwchan=',ggsum,iwchan,ch10,&
&         'gnorm=',gnorm,'.'
         MSG_BUG(msg)
       end if

!      Check if the density of states changed by more than dostol

       gerr=zero
       if (ngrids>1) then
         do ichan=1,nchan2(iwchan)
           gerr=gerr+abs(gg(ichan,iwchan)/ggsum-gdos(ichan,iwchan))
         end do
       end if

       if(gerr>anaddb_dtset%dostol.and.ngrids>1) then
         wgcnv(iwchan)=.false.
       else
         wgcnv(iwchan)=.true.
       end if

!      g(w) is updated
       do ichan=1,nchan2(iwchan)
         gdos(ichan,iwchan)=gg(ichan,iwchan)/ggsum
         gg_sum(ichan,iwchan)=gg_sum(ichan,iwchan)/ggsumsum
         gg_rest(ichan,iwchan)=gg_rest(ichan,iwchan)/ggrestsum
       end do
       do ichan=1,nchan2(iwchan)
         gw(ichan,iwchan)=gw(ichan,iwchan)/ggsum
       end do

!      Write gerr for each q sampling and w width
       write(msg,'(a,a,i3,3i6,f10.1,f10.5)') ch10, &
&       'harmonic_thermo: iwchan,igqpt(:),norm,error=',iwchan,igqpt2(1),igqpt2(2),igqpt2(3),ggsum+tol8,gerr+tol10
       call wrtout(std_out,msg,'COLL')

!      If the DOS with a channel width is newly converged,
!      print it out and calculate the thermodynamic functions.
       convth=0
       if(wgcnv(iwchan)) then
         if (ngrids==1) then
           if (anaddb_dtset%dossum /= 0 ) then
             write(msg,'(a65,i5,a16)') ' DOS, SUM AND DIFFERENCE OF TWO PHONON DOS with channel width= ',iwchan,':'
           else
             write(msg,'(a25,i5,a16)') ' DOS  with channel width=  ',iwchan,':'
           end if
         else
           if (anaddb_dtset%dossum /= 0 ) then
             write(msg,'(a65,i5,a16)')&
&             ' DOS, SUM AND DIFFERENCE OF TWO PHONON DOS  with channel width= ',iwchan,' newly converged'
           else
             write(msg,'(a25,i5,a16)') ' DOS  with channel width=  ',iwchan,' newly converged'
           end if
         end if

         call wrtout(std_out,msg,'COLL')
         call wrtout(iout,msg,'COLL')
         do ichan=1,nchan2(iwchan)
           if (anaddb_dtset%dossum /= 0 ) then
             write(msg,'(i8,f11.1,3(f12.5,3x))') ichan,gg(ichan,iwchan)+tol10,&
&             gdos(ichan,iwchan)+tol10, gg_sum(ichan,iwchan)*(3.0*natom*3.0*natom)+tol10, &
&             gg_rest(ichan,iwchan)*(3.0*natom*(3.0*natom-1))+tol10
           else
             write(msg,'(i8,f11.1,1x,f12.5)') ichan,gg(ichan,iwchan)+tol10,&
&             gdos(ichan,iwchan)+tol10
           end if
           call wrtout(std_out,msg,'COLL')
         end do

         if (ngrids>1) then
           write(msg,'(a24,f10.5)')'   with maximal error = ',gerr+tol10
           call wrtout(std_out,msg,'COLL')
           call wrtout(iout,msg,'COLL')
         end if

         nomega = nchan2(iwchan)
         dosinc=dble(iwchan)

         ABI_ALLOCATE(phon_dos,(nomega))
         phon_dos = gdos(1:nomega,iwchan)

!Put zeroes for F, E, S, Cv
         free(:)=zero
         energy(:)=zero
         entropy(:)=zero
         spheat(:)=zero
         wme(:)=zero

         do itemper=1,ntemper

           tmp=anaddb_dtset%tempermin+anaddb_dtset%temperinc*dble(itemper-1)
!          The temperature (tmp) is given in Kelvin
           if (tmp < tol6) cycle

           do iomega=1,nomega

!            wovert= hbar*w / 2kT dimensionless
             wovert=dosinc*(dble(iomega)-0.5_dp)/Ha_cmm1/(2._dp*kb_HaK*tmp)
             expm2x=exp(-2.0_dp*wovert)
             ln2shx=wovert+log(1.0_dp-expm2x)
             cothx=(1.0_dp+expm2x)/(1.0_dp-expm2x)
             factor=dble(3*natom)*phon_dos(iomega)
             factorw=3*natom*gw(iomega,iwchan)

!            This matches the equations published in Lee & Gonze, PRB 51, 8610 (1995) [[cite:Lee1995]]
             free(itemper)=free(itemper) +factor*kb_HaK*tmp*ln2shx
             energy(itemper)=energy(itemper) + factor*kb_HaK*tmp*wovert*cothx
             entropy(itemper)=entropy(itemper) + factor*kb_HaK*(wovert*cothx - ln2shx)

!            The contribution is much lower than 1.0d-16 when wovert<100.0_dp
             if(wovert<100.0_dp)then
               spheat(itemper)=spheat(itemper)+factor*kb_HaK*wovert**2/sinh(wovert)**2
             end if
             wme(itemper)=wme(itemper)+factorw*kb_HaK*wovert**2/sinh(wovert)**2

           end do ! iomega

           if (abs(spheat(itemper))>tol8) wme(itemper)=wme(itemper)/spheat(itemper)
         end do ! itemper
         ABI_DEALLOCATE(phon_dos)

!        Check if the thermodynamic functions change within tolerance,
         if (ngrids>1) then
           write(msg,'(a,a,a)')&
&           ' harmonic_thermo : check if the thermodynamic functions',ch10,&
&           '    change within tolerance.'
           call wrtout(std_out,msg,'COLL')
           convth=1

           do itemper=1,ntemper
             change=free(itemper)-free0(itemper)
             relchg=change/free0(itemper)
             if(change>1d-14*dble(mqpt2) .and. relchg>thmtol)then
               write(msg,'(a,es14.4,a,a,es14.4)' )&
&               ' harmonic_thermo : free energy relative changes ',relchg,ch10,&
&               '        are larger than thmtol ',thmtol
               call wrtout(std_out,msg,'COLL')
               convth=0
             end if
             change=energy(itemper)-energy0(itemper)
             relchg=change/energy0(itemper)
             if(change>1d-14*dble(mqpt2) .and. relchg>thmtol)then
               write(msg,'(a,es14.4,a,a,es14.4)' )&
&               ' harmonic_thermo : energy relative changes ',relchg,ch10,&
&               '        are larger than thmtol ',thmtol
               call wrtout(std_out,msg,'COLL')
               convth=0
             end if
             change=entropy(itemper)-entropy0(itemper)
             relchg=change/entropy0(itemper)
             if(change>1d-14*dble(mqpt2) .and. relchg>thmtol)then
               write(msg,'(a,es14.4,a,a,es14.4)' )&
&               ' harmonic_thermo : entropy relative changes ',relchg,ch10,&
&               '        are larger than thmtol ',thmtol
               call wrtout(std_out,msg,'COLL')
               convth=0
             end if
             change=spheat(itemper)-spheat0(itemper)
             relchg=change/spheat0(itemper)
             if(change>1d-14*dble(mqpt2) .and. relchg>thmtol)then
               write(msg,'(a,es14.4,a,a,es14.4)' )&
&               ' harmonic_thermo : specific heat relative changes ',relchg,ch10,&
&               '        are larger than thmtol ',thmtol
               call wrtout(std_out,msg,'COLL')
               convth=0
             end if

             if(convth==0)exit
           end do ! End of check if the thermodynamic functions change within tolerance

         else
           convth=1
         end if

!        Update F,E,S,C and eventually write them if converged
         if(convth==1)then
           part1=.true.
           write(msg,'(a,a,a)') ch10,&
&           ' # At  T     F(J/mol-c)     E(J/mol-c)     S(J/(mol-c.K)) C(J/(mol-c.K)) Omega_mean(cm-1)'
           call wrtout(iout,msg,'COLL')
           call wrtout(thermal_unit,msg,'COLL')
           msg = ' # (A mol-c is the abbreviation of a mole-cell, that is, the'
           call wrtout(iout,msg,'COLL')
           call wrtout(thermal_unit,msg,'COLL')
           msg = ' #  number of Avogadro times the atoms in a unit cell)'
           call wrtout(iout,msg,'COLL')
           call wrtout(thermal_unit,msg,'COLL')

           write(msg, '(a,a,a)' )&
&           ' harmonic_thermo : thermodynamic functions have converged',ch10,&
&           '     see main output file ...'
           call wrtout(std_out,msg,'COLL')
         end if

         do itemper=1,ntemper
           free0(itemper)=free(itemper)
           energy0(itemper)=energy(itemper)
           entropy0(itemper)=entropy(itemper)
           spheat0(itemper)=spheat(itemper)

           if(convth==1)then
             tmp=anaddb_dtset%tempermin+anaddb_dtset%temperinc*dble(itemper-1)
             write(msg,'(es11.3,5es15.7)') tmp+tol8,&
&             Ha_eV*e_Cb*Avogadro*free(itemper),&
&             Ha_eV*e_Cb*Avogadro*energy(itemper),&
&             Ha_eV*e_Cb*Avogadro*entropy(itemper),&
&             Ha_eV*e_Cb*Avogadro*spheat(itemper),&
&             wme(itemper)
             call wrtout(iout,msg,'COLL')
             call wrtout(thermal_unit,msg,'COLL')
           end if
         end do
       end if

       if(convth==1)exit
     end do
   end if

   if(.not.part2)then
     ! Atomic temperature factor calculation
     do iwchan=nwchan,1,-1
       if (wgijcnv(iwchan))cycle

       ! Calculate gij(k|w) and Bij(k)
       ! Check if the density of states changed by more than dostol
       gijsum =zero
       wgijcnv(iwchan)=.true.
       if (ngrids>1) then
         do iatom=1,natom
           do ij=1,6
             gijerr=zero
             do ichan=1,nchan2(iwchan)
               gijsum = gijsum + gij(ij,iatom,ichan,iwchan)
               gijerr=gijerr&
&               +abs(ggij(ij,iatom,ichan,iwchan)/gnorm&
&               -     gij(ij,iatom,ichan,iwchan))
             end do
             if(gijerr>anaddb_dtset%dostol) then
               wgijcnv(iwchan)=.false.
               exit
             end if
           end do
         end do
       else
         gijerr=0.d0
       end if

!      gij(k|w) is updated

       do ichan=1,nchan2(iwchan)
         do iatom=1,natom
           do ij=1,6
             gij(ij,iatom,ichan,iwchan)=ggij(ij,iatom,ichan,iwchan)/(gnorm/(3*natom))
           end do
!if (iwchan==1) write (200+iatom,'(I6,6(E20.10,2x))') ichan, gij(1:6,iatom,ichan,iwchan)
         end do
       end do

!      Write gijerr for each q sampling and w width

       write(msg,'(a,a,i3,3i6,f10.5,f10.5)') ch10,&
&       ' iwchan,igqpt(i),gijsum, gij error= ',&
&       iwchan,igqpt2(1),igqpt2(2),igqpt2(3),gijsum,gijerr+tol10
       call wrtout(std_out,msg,'COLL')

!      If the generalized DOS with a channel width is newly converged,
!      print it out and calculate Bij(k).
       if(wgijcnv(iwchan)) then

         if (ngrids==1) then
           write(msg,'(a,i5,a)') ' gij with channel width=  ',iwchan,':'
         else
           write(msg,'(a,i5,a)') ' gij with channel width=  ',iwchan,' newly converged'
         end if
         call wrtout(iout,msg,'COLL')

         write(msg,'(a,2i3,3i6,f10.5)')'iatom,iwchan,igqpt2(i),gij error= ',&
&         iatom,iwchan,igqpt2(1),igqpt2(2),igqpt2(3),gijerr+tol10
         call wrtout(iout,msg,'COLL')

         do itemper=1,ntemper

!          Put zeroes for Bij(k)
           do iatom=1,natom
             do ij=1,6
               bbij(ij,iatom,itemper)=0._dp
               vij(ij,iatom,itemper)=0._dp
             end do
           end do

           tmp=anaddb_dtset%tempermin+anaddb_dtset%temperinc*dble(itemper-1)
!          tmp in K
           if (tmp < tol6) cycle

           dosinc=dble(iwchan)
!
           do ichan=1,nchan2(iwchan)
!
!$wovert= \hbar*w / 2kT$, dimensionless
             wovert=dosinc*(dble(ichan)-half)/Ha_cmm1/(two*kb_HaK*tmp)
             expm2x=exp(-two*wovert)
             do iatom=1,natom
!   factor contains 1 / (2 omega)
               factor=Ha_cmm1/(two*dosinc*(dble(ichan)-half))    &
&               *(one+expm2x)/(one-expm2x) /amu(Crystal%typat(iatom))/amu_emass

!   this becomes * 0.5 * omega for the velocities
               factorv=(half*dosinc*(dble(ichan)-half)/Ha_cmm1)    &
&               *(one+expm2x)/(one-expm2x) /amu(Crystal%typat(iatom))/amu_emass

               do ij=1,6
                 bbij(ij,iatom,itemper)=bbij(ij,iatom,itemper) + factor*gij(ij,iatom,ichan,iwchan)
                 vij(ij,iatom,itemper)=vij(ij,iatom,itemper) + factorv*gij(ij,iatom,ichan,iwchan)
               end do
             end do

           end do

         end do

!        B matrix is now in atomic unit in the Cartesian coordinates.
!        Check if Bij(k) changed within tolerance.
         convth=1
         if (ngrids>1) then
           do itemper=1,ntemper
             do iatom=1,natom
               do ij=1,6
                 diffbb=bbij(ij,iatom,itemper)-bij(ij,iatom,itemper)
                 if (diffbb > 1d-10  .and. diffbb/bij(ij,iatom,itemper) > thmtol) then
                   write(msg,'(a)' )' harmonic_thermo : Bij changes are larger than thmtol '
                   call wrtout(std_out,msg,'COLL')
                   convth=0
                 end if
                 if(convth==0)exit
               end do
               if(convth==0)exit
             end do
             if(convth==0)exit
           end do
         end if

         bij=bbij ! save for next iteration

         !Update Bij(k) and write them. B matrix printed in angstrom^2
         !TODO : get rid of this version in the log and output file. Prefer
         !external files
         if (convth==1) then
           write(msg, '(a,a,a)' )&
&           ' B matrix elements as a function of T',ch10,&
&           '    Angstrom^2, cartesian coordinates'
           call wrtout(std_out,msg,'COLL')
           call wrtout(iout,msg,'COLL')

           do itemper=1,ntemper
!            tmp in K
             tmp=anaddb_dtset%tempermin+anaddb_dtset%temperinc*dble(itemper-1)
             do iatom=1,natom
               write(iout,'(2i3,es11.3,6es12.4)')&
&               iwchan,iatom,tmp+tol10,&
&               Bohr_Ang**2*bij(1,iatom,itemper)+tol10,&
&               Bohr_Ang**2*bij(2,iatom,itemper)+tol10,&
&               Bohr_Ang**2*bij(3,iatom,itemper)+tol10,&
&               Bohr_Ang**2*bij(4,iatom,itemper)+tol10,&
&               Bohr_Ang**2*bij(5,iatom,itemper)+tol10,&
&               Bohr_Ang**2*bij(6,iatom,itemper)+tol10
             end do ! end loop over natom
           end do ! end loop over ntemper

!        Mean square velocity matrix printed in angstrom^2/picosec^2
           write(msg, '(a,a,a)' )&
&           ' <vel^2> matrix elements as a function of T',ch10,&
&           '    Angstrom^2/(picosec)^2, cartesian coordinates'
           call wrtout(std_out,msg,'COLL')
           call wrtout(iout,msg,'COLL')

           do itemper=1,ntemper
!            tmp in K
             tmp=anaddb_dtset%tempermin+anaddb_dtset%temperinc*float(itemper-1)
             do iatom=1,natom
               vij(:,iatom,itemper)=Bohr_Ang**2*vij(:,iatom,itemper)/(Time_Sec*1.0e12)**2
!              The following check zeros out <v^2> if it is very small, in order to
!              avoid numerical noise being interpreted by the automatic tests as
!              something real. Note also that we compare it in
!              absolute value, that's because if any of the phonon frequencies are
!              computed as negative, <v^2> can take a negative value.
               do icomp=1, 6
                 if (abs(vij(icomp,iatom,itemper)) < 1.0e-12) vij(icomp,iatom,itemper)=zero
               end do
               write(iout,'(2i3,es11.3,6es12.4)')&
&               iwchan,iatom,tmp+tol10,&
&               vij(1,iatom,itemper),&
&               vij(2,iatom,itemper),&
&               vij(3,iatom,itemper),&
&               vij(4,iatom,itemper),&
&               vij(5,iatom,itemper),&
&               vij(6,iatom,itemper)
             end do ! end loop over natom
           end do ! end loop over ntemper
         end if ! end check on convergence


         ! keep this one !!!!!!!!!!!!!!!!!!
         if (convth==1) then
           write(msg, '(a,a,a)' )&
&           '# B matrix elements as a function of T, for each atom, and smallest omega channel width',ch10,&
&           '#    Angstrom^2, cartesian coordinates'
           call wrtout(bij_unit,msg,'COLL')
           do iatom=1,natom
             write(msg, '(2a,i10)' ) ch10, '# for atom ', iatom
             call wrtout(bij_unit,msg,'COLL')
             do itemper=1,ntemper
!              tmp in K
               tmp=anaddb_dtset%tempermin+anaddb_dtset%temperinc*dble(itemper-1)
               write(msg,'(es11.3,6es12.4)')&
&               tmp,&
&               Bohr_Ang**2*bij(1,iatom,itemper),&
&               Bohr_Ang**2*bij(2,iatom,itemper),&
&               Bohr_Ang**2*bij(3,iatom,itemper),&
&               Bohr_Ang**2*bij(4,iatom,itemper),&
&               Bohr_Ang**2*bij(5,iatom,itemper),&
&               Bohr_Ang**2*bij(6,iatom,itemper)
               call wrtout(bij_unit,msg,'COLL')
             end do ! end loop over ntemper
           end do ! end loop over natom

!        Mean square velocity matrix printed in angstrom^2/picosec^2
           write(msg, '(a,a,a)' )&
&           '# <vel^2> matrix elements as a function of T, for each atom, and smallest channel width',ch10,&
&           '#    Angstrom^2/(picosec)^2, cartesian coordinates'
           call wrtout(vij_unit,msg,'COLL')

           do iatom=1,natom
             write(msg, '(2a,i10)' ) ch10, '# for atom ', iatom
             call wrtout(vij_unit,msg,'COLL')
             do itemper=1,ntemper
!            tmp in K
               tmp=anaddb_dtset%tempermin+anaddb_dtset%temperinc*float(itemper-1)
               vij(:,iatom,itemper)=Bohr_Ang**2*vij(:,iatom,itemper)/(Time_Sec*1.0e12)**2

!            The following check zeros out <v^2> if it is very small, in order to
!            avoid numerical noise being interpreted by the automatic tests as
!            something real. Note also that we compare it in
!            absolute value, that's because if any of the phonon frequencies are
!            computed as negative, <v^2> can take a negative value.
               do icomp=1, 6
                 if (abs(vij(icomp,iatom,itemper)) < 1.0e-12) vij(icomp,iatom,itemper)=zero
               end do
               write(vij_unit,'(es11.3,6es12.4)')&
&               tmp,&
&               vij(1,iatom,itemper),&
&               vij(2,iatom,itemper),&
&               vij(3,iatom,itemper),&
&               vij(4,iatom,itemper),&
&               vij(5,iatom,itemper),&
&               vij(6,iatom,itemper)
             end do ! end loop over ntemper
           end do ! end loop over natom
         end if ! end check on convergence

         if(convth==1)part2=.true.

       end if ! End of test on wgijcnv
     end do ! End of loop over iwchan
   end if ! End of part2

   if(part1.and.part2)exit

   ABI_DEALLOCATE(indqpt1)
   ABI_DEALLOCATE(qpt2)
   ABI_DEALLOCATE(spqpt2)
   ABI_DEALLOCATE(wtq)
   ABI_DEALLOCATE(wtq2)

 end do ! End of the Loop on the q point grids

 ABI_DEALLOCATE(bbij)
 ABI_DEALLOCATE(bij)
 ABI_DEALLOCATE(energy)
 ABI_DEALLOCATE(energy0)
 ABI_DEALLOCATE(entropy)
 ABI_DEALLOCATE(entropy0)
 ABI_DEALLOCATE(free)
 ABI_DEALLOCATE(free0)
 ABI_DEALLOCATE(gdos)
 ABI_DEALLOCATE(gg)
 ABI_DEALLOCATE(gg_sum)
 ABI_DEALLOCATE(gg_rest)
 ABI_DEALLOCATE(ggij)
 ABI_DEALLOCATE(gij)
 ABI_DEALLOCATE(nchan2)
 ABI_DEALLOCATE(spheat)
 ABI_DEALLOCATE(spheat0)
 ABI_DEALLOCATE(vij)
 ABI_DEALLOCATE(wgcnv)
 ABI_DEALLOCATE(wgijcnv)
 if(allocated(indqpt1)) then
   ABI_DEALLOCATE(indqpt1)
 end if
 if(allocated(qpt2)) then
   ABI_DEALLOCATE(qpt2)
 end if
 if(allocated(spqpt2)) then
   ABI_DEALLOCATE(spqpt2)
 end if
 if(allocated(wtq)) then
   ABI_DEALLOCATE(wtq)
 end if
 if(allocated(wtq2)) then
   ABI_DEALLOCATE(wtq2)
 end if
 ABI_DEALLOCATE(gw)
 ABI_DEALLOCATE(wme)

 if(.not.part1)then
   write(msg, '(a,a,a,a,a,a,a,a,a)' )&
&   'No thermodynamical function is printed out :',ch10,&
&   'the tolerance level that was asked ',ch10,&
&   'has not been match with the grids specified.',ch10,&
&   'Action: in the input file, increase the resolution',ch10,&
&   'of grids ng2qpt, or decrease the accuracy requirement thmtol.'
   MSG_ERROR(msg)
 end if

 if(.not.part2)then
   write(msg,'(a,a,a,a,a,a,a,a,a)')&
&   'No atomic factor tensor is printed out :',ch10,&
&   'the tolerance level that was asked ',ch10,&
&   'has not been match with the grids specified.',ch10,&
&   'Action: in the input file, increase the resolution',ch10,&
&   'of grids ng2qpt, or decrease the accuracy requirement thmtol.'
   MSG_WARNING(msg)
 end if

 close (thermal_unit)
 close (bij_unit)
 close (vij_unit)

end subroutine harmonic_thermo
!!***

end module m_harmonic_thermo
!!***
