!{\src2tex{textfont=tt}}
!!****f* ABINIT/initorbmag
!! NAME
!! initorbmag
!!
!! FUNCTION
!! Initialization of orbital magnetization calculation; similar to initberry
!!
!! COPYRIGHT
!! Copyright (C) 2004-2017 ABINIT group.
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables in this dataset
!!  gmet(3,3) = reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3) = primitive translations in recip space
!!  kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!!  npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!!  occ(mband*nkpt*nsppol) = occup number for each band at each k point
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3) = dimensional primitive vectors
!!  symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!    reciprocal space primitive translations
!!  xred(3,natom) = location of atoms in reduced units
!!
!! OUTPUT
!!  dtorbmag <type(orbmag_type)> = variables related to orbital magnetization
!!
!! SIDE EFFECTS
!!  mpi_enreg = information about MPI parallelization
!!
!! TO DO
!!
!! NOTES
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

subroutine initorbmag(dtorbmag,dtset,gmet,gprimd,kg,mpi_enreg,npwarr,occ,&
&                     pawang,pawrad,pawtab,psps,pwind,pwind_alloc,pwnsfac,&
&                     rprimd,symrec,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_orbmag
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_fftcore, only : kpgsph
 use m_pawang,  only : pawang_type
 use m_pawrad,  only : pawrad_type
 use m_pawtab,  only : pawtab_type
 use m_pawcprj, only : pawcprj_alloc, pawcprj_getdim

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initorbmag'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_56_recipspace
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 !scalars
 integer,intent(out) :: pwind_alloc
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(inout) :: dtset
 type(orbmag_type),intent(inout) :: dtorbmag
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 !arrays
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),npwarr(dtset%nkpt)
 integer,intent(in) :: symrec(3,3,dtset%nsym)
 integer,pointer :: pwind(:,:,:)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: rprimd(3,3),xred(3,dtset%natom)
 real(dp),pointer :: pwnsfac(:,:)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer :: brav,exchn2n3d,fnkpt_computed
 integer :: iband,icg,icprj,idir,idum,idum1,ierr,ifor,ikg,ikg1
 integer :: ikpt,ikpt_loc,ikpti,ikpt1,ikpt1f,ikpt1i
 integer :: index,ipw,ipwnsfac,isign,isppol,istwf_k,isym,isym1,itrs,itypat
 integer :: jpw,lmax,lmn2_size_max
 integer :: mband_occ_k,me,me_g0,mkmem_,mkpt,my_nspinor,nband_k,nkptlatt,nproc,npw_k,npw_k1
 integer :: option,spaceComm
 real(dp) :: diffk1,diffk2,diffk3,ecut_eff,kpt_shifted1,kpt_shifted2,kpt_shifted3,rdum
 character(len=500) :: message
 !arrays
 integer :: iadum(3),iadum1(3),dg(3)
 integer,allocatable :: kg1_k(:,:)
 real(dp) :: diffk(3),dk(3),dum33(3,3),kpt1(3),tsec(2)
 real(dp),allocatable :: spkpt(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(1001,1,tsec)
 call timab(1002,1,tsec)

 !save the current value of nspinor
 dtorbmag%nspinor = dtset%nspinor

!----------------------------------------------------------------------------
!-------------------- Obtain k-point grid in the full BZ --------------------
!----------------------------------------------------------------------------

 if(dtset%kptopt==1 .or. dtset%kptopt==2 .or. dtset%kptopt==4)then
!  Compute the number of k points in the G-space unit cell
   nkptlatt=dtset%kptrlatt(1,1)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,3) &
&   +dtset%kptrlatt(1,2)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,1) &
&   +dtset%kptrlatt(1,3)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,2) &
&   -dtset%kptrlatt(1,2)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,3) &
&   -dtset%kptrlatt(1,3)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,1) &
&   -dtset%kptrlatt(1,1)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,2)

!  Call smpbz to obtain the list of k-point in the full BZ - without symmetry reduction
   option = 0
   brav = 1
   mkpt=nkptlatt*dtset%nshiftk
   ABI_ALLOCATE(spkpt,(3,mkpt))
   call smpbz(1,ab_out,dtset%kptrlatt,mkpt,fnkpt_computed,dtset%nshiftk,option,dtset%shiftk,spkpt)
   dtorbmag%fnkpt = fnkpt_computed
   ABI_ALLOCATE(dtorbmag%fkptns,(3,dtorbmag%fnkpt))
   dtorbmag%fkptns(:,:)=spkpt(:,1:dtorbmag%fnkpt)
   ABI_DEALLOCATE(spkpt)
 else if(dtset%kptopt==3.or.dtset%kptopt==0)then
   dtorbmag%fnkpt=dtset%nkpt
   ABI_ALLOCATE(dtorbmag%fkptns,(3,dtorbmag%fnkpt))
   dtorbmag%fkptns(1:3,1:dtorbmag%fnkpt)=dtset%kpt(1:3,1:dtorbmag%fnkpt)
   if(dtset%kptopt==0)then
     write(message,'(10a)') ch10,&
&     ' initorbmag : WARNING -',ch10,&
&     '  you have defined manually the k-point grid with kptopt = 0',ch10,&
&     '  the orbital magnetization calculation works only with a regular k-points grid,',ch10,&
&     '  abinit doesn''t check if your grid is regular...'
     call wrtout(std_out,message,'PERS')
   end if
 end if

!call listkk to get mapping from FBZ to IBZ
 rdum=1.0d-5  ! cutoff distance to decide when two k points match
 ABI_ALLOCATE(dtorbmag%indkk_f2ibz,(dtorbmag%fnkpt,6))

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

!JWZ: The following may need modification in the future
!**** no spin-polarization doubling ; do not allow use of time reversal symmetry ****

 call timab(1002,2,tsec)
 call timab(1003,1,tsec)

 call listkk(rdum,gmet,dtorbmag%indkk_f2ibz,dtset%kptns,dtorbmag%fkptns,dtset%nkpt,&
& dtorbmag%fnkpt,dtset%nsym,1,dtset%symafm,symrec,0,use_symrec=.True.)

 call timab(1003,2,tsec)
 call timab(1004,1,tsec)

!Construct i2fbz and f2ibz
 ABI_ALLOCATE(dtorbmag%i2fbz,(dtset%nkpt))
 idum=0
 do ikpt=1,dtorbmag%fnkpt
   if (dtorbmag%indkk_f2ibz(ikpt,2)==1 .and. &
&   dtorbmag%indkk_f2ibz(ikpt,6) == 0 .and. &
&   maxval(abs(dtorbmag%indkk_f2ibz(ikpt,3:5))) == 0 ) then
     dtorbmag%i2fbz(dtorbmag%indkk_f2ibz(ikpt,1))=ikpt
     idum=idum+1
   end if
 end do
 if (idum/=dtset%nkpt)then
   message = ' Found wrong number of k-points in IBZ'
   MSG_ERROR(message)
 end if

!----------------------------------------------------------------------------
!------------- Allocate PAW space as necessary ------------------------------
!----------------------------------------------------------------------------

 dtorbmag%usepaw   = psps%usepaw
 dtorbmag%natom    = dtset%natom
 dtorbmag%my_natom = mpi_enreg%my_natom

 ABI_ALLOCATE(dtorbmag%lmn_size,(dtset%ntypat))
 ABI_ALLOCATE(dtorbmag%lmn2_size,(dtset%ntypat))
 do itypat = 1, dtset%ntypat
    dtorbmag%lmn_size(itypat) = pawtab(itypat)%lmn_size
    dtorbmag%lmn2_size(itypat) = pawtab(itypat)%lmn2_size
 end do

 lmn2_size_max = psps%lmnmax*(psps%lmnmax+1)/2
 dtorbmag%lmn2max = lmn2_size_max

! expibi and qijb_kk are NOT parallelized over atoms
! this may change in the future (JZwanziger 18 March 2014)
 ABI_ALLOCATE(dtorbmag%qijb_kk,(2,lmn2_size_max,dtorbmag%natom,3))
 ABI_ALLOCATE(dtorbmag%expibi,(2,dtorbmag%natom,3))
 dtorbmag%has_expibi = 1
 dtorbmag%has_qijb = 1

 ABI_ALLOCATE(dtorbmag%cprjindex,(dtset%nkpt,dtset%nsppol))
 dtorbmag%cprjindex(:,:) = 0

 if (dtset%kptopt /= 3) then
    ABI_ALLOCATE(dtorbmag%atom_indsym,(4,dtset%nsym,dtorbmag%natom))
    call symatm(dtorbmag%atom_indsym,dtorbmag%natom,dtset%nsym,symrec,dtset%tnons,tol8,dtset%typat,xred)
    lmax = psps%mpsang - 1
    ABI_ALLOCATE(dtorbmag%zarot,(2*lmax+1,2*lmax+1,lmax+1,dtset%nsym))
    call setsymrhoij(gprimd,lmax,dtset%nsym,1,rprimd,symrec,dtorbmag%zarot)
    dtorbmag%nsym = dtset%nsym
    dtorbmag%lmax = lmax
    dtorbmag%lmnmax = psps%lmnmax
 end if

! !------------------------------------------------------------------------------
! !------------------- Compute variables related to MPI // ----------------------
! !------------------------------------------------------------------------------
 spaceComm=mpi_enreg%comm_cell
 nproc=xmpi_comm_size(spaceComm)
 me=xmpi_comm_rank(spaceComm)

 if (nproc==1) then
   dtorbmag%fmkmem = dtorbmag%fnkpt
   dtorbmag%fmkmem_max = dtorbmag%fnkpt
   dtorbmag%mkmem_max = dtset%nkpt
 else
   dtorbmag%fmkmem = 0
   do ikpt = 1, dtorbmag%fnkpt
     ikpti = dtorbmag%indkk_f2ibz(ikpt,1)
     nband_k = dtset%nband(ikpti)
     if (.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,-1,me))) &
&     dtorbmag%fmkmem = dtorbmag%fmkmem + 1
   end do
!  Maximum value of mkmem and fmkmem
   call xmpi_max(dtorbmag%fmkmem,dtorbmag%fmkmem_max,spaceComm,ierr)
!  I have to use the dummy variable mkmem_ because
!  mkmem is declared as intent(in) while the first
!  argument of xmpi_max must be intent(inout)
   mkmem_ = dtset%mkmem
   call xmpi_max(mkmem_,dtorbmag%mkmem_max,spaceComm,ierr)
 end if

 ABI_ALLOCATE(mpi_enreg%kpt_loc2fbz_sp,(0:nproc-1,1:dtorbmag%fmkmem_max*dtset%nsppol, 1:2))
 ABI_ALLOCATE(mpi_enreg%kpt_loc2ibz_sp,(0:nproc-1,1:dtorbmag%mkmem_max*dtset%nsppol, 1:2))
 ABI_ALLOCATE(mpi_enreg%kptdstrb,(nproc,6,dtorbmag%fmkmem_max*dtset%nsppol*2))
 ABI_ALLOCATE(mpi_enreg%mkmem,(0:nproc-1))
 mpi_enreg%kpt_loc2fbz_sp(:,:,:) = 0
 mpi_enreg%kpt_loc2ibz_sp(:,:,:) = 0
 mpi_enreg%kptdstrb(:,:,:)       = 0
 mpi_enreg%mkmem(:)              = 0

 pwind_alloc = dtset%mpw*dtorbmag%fmkmem_max
 ABI_ALLOCATE(pwind,(pwind_alloc,2,3))
 ABI_ALLOCATE(pwnsfac,(2,pwind_alloc))

! !------------------------------------------------------------------------------
! !---------------------- Compute orbmag_type variables -------------------------
! !------------------------------------------------------------------------------

 !Initialization of orbmag_type variables
 dtorbmag%dkvecs(:,:) = zero
 ABI_ALLOCATE(dtorbmag%ikpt_dk,(dtorbmag%fnkpt,2,3))
 ABI_ALLOCATE(dtorbmag%cgindex,(dtset%nkpt,dtset%nsppol))
 ABI_ALLOCATE(dtorbmag%kgindex,(dtset%nkpt))
 ABI_ALLOCATE(dtorbmag%fkgindex,(dtorbmag%fnkpt))
 dtorbmag%ikpt_dk(:,:,:) = 0
 dtorbmag%cgindex(:,:) = 0
 dtorbmag%mband_occ = 0
 ABI_ALLOCATE(dtorbmag%nband_occ,(dtset%nsppol))
 dtorbmag%kgindex(:) = 0
 dtorbmag%fkgindex(:) = 0

!Compute spin degeneracy
 if (dtset%nsppol == 1 .and. dtset%nspinor == 1) then
   dtorbmag%sdeg = two
 else if (dtset%nsppol == 2 .or. my_nspinor == 2) then
   dtorbmag%sdeg = one
 end if

!Compute the number of occupied bands and check that
!it is the same for each k-point

 index = 0
 do isppol = 1, dtset%nsppol
   dtorbmag%nband_occ(isppol) = 0
   do ikpt = 1, dtset%nkpt

     mband_occ_k = 0
     nband_k = dtset%nband(ikpt + (isppol - 1)*dtset%nkpt)

     do iband = 1, nband_k
       index = index + 1
       if (abs(occ(index) - dtorbmag%sdeg) < tol8) mband_occ_k = mband_occ_k + 1
     end do

     if (ikpt > 1) then
       if (dtorbmag%nband_occ(isppol) /= mband_occ_k) then
         message = "The number of valence bands is not the same for every k-point of present spin channel"
         MSG_ERROR(message)
       end if
     else
       dtorbmag%mband_occ         = max(dtorbmag%mband_occ, mband_occ_k)
       dtorbmag%nband_occ(isppol) = mband_occ_k
     end if

   end do                ! close loop over ikpt
 end do                ! close loop over isppol

!Compute the location of each wavefunction

 icg = 0
 icprj = 0
!ikg = 0
 do isppol = 1, dtset%nsppol
   do ikpt = 1, dtset%nkpt

     nband_k = dtset%nband(ikpt + (isppol-1)*dtset%nkpt)

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle
     
     dtorbmag%cgindex(ikpt,isppol) = icg
     npw_k = npwarr(ikpt)
     icg = icg + npw_k*dtorbmag%nspinor*nband_k

     if (psps%usepaw == 1) then
       dtorbmag%cprjindex(ikpt,isppol) = icprj
       icprj = icprj + dtorbmag%nspinor*nband_k
     end if

   end do
 end do

 ikg = 0
 do ikpt = 1, dtset%nkpt
   if ((proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,1,me)).and.&
&   (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,dtset%nsppol,me))) cycle
   
   npw_k = npwarr(ikpt)
   dtorbmag%kgindex(ikpt) = ikg
   ikg = ikg + npw_k
 end do

 call timab(1004,2,tsec)

!------------------------------------------------------------------------------
!---------------------- Compute dk --------------------------------------------
!------------------------------------------------------------------------------

 call timab(1005,1,tsec)

 do idir = 1, 3

    !    Compute dk(:), the vector between a k-point and its nearest
    !    neighbour along the direction idir

    dk(:) = zero
    dk(idir) = 1._dp   ! 1 mean there is no other k-point un the direction idir
    do ikpt = 2, dtorbmag%fnkpt
       diffk(:) = abs(dtorbmag%fkptns(:,ikpt) - dtorbmag%fkptns(:,1))
       if ((diffk(1) < dk(1)+tol8).and.(diffk(2) < dk(2)+tol8).and.&
&          (diffk(3) < dk(3)+tol8)) dk(:) = diffk(:)
    end do
    dtorbmag%dkvecs(:,idir) = dk(:)
    !    DEBUG
    !    write(std_out,*)' initorbmag : idir, dk', idir, dk
    !    ENDDEBUG

    !    For each k point, find k_prim such that k_prim= k + dk mod(G)
    !    where G is a vector of the reciprocal lattice

    do ikpt = 1, dtorbmag%fnkpt

       !      First k+dk, then k-dk
       do isign=-1,1,2
          kpt_shifted1=dtorbmag%fkptns(1,ikpt)- isign*dk(1)
          kpt_shifted2=dtorbmag%fkptns(2,ikpt)- isign*dk(2)
          kpt_shifted3=dtorbmag%fkptns(3,ikpt)- isign*dk(3)
          !        Note that this is still a order fnkpt**2 algorithm.
          !        It is possible to implement a order fnkpt algorithm, see listkk.F90.
          do ikpt1 = 1, dtorbmag%fnkpt
             diffk1=dtorbmag%fkptns(1,ikpt1) - kpt_shifted1
             if(abs(diffk1-nint(diffk1))>tol8)cycle
             diffk2=dtorbmag%fkptns(2,ikpt1) - kpt_shifted2
             if(abs(diffk2-nint(diffk2))>tol8)cycle
             diffk3=dtorbmag%fkptns(3,ikpt1) - kpt_shifted3
             if(abs(diffk3-nint(diffk3))>tol8)cycle
             dtorbmag%ikpt_dk(ikpt,(isign+3)/2,idir) = ikpt1
             exit
          end do   ! ikpt1
       end do     ! isign

    end do     ! ikpt

 end do     ! close loop over idir

 call timab(1005,2,tsec)
 call timab(1006,1,tsec)

!------------------------------------------------------------------------------
!------------ Compute PAW on-site terms if necessary --------------------------
!------------------------------------------------------------------------------

 if (dtorbmag%usepaw == 1 .and. dtorbmag%has_expibi == 1) then
   call expibi(dtorbmag%expibi,dtorbmag%dkvecs,dtset%natom,xred)
   dtorbmag%has_expibi = 2
 end if

 if (dtorbmag%usepaw == 1 .and. dtorbmag%has_qijb == 1) then

   call qijb_kk(dtorbmag%qijb_kk,dtorbmag%dkvecs,dtorbmag%expibi,&
&   gprimd,dtorbmag%lmn2max,dtset%natom,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat)

 end if
 
 call timab(1007,2,tsec)
 call timab(1008,1,tsec)

!------------------------------------------------------------------------------
!------------ Build the array pwind that is needed to compute the -------------
!------------ overlap matrices at k +- dk                         -------------
!------------------------------------------------------------------------------

 ecut_eff = dtset%ecut*(dtset%dilatmx)**2
 exchn2n3d = 0 ; istwf_k = 1 ; ikg1 = 0
 pwind(:,:,:) = 0
 pwnsfac(1,:) = 1.0_dp
 pwnsfac(2,:) = 0.0_dp
 ABI_ALLOCATE(kg1_k,(3,dtset%mpw))

 ipwnsfac = 0

 do idir = 1, 3

    dk(:) = dtorbmag%dkvecs(:,idir)

    do ifor = 1, 2

       if (ifor == 2) dk(:) = -1._dp*dk(:)

       !      Build pwind and kgindex
       !      NOTE: The array kgindex is important for parallel execution.
       !      In case nsppol = 2, it may happen that a particular processor
       !      treats k-points at different spin polarizations.
       !      In this case, it is not possible to address the elements of
       !      pwind correctly without making use of the kgindex array.
       
       ikg = 0 ; ikpt_loc = 0 ; isppol = 1
       do ikpt = 1, dtorbmag%fnkpt

          ikpti = dtorbmag%indkk_f2ibz(ikpt,1)
          nband_k = dtset%nband(ikpti)
          ikpt1f = dtorbmag%ikpt_dk(ikpt,ifor,idir)
          ikpt1i = dtorbmag%indkk_f2ibz(ikpt1f,1)

          if ((proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,1,me)).and.&
&             (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,dtset%nsppol,me))) cycle

          ikpt_loc = ikpt_loc + 1

          !        Build basis sphere of plane waves for the nearest neighbour of
          !        the k-point (important for MPI //)

          kg1_k(:,:) = 0
          kpt1(:) = dtset%kptns(:,ikpt1i)
          call kpgsph(ecut_eff,exchn2n3d,gmet,ikg1,ikpt,istwf_k,kg1_k,kpt1,&
&                     1,mpi_enreg,dtset%mpw,npw_k1)
          me_g0=mpi_enreg%me_g0


          !        ji: fkgindex is defined here !
          dtorbmag%fkgindex(ikpt) = ikg

          !        
          !        Deal with symmetry transformations
          !        

          !        bra k-point k(b) and IBZ k-point kIBZ(b) related by
          !        k(b) = alpha(b) S(b)^t kIBZ(b) + G(b)
          !        where alpha(b), S(b) and G(b) are given by indkk_f2ibz
          !        
          !        For the ket k-point:
          !        k(k) = alpha(k) S(k)^t kIBZ(k) + G(k) - GBZ(k)
          !        where GBZ(k) takes k(k) to the BZ
          !        

          isym  = dtorbmag%indkk_f2ibz(ikpt,2)
          isym1 = dtorbmag%indkk_f2ibz(ikpt1f,2)

          !        Construct transformed G vector that enters the matching condition:
          !        alpha(k) S(k)^{t,-1} ( -G(b) - GBZ(k) + G(k) )

          dg(:) = -dtorbmag%indkk_f2ibz(ikpt,3:5) &
&                 -nint(-dtorbmag%fkptns(:,ikpt) - dk(:) - tol10 &
&                 +dtorbmag%fkptns(:,ikpt1f)) &
&                 +dtorbmag%indkk_f2ibz(ikpt1f,3:5)

          iadum(:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),dg(:))

          dg(:) = iadum(:)

          if ( dtorbmag%indkk_f2ibz(ikpt1f,6) == 1 ) dg(:) = -dg(:)

          !        Construct S(k)^{t,-1} S(b)^{t}

          dum33(:,:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),symrec(:,:,isym))

          !        Construct alpha(k) alpha(b)

          if (dtorbmag%indkk_f2ibz(ikpt,6) == dtorbmag%indkk_f2ibz(ikpt1f,6)) then
             itrs=0
          else
             itrs=1
          end if


          npw_k  = npwarr(ikpti)
          !        npw_k1 = npwarr(ikpt1i)

          !        loop over bra G vectors
          do ipw = 1, npw_k

             !          NOTE: the bra G vector is taken for the sym-related IBZ k point,
             !          not for the FBZ k point
             iadum(:) = kg(:,dtorbmag%kgindex(ikpti) + ipw)

             !          Store non-symmorphic operation phase factor exp[i2\pi \alpha G \cdot t]

             if ( ipwnsfac == 0 ) then
                rdum=0.0_dp
                do idum=1,3
                   rdum=rdum+dble(iadum(idum))*dtset%tnons(idum,isym)
                end do
                rdum=two_pi*rdum
                if ( dtorbmag%indkk_f2ibz(ikpt,6) == 1 ) rdum=-rdum
                pwnsfac(1,ikg+ipw) = cos(rdum)
                pwnsfac(2,ikg+ipw) = sin(rdum)
             end if

             !          to determine r.l.v. matchings, we transformed the bra vector
             !          Rotation
             iadum1(:)=0
             do idum1=1,3
                iadum1(:)=iadum1(:)+dum33(:,idum1)*iadum(idum1)
             end do
             iadum(:)=iadum1(:)
             !          Time reversal
             if (itrs==1) iadum(:)=-iadum(:)
             !          Translation
             iadum(:) = iadum(:) + dg(:)

             do jpw = 1, npw_k1
                iadum1(1:3) = kg1_k(1:3,jpw)
                if ( (iadum(1) == iadum1(1)).and. &
&                    (iadum(2) == iadum1(2)).and. &
&                    (iadum(3) == iadum1(3)) ) then
                   pwind(ikg + ipw,ifor,idir) = jpw
                   !              write(std_out,'(a,2x,3i4,2x,i4)') 'Found !:',iadum1(:),jpw
                   exit
                end if
             end do
          end do

          ikg  = ikg + npw_k

       end do    ! close loop over ikpt

       ipwnsfac = 1

    end do    ! close loop over ifor

 end do        ! close loop over idir


 call timab(1008,2,tsec)
 call timab(1009,1,tsec)

!Build mpi_enreg%kptdstrb
!array required to communicate the WFs between cpus
!(MPI // over k-points)
 if (nproc>1) then
    do idir = 1, 3
       do ifor = 1, 2

          ikpt_loc = 0
          do isppol = 1, dtset%nsppol

             do ikpt = 1, dtorbmag%fnkpt

                ikpti = dtorbmag%indkk_f2ibz(ikpt,1)
                nband_k = dtset%nband(ikpti)
                ikpt1f = dtorbmag%ikpt_dk(ikpt,ifor,idir)
                ikpt1i = dtorbmag%indkk_f2ibz(ikpt1f,1)

                if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,isppol,me)) cycle
             
                ikpt_loc = ikpt_loc + 1
                mpi_enreg%kptdstrb(me + 1,ifor+2*(idir-1),ikpt_loc) = &
&                 ikpt1i + (isppol - 1)*dtset%nkpt

                mpi_enreg%kptdstrb(me+1,ifor+2*(idir-1),ikpt_loc+dtorbmag%fmkmem_max*dtset%nsppol) = &
&                 ikpt1f + (isppol - 1)*dtorbmag%fnkpt

             end do   ! ikpt
          end do     ! isppol
       end do       ! ifor
    end do           ! idir
 end if             ! nproc>1

!build mpi_enreg%kpt_loc2fbz_sp 
 ikpt_loc = 0
 do isppol = 1, dtset%nsppol
    do ikpt = 1, dtorbmag%fnkpt

       ikpti = dtorbmag%indkk_f2ibz(ikpt,1)
       nband_k = dtset%nband(ikpti)

       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,isppol,me)) cycle
     
       ikpt_loc = ikpt_loc + 1

       mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc, 1) = ikpt
       mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc, 2) = isppol

    end do
 end do

!should be temporary
!unassigned mpi_enreg%kpt_loc2fbz_sp are empty ; inform other cpu (there are better ways...)
 mpi_enreg%mkmem(me) = dtset%mkmem
!do ii=ikpt_loc+1,dtefield%fmkmem_max
!mpi_enreg%kpt_loc2fbz_sp(me, ii, 1) = -1
!end do

 call xmpi_sum(mpi_enreg%kptdstrb,spaceComm,ierr)
 call xmpi_sum(mpi_enreg%kpt_loc2fbz_sp,spaceComm,ierr)

 ABI_DEALLOCATE(kg1_k)

 call timab(1009,2,tsec)
 call timab(1001,2,tsec)

 DBG_EXIT("COLL")

end subroutine initorbmag
!!***
