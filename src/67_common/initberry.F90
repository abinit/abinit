!{\src2tex{textfont=tt}}
!!****f* ABINIT/initberry
!! NAME
!! initberry
!!
!! FUNCTION
!! Initialization of Berryphase calculation of the polarization, the
!! ddk and the response of an insulator to a homogenous electric field.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2018 ABINIT group (MVeithen).
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
!!  mband = maximum number of bands
!!  mkmem = maximum number of k-points in core memory
!!  mpw = maximum number of plane waves
!!  natom = number of atoms in unit cell
!!  nkpt = number of k points
!!  npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!!  nsppol = 1 for unpolarized, 2 for spin-polarized
!!  nsym = number of symmetry operations
!!  ntypat = number of types of atoms in unit cell
!!  occ(mband*nkpt*nsppol) = occup number for each band at each k point
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3) = dimensional primitive vectors
!!  symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!    reciprocal space primitive translations
!!  typat = typat(natom) list of atom types
!!  usepaw = flag for PAW (1 PAW, 0 NCPP)
!!  xred(3,natom) = location of atoms in reduced units
!!
!! OUTPUT
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations
!!  pwind(pwind_alloc,2,3) = array used to compute the overlap matrix smat
!!                         between k-points k and k +- dk where dk is
!!                         parallel to the direction idir
!!    jpw = pwind(ipw,ifor,idir)
!!      * ipw = index of plane wave vector G for a given k-point k
!!      * ifor = 1: k + dk
!!               2: k - dk
!!      * idir = direction of the polarization/ddk calculation [dk(idir)
!!               is the only non-zero element of dk(:)]
!!      * jpw = index of plane wave vector G (+dG) at k +- dk
!!              where dG is a shift of one reciprocal lattice vector
!!              (required to close the strings of k-points using the
!!               periodic gauge condition)
!!    In case a G-vector of the basis sphere of plane waves at k
!!    does not belong to the basis sphere of plane waves at k+dk, jpw = 0.
!!   pwind_alloc = first dimension of pwind and pwnsfac
!!   pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!
!! SIDE EFFECTS
!!  mpi_enreg = informations about MPI parallelization
!!    kptdstrb(nproc,nneighbour,fmkmem_max*nsppol) : Array required
!!      by berryphase_new.f for MPI // over k-points. Defined
!!      for k-points in the fBZ
!!      but for k-points in the iBZ. Used by vtorho.f
!!           nproc = number of cpus
!!           nneighbour = number of neighbours for each k-point (= 6)
!!
!! TO DO
!!
!! NOTES
!!
!! PARENTS
!!      init_e_field_vars
!!
!! CHILDREN
!!      expibi,kpgsph,listkk,metric,pawcprj_alloc,pawcprj_getdim,qijb_kk
!!      setsymrhoij,smpbz,symatm,timab,wrtout,xmpi_max,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initberry(dtefield,dtset,gmet,gprimd,kg,mband,&
&              mkmem,mpi_enreg,mpw,natom,nkpt,npwarr,nsppol,&
&              nsym,ntypat,occ,pawang,pawrad,pawtab,psps,&
&              pwind,pwind_alloc,pwnsfac,&
&              rprimd,symrec,typat,usepaw,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_efield

 use m_fftcore, only : kpgsph
 use m_pawang,  only : pawang_type
 use m_pawrad,  only : pawrad_type
 use m_pawtab,  only : pawtab_type
 use m_pawcprj, only : pawcprj_alloc, pawcprj_getdim

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initberry'
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
 integer,intent(in) :: mband,mkmem,mpw,natom,nkpt,nsppol,nsym,ntypat,usepaw
 integer,intent(out) :: pwind_alloc
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(inout) :: dtset
 type(efield_type),intent(inout) :: dtefield !vz_i
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt)
 integer,intent(in) :: symrec(3,3,nsym),typat(natom)
 integer,pointer :: pwind(:,:,:)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 real(dp),pointer :: pwnsfac(:,:)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: exchn2n3d,flag,flag_kpt,fnkpt_computed,iband,icg,icprj
 integer :: idir,idum,idum1,ierr,ifor,ikg,ikg1,ikpt,ikpt1,ikpt1f
 integer :: ikpt1i,ikpt2,ikpt_loc,ikptf,ikpti,ikstr,index,ineigh,ipw,ipwnsfac
 integer :: isppol,istr,istwf_k,isym,isym1,itrs,itypat,iunmark,jpw,klmn,lmax,lmn2_size_max
 integer :: me,me_g0,mkmem_,my_nspinor,nband_k,mband_occ_k,ncpgr,nkstr,nproc,npw_k,npw_k1,spaceComm
 integer :: option, brav, mkpt, nkptlatt
 integer :: jstr,ii,jj,isign
 integer :: dk_flag, coord1, coord2
 integer :: mult
 real(dp) :: c1,ecut_eff,eg,eg_ev,rdum,diffk1,diffk2,diffk3
 real(dp) :: dist_, max_dist, last_dist, dist,kpt_shifted1,kpt_shifted2,kpt_shifted3
 real(dp) :: gprimdlc(3,3),rmetllc(3,3),gmetllc(3,3),ucvol_local
! gprimd(3,3) = inverse of rprimd
! rmetlcl(3,3)=real-space metric (same as rmet in metric.F90)
! gmetlcl(3,3)= same as gmet in metric.F90 
! ucvol = volume of the unit cell in Bohr**3

 character(len=500) :: message
 logical :: calc_epaw3_force,calc_epaw3_stress,fieldflag
!arrays
 integer :: dg(3),iadum(3),iadum1(3),neigh(6)
 integer,allocatable :: dimlmn(:),kg1_k(:,:),kpt_mark(:),nattyp_dum(:)
 real(dp) :: diffk(3),dk(3),dum33(3,3),eg_dir(3)
 real(dp) :: kpt1(3)
 real(dp) :: delta_str3(2), dstr(2),dk_str(2,2,3)
 real(dp) :: tsec(2)
 real(dp),allocatable :: calc_expibi(:,:),calc_qijb(:,:,:),spkpt(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(1001,1,tsec)
 call timab(1002,1,tsec)

!save the current value of berryopt
 dtefield%berryopt = dtset%berryopt
!save the current value of nspinor
 dtefield%nspinor = dtset%nspinor

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
   dtefield%fnkpt = fnkpt_computed
   ABI_ALLOCATE(dtefield%fkptns,(3,dtefield%fnkpt))
   dtefield%fkptns(:,:)=spkpt(:,1:dtefield%fnkpt)
   ABI_DEALLOCATE(spkpt)
 else if(dtset%kptopt==3.or.dtset%kptopt==0)then
   dtefield%fnkpt=nkpt
   ABI_ALLOCATE(dtefield%fkptns,(3,dtefield%fnkpt))
   dtefield%fkptns(1:3,1:dtefield%fnkpt)=dtset%kpt(1:3,1:dtefield%fnkpt)
   if(dtset%kptopt==0)then
     write(message,'(10a)') ch10,&
&     ' initberry : WARNING -',ch10,&
&     '  you have defined manually the k-point grid with kptopt = 0',ch10,&
&     '  the berry phase calculation works only with a regular k-points grid,',ch10,&
&     '  abinit doesn''t check if your grid is regular...'
     call wrtout(std_out,message,'PERS')
   end if
 end if

!call listkk to get mapping from FBZ to IBZ
 rdum=1.0d-5  ! cutoff distance to decide when two k points match
 ABI_ALLOCATE(dtefield%indkk_f2ibz,(dtefield%fnkpt,6))

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

!ji: The following may need modification in the future
!**** no spin-polarization doubling ; allow use of time reversal symmetry ****

!Here is original call
!
!call listkk(rdum,gmet,dtefield%indkk_f2ibz,dtset%kptns,dtefield%fkptns,nkpt,&
!& dtefield%fnkpt,dtset%nsym,1,dtset%symafm,dtset%symrel,1)

 call timab(1002,2,tsec)
 call timab(1003,1,tsec)

 call listkk(rdum,gmet,dtefield%indkk_f2ibz,dtset%kptns,dtefield%fkptns,nkpt,&
& dtefield%fnkpt,dtset%nsym,1,dtset%symafm,symrec,1,use_symrec=.True.)

 call timab(1003,2,tsec)
 call timab(1004,1,tsec)

!Construct i2fbz and f2ibz
 ABI_ALLOCATE(dtefield%i2fbz,(nkpt))
 idum=0
 do ikpt=1,dtefield%fnkpt
   if (dtefield%indkk_f2ibz(ikpt,2)==1 .and. &
&   dtefield%indkk_f2ibz(ikpt,6) == 0 .and. &
&   maxval(abs(dtefield%indkk_f2ibz(ikpt,3:5))) == 0 ) then
     dtefield%i2fbz(dtefield%indkk_f2ibz(ikpt,1))=ikpt
     idum=idum+1
   end if
 end do
 if (idum/=nkpt)then
   message = ' Found wrong number of k-points in IBZ'
   MSG_ERROR(message)
 end if

!set flags for fields, forces, stresses
 fieldflag = ( (dtset%berryopt== 4) .or. (dtset%berryopt== 6) .or. (dtset%berryopt== 7)  & 
& .or. (dtset%berryopt==14) .or. (dtset%berryopt==16) .or. (dtset%berryopt==17) )
! following two flags activates computation of projector gradient contributions to force and
! stress in finite field PAW calculations
 calc_epaw3_force = (fieldflag .and. usepaw == 1 .and. dtset%optforces /= 0)
 calc_epaw3_stress = (fieldflag .and. usepaw == 1 .and. dtset%optstress /= 0)



!----------------------------------------------------------------------------
!------------- Allocate PAW space if necessary ------------------------------
!----------------------------------------------------------------------------

 if (usepaw == 1) then

   dtefield%usepaw   = usepaw
   dtefield%natom    = natom
   dtefield%my_natom = mpi_enreg%my_natom

   ABI_ALLOCATE(dtefield%lmn_size,(ntypat))
   ABI_ALLOCATE(dtefield%lmn2_size,(ntypat))
   do itypat = 1, ntypat
     dtefield%lmn_size(itypat) = pawtab(itypat)%lmn_size
     dtefield%lmn2_size(itypat) = pawtab(itypat)%lmn2_size
   end do

   lmn2_size_max = psps%lmnmax*(psps%lmnmax+1)/2
   dtefield%lmn2max = lmn2_size_max

! expibi and qijb_kk are NOT parallelized over atoms
! this may change in the future (JZwanziger 18 March 2014)
   ABI_ALLOCATE(dtefield%qijb_kk,(2,lmn2_size_max,dtefield%natom,3))
   ABI_ALLOCATE(dtefield%expibi,(2,dtefield%natom,3))
   dtefield%has_expibi = 1
   dtefield%has_qijb = 1

   if ( fieldflag .and. dtefield%has_rij==0) then  
     lmn2_size_max = psps%lmnmax*(psps%lmnmax+1)/2
     ABI_ALLOCATE(dtefield%rij,(lmn2_size_max,ntypat,3))
     dtefield%has_rij = 1
   end if

! additional F3-type force term for finite electric field with PAW. Same term
! might also apply for other displacement-type field calculations, but not sure yet
! JZwanziger 4 April 2014
   if ( calc_epaw3_force ) then
     ABI_ALLOCATE(dtefield%epawf3,(dtefield%natom,3,3))
     dtefield%has_epawf3 = 1
   end if
   if ( calc_epaw3_stress ) then
     ABI_ALLOCATE(dtefield%epaws3,(dtefield%natom,3,6))
     dtefield%has_epaws3 = 1
   end if

   ncpgr = 0 
   if ( fieldflag .and. dtefield%usecprj == 0) then
     ABI_ALLOCATE(dimlmn,(natom))
     call pawcprj_getdim(dimlmn,natom,nattyp_dum,ntypat,typat,pawtab,'R')
!    allocate space for cprj at kpts in BZ (IBZ or FBZ)
     ABI_DATATYPE_ALLOCATE(dtefield%cprj,(natom, mband*dtset%nspinor*dtset%nkpt*nsppol))
!    write(std_out,*) "initberry alloc of cprj ", shape(dtefield%cprj)
     if (calc_epaw3_force .and. .not. calc_epaw3_stress) ncpgr = 3
     if (.not. calc_epaw3_force .and. calc_epaw3_stress) ncpgr = 6
     if (calc_epaw3_force .and. calc_epaw3_stress) ncpgr = 9
     call pawcprj_alloc(dtefield%cprj,ncpgr,dimlmn)
     dtefield%usecprj = 1
     ABI_DEALLOCATE(dimlmn)
   end if

   ABI_ALLOCATE(dtefield%cprjindex,(nkpt,nsppol))
   dtefield%cprjindex(:,:) = 0

   if (dtset%kptopt /= 3) then
     ABI_ALLOCATE(dtefield%atom_indsym,(4,nsym,natom))
     call symatm(dtefield%atom_indsym,natom,nsym,symrec,dtset%tnons,tol8,typat,xred)
     lmax = psps%mpsang - 1
     ABI_ALLOCATE(dtefield%zarot,(2*lmax+1,2*lmax+1,lmax+1,nsym))
     call setsymrhoij(gprimd,lmax,nsym,1,rprimd,symrec,dtefield%zarot)
     dtefield%nsym = nsym
     dtefield%lmax = lmax
     dtefield%lmnmax = psps%lmnmax
   end if

 end if

!------------------------------------------------------------------------------
!------------------- Compute variables related to MPI // ----------------------
!------------------------------------------------------------------------------
 spaceComm=mpi_enreg%comm_cell
 nproc=xmpi_comm_size(spaceComm)
 me=xmpi_comm_rank(spaceComm)

 if (nproc==1) then
   dtefield%fmkmem = dtefield%fnkpt
   dtefield%fmkmem_max = dtefield%fnkpt
   dtefield%mkmem_max = nkpt
 else
   dtefield%fmkmem = 0
   do ikpt = 1, dtefield%fnkpt
     ikpti = dtefield%indkk_f2ibz(ikpt,1)
     nband_k = dtset%nband(ikpti)
     if (.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,-1,me))) &
&     dtefield%fmkmem = dtefield%fmkmem + 1
   end do
!  Maximum value of mkmem and fmkmem
   call xmpi_max(dtefield%fmkmem,dtefield%fmkmem_max,spaceComm,ierr)
!  I have to use the dummy variable mkmem_ because
!  mkmem is declared as intent(in) while the first
!  argument of xmpi_max must be intent(inout)
   mkmem_ = mkmem
   call xmpi_max(mkmem_,dtefield%mkmem_max,spaceComm,ierr)
 end if

 ABI_ALLOCATE(mpi_enreg%kpt_loc2fbz_sp,(0:nproc-1,1:dtefield%fmkmem_max*nsppol, 1:2))
 ABI_ALLOCATE(mpi_enreg%kpt_loc2ibz_sp,(0:nproc-1,1:dtefield%mkmem_max*nsppol, 1:2))
 ABI_ALLOCATE(mpi_enreg%kptdstrb,(nproc,6,dtefield%fmkmem_max*nsppol*2))
 ABI_ALLOCATE(mpi_enreg%mkmem,(0:nproc-1))
 mpi_enreg%kpt_loc2fbz_sp(:,:,:) = 0
 mpi_enreg%kpt_loc2ibz_sp(:,:,:) = 0
 mpi_enreg%kptdstrb(:,:,:)       = 0
 mpi_enreg%mkmem(:)              = 0

 if (fieldflag) then
   ABI_ALLOCATE(dtefield%cgqindex,(3,6,nkpt*nsppol))
   ABI_ALLOCATE(dtefield%nneigh,(nkpt))
   dtefield%cgqindex(:,:,:) = 0 ; dtefield%nneigh(:) = 0
 end if

 pwind_alloc = mpw*dtefield%fmkmem_max
 ABI_ALLOCATE(pwind,(pwind_alloc,2,3))
 ABI_ALLOCATE(pwnsfac,(2,pwind_alloc))

!------------------------------------------------------------------------------
!---------------------- Compute efield_type variables -------------------------
!------------------------------------------------------------------------------

!Initialization of efield_type variables
 mult=dtset%useria+1
 dtefield%efield_dot(:) = zero
 dtefield%dkvecs(:,:) = zero
 dtefield%maxnstr = 0    ; dtefield%maxnkstr  = 0
 dtefield%nstr(:) = 0    ; dtefield%nkstr(:) = 0
 ABI_ALLOCATE(dtefield%ikpt_dk,(dtefield%fnkpt,2,3))
 ABI_ALLOCATE(dtefield%cgindex,(nkpt,nsppol))
 ABI_ALLOCATE(dtefield%kgindex,(nkpt))
 ABI_ALLOCATE(dtefield%fkgindex,(dtefield%fnkpt))
 dtefield%ikpt_dk(:,:,:) = 0
 dtefield%cgindex(:,:) = 0
 dtefield%mband_occ = 0
 ABI_ALLOCATE(dtefield%nband_occ,(nsppol))
 dtefield%kgindex(:) = 0
 dtefield%fkgindex(:) = 0

 if (fieldflag) then
   dtset%rfdir(1:3) = 1
 end if


!Compute spin degeneracy
 if (nsppol == 1 .and. dtset%nspinor == 1) then
   dtefield%sdeg = two
 else if (nsppol == 2 .or. my_nspinor == 2) then
   dtefield%sdeg = one
 end if

!Compute the number of occupied bands and check that
!it is the same for each k-point

 index = 0
 do isppol = 1, nsppol
   dtefield%nband_occ(isppol) = 0
   do ikpt = 1, nkpt

     mband_occ_k = 0
     nband_k = dtset%nband(ikpt + (isppol - 1)*nkpt)

     do iband = 1, nband_k
       index = index + 1
       if (abs(occ(index) - dtefield%sdeg) < tol8) mband_occ_k = mband_occ_k + 1
     end do

     if (fieldflag) then
       if (nband_k /= mband_occ_k) then
         write(message,'(a,a,a)')&
&         '  In a finite electric field, nband must be equal ',ch10,&
&         '  to the number of valence bands.'
         MSG_ERROR(message)
       end if
     end if

     if (ikpt > 1) then
       if (dtefield%nband_occ(isppol) /= mband_occ_k) then
         message = "The number of valence bands is not the same for every k-point of present spin channel"
         MSG_ERROR(message)
       end if
     else
       dtefield%mband_occ         = max(dtefield%mband_occ, mband_occ_k)
       dtefield%nband_occ(isppol) = mband_occ_k
     end if

   end do                ! close loop over ikpt
 end do                ! close loop over isppol

 if (fieldflag) then
   ABI_ALLOCATE(dtefield%smat,(2,dtefield%mband_occ,dtefield%mband_occ,nkpt*nsppol,2,3))

   dtefield%smat(:,:,:,:,:,:) = zero
 end if

 ABI_ALLOCATE(dtefield%sflag,(dtefield%mband_occ,nkpt*nsppol,2,3))
 dtefield%sflag(:,:,:,:) = 0

!Compute the location of each wavefunction

 icg = 0
 icprj = 0
!ikg = 0
 do isppol = 1, nsppol
   do ikpt = 1, nkpt

     nband_k = dtset%nband(ikpt + (isppol-1)*nkpt)

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle
     
     dtefield%cgindex(ikpt,isppol) = icg
     npw_k = npwarr(ikpt)
     icg = icg + npw_k*dtefield%nspinor*nband_k

     if (usepaw == 1) then
       dtefield%cprjindex(ikpt,isppol) = icprj
       icprj = icprj + dtefield%nspinor*nband_k
     end if

   end do
 end do

 ikg = 0
 do ikpt = 1, nkpt
   if ((proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,1,me)).and.&
&   (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,nsppol,me))) cycle
   
   npw_k = npwarr(ikpt)
   dtefield%kgindex(ikpt) = ikg
   ikg = ikg + npw_k
 end do

!Need to use dtset%red_efieldbar in the whole code
!Compute the reciprocal lattice coordinates of the electric field  
 if (fieldflag) then
   
   call  metric(gmetllc,gprimdlc,-1,rmetllc,rprimd,ucvol_local)

   if (dtset%berryopt == 4 .or. dtset%berryopt == 6 .or. dtset%berryopt == 7) then
     
     do ii=1,3
       dtset%red_efieldbar(ii) = dot_product(dtset%efield(:),rprimd(:,ii))
       dtefield%efield_dot(ii) =  dtset%red_efieldbar(ii)
     end do

!    dtefield%efield_dot(1) = dot_product(dtset%efield(:),rprimd(:,1))
!    dtefield%efield_dot(2) = dot_product(dtset%efield(:),rprimd(:,2))
!    dtefield%efield_dot(3) = dot_product(dtset%efield(:),rprimd(:,3))

     write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
&     ' initberry: Reduced electric field (ebar)',ch10,&
&     '  red_efieldbar(1:3) = ',dtset%red_efieldbar(1:3),ch10
     call wrtout(std_out,message,'COLL')

   end if

   if (dtset%berryopt == 6 .or. dtset%berryopt ==7 ) then

     do ii=1,3 
       dtset%red_dfield(ii)= (dot_product(dtset%dfield(:),gprimdlc(:,ii)))*ucvol_local/(4.d0*pi)
     end do

     write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
&     ' initberry: Reduced electric displacement field',ch10,&
&     '  red_dfield(1:3) = ',dtset%red_dfield(1:3),ch10
     call wrtout(std_out,message,'COLL')

   end if


   if (  dtset%berryopt == 14 ) then
!    transfer to unreduced electric field.
     do idir=1,3
       dtset%efield(idir)= dot_product(dtset%red_efieldbar(:),gprimdlc(:,idir))
       dtefield%efield_dot(idir) = dtset%red_efieldbar(idir) 
!      dtefield%efield2(idir)=dtset%red_efieldbar(idir) 
     end do

!    dtefield%efield_dot(1) = dtset%red_efieldbar(1) 
!    dtefield%efield_dot(2) = dtset%red_efieldbar(2)
!    dtefield%efield_dot(3) = dtset%red_efieldbar(3)

     write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
&     ' initberry: Unreduced electric field (a.u.)',ch10,&
&     '  efield(1:3) = ',dtset%efield(1:3),ch10
     call wrtout(std_out,message,'COLL')

     write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
&     ' initberry: Reduced electric field (ebar)',ch10,&
&     '  red_efieldbar(1:3) = ',dtset%red_efieldbar(1:3),ch10
     call wrtout(std_out,message,'COLL')

   end if


   if ( dtset%berryopt == 16 ) then

!    to calculate D
     do ii=1,3
       dtset%dfield(ii)  =(4*pi/ucvol_local)*dot_product(dtset%red_dfield(:),rprimd(:,ii))
     end do

     do idir=1,3
       dtset%efield(idir)= (4*pi/ucvol_local)*dot_product(dtset%red_efield(:),rprimd(:,idir))
     end do

     do idir=1,3
       dtset%red_efieldbar(idir)= (4*pi/ucvol_local)*dot_product(dtset%red_efield(:),rmetllc(:,idir))
       dtefield%efield_dot(idir) = dtset%red_efieldbar(idir) 
     end do
     
!    dtefield%efield_dot(1) = dtset%red_efieldbar(1) 
!    dtefield%efield_dot(2) = dtset%red_efieldbar(2)
!    dtefield%efield_dot(3) = dtset%red_efieldbar(3)


     write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
&     ' initberry: Unreduced electric displacement field (a.u.)',ch10,&
&     '  dfield(1:3) = ',dtset%dfield(1:3),ch10
     call wrtout(std_out,message,'COLL')

     write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
&     ' initberry: Unreduced electric field (a.u.)',ch10,&
&     '  efield(1:3) = ',dtset%efield(1:3),ch10
     call wrtout(std_out,message,'COLL')

     write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
&     ' initberry: Reduced electric field (ebar)',ch10,&
&     '  red_efieldbar(1:3) = ',dtset%red_efieldbar(1:3),ch10
     call wrtout(std_out,message,'COLL')

   end if

   if ( dtset%berryopt ==17) then

!    to calculate D
     
     do idir=1,3
       dtset%efield(idir)= dot_product(dtset%red_efieldbar(:),gprimdlc(:,idir))  ! from ebar 
       dtset%dfield(idir)  =(4*pi/ucvol_local)*dot_product(dtset%red_dfield(:),rprimd(:,idir))
!      dtset%red_efield(idir) = (ucvol_local/(4*pi))*dot_product(dtset%red_efieldbar(:),gmetllc(:,idir))
       dtefield%efield_dot(idir) = dtset%red_efieldbar(idir) 
     end do

     write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
&     ' initberry: Reduced electric field (ebar)',ch10,&
&     '  red_efieldbar(1:3) = ',dtset%red_efieldbar(1:3),ch10
     call wrtout(std_out,message,'COLL')


     write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
&     ' initberry: Unreduced electric field (a.u.)',ch10,&
&     '  efield(1:3) = ',dtset%efield(1:3),ch10
     call wrtout(std_out,message,'COLL')

     write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
&     ' initberry: Reduced electric displacement field (a.u.)',ch10,&
&     '  red_dfield(1:3) = ',dtset%red_dfield(1:3),ch10
     call wrtout(std_out,message,'COLL')

     write(message,'(a,a,a,a,3(2x,f16.9),a)')ch10,&
&     ' initberry: Unreduced electric displacement field (a.u.)',ch10,&
&     '  dfield(1:3) = ',dtset%dfield(1:3),ch10
     call wrtout(std_out,message,'COLL')


   end if



 end if

 call timab(1004,2,tsec)

!------------------------------------------------------------------------------
!---------------------- Compute dk --------------------------------------------
!------------------------------------------------------------------------------

 call timab(1005,1,tsec)

 do idir = 1, 3

   if (dtset%rfdir(idir) == 1) then

!    Compute dk(:), the vector between a k-point and its nearest
!    neighbour along the direction idir

     dk(:) = zero
     dk(idir) = 1._dp   ! 1 mean there is no other k-point un the direction idir
     do ikpt = 2, dtefield%fnkpt
       diffk(:) = abs(dtefield%fkptns(:,ikpt) - dtefield%fkptns(:,1))
       if ((diffk(1) < dk(1)+tol8).and.(diffk(2) < dk(2)+tol8).and.&
&       (diffk(3) < dk(3)+tol8)) dk(:) = diffk(:)
     end do
     dtefield%dkvecs(:,idir) = dk(:)
!    DEBUG
!    write(std_out,*)' initberry : idir, dk', idir, dk
!    ENDDEBUG

!    For each k point, find k_prim such that k_prim= k + dk mod(G)
!    where G is a vector of the reciprocal lattice

     do ikpt = 1, dtefield%fnkpt

!      First k+dk, then k-dk
       do isign=-1,1,2
         kpt_shifted1=dtefield%fkptns(1,ikpt)- isign*dk(1)
         kpt_shifted2=dtefield%fkptns(2,ikpt)- isign*dk(2)
         kpt_shifted3=dtefield%fkptns(3,ikpt)- isign*dk(3)
!        Note that this is still a order fnkpt**2 algorithm.
!        It is possible to implement a order fnkpt algorithm, see listkk.F90.
         do ikpt1 = 1, dtefield%fnkpt
           diffk1=dtefield%fkptns(1,ikpt1) - kpt_shifted1
           if(abs(diffk1-nint(diffk1))>tol8)cycle
           diffk2=dtefield%fkptns(2,ikpt1) - kpt_shifted2
           if(abs(diffk2-nint(diffk2))>tol8)cycle
           diffk3=dtefield%fkptns(3,ikpt1) - kpt_shifted3
           if(abs(diffk3-nint(diffk3))>tol8)cycle
           dtefield%ikpt_dk(ikpt,(isign+3)/2,idir) = ikpt1
           exit
         end do   ! ikpt1
       end do     ! isign

!      OLD CODING
!      First: k + dk 
!      do ikpt1 = 1, dtefield%fnkpt
!      diffk(:) = abs(dtefield%fkptns(:,ikpt1) - &
!      &         dtefield%fkptns(:,ikpt) - dk(:))
!      if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
!      dtefield%ikpt_dk(ikpt,1,idir) = ikpt1
!      exit
!      end if
!      end do
       
!      Second: k - dk
!      do ikpt1 = 1, dtefield%fnkpt
!      diffk(:) = abs(dtefield%fkptns(:,ikpt1) - &
!      &         dtefield%fkptns(:,ikpt) + dk(:))
!      if(sum(abs(diffk(:) - nint(diffk(:)))) < 3*tol8) then
!      dtefield%ikpt_dk(ikpt,2,idir) = ikpt1
!      exit
!      end if
!      end do 

     end do     ! ikpt

!    Find the string length, starting from k point 1
!    (all strings must have the same number of points)

     nkstr = 1
     ikpt1 = 1
     do ikpt = 1, dtefield%fnkpt
       ikpt1 = dtefield%ikpt_dk(ikpt1,1,idir)
       if (ikpt1 == 1) exit
       nkstr = nkstr + 1
     end do

!    Check that the string length is a divisor of nkpt
     if(mod(dtefield%fnkpt,nkstr) /= 0) then
       write(message,'(a,i5,a,i7)')&
&       ' The string length = ',nkstr,&
&       ', is not a divisor of fnkpt =',dtefield%fnkpt
       MSG_BUG(message)
     end if

     dtefield%nkstr(idir) = nkstr
     dtefield%nstr(idir)  = dtefield%fnkpt/nkstr

   end if      ! dtset%rfdir(idir) == 1

   write(message,'(a,i1,a,i3,a,i6)')&
&   '  initberry: for direction ',idir,', nkstr = ',dtefield%nkstr(idir),&
&   ', nstr = ',dtefield%nstr(idir)
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

 end do     ! close loop over idir

 call timab(1005,2,tsec)
 call timab(1006,1,tsec)

 dtefield%maxnstr  = maxval(dtefield%nstr(:))
 dtefield%maxnkstr = maxval(dtefield%nkstr(:))
 ABI_ALLOCATE(dtefield%idxkstr,(dtefield%maxnkstr,dtefield%maxnstr,3))
 dtefield%idxkstr(:,:,:) = 0

!for the geometry of the string space :
 ABI_ALLOCATE(dtefield%coord_str,(2,dtefield%maxnstr,3))
 ABI_ALLOCATE(dtefield%str_neigh,(-2:2,dtefield%maxnstr,3))
 ABI_ALLOCATE(dtefield%strg_neigh,(-2:2,dtefield%maxnstr,2,3))
 dtefield%coord_str(:,:,:) = 0.d0
 dtefield%str_neigh(:,:,:)=0
 dtefield%strg_neigh(:,:,:,:)=0
 dtefield%gmet_str(:,:,:)=0.d0

!------------------------------------------------------------------------------
!---------------------- Build the strings -------------------------------------
!------------------------------------------------------------------------------

 ABI_ALLOCATE(kpt_mark,(dtefield%fnkpt))
 do idir = 1, 3

   if (dtset%rfdir(idir) == 1) then

     iunmark = 1
     kpt_mark(:) = 0
     do istr = 1, dtefield%nstr(idir)

       do while(kpt_mark(iunmark) /= 0)
         iunmark = iunmark + 1
       end do
       dtefield%idxkstr(1,istr,idir) = iunmark
       kpt_mark(iunmark) = 1
       do ikstr = 2, dtefield%nkstr(idir)
         ikpt1 = dtefield%idxkstr(ikstr-1,istr,idir)
         ikpt2 = dtefield%ikpt_dk(ikpt1,1,idir)
         dtefield%idxkstr(ikstr,istr,idir) = ikpt2
         kpt_mark(ikpt2) = 1
       end do

     end do    ! istr

!    compute distance between strings
!    compute the metric matrix of the strings space in the direction idir
     do ii = 1,3
       do jj = 1,3
         if (ii<idir.and.jj<idir) dtefield%gmet_str(ii  ,jj  ,idir) = &
&         gmet(ii,jj) - gmet(ii,idir)*gmet(jj,idir)/gmet(idir,idir)
         if (ii<idir.and.jj>idir) dtefield%gmet_str(ii  ,jj-1,idir) = &
&         gmet(ii,jj) - gmet(ii,idir)*gmet(jj,idir)/gmet(idir,idir)
         if (ii>idir.and.jj<idir) dtefield%gmet_str(ii-1,jj  ,idir) = &
&         gmet(ii,jj) - gmet(ii,idir)*gmet(jj,idir)/gmet(idir,idir)
         if (ii>idir.and.jj>idir) dtefield%gmet_str(ii-1,jj-1,idir) = &
&         gmet(ii,jj) - gmet(ii,idir)*gmet(jj,idir)/gmet(idir,idir)
       end do
     end do
!    DEBUG
!    write(std_out,*)'gmet'
!    do ii=1,3
!    write(std_out,*)gmet(ii,:)
!    end do
!    write(std_out,*)'gmet_str'
!    do ii=1,2
!    write(std_out,*)dtefield%gmet_str(ii,:,idir)
!    end do
!    ENDDEBUG
     do istr = 1, dtefield%nstr(idir)
       do ii = 1,3
         if (ii<idir) dtefield%coord_str(ii,istr,idir)=dtefield%fkptns(ii,dtefield%idxkstr(1,istr,idir))
         if (ii>idir) dtefield%coord_str(ii-1,istr,idir)=dtefield%fkptns(ii,dtefield%idxkstr(1,istr,idir))
       end do
     end do

!    the following is very similar to getshell
     dist_ = 0._dp
     do ii = 1,2
       dist_ = dist_ + dtefield%gmet_str(ii,ii,idir)
     end do
     max_dist = 2._dp * dist_ * 2._dp

     dk_str(:,:,idir) = 0._dp
     last_dist = 0._dp
!    ishell = 0
!    dtefield%str_neigh(:,:,:) = 0
     dk_flag = 0
     do while (dk_flag /= 2)
!      Advance shell counter
!      ishell = ishell + 1

!      Search the smallest distance between two strings
       dist = max_dist
       do istr = 1,dtefield%nstr(idir)
         delta_str3(:) = dtefield%coord_str(:,1,idir) - dtefield%coord_str(:,istr,idir)
         do coord1 = -1,1  !two loop to search also on the border of the BZ
           do coord2 = -1,1
             dist_ = 0._dp
             dstr(:) = delta_str3(:) - nint(delta_str3(:))
             dstr(1) = dstr(1) + real(coord1,dp)
             dstr(2) = dstr(2) + real(coord2,dp)
             do ii = 1,2
               do jj = 1,2
                 dist_ = dist_ + dstr(ii)*dtefield%gmet_str(ii,jj,idir)*dstr(jj)
               end do
             end do
             if ((dist_ < dist).and.(dist_ - last_dist > tol8)) then
               dist = dist_
             end if
           end do
         end do
       end do

       last_dist = dist

!      search the connecting vectors for that distance
       do istr = 1,dtefield%nstr(idir)
         delta_str3(:) = dtefield%coord_str(:,istr,idir) - dtefield%coord_str(:,1,idir)
         do coord1 = -1,1
           do coord2 = -1,1
             dist_ = 0._dp
             dstr(:) = delta_str3(:) - nint(delta_str3(:))
             dstr(1) = dstr(1) + real(coord1,dp)
             dstr(2) = dstr(2) + real(coord2,dp)
             do ii = 1,2
               do jj = 1,2
                 dist_ = dist_ + dstr(ii)*dtefield%gmet_str(ii,jj,idir)*dstr(jj)
               end do
             end do
             if (abs(dist_ - dist) < tol8) then
               if (dk_flag == 0) then
                 dk_str(:,1,idir) = dstr(:)
                 dk_flag = 1
!                DEBUG
!                write(std_out,'(a,i4,2e15.4)')'1st connect', istr, dstr
!                ENDDEBUG
               elseif (dk_str(1,1,idir)*dstr(2)-dk_str(2,1,idir)*dstr(1) > tol8) then
                 dk_str(:,2,idir) = dstr(:)
                 dk_flag = 2
!                DEBUG
!                write(std_out,'(a,i4,2e15.4)')'2nd connect', istr, dstr
!                ENDDEBUG
                 exit
               end if
             end if
           end do
           if (dk_flag == 2) exit
         end do
         if (dk_flag == 2) exit
       end do

     end do ! do while

!    search the two neighbours for each string
     do istr = 1,dtefield%nstr(idir)
       dtefield%str_neigh(0,istr,idir) = istr
       dtefield%strg_neigh(0,istr,:,idir) = 0
       do jstr = 1,dtefield%nstr(idir)
         delta_str3(:) = dtefield%coord_str(:,jstr,idir) - dtefield%coord_str(:,istr,idir)
         do coord1 = -1,1
           do coord2 = -1,1
             dist_ = 0._dp
             dstr(:) = delta_str3(:) - nint(delta_str3(:))
             dstr(1) = dstr(1) + real(coord1,dp)
             dstr(2) = dstr(2) + real(coord2,dp)
             do ii = 1,2
               if (sum(abs(dstr(:)-dk_str(:,ii,idir)))<tol8) then
                 dtefield%str_neigh(ii,istr,idir) = jstr
                 dtefield%strg_neigh(ii,istr,1,idir) = coord1
                 dtefield%strg_neigh(ii,istr,2,idir) = coord2
               elseif (sum(abs(dstr(:)+dk_str(:,ii,idir)))<tol8) then
                 dtefield%str_neigh(-ii,istr,idir) = jstr
                 dtefield%strg_neigh(-ii,istr,1,idir) = coord1
                 dtefield%strg_neigh(-ii,istr,2,idir) = coord2
               end if
             end do
           end do
         end do
       end do
     end do

!    DEBUG
!    write(std_out,'(a,e15.4,e15.4,e15.4,e15.4)')'dk_str',dk_str(1,1,idir),dk_str(2,1,idir),dk_str(1,2,idir),dk_str(2,2,idir)
!    write(std_out,*)'istr, neigh1, strg(1,:), neigh2, strg(2,:),neigh-1, strg(-1,:), neigh-2, strg(-2,:)'
!    do istr=1,dtefield%nstr(idir)
!    write(std_out,'(13i4)')istr, &
!    &       dtefield%str_neigh(1,istr,idir), dtefield%strg_neigh(1,istr,:,idir),&
!    &       dtefield%str_neigh(2,istr,idir), dtefield%strg_neigh(2,istr,:,idir),&
!    &       dtefield%str_neigh(-1,istr,idir), dtefield%strg_neigh(-1,istr,:,idir),&
!    &       dtefield%str_neigh(-2,istr,idir), dtefield%strg_neigh(-2,istr,:,idir)
!    end do
!    ENDDEBUG


   end if         ! rfdir(idir) == 1

 end do           ! close loop over idir

 ABI_DEALLOCATE(kpt_mark)

 call timab(1006,2,tsec)
 call timab(1007,1,tsec)

!------------------------------------------------------------------------------
!------------ Compute PAW on-site terms if necessary --------------------------
!------------------------------------------------------------------------------

 if (usepaw == 1 .and. dtefield%has_expibi == 1) then
    ABI_ALLOCATE(calc_expibi,(2,natom))
    do idir = 1, 3
       dk = dtefield%dkvecs(1:3,idir)
       calc_expibi = zero
       call expibi(calc_expibi,dk,natom,xred)
       dtefield%expibi(1:2,1:natom,idir) = calc_expibi
    end do
!   call expibi(dtefield%expibi,dtefield%dkvecs,natom,xred)
   dtefield%has_expibi = 2
   ABI_DEALLOCATE(calc_expibi)
 end if

 if (usepaw == 1 .and. dtefield%has_qijb == 1) then
    ABI_ALLOCATE(calc_qijb,(2,dtefield%lmn2max,natom))

    do idir = 1, 3
       dk = dtefield%dkvecs(1:3,idir)
       calc_qijb = zero
       call qijb_kk(calc_qijb,dk,dtefield%expibi(1:2,1:natom,idir),&
&        gprimd,dtefield%lmn2max,natom,ntypat,pawang,pawrad,pawtab,typat)
       dtefield%qijb_kk(1:2,1:dtefield%lmn2max,1:natom,idir) = calc_qijb
!    call qijb_kk(dtefield%qijb_kk,dtefield%dkvecs,dtefield%expibi,&
! &   gprimd,dtefield%lmn2max,natom,ntypat,pawang,pawrad,pawtab,typat)
    end do
    dtefield%has_qijb = 2
    ABI_DEALLOCATE(calc_qijb)
 end if
 
 if (usepaw == 1 .and. dtefield%has_rij == 1) then
   c1=sqrt(four_pi/three)
   do itypat = 1, ntypat
     do klmn = 1, pawtab(itypat)%lmn2_size
       dtefield%rij(klmn,itypat,1) = c1*pawtab(itypat)%qijl(4,klmn) ! S_{1,1} ~ x
       dtefield%rij(klmn,itypat,2) = c1*pawtab(itypat)%qijl(2,klmn) ! S_{1,-1} ~ y
       dtefield%rij(klmn,itypat,3) = c1*pawtab(itypat)%qijl(3,klmn) ! S_{1,0} ~ z
     end do ! end loop over klmn
   end do ! end loop over itypat
   dtefield%has_rij = 2
 end if !

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
 ABI_ALLOCATE(kg1_k,(3,mpw))

 ipwnsfac = 0

 do idir = 1, 3

   if (dtset%rfdir(idir) == 1) then

     dk(:) = dtefield%dkvecs(:,idir)

     do ifor = 1, 2

       if (ifor == 2) dk(:) = -1._dp*dk(:)

!      Build pwind and kgindex
!      NOTE: The array kgindex is important for parallel execution.
!      In case nsppol = 2, it may happent that a particular processor
!      treats k-points at different spin polarizations.
!      In this case, it is not possible to address the elements of
!      pwind correctly without making use of the kgindex array.

       ikg = 0 ; ikpt_loc = 0 ; isppol = 1
       do ikpt = 1, dtefield%fnkpt

         ikpti = dtefield%indkk_f2ibz(ikpt,1)
         nband_k = dtset%nband(ikpti)
         ikpt1f = dtefield%ikpt_dk(ikpt,ifor,idir)
         ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

         if ((proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,1,me)).and.&
&         (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,nsppol,me))) cycle

         ikpt_loc = ikpt_loc + 1

!        Build basis sphere of plane waves for the nearest neighbour of
!        the k-point (important for MPI //)

         kg1_k(:,:) = 0
         kpt1(:) = dtset%kptns(:,ikpt1i)
         call kpgsph(ecut_eff,exchn2n3d,gmet,ikg1,ikpt,istwf_k,kg1_k,kpt1,&
&         1,mpi_enreg,mpw,npw_k1)
         me_g0=mpi_enreg%me_g0


!        ji: fkgindex is defined here !
         dtefield%fkgindex(ikpt) = ikg

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

         isym  = dtefield%indkk_f2ibz(ikpt,2)
         isym1 = dtefield%indkk_f2ibz(ikpt1f,2)

!        Construct transformed G vector that enters the matching condition:
!        alpha(k) S(k)^{t,-1} ( -G(b) - GBZ(k) + G(k) )

         dg(:) = -dtefield%indkk_f2ibz(ikpt,3:5) &
&         -nint(-dtefield%fkptns(:,ikpt) - dk(:) - tol10 + &
&         dtefield%fkptns(:,ikpt1f)) &
&         +dtefield%indkk_f2ibz(ikpt1f,3:5)

!        old code
!        iadum(:)=0
!        do idum=1,3
!        iadum(:)=iadum(:)+ symrec(:,idum,isym1)*dg(idum)
!        end do

!        new code
         iadum(:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),dg(:))

         dg(:) = iadum(:)

         if ( dtefield%indkk_f2ibz(ikpt1f,6) == 1 ) dg(:) = -dg(:)

!        Construct S(k)^{t,-1} S(b)^{t}

         dum33(:,:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),symrec(:,:,isym))

!        Construct alpha(k) alpha(b)

         if (dtefield%indkk_f2ibz(ikpt,6) == dtefield%indkk_f2ibz(ikpt1f,6)) then
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
           iadum(:) = kg(:,dtefield%kgindex(ikpti) + ipw)

!          Store non-symmorphic operation phase factor exp[i2\pi \alpha G \cdot t]

           if ( ipwnsfac == 0 ) then
!            old code
             rdum=0.0_dp
             do idum=1,3
               rdum=rdum+dble(iadum(idum))*dtset%tnons(idum,isym)
             end do
             rdum=two_pi*rdum
             if ( dtefield%indkk_f2ibz(ikpt,6) == 1 ) rdum=-rdum
             pwnsfac(1,ikg+ipw) = cos(rdum)
             pwnsfac(2,ikg+ipw) = sin(rdum)
!            
!            new code
!            rdum = DOT_PRODUCT(dble(iadum(:)),dtset%tnons(:,isym))
!            rdum= two_pi*rdum
!            if ( dtefield%indkk_f2ibz(ikpt,6) == 1 ) rdum=-rdum
!            pwnsfac(1,ikg+ipw) = cos(rdum)
!            pwnsfac(2,ikg+ipw) = sin(rdum)

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
&             (iadum(2) == iadum1(2)).and. &
&             (iadum(3) == iadum1(3)) ) then
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

   end if      ! rfdir(idir) == 1

 end do        ! close loop over idir


 call timab(1008,2,tsec)
 call timab(1009,1,tsec)

!Build mpi_enreg%kptdstrb
!array required to communicate the WFs between cpus in berryphase_new.f
!(MPI // over k-points)
 if (nproc>1) then
   do idir = 1, 3
     if (dtset%rfdir(idir) == 1) then
       do ifor = 1, 2

         ikpt_loc = 0
         do isppol = 1, nsppol

           do ikpt = 1, dtefield%fnkpt

             ikpti = dtefield%indkk_f2ibz(ikpt,1)
             nband_k = dtset%nband(ikpti)
             ikpt1f = dtefield%ikpt_dk(ikpt,ifor,idir)
             ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

             if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,isppol,me)) cycle
             
             ikpt_loc = ikpt_loc + 1
             mpi_enreg%kptdstrb(me + 1,ifor+2*(idir-1),ikpt_loc) = &
&             ikpt1i + (isppol - 1)*nkpt

             mpi_enreg%kptdstrb(me+1,ifor+2*(idir-1),&
&             ikpt_loc+dtefield%fmkmem_max*nsppol) = &
&             ikpt1f + (isppol - 1)*dtefield%fnkpt

           end do   ! ikpt
         end do     ! isppol
       end do       ! ifor
     end if         ! dtset%rfdir(idir) == 1
   end do           ! idir
 end if             ! nproc>1

!build mpi_enreg%kpt_loc2fbz_sp 
 ikpt_loc = 0
 do isppol = 1, nsppol
   do ikpt = 1, dtefield%fnkpt

     ikpti = dtefield%indkk_f2ibz(ikpt,1)
     nband_k = dtset%nband(ikpti)

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,isppol,me)) cycle
     
     ikpt_loc = ikpt_loc + 1

     mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc, 1) = ikpt
     mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc, 2) = isppol

   end do
 end do


!parallel case only :
!build mpi_enreg%kpt_loc2ibz_sp, dtefield%cgqindex and dtefield%nneigh
 if ((fieldflag).and.(nproc>1)) then
   ikpt_loc = 0
   do isppol = 1, nsppol
     do ikpt = 1, nkpt

       ikptf = dtefield%i2fbz(ikpt)
       nband_k = dtset%nband(ikpti)

       neigh(:) = 0 ; icg = 0 ; ikg = 0 ; flag_kpt = 0; icprj = 0
       do idir=1, 3

!        skip idir values for which efield_dot(idir) = 0
         if (abs(dtefield%efield_dot(idir)) < tol12) cycle

         do ifor = 1, 2

           flag = 0

           ikpt1f = dtefield%ikpt_dk(ikptf,ifor,idir)
           ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

           dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ikg
           ikg = ikg + npwarr(ikpt1i)

!          check if this neighbour is also a previous neighbour
           do ineigh = 1, (ifor+2*(idir-1))
             if (neigh(ineigh) == ikpt1i) then
               flag = 1
               dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ineigh
               dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
&               dtefield%cgqindex(2,ineigh,ikpt+(isppol-1)*nkpt)
               exit
             end if
           end do
!          create the cgqindex of the neighbour if necessary
           if (flag == 0) then
             neigh(ifor+2*(idir-1)) = ikpt1i
             dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
&             ifor+2*(idir-1)
             dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = icg
             if (isppol == 1) dtefield%nneigh(ikpt) = dtefield%nneigh(ikpt) + 1
             icg = icg + npwarr(ikpt1i)*dtefield%nspinor*nband_k
           end if
         end do !ifor
       end do !idir

       if (.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me))) then
!        ikpt is one of my kpt_loc
         ikpt_loc = ikpt_loc + 1
         mpi_enreg%kpt_loc2ibz_sp(me, ikpt_loc, 1) = ikpt
         mpi_enreg%kpt_loc2ibz_sp(me, ikpt_loc, 2) = isppol
       end if

     end do !ikpt
   end do !isppol
 end if !nproc>1

!should be temporary
!unassigned mpi_enreg%kpt_loc2fbz_sp are empty ; inform other cpu (there are better ways...)
 mpi_enreg%mkmem(me) = mkmem
!do ii=ikpt_loc+1,dtefield%fmkmem_max
!mpi_enreg%kpt_loc2fbz_sp(me, ii, 1) = -1
!end do


!(same as mpi_enreg%kptdstrb but for k-points in the iBZ),
!dtefield%cgqindex and dtefield%nneigh

 if ((fieldflag).and.(nproc>1)) then

   ikpt_loc = 1
   do isppol = 1, nsppol
     do ikpt = 1, nkpt

       nband_k = dtset%nband(ikpt)
       ikptf = dtefield%i2fbz(ikpt)

       neigh(:) = 0 ; icg = 0 ; ikg = 0 ; flag_kpt = 0; icprj = 0
       do idir = 1, 3

!        Skip idir values for which efield_dot(idir) = 0
         if (abs(dtefield%efield_dot(idir)) < tol12 .and. (fieldflag)) cycle

         do ifor = 1, 2

           ikpt1f = dtefield%ikpt_dk(ikptf,ifor,idir)
           ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)

!          dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ikg
           ikg = ikg + npwarr(ikpt1i)

           flag = 0
           do ineigh = 1, (ifor+2*(idir-1))
             if (neigh(ineigh) == ikpt1i) then
               flag = 1
!              dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = ineigh
!              dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
!              &               dtefield%cgqindex(2,ineigh,ikpt+(isppol-1)*nkpt)
               exit
             end if
           end do
           if (flag == 0) then
!            neigh(ifor+2*(idir-1)) = ikpt1i
!            dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = &
!            &             ifor+2*(idir-1)
!            dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt) = icg
!            if (isppol == 1) dtefield%nneigh(ikpt) = dtefield%nneigh(ikpt) + 1
!            icg = icg + npwarr(ikpt1i)*dtset%nspinor*nband_k
           end if

           if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle

           flag_kpt = 1

!          MVeithen: the if condition allows to avoid that the same wavefunction
!          is send several times to a particular cpu

         end do    ! ifor
       end do    ! idir

       if (flag_kpt == 1) ikpt_loc = ikpt_loc + 1

     end do    ! ikpt
   end do    ! isppol

 end if   ! fieldflag

 call xmpi_sum(mpi_enreg%kptdstrb,spaceComm,ierr)
 call xmpi_sum(mpi_enreg%kpt_loc2fbz_sp,spaceComm,ierr)
 if (fieldflag) then
   call xmpi_sum(mpi_enreg%kpt_loc2ibz_sp,spaceComm,ierr)
   call xmpi_sum(mpi_enreg%mkmem,spaceComm,ierr)
 end if

!------------------------------------------------------------------------------
!------------------------ Estimate critical field -----------------------------
!------------------------------------------------------------------------------

!Compute the minimal value of the bandgap required to be below
!the critical field as defined by the relation
!| E_i*a_i | < E_g/n_i

 if (fieldflag) then

   do idir = 1, 3
!    eg_dir(idir) = abs(dtefield%efield_dot(idir))*dtefield%nkstr(idir)
     eg_dir(idir) = abs(dtset%red_efieldbar(idir))*dtefield%nkstr(idir)
   end do

   
   eg = maxval(eg_dir)
   eg_ev = eg*Ha_eV

   if (dtset%optcell ==0 .and. (dtset%berryopt == 4 .or. dtset%berryopt == 14)) then
     write(message,'(a,a,a,a,a,a,a,a,f7.2,a,a)')ch10,&
&     ' initberry: COMMENT - ',ch10,&
&     '  As a rough estimate,',ch10,&
&     '  to be below the critical field, the bandgap of your system',ch10,&
&     '  should be larger than ',eg_ev,' eV.',ch10
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

   else

     write(message,'(a,a,a,a,a,a,a)') ch10,&
&     ' initberry: COMMENT - ',ch10,&
&     '  The estimation of critical electric field should be checked after calculation.',ch10,&
&     '  It is printed out just after total energy.' ,ch10

     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

   end if 

 end if

 ABI_DEALLOCATE(kg1_k)

 call timab(1009,2,tsec)
 call timab(1001,2,tsec)

 DBG_EXIT("COLL")

end subroutine initberry
!!***
