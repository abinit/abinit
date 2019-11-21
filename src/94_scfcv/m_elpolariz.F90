!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_elpolariz
!! NAME
!! m_elpolariz
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2019 ABINIT group (XG, NSAI, MKV)
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

module m_elpolariz

 use defs_basis
 use m_abicore
 use m_errors
 use m_efield
 use m_xmpi
 use m_wffile
 use m_hdr
 use m_dtset
 use m_dtfil

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes, only : MPI_type
 use m_geometry, only : metric
 use m_symtk,    only : matr3inv
 use m_hide_lapack,  only : dzgedi, dzgefa
 use m_rwwf,     only : rwwf
 use m_pawtab,   only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type
 use m_pawcprj,  only : pawcprj_type
 use m_berryphase, only : berryphase
 use m_berryphase_new, only : berryphase_new

 implicit none

 private
!!***

 public :: elpolariz
!!***

contains
!!***

!!****f* ABINIT/elpolariz
!! NAME
!! elpolariz
!!
!! FUNCTION
!! Calculate corrections to total energy from polarising
!! electric field with or without Berry phases (berryopt keyword)
!!
!! INPUTS
!! atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!! cg(2,mcg)=planewave coefficients of wavefunctions
!! cprj(natom,mcprj*usecprj)=<p_lmn|Cnk> coefficients for each WF |Cnk>
!!                           and each |p_lmn> non-local projector
!! dtfil <type(datafiles_type)>=variables related to files
!! dtset <type(dataset_type)>=all input variables in this dataset
!! gprimd(3,3)=reciprocal space dimensional primitive translations
!! hdr <type(hdr_type)>=the header of wf, den and pot files
!! kg(3,mpw*mkmem)=reduced planewave coordinates
!! mband=maximum number of bands
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mkmem=number of k points treated by this node.
!! mpi_enreg=information about MPI parallelization
!! mpw=maximum dimensioned size of npw
!! my_natom=number of atoms treated by current processor
!! natom=number of atoms in cell
!! nattyp(ntypat)= # atoms of each type.
!! nkpt=number of k points
!! npwarr(nkpt)=number of planewaves in basis at this k point
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! ntypat=number of types of atoms in unit cell
!! nkpt=number of k-points
!! option = 1: compute Berryphase polarization
!!          2: compute finite difference expression of the ddk
!!          3: compute polarization & ddk
!! pawrhoij(my_natom*usepaw) <type(pawrhoij_type)> atomic occupancies
!! pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! pel_cg(3) = reduced coordinates of the electronic polarization (a. u.)
!!             computed in the SCF loop
!! pelev(3)= expectation value polarization term (PAW only) in cartesian coordinates
!! psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!! pwind_alloc = first dimension of pwind
!! pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! ucvol=unit cell volume in bohr**3.
!! usecprj=1 if cprj datastructure has been allocated
!! xred(3,natom)=reduced atomic coordinates
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! dtefield <type(efield_type)> = variables related to Berry phase
!!       and electric field calculations (see initberry.f).
!!       In case berryopt = 4/6/7/14/16/17, the overlap matrices computed
!!       in this routine are stored in dtefield%smat in order
!!       to be used in the electric field calculation.
!! enefield=field energy
!! etotal=total energy, might be correct by improved polarization computation
!! pel(3) = reduced coordinates of the electronic polarization (a. u.)
!! pion(3)= reduced coordinates of the ionic polarization (a. u.)
!!
!! PARENTS
!!      afterscfloop
!!
!! CHILDREN
!!      berryphase,berryphase_new,metric,uderiv,wrtout
!!
!! SOURCE

subroutine elpolariz(atindx1,cg,cprj,dtefield,dtfil,dtset,etotal,enefield,gprimd,hdr,&
& kg,mband,mcg,mcprj,mkmem,mpi_enreg,mpw,my_natom,natom,nattyp,nkpt,&
& npwarr,nsppol,ntypat,pawrhoij,pawtab,&
& pel,pel_cg,pelev,pion,psps,pwind,pwind_alloc,&
& pwnsfac,rprimd,ucvol,usecprj,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mcprj,mkmem,mpw,my_natom,natom,nkpt,nsppol,ntypat
 integer,intent(in) :: pwind_alloc,usecprj
 real(dp),intent(in) :: ucvol
 real(dp),intent(inout) :: enefield,etotal
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(efield_type),intent(inout) :: dtefield
 type(hdr_type),intent(inout) :: hdr
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx1(natom),kg(3,mpw*mkmem),nattyp(ntypat)
 integer,intent(in) :: npwarr(nkpt),pwind(pwind_alloc,2,3)
 real(dp),intent(in) :: cg(2,mcg),gprimd(3,3)
 real(dp),intent(in) :: pel_cg(3),pwnsfac(2,pwind_alloc),rprimd(3,3)
 real(dp),intent(inout) :: pel(3),pelev(3),pion(3),xred(3,natom)
 type(pawcprj_type),intent(in) :: cprj(natom,mcprj*usecprj)
 type(pawrhoij_type), intent(in) :: pawrhoij(my_natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: my_nspinor,option,unit_out,iir,jjr,kkr
 real(dp) :: pdif_mod,eenth,ucvol_local
 character(len=500) :: message
!arrays
 real(dp) :: gmet(3,3),gprimdlc(3,3),pdif(3),ptot(3),red_ptot(3),rmet(3,3)
!! ptot(3) = total polarization (not reduced) REC
!! red_ptot(3) = internal reduced total polarization REC
 real(dp) ::  A(3,3),A1(3,3),A_new(3,3),efield_new(3)

! *************************************************************************

 DBG_ENTER("COLL")

 if (usecprj==0.and.psps%usepaw==1) then
   write (message,'(3a)')&
&   'cprj datastructure must be allocated !',ch10,&
&   'Action: change pawusecp input keyword.'
   MSG_ERROR(message)
 end if

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

 if(dtset%berryopt>0 .and. dtset%berryopt/=4 .and. dtset%berryopt/=6 .and. dtset%berryopt/=7 .and. &
& dtset%berryopt/=14 .and. dtset%berryopt/=16 .and. dtset%berryopt/=17)then  !!HONG

   if (dtset%berryopt==1 .or. dtset%berryopt==3) then
     call berryphase(atindx1,dtset%bdberry,cg,gprimd,dtset%istwfk,&
&     dtset%kberry,kg,dtset%kptns,dtset%kptopt,dtset%kptrlatt,&
&     mband,mcg,mkmem,mpw,natom,nattyp,dtset%nband,dtset%nberry,npwarr,&
&     my_nspinor,nsppol,psps%ntypat,nkpt,rprimd,ucvol,&
&     xred,psps%ziontypat)
   end if

   if (dtset%berryopt==2 .or. dtset%berryopt==3) then
     call uderiv(dtset%bdberry,cg,gprimd,hdr,dtset%istwfk,&
&     dtset%kberry,kg,dtset%kptns,dtset%kptopt,&
&     dtset%kptrlatt,mband,mcg,mkmem,mpi_enreg,mpw,&
&     natom,dtset%nband,dtset%nberry,npwarr,my_nspinor,nsppol,&
&     nkpt,dtfil%unddk,dtfil%fnameabo_1wf)
   end if

 else if(dtset%berryopt<0 .or. dtset%berryopt==4 .or. dtset%berryopt==6 .or. dtset%berryopt==7 .or.  &
&   dtset%berryopt==14 .or. dtset%berryopt==16 .or. dtset%berryopt==17)then   !!HONG

   select case (dtset%berryopt)
   case (-5)
     option = 2
   case (-3)
     option = 3
   case (-2)
     option = 2
   case (-1)
     option = 1
   case (4)
     option = 1
     pel(:) = zero
     pelev(:) = zero
   case (6)                !!HONG
     option = 1
     pel(:) = zero
     pelev(:) = zero
   case (7)                !!HONG
     option = 1
     pel(:) = zero
     pelev(:) = zero
   case (14)                !!HONG
     option = 1
     pel(:) = zero
     pelev(:) = zero
   case (16)                !!HONG
     option = 1
     pel(:) = zero
     pelev(:) = zero
   case (17)                !!HONG
     option = 1
     pel(:) = zero
     pelev(:) = zero
   end select

   unit_out = ab_out
   call berryphase_new(atindx1,cg,cprj,dtefield,dtfil,dtset,psps,&
&   gprimd,hdr,psps%indlmn,kg,&
&   psps%lmnmax,mband,mcg,mcprj,mkmem,mpi_enreg,mpw,my_natom,natom,npwarr,&
&   nsppol,psps%ntypat,nkpt,option,pawrhoij,&
&   pawtab,pel,pelev,pion,ptot,red_ptot,pwind,&                            !!REC
&  pwind_alloc,pwnsfac,rprimd,dtset%typat,ucvol,&
&   unit_out,usecprj,psps%usepaw,xred,psps%ziontypat)

   dtefield%red_ptot1(:)=red_ptot(:)

   if (dtset%berryopt == 4 .or. dtset%berryopt == 6 .or. dtset%berryopt == 7 .or.  &
&   dtset%berryopt == 14 .or. dtset%berryopt == 16 .or. dtset%berryopt == 17 ) then   !!HONG

!    Check if pel has the same value as pel_cg
!    if (psps%usepaw == 1) pel(:) = pel(:) + pelev(:) ! add on-site term for PAW
!    if (psps%usepaw == 1) red_ptot(:) = red_ptot(:) + pelev(:) ! add on-site term for PAW  !! REC
!    above line suppressed because in the PAW case, pel already includes all on-site
!    terms and pelev should not be added in additionally. We are computing pelev separately for
!    reporting purposes only.
!    13 June 2012 J Zwanziger

     pdif(:) = pel_cg(:) - pel(:)
     pdif_mod = pdif(1)**2 + pdif(2)**2 + pdif(3)**2

     if (pdif_mod > tol8) then
       write(message,'(11(a),e16.9)')ch10,&
&       ' scfcv (electric field calculation) : WARNING -',ch10,&
&       '   The difference between pel (electronic Berry phase updated ',ch10,&
&       '   at each SCF cycle)',ch10,&
&       '   and pel_cg (electronic Berryphase computed using the ',&
&       'berryphase routine) is',ch10,&
&       '   pdif_mod = ',pdif_mod
       call wrtout(std_out,message,'COLL')
       write(message,'(a,6(a,e16.9,a))') ch10,&
&       'pel_cg(1) = ',pel_cg(1),ch10,&
&       'pel_cg(2) = ',pel_cg(2),ch10,&
&       'pel_cg(3) = ',pel_cg(3),ch10,&
&       'pel(1) = ',pel(1),ch10,&
&       'pel(2) = ',pel(2),ch10,&
&       'pel(3) = ',pel(3),ch10
       MSG_ERROR(message)
     end if

!    Use this (more accurate) value of P to recompute enefield
     if (dtset%berryopt == 4 .or. dtset%berryopt == 14 ) then             !!HONG
       etotal = etotal - enefield

       enefield = -dot_product(dtset%red_efieldbar,red_ptot)
       call metric(gmet,gprimdlc,-1,rmet,rprimd,ucvol_local)
       eenth = zero
       do iir=1,3
         do jjr=1,3
           eenth= eenth+gmet(iir,jjr)*dtset%red_efieldbar(iir)*dtset%red_efieldbar(jjr)         !! HONG g^{-1})_ij ebar_i ebar_j
         end do
       end do
       eenth=-1_dp*(ucvol_local/(8.0d0*pi))*eenth
       enefield=enefield+eenth

       etotal = etotal + enefield

       write(message,'(a,a)')ch10,&
&       ' Stress tensor under a constant electric field:'
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')

     end if

!    ! In finite D-field case, turn it into internal energy    !!HONG
     if (dtset%berryopt == 6 .or. dtset%berryopt == 16 )  then
       etotal = etotal - enefield

       enefield=zero
       call metric(gmet,gprimdlc,-1,rmet,rprimd,ucvol_local)
       do iir=1,3
         do jjr=1,3
           enefield= enefield+gmet(iir,jjr)*dtset%red_efieldbar(iir)*dtset%red_efieldbar(jjr)         !! HONG g^{-1})_ij ebar_i ebar_j
         end do
       end do
       enefield= ucvol_local/(8.0d0*pi)*enefield

       etotal = etotal + enefield

       write(message,'(a,a)')ch10,&
&       ' Stress tensor under a constant electric displacement field:'
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')

     end if

!    HONG  calculate internal energy and electric enthalpy for mixed BC case.
     if ( dtset%berryopt == 17 ) then
       etotal = etotal - enefield
       enefield = zero

       call metric(gmet,gprimdlc,-1,rmet,rprimd,ucvol_local)
       A(:,:)=(4*pi/ucvol_local)*rmet(:,:)
       A1(:,:)=A(:,:)
       A_new(:,:)=A(:,:)
       efield_new(:)=dtset%red_efield(:)
       eenth = zero

       do kkr=1,3
         if (dtset%jfielddir(kkr)==1) then    ! fixed ebar direction
!          step 1 add -ebar*p
           eenth=eenth - dtset%red_efieldbar(kkr)*red_ptot(kkr)

!          step 2  chang to e_new (change e to ebar)
           efield_new(kkr)=dtset%red_efieldbar(kkr)

!          step 3  chang matrix A to A1

           do iir=1,3
             do jjr=1,3
               if (iir==kkr .and. jjr==kkr) A1(iir,jjr)=-1.0/A(kkr,kkr)
               if ((iir==kkr .and. jjr/=kkr) .or.  (iir/=kkr .and.  jjr==kkr)) &
&               A1(iir,jjr)=-1.0*A(iir,jjr)/A(kkr,kkr)
               if (iir/=kkr .and. jjr/=kkr) A1(iir,jjr)=A(iir,jjr)-A(iir,kkr)*A(kkr,jjr)/A(kkr,kkr)
             end do
           end do

           A(:,:)=A1(:,:)
           A_new(:,:)=A1(:,:)
         end if

       end do  ! end fo kkr


       do iir=1,3
         do jjr=1,3
           eenth= eenth+(1/2.0)*A_new(iir,jjr)*efield_new(iir)*efield_new(jjr)
         end do
       end do

       enefield=eenth
       etotal = etotal + enefield

       write(message,'(a,a)')ch10,&
&       ' Stress tensor under a constant (mixed) electric and electric displacement field:'
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')

     end if   ! berryopt==17


!    MVeithen: to clarify
!    Which stress tensor should be used in structural optimizations?
!    The one at constant electric field or at constant potential drop.
!    write(message,'(a,a)')ch10,&
!    &     ' Stress tensor imposing a constant electric field:'
!    call wrtout(std_out,message,'COLL')
!    call wrtout(ab_out,message,'COLL')

   end if ! dtset%berryopt == 4/6/7/14/16/17

 end if ! dtset%berryopt>0 or dtset%berryopt/=4/6/7/14/16/17

 DBG_EXIT("COLL")

end subroutine elpolariz
!!***

!!****f* ABINIT/uderiv
!! NAME
!! uderiv
!!
!! FUNCTION
!! This routine is called computes the derivative of
!! ground-state wavefunctions with respect to k (du/dk) by finite differencing
!! on neighbouring k points
!! Work for nsppol=1 or 2, but only accept nspinor=1,
!!
!! INPUTS
!!  bdberry(4)=band limits for Berry phase contributions (or du/dk)
!!   spin up and spin down (bdberry(3:4) is irrelevant when nsppol=1)
!!  cg(2,mcg)=planewave coefficients of wavefunctions
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  istwfk(nkpt_)=input option parameter that describes the storage of wfs
!!  kberry(3,20)= different delta k for Berry phases(or du/dk),
!!   in unit of kptrlatt only kberry(1:3,1:nberry) is relevant
!!  kg(3,mpw*mkmem)=reduced planewave coordinates
!!  kpt_(3,nkpt_)=reduced coordinates of k points generated by ABINIT,
!!               kpt_ samples half the BZ if time-reversal symetrie is used
!!  kptopt=2 when time-reversal symmetry is used
!!  kptrlatt(3,3)=k-point lattice specification
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mkmem=number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw
!!  natom=number of atoms in cell
!!  nband(nkpt*nsppol)=number of bands at each k point, for each polarization
!!  nberry=number of Berry phases(or du/dk) to be computed
!!  nkpt=number of k points
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  unddk=unit number for ddk file
!!
!! OUTPUT
!!  (the ddk wavefunctions are written on disk)
!!
!! SIDE EFFECTS
!!
!! TODO
!!  Cleaning, checking for rules
!!  Should allow for time-reversal symmetry (istwfk)
!!  WARNING : the use of nspinor is completely erroneous
!!
!! NOTES
!! Local Variables:
!!  cmatrix(:,:,:)= overlap matrix of size maxband*maxband
!!  cg_index(:,:,:)= unpacked cg index array for specific band,
!!   k point and polarization.
!!  det(2,2)= intermediate output of Lapack routine zgedi.f
!!  dk(3)= step taken to the next k mesh point along the kberry direction
!!  gpard(3)= dimensionalreciprocal lattice vector G along which the
!!          polarization is computed
!!  kg_kpt(:,:,:)= unpacked reduced planewave coordinates with subscript of
!!          planewave and k point
!!  kpt(3,nkpt)=reduced coordinates of k-point grid that samples the whole BZ
!!  kpt_flag(nkpt)=kpt_flag(ikpt)=0 when the wf was generated by the ABINIT
!!                 code
!!                 kpt_flag(ikpt) gives the indices of the k-point
!!                 related to ikpt by time revers
!!  maxband/minband= control the minimum and maximum band calculated in the
!!           overlap matrix
!!  npw_k= npwarr(ikpt), number of planewaves in basis at this k point
!!  shift_g_2(nkpt,nkpt)= .true. if the k point should be shifted by a G vector;
!!  .false. if not
!!  tr(2)=variable that changes k to -k
!!                              G to -G
!!                     $c_g$ to $c_g^*$ when time-reversal symetrie is used
!!
!! PARENTS
!!      elpolariz
!!
!! CHILDREN
!!      appdig,dzgedi,dzgefa,hdr_io,matr3inv,rwwf,waveformat,wffclose,wffopen
!!      wrtout,xdefineoff
!!
!! SOURCE

subroutine uderiv(bdberry,cg,gprimd,hdr,istwfk,kberry,kg,kpt_,kptopt,kptrlatt,&
& mband,mcg,mkmem,mpi_enreg,mpw,natom,nband,nberry,npwarr,nspinor,nsppol,nkpt_,&
& unddk,fnameabo_1wf)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kptopt,mband,mcg,mkmem,mpw,natom,nberry,nkpt_,nspinor
 integer,intent(in) :: nsppol,unddk
 type(MPI_type),intent(in) :: mpi_enreg
 type(hdr_type),intent(inout) :: hdr
!arrays
 integer,intent(in) :: bdberry(4),istwfk(nkpt_),kberry(3,20),kg(3,mpw*mkmem)
 integer,intent(in) :: kptrlatt(3,3),nband(nkpt_*nsppol),npwarr(nkpt_)
 real(dp),intent(in) :: cg(2,mcg),gprimd(1:3,1:3)
 real(dp),intent(in) :: kpt_(3,nkpt_)
 character(len=fnlen),intent(in) :: fnameabo_1wf

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 integer :: iomode,band_in,cg_index_iband,fform,flag1
 integer :: formeig,iband,iberry,icg,idir,ierr,ifor,ii,ikpt,ikpt2,ikpt_
 integer :: index,index1,info,ipert,ipw,isgn,isppol,jband,jj,jkpt,jkpt_
 integer :: maxband,mcg_disk,me,minband,nband_diff,nband_k
 integer :: nkpt,npw_k,pertcase,rdwr,read_k,spaceComm
 integer :: tim_rwwf
 real(dp) :: gmod,twodk
 character(len=500) :: message
 character(len=fnlen) :: fiwf1o
 type(wffile_type) :: wffddk
!arrays
 integer :: kg_jl(0,0,0),kpt_flag(2*nkpt_)
 integer,allocatable :: cg_index(:,:,:),ikpt_dk(:,:),ipvt(:)
 integer,allocatable :: kg_kpt(:,:,:)
 real(dp) :: det(2,2),diffk(3),diffk2(3),dk(3),gpard(3),klattice(3,3)
 real(dp) :: kptrlattr(3,3),tr(2)
 real(dp) :: cg_disk(0,0,0)
 real(dp),allocatable :: cmatrix(:,:,:),dudk(:,:)
 real(dp),allocatable :: eig_dum_2(:),kpt(:,:)
 real(dp),allocatable :: occ_dum_2(:),phi(:,:,:),u_tilde(:,:,:,:),zgwork(:,:)
 logical,allocatable :: shift_g_2(:,:)

! *************************************************************************

 if(min(2,(1+mpi_enreg%paral_spinor)*nspinor)==2)then
   MSG_ERROR('uderiv: does not yet work for nspinor=2')
 end if

 if(maxval(istwfk(:))/=1)then
   write(message,'(3a)')&
&   'Sorry, this routine does not work yet with istwfk/=1.',ch10,&
&   'This should have been tested previously ...'
   MSG_BUG(message)
 end if

 if (kptopt==3) then
   nkpt = nkpt_
   ABI_ALLOCATE(kpt,(3,nkpt))
   kpt(:,:)=kpt_(:,:)
 else if (kptopt==2) then
   nkpt = nkpt_*2
   ABI_ALLOCATE(kpt,(3,nkpt))
   do ikpt = 1,nkpt/2
     kpt_flag(ikpt) = 0
     kpt(:,ikpt)=kpt_(:,ikpt)
   end do
   index = 0
   do ikpt = (nkpt/2+1),nkpt
     flag1 = 0
     do jkpt = 1, nkpt/2
       if (((abs(kpt_(1,ikpt-nkpt/2)+kpt_(1,jkpt))<1.0d-8).or.&
&       (abs(1-abs(kpt_(1,ikpt-nkpt/2)+kpt_(1,jkpt)))<1.0d-8))&
&       .and.((abs(kpt_(2,ikpt-nkpt/2)+kpt_(2,jkpt))<1.0d-8).or.&
&       (abs(1-abs(kpt_(2,ikpt-nkpt/2)+kpt_(2,jkpt)))<1.0d-8))&
&       .and.((abs(kpt_(3,ikpt-nkpt/2)+kpt_(3,jkpt))<1.0d-8).or.&
&       (abs(1-abs(kpt_(3,ikpt-nkpt/2)+kpt_(3,jkpt)))<1.0d-8))) then
         flag1 = 1
         index = index + 1
         exit
       end if
     end do
     if (flag1==0) then
       kpt_flag(ikpt-index)=ikpt-nkpt/2
       kpt(:,ikpt-index)=-kpt_(:,ikpt-nkpt/2)
     end if
   end do
   nkpt = nkpt - index
 end if

!DEBUG
!write(101,*) 'beginning write kpt'
!do ikpt=1,nkpt
!write(101,*) kpt(:,ikpt)
!end do
!ENDDEBUG

 ABI_ALLOCATE(shift_g_2,(nkpt,nkpt))

!Compute primitive vectors of the k point lattice
!Copy to real(dp)
 kptrlattr(:,:)=kptrlatt(:,:)
!Go to reciprocal space (in reduced coordinates)
 call matr3inv(kptrlattr,klattice)

 do iberry=1,nberry

!  **************************************************************************
!  Determine the appended index for ddk 1WF files

   do idir=1,3
     if (kberry(idir,iberry) ==1) then
       ipert=natom+1
       pertcase=idir+(ipert-1)*3
     end if
   end do

!  open ddk 1WF file
   formeig=1

   call appdig(pertcase,fnameabo_1wf,fiwf1o)
   !call wfk_open_read(wfk, fiwf1o, formeig, iomode, unddk, spaceComm)

   spaceComm=xmpi_comm_self; me=0 ; iomode=IO_MODE_FORTRAN
   call WffOpen(iomode,spaceComm,fiwf1o,ierr,wffddk,master,me,unddk)

   rdwr=2 ; fform=2
   call hdr_io(fform,hdr,rdwr,wffddk)

!  Define offsets, in case of MPI I/O
   call xdefineOff(formeig,wffddk,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt_)

!  *****************************************************************************
!  Calculate dimensional recip lattice vector along which P is calculated
!  dk =  step to the nearest k point along that direction
!  in reduced coordinates

   dk(:)=dble(kberry(1,iberry))*klattice(:,1)+&
&   dble(kberry(2,iberry))*klattice(:,2)+&
&   dble(kberry(3,iberry))*klattice(:,3)

   do idir=1,3
     if (dk(idir)/=0) then
       twodk=2*dk(idir)
     end if
   end do

   gpard(:)=dk(1)*gprimd(:,1)+dk(2)*gprimd(:,2)+dk(3)*gprimd(:,3)
   gmod=sqrt(dot_product(gpard,gpard))

!  ******************************************************************************
!  Select the k grid  points along the direction to compute dudk
!  dk =  step to the nearest k point along that direction

!  For each k point, find k_prim such that k_prim= k + dk mod(G)
!  where G is a vector of the reciprocal lattice
   ABI_ALLOCATE(ikpt_dk,(2,nkpt))
   ikpt_dk(1:2,1:nkpt)=0
   shift_g_2(:,:)= .false.

   do ikpt=1,nkpt
     do ikpt2=1,nkpt
       diffk(:)=abs(kpt(:,ikpt2)-kpt(:,ikpt)-dk(:))
       diffk2(:)=abs(kpt(:,ikpt2)-kpt(:,ikpt)+dk(:))
       if (sum(abs(diffk(:)-nint(diffk(:))))<3*tol8)then
         ikpt_dk(1,ikpt)=ikpt2
         if(sum(diffk(:))>=3*tol8) shift_g_2(ikpt,ikpt2) = .true.
       end if
       if (sum(abs(diffk2(:)-nint(diffk2(:))))<3*tol8)then
         ikpt_dk(2,ikpt)=ikpt2
         if(sum(diffk2(:))>=3*tol8) shift_g_2(ikpt,ikpt2) = .true.
       end if
     end do
   end do

   write(message,'(a,a,a,3f9.5,a,a,3f9.5,a)')ch10,&
&   ' Computing the derivative for reciprocal vector:',ch10,&
&   dk(:),' (in reduced coordinates)',ch10,&
&   gpard(1:3),' (in cartesian coordinates - atomic units)'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   if(nsppol==1)then
     write(message, '(a,i5,a,i5)')&
&     ' From band number',bdberry(1),'  to band number',bdberry(2)
   else
     write(message, '(a,i5,a,i5,a,a,a,i5,a,i5,a)')&
&     ' From band number',bdberry(1),'  to band number',bdberry(2),' for spin up,',&
&     ch10,&
&     ' from band number',bdberry(3),'  to band number',bdberry(4),' for spin down.'
   end if
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

!  *****************************************************************************
   ABI_ALLOCATE(dudk,(2,mpw*nspinor*mband*nsppol))
   ABI_ALLOCATE(eig_dum_2,((2*mband)**formeig*mband))
   ABI_ALLOCATE(occ_dum_2,((2*mband)**formeig*mband))
   dudk(1:2,:)=0.0_dp
   eig_dum_2=0.0_dp
   occ_dum_2=0.0_dp

   if (mkmem/=0) then

!    Find the location of each wavefunction

     ABI_ALLOCATE(cg_index,(mband,nkpt_,nsppol))
     icg = 0
     do isppol=1,nsppol
       do ikpt=1,nkpt_
         nband_k=nband(ikpt+(isppol-1)*nkpt_)
         npw_k=npwarr(ikpt)
         do iband=1,nband_k
           cg_index(iband,ikpt,isppol)=(iband-1)*npw_k*nspinor+icg
         end do
         icg=icg+npw_k*nspinor*nband(ikpt)
       end do
     end do

!    Find the planewave vectors for each k point
!    SHOULD BE REMOVED WHEN ANOTHER INDEXING TECHNIQUE WILL BE USED FOR kg
     ABI_ALLOCATE(kg_kpt,(3,mpw*nspinor,nkpt_))
     kg_kpt(:,:,:) = 0
     index1 = 0
     do ikpt=1,nkpt_
       npw_k=npwarr(ikpt)
       do ipw=1,npw_k*nspinor
         kg_kpt(1:3,ipw,ikpt)=kg(1:3,ipw+index1)
       end do
       index1=index1+npw_k*nspinor
     end do
   end if

!  *************************************************************************
!  Loop over spins
   do isppol=1,nsppol

     minband=bdberry(2*isppol-1)
     maxband=bdberry(2*isppol)

     if(minband<1)then
       write(message,'(a,i0,a)')'  The band limit minband= ',minband,', is lower than 0.'
       MSG_BUG(message)
     end if

     if(maxband<1)then
       write(message,'(a,i0,a)')' The band limit maxband= ',maxband,', is lower than 0.'
       MSG_BUG(message)
     end if

     if(maxband<minband)then
       write(message,'(a,i0,a,i0)')' maxband= ',maxband,', is lower than minband= ',minband
       MSG_BUG(message)
     end if

!    Loop over k points
     do ikpt_=1,nkpt_

       read_k = 0

       ikpt=ikpt_
       tr(1) = 1.0_dp

       if (kptopt==2) then
         if (read_k == 0) then
           if (kpt_flag(ikpt_)/=0) then
             tr(1) = -1.0_dp
             ikpt= kpt_flag(ikpt_)
           end if
         else           !read_k
           if (kpt_flag(ikpt_)/=0) then
             tr(-1*read_k+3) = -1.0_dp
             ikpt= kpt_flag(ikpt_)
           end if
         end if       !read_k
       end if           !kptopt

       nband_k=nband(ikpt+(isppol-1)*nkpt_)

       if(nband_k<maxband)then
         write(message,'(a,i0,a,i0)')'  maxband=',maxband,', is larger than nband(i,isppol)=',nband_k
         MSG_BUG(message)
       end if

       npw_k=npwarr(ikpt)

       ABI_ALLOCATE(u_tilde,(2,npw_k,maxband,2))
       u_tilde(1:2,1:npw_k,1:maxband,1:2)=0.0_dp

!      ifor = 1,2 represents forward and backward neighbouring k points of ikpt
!      respectively along dk direction

       do ifor=1,2

         ABI_ALLOCATE(phi,(2,mpw,mband))
         ABI_ALLOCATE(cmatrix,(2,maxband,maxband))
         phi(1:2,1:mpw,1:mband)=0.0_dp; cmatrix(1:2,1:maxband,1:maxband)=0.0_dp

         isgn=(-1)**ifor
         jkpt_= ikpt_dk(ifor,ikpt_)

         tr(2) = 1.0_dp

         jkpt=jkpt_

         if (kptopt==2) then
           if (read_k == 0) then
             if (kpt_flag(jkpt_)/=0) then
               tr(2) = -1.0_dp
               jkpt= kpt_flag(jkpt_)
             end if
           else           !read_k
             if (kpt_flag(jkpt_)/=0) then
               tr(read_k) = -1.0_dp
               jkpt= kpt_flag(jkpt_)
             end if
           end if       !read_k
         end if           !kptopt

         if (ifor==1) read_k = 2

         jj = read_k
         ii = -1*read_k+3

         call waveformat(cg,cg_disk,cg_index,phi,dk,ii,ikpt,&
&         ikpt_,isgn,isppol,jj,jkpt,jkpt_,kg_kpt,kpt,kg_jl,maxband,mband,mcg,mcg_disk,&
&         minband,mkmem,mpw,nkpt,nkpt_,npwarr,nsppol,nspinor,shift_g_2,tr)

!        Compute the overlap matrix <u_k|u_k+b>

         do iband=minband,maxband
           cg_index_iband=cg_index(iband,ikpt,isppol)
           do jband=minband,maxband
             do ipw=1,npwarr(ikpt)
               cmatrix(1,iband,jband)=cmatrix(1,iband,jband)+&
&               cg(1,ipw+cg_index_iband)*phi(1,ipw,jband)+&
&               tr(ii)*cg(2,ipw+cg_index_iband)*tr(jj)*phi(2,ipw,jband)

               cmatrix(2,iband,jband)=cmatrix(2,iband,jband)+&
&               cg(1,ipw+cg_index_iband)*tr(jj)*phi(2,ipw,jband)-&
&               tr(ii)*cg(2,ipw+cg_index_iband)*phi(1,ipw,jband)
             end do
           end do
         end do

!        Compute the inverse of cmatrix(1:2,minband:maxband, minband:maxband)

         band_in = maxband - minband + 1
         ABI_ALLOCATE(ipvt,(maxband))
         ABI_ALLOCATE(zgwork,(2,1:maxband))

!        Last argument of zgedi means calculate inverse only
         call dzgefa(cmatrix(1,minband,minband),maxband, band_in,ipvt,info)
         call dzgedi(cmatrix(1,minband,minband),maxband, band_in,ipvt,det,zgwork,01)

         ABI_DEALLOCATE(zgwork)
         ABI_DEALLOCATE(ipvt)

!        Compute the product of Inverse overlap matrix with the wavefunction

         do iband=minband,maxband
           do ipw=1,npwarr(ikpt)
             u_tilde(1,ipw,iband,ifor)= &
&             dot_product(cmatrix(1,minband:maxband,iband),&
&             phi(1,ipw,minband:maxband))-&
&             dot_product(cmatrix(2,minband:maxband,iband),&
&             tr(jj)*phi(2,ipw,minband:maxband))
             u_tilde(2,ipw,iband,ifor)= &
&             dot_product(cmatrix(1,minband:maxband,iband),&
&             tr(jj)*phi(2,ipw,minband:maxband))+&
&             dot_product(cmatrix(2,minband:maxband,iband),&
&             phi(1,ipw,minband:maxband))
           end do
         end do
         ABI_DEALLOCATE(cmatrix)
         ABI_DEALLOCATE(phi)

       end do !ifor

!      Compute dudk for ikpt

       npw_k=npwarr(ikpt)

       do iband=minband,maxband

         icg=(iband-minband)*npw_k

         dudk(1,1+icg:npw_k+icg)=(u_tilde(1,1:npw_k,iband,1)-&
&         u_tilde(1,1:npw_k,iband,2))/twodk

         dudk(2,1+icg:npw_k+icg)=(u_tilde(2,1:npw_k,iband,1)-&
&         u_tilde(2,1:npw_k,iband,2))/twodk

       end do

       tim_rwwf=0
       mcg_disk=mpw*nspinor*mband
       nband_diff=maxband-minband+1
       call rwwf(dudk,eig_dum_2,formeig,0,0,ikpt,isppol,kg_kpt(:,:,ikpt),&
&       mband,mcg_disk,mpi_enreg,nband_diff,nband_diff,&
&       npw_k,nspinor,occ_dum_2,2,1,tim_rwwf,wffddk)

       !call wfk_read_band_block(wfk, band_block, ikpt, isppol, sc_mode,
       !  kg_k=kg_kpt(:,:,ikpt), cg_k=dudk, eig_k=eig_dum, occ_k=occ_dum)

       ABI_DEALLOCATE(u_tilde)

     end do !ikpt
   end do  !isppol

   ABI_DEALLOCATE(eig_dum_2)
   ABI_DEALLOCATE(occ_dum_2)
   ABI_DEALLOCATE(dudk)

   call WffClose(wffddk,ierr)
   !call wfk_close(wfk)

   ABI_DEALLOCATE(kg_kpt)
   ABI_DEALLOCATE(cg_index)
   ABI_DEALLOCATE(ikpt_dk)

 end do ! iberry

 ABI_DEALLOCATE(shift_g_2)
 ABI_DEALLOCATE(kpt)

 write(std_out,*) 'uderiv:  exit '

end subroutine uderiv
!!***


!!****f* ABINIT/waveformat
!! NAME
!! waveformat
!!
!! FUNCTION
!! This routine is to find the matched pairs of plane waves between
!! two neighbouring k points and load a new pw coefficients array cg_new
!! Was written first by Na Sai (thanks), but unfortunately without
!! any comment ...
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)= input planewave coefficients, in case mkmem/=0
!!  cg_disk(2,mpw*nspinor*mband,2)= input planewave coefficients, in case mkmem==0
!!  cg_index(mband,nkpt_,nsppol)=index of wavefunction iband,ikpt,isppol in the array cg.
!!  dk(3)= step taken to the next k mesh point along the kberry direction (see also isgn)
!!  ii=(to be documented)
!!  ikpt=index of the first k-point in the reduced Brillouin zone
!!  ikpt_=index of the first k-point in the full Brillouin zone
!!  isgn=1 if dk(3) is connecting the k-points (ikpt_ and jkpt)
!!      =-1 if -dk(3) is connecting the k-points
!!  isppol=1 if spin-up, =2 if spin-down
!!  jj=(to be documented)
!!  jkpt=index of the second k-point in the reduced Brillouin zone
!!  jkpt_=index of the second k-point in the full Brillouin zone
!!  kg_kpt(:,:,:)= unpacked reduced planewave coordinates with subscript of
!!          planewave and k point
!!  kpt(3,nkpt)=reduced coordinates of k-point grid that samples the whole BZ
!!  kg_jl(3,mpw,2)=(to be documented)
!!  maxband/minband= control the minimum and maximum band calculated in the
!!           overlap matrix
!!  mband=maximum number of bands (dimension of several cg* arrays)
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcg_disk=size of wave-functions array (cg_disk) =mpw*nspinor*mband
!!  mkmem= if 0, the wavefunctions are input in cg_disk, otherwise in cg
!!  mpw=maximum number of planewaves (dimension of several cg* arrays)
!!  nkpt=number of k points (full Brillouin zone !?!)
!!  nkpt_=number of k points (reduced Brillouin zone !?!)
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  shift_g_2(nkpt,nkpt)=non-zero if a G vector along kberry is needed to connect k points
!!  tr(2)=variable that changes k to -k
!!                              G to -G
!!                     $c_g$ to $c_g^*$ when time-reversal symetrie is used
!!
!! OUTPUT
!!  cg_new(2,mpw,maxband)=planewave coefficients transferred onto the
!!   set of planewaves at k
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      uderiv
!!
!! CHILDREN
!!
!! SOURCE

subroutine waveformat(cg,cg_disk,cg_index,cg_new,dk,ii,ikpt,&
& ikpt_,isgn,isppol,jj,jkpt,jkpt_,kg_kpt,kpt,kg_jl,maxband,mband,mcg,mcg_disk,&
& minband,mkmem,mpw,nkpt,nkpt_,npwarr,nsppol,nspinor,shift_g_2,tr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ii,ikpt,ikpt_,isgn,isppol,jj,jkpt,jkpt_,maxband,mband,mcg,mcg_disk
 integer,intent(in) :: minband,mkmem,mpw,nkpt,nkpt_,nspinor,nsppol
!arrays
 integer,intent(in) :: cg_index(mband,nkpt_,nsppol),kg_jl(3,mpw,2)
 integer,intent(in) :: kg_kpt(3,mpw*nspinor,nkpt_),npwarr(nkpt_)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: cg_disk(2,mcg_disk,2),dk(3),kpt(3,nkpt),tr(2)
 real(dp),intent(out) :: cg_new(2,mpw,maxband)
 logical,intent(in) :: shift_g_2(nkpt,nkpt)

!Local variables -------------------------
!scalars
 integer :: cg_index_iband,iband,ipw,jpw,nomatch,npw_k
 logical :: found_match
!arrays
 integer :: dg(3)

! ***********************************************************************

 npw_k=npwarr(ikpt)


 nomatch=0

!If there is no shift of G-vector between ikpt_ and jkpt_
 if(shift_g_2(ikpt_,jkpt_) .eqv. .false.) then

!  DEBUG
!  write(111,*)'pair', ikpt_,jkpt_,'noshift'
!  ENDDEBUG

!  If the original wavefunction is contained in cg_disk
   if(mkmem==0) then

     do ipw=1,npw_k

       found_match = .false.

       do jpw=1,npwarr(jkpt)
         if (sum(abs(tr(ii)*kg_jl(:,ipw,ii)-tr(jj)*kg_jl(:,jpw,jj)))<3*tol8)then
           do iband=minband, maxband
             cg_index_iband=(iband-1)*npwarr(jkpt)
             cg_new(1:2,ipw,iband)=cg_disk(1:2,jpw+cg_index_iband,jj)
           end do
           found_match = .true.
           exit
         end if
       end do

       if (found_match .eqv. .false.) then
         do iband=minband,maxband
           cg_new(1:2,ipw,iband)=zero
         end do
         nomatch = nomatch + 1
       end if

     end do

!    Here, the wavefunctions are contained in cg
   else

     do ipw=1,npw_k

       found_match = .false.

       do jpw=1,npwarr(jkpt)
         if (sum(abs(tr(ii)*kg_kpt(:,ipw,ikpt)-tr(jj)*kg_kpt(:,jpw,jkpt)))<3*tol8)then
           do iband=minband, maxband
             cg_index_iband=cg_index(iband,jkpt,isppol)
             cg_new(1:2,ipw,iband)=cg(1:2,jpw+cg_index_iband)
           end do
           found_match = .true.
           exit
         end if
       end do

       if (found_match .eqv. .false.) then
         do iband=minband,maxband
           cg_new(1:2,ipw,iband)=(0.0_dp,0.0_dp)
         end do
         nomatch = nomatch + 1
       end if
     end do

   end if

!  DEBUG
!  write(111,*) 'normal pair nomatch=',nomatch
!  ENDDEBUG

!  If there is a G-vector shift between ikpt_ and jkpt_
 else

!  DEBUG
!  write(111,*) 'pair',ikpt_,jkpt_,' need shift'
!  ENDDEBUG

   dg(:) = -1*nint(tr(jj)*kpt(:,jkpt)-tr(ii)*kpt(:,ikpt)+isgn*dk(:))

!  If the original wavefunction is contained in cg_disk
   if(mkmem==0) then

     do ipw=1,npw_k

       found_match = .false.

       do jpw=1,npwarr(jkpt)
         if (sum(abs(tr(ii)*kg_jl(:,ipw,ii)-(tr(jj)*kg_jl(:,jpw,jj)-&
&         dg(:))))<3*tol8)then

           do iband=minband, maxband
             cg_index_iband=(iband-1)*npwarr(jkpt)
             cg_new(1:2,ipw,iband)=cg_disk(1:2,jpw+cg_index_iband,jj)
           end do
           found_match = .true.
           exit
         end if
       end do

       if (found_match .eqv. .false.) then
         do iband=minband,maxband
           cg_new(1:2,ipw,iband)=(0.0_dp,0.0_dp)
         end do
         nomatch = nomatch + 1
       end if
     end do

!    Here, the wavefunctions are contained in cg
   else

     do ipw=1,npw_k

       found_match = .false.

       do jpw=1,npwarr(jkpt)
         if (sum(abs(tr(ii)*kg_kpt(:,ipw,ikpt)-(tr(jj)*kg_kpt(:,jpw,jkpt)-&
&         dg(:))))<3*tol8)then
           do iband=minband, maxband
             cg_index_iband=cg_index(iband,jkpt,isppol)
             cg_new(1:2,ipw,iband)=cg(1:2,jpw+cg_index_iband)
           end do
           found_match = .true.
           exit
         end if
       end do

       if (found_match .eqv. .false.) then
         do iband=minband,maxband
           cg_new(1:2,ipw,iband)=zero
         end do
         nomatch = nomatch + 1
       end if
     end do

   end if

!  DEBUG
!  write(111,*) 'special pair nomatch=',nomatch
!  ENDDEBUG

 end if

end subroutine waveformat
!!***

end module m_elpolariz
!!***
