!{\src2tex{textfont=tt}}
!!****f* ABINIT/elpolariz
!! NAME
!! elpolariz
!!
!! FUNCTION
!! Calculate corrections to total energy from polarising
!! electric field with or without Berry phases (berryopt keyword)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine elpolariz(atindx1,cg,cprj,dtefield,dtfil,dtset,etotal,enefield,gprimd,hdr,&
& kg,mband,mcg,mcprj,mkmem,mpi_enreg,mpw,my_natom,natom,nattyp,nkpt,&
& npwarr,nsppol,ntypat,pawrhoij,pawtab,&
& pel,pel_cg,pelev,pion,psps,pwind,pwind_alloc,&
& pwnsfac,rprimd,ucvol,usecprj,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_efield

 use m_pawtab,   only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type
 use m_pawcprj,  only : pawcprj_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elpolariz'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_67_common
!End of the abilint section

 implicit none

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
