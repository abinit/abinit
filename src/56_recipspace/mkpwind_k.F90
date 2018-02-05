!{\src2tex{textfont=tt}}
!****f* ABINIT/mkpwind_k
!! NAME
!! mkpwind_k
!!
!! FUNCTION
!! Make plane wave index at k point for basis at second k point,
!! needed to compute overlaps $\langle u_{k,n}|u_{k+b,n}\rangle$
!! as appear in Berry phase derived quantities
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! dk(3)=real vector difference of ket kpt - bra kpt
!! dtset <type(dataset_type)>=all input variables in this dataset
!! fnkpt=number of kpts in full BZ
!! fkptns=kpts in full BZ
!! gmet(3,3)=metric in reciprocal space
!! indkk_f2ibz(fnkpt,6)=information on folding from FBZ to IBZ (see initberry or initorbmag)
!! ikpt=index of bra k pt
!! ikpt1=index of neighbour ket k pt
!! kg(3,dtset%mpw*dtset%mkmem)=planewave basis data
!! kgindex(dtset%nkpt)= index of kg per kpt
!! mpi_enreg=information about MPI parallelization
!! npw_k=number of planewaves at k
!! pwind_k1(dtset%mpw)=output index of ikpt1 basis states refered to ikpt
!! symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!   reciprocal space primitive translations
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
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

subroutine mkpwind_k(dk,dtset,fnkpt,fkptns,gmet,indkk_f2ibz,ikpt,ikpt1,&
     & kg,kgindex,mpi_enreg,npw_k,pwind_k1,symrec)

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_xmpi
 use m_errors
 use m_profiling_abi

 use m_fftcore, only : kpgsph

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkpwind_k'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 !scalars
 integer,intent(in) :: fnkpt,ikpt,ikpt1,npw_k
 type(dataset_type),intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg

 !arrays
 integer,intent(in) :: indkk_f2ibz(fnkpt,6),kg(3,dtset%mpw*dtset%mkmem),kgindex(dtset%nkpt)
 integer,intent(in) :: symrec(3,3,dtset%nsym)
 integer,intent(out) :: pwind_k1(dtset%mpw)
 real(dp),intent(in) :: dk(3),fkptns(3,fnkpt),gmet(3,3)
 
 !Local variables -------------------------
 !scalars
 integer :: exchn2n3d,idum1,ikg1,ipw,istwf_k,isym,isym1,jpw,npw_k1
 real(dp) :: ecut_eff
 
 !arrays
 integer,allocatable :: kg1_k(:,:)
 real(dp) :: dg(3),dum33(3,3),kpt1(3),iadum(3),iadum1(3)

 ! ***********************************************************************

 ABI_ALLOCATE(kg1_k,(3,dtset%mpw))
 
 ecut_eff = dtset%ecut*(dtset%dilatmx)**2
 exchn2n3d = 0 ; istwf_k = 1 ; ikg1 = 0

 ! Build basis sphere of plane waves for the nearest neighbour of the k-point 

 kg1_k(:,:) = 0
 kpt1(:) = dtset%kptns(:,ikpt1)
 call kpgsph(ecut_eff,exchn2n3d,gmet,ikg1,ikpt1,istwf_k,kg1_k,kpt1,1,mpi_enreg,dtset%mpw,npw_k1)

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
                   
 isym  = indkk_f2ibz(ikpt,2)
 isym1 = indkk_f2ibz(ikpt1,2)

 !        Construct transformed G vector that enters the matching condition:
 !        alpha(k) S(k)^{t,-1} ( -G(b) - GBZ(k) + G(k) )

 dg(:) = -indkk_f2ibz(ikpt,3:5) &
      &         -nint(-fkptns(:,ikpt) - dk(:) - tol10 + &
      &         fkptns(:,ikpt1)) &
      &         +indkk_f2ibz(ikpt1,3:5)

 iadum(:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),dg(:))

 dg(:) = iadum(:)
 
 !        Construct S(k)^{t,-1} S(b)^{t}

 dum33(:,:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),symrec(:,:,isym))

 !        Construct alpha(k) alpha(b)
 
 pwind_k1(:) = 0
 do ipw = 1, npw_k

    !          NOTE: the bra G vector is taken for the sym-related IBZ k point,
    !          not for the FBZ k point
    iadum(:) = kg(:,kgindex(ikpt) + ipw)
    
    !          to determine r.l.v. matchings, we transformed the bra vector
    !          Rotation
    iadum1(:)=0
    do idum1=1,3
       iadum1(:)=iadum1(:)+dum33(:,idum1)*iadum(idum1)
    end do
    iadum(:)=iadum1(:)
    iadum(:) = iadum(:) + dg(:)

    do jpw = 1, npw_k1
       iadum1(1:3) = kg1_k(1:3,jpw)
       if ( (iadum(1) == iadum1(1)).and. &
            &             (iadum(2) == iadum1(2)).and. &
            &             (iadum(3) == iadum1(3)) ) then
          pwind_k1(ipw) = jpw
          ! write(std_out,'(a,2i4)')'JWZ debug : bg ipw == jpw ',ipw,jpw
          exit
       end if
    end do
 end do

 ABI_DEALLOCATE(kg1_k)

end subroutine mkpwind_k
!!***
