!!****m* ABINIT/m_wfd_optic
!! NAME
!!  m_wfd_optic
!!
!! FUNCTION
!!  Functions to compute optical matrix elements using the wavefunction descriptor.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MG)
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

module m_wfd_optic

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi

 use defs_datatypes,      only : ebands_t, pseudopotential_type
 use m_hide_lapack,       only : matrginv
 use m_bz_mesh,           only : kmesh_t, get_BZ_item
 use m_crystal,           only : crystal_t
 use m_vkbr,              only : vkbr_t, vkbr_free, vkbr_init, nc_ihr_comm
 use m_wfd,               only : wfd_t, wave_t
 use m_pawtab,            only : pawtab_type
 use m_pawcprj,           only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_paw_hr,            only : pawhur_t, paw_ihr

 implicit none

 private
!!***

 public :: calc_optical_mels
!!***

contains
!!***

!!****f* ABINIT/calc_optical_mels
!! NAME
!!  calc_optical_mels
!!
!! FUNCTION
!!  Calculate all optical matrix elements in the BZ.
!!
!! INPUTS
!! lomo_spin(Wfd%nsppol)=Index of the lomo band for the different spins.
!! lomo_min,max_band=minimum and max band index to be calculated.
!! nkbz=Number of points in the full Brillouin zone.
!! inclvkb=if different from 0, [Vnl,r] is included in the calculation of the
!!   matrix element of the velocity operator. No meaning for PAW (except for DFT+U)
!! qpt(3)
!! Kmesh<kmesh_t>=Info on the k-point sampling for wave functions.
!! Cryst<crystal_t>=Structure defining the crystalline structure.
!! KS_Bst<ebands_t>
!! Pawtab(Cryst%ntypat*usepaw)<pawtab_type>=PAW tabulated starting data
!! Psps <pseudopotential_type>=variables related to pseudopotentials.
!! Hur(Cryst%natom*usepaw)<pawhur_t>=Only for PAW and DFT+U, quantities used to evaluate the commutator [H_u,r].
!! Wfd<wfd_t>=Handler for the wavefunctions.
!!
!! OUTPUT
!! opt_cvk(lomo_min:max_band,lomo_min:max_band,nkbz,nsppol)=Matrix elements <c k|e^{+iqr}|v k>
!!
!! PARENTS
!!      m_exc_spectra,m_haydock
!!
!! CHILDREN
!!      get_bz_item,matrginv,pawcprj_alloc,pawcprj_free,vkbr_free,vkbr_init
!!      wfd%distribute_bbp,wfd%get_cprj,wrtout,xmpi_barrier,xmpi_sum
!!
!! SOURCE

subroutine calc_optical_mels(Wfd,Kmesh,KS_Bst,Cryst,Psps,Pawtab,Hur,&
&  inclvkb,lomo_spin,lomo_min,max_band,nkbz,qpoint,opt_cvk)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkbz,inclvkb,lomo_min,max_band
 type(kmesh_t),intent(in) :: Kmesh
 type(crystal_t),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
 type(ebands_t),intent(in) :: KS_Bst
 type(wfd_t),target,intent(inout) :: Wfd
!arrays
 integer,intent(in) :: lomo_spin(Wfd%nsppol)
 real(dp),intent(in) :: qpoint(3)
 complex(dpc),intent(out) :: opt_cvk(lomo_min:max_band,lomo_min:max_band,nkbz,Wfd%nsppol)
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(pawhur_t),intent(in) :: Hur(Cryst%natom*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer :: nsppol,usepaw,nspinor,comm,spin,npw_k,istwf_k,my_nbbp
 integer :: ik_bz,ik_ibz,itim_k,isym_k,ib_c,ib_v,ierr,my_rank
 real(dp) :: ediff
 complex(dpc) :: emcvk
 character(len=500) :: msg
 type(vkbr_t) :: vkbr
 type(wave_t),pointer :: wave_v, wave_c
!arrays
 integer,allocatable :: bbp_distrb(:,:)
 integer,ABI_CONTIGUOUS pointer :: kg_k(:,:)
 real(dp) :: mat_dp(3,3),qrot(3),b1(3),b2(3),b3(3),kbz(3)
 complex(dpc),allocatable :: ir_kibz(:,:,:,:,:)
 complex(gwpc), ABI_CONTIGUOUS pointer :: ug_c(:),ug_v(:)
 complex(gwpc) :: ihrc(3,Wfd%nspinor**2)
 logical :: bbp_mask(Wfd%mband,Wfd%mband)
 type(pawcprj_type),allocatable :: Cp_v(:,:),Cp_c(:,:)

!************************************************************************

 call wrtout(std_out," Calculating optical matrix elements in the IBZ","COLL")
 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")

 comm = Wfd%comm
 my_rank = Wfd%my_rank

 nsppol  = Wfd%nsppol
 nspinor = Wfd%nspinor
 usepaw  = Wfd%usepaw

 if (usepaw==1) then
   ABI_MALLOC(Cp_v,(Wfd%natom,nspinor))
   call pawcprj_alloc(Cp_v,0,Wfd%nlmn_atm)
   ABI_MALLOC(Cp_c,(Wfd%natom,nspinor))
   call pawcprj_alloc(Cp_c,0,Wfd%nlmn_atm)
 end if

 if (inclvkb==1.and.usepaw==0) then
   MSG_ERROR("inclvkb==1 not coded,using inclvkb==2")
 end if
 !
 ! Calculate the matrix elements of ir in the IBZ.
 ABI_MALLOC(ir_kibz,(3,lomo_min:max_band,lomo_min:max_band,Wfd%nkibz,nsppol))
 ir_kibz=czero

 ABI_MALLOC(bbp_distrb, (Wfd%mband,Wfd%mband))

 do spin=1,nsppol
   do ik_ibz=1,Wfd%nkibz
    !
    ! Distribute the (b,b') entries.
    bbp_mask=.FALSE.; bbp_mask(lomo_spin(spin):max_band,lomo_spin(spin):max_band)=.TRUE.
    call wfd%distribute_bbp(ik_ibz,spin,"All",my_nbbp,bbp_distrb,bbp_mask=bbp_mask)
    if (ALL(bbp_distrb/=my_rank)) CYCLE

    istwf_k = Wfd%istwfk(ik_ibz)
    ABI_CHECK(istwf_k==1,"istwf_k/=1 not coded") ! KB stuff is missing.
    npw_k = Wfd%npwarr(ik_ibz)
    kg_k  => Wfd%Kdata(ik_ibz)%kg_k

    if (inclvkb/=0.and.usepaw==0) then
      ! Prepare term i <n,k|[Vnl,r]|n"k>
      call vkbr_init(vkbr,Cryst,Psps,inclvkb,istwf_k,npw_k,Kmesh%ibz(:,ik_ibz),kg_k)
    end if

    ! Note: spinorial case is not coded therefore we work with ihrc(:,1).
    ! TODO: The lower triangle can be Reconstructed by symmetry.
    do ib_v=lomo_spin(spin),max_band ! Loop over bands
      if ( ALL(bbp_distrb(ib_v,:)/=my_rank) ) CYCLE

      ABI_CHECK(wfd%get_wave_ptr(ib_v, ik_ibz, spin, wave_v, msg) == 0, msg)
      ug_v => wave_v%ug
      if (usepaw==1) call wfd%get_cprj(ib_v,ik_ibz,spin,Cryst,Cp_v,sorted=.FALSE.)

      do ib_c=lomo_spin(spin),max_band
       if (bbp_distrb(ib_v,ib_c)/=my_rank) CYCLE
       ABI_CHECK(wfd%get_wave_ptr(ib_c, ik_ibz, spin, wave_c, msg) == 0, msg)
       ug_c => wave_c%ug

       if (usepaw==0) then
         ! Calculate matrix elements of i[H,r] for NC pseudopotentials.
         ihrc = nc_ihr_comm(vkbr,cryst,psps,npw_k,nspinor,istwf_k,inclvkb,Kmesh%ibz(:,ik_ibz),ug_c,ug_v,kg_k)

       else
         ! Matrix elements of i[H,r] for PAW.
         call wfd%get_cprj(ib_c,ik_ibz,spin,Cryst,Cp_c,sorted=.FALSE.)

         ihrc = paw_ihr(spin,nspinor,npw_k,istwf_k,Kmesh%ibz(:,ik_ibz),Cryst,Pawtab,ug_c,ug_v,kg_k,Cp_c,Cp_v,HUr)
       end if
       !
       ! Save matrix elements of i*r in the IBZ
       ediff = KS_Bst%eig(ib_c,ik_ibz,spin) - KS_BSt%eig(ib_v,ik_ibz,spin)
       if (ABS(ediff)<tol16) ediff=tol6  ! Treat a possible degeneracy between v and c.
       ir_kibz(:,ib_c,ib_v,ik_ibz,spin) = ihrc(:,1)/ediff

      end do !ib_c
    end do !ib_v

    call vkbr_free(vkbr)
   end do !spin
 end do !ik_ibz

 ! Collect results on each node.
 call xmpi_sum(ir_kibz,comm,ierr)

 ABI_FREE(bbp_distrb)

 if (usepaw==1) then
   call pawcprj_free(Cp_v)
   ABI_FREE(Cp_v)
   call pawcprj_free(Cp_c)
   ABI_FREE(Cp_c)
 end if
 !
 ! ======================================================
 ! ==== Calculate Fcv(kBZ) in the full Brilouin zone ====
 ! ======================================================
 !
 ! Symmetrization of the matrix elements of the position operator.
 ! <Sk b|r|Sk b'> = R <k b|r|k b'> + \tau \delta_{bb'}
 !   where S is one of the symrec operations in reciprocal space, R is the
 !   corresponding operation in real space, \tau being the associated fractional translations.
 !
 ! q.Mcv( Sk) =  S^{-1}q. Mcv(k)
 ! q.Mcv(-Sk) = -S^{-1}q. CONJG(Mcv(k)) if time-reversal is used.

 b1=Cryst%gprimd(:,1)*two_pi
 b2=Cryst%gprimd(:,2)*two_pi
 b3=Cryst%gprimd(:,3)*two_pi

 opt_cvk = czero
 do spin=1,nsppol
   do ik_bz=1,nkbz
    !
    ! Get ik_ibz, and symmetries index from ik_bz.
    call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k)

    mat_dp = DBLE(Cryst%symrec(:,:,isym_k))
    call matrginv(mat_dp,3,3) ! Invert
    qrot = (3-2*itim_k) * MATMUL(mat_dp,qpoint)

    do ib_v=lomo_spin(spin),max_band !  Loops over the bands C and V start
      do ib_c=lomo_spin(spin),max_band
        !if (ib_c==ib_v) CYCLE
        emcvk = pdtqrc(qrot,ir_kibz(:,ib_c,ib_v,ik_ibz,spin),b1,b2,b3)
        if (itim_k==2) emcvk = CONJG(emcvk)
        opt_cvk(ib_c,ib_v,ik_bz,spin) = emcvk
      end do !ib_c
    end do !ib_v

   end do !ik_bz
 end do !spin

 ABI_FREE(ir_kibz)

 call xmpi_barrier(comm)

contains
!!***

!!****f* ABINIT/pdtqrc
!! NAME
!!  pdtqrc
!!
!! FUNCTION
!!  Calculate the dot product of a real vector with a complex vector, where each is in terms of b1-b3
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! SOURCE

pure function pdtqrc(R,C,b1,b2,b3)

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: R(3),b1(3),b2(3),b3(3)
 complex(dpc),intent(in) :: C(3)
 complex(dpc) :: pdtqrc

!Local variables ------------------------------
!scalars
 integer :: ii

!************************************************************************

 pdtqrc=czero
 do ii=1,3
   pdtqrc = pdtqrc + (R(1)*b1(ii)+R(2)*b2(ii)+R(3)*b3(ii)) * &
&                    (C(1)*b1(ii)+C(2)*b2(ii)+C(3)*b3(ii))
 end do

end function pdtqrc
!!***

end subroutine calc_optical_mels
!!***

end module m_wfd_optic
!!***
