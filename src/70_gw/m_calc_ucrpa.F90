!!****m* ABINIT/m_calc_ucrpa
!! NAME
!!  m_calc_ucrpa
!!
!! FUNCTION
!! Calculate the effective interaction in the correlated orbitals
!!
!! COPYRIGHT
!! Copyright (C) 2006-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

MODULE m_calc_ucrpa

#ifdef FC_INTEL
!DEC$ NOOPTIMIZE
#endif


  use defs_basis
 implicit none

 private

 public :: calc_ucrpa
!!***

contains
!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_ucrpa
!! NAME
!! calc_ucrpa
!!
!! FUNCTION
!! Calculate the effective interaction in the correlated orbitals
!!
!! COPYRIGHT
!! Copyright (C) 1999-2019 ABINIT group (TApplencourt,BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! npwe : number of plane wave for the dielectric constant
!! npw : number of plane wave
!! nomega  : number of frequencis
!! bandinf,bandsup : kohn sham band
!! optimisation : string for the optimisation
!! Wfd:: MPI communicator
!! mesh <kmesh_t>
!!    %nbz=Number of points in the BZ
!!    %nibz=Number of points in IBZ
!!    %kibz(,nibz)=k-point coordinates, irreducible Brillouin zone
!!    %kbz(3,nbz)=k-point coordinates, full Brillouin zone
!!    %ktab(nbz)= table giving for each k-point in the BZ (kBZ), the corresponding
!!    %ktabi(nbz)= for each k-point in the BZ defines whether inversion has to be considered
!!    %ktabp(nbz)= phase factor associated to tnons
!! M1_q_m(bandinf:bandsup,bandinf:bandsup,npw,Qmesh%nibz): Oscillator strengh in Wannier basis
!! rhot1_q_m(bandinf:bandsup,bandinf:bandsup,npw,Qmesh%nibz): Oscillator strengh
!!                                          multiplied by coulomb potential in Wannier basis
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      affichage,checkk,cpu_time,get_bz_item,read_screening,wrtout
!!      xmpi_barrier,xmpi_sum
!!
!! SOURCE

 subroutine calc_ucrpa(itypatcor,cryst,Kmesh,lpawu,M1_q_m,Qmesh,npwe,&
& npw,nsym,nomega,omegamin,omegamax,bandinf,bandsup,optimisation,ucvol,Wfd,fname,plowan_compute,rhot1,wanbz)

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors

 use m_io_tools,      only : open_file
 use m_wfd,           only : wfd_t
 use m_io_screening,  only : read_screening, em1_ncname
 use m_bz_mesh,       only : kmesh_t, get_BZ_item
 use m_crystal,       only : crystal_t
 use m_plowannier,    only : operwan_realspace_type,plowannier_type
 implicit none
!   _____                   _
!  |_   _|                 | |
!    | |  _ __  _ __  _   _| |_
!    | | | '_ \| '_ \| | | | __|
!   _| |_| | | | |_) | |_| | |_
!  |_____|_| |_| .__/ \__,_|\__|
!              | |
!              |_|

!Arguments ------------------------------------
 integer, intent(in)   :: itypatcor,lpawu,npw,npwe,nsym
 integer, intent(in)   :: nomega
 integer, intent(in)   :: bandinf
 integer, intent(in)   :: bandsup
 integer, intent(in)   :: plowan_compute 
 character(len=fnlen), intent(in) :: fname
 character(len=*), intent(in) :: optimisation
 real(dp), intent(in) :: ucvol,omegamin,omegamax

 type(wfd_t),intent(inout) :: Wfd
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(crystal_t),intent(in) :: Cryst
 type(operwan_realspace_type),intent(in) :: rhot1(npw,Qmesh%nibz)
 type(plowannier_type),intent(in) :: wanbz
 complex(dpc), intent(in) :: M1_q_m(cryst%nattyp(itypatcor),Wfd%nspinor,Wfd%nspinor,2*lpawu+1,2*lpawu+1,npw,Qmesh%nibz)

!Local variables ------------------------------
!scalars
 real(dp) :: x
 real(dp) :: t1,t2
 real(dp):: tol
 complex(dpc) :: uu,jj

 complex :: nC,ualter

 integer :: iatom1,iatom2,iatom3,iatom4,pos1,pos2,pos3,pos4,il1,il2,il3,il4,ispin,one_orbital
 integer :: im_paral,iqalloc,ib1,ib2,m1,m2,m3,m4,iqibz,mbband1,mbband2,mbband3,mbband4,spin1,spin2
 integer :: ierr,ik_bz,ik_ibz,iq_ibz,i,iG1,iG2,iG,iiG,iomega,iomega1,ispinor1,ispinor2,ispinor3,ispinor4
 integer :: lpawu_read,nkibz,nbband,nkbz,nprocs,nqalloc,nqibz,ms1,ms2,ms3,ms4,mbband,nspinor
 integer :: isym_kgw,iik,unt, cp_paral
 complex(dpc) ::ph_mkt

 logical :: wannier=.TRUE.
 logical :: verbose=.FALSE.
 logical :: bug=.FALSE.
 logical :: lscr_one

 character(len=500) :: message

!arrays
 complex(dpc), allocatable :: V_m(:,:,:,:)
 complex(dpc), allocatable :: U_m(:,:,:,:)
 complex(dpc),allocatable :: uspin(:,:),jspin(:,:)
! complex(dpc), allocatable :: coeffW_BZ(:,:,:),coeffW_IBZ(:,:,:)
 complex(dpc), allocatable :: rhot_q_m1m3(:,:,:,:,:,:),rhot_q_m2m4(:,:,:,:,:,:)
 complex(dpc), allocatable :: rhot_q_m1m3_npwe(:,:,:,:,:,:),rhot_q_m2m4_npwe(:,:,:,:,:,:)
 complex(dpc),allocatable :: trrho(:,:),sumrhorhoeps(:)
 complex(gwpc), allocatable :: scr(:,:,:,:)

 real(dp),allocatable :: k_coord(:,:)!,k_coordIBZ(:,:)
 real(dp),allocatable :: q_coord(:,:)
 real(dp),allocatable:: normG(:)
 complex(dpc),allocatable:: uomega(:),jomega(:)
 real(dp),allocatable:: omega(:)

 integer,allocatable:: ikmq_bz_t(:,:)

 logical,allocatable :: bijection(:)
!************************************************************************

 write(message,*) ch10, '==== Calculation of the screened interaction ===='
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 write(message,*) ""
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 nkbz = Kmesh%nbz
 nqibz= Qmesh%nibz
 nspinor=Wfd%nspinor

 nbband=1+bandsup-bandinf
!  _  __            ____
! | |/ /   ___     / __ \
! | ' /   ( _ )   | |  | |
! |  <    / _ \/\ | |  | |
! | . \  | (_>  < | |__| |
! |_|\_\  \___/\/  \___\_\

 write(message,*) "Read K and Q mesh"
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 ABI_ALLOCATE(k_coord,(nkbz,3))
 ABI_ALLOCATE(q_coord,(nqibz,4))


!==Read k and q==!
!open(unit=2012,file='ikbz_COORD',form='formatted',status='unknown')
!read(2012,*) (ik_bz,k_coord(ik_bz,:),i=1,nkbz)
!close(2012)

 do ik_bz=1,nkbz
   call get_BZ_item(Kmesh,ik_bz,k_coord(ik_bz,:),ik_ibz,isym_kgw,iik,ph_mkt)
 end do

! open(unit=2012,file='iqbz_COORD',form='formatted',status='unknown')
! read(2012,*)
 do i=1,nqibz
!   read(2012,*) iq_ibz,q_coord(iq_ibz,:)
   q_coord(i,1)=Qmesh%ibz(1,i)
   q_coord(i,2)=Qmesh%ibz(2,i)
   q_coord(i,3)=Qmesh%ibz(3,i)
   q_coord(i,4)=Qmesh%wt(i)
!   if (iq_ibz > nqibz) then
!     write(message,*) iq_ibz,nqibz," Error on line",i,"Are you in iBZ ?"
!     call wrtout(std_out,message,'COLL')
!   end if
 end do
! close(2012)

!==Bijection and array for k-q==!
 ABI_ALLOCATE(bijection,(nkbz))
 ABI_ALLOCATE(ikmq_bz_t,(nkbz,nqibz))
 bijection(:)=.FALSE.
 if (nsym==1) then
 do ik_bz=1,nkbz
   do iq_ibz=1,nqibz
     ikmq_bz_t(ik_bz,iq_ibz)=findkmq(ik_bz,k_coord,q_coord(iq_ibz,:),nkbz)
     if (ikmq_bz_t(ik_bz,iq_ibz)>nkbz.and.nsym==1) then
       BUG=.TRUE.
       write(message,*) "No K-Q for K/Q =",ik_bz,iq_ibz
       MSG_ERROR(message)
     end if
     bijection(ikmq_bz_t(ik_bz,iq_ibz))=.TRUE.
   end do

   if (count(bijection).NE.nqibz.and.nsym==1) then
   BUG=.TRUE.
   write(message,*) 'No bijection ',ik_bz
   MSG_ERROR(message)
   end if

   bijection(:)=.FALSE.
 end do

 if (.NOT.BUG.and.nsym==1) then
   write(message,*)  "Bijection Ok."
   call wrtout(std_out,message,'COLL')
 end if
endif
!                                           _____
!                                          / ____|
!   _ __   ___  _ __ _ __ ___   ___       | |  __
!  | '_ \ / _ \| '__| '_ ` _ \ / _ \      | | |_ |
!  | | | | (_) | |  | | | | | |  __/      | |__| |
!  |_| |_|\___/|_|  |_| |_| |_|\___|       \_____|

 ABI_ALLOCATE(normG,(npw))
 if (verbose) then
   write(message,*) 'Read the potential and G norm'
   call wrtout(std_out,message,'COLL')
   if (open_file('normeG',message,newunit=unt,form='formatted',status='unknown') /= 0) then
     MSG_ERROR(message)
   end if
   read(unt,*) (iiG,x,normG(iiG),iG=1,npw)
   close(unt)
   !!False norme for G=0 idd G is the inverse of the potential inverse du potentiel (q=0)
   normG(1)=0
 end if

!========================================================================
!------------------------------------------------------------------------

!  FIRST PART OF THE ROUTINE: USE M_G^(nn')(q,k) to do formal checks.
!                             USE rhot_q_n to do compute bare interaction

!------------------------------------------------------------------------
!========================================================================

!                               _
!                              ( )
!   _ __ ___        _ __  _ __ |/
!  | '_ ` _ \      | '_ \| '_ \
!  | | | | | |     | | | | | | |
!  |_| |_| |_|     |_| |_|_| |_|

!==========================================================
!==========================================================
! tol=1E-1
! tolerance for the normalization of wfc: should be around 0.01.
 tol = 1 ! very large for test.
 write(message,*) 'Check the norm of M'
 call wrtout(std_out,message,'COLL')
 write(message,*) 'Tolerance :',tol
 call wrtout(std_out,message,'COLL')
!  __      __
!  \ \    / /
!   \ \  / /      _ __
!    \ \/ /      | '_ \
!     \  /       | | | |
!      \/        |_| |_|
!==========================================================================
!==Compute V_{n,n'}: bare interaction in the KS basis <- rhot_q_n -> V_n
!==========================================================================
!==========================================================
!==Compute V_{n,n'}
!==========================================================
 if(verbose) then
   write(message,*) ""
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   write(message,*)  "==Calcul of the bare kohn-sham interaction V n=="
   call wrtout(std_out,message,'COLL')
 endif
 tol=1E+1

 if (.NOT.wannier) RETURN

!========================================================================
!------------------------------------------------------------------------

!  SECOND PART OF THE ROUTINE: Read Wannier projections

!------------------------------------------------------------------------
!========================================================================
!
! \ \        / /       (_)
!  \ \  /\  / /_ _ _ __  _ __  _  ___ _ __
!   \ \/  \/ / _` | '_ \| '_ \| |/ _ \ '__|
!    \  /\  / (_| | | | | | | | |  __/ |
!     \/  \/ \__,_|_| |_|_| |_|_|\___|_|

!==========================================================
!==========================================================
!== Read Wannier coefficient in forlb.ovlp
!==========================================================
 nkibz=Kmesh%nibz

!Read "l"
 if (plowan_compute<10) then
   write(message,*) ""
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   write(message,*) "Read wannier in iBZ"
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   if (open_file('forlb.ovlp',message,newunit=unt,form='formatted',status='unknown') /= 0) then
     MSG_ERROR(message)
   end if
   rewind(unt)
 read(unt,*) message
   read(unt,*) message, lpawu_read
   read(unt,*) message, ib1, ib2
   close(unt)
   mbband=2*lpawu_read+1
 else 
   ib1=wanbz%bandi_wan
   ib2=wanbz%bandf_wan
   mbband=2*wanbz%latom_wan(1)%lcalc(1)+1
   write(message,*)"Read l and bands from wanbz",ib1,ib2,mbband
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 endif
 if(ib1/=bandinf.and.ib2/=bandsup) then
   write(message,*) "Error with bands",ib1,bandinf,ib2,bandsup
   MSG_ERROR(message)
 endif
!!Read the bandinf, bandinf redondance information
  
  
!*******************************************************
!if (3==4) then
!!  USELESS START
!!*******************************************************
! ABI_ALLOCATE(coeffW_IBZ,(bandinf:bandsup,nkibz,mbband))
! coeffW_IBZ=czero
!
!
! if(Wfd%my_rank==0) then
!   do ik_ibz=1,nkibz
!     !read k
!      read(2012,*)
!      do iband=bandinf,bandsup
!     !read band
!        read(2012,*)
!     !read projection
!        do m1=1,mbband
!          read(2012,*) binR,binR,binR,x,y
!          !write(message,*)  binR,binR,binR,x,y
!          coeffW_IBZ(iband,ik_ibz,m1)=cmplx(x,y)
!        end do
!      end do
!   end do
! endif
! call xmpi_barrier(Wfd%comm)
! call xcast_mpi(coeffW_IBZ,0,Wfd%comm,ierr)
! call xmpi_barrier(Wfd%comm)
! close(2012)
!
! ABI_ALLOCATE(coeffW_BZ,(bandinf:bandsup,nkbz,mbband))
!
! if (nkbz==nkibz) then
!   coeffW_BZ=coeffW_IBZ
! else
!   write(message,*) "Reconstruct in full BZ"
!   call wrtout(std_out,message,'COLL')
!   call wrtout(ab_out,message,'COLL')
!   ABI_ALLOCATE(k_coordIBZ,(nkibz,3))
!
!!   k_coordIBZ(:,:)=q_coord(:,1:3)
!
!   bijection(:)=.FALSE.
!   write(message,*) "Indice in iBZ | Indice in BZ | Inverse in BZ"
!   call wrtout(std_out,message,'COLL')
!   do ik_bz=1,Kmesh%nbz
!     write(6,*) "ik",ik_bz,Kmesh%tab(ik_bz),Kmesh%tabi(ik_bz),Kmesh%tabo(ik_bz)
!     if(Kmesh%tabi(ik_bz)==1) then
!       coeffW_BZ(:,ik_bz,:)=coeffW_IBZ(:,Kmesh%tab(ik_bz),:)
!       bijection(ik_bz)=.TRUE.
!     else if(Kmesh%tabi(ik_bz)==-1) then
!       coeffW_BZ(:,ik_bz,:)=conjg(coeffW_IBZ(:,Kmesh%tab(ik_bz),:))
!       write(message,*) Kmesh%tab(ik_bz),ik_bz
!       inverse_ik_bz=Kmesh%tab(ik_bz)
!       bijection(ik_bz)=.TRUE.
!!       bijection(inverse_ik_bz)=.TRUE.
!     endif
!   enddo
!
!   if (count(bijection).NE.nkbz) then
!    BUG=.TRUE.
!    write(message,*) 'Miss somme K point for the Wannier',count(bijection),"/",nkbz
!    MSG_ERROR(message)
!   end if
!
!   if (.NOT.BUG) then
!     write(message,*) "Reconstruction Success"
!     call wrtout(std_out,message,'COLL')
!     call wrtout(ab_out,message,'COLL')
!   end if
!   ABI_DEALLOCATE(k_coordIBZ)
! end if
!
! ABI_DEALLOCATE(coeffW_IBZ)
!
! wk=1.0/nkbz
!
! write(message,*) 'Orthogonality check'
! call wrtout(std_out,message,'COLL')
! write(message,*)  'Sum on all the k point ,on all the Kohn-Sham band of C_(m1)*C_m(2)'
! call wrtout(std_out,message,'COLL')
!
!! tolerance for the sum over k-points of the Wannier functions (orthogonality)
! tol=1E-5
!
!! Sum for one k-point (should be around 0.1).
! tol2=10E0
!
! write(message,*) 'Tolerance : k',tol2,'m',tol
! call wrtout(std_out,message,'COLL')
!
!
! ! Here checks on Wannier coeff are done (ortho, norm..)
! nC=cmplx(0,0)
! BUG=.FALSE.
!
! do m1=1,mbband
!   do m2=1,mbband
!     do ik_bz=1,nkbz
!       nCt=sum(conjg(coeffW_BZ(:,ik_bz,m1))*coeffW_BZ(:,ik_bz,m2))
!       if (  ((m1==m2).and.(abs(abs(ncT)-1)>tol2)).OR.&
!       ((m1.NE.m2).and.(abs(ncT)>tol2)) ) then
!         BUG=.TRUE.
!         write(message,*)  "No orthogonality for m1,m2",m1,m2,"kpt",ik_bz,abs(nCt)
!         MSG_ERROR(message)
!       end if
!       nC=nC+wk*nCt
!     end do
!
!     if (  ((m1==m2).and.(abs(abs(nC)-1)>tol)).OR.&
!     ((m1.NE.m2).and.(abs(nC)>tol)) ) then
!       bug=.TRUE.
!       write(message,*) "No orthogonality for",m1,m2,abs(nC)
!       MSG_ERROR(message)
!     end if
!     write(message,*)  m1,m2,abs(nC)
!     nC=cmplx(0,0)
!   end do
! end do
! if (.NOT.bug) then
!   write(message,*) "Orthogonality check"
!   call wrtout(std_out,message,'COLL')
! end if
!
!!*******************************************************
!end if
!  USELESS END
!*******************************************************

! do iq_ibz=1,nqibz
!   do m1=1,mbband
!     do m2=1,mbband
!       write(6,*) "M1_q_m",M1_q_m(m1,m2,1,iq_ibz)
!     end do
!   end do
! end do

 if(real(M1_q_m(1,1,1,1,1,1,1))>0) then
 endif

!                                       _
!                                      ( )
!   _ __ ___        _ __ ___  _ __ ___ |/
!  | '_ ` _ \      | '_ ` _ \| '_ ` _ \
!  | | | | | |     | | | | | | | | | | |
!  |_| |_| |_|     |_| |_| |_|_| |_| |_|


!!!==Calculation of M_G^(mm')(q)==!

 write(message,*) 'Calculation of M  m'
 call wrtout(std_out,message,'COLL')

!=================================================================!
!==Compute V_{m,m'}(q,z): Oscillator strengh in the Wannier basis
!=================================================================!

! Sum over k-points for the oscillator strengh in the Wannier basis

 if (verbose) then
!   call  Sauvegarde_M_q_m(M1_q_m,normG,nqibz,npw,mbband)
 end if

! write(message,*)  "M1_q_m,m for iG=53 and iq=8"
!    call wrtout(std_out,message,'COLL')

!        _                                     _
!       | |                                   ( )
!   _ __| |__   ___        _ __ ___  _ __ ___ |/
!  | '__| '_ \ / _ \      | '_ ` _ \| '_ ` _ \
!  | |  | | | | (_) |     | | | | | | | | | | |
!  |_|  |_| |_|\___/      |_| |_| |_|_| |_| |_|

 write(message,*) ""
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
!==Calcul de Ro_G^(mm')(q)==!
 write(message,*) 'Calculation of rhotwilde  m'
 call wrtout(std_out,message,'COLL')

!===============================================================!
!===============================================================!
!==Compute V_{m,m'}(q,z): bare interaction in the Wannier basis
!===============================================================!
!===============================================================!
 ABI_DEALLOCATE(ikmq_bz_t)

 !ABI_ALLOCATE(rhot_q_m,(nspinor,nspinor,mbband,mbband,npw,nqibz))
 
 ! if (plowan_compute>=10)then
 !   write(message,*)" cRPA calculation using plowan module"
 !   call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
 !   do iG=1,npw
 !     do iqibz=1,nqibz
 !       do ispinor1=1,nspinor
 !         do ispinor2=1,nspinor
 !           do m1=1,2*wanbz%latom_wan(1)%lcalc(1)+1
 !             do m2=1,2*wanbz%latom_wan(1)%lcalc(1)+1
 !               rhot_q_m(ispinor1,ispinor2,m1,m2,iG,iqibz)=&
 !                 &rhot1(iG,iqibz)%atom_index(1,1)%position(1,1)%atom(1,1)%matl(m1,m2,1,ispinor1,ispinor2)
 !               !write(6,*) rhot_q_m(ispinor1,ispinor2,m1,m2,iG,iqbz)
 !             enddo
 !           enddo
 !         enddo
 !       enddo
 !     enddo
 !   enddo
 ! endif
  do iatom1=1,wanbz%natom_wan
  do iatom2=1,wanbz%natom_wan
  do iatom3=1,wanbz%natom_wan
  do iatom4=1,wanbz%natom_wan
    if (iatom1/=iatom3 .or. iatom2/=iatom4)cycle
    do pos1=1,size(wanbz%nposition(iatom1)%pos,1)
    do pos2=1,size(wanbz%nposition(iatom2)%pos,1)
    do pos3=1,size(wanbz%nposition(iatom3)%pos,1)
    do pos4=1,size(wanbz%nposition(iatom4)%pos,1)
      if (pos1/=pos3 .or. pos2/=pos4)cycle
      do il1=1,wanbz%nbl_atom_wan(iatom1)
      do il2=1,wanbz%nbl_atom_wan(iatom2)
      do il3=1,wanbz%nbl_atom_wan(iatom3)
      do il4=1,wanbz%nbl_atom_wan(iatom4)
      if(il1/=il3 .or. il2/=il4)cycle
      if (wanbz%nsppol/=1)then
        ABI_ALLOCATE(uspin,(4,nomega))
        ABI_ALLOCATE(jspin,(4,nomega))
      endif
      if (iatom1==iatom2 .and. pos1==pos2 .and. il1 == il2)then
        one_orbital=1
      else
        one_orbital=0
      endif
      ABI_ALLOCATE(omega,(nomega))
      do spin1=1,wanbz%nsppol
      do spin2=1,wanbz%nsppol
        cp_paral=0
        mbband1=2*wanbz%latom_wan(iatom1)%lcalc(il1)+1
        mbband2=2*wanbz%latom_wan(iatom2)%lcalc(il2)+1
        mbband3=2*wanbz%latom_wan(iatom3)%lcalc(il3)+1
        mbband4=2*wanbz%latom_wan(iatom4)%lcalc(il4)+1
        ABI_ALLOCATE(rhot_q_m1m3,(nspinor,nspinor,mbband1,mbband3,npw,nqibz))
        ABI_ALLOCATE(rhot_q_m2m4,(nspinor,nspinor,mbband2,mbband4,npw,nqibz))
        do ispinor1=1,wanbz%nspinor
        do ispinor2=1,wanbz%nspinor
        do ispinor3=1,wanbz%nspinor
        do ispinor4=1,wanbz%nspinor
          do m1=1,mbband1
          do m2=1,mbband2
          do m3=1,mbband3
          do m4=1,mbband4
            do iqibz=1,nqibz
            do iG=1,npw
              ! cp_paral=cp_paral+1
              ! if(mod(cp_paral-1,nprocs)==Wfd%my_rank) then
                rhot_q_m1m3(ispinor1,ispinor3,m1,m3,iG,iqibz)=&
                  &rhot1(iG,iqibz)%atom_index(iatom1,iatom3)%position(pos1,pos3)%atom(il1,il3)%matl(m1,m3,spin1,ispinor1,ispinor3)
          
                rhot_q_m2m4(ispinor2,ispinor4,m2,m4,iG,iqibz)=&
                  &rhot1(iG,iqibz)%atom_index(iatom2,iatom4)%position(pos2,pos4)%atom(il2,il4)%matl(m2,m4,spin2,ispinor2,ispinor4)
              !endif
            enddo!iG
            enddo!iqibz
          enddo!m4
          enddo!m3
          enddo!m2
          enddo!m1
        enddo!ispinor4
        enddo!ispinor3
        enddo!ispinor2
        enddo!ispinor1 
        ! call xmpi_barrier(Wfd%comm)  ! First make sure that all processors are here
        ! call xmpi_sum(rhot_q_m1m3,Wfd%comm,ierr)
        ! call xmpi_sum(rhot_q_m2m4,Wfd%comm,ierr)
        ! call xmpi_barrier(Wfd%comm)  ! First make sure that all processors are here

!    __      __
!    \ \    / /
!     \ \  / /      _ __ ___
!      \ \/ /      | '_ ` _ \
!       \  /       | | | | | |
!        \/        |_| |_| |_|

        
        ABI_ALLOCATE(V_m,(mbband1*nspinor,mbband2*nspinor,mbband3*nspinor,mbband4*nspinor))
        write(message,*) "";call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
        write(message,*)  "==Calculation of the bare interaction  V m=="
        call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
        
        call print_orbitals(spin1,spin2,iatom1,iatom2,iatom3,iatom4,pos1,pos2,pos3,pos4,il1,il2,il3,il4,wanbz,1)
        V_m=czero
        call cpu_time ( t1 )
        im_paral=0
        nprocs  = xmpi_comm_size(Wfd%comm)
   
   




        do ispinor1=1,nspinor
        do ispinor2=1,nspinor
        do ispinor3=1,nspinor
        do ispinor4=1,nspinor
          do m1=1,mbband1
          do m2=1,mbband2
          do m3=1,mbband3
          do m4=1,mbband4
            ms1=m1+(ispinor1-1)*mbband1
            ms2=m2+(ispinor2-1)*mbband2
            ms3=m3+(ispinor3-1)*mbband3
            ms4=m4+(ispinor4-1)*mbband4
            im_paral=im_paral+1
            if(mod(im_paral-1,nprocs)==Wfd%my_rank) then
!!somme interne sur iG, puis somme externe sur iq_ibz
!! Sum_(iq_ibz) wi(iq_ibz)*Sum_ig  Rho(m3,m1,iG,iq)cong*Rho(m2,m4,ig,iq)
              V_m(ms1,ms2,ms3,ms4)=sum(q_coord(:,4)*sum(conjg(rhot_q_m1m3(ispinor3,ispinor1,m3,m1,:,:))* &
              &rhot_q_m2m4(ispinor2,ispinor4,m2,m4,:,:),dim = 1))*Ha_eV/(ucvol)
            endif!paral
          end do!m4
          end do!m3
          end do!m2
          end do!m1
        end do!ispinor4
        end do!ispinor3
        end do!ispinor2
        end do!ispinor1
        call xmpi_barrier(Wfd%comm)  ! First make sure that all processors are here
        call xmpi_sum(V_m,Wfd%comm,ierr)
        call xmpi_barrier(Wfd%comm)  ! First make sure that all processors are here
        call cpu_time ( t2 )
        write(message,*)  "in ",t2-t1,"sec"
        call wrtout(std_out,message,'COLL')
!  !==Check if calculation is correct
        tol=1E-2
      
        write(message,*)  "BARE INTERACTION"
        call wrtout(std_out,message,'COLL')
        call checkk(V_m,1,mbband*nspinor,1,0,uu,jj,"bare interaction",mbband1,mbband2,mbband3,mbband4,nspinor,one_orbital)
        !call print_U(mbband1,mbband2,mbband3,mbband4,nspinor,V_m)
!  ========================================================================
!  ------------------------------------------------------------------------

!    THIRD PART OF THE ROUTINE: Read dielectric function

!  ------------------------------------------------------------------------
!  ========================================================================

!      _____                          _
!     / ____|                        (_)
!    | (___   ___ _ __ ___  ___ _ __  _ _ __   __ _
!     \___ \ / __| '__/ _ \/ _ \ '_ \| | '_ \ / _` |
!     ____) | (__| | |  __/  __/ | | | | | | | (_| |
!    |_____/ \___|_|  \___|\___|_| |_|_|_| |_|\__, |
!                                              __/ |
!                                             |___/
!  ==========================================================
!  == Read Dielectric Matrix for _SCR file
!  ==========================================================
        write(message,*) "";call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
        call wrtout(std_out,message,'COLL')
        call wrtout(ab_out,message,'COLL')
 !  write(message,*) "==Read the dielectric matrix=="
 !  call wrtout(ab_out,message,'COLL'); call wrtout(std_out,message,'COLL')

        lscr_one=.true.

   ! if symetry is not enabled then a large number of q-point are used
   ! and they are read with direct access to avoid having too large
   ! memory.
   !-------------------------------------------------------------------
        if(nsym>1.and..not.lscr_one) nqalloc=nqibz
        if(nsym==1.or.lscr_one) nqalloc=1

        ABI_ALLOCATE(scr,(npwe,npwe,nomega,nqalloc))
        scr=czero
        if(nsym>1.and..not.lscr_one) then
          write(message,*) "==Read the dielectric matrix=="
          call wrtout(ab_out,message,'COLL'); call wrtout(std_out,message,'COLL')
          call read_screening(em1_ncname,fname,npwe,nqibz,nomega,scr,IO_MODE_MPI,Wfd%comm)
        endif

!   if (verbose) then
!     open(unit=2211,file='Screening',form='formatted',status='unknown')
!     do iG1=1,npwe
!       do iG2=1,npwe
!         write(2211,*) iG1,iG2,normG(iG1),normG(iG2),abs(scr(iG1,iG2,1,1)),abs(scr(iG1,iG2,2,1))
!       end do
!     end do
!     close(2211)
!   end if
!
        if(nsym>1.and..not.lscr_one) then
          write(message,*) "Check the hermiticity"
          call wrtout(ab_out,message,'COLL')
          call wrtout(std_out,message,'COLL')
          tol = 1E-2
          do iG1=1,npwe
            if (modulo(iG1,100).EQ.1) then
              write(message,*)  iG1,"/",npw
              call wrtout(std_out,message,'COLL')
            end if
            do iG2=iG1,npwe
              if (ANY(abs(scr(iG1,iG2,:,:)-scr(iG2,iG1,:,:))>tol)) then
                write(message,*) iG1,iG2,"False"
                MSG_ERROR(message)
                do iomega1=1,nomega
                  if(abs(scr(iG1,iG2,iomega1,1)-scr(iG2,iG1,iomega1,1))>tol) then
                    write(message,*) iG1,iG2,"False",scr(iG1,iG2,iomega1,1),scr(iG2,iG1,iomega1,1)
                    call wrtout(std_out,message,'COLL')
                  endif
                enddo
              end if
              if(iG1==iG2) then
                scr(iG1,iG1,:,:)=scr(iG1,iG1,:,:)-one ! unscreened part of
              endif
            end do
          end do
          write(message,*)  "Done: Hermiticity of dielectric matrix checked"
          call wrtout(std_out,message,'COLL')
          call wrtout(ab_out,message,'COLL')
        endif

!  ========================================================================
!  ------------------------------------------------------------------------

!    FOURTH PART OF THE ROUTINE: Use dielectric function and oscillator
!    strengh in Wannier basis to compute Screened cRPA interactions.

!  ------------------------------------------------------------------------
!  ========================================================================
!     _    _
!    | |  | |
!    | |  | |      _ __ ___
!    | |  | |     | '_ ` _ \
!    | |__| |     | | | | | |
!     \____/      |_| |_| |_|
!
        write(message,*) "";call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
        write(message,*)  "== Calculation of the screened interaction on the correlated orbital U m =="
        call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
        write(message,*)ch10,  " = Start loop over frequency "
        call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
        
        ABI_ALLOCATE(U_m,(mbband1*nspinor,mbband2*nspinor,mbband3*nspinor,mbband4*nspinor))
        ABI_ALLOCATE(rhot_q_m1m3_npwe,(nspinor,nspinor,mbband1,mbband3,npwe,nqibz))
        ABI_ALLOCATE(rhot_q_m2m4_npwe,(nspinor,nspinor,mbband2,mbband4,npwe,nqibz))
        ABI_ALLOCATE(trrho,(npwe,nqibz))
        ABI_ALLOCATE(sumrhorhoeps,(nqibz))

   ! Following lines are only for debug
   !--------------------------------------------
        trrho(:,:)=czero
        do iG1=1,npwe
          rhot_q_m1m3_npwe(:,:,:,:,iG1,:)=rhot_q_m1m3(:,:,:,:,iG1,:)
          rhot_q_m2m4_npwe(:,:,:,:,iG1,:)=rhot_q_m2m4(:,:,:,:,iG1,:)
        enddo
        !  do  ispinor1=1,nspinor
        !     do m1=1,mbband
        !       trrho(iG1,:)= trrho(iG1,:)+rhot_q_m(ispinor1,ispinor1,m1,m1,iG1,:)
        !     enddo
        !   enddo
        ! enddo
!   write(6,*) "trrho"
!   do iG1=1,npwe
!     do iq_ibz=1,nqibz
!        write(6,*) iG1,iq_ibz,trrho(iG1,iq_ibz)
!     enddo
!   enddo
!   write(6,*)

   ! Loop over frequencies to compute cRPA U(w)
   !--------------------------------------------
        ABI_ALLOCATE(uomega,(nomega))
        ABI_ALLOCATE(jomega,(nomega))
        
        iomega=1
        do iomega=1,nomega
          write(message,'(2a,i4,2a)')ch10,  " --- For frequency w =",iomega, "  -------------",ch10
          call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
          sumrhorhoeps(:)=czero
          ualter=czero
          
          write(std_out,*) Optimisation
          U_m=cmplx(0,0)
          SELECT CASE(trim(Optimisation))
          CASE("naif")
          CASE("onlyG")
          CASE("Gsum2")
          CASE("Gsum")
! write(message,*)  "Optimisation on G and sum"
! call wrtout(std_out,message,'COLL')
            U_m=cmplx(0,0)
            nc=cmplx(0,0)
            im_paral=0
            call cpu_time ( t1 )
            nprocs  = xmpi_comm_size(Wfd%comm)
            do iq_ibz=1,nqibz
              
              if(nsym>1.and..not.lscr_one)  iqalloc=iq_ibz
              if(nsym==1.or.lscr_one) iqalloc=1
              if(nsym==1.or.lscr_one) then
#ifdef HAVE_MPI_IO
                call read_screening(em1_ncname,fname,npwe,1,nomega,scr,IO_MODE_MPI,Wfd%comm,iqiA=iq_ibz)
#else
                call read_screening(em1_ncname,fname,npwe,1,nomega,scr,IO_MODE_FORTRAN,Wfd%comm,iqiA=iq_ibz)
#endif
                write(message,*) "Check the hermiticity"
                call wrtout(std_out,message,'COLL')
                tol = 0.01_dp
                do iG1=1,npwe
                  if (modulo(iG1,100).EQ.1) then
                    write(message,*)  iG1,"/",npw
                    call wrtout(std_out,message,'COLL')
                  end if
                  do iG2=iG1,npwe
                    if (ANY(abs(scr(iG1,iG2,:,:)-scr(iG2,iG1,:,:))>tol)) then
                      do iomega1=1,nomega
                        if(abs(scr(iG1,iG2,iomega1,1)-scr(iG2,iG1,iomega1,1))>tol) then
                          write(message,*) iG1,iG2,"False",scr(iG1,iG2,iomega1,1),scr(iG2,iG1,iomega1,1)
                          call wrtout(std_out,message,'COLL')
                        endif
                      enddo
                      write(message,*) "CHECK THE HERMITICITY"
                      MSG_WARNING(message)
                    end if
                    if(iG1==iG2) then
                      scr(iG1,iG1,:,:)=scr(iG1,iG1,:,:)-one ! unscreened part of
!                 interaction is computed before
!                 scr(iG1,iG1,:,:)=one
!               else
!                 scr(iG1,iG2,:,:)=zero
                    endif
                  end do
                end do
              endif
              write(message,*)  "Done: Hermiticity of dielectric matrix checked"
              
              ! do iG1=1,npwe
!                 do iG2=1,npwe
!                   sumrhorhoeps(iq_ibz)=sumrhorhoeps(iq_ibz)+conjg(trrho(iG1,iq_ibz))*trrho(iG2,iq_ibz)*scr(iG2,iG1,iomega,1)
!                   sumrhorhoeps(iq_ibz)=sumrhorhoeps(iq_ibz) + scr(iG2,iG1,iomega,iqalloc)
!                 enddo
!               enddo
!             ualter=ualter+sumrhorhoeps(iq_ibz)*q_coord(iq_ibz,4)*Ha_eV/ucvol/(mbband)**2
!         write(6,*) "sumrhorhoeps",iq_ibz,sumrhorhoeps(iq_ibz)
!         write(6,*) "ualter",iq_ibz,ualter
              do ispinor1=1,nspinor
              do ispinor2=1,nspinor
              do ispinor3=1,nspinor
              do ispinor4=1,nspinor
                do m1=1,mbband1
                do m2=1,mbband2
                do m3=1,mbband3
                do m4=1,mbband4
                  ms1=m1+(ispinor1-1)*mbband1
                  ms2=m2+(ispinor2-1)*mbband2
                  ms3=m3+(ispinor3-1)*mbband3
                  ms4=m4+(ispinor4-1)*mbband4
                  im_paral=im_paral+1
                  if(mod(im_paral-1,nprocs)==Wfd%my_rank) then

!     sum q sum G1 sum G2 f(G1,q)f(G2,q)G(G1,G2,q)
!     sum q sum G1 f(g1,q) sum G2 f(G2,q)G(G1,G2,q)
                    do iG1=1,npwe
                      nc=nc+conjg(rhot_q_m1m3_npwe(ispinor1,ispinor3,m1,m3,iG1,iq_ibz))*&
                      &sum(rhot_q_m2m4_npwe(ispinor2,ispinor4,m2,m4,:,iq_ibz)*scr(:,iG1,iomega,iqalloc))
                    end do
                    
                    U_m(ms1,ms2,ms3,ms4)=U_m(ms1,ms2,ms3,ms4)+nc*q_coord(iq_ibz,4)
!                   if(m1==1.and.m2==1.and.m3==1.and.m4==1) then
!                     write(6,*) "TEST2"
!                     write(6,*) iq_ibz,nc
!                     write(6,*) q_coord(iq_ibz,4)
!                   endif
                    nc=cmplx(0,0)
                  endif!paral
                end do!m4
                end do!m3
                end do!m2
                end do!m1
              end do!ispinor4
              end do!ispinor3
              end do!ispinor2
              end do!ispinor1
            end do!iq_ibz
            U_m(:,:,:,:)=U_m(:,:,:,:)*Ha_eV/(ucvol)
            call xmpi_barrier(Wfd%comm)  ! First make sure that all processors are here
            call xmpi_sum(U_m,Wfd%comm,ierr)
            call xmpi_barrier(Wfd%comm)  ! First make sure that all processors are here
            call cpu_time ( t2 )
            write(message,*)  "in ",t2-t1, "sec"
            call wrtout(std_out,message,'COLL')
          END SELECT
    ! tolerance of the symetry of screened U.
          tol=1E-2
          call checkk(U_m,1,mbband*nspinor,0,iomega,uu,jj,"UminusVbare",mbband1,mbband2,mbband3,mbband4,nspinor,one_orbital)
          U_m=V_m+U_m
          write(message,*)  "UCRPA interaction"
          call wrtout(std_out,message,'COLL')
          call checkk(U_m,1,mbband*nspinor,1,iomega,uomega(iomega),jomega(iomega),&
            &"cRPA interaction",mbband1,mbband2,mbband3,mbband4,nspinor,one_orbital)
          if (spin1==1 .and. spin2==1) then
            ispin=1
          else if(spin1==1 .and. spin2==2) then
            ispin=2
          else if (spin1==2 .and. spin2==1)then
            ispin=3
          else
            ispin=4
          endif
          if (wanbz%nsppol/=1)then
            uspin(ispin,iomega)=uomega(iomega)
            jspin(ispin,iomega)=jomega(iomega)
          endif
          !write(message,*)ch10,"SCREENED INTERACTION"
          !call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
          !call print_U(mbband1,mbband2,mbband3,mbband4,nspinor,U_m)
        enddo

!   Summarize the calculation if nsppol==1
        if (wanbz%nsppol /=1 .and. spin1==2 .and. spin2==2)then
          uomega(:)=sum(uspin,dim=1)/4
          jomega(:)=sum(jspin,dim=1)/4
        endif
        call print_orbitals(spin1,spin2,iatom1,iatom2,iatom3,iatom4,pos1,pos2,pos3,pos4,il1,il2,il3,il4,wanbz,0)
        write(message,*)ch10,"  -------------------------------------------------------------"
        call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
        write(message,*)"           Average U and J as a function of frequency   "
        call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
        write(message,*)"  -------------------------------------------------------------"
        call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
        write(message,*)"        omega           U(omega)            J(omega)"
        call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
        do iomega=1,nomega
          if(nomega==1) then
            omega(iomega)=omegamin
          else
            omega(iomega)=(omegamax-omegamin)/(nomega-1)*(iomega-1)+omegamin
          endif
          if(nomega==1)  omega(iomega)=omegamin
          write(message,'(2x,f11.3,2x,2f10.4,2x,2f10.4)')  omega(iomega)*Ha_eV, uomega(iomega),jomega(iomega)
          call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
        enddo
        write(message,*)"  -------------------------------------------------------------"
        call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
        
!  END OF LOOP ON ATOMS !!
        
        ABI_DEALLOCATE(uomega)
        ABI_DEALLOCATE(jomega)
        ABI_DEALLOCATE(V_m)
        ABI_DEALLOCATE(rhot_q_m1m3_npwe)
        ABI_DEALLOCATE(rhot_q_m2m4_npwe)
        ABI_DEALLOCATE(trrho)
        ABI_DEALLOCATE(sumrhorhoeps)
        ABI_DEALLOCATE(scr)
        ABI_DEALLOCATE(U_m)
        ABI_DEALLOCATE(rhot_q_m1m3)
        ABI_DEALLOCATE(rhot_q_m2m4)
!    Print dielectric matrix
!   do iq_ibz=1,nqibz
!       call read_screening(em1_ncname,fname,npwe,1,nomega,scr,IO_MODE_FORTRAN,Wfd%comm,iqiA=iq_ibz)
!   enddo
      enddo!spin2
      enddo!spin1
      if (wanbz%nsppol/=1)then
        do iomega=1,nomega
          if(nomega==1) then
            omega(iomega)=omegamin
          else
            omega(iomega)=(omegamax-omegamin)/(nomega-1)*(iomega-1)+omegamin
          endif
          if(nomega==1)  omega(iomega)=omegamin
        enddo
        call print_orbitals(1,1,iatom1,iatom2,iatom3,iatom4,pos1,pos2,pos3,pos4,il1,il2,il3,il4,wanbz,2)
        call print_uj_spin(nomega,uspin,jspin,omega,one_orbital)
      endif
      if (wanbz%nsppol/=1)then
        ABI_DEALLOCATE(uspin)
        ABI_DEALLOCATE(jspin)
      endif
      ABI_DEALLOCATE(omega)
      enddo!il4
      enddo!il3
      enddo!il2
      enddo!il1
    enddo!pos4
    enddo!pos3
    enddo!pos2
    enddo!pos1
  enddo!iatom4
  enddo!iatom3
  enddo!iatom2
  enddo!iatom1
 ABI_DEALLOCATE(k_coord)
 ABI_DEALLOCATE(q_coord)
 ABI_DEALLOCATE(bijection)
 ABI_DEALLOCATE(normG)
! ABI_DEALLOCATE(coeffW_BZ)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 CONTAINS

!!    __                 _   _
!!   / _|               | | (_)
!!  | |_ ___  _ __   ___| |_ _  ___  _ __
!!  |  _/ _ \| '_ \ / __| __| |/ _ \| '_ \
!!  | || (_) | | | | (__| |_| | (_) | | | |
!!  |_| \___/|_| |_|\___|\__|_|\___/|_| |_|
!
 integer FUNCTION fi(nkbz,k_coord,kprime_coord)

      implicit none
       integer,intent(in) :: nkbz
      real(dp),dimension(nkbz,3),intent(in) ::k_coord
      real(dp),dimension(3),intent(in) :: kprime_coord(3)

 do fi=1,nkbz
   if (ALL(abs(kprime_coord(:)-k_coord(fi,:))<0.001)) then
     exit
   end if
 end do
 END FUNCTION fi

 integer FUNCTION findkmq(ik_bz,k_coord,q_coord,nkbz)

      implicit none
      integer,intent(in) :: ik_bz,nkbz
      real(dp),dimension(nkbz,3),intent(in) ::k_coord
      real(dp),dimension(4),intent(in) ::q_coord
      real(dp),dimension(3) :: kprime_coord
      integer :: i,j,k


 kprime_coord(:)=k_coord(ik_bz,:)-q_coord(1:3)

 where (kprime_coord > 0.5)
   kprime_coord(:)= kprime_coord(:)-1
 elsewhere(kprime_coord < -0.5)
   kprime_coord(:)= kprime_coord(:)+1
 end where

!! indice of k -q
 findkmq=fi(nkbz,k_coord,kprime_coord)

!!Test if k-q exists
 if (findkmq.EQ.(nkbz+1)) then
!!The prb comes from PBC included born_inf,born_sup]
!! One test all combination over boundaries
   do i=1,2
     if (abs(abs(kprime_coord(1))-0.5)<0.001) kprime_coord(1)=(-1)**i*0.5
     do j=1,2
       if (abs(abs(kprime_coord(2))-0.5)<0.001) kprime_coord(2)=(-1)**j*0.5
       do k=1,2
         if (abs(abs(kprime_coord(3))-0.5)<0.01) kprime_coord(3)=(-1)**k*0.5
         findkmq=fi(nkbz,k_coord,kprime_coord)
        !!Quand on a trouver la bonne valeur on part
         if (findkmq.NE.(nkbz+1)) return
       end do
     end do
   end do
 end if
 END FUNCTION findkmq

 SUBROUTINE checkk(Interaction,m_inf,m_sup,prtopt,ifreq,uu,jj,utype,mbband1,mbband2,mbband3,mbband4,nspinor,one_orbital)

 implicit none
 integer, intent(in) :: m_inf,m_sup,ifreq,mbband1,mbband2,mbband3,mbband4,nspinor,one_orbital
 complex(dpc), intent(in) :: Interaction(mbband1*nspinor,mbband2*nspinor,mbband3*nspinor,mbband4*nspinor)
 complex(dpc), intent(out)    :: uu,jj
 character(len=*), intent(in) :: utype
 integer :: prtopt


!==Check correctness
 ! write(message,*)  "== Check == "
 ! call wrtout(std_out,message,'COLL')
 ! write(message,*) 'Tolerance :',tol
 ! call wrtout(std_out,message,'COLL')
 ! do i=m_inf,m_sup
 !   do j=i+1,m_sup
      !if (abs(abs(Interaction(i,i,i,i)-Interaction(j,j,j,j)))>tol) then
      !     BUG=.TRUE.
      !     write(message,*) "Problem in the interband calculation"&
!&      ,i,j,abs(Interaction(i,i,i,i)),abs(Interaction(j,j,j,j))
      !     call wrtout(std_out,message,'COLL')
      !end if

!      if (abs(Interaction(i,j,i,j)-Interaction(j,i,j,i))>tol) then
!        BUG=.TRUE.
!        write(message,*) "Warning in the symetry of U'",i,j,abs(Interaction(i,j,i,j)),&
! &       abs(Interaction(j,i,j,i)),abs(Interaction(i,j,i,j)-Interaction(j,i,j,i))
!        call wrtout(std_out,message,'COLL')
!      end if

!      if (abs(Interaction(i,i,j,j)-Interaction(j,j,i,i))>tol) then
!        BUG=.TRUE.
!        write(message,*) "Warning in the symetry of J'",i,j,abs(Interaction(i,i,j,j)),&
! &       abs(Interaction(j,j,i,i)),abs(Interaction(j,j,i,i)-Interaction(i,i,j,j))
!        call wrtout(std_out,message,'COLL')
!      end if
!    end do
!  end do


! do i=m_inf,m_sup
!   do j=m_inf,m_sup
!     if (i.EQ.j) cycle
!     if (abs(Interaction(i,j,j,j))>tol) then
!       BUG=.TRUE.
!       write(message,*) "Warning in the symetry U(i,j,j,j) should vanish (in the Ylm basis) for",&
!&       i,j,abs(Interaction(i,j,j,j))
!       call wrtout(std_out,message,'COLL')
!     end if
!   end do
! end do

! if (.not.BUG) then
!    call wrtout(std_out,'Calcul is possibly correct','COLL')
!    call Affichage(Interaction,m_inf,m_sup,2)
! else
!     call wrtout(std_out,'Maybe somme error','COLL')
!     call Affichage(Interaction,m_inf,m_sup,1)
! end if

 if(prtopt>0)  call Affichage(Interaction,m_inf,m_sup,1,ifreq,uu,jj,utype,mbband1,mbband2,mbband3,mbband4,nspinor,one_orbital)

 END SUBROUTINE checkk

 SUBROUTINE Affichage(Interaction,m_inf,m_sup,option,ifreq,uu,jj,utype,mbband1,mbband2,mbband3,mbband4,nspinor,one_orbital)

  implicit none
  integer, intent(in) :: m_inf,m_sup,option,ifreq,mbband1,mbband2,mbband3,mbband4,nspinor,one_orbital
  complex(dpc), intent(in) :: Interaction(mbband1*nspinor,mbband2*nspinor,mbband3*nspinor,mbband4*nspinor)
  complex(dpc),intent(out) :: UU,JJ
  character(len=*), intent(in) :: utype
  complex(dpc) :: UU1,UUmJJ,JJ1,JJ2
  integer :: m1,m2
  logical :: lprint
  character(len=500) :: message


  if(utype=="cRPA interaction".or.utype=="bare interaction") then
   lprint=.true.
  else
   lprint=.false.
  endif
  write(message,*) "";call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')

  if (one_orbital==1) then
    if(lprint) then
      write(message,*)" Diagonal ",utype
      call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
    endif
    
    do m1=m_inf,m_sup
      if (option.EQ.1) then
        write(message,'(a,i3,14f7.3)') " ",m1,real(Interaction(m1,m1,m1,m1))
        call wrtout(std_out,message,'COLL')
        call wrtout(ab_out,message,'COLL')
      end if
    end do
  endif
  
  if(lprint) then
    write(message,*) "";call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
    write(message,*)" U'=U(m1,m2,m1,m2) for the ",utype
    call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')

    write(message,'(a,14i7)') " -",(m2,m2=1,mbband2)
    call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
    do m1=m_inf,m_sup
      if (option.EQ.1) then
        write(message,'(a,i3,14f7.3)') " ",m1,(real(Interaction(m1,m2,m1,m2)),m2=1,mbband2)
        call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
      end if
    end do
  endif

 write(message,*) ""
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 UU=czero
 do m1=1,mbband1
   do m2=1,mbband2
     UU=UU+Interaction(m1,m2,m1,m2)
   enddo
 enddo
 UU=UU/float(mbband1*mbband2)
 if(ifreq/=0) write(message,'(3a,i3,a,2f10.4,a)')'  Hubbard ',utype,' for w =',ifreq,', U=1/(2l+1)**2 \sum U(m1,m2,m1,m2)=',UU,ch10
 if(ifreq==0) write(message,'(3a,2f10.4,a)') '  Hubbard ',utype, ' U=1/(2l+1)**2 \sum U(m1,m2,m1,m2)=',UU,ch10
 call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
 

if (one_orbital==1)then
  UU1=czero
  do m1=1,mbband1
    UU1=UU1+Interaction(m1,m1,m1,m1)
  enddo
  UU1=UU1/((mbband1))
  if(ifreq/=0) write(message,'(3a,i4,a,2f10.4,2a)')&
    & '  (Hubbard ',utype,' for w =',ifreq,', U=1/(2l+1) \sum U(m1,m1,m1,m1)=',UU1,')',ch10
  if(ifreq==0) write(message,'(3a,2f10.4,2a)')' (Hubbard ',utype,' U=1/(2l+1) \sum U(m1,m1,m1,m1)=',UU1,')',ch10
  call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
endif

 if(lprint .and. one_orbital==1) then
   write(message,*)' Hund coupling J=U(m1,m1,m2,m2) for the ', utype
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')

   write(message,'(a,14i7)') " -",(m2,m2=1,mbband1)
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   do m1=1,mbband1
     if (option.EQ.1) then
        write(message,'(a,i3,14f7.3)') " ",m1,(real(Interaction(m1,m1,m2,m2)),m2=1,mbband1)
        call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
     end if
   end do
 endif
 if (one_orbital==1)then
   UUmJJ=czero
   do m1=1,mbband1
     do m2=1,mbband1
       UUmJJ=UUmJJ+Interaction(m1,m2,m1,m2)-Interaction(m1,m2,m2,m1)
     enddo
   enddo
   if (mbband1/=1) then
     UUmJJ=UUmJJ/float((mbband1)*(mbband1-1))
     JJ1=UU-UUmJJ
   endif
   
   
   JJ=czero
   do m1=1,mbband1
     do m2=1,mbband1
       if(m1/=m2) JJ=JJ+Interaction(m1,m2,m2,m1)
     enddo
   enddo
   if (mbband1/=1) then
     JJ=JJ/float((mbband1)*(mbband1-1))
   endif

   write(message,'(a,3x,2a,2f10.4,a)')ch10,utype,&
     & ' value of J=U-1/((2l+1)(2l)) \sum_{m1,m2} (U(m1,m2,m1,m2)-U(m1,m2,m2,m1))=',JJ1,ch10
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   
   if (mbband1/=1) then
     UUmJJ=czero
     do m1=1,mbband1
       do m2=1,mbband1
         UUmJJ=UUmJJ+Interaction(m1,m2,m1,m2)-Interaction(m1,m1,m2,m2)
       enddo
     enddo
     UUmJJ=UUmJJ/float((mbband1)*(mbband1-1))
     JJ2=UU-UUmJJ
     if(abs(JJ1-JJ2)<0.0001) then
       JJ=JJ1
     else
     
     write(message,'(a,3x,2a,2f10.4,a)')ch10,utype,&
       &   ' value of J=U-1/((2l+1)(2l)) \sum_{m1,m2} (U(m1,m2,m1,m2)-U(m1,m1,m2,m2))=',JJ2,ch10
     call wrtout(std_out,message,'COLL')
     stop
   endif
 endif
endif
 
 if(lprint .and. one_orbital==1) then
   write(message,*)ch10,' Hund coupling J2=U(m1,m2,m2,m1) for the ', utype
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')

   write(message,'(a,14i7)') " -",(m2,m2=1,mbband1)
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   do m1=1,mbband1
     if (option.EQ.1) then
        write(message,'(a,i3,14f7.3)') " ",m1,(real(Interaction(m1,m2,m2,m1)),m2=1,mbband1)
        call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
     end if
   end do
! write(message,*) "";call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
! write(message,*) "U'=U-2J for the t2g should be checked"
! call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
 endif

 END SUBROUTINE Affichage

 SUBROUTINE Sauvegarde_M_q_m(M_q_m,normG,nqibz,npw,mbband)

 implicit none
 integer, intent(in) :: nqibz,npw,mbband
 complex(dpc), intent(in) :: M_q_m(mbband,mbband,npw,nqibz)
 real(dp), intent(in) :: normG(npw)
 integer :: i,j,iq_ibz,iG,unt
 character(len=500) :: msg

!==Ecriture de M_(G=0)^(mm')(q) ==!
 if (open_file('M_mimj(n=1_2_3)(q,G=0)',msg,newunit=unt,form='formatted',status='unknown') /=0) then
   MSG_ERROR(msg)
 end if
 do iq_ibz=1,nqibz
    write(unt,*) iq_ibz,((abs(M_q_m(i,j,1,iq_ibz)),i=1,mbband),j=1,mbband)
 end do
 close(unt)

!==Ecriture de M_G^(mm')(q=0) ==!
 if (open_file('M_mm(m=1..mbband)(q=0)',msg,newunit=unt,form='formatted',status='unknown') /= 0) then
   MSG_ERROR(msg)
 end if
 do iG=1,npw
    write(unt,*) normG(iG),(abs(M_q_m(i,i,iG,1)),i=1,mbband)
 end do
 close(unt)
 END SUBROUTINE Sauvegarde_M_q_m

 subroutine print_orbitals(spin1,spin2,iatom1,iatom2,iatom3,iatom4,pos1,pos2,pos3,pos4,il1,il2,il3,il4,wanbz,opt)
   implicit none
   integer, intent(in) :: spin1,spin2,iatom1,iatom2,iatom3,iatom4,pos1,pos2,pos3,pos4,il1,il2,il3,il4,opt
   type(plowannier_type),intent(in) :: wanbz
   character(len=5000) ::message
   character(len=10):: print_spin

   if (spin1==spin2 .and. spin1==1)then
     print_spin=" Up-Up"
   else if(spin1==1 .and. spin2==2)then
     print_spin=" Up-Down"
   else if(spin1==2 .and. spin2==1) then
     print_spin=" Down-Up"
   else 
     print_spin=" Down-Down"
   endif
   if (opt==1)then
     write(message,*)ch10,"==Definition of the orbitals=="
   else if (opt==0) then 
     write(message,*)ch10,"==Reminder of the orbitals=="
   else if (opt==2) then
     write(message,*)ch10,"==Reminder of the orbitals=="
     write(print_spin,*)" summary"
   else if (opt==3) then
     write(message,*)ch10,"==Reminder of the orbitals=="
     write(print_spin,*)" average"
   endif
   call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
   if (iatom1==iatom2.and. iatom3==iatom4 .and. iatom3==iatom2)then
     if(pos1==pos2 .and. pos3==pos4 .and. pos1==pos3)then
       if (il1==il2 .and. il3==il4 .and. il3==il1)then
         write(message,*)" Only one orbital"
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)" Orbital with l=",wanbz%latom_wan(iatom1)%lcalc(il1),&
           &"on atom",wanbz%iatom_wan(iatom1),"with spin's orientations",print_spin
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
       else
         write(message,*)"Different orbitals on atom",wanbz%iatom_wan(iatom1),"with spin's orientations",print_spin
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)"  orbital 1 with l=",wanbz%latom_wan(iatom1)%lcalc(il1)
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)"  orbital 2 with l=",wanbz%latom_wan(iatom1)%lcalc(il2)
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)"  orbital 3 with l=",wanbz%latom_wan(iatom1)%lcalc(il3)
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)"  orbital 4 with l=",wanbz%latom_wan(iatom1)%lcalc(il4)
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
       endif
     else
       if(il1==il2 .and. il3==il4 .and. il1==il3)then
         write(message,*)" Different position of the same orbital"
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)" Orbitals with l=",wanbz%latom_wan(iatom1)%lcalc(il1),&
           &"on atom",wanbz%iatom_wan(iatom1),"with spin's orientations",print_spin
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)"  orbital 1 at postion ",wanbz%nposition(iatom1)%pos(pos1,:)
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)"  orbital 2 at postion ",wanbz%nposition(iatom1)%pos(pos2,:)
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)"  orbital 3 at postion ",wanbz%nposition(iatom1)%pos(pos3,:)
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)"  orbital 4 at postion ",wanbz%nposition(iatom1)%pos(pos4,:)
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
       else 
         write(message,*)" Different orbitals on the same atom, with different positions"
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)" Orbitals on atom",wanbz%iatom_wan(iatom1),"with spin's orientations",print_spin
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)"  orbital 1 with l=",wanbz%latom_wan(iatom1)%lcalc(il1),"at position ",wanbz%nposition(iatom1)%pos(pos1,:)
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)"  orbital 2 with l=",wanbz%latom_wan(iatom1)%lcalc(il2),"at position ",wanbz%nposition(iatom1)%pos(pos2,:)
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)"  orbital 3 with l=",wanbz%latom_wan(iatom1)%lcalc(il3),"at position ",wanbz%nposition(iatom1)%pos(pos3,:)
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
         write(message,*)"  orbital 4 with l=",wanbz%latom_wan(iatom1)%lcalc(il4),"at position ",wanbz%nposition(iatom1)%pos(pos4,:)
         call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
       endif
     endif
   else
     if(pos1==pos2 .and. pos3==pos4 .and. pos1==pos3)then
       write(message,*)"Different atoms, with spin orientation",print_spin
       call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
       write(message,*)"  orbital 1 with l=",wanbz%latom_wan(iatom1)%lcalc(il1),"on atom",wanbz%iatom_wan(iatom1)
       call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
       write(message,*)"  orbital 2 with l=",wanbz%latom_wan(iatom2)%lcalc(il2),"on atom",wanbz%iatom_wan(iatom2)
       call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
       write(message,*)"  orbital 3 with l=",wanbz%latom_wan(iatom3)%lcalc(il3),"on atom",wanbz%iatom_wan(iatom3)
       call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
       write(message,*)"  orbital 4 with l=",wanbz%latom_wan(iatom4)%lcalc(il4),"on atom",wanbz%iatom_wan(iatom4)
       call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
     else
       write(message,*)"Different atoms, in different postions with spin's orientations",print_spin
       call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
       write(message,*)"  orbital 1 with l=",wanbz%latom_wan(iatom1)%lcalc(il1),&
         &"on atom",wanbz%iatom_wan(iatom1),"at position ",wanbz%nposition(iatom1)%pos(pos1,:)
       call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
       write(message,*)"  orbital 2 with l=",wanbz%latom_wan(iatom2)%lcalc(il2),&
         &"on atom",wanbz%iatom_wan(iatom2),"at position ",wanbz%nposition(iatom2)%pos(pos2,:)
       call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
       write(message,*)"  orbital 3 with l=",wanbz%latom_wan(iatom3)%lcalc(il3),&
         &"on atom",wanbz%iatom_wan(iatom3),"at position ",wanbz%nposition(iatom3)%pos(pos3,:)
       call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
       write(message,*)"  orbital 4 with l=",wanbz%latom_wan(iatom4)%lcalc(il4),&
         &"on atom",wanbz%iatom_wan(iatom4),"at position ",wanbz%nposition(iatom4)%pos(pos4,:)
       call wrtout(std_out,message,'COLL');call wrtout(ab_out,message,'COLL')
     endif
   endif
 end subroutine print_orbitals


 subroutine print_uj_spin(nomega,uspin,jspin,omega,one_orbital)
   implicit none
   integer,intent(in) :: nomega,one_orbital
   complex(dpc),intent(in) :: uspin(4,nomega)
   complex(dpc),intent(in) :: jspin(4,nomega)
   real(dp),intent(in) :: omega(nomega)
   integer :: iomega,ispin
   real(dp) :: uomega(nomega),jomega(nomega)
   character(len=500)::message
   
   write(message,*)ch10," --------------------------------------------------------------------" 
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   write(message,*)" Sum up of the calcul for different spin polarization and frequencies"
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   write(message,*)" --------------------------------------------------------------------" 
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   
   
   write(message,*)ch10,"Sum up of U",ch10
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   
   write(message,'(6a)')" -omega (eV)- "," Up-Up "," Up-Down "," Down-Up "," Down-Down "," Average "
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   do iomega=1,nomega
     write(message,'(a,f7.2,5f9.3)')"    ",omega(iomega)*Ha_eV,(real(uspin(ispin,iomega)),ispin=1,4),sum(real(uspin(:,iomega)))/4
     call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   end do
   
   if (one_orbital==1)then
     write(message,*)ch10,"Sum up of J",ch10
     call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
     
     write(message,'(6a)')" -omega (eV)- "," Up-Up "," Up-Down "," Down-Up "," Down-Down "," Average "
     call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
     do iomega=1,nomega
       write(message,'(a,f7.2,5f9.3)')"    ",omega(iomega)*Ha_eV,(real(jspin(ispin,iomega)),ispin=1,4),sum(real(jspin(:,iomega)))/4
       call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
     end do
   endif
   uomega(:)=sum(uspin,dim=1)/4
   jomega(:)=sum(jspin,dim=1)/4
   call print_orbitals(spin1,spin2,iatom1,iatom2,iatom3,iatom4,pos1,pos2,pos3,pos4,il1,il2,il3,il4,wanbz,3)
   write(message,*)ch10,"  -------------------------------------------------------------"
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   write(message,*)"           Average U and J as a function of frequency   "
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   write(message,*)"  -------------------------------------------------------------"
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   write(message,*)"        omega           U(omega)            J(omega)"
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   do iomega=1,nomega
     write(message,'(2x,f11.3,2x,2f10.4,2x,2f10.4)')  omega(iomega)*Ha_eV, uomega(iomega),jomega(iomega)
     call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
   enddo
   write(message,*)"  -------------------------------------------------------------"
   call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')
 end subroutine print_uj_spin
           
 end subroutine calc_ucrpa
!!***

END MODULE m_calc_ucrpa
!!***
