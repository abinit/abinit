!!****m* m_paw_correlations/m_paw_correlations
!! NAME
!!  m_paw_correlations
!!
!! FUNCTION
!!  This module contains several routines related to the treatment of electronic
!!    correlations in the PAW approach (DFT+U, exact-exchange, ...).
!!
!! COPYRIGHT
!! Copyright (C) 2018-2019 ABINIT group (BA,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_correlations

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_dtset
 use m_linalg_interfaces
 use m_special_funcs

 use m_io_tools,    only : open_file
 use m_pawang,      only : pawang_type,pawang_init,pawang_free
 use m_pawrad,      only : pawrad_type,simp_gen,nderiv_gen,pawrad_ifromr,poisson
 use m_pawtab,      only : pawtab_type
 use m_pawrhoij,    only : pawrhoij_type
 use m_paw_ij,      only : paw_ij_type
 use m_paw_sphharm, only : mat_mlms2jmj,mat_slm2ylm
 use m_paw_io,      only : pawio_print_ij
 use m_paral_atom,  only : get_my_atmtab,free_my_atmtab

 implicit none

 private

!public procedures.
 public :: pawpuxinit   ! Initialize some data for PAW+U/PAW+LocalExactExchange/PAW+DMFT
 public :: pawuenergy   ! Compute contributions to energy for PAW+U
 public :: pawxenergy   ! Compute contributions to energy for PAW+[local exact exchange]
 public :: setnoccmmp   ! Compute LDA+U density matrix nocc_{m,m_prime} or impose it
 public :: setrhoijpbe0 ! Impose value of rhoij for using an auxiliairy file (PBE0 only)
 public :: calc_ubare   ! Calculate the bare interaction on atomic orbitals

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_correlations/pawpuxinit
!! NAME
!! pawpuxinit
!!
!! FUNCTION
!! Initialize some starting values of several arrays used in
!! PAW+U/+DMFT or local exact-exchange calculations
!!
!! A-define useful indices for LDA+U/local exact-exchange
!! B-Compute overlap between atomic wavefunction
!! C-Compute matrix elements of coulomb interaction (see PRB vol.52 5467) [[cite:Liechenstein1995]]
!!    (angular part computed from Gaunt coefficients)
!!
!! INPUTS
!!  dmatpuopt= select expression for the density matrix
!!  exchmix= mixing factor for local exact-exchange
!!  is_dfpt=true if we are running a DFPT calculation
!!  jpawu(ntypat)= value of J
!!  llexexch(ntypat)= value of l on which local exact-exchange applies
!!  llpawu(ntypat)= value of l on which LDA+U applies
!!  ntypat=number of types of atoms in unit cell.
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!     %lmax=Maximum value of angular momentum l+1
!!     %gntselect((2*l_max-1)**2,l_max**2,l_max**2)=
!!                     selection rules for Gaunt coefficients
!!  pawprtvol=output printing level for PAW
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data:
!!  upawu(ntypat)= value of U
!!  use_dmft = 0 no PAW+DMFT, =1 PAW+DMFT
!!  useexexch= 0 if no local exact-exchange; 1 if local exact-exchange
!!  usepawu= 0 if no LDA+U; /=0 if LDA+U
!!
!! OUTPUT
!!  pawtab <type(pawtab_type)>=paw tabulated data read at start:
!!     %euijkl=(2,2,lmn2_size,lmn2_size)= array for computing LDA+U terms without occupancies
!!     %ij_proj= nproj*(nproju+1)/2
!!     %klmntomn(4,lmn2_size)= Array giving im, jm ,in, and jn for each klmn=(ilmn,jlmn)
!!     %lnproju(nproj)= value of ln for projectors on which paw+u/local exact-exchange acts.
!!     %nproju=number of projectors for orbitals on which paw+u/local exact-exchange acts.
!!     %phiphjint(pawtabitypat%ij_proj)=Integral of Phi(:,i)*Phi(:,j) for correlated orbitals.
!!     %usepawu=0 if no LDA+U; /=0 if LDA+U
!!     %useexexch=0 if no local exact-exchange; 1 if local exact-exchange
!!     === if usepawu/=0
!!     %jpawu= value of J
!!     %upawu= value of U
!!     %vee(2*lpawu+1,:,:,:)=matrix of the screened interaction for correlated orbitals
!!     === if useexexch/=0
!!     %fk
!!     %vex(2*lpawu+1,:,:,:)=matrix of the screened interaction for correlated orbitals
!!
!! PARENTS
!!      bethe_salpeter,gstate,m_entropyDMFT,respfn,screening,sigma
!!
!! CHILDREN
!!      calc_ubare,poisson,simp_gen,wrtout
!!
!! SOURCE

 subroutine pawpuxinit(dmatpuopt,exchmix,f4of2_sla,f6of2_sla,is_dfpt,jpawu,llexexch,llpawu,&
&           ntypat,pawang,pawprtvol,pawrad,pawtab,upawu,use_dmft,useexexch,usepawu,&
&           ucrpa) ! optional argument

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: dmatpuopt,ntypat,pawprtvol,use_dmft,useexexch,usepawu
 logical :: is_dfpt
 real(dp),intent(in) :: exchmix
 type(pawang_type), intent(in) :: pawang
 integer,optional, intent(in) :: ucrpa
!arrays
 integer,intent(in) :: llexexch(ntypat),llpawu(ntypat)
 real(dp),intent(in) :: jpawu(ntypat),upawu(ntypat)
 real(dp),intent(in) :: f4of2_sla(ntypat),f6of2_sla(ntypat)
 type(pawrad_type),intent(inout) :: pawrad(ntypat)
 type(pawtab_type),target,intent(inout) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: icount,il,ilmn,ilmnp,isela,iselb,itemp,itypat,iu,iup,j0lmn,jl,jlmn,jlmnp,ju,jup
 integer :: klm0u,klm0x,klma,klmb,klmn,klmna,klmnb,kln,kln1,kln2,kyc,lcur,lexexch,lkyc,ll,ll1
 integer :: lmexexch,lmkyc,lmn_size,lmn2_size,lmpawu,lpawu
 integer :: m1,m11,m2,m21,m3,m31,m4,m41
 integer :: mesh_size,int_meshsz,mkyc,sig,sigp,sz1
 logical :: compute_euijkl,compute_euij_fll
 real(dp) :: ak,f4of2,f6of2,int1,intg,phiint_ij,phiint_ipjp,vee1,vee2
 character(len=500) :: message
!arrays
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp) :: euijkl_temp(2,2),euijkl_temp2(2,2),euijkl_dc(2,2)
 real(dp),allocatable :: ff(:),fk(:),gg(:)

! *************************************************************************

 DBG_ENTER("COLL")

!No correlations= nothing to do
 if(useexexch==0.and.usepawu==0.and.use_dmft==0) then
   do itypat=1,ntypat
     pawtab(itypat)%usepawu=0;pawtab(itypat)%useexexch=0;pawtab(itypat)%exchmix=zero
   end do
   return
 end if

!PAW+U and local exact-exchange restriction
 if(useexexch/=0.and.usepawu/=0)then
   do itypat=1,ntypat
     if (llpawu(itypat)/=llexexch(itypat).and.llpawu(itypat)/=-1.and.llexexch(itypat)/=-1) then
       write(message, '(5a,i2,3a)' )&
&       '  When PAW+U (usepawu/=0) and local exact-exchange (exexch/=0)',ch10,&
&       '  are selected together, they must apply on the same',ch10,&
&       '  angular momentum (lpawu/=lexexch forbidden, here for typat=',itypat,') !',ch10,&
&       '  Action: correct your input file.'
       MSG_ERROR(message)
     end if
   end do
 end if

!Print title
 if((abs(usepawu)>=1.and.abs(usepawu)<=4).or.useexexch/=0) &
&  write(message, '(3a)' ) ch10,ch10," ******************************************"
 if(usepawu==1) then
   write(message, '(3a)' ) trim(message),ch10," LDA+U Method used: FLL"
 else if(usepawu==2) then
   write(message, '(3a)' ) trim(message),ch10," LDA+U Method used: AMF"
 else if(usepawu==3) then
   write(message, '(3a)' ) trim(message),ch10," LDA+U Method used: AMF (alternative)"
 else if(usepawu==4) then
   write(message, '(3a)' ) trim(message),ch10," LDA+U Method used: FLL with no spin polarization in the xc functional"
 else if(usepawu==-1) then
   write(message, '(3a)' ) trim(message),ch10," LDA+U Method used: FLL (no use of occupation matrix) - experimental"
 else if(usepawu==-2) then
   write(message, '(3a)' ) trim(message),ch10," LDA+U Method used: AMF (no use of occupation matrix) - experimental"
 end if
 if(useexexch/=0) write(message, '(3a)' ) trim(message),ch10," PAW Local Exact exchange: PBE0"
 if((abs(usepawu)>=1.and.abs(usepawu)<=4).or.useexexch/=0) &
&  write(message, '(3a)' ) trim(message),ch10," ******************************************"
 if(use_dmft==0) then
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if
!if(use_dmft>0) then
!write(message, '(3a)' ) ch10, " (see DMFT data in log file) "
!call wrtout(ab_out,message,'COLL')
!endif

!Loop on atom types
 do itypat=1,ntypat
   indlmn => pawtab(itypat)%indlmn
   lmn_size=pawtab(itypat)%lmn_size
   lmn2_size=pawtab(itypat)%lmn2_size
   mesh_size=pawtab(itypat)%mesh_size
   int_meshsz=pawrad(itypat)%int_meshsz
   lcur=-1

!  PAW+U data
   if (usepawu/=0.or.use_dmft>0) then
     lcur=llpawu(itypat)
     pawtab(itypat)%lpawu=lcur
     if(lcur/=-1) then
       pawtab(itypat)%usepawu=usepawu
       pawtab(itypat)%upawu=upawu(itypat)
       pawtab(itypat)%jpawu=jpawu(itypat)
       pawtab(itypat)%f4of2_sla=f4of2_sla(itypat)
       pawtab(itypat)%f6of2_sla=f6of2_sla(itypat)
     else
       pawtab(itypat)%usepawu=0
       pawtab(itypat)%upawu=zero
       pawtab(itypat)%jpawu=zero
       pawtab(itypat)%f4of2_sla=zero
       pawtab(itypat)%f6of2_sla=zero
     end if
   end if

!  Local exact-echange data
   if (useexexch/=0) then
     lcur=llexexch(itypat)
     pawtab(itypat)%lexexch=lcur
     pawtab(itypat)%exchmix=exchmix
     if(pawtab(itypat)%lexexch==-1) pawtab(itypat)%useexexch=0
     if(pawtab(itypat)%lexexch/=-1) pawtab(itypat)%useexexch=useexexch
   end if

!  Select only atoms with +U
   if(lcur/=-1) then

!    Compute number of projectors for LDA+U/local exact-exchange/LDA+DMFT
     icount=count(indlmn(1,1:lmn_size)==lcur)
     pawtab(itypat)%nproju=icount/(2*lcur+1)
     if(useexexch/=0.and.pawtab(itypat)%nproju>2)  then
       write(message, '(a,a,a)' )&
&       '  Error on the number of projectors ',ch10,&
&       '  more than 2 projectors is not allowed for local exact-exchange'
       MSG_ERROR(message)
     end if
     if(pawtab(itypat)%nproju*(2*lcur+1)/=icount)  then
       message = 'pawpuxinit: Error on the number of projectors '
       MSG_BUG(message)
     end if
     write(message, '(a,a,i4,a,a,i4)' ) ch10,&
&     ' pawpuxinit : for species ',itypat,ch10,&
&     '   number of projectors is',pawtab(itypat)%nproju
     call wrtout(std_out,message,'COLL')

     pawtab(itypat)%ij_proj=pawtab(itypat)%nproju*(pawtab(itypat)%nproju+1)/2

!    ==================================================
!    A-define useful indexes
!    --------------------------------------------------
     if (allocated(pawtab(itypat)%lnproju)) then
       ABI_DEALLOCATE(pawtab(itypat)%lnproju)
     end if
     ABI_ALLOCATE(pawtab(itypat)%lnproju,(pawtab(itypat)%nproju))
     icount=0
     do ilmn=1,lmn_size
       if(indlmn(1,ilmn)==lcur) then
         itemp=icount/(2*lcur+1)
         if (itemp*(2*lcur+1)==icount) then
           pawtab(itypat)%lnproju(itemp+1)=indlmn(5,ilmn)
         end if
         icount=icount+1
       end if
     end do

     if (allocated(pawtab(itypat)%klmntomn)) then
       ABI_DEALLOCATE(pawtab(itypat)%klmntomn)
     end if
     ABI_ALLOCATE(pawtab(itypat)%klmntomn,(4,lmn2_size))
     do jlmn=1,lmn_size
       jl= indlmn(1,jlmn)
       j0lmn=jlmn*(jlmn-1)/2
       do ilmn=1,jlmn
         il= indlmn(1,ilmn)
         klmn=j0lmn+ilmn
         pawtab(itypat)%klmntomn(1,klmn)=indlmn(2,ilmn)+il+1
         pawtab(itypat)%klmntomn(2,klmn)=indlmn(2,jlmn)+jl+1
         pawtab(itypat)%klmntomn(3,klmn)=indlmn(3,ilmn)
         pawtab(itypat)%klmntomn(4,klmn)=indlmn(3,jlmn)
       end do
     end do

!    ==================================================
!    B-PAW+U: overlap between atomic wavefunctions
!    --------------------------------------------------
     if(dmatpuopt==1) then
       write(message, '(4a)' ) ch10,&
&       ' pawpuxinit : dmatpuopt=1 ',ch10,&
&       '   PAW+U: dens. mat. constructed by projection on atomic wfn inside PAW augm. region(s)'
       call wrtout(std_out,message,'COLL')
       write(message, '(8a)' ) ch10,&
&       ' pawpuxinit: WARNING: Check that the first partial wave for lpawu:', ch10, &
&       '                      - Is an atomic eigenfunction  ',ch10, &
&       '                      - Is normalized ',ch10, &
&       '                      In other cases, choose dmatpuopt=2'
       call wrtout(std_out,message,'COLL')
     else if(dmatpuopt==2) then
       write(message, '(6a)' ) ch10,&
&       ' pawpuxinit : dmatpuopt=2 ',ch10,&
&       '   PAW+U: dens. mat. constructed by selecting contribution',ch10,&
&       '          for each angular momentum to the density (inside PAW augm. region(s))'
       call wrtout(std_out,message,'COLL')
     else if(dmatpuopt==3) then
       write(message, '(a,a,a,a,a,a)' ) ch10,&
&       ' pawpuxinit : dmatpuopt=3 ',ch10,&
&       '    PAW+U: dens. mat. constructed by projection on atomic wfn inside PAW augm. region(s)',ch10,&
&       '           and normalized inside PAW augm. region(s)'
       call wrtout(std_out,message,'COLL')
       write(message, '(6a)' ) ch10,&
&       ' pawpuxinit: WARNING: Check that the first partial wave for lpawu:', ch10, &
&       '                     is an atomic eigenfunction',ch10, &
&       '                     In the other case, choose dmatpuopt=2'
       call wrtout(std_out,message,'COLL')
     end if

     ABI_ALLOCATE(ff,(mesh_size))
     ff(:)=zero

     if (allocated(pawtab(itypat)%ph0phiint)) then
       ABI_DEALLOCATE(pawtab(itypat)%ph0phiint)
     end if
     if (allocated(pawtab(itypat)%zioneff)) then
       ABI_DEALLOCATE(pawtab(itypat)%zioneff)
     end if
     ABI_ALLOCATE(pawtab(itypat)%ph0phiint,(pawtab(itypat)%nproju))
     ABI_ALLOCATE(pawtab(itypat)%zioneff,(pawtab(itypat)%nproju))

     icount=0
     do iu=1,pawtab(itypat)%nproju
!      write(std_out,*)'DJA iu',iu,' mesh_size',pawtab(itypat)%mesh_size
!      do ju=2,pawtab(itypat)%mesh_size
!      ff(ju)=pawtab(itypat)%phi(ju,pawtab(itypat)%lnproju(iu))/pawrad(itypat)%rad(ju)
!      write(std_out,fmt='(i5,3e15.5)')ju,pawrad(itypat)%rad(ju),ff(ju),&
!      &         RadFnH(pawrad(itypat)%rad(ju),4,3,15.0_dp)
!      end do
!      ff(1:mesh_size)=pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(iu))**2
!      call simp_gen(int1,ff,pawrad(itypat))
!      write(std_out,*)'DJA iu',iu,'int1 ',int1
!      write(std_out,*)'DJA int1,IRadFnH',int1,IRadFnH(0.0_dp,pawrad(itypat)%rmax,4,3,12)
!      Calculation of zioneff
       ju=pawtab(itypat)%mesh_size-1
       ak=pawtab(itypat)%phi(ju,pawtab(itypat)%lnproju(iu))/pawtab(itypat)%phi(ju+1,pawtab(itypat)%lnproju(iu))
       ak=ak*(pawrad(itypat)%rad(ju+1)/pawrad(itypat)%rad(ju))**(pawtab(itypat)%lpawu-1)
       pawtab(itypat)%zioneff(iu)=log(ak)/(pawrad(itypat)%rad(ju+1)-pawrad(itypat)%rad(ju))
!      Calculation of ph0phiint
       ff(1:mesh_size)=pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(1))&
&       *pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(iu))
       call simp_gen(int1,ff,pawrad(itypat))
       pawtab(itypat)%ph0phiint(iu)=int1
     end do
     if(abs(pawprtvol)>=2) then
       do icount=1,pawtab(itypat)%nproju
         write(message, '(a,a,i2,f9.5,a)' ) ch10,&
&         '  pawpuxinit: icount, ph0phiint(icount)=',icount,pawtab(itypat)%ph0phiint(icount)
         call wrtout(std_out,message,'COLL')
         write(message, '(a,f15.5)' ) &
&         '  pawpuxinit: zioneff=',pawtab(itypat)%zioneff(icount)
         call wrtout(std_out,message,'COLL')
       end do
       write(message, '(a)' ) ch10
       call wrtout(std_out,message,'COLL')
     end if

     if (allocated(pawtab(itypat)%phiphjint)) then
       ABI_DEALLOCATE(pawtab(itypat)%phiphjint)
     end if
     ABI_ALLOCATE(pawtab(itypat)%phiphjint,(pawtab(itypat)%ij_proj))

     icount=0
     do ju=1,pawtab(itypat)%nproju
       do iu=1,ju
         icount=icount+1
         if ((dmatpuopt==1).and.(useexexch==0)) then
           pawtab(itypat)%phiphjint(icount)=pawtab(itypat)%ph0phiint(iu)*&
&           pawtab(itypat)%ph0phiint(ju)
         else if((dmatpuopt==2).or.(useexexch/=0)) then
           ff(1:mesh_size)=pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(iu))&
&           *pawtab(itypat)%phi(1:mesh_size,pawtab(itypat)%lnproju(ju))
           call simp_gen(int1,ff,pawrad(itypat))
           pawtab(itypat)%phiphjint(icount)=int1
         else if((dmatpuopt>=3).and.(useexexch==0)) then
           pawtab(itypat)%phiphjint(icount)=pawtab(itypat)%ph0phiint(iu)* &
&           pawtab(itypat)%ph0phiint(ju)/pawtab(itypat)%ph0phiint(1)**(dmatpuopt-2)
         else
           write(message, '(3a)' )&
&           '  PAW+U: dmatpuopt has a wrong value !',ch10,&
&           '  Action : change value in input file'
           MSG_ERROR(message)
         end if
       end do
     end do
     if(pawtab(itypat)%ij_proj/=icount)  then
       message = ' Error in the loop for calculating phiphjint '
       MSG_ERROR(message)
     end if
     ABI_DEALLOCATE(ff)
     if(abs(pawprtvol)>=2) then
       do icount=1,pawtab(itypat)%ij_proj
         write(message, '(a,a,i2,f9.5,a)' ) ch10,&
&         '  PAW+U: icount, phiphjint(icount)=',icount,pawtab(itypat)%phiphjint(icount)
         call wrtout(std_out,message,'COLL')
       end do
     end if
!    end if

!    ======================================================================
!    C-PAW+U: Matrix elements of coulomb interaction (see PRB vol.52 5467) [[cite:Liechenstein1995]]
!    1. angular part computed from Gaunt coefficients
!    --------------------------------------------------------------------
     if (usepawu/=0) then
       lpawu=lcur

!      a. compute F(k)
!      ---------------------------------------------
       ABI_ALLOCATE(fk,(lpawu+1))
       fk(1)=pawtab(itypat)%upawu
!      cf Slater Physical Review 165, p 665 (1968) [[cite:Slater1958]]
!      write(std_out,*) "f4of2_sla",pawtab(itypat)%f4of2_sla
       if(lpawu==0) then
         fk(1)=fk(1)
       else if(lpawu==1) then
         fk(2)=pawtab(itypat)%jpawu*5._dp
       else if(lpawu==2) then
!        f4of2=0._dp
         if(pawtab(itypat)%f4of2_sla<-0.1_dp)  then
           f4of2=0.625_dp
           pawtab(itypat)%f4of2_sla=f4of2
         else
           f4of2=pawtab(itypat)%f4of2_sla
         end if
         fk(2)=pawtab(itypat)%jpawu*14._dp/(One+f4of2)
         fk(3)=fk(2)*f4of2
!        if(abs(pawprtvol)>=2) then
         write(message,'(a,3x,a,f9.4,f9.4,f9.4,f9.4)') ch10,&
&         "Slater parameters F^0, F^2, F^4 are",fk(1),fk(2),fk(3)
         call wrtout(std_out,message,'COLL')
!        end if
       else if(lpawu==3) then
         f4of2=0.6681_dp
         f6of2=0.4943_dp
         if(pawtab(itypat)%f4of2_sla<-0.1_dp)  then
           f4of2=0.6681_dp
           pawtab(itypat)%f4of2_sla=f4of2
         else
           f4of2=pawtab(itypat)%f4of2_sla
         end if
         if(pawtab(itypat)%f6of2_sla<-0.1_dp)  then
           f6of2=0.4943_dp
           pawtab(itypat)%f6of2_sla=f6of2
         else
           f6of2=pawtab(itypat)%f6of2_sla
         end if
         fk(2)=pawtab(itypat)%jpawu*6435._dp/(286._dp+195._dp*f4of2+250._dp*f6of2)
         fk(3)=fk(2)*f4of2
         fk(4)=fk(2)*f6of2
         write(std_out,'(a,3x,a,f9.4,f9.4,f9.4,f9.4)') ch10,&
&         "Slater parameters F^0, F^2, F^4, F^6 are",fk(1),fk(2),fk(3),fk(4)
       else
         write(message, '(a,i0,2a)' )&
&         ' lpawu=',lpawu,ch10,&
&         ' lpawu not equal to 0 ,1 ,2 or 3 is not allowed'
         MSG_ERROR(message)
       end if

!      b. Compute ak and vee.
!      ---------------------------------------------
       if (allocated(pawtab(itypat)%vee)) then
         ABI_DEALLOCATE(pawtab(itypat)%vee)
       end if
       sz1=2*lpawu+1
       ABI_ALLOCATE(pawtab(itypat)%vee,(sz1,sz1,sz1,sz1))
       pawtab(itypat)%vee=zero
       lmpawu=(lpawu-1)**2+2*(lpawu-1)+1  ! number of m value below correlated orbitals
       klm0u=lmpawu*(lmpawu+1)/2          ! value of klmn just below correlated orbitals
!      --------- 4 loops for interaction matrix
       do m1=-lpawu,lpawu
         m11=m1+lpawu+1
         do m2=-lpawu,m1
           m21=m2+lpawu+1
!          klma= number of pair before correlated orbitals +
!          number of pair for m1 lower than correlated orbitals
!          (m1+lpawu+1)*(lpawu-1) + number of pairs for correlated orbitals
!          before (m1,m2) + number of pair for m2 lower than current value
           klma=klm0u+m11*lmpawu+(m11-1)*m11/2+m21
           do m3=-lpawu,lpawu
             m31=m3+lpawu+1
             do m4=-lpawu,m3
               m41=m4+lpawu+1
               klmb=klm0u+m31*lmpawu+(m31-1)*m31/2+m41
!              --------- loop on k=1,2,3 (4 if f orbitals)
               do kyc=1,2*lpawu+1,2
                 lkyc=kyc-1
                 lmkyc=(lkyc+1)*(lkyc)+1
                 ak=zero
                 do mkyc=-lkyc,lkyc,1
                   isela=pawang%gntselect(lmkyc+mkyc,klma)
                   iselb=pawang%gntselect(lmkyc+mkyc,klmb)
                   if (isela>0.and.iselb>0) ak=ak +pawang%realgnt(isela)*pawang%realgnt(iselb)
                 end do
!                ----- end loop on k=1,2,3 (4 if f orbitals)
                 ak=ak/(two*dble(lkyc)+one)
                 pawtab(itypat)%vee(m11,m31,m21,m41)=ak*fk(lkyc/2+1)+pawtab(itypat)%vee(m11,m31,m21,m41)
               end do  !kyc
               pawtab(itypat)%vee(m11,m31,m21,m41)=pawtab(itypat)%vee(m11,m31,m21,m41)*four_pi
               pawtab(itypat)%vee(m21,m31,m11,m41)=pawtab(itypat)%vee(m11,m31,m21,m41)
               pawtab(itypat)%vee(m11,m41,m21,m31)=pawtab(itypat)%vee(m11,m31,m21,m41)
               pawtab(itypat)%vee(m21,m41,m11,m31)=pawtab(itypat)%vee(m11,m31,m21,m41)
             end do
           end do
         end do
       end do
       ABI_DEALLOCATE(fk)

     !  testu=0
     !  write(std_out,*) " Matrix of interaction vee(m1,m2,m1,m2)"
     !  do m1=1,2*lpawu+1
     !    write(std_out,'(2x,14(f12.6,2x))') (pawtab(itypat)%vee(m1,m2,m1,m2),m2=1,2*lpawu+1)
     !    do m2=1,2*lpawu+1
     !      testu=testu+ pawtab(itypat)%vee(m1,m2,m1,m2)
     !   enddo
     !  enddo
     !  testu=testu/((two*lpawu+one)**2)
     !  write(std_out,*) "------------------------"
     !  write(std_out,'(a,f12.6)') " U=", testu
     !  write(std_out,*) "------------------------"
     !  write(std_out,*) " Matrix of interaction vee(m1,m2,m1,m2)-vee(m1,m2,m2,m1)"
     !  do m1=1,2*lpawu+1
     !    write(std_out,'(2x,14(f12.6,2x))') ((pawtab(itypat)%vee(m1,m2,m1,m2)-pawtab(itypat)%vee(m1,m2,m2,m1)),m2=1,2*lpawu+1)
     !    do m2=1,2*lpawu+1
     !    if(m1/=m2) testumj=testumj+ pawtab(itypat)%vee(m1,m2,m1,m2)-pawtab(itypat)%vee(m1,m2,m2,m1)
     !   enddo
     !  enddo
     !  testumj=testumj/((two*lpawu)*(two*lpawu+one))
     !  write(std_out,*) "------------------------"
     !  write(std_out,'(a,f12.6)') " U-J=", testumj
     !  write(std_out,*) "------------------------"
     !  write(std_out,*) "------------------------"
     !  write(std_out,'(a,f12.6)')  " J=", testu-testumj
     !  write(std_out,*) "------------------------"

!      c. For DFT (or with exp. values usepawu=-1 or -2), compute euijkl
!      ---------------------------------------------
       compute_euijkl=(is_dfpt.or.usepawu<0)
       if (compute_euijkl) then
         if (allocated(pawtab(itypat)%euijkl)) then
           ABI_DEALLOCATE(pawtab(itypat)%euijkl)
         end if
         ABI_ALLOCATE(pawtab(itypat)%euijkl,(2,2,lmn_size,lmn_size,lmn_size,lmn_size))
         pawtab(itypat)%euijkl = zero
         compute_euij_fll = .false.
         euijkl_temp2=zero
         if (abs(usepawu)==1.or.abs(usepawu)==4) then ! Only for FLL
           if (allocated(pawtab(itypat)%euij_fll)) then ! allocate euij_fll for FLL
             ABI_DEALLOCATE(pawtab(itypat)%euij_fll)
           end if
           ABI_ALLOCATE(pawtab(itypat)%euij_fll,(lmn2_size))
           pawtab(itypat)%euij_fll = zero
           compute_euij_fll = .true.
         end if

!        loop on i,j
         do klmna=1,lmn2_size
           ilmn=pawtab(itypat)%indklmn(7,klmna) ! i
           jlmn=pawtab(itypat)%indklmn(8,klmna) ! j
           if (pawtab(itypat)%indlmn(1,ilmn)==lpawu.and.pawtab(itypat)%indlmn(1,jlmn)==lpawu) then ! only correlated orbitals
             iu = pawtab(itypat)%indlmn(3,ilmn) ! ni
             ju = pawtab(itypat)%indlmn(3,jlmn) ! nj
             phiint_ij = pawtab(itypat)%phiphjint(iu+(ju*(ju-1))/2) ! iu <= ju by construction (ilmn<=jlmn)
             m2 = pawtab(itypat)%indlmn(2,ilmn) ! mi
             m21=m2+lpawu+1
             m1 = pawtab(itypat)%indlmn(2,jlmn) ! mj
             m11=m1+lpawu+1

             if (compute_euij_fll.and.m1==m2) then ! FLL
               pawtab(itypat)%euij_fll(klmna) = - half * phiint_ij * ( pawtab(itypat)%jpawu - pawtab(itypat)%upawu )
             end if

!            loop on ip,jp (=k,l)
             do klmnb=1,lmn2_size
               ilmnp=pawtab(itypat)%indklmn(7,klmnb) ! ip (=k)
               jlmnp=pawtab(itypat)%indklmn(8,klmnb) ! jp (=l)
               if (pawtab(itypat)%indlmn(1,ilmnp)==lpawu.and.pawtab(itypat)%indlmn(1,jlmnp)==lpawu) then ! correlated orbitals
                 iup = pawtab(itypat)%indlmn(3,ilmnp) ! nip
                 jup = pawtab(itypat)%indlmn(3,jlmnp) ! njp
                 phiint_ipjp = pawtab(itypat)%phiphjint(iup+(jup*(jup-1))/2) ! iup <= jup by construction (ilmnp<=jlmnp)
                 m4 = pawtab(itypat)%indlmn(2,ilmnp) ! mip
                 m41=m4+lpawu+1
                 m3 = pawtab(itypat)%indlmn(2,jlmnp) ! mjp
                 m31=m3+lpawu+1

                 euijkl_dc = zero
!                Compute the double-counting part of euijkl (invariant when exchanging i<-->j or ip<-->jp)
                 if (m1==m2.and.m3==m4) then ! In that case, we have to add the double-counting term

                     do sig=1,2
                       do sigp=1,2

                         if (abs(usepawu)==1.or.abs(usepawu)==4) then ! FLL

                           if (sig==sigp) then
                             euijkl_dc(sig,sigp) = &
&                             phiint_ij * phiint_ipjp * ( pawtab(itypat)%upawu - pawtab(itypat)%jpawu )
                           else
                             euijkl_dc(sig,sigp) = &
&                             phiint_ij * phiint_ipjp * pawtab(itypat)%upawu
                           end if

                         else if (abs(usepawu)==2) then ! AMF

                           if (sig==sigp) then
                             euijkl_dc(sig,sigp) = &
&                             two*lpawu/(two*lpawu+one) * phiint_ij * phiint_ipjp * ( pawtab(itypat)%upawu - pawtab(itypat)%jpawu )
                           else
                             euijkl_dc(sig,sigp) = &
&                             phiint_ij * phiint_ipjp * pawtab(itypat)%upawu
                           end if

                         end if

                       end do ! sigp
                     end do ! sig

                 end if ! double-counting term

                 vee1 = pawtab(itypat)%vee(m11,m31,m21,m41)
!                Note : vee(13|24) = vee(23|14) ( so : i    <--> j     )
!                       vee(13|24) = vee(14|23) ( so : ip   <--> jp    )
!                       vee(13|24) = vee(24|13) ( so : i,ip <--> j,jp  )
!                Also : vee(13|24) = vee(31|42) ( so : i,j  <--> ip,jp )
!                ==> vee1 is invariant with respect to the permutations i <--> j , ip <--> jp and i,ip <--> j,jp
!                ( The term 'phiint_ij * phiint_ipjp' has the same properties)
                 do sig=1,2
                   do sigp=1,2
                     euijkl_temp(sig,sigp) = phiint_ij * phiint_ipjp * vee1
                   end do
                 end do

                 vee2 = pawtab(itypat)%vee(m11,m31,m41,m21)
!                Note : vee(13|42) = vee(43|12) ( so : ip   <--> j     )
!                       vee(13|42) = vee(12|43) ( so : i    <--> jp    )
!                       vee(13|42) = vee(42|13) ( so : i,ip <--> jp,j  )
!                Also : vee(13|42) = vee(31|24) ( so : i,j  <--> ip,jp )
!                Combining the third and fourth rule we get:
!                       vee(13|42) = vee(42|13) = vee(24|31) ( so : i,ip  <--> j,jp )
!                ==> vee2 is invariant only with respect to the permutation i,ip <--> j,jp

!                Terms i,j,ip,jp (m2,m1,m4,m3) and j,i,jp,ip (m1,m2,m3,m4)
                 do sig=1,2
                   euijkl_temp2(sig,sig) = phiint_ij * phiint_ipjp * vee2
                 end do
                 pawtab(itypat)%euijkl(:,:,ilmn,jlmn,ilmnp,jlmnp) = euijkl_temp(:,:) - euijkl_temp2(:,:) - euijkl_dc(:,:)
                 pawtab(itypat)%euijkl(:,:,jlmn,ilmn,jlmnp,ilmnp) = pawtab(itypat)%euijkl(:,:,ilmn,jlmn,ilmnp,jlmnp)

!                Term j,i,ip,jp (m1,m2,m4,m3)
                 vee2 = pawtab(itypat)%vee(m21,m31,m41,m11)
                 do sig=1,2
                   euijkl_temp2(sig,sig) = phiint_ij * phiint_ipjp * vee2
                 end do
                 pawtab(itypat)%euijkl(:,:,jlmn,ilmn,ilmnp,jlmnp) = euijkl_temp(:,:) - euijkl_temp2(:,:) - euijkl_dc(:,:)

!                Term i,j,jp,ip (m2,m1,m3,m4)
                 vee2 = pawtab(itypat)%vee(m11,m41,m31,m21)
                 do sig=1,2
                   euijkl_temp2(sig,sig) = phiint_ij * phiint_ipjp * vee2
                 end do
                 pawtab(itypat)%euijkl(:,:,ilmn,jlmn,jlmnp,ilmnp) = euijkl_temp(:,:) - euijkl_temp2(:,:) - euijkl_dc(:,:)

               end if ! correlated orbitals
             end do ! klmnb
           end if ! correlated orbitals
         end do ! klmna

       end if ! compute_euijkl
     end if ! usepawu

!    ======================================================================
!    D-Local ex-exchange: Matrix elements of coulomb interaction and Fk
!    ----------------------------------------------------------------------
     if (useexexch/=0) then
       lexexch=lcur

!      a. compute F(k)
!      ---------------------------------------------
       if (allocated(pawtab(itypat)%fk)) then
         ABI_DEALLOCATE(pawtab(itypat)%fk)
       end if
       ABI_ALLOCATE(pawtab(itypat)%fk,(6,4))
       pawtab(itypat)%fk=zero
       ABI_ALLOCATE(ff,(mesh_size))
       ABI_ALLOCATE(gg,(mesh_size))
       ff(:)=zero;gg(:)=zero
       kln=(pawtab(itypat)%lnproju(1)*( pawtab(itypat)%lnproju(1)+1)/2)
       do ll=1,lexexch+1
         ll1=2*ll-2
         if (int_meshsz<mesh_size) ff(int_meshsz+1:mesh_size)=zero
         ff(1:int_meshsz)=pawtab(itypat)%phiphj(1:int_meshsz,kln)
         call poisson(ff,ll1,pawrad(itypat),gg)
         ff(1)=zero
         ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln)*gg(2:mesh_size))&
&         /pawrad(itypat)%rad(2:mesh_size)
         call simp_gen(intg,ff,pawrad(itypat))
         pawtab(itypat)%fk(1,ll)=intg*(two*ll1+one)
       end do
       if (pawtab(itypat)%nproju==2) then
         kln1=kln+pawtab(itypat)%lnproju(1)
         kln2=kln1+1
         do ll=1,lexexch+1
           ll1=2*ll-2
           if (int_meshsz<mesh_size) ff(int_meshsz+1:mesh_size)=zero
           ff(1:int_meshsz)=pawtab(itypat)%phiphj(1:int_meshsz,kln1)
           call poisson(ff,ll1,pawrad(itypat),gg)
           ff(1)=zero
           ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln1)*gg(2:mesh_size))&
&           /pawrad(itypat)%rad(2:mesh_size)
           call simp_gen(intg,ff,pawrad(itypat))
           pawtab(itypat)%fk(2,ll)=intg*(two*ll1+one)
         end do
         do ll=1,lexexch+1
           ll1=2*ll-2
           if (int_meshsz<mesh_size) ff(int_meshsz+1:mesh_size)=zero
           ff(1:int_meshsz)=pawtab(itypat)%phiphj(1:int_meshsz,kln2)
           call poisson(ff,ll1,pawrad(itypat),gg)
           ff(1)=zero
           ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln2)*gg(2:mesh_size))&
&           /pawrad(itypat)%rad(2:mesh_size)
           call simp_gen(intg,ff,pawrad(itypat))
           pawtab(itypat)%fk(3,ll)=intg*(two*ll1+one)
         end do
         do ll=1,lexexch+1
           ll1=2*ll-2
           if (int_meshsz<mesh_size) ff(int_meshsz+1:mesh_size)=zero
           ff(1:int_meshsz)=pawtab(itypat)%phiphj(1:int_meshsz,kln)
           call poisson(ff,ll1,pawrad(itypat),gg)
           ff(1)=zero
           ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln1)*gg(2:mesh_size))&
&           /pawrad(itypat)%rad(2:mesh_size)
           call simp_gen(intg,ff,pawrad(itypat))
           pawtab(itypat)%fk(4,ll)=intg*(two*ll1+one)
         end do
         do ll=1,lexexch+1
           ll1=2*ll-2
           if (int_meshsz<mesh_size) ff(int_meshsz+1:mesh_size)=zero
           ff(1:int_meshsz)=pawtab(itypat)%phiphj(1:int_meshsz,kln)
           call poisson(ff,ll1,pawrad(itypat),gg)
           ff(1)=zero
           ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln2)*gg(2:mesh_size))&
&           /pawrad(itypat)%rad(2:mesh_size)
           call simp_gen(intg,ff,pawrad(itypat))
           pawtab(itypat)%fk(5,ll)=intg*(two*ll1+one)
         end do
         do ll=1,lexexch+1
           ll1=2*ll-2
           if (int_meshsz<mesh_size) ff(int_meshsz+1:mesh_size)=zero
           ff(1:int_meshsz)=pawtab(itypat)%phiphj(1:int_meshsz,kln1)
           call poisson(ff,ll1,pawrad(itypat),gg)
           ff(1)=zero
           ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln2)*gg(2:mesh_size))&
&           /pawrad(itypat)%rad(2:mesh_size)
           call simp_gen(intg,ff,pawrad(itypat))
           pawtab(itypat)%fk(6,ll)=intg*(two*ll1+one)
         end do
         f4of2=0.6681_dp
         f6of2=0.4943_dp
       end if
       ABI_DEALLOCATE(ff)
       ABI_DEALLOCATE(gg)

!      b. Compute vex.
!      ---------------------------------------------
       if (allocated(pawtab(itypat)%vex)) then
         ABI_DEALLOCATE(pawtab(itypat)%vex)
       end if
       sz1=2*lexexch+1
       ABI_ALLOCATE(pawtab(itypat)%vex,(sz1,sz1,sz1,sz1,4))
       pawtab(itypat)%vex=zero
       lmexexch=(lexexch-1)**2+2*(lexexch-1)+1  ! number of m value below correlated orbitals
       klm0x=lmexexch*(lmexexch+1)/2            ! value of klmn just below correlated orbitals
!      --------- 4 loops for interaction matrix
       do m1=-lexexch,lexexch
         m11=m1+lexexch+1
         do m2=-lexexch,m1
           m21=m2+lexexch+1
!          klma= number of pair before correlated orbitals +
!          number of pair for m1 lower than correlated orbitals
!          (m1+lexexch+1)*(lexexch-1) + number of pairs for correlated orbitals
!          before (m1,m2) + number of pair for m2 lower than current value
           klma=klm0x+m11*lmexexch+(m11-1)*m11/2+m21
           do m3=-lexexch,lexexch
             m31=m3+lexexch+1
             do m4=-lexexch,m3
               m41=m4+lexexch+1
               klmb=klm0x+m31*lmexexch+(m31-1)*m31/2+m41
!              --------- loop on k=1,2,3 (4 if f orbitals)
               do kyc=1,2*lexexch+1,2
                 lkyc=kyc-1
                 ll=(kyc+1)/2
                 lmkyc=(lkyc+1)*(lkyc)+1
                 ak=zero
                 do mkyc=-lkyc,lkyc,1
                   isela=pawang%gntselect(lmkyc+mkyc,klma)
                   iselb=pawang%gntselect(lmkyc+mkyc,klmb)
                   if (isela>0.and.iselb>0) ak=ak +pawang%realgnt(isela)*pawang%realgnt(iselb)
                 end do
!                ----- end loop on k=1,2,3 (4 if f orbitals)
                 pawtab(itypat)%vex(m11,m31,m21,m41,ll)=ak/(two*dble(lkyc)+one)
               end do  !kyc
               do ll=1,4
                 pawtab(itypat)%vex(m11,m31,m21,m41,ll)=pawtab(itypat)%vex(m11,m31,m21,m41,ll)*four_pi
                 pawtab(itypat)%vex(m21,m31,m11,m41,ll)=pawtab(itypat)%vex(m11,m31,m21,m41,ll)
                 pawtab(itypat)%vex(m11,m41,m21,m31,ll)=pawtab(itypat)%vex(m11,m31,m21,m41,ll)
                 pawtab(itypat)%vex(m21,m41,m11,m31,ll)=pawtab(itypat)%vex(m11,m31,m21,m41,ll)
               end do
             end do
           end do
         end do
       end do

     end if !useexexch/=0

     if (present(ucrpa)) then
       if (ucrpa>=1) then
         call calc_ubare(itypat,lcur,pawang,pawrad(itypat),pawtab(itypat))
         call calc_ubare(itypat,lcur,pawang,pawrad(itypat),pawtab(itypat),pawtab(itypat)%rpaw)
       end if
     end if
   end if !lcur/=-1
 end do !end loop on typat

 DBG_EXIT("COLL")

 end subroutine pawpuxinit
!!***

!----------------------------------------------------------------------

!!****f* m_paw_correlations/pawuenergy
!! NAME
!! pawuenergy
!!
!! FUNCTION
!! Compute contributions to energy for PAW+U calculations
!!
!! INPUTS
!!  iatom=index of current atom (absolute index, the index on current proc)
!!  noccmmp(2*lpawu+1,2*lpawu+1,nspden)=density matrix in the PAW augm. region
!!  nocctot(nspden)=number of electrons in the correlated subspace
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawtab <type(pawtab_type)>=paw tabulated starting data:
!!     %lpawu=l used for lda+u
!!     %vee(2*lpawu+1*4)=screened coulomb matrix
!!  dmft_dc,e_ee,e_dc,e_dcdc,u_dmft,j_dmft= optional arguments for DMFT
!!
!! OUTPUT
!!  eldaumdc= PAW+U contribution to total energy
!!  eldaumdcdc= PAW+U contribution to double-counting total energy
!!
!! PARENTS
!!      m_energy,pawdenpot
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

 subroutine pawuenergy(iatom,eldaumdc,eldaumdcdc,noccmmp,nocctot,pawprtvol,pawtab,&
 &                     dmft_dc,e_ee,e_dc,e_dcdc,u_dmft,j_dmft) ! optional arguments (DMFT)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: iatom,pawprtvol
 integer,optional,intent(in) :: dmft_dc
 real(dp),intent(in) :: noccmmp(:,:,:,:),nocctot(:)
 real(dp),intent(inout) :: eldaumdc,eldaumdcdc
 real(dp),optional,intent(inout) :: e_ee,e_dc,e_dcdc
 real(dp),optional,intent(in) :: j_dmft,u_dmft
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
!Option for interaction energy in case of non-collinear magnetism:
!           1: E_int=-U/4.N.(N-2)
!           2: E_int=-U/2.(Nup.(Nup-1)+Ndn.(Ndn-1))
 integer,parameter :: option_interaction=1
 integer :: cplex_occ,dmftdc,ispden,jspden,lpawu,m1,m11,m2,m21,m3,m31,m4,m41,nspden
 real(dp) :: eks_opt3,edcdc_opt3,edcdctemp,edctemp,eldautemp,jpawu,jpawu_dc,mnorm,mx,my,mz
 real(dp) :: n_sig,n_sigs,n_msig,n_msigs,n_dndn,n_tot,n_upup
 real(dp) :: n12_ud_im,n12_du_im
 real(dp) :: n12_ud_re,n12_du_re
 real(dp) :: n34_ud_im,n34_du_im
 real(dp) :: n34_ud_re,n34_du_re
 real(dp) :: upawu
 real(dp),allocatable :: n12_sig(:),n34_msig(:),n34_sig(:)
 character(len=500) :: message

! *****************************************************

 nspden=size(nocctot)
 cplex_occ=size(noccmmp,1)

 if (size(noccmmp,4)/=nspden) then
   message='size of nocctot and noccmmp are inconsistent!'
   MSG_BUG(message)
 end if
 if (pawtab%usepawu<0) then
   message='not allowed for usepawu<0!'
   MSG_BUG(message)
 end if
 if(present(dmft_dc))  then
   dmftdc=dmft_dc
   if(pawtab%usepawu<10) then
     write(message,'(a,i5)') "usepawu should be =10 if dmft_dc is present ",pawtab%usepawu
     MSG_BUG(message)
   end if
 else
   dmftdc=0
 end if

 DBG_ENTER("COLL")

 lpawu=pawtab%lpawu
 upawu=pawtab%upawu;if(present(u_dmft)) upawu=u_dmft
 jpawu=pawtab%jpawu;if(present(j_dmft)) jpawu=j_dmft

!======================================================
!Compute LDA+U Energy
!-----------------------------------------------------

 eldautemp=zero
 edcdc_opt3=zero
 eks_opt3=zero

 ABI_ALLOCATE(n12_sig,(cplex_occ))
 ABI_ALLOCATE(n34_msig,(cplex_occ))
 ABI_ALLOCATE(n34_sig,(cplex_occ))
 do ispden=1,min(nspden,2)
   jspden=min(nspden,2)-ispden+1

!  Compute n_sigs and n_msigs for pawtab%usepawu=3
   if (nspden<=2) then
     n_sig =nocctot(ispden)
     n_msig=nocctot(jspden)
     n_tot=n_sig+n_msig
   else
     n_tot=nocctot(1)
     mx=nocctot(2)
     my=nocctot(3)
     mz=nocctot(4)
     mnorm=sqrt(mx*mx+my*my+mz*mz)
     if (ispden==1) then
!      n_sig =half*(n_tot+mnorm)
!      n_msig=half*(n_tot-mnorm)
       n_sig =half*(n_tot+sign(mnorm,mz))
       n_msig=half*(n_tot-sign(mnorm,mz))
     else
!      n_sig =half*(n_tot-mnorm)
!      n_msig=half*(n_tot+mnorm)
       n_sig =half*(n_tot-sign(mnorm,mz))
       n_msig=half*(n_tot+sign(mnorm,mz))
     end if
   end if
   n_sigs =n_sig/(float(2*lpawu+1))
   n_msigs =n_msig/(float(2*lpawu+1))
!  if(pawtab%usepawu==3) then
!    write(message,fmt=12) "noccmmp11 ",ispden,noccmmp(1,1,1,ispden)
!    call wrtout(std_out,message,'COLL')
!    write(message,fmt=12) "noccmmp11 ",jspden,noccmmp(1,1,1,jspden)
!    call wrtout(std_out,message,'COLL')
!    write(message,fmt=12) "n_sig      ",ispden,n_sig
!    call wrtout(std_out,message,'COLL')
!    write(message,fmt=12) "n_msig     ",jspden,n_msig
!    call wrtout(std_out,message,'COLL')
!    write(message,fmt=12) "n_sigs     ",ispden,n_sigs
!    call wrtout(std_out,message,'COLL')
!    write(message,fmt=12) "n_msigs    ",jspden,n_msigs
!    call wrtout(std_out,message,'COLL')
!  endif
!  12 format(a,i4,e20.10)

!  Compute interaction energy E_{ee}
   do m1=-lpawu,lpawu
     m11=m1+lpawu+1
     do m2=-lpawu,lpawu
       m21=m2+lpawu+1
       n12_sig(:)=noccmmp(:,m11,m21,ispden)
       if(m21==m11.and.(pawtab%usepawu==3.or.dmftdc==3)) n12_sig(1)=n12_sig(1)-n_sigs
       do m3=-lpawu,lpawu
         m31=m3+lpawu+1
         do m4=-lpawu,lpawu
           m41=m4+lpawu+1
           n34_sig(:) =noccmmp(:,m31,m41,ispden)
           n34_msig(:)=noccmmp(:,m31,m41,jspden)
           if(m31==m41.and.(pawtab%usepawu==3.or.dmftdc==3)) then
             n34_sig(1)= n34_sig(1) - n_sigs
             n34_msig(1)= n34_msig(1) - n_msigs
           end if
           eldautemp=eldautemp &
&           + n12_sig(1)*n34_msig(1)*pawtab%vee(m11,m31,m21,m41) &
&           + n12_sig(1)*n34_sig(1) *(pawtab%vee(m11,m31,m21,m41)-pawtab%vee(m11,m31,m41,m21))
           if(cplex_occ==2) then
             eldautemp=eldautemp &
&             - n12_sig(2)*n34_msig(2)*pawtab%vee(m11,m31,m21,m41) &
&             - n12_sig(2)*n34_sig(2) *(pawtab%vee(m11,m31,m21,m41)-pawtab%vee(m11,m31,m41,m21))
           end if
           if (pawtab%usepawu==3.or.dmftdc==3) then
             edcdc_opt3=edcdc_opt3 &
&             + n_sigs*n34_msig(1)*pawtab%vee(m11,m31,m21,m41) &
&             + n_sigs*n34_sig(1) *(pawtab%vee(m11,m31,m21,m41)-pawtab%vee(m11,m31,m41,m21))
             eks_opt3=eks_opt3 &
&             + noccmmp(1,m11,m21,ispden)*n34_msig(1)*pawtab%vee(m11,m31,m21,m41) &
&             + noccmmp(1,m11,m21,ispden)*n34_sig(1) *(pawtab%vee(m11,m31,m21,m41)-pawtab%vee(m11,m31,m41,m21))
             if(cplex_occ==2) then
               eks_opt3=eks_opt3 &
&               - noccmmp(2,m11,m21,ispden)*n34_msig(2)*pawtab%vee(m11,m31,m21,m41) &
&               - noccmmp(2,m11,m21,ispden)*n34_sig(2) *(pawtab%vee(m11,m31,m21,m41)-pawtab%vee(m11,m31,m41,m21))
             end if
           end if
         end do ! m4
       end do ! m3
     end do ! m2
   end do ! m1

 end do ! ispden
 if (nspden==1) eldautemp=two*eldautemp ! Non-magn. system: sum up and dn energies
 ABI_DEALLOCATE(n12_sig)
 ABI_DEALLOCATE(n34_msig)
 ABI_DEALLOCATE(n34_sig)

!Non-collinear magnetism: add non-diagonal term; see (Eq 3) in PRB 72, 024458 (2005) [[cite:Shurikov2005]]
 if (nspden==4) then
   do m1=-lpawu,lpawu
     m11=m1+lpawu+1
     do m2=-lpawu,lpawu
       m21=m2+lpawu+1
       n12_ud_re=noccmmp(1,m11,m21,3) ! updn
       n12_ud_im=noccmmp(2,m11,m21,3) ! updn
       n12_du_re=noccmmp(1,m11,m21,4) ! dnup
       n12_du_im=noccmmp(2,m11,m21,4) ! dnup
       do m3=-lpawu,lpawu
         m31=m3+lpawu+1
         do m4=-lpawu,lpawu
           m41=m4+lpawu+1
           n34_ud_re=noccmmp(1,m31,m41,3)  ! updn
           n34_ud_im=noccmmp(2,m31,m41,3)  ! updn
           n34_du_re=noccmmp(1,m31,m41,4)  ! dnup
           n34_du_im=noccmmp(2,m31,m41,4)  ! dnup
           eldautemp=eldautemp-pawtab%vee(m11,m31,m41,m21) &
&           *(n12_ud_re*n34_du_re-n12_ud_im*n34_du_im &
&           +n12_du_re*n34_ud_re-n12_du_im*n34_ud_im)
           if (pawtab%usepawu==3.or.dmftdc==3) then
             eks_opt3=eks_opt3-pawtab%vee(m11,m31,m41,m21) &
&             *(n12_ud_re*n34_du_re-n12_ud_im*n34_du_im &
&             +n12_du_re*n34_ud_re-n12_du_im*n34_ud_im)
           end if
         end do ! m4
       end do ! m3
     end do ! m2
   end do ! m1
 end if

!Divide eldautemp by 2; see (Eq 1) in PRB 77, 155104 (2008) [[cite:Amadon2008a]]
 eldautemp=half*eldautemp

!if (nspden==1) then
!n_tot=two*nocctot(1)
!n_upup=nocctot(1)
!n_dndn=nocctot(1)
!else if (nspden==2) then
!n_tot=nocctot(1)+nocctot(2)
!n_upup=nocctot(1)
!n_dndn=nocctot(2)
!else if (nspden==4) then
!n_tot=nocctot(1)
!mx=nocctot(2)
!my=nocctot(3)
!mz=nocctot(4)
!mnorm=sqrt(mx*mx+my*my+mz*mz)
!n_upup=half*(n_tot+mnorm)
!n_dndn=half*(n_tot-mnorm)
!end if
 n_upup=n_sig
 n_dndn=n_msig

 edcdctemp=zero;edctemp=zero

!Full localized limit
 if((pawtab%usepawu==1.or.pawtab%usepawu==4).or.(dmftdc==1.or.dmftdc==4.or.dmftdc==5)) then
   jpawu_dc=jpawu
   if(dmftdc==4)  then
     jpawu_dc=zero
   end if
   edcdctemp=edcdctemp-half*upawu*n_tot**2
   edctemp  =edctemp  +half*upawu*(n_tot*(n_tot-one))
   if (nspden/=4.or.option_interaction==2) then
     if(dmftdc/=5.and.pawtab%usepawu/=4) then
       edcdctemp=edcdctemp+half*jpawu_dc*(n_upup**2+n_dndn**2)
       edctemp  =edctemp  -half*jpawu_dc*(n_upup*(n_upup-one)+n_dndn*(n_dndn-one))
     else if(dmftdc==5.or.pawtab%usepawu==4)  then
       edcdctemp=edcdctemp+quarter*jpawu_dc*n_tot**2
       edctemp  =edctemp  -quarter*jpawu_dc*(n_tot*(n_tot-two))
     end if
   else if (nspden==4.and.option_interaction==1) then
!    write(message,'(a)') "  warning: option_interaction==1 for test         "
!    call wrtout(std_out,message,'COLL')
     edcdctemp=edcdctemp+quarter*jpawu_dc*n_tot**2
     edctemp  =edctemp  -quarter*jpawu_dc*(n_tot*(n_tot-two))
   else if (nspden==4.and.option_interaction==3) then
!    edcdctemp= \frac{J}/{4}[ N(N) + \vect{m}.\vect{m}]
     edcdctemp=edcdctemp+quarter*jpawu_dc*(n_tot**2 + &
&     mx**2+my**2+mz**2)  ! +\frac{J}/{4}\vect{m}.\vect{m}
!    edctemp= -\frac{J}/{4}[ N(N-2) + \vect{m}.\vect{m}]
     edctemp  =edctemp  -quarter*jpawu_dc*(  &
&     (n_tot*(n_tot-two)) +   &
&     mx**2+my**2+mz**2)  ! -\frac{J}/{4}\vect{m}.\vect{m}
   end if

!  Around mean field
 else if(pawtab%usepawu==2.or.dmftdc==2) then
   edctemp=edctemp+upawu*(n_upup*n_dndn)&
&   +half*(upawu-jpawu)*(n_upup**2+n_dndn**2) &
&   *(dble(2*lpawu)/dble(2*lpawu+1))
   edcdctemp=-edctemp
 else if(pawtab%usepawu==3.or.dmftdc==3) then
   edcdctemp=edcdc_opt3
   if(abs(pawprtvol)>=3) then
     write(message,fmt=11) "edcdc_opt3          ",edcdc_opt3
     call wrtout(std_out,message,'COLL')
     write(message,fmt=11) "eks_opt3            ",eks_opt3
     call wrtout(std_out,message,'COLL')
     write(message,fmt=11) "eks+edcdc_opt3      ",eks_opt3+edcdc_opt3
     call wrtout(std_out,message,'COLL')
     write(message,fmt=11) "(eks+edcdc_opt3)/2  ",(eks_opt3+edcdc_opt3)/2.d0
     call wrtout(std_out,message,'COLL')
   end if
 end if

 eldaumdc  =eldaumdc  +eldautemp-edctemp
 eldaumdcdc=eldaumdcdc-eldautemp-edcdctemp

!if(pawtab%usepawu/=10.or.pawprtvol>=3) then
 if(abs(pawprtvol)>=3) then
   if(pawtab%usepawu<10) then
     write(message, '(5a,i4)')ch10,'======= LDA+U Energy terms (in Hartree) ====',ch10,&
&     ch10,' For Atom ',iatom
   else if (pawtab%usepawu >= 10) then
     write(message, '(5a,i4)')ch10,'  ===   LDA+U Energy terms for the DMFT occupation matrix ==',ch10,&
&     ch10,' For Atom ',iatom
   end if

   call wrtout(std_out,message,'COLL')
   write(message, '(a)' )"   Contributions to the direct expression of energy:"
   call wrtout(std_out,  message,'COLL')
   write(message,fmt=11) "     Double counting  correction   =",edctemp
   call wrtout(std_out,  message,'COLL')
   write(message,fmt=11) "     Interaction energy            =",eldautemp
   call wrtout(std_out,  message,'COLL')
   write(message,fmt=11) "     Total LDA+U Contribution      =",eldautemp-edctemp
   call wrtout(std_out,  message,'COLL')
   write(message, '(a)' )' '
   call wrtout(std_out,  message,'COLL')
   write(message, '(a)' )"   For the ""Double-counting"" decomposition:"
   call wrtout(std_out,  message,'COLL')
   write(message,fmt=11) "     LDA+U Contribution            =",-eldautemp-edcdctemp
   call wrtout(std_out,  message,'COLL')
   11 format(a,e20.10)
   if(abs(pawprtvol)>=2) then
     write(message,fmt=11)"     edcdctemp                     =",edcdctemp
     call wrtout(std_out,  message,'COLL')
     write(message,fmt=11)"     eldaumdcdc for current atom   =",-eldautemp-edcdctemp
     call wrtout(std_out,  message,'COLL')
     write(message, '(a)' )' '
     call wrtout(std_out,  message,'COLL')
     write(message,fmt=11)"   pawuenergy: -VUKS pred          =",eldaumdcdc-eldaumdc
     call wrtout(std_out,  message,'COLL')
   end if
   write(message, '(a)' )' '
   call wrtout(std_out,  message,'COLL')
 end if

!For DMFT calculation
 if(present(e_ee))   e_ee=e_ee+eldautemp
 if(present(e_dc))   e_dc=e_dc+edctemp
 if(present(e_dcdc)) e_dcdc=e_dcdc+edcdctemp

 DBG_EXIT("COLL")

 end subroutine pawuenergy
!!***

!----------------------------------------------------------------------

!!****f* m_paw_correlations/pawxenergy
!! NAME
!! pawxenergy
!!
!! FUNCTION
!! Compute contributions to energy for PAW+ local exact exchange calculations
!!
!! INPUTS
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrhoij <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab <type(pawtab_type)>=paw tabulated starting data:
!!     %lexexch=l used for local exact-exchange
!!     %vex(2*lexexch+1*4)=screened coulomb matrix
!!
!! SIDE EFFECTS
!!  eexex=energy is updated with the contribution of the cuyrrent atom
!!
!! PARENTS
!!      pawdenpot
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

 subroutine pawxenergy(eexex,pawprtvol,pawrhoij,pawtab)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: pawprtvol
 real(dp),intent(inout) :: eexex
 type(pawrhoij_type),intent(in) :: pawrhoij
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: irhoij,irhoij1,ispden,jrhoij,jrhoij1,klmn,klmn1,lexexch,ll,m11,m21,m31,m41,n1
 integer :: n2,n3,n4,nk,nn1,nn2
 real(dp) :: eexextemp
 character(len=500) :: message
!arrays
 integer :: indn(3,3)
 real(dp) :: factnk(6)

! *****************************************************

 DBG_ENTER("COLL")

 if (pawrhoij%qphase==2) then
   message='pawxenergy: local exact-exchange not compatible with qphase=2!'
   MSG_ERROR(message)
 end if

 lexexch=pawtab%lexexch
 if (pawtab%nproju==1) nk=1
 if (pawtab%nproju==2) nk=6
 factnk(1)=one;factnk(2)=one;factnk(3)=one
 factnk(4)=two;factnk(5)=two;factnk(6)=two
 indn(1,1)=1;indn(1,2)=4;indn(1,3)=5
 indn(2,1)=4;indn(2,2)=2;indn(2,3)=6
 indn(3,1)=5;indn(3,2)=6;indn(3,3)=3

!======================================================
!Compute local exact exchange Energy
!-----------------------------------------------------
 eexextemp=zero

 do ispden=1,pawrhoij%nspden
   jrhoij=1
   do irhoij=1,pawrhoij%nrhoijsel
     klmn=pawrhoij%rhoijselect(irhoij)
     if(pawtab%indklmn(3,klmn)==0.and.pawtab%indklmn(4,klmn)==2*lexexch) then
       m11=pawtab%klmntomn(1,klmn);m21=pawtab%klmntomn(2,klmn)
       n1=pawtab%klmntomn(3,klmn);n2=pawtab%klmntomn(4,klmn)
       nn1=(n1*n2)/2+1
       jrhoij1=1
       do irhoij1=1,pawrhoij%nrhoijsel
         klmn1=pawrhoij%rhoijselect(irhoij1)
         if(pawtab%indklmn(3,klmn1)==0.and.pawtab%indklmn(4,klmn1)==2*lexexch) then
           m31=pawtab%klmntomn(1,klmn1);m41=pawtab%klmntomn(2,klmn1)
           n3=pawtab%klmntomn(3,klmn1);n4=pawtab%klmntomn(4,klmn1)
           nn2=(n3*n4)/2+1
           do ll=1,lexexch+1
             eexextemp=eexextemp-pawtab%vex(m11,m31,m41,m21,ll)*pawtab%dltij(klmn)*pawtab%fk(indn(nn1,nn2),ll)&
&             *pawtab%dltij(klmn1)*pawrhoij%rhoijp(jrhoij,ispden)*pawrhoij%rhoijp(jrhoij1,ispden)
           end do
         end if
         jrhoij1=jrhoij1+pawrhoij%cplex_rhoij
       end do !irhoij1
     end if
     jrhoij=jrhoij+pawrhoij%cplex_rhoij
   end do !irhoij
 end do ! ispden
 eexextemp=eexextemp/two
 eexex=eexex+eexextemp*pawtab%exchmix

 if (abs(pawprtvol)>=2) then
   write(message, '(a)' )"   Contributions to the direct expression of energy:"
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a,e20.10,a)') "     HF exchange energy  =",eexextemp,ch10
   call wrtout(std_out,message,'COLL')
 end if

 DBG_EXIT("COLL")

 end subroutine pawxenergy
!!***

!----------------------------------------------------------------------

!!****f* m_paw_correlations/setnoccmmp
!! NAME
!! setnoccmmp
!!
!! FUNCTION
!! PAW+U only:
!!   Compute density matrix nocc_{m,m_prime}
!!   or
!!   Impose value of density matrix using dmatpawu input array, then symetrize it.
!!
!! noccmmp^{\sigma}_{m,m'}=\sum_{ni,nj}[\rho^{\sigma}_{ni,nj}*phiphjint_{ni,nj}]
!!
!! INPUTS
!!  compute_dmat= flag: if 1, nocc_{m,mp} is computed
!!  dimdmat=first dimension of dmatpawu array
!!  dmatpawu(dimdmat,dimdmat,nsppol*nspinor,natpawu)=input density matrix to be copied into noccmpp
!!  dmatudiag= flag controlling the use of diagonalization:
!!             0: no diagonalization of nocc_{m,mp}
!!             1: diagonalized nocc_{m,mp} matrix is printed
!!             2: dmatpawu matrix is expressed in the basis where nocc_(m,mp} is diagonal
!!  impose_dmat= flag: if 1, nocc_{m,mp} is replaced by dmatpawu
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell
!!  natpawu=number of atoms on which PAW+U is applied
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=number of independant spin components
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of atom types
!!  paw_ij(my_natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  spinat(3,matom)=initial spin of each atom, in unit of hbar/2
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  typat(natom)=type for each atom
!!  useexexch=1 if local-exact-exchange is activated
!!  usepawu= /=0 if PAW+U is activated
!!
!! OUTPUT
!!   paw_ij(natom)%noccmmp(cplex_dij,2*lpawu+1,2*lpawu+1,nsppol or ndij)=density matrix
!!
!! NOTES
!! For non-collinear magnetism,
!! - nocc_{m,mp} is computed as:noccmmp(:,:,:,1)=   n{m,mp}
!!                              noccmmp(:,:,:,2)=   m_x{m,mp}
!!                              noccmmp(:,:,:,3)=   m_y{m,mp}
!!                              noccmmp(:,:,:,4)=   m_z{m,mp}
!! - but nocc_{m,mp} is stored as: noccmmp(:,:,:,1)=   n^{up,up}_{m,mp}
!!                                 noccmmp(:,:,:,2)=   n^{dn,dn}_{m,mp}
!!                                 noccmmp(:,:,:,3)=   n^{up,dn}_{m,mp}
!!                                 noccmmp(:,:,:,4)=   n^{dn,up}_{m,mp}
!!   We choose to have noccmmp complex when ndij=4 (ie nspinor=2)
!!    If ndij=4 and pawspnorb=0, one could keep noccmmp real
!!    with the n11, n22, Re(n12), Im(n21) representation, but it would
!!    less clear to change the representation when pawspnorb is activated.
!!   If ndij=4, nocc_{m,mp} is transformed to the Ylm basis
!!    and then to the J, M_J basis (if cplex_dij==2)
!!
!!  Note that n_{m,mp}=<mp|hat(n)|m> because rhoij=<p_j|...|p_i>
!!
!! PARENTS
!!      afterscfloop,pawdenpot,pawprt,scfcv,vtorho
!!
!! CHILDREN
!!      dgemm,dsyev,free_my_atmtab,get_my_atmtab,mat_mlms2jmj,mat_slm2ylm
!!      wrtout,zgemm,zheev
!!
!! SOURCE

subroutine setnoccmmp(compute_dmat,dimdmat,dmatpawu,dmatudiag,impose_dmat,indsym,my_natom,natom,&
&                     natpawu,nspinor,nsppol,nsym,ntypat,paw_ij,pawang,pawprtvol,pawrhoij,pawtab,&
&                     spinat,symafm,typat,useexexch,usepawu, &
&                     mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: compute_dmat,dimdmat,dmatudiag,impose_dmat,my_natom,natom,natpawu
 integer,intent(in) :: nspinor,nsppol,nsym,ntypat,useexexch,usepawu
 integer,optional,intent(in) :: comm_atom
 type(pawang_type),intent(in) :: pawang
 integer,intent(in) :: pawprtvol
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symafm(nsym),typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: dmatpawu(dimdmat,dimdmat,nspinor*nsppol,natpawu*impose_dmat)
 !real(dp),intent(in) :: dmatpawu(:,:,:,:)
 real(dp),intent(in) :: spinat(3,natom)
 type(paw_ij_type),intent(inout) :: paw_ij(my_natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer,parameter :: limp=0 ! could become an input variable
 integer :: at_indx,cplex_dij,cplex_rhoij,dmatudiag_loc,iafm,iatom,iatom_tot,iatpawu,icount
 integer :: ilm,im1,im2,in1,in2,info,iplex,irot,ispden, irhoij,itypat,jlm,jrhoij
 integer :: jspden,klmn,kspden,lcur,ldim,lmax,lmin,lpawu,lwork,my_comm_atom,ndij,nmat,nspden,nsploop
 logical,parameter :: afm_noncoll=.true.  ! TRUE if antiferro symmetries are used with non-collinear magnetism
 logical :: antiferro,my_atmtab_allocated,noccsym_error,paral_atom,use_afm
 real(dp),parameter :: invsqrt2=one/sqrt2
 real(dp) :: factafm,mnorm,mx,my,mz,ntot,nup,ndn,snorm,sx,sy,szp,szm
 character(len=4) :: wrt_mode
 character(len=500) :: message
!arrays
 integer :: nsym_used(2)
 integer,pointer :: my_atmtab(:)
 real(dp) :: ro(2),sumocc(2)
 real(dp),allocatable :: eig(:),hdp(:,:,:),hdp2(:,:),noccmmptemp(:,:,:,:),noccmmp_tmp(:,:,:,:)
 real(dp),allocatable :: rwork(:),noccmmp2(:,:,:,:),nocctot2(:)
 complex(dpc),allocatable :: noccmmp_ylm(:,:,:),noccmmp_jmj(:,:),noccmmp_slm(:,:,:)
 complex(dpc),allocatable :: zhdp(:,:),zhdp2(:,:),znoccmmp_tmp(:,:),zwork(:)
 character(len=9),parameter :: dspin(6)=  (/"up       ","down     ","up-up    ","down-down","Re[up-dn]","Im[up-dn]"/)
 character(len=9),parameter :: dspinc(6)= (/"up       ","down     ","up-up    ","down-down","up-dn    ","dn-up    "/)
 character(len=9),parameter :: dspinc2(6)=(/"up       ","down     ","dn-dn    ","up-up    ","dn-up    ","up-dn    "/)
 character(len=9),parameter :: dspinm(6)= (/"dn       ","up i     ","n        ","mx       ","my       ","mz       "/)
 type(coeff4_type),allocatable :: tmp_noccmmp(:)

!*********************************************************************

 DBG_ENTER("COLL")

!Tests
 if (my_natom>0) then
   if (nsppol/=paw_ij(1)%nsppol) then
     message='inconsistent values for nsppol!'
     MSG_BUG(message)
   end if
   if (compute_dmat>0) then
     if (pawrhoij(1)%nspden/=paw_ij(1)%nspden.and.&
&        pawrhoij(1)%nspden/=4.and.paw_ij(1)%nspden/=1) then
       message=' inconsistent values for nspden!'
       MSG_BUG(message)
     end if
   end if
   if (pawrhoij(1)%qphase==2) then
     message='setnoccmmp not compatible with qphase=2!'
     MSG_BUG(message)
   end if
 end if
 if (usepawu/=0.and.useexexch/=0) then
   message='usepawu/=0 and useexexch>0 not allowed!'
   MSG_BUG(message)
 end if
 if (impose_dmat/=0.and.dimdmat==0) then
   message='dmatpawu must be allocated when impose_dmat/=0!'
   MSG_BUG(message)
 end if
 if (usepawu>0.and.compute_dmat/=0.and.impose_dmat/=0.and.pawang%nsym==0) then
   message='pawang%zarot must be allocated!'
   MSG_BUG(message)
 end if

!Some inits
 if (usepawu==0.and.useexexch==0) return
 nspden=1;ndij=1;cplex_dij=1
 if (my_natom>0) then
   nspden=paw_ij(1)%nspden
   ndij=paw_ij(1)%ndij
   cplex_dij=paw_ij(1)%cplex_dij
 end if
 antiferro=(nspden==2.and.nsppol==1)
 use_afm=((antiferro).or.((nspden==4).and.afm_noncoll))
 dmatudiag_loc=dmatudiag
 if (dmatudiag==2.and.(dimdmat==0.or.impose_dmat==0)) dmatudiag_loc=1

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom) !vz_d
 wrt_mode='COLL';if (paral_atom) wrt_mode='PERS'

!If needed, store dmatpu in suitable format in tmp_noccmmp
 if (usepawu/=0.and.impose_dmat/=0) then
   iatpawu=0
   ABI_DATATYPE_ALLOCATE(tmp_noccmmp,(natom))
   do iatom_tot=1,natom
     itypat=typat(iatom_tot)
     lpawu=pawtab(itypat)%lpawu
     if (lpawu/=-1) then
       iatpawu=iatpawu+1
       if (ndij/=4) then
         ABI_ALLOCATE(tmp_noccmmp(iatom_tot)%value,(cplex_dij,2*lpawu+1,2*lpawu+1,nsppol))
         tmp_noccmmp(iatom_tot)%value(1,1:2*lpawu+1,1:2*lpawu+1,1:nsppol)=&
&         dmatpawu(1:2*lpawu+1,1:2*lpawu+1,1:nsppol,iatpawu)
       else
         ABI_ALLOCATE(tmp_noccmmp(iatom_tot)%value,(cplex_dij,2*lpawu+1,2*lpawu+1,ndij))
         tmp_noccmmp(iatom_tot)%value=zero
         if(limp==0) then ! default reading
           snorm=sqrt(spinat(1,iatom_tot)**2+spinat(1,iatom_tot)**2+spinat(3,iatom_tot)**2)
           if (snorm>tol12) then
             sx=half*spinat(1,iatom_tot)/snorm
             sy=half*spinat(2,iatom_tot)/snorm
             szp=half*(one+spinat(3,iatom_tot)/snorm)
             szm=half*(one-spinat(3,iatom_tot)/snorm)
           else
             sx=zero;sy=zero
             szp=one;szm=zero
           end if
           do im2=1,2*lpawu+1
             do im1=1,2*lpawu+1
               nup=dmatpawu(im1,im2,1,iatpawu);ndn=dmatpawu(im1,im2,2,iatpawu)
               tmp_noccmmp(iatom_tot)%value(1,im1,im2,1)=nup*szp+ndn*szm
               tmp_noccmmp(iatom_tot)%value(1,im1,im2,2)=nup*szm+ndn*szp
               tmp_noccmmp(iatom_tot)%value(1,im1,im2,3)=(nup-ndn)*sx
               tmp_noccmmp(iatom_tot)%value(1,im1,im2,4)=(ndn-nup)*sy
             end do
           end do
         else if(limp>=1) then
           ABI_ALLOCATE(noccmmp_ylm,(2*lpawu+1,2*lpawu+1,ndij))
           noccmmp_ylm=czero
           ABI_ALLOCATE(noccmmp_slm,(2*lpawu+1,2*lpawu+1,ndij))
           noccmmp_slm=czero
           ABI_ALLOCATE(noccmmp_jmj,(2*(2*lpawu+1),2*(2*lpawu+1)))
           noccmmp_jmj=czero
           if(limp==1) then ! read input matrix in J,M_J basis (l-1/2, then l+1/2)
             noccmmp_jmj=czero
             do im1=1,2*lpawu+1
               noccmmp_jmj(im1,im1)=cmplx(dmatpawu(im1,im1,1,iatpawu),zero,kind=dp)
               noccmmp_jmj(im1+lpawu,im1+lpawu)=cmplx(dmatpawu(im1+lpawu,im1+lpawu,2,iatpawu),zero,kind=dp)
             end do
             write(message,'(2a,i3,a)') ch10,'== Atom ',iatom_tot,&
&             ' == Imposed occupation matrix (in the J M_J basis: L-1/2 and L+1/2 states)'
             call wrtout(std_out,message,wrt_mode)
             call mat_mlms2jmj(lpawu,noccmmp_ylm,noccmmp_jmj,ndij,&
&             2,2,pawprtvol,std_out,wrt_mode) !  optspin=1: up spin are first
           end if
           if(limp==2) then ! read input matrix in Ylm basis
             noccmmp_ylm=czero
             do im1=1,2*lpawu+1
               noccmmp_ylm(im1,im1,1)=cmplx(dmatpawu(im1,im1,1,iatpawu),zero,kind=dp)
               noccmmp_ylm(im1,im1,2)=cmplx(dmatpawu(im1,im1,2,iatpawu),zero,kind=dp)
             end do
             write(message,'(2a,i3,a)') ch10,'== Atom ',iatom_tot,&
&             ' == Imposed occupation matrix (in the Ylm basis), for dn and up spin'
             call wrtout(std_out,message,wrt_mode)
           end if
           call mat_slm2ylm(lpawu,noccmmp_ylm,noccmmp_slm,ndij,&
&           2,2,pawprtvol,std_out,wrt_mode) ! optspin=1 because up spin are first
!          interchange upup and dndn
           if(limp>=1) then
             tmp_noccmmp(iatom_tot)%value(1,:,:,1)=real(noccmmp_slm(:,:,2))
             tmp_noccmmp(iatom_tot)%value(2,:,:,1)=aimag(noccmmp_slm(:,:,2))
             tmp_noccmmp(iatom_tot)%value(1,:,:,2)=real(noccmmp_slm(:,:,1))
             tmp_noccmmp(iatom_tot)%value(2,:,:,2)=aimag(noccmmp_slm(:,:,1))
             tmp_noccmmp(iatom_tot)%value(1,:,:,3)=real(noccmmp_slm(:,:,4))
             tmp_noccmmp(iatom_tot)%value(2,:,:,3)=aimag(noccmmp_slm(:,:,4))
             tmp_noccmmp(iatom_tot)%value(1,:,:,4)=real(noccmmp_slm(:,:,3))
             tmp_noccmmp(iatom_tot)%value(2,:,:,4)=aimag(noccmmp_slm(:,:,3))
           end if
           if(abs(pawprtvol)>2) then
             write(message, '(2a)' ) ch10,&
&             " Check Imposed density matrix in different basis"
             call wrtout(std_out,message,wrt_mode)
             call mat_slm2ylm(lpawu,noccmmp_slm,noccmmp_ylm,ndij,&
&             1,2,pawprtvol,std_out,wrt_mode) ! optspin=1 because up spin are first
             call mat_mlms2jmj(lpawu,noccmmp_ylm,noccmmp_jmj,ndij,1,2,&
&             pawprtvol,std_out,wrt_mode) !  optspin=1: up spin are first
           end if
           ABI_DEALLOCATE(noccmmp_ylm)
           ABI_DEALLOCATE(noccmmp_jmj)
           ABI_DEALLOCATE(noccmmp_slm)
         end if
       end if
     end if
   end do
 end if  ! impose_dmat/=0

!Print message
 if (usepawu/=0.and.impose_dmat/=0) then
   if (dmatudiag_loc/=2) then
     write(message,'(6a)') ch10,'Occupation matrix for correlated orbitals is kept constant',ch10,&
&     'and equal to dmatpawu from input file !',ch10,&
&     '----------------------------------------------------------'
   else
     write(message,'(6a)') ch10,'Occupation matrix for correlated orbitals is imposed',ch10,&
&     'and equal to dmatpawu in the diagonal basis !',ch10,&
&     '----------------------------------------------------------'
   end if
   call wrtout(std_out,message,'COLL')
 end if

 if (usepawu/=0.and.dmatudiag_loc/=0.and.my_natom>0) then
   write(message,'(4a)') ch10,'Diagonalized occupation matrix "noccmmp" is printed !',ch10,&
&   '-------------------------------------------------------------'
   call wrtout(std_out,message,wrt_mode)
 end if

!Loops over atoms
 do iatom=1,my_natom
   iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)
   itypat=pawrhoij(iatom)%itypat
   cplex_rhoij=pawrhoij(iatom)%cplex_rhoij

   if (useexexch/=0) then
     lcur=pawtab(itypat)%lexexch
   else if (usepawu/=0) then
     lcur=pawtab(itypat)%lpawu
   end if
   if (lcur/=-1) then

!    ########################################################################################
!    # Compute nocc_mmp
!    ########################################################################################
     if ((usepawu/=0.and.compute_dmat/=0).or.useexexch/=0) then


       paw_ij(iatom)%noccmmp(:,:,:,:)=zero

!      Loop over spin components
       ABI_ALLOCATE(noccmmptemp,(cplex_dij,2*lcur+1,2*lcur+1,ndij))
       noccmmptemp(:,:,:,:)=zero
       if(ndij==4)  then
         ABI_ALLOCATE(noccmmp2,(cplex_dij,2*lcur+1,2*lcur+1,ndij))
       end if
       if(ndij==4)  then
         if(allocated(nocctot2)) then
           ABI_DEALLOCATE(nocctot2)
         end if
         ABI_ALLOCATE(nocctot2,(ndij))
       end if
       do ispden=1,ndij
         jrhoij=1
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=pawrhoij(iatom)%rhoijselect(irhoij)
           im1=pawtab(itypat)%klmntomn(1,klmn)
           im2=pawtab(itypat)%klmntomn(2,klmn)
           in1=pawtab(itypat)%klmntomn(3,klmn)
           in2=pawtab(itypat)%klmntomn(4,klmn)
           lmin=pawtab(itypat)%indklmn(3,klmn)
           lmax=pawtab(itypat)%indklmn(4,klmn)

           ro(1:2)=zero
           ro(1:cplex_rhoij)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+cplex_rhoij-1,ispden)
           if (ndij==1) ro(1:2)=half*ro(1:2)
!          Non-collinear magnetism: keep n, m storage because
!            it is easier for the computation of noccmmp from rhoij)

           if(lmin==0.and.lmax==2*lcur) then
             icount=in1+(in2*(in2-1))/2
             if(pawtab(itypat)%ij_proj<icount)  then
               message='PAW+U: Problem in the loop calculating noccmmp!'
               MSG_BUG(message)
             end if
             if(in1/=in2) then
               if(im2<=im1) then
                 noccmmptemp(1:cplex_dij,im1,im2,ispden)=noccmmptemp(1:cplex_dij,im1,im2,ispden) &
&                           +ro(1:cplex_dij)*pawtab(itypat)%phiphjint(icount)
               end if
             end if
             if(im2>=im1) then
               paw_ij(iatom)%noccmmp(1:cplex_dij,im1,im2,ispden)=paw_ij(iatom)%noccmmp(1:cplex_dij,im1,im2,ispden) &
&                           +ro(1:cplex_dij)*pawtab(itypat)%phiphjint(icount)
             end if
           end if
           jrhoij=jrhoij+cplex_rhoij
         end do ! irhoij
         do im2=1,2*lcur+1
           do im1=1,im2
             paw_ij(iatom)%noccmmp(1,im1,im2,ispden)=paw_ij(iatom)%noccmmp(1,im1,im2,ispden) &
&             +noccmmptemp(1,im2,im1,ispden)
             if(cplex_dij==2) paw_ij(iatom)%noccmmp(2,im1,im2,ispden)=paw_ij(iatom)%noccmmp(2,im1,im2,ispden) &
&             -noccmmptemp(2,im2,im1,ispden)
           end do
         end do
         do im1=1,2*lcur+1
           do im2=1,im1
             paw_ij(iatom)%noccmmp(1,im1,im2,ispden)=paw_ij(iatom)%noccmmp(1,im2,im1,ispden)
             if(cplex_dij==2) paw_ij(iatom)%noccmmp(2,im1,im2,ispden)=-paw_ij(iatom)%noccmmp(2,im2,im1,ispden)
           end do
         end do
       end do ! ispden
       ABI_DEALLOCATE(noccmmptemp)
!      Compute noccmmp2, occupation matrix in the spin basis (upup, dndn, updn, dnup)
       if(ndij==4) then
         noccmmp2(:,:,:,:)=zero
         do im1=1,2*lcur+1
           do im2=1,2*lcur+1
             noccmmp2(1,im1,im2,1)=half*(paw_ij(iatom)%noccmmp(1,im1,im2,1)+paw_ij(iatom)%noccmmp(1,im1,im2,4))
             noccmmp2(2,im1,im2,1)=half*(paw_ij(iatom)%noccmmp(2,im1,im2,1)+paw_ij(iatom)%noccmmp(2,im1,im2,4))
             noccmmp2(1,im1,im2,2)=half*(paw_ij(iatom)%noccmmp(1,im1,im2,1)-paw_ij(iatom)%noccmmp(1,im1,im2,4))
             noccmmp2(2,im1,im2,2)=half*(paw_ij(iatom)%noccmmp(2,im1,im2,1)-paw_ij(iatom)%noccmmp(2,im1,im2,4))
             noccmmp2(1,im1,im2,3)=half*(paw_ij(iatom)%noccmmp(1,im1,im2,2)+paw_ij(iatom)%noccmmp(2,im1,im2,3))
             noccmmp2(2,im1,im2,3)=half*(paw_ij(iatom)%noccmmp(2,im1,im2,2)-paw_ij(iatom)%noccmmp(1,im1,im2,3))
             noccmmp2(1,im1,im2,4)=half*(paw_ij(iatom)%noccmmp(1,im1,im2,2)-paw_ij(iatom)%noccmmp(2,im1,im2,3))
             noccmmp2(2,im1,im2,4)=half*(paw_ij(iatom)%noccmmp(2,im1,im2,2)+paw_ij(iatom)%noccmmp(1,im1,im2,3))
           end do
         end do
         if(abs(pawprtvol)>=1) then
           write(message,'(2a)') ch10,"== Calculated occupation matrix for correlated orbitals in the n, m basis :"
           call wrtout(std_out,message,wrt_mode)
           do ispden=1,ndij
             write(message,'(3a)') ch10,"Calculated occupation matrix for component ",trim(dspinm(ispden+2*(ndij/4)))
             call wrtout(std_out,message,wrt_mode)
             do im1=1,lcur*2+1  ! ( order of indices in noccmmp is exchanged in order to have the same convention as rhoij: transposition is done after )
               if(cplex_dij==1)&
&               write(message,'(12(1x,9(1x,f10.5)))')&
&               (paw_ij(iatom)%noccmmp(1,im2,im1,ispden),im2=1,lcur*2+1)
               if(cplex_dij==2)&
!              &               write(message,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))')&
&               write(message,'(12(1x,9(1x,"(",f10.5,",",f10.5,")")))')&
&               (paw_ij(iatom)%noccmmp(:,im2,im1,ispden),im2=1,lcur*2+1)
               call wrtout(std_out,message,wrt_mode)
             end do
           end do
         end if ! pawprtvol >=1
       end if

!      Compute total number of electrons per spin
       paw_ij(iatom)%nocctot(:)=zero ! contains nmmp in the n m representation
       if(ndij==4) nocctot2(:)=zero ! contains nmmp in the upup dndn updn dnup  representation
       do ispden=1,ndij
         do im1=1,2*lcur+1
           if(ndij==4) then
             paw_ij(iatom)%nocctot(ispden)=paw_ij(iatom)%nocctot(ispden)+paw_ij(iatom)%noccmmp(1,im1,im1,ispden)
             nocctot2(ispden)=nocctot2(ispden)+noccmmp2(1,im1,im1,ispden)
           else
             paw_ij(iatom)%nocctot(ispden)=paw_ij(iatom)%nocctot(ispden)+paw_ij(iatom)%noccmmp(1,im1,im1,ispden)
           end if
         end do
       end do
!      noccmmp will now be in the up up , dn dn... representation and now n_mmp=<m|n|mp> instead of <mp|n|m> !
       if(ndij==4) then
         do ispden=1,ndij
           do iplex=1,cplex_dij
             do im1=1,2*lcur+1
               do im2=1,2*lcur+1
                 paw_ij(iatom)%noccmmp(iplex,im1,im2,ispden)=noccmmp2(iplex,im2,im1,ispden) ! now, noccmmp is in the upup dndn updn dnup representation
               end do
             end do
           end do
         end do
         ABI_DEALLOCATE(noccmmp2)
       end if
!      Printing of new nocc_mmp
       if ((usepawu/=0.and.abs(usepawu)<10).or.(usepawu>=10.and.pawprtvol>=3)) then
         write(message, '(2a)' )  ch10, &
&         '========== LDA+U DATA =================================================== '
       end if
       if (useexexch/=0) then
         write(message, '(2a)' ) ch10, &
&         '======= Local ex-exchange (PBE0) DATA =================================== '
       end if
       if (((usepawu/=0.and.abs(usepawu)<10).or.(usepawu>=10.and.pawprtvol>=3)).or.useexexch/=0) then
         call wrtout(std_out,message,wrt_mode)
       end if
       if (usepawu>=10.and.pawprtvol>=3) then
         write(message, '(6a)' )  ch10,'    ( A DFT+DMFT calculation is carried out                              ',&
         ch10,'      Thus, the following LDA+U occupation matrices are not physical     ',&
         ch10,'      and just informative )'
         call wrtout(std_out,message,wrt_mode)
       end if
       if(abs(usepawu)<10.or.pawprtvol>=3) then ! Always write except if DMFT and pawprtvol low
         write(message,'(2a,i5,a,i4,a)') ch10,"====== For Atom", iatom_tot,&
&         ", occupations for correlated orbitals. l =",lcur,ch10
         call wrtout(std_out,message,wrt_mode)
         if(ndij==2) then
           do ispden=1,2
             write(message,'(a,i4,3a,f10.5)') "Atom", iatom_tot,". Occupations for spin ",&
&             trim(dspin(ispden))," =",paw_ij(iatom)%nocctot(ispden)
             call wrtout(std_out,message,wrt_mode)
           end do
           write(message,'(a,i4,a,2x,e16.8)') "=> On atom",iatom_tot,", local Mag. is  ",&
&           paw_ij(iatom)%nocctot(2)-paw_ij(iatom)%nocctot(1)
           call wrtout(std_out,message,wrt_mode)
         end if
         if(ndij==4) then
           ntot=paw_ij(iatom)%nocctot(1)
           mx=paw_ij(iatom)%nocctot(2)
           my=paw_ij(iatom)%nocctot(3)
           mz=paw_ij(iatom)%nocctot(4)
           mnorm=sqrt(mx*mx+my*my+mz*mz)
           nup=nocctot2(1)
           ndn=nocctot2(2)
           write(message,'(a,i4,a,2x,e16.8)') "=> On atom",iatom_tot,", local Mag. x is ",mx
           call wrtout(std_out,message,wrt_mode)
           write(message,'(14x,a,2x,e16.8)') "  local Mag. y is ",my
           call wrtout(std_out,message,wrt_mode)
           write(message,'(14x,a,2x,e16.8)') "  local Mag. z is ",mz
           call wrtout(std_out,message,wrt_mode)
           write(message,'(14x,a,2x,e16.8)') "  norm of Mag. is ",mnorm
           call wrtout(std_out,message,wrt_mode)
           write(message,'(14x,a,2x,f10.5)') "  occ. of majority spin is ",half*(ntot+mnorm)  ! to be checked versus direct calc from noccmmp
           call wrtout(std_out,message,wrt_mode)
           if(abs(pawprtvol)>=1) write(message,'(14x,a,2x,f10.5)') "  occ. for spin up (along z) ",nup
           if(abs(pawprtvol)>=1) then
             call wrtout(std_out,message,wrt_mode)
           end if
           write(message,'(14x,a,2x,f10.5)') "  occ. of minority spin is ",half*(ntot-mnorm)
           call wrtout(std_out,message,wrt_mode)
           if(abs(pawprtvol)>=1) write(message,'(14x,a,2x,f10.5)') "  occ. for spin dn (along z) ",ndn
           if(abs(pawprtvol)>=1) then
             call wrtout(std_out,message,wrt_mode)
           end if
           if(ndij==4)  then
             ABI_DEALLOCATE(nocctot2)
           end if
         end if
         write(message,'(2a)') ch10,"== Calculated occupation matrix for correlated orbitals:"
         call wrtout(std_out,message,wrt_mode)
         do ispden=1,ndij
           write(message,'(3a)') ch10,"Calculated occupation matrix for component ",trim(dspinc(ispden+2*(ndij/4)))
           call wrtout(std_out,message,wrt_mode)
           do im1=1,lcur*2+1
             if(cplex_dij==1)&
&             write(message,'(12(1x,9(1x,f10.5)))')&
&             (paw_ij(iatom)%noccmmp(1,im1,im2,ispden),im2=1,lcur*2+1)
             if(cplex_dij==2)&
&             write(message,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))')&
&             (paw_ij(iatom)%noccmmp(:,im1,im2,ispden),im2=1,lcur*2+1)
             call wrtout(std_out,message,wrt_mode)
           end do
         end do
       end if

!      Transformation matrices: real->complex spherical harmonics (for test)
       if(ndij==4.and.abs(pawprtvol)>=0) then
         ABI_ALLOCATE(noccmmp_ylm,(2*lcur+1,2*lcur+1,ndij))
         noccmmp_ylm=czero
         ABI_ALLOCATE(noccmmp_slm,(2*lcur+1,2*lcur+1,ndij))
         noccmmp_slm=czero
         ABI_ALLOCATE(noccmmp_jmj,(2*(2*lcur+1),2*(2*lcur+1)))
         noccmmp_jmj=czero
!        go from real notation for complex noccmmp to complex notation in noccmmp_slm
         noccmmp_slm(:,:,:)=cmplx(paw_ij(iatom)%noccmmp(1,:,:,:)&
&         ,paw_ij(iatom)%noccmmp(2,:,:,:),kind=dp)
         call mat_slm2ylm(lcur,noccmmp_slm,noccmmp_ylm,ndij,1,1,pawprtvol,std_out,wrt_mode) ! optspin=1: up spin are first

         do ispden=1,ndij
           write(message,'(3a)') ch10,"Calculated Ylm occupation matrix for component ",trim(dspinc(ispden+2*(ndij/4)))
           call wrtout(std_out,message,wrt_mode)
           do im1=1,lcur*2+1
             write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))') (noccmmp_ylm(im1,im2,ispden),im2=1,lcur*2+1)
             call wrtout(std_out,message,wrt_mode)
           end do
         end do
         call mat_mlms2jmj(lcur,noccmmp_ylm,noccmmp_jmj,ndij,1,1,pawprtvol,std_out,wrt_mode) !  optspin=1: up spin are first
         ABI_DEALLOCATE(noccmmp_ylm)
         ABI_DEALLOCATE(noccmmp_jmj)
         ABI_DEALLOCATE(noccmmp_slm)
       end if !ndij==4

     end if ! impose_dmat==0

!    ########################################################################################
!    # Diagonalize nocc_mmp
!    ########################################################################################
     if(usepawu/=0.and.dmatudiag_loc>0) then

       lpawu=lcur;ldim=2*lpawu+1
       ABI_ALLOCATE(noccmmp_tmp,(1,ldim,ldim,ndij))
       if (ndij==4)  then
         ABI_ALLOCATE(znoccmmp_tmp,(2*ldim,2*ldim))
       end if

!      Select noccmmp for this atom
       do ispden=1,ndij
         noccmmp_tmp(1,:,:,ispden)=paw_ij(iatom)%noccmmp(1,:,:,ispden)
       end do
       if (ndij==4) then
         do im2=1,ldim
           do im1=1,ldim
             znoccmmp_tmp(im1     ,     im2)=cmplx(paw_ij(iatom)%noccmmp(1,im1,im2,1)&
&             ,paw_ij(iatom)%noccmmp(2,im1,im2,1),kind=dp)
             znoccmmp_tmp(ldim+im1,ldim+im2)=cmplx(paw_ij(iatom)%noccmmp(1,im1,im2,2)&
&             ,paw_ij(iatom)%noccmmp(2,im1,im2,2),kind=dp)
             znoccmmp_tmp(     im1,ldim+im2)=cmplx(paw_ij(iatom)%noccmmp(1,im1,im2,3)&
&             ,paw_ij(iatom)%noccmmp(2,im1,im2,3),kind=dp)
             znoccmmp_tmp(ldim+im1,     im2)=cmplx(paw_ij(iatom)%noccmmp(1,im1,im2,4)&
&             ,paw_ij(iatom)%noccmmp(2,im1,im2,4),kind=dp)
           end do
         end do
       end if

!      Diagonalize nocc_mmp
       if (ndij/=4) then
         ABI_ALLOCATE(hdp,(ldim,ldim,ndij))
         hdp=zero
         lwork=3*ldim-1
         ABI_ALLOCATE(rwork,(lwork))
         ABI_ALLOCATE(eig,(ldim))
         do ispden=1,ndij
           call dsyev('v','u',ldim,noccmmp_tmp(1,:,:,ispden),ldim,eig,rwork,lwork,info)
           if(info/=0) then
             message=' Error in diagonalization of noccmmp (DSYEV)!'
             MSG_ERROR(message)
           end if
           do ilm=1,ldim
             hdp(ilm,ilm,ispden)=eig(ilm)
           end do
         end do ! ispden
         ABI_DEALLOCATE(rwork)
         ABI_DEALLOCATE(eig)
       else
         ABI_ALLOCATE(hdp,(2*ldim,2*ldim,1))
         hdp=zero
         lwork=4*ldim-1
         ABI_ALLOCATE(rwork,(6*ldim-2))
         ABI_ALLOCATE(zwork,(lwork))
         ABI_ALLOCATE(eig,(2*ldim))
         call zheev('v','u',2*ldim,znoccmmp_tmp,2*ldim,eig,zwork,lwork,rwork,info)
         if(info/=0) then
           message=' Error in diagonalization of znoccmmp_tmp (zheev) !'
           MSG_ERROR(message)
         end if
         do ilm=1,2*ldim
           hdp(ilm,ilm,1)=eig(ilm)
         end do
         ABI_DEALLOCATE(rwork)
         ABI_DEALLOCATE(zwork)
         ABI_DEALLOCATE(eig)
       end if

!      Print diagonalized matrix and eigenvectors
       do ispden=1,size(hdp,3)
         write(message,'(2a,i3,a)') ch10,'== Atom ',iatom_tot,' == Diagonalized Occupation matrix'
         if (ndij==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
         if (ndij==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
         if (ndij==4) write(message,fmt='(2a,i3,a)')trim(message)," =="
         call wrtout(std_out,message,wrt_mode)
         do ilm=1,size(hdp,1)
           write(message,'(12(1x,9(1x,f10.5)))') (hdp(ilm,jlm,ispden),jlm=1,size(hdp,2))
           call wrtout(std_out,message,wrt_mode)
         end do
       end do ! ispden
       if(abs(pawprtvol)>=1) then
         if (ndij/=4) then
           do ispden=1,ndij
             write(message,'(2a,i3,a)') ch10,'== Atom ',iatom_tot,' == Eigenvectors'
             if (ndij==1) write(message,fmt='(2a)')     trim(message),' for spin up =='
             if (ndij==2) write(message,fmt='(2a,i3,a)')trim(message),' for spin ',ispden,' =='
             call wrtout(std_out,message,wrt_mode)
             do ilm=1,ldim
               write(message,'(12(1x,9(1x,f10.5)))') (noccmmp_tmp(1,ilm,jlm,ispden),jlm=1,ldim)
               call wrtout(std_out,message,wrt_mode)
             end do
           end do
         else
           write(message,'(2a,i3,a)') ch10,'== Atom ',iatom_tot,' == Eigenvectors (spinors) in the real harmonics basis =='
           call wrtout(std_out,message,wrt_mode)
           do ilm=1,2*ldim
             write(message,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))') (znoccmmp_tmp(ilm,jlm),jlm=1,2*ldim)
             call wrtout(std_out,message,wrt_mode)
           end do
         end if
       end if

!      Back rotation of diagonalized matrix and printing
       if(abs(pawprtvol)>=1) then
         if (ndij/=4) then
           ABI_ALLOCATE(hdp2,(ldim,ldim))
           do ispden=1,ndij
             call dgemm('n','t',ldim,ldim,ldim,one,hdp(:,:,ispden),ldim,noccmmp_tmp(1,:,:,ispden),ldim,zero,hdp2,ldim)
             call dgemm('n','n',ldim,ldim,ldim,one,noccmmp_tmp(1,:,:,ispden),ldim,hdp2,ldim,zero,hdp(:,:,ispden),ldim)
             noccmmp_tmp(1,:,:,ispden)=hdp(:,:,ispden)
           end do ! ispden
           ABI_DEALLOCATE(hdp2)
         else
           ABI_ALLOCATE(zhdp,(2*ldim,2*ldim))
           ABI_ALLOCATE(zhdp2,(2*ldim,2*ldim))
           zhdp(:,:)=cmplx(hdp(:,:,1),zero,kind=dp)
           zhdp2(:,:)=cmplx(zero,zero,kind=dp)
           call zgemm('n','c',2*ldim,2*ldim,2*ldim,cone,zhdp,2*ldim,znoccmmp_tmp,2*ldim,czero,zhdp2,2*ldim)
           zhdp(:,:)=cmplx(zero,zero,kind=dp)
           call zgemm('n','n',2*ldim,2*ldim,2*ldim,cone,znoccmmp_tmp,2*ldim,zhdp2,2*ldim,czero,zhdp,2*ldim)
           znoccmmp_tmp=zhdp
           ABI_DEALLOCATE(zhdp)
           ABI_DEALLOCATE(zhdp2)
         end if
         nmat=ndij ; if(ndij==4.and.cplex_dij==2) nmat=1
         do ispden=1,nmat
           write(message,'(2a,i3,a)') ch10,'== Atom ',iatom_tot,&
&           ' == Rotated back diagonalized matrix'
           if (ndij==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
           if (ndij==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
           if (ndij==4.and.cplex_dij==2) write(message,fmt='(4a)')     trim(message)," for all component "
           call wrtout(std_out,message,wrt_mode)
           do ilm=1,ldim*cplex_dij
             if(ndij==1.or.ndij==2)&
&             write(message,'(12(1x,9(1x,f10.5)))')&
&             (noccmmp_tmp(1,ilm,jlm,ispden),jlm=1,ldim)
             if(ndij==4.and.cplex_dij==2)&
&             write(message,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))')&
&             (znoccmmp_tmp(ilm,jlm),jlm=1,ldim*cplex_dij)
             call wrtout(std_out,message,wrt_mode)
           end do
         end do ! ispden
       end if
       ABI_DEALLOCATE(hdp)

     end if ! dmatudiag_loc

!    ########################################################################################
!    # Impose value of nocc_mmp from dmatpu; symetrize it
!    ########################################################################################
     if (usepawu/=0.and.impose_dmat/=0) then

       lpawu=lcur
       nsploop=nsppol;if (ndij==4) nsploop=4
       noccsym_error=.false.

!      Loop over spin components
       do ispden=1,nsploop
         if (ndij/=4) then
           jspden=min(3-ispden,paw_ij(iatom)%nsppol)
         else if (ispden<=2) then
           jspden=3-ispden
         else
           jspden=ispden
         end if

!        Loops over components of nocc_mmp
         do jlm=1,2*lpawu+1
           do ilm=1,2*lpawu+1

             if(nsym>1.and.ndij<4) then

               nsym_used(1:2)=0
               sumocc(1:2)=zero

!              Accumulate values of nocc_mmp over symmetries
               do irot=1,nsym
                 if ((symafm(irot)/=1).and.(.not.use_afm)) cycle
                 kspden=ispden;if (symafm(irot)==-1) kspden=jspden
                 factafm=one;if (ispden>3) factafm=dble(symafm(irot))
                 iafm=1;if ((antiferro).and.(symafm(irot)==-1)) iafm=2
                 nsym_used(iafm)=nsym_used(iafm)+1
                 at_indx=indsym(4,irot,iatom_tot)
                 do im2=1,2*lpawu+1
                   do im1=1,2*lpawu+1
!                    Be careful: use here R_rel^-1 in term of spherical harmonics
!                    which is tR_rec in term of spherical harmonics
!                    so, use transpose[zarot]
                     sumocc(iafm)=sumocc(iafm)+factafm*tmp_noccmmp(at_indx)%value(1,im1,im2,kspden) &
&                     *pawang%zarot(im1,ilm,lpawu+1,irot)&
&                     *pawang%zarot(im2,jlm,lpawu+1,irot)
!                    sumocc(iafm)=sumocc(iafm)+factafm*tmp_noccmmp(at_indx)%value(im1,im2,kspden) &
!                    &                     *pawang%zarot(ilm,im1,lpawu+1,irot)&
!                    &                     *pawang%zarot(jlm,im2,lpawu+1,irot)
                   end do
                 end do
               end do ! End loop over symmetries

!              Store new values of nocc_mmp
               paw_ij(iatom)%noccmmp(1,ilm,jlm,ispden)=sumocc(1)/nsym_used(1)
               if (.not.noccsym_error)&
&               noccsym_error=(abs(paw_ij(iatom)%noccmmp(1,ilm,jlm,ispden) &
&               -tmp_noccmmp(iatom_tot)%value(1,ilm,jlm,ispden))>tol5)

!              Antiferromagnetic case: has to fill up "down" component of nocc_mmp
               if (antiferro.and.nsym_used(2)>0) paw_ij(iatom)%noccmmp(1,ilm,jlm,2)=sumocc(2)/nsym_used(2)

             else  ! nsym=1

!              Case without symetries
               paw_ij(iatom)%noccmmp(:,ilm,jlm,ispden)= tmp_noccmmp(iatom_tot)%value(:,ilm,jlm,ispden)
             end if

           end do !ilm
         end do !jlm
       end do ! ispden
       do ispden=1,nsploop
         paw_ij(iatom)%nocctot(ispden)=zero ! contains nmmp in the n m representation
         do im1=1,2*lcur+1
           if(ndij==4.and.ispden==1) then
!            in this case, on computes total number or electron for double counting correction
             paw_ij(iatom)%nocctot(ispden)=paw_ij(iatom)%nocctot(ispden)+&
&             paw_ij(iatom)%noccmmp(1,im1,im1,1)+paw_ij(iatom)%noccmmp(1,im1,im1,2)
           else if(ndij==4.and.ispden==2) then
             paw_ij(iatom)%nocctot(ispden)=paw_ij(iatom)%nocctot(ispden)+&
&             paw_ij(iatom)%noccmmp(1,im1,im1,3)+paw_ij(iatom)%noccmmp(1,im1,im1,4)
           else if(ndij==4.and.ispden==3) then
             paw_ij(iatom)%nocctot(ispden)=paw_ij(iatom)%nocctot(ispden)-&
&             paw_ij(iatom)%noccmmp(2,im1,im1,3)+paw_ij(iatom)%noccmmp(2,im1,im1,4)
           else if(ndij==4.and.ispden==4) then
             paw_ij(iatom)%nocctot(ispden)=paw_ij(iatom)%nocctot(ispden)+&
&             paw_ij(iatom)%noccmmp(2,im1,im1,1)-paw_ij(iatom)%noccmmp(2,im1,im1,2)
           else
             paw_ij(iatom)%nocctot(ispden)=paw_ij(iatom)%nocctot(ispden)+&
&             paw_ij(iatom)%noccmmp(1,im1,im1,ispden)
           end if
         end do
       end do ! ispden

!      Printing of new nocc_mmp
       do ispden=1,ndij
         if(dmatudiag_loc==2) then
           write(message,'(2a,i3,a)') ch10,'== Atom ',iatom_tot,&
&           ' == Imposed occupation matrix (in the basis of diagonalization!!)'
         else
           write(message,'(2a,i3,a)') ch10,'== Atom ',iatom_tot,&
&           ' == Imposed occupation matrix'
         end if
         if (ndij==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
         if (ndij==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
         if (ndij==4) write(message,fmt='(4a)')     trim(message)," for component ", &
&         trim(dspinc(ispden+2*(ndij/4)))," =="
         call wrtout(std_out,message,wrt_mode)
         do ilm=1,2*lpawu+1
           if(cplex_dij==1)&
&           write(message,'(12(1x,9(1x,f10.5)))')&
&           (paw_ij(iatom)%noccmmp(1,ilm,jlm,ispden),jlm=1,2*lpawu+1)
           if(cplex_dij==2)&
&           write(message,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))')&
&           (paw_ij(iatom)%noccmmp(:,ilm,jlm,ispden),jlm=1,2*lpawu+1)
           call wrtout(std_out,message,wrt_mode)
         end do
       end do

!      WARNING if symmetrization changes the matrix
       if (noccsym_error) then
         write(message, '(a,i4,6a)' ) &
         '   After symmetrization, imposed occupation matrix for atom ',iatom_tot,ch10,&
&         '   is different from dmatpawu value set in input file !',ch10,&
&         '   It is likely that dmatpawu does not match the symmetry operations of the system.',ch10,&
&         '   Action: change dmatpawu in input file or increase precision until 0.00001'
         MSG_WARNING(message)
       end if

     end if ! impose_dmat/=0

!    ########################################################################################
!    # Rotate imposed occupation matrix in the non-diagonal basis
!    ########################################################################################
     if (usepawu/=0.and.impose_dmat/=0.and.dmatudiag_loc==2) then

       lpawu=lcur;ldim=2*lpawu+1

!      Rotation of imposed nocc_mmp
       if (ndij/=4) then
         ABI_ALLOCATE(hdp2,(ldim,ldim))
         do ispden=1,ndij
           call dgemm('n','t',ldim,ldim,ldim,one,&
&           paw_ij(iatom)%noccmmp(1,:,:,ispden),ldim,noccmmp_tmp(1,:,:,ispden),ldim,zero,hdp2,ldim)
           call dgemm('n','n',ldim,ldim,ldim,one,&
&           noccmmp_tmp(1,:,:,ispden),ldim,hdp2,ldim,zero,paw_ij(iatom)%noccmmp(1,:,:,ispden),ldim)
         end do ! ispden
         ABI_DEALLOCATE(hdp2)
       else
         ABI_ALLOCATE(zhdp,(2*ldim,2*ldim))
         ABI_ALLOCATE(zhdp2,(2*ldim,2*ldim))
         do im2=1,ldim
           do im1=1,ldim
             zhdp(     im1,     im2)=cmplx(paw_ij(iatom)%noccmmp(1,im1,im2,1),zero,kind=dp)  ! to be checked
             zhdp(ldim+im1,ldim+im2)=cmplx(paw_ij(iatom)%noccmmp(1,im1,im2,2),zero,kind=dp)  ! to be checked
             zhdp(     im1,ldim+im2)=&
&             cmplx(paw_ij(iatom)%noccmmp(1,im1,im2,3),+paw_ij(iatom)%noccmmp(1,im1,im2,4),kind=dp)  ! to be checked
             zhdp(ldim+im1,     im2)=&
&             cmplx(paw_ij(iatom)%noccmmp(1,im2,im1,3),-paw_ij(iatom)%noccmmp(1,im2,im1,4),kind=dp)  ! to be checked
           end do
         end do
         call zgemm('n','c',2*ldim,2*ldim,2*ldim,cone,zhdp,2*ldim,znoccmmp_tmp,2*ldim,czero,zhdp2,2*ldim)
         call zgemm('n','n',2*ldim,2*ldim,2*ldim,cone,znoccmmp_tmp,2*ldim,zhdp2,2*ldim,czero,zhdp,2*ldim)
         do jlm=1,ldim
           do ilm=1,ldim
             paw_ij(iatom)%noccmmp(1,ilm,jlm,1)= real(znoccmmp_tmp(     ilm,     jlm))  ! to be checked
             paw_ij(iatom)%noccmmp(1,ilm,jlm,2)= real(znoccmmp_tmp(ldim+ilm,ldim+jlm))  ! to be checked
             paw_ij(iatom)%noccmmp(1,ilm,jlm,3)= real(znoccmmp_tmp(     ilm,ldim+jlm))  ! to be checked
             paw_ij(iatom)%noccmmp(1,ilm,jlm,4)=aimag(znoccmmp_tmp(     ilm,ldim+jlm))  ! to be checked
           end do
         end do
         ABI_DEALLOCATE(zhdp)
         ABI_DEALLOCATE(zhdp2)
       end if

!      Printing of rotated imposed matrix
       do ispden=1,ndij
         write(message,'(2a,i3,a)') ch10,'== Atom ',iatom_tot,&
&         ' == Imposed density matrix in original basis'
         if (ndij==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
         if (ndij==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
         if (ndij==4) write(message,fmt='(4a)')     trim(message)," for component ", &
&         trim(dspin(ispden+2*(ndij/4)))," =="
         call wrtout(std_out,message,wrt_mode)
         do ilm=1,2*lpawu+1
           write(message,'(12(1x,9(1x,f10.5)))') (paw_ij(iatom)%noccmmp(1,ilm,jlm,ispden),jlm=1,2*lpawu+1)  ! to be checked
           call wrtout(std_out,message,wrt_mode)
         end do
       end do ! ispden

     end if ! dmatudiag_loc==2

     if (usepawu/=0.and.dmatudiag_loc>0) then
       ABI_DEALLOCATE(noccmmp_tmp)
       if (ndij==4)  then
         ABI_DEALLOCATE(znoccmmp_tmp)
       end if
     end if

     paw_ij(iatom)%has_pawu_occ=2

   end if ! lcur
 end do ! iatom

!Memory deallocation
 if (usepawu/=0.and.impose_dmat/=0) then
   do iatom_tot=1,natom
     lpawu=pawtab(typat(iatom_tot))%lpawu
     if (lpawu/=-1)  then
       ABI_DEALLOCATE(tmp_noccmmp(iatom_tot)%value)
     end if
   end do
   ABI_DATATYPE_DEALLOCATE(tmp_noccmmp)
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 DBG_EXIT("COLL")

end subroutine setnoccmmp
!!***

!----------------------------------------------------------------------

!----------------------------------------------------------------------

!!****f* m_paw_correlations/setrhoijpbe0
!! NAME
!! setrhoijpbe0
!!
!! FUNCTION
!! PAW local exact exchange only:
!! Impose value of rhoij for f electrons using an auxiliairy file
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  initialized= if 0, the initialization of the gstate run is not yet finished
!!  istep=index of the number of steps in the routine scfcv
!!  istep_mix=index of the number of steps for the SCF mixing (can be <istep)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  mpi_comm_read=MPI communicator containing all the processes reading the PBE0 file
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell
!!  ntypat=number of types of atoms in unit cell
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  typat(natom)=type integer for each atom in cell
!!
!! SIDE EFFECTS
!!  istep_mix=index of the number of steps for the SCF mixing (can be <istep)
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!
!! NOTES
!!  Only valid for f electrons !!!
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawio_print_ij,wrtout,xmpi_sum
!!
!! SOURCE

subroutine setrhoijpbe0(dtset,initialized,istep,istep_mix,&
&                       mpi_comm_read,my_natom,natom,ntypat,pawrhoij,pawtab,typat,&
&                       mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: initialized,istep,mpi_comm_read,my_natom,natom,ntypat
 integer,intent(inout) :: istep_mix
 integer,optional,intent(in) :: comm_atom
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer,parameter :: ll=3
 integer :: cplex_rhoij,iatom,iatom_tot,ierr,ii,ios,iread,irhoij,ispden,itypat,jj
 integer :: klmn,my_comm_atom,my_rank,nselect,nstep1,nstep1_abs,rhoijshft,rhoijsz
 logical :: my_atmtab_allocated,paral_atom,test0
 character(len=9),parameter :: filnam='rhoijpbe0'
 character(len=9),parameter :: dspin(6)=(/"up       ","down     ","up-up    ","down-down","Re[up-dn]","Im[up-dn]"/)
 character(len=500) :: strg, message
!arrays
 integer, allocatable :: nspden_tmp(:)
 integer,pointer :: my_atmtab(:)
 real(dp),allocatable :: rhoijtmp(:,:),rhoijtmp1(:,:),rhoijtmp2(:,:,:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Some limitation
 if (my_natom>0) then
   if (pawrhoij(1)%qphase==2) then
     message='setrhoijpbe0 not compatible with qphase=2!'
     MSG_BUG(message)
   end if
 end if

!Test existence of file and open it
 inquire(file=filnam,iostat=ios,exist=test0)
 if(.not.test0) return

!Look for parallelisation over atomic sites
 paral_atom=(present(comm_atom).and.(my_natom/=natom))

!Test if exact-exch. is on f electrons
 test0=.false.
 do itypat=1,ntypat
   if (pawtab(itypat)%useexexch/=0.and.pawtab(itypat)%lexexch/=ll) test0=.true.
 end do
 if (test0) then
   write(message, '(3a,i1,a)' ) &
&   ' Local exact exchange: occ. matrix can only be imposed for l=',ll,' !'
   MSG_ERROR(message)
 end if

!============================================================
!===== First case: no parallelisation over atomic sites =====
!============================================================

 if (.not.paral_atom) then

!  Open file
   if (open_file(filnam,message,unit=77,form='formatted') /= 0) then
     MSG_ERROR(message)
   end if

!  Read step number and eventually exit
   nstep1=0;test0=.false.
   do while (.not.test0)
     read(77,'(A)') strg
     test0=(strg(1:1)/="#")
     if (test0) read(unit=strg,fmt=*) nstep1
   end do
   nstep1_abs=abs(nstep1)
   if (nstep1_abs==0.or.istep>nstep1_abs.or.(nstep1>0.and.initialized/=0)) then
     close(77)
!    Reinitalize mixing when rhoij is allowed to change; for experimental purpose...
     if (dtset%userib==1234.and.istep==1+nstep1_abs.and.(nstep1<0.or.initialized==0)) istep_mix=1
     return
   end if

!  Loop on atoms
   do iatom=1,natom
     itypat=typat(iatom)
     cplex_rhoij=pawrhoij(iatom)%cplex_rhoij

     if (pawtab(itypat)%useexexch/=0) then

!      Set sizes depending on ll
       rhoijsz=4*ll+2
       rhoijshft=2*ll*ll

!      Uncompress rhoij
       ABI_ALLOCATE(rhoijtmp,(pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
       do ispden=1,pawrhoij(iatom)%nspden
         rhoijtmp=zero
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=pawrhoij(iatom)%rhoijselect(irhoij)
           rhoijtmp(klmn,ispden)=pawrhoij(iatom)%rhoijp(irhoij,ispden)
         end do
       end do
!      Read rhoij from file
       ABI_ALLOCATE(rhoijtmp1,(rhoijsz,rhoijsz))
       do ispden=1,pawrhoij(iatom)%nspden
         do ii=1,rhoijsz
           test0=.false.
           do while (.not.test0)
             read(77,'(A)') strg
             test0=(strg(1:1)/="#")
             if (test0)  read(unit=strg,fmt=*) (rhoijtmp1(ii,jj), jj=1,rhoijsz)
           end do
         end do

!        Impose rhoij
         do jj=1,rhoijsz
           do ii=1,jj
             rhoijtmp((jj+rhoijshft)*((jj+rhoijshft)-1)/2+ii+rhoijshft,ispden)=rhoijtmp1(ii,jj)
           end do
         end do

       end do
       ABI_DEALLOCATE(rhoijtmp1)

!      Compress rhoij
       nselect=0 ; pawrhoij(iatom)%rhoijselect=0
       pawrhoij(iatom)%rhoijp=zero
       do klmn=1,pawrhoij(iatom)%lmn2_size
         if (any(abs(rhoijtmp(klmn,:))>tol10)) then
           nselect=nselect+1 ; ii=cplex_rhoij*(nselect-1)+1
           do ispden=1,pawrhoij(iatom)%nspden
             pawrhoij(iatom)%rhoijp(ii,ispden)=rhoijtmp(klmn,ispden)
           end do
           pawrhoij(iatom)%rhoijselect(nselect)=klmn
         end if
       end do
       pawrhoij(iatom)%nrhoijsel=nselect
       ABI_DEALLOCATE(rhoijtmp)

!      Print new rhoij
       do ispden=1,pawrhoij(iatom)%nspden
         write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,&
&         ' == Imposed occupation matrix'
         if (pawrhoij(iatom)%nspden==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
         if (pawrhoij(iatom)%nspden==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
         if (pawrhoij(iatom)%nspden==4) write(message,fmt='(4a)')     trim(message)," for component ", &
&         trim(dspin(ispden+2*(pawrhoij(iatom)%nspden/4)))," =="
         call wrtout(std_out,message,'COLL')
         call pawio_print_ij(std_out,pawrhoij(iatom)%rhoijp(:,ispden),pawrhoij(iatom)%nrhoijsel,&
&         pawrhoij(iatom)%cplex_rhoij,pawrhoij(iatom)%lmn_size,ll,&
&         pawtab(itypat)%indlmn(1,1:pawtab(itypat)%lmn_size),&
&         1,-1,pawrhoij(iatom)%rhoijselect(:),-1.d0,1,mode_paral='COLL')
       end do

!      End loop on atoms
     end if
   end do

!  Close file
   close (77)

 else

!  ============================================================
!  ====== 2nd case: no parallelisation over atomic sites =====
!  ============================================================

   my_rank=xmpi_comm_rank(mpi_comm_read)

!  Read step number and eventually exit
   iread=0
   if (my_rank==0) then
     if (open_file(filnam,message,unit=77,form='formatted') /=0 ) then
       MSG_ERROR(message)
     end if
     nstep1=0;test0=.false.
     do while (.not.test0)
       read(77,'(A)') strg
       test0=(strg(1:1)/="#")
       if (test0) read(unit=strg,fmt=*) nstep1
     end do
     nstep1_abs=abs(nstep1)
     if (nstep1_abs==0.or.istep>nstep1_abs.or.(nstep1>0.and.initialized/=0)) then
       close(77)
!      Reinitalize mixing when rhoij is allowed to change; for experimental purpose...
       if (dtset%userib==1234.and.istep==1+nstep1_abs.and.(nstep1<0.or.initialized==0)) istep_mix=1
       iread=1
     end if
   end if
   call xmpi_sum(iread,mpi_comm_read,ierr)
   if (iread/=0) return

!  Set up parallelism over atoms
   nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
   my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
   call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!  Store number of component for rhoij
   ABI_ALLOCATE(nspden_tmp,(natom))
   nspden_tmp(:)=zero
   do iatom=1,my_natom
     iatom_tot=my_atmtab(iatom)
     nspden_tmp(iatom_tot)=pawrhoij(iatom)%nspden
   end do
   call xmpi_sum(nspden_tmp,mpi_comm_read,ierr)

!  To be improve if too much memory
   ABI_ALLOCATE(rhoijtmp2,(natom,4,rhoijsz,rhoijsz))
   rhoijtmp2=zero

!  Read rhoij from file
   if (my_rank==0) then
     do iatom=1,natom
       itypat=typat(iatom)
       if (pawtab(itypat)%useexexch/=0) then
         rhoijsz=4*ll+2
         do ispden=1,nspden_tmp(iatom)
           do ii=1,rhoijsz
             test0=.false.
             do while (.not.test0)
               read(77,'(A)') strg
               test0=(strg(1:1)/="#")
               if (test0)  read(unit=strg,fmt=*) (rhoijtmp2(iatom,ispden,ii,jj),jj=1,rhoijsz)
             end do
           end do
         end do
       end if
     end do
   end if
   call xmpi_sum(rhoijtmp2,mpi_comm_read,ierr)

!  Now, distribute rhoij
   do iatom=1,my_natom
     iatom_tot=my_atmtab(iatom)
     itypat=pawrhoij(iatom)%itypat
     cplex_rhoij=pawrhoij(iatom)%cplex_rhoij

     if (pawtab(itypat)%useexexch/=0) then

!      Set sizes depending on ll
       rhoijsz=4*ll+2
       rhoijshft=2*ll*ll

!      Uncompress rhoij
       ABI_ALLOCATE(rhoijtmp,(pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
       do ispden=1,pawrhoij(iatom)%nspden
         rhoijtmp=zero
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=pawrhoij(iatom)%rhoijselect(irhoij)
           rhoijtmp(klmn,ispden)=pawrhoij(iatom)%rhoijp(irhoij,ispden)
         end do

!        Impose rhoij
         do jj=1,rhoijsz
           do ii=1,jj
             rhoijtmp((jj+rhoijshft)*((jj+rhoijshft)-1)/2+ii+rhoijshft,ispden)=rhoijtmp2(iatom_tot,ispden,ii,jj)
           end do
         end do

       end do

!      Compress rhoij
       nselect=0 ; pawrhoij(iatom)%rhoijselect=0
       pawrhoij(iatom)%rhoijp=zero
       do klmn=1,pawrhoij(iatom)%lmn2_size
         if (any(abs(rhoijtmp(klmn,:))>tol10)) then
           nselect=nselect+1 ; ii=cplex_rhoij*(nselect-1)+1
           do ispden=1,pawrhoij(iatom)%nspden
             pawrhoij(iatom)%rhoijp(ii,ispden)=rhoijtmp(klmn,ispden)
           end do
           pawrhoij(iatom)%rhoijselect(nselect)=klmn
         end if
       end do
       pawrhoij(iatom)%nrhoijsel=nselect
       ABI_DEALLOCATE(rhoijtmp)

     end if ! useexexch/=0

!    Print new rhoij
     do ispden=1,pawrhoij(iatom)%nspden
       write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,' == Imposed occupation matrix'
       if (pawrhoij(iatom)%nspden==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
       if (pawrhoij(iatom)%nspden==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
       if (pawrhoij(iatom)%nspden==4) write(message,fmt='(4a)')     trim(message)," for component ", &
&       trim(dspin(ispden+2*(pawrhoij(iatom)%nspden/4)))," =="
       call wrtout(std_out,message,'PERS')
       call pawio_print_ij(std_out,pawrhoij(iatom)%rhoijp(:,ispden),pawrhoij(iatom)%nrhoijsel,&
&       pawrhoij(iatom)%cplex_rhoij,pawrhoij(iatom)%lmn_size,ll,&
&       pawtab(itypat)%indlmn(1,1:pawtab(itypat)%lmn_size),&
&       1,-1,pawrhoij(iatom)%rhoijselect(:),-1.d0,1,mode_paral='PERS')
     end do

!    end loop on atoms
   end do

   ABI_DEALLOCATE(nspden_tmp)
   ABI_DEALLOCATE(rhoijtmp2)

!  Destroy atom table used for parallelism
   call free_my_atmtab(my_atmtab,my_atmtab_allocated)

!  ============================================================
 end if ! paral_atom

 DBG_EXIT("COLL")

end subroutine setrhoijpbe0
!!***

!----------------------------------------------------------------------

!!****f* m_paw_correlations/calc_ubare
!! NAME
!! calc_ubare
!!
!! FUNCTION
!! Calculate the bare interaction on atomic orbitals
!!
!! INPUTS
!!  itypatcor = value of itypat for correlated species
!!  lpawu = angular momentum for correlated species
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data:
!!  pawang
!!     %lmax=Maximum value of angular momentum l+1
!!     %gntselect((2*l_max-1)**2,l_max**2,l_max**2)=
!!                     selection rules for Gaunt coefficients
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data:
!!     %mesh_size=Dimension of radial mesh
!!     %rad(mesh_size)=The coordinates of all the points of the radial mesh
!!     %radfact(mesh_size)=Factor used to compute radial integrals
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      pawpuxinit
!!
!! CHILDREN
!!      poisson,simp_gen,wrtout
!!
!! SOURCE

 subroutine calc_ubare(itypatcor,lpawu,pawang,pawrad,pawtab,rmax)

!Arguments ------------------------------------
 integer, intent(in)   :: itypatcor,lpawu
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),target,intent(in) :: pawtab
 real(dp), optional, intent(in) :: rmax

!Local variables ------------------------------
!scalars
 integer :: ilmn,ilmn1,iln,iln1,isel,isel1,itypat,jlmn,jlmn1,jln,jln1
 integer :: klm,klm1,klmn,klmn1,ll,lm0
 integer :: lmin,lmax,lmn2_size,mesh_size,meshsz,mm
 real(dp) :: norm,r_for_intg,rg,rg1,ubare,uint,uint_tmp
 character(len=800) :: message
!arrays
 real(dp),allocatable :: ff(:),gg(:),phiphj(:),phiphj1(:)

!************************************************************************

 itypat=itypatcor
 write(message,'(11a,f12.4,2a,i7,2a,f12.4,2a,i7,2a,f12.4)') &
& ch10," =======================================================================",ch10, &
& "  == Calculation of diagonal bare Coulomb interaction on ATOMIC orbitals ",ch10, &
& "     (it is assumed that the wavefunction for the first reference ",ch10, &
& "             energy in PAW atomic data is an atomic eigenvalue)",ch10,ch10, &
& " Max value of the radius in atomic data file   =", pawrad%rmax ,ch10, &
& " Max value of the mesh   in atomic data file   =", pawrad%mesh_size,ch10, &
& " PAW radius is                                 =", pawtab%rpaw,ch10, &
& " PAW value of the mesh for integration is      =", pawrad%int_meshsz,ch10, &
& " Integral of atomic wavefunction until rpaw    =", pawtab%ph0phiint(1)
 if(.not.present(rmax)) then
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 mesh_size=pawrad%mesh_size

!  Definition of the mesh used for integration.
 if(present(rmax)) then
   if(rmax>pawrad%rmax)  then
     write(message, '(a)' ) 'calc_ubare: the radius cannot be larger than the maximum radius of the mesh'
     MSG_ERROR(message)
   end if
   meshsz=pawrad_ifromr(pawrad,rmax)+5
   r_for_intg=rmax
 else
   meshsz=pawtab%partialwave_mesh_size
   r_for_intg=pawrad%rad(meshsz)  ! (we could use r_for_intg=-1)
 end if

 lmn2_size=pawtab%lmn2_size
 ABI_ALLOCATE(ff,(mesh_size))
 ABI_ALLOCATE(gg,(mesh_size))
 ABI_ALLOCATE(phiphj,(mesh_size))
 ABI_ALLOCATE(phiphj1,(mesh_size))
 do klmn=1,lmn2_size
   ilmn=pawtab%indklmn(7,klmn);jlmn=pawtab%indklmn(8,klmn)
   ! Select lpawu and first projectors il=jl=lpawu and first proj only
   if (( pawtab%indklmn(3,klmn)+pawtab%indklmn(4,klmn)==2*lpawu).and. &
&   (-pawtab%indklmn(3,klmn)+pawtab%indklmn(4,klmn)==2*lpawu).and. &
&   (pawtab%indlmn(3,ilmn)==1).and.(pawtab%indlmn(3,jlmn)==1) ) then
     klm=pawtab%indklmn(1,klmn);iln=pawtab%indlmn(5,ilmn);jln=pawtab%indlmn(5,jlmn)
     lmin=pawtab%indklmn(3,klmn);lmax=pawtab%indklmn(4,klmn)
     phiphj(1:meshsz)=pawtab%phi(1:meshsz,iln)*pawtab%phi(1:meshsz,jln)
     !write(6,*) "A",klmn,pawtab%klmntomn(1,klmn),pawtab%klmntomn(2,klmn),&
     !&pawtab%indklmn(7,klmn),pawtab%indklmn(8,klmn),pawtab%klmntomn(3,klmn),pawtab%klmntomn(4,klmn)
     do ll=lmin,lmin,2
       lm0=ll*ll+ll+1
       ff(1:meshsz)=phiphj(1:meshsz)
       call simp_gen(norm,ff,pawrad,r_for_intg=r_for_intg)
       call poisson(ff,ll,pawrad,gg)
       do klmn1=klmn,lmn2_size
         ilmn1=pawtab%indklmn(7,klmn);jlmn1=pawtab%indklmn(8,klmn)
         ! Select lpawu and first projectors il=jl=lpawu and first proj only
         if (( pawtab%indklmn(3,klmn1)+pawtab%indklmn(4,klmn1)==2*lpawu).and. &
&         (-pawtab%indklmn(3,klmn1)+pawtab%indklmn(4,klmn1)==2*lpawu).and. &
&         (pawtab%indlmn(3,ilmn1)==1).and.(pawtab%indlmn(3,jlmn1)==1) ) then
           !write(6,*) "A1",klmn1,pawtab%klmntomn(1,klmn1),pawtab%klmntomn(2,klmn1),&
           !&pawtab%indklmn(7,klmn1),pawtab%indklmn(8,klmn1),pawtab%klmntomn(3,klmn1),pawtab%klmntomn(4,klmn1)
           klm1=pawtab%indklmn(1,klmn1);iln1=pawtab%indlmn(5,ilmn1);jln1=pawtab%indlmn(5,jlmn1)
           phiphj1(1:meshsz)=pawtab%phi(1:meshsz,iln1)*pawtab%phi(1:meshsz,jln1)
           uint_tmp=zero
           if ((ll==lmin)) then
             ff(1)=zero
             ff(2:meshsz)=phiphj1(2:meshsz)*gg(2:meshsz)*two/pawrad%rad(2:meshsz)
             call simp_gen(uint_tmp,ff,pawrad,r_for_intg=r_for_intg)
           end if
           uint=zero
           do mm=-ll,ll
             isel =pawang%gntselect(lm0+mm,klm)
             isel1=pawang%gntselect(lm0+mm,klm1)
             if (isel>0.and.isel1>0) then
               rg =pawang%realgnt(isel)
               rg1=pawang%realgnt(isel1)
               uint=uint+uint_tmp*rg*rg1*two_pi
             end if
           end do
           if((pawtab%indklmn(5,klmn)==pawtab%indklmn(6,klmn)).and.&
&           (pawtab%indklmn(5,klmn1)==pawtab%indklmn(6,klmn1)).and.&
&           (pawtab%indklmn(5,klmn)==pawtab%indklmn(5,klmn1))) then
             ubare=uint*Ha_eV
           end if
         end if
       end do
     end do
   end if
 end do
 ABI_DEALLOCATE(gg)
 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(phiphj)
 ABI_DEALLOCATE(phiphj1)

 write(message,'(a,3(a,f12.4,a),2a,f12.4,a)') ch10," For an atomic wfn truncated at rmax =",r_for_intg,ch10,&
& "     The norm of the wfn is                    =",norm,ch10,&
& "     The bare interaction (no renormalization) =",ubare," eV",ch10,&
& "     The bare interaction (for a renorm. wfn ) =",ubare/norm/norm," eV"
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 if(r_for_intg < 10_dp .and. .not.present(rmax)) then
   write(message,'(a,f6.2,4a)') '   ( WARNING: The radial mesh in the atomic data file is cut at',r_for_intg,ch10,&
&   '   Use XML atomic data files to compute the bare Coulomb interaction',ch10,&
&   '   on a true normalized atomic wavefunction )'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if
 if(present(rmax)) then
   write(message,'(2a)')  " =======================================================================",ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 end subroutine calc_ubare
!!***

!----------------------------------------------------------------------

END MODULE m_paw_correlations
!!***
