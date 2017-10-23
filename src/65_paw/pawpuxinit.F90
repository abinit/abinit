!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawpuxinit
!! NAME
!! pawpuxinit
!!
!! FUNCTION
!! Initialize some starting values of several arrays used in
!! PAW+U/+DMFT or local exact-exchange calculations
!!
!! A-define useful indices for LDA+U/local exact-exchange
!! B-Compute overlap between atomic wavefunction
!! C-Compute matrix elements of coulomb interaction (see PRB vol.52 5467)
!!    (angular part computed from Gaunt coefficients)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (BA,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors.
!!
!! INPUTS
!!  dmatpuopt= select expression for the density matrix
!!  exchmix= mixing factor for local exact-exchange
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
!!  usepawu= 0 if no LDA+U; 1 if LDA+U
!!
!! OUTPUT
!!  pawtab <type(pawtab_type)>=paw tabulated data read at start:
!!     %ij_proj=nproj*(nproju+1)/2
!!     %klmntomn(4,lmn2_size) = Array giving im, jm ,in, and jn for each klmn=(ilmn,jlmn)
!!     %lnproju(nproj)= value of ln for projectors on which paw+u/local exact-exchange acts.
!!     %nproju=number of projectors for orbitals on which paw+u/local exact-exchange acts.
!!     %phiphjint(pawtabitypat%ij_proj)=Integral of Phi(:,i)*Phi(:,j) for correlated orbitals.
!!     %usepawu=0 if no LDA+U; 1 if LDA+U
!!     %useexexch=0 if no local exact-exchange; 1 if local exact-exchange
!!     === if usepawu>0
!!     %jpawu= value of J
!!     %upawu= value of U
!!     %vee(2*lpawu+1,:,:,:)=matrix of the screened interaction for correlated orbitals
!!     === if useexexch>0
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawpuxinit(dmatpuopt,exchmix,f4of2_sla,f6of2_sla,jpawu,llexexch,llpawu,&
&           ntypat,pawang,pawprtvol,pawrad,pawtab,upawu,use_dmft,useexexch,usepawu,ucrpa)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_special_funcs

 use m_pawang, only : pawang_type
 use m_pawrad, only : pawrad_type, simp_gen, poisson
 use m_pawtab, only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawpuxinit'
 use interfaces_14_hidewrite
 use interfaces_65_paw, except_this_one => pawpuxinit
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: dmatpuopt,ntypat,pawprtvol,use_dmft,useexexch,usepawu
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
 integer :: icount,il,ilmn,isela,iselb,itemp,itypat,iu,j0lmn,jl,jlmn,ju,klm0u
 integer :: klm0x,klma,klmb,klmn,kln,kln1,kln2,kyc,lcur,lexexch,lkyc,ll,ll1
 integer :: lmexexch,lmkyc,lmn_size,lmn2_size,lmpawu,lpawu
 integer :: m1,m11,m2,m21,m3,m31,m4,m41
 integer :: mesh_size,int_meshsz,mkyc,sz1
 real(dp) :: ak,f4of2,f6of2,int1,intg,!testj,testu,testumj
 character(len=500) :: message
!arrays
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp),allocatable :: ff(:),fk(:),gg(:)

! *************************************************************************

 DBG_ENTER("COLL")

 if(useexexch==0.and.usepawu==0.and.use_dmft==0) return

!PAW+U and local exact-exchange restriction
 if(useexexch>0.and.usepawu>0)then
   do itypat=1,ntypat
     if (llpawu(itypat)/=llexexch(itypat).and.llpawu(itypat)/=-1.and.llexexch(itypat)/=-1) then
       write(message, '(5a,i2,3a)' )&
&       '  When PAW+U (usepawu>0) and local exact-exchange (exexch>0)',ch10,&
&       '  are selected together, they must apply on the same',ch10,&
&       '  angular momentum (lpawu/=lexexch forbidden, here for typat=',itypat,') !',ch10,&
&       '  Action: correct your input file.'
       MSG_ERROR(message)
     end if
   end do
 end if

!Print title
 if((usepawu>=1.and.usepawu<=4).or.useexexch>0) write(message, '(3a)' ) ch10,ch10," ******************************************"
 if(usepawu==1) then
   write(message, '(3a)' ) trim(message),ch10," LDA+U Method used: FLL"
 else if(usepawu==2) then
   write(message, '(3a)' ) trim(message),ch10," LDA+U Method used: AMF"
 else if(usepawu==3) then
   write(message, '(3a)' ) trim(message),ch10," LDA+U Method used: AMF (alternative)"
 else if(usepawu==4) then
   write(message, '(3a)' ) trim(message),ch10," LDA+U Method used: FLL with no spin polarization in the xc functional"
 end if
 if(useexexch>0) write(message, '(3a)' ) trim(message),ch10," PAW Local Exact exchange: PBE0"
 if((usepawu>=1.and.usepawu<=4).or.useexexch>0) &
 write(message, '(3a)' ) trim(message),ch10," ******************************************"
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
   if (usepawu>0.or.use_dmft>0) then
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
   if (useexexch>0) then
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
     if(useexexch>0.and.pawtab(itypat)%nproju>2)  then
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
         icount=icount+1
         itemp=(icount-1)/(2*lcur+1)
         if (itemp*(2*lcur+1)==icount-1) then
           pawtab(itypat)%lnproju(itemp+1)=indlmn(5,ilmn)
         end if
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
!    if (usepawu>0) then
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
         else if((dmatpuopt==2).or.(useexexch>0)) then
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
!    C-PAW+U: Matrix elements of coulomb interaction (see PRB vol.52 5467)
!    1. angular part computed from Gaunt coefficients
!    --------------------------------------------------------------------
     if (usepawu>0) then
       lpawu=lcur

!      a. compute F(k)
!      ---------------------------------------------
       ABI_ALLOCATE(fk,(lpawu+1))
       fk(1)=pawtab(itypat)%upawu
!      cf Slater Physical Review 165, p 665 (1968)
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

!!!!  pawtab(itypat)%vee(m11,m31,m21,m41)= <m11 m31| vee| m21 m41 >
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
     end if ! usepawu

!    ======================================================================
!    D-Local ex-exchange: Matrix elements of coulomb interaction and Fk
!    ----------------------------------------------------------------------
     if (useexexch>0) then
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

     end if !useexexch>0

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
