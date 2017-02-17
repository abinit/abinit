!{\src2tex{textfont=tt}}
!!****f* ABINIT/newocc
!! NAME
!! newocc
!!
!! FUNCTION
!! Compute new occupation numbers at each k point,
!! from eigenenergies eigen, according to the
!! smearing scheme defined by occopt (smearing width tsmear and
!! physical temperature tphysel),
!! with the constraint of number of valence electrons per unit cell nelect.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  eigen(mband*nkpt*nsppol)=eigenvalues (input or init to large number), hartree
!!  spinmagntarget=if differ from -99.99_dp, fix the magnetic moment (in Bohr magneton)
!!  mband=maximum number of bands
!!  nband(nkpt)=number of bands at each k point
!!  nelect=number of electrons per unit cell
!!  nkpt=number of k points
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occopt=option for occupancies
!!  prtvol=control print volume and debugging output
!!  stmbias=if non-zero, compute occupation numbers for STM (non-zero around the Fermi energy)
!!   NOTE : in this case, only fermie and occ are meaningful outputs.
!!  tphysel="physical" electronic temperature with FD occupations
!!  tsmear=smearing width (or temperature)
!!  wtk(nkpt)=k point weights
!!
!! OUTPUT
!!  doccde(maxval(nband(:))*nkpt*nsppol)=derivative of occupancies wrt
!!           the energy for each band and k point
!!  entropy= entropy associated with the smearing (adimensional)
!!  fermie= fermi energy (Hartree)
!!  occ(maxval(nband(:))*nkpt*nsppol)=occupancies for each band and k point
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      gstate,m_ebands,respfn,vtorho
!!
!! CHILDREN
!!      getnel,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine newocc(doccde,eigen,entropy,fermie,spinmagntarget,mband,nband,&
&  nelect,nkpt,nspinor,nsppol,occ,occopt,prtvol,stmbias,tphysel,tsmear,wtk)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'newocc'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_61_occeig, except_this_one => newocc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkpt,nspinor,nsppol,occopt,prtvol
 real(dp),intent(in) :: spinmagntarget,nelect,stmbias,tphysel,tsmear
 real(dp),intent(out) :: entropy,fermie
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),wtk(nkpt)
 real(dp),intent(out) :: doccde(mband*nkpt*nsppol)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol) !vz_i

!Local variables-------------------------------
 integer,parameter :: niter_max=120,nkpt_max=50,fake_unit=-666,option1=1
 integer :: cnt,cnt2,cnt3,ib,ii,ik,ikpt,is,isppol,nkpt_eff
 integer :: sign
 integer,allocatable :: nbandt(:)
 real(dp) :: dosdeltae,entropy_tmp,fermihi,fermilo,fermimid,fermimid_tmp
 real(dp) :: fermi_biased,maxocc
 real(dp) :: nelect_tmp,nelecthi,nelectlo,nelectmid,nelect_biased
 real(dp) :: entropyt(2),fermihit(2),fermilot(2),fermimidt(2),nelecthit(2)
 real(dp) :: nelectlot(2),nelectt(2),tsec(2)
 real(dp),allocatable :: doccdet(:),eigent(:),occt(:)
 character(len=500) :: message

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(74,1,tsec)

!Here treat the case where occopt does not correspond to a metallic occupation scheme
 if (occopt<3 .or. occopt>8) then
   write(message,'(a,i0,a)')' occopt= ',occopt,', a value not allowed in newocc.'
   MSG_BUG(message)
 end if ! test of metallic occopt

!Check whether nband is a constant for all k point and spin-pol
 do isppol=1,nsppol
   do ikpt=1,nkpt
     if(nband(ikpt+(isppol-1)*nkpt)/=nband(1))then
       write(message,'(3a,i0,a,i0,a,i0,a)')&
&       'The number of bands must be the same for all k-points ',ch10,&
&       'but nband(1)= ',nband(1),' is different of nband(',&
&       ikpt+(isppol-1)*nkpt,') = ',nband(ikpt+(isppol-1)*nkpt),'.'
       MSG_BUG(message)
     end if
   end do
 end do

!Check whether nelect is strictly positive
 if(nelect<=zero)then
   write(message,'(3a,es16.8,a)')&
&   'nelect must be a positive number, while ',ch10,&
&   'the calling routine ask nelect=',nelect,'.'
   MSG_BUG(message)
 end if

 maxocc=two/(nsppol*nspinor)
!Check whether nelect is coherent with nband (nband(1) is enough,
!since it was checked that nband is independent of k-point and spin-pol
 if( nelect > nband(1)*nsppol*maxocc )then
   write(message,'(3a,es16.8,a,i0,a,es16.8,a)' )&
&   'nelect must be smaller than nband*maxocc, while ',ch10,&
&   'the calling routine gives nelect= ',nelect,', nband= ',nband(1),' and maxocc= ',maxocc,'.'
   MSG_BUG(message)
 end if

!Use bissection algorithm to find fermi energy
!This choice is due to the fact that it will always give sensible
!result (because the answer is bounded, even if the smearing function
!is non-monotonic (which is the case for occopt=4 or 6)
!Might speed up it, if needed !

!Lowest and largest trial fermi energies, and corresponding number of electrons
!They are obtained from the smallest or largest eigenenergy, plus a range of
!energy that allows for complete occupation of all bands, or, on the opposite,
!for zero occupation of all bands (see getnel.f)
 dosdeltae=zero  ! the DOS is not computed, with option=1
 fermilo=minval(eigen(1:nband(1)*nkpt*nsppol))-6.001_dp*tsmear
 if(occopt==3)fermilo=fermilo-24.0_dp*tsmear

 call getnel(doccde,dosdeltae,eigen,entropy,fermilo,maxocc,mband,nband,&
& nelectlo,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk)

 fermihi=maxval(eigen(1:nband(1)*nkpt*nsppol))+6.001_dp*tsmear
!safety value
 fermihi = min(fermihi, 1.e6_dp)
 if(occopt==3)fermihi=fermihi+24.0_dp*tsmear

 call getnel(doccde,dosdeltae,eigen,entropy,fermihi,maxocc,mband,nband,&
& nelecthi,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk)

!Prepare fixed moment calculation
 if(abs(spinmagntarget+99.99_dp)>1.0d-10)then
   sign = 1
   do is = 1, nsppol
     fermihit(is) = fermihi
     fermilot(is) = fermilo
     nelectt(is) = half*(nelect+sign*spinmagntarget)
     sign = -sign
     nelecthit(is) = nelecthi
     nelectlot(is) = nelectlo
   end do
 end if

!If the target nelect is not between nelectlo and nelecthi, exit
 if(nelect<nelectlo .or. nelect>nelecthi)then
   write(message, '(a,a,a,a,d16.8,a,a,d16.8,a,d16.8,a,a,d16.8,a,d16.8)') ch10,&
&   ' newocc : ',ch10,&
&   '  The calling routine gives nelect=',nelect,ch10,&
&   '  The lowest bound is ',fermilo,', with nelect=',nelectlo,ch10,&
&   '  The highest bound is ',fermihi,', with nelect=',nelecthi
   call wrtout(std_out,message,'COLL')

   write(message, '(11a)' )&
&   'In order to get the right number of electrons,',ch10,&
&   'it seems that the Fermi energy must be outside the range',ch10,&
&   'of eigenenergies, plus 6 or 30 times the smearing, which is strange.',ch10,&
&   'It might be that your number of bands (nband) corresponds to the strictly',ch10,&
&   'minimum number of bands to accomodate your electrons (so, OK for an insulator),',ch10,&
&   'while you are trying to describe a metal. In this case, increase nband, otherwise ...'
   MSG_BUG(message)
 end if

 if( abs(spinmagntarget+99.99_dp) < tol10 ) then

!  Usual bissection loop
   do ii=1,niter_max
     fermimid=(fermihi+fermilo)*half
!    Produce nelectmid from fermimid
     call getnel(doccde,dosdeltae,eigen,entropy,fermimid,maxocc,mband,nband,&
&     nelectmid,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk)
!    write(std_out,'(a,es24.16,a,es24.16)' )' newocc : from fermi=',fermimid,', getnel gives nelect=',nelectmid
     if(nelectmid>nelect*(one-tol14))then
       fermihi=fermimid
       nelecthi=nelectmid
     end if
     if(nelectmid<nelect*(one+tol14))then
       fermilo=fermimid
       nelectlo=nelectmid
     end if
     if( abs(nelecthi-nelectlo) <= nelect*two*tol14 .or. &
&     abs(fermihi-fermilo) <= tol14*abs(fermihi+fermilo) ) exit
     if(ii==niter_max)then
       write(message,'(a,i0,3a,es22.14,a,es22.14,a)')&
&       'It was not possible to find Fermi energy in ',niter_max,' bissections.',ch10,&
&       'nelecthi = ',nelecthi,', and nelectlo = ',nelectlo,'.'
       MSG_BUG(message)
     end if
   end do ! End of bissection loop

   fermie=fermimid
   write(message, '(a,f14.6,a,f14.6,a,a,i4)' ) &
&   ' newocc: new Fermi energy is ',fermie,' , with nelect=',nelectmid,ch10,&
&   '  Number of bissection calls =',ii
   call wrtout(std_out,message,'COLL')

!  Compute occupation numbers for prtstm/=0, close to the Fermi energy
   if(abs(stmbias)>tol10)then
     fermi_biased=fermie-stmbias
     ABI_ALLOCATE(occt,(mband*nkpt*nsppol))
     call getnel(doccde,dosdeltae,eigen,entropy,fermi_biased,maxocc,mband,nband,&
&     nelect_biased,nkpt,nsppol,occt,occopt,option1,tphysel,tsmear,fake_unit,wtk)
     occ(:)=occ(:)-occt(:)
     nelect_biased=abs(nelectmid-nelect_biased)
!    Here, arrange to have globally positive occupation numbers,
!    irrespective of the stmbias sign
     if(-stmbias>tol10)occ(:)=-occ(:)
     ABI_DEALLOCATE(occt)

     write(message,'(a,f14.6)')' newocc : the number of electrons in the STM range is nelect_biased=',nelect_biased
     call wrtout(std_out,message,'COLL')

   end if

 else ! Calculations with a specified moment

!  Bissection loop
   cnt2=0
   cnt3=0
   entropy=zero
   maxocc=one
   ABI_ALLOCATE(doccdet,(nkpt*mband))
   ABI_ALLOCATE(eigent,(nkpt*mband))
   ABI_ALLOCATE(occt,(nkpt*mband))
   ABI_ALLOCATE(nbandt,(nkpt))

   do is = 1, nsppol
     nelect_tmp = nelectt(is)
     fermihi = fermihit(is)
     fermilo = fermilot(is)
     nelecthi = nelecthit(is)
     nelectlo = nelectlot(is)
!    DEBUG
!    write(std_out,'(a,i1,3(f8.4,1x))') "Spin, N(spin):", is, nelect, fermihi, fermilo
!    write(std_out,'(a,2(f8.4,1x))') "Hi, lo:", nelecthi, nelectlo
!    ENDDEBUG

     do ii=1,niter_max
       fermimid_tmp=(fermihi+fermilo)/2.0_dp
!      temporary arrays
       cnt = 0
       do ik = 1, nkpt
         nbandt(ik) = mband
         do ib = 1, mband
           cnt = cnt + 1
           eigent(cnt) = eigen(cnt+cnt2)
           occt(cnt) = occ(cnt+cnt2)
           doccdet(cnt) = doccde(cnt+cnt2)
         end do
       end do

!      Produce nelectmid from fermimid
       call getnel(doccdet,dosdeltae,eigent,entropy_tmp,fermimid_tmp,maxocc,mband,nbandt,&
&       nelectmid,nkpt,1,occt,occopt,option1,tphysel,tsmear,fake_unit,wtk)
       entropyt(is) = entropy_tmp
       fermimidt(is) = fermimid_tmp
       fermimid = fermimidt(is)
!      temporary arrays
       cnt = 0
       do ik = 1, nkpt
         do ib = 1, mband
           cnt = cnt + 1
           occ(cnt+cnt2) = occt(cnt)
           doccde(cnt+cnt2) = doccdet(cnt)
         end do
       end do

!      DEBUG
!      write(std_out,'(a,es24.16,a,es24.16)' )&
!      &    ' newocc : from fermi=',fermimid,', getnel gives nelect=',nelectmid
!      ENDDEBUG

       if(nelectmid>=nelect_tmp)then
         fermihi=fermimid_tmp
         nelecthi=nelectmid
       else
         fermilo=fermimid_tmp
         nelectlo=nelectmid
       end if
       if( abs(nelecthi-nelectlo) <= 1.0d-13 .or. abs(fermihi-fermilo) <= 0.5d-14*abs(fermihi+fermilo) ) exit

       if(ii==niter_max)then
         write(message,'(a,i3,3a,es22.14,a,es22.14,a)')&
&         '  It was not possible to find Fermi energy in ',niter_max,' bissections.',ch10,&
&         '  nelecthi=',nelecthi,', and nelectlo=',nelectlo,'.'
         MSG_BUG(message)
       end if
     end do ! End of bissection loop

     cnt2 = cnt2 + nkpt*mband
     entropy = entropy + entropyt(is)
     fermie=fermimid
     write(message, '(a,i2,a,f14.6,a,f14.6,a,a,i4)' ) &
&     ' newocc : new Fermi energy for spin ', is, ' is ',fermie,' , with nelect=',nelectmid,ch10,&
&     '  Number of bissection calls =',ii
     call wrtout(std_out,message,'COLL')

   end do ! spin

   ABI_DEALLOCATE(doccdet)
   ABI_DEALLOCATE(eigent)
   ABI_DEALLOCATE(nbandt)
   ABI_DEALLOCATE(occt)

 end if !  End of logical on fixed moment calculations

!write(std_out,*) "kT*Entropy:", entropy*tsmear

 nkpt_eff=nkpt
 if(prtvol==0)nkpt_eff=min(nkpt_max,nkpt)

 if(nsppol==1)then
   write(message, '(a,i0,a)' ) &
&   ' newocc : computed new occ. numbers for occopt= ',occopt,&
&   ' , spin-unpolarized case. '
   call wrtout(std_out,message,'COLL')
   do ikpt=1,nkpt_eff
     write(message,'(a,i4,a)' ) ' k-point number ',ikpt,' :'
     do ii=0,(nband(1)-1)/12
       write(message,'(12f6.3)') &
&       occ(1+ii*12+(ikpt-1)*nband(1):min(12+ii*12,nband(1))+(ikpt-1)*nband(1))
       call wrtout(std_out,message,'COLL')
     end do
   end do
   if(nkpt/=nkpt_eff)then
     call wrtout(std_out,' newocc: prtvol=0, stop printing more k-point information','COLL')
   end if

!  DEBUG
!  write(message, '(a)' ) &
!  &   ' newocc : corresponding derivatives are '
!  call wrtout(std_out,message,'COLL')
!  do ikpt=1,nkpt_eff
!  write(message,'(a,i4,a)' ) ' k-point number ',ikpt,' :'
!  do ii=0,(nband(1)-1)/12
!  write(message,'(12f6.1)') &
!  &    doccde(1+ii*12+(ikpt-1)*nband(1):min(12+ii*12,nband(1))+(ikpt-1)*nband(1))
!  call wrtout(std_out,message,'COLL')
!  end do
!  end do
!  if(nkpt/=nkpt_eff)then
!  write(message,'(a)') &
!  &    ' newocc : prtvol=0, stop printing more k-point informations'
!  call wrtout(std_out,message,'COLL')
!  end if
!  ENDDEBUG
 else
   write(message, '(a,i0,a,a)' ) &
&   ' newocc : computed new occupation numbers for occopt= ',occopt,&
&   ch10,'  (1) spin up   values  '
   call wrtout(std_out,message,'COLL')
   do ikpt=1,nkpt_eff
     write(message,'(a,i0,a)' ) ' k-point number ',ikpt,':'
     do ii=0,(nband(1)-1)/12
       write(message,'(12f6.3)') &
&       occ(1+ii*12+(ikpt-1)*nband(1):min(12+ii*12,nband(1))+(ikpt-1)*nband(1))
       call wrtout(std_out,message,'COLL')
     end do
   end do
   if(nkpt/=nkpt_eff)then
     call wrtout(std_out,'newocc: prtvol=0, stop printing more k-point information','COLL')
   end if

   call wrtout(std_out,'  (2) spin down values  ','COLL')
   do ikpt=1,nkpt_eff
     do ii=0,(nband(1)-1)/12
       write(message,'(12f6.3)') &
&       occ( 1+ii*12+(ikpt-1+nkpt)*nband(1) : &
&       min(12+ii*12,nband(1))+(ikpt-1+nkpt)*nband(1) )
       call wrtout(std_out,message,'COLL')
     end do
   end do
   if(nkpt/=nkpt_eff)then
     call wrtout(std_out,' newocc: prtvol=0, stop printing more k-point information','COLL')
   end if

 end if !  End choice based on spin

 call timab(74,2,tsec)

 DBG_EXIT("COLL")

end subroutine newocc
!!***
