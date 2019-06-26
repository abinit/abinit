!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_occ
!! NAME
!! m_occ
!!
!! FUNCTION
!!  Low-level functions for occupation factors.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (XG, AF)
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

module m_occ

 use defs_basis
 use m_errors
 use m_abicore
 use m_splines
 use m_xmpi
 use m_hightemp

 use m_time,         only : timab
 use m_fstrings,     only : sjoin, itoa
 use defs_abitypes,  only : MPI_type
 use m_mpinfo,       only : proc_distrb_cycle

 implicit none

 private
!!***

 real(dp),parameter :: huge_tsmearinv = 1e50_dp
 real(dp),parameter :: maxFDarg=500.0_dp
 real(dp),parameter :: maxDFDarg=200.0_dp
 real(dp),parameter :: maxBEarg=600.0_dp
 real(dp),parameter :: maxDBEarg=200.0_dp

 public :: getnel        ! Compute total number of electrons from efermi or DOS
 public :: newocc        ! Compute new occupation numbers at each k point,
 public :: occeig        ! (occ_{k,q}(m)-occ_k(n))/(eig0_{k,q}(m)-eig0_k(n))$,
 public :: occ_fd        ! Fermi-Dirac statistic 1 / [(exp((e - mu)/ KT) + 1]
 public :: occ_dfd       ! Derivative of Fermi-Dirac statistic: (exp((e - mu)/ KT) / KT[(exp((e - mu)/ KT) + 1]^2
 public :: occ_be        ! Bose-Einstein statistic  1 / [(exp((e - mu)/ KT) - 1]
 public :: occ_dbe       ! Derivative of Bose-Einstein statistic  (exp((e - mu)/ KT) / KT[(exp((e - mu)/ KT) - 1]^2
 public :: dos_hdr_write
 public :: pareigocc

contains
!!***

!!****f* m_abinit/getnel
!! NAME
!! getnel
!!
!! FUNCTION
!! Option=1:
!!   Get the total number of electrons nelect, given a trial fermienergy fermie.
!!   For this, compute new occupation numbers at each k point,
!!   from eigenenergies eigen, according to the
!!   smearing scheme defined by occopt (and smearing width tsmear or tphysel).
!!
!! Option=2:
!!   Compute and output the smeared density of states, and the integrated density
!!   of states, then write these data
!!
!! Warning: this routine assumes checks have been done in the calling
!! routine, and that the values of the arguments are sensible
!!
!! NOTE
!! in order to speed the calculation, it would be easy to
!! compute the entropy only when the fermi energy is well converged
!!
!! INPUTS
!! dosdeltae= DOS delta of Energy (needed if Option=2)
!! eigen(mband*nkpt*nsppol)=eigenvalues (input or init to large number), hartree
!! fermie= fermi energy (Hartree)
!! maxocc=asymptotic maximum occupation number per band
!! mband=maximum number of bands
!! nband(nkpt*nsppol)=number of bands at each k point
!! nkpt=number of k points
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! occopt=option for occupancies, or re-smearing scheme if dblsmr /= 0
!! option=see above
!! tphysel="physical" electronic temperature with FD occupations
!! tsmear=smearing width (or temperature)
!! unitdos=unit number of output of the DOS. Not needed if option==1
!! wtk(nkpt)=k point weights
!!
!! OUTPUT
!! doccde(mband*nkpt*nsppol)=derivative of occupancies wrt the energy for each band and k point.
!! entropy= entropy associated with the smearing (adimensional)
!! nelect=number of electrons per unit cell
!! occ(mband*nkpt*nsppol)=occupancies for each band and k point.
!!
!! NOTES
!! Modified beginning 23/11/2000 by MV
!! Add an additional smearing on top of a FD type, in order to improve k-point
!! convergence: tsmear = 0 and tphysel ~= 2.e-3 corresponds to a small (300K)
!! temperature on the electrons insufficient for convergence purposes.
!! Feed re-smeared "Dirac delta" to the rest of ABINIT with only one parameter,
!! tphysel, which is the physical temperature.
!! encorr = correction to energy for terms of order tsmear^2:
!!       $  E_{phys} = E_{free} - encorr*(E_{int}-E_{free}) + O(tsmear^3)  $
!!
!! PARENTS
!!      cchi0q0_intraband,clnup1,conducti_nc,dfpt_looppert,m_ebands,newocc
!!
!! CHILDREN
!!      dos_hdr_write,init_occ_ent,splfit,wrtout
!!
!! SOURCE

subroutine getnel(doccde,dosdeltae,eigen,entropy,fermie,maxocc,mband,nband,&
&  nelect,nkpt,nsppol,occ,occopt,option,tphysel,tsmear,unitdos,wtk)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkpt,nsppol,occopt,option,unitdos
 real(dp),intent(in) :: dosdeltae,fermie,maxocc,tphysel,tsmear
 real(dp),intent(out) :: entropy,nelect
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),wtk(nkpt)
 real(dp),intent(out) :: doccde(mband*nkpt*nsppol) !vz_i
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol) !vz_i

!Local variables-------------------------------
! nptsdiv2 is the number of integration points, divided by 2.
! tratio  = ratio tsmear/tphysel for convoluted smearing function
! save values so we can impose recalculation of smdfun when
! the smearing or electronic temperature change between
! datasets
! corresponds roughly to delta_FD (maxFDarg) = 1.0d-100
!
! return fermi-dirac smearing function analytically
!
! real(dp) :: smdFD
! smdFD (tt) = 1.0_dp / (exp(-tt/2.0_dp) + exp(tt/2.0_dp))**2
!scalars
! TODO: This parameter is defined in init_occ_ent but we cannot call the
! routine to get this value since the same variable is used to dimension the
! arrays! This Constants should be stored somewhere in a module.
 integer,parameter :: nptsdiv2_def=6000,  prtdos1=1
 integer :: bantot,iband,iene,ikpt,index,index_start,isppol, nene,nptsdiv2
 real(dp) :: buffer,deltaene,dosdbletot,doshalftot,dostot
 real(dp) :: enemax,enemin,enex,intdostot,limit,tsmearinv
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: entfun(:,:),occfun(:,:)
 real(dp),allocatable :: smdfun(:,:),xgrid(:)
 real(dp),allocatable :: arg(:),derfun(:),dos(:),dosdble(:),doshalf(:),ent(:)
 real(dp),allocatable :: intdos(:)

! *************************************************************************

 DBG_ENTER("COLL")

 if(option/=1 .and. option/=2)then
   MSG_BUG(sjoin('Option must be either 1 or 2. It is:', itoa(option)))
 end if

!Initialize the occupation function and generalized entropy function,
!at the beginning, or if occopt changed

!Just get the number nptsdiv2 and allocate entfun, occfun, smdfun and xgrid accordingly
 nptsdiv2 = nptsdiv2_def

! call init_occ_ent(entfun, limit, &
!& nptsdiv2, occfun, occopt, -1, smdfun, tphysel, &
!& tsmear, tsmearinv, xgrid)

 ABI_ALLOCATE(entfun,(-nptsdiv2:nptsdiv2,2))
 ABI_ALLOCATE(occfun,(-nptsdiv2:nptsdiv2,2))
 ABI_ALLOCATE(smdfun,(-nptsdiv2:nptsdiv2,2))
 ABI_ALLOCATE(xgrid,(-nptsdiv2:nptsdiv2))

!Call to init_occ_ent
 call init_occ_ent(entfun, limit, nptsdiv2, occfun, occopt, option, smdfun, tphysel, &
& tsmear, tsmearinv, xgrid)

!The initialisation of occfun and entfun is done

!---------------------------------------------------------------------

!write(std_out,*)' getnel : debug  tphysel, tsmear = ', tphysel, tsmear
 bantot=sum(nband(:))

 ABI_ALLOCATE(arg,(bantot))
 ABI_ALLOCATE(derfun,(bantot))
 ABI_ALLOCATE(ent,(bantot))

 if(option==1)then
   ! normal evaluation of occupations and entropy

!  Compute the arguments of the occupation and entropy functions
!  HM 20/08/2018 Treat the T --> 0 limit
   if (tsmear==0) then
     arg(:)=sign(huge_tsmearinv,fermie-eigen(1:bantot))
   else
     arg(:)=(fermie-eigen(1:bantot))*tsmearinv
   endif

!  Compute the values of the occupation function, and the entropy function
!  Note : splfit also takes care of the points outside of the interval,
!  and assign to them the value of the closest extremal point,
!  which is what is needed here.

   call splfit(xgrid,doccde,occfun,1,arg,occ,(2*nptsdiv2+1),bantot)
   call splfit(xgrid,derfun,entfun,0,arg,ent,(2*nptsdiv2+1),bantot)

!  Normalize occ and ent, and sum number of electrons and entropy
   nelect=zero; entropy=zero
   index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       do iband=1,nband(ikpt+nkpt*(isppol-1))
         index=index+1
         ent(index)=ent(index)*maxocc
         occ(index)=occ(index)*maxocc
         doccde(index)=-doccde(index)*maxocc*tsmearinv
         entropy=entropy+wtk(ikpt)*ent(index)
         nelect=nelect+wtk(ikpt)*occ(index)
       end do
     end do
   end do

!  write(std_out,*) ' getnel : debug   wtk, occ, eigen = ', wtk, occ, eigen
!  write(std_out,*)xgrid(-nptsdiv2),xgrid(nptsdiv2)
!  write(std_out,*)'fermie',fermie
!  do ii=1,bantot
!  write(std_out,*)ii,arg(ii),doccde(ii)
!  end do
!  write(std_out,*)'eigen',eigen(:)
!  write(std_out,*)'arg',arg(:)
!  write(std_out,*)'occ',occ(:)
!  write(std_out,*)'nelect',nelect

 else if(option==2)then
  ! evaluate DOS for smearing, half smearing, and double.

   buffer=limit/tsmearinv*.5_dp

  ! A Similar section is present is dos_calcnwrite. Should move all DOS stuff to m_ebands
  ! Choose the lower and upper energies
   enemax=maxval(eigen(1:bantot))+buffer
   enemin=minval(eigen(1:bantot))-buffer

  ! Extend the range to a nicer value
   enemax=0.1_dp*ceiling(enemax*10._dp)
   enemin=0.1_dp*floor(enemin*10._dp)

  ! Choose the energy increment
   if(abs(dosdeltae)<tol10)then
     deltaene=0.001_dp
     if(prtdos1>=2)deltaene=0.0005_dp ! Higher resolution possible (and wanted) for tetrahedron
   else
     deltaene=dosdeltae
   end if
   nene=nint((enemax-enemin)/deltaene)+1

!  Write the header of the DOS file, and also decides the energy range and increment
   call dos_hdr_write(deltaene,eigen,enemax,enemin,fermie,mband,nband,nene,&
&   nkpt,nsppol,occopt,prtdos1,tphysel,tsmear,unitdos)

   ABI_ALLOCATE(dos,(bantot))
   ABI_ALLOCATE(dosdble,(bantot))
   ABI_ALLOCATE(doshalf,(bantot))
   ABI_ALLOCATE(intdos,(bantot))

   do isppol=1,nsppol

     if (nsppol==2) then
       if(isppol==1) write(msg,'(a,16x,a)')  '#','Spin-up DOS'
       if(isppol==2) write(msg,'(2a,16x,a)')  ch10,'#','Spin-dn DOS '
       call wrtout(unitdos,msg,'COLL')
     end if
     index_start=0
     if(isppol==2)then
       do ikpt=1,nkpt
         index_start=index_start+nband(ikpt)
       end do
     end if

     enex=enemin
     do iene=1,nene

!      Compute the arguments of the dos and occupation function
       arg(:)=(enex-eigen(1:bantot))*tsmearinv

       call splfit(xgrid,derfun,smdfun,0,arg,dos,(2*nptsdiv2+1),bantot)
       call splfit(xgrid,derfun,occfun,0,arg,intdos,(2*nptsdiv2+1),bantot)
!      Also compute the dos with tsmear halved and doubled
       arg(:)=arg(:)*2.0_dp
       call splfit(xgrid,derfun,smdfun,0,arg,doshalf,(2*nptsdiv2+1),bantot)
!      Since arg was already doubled, must divide by four
       arg(:)=arg(:)*0.25_dp
       call splfit(xgrid,derfun,smdfun,0,arg,dosdble,(2*nptsdiv2+1),bantot)

!      Now, accumulate the contribution from each eigenenergy
       dostot=zero
       intdostot=zero
       doshalftot=zero
       dosdbletot=zero
       index=index_start

!      write(std_out,*)' eigen, arg, dos, intdos, doshalf, dosdble'
       do ikpt=1,nkpt
         do iband=1,nband(ikpt+nkpt*(isppol-1))
           index=index+1
           dostot=dostot+wtk(ikpt)*maxocc*dos(index)*tsmearinv
           intdostot=intdostot+wtk(ikpt)*maxocc*intdos(index)
           doshalftot=doshalftot+wtk(ikpt)*maxocc*doshalf(index)*tsmearinv*2.0_dp
           dosdbletot=dosdbletot+wtk(ikpt)*maxocc*dosdble(index)*tsmearinv*0.5_dp
         end do
       end do

!      Print the data for this energy
       write(unitdos, '(f8.3,2f14.6,2f14.3)' )enex,dostot,intdostot,doshalftot,dosdbletot

       enex=enex+deltaene
     end do ! iene
   end do ! isppol

   ABI_DEALLOCATE(dos)
   ABI_DEALLOCATE(dosdble)
   ABI_DEALLOCATE(doshalf)
   ABI_DEALLOCATE(intdos)

!  MG: It does not make sense to close the unit here since the routines
!  did not open the file here!
!  Close the DOS file
   close(unitdos)
 end if

 ABI_DEALLOCATE(arg)
 ABI_DEALLOCATE(derfun)
 ABI_DEALLOCATE(ent)
 ABI_DEALLOCATE(entfun)
 ABI_DEALLOCATE(occfun)
 ABI_DEALLOCATE(smdfun)
 ABI_DEALLOCATE(xgrid)

 DBG_EXIT("COLL")

end subroutine getnel
!!***

!!****f* m_occ/newocc
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
!!   NOTE: in this case, only fermie and occ are meaningful outputs.
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
!! PARENTS
!!      gstate,m_ebands,respfn,vtorho
!!
!! CHILDREN
!!      getnel,timab,wrtout
!!
!! SOURCE

subroutine newocc(doccde,eigen,entropy,fermie,spinmagntarget,mband,nband,&
&  nelect,nkpt,nspinor,nsppol,occ,occopt,prtvol,stmbias,tphysel,tsmear,wtk,hightemp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkpt,nspinor,nsppol,occopt,prtvol
 real(dp),intent(in) :: spinmagntarget,nelect,stmbias,tphysel,tsmear
 real(dp),intent(out) :: entropy,fermie
 type(hightemp_type),intent(inout),optional :: hightemp
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),wtk(nkpt)
 real(dp),intent(out) :: doccde(mband*nkpt*nsppol)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol) !vz_i

!Local variables-------------------------------
 integer,parameter :: niter_max=120,nkpt_max=2,fake_unit=-666,option1=1
 integer :: cnt,cnt2,cnt3,ib,ii,ik,ikpt,is,isppol,nkpt_eff
 integer :: sign
 integer,allocatable :: nbandt(:)
 real(dp) :: dosdeltae,entropy_tmp,fermihi,fermilo,fermimid,fermimid_tmp
 real(dp) :: fermi_biased,maxocc
 real(dp) :: nelect_tmp,nelecthi,nelectlo,nelectmid,nelect_biased
 real(dp) :: entropyt(2),fermihit(2),fermilot(2),fermimidt(2),nelecthit(2)
 real(dp) :: nelectlot(2),nelectt(2),tsec(2)
 real(dp),allocatable :: doccdet(:),eigent(:),occt(:)
 character(len=500) :: msg

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(74,1,tsec)

 ! Here treat the case where occopt does not correspond to a metallic occupation scheme
 if (occopt<3 .or. occopt>8) then
   MSG_BUG(sjoin(' occopt= ',itoa(occopt),', a value not allowed in newocc.'))
 end if

 ! Check whether nband is a constant for all k point and spin-pol
 do isppol=1,nsppol
   do ikpt=1,nkpt
     if(nband(ikpt+(isppol-1)*nkpt)/=nband(1))then
       write(msg,'(3a,i0,a,i0,a,i0,a)')&
        'The number of bands must be the same for all k-points ',ch10,&
        'but nband(1)= ',nband(1),' is different of nband(',ikpt+(isppol-1)*nkpt,') = ',nband(ikpt+(isppol-1)*nkpt),'.'
       MSG_BUG(msg)
     end if
   end do
 end do

 ! Check whether nelect is strictly positive
 if(nelect <= zero)then
   write(msg,'(3a,es16.8,a)')&
&   'nelect must be a positive number, while ',ch10,&
&   'the calling routine asks nelect= ',nelect,'.'
   MSG_BUG(msg)
 end if

 maxocc=two/(nsppol*nspinor)
!Check whether nelect is coherent with nband (nband(1) is enough,
!since it was checked that nband is independent of k-point and spin-pol
 if (nelect > nband(1)*nsppol*maxocc) then
   write(msg,'(3a,es16.8,a,i0,a,es16.8,a)' )&
   'nelect must be smaller than nband*maxocc, while ',ch10,&
   'the calling routine gives nelect= ',nelect,', nband= ',nband(1),' and maxocc= ',maxocc,'.'
   MSG_BUG(msg)
 end if

!Use bisection algorithm to find fermi energy
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

!Blanchet - Compute the number of free electrons with corresponding chemical
!potential and add to nelect bounds.
 if(hightemp%enabled) then
   call hightemp_getnfreeel(hightemp%ebcut,fermilo,1024,nelect_tmp,tsmear,hightemp%u0,hightemp%ucvol)
   nelectlo=nelectlo+nelect_tmp
 end if

 fermihi=maxval(eigen(1:nband(1)*nkpt*nsppol))+6.001_dp*tsmear
!safety value
 fermihi = min(fermihi, 1.e6_dp)
 if(occopt==3)fermihi=fermihi+24.0_dp*tsmear

 call getnel(doccde,dosdeltae,eigen,entropy,fermihi,maxocc,mband,nband,&
& nelecthi,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk)

!Blanchet - Compute the number of free electrons with corresponding chemical
!potential and add to nelect bounds.
 if(hightemp%enabled) then
   call hightemp_getnfreeel(hightemp%ebcut,fermihi,1024,nelect_tmp,tsmear,hightemp%u0,hightemp%ucvol)
   nelecthi=nelecthi+nelect_tmp
 end if

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
   write(msg, '(a,a,a,a,d16.8,a,a,d16.8,a,d16.8,a,a,d16.8,a,d16.8)') ch10,&
&   ' newocc: ',ch10,&
&   '  The calling routine gives nelect=',nelect,ch10,&
&   '  The lowest bound is ',fermilo,', with nelect=',nelectlo,ch10,&
&   '  The highest bound is ',fermihi,', with nelect=',nelecthi
   call wrtout(std_out,msg,'COLL')

   write(msg, '(11a)' )&
&   'In order to get the right number of electrons,',ch10,&
&   'it seems that the Fermi energy must be outside the range',ch10,&
&   'of eigenenergies, plus 6 or 30 times the smearing, which is strange.',ch10,&
&   'It might be that your number of bands (nband) corresponds to the strictly',ch10,&
&   'minimum number of bands to accomodate your electrons (so, OK for an insulator),',ch10,&
&   'while you are trying to describe a metal. In this case, increase nband, otherwise ...'
   MSG_BUG(msg)
 end if

 if( abs(spinmagntarget+99.99_dp) < tol10 ) then

!  Usual bisection loop
   do ii=1,niter_max
     fermimid=(fermihi+fermilo)*half
!    Produce nelectmid from fermimid
     call getnel(doccde,dosdeltae,eigen,entropy,fermimid,maxocc,mband,nband,&
&     nelectmid,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk)
     !Blanchet - Compute the number of free electrons with corresponding chemical
     !potential and add to nelect bounds.
     if(hightemp%enabled) then
       call hightemp_getnfreeel(hightemp%ebcut,fermimid,1024,nelect_tmp,tsmear,hightemp%u0,hightemp%ucvol)
       nelectmid=nelectmid+nelect_tmp
     end if

!    write(std_out,'(a,es24.16,a,es24.16)' )' newocc: from fermi=',fermimid,', getnel gives nelect=',nelectmid
     if(nelectmid>nelect*(one-tol14))then
       fermihi=fermimid
       nelecthi=nelectmid
     end if
     if(nelectmid<nelect*(one+tol14))then
       fermilo=fermimid
       nelectlo=nelectmid
     end if
     if( abs(nelecthi-nelectlo) <= nelect*two*tol14 .or. abs(fermihi-fermilo) <= tol14*abs(fermihi+fermilo) ) exit
     if(ii==niter_max)then
       write(msg,'(a,i0,3a,es22.14,a,es22.14,a)')&
&       'It was not possible to find Fermi energy in ',niter_max,' bisections.',ch10,&
&       'nelecthi = ',nelecthi,', and nelectlo = ',nelectlo,'.'
       MSG_BUG(msg)
     end if
   end do ! End of bisection loop

   fermie=fermimid
   write(msg, '(2(a,f14.6),a,i0)' ) &
&   ' newocc: new Fermi energy is ',fermie,' , with nelect=',nelectmid,', Number of bisection calls: ',ii
   call wrtout(std_out,msg,'COLL')

!  Compute occupation numbers for prtstm/=0, close to the Fermi energy
   if(abs(stmbias)>tol10)then
     fermi_biased=fermie-stmbias
     ABI_ALLOCATE(occt,(mband*nkpt*nsppol))
     call getnel(doccde,dosdeltae,eigen,entropy,fermi_biased,maxocc,mband,nband,&
&     nelect_biased,nkpt,nsppol,occt,occopt,option1,tphysel,tsmear,fake_unit,wtk)
     occ(:)=occ(:)-occt(:)
     nelect_biased=abs(nelectmid-nelect_biased)
!    Here, arrange to have globally positive occupation numbers, irrespective of the stmbias sign
     if(-stmbias>tol10)occ(:)=-occ(:)
     ABI_DEALLOCATE(occt)

     write(msg,'(a,f14.6)')' newocc: the number of electrons in the STM range is nelect_biased=',nelect_biased
     call wrtout(std_out,msg,'COLL')
   end if

 else ! Calculations with a specified moment

!  Bisection loop
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
!    write(std_out,'(a,i1,3(f8.4,1x))') "Spin, N(spin):", is, nelect, fermihi, fermilo
!    write(std_out,'(a,2(f8.4,1x))') "Hi, lo:", nelecthi, nelectlo

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
!      write(std_out,'(a,es24.16,a,es24.16)' )' newocc: from fermi=',fermimid,', getnel gives nelect=',nelectmid

       if(nelectmid>=nelect_tmp)then
         fermihi=fermimid_tmp
         nelecthi=nelectmid
       else
         fermilo=fermimid_tmp
         nelectlo=nelectmid
       end if
       if( abs(nelecthi-nelectlo) <= 1.0d-13 .or. abs(fermihi-fermilo) <= 0.5d-14*abs(fermihi+fermilo) ) exit

       if(ii==niter_max)then
         write(msg,'(a,i3,3a,es22.14,a,es22.14,a)')&
&         'It was not possible to find Fermi energy in ',niter_max,' bisections.',ch10,&
&         'nelecthi= ',nelecthi,', and nelectlo= ',nelectlo,'.'
         MSG_BUG(msg)
       end if
     end do ! End of bisection loop

     cnt2 = cnt2 + nkpt*mband
     entropy = entropy + entropyt(is)
     fermie=fermimid
     write(msg, '(a,i2,a,f14.6,a,f14.6,a,a,i4)' ) &
&     ' newocc: new Fermi energy for spin ', is, ' is ',fermie,' , with nelect=',nelectmid,ch10,&
&     '  Number of bisection calls =',ii
     call wrtout(std_out,msg,'COLL')

   end do ! spin

   ABI_DEALLOCATE(doccdet)
   ABI_DEALLOCATE(eigent)
   ABI_DEALLOCATE(nbandt)
   ABI_DEALLOCATE(occt)

 end if !  End of logical on fixed moment calculations

!write(std_out,*) "kT*Entropy:", entropy*tsmear

 nkpt_eff=nkpt
 if(prtvol==0)nkpt_eff=min(nkpt_max,nkpt)

 if (nsppol==1)then
   write(msg, '(a,i0,a)' )' newocc: computed new occ. numbers for occopt= ',occopt,' , spin-unpolarized case. '
   call wrtout(std_out,msg,'COLL')
   do ikpt=1,nkpt_eff
     write(msg,'(a,i4,a)' ) ' k-point number ',ikpt,' :'
     do ii=0,(nband(1)-1)/12
       if (ii == 3 .and. prtvol /= 0) exit
       write(msg,'(12f6.3)') occ(1+ii*12+(ikpt-1)*nband(1):min(12+ii*12,nband(1))+(ikpt-1)*nband(1))
       call wrtout(std_out,msg,'COLL')
     end do
   end do
   if (nkpt/=nkpt_eff) call wrtout(std_out,' newocc: prtvol=0, stop printing more k-point information','COLL')

!  DEBUG
!  call wrtout(std_out,' newocc: corresponding derivatives are ','COLL')
!  do ikpt=1,nkpt_eff
!  write(msg,'(a,i4,a)' ) ' k-point number ',ikpt,' :'
!  do ii=0,(nband(1)-1)/12
!  write(msg,'(12f6.1)') doccde(1+ii*12+(ikpt-1)*nband(1):min(12+ii*12,nband(1))+(ikpt-1)*nband(1))
!  call wrtout(std_out,msg,'COLL')
!  end do
!  end do
!  if(nkpt/=nkpt_eff)then
!    call wrtout(std_out,'newocc: prtvol=0, stop printing more k-point information','COLL')
!  end if
!  ENDDEBUG
 else
   write(msg, '(a,i0,2a)' )' newocc: computed new occupation numbers for occopt= ',occopt,ch10,'  (1) spin up   values  '
   call wrtout(std_out,msg,'COLL')
   do ikpt=1,nkpt_eff
     write(msg,'(a,i0,a)' ) ' k-point number ',ikpt,':'
     do ii=0,(nband(1)-1)/12
       if (ii == 3 .and. prtvol /= 0) exit
       write(msg,'(12f6.3)') occ(1+ii*12+(ikpt-1)*nband(1):min(12+ii*12,nband(1))+(ikpt-1)*nband(1))
       call wrtout(std_out,msg,'COLL')
     end do
   end do
   if (nkpt/=nkpt_eff) call wrtout(std_out,'newocc: prtvol=0, stop printing more k-point information','COLL')

   call wrtout(std_out,'  (2) spin down values  ','COLL')
   do ikpt=1,nkpt_eff
     do ii=0,(nband(1)-1)/12
       if (ii == 3 .and. prtvol /= 0) exit
       write(msg,'(12f6.3)') occ( 1+ii*12+(ikpt-1+nkpt)*nband(1):min(12+ii*12,nband(1))+(ikpt-1+nkpt)*nband(1) )
       call wrtout(std_out,msg,'COLL')
     end do
   end do
   if(nkpt/=nkpt_eff) call wrtout(std_out,' newocc: prtvol=0, stop printing more k-point information','COLL')
 end if ! End choice based on spin

 call timab(74,2,tsec)

 DBG_EXIT("COLL")

end subroutine newocc
!!***

!!****f* m_occ/init_occ_ent
!! NAME
!! init_occ_ent
!!
!! FUNCTION
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! PARENTS
!!      getnel
!!
!! CHILDREN
!!      spline
!!
!! SOURCE

subroutine init_occ_ent(entfun,limit,nptsdiv2,occfun,occopt,option,smdfun,tphysel,tsmear,tsmearinv,xgrid)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: occopt,option
 real(dp),intent(in) :: tphysel,tsmear
 integer,intent(inout) :: nptsdiv2
 real(dp),intent(out) :: limit,tsmearinv
 real(dp),intent(inout) :: entfun(-nptsdiv2:nptsdiv2,2),occfun(-nptsdiv2:nptsdiv2,2)
 real(dp),intent(inout) :: smdfun(-nptsdiv2:nptsdiv2,2),xgrid(-nptsdiv2:nptsdiv2)


!Local variables-------------------------------
!scalars
 integer :: algo,ii,jj,nconvd2
 integer :: nmaxFD,nminFD
 integer,parameter :: nptsdiv2_def=6000
 integer,save :: dblsmr,occopt_prev=-9999
 real(dp),save :: convlim,incconv,limit_occ,tphysel_prev=-9999,tsmear_prev=-9999
 real(dp) :: aa,dsqrpi,encorr,factor
 real(dp) :: expinc,expx22,expxo2,gauss,increm
 real(dp) :: resFD1,resFD2,resFD3,resFD4,resmom,resmom1,resmom2
 real(dp) :: resmom3,resmom4,secmom,smom1,smom2,thdmom,tmom1,tmom2,tmpexpsum
 real(dp) :: tmpsmdfun,tratio,tt,xx,yp1,ypn
 character(len=500) :: msg
!arrays
 real(dp),save :: entfun_prev(-nptsdiv2_def:nptsdiv2_def,2),occfun_prev(-nptsdiv2_def:nptsdiv2_def,2)
 real(dp),save :: smdfun_prev(-nptsdiv2_def:nptsdiv2_def,2),xgrid_prev(-nptsdiv2_def:nptsdiv2_def)
 real(dp),allocatable :: entder(:),occder(:),smd1(:),smd2(:)
 real(dp),allocatable :: smdder(:),tgrid(:),work(:),workfun(:)

! *************************************************************************

!Initialize the occupation function and generalized entropy function,
!at the beginning, or if occopt changed

 if(option==-1)then
   nptsdiv2 = nptsdiv2_def
   return
 end if

 if (occopt_prev/=occopt .or. abs(tsmear_prev-tsmear)  >tol12 .or. abs(tphysel_prev-tphysel)>tol12) then
   occopt_prev=occopt
   tsmear_prev=tsmear
   tphysel_prev=tphysel

!  Check whether input values of tphysel tsmear and occopt are consistent
   dblsmr = 0
   if (abs(tphysel)>tol12) then
!    Use re-smearing scheme
     if (abs(tsmear)>tol12) then
       dblsmr = 1
!      Use FD occupations (one smearing) only with "physical" temperature tphysel
     else if (occopt /= 3) then
       write(msg, '(a,i6,a)' )' tphysel /= 0, tsmear == 0, but occopt is not = 3, but ',occopt,'.'
       MSG_ERROR(msg)
     end if
   end if

   ABI_ALLOCATE(entder,(-nptsdiv2_def:nptsdiv2_def))
   ABI_ALLOCATE(occder,(-nptsdiv2_def:nptsdiv2_def))
   ABI_ALLOCATE(smdder,(-nptsdiv2_def:nptsdiv2_def))
   ABI_ALLOCATE(workfun,(-nptsdiv2_def:nptsdiv2_def))
   ABI_ALLOCATE(work,(-nptsdiv2_def:nptsdiv2_def))

!  Prepare the points on the grid
!  limit is the value of the argument that will give 0.0 or 1.0 , with
!  less than about 1.0d-15 error for 4<=occopt<=8, and less than about 1.0d-12
!  error for occopt==3. It is not worth to compute the function beyond
!  that point. Even with a less severe requirement, it is significantly
!  larger for occopt==3, with an exponential
!  tail, than for the other occupation functions, with a Gaussian tail.
!  Note that these values are useful in newocc.f also.
   limit_occ=6.0_dp
   if(occopt==3)limit_occ=30.0_dp
   if(dblsmr /= 0) then
     tratio = tsmear / tphysel
     limit_occ=30.0_dp + 6.0_dp*tratio
   end if

!  With nptsdiv2_def=6000 (thus increm=0.001 for 4<=occopt<=8,
!  and increm=0.005 for occopt==3, the O(1/N4) algorithm gives 1.0d-12
!  accuracy on the stored values occfun and entfun. These, together
!  with smdfun and xgrid_prev, need permanently about 0.67 MB, which is affordable.
   increm=limit_occ/nptsdiv2_def
   do ii=-nptsdiv2_def,nptsdiv2_def
     xgrid_prev(ii)=ii*increm
   end do

!  ---------------------------------------------------------
!  Ordinary (unique) smearing function
!  ---------------------------------------------------------
   if (dblsmr == 0) then

!    Compute the unnormalized smeared delta function between -limit_occ and +limit_occ
!    (well, they are actually normalized ...)
     if(occopt==3)then

!      Fermi-Dirac
       do ii=0,nptsdiv2_def
         xx=xgrid_prev(ii)
         smdfun_prev( ii,1)=0.25_dp/(cosh(xx/2.0_dp)**2)
         smdfun_prev(-ii,1)=smdfun_prev(ii,1)
       end do

     else if(occopt==4 .or. occopt==5)then

!      Cold smearing of Marzari, two values of the "a" parameter being possible
!      first value gives minimization of the bump
       if(occopt==4)aa=-.5634
!      second value gives monotonic occupation function
       if(occopt==5)aa=-.8165

       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nptsdiv2_def
         xx=xgrid_prev(ii)
         gauss=dsqrpi*exp(-xx**2)
         smdfun_prev( ii,1)=gauss*(1.5_dp+xx*(-aa*1.5_dp+xx*(-1.0_dp+aa*xx)))
         smdfun_prev(-ii,1)=gauss*(1.5_dp+xx*( aa*1.5_dp+xx*(-1.0_dp-aa*xx)))
       end do

     else if(occopt==6)then

!      First order Hermite-Gaussian of Paxton and Methfessel
       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nptsdiv2_def
         xx=xgrid_prev(ii)
         smdfun_prev( ii,1)=dsqrpi*(1.5_dp-xx**2)*exp(-xx**2)
         smdfun_prev(-ii,1)=smdfun_prev(ii,1)
       end do

     else if(occopt==7)then

!      Gaussian smearing
       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nptsdiv2_def
         xx=xgrid_prev(ii)
         smdfun_prev( ii,1)=dsqrpi*exp(-xx**2)
         smdfun_prev(-ii,1)=smdfun_prev(ii,1)
       end do

     else if(occopt==8)then

!      Constant value of the delta function over the smearing interval, for testing purposes only.
       do ii=0,nptsdiv2_def
         xx=xgrid_prev(ii)
         if(xx>half+tol8)then
           smdfun_prev( ii,1)=zero
         else if(xx<half-tol8)then
           smdfun_prev( ii,1)=one
         else
           smdfun_prev( ii,1)=half
         end if
         smdfun_prev(-ii,1)=smdfun_prev(ii,1)
       end do

     else
       write(msg, '(a,i0,a)' )' Occopt=',occopt,' is not allowed in getnel. '
       MSG_BUG(msg)
     end if

!    ---------------------------------------------------------
!    smear FD delta with occopt delta calculated in smdfun_prev
!    ---------------------------------------------------------
   else if (dblsmr /= 0) then

     nconvd2 = 6000
     convlim = 10.0_dp
     incconv = convlim / nconvd2

!    store smearing functions in smd1 and smd2
     ABI_ALLOCATE(smd1,(-nconvd2:nconvd2))
     ABI_ALLOCATE(smd2,(-nconvd2:nconvd2))
     ABI_ALLOCATE(tgrid,(-nconvd2:nconvd2))

!    FD function in smd1( ii) and second smearing delta in smd2( ii)
!
!    smd1(:) contains delta_FD ( x )
     do ii=0,nconvd2
       tgrid(ii)=ii*incconv
       tgrid(-ii)=-tgrid(ii)
       tt=tgrid(ii)
       smd1( ii)=0.25_dp/(cosh(tt/2.0_dp)**2)
       smd1(-ii)=smd1(ii)
     end do

!    check input values of occopt and fill smd2(:) with appropriate data:
!    smd2(:) contains delta_resmear ( x )
     if(occopt == 3) then
       write(msg, '(a,a)' )&
&       'Occopt=3 is not allowed as a re-smearing.', &
&       'Use a single FD, or re-smear with a different delta type (faster cutoff). '
       MSG_ERROR(msg)
     else if(occopt==4 .or. occopt==5)then
!      Cold smearing of Marzari, two values of the "a" parameter being possible
!      first value gives minimization of the bump
       if(occopt==4)aa=-.5634
!      second value gives monotonic occupation function
       if(occopt==5)aa=-.8165

       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nconvd2
         tt=tgrid(ii)
         gauss=dsqrpi*exp(-tt**2)
         smd2( ii)=gauss*(1.5_dp+tt*(-aa*1.5_dp+tt*(-1.0_dp+aa*tt)))
         smd2(-ii)=gauss*(1.5_dp+tt*( aa*1.5_dp+tt*(-1.0_dp-aa*tt)))
       end do
     else if(occopt==6)then
       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nconvd2
         tt=tgrid(ii)
         smd2( ii)=dsqrpi*(1.5_dp-tt**2)*exp(-tt**2)
         smd2(-ii)=smd2(ii)
       end do
     else if(occopt==7)then
       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nconvd2
         tt=tgrid(ii)
         smd2( ii)=dsqrpi*exp(-tt**2)
         smd2(-ii)=smd2(ii)
       end do
     else if(occopt==8)then
       do ii=0,nconvd2
         tt=tgrid(ii)
         if(tt>half+tol8)then
           smd2( ii)=zero
         else if(tt<half-tol8)then
           smd2( ii)=one
         else
           smd2( ii)=half
         end if
         smd2(-ii)=smd2(ii)
       end do
     else
       write(msg, '(a,i0,a)' )' Occopt= ',occopt,' is not allowed in getnel. '
       MSG_BUG(msg)
     end if

!    Use O(1/N4) algorithm from Num Rec (see below)
!
!    The grid for the convoluted delta is taken (conservatively)
!    to be that for the FD delta ie 6000 pts in [-limit_occ;limit_occ]
!    Smearing functions are given on [-dbllim;dbllim] and the grid must
!    superpose the normal grid on [-limit_occ:limit_occ]
!    The maximal interval for integration of the convolution is
!    [-dbllim+limit_occ+lim(delta2);dbllim-limit_occ-lim(delta2)] =
!    [-dbllim+36;dbllim-36]

!    test the smdFD function for extreme values:
!    do jj=-nptsdiv2_def,-nptsdiv2_def
!    do ii=-nconvd2+4,nconvd2
!    call smdFD(xgrid_prev(jj) - tgrid(ii)*tratio, resFD)
!    write(std_out,*) 'ii jj = ', ii,jj, ' smdFD (', xgrid_prev(jj) - tgrid(ii)*tratio, ') ', resFD
!    end do
!    end do

     expinc = exp(half*incconv*tratio)

!    jj = position of point at which we are calculating smdfun_prev
     do jj=-nptsdiv2_def,nptsdiv2_def
!      Do not care about the 8 boundary points,
!      where the values should be extremely small anyway
       smdfun_prev(jj,1)=0.0_dp
!      only add contribution with delta_FD > 1.0d-100
       nmaxFD = floor  (( maxFDarg+xgrid_prev(jj)) / tratio / incconv )
       nmaxFD = min (nmaxFD, nconvd2)
       nminFD = ceiling((-maxFDarg+xgrid_prev(jj)) / tratio / incconv )
       nminFD = max (nminFD, -nconvd2)

!      Calculate the Fermi-Dirac distrib at point xgrid_prev(jj)-tgrid(ii)*tratio
       expxo2 = exp (-half*(xgrid_prev(jj) - (nminFD)*incconv*tratio))
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD4 = tmpexpsum * tmpexpsum
       expxo2 = expxo2*expinc
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD3 = tmpexpsum * tmpexpsum
       expxo2 = expxo2*expinc
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD2 = tmpexpsum * tmpexpsum
       expxo2 = expxo2*expinc
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD1 = tmpexpsum * tmpexpsum

!      core contribution to the integral with constant weight (48)
       tmpsmdfun = 0.0_dp
       do ii=nminFD+4,nmaxFD-4
         expxo2 = expxo2*expinc
!        tmpexpsum = 1.0_dp / (expxo2 + 1.0_dp / expxo2 )
         expx22 = expxo2*expxo2
         tmpexpsum = expxo2 / (expx22 + 1.0_dp)
         tmpsmdfun = tmpsmdfun + smd2(ii) * tmpexpsum * tmpexpsum
       end do

!      Add on end contributions for show (both functions smd and smdFD are very small
       smdfun_prev(jj,1)=smdfun_prev(jj,1)       +48.0_dp*tmpsmdfun             &
&       + 31.0_dp*smd2(nminFD+3)*resFD1 -11.0_dp*smd2(nminFD+2)*resFD2 &
&       +  5.0_dp*smd2(nminFD+1)*resFD3 -       smd2(nminFD)*resFD4

       expxo2 = expxo2*expinc
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD1 = tmpexpsum * tmpexpsum
       expxo2 = expxo2*expinc
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD2 = tmpexpsum * tmpexpsum
       expxo2 = expxo2*expinc
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD3 = tmpexpsum * tmpexpsum
       expxo2 = expxo2*expinc
       expx22 = expxo2*expxo2
       tmpexpsum = expxo2 / (expx22 + 1.0_dp)
       resFD4 = tmpexpsum * tmpexpsum

!      Contribution above
       smdfun_prev(jj,1)=smdfun_prev(jj,1)                                      &
&       + 31.0_dp*smd2(nmaxFD-3)*resFD1  -11.0_dp*smd2(nmaxFD-2)*resFD2 &
&       +  5.0_dp*smd2(nmaxFD-1)*resFD3  -       smd2(nmaxFD)*resFD4
       smdfun_prev(jj,1)=incconv*smdfun_prev(jj,1)/48.0_dp
     end do

     secmom = 0.0_dp
     thdmom = 0.0_dp
     resmom4 = xgrid_prev(-nptsdiv2_def  )*xgrid_prev(-nptsdiv2_def  )*smdfun_prev(-nptsdiv2_def  ,  1)
     resmom3 = xgrid_prev(-nptsdiv2_def+1)*xgrid_prev(-nptsdiv2_def+1)*smdfun_prev(-nptsdiv2_def+1,  1)
     resmom2 = xgrid_prev(-nptsdiv2_def+2)*xgrid_prev(-nptsdiv2_def+2)*smdfun_prev(-nptsdiv2_def+2,  1)
     resmom1 = xgrid_prev(-nptsdiv2_def+3)*xgrid_prev(-nptsdiv2_def+3)*smdfun_prev(-nptsdiv2_def+3,  1)
     resmom  = xgrid_prev(-nptsdiv2_def+4)*xgrid_prev(-nptsdiv2_def+4)*smdfun_prev(-nptsdiv2_def+4,  1)
     do ii=-nptsdiv2_def+4,nptsdiv2_def-1
       secmom = secmom +                                   &
&       ( 17.0_dp*xgrid_prev(ii)  *xgrid_prev(ii)  *smdfun_prev(ii,  1)   &
&       +42.0_dp*xgrid_prev(ii-1)*xgrid_prev(ii-1)*smdfun_prev(ii-1,1)   &
&       -16.0_dp*xgrid_prev(ii-2)*xgrid_prev(ii-2)*smdfun_prev(ii-2,1)   &
&       + 6.0_dp*xgrid_prev(ii-3)*xgrid_prev(ii-3)*smdfun_prev(ii-3,1)   &
&       -       xgrid_prev(ii-4)*xgrid_prev(ii-4)*smdfun_prev(ii-4,1)  )
       resmom4 = resmom3
       resmom3 = resmom2
       resmom2 = resmom1
       resmom1 = resmom
       resmom  = xgrid_prev(ii+1)  *xgrid_prev(ii+1)  *smdfun_prev(ii+1,  1)
     end do
     secmom=increm * secmom / 48.0_dp
!    thdmom=increm * thdmom / 48.0_dp
!
!    smom1  = second moment of delta in smd1(:)
!    smom2  = second moment of delta in smd2(:)
!
     smom1  = 0.0_dp
     smom2  = 0.0_dp
     tmom1  = 0.0_dp
     tmom2  = 0.0_dp
     do ii=-nconvd2+4,nconvd2
       smom1 = smom1+                                       &
&       ( 17.0_dp*tgrid(ii)  *tgrid(ii)  *smd1(ii)         &
&       +42.0_dp*tgrid(ii-1)*tgrid(ii-1)*smd1(ii-1)       &
&       -16.0_dp*tgrid(ii-2)*tgrid(ii-2)*smd1(ii-2)       &
&       + 6.0_dp*tgrid(ii-3)*tgrid(ii-3)*smd1(ii-3)       &
&       -       tgrid(ii-4)*tgrid(ii-4)*smd1(ii-4)  )
       smom2 = smom2+                                       &
&       ( 17.0_dp*tgrid(ii)  *tgrid(ii)  *smd2(ii  )     &
&       +42.0_dp*tgrid(ii-1)*tgrid(ii-1)*smd2(ii-1)     &
&       -16.0_dp*tgrid(ii-2)*tgrid(ii-2)*smd2(ii-2)     &
&       + 6.0_dp*tgrid(ii-3)*tgrid(ii-3)*smd2(ii-3)     &
&       -       tgrid(ii-4)*tgrid(ii-4)*smd2(ii-4)  )
     end do
     smom1 =incconv * smom1  / 48.0_dp
     smom2 =incconv * smom2  / 48.0_dp
!    tmom1 =incconv * tmom1  / 48.0_dp
!    tmom2 =incconv * tmom2  / 48.0_dp

     encorr =  smom2*tratio*tratio/secmom

     ABI_DEALLOCATE(tgrid)
     ABI_DEALLOCATE(smd1)
     ABI_DEALLOCATE(smd2)

   end if

!  --------------------------------------------------------
!  end of smearing function initialisation, dblsmr case
!  --------------------------------------------------------


!  Now that the smeared delta function has been initialized, compute the
!  occupation function
   occfun_prev(-nptsdiv2_def,1)=zero
   entfun_prev(-nptsdiv2_def,1)=zero

!  Different algorithms are possible, corresponding to the formulas
!  (4.1.11), (4.1.12) and (4.1.14) in Numerical recipes (pp 107 and 108),
!  with respective O(1/N2), O(1/N3), O(1/N4) convergence, where N is the
!  number of points in the interval.
   algo=4

   if(algo==2)then

!    Extended trapezoidal rule (4.1.11), taken in a cumulative way
     do ii=-nptsdiv2_def+1,nptsdiv2_def
       occfun_prev(ii,1)=occfun_prev(ii-1,1)+increm*(smdfun_prev(ii,1)+smdfun_prev(ii-1,1))/2.0_dp
       entfun_prev(ii,1)=entfun_prev(ii-1,1)+increm*&
&       ( -xgrid_prev(ii)*smdfun_prev(ii,1) -xgrid_prev(ii-1)*smdfun_prev(ii-1,1) )/2.0_dp
     end do

   else if(algo==3)then

!    Derived from (4.1.12). Converges as O(1/N3).
!    Do not care about the following points,
!    where the values are extremely small anyway
     occfun_prev(-nptsdiv2_def+1,1)=0.0_dp ;   entfun_prev(-nptsdiv2_def+1,1)=0.0_dp
     do ii=-nptsdiv2_def+2,nptsdiv2_def
       occfun_prev(ii,1)=occfun_prev(ii-1,1)+increm*&
&       ( 5.0_dp*smdfun_prev(ii,1) + 8.0_dp*smdfun_prev(ii-1,1) - smdfun_prev(ii-2,1) )/12.0_dp
       entfun_prev(ii,1)=entfun_prev(ii-1,1)+increm*&
&       ( 5.0_dp*(-xgrid_prev(ii)  )*smdfun_prev(ii,1)  &
&       +8.0_dp*(-xgrid_prev(ii-1))*smdfun_prev(ii-1,1)&
&       -      (-xgrid_prev(ii-2))*smdfun_prev(ii-2,1) )/12.0_dp
     end do

   else if(algo==4)then

!    Derived from (4.1.14)- alternative extended Simpsons rule. Converges as O(1/N4).
!    Do not care about the following points,
!    where the values are extremely small anyway
     occfun_prev(-nptsdiv2_def+1,1)=0.0_dp ;   entfun_prev(-nptsdiv2_def+1,1)=0.0_dp
     occfun_prev(-nptsdiv2_def+2,1)=0.0_dp ;   entfun_prev(-nptsdiv2_def+2,1)=0.0_dp
     occfun_prev(-nptsdiv2_def+3,1)=0.0_dp ;   entfun_prev(-nptsdiv2_def+3,1)=0.0_dp
     do ii=-nptsdiv2_def+4,nptsdiv2_def
       occfun_prev(ii,1)=occfun_prev(ii-1,1)+increm*&
&       ( 17.0_dp*smdfun_prev(ii,1)  &
&       +42.0_dp*smdfun_prev(ii-1,1)&
&       -16.0_dp*smdfun_prev(ii-2,1)&
&       + 6.0_dp*smdfun_prev(ii-3,1)&
&       -       smdfun_prev(ii-4,1) )/48.0_dp
       entfun_prev(ii,1)=entfun_prev(ii-1,1)+increm*&
&       ( 17.0_dp*(-xgrid_prev(ii)  )*smdfun_prev(ii,1)  &
&       +42.0_dp*(-xgrid_prev(ii-1))*smdfun_prev(ii-1,1)&
&       -16.0_dp*(-xgrid_prev(ii-2))*smdfun_prev(ii-2,1)&
&       + 6.0_dp*(-xgrid_prev(ii-3))*smdfun_prev(ii-3,1)&
&       -       (-xgrid_prev(ii-4))*smdfun_prev(ii-4,1) )/48.0_dp
     end do

   end if ! End of choice between different algorithms for integration

!  Normalize the functions (actually not needed for occopt=3..7)
   factor=1.0_dp/occfun_prev(nptsdiv2_def,1)
   smdfun_prev(:,1)=smdfun_prev(:,1)*factor
   occfun_prev(:,1)=occfun_prev(:,1)*factor
   entfun_prev(:,1)=entfun_prev(:,1)*factor

!  Compute the cubic spline fitting of the smeared delta function
   yp1=0.0_dp ; ypn=0.0_dp
   workfun(:)=smdfun_prev(:,1)
   call spline(xgrid_prev, workfun, (2*nptsdiv2_def+1), yp1, ypn, smdder)
   smdfun_prev(:,2)=smdder(:)

!  Compute the cubic spline fitting of the occupation function
   yp1=0.0_dp ; ypn=0.0_dp
   workfun(:)=occfun_prev(:,1)
   call spline(xgrid_prev, workfun, (2*nptsdiv2_def+1), yp1, ypn, occder)
   occfun_prev(:,2)=occder(:)

!  Compute the cubic spline fitting of the entropy function
   yp1=0.0_dp ; ypn=0.0_dp
   workfun(:)=entfun_prev(:,1)
   call spline(xgrid_prev, workfun, (2*nptsdiv2_def+1), yp1, ypn, entder)
   entfun_prev(:,2)=entder(:)

   ABI_DEALLOCATE(entder)
   ABI_DEALLOCATE(occder)
   ABI_DEALLOCATE(smdder)
   ABI_DEALLOCATE(work)
   ABI_DEALLOCATE(workfun)

 end if

 if (abs(tphysel)<tol12) then
   if (tsmear == zero) then
     tsmearinv = huge_tsmearinv
   else
     tsmearinv=one/tsmear
   end if
 else
   tsmearinv=one/tphysel
 end if

 entfun(:,:) = entfun_prev(:,:)
 occfun(:,:) = occfun_prev(:,:)
 smdfun(:,:) = smdfun_prev(:,:)
 xgrid(:) = xgrid_prev(:)
 limit = limit_occ
 nptsdiv2 = nptsdiv2_def

end subroutine init_occ_ent
!!***

!!****f* m_occ/occeig
!! NAME
!! occeig
!!
!! FUNCTION
!! For each pair of active bands (m,n), generates ratios
!! that depend on the difference between occupation numbers
!! and eigenvalues.
!!
!! INPUTS
!!  doccde_k(nband_k)=derivative of occ_k wrt the energy
!!  doccde_kq(nband_k)=derivative of occ_kq wrt the energy
!!  eig0_k(nband_k)=GS eigenvalues at k
!!  eig0_kq(nband_k)=GS eigenvalues at k+q
!!  nband_k=number of bands
!!  occopt=option for occupancies
!!  occ_k(nband_k)=occupation number for each band at k
!!  occ_kq(nband_k)=occupation number for each band at k+q
!!
!! OUTPUT
!!  rocceig(nband_k,nband_k)$= (occ_{k,q}(m)-occ_k(n))/(eig0_{k,q}(m)-eig0_k(n))$,
!!   if this ratio has been attributed to the band n, 0.0_dp otherwise
!!
!! NOTES
!! Supposing the occupations numbers differ:
!! if $abs(occ_{k,q}(m)) < abs(occ_k(n))$
!!  $rocceig(m,n)=(occ_{k,q}(m)-occ_k(n))/(eig0_{k,q}(m)-eig0_k(n)) $
!! if $abs(occ_{k,q}(m))>abs(occ_k(n))$
!!  rocceig(m,n)=0.0_dp
!!
!! If the occupation numbers are close enough, then
!! if the eigenvalues are also close, take the derivative
!!  $ rocceig(m,n)=\frac{1}{2}*docc/deig0 $
!! otherwise,
!!  $ rocceig(m,n)=\frac{1}{2}*(occ_{k,q}(m)-occ_k(n))/(eig0_{k,q}(m)-eig0_k(n))$
!!
!! PARENTS
!!      dfpt_nstpaw,dfpt_rhofermi,dfpt_vtorho
!!
!! CHILDREN
!!
!! SOURCE

subroutine occeig(doccde_k,doccde_kq,eig0_k,eig0_kq,nband_k,occopt,occ_k,occ_kq,rocceig)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nband_k,occopt
!arrays
 real(dp),intent(in) :: doccde_k(nband_k),doccde_kq(nband_k),eig0_k(nband_k)
 real(dp),intent(in) :: eig0_kq(nband_k),occ_k(nband_k),occ_kq(nband_k)
 real(dp),intent(out) :: rocceig(nband_k,nband_k)

!Local variables-------------------------------
!scalars
 integer :: ibandk,ibandkq
 real(dp) :: diffabsocc,diffeig,diffocc,ratio,sumabsocc
 character(len=500) :: msg

! *************************************************************************

!The parameter tol5 defines the treshhold for degeneracy, and the width of the step function

 rocceig(:,:)=0.0_dp

 do ibandk=1,nband_k
   do ibandkq=1,nband_k
     diffeig=eig0_kq(ibandkq)-eig0_k(ibandk)
     diffocc=occ_kq(ibandkq)-occ_k(ibandk)

     if( abs(diffeig) > tol5 ) then
       ratio=diffocc/diffeig
     else
       if(occopt<3)then
!        In a non-metallic case, if the eigenvalues are degenerate,
!        the occupation numbers must also be degenerate, in which
!        case there is no contribution from this pair of bands
         if( abs(diffocc) > tol5 ) then
           write(msg,'(a,a,a,a,a,a,a,2(a,i4,a,es16.6,a,es16.6,a,a),a)' ) &
&           'In a non-metallic case (occopt<3), for a RF calculation,',ch10,&
&           'if the eigenvalues are degenerate,',' the occupation numbers must also be degenerate.',ch10,&
&           'However, the following pair of states gave :',ch10,&
&           'k -state, band number',ibandk,', occ=',occ_k(ibandk),'eigenvalue=',eig0_k(ibandk),',',ch10,&
&           ' kq-state, band number',ibandkq,', occ=',occ_kq(ibandkq),', eigenvalue=',eig0_kq(ibandkq),'.',ch10,&
&           'Action: change occopt, consistently, in GS and RF calculations.'
           MSG_ERROR(msg)
         end if
         ratio=0.0_dp
       else
!        In the metallic case, one can compute a better approximation of the
!        ratio by using derivatives doccde
         ratio=0.5_dp*(doccde_kq(ibandkq)+doccde_k(ibandk))
!        write(std_out,*)' occeig : ibandkq,doccde_kq(ibandkq)',ibandkq,doccde_kq(ibandkq)
!        write(std_out,*)'          ibandk ,doccde_k (ibandk )',ibandk,doccde_k(ibandk)
       end if
     end if

!    Here, must pay attention to the smallness of some coefficient
     diffabsocc=abs(occ_k(ibandk))-abs(occ_kq(ibandkq))
     sumabsocc=abs(occ_k(ibandk))+abs(occ_kq(ibandkq))
     if(sumabsocc>tol8)then
       if( diffabsocc > sumabsocc*tol5 ) then
         rocceig(ibandkq,ibandk)=ratio
       else if ( diffabsocc >= -sumabsocc*tol5 ) then
         rocceig(ibandkq,ibandk)=0.5_dp*ratio
       else
         rocceig(ibandkq,ibandk)=0.0_dp
       end if
     end if

   end do ! ibandkq
 end do ! ibandk

end subroutine occeig
!!***

!----------------------------------------------------------------------

!!****f* m_occ/occ_fd
!! NAME
!!  occ_fd
!!
!! FUNCTION
!!  Fermi-Dirac statistic: 1 / [(exp((e - mu)/ KT) + 1]
!!  Note that occ_fs in [0, 1] so the spin factor is not included, unlike the
!!  occupations stored in ebands%occ.
!!
!! INPUTS
!!   ee=Single particle energy in Ha
!!   kT=Value of K_Boltzmann x T in Ha.
!!   mu=Chemical potential in Ha.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental real(dp) function occ_fd(ee, kT, mu)

!Arguments ------------------------------------
 real(dp),intent(in) :: ee, kT, mu

!Local variables ------------------------------
 real(dp) :: ee_mu,arg
! *************************************************************************

 ee_mu = ee - mu

 !TODO: Find good tols.
 ! 1 kelvin [K] = 3.16680853419133E-06 Hartree
 if (kT > tol6) then
   arg = ee_mu / kT
   if (arg > maxFDarg) then
     occ_fd = zero
   else if (arg < -maxFDarg) then
     occ_fd = one
   else
     occ_fd = one / (exp(arg) + one)
   end if
 else
   ! Heaviside
   if (ee_mu > zero) then
     occ_fd = zero
   else if (ee_mu < zero) then
     occ_fd = one
   else
     occ_fd = half
   end if
 end if

end function occ_fd
!!***

!----------------------------------------------------------------------

!!****f* m_occ/occ_dfd
!! NAME
!!  occ_dfd
!!
!! FUNCTION
!!  Derivative of Fermi-Dirac statistic: (exp((e - mu)/ KT) / KT[(exp((e - mu)/ KT) + 1]^2
!!  Note that kT is given in Hartree so the derivative as well
!!
!! INPUTS
!!   ee=Single particle energy in Ha
!!   kT=Value of K_Boltzmann x T in Ha.
!!   mu=Chemical potential in Ha.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental real(dp) function occ_dfd(ee, kT, mu)

!Arguments ------------------------------------
 real(dp),intent(in) :: ee, kT, mu

!Local variables ------------------------------
 real(dp) :: ee_mu,arg
! *************************************************************************

 ee_mu = ee - mu

 !TODO: Find good tols.
 ! 1 kelvin [K] = 3.16680853419133E-06 Hartree
 if (kT > tol6) then
   arg = ee_mu / kT
   if (arg > maxDFDarg) then
     occ_dfd = zero
   else if (arg < -maxDFDarg) then
     occ_dfd = zero
   else
     occ_dfd = exp(arg) / (exp(arg) + one)**2 / kT
   end if
 else
   occ_dfd = zero
 end if

end function occ_dfd
!!***

!----------------------------------------------------------------------

!!****f* m_occ/occ_be
!! NAME
!!  occ_be
!!
!! FUNCTION
!!   Bose-Einstein statistic  1 / [(exp((e - mu)/ KT) - 1]
!!
!! INPUTS
!!   ee=Single particle energy in Ha
!!   kT=Value of K_Boltzmann x T in Ha.
!!   mu=Chemical potential in Ha (usually zero)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental real(dp) function occ_be(ee, kT, mu)

!Arguments ------------------------------------
 real(dp),intent(in) :: ee, kT, mu

!Local variables ------------------------------
 real(dp) :: ee_mu, arg
! *************************************************************************

 ee_mu = ee - mu

 !TODO: Find good tols.
 ! 1 kelvin [K] = 3.16680853419133E-06 Hartree
 if (kT > tol12) then
   arg = ee_mu / kT
   if (arg > tol12 .and. arg < maxBEarg) then
     occ_be = one / (exp(arg) - one)
   else
     occ_be = zero
   end if
 else
   ! No condensate for T --> 0
   occ_be = zero
 end if

end function occ_be
!!***

!----------------------------------------------------------------------

!!****f* m_occ/occ_dbe
!! NAME
!!  occ_dbe
!!
!! FUNCTION
!!   Derivative of Bose-Einstein statistic  (exp((e - mu)/ KT) / KT[(exp((e - mu)/ KT) - 1]^2
!!   Note that kT is given in Hartree so the derivative as well
!!
!! INPUTS
!!   ee=Single particle energy in Ha
!!   kT=Value of K_Boltzmann x T in Ha.
!!   mu=Chemical potential in Ha (usually zero)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental real(dp) function occ_dbe(ee, kT, mu)

!Arguments ------------------------------------
 real(dp),intent(in) :: ee, kT, mu

!Local variables ------------------------------
 real(dp) :: ee_mu, arg
! *************************************************************************

 ee_mu = ee - mu

 !TODO: Find good tols.
 ! 1 kelvin [K] = 3.16680853419133E-06 Hartree
 if (kT > tol12) then
   arg = ee_mu / kT
   if (arg > tol12 .and. arg < maxDBEarg) then
     occ_dbe = exp(arg) / (kT * (exp(arg) - one)**2)
   else
     occ_dbe = zero
   end if
 else
   ! No condensate for T --> 0
   occ_dbe = zero
 end if

end function occ_dbe
!!***

!----------------------------------------------------------------------


!!****f* m_occ/dos_hdr_write
!!
!! NAME
!! dos_hdr_write
!!
!! FUNCTION
!! Write the header of the DOS files, for both smearing and tetrahedron methods.
!!
!! INPUTS
!! deltaene=increment of DOS energy arguments
!! enemax=maximal value of the DOS energy argument
!! enemin=minimal value of the DOS energy argument
!! nene=number of DOS energy argument
!! eigen(mband*nkpt*nsppol)=eigenvalues (input or init to large number), hartree
!! fermie=fermi energy useful for band alignment...
!! mband=maximum number of bands
!! nband(nkpt*nsppol)=number of bands at each k point
!! nkpt=number of k points
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! occopt=option for occupancies, or re-smearing scheme if dblsmr /= 0
!! prtdos=1 for smearing technique, 2 or 3 for tetrahedron technique
!! tphysel="physical" electronic temperature with FD occupations
!! tsmear=smearing width (or temperature)
!! unitdos=unit number of output of the DOS.
!!
!! OUTPUT
!!   Only writing.
!!
!! PARENTS
!!      getnel,m_epjdos
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine dos_hdr_write(deltaene,eigen,enemax,enemin,fermie,mband,nband,nene,&
&  nkpt,nsppol,occopt,prtdos,tphysel,tsmear,unitdos)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkpt,nsppol,occopt,prtdos,unitdos,nene
 real(dp),intent(in) :: fermie,tphysel,tsmear
 real(dp),intent(in) :: deltaene,enemax,enemin
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *************************************************************************

!Write the DOS file
 write(msg, '(7a,i2,a,i5,a,i4)' ) "#",ch10, &
& '# ABINIT package : DOS file  ',ch10,"#",ch10,&
& '# nsppol =',nsppol,', nkpt =',nkpt,', nband(1)=',nband(1)
 call wrtout(unitdos,msg,'COLL')

 if (any(prtdos== [1,4])) then
   write(msg, '(a,i2,a,f6.3,a,f6.3,a)' )  &
&   '# Smearing technique, occopt =',occopt,', tsmear=',tsmear,' Hartree, tphysel=',tphysel,' Hartree'
 else
   write(msg, '(a)' ) '# Tetrahedron method '
 end if
 call wrtout(unitdos,msg,'COLL')

 if (mband*nkpt*nsppol>=3) then
   write(msg, '(a,3f8.3,2a)' )'# For identification : eigen(1:3)=',eigen(1:3),ch10,"#"
 else
   write(msg, '(a,3f8.3)' ) '# For identification : eigen=',eigen
   write(msg, '(3a)')trim(msg),ch10,"#"
 end if
 call wrtout(unitdos,msg,'COLL')

 write(msg, '(a,f16.8)' ) '# Fermi energy : ', fermie
 call wrtout(unitdos,msg,'COLL')

 if (prtdos==1) then
   write(msg, '(5a)' ) "#",ch10,&
&   '# The DOS (in electrons/Hartree/cell) and integrated DOS (in electrons/cell),',&
&   ch10,'# as well as the DOS with tsmear halved and doubled, are computed,'

 else if (prtdos==2)then
   write(msg, '(3a)' ) "#",ch10,&
&   '# The DOS (in electrons/Hartree/cell) and integrated DOS (in electrons/cell) are computed,'

 else if (any(prtdos == [3, 4])) then
   write(msg, '(5a)' ) "#",ch10,&
&   '# The local DOS (in electrons/Hartree for one atomic sphere)',ch10,&
&   '# and integrated local DOS (in electrons for one atomic sphere) are computed.'

 else if (prtdos==5)then
   write(msg, '(9a)' ) "#",ch10,&
&   '# The spin component DOS (in electrons/Hartree/cell)',ch10,&
&   '# and integrated spin component DOS (in electrons/cell) are computed.',ch10,&
&   '# Remember that the wf are eigenstates of S_z and S^2, not S_x and S_y',ch10,&
&   '#   so the latter will not always sum to 0 for paired electronic states.'
 end if
 call wrtout(unitdos,msg,'COLL')

 write(msg, '(a,i5,a,a,a,f9.4,a,f9.4,a,f8.5,a,a,a)' )&
& '# at ',nene,' energies (in Hartree) covering the interval ',ch10,&
& '# between ',enemin,' and ',enemax,' Hartree by steps of ',deltaene,' Hartree.',ch10,"#"
 call wrtout(unitdos,msg,'COLL')

 if (prtdos==1) then
   write(msg, '(a,a)' )&
&   '#       energy        DOS       Integr. DOS   ','     DOS           DOS    '
   call wrtout(unitdos,msg,'COLL')

   write(msg, '(a)' )&
&   '#                                              (tsmear/2)    (tsmear*2) '
   call wrtout(unitdos,msg,'COLL')
 else
   write(msg, '(a)' ) '#       energy        DOS '
 end if

end subroutine dos_hdr_write
!!***

!!****f* m_occ/pareigocc
!! NAME
!! pareigocc
!!
!! FUNCTION
!! This subroutine transmit to all processors, using MPI:
!!   - the eigenvalues and,
!!   - if ground-state, the occupation numbers
!!     (In fact, in the present status of the routine,
!!     occupation numbers are NOT transmitted)
!!     transmit_occ = 2 is used in case the occ should be transmitted.
!!     Yet the code is not already written.
!!
!! INPUTS
!!  formeig=format of eigenvalues (0 for GS, 1 for RF)
!!  localrdwf=(for parallel case) if 1, the eig and occ initial values
!!            are local to each machine, if 0, they are on proc me=0.
!!  mband=maximum number of bands of the output wavefunctions
!!  mpi_enreg=information about MPI parallelization
!!  nband(nkpt*nsppol)=desired number of bands at each k point
!!  nkpt=number of k points
!!  nsppol=1 for unpolarized, 2 for spin-polarized, output wf file processors,
!!         Warning : defined only when paralbd=1
!!  transmit_occ/=2 transmit only eigenvalues, =2 for transmission of occ also
!!         (yet transmit_occ=2 is not safe or finished at all)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  eigen(mband*nkpt*nsppol)=eigenvalues (input or init to large number), (Ha)
!!  occ(mband*nkpt*nsppol)=occupation (input or init to 0.0)  NOT USED NOW
!!
!! NOTES
!! * The case paralbd=1 with formeig=0 is implemented, but not yet used.
!!
!! * The transmission of occ is not activated yet !
!!
!! * The routine takes the eigenvalues in the eigen array on one of the
!!   processors that possess the wavefunctions, and transmit it to all procs.
!!   If localrdwf==0, me=0 has the full array at start,
!!   If localrdwf==1, the transfer might be more complex.
!!
!! * This routine should not be used for RF wavefunctions, since
!!   it does not treat the eigenvalues as a matrix.
!!
!! PARENTS
!!      newkpt,wfsinp
!!
!! CHILDREN
!!      timab,xmpi_bcast,xmpi_sum
!!
!! SOURCE

subroutine pareigocc(eigen,formeig,localrdwf,mpi_enreg,mband,nband,nkpt,nsppol,occ,transmit_occ)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: formeig,localrdwf,mband,nkpt,nsppol,transmit_occ
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(inout) :: eigen(mband*(2*mband)**formeig*nkpt*nsppol)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer :: band_index,iband,ierr,ikpt,isppol,me,nbks,spaceComm
 !character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: buffer1(:),buffer2(:)

! *************************************************************************

 if(xmpi_paral==1)then

!  Init mpi_comm
   spaceComm=mpi_enreg%comm_cell
   if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_kpt
   if(mpi_enreg%paral_hf==1) spaceComm=mpi_enreg%comm_kpt
!  Init me
   me=mpi_enreg%me_kpt

   if(localrdwf==0)then
     call xmpi_bcast(eigen,0,spaceComm,ierr)

   else if(localrdwf==1)then

!    Prepare transmission of eigen (and occ)
     ABI_ALLOCATE(buffer1,(2*mband**(formeig+1)*nkpt*nsppol))
     ABI_ALLOCATE(buffer2,(2*mband**(formeig+1)*nkpt*nsppol))
     buffer1(:)=zero
     buffer2(:)=zero

     band_index=0
     do isppol=1,nsppol
       do ikpt=1,nkpt
         nbks=nband(ikpt+(isppol-1)*nkpt)

         if(mpi_enreg%paralbd==0)then

           if(formeig==0)then
             buffer1(2*band_index+1:2*band_index+nbks) = eigen(band_index+1:band_index+nbks)
             if(transmit_occ==2) then
               buffer1(2*band_index+nbks+1:2*band_index+2*nbks) = occ(band_index+1:band_index+nbks)
             end if
             band_index=band_index+nbks
           else if(formeig==1)then
             buffer1(band_index+1:band_index+2*nbks**2) = eigen(band_index+1:band_index+2*nbks**2)
             band_index=band_index+2*nbks**2
           end if

         else if(mpi_enreg%paralbd==1)then

!          Skip this k-point if not the proper processor
           if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nbks,isppol,me)) then
             if(formeig==0) then
               band_index=band_index+nbks
             else
               band_index=band_index+2*nbks**2
             end if
             cycle
           end if
!          Loop on bands
           do iband=1,nbks
             if(mpi_enreg%proc_distrb(ikpt, iband,isppol) /= me)cycle
             if(formeig==0)then
               buffer1(2*band_index+iband)=eigen(band_index+iband)
!              if(transmit_occ==2) buffer1(2*band_index+iband+nbdks)=occ(band_index+iband)
             else if (formeig==1)then
               buffer1(band_index+(iband-1)*2*nbks+1:band_index+(iband-1)*2*nbks+2*nbks) = &
&               eigen(band_index+(iband-1)*2*nbks+1:band_index+(iband-1)*2*nbks+2*nbks)
             end if
           end do
           if(formeig==0)then
             band_index=band_index+nbks
           else
             band_index=band_index+2*nbks**2
           end if
         end if

       end do
     end do

!    Build sum of everything
     call timab(48,1,tsec)
     if(formeig==0)band_index=band_index*2
     call xmpi_sum(buffer1,buffer2,band_index,spaceComm,ierr)
     call timab(48,2,tsec)

     band_index=0
     do isppol=1,nsppol
       do ikpt=1,nkpt
         nbks=nband(ikpt+(isppol-1)*nkpt)
         if(formeig==0)then
           eigen(band_index+1:band_index+nbks) = buffer2(2*band_index+1:2*band_index+nbks)
           if(transmit_occ==2) then
             occ(band_index+1:band_index+nbks) = buffer2(2*band_index+nbks+1:2*band_index+2*nbks)
           end if
           band_index=band_index+nbks
         else if(formeig==1)then
           eigen(band_index+1:band_index+2*nbks**2) = buffer1(band_index+1:band_index+2*nbks**2)
           band_index=band_index+2*nbks**2
         end if
       end do
     end do

     ABI_DEALLOCATE(buffer1)
     ABI_DEALLOCATE(buffer2)
   end if
 end if

end subroutine pareigocc
!!***

end module m_occ
!!***
