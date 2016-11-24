!{\src2tex{textfont=tt}}
!!****f* ABINIT/getnel
!! NAME
!! getnel
!!
!! FUNCTION
!! Option=1 :
!! Get the total number of electrons nelect, given a trial fermienergy fermie.
!! For this, compute new occupation numbers at each k point,
!! from eigenenergies eigen, according to the
!! smearing scheme defined by occopt (and smearing width tsmear or tphysel).
!!
!! Option=2 :
!! Compute and output the smeared density of states, and the integrated density
!! of states, then write these data
!!
!! Warning : this routine assumes checks have been done in the calling
!! routine, and that the values of the arguments are sensible
!!
!! NOTE
!! in order to speed the calculation, it would be easy to
!! compute the entropy only when the fermi energy is well converged
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG, AF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!       $  E_{phys} = E_{free} - encorr*(E_{int}-E_{free}) + O(tseamr^3)  $
!!
!! PARENTS
!!      cchi0q0_intraband,clnup1,conducti_nc,dfpt_looppert,m_ebands,newocc
!!
!! CHILDREN
!!      dos_hdr_write,init_occ_ent,splfit,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine getnel(doccde,dosdeltae,eigen,entropy,fermie,maxocc,mband,nband,&
&  nelect,nkpt,nsppol,occ,occopt,option,tphysel,tsmear,unitdos,wtk)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_splines

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getnel'
 use interfaces_14_hidewrite
 use interfaces_61_occeig, except_this_one => getnel
!End of the abilint section

 implicit none

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
 integer,parameter :: nptsdiv2_def=6000
 integer,parameter :: prtdos1=1
 integer :: bantot,iband,iene,ikpt,index,index_start,isppol
 integer :: nene,nptsdiv2
 real(dp) :: buffer,deltaene,dosdbletot,doshalftot,dostot
 real(dp) :: enemax,enemin,enex,intdostot,limit,tsmearinv
 character(len=500) :: message
!arrays
 real(dp),allocatable :: entfun(:,:),occfun(:,:)
 real(dp),allocatable :: smdfun(:,:),xgrid(:)
 real(dp),allocatable :: arg(:),derfun(:),dos(:),dosdble(:),doshalf(:),ent(:)
 real(dp),allocatable :: intdos(:)

! *************************************************************************

 DBG_ENTER("COLL")

 if(option/=1 .and. option/=2)then
   write(message,'(a,i0,a)')' Option must be either 1 or 2. It is ',option,'.'
   MSG_BUG(message)
 end if

!Initialize the occupation function and generalized entropy function,
!at the beginning, or if occopt changed

!Just get the number nptsdiv2 and allocate entfun, occfun,
!smdfun and xgrid accordingly
 nptsdiv2 = nptsdiv2_def

! call init_occ_ent(entfun, limit, &
!& nptsdiv2, occfun, occopt, -1, smdfun, tphysel, &
!& tsmear, tsmearinv, xgrid)

 ABI_ALLOCATE(entfun,(-nptsdiv2:nptsdiv2,2))
 ABI_ALLOCATE(occfun,(-nptsdiv2:nptsdiv2,2))
 ABI_ALLOCATE(smdfun,(-nptsdiv2:nptsdiv2,2))
 ABI_ALLOCATE(xgrid,(-nptsdiv2:nptsdiv2))

!Call to init_occ_ent
 call init_occ_ent(entfun, limit, &
& nptsdiv2, occfun, occopt, option, smdfun, tphysel, &
& tsmear, tsmearinv, xgrid)

!The initialisation of occfun and entfun is done

!---------------------------------------------------------------------

!write(std_out,*)' getnel : debug  tphysel, tsmear = ', tphysel, tsmear
 bantot=sum(nband(:))

 ABI_ALLOCATE(arg,(bantot))
 ABI_ALLOCATE(derfun,(bantot))
 ABI_ALLOCATE(ent,(bantot))

 if(option==1)then
   !normal evaluation of occupations and entropy

!  Compute the arguments of the occupation and entropy functions
   arg(:)=(fermie-eigen(1:bantot))*tsmearinv

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
   call dos_hdr_write(buffer,deltaene,dosdeltae,eigen,enemax,enemin,fermie,mband,nband,nene,&
&   nkpt,nsppol,occopt,prtdos1,tphysel,tsmear,unitdos)

   ABI_ALLOCATE(dos,(bantot))
   ABI_ALLOCATE(dosdble,(bantot))
   ABI_ALLOCATE(doshalf,(bantot))
   ABI_ALLOCATE(intdos,(bantot))

   do isppol=1,nsppol

     if (nsppol==2) then
       if(isppol==1) write(message,'(a,16x,a)')  '#','Spin-up DOS'
       if(isppol==2) write(message,'(2a,16x,a)')  ch10,'#','Spin-dn DOS '
       call wrtout(unitdos,message,'COLL')
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
