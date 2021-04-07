! CP modified to include occopt 9 option
!!****m* ABINIT/m_occ
!! NAME
!! m_occ
!!
!! FUNCTION
!!  Low-level functions for occupation factors.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2021 ABINIT group (XG, AF)
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

 use m_time,         only : timab, cwtime, cwtime_report
 use m_fstrings,     only : sjoin, itoa

 implicit none

 private
!!***

 public :: getnel        ! Compute total number of electrons from efermi or DOS
 public :: newocc        ! Compute new occupation numbers at each k point,
 public :: occeig        ! (occ_{k,q}(m)-occ_k(n))/(eig0_{k,q}(m)-eig0_k(n))$,
 public :: occ_fd        ! Fermi-Dirac statistics 1 / [(exp((e - mu)/ KT) + 1]
 public :: occ_dfde      ! Derivative of Fermi-Dirac statistics wrt e: (exp((e - mu)/ KT) / KT[(exp((e - mu)/ KT) + 1]^2
 public :: occ_be        ! Bose-Einstein statistics  1 / [(exp((e - mu)/ KT) - 1]
 public :: occ_dbe       ! Derivative of Bose-Einstein statistics  (exp((e - mu)/ KT) / KT[(exp((e - mu)/ KT) - 1]^2
 public :: dos_hdr_write


 integer,parameter :: nptsdiv2_def=6000
 ! This parameter is used in init_occ_ent and getnel
 ! nptsdiv2 is the number of integration points, divided by 2.

 real(dp),parameter :: huge_tsmearinv = 1e50_dp
 real(dp),parameter :: maxFDarg = 500.0_dp
 real(dp),parameter :: maxDFDarg = 200.0_dp
 real(dp),parameter :: maxBEarg = 600.0_dp
 real(dp),parameter :: maxDBEarg = 200.0_dp


contains
!!***

!!****f* m_abinit/getnel
!! NAME
!! getnel
!!
!! FUNCTION
!! Option = 1:
!!   Get the total number of electrons nelect, given a trial fermienergy fermie.
!!   For this, compute new occupation numbers at each k point,
!!   from eigenenergies eigen, according to the
!!   smearing scheme defined by occopt (and smearing width tsmear or tphysel).
!!
!! Option = 2:
!!   Compute and output the smeared density of states, and the integrated density
!!   of states, then write these data
!!
!! Warning: this routine assumes checks have been done in the calling
!! routine, and that the values of the arguments are sensible
!!
!! NOTE
!! In order to speed the calculation, it would be easy to
!! compute the entropy only when the fermi energy is well converged
!!
!! INPUTS
!! dosdeltae= DOS delta of Energy (needed if Option=2)
!! eigen(mband*nkpt*nsppol)=eigenvalues (input or init to large number), hartree
!! fermie= fermi energy/ fermi energy for excited electrons if occopt = 9 (Hartree) ! CP description modified
!! fermih= fermi energy for excited holes (Hartree) ! CP added
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
!! iB1, iB2 = band min and max between which to calculate the number of electrons
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
!!
!!       $  E_{phys} = E_{free} - encorr*(E_{int}-E_{free}) + O(tsmear^3)  $
!!
!! PARENTS
!!      m_chi0,m_conducti,m_dfpt_looppert,m_ebands,m_gstate,m_occ
!!
!! CHILDREN
!!      timab,xmpi_bcast,xmpi_sum
!!
!! SOURCE

subroutine getnel(doccde, dosdeltae, eigen, entropy, fermie, fermih, maxocc, mband, nband, &
                  nelect, nkpt, nsppol, occ, occopt, option, tphysel, tsmear, unitdos, wtk, iB1, iB2)
          ! CP added fermih, iB1 and iB2 in the list of arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkpt,nsppol,occopt,option,unitdos
 real(dp),intent(in) :: dosdeltae,fermie,fermih,maxocc,tphysel,tsmear ! CP added fermih
 real(dp),intent(out) :: entropy,nelect
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),wtk(nkpt)
 real(dp),intent(out) :: doccde(mband*nkpt*nsppol)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
 integer, intent(in), optional:: iB1, iB2 !! CP: added optional arguments to get number of electrons between bands iB1 and iB2
 !! Used only when occopt = 9

!Local variables-------------------------------
! nptsdiv2 is the number of integration points, divided by 2.
! tratio  = ratio tsmear/tphysel for convoluted smearing function
! save values so we can impose recalculation of smdfun when
! the smearing or electronic temperature change between datasets
! corresponds roughly to delta_FD (maxFDarg) = 1.0d-100
!
! return fermi-dirac smearing function analytically
! real(dp) :: smdFD
! smdFD (tt) = 1.0_dp / (exp(-tt/2.0_dp) + exp(tt/2.0_dp))**2
!scalars
 integer,parameter :: prtdos1=1
 integer :: iband,iene,ikpt,index,index_tot,index_start,isppol, nene,nptsdiv2 ! CP added index_tot, removed bantot
 real(dp) :: buffer,deltaene,dosdbletot,doshalftot,dostot, wk
 real(dp) :: enemax,enemin,enex,intdostot,limit,tsmearinv
 !real(dp) :: cpu, wall, gflops
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: entfun(:,:),occfun(:,:)
 real(dp),allocatable :: smdfun(:,:),xgrid(:)
 real(dp),allocatable :: arg(:),derfun(:),dos(:),dosdble(:),doshalf(:),ent(:)
 real(dp),allocatable :: intdos(:)
 ! CP added: variables
 real(dp),allocatable :: occ_tmp(:),ent_tmp(:),doccde_tmp(:)
 integer:: low_band_index, high_band_index, number_of_bands
 ! End CP added: variables

! *************************************************************************

 !call cwtime(cpu, wall, gflops, "start")

 if (option/=1 .and. option/=2)then
   ABI_BUG(sjoin('Option must be either 1 or 2. It is:', itoa(option)))
 end if

 ! Initialize the occupation function and generalized entropy function,
 ! at the beginning, or if occopt changed

 ! CP added
 if (occopt==9) then
    low_band_index  = iB1
    high_band_index = iB2
    number_of_bands = (iB2-iB1+1)*nkpt*nsppol
 else
    low_band_index  = 1
    high_band_index = nband(1)
    number_of_bands = sum(nband(:))
 end if
 ABI_MALLOC(occ_tmp,(number_of_bands))
 ABI_MALLOC(ent_tmp,(number_of_bands))
 ABI_MALLOC(doccde_tmp,(number_of_bands))
 ! End CP addition

 ! Just get the number nptsdiv2 and allocate entfun, occfun, smdfun and xgrid accordingly
 nptsdiv2 = nptsdiv2_def

 ABI_MALLOC(entfun,(-nptsdiv2:nptsdiv2,2))
 ABI_MALLOC(occfun,(-nptsdiv2:nptsdiv2,2))
 ABI_MALLOC(smdfun,(-nptsdiv2:nptsdiv2,2))
 ABI_MALLOC(xgrid,(-nptsdiv2:nptsdiv2))

 call init_occ_ent(entfun, limit, nptsdiv2, occfun, occopt, option, smdfun, tphysel, tsmear, tsmearinv, xgrid)
 ! The initialisation of occfun and entfun is done

!---------------------------------------------------------------------

 ! write(std_out,*)' getnel : debug  tphysel, tsmear = ', tphysel, tsmear
 ! CP delete
 !bantot=sum(nband(:))
 ! End CP delete
 ! CP modified
 !ABI_MALLOC(arg,(bantot))
 !ABI_MALLOC(derfun,(bantot))
 !ABI_MALLOC(ent,(bantot))
 ABI_MALLOC(arg,(number_of_bands))
 ABI_MALLOC(derfun,(number_of_bands))
 ABI_MALLOC(ent,(number_of_bands))
! Enc CP modified
 if (option==1) then
   ! normal evaluation of occupations and entropy

   ! CP modified for occopt = 9 option, to account for subset of bands bw iB1 and iB2
   !! Compute the arguments of the occupation and entropy functions
   !! HM 20/08/2018 Treat the T --> 0 limit
   !if (tsmear==0) then
   !  arg(:)=sign(huge_tsmearinv,fermie-eigen(1:bantot))
   !else
   !  arg(:)=(fermie-eigen(1:bantot))*tsmearinv
   !endif
   index = 0
   index_tot = 0
   do isppol=1,nsppol
      do ikpt=1,nkpt
         if (occopt == 2) high_band_index=nband(ikpt+nkpt*(isppol-1))
         do iband=low_band_index,high_band_index
            index = index + 1
            if (tsmear==0) then
               arg(index) = sign(huge_tsmearinv,fermie-eigen(index_tot + iband))
            else
               arg(index)=(fermie-eigen(index_tot + iband))*tsmearinv
            end if
         end do
         index_tot = index_tot + nband(ikpt+nkpt*(isppol-1))
      end do
   end do
   ! End CP modified

   ! MG TODO: This part is expensive for dense k-meshes
   ! Compute the values of the occupation function, and the entropy function
   ! Note: splfit also takes care of the points outside of the interval,
   ! and assign to them the value of the closest extremal point,
   ! which is what is needed here.

   ! CP modified
   !call splfit(xgrid, doccde, occfun, 1, arg, occ, (2*nptsdiv2+1), bantot)
   !call splfit(xgrid, derfun, entfun, 0, arg, ent, (2*nptsdiv2+1), bantot)
   call splfit(xgrid, doccde_tmp, occfun, 1, arg, occ_tmp, (2*nptsdiv2+1), number_of_bands)
   call splfit(xgrid, derfun, entfun, 0, arg, ent, (2*nptsdiv2+1), number_of_bands)
   ! End CP modified

   ! Normalize occ and ent, and sum number of electrons and entropy
   ! Use different loops for nelect and entropy because bantot may be quite large in the EPH code
   ! when we use very dense k-meshes.

   ! CP modified
   !ent = ent * maxocc
   !occ = occ * maxocc
   !doccde = -doccde * maxocc * tsmearinv

   nelect=zero; entropy=zero
   index=0
   !do isppol=1,nsppol
   !  do ikpt=1,nkpt
   !    wk = wtk(ikpt)
   !    do iband=1,nband(ikpt+nkpt*(isppol-1))
   !      index = index  +1
   !      entropy = entropy + wk * ent(index)
   !    end do
   !  end do
   !end do

   !index=0
   !do isppol=1,nsppol
   !  do ikpt=1,nkpt
   !    wk = wtk(ikpt)
   !    do iband=1,nband(ikpt+nkpt*(isppol-1))
   !      index = index  +1
   !      nelect = nelect + wk * occ(index)
   !    end do
   !  end do
   !end do
   index_tot = 0
   do isppol=1,nsppol
      do ikpt=1,nkpt
         wk = wtk(ikpt)
         if (occopt == 2) high_band_index=nband(ikpt+nkpt*(isppol-1))
         do iband=low_band_index,high_band_index
            index = index + 1
            ent(index)                = ent(index)*maxocc
            occ(iband + index_tot)    = occ_tmp(index)*maxocc
            doccde(iband + index_tot) = -doccde_tmp(index)*maxocc*tsmearinv
            entropy                   = entropy + wk*ent(index)
            nelect                    = nelect + wk*occ(iband + index_tot)
         end do
         index_tot = index_tot + nband(ikpt+nkpt*(isppol-1))
      end do
   end do
   ! End CP modified

   !write(std_out,*) ' getnel : debug   wtk, occ, eigen = ', wtk, occ, eigen
   !write(std_out,*)xgrid(-nptsdiv2),xgrid(nptsdiv2)
   !write(std_out,*)'fermie',fermie
   !do ii=1,bantot
   !write(std_out,*)ii,arg(ii),doccde(ii)
   !end do
   !write(std_out,*)'eigen',eigen(:)
   !write(std_out,*)'arg',arg(:)
   !write(std_out,*)'occ',occ(:)
   !write(std_out,*)'nelect',nelect

 else if (option==2) then
   ! evaluate DOS for smearing, half smearing, and double.

   buffer=limit/tsmearinv*.5_dp

   ! A Similar section is present is dos_calcnwrite. Should move all DOS stuff to m_ebands
   ! Choose the lower and upper energies
   ! CP modified
   !enemax=maxval(eigen(1:bantot))+buffer
   !enemin=minval(eigen(1:bantot))-buffer
   enemax=maxval(eigen(1:number_of_bands))+buffer
   enemin=minval(eigen(1:number_of_bands))-buffer
   ! End CP modified

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

   ! Write the header of the DOS file, and also decides the energy range and increment
   ! CP modified to fit new argument list of dos_hdr_write
   !call dos_hdr_write(deltaene,eigen,enemax,enemin,fermie,mband,nband,nene,&
   !  nkpt,nsppol,occopt,prtdos1,tphysel,tsmear,unitdos)
   call dos_hdr_write(deltaene,eigen,enemax,enemin,fermie,fermih,mband,nband,nene,&
           nkpt,nsppol,occopt,prtdos1,tphysel,tsmear,unitdos)
   !ABI_MALLOC(dos,(bantot))
   !ABI_MALLOC(dosdble,(bantot))
   !ABI_MALLOC(doshalf,(bantot))
   !ABI_MALLOC(intdos,(bantot))
   ABI_MALLOC(dos,(number_of_bands))
   ABI_MALLOC(dosdble,(number_of_bands))
   ABI_MALLOC(doshalf,(number_of_bands))
   ABI_MALLOC(intdos,(number_of_bands))

   do isppol=1,nsppol

     if (nsppol==2) then
       if(isppol==1) write(msg,'(a,16x,a)')  '#','Spin-up DOS'
       if(isppol==2) write(msg,'(2a,16x,a)')  ch10,'#','Spin-dn DOS '
       call wrtout(unitdos,msg)
     end if
     index_start=0
     if(isppol==2)then
       do ikpt=1,nkpt
   !      index_start=index_start+nband(ikpt)
          if (occopt == 2) high_band_index = nband(ikpt + nkpt*(isppol - 1))
          index_start=index_start + high_band_index - low_band_index + 1
       end do
     end if
   ! End CP modified

     enex=enemin
     do iene=1,nene

       ! Compute the arguments of the dos and occupation function
       ! CP modified
       ! arg(:)=(enex-eigen(1:bantot))*tsmearinv
       arg(:)=(enex-eigen(1:number_of_bands))*tsmearinv
       ! End CP modified

       !call splfit(xgrid,derfun,smdfun,0,arg,dos,(2*nptsdiv2+1),bantot)
       !call splfit(xgrid,derfun,occfun,0,arg,intdos,(2*nptsdiv2+1),bantot)
       call splfit(xgrid,derfun,smdfun,0,arg,dos,(2*nptsdiv2+1),number_of_bands)
       call splfit(xgrid,derfun,occfun,0,arg,intdos,(2*nptsdiv2+1),number_of_bands)

       ! Also compute the dos with tsmear halved and doubled
       arg(:)=arg(:)*2.0_dp
       !call splfit(xgrid,derfun,smdfun,0,arg,doshalf,(2*nptsdiv2+1),bantot)
       call splfit(xgrid,derfun,smdfun,0,arg,doshalf,(2*nptsdiv2+1),number_of_bands)

       ! Since arg was already doubled, must divide by four
       arg(:)=arg(:)*0.25_dp
       !call splfit(xgrid,derfun,smdfun,0,arg,dosdble,(2*nptsdiv2+1),bantot)
       call splfit(xgrid,derfun,smdfun,0,arg,dosdble,(2*nptsdiv2+1),number_of_bands)

       ! Now, accumulate the contribution from each eigenenergy
       dostot=zero
       intdostot=zero
       doshalftot=zero
       dosdbletot=zero
       index=index_start

     ! CP modified
       ! write(std_out,*)' eigen, arg, dos, intdos, doshalf, dosdble'
       !do ikpt=1,nkpt
       !  do iband=1,nband(ikpt+nkpt*(isppol-1))
       !    index=index+1
       !    dostot=dostot+wtk(ikpt)*maxocc*dos(index)*tsmearinv
       !    intdostot=intdostot+wtk(ikpt)*maxocc*intdos(index)
       !    doshalftot=doshalftot+wtk(ikpt)*maxocc*doshalf(index)*tsmearinv*2.0_dp
       !    dosdbletot=dosdbletot+wtk(ikpt)*maxocc*dosdble(index)*tsmearinv*0.5_dp
       !  end do
       !end do

       do ikpt=1,nkpt
          if (occopt == 2) high_band_index=nband(ikpt+nkpt*(isppol-1))
          do iband=low_band_index,high_band_index
             index=index+1
             dostot=dostot+wtk(ikpt)*maxocc*dos(index)*tsmearinv
             intdostot=intdostot+wtk(ikpt)*maxocc*intdos(index)
             doshalftot=doshalftot+wtk(ikpt)*maxocc*doshalf(index)*tsmearinv*2.0_dp
             dosdbletot=dosdbletot+wtk(ikpt)*maxocc*dosdble(index)*tsmearinv*0.5_dp
          end do
       end do
   ! End CP modified

       ! Print the data for this energy
       write(unitdos, '(f8.3,2f14.6,2f14.3)' )enex,dostot,intdostot,doshalftot,dosdbletot

       enex=enex+deltaene
     end do ! iene
   end do ! isppol

   ABI_FREE(dos)
   ABI_FREE(dosdble)
   ABI_FREE(doshalf)
   ABI_FREE(intdos)

   ! MG: It does not make sense to close the unit here since the routines
   ! did not open the file here!
   ! Close the DOS file
   close(unitdos)
 end if

 ABI_FREE(arg)
 ABI_FREE(derfun)
 ABI_FREE(ent)
 ABI_FREE(entfun)
 ABI_FREE(occfun)
 ABI_FREE(smdfun)
 ABI_FREE(xgrid)
! CP added
 ABI_FREE(occ_tmp)
 ABI_FREE(doccde_tmp)
 ABI_FREE(ent_tmp)
! End CP added

 !call cwtime_report(" getnel", cpu, wall, gflops, end_str=ch10)

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
!!  ne_qFD, nh_qFD=number of thermalized excited electrons (resp. holes) in bands > ivalence (resp. <= ivalence) ! CP added for
!occopt 9 case
!!  ivalence= band index of the last valence band ! CP added for occopt 9 case
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
!!  fermie= fermi energy (Hartree)/fermi level for thermalized excited electrons in bands > ivalence when occopt=9
!!  fermih= fermi level for thermalized excited holes in bands <= ivalence ! CP added for occopt 9 case
!!  occ(maxval(nband(:))*nkpt*nsppol)=occupancies for each band and k point
!!
!! PARENTS
!!      m_ebands,m_gstate,m_respfn_driver,m_vtorho
!!
!! CHILDREN
!!      timab,xmpi_bcast,xmpi_sum
!!
!! SOURCE

subroutine newocc(doccde, eigen, entropy, fermie, fermih, ivalence, spinmagntarget, mband, nband, &
  nelect, ne_qFD, nh_qFD, nkpt, nspinor, nsppol, occ, occopt, prtvol, stmbias, tphysel, tsmear, wtk) ! CP modified:
!  added fermih, ivalence, ne_qFD, nh_qFD for occopt 9 case

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkpt,nspinor,nsppol,occopt,prtvol
 real(dp),intent(in) :: spinmagntarget,nelect,stmbias,tphysel,tsmear
 real(dp),intent(out) :: entropy,fermie,fermih ! CP added fermih
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),wtk(nkpt)
 real(dp),intent(out) :: doccde(mband*nkpt*nsppol)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
 !real(dp),optional,intent(in)
 ! CP added
 real(dp),intent(in) :: ne_qFD, nh_qFD
 integer, intent(in) :: ivalence
 ! End CP added

!Local variables-------------------------------
 integer,parameter :: niter_max = 120, nkpt_max = 2, fake_unit=-666, option1 = 1
 integer :: cnt,cnt2,cnt3,ib,ii,ik,ikpt,is,isppol,nkpt_eff,  sign
 integer,allocatable :: nbandt(:)
 real(dp),parameter :: tol = tol14
 !real(dp),parameter :: tol = tol10
 real(dp) :: dosdeltae,entropy_tmp,fermie_hi,fermie_lo,fermie_mid,fermie_mid_tmp ! CP modified
 real(dp) :: fermih_lo,fermih_mid,fermih_hi ! CP added
 real(dp) :: fermie_biased,maxocc
 real(dp) :: nelect_tmp,nelecthi,nelectlo,nelectmid,nelect_biased
 real(dp) :: nholeshi,nholeslo,nholesmid ! CP added
 real(dp) :: entropyet(2),fermie_hit(2),fermie_lot(2),fermie_midt(2),nelecthit(2) ! CP modified
 real(dp) :: nelectlot(2),nelectt(2),tsec(2)
 real(dp) :: entropye, entropyh ! CP added
 real(dp),allocatable :: doccdet(:),eigent(:),occt(:)
 character(len=500) :: msg
 ! CP added
 logical::not_enough_bands=.false.
 ! End CP added

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(74,1,tsec)

 ! Here treat the case where occopt does not correspond to a metallic occupation scheme
 if (occopt < 3 .or. occopt > 9) then ! CP modified
   ABI_BUG(sjoin(' occopt= ',itoa(occopt),', a value not allowed in newocc.'))
 end if

 ! Check whether nband is a constant for all k point and spin-pol
 do isppol=1,nsppol
   do ikpt=1,nkpt
     if(nband(ikpt+(isppol-1)*nkpt)/=nband(1))then
       write(msg,'(3a,i0,a,i0,a,i0,a)')&
        'The number of bands must be the same for all k-points ',ch10,&
        'but nband(1)= ',nband(1),' is different of nband(',ikpt+(isppol-1)*nkpt,') = ',nband(ikpt+(isppol-1)*nkpt),'.'
       ABI_BUG(msg)
     end if
   end do
 end do

 ! Check whether nelect is strictly positive
 if (nelect <= zero) then
   write(msg,'(3a,es16.8,a)')&
   'nelect must be a positive number, while ',ch10, 'the calling routine asks nelect= ',nelect,'.'
   ABI_BUG(msg)
 end if

 ! CP added: Check whether the number of holes and electrons if positive
 if (occopt == 9) then
    if ( (ne_qFD < zero) .or. (nh_qFD < zero) ) then
       write(msg,'(3a,es16.8,a,es16.8,a)')&
&   'ne_qFD or nh_qFD must be positive numbers, while ',ch10,&
&   'the calling routine asks ne_qFD= ',ne_qFD,' and nh_qFD= ',nh_qFD, '.'
   ABI_BUG(msg)
    end if
 end if
 ! End CP added

 maxocc = two / (nsppol * nspinor)

 ! Check whether nelect is coherent with nband (nband(1) is enough,
 ! since it was checked that nband is independent of k-point and spin-pol
 if (nelect > nband(1) * nsppol * maxocc) then
   write(msg,'(3a,es16.8,a,i0,a,es16.8,a)' )&
   'nelect must be smaller than nband*maxocc, while ',ch10,&
   'the calling routine gives nelect= ',nelect,', nband= ',nband(1),' and maxocc= ',maxocc,'.'
   ABI_BUG(msg)
 end if

! CP added: Providing additional checks to ensure that there are enough valence and conduction bands to accomodate ne_qFD and nh_qFD
 if( occopt==9 .and. ne_qFD > (nband(1)-ivalence)*nsppol*maxocc )then
   write(msg,'(a,es16.8,2a,es16.8,a)') 'ne_qFD = ', ne_qFD ,ch10, &
&   'must be smaller than (nband-ivalence)*maxocc*nsppol = ', &
&   (nband(1)-ivalence)*nsppol*maxocc,'.'
   ABI_BUG(msg)
 else if( occopt==9 .and. (nh_qFD > ivalence*nsppol*maxocc .or. &
&                          nelect - nh_qFD > ivalence*nsppol*maxocc ) )then
   write(msg,'(a,es16.8,2a,es16.8,2a,es16.8,a)') 'nh_qFD = ', nh_qFD ,ch10, &
&   'and nelect-nh_qFD = ', nelect - nh_qFD,ch10, ' must be smaller than ivalence*maxocc*nsppol = ', &
&   ivalence*nsppol*maxocc,'.'
   ABI_BUG(msg)
  end if
! End CP added

 ! Use bisection algorithm to find fermi energy
 ! This choice is due to the fact that it will always give sensible
 ! result (because the answer is bounded, even if the smearing function
 ! is non-monotonic (which is the case for occopt=4 or 6)
 ! Might speed up it, if needed !

 ! Lowest and largest trial fermi energies, and corresponding number of electrons
 ! They are obtained from the smallest or largest eigenenergy, plus a range of
 ! energy that allows for complete occupation of all bands, or, on the opposite,
 ! for zero occupation of all bands (see getnel.f)

 dosdeltae = zero  ! the DOS is not computed, with option=1
 fermie_lo = minval(eigen(1:nband(1)*nkpt*nsppol)) - 6.001_dp * tsmear ! CP modified fermi_lo ->fermie_lo
 if (occopt == 3) fermie_lo = fermie_lo - 24.0_dp * tsmear ! CP modified
 if(occopt==9)fermih_lo = fermie_lo ! CP added to take into account holes
 !if (present(ef_range) fermilo = ef_range(1)

 ! CP modified
 !call getnel(doccde,dosdeltae,eigen,entropy,fermilo,maxocc,mband,nband,&
 ! nelectlo,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk)
 if(occopt >= 3 .and. occopt <= 8) then
    call getnel(doccde,dosdeltae,eigen,entropye,fermie_lo,fermie_lo,maxocc,mband,nband,&
& nelectlo,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk,1,nband(1))
 else if (occopt == 9) then
    call getnel(doccde,dosdeltae,eigen,entropye,fermie_lo,fermie_lo,maxocc,mband,nband,&
& nelectlo,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk, ivalence+1, nband(1)) ! Excited electrons
    call getnel(doccde,dosdeltae,eigen,entropyh,fermih_lo,fermih_lo,maxocc,mband,nband,&
& nholeslo,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk,1, ivalence)
 end if
 !
 !fermihi = maxval(eigen(1:nband(1)*nkpt*nsppol)) + 6.001_dp * tsmear
 fermie_hi = maxval(eigen(1:nband(1)*nkpt*nsppol)) + 6.001_dp * tsmear
 !! Safety value
 !fermihi = min(fermihi, 1.e6_dp)
 fermie_hi = min(fermie_hi, 1.e6_dp)
 !if(occopt == 3) fermihi = fermihi + 24.0_dp * tsmear
 if(occopt == 3 .or. occopt == 9) fermie_hi = fermie_hi + 24.0_dp * tsmear
 if(occopt == 9) fermih_hi=fermie_hi
 !!if (present(ef_range) fermihi = ef_range(2)
 !
 !call getnel(doccde,dosdeltae,eigen,entropy,fermihi,maxocc,mband,nband,&
 ! nelecthi,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk)
 if (occopt >= 3 .and. occopt <= 8) then
    call getnel(doccde,dosdeltae,eigen,entropye,fermie_hi,fermie_hi,maxocc,mband,nband,&
& nelecthi,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk,1,nband(1))
 else if (occopt == 9) then
    call getnel(doccde,dosdeltae,eigen,entropye,fermie_hi,fermie_hi,maxocc,mband,nband,&
& nelecthi,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk, ivalence+1, nband(1)) ! Excited electrons
    fermih_hi=fermie_hi
    call getnel(doccde,dosdeltae,eigen,entropyh,fermih_hi,fermih_hi,maxocc,mband,nband,&
& nholeshi,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk,1, ivalence)
 end if
 !
 !!write(std_out,'(2(a, es16.8))' )' newocc: initial nelect_lo: ',nelectlo, " nelect_hi: ", nelecthi

 ! End CP modified

 ! Prepare fixed moment calculation
 if(abs(spinmagntarget+99.99_dp)>1.0d-10)then
 ! CP added
   if (occopt==9)then
      write(msg,'(a)') 'occopt=9 and spinmagntarget not implemented.'
      ABI_ERROR(msg)
   end if
 ! End CP added
   sign = 1
   do is = 1, nsppol
     fermie_hit(is) = fermie_hi ! CP modified (name only)
     fermie_lot(is) = fermie_lo ! CP modified (name only)
     nelectt(is) = half*(nelect+sign*spinmagntarget)
     sign = -sign
     nelecthit(is) = nelecthi
     nelectlot(is) = nelectlo
   end do
 end if

 ! If the target nelect is not between nelectlo and nelecthi, exit
 if (nelect < nelectlo .or. nelect > nelecthi) then
   not_enough_bands = .true.
   write(msg, '(a,a,a,a,d16.8,a,a,d16.8,a,d16.8,a,a,d16.8,a,d16.8)') ch10,&
    ' newocc: ',ch10,&
    '  The calling routine gives nelect= ',nelect,ch10,&
    '  The lowest bound is ',fermie_lo,', with nelect=',nelectlo,ch10,&
    '  The highest bound is ',fermie_hi,', with nelect=',nelecthi
   call wrtout(std_out, msg)
   ABI_BUG(msg)
 end if

  ! CP added special test for occopt == 9
 if( occopt==9 ) then
    if ((nelect-nh_qFD)<nholeslo .or. (nelect-nh_qFD)>nholeshi) then
       not_enough_bands = .true.
       write(msg,'(a,a,a,d16.8,a,a,d16.8,a,d16.8,a)') 'newocc : ',ch10, &
&      'The calling routine gives nelect-nh_qFD = ', nelect-nh_qFD, ch10, &
&       'The lowest (highest resp.) bound for nelect-nh_qFD is ', &
&   nholeslo, ' ( ', nholeshi, ' ).'
       ABI_BUG(msg)
    endif
    if ((ne_qFD < nelectlo) .or. (ne_qFD > nelecthi) ) then
       not_enough_bands = .true.
       write(msg,'(a,a,a,d16.8,a,a,d16.8,a,d16.8,a)') 'newocc : ',ch10, &
&   'The calling routine gives ne_qFD = ', ne_qFD, ch10, 'The lowest (highest resp.) bound for ne_qFD are ',&
&   nelectlo, ' ( ', nelecthi, ' ) .'
       ABI_BUG(msg)
    endif

   if (not_enough_bands) then
      write(msg, '(11a)' )&
&      'In order to get the right number of carriers,',ch10,&
&      'it seems that the Fermi energies must be outside the range',ch10,&
&      'of eigenenergies, plus 6 or 30 times the smearing, which is strange.',ch10,&
&      'It might be that your number of bands (nband) corresponds to the strictly',ch10,&
&      'minimum number of bands to accomodate your electrons (so, OK for an insulator),',ch10,&
&      'while you are trying to describe a metal. In this case, increase nband, otherwise ...'
      ABI_BUG(msg)
   end if
 end if
 ! End CP added

 if( abs(spinmagntarget+99.99_dp) < tol10) then

   ! Usual bisection loop
   do ii=1,niter_max
     fermie_mid = (fermie_hi + fermie_lo) * half ! CP modified (name only)
     ! CP modified
     if (occopt == 9) fermih_mid=(fermih_hi+fermih_lo)*half
     ! Produce nelectmid from fermimid
     !call getnel(doccde,dosdeltae,eigen,entropy,fermimid,maxocc,mband,nband,&
     !  nelectmid,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk)
     if (occopt /= 9) then

       call getnel(doccde,dosdeltae,eigen,entropye,fermie_mid,fermie_mid,maxocc,mband,nband,&
&     nelectmid,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk, 1, nband(1))
      !write(std_out,'(a,i0,1x, 3(a,es13.5))' ) " iter: ", ii, &
      !  ' fermi_mid: ',fermimid * Ha_eV, ', n_mid: ',nelectmid, &
      !  ", (n_mid-nelect)/nelect: ", (nelectmid - nelect) / nelect

      !if (nelectmid > nelect * (one - tol)) then
      !  fermihi = fermimid
      !  nelecthi = nelectmid
      !end if
      !if (nelectmid < nelect * (one + tol)) then
      !  fermilo = fermimid
      !  nelectlo = nelectmid
      !end if
       if(nelectmid>nelect*(one-tol14))then
         fermie_hi=fermie_mid
         nelecthi=nelectmid
       end if
       if(nelectmid<nelect*(one+tol14))then
         fermie_lo=fermie_mid
         nelectlo=nelectmid
       end if

     else

       call getnel(doccde,dosdeltae,eigen,entropye,fermie_mid,fermie_mid,maxocc,mband,nband,&
&     nelectmid,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk, ivalence+1, nband(1))
       call getnel(doccde,dosdeltae,eigen,entropyh,fermih_mid,fermih_mid,maxocc,mband,nband,&
&     nholesmid,nkpt,nsppol,occ,occopt,option1,tphysel,tsmear,fake_unit,wtk,1,ivalence)
       if(nelectmid>ne_qFD*(one-tol14))then
         fermie_hi = fermie_mid
         nelecthi  = nelectmid
       else if (nelectmid<ne_qFD*(one-tol14))then
         fermie_lo = fermie_mid
         nelectlo  = nelectmid
       end if
       if(nholesmid>(nelect-nh_qFD)*(one-tol14))then
         fermih_hi = fermih_mid
         nholeshi  = nholesmid
       else if(nholesmid<(nelect-nh_qFD)*(one+tol14))then
         fermih_lo = fermih_mid
         nholeslo  = nholesmid
       end if

     end if

     !if (abs(nelectmid - nelect) <= nelect*two*tol) exit
     !write(std_out,'(2(a,es13.5))' )' bisection move: fermi_lo: ',fermilo * Ha_eV,", fermi_hi: ", fermihi * Ha_eV

!     if (abs(nelecthi - nelectlo) <= nelect*two*tol .or. &
!         abs(fermihi - fermilo) <= tol * abs(fermihi + fermilo) ) exit
     if (occopt /= 9) then
        if( abs(nelecthi-nelectlo) <= nelect*two*tol14 .or. abs(fermie_hi-fermie_lo) <= tol14*abs(fermie_hi+fermie_lo) ) exit
     else
        if( ( abs(nelecthi-nelectlo) <= ne_qFD*two*tol14 .or. &
&             abs(fermie_hi-fermie_lo) <= tol14*abs(fermie_hi+fermie_lo) ) .and. &
            ( abs(nholeshi-nholeslo) <= (nelect-nh_qFD)*two*tol14 .or. &
&             abs(fermih_hi-fermih_lo) <= tol14*abs(fermih_hi+fermih_lo) ) ) exit
     end if

     if (ii == niter_max) then
       write(msg,'(a,i0,3a,es22.14,a,es22.14,a)')&
        'It was not possible to find Fermi energy in ',niter_max,' max bisections.',ch10,&
        'nelecthi: ',nelecthi,', and nelectlo: ',nelectlo,'.'
       ABI_BUG(msg)
       if (occopt == 9) then
          write(msg,'(a,es22.14,a,es22.14,a)')&
          'nholesi = ',nholeshi,', and holeslo = ',nholeslo,'.'
       end if
     end if
   ! End CP modified
   end do ! End of bisection loop

   fermie = fermie_mid ! CP modified (name only)
   entropy= entropye ! CP added (necessary because now we have entropy for the electron and hole subsystems)

   ! CP modified
   !write(msg, '(2(a,f14.6),3a,f14.6,a,i0)' ) &
   !' newocc: new Fermi energy is ',fermie," (Ha),", fermie * Ha_eV, " (eV)", ch10,&
   !'         with nelect: ',nelectmid,', after number of bisections: ',ii
   if (occopt /= 9) then
      write(msg, '(2(a,f14.6),a,i0)' ) &
&   ' newocc: new Fermi energy is ',fermie,' , with nelect=',nelectmid,', Number of bisection calls: ',ii
   else
      fermih=fermih_mid
      entropy = entropy + entropyh ! CP: adding entropy of the holes subsystem
      write(msg, '(2(a,f14.6),a,i0)' ) &
&   ' newocc: new Fermi energy for excited electrons is ',fermie,' , with ne_qFD=',nelectmid,', Number of bisection calls: ',ii
      call wrtout(std_out,msg,'COLL')
      write(msg, '(2(a,f14.6),a,i0)' ) &
&   ' newocc: new Fermi energy for excited holes     is ',fermih,' , with nh_qFD=',nelect-nholesmid,&
&   ', Number of bisection calls: ',ii
   end if
   ! End CP modified
   call wrtout(std_out,msg)

   !  Compute occupation numbers for prtstm/=0, close to the Fermi energy
   if (abs(stmbias) > tol10) then
      ! CP added to prevent use with occopt = 9 so far
      if (occopt == 9) then
         write(msg,'(a)') 'Occopt 9 and prtstm /=0 not implemented together. Change occopt or prtstm.'
         ABI_ERROR(msg)
      end if
      ! End CP added
     fermie_biased = fermie - stmbias ! CP modify name
     ABI_MALLOC(occt,(mband*nkpt*nsppol))

     ! CP modify
     !call getnel(doccde,dosdeltae,eigen,entropy,fermi_biased,maxocc,mband,nband,&
     !  nelect_biased,nkpt,nsppol,occt,occopt,option1,tphysel,tsmear,fake_unit,wtk)
     call getnel(doccde,dosdeltae,eigen,entropy,fermie_biased,fermie_biased,maxocc,mband,nband,&
&       nelect_biased,nkpt,nsppol,occt,occopt,option1,tphysel,tsmear,fake_unit,wtk,1,nband(1))
     ! End CP modify
     occ(:)=occ(:)-occt(:)
     nelect_biased = abs(nelectmid - nelect_biased)
     ! Here, arrange to have globally positive occupation numbers, irrespective of the stmbias sign
     if (-stmbias > tol10) occ(:) = -occ(:)
     ABI_FREE(occt)

     write(msg,'(a,f14.6)')' newocc: the number of electrons in the STM range is nelect_biased=',nelect_biased
     call wrtout(std_out,msg)
   end if


 else
   ! Calculations with a specified moment
   ! Bisection loop
   cnt2=0
   cnt3=0
   entropy=zero
   maxocc=one
   ABI_MALLOC(doccdet,(nkpt*mband))
   ABI_MALLOC(eigent,(nkpt*mband))
   ABI_MALLOC(occt,(nkpt*mband))
   ABI_MALLOC(nbandt,(nkpt))

   do is = 1, nsppol
     nelect_tmp = nelectt(is)
     fermie_hi = fermie_hit(is) ! CP modify name
     fermie_lo = fermie_lot(is) ! CP modify name
     nelecthi = nelecthit(is)
     nelectlo = nelectlot(is)
     ! write(std_out,'(a,i1,3(f8.4,1x))') "Spin, N(spin):", is, nelect, fermihi, fermilo
     ! write(std_out,'(a,2(f8.4,1x))') "Hi, lo:", nelecthi, nelectlo

     do ii=1,niter_max
       fermie_mid_tmp=(fermie_hi+fermie_lo)/2.0_dp ! CP modify name
       ! temporary arrays
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

       ! Produce nelectmid from fermimid
       ! CP modify
       ! call getnel(doccdet,dosdeltae,eigent,entropy_tmp,fermimid_tmp,maxocc,mband,nbandt,&
       !  nelectmid,nkpt,1,occt,occopt,option1,tphysel,tsmear,fake_unit,wtk)
       call getnel(doccdet,dosdeltae,eigent,entropy_tmp,fermie_mid_tmp,fermie_mid_tmp,maxocc,mband,nbandt,&
         nelectmid,nkpt,1,occt,occopt,option1,tphysel,tsmear,fake_unit,wtk,1,nband(1))

       !entropyt(is) = entropy_tmp
       !fermimidt(is) = fermimid_tmp
       !fermimid = fermimidt(is)
       entropyet(is) = entropy_tmp
       fermie_midt(is) = fermie_mid_tmp
       fermie_mid = fermie_midt(is)
       ! End CP modify

       ! temporary arrays
       cnt = 0
       do ik = 1, nkpt
         do ib = 1, mband
           cnt = cnt + 1
           occ(cnt+cnt2) = occt(cnt)
           doccde(cnt+cnt2) = doccdet(cnt)
         end do
       end do
       ! write(std_out,'(a,es24.16,a,es24.16)' )' newocc: from fermi=',fermimid,', getnel gives nelect=',nelectmid

       if(nelectmid>=nelect_tmp)then
         fermie_hi=fermie_mid_tmp ! CP modify name
         nelecthi=nelectmid
       else
         fermie_lo=fermie_mid_tmp ! CP modify name
         nelectlo=nelectmid
       end if
       if( abs(nelecthi-nelectlo) <= 1.0d-13 .or. abs(fermie_hi-fermie_lo) <= 0.5d-14*abs(fermie_hi+fermie_lo) ) exit

       if(ii==niter_max)then
         write(msg,'(a,i3,3a,es22.14,a,es22.14,a)')&
          'It was not possible to find Fermi energy in ',niter_max,' bisections.',ch10,&
          'nelecthi: ',nelecthi,', and nelectlo: ',nelectlo,'.'
         ABI_BUG(msg)
       end if
     end do ! End of bisection loop

     cnt2 = cnt2 + nkpt*mband
     entropy = entropy + entropyet(is)
     fermie=fermie_mid ! CP modify name
     write(msg, '(a,i2,a,f14.6,a,f14.6,a,a,i4)' ) &
       ' newocc: new Fermi energy for spin ', is, ' is ',fermie,' , with nelect: ',nelectmid,ch10,&
       '  Number of bisection calls =',ii
     call wrtout(std_out,msg)

   end do ! spin

   ABI_FREE(doccdet)
   ABI_FREE(eigent)
   ABI_FREE(nbandt)
   ABI_FREE(occt)

 end if ! End of logical on fixed moment calculations

 !write(std_out,*) "kT*Entropy:", entropy*tsmear

 ! MG: If you are wondering why this part is npw disabled by default consider that this output
 ! is produced many times in the SCF cycle and in EPH we have to call this routine for
 ! several temperature and the log becomes unreadable.
 ! If you really need to look at the occupation factors use prtvol > 0.
 nkpt_eff = nkpt
 if (prtvol == 0) nkpt_eff = 0
 if (prtvol == 1) nkpt_eff = min(nkpt_max, nkpt)

 if (nsppol == 1)then

   if (nkpt_eff /= 0) then
     write(msg, '(a,i0,a)' )' newocc: computed new occ. numbers for occopt= ',occopt,' , spin-unpolarized case. '
     call wrtout(std_out,msg)
     do ikpt=1,nkpt_eff
       write(msg,'(a,i4,a)' ) ' k-point number ',ikpt,' :'
       do ii=0,(nband(1)-1)/12
         if (ii == 3 .and. prtvol /= 0) exit
         write(msg,'(12f6.3)') occ(1+ii*12+(ikpt-1)*nband(1):min(12+ii*12,nband(1))+(ikpt-1)*nband(1))
         call wrtout(std_out,msg)
       end do
     end do
     if (nkpt /= nkpt_eff) call wrtout(std_out,' newocc: prtvol=0, stop printing more k-point information')

     !call wrtout(std_out,' newocc: corresponding derivatives are ')
     !do ikpt=1,nkpt_eff
     !write(msg,'(a,i4,a)' ) ' k-point number ',ikpt,' :'
     !do ii=0,(nband(1)-1)/12
     !write(msg,'(12f6.1)') doccde(1+ii*12+(ikpt-1)*nband(1):min(12+ii*12,nband(1))+(ikpt-1)*nband(1))
     !call wrtout(std_out,msg)
     !end do
     !end do
     !if(nkpt/=nkpt_eff)then
     !  call wrtout(std_out,'newocc: prtvol=0, stop printing more k-point information')
     !end if
   end if

 else

   if (nkpt_eff /= 0) then
     write(msg, '(a,i0,2a)' )' newocc: computed new occupation numbers for occopt= ',occopt,ch10,'  (1) spin up   values  '
     call wrtout(std_out, msg)
     do ikpt=1,nkpt_eff
       write(msg,'(a,i0,a)' ) ' k-point number ',ikpt,':'
       do ii=0,(nband(1)-1)/12
         if (ii == 3 .and. prtvol /= 0) exit
         write(msg,'(12f6.3)') occ(1+ii*12+(ikpt-1)*nband(1):min(12+ii*12,nband(1))+(ikpt-1)*nband(1))
         call wrtout(std_out,msg)
       end do
     end do
     if (nkpt/=nkpt_eff) call wrtout(std_out,' newocc: prtvol=0, stop printing more k-point information')

     call wrtout(std_out,'  (2) spin down values  ')
     do ikpt=1,nkpt_eff
       do ii=0,(nband(1)-1)/12
         if (ii == 3 .and. prtvol /= 0) exit
         write(msg,'(12f6.3)') occ( 1+ii*12+(ikpt-1+nkpt)*nband(1):min(12+ii*12,nband(1))+(ikpt-1+nkpt)*nband(1) )
         call wrtout(std_out,msg)
       end do
     end do
     if(nkpt/=nkpt_eff) call wrtout(std_out,' newocc: prtvol=0, stop printing more k-point information')
   end if

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
!!      m_occ
!!
!! CHILDREN
!!      timab,xmpi_bcast,xmpi_sum
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

 ! Initialize the occupation function and generalized entropy function,
 ! at the beginning, or if occopt changed

 if(option==-1)then
   nptsdiv2 = nptsdiv2_def
   return
 end if

 if (occopt_prev/=occopt .or. abs(tsmear_prev-tsmear)  >tol12 .or. abs(tphysel_prev-tphysel)>tol12) then
   occopt_prev=occopt
   tsmear_prev=tsmear
   tphysel_prev=tphysel

   ! Check whether input values of tphysel tsmear and occopt are consistent
   dblsmr = 0
   if (abs(tphysel)>tol12) then
     ! Use re-smearing scheme
     if (abs(tsmear)>tol12) then
       dblsmr = 1
       ! Use FD occupations (one smearing) only with "physical" temperature tphysel
     ! CP modify
     !else if (occopt /= 3) then
     !  write(msg, '(a,i6,a)' )' tphysel /= 0, tsmear == 0, but occopt is not = 3, but ',occopt,'.'
     else if (occopt /= 3 .and. occopt/=9) then
       write(msg, '(a,i6,a)' )' tphysel /= 0, tsmear == 0, but occopt is not = 3 or 9, but ',occopt,'.'
     ! End CP modify
       ABI_ERROR(msg)
     end if
   end if

   ABI_MALLOC(entder,(-nptsdiv2_def:nptsdiv2_def))
   ABI_MALLOC(occder,(-nptsdiv2_def:nptsdiv2_def))
   ABI_MALLOC(smdder,(-nptsdiv2_def:nptsdiv2_def))
   ABI_MALLOC(workfun,(-nptsdiv2_def:nptsdiv2_def))
   ABI_MALLOC(work,(-nptsdiv2_def:nptsdiv2_def))

   ! Prepare the points on the grid
   ! limit is the value of the argument that will give 0.0 or 1.0 , with
   ! less than about 1.0d-15 error for 4<=occopt<=8, and less than about 1.0d-12
   ! error for occopt==3. It is not worth to compute the function beyond
   ! that point. Even with a less severe requirement, it is significantly
   ! larger for occopt==3, with an exponential
   ! tail, than for the other occupation functions, with a Gaussian tail.
   ! Note that these values are useful in newocc.f also.
   limit_occ=6.0_dp
   ! CP modify
   !if(occopt==3)limit_occ=30.0_dp
   if(occopt==3 .or. occopt==9)limit_occ=30.0_dp
   ! End CP modify
   if(dblsmr /= 0) then
     tratio = tsmear / tphysel
     limit_occ=30.0_dp + 6.0_dp*tratio
   end if

   ! With nptsdiv2_def=6000 (thus increm=0.001 for 4<=occopt<=8,
   ! and increm=0.005 for occopt==3, the O(1/N4) algorithm gives 1.0d-12
   ! accuracy on the stored values occfun and entfun. These, together
   ! with smdfun and xgrid_prev, need permanently about 0.67 MB, which is affordable.
   increm=limit_occ/nptsdiv2_def
   do ii=-nptsdiv2_def,nptsdiv2_def
     xgrid_prev(ii)=ii*increm
   end do

   !  ---------------------------------------------------------
   !  Ordinary (unique) smearing function
   !  ---------------------------------------------------------
   if (dblsmr == 0) then

     ! Compute the unnormalized smeared delta function between -limit_occ and +limit_occ
     ! (well, they are actually normalized ...)

     ! CP modify
     !if(occopt==3)then
     if(occopt==3 .or. occopt==9)then
     ! End CP modify
       ! Fermi-Dirac
       do ii=0,nptsdiv2_def
         xx=xgrid_prev(ii)
         smdfun_prev( ii,1)=0.25_dp/(cosh(xx/2.0_dp)**2)
         smdfun_prev(-ii,1)=smdfun_prev(ii,1)
       end do

     else if(occopt==4 .or. occopt==5)then
       ! Cold smearing of Marzari, two values of the "a" parameter being possible
       ! first value gives minimization of the bump
       if(occopt==4)aa=-.5634
       ! second value gives monotonic occupation function
       if(occopt==5)aa=-.8165

       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nptsdiv2_def
         xx=xgrid_prev(ii)
         gauss=dsqrpi*exp(-xx**2)
         smdfun_prev( ii,1)=gauss*(1.5_dp+xx*(-aa*1.5_dp+xx*(-1.0_dp+aa*xx)))
         smdfun_prev(-ii,1)=gauss*(1.5_dp+xx*( aa*1.5_dp+xx*(-1.0_dp-aa*xx)))
       end do

     else if(occopt==6)then

       ! First order Hermite-Gaussian of Paxton and Methfessel
       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nptsdiv2_def
         xx=xgrid_prev(ii)
         smdfun_prev( ii,1)=dsqrpi*(1.5_dp-xx**2)*exp(-xx**2)
         smdfun_prev(-ii,1)=smdfun_prev(ii,1)
       end do

     else if(occopt==7)then

       ! Gaussian smearing
       dsqrpi=1.0_dp/sqrt(pi)
       do ii=0,nptsdiv2_def
         xx=xgrid_prev(ii)
         smdfun_prev( ii,1)=dsqrpi*exp(-xx**2)
         smdfun_prev(-ii,1)=smdfun_prev(ii,1)
       end do

     else if(occopt==8)then

       ! Constant value of the delta function over the smearing interval, for testing purposes only.
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
       ABI_BUG(sjoin('Occopt: ', itoa(occopt),' is not allowed in getnel.'))
     end if

   else if (dblsmr /= 0) then
     !    ---------------------------------------------------------
     !    smear FD delta with occopt delta calculated in smdfun_prev
     !    ---------------------------------------------------------

     nconvd2 = 6000
     convlim = 10.0_dp
     incconv = convlim / nconvd2

     ! store smearing functions in smd1 and smd2
     ABI_MALLOC(smd1,(-nconvd2:nconvd2))
     ABI_MALLOC(smd2,(-nconvd2:nconvd2))
     ABI_MALLOC(tgrid,(-nconvd2:nconvd2))

     ! FD function in smd1( ii) and second smearing delta in smd2( ii)
     !
     ! smd1(:) contains delta_FD ( x )
     do ii=0,nconvd2
       tgrid(ii)=ii*incconv
       tgrid(-ii)=-tgrid(ii)
       tt=tgrid(ii)
       smd1( ii)=0.25_dp/(cosh(tt/2.0_dp)**2)
       smd1(-ii)=smd1(ii)
     end do

     ! check input values of occopt and fill smd2(:) with appropriate data:
     ! smd2(:) contains delta_resmear ( x )
     ! CP modify
     !if(occopt == 3) then
     if(occopt == 3 .or. occopt==9) then
     ! End CP modify
       write(msg, '(a,a)' )&
        'Occopt=3 is not allowed as a re-smearing.', &
        'Use a single FD, or re-smear with a different delta type (faster cutoff). '
       ABI_ERROR(msg)
     else if(occopt==4 .or. occopt==5)then
       ! Cold smearing of Marzari, two values of the "a" parameter being possible
       ! first value gives minimization of the bump
       if(occopt==4)aa=-.5634
       ! second value gives monotonic occupation function
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
       ABI_BUG(sjoin('Occopt: ', itoa(occopt),' is not allowed in getnel.'))
     end if

     ! Use O(1/N4) algorithm from Num Rec (see below)
     !
     ! The grid for the convoluted delta is taken (conservatively)
     ! to be that for the FD delta ie 6000 pts in [-limit_occ;limit_occ]
     ! Smearing functions are given on [-dbllim;dbllim] and the grid must
     ! superpose the normal grid on [-limit_occ:limit_occ]
     ! The maximal interval for integration of the convolution is
     ! [-dbllim+limit_occ+lim(delta2);dbllim-limit_occ-lim(delta2)] =
     ! [-dbllim+36;dbllim-36]

     ! test the smdFD function for extreme values:
     ! do jj=-nptsdiv2_def,-nptsdiv2_def
     ! do ii=-nconvd2+4,nconvd2
     ! call smdFD(xgrid_prev(jj) - tgrid(ii)*tratio, resFD)
     ! write(std_out,*) 'ii jj = ', ii,jj, ' smdFD (', xgrid_prev(jj) - tgrid(ii)*tratio, ') ', resFD
     ! end do
     ! end do

     expinc = exp(half*incconv*tratio)

     ! jj = position of point at which we are calculating smdfun_prev
     do jj=-nptsdiv2_def,nptsdiv2_def
       ! Do not care about the 8 boundary points,
       ! where the values should be extremely small anyway
       smdfun_prev(jj,1)=0.0_dp
       ! only add contribution with delta_FD > 1.0d-100
       nmaxFD = floor  (( maxFDarg+xgrid_prev(jj)) / tratio / incconv )
       nmaxFD = min (nmaxFD, nconvd2)
       nminFD = ceiling((-maxFDarg+xgrid_prev(jj)) / tratio / incconv )
       nminFD = max (nminFD, -nconvd2)

       ! Calculate the Fermi-Dirac distrib at point xgrid_prev(jj)-tgrid(ii)*tratio
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

       ! core contribution to the integral with constant weight (48)
       tmpsmdfun = 0.0_dp
       do ii=nminFD+4,nmaxFD-4
         expxo2 = expxo2*expinc
         ! tmpexpsum = 1.0_dp / (expxo2 + 1.0_dp / expxo2 )
         expx22 = expxo2*expxo2
         tmpexpsum = expxo2 / (expx22 + 1.0_dp)
         tmpsmdfun = tmpsmdfun + smd2(ii) * tmpexpsum * tmpexpsum
       end do

       ! Add on end contributions for show (both functions smd and smdFD are very small
       smdfun_prev(jj,1)=smdfun_prev(jj,1)       +48.0_dp*tmpsmdfun             &
         + 31.0_dp*smd2(nminFD+3)*resFD1 -11.0_dp*smd2(nminFD+2)*resFD2 &
         +  5.0_dp*smd2(nminFD+1)*resFD3 -       smd2(nminFD)*resFD4

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

       ! Contribution above
       smdfun_prev(jj,1)=smdfun_prev(jj,1)                                      &
         + 31.0_dp*smd2(nmaxFD-3)*resFD1  -11.0_dp*smd2(nmaxFD-2)*resFD2 &
         +  5.0_dp*smd2(nmaxFD-1)*resFD3  -       smd2(nmaxFD)*resFD4
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
     ! thdmom=increm * thdmom / 48.0_dp
     !
     ! smom1  = second moment of delta in smd1(:)
     ! smom2  = second moment of delta in smd2(:)
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

     ABI_FREE(tgrid)
     ABI_FREE(smd1)
     ABI_FREE(smd2)

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

     ! Extended trapezoidal rule (4.1.11), taken in a cumulative way
     do ii=-nptsdiv2_def+1,nptsdiv2_def
       occfun_prev(ii,1)=occfun_prev(ii-1,1)+increm*(smdfun_prev(ii,1)+smdfun_prev(ii-1,1))/2.0_dp
       entfun_prev(ii,1)=entfun_prev(ii-1,1)+increm*&
&       ( -xgrid_prev(ii)*smdfun_prev(ii,1) -xgrid_prev(ii-1)*smdfun_prev(ii-1,1) )/2.0_dp
     end do

   else if(algo==3)then

     ! Derived from (4.1.12). Converges as O(1/N3).
     ! Do not care about the following points,
     ! where the values are extremely small anyway
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

     ! Derived from (4.1.14)- alternative extended Simpsons rule. Converges as O(1/N4).
     ! Do not care about the following points,
     ! where the values are extremely small anyway
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

   ! Normalize the functions (actually not needed for occopt=3..7)
   factor=1.0_dp/occfun_prev(nptsdiv2_def,1)
   smdfun_prev(:,1)=smdfun_prev(:,1)*factor
   occfun_prev(:,1)=occfun_prev(:,1)*factor
   entfun_prev(:,1)=entfun_prev(:,1)*factor

   !  Compute the cubic spline fitting of the smeared delta function
   yp1=0.0_dp ; ypn=0.0_dp
   workfun(:)=smdfun_prev(:,1)
   call spline(xgrid_prev, workfun, (2*nptsdiv2_def+1), yp1, ypn, smdder)
   smdfun_prev(:,2)=smdder(:)

   ! Compute the cubic spline fitting of the occupation function
   yp1=0.0_dp ; ypn=0.0_dp
   workfun(:)=occfun_prev(:,1)
   call spline(xgrid_prev, workfun, (2*nptsdiv2_def+1), yp1, ypn, occder)
   occfun_prev(:,2)=occder(:)

   ! Compute the cubic spline fitting of the entropy function
   yp1=0.0_dp ; ypn=0.0_dp
   workfun(:)=entfun_prev(:,1)
   call spline(xgrid_prev, workfun, (2*nptsdiv2_def+1), yp1, ypn, entder)
   entfun_prev(:,2)=entder(:)

   ABI_FREE(entder)
   ABI_FREE(occder)
   ABI_FREE(smdder)
   ABI_FREE(work)
   ABI_FREE(workfun)

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
!! that depend on the difference between occupation numbers and eigenvalues.
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
!!      m_dfpt_nstwf,m_dfpt_scfcv,m_dfpt_vtorho
!!
!! CHILDREN
!!      timab,xmpi_bcast,xmpi_sum
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

 rocceig(:,:) = zero

 do ibandk=1,nband_k
   do ibandkq=1,nband_k
     diffeig=eig0_kq(ibandkq)-eig0_k(ibandk)
     diffocc=occ_kq(ibandkq)-occ_k(ibandk)

     if( abs(diffeig) > tol5 ) then
       ratio=diffocc/diffeig
     else
       if(occopt<3)then
         ! In a non-metallic case, if the eigenvalues are degenerate,
         ! the occupation numbers must also be degenerate, in which
         ! case there is no contribution from this pair of bands
         if( abs(diffocc) > tol5 ) then
           write(msg,'(a,a,a,a,a,a,a,2(a,i4,a,es16.6,a,es16.6,a,a),a)' ) &
           'In a non-metallic case (occopt<3), for a RF calculation,',ch10,&
           'if the eigenvalues are degenerate,',' the occupation numbers must also be degenerate.',ch10,&
           'However, the following pair of states gave :',ch10,&
           'k -state, band number',ibandk,', occ=',occ_k(ibandk),'eigenvalue=',eig0_k(ibandk),',',ch10,&
           ' kq-state, band number',ibandkq,', occ=',occ_kq(ibandkq),', eigenvalue=',eig0_kq(ibandkq),'.',ch10,&
           'Action: change occopt, consistently, in GS and RF calculations.'
           ABI_ERROR(msg)
         end if
         ratio=0.0_dp
       else
         ! In the metallic case, one can compute a better approximation of the
         ! ratio by using derivatives doccde
         ratio=0.5_dp*(doccde_kq(ibandkq)+doccde_k(ibandk))
         ! write(std_out,*)' occeig : ibandkq,doccde_kq(ibandkq)',ibandkq,doccde_kq(ibandkq)
         ! write(std_out,*)'          ibandk ,doccde_k (ibandk )',ibandk,doccde_k(ibandk)
       end if
     end if

     ! Here, must pay attention to the smallness of some coefficient
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
!!  Fermi-Dirac statistics: 1 / [(exp((e - mu)/ KT) + 1]
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

!!****f* m_occ/occ_dfde
!! NAME
!!  occ_dfde
!!
!! FUNCTION
!!  Derivative of Fermi-Dirac statistics: - (exp((e - mu)/ KT) / KT[(exp((e - mu)/ KT) + 1]^2
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

elemental real(dp) function occ_dfde(ee, kT, mu)

!Arguments ------------------------------------
 real(dp),intent(in) :: ee, kT, mu

!Local variables ------------------------------
 real(dp) :: ee_mu,arg
! *************************************************************************

 ee_mu = ee - mu

 ! 1 kelvin [K] = 3.16680853419133E-06 Hartree
 if (kT > tol6) then
   arg = ee_mu / kT
   if (arg > maxDFDarg) then
     occ_dfde = zero
   else if (arg < -maxDFDarg) then
     occ_dfde = zero
   else
     occ_dfde = - exp(arg) / (exp(arg) + one)**2 / kT
   end if
 else
   occ_dfde = zero
 end if

end function occ_dfde
!!***

!----------------------------------------------------------------------

!!****f* m_occ/occ_be
!! NAME
!!  occ_be
!!
!! FUNCTION
!!   Bose-Einstein statistics  1 / [(exp((e - mu)/ KT) - 1]
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
!!   Derivative of Bose-Einstein statistics  (exp((e - mu)/ KT) / KT[(exp((e - mu)/ KT) - 1]^2
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
!! fermih= fermi energy of thermalized excited holes when occopt = 9 ! CP added
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
!!      m_epjdos,m_occ
!!
!! CHILDREN
!!      timab,xmpi_bcast,xmpi_sum
!!
!! SOURCE

subroutine dos_hdr_write(deltaene,eigen,enemax,enemin,fermie,fermih,mband,nband,nene,&
                         nkpt,nsppol,occopt,prtdos,tphysel,tsmear,unitdos)
! CP modified arguments list and added fermih

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkpt,nsppol,occopt,prtdos,unitdos,nene
 ! CP modify
 real(dp),intent(in) :: fermie,fermih,tphysel,tsmear
 ! End CP modify
 real(dp),intent(in) :: deltaene,enemax,enemin
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *************************************************************************

 ! Write the DOS file
 write(msg, '(7a,i2,a,i5,a,i4)' ) "#",ch10, &
  '# ABINIT package : DOS file  ',ch10,"#",ch10,&
  '# nsppol =',nsppol,', nkpt =',nkpt,', nband(1)=',nband(1)
 call wrtout(unitdos, msg)

 if (any(prtdos== [1,4])) then
   write(msg, '(a,i2,a,f6.3,a,f6.3,a)' )  &
    '# Smearing technique, occopt =',occopt,', tsmear=',tsmear,' Hartree, tphysel=',tphysel,' Hartree'
 else
   write(msg, '(a)' ) '# Tetrahedron method '
 end if
 call wrtout(unitdos, msg)

 if (mband*nkpt*nsppol>=3) then
   write(msg, '(a,3f8.3,2a)' )'# For identification : eigen(1:3)=',eigen(1:3),ch10,"#"
 else
   write(msg, '(a,3f8.3)' ) '# For identification : eigen=',eigen
   write(msg, '(3a)')trim(msg),ch10,"#"
 end if
 call wrtout(unitdos, msg)

 ! CP modified
 !write(msg, '(a,f16.8)' ) '# Fermi energy : ', fermie
 if (occopt == 9) then
    write(msg, '(a,f16.8, f16.8)' ) '# Fermi energy for electrons and holes ', fermie, fermih
 else
    write(msg, '(a,f16.8)' ) '# Fermi energy : ', fermie
 end if
 ! End CP modify
 call wrtout(unitdos, msg)

 if (prtdos==1) then
   write(msg, '(5a)' ) "#",ch10,&
    '# The DOS (in electrons/Hartree/cell) and integrated DOS (in electrons/cell),',&
    ch10,'# as well as the DOS with tsmear halved and doubled, are computed,'

 else if (prtdos==2)then
   write(msg, '(3a)' ) "#",ch10,&
    '# The DOS (in electrons/Hartree/cell) and integrated DOS (in electrons/cell) are computed,'

 else if (any(prtdos == [3, 4])) then
   write(msg, '(5a)' ) "#",ch10,&
    '# The local DOS (in electrons/Hartree for one atomic sphere)',ch10,&
    '# and integrated local DOS (in electrons for one atomic sphere) are computed.'

 else if (prtdos==5)then
   write(msg, '(9a)' ) "#",ch10,&
   '# The spin component DOS (in electrons/Hartree/cell)',ch10,&
   '# and integrated spin component DOS (in electrons/cell) are computed.',ch10,&
   '# Remember that the wf are eigenstates of S_z and S^2, not S_x and S_y',ch10,&
   '#   so the latter will not always sum to 0 for paired electronic states.'
 end if
 call wrtout(unitdos, msg)

 write(msg, '(a,i5,a,a,a,f9.4,a,f9.4,a,f8.5,a,a,a)' )&
  '# at ',nene,' energies (in Hartree) covering the interval ',ch10,&
  '# between ',enemin,' and ',enemax,' Hartree by steps of ',deltaene,' Hartree.',ch10,"#"
 call wrtout(unitdos, msg)

 if (prtdos==1) then
   write(msg, '(a,a)' )&
    '#       energy        DOS       Integr. DOS   ','     DOS           DOS    '
   call wrtout(unitdos,msg)

   write(msg, '(a)' )&
    '#                                              (tsmear/2)    (tsmear*2) '
   call wrtout(unitdos,msg)
 else
   write(msg, '(a)' ) '#       energy        DOS '
 end if

end subroutine dos_hdr_write
!!***

end module m_occ
!!***
