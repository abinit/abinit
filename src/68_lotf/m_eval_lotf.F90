!!****m* ABINIT/eval_lotf
!! NAME
!! eval_lotf
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2005-2020 ABINIT group (MMancini)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module eval_lotf

 use defs_basis
 use m_errors

 implicit none
 public  

 public ::                 &
   eval_forces_u_n,        &
   phi_n_calc,             &
   calc_coord2,            &
   eval_forces_u_n_2,      &
   eval_force_devs_new_d,  &
   upd_lis0,               &
   tuneparms

contains
!!***

 !!****f* eval_lotf/phi_n_calc
 !! NAME
 !! phi_n_calc
 !!
 !! FUNCTION
 !!
 !! INPUTS
 !! PARENTS
!!      m_lotf
!!
 !! CHILDREN
!!      dist_pbc,wrtout
!!
 !! SOURCE
 !--Similar to SWCALC but without the triplets (used for the pair potential part and for the coordination) 
 subroutine phi_n_calc(alpha_dum,nneig,nlist,r0,rv,epot_dum,&
   &                        forcv,coordatom_dum,alpha_fdum)

  use defs_param_lotf,only : lotfvar
  use bond_lotf,only : nbondex,tafit,imat,ibmat,ibmat_large
  use glue_lotf,only : dphi,rcphi,rmrho,glue_pair,glue_pair_devs,calc_coord,calc_coord_new_d
  use pbc_lotf,only : dist_pbc

  ! INPUT/OUTPUT
  real(dp),intent(in) :: alpha_dum(3,nbondex)
  integer ::  nneig,  iat, jat
  integer ::  nlist(0:nneig)
  real(dp),intent(in) :: r0(3)
  real(dp) ::  forcv(3,0:nneig)
  real(dp) ::  RDV(3,nneig)
  real(dp) ::  r_au(nneig)
  real(dp) ::  rv(3,nneig)
  real(dp) ::  coordatom_dum
  real(dp) ::  alpha_fdum(3,2,nbondex,1)

  !Local ---------------------------
  integer ::   i1,  i0, i3, nbond
  real(dp) ::  epot_dum, epot_2, drcphi
  real(dp) ::  rcut2_pair, rcut2_rho
  real(dp) ::  rref(3), fp(3), dfp(3,2)

! *************************************************************************

  drcphi = rcphi - dphi
  rref(:) = zero
  forcv(:,0:nneig) = zero
  epot_dum = zero
  coordatom_dum = zero

  iat  = nlist(0)


  run_over_neighbours: do i1 = 1, nneig  ! Run over neighbours
    !--PAIR: POTENTIAL, FORCES, DERIVATIVES
    call dist_pbc(rv(:,i1),r0,r_au(i1),RDV(:,i1))

    jat = nlist(i1)
    if (tafit(jat).and.(jat > iat)) then ! Only for fit pairs + NO doUBLE COUNTING

      i0 = imat(iat) ! Index of atom iat in the fitting indexing (1,nfit)
      nbond = ibmat(i1,i0)

      !--Cutoff radius for the pair potential
      rcut2_pair = (alpha_dum(1,nbond) + drcphi)**2 

      if ((r_au(i1) < rcut2_pair)) then 

        call glue_pair_devs(alpha_dum(:,nbond),RDV(:,i1),r_au(i1),epot_2,fp,dfp)

        epot_dum = epot_dum + epot_2

        forcv(:,0)  = forcv(:,0)  + fp(:)
        forcv(:,i1) = forcv(:,i1) - fp(:)

        alpha_fdum(:,1,nbond,1) =  dfp(:,1)
        alpha_fdum(:,2,nbond,1) =  dfp(:,2)
      endif

    else !--We have to count the forces on the fit atoms by the non-fit neighbours

      rcut2_pair = rcphi**2
      if ( (r_au(i1) < rcut2_pair).and.(jat > iat) ) then 
        !--Cutoff pair potential + NO MULTIPLE COUNTING
        call glue_pair(RDV(1,i1),r_au(i1),epot_2,fp)
        epot_dum = epot_dum + epot_2

        forcv(:,0)  = forcv(:,0)  + fp(:)
        forcv(:,i1) = forcv(:,i1) - fp(:)
      endif
    endif ! Pair potential and pair forces and derivs only for fit pairs

    !--COORDINATION 
    select case(lotfvar%classic)
    case(5)
      !--Cutoff radius for the density 
      rcut2_rho  = rmrho**2 

      !--Cutoff density
      if (r_au(i1) < rcut2_rho) then
        call calc_coord(r_au(i1),coordatom_dum)
      end if

    case(6) 
      if (tafit(jat)) then

        i0 = imat(iat) ! Index of atom iat in the fitting indexing (1,nfit)
        i3 = imat(jat)
        nbond = ibmat_large(i3,i0)

        !--Cutoff radius for the density
        rcut2_rho  = (rmrho+alpha_dum(1,nbond)-dphi)**2 

        !--Cutoff density
        if (r_au(i1) < rcut2_rho) then
          call calc_coord_new_d(r_au(i1),alpha_dum(1,nbond),coordatom_dum)
        end if

      else ! tafit = .false
        !--Cutoff radius for the density
        rcut2_rho  = rmrho**2 
        if (r_au(i1) < rcut2_rho) then
          call calc_coord(r_au(i1),coordatom_dum)
        end if

      endif ! Fit/NonFit
    end select
  enddo run_over_neighbours

 end subroutine phi_n_calc
 !!***


 !!****f* eval_lotf/calc_coord2
 !! NAME
 !! calc_coord2
 !!
 !! FUNCTION
 !!
 !! INPUTS
 !! PARENTS
!!      m_lotf
!!
 !! CHILDREN
!!      dist_pbc,wrtout
!!
 !! SOURCE
 subroutine calc_coord2(nneig,r0,rv,coordatom_dum)

  USE pbc_lotf,only : dist_pbc
  USE GLUE_LOTF,only : calc_coord,rmrho
  integer,intent(in) :: nneig
  real(dp),intent(out) :: coordatom_dum
  real(dp),intent(in) :: r0(3),rv(3,nneig)

  !Local ---------------------------
  integer :: ii
  real(dp) :: r_au(nneig),RDV(3,nneig),rcut2_rho

! *************************************************************************

  !--Cutoff radius for the density
  rcut2_rho  = rmrho**2 

  !--Run over neighbours
  do ii=1,nneig  
    call dist_pbc(rv(:,ii),r0,r_au(ii),RDV(:,ii))

    !--Cutoff density
    if (r_au(ii) < rcut2_rho) then
      call calc_coord(r_au(ii),coordatom_dum)
    end if
  enddo
 end subroutine calc_coord2
 !!***

 !!****f* eval_lotf/eval_forces_U_n
 !! NAME
 !! eval_forces_U_n
 !!
 !! FUNCTION
 !!
 !! INPUTS
 !! PARENTS
!!      m_lotf
!!
 !! CHILDREN
!!      dist_pbc,wrtout
!!
 !! SOURCE
 subroutine eval_forces_U_n(nneig,nlist,r0,rv,up_list,forc_dum2)
  use pbc_lotf,only : dist_pbc
  use defs_param_lotf,only : lotfvar
  use GLUE_LOTF, only :  rmrho,calc_rhop
  integer,intent(in) ::  nneig,nlist(0:lotfvar%nneigx)
  real(dp),intent(in) :: r0(3),rv(3,nneig),up_list(lotfvar%natom)
  real(dp),intent(out):: forc_dum2(3)

  !Local ---------------------------
  integer :: jat,iat
  real(dp) :: U_tot,r_au,RDV(3),r_st,rcut2_rho,rhop_dum

! *************************************************************************

  !--Cutoffradius for the density
  rcut2_rho  = rmrho**2 
  U_tot = zero
  forc_dum2(:) = zero

  iat = nlist(0)
  do jat=1,nneig
    call dist_pbc(rv(:,jat),r0,r_au,RDV)
    if(r_au < rcut2_rho) then

      U_tot = up_list(iat) + up_list(nlist(jat))
      r_st = sqrt(r_au)
      call calc_rhop(r_st,rhop_dum)

      U_tot = U_tot * rhop_dum
      forc_dum2(:) = forc_dum2(:) + U_tot * RDV(:) / r_st
    endif
  enddo
 end subroutine eval_forces_U_n
 !!***


 !!****f* eval_lotf/eval_forces_U_n_2
 !! NAME
 !! eval_forces_U_n_2
 !!
 !! FUNCTION
 !!
 !! INPUTS
 !! PARENTS
!!      m_lotf
!!
 !! CHILDREN
!!      dist_pbc,wrtout
!!
 !! SOURCE
 subroutine eval_forces_U_n_2(alpha_dum,nneig,nlist,&
   r0,rv,up_list,rho_p_sum,forc_dum2)

  use pbc_lotf,only : dist_pbc
  use GLUE_LOTF,only : rmrho,dphi,rhop_value
  use defs_param_lotf,only : lotfvar
  use bond_lotf,only : nbondex,tafit,imat,ibmat_large

  ! Input/Output variables
  real(dp),intent(in) ::  alpha_dum(3,nbondex)
  integer,intent(in) ::  nneig,nlist(0:lotfvar%nneigx)
  real(dp),intent(in) ::  r0(3),rv(3,nneig),up_list(lotfvar%natom)
  real(dp),intent(out) ::  rho_p_sum(3),forc_dum2(3)

  !Local ---------------------------
  integer :: iat,ii,jat,i0,nbond,i2
  real(dp) :: U_tot,r_au,RDV(3),r_st,rcut2_rho
  real(dp) :: rho_p_esc

! *************************************************************************

  U_tot = zero
  rho_p_sum(:) = zero

  iat = nlist(0)

  !--(0) RHOP_SUM(iat) and FORC_DUM2(iat) (both need to be summed up before calculating FORCE DERIVATIVES)
  do ii=1,nneig

    jat=nlist(ii)
    call dist_pbc(rv(:,ii),r0,r_au,RDV)
    r_st = sqrt(r_au)

    if (tafit(jat)) then

      !--(0a) RHOP_SUM
      i0 = imat(iat) 
      i2 = imat(jat)
      nbond = ibmat_large(i2,i0)

      !--Cutoff radius for the density
      rcut2_rho  = (rmrho+alpha_dum(1,nbond)-dphi)**2 
      if(r_au < rcut2_rho) then
        call rhop_value(r_st,alpha_dum(1,nbond),rho_p_esc)
      end if
    else

      rcut2_rho = rmrho**2
      if(r_au < rcut2_rho) then
        call rhop_value(r_st,dphi,rho_p_esc)
      end if
    endif

    if (r_au < rcut2_rho) then 

      rho_p_sum(:) = rho_p_sum(:) + rho_p_esc * RDV(:) / r_st 

      !--(0b) FORCES 
      U_tot = up_list(iat) + up_list(jat)
      U_tot = U_tot * rho_p_esc

      !--2-body + many-body
      forc_dum2(:) = forc_dum2(:) + U_tot * RDV(:) / r_st
    endif
  enddo
 end subroutine eval_forces_U_n_2
 !!***



 !!****f* eval_lotf/eval_force_devs_new_d
 !! NAME
 !! eval_force_devs_new_d
 !!
 !! FUNCTION
 !!
 !! INPUTS
 !! PARENTS
!!      m_lotf
!!
 !! CHILDREN
!!      dist_pbc,wrtout
!!
 !! SOURCE
 subroutine eval_force_devs_new_d(alpha_dum,nneig,nlist,neig2,nlist2,&
   r0,rv,rv2,up_list,upp_list,fact2,ffit,&
   forc_dum2,rho_p_sum,dcost_dalpha)

  use pbc_lotf,only : dist_pbc
  use GLUE_LOTF,only : rmrho,dphi,rho_devs
  use defs_param_lotf,only : lotfvar
  use bond_lotf,only : nbondex,nfitmax,tafit,imat,ibmat_large
  implicit none

  !Arguments ------------------------
  real(dp),intent(in) ::  r0(3)
  real(dp) ::  alpha_dum(3,nbondex)
  integer ::  nneig,nlist(0:lotfvar%nneigx)
  integer ::  neig2(lotfvar%nneigx),nlist2(lotfvar%nneigx,lotfvar%nneigx)
  real(dp) ::  rv(3,lotfvar%nneigx),rv2(3,lotfvar%nneigx,lotfvar%nneigx)
  real(dp) ::  up_list(lotfvar%natom),upp_list(lotfvar%natom)
  real(dp) ::  fact2(nfitmax),ffit(3,nfitmax)
  real(dp) ::  forc_dum2(3,0:nfitmax)
  real(dp) ::  dcost_dalpha(3,nbondex)
  !Local ---------------------------
  integer :: iat,ii,jat,i0,nbond,ixyz
  real(dp) :: r_au,RDV(3),r_st,rcut2_rho,U_tot0
  real(dp) :: rho_neigh_d,rho_neigh_p,rho_neigh_pd
  real(dp) :: rho_p_sum(3,nfitmax) 
  real(dp) :: dFi_dbond_ij(3)
  real(dp) :: r_au_ik,r_st_ik,RDV_ik(3)
  real(dp) :: r_au_jk,r_st_jk,RDV_jk(3)
  real(dp) :: rcut2_rho_ik
  integer :: nbond_ik
  real(dp) :: rho_d_ik,rho_p_ik,rho_pd_ik
  real(dp) :: rcut2_rho_jk
  integer :: nbond_jk
  real(dp) :: rho_d_jk,rho_p_jk,rho_pd_jk
  real(dp):: dFi_dbond_jk(3)
  integer :: i2,i3,i4,kat
  integer :: iunique(lotfvar%nneigx*lotfvar%nneigx),indunique

! *************************************************************************

  U_tot0 = zero
  dFi_dbond_ij(:) = zero
  dFi_dbond_jk(:) = zero
  iat = nlist(0)

  indunique = 0
  iunique(:) = 0

  !--1st NEIGHBOURS
  first_neighbours : do ii=1,nneig 
    jat=nlist(ii)

    call dist_pbc(rv(:,ii),r0,r_au,RDV)
    r_st = sqrt(r_au)
    i0 = imat(iat)
    i2 = imat(jat)

    rho_neigh_d = zero
    rho_neigh_p = zero
    rho_neigh_pd = zero

    !--Fit
    if (tafit(jat)) then  

      nbond = ibmat_large(i2,i0)
      rcut2_rho  = (rmrho+alpha_dum(1,nbond)-dphi)**2 

      !--Cutoff 
      if (r_au < rcut2_rho) then 
        call rho_devs(r_au,alpha_dum(1,nbond),rho_neigh_d,&
          rho_neigh_p,rho_neigh_pd)

        U_tot0 = up_list(iat) + up_list(jat)

        dFi_dbond_ij(:) = upp_list(iat)*rho_neigh_d*rho_p_sum(:,i0) + &
          upp_list(jat)*rho_neigh_d*rho_neigh_p*RDV(:)/r_st+&
          U_tot0*rho_neigh_pd*RDV(:)/r_st

        do ixyz=1,3
          !--Atom iat
          dcost_dalpha(1,nbond) = dcost_dalpha(1,nbond) - &
            fact2(i0) * (forc_dum2(ixyz,i0) - ffit(ixyz,i0)) *&
            dFi_dbond_ij(ixyz)
        enddo

      endif
    endif !--jat FIT atom

    !--Debug
    go to 1002

    !--2nd NEIGHBOURS
    second_neighbours : do i3=1,neig2(ii) 

      kat=nlist2(i3,ii)

      ! Make sure that kat is 'NEW' (not in the neigbours of iat neither
      ! in previously explored neigbours lists...)

      if (tafit(jat).AND.tafit(kat).AND.kat /= iat) then
        call dist_pbc(rv2(:,i3,ii),r0,r_au_ik,RDV_ik)
        r_st_ik = sqrt(r_au_ik)
        call dist_pbc(rv2(:,i3,ii),rv(:,ii),r_au_jk,RDV_jk)
        r_st_jk = sqrt(r_au_jk)

        i4 = imat(kat)

        rho_d_ik = zero
        rho_p_ik = zero
        rho_pd_ik = zero

        !--FIT pair iat-kat 
        if (tafit(kat)) then 

          nbond_ik = ibmat_large(i4,i0)
          rcut2_rho_ik  = (rmrho+alpha_dum(1,nbond_ik)-dphi)**2
          !--Cutoff 
          if (r_au_ik < rcut2_rho_ik) then
            call rho_devs(r_au_ik,alpha_dum(1,nbond_ik),rho_d_ik,rho_p_ik,rho_pd_ik)
          endif ! cutoffs

        endif

        rho_d_jk = zero
        rho_p_jk = zero
        rho_pd_jk = zero

        !--FIT pair jat-kat 
        if (tafit(kat).AND.tafit(jat)) then 

          nbond_jk = ibmat_large(i4,i0)
          rcut2_rho_jk  = (rmrho+alpha_dum(1,nbond_jk)-dphi)**2

          !--Cutoff + no multiple counting
          if (r_au_jk < rcut2_rho_jk.AND.kat > iat) then 
            call rho_devs(r_au_jk,alpha_dum(1,nbond_jk),rho_d_jk,&
              rho_p_jk,rho_pd_jk)
          endif
        endif ! FIT pair jat-kat

        if(r_au_jk < rcut2_rho_jk) then
          if(r_au_ik < rcut2_rho_ik) then
            dFi_dbond_jk(:) = dFi_dbond_jk(:) + upp_list(kat) * rho_d_jk * rho_p_ik *&
              RDV_ik(:)/r_st_ik 
          endif ! ik cutoff

          if(r_au < rcut2_rho) then
            dFi_dbond_jk(:) = dFi_dbond_jk(:) + upp_list(jat) * rho_d_jk * rho_neigh_p *&
              RDV(:)/r_st 
          endif ! ij cutoff
        endif ! jk cutoff 

        !--Cutoff 
        if (r_au_jk < rcut2_rho_jk) then
          do ixyz=1,3
            !--Atom iat
            dcost_dalpha(1,nbond_jk) = dcost_dalpha(1,nbond_jk) + &
              fact2(i0) * (forc_dum2(ixyz,i0) - ffit(ixyz,i0)) *&
              dFi_dbond_jk(ixyz)
          enddo
        endif

      endif ! atoms jat and kat are FIT and kat is not iat again 

    enddo second_neighbours

1002 continue

  enddo first_neighbours

 end subroutine eval_force_devs_new_d
 !!***

 !!****f* eval_lotf/upd_lis0
 !! NAME
 !! upd_lis0
 !!
 !! FUNCTION
 !!  update neigbours relations
 !!
 !! lotfvar%nneigx : max number of neighbours
 !! tau0(3,natom) : atomic positions 
 !! neighl(lotfvar%nneigx,natom) : list of neighbours 
 !! nneig(natom) : number of neighbours 
 !! niter  : iteration number (itime) 
 !!
 !! upgraded to use the linked cell method (Knuth) 
 !! 
 !! INPUTS
 !! PARENTS
!!      m_lotf
!!
 !! CHILDREN
!!      dist_pbc,wrtout
!!
 !! SOURCE
 subroutine  upd_lis0(tau0,neighl,nneig,niter)   
  use defs_param_lotf,only : lotfvar
  use work_var_lotf,only : rcut_nbl
  use tools_lotf,only : icf
  use pbc_lotf,only : dist_pbc,r2,rd,pbc_bb_proj,pbc_bb_contract,pbc_aa_contract
  implicit none

  !Arguments ------------------------
  integer,intent(in)   :: niter  
  integer,intent(out)  :: nneig(lotfvar%natom), neighl(lotfvar%nneigx,lotfvar%natom)
  real(dp),intent(inout) :: tau0(3,lotfvar%natom)
  !Local --------------------------- 
  integer,parameter :: icellx_a=9500
  integer  :: icell_tot,icell
  integer  :: icnr, icnr2, icn, icn2
  integer  :: ic
  integer  :: ix, iy, iz
  integer  :: n_count, iat,i 
  integer  :: ica(3) 
  integer  :: head(icellx_a)  
  integer  :: list(lotfvar%natom)  
  integer  :: neigh_cel(icellx_a,13) 
  real(dp) :: r3(3), rcut2 
  real(dp) :: rdum(3)
  real(dp) :: raa(3) 
  character(len=500) :: message

! *************************************************************************

  n_count = niter - lotfvar%n0
  write(message,'(a,4i8)')&
    &  'Checking Input',n_count,niter,mod(n_count,lotfvar%nitex),lotfvar%nitex
  call wrtout(std_out,message,'COLL')

  rcut2   = rcut_nbl*rcut_nbl     

  !--(1) gets all atoms within the "same" repeated cell zone
  !     (from -0.5 to 0.5 in relative units) 
  !at this part is changed to take into account symetries other than orthorombic.
  r3(:) = zero

  !--length of the cell vectors :   
  raa(:) = pbc_aa_contract()

  !--put tau0 and taum in the cell -0.5/+0.5 for tau0
  if(lotfvar%version==1) then 
    do i =1, lotfvar%natom
      !--compute rd
      call dist_pbc(tau0(:,i),r3)
      rdum(:)  =  rd(:) - tau0(:,i)
      tau0(:,i) =  tau0(:,i) + rdum(:)
      ! taum(:,i) =  taum(:,i) + rdum(:)
    enddo
  elseif(lotfvar%version==2) then  
    do i =1, lotfvar%natom
      !--compute rd
      call dist_pbc(tau0(:,i),r3)
      rdum(:)  =  rd(:) - tau0(:,i)
      tau0(:,i) =  tau0(:,i) + rdum(:)
    enddo
  endif ! lotfvar%version says if taum is touched

  !--clears neighbour lists
  nneig(:) = 0 
  neighl(:,:) = 0


  !--OK, NOW KNUTH TRICK 
  ! determines the cell partition of the run 
  ! it means we cut volumes with edges // to the cell vectors   
  ! but with length divided by an integer ica(1:3)

  ica(:) = int(raa(:)/sqrt(4*rcut2)) 


  write(message,'(a,3f12.8,4a,2f12.8)')&
    &     ' raa = ', raa(:),ch10,&
    &     ' Neighbour List Cutoff: ',ch10,&
    &     '   rcut2, sqrt(4*rcut2) ',rcut2,sqrt(4*rcut2) 
  call wrtout(std_out,message,'COLL')

  where(ica < 1) ica = 1
  icell_tot =  product(ica)

  write(message,'(3a,2i8,2a,3i8)')&
    & 'UPDLIS0 : SYSTEM SUBCELL DISTRIBUTION: ',ch10,&
    & 'TOT. NO. OF CELLS (& max.) : ',icell_tot, icellx_a,ch10,&
    & 'ICX,ICY,ICZ = ',ica(1),ica(2),ica(3)
  call wrtout(std_out,message,'COLL')


  if(icell_tot > icellx_a) then 
    write(message,'(a,i8,2a)')&
      & 'ERROR IN UPD_LIS0: ICELL_TOT = ',icell_tot,ch10,&
      & 'ERROR IN UPD_LIS0: RAISE ICELLX_A'
    MSG_ERROR(message)
  endif

  !-- clears head vector & constructs head & list      
  head(:icellx_a) = 0

  do i=1,lotfvar%natom
    rdum = pbc_bb_proj(tau0(:,i))
    icell = 1 + mod(int( (rdum(1)+0.5)* float(ica(1)) ),ica(1)) &
      + mod(int( (rdum(2)+0.5)* float(ica(2)) ),ica(2)) * ica(1) & 
      + mod(int( (rdum(3)+0.5)* float(ica(3)) ),ica(3)) * ica(1) * ica(2)

    list(i) = head(icell) 
    head(icell) = i 
  enddo

  !--Constructs the 13 neighbours of each cell 
  do  IZ = 1, ica(3)
    do  IY = 1, ica(2)
      do  IX = 1, ica(1)
        ic               = icf(ix  ,iy,iz    ,ica(1),ica(2),ica(3))
        neigh_cel(ic,1)  = icf(ix+1,iy,iz    ,ica(1),ica(2),ica(3)) 
        neigh_cel(ic,2)  = icf(ix+1,iy+1,iz  ,ica(1),ica(2),ica(3)) 
        neigh_cel(ic,3)  = icf(ix  ,iy+1,iz  ,ica(1),ica(2),ica(3)) 
        neigh_cel(ic,4)  = icf(ix-1,iy+1,iz  ,ica(1),ica(2),ica(3)) 
        neigh_cel(ic,5)  = icf(ix+1,iy  ,iz-1,ica(1),ica(2),ica(3)) 
        neigh_cel(ic,6)  = icf(ix+1,iy+1,iz-1,ica(1),ica(2),ica(3)) 
        neigh_cel(ic,7)  = icf(ix  ,iy+1,iz-1,ica(1),ica(2),ica(3)) 
        neigh_cel(ic,8)  = icf(ix-1,iy+1,iz-1,ica(1),ica(2),ica(3)) 
        neigh_cel(ic,9)  = icf(ix+1,iy  ,iz+1,ica(1),ica(2),ica(3)) 
        neigh_cel(ic,10) = icf(ix+1,iy+1,iz+1,ica(1),ica(2),ica(3)) 
        neigh_cel(ic,11) = icf(ix  ,iy+1,iz+1,ica(1),ica(2),ica(3)) 
        neigh_cel(ic,12) = icf(ix-1,iy+1,iz+1,ica(1),ica(2),ica(3)) 
        neigh_cel(ic,13) = icf(ix  ,iy  ,iz+1,ica(1),ica(2),ica(3)) 
      enddo
    enddo
  enddo

  !--Safety loops to avoid repetitions

  !--(1) to avoid having twice the same neigh. cell 
  do icell = 1,icell_tot
    do icnr = 1,13
      icn =  neigh_cel(icell,icnr)
      do icnr2 = icnr+1,13
        icn2 =  neigh_cel(icell,icnr2)
        if(icn2==icn)  neigh_cel(icell,icnr2) = icell 
      enddo
    enddo
  enddo

  !--(2) to avoid counting twice the interaction between the same 
  !     couple of  neigh. cells 
  do icell = 1,icell_tot
    do icnr = 1,13
      icn =  neigh_cel(icell,icnr)
      do icnr2 = 1,13
        icn2 =  neigh_cel(icn,icnr2)
        if(icn2==icell) neigh_cel(icn,icnr2) = icn 
      enddo
    enddo
  enddo

  !--(3) constructs neighbour list looking through neighbour cells only
  do icell = 1,icell_tot
    iat = head(icell) 
    do while(iat > 0)
      i = list(iat) 
      do while(i  > 0) 
        if(i==iat) MSG_ERROR("i==iat")
        call dist_pbc(tau0(:,i),tau0(:,iat))
        if(r2 < rcut2) then 
          nneig(iat) = nneig(iat) + 1
          nneig(i)   = nneig(i)   + 1
          neighl(nneig(iat),iat) = i 
          neighl(nneig(i)  ,i  ) = iat 
        endif  ! distance check 
        i = list(i) 
      enddo

      if(ANY(nneig(:) > lotfvar%nneigx))  then 
        write(message,'(a,i8,a)')&
          'UPD_LIS CLASSIC: max no. of neighbours: ',lotfvar%nneigx,' is too small'
        MSG_ERROR(message)
      endif


      !   mask to avoid self interaction within the same cell
      !   considered as a neighbor and with different cells more than once
      do icnr = 1,13
        icn = neigh_cel(icell,icnr)
        if(icn==icell) cycle
        i = head(icn) 

        do while(i>0)
          if(i==iat) MSG_ERROR("i==iat")
          call dist_pbc(tau0(:,i),tau0(:,iat))
          if(r2 < rcut2) then 
            nneig(iat) = nneig(iat) + 1
            nneig(i)   = nneig(i)   + 1
            neighl(nneig(iat),iat) = i 
            neighl(nneig(i)  ,i  ) = iat 
          endif ! distance check  
          i = list(i) 
        enddo
        if(ANY(nneig(:) > lotfvar%nneigx))  then 
          write(message,'(a,i8,a)')&
            'UPD_LIS CLASSIC: max no. of neighbours: ',lotfvar%nneigx,' is too small'
          MSG_ERROR(message)
        endif
      enddo

      iat = list(iat) 
    enddo
  enddo
 end subroutine upd_lis0
 !!***



 !!****f* eval_lotf/tuneparms
 !! NAME
 !! tuneparms
 !!
 !! FUNCTION
 !!
 !! INPUTS
 !!  tau0(3,natom)=atomic positions 
 !!  rcf2_int=square of the cut off radius for this interaction 
 !!
 !! OUTPUT
 !!  tfit_int(3,nbondex)=logical that determines which parameters
 !!  will be optimized ( if true the parameter is optimized)  
 !!
 !! NOTES
 !!  simple SW version (eventually modify for Tersoff or others             
 !!  all the bonds considered here are already 100% in the fitting zone.    
 !!  Here we just need to eliminate bonds or triplets that are too long..
 !! 
 !! PARENTS
!!      m_lotf
!!
 !! CHILDREN
!!      dist_pbc,wrtout
!!
 !! SOURCE
 subroutine tuneparms(tau0,tfit_int,rcf2_int)

  use defs_param_lotf,only : lotfvar
  use bond_lotf,only : ibn_tot,ibnd_mat,nbondex
  USE pbc_lotf,only : dist_pbc
  implicit none

  !Arguments ------------------------
  real(dp),intent(in) :: rcf2_int
  real(dp),intent(in) :: tau0(3,lotfvar%natom)
  logical,intent(out) :: tfit_int(3,nbondex)
  !Local variables ------------------------------
  integer :: il_fit_int,ibn,iat,i,il_fit_par,ntot
  real(dp) :: r,r12,R1(3)
  character(len=500) :: message
  integer :: nbnds(lotfvar%natom),npract(lotfvar%natom)

! *************************************************************************

  tfit_int = .false.
  nbnds(:) = 0
  npract(:) = 0

  il_fit_int = 0
  il_fit_par = 0
  do ibn = 1,ibn_tot
    iat = ibnd_mat(1,ibn)
    i   = ibnd_mat(2,ibn)
    call dist_pbc(tau0(:,i),tau0(:,iat),r,R1)
    r12 = sqrt(r) 
    if(r < rcf2_int) then
      il_fit_int = il_fit_int + 1
      nbnds(iat) = nbnds(iat) + 1 
      nbnds(i)   = nbnds(i) + 1 

      tfit_int(:2,ibn) = .true.
      il_fit_par = il_fit_par + 2
      npract(iat) = npract(iat) + 2
      npract(i)   = npract(i) + 2
    endif
  enddo

  ntot = sum(npract)

  do ibn = 1,ibn_tot
    iat = ibnd_mat(1,ibn)
    i   = ibnd_mat(2,ibn)
    if( (npract(iat) > 50).and.(npract(i) > 50) )then  
      MSG_ERROR('LOTF: NPRACT 50 (A) ') 
    endif
  enddo

  write(message,'(2(a,i8,a))')&
    &    ' Total Active Bonds:          ',il_fit_int,ch10,&
    &    ' Total Active Bondparms:      ',il_fit_par,ch10
  call wrtout(std_out,message,'COLL')


 end subroutine tuneparms

end module eval_lotf
!!***
