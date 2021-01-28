!!****m* ABINIT/lotfpath
!! NAME
!! lotfpath
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

module LOTFPATH

 use defs_basis
 use m_errors
 use m_abicore

 use defs_param_lotf, only : lotfvar

 implicit none
 public

 !--Lotf variables used in non-lotf procedure (ex.moldyn)
 real(dp),allocatable,dimension(:,:) :: alpha,alpha_in,alpha_end

 real(dp),private :: epotlotf
 integer,private,allocatable,dimension(:)   :: nneig,nneig_old
 integer,private,allocatable,dimension(:,:) :: neighl,neighl_old
 real(dp),private,allocatable,dimension(:,:) :: xcartfit,velfit,fcartfit
 character(len=500),private :: message

 public  :: init_lotf
 public  :: end_lotf
 public  :: fitclus
 public  :: upd_lis
 public  :: force0
 public  :: intparms
 public  :: lotf_extrapolation
 public  :: vel_to_gauss
 public  :: force_to_vel
 public  :: vel_rescale
 public  :: extrapolation_loop

 private :: alpha_update

contains
!!***

!!****f* lotfpath/init_lotf
!! NAME
!! init_lotf
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_pred_lotf
!!
!! CHILDREN
!!      force0,force_to_vel,vel_to_gauss,wrtout
!!
!! SOURCE

 subroutine init_lotf(itime,natom,acell,rprimd,xcart)
  ! natom : number of atoms, abinit notation
  ! acell : cell length, abinit notation

  use glue_lotf,only  : glue_init
  use pbc_lotf,only   : pbc_init
  use eval_lotf,only  : upd_lis0
  use work_var_lotf,only : work_var_set,smallfit,cutoff_init

  use bond_lotf,only : ibn_tot,ibn_tots,nbondex,ibn_tot2,nfitmax,&
&                  bond_fit_set,bond_tafit_init,&
&                  bond_atom_init,     &
&                  bond_matrix_alloc,  &
&                  bond_matrix_set
  integer,intent(in) :: itime
  integer,intent(in) :: natom
  real(dp),intent(in) :: acell(3),rprimd(3,3)
  real(dp),intent(out) :: xcart(3,natom)


  integer :: nfitdum

! *************************************************************************

  ! I should modify units : LOTF uses atomic units
  !                         ABINIT uses  ???

  ! !-----------------------------------------
  ! Transfered in gstate
  ! call lotfvar_init(natom,&
  !   &               2, & !--version: set type of MD algo
  !   &           itime, & !--nstart: initial step
  !   &           nitex,& !--nitex: number of LOTF steps
  !   &               40,& !--nneigx: roughly decide the number of neighbours
  !   &               5, & !--classic: stick with the adaptable Glue model (rough version)
  !   &               1,1& !--me,nproc : disabled parallel LOTF
  !   )

  !-----------------------------------------
  !--Prepare LOTF variables used in moldyn.f90

  !--initialize potential energy :
   epotlotf = zero

   ABI_MALLOC(fcartfit,(3,natom))
   ABI_MALLOC(xcartfit,(3,natom))
   ABI_MALLOC(velfit,(3,natom))

   ABI_MALLOC(neighl,(lotfvar%nneigx,natom))
   ABI_MALLOC(nneig,(lotfvar%nneigx))
   ABI_MALLOC(neighl_old,(lotfvar%nneigx,natom))
   ABI_MALLOC(nneig_old,(lotfvar%nneigx))
   neighl = 0
   nneig = 0
   neighl_old = 0
   nneig_old = 0

  !--set work variables of LOTF
   call work_var_set()

  !--control on lotfvar%classic
   if(lotfvar%classic/=5 .and. lotfvar%classic/=6) then
     write(message,'(3a,3f12.6,a)')&
&     'LOTF: INIT_LIST: wrong value for lotfvar%classic = ',&
&     lotfvar%classic,ch10,&
&     'change lotfvar%classic 5 or 6 '
     ABI_ERROR(message)
   end if

  !--Init cell and pbc_lotf
   call pbc_init(rprimd)

  !--Initializes Glue, here we suppose we have atomic units.
   call Glue_INIT()

  !--Initializes cut-off radius
   call cutoff_init()

  !--last argument is to force the search for neighbors in upd_lis0
   call upd_lis0(xcart,neighl,nneig,itime)

  ! SMALL_FIT FINDS A RESTRICTED REGION WHERE THE FIT WILL BE
  ! PERFOMED. dimensionS WILL BE SET ACCORDING THE NUMBER OF ATOMS
  ! IN THIS SMALL REGION
  ! if niter=lotfvar%n0 it will set some dimensions for later use

  !--set tafit
   call bond_tafit_init(lotfvar%natom)
   call SMALLFIT(xcart,nfitdum)

  !--set nfitmax,ifit,nfit
   call bond_fit_set(lotfvar%natom,nfitdum)

  !--Initialize nbondex,neighl,neeig
   call bond_atom_init(lotfvar%nneigx,nneig,neighl)

  !---------------------------------------------
  !--Initialize/creates arrays needed by updlis :

  !--Allocate bond matrix
   call bond_matrix_alloc(lotfvar%natom,lotfvar%nneigx)

  !--updates bond matrix, associates the bond to the atom neighlists
   call bond_matrix_set(nneig,neighl)

  !--initialize bond parms alpha
   ABI_MALLOC(alpha,(3,nbondex))
   ABI_MALLOC(alpha_in,(3,nbondex))
   ABI_MALLOC(alpha_end,(3,nbondex))
   alpha = zero
   alpha_in = zero
   alpha_end = zero

   write(message,'(2a,i6,a,i8,2a,i8,2a,i8,i8)')ch10,&
&   'ITERATION:  ',itime,&
&   ' NO. OF BONDS BETWEEN FITTED ATOMS: ', ibn_tot,ch10,&
&   'TOTAL N.OF BOUNDS : ', ibn_tots,ch10,&
&   'BORDER BOUNDS (&max) : ', ibn_tot2,6*nfitmax
   call wrtout(std_out,message,'COLL')

 end subroutine init_lotf
 !!***

!!****f* lothpath/end_lotf
!! NAME
!! end_lotf
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!      force0,force_to_vel,vel_to_gauss,wrtout
!!
!! SOURCE

 subroutine end_LOTF()

  use work_var_lotf
  use bond_lotf

  !--init_lotf
! *************************************************************************
   ABI_FREE(alpha)
   ABI_FREE(alpha_in)
   ABI_FREE(alpha_end)
   ABI_FREE(nneig)
   ABI_FREE(neighl)
   ABI_FREE(fcartfit)
   ABI_FREE(xcartfit)
   ABI_FREE(velfit)

  !--deallocate LOTF internal variables
   call work_var_dealloc()
  !--deallocate LOTF bond variables
   call bond_dealloc()

 end subroutine end_LOTF
 !!***


!!****f* lothpath/fitclus
!! NAME
!! fitclus
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_pred_lotf
!!
!! CHILDREN
!!      force0,force_to_vel,vel_to_gauss,wrtout
!!
!! SOURCE

 subroutine fitclus(tfor,forc_in,tau0,alpha_tr,nqmin,nqmax)

  ! tfor : true if the forces are given in input
  !        false if fitclus calculates the forces
  ! forc_in(3,natom) : input forces to be fitted
  ! tau0(3,natom) : atomic positions
  ! alpha_tr(3,nbondex) : two body potential parameters
  ! neighl(lotfvar%nneigx,natom) : list of neighbours
  ! nneig(natom) : number of neighbours
  ! nqmin : lowest  index for the quantum atoms (nqmin = 1 ) !!!!NOT USED!!!
  ! nqmax : maximum index for the quantum atoms ( nqmax = natom)!!!!NOT USED!!!
  ! niter  : iteration number (itime)

  use work_var_lotf,only : rcut,ifixed
  use bond_lotf,only : nbondex,ifit,ibn_tots,nfitmax,imat,tafit,&
&                  ibnd_mat,ibn_tot,nfit
  use tools_lotf,only : dlvsum
  use glue_lotf, only : dphi,calc_coord_new_d,eval_u_n,eval_Upp_n
  use eval_lotf,only : tuneparms,phi_n_calc,calc_coord2,eval_forces_U_n,eval_forces_U_n_2,eval_force_devs_new_d

  !--Evaluates "real" forces from external source, and computes
  ! the istantaneous best fit of the current model.
  implicit none
  !Arguments ------------------------------------
  logical,intent(in) :: tfor
  integer,intent(in) :: nqmin, nqmax
  real(dp),intent(in) :: forc_in(3,lotfvar%natom)
  real(dp) :: tau0(3,lotfvar%natom)
  real(dp),intent(out) :: alpha_tr(3,nbondex)
  !Local ----------------------------------------
  real(dp), allocatable :: fspare(:,:),alpha_fce(:,:,:,:)
  real(dp)              :: rcut_fit_int, rcf2_int
  integer  :: i,j,k, iprecmass, iat
  integer  :: jat
  real(dp)  :: dcost,dcost_rms,dcost_old
  integer :: ibn,nitmax,icheck,nwarn,ibn_count,ifiat
  integer :: ifi
  real(dp)  :: epotlotf_dum,f_dum,dampfact,dcost_rm0,dcost_rm1
  real(dp)  :: d_dcost,dcost_check,dummy,vel,force
  real(dp)  :: dtm,toeva,dtest
  integer :: n,iem

  logical    :: tfit(3,nbondex),t_2nd(nbondex)
  real(dp)   :: fact(nfitmax),fact2(nfitmax),ffit(3,nfitmax)

  real(dp), allocatable :: alpha_dum(:,:),alpha_old(:,:)
  real(dp)   :: dcost_dalpha(3,nbondex),fcold(3,nbondex)

  real(dp) :: forc0_dum(3,0:nfitmax)
  integer  ::  iwrdri
  real(dp) :: dt_par(3),dmaspars(3)
  real(dp),parameter :: prec_lotf = 1.0d-5

  !--Glue variables
  integer :: istride, imin,imax,kat
  integer :: nlist(0:lotfvar%nneigx)
  integer :: neig2(lotfvar%nneigx),nlist2(lotfvar%nneigx,lotfvar%nneigx)
  real(dp) :: tauv(3,lotfvar%nneigx),forcv(3,0:lotfvar%nneigx)
  real(dp) :: stress(3,3)
  real(dp) :: coordatom(lotfvar%natom),up_list(lotfvar%natom),upp_list(lotfvar%natom),forc_dum2(3)
  real(dp) :: alpha_avg(3), alpha_disp(3)
  real(dp) :: tauv2(3,lotfvar%nneigx,lotfvar%nneigx)
  real(dp) :: rho_p_sum(3,nfitmax)

! *************************************************************************

  ! #####################  end DECLARATION #########################
   call wrtout(std_out,' LOTF : FITCLUS','COLL')

   iwrdri = 0

   ABI_MALLOC(fspare,(3,nbondex))
   ABI_MALLOC(alpha_fce,(3,2,nbondex,1))

  !--(0,-2) initialises a few parameters
   rcut_fit_int = rcut
   rcf2_int     = rcut_fit_int * rcut_fit_int

   dmaspars = (/ 0.01_dp, 0.1_dp,  0.1_dp /)
   dt_par   = (/ 0.4_dp, 0.4_dp, 0.0_dp /)

  !--Not needed in this case... see if we can safely get rid of it
  !     dt_ang          = zero
  ! (0,-1) sets (optional: preconditions) masses of bond parameters

   iprecmass = 0

  !--Iitialize variables  :
   t_2nd = .true.
   tfit = .false.
   ffit = zero
   dcost_dalpha = zero

  !--set constant array for speed-up : fact(i)
   do i = 1,nfit
     iat = ifit(i)
     fact2(i) =        real(ifixed(iat),dp)
     fact(i)  = half * real(ifixed(iat),dp)
   end do

  !--(0,0) decides which parameters will be fitted
   call tuneparms(tau0,tfit,rcf2_int)

  !--(1,0) computes true forces: needs the external "expensive" routine.
   fcold = zero

   if (tfor) then ! tfor = true <-> WE CALCULATE THE FORCES HERE !
     do i =1, nfit
       iat = ifit(i)
      !write(97,*) iat,forc_in(1,iat)
       ffit(:,i) = forc_in(:,iat)
     end do
   else
     ABI_ERROR('LOTF : HERE WE SHOULD HAVE THE FORCES ALREADY !! ')
   end if ! TFOR

  !--THE REST OF THE ROUTINE IS THE FIT
   write(message,'(2(a,i8,a))')&
   ' nbondex',  nbondex,ch10,' ibn_tots : ', ibn_tots,ch10
   call wrtout(std_out,message,'COLL')


   ABI_MALLOC(alpha_dum,(3,nbondex))
   ABI_MALLOC(alpha_old,(3,nbondex))
   alpha_dum = zero
   alpha_old = zero

  !--(1,2) initialises the parameter set Alpha_dum, which IS modified
  !       in the optimisation section
   forall(ibn=1:ibn_tots)  alpha_dum(:,ibn) = (/ dphi,one, zero /)

  !--(2,0) computes derivatives of the ionic forces wrt parameters
  !       at fixed ionic positions.

  !--Debug.......................
   nitmax = 1000000

   icheck     = 0
   nwarn      = 0
   dcost      = zero
   dcost_old  = zero
   alpha_fce  = zero
   dcost_rms  = one

  !###################################
  !--MAIN MINIMISATION LOOP BEGINS
  !###################################
   main_minimization: do while(dcost_rms >  prec_lotf)
    ! if(mod(icheck,20)==1.and.(lotfvar%me==1)) write(*,*)  icheck&
    !   &,dcost,dcost_rms , prec_lotf,dcost_rms > prec_lotf

    !---------
    !       derivatives of two-body forces w.r.t. parameters
    !---------

    !--I choose to duplicate the code, can be made better :
     if(lotfvar%classic==5) then
       forc0_dum = zero
      !--Testing........................................
       alpha_fce = zero
      !................................................
      !--Test_5
       dcost_dalpha = zero
      ! End Test_5

      !--I cannot run only over fit pairs because of the term U(n)
      ! so I run over fit atoms and ALL its neighbours twice, as I did in Force_Glue...

      !--Also: I only consider 1 parameter as variable, alpha(1,:)
      ! therefore I don't use t_2nd nor tfit(:,ibn_count)

      !--(I) Pair potential (and its forces), coordination,  U, and  Up
       nlist(:) = 0
       stress(:,:) = zero
       coordatom(:) = zero
       up_list(:) = zero

      !--Parallel version
       istride = lotfvar%natom/lotfvar%nproc
       if(istride*lotfvar%nproc < lotfvar%natom) istride = istride + 1
       imin = (lotfvar%me-1)*istride + 1
       imax = lotfvar%me*istride
       if(imax > lotfvar%natom) imax = lotfvar%natom

      !--First run over pairs: pair potential, its forces, coordination, U, total E, and Up
       overatom1 : do iat = imin,imax
         nlist(0) = iat
        !   write(*,*) 'number of neighbours',nneig(iat)

         do j = 1,nneig(iat)
           i = neighl(j,iat)
           tauv(:,j) = tau0(:,i)
           nlist(j) = i
         end do

         if (tafit(iat)) then !--Only for fit atoms

          !--Pair potential (and its forces) & coordination(i)
           call phi_n_calc(alpha_dum,nneig(iat),nlist,tau0(:,iat),&
&           tauv,epotlotf_dum,forcv,coordatom(iat),alpha_fce)

           forc0_dum(:,imat(iat)) = forc0_dum(:,imat(iat)) + forcv(:,0)

           do j = 1,nneig(iat)
             i = neighl(j,iat)
             forc0_dum(:,imat(i)) = forc0_dum(:,imat(i)) + forcv(:,j)
           end do

           call eval_U_n(coordatom(iat),epotlotf_dum,up_list(iat))
         else   !--non fit atoms

          !--Evaluate coordination for ALL the non fit atoms ... SHOULD BE IMPROVED to include only neighbours of neighbours of fit atoms...
           call calc_coord2(nneig(iat),tau0(1,iat),tauv,coordatom(iat))

           call eval_U_n(coordatom(iat),epotlotf_dum,up_list(iat))

         end if !--fit / non fit atoms
       end do overatom1

       call dlvsum(lotfvar%me-1,lotfvar%nproc,alpha_fce(1,1,1,1),6*nbondex)
       call dlvsum(lotfvar%me-1,lotfvar%nproc,up_list(1),lotfvar%natom)


      !--(II) Glue forces (coming only from the embedding function)
       do iat = imin,imax
         if (tafit(iat)) then
           nlist(0) = iat
           do j = 1,nneig(iat)
             i = neighl(j,iat)
             tauv(:,j) = tau0(:,i)
             nlist(j) = i
           end do

           call eval_forces_U_n(nneig(iat),nlist,tau0(1,iat),tauv,up_list,&
&           forc_dum2)
           forc0_dum(:,imat(iat)) = forc0_dum(:,imat(iat)) + forc_dum2(:)
         end if
       end do
     elseif (lotfvar%classic==6) then !-----------------------------------------------------------
       forc0_dum = zero
      !--Testing
       alpha_fce = zero

      !--Test_5
       dcost_dalpha = zero
      !--End Test_5

      !--(I) Pair potential (and its forces), coordination,  U, and  Up
       nlist(:) = 0
       stress(:,:) = zero
       coordatom(:) = zero
       up_list(:) = zero
       upp_list(:) = zero

      !--Parallel version
       istride = lotfvar%natom/lotfvar%nproc
       if(istride*lotfvar%nproc < lotfvar%natom) istride = istride + 1
       imin = (lotfvar%me-1)*istride + 1
       imax = lotfvar%me*istride
       if(imax > lotfvar%natom) imax = lotfvar%natom

      !--First run over pairs: pair potential, its forces, coordination, U, total E, and Up
       overatom2 : do iat = imin,imax
         nlist(0) = iat
         do j = 1,nneig(iat)
           i = neighl(j,iat)
           tauv(:,j) = tau0(:,i)
           nlist(j) = i
         end do

         if(tafit(iat)) then ! Only for fit atoms

          !--Pair potential (and its forces) & coordination(i)
           call phi_n_calc(alpha_dum,nneig(iat),nlist,tau0(:,iat),&
&           tauv,epotlotf_dum,forcv,coordatom(iat),alpha_fce) ! COORDATOM is OK (lotfvar%classic dependence)
           forc0_dum(:,imat(iat)) = forc0_dum(:,imat(iat)) + forcv(:,0)
           do j = 1,nneig(iat)
             i = neighl(j,iat)
             forc0_dum(:,imat(i)) = forc0_dum(:,imat(i)) + forcv(:,j)
           end do

          !--Up(coord(i)) & Upp(coord(i))
           call eval_Upp_n(coordatom(iat),up_list(iat),upp_list(iat)) ! Up and Upp
         else   ! non fit atoms

          !--Evaluate coordination for ALL the non fit atoms ... CAN BE IMPROVED to include only neighbours of neighbours of fit atoms...
           call calc_coord2(nneig(iat),tau0(:,iat),tauv,&
&           coordatom(iat))                   ! OK (non-fit)

           call eval_Upp_n(coordatom(iat),up_list(iat),upp_list(iat)) ! OK (non-fit)

         end if ! fit / non fit atoms
       end do overatom2

       call dlvsum(lotfvar%me-1,lotfvar%nproc,up_list(1),lotfvar%natom)
       call dlvsum(lotfvar%me-1,lotfvar%nproc,upp_list(1),lotfvar%natom)
       call dlvsum(lotfvar%me-1,lotfvar%nproc,alpha_fce(1,1,1,1),6*nbondex)

      !--(II) Glue forces (coming only from the embedding function)
       do iat = imin,imax
         if (tafit(iat)) then
           nlist(0) = iat
           do j = 1,nneig(iat)
             jat = neighl(j,iat)
             tauv(:,j) = tau0(:,jat)
             nlist(j) = jat
           end do

           call eval_forces_U_n_2(alpha_dum,nneig(iat),nlist,&
&           tau0(:,iat),tauv,up_list,&
&           rho_p_sum(:,imat(iat)),&
&           forc0_dum(:,imat(iat)))
         end if
       end do

       call dlvsum(lotfvar%me-1,lotfvar%nproc,rho_p_sum,3*nfitmax)
       call dlvsum(lotfvar%me-1,lotfvar%nproc,forc0_dum,(nfitmax+1)*3)

      !--(III) Derivatives of the forces (inside DCOST_DALPHA)
       neig2(:) = 0
       nlist2(:,:) = 0
       tauv(:,:) = zero
       tauv2(:,:,:) = zero

       do iat = imin,imax
         if (tafit(iat)) then
           nlist(0) = iat
           do j = 1,nneig(iat)
             jat = neighl(j,iat)
             tauv(:,j) = tau0(:,jat)
             nlist(j) = jat
             do k=1,nneig(jat)
               neig2(j)=nneig(jat)
               kat = neighl(k,jat)
               nlist2(k,j)=kat
               tauv2(:,k,j)=tau0(:,kat)
             end do
           end do

           call eval_force_devs_new_d(alpha_dum,nneig(iat),nlist,neig2,nlist2,&
&           tau0(:,iat),tauv,tauv2,up_list,upp_list,&
&           fact2,ffit,forc0_dum,rho_p_sum,dcost_dalpha)
         end if
       end do

      ! In paralell all this is going to cause trouble...........................
      !      call dlvsum(lotfvar%me-1,lotfvar%nproc,dcost_dalpha,3*nbondex)
      !......... check later more operations over dcost_dalpha...................

      !--------------------------------------------------------------------------------------------------

     end if ! lotfvar%classic

    ! if(lotfvar%classic /= 5) then
    !   call dlvsum(lotfvar%me-1,lotfvar%nproc,forc0_dum,(nfitmax+1)*3)
    ! end if

     if (lotfvar%classic /= 6)  then
       call dlvsum(lotfvar%me-1,lotfvar%nproc,forc0_dum,(nfitmax+1)*3)
     end if

    !--Check1
    !--(1) evaluate cost function: quadratic deviation from true forces
     dcost = zero
     do i = 1,nfit
       dcost = dcost+fact(i)*sum((forc0_dum(:,i)-ffit(:,i))**2,dim=1)
     end do
     dcost_rms = (sqrt((two*dcost)/nfit))*two*13.6058d0/0.529177d0

    !--minimisation is achived
     if(dcost_rms < prec_lotf) exit

    !--(2) evaluate its derivative downhill w.r.t two body parms
     do ibn_count = 1,ibn_tot

       if( (lotfvar%me == (mod(ibn_count-1,lotfvar%nproc)+1)) ) then

         iat   = ibnd_mat(1,ibn_count)
         i     = ibnd_mat(2,ibn_count)
         ifiat = imat(iat)
         ifi   = imat(i)

         do n=1,2
           f_dum = zero
           if(tfit(n,ibn_count)) then
            !--note : reduction here does not give good results so the
            !loop in k has to be left
             do k=1,3
               f_dum = f_dum +  fact2(ifi)&
&               *  (forc0_dum(k,ifi)   - ffit(k,ifi))&
&               *  alpha_fce(k,n,ibn_count,1)

               f_dum = f_dum -  fact2(ifiat)&
&               *  (forc0_dum(k,ifiat) - ffit(k,ifiat))&
&               *  alpha_fce(k,n,ibn_count,1)
             end do
           end if

          ! Debug
           dcost_dalpha(n,ibn_count) = dcost_dalpha(n,ibn_count) + f_dum
         end do
       end if
     end do


    !--parameter OPTIMIZATION STEP --------------------------- start
     icheck = icheck + 1

     if(icheck ==1) alpha_old(:,:ibn_tot) = alpha_dum(:,:ibn_tot)

    !-----------------------------------------
    !--damped dynamics

    !--(1) damping factor
     if(icheck > 40) then
       dampfact = 0.6
       if(icheck > 200)  dampfact = 0.6
     else
       dampfact = 0.6
     end if

     if(dcost  >  dcost_old)  nwarn = nwarn + 1

     dcost_rm0 = (sqrt((two*dcost_old)/nfit))*two*13.6058d0/0.529177d0
     dcost_rm1 = (sqrt((two*dcost)/nfit))*two*13.6058d0/0.529177d0
     d_dcost = dcost_rm1 - dcost_rm0
     dcost_old = dcost


     if(icheck < nitmax) then
       if(icheck==nitmax/2) dcost_check = dcost_rms

      !--II BODY
       bodyii: do  ibn_count = 1, ibn_tot
         if( (lotfvar%me == (mod(ibn_count-1,lotfvar%nproc)+1)) ) then
           if(.not.t_2nd(ibn_count)) cycle

          !--(2) updating strategy
           do n=1,3
             dummy = alpha_dum(n,ibn_count)
             if(tfit(n,ibn_count)) then
               vel   = (alpha_dum(n,ibn_count)-alpha_old(n,ibn_count))
               force = dcost_dalpha(n,ibn_count)

               fcold(n,ibn_count) = force

              !--(3) Verlet step
               alpha_dum(n,ibn_count) = alpha_dum(n,ibn_count)&
&               + dampfact * vel + dt_par(n) * dt_par(n) * force/dmaspars(n)

              !--(4) bound control for the parameters
               if(alpha_dum(1,ibn_count)  >  6.1d0) then
                 alpha_dum(1,ibn_count) = 6.1d0
               elseif(alpha_dum(1,ibn_count)  <  two ) then
                 alpha_dum(1,ibn_count) = two
                !ABI_ERROR('LOTF: Alpha1 reaches 2 au... Too small value!')
               end if

             end if   ! tfit(n)

             alpha_old(n,ibn_count) = dummy
           end do

         else
           alpha_old(:,ibn_count) = zero
           alpha_dum(:,ibn_count) = zero
         end if

       end do bodyii

       call dlvsum(lotfvar%me-1,lotfvar%nproc,alpha_dum,3*ibn_tot)

       dcost_rms=(sqrt((two*dcost)/nfit))*two*13.6058d0/0.529177d0

      !--Minimization is not achived
      !if(dcost_rms >  prec_lotf)  cycle
     end if   ! if icheck < ##
    !!qui!!
    !--Minimization is not achived (dcost_rms > prec_lotf)
   end do main_minimization

   if(dcost_rms >  prec_lotf) then
     ABI_ERROR('LOTF: ACHTUNG: REQD.TOLERANCE NOT ACHIEVED IN THE FIT')
   end if

   iwrdri = iwrdri + 1

  !--To have a look at what is going on with alpha...
   if(lotfvar%me==1) then

     alpha_avg(:) = sum(alpha_dum(:,:ibn_tots),dim=2)
     alpha_avg(:) = (one/ibn_tots)*alpha_avg(:)

     alpha_disp(:) = zero
     do i=1,ibn_tots
       alpha_disp(:) = alpha_disp(:) + (alpha_dum(:,i) - alpha_avg(:))**2
     end do
     alpha_disp(:) = sqrt((one/ibn_tots)*alpha_disp(:))

     write(message,'(2(2a,2f12.8))')ch10,&
&     'Alpha_avg =',(alpha_avg(i),i=1,2),ch10,&
&     'Alpha_disp =',(alpha_disp(i),i=1,2)
     call wrtout(std_out,message,'COLL')


     dtm = zero

     write(message,'(a,2f12.8,a)')'FIT.PREC. : ',dcost, dcost_rms,' EV/A '
     call wrtout(std_out,message,'COLL')

     toeva = 2.d0*13.6058d0/0.529177d0
     do i = 1,nfit
       dtest = zero
       do k=1,3
         dtest = dtest+sqrt(fact2(i)*(forc0_dum(k,i)-ffit(k,i))**2)*toeva
       end do
       if(dtest > dtm) then
         dtm = dtest
         iem = ifit(i)
       end if
     end do
     write(message,'(a,f12.8,a,i8)')'Max Error : ',dtm, ' on atom : ',iem
     call wrtout(std_out,message,'COLL')

   end if

  !--PARAMETER OPTIMIZATION STEP ---------------------------end
  !if (lotfvar%classic /= 5 .AND. lotfvar%classic /= 6) then
  !  call dlvsum(lotfvar%me-1,lotfvar%nproc,alpha_dum,3*ibn_tot)
  !end if

  !--Final result for alpha............................
   if (lotfvar%me==1) then
     call wrtout(std_out,'alpha(1)   alpha(2)','COLL')
     do ibn_count = 1,ibn_tots
       write(message,'(2f12.8)') (alpha_dum(n,ibn_count),n=1,2)
       call wrtout(std_out,message,'COLL')
     end do
   end if

  !--prepares "true" updated parameters
   alpha_tr(:,:ibn_tots) = alpha_dum(:,:ibn_tots)

   ABI_FREE(alpha_fce)
   ABI_FREE(alpha_dum)
   ABI_FREE(alpha_old)
  !-----------------------------------------------------------

 end subroutine fitclus
 !!***


!!****f* lothpath/upd_lis
!! NAME
!! upd_lis
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!      force0,force_to_vel,vel_to_gauss,wrtout
!!
!! SOURCE

 subroutine  upd_lis(tau0,niter,alpha_dum)
  ! lotfvar%nneigx : max number of neighbours
  ! tau0(3,natom) : atomic positions
  ! neighl(lotfvar%nneigx,natom) : list of neighbours
  ! neighl_old(lotfvar%nneigx,natom) : old list of neighbours
  ! nneig_old(natom) : old number of neighbours
  ! niter  : iteration number (itime)

  !  updates the neighbour list for each atom, and updates the active
  !  parameter lists

  USE GLUE_LOTF,only : dphi
  use work_var_lotf,only : smallfit
  use bond_lotf,only : ibn_tot,ibn_tot2,ibn_tots,ibmat,ibnd_mat,&
&                  nfitmax,imat,bond_matrix_set,bond_compute,&
&                  bond_fit_set,nbondex
  use eval_lotf,only : upd_lis0
  implicit none

  !Arguments ------------------------------------
  integer,intent(in)   :: niter
  real(dp),intent(inout) :: tau0(3,lotfvar%natom)
  real(dp),intent(inout) :: alpha_dum(3,nbondex)
  !Local ----------------------------------------
  integer :: j,ibn_count
  integer :: ifo,jato,jb_old
  integer :: iat,jat,nfitdum
  integer :: jbo(nbondex)
  !take care nsurf_at*(nbond_at)/2

! *************************************************************************

  !--INITIALIZATION AND DECISION :
   jbo(:) = 0

  !--neighbours update
   call upd_lis0(tau0,neighl,nneig,niter)

  !--set tafit and compute nfitdum
   call SMALLFIT(tau0,nfitdum)

  !--control nfitmax and set ifit,nfit
   call bond_fit_set(lotfvar%natom,nfitdum)

  !--Updates bond matrix, associates the bond to the atom neighlists
   call bond_compute(nneig,neighl)

  !--CREATES THE REORDERING ARRAY (ibmat)
  !  ONLY variational BONDS ARE CONSIDERED IN HISTORY
   do ibn_count = 1, ibn_tot  ! only on fitted atoms pairs
    !--finds to which atoms it corresponds
     iat = ibnd_mat(1,ibn_count)
     jat = ibnd_mat(2,ibn_count)

     ifo = imat(iat)    !--imat finds the old atomic 'fitted number'
     if(iat > jat) ABI_ERROR('UPDLIS 177')

    !--Set to 0, finds to which old bond (if any) these two correspond
     jb_old = 0               !--atom jat is a new neighbour of atom iat
     do j =1,nneig_old(iat)
       jato = neighl_old(j,iat)
       if(jato==jat) then
         jb_old = ibmat(j,ifo)  !--atom jat was already a neighbour of atom iat
         exit
       end if
     end do
     jbo(ibn_count) = jb_old  !--index array for reordering
   end do

  !--updates bond matrix, associates the bond to the atom neighlists
   call bond_matrix_set(nneig,neighl)


   write(message,'(2a,2(a,i8,a),(a,2i8,a),a,i8,a)') &
&   ' LOTF : UPD_LIS',ch10,&
&   ' ITERATION:  ',niter,ch10,&
&   ' NO. OF BONDS IN ACTIVE ZONE : ',ibn_tot,ch10,&
&   ' BORDER BOUNDS (&max) : ', ibn_tot2,6*nfitmax ,ch10,&
&   ' TOTAL N.OF BOUNDS : ', ibn_tots,ch10
   call wrtout(std_out,message,'COLL')

  !--updates old neighbour lists
   do iat =1,lotfvar%natom
     nneig_old(iat) = nneig(iat)
     neighl_old(:nneig(iat),iat) = neighl(:nneig(iat),iat)
    !write(6,*) iat,nneig(iat),(neighl(j,iat),j=1,nneig(iat))
   end do

  !--updates parameter list
   if(niter /= lotfvar%n0) then
     call alpha_update(dphi,jbo,alpha_dum)
   end if

   if (ibn_tots > nbondex) then
     write(message,'(2a,2(a,i8))') 'LOTF: ibn_tots > nbondex  ! ',ch10,&
&     'UPDLIS  stop : IBNTOTS = ',ibn_tots,' NBONDEX = ',nbondex
     ABI_ERROR(message)
   end if


 end subroutine upd_lis
!!***


!!****f* lothpath/force0
!! NAME
!! force0
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!      force0,force_to_vel,vel_to_gauss,wrtout
!!
!! SOURCE

 subroutine force0(tau0,forc0,alpha_tr,amass)
  ! tau0(3,natom) : atomic positions (input)
  ! forc0(3,natom) : atomic forces(output)
  ! alpha_tr(3,nbondex) : two body potential parameters
  ! neighl(lotfvar%nneigx,natom) : list of neighbours
  ! nneig(natom) : number of neighbours
  ! epotlotf : potential energy (parameter dependent)
  use glue_lotf,only : eval_U_n
  use tools_lotf,only : dlvsum
  use bond_lotf,only : ibn_tot,nbondex,tafit
  use eval_lotf, only : phi_n_calc,calc_coord2,eval_forces_u_n
  implicit none

  !Arguments ------------------------------------
  real(dp),intent(in) :: amass(lotfvar%natom)
  real(dp),intent(in) :: tau0(3,lotfvar%natom)
  real(dp),intent(out) :: forc0(3,lotfvar%natom)
  real(dp),intent(in) :: alpha_tr(3,nbondex)
  !Local ----------------------------------------
  integer :: i, j
  integer :: istride,ibmin,ibmax,iat
  logical  :: tcalc(3)
  real(dp) :: epotlotf_dum
  real(dp) :: epotlotf_2
  real(dp) :: stress(3,3)

  !--Glue variables
  integer ::  nlist(0:lotfvar%nneigx)
  real(dp) ::  tauv(3,lotfvar%nneigx),forcv(3,0:lotfvar%nneigx)
  real(dp) ::  coordatom(lotfvar%natom),alpha_fdum_v(3,2,nbondex,1)
  real(dp) ::  up_list(lotfvar%natom),epotlotf_dum2,epotlotf_dum3
  real(dp) ::  forc_dum2(3)

! *************************************************************************

   epotlotf_2 = zero
   tcalc(:) = .false.

   forc0(:,:) = zero    ! MODifIED FOR ABINITttime 23/07/2008

  !--parallel version
   istride = ibn_tot/lotfvar%nproc
   if(istride*lotfvar%nproc < ibn_tot) istride = istride + 1
   ibmin = (lotfvar%me-1) * istride + 1
   ibmax =  lotfvar%me * istride
   if(ibmax > ibn_tot) ibmax = ibn_tot


  !--All but glue forces
   nlist(:) = 0
   epotlotf_dum = zero
   stress(:,:) = zero
   coordatom(:) = zero
   up_list(:) = zero
   epotlotf_dum2 = zero

   do iat=1,lotfvar%natom

     nlist(0) = iat
     do j = 1,nneig(iat)
       i = neighl(j,iat)
       tauv(:,j) = tau0(:,i)
       nlist(j) = i
     end do

     if (tafit(iat)) then
       call phi_n_calc(alpha_tr,nneig(iat),nlist,tau0(1,iat),&
       tauv,epotlotf_dum,forcv,coordatom(iat),alpha_fdum_v)
      !--PAIR energy: Fit atoms and NEIGHBOURS --> OK
       epotlotf_2 = epotlotf_2 + epotlotf_dum

       forc0(:,iat) = forc0(:,iat) + forcv(:,0)
       do j = 1,nneig(iat)
         i = neighl(j,iat)
         forc0(:,i) = forc0(:,i) + forcv(:,j)
       end do

       call eval_U_n(coordatom(iat),epotlotf_dum2,up_list(iat))
       epotlotf_2 = epotlotf_2 + epotlotf_dum2 ! GLUE energy: Fit atoms ONLY --> OK

     else

      !--Evaluate coordination and up_list for ALL the non fit atoms ... CAN BE IMPROVED to include only neighbours of neighbours of fit atoms...
      ! For example: run over neigbours, in case one of them is FIT atom, proceed, otherwise, do nothing
       call calc_coord2(nneig(iat),tau0(1,iat),tauv,coordatom(iat))

       call eval_U_n(coordatom(iat),epotlotf_dum3,up_list(iat))
      !--We do NOT accumulate epotlotf_dum3
     end if

   end do

  !--Glue forces
   do iat=1,lotfvar%natom
     if (tafit(iat)) then
       nlist(0) = iat
       do j = 1,nneig(iat)
         i = neighl(j,iat)
         tauv(:,j) = tau0(:,i)
         nlist(j) = i
       end do
       call eval_forces_U_n(nneig(iat),nlist,tau0(1,iat),tauv,up_list,forc_dum2)

       forc0(:,iat) = forc0(:,iat) + forc_dum2(:)
     end if

    !--Renomalization of the force (to get acceleration)
     forc0(:,iat) = forc0(:,iat)/amass(iat)
   end do

   epotlotf = epotlotf + epotlotf_2
  !--ends parallelisation

 end subroutine force0
 !!***

!!****f* lothpath/intparms
!! NAME
!! intparms
!!
!! FUNCTION
!!
!! INPUTS
!! PARENTS
!!      m_pred_lotf
!!
!! CHILDREN
!!      force0,force_to_vel,vel_to_gauss,wrtout
!!
!! SOURCE

 subroutine intparms(itime)
  use bond_lotf,only : ibn_tot
  use tools_lotf,only : pinterp,pinterp_nolinear
  implicit none
  !Arguments ------------------------------------
  integer,intent(in) :: itime!,nitex
  !Local ----------------------------------------
  integer :: ibn,interp_type,nitdu

! *************************************************************************

   nitdu = mod(itime-lotfvar%n0,lotfvar%nitex) + 1

  !--Reset alpha
   alpha(:,:) = zero

  !--Select here the type of interpolation....................
   interp_type = 1
  !   interp_type = 2
  !--II body
   if(interp_type==1) then
     do ibn = 1, ibn_tot
       call pinterp(alpha_in(:,ibn),alpha_end(:,ibn),alpha(:,ibn),3,lotfvar%nitex,nitdu)
     end do
   end if
 end subroutine intparms
 !!***

!!****f* lothpath/alpha_update
!! NAME
!! alpha_update
!!
!! FUNCTION
!!  updates parameter list
!! INPUTS
!!  dphi=parameter to reinatialise bond parameters
!!  jbo= index array for reordering
!!
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!      force0,force_to_vel,vel_to_gauss,wrtout
!!
 !! SOURCE
 subroutine alpha_update(dphi,jbo,alpha_dum)
  use bond_lotf,only : ibn_tot
  implicit none

  !Arguments ------------------------
  real(dp),intent(in) :: dphi
  integer,intent(in) :: jbo(:)
  real(dp),intent(inout) :: alpha_dum(:,:)
  !Local ---------------------------
  integer :: ibn,jb_old
  real(dp) :: alphawork(3,ibn_tot)

! *************************************************************************

   do ibn = 1,ibn_tot
     jb_old = jbo(ibn)
     if(jb_old /= 0) then
      !--swaps old parms
       alphawork(:,ibn) = alpha_dum(:,jb_old)
     else
      !--or reinitialise bond parms
       alphawork(:,ibn) = (/ dphi, one, zero /)
     end if
   end do

   alpha_dum(:,:ibn_tot) = alphawork(:,:)

 end subroutine alpha_update
 !!***

!!****f* lothpath/lotf_extrapolation
!! NAME
!! lotf_extrapolation
!!
!! FUNCTION
!!  return true if mod(itime,nitex) == 0
!! INPUTS
!! CHILDREN
!!
!! SOURCE

 function lotf_extrapolation(itime)

  implicit none

  !Arguments ------------------------------------
  integer,intent(in) :: itime
  logical :: lotf_extrapolation

! *************************************************************************

   if(itime == 1) then
     lotf_extrapolation = .true.
     return
   end if
   if(lotfvar%nitex-lotfvar%n0==0) then
     lotf_extrapolation = .true.
   else
     lotf_extrapolation = mod(itime-lotfvar%n0,lotfvar%nitex)==0
   end if
   return
 end function lotf_extrapolation
 !!***

!!****f* lothpath/force_to_vel
!! NAME
!! force_to_vel
!!
!! FUNCTION
!!  Compute velocity starting from : vel_in,v2gauss and forces
!! INPUTS
!!  v2gauss=gauss factor (2*kinetic energy)
!!  dtion=time step for Molecular Dynamics
!!  amass(natom)=masse of the ions
!!  vel_in(3,natom)=initial velocity of ions
!!  fcart(3,natom)=force on ions
!!
!! OUTPUTS
!!  vel_out(3,natom)=new velocity of ions
!!
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!      force0,force_to_vel,vel_to_gauss,wrtout
!!
!! SOURCE

 subroutine force_to_vel(v2gauss,dtion,amass,vel_in,fcart,vel_out)

  implicit none

  !Arguments ------------------------------------
  real(dp),intent(in) :: v2gauss,dtion
  real(dp),intent(in) :: amass(:)
  real(dp),intent(in) :: vel_in(:,:)
  real(dp),intent(in) :: fcart(:,:)
  real(dp),intent(out) :: vel_out(:,:)
  !Local ---------------------------
  integer :: iatom,idim
  real(dp) :: a,b,sqb,as,s1,s2,s,scdot

! *************************************************************************

  !--Compute a and b (4.13 de Ref.1)
   a = zero
   b = zero
   do iatom=1,size(amass)
     do idim=1,3
       a=a+fcart(idim,iatom)*vel_in(idim,iatom)*amass(iatom)
       b=b+fcart(idim,iatom)*fcart(idim,iatom)*amass(iatom)
     end do
   end do

   a= a/v2gauss
   b= b/v2gauss

  !--Campute  s et scdot
   sqb = sqrt(b)
   as = sqb*dtion/2.
   s1 = cosh(as)
   s2 = sinh(as)
   s = a*(s1-1.)/b+s2/sqb
   scdot = a*s2/sqb+s1
   vel_out(:,:) = (vel_in(:,:)+fcart(:,:)*s)/scdot
 end subroutine force_to_vel
 !!***

!!****f* lothpath/vel_to_gauss
!! NAME
!! vel_to_gauss
!!
!! FUNCTION
!!  Compute gauss factor sum(v**2*m) the double of kinetic energy
!!  If present vtest compute also sum(v)/(3*sum(m))
!! INPUTS
!!  vel_in(3,natom)=velocity to use
!!  amass(natom)=masse of the ions
!!
!! OUTPUT
!!  v2gauss=2*kinetic energy
!!  vtest=pick velocity
!!
!! PARENTS
!!      m_lotf,m_pred_lotf
!!
!! CHILDREN
!!      force0,force_to_vel,vel_to_gauss,wrtout
!!
!! SOURCE

 subroutine vel_to_gauss(vel_in,amass,v2gauss,vtest)

  implicit none

  !Arguments ------------------------------------
  real(dp),intent(out) :: v2gauss
  real(dp),intent(out),optional :: vtest
  real(dp),intent(in) :: vel_in(:,:)
  real(dp),intent(in) :: amass(:)
  !Local ---------------------------
  integer :: iatom,idim

! *************************************************************************

   if(present(vtest)) then
     v2gauss = zero
     vtest = zero
     do iatom=1,size(amass)
       do idim=1,3
         v2gauss = v2gauss+vel_in(idim,iatom)*vel_in(idim,iatom)*amass(iatom)
         vtest = vtest+vel_in(idim,iatom)
       end do
     end do
     vtest = vtest/(3._dp*size(amass))
   else
    !--Recompute v2gauss with the rescaled velocities
     v2gauss = zero
     do iatom=1,size(amass)
       do idim=1,3
         v2gauss = v2gauss+vel_in(idim,iatom)*vel_in(idim,iatom)*amass(iatom)
       end do
     end do
   end if
 end subroutine vel_to_gauss
 !!***

!!****f* lothpath/vel_rescale
!! NAME
!! vel_rescale
!!
!! FUNCTION
!!  Starting from the velocity, it recompute the velocities in the
!!  center of mass. Then compute the gauss factor and renormalises
!!  the velocities with respect the gauss distribution. Then is
!!  recompute the gauss factor
!!
!! INPUTS
!!  mditemp=temperature of ions
!!  amass(natom)=masse of the ions
!!
!! OUTPUTS
!!  v2gauss=2*kinetic energy
!!
!! SIDE EFFECTS
!!  vel(3,natom)=velocity of ions
!!
!! PARENTS
!!      m_lotf,m_pred_lotf
!!
!! CHILDREN
!!      force0,force_to_vel,vel_to_gauss,wrtout
!!
!! SOURCE

 subroutine vel_rescale(mditemp,vel,amass,v2gauss)

  real(dp),intent(in) :: mditemp
  real(dp),intent(out) :: v2gauss
  real(dp),intent(in) :: amass(:)
  real(dp),intent(inout) :: vel(:,:)

  integer :: idim,natom
  real(dp) :: mass_tot,momentum,rescale_vel,vtest,sigma2

! *************************************************************************

   natom = size(amass)

  !--Get rid of center-of-mass velocity
   mass_tot = sum(amass(:))
   do idim=1,3
     momentum = sum(amass(:)*vel(idim,:))
     vel(idim,:) = vel(idim,:)-momentum/mass_tot
   end do

  !--Recompute v2gauss
   call vel_to_gauss(vel,amass,v2gauss,vtest)

  !--Now rescale the velocities to give the exact temperature
   rescale_vel = sqrt(three*natom*kb_HaK*mditemp/v2gauss)
   vel(:,:) = vel(:,:)*rescale_vel

  !--Recompute v2gauss
   call vel_to_gauss(vel,amass,v2gauss)

  !--Compute the variance and print
   sigma2 = (v2gauss/(3._dp*natom)-amass(1)*vtest**2)/kb_HaK

   write(message, '(a,d12.5,a,D12.5)' )&
&   ' --- LOTF STEP : Effective temperature',&
&   v2gauss/(3*natom*kb_HaK),' From variance', sigma2
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')


 end subroutine vel_rescale
!!***

!!****f* lothpath/extrapolation_loop
!! NAME
!! extrapolation_loop
!!
!! FUNCTION
!!  Compute the LOTF extrapolation:
!!  Starting from xcart_0,vel_0,fcart_0, it computes in first the new
!!  bond, then the forces on the atoms.
!!  At this point if compute the position,speed and forces on atoms
!!  upto the step nitex:  xcart_nitex, vel_nitex, fcart_nitex
!!  This will be used in a SCF to compute the forces in the final
!!  point of the LOTF approximation.
!!  In output the positions (extrapoled) in the step ntitex.
!!
!! INPUTS
!!  itime=numeber of MD step
!!  mditemp=temperature of ions
!!  dtion=time step for Molecular Dynamics
!!  amass(natom)=masse of the ions
!!  xcart(3,natom)=position of the ions initial
!!  vel(3,natom)=speed of the ions initial
!!
!! OUT
!!  xcart_next(3,natom)=positions of the ions final step
!!  the following variable are used as work variables:
!!  vel_nexthalf(3,natom)=velocity of ions in next half step
!!
!! PARENTS
!!      m_pred_lotf
!!
!! CHILDREN
!!      force0,force_to_vel,vel_to_gauss,wrtout
!!
!! SOURCE

 subroutine extrapolation_loop(itime,mditemp,dtion,amass,&
&                           xcart,vel,xcart_next,vel_nexthalf)

  implicit none

  !Arguments ------------------------------------
  integer,intent(in) :: itime
  real(dp),intent(in) :: mditemp,dtion
  real(dp),intent(in) :: amass(:)
  real(dp),intent(in) :: vel(:,:)
  real(dp),intent(in) :: xcart(:,:)
  real(dp),intent(out) :: vel_nexthalf(:,:)
  real(dp),intent(out) :: xcart_next(:,:)
  !Local ---------------------------
  integer :: itex
  real(dp) :: v2gauss

! *************************************************************************

  !--inital values for position and forces
   xcartfit(:,:) = xcart(:,:)
   velfit(:,:) = vel(:,:)

  !--Update bond variables and alpha
   call upd_lis(xcartfit,itime,alpha_in)

  !print *,'b-end uplis',itime,sum(sum(alpha,dim=2)),sum(sum(alpha_in,dim=2)),sum(sum(alpha_end,dim=2))

  !--Store alpha in alpha_in
   alpha(:,:) = alpha_in(:,:)

  !--compute fcartfit (fcartfit is the accelleration)
   call force0(xcartfit,fcartfit,alpha,amass)

   do_itex : do itex = 1, lotfvar%nitex
    !--Compute rescaled velfit (center of mass speed and gaussian)
     call vel_rescale(mditemp,velfit,amass,v2gauss)

    ! start  verletvel here
     if(itime==1.AND.itex==1) then  ! check this
       vel_nexthalf(:,:) = velfit(:,:)
       xcart_next(:,:) = xcartfit(:,:)
     else
      !--Computation of vel_nexthalf (4.16 de Ref.1)
       call force_to_vel(v2gauss,dtion,amass,velfit,fcartfit,vel_nexthalf)

      !--Computation of the next positions
       xcart_next(:,:) = xcartfit(:,:)+vel_nexthalf(:,:)*dtion


      !--compute fcartfit and alpha (fcartfit is the accelleration)
       call force0(xcart_next,fcartfit,alpha,amass)


      !--Computation of vel(:,:) at the next positions
      !--Computation of v2gauss
       call vel_to_gauss(vel_nexthalf,amass,v2gauss)

      !--Compute velocity from force
       call force_to_vel(v2gauss,dtion,amass,vel_nexthalf,fcartfit,velfit)
     end if
     xcartfit = xcart_next

   end do do_itex

 end subroutine extrapolation_loop
 !!***

!!****f* lothpath/lotf_interpolation
!! NAME
!! lotf_interpolation
!!
!! FUNCTION
!!  Compute the LOTF interpolation:
!!
!! INPUTS
!!  itime=numeber of MD step
!!  dtion=time step for Molecular Dynamics
!!  amass(natom)=masse of the ions
!!  xcart(3,natom)=position of the ions initial
!!
!! OUTPUTS
!!  xcart_next(3,natom)=positions of the ions final step
!!  vel_nexthalf(3,natom)=velocity of ions in next half step
!!
!! SIDE EFFECTS
!!  v2gauss=gauss factor (twice the kinetic energy) (initial and final)
!!  vel(3,natom)=velocity of ions
!!  fcart_m(3,natom)=forces on ions in
!!
!! PARENTS
!!      m_pred_lotf
!!
!! CHILDREN
!!      force0,force_to_vel,vel_to_gauss,wrtout
!!
!! SOURCE

 subroutine lotf_interpolation(itime,dtion,v2gauss,amass,xcart,vel,&
&                           fcart_m,xcart_next,vel_nexthalf)

  implicit none

  !Arguments ------------------------------------
  integer,intent(in) :: itime
  real(dp),intent(in) :: dtion
  real(dp),intent(inout) :: v2gauss
  real(dp),intent(in) :: amass(:)
  real(dp),intent(in) :: xcart(:,:)
  real(dp),intent(inout) :: vel(:,:)
  real(dp),intent(inout) :: fcart_m(:,:)
  real(dp),intent(out) :: vel_nexthalf(:,:)
  real(dp),intent(out) :: xcart_next(:,:)

! *************************************************************************

   write(message,'(a,i8)') ' ---LOTF interpolation: itime=',itime
   call wrtout(std_out,message,'COLL')

  !--VERLETVEL (interpolation)
   if(itime==lotfvar%n0) then  ! check this
    !--itime=0 fist step
     vel_nexthalf(:,:) = vel(:,:)
     xcart_next(:,:)  = xcart(:,:)
   else
    !--Computation of vel_nexthalf (4.16 de Ref.1)
     call force_to_vel(v2gauss,dtion,amass,vel,fcart_m,vel_nexthalf)

    !--Computation of the next positions
     xcart_next(:,:) = xcart(:,:)+vel_nexthalf(:,:)*dtion

    !--compute fcart_m and alpha (fcart_m is the accelleration)
     call force0(xcart_next,fcart_m,alpha,amass)

    !--Computation of vel(:,:) at the next positions
     call vel_to_gauss(vel_nexthalf,amass,v2gauss)
     call force_to_vel(v2gauss,dtion,amass,vel_nexthalf,fcart_m,vel)

   end if

 end subroutine lotf_interpolation

end module LOTFPATH
!!***
