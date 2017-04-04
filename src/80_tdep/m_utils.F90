#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_utils
  
  use defs_basis
  use m_latt,        only : Lattice_Variables_type, inbox
  use m_readwrite,   only : Input_Variables_type
  use m_sym,         only : Symetries_Variables_type, SearchMatR_1at, SearchMatR_2at

 implicit none

 public :: write_d1NN
 public :: calc_MoorePenrose
 public :: MatchIdeal2Average
 public :: model

contains

!====================================================================================================
 subroutine write_d1NN(distance_average,FromIdeal2Average,InVar,bond_ref,Sym,xcart,xcart_ideal)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_d1NN'
!End of the abilint section

  implicit none 

  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  integer,intent(in) :: FromIdeal2Average(InVar%natom)
  integer,intent(in) :: bond_ref(InVar%natom,InVar%natom,3)
  double precision, intent(in) :: xcart_ideal(3,InVar%natom)
  double precision, intent(in) :: distance_average(InVar%natom,InVar%natom,4)
  double precision, intent(in) :: xcart(3,InVar%natom,InVar%nstep)

  integer :: ii,jj,kk,iatom,iatcell,jatom,katom,latom,istep,isym
  integer :: total_point,index_proj_1d
  integer :: datagrid3d(3),index_proj(3)
  integer :: index_ideal(3)
  integer, allocatable :: proba_3d(:,:,:), proba_1d(:)
  double precision :: ratio,dist,dist_ideal,dist_mean,dist_average
  double precision :: tmp(3),tmp_ideal(3),tmp_mean(3),tmp_ave(3),tmp2(3),tmp3(3)
  double precision :: maxdist(3)

! Define quantities
  maxdist(:)=10
  datagrid3d(:)=50
  allocate(proba_3d(datagrid3d(1),datagrid3d(2),datagrid3d(3))); proba_3d(:,:,:)=0
  allocate(proba_1d(int(2*InVar%Rcut*datagrid3d(1)))); proba_1d(:)=0

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '####### Compute the distance for the first nearest neighbour ################'
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,'(a,i4)') ' Write 1NN density of probability of the shell number',InVar%d1NN
  write(InVar%stdout,'(a)')    ' See the Xcrysden "d1NN.xsf", XYZ "d1NN.xyz" and 1D "d1NN.1d" files'
  write(InVar%stdout,'(a)') ' '


! Compute the density of probability in 3D and 1D
! Search the couple of atoms linked using a symetry relation
  open(unit=15,file='d1NN.xyz')
  dist_average=zero
  dist_mean   =zero
  tmp_mean(:)=zero
  tmp_ave (:)=zero
  do iatcell=1,InVar%natom_unitcell
    do jatom=1,InVar%natom
      if (jatom==iatcell.or.bond_ref(jatom,iatcell,3).ne.InVar%d1NN) cycle
      tmp_ideal(:)=xcart_ideal(:,jatom)-xcart_ideal(:,iatcell)
      dist_ideal=dsqrt(sum(tmp_ideal(1:3)*tmp_ideal(1:3)))

      do katom=1,InVar%natom
        do latom=1,InVar%natom
          if ((bond_ref(katom,latom,1)==iatcell).and.(bond_ref(katom,latom,2)==jatom)) then
            isym=abs(Sym%matR(katom,latom))
            do istep=1,InVar%nstep_remain

!             Verify that the distance is "inside the sphere of NN" 
!             The distance could be badly computed (higher) due to PBC
              tmp2(:)=zero
              do jj=1,3
                tmp2(jj)=InVar%xred(jj,FromIdeal2Average(latom),InVar%nstep1+istep-1)&
&                       -InVar%xred(jj,FromIdeal2Average(katom),InVar%nstep1+istep-1)
              end do
              call inbox(tmp2,1,1d-3)

!             Transform the rotated bond in cartesian coordinates
              tmp3(:)=zero
              call DGEMV('T',3,3,1.d0,InVar%rprimd_MD(:,:),3,tmp2(:),1,0.d0,tmp3(:),1)
              dist=dsqrt(sum(tmp3(1:3)*tmp3(1:3)))
              ratio=dist/dist_ideal
              if ((ratio.gt.2.0).or.(ratio.lt.0.5)) then
                write(InVar%stdout,'(a)') ' WARNING: a distance is out of range.  Something wrong!'
                cycle
              end if

!             Transform the "present bond" toward the "reference bond"
              tmp(:)=zero
              do ii=1,3
                do jj=1,3
!                 Transpose MatR
                  tmp(ii)=tmp(ii)+Sym%matR_ref(jj,ii,isym,1)*tmp3(jj)
                end do
              end do 
              if (Sym%matR(katom,latom).lt.0) tmp(:)=-tmp(:)
!FB              write(InVar%stdout,'(a,3f15.10)') ' Vector after transformation:',tmp(:)

              do ii=1,3
                if (abs(tmp(ii)).gt.maxdist(ii)) then
                  write(InVar%stdout,'(a)') ' STOP: the coordinate is larger than the grid!'
                  write(InVar%stdout,'(i2,x,2f15.10)') ii,tmp(ii),maxdist(ii)
                  stop
                end if  
              end do  

!             Compute the mean value of the bond distance
              tmp_mean(:)=tmp_mean(:)+tmp(:)
!USE MATREF BEFORE tmp_ave (:)=tmp_ave (:)+distance_average(katom,latom,2:4)
              dist_mean   =dist_mean   +dist
              dist_average=dist_average+distance_average(katom,latom,1)

!             Compute the density of probability in 1D and 3D
              index_proj(:)=int((tmp(:)/maxdist(:)+1)/2.d0*datagrid3d(:))
              index_proj_1d=int(dist*datagrid3d(1))
              do ii=1,3
                if (index_proj(ii).le.0) then
                  write(InVar%stdout,'(a)') ' STOP: the index is lower than 0'
                  stop
                end if
              end do
              if ((index_proj_1d.le.0).or.(minval(index_proj(:)).le.0)&
&             .or.(index_proj_1d.gt.int(2*InVar%Rcut*datagrid3d(1)))&
&             .or.(index_proj(1).gt.datagrid3d(1))&
&             .or.(index_proj(2).gt.datagrid3d(2))&
&             .or.(index_proj(3).gt.datagrid3d(3))) then
                write(InVar%stdout,'(a)') ' BUG: the code stop'
                write(InVar%stdout,'(a,i8)') ' index_proj_1d=',index_proj_1d
                write(InVar%stdout,'(a,i8)') ' int(2*InVar%Rcut*datagrid3d(1))=',int(2*InVar%Rcut*datagrid3d(1))
                write(InVar%stdout,'(a,3i8)') ' index_proj(:)=',index_proj(:)
                write(InVar%stdout,'(a,3i8)') ' datagrid3d(:)=',datagrid3d(:)
                stop
              end if
!             Plot only a slice of the total distribution
              if (InVar%bravais(2).eq.-1) then
                if (abs((tmp(1)-tmp_ideal(1))+(tmp(2)-tmp_ideal(2))).lt.0.1d0) then
                  write(15,'(3(f15.10))')tmp(:)
                  proba_1d(index_proj_1d)=proba_1d(index_proj_1d)+1
                end if  
              else if (InVar%bravais(2).eq.-3) then
                if (abs( tmp(3)-tmp_ideal(3))                       .lt.0.1d0) then
                  write(15,'(3(f15.10))')tmp(:)
                  proba_1d(index_proj_1d)=proba_1d(index_proj_1d)+1
                end if  
              else
                write(15,'(3(f15.10))')tmp(:)
                proba_1d(index_proj_1d)=proba_1d(index_proj_1d)+1
              end if
              proba_3d(index_proj(1),index_proj(2),index_proj(3))=proba_3d(index_proj(1),index_proj(2),index_proj(3))+1
            end do !istep
          end if  !(iatcell,jatom)<->(katom,latom)
        end do !latom
      end do !katom
    end do !jatom
  end do !iatcell
  close(15)

! Write the density of probability in 3D
  open(unit=16,file='d1NN.xsf')
  write(16,'(a)') 'DIM-GROUP'
  write(16,'(a)') '3  1'
  write(16,'(a)') 'PRIMVEC'
  do ii=1,3
    write(16,'(3f15.10)') InVar%rprimd_MD(ii,:)
  end do  
  write(16,'(a)') 'PRIMCOORD'
  write(16,'(i4,a4)') InVar%natom,' 1'
  do iatom=1,InVar%natom
    if (iatom.le.InVar%natom_unitcell) then
      write(16,'(a4,3f15.10)') ' 22',xcart_ideal(:,iatom)
    else  
      write(16,'(a4,3f15.10)') ' 13',xcart_ideal(:,iatom)
    end if  
  end do  
  write(16,'(a)') 'ATOMS'
  do iatom=1,InVar%natom
    if (iatom.le.InVar%natom_unitcell) then
      write(16,'(a4,3f15.10)') ' 22',xcart_ideal(:,iatom)
    else  
      write(16,'(a4,3f15.10)') ' 13',xcart_ideal(:,iatom)
    end if  
  end do  
  write(16,'(a)') 'BEGIN_BLOCK_DATAGRID3D'
  write(16,'(a)') 'datagrids'
  write(16,'(a)') 'DATAGRID_3D_DENSITY'
  write(16,'(3i4)') datagrid3d(1),datagrid3d(2),datagrid3d(3)
  write(16,'(3f15.10)') -maxdist(1),-maxdist(2),-maxdist(3)
  write(16,'(3f15.10)') 2*maxdist(1),0.d0,0.d0
  write(16,'(3f15.10)') 0.d0,2*maxdist(2),0.d0
  write(16,'(3f15.10)') 0.d0,0.d0,2*maxdist(3)

  total_point=zero
  do kk=1,datagrid3d(3)
    do jj=1,datagrid3d(2)
      do ii=1,datagrid3d(1)
        write(16,'(i8)') proba_3d(ii,jj,kk)
        total_point=total_point + proba_3d(ii,jj,kk)
      end do
    end do
  end do
  write(16,'(a)') ' ' 
  write(16,'(a)') 'END_DATAGRID_3D'
  write(16,'(a)') 'END_BLOCK_DATAGRID3D'
  close(16)

! Write the density of probability in 1D
  open(unit=17,file='d1NN.1d')
  do ii=1,int(2*InVar%Rcut*datagrid3d(1))
    write(17,'(f15.10,i8)') dfloat(ii)/dfloat(datagrid3d(1)),proba_1d(ii)
  end do
  close(17)
  deallocate(proba_3d,proba_1d)

  write(InVar%stdout,'(a,f15.10)') ' The ideal   distance is=',dist_ideal 
  write(InVar%stdout,'(a,f15.10)') ' The mean    distance is=',dist_mean/dfloat(total_point) 
  write(InVar%stdout,'(a,f15.10)') ' The distance between mean positions is=',dist_average/dfloat(total_point) 
  write(InVar%stdout,'(a,3f15.10)') ' The coordinates of the ideal vector are=',tmp_ideal(:)
!NO MEANING AFTER A TRANSFORMATION x<->y  write(InVar%stdout,'(a,3f15.10)') ' The coordinates of the mean  vector are=',tmp_mean(:)/dfloat(total_point)
!FB  write(InVar%stdout,'(a,3f15.10)') ' The coordinates of the vector between mean postitions are=',tmp_ave(:)/dfloat(total_point)

 end subroutine write_d1NN

!====================================================================================================
 subroutine calc_MoorePenrose(fcartij,fcoeff,InVar,ntotcoeff,Phij_coeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_MoorePenrose'
!End of the abilint section

  implicit none 

  type(Input_Variables_type),intent(in) :: InVar
  integer,intent(in) :: ntotcoeff
  double precision, intent(in)  :: fcoeff(3*InVar%natom*InVar%nstep_remain,ntotcoeff)
  double precision, intent(in)  :: fcartij(3*InVar%natom*InVar%nstep_remain)
  double precision, intent(out)  :: Phij_coeff(ntotcoeff,1)
  
  integer :: ii,LWORK,INFO
  integer, allocatable :: IWORK(:)
  double precision, allocatable :: WORK(:),pseudo_inverse(:,:),sigma(:),matU(:,:)
  double precision, allocatable :: pseudo_sigma(:,:),transmatV(:,:),tmp1(:,:),fcartij_tmp(:,:)


  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '########################## Compute the pseudo-inverse #######################'
  write(InVar%stdout,*) '#############################################################################'

  allocate(pseudo_inverse(ntotcoeff,3*InVar%natom*InVar%nstep_remain)) ; pseudo_inverse(:,:)=0.d0
  allocate(sigma(ntotcoeff)) ; sigma(:)=0.d0 
  allocate(matU(3*InVar%natom*InVar%nstep_remain,ntotcoeff)) ; matU(:,:)=0.d0
  allocate(transmatV(ntotcoeff,ntotcoeff)) ; transmatV(:,:)=0.d0
  LWORK=3*(ntotcoeff)**2+max(3*InVar%natom*InVar%nstep_remain,4*(ntotcoeff)**2+4*(ntotcoeff))
  allocate(WORK(LWORK),IWORK(8*(ntotcoeff))) ; WORK(:)=0.d0
  call DGESDD('S',3*InVar%natom*InVar%nstep_remain,ntotcoeff,fcoeff,3*InVar%natom*InVar%nstep_remain,sigma,matU,3*InVar%natom*InVar%nstep_remain,transmatV,ntotcoeff,WORK,LWORK,IWORK,INFO) 
  deallocate(WORK,IWORK)

  allocate(pseudo_sigma(ntotcoeff,ntotcoeff)) ; pseudo_sigma(:,:)=0.d0
  sigma(:)=1.d0/sigma(:)
  write(InVar%stdout,*) 'The eigenvalues are:'
  do ii=1,ntotcoeff
    if (sigma(ii).lt.1.d8) then
      pseudo_sigma(ii,ii)=sigma(ii)
    end if  
    write(InVar%stdout,'(x,i4,x,f15.10)') ii,sigma(ii)
  end do
  write(InVar%stdout,'(a,x,f15.10)')'  condition number=',maxval(sigma(:))/minval(sigma(:))
  deallocate(sigma)
  
  allocate(tmp1(ntotcoeff,3*InVar%natom*InVar%nstep_remain)) ; tmp1(:,:)=0.d0
  call DGEMM('N','T',ntotcoeff,3*InVar%natom*InVar%nstep_remain,ntotcoeff,1.d0,pseudo_sigma,ntotcoeff,matU,3*InVar%natom*InVar%nstep_remain,1.d0,tmp1,ntotcoeff)
  call DGEMM('T','N',ntotcoeff,3*InVar%natom*InVar%nstep_remain,ntotcoeff,1.d0,transmatV,ntotcoeff,tmp1,ntotcoeff,1.d0,pseudo_inverse,ntotcoeff)
  deallocate(tmp1)
  deallocate(matU,transmatV,pseudo_sigma)
  allocate(fcartij_tmp(3*InVar%natom*InVar%nstep_remain,1)); fcartij_tmp(:,:)=zero  
  fcartij_tmp(:,1)=fcartij(:)
! NOTE, we have to solve F_ij = -\sum_j \Phi_ij u_j, so we add a minus sign to the pseudo_inverse  
  call DGEMM('N','N',ntotcoeff,1,3*InVar%natom*InVar%nstep_remain,1.d0,-pseudo_inverse,ntotcoeff,fcartij_tmp,3*InVar%natom*InVar%nstep_remain,1.d0,Phij_coeff,ntotcoeff)
  deallocate(fcartij_tmp)
  deallocate(pseudo_inverse)


 end subroutine calc_MoorePenrose

!====================================================================================================
 subroutine MatchIdeal2Average(bond_ref,distance,fcartij,InVar,Lattice,nshell,&
&                              Rlatt_cart,Rlatt4dos,Sym,ucart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'MatchIdeal2Average'
!End of the abilint section

  implicit none 

  type(Input_Variables_type),intent(inout) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Symetries_Variables_type),intent(inout) :: Sym
  integer, intent(out)  :: nshell
  integer, intent(out)  :: bond_ref(InVar%natom,InVar%natom,3)
  double precision, intent(out)  :: distance(InVar%natom,InVar%natom,4)
  double precision, intent(out)  :: fcartij(3*InVar%natom*InVar%nstep_remain)
  double precision, intent(out)  :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)
  double precision, intent(out)  :: Rlatt4dos (3,InVar%natom_unitcell,InVar%natom)
  double precision, intent(out)  :: ucart(3,InVar%natom,InVar%nstep_remain)
  
  integer :: ii,jj,kk,max_ijk,iatcell,jatcell,iatom,jatom,katom,latom,istep,jstep
  integer :: foo,foo2,atom_ref,ishell
  double precision :: tmp(3),tmp1(3),tmp2(3),Rlatt(3),xred_tmp(3),rprimd_MD_tmp(3,3)
  double precision, allocatable :: dist_unitcell(:,:,:),xcart_average(:,:)
  double precision, allocatable :: fcart_tmp(:,:,:),ucart_tmp(:,:,:)
  double precision, allocatable  :: xred_average(:,:)
  double precision, allocatable  :: xred_center(:,:)
  double precision, allocatable  :: Rlatt_red (:,:,:)
  double precision, allocatable  :: xred_ideal(:,:)
  double precision, allocatable  :: distance_average(:,:,:)
  integer, allocatable  :: FromIdeal2Average(:)
  double precision, allocatable  :: xcart(:,:,:)
  double precision, allocatable  :: xcart_ideal(:,:)
  logical :: ok,ok1

!==========================================================================================
!======== 1/ Determine ideal positions and distances ======================================
!==========================================================================================
!Check that atoms (defined in the input.in file) are set correctly in the unitcell
  do ii=1,3
    do iatcell=1,InVar%natom_unitcell
      if ((InVar%xred_unitcell(ii,iatcell).le.(-0.5)).or.(InVar%xred_unitcell(ii,iatcell).gt.(0.5))) then
        write(InVar%stdout,*) 'xred_unitcell='
        write(InVar%stdout,*)  InVar%xred_unitcell(:,1:InVar%natom_unitcell)
        write(InVar%stdout,*) 'Please put the atoms in the ]-0.5;0.5] range'
        stop
      endif
    end do
  end do

! Define the bigbox with ideal positions
  allocate(Rlatt_red (3,InVar%natom_unitcell,InVar%natom)); Rlatt_red (:,:,:)=0.d0
  allocate(xred_ideal(3,InVar%natom))                     ; xred_ideal(:,:)=0.d0
  max_ijk=40
  iatom=1
  do ii=-max_ijk,max_ijk
    do jj=-max_ijk,max_ijk
      do kk=-max_ijk,max_ijk
        do iatcell=1,InVar%natom_unitcell
          if (iatcell==1) ok=.false.
          Rlatt(1)=dfloat(ii-1)
          Rlatt(2)=dfloat(jj-1)
          Rlatt(3)=dfloat(kk-1)
!         Then compute the reduced positions
          tmp(:)=Rlatt(:)+InVar%xred_unitcell(:,iatcell)
          call DGEMV('T',3,3,1.d0,Lattice%multiplicitym1(:,:),3,tmp(:),1,0.d0,xred_tmp(:),1)

!         If the first atom of the pattern is in the [0;1[ range then keep all the
!         atoms of the pattern (even if the others are outside the box). Elsewhere, 
!         none are taken.
          if (iatcell==1) then
            if (minval(xred_tmp(:)).lt.0.d0.or.maxval(xred_tmp(:)).ge.(1.d0-1.d-12)) then
              cycle
            else
              ok=.true.
            end if  
          else 
            if (.not.ok) cycle
          end if
          if (iatom.gt.(InVar%natom+1)) then
            write(InVar%stdout,*) 'The number of atoms found in the bigbox exceeds natom' 
            stop
          end if  
          xred_ideal(:,iatom)=xred_tmp(:)
          call DGEMV('T',3,3,1.d0,Lattice%multiplicitym1(:,:),3,Rlatt(:),1,0.d0,Rlatt_red(:,1,iatom),1)
          iatom=iatom+1
        end do   
      end do
    end do
  end do

  if (iatom.lt.InVar%natom) then
    write(InVar%stdout,*) 'STOP : The number of atoms found in the bigbox is lower than natom' 
    stop
  end if  

! Compute the distances between ideal positions in the SUPERcell
  do katom=1,InVar%natom
    do latom=1,InVar%natom
      tmp(:)=xred_ideal(:,latom)-xred_ideal(:,katom)
      call inbox(tmp,1,1d-3)
      call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,tmp(:),1,0.d0,distance(katom,latom,2:4),1)
      do ii=1,3
!       Remove the rounding errors before writing (for non regression testing purposes)
        if (abs(distance(katom,latom,ii+1)).lt.tol8) distance(katom,latom,ii+1)=zero
        distance(katom,latom,1)=distance(katom,latom,1)+(distance(katom,latom,ii+1))**2
      end do
      distance(katom,latom,1)=distance(katom,latom,1)**0.5
    end do  
  end do  

! Compute the distances between ideal positions in the UNITcell
  allocate(dist_unitcell(InVar%natom_unitcell,InVar%natom_unitcell,3)); dist_unitcell(:,:,:)=zero
  do iatcell=1,InVar%natom_unitcell
    do jatcell=1,InVar%natom_unitcell
      tmp(:)=xred_ideal(:,jatcell)-xred_ideal(:,iatcell)
      call inbox(tmp,1,tol8)
      dist_unitcell(iatcell,jatcell,:)=tmp(:)
    end do
  end do  

!==========================================================================================
!======== 2/ Find the matching between the ideal and average ==============================
!========   (from the MD simulations) positions. ==========================================
!======== NOTE: - xred_center is used to find the matching with the ideal positions ======= 
!========       - xred_average is used to compute the displacements (from MD trajectories) 
!==========================================================================================
  allocate(xred_average(3,InVar%natom))             ; xred_average(:,:)=0.d0
  allocate(xred_center(3,InVar%natom))              ; xred_center(:,:)=0.d0
! Average positions from MD (on nstep_remain steps)
  do istep=1,InVar%nstep_remain
    do iatom=1,InVar%natom
      xred_average(:,iatom)=xred_average(:,iatom)+InVar%xred(:,iatom,InVar%nstep1+istep-1)
    end do
  end do
  xred_average(:,:)=xred_average(:,:)/dfloat(InVar%nstep_remain)

! Search the basis of atoms in the supercell
  ok=.true.
  xred_center(:,:)=xred_average(:,:)
  do iatom=1,InVar%natom
    foo2=0
    do jatom=1,InVar%natom
      tmp(:)=xred_center(:,jatom)-xred_center(:,iatom)
      call inbox(tmp,1,1d-3)
      iatcell=1
      do jatcell=1,InVar%natom_unitcell
        foo=0
        do ii=1,3
          if ((abs(tmp(ii)-dist_unitcell(iatcell,jatcell,ii)).le.InVar%tolmotif).and.&
&               InVar%typat(iatom).eq.InVar%typat_unitcell(iatcell).and.&
&               InVar%typat(jatom).eq.InVar%typat_unitcell(jatcell)) then
            foo=foo+1
          end if
        end do
        if (foo==3) then
          foo2=foo2+1
          exit
        end if
      end do
    end do
    if (foo2.eq.InVar%natom_unitcell) then
      atom_ref=iatom
      write(InVar%stdout,*) 'ATOM REF=',atom_ref
      ok=.false.
      exit
    else if (foo2.gt.InVar%natom_unitcell) then
      write(InVar%stdout,*) ' Something wrong: foo2=',foo2
      stop
    endif
  end do  
  if (ok) then
    open(unit=31,file='xred_average.xyz')
    do iatom=1,InVar%natom
      write(31,'(a,x,3(f10.6,x))') 'C',xred_center(:,iatom)
      write(31,'(a,x,3(f10.6,x))') 'I',xred_ideal (:,iatom)
    end do
    close(31)
    write(InVar%stdout,*) 'The basis of atoms written in input.in file does not appear in the MD trajectory'
    write(InVar%stdout,*) 'Perhaps, you can adjust the tolerance (tolmotif)'
    stop
  end if  
  deallocate(dist_unitcell)

! Modification of xred and Rlatt tabs
! For averaged quantities --> kk=1: xred_center, xred_average et xred 
! For ideal quantities    --> kk=2: Rlatt_red et xred_ideal 
  do kk=1,2
!   1/ The "atom_ref" (kk=1) or iatom=1 (kk=2) atom is put in (0.0;0.0;0.0)
    if (kk==1) then
      tmp(:)=xred_center(:,atom_ref)
    else if (kk==2) then
      tmp1(:)=xred_ideal(:,1)
      tmp2(:)=Rlatt_red(:,1,1)
    end if  
    do jatom=1,InVar%natom
      if (kk==1) xred_center(:,jatom)=xred_center(:,jatom)-tmp(:)
      if (kk==2) then
        xred_ideal(:,jatom)=  xred_ideal(:,jatom)  -tmp1(:)
        Rlatt_red (:,1,jatom)=Rlatt_red (:,1,jatom)-tmp2(:)
      end if        
    end do  
!   2/ All the atoms are put in the range [-0.5;0.5[ (use of PBC)
    do jatom=1,InVar%natom
      if (kk==1) then
        tmp(:)=xred_center(:,jatom)
        call inbox(tmp,1,InVar%tolinbox,xred_center (:,jatom))
        call inbox(tmp,1,InVar%tolinbox,xred_average(:,jatom))
        do istep=1,InVar%nstep_remain
          call inbox(tmp,1,InVar%tolinbox,InVar%xred(:,jatom,InVar%nstep1+istep-1))
        end do  
      else if (kk==2) then
        tmp(:)=xred_ideal(:,jatom)
        call inbox(tmp,1,tol8,xred_ideal(:,jatom))
        call inbox(tmp,1,tol8,Rlatt_red(:,1,jatom))
!FB        call inbox(Rlatt_red(:,1,jatom),1,tol8)
      end if  
    end do  
  end do  

! When the multiplicity equals 1 along one direction, there is some trouble
! To clean!!!!!!!
  do ii=1,3
    if ((InVar%multiplicity(ii,ii).eq.1).and.(InVar%multiplicity(ii,mod(ii  ,3)+1).eq.0)&
&                                 .and.(InVar%multiplicity(ii,mod(ii+1,3)+1).eq.0)) then
      Rlatt_red(ii,1,:)=0.d0
      write(InVar%stdout,*) 'WARNING: multiplicity=1 for ii=',ii
    end if
  end do

! Define Rlatt for all the atoms in the basis (Rlatt_red varies as a function of iatcell)
  if (InVar%natom_unitcell.gt.1) then
    do iatcell=2,InVar%natom_unitcell
      Rlatt_red(:,iatcell,:)=Rlatt_red(:,1,:)
    end do
  end if  
  do iatom=1,InVar%natom
    do iatcell=1,InVar%natom_unitcell
      tmp(:)=xred_ideal(:,iatom)-xred_ideal(:,iatcell)
      call inbox(tmp,1,tol8,Rlatt_red(:,iatcell,iatom))
    end do
  end do
  if (InVar%debug) then
    do iatcell=1,InVar%natom_unitcell
      write(InVar%stdout,*) 'For iatcell=',iatcell
      do jatom=1,InVar%natom
        write(InVar%stdout,'(a,i4,a,3(f16.10,x))') 'For jatom=',jatom,', Rlatt=',Rlatt_red(1:3,iatcell,jatom)
      end do  
    end do  
  end if

! Matching between Ideal and Average positions: xred_ideal and xred_center
! Then, write them in the xred_average.xyz file.
  allocate(FromIdeal2Average(InVar%natom))             ; FromIdeal2Average(:)=0
  open(unit=31,file='xred_average.xyz')
  write(31,'(i4)') InVar%natom*2
  write(31,'(i4)') 1
  do iatom=1,InVar%natom
    do jatom=1,InVar%natom
      ok =.true.
      ok1=.true.
      foo=0
      do ii=1,3
        if ((abs(xred_center(ii,iatom)-xred_ideal(ii,jatom)     ).le.InVar%tolmatch.and.(InVar%typat(iatom).eq.InVar%typat_unitcell(mod(jatom-1,InVar%natom_unitcell)+1)))) then
          foo=foo+1
        else if ((abs(xred_center(ii,iatom)-xred_ideal(ii,jatom)-1.d0).le.InVar%tolmatch.and.(InVar%typat(iatom).eq.InVar%typat_unitcell(mod(jatom-1,InVar%natom_unitcell)+1))).or.&
&           (abs(xred_center(ii,iatom)-xred_ideal(ii,jatom)+1.d0).le.InVar%tolmatch.and.(InVar%typat(iatom).eq.InVar%typat_unitcell(mod(jatom-1,InVar%natom_unitcell)+1)))) then
          foo=foo+1
          ok1=.false.
        endif
      end do
      if (foo==3.and.ok1) then
        FromIdeal2Average(jatom)=iatom
        ok=.false.
        exit
      else if (foo==3.and..not.ok1) then
        write(InVar%stdout,*) '  THE CODE STOPS'
        write(InVar%stdout,*) '  Some positions are outside the [-0.5;0.5[ range:'
        write(InVar%stdout,*) '  xred_center(:,',iatom,')='
        write(InVar%stdout,*) xred_center(:,iatom)
        write(InVar%stdout,*) '  xred_ideal(:,',jatom,')='
        write(InVar%stdout,*) xred_ideal(:,jatom)
        write(InVar%stdout,*) 'Perhaps, you can adjust the tolerance (tolinbox)'
        stop
      end if  
    end do  
    if (ok) then
      write(InVar%stdout,*) 'Problem to find the average position for iatom=',iatom
      write(InVar%stdout,*) '  Reasons:'
      write(InVar%stdout,*) '    1/ One atom jump to another equilibrium position'
      write(InVar%stdout,*) '    2/ The system is no more solid'
      write(InVar%stdout,*) '    3/ Perhaps, you can adjust the tolerance (tolmatch)'
      write(InVar%stdout,*) '  xred_center=',(xred_center(ii,iatom),ii=1,3)
      close(31)
      do katom=1,InVar%natom
        write(31,'(a,x,3(f10.6,x))') 'I',xred_ideal (:,katom)
        write(31,'(a,x,3(f10.6,x))') 'C',xred_center(:,katom)
      end do
      stop
    end if  
  end do

! WARNING: VERY IMPORTANT: The positions are displayed/sorted (and used in the following) 
! according to ideal positions xred_ideal. 
! --> In reduced coordinates
  do iatom=1,InVar%natom
    write(31,'(a,x,3(f10.6,x))') 'Ired',xred_ideal (:,iatom)
    write(31,'(a,x,3(f10.6,x))') 'Cred',xred_center(:,FromIdeal2Average(iatom))
  end do  
! --> In cartesian coordinates  
  do iatom=1,InVar%natom
    tmp(:)=zero
    call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,xred_ideal (:,iatom),1,0.d0,tmp(:),1)
    write(31,'(a,x,3(f10.6,x))') 'Icart',tmp(:)
    tmp(:)=zero
    call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,xred_center(:,FromIdeal2Average(iatom)),1,0.d0,tmp(:),1)
    write(31,'(a,x,3(f10.6,x))') 'Ccart',tmp(:)
  end do  
  close(31)

! Average distances between atoms --> distance_average
  allocate(distance_average(InVar%natom,InVar%natom,4))      ; distance_average(:,:,:)=0.d0
  do katom=1,InVar%natom
    do latom=1,InVar%natom
      tmp(:)=xred_center(:,FromIdeal2Average(latom))-xred_center(:,FromIdeal2Average(katom))
      call inbox(tmp,1,1d-3)
      call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,tmp(:),1,0.d0,distance_average(katom,latom,2:4),1)
      do ii=1,3
!       Remove the rounding errors before writing (for non regression testing purposes)
        if (abs(distance_average(katom,latom,ii+1)).lt.tol8) distance_average(katom,latom,ii+1)=zero
        distance_average(katom,latom,1)=distance_average(katom,latom,1)+(distance_average(katom,latom,ii+1))**2
      end do
      distance_average(katom,latom,1)=distance_average(katom,latom,1)**0.5
    end do  
  end do  
  deallocate(xred_center)

!====================================================================================
!====================== END OF REDUCED COORDINATES ==================================
!====================================================================================
! a/ Get cartesian coordinates from reduced ones
! b/ Compute ucart and fcart tabs
! c/ The atoms are sorted according the IDEAL arrangement
!    The correspondance function is contained in: FromIdeal2Average
!    WARNING : Consequently the arrangement of the xcart* tabs is not modified. 
  allocate(xcart        (3,InVar%natom,InVar%nstep))       ; xcart(:,:,:)=0.d0
  allocate(xcart_ideal  (3,InVar%natom))                   ; xcart_ideal(:,:)=0.d0
  allocate(xcart_average(3,InVar%natom))                   ; xcart_average(:,:)=0.d0
  allocate(ucart_tmp    (3,InVar%natom,InVar%nstep_remain)); ucart_tmp(:,:,:)=0.d0
  do iatom=1,InVar%natom
    call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,xred_ideal  (:,iatom),1,0.d0,xcart_ideal  (:,iatom),1)
    call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,xred_average(:,iatom),1,0.d0,xcart_average(:,iatom),1)
    do iatcell=1,InVar%natom_unitcell
      call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,Rlatt_red(:,iatcell,iatom),1,0.d0,Rlatt_cart(:,iatcell,iatom),1)
    end do  
  end do
  do istep=1,InVar%nstep_remain
    do iatom=1,InVar%natom
      call DGEMV('T',3,3,1.d0,Lattice%rprimd_MD(:,:),3,InVar%xred(:,FromIdeal2Average(iatom),InVar%nstep1+istep-1),1,0.d0,xcart(:,FromIdeal2Average(iatom),InVar%nstep1+istep-1),1)
      if (InVar%Use_ideal_positions.eq.0) then
        ucart_tmp(:,iatom,istep)=xcart(:,FromIdeal2Average(iatom),InVar%nstep1+istep-1)-xcart_average(:,FromIdeal2Average(iatom))
      else
        ucart_tmp(:,iatom,istep)=xcart(:,FromIdeal2Average(iatom),InVar%nstep1+istep-1)-xcart_ideal  (:,iatom)
      end if  
    end do
  end do
  deallocate(xred_average)

! Rearrangement of the fcart tabs in column --> fcartij
  allocate(fcart_tmp(3,InVar%natom,InVar%nstep)); fcart_tmp(:,:,:)=0.d0
  do istep=1,InVar%nstep
    do iatom=1,InVar%natom
      fcart_tmp(:,iatom,istep)=InVar%fcart(:,FromIdeal2Average(iatom),istep)
    end do  
  end do  
  jstep=0
  do istep=1,InVar%nstep_remain
    jstep=jstep+1
    do jatom=1,InVar%natom
      do ii=1,3 
        fcartij(ii+3*(jatom-1)+3*InVar%natom*(jstep-1))=fcart_tmp(ii,jatom,InVar%nstep1+istep-1)
        ucart(ii,jatom,jstep)=ucart_tmp(ii,jatom,istep)
      enddo  
    enddo
  enddo  
  deallocate(ucart_tmp)
  deallocate(xcart_average)
  deallocate(fcart_tmp)

! Define Rlatt4dos, fulfilling the definition of mkphdos (ABINIT routine)
  do ii=1,3
    rprimd_MD_tmp(ii,:)=Lattice%rprimd_MD(ii,:)/Lattice%acell_unitcell(ii)
  end do  
  do iatom=1,InVar%natom
    do iatcell=1,InVar%natom_unitcell
      call DGEMV('T',3,3,1.d0,rprimd_MD_tmp(:,:),3,Rlatt_red(:,iatcell,iatom),1,0.d0,Rlatt4dos(:,iatcell,iatom),1)
    end do  
  end do
  deallocate(Rlatt_red)
  deallocate(InVar%fcart)

!==========================================================================================
!======== 3/ Find the symetry operation between the reference and image bonds =============
!==========================================================================================
  call SearchMatR_1at(InVar,Lattice,Sym,xred_ideal)

  allocate(Sym%matR(InVar%natom,InVar%natom)) ; Sym%matR(:,:)=zero
  ishell=0
  do iatcell=1,InVar%natom_unitcell
    do jatom=1,InVar%natom
!     Interactions are only computed until Rcut<acell/2 in order to have complete shell of neighbours.
!     Otherwise the symetries are broken.
      if ((bond_ref(iatcell,jatom,1).ne.0).or.(distance(iatcell,jatom,1).gt.(InVar%Rcut*0.99))) cycle
      ishell=ishell+1
      do katom=1,InVar%natom
        do latom=1,InVar%natom
          if ((bond_ref(katom,latom,1).eq.0).and.(bond_ref(katom,latom,2).eq.0).and.&
&            (abs(distance(iatcell,jatom,1)-distance(katom,latom,1)).lt.1.d-3)) then
            
            call SearchMatR_2at(InVar,iatcell,jatom,katom,latom,Sym,xred_ideal)
            if (Sym%matR(katom,latom).ne.0) then
                bond_ref(katom,latom,1)=iatcell
                bond_ref(katom,latom,2)=jatom
                bond_ref(katom,latom,3)=ishell
            end if  
          end if
        end do !latom
      end do !katom
    end do !jatom 
  end do !iatcell 
  nshell=ishell
  deallocate(xred_ideal)

!==========================================================================================
!======== 4/ Write output quantities needed to visualize the neighbouring distances =======
!==========================================================================================

  if (InVar%d1NN.gt.0) then
    call write_d1NN(distance_average,FromIdeal2Average,InVar,bond_ref,Sym,xcart,xcart_ideal)
  end if
  deallocate(InVar%xred)
  deallocate(distance_average)
  deallocate(FromIdeal2Average)
  deallocate(xcart)
  deallocate(xcart_ideal)


 end subroutine MatchIdeal2Average

!====================================================================================================
 subroutine model(DeltaFree_AH2,Forces_MD,InVar,Phij_NN,ucart,U0)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'model'
!End of the abilint section

  implicit none 

  type(Input_Variables_type),intent(in) :: InVar
  double precision, intent(in)  :: Phij_NN(3*InVar%natom,         3*InVar%natom)
  double precision, intent(in)  :: ucart(3,InVar%natom,InVar%nstep_remain)
  double precision, intent(in)  :: Forces_MD(3*InVar%natom*InVar%nstep_remain)
  double precision, intent(out)  :: U0,DeltaFree_AH2
  
  integer :: ii,jj,istep,iatom,jatom,islice,nslice,istepmin,istepmax
  double precision :: force_i,MinDeltaForces,DeltaFree_AH,tmp0,tmp1,tmp2,tmp3
  double precision, allocatable :: U_MD(:),U_TDEP(:),PhijUiUj(:),Forces_TDEP(:)

  write(InVar%stdout,*)' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '######################### Energies, errors,...  #############################'
  write(InVar%stdout,*) '#############################################################################'
  allocate(U_MD(InVar%nstep_remain)); U_MD(:)=0.d0
  do istep=1,InVar%nstep_remain
    U_MD(istep)=InVar%etot(InVar%nstep1+istep-1)
  end do  
  
! Compute Forces of the model TDEP
  allocate(PhijUiUj(InVar%nstep_remain)); PhijUiUj(:)=0.d0 
  allocate(Forces_TDEP(3*InVar%natom*InVar%nstep_remain)); Forces_TDEP(:)=0.d0 
  do istep=1,InVar%nstep_remain
    do iatom=1,InVar%natom
      do ii=1,3
        force_i=0.d0
        do jatom=1,InVar%natom
          do jj=1,3
            force_i=force_i-Phij_NN((iatom-1)*3+ii,3*(jatom-1)+jj)*ucart(jj,jatom,istep)
            PhijUiUj(istep)=PhijUiUj(istep)+Phij_NN((iatom-1)*3+ii,3*(jatom-1)+jj)*ucart(ii,iatom,istep)*ucart(jj,jatom,istep)/2.d0
          end do  
        end do
        Forces_TDEP(ii+3*(iatom-1)+3*InVar%natom*(istep-1))=force_i
      end do
    end do
  end do

! Compute U0, U_TDEP, DeltaFree_AH and write them in the data.out file
  allocate(U_TDEP(InVar%nstep_remain)); U_TDEP(:)=0.d0
  write(InVar%stdout,'(a)') ' Average quantities highlighting the convergence (in a.u./atom) :'
  write(InVar%stdout,'(a)') ' Istep                 U_0           DeltaFree_AH          DeltaFree_AH2        (F_MD-F_TDEP)**2'
  tmp0=zero
  tmp3=zero
  nslice=10
  do islice=1,nslice
    istepmin=int(InVar%nstep_remain*(islice-1)/nslice+1)
    istepmax=int(InVar%nstep_remain*islice/nslice)
!   Compute U0 = < U_MD - \Sum_ij \Phi_ij u_i u_j >   
    do istep=istepmin,istepmax
      tmp0=tmp0+(U_MD(istep)-PhijUiUj(istep))
    end do
!   Compute U_TDEP = U0 + PhijUiUj  
    U_TDEP(:)=tmp0/dfloat(istepmax)+PhijUiUj(:)
!   Compute DeltaFree_AH = < U_MD - U_TDEP >
    tmp1=zero
    do istep=1,istepmax
      tmp1 =tmp1 +(U_MD(istep)-U_TDEP(istep))
    end do
!   Compute DeltaFree_AH2 = - <( (U_MD-U_TDEP) - <U_MD-U_TDEP> )**2> / (2.k_B.T)
    tmp2=zero
    do istep=1,istepmax
      tmp2=tmp2-( (U_MD(istep)-U_TDEP(istep)) - tmp1/dfloat(istepmax))**2
    end do
!   Compute eucledian distance for forces    
    do istep=istepmin,istepmax
      do iatom=1,InVar%natom
        do ii=1,3
          tmp3=tmp3+(Forces_MD(ii+3*(iatom-1)+3*InVar%natom*(istep-1))&
&                 -Forces_TDEP(ii+3*(iatom-1)+3*InVar%natom*(istep-1)))**2
        end do
      end do  
    end do
    U0            =tmp0/dfloat(istepmax*InVar%natom)
    DeltaFree_AH  =tmp1/dfloat(istepmax*InVar%natom)
    DeltaFree_AH2 =tmp2/dfloat(istepmax*InVar%natom)/(2.d0*kb_HaK*InVar%temperature)
    MinDeltaForces=tmp3/dfloat(istepmax*InVar%natom*3)
    write(InVar%stdout,'(i5,4(10x,e12.5))') istepmax,U0,DeltaFree_AH,DeltaFree_AH2,MinDeltaForces
  end do  
  deallocate(PhijUiUj)

! Write (U_TDEP vs U_MD) and (Forces_TDEP vs Forces_MD) 
  write(InVar%stdout,'(a)') ' '
  write(InVar%stdout,'(a)') ' See the etotMDvsTDEP.dat file and (if DEBUG) the fcartMDvsTDEP.dat file'
  open(unit=32,file='etotMDvsTDEP.dat')
  open(unit=33,file='fcartMDvsTDEP.dat')
  do istep=1,InVar%nstep_remain
    write(32,'(i6,x,2(e17.10,x))') istep,U_MD(istep),U_TDEP(istep)
    if (InVar%debug) then
      do iatom=1,InVar%natom
        do ii=1,3
          write(33,'(2(e17.10,x))') Forces_MD  (ii+3*(iatom-1)+3*InVar%natom*(istep-1)),&
&                                   Forces_TDEP(ii+3*(iatom-1)+3*InVar%natom*(istep-1))
        end do
      end do  
    end if  
  end do  
  close(32)
  close(33)
  deallocate(U_MD,U_TDEP)
  deallocate(Forces_TDEP)

 end subroutine model
!====================================================================================================

end module m_utils
