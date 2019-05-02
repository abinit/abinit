
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_sym

 use defs_basis
 use m_abicore
 use m_errors

 use m_symtk,            only : mati3inv, symatm
 use m_symfind,          only : symfind, symanal, symlatt
 use m_tdep_latt,        only : Lattice_Variables_type, tdep_make_inbox
 use m_tdep_readwrite,   only : Input_Variables_type

 implicit none

 type,public :: Symetries_Variables_type

   integer :: msym
   integer :: nsym
   integer :: nptsym
   integer :: spgroup
   integer, allocatable :: ptsymrel(:,:,:)
   integer, allocatable :: symrec(:,:,:)
   integer, allocatable :: symafm(:)
   integer, allocatable :: indsym(:,:,:)
   double precision, allocatable :: S_ref(:,:,:,:)
   double precision, allocatable :: S_inv(:,:,:,:)
   double precision, allocatable :: tnons(:,:)
   double precision, allocatable :: xred_zero(:,:)

 end type Symetries_Variables_type

 public :: tdep_make_sym
 public :: tdep_SearchS_1at
 public :: tdep_SearchS_2at
 public :: tdep_SearchS_3at
 public :: tdep_calc_indsym2

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_make_sym(Invar,Lattice,Sym)

  implicit none

  integer :: isym,ii,jj
  double precision, allocatable :: tmp1(:,:)
  type(Symetries_Variables_type),intent(out) :: Sym
  type(Input_Variables_type),intent(inout) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice

! Compute all the symetries coming from the bravais lattice
! The routine used is symlatt (from Abinit code)
  Sym%msym=1000 !msym needs to be very large due to non-primitive cell calculations
  ABI_MALLOC(Sym%ptsymrel,(3,3,Sym%msym)) ; Sym%ptsymrel(:,:,:)=0
  call symlatt(InVar%bravais,Sym%msym,Sym%nptsym,Sym%ptsymrel,Lattice%rprimdt,tol8)
  write(Invar%stdout,'(a,1x,11(i4,1x))')' bravais=',InVar%bravais(:)
  Sym%nsym=Sym%nptsym

! Transform S_ref in cartesian coordinates
  ABI_MALLOC(Sym%S_ref,(3,3,Sym%nsym,2)) ; Sym%S_ref(:,:,:,1)=real(Sym%ptsymrel(:,:,1:Sym%nsym))
  ABI_MALLOC(Sym%S_inv,(3,3,Sym%nsym,2)) ; Sym%S_inv(:,:,:,1)=zero
  ABI_MALLOC(tmp1,(3,3)); tmp1(:,:)=0.d0
  open(unit=75,file=trim(InVar%output_prefix)//'sym.dat')
  do isym=1,Sym%nsym
    write(75,*) ' '
    write(75,*) 'For isym=',isym
    write(75,*) 'In reduced coordinates:'
    write(75,'(3(f5.2,1x))') Sym%S_ref(1,1,isym,1),Sym%S_ref(1,2,isym,1),Sym%S_ref(1,3,isym,1)
    write(75,'(3(f5.2,1x))') Sym%S_ref(2,1,isym,1),Sym%S_ref(2,2,isym,1),Sym%S_ref(2,3,isym,1)
    write(75,'(3(f5.2,1x))') Sym%S_ref(3,1,isym,1),Sym%S_ref(3,2,isym,1),Sym%S_ref(3,3,isym,1)
    call DGEMM('N','N',3,3,3,1.d0,Lattice%rprimdt,3,Sym%S_ref(:,:,isym,1),3,0.d0,tmp1,3)
    call DGEMM('N','N',3,3,3,1.d0,tmp1,3,Lattice%rprimdtm1,3,0.d0,Sym%S_ref(:,:,isym,1),3)
    write(75,*) 'In cartesian coordinates:'
    write(75,'(3(f5.2,1x))') Sym%S_ref(1,1,isym,1),Sym%S_ref(1,2,isym,1),Sym%S_ref(1,3,isym,1)
    write(75,'(3(f5.2,1x))') Sym%S_ref(2,1,isym,1),Sym%S_ref(2,2,isym,1),Sym%S_ref(2,3,isym,1)
    write(75,'(3(f5.2,1x))') Sym%S_ref(3,1,isym,1),Sym%S_ref(3,2,isym,1),Sym%S_ref(3,3,isym,1)

    do ii=1,3
      do jj=1,3
        if (abs(Sym%S_ref(ii,jj,isym,1)).lt.tol8) Sym%S_ref(ii,jj,isym,1)=zero
      end do
    end do

!   Inversion of S_ref (which is equivalent to the transposition)
    do ii=1,3
      do jj=1,3
        Sym%S_inv(ii,jj,isym,1)=Sym%S_ref(jj,ii,isym,1)
      end do
    end do

  end do
  write(InVar%stdout,'(a)') ' See the sym.dat file'

! We verify that S.S^T = Id
  tmp1(:,:)=zero
  do isym=1,Sym%nsym
    call DGEMM('T','N',3,3,3,1.d0,Sym%S_ref(:,:,isym,1),3,Sym%S_ref(:,:,isym,1),3,0.d0,tmp1,3)
    do ii=1,3
      do jj=1,3
        if ((ii/=jj.and.abs(tmp1(ii,jj)).gt.tol8).or.(ii==jj.and.abs(tmp1(ii,jj)-1.d0).gt.tol8)) then
          write(InVar%stdout,'(a)') ' STOP : the matrix is not orthogonal : '
          write(InVar%stdout,'(a,1x,i3,1x,i3,1x,f15.10)') ' ii, jj, sym(ii,jj)=',ii,jj,tmp1(ii,jj)
          MSG_ERROR('The matrix is not orthogonal')
        end if
      end do
    end do
  end do

  ABI_FREE(tmp1)

 end subroutine tdep_make_sym

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine tdep_SearchS_1at(InVar,Lattice,Sym,xred_ideal)

  implicit none

  type(Input_Variables_type), intent(in) :: InVar
  type(Lattice_Variables_type), intent(in) :: Lattice
  type(Symetries_Variables_type), intent(inout) :: Sym
  double precision, intent(in) :: xred_ideal(3,InVar%natom)

  integer :: ptgroupma,isym,jatom,ii,jj,mu
  integer :: nsym,iatom_unitcell,jatom_unitcell,iatom
  integer :: berryopt,jellslab,noncoll,use_inversion,chkprim,counter,nzchempot
  integer :: bravais(11)
  integer :: vecti(3),vectj(3),vectsym(4,Sym%nptsym)
  integer, allocatable :: symrel(:,:,:),symrel_tmp(:,:,:)
  integer, allocatable :: indsym2(:,:,:,:)
  double precision :: genafm(3),efield(3)
  double precision :: tmpi(3,InVar%natom),tmpj(3,InVar%natom),temp3(3,1),tmp_store(3,InVar%natom_unitcell)
  double precision, allocatable :: spinat(:,:),xred_temp(:,:),tnons_tmp(:,:)

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '###################### Find the symetry operations ##########################'
  write(InVar%stdout,*) '#################### (connecting the atoms together) ########################'
  write(InVar%stdout,*) '#############################################################################'

!  write(InVar%stdout,'(a)') 'Begining SearchS'
! TODO: Initialize the local variables at 0

! Initialize the all the symetries using the first atom at (0.0;0.0;0.0)
  ABI_MALLOC(    xred_temp,(3,InVar%natom_unitcell))        ;     xred_temp(:,:)=0.d0
  ABI_MALLOC(Sym%xred_zero,(3,InVar%natom_unitcell))        ; Sym%xred_zero(:,:)=0.d0
  if (InVar%natom_unitcell.gt.1) then
    do iatom_unitcell=2,InVar%natom_unitcell
      temp3(:,1)=xred_ideal(:,iatom_unitcell)-xred_ideal(:,1)
      call DGEMV('T',3,3,1.d0,InVar%multiplicity(:,:),3,temp3(:,1),1,0.d0,xred_temp(:,iatom_unitcell),1)
      Sym%xred_zero(:,iatom_unitcell)=xred_temp(:,iatom_unitcell)
    end do
  end if

!  write(InVar%stdout,'(a)') 'Before tnons'

! Calcul des tnons
  berryopt=0
  efield(:)=0.d0
  jellslab=0
  noncoll=0
  nzchempot=0
  ABI_MALLOC(spinat,(3,InVar%natom_unitcell)); spinat(:,:)=0.d0
  use_inversion=1
  ABI_MALLOC(symrel    ,(3,3,Sym%msym)) ; symrel    (:,:,:)=0
  ABI_MALLOC(symrel_tmp,(3,3,Sym%msym)) ; symrel_tmp(:,:,:)=0
  ABI_MALLOC(Sym%tnons   ,(3,Sym%msym)) ; Sym%tnons   (:,:)=0.d0
  ABI_MALLOC(Sym%symafm    ,(Sym%msym)) ; Sym%symafm(:)    =1
  ABI_MALLOC(tnons_tmp   ,(3,Sym%msym)) ; tnons_tmp   (:,:)=0.d0
  call symfind(berryopt,efield,Lattice%gprimd,jellslab,Sym%msym,InVar%natom_unitcell,noncoll,Sym%nptsym,nsym,&
&      nzchempot,0,Sym%ptsymrel,spinat,Sym%symafm,symrel_tmp,tnons_tmp,tol8,InVar%typat_unitcell,use_inversion,xred_temp)
  if (Sym%nptsym.eq.nsym) then
    symrel(:,:,:) =symrel_tmp(:,:,:)
    Sym%tnons(:,:)=tnons_tmp   (:,:)
  else
    write(InVar%stdout,'(a,i4,a,i4)') 'WARNING: nsym=',nsym,' is not equal to nptsym=',Sym%nptsym
!   Warning: EXPERIMENTAL
!   This case is only useful for testing purpose
!   when the unitcell is not primitive (chkprim=0)
    if (nsym.eq.Sym%nptsym*InVar%natom_unitcell) then
      write(InVar%stdout,'(a)') 'The unitcell is probably non-primitive, with probably nsym=nptsym*natom_unitcell'
      do isym=1,Sym%nptsym
        symrel(:,:,isym)=symrel_tmp(:,:,(isym-1)*InVar%natom_unitcell+1)
        Sym%tnons(:,isym)=tnons_tmp(:,(isym-1)*InVar%natom_unitcell+1)
      end do
    else
      MSG_ERROR('There is no relationship between nsym and nptsym')
    end if
  end if
  ABI_FREE(symrel_tmp)
  ABI_FREE(tnons_tmp)
  do isym=1,Sym%nptsym
    counter=0
    do ii=1,3
      do jj=1,3
        counter=counter+abs(Sym%ptsymrel(ii,jj,isym)-symrel(ii,jj,isym))
      end do
    end do
    if (counter.ne.0) then
      exit
    end if
  end do

  if (counter.ne.0) then
    do isym=1,nsym
      write(InVar%stdout,*) ' For sym=',isym
      write(InVar%stdout,'(a)') ' Symrel='
      write(InVar%stdout,*) symrel(1,:,isym)
      write(InVar%stdout,*) symrel(2,:,isym)
      write(InVar%stdout,*) symrel(3,:,isym)
      write(InVar%stdout,'(a)') ' Ptsymrel='
      write(InVar%stdout,*) Sym%ptsymrel(1,:,isym)
      write(InVar%stdout,*) Sym%ptsymrel(2,:,isym)
      write(InVar%stdout,*) Sym%ptsymrel(3,:,isym)
    end do
    MSG_ERROR('The symrel and ptsymrel tabs do not correspond')
  end if

!Calcul des symrec
  ABI_MALLOC(Sym%symrec,(3,3,Sym%nptsym)); Sym%symrec(:,:,:)=0
  do isym=1,Sym%nptsym
    call mati3inv(symrel(:,:,isym),Sym%symrec(:,:,isym))
  end do
! Calcul du spgroup
!FB  chkprim=1
  chkprim=0
  genafm(:)=0.d0
  ptgroupma=0
  Sym%spgroup=0
  call symanal(bravais,chkprim,genafm,Sym%msym,Sym%nptsym,ptgroupma,Lattice%rprimdt,Sym%spgroup,Sym%symafm,symrel,Sym%tnons,tol8)
!FB  call symanal(bravais,chkprim,genafm,Sym%msym,Sym%nptsym,ptgroupma,Lattice%rprimd,Sym%spgroup,Sym%symafm,symrel,Sym%tnons,tol8)
  ABI_FREE(symrel)

! Compute the indsym fundamental quantity
! === Obtain a list of rotated atoms ===
! $ R^{-1} (xred(:,iat)-\tau) = xred(:,iat_sym) + R_0 $
! * indsym(4,  isym,iat) gives iat_sym in the original unit cell.
! * indsym(1:3,isym,iat) gives the lattice vector $R_0$.
  write(InVar%stdout,'(a)') ' Search the matrix transformation going from (k) to (i)...'
  ABI_MALLOC(Sym%indsym,(4,Sym%nptsym,InVar%natom)); Sym%indsym(:,:,:)=zero
  call symatm(Sym%indsym(:,:,1:InVar%natom_unitcell),InVar%natom_unitcell,Sym%nptsym,&
&   Sym%symrec,Sym%tnons,tol8,InVar%typat_unitcell,xred_temp)

! Store the positions of the atoms in the motif
  do iatom=1,InVar%natom_unitcell
    call DGEMV('T',3,3,1.d0,InVar%multiplicity(:,:),3,xred_ideal(:,iatom),1,0.d0,tmp_store(:,iatom),1)
  end do
! Write the Indsym of the atoms included in the (reference) unitcell (i.e.: the motif)
  if (InVar%debug) then
    open(unit=40,file=trim(InVar%output_prefix)//'Indsym-unitcell.dat')
    do iatom=1,InVar%natom_unitcell
      write(40,*) '=========================================='
      write(40,'(a,i4,a,3(f10.5,1x))') 'For iatom=',iatom,' with xred (supercell)=',xred_ideal(:,iatom)
      do isym=1,Sym%nptsym
        write(40,'(a,i2,a,i4,a,3(i4,1x),a,i2,a,3(f10.5,1x))') '  indsym(isym=',isym,',',iatom,')=',&
&         Sym%indsym(1:3,isym,iatom),'|iat=',Sym%indsym(4,isym,iatom),'| with tnons=',Sym%tnons(:,isym)
      end do
    end do
    close(40)
  end if

! Search the matrix transformation going from (k,l) to (i,j)
  write(InVar%stdout,'(a)') ' Search the matrix transformation going from (k,l) to (i,j)...'
  ABI_MALLOC(indsym2,(8,Sym%nptsym,InVar%natom,InVar%natom)); indsym2(:,:,:,:)=0
  tmpi(:,:)=0.d0
  tmpj(:,:)=0.d0
  if (InVar%debug) open(unit=40,file=trim(InVar%output_prefix)//'Indsym-supercell.dat')
  do iatom=1,InVar%natom
    if (InVar%debug) write(40,*) '=========================================='
    if (InVar%debug) write(40,'(a,i4,a,3(f10.5,1x))') 'For iatom=',iatom,' with xred (supercell)=',xred_ideal(:,iatom)
!   For a single iatom 
    call DGEMV('T',3,3,1.d0,InVar%multiplicity(:,:),3,xred_ideal(:,iatom),1,0.d0,tmpi(:,iatom),1)
    iatom_unitcell=mod(iatom-1,InVar%natom_unitcell)+1
    vecti(:)=nint(tmpi(:,iatom)-tmp_store(:,iatom_unitcell))
    do isym=1,Sym%nptsym
      vectsym(:,:)=0
      do mu=1,3 ! Apply inverse transformation to original coordinates. Note transpose of symrec.
        vectsym(mu,isym) = Sym%symrec(1,mu,isym)*vecti(1)+Sym%symrec(2,mu,isym)*vecti(2)+Sym%symrec(3,mu,isym)*vecti(3)
      end do
      Sym%indsym(1:4,isym,iatom)=Sym%indsym(1:4,isym,iatom_unitcell)+vectsym(1:4,isym)
      if (InVar%debug) then
        write(40,'(a,i2,a,i4,a,3(i4,1x),a,i2,a,3(f10.5,1x))') '  indsym(isym=',isym,',',iatom,')=',&
&         Sym%indsym(1:3,isym,iatom),'|iat=',Sym%indsym(4,isym,iatom),'| with tnons=',Sym%tnons(:,isym)
      end if
      do jatom=1,InVar%natom
        indsym2(1:4,isym,iatom,jatom)=Sym%indsym(1:4,isym,iatom_unitcell)+vectsym(1:4,isym)
      end do
    end do
    if (InVar%debug) write(40,'(a,i4)') ' '
  end do
  if (InVar%debug) close(40)

! For a couple of (iatom,jatom). The (iatom,jatom) vector depends on the position of iatom (due to PBC)
  if (InVar%debug) open(unit=40,file=trim(InVar%output_prefix)//'Indsym-2atoms.dat')
  do iatom=1,InVar%natom
    if (InVar%debug) write(40,*) '=========================================='
    if (InVar%debug) write(40,'(a,i4,a,3(f10.5,1x))') 'For iatom=',iatom,' with xred (supercell)=',xred_ideal(:,iatom)
    do jatom=1,InVar%natom
      if (InVar%debug) write(40,*) '  =========================================='
      if (InVar%debug) write(40,'(a,i4,a,3(f10.5,1x))') '  For jatom=',jatom,' with xred (supercell)=',xred_ideal(:,jatom)
      temp3(:,1)=xred_ideal(:,jatom)-xred_ideal(:,iatom)
      call tdep_make_inbox(temp3(:,1),1,1d-3,temp3(:,1))
      temp3(:,1)=xred_ideal(:,iatom)+temp3(:,1)
      call DGEMV('T',3,3,1.d0,InVar%multiplicity(:,:),3,temp3(:,1),1,0.d0,tmpj(:,jatom),1)
      jatom_unitcell=mod(jatom-1,InVar%natom_unitcell)+1
      vectj(:)=nint(tmpj(:,jatom)-tmp_store(:,jatom_unitcell))
      do isym=1,Sym%nptsym
        vectsym(:,:)=0
        do mu=1,3 ! Apply inverse transformation to original coordinates. Note transpose of symrec.
          vectsym(mu,isym) = Sym%symrec(1,mu,isym)*vectj(1)+Sym%symrec(2,mu,isym)*vectj(2)+Sym%symrec(3,mu,isym)*vectj(3)
        end do
        indsym2(5:8,isym,iatom,jatom)=Sym%indsym(1:4,isym,jatom_unitcell)+vectsym(1:4,isym)
        if (InVar%debug) then
          write(40,'(a,i2,a,i4,a,i4,a,3(i4,1x),a,i2,a,3(i4,1x),a,i2,a)') '  indsym2(isym=',isym,',',iatom,',',jatom,')=',&
&           indsym2(1:3,isym,iatom,jatom),&
&           '|iat=',indsym2(4,isym,iatom,jatom),'|',indsym2(5:7,isym,iatom,jatom),'|iat=',indsym2(8,isym,iatom,jatom),'|'
        end if
      end do
    end do
    if (InVar%debug) write(40,'(a,i4)') ' '
  end do
  if (InVar%debug) close(40)
  ABI_FREE(indsym2)
  write(InVar%stdout,'(a)') ' See the Indsym*.dat files (if debug)'

 end subroutine tdep_SearchS_1at

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_SearchS_2at(InVar,iatom,jatom,eatom,fatom,Isym2at,Sym,xred_ideal)

  implicit none

  integer, intent(in) :: iatom,jatom,eatom,fatom
  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  integer, intent(inout) :: Isym2at(InVar%natom,InVar%natom,2)
  double precision, intent(in) :: xred_ideal(3,InVar%natom)

  integer :: isym,ee,ff,ii
  integer :: iatom_unitcell,jatom_unitcell
  integer :: vecti(3),vectj(3),latte(3),lattf(3),indsym2(8)
  double precision :: temp(3),tmp_store(3,InVar%natom_unitcell)
  double precision :: tmp(3,InVar%natom)
  logical :: ok

! Search the couple of atoms (e,f) obtained through the transformation (R+t)
! starting from (i,j)
! In that case, note that two cases are possible:
! - The bond is just transformed, so indsym(4,e)=i and indsym(4,f)=j
! - The bond is transformed and reversed, so indsym(4,e)=j and indsym(4,f)=i

! Store the positions of the atoms in the motif
  do ii=1,InVar%natom_unitcell
    call DGEMV('T',3,3,1.d0,InVar%multiplicity(:,:),3,xred_ideal(:,ii),1,0.d0,tmp_store(:,ii),1)
  end do

! Search the atom equivalent to iatom in the (reference) unitcell
  do ii=1,InVar%natom
    call DGEMV('T',3,3,1.d0,InVar%multiplicity(:,:),3,xred_ideal(:,ii),1,0.d0,tmp(:,ii),1)
  end do
! Note that in the (present) particular case : iatom_unitcell=iatom and vecti(:)=zero
  iatom_unitcell=mod(iatom-1,InVar%natom_unitcell)+1
  vecti(:)=nint(tmp(:,iatom)-tmp(:,iatom_unitcell))

! Search the atom equivalent to jatom in the (reference) unitcell.
! jatom can be outside the box (according to the value of iatom)
! so we use the inbox procedure to put the distance within [-0.5,0.5[
  temp(:)=xred_ideal(:,jatom)-xred_ideal(:,iatom)
  call tdep_make_inbox(temp,1,1d-3,temp)
  temp(:)=xred_ideal(:,iatom)+temp(:)
  call DGEMV('T',3,3,1.d0,InVar%multiplicity(:,:),3,temp(:),1,0.d0,tmp(:,jatom),1)
  jatom_unitcell=mod(jatom-1,InVar%natom_unitcell)+1
  vectj(:)=nint(tmp(:,jatom)-tmp_store(:,jatom_unitcell))

! To understand the meaning of "latt", see SearchS_1at
  ok=.false.
  do isym=1,Sym%nsym
    indsym2(:)=0
    call tdep_calc_indsym2(InVar,eatom,fatom,indsym2,isym,Sym,xred_ideal)
    ee=indsym2(4)
    ff=indsym2(8)
    latte(:)=indsym2(1:3)
    lattf(:)=indsym2(5:7)
    if (ee==iatom_unitcell.and.ff==jatom_unitcell) then
!FB      if (InVar%debug) then
!FB        write(InVar%stdout,*) 'For ee=iatom and ff=jatom, isym=',isym
!FB        write(InVar%stdout,'(4(a,i4))') ' + ',vecti(1),' - ',vectj(1),' - ',latte(1),' + ',lattf(1)
!FB        write(InVar%stdout,'(4(a,i4))') ' + ',vecti(2),' - ',vectj(2),' - ',latte(2),' + ',lattf(2)
!FB        write(InVar%stdout,'(4(a,i4))') ' + ',vecti(3),' - ',vectj(3),' - ',latte(3),' + ',lattf(3)
!FB      end if
      if (sum(abs((vecti(:)-latte(:))-(vectj(:)-lattf(:)))).lt.tol8) then
        Isym2at(eatom,fatom,1)=isym
        Isym2at(eatom,fatom,2)=1
        ok=.true.
      end if
    end if
    if (ok) exit
  end do

 end subroutine tdep_SearchS_2at

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_SearchS_3at(InVar,iatom,jatom,katom,eatom,fatom,gatom,Isym3at,Sym,xred_ideal)

  implicit none

  integer, intent(in) :: iatom,jatom,katom,eatom,fatom,gatom
  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  integer, intent(inout) :: Isym3at(2)
  double precision, intent(in) :: xred_ideal(3,InVar%natom)

  integer :: isym,ee,ff,gg,ii
  integer :: iatom_unitcell,jatom_unitcell,katom_unitcell
  integer :: vecti(3),vectj(3),vectk(3),indsym2(8)
  integer :: lattef(3),lattfg(3),lattge(3),lattfe(3),lattgf(3),latteg(3)
  double precision :: temp(3),tmp_store(3,InVar%natom_unitcell)
  double precision :: tmp(3,InVar%natom)
  logical :: ok

! Search the couple of atoms (e,f) obtained through the transformation (R+t)
! starting from (i,j)
! In that case, note that two cases are possible:
! - The bond is just transformed, so indsym(4,e)=i and indsym(4,f)=j
! - The bond is transformed and reversed, so indsym(4,e)=j and indsym(4,f)=i

! Store the positions of the atoms in the motif
  do ii=1,InVar%natom_unitcell
    call DGEMV('T',3,3,1.d0,InVar%multiplicity(:,:),3,xred_ideal(:,ii),1,0.d0,tmp_store(:,ii),1)
  end do

! Search the atom equivalent to iatom in the (reference) unitcell
  do ii=1,InVar%natom
    call DGEMV('T',3,3,1.d0,InVar%multiplicity(:,:),3,xred_ideal(:,ii),1,0.d0,tmp(:,ii),1)
  end do
! Note that in the (present) particular case : iatom_unitcell=iatom and vecti(:)=zero
  iatom_unitcell=mod(iatom-1,InVar%natom_unitcell)+1
  vecti(:)=nint(tmp(:,iatom)-tmp(:,iatom_unitcell))

! Search the atom equivalent to jatom in the (reference) unitcell.
! jatom can be outside the box (according to the value of iatom)
! so we use the inbox procedure to put the distance within [-0.5,0.5[
  temp(:)=xred_ideal(:,jatom)-xred_ideal(:,iatom)
  call tdep_make_inbox(temp,1,1d-3,temp)
  temp(:)=xred_ideal(:,iatom)+temp(:)
  call DGEMV('T',3,3,1.d0,InVar%multiplicity(:,:),3,temp(:),1,0.d0,tmp(:,jatom),1)
  jatom_unitcell=mod(jatom-1,InVar%natom_unitcell)+1
  vectj(:)=nint(tmp(:,jatom)-tmp_store(:,jatom_unitcell))

! Search the atom equivalent to katom in the (reference) unitcell.
! katom can be outside the box (according to the value of iatom)
! so we use the inbox procedure to put the distance within [-0.5,0.5[
  temp(:)=xred_ideal(:,katom)-xred_ideal(:,iatom)
  call tdep_make_inbox(temp,1,1d-3,temp)
  temp(:)=xred_ideal(:,iatom)+temp(:)
  call DGEMV('T',3,3,1.d0,InVar%multiplicity(:,:),3,temp(:),1,0.d0,tmp(:,katom),1)
  katom_unitcell=mod(katom-1,InVar%natom_unitcell)+1
  vectk(:)=nint(tmp(:,katom)-tmp_store(:,katom_unitcell))

! To understand the meaning of "latt", see SearchS_1at
  ok=.false.
  do isym=1,Sym%nsym
!   TODO : A CHECKER!!!!!!!!!!!!!!!!!  
    indsym2(:)=zero
    call tdep_calc_indsym2(InVar,eatom,fatom,indsym2,isym,Sym,xred_ideal)
    ee=indsym2(4)
    lattef(:)=indsym2(1:3)
    lattfe(:)=indsym2(5:7)

    indsym2(:)=zero
    call tdep_calc_indsym2(InVar,fatom,gatom,indsym2,isym,Sym,xred_ideal)
    ff=indsym2(4)
    lattfg(:)=indsym2(1:3)
    lattgf(:)=indsym2(5:7)

    indsym2(:)=zero
    call tdep_calc_indsym2(InVar,gatom,eatom,indsym2,isym,Sym,xred_ideal)
    gg=indsym2(4)
    lattge(:)=indsym2(1:3)
    latteg(:)=indsym2(5:7)
    if (ee==iatom_unitcell.and.ff==jatom_unitcell.and.gg==katom_unitcell) then
!FB      if (InVar%debug) then
!FB        write(InVar%stdout,*) 'For ee=iatom, ff=jatom and gg=katom, isym=',isym
!FB        write(InVar%stdout,'(4(a,i4))') ' + ',vecti(1),' - ',vectj(1),' - ',lattef(1),' + ',lattfe(1)
!FB        write(InVar%stdout,'(4(a,i4))') ' + ',vecti(2),' - ',vectj(2),' - ',lattef(2),' + ',lattfe(2)
!FB        write(InVar%stdout,'(4(a,i4))') ' + ',vecti(3),' - ',vectj(3),' - ',lattef(3),' + ',lattfe(3)
!FB        write(InVar%stdout,*) ' '
!FB        write(InVar%stdout,'(4(a,i4))') ' + ',vectj(1),' - ',vectk(1),' - ',lattfg(1),' + ',lattgf(1)
!FB        write(InVar%stdout,'(4(a,i4))') ' + ',vectj(2),' - ',vectk(2),' - ',lattfg(2),' + ',lattgf(2)
!FB        write(InVar%stdout,'(4(a,i4))') ' + ',vectj(3),' - ',vectk(3),' - ',lattfg(3),' + ',lattgf(3)
!FB        write(InVar%stdout,*) ' '
!FB        write(InVar%stdout,'(4(a,i4))') ' + ',vectk(1),' - ',vecti(1),' - ',lattge(1),' + ',latteg(1)
!FB        write(InVar%stdout,'(4(a,i4))') ' + ',vectk(2),' - ',vecti(2),' - ',lattge(2),' + ',latteg(2)
!FB        write(InVar%stdout,'(4(a,i4))') ' + ',vectk(3),' - ',vecti(3),' - ',lattge(3),' + ',latteg(3)
!FB        write(InVar%stdout,*) ' '
!FB      end if
      if ((sum(abs((vecti(:)-lattef(:))-(vectj(:)-lattfe(:)))).lt.tol8).and.&
&         (sum(abs((vectj(:)-lattfg(:))-(vectk(:)-lattgf(:)))).lt.tol8).and.&
&         (sum(abs((vectk(:)-lattge(:))-(vecti(:)-latteg(:)))).lt.tol8)) then
        Isym3at(1)=isym
        Isym3at(2)=1
        ok=.true.
      end if
    end if
    if (ok) exit
  end do

 end subroutine tdep_SearchS_3at

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine tdep_calc_indsym2(InVar,iatom,jatom,indsym2,isym,Sym,xred_ideal)

  implicit none

  integer, intent(in) :: iatom,jatom
  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  integer, intent(out) :: indsym2(8)
  double precision, intent(in) :: xred_ideal(3,InVar%natom)

  integer :: isym,mu,ii
  integer :: iatom_unitcell,jatom_unitcell
  integer :: vecti(3),vectj(3),vectsym(4,Sym%nptsym)
  double precision :: tmpi(3,InVar%natom),tmpj(3,InVar%natom),temp3(3,1),tmp_store(3,InVar%natom_unitcell)

! Store the positions of the atoms in the motif
  do ii=1,InVar%natom_unitcell
    call DGEMV('T',3,3,1.d0,InVar%multiplicity(:,:),3,xred_ideal(:,ii),1,0.d0,tmp_store(:,ii),1)
  end do

! Search the matrix transformation going from (k,l) to (i,j)
  tmpi(:,:)=0.d0
  tmpj(:,:)=0.d0
! For a single iatom
  call DGEMV('T',3,3,1.d0,InVar%multiplicity(:,:),3,xred_ideal(:,iatom),1,0.d0,tmpi(:,iatom),1)
  iatom_unitcell=mod(iatom-1,InVar%natom_unitcell)+1
  vecti(:)=nint(tmpi(:,iatom)-tmp_store(:,iatom_unitcell))
  vectsym(:,:)=0
  do mu=1,3 ! Apply inverse transformation to original coordinates. Note transpose of symrec.
    vectsym(mu,isym) = Sym%symrec(1,mu,isym)*vecti(1)+Sym%symrec(2,mu,isym)*vecti(2)+Sym%symrec(3,mu,isym)*vecti(3)
  end do
  indsym2(1:4)=Sym%indsym(1:4,isym,iatom_unitcell)+vectsym(1:4,isym)

! For a couple of (iatom,jatom). The (iatom,jatom) vector depends on the position of iatom (due to PBC)
!FB  if (InVar%debug) write(InVar%stdout,*) '=========================================='
!FB  if (InVar%debug) write(InVar%stdout,'(a,i4,a,3(f10.5,x))') 'For iatom=',iatom,' with xred=',xred_ideal(:,iatom)
!FB  if (InVar%debug) write(InVar%stdout,*) '  =========================================='
!FB  if (InVar%debug) write(InVar%stdout,'(a,i4,a,3(f10.5,x))') '  For jatom=',jatom,' with xred=',xred_ideal(:,jatom)
  temp3(:,1)=xred_ideal(:,jatom)-xred_ideal(:,iatom)
  call tdep_make_inbox(temp3(:,1),1,1d-3,temp3(:,1))
  temp3(:,1)=xred_ideal(:,iatom)+temp3(:,1)
  call DGEMV('T',3,3,1.d0,InVar%multiplicity(:,:),3,temp3(:,1),1,0.d0,tmpj(:,jatom),1)
  jatom_unitcell=mod(jatom-1,InVar%natom_unitcell)+1
  vectj(:)=nint(tmpj(:,jatom)-tmp_store(:,jatom_unitcell))
  vectsym(:,:)=0
  do mu=1,3 ! Apply inverse transformation to original coordinates. Note transpose of symrec.
    vectsym(mu,isym) = Sym%symrec(1,mu,isym)*vectj(1)+Sym%symrec(2,mu,isym)*vectj(2)+Sym%symrec(3,mu,isym)*vectj(3)
  end do
  indsym2(5:8)=Sym%indsym(1:4,isym,jatom_unitcell)+vectsym(1:4,isym)
!FB  if (InVar%debug) then
!FB    write(InVar%stdout,'(a,i2,a,i4,a,i4,a,3(i4,x),a,i2,a,3(i4,x),a,i2,a)') '  indsym2(isym=',isym,',',iatom,',',jatom,')=',indsym2(1:3),&
!FB&     '|iat=',indsym2(4),'|',indsym2(5:7),'|iat=',indsym2(8),'|'
!FB  end if

 end subroutine tdep_calc_indsym2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module m_tdep_sym
