!{\src2tex{textfont=tt}}
!!****f* ABINIT/ddb_io_out
!!
!! NAME
!! ddb_io_out
!!
!! FUNCTION
!! Open Derivative DataBase, then
!! write Derivative DataBase preliminary information.
!! Note: only one processor writes the DDB.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! amu(mtypat)=mass of the atoms (atomic mass unit)
!! dilatmx=the maximal dilatation factor
!! character(len=fnlen) dscrpt:string that describe the output database
!! ecut=kinetic energy planewave cutoff (hartree)
!! ecutsm=smearing energy for plane wave kinetic energy (Ha)
!! character(len=fnlen) filnam: name of output file
!! intxc=control xc quadrature
!! iscf=parameter controlling scf or non-scf choice
!! ixc=exchange-correlation choice parameter
!! kpt(3,mkpt)=k point set (reduced coordinates)
!! kptnrm=normalisation of k points
!! matom=maximum number of atoms
!! mband=maximum number of bands
!! mkpt=maximum number of special points
!! msym=maximum number of symetries
!! mtypat=maximum number of atom types
!! natom=number of atoms in the unit cell
!! nband(mkpt)=number of bands at each k point, for each polarization
!! ngfft(18)=contain all needed information about 3D FFT,
!!        see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nkpt=number of k points
!! nspden=number of spin-density components
!! nspinor=number of spinorial components of the wavefunctions
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! nsym=number of symmetry elements in space group
!! ntypat=number of atom types
!! occ(mband*mkpt)=occupation number for each band and k
!! occopt=option for occupancies
!! pawecutdg=cut-off for fine "double grid" used in PAW calculations (unused for NCPP)
!! rprim(3,3)=dimensionless primitive translations in real space
!! dfpt_sciss=scissor shift (Ha)
!! spinat(3,matom)=initial spin of each atom, in unit of hbar/2
!! symafm(msym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,msym)=symmetry operations in real space
!! tnons(3,msym)=nonsymmorphic translations for symmetry operations
!! tolwfr=tolerance on largest wf residual
!! tphysel="physical" electronic temperature with FD occupations
!! tsmear=smearing width (or temperature) in Hartree
!! typat(matom)=type of each atom
!! unddb=unit number for output
!! usepaw=flag for PAW
!! vrsddb=6 digit integer giving date, in form yymmdd for month=mm(1-12),
!!  day=dd(1-31), and year=yy(90-99 for 1990 to 1999,00-89 for 2000 to 2089),
!!  of current DDB version.
!! wtk(mkpt)=weight assigned to each k point
!! xred(3,matom)=reduced atomic coordinates
!! zion(mtypat)=valence charge of each type of atom
!! znucl(mtypat)=atomic number of atom type
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      dfpt_looppert,eig2tot,gstate,mblktyp1,mblktyp5,nonlinear,respfn
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ddb_io_out (dscrpt,filnam,matom,mband,&
&  mkpt,msym,mtypat,unddb,vrsddb,&
&  acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
&  natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
&  pawecutdg,rprim,dfpt_sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,&
&  typat,usepaw,wtk,xred,zion,znucl)


 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_io_tools,     only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_io_out'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: matom,mband,mkpt,msym,mtypat,unddb,vrsddb
 integer,intent(in) :: intxc,iscf,ixc,natom,nkpt,nspden,nspinor,nsppol,nsym
 integer,intent(in) :: ntypat,occopt,usepaw
 real(dp),intent(in) :: dilatmx,ecut,ecutsm,kptnrm,pawecutdg,dfpt_sciss,tolwfr,tphysel
 real(dp),intent(in) :: tsmear
 character(len=fnlen),intent(in) :: dscrpt,filnam
!arrays
 integer,intent(in) :: nband(mkpt),ngfft(18),symafm(msym),symrel(3,3,msym)
 integer,intent(in) :: typat(matom)
 real(dp),intent(in) :: acell(3),amu(mtypat),kpt(3,mkpt),occ(mband*mkpt)
 real(dp),intent(in) :: rprim(3,3),spinat(3,matom),tnons(3,msym),wtk(mkpt)
 real(dp),intent(in) :: xred(3,matom),zion(mtypat),znucl(mtypat)

!Local variables -------------------------
!Set routine version number here:
!scalars
 integer,parameter :: vrsio8=100401,vrsio8_old=010929,vrsio8_old_old=990527
 integer :: bantot,ii,ij,ikpt,iline,im
 character(len=500) :: message
!arrays
 character(len=9) :: name(9)

! *********************************************************************

 DBG_ENTER("COLL")

!Check ioddb8 version number (vrsio8) against mkddb version number
!(vrsddb)
 if (vrsio8/=vrsddb) then
   write(message, '(a,a,a,i10,a,a,i10,a)' )&
&   ' ddb_io_out: WARNING -',ch10,&
&   '  The input/output DDB version number=',vrsio8,ch10,&
&   '  is not equal to the DDB version number=',vrsddb,'.'
   call wrtout(std_out,message,'COLL')
 end if

!Open the output derivative database.
!(version 2.1. : changed because of a bug in a Perl script
!should set up a name checking procedure, with change of name
!like for the output file)
 if (open_file(filnam,message,unit=unddb,status='unknown',form='formatted') /= 0) then
   MSG_ERROR(message)
 end if

!Write the heading
 write(unddb, '(/,a,/,a,i10,/,/,a,a,/)' ) &
& ' **** DERIVATIVE DATABASE ****    ',&
& '+DDB, Version number',vrsddb,' ',trim(dscrpt)

!Write the descriptive data
!1. usepaw
 write(unddb, '(1x,a9,i10)' )'   usepaw',usepaw
!2. natom
 write(unddb, '(1x,a9,i10)' )'    natom',natom
!3. nkpt
 write(unddb, '(1x,a9,i10)' )'     nkpt',nkpt
!4. nsppol
 write(unddb, '(1x,a9,i10)' )'   nsppol',nsppol
!5. nsym
 write(unddb, '(1x,a9,i10)' )'     nsym',nsym
!6. ntypat
 write(unddb, '(1x,a9,i10)' )'   ntypat',ntypat
!7. occopt
 write(unddb, '(1x,a9,i10)' )'   occopt',occopt
!8. nband
 if(occopt==2)then
   im=12
   name(1)='    nband'
   do iline=1,(nkpt+11)/12
     if(iline==(nkpt+11)/12)im=nkpt-12*(iline-1)
     write(unddb, '(1x,a9,5x,12i5)' )name(1),&
&     (nband((iline-1)*12+ii),ii=1,im)
     name(1)='         '
   end do
   bantot=0
   do ikpt=1,nkpt
     bantot=bantot+nband(ikpt)
   end do
 else
   write(unddb, '(1x,a9,i10)' )'    nband',nband(1)
   bantot=nkpt*nband(1)
 end if

!9. acell
 write(unddb, '(1x,a9,3d22.14)' )'    acell',acell
!10. amu
 im=3
 name(1)='      amu'
 do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   write (unddb, '(1x,a9,3d22.14)' )name(1),&
&   (amu((iline-1)*3+ii),ii=1,im)
   name(1)='         '
 end do
!11. dilatmx
 write(unddb, '(1x,a9,d22.14)' )'  dilatmx',dilatmx
!12. ecut
 write(unddb, '(1x,a9,d22.14)' )'     ecut',ecut
!12b. pawecutdg (PAW)
 if (usepaw==1) then
   write(unddb, '(1x,a9,d22.14)' )'pawecutdg',pawecutdg
 end if
!13. ecutsm
 write(unddb, '(1x,a9,d22.14)' )'   ecutsm',ecutsm
!14. intxc
 write(unddb, '(1x,a9,i10)' )'    intxc',intxc
!15. iscf
 write(unddb, '(1x,a9,i10)' )'     iscf',iscf
!16. ixc
 write(unddb, '(1x,a9,i10)' )'      ixc',ixc
!17. kpt
 name(1)='      kpt'
 do iline=1,nkpt
   write (unddb, '(1x,a9,3d22.14)' )name(1),&
&   (kpt(ii,iline),ii=1,3)
   name(1)='      '
 end do
!18. kptnrm
 write(unddb, '(1x,a9,d22.14)' )'   kptnrm',kptnrm
!19. ngfft
 write(unddb, '(1x,a9,5x,3i5)' )'    ngfft',ngfft(1:3)
!20. nspden
 write(unddb, '(1x,a9,i10)' )'   nspden',nspden
!21. nspinor
 write(unddb, '(1x,a9,i10)' )'  nspinor',nspinor
!22. occ
 if(occopt==2)then
   im=3
   name(1)='      occ'
   do iline=1,(bantot+2)/3
     if(iline==(bantot+2)/3)im=bantot-3*(iline-1)
     write(unddb, '(1x,a9,3d22.14)' )name(1),&
&     (occ((iline-1)*3+ii),ii=1,im)
     name(1)='         '
   end do
 else
   im=3
   name(1)='      occ'
   do iline=1,(nband(1)+2)/3
     if(iline==(nband(1)+2)/3)im=nband(1)-3*(iline-1)
     write(unddb, '(1x,a9,3d22.14)' )name(1),&
&     (occ((iline-1)*3+ii),ii=1,im)
     name(1)='         '
   end do
 end if
!23. rprim
 name(1)='    rprim'
 do iline=1,3
   write(unddb, '(1x,a9,3d22.14)' )name(1),&
&   (rprim(ii,iline),ii=1,3)
   name(1)='      '
 end do
!24. dfpt_sciss
 write(unddb, '(1x,a11,d22.14)' )' dfpt_sciss',dfpt_sciss
!25. spinat
 name(1)='   spinat'
 do iline=1,natom
   write(unddb, '(1x,a9,3d22.14)' )name(1),&
&   (spinat(ii,iline),ii=1,3)
   name(1)='         '
 end do
!26. symafm
 im=12
 name(1)='   symafm'
 do iline=1,(nsym+11)/12
   if(iline==(nsym+11)/12)im=nsym-12*(iline-1)
   write(unddb, '(1x,a9,5x,12i5)' )name(1),&
&   (symafm((iline-1)*12+ii),ii=1,im)
   name(1)='         '
 end do
!27. symrel
 name(1)='   symrel'
 do iline=1,nsym
   write(unddb, '(1x,a9,5x,9i5)' )name(1),&
&   ((symrel(ii,ij,iline),ii=1,3),ij=1,3)
   name(1)='         '
 end do
!28. tnons
 name(1)='    tnons'
 do iline=1,nsym
   write(unddb, '(1x,a9,3d22.14)' )name(1),&
&   (tnons(ii,iline),ii=1,3)
   name(1)='         '
 end do
!29. tolwfr
 write(unddb, '(1x,a9,d22.14)' )'   tolwfr',tolwfr
!30. tphysel
 write(unddb, '(1x,a9,d22.14)' )'  tphysel',tphysel
!31. tsmear
 write(unddb, '(1x,a9,d22.14)' )'   tsmear',tsmear
!32. typat
 im=12
 name(1)='    typat'
 do iline=1,(natom+11)/12
   if(iline==(natom+11)/12)im=natom-12*(iline-1)
   write(unddb, '(1x,a9,5x,12i5)' )name(1),&
&   (typat((iline-1)*12+ii),ii=1,im)
   name(1)='         '
 end do
!33. wtk
 name(1)='      wtk'
 im=3
 do iline=1,(nkpt+2)/3
   if(iline==(nkpt+2)/3)im=nkpt-3*(iline-1)
   write(unddb, '(1x,a9,3d22.14)' )name(1),&
&   (wtk((iline-1)*3+ii),ii=1,im)
   name(1)='         '
 end do
!34. xred
 name(1)='     xred'
 do iline=1,natom
   write(unddb, '(1x,a9,3d22.14)' )name(1),&
&   (xred(ii,iline),ii=1,3)
   name(1)='         '
 end do
!35. znucl
 name(1)='    znucl'
 im=3
 do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   write(unddb, '(1x,a9,3d22.14)' )name(1),&
&   (znucl((iline-1)*3+ii),ii=1,im)
   name(1)='         '
 end do
!36. zion
 name(1)='     zion'
 im=3
 do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   write(unddb, '(1x,a9,3d22.14)' )name(1),&
&   (zion((iline-1)*3+ii),ii=1,im)
   name(1)='         '
 end do

 DBG_EXIT("COLL")

end subroutine ddb_io_out
!!***
