!{\src2tex{textfont=tt}}
!!****p* ABINIT/band2eps
!! NAME
!! band2eps
!!
!! FUNCTION
!! Draws the phonon dispersion curve in Encapsuled PostScript (EPS)
!! in black and white or in color according to the displacement participation
!! of each atom.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (FDortu,MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


program band2eps

 use defs_basis
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'band2eps'
!End of the abilint section

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
!no_abirules
 character(len=fnlen) :: filnam(4)
 real(dp) :: E,Emax,Emin,deltaE
 integer :: EmaxN,EminN,gradRes,kmaxN,kminN,lastPos,pos,posk
 integer :: iatom,ii,imode,iqpt,jj,natom,nqpt
 real(dp),allocatable :: phfrq(:),phfrqqm1(:)
 real(dp),allocatable :: color(:,:)
 real(dp) :: facUnit,norm,renorm
 real(dp),allocatable :: scale(:)
 character(len=6),allocatable :: qname(:)
 integer,allocatable :: nqptl(:)
 real(dp),allocatable :: colorAtom(:,:)
 real(dp),allocatable :: displ(:,:)
  !scale : hold the scale for each line (dimension=nlines)
  !qname : hold the name (gamma,R,etc..) for each extremity of line (dimension=nlines+1)
  !nqptl : =nqpt by line (dimension=nlines)
 integer :: cunits,nlines
  !nlines : number of lines
  !Emin is the minimum energy of the vertical axe
  !Emax is the maximum energy of the vertical axe
  !EminN is the minimum value of the vertical axe(in point)
  !EmaxN is the maximum value of the vertical axe(in point)
  !kminN is the minimum value of the horizontal axe(in point)
  !kmaxN is the maximum value of the horizontal axe(in point)
  !E,deltaE,pos are a work variables
  !gradRes is the number of intervals for the graduation along vertical axe

! *********************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 call xmpi_init()

!Initialize memory profiling if it is activated
!if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

!read the .file file
!File names refer to following files, in order:
!(1) Formatted input file
!(2) EPS graphic
!(3) Input phonon energies (from sortph.f)
!(4) Input displacements (from sortph.f)
 read(std_in, '(a)' ) filnam(1)
 if(len_trim(filnam(1))>132)then
   write(std_out,*)' band2eps : WARNING -'
   write(std_out,*)'  The length of this filename exceeds 132 characters. This might induce trouble later.'
 end if
 read(std_in, '(a)' ) filnam(2)
 if(len_trim(filnam(1))>132)then
   write(std_out,*)' band2eps : WARNING -'
   write(std_out,*)'  The length of this filename exceeds 132 characters. This might induce trouble later.'
 end if
 read(std_in, '(a)' ) filnam(3)
 if(len_trim(filnam(1))>132)then
   write(std_out,*)' band2eps : WARNING -'
   write(std_out,*)'  The length of this filename exceeds 132 characters. This might induce trouble later.'
 end if
 read(std_in, '(a)' ) filnam(4)
 if(len_trim(filnam(1))>132)then
   write(std_out,*)' band2eps : WARNING -'
   write(std_out,*)'  The length of this filename exceeds 132 characters. This might induce trouble later.'
 end if
!open(18,FILE='band.eps',STATUS='replace',ACCESS='sequential',ACTION='write')
!Open the '.eps' file for write
 write(std_out,'(a,a)') 'Opening file ', filnam(2)
 open(18, file=filnam(2))
!Open the phonon energies file
 write(std_out,'(a,a)') 'Opening file ', filnam(3)
 open(19, file=filnam(3))
 if(filnam(4)/='no') then
!  Open the displacements file
   write(std_out,'(a,a)') 'Opening file ', filnam(4)
   open(20, file=filnam(4))
 end if
!Open the input file
 write(std_out,'(a,a)') 'Opening file ', filnam(1)
 open(21, file=filnam(1))
!Boundings of the plot (only the plot and not what is around)
 EminN=6900
 EmaxN=2400
 kminN=2400
 kmaxN=9600
!Read input file (input.band) to know how to format the graph
 read(21,*)
 read(21,*)
 read(21,*) natom
 read(21,*)
 read(21,*) Emin,Emax,gradRes
 read(21,*)
 read(21,*) cunits
 read(21,*)
 read(21,*) nlines
!Allocate dynamique variables
 ABI_ALLOCATE(phfrqqm1,(3*natom))
 ABI_ALLOCATE(phfrq,(3*natom))
 ABI_ALLOCATE(color,(3,3*natom))
 ABI_ALLOCATE(qname,(nlines+1))
 ABI_ALLOCATE(scale,(nlines))
 ABI_ALLOCATE(nqptl,(nlines))
 ABI_ALLOCATE(colorAtom,(3,natom))
!colorAtom(1,1:5) : atoms contributing to red (ex : [1 0 0 0 0])
!colorAtom(2,1:5) : atoms contributing to green (ex : [0 1 0 0 0])
!colorAtom(3,1:5) : atoms contributing to blue (ex : [0 0 1 1 1])
 ABI_ALLOCATE(displ,(natom,3*natom))
!Read end of input file
 read(21,*)
 read(21,*) (qname(ii),ii=1,nlines+1)
 read(21,*)
 read(21,*) (nqptl(ii),ii=1,nlines)
 read(21,*)
 read(21,*) (scale(ii),ii=1,nlines)
 read(21,*)
 read(21,*)
 read(21,*)
 read(21,*) (colorAtom(1,ii),ii=1,natom)
 read(21,*)
 read(21,*) (colorAtom(2,ii),ii=1,natom)
 read(21,*)
 read(21,*) (colorAtom(3,ii),ii=1,natom)
!Multiplication factor for units (from Hartree to cm-1 or THz)
 if(cunits==1) then
   facUnit=Ha_cmm1
 elseif(cunits==2) then
   facUnit=Ha_THz
 else
 end if
!calculate nqpt
 nqpt=0
 do ii=1,nlines
   nqpt=nqpt+nqptl(ii)
 end do
!compute normalisation factor
 renorm=0
 do ii=1,nlines
   renorm=renorm+nqptl(ii)*scale(ii)
 end do
 renorm=renorm/nqpt
!Calculate Emin and Emax
 Emin=Emin/FacUnit
 Emax=Emax/FacUnit

!*******************************************************
!Begin to write some comments in the eps file
!This is based to 'xfig'

 write(18,'(a)') '% !PS-Adobe-2.0 EPSF-2.0'
 write(18,'(a)') '%%Title: band.ps'
 write(18,'(a)') '%%BoundingBox: 0 0 581 310'
 write(18,'(a)') '%%Magnification: 1.0000'

 write(18,'(a)') '/$F2psDict 200 dict def'
 write(18,'(a)') '$F2psDict begin'
 write(18,'(a)') '$F2psDict /mtrx matrix put'
 write(18,'(a)') '/col-1 {0 setgray} bind def'
 write(18,'(a)') '/col0 {0.000 0.000 0.000 srgb} bind def'
 write(18,'(a)') 'end'
 write(18,'(a)') 'save'
 write(18,'(a)') 'newpath 0 310 moveto 0 0 lineto 581 0 lineto 581 310 lineto closepath clip newpath'
 write(18,'(a)') '-36.0 446.0 translate'
 write(18,'(a)') '1 -1 scale'

 write(18,'(a)') '/cp {closepath} bind def'
 write(18,'(a)') '/ef {eofill} bind def'
 write(18,'(a)') '/gr {grestore} bind def'
 write(18,'(a)') '/gs {gsave} bind def'
 write(18,'(a)') '/sa {save} bind def'
 write(18,'(a)') '/rs {restore} bind def'
 write(18,'(a)') '/l {lineto} bind def'
 write(18,'(a)') '/m {moveto} bind def'
 write(18,'(a)') '/rm {rmoveto} bind def'
 write(18,'(a)') '/n {newpath} bind def'
 write(18,'(a)') '/s {stroke} bind def'
 write(18,'(a)') '/sh {show} bind def'
 write(18,'(a)') '/slc {setlinecap} bind def'
 write(18,'(a)') '/slj {setlinejoin} bind def'
 write(18,'(a)') '/slw {setlinewidth} bind def'
 write(18,'(a)') '/srgb {setrgbcolor} bind def'
 write(18,'(a)') '/rot {rotate} bind def'
 write(18,'(a)') '/sc {scale} bind def'
 write(18,'(a)') '/sd {setdash} bind def'
 write(18,'(a)') '/ff {findfont} bind def'
 write(18,'(a)') '/sf {setfont} bind def'
 write(18,'(a)') '/scf {scalefont} bind def'
 write(18,'(a)') '/sw {stringwidth} bind def'
 write(18,'(a)') '/tr {translate} bind def'
 write(18,'(a)') '/tnt {dup dup currentrgbcolor'

 write(18,'(a)') '4 -2 roll dup 1 exch sub 3 -1 roll mul add'
 write(18,'(a)') '4 -2 roll dup 1 exch sub 3 -1 roll mul add'
 write(18,'(a)') '4 -2 roll dup 1 exch sub 3 -1 roll mul add srgb}'
 write(18,'(a)') 'bind def'
 write(18,'(a)') '/shd {dup dup currentrgbcolor 4 -2 roll mul 4 -2 roll mul'
 write(18,'(a)') ' 4 -2 roll mul srgb} bind def'
 write(18,'(a)') '/$F2psBegin {$F2psDict begin /$F2psEnteredState save def} def'
 write(18,'(a)') '/$F2psEnd {$F2psEnteredState restore end} def'
 write(18,'(a)') '$F2psBegin'
 write(18,'(a)') '%%Page: 1 1'
 write(18,'(a)') '10 setmiterlimit'
 write(18,'(a)') '0.06000 0.06000 sc'

!****************************************************************
!Begin of the intelligible part of the postcript document

 write(18,'(a)') '%**************************************'
!****************************************************************
!Draw the box containing the plot
 write(18,'(a)') '%****Big Box****'
 write(18,'(a)') '7.500 slw'
 write(18,'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)') 'n ', kminN,' ', EmaxN,&
& ' m ', kmaxN,' ', EmaxN, ' l ', &
& kmaxN,' ', EminN, ' l ', kminN,' ', EminN, ' l'
 write(18,'(a)') 'cp gs col0 s gr'

!****************************************************************
!Write unit on the middle left of the vertical axe
 write(18,'(a)') '%****Units****'



 if(cunits==1) then
!  1/lambda
   write(18,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(18,'(a)') '1425 5650 m'
   write(18,'(3a)') 'gs 1 -1 sc  90.0 rot (Frequency ',achar(92),'(cm) col0 sh gr'
!  cm-1
   write(18,'(a)') '/Times-Roman ff 200.00 scf sf'
   write(18,'(a)') '1325 4030 m'
   write(18,'(a)') 'gs 1 -1 sc 90.0 rot  (-1) col0 sh gr'
   write(18,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(18,'(a)') '1425 3850 m'
   write(18,'(3a)') 'gs 1 -1 sc  90.0 rot (',achar(92),')) col0 sh gr'
 else
!  Freq
   write(18,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(18,'(a)') '825 4850 m'
   write(18,'(a)') 'gs 1 -1 sc  90.0 rot (Freq) col0 sh gr'
!  THz
   write(18,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(18,'(a)') '825 4350 m'
   write(18,'(a)') 'gs 1 -1 sc 90.0 rot  (THz) col0 sh gr'
 end if
!*****************************************************************
!Write graduation on the vertical axe
 write(18,'(a)') '%****Vertical graduation****'
 deltaE=(Emax-Emin)/gradRes

!Replacing do loop with real variables with standard g95 do loop
 E=Emin
 do
!  do E=Emin,(Emax-deltaE/2),deltaE
   if (E >= (Emax-deltaE/2)-tol6) exit
   pos=int(((EminN-EmaxN)*E &
&   +EmaxN*Emin -EminN*Emax)/(Emin-Emax))

!  write the value of energy(or frequence)
   write(18,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(18,'(i4,a,i4,a)') kminN-800,' ',pos+60,' m'        !-1300 must be CHANGED
!  as a function of the width of E
   write(18,'(a,i6,a)') 'gs 1 -1 sc (', nint(E*facUnit),') col0 sh gr'

!  write a little bar
   write(18,'(a,i4,a,i4,a,i4,a,i4,a)') 'n ', kminN,' ',pos ,' m ', kminN+100,' ', pos, ' l'
   write(18,'(a)') 'gs col0 s gr '

   E = E+deltaE
 end do

!do the same thing for E=Emax (floating point error)
 write(18,'(a)') '/Times-Roman ff 270.00 scf sf'
 write(18,'(i4,a,i4,a)') kminN-800,' ',EmaxN+60,' m'        !-1300 must be changed as E
 write(18,'(a,i6,a)') 'gs 1 -1 sc (', nint(Emax*facUnit),') col0 sh gr'


!draw zero line
 E=0
 pos=int(((EminN-EmaxN)*E &
& +EmaxN*Emin -EminN*Emax)/(Emin-Emax))
 write(18,'(a,i4,a,i4,a,i4,a,i4,a)') 'n ', kminN,' ',pos ,' m ', kmaxN,' ', pos, ' l'
 write(18,'(a)') 'gs col0 s gr '


!******************************************************
!draw legend of horizontal axe
!+vertical line

 write(18,'(a)') '%****Horizontal graduation****'



 lastPos=kminN

 do ii=0,nlines

   if(ii/=0) then
     posk=int(((kminN-kmaxN)*(nqptl(ii))) &
&     *scale(ii)/renorm/(-nqpt))
   else
     posk=0
   end if

   posk=posk+lastPos
   lastPos=posk

   if(qname(ii+1)=='gamma') then             !GAMMA
     write(18,'(a)') '/Symbol ff 270.00 scf sf'
     write(18,'(i4,a,i4,a)') posk-100,' ', 7150, ' m'
     write(18,'(a)') 'gs 1 -1 sc (G) col0 sh gr'
   elseif(qname(ii+1)=='lambda') then              !LAMBDA
     write(18,'(a)') '/Symbol ff 270.00 scf sf'
     write(18,'(i4,a,i4,a)') posk-100,' ', 7150, ' m'
     write(18,'(a)') 'gs 1 -1 sc (L) col0 sh gr'
   else                                     !autre
     write(18,'(a)') '/Times-Roman ff 270.00 scf sf'
     write(18,'(i4,a,i4,a)') posk-100,' ', 7150, ' m'
     write(18,'(a,a1,a)') 'gs 1 -1 sc (',qname(ii+1),') col0 sh gr'
   end if


!  draw vertical line
   write(18,'(a,i4,a,i4,a,i4,a,i4,a)') 'n ', posk,' ',EminN ,' m ', posk,' ', EmaxN, ' l'
   write(18,'(a)') 'gs col0 s gr '


 end do




!***********************************************************
!Write the bands (the most important part actually)

 write(18,'(a)') '%****Write Bands****'

 read(19,*) (phfrqqm1(ii),ii=1,3*natom)

 lastPos=kminN
 do jj=1,nlines
   do iqpt=1,nqptl(jj)
     read(19,*) (phfrq(ii),ii=1,3*natom)
     do imode=1,3*natom


       if(filnam(4)/='no') then       !calculate the color else in black and white
         read(20,*) (displ(iatom,imode),iatom=1,natom)


!        normalize displ
         norm=0
         do iatom=1,natom
           norm=norm+displ(iatom,imode)
         end do

         do iatom=1,natom
           displ(iatom,imode)=displ(iatom,imode)/norm
         end do

!        Treat color
         color(:,imode)=0
         do ii=1,natom
!          Red
           color(1,imode)=color(1,imode)+displ(ii,imode)*colorAtom(1,ii)
!          Green
           color(2,imode)=color(2,imode)+displ(ii,imode)*colorAtom(2,ii)
!          Blue
           color(3,imode)=color(3,imode)+displ(ii,imode)*colorAtom(3,ii)
         end do

       end if

       pos=int(((EminN-EmaxN)*phfrqqm1(imode) &
&       +EmaxN*Emin -EminN*Emax)/(Emin-Emax))

       posk=int(((kminN-kmaxN)*(iqpt-1) &
&       *scale(jj)/renorm/(-nqpt)))
       posk=posk+lastPos


       write(18,'(a,i4,a,i4,a)') 'n ',posk,' ',pos,' m'


       pos=int(((EminN-EmaxN)*phfrq(imode) &
&       +EmaxN*Emin -EminN*Emax)/(Emin-Emax))
       posk=int(((kminN-kmaxN)*(iqpt) &
&       *scale(jj)/renorm/(-nqpt)))
       posk=posk+lastPos
       write(18,'(i4,a,i4,a)') posk,' ',pos,' l gs'



       if(filnam(4)/='no') then     !(in color)
         write(18,'(f6.3,a,f6.3,a,f6.3,a)') color(1,imode),' ', &
&         color(2,imode),' ',color(3,imode), ' srgb s gr'
       else
         write(18,'(f6.3,a,f6.3,a,f6.3,a)') 0.0,' ', &
&         0.0,' ',0.0, ' srgb s gr'
       end if


     end do


     phfrqqm1=phfrq

   end do
   lastPos=posk

 end do


!**********************************************************
!Ending the poscript document
 write(18,'(a)') '$F2psEnd'
 write(18,'(a)') 'rs'

 !call abinit_doctor("__band2eps")

 call xmpi_end()

 end program band2eps
!!***
