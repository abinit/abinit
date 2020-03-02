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
!! Copyright (C) 1999-2020 ABINIT group (FDortu,MVeithen)
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
!!      abi_io_redirect,abimem_init,band2eps_dtset_free,instrng,inupper
!!      invars11,outvars_band2eps,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


program band2eps

 use defs_basis
 use m_abimover
 use m_build_info
 use m_xmpi
 use m_abicore
 use m_errors
 use m_effective_potential
 use m_multibinit_dataset
 use m_effective_potential_file
 use m_band2eps_dataset

 use m_io_tools,      only : open_file
 use m_fstrings,      only : int2char4, tolower, inupper
 use m_time,          only : asctime
 use m_parser,        only : instrng

 implicit none

!Arguments -----------------------------------
!Local variables-------------------------------
 integer,parameter :: master=0
 character(len=fnlen) :: filnam(4)
 real(dp) :: E,deltaE
 integer :: comm,EmaxN,EminN,kmaxN,kminN,lastPos,lenstr,pos,posk
 integer :: iatom,ii,imode,io,iqpt,jj,nqpt
 integer :: nproc,my_rank
 integer :: option,unt1,unt2,unt3
 logical :: iam_master
!array
 real(dp),allocatable :: phfrq(:),phfrqqm1(:)
 real(dp),allocatable :: color(:,:)
 real(dp) :: facUnit,norm,renorm
 real(dp),allocatable :: colorAtom(:,:)
 real(dp),allocatable :: displ(:,:)
 type(band2eps_dataset_type) :: inp
 character(len=500) :: message
 character(len=strlen) :: string
  !scale : hold the scale for each line (dimension=nlines)
  !qname : hold the name (gamma,R,etc..) for each extremity of line (dimension=nlines+1)
  !nqptl : =nqpt by line (dimension=nlines)
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

!Initialize MPI
 call xmpi_init()
 comm = xmpi_world

!MPI variables
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

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
 write(std_out,*)' Give name for formatted input file : '
 read(std_in, '(a)',IOSTAT=io) filnam(1)
 write(std_out,'(a,a)' )'-   ',trim(filnam(1))

 write(std_out,*)' Give name for formatted output eps file : '
 read(std_in, '(a)',IOSTAT=io) filnam(2)
 write(std_out,'(a,a)' )'-   ',trim(filnam(2))

 write(std_out,*)' Give name for formatted phonon frequency file : '
 read(std_in, '(a)',IOSTAT=io) filnam(3)
 write(std_out,'(a,a)' )'-   ',trim(filnam(3))

 write(std_out,*)' Give name for formatted displacements file : '
 read(std_in, '(a)',IOSTAT=io) filnam(4)
 write(std_out,'(a,a)' )'-   ',trim(filnam(4))

!Read the input file, and store the information in a long string of characters
!strlen from defs_basis module
 write(std_out,'(a,a)') 'Opening and reading input file: ', filnam(1)
 option=1
 call instrng (filnam(1),lenstr,option,strlen,string)
 !To make case-insensitive, map characters to upper case:
 call inupper(string(1:lenstr))

!Read the input file
 call invars11(inp,lenstr,string)
 if(inp%prtout == 1) call outvars_band2eps(inp,std_out)

!Open the '.eps' file for write
 write(std_out,'(a,a)') 'Creation of file ', filnam(2)
 if (open_file(filnam(2),message,newunit=unt1,form="formatted",status="unknown",action="write") /= 0) then
   MSG_ERROR(message)
 end if
!Open the phonon energies file
 if (open_file(filnam(3),message,newunit=unt2,form="formatted") /= 0) then
   MSG_ERROR(message)
 end if
 if(filnam(4)/='no') then
!  Open the displacements file
   if (open_file(filnam(4),message,newunit=unt3,form="formatted",status="old",action='read') /= 0) then
     MSG_ERROR(message)
   end if
 end if


!Boundings of the plot (only the plot and not what is around)
 EminN=6900
 EmaxN=2400
 kminN=2400
 kmaxN=9600

!Allocate dynamique variables
 ABI_ALLOCATE(phfrqqm1,(3*inp%natom))
 ABI_ALLOCATE(phfrq,(3*inp%natom))
 ABI_ALLOCATE(color,(3,3*inp%natom))
 ABI_ALLOCATE(colorAtom,(3,inp%natom))
!colorAtom(1,1:5) : atoms contributing to red (ex : [1 0 0 0 0])
!colorAtom(2,1:5) : atoms contributing to green (ex : [0 1 0 0 0])
!colorAtom(3,1:5) : atoms contributing to blue (ex : [0 0 1 1 1])
!tranfert color from input
 colorAtom(1,:) = inp%red
 colorAtom(2,:) = inp%green
 colorAtom(3,:) = inp%blue
 ABI_ALLOCATE(displ,(inp%natom,3*inp%natom))
!Read end of input file

!Multiplication factor for units (from Hartree to cm-1 or THz)
 if(inp%cunit==1) then
   facUnit=Ha_cmm1
 elseif(inp%cunit==2) then
   facUnit=Ha_THz
 else
 end if
!calculate nqpt
 nqpt=0
 do ii=1,inp%nlines
   nqpt=nqpt+inp%nqline(ii)
 end do
!compute normalisation factor
 renorm=0
 do ii=1,inp%nlines
   renorm=renorm+inp%nqline(ii)*inp%scale(ii)
 end do
 renorm=renorm/nqpt
!Calculate inp%min and inp%max
 inp%min=inp%min/FacUnit
 inp%max=inp%max/FacUnit

!*******************************************************
!Begin to write some comments in the eps file
!This is based to 'xfig'

 write(unt1,'(a)') '% !PS-Adobe-2.0 EPSF-2.0'
 write(unt1,'(a)') '%%Title: band.ps'
 write(unt1,'(a)') '%%BoundingBox: 0 0 581 310'
 write(unt1,'(a)') '%%Magnification: 1.0000'

 write(unt1,'(a)') '/$F2psDict 200 dict def'
 write(unt1,'(a)') '$F2psDict begin'
 write(unt1,'(a)') '$F2psDict /mtrx matrix put'
 write(unt1,'(a)') '/col-1 {0 setgray} bind def'
 write(unt1,'(a)') '/col0 {0.000 0.000 0.000 srgb} bind def'
 write(unt1,'(a)') 'end'
 write(unt1,'(a)') 'save'
 write(unt1,'(a)') 'newpath 0 310 moveto 0 0 lineto 581 0 lineto 581 310 lineto closepath clip newpath'
 write(unt1,'(a)') '-36.0 446.0 translate'
 write(unt1,'(a)') '1 -1 scale'

 write(unt1,'(a)') '/cp {closepath} bind def'
 write(unt1,'(a)') '/ef {eofill} bind def'
 write(unt1,'(a)') '/gr {grestore} bind def'
 write(unt1,'(a)') '/gs {gsave} bind def'
 write(unt1,'(a)') '/sa {save} bind def'
 write(unt1,'(a)') '/rs {restore} bind def'
 write(unt1,'(a)') '/l {lineto} bind def'
 write(unt1,'(a)') '/m {moveto} bind def'
 write(unt1,'(a)') '/rm {rmoveto} bind def'
 write(unt1,'(a)') '/n {newpath} bind def'
 write(unt1,'(a)') '/s {stroke} bind def'
 write(unt1,'(a)') '/sh {show} bind def'
 write(unt1,'(a)') '/slc {setlinecap} bind def'
 write(unt1,'(a)') '/slj {setlinejoin} bind def'
 write(unt1,'(a)') '/slw {setlinewidth} bind def'
 write(unt1,'(a)') '/srgb {setrgbcolor} bind def'
 write(unt1,'(a)') '/rot {rotate} bind def'
 write(unt1,'(a)') '/sc {scale} bind def'
 write(unt1,'(a)') '/sd {setdash} bind def'
 write(unt1,'(a)') '/ff {findfont} bind def'
 write(unt1,'(a)') '/sf {setfont} bind def'
 write(unt1,'(a)') '/scf {scalefont} bind def'
 write(unt1,'(a)') '/sw {stringwidth} bind def'
 write(unt1,'(a)') '/tr {translate} bind def'
 write(unt1,'(a)') '/tnt {dup dup currentrgbcolor'

 write(unt1,'(a)') '4 -2 roll dup 1 exch sub 3 -1 roll mul add'
 write(unt1,'(a)') '4 -2 roll dup 1 exch sub 3 -1 roll mul add'
 write(unt1,'(a)') '4 -2 roll dup 1 exch sub 3 -1 roll mul add srgb}'
 write(unt1,'(a)') 'bind def'
 write(unt1,'(a)') '/shd {dup dup currentrgbcolor 4 -2 roll mul 4 -2 roll mul'
 write(unt1,'(a)') ' 4 -2 roll mul srgb} bind def'
 write(unt1,'(a)') '/$F2psBegin {$F2psDict begin /$F2psEnteredState save def} def'
 write(unt1,'(a)') '/$F2psEnd {$F2psEnteredState restore end} def'
 write(unt1,'(a)') '$F2psBegin'
 write(unt1,'(a)') '%%Page: 1 1'
 write(unt1,'(a)') '10 setmiterlimit'
 write(unt1,'(a)') '0.06000 0.06000 sc'

!****************************************************************
!Begin of the intelligible part of the postcript document

 write(unt1,'(a)') '%**************************************'
!****************************************************************
!Draw the box containing the plot
 write(unt1,'(a)') '%****Big Box****'
 write(unt1,'(a)') '12 slw'
 write(unt1,'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)') 'n ', kminN,' ', EmaxN,&
& ' m ', kmaxN,' ', EmaxN, ' l ', &
& kmaxN,' ', EminN, ' l ', kminN,' ', EminN, ' l'
 write(unt1,'(a)') 'cp gs col0 s gr'

!****************************************************************
!Write unit on the middle left of the vertical axe
 write(unt1,'(a)') '%****Units****'

 if(inp%cunit==1) then
!  1/lambda
   write(unt1,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(unt1,'(a)') '1425 5650 m'
   write(unt1,'(3a)') 'gs 1 -1 sc  90.0 rot (Frequency ',achar(92),'(cm) col0 sh gr'
!  cm-1
   write(unt1,'(a)') '/Times-Roman ff 200.00 scf sf'
   write(unt1,'(a)') '1325 4030 m'
   write(unt1,'(a)') 'gs 1 -1 sc 90.0 rot  (-1) col0 sh gr'
   write(unt1,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(unt1,'(a)') '1425 3850 m'
   write(unt1,'(3a)') 'gs 1 -1 sc  90.0 rot (',achar(92),')) col0 sh gr'
 else
!  Freq
   write(unt1,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(unt1,'(a)') '825 4850 m'
   write(unt1,'(a)') 'gs 1 -1 sc  90.0 rot (Freq) col0 sh gr'
!  THz
   write(unt1,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(unt1,'(a)') '825 4350 m'
   write(unt1,'(a)') 'gs 1 -1 sc 90.0 rot  (THz) col0 sh gr'
 end if
!*****************************************************************
!Write graduation on the vertical axe
 write(unt1,'(a)') '%****Vertical graduation****'
 deltaE=(inp%max-inp%min)/inp%ngrad

 E=inp%min
 do
!  do E=inp%min,(inp%max-deltaE/2),deltaE
   if (E >= (inp%max-deltaE/2)-tol6) exit
   pos=int(((EminN-EmaxN)*E &
&   +EmaxN*inp%min -EminN*inp%max)/(inp%min-inp%max))

!  write the value of energy(or frequence)
   write(unt1,'(a)') '/Times-Roman ff 270.00 scf sf'
   write(unt1,'(i4,a,i4,a)') kminN-800,' ',pos+60,' m'        !-1300 must be CHANGED
!  as a function of the width of E
   write(unt1,'(a,i6,a)') 'gs 1 -1 sc (', nint(E*facUnit),') col0 sh gr'

!  write a little bar
   write(unt1,'(a,i4,a,i4,a,i4,a,i4,a)') 'n ', kminN,' ',pos ,' m ', kminN+100,' ', pos, ' l'
   write(unt1,'(a)') 'gs col0 s gr '

   E = E+deltaE
 end do

!do the same thing for E=inp%max (floating point error)
 write(unt1,'(a)') '/Times-Roman ff 270.00 scf sf'
 write(unt1,'(i4,a,i4,a)') kminN-800,' ',EmaxN+60,' m'        !-1300 must be changed as E
 write(unt1,'(a,i6,a)') 'gs 1 -1 sc (', nint(inp%max*facUnit),') col0 sh gr'


!draw zero line
 E=0
 pos=int(((EminN-EmaxN)*E &
& +EmaxN*inp%min -EminN*inp%max)/(inp%min-inp%max))
 write(unt1,'(a,i4,a,i4,a,i4,a,i4,a)') 'n ', kminN,' ',pos ,' m ', kmaxN,' ', pos, ' l'
 write(unt1,'(a)') 'gs col0 s gr '


!******************************************************
!draw legend of horizontal axe
!+vertical line

 write(unt1,'(a)') '%****Horizontal graduation****'

 lastPos=kminN

 do ii=0,inp%nlines

   if(ii/=0) then
     posk=int(((kminN-kmaxN)*(inp%nqline(ii))) &
&     *inp%scale(ii)/renorm/(-nqpt))
   else
     posk=0
   end if

   posk=posk+lastPos
   lastPos=posk

   if(tolower(inp%qpoint_name(ii+1))=='gamma') then             !GAMMA
     write(unt1,'(a)') '/Symbol ff 270.00 scf sf'
     write(unt1,'(i4,a,i4,a)') posk-100,' ', 7150, ' m'
     write(unt1,'(a)') 'gs 1 -1 sc (G) col0 sh gr'
   elseif(tolower(inp%qpoint_name(ii+1))=='lambda') then              !LAMBDA
     write(unt1,'(a)') '/Symbol ff 270.00 scf sf'
     write(unt1,'(i4,a,i4,a)') posk-100,' ', 7150, ' m'
     write(unt1,'(a)') 'gs 1 -1 sc (L) col0 sh gr'
   else                                     !autre
     write(unt1,'(a)') '/Times-Roman ff 270.00 scf sf'
     write(unt1,'(i4,a,i4,a)') posk-100,' ', 7150, ' m'
     write(unt1,'(a,a1,a)') 'gs 1 -1 sc (',inp%qpoint_name(ii+1),') col0 sh gr'
   end if

!  draw vertical line
   write(unt1,'(a,i4,a,i4,a,i4,a,i4,a)') 'n ', posk,' ',EminN ,' m ', posk,' ', EmaxN, ' l'
   write(unt1,'(a)') 'gs col0 s gr '

 end do

!***********************************************************
!Write the bands (the most important part actually)

 write(unt1,'(a)') '%****Write Bands****'

 lastPos=kminN

 read(unt2,*) (phfrqqm1(ii),ii=1,3*inp%natom)

 do jj=1,inp%nlines
   do iqpt=1,inp%nqline(jj)
     read(unt2,*) (phfrq(ii),ii=1,3*inp%natom)
     do imode=1,3*inp%natom

       if(filnam(4)/='no') then       !calculate the color else in black and white
         do iatom=1,inp%natom
           read(unt3,*) displ(iatom,imode)
         end do
!        normalize displ
         norm=0
         do iatom=1,inp%natom
           norm=norm+displ(iatom,imode)
         end do

         do iatom=1,inp%natom
           displ(iatom,imode)=displ(iatom,imode)/norm
         end do

!        Treat color
         color(:,imode)=0
         do ii=1,inp%natom
!          Red
           color(1,imode)=color(1,imode)+displ(ii,imode)*colorAtom(1,ii)
!          Green
           color(2,imode)=color(2,imode)+displ(ii,imode)*colorAtom(2,ii)
!          Blue
           color(3,imode)=color(3,imode)+displ(ii,imode)*colorAtom(3,ii)
         end do
       end if

       pos=int(((EminN-EmaxN)*phfrqqm1(imode) &
&       +EmaxN*inp%min -EminN*inp%max)/(inp%min-inp%max))

       posk=int(((kminN-kmaxN)*(iqpt-1) &
&       *inp%scale(jj)/renorm/(-nqpt)))
       posk=posk+lastPos

       write(unt1,'(a,i5,a,i5,a)') 'n ',posk,' ',pos,' m'

       pos=int(((EminN-EmaxN)*phfrq(imode) &
&       +EmaxN*inp%min -EminN*inp%max)/(inp%min-inp%max))
       posk=int(((kminN-kmaxN)*(iqpt) &
&       *inp%scale(jj)/renorm/(-nqpt)))
       posk=posk+lastPos
       write(unt1,'(i5,a,i5,a)') posk,' ',pos,' l gs'


       if(filnam(4)/='no') then     !(in color)
         write(unt1,'(f6.3,a,f6.3,a,f6.3,a)') color(1,imode),' ', &
&         color(2,imode),' ',color(3,imode), ' srgb s gr'
       else
         write(unt1,'(f6.3,a,f6.3,a,f6.3,a)') 0.0,' ', &
&         0.0,' ',0.0, ' srgb s gr'
       end if


     end do

     phfrqqm1=phfrq

   end do
   lastPos=posk

 end do


!**********************************************************
!Ending the poscript document
 write(unt1,'(a)') '$F2psEnd'
 write(unt1,'(a)') 'rs'

 !call abinit_doctor("__band2eps")

 call band2eps_dtset_free(inp)

 call xmpi_end()

 end program band2eps
!!***
