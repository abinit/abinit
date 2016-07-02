program ebands

 implicit none

 !Variable declarations
 character(len=20) :: in_file,out_file,band_file,plt_file
 integer :: in,nband,kptoption
 real :: acell(3)
 real :: xrange
 character(len=2) :: band_dataset
 logical :: flag

 write(*,*)"Energy bands reader tool for ABINIT by Ersen Mete,Ph.D."
 write(*,*)"This utility generates a text formatted bands file and a&
 & raw gnuplot script."
 write(*,*)"Enter the ABINIT output file name :"
 read(*,*)out_file
 out_file=trim(out_file)
 call check_out(out_file,nband,band_dataset,kptoption)
 if (index(out_file,".out").gt.0) then
  band_file=out_file(1:index(out_file,".out"))//"bnd"
 else
  band_file=out_file//".bnd"
 end if
 band_file=trim(band_file)
 if (index(out_file,".out").gt.0) then
  plt_file=out_file(1:index(out_file,".out"))//"plt"
 else
  plt_file=out_file//".plt"
 end if
 inquire(file=band_file,exist=flag)
 if(flag)then
  in=65
  do while(flag)
   if(in.gt.90)then
    write(*,*)"Error : cannot assign a name for output files"
    stop
   end if
   if (index(out_file,".out").gt.0) then
    band_file=out_file(1:index(out_file,".out")-1)//char(in)//".bnd"
   else
    band_file=out_file//char(in)//".bnd"
   end if
   band_file=trim(band_file)
   if (index(out_file,".out").gt.0) then
    plt_file=out_file(1:index(out_file,".out")-1)//char(in)//".plt"
   else
    plt_file=out_file//char(in)//".plt"
   end if
   inquire(file=band_file,exist=flag)
   in=in+1
  end do
 end if
 call process_ebands(out_file,band_file,band_dataset,xrange)
 call write_plt(band_file,plt_file,nband,xrange,kptoption)
 write(*,*)"Output successfully generated!"
 write(*,*)"You may need to make modifications in the gnuplot script file."

contains

 subroutine check_out(out_file,nband,band_dataset,kptoption)

  implicit none

  character(len=20),intent(in) :: out_file
  integer,intent(inout) :: nband
  character(len=2),intent(inout) :: band_dataset
  integer,intent(inout) :: kptoption

  character(len=80) :: buffer
  character(len=2) :: ndtset,iscf,kptopt
  integer :: i,ios

  ! Check if the file exists
  open(5,file=out_file,status="old",iostat=ios)
  if (ios/=0) then
   write(*,*)"Error opening file : ",out_file
   stop
  end if
  ! Is this an ABINIT output?
  read(5,'(A80))')buffer
  do while(index(buffer,"ABINIT").lt.1.and.ios==0)
   read(5,'(A80)',iostat=ios)buffer
  end do
  if (index(buffer,"ABINIT").lt.1) then
   write(*,*)"This does not seem to be an ABINIT output file!"
   stop
  end if
  ! Check how many datasets exist
  do while(index(buffer,"ndtset").lt.1)
   read(5,'(A80)')buffer
  end do
  buffer=buffer(6+index(buffer,"ndtset"):80)
  ndtset=adjustl(buffer)
  ! Check which dataset, bands are computed in
  rewind(5)
  do while(index(buffer,"-outvars").lt.1)
   read(5,'(A80)')buffer
  end do
  do while(index(buffer,"iscf").lt.1)
   read(5,'(A80)')buffer
  end do
  do while(index(buffer,"iscf").gt.0.and.index(buffer,"-2").lt.1)
   read(5,'(A80)')buffer
  end do
  buffer=buffer(4+index(buffer,"iscf"):80)
  iscf=adjustl(buffer)
  rewind(5)
  do while(index(buffer,"-outvars").lt.1)
   read(5,'(A80)')buffer
  end do
  do while(index(buffer,"kptopt").lt.1)
   read(5,'(A80)')buffer
  end do
  do while(index(buffer,"kptopt").gt.0.and.index(buffer,"-").lt.1)
   read(5,'(A80)')buffer
  end do
  buffer=buffer(6+index(buffer,"kptopt"):80)
  kptopt=adjustl(buffer)
  if (iscf/=kptopt) then
   write(*,*)"This does not seem to be a proper band calculation output!"
   stop
  end if
  band_dataset=kptopt
  backspace(5)
  read(5,'(a12,i2)')buffer,kptoption
  ! Extract nband
  do while(index(buffer,"nband"//band_dataset).lt.1)
   read(5,'(A80)')buffer
  end do
  backspace(5)
  read(5,'(A13,I5)')buffer,nband
  close(5)

 end subroutine check_out

 subroutine process_ebands(out_file,band_file,band_dataset,xrange)

  implicit none

  ! Arguments
  character(len=20),intent(in) :: out_file,band_file
  character(len=2),intent(in) :: band_dataset
  real,intent(inout) :: xrange

  character(len=2) :: appen
  character(len=80) :: buffer
  real,allocatable :: kpts(:,:)
  real :: acell(3)
  real :: G1(3),G2(3),G3(3)
  character(len=*), parameter :: format01150a="(1x,a9,a,1x,(t13,3es16.8))" 
  character(len=*), parameter :: format01155a="(1x,a9,a,1x,(t13,i10))"
  integer :: i,ios,nkpt
  real :: x,kx,ky,kz,old_increment,increment
  real :: pi=3.1415926535897932384626433832795

  ! Open files for read and write
  open(5,file=out_file,status="old",iostat=ios)
  if (ios/=0) then
   write(*,*)"Error opening file : ",out_file
   stop
  end if
  open(8,file=band_file,status="unknown",iostat=ios)
  if (ios/=0) then
   write(*,*)"Error opening file : ",band_file
   stop
  end if
  ! Get nkpt
  do while(index(buffer,trim("nkpt"//band_dataset)).lt.1)
   read(5,'(a80)')buffer
  end do
  backspace(5)
  read(5,format01155a)buffer,appen,nkpt
  ! Get k-points
  rewind(5)
  kpts(:,:)=0.0d0
  do while(index(buffer,trim("kpt"//band_dataset)).lt.1)
   read(5,'(A80)')buffer
  end do
  backspace(5)
  allocate(kpts(3,nkpt))
  do i=1,nkpt
   read(5,format01150a)buffer,appen,kpts(1,i),kpts(2,i),kpts(3,i)
  end do
  ! Get acell
  rewind(5)
  acell(:)=1
  do while(index(buffer,"acell").lt.1)
   read(5,'(A80)')buffer
  end do
  backspace(5)
  read(5,'(1x,a10,1x,(t13,3es18.10))')buffer,acell(1),acell(2),acell(3)
  write(*,*)'Cell dimensions :'
  write(*,*)'acell = ',acell(1),acell(2),acell(3)
  ! Get reciprocal primitives
  rewind(5)
  do while(index(buffer,"R(1)").lt.1)
   read(5,'(A80)')buffer
  end do
  backspace(5)
  read(5,'(a6,3f11.7,a7,3f11.7)')buffer,G1(1),G1(2),G1(3),buffer,G1(1),G1(2),G1(3)
  read(5,'(a6,3f11.7,a7,3f11.7)')buffer,G2(1),G2(2),G2(3),buffer,G2(1),G2(2),G2(3)
  read(5,'(a6,3f11.7,a7,3f11.7)')buffer,G3(1),G3(2),G3(3),buffer,G3(1),G3(2),G3(3)
  write(*,*)'Reciprocal vectors (to be multiplied by 2*pi):'
  write(*,*)'G1= ',G1(1),G1(2),G1(3)
  write(*,*)'G2= ',G2(1),G2(2),G2(3)
  write(*,*)'G3= ',G3(1),G3(2),G3(3)
  do i=1,3
   G1(i)=G1(i)*2*pi
   G2(i)=G2(i)*2*pi
   G3(i)=G3(i)*2*pi
  end do
  ! Convert kpoints to cartesian coordinates
  do i=1,nkpt
   kx=kpts(1,i)*G1(1)+kpts(2,i)*G2(1)+kpts(3,i)*G3(1)
   ky=kpts(1,i)*G1(2)+kpts(2,i)*G2(2)+kpts(3,i)*G3(2)
   kz=kpts(1,i)*G1(3)+kpts(2,i)*G2(3)+kpts(3,i)*G3(3)
   kpts(1,i)=kx
   kpts(2,i)=ky
   kpts(3,i)=kz
  end do
  ! Read & write bands
  do while(.not.(index(buffer,"DATASET").gt.0.and.index(buffer,band_dataset).gt.0))
   read(5,'(a80)')buffer
  end do
  x=0
  do while(index(buffer,"kpt# ").lt.1)
   read(5,'(a80)')buffer
  end do
  old_increment=1.0
  do i=1,nkpt
   buffer=''
   write(8,'(i3,1x,f9.6,1x)',advance='no')i,x
   do while(index(buffer,"kpt# ").lt.1.and.index(buffer,"Min").lt.1)
    buffer=''
    read(5,'(a)')buffer
    buffer=trim(buffer)
    if (index(buffer,"kpt# ").lt.1) write(8,'(a)',advance='no')buffer(1:len_trim(buffer))
   end do
   write(8,'(a)')''
   increment=sqrt((kpts(1,i)-kpts(1,i+1))**2+(kpts(2,i)-kpts(2,i+1))**2+(kpts(3,i)-kpts(3,i+1))**2)
   if (7*old_increment>increment) then 
    x=x+increment
   else 
    increment=sqrt((kpts(1,i+1)-kpts(1,i+2))**2+(kpts(2,i+1)-kpts(2,i+2))**2+(kpts(3,i+1)-kpts(3,i+2))**2)
    x=x+increment
   end if
   old_increment=increment
   if (i==nkpt-1) xrange=x
  end do
  close(5)
  close(8)
  deallocate(kpts)

 end subroutine process_ebands

 subroutine write_plt(band_file,plt_file,nband,xrange,kptoption)

  implicit none

  ! Arguments
  character(len=20),intent(in)::band_file,plt_file
  integer,intent(in)::nband
  real,intent(in)::xrange
  integer,intent(in)::kptoption

  integer::i,ios

  open(9,file=plt_file,status="unknown",iostat=ios)
  if (ios/=0) then
   write(*,*)"Error opening file :",plt_file
   stop
  end if
  write(9,*)"#Gnuplot script generated by ebands utility."
  write(9,*)"y1=-15"
  write(9,*)"y2=15"
  write(9,*)"set style data lines"
  write(9,*)"set nokey"
  write(9,*)"set noxtics"
  do i=1,abs(kptoption)-1
   write(9,'(a,i2,a)')"set arrow",i," from 0,y1 to 0,y2 lt 9 nohead"
  end do
  write(9,'(1x,a,f7.4,a)')"set xrange [0:",xrange,"]"
  write(9,*)"set yrange[y1:y2]"
  write(9,'(1x,a,a,a)',advance='no')'plot ''',trim(band_file),''' u 2:3'
  do i=2,nband
   write(9,'(a,i3)',advance='no')','''' u 2:',i+2
  end do
  close(9)

 end subroutine write_plt

end program ebands
