program crea_descrip_mms
!
!****************************************************************************!
!                                                                            !
! Creates OpenVMS descrip.mms from templates in ABINIT                       !
!                                                                            !
! Author : J.Jansen                                                          !
! Date : 13 August 2003                                                      !
! Revision 1.0 , 3 September 2004                                            !
! Revision 2.0 , 11 January 2006                                             !
!                                                                            !
!****************************************************************************!
!
#if defined VMS
!DEC$ ATTRIBUTES ALIAS:'LIB$FIND_FILE' :: lib$find_file
!
  implicit none
!
  include '($rmsdef)/list'
!
  character ( len = 256 ) :: card , card1
  integer ( kind = 4 ) :: ipos , ipos0 , icount = 0 , numobj = 1 , i
  integer ( kind = 4 ) :: lib$find_file , istat , ic
!
  open( unit=2 , file='descrip.mms' , status='new' , carriagecontrol='list' )
  open( unit=3 , status='scratch' )
!
! Write general rules
  write( 2 , '( a )' ) '.f.obj :'
  write( 2 , '( a )' ) '	cc/define=(VMS,HAVE_NETCDF,HAVE_FOX)/comment=as_is/prep=$(MMS$TARGET_NAME).f_\'
  write( 2 , '( a )' ) '	/include=([--.src.10_defs]) $(MMS$TARGET_NAME).F'
  write( 2 , '( a )' ) '	f90/source=fixed/convert=big_endian/name=lowercase\'
  write( 2 , '( a )' ) '	/opti=(tune=host,nopipe)/arch=host/debug/include=([--.src.10_defs],[--],sys$common:[syslib])\'
  write( 2 , '( a )' ) '	/check=(noarg,bo,form,fp,noout,over,powe,under) $(MMS$TARGET_NAME).f_'
  write( 2 , '( a )' ) '	delete $(MMS$TARGET_NAME).f_;*'
  write( 2 , '( )' )
  write( 2 , '( a )' ) '.f90.obj :'
  write( 2 , '( a )' ) '	cc/define=(VMS,HAVE_NETCDF,HAVE_FOX)/comment=as_is/prep=$(MMS$TARGET_NAME).f_\'
  write( 2 , '( a )' ) '	/include=([--.src.10_defs]) $(MMS$TARGET_NAME).F90'
  write( 2 , '( a )' ) '	f90/source=free/convert=big_endian/name=lowercase\'
  write( 2 , '( a )' ) '	/opti=(tune=host,nopipe)/arch=host/debug/include=([--.src.10_defs],[--],sys$common:[syslib])\'
  write( 2 , '( a )' ) '	/check=(noarg,bo,form,fp,noout,over,powe,under) $(MMS$TARGET_NAME).f_'
  write( 2 , '( a )' ) '	delete $(MMS$TARGET_NAME).f_;*'
  write( 2 , '( )' )
!
  ic = 0
  card = 'NL:'
  write( 2 , '( a )' ) 'OBJS =\'
  istat = lib$find_file( '*.F' , card1 , ic )
  if ( index( card1 , ']*' ) == 0 ) then
    card = card1
    istat = lib$find_file( '*.F' , card1 , ic )
    do while ( istat .ne. rms$_nmf )
      ipos0 = index( card , ']' ) + 1
      ipos = index( card , '.F' )
      write( 3 , '( a )' ) card( ipos0:ipos ) // 'obj : ' // card( ipos0:ipos+1 )
      if ( icount == 44 ) then
        write( 2 , '(a)' ) char( 9 ) // card( ipos0:ipos ) // 'obj'
        write( 2 , '( a , i1 , a )' ) 'OBJS' , numobj , ' =\'
        numobj = numobj + 1
        icount = -1
      else
        write( 2 , '(a)' ) char( 9 ) // card( ipos0:ipos ) // 'obj,\'
      end if
      icount = icount + 1
      card = card1
      istat = lib$find_file( '*.F' , card1 , ic )
    end do
  end if
  ic = 0
  istat = lib$find_file( '*.F90' , card1 , ic )
  if ( index( card1 , ']*' ) == 0 ) then
    if ( card( 1:3 ) .ne. 'NL:' ) then
      ipos0 = index( card , ']' ) + 1
      ipos = index( card , '.F' )
      write( 3 , '( a )' ) card( ipos0:ipos ) // 'obj : ' // card( ipos0:ipos+1 )
      if ( icount == 44 ) then
        write( 2 , '(a)' ) char( 9 ) // card( ipos0:ipos ) // 'obj'
        write( 2 , '( a , i1 , a )' ) 'OBJS' , numobj , ' =\'
        numobj = numobj + 1
        icount = -1
      else
        write( 2 , '(a)' ) char( 9 ) // card( ipos0:ipos ) // 'obj,\'
      end if
      icount = icount + 1
    end if
    card = card1
    istat = lib$find_file( '*.F90' , card1 , ic )
    do while ( istat .ne. rms$_nmf )
      ipos0 = index( card , ']' ) + 1
      ipos = index( card , '.F90' )
      write( 3 , '( a )' ) card( ipos0:ipos ) // 'obj : ' // card( ipos0:ipos+3 )
      if ( icount == 44 ) then
        write( 2 , '(a)' ) char( 9 ) // card( ipos0:ipos ) // 'obj'
        write( 2 , '( a , i1 , a )' ) 'OBJS' , numobj , ' =\'
        numobj = numobj + 1
        icount = -1
      else
        write( 2 , '(a)' ) char( 9 ) // card( ipos0:ipos ) // 'obj,\'
      end if
      icount = icount + 1
      card = card1
      istat = lib$find_file( '*.F90' , card1 , ic )
    end do
  end if
  ipos0 = index( card , ']' ) + 1
  ipos = index( card , '.F' )
  write( 3 , '( a )' ) card( ipos0:ipos ) // 'obj : ' // card( ipos0:ipos+3 )
  write( 2 , '(a)' ) char( 9 ) // card( ipos0:ipos ) // 'obj'
  ipos0 = ipos0 - 2
  ipos = ipos0
  do while ( card( ipos:ipos ) .ne. '.' )
    ipos = ipos - 1
  end do
  ipos = ipos + 1
  write( 2 , '( )' )
  select case ( numobj )
    case( 1 )
      write( 2 , '( a )' ) 'all : $(OBJS)'
    case( 2 )
      write( 2 , '( a )' ) 'all : $(OBJS) $(OBJS1)'
    case( 3 )
      write( 2 , '( a )' ) 'all : $(OBJS) $(OBJS1) $(OBJS2)'
    case( 4 )
      write( 2 , '( a )' ) 'all : $(OBJS) $(OBJS1) $(OBJS2) $(OBJS3)'
  end select
  write( 2 , '( a , a , a )' ) '	library/create lib' , card( &
&   ipos:ipos0 ) , '.olb $(OBJS)'
  do i = 1 , numobj - 1
    write( 2 , '( a , a , a , i1 , a )' ) '	library lib' , card( ipos:ipos0 ) , &
&     '.olb $(OBJS' , i , ')'
  end do
  write( 2 , '( )' )
  rewind 3
  53 read( 3 , '( a )' , end=51 ) card
  ipos = lastns( card )
  write( 2 , '(a)' ) card( 1:ipos )
  goto 53
!
  51 stop 'descrip.mms created'
  
contains

!
function lastns( card )
!
  integer ( kind = 4 ) :: lastns , result
  character( len = 80 ) :: card
!
  result = 80
  do while( card( result:result ) == ' ' )
    result = result - 1
    if ( result == 0 ) exit
  end do
!
  lastns = result
  return
end function lastns
!
#else
  stop 'VMS support has not been activated'
#endif
end program crea_descrip_mms
