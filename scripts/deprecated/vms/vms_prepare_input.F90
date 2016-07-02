program vms_prepare_input
!
!******************************************************************************
!                                                                             *
!            Author : J.Jansen               Date : 14 August 2003            *
!                                                                             *
!-----------------------------------------------------------------------------*
!                                                                             *
! Manipulate input files for OpenVMS                                          *
!                                                                             *
!******************************************************************************
!
#if defined VMS
!
  implicit none
!
  character ( len = 80 ) :: file_name , string , out
  integer ( kind = 4 ) :: ipos , ipos1 , iout
!
  call getarg( 1 , file_name )
  ipos = lastns( file_name )
!
  open( unit=1 , file=file_name(1:ipos) , status='old' )
  open( unit=2 , file=file_name(1:ipos) , status='new' , carriagecontrol='list' )
  50 read( 1 , '(a)' , end=51 ) string
  ipos = lastns( string )
  out = '['
  iout = 1
  ipos1 = index( string , '/' )
  do while( ipos1 > 0 )
    if ( string( 1:ipos1-1 ) == '..' ) then
      iout = iout + 1
      out( iout:iout ) = '-'
    else
      out( iout+1: ) = '.' // string
      iout = iout + ipos1
    endif
    string = string( ipos1+1: )
    ipos = ipos - ipos1
    ipos1 = index( string , '/' )
  end do
  if ( iout == 1 ) then
    out = string( 1:ipos )
    iout = ipos
  else
    out( iout+1: ) = ']' // string( 1:ipos )
    iout = iout + ipos + 1
  end if
  write( 2 , '(a)' ) out( 1:iout )
  goto 50
!
  51 stop 'input updated'

contains

!
      subroutine getarg( numpar , param )
!
!*******************************************************************************
!                                                                              *
!            Author : J.Jansen               Date : 3 May 1993                 *
!   Revision : 1.0          Author : J.Jansen      Date : 25 September 2000    !
!   Revision : 2.0          Author : J.Jansen      Date : 6 June 2002          !
!   Revision : 3.0          Author : J.Jansen      Date : 13 May 2003          !
!                                                                              *
!------------------------------------------------------------------------------*
!                                                                              *
! Subroutine to get command parameters                                         *
!                                                                              *
! VMS(VAX,AXP) version                                                         *
!                                                                              *
!------------------------------------------------------------------------------*
!                                                                              *
!    Parameter         Type                   Purpose                          *
!    numpar           integer      Number of parameters                        *
!    param            char*(*)     parameter value                             *
!    para_name        char*2       Parameter name                              *
!                                                                              *
!*******************************************************************************
!
!DEC$ ATTRIBUTES ALIAS:'LIB$GET_FOREIGN' :: lib$get_foreign
!
        implicit none
!
        character*(*) param
        character ( len = 1024 ) string
        integer ( kind = 4 ) :: numpar , len_string , ipos , iend , i
!
      call lib$get_foreign( string , , len_string )
      if ( numpar == 0 ) then
        do while ( string( 1:1 ) == ' ' )
          string = string( 2: )
        end do
        ipos = index( string , ' ' )
        param = string( 1:ipos )
      else
        ipos = 1
        do i = 1 , numpar - 1
          do while ( string( ipos:ipos ) == ' ' )
            ipos = ipos + 1
          end do
          do while ( string( ipos:ipos ) /= ' ' )
            ipos = ipos + 1
          end do
        end do
        do while ( string( ipos:ipos ) == ' ' )
          ipos = ipos + 1
        end do
        iend = ipos + index( string( ipos: ) , ' ' ) - 1
        param = string( ipos:iend )
      end if
!
      return
      end subroutine getarg
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

#else
  stop 'VMS support has not been activated'
#endif

end program vms_prepare_input
