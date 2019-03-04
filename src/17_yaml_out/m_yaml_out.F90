#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
#define SET_DEFAULT(v, optv, defv) v = defv; if(present(optv)) v = optv
#define ERROR(msg) write(0,*) msg; stop
#define ERROR_NO_OUT ERROR("No output medium have been provided.")
#define ASSERT(cond, msg) if(.not. cond) then; ERROR(msg); end if
module m_yaml_out


  use m_pair_list
  use m_stream_string

  implicit none

  private

  integer,parameter :: dp=kind(1.0D0)
  character(len=1),parameter :: eol=char(10)
  character(len=9),parameter :: default_rfmt='(ES25.17)'
  character(len=4),parameter :: default_ifmt='(I8)'
  character(len=13),parameter :: default_kfmt="(A)"
  character(len=13),parameter :: default_sfmt="(A)"
  integer,parameter :: default_keysize=30
  integer,parameter :: default_stringsize=500

  integer :: doclock


  public :: yaml_open_doc, yaml_close_doc, yaml_single_dict, yaml_iterstart
  public :: yaml_add_realfield, yaml_add_intfield, yaml_add_stringfield
  public :: yaml_add_real1d, yaml_add_real2d
  public :: yaml_add_dict,  yaml_add_dictlist
  public :: yaml_add_int1d, yaml_add_int2d

  contains

! private
  subroutine string_clear(string)
    character(len=*),intent(inout) :: string
    string = repeat(' ', len(string))
  end subroutine string_clear

  pure function yaml_quote_string(string) result(quoted)
    character(len=*),intent(in) :: string
    character(len=len(string)+2) :: quoted
    logical :: multiline, spec_char, quote
    spec_char = index(string, ':', back=.true.) /= 0
    spec_char = spec_char .or. index(string, '{') /= 0
    spec_char = spec_char .or. index(string, '}') /= 0
    spec_char = spec_char .or. index(string, '[') /= 0
    spec_char = spec_char .or. index(string, ']') /= 0
    spec_char = spec_char .or. index(string, ',') /= 0
    spec_char = spec_char .or. index(string, '&') /= 0
    spec_char = spec_char .or. index(string, '*') /= 0
    spec_char = spec_char .or. index(string, '#') /= 0
    spec_char = spec_char .or. index(string, '?') /= 0
    spec_char = spec_char .or. index(string, '|') /= 0
    spec_char = spec_char .or. index(string, '-') /= 0
    spec_char = spec_char .or. index(string, '<') /= 0
    spec_char = spec_char .or. index(string, '>') /= 0
    spec_char = spec_char .or. index(string, '=') /= 0
    spec_char = spec_char .or. index(string, '!') /= 0
    spec_char = spec_char .or. index(string, '%') /= 0
    spec_char = spec_char .or. index(string, '@') /= 0
    spec_char = spec_char .or. index(string, '`') /= 0

    quote = index(string, "'") /= 0
    multiline = index(string, eol, back=.true.) /= 0


    if (quote) then
      quoted='"'//string//'"'
    else if(multiline .or. spec_char) then
      quoted="'"//string//"'"
    else
      quoted=string
    endif
  end function yaml_quote_string

  subroutine yaml_start_field(stream, label, tag, width)
    type(stream_string),intent(inout) :: stream
    character(len=*),intent(in) :: label
    integer,optional :: width
    character(len=*),intent(in),optional :: tag
    character(len=len_trim(label)+2) :: quoted

    quoted = yaml_quote_string(label)
    if(present(width)) then
      if(width > len_trim(label)) then
        call stream%write(trim(quoted)//repeat(' ', width-len_trim(quoted))//':')
      else
        call stream%write(trim(quoted)//':')
      end if
    else
      call stream%write(trim(quoted)//':')
    end if
    if (present(tag)) then
      call stream%write(' !'//trim(tag))
    end if
  end subroutine yaml_start_field

  subroutine yaml_print_real1d(stream, length, arr, rfmt, vmax)
    type(stream_string),intent(inout) :: stream
    integer,intent(in) :: vmax
    integer,intent(in) :: length
    real(kind=dp),intent(in) :: arr(length)
    character(len=*),intent(in) :: rfmt
    character(len=50) :: tmp_r

    integer :: i

    if (length > vmax) then
      call stream%write(' ['//eol//'    ')
    else
      call stream%write(' [')
    end if

    do i=1,length
      call string_clear(tmp_r)
      write(tmp_r, rfmt) arr(i)
      call stream%write(trim(tmp_r))
      if (i > 0 .and. mod(i, vmax) == 0 .and. i /= length) then
        call stream%write(', '//eol//'    ')
      else
        call stream%write(', ')
      end if
    end do
    if (length > vmax) then
      call stream%write(eol)
    end if
    call stream%write(']')
  end subroutine yaml_print_real1d

  subroutine yaml_print_int1d(stream, length, arr, ifmt, vmax)
    type(stream_string),intent(inout) :: stream
    integer,intent(in) :: vmax
    integer,intent(in) :: length
    integer,intent(in) :: arr(length)
    character(len=*),intent(in) :: ifmt
    character(len=50) :: tmp_i

    integer :: i

    if (length > vmax) then
      call stream%write(' ['//eol//'    ')
    else
      call stream%write(' [')
    end if

    do i=1,length
      call string_clear(tmp_i)
      write(tmp_i, ifmt) arr(i)
      call stream%write(trim(tmp_i))
      if (i > 0 .and. mod(i, vmax) == 0 .and. i /= length) then
        call stream%write(', '//eol//'    ')
      else
        call stream%write(', ')
      end if
    end do
    if (length > vmax) then
      call stream%write(eol)
    end if
    call stream%write(']')
  end subroutine yaml_print_int1d

  subroutine yaml_print_dict(stream, pl, key_size, s_size, kfmt, ifmt, rfmt, sfmt, vmax)
    type(stream_string),intent(inout) :: stream
    integer,intent(in) :: vmax
    type(pair_list),intent(inout) :: pl
    character(len=*),intent(in) :: ifmt, rfmt, kfmt, sfmt
    integer,intent(in) :: key_size, s_size
    character(len=key_size) :: key
    character(len=key_size+5) :: tmp_key
    character(len=100) :: tmp_r
    character(len=100) :: tmp_i
    character(len=s_size) :: tmp_s

    integer :: i, vi, type_code
    real(kind=dp) :: vr
    character(len=s_size) :: vs

    if (pl%length() > vmax) then
      call stream%write(' {'//eol//'    ')
    else
      call stream%write(' {')
    end if

    call pl%restart()
    do i=1,pl%length()
      call pl%iter(key, type_code, vi, vr, vs)

      call string_clear(tmp_key)
      write(tmp_key, kfmt) '"'//trim(key)//'"'
      call stream%write(trim(tmp_key)//': ')
      if(type_code == TC_INT) then
        call string_clear(tmp_i)
        write(tmp_i, ifmt) vi
        call stream%write(trim(tmp_i))
      else if(type_code == TC_REAL) then
        call string_clear(tmp_r)
        write(tmp_r, rfmt) vr
        call stream%write(trim(tmp_r))
      else if(type_code == TC_STRING) then
        call string_clear(tmp_s)
        write(tmp_s, sfmt) vs
        call yaml_print_string(stream, trim(tmp_s), 100)
      end if
      if (i > 0 .and. mod(i, vmax) == 0 .and. i /= pl%length()) then
        call stream%write(', '//eol//'    ')
      else
        call stream%write(', ')
      end if
    end do

    if (pl%length() > vmax) then
      call stream%write(eol)
    end if
    call stream%write('}')
  end subroutine yaml_print_dict

  subroutine yaml_print_string(stream, string, vmax)
    type(stream_string),intent(inout) :: stream
    integer,intent(in) :: vmax
    character(len=*),intent(in) :: string
    character(len=len_trim(string)+2) :: quoted

    logical :: auto_wrap

    if(len(string) > vmax) then
      auto_wrap = .true.
    else
      auto_wrap = .false.
    end if


    quoted = yaml_quote_string(string)
    call stream%write(trim(quoted))
  end subroutine yaml_print_string

! public
  subroutine yaml_iterstart(label, val, file_d, string, stream, newline, width)
    integer,intent(in) :: val
    character(len=*),intent(in) :: label
    integer,intent(in), optional :: file_d
    type(stream_string),intent(inout),optional :: stream
    character(len=*),intent(out),optional :: string
    logical,intent(in),optional :: newline
    integer,intent(in),optional :: width

    integer :: w
    character(len=6) :: tmp_i
    logical :: nl

    SET_DEFAULT(nl, newline, .true.)
    SET_DEFAULT(w, width, 0)
    write(tmp_i, '(I6)') val

    if(present(stream)) then
      call stream%write('--- !IterStart'//eol//label//':'//tmp_i//eol//'...')
      if(nl) call stream%write(eol)
    else if(present(string)) then
      if(nl) then
        write(string, '(A)') '--- !IterStart'//eol//label//':'//tmp_i//eol//'...'//eol
      else
        write(string, '(A)') '--- !IterStart'//eol//label//':'//tmp_i//eol//'...'
      end if
    else if(present(file_d)) then
      if(nl) then
        write(file_d, '(A)') '--- !IterStart'//eol//label//':'//tmp_i//eol//'...'
      else
        write(file_d, '(A)', advance='no') '--- !IterStart'//eol//label//':'//tmp_i//eol//'...'
      end if
    else
      ERROR_NO_OUT
    end if
  end subroutine yaml_iterstart

  subroutine yaml_open_doc(label, comment, tag, file_d, string, stream, newline, width)
    character(len=*),intent(in) :: label
    character(len=*),intent(in) :: comment
    character(len=*),intent(in),optional :: tag
    integer,intent(in), optional :: file_d
    type(stream_string),intent(inout),optional :: stream
    character(len=*),intent(out),optional :: string
    logical,intent(in),optional :: newline
    integer,intent(in),optional :: width

    integer :: w
    type(stream_string) :: interm
    logical :: nl

    SET_DEFAULT(nl, newline, .true.)
    SET_DEFAULT(w, width, 0)
  
    if (doclock == 1) then
     call interm%write('...')
    end if
    doclock = 1
    
    if(present(tag)) then
      call interm%write('---'//' !'//trim(tag)//eol//'label')
    else
      call interm%write('---'//eol//'label')
    end if
    if(present(width)) then
      if(width > 5) then
        call interm%write(repeat(' ', width - 5))
      end if
    end if
    call interm%write(': '//trim(label))

    if (comment /= '') then
      call interm%write(eol//'comment')
      if(present(width)) then
        if(width > 7) then
          call interm%write(repeat(' ', width - 7))
        end if
      end if
      call interm%write(': ')
      call yaml_print_string(interm, comment, 70)
    end if
    if(nl) call interm%write(eol)

    if(present(stream)) then
      call interm%transfer(stream)
    else if(present(string)) then
      call interm%to_string(string)
    else if(present(file_d)) then
      call interm%to_file(file_d)
    else
      ERROR_NO_OUT
    end if
  end subroutine yaml_open_doc

  subroutine yaml_add_realfield(label, val, file_d, string, stream, tag, &
&                               real_fmt, newline, width)
    real(kind=dp) :: val
    character(len=*),intent(in) :: label
    character(len=*),intent(in),optional :: tag, real_fmt
    character(len=30) :: rfmt
    integer,intent(in), optional :: file_d
    type(stream_string),intent(inout),optional :: stream
    character(len=*),intent(out),optional :: string
    logical,intent(in),optional :: newline
    integer,intent(in),optional :: width

    integer :: w
    character(len=50) :: tmp_r
    type(stream_string) :: interm
    logical :: nl

    SET_DEFAULT(nl, newline, .true.)
    SET_DEFAULT(w, width, 0)

    ASSERT(doclock == 1, "No document is opened yet.")

    rfmt = '                              '
    SET_DEFAULT(rfmt, real_fmt, default_rfmt)
    if(present(tag)) then
      call yaml_start_field(interm, label, width=w, tag=tag)
    else
      call yaml_start_field(interm, label, width=w)
    end if

    call interm%write(' ')
    write(tmp_r, trim(rfmt)) val
    call interm%write(trim(tmp_r))
    if(nl) call interm%write(eol)

    if(present(stream)) then
      call interm%transfer(stream)
    else if(present(string)) then
      call interm%to_string(string)
    else if(present(file_d)) then
      call interm%to_file(file_d)
    else
      ERROR_NO_OUT
    end if
  end subroutine yaml_add_realfield

  subroutine yaml_add_intfield(label, val, file_d, string, stream, tag, &
&                              int_fmt, newline, width)
    integer :: val
    character(len=*),intent(in) :: label
    character(len=*),intent(in),optional :: tag, int_fmt
    character(len=30) :: ifmt
    integer,intent(in), optional :: file_d
    type(stream_string),intent(inout),optional :: stream
    character(len=*),intent(out),optional :: string
    logical,intent(in),optional :: newline
    integer,intent(in),optional :: width

    integer :: w
    character(50) :: tmp_i
    type(stream_string) :: interm
    logical :: nl

    SET_DEFAULT(nl, newline, .true.)
    SET_DEFAULT(w, width, 0)

    ASSERT(doclock == 1, "No document is opened yet.")

    ifmt = '                              '
    SET_DEFAULT(ifmt, int_fmt, default_ifmt)
    if(present(tag)) then
      call yaml_start_field(interm, label, width=w, tag=tag)
    else
      call yaml_start_field(interm, label, width=w)
    end if

    call interm%write(' ')
    write(tmp_i, trim(ifmt)) val
    call interm%write(trim(tmp_i))
    if(nl) call interm%write(eol)

    if(present(stream)) then
      call interm%transfer(stream)
    else if(present(string)) then
      call interm%to_string(string)
    else if(present(file_d)) then
      call interm%to_file(file_d)
    else
      ERROR_NO_OUT
    end if
  end subroutine yaml_add_intfield

  subroutine yaml_add_stringfield(label, val, file_d, string, stream, tag, newline, width)
    character(len=*) :: val
    character(len=*),intent(in) :: label
    character(len=*),intent(in),optional :: tag
    integer,intent(in), optional :: file_d
    type(stream_string),intent(inout),optional :: stream
    character(len=*),intent(out),optional :: string
    logical,intent(in),optional :: newline
    integer,intent(in),optional :: width

    integer :: w
    type(stream_string) :: interm
    logical :: nl

    SET_DEFAULT(nl, newline, .true.)
    SET_DEFAULT(w, width, 0)

    ASSERT(doclock == 1, "No document is opened yet.")

    if(present(tag)) then
      call yaml_start_field(interm, label, width=w, tag=tag)
    else
      call yaml_start_field(interm, label, width=w)
    end if

    call interm%write(' ')
    call yaml_print_string(interm, val, 70)
    if(nl) call interm%write(eol)

    if(present(stream)) then
      call interm%transfer(stream)
    else if(present(string)) then
      call interm%to_string(string)
    else if(present(file_d)) then
      call interm%to_file(file_d)
    else
      ERROR_NO_OUT
    end if
  end subroutine yaml_add_stringfield

  subroutine yaml_add_real1d(label, length, arr, file_d, string, stream, tag, real_fmt, multiline_trig, newline, width)
    integer,intent(in) :: length
    integer,intent(in),optional :: multiline_trig
    real(kind=dp),intent(in) :: arr(length)
    character(len=*),intent(in) :: label
    character(len=*),intent(in),optional :: tag, real_fmt
    integer,intent(in), optional :: file_d
    type(stream_string),intent(inout),optional :: stream
    character(len=*),intent(out),optional :: string
    logical,intent(in),optional :: newline
    integer,intent(in),optional :: width

    type(stream_string) :: interm
    integer :: w
    character(len=30) :: rfmt
    integer :: vmax
    logical :: nl

    SET_DEFAULT(nl, newline, .true.)
    SET_DEFAULT(w, width, 0)

    ASSERT(doclock == 1, "No document is opened yet.")

    rfmt = '                              '
    SET_DEFAULT(rfmt, real_fmt, default_rfmt)
    SET_DEFAULT(vmax, multiline_trig, 5)

    if(present(tag)) then
      call yaml_start_field(interm, label, width=w, tag=tag)
    else
      call yaml_start_field(interm, label, width=w)
    end if

    call yaml_print_real1d(interm, length, arr, trim(rfmt), vmax)
    if(nl) call interm%write(eol)

    if(present(stream)) then
      call interm%transfer(stream)
    else if(present(string)) then
      call interm%to_string(string)
    else if(present(file_d)) then
      call interm%to_file(file_d)
    else
      ERROR_NO_OUT
    end if
  end subroutine yaml_add_real1d

  subroutine yaml_add_int1d(label, length, arr, file_d, string, stream, tag, int_fmt, multiline_trig, newline, width)
    integer,intent(in) :: length
    integer,intent(in),optional :: multiline_trig
    integer,intent(in) :: arr(length)
    character(len=*),intent(in) :: label
    character(len=*),intent(in),optional :: tag, int_fmt
    integer,intent(in), optional :: file_d
    type(stream_string),intent(inout),optional :: stream
    character(len=*),intent(out),optional :: string
    logical,intent(in),optional :: newline
    integer,intent(in),optional :: width

    type(stream_string) :: interm
    character(len=30) :: ifmt
    integer :: w
    integer :: vmax
    logical :: nl

    SET_DEFAULT(nl, newline, .true.)
    SET_DEFAULT(w, width, 0)

    ASSERT(doclock == 1, "No document is opened yet.")

    ifmt = '                              '
    SET_DEFAULT(ifmt, int_fmt, default_ifmt)
    SET_DEFAULT(vmax, multiline_trig, 5)

    if(present(tag)) then
      call yaml_start_field(interm, label, width=w, tag=tag)
    else
      call yaml_start_field(interm, label, width=w)
    end if

    call yaml_print_int1d(interm, length, arr, trim(ifmt), vmax)
    if(nl) call interm%write(eol)

    if(present(stream)) then
      call interm%transfer(stream)
    else if(present(string)) then
      call interm%to_string(string)
    else if(present(file_d)) then
      call interm%to_file(file_d)
    else
      ERROR_NO_OUT
    end if
  end subroutine yaml_add_int1d

  subroutine yaml_add_dict(label, pl, file_d, string, stream, tag, key_size, string_size, key_fmt, int_fmt, real_fmt, string_fmt, &
&                          multiline_trig, newline, width)
    type(pair_list),intent(inout) :: pl
    character(len=*),intent(in) :: label
    integer,intent(in),optional :: string_size, key_size, multiline_trig
    character(len=*),intent(in),optional :: tag, key_fmt, int_fmt, real_fmt, string_fmt
    integer,intent(in), optional :: file_d
    type(stream_string),intent(inout),optional :: stream
    character(len=*),intent(out),optional :: string
    logical,intent(in),optional :: newline
    integer,intent(in),optional :: width

    type(stream_string) :: interm
    integer :: w
    character(len=30) :: kfmt, ifmt, rfmt, sfmt
    integer :: vmax, ks, ss
    logical :: nl

    SET_DEFAULT(nl, newline, .true.)
    SET_DEFAULT(w, width, 0)
    SET_DEFAULT(ks, key_size, default_keysize)
    SET_DEFAULT(ss, string_size, default_stringsize)

    ASSERT(doclock == 1, "No document is opened yet.")

    kfmt = '                              '
    rfmt = '                              '
    ifmt = '                              '
    sfmt = '                              '
    SET_DEFAULT(kfmt, key_fmt, default_kfmt)
    SET_DEFAULT(rfmt, real_fmt, default_rfmt)
    SET_DEFAULT(ifmt, int_fmt, default_ifmt)
    SET_DEFAULT(sfmt, string_fmt, default_sfmt)
    SET_DEFAULT(vmax, multiline_trig, 5)

    if(present(tag)) then
      call yaml_start_field(interm, label, width=w, tag=tag)
    else
      call yaml_start_field(interm, label, width=w)
    end if

    call yaml_print_dict(interm, pl, ks, ss, trim(kfmt), trim(ifmt), trim(rfmt), trim(sfmt), vmax)
    if(nl) call interm%write(eol)

    if(present(stream)) then
      call interm%transfer(stream)
    else if(present(string)) then
      call interm%to_string(string)
    else if(present(file_d)) then
      call interm%to_file(file_d)
    else
      ERROR_NO_OUT
    end if
  end subroutine yaml_add_dict

  subroutine yaml_add_real2d(label, m, n, arr, file_d, string, stream, tag, real_fmt, multiline_trig, newline, width)
    integer,intent(in) :: m, n
    real(kind=dp),intent(in) :: arr(m, n)
    character(len=*),intent(in) :: label
    character(len=*),intent(in),optional :: tag, real_fmt
    integer,intent(in),optional :: multiline_trig
    integer,intent(in), optional :: file_d
    type(stream_string),intent(inout),optional :: stream
    character(len=*),intent(out),optional :: string
    logical,intent(in),optional :: newline
    integer,intent(in),optional :: width

    type(stream_string) :: interm
    integer :: w
    integer :: i, vmax
    real(kind=dp) :: line(n)
    character(len=30) :: rfmt
    logical :: nl

    SET_DEFAULT(nl, newline, .true.)
    SET_DEFAULT(w, width, 0)

    ASSERT(doclock == 1, "No document is opened yet.")
    
    rfmt = '                              '
    SET_DEFAULT(rfmt, real_fmt, default_rfmt)
    SET_DEFAULT(vmax, multiline_trig, 5)

    if(present(tag)) then
      call yaml_start_field(interm, label, width=w, tag=tag)
    else
      call yaml_start_field(interm, label, width=w)
    end if

    do i=1,m
      call interm%write(eol//'-')
      line = arr(i,:)
      call yaml_print_real1d(interm, n, line, rfmt, vmax)
    end do

    if(nl) call interm%write(eol)

    if(present(stream)) then
      call interm%transfer(stream)
    else if(present(string)) then
      call interm%to_string(string)
    else if(present(file_d)) then
      call interm%to_file(file_d)
    else
      ERROR_NO_OUT
    end if
  end subroutine yaml_add_real2d

  subroutine yaml_add_int2d(interm, label, m, n, arr, file_d, string, stream, tag, int_fmt, multiline_trig, newline, width)
    integer,intent(in) :: m, n
    integer,intent(in) :: arr(m, n)
    character(len=*),intent(in) :: label
    character(len=*),intent(in),optional :: tag, int_fmt
    integer,intent(in),optional :: multiline_trig
    integer,intent(in), optional :: file_d
    type(stream_string),intent(inout),optional :: stream
    character(len=*),intent(out),optional :: string
    logical,intent(in),optional :: newline
    integer,intent(in),optional :: width

    type(stream_string) :: interm
    integer :: w
    integer :: i, vmax
    integer :: line(n)
    character(len=30) :: ifmt
    logical :: nl

    SET_DEFAULT(nl, newline, .true.)
    SET_DEFAULT(w, width, 0)

    ASSERT(doclock == 1, "No document is opened yet.")
    
    ifmt = '                              '
    SET_DEFAULT(ifmt, int_fmt, default_ifmt)
    SET_DEFAULT(vmax, multiline_trig, 5)

    if(present(tag)) then
      call yaml_start_field(interm, label, width=w, tag=tag)
    else
      call yaml_start_field(interm, label, width=w)
    end if

    do i=1,m
      call interm%write(eol//'-')
      line = arr(i,:)
      call yaml_print_int1d(interm, n, line, ifmt, vmax)
    end do

    if(nl) call interm%write(eol)

    if(present(stream)) then
      call interm%transfer(stream)
    else if(present(string)) then
      call interm%to_string(string)
    else if(present(file_d)) then
      call interm%to_file(file_d)
    else
      ERROR_NO_OUT
    end if
  end subroutine yaml_add_int2d

  subroutine yaml_add_dictlist(label, n, plarr, file_d, string, stream, tag, key_size, string_size, key_fmt, int_fmt, &
&                              real_fmt, string_fmt, multiline_trig, newline, width)
    integer,intent(in) :: n
    type(pair_list),intent(inout) :: plarr(n)
    character(len=*),intent(in) :: label
    integer,intent(in),optional :: key_size, string_size
    integer,intent(in),optional :: multiline_trig
    character(len=*),intent(in),optional :: tag, key_fmt, int_fmt, real_fmt, string_fmt
    integer,intent(in), optional :: file_d
    type(stream_string),intent(inout),optional :: stream
    character(len=*),intent(out),optional :: string
    logical,intent(in),optional :: newline
    integer,intent(in),optional :: width

    type(stream_string) :: interm
    integer :: w
    character(len=30) :: kfmt, ifmt, rfmt, sfmt
    integer :: vmax, ks, i, ss
    logical :: nl

    ASSERT(doclock == 1, "No document is opened yet.")
    
    SET_DEFAULT(nl, newline, .true.)
    SET_DEFAULT(w, width, 0)

    kfmt = '                              '
    rfmt = '                              '
    ifmt = '                              '
    SET_DEFAULT(kfmt, key_fmt, default_kfmt)
    SET_DEFAULT(rfmt, real_fmt, default_rfmt)
    SET_DEFAULT(ifmt, int_fmt, default_ifmt)
    SET_DEFAULT(sfmt, string_fmt, default_sfmt)
    SET_DEFAULT(vmax, multiline_trig, 5)
    SET_DEFAULT(ks, key_size, default_keysize)
    SET_DEFAULT(ss, string_size, default_keysize)

    if(present(tag)) then
      call yaml_start_field(interm, label, width=w, tag=tag)
    else
      call yaml_start_field(interm, label, width=w)
    end if
    call interm%write(eol)

    do i=1,n
      call interm%write('- ')
      call yaml_print_dict(interm, plarr(i), ks, ss, trim(kfmt), trim(ifmt), trim(rfmt), trim(sfmt), vmax)
      if(nl .or. i/=n) then
        call interm%write(eol)
      end if
    end do


    if(present(stream)) then
      call interm%transfer(stream)
    else if(present(string)) then
      call interm%to_string(string)
    else if(present(file_d)) then
      call interm%to_file(file_d)
    else
      ERROR_NO_OUT
    end if
  end subroutine yaml_add_dictlist

  subroutine yaml_single_dict(label, comment, pl, key_size, string_size, file_d, string, stream, tag, &
&                             int_fmt, real_fmt, string_fmt, newline, width)
    type(pair_list),intent(inout) :: pl
    character(len=*),intent(in) :: label
    character(len=*),intent(in) :: comment
    integer,intent(in) :: key_size, string_size
    character(len=*),intent(in),optional :: tag, int_fmt, real_fmt, string_fmt
    integer,intent(in), optional :: file_d, width
    type(stream_string),intent(inout),optional :: stream
    character(len=*),intent(out),optional :: string
    logical,intent(in),optional :: newline

    type(stream_string) :: interm
    character(len=30) :: kfmt, ifmt, rfmt, sfmt
    character(len=string_size) :: vs, tmp_s
    character(len=key_size) :: key
    integer :: vi, k, type_code, w
    character(len=50) :: tmp_i, tmp_r
    real(kind=dp) :: vr
    logical :: nl
    
    SET_DEFAULT(nl, newline, .true.)

    kfmt = '                              '
    rfmt = '                              '
    ifmt = '                              '
    SET_DEFAULT(rfmt, real_fmt, default_rfmt)
    SET_DEFAULT(ifmt, int_fmt, default_ifmt)
    SET_DEFAULT(sfmt, string_fmt, default_sfmt)
    SET_DEFAULT(w, width, 0)

    if (doclock == 1) then
     call interm%write('...'//eol)
    end if
    doclock = 1
    
    call interm%write('---')
    if(present(tag)) then
      call interm%write(' !'//tag)
    end if
    call interm%write(eol)
    call yaml_start_field(interm, 'label', width=width)
    call interm%write(' '//label)

    if (comment /= '') then
      call interm%write(eol)
      call yaml_start_field(interm, 'comment', width=width)
      call yaml_print_string(interm, comment, 70)
    end if
    call interm%write(eol)

    call pl%restart()
    do k=1,pl%length()
      call string_clear(key)
      call string_clear(vs)
      call pl%iter(key, type_code, vi, vr, vs)
      write(*,*) key, type_code

      call yaml_start_field(interm, trim(key), width=width)
      call interm%write(' ')
      if(type_code == TC_INT) then
        call string_clear(tmp_i)
        write(tmp_i, ifmt) vi
        call interm%write(trim(tmp_i))
      else if(type_code == TC_REAL) then
        call string_clear(tmp_r)
        write(tmp_r, rfmt) vr
        call interm%write(trim(tmp_r))
      else if(type_code == TC_STRING) then
        call string_clear(tmp_s)
        write(tmp_s, sfmt) vs
        call yaml_print_string(interm, trim(tmp_s), 100)
      end if
      call interm%write(eol)
    end do

    call yaml_close_doc(stream=interm, newline=nl)

    if(present(stream)) then
      call interm%transfer(stream)
    else if(present(string)) then
      call interm%to_string(string)
    else if(present(file_d)) then
      call interm%to_file(file_d)
    else
      ERROR_NO_OUT
    end if
  end subroutine yaml_single_dict

  subroutine yaml_close_doc(file_d, string, stream, newline)
    integer,intent(in), optional :: file_d
    type(stream_string),intent(inout),optional :: stream
    character(len=*),intent(out),optional :: string
    logical,intent(in),optional :: newline
    logical :: nl

    SET_DEFAULT(nl, newline, .true.)

    if (doclock == 1) then
      if(present(stream)) then
        call stream%write('...')
        if(nl) call stream%write(eol)
      else if(present(string)) then
        if(nl) then
          write(string, '(A)') '...'//eol
        else
          write(string, '(A)') '...'
        end if
      else if(present(file_d)) then
        if(nl) then
          write(file_d, '(A)') '...'
        else
          write(file_d, '(A)', advance='no') '...'
        end if
      else
        ERROR_NO_OUT
      end if
    end if
    doclock = 0
  end subroutine yaml_close_doc
end module m_yaml_out
