#set(CMAKE_POSITION_INDEPENDENT_CODE ON)

##################################################################
# Fortran
##################################################################

# full list of compilers supported by cmake :
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html

if(CMAKE_Fortran_COMPILER_ID MATCHES "Absoft")

  # TODO

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "ARMClang")

  # debug flags
  if (ABI_DEBUG_FLAVOR MATCHES "enhanced")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fbacktrace;-finit-real=nan;-Wimplicit>")
  elseif (ABI_DEBUG_FLAVOR MATCHES "paranoid")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fbacktrace;-finit-real=nan;-Wimplicit;-Wall;-Wextra>")
  elseif (ABI_DEBUG_FLAVOR MATCHES "naughty")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fbacktrace;-finit-real=nan;-Wimplicit;-Wall;-Wextra>")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Clang")

  # debug flags
  if (ABI_DEBUG_FLAVOR MATCHES "enhanced")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fbacktrace;-finit-real=nan;-Wimplicit>")
  elseif (ABI_DEBUG_FLAVOR MATCHES "paranoid")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fbacktrace;-finit-real=nan;-Wimplicit;-Wall;-Wextra>")
  elseif (ABI_DEBUG_FLAVOR MATCHES "naughty")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fbacktrace;-finit-real=nan;-Wimplicit;-Wall;-Wextra>")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Flang")

  # TODO

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

  if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    add_compile_options(
      "$<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none;-fallow-argument-mismatch>"
      )
  else()
    add_compile_options(
      "$<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>"
      )
  endif()

  # debug flags
  if (ABI_DEBUG_FLAVOR MATCHES "enhanced")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fbacktrace;-Werror=array-bounds;-finit-real=nan;-Wimplicit-interface;-Wno-maybe-uninitialized;-Wtabs>"
      )
  elseif(ABI_DEBUG_FLAVOR MATCHES "paranoid")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fbacktrace;-Werror=array-bounds;-finit-real=nan;-Wimplicit-interface;-Wno-maybe-uninitialized;-Wtabs;-ffpe-trap=invalid,zero,overflow;-Wall;-Wextra>")
  elseif(ABI_DEBUG_FLAVOR MATCHES "naughty")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fbacktrace;-Werror=array-bounds;-finit-real=nan;-Wimplicit-interface;-Wno-maybe-uninitialized;-Wtabs;-ffpe-trap=invalid,zero,overflow;-Wall;-Wextra;-fcheck=all;-pedantic>")
  endif()

  # release / optim flags
  add_compile_options(
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-fno-backtrace;-Wno-maybe-uninitialized>"
    )

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")

  add_compile_options(
    "$<$<COMPILE_LANGUAGE:Fortran>:-traceback;-heap-arrays>"
    )

  # debug flags
  if (ABI_DEBUG_FLAVOR MATCHES "enhanced")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-check uninit;-debug all;-fp-model source;-ftrapuv;-traceback>")
  elseif (ABI_DEBUG_FLAVOR MATCHES "paranoid")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-check uninit;-debug all;-fp-model source;-ftrapuv;-traceback;-fp-stack-check;-implicitnone;-init=snan;-warn all>")
  elseif (ABI_DEBUG_FLAVOR MATCHES "naughty")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-check uninit;-debug all;-fp-model source;-ftrapuv;-traceback-fp-stack-check;-implicitnone;-init=snan;-warn all;-check all;-WB>")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")

  # https://www.nag.co.uk/nagware/np/r70_doc/manual/compiler_2_4.html#OPTIONS
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:Fortran>:-f2018;-C;-colour;-gline;-u>"
    )

  # debug flags
  if (ABI_DEBUG_FLAVOR MATCHES "enhanced")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-mtrace=verbose;-nan>")
  elseif (ABI_DEBUG_FLAVOR MATCHES "paranoid")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-mtrace=verbose;-nan;-info>")
  elseif (ABI_DEBUG_FLAVOR MATCHES "naughty")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-mtrace=verbose;-nan;-info;-C>")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC") # NVFORTRAN

  # TODO : improve
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:Fortran>:-Mextend;-Mpreprocess;-Mfree;-tp=px;-traceback;-Minfo=mp,accel,par,pfo>")

  # debug flags
  add_compile_options(
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-traceback>")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")

  # TODO

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "XL")

  # debug flags
  if (ABI_DEBUG_FLAVOR MATCHES "enhanced")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-qnooptimize;-qextcheck;-qflag=i:i;-qfloat=nans;-qinitauto=7FBFFFFF>")
  elseif (ABI_DEBUG_FLAVOR MATCHES "paranoid")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-qnooptimize;-qextcheck;-qflag=i:i;-qfloat=nans;-qinitauto=7FBFFFFF;-qflttrap=overflow:underflow:zerodivide:invalid:enable;-qsigtrap>")
  elseif (ABI_DEBUG_FLAVOR MATCHES "naughty")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-qnooptimize;-qextcheck;-qflag=i:i;-qfloat=nans;-qinitauto=7FBFFFFF;-qflttrap=overflow:underflow:zerodivide:invalid:enable;-qsigtrap;-C;-qcheck>")
  endif()

endif()

##################################################################
# C
##################################################################

# full list of compilers supported by cmake :
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html

if(CMAKE_C_COMPILER_ID MATCHES "Clang")

  # release / optim flags
  if (ABI_OPTIM_FLAVOR MATCHES "safe")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Release>>:-O2>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "standard")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Release>>:-O2;-mtune=native;-mcpu=native>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "aggressive")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Release>>:-O3;-mtune=native;-mcpu=native;-ffp-contract=fast>"
      )
  endif()

elseif(CMAKE_C_COMPILER_ID MATCHES "ARMClang")

  # TODO

elseif(CMAKE_C_COMPILER_ID MATCHES "GNU")

  # debug flags
  if (ABI_DEBUG_FLAVOR MATCHES "enhanced")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Debug>>:-g3;-ggdb>"
      )
  elseif (ABI_DEBUG_FLAVOR MATCHES "paranoid")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Debug>>:-g3;-ggdb;-Wall;-Wextra>"
      )
  elseif (ABI_DEBUG_FLAVOR MATCHES "naughty")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Debug>>:-g3;-ggdb;-Wall;-Wextra>"
      )
  endif()

  # release / optim flags
  if (ABI_OPTIM_FLAVOR MATCHES "safe")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Release>>:-fno-backtrace;-Wno-maybe-uninitialized;-O2>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "standard")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Release>>:-fno-backtrace;-Wno-maybe-uninitialized;-O2;-mtune=native;-march=native>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "aggressive")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Release>>:-fno-backtrace;-Wno-maybe-uninitialized;-O3;-mtune=native;-march=native>"
      )
  endif()

elseif(CMAKE_C_COMPILER_ID MATCHES "^Intel")

  # release / optim flags
  if (ABI_OPTIM_FLAVOR MATCHES "safe")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Release>>:-O2>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "standard")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Release>>:-O2>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "aggressive")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Release>>:-O3>"
      )
  endif()

elseif(CMAKE_C_COMPILER_ID MATCHES "NVHPC") # NVC

  # TODO

elseif(CMAKE_C_COMPILER_ID MATCHES "PGI")

  # TODO

elseif(CMAKE_C_COMPILER_ID MATCHES "XL")

  # release / optim flags
  if (ABI_OPTIM_FLAVOR MATCHES "safe")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Release>>:-O2;-qarch=auto;-qtune=auto;-qstrict;-qspill=2000>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "standard")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Release>>:-O3;-qarch=auto;-qtune=auto;-qstrict;-qspill=2000>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "aggressive")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:Release>>:-O4;-qarch=auto;-qtune=auto;-qstrict;-qspill=2000>"
      )
  endif()

endif()

##################################################################
# CXX
##################################################################

# full list of compilers supported by cmake :
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")

  # release / optim flags
  if (ABI_OPTIM_FLAVOR MATCHES "safe")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:-O2>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "standard")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:-O2;-mtune=native;-mcpu=native>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "aggressive")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:-O3;-mtune=native;-mcpu=native;-ffp-contract=fast>"
      )
  endif()

elseif(CMAKE_CXX_COMPILER_ID MATCHES "ARMClang")

  # TODO

elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")

  # debug flags
  if (ABI_DEBUG_FLAVOR MATCHES "enhanced")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>>:-g3;-ggdb>"
      )
  elseif (ABI_DEBUG_FLAVOR MATCHES "paranoid")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>>:-g3;-ggdb;-Wall;-Wextra>"
      )
  elseif (ABI_DEBUG_FLAVOR MATCHES "naughty")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>>:-g3;-ggdb;-Wall;-Wextra>"
      )
  endif()

  # release / optim flags
  if (ABI_OPTIM_FLAVOR MATCHES "safe")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:-fno-backtrace;-Wno-maybe-uninitialized;-O2>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "standard")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:-fno-backtrace;-Wno-maybe-uninitialized;-O2;-mtune=native;-march=native>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "aggressive")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:-fno-backtrace;-Wno-maybe-uninitialized;-O3;-mtune=native;-march=native>"
      )
  endif()

elseif(CMAKE_CXX_COMPILER_ID MATCHES "^Intel")

  # release / optim flags
  if (ABI_OPTIM_FLAVOR MATCHES "safe")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:-O2>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "standard")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:-O2>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "aggressive")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:-O3>"
      )
  endif()

elseif(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC") # NVC++

  # TODO

elseif(CMAKE_CXX_COMPILER_ID MATCHES "PGI")

  # TODO

elseif(CMAKE_CXX_COMPILER_ID MATCHES "XL")

  # release / optim flags
  if (ABI_OPTIM_FLAVOR MATCHES "safe")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:-O2;-qarch=auto;-qtune=auto;-qstrict;-qspill=2000;-qessl>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "standard")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:-O3;-qarch=auto;-qtune=auto;-qstrict;-qspill=2000;-qessl>"
      )
  elseif (ABI_OPTIM_FLAVOR MATCHES "aggressive")
    add_compile_options(
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>:-O4;-qarch=auto;-qtune=auto;-qstrict;-qspill=2000;-qessl>"
      )
  endif()

endif()
