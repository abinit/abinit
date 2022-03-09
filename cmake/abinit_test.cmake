#
# Some test target
#

add_custom_target(check-local)
add_custom_command(
  TARGET check-local
  POST_BUILD
  COMMAND cat ${CMAKE_SOURCE_DIR}/tests/built-in/README)

add_custom_target(test_fast)
add_custom_command(
  TARGET test_fast
  POST_BUILD
  COMMAND cp ${CMAKE_BINARY_DIR}/src/98_main/abinit ${CMAKE_SOURCE_DIR}/tests/built-in/Input/abinit && cd ${CMAKE_SOURCE_DIR}/tests/built-in/Input && export ABI_PSPDIR="${CMAKE_SOURCE_DIR}/tests/Psps_for_tests" && ./abinit testin_fast.abi >& testin_fast.stdout && cat testin_fastt_STATUS && rm -f abinit *DDB *EIG *out* *nc *WFK *abo* *o_* *t_STATUS*
  )

add_custom_target(test_v1)
add_custom_command(
  TARGET test_v1
  POST_BUILD
  COMMAND cp ${CMAKE_BINARY_DIR}/src/98_main/abinit ${CMAKE_SOURCE_DIR}/tests/built-in/Input/abinit && cd ${CMAKE_SOURCE_DIR}/tests/built-in/Input && export ABI_PSPDIR="${CMAKE_SOURCE_DIR}/tests/Psps_for_tests" && ./abinit testin_v1.abi >& testin_v1.stdout && cat testin_v1t_STATUS && rm -f abinit *DDB *EIG *out* *nc *WFK *abo* *o_* *t_STATUS*
  )

add_custom_target(test_v5)
add_custom_command(
  TARGET test_v5
  POST_BUILD
  COMMAND cp ${CMAKE_BINARY_DIR}/src/98_main/abinit ${CMAKE_SOURCE_DIR}/tests/built-in/Input/abinit && cd ${CMAKE_SOURCE_DIR}/tests/built-in/Input && export ABI_PSPDIR="${CMAKE_SOURCE_DIR}/tests/Psps_for_tests" && ./abinit testin_v5.abi >& testin_v5.stdout && cat testin_v5t_STATUS && rm -f abinit *DDB *EIG *out* *nc *WFK *abo* *o_* *t_STATUS*
  )

if (TARGET abinit::libbigdft)
  add_custom_target(test_bigdft)
  add_custom_command(
    TARGET test_bigdft
    POST_BUILD
    COMMAND cp ${CMAKE_BINARY_DIR}/src/98_main/abinit ${CMAKE_SOURCE_DIR}/tests/built-in/Input/abinit && cd ${CMAKE_SOURCE_DIR}/tests/built-in/Input && export ABI_PSPDIR="${CMAKE_SOURCE_DIR}/tests/Psps_for_tests" && ./abinit testin_bigdft.abi >& testin_bigdft.stdout && cat testin_bigdftt_STATUS && rm -f abinit *DDB *EIG *out* *nc *WFK *abo* *o_* *t_STATUS*)
endif()

if (TARGET abinit::libxc)
  add_custom_target(test_libxc)
  add_custom_command(
    TARGET test_libxc
    POST_BUILD
    COMMAND cp ${CMAKE_BINARY_DIR}/src/98_main/abinit ${CMAKE_SOURCE_DIR}/tests/built-in/Input/abinit && cd ${CMAKE_SOURCE_DIR}/tests/built-in/Input && export ABI_PSPDIR="${CMAKE_SOURCE_DIR}/tests/Psps_for_tests" && ./abinit testin_libxc.abi >& testin_libxc.stdout && cat testin_libxct_STATUS && rm -f abinit *DDB *EIG *out* *nc *WFK *abo* *o_* *t_STATUS*)
endif()

if (TARGET abinit::wannier)
  add_custom_target(test_wannier90)
  add_custom_command(
    TARGET test_wannier90
    POST_BUILD
    COMMAND cp ${CMAKE_BINARY_DIR}/src/98_main/abinit ${CMAKE_SOURCE_DIR}/tests/built-in/Input/abinit && cd ${CMAKE_SOURCE_DIR}/tests/built-in/Input && export ABI_PSPDIR="${CMAKE_SOURCE_DIR}/tests/Psps_for_tests" && ./abinit testin_wannier90.abi >& testin_wannier90.stdout && cat testin_wannier90t_STATUS  && rm -f abinit *DDB *EIG *out* *nc *WFK *DEN *chk *eig *mmn *amn *abo* *o_* *t_STATUS*)
endif()

add_custom_target(tests_in DEPENDS
  test_fast test_v1 test_v5 test_libxc test_wannier90)
