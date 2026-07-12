file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/schemes.txt" SCHEMES)
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/fsize.txt"   FSIZES)
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/blsize.txt"  BLSIZES)

foreach(SCHEME ${SCHEMES})
  set(RUN_TESTS "")

  foreach(FSIZE ${FSIZES})
    foreach(BLSIZE ${BLSIZES})
      set(TEST_NAME "test_error_heat_flux_hybrid_${SCHEME}_${FSIZE}_${BLSIZE}")

      add_test(NAME ${TEST_NAME}
        COMMAND bash run_case.sh ${SCHEME} ${FSIZE} ${BLSIZE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

      list(APPEND RUN_TESTS ${TEST_NAME})
    endforeach()
  endforeach()

  set(VISU_TEST "test_error_heat_flux_hybrid_${SCHEME}_visu")

  add_test(NAME ${VISU_TEST}
    COMMAND bash visu.sh ${SCHEME}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  set_tests_properties(${VISU_TEST}
    PROPERTIES DEPENDS "${RUN_TESTS}")

endforeach()
