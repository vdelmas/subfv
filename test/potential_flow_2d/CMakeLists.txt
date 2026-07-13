file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/../schemes.txt" SCHEMES)
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/methods.txt" METHODS)

set(SECOND ".false.")
set(METHOD "0")
foreach(SCHEME ${SCHEMES})
  set(TEST_NAME "potential_flow_2d_${SCHEME}")

  add_test(NAME ${TEST_NAME}
    COMMAND bash run.sh ${SCHEME} ${SECOND} ${METHOD}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_test(NAME ${TEST_NAME}_visu
    COMMAND bash visu.sh ${SCHEME} ${SECOND} ${METHOD}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  set_tests_properties(${TEST_NAME}_visu PROPERTIES DEPENDS ${TEST_NAME})

endforeach()

set(SECOND ".true.")
foreach(METHOD ${METHODS})
  foreach(SCHEME ${SCHEMES})
    set(TEST_NAME "potential_flow_2d_${SCHEME}_o2_${METHOD}")

    add_test(NAME ${TEST_NAME}
      COMMAND bash run.sh ${SCHEME} ${SECOND} ${METHOD}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    add_test(NAME ${TEST_NAME}_visu
      COMMAND bash visu.sh ${SCHEME} ${SECOND} ${METHOD}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    set_tests_properties(${TEST_NAME}_visu PROPERTIES DEPENDS ${TEST_NAME})

  endforeach()
endforeach()
