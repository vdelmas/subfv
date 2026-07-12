file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/../../schemes.txt" SCHEMES)

foreach(SCHEME ${SCHEMES})
  set(TEST_NAME "test_6_wall_${SCHEME}")

  add_test(NAME ${TEST_NAME}
    COMMAND bash run.sh ${SCHEME}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_test(NAME ${TEST_NAME}_visu
    COMMAND bash plot_all_schemes.sh ${SCHEME}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(${TEST_NAME}_visu PROPERTIES DEPENDS ${TEST_NAME})
  set_tests_properties(${TEST_NAME} PROPERTIES LABELS "quick_test")
  set_tests_properties(${TEST_NAME}_visu PROPERTIES LABELS "quick_test")

endforeach()
