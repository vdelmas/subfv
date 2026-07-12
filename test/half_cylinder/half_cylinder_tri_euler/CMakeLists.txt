file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/../../schemes.txt" SCHEMES)
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/methods.txt" METHODS)

set(SECOND ".false.")
set(METHOD "0")
foreach(SCHEME ${SCHEMES})
  set(TEST_NAME "half_cylinder_tri_euler_${SCHEME}")

  add_test(NAME ${TEST_NAME}
    COMMAND bash run.sh ${SCHEME} ${SECOND} ${METHOD}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_test(NAME ${TEST_NAME}_visu
    COMMAND bash visu.sh ${SCHEME} ${SECOND} ${METHOD}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  set_tests_properties(${TEST_NAME}_visu PROPERTIES DEPENDS ${TEST_NAME})
  set_tests_properties(${TEST_NAME} PROPERTIES LABELS "quick_test")
  set_tests_properties(${TEST_NAME}_visu PROPERTIES LABELS "quick_test")

  list(APPEND half_cylinder_tri_euler_VISU_TESTS_O1 ${TEST_NAME}_visu)

endforeach()

add_test(NAME half_cylinder_tri_euler_plot_line
  COMMAND bash plot_all_lines.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set_tests_properties(half_cylinder_tri_euler_plot_line
  PROPERTIES DEPENDS "half_cylinder_tri_euler_${VISU_TESTS_O1}")

set(SECOND ".true.")
foreach(METHOD ${METHODS})
  foreach(SCHEME ${SCHEMES})
    set(TEST_NAME "half_cylinder_tri_euler_${SCHEME}_o2_${METHOD}")

    add_test(NAME ${TEST_NAME}
      COMMAND bash run.sh ${SCHEME} ${SECOND} ${METHOD}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    add_test(NAME ${TEST_NAME}_visu
      COMMAND bash visu.sh ${SCHEME} ${SECOND} ${METHOD}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    set_tests_properties(${TEST_NAME}_visu PROPERTIES DEPENDS ${TEST_NAME})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "quick_test")
    set_tests_properties(${TEST_NAME}_visu PROPERTIES LABELS "quick_test")

    list(APPEND half_cylinder_tri_euler_o2_VISU_TESTS_O2 ${TEST_NAME}_visu)

  endforeach()
endforeach()

add_test(NAME half_cylinder_tri_euler_o2_plot_line
  COMMAND bash plot_all_lines_o2.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set_tests_properties(half_cylinder_tri_euler_o2_plot_line
  PROPERTIES DEPENDS "half_cylinder_tri_euler_o2_${VISU_TESTS_O2}")
