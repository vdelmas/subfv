file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/../schemes.txt" SCHEMES)

set(SUMMARY_DEPS "")
foreach(SCHEME ${SCHEMES})
  set(TEST_NAME "gresho_summary_${SCHEME}")

  add_test(NAME ${TEST_NAME}
    COMMAND bash summary.sh ${SCHEME}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  set_tests_properties(${TEST_NAME} PROPERTIES
    DEPENDS "gresho_quad_${SCHEME}_visu;gresho_tri_${SCHEME}_visu;convergence_gresho_quad_${SCHEME}_visu;convergence_gresho_tri_${SCHEME}_visu"
    LABELS "summary")

  list(APPEND SUMMARY_DEPS ${TEST_NAME})
endforeach()

add_test(NAME gresho_summary_combine
  COMMAND bash combine.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set_tests_properties(gresho_summary_combine PROPERTIES
  DEPENDS "${SUMMARY_DEPS}"
  LABELS "summary")
