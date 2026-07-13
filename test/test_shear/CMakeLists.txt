file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/../schemes.txt" SCHEMES)
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/meshes.txt" MESHES)

set(SHEAR_GNUPLOT_TESTS "")

foreach(SCHEME ${SCHEMES})
  foreach(MESH ${MESHES})
    set(TEST_NAME "test_shear_${SCHEME}_${MESH}")

    add_test(NAME ${TEST_NAME}
      COMMAND bash run.sh ${SCHEME} ${MESH}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    add_test(NAME ${TEST_NAME}_gnuplot
      COMMAND bash plot_shear.sh ${SCHEME} ${MESH}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    set_tests_properties(${TEST_NAME}_gnuplot
      PROPERTIES DEPENDS ${TEST_NAME})

    list(APPEND SHEAR_GNUPLOT_TESTS ${TEST_NAME}_gnuplot)

  endforeach()
endforeach()

add_test(NAME test_shear_latex
  COMMAND bash compile_shear.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set_tests_properties(test_shear_latex
  PROPERTIES DEPENDS "${SHEAR_GNUPLOT_TESTS}")

