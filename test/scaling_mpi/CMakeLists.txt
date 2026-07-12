set(TEST_NAME "scaling_mpi")

add_test(NAME ${TEST_NAME}
  COMMAND bash run.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
