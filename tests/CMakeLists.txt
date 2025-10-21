# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The Pipal project is distributed under the MIT License.                                         #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# e-mail: davide.stocco@unitn.it                               e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include(PipalGTest)

file(GLOB_RECURSE TEST_BASIC "${CMAKE_CURRENT_SOURCE_DIR}/test_basic.cc")
add_executable(test_basic ${TEST_BASIC})
target_link_libraries(test_basic PRIVATE Pipal GTest::gtest_main)

file(GLOB_RECURSE TEST_LINEAR "${CMAKE_CURRENT_SOURCE_DIR}/test_linear.cc")
add_executable(test_linear ${TEST_LINEAR})
target_link_libraries(test_linear PRIVATE Pipal GTest::gtest_main)
