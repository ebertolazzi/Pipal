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

file(GLOB_RECURSE TEST_ROSENBROCK_BOX "${CMAKE_CURRENT_SOURCE_DIR}/test_rosenbrock_box.cc")
add_executable(test_rosenbrock_box ${TEST_ROSENBROCK_BOX})
target_link_libraries(test_rosenbrock_box PRIVATE Pipal GTest::gtest_main)

file(GLOB_RECURSE TEST_ROSENBROCK_SUZUKI "${CMAKE_CURRENT_SOURCE_DIR}/test_rosenbrock_suzuki.cc")
add_executable(test_rosenbrock_suzuki ${TEST_ROSENBROCK_SUZUKI})
target_link_libraries(test_rosenbrock_suzuki PRIVATE Pipal GTest::gtest_main)
