cmake_minimum_required(VERSION 3.15)
project(decomposition_library_test)

set(CMAKE_CXX_STANDARD 20)


add_subdirectory(lib/googletest)
include_directories(lib/googletest/googletest/include ../src)


file(GLOB_RECURSE CPP_TEST_SOURCE_FILES "src/*.cpp")
add_executable(decomposition_library_test ${CPP_TEST_SOURCE_FILES})

target_link_libraries(decomposition_library_test decomposition_library gtest gtest_main)