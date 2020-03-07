cmake_minimum_required(VERSION 3.15)
project(decomposition_library)

set(CMAKE_CXX_STANDARD 20)

get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(NOT HAS_PARENT)
    add_subdirectory(test) # Build tests only for root project
endif()


file(GLOB_RECURSE CPP_SOURCE_FILES "src/*.cpp")
add_library(decomposition_library ${CPP_SOURCE_FILES})

target_include_directories(decomposition_library PUBLIC include PUBLIC lib/glm/glm)