# cmake version to be used
cmake_minimum_required( VERSION 3.0 )

# project name
project(quopt)

# output directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ./lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ./lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin)

# target
add_executable( main src/main.cpp )
add_executable( test1 tests/test1.cpp )

# external libs
find_package(Ceres REQUIRED)
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
)

FetchContent_MakeAvailable(googletest)

target_include_directories(main
  PRIVATE
    ${CERES_INCLUDE_DIRS}
    ${MPI_INCLUDE_DIRS}
)

target_link_libraries(main
  PRIVATE
    ${CERES_LIBRARIES}
    ${MPI_C_LIBRARIES}
)

# GoogleTest testing setup
enable_testing()

target_link_libraries(test1
  gtest_main
)

include(GoogleTest)
gtest_discover_tests(test1)
