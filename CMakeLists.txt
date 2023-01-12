cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project (parallel-equitruss CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#use older method for older cmake versions
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
#link_libraries(OpenMP::OpenMP_CXX)

enable_testing()

add_subdirectory(gapbs-agile)

add_executable(
  EquiTruss
  EquiTruss.cpp
  )

target_include_directories(EquiTruss PRIVATE ${CMAKE_SOURCE_DIR})
#target_compile_definitions(EquiTruss PRIVATE _GLIBCXX_PARALLEL)
