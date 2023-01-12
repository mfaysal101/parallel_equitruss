add_executable(
  tc
  ${CMAKE_SOURCE_DIR}/gapbs-agile/src/tc.cc)

target_include_directories(tc PRIVATE ${CMAKE_SOURCE_DIR}/gapbs-agile/src)
target_compile_definitions(tc PRIVATE _GLIBCXX_PARALLEL)

add_executable(
  converter
  ${CMAKE_SOURCE_DIR}/gapbs-agile/src/converter.cc)

target_include_directories(converter PRIVATE ${CMAKE_SOURCE_DIR}/gapbs-agile/src)
target_compile_definitions(converter PRIVATE _GLIBCXX_PARALLEL)
