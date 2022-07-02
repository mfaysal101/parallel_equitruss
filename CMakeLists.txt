add_executable(
  tc
  ${CMAKE_SOURCE_DIR}/gapbs-agile/src/tc.cc)

target_include_directories(tc PRIVATE ${asa_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/gapbs-agile/src)
target_link_libraries(tc asa::asa)
target_compile_definitions(tc PRIVATE _GLIBCXX_PARALLEL)
