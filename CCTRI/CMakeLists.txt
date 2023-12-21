cmake_minimum_required(VERSION 3.10)
project(TRIRG LANGUAGES CXX VERSION 0.1)
add_executable(TRIRG TRIRG.cpp)
add_library(randNums randNums.h)
set_target_properties(randNums PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(TRIRG PUBLIC randNums)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set(SEED 12)
target_include_directories(TRIRG PUBLIC
                            "${PROJECT_BINARY_DIR}"
                            )
target_compile_options(TRIRG PRIVATE -O3)
enable_testing()