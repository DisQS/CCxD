cmake_minimum_required(VERSION 3.10)
project(TRIRG VERSION 0.1)
add_executable(TRIRG TRIRG.cpp)
add_library(randNums randNums.cpp)
target_link_libraries(TRIRG PUBLIC randNums)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED true)
configure_file(TRIRGConfig.h.in TRIRGConfig.h)
target_include_directories(TRIRG PUBLIC
                            "${PROJECT_BINARY_DIR}"
                            )
target_compile_options(TRIRG PRIVATE -O3)