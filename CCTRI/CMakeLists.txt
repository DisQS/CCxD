cmake_minimum_required(VERSION 3.10)
project(TRIRG LANGUAGES CXX VERSION 1.0)
add_executable(TRIRG TRIRG.cpp generalFunctions.cpp)
add_library(randNums randNumGen.h)
set_target_properties(randNums PROPERTIES LINKER_LANGUAGE CXX)
add_library(histogramOperations histogramOperations.h)
set_target_properties(histogramOperations PROPERTIES LINKER_LANGUAGE CXX)
add_library(RGStep RGStep.h)
add_custom_target(build-time-data-directory ALL COMMAND ${CMAKE_COMMAND} -E make_directory "./Data")
add_custom_target(build-time-analysis-directory ALL COMMAND ${CMAKE_COMMAND} -E make_directory "./Analysis")
configure_file(analyse.py "./Analysis/analyse.py" COPYONLY)
configure_file(phasemap.py "./Analysis/phasemap.py" COPYONLY)
configure_file(phaseMapMaxHeight.py "./Analysis/phaseMapMaxHeight.py" COPYONLY)
configure_file(analyseSingleAngle.py "./Analysis/analyseSingleAngle.py" COPYONLY)
configure_file(runScript.sh runScript.sh COPYONLY)
set_target_properties(RGStep PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(TRIRG PUBLIC randNums RGStep histogramOperations)
set(CMAKE_C.XX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set(SEED 11554478654588)
configure_file (
    "${PROJECT_SOURCE_DIR}/TRIRGConfig.h.in"
    "${PROJECT_BINARY_DIR}/TRIRGConfig.h"
)
target_include_directories(TRIRG PUBLIC
                            "${PROJECT_BINARY_DIR}"
                            )
target_compile_options(TRIRG PRIVATE -O3 -g)
