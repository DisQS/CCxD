Cmake_minimum_required(VERSION 3.10)
project(TRIRG LANGUAGES CXX)

set(TRIRG_VERSION_MAJOR 2)
set(TRIRG_VERSION_MINOR 0)
find_package(MPI REQUIRED)
add_library(randNums CCTRI/randNumGen.hpp)
set_target_properties(randNums PROPERTIES LINKER_LANGUAGE CXX)
add_library(histogramOperations CCTRI/histogramOperations.hpp)
set_target_properties(histogramOperations PROPERTIES LINKER_LANGUAGE CXX)
add_library(RGStep CCTRI/RGStep.hpp)
set_target_properties(RGStep PROPERTIES LINKER_LANGUAGE CXX)


add_custom_target(build-time-data-directory ALL COMMAND ${CMAKE_COMMAND} -E make_directory "./Data")
add_custom_target(build-time-analysis-directory ALL COMMAND ${CMAKE_COMMAND} -E make_directory "./Analysis")
configure_file(analysis/analyse.py "./Analysis/analyse.py" COPYONLY)
configure_file(analysis/phasemap.py "./Analysis/phasemap.py" COPYONLY)
configure_file(analysis/phaseMapMaxHeight.py "./Analysis/phaseMapMaxHeight.py" COPYONLY)
configure_file(analysis/analyseSingleAngle.py "./Analysis/analyseSingleAngle.py" COPYONLY)
configure_file(scripts/2DPlotScript.sh 2DPlotScript.sh COPYONLY)


add_executable(TRIRG CCTRI/TRIRG.cpp CCTRI/generalFunctions.cpp)
add_executable(CC2D CCTRI/CC2D.cpp CCTRI/generalFunctions.cpp)

target_include_directories(TRIRG PUBLIC "${PROJECT_SOURCE_DIR}/CCTRI")
target_include_directories(CC2D PUBLIC "${PROJECT_SOURCE_DIR}/CCTRI")

target_link_libraries(TRIRG PUBLIC randNums RGStep histogramOperations MPI::MPI_CXX)
target_link_libraries(CC2D PUBLIC randNums RGStep histogramOperations MPI::MPI_CXX)

set(CMAKE_C.XX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED true)
set(SEED 27)


configure_file (
    "${PROJECT_SOURCE_DIR}/CCTRI/TRIRGConfig.hpp.in"
    "${PROJECT_BINARY_DIR}/CCTRI/TRIRGConfig.hpp"
)
target_include_directories(TRIRG PUBLIC
                            "${PROJECT_BINARY_DIR}"
                            )
target_compile_options(TRIRG PRIVATE -O3 -g)
