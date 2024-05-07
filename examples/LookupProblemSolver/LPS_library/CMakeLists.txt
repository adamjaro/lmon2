
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(LPS)

#BOOST
find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

#ROOT
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree MathCore EG Foam)
include(${ROOT_USE_FILE})

#top directory for LPS codes
set(LPS_TOP ${CMAKE_CURRENT_SOURCE_DIR}/../../../reconstruction)

#include and sources for LPS library
include_directories(${LPS_TOP}/include)

#sources for local library (LIB)
file(GLOB SRC
  ${LPS_TOP}/src/FoamCellFinder.cxx
  ${LPS_TOP}/src/MedCellFinder.cxx
  ${LPS_TOP}/src/SumCellFinder.cxx
  ${LPS_TOP}/src/LookupProblemSolver.cxx
)

add_library (LPS SHARED ${SRC})
target_link_libraries(LPS ${Boost_LIBRARIES} ${ROOT_LIBRARIES})

