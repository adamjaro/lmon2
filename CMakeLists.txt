
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(lmon2)

#executable
set (BIN run_lmon2)

#ROOT
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree MathCore EG)
include(${ROOT_USE_FILE})

#Geant
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

#BOOST
find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

#HepMC3
option(USE_HEPMC3 "Use HepMC3 input" OFF)
set(HepMC3_DIR "" CACHE PATH "HepMC3 location")
configure_file(base/include/use_hepmc3.h.in use_hepmc3.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

#set the HepMC3 if requested
if( USE_HEPMC3 )
  find_package(HepMC3 REQUIRED)
  include_directories(${HEPMC3_INCLUDE_DIR})
endif()

#subdirectory libraries and include
set(SUBDIR_LIBS)
set(SUBDIR_INCLUDE)

#subdirectories
add_subdirectory(base)
include_directories(${SUBDIR_INCLUDE})

add_subdirectory(calorimeters)
add_subdirectory(beamline)
add_subdirectory(trackers)

add_subdirectory(reconstruction)

#binary executable
add_executable(${BIN} run_lmon2.cxx)
target_link_libraries(${BIN} ${Geant4_LIBRARIES} ${Boost_LIBRARIES} ${SUBDIR_LIBS})

