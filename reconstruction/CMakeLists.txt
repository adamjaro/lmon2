
set(LIB lmon2Reco)

#local and base includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_PROJECT_SOURCE_DIR}/base/include)

#specific subdirectories include
include_directories(${SUBDIR_INCLUDE})

#sources for local library (LIB)
file(GLOB SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx)

add_library (${LIB} SHARED ${SRC})
target_link_libraries(${LIB} ${Geant4_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} ${SUBDIR_LIBS})

