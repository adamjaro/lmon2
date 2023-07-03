
set(LIB lmon2Base)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

file(GLOB SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx)

add_library (${LIB} SHARED ${SRC})
target_link_libraries(${LIB} ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${HEPMC3_LIB})

set(SUBDIR_LIBS ${SUBDIR_LIBS} ${LIB} PARENT_SCOPE)
set(SUBDIR_INCLUDE ${SUBDIR_INCLUDE} ${CMAKE_CURRENT_SOURCE_DIR}/include PARENT_SCOPE)

