cmake_minimum_required(VERSION 2.4.7)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/")
set(BIN_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/bin")

Include (CheckCXXSourceCompiles)

#SET(CMAKE_CXX_FLAGS "-Wall" )
#SET(CMAKE_CXX_FLAGS_DEBUG "-Werror -O0 -g -DDEBUG_FITS" )
#SET(CMAKE_C_FLAGS "-Wall" )
#SET(CMAKE_C_FLAGS_DEBUG "-Werror -O0 -g -DDEBUG_FITS" )


find_package(INDI REQUIRED)
find_package(Nova REQUIRED)
find_package(ZLIB REQUIRED)
#find_package(Threads REQUIRED)

set(EQMOD_VERSION_MAJOR 0)
set(EQMOD_VERSION_MINOR 1)
option(WITH_SIMULATOR "Enable Simulator" ON)
option(WITH_LOGGER "Enable File Logging" OFF)
option(WITH_NOFMANY "Use NOFMANY switches" OFF)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h )

include_directories( ${CMAKE_CURRENT_BINARY_DIR})
include_directories( ${CMAKE_SOURCE_DIR})
include_directories( ${INDI_INCLUDE_DIR})
include_directories( ${NOVA_INCLUDE_DIR})

########### EQMod ###############
set(eqmod_SRCS
   ${CMAKE_SOURCE_DIR}/eqmod.cpp ${CMAKE_SOURCE_DIR}/eqmoderror.cpp ${CMAKE_SOURCE_DIR}/skywatcher.cpp 
   ${CMAKE_SOURCE_DIR}/logger/Logger.cpp
   ${CMAKE_SOURCE_DIR}/simulator/simulator.cpp    ${CMAKE_SOURCE_DIR}/simulator/skywatcher-simulator.cpp
   ${CMAKE_SOURCE_DIR}/align/align.cpp ${CMAKE_SOURCE_DIR}/align/pointset.cpp ${CMAKE_SOURCE_DIR}/align/htm.c)

add_executable(indi_eqmod_telescope ${eqmod_SRCS})

target_link_libraries(indi_eqmod_telescope ${INDI_LIBRARIES} ${INDI_DRIVER_LIBRARIES} ${NOVA_LIBRARIES} z)

install(TARGETS indi_eqmod_telescope RUNTIME DESTINATION bin )

if (WITH_SIMULATOR)
  install( FILES  indi_eqmod.xml indi_eqmod_sk.xml align/indi_align_sk.xml simulator/indi_eqmod_simulator_sk.xml DESTINATION ${INDI_DATA_DIR})
else(WITH_SIMULATOR)
  install( FILES  indi_eqmod.xml indi_eqmod_sk.xml align/indi_align_sk.xml DESTINATION ${INDI_DATA_DIR})
endif(WITH_SIMULATOR)