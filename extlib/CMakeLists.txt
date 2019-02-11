################################################################################
# CMakeLists.txt
#
# Build Silo & Topo (they only ship makefiles)
#
# Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
################################################################################

# Topo
file(GLOB TOPO_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  Topo/include/*.h
  Topo/source/*.cpp
  )

add_library(topo STATIC ${TOPO_SRCS})
set_target_properties(topo PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(topo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Topo/include)


# Silo
file(GLOB SILO_SRCS_GENERIC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  Silo/include/*.h
  Silo/include/silo/*.h
  Silo/source/consume.cpp
  Silo/source/pointermap.cpp
  Silo/source/osmemory.cpp
  Silo/source/silo.cpp
  )

if (WIN32)
  set(SILO_SRCS ${SILO_SRCS_GENERIC} Silo/source/osmemory-windows.cpp)
else()
  set(SILO_SRCS ${SILO_SRCS_GENERIC} Silo/source/osmemory-linux.cpp)
endif()

add_library(silo STATIC ${SILO_SRCS})
set_target_properties(silo PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(silo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Silo/include/silo ${CMAKE_CURRENT_SOURCE_DIR}/Topo/include)
