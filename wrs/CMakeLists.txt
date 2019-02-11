################################################################################
# wrs/CMakeLists.txt
#
# CMake file for weighted random sampling
#
# Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
################################################################################

# glob general sources
file(GLOB WRS_SRCS
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  *.[ch]pp
  generators/*.[ch]pp
  )

# if we have Intel MKL, add generators/mkl.hpp
if(MKL_FOUND)
#  set(WRS_SRCS ${WRS_SRCS} generators/mkl.hpp)
endif()

add_library(wrs STATIC ${WRS_SRCS})
set_target_properties(wrs PROPERTIES LINKER_LANGUAGE CXX)
target_compile_definitions(wrs PUBLIC ${WRS_DEFINITIONS})
target_include_directories(wrs PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(wrs SYSTEM PUBLIC ${WRS_INCLUDE_DIRS})
target_link_libraries(wrs ${WRS_LINK_LIBRARIES})

################################################################################
