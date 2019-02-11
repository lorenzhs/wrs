################################################################################
# CMakeLists.txt
#
# Root CMake build script for sampler
#
# Copyright (C) 2015-2017 Timo Bingmann <tb@panthema.net>
# Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
################################################################################

cmake_minimum_required(VERSION 3.1...3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# custom cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project(wrs LANGUAGES CXX)
#set_target_properties(wrs PROPERTIES LINKER_LANGUAGE CXX)

# prohibit in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif()

# default to Debug building for single-config generators
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message("Defaulting CMAKE_BUILD_TYPE to Debug")
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type")
endif()


option(WRS_USE_LTO
  "Compile with -flto (link-time optimization)." OFF)

option(WRS_USE_MKL
  "Build with Intel MKL random generator if available." ON)

option(WRS_USE_NUMA
  "Build with NUMA awareness if available." ON)

option(WRS_USE_GSL
  "Build with GSL (GNU Scientific Library) if available." ON)


################################################################################

# variables to collect compile-time definitions, include dirs, and libraries
set(WRS_DEFINITIONS "")
set(WRS_INCLUDE_DIRS "")
set(WRS_LINK_LIBRARIES "")
set(EXEC_LINK_LIBRARIES "wrs")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT MSVC)
  # Debug builds should use -Og instead of -O0, if available
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-Og" WRS_HAS_OG)
  if(WRS_HAS_OG)
    set(CMAKE_CXX_FLAGS_DEBUG "-Og -g")
  endif()

  # Require pthread
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

  # enable -g on Release builds to aid debugging
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

  # enable warnings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra")

  # enable more warnings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic")

  # enable -march=native on Debug and Release builds
  check_cxx_compiler_flag("-march=native" WRS_HAS_MARCH_NATIVE)
  if(WRS_HAS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  endif()

  # remove -rdynamic from linker flags (smaller binaries which cannot be loaded
  # with dlopen() -- something no one needs)
  string(REGEX REPLACE "-rdynamic" ""
    CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_C_FLAGS}")
  string(REGEX REPLACE "-rdynamic" ""
    CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS}")

  # warn on conversions
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wconversion -Werror")

  # enable AddressSanitizer
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

  # enable ThreadSanitizer
  if(OFF)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=thread -pie -fPIC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DWRS_HAVE_THREAD_SANITIZER=1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread -pie -fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWRS_HAVE_THREAD_SANITIZER=1")
  endif()

  # enable UndefinedBehaviorSanitizer
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined")

  # enable extra warnings on gcc
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcast-qual")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Winit-self -Wnoexcept")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual -Wredundant-decls")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wstrict-null-sentinel -Wstrict-overflow=5")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wtautological-compare")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fipa-pure-const -Wsuggest-attribute=const")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-conversion")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-promo")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wold-style-cast")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override")
  endif()
  # enable extra warnings on clang
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wdeprecated")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wabstract-vbase-init")
  endif()

  # Clang < 3.6 0 (?) does not support debug output for auto return types yet.
  # try compiling a platform test for auto return types
  if(ON)
    include(CheckCXXSourceCompiles)
    set(OLD_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS} -g")

    check_cxx_source_compiles(
      "template <typename T> struct A { auto func(int i) { return 42 + i; } };
     int main() { A<int> a; return 0; }"
      WRS_CLANG_AUTO_RETURN_DEBUG_INFO)

    if (NOT WRS_CLANG_AUTO_RETURN_DEBUG_INFO)
      message(STATUS "compiler does not support -g debug info with auto returns")
      string(REPLACE "-g" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
      string(REPLACE "-g" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
    endif()

    set(CMAKE_REQUIRED_FLAGS "${OLD_CMAKE_REQUIRED_FLAGS}")
  endif()
elseif(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
  # raise warnings as errors
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX")

  ### disable verbose warnings:
  # warning C4589: Constructor of abstract class '...' ignores initializer for
  # virtual base class '...' (false positive warnings)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4589")
  # warning C4127: conditional expression is constant
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4127")
  # warning C4458: declaration of '...' hides class member
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4458")
  # warning C4459: declaration of '...' hides global declaration
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4459")
  # warning C4702: unreachable code
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4702")
  # warning C4250: ABC inherits XYZ via dominance
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4250")
  # warning C4503: decorated name length exceeded, name was truncated
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4503")
  # disable lots of warnings about "unsecure" C runtime function
  list(APPEND WRS_DEFINITIONS "_CRT_SECURE_NO_WARNINGS")
  # disable "The POSIX name for this item is deprecated. Instead, use the ISO C
  # and C++ conformant name.", Nope. We will not.
  list(APPEND WRS_DEFINITIONS "_CRT_NONSTDC_NO_DEPRECATE")
  # disable lots of warnings about "unsecure" STL functions
  list(APPEND WRS_DEFINITIONS "_SCL_SECURE_NO_WARNINGS")
  # windef.h bizzarly defines min and max as macros, unless this is defined.
  list(APPEND WRS_DEFINITIONS "NOMINMAX")
endif()


if(WRS_USE_LTO)
  # build with link-time optimization
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-flto CXX_HAS_LTO_FLAG)

  if(CMAKE_BUILD_TYPE MATCHES Release AND CXX_HAS_LTO_FLAG)
    find_program(CMAKE_GCC_AR
      NAMES ${_CMAKE_TOOLCHAIN_PREFIX}gcc-ar${_CMAKE_TOOLCHAIN_SUFFIX}
      HINTS ${_CMAKE_TOOLCHAIN_LOCATION})

    find_program(CMAKE_GCC_NM
      NAMES ${_CMAKE_TOOLCHAIN_PREFIX}gcc-nm
      HINTS ${_CMAKE_TOOLCHAIN_LOCATION})

    find_program(CMAKE_GCC_RANLIB
      NAMES ${_CMAKE_TOOLCHAIN_PREFIX}gcc-ranlib
      HINTS ${_CMAKE_TOOLCHAIN_LOCATION})

    if(CMAKE_GCC_AR AND CMAKE_GCC_NM AND CMAKE_GCC_RANLIB)
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto")
      set(CMAKE_AR "${CMAKE_GCC_AR}")
      set(CMAKE_NM "${CMAKE_GCC_NM}")
      set(CMAKE_RANLIB "${CMAKE_GCC_RANLIB}")
    else()
      message(WARNING "GCC indicates LTO support, but binutils wrappers could not be found. Disabling LTO.")
    endif()
  endif()
endif()

if(APPLE)
  # disable warnings about "ranlib: file: libsampling.a(...cpp.o) has no symbols"
  set(CMAKE_C_ARCHIVE_FINISH   "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
  set(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
endif()

message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

###############################################################################

# set_join(VAR "foo" "bar" "abc") sets VAR="foo bar abc"
macro(set_join var)
  set(${var})
  foreach(v ${ARGN})
    set(${var} "${${var}} ${v}")
  endforeach()
  string(STRIP ${var} "${${var}}")
endmacro(set_join)

################################################################################

# find pthreads

find_package(Threads REQUIRED)
set(WRS_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT} ${WRS_LINK_LIBRARIES})
if(CMAKE_USE_PTHREADS_INIT)
  set(WRS_LINK_LIBRARIES pthread atomic ${WRS_LINK_LIBRARIES})
endif()

# use tlx

add_subdirectory(extlib/tlx)
set(WRS_LINK_LIBRARIES tlx ${WRS_LINK_LIBRARIES})


# try to find libnuma
if(WRS_USE_NUMA)
  find_package(Numa)
  if(NUMA_FOUND)
    list(APPEND WRS_DEFINITIONS "WRS_HAVE_NUMA")

    # use Silo & Topo
    add_subdirectory(extlib)
    list(APPEND WRS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/extlib/Silo/include)
    list(APPEND WRS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/extlib/Topo/include)
    set(WRS_LINK_LIBRARIES silo topo hwloc numa ${WRS_LINK_LIBRARIES})
  else()
    message(STATUS "Could not find libnuma, disabling NUMA-awareness")
  endif()
endif()

# try to find GNU Scienific Library (gsl)
if(WRS_USE_GSL)
  find_package(GSL)
  if(GSL_FOUND)
    list(APPEND WRS_DEFINITIONS "WRS_HAVE_GSL")
    set(WRS_LINK_LIBRARIES GSL::gsl ${WRS_LINK_LIBRARIES})
  else()
    message(STATUS "Could not find GSL, no problem, not building that")
  endif()
endif()

# try to find MKL
if(WRS_USE_MKL)
  find_package(MKL)
  if(MKL_FOUND)
    #message(STATUS "Found MKL: ${MKL_ROOT_DIR}")
    list(APPEND WRS_DEFINITIONS "WRS_HAVE_MKL")
    list(APPEND WRS_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
    list(APPEND EXEC_LINK_LIBRARIES ${MKL_LP_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_SEQUENTIAL_LIBRARY})
  else()
    message(STATUS "Could not find MKL, no problem, using slower generators")
  endif()
endif()

################################################################################
### Build Macros

# macro for building main programs with correct libraries
macro(do_build TARGETNAME)
  add_executable(${TARGETNAME} ${ARGN})
  target_link_libraries(${TARGETNAME} ${EXEC_LINK_LIBRARIES})
endmacro(do_build)

# macro for building a single-source program
macro(build_prog PROGNAME)
  string(REPLACE "/" "_" TARGETNAME "${PROGNAME}") # replace slashes
  do_build(${TARGETNAME} ${PROGNAME}.cpp ${ARGN})
endmacro(build_prog)

################################################################################

add_subdirectory(wrs)

add_subdirectory(benchmark)

################################################################################
