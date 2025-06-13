cmake_minimum_required(VERSION 3.10)
project(flashprediction VERSION 1.0.0)

#---Install location by default is in the build directory
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/installed" CACHE PATH "${CMAKE_BUILD_DIR}" FORCE )
endif()

#---Module directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

## Dependencies

#---CERN ROOT
find_package(ROOT CONFIG REQUIRED)
# we set the compiler version based on the one used for ROOT
execute_process(
    COMMAND bash "-c" "root-config --features"
    OUTPUT_VARIABLE TEST_ROOT_CXX_STANDARD )
if (${TEST_ROOT_CXX_STANDARD} MATCHES "cxx17")
  message( STATUS "using cxx17")
  set(CMAKE_CXX_STANDARD 17)
elseif (${TEST_ROOT_CXX_STANDARD} MATCHES "cxx14")
  message( STATUS "using cxx14")
  set(CMAKE_CXX_STANDARD 14)
else ()
  message( STATUS "using cxx11")  
  set(CMAKE_CXX_STANDARD 11)
endif()

#---Eigen
#find_package(Eigen3 REQUIRED)

#---PYTHON
#find_package( PythonInterp 3 REQUIRED )
#find_package( PythonLibs 3 REQUIRED )
#execute_process( COMMAND python3 -c "import numpy; print(numpy.get_include())" OUTPUT_VARIABLE NUMPY_INCLUDES OUTPUT_STRIP_TRAILING_WHITESPACE )

# ========================
# UBDL MODULES
# ========================

#---larbys GEO2D
# use it if env variables defined
#if (DEFINED ENV{GEO2D_BASEDIR})
#  set(LARCV_USE_GEO2D 1)
#else()
#  set(LARCV_USE_GEO2D 0)
#endif()

set( larlite_DIR "$ENV{LARLITE_LIBDIR}/cmake/larlite" )
find_package( larlite REQUIRED )

set( larcv_DIR "$ENV{LARCV_LIBDIR}/cmake/larcv")
find_package( larcv REQUIRED )

set( ublarcvapp_DIR "$ENV{UBLARCVAPP_LIBDIR}/cmake/ublarcvapp")
find_package( ublarcvapp REQUIRED )

set( larflow_DIR "$ENV{LARFLOW_LIBDIR}/cmake/larflow")
find_package( larflow REQUIRED )


# ## Git Submodules
# find_package(Git QUIET)
# if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
#   # Update submodules as needed
#   option(GIT_SUBMODULE "Check submodules during build" ON)
#   if(GIT_SUBMODULE)
#     message(STATUS "Submodule update")
#     execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
#       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#       RESULT_VARIABLE GIT_SUBMOD_RESULT)
#     if(NOT GIT_SUBMOD_RESULT EQUAL "0")
#       message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
#     endif()
#   endif()
# endif()

#---check json
# if(EXISTS "$ENV{NLOHMANN_JSON_DIR}/share/cmake/nlohmann_json/nlohmann_jsonConfig.cmake")
#   message( STATUS "Use cmvfs nlohmann/json")
#   set( nlohmann_json_DIR )
#   find_package( nlohmann_json REQUIRED )
#   set(LARCV_JSON_INCLUDE_DIR $ENV{NLOHMANN_JSON_INC})
# else()
#   if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/json/CMakeLists.txt")
#     message(FATAL_ERROR "The json submodule was not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
#   else()
#     ## External repos (from submodules)
#     message( STATUS "Use extern copy of json included in repo")
#     include_directories(extern/json)
#     add_subdirectory(extern/json)
#     set(LARCV_JSON_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/json/single_include)
#   endif()
# endif()

## Setup

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

## Add Source

#---Add base of source
add_subdirectory(flashprediction)

#---Configuration export
include(CMakePackageConfigHelpers)
set(FLASHPREDICTION_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")
set(FLASHPREDICTION_LIBRARY_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(TARGETS_EXPORT_NAME "flashpredictionTargets")
configure_package_config_file("cmake/Config.cmake.in"
  "${CMAKE_BINARY_DIR}/cmake/flashpredictionConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/flashprediction"
  PATH_VARS FLASHPREDICTION_INCLUDE_DIR FLASHPREDICTION_LIBRARY_DIR)
write_basic_package_version_file("${CMAKE_BINARY_DIR}/cmake/flashpredictionConfigVersion.cmake"
  COMPATIBILITY SameMajorVersion)


## Export Targets

#---Create info to include Flashprediction in other cmake files
install(EXPORT flashprediction
        FILE flashpredictionTargets.cmake
        DESTINATION "${CMAKE_BINARY_DIR}/installed/lib/cmake/flashprediction")

#---Create info to include Flashprediction in other cmake files
install(FILES "${CMAKE_BINARY_DIR}/cmake/flashpredictionConfig.cmake"
  "${CMAKE_BINARY_DIR}/cmake/flashpredictionConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/flashprediction")
