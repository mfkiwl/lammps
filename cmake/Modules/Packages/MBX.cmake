# MBX v1.3.0+ support for MBX package
# Based off of PLUMED.cmake

# set policy to silence warnings about timestamps of downloaded files. review occasionally if it may be set to NEW
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 OLD)
endif()

# for supporting multiple concurrent mbx installations for debugging and testing
set(MBX_SUFFIX "" CACHE STRING "Suffix for MBX library")
mark_as_advanced(MBX_SUFFIX)

set(MBX_CONFIG_CC  ${CMAKE_C_COMPILER})
set(MBX_CONFIG_CXX  ${CMAKE_CXX_COMPILER})
if(BUILD_MPI)
  set(MBX_CONFIG_MPI "--enable-mpi")
else()
  set(MBX_CONFIG_MPI "--disable-mpi")
endif()

set(MBX_MODE "static" CACHE STRING "Linkage mode for MBX library")
set(MBX_MODE_VALUES static shared runtime)
set_property(CACHE MBX_MODE PROPERTY STRINGS ${MBX_MODE_VALUES})
validate_option(MBX_MODE MBX_MODE_VALUES)
string(TOUPPER ${MBX_MODE} MBX_MODE)

set(MBX_LINK_LIBS)
if(MBX_MODE STREQUAL "STATIC")
  find_package(FFTW3 QUIET)
  if(FFTW3_FOUND)
    list(APPEND MBX_LINK_LIBS FFTW3::FFTW3)
  endif()
endif()

find_package(PkgConfig QUIET)
set(DOWNLOAD_MBX_DEFAULT ON)
if(PKG_CONFIG_FOUND AND NOT DEFINED MBX_GIT_TAG)
  pkg_check_modules(MBX QUIET mbx)
  if(MBX_FOUND)
    set(DOWNLOAD_MBX_DEFAULT OFF)
  endif()
endif()

option(DOWNLOAD_MBX "Download MBX package instead of using an already installed one" ${DOWNLOAD_MBX_DEFAULT})
if(DOWNLOAD_MBX)
  message(STATUS "MBX download requested - we will build our own")
  
  if(MBX_GIT_TAG)
    message(STATUS "Using user-specified MBX branch/tag: ${MBX_GIT_TAG}")
    set(MBX_GIT_TAG ${MBX_GIT_TAG})
  else()
    message(STATUS "Using default MBX tag: v1.3.2")
    set(MBX_GIT_TAG v1.3.2)
  endif()

  set(MBX_BUILD_BYPRODUCTS "<INSTALL_DIR>/lib/${CMAKE_STATIC_LIBRARY_PREFIX}mbx${CMAKE_STATIC_LIBRARY_SUFFIX}")

  # MBX perform autoreconf -fi

  message(STATUS "MBX_CONFIG_MPI: ${MBX_CONFIG_MPI}")
  message(MBX_CONFIG_CXX: ${MBX_CONFIG_CXX})
  message(CPPFLAGS: ${MBX_CONFIG_CPPFLAGS})

  include(ExternalProject)
  ExternalProject_Add(mbx_build
    GIT_REPOSITORY https://github.com/paesanilab/MBX.git
    GIT_TAG ${MBX_GIT_TAG}
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
                      ${MBX_CONFIG_MPI}
                      CXX=${MBX_CONFIG_CXX}
                      CC=${MBX_CONFIG_CC}
                      # CPPFLAGS=${MBX_CONFIG_CPPFLAGS}
    BUILD_BYPRODUCTS ${MBX_BUILD_BYPRODUCTS}
  )
  ExternalProject_get_property(mbx_build INSTALL_DIR)
  add_library(LAMMPS::MBX UNKNOWN IMPORTED)
  add_dependencies(LAMMPS::MBX mbx_build)
  set_target_properties(LAMMPS::MBX PROPERTIES
    IMPORTED_LOCATION ${INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}mbx${CMAKE_STATIC_LIBRARY_SUFFIX}
    INTERFACE_LINK_LIBRARIES "${MBX_LINK_LIBS};${CMAKE_DL_LIBS}"
  )

  set_target_properties(LAMMPS::MBX PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_DIR}/include)
  file(MAKE_DIRECTORY ${INSTALL_DIR}/include)
  if(CMAKE_PROJECT_NAME STREQUAL "lammps")
    target_link_libraries(lammps PRIVATE LAMMPS::MBX)
  endif()
else()
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(MBX REQUIRED mbx${MBX_SUFFIX})
  add_library(LAMMPS::MBX INTERFACE IMPORTED)

  file(MAKE_DIRECTORY ${MBX_INCLUDE_DIRS})

  target_include_directories(LAMMPS::MBX INTERFACE ${MBX_INCLUDE_DIRS})
  target_link_directories(LAMMPS::MBX INTERFACE ${MBX_LIBRARY_DIRS})
  target_link_libraries(LAMMPS::MBX INTERFACE "${MBX_LINK_LIBS};${MBX_LIBRARIES};${CMAKE_DL_LIBS}")

  if(CMAKE_PROJECT_NAME STREQUAL "lammps")
    target_link_libraries(lammps PUBLIC LAMMPS::MBX)
  endif()
endif()

