# Install script for directory: /exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/src/CGAL

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/usr")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64" TYPE SHARED_LIBRARY FILES
    "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/src/CGAL/CMakeFiles/CMakeRelink.dir/libCGAL.so.12.0.0"
    "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/src/CGAL/CMakeFiles/CMakeRelink.dir/libCGAL.so.12"
    "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/src/CGAL/CMakeFiles/CMakeRelink.dir/libCGAL.so"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/CGAL/CGALExports.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/CGAL/CGALExports.cmake"
         "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/src/CGAL/CMakeFiles/Export/lib64/CGAL/CGALExports.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/CGAL/CGALExports-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib64/CGAL/CGALExports.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/CGAL" TYPE FILE FILES "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/src/CGAL/CMakeFiles/Export/lib64/CGAL/CGALExports.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/CGAL" TYPE FILE FILES "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/src/CGAL/CMakeFiles/Export/lib64/CGAL/CGALExports-release.cmake")
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/CGAL" TYPE FILE FILES "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/CGALLibConfig.cmake")
endif()

