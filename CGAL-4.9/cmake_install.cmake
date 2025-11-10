# Install script for directory: /exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9

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
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/usr/share/doc/CGAL-4.9/AUTHORS;/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/usr/share/doc/CGAL-4.9/CHANGES;/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/usr/share/doc/CGAL-4.9/LICENSE;/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/usr/share/doc/CGAL-4.9/LICENSE.FREE_USE;/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/usr/share/doc/CGAL-4.9/LICENSE.GPL;/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/usr/share/doc/CGAL-4.9/LICENSE.LGPL")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/usr/share/doc/CGAL-4.9" TYPE FILE FILES
    "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/AUTHORS"
    "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/CHANGES"
    "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/LICENSE"
    "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/LICENSE.FREE_USE"
    "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/LICENSE.GPL"
    "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/LICENSE.LGPL"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/include/CGAL" REGEX "/\\.svn$" EXCLUDE)
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/include/CGAL" REGEX "/\\.svn$" EXCLUDE)
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES
    "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/scripts/cgal_create_CMakeLists"
    "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/scripts/cgal_create_cmake_script"
    "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/scripts/cgal_make_macosx_app"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/CGAL" TYPE DIRECTORY FILES "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/cmake/modules/")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/CGAL" TYPE FILE FILES "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/cmake/modules/UseCGAL.cmake")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/CGAL" TYPE FILE FILES "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/config/CGALConfig.cmake")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/usr/share/man/man1/cgal_create_cmake_script.1")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/usr/share/man/man1" TYPE FILE FILES "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/auxiliary/cgal_create_cmake_script.1")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/src/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/exports/csce/datastore/geos/groups/eddie_geos_eps_eip/geos_eps_eip_backedup/Xin/CGAL-4.9/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
