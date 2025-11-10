set(WITH_CGAL "ON")

# The else condition of this code is never used in an installed
# version since it cannot happen there. Note also that for
# CMake<=2.8.11 (detected by the absence of CMP0024), the else()
# condition is never used.
if(NOT POLICY CMP0024 OR NOT CGAL_BUILDING_LIBS)
  if(NOT MSVC AND NOT CGAL_HEADER_ONLY)
    get_property(CGAL_LIBRARY TARGET CGAL::CGAL PROPERTY LOCATION)
  else()
    set(CGAL_LIBRARY "")
  endif()
else()
  # We are currently in a CGAL Build and CGALExports.cmake has not
  # necessarily been created yet. Just alias the targets. Also don't
  # access the LOCATION property here to set lib_LIBRARY, since those
  # targets are not imported and this is disallowed by CMP0026. Just
  # set it to the target name.
  if(TARGET CGAL AND NOT TARGET CGAL::CGAL AND NOT CGAL_HEADER_ONLY)
    add_library(CGAL::CGAL ALIAS CGAL)
    set(CGAL_LIBRARY CGAL::CGAL)
  else()
    set(CGAL_LIBRARY "")
  endif()
endif()


# 3RD_PARTY variables.
set(CGAL_3RD_PARTY_INCLUDE_DIRS   "/exports/igmm/software/pkg/el7/libs/boost/1.59.0/include")
set(CGAL_3RD_PARTY_DEFINITIONS    "")
set(CGAL_3RD_PARTY_LIBRARIES_DIRS "/exports/igmm/software/pkg/el7/libs/boost/1.59.0/lib")
set(CGAL_3RD_PARTY_LIBRARIES      "/exports/igmm/software/pkg/el7/libs/boost/1.59.0/lib/libboost_thread.so;/exports/igmm/software/pkg/el7/libs/boost/1.59.0/lib/libboost_system.so")
