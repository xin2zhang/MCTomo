#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "CGAL::CGAL_ImageIO" for configuration "Release"
set_property(TARGET CGAL::CGAL_ImageIO APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(CGAL::CGAL_ImageIO PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "/usr/lib64/libmpfr.so;/usr/lib64/libgmp.so;CGAL::CGAL;/exports/igmm/software/pkg/el7/libs/boost/1.59.0/lib/libboost_thread.so;/exports/igmm/software/pkg/el7/libs/boost/1.59.0/lib/libboost_system.so;/usr/lib64/libz.so"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libCGAL_ImageIO.so.12.0.0"
  IMPORTED_SONAME_RELEASE "libCGAL_ImageIO.so.12"
  )

list(APPEND _IMPORT_CHECK_TARGETS CGAL::CGAL_ImageIO )
list(APPEND _IMPORT_CHECK_FILES_FOR_CGAL::CGAL_ImageIO "${_IMPORT_PREFIX}/lib64/libCGAL_ImageIO.so.12.0.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
