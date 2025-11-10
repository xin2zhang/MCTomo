set(WITH_CGAL_ImageIO "ON")

# The else condition of this code is never used in an installed
# version since it cannot happen there. Note also that for
# CMake<=2.8.11 (detected by the absence of CMP0024), the else()
# condition is never used.
if(NOT POLICY CMP0024 OR NOT CGAL_BUILDING_LIBS)
  if(NOT MSVC AND NOT CGAL_HEADER_ONLY)
    get_property(CGAL_ImageIO_LIBRARY TARGET CGAL::CGAL_ImageIO PROPERTY LOCATION)
  else()
    set(CGAL_ImageIO_LIBRARY "")
  endif()
else()
  # We are currently in a CGAL Build and CGALExports.cmake has not
  # necessarily been created yet. Just alias the targets. Also don't
  # access the LOCATION property here to set lib_LIBRARY, since those
  # targets are not imported and this is disallowed by CMP0026. Just
  # set it to the target name.
  if(TARGET CGAL_ImageIO AND NOT TARGET CGAL::CGAL_ImageIO AND NOT CGAL_HEADER_ONLY)
    add_library(CGAL::CGAL_ImageIO ALIAS CGAL_ImageIO)
    set(CGAL_ImageIO_LIBRARY CGAL::CGAL_ImageIO)
  else()
    set(CGAL_ImageIO_LIBRARY "")
  endif()
endif()


# 3RD_PARTY variables.
set(CGAL_ImageIO_3RD_PARTY_INCLUDE_DIRS   "/usr/include")
set(CGAL_ImageIO_3RD_PARTY_DEFINITIONS    "-DCGAL_USE_ZLIB")
set(CGAL_ImageIO_3RD_PARTY_LIBRARIES_DIRS "")
set(CGAL_ImageIO_3RD_PARTY_LIBRARIES      "/usr/lib64/libz.so")
