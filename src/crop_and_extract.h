#ifndef _crop_and_extract_
#define _crop_and_extract_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;

/*
 * functions
 *
 */
bool ray_bbox_intersector( const K::Ray_3, const CGAL::Bbox_3, K::Point_3&);

bool crop_and_extract_points(const K::Ray_3, const K::Ray_3, const CGAL::Bbox_3,
        std::vector<K::Point_3>& );

#endif
