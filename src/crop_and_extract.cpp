#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Object.h>
#include <CGAL/intersections.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef K::Ray_3 Ray;
typedef K::Vector_3 Vector;
typedef K::Plane_3 Plane;

bool intersection_plane_bbox(const Plane plane, const CGAL::Bbox_3 bbox,
        std::vector<Point>& points) {

    // along x axis
    Segment segment(Point(bbox.xmin(),bbox.ymin(),bbox.zmin()),Point(bbox.xmax(),bbox.ymin(),bbox.zmin()));
    if( do_intersect(plane, segment) ){
        CGAL::cpp11::result_of<K::Intersect_3(Plane, Segment)>::type
            result = intersection(plane, segment);
        if( const Point* s = boost::get<Point>(&*result) ){
            points.push_back(*s);
        }            
        else{
            std::cout << "Warning: intersection of a segment and the plane is not a point" << std::endl;
        }
    }

    segment = Segment(Point(bbox.xmin(),bbox.ymax(),bbox.zmin()),Point(bbox.xmax(),bbox.ymax(),bbox.zmin()));
    if( do_intersect(plane, segment) ){
        CGAL::cpp11::result_of<K::Intersect_3(Plane, Segment)>::type
            result = intersection(plane, segment);
        if( const Point* s = boost::get<Point>(&*result) ){
            points.push_back(*s);
        }            
        else{
            std::cout << "Warning: intersection of a segment and the plane is not a point" << std::endl;
        }
    }

    segment = Segment(Point(bbox.xmin(),bbox.ymin(),bbox.zmax()),Point(bbox.xmax(),bbox.ymin(),bbox.zmax()));
    if( do_intersect(plane, segment) ){
        CGAL::cpp11::result_of<K::Intersect_3(Plane, Segment)>::type
            result = intersection(plane, segment);
        if( const Point* s = boost::get<Point>(&*result) ){
            points.push_back(*s);
        }            
        else{
            std::cout << "Warning: intersection of a segment and the plane is not a point" << std::endl;
        }
    }

    segment = Segment(Point(bbox.xmin(),bbox.ymax(),bbox.zmax()),Point(bbox.xmax(),bbox.ymax(),bbox.zmax()));
    if( do_intersect(plane, segment) ){
        CGAL::cpp11::result_of<K::Intersect_3(Plane, Segment)>::type
            result = intersection(plane, segment);
        if( const Point* s = boost::get<Point>(&*result) ){
            points.push_back(*s);
        }            
        else{
            std::cout << "Warning: intersection of a segment and the plane is not a point" << std::endl;
        }
    }

    // along y axis
    segment = Segment(Point(bbox.xmin(),bbox.ymin(),bbox.zmin()),Point(bbox.xmin(),bbox.ymax(),bbox.zmin()));
    if( do_intersect(plane, segment) ){
        CGAL::cpp11::result_of<K::Intersect_3(Plane, Segment)>::type
            result = intersection(plane, segment);
        if( const Point* s = boost::get<Point>(&*result) ){
            points.push_back(*s);
        }            
        else{
            std::cout << "Warning: intersection of a segment and the plane is not a point" << std::endl;
        }
    }

    segment = Segment(Point(bbox.xmax(),bbox.ymin(),bbox.zmin()),Point(bbox.xmax(),bbox.ymax(),bbox.zmin()));
    if( do_intersect(plane, segment) ){
        CGAL::cpp11::result_of<K::Intersect_3(Plane, Segment)>::type
            result = intersection(plane, segment);
        if( const Point* s = boost::get<Point>(&*result) ){
            points.push_back(*s);
        }            
        else{
            std::cout << "Warning: intersection of a segment and the plane is not a point" << std::endl;
        }
    }

    segment = Segment(Point(bbox.xmin(),bbox.ymin(),bbox.zmax()),Point(bbox.xmin(),bbox.ymax(),bbox.zmax()));
    if( do_intersect(plane, segment) ){
        CGAL::cpp11::result_of<K::Intersect_3(Plane, Segment)>::type
            result = intersection(plane, segment);
        if( const Point* s = boost::get<Point>(&*result) ){
            points.push_back(*s);
        }            
        else{
            std::cout << "Warning: intersection of a segment and the plane is not a point" << std::endl;
        }
    }

    segment = Segment(Point(bbox.xmax(),bbox.ymin(),bbox.zmax()),Point(bbox.xmax(),bbox.ymax(),bbox.zmax()));
    if( do_intersect(plane, segment) ){
        CGAL::cpp11::result_of<K::Intersect_3(Plane, Segment)>::type
            result = intersection(plane, segment);
        if( const Point* s = boost::get<Point>(&*result) ){
            points.push_back(*s);
        }            
        else{
            std::cout << "Warning: intersection of a segment and the plane is not a point" << std::endl;
        }
    }

    // along z axis
    segment = Segment(Point(bbox.xmin(),bbox.ymin(),bbox.zmin()),Point(bbox.xmin(),bbox.ymin(),bbox.zmax()));
    if( do_intersect(plane, segment) ){
        CGAL::cpp11::result_of<K::Intersect_3(Plane, Segment)>::type
            result = intersection(plane, segment);
        if( const Point* s = boost::get<Point>(&*result) ){
            points.push_back(*s);
        }            
        else{
            std::cout << "Warning: intersection of a segment and the plane is not a point" << std::endl;
        }
    }

    segment = Segment(Point(bbox.xmax(),bbox.ymin(),bbox.zmin()),Point(bbox.xmax(),bbox.ymin(),bbox.zmax()));
    if( do_intersect(plane, segment) ){
        CGAL::cpp11::result_of<K::Intersect_3(Plane, Segment)>::type
            result = intersection(plane, segment);
        if( const Point* s = boost::get<Point>(&*result) ){
            points.push_back(*s);
        }            
        else{
            std::cout << "Warning: intersection of a segment and the plane is not a point" << std::endl;
        }
    }

    segment = Segment(Point(bbox.xmin(),bbox.ymax(),bbox.zmin()),Point(bbox.xmin(),bbox.ymax(),bbox.zmax()));
    if( do_intersect(plane, segment) ){
        CGAL::cpp11::result_of<K::Intersect_3(Plane, Segment)>::type
            result = intersection(plane, segment);
        if( const Point* s = boost::get<Point>(&*result) ){
            points.push_back(*s);
        }            
        else{
            std::cout << "Warning: intersection of a segment and the plane is not a point" << std::endl;
        }
    }

    segment = Segment(Point(bbox.xmax(),bbox.ymax(),bbox.zmin()),Point(bbox.xmax(),bbox.ymax(),bbox.zmax()));
    if( do_intersect(plane, segment) ){
        CGAL::cpp11::result_of<K::Intersect_3(Plane, Segment)>::type
            result = intersection(plane, segment);
        if( const Point* s = boost::get<Point>(&*result) ){
            points.push_back(*s);
        }            
        else{
            std::cout << "Warning: intersection of a segment and the plane is not a point" << std::endl;
        }
    }

    if(points.empty()){
        return false;
    }
    else {
        return true;
    }
}
    

/*
    Point get_intersector( const Point p0, const Point p1, const Point point0, const Point point1) {

        Point inter_x = (p1.x()-p0.x())>0 ? p0+(p1-p0)*(point1.x()-p0.x())/(p1.x()-p0.x()) : p0+(p1-p0)*(point0.x()-p0.x())/(p1.x()-p0.x());
        Point inter_y = (p1.y()-p0.y())>0 ? p0+(p1-p0)*(point1.y()-p0.y())/(p1.y()-p0.y()) : p0+(p1-p0)*(point0.y()-p0.y())/(p1.y()-p0.y());
        Point inter_z = (p1.z()-p0.z())>0 ? p0+(p1-p0)*(point1.z()-p0.z())/(p1.z()-p0.z()) : p0+(p1-p0)*(point0.z()-p0.z())/(p1.z()-p0.z());

        Point small = CGAL::has_smaller_distance_to_point(p0, inter_x, inter_y) ? inter_x : inter_y;
        small = CGAL::has_smaller_distance_to_point(p0, small, inter_z) ? small : inter_z;

        return small;
    }

    bool get_corner( const Point p0, const Point p1, const Point p2, 
            const Point point0, const Point point1, Point& crn_ptr) {
        return true;
    }
    */

    bool ray_bbox_intersector( const Ray ray, const CGAL::Bbox_3 bbox, Point &bnd_ptr) {

        if( do_intersect(bbox, ray) ) {
            CGAL::cpp11::result_of<K::Intersect_3(CGAL::Bbox_3, K::Ray_3)>::type
             result = intersection(bbox,ray);
            if( const K::Segment_3* s = boost::get<K::Segment_3>(&*result) ) {
                bnd_ptr = s->target();
                return true;
            }
            else {
                const Point* p = boost::get<Point> (&*result);
                return false;
            }

        }
        else {
            //std::cerr << "Warning ray does not intersect with bbox!" << std::endl;
            return false;
        }

    }

    /* tell a point is in a polygon with an infinite vertex
     */
    bool within_polygon(const Point point, const Ray ray1, const Ray ray2){

        // in the same side of ray1
        Vector vector_1( ray1 );
        Vector vector_1_2( ray1.point(0),ray2.point(1) );
        Vector vector_1_p( ray1.point(0), point );
        // in the same side of ray2
        Vector vector_2( ray2 );
        Vector vector_2_1( ray2.point(0), ray1.point(1) );
        Vector vector_2_p( ray2.point(0), point );
        // in the same side of segment between 1 and 2
        bool same_side_of_ray12 = true;
        if( ray1.point(0) != ray2.point(0) ){
            Vector vector_3( ray1.point(0), ray2.point(0) );
            if( CGAL::cross_product(vector_3, vector_1_2) * CGAL::cross_product(vector_3, vector_1_p) > 0 ){
                same_side_of_ray12 = true;
            }
            else {
                same_side_of_ray12 = false;
            }
        }

        if( CGAL::cross_product(vector_1, vector_1_2) * CGAL::cross_product(vector_1, vector_1_p) > 0 &&
            CGAL::cross_product(vector_2, vector_2_1) * CGAL::cross_product(vector_2, vector_2_p) > 0 &&
            same_side_of_ray12 )
            return true;

        return false;
    }
    /* crop one facet of Voronoi polyhedra with a axis-lined bouding box
     * and extract its boundary points
     */
    bool crop_and_extract_points(const Ray ray1, const Ray ray2, const CGAL::Bbox_3 bbox,
                              std::vector<Point>& bpts) {

        /*if( !CGAL::coplanar(ray1.point(0),ray1.point(1),ray2.point(0),ray2.point(1)) ){
            std::cout << ray1 << ray2 << std::endl;
            std::cerr << "Error: two rays are supposed to be coplanar from one facet of polyhedra!" <<std::endl;
            //exit(1);
        }
        */

        std::vector<Point> ppoints;
        Plane plane(ray1, ray2.point(1));
        if( intersection_plane_bbox(plane, bbox, ppoints) ){
            for(std::vector<Point>::iterator it=ppoints.begin(); it!=ppoints.end(); ++it){
                if(within_polygon(*it,ray1,ray2)) bpts.push_back(*it);
            }
        }

        Point rpoint;
        if( ray_bbox_intersector(ray1,bbox,rpoint) ) bpts.push_back(rpoint);
        if( ray_bbox_intersector(ray2,bbox,rpoint) ) bpts.push_back(rpoint);
        //bpts.insert( bpts.end(), ppoints.begin(), ppoints.end() );

        if(bpts.empty()){
            return false;
        }
        else {
            return true;
        }

    }
