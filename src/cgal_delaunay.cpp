#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <cassert>
#include <math.h>
//#include <omp.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Object.h>
#include <CGAL/intersections.h>
#include <CGAL/bounding_box.h>

#include "crop_and_extract.h"

typedef struct {double vp,vs,rho;} p3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<p3, K>          Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay_fast;
typedef Delaunay::Point 					    Point;
typedef Delaunay::Segment                                           Segment;
typedef Delaunay::Ray                                               Ray;
typedef Delaunay::Vertex_handle					    Vertex_handle;
typedef Delaunay::Cell_handle					    Cell_handle;
typedef Delaunay::Edge                                              Edge;

typedef struct {double x,y,z;} d3;
/*typedef struct {
    d3 p0;
    d3 p1;
    int type;
} line3; */
//enum Bnds { PLANE_X0, PLANE_X1, PLANE_Y0, PLANE_Y1, PLANE_Z0, PLANE_Z1 };

extern "C" {

    /* build delaunay triangulation
     * input: points, parameters, and their number n
     * output: delaunay pointer
     */
    Delaunay* delaunay_build (const d3 *ppoint, const p3 *parameter, const int n){

	    std::vector<std::pair<Point,p3> > P;
	    
	    for(int i=0 ; i<n ; i++){
	        P.push_back( std::make_pair(Point(ppoint->x,ppoint->y,ppoint->z),*(parameter+i)) );
	        ppoint++;
	    }
	    // building their delaunay triangulation
	    Delaunay* pt = new Delaunay(P.begin(), P.end());
	    
	    return pt;
    }

    size_t delaunay_size(const Delaunay* pt) {
        size_t n = pt->number_of_vertices();
        return n;
    }

    void get_vertices(const Delaunay* pt, double* points, double* parameters) {

        typedef double (*arr2d_t)[ 3 ];
        arr2d_t pt2d = (arr2d_t) points;
        arr2d_t pm2d = (arr2d_t) parameters;

        int i=0;
        for(Delaunay::Finite_vertices_iterator vit= pt->finite_vertices_begin();
                vit != pt->finite_vertices_end(); vit++, i++) {
            pt2d[i][0] = vit->point().x();
            pt2d[i][1] = vit->point().y();
            pt2d[i][2] = vit->point().z();
            pm2d[i][0] = vit->info().vp;
            pm2d[i][1] = vit->info().vs;
            pm2d[i][2] = vit->info().rho;
        }

        return;
    }

    void reset_values(const Delaunay* pt, double* parameters) {

        typedef double (*arr2d_t)[ 3 ];
        arr2d_t pm2d = (arr2d_t) parameters;

        int i=0;
        for(Delaunay::Finite_vertices_iterator vit= pt->finite_vertices_begin();
                vit != pt->finite_vertices_end(); vit++, i++) {
            vit->info().vp = pm2d[i][0];
            vit->info().vs = pm2d[i][1];
            vit->info().rho = pm2d[i][2];
        }

        return;
    }

    Vertex_handle delaunay_vertex(const Delaunay* pt, const int index) {

	    Delaunay::Finite_vertices_iterator vit = pt->finite_vertices_begin();

        int i = 0;
	    for(i=1 ; i < index && vit != pt->finite_vertices_end(); i++)
	        vit++;
	    if( i != index ) {
 	        std::cerr << "Error vertex index" << std::endl;
	        exit(1);
	    }

        return vit;
    }

    Delaunay* delaunay_copy(const Delaunay* pt){
        Delaunay* pt_copy = new Delaunay(*pt); 
        return pt_copy;
    }

    Vertex_handle delaunay_vertex2(const Delaunay* pt, const Point query) {

        Vertex_handle vit;
        if( pt->is_vertex(query, vit) ) return vit;

        std::cout << "Query point is not a vertex, using nearest vertex instead!" << std::endl;

	    vit = pt->nearest_vertex(query);
        return vit;
    }

    /*
     *  convert delaunay triangulation to voronoi diagram
     *  input: delaunay, vertex_handle, boundary point
     *  output: voronoi cell which is stored as a serial of segments/rays
     *
     */
    void delaunay_voronoi(const Delaunay* pt, const Vertex_handle vt,
                        const Point p0, const Point p1, std::vector<Point>& vertices) {
        // get the  incident edges of vt
        std::vector<Edge> edges;
        pt->finite_incident_edges( vt, std::back_inserter(edges) );

        std::vector< std::vector<Segment> > voronoi;
        for( size_t i=0; i< edges.size(); i++) {
            Delaunay::Facet_circulator facets = pt->incident_facets(edges[i]);

            std::vector<Segment> voronoi_plane;
            Delaunay::Facet_circulator done = facets;
            do {
                Delaunay::Object obj = pt->dual( *facets );
                if( const Delaunay::Segment* segment = CGAL::object_cast<Delaunay::Segment>(&obj) ) {
                    voronoi_plane.push_back(*segment);
                }
                else if( const Delaunay::Ray* ray = CGAL::object_cast<Delaunay::Ray>(&obj) ) {
                    Segment segment(ray->point(0), ray->point(1));
                    voronoi_plane.push_back(segment);
                   // rays.push_back(*ray);
                   // irays.push_back( voronoi_plane.size() );
                }
                else {
                    std::cerr << "Error: the dual of facet is not recognized!" << std::endl;
                    exit(1);
                }

                facets++;
            } while( facets != done);

            // get the vertices of this facet of the voronoi cell
        }

        return;
    }
    /* 
     * convert delaunay triangulation to voronoi diagram
     * input: delaunay, vertex_handle, boudary point
     * output: voronoi which is stored asiplane in a 2d array, each row is a plane (points)
     *
     */
    void delaunay_voronoi2(const Delaunay* pt, const Vertex_handle vt, 
                          const Point point0, const Point point1, std::vector<Point>& voronoi) {
	
        // get the incident edges of vt
        std::vector<Edge> edges;
        pt->finite_incident_edges( vt, std::back_inserter(edges) );

        CGAL::Bbox_3 bbox(point0.x(),point0.y(),point0.z(),point1.x(),point1.y(),point1.z());
        for(size_t i=0; i < edges.size(); i++) {
            Delaunay::Cell_circulator cells = pt->incident_cells(edges[i]);
            while( pt->is_infinite(cells) ) cells++;

            //std::vector<Point> facets;
            std::vector<Ray> rays;
            Delaunay::Cell_circulator done = cells;
            do  {
                if( !pt->is_infinite(cells) ) {
                    voronoi.push_back( pt->dual(cells) );
                }
                else {
                    int j = 0;
                    for(j=0; j<4; j++) {
                        Vertex_handle vit = cells->vertex(j);
                        if(pt->is_infinite(vit) ) break;
                    }
                    Delaunay::Facet facet(cells,j);
                    Delaunay::Object obj= pt->dual( facet );
                    if( const Delaunay::Ray* ray = CGAL::object_cast<Delaunay::Ray>(&obj) ) {
                        rays.push_back(*ray);
                    } 

                    if(rays.size()==2){
                        std::vector<Point> bpts;
                        if(crop_and_extract_points(rays[0],rays[1],bbox,bpts)){
                            voronoi.insert( voronoi.end(), bpts.begin(), bpts.end() );
                        }
                    }

                }

                cells++;

            }while( cells != done);
            
            //voronoi.push_back(facets);
        }

        return;

    } 

    /* This function is to get voronoi vertices from delaunay
     * there is an error when dealing with boudary
     * use delaunay_voronoi2 instead
     */

    void voronoi_vertex ( const Delaunay* pt, Vertex_handle vt,const Point point0, const Point point1, 
                            std::vector<Point>& voronoi) {
        // get the incident cells of vt
        std::vector<Cell_handle> cells;
        pt->incident_cells( vt, std::back_inserter(cells) );

        CGAL::Bbox_3 bbox(point0.x(),point0.y(),point0.z(),point1.x(),point1.y(),point1.z());
        for( size_t i=0; i < cells.size(); i++) {
            if(!pt->is_infinite(cells[i])) {
                Point circumcenter = pt->dual(cells[i]);
                voronoi.push_back(circumcenter);
            }
            else {
                int j = 0;
                for(j=0; j<4; j++) {
                    Vertex_handle vit = cells[i]->vertex(j);
                    if(pt->is_infinite(vit) ) break;
                }
                Delaunay::Facet facet(cells[i],j);
                Delaunay::Object obj= pt->dual( facet );
                if( const Delaunay::Ray* ray = CGAL::object_cast<Delaunay::Ray>(&obj) ) {
		    Point bnd_ptr;
                    if (ray_bbox_intersector( *ray, bbox, bnd_ptr ) )
                        voronoi.push_back(bnd_ptr);
                }
                else {
                    std::cerr << "Error: the return of dual of one hull facet is not a ray" << std::endl;
                    exit(1);
                }
            }
        }

        return;
    }

    /* insert a point to the delaunay triangulation
     * input: delaunay, the point and parameters to be inserted
     * output: inserted voronoi
     */
    void delaunay_insert (Delaunay* pt, const d3 *pinsert, const p3 *parameter, const d3 *p0, const d3 *p1, d3 *bnd_box) {

	    Vertex_handle vt;
	    // insert the point to the delaunay triangulation
	    vt = pt->insert( Point(pinsert->x,pinsert->y,pinsert->z) );
        vt->info() = *parameter;
	
	    // get the vertex of voronoi cell
	    Point point0(p0->x,p0->y,p0->z);
	    Point point1(p1->x,p1->y,p1->z);
        std::vector<Point> circumcenter;
        //voronoi_vertex(pt, vt, point0, point1, circumcenter);
        delaunay_voronoi2(pt, vt, point0, point1, circumcenter);

        // get the bouding box of this voronoi cell
        K::Iso_cuboid_3 c3 = CGAL::bounding_box( circumcenter.begin(), circumcenter.end() );
        bnd_box->x = c3.xmin() > p0->x ? c3.xmin() : p0->x;
        bnd_box->y = c3.ymin() > p0->y ? c3.ymin() : p0->y;
        bnd_box->z = c3.zmin() > p0->z ? c3.zmin() : p0->z;
        bnd_box++;                                        
        bnd_box->x = c3.xmax() < p1->x ? c3.xmax() : p1->x;
        bnd_box->y = c3.ymax() < p1->y ? c3.ymax() : p1->y;
        bnd_box->z = c3.zmax() < p1->z ? c3.zmax() : p1->z;

	    return;
    }

    /* remove a point from the delaunay triangulaiton
     * input: delaunay class, remove point
     * output: the removed voronoi
     */
    void delaunay_remove (Delaunay* pt, const d3 *removal, const d3 *p0, const d3 *p1, p3* removed_parameter, d3 *bnd_box) {
	// find the vertex handle
        Point removed_point(removal->x,removal->y,removal->z);
        Vertex_handle vit = delaunay_vertex2(pt, removed_point);
        *removed_parameter = vit->info();

        // get the vertex of voronoi cell
	    Point point0(p0->x,p0->y,p0->z);
	    Point point1(p1->x,p1->y,p1->z);
        std::vector<Point> circumcenter;
        delaunay_voronoi2(pt, vit, point0, point1, circumcenter);
        //voronoi_vertex(pt, vit, point0, point1, circumcenter);

        // get the bouding box of this voronoi cell
        K::Iso_cuboid_3 c3 = CGAL::bounding_box( circumcenter.begin(), circumcenter.end() );
        bnd_box->x = c3.xmin() > p0->x ? c3.xmin() : p0->x;
        bnd_box->y = c3.ymin() > p0->y ? c3.ymin() : p0->y;
        bnd_box->z = c3.zmin() > p0->z ? c3.zmin() : p0->z;
        bnd_box++;
        bnd_box->x = c3.xmax() < p1->x ? c3.xmax() : p1->x;
        bnd_box->y = c3.ymax() < p1->y ? c3.ymax() : p1->y;
        bnd_box->z = c3.zmax() < p1->z ? c3.zmax() : p1->z;

        // remove the point
	    pt->remove(vit);

	    return;
    }

    /* move a point in the delaunay triangulation
     * input: delaunay class, source point, destination point
     * output: the area which has been changed
     */
    void delaunay_move (Delaunay *pt, const d3 *p_src, const d3 *p_dst, 
            const d3 *p0, const d3 *p1, d3 *bnd_box, int* verbose) {
        
	    // find the vertex handle
        Point src( p_src->x, p_src->y, p_src->z);
        Vertex_handle vit = delaunay_vertex2(pt, src);

    	// get the voronoi cell for old location p_src
    	Point point0(p0->x,p0->y,p0->z);
    	Point point1(p1->x,p1->y,p1->z);
        std::vector<Point> circumcenter;
        //voronoi_vertex(pt, vit, point0, point1, circumcenter);
        delaunay_voronoi2(pt, vit, point0, point1, circumcenter);
	
        p3 pm = vit->info();

        Point dst(p_dst->x, p_dst->y, p_dst->z);

        // if outside the boudary, return
        if(dst.x()<p0->x || dst.x()>p1->x || dst.y()<p0->y || dst.y()>p1->y ||
            dst.z()<p0->z || dst.z()>p1->z){
            *verbose = 0;
            return;
        }

    	Vertex_handle vt_dst;
    	vt_dst = pt->move_if_no_collision( vit, dst );
    	if(vt_dst!=vit){
    	    *verbose = 0;
            return;
        }

        assert(vt_dst->info().vs == pm.vs);

	    // get the voronoi cell for new location p_dst
        //voronoi_vertex(pt, vt_dst, point0, point1, circumcenter);
        delaunay_voronoi2(pt, vt_dst, point0, point1, circumcenter);

	    // calculate the area affected by move
        // get the bouding box of this voronoi cell
        //std::cout << circumcenter.size() << std::endl;
        K::Iso_cuboid_3 c3 = CGAL::bounding_box( circumcenter.begin(), circumcenter.end() );
        bnd_box->x = c3.xmin() > p0->x ? c3.xmin() : p0->x;
        bnd_box->y = c3.ymin() > p0->y ? c3.ymin() : p0->y;
        bnd_box->z = c3.zmin() > p0->z ? c3.zmin() : p0->z;
        bnd_box++;                                        
        bnd_box->x = c3.xmax() < p1->x ? c3.xmax() : p1->x;
        bnd_box->y = c3.ymax() < p1->y ? c3.ymax() : p1->y;
        bnd_box->z = c3.zmax() < p1->z ? c3.zmax() : p1->z;

        *verbose = 1;
	    return;
    }

    /* delaunay_value
     * change the parameter values of one delaunay cell
     * input: vertex_index, value disturbation
     * output: new delaunay, boudary_box
     */
    void delaunay_value(Delaunay *pt, const d3 *pt_src, const p3 *pm_dst, const d3 *p0, const d3 *p1, d3 *bnd_box, int* verbose){

	    // find the vertex handle
        Point p_src(pt_src->x, pt_src->y, pt_src->z);
        Vertex_handle vit = delaunay_vertex2(pt, p_src);
        // change the values
        vit->info() = *pm_dst;

        // get the vertex of voronoi cell
	    Point point0(p0->x,p0->y,p0->z);
	    Point point1(p1->x,p1->y,p1->z);
        std::vector<Point> circumcenter;
        //voronoi_vertex(pt, vit, point0, point1, circumcenter);
        delaunay_voronoi2(pt, vit, point0, point1, circumcenter);

        // get the bouding box of this voronoi cell
        K::Iso_cuboid_3 c3 = CGAL::bounding_box( circumcenter.begin(), circumcenter.end() );
        bnd_box->x = c3.xmin() > p0->x ? c3.xmin() : p0->x;
        bnd_box->y = c3.ymin() > p0->y ? c3.ymin() : p0->y;
        bnd_box->z = c3.zmin() > p0->z ? c3.zmin() : p0->z;
        bnd_box++;                                        
        bnd_box->x = c3.xmax() < p1->x ? c3.xmax() : p1->x;
        bnd_box->y = c3.ymax() < p1->y ? c3.ymax() : p1->y;
        bnd_box->z = c3.zmax() < p1->z ? c3.zmax() : p1->z;

        *verbose = 1;
        return;
    }

    void delaunay_box(Delaunay *pt, const d3 *pt_src, const d3 *p0, const d3 *p1, d3 *bnd_box ){

	    // find the vertex handle
        Point p_src(pt_src->x, pt_src->y, pt_src->z);
        Vertex_handle vit = delaunay_vertex2(pt, p_src);

        // get the vertex of voronoi cell
	    Point point0(p0->x,p0->y,p0->z);
	    Point point1(p1->x,p1->y,p1->z);
        std::vector<Point> circumcenter;
        //voronoi_vertex(pt, vit, point0, point1, circumcenter);
        delaunay_voronoi2(pt, vit, point0, point1, circumcenter);

        // get the bouding box of this voronoi cell
        K::Iso_cuboid_3 c3 = CGAL::bounding_box( circumcenter.begin(), circumcenter.end() );
        bnd_box->x = c3.xmin() > p0->x ? c3.xmin() : p0->x;
        bnd_box->y = c3.ymin() > p0->y ? c3.ymin() : p0->y;
        bnd_box->z = c3.zmin() > p0->z ? c3.zmin() : p0->z;
        bnd_box++;                                        
        bnd_box->x = c3.xmax() < p1->x ? c3.xmax() : p1->x;
        bnd_box->y = c3.ymax() < p1->y ? c3.ymax() : p1->y;
        bnd_box->z = c3.zmax() < p1->z ? c3.zmax() : p1->z;

        return;
    }

    void delaunay_locate (const Delaunay *pt, const d3 *ppoint, p3 *parameter ) {
    	Vertex_handle vt;
    	// locate the nearest vertex
    	vt = pt->nearest_vertex( Point(ppoint->x,ppoint->y,ppoint->z) );
    	*parameter = vt->info();
    
    	return;

    }

    /* convert delaunay triangulation to regular grid point
     * input: delaunay, boundary point, the grid number nx,ny,nz
     * output: vp, vs and rho array of each grid point
     *          vertex array of each grid point
     */
    // fast location
    void delaunay_fast_to_grid(const Delaunay *pt, const d3 *point0, const d3 *bnd0, 
                          const d3 *bnd1, const double dx, const double dy, const double dz, 
		                  const int nx, const int ny, const int nz, double *vp,
			              double *vs, double *rho) {

        typedef double (*arr3d_t)[ ny ] [ nz ];
        arr3d_t vp3d = (arr3d_t) vp;
        arr3d_t vs3d = (arr3d_t) vs;
        arr3d_t rho3d = (arr3d_t) rho;
        
        int ix0 = static_cast<int> ( std::floor(fabs(bnd0->x - point0->x)/dx) );
        int ix1 = static_cast<int> ( std::ceil(fabs(bnd1->x - point0->x)/dx) );
        int iy0 = static_cast<int> ( std::floor(fabs(bnd0->y - point0->y)/dy) );
        int iy1 = static_cast<int> ( std::ceil(fabs(bnd1->y - point0->y)/dy) );
        int iz0 = static_cast<int> ( std::floor(fabs(bnd0->z - point0->z)/dz) );
        int iz1 = static_cast<int> ( std::ceil(fabs(bnd1->z - point0->z)/dz) );

        ix1 = ix1>=nx ? nx-1 : ix1;
        iy1 = iy1>=ny ? ny-1 : iy1;
        iz1 = iz1>=nz ? nz-1 : iz1;

        /*
        double dx = (point1->x - point0->x)/(nx-1);
        double dy = (point1->y - point0->y)/(ny-1);
        double dz = (point1->z - point0->z)/(nz-1);
        */

        // loop over the grids
//        #pragma omp parallel for
        for(int i= ix0; i <= ix1; i++) 
            for(int j=iy0; j <= iy1; j++) 
                for(int k=iz0; k <= iz1; k++) {
                    double x = point0->x + i*dx;
                    double y = point0->y + j*dy;
                    double z = point0->z + k*dz;
                    Vertex_handle vit = pt->nearest_vertex( Point(x,y,z) );

                    p3 parameter = vit->info();
                    vp3d[i][j][k] = parameter.vp;
                    vs3d[i][j][k] = parameter.vs;
                    rho3d[i][j][k] = parameter.rho;
                }


        return;
    }

    // compact location
    void delaunay_to_grid(const Delaunay *pt, const d3 *point0, const d3 *bnd0, const d3 *bnd1, 
                          const double dx, const double dy, const double dz, 
                          const int nx, const int ny, const int nz,
                          double *vp, double *vs, double *rho ) {

        typedef double (*arr3d_t)[ ny ] [ nz ];
        arr3d_t vp3d = (arr3d_t) vp;
        arr3d_t vs3d = (arr3d_t) vs;
        arr3d_t rho3d = (arr3d_t) rho;
        
        int ix0 = static_cast<int> ( std::floor(fabs(bnd0->x - point0->x)/dx) );
        int ix1 = static_cast<int> ( std::ceil(fabs(bnd1->x - point0->x)/dx) );
        int iy0 = static_cast<int> ( std::floor(fabs(bnd0->y - point0->y)/dy) );
        int iy1 = static_cast<int> ( std::ceil(fabs(bnd1->y - point0->y)/dy) );
        int iz0 = static_cast<int> ( std::floor(fabs(bnd0->z - point0->z)/dz) );
        int iz1 = static_cast<int> ( std::ceil(fabs(bnd1->z - point0->z)/dz) );

        ix1 = ix1>=nx ? nx-1 : ix1;
        iy1 = iy1>=ny ? ny-1 : iy1;
        iz1 = iz1>=nz ? nz-1 : iz1;

        Cell_handle cell = pt->locate( Point(point0->x+ix0*dx,point0->y+iy0*dy,point0->z+iz0*dz) );
        Cell_handle cell_y0 = cell;
        Cell_handle cell_z0 = cell;

        // loop over the grids
        for(int i=ix0 ; i <= ix1; i++) {
            // initialize 
            cell = cell_y0;
            for(int j=iy0; j <= iy1; j++) {
                // initialize
                cell = cell_z0;
                if( j==iy0 ) cell = cell_y0;

                for(int k=iz0; k <= iz1; k++) {
                    double x = point0->x + i*dx;
                    double y = point0->y + j*dy;
                    double z = point0->z + k*dz;
                    cell = pt->locate( Point(x,y,z),cell );
                    Vertex_handle vit = pt->nearest_vertex( Point(x,y,z), cell );

                    p3 parameter = vit->info();
                    vp3d[i][j][k] = parameter.vp;
                    vs3d[i][j][k] = parameter.vs;
                    rho3d[i][j][k] = parameter.rho;
                    // if k==0, store the first z point
                    if(k==iz0) cell_z0 = cell;
                }
                // if j==0, store the first y point
                if(j==iy0) cell_y0 = cell_z0;
            }

        }

        return;
    }

    /* read delaunay triangulation from file
     */
    Delaunay* delaunay_read(const char* file) {

        std::ifstream iFileT(file,std::ios::in);
        //reading delaunay triangulation from file
        Delaunay dt;
        iFileT >> dt;
        Delaunay* pt = new Delaunay(dt);

        return pt;
    }

    /* output delaunay triangulation to file
     */
    void delaunay_write(Delaunay *pt, const char* file) {
        std::ofstream oFileT(file,std::ios::out);
        // writing file output;
        oFileT << *pt;

        return;
    }

    /* output vertex and its info of delaunay triangulation to file
     */
    void vertex_write(Delaunay *pt, const char* file) {
        std::ofstream oFileT(file,std::ios::out);
        // writing file output;
        oFileT << pt->number_of_vertices() << std::endl;
        for(Delaunay::Finite_vertices_iterator vit= pt->finite_vertices_begin();
                vit != pt->finite_vertices_end(); vit++) {
            oFileT << vit->point().x() << ' ' << vit->point().y() << ' ' << vit->point().z() << ' ';
            oFileT << vit->info().vp << ' ' << vit->info().vs << ' ' << vit->info().rho << std::endl;
        }

        return;
    }

    /* free the memory in heap
     */
    void delaunay_delete (Delaunay *pt) {
	delete pt; pt = NULL;
    }

}
