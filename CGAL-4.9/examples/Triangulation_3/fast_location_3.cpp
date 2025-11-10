#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>

#include <vector>
#include <cassert>

using namespace std;

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Delaunay_fast;
typedef Delaunay::Point Point;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Cell_handle Cell_handle;

int main()
{
  // generating points on a grid.
  std::vector<Point> P;

  //for (int z=0 ; z<20 ; z++)
  //  for (int y=0 ; y<10 ; y++)
   //   for (int x=0 ; x<10 ; x++)
  for (int i=0; i<3000; ++i)
	  //P.push_back(Point(x,y,z));
          P.push_back(Point(CGAL::get_default_random().get_double(-9, 3),
			   CGAL::get_default_random().get_double(48, 61),
			   CGAL::get_default_random().get_double(0, 60)));

  clock_t t;
  t = clock();
  // building their Delaunay triangulation.
  Delaunay T(P.begin(), P.end());
  T.insert(Point(2,30,20));

  clock_t t1 = clock();
  cout << "delaunay:" << (float(t1-t))/CLOCKS_PER_SEC << endl;

  assert( T.number_of_vertices() == 3001 );

  // performing nearest vertex queries to a series of random points,
  // which is a case where the Fast_location policy is beneficial.
  double dx = 12.0/208;
  double dy = (61.0-48.0)/192;
  double dz = 60.0/119;
  Cell_handle cell = T.locate(Point(-9,48,0));
  Cell_handle first = cell;
  Cell_handle firsti = cell;
  for (int i=0 ; i<209 ; i++){
      cell = firsti;
      double x = -9 + i*dx;
     for (int j=0 ; j<193 ; j++) {
	cell = first;
        if(j==0) cell = firsti;
	double y = 48 + j*dy;
        for (int k=0 ; k<120 ; k++){
	  double z = 0.0 + k*dz;
  	  cell = T.locate(Point(x,y,z),cell);
	  Vertex_handle vit = T.nearest_vertex(Point(x,y,z), cell);
	  //std::vector<Cell_handle> cells;
          //T.incident_cells( vit, std::back_inserter(cells) );
	  //cell = cells[0];
 	  if(k==0) first = cell;
        }
        if(j == 0) firsti = first;
     }
  }
   //T.nearest_vertex(Point(CGAL::get_default_random().get_double(-9, 3),
   //			   CGAL::get_default_random().get_double(48, 61),
   //			   CGAL::get_default_random().get_double(0, 60)));

  t = clock() - t1;
  cout << "compact location:" << (float(t))/CLOCKS_PER_SEC << endl;

  // fast location
  Delaunay_fast T1(P.begin(), P.end());

  assert( T1.number_of_vertices() == 3000 );

  t1 = clock();
  // performing nearest vertex queries to a series of random points,
  // which is a case where the Fast_location policy is beneficial.
  for (int i=0 ; i<209 ; i++)
     for (int j=0 ; j<193 ; j++) {
        for (int k=0 ; k<120 ; k++){
	  double x = -9 + i*dx;
	  double y = 48 + j*dy;
	  double z = 0.0 + k*dz;
	  Delaunay_fast::Vertex_handle vit = T1.nearest_vertex(Point(x,y,z));
        }
     }
  t = clock() - t1;
  cout << "fast location:" << (float(t))/CLOCKS_PER_SEC << endl;
  return 0;
}
